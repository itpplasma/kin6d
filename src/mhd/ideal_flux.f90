module milca_ideal_flux
    ! Fixed-boundary ideal-MHD stiffness and mass on nested flux surfaces
    ! (docs/stability.md, P1e.3). Coordinates (rho, theta, zeta), rho the
    ! radial label with the axis at 0 and the wall at 1. Unknowns are the
    ! contravariant densities D^u = sqrt(g) xi^u of one mode family,
    ! D^u = sum_h D^u_h(rho) exp(i(m_h theta + n_h zeta)): D^rho_h in S_p and
    ! D^theta_h, D^zeta_h in the derivative space S_(p-1), the hybrid pair of
    ! milca_ideal_cylinder. Harmonics couple through the angular Fourier
    ! coefficients of the pointwise form, evaluated on a uniform angle grid.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortnum_quadrature, only: gauss_legendre_ab
    use milca_real_identity, only: same_real_value
    use fortfem_bspline_feec, only: evaluate_bspline_basis
    use fortsparse, only: csc_z_t, csc_from_triplet, fortsparse_status_t, &
        FORTSPARSE_OK, status_set, FORTSPARSE_INVALID_MATRIX
    use milca_ideal_flux_leaves, only: flux_jets_array, flux_form_matrix_array
    implicit none
    private
    public :: flux_space_t, flux_space_init, flux_assemble, flux_equilibrium_t
    public :: flux_point_t, flux_point

    type :: flux_space_t
        integer :: degree = 0, nh = 0, nr = 0, ndof = 0, ntheta = 0, nzeta = 0
        integer :: nfp = 1
        integer, allocatable :: m(:), n(:)
        real(dp), allocatable :: knots(:), edges(:)
        integer, allocatable :: radial_dof(:, :), theta_dof(:, :), zeta_dof(:, :)
    end type flux_space_t

    !> Pointwise equilibrium data at (rho, theta, zeta): Cartesian chart jets,
    !> toroidal flux per radian psi(rho), iota, lambda jets, pressure, density.
    type :: flux_point_t
        real(dp) :: xu(3, 3) = 0, xuu(3, 3, 3) = 0
        real(dp) :: psi1 = 0, psi2 = 0, iota = 0, iota1 = 0
        !> lambda_theta, lambda_zeta, then d_rho, d_theta, d_zeta of
        !> lambda_theta and d_rho, d_zeta of lambda_zeta.
        real(dp) :: lam(7) = 0
        real(dp) :: p = 0, pp = 0, density = 0
    end type flux_point_t

    !> An equilibrium on nested surfaces, evaluated pointwise.
    type, abstract :: flux_equilibrium_t
    contains
        procedure(flux_evaluate), deferred :: evaluate
    end type flux_equilibrium_t

    abstract interface
        subroutine flux_evaluate(self, rho, theta, zeta, point)
            import :: dp, flux_point_t, flux_equilibrium_t
            class(flux_equilibrium_t), intent(in) :: self
            real(dp), intent(in) :: rho, theta, zeta
            type(flux_point_t), intent(out) :: point
        end subroutine flux_evaluate
    end interface
contains
    subroutine flux_space_init(space, edges, degree, m, n, nfp, ntheta, nzeta, &
            status)
        ! D^rho vanishes on the axis and the wall, D^zeta on the axis, D^theta
        ! on the axis unless |m| = 1, as for the cylinder densities.
        type(flux_space_t), intent(out) :: space
        real(dp), intent(in) :: edges(0:)
        integer, intent(in) :: degree, m(:), n(:), nfp, ntheta, nzeta
        integer, intent(out) :: status
        integer :: ne, k, h, at, next
        status = 1
        ne = size(edges) - 1
        if (degree < 1 .or. ne < 1 .or. size(m) < 1) return
        if (size(n) /= size(m) .or. nfp < 1 .or. ntheta < 1 .or. nzeta < 1) return
        if (.not. same_real_value(edges(0), 0.0_dp)) return
        if (.not. same_real_value(edges(ne), 1.0_dp)) return
        if (any(edges(1:) <= edges(:ne-1))) return
        space%degree = degree; space%nh = size(m); space%nfp = nfp
        space%m = m; space%n = n; space%edges = edges
        space%ntheta = ntheta; space%nzeta = nzeta
        allocate(space%knots(2*(degree + 1) + (ne - 1)*degree))
        space%knots(1:degree + 1) = 0
        at = degree + 1
        do k = 1, ne - 1
            space%knots(at + 1:at + degree) = edges(k)
            at = at + degree
        end do
        space%knots(at + 1:at + degree + 1) = 1
        space%nr = size(space%knots) - degree - 1
        allocate(space%radial_dof(space%nr, space%nh), &
            space%theta_dof(space%nr - 1, space%nh), &
            space%zeta_dof(space%nr - 1, space%nh))
        space%radial_dof = 0; space%theta_dof = 0; space%zeta_dof = 0
        ! Radial-major numbering keeps the matrix banded: all harmonics and
        ! components of one radial basis index are adjacent.
        next = 0
        do k = 1, space%nr
            do h = 1, space%nh
                if (k >= 2 .and. k <= space%nr - 1) then
                    next = next + 1; space%radial_dof(k, h) = next
                end if
                if (k <= space%nr - 1) then
                    if (k /= 1 .or. abs(m(h)) == 1) then
                        next = next + 1; space%theta_dof(k, h) = next
                    end if
                    if (k >= 2) then
                        next = next + 1; space%zeta_dof(k, h) = next
                    end if
                end if
            end do
        end do
        space%ndof = next
        status = 0
    end subroutine flux_space_init

    subroutine flux_point(equilibrium, rho, theta, zeta, gam, mu, weight_form, &
            weight_mass, jac, status)
        ! Pointwise coefficient matrices in the density-jet order
        ! (D(3), D_rho(3), D_theta(3), D_zeta(3)): T^T C T and rho T0^T T0.
        class(flux_equilibrium_t), intent(in) :: equilibrium
        real(dp), intent(in) :: rho, theta, zeta, gam, mu
        real(dp), intent(out) :: weight_form(12, 12), weight_mass(3, 3), jac
        integer, intent(out) :: status
        type(flux_point_t) :: e
        real(dp) :: b(3), g(3, 3), q(3), t(12, 12), c(12, 12)
        status = 1
        call equilibrium%evaluate(rho, theta, zeta, e)
        call flux_jets_array(e%xu, e%xuu, e%psi1, e%psi2, e%iota, e%iota1, &
            e%pp, e%lam, jac, b, g, q, t)
        call flux_form_matrix_array(b, g, q, e%p, gam, mu, c)
        weight_form = matmul(transpose(t), matmul(c, t))
        weight_mass = e%density*matmul(transpose(t(1:3, 1:3)), t(1:3, 1:3))
        if (.not. all(ieee_is_finite(weight_form))) return
        if (.not. all(ieee_is_finite(weight_mass))) return
        if (.not. ieee_is_finite(jac)) return
        if (abs(jac) <= tiny(jac)) return
        status = 0
    end subroutine flux_point

    subroutine flux_assemble(space, equilibrium, gam, mu, stiffness, mass, status)
        ! Weights: Gauss in rho, trapezoidal in theta and zeta (one field
        ! period), times |sqrt(g)|. The common factor nfp cancels. Harmonic
        ! pairs couple through the angular Fourier coefficient of T^T C T at
        ! their difference (dm, dn); each basis jet has four nonzero entries.
        type(flux_space_t), intent(in) :: space
        class(flux_equilibrium_t), intent(in) :: equilibrium
        real(dp), intent(in) :: gam, mu
        type(csc_z_t), intent(out) :: stiffness, mass
        type(fortsparse_status_t), intent(out) :: status
        real(dp), parameter :: twopi = 2*acos(-1.0_dp)
        integer, allocatable :: rows(:), cols(:), pair(:, :), dlist(:, :)
        complex(dp), allocatable :: kval(:), mval(:), wave(:, :)
        real(dp), allocatable :: angle(:, :)
        real(dp) :: wa
        integer :: p_, nq, ne, na, nb, nloc, ntrip, nd, e, i, j, h1, h2, d, at, mr
        integer :: bad
        call status_set(status, FORTSPARSE_INVALID_MATRIX, 'invalid flux input')
        if (space%ndof < 1) return
        if (.not. ieee_is_finite(mu) .or. mu <= 0) return
        if (.not. ieee_is_finite(gam) .or. gam < 0) return
        p_ = space%degree; nq = p_ + 3; ne = size(space%edges) - 1
        na = space%ntheta*space%nzeta
        nb = 3*p_ + 1; nloc = nb*space%nh
        ntrip = nloc*nloc*ne
        allocate(rows(ntrip), cols(ntrip), kval(ntrip), mval(ntrip), angle(2, na), &
            pair(space%nh, space%nh), dlist(2, space%nh*space%nh))
        rows = 1; cols = 1; kval = 0; mval = 0
        do j = 1, space%nzeta
            do i = 1, space%ntheta
                angle(:, i + space%ntheta*(j - 1)) = [twopi*(i - 1)/space%ntheta, &
                    twopi*(j - 1)/(space%nfp*space%nzeta)]
            end do
        end do
        ! Distinct differences (m2 - m1, n2 - n1) and their phase factors.
        nd = 0
        do h2 = 1, space%nh
            do h1 = 1, space%nh
                pair(h1, h2) = 0
                do d = 1, nd
                    if (dlist(1, d) == space%m(h2) - space%m(h1) .and. &
                        dlist(2, d) == space%n(h2) - space%n(h1)) pair(h1, h2) = d
                end do
                if (pair(h1, h2) == 0) then
                    nd = nd + 1
                    dlist(:, nd) = [space%m(h2) - space%m(h1), &
                        space%n(h2) - space%n(h1)]
                    pair(h1, h2) = nd
                end if
            end do
        end do
        allocate(wave(na, nd))
        do d = 1, nd
            wave(:, d) = exp(cmplx(0, dlist(1, d)*angle(1, :) &
                + dlist(2, d)*angle(2, :), dp))
        end do
        wa = (twopi/space%ntheta)*(twopi/(space%nfp*space%nzeta))
        bad = 0
        !$omp parallel do default(shared) private(e) reduction(+:bad) schedule(dynamic)
        do e = 1, ne
            call element(e, bad)
        end do
        !$omp end parallel do
        if (bad /= 0) return
        if (.not. all(ieee_is_finite(real(kval)))) return
        if (.not. all(ieee_is_finite(aimag(kval)))) return
        do e = 1, ne
            do j = 1, nloc
                do i = 1, j
                    at = (e - 1)*nloc*nloc + (j - 1)*nloc + i
                    mr = (e - 1)*nloc*nloc + (i - 1)*nloc + j
                    kval(at) = (kval(at) + conjg(kval(mr)))/2
                    kval(mr) = conjg(kval(at))
                    mval(at) = (mval(at) + conjg(mval(mr)))/2
                    mval(mr) = conjg(mval(at))
                end do
            end do
        end do
        call csc_from_triplet(space%ndof, space%ndof, rows, cols, kval, &
            stiffness, status)
        if (status%code /= FORTSPARSE_OK) return
        call csc_from_triplet(space%ndof, space%ndof, rows, cols, mval, mass, status)
    contains
        subroutine element(e, bad)
            ! Element-local work with private scratch; writes the disjoint
            ! triplet block of element e.
            integer, intent(in) :: e
            integer, intent(inout) :: bad
            complex(dp) :: fk(12, 12, nd), fm(3, 3, nd), jv(4, nloc), kij, mij
            real(dp) :: wf(12, 12, na), wm(3, 3, na), xq(nq), wq(nq), wf1(12, 12), &
                wm1(3, 3), jac
            real(dp), allocatable :: sv(:), sd(:), tv(:), td(:)
            integer :: dofs(nloc), harm(nloc), comp(nloc), q, a, l, c, k, first, ok
            integer :: ii, jj, ia, ib, dd, base, hh, at
            call gauss_legendre_ab(nq, space%edges(e - 1), space%edges(e), xq, wq)
            base = (e - 1)*nloc*nloc
            do q = 1, nq
                do a = 1, na
                    call flux_point(equilibrium, xq(q), angle(1, a), angle(2, a), &
                        gam, mu, wf1, wm1, jac, ok)
                    if (ok /= 0) then
                        bad = bad + 1; return
                    end if
                    wf(:, :, a) = wa*abs(jac)*wf1
                    wm(:, :, a) = wa*abs(jac)*wm1
                end do
                do dd = 1, nd
                    fk(:, :, dd) = 0; fm(:, :, dd) = 0
                    do a = 1, na
                        fk(:, :, dd) = fk(:, :, dd) + wf(:, :, a)*wave(a, dd)
                        fm(:, :, dd) = fm(:, :, dd) + wm(:, :, a)*wave(a, dd)
                    end do
                end do
                call evaluate_bspline_basis(space%knots, p_, xq(q), sv, sd, ok)
                if (ok /= 0) then
                    bad = bad + 1; return
                end if
                call evaluate_bspline_basis(space%knots(2:size(space%knots) - 1), &
                    p_ - 1, xq(q), tv, td, ok)
                if (ok /= 0) then
                    bad = bad + 1; return
                end if
                first = (e - 1)*p_
                l = 0
                do hh = 1, space%nh
                    do c = 1, 3
                        do k = first + 1, first + p_ + merge(1, 0, c == 1)
                            l = l + 1
                            harm(l) = hh; comp(l) = c
                            if (c == 1) then
                                dofs(l) = space%radial_dof(k, hh)
                                call set_jet(jv(:, l), sv(k), sd(k), hh)
                            else
                                if (c == 2) dofs(l) = space%theta_dof(k, hh)
                                if (c == 3) dofs(l) = space%zeta_dof(k, hh)
                                call set_jet(jv(:, l), tv(k), td(k), hh)
                            end if
                        end do
                    end do
                end do
                do jj = 1, nloc
                    if (dofs(jj) == 0) cycle
                    do ii = 1, nloc
                        if (dofs(ii) == 0) cycle
                        dd = pair(harm(ii), harm(jj))
                        kij = 0
                        do ib = 1, 4
                            do ia = 1, 4
                                kij = kij + conjg(jv(ia, ii)) &
                                    *fk(comp(ii) + 3*(ia - 1), comp(jj) + 3*(ib - 1), dd) &
                                    *jv(ib, jj)
                            end do
                        end do
                        mij = conjg(jv(1, ii))*fm(comp(ii), comp(jj), dd)*jv(1, jj)
                        at = base + (jj - 1)*nloc + ii
                        rows(at) = dofs(ii); cols(at) = dofs(jj)
                        kval(at) = kval(at) + wq(q)*kij
                        mval(at) = mval(at) + wq(q)*mij
                    end do
                end do
            end do
        end subroutine element

        subroutine set_jet(jet, value, slope, h)
            ! Nonzero jet entries of one basis function: D, D_rho, D_theta,
            ! D_zeta of its component, for the harmonic h.
            complex(dp), intent(out) :: jet(4)
            real(dp), intent(in) :: value, slope
            integer, intent(in) :: h
            jet = [cmplx(value, 0, dp), cmplx(slope, 0, dp), &
                cmplx(0, space%m(h)*value, dp), cmplx(0, space%n(h)*value, dp)]
        end subroutine set_jet
    end subroutine flux_assemble
end module milca_ideal_flux
