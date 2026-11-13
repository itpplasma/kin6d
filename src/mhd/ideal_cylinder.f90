module milca_ideal_cylinder
    ! Fixed-boundary ideal-MHD stiffness and mass of a periodic screw pinch for
    ! one Fourier pair exp(i(m theta + kz z)) (docs/stability.md, P1e.2).
    ! Unknowns are the contravariant densities X = r xi_r in the B-spline space
    ! S_p and T = xi_theta, Z = r xi_z in its derivative space S_(p-1), so that
    ! r div xi = X' + i (m T + kz Z) is exact in S_(p-1). Interior knots of
    ! multiplicity p make S_p continuous and S_(p-1) discontinuous: the hybrid
    ! element of the plan.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortnum_quadrature, only: gauss_legendre_ab
    use fortnum_special_bessel, only: bessel_kn, bessel_in
    use milca_real_identity, only: same_real_value
    use fortfem_bspline_feec, only: evaluate_bspline_basis
    use fortsparse, only: csc_z_t, csc_from_triplet, fortsparse_status_t, &
        FORTSPARSE_OK, status_set, FORTSPARSE_INVALID_MATRIX
    use milca_ideal_xi, only: xi_gradient_n, XI_OK
    use milca_ideal_dw, only: dw_bernstein_form
    use milca_ideal_dw_leaves, only: dw_screw_pinch_field, dw_cylinder_density_jets
    use milca_s0_leaves, only: s0_mass_dp
    implicit none
    private
    public :: cylinder_space_t, cylinder_space_init, cylinder_assemble
    public :: cylinder_equilibrium, cylinder_mode

    type :: cylinder_space_t
        integer :: degree = 0, m = 0, nr = 0, ndof = 0
        real(dp) :: kz = 0
        !> Free plasma boundary at the last edge; vacuum out to a conducting
        !> wall at radius wall (<= 0: no wall).
        logical :: free = .false.
        real(dp) :: wall = 0
        real(dp), allocatable :: knots(:), edges(:)
        !> dof of each S_p basis (0: constrained); S_(p-1) dofs follow.
        integer, allocatable :: radial_dof(:), theta_dof(:), axial_dof(:)
    end type cylinder_space_t

    abstract interface
        pure subroutine cylinder_equilibrium(r, bt, btp, bz, bzp, p, pp, rho)
            import :: dp
            real(dp), intent(in) :: r
            real(dp), intent(out) :: bt, btp, bz, bzp, p, pp, rho
        end subroutine cylinder_equilibrium
    end interface
contains
    subroutine cylinder_space_init(space, edges, degree, m, kz, status, wall)
        ! edges(0:n): 0 = axis < ... < plasma edge. X vanishes on the axis, Z on
        ! the axis (finite xi_z), T on the axis unless |m| = 1. Without wall,
        ! X also vanishes at the edge (fixed boundary); with wall present the
        ! edge is free and a vacuum fills the region up to the conducting wall
        ! at radius wall (wall <= 0: vacuum to infinity).
        type(cylinder_space_t), intent(out) :: space
        real(dp), intent(in) :: edges(0:)
        integer, intent(in) :: degree, m
        real(dp), intent(in) :: kz
        integer, intent(out) :: status
        real(dp), intent(in), optional :: wall
        integer :: n, k, at, next
        status = 1
        n = size(edges) - 1
        if (degree < 1 .or. n < 1) return
        if (.not. same_real_value(edges(0), 0.0_dp)) return
        if (any(edges(1:) <= edges(:n-1))) return
        space%degree = degree; space%m = m; space%kz = kz
        space%edges = edges
        if (present(wall)) then
            space%free = .true.; space%wall = wall
            if (wall > 0 .and. wall <= edges(n)) return
            if (m == 0 .and. abs(kz) <= 0) return
        end if
        allocate(space%knots(2*(degree + 1) + (n - 1)*degree))
        space%knots(1:degree + 1) = edges(0)
        at = degree + 1
        do k = 1, n - 1
            space%knots(at + 1:at + degree) = edges(k)
            at = at + degree
        end do
        space%knots(at + 1:at + degree + 1) = edges(n)
        space%nr = size(space%knots) - degree - 1
        allocate(space%radial_dof(space%nr), space%theta_dof(space%nr - 1), &
            space%axial_dof(space%nr - 1))
        next = 0
        space%radial_dof = 0; space%theta_dof = 0; space%axial_dof = 0
        do k = 2, space%nr - merge(0, 1, space%free)
            next = next + 1; space%radial_dof(k) = next
        end do
        do k = 1, space%nr - 1
            if (k == 1 .and. abs(m) /= 1) cycle
            next = next + 1; space%theta_dof(k) = next
        end do
        do k = 2, space%nr - 1
            next = next + 1; space%axial_dof(k) = next
        end do
        space%ndof = next
        status = 0
    end subroutine cylinder_space_init

    subroutine cylinder_assemble(space, equilibrium, gam, mu, stiffness, mass, &
            status)
        ! Weights are r dr; the common 4 pi^2 R0 volume factor cancels. The
        ! Eulerian form needs only div xi and tangential derivatives, both
        ! defined for the discontinuous tangential densities; the Lagrangian
        ! form needs H1 fields and is not used here.
        type(cylinder_space_t), intent(in) :: space
        procedure(cylinder_equilibrium) :: equilibrium
        real(dp), intent(in) :: gam, mu
        type(csc_z_t), intent(out) :: stiffness, mass
        type(fortsparse_status_t), intent(out) :: status
        integer, allocatable :: rows(:), cols(:)
        complex(dp), allocatable :: kval(:), mval(:)
        real(dp), allocatable :: xq(:), wq(:)
        complex(dp) :: v(3, 3*space%degree + 1), a(3, 3, 3*space%degree + 1)
        integer :: dofs(3*space%degree + 1)
        real(dp) :: b(3), g(3, 3), gp(3), bt, btp, bz, bzp, p, pp, rho, w
        integer :: p_, nb, nq, e, q, i, j, at, mr, ntrip, nspan, ok
        call status_set(status, FORTSPARSE_INVALID_MATRIX, 'invalid cylinder input')
        if (space%ndof < 1) return
        if (.not. ieee_is_finite(mu) .or. mu <= 0) return
        if (.not. ieee_is_finite(gam) .or. gam < 0) return
        p_ = space%degree; nb = 3*p_ + 1; nq = p_ + 3
        nspan = size(space%edges) - 1
        ntrip = nb*nb*nspan
        allocate(rows(ntrip), cols(ntrip), kval(ntrip), mval(ntrip), xq(nq), wq(nq))
        rows = 1; cols = 1; kval = 0; mval = 0
        do e = 1, nspan
            call gauss_legendre_ab(nq, space%edges(e - 1), space%edges(e), xq, wq)
            do q = 1, nq
                call span_basis(space, e, xq(q), dofs, v, a, ok)
                if (ok /= 0) return
                call equilibrium(xq(q), bt, btp, bz, bzp, p, pp, rho)
                call dw_screw_pinch_field(xq(q), bt, btp, bz, bzp, b(1), b(2), &
                    b(3), g(1,1), g(2,1), g(3,1), g(1,2), g(2,2), g(3,2), &
                    g(1,3), g(2,3), g(3,3))
                gp = [pp, 0.0_dp, 0.0_dp]
                w = wq(q)*xq(q)
                do j = 1, nb
                    if (dofs(j) == 0) cycle
                    do i = 1, nb
                        if (dofs(i) == 0) cycle
                        at = (e - 1)*nb*nb + (j - 1)*nb + i
                        rows(at) = dofs(i); cols(at) = dofs(j)
                        kval(at) = kval(at) + w*dw_bernstein_form(b, g, p, gp, &
                            gam, mu, v(:, j), a(:, :, j), v(:, i), a(:, :, i))
                        mval(at) = mval(at) + w*full_mass(rho, v(:, j), v(:, i))
                    end do
                end do
            end do
        end do
        if (space%free) then
            ! Vacuum: 2 delta W_V/(4 pi^2 R0) = -(a/mu) F_a^2 G |xi_r(a)|^2 with
            ! G = phi/phi' at the edge, xi_r = X/a, on the edge S_p function.
            call equilibrium(space%edges(nspan), bt, btp, bz, bzp, p, pp, rho)
            w = vacuum_ratio(space%m, abs(space%kz), space%edges(nspan), space%wall)
            w = -space%edges(nspan)/mu*(space%m*bt/space%edges(nspan) &
                + space%kz*bz)**2*w/space%edges(nspan)**2
            at = (nspan - 1)*nb*nb + p_*nb + p_ + 1
            i = space%radial_dof(space%nr)
            rows(at) = i; cols(at) = i
            kval(at) = kval(at) + w
        end if
        if (.not. all(ieee_is_finite(real(kval)))) return
        if (.not. all(ieee_is_finite(aimag(kval)))) return
        ! Exact Hermitian symmetry for the native LDL^H contract.
        do e = 1, nspan
            do j = 1, nb
                do i = 1, j
                    at = (e - 1)*nb*nb + (j - 1)*nb + i
                    mr = (e - 1)*nb*nb + (i - 1)*nb + j
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
    end subroutine cylinder_assemble

    subroutine cylinder_mode(space, coefficients, r, xi, status)
        ! Physical (xi_r, xi_theta, xi_z) amplitudes of a coefficient vector.
        type(cylinder_space_t), intent(in) :: space
        complex(dp), intent(in) :: coefficients(:)
        real(dp), intent(in) :: r
        complex(dp), intent(out) :: xi(3)
        integer, intent(out) :: status
        complex(dp) :: v(3, 3*space%degree + 1), a(3, 3, 3*space%degree + 1)
        integer :: dofs(3*space%degree + 1), e, i
        xi = 0; status = 1
        if (size(coefficients) /= space%ndof) return
        if (r < space%edges(0) .or. r > space%edges(size(space%edges) - 1)) return
        e = max(1, min(size(space%edges) - 1, count(space%edges(1:) <= r) + 1))
        if (r <= 0) return
        call span_basis(space, e, r, dofs, v, a, status)
        if (status /= 0) return
        do i = 1, size(dofs)
            if (dofs(i) /= 0) xi = xi + coefficients(dofs(i))*v(:, i)
        end do
    end subroutine cylinder_mode

    subroutine span_basis(space, e, r, dofs, v, a, status)
        ! Local fields on span e: S_p basis k = (e-1)p+1 .. ep+1 in xi_r and
        ! S_(p-1) basis (e-1)p+1 .. ep in xi_theta and xi_z.
        type(cylinder_space_t), intent(in) :: space
        integer, intent(in) :: e
        real(dp), intent(in) :: r
        integer, intent(out) :: dofs(:), status
        complex(dp), intent(out) :: v(:, :), a(:, :, :)
        real(dp), allocatable :: sv(:), sd(:), tv(:), td(:)
        complex(dp) :: dz(3), dr(3)
        integer :: p_, first, k, l, c, gstat
        p_ = space%degree
        v = 0; a = 0; dofs = 0
        call evaluate_bspline_basis(space%knots, p_, r, sv, sd, status)
        if (status /= 0) return
        call evaluate_bspline_basis(space%knots(2:size(space%knots) - 1), p_ - 1, &
            r, tv, td, status)
        if (status /= 0) return
        first = (e - 1)*p_
        l = 0
        do k = first + 1, first + p_ + 1
            l = l + 1
            dofs(l) = space%radial_dof(k)
            call density(l, sv(k), sd(k), 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp)
        end do
        do c = 2, 3
            do k = first + 1, first + p_
                l = l + 1
                if (c == 2) then
                    dofs(l) = space%theta_dof(k)
                    call density(l, 0.0_dp, 0.0_dp, tv(k), td(k), 0.0_dp, 0.0_dp)
                else
                    dofs(l) = space%axial_dof(k)
                    call density(l, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, tv(k), td(k))
                end if
            end do
        end do
    contains
        subroutine density(l, xd, xdp, td_, tdp, zd, zdp)
            integer, intent(in) :: l
            real(dp), intent(in) :: xd, xdp, td_, tdp, zd, zdp
            real(dp) :: f(6)
            call dw_cylinder_density_jets(r, xd, xdp, td_, tdp, zd, zdp, f(1), &
                f(2), f(3), f(4), f(5), f(6))
            v(:, l) = cmplx(f(1:5:2), 0, dp)
            dr = cmplx(f(2:6:2), 0, dp)
            dz = cmplx(0, space%kz, dp)*v(:, l)
            call xi_gradient_n(space%m, r, v(:, l), dr, dz, a(:, :, l), gstat)
            if (gstat /= XI_OK) status = 1
        end subroutine density
    end subroutine span_basis

    real(dp) function vacuum_ratio(m, k, a, b) result(g)
        ! phi/phi' at r = a for the vacuum potential with phi' = 0 at r = b
        ! (b <= 0: decaying to infinity); exp(i(m theta + k z)), k >= 0.
        integer, intent(in) :: m
        real(dp), intent(in) :: k, a, b
        real(dp) :: x, xb, c, ratio
        integer :: mm
        mm = abs(m)
        if (k <= 0) then
            ratio = 0
            if (b > 0) ratio = (a/b)**(2*mm)
            g = -a*(1 + ratio)/(mm*(1 - ratio))
            return
        end if
        x = k*a
        c = 0
        if (b > 0) then
            xb = k*b
            c = kprime(mm, xb)/iprime(mm, xb)
        end if
        g = (bessel_kn(mm, x) - c*bessel_in(mm, x)) &
            /(k*(kprime(mm, x) - c*iprime(mm, x)))
    contains
        real(dp) function kprime(n, y)
            integer, intent(in) :: n
            real(dp), intent(in) :: y
            kprime = -(bessel_kn(abs(n - 1), y) + bessel_kn(n + 1, y))/2
        end function kprime
        real(dp) function iprime(n, y)
            integer, intent(in) :: n
            real(dp), intent(in) :: y
            iprime = (bessel_in(abs(n - 1), y) + bessel_in(n + 1, y))/2
        end function iprime
    end function vacuum_ratio

    pure complex(dp) function full_mass(rho, x, y) result(value)
        real(dp), intent(in) :: rho
        complex(dp), intent(in) :: x(3), y(3)
        real(dp) :: rr, ii, ir, ri
        call s0_mass_dp(rho, real(x(1)), real(x(2)), real(x(3)), real(y(1)), &
            real(y(2)), real(y(3)), rr)
        call s0_mass_dp(rho, aimag(x(1)), aimag(x(2)), aimag(x(3)), aimag(y(1)), &
            aimag(y(2)), aimag(y(3)), ii)
        call s0_mass_dp(rho, aimag(x(1)), aimag(x(2)), aimag(x(3)), real(y(1)), &
            real(y(2)), real(y(3)), ir)
        call s0_mass_dp(rho, real(x(1)), real(x(2)), real(x(3)), aimag(y(1)), &
            aimag(y(2)), aimag(y(3)), ri)
        value = cmplx(rr + ii, ir - ri, dp)
    end function full_mass
end module milca_ideal_cylinder
