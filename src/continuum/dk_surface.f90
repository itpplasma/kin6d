!> Transport matrix of a 3D surface with the energy-coupled linearised
!! collision operator in real arithmetic (the n = 0 problem has no complex
!! part), the faster solver of the operator of ntv_surface_coefficients: same
!! drives, responses and units (see there), same collision rows
!! (ntv_collision_rows), same Legendre x Fourier surface discretisation, with
!! the density null space of each speed node bordered as in dk_mono (a
!! density source and <w h_0> = 0 per node instead of the unweighted mean).
!!
!! Unknowns h(i, l, m): surface point i, Legendre degree l, speed-species
!! index m = (a - 1) n_x + k. After the row scaling the operator is
!! S x R - I x C_l (S streaming and mirror force per unit speed, R the
!! streaming coefficient of each node, C_l the collision rows of degree l);
!! S acts on all m at once as two dense surface matrices per degree.
!! FGMRES is right-preconditioned by one dk_mono factorisation per eigenmode
!! of R^-1 C_4, with the diagonal of V^-1 R^-1 C_l V as the rate of degree l:
!! exact for pitch-angle scattering. Memory is that of n_x n_s monoenergetic
!! surface problems; the modes are factorised and solved in parallel.
module milca_dk_surface
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_krylov, only: real_fgmres_operator, KRYLOV_OK
    use milca_dk_mono, only: mono_factor_t, mono_surface_point, mono_factor_rates, &
        mono_solve, mono_coupling, MONO_OK
    use milca_dk_ntv, only: ntv_collision_rows, NTV_OK
    use milca_dk_collop, only: collop_species_t, collop_grid_t
    implicit none
    private
    public :: surface_transport
    integer, parameter, public :: SURFACE_OK = 0, SURFACE_BAD_INPUT = 1, &
        SURFACE_SINGULAR = 2, SURFACE_NO_CONVERGENCE = 3

    real(dp), parameter :: pi = 4*atan(1.0_dp)
    ! Relative residual of the FGMRES solves.
    real(dp), parameter :: tolerance = 1e-10_dp
    integer, parameter :: restart = 100, max_iterations = 5000

    interface
        subroutine dgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
            import :: dp
            character, intent(in) :: transa, transb
            integer, intent(in) :: m, n, k, lda, ldb, ldc
            real(dp), intent(in) :: alpha, beta, a(lda, *), b(ldb, *)
            real(dp), intent(inout) :: c(ldc, *)
        end subroutine dgemm
        subroutine dgeev(jobvl, jobvr, n, a, lda, wr, wi, vl, ldvl, vr, ldvr, work, &
                lwork, info)
            import :: dp
            character, intent(in) :: jobvl, jobvr
            integer, intent(in) :: n, lda, ldvl, ldvr, lwork
            real(dp), intent(inout) :: a(lda, *)
            real(dp), intent(out) :: wr(*), wi(*), vl(ldvl, *), vr(ldvr, *), work(*)
            integer, intent(out) :: info
        end subroutine dgeev
        subroutine dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
            import :: dp
            integer, intent(in) :: n, nrhs, lda, ldb
            real(dp), intent(inout) :: a(lda, *), b(ldb, *)
            integer, intent(out) :: ipiv(*), info
        end subroutine dgesv
    end interface

contains

    !> Arguments and d as in ntv_surface_coefficients; iterations the FGMRES
    !> iterations summed over the 3 n_s drives.
    subroutine surface_transport(b, nper, q, i_over_g, kappa, species, lnlambda, grid, &
            model, n_leg, d, iterations, status)
        real(dp), intent(in) :: b(:, :), q, i_over_g, kappa
        integer, intent(in) :: nper, n_leg
        type(collop_species_t), intent(in) :: species(:)
        real(dp), intent(in) :: lnlambda(:, :)
        type(collop_grid_t), intent(in) :: grid
        integer, intent(in) :: model
        real(dp), intent(out) :: d(:, :, :, :)
        integer, intent(out) :: iterations, status
        type(mono_factor_t) :: geo
        type(mono_factor_t), allocatable :: modes(:)
        real(dp), allocatable :: rows(:, :, :, :, :), ratio(:, :), d3(:), d4(:), d6(:), &
            vt(:), tz(:), coll(:, :, :), lo(:, :, :), up(:, :, :), v(:, :), vu(:, :), &
            rinv(:), rhs(:, :, :), x(:, :, :), radial(:), measure(:), bv(:)
        real(dp) :: j0
        integer :: n, lmax, ns, nx, nm, a, sb, jf, k, l, st, its

        status = SURFACE_BAD_INPUT
        iterations = 0
        ns = size(species)
        if (any(shape(d) /= [ns, ns, 3, 3]) .or. mod(n_leg, 2) /= 0) return
        call mono_surface_point(b, nper, q, i_over_g, kappa, n_leg, geo, st)
        if (st /= MONO_OK) return
        call ntv_collision_rows(species, lnlambda, grid, model, n_leg, rows, ratio, d3, &
            d4, d6, vt, tz, st)
        if (st /= NTV_OK) return
        n = geo%n
        lmax = geo%lmax
        nx = size(ratio, 1)
        nm = nx*ns
        ! Collision rows over (node, species), node fastest.
        allocate (coll(nm, nm, 0:lmax), lo(n, n, 0:lmax), up(n, n, 0:lmax))
        do l = 0, lmax
            do sb = 1, ns
                do a = 1, ns
                    do k = 1, nx
                        coll((a - 1)*nx + k, (sb - 1)*nx + 1:sb*nx, l) = rows(k, :, a, sb, l)
                    end do
                end do
            end do
            call mono_coupling(geo, l, lo(:, :, l), up(:, :, l))
        end do
        rinv = 1/reshape(ratio, [nm])
        status = SURFACE_SINGULAR
        call modal_factors(st)
        if (st /= SURFACE_OK) return

        bv = reshape(b, [n])
        ! Radial drive (1 + xi^2) S/b, S = b_drift/q, in terms of the dk_mono
        ! drift: P_0 and P_2 components -2 drift0/(kappa q), -2 drift2/(kappa q).
        radial = -2/(kappa*q)*geo%p%drift0
        measure = 2/bv**3
        j0 = sum(1/bv**2)
        allocate (rhs(n, 0:lmax, nm))
        do sb = 1, ns
            do jf = 1, 3
                rhs = 0
                do k = 1, nx
                    associate (m => (sb - 1)*nx + k, row => ratio(k, sb)/(vt(sb)*d3(k)))
                        if (jf < 3) then
                            rhs(:, 0, m) = row*tz(sb)*merge(d4(k), d6(k), jf == 1)*radial
                            rhs(:, 2, m) = row*tz(sb)*merge(d4(k), d6(k), jf == 1) &
                                *(-2/(kappa*q))*geo%p%drift2
                        else
                            rhs(:, 1, m) = row*vt(sb)*d3(k)*geo%p%drive_e
                        end if
                    end associate
                end do
                call solve(rhs, x, its, st)
                iterations = iterations + its
                if (st /= SURFACE_OK) then
                    status = st
                    return
                end if
                do a = 1, ns
                    d(a, sb, :, jf) = 0
                    do k = 1, nx
                        associate (m => (a - 1)*nx + k)
                            ! flux: 2 pi/n sum measure S (8/3 h_0 + 4/15 h_2);
                            ! flow: 2 pi/n sum measure b^2 (2/3) h_1; over 2 J_0.
                            d(a, sb, 1, jf) = d(a, sb, 1, jf) + d4(k)*flux(x(:, :, m))
                            d(a, sb, 2, jf) = d(a, sb, 2, jf) + d6(k)*flux(x(:, :, m))
                            d(a, sb, 3, jf) = d(a, sb, 3, jf) + d3(k) &
                                *sum(measure*bv**2*2/3.0_dp*x(:, 1, m))
                        end associate
                    end do
                    d(a, sb, 1:2, jf) = tz(a)/(sqrt(pi)*2*j0)*d(a, sb, 1:2, jf)
                    d(a, sb, 3, jf) = vt(a)/(sqrt(pi)*2*j0)*d(a, sb, 3, jf)
                end do
            end do
        end do
        status = SURFACE_OK
    contains
        !> sum measure S (8/3 h_0 + 4/15 h_2), S = b_drift/q = -3 b drift0/(2 kappa q).
        real(dp) function flux(h)
            real(dp), intent(in) :: h(:, 0:)
            flux = sum(measure*(-1.5_dp*bv*geo%p%drift0/(kappa*q)) &
                *(8/3.0_dp*h(:, 0) + 4/15.0_dp*h(:, 2)))
        end function flux

        !> V, V^-1 and one dk_mono factorisation per eigenmode of R^-1 C_4.
        subroutine modal_factors(status)
            integer, intent(out) :: status
            real(dp), allocatable :: w(:, :), wr(:), wi(:), work(:), vl(:, :), rate(:, :)
            integer, allocatable :: piv(:), mode_status(:)
            integer :: m, info
            status = SURFACE_SINGULAR
            allocate (w(nm, nm), wr(nm), wi(nm), work(8*nm), vl(1, 1), v(nm, nm), &
                vu(nm, nm), piv(nm), rate(0:lmax, nm), mode_status(nm))
            ! Eigenvectors of degree 4: fewest FGMRES iterations among degrees 1
            ! to 9 on the DIII-D helical core (2402 for degree 1, 857 for 4).
            do m = 1, nm
                w(m, :) = rinv(m)*coll(m, :, min(4, lmax))
            end do
            call dgeev('N', 'V', nm, w, nm, wr, wi, vl, 1, v, nm, work, size(work), info)
            if (info /= 0) return
            if (any(abs(wi) > 1e-10_dp*maxval(abs(wr)))) return
            w = v
            vu = 0
            do m = 1, nm
                vu(m, m) = 1
            end do
            call dgesv(nm, nm, w, nm, piv, vu, nm, info)
            if (info /= 0) return
            do l = 0, lmax
                do m = 1, nm
                    w(m, :) = rinv(m)*coll(m, :, l)
                end do
                w = matmul(vu, matmul(w, v))
                do m = 1, nm
                    rate(l, m) = -w(m, m)
                end do
            end do
            allocate (modes(nm))
            !$omp parallel do schedule(dynamic, 1) default(shared) private(m)
            do m = 1, nm
                modes(m) = geo
                call mono_factor_rates(modes(m), rate(:, m), mode_status(m))
            end do
            !$omp end parallel do
            if (any(mode_status /= MONO_OK)) return
            status = SURFACE_OK
        end subroutine modal_factors

        !> FGMRES on the bordered operator: unknowns h(n, 0:lmax, nm) and the
        !> density sources mu(nm).
        subroutine solve(rhs, x, its, status)
            real(dp), intent(in) :: rhs(:, 0:, :)
            real(dp), allocatable, intent(out) :: x(:, :, :)
            integer, intent(out) :: its, status
            real(dp), allocatable :: bvec(:), u(:)
            real(dp) :: scale, residual
            integer :: ntot, info
            ntot = n*(lmax + 1)*nm + nm
            allocate (bvec(ntot), u(ntot))
            bvec = 0
            bvec(1:ntot - nm) = reshape(rhs, [ntot - nm])
            scale = norm2(bvec)
            its = 0
            allocate (x(n, 0:lmax, nm))
            x = 0
            status = SURFACE_OK
            if (scale == 0) return
            bvec = bvec/scale
            u = 0
            call real_fgmres_operator(apply, precondition, bvec, u, tolerance, &
                max_iterations, restart, info, its, residual)
            status = SURFACE_NO_CONVERGENCE
            if (info /= KRYLOV_OK) return
            x = scale*reshape(u(1:ntot - nm), [n, lmax + 1, nm])
            status = SURFACE_OK
        end subroutine solve

        !> y = A x with the border: h-rows plus mu at l = 0, then <w h_0> rows.
        subroutine apply(input, output)
            real(dp), intent(in) :: input(:)
            real(dp), intent(out) :: output(:)
            real(dp), allocatable :: h(:, :, :), y(:, :, :)
            integer :: l, m, nh
            nh = n*(lmax + 1)*nm
            h = reshape(input(1:nh), [n, lmax + 1, nm])
            allocate (y(n, 0:lmax, nm))
            !$omp parallel do schedule(static) default(shared) private(l, m)
            do l = 0, lmax
                y(:, l, :) = 0
                if (l > 0) call dgemm('N', 'N', n, nm, n, 1.0_dp, lo(:, :, l), n, &
                    h(:, l, :), n, 1.0_dp, y(:, l, :), n)
                if (l < lmax) call dgemm('N', 'N', n, nm, n, 1.0_dp, up(:, :, l), n, &
                    h(:, l + 2, :), n, 1.0_dp, y(:, l, :), n)
                do m = 1, nm
                    y(:, l, m) = ratio_flat(m)*y(:, l, m)
                end do
                call dgemm('N', 'T', n, nm, nm, -1.0_dp, h(:, l + 1, :), n, coll(:, :, l), &
                    nm, 1.0_dp, y(:, l, :), n)
            end do
            !$omp end parallel do
            do m = 1, nm
                y(:, 0, m) = y(:, 0, m) + input(nh + m)
                output(nh + m) = sum(geo%p%weight*h(:, 1, m))
            end do
            output(1:nh) = reshape(y, [nh])
        end subroutine apply

        pure real(dp) function ratio_flat(m)
            integer, intent(in) :: m
            ratio_flat = 1/rinv(m)
        end function ratio_flat

        !> x = P^-1 input: modal coordinates V^-1 R^-1 for the rows and V^-1
        !> for the constraints, one bordered dk_mono solve per mode, back by V
        !> (h) and R V (mu).
        subroutine precondition(input, output)
            real(dp), intent(in) :: input(:)
            real(dp), intent(out) :: output(:)
            real(dp), allocatable :: r(:, :, :), z(:, :, :), c(:), cm(:), mu(:)
            integer :: l, m, nh
            nh = n*(lmax + 1)*nm
            r = reshape(input(1:nh), [n, lmax + 1, nm])
            allocate (z(n, lmax + 1, nm), mu(nm))
            do l = 1, lmax + 1
                do m = 1, nm
                    r(:, l, m) = rinv(m)*r(:, l, m)
                end do
                z(:, l, :) = matmul(r(:, l, :), transpose(vu))
            end do
            c = input(nh + 1:nh + nm)
            cm = matmul(vu, c)
            !$omp parallel do schedule(dynamic, 1) default(shared) private(m)
            do m = 1, nm
                block
                    real(dp), allocatable :: xm(:, :, :)
                    real(dp) :: mum(1)
                    integer :: st
                    call mono_solve(modes(m), reshape(z(:, :, m), [n, lmax + 1, 1]), &
                        [cm(m)], xm, mum, st)
                    z(:, :, m) = xm(:, :, 1)
                    mu(m) = mum(1)
                end block
            end do
            !$omp end parallel do
            do l = 1, lmax + 1
                r(:, l, :) = matmul(z(:, l, :), transpose(v))
            end do
            output(1:nh) = reshape(r, [nh])
            output(nh + 1:nh + nm) = matmul(v, mu)/rinv
        end subroutine precondition
    end subroutine surface_transport
end module milca_dk_surface
