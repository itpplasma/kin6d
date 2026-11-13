program test_ideal_cylinder
    ! Rungs R1 and R2 of docs/stability.md (P1e.2) for the radial hybrid space.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortsparse, only: csc_z_t, fortsparse_status_t, FORTSPARSE_OK, &
        FORTSPARSE_INERTIA_RELIABLE
    use fortnum_eigen, only: eigen_options_t, EIGEN_OK, EIGEN_NEAREST
    use milca_ideal_cylinder, only: cylinder_space_t, cylinder_space_init, &
        cylinder_assemble, cylinder_equilibrium
    use milca_s0_spectrum, only: s0_pencil_t, s0_report_t, s0_prepare, s0_solve
    implicit none
    interface
        subroutine zhegv(itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, &
                rwork, info)
            import :: dp
            integer, intent(in) :: itype, n, lda, ldb, lwork
            character, intent(in) :: jobz, uplo
            complex(dp), intent(inout) :: a(lda, *), b(ldb, *)
            real(dp), intent(out) :: w(*), rwork(*)
            complex(dp), intent(inout) :: work(*)
            integer, intent(out) :: info
        end subroutine zhegv
    end interface
    real(dp), parameter :: kink_reference = -0.1949752661_dp
    real(dp) :: ext_bz = 5
    character(32) :: mode
    mode = 'all'
    if (command_argument_count() > 0) call get_command_argument(1, mode)
    if (mode == 'all' .or. mode == 'uniform') call uniform_cylinder()
    if (mode == 'all' .or. mode == 'kink') call fixed_wall_kink()
    if (mode == 'all' .or. mode == 'limit') call unstable_limit()
    if (mode == 'all' .or. mode == 'free') call external_kink()
    print '(a)', 'PASS ideal cylinder ('//trim(mode)//')'
contains
    subroutine uniform_cylinder()
        ! Uniform Bz = 1, p = 0.3, gamma = 5/3, rho = mu = 1, a = 1, m = 2, kz = 1.
        ! Exact spectrum: the Alfven point kz^2 vA^2 = 1, and fast and slow
        ! roots of omega^4 - omega^2 k^2 (vA^2+cs^2) + k^2 kz^2 vA^2 cs^2 = 0
        ! with k^2 = kz^2 + j^2, J_m'(j) = 0 (xi_r ~ J_m'(j r) vanishes at a).
        real(dp), allocatable :: lam(:)
        real(dp) :: zeros(4), fast(4), slow(4), k2, s, d, err
        integer :: nu, nlow
        call bessel_prime_zeros(2, zeros)
        do nu = 1, 4
            k2 = 1 + zeros(nu)**2
            s = k2*(1 + 0.5_dp); d = sqrt(s**2 - 4*k2*0.5_dp)
            fast(nu) = (s + d)/2; slow(nu) = 2*k2*0.5_dp/(s + d)
        end do
        call spectrum(uniform_field, 1.0_dp, 2, 1.0_dp, 5.0_dp/3, 24, 3, lam)
        do nu = 1, 3
            err = minval(abs(lam - fast(nu)))/fast(nu)
            print '(a,i1,2es14.6,es10.2)', 'fast ', nu, fast(nu), &
                lam(minloc(abs(lam - fast(nu)), 1)), err
            if (err > 1e-4_dp) error stop 'fast magnetosonic eigenvalue'
            err = minval(abs(lam - slow(nu)))/slow(nu)
            print '(a,i1,2es14.6,es10.2)', 'slow ', nu, slow(nu), &
                lam(minloc(abs(lam - slow(nu)), 1)), err
            if (err > 1e-4_dp) error stop 'slow magnetosonic eigenvalue'
        end do
        ! O5: the gaps (slow_1, 1) and (1, fast_1) of the exact spectrum.
        nlow = count(lam > slow(1)*(1 + 1e-4_dp) .and. lam < 1 - 1e-6_dp)
        nlow = nlow + count(lam > 1 + 1e-6_dp .and. lam < fast(1)*(1 - 1e-4_dp))
        print '(a,i0)', 'eigenvalues inside exact gaps: ', nlow
        if (nlow /= 0) error stop 'spectral pollution'
        if (any(lam < slow(size(slow))*(1 - 1e-2_dp) .and. lam < 1.0_dp/3)) &
            error stop 'eigenvalue below the slow accumulation point'
    end subroutine uniform_cylinder

    subroutine fixed_wall_kink()
        ! docs/s0-spectrum.md (S2.1-3): a = 2, Btheta = r, Bz = sqrt(1/4 + r^2),
        ! p = 3/2 (a^2 - r^2), m = 1, kz = -0.4, gamma = 0. The independent P2
        ! oracle gives omega^2 = -0.1949752661 (converged to 1e-9).
        real(dp), allocatable :: lam(:)
        real(dp) :: values(3, 3)
        integer :: level, deg
        do deg = 2, 3
            do level = 1, 3
                call spectrum(kink_field, 2.0_dp, 1, -0.4_dp, 0.0_dp, &
                    4*2**level, deg, lam)
                values(level, deg) = lam(1)
                print '(a,i1,a,i3,es20.11,es10.2)', 'kink p=', deg, ' n=', &
                    4*2**level, lam(1), abs(lam(1)/kink_reference - 1)
            end do
        end do
        if (abs(values(3, 3)/kink_reference - 1) > 1e-6_dp) &
            error stop 'fixed-wall kink differs from the independent oracle'
        if (count(lam < -1e-4_dp) /= 1) error stop 'unstable count (O1)'
        ! Compression adds gamma p |div xi|^2 >= 0 and can only stabilize.
        call spectrum(kink_field, 2.0_dp, 1, -0.4_dp, 5.0_dp/3, 32, 3, lam)
        print '(a,es20.11)', 'kink gamma = 5/3        ', lam(1)
        if (lam(1) < values(3, 3) - 1e-12_dp) error stop 'compression destabilized'
        call production_inertia()
    end subroutine fixed_wall_kink

    subroutine production_inertia()
        ! Sparse LDL^H inertia of A + delta M (O1) and shift-invert Lanczos.
        type(cylinder_space_t) :: space
        type(csc_z_t) :: stiffness, mass
        type(fortsparse_status_t) :: status
        type(s0_pencil_t) :: pencil
        type(s0_report_t) :: report
        type(eigen_options_t) :: options
        complex(dp), allocatable :: start(:), vectors(:, :), omega2(:)
        real(dp), allocatable :: residuals(:)
        integer :: ok, i
        call cylinder_space_init(space, edges(32, 2.0_dp), 3, 1, -0.4_dp, ok)
        if (ok /= 0) error stop 'space'
        call cylinder_assemble(space, kink_field, 0.0_dp, 1.0_dp, stiffness, mass, &
            status)
        if (status%code /= FORTSPARSE_OK) error stop 'assembly'
        call s0_prepare(pencil, stiffness, mass, cmplx(-1e-4_dp, 0, dp), report)
        if (report%status /= EIGEN_OK) error stop 'shifted factor'
        print '(a,3i8)', 'inertia below/at/above -1e-4:', report%below, report%at, &
            report%above
        if (report%inertia_status /= FORTSPARSE_INERTIA_RELIABLE) &
            error stop 'inertia reliability'
        if (report%below /= 1) error stop 'sparse unstable count (O1)'
        call s0_prepare(pencil, stiffness, mass, cmplx(-0.25_dp, 0, dp), report)
        allocate(start(space%ndof), vectors(space%ndof, 1), omega2(1), residuals(1))
        do i = 1, space%ndof
            start(i) = cmplx(sin(1.3_dp*i), cos(0.7_dp*i), dp)
        end do
        options%max_basis = 20; options%max_restarts = 40
        options%which = EIGEN_NEAREST
        call s0_solve(pencil, start, 1e-11_dp, options, omega2, vectors, &
            residuals, report)
        if (report%status /= EIGEN_OK) error stop 'shift-invert Lanczos'
        print '(a,es20.11)', 'kink shift-invert       ', real(omega2(1))
        if (abs(real(omega2(1))/kink_reference - 1) > 1e-6_dp) &
            error stop 'sparse eigenvalue differs'
    end subroutine production_inertia

    subroutine unstable_limit()
        ! Btheta = r, Bz -> 0, |kz| -> infinity: omega^2 -> -3 at order kz^-2
        ! (docs/s0-spectrum.md, S2.1-5); a = 2, p = a^2 - r^2, gamma = 0.
        real(dp), allocatable :: lam(:)
        real(dp) :: values(3), extrapolated
        real(dp), parameter :: kz(3) = [-40.0_dp, -80.0_dp, -160.0_dp]
        integer :: i
        do i = 1, 3
            call spectrum(pinch_field, 2.0_dp, 1, kz(i), 0.0_dp, 64, 4, lam)
            values(i) = lam(1)
            print '(a,f8.1,es20.11)', 'limit kz=', kz(i), values(i)
        end do
        extrapolated = (4*values(3) - values(2))/3
        print '(a,es20.11)', 'extrapolated            ', extrapolated
        if (abs(extrapolated/(-3) - 1) > 1e-4_dp) error stop 'analytic limit -3'
    end subroutine unstable_limit

    subroutine external_kink()
        ! docs/ideal-free-boundary.md: B = (-y, x, Bz), p = 1 - r^2, a = mu = rho
        ! = 1, m = 1, kz = -0.1, q = -kz Bz = 0.5. The rigid displacement
        ! (X = r, T = i, Z = 0) has omega^2 = q^2 - 1 + f (q - 1)^2 with the
        ! vacuum factor f = K_1/(K_1 + ka K_0); for ka = 0.1, -0.50601.
        use fortnum_special_bessel, only: bessel_kn
        type(cylinder_space_t) :: space
        type(csc_z_t) :: stiffness, mass
        type(fortsparse_status_t) :: status
        complex(dp), allocatable :: v(:), kv(:), mv(:)
        real(dp), allocatable :: lam(:), greville(:)
        real(dp) :: f, exact, rq, previous
        integer :: ok, k, j
        real(dp), parameter :: walls(4) = [3.0_dp, 1.6_dp, 1.3_dp, 1.1_dp]
        real(dp), parameter :: marginal(2) = [1.40_dp, 1.43_dp]
        f = bessel_kn(1, 0.1_dp)/(bessel_kn(1, 0.1_dp) + 0.1_dp*bessel_kn(0, 0.1_dp))
        exact = 0.25_dp - 1 + f*0.25_dp
        call cylinder_space_init(space, edges(16, 1.0_dp), 3, 1, -0.1_dp, ok, &
            wall=0.0_dp)
        if (ok /= 0) error stop 'free-boundary space'
        call cylinder_assemble(space, external_field, 0.0_dp, 1.0_dp, stiffness, &
            mass, status)
        if (status%code /= FORTSPARSE_OK) error stop 'free-boundary assembly'
        allocate(v(space%ndof)); v = 0
        greville = [(sum(space%knots(k + 1:k + space%degree))/space%degree, &
            k = 1, space%nr)]
        do k = 1, space%nr
            if (space%radial_dof(k) > 0) v(space%radial_dof(k)) = greville(k)
        end do
        do k = 1, space%nr - 1
            if (space%theta_dof(k) > 0) v(space%theta_dof(k)) = (0.0_dp, 1.0_dp)
        end do
        kv = matvec(stiffness, v); mv = matvec(mass, v)
        rq = real(dot_product(v, kv))/real(dot_product(v, mv))
        print '(a,2es18.9,es10.2)', 'rigid kink quotient, exact', rq, exact, &
            abs(rq/exact - 1)
        if (abs(rq/exact - 1) > 1e-10_dp) error stop 'rigid external kink energy'
        if (abs(exact + 0.50601_dp) > 1e-5_dp) error stop 'vacuum factor'
        call spectrum_space(space, stiffness, mass, lam)
        print '(a,es18.9)', 'lowest eigenvalue, no wall', lam(1)
        if (lam(1) > rq + 1e-12_dp) error stop 'Ritz bound violated'
        ! Compression: the rigid displacement is divergence free, so it still
        ! bounds the eigenvalue, which the adiabatic term can only raise.
        block
            real(dp), allocatable :: adiabatic(:)
            call cylinder_assemble(space, external_field, 5.0_dp/3, 1.0_dp, &
                stiffness, mass, status)
            if (status%code /= FORTSPARSE_OK) error stop 'adiabatic assembly'
            call spectrum_space(space, stiffness, mass, adiabatic)
            print '(a,es18.9)', 'lowest eigenvalue, no wall, gamma = 5/3', &
                adiabatic(1)
            if (adiabatic(1) > rq + 1e-12_dp) error stop 'adiabatic Ritz bound'
            if (adiabatic(1) < lam(1)) error stop 'compression destabilized'
        end block
        ! Long-wavelength limit ka -> 0: omega^2 -> 2 q (q - 1) = -0.5 (O2).
        ext_bz = 500
        call cylinder_space_init(space, edges(16, 1.0_dp), 3, 1, -1e-3_dp, ok, &
            wall=0.0_dp)
        call cylinder_assemble(space, external_field, 0.0_dp, 1.0_dp, stiffness, &
            mass, status)
        ext_bz = 5
        block
            real(dp), allocatable :: low(:)
            call spectrum_space(space, stiffness, mass, low)
            print '(a,es18.9)', 'lowest eigenvalue, ka = 1e-3', low(1)
            if (abs(low(1)/(-0.5_dp) - 1) > 1e-4_dp) error stop 'kink limit 2q(q-1)'
        end block
        if (count(lam < -1e-4_dp) /= 1) error stop 'external kink count (O1)'
        previous = lam(1)
        do j = 1, size(walls)
            call cylinder_space_init(space, edges(16, 1.0_dp), 3, 1, -0.1_dp, ok, &
                wall=walls(j))
            call cylinder_assemble(space, external_field, 0.0_dp, 1.0_dp, &
                stiffness, mass, status)
            call spectrum_space(space, stiffness, mass, lam)
            print '(a,f5.2,es18.9)', 'lowest eigenvalue, wall at', walls(j), lam(1)
            if (lam(1) < previous) error stop 'a closer wall destabilized'
            previous = lam(1)
        end do
        if (lam(1) < -1e-4_dp) error stop 'close wall did not stabilize'
        ! Marginal wall at ka -> 0: q^2 - 1 + (q - 1)^2 (1 + c)/(1 - c) = 0 with
        ! c = (a/b)^2 gives b/a = sqrt(2) for q = 0.5 (O3).
        ext_bz = 500
        do j = 1, 2
            call cylinder_space_init(space, edges(16, 1.0_dp), 3, 1, -1e-3_dp, ok, &
                wall=marginal(j))
            call cylinder_assemble(space, external_field, 0.0_dp, 1.0_dp, &
                stiffness, mass, status)
            call spectrum_space(space, stiffness, mass, lam)
            print '(a,f5.2,es18.9)', 'ka = 1e-3, wall at', marginal(j), lam(1)
            if (j == 1 .and. lam(1) < -1e-4_dp) error stop 'wall inside marginal'
            if (j == 2 .and. lam(1) > -1e-4_dp) error stop 'wall outside marginal'
        end do
        ext_bz = 5
    end subroutine external_kink

    pure subroutine external_field(r, bt, btp, bz, bzp, p, pp, rho)
        real(dp), intent(in) :: r
        real(dp), intent(out) :: bt, btp, bz, bzp, p, pp, rho
        bt = r; btp = 1; bz = ext_bz; bzp = 0; p = 1 - r*r; pp = -2*r; rho = 1
    end subroutine external_field

    function matvec(a, x) result(y)
        type(csc_z_t), intent(in) :: a
        complex(dp), intent(in) :: x(:)
        complex(dp) :: y(size(x))
        integer :: j, k
        y = 0
        do j = 1, a%ncol
            do k = a%col_ptr(j), a%col_ptr(j + 1) - 1
                y(a%row_idx(k)) = y(a%row_idx(k)) + a%val(k)*x(j)
            end do
        end do
    end function matvec

    subroutine spectrum_space(space, stiffness, mass, lam)
        type(cylinder_space_t), intent(in) :: space
        type(csc_z_t), intent(in) :: stiffness, mass
        real(dp), allocatable, intent(out) :: lam(:)
        complex(dp), allocatable :: ka(:, :), ma(:, :), work(:)
        real(dp), allocatable :: rwork(:)
        integer :: j, k, info
        allocate(ka(space%ndof, space%ndof), ma(space%ndof, space%ndof))
        ka = 0; ma = 0
        do j = 1, space%ndof
            do k = stiffness%col_ptr(j), stiffness%col_ptr(j + 1) - 1
                ka(stiffness%row_idx(k), j) = stiffness%val(k)
                ma(mass%row_idx(k), j) = mass%val(k)
            end do
        end do
        allocate(lam(space%ndof), work(4*space%ndof), rwork(3*space%ndof))
        call zhegv(1, 'N', 'U', space%ndof, ka, space%ndof, ma, space%ndof, lam, &
            work, size(work), rwork, info)
        if (info /= 0) error stop 'dense generalized eigensolve'
    end subroutine spectrum_space

    subroutine spectrum(field, a, m, kz, gam, n, deg, lam)
        procedure(cylinder_equilibrium) :: field
        real(dp), intent(in) :: a, kz, gam
        integer, intent(in) :: m, n, deg
        real(dp), allocatable, intent(out) :: lam(:)
        type(cylinder_space_t) :: space
        type(csc_z_t) :: stiffness, mass
        type(fortsparse_status_t) :: status
        complex(dp), allocatable :: ka(:, :), ma(:, :), work(:)
        real(dp), allocatable :: rwork(:)
        integer :: ok, j, k, info
        call cylinder_space_init(space, edges(n, a), deg, m, kz, ok)
        if (ok /= 0) error stop 'space'
        call cylinder_assemble(space, field, gam, 1.0_dp, stiffness, mass, status)
        if (status%code /= FORTSPARSE_OK) error stop 'assembly'
        allocate(ka(space%ndof, space%ndof), ma(space%ndof, space%ndof))
        ka = 0; ma = 0
        do j = 1, space%ndof
            do k = stiffness%col_ptr(j), stiffness%col_ptr(j + 1) - 1
                ka(stiffness%row_idx(k), j) = stiffness%val(k)
                ma(mass%row_idx(k), j) = mass%val(k)
            end do
        end do
        allocate(lam(space%ndof), work(4*space%ndof), rwork(3*space%ndof))
        call zhegv(1, 'N', 'U', space%ndof, ka, space%ndof, ma, space%ndof, lam, &
            work, size(work), rwork, info)
        if (info /= 0) error stop 'dense generalized eigensolve'
    end subroutine spectrum

    function edges(n, a) result(r)
        integer, intent(in) :: n
        real(dp), intent(in) :: a
        real(dp) :: r(0:n)
        integer :: i
        r = [(a*i/n, i = 0, n)]
    end function edges

    pure subroutine uniform_field(r, bt, btp, bz, bzp, p, pp, rho)
        real(dp), intent(in) :: r
        real(dp), intent(out) :: bt, btp, bz, bzp, p, pp, rho
        bt = 0*r; btp = 0; bz = 1; bzp = 0; p = 0.3_dp; pp = 0; rho = 1
    end subroutine uniform_field

    pure subroutine kink_field(r, bt, btp, bz, bzp, p, pp, rho)
        real(dp), intent(in) :: r
        real(dp), intent(out) :: bt, btp, bz, bzp, p, pp, rho
        bt = r; btp = 1; bz = sqrt(0.25_dp + r*r); bzp = r/bz
        p = 1.5_dp*(4 - r*r); pp = -3*r; rho = 1
    end subroutine kink_field

    pure subroutine pinch_field(r, bt, btp, bz, bzp, p, pp, rho)
        real(dp), intent(in) :: r
        real(dp), intent(out) :: bt, btp, bz, bzp, p, pp, rho
        bt = r; btp = 1; bz = 1e-8_dp; bzp = 0
        p = 4 - r*r; pp = -2*r; rho = 1
    end subroutine pinch_field

    subroutine bessel_prime_zeros(m, zeros)
        ! Positive zeros of J_m' by sign scan and bisection.
        integer, intent(in) :: m
        real(dp), intent(out) :: zeros(:)
        real(dp) :: x0, x1, xm, f0, fm
        integer :: found, it
        found = 0; x0 = 0.5_dp; f0 = jp(m, x0)
        do while (found < size(zeros))
            x1 = x0 + 0.01_dp
            if (f0*jp(m, x1) < 0) then
                do it = 1, 200
                    xm = (x0 + x1)/2; fm = jp(m, xm)
                    if (f0*fm <= 0) then
                        x1 = xm
                    else
                        x0 = xm; f0 = fm
                    end if
                end do
                found = found + 1; zeros(found) = (x0 + x1)/2
                x0 = x1 + 1e-6_dp
            else
                x0 = x1
            end if
            f0 = jp(m, x0)
        end do
    end subroutine bessel_prime_zeros

    real(dp) function jp(m, x)
        integer, intent(in) :: m
        real(dp), intent(in) :: x
        jp = (bessel_jn(m - 1, x) - bessel_jn(m + 1, x))/2
    end function jp
end program test_ideal_cylinder
