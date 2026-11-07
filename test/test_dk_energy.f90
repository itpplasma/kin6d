!> Axisymmetric drift-kinetic solve on the energy grid (#216, #184).
!!   rigid_shift   with the momentum-conserving operator the D31 drive has
!!                 the exact solution h = x xi/b, so D31 = sum w x^2/3; the
!!                 test-particle operator misses it;
program test_dk_energy
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use kin6d_dk_energy
    use kin6d_dk_collop, only: collop_species_t, collop_grid_t, collop_make_grid, &
        collop_t, collop_build, COLLOP_FULL, COLLOP_LORENTZ, collop_qe, COLLOP_OK
    implicit none
    real(dp), parameter :: pi = 4*atan(1.0_dp)
    character(64) :: check
    logical :: ok
    call get_command_argument(1, check)
    select case (trim(check))
    case ('rigid_shift'); ok = rigid_shift()
    case ('jvp'); ok = jvp()
    case ('operator_cache'); ok = operator_cache()
    case ('batch'); ok = batch()
    case ('pfirsch_schluter'); ok = pfirsch_schluter()
    case ('profile'); ok = profile()
    case default
        print '(a)', 'unknown check '//trim(check)
        ok = .false.
    end select
    if (.not. ok) error stop 1
contains
    logical function profile() result(ok)
        type(collop_grid_t) :: grid
        type(energy_axi_t) :: value
        integer :: status
        integer(int64) :: start, finish, rate
        call collop_make_grid(6.0_dp, 2, 5, grid, status)
        call system_clock(start, rate)
        call energy_axi_coefficients(field(17, 0.1_dp), 1.0_dp, deuterons(1e20_dp), &
            17.0_dp, grid, 16, COLLOP_FULL, value, status)
        call system_clock(finish)
        if (status /= ENERGY_OK) error stop 'energy profile failed'
        print '(a)', 'synthetic: theta=17, speed grid=(xmax=6,elements=2,order=5), Legendre=16'
        print '(a,es18.10)', 'elapsed [s]: ', real(finish - start, dp)/real(rate, dp)
        print '(a,3es18.10)', 'D31, D32, k: ', value%d31, value%d32, value%k
        print '(a,3es18.10)', 's tangents: ', value%dd31_ds, value%dd32_ds, value%dk_ds
        print '(a,2es18.10)', 'backward errors: ', value%residual
        ok = .true.
    end function profile

    !> Landreman (2013), equations (2) and (5), for full ion self-collisions:
    !> https://terpconnect.umd.edu/~mattland/assets/notes/20131019-03-PfirschSchluterRegimeIonNeoclassicalFlow.pdf
    logical function pfirsch_schluter() result(ok)
        type(collop_grid_t) :: grid
        type(energy_axi_t) :: low, high
        real(dp), parameter :: eps = 0.01_dp, k0 = -1.97213_dp, k1 = -0.09361_dp
        real(dp) :: ratio, reference, limit
        integer :: status
        call collop_make_grid(6.0_dp, 4, 7, grid, status)
        if (status /= COLLOP_OK) error stop 'PS grid failed'
        call energy_axi_coefficients(field(9, eps), 1.0_dp, deuterons(1e24_dp), &
            17.0_dp, grid, 8, COLLOP_FULL, low, status)
        if (status /= ENERGY_OK) error stop 'PS lower-collisionality solve failed'
        call energy_axi_coefficients(field(9, eps), 1.0_dp, deuterons(1e25_dp), &
            17.0_dp, grid, 8, COLLOP_FULL, high, status)
        if (status /= ENERGY_OK) error stop 'PS higher-collisionality solve failed'
        ! Analytic surface integrals for b=1/(1+eps*cos(theta)), parallel
        ! derivative kappa*b*d_theta and surface measure 1/b**2.
        ratio = 2*(1 - eps**2)**2/((1 + eps**2/2)*(1 + sqrt(1 - eps**2)))
        reference = k0 + k1*ratio
        ! The leading finite-collisionality correction is proportional to
        ! 1/nu**2; extrapolate between rates differing by a factor of ten.
        limit = (100*high%k - low%k)/99
        print '(a,3es18.10)', 'PS low, high, extrapolated k: ', low%k, high%k, limit
        print '(a,es18.10)', 'PS geometry-dependent reference: ', reference
        ok = abs(low%k - reference) < 6e-3_dp &
            .and. abs(high%k - reference) < 1e-4_dp &
            .and. abs(limit - reference) < 2e-5_dp &
            .and. max(abs(low%d31 - 0.5_dp), abs(high%d31 - 0.5_dp)) < 1e-4_dp
    end function pfirsch_schluter

    logical function batch() result(ok)
        type(collop_grid_t) :: grid
        type(collop_species_t) :: sp(3)
        type(energy_axi_t) :: value(3), fresh
        real(dp) :: b(9, 3), lnl(3)
        integer :: status(3), st, j
        call collop_make_grid(6.0_dp, 1, 3, grid, st)
        sp = [deuterons(-1.0_dp), deuterons(1e20_dp), deuterons(2e20_dp)]
        do j = 1, 3
            b(:, j) = field(9, 0.1_dp)
        end do
        lnl = 17
        call energy_axi_batch(b, [1.0_dp, 1.0_dp, 1.0_dp], sp, lnl, grid, 8, &
            COLLOP_FULL, value, status)
        if (status(1) /= ENERGY_BAD_INPUT .or. any(status(2:) /= ENERGY_OK)) &
            error stop 'bad surface discarded valid batch members'
        ok = .true.
        do j = 2, 3
            call energy_axi_coefficients(b(:, j), 1.0_dp, sp(j), lnl(j), grid, 8, &
                COLLOP_FULL, fresh, st)
            if (st /= ENERGY_OK) error stop 'fresh batch oracle failed'
            ok = ok .and. maxval(abs([value(j)%d31 - fresh%d31, value(j)%d32 - fresh%d32, &
                value(j)%k - fresh%k, value(j)%dk_ds - fresh%dk_ds])) < 1e-10_dp
        end do
        b(:, 2) = 0
        call energy_axi_batch(b, [1.0_dp, 1.0_dp, 1.0_dp], sp, lnl, grid, 8, &
            COLLOP_FULL, value, status)
        ok = ok .and. all(status == [ENERGY_BAD_INPUT, ENERGY_BAD_INPUT, ENERGY_OK])
    end function batch

    pure function field(n, eps) result(b)
        integer, intent(in) :: n
        real(dp), intent(in) :: eps
        real(dp) :: b(n)
        integer :: j
        b = [(1/(1 + eps*cos(2*pi*(j - 1)/n)), j=1, n)]
    end function field

    !> Deuterons at 1 keV with density n [m^-3].
    function deuterons(n) result(s)
        real(dp), intent(in) :: n
        type(collop_species_t) :: s
        s = collop_species_t(3.3436e-27_dp, 1.0_dp, n, 1e3_dp*collop_qe)
    end function deuterons

    logical function rigid_shift() result(ok)
        type(collop_grid_t) :: grid
        type(energy_axi_t) :: full, tp
        integer :: status
        call collop_make_grid(6.0_dp, 4, 5, grid, status)
        call energy_axi_coefficients(field(33, 0.1_dp), 1.0_dp, deuterons(1e20_dp), 17.0_dp, &
            grid, 24, COLLOP_FULL, full, status)
        if (status /= ENERGY_OK) error stop 'full solve failed'
        call energy_axi_coefficients(field(33, 0.1_dp), 1.0_dp, deuterons(1e20_dp), 17.0_dp, &
            grid, 24, COLLOP_LORENTZ, tp, status)
        if (status /= ENERGY_OK) error stop 'test-particle solve failed'
        print '(a,3es16.8)', 'D31 full, test particle, rigid shift:', full%d31, tp%d31, &
            full%d31_rigid
        print '(a,2f12.6,a,2es9.1)', 'k full, test particle:', full%k, tp%k, '  residuals', &
            full%residual
        ok = abs(full%d31/full%d31_rigid - 1) < 1e-8_dp .and. all(full%residual < 1e-8_dp)
    end function rigid_shift

    logical function jvp() result(ok)
        type(collop_grid_t) :: grid
        type(collop_species_t) :: sp, up, down
        type(energy_axi_t) :: value, plus, minus
        real(dp), parameter :: h = 1e-4_dp
        real(dp) :: exact(3), fd(3)
        integer :: status
        call collop_make_grid(6.0_dp, 1, 3, grid, status)
        sp = deuterons(1e20_dp)
        call energy_axi_coefficients(field(9, 0.1_dp), 1.0_dp, sp, 17.0_dp, &
            grid, 8, COLLOP_FULL, value, status)
        if (status /= ENERGY_OK) error stop 'tangent solve failed'
        up = sp; down = sp
        up%n = sp%n*exp(h); down%n = sp%n*exp(-h)
        call energy_axi_coefficients(field(9, 0.1_dp), 1.0_dp, up, 17.0_dp, &
            grid, 8, COLLOP_FULL, plus, status)
        if (status /= ENERGY_OK) error stop 'positive density solve failed'
        call energy_axi_coefficients(field(9, 0.1_dp), 1.0_dp, down, 17.0_dp, &
            grid, 8, COLLOP_FULL, minus, status)
        if (status /= ENERGY_OK) error stop 'negative density solve failed'
        fd = [plus%d31 - minus%d31, plus%d32 - minus%d32, plus%k - minus%k]/(2*h)
        exact = [value%dd31_ds, value%dd32_ds, value%dk_ds]
        print '(a,3es14.5)', 'log-density tangent errors:', abs(exact - fd)
        ok = all(abs(exact - fd) < 1e-7_dp*(1 + abs(fd)))
        up = sp; down = sp
        up%t = sp%t*exp(h); down%t = sp%t*exp(-h)
        call energy_axi_coefficients(field(9, 0.1_dp), 1.0_dp, up, 17.0_dp, &
            grid, 8, COLLOP_FULL, plus, status)
        if (status /= ENERGY_OK) error stop 'positive temperature solve failed'
        call energy_axi_coefficients(field(9, 0.1_dp), 1.0_dp, down, 17.0_dp, &
            grid, 8, COLLOP_FULL, minus, status)
        if (status /= ENERGY_OK) error stop 'negative temperature solve failed'
        fd = [plus%d31 - minus%d31, plus%d32 - minus%d32, plus%k - minus%k]/(2*h)
        print '(a,3es14.5)', 'log-temperature tangent errors:', abs(-2*exact - fd)
        ok = ok .and. all(abs(-2*exact - fd) < 1e-7_dp*(1 + abs(fd)))
    end function jvp

    logical function operator_cache() result(ok)
        type(collop_grid_t) :: grid, other
        type(collop_species_t) :: sp(1)
        type(collop_t) :: op
        type(energy_axi_t) :: cached, fresh
        real(dp) :: lnl(1, 1)
        integer :: status
        call collop_make_grid(6.0_dp, 1, 3, grid, status)
        sp = [deuterons(2e20_dp)]
        sp%t = 2*sp%t
        lnl = 17
        call collop_build(sp, lnl, grid, COLLOP_FULL, 7, op, status)
        if (status /= COLLOP_OK) error stop 'cache build failed'
        call energy_axi_coefficients(field(9, 0.1_dp), 1.0_dp, deuterons(1e20_dp), &
            17.0_dp, grid, 8, COLLOP_FULL, cached, status, operator=op)
        if (status /= ENERGY_OK) error stop 'cached solve failed'
        call energy_axi_coefficients(field(9, 0.1_dp), 1.0_dp, deuterons(1e20_dp), &
            17.0_dp, grid, 8, COLLOP_FULL, fresh, status)
        if (status /= ENERGY_OK) error stop 'fresh solve failed'
        print '(a,3es14.5)', 'cached operator differences:', cached%d31 - fresh%d31, &
            cached%d32 - fresh%d32, cached%k - fresh%k
        ok = maxval(abs([cached%d31 - fresh%d31, cached%d32 - fresh%d32, &
            cached%k - fresh%k])) < 1e-10_dp
        call collop_make_grid(7.0_dp, 1, 3, other, status)
        call energy_axi_coefficients(field(9, 0.1_dp), 1.0_dp, deuterons(1e20_dp), &
            17.0_dp, other, 8, COLLOP_FULL, cached, status, operator=op)
        print '(a,i3)', 'incompatible grid status:', status
        ok = ok .and. status == ENERGY_BAD_INPUT
    end function operator_cache
end program test_dk_energy
