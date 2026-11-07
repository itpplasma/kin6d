program test_dk_force_balance
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use kin6d_dk_force_balance
    use kin6d_dk_energy, only: ENERGY_OK, ENERGY_BAD_INPUT
    use kin6d_dk_energy_leaves, only: dk_energy_force_balance_vjp
    implicit none
    real(dp), parameter :: p(12) = [1.3_dp, 3e13_dp, 3e-9_dp, -4e11_dp, -5e-11_dp, &
        2.1e4_dp, -2e4_dp, 15.0_dp, 2e6_dp, -2.5_dp, 165.0_dp, 4.8031e-10_dp]
    character(64) :: check
    call get_command_argument(1, check)
    select case (trim(check))
    case ('physics'); call physics()
    case ('derivatives'); call derivatives()
    case ('refusals'); call refusals()
    case default; error stop 'unknown force-balance check'
    end select
contains
    function input(x) result(v)
        real(dp), intent(in) :: x(12)
        type(energy_balance_input_t) :: v
        v = energy_balance_input_t(x(1), x(2), x(3), x(4), x(5), x(6), &
            x(7), x(8), x(9), x(10), x(11), x(12))
    end function input

    subroutine physics()
        type(energy_balance_t) :: v
        integer :: status
        call energy_force_balance_cgs(input(p), v, status)
        if (status /= ENERGY_OK) error stop 'valid force balance refused'
        ! KAMEL 266210c5, python/profile_processor/profile_processor.py,
        ! calc_Er_prof evaluated on the synthetic CGS inputs above.
        if (abs(v%er/(-0.001105383351174903_dp) - 1) > 1e-13_dp) &
            error stop 'KAMEL radial electric field differs'
        if (abs(v%vth/183994.07459402047_dp - 1) > 1e-13_dp) &
            error stop 'KAMEL poloidal flow differs'
        if (v%certification /= 'unavailable') error stop 'unjustified certificate'
        print '(a,es18.10,a,es18.10)', 'Er [statV/cm]: ', v%er, ' Vtheta [cm/s]: ', v%vth
    end subroutine physics

    subroutine derivatives()
        type(energy_balance_t) :: v, plus, minus
        real(dp), parameter :: h = 1e-5_dp, weights(2) = [0.7_dp, -0.3_dp]
        real(dp) :: d(12), g(12), fd(2), tangent(2), pair, exact
        integer :: j, status
        do j = 1, 12
            d = 0; d(j) = p(j)
            call energy_force_balance_cgs(input(p), v, status, input(d))
            if (status /= ENERGY_OK) error stop 'force-balance tangent refused'
            call energy_force_balance_cgs(input(p + h*d), plus, status)
            if (status /= ENERGY_OK) error stop 'positive perturbation refused'
            call energy_force_balance_cgs(input(p - h*d), minus, status)
            if (status /= ENERGY_OK) error stop 'negative perturbation refused'
            fd = [plus%er - minus%er, plus%vth - minus%vth]/(2*h)
            tangent = [v%der, v%dvth]
            if (any(abs(fd - tangent) > 2e-8_dp*max(abs(tangent), [1e-3_dp, 1e5_dp]))) &
                error stop 'force-balance finite-difference mismatch'
        end do
        d = p*[(real(j, dp)/12, j=1, 12)]
        call energy_force_balance_cgs(input(p), v, status, input(d))
        if (status /= ENERGY_OK) error stop 'combined tangent refused'
        call dk_energy_force_balance_vjp(p(1), p(2), p(3), p(4), p(5), p(6), &
            p(7), p(8), p(9), p(10), p(11), p(12), weights(1), weights(2), &
            g(1), g(2), g(3), g(4), g(5), g(6), g(7), g(8), g(9), g(10), g(11), g(12))
        pair = dot_product(g, d)
        exact = dot_product(weights, [v%der, v%dvth])
        if (abs(pair - exact) > 1e-13_dp*max(1.0_dp, abs(exact))) &
            error stop 'force-balance adjoint pairing differs'
    end subroutine derivatives

    subroutine refusals()
        type(energy_balance_t) :: v
        real(dp) :: x(12), direction(12)
        integer :: j, status
        do j = 1, 12
            x = p; x(j) = ieee_value(0.0_dp, ieee_quiet_nan)
            call energy_force_balance_cgs(input(x), v, status)
            if (status /= ENERGY_BAD_INPUT) error stop 'nonfinite input accepted'
            direction = 0; direction(j) = x(j)
            call energy_force_balance_cgs(input(p), v, status, input(direction))
            if (status /= ENERGY_BAD_INPUT) error stop 'nonfinite direction accepted'
        end do
        do j = 1, 12
            if (.not. any(j == [2, 3, 6, 10, 11, 12])) cycle
            x = p; x(j) = 0
            call energy_force_balance_cgs(input(x), v, status)
            if (status /= ENERGY_BAD_INPUT) error stop 'singular input accepted'
        end do
        x = p; x(8) = -1
        call energy_force_balance_cgs(input(x), v, status)
        if (status /= ENERGY_BAD_INPUT) error stop 'negative radius accepted'
        call energy_force_balance_cgs(input(p), v, status)
        if (status /= ENERGY_OK) error stop 'refusal poisoned later valid call'
    end subroutine refusals
end program test_dk_force_balance
