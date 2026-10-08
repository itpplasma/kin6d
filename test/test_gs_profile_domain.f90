program test_gs_profile_domain
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
    use kin6d_grad_shafranov, only: gs_state_t, initialize_gs, solve_gs_profiles
    implicit none
    type(gs_state_t) :: state
    real(dp) :: contour(2, 4), mu0
    integer :: status, order
    character(32) :: mode, order_arg
    call get_command_argument(1, mode)
    call get_command_argument(2, order_arg)
    order = 1
    if (len_trim(order_arg) > 0) read (order_arg, *) order
    contour = reshape([5.8_dp, -0.4_dp, 6.6_dp, -0.4_dp, &
                         6.6_dp, 0.4_dp, 5.8_dp, 0.4_dp], [2, 4])
    call initialize_gs(state, contour, 0.02_dp, status, order)
    if (status /= 0) error stop 'domain oracle mesh failed'
    select case (trim(mode))
    case ('valid')
        mu0 = 1.0_dp
    case ('zero')
        mu0 = 0.0_dp
    case ('negative')
        mu0 = -1.0_dp
    case ('nan')
        mu0 = ieee_value(0.0_dp, ieee_quiet_nan)
    case default
        error stop 'unknown domain oracle mode'
    end select
    call solve_gs_profiles(state, mu0, [-1.0_dp, 0.0_dp], &
        [1.0_dp, 1.0_dp], [0.0_dp, 0.0_dp], status)
    print *, 'mu0 domain mode/status:', trim(mode), status
    if (trim(mode) == 'valid') then
        if (status /= 0) error stop 'positive finite mu0 rejected'
        if (.not. all(ieee_is_finite(state%psi))) error stop 'nonfinite valid solution'
        if (maxval(state%psi) <= 0) error stop 'nontrivial valid source not solved'
        if (state%residual > 1.0e-8_dp) error stop 'valid weak residual failed'
    else
        if (status /= 1) error stop 'invalid mu0 was not rejected at entry'
        if (maxval(abs(state%psi)) > 0) error stop 'invalid input changed the state'
    end if
    print *, 'PASS: independent finite-positive GS profile domain', trim(mode)
end program test_gs_profile_domain
