program test_gs_device_selection
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use kin6d_grad_shafranov, only: gs_state_t, initialize_gs, solve_gs, &
        gs_requested_device, GS_DEVICE_FAILURE, GS_UNSUPPORTED_DEVICE
    use kin6d_gs_options, only: gs_set_options
    use kin6d_gs_profile_law, only: gs_profile_law_t
    use kin6d_gs_nonlinear, only: solve_gs_cubic_profiles, &
        gs_profile_jacobian_t, cubic_response => gs_profile_response
    use kin6d_gs_free_boundary, only: gs_fb_system_t, gs_fb_config_t, &
        gs_fb_eval_t, gs_fb_factor_t, gs_fb_setup, gs_fb_newton, gs_fb_factor, &
        GS_FB_UNSUPPORTED_DEVICE
    implicit none
    type(gs_state_t) :: cpu, requested, invalid
    type(gs_profile_law_t) :: law
    type(gs_profile_jacobian_t) :: jacobian
    type(gs_fb_system_t) :: system
    type(gs_fb_config_t) :: config
    type(gs_fb_eval_t) :: evaluation
    type(gs_fb_factor_t) :: factor
    real(dp) :: contour(2, 4)
    real(dp), allocatable :: rhs(:), response(:)
    integer :: status

    contour = reshape([2.0_dp, -1.0_dp, 4.0_dp, -1.0_dp, &
        4.0_dp, 1.0_dp, 2.0_dp, 1.0_dp], [2, 4])
    call initialize_gs(invalid, contour, 0.5_dp, status, 2, device=-2)
    call require(status == GS_DEVICE_FAILURE, 'invalid device rejected')
    call require(.not. invalid%ready, 'invalid device has no ready state')
    call initialize_gs(invalid, contour, 0.5_dp, status, 2, &
        device=0, device_memory_limit_bytes=0_int64)
    call require(status == GS_DEVICE_FAILURE, 'nonpositive device cap rejected')

    call initialize_gs(requested, contour, 0.5_dp, status, 2, &
        device=0, device_memory_limit_bytes=1024_int64**3)
    call require(status == 0, 'explicit device state initialized')
    call require(gs_requested_device(requested) == 0, 'device policy retained')
    ! Device availability is intentionally not needed for this control. A
    ! conflicting iterative-only selection must not silently run CPU CG.
    call gs_set_options(.true., .false.)
    call solve_gs(requested, 1.0_dp, 0.1_dp, 0.2_dp, status)
    call require(status == GS_DEVICE_FAILURE, 'explicit request has no fallback')
    call require(maxval(abs(requested%psi)) == 0, 'failed request keeps flux')
    call require(.not. requested%direct%factored, 'no hidden CPU factor')

    ! CPU default retains its existing selectable iterative path.
    call initialize_gs(cpu, contour, 0.5_dp, status, 2)
    call require(status == 0, 'default CPU state')
    call require(gs_requested_device(cpu) == -1, 'default remains CPU')
    call solve_gs(cpu, 1.0_dp, 0.1_dp, 0.2_dp, status)
    call require(status == 0, 'default CPU iterative solve')
    call require(cpu%residual <= 1.0e-8_dp, 'unchanged CPU residual gate')
    call require(.not. cpu%direct%factored, 'iterative control used no factor')
    call gs_set_options(.true., .true.)
    call solve_gs(cpu, 1.0_dp, 0.1_dp, 0.2_dp, status)
    call require(status == 0, 'default CPU direct solve')
    call require(cpu%direct%stats%active_device == -1, 'CPU factor stays CPU')

    ! These Newton paths need indefinite/updated operators and are not part
    ! of the supported fixed-source SPD device contract.
    call solve_gs_cubic_profiles(requested, 1.0_dp, law, status)
    call require(status == GS_UNSUPPORTED_DEVICE, 'cubic Newton not relabelled')
    allocate (rhs(size(requested%psi)), response(size(requested%psi)))
    rhs = 0
    call cubic_response(requested, jacobian, rhs, response, status)
    call require(status == GS_UNSUPPORTED_DEVICE, 'cubic response not relabelled')
    call gs_fb_setup(requested, system, law, config, status=status)
    call require(status == GS_FB_UNSUPPORTED_DEVICE, 'free-boundary setup')
    call gs_fb_newton(requested, system, [real(dp) ::], status)
    call require(status == GS_FB_UNSUPPORTED_DEVICE, 'free-boundary Newton')
    call gs_fb_factor(requested, system, evaluation, factor, status)
    call require(status == GS_FB_UNSUPPORTED_DEVICE, 'free-boundary factor')
    print '(a)', 'PASS: GS per-state device contract; no GPU required'
contains
    subroutine require(ok, message)
        logical, intent(in) :: ok
        character(*), intent(in) :: message
        if (ok) return
        print '(a)', message
        error stop 'GS device selection contract'
    end subroutine
end program test_gs_device_selection
