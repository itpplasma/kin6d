program test_gs_cuda_picard
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use kin6d_grad_shafranov, only: gs_state_t, initialize_gs, solve_gs_profiles, &
        gs_profile_load, gs_profile_value, gs_profile_integral, apply_gs, &
        gs_requested_device, gs_factor_ready
    use kin6d_grad_shafranov_leaves, only: gs_source
    use kin6d_gs_direct, only: gs_direct_free
    use kin6d_gs_options, only: gs_set_options
    use kin6d_gs_timers, only: gs_wall_s, gs_timer_reset, gs_timer_total, &
        gs_timer_calls, GS_T_FACTOR, GS_T_DIRECT
    implicit none
    ! The existing affine P2 oracle in test_gs_p2: psi=R^2 cancels its
    ! profile-dependent source. This target PDE is linear in psi, with no axis.
    integer(int64), parameter :: budget = 1024_int64**3
    real(dp), parameter :: tolerance = 1.0e-11_dp, load_tolerance = 1.0e-13_dp
    real(dp), parameter :: flux(2) = [0.0_dp, 50.0_dp]
    real(dp), parameter :: pp(2) = [-0.01_dp, -0.01_dp]
    real(dp), parameter :: ff(2) = [0.0_dp, 0.5_dp], slope = 0.01_dp
    ! Integral_rectangle |grad(R^2)|^2/R dR dZ, not total magnetic energy.
    real(dp), parameter :: gradient_norm_squared = 15.872_dp
    type :: metrics_t
        real(dp) :: init_s, solve_s, validation_s, factor_s, direct_s
        real(dp) :: element_residual, free_residual, boundary_residual
        real(dp) :: flux_error, load_estimate, source_defect
        real(dp) :: profile_errors(3)
        integer(int64) :: factors, backsolves
        integer :: requested_device, active_device
    end type
    type(gs_state_t) :: cpu, gpu
    type(metrics_t) :: cm, gm
    real(dp), allocatable :: cpu_initial(:), gpu_initial(:), delta(:), action(:)
    real(dp) :: clock, pair_validation_s, cleanup_s, flux_parity, gradient_parity
    integer(int64) :: device_payload
    real(dp) :: profile_parity(3), change, scale, validated_total
    real(dp) :: quadratic, quadratic_roundoff
    integer :: j

    call run_case(cpu, -1, cm, cpu_initial)
    call run_case(gpu, 0, gm, gpu_initial)
    clock = gs_wall_s()
    call require(size(cpu%psi) == size(gpu%psi), 'matched dimension')
    call require(all(shape(cpu%points) == shape(gpu%points)), 'matched point count')
    call require(all(cpu%points == gpu%points), 'matched physical mesh')
    call require(all(shape(cpu%triangles) == shape(gpu%triangles)), &
        'matched topology size')
    call require(all(cpu%triangles == gpu%triangles), 'matched topology')
    call require(all(cpu%boundary .eqv. gpu%boundary), 'matched boundary mask')
    call require(all(cpu_initial == gpu_initial), 'matched initial flux')
    call require(.not. associated(cpu%direct, gpu%direct), 'independent caches')
    delta = gpu%psi - cpu%psi
    call require(all(pack(delta, gpu%boundary) == 0), 'exact boundary parity')
    allocate (action(size(delta)))
    call gs_set_options(.false., .true.)
    call apply_gs(cpu, delta, action)
    quadratic = dot_product(delta, action)
    quadratic_roundoff = 64*epsilon(1.0_dp)*norm2(delta)*norm2(action)
    call require(ieee_is_finite(quadratic), 'finite gradient quadratic')
    call require(ieee_is_finite(quadratic_roundoff), 'finite quadratic roundoff')
    call require(quadratic >= -quadratic_roundoff, 'nonnegative gradient quadratic')
    gradient_parity = sqrt(max(quadratic, 0.0_dp)/gradient_norm_squared)
    flux_parity = maxval(abs(delta))/maxval(abs(cpu%psi))
    call require(ieee_is_finite(gradient_parity), 'finite gradient parity')
    call require(ieee_is_finite(flux_parity), 'finite flux parity')
    call require(max(gradient_parity, flux_parity) <= 1.0e-9_dp, &
        'CPU/GPU flux and weighted-gradient parity')
    profile_parity = 0
    do j = 1, size(delta)
        profile_parity = max(profile_parity, &
            abs(profiles(gpu%psi(j)) - profiles(cpu%psi(j))))
    end do
    change = maxval(abs(delta))
    scale = max(maxval(abs(cpu%psi)), maxval(abs(gpu%psi)))
    call profile_close(profile_parity(1), 0.0_dp, slope*change, 'FFprime parity')
    call profile_close(profile_parity(2), 0.0_dp, slope*change, &
        'pressure primitive parity')
    call profile_close(profile_parity(3), 0.0_dp, slope*scale*change, &
        'FFprime primitive parity')
    call gs_set_options(.true., .true.)
    pair_validation_s = gs_wall_s() - clock
    device_payload = gpu%direct%stats%device_peak_bytes
    clock = gs_wall_s()
    call gs_direct_free(cpu%direct)
    call gs_direct_free(gpu%direct)
    cleanup_s = gs_wall_s() - clock
    validated_total = cm%init_s + cm%solve_s + cm%validation_s &
        + gm%init_s + gm%solve_s + gm%validation_s + pair_validation_s + cleanup_s
    print '(a)', 'observation=single_CPU_first_matched_pair; no speedup claim'
    print '(a)', 'case=affine_profile_P2_Picard; target_PDE=linear_in_psi'
    print '(a)', 'scope=fixed_boundary; no_closed_surfaces; certification=unavailable'
    print '(a)', 'solve_cost=whole_outer_loop; two_adaptive_profile_loads_per_outer'
    print '(a,i0)', 'dofs=', size(delta)
    call report('CPU', cpu, cm)
    call report('GPU', gpu, gm)
    print '(a,2es17.9)', 'CPU_GPU_relative_flux_weighted_gradient_difference=', &
        flux_parity, gradient_parity
    print '(a,2es17.9)', 'difference_gradient_quadratic_roundoff_allowance=', &
        quadratic, quadratic_roundoff
    print '(a,es17.9)', 'analytic_weighted_gradient_norm_squared=', &
        gradient_norm_squared
    print '(a,3es17.9)', 'CPU_GPU_FFprime_pressure_primitive_FFprime_primitive=', &
        profile_parity
    print '(a,3es17.9)', 'pair_validation_cleanup_aggregate_validated_total_s=', &
        pair_validation_s, cleanup_s, validated_total
    print '(a,2(i0,1x))', 'GPU_device_payload_cap_bytes=', device_payload, budget
    print '(a)', 'PASS: complete affine-profile public CPU/GPU Picard fixed point'
contains
    subroutine run_case(state, device, metric, initial)
        type(gs_state_t), intent(out) :: state
        integer, intent(in) :: device
        type(metrics_t), intent(out) :: metric
        real(dp), allocatable, intent(out) :: initial(:)
        real(dp), parameter :: contour(2, 4) = reshape( &
            [5.8_dp, -0.4_dp, 6.6_dp, -0.4_dp, &
             6.6_dp, 0.4_dp, 5.8_dp, 0.4_dp], [2, 4])
        real(dp), allocatable :: rhs(:), applied(:)
        real(dp) :: started, exact, actual_profiles(3), exact_profiles(3), source
        real(dp) :: expected_source, primitive_budget
        integer :: status, i
        call gs_set_options(.true., .true.)
        started = gs_wall_s()
        if (device < 0) then
            call initialize_gs(state, contour, 0.012_dp, status, 2)
        else
            call initialize_gs(state, contour, 0.012_dp, status, 2, device=device, &
                device_memory_limit_bytes=budget)
        end if
        call require(status == 0, 'affine P2 initialization')
        call require(.not. state%mapped, 'affine profile geometry')
        state%psi = 0
        where (state%boundary) state%psi = state%points(1, :)**2
        initial = state%psi
        metric%init_s = gs_wall_s() - started
        call gs_timer_reset()
        started = gs_wall_s()
        call solve_gs_profiles(state, 1.0_dp, flux, pp, ff, status, tolerance)
        metric%solve_s = gs_wall_s() - started
        started = gs_wall_s()
        metric%factors = gs_timer_calls(GS_T_FACTOR)
        metric%backsolves = gs_timer_calls(GS_T_DIRECT)
        metric%factor_s = gs_timer_total(GS_T_FACTOR)
        metric%direct_s = gs_timer_total(GS_T_DIRECT)
        call require(status == 0, 'complete public Picard solve')
        call require(all(ieee_is_finite(state%psi)), 'finite Picard solution')
        call require(ieee_is_finite(state%residual), 'finite stored residual')
        call require(state%residual <= tolerance, 'original outer residual')
        call require(state%iterations > 1, 'nontrivial outer iteration')
        call require(state%iterations <= 1000, 'original outer cap')
        call require(metric%factors == 1_int64, 'one factor across outer solves')
        call require(metric%backsolves >= int(state%iterations, int64), &
            'observed cached backsolves across outer iterations')
        call require(associated(state%direct), 'owned direct cache')
        call require(state%direct%factored, 'factor persists after whole solve')
        call require(gs_factor_ready(state), 'public factor cache ready')
        call require(gs_timer_calls(GS_T_FACTOR) == metric%factors, &
            'ready query reuses factor')
        metric%requested_device = gs_requested_device(state)
        metric%active_device = state%direct%stats%active_device
        call require(metric%requested_device == device, 'state device persists')
        call require(metric%active_device == device, &
            'actual native factor device; no fallback')
        if (device >= 0) then
            call require(state%direct%stats%device_peak_bytes > 0, &
                'positive GPU payload')
            call require(state%direct%stats%device_peak_bytes <= budget, &
                'positive GPU cap preserved')
        end if
        call require(all(pack(state%psi, state%boundary) &
            == pack(initial, state%boundary)), 'fixed boundary preserved')
        metric%flux_error = maxval(abs(state%psi - state%points(1, :)**2))
        call require(metric%flux_error <= 2.0e-7_dp, 'original exact R^2 oracle')
        allocate (rhs(size(state%psi)), applied(size(state%psi)))
        call gs_profile_load(state, 1.0_dp, flux, pp, ff, rhs, status, &
            relative_tolerance=load_tolerance, error_estimate=metric%load_estimate)
        call require(status == 0, 'independent profile load and estimator budget')
        call require(ieee_is_finite(metric%load_estimate), 'finite load estimate')
        call gs_set_options(.false., .true.)
        call apply_gs(state, state%psi, applied)
        metric%element_residual = norm2(applied - rhs) &
            /max(norm2(rhs), tiny(1.0_dp))
        metric%free_residual = norm2(pack(applied - rhs, .not. state%boundary)) &
            /max(norm2(pack(rhs, .not. state%boundary)), tiny(1.0_dp))
        metric%boundary_residual = maxval(abs(pack(applied - rhs, state%boundary)))
        call require(ieee_is_finite(metric%element_residual), &
            'finite independent element residual')
        call require(metric%element_residual <= tolerance, &
            'independent full Picard equation residual')
        metric%profile_errors = 0
        metric%source_defect = 0
        do i = 1, size(state%psi)
            exact = state%points(1, i)**2
            actual_profiles = profiles(state%psi(i))
            exact_profiles = [slope*exact, slope*(50 - exact), &
                0.5_dp*slope*(exact**2 - 50.0_dp**2)]
            metric%profile_errors = max(metric%profile_errors, &
                abs(actual_profiles - exact_profiles))
            call profile_close(gs_profile_value(flux, pp, state%psi(i)), &
                -slope, 0.0_dp, 'constant pressure derivative')
            call profile_close(actual_profiles(1), exact_profiles(1), &
                slope*metric%flux_error, 'analytic FFprime profile')
            call profile_close(actual_profiles(2), exact_profiles(2), &
                slope*metric%flux_error, 'analytic pressure primitive')
            primitive_budget = slope*(abs(exact) + 0.5_dp*metric%flux_error) &
                *metric%flux_error
            call profile_close(actual_profiles(3), exact_profiles(3), &
                primitive_budget, 'analytic FFprime primitive')
            call gs_source(state%points(1, i), 1.0_dp, &
                gs_profile_value(flux, pp, state%psi(i)), actual_profiles(1), source)
            expected_source = slope*(state%psi(i) - exact)/state%points(1, i)
            call profile_close(source, expected_source, 0.0_dp, &
                'independent R^2 source cancellation')
            metric%source_defect = max(metric%source_defect, abs(source))
        end do
        call gs_set_options(.true., .true.)
        metric%validation_s = gs_wall_s() - started
    end subroutine run_case

    function profiles(psi) result(values)
        real(dp), intent(in) :: psi
        real(dp) :: values(3)
        values = [gs_profile_value(flux, ff, psi), &
            gs_profile_integral(flux, pp, psi), gs_profile_integral(flux, ff, psi)]
    end function profiles

    subroutine report(label, state, metric)
        character(*), intent(in) :: label
        type(gs_state_t), intent(in) :: state
        type(metrics_t), intent(in) :: metric
        print '(a,a,i0)', label, '_outer_iterations=', state%iterations
        print '(a,a,3es17.9)', label, '_stored_element_residual_exact_flux_error=', &
            state%residual, metric%element_residual, metric%flux_error
        print '(a,a,2es17.9)', label, &
            '_lifted_free_relative_boundary_absolute_residual=', &
            metric%free_residual, metric%boundary_residual
        print '(a,a,2es17.9)', label, '_load_estimate_nodal_source_defect=', &
            metric%load_estimate, metric%source_defect
        print '(a,a,3es17.9)', label, '_FFprime_pressure_FFprime_primitive_errors=', &
            metric%profile_errors
        print '(a,a,3es17.9)', label, '_init_whole_Picard_init_plus_solve_s=', &
            metric%init_s, metric%solve_s, metric%init_s + metric%solve_s
        print '(a,a,2es17.9)', label, '_validation_validated_route_total_s=', &
            metric%validation_s, metric%init_s + metric%solve_s + metric%validation_s
        print '(a,a,2(i0,1x))', label, '_requested_active_device=', &
            metric%requested_device, metric%active_device
        print '(a,a,2(i0,1x))', label, '_factor_backsolve_calls=', &
            metric%factors, metric%backsolves
        print '(a,a,2es17.9)', label, '_factor_backsolve_s=', &
            metric%factor_s, metric%direct_s
    end subroutine report

    subroutine profile_close(actual, expected, propagated, message)
        real(dp), intent(in) :: actual, expected, propagated
        character(*), intent(in) :: message
        real(dp) :: roundoff
        roundoff = 64*epsilon(1.0_dp)*max(1.0_dp, abs(expected))
        call require(ieee_is_finite(actual), message)
        call require(abs(actual - expected) <= propagated + roundoff, message)
    end subroutine profile_close

    subroutine require(ok, message)
        logical, intent(in) :: ok
        character(*), intent(in) :: message
        if (ok) return
        print '(a)', message
        error stop 'public affine-profile GS CUDA Picard'
    end subroutine require
end program test_gs_cuda_picard
