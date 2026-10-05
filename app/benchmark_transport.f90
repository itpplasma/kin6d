program benchmark_transport
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use kin6d_transport_laboratory, only: transport_state_t, initialize_transport, &
        solve_transport, solve_transport_tangent
    use kin6d_transport_certificate, only: transport_certificate_t, certify_transport, &
        certify_transport_tangent, transport_certificate_accepts, &
        certify_transport_sample, certify_transport_tangent_sample, &
        transport_observable_band
    use kin6d_transport_leaves, only: transport_source_affine, transport_plot_x, &
        transport_cosine_moment_cell_interval, transport_scalar_add_interval, &
        transport_accuracy_interval
    use fortnum_interval, only: interval_t, interval, rad
    implicit none
    integer, parameter :: modes = 2, points = 96, trials = 5
    integer, parameter :: solves = 100, checks = 3, renders = 3
    integer, parameter :: sweep(4) = [32, 64, 256, 4096]
    integer :: cells, parsed, case_index, first, last
    character(len=64) :: argument

    first = 0; last = 0; cells = 256
    if (command_argument_count() > 1) &
        error stop 'usage: benchmark_transport [cells|sweep]'
    if (command_argument_count() == 1) then
        call get_command_argument(1, argument)
        if (trim(argument) == 'sweep') then
            first = 1; last = size(sweep)
        else
            read (argument, *, iostat=parsed) cells
            if (parsed /= 0) error stop 'cells must be an integer or sweep'
        end if
    end if
    if (cells < 4) error stop 'at least four velocity cells required'
    if (iand(cells, cells - 1) /= 0) &
        error stop 'power-of-two cells required for exactly matched benchmark source'
    print '(a)', 'benchmark=general_damped_transport; nu=1; modes=2; E=0'
    print '(a)', 'source=exact_dyadic_Fourier_x_P1_v; certification=continuous_sup'
    print '(a)', 'goal=normalized_cosine_density; tangent=nu_at_fixed_source'
    print '(a)', 'goal_scale=0.25; relative_tolerance=0.01; absolute_floor=1e-8'
    print '(a)', 'postprocessing=96 whole_x_cell state/tangent bands plus goals'
    print '(a)', 'adjoint_goal_gradient=unavailable; external_speed_rank=unestablished'
    print '(a)', 'standalone primal prepare/reprepare are excluded from combined totals'
    do case_index = first, last
        if (case_index > 0) cells = sweep(case_index)
        call measure(cells)
    end do
contains
    subroutine measure(cell_count)
        integer, intent(in) :: cell_count
        type(transport_state_t) :: state
        type(transport_certificate_t) :: certificate
        real(dp), allocatable :: dsource_cos(:, :), dsource_sin(:, :)
        real(dp), allocatable :: dcos(:, :), dsin(:, :)
        real(dp) :: x(points + 1), state_lo(points), state_hi(points)
        real(dp) :: tangent_lo(points), tangent_hi(points)
        real(dp) :: times(7), start, finish, setup_time, combined_time, first_pipeline
        real(dp) :: goal_lo, goal_hi, derivative_lo, derivative_hi
        integer :: j, trial, repeat, evaluator, attributed_bytes
        logical :: ok, strict, accepted, bands_valid, targets_met
        character(len=16) :: evaluator_name

        call cpu_time(start)
        call initialize_transport(state, cell_count, modes, 1.0_dp, ok)
        if (.not. ok) error stop 'benchmark initialization failed'
        allocate (dsource_cos(0:modes, cell_count + 1), &
                  dsource_sin(0:modes, cell_count + 1), &
                  dcos(0:modes, cell_count + 1), dsin(0:modes, cell_count + 1))
        dsource_cos = 0.0_dp; dsource_sin = 0.0_dp
        do j = 1, cell_count + 1
            call transport_source_affine(state%v(j), 1.0_dp, 0.125_dp, &
                state%source_cos(0, j))
            call transport_source_affine(state%v(j), 0.25_dp, -0.25_dp, &
                state%source_cos(1, j))
            call transport_source_affine(state%v(j), 0.0_dp, 0.125_dp, &
                state%source_sin(2, j))
        end do
        do j = 1, points + 1
            call transport_plot_x(real(j - 1, dp), real(points, dp), x(j))
        end do
        call cpu_time(finish)
        setup_time = finish - start
        attributed_bytes = (size(state%v) + size(state%source_cos) + &
            size(state%source_sin) + size(state%cosine) + size(state%sine) + &
            size(dsource_cos) + size(dsource_sin) + size(dcos) + size(dsin) + &
            size(x) + size(state_lo) + size(state_hi) + size(tangent_lo) + &
            size(tangent_hi))*storage_size(start)/8
        print '(a,i0,a,i0,a,es24.16)', 'cells=', cell_count, &
            ' attributed_real_buffer_bytes=', attributed_bytes, &
            ' setup_cpu_s=', setup_time

        do evaluator = 1, 2
            strict = evaluator == 2
            evaluator_name = 'bernstein'
            if (strict) evaluator_name = 'strict_box'
            ! Time the actual first combined pipeline, without a reference solution.
            call cpu_time(start)
            call solve_transport(state, ok)
            call solve_transport_tangent(state, 1.0_dp, dsource_cos, dsource_sin, &
                dcos, dsin, ok)
            call certify_transport_tangent(state, 1.0_dp, dsource_cos, dsource_sin, &
                dcos, dsin, certificate, ok, strict=strict)
            call postprocess(state, certificate, dcos, dsin, x, &
                state_lo, state_hi, tangent_lo, tangent_hi, &
                goal_lo, goal_hi, derivative_lo, derivative_hi, bands_valid)
            targets_met = prepared_targets_met(goal_lo, goal_hi, &
                derivative_lo, derivative_hi)
            call cpu_time(finish)
            first_pipeline = finish - start
            ! Verification precedes repeated trials and stays outside timed phases.
            if (.not. ok) error stop 'warm state/tangent certificate failed'
            if (.not. bands_valid) error stop 'warm postprocessing failed'
            accepted = transport_certificate_accepts(state, certificate%state_radius, &
                strict=strict)
            if (.not. accepted) error stop 'warm certificate rejected'
            if (transport_certificate_accepts(state, 0.0_dp, strict=strict)) &
                error stop 'zero radius unexpectedly accepted for curved state'
            print '(a,a,a,es24.16,a,es24.16)', 'evaluator=', trim(evaluator_name), &
                ' state_radius=', certificate%state_radius, &
                ' nu_tangent_radius=', certificate%tangent_radius
            print '(a,es24.16,a,es24.16,a,l1)', &
                'first_combined_pipeline_cpu_s=', first_pipeline, &
                ' setup_plus_first_pipeline_cpu_s=', setup_time + first_pipeline, &
                ' requested_targets_met=', targets_met
            print '(a,2es24.16)', 'Q_interval=', goal_lo, goal_hi
            print '(a,2es24.16)', 'dQ_dnu_interval=', derivative_lo, derivative_hi
            call report_accuracy(goal_lo, goal_hi, 'Q')
            call report_accuracy(derivative_lo, derivative_hi, 'dQ_dnu')
            print '(a)', 'trial,cells,evaluator,solve_cpu_s,prepare_cpu_s,'// &
                'accept_reprepare_cpu_s,tangent_solve_cpu_s,tangent_prepare_cpu_s,'// &
                'band_postprocess_cpu_s,prepared_target_accept_cpu_s,'// &
                'combined_warm_cpu_s,setup_plus_warm_cpu_s'
            do trial = 1, trials
                call cpu_time(start)
                do repeat = 1, solves
                    call solve_transport(state, ok)
                end do
                call cpu_time(finish)
                times(1) = (finish - start)/solves
                call cpu_time(start)
                do repeat = 1, checks
                    call certify_transport(state, certificate, ok, strict=strict)
                end do
                call cpu_time(finish)
                times(2) = (finish - start)/checks
                call cpu_time(start)
                do repeat = 1, checks
                    accepted = transport_certificate_accepts(state, &
                        certificate%state_radius, strict=strict)
                end do
                call cpu_time(finish)
                times(3) = (finish - start)/checks
                call cpu_time(start)
                do repeat = 1, solves
                    call solve_transport_tangent(state, 1.0_dp, &
                        dsource_cos, dsource_sin, dcos, dsin, ok)
                end do
                call cpu_time(finish)
                times(4) = (finish - start)/solves
                call cpu_time(start)
                do repeat = 1, checks
                    call certify_transport_tangent(state, 1.0_dp, dsource_cos, &
                        dsource_sin, dcos, dsin, certificate, ok, strict=strict)
                end do
                call cpu_time(finish)
                times(5) = (finish - start)/checks
                call cpu_time(start)
                do repeat = 1, renders
                    call postprocess(state, certificate, dcos, dsin, x, &
                        state_lo, state_hi, tangent_lo, tangent_hi, &
                        goal_lo, goal_hi, derivative_lo, derivative_hi, bands_valid)
                end do
                call cpu_time(finish)
                times(6) = (finish - start)/renders
                call cpu_time(start)
                do repeat = 1, solves
                    targets_met = prepared_targets_met(goal_lo, goal_hi, &
                        derivative_lo, derivative_hi)
                end do
                call cpu_time(finish)
                times(7) = (finish - start)/solves
                ! Verification and output are outside the raw timed phases.
                if (.not. ok) error stop 'timed preparation failed'
                if (.not. accepted) error stop 'timed acceptance failed'
                if (.not. bands_valid) error stop 'timed postprocessing failed'
                if (.not. all(ieee_is_finite(times))) error stop 'invalid timer trial'
                combined_time = times(1) + sum(times(4:7))
                print '(i0,",",i0,",",a,9(",",es24.16))', trial, cell_count, &
                    trim(evaluator_name), times, combined_time, &
                    setup_time + combined_time
                if (.not. targets_met) print '(a)', 'requested_goal_target=not_met'
            end do
        end do
    end subroutine
    subroutine postprocess(state, certificate, dcos, dsin, x, &
        state_lo, state_hi, tangent_lo, tangent_hi, &
        goal_lo, goal_hi, derivative_lo, derivative_hi, valid)
        type(transport_state_t), intent(in) :: state
        type(transport_certificate_t), intent(in) :: certificate
        real(dp), intent(in) :: dcos(0:, :), dsin(0:, :), x(:)
        real(dp), intent(out) :: state_lo(:), state_hi(:), tangent_lo(:), tangent_hi(:)
        real(dp), intent(out) :: goal_lo, goal_hi, derivative_lo, derivative_hi
        logical, intent(out) :: valid
        type(interval_t) :: cell
        integer :: j
        logical :: ok

        valid = .true.
        do j = 1, size(x) - 1
            cell = interval(x(j), x(j + 1))
            call certify_transport_sample(state, certificate, cell, &
                interval(0.3125_dp), state_lo(j), state_hi(j), ok)
            valid = valid .and. ok
            call certify_transport_tangent_sample(state, certificate, dcos, dsin, &
                cell, interval(0.3125_dp), tangent_lo(j), tangent_hi(j), ok)
            valid = valid .and. ok
        end do
        call goal_band(state, state%cosine, certificate%state_radius, &
            goal_lo, goal_hi, ok)
        valid = valid .and. ok
        call goal_band(state, dcos, certificate%tangent_radius, &
            derivative_lo, derivative_hi, ok)
        valid = valid .and. ok
    end subroutine

    subroutine goal_band(state, coefficients, radius, lower, upper, ok)
        type(transport_state_t), intent(in) :: state
        real(dp), intent(in) :: coefficients(0:, :), radius
        real(dp), intent(out) :: lower, upper
        logical, intent(out) :: ok
        type(interval_t) :: total, local, updated
        integer :: j

        total = interval(0.0_dp)
        do j = 1, size(state%v) - 1
            call transport_cosine_moment_cell_interval(interval(state%v(j)), &
                interval(state%v(j + 1)), interval(coefficients(1, j)), &
                interval(coefficients(1, j + 1)), local)
            call transport_scalar_add_interval(total, local, updated)
            total = updated
        end do
        call transport_observable_band(total, interval(2.0_dp), &
            radius, lower, upper, ok)
    end subroutine

    logical function prepared_targets_met(lower, upper, dlower, dupper) result(met)
        real(dp), intent(in) :: lower, upper, dlower, dupper
        type(interval_t) :: radius, relative, target, derivative_radius

        radius = interval(rad(interval(lower, upper)))
        derivative_radius = interval(rad(interval(dlower, dupper)))
        call transport_accuracy_interval(radius, interval(0.25_dp), interval(0.01_dp), &
            interval(1.0e-8_dp), relative, target)
        met = radius%hi <= target%lo
        met = met .and. derivative_radius%hi <= target%lo
    end function

    subroutine report_accuracy(lower, upper, label)
        real(dp), intent(in) :: lower, upper
        character(len=*), intent(in) :: label
        type(interval_t) :: radius, relative, target

        radius = interval(rad(interval(lower, upper)))
        call transport_accuracy_interval(radius, interval(0.25_dp), interval(0.01_dp), &
            interval(1.0e-8_dp), relative, target)
        print '(a,a,es24.16,a,es24.16,a,l1)', trim(label), &
            ' normalized_radius_upper=', relative%hi, ' target_lower=', target%lo, &
            ' requested_radius_met=', radius%hi <= target%lo
    end subroutine
end program benchmark_transport
