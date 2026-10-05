program benchmark_streaming_comparison_main
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_streaming_comparison, only: streaming_state_t, init_streaming, &
        advance_streaming, evaluate_streaming, streaming_radius
    use kin6d_streaming_leaves, only: streaming_velocity_coordinate
    use kin6d_result_io, only: run_t, result_file_t
    implicit none
    integer :: cells = 128, repeats = 50, trials = 5
    type(streaming_state_t) :: state
    type(run_t) :: run
    type(result_file_t) :: result
    real(dp), allocatable :: velocity(:), a0(:, :), b0(:, :), elapsed(:, :), sorted(:, :)
    real(dp) :: start, finish, scratch
    real(dp) :: value, lower, upper, radius
    integer :: j, trial, order, variant, repeat, k, unit, ios, mid
    character(len=256) :: msg
    logical :: ok
    namelist /benchmark_streaming_comparison/ cells, repeats, trials

    call run%start('benchmark_streaming_comparison')
    if (run%argument_count() > 0) &
        error stop 'usage: benchmark_streaming_comparison [--config FILE]'
    ios = 0; msg = ''
    call run%open_config(unit)
    if (run%has_config()) read (unit, nml=benchmark_streaming_comparison, iostat=ios, &
        iomsg=msg)
    call run%close_config(unit, 'benchmark_streaming_comparison', ios, msg)
    if (cells < 1 .or. repeats < 1 .or. trials < 1) error stop 'invalid benchmark sizes'
    call run%open_resolved(unit)
    write (unit, nml=benchmark_streaming_comparison)
    call run%close_resolved(unit)
    allocate (velocity(cells + 1), a0(0:1, cells + 1), b0(0:1, cells + 1), &
        elapsed(trials, 2), sorted(trials, 2))
    do j = 1, cells + 1
        call streaming_velocity_coordinate(real(j - 1, dp), real(cells, dp), &
            velocity(j))
    end do
    a0(0, :) = 1.0_dp
    a0(1, :) = 0.5_dp
    b0 = 0.0_dp
    call init_streaming(state, velocity, a0, b0, ok)
    if (.not. ok) error stop 'benchmark initialization failed'
    call advance_streaming(state, 12.0_dp, ok)
    if (.not. ok) error stop 'benchmark warm certificate failed'

    print '(a)', 'trial,variant,warm_seconds_per_result'
    do trial = 1, trials
        do order = 1, 2
            variant = order
            if (mod(trial, 2) == 0) variant = 3 - order
            call cpu_time(start)
            do repeat = 1, repeats
                call advance_streaming(state, 12.0_dp, ok, certify=(variant == 2))
                if (.not. ok) error stop 'benchmark evolution failed'
            end do
            call cpu_time(finish)
            elapsed(trial, variant) = (finish - start)/real(repeats, dp)
            print '(i0,a,i0,a,es24.16)', trial, ',', variant, ',', &
                elapsed(trial, variant)
        end do
    end do
    sorted = elapsed
    do variant = 1, 2
        do j = 2, trials
            do k = j, 2, -1
                if (sorted(k, variant) >= sorted(k - 1, variant)) exit
                scratch = sorted(k, variant)
                sorted(k, variant) = sorted(k - 1, variant)
                sorted(k - 1, variant) = scratch
            end do
        end do
    end do
    mid = (trials + 1)/2
    print '(a)', 'variant,median_seconds,min_seconds,max_seconds'
    do variant = 1, 2
        print '(i0,3(a,es24.16))', variant, ',', sorted(mid, variant), ',', &
            sorted(1, variant), ',', sorted(trials, variant)
    end do
    print '(a,es24.16)', 'extra_certificate_seconds=', &
        sorted(mid, 2) - sorted(mid, 1)
    print '(a,f12.4)', 'certified_evolution_over_solve_only=', &
        sorted(mid, 2)/sorted(mid, 1)
    call advance_streaming(state, 12.0_dp, ok)
    if (.not. ok) error stop 'benchmark final certificate failed'
    call evaluate_streaming(state, 0.0_dp, 0.3125_dp, value, lower, upper, ok)
    if (.not. ok) error stop 'benchmark final enclosure failed'
    call streaming_radius(state, radius, ok)
    if (.not. ok) error stop 'benchmark radius failed'
    print '(a,es24.16)', 'velocity_interpolation_radius=', radius
    print '(a,3(es24.16,1x))', 'checked_sample_center_lower_upper=', value, lower, upper
    print '(i0,a)', cells, ' velocity cells; one harmonic; exact binary64 t=12.'
    print '(a)', 'Variant 1: primal nodal evolution. Variant 2: same + interval packet.'
    print '(a)', 'Warm reuse; observable sampling/rendering/encoding and derivatives excluded.'
    print '(a)', 'Cold setup and particle-density overhead: unavailable in this microbenchmark.'
    call execute_command_line('mkdir -p output')
    call result%create(run)
    call result%write('warm_seconds_per_result', elapsed, 's', &
        'per-trial warm time; variant 1=primal, 2=primal+interval packet', &
        dims='trial,variant', status='measurement')
    call result%write('velocity_interpolation_radius', radius, '1', &
        'velocity interpolation radius', status='certified-continuous streaming radius')
    call result%write('sample_center', value, '1', 'f(x=0,v=0.3125,t=12) center', &
        status='candidate center')
    call result%write('sample_lower', lower, '1', 'f sample lower bound', &
        status='certified-continuous enclosure')
    call result%write('sample_upper', upper, '1', 'f sample upper bound', &
        status='certified-continuous enclosure')
    call result%close()
end program benchmark_streaming_comparison_main