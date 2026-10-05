program benchmark_streaming_comparison
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_streaming_comparison, only: streaming_state_t, init_streaming, &
        advance_streaming, evaluate_streaming, streaming_radius
    use kin6d_streaming_leaves, only: streaming_velocity_coordinate
    implicit none
    integer, parameter :: cells = 128, repeats = 50, trials = 5
    type(streaming_state_t) :: state
    real(dp) :: velocity(cells + 1), a0(0:1, cells + 1), b0(0:1, cells + 1)
    real(dp) :: start, finish, elapsed(trials, 2), sorted(trials, 2), scratch
    real(dp) :: value, lower, upper, radius
    integer :: j, trial, order, variant, repeat, k
    logical :: ok

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
    print '(a)', 'variant,median_seconds,min_seconds,max_seconds'
    do variant = 1, 2
        print '(i0,3(a,es24.16))', variant, ',', sorted(3, variant), ',', &
            sorted(1, variant), ',', sorted(5, variant)
    end do
    print '(a,es24.16)', 'extra_certificate_seconds=', &
        sorted(3, 2) - sorted(3, 1)
    print '(a,f12.4)', 'certified_evolution_over_solve_only=', &
        sorted(3, 2)/sorted(3, 1)
    call advance_streaming(state, 12.0_dp, ok)
    if (.not. ok) error stop 'benchmark final certificate failed'
    call evaluate_streaming(state, 0.0_dp, 0.3125_dp, value, lower, upper, ok)
    if (.not. ok) error stop 'benchmark final enclosure failed'
    call streaming_radius(state, radius, ok)
    if (.not. ok) error stop 'benchmark radius failed'
    print '(a,es24.16)', 'velocity_interpolation_radius=', radius
    print '(a,3(es24.16,1x))', 'checked_sample_center_lower_upper=', value, lower, upper
    print '(a)', '128 velocity cells; one harmonic; exact binary64 t=12.'
    print '(a)', 'Variant 1: primal nodal evolution. Variant 2: same + interval packet.'
    print '(a)', 'Warm reuse; observable sampling/rendering/encoding and derivatives excluded.'
    print '(a)', 'Cold setup and particle-density overhead: unavailable in this microbenchmark.'
end program benchmark_streaming_comparison
