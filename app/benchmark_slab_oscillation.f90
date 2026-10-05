program benchmark_slab_oscillation
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_slab_oscillation, only: slab_oscillation_state_t, &
         initialize_slab_oscillation, solve_slab_oscillation, slab_oscillation_radius, &
                                      sample_slab_oscillation
    use kin6d_slab_oscillation_leaves, only: slab_oscillation_period
    implicit none
    integer, parameter :: steps = 128, samples = 15, repetitions = 160
    type(slab_oscillation_state_t) :: state
    real(dp) :: duration, start, finish, candidate, complete, radius
    real(dp) :: q, p, qlo, qhi, plo, phi
    integer :: sample, repeat
    logical :: ok

    call slab_oscillation_period(2.0_dp, duration)
    call initialize_slab_oscillation(state, steps, duration, 1.0_dp, 0.0_dp, ok)
    if (.not. ok) error stop 'initialization failed'
    call slab_oscillation_radius(state, radius, ok)
    if (.not. ok .or. radius > 0.01_dp) error stop 'one-percent target failed'
    print '(a)', 'model=cold_planar_slab; tau=omega_pe*t; q=xi/xi0; p=dq/dtau'
    print '(a)', 'solver=FortNum_fixed_Cash_Karp_RK5(4); certificate=Hermite_residual'
    print '(a)', 'validation_outside_timing=true; setup_and_rendering_excluded=true'
    print '(a,i0)', 'steps=', steps
    print '(a,i0)', 'repetitions_per_sample=', repetitions
    print '(a,i0)', 'trace_payload_bytes=', 8*(5*steps + 4)
    print '(a)', 'sample,candidate_s,complete_s,certificate_extra_s,factor,radius'
    do sample = 1, 3
        call solve_slab_oscillation(state, ok)
        if (.not. ok) error stop 'warmup failed'
    end do
    do sample = 1, samples
        call cpu_time(start)
        do repeat = 1, repetitions
            call solve_slab_oscillation(state, ok, certify=.false.)
        end do
        call cpu_time(finish)
        candidate = (finish - start)/repetitions
        if (.not. ok) error stop 'candidate solve failed'
        call cpu_time(start)
        do repeat = 1, repetitions
            call solve_slab_oscillation(state, ok)
        end do
        call cpu_time(finish)
        complete = (finish - start)/repetitions
        if (.not. ok) error stop 'certificate failed'
        call slab_oscillation_radius(state, radius, ok)
        call sample_slab_oscillation(state, 0.0_dp, q, p, qlo, qhi, plo, phi, ok)
        if (.not. ok) error stop 'validated start unavailable'
        if (qlo > 1.0_dp .or. qhi < 1.0_dp .or. plo > 0.0_dp .or. phi < 0.0_dp) &
            error stop 'initial-value identity excluded'
        if (candidate <= 0.0_dp) error stop 'timer resolution insufficient'
        print '(i0,5(a,es24.16))', sample, ',', candidate, ',', complete, &
            ',', complete - candidate, ',', complete/candidate, ',', radius
    end do
end program benchmark_slab_oscillation
