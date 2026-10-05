program benchmark_slab_oscillation_main
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_slab_oscillation, only: slab_oscillation_state_t, &
         initialize_slab_oscillation, solve_slab_oscillation, slab_oscillation_radius, &
                                      sample_slab_oscillation
    use kin6d_slab_oscillation_leaves, only: slab_oscillation_period
    use kin6d_result_io, only: run_t, result_file_t
    implicit none
    integer :: steps = 128, samples = 15, repetitions = 160
    type(run_t) :: run
    type(result_file_t) :: result
    real(dp), allocatable :: candidate_s(:), complete_s(:), radius_s(:)
    integer :: unit, ios
    character(len=256) :: msg
    namelist /benchmark_slab_oscillation/ steps, samples, repetitions
    type(slab_oscillation_state_t) :: state
    real(dp) :: duration, start, finish, candidate, complete, radius
    real(dp) :: q, p, qlo, qhi, plo, phi
    integer :: sample, repeat
    logical :: ok

    call run%start('benchmark_slab_oscillation')
    if (run%argument_count() > 0) &
        error stop 'usage: benchmark_slab_oscillation [--config FILE]'
    ios = 0; msg = ''
    call run%open_config(unit)
    if (run%has_config()) read (unit, nml=benchmark_slab_oscillation, iostat=ios, &
        iomsg=msg)
    call run%close_config(unit, 'benchmark_slab_oscillation', ios, msg)
    if (steps < 1 .or. samples < 1 .or. repetitions < 1) &
        error stop 'invalid benchmark_slab_oscillation sizes'
    call run%open_resolved(unit)
    write (unit, nml=benchmark_slab_oscillation)
    call run%close_resolved(unit)
    allocate (candidate_s(samples), complete_s(samples), radius_s(samples))
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
        candidate_s(sample) = candidate; complete_s(sample) = complete
        radius_s(sample) = radius
    end do
    call execute_command_line('mkdir -p output')
    call result%create(run)
    call result%write('candidate_s', candidate_s, 's', &
        'candidate solve CPU time per solve', dims='sample', status='measurement')
    call result%write('complete_s', complete_s, 's', &
        'candidate plus Hermite residual certificate CPU time per solve', &
        dims='sample', status='measurement')
    call result%write('radius', radius_s, '1', 'certified trajectory radius', &
        dims='sample', status='certified Hermite residual radius (<= 0.01 required)')
    call result%close()
end program benchmark_slab_oscillation_main