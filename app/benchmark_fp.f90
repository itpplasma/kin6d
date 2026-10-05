program benchmark_fp_main
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_fp_laboratory, only: fp_state_t, initialize_fp, solve_fp, &
        fp_certificate_radius, fp_certificate_accepts, sample_fp
    use kin6d_result_io, only: run_t, result_file_t
    implicit none
    type(fp_state_t) :: state
    type(run_t) :: run
    type(result_file_t) :: result
    integer :: trials = 5, solves = 100, checks = 10
    integer :: cells = 4096
    real(dp), parameter :: exact_zero = 0.41795955023513457256435545185990697_dp
    real(dp) :: observable_scale = 1.0_dp, requested_radius = 0.01_dp
    real(dp), allocatable :: solve_time(:), check_time(:), accept_time(:)
    real(dp) :: start, finish, residual, radius, value, temporary
    integer :: trial, repeat, i, j, parse_status, unit, mid
    character(len=64) :: argument
    character(len=256) :: msg
    logical :: accepted
    namelist /benchmark_fp/ cells, trials, solves, checks, observable_scale, &
        requested_radius
    call run%start('benchmark_fp')
    parse_status = 0; msg = ''
    call run%open_config(unit)
    if (run%has_config()) read (unit, nml=benchmark_fp, iostat=parse_status, iomsg=msg)
    call run%close_config(unit, 'benchmark_fp', parse_status, msg)
    if (run%argument_count() > 1) error stop 'usage: benchmark_fp [--config FILE] [cells]'
    if (run%argument_count() == 1) then
        call run%argument(1, argument)
        read (argument, *, iostat=parse_status) cells
        if (parse_status /= 0) error stop 'cells must be an integer'
    end if
    if (cells < 8) error stop 'benchmark requires at least eight cells'
    if (mod(cells, 2) /= 0) error stop 'even cells required for f(0) node oracle'
    if (trials < 1 .or. solves < 1 .or. checks < 1) error stop 'invalid repetition counts'
    if (.not. (observable_scale > 0.0_dp .and. requested_radius > 0.0_dp)) &
        error stop 'invalid requested radius'
    call run%open_resolved(unit)
    write (unit, nml=benchmark_fp)
    call run%close_resolved(unit)
    allocate (solve_time(trials), check_time(trials), accept_time(trials))
    call initialize_fp(state, cells)
    call solve_fp(state, residual)
    radius = fp_certificate_radius(state)
    do trial = 1, trials
        call cpu_time(start)
        do repeat = 1, solves
            call solve_fp(state, residual)
        end do
        call cpu_time(finish)
        solve_time(trial) = (finish - start)/solves
        call cpu_time(start)
        do repeat = 1, checks
            radius = fp_certificate_radius(state)
        end do
        call cpu_time(finish)
        check_time(trial) = (finish - start)/checks
        ! Acceptance currently recomputes the certificate; measure that cost.
        call cpu_time(start)
        do repeat = 1, checks
            accepted = fp_certificate_accepts(state, radius)
        end do
        call cpu_time(finish)
        accept_time(trial) = (finish - start)/checks
        if (.not. accepted) error stop 'timed FP acceptance failed'
    end do
    ! Exact finite-domain erf fixture independently checks the large case.
    call sample_fp(state, 0.0_dp, value)
    if (abs(value - exact_zero) > radius) error stop 'large FP oracle failed'
    if (.not. fp_certificate_accepts(state, radius)) error stop 'large FP bound failed'
    if (fp_certificate_accepts(state, radius/2.0_dp)) &
        error stop 'undersized FP bound accepted'
    print '(a,i0,a,es24.16)', 'FP cells=', cells, ' achieved_sup_radius=', radius
    print '(a)', 'certificate_kind=analytic_reference_continuous_sup'
    print '(a)', 'arithmetic_mode=strict_interval derivative_status=not_benchmarked'
    print '(a,es24.16,a,es24.16)', 'observable_scale=', observable_scale, &
        ' requested_radius=', requested_radius
    print '(a,l1)', 'requested_radius_met=', radius <= requested_radius
    print '(a,es24.16,a,es24.16)', 'observable f(0)=', value, ' flux_residual=', residual
    print '(a,i0)', 'state_bytes=', (size(state%v) + size(state%f))*storage_size(value)/8
    print '(a)', 'trial,solve_cpu_seconds,certificate_cpu_seconds,'// &
        'acceptance_cpu_seconds,total_cpu_seconds'
    do trial = 1, trials
        print '(i0,4(",",es24.16))', trial, solve_time(trial), check_time(trial), &
            accept_time(trial), solve_time(trial) + check_time(trial) + accept_time(trial)
    end do
    ! Sorting five timer samples is instrumentation, outside the physics kernels.
    do i = 2, trials
        do j = i, 2, -1
            if (solve_time(j) < solve_time(j - 1)) then
                temporary = solve_time(j)
                solve_time(j) = solve_time(j - 1); solve_time(j - 1) = temporary
            end if
            if (check_time(j) < check_time(j - 1)) then
                temporary = check_time(j)
                check_time(j) = check_time(j - 1); check_time(j - 1) = temporary
            end if
            if (accept_time(j) < accept_time(j - 1)) then
                temporary = accept_time(j)
                accept_time(j) = accept_time(j - 1); accept_time(j - 1) = temporary
            end if
        end do
    end do
    mid = (trials + 1)/2
    print '(a,3es16.7)', 'solve median/min/max=', solve_time(mid), &
        solve_time(1), solve_time(trials)
    print '(a,3es16.7)', 'certificate median/min/max=', check_time(mid), &
        check_time(1), check_time(trials)
    print '(a,3es16.7)', 'acceptance median/min/max=', accept_time(mid), &
        accept_time(1), accept_time(trials)
    call execute_command_line('mkdir -p output')
    call result%create(run)
    call result%write('achieved_sup_radius', radius, '1', &
        'continuous sup-norm radius of the FP profile', &
        status='certified-continuous analytic reference sup bound (strict interval)')
    call result%write('f0', value, '1', 'candidate observable f(0)', &
        status='candidate; |f(0)-exact| <= achieved_sup_radius checked')
    call result%write('f0_exact', exact_zero, '1', 'exact finite-domain erf oracle f(0)')
    call result%write('flux_residual', residual, '1', 'final flux residual', &
        status='diagnostic')
    call result%write('solve_cpu_seconds', solve_time, 's', &
        'sorted per-solve CPU time', dims='trial', status='measurement')
    call result%write('certificate_cpu_seconds', check_time, 's', &
        'sorted certificate CPU time', dims='trial', status='measurement')
    call result%write('acceptance_cpu_seconds', accept_time, 's', &
        'sorted acceptance CPU time', dims='trial', status='measurement')
    call result%close()
end program benchmark_fp_main