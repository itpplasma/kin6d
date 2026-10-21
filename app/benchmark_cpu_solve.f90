program benchmark_cpu_solve
    ! Prototype consumer of the owning fwprec API; no private solver.
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use kin6d_grad_shafranov, only: gs_state_t, initialize_gs, solve_gs, apply_gs
    use kin6d_gs_csr, only: gs_csr_spd_matrix
    use fortnum_precision, only: precision_policy_t, precision_pair_t, &
                            precision_preset, FP64, FP32, FP16, BF16, STAGE_PRECONDITION
    use fortsparse_precision_csr, only: precision_csr_t, csr_gmres_ir
    implicit none
    type(gs_state_t) :: state
    type(precision_csr_t) :: matrix
    type(precision_policy_t) :: policy
    type(precision_pair_t) :: final_precision
    integer, parameter :: formats(4) = [FP64, FP32, FP16, BF16]
    integer, allocatable :: row(:), column(:)
    real(dp), allocatable :: values(:), rhs(:), x(:), ax(:), contour(:, :)
    real(dp) :: theta, h, t, seconds, elapsed, residual, true_residual, error, scale
    integer :: nr, i, f, status, info, steps, rep, repeats
    logical :: ok
    character(32) :: arg
    nr = 24; repeats = 3
    if (command_argument_count() > 0) then
        call get_command_argument(1, arg); read (arg, *) nr
    end if
    if (command_argument_count() > 1) then
        call get_command_argument(2, arg); read (arg, *) repeats
    end if
    if (repeats < 1) error stop 'repeats >= 1'
    if (nr < 8) error stop 'nr >= 8'
    allocate (contour(2, nr))
    do i = 1, nr
        theta = 2*acos(-1.0_dp)*(i - 1)/nr
        contour(:, i) = [3 + cos(theta), sin(theta)]
    end do
    h = 2*acos(-1.0_dp)/nr
    call initialize_gs(state, contour, h*h/2, status, 2)
    if (status /= 0) error stop 'GS mesh'
    call solve_gs(state, 1.0_dp, -0.1_dp, 0.0_dp, status)
    if (status /= 0) error stop 'native GS reference'
    allocate (rhs(size(state%psi)), x(size(state%psi)), ax(size(state%psi)))
    ! Independent exchange-form production operator supplies RHS and checks
    ! the solve; the mixed solver receives an ordinary masked CSR matrix.
    call apply_gs(state, state%psi, rhs)
    call gs_csr_spd_matrix(state%op, state%boundary, row, column, values)
    call matrix%store(size(x), row, column, values, FP64, ok)
    if (.not. ok) error stop 'original CSR storage'
    scale = maxval(abs(state%psi))
    print '(a)', '# [V] discrete GS comparison; continuous certificate unavailable'
print '(a)', '# n,storage,compute,seconds,steps,residual,exchange_residual,state_defect'
    do f = 1, 4
        seconds = huge(1.0_dp)
        do rep = 1, repeats
            policy = precision_preset('cpu-mixed', ok)
            if (.not. ok) error stop 'precision preset'
            call policy%set(STAGE_PRECONDITION, formats(f), &
                            merge(FP64, FP32, formats(f) == FP64), ok)
            if (.not. ok) error stop 'precision override'
            x = 0
            t = wall()
            call csr_gmres_ir(matrix, rhs, x, policy, 1e-5_dp, 30, 200, 40, &
                              info, steps, residual)
            elapsed = wall() - t
            seconds = min(seconds, elapsed)
            if (info /= 0) error stop 'shared GMRES-IR failed'
        end do
        call apply_gs(state, x, ax)
        true_residual = norm2(rhs - ax)/norm2(rhs)
        error = maxval(abs(x - state%psi))/scale
        final_precision = policy%get(STAGE_PRECONDITION)
        print '(i0,2(a,i0),a,es23.15,a,i0,3(a,es23.15))', size(x), ',', &
            final_precision%storage, ',', final_precision%compute, ',', seconds, &
            ',', steps, ',', residual, ',', true_residual, ',', error
        if (info /= 0) error stop 'shared GMRES-IR failed'
        if (true_residual > 1.1e-5_dp) error stop 'independent GS residual'
        if (error > 1e-2_dp) error stop 'GS state error exceeds 1%'
    end do
contains
    function wall() result(t)
        real(dp) :: t
        integer(int64) :: count, rate
        call system_clock(count, rate)
        t = real(count, dp)/real(rate, dp)
    end function
end program benchmark_cpu_solve
