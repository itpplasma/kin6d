program test_gs_native_adapter
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_gs_direct, only: gs_direct_t, gs_direct_available, &
        gs_direct_factor, gs_direct_solve, gs_direct_solve_block, &
        gs_direct_free, GS_DIRECT_INDEFINITE
    use fortsparse, only: FORTSPARSE_BACKEND_NATIVE
    implicit none
    type(gs_direct_t) :: factor
    integer, parameter :: ptr(4) = [1, 3, 6, 8]
    integer, parameter :: col(7) = [1, 2, 1, 2, 3, 2, 3]
    real(dp), parameter :: val(7) = [4.0_dp, 1.0_dp, 1.0_dp, 3.0_dp, &
        1.0_dp, 1.0_dp, 2.0_dp]
    real(dp) :: x(3), rhs(3, 2), block(3, 2), saddle(2), saved(2)
    logical :: ok

    if (.not. gs_direct_available()) error stop 'native backend required'
    ! Unsupported device/symmetry requests fail before touching any GPU.
    call gs_direct_factor(factor, 3, ptr, col, val, ok, device=-2)
    if (ok .or. factor%factored .or. .not. factor%failed) &
        error stop 'invalid device was accepted'
    call gs_direct_factor(factor, 2, [1, 3, 5], [1, 2, 1, 2], &
        [0.0_dp, 1.0_dp, 1.0_dp, 0.0_dp], ok, &
        symmetry=GS_DIRECT_INDEFINITE, device=0)
    if (ok .or. factor%factored .or. .not. factor%failed) &
        error stop 'unsupported CUDA indefinite factor was accepted'
    call gs_direct_factor(factor, 3, ptr, col, val, ok)
    if (.not. ok) error stop 'SPD factor'
    if (factor%stats%active_device /= -1) error stop 'default CPU device changed'
    if (factor%solver%backend_id /= FORTSPARSE_BACKEND_NATIVE) &
        error stop 'explicit native backend'
    call gs_direct_solve(factor, [5.0_dp, 5.0_dp, 3.0_dp], x, ok)
    if (.not. ok) error stop 'SPD solve'
    if (maxval(abs(x - 1.0_dp)) > 1.0e-13_dp) error stop 'SPD known solution'
    rhs(:, 1) = [5.0_dp, 5.0_dp, 3.0_dp]
    rhs(:, 2) = [2.0_dp, -4.0_dp, 0.0_dp]
    call gs_direct_solve_block(factor, rhs, block, ok)
    if (.not. ok) error stop 'retained block solve'
    if (maxval(abs(block(:, 1) - 1.0_dp)) > 1.0e-13_dp) error stop 'block one'
    if (maxval(abs(block(:, 2) - [1.0_dp, -2.0_dp, 1.0_dp])) &
        > 1.0e-13_dp) error stop 'block two'

    ! A zero diagonal needs an indefinite 2-by-2 pivot, and must fail SPD.
    call gs_direct_factor(factor, 2, [1, 3, 5], [1, 2, 1, 2], &
        [0.0_dp, 1.0_dp, 1.0_dp, 0.0_dp], ok)
    if (ok .or. factor%factored .or. .not. factor%failed) &
        error stop 'SPD rejection status'
    call gs_direct_factor(factor, 2, [1, 3, 5], [1, 2, 1, 2], &
        [0.0_dp, 1.0_dp, 1.0_dp, 0.0_dp], ok, GS_DIRECT_INDEFINITE)
    if (.not. ok) error stop 'indefinite recovery factor'
    saved = [2.0_dp, -3.0_dp]
    call gs_direct_solve(factor, [-3.0_dp, 2.0_dp], saddle, ok)
    if (.not. ok) error stop 'indefinite solve'
    if (maxval(abs(saddle - saved)) > 1.0e-13_dp) &
        error stop 'indefinite known solution'
    call gs_direct_free(factor)
    call gs_direct_free(factor)
    if (factor%factored .or. factor%failed) error stop 'free state'
    call gs_direct_solve(factor, [1.0_dp, 1.0_dp], saddle, ok)
    if (ok) error stop 'solve after free'
end program
