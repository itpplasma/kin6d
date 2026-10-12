program test_gs_solver_work
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortnum_krylov, only: KRYLOV_OK, KRYLOV_MAX_ITERATIONS
    use kin6d_grad_shafranov, only: gs_state_t, initialize_gs, solve_gs, apply_gs, &
        gs_linear_response
    implicit none
    type(gs_state_t) :: state
    real(dp), allocatable :: contour(:, :), input(:), action(:)
    logical, allocatable :: boundary(:)
    real(dp) :: h
    integer :: mesh, n, i, status, previous

    previous = 0
    do mesh = 1, 2
        n = 192*2**(mesh - 1)
        allocate (contour(2, n))
        do i = 1, n
            call curve(real(i - 1, dp)/n, contour(:, i))
        end do
        h = 2*acos(-1.0_dp)*2/n
        call initialize_gs(state, contour, h*h/2, status, 2, curve=curve)
        if (status /= 0) error stop 'shaped mesh'
        allocate (input(size(state%psi)), action(size(state%psi)))
        boundary = state%boundary
        state%boundary = .false.
        input = 1
        call apply_gs(state, input, action)
        if (any(action /= 0)) error stop 'constant flux has a spurious weak action'
        state%boundary = boundary
        call solve_gs(state, 1.0_dp, -0.08_dp, -0.7_dp, status)
        if (status /= 0) error stop 'shaped solve stagnation'
        print *, n, size(state%psi), state%iterations, state%residual
        if (state%residual > 1.0e-10_dp) error stop 'true residual target'
        if (state%iterations > 8*n) error stop 'unbounded shaped solver work'
        if (mesh > 1) then
            if (state%iterations > 3*previous) error stop 'solver work scaling'
        end if
        previous = state%iterations
        deallocate (contour, input, action, boundary)
    end do
    call check_capped_linear_response()
    call check_failed_solve_residual()
contains
    subroutine build_mesh(state, n)
        type(gs_state_t), intent(out) :: state
        integer, intent(in) :: n
        real(dp), allocatable :: contour(:, :)
        real(dp) :: h
        integer :: i, status
        allocate (contour(2, n))
        do i = 1, n
            call curve(real(i - 1, dp)/n, contour(:, i))
        end do
        h = 2*acos(-1.0_dp)*2/n
        call initialize_gs(state, contour, h*h/2, status, 2, curve=curve)
        if (status /= 0) error stop 'regression mesh'
    end subroutine

    subroutine check_capped_linear_response()
        ! Forces the capped CG path (MAX_ITER with stagnation_ok) and checks
        ! that the reported residual equals the independently recomputed one.
        type(gs_state_t) :: state
        real(dp), allocatable :: source(:), response(:), action(:)
        real(dp) :: residual_norm, expected
        integer :: i, status
        call build_mesh(state, 192)
        allocate (source(size(state%psi)), action(size(state%psi)))
        source = [(sin(0.37_dp*real(i, dp)), i=1, size(source))]
        call gs_linear_response(state, source, response, status, &
            relative_tolerance=1.0e-300_dp, stagnation_ok=.true., &
            residual_norm=residual_norm)
        if (status /= KRYLOV_MAX_ITERATIONS) error stop 'capped response status'
        if (.not. all(ieee_is_finite(response))) error stop 'capped response finite'
        call apply_gs(state, response, action)
        expected = norm2(source - action)
        if (.not. ieee_is_finite(residual_norm)) error stop 'residual_norm finite'
        if (abs(residual_norm - expected) > 1.0e-10_dp*expected) then
            print *, residual_norm, expected
            error stop 'residual_norm recomputed'
        end if
    end subroutine

    subroutine check_failed_solve_residual()
        ! With no Dirichlet boundary the stiffness operator is singular and the
        ! inconsistent source leaves CG without a solution, so the solve fails.
        ! The failed solve must still report a finite, non-sentinel residual.
        type(gs_state_t) :: state
        integer :: status
        call build_mesh(state, 192)
        state%boundary = .false.
        call solve_gs(state, 1.0_dp, -0.08_dp, -0.7_dp, status)
        if (status == KRYLOV_OK) error stop 'floating solve must not converge'
        if (.not. ieee_is_finite(state%residual)) error stop 'solve residual finite'
        if (state%residual >= huge(1.0_dp)) error stop 'solve residual not sentinel'
    end subroutine

    subroutine curve(t, point)
        real(dp), intent(in) :: t
        real(dp), intent(out) :: point(2)
        real(dp) :: angle
        angle = 2*acos(-1.0_dp)*t
        point = [6.2_dp + 2*cos(angle + asin(0.33_dp)*sin(angle)), &
            3.4_dp*sin(angle)]
    end subroutine
end program
