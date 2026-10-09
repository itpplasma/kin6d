program test_gs_solver_work
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_grad_shafranov, only: gs_state_t, initialize_gs, solve_gs, apply_gs
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
contains
    subroutine curve(t, point)
        real(dp), intent(in) :: t
        real(dp), intent(out) :: point(2)
        real(dp) :: angle
        angle = 2*acos(-1.0_dp)*t
        point = [6.2_dp + 2*cos(angle + asin(0.33_dp)*sin(angle)), &
            3.4_dp*sin(angle)]
    end subroutine
end program
