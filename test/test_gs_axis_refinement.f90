program test_gs_axis_refinement
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_grad_shafranov, only: gs_state_t, initialize_gs, solve_gs
    use kin6d_gs_axis_recovery, only: gs_magnetic_axis
    implicit none
    real(dp), parameter :: r0 = 6.2_dp, a = r0/3, c = 32.86_dp/(12*r0**3)
    type(gs_state_t) :: state
    real(dp), allocatable :: contour(:, :)
    real(dp) :: axis(2), error(3), rate, h
    integer :: mesh, n, i, status
    do mesh = 1, 3
        n = 96*2**(mesh - 1)
        allocate (contour(2, n))
        do i = 1, n
            call boundary(real(i - 1, dp)/n, contour(:, i))
        end do
        h = 2*acos(-1.0_dp)*a/n
        call initialize_gs(state, contour, h*h/2, status, 2, curve=boundary)
        if (status /= 0) error stop 'axis mesh'
        call solve_gs(state, 1.0_dp, -16*c, 0.0_dp, status)
        if (status /= 0) error stop 'axis solve'
        call gs_magnetic_axis(state, axis, status)
        if (status /= 0) error stop 'axis candidate'
        error(mesh) = norm2(axis - [r0, 0.0_dp])/r0
        print *, n, error(mesh)
        deallocate (contour)
    end do
    do mesh = 2, 3
        rate = log(error(mesh - 1)/error(mesh))/log(2.0_dp)
        print *, 'axis rate', rate
        if (rate < 1.7_dp) error stop 'axis below field order'
    end do
contains
    subroutine boundary(t, point)
        real(dp), intent(in) :: t
        real(dp), intent(out) :: point(2)
        real(dp) :: angle
        angle = 2*acos(-1.0_dp)*t
        point(1) = sqrt(r0*r0 + 2*r0*a*cos(angle))
        point(2) = r0*a*sin(angle)/point(1)
    end subroutine
end program
