program test_gs_readback_work
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_grad_shafranov, only: gs_state_t, initialize_gs, solve_gs, &
        gs_axis_candidate
    use kin6d_gs_coarea, only: gs_level_coarea
    use kin6d_gs_coarea_leaves, only: gs_cubic_flux, gs_cubic_flux_jvp, &
        gs_cubic_flux_vjp
    implicit none
    type(gs_state_t) :: state
    real(dp), parameter :: r0 = 6.2_dp, a = 0.62_dp
    real(dp), allocatable :: contour(:, :)
    real(dp) :: h, axis(2), axis_flux, hessian(2, 2), levels(5), coarea(5)
    real(dp) :: work_per_contour(3), ratio, reference(5)
    integer :: mesh, n, j, status, element, work(5)
    call check_cubic_leaf()
    do mesh = 1, 3
        n = 48*2**(mesh - 1)
        allocate (contour(2, n))
        do j = 1, n
            call boundary(real(j - 1, dp)/n, contour(:, j))
        end do
        h = 2*acos(-1.0_dp)*a/n
        call initialize_gs(state, contour, h*h/2, status, 3, curve=boundary)
        if (status /= 0) error stop 'work mesh'
        call solve_gs(state, 1.0_dp, -16.0_dp, 0.0_dp, status)
        if (status /= 0) error stop 'work solve'
        call gs_axis_candidate(state, axis, axis_flux, hessian, element, status)
        if (status /= 0) error stop 'work axis'
        levels = axis_flux*[1.0_dp, 0.8_dp, 0.6_dp, 0.4_dp, 0.2_dp]
        call gs_level_coarea(state, levels, 80, coarea, axis, status, work=work)
        if (status /= 0) error stop 'work contours'
        if (mesh == 1) then
            call gs_level_coarea(state, levels, 80, reference, axis, status, &
                reference_cubic=.true.)
            if (status /= 0) error stop 'reference schedule'
            if (maxval(abs(coarea/reference - 1)) > 1.0e-11_dp) &
                error stop 'cached cubic differs from native basis'
        end if
        work_per_contour(mesh) = sum(real(work(2:), dp))/4
        print '(a,i4,a,f12.1)', 'n=', n, ' polynomial/location work per contour=', &
            work_per_contour(mesh)
        if (mesh > 1) then
            ratio = work_per_contour(mesh)/work_per_contour(mesh - 1)
            ! Doubling resolution doubles the number of crossed cells. It
            ! must not quadruple the number of edge roots or location trials.
            if (ratio > 2.6_dp) error stop 'readback work scales with mesh area'
        end if
        deallocate (contour)
    end do
contains
    subroutine check_cubic_leaf()
        real(dp) :: x(12), dx(12), xp(12), xm(12), dy(3), plus(3), minus(3)
        real(dp) :: dual(3), bar(12), step
        integer :: j
        x = [(sin(real(j, dp)), j=1,12)]
        dx = [(cos(real(j, dp)), j=1,12)]
        dual = [0.3_dp, -0.7_dp, 0.2_dp]
        step = 1.0e-5_dp
        xp = x + step*dx; xm = x - step*dx
        call gs_cubic_flux(xp(1), xp(2), xp(3), xp(4), xp(5), xp(6), xp(7), &
            xp(8), xp(9), xp(10), xp(11), xp(12), plus(1), plus(2), plus(3))
        call gs_cubic_flux(xm(1), xm(2), xm(3), xm(4), xm(5), xm(6), xm(7), &
            xm(8), xm(9), xm(10), xm(11), xm(12), minus(1), minus(2), minus(3))
        call gs_cubic_flux_jvp(x(1), x(2), x(3), x(4), x(5), x(6), x(7), x(8), &
            x(9), x(10), x(11), x(12), dx(1), dx(2), dx(3), dx(4), dx(5), &
            dx(6), dx(7), dx(8), dx(9), dx(10), dx(11), dx(12), dy(1), dy(2), dy(3))
        if (maxval(abs(dy - (plus - minus)/(2*step))) > 1.0e-8_dp) &
            error stop 'cubic cache JVP'
        call gs_cubic_flux_vjp(x(1), x(2), x(3), x(4), x(5), x(6), x(7), x(8), &
            x(9), x(10), x(11), x(12), dual(1), dual(2), dual(3), bar(1), bar(2), &
            bar(3), bar(4), bar(5), bar(6), bar(7), bar(8), bar(9), bar(10), &
            bar(11), bar(12))
        if (abs(dot_product(dy, dual) - dot_product(dx, bar)) > 1.0e-12_dp) &
            error stop 'cubic cache VJP pairing'
    end subroutine

    subroutine boundary(t, point)
        real(dp), intent(in) :: t
        real(dp), intent(out) :: point(2)
        real(dp) :: theta
        theta = 2*acos(-1.0_dp)*t
        point(1) = sqrt(r0*r0 + 2*r0*a*cos(theta))
        point(2) = r0*a*sin(theta)/point(1)
    end subroutine
end program
