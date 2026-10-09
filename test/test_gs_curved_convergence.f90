program test_gs_curved_convergence
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use kin6d_grad_shafranov, only: gs_state_t, initialize_gs, solve_gs, &
        gs_cell_data, gs_element_jet
    use kin6d_gs_majorant, only: gs_functional_majorant
    use fortfem_triangle_duffy_quadrature, only: triangle_duffy_quadrature
    implicit none
    real(dp), parameter :: r0 = 6.2_dp, a = 0.62_dp, f = 32.86_dp
    real(dp), parameter :: c = f/(12*r0**3), level = 4*r0*r0*a*a
    type(gs_state_t) :: state
    real(dp), allocatable :: boundary(:, :), x(:), y(:), w(:)
    real(dp) :: errors(3, 3, 2), rates(2), point(2), det, shape(6)
    real(dp) :: gradients(2, 6), hessians(2, 2, 6), grad(2), exact_grad(2)
    real(dp) :: psi, exact, h, elapsed, physical_h(2, 2), value, read_grad(2)
    real(dp) :: majorant, flux_defect, balance_defect, effectivity
    integer :: mode, mesh, n, j, e, q, status
    integer(int64) :: start, finish, clock_rate
    call triangle_duffy_quadrature(10, x, y, w, status)
    if (status /= 0) error stop 'quadrature'
    print '(a)', 'geometry n dof wall_s psi_L2 Bpol_L2 energy'
    do mode = 1, 2
        do mesh = 1, 3
            n = 48*2**(mesh - 1)
            allocate (boundary(2, n))
            call system_clock(start, clock_rate)
            do j = 1, n
                call exact_boundary(real(j - 1, dp)/n, boundary(:, j))
            end do
            h = 2*acos(-1.0_dp)*a/n
            if (mode == 1) then
                call initialize_gs(state, boundary, h*h/2, status, 2, &
                    curve=exact_boundary)
            else
                call initialize_gs(state, boundary, h*h/2, status, 2)
            end if
            if (status /= 0) error stop 'mesh initialization'
            call solve_gs(state, 1.0_dp, -16*c, 0.0_dp, status)
            if (status /= 0) error stop 'linear solve'
            call system_clock(finish)
            elapsed = real(finish - start, dp)/clock_rate
            if (elapsed > 300) error stop 'solve time limit'
            errors(:, mesh, mode) = 0
            do e = 1, size(state%triangles, 2)
                do q = 1, size(w)
                    call gs_cell_data(state, e, [x(q), y(q)], point, det, shape, &
                        gradients, hessians, status)
                    if (status /= 0) error stop 'physical cell'
                    psi = dot_product(shape, state%psi(state%triangles(:, e)))
                    grad = matmul(gradients, state%psi(state%triangles(:, e)) &
                        - state%psi(state%triangles(1, e)))
                    exact = c*((point(1)**2 - r0*r0)**2 &
                        + 4*point(1)**2*point(2)**2 - level)
                    exact_grad = [4*c*point(1)*(point(1)**2 - r0*r0 &
                        + 2*point(2)**2), 8*c*point(1)**2*point(2)]
                    errors(1, mesh, mode) = errors(1, mesh, mode) &
                        + w(q)*det*(psi - exact)**2
                    errors(2, mesh, mode) = errors(2, mesh, mode) &
                        + w(q)*det*sum((grad - exact_grad)**2)/point(1)**2
                    errors(3, mesh, mode) = errors(3, mesh, mode) &
                        + w(q)*det*sum((grad - exact_grad)**2)/point(1)
                    if (q /= 1) cycle
                    call gs_element_jet(state, e, point(1), point(2), value, &
                        read_grad, physical_h, status)
                    if (status /= 0) error stop 'physical readback'
                    if (abs(value - psi) > 1.0e-11_dp .or. &
                        norm2(read_grad - grad) > 1.0e-10_dp) error stop 'readback jet'
                end do
            end do
            errors(:, mesh, mode) = sqrt(errors(:, mesh, mode))
            print '(i1,1x,i4,1x,i7,4(1x,es13.5))', mode, n, &
                count(.not. state%boundary), elapsed, errors(:, mesh, mode)
            if (mode == 1) then
                call gs_functional_majorant(state, 1.0_dp, -16*c, 0.0_dp, &
                    majorant, flux_defect, balance_defect, status)
                if (status /= 0) error stop 'functional majorant recovery'
                effectivity = majorant/errors(3, mesh, mode)
                print '(a,4es14.5)', 'majorant, flux, balance, effectivity: ', &
                    majorant, flux_defect, balance_defect, effectivity
                if (effectivity < 1 .or. effectivity > 3) error stop 'effectivity'
                if (mesh == 1) then
                    state%psi = 0
                    call gs_functional_majorant(state, 1.0_dp, 0.0_dp, 0.0_dp, &
                        majorant, flux_defect, balance_defect, status)
                    if (status /= 0 .or. majorant > tiny(1.0_dp)) &
                        error stop 'zero-source exact majorant'
                end if
            end if
            deallocate (boundary)
        end do
        do mesh = 2, 3
            rates = log(errors(1:2, mesh - 1, mode)/errors(1:2, mesh, mode))/log(2.0_dp)
            print '(a,i1,a,i1,2f10.5)', 'rates geometry=', mode, ' pair=', mesh - 1, rates
            if (mode == 1) then
                if (abs(rates(1) - 3) > 0.3_dp) error stop 'curved psi L2 rate'
                if (abs(rates(2) - 2) > 0.3_dp) error stop 'curved B L2 rate'
            else
                if (rates(1) > 2.5_dp) error stop 'affine boundary degradation missing'
                if (rates(2) > 1.85_dp) error stop 'affine field degradation missing'
            end if
        end do
    end do
contains
    subroutine exact_boundary(t, point)
        real(dp), intent(in) :: t
        real(dp), intent(out) :: point(2)
        real(dp) :: angle
        angle = 2*acos(-1.0_dp)*t
        point(1) = sqrt(r0*r0 + 2*r0*a*cos(angle))
        point(2) = r0*a*sin(angle)/point(1)
    end subroutine
end program
