program test_gs_nonlinear
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use kin6d_grad_shafranov, only: gs_state_t, initialize_gs, &
                                    gs_cell_data, gs_element_jet
    use kin6d_gs_majorant, only: gs_functional_majorant
    use kin6d_gs_profile_law, only: gs_profile_law_t
    use kin6d_gs_nonlinear, only: solve_gs_cubic_profiles
    use fortfem_triangle_duffy_quadrature, only: triangle_duffy_quadrature
    implicit none
    real(dp), parameter :: r0 = 6.2_dp, f = 32.86_dp
    real(dp) :: a, level
    real(dp), parameter :: c = f/(12*r0**3)
    type(gs_state_t) :: state
    type(gs_profile_law_t) :: law, absolute_law
    real(dp), allocatable :: boundary(:, :), x(:), y(:), w(:)
    real(dp) :: errors(3, 3), rates(2), point(2), det, shape(10)
    real(dp) :: gradients(2, 10), hessians(2, 2, 10), grad(2), exact_grad(2)
    real(dp) :: psi, exact, h, elapsed, physical_h(2, 2), value, read_grad(2)
    real(dp) :: s, base, majorant, flux_defect, balance_defect, effectivity
    integer :: counts(3)
    integer :: mesh, n, j, e, q, status, degree, ndof
    character(8) :: argument
    integer(int64) :: start, finish, clock_rate
    degree = 3
    if (command_argument_count() > 0) then
        call get_command_argument(1, argument)
        read (argument, *) degree
    end if
    a = r0/3
    level = 4*r0*r0*a*a
    law%knots = [0.0_dp, 0.31_dp, 0.77_dp, 1.0_dp]
    allocate (law%pp(4, 3), law%ff(4, 3))
    law%normalized = .true.; law%axis = -c*level; law%p_edge = 0; law%f_edge = f
    base = -16*c
    do j = 1, 3
        s = law%knots(j)
        law%pp(:, j) = base*[1 + 0.2_dp*s*s, 0.4_dp*s, 0.2_dp, 0.0_dp]
    end do
    law%ff = 0
    ndof = (degree + 1)*(degree + 2)/2
    call triangle_duffy_quadrature(10, x, y, w, status)
    if (status /= 0) error stop 'quadrature'
    print '(a)', 'n dof wall_s psi_L2 Bpol_L2 energy'
    do mesh = 1, 3
        n = 32*2**(mesh - 1)
        allocate (boundary(2, n))
        call system_clock(start, clock_rate)
        do j = 1, n
            call exact_boundary(real(j - 1, dp)/n, boundary(:, j))
        end do
        h = 2*acos(-1.0_dp)*a/n
        call initialize_gs(state, boundary, h*h/2, status, degree, &
                           curve=exact_boundary)
        if (status /= 0) error stop 'mesh initialization'
        call solve_gs_cubic_profiles(state, 1.0_dp, law, status, &
                             relative_tolerance=1.0e-10_dp, forcing=manufactured_source)
        print *, 'nonlinear status/iterations/residual/axis:', status, &
            state%iterations, state%residual, law%axis
        if (status /= 0) error stop 'nonlinear solve'
        counts(mesh) = state%iterations
        if (counts(mesh) > 8) error stop 'Newton iteration bound'
        call system_clock(finish)
        elapsed = real(finish - start, dp)/clock_rate
        if (elapsed > 300) error stop 'solve time limit'
        errors(:, mesh) = 0
        do e = 1, size(state%triangles, 2)
            do q = 1, size(w)
                call gs_cell_data(state, e, [x(q), y(q)], point, det, &
                       shape(:ndof), gradients(:, :ndof), hessians(:, :, :ndof), status)
                if (status /= 0) error stop 'physical cell'
                psi = dot_product(shape(:ndof), state%psi(state%triangles(:, e)))
                grad = matmul(gradients(:, :ndof), &
                    state%psi(state%triangles(:, e)) - state%psi(state%triangles(1, e)))
                exact = c*((point(1)**2 - r0*r0)**2 &
                           + 4*point(1)**2*point(2)**2 - level)
                exact_grad = [4*c*point(1)*(point(1)**2 - r0*r0 &
                                            + 2*point(2)**2), 8*c*point(1)**2*point(2)]
                errors(1, mesh) = errors(1, mesh) &
                                  + w(q)*det*(psi - exact)**2
                errors(2, mesh) = errors(2, mesh) &
                                  + w(q)*det*sum((grad - exact_grad)**2)/point(1)**2
                errors(3, mesh) = errors(3, mesh) &
                                  + w(q)*det*sum((grad - exact_grad)**2)/point(1)
                if (q /= 1) cycle
                call gs_element_jet(state, e, point(1), point(2), value, &
                                    read_grad, physical_h, status)
                if (status /= 0) error stop 'physical readback'
                if (abs(value - psi) > 1.0e-11_dp .or. &
                    norm2(read_grad - grad) > 1.0e-10_dp) error stop 'readback jet'
            end do
        end do
        errors(:, mesh) = sqrt(errors(:, mesh))
        print '(i4,1x,i7,4(1x,es13.5))', n, &
            count(.not. state%boundary), elapsed, errors(:, mesh)
        if (mesh == 1) then
            ! The exact manufactured problem also has an absolute-flux law
            ! with nonpositive source derivative: alpha=1 by monotonicity.
            absolute_law = law
            absolute_law%normalized = .false.
            absolute_law%knots = -c*level*(1 - law%knots)
            do j = 1, 4
                absolute_law%pp(j, :) = law%pp(j, :)/(c*level)**(j - 1)
            end do
            call gs_functional_majorant(state, 1.0_dp, 0.0_dp, 0.0_dp, &
                          majorant, flux_defect, balance_defect, status, law=absolute_law, &
                              stability_lower_bound=1.0_dp, forcing=manufactured_source)
            if (status /= 0) error stop 'nonlinear recovered-flux majorant'
            effectivity = majorant/errors(3, mesh)
            print *, 'nonlinear monotone majorant, effectivity:', majorant, effectivity
            if (effectivity < 1 .or. effectivity > 3) error stop 'nonlinear effectivity'
        end if
        deallocate (boundary)
    end do
    do mesh = 2, 3
        rates = log(errors(1:2, mesh - 1)/errors(1:2, mesh))/log(2.0_dp)
        print '(a,i1,2f10.5)', 'rates pair=', mesh - 1, rates
        if (abs(rates(1) - (degree + 1)) > 0.5_dp) &
            error stop 'curved psi L2 rate'
        if (abs(rates(2) - degree) > 0.5_dp) error stop 'curved B L2 rate'
    end do
    if (maxval(counts) - minval(counts) > 2) error stop 'Newton mesh dependence'
contains
    function manufactured_source(r, z) result(source)
        real(dp), intent(in) :: r, z
        real(dp) :: source, exact_s
        exact_s = ((r*r - r0*r0)**2 + 4*r*r*z*z)/level
        source = 16*c*r*0.2_dp*exact_s**2
    end function
    subroutine exact_boundary(t, point)
        real(dp), intent(in) :: t
        real(dp), intent(out) :: point(2)
        real(dp) :: angle
        angle = 2*acos(-1.0_dp)*t
        point(1) = sqrt(r0*r0 + 2*r0*a*cos(angle))
        point(2) = r0*a*sin(angle)/point(1)
    end subroutine
end program
