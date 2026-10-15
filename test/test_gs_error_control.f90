program test_gs_error_control
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_grad_shafranov, only: gs_state_t, initialize_gs, &
                                    gs_linear_response, gs_cell_data
    use kin6d_gs_profile_law, only: gs_profile_law_t
    use kin6d_gs_error_control, only: gs_stability_t, gs_linear_goal_t, &
                  gs_fixed_stability, gs_fixed_source_stability, gs_fixed_linear_goal, &
                              gs_state_error_bound, GS_ERROR_OK, GS_ERROR_UNAVAILABLE, &
                                      GS_ERROR_NONCOERCIVE
    use fortfem_triangle_duffy_quadrature, only: triangle_duffy_quadrature
    implicit none
    type(gs_state_t) :: state
    type(gs_profile_law_t) :: law
    type(gs_stability_t) :: stability
    type(gs_linear_goal_t) :: goal, perturbed, reversed
    real(dp), parameter :: contour(2, 4) = reshape([2.0_dp, -1.0_dp, &
                               4.0_dp, -1.0_dp, 4.0_dp, 1.0_dp, 2.0_dp, 1.0_dp], [2, 4])
    real(dp), allocatable :: rhs(:), psi(:)
    real(dp) :: cf2, expected, radius
    integer :: status

    call initialize_gs(state, contour, 0.08_dp, status, 2)
    if (status /= 0) error stop 'mesh'
    call gs_fixed_source_stability(state, stability, status)
    if (status /= GS_ERROR_OK) error stop 'fixed source stability'
    if (abs(stability%alpha - 1) > 1e-14_dp) error stop 'fixed source alpha'
    cf2 = 8.0_dp/acos(-1.0_dp)**2
    if (abs(stability%friedrichs_squared - cf2) > 1e-13_dp) &
        error stop 'continuous rectangle Friedrichs bound'
    if (stability%arithmetic_rigorous) error stop 'unchecked arithmetic promoted'
    allocate (law%knots(2), law%pp(4, 1), law%ff(4, 1))
    law%knots = [-1.0_dp, 1.0_dp]
    law%pp = 0; law%ff = 0; law%pp(2, 1) = 0.01_dp
    call gs_fixed_stability(state, law, 1.0_dp, stability, status)
    expected = 1 - 0.04_dp*cf2
    if (status /= GS_ERROR_OK) error stop 'absolute affine stability'
    if (abs(stability%alpha - expected) > 1e-13_dp) &
        error stop 'independent analytical coercivity'
    call gs_state_error_bound(0.2_dp, stability, radius, status)
    if (status /= GS_ERROR_OK) error stop 'state error composition'
    if (abs(radius - 0.2_dp/expected) > 1e-14_dp) error stop 'state radius'
    law%pp(2, 1) = 10
    call gs_fixed_stability(state, law, 1.0_dp, stability, status)
    if (status /= GS_ERROR_NONCOERCIVE .or. stability%available) &
        error stop 'noncoercive profile accepted'
    call gs_state_error_bound(0.2_dp, stability, radius, status)
    if (status /= GS_ERROR_UNAVAILABLE) error stop 'invalid state band'
    law%pp = 0; law%pp(3, 1) = 0.01_dp
    call gs_fixed_stability(state, law, 1.0_dp, stability, status)
    if (status /= GS_ERROR_UNAVAILABLE) error stop 'unbounded extension accepted'
    law%pp = 0
    law%pp(:, 1) = [0.0_dp, -0.1_dp, 0.3_dp, -0.2_dp]
    call gs_fixed_stability(state, law, 1.0_dp, stability, status)
    if (status /= GS_ERROR_OK) error stop 'globally bounded cubic slope refused'
    expected = 1 - 4*0.05_dp*cf2
    if (abs(stability%alpha - expected) > 1e-13_dp) &
        error stop 'cubic derivative interior maximum missed'
    law%pp(:, 1) = [0.0_dp, -0.1_dp, 0.0_dp, -0.1_dp]
    call gs_fixed_stability(state, law, 1.0_dp, stability, status)
    if (status /= GS_ERROR_OK) error stop 'dissipative cubic refused'
    if (abs(stability%alpha - 1) > 1e-13_dp) &
        error stop 'dissipative slope treated as positive Lipschitz'
    law%pp = 0; law%normalized = .true.; law%knots = [0.0_dp, 1.0_dp]
    call gs_fixed_stability(state, law, 1.0_dp, stability, status)
    if (status /= GS_ERROR_UNAVAILABLE) error stop 'moving axis bound invented'

    ! Positive leading derivative is admissible on a finite interior piece;
    ! constant endpoint extensions keep the global upper bound finite.
    deallocate (law%knots, law%pp, law%ff)
    allocate (law%knots(4), law%pp(4, 3), law%ff(4, 3))
    law%normalized = .false.
    law%knots = [-2.0_dp, -1.0_dp, 1.0_dp, 2.0_dp]
    law%pp = 0; law%ff = 0
    law%pp(:, 2) = [0.0_dp, 0.02_dp, -0.03_dp, 0.01_dp]
    call gs_fixed_stability(state, law, 1.0_dp, stability, status)
    if (status /= GS_ERROR_OK) error stop 'finite interior cubic refused'
    if (abs(stability%alpha - (1 - 0.08_dp*cf2)) > 1e-13_dp) &
        error stop 'finite piece endpoint slope maximum'

    call load(rhs)
    call gs_linear_response(state, rhs, psi, status, relative_tolerance=1e-13_dp)
    if (status > 1) error stop 'primal solve'
    call gs_fixed_linear_goal(state, psi, source, weight, goal, status)
    if (status /= GS_ERROR_OK .or. .not. goal%available) error stop 'linear goal'
    if (.not. goal%approximate) error stop 'quadrature status lost'
    ! Independent exact oracle: int (R-2)(4-R)(1-Z^2) dR dZ = 16/9.
    expected = 16.0_dp/9
    if (expected < goal%lower .or. expected > goal%upper) &
        error stop 'manufactured scalar observable outside product band'
    if (abs(goal%correction) > 1e-9_dp) error stop 'Galerkin cancellation'
    call gs_fixed_linear_goal(state, 0.8_dp*psi, source, weight, perturbed, status)
    if (status /= GS_ERROR_OK) error stop 'perturbed candidate'
    if (abs(perturbed%correction) < 0.1_dp) error stop 'residual correction absent'
    if (abs(perturbed%value + perturbed%correction &
            - goal%value - goal%correction) > 1e-9_dp) &
        error stop 'adjoint residual pairing has wrong direction or sign'
    if (expected < perturbed%lower .or. expected > perturbed%upper) &
        error stop 'perturbed observable containment'
    call gs_fixed_linear_goal(state, psi, source, negative_weight, reversed, status)
    if (status /= GS_ERROR_OK) error stop 'negative goal'
    if (abs(reversed%value + goal%value) > 1e-12_dp) error stop 'adjoint direction'
    if (abs(reversed%radius - goal%radius) > 1e-10_dp) &
        error stop 'adjoint signed goal radius'
    print '(a,3es12.4)', 'goal, exact error, product radius: ', goal%value, &
        abs(expected - goal%value), goal%radius
contains
    function source(cell, r, z, psi) result(f)
        integer, intent(in) :: cell
        real(dp), intent(in) :: r, z, psi
        real(dp) :: f, p
        p = (r - 2)*(4 - r)
        f = (2*(1 - z*z) + 2*p)/r + (6 - 2*r)*(1 - z*z)/(r*r) &
            + 0.0_dp*(cell + psi)
    end function
    function weight(r, z) result(q)
        real(dp), intent(in) :: r, z
        real(dp) :: q
        q = 1 + 0.0_dp*(r + z)
    end function
    function negative_weight(r, z) result(q)
        real(dp), intent(in) :: r, z
        real(dp) :: q
        q = -weight(r, z)
    end function
    subroutine load(rhs)
        real(dp), allocatable, intent(out) :: rhs(:)
        real(dp), allocatable :: x(:), y(:), w(:)
        real(dp) :: point(2), det, shape(6), grad(2, 6), hess(2, 2, 6)
        integer :: e, q, stat
        allocate (rhs(size(state%psi)))
        rhs = 0
        call triangle_duffy_quadrature(14, x, y, w, stat)
        if (stat /= 0) error stop 'quadrature'
        do e = 1, size(state%triangles, 2)
            do q = 1, size(w)
                call gs_cell_data(state, e, [x(q), y(q)], point, det, &
                                  shape, grad, hess, stat)
                if (stat /= 0) error stop 'cell'
                rhs(state%triangles(:, e)) = rhs(state%triangles(:, e)) &
                                  + w(q)*det*source(e, point(1), point(2), 0.0_dp)*shape
            end do
        end do
        where (state%boundary) rhs = 0
    end subroutine
end program
