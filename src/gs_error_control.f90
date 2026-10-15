module kin6d_gs_error_control
    !! Continuous fixed-Dirichlet GS error-control composition. All numerical
    !! evaluations here use FP64: analytical coercivity is not an arithmetic
    !! certificate. Physical boundary/model errors remain separate.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use kin6d_grad_shafranov, only: gs_state_t, gs_cell_data, apply_gs, &
                                    gs_linear_response
    use kin6d_triangle_basis, only: triangle_bernstein_controls
    use kin6d_gs_profile_law, only: gs_profile_law_t
    use kin6d_grad_shafranov_leaves, only: gs_cubic_jvp
    use kin6d_gs_majorant, only: gs_point_source_fn, gs_goal_weight_fn
    use kin6d_gs_equilibrated, only: gs_equilibrated_majorant, &
                                     equilibrated_flux_estimate_t, GS_EQ_OK
    use kin6d_gs_error_leaves, only: gs_error_friedrichs_squared, &
                        gs_error_source_lipschitz, gs_error_coercivity, gs_error_goal, &
                                     gs_error_goal_corrected, gs_error_state, &
                                     gs_error_slope_vertex
    use fortfem_triangle_duffy_quadrature, only: triangle_duffy_quadrature
    implicit none
    private
    public :: gs_stability_t, gs_linear_goal_t, gs_fixed_stability
    public :: gs_fixed_source_stability, gs_coercivity_bound
    public :: gs_fixed_linear_goal, gs_state_error_bound
    integer, parameter, public :: GS_ERROR_OK = 0, GS_ERROR_BAD_INPUT = 1
    integer, parameter, public :: GS_ERROR_UNAVAILABLE = 2, GS_ERROR_NONCOERCIVE = 3
    integer, parameter, public :: GS_ERROR_ESTIMATOR_FAILED = 4
    integer, parameter, public :: GS_ERROR_SOLVE_FAILED = 5

    type :: gs_stability_t
        ! source_lipschitz is a nonnegative one-sided derivative upper bound.
        real(dp) :: alpha = 0, source_lipschitz = 0, friedrichs_squared = 0
        logical :: available = .false., arithmetic_rigorous = .false.
    end type

    type :: gs_linear_goal_t
        real(dp) :: value = 0, correction = 0, radius = huge(1.0_dp)
        real(dp) :: lower = -huge(1.0_dp), upper = huge(1.0_dp)
        real(dp) :: primal_majorant = huge(1.0_dp)
        real(dp) :: adjoint_majorant = huge(1.0_dp)
        logical :: approximate = .true., available = .false.
    end type
contains
    subroutine gs_coercivity_bound(cf2, source_lipschitz, alpha, status)
        !! ||v||_L2^2 <= cf2*a(v,v), f'(psi) <= source_lipschitz implies
        !! <F'(psi)v,v> >= alpha*a(v,v), alpha=1-source_lipschitz*cf2.
        real(dp), intent(in) :: cf2, source_lipschitz
        real(dp), intent(out) :: alpha
        integer, intent(out) :: status
        alpha = 0
        status = GS_ERROR_BAD_INPUT
        if (.not. all(ieee_is_finite([cf2, source_lipschitz]))) return
        if (cf2 <= 0 .or. source_lipschitz < 0) return
        call gs_error_coercivity(source_lipschitz, cf2, alpha)
        status = GS_ERROR_NONCOERCIVE
        if (.not. ieee_is_finite(alpha)) return
        if (alpha <= 0) return
        status = GS_ERROR_OK
    end subroutine

    subroutine gs_fixed_source_stability(state, stability, status)
        type(gs_state_t), intent(in) :: state
        type(gs_stability_t), intent(out) :: stability
        integer, intent(out) :: status
        real(dp) :: lo(2), hi(2)
        call domain_bounds(state, lo, hi, status)
        if (status /= GS_ERROR_OK) return
        call gs_error_friedrichs_squared(hi(1), hi(1) - lo(1), &
                                         hi(2) - lo(2), stability%friedrichs_squared)
        call gs_coercivity_bound(stability%friedrichs_squared, 0.0_dp, &
                                 stability%alpha, status)
        stability%available = status == GS_ERROR_OK
    end subroutine

    subroutine gs_fixed_stability(state, law, mu0, stability, status)
        !! Global monotonicity bound for absolute-flux polynomial laws, including
        !! endpoint-polynomial extensions. One-sided derivative maxima cover
        !! every finite interval and both half-lines. Unbounded positive
        !! derivatives and normalized moving-axis laws are refused.
        type(gs_state_t), intent(in) :: state
        type(gs_profile_law_t), intent(in) :: law
        real(dp), intent(in) :: mu0
        type(gs_stability_t), intent(out) :: stability
        integer, intent(out) :: status
        real(dp) :: lo(2), hi(2), pp_upper, ff_upper
        status = GS_ERROR_BAD_INPUT
        if (.not. law%valid()) return
        if (.not. ieee_is_finite(mu0)) return
        if (mu0 <= 0) return
        status = GS_ERROR_UNAVAILABLE
        if (law%normalized) return
        call profile_slope_upper(law%knots, law%pp, pp_upper, status)
        if (status /= GS_ERROR_OK) return
        call profile_slope_upper(law%knots, law%ff, ff_upper, status)
        if (status /= GS_ERROR_OK) return
        call domain_bounds(state, lo, hi, status)
        if (status /= GS_ERROR_OK) return
        call gs_error_friedrichs_squared(hi(1), hi(1) - lo(1), &
                                         hi(2) - lo(2), stability%friedrichs_squared)
        call gs_error_source_lipschitz(mu0, hi(1), lo(1), &
                                       max(0.0_dp, pp_upper), max(0.0_dp, ff_upper), &
                                       stability%source_lipschitz)
        call gs_coercivity_bound(stability%friedrichs_squared, &
                                 stability%source_lipschitz, stability%alpha, status)
        stability%available = status == GS_ERROR_OK
    end subroutine

    subroutine profile_slope_upper(knots, coefficients, upper, status)
        !! Maxima of a quadratic derivative occur at endpoints or its vertex;
        !! the endpoint-polynomial extensions also require tail sign checks.
        real(dp), intent(in) :: knots(:), coefficients(:, :)
        real(dp), intent(out) :: upper
        integer, intent(out) :: status
        real(dp) :: width, vertex, value
        integer :: i, n
        logical :: first, last
        upper = -huge(1.0_dp)
        status = GS_ERROR_UNAVAILABLE
        n = size(coefficients, 2)
        do i = 1, n
            first = i == 1; last = i == n
            if (first .or. last) then
                if (coefficients(4, i) > 0) return
                if (abs(coefficients(4, i)) <= 0) then
                    if (first .and. coefficients(3, i) < 0) return
                    if (last .and. coefficients(3, i) > 0) return
                end if
            end if
            width = knots(i + 1) - knots(i)
            call slope(0.0_dp, i, value)
            if (.not. ieee_is_finite(value)) return
            upper = max(upper, value)
            call slope(width, i, value)
            if (.not. ieee_is_finite(value)) return
            upper = max(upper, value)
            if (coefficients(4, i) < 0) then
                call gs_error_slope_vertex(coefficients(3, i), &
                                           coefficients(4, i), vertex)
                if (.not. ieee_is_finite(vertex)) return
                if (first .or. vertex >= 0) then
                    if (last .or. vertex <= width) then
                        call slope(vertex, i, value)
                        if (.not. ieee_is_finite(value)) return
                        upper = max(upper, value)
                    end if
                end if
            end if
        end do
        if (.not. ieee_is_finite(upper)) return
        status = GS_ERROR_OK
    contains
        subroutine slope(x, interval, value)
            real(dp), intent(in) :: x
            integer, intent(in) :: interval
            real(dp), intent(out) :: value
            call gs_cubic_jvp(x, coefficients(1, interval), &
                              coefficients(2, interval), coefficients(3, interval), &
                              coefficients(4, interval), 1.0_dp, 0.0_dp, 0.0_dp, &
                              0.0_dp, 0.0_dp, value)
        end subroutine
    end subroutine

    subroutine gs_state_error_bound(residual_majorant, stability, radius, status)
        real(dp), intent(in) :: residual_majorant
        type(gs_stability_t), intent(in) :: stability
        real(dp), intent(out) :: radius
        integer, intent(out) :: status
        radius = huge(1.0_dp)
        status = GS_ERROR_UNAVAILABLE
        if (.not. stability%available) return
        status = GS_ERROR_BAD_INPUT
        if (.not. ieee_is_finite(residual_majorant)) return
        if (residual_majorant < 0 .or. stability%alpha <= 0) return
        call gs_error_state(residual_majorant, stability%alpha, radius)
        if (.not. ieee_is_finite(radius)) return
        status = GS_ERROR_OK
    end subroutine

    subroutine domain_bounds(state, lo, hi, status)
        type(gs_state_t), intent(in) :: state
        real(dp), intent(out) :: lo(2), hi(2)
        integer, intent(out) :: status
        real(dp) :: controls(2, 10)
        integer :: e, k, nd
        status = GS_ERROR_BAD_INPUT
        lo = huge(1.0_dp); hi = -huge(1.0_dp)
        if (.not. state%ready) return
        nd = size(state%triangles, 1)
        if (nd /= 3 .and. nd /= 6 .and. nd /= 10) return
        ! Coordinate Bezier hulls enclose curved cells, unlike Lagrange nodes.
        do e = 1, size(state%triangles, 2)
            do k = 1, 2
                call triangle_bernstein_controls( &
                    state%points(k, state%triangles(:, e)), controls(k, :nd))
            end do
            lo = min(lo, minval(controls(:, :nd), dim=2))
            hi = max(hi, maxval(controls(:, :nd), dim=2))
        end do
        if (.not. all(ieee_is_finite([lo, hi]))) return
        if (lo(1) <= 0 .or. any(hi <= lo)) return
        status = GS_ERROR_OK
    end subroutine

    subroutine gs_fixed_linear_goal(state, psi, source, weight, goal, status)
        !! Q(u)=int weight*u for FIXED source -div(grad u/R)=source(x).
        !! The callback MUST be independent of its psi argument; this API does
        !! not cover nonlinear states/goals. Independently solve a(z_h,v)=Q(v),
        !! then |Q(u)-Q(psi)-r(z_h)| <= M_primal*M_adjoint. Same-space Galerkin
        !! cancellation is retained; a discrete residual alone is no bound.
        type(gs_state_t), intent(in) :: state
        real(dp), intent(in) :: psi(:)
        procedure(gs_point_source_fn) :: source
        procedure(gs_goal_weight_fn) :: weight
        type(gs_linear_goal_t), intent(out) :: goal
        integer, intent(out) :: status
        type(equilibrated_flux_estimate_t) :: primal, adjoint
        real(dp), allocatable :: x(:), y(:), w(:), rhs(:, :), z(:), action(:)
        real(dp) :: point(2), det, shape(10), grad(2, 10), hess(2, 2, 10)
        real(dp) :: ph, f, q, center, value, local_rhs(10, 2)
        integer :: stat, e, j, nd, n, a, failed
        status = GS_ERROR_BAD_INPUT
        if (.not. state%ready) return
        n = size(state%psi)
        if (size(psi) /= n) return
        if (.not. all(ieee_is_finite(psi))) return
        if (any(abs(pack(psi, state%boundary)) > tiny(1.0_dp))) return
        nd = size(state%triangles, 1)
        allocate (rhs(n, 2), action(n))
        rhs = 0
        call triangle_duffy_quadrature(14, x, y, w, stat)
        if (stat /= 0) return
        value = 0; failed = 0
        !$omp parallel do default(shared) schedule(static) &
        !$omp private(j, a, point, det, shape, grad, hess, ph, f, q, stat, local_rhs) &
        !$omp reduction(+:value) reduction(max:failed)
        do e = 1, size(state%triangles, 2)
            local_rhs = 0
            do j = 1, size(w)
                call gs_cell_data(state, e, [x(j), y(j)], point, det, &
                                  shape(:nd), grad(:, :nd), hess(:, :, :nd), stat)
                if (stat /= 0) then
                    failed = 1
                    cycle
                end if
                ph = dot_product(shape(:nd), psi(state%triangles(:, e)))
                f = source(e, point(1), point(2), ph)
                q = weight(point(1), point(2))
                if (.not. all(ieee_is_finite([f, q]))) then
                    failed = 1
                    cycle
                end if
                local_rhs(:nd, 1) = local_rhs(:nd, 1) + w(j)*det*f*shape(:nd)
                local_rhs(:nd, 2) = local_rhs(:nd, 2) + w(j)*det*q*shape(:nd)
                value = value + w(j)*det*q*ph
            end do
            do a = 1, nd
                !$omp atomic update
                rhs(state%triangles(a, e), 1) = rhs(state%triangles(a, e), 1) &
                    + local_rhs(a, 1)
                !$omp atomic update
                rhs(state%triangles(a, e), 2) = rhs(state%triangles(a, e), 2) &
                    + local_rhs(a, 2)
            end do
        end do
        !$omp end parallel do
        if (failed /= 0) return
        goal%value = value
        where (spread(state%boundary, 2, 2)) rhs = 0
        status = GS_ERROR_SOLVE_FAILED
        call gs_linear_response(state, rhs(:, 2), z, stat, &
                                relative_tolerance=1.0e-12_dp)
        if (stat > 1) return
        call apply_gs(state, psi, action)
        goal%correction = dot_product(z, rhs(:, 1) - action)
        status = GS_ERROR_ESTIMATOR_FAILED
        call gs_equilibrated_majorant(state, psi, source, primal, stat)
        if (stat /= GS_EQ_OK) return
        call gs_equilibrated_majorant(state, z, adjoint_source, adjoint, stat)
        if (stat /= GS_EQ_OK) return
        goal%primal_majorant = primal%bound
        goal%adjoint_majorant = adjoint%bound
        call gs_error_goal(primal%bound, adjoint%bound, 1.0_dp, goal%radius)
        call gs_error_goal_corrected(goal%value, goal%correction, center)
        goal%lower = center - goal%radius
        goal%upper = center + goal%radius
        if (.not. all(ieee_is_finite([center, goal%radius]))) return
        goal%available = .true.
        status = GS_ERROR_OK
    contains
        function adjoint_source(cell, r, zz, p) result(f)
            integer, intent(in) :: cell
            real(dp), intent(in) :: r, zz, p
            real(dp) :: f
            f = weight(r, zz) + 0.0_dp*(cell + p)
        end function
    end subroutine
end module kin6d_gs_error_control
