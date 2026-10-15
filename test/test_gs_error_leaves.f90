program test_gs_error_leaves
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_interval, only: interval_t, interval, contains
    use kin6d_grad_shafranov_leaves, only: gs_source_interval
    use kin6d_gs_error_leaves, only: &
        gs_error_coefficient, gs_error_coefficient_jvp, &
        gs_error_coefficient_vjp, gs_error_flux, &
        gs_error_flux_jvp, gs_error_flux_vjp, &
        gs_error_friedrichs, gs_error_friedrichs_jvp, &
        gs_error_friedrichs_vjp, gs_error_friedrichs_squared, &
        gs_error_friedrichs_squared_jvp, gs_error_friedrichs_squared_vjp, &
        gs_error_coercivity, gs_error_coercivity_jvp, &
        gs_error_coercivity_vjp, gs_error_source_lipschitz, &
        gs_error_source_lipschitz_jvp, gs_error_source_lipschitz_vjp, &
        gs_error_state, gs_error_state_jvp, &
        gs_error_state_vjp, gs_error_goal, &
        gs_error_goal_jvp, gs_error_goal_vjp, &
        gs_error_goal_corrected, gs_error_goal_corrected_jvp, &
        gs_error_goal_corrected_vjp
    implicit none
    real(dp) :: x(5), direction(5), weights(2), value(2), tangent(2), gradient(5)
    real(dp) :: plus(2), minus(2), oracle(2), h, cf, cf2, pi
    type(interval_t) :: source_enclosure
    integer :: leaf, trial, ir, ip, iff
    real(dp) :: radius, pprime, ffprime, exact_source

    pi = acos(-1.0_dp)
    h = 2.0e-5_dp
    direction = [0.13_dp, -0.21_dp, 0.17_dp, 0.09_dp, -0.07_dp]
    weights = [-0.7_dp, 1.3_dp]
    do leaf = 1, 9
        do trial = 1, 8
            x = [2.0_dp, 1.1_dp, 0.8_dp, 0.3_dp, 0.4_dp] &
                + 0.03_dp*real(trial, dp)
            call evaluate(leaf, x, value)
            call derivative(leaf, x, direction, weights, tangent, gradient)
            call evaluate(leaf, x + h*direction, plus)
            call evaluate(leaf, x - h*direction, minus)
            if (maxval(abs(tangent - (plus - minus)/(2*h))) > 2.0e-9_dp) &
                error stop 'GS error leaf independent centered difference'
            if (abs(dot_product(weights, tangent) &
                - dot_product(direction, gradient)) > 1.0e-13_dp) &
                error stop 'GS error leaf JVP/VJP Euclidean pairing'
            oracle = 0
            select case (leaf)
            case (1)
                oracle(1) = 1/x(1)
            case (2)
                oracle = x(2:3)/x(1)
            case (3)
                oracle(1) = sqrt(x(1)/(pi*pi*(x(2)**2 + x(3)**2))) &
                    *x(2)*x(3)
            case (4)
                oracle(1) = x(1)*x(2)**2*x(3)**2 &
                    /(pi*pi*(x(2)**2 + x(3)**2))
            case (5)
                oracle(1) = 1 - x(1)*x(2)
            case (6)
                oracle(1) = x(1)*x(2)*x(4) + x(5)/x(3)
            case (7)
                oracle(1) = x(1)/x(2)
            case (8)
                oracle(1) = x(1)*x(2)/x(3)
            case (9)
                oracle(1) = x(1) + x(2)
            end select
            if (maxval(abs(value - oracle)) > 5.0e-14_dp) &
                error stop 'GS error leaf independent primal formula'
        end do
    end do
    ! Domain monotonicity: the squared rectangle constant equals the inverse
    ! lowest Dirichlet Laplacian eigenvalue times the maximum major radius.
    call gs_error_friedrichs(3.0_dp, 2.0_dp, 1.0_dp, cf)
    call gs_error_friedrichs_squared(3.0_dp, 2.0_dp, 1.0_dp, cf2)
    if (abs(cf**2 - cf2) > 1.0e-15_dp) error stop 'Friedrichs squared consistency'
    if (abs(cf2*(pi**2/4 + pi**2) - 3) > 1.0e-14_dp) &
        error stop 'rectangle Dirichlet eigenvalue identity'
    ! Enclose a whole input box; exact arithmetic source values are sampled
    ! independently, with both source signs included.
    call gs_source_interval(interval(2.0_dp, 3.0_dp), interval(0.25_dp), &
        interval(-2.0_dp, 1.0_dp), interval(-3.0_dp, 2.0_dp), source_enclosure)
    do ir = 0, 4
        radius = 2.0_dp + real(ir, dp)/4
        do ip = 0, 4
            pprime = -2.0_dp + 3*real(ip, dp)/4
            do iff = 0, 4
                ffprime = -3.0_dp + 5*real(iff, dp)/4
                exact_source = 0.25_dp*radius*pprime + ffprime/radius
                if (.not. contains(source_enclosure, exact_source)) &
                    error stop 'GS interval source box containment'
            end do
        end do
    end do
    print *, 'PASS: GS error leaf formulas, directional differences and adjoint pairing'
contains
    subroutine evaluate(leaf, x, y)
        integer, intent(in) :: leaf
        real(dp), intent(in) :: x(5)
        real(dp), intent(out) :: y(2)
        y = 0
        select case (leaf)
        case (1)
            call gs_error_coefficient(x(1), y(1))
        case (2)
            call gs_error_flux(x(1), x(2), x(3), y(1), y(2))
        case (3)
            call gs_error_friedrichs(x(1), x(2), x(3), y(1))
        case (4)
            call gs_error_friedrichs_squared(x(1), x(2), x(3), y(1))
        case (5)
            call gs_error_coercivity(x(1), x(2), y(1))
        case (6)
            call gs_error_source_lipschitz(x(1), x(2), x(3), x(4), x(5), y(1))
        case (7)
            call gs_error_state(x(1), x(2), y(1))
        case (8)
            call gs_error_goal(x(1), x(2), x(3), y(1))
        case (9)
            call gs_error_goal_corrected(x(1), x(2), y(1))
        end select
    end subroutine evaluate

    subroutine derivative(leaf, x, d, b, tangent, gradient)
        integer, intent(in) :: leaf
        real(dp), intent(in) :: x(5), d(5), b(2)
        real(dp), intent(out) :: tangent(2), gradient(5)
        tangent = 0; gradient = 0
        select case (leaf)
        case (1)
            call gs_error_coefficient_jvp(x(1), d(1), tangent(1))
            call gs_error_coefficient_vjp(x(1), b(1), gradient(1))
        case (2)
            call gs_error_flux_jvp(x(1), x(2), x(3), d(1), d(2), d(3), &
                tangent(1), tangent(2))
            call gs_error_flux_vjp(x(1), x(2), x(3), b(1), b(2), &
                gradient(1), gradient(2), gradient(3))
        case (3)
            call gs_error_friedrichs_jvp(x(1), x(2), x(3), d(1), d(2), d(3), tangent(1))
            call gs_error_friedrichs_vjp(x(1), x(2), x(3), b(1), &
                gradient(1), gradient(2), gradient(3))
        case (4)
            call gs_error_friedrichs_squared_jvp(x(1), x(2), x(3), &
                d(1), d(2), d(3), tangent(1))
            call gs_error_friedrichs_squared_vjp(x(1), x(2), x(3), b(1), &
                gradient(1), gradient(2), gradient(3))
        case (5)
            call gs_error_coercivity_jvp(x(1), x(2), d(1), d(2), tangent(1))
            call gs_error_coercivity_vjp(x(1), x(2), b(1), gradient(1), gradient(2))
        case (6)
            call gs_error_source_lipschitz_jvp(x(1), x(2), x(3), x(4), x(5), &
                d(1), d(2), d(3), d(4), d(5), tangent(1))
            call gs_error_source_lipschitz_vjp(x(1), x(2), x(3), x(4), x(5), b(1), &
                gradient(1), gradient(2), gradient(3), gradient(4), gradient(5))
        case (7)
            call gs_error_state_jvp(x(1), x(2), d(1), d(2), tangent(1))
            call gs_error_state_vjp(x(1), x(2), b(1), gradient(1), gradient(2))
        case (8)
            call gs_error_goal_jvp(x(1), x(2), x(3), d(1), d(2), d(3), tangent(1))
            call gs_error_goal_vjp(x(1), x(2), x(3), b(1), &
                gradient(1), gradient(2), gradient(3))
        case (9)
            call gs_error_goal_corrected_jvp(x(1), x(2), d(1), d(2), tangent(1))
            call gs_error_goal_corrected_vjp(x(1), x(2), b(1), gradient(1), gradient(2))
        end select
    end subroutine derivative
end program test_gs_error_leaves
