module kin6d_slab_oscillation
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortnum_interval, only: interval_t, interval, hull, intersect, is_empty, mag, &
                                linear_interpolation_remainder
    use fortnum_ode_cash_karp, only: cash_karp_step
    use fortnum_status, only: fortnum_status_t, FORTNUM_OK
    use fortnum_ode_residual_certificate, only: certify_linear_trajectory
    use kin6d_slab_oscillation_leaves, only: slab_oscillation_rhs, &
               slab_oscillation_rhs_interval, &
               slab_oscillation_rhs_vjp, slab_oscillation_time, slab_oscillation_span, &
                                             slab_oscillation_hermite, &
                slab_oscillation_hermite_interval, &
           slab_oscillation_initial_radius_interval, &
            slab_oscillation_enlarge_interval, slab_oscillation_error_curvature_interval
    implicit none
    private
    public :: slab_oscillation_state_t, initialize_slab_oscillation
    public :: solve_slab_oscillation, certify_slab_oscillation_trace
    public :: sample_slab_oscillation, sample_slab_oscillation_candidate_interval
    public :: sample_slab_oscillation_interval, slab_oscillation_radius
    public :: slab_oscillation_chord, get_slab_oscillation_trace
    public :: initialize_slab_oscillation_tangent, slab_oscillation_terminal_vjp

    ! Cold rigid electron slab against fixed planar ions, not periodic 3D OCP.
    ! tau=omega_pe*t, q=xi/xi0, p=dq/dtau; q'=p, p'=-q.
    ! The untrusted generic RK candidate is reconstructed componentwise as C1
    ! cubic Hermite q_h (slopes p_k) and p_h (slopes -q_k); FortNum's
    ! fortnum_ode_residual_certificate bounds sup||(q_h'-p_h, p_h'+q_h)||2 by
    ! Bernstein controls and, for the skew generator, accumulates
    ! ||e(t)||2 <= ||e(0)||2 + integral of the residual bound. The bound is
    ! independent of the analytical oracle.
    ! Displacement/density-edge events and fixed integration-mesh leaf derivatives
    ! are separate from an adaptive-schedule derivative or sensitivity certificate.
    type :: slab_oscillation_state_t
        private
        real(dp), allocatable :: time(:), q(:), p(:), radius(:), residual(:)
        real(dp) :: qinitial = 1.0_dp, pinitial = 0.0_dp
        logical :: initialized = .false., prepared = .false.
    end type slab_oscillation_state_t

contains

    subroutine initialize_slab_oscillation(state, steps, duration, &
                                           qinitial, pinitial, ok)
        type(slab_oscillation_state_t), intent(out) :: state
        integer, intent(in) :: steps
        real(dp), intent(in) :: duration, qinitial, pinitial
        logical, intent(out) :: ok
        integer :: j

        ok = .false.
        if (steps < 1 .or. steps >= huge(steps)) return
        if (.not. ieee_is_finite(duration)) return
        if (.not. ieee_is_finite(qinitial) .or. .not. ieee_is_finite(pinitial)) return
        if (duration <= 0.0_dp) return
        allocate (state%time(steps + 1), state%q(steps + 1), state%p(steps + 1), &
                  state%radius(steps + 1), state%residual(steps))
        do j = 1, steps + 1
            call slab_oscillation_time(real(j - 1, dp), real(steps, dp), &
                                       duration, state%time(j))
        end do
        state%time(1) = 0.0_dp
        state%time(steps + 1) = duration
        if (.not. all(ieee_is_finite(state%time))) return
        do j = 1, steps
            if (state%time(j + 1) <= state%time(j)) return
        end do
        state%qinitial = qinitial
        state%pinitial = pinitial
        state%initialized = .true.
        call solve_slab_oscillation(state, ok)
    end subroutine initialize_slab_oscillation

    subroutine solve_slab_oscillation(state, ok, certify)
        type(slab_oscillation_state_t), intent(inout) :: state
        logical, intent(out) :: ok
        logical, intent(in), optional :: certify
        real(dp) :: y(2), ynext(2), k1(2), k2(2), k3(2), k4(2), k5(2), k6(2)
        real(dp) :: scratch(2), embedded_error(2), h
        integer :: j, nfev

        ok = .false.
        state%prepared = .false.
        if (.not. state%initialized) return
        state%q(1) = state%qinitial
        state%p(1) = state%pinitial
        y = [state%qinitial, state%pinitial]
        nfev = 0
        do j = 1, size(state%time) - 1
            call slab_oscillation_span(state%time(j), state%time(j + 1), h)
            call cash_karp_step(rhs, state%time(j), y, h, .false., &
                           k1, k2, k3, k4, k5, k6, scratch, ynext, embedded_error, nfev)
            if (.not. all(ieee_is_finite(ynext))) return
            state%q(j + 1) = ynext(1)
            state%p(j + 1) = ynext(2)
            y = ynext
        end do
        ok = .true.
        if (present(certify)) then
            if (.not. certify) return
        end if
        call certify_slab_oscillation_trace(state%time, state%q, state%p, &
                       state%qinitial, state%pinitial, state%radius, state%residual, ok)
        state%prepared = ok
    end subroutine solve_slab_oscillation

    subroutine rhs(tau, y, derivative, ctx)
        real(dp), intent(in) :: tau, y(:)
        real(dp), intent(out) :: derivative(:)
        class(*), intent(in), optional :: ctx

        call slab_oscillation_rhs(tau, y(1), y(2), derivative(1), derivative(2))
    end subroutine rhs

    subroutine initialize_slab_oscillation_tangent(primal, dq0, dp0, tangent, ok)
        type(slab_oscillation_state_t), intent(in) :: primal
        real(dp), intent(in) :: dq0, dp0
        type(slab_oscillation_state_t), intent(out) :: tangent
        logical, intent(out) :: ok

        ok = .false.
        if (.not. primal%initialized) return
        call initialize_slab_oscillation(tangent, size(primal%time) - 1, &
                                         primal%time(size(primal%time)), dq0, dp0, ok)
    end subroutine initialize_slab_oscillation_tangent

    subroutine slab_oscillation_terminal_vjp(state, bq, bp, gq0, gp0, ok)
        type(slab_oscillation_state_t), intent(in) :: state
        real(dp), intent(in) :: bq, bp
        real(dp), intent(out) :: gq0, gp0
        logical, intent(out) :: ok
        real(dp) :: y(2), ynext(2), k1(2), k2(2), k3(2), k4(2), k5(2), k6(2)
        real(dp) :: scratch(2), embedded_error(2), h
        integer :: j, nfev

        ok = .false.
        gq0 = 0.0_dp; gp0 = 0.0_dp
        if (.not. state%initialized) return
        if (.not. ieee_is_finite(bq) .or. .not. ieee_is_finite(bp)) return
        y = [bq, bp]
        nfev = 0
        ! Continuous adjoint with reverse-time parameter: lambda'=J^T lambda.
        ! For constant J, R(hJ)^T=R(hJ^T) in real RK arithmetic. Floating
        ! rounding is checked separately; no continuous VJP bound is claimed.
        do j = size(state%time) - 1, 1, -1
            call slab_oscillation_span(state%time(j), state%time(j + 1), h)
            call cash_karp_step(adjoint_rhs, state%time(j), y, h, .false., &
                           k1, k2, k3, k4, k5, k6, scratch, ynext, embedded_error, nfev)
            if (.not. all(ieee_is_finite(ynext))) return
            y = ynext
        end do
        gq0 = y(1); gp0 = y(2)
        ok = .true.
    end subroutine slab_oscillation_terminal_vjp

    subroutine adjoint_rhs(tau, y, derivative, ctx)
        real(dp), intent(in) :: tau, y(:)
        real(dp), intent(out) :: derivative(:)
        class(*), intent(in), optional :: ctx
        real(dp) :: gtime

        call slab_oscillation_rhs_vjp(tau, 0.0_dp, 0.0_dp, y(1), y(2), &
                                      gtime, derivative(1), derivative(2))
    end subroutine adjoint_rhs

    subroutine certify_slab_oscillation_trace(time, q, p, qinitial, pinitial, &
                                              radius, residual, ok)
        real(dp), intent(in) :: time(:), q(:), p(:), qinitial, pinitial
        real(dp), intent(out) :: radius(:), residual(:)
        logical, intent(out) :: ok
        real(dp), allocatable :: y(:, :), f(:, :)
        type(interval_t) :: c(2, 0:3), w(2, 0:3), acc
        type(fortnum_status_t) :: status
        integer :: j, n

        ok = .false.
        n = size(time)
        if (n < 2) return
        if (size(q) /= n .or. size(p) /= n .or. size(radius) /= n) return
        if (size(residual) /= n - 1) return
        if (.not. all(ieee_is_finite(q)) .or. .not. all(ieee_is_finite(p))) return
        if (.not. ieee_is_finite(qinitial) .or. .not. ieee_is_finite(pinitial)) return
        if (time(1) < 0.0_dp .or. time(1) > 0.0_dp) return
        call slab_oscillation_initial_radius_interval(interval(q(1)), interval(p(1)), &
                                            interval(qinitial), interval(pinitial), acc)
        if (.not. finite_interval(acc)) return
        allocate (y(2, n), f(2, n))
        do j = 1, n
            y(1, j) = q(j); y(2, j) = p(j)
            call slab_oscillation_rhs(time(j), q(j), p(j), f(1, j), f(2, j))
        end do
        call certify_linear_trajectory(slab_action, time, y, f, interval(0.0_dp), &
                                       acc%hi, c, w, residual, radius, status)
        ok = status%code == FORTNUM_OK
    end subroutine certify_slab_oscillation_trace

    ! Strict interval action of the constant generator [0, 1; -1, 0].
    subroutine slab_action(x, ax, ok, ctx)
        type(interval_t), intent(in) :: x(:)
        type(interval_t), intent(out) :: ax(:)
        logical, intent(out) :: ok
        class(*), intent(in), optional :: ctx

        call slab_oscillation_rhs_interval(interval(0.0_dp), x(1), x(2), ax(1), ax(2))
        ok = finite_interval(ax(1)) .and. finite_interval(ax(2))
    end subroutine slab_action

    subroutine sample_slab_oscillation(state, tau, q, p, qlower, qupper, &
                                       plower, pupper, ok)
        type(slab_oscillation_state_t), intent(in) :: state
        real(dp), intent(in) :: tau
        real(dp), intent(out) :: q, p, qlower, qupper, plower, pupper
        logical, intent(out) :: ok
        type(interval_t) :: qi, pi
        real(dp) :: qpp, pdot
        integer :: j

        ok = .false.
        q = 0.0_dp; p = 0.0_dp
        qlower = 0.0_dp; qupper = 0.0_dp; plower = 0.0_dp; pupper = 0.0_dp
        if (.not. state%prepared) return
        if (.not. ieee_is_finite(tau)) return
        if (tau < state%time(1) .or. tau > state%time(size(state%time))) return
        call sample_slab_oscillation_interval(state, interval(tau), qi, pi, ok)
        if (.not. ok) return
        j = locate(state%time, tau)
        call slab_oscillation_hermite(tau, state%time(j), state%time(j + 1), &
                      state%q(j), state%q(j + 1), state%p(j), state%p(j + 1), q, pdot, qpp)
        ! p_h is its own Hermite with slopes -q_k (exact negation).
        call slab_oscillation_hermite(tau, state%time(j), state%time(j + 1), &
                    state%p(j), state%p(j + 1), -state%q(j), -state%q(j + 1), p, pdot, qpp)
        if (.not. ieee_is_finite(q) .or. .not. ieee_is_finite(p)) then
            ok = .false.
            return
        end if
        qi = hull(qi, interval(q)); pi = hull(pi, interval(p))
        qlower = qi%lo; qupper = qi%hi
        plower = pi%lo; pupper = pi%hi
    end subroutine sample_slab_oscillation

    subroutine sample_slab_oscillation_candidate_interval(state, tau, q, p, ok)
        type(slab_oscillation_state_t), intent(in) :: state
        type(interval_t), intent(in) :: tau
        type(interval_t), intent(out) :: q, p
        logical, intent(out) :: ok

        call sample_packet(state, tau, q, p, .false., ok)
    end subroutine sample_slab_oscillation_candidate_interval

    subroutine sample_slab_oscillation_interval(state, tau, q, p, ok)
        type(slab_oscillation_state_t), intent(in) :: state
        type(interval_t), intent(in) :: tau
        type(interval_t), intent(out) :: q, p
        logical, intent(out) :: ok

        call sample_packet(state, tau, q, p, .true., ok)
    end subroutine sample_slab_oscillation_interval

    subroutine sample_packet(state, tau, q, p, widen, ok)
        type(slab_oscillation_state_t), intent(in) :: state
        type(interval_t), intent(in) :: tau
        type(interval_t), intent(out) :: q, p
        logical, intent(in) :: widen
        logical, intent(out) :: ok
        type(interval_t) :: part, qi, pi, pdot, qpp, lo, hi
        integer :: first, last, j

        ok = .false.
        q = interval(0.0_dp); p = interval(0.0_dp)
        if (.not. state%prepared) return
        if (.not. finite_interval(tau)) return
        if (tau%lo < state%time(1)) return
        if (tau%hi > state%time(size(state%time))) return
        first = locate(state%time, tau%lo)
        last = locate(state%time, tau%hi)
        do j = first, last
            part = intersect(tau, interval(state%time(j), state%time(j + 1)))
            call slab_oscillation_hermite_interval(part, interval(state%time(j)), &
                                    interval(state%time(j + 1)), interval(state%q(j)), &
                                       interval(state%q(j + 1)), interval(state%p(j)), &
                                                  interval(state%p(j + 1)), qi, pdot, qpp)
            call slab_oscillation_hermite_interval(part, interval(state%time(j)), &
                                    interval(state%time(j + 1)), interval(state%p(j)), &
                                      interval(state%p(j + 1)), interval(-state%q(j)), &
                                                interval(-state%q(j + 1)), pi, pdot, qpp)
            if (widen) then
             call slab_oscillation_enlarge_interval(qi, interval(state%radius(j + 1)), &
                                                       lo, hi)
                qi = hull(lo, hi)
             call slab_oscillation_enlarge_interval(pi, interval(state%radius(j + 1)), &
                                                       lo, hi)
                pi = hull(lo, hi)
            end if
            if (.not. finite_interval(qi) .or. .not. finite_interval(pi)) return
            if (j == first) then
                q = qi; p = pi
            else
                q = hull(q, qi); p = hull(p, pi)
            end if
        end do
        ok = .true.
    end subroutine sample_packet

    subroutine slab_oscillation_radius(state, radius, ok)
        type(slab_oscillation_state_t), intent(in) :: state
        real(dp), intent(out) :: radius
        logical, intent(out) :: ok

        ok = state%prepared
        radius = 0.0_dp
        if (.not. ok) return
        radius = state%radius(size(state%radius))
    end subroutine slab_oscillation_radius

    subroutine slab_oscillation_chord(state, tl, tr, padding, ok, error_curve)
        type(slab_oscillation_state_t), intent(in) :: state
        real(dp), intent(in) :: tl, tr
        real(dp), intent(out) :: padding
        logical, intent(out) :: ok
        logical, intent(in), optional :: error_curve
        type(interval_t) :: part, q, p, qpp, curvature, next
        integer :: first, last, j

        ok = .false.
        padding = 0.0_dp
        if (.not. state%prepared) return
        if (.not. ieee_is_finite(tl) .or. .not. ieee_is_finite(tr)) return
        if (tl < state%time(1) .or. tr > state%time(size(state%time))) return
        if (tr < tl) return
        first = locate(state%time, tl)
        last = locate(state%time, tr)
        curvature = interval(0.0_dp)
        do j = first, last
            part = intersect(interval(tl, tr), &
                             interval(state%time(j), state%time(j + 1)))
            call slab_oscillation_hermite_interval(part, interval(state%time(j)), &
                                    interval(state%time(j + 1)), interval(state%q(j)), &
                                       interval(state%q(j + 1)), interval(state%p(j)), &
                                                   interval(state%p(j + 1)), q, p, qpp)
            if (.not. finite_interval(qpp)) return
            curvature = hull(curvature, interval(0.0_dp, mag(qpp)))
        end do
        if (present(error_curve)) then
            if (error_curve) then
                call slab_oscillation_error_curvature_interval(curvature, &
                               interval(state%qinitial), interval(state%pinitial), next)
                if (.not. finite_interval(next)) return
                curvature = interval(0.0_dp, next%hi)
            end if
        end if
        next = linear_interpolation_remainder(interval(tl), interval(tr), curvature)
        if (.not. finite_interval(next)) return
        padding = next%hi
        ok = .true.
    end subroutine slab_oscillation_chord

    subroutine get_slab_oscillation_trace(state, time, q, p, radius, ok)
        type(slab_oscillation_state_t), intent(in) :: state
        real(dp), intent(out) :: time(:), q(:), p(:), radius(:)
        logical, intent(out) :: ok
        integer :: n

        ok = .false.
        if (.not. state%prepared) return
        n = size(state%time)
        if (size(time) /= n .or. size(q) /= n .or. size(p) /= n .or. &
            size(radius) /= n) return
        time = state%time; q = state%q; p = state%p; radius = state%radius
        ok = .true.
    end subroutine get_slab_oscillation_trace

    pure logical function finite_interval(a) result(ok)
        type(interval_t), intent(in) :: a

        ok = .false.
        if (is_empty(a)) return
        if (.not. ieee_is_finite(a%lo) .or. .not. ieee_is_finite(a%hi)) return
        ok = .true.
    end function finite_interval

    integer function locate(grid, value) result(cell)
        real(dp), intent(in) :: grid(:), value
        integer :: left, right, middle

        left = 1; right = size(grid)
        do while (right - left > 1)
            middle = (left + right)/2
            if (grid(middle) <= value) then
                left = middle
            else
                right = middle
            end if
        end do
        cell = min(left, size(grid) - 1)
    end function locate

end module kin6d_slab_oscillation
