module kin6d_gs_prescribed_q
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use kin6d_grad_shafranov, only: gs_state_t, solve_gs, gs_axis_candidate
    use kin6d_gs_profile_law, only: gs_profile_law_t
    use kin6d_gs_nonlinear, only: solve_gs_cubic_profiles
    use kin6d_gs_readback, only: gs_readback_t, read_gs_equilibrium, &
                                 gs_readback_resolution
    use fortnum_bspline, only: bspline_workspace_t, bspline_init, &
                               bspline_set_knots, bspline_eval_basis, bspline_eval_deriv
    use fortnum_multiroot, only: multiroot_hybrids
    use fortnum_polynomial, only: lagrange_weights
    use fortnum_status, only: fortnum_status_t, FORTNUM_OK
    implicit none
    private
    public :: solve_gs_prescribed_q, gs_q_table_value
contains
    real(dp) function gs_q_table_value(s, q, x) result(value)
        real(dp), intent(in) :: s(:), q(:), x
        real(dp) :: weights(4)
        integer :: first
        first = max(1, min(count(s < x) - 1, size(s) - 3))
        call lagrange_weights(4, x, s(first:first + 3), weights)
        value = dot_product(weights, q(first:first + 3))
    end function

    subroutine solve_gs_prescribed_q(state, mu0, pprime, s, q, phi, law, &
                          status, intervals, radial_points, angular_points, tolerance, &
                                     q_error, evaluations, q_tolerance)
        type(gs_state_t), intent(inout) :: state
        real(dp), intent(in) :: mu0, pprime, s(:), q(:), phi
        type(gs_profile_law_t), intent(out) :: law
        integer, intent(out) :: status
        integer, intent(in), optional :: intervals, radial_points, angular_points
        real(dp), intent(in), optional :: tolerance
        real(dp), intent(in), optional :: q_tolerance
        real(dp), intent(out), optional :: q_error
        integer, intent(out), optional :: evaluations
        type(gs_state_t) :: seed, candidate
        type(gs_profile_law_t) :: candidate_law
        type(gs_readback_t) :: readback
        type(bspline_workspace_t) :: spline
        type(fortnum_status_t) :: stat
        real(dp), allocatable :: x0(:), solution(:), residual(:), breaks(:)
        real(dp), allocatable :: basis(:, :), conversion(:, :, :), target(:), weights(:)
        real(dp), allocatable :: derivative(:, :)
        real(dp) :: axis(2), hessian(2, 2), pressure_axis, unit_axis
        real(dp) :: axis_target, integral, h, x, ffscale, fscale, tol, ffseed
        real(dp) :: left, right, middle, qtol, defect
        integer :: n, nb, nr, na, i, j, element, calls, failure
        status = 1; calls = 0; failure = 0
        if (present(q_error)) q_error = huge(1.0_dp)
        if (present(evaluations)) evaluations = 0
        if (size(s) < 4 .or. size(q) /= size(s)) return
        if (.not. all(ieee_is_finite(s)) .or. .not. all(ieee_is_finite(q))) return
        if (.not. all(ieee_is_finite([mu0, pprime, phi]))) return
        if (mu0 <= 0 .or. abs(phi) <= tiny(1.0_dp)) return
        if (any(q*sign(1.0_dp, phi) <= 0)) return
        if (abs(s(1)) > 0 .or. abs(s(size(s)) - 1) > 0) return
        if (any(s(2:) <= s(:size(s) - 1))) return
        if (.not. state%ready .or. state%element_order < 2) return
        if (any(abs(pack(state%psi, state%boundary)) > 0)) return
        nb = 3; nr = 0; na = 0; tol = 1.0e-9_dp; qtol = 1.0e-3_dp
        if (present(intervals)) nb = intervals + 1
        if (present(radial_points)) nr = radial_points
        if (present(angular_points)) na = angular_points
        if (present(tolerance)) tol = tolerance
        if (present(q_tolerance)) qtol = q_tolerance
        if (nb < 2 .or. nb > 17) return
        if (.not. all(ieee_is_finite([tol, qtol]))) return
        if (tol <= 0 .or. tol >= 1) return
        if (qtol <= 0 .or. qtol >= 1) return
        call gs_readback_resolution(state, nr, na)
        if (nr < 2*(nb + 2) .or. nr > 2049 .or. na < 16 .or. na > 8192) return
        integral = 0
        do i = 1, size(s) - 1
            h = s(i + 1) - s(i)
            middle = gs_q_table_value(s, q, (s(i) + s(i + 1))/2)
            integral = integral + h*(q(i) + 4*middle + q(i + 1))/6
        end do
        axis_target = -phi/(2*acos(-1.0_dp)*integral)
        if (axis_target >= 0 .or. pprime > 0) return
        seed = state
        call solve_gs(seed, mu0, 0.0_dp, -1.0_dp, status)
        if (status /= 0) return
        call gs_axis_candidate(seed, axis, unit_axis, hessian, element, status)
        if (status /= 0) return
        pressure_axis = 0
        if (abs(pprime) > 0) then
            call solve_gs(seed, mu0, pprime, 0.0_dp, status)
            if (status /= 0) return
            call gs_axis_candidate(seed, axis, pressure_axis, hessian, element, status)
            if (status /= 0) return
        end if
        ffseed = -(axis_target - pressure_axis)/unit_axis
        ffscale = abs(axis_target/unit_axis)
        call solve_gs(seed, mu0, pprime, ffseed, status)
        if (status /= 0) return
        call read_gs_equilibrium(seed, 1.0_dp, 0.0_dp, nr, na, readback, status)
        if (status /= 0) return
        fscale = abs(phi/readback%phi_edge)
        status = 1
        call bspline_init(spline, 4, nb, stat)
        if (stat%code /= FORTNUM_OK) return
        n = spline%ncoef
        allocate (breaks(nb), conversion(4, nb - 1, n), derivative(0:3, n), &
                  basis(n, nr), target(nr), weights(nr), x0(n + 1), &
                  solution(n + 1), residual(n + 1))
        breaks = [(real(i - 1, dp)/(nb - 1), i=1, nb)]
        call bspline_set_knots(spline, breaks, stat)
        if (stat%code /= FORTNUM_OK) return
        do i = 1, nb - 1
            call bspline_eval_deriv(spline, breaks(i), 3, derivative, stat)
            if (stat%code /= FORTNUM_OK) return
            conversion(:, i, :) = derivative
            conversion(3, i, :) = conversion(3, i, :)/2
            conversion(4, i, :) = conversion(4, i, :)/6
        end do
        do i = 1, nr
            x = (real(i - 1, dp)/(nr - 1))**2
            call bspline_eval_basis(spline, x, basis(:, i), stat)
            if (stat%code /= FORTNUM_OK) return
            target(i) = gs_q_table_value(s, q, x)
            left = (real(max(0, i - 2), dp)/(nr - 1))**2
            right = (real(min(nr - 1, i), dp)/(nr - 1))**2
            ! Interior weak q constraints avoid assigning finite measure to
            ! the pointwise axis Hessian or the one-sided boundary gradient.
            weights(i) = (right - left)*x*(1 - x)/2
        end do
        do j = 1, n
            basis(j, :) = basis(j, :)*weights
            basis(j, :) = basis(j, :)/sum(basis(j, :))
        end do
        law%knots = breaks
        allocate (law%pp(4, nb - 1), law%ff(4, nb - 1))
        law%pp = 0; law%pp(1, :) = pprime
        law%normalized = .true.; law%axis = axis_target
        law%f_edge = sign(fscale, phi)
        x0(:n) = ffseed/ffscale; x0(n + 1) = 0
        ! A numerical root solver proposes coefficients. No iteration or finite-
        ! difference solver Jacobian is exposed as an equilibrium sensitivity.
        call multiroot_hybrids(equations, n + 1, x0, solution, stat, &
                               xtol=tol*0.1_dp, ftol=tol, max_iter=20)
        call equations(solution, residual)
        if (present(evaluations)) evaluations = calls
        status = 2
        if (failure /= 0) then
            status = failure
            return
        end if
        if (maxval(abs(residual)) > tol) return
        state = candidate; law = candidate_law
        ! Independent interleaved radial probes and doubled arc quadrature
        ! check the q profile, beyond the finite set of fitted moments.
        call read_gs_equilibrium(state, law%f_edge, 0.0_dp, 2*nr - 1, &
                                 2*na, readback, status, law)
        if (status /= 0) return
        defect = 0
        do i = 1, size(readback%q)
            x = (real(i - 1, dp)/(size(readback%q) - 1))**2
            defect = max(defect, abs(readback%q(i)/gs_q_table_value(s, q, x) - 1))
        end do
        if (present(q_error)) q_error = defect
        if (defect > qtol) then
            status = 9
            return
        end if
        status = 0
    contains
        subroutine equations(coefficients, values, context)
            real(dp), intent(in) :: coefficients(:)
            real(dp), intent(out) :: values(:)
            class(*), intent(in), optional :: context
            integer :: k, stat
            calls = calls + 1; values = 1.0e3_dp; failure = 0
            if (present(context)) error stop 'unexpected inverse solver context'
            if (.not. all(ieee_is_finite(coefficients))) then
                failure = 3
                return
            end if
            if (abs(coefficients(n + 1)) > 100) then
                failure = 3
                return
            end if
            candidate = seed; candidate_law = law
            candidate_law%ff = 0
            do k = 1, n
                candidate_law%ff = candidate_law%ff &
                                   + ffscale*coefficients(k)*conversion(:, :, k)
            end do
            candidate_law%f_edge = sign(fscale*exp(coefficients(n + 1)), phi)
            call solve_gs_cubic_profiles(candidate, mu0, candidate_law, stat, &
                           min(1.0e-11_dp, tol*0.01_dp), quadrature_tolerance=1.0e-9_dp)
            if (stat /= 0) then
                failure = 100 + stat
                return
            end if
            call read_gs_equilibrium(candidate, candidate_law%f_edge, 0.0_dp, &
                                     nr, na, readback, stat, candidate_law)
            if (stat /= 0) then
                failure = 200 + stat
                return
            end if
            values(:n) = matmul(basis, readback%q/target - 1)
            values(n + 1) = readback%phi_edge/phi - 1
        end subroutine
    end subroutine
end module
