module kin6d_gs_prescribed_q
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use kin6d_grad_shafranov, only: gs_state_t, solve_gs, gs_axis_candidate
    use kin6d_gs_profile_law, only: gs_profile_law_t, gs_law_primitives
    use kin6d_gs_nonlinear, only: solve_gs_cubic_profiles, gs_forcing_t
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
                                     q_error, evaluations, q_tolerance, &
                                     pressure_law, forcing)
        ! Pressure source dp/dpsi: the scalar pprime on every interval, or, when
        ! pressure_law is present, its piecewise-cubic pp(s) in the normalized
        ! flux s (0 axis, 1 edge; pprime is then ignored). Only the knots and pp
        ! of pressure_law are used. The supported branch has psi_axis < 0,
        ! psi_edge = 0 and p(s) = psi_axis*integral_s^1 pp ds' >= 0 everywhere
        ! (equivalently pprime <= 0 for a constant law). forcing is an optional
        ! extra source for manufactured-solution verification.
        type(gs_state_t), intent(inout) :: state
        real(dp), intent(in) :: mu0, pprime, s(:), q(:), phi
        type(gs_profile_law_t), intent(in), optional :: pressure_law
        procedure(gs_forcing_t), optional :: forcing
        type(gs_profile_law_t), intent(out) :: law
        integer, intent(out) :: status
        integer, intent(in), optional :: intervals, radial_points, angular_points
        real(dp), intent(in), optional :: tolerance
        real(dp), intent(in), optional :: q_tolerance
        real(dp), intent(out), optional :: q_error
        integer, intent(out), optional :: evaluations
        type(gs_state_t) :: seed, candidate
        type(gs_profile_law_t) :: candidate_law, plaw
        type(gs_readback_t) :: readback
        type(bspline_workspace_t) :: spline
        type(fortnum_status_t) :: stat
        real(dp), allocatable :: x0(:), solution(:), residual(:), breaks(:)
        real(dp), allocatable :: basis(:, :), conversion(:, :, :), target(:), weights(:)
        real(dp), allocatable :: derivative(:, :), ukn(:), conv(:, :, :)
        integer, allocatable :: ipf(:), ipp(:)
        real(dp) :: axis(2), hessian(2, 2), pressure_axis, unit_axis
        real(dp) :: axis_target, integral, h, x, ffscale, fscale, tol, ffseed
        real(dp) :: left, right, middle, qtol, defect, ppeff, psi_probe
        real(dp) :: pressure, fprobe
        integer :: n, nb, nr, na, nu, i, j, k, element, calls, failure
        status = 1; calls = 0; failure = 0
        if (present(q_error)) q_error = huge(1.0_dp)
        if (present(evaluations)) evaluations = 0
        if (size(s) < 4 .or. size(q) /= size(s)) return
        if (.not. all(ieee_is_finite(s)) .or. .not. all(ieee_is_finite(q))) return
        if (.not. all(ieee_is_finite([mu0, phi]))) return
        if (present(pressure_law)) then
            if (.not. allocated(pressure_law%knots)) return
            if (.not. allocated(pressure_law%pp)) return
            if (.not. pressure_law%normalized) return
            if (size(pressure_law%knots) < 2) return
            if (any(shape(pressure_law%pp) /= [4, size(pressure_law%knots) - 1])) return
            plaw%knots = pressure_law%knots; plaw%pp = pressure_law%pp
            plaw%axis = pressure_law%axis
        else
            if (.not. ieee_is_finite(pprime)) return
            plaw%knots = [0.0_dp, 1.0_dp]
            allocate (plaw%pp(4, 1)); plaw%pp = 0; plaw%pp(1, 1) = pprime
        end if
        allocate (plaw%ff(4, size(plaw%knots) - 1)); plaw%ff = 0
        plaw%normalized = .true.; plaw%f_edge = 1; plaw%p_edge = 0
        if (abs(plaw%axis) <= tiny(1.0_dp) .or. plaw%axis > 0) plaw%axis = -1
        if (.not. plaw%valid()) return
        if (abs(plaw%knots(1)) > 0) return
        if (abs(plaw%knots(size(plaw%knots)) - 1) > 0) return
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
        if (axis_target >= 0) return
        ! Per-point pressure sign: p(s) >= 0 on the supported branch, checked
        ! with the law's own primitive on a dense probe of every interval.
        plaw%axis = axis_target
        do j = 1, size(plaw%knots) - 1
            do k = 0, 32
                left = plaw%knots(j) + (plaw%knots(j + 1) - plaw%knots(j))*k/32
                psi_probe = axis_target*(1 - left)
                call gs_law_primitives(plaw, psi_probe, pressure, fprobe, status)
                if (status /= 0) then
                    status = 1
                    return
                end if
            end do
        end do
        ppeff = 0
        do j = 1, size(plaw%knots) - 1
            h = plaw%knots(j + 1) - plaw%knots(j)
            ppeff = ppeff + h*(plaw%pp(1, j) + h*(plaw%pp(2, j)/2 &
                                + h*(plaw%pp(3, j)/3 + h*plaw%pp(4, j)/4)))
        end do
        seed = state
        call solve_gs(seed, mu0, 0.0_dp, -1.0_dp, status)
        if (status /= 0) return
        call gs_axis_candidate(seed, axis, unit_axis, hessian, element, status)
        if (status /= 0) return
        pressure_axis = 0
        if (abs(ppeff) > 0) then
            call solve_gs(seed, mu0, ppeff, 0.0_dp, status)
            if (status /= 0) return
            call gs_axis_candidate(seed, axis, pressure_axis, hessian, element, status)
            if (status /= 0) return
        end if
        ffseed = -(axis_target - pressure_axis)/unit_axis
        ffscale = abs(axis_target/unit_axis)
        call solve_gs(seed, mu0, ppeff, ffseed, status)
        if (status /= 0) return
        call read_gs_equilibrium(seed, 1.0_dp, 0.0_dp, nr, na, readback, status)
        if (status /= 0) return
        fscale = abs(phi/readback%phi_edge)
        status = 1
        call bspline_init(spline, 4, nb, stat)
        if (stat%code /= FORTNUM_OK) return
        n = spline%ncoef
        allocate (breaks(nb), conv(4, nb - 1, n), derivative(0:3, n), &
                  basis(n, nr), target(nr), weights(nr), x0(n + 1), &
                  solution(n + 1), residual(n + 1))
        breaks = [(real(i - 1, dp)/(nb - 1), i=1, nb)]
        call bspline_set_knots(spline, breaks, stat)
        if (stat%code /= FORTNUM_OK) return
        do i = 1, nb - 1
            call bspline_eval_deriv(spline, breaks(i), 3, derivative, stat)
            if (stat%code /= FORTNUM_OK) return
            conv(:, i, :) = derivative
            conv(3, i, :) = conv(3, i, :)/2
            conv(4, i, :) = conv(4, i, :)/6
        end do
        call union_knots(breaks, plaw%knots, ukn, ipf, ipp)
        nu = size(ukn)
        allocate (conversion(4, nu - 1, n))
        do j = 1, nu - 1
            do k = 1, n
                conversion(:, j, k) = shifted(conv(:, ipf(j), k), &
                                              ukn(j) - breaks(ipf(j)))
            end do
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
        law%knots = ukn
        allocate (law%pp(4, nu - 1), law%ff(4, nu - 1))
        law%ff = 0
        do j = 1, nu - 1
            law%pp(:, j) = shifted(plaw%pp(:, ipp(j)), ukn(j) - plaw%knots(ipp(j)))
        end do
        law%normalized = .true.; law%axis = axis_target
        law%f_edge = sign(fscale, phi)
        if (.not. law%valid()) return
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
                           min(1.0e-11_dp, tol*0.01_dp), forcing=forcing, &
                                         quadrature_tolerance=1.0e-9_dp)
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

    pure function shifted(c, d) result(r)
        ! Cubic coefficients about a point displaced by d from the expansion point.
        real(dp), intent(in) :: c(4), d
        real(dp) :: r(4)
        r(1) = c(1) + d*(c(2) + d*(c(3) + d*c(4)))
        r(2) = c(2) + d*(2*c(3) + 3*d*c(4))
        r(3) = c(3) + 3*d*c(4)
        r(4) = c(4)
    end function

    subroutine union_knots(first, second, knots, parent_first, parent_second)
        ! Sorted union of two knot sets (merged within 1e-12) and the interval
        ! of each set containing every union interval.
        real(dp), intent(in) :: first(:), second(:)
        real(dp), allocatable, intent(out) :: knots(:)
        integer, allocatable, intent(out) :: parent_first(:), parent_second(:)
        real(dp), allocatable :: merged(:)
        real(dp) :: middle
        integer :: i, j, m
        allocate (merged(size(first) + size(second)))
        merged(:size(first)) = first
        m = size(first)
        do i = 1, size(second)
            if (all(abs(first - second(i)) > 1.0e-12_dp)) then
                m = m + 1; merged(m) = second(i)
            end if
        end do
        knots = merged(:m)
        do i = 2, m
            middle = knots(i)
            j = i - 1
            do while (j >= 1)
                if (knots(j) <= middle) exit
                knots(j + 1) = knots(j); j = j - 1
            end do
            knots(j + 1) = middle
        end do
        allocate (parent_first(m - 1), parent_second(m - 1))
        do i = 1, m - 1
            middle = (knots(i) + knots(i + 1))/2
            parent_first(i) = min(size(first) - 1, max(1, count(first <= middle)))
            parent_second(i) = min(size(second) - 1, max(1, count(second <= middle)))
        end do
    end subroutine
end module
