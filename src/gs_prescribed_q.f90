module kin6d_gs_prescribed_q
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64, output_unit
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use kin6d_grad_shafranov, only: gs_state_t, solve_gs, gs_axis_candidate, &
                                    gs_factor_ready
    use kin6d_gs_profile_law, only: gs_profile_law_t, gs_law_primitives
    use kin6d_gs_nonlinear, only: solve_gs_cubic_profiles, gs_forcing_t, &
                                  gs_cubic_load, gs_profile_jacobian_t, &
                                  gs_profile_response, gs_profile_axis
    use kin6d_gs_readback, only: gs_readback_t, read_gs_equilibrium, &
                                 gs_readback_resolution
    use fortnum_bspline, only: bspline_workspace_t, bspline_init, &
                               bspline_set_knots, bspline_eval_basis, bspline_eval_deriv
    use fortnum_multiroot, only: multiroot_hybrids
    use fortnum_linalg, only: dense_solve
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
                                     pressure_law, forcing, knots, jacobian, &
                                     newton_iterations, verbose, probe_x, &
                                     probe_residual, probe_jacobian, probe_seed)
        ! Pressure source dp/dpsi: the scalar pprime on every interval, or, when
        ! pressure_law is present, its piecewise-cubic pp(s) in the normalized
        ! flux s (0 axis, 1 edge; pprime is then ignored). Only the knots and pp
        ! of pressure_law are used. The supported branch has psi_axis < 0,
        ! psi_edge = 0 and p(s) = psi_axis*integral_s^1 pp ds' >= 0 everywhere
        ! (equivalently pprime <= 0 for a constant law). forcing is an optional
        ! extra source for manufactured-solution verification.
        !
        ! FF'(s) is a continuous cubic B-spline. Its break points are the
        ! optional knots (strictly increasing from 0 to 1, any grading) or else
        ! intervals uniform intervals. jacobian selects the outer Newton
        ! Jacobian: 'exact' (default) differentiates the discrete residual, i.e.
        ! one linear solve with the converged weak-equation Jacobian per FF'
        ! coefficient plus a contour readback difference along the resulting
        ! field tangent, and never repeats a nonlinear forward solve per column;
        ! 'secant' is the central-difference fallback of full forward solves.
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
        real(dp), intent(in), optional :: knots(:)
        character(*), intent(in), optional :: jacobian
        integer, intent(out), optional :: newton_iterations
        logical, intent(in), optional :: verbose
        ! Diagnostic evaluation: with probe_x (n+1 outer unknowns: FF' spline
        ! coefficients in units of the seed scale, then log F_edge/seed) the
        ! solver stops after one converged forward solve and returns its outer
        ! residual and, if requested, the exact Jacobian there. The state is
        ! not replaced; status 0 means both were produced.
        real(dp), intent(in), optional :: probe_x(:)
        real(dp), intent(out), optional :: probe_residual(:)
        real(dp), intent(out), optional :: probe_jacobian(:, :)
        real(dp), intent(out), optional :: probe_seed(:)
        type(gs_state_t) :: seed, candidate, accepted_state
        type(gs_profile_law_t) :: candidate_law, plaw, accepted_law
        type(gs_readback_t) :: readback
        type(bspline_workspace_t) :: spline
        type(fortnum_status_t) :: stat
        real(dp), allocatable :: x0(:), solution(:), residual(:), breaks(:)
        real(dp), allocatable :: basis(:, :), conversion(:, :, :), target(:), weights(:)
        real(dp), allocatable :: derivative(:, :), ukn(:), conv(:, :, :)
        real(dp), allocatable :: rho(:), rho_check(:), warm_psi(:)
        integer, allocatable :: ipf(:), ipp(:)
        real(dp) :: axis(2), hessian(2, 2), pressure_axis, unit_axis
        real(dp) :: axis_target, integral, h, x, ffscale, fscale, tol, ffseed
        real(dp) :: left, right, middle, qtol, defect, ppeff, psi_probe
        real(dp) :: pressure, fprobe, jacobian_seconds, solve_seconds
        real(dp) :: t_load, t_response, t_readback
        integer :: n, nb, nr, na, nu, i, j, k, element, calls, failure, newton
        logical :: exact, talk, uniform_grid, good
        status = 1; calls = 0; failure = 0; newton = 0
        jacobian_seconds = 0; solve_seconds = 0
        t_load = 0; t_response = 0; t_readback = 0
        exact = .true.; talk = .false.
        if (present(jacobian)) then
            if (jacobian /= 'exact' .and. jacobian /= 'secant') return
            exact = jacobian == 'exact'
        end if
        if (present(verbose)) talk = verbose
        if (present(newton_iterations)) newton_iterations = 0
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
        if (present(knots)) then
            nb = size(knots)
            if (nb < 2) return
            if (.not. all(ieee_is_finite(knots))) return
            if (abs(knots(1)) > 0 .or. abs(knots(nb) - 1) > 0) return
            if (any(knots(2:) <= knots(:nb - 1))) return
        end if
        if (nb < 2) return
        if (.not. all(ieee_is_finite([tol, qtol]))) return
        if (tol <= 0 .or. tol >= 1) return
        if (qtol <= 0 .or. qtol >= 1) return
        call gs_readback_resolution(state, nr, na)
        if (nr < 5 .or. nr > 2049 .or. na < 16 .or. na > 8192) return
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
        allocate (breaks(nb))
        breaks = [(real(i - 1, dp)/(nb - 1), i=1, nb)]
        if (present(knots)) breaks = knots
        call bspline_set_knots(spline, breaks, stat)
        if (stat%code /= FORTNUM_OK) return
        call sample_grid(nr, breaks, rho, uniform_grid)
        nr = size(rho)
        if (nr > 2049) return
        allocate (conv(4, nb - 1, n), derivative(0:3, n), &
                  basis(n, nr), target(nr), weights(nr), x0(n + 1), &
                  solution(n + 1), residual(n + 1))
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
            x = rho(i)**2
            call bspline_eval_basis(spline, x, basis(:, i), stat)
            if (stat%code /= FORTNUM_OK) return
            target(i) = gs_q_table_value(s, q, x)
            left = rho(max(1, i - 1))**2
            right = rho(min(nr, i + 1))**2
            ! Interior weak q constraints avoid assigning finite measure to
            ! the pointwise axis Hessian or the one-sided boundary gradient.
            weights(i) = (right - left)*x*(1 - x)/2
        end do
        do j = 1, n
            basis(j, :) = basis(j, :)*weights
            ! Every coefficient needs a sample inside its support.
            if (sum(basis(j, :)) <= 0) return
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
        if (present(probe_seed)) then
            if (size(probe_seed) == n + 1) probe_seed = x0
        end if
        if (present(probe_seed) .and. .not. present(probe_x)) then
            status = 0
            return
        end if
        if (present(probe_x)) then
            status = 2
            if (size(probe_x) /= n + 1) return
            if (.not. present(probe_residual)) return
            if (size(probe_residual) /= n + 1) return
            allocate (warm_psi(size(seed%psi)))
            warm_psi = seed%psi
            exact = .false.
            call equations(probe_x, probe_residual)
            if (failure /= 0) then
                status = failure
                return
            end if
            if (present(probe_jacobian)) then
                if (any(shape(probe_jacobian) /= [n + 1, n + 1])) return
                call accept(probe_x, probe_residual)
                call exact_jacobian(probe_jacobian, good)
                if (.not. good) then
                    status = 4
                    return
                end if
            end if
            status = 0
            return
        end if
        if (exact) then
            allocate (warm_psi(size(seed%psi)))
            warm_psi = seed%psi
            call newton_exact(x0)
            if (failure /= 0 .and. newton == 0 .and. any(abs(law%pp(:, 2:)) > 0)) &
                call pressure_homotopy()
        else
            ! The secant fallback proposes coefficients from central differences
            ! of full forward solves; nothing it computes is a sensitivity.
            call multiroot_hybrids(equations, n + 1, x0, solution, stat, &
                                   xtol=tol*0.1_dp, ftol=tol, max_iter=20)
            call equations(solution, residual)
        end if
        if (present(evaluations)) evaluations = calls
        if (present(newton_iterations)) newton_iterations = newton
        status = 2
        if (failure /= 0) then
            status = failure
            return
        end if
        if (maxval(abs(residual)) > tol) return
        state = candidate; law = candidate_law
        ! Independent interleaved radial probes and doubled arc quadrature
        ! check the q profile, beyond the finite set of fitted moments.
        call check_grid(rho, uniform_grid, rho_check)
        call read_gs_equilibrium(state, law%f_edge, 0.0_dp, size(rho_check), &
                                 2*na, readback, status, law, rho_grid=rho_check)
        if (status /= 0) return
        defect = 0
        do i = 1, size(readback%q)
            x = rho_check(i)**2
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
            if (exact) candidate%psi = warm_psi
            candidate_law%ff = 0
            do k = 1, n
                candidate_law%ff = candidate_law%ff &
                                   + ffscale*coefficients(k)*conversion(:, :, k)
            end do
            call join_cubics(candidate_law%knots, candidate_law%ff)
            candidate_law%f_edge = sign(fscale*exp(coefficients(n + 1)), phi)
            call solve_gs_cubic_profiles(candidate, mu0, candidate_law, stat, &
                           min(1.0e-11_dp, tol*0.01_dp), forcing=forcing, &
                                         quadrature_tolerance=1.0e-9_dp)
            if (stat /= 0) then
                failure = 100 + stat
                if (talk) then
                    print '(a,i0,a,es10.3)', '  forward solve status ', stat, &
                        ' residual ', candidate%residual
                    flush (output_unit)
                end if
                return
            end if
            call read_gs_equilibrium(candidate, candidate_law%f_edge, 0.0_dp, &
                              nr, na, readback, stat, candidate_law, rho_grid=rho)
            if (stat /= 0) then
                failure = 200 + stat
                return
            end if
            values(:n) = matmul(basis, readback%q/target - 1)
            values(n + 1) = readback%phi_edge/phi - 1
        end subroutine

        subroutine pressure_homotopy()
            ! The seed field is built with the mean p'. With a strongly
            ! structured p'(s) (pedestal) the first forward field can have no
            ! readable q profile, so the outer Newton cannot start. Solve the
            ! inverse problem for p' blended from the mean to the full law,
            ! each stage warm-started from the previous one, halving the
            ! blend step on failure.
            real(dp), allocatable :: pp_full(:, :), pp_mean(:, :), x_good(:), psi_good(:)
            real(dp) :: blend, blend_ok, step
            integer :: tries
            logical :: stage_ok
            allocate (pp_full(size(law%pp, 1), size(law%pp, 2)))
            allocate (pp_mean, mold=pp_full)
            pp_full = law%pp
            pp_mean = 0; pp_mean(1, :) = ppeff
            allocate (x_good(n + 1), psi_good(size(seed%psi)))
            x_good = x0; psi_good = seed%psi
            blend = 0; blend_ok = -1; step = 0.5_dp
            do tries = 1, 64
                law%pp = (1 - blend)*pp_mean + blend*pp_full
                failure = 0; newton = 0
                warm_psi = psi_good
                call newton_exact(x_good)
                stage_ok = failure == 0
                if (stage_ok) stage_ok = maxval(abs(residual)) <= tol
                if (stage_ok) then
                    x_good = solution; psi_good = accepted_state%psi
                    blend_ok = blend
                    if (talk) print '(a,f8.4)', '  pressure homotopy blend ', blend
                    if (blend >= 1) exit
                    blend = min(1.0_dp, blend + step)
                else
                    if (blend_ok < 0) exit
                    step = step/2
                    if (step < 1.0_dp/32) exit
                    blend = min(1.0_dp, blend_ok + step)
                end if
            end do
            law%pp = pp_full
            if (blend_ok < 1) then
                if (failure == 0) failure = 105
            end if
        end subroutine

        subroutine newton_exact(xstart)
            ! Levenberg-Marquardt-damped Newton with the equation-level Jacobian,
            ! evaluated only at accepted iterates. A trial is the Newton step
            ! (mu = 0) or the minimizer of |r + J dx| + mu |dx| (mu relative to
            ! the largest diagonal of J^T J). Halving the Newton step instead
            ! kept its direction: on TC24 the edge FF' coefficients are weakly
            ! determined by q, the Newton step there is huge, and every scaled
            ! copy of it left the forward problem's solution set. The forward
            ! solve of each trial starts from the last accepted field. Near
            ! the solution mu returns to zero, so convergence is quadratic.
            real(dp), intent(in) :: xstart(:)
            real(dp), allocatable :: x(:), xtrial(:), r(:), rtrial(:), dx(:), &
                                     jac(:, :), normal(:, :), gradient(:), system(:, :)
            real(dp) :: mu, merit, merit_trial, predicted, diagonal
            integer(int64) :: tick, tock, rate
            integer :: iteration, trial, info, k
            logical :: ok, moved
            allocate (x(n + 1), xtrial(n + 1), r(n + 1), rtrial(n + 1), dx(n + 1), &
                      jac(n + 1, n + 1), normal(n + 1, n + 1), gradient(n + 1), &
                      system(n + 1, n + 1))
            x = xstart
            mu = 0
            call equations(x, r)
            solution = x; residual = r
            if (failure /= 0) return
            call accept(x, r)
            do iteration = 1, 25
                if (maxval(abs(r)) <= tol) exit
                call system_clock(tick, rate)
                call exact_jacobian(jac, ok)
                call system_clock(tock)
                jacobian_seconds = real(tock - tick, dp)/rate
                if (.not. ok) then
                    failure = 4
                    exit
                end if
                newton = iteration
                merit = 0.5_dp*dot_product(r, r)
                normal = matmul(transpose(jac), jac)
                gradient = matmul(transpose(jac), r)
                diagonal = maxval([(normal(k, k), k=1, n + 1)])
                moved = .false.
                call system_clock(tick)
                do trial = 1, 16
                    if (mu > 0) then
                        system = normal
                        do k = 1, n + 1
                            system(k, k) = system(k, k) + mu*diagonal
                        end do
                        call dense_solve(system, -gradient, dx, info)
                    else
                        call dense_solve(jac, -r, dx, info)
                    end if
                    if (info /= 0 .or. .not. all(ieee_is_finite(dx))) then
                        mu = max(10*mu, 1.0e-6_dp)
                        cycle
                    end if
                    xtrial = x + dx
                    call equations(xtrial, rtrial)
                    if (failure == 0) then
                        merit_trial = 0.5_dp*dot_product(rtrial, rtrial)
                        predicted = merit - 0.5_dp*sum((r + matmul(jac, dx))**2)
                        if (merit_trial < merit .and. &
                            merit - merit_trial >= 1.0e-4_dp*predicted) then
                            moved = .true.
                            exit
                        end if
                    end if
                    mu = max(10*mu, 1.0e-6_dp)
                end do
                call system_clock(tock)
                solve_seconds = real(tock - tick, dp)/rate
                if (talk) print '(a,i3,a,es10.3,a,f9.2,a,f9.2,a,es9.2)', &
                    'inverse newton ', iteration, ' max|r|=', maxval(abs(r)), &
                    ' jacobian_s=', jacobian_seconds, ' forward_s=', solve_seconds, &
                    ' mu=', mu
                if (talk) flush (output_unit)
                if (.not. moved) then
                    failure = 0
                    exit
                end if
                mu = mu/10
                if (mu < 1.0e-6_dp) mu = 0
                x = xtrial; r = rtrial
                call accept(x, r)
            end do
            ! Report the last accepted iterate.
            if (failure == 0 .or. failure == 4) then
                candidate = accepted_state; candidate_law = accepted_law
            end if
            solution = x; residual = r
            if (failure == 4) return
            failure = 0
        end subroutine

        subroutine perturbed(step, kk, base, dpsi, moved, grad, law_out, readout, st)
            ! Readback of the field moved along dpsi with FF' moved along the
            ! k-th basis function; the normalizing axis flux follows the field.
            real(dp), intent(in) :: step, dpsi(:)
            integer, intent(in) :: kk
            type(gs_profile_law_t), intent(in) :: base
            type(gs_state_t), intent(inout) :: moved
            real(dp), intent(inout) :: grad(:)
            type(gs_profile_law_t), intent(out) :: law_out
            type(gs_readback_t), intent(out) :: readout
            integer, intent(out) :: st
            real(dp) :: axis_value
            law_out = base
            law_out%ff = base%ff + step*ffscale*conversion(:, :, kk)
            call join_cubics(law_out%knots, law_out%ff)
            moved%psi = accepted_state%psi + step*dpsi
            call gs_profile_axis(moved, axis_value, grad, st)
            if (st /= 0) return
            law_out%axis = axis_value
            call read_gs_equilibrium(moved, law_out%f_edge, 0.0_dp, nr, na, &
                                     readout, st, law_out, rho_grid=rho)
        end subroutine

        subroutine accept(x, r)
            real(dp), intent(in) :: x(:), r(:)
            accepted_state = candidate; accepted_law = candidate_law
            warm_psi = candidate%psi
            solution = x; residual = r
        end subroutine

        subroutine exact_jacobian(jac, ok)
            ! Jacobian of the discrete outer residual at the accepted iterate.
            !   column k (FF' coefficient): J dpsi = dL/dc_k with J the
            !     converged weak-equation Jacobian (profile load, axis rank-one
            !     term) and dL/dc_k the load of the k-th basis law at fixed psi;
            !     the residual derivative is the central difference of the
            !     contour readback along (dpsi, dFF'_k) with the axis flux
            !     following the field. No nonlinear forward solve is made.
            !     Loads and readbacks of the columns run on OpenMP threads in
            !     blocks; the response solves between them run one after the
            !     other, because they share the cached MUMPS Cholesky factor
            !     of A (a MUMPS instance is not reentrant) and the phase timers.
            !   column n+1: F_edge enters only the F primitive, not psi.
            real(dp), intent(out) :: jac(:, :)
            logical, intent(out) :: ok
            type(gs_state_t) :: work
            type(gs_profile_law_t) :: base, plus, minus
            type(gs_profile_jacobian_t) :: operator
            type(gs_readback_t) :: rplus, rminus
            real(dp), allocatable :: rhs0(:), rhs_block(:, :), dpsi_block(:, :)
            real(dp) :: load_tol, tl, tr, tb
            integer :: kk, st, first, last, b, nblock
            integer(int64) :: t1, t2, rate
            logical :: column_ok(n), load_ok(n), direct_ready
            integer, parameter :: block_size = 64
            ok = .false.
            jac = 0
            work = accepted_state; base = accepted_law
            allocate (rhs0(size(work%psi)))
            load_tol = 1.0e-6_dp
            call gs_cubic_load(work, mu0, base, rhs0, operator, st, load_tol, &
                               jacobian_tolerance=1.0e-6_dp)
            if (st /= 0) return
            ! The first use creates the factor; the parallel phases only read it.
            direct_ready = gs_factor_ready(accepted_state)
            column_ok = .false.; load_ok = .false.
            tl = 0; tr = 0; tb = 0
            nblock = min(n, block_size)
            allocate (rhs_block(size(work%psi), nblock), dpsi_block(size(work%psi), nblock))
            do first = 1, n, nblock
                last = min(n, first + nblock - 1)
                !$omp parallel do default(shared) schedule(dynamic) &
                !$omp reduction(+:tl) private(kk, b)
                do kk = first, last
                    b = kk - first + 1
                    call column_load(kk, base, rhs_block(:, b), load_ok(kk), tl)
                end do
                !$omp end parallel do
                call system_clock(t1, rate)
                do kk = first, last
                    b = kk - first + 1
                    dpsi_block(:, b) = 0
                    if (.not. load_ok(kk)) cycle
                    call gs_profile_response(accepted_state, operator, &
                                             rhs_block(:, b), dpsi_block(:, b), st)
                    load_ok(kk) = st == 0
                end do
                call system_clock(t2)
                tr = tr + real(t2 - t1, dp)/rate
                !$omp parallel do default(shared) schedule(dynamic) &
                !$omp reduction(+:tb) private(kk, b)
                do kk = first, last
                    b = kk - first + 1
                    if (load_ok(kk)) call column_readback(kk, base, dpsi_block(:, b), &
                                                          jac(:, kk), column_ok(kk), tb)
                end do
                !$omp end parallel do
            end do
            t_load = t_load + tl; t_response = t_response + tr
            t_readback = t_readback + tb
            if (.not. all(column_ok)) return
            ! F_edge direction: the field is untouched.
            plus = base; minus = base
            plus%f_edge = base%f_edge*exp(1.0e-3_dp)
            minus%f_edge = base%f_edge*exp(-1.0e-3_dp)
            allocate (rplus%q(nr), rminus%q(nr))
            call read_gs_equilibrium(accepted_state, plus%f_edge, 0.0_dp, nr, na, &
                                     rplus, st, plus, rho_grid=rho)
            if (st /= 0) return
            call read_gs_equilibrium(accepted_state, minus%f_edge, 0.0_dp, nr, na, &
                                     rminus, st, minus, rho_grid=rho)
            if (st /= 0) return
            jac(:n, n + 1) = matmul(basis, (rplus%q - rminus%q)/(2.0e-3_dp*target))
            jac(n + 1, n + 1) = (rplus%phi_edge - rminus%phi_edge)/(2.0e-3_dp*phi)
            ok = all(ieee_is_finite(jac))
            if (talk) print '(a,3f9.2,a,l2)', '  jacobian load/response/readback thread-s:', &
                t_load, t_response, t_readback, '  direct factor:', direct_ready
            if (talk) flush (output_unit)
        end subroutine

        subroutine column_load(kk, base, rhs_k, good_col, tl)
            ! Load of the k-th basis law at fixed psi (OpenMP-safe: all state
            ! is local, the accepted iterate is read-only).
            integer, intent(in) :: kk
            type(gs_profile_law_t), intent(in) :: base
            real(dp), intent(out) :: rhs_k(:)
            logical, intent(out) :: good_col
            real(dp), intent(inout) :: tl
            type(gs_state_t) :: work
            type(gs_profile_law_t) :: column
            type(gs_profile_jacobian_t) :: unused
            real(dp) :: load_tol
            integer :: jmin, jmax, st, m, nk, i, j
            integer(int64) :: t1, t2, rate
            good_col = .false.
            rhs_k = 0
            work = accepted_state
            load_tol = 1.0e-6_dp
            jmin = 0; jmax = 0
            do j = 1, nu - 1
                if (maxval(abs(conversion(:, j, kk))) <= 0) cycle
                if (jmin == 0) jmin = j
                jmax = j
            end do
            if (jmin == 0) return
            nk = jmax - jmin + 1
            m = nk + 1
            column%normalized = .true.; column%axis = base%axis
            column%f_edge = base%f_edge; column%p_edge = base%p_edge
            if (jmin > 1) m = m + 1
            if (jmax < nu - 1) m = m + 1
            allocate (column%knots(m), column%pp(4, m - 1), column%ff(4, m - 1))
            column%pp = 0; column%ff = 0
            i = 0
            if (jmin > 1) then
                i = i + 1; column%knots(i) = 0
            end if
            column%knots(i + 1:i + nk + 1) = ukn(jmin:jmax + 1)
            column%ff(:, i + 1:i + nk) = conversion(:, jmin:jmax, kk)
            i = i + nk + 1
            if (jmax < nu - 1) column%knots(i + 1) = 1
            call system_clock(t1, rate)
            call gs_cubic_load(work, mu0, column, rhs_k, unused, st, load_tol, &
                               load_only=.true.)
            if (st /= 0) return
            rhs_k = ffscale*rhs_k
            call system_clock(t2)
            tl = tl + real(t2 - t1, dp)/rate
            good_col = .true.
        end subroutine

        subroutine column_readback(kk, base, dpsi, col, good_col, tb)
            ! Central-difference readback of the k-th column along dpsi.
            integer, intent(in) :: kk
            type(gs_profile_law_t), intent(in) :: base
            real(dp), intent(in) :: dpsi(:)
            real(dp), intent(out) :: col(:)
            logical, intent(out) :: good_col
            real(dp), intent(inout) :: tb
            type(gs_state_t) :: moved
            type(gs_profile_law_t) :: plus, minus
            type(gs_readback_t) :: rplus, rminus
            real(dp), allocatable :: grad(:)
            real(dp) :: amplitude, eps
            integer :: attempt, st
            integer(int64) :: t3, t4, rate
            logical :: good
            good_col = .false.
            col = 0
            allocate (rplus%q(nr), rminus%q(nr))
            rplus%q = 0; rminus%q = 0
            moved = accepted_state
            allocate (grad(size(moved%psi)))
            call system_clock(t3, rate)
            ! The step moves the FF' coefficient (unit scale) and the field
            ! together: it is bounded in both. A column whose field response is
            ! small (edge FF' basis) took eps up to 1 and a 1e-2 difference
            ! error on TC24.
            amplitude = maxval(abs(dpsi))
            eps = 1.0e-3_dp
            if (amplitude > 0) eps = min(1.0e-3_dp, 1.0e-3_dp*abs(base%axis)/amplitude)
            good = .false.
            do attempt = 1, 8
                call perturbed(eps, kk, base, dpsi, moved, grad, plus, rplus, st)
                if (st == 0) call perturbed(-eps, kk, base, dpsi, moved, grad, &
                                            minus, rminus, st)
                if (st == 0) then
                    good = .true.
                    exit
                end if
                eps = eps/2
            end do
            if (.not. good) return
            call system_clock(t4)
            tb = tb + real(t4 - t3, dp)/rate
            col(:n) = matmul(basis, (rplus%q - rminus%q)/(2*eps*target))
            col(n + 1) = (rplus%phi_edge - rminus%phi_edge)/(2*eps*phi)
            good_col = .true.
        end subroutine
    end subroutine

    subroutine join_cubics(knots, ff)
        ! The power-basis conversion of a B-spline carries rounding errors of
        ! a few epsilon, which accumulate with the number of knots. The profile
        ! law demands exact continuity, so every interval starts at the value
        ! its predecessor ends with; the change is rounding-sized.
        real(dp), intent(in) :: knots(:)
        real(dp), intent(inout) :: ff(:, :)
        real(dp) :: h
        integer :: i
        do i = 1, size(knots) - 2
            h = knots(i + 1) - knots(i)
            ff(1, i + 1) = ff(1, i) + h*(ff(2, i) + h*(ff(3, i) + h*ff(4, i)))
        end do
    end subroutine

    subroutine sample_grid(nr, breaks, rho, uniform)
        ! Radial sampling in sqrt(s_pol): a uniform grid of nr points plus, in
        ! every FF' spline interval holding fewer than two of them, the points
        ! that bring it to two. A graded knot vector is thereby always sampled
        ! inside every coefficient's support.
        integer, intent(in) :: nr
        real(dp), intent(in) :: breaks(:)
        real(dp), allocatable, intent(out) :: rho(:)
        logical, intent(out) :: uniform
        real(dp), allocatable :: base(:), extra(:)
        real(dp) :: lo, hi, s_point
        integer :: i, m, have, add, t, nx, a, b
        allocate (base(nr), extra(2*size(breaks)))
        base = [(real(i - 1, dp)/(nr - 1), i=1, nr)]
        nx = 0
        do m = 1, size(breaks) - 1
            lo = breaks(m); hi = breaks(m + 1)
            have = count(base**2 > lo .and. base**2 < hi)
            add = max(0, 2 - have)
            do t = 1, add
                s_point = lo + (hi - lo)*t/(add + 1)
                nx = nx + 1
                extra(nx) = sqrt(s_point)
            end do
        end do
        uniform = nx == 0
        allocate (rho(nr + nx))
        a = 1; b = 1
        do i = 1, nr + nx
            if (b > nx) then
                rho(i) = base(a); a = a + 1
            else if (a > nr) then
                rho(i) = extra(b); b = b + 1
            else if (base(a) <= extra(b)) then
                rho(i) = base(a); a = a + 1
            else
                rho(i) = extra(b); b = b + 1
            end if
        end do
    end subroutine

    subroutine check_grid(rho, uniform, finer)
        ! Interleaved probes between the fitted levels. The uniform grid keeps
        ! the exact uniform spacing of twice the points minus one.
        real(dp), intent(in) :: rho(:)
        logical, intent(in) :: uniform
        real(dp), allocatable, intent(out) :: finer(:)
        integer :: i, nr
        nr = size(rho)
        allocate (finer(2*nr - 1))
        if (uniform) then
            finer = [(real(i - 1, dp)/(2*nr - 2), i=1, 2*nr - 1)]
        else
            finer(1:2*nr - 1:2) = rho
            finer(2:2*nr - 2:2) = (rho(:nr - 1) + rho(2:))/2
        end if
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
