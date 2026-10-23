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
    interface
        subroutine dgels(trans, m, n, nrhs, a, lda, b, ldb, work, lwork, info)
            import :: dp
            character, intent(in) :: trans
            integer, intent(in) :: m, n, nrhs, lda, ldb, lwork
            real(dp), intent(inout) :: a(lda, *), b(ldb, *)
            real(dp), intent(out) :: work(*)
            integer, intent(out) :: info
        end subroutine
    end interface
    type :: solve_gs_prescribed_q_equations_context_t
        real(dp), pointer :: basis(:,:) => null()
        integer, pointer :: calls => null()
        type(gs_state_t), pointer :: candidate => null()
        type(gs_profile_law_t), pointer :: candidate_law => null()
        real(dp), pointer :: conversion(:,:,:) => null()
        logical, pointer :: exact => null()
        integer, pointer :: failure => null()
        real(dp), pointer :: ffscale => null()
        procedure(gs_forcing_t), pointer, nopass :: forcing => null()
        real(dp), pointer :: fscale => null()
        type(gs_profile_law_t), pointer :: law => null()
        real(dp), pointer :: mu0 => null()
        integer, pointer :: n => null()
        integer, pointer :: na => null()
        integer, pointer :: nr => null()
        real(dp), pointer :: phi => null()
        type(gs_readback_t), pointer :: readback => null()
        real(dp), pointer :: rho(:) => null()
        real(dp), pointer :: sample_now(:) => null()
        type(gs_state_t), pointer :: seed => null()
        logical, pointer :: talk => null()
        real(dp), pointer :: target(:) => null()
        real(dp), pointer :: tol => null()
        real(dp), pointer :: warm_psi(:) => null()
        real(dp), pointer :: wsqrt(:) => null()
    end type
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
                                     probe_residual, probe_jacobian, probe_seed, &
                                     regularization, regularization_weight, &
                                     regularization_used, q_accuracy)
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
        real(dp), intent(in), target :: mu0, pprime, s(:), q(:), phi
        type(gs_profile_law_t), intent(in), optional :: pressure_law
        procedure(gs_forcing_t), optional :: forcing
        type(gs_profile_law_t), intent(out), target :: law
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
        ! Smoothness regularization of FF' (default 'none': the n q moments
        ! and the flux equation are solved exactly). The exact moments ask the
        ! discrete forward model to reproduce q to better than its own
        ! discretization error; FF' modes that q barely sees (an alternating
        ! edge mode on edge-graded knots, Jacobian cond 1.7e5 on TC24) then
        ! absorb that error and the iteration ends at a fold. Regularized,
        ! the weighted q samples themselves are fitted (overdetermined least
        ! squares, w = sample weights normalized to unit sum) together with
        ! the flux equation, and
        !   |r|^2 + lambda^2 integral_0^1 (d^2 FF'/ds^2)^2 ds
        ! is minimized (FF' in units of the seed scale; F_edge free) by
        ! Levenberg-Marquardt-damped Gauss-Newton. 'discrepancy' (Morozov)
        ! chooses lambda so that the weighted rms q misfit equals q_accuracy,
        ! the stated accuracy of the q data relative to the discrete model.
        ! It is re-chosen on each linearization (goal max(q_accuracy, half the
        ! current misfit), a continuation) and frozen once the goal is
        ! q_accuracy and two choices agree within 25 %; the result is the
        ! Tikhonov minimizer at that lambda. A positive regularization_weight
        ! fixes lambda instead. Convergence: the Gauss-Newton step changes no
        ! residual or penalty component by more than tolerance at the final
        ! lambda (regularization_used).
        character(*), intent(in), optional :: regularization
        real(dp), intent(in), optional :: regularization_weight, q_accuracy
        real(dp), intent(out), optional :: regularization_used
        type(gs_state_t), target :: seed, candidate, accepted_state
        type(gs_profile_law_t), target :: candidate_law, plaw, accepted_law
        type(gs_readback_t), target :: readback
        type(bspline_workspace_t) :: spline
        type(fortnum_status_t) :: stat
        real(dp), allocatable :: x0(:), solution(:), residual(:), breaks(:)
        real(dp), allocatable , target :: basis(:, :), conversion(:, :, :), &
            target(:), weights(:)
        real(dp), allocatable :: derivative(:, :), ukn(:), conv(:, :, :)
        real(dp), allocatable , target :: rho(:), rho_check(:), warm_psi(:)
        integer, allocatable :: ipf(:), ipp(:)
        real(dp) :: axis(2), hessian(2, 2), pressure_axis, unit_axis
        real(dp), target :: axis_target, integral, h, x, ffscale, fscale, tol, ffseed
        real(dp) :: left, right, middle, qtol, defect, ppeff, psi_probe
        real(dp) :: pressure, fprobe, jacobian_seconds, solve_seconds
        real(dp) :: t_load, t_response, t_readback
        integer , target :: n, nb, nr, na, nu, i, j, k, element, calls, failure, newton
        logical , target :: exact, talk, uniform_grid, good, regularized, adaptive, &
            converged
        real(dp), allocatable , target :: smooth_rows(:, :), wsqrt(:), sample_now(:)
        real(dp), allocatable :: jac_samples(:, :)
        real(dp) :: lambda_fixed, lambda, misfit_target
        type(solve_gs_prescribed_q_equations_context_t) :: equations_context
        status = 1; calls = 0; failure = 0; newton = 0
        jacobian_seconds = 0; solve_seconds = 0
        t_load = 0; t_response = 0; t_readback = 0
        exact = .true.; talk = .false.
        if (present(jacobian)) then
            if (jacobian /= 'exact' .and. jacobian /= 'secant') return
            exact = jacobian == 'exact'
        end if
        if (present(verbose)) talk = verbose
        regularized = .false.; adaptive = .false.; lambda_fixed = 0; lambda = 0
        misfit_target = 0
        converged = .false.
        if (present(regularization_used)) regularization_used = 0
        if (present(regularization)) then
            if (regularization /= 'none' .and. regularization /= 'discrepancy') return
            adaptive = regularization == 'discrepancy'
        end if
        if (present(regularization_weight)) then
            if (.not. ieee_is_finite(regularization_weight)) return
            if (regularization_weight < 0) return
            lambda_fixed = regularization_weight
        end if
        if (adaptive) then
            if (lambda_fixed > 0) return
            if (.not. present(q_accuracy)) return
            if (.not. ieee_is_finite(q_accuracy)) return
            if (q_accuracy <= 0 .or. q_accuracy >= 1) return
            misfit_target = q_accuracy
        end if
        regularized = adaptive .or. lambda_fixed > 0
        if (regularized .and. .not. exact) return
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
        ! Rows of the FF'' seminorm: two Gauss points per break interval
        ! integrate the squared (piecewise linear) second derivative exactly.
        allocate (smooth_rows(2*(nb - 1), n + 1))
        smooth_rows = 0
        do i = 1, nb - 1
            h = breaks(i + 1) - breaks(i)
            do k = 1, 2
                x = (breaks(i) + breaks(i + 1))/2 + (2*k - 3)*h/(2*sqrt(3.0_dp))
                call bspline_eval_deriv(spline, x, 3, derivative, stat)
                if (stat%code /= FORTNUM_OK) return
                smooth_rows(2*(i - 1) + k, :n) = sqrt(h/2)*derivative(2, :)
            end do
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
        ! Regularized mode fits the samples themselves (an overdetermined
        ! least-squares problem with a data/model floor) with weights uniform
        ! in rho = sqrt(s), reaching the axis and the edge: the q check bounds
        ! q there too, and the regularization, not a vanishing weight, keeps
        ! the less accurate axis and edge readback from being over-fitted.
        allocate (wsqrt(nr), sample_now(nr + 1), jac_samples(nr + 1, n + 1))
        do i = 1, nr
            wsqrt(i) = (rho(min(nr, i + 1)) - rho(max(1, i - 1)))/2
        end do
        wsqrt = sqrt(wsqrt/sum(wsqrt))
        sample_now = 0; jac_samples = 0
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
        nullify(equations_context%basis)
        if (allocated(basis)) equations_context%basis => basis
        equations_context%calls => calls
        equations_context%candidate => candidate
        equations_context%candidate_law => candidate_law
        nullify(equations_context%conversion)
        if (allocated(conversion)) equations_context%conversion => conversion
        equations_context%exact => exact
        equations_context%failure => failure
        equations_context%ffscale => ffscale
        if (present(forcing)) equations_context%forcing => forcing
        equations_context%fscale => fscale
        equations_context%law => law
        equations_context%mu0 => mu0
        equations_context%n => n
        equations_context%na => na
        equations_context%nr => nr
        equations_context%phi => phi
        equations_context%readback => readback
        nullify(equations_context%rho)
        if (allocated(rho)) equations_context%rho => rho
        nullify(equations_context%sample_now)
        if (allocated(sample_now)) equations_context%sample_now => sample_now
        equations_context%seed => seed
        equations_context%talk => talk
        nullify(equations_context%target)
        if (allocated(target)) equations_context%target => target
        equations_context%tol => tol
        nullify(equations_context%warm_psi)
        if (allocated(warm_psi)) equations_context%warm_psi => warm_psi
        nullify(equations_context%wsqrt)
        if (allocated(wsqrt)) equations_context%wsqrt => wsqrt
            call solve_gs_prescribed_q_equations_callback(probe_x, probe_residual, &
                context=equations_context)
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
        nullify(equations_context%basis)
        if (allocated(basis)) equations_context%basis => basis
        equations_context%calls => calls
        equations_context%candidate => candidate
        equations_context%candidate_law => candidate_law
        nullify(equations_context%conversion)
        if (allocated(conversion)) equations_context%conversion => conversion
        equations_context%exact => exact
        equations_context%failure => failure
        equations_context%ffscale => ffscale
        if (present(forcing)) equations_context%forcing => forcing
        equations_context%fscale => fscale
        equations_context%law => law
        equations_context%mu0 => mu0
        equations_context%n => n
        equations_context%na => na
        equations_context%nr => nr
        equations_context%phi => phi
        equations_context%readback => readback
        nullify(equations_context%rho)
        if (allocated(rho)) equations_context%rho => rho
        nullify(equations_context%sample_now)
        if (allocated(sample_now)) equations_context%sample_now => sample_now
        equations_context%seed => seed
        equations_context%talk => talk
        nullify(equations_context%target)
        if (allocated(target)) equations_context%target => target
        equations_context%tol => tol
        nullify(equations_context%warm_psi)
        if (allocated(warm_psi)) equations_context%warm_psi => warm_psi
        nullify(equations_context%wsqrt)
        if (allocated(wsqrt)) equations_context%wsqrt => wsqrt
            call newton_exact(x0, tol)
            if (talk .and. regularized) print '(a,es10.3)', &
                '  regularization lambda ', lambda
        nullify(equations_context%basis)
        if (allocated(basis)) equations_context%basis => basis
        equations_context%calls => calls
        equations_context%candidate => candidate
        equations_context%candidate_law => candidate_law
        nullify(equations_context%conversion)
        if (allocated(conversion)) equations_context%conversion => conversion
        equations_context%exact => exact
        equations_context%failure => failure
        equations_context%ffscale => ffscale
        if (present(forcing)) equations_context%forcing => forcing
        equations_context%fscale => fscale
        equations_context%law => law
        equations_context%mu0 => mu0
        equations_context%n => n
        equations_context%na => na
        equations_context%nr => nr
        equations_context%phi => phi
        equations_context%readback => readback
        nullify(equations_context%rho)
        if (allocated(rho)) equations_context%rho => rho
        nullify(equations_context%sample_now)
        if (allocated(sample_now)) equations_context%sample_now => sample_now
        equations_context%seed => seed
        equations_context%talk => talk
        nullify(equations_context%target)
        if (allocated(target)) equations_context%target => target
        equations_context%tol => tol
        nullify(equations_context%warm_psi)
        if (allocated(warm_psi)) equations_context%warm_psi => warm_psi
        nullify(equations_context%wsqrt)
        if (allocated(wsqrt)) equations_context%wsqrt => wsqrt
            if (failure /= 0 .and. newton == 0 .and. any(abs(law%pp(:, 2:)) > 0)) &
                call pressure_homotopy()
        else
            ! The secant fallback proposes coefficients from central differences
            ! of full forward solves; nothing it computes is a sensitivity.
        nullify(equations_context%basis)
        if (allocated(basis)) equations_context%basis => basis
        equations_context%calls => calls
        equations_context%candidate => candidate
        equations_context%candidate_law => candidate_law
        nullify(equations_context%conversion)
        if (allocated(conversion)) equations_context%conversion => conversion
        equations_context%exact => exact
        equations_context%failure => failure
        equations_context%ffscale => ffscale
        if (present(forcing)) equations_context%forcing => forcing
        equations_context%fscale => fscale
        equations_context%law => law
        equations_context%mu0 => mu0
        equations_context%n => n
        equations_context%na => na
        equations_context%nr => nr
        equations_context%phi => phi
        equations_context%readback => readback
        nullify(equations_context%rho)
        if (allocated(rho)) equations_context%rho => rho
        nullify(equations_context%sample_now)
        if (allocated(sample_now)) equations_context%sample_now => sample_now
        equations_context%seed => seed
        equations_context%talk => talk
        nullify(equations_context%target)
        if (allocated(target)) equations_context%target => target
        equations_context%tol => tol
        nullify(equations_context%warm_psi)
        if (allocated(warm_psi)) equations_context%warm_psi => warm_psi
        nullify(equations_context%wsqrt)
        if (allocated(wsqrt)) equations_context%wsqrt => wsqrt
            call multiroot_hybrids(solve_gs_prescribed_q_equations_callback, n + 1, &
                x0, solution, stat, &
                                   xtol=tol*0.1_dp, ftol=tol, max_iter=20, &
                                       ctx=equations_context)
        nullify(equations_context%basis)
        if (allocated(basis)) equations_context%basis => basis
        equations_context%calls => calls
        equations_context%candidate => candidate
        equations_context%candidate_law => candidate_law
        nullify(equations_context%conversion)
        if (allocated(conversion)) equations_context%conversion => conversion
        equations_context%exact => exact
        equations_context%failure => failure
        equations_context%ffscale => ffscale
        if (present(forcing)) equations_context%forcing => forcing
        equations_context%fscale => fscale
        equations_context%law => law
        equations_context%mu0 => mu0
        equations_context%n => n
        equations_context%na => na
        equations_context%nr => nr
        equations_context%phi => phi
        equations_context%readback => readback
        nullify(equations_context%rho)
        if (allocated(rho)) equations_context%rho => rho
        nullify(equations_context%sample_now)
        if (allocated(sample_now)) equations_context%sample_now => sample_now
        equations_context%seed => seed
        equations_context%talk => talk
        nullify(equations_context%target)
        if (allocated(target)) equations_context%target => target
        equations_context%tol => tol
        nullify(equations_context%warm_psi)
        if (allocated(warm_psi)) equations_context%warm_psi => warm_psi
        nullify(equations_context%wsqrt)
        if (allocated(wsqrt)) equations_context%wsqrt => wsqrt
            call solve_gs_prescribed_q_equations_callback(solution, residual, &
                context=equations_context)
            converged = failure == 0 .and. maxval(abs(residual)) <= tol
        end if
        if (present(evaluations)) evaluations = calls
        if (present(newton_iterations)) newton_iterations = newton
        status = 2
        if (failure /= 0) then
            status = failure
            return
        end if
        if (.not. converged) return
        if (present(regularization_used)) regularization_used = lambda
        state = candidate; law = candidate_law
        ! Independent interleaved radial probes and doubled arc quadrature
        ! check the q profile, beyond the finite set of fitted moments.
        call check_grid(rho, uniform_grid, rho_check)
        call read_gs_equilibrium(state, law%f_edge, 0.0_dp, size(rho_check), &
                                 2*na, readback, status, law, rho_grid=rho_check)
        if (status /= 0) return
        defect = 0; j = 1
        do i = 1, size(readback%q)
            x = rho_check(i)**2
            if (abs(readback%q(i)/gs_q_table_value(s, q, x) - 1) > defect) j = i
            defect = max(defect, abs(readback%q(i)/gs_q_table_value(s, q, x) - 1))
        end do
        if (talk) print '(a,es10.3,a,f9.6,a,2es14.6)', '  q check: max defect ', &
            defect, &
            ' at s ', rho_check(j)**2, ' q, target ', readback%q(j), &
            gs_q_table_value(s, q, rho_check(j)**2)
        if (present(q_error)) q_error = defect
        if (defect > qtol) then
            status = 9
            return
        end if
        status = 0
    contains
        subroutine pressure_homotopy()
            ! The seed field is built with the mean p'. With a strongly
            ! structured p'(s) (pedestal) the first forward field can have no
            ! readable q profile, so the outer Newton cannot start. Solve the
            ! inverse problem for p' blended from the mean to the full law,
            ! each stage warm-started from the previous one, halving the
            ! blend step on failure.
            real(dp), allocatable :: pp_full(:, :), pp_mean(:, :), x_good(:), psi_good(:)
            real(dp) :: blend, blend_ok, step, stage_target
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
                ! An intermediate stage is only the warm start of the next:
                ! it is accepted at max|r| <= 1e-3. Only the full law must
                ! meet tol. On TC24 (16 graded intervals) the q moments of the
                ! blended laws carry a component along the weakest Jacobian
                ! direction (sigma_min/sigma_max ~ 6e-6, an alternating edge
                ! FF' mode) that no nearby law removes; the stages stalled at
                ! max|r| 2e-5..1e-4 and the homotopy gave up at blend 0.25.
                stage_target = tol
                if (blend < 1) stage_target = max(tol, 1.0e-3_dp)
                call newton_exact(x_good, stage_target)
                stage_ok = failure == 0 .and. converged
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
        subroutine newton_exact(xstart, target)
            ! Levenberg-Marquardt-damped Newton with the equation-level Jacobian,
            ! evaluated only at accepted iterates. A trial is the Newton step
            ! (mu = 0) or the minimizer of |r + J dx| + mu |dx| (mu relative to
            ! the largest diagonal of J^T J). Halving the Newton step instead
            ! kept its direction: on TC24 the edge FF' coefficients are weakly
            ! determined by q, the Newton step there is huge, and every scaled
            ! copy of it left the forward problem's solution set. The forward
            ! solve of each trial starts from the last accepted field. Near
            ! the solution mu returns to zero, so convergence is quadratic.
            ! Regularized: the same iteration on |r|^2 + lambda^2 |S x|^2
            ! (S = smooth_rows), Gauss-Newton steps from stacked least squares,
            ! converged when the undamped step is below target at frozen lambda.
            real(dp), intent(in) :: xstart(:), target
            real(dp), allocatable :: x(:), xtrial(:), r(:), rtrial(:), dx(:), &
                                     jac(:, :), normal(:, :), gradient(:), system(:, :)
            real(dp) :: mu, merit, merit_trial, predicted, diagonal, stall_merit
            real(dp) :: lambda_new, change
            logical :: settled
            real(dp), allocatable :: rs(:), rstrial(:), js(:, :)
            integer(int64) :: tick, tock, rate
            integer :: iteration, trial, info, k, stalled, updates
            logical :: ok, moved, frozen
            allocate (x(n + 1), xtrial(n + 1), r(n + 1), rtrial(n + 1), dx(n + 1), &
                      jac(n + 1, n + 1), normal(n + 1, n + 1), gradient(n + 1), &
                      system(n + 1, n + 1))
            x = xstart
            mu = 0; stalled = 0; stall_merit = 0
            converged = .false.
            frozen = .not. adaptive; updates = 0
            if (lambda_fixed > 0) lambda = lambda_fixed
            if (adaptive) lambda = 0
            call solve_gs_prescribed_q_equations_callback(x, r, &
                context=equations_context)
            solution = x; residual = r
            if (failure /= 0) return
            rs = sample_now; rstrial = rs
            call accept(x, r)
            do iteration = 1, merge(80, 25, adaptive)
                if (.not. regularized .and. maxval(abs(r)) <= target) then
                    converged = .true.
                    exit
                end if
                call system_clock(tick, rate)
                call exact_jacobian(jac, ok)
                call system_clock(tock)
                jacobian_seconds = real(tock - tick, dp)/rate
                if (.not. ok) then
                    failure = 4
                    exit
                end if
                js = jac_samples
                if (regularized) then
                    ! Discrepancy (Morozov) choice on each linearization until
                    ! it settles: the step's linearized q misfit is the larger
                    ! of q_accuracy and half the current misfit (a continuation
                    ! that keeps early steps smooth). lambda is frozen once the
                    ! goal is q_accuracy and two choices agree within 25 %.
                    if (.not. frozen) then
                        call discrepancy_lambda(js, rs, x, lambda_new, settled)
                        if (settled .and. lambda > 0) frozen = lambda_new <= &
                            1.25_dp*lambda .and. lambda <= 1.25_dp*lambda_new
                        updates = updates + 1
                        if (updates >= 40) frozen = .true.
                        if (.not. frozen) lambda = lambda_new
                    end if
                    call stacked_step(js, rs, x, lambda, 0.0_dp, dx, info)
                    if (info /= 0) then
                        failure = 4
                        exit
                    end if
                    ! Largest change the Gauss-Newton step makes in a fitted
                    ! residual or penalty component (relative q units).
                    change = max(maxval(abs(matmul(js, dx))), &
                                 lambda*maxval(abs(matmul(smooth_rows, dx))))
                    if (talk) print '(a,es10.3,a,es10.3,a,es10.3,a,l2)', &
                        '  Gauss-Newton step change ', change, ' lambda ', lambda, &
                        ' rms q misfit ', norm2(rs(:nr)), ' frozen', frozen
                    if (frozen .and. change <= target) then
                        converged = .true.
                        exit
                    end if
                end if
                newton = iteration
                merit = 0.5_dp*dot_product(r, r)
                if (regularized) merit = 0.5_dp*(dot_product(rs, rs) &
                                         + lambda**2*sum(matmul(smooth_rows, x)**2))
                normal = matmul(transpose(jac), jac)
                if (regularized) normal = matmul(transpose(js), js)
                gradient = matmul(transpose(jac), r)
                diagonal = maxval([(normal(k, k), k=1, n + 1)])
                moved = .false.
                call system_clock(tick)
                do trial = 1, 16
                    if (regularized) then
                        call stacked_step(js, rs, x, lambda, mu*diagonal, dx, info)
                    else if (mu > 0) then
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
                    call solve_gs_prescribed_q_equations_callback(xtrial, rtrial, &
                        context=equations_context)
                    if (failure == 0) then
                        merit_trial = 0.5_dp*dot_product(rtrial, rtrial)
                        predicted = merit - 0.5_dp*sum((r + matmul(jac, dx))**2)
                        if (regularized) then
                            rstrial = sample_now
                            merit_trial = 0.5_dp*(dot_product(rstrial, rstrial) + &
                                lambda**2*sum(matmul(smooth_rows, xtrial)**2))
                            predicted = merit - 0.5_dp*(sum((rs + matmul(js, dx))**2) &
                                        + lambda**2*sum(matmul(smooth_rows, xtrial)**2))
                        end if
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
                if (talk) print '(a,i3,a,es10.3,a,es10.3,a,f9.2,a,f9.2,a,es9.2)', &
                    'inverse newton ', iteration, ' max|r|=', maxval(abs(r)), &
                    ' max|dx|=', maxval(abs(dx)), &
                    ' jacobian_s=', jacobian_seconds, ' forward_s=', solve_seconds, &
                    ' mu=', mu
                if (talk) flush (output_unit)
                if (.not. moved) then
                    failure = 0
                    exit
                end if
                ! Three consecutive damped steps that together reduce the
                ! merit by less than 1%: a nonzero residual minimum (singular
                ! J, e.g. the unregularized TC24 floor max|r| 5.5e-4); stop.
                if (mu > 0) then
                    if (stalled == 0) stall_merit = merit
                    stalled = stalled + 1
                else
                    stalled = 0
                end if
                mu = mu/10
                if (mu < 1.0e-6_dp) mu = 0
                x = xtrial; r = rtrial; rs = rstrial
                call accept(x, r)
                if (stalled == 3) then
                    if (merit_trial > 0.99_dp*stall_merit) exit
                    stalled = 0
                end if
            end do
            ! Report the last accepted iterate.
            if (failure == 0 .or. failure == 4) then
                candidate = accepted_state; candidate_law = accepted_law
            end if
            solution = x; residual = r
            if (failure == 4) return
            failure = 0
        end subroutine
        subroutine stacked_step(jac, r, x, weight, damping, dx, info)
            ! Least-squares step of [J; weight S; sqrt(damping) I] dx =
            ! [-r; -weight S x; 0] (LAPACK QR; no normal equations).
            real(dp), intent(in) :: jac(:, :), r(:), x(:), weight, damping
            real(dp), intent(out) :: dx(:)
            integer, intent(out) :: info
            real(dp), allocatable :: a(:, :), b(:), work(:)
            real(dp) :: query(1)
            integer :: m, p, k, q
            p = size(smooth_rows, 1); q = size(jac, 1); m = q + p + (n + 1)
            allocate (a(m, n + 1), b(m))
            a = 0; b = 0
            a(:q, :) = jac; b(:q) = -r
            a(q + 1:q + p, :) = weight*smooth_rows
            b(q + 1:q + p) = -weight*matmul(smooth_rows, x)
            do k = 1, n + 1
                a(q + p + k, k) = sqrt(damping)
            end do
            call dgels('N', m, n + 1, 1, a, m, b, m, query, -1, info)
            allocate (work(max(1, int(query(1)))))
            call dgels('N', m, n + 1, 1, a, m, b, m, work, size(work), info)
            dx = b(:n + 1)
            if (info == 0 .and. .not. all(ieee_is_finite(dx))) info = 1
        end subroutine
        subroutine discrepancy_lambda(jac, r, x, weight, settled)
            ! lambda whose linearized step at x has q misfit |r_q + J_q dx| =
            ! max(q_accuracy, |r_q|/2) (at least 1.01 times the unregularized
            ! floor), by bisection in log lambda; the misfit grows with lambda.
            ! settled: the goal is q_accuracy or, when the data cannot be met
            ! that closely (an intermediate homotopy law), the floor.
            real(dp), intent(in) :: jac(:, :), r(:), x(:)
            real(dp), intent(out) :: weight
            logical, intent(out) :: settled
            real(dp) :: lo, hi, mid, goal, floor_misfit, step(n + 1), base
            integer :: k, st
            base = norm2(jac(:, :n))/max(norm2(smooth_rows), tiny(1.0_dp))
            hi = log(base) + 6*log(10.0_dp); lo = log(base) - 12*log(10.0_dp)
            call stacked_step(jac, r, x, exp(lo), 0.0_dp, step, st)
            floor_misfit = norm2(r(:nr) + matmul(jac(:nr, :), step))
            goal = max(misfit_target, 0.5_dp*norm2(r(:nr)), 1.01_dp*floor_misfit)
            settled = 0.5_dp*norm2(r(:nr)) <= max(misfit_target, 1.01_dp*floor_misfit)
            do k = 1, 80
                mid = (lo + hi)/2
                call stacked_step(jac, r, x, exp(mid), 0.0_dp, step, st)
                if (norm2(r(:nr) + matmul(jac(:nr, :), step)) > goal) then
                    hi = mid
                else
                    lo = mid
                end if
                if (hi - lo < 1.0e-3_dp) exit
            end do
            weight = exp(lo)
            if (talk) print '(a,es10.3,a,es10.3,a,l2)', '  discrepancy lambda ', &
                weight, ' goal ', goal, ' settled', settled
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
            !     other, because they share the cached native Cholesky factor
            !     of A (retained factors own mutable workspaces) and phase timers.
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
                    if (load_ok(kk)) call column_readback(kk, base, &
                        dpsi_block(:, b), jac(:, kk), jac_samples(:, kk), &
                        column_ok(kk), tb)
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
            jac_samples(:nr, n + 1) = wsqrt*(rplus%q - rminus%q)/(2.0e-3_dp*target)
            jac_samples(nr + 1, n + 1) = jac(n + 1, n + 1)
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
        subroutine column_readback(kk, base, dpsi, col, col_samples, good_col, tb)
            ! Central-difference readback of the k-th column along dpsi.
            integer, intent(in) :: kk
            type(gs_profile_law_t), intent(in) :: base
            real(dp), intent(in) :: dpsi(:)
            real(dp), intent(out) :: col(:), col_samples(:)
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
            col = 0; col_samples = 0
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
            col_samples(:nr) = wsqrt*(rplus%q - rminus%q)/(2*eps*target)
            col_samples(nr + 1) = col(n + 1)
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
        subroutine solve_gs_prescribed_q_equations_callback(coefficients, values, &
            context)
            real(dp), intent(in) :: coefficients(:)
            real(dp), intent(out) :: values(:)
            class(*), intent(in), optional  :: context
            integer :: k, stat
        type(solve_gs_prescribed_q_equations_context_t) :: callback_context
        if (.not. present(context)) &
            error stop "missing private callback context"
        select type(selected_context => context)
        type is (solve_gs_prescribed_q_equations_context_t)
            callback_context = selected_context
        class default
            error stop "invalid private callback context"
        end select
            callback_context%calls = callback_context%calls + 1; values = 1.0e3_dp; &
                callback_context%failure = 0
            if (.not. all(ieee_is_finite(coefficients))) then
                callback_context%failure = 3
                return
            end if
            if (abs(coefficients(callback_context%n + 1)) > 100) then
                callback_context%failure = 3
                return
            end if
            callback_context%candidate = callback_context%seed; &
                callback_context%candidate_law = callback_context%law
            if (callback_context%exact) callback_context%candidate%psi = &
                callback_context%warm_psi
            callback_context%candidate_law%ff = 0
            do k = 1, callback_context%n
                callback_context%candidate_law%ff = callback_context%candidate_law%ff &
                                   + callback_context%ffscale*coefficients(k)*callback_context%conversion(:, :, k)
            end do
            call join_cubics(callback_context%candidate_law%knots, &
                callback_context%candidate_law%ff)
            callback_context%candidate_law%f_edge = &
                sign(callback_context%fscale*exp(coefficients(callback_context%n + &
                1)), callback_context%phi)
            call solve_gs_cubic_profiles(callback_context%candidate, &
                callback_context%mu0, callback_context%candidate_law, stat, &
                           min(1.0e-11_dp, callback_context%tol*0.01_dp), &
                               forcing=callback_context%forcing, &
                                         quadrature_tolerance=1.0e-9_dp)
            if (stat /= 0) then
                callback_context%failure = 100 + stat
                if (callback_context%talk) then
                    print '(a,i0,a,es10.3)', '  forward solve status ', stat, &
                        ' residual ', callback_context%candidate%residual
                    flush (output_unit)
                end if
                return
            end if
            call read_gs_equilibrium(callback_context%candidate, &
                callback_context%candidate_law%f_edge, 0.0_dp, &
                              callback_context%nr, callback_context%na, &
                                  callback_context%readback, stat, &
                                  callback_context%candidate_law, &
                                  rho_grid=callback_context%rho)
            if (stat /= 0) then
                callback_context%failure = 200 + stat
                return
            end if
            values(:callback_context%n) = matmul(callback_context%basis, &
                callback_context%readback%q/callback_context%target - 1)
            values(callback_context%n + 1) = &
                callback_context%readback%phi_edge/callback_context%phi - 1
            callback_context%sample_now(:callback_context%nr) = &
                callback_context%wsqrt*(callback_context%readback%q/callback_context%target - 1)
            callback_context%sample_now(callback_context%nr + 1) = &
                values(callback_context%n + 1)
        end subroutine
end module
