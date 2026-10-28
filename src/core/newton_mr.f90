module kin6d_newton_mr
    ! Matrix-free inexact Newton for a smooth objective with a symmetric,
    ! possibly indefinite or singular Hessian: Newton-MR (Liu and Roosta,
    ! "A Newton-MR algorithm with complexity guarantees for nonconvex smooth
    ! unconstrained optimization", 2022), with preconditioned MINRES as the
    ! inner solver (Paige and Saunders 1975; residual recurrence and the
    ! nonpositive-curvature test of Liu and Roosta, "MINRES: from negative
    ! curvature detection to monotonicity properties", SIAM J. Optim. 2022).
    !
    ! MINRES solves H d = -g to the relative residual eta (Euclidean norm of
    ! the true residual, carried by its recurrence and checked once at the
    ! end). With npc_exit it stops at the first Krylov residual r with
    ! <M r, H M r> <= 0 and returns d = M r instead (an NPC direction, a
    ! descent direction for the Newton system it came from). A Newton-MR
    ! step on a "solution" (SOL) direction backtracks from alpha = 1 under
    ! the Armijo rule on the objective; on an NPC direction it tracks forward
    ! (doubling alpha while the Armijo rule holds) or backward. When the
    ! objective change is below its roundoff band, steps are ranked by the
    ! trapezoidal change (g(x) + g(x + s)) . s / 2, as in kin6d_ideal_map.
    !
    ! With merit = MR_MERIT_GRADIENT the line search is on |g|^2 / 2 instead
    ! (the original Newton-MR of Roosta, Liu, Xu and Mahoney, "Newton-MR:
    ! inexact Newton method with minimum residual sub-problem solver",
    ! EURO J. Comput. Optim. 2022, for invex objectives). This is the merit
    ! for stationary points that are not minimisers, e.g. the saddle points
    ! of a Lagrangian with multipliers: every MINRES iterate d satisfies
    ! g . H d = (|H d + g|^2 - |H d|^2 - |g|^2)/2 < 0 (MINRES residuals do not
    ! increase), so d descends |g|^2 / 2 without any curvature condition, and
    ! the NPC exit is not used.
    !
    ! Callbacks are type-bound procedures of caller-owned objects (no
    ! internal-procedure actuals, hence no executable-stack trampolines).
    ! The generic Krylov kernel is a candidate for FortNum (fortnum_krylov
    ! has CG, GMRES and FGMRES but no MINRES); it lives here until moved.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    private
    public :: mr_operator_t, mr_problem_t, mr_options_t, mr_report_t
    public :: minres_solve, newton_mr_solve
    integer, parameter, public :: MR_OK = 0, MR_MAX_ITERATIONS = 1, &
        MR_BREAKDOWN = 2, MR_NPC = 3, MR_LINE_SEARCH = 4, MR_BAD_INPUT = 5, &
        MR_EVALUATION = 6
    ! Forcing: eta = min(forcing_max, (|g|/|g0|)^(1/2)) as the kin6d_ideal_map
    ! krylov_step (Eisenstat-Walker type), or the fixed eta = forcing_max.
    integer, parameter, public :: MR_FORCING_SQRT = 0, MR_FORCING_FIXED = 1
    ! Line-search merit: the objective f, or |g|^2 / 2 (stationary points).
    integer, parameter, public :: MR_MERIT_OBJECTIVE = 0, MR_MERIT_GRADIENT = 1

    type, abstract :: mr_operator_t
        ! Symmetric operator y = A x and SPD preconditioner z = M r ~ A^-1 r.
        ! The default preconditioner is the identity.
    contains
        procedure(mr_apply_i), deferred :: apply
        procedure :: precondition => identity_precondition
    end type

    type, abstract, extends(mr_operator_t) :: mr_problem_t
        ! Objective f, gradient g and, after linearize(x), the Hessian
        ! action apply and its preconditioner at x.
    contains
        procedure(mr_evaluate_i), deferred :: evaluate
        procedure :: linearize => no_linearize
    end type

    abstract interface
        subroutine mr_apply_i(self, x, y)
            import :: mr_operator_t, dp
            class(mr_operator_t), intent(inout) :: self
            real(dp), intent(in) :: x(:)
            real(dp), intent(out) :: y(:)
        end subroutine
        subroutine mr_evaluate_i(self, x, f, g, status)
            import :: mr_problem_t, dp
            class(mr_problem_t), intent(inout) :: self
            real(dp), intent(in) :: x(:)
            real(dp), intent(out) :: f, g(:)
            integer, intent(out) :: status
        end subroutine
    end interface

    type :: mr_options_t
        ! Converged when |g| <= max(tolerance |g0|, absolute_tolerance).
        real(dp) :: tolerance = 1.0e-10_dp, absolute_tolerance = 0
        integer :: max_iterations = 50, max_minres = 1000
        integer :: forcing = MR_FORCING_SQRT
        real(dp) :: forcing_max = 0.5_dp
        ! Armijo constant and the maximum line-search trials per step.
        real(dp) :: armijo = 1.0e-4_dp
        integer :: max_line_search = 40
        logical :: npc_exit = .true.
        integer :: merit = MR_MERIT_OBJECTIVE
    end type

    type :: mr_report_t
        integer :: status = MR_BAD_INPUT
        integer :: iterations = 0, minres_iterations = 0, minres_max = 0
        integer :: npc_steps = 0, backtracks = 0, evaluations = 0
        real(dp) :: energy = 0, initial_gradient = 0, gradient = 0
        ! Wall seconds: linearize (Hessian point and preconditioner setup),
        ! MINRES (operator and preconditioner applications), evaluations.
        real(dp) :: linearize_s = 0, minres_s = 0, evaluate_s = 0
    end type
contains
    subroutine identity_precondition(self, r, z, status)
        class(mr_operator_t), intent(inout) :: self
        real(dp), intent(in) :: r(:)
        real(dp), intent(out) :: z(:)
        integer, intent(out) :: status
        z = r
        status = 0
    end subroutine

    subroutine no_linearize(self, x, status)
        class(mr_problem_t), intent(inout) :: self
        real(dp), intent(in) :: x(:)
        integer, intent(out) :: status
        status = 0
    end subroutine

    subroutine minres_solve(op, b, x, rtol, max_iterations, info, iterations, &
                            residual, npc_exit)
        ! Preconditioned MINRES for A x = b from x = 0. Stops when the
        ! recurred true residual satisfies |b - A x| <= rtol |b|. With npc_exit
        ! a residual r_(k-1) of nonpositive curvature, <M r, A M r> <= 0 (the
        ! test c_(k-1) gamma1_k >= 0 of the QR recurrence), ends the solve
        ! with info = MR_NPC and x = M r_(k-1). residual is |b - A x| / |b|,
        ! recomputed explicitly for the returned solution (not for an NPC
        ! direction). A preconditioner with r . M r < 0 or a failed apply is
        ! MR_BREAKDOWN. A singular A with consistent b is admitted: the
        ! Krylov space stays in the range and x is the M-minimum-norm
        ! solution component.
        class(mr_operator_t), intent(inout) :: op
        real(dp), intent(in) :: b(:), rtol
        real(dp), intent(out) :: x(:)
        integer, intent(in) :: max_iterations
        integer, intent(out) :: info, iterations
        real(dp), intent(out) :: residual
        logical, intent(in), optional :: npc_exit
        real(dp), allocatable :: v(:), y(:), r1(:), r2(:), w(:), w1(:), w2(:)
        real(dp), allocatable :: res(:)
        real(dp) :: beta, oldb, alpha, delta, gbar, epsln, oldeps, dbar
        real(dp) :: gamma, cs, sn, phi, phibar, bnorm, rnorm, target
        integer :: n, k, st
        logical :: npc_test

        n = size(b)
        npc_test = .false.
        if (present(npc_exit)) npc_test = npc_exit
        x = 0; iterations = 0; residual = 0; info = MR_BAD_INPUT
        if (size(x) /= n .or. max_iterations < 0 .or. .not. rtol >= 0) return
        bnorm = norm2(b)
        info = MR_OK
        if (.not. bnorm > 0) return
        allocate (v(n), y(n), r1(n), r2(n), w(n), w1(n), w2(n), res(n))
        r1 = 0; r2 = b; res = b; w = 0; w1 = 0
        call op%precondition(r2, y, st)
        beta = dot_product(r2, y)
        info = MR_BREAKDOWN
        if (st /= 0 .or. .not. beta >= 0 .or. .not. ieee_is_finite(beta)) return
        beta = sqrt(beta)
        oldb = 0; dbar = 0; epsln = 0; phibar = beta; cs = -1; sn = 0
        target = rtol*bnorm
        info = MR_MAX_ITERATIONS
        do k = 1, max_iterations
            if (.not. beta > 0) then
                ! Invariant Krylov space: the last iterate solves the
                ! (consistent) system exactly.
                info = MR_OK; exit
            end if
            v = y/beta
            call op%apply(v, y)
            if (k > 1) y = y - (beta/oldb)*r1
            alpha = dot_product(v, y)
            y = y - (alpha/beta)*r2
            r1 = r2; r2 = y
            call op%precondition(r2, y, st)
            oldb = beta
            beta = dot_product(r2, y)
            if (st /= 0 .or. .not. beta >= 0 .or. .not. ieee_is_finite(beta)) then
                info = MR_BREAKDOWN; exit
            end if
            beta = sqrt(beta)
            oldeps = epsln
            delta = cs*dbar + sn*alpha
            gbar = sn*dbar - cs*alpha
            epsln = sn*beta
            dbar = -cs*beta
            if (npc_test .and. cs*gbar >= 0) then
                ! The residual of x_(k-1) has nonpositive curvature.
                call op%precondition(res, x, st)
                info = merge(MR_NPC, MR_BREAKDOWN, st == 0)
                iterations = k - 1
                return
            end if
            gamma = sqrt(gbar**2 + beta**2)
            if (.not. gamma > 0) then
                info = MR_BREAKDOWN; exit
            end if
            cs = gbar/gamma; sn = beta/gamma
            phi = cs*phibar; phibar = sn*phibar
            w2 = w1; w1 = w
            w = (v - oldeps*w2 - delta*w1)/gamma
            x = x + phi*w
            iterations = k
            if (beta > 0) then
                res = sn**2*res - (cs*phibar/beta)*r2
            else
                res = 0
            end if
            rnorm = norm2(res)
            if (.not. ieee_is_finite(rnorm)) then
                info = MR_BREAKDOWN; exit
            end if
            if (rnorm <= target) then
                info = MR_OK; exit
            end if
        end do
        ! Explicit residual of the returned iterate (one extra product).
        call op%apply(x, y)
        residual = norm2(b - y)/bnorm
    end subroutine

    subroutine newton_mr_solve(problem, x, options, report)
        ! Newton-MR from x. Each iteration linearizes the problem at x,
        ! solves H d = -g by MINRES to the forcing eta and takes a line search
        ! along d. report%status: MR_OK (gradient tolerance met), MR_MAX_
        ! ITERATIONS, MR_LINE_SEARCH (no acceptable step), MR_BREAKDOWN
        ! (MINRES/preconditioner failure or a non-descent direction that the
        ! preconditioned gradient could not replace) or MR_EVALUATION.
        class(mr_problem_t), intent(inout) :: problem
        real(dp), intent(inout) :: x(:)
        type(mr_options_t), intent(in) :: options
        type(mr_report_t), intent(out) :: report
        real(dp), allocatable :: g(:), gt(:), d(:), xt(:), mg(:)
        real(dp) :: f, ft, gnorm, g0, eta, gd, alpha, gtol, residual, t0
        real(dp) :: f_best, alpha_best
        integer :: it, info, its, st, ls, n, forward
        logical :: npc, accepted, gradient_merit

        n = size(x)
        gradient_merit = options%merit == MR_MERIT_GRADIENT
        allocate (g(n), gt(n), d(n), xt(n), mg(n))
        report%status = MR_EVALUATION
        t0 = clock()
        call problem%evaluate(x, f, g, st)
        report%evaluate_s = report%evaluate_s + clock() - t0
        report%evaluations = 1
        if (st /= 0 .or. .not. ieee_is_finite(f)) return
        g0 = norm2(g); gnorm = g0
        report%initial_gradient = g0; report%energy = f; report%gradient = g0
        gtol = max(options%tolerance*g0, options%absolute_tolerance)
        report%status = MR_MAX_ITERATIONS
        do it = 1, options%max_iterations
            if (gnorm <= gtol) then
                report%status = MR_OK; exit
            end if
            report%iterations = it
            t0 = clock()
            call problem%linearize(x, st)
            report%linearize_s = report%linearize_s + clock() - t0
            if (st /= 0) then
                report%status = MR_BREAKDOWN; exit
            end if
            eta = options%forcing_max
            if (options%forcing == MR_FORCING_SQRT) &
                eta = min(eta, sqrt(gnorm/max(g0, tiny(1.0_dp))))
            ! A forcing far below what the stopping test needs only adds
            ! iterations at the roundoff floor.
            eta = max(eta, 0.5_dp*gtol/gnorm)
            t0 = clock()
            call minres_solve(problem, -g, d, eta, options%max_minres, info, its, &
                residual, options%npc_exit .and. .not. gradient_merit)
            report%minres_s = report%minres_s + clock() - t0
            report%minres_iterations = report%minres_iterations + its
            report%minres_max = max(report%minres_max, its)
            if (info == MR_BREAKDOWN .or. info == MR_BAD_INPUT) then
                report%status = MR_BREAKDOWN; exit
            end if
            ! MR_MAX_ITERATIONS: the last MINRES iterate still has a smaller
            ! residual than d = 0 and is used as an inexact step.
            npc = info == MR_NPC
            if (npc) report%npc_steps = report%npc_steps + 1
            if (gradient_merit) then
                ! Slope of |g|^2 / 2 along d: g . H d (one extra product).
                call problem%apply(d, mg)
                gd = dot_product(g, mg)
                if (.not. gd < 0) then
                    report%status = MR_BREAKDOWN; exit
                end if
            else
                gd = dot_product(g, d)
            end if
            if (.not. gd < 0) then
                ! Not a descent direction (only possible through roundoff or
                ! an unconverged NPC exit at k = 1): preconditioned gradient.
                call problem%precondition(-g, mg, st)
                gd = dot_product(g, mg)
                if (st /= 0 .or. .not. gd < 0) then
                    report%status = MR_BREAKDOWN; exit
                end if
                d = mg; npc = .true.
            end if
            accepted = .false.
            alpha = 1
            f_best = huge(1.0_dp); alpha_best = 0
            do ls = 1, options%max_line_search
                xt = x + alpha*d
                t0 = clock()
                call problem%evaluate(xt, ft, gt, st)
                report%evaluate_s = report%evaluate_s + clock() - t0
                report%evaluations = report%evaluations + 1
                if (st == 0 .and. armijo_ok(f, ft, alpha*gd, g, gt, alpha*d)) then
                    accepted = .true.
                    f_best = ft; alpha_best = alpha
                    ! NPC: forward tracking while the Armijo rule holds.
                    if (.not. (npc .and. ls == 1)) exit
                    do forward = 1, options%max_line_search
                        alpha = 2*alpha
                        xt = x + alpha*d
                        t0 = clock()
                        call problem%evaluate(xt, ft, gt, st)
                        report%evaluate_s = report%evaluate_s + clock() - t0
                        report%evaluations = report%evaluations + 1
                        if (st /= 0) exit
                        if (.not. armijo_ok(f, ft, alpha*gd, g, gt, alpha*d)) exit
                        if (.not. ft < f_best) exit
                        f_best = ft; alpha_best = alpha
                    end do
                    exit
                end if
                report%backtracks = report%backtracks + 1
                alpha = alpha/2
            end do
            if (.not. accepted) then
                report%status = MR_LINE_SEARCH; exit
            end if
            x = x + alpha_best*d
            t0 = clock()
            call problem%evaluate(x, f, g, st)
            report%evaluate_s = report%evaluate_s + clock() - t0
            report%evaluations = report%evaluations + 1
            if (st /= 0 .or. .not. ieee_is_finite(f)) then
                report%status = MR_EVALUATION; exit
            end if
            gnorm = norm2(g)
            report%energy = f; report%gradient = gnorm
        end do
        if (gnorm <= gtol) report%status = MR_OK
        report%energy = f; report%gradient = gnorm
    contains
        logical function armijo_ok(f0, f1, slope, g0v, g1v, s) result(ok)
            ! Armijo sufficient decrease; inside the roundoff band of f the
            ! trapezoidal change ranks the step instead. Gradient merit:
            ! |g1|^2 <= |g0|^2 + 2 armijo slope.
            real(dp), intent(in) :: f0, f1, slope, g0v(:), g1v(:), s(:)
            real(dp) :: band
            ok = .false.
            if (.not. ieee_is_finite(f1)) return
            if (gradient_merit) then
                ok = 0.5_dp*sum(g1v**2) <= 0.5_dp*sum(g0v**2) + options%armijo*slope
                return
            end if
            band = 64*epsilon(1.0_dp)*max(abs(f0), abs(f1))
            if (abs(f1 - f0) > band) then
                ok = f1 <= f0 + options%armijo*slope
            else
                ok = 0.5_dp*dot_product(g0v + g1v, s) <= options%armijo*slope
            end if
        end function
    end subroutine

    real(dp) function clock()
        integer(8) :: count, rate
        call system_clock(count, rate)
        clock = real(count, dp)/real(rate, dp)
    end function
end module kin6d_newton_mr
