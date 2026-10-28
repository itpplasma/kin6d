submodule(kin6d_nonnested_beltrami) kin6d_nonnested_beltrami_mr
    ! Opt-in matrix-free MINRES path of the Taylor closure
    ! (docs/equilibrium-beltrami.md, docs/nonnested-newton-mr.md). After
    ! lifting the boundary trace the fixed-mu equation is the symmetric
    ! indefinite saddle system
    !
    !     [ S - mu H_s   G ] [A]   [b(mu)]        S   = C^T M2(K) C
    !     [ G^T          0 ] [p] = [c    ],       H_s = (P C + C^T P^T)/2
    !                                             G   = M1 grad_0 (n = 0 only)
    !
    ! on the free dofs (P C is symmetric on them to roundoff; the factored
    ! path symmetrises the assembled block the same way, while the lifted
    ! trace enters through P C itself). S is never assembled: its action is
    ! FortFEM's sum-factorised meridian_fourier_curl_curl_apply (FFT in phi per
    ! meridional quadrature point), as in kin6d_nonnested_newton_mr. C, P and
    ! G are mode-block-diagonal meridional matrices and stay sparse.
    !
    ! Preconditioner (BELTRAMI_PC_BLOCK): the block-diagonal SPD
    !
    !     P^-1 = blockdiag( Sbar + tau M1,  L0 / tau ),   L0 = grad_0^T M1 grad_0,
    !
    ! with Sbar the phi-averaged (metric-lumped) curl-curl. Both Sbar and M1
    ! commute with the Fourier basis, so the A block is one meridional
    ! factor per mode n on its (cos, sin) pair. Justification:
    ! - The multiplier Schur complement is exact for every tau > 0: grad V0
    !   is the kernel of C, so (S + tau M1) grad q = tau M1 grad q and
    !   G^T (S + tau M1)^-1 G = L0 / tau, the scalar Laplacian (Greif and
    !   Schoetzau, Numer. Linear Algebra Appl. 14 (2007) 281, for the mixed
    !   Maxwell saddle system). On grad V0 x V0 the preconditioned saddle has
    !   the eigenvalues +-1.
    ! - On the Coulomb complement the A block K = S - mu H_s satisfies
    !   |a.H_s a| <= |a| |curl a| <= C_P |curl a|^2 (discrete Poincare,
    !   uniform in h for Nedelec), so S^-1 K = I - mu S^-1 H_s is the
    !   identity plus a compact, h-uniformly bounded perturbation (mapped:
    !   plus the bounded lumping error of Sbar): it clusters at 1, the negative
    !   part consists of the finitely many Beltrami eigenvalues below |mu|
    !   (the inertia of K). MINRES then needs iterations bounded in h and n,
    !   growing with the number of those eigenvalues and with the inverse
    !   distance of mu to the nearest one (a resonance is the only way to
    !   approach zero; no inertia is available on this path).
    ! - tau only regularises the gradient kernel of the n = 0 block, so it is
    !   small and mu-independent: tau = 0.01/a^2 (a the meridional half-
    !   width). The H(curl) choice tau ~ mu^2 squeezes the low curl modes
    !   towards zero and was measured to need 2-5 times more iterations; the
    !   factors are reused unchanged when mu moves (Newton on mu).
    ! Mapped domains: Sbar ~ S with the constants of sup K/Kbar and
    ! sup Kbar/K (the metric lumping of kin6d_nonnested_newton_mr).
    !
    ! Helicity constraint: H(A) = H0 makes mu an unknown. The Taylor
    ! Lagrangian
    !     L(A, p, mu) = W(A) + p.G^T A - mu (H(A) - H0) / (4 pi)
    ! (W the energy form of the equation, H = 2 pi (A.P C A + 2 A.P h)) has
    ! the stationary equations "fixed-mu system" and "H = H0", and the
    ! symmetric indefinite bordered Hessian [[K, G, -g], [G^T, 0, 0],
    ! [-g^T, 0, 0]], g = H_s A + P h. Its stationary point is solved by
    ! Newton-MR with the |grad L|^2 merit (MR_MERIT_GRADIENT); the mu block
    ! of the preconditioner is the Schur approximation 1 / (g^T P_A g).
    use fortfem_meridian_feec_n, only: meridian_lagrange_nodes
    use fortfem_meridian_fourier_apply, only: meridian_fourier_op_init, &
        meridian_fourier_curl_curl_apply, meridian_fourier_op_bytes
    use fortfem_meridian_fourier_curl, only: meridian_fourier_mass
    use kin6d_newton_mr, only: mr_problem_t, mr_report_t, minres_solve, &
        newton_mr_solve, MR_OK
    use kin6d_equilibrium_nonnested, only: nn_factor_t, nn_factor, &
        nn_factor_solve_real, nn_factor_free
    use kin6d_nonnested_newton_mr, only: lumped_metric_t, mode_curl, mode_mass
    implicit none

    type, extends(mr_problem_t) :: saddle_t
        ! x = [A on the free dofs (system%index order), p on the interior
        ! V0 nodes, mu if bordered].
        type(beltrami_system_t), pointer :: sys => null()
        integer :: kind = BELTRAMI_PC_BLOCK, nf = 0, ng = 0, n = 0
        real(dp) :: mu = 0, tau = 1, schur = 1, target = 0
        logical :: bordered = .false.
        integer, allocatable :: free(:), nodes(:)
        ! Equation load in the full layout: b(A) = load0 + mu load1.
        real(dp), allocatable :: load0(:), load1(:), border(:), diag(:)
        type(nn_factor_t), allocatable :: block(:)
        logical, allocatable :: active(:)
        type(nn_factor_t) :: laplace
    contains
        procedure :: apply => saddle_apply
        procedure :: precondition => saddle_precondition
        procedure :: evaluate => saddle_evaluate
        procedure :: linearize => saddle_linearize
    end type
contains
    module subroutine beltrami_mr_prepare(system, field, metric, nphi, status)
        ! Energy action data, transposes, and the per-mode lumped curl-curl
        ! Sbar_n, V1 mass M_n and gauge Laplacian L0 of the preconditioner.
        type(beltrami_system_t), intent(inout) :: system
        type(free_b_field_t), intent(in) :: field
        type(map_metric_t), intent(in), target :: metric
        integer, intent(in) :: nphi
        integer, intent(out) :: status
        type(lumped_metric_t) :: lumped
        type(fortsparse_status_t) :: fst
        type(csc_t) :: m0, cn, cnt, mn, mc, gt
        real(dp), allocatable :: nodes(:, :)
        real(dp) :: t0
        integer :: n, st
        status = 1
        t0 = wall()
        call meridian_fourier_op_init(field%fs, 2, system%modes, nphi, metric, &
            system%op, st)
        if (st /= 0) return
        system%operator_bytes = meridian_fourier_op_bytes(system%op)
        call csc_transpose(system%curl, system%curl_t, fst)
        if (.not. status_ok(fst)) return
        call csc_transpose(system%helicity, system%helicity_t, fst)
        if (.not. status_ok(fst)) return
        call csc_transpose(system%gauge, system%gauge_t, fst)
        if (.not. status_ok(fst)) return
        call csc_transpose(field%fs%grad%a, gt, fst)
        if (.not. status_ok(fst)) return
        call csc_matmul(gt, system%gauge, system%laplace, fst)
        if (.not. status_ok(fst)) return
        lumped%base => metric; lumped%nphi = nphi
        call meridian_fourier_mass(field%fs, 2, 0, 1, lumped, m0, st)
        if (st /= 0) return
        allocate (system%lumped(0:system%modes), system%mass(0:system%modes))
        do n = 0, system%modes
            if (mod(n, system%nfp) /= 0) cycle
            call mode_curl(field%fs, n, cn, st)
            if (st /= 0) return
            call mode_mass(m0, n, mn, st)
            if (st /= 0) return
            call csc_matmul(mn, cn, mc, fst)
            if (.not. status_ok(fst)) return
            call csc_transpose(cn, cnt, fst)
            if (.not. status_ok(fst)) return
            call csc_matmul(cnt, mc, system%lumped(n), fst)
            if (.not. status_ok(fst)) return
            call mode_mass(field%fs%mass(1), n, system%mass(n), st)
            if (st /= 0) return
        end do
        call meridian_lagrange_nodes(field%fs, nodes)
        system%length2 = (0.5_dp*max(maxval(nodes(1, :)) - minval(nodes(1, :)), &
            maxval(nodes(2, :)) - minval(nodes(2, :))))**2
        if (.not. system%length2 > 0) return
        system%setup_s = wall() - t0
        status = 0
    end subroutine

    module subroutine beltrami_mr_solve(system, mu, field, report, source, &
        potential, dmu)
        ! Fixed-mu solve by preconditioned MINRES; same outputs and
        ! acceptance (relative residual <= 1e-9) as the factored path. A
        ! MINRES failure is reported as BELTRAMI_RESONANT without changing B.
        type(beltrami_system_t), intent(in), target :: system
        real(dp), intent(in) :: mu
        type(free_b_field_t), intent(inout) :: field
        type(beltrami_report_t), intent(inout) :: report
        real(dp), optional, intent(in) :: source(:)
        real(dp), allocatable, optional, intent(out) :: potential(:), dmu(:)
        type(saddle_t) :: problem
        real(dp), allocatable :: x(:), full(:), rhs(:), dx(:), bh(:), sx(:)
        real(dp) :: residual
        integer :: j, info, st
        report%solver = 'minres'
        if (present(source)) then
            if (size(source) /= size(system%lift)) return
            if (.not. all(ieee_is_finite(source))) return
            do j = 1, system%modes
                if (mod(j, system%nfp) == 0) cycle
                if (any(abs(source((2*j - 1)*system%n1 + 1: &
                    (2*j + 1)*system%n1)) > 0)) return
            end do
        end if
        call saddle_init(problem, system, mu, .false., report, st, source)
        if (st /= 0) then
            call saddle_free(problem); return
        end if
        call fixed_solve(problem, mu, x, info, report)
        if (info /= MR_OK .or. .not. report%residual <= max(1.0e-9_dp, &
            10*system%mr%tolerance)) then
            report%status = BELTRAMI_RESONANT
            report%message = 'MINRES did not reach the tolerance (near-resonant '// &
                'or ill-conditioned Beltrami operator)'
            call saddle_free(problem); return
        end if
        full = system%lift
        full(problem%free) = full(problem%free) + x(:problem%nf)
        if (present(dmu)) then
            ! Differentiated equation K dA/dmu = H_s A + P h (same operator).
            allocate (rhs(problem%n)); rhs = 0
            rhs(:problem%nf) = pack_free(problem, csc_matvec(system%helicity, full) &
                + system%helicity_h)
            residual = report%residual
            call krylov(problem, rhs, dx, info, report)
            report%residual = residual
            if (info /= MR_OK) then
                report%status = BELTRAMI_RESONANT
                report%message = 'MINRES failed on the mu tangent'
                call saddle_free(problem); return
            end if
            allocate (dmu(size(full))); dmu = 0
            dmu(problem%free) = dx(:problem%nf)
        end if
        call saddle_free(problem)
        report%helicity = 2*pi*(dot_product(full, csc_matvec(system%helicity, full)) &
            + 2*dot_product(full, system%helicity_h))
        allocate (sx(size(full)))
        call meridian_fourier_curl_curl_apply(system%op, system%curl, system%curl_t, &
            full, sx)
        report%energy = pi*(dot_product(full, sx) + 2*dot_product(full, &
            system%energy_h) + system%h_energy)/field%mu0
        report%current = mu*system%flux/field%mu0
        if (present(source)) report%current = report%current &
            + sum(source(system%n1 - size(system%gauge_index) + 1: &
            system%n1))/field%mu0
        bh = csc_matvec(system%curl, full)
        field%b(:, 0) = cmplx(bh(:system%n2), 0.0_dp, dp)
        do j = 1, system%modes
            field%b(:, j) = cmplx(bh((2*j - 1)*system%n2 + 1:2*j*system%n2), &
                -bh(2*j*system%n2 + 1:(2*j + 1)*system%n2), dp)/2
        end do
        field%psi = full(field%fs%nn + 1:system%n1)
        field%f_vac = system%harmonic
        call nn_divergence(field, report%divergence, st)
        if (st /= 0) return
        if (present(potential)) potential = full
        report%status = BELTRAMI_OK
        report%message = 'Beltrami candidate (MINRES); bounds unavailable'
    end subroutine

    module subroutine beltrami_mr_helicity(system, target, bracket, tolerance, &
        field, report)
        ! Newton-MR on the stationary point of the Taylor Lagrangian from
        ! mu0 = mid-bracket and A0 = the fixed-mu0 solution. The converged
        ! mu must lie in the bracket; the state is then re-solved at that mu
        ! (published only if |H - H0| <= tolerance max(1, |H0|)).
        type(beltrami_system_t), intent(in), target :: system
        real(dp), intent(in) :: target, bracket(2), tolerance
        type(free_b_field_t), intent(inout) :: field
        type(beltrami_report_t), intent(out) :: report
        type(saddle_t) :: problem
        type(free_b_field_t) :: candidate
        type(beltrami_report_t) :: work
        type(mr_report_t) :: newton
        type(mr_options_t) :: options
        real(dp), allocatable :: x(:), y(:)
        real(dp) :: mu0, t0
        integer :: info, st
        report%solver = 'newton-mr'; report%nfp = system%nfp
        mu0 = 0.5_dp*(bracket(1) + bracket(2))
        call saddle_init(problem, system, mu0, .false., work, st)
        if (st /= 0) then
            call saddle_free(problem); report = work; return
        end if
        call fixed_solve(problem, mu0, x, info, work)
        if (info /= MR_OK) then
            report = work; report%status = BELTRAMI_RESONANT
            report%message = 'MINRES failed at the initial multiplier'
            call saddle_free(problem); return
        end if
        ! Border the system: x = [A, p, mu]; keep the factors of mu0.
        problem%bordered = .true.; problem%target = target
        problem%n = problem%n + 1
        allocate (y(problem%n)); y(:problem%n - 1) = x; y(problem%n) = mu0
        options = system%mr%newton
        options%absolute_tolerance = max(options%absolute_tolerance, &
            0.1_dp*tolerance*max(1.0_dp, abs(target))/(4*pi))
        t0 = wall()
        call newton_mr_solve(problem, y, options, newton)
        report%minres_s = work%minres_s + wall() - t0
        report%minres_iterations = work%minres_iterations + newton%minres_iterations
        report%newton_iterations = newton%iterations
        report%preconditioner_s = work%preconditioner_s
        call saddle_free(problem)
        report%mu = y(problem%n)
        if (newton%status /= MR_OK) then
            report%status = BELTRAMI_INVALID
            report%message = 'helicity Newton-MR did not converge'
            return
        end if
        if (report%mu < bracket(1) .or. report%mu > bracket(2)) then
            report%status = BELTRAMI_INVALID
            report%message = 'helicity Newton-MR left the multiplier bracket'
            return
        end if
        candidate = field
        call beltrami_mr_solve(system, report%mu, candidate, work)
        work%minres_iterations = work%minres_iterations + report%minres_iterations
        work%newton_iterations = report%newton_iterations
        work%minres_s = work%minres_s + report%minres_s
        work%preconditioner_s = work%preconditioner_s + report%preconditioner_s
        work%solver = report%solver; work%mu = report%mu; work%nfp = system%nfp
        report = work
        if (report%status /= BELTRAMI_OK) return
        if (abs(report%helicity - target) > tolerance*max(1.0_dp, abs(target))) then
            report%status = BELTRAMI_INVALID
            report%message = 'helicity constraint did not converge'
            return
        end if
        field = candidate
    end subroutine

    subroutine saddle_init(problem, system, mu, bordered, report, status, source)
        ! Index lists, loads and the preconditioner for the multiplier mu.
        type(saddle_t), intent(out) :: problem
        type(beltrami_system_t), intent(in), target :: system
        real(dp), intent(in) :: mu
        logical, intent(in) :: bordered
        type(beltrami_report_t), intent(inout) :: report
        integer, intent(out) :: status
        real(dp), optional, intent(in) :: source(:)
        integer :: j
        real(dp) :: t0
        status = 1
        problem%sys => system
        problem%mu = mu; problem%bordered = bordered
        problem%kind = system%mr%preconditioner
        problem%nf = system%nf; problem%ng = system%ng
        problem%n = system%nf + system%ng + merge(1, 0, bordered)
        problem%free = pack([(j, j=1, size(system%index))], system%index > 0)
        problem%nodes = pack([(j, j=1, size(system%gauge_index))], &
            system%gauge_index > 0)
        problem%load0 = -system%energy_h
        if (present(source)) problem%load0 = problem%load0 + source
        problem%load1 = system%helicity_h
        problem%tau = system%mr%shift
        if (.not. problem%tau > 0) problem%tau = 0.01_dp/system%length2
        t0 = wall()
        call saddle_factor(problem, status)
        report%preconditioner_s = report%preconditioner_s + wall() - t0
    end subroutine

    subroutine saddle_factor(problem, status)
        ! Per-mode factors of Sbar_n + tau M_n (or their diagonals) and the
        ! factor of the gauge Laplacian L0.
        type(saddle_t), intent(inout) :: problem
        integer, intent(out) :: status
        type(csc_t) :: k
        logical, allocatable :: fixed(:)
        integer :: n, lo, hi, st
        status = 0
        if (problem%kind == BELTRAMI_PC_NONE) return
        status = 1
        if (problem%kind /= BELTRAMI_PC_JACOBI .and. &
            problem%kind /= BELTRAMI_PC_BLOCK) return
        associate (sys => problem%sys)
            allocate (problem%block(0:sys%modes), problem%active(0:sys%modes))
            problem%active = .false.
            if (problem%kind == BELTRAMI_PC_JACOBI) then
                allocate (problem%diag(problem%n)); problem%diag = 1
            end if
            do n = 0, sys%modes
                call mode_range(sys, n, lo, hi)
                fixed = sys%index(lo:hi) == 0
                if (all(fixed)) cycle
                call csc_add(sys%lumped(n), sys%mass(n), problem%tau, k, st)
                if (st /= 0) return
                problem%active(n) = .true.
                if (problem%kind == BELTRAMI_PC_JACOBI) then
                    problem%diag(pack(sys%index(lo:hi), .not. fixed)) = &
                        pack(diagonal(k), .not. fixed)
                    cycle
                end if
                call nn_factor(k, fixed, problem%block(n), st)
                if (st /= 0) return
            end do
            if (problem%ng > 0) then
                fixed = sys%gauge_index == 0
                if (problem%kind == BELTRAMI_PC_JACOBI) then
                    problem%diag(problem%nf + 1:problem%nf + problem%ng) = &
                        pack(diagonal(sys%laplace), .not. fixed)/problem%tau
                else
                    call nn_factor(sys%laplace, fixed, problem%laplace, st)
                    if (st /= 0) return
                end if
            end if
        end associate
        status = 0
    end subroutine

    subroutine fixed_solve(problem, mu, x, info, report)
        ! MINRES for the lifted fixed-mu system; report%residual is
        ! |K x - b| / max(1, |b|), the normalisation of the factored path.
        type(saddle_t), intent(inout) :: problem
        real(dp), intent(in) :: mu
        real(dp), allocatable, intent(out) :: x(:)
        integer, intent(out) :: info
        type(beltrami_report_t), intent(inout) :: report
        real(dp), allocatable :: rhs(:), ya(:), yp(:), zero(:)
        allocate (rhs(problem%n), zero(size(problem%sys%gauge_index)))
        zero = 0
        call saddle_core(problem%sys, problem%sys%lift, zero, mu, ya, yp, .true.)
        rhs(:problem%nf) = pack_free(problem, problem%load0 + mu*problem%load1 - ya)
        rhs(problem%nf + 1:) = -yp(problem%nodes)
        call krylov(problem, rhs, x, info, report)
    end subroutine

    subroutine krylov(problem, rhs, x, info, report)
        type(saddle_t), intent(inout) :: problem
        real(dp), intent(in) :: rhs(:)
        real(dp), allocatable, intent(out) :: x(:)
        integer, intent(out) :: info
        type(beltrami_report_t), intent(inout) :: report
        real(dp) :: bnorm, rtol, residual, t0
        integer :: its
        allocate (x(size(rhs)))
        bnorm = norm2(rhs)
        rtol = problem%sys%mr%tolerance*max(1.0_dp, bnorm)/max(bnorm, tiny(1.0_dp))
        t0 = wall()
        call minres_solve(problem, rhs, x, rtol, problem%sys%mr%max_minres, info, &
            its, residual, npc_exit=.false.)
        report%minres_s = report%minres_s + wall() - t0
        report%minres_iterations = report%minres_iterations + its
        report%residual = residual*bnorm/max(1.0_dp, bnorm)
    end subroutine

    subroutine saddle_core(sys, a, p, mu, ya, yp, lifted)
        ! Full-layout action: ya = (S - mu H_s) a + G p, yp = G^T a; with
        ! lifted, P C replaces H_s (the equation rows applied to the trace).
        type(beltrami_system_t), intent(in) :: sys
        real(dp), intent(in) :: a(:), p(:), mu
        real(dp), allocatable, intent(out) :: ya(:), yp(:)
        logical, intent(in) :: lifted
        allocate (ya(size(a)))
        call meridian_fourier_curl_curl_apply(sys%op, sys%curl, sys%curl_t, a, ya)
        if (lifted) then
            ya = ya - mu*csc_matvec(sys%helicity, a)
        else
            ya = ya - mu*symmetric_helicity(sys, a)
        end if
        ya(:sys%n1) = ya(:sys%n1) + csc_matvec(sys%gauge, p)
        yp = csc_matvec(sys%gauge_t, a(:sys%n1))
    end subroutine

    function symmetric_helicity(sys, a) result(y)
        type(beltrami_system_t), intent(in) :: sys
        real(dp), intent(in) :: a(:)
        real(dp), allocatable :: y(:)
        y = 0.5_dp*(csc_matvec(sys%helicity, a) + csc_matvec(sys%helicity_t, a))
    end function

    function pack_free(problem, full) result(x)
        type(saddle_t), intent(in) :: problem
        real(dp), intent(in) :: full(:)
        real(dp), allocatable :: x(:)
        x = full(problem%free)
    end function

    subroutine unpack(problem, x, a, p)
        ! Free-dof vector to zero-extended full A and V0 multiplier vectors.
        type(saddle_t), intent(in) :: problem
        real(dp), intent(in) :: x(:)
        real(dp), allocatable, intent(out) :: a(:), p(:)
        allocate (a(size(problem%sys%index)), p(size(problem%sys%gauge_index)))
        a = 0; p = 0
        a(problem%free) = x(:problem%nf)
        p(problem%nodes) = x(problem%nf + 1:problem%nf + problem%ng)
    end subroutine

    subroutine saddle_apply(self, x, y)
        class(saddle_t), intent(inout) :: self
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: y(:)
        real(dp), allocatable :: a(:), p(:), ya(:), yp(:)
        call unpack(self, x, a, p)
        call saddle_core(self%sys, a, p, self%mu, ya, yp, .false.)
        y(:self%nf) = ya(self%free)
        y(self%nf + 1:self%nf + self%ng) = yp(self%nodes)
        if (self%bordered) then
            y(:self%nf) = y(:self%nf) - self%border*x(self%n)
            y(self%n) = -dot_product(self%border, x(:self%nf))
        end if
    end subroutine

    subroutine saddle_evaluate(self, x, f, g, status)
        ! Lagrangian and its gradient at x = [A, p, mu] (bordered) or the
        ! fixed-mu equation residual at x = [A, p].
        class(saddle_t), intent(inout) :: self
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: f, g(:)
        integer, intent(out) :: status
        real(dp), allocatable :: a(:), p(:), ya(:), yp(:), yl(:), ypl(:)
        real(dp) :: mu, hq, ah
        mu = self%mu
        if (self%bordered) mu = x(self%n)
        call unpack(self, x, a, p)
        call saddle_core(self%sys, a, p, mu, ya, yp, .false.)
        if (any(abs(self%sys%lift) > 0)) then
            ! Open wall only (the bordered problem requires a closed wall).
            call saddle_core(self%sys, self%sys%lift, 0*p, mu, yl, ypl, .true.)
            ya = ya + yl; yp = yp + ypl
            a = a + self%sys%lift
        end if
        g(:self%nf) = pack_free(self, ya - self%load0 - mu*self%load1)
        g(self%nf + 1:self%nf + self%ng) = yp(self%nodes)
        ah = dot_product(a, self%load1)
        hq = dot_product(a, symmetric_helicity(self%sys, a)) + 2*ah
        ! 0.5 a.S a = 0.5 (a.ya + mu a.H_s a - a.G p); a.G p = p.yp.
        f = 0.5_dp*(dot_product(a, ya) + mu*(hq - 2*ah) - dot_product(p, yp)) &
            - dot_product(a, self%load0) + dot_product(p, yp)
        if (self%bordered) then
            f = f - mu*(2*pi*hq - self%target)/(4*pi)
            g(self%n) = -(2*pi*hq - self%target)/(4*pi)
        else
            f = f - mu*0.5_dp*hq
        end if
        status = merge(0, 1, ieee_is_finite(f))
    end subroutine

    subroutine saddle_linearize(self, x, status)
        ! Bordered Hessian at x: multiplier and border g = H_s A + P h; the
        ! mu block of the preconditioner is 1 / (g^T P_A g).
        class(saddle_t), intent(inout) :: self
        real(dp), intent(in) :: x(:)
        integer, intent(out) :: status
        real(dp), allocatable :: a(:), p(:), r(:), z(:)
        status = 0
        if (.not. self%bordered) return
        self%mu = x(self%n)
        if (self%kind == BELTRAMI_PC_JACOBI .and. size(self%diag) < self%n) &
            self%diag = [self%diag, 1.0_dp]
        call unpack(self, x, a, p)
        a = a + self%sys%lift
        self%border = pack_free(self, symmetric_helicity(self%sys, a) + self%load1)
        allocate (r(self%n), z(self%n)); r = 0
        r(:self%nf) = self%border
        self%schur = 1
        call self%precondition(r, z, status)
        if (status /= 0) return
        self%schur = dot_product(r(:self%nf), z(:self%nf))
        if (.not. self%schur > 0) self%schur = max(1.0_dp, norm2(self%border)**2)
        if (self%kind == BELTRAMI_PC_JACOBI) self%diag(self%n) = self%schur
    end subroutine

    subroutine saddle_precondition(self, r, z, status)
        class(saddle_t), intent(inout) :: self
        real(dp), intent(in) :: r(:)
        real(dp), intent(out) :: z(:)
        integer, intent(out) :: status
        real(dp), allocatable :: a(:), p(:), za(:), zp(:)
        integer :: n, lo, hi
        status = 0
        select case (self%kind)
        case (BELTRAMI_PC_NONE)
            z = r; return
        case (BELTRAMI_PC_JACOBI)
            z = r/self%diag; return
        end select
        call unpack(self, r, a, p)
        allocate (za(size(a)), zp(size(p))); za = 0; zp = 0
        do n = 0, self%sys%modes
            if (.not. self%active(n)) cycle
            call mode_range(self%sys, n, lo, hi)
            call nn_factor_solve_real(self%block(n), a(lo:hi), za(lo:hi), status)
            if (status /= 0) return
        end do
        if (self%ng > 0) then
            call nn_factor_solve_real(self%laplace, p, zp, status)
            if (status /= 0) return
        end if
        z(:self%nf) = za(self%free)
        z(self%nf + 1:self%nf + self%ng) = self%tau*zp(self%nodes)
        if (self%bordered) z(self%n) = r(self%n)/self%schur
    end subroutine

    subroutine saddle_free(problem)
        type(saddle_t), intent(inout) :: problem
        integer :: n
        if (allocated(problem%block)) then
            do n = lbound(problem%block, 1), ubound(problem%block, 1)
                call nn_factor_free(problem%block(n))
            end do
            deallocate (problem%block)
        end if
        call nn_factor_free(problem%laplace)
    end subroutine

    subroutine mode_range(sys, n, lo, hi)
        ! Full-layout V1 range of mode n: mean, or the (cos, sin) pair.
        type(beltrami_system_t), intent(in) :: sys
        integer, intent(in) :: n
        integer, intent(out) :: lo, hi
        lo = merge(1, (2*n - 1)*sys%n1 + 1, n == 0)
        hi = merge(sys%n1, (2*n + 1)*sys%n1, n == 0)
    end subroutine

    subroutine csc_add(a, b, beta, c, status)
        ! c = a + beta b (same shape).
        type(csc_t), intent(in) :: a, b
        real(dp), intent(in) :: beta
        type(csc_t), intent(out) :: c
        integer, intent(out) :: status
        type(fortsparse_status_t) :: fst
        integer, allocatable :: row(:), col(:)
        real(dp), allocatable :: val(:)
        integer :: j, q, k
        allocate (row(a%nnz + b%nnz), col(a%nnz + b%nnz), val(a%nnz + b%nnz))
        k = 0
        do j = 1, a%ncol
            do q = a%col_ptr(j), a%col_ptr(j + 1) - 1
                k = k + 1; row(k) = a%row_idx(q); col(k) = j; val(k) = a%val(q)
            end do
            do q = b%col_ptr(j), b%col_ptr(j + 1) - 1
                k = k + 1; row(k) = b%row_idx(q); col(k) = j; val(k) = beta*b%val(q)
            end do
        end do
        call csc_from_triplet(a%nrow, a%ncol, row(:k), col(:k), val(:k), c, fst)
        status = merge(0, 1, status_ok(fst))
    end subroutine

    function diagonal(a) result(d)
        type(csc_t), intent(in) :: a
        real(dp), allocatable :: d(:)
        integer :: j, q
        allocate (d(a%ncol)); d = 0
        do j = 1, a%ncol
            do q = a%col_ptr(j), a%col_ptr(j + 1) - 1
                if (a%row_idx(q) == j) d(j) = d(j) + a%val(q)
            end do
        end do
    end function

    real(dp) function wall()
        integer(i8) :: count, rate
        call system_clock(count, rate)
        wall = real(count, dp)/real(rate, dp)
    end function
end submodule kin6d_nonnested_beltrami_mr
