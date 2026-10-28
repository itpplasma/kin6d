module newton_mr_dense_fixture
    ! Dense operators and objectives for the kin6d_newton_mr unit checks.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_newton_mr, only: mr_operator_t, mr_problem_t
    implicit none

    type, extends(mr_operator_t) :: dense_op_t
        real(dp), allocatable :: a(:, :), m(:) ! m: diagonal preconditioner
    contains
        procedure :: apply => dense_apply
        procedure :: precondition => dense_precondition
    end type

    type, extends(mr_problem_t) :: quartic_t
        ! f = x^T D x / 2 + sum x^4 / 4 - c^T x, D symmetric indefinite.
        real(dp), allocatable :: d(:, :), c(:), x0(:)
    contains
        procedure :: apply => quartic_apply
        procedure :: evaluate => quartic_evaluate
        procedure :: linearize => quartic_linearize
    end type
contains
    subroutine dense_apply(self, x, y)
        class(dense_op_t), intent(inout) :: self
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: y(:)
        y = matmul(self%a, x)
    end subroutine

    subroutine dense_precondition(self, r, z, status)
        class(dense_op_t), intent(inout) :: self
        real(dp), intent(in) :: r(:)
        real(dp), intent(out) :: z(:)
        integer, intent(out) :: status
        z = r
        if (allocated(self%m)) z = self%m*r
        status = 0
    end subroutine

    subroutine quartic_apply(self, x, y)
        class(quartic_t), intent(inout) :: self
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: y(:)
        y = matmul(self%d, x) + 3*self%x0**2*x
    end subroutine

    subroutine quartic_linearize(self, x, status)
        class(quartic_t), intent(inout) :: self
        real(dp), intent(in) :: x(:)
        integer, intent(out) :: status
        self%x0 = x
        status = 0
    end subroutine

    subroutine quartic_evaluate(self, x, f, g, status)
        class(quartic_t), intent(inout) :: self
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: f, g(:)
        integer, intent(out) :: status
        g = matmul(self%d, x)
        f = 0.5_dp*dot_product(x, g) + 0.25_dp*sum(x**4) - dot_product(self%c, x)
        g = g + x**3 - self%c
        status = 0
    end subroutine

    function dense_solve(a, b) result(x)
        ! Gaussian elimination with partial pivoting (test oracle).
        real(dp), intent(in) :: a(:, :), b(:)
        real(dp) :: x(size(b)), m(size(b), size(b) + 1), row(size(b) + 1)
        integer :: n, k, p, i
        n = size(b)
        m(:, :n) = a; m(:, n + 1) = b
        do k = 1, n
            p = k - 1 + maxloc(abs(m(k:, k)), 1)
            row = m(k, :); m(k, :) = m(p, :); m(p, :) = row
            do i = k + 1, n
                m(i, k:) = m(i, k:) - m(i, k)/m(k, k)*m(k, k:)
            end do
        end do
        do k = n, 1, -1
            x(k) = (m(k, n + 1) - dot_product(m(k, k + 1:n), x(k + 1:n)))/m(k, k)
        end do
    end function

    function orthogonal(n, seed) result(q)
        ! Deterministic orthogonal matrix by modified Gram-Schmidt.
        integer, intent(in) :: n, seed
        real(dp) :: q(n, n)
        integer :: i, j
        do j = 1, n
            do i = 1, n
                q(i, j) = sin(real(seed*i + 7*j*j + i*j, dp))
            end do
        end do
        do j = 1, n
            do i = 1, j - 1
                q(:, j) = q(:, j) - dot_product(q(:, i), q(:, j))*q(:, i)
            end do
            q(:, j) = q(:, j)/norm2(q(:, j))
        end do
    end function

end module newton_mr_dense_fixture

program test_newton_mr
    ! Usage: minres | newton
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use newton_mr_dense_fixture
    use kin6d_newton_mr, only: minres_solve, newton_mr_solve, mr_options_t, &
        mr_report_t, MR_OK, MR_NPC
    implicit none
    character(32) :: mode
    call get_command_argument(1, mode)
    select case (trim(mode))
    case ('minres')
        call minres_checks()
    case ('newton')
        call newton_checks()
    case default
        error stop 'usage: minres | newton'
    end select
contains
    subroutine minres_checks()
        ! MINRES against a dense solve: SPD (with and without a diagonal
        ! preconditioner), indefinite, singular consistent (minimum norm),
        ! and the nonpositive-curvature exit.
        integer, parameter :: n = 40
        type(dense_op_t) :: op
        real(dp) :: q(n, n), lam(n), b(n), x(n), xref(n), res, curv
        integer :: info, its, k
        q = orthogonal(n, 3)
        do k = 1, n
            lam(k) = 10.0_dp**(3*real(k - 1, dp)/(n - 1)) ! cond 1e3
            b(k) = cos(real(k, dp))
        end do
        op%a = matmul(q, matmul(diag(lam), transpose(q)))
        xref = dense_solve(op%a, b)
        call minres_solve(op, b, x, 1.0e-12_dp, 500, info, its, res)
        write (*, '(a,i4,2es11.3)') 'minres SPD: its, residual, error', its, res, &
            norm2(x - xref)/norm2(xref)
        if (info /= MR_OK .or. res > 2.0e-12_dp) error stop 'minres SPD residual'
        if (norm2(x - xref) > 1.0e-8_dp*norm2(xref)) error stop 'minres SPD error'
        ! Diagonal preconditioner on a badly scaled SPD matrix.
        op%a = matmul(diag(lam), matmul(op%a, diag(lam)))
        xref = dense_solve(op%a, b)
        allocate (op%m(n))
        do k = 1, n
            op%m(k) = 1/op%a(k, k)
        end do
        call minres_solve(op, b, x, 1.0e-12_dp, 500, info, its, res)
        write (*, '(a,i4,2es11.3)') 'minres PSPD: its, residual, error', its, res, &
            norm2(x - xref)/norm2(xref)
        ! cond(A) ~ 1e9: the recurred residual drifts from the true one by
        ! about cond eps, so the explicit residual gate is looser.
        if (info /= MR_OK .or. res > 1.0e-10_dp) error stop 'minres preconditioned'
        if (norm2(x - xref) > 1.0e-6_dp*norm2(xref)) error stop 'minres PSPD error'
        deallocate (op%m)
        ! Indefinite: half of the spectrum negative.
        lam(1:n:2) = -lam(1:n:2)
        op%a = matmul(q, matmul(diag(lam), transpose(q)))
        xref = dense_solve(op%a, b)
        call minres_solve(op, b, x, 1.0e-12_dp, 500, info, its, res)
        write (*, '(a,i4,2es11.3)') 'minres indefinite: its, residual, error', its, &
            res, norm2(x - xref)/norm2(xref)
        if (info /= MR_OK .or. res > 2.0e-12_dp) error stop 'minres indefinite'
        if (norm2(x - xref) > 1.0e-8_dp*norm2(xref)) error stop 'minres indef error'
        ! NPC exit: the returned direction has nonpositive curvature and is a
        ! descent direction of the Newton system A d = -g (b = -g).
        call minres_solve(op, b, x, 1.0e-12_dp, 500, info, its, res, npc_exit=.true.)
        curv = dot_product(x, matmul(op%a, x))
        write (*, '(a,i4,2es11.3)') 'minres NPC: its, curvature, g.d', its, curv, &
            -dot_product(b, x)
        if (info /= MR_NPC .or. curv > 0) error stop 'minres NPC exit'
        if (.not. dot_product(b, x) > 0) error stop 'minres NPC not descent'
        ! Singular PSD with consistent b: the minimum-norm solution.
        lam = abs(lam); lam(n - 4:) = 0
        op%a = matmul(q, matmul(diag(lam), transpose(q)))
        b = matmul(op%a, [(sin(real(k, dp)), k=1, n)])
        xref = matmul(q, matmul(diag(merge(1/max(lam, tiny(1.0_dp)), 0.0_dp, &
            lam > 0)), matmul(transpose(q), b)))
        call minres_solve(op, b, x, 1.0e-12_dp, 500, info, its, res, npc_exit=.true.)
        write (*, '(a,i4,2es11.3)') 'minres singular: its, residual, error', its, res, &
            norm2(x - xref)/norm2(xref)
        if (info /= MR_OK .or. res > 2.0e-12_dp) error stop 'minres singular'
        if (norm2(x - xref) > 1.0e-8_dp*norm2(xref)) error stop 'minres min norm'
    end subroutine

    function diag(v) result(d)
        real(dp), intent(in) :: v(:)
        real(dp) :: d(size(v), size(v))
        integer :: k
        d = 0
        do k = 1, size(v)
            d(k, k) = v(k)
        end do
    end function

    subroutine newton_checks()
        ! Newton-MR on a convex quartic (no NPC step) and from a saddle of a
        ! nonconvex one (NPC steps, then a local minimiser with a positive
        ! Hessian).
        integer, parameter :: n = 20
        type(quartic_t) :: p
        type(mr_options_t) :: opt
        type(mr_report_t) :: rep
        real(dp) :: q(n, n), lam(n), x(n), h(n, n), hl(n, n)
        integer :: k, j
        q = orthogonal(n, 5)
        do k = 1, n
            lam(k) = 1 + k
        end do
        p%d = matmul(q, matmul(diag(lam), transpose(q)))
        p%c = [(cos(real(k, dp)), k=1, n)]
        allocate (p%x0(n))
        ! Convex (positive D): no NPC step, gradient tolerance met.
        x = 0
        opt%tolerance = 1.0e-12_dp
        call newton_mr_solve(p, x, opt, rep)
        write (*, '(a,i3,i5,i3,es11.3)') 'newton convex: steps, minres, npc, |g|/|g0|', &
            rep%iterations, rep%minres_iterations, rep%npc_steps, &
            rep%gradient/rep%initial_gradient
        if (rep%status /= MR_OK .or. rep%npc_steps /= 0) error stop 'newton convex'
        ! Nonconvex: D with negative eigenvalues; start near the saddle 0.
        lam(1:n:3) = -2*lam(1:n:3)
        p%d = matmul(q, matmul(diag(lam), transpose(q)))
        p%c = 1.0e-3_dp*p%c
        x = 1.0e-4_dp*[(sin(real(k, dp)), k=1, n)]
        ! |g0| is small near the saddle; the roundoff floor of |g| at the
        ! minimiser (|x| ~ 5) is about 1e-13, i.e. 3e-11 relative.
        opt%tolerance = 1.0e-9_dp
        call newton_mr_solve(p, x, opt, rep)
        h = p%d
        do k = 1, n
            h(k, k) = h(k, k) + 3*x(k)**2
        end do
        ! Positive definiteness by an unpivoted Cholesky attempt.
        hl = h
        do j = 1, n
            if (.not. hl(j, j) - sum(hl(j, :j - 1)**2) > 0) &
                error stop 'newton nonconvex: Hessian not positive'
            hl(j, j) = sqrt(hl(j, j) - sum(hl(j, :j - 1)**2))
            do k = j + 1, n
                hl(k, j) = (hl(k, j) - dot_product(hl(k, :j - 1), hl(j, :j - 1)))/hl(j, j)
            end do
        end do
        write (*, '(a,i3,i5,i3,i3,2es11.3)') 'newton nonconvex: steps, minres, npc, ' &
            //'backtracks, |g|/|g0|, f', rep%iterations, rep%minres_iterations, &
            rep%npc_steps, rep%backtracks, rep%gradient/rep%initial_gradient, rep%energy
        if (rep%status /= MR_OK) error stop 'newton nonconvex status'
        if (rep%npc_steps < 1) error stop 'newton nonconvex: no NPC step taken'
        if (.not. rep%energy < 0) error stop 'newton nonconvex: saddle not left'
    end subroutine

end program test_newton_mr
