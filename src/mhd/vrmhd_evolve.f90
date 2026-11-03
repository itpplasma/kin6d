module kin6d_vrmhd_evolve
    !! Weyl-gauge time evolution of the invariant MHD section as a FortNum
    !! BDF problem; see docs/vrmhd-section.md#time-evolution.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_is_finite
    use kin6d_vrmhd_parameters, only: vrmhd_parameters_t
    use kin6d_vrmhd_section, only: vrmhd_section_t, vrmhd_section_residual, &
        vrmhd_section_jacobian
    use fortnum_ode_bdf, only: bdf_problem_t
    use fortsparse, only: csc_t, csc_matvec, csc_transpose, csc_matmul, &
        csc_from_triplet, sparse_solver_t, sparse_factor, sparse_solve, &
        sparse_free, sparse_solve_once, fortsparse_status_t, status_ok
    implicit none
    private

    !> The state is u = lift + P y. Rows are tested with P^T, so a periodic
    !> identification or a Dirichlet elimination is a choice of P and lift.
    type, extends(bdf_problem_t), public :: vrmhd_evolve_t
        type(vrmhd_section_t), pointer :: op => null()
        type(vrmhd_parameters_t) :: p
        real(dp), allocatable :: source(:, :, :), load(:), lift(:)
        type(csc_t) :: reduction, restriction
        logical, allocatable :: algebraic(:)
        integer :: n = 0
        !> Linearizations that reuse one factorization; 1 is exact Newton.
        integer :: lag = 1
        !> Relative change of the BDF shift alpha that forces a new factor.
        real(dp) :: shift_tolerance = 0.05_dp
        integer :: since = 0, factorizations = 0
        real(dp), allocatable :: u(:), ud(:)
        real(dp) :: alpha = 0, factored_alpha = 0
        logical :: factored = .false.
        type(sparse_solver_t) :: solver
    contains
        procedure :: evaluate => evolve_evaluate
        procedure :: linearize => evolve_linearize
        procedure :: apply => evolve_apply
        procedure :: precondition => evolve_precondition
    end type
    public :: vrmhd_evolve_build, vrmhd_evolve_expand, vrmhd_evolve_consistent
contains
    subroutine vrmhd_evolve_build(problem, op, p, source, lift, reduction, status, &
            load, lag)
        type(vrmhd_evolve_t), intent(out) :: problem
        type(vrmhd_section_t), intent(in), target :: op
        type(vrmhd_parameters_t), intent(in) :: p
        real(dp), intent(in) :: source(:, :, :), lift(:)
        type(csc_t), intent(in) :: reduction
        integer, intent(out) :: status
        !> Voltage load of vrmhd_drive_build; the residual subtracts it.
        real(dp), intent(in), optional :: load(:)
        integer, intent(in), optional :: lag
        type(fortsparse_status_t) :: st
        integer :: j, k, current
        logical :: any_current, any_other
        status = 1
        if (reduction%nrow /= op%ndof .or. reduction%ncol < 1) return
        if (size(lift) /= op%ndof) return
        if (size(source, 1) /= 6 .or. size(source, 2) /= op%fs%nq .or. &
            size(source, 3) /= op%fs%mm%topology%n_triangles) return
        problem%op => op
        problem%p = p
        problem%source = source
        problem%lift = lift
        allocate (problem%load(op%ndof), source=0.0_dp)
        if (present(load)) then
            if (size(load) /= op%ndof) return
            problem%load = load
        end if
        if (present(lag)) then
            if (lag < 1) return
            problem%lag = lag
        end if
        problem%reduction = reduction
        call csc_transpose(reduction, problem%restriction, st)
        if (.not. status_ok(st)) return
        problem%n = reduction%ncol
        allocate (problem%algebraic(problem%n), problem%u(op%ndof), &
            problem%ud(op%ndof))
        ! Weak Ampere rows have zero mass; a reduced coordinate must not mix them.
        do j = 1, problem%n
            any_current = .false.; any_other = .false.
            do k = reduction%col_ptr(j), reduction%col_ptr(j + 1) - 1
                current = current_row(op, reduction%row_idx(k))
                any_current = any_current .or. current == 1
                any_other = any_other .or. current == 0
            end do
            if (any_current .eqv. any_other) return
            problem%algebraic(j) = any_current
        end do
        status = 0
    end subroutine

    pure integer function current_row(op, i) result(current)
        type(vrmhd_section_t), intent(in) :: op
        integer, intent(in) :: i
        integer :: stride, local
        stride = op%offset(9)
        local = modulo(i - 1, stride) + 1
        current = 0
        if (local > op%offset(8)) current = 1
    end function

    subroutine vrmhd_evolve_expand(problem, y, u)
        type(vrmhd_evolve_t), intent(in) :: problem
        real(dp), intent(in) :: y(:)
        real(dp), intent(out) :: u(:)
        u = problem%lift + csc_matvec(problem%reduction, y)
    end subroutine

    subroutine evolve_evaluate(self, t, y, ydot, residual, status)
        class(vrmhd_evolve_t), intent(inout) :: self
        real(dp), intent(in) :: t, y(:), ydot(:)
        real(dp), intent(out) :: residual(:)
        integer, intent(out) :: status
        real(dp), allocatable :: u(:), ud(:), zero(:), r(:), jv(:)
        status = 1
        residual = 0
        if (.not. ieee_is_finite(t)) return
        allocate (u(self%op%ndof), ud(self%op%ndof), zero(self%op%ndof), &
            r(self%op%ndof), jv(self%op%ndof))
        call vrmhd_evolve_expand(self, y, u)
        ud = csc_matvec(self%reduction, ydot)
        zero = 0
        call vrmhd_section_residual(self%op, self%p, u, ud, zero, zero, &
            self%source, r, jv, status)
        residual = 0
        if (status /= 0) return
        residual = csc_matvec(self%restriction, r - self%load)
    end subroutine

    subroutine evolve_linearize(self, t, y, ydot, alpha, status)
        class(vrmhd_evolve_t), intent(inout) :: self
        real(dp), intent(in) :: t, y(:), ydot(:), alpha
        integer, intent(out) :: status
        type(csc_t) :: reduced
        type(fortsparse_status_t) :: st
        status = 1
        if (.not. ieee_is_finite(t)) return
        call vrmhd_evolve_expand(self, y, self%u)
        self%ud = csc_matvec(self%reduction, ydot)
        self%alpha = alpha
        status = 0
        ! A lagged factor stays a good preconditioner while the shift is close.
        if (self%factored .and. self%since < self%lag .and. &
            abs(alpha - self%factored_alpha) <= &
            self%shift_tolerance*abs(self%factored_alpha)) then
            self%since = self%since + 1
            return
        end if
        status = 1
        call reduced_jacobian(self, self%u, self%ud, alpha, reduced, status)
        if (status /= 0) return
        if (self%factored) call sparse_free(self%solver)
        self%factored = .false.
        call sparse_factor(self%solver, reduced, st)
        status = 1
        if (.not. status_ok(st)) return
        self%factored = .true.
        self%factored_alpha = alpha
        self%factorizations = self%factorizations + 1
        self%since = 1
        status = 0
    end subroutine

    subroutine evolve_apply(self, input, output)
        class(vrmhd_evolve_t), intent(inout) :: self
        real(dp), intent(in) :: input(:)
        real(dp), intent(out) :: output(:)
        real(dp), allocatable :: du(:), dud(:), r(:), jv(:)
        integer :: status
        allocate (du(self%op%ndof), dud(self%op%ndof), r(self%op%ndof), &
            jv(self%op%ndof))
        du = csc_matvec(self%reduction, input)
        dud = self%alpha*du
        call vrmhd_section_residual(self%op, self%p, self%u, self%ud, du, dud, &
            self%source, r, jv, status)
        if (status /= 0) then
            output = ieee_value(1.0_dp, ieee_quiet_nan)
            return
        end if
        output = csc_matvec(self%restriction, jv)
    end subroutine

    subroutine evolve_precondition(self, r, z, status)
        class(vrmhd_evolve_t), intent(inout) :: self
        real(dp), intent(in) :: r(:)
        real(dp), intent(out) :: z(:)
        integer, intent(out) :: status
        real(dp), allocatable :: b(:), x(:)
        type(fortsparse_status_t) :: st
        status = 1
        z = 0
        if (.not. self%factored) return
        b = r
        allocate (x(size(r)))
        call sparse_solve(self%solver, b, x, st)
        if (.not. status_ok(st)) return
        z = x
        status = 0
    end subroutine

    subroutine reduced_jacobian(self, u, ud, shift, reduced, status)
        !! P^T (R_U + (M ud)_U + shift M) P from the assembled section Jacobian.
        type(vrmhd_evolve_t), intent(in) :: self
        real(dp), intent(in) :: u(:), ud(:), shift
        type(csc_t), intent(out) :: reduced
        integer, intent(out) :: status
        type(csc_t) :: jac, half
        type(fortsparse_status_t) :: st
        call vrmhd_section_jacobian(self%op, self%p, u, ud, self%source, shift, &
            jac, status)
        if (status /= 0) return
        status = 1
        call csc_matmul(jac, self%reduction, half, st)
        if (.not. status_ok(st)) return
        call csc_matmul(self%restriction, half, reduced, st)
        if (status_ok(st)) status = 0
    end subroutine

    subroutine vrmhd_evolve_consistent(problem, t, y, ydot, status)
        !! Solves weak Ampere for the current at the given vector potential and
        !! returns ydot from M ydot = -R on differential rows and the time
        !! derivative of the algebraic rows, d/dt R_J = R_J,U ydot = 0.
        type(vrmhd_evolve_t), intent(inout) :: problem
        real(dp), intent(in) :: t
        real(dp), intent(inout) :: y(:)
        real(dp), intent(out) :: ydot(:)
        integer, intent(out) :: status
        type(csc_t) :: steady, shifted, ampere, combined
        type(fortsparse_status_t) :: st
        real(dp), allocatable :: zero(:), r(:), rhs(:), step(:), u(:), uzero(:)
        integer, allocatable :: index(:), rows(:), cols(:)
        real(dp), allocatable :: vals(:)
        integer :: i, j, k, na, pos
        status = 1
        ydot = 0
        if (size(y) /= problem%n .or. size(ydot) /= problem%n) return
        allocate (zero(problem%n), r(problem%n), u(problem%op%ndof), &
            uzero(problem%op%ndof), index(problem%n))
        zero = 0; uzero = 0
        na = count(problem%algebraic)
        index = 0
        index(pack([(i, i=1, problem%n)], problem%algebraic)) = [(i, i=1, na)]
        status = 0
        if (na > 0) call solve_current(status)
        if (status /= 0) return
        call problem%evaluate(t, y, zero, r, status)
        if (status /= 0) return
        call vrmhd_evolve_expand(problem, y, u)
        call reduced_jacobian(problem, u, uzero, 0.0_dp, steady, status)
        if (status /= 0) return
        call reduced_jacobian(problem, u, uzero, 1.0_dp, shifted, status)
        if (status /= 0) return
        status = 1
        ! Differential rows take the mass (shift 1 minus shift 0); algebraic
        ! rows keep the steady Jacobian.
        allocate (rows(size(steady%val) + size(shifted%val)), &
            cols(size(steady%val) + size(shifted%val)), &
            vals(size(steady%val) + size(shifted%val)))
        pos = 0
        do j = 1, problem%n
            do k = shifted%col_ptr(j), shifted%col_ptr(j + 1) - 1
                i = shifted%row_idx(k)
                if (problem%algebraic(i)) cycle
                pos = pos + 1
                rows(pos) = i; cols(pos) = j; vals(pos) = shifted%val(k)
            end do
            do k = steady%col_ptr(j), steady%col_ptr(j + 1) - 1
                i = steady%row_idx(k)
                pos = pos + 1
                rows(pos) = i; cols(pos) = j
                vals(pos) = steady%val(k)
                if (.not. problem%algebraic(i)) vals(pos) = -steady%val(k)
            end do
        end do
        call csc_from_triplet(problem%n, problem%n, rows(:pos), cols(:pos), &
            vals(:pos), combined, st)
        if (.not. status_ok(st)) return
        rhs = -r
        where (problem%algebraic) rhs = 0
        call sparse_solve_once(combined, rhs, ydot, st)
        if (status_ok(st)) status = 0
    contains
        subroutine solve_current(status)
            !! Ampere is linear in the current, so one direct solve suffices.
            integer, intent(out) :: status
            call problem%evaluate(t, y, zero, r, status)
            if (status /= 0) return
            call vrmhd_evolve_expand(problem, y, u)
            call reduced_jacobian(problem, u, uzero, 0.0_dp, steady, status)
            if (status /= 0) return
            status = 1
            allocate (rows(size(steady%val)), cols(size(steady%val)), &
                vals(size(steady%val)))
            pos = 0
            do j = 1, problem%n
                if (index(j) == 0) cycle
                do k = steady%col_ptr(j), steady%col_ptr(j + 1) - 1
                    i = steady%row_idx(k)
                    if (index(i) == 0) cycle
                    pos = pos + 1
                    rows(pos) = index(i); cols(pos) = index(j); vals(pos) = steady%val(k)
                end do
            end do
            call csc_from_triplet(na, na, rows(:pos), cols(:pos), vals(:pos), ampere, st)
            deallocate (rows, cols, vals)
            if (.not. status_ok(st)) return
            rhs = -pack(r, problem%algebraic)
            allocate (step(na))
            call sparse_solve_once(ampere, rhs, step, st)
            if (.not. status_ok(st)) return
            y = unpack(pack(y, problem%algebraic) + step, problem%algebraic, y)
            status = 0
        end subroutine
    end subroutine
end module kin6d_vrmhd_evolve
