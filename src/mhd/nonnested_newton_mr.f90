module kin6d_nonnested_newton_mr
    ! Opt-in matrix-free Newton-MR for the mapped NN1 vacuum energy
    ! (docs/nn1-mapped-3d.md, docs/nonnested-newton-mr.md):
    !     W(A) = int |h + curl A|^2 / (2 mu0),   S = C^T M_2(K) C,
    ! with the boundary traces of data fixed. The default path,
    ! kin6d_nonnested_mapped:nn_mapped_solve, assembles S (mode-dense element
    ! blocks) and solves by FGMRES with a factored mean-metric preconditioner.
    ! This path never assembles S: the Hessian action is FortFEM's
    ! sum-factorised meridian_fourier_curl_curl_apply (FFT in phi per
    ! meridional quadrature point), the step is Newton-MR with MINRES
    ! (kin6d_newton_mr), and the preconditioner is chosen by options:
    !   NN_MR_PC_NONE   identity on the free dofs;
    !   NN_MR_PC_JACOBI diagonal of the metric-lumped operator below;
    !   NN_MR_PC_LUMPED metric lumping: K is replaced by its phi average
    !                   Kbar(R, Z) at every meridional quadrature point. The
    !                   lumped operator commutes with the Fourier basis, so it
    !                   is a Kronecker sum of the meridional forms with the
    !                   diagonal 1-D phi mass (1, 1/2, ...) and stiffness n^2;
    !                   in phi the fast diagonalisation is the Fourier basis
    !                   itself. Each mode n is one meridional system
    !                   C_n^T M_2,n(Kbar) C_n on the (cos, sin) pair, factored
    !                   once by nn_factor. The meridional mesh is unstructured
    !                   (triangles), so no further Kronecker split in (R, Z)
    !                   exists; the per-mode factor is the exact inverse of the
    !                   lumped model instead of a 1-D fast diagonalisation.
    !                   Triangles have no polar logical axis, so no dense axis
    !                   block is needed either.
    ! The n = 0 Nedelec block of S is singular on grad V0; the right-hand
    ! side is consistent (C grad = 0), MINRES stays in the range, and the
    ! lumped n = 0 factor is regularised by eps M_1 exactly as the default
    ! preconditioner. B = C A is unique. Geometric multigrid is not offered:
    ! the meridional ring meshes are not nested under refinement.
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_meridian_feec_n, only: meridian_feec_t
    use fortfem_meridian_fourier_scalar, only: meridian_fourier_tensor_t
    use fortfem_meridian_fourier_curl, only: meridian_fourier_mass, &
        meridian_fourier_curl_curl, meridian_fourier_curl
    use fortfem_meridian_fourier_apply, only: meridian_fourier_op_t, &
        meridian_fourier_op_init, meridian_fourier_curl_curl_apply, &
        meridian_fourier_op_bytes
    use fortsparse, only: csc_t, csc_matvec, csc_matmul, csc_transpose, &
        csc_from_triplet, fortsparse_status_t, status_ok
    use kin6d_newton_mr, only: mr_problem_t, mr_options_t, mr_report_t, &
        newton_mr_solve, MR_OK
    use kin6d_equilibrium_nonnested, only: free_b_field_t, nn_factor_t, nn_factor, &
        nn_factor_solve_real, nn_factor_free, nn_divergence
    use kin6d_nonnested_mapped, only: nn_map_t, map_metric_t, nn_harmonic_forms, &
        nn_mapped_report_t
    implicit none
    private
    public :: nn_mr_options_t, nn_mr_report_t, nn_mapped_solve_mr
    public :: lumped_metric_t, mapped_vacuum_t, nn_mr_prepare, nn_mr_free
    ! Per-mode lumped pieces, shared with the Beltrami MINRES path.
    public :: mode_curl, mode_mass
    integer, parameter, public :: NN_MR_PC_NONE = 0, NN_MR_PC_JACOBI = 1, &
        NN_MR_PC_LUMPED = 2
    real(dp), parameter :: pi = acos(-1.0_dp)

    type :: nn_mr_options_t
        integer :: preconditioner = NN_MR_PC_LUMPED
        ! False: apply the assembled S instead (comparison only).
        logical :: matrix_free = .true.
        ! Angles of the phi average of K in the lumped model.
        integer :: nphi_lump = 0 ! 0: the assembly rule nphi
        type(mr_options_t) :: newton
    end type

    type :: nn_mr_report_t
        type(mr_report_t) :: newton
        ! Wall seconds of the operator setup (FortFEM apply data or the
        ! assembled S) and of the preconditioner setup (lumped assembly and
        ! per-mode factors); bytes of the matrix-free operator data.
        real(dp) :: operator_s = 0, preconditioner_s = 0, total_s = 0
        integer(int64) :: operator_bytes = 0
    end type

    type, extends(meridian_fourier_tensor_t) :: lumped_metric_t
        ! Phi average of a mapped metric on nphi uniform angles.
        type(map_metric_t), pointer :: base => null()
        integer :: nphi = 0
    contains
        procedure :: tensor => lumped_tensor
    end type

    type, extends(mr_problem_t) :: mapped_vacuum_t
        ! f(A) = A^T S A / 2 - rhs^T A on the free dofs; fixed dofs keep data.
        integer :: modes = 0, nb = 1, n1 = 0, nn = 0, kind = NN_MR_PC_LUMPED
        logical :: matrix_free = .true., failed = .false.
        type(meridian_fourier_op_t) :: op
        type(csc_t) :: c, ct, s
        logical, allocatable :: fixed(:)
        real(dp), allocatable :: rhs(:), diag(:)
        type(nn_factor_t), allocatable :: block(:) ! modes 0..N
    contains
        procedure :: apply => vacuum_apply
        procedure :: precondition => vacuum_precondition
        procedure :: evaluate => vacuum_evaluate
    end type
contains
    subroutine nn_mapped_solve_mr(field, map, nphi, data, options, report, mr, status)
        ! Same problem, data, fixed traces and report as nn_mapped_solve.
        ! report%iterations counts MINRES iterations of all Newton steps;
        ! report%info is the Newton-MR status (MR_OK = 0).
        type(free_b_field_t), intent(inout), target :: field
        type(nn_map_t), intent(in), target :: map
        integer, intent(in) :: nphi
        real(dp), intent(in) :: data(:)
        type(nn_mr_options_t), intent(in) :: options
        type(nn_mapped_report_t), intent(out) :: report
        type(nn_mr_report_t), intent(out) :: mr
        integer, intent(out) :: status
        type(mapped_vacuum_t) :: problem
        type(map_metric_t), target :: metric
        real(dp), allocatable :: x(:), sx(:), r(:), bh(:)
        real(dp) :: energy_h, scale, t0
        integer :: n, n2, st
        status = 1
        t0 = clock()
        metric%fs => field%fs; metric%map => map
        call nn_mr_prepare(field, metric, nphi, options, problem, energy_h, mr, st)
        if (st /= 0) then
            call nn_mr_free(problem); return
        end if
        if (size(data) /= size(problem%rhs)) then
            call nn_mr_free(problem); return
        end if
        report%ndof = count(.not. problem%fixed)
        x = merge(data, 0.0_dp, problem%fixed)
        allocate (sx(size(x)))
        call full_apply(problem, x, sx)
        r = merge(problem%rhs - sx, 0.0_dp, .not. problem%fixed)
        scale = norm2(r)
        call newton_mr_solve(problem, x, options%newton, mr%newton)
        report%iterations = mr%newton%minres_iterations
        report%info = mr%newton%status
        call nn_mr_free(problem)
        if (mr%newton%status /= MR_OK) then
            call nn_mr_free(problem); return
        end if
        call full_apply(problem, x, sx)
        r = merge(problem%rhs - sx, 0.0_dp, .not. problem%fixed)
        ! The normalisation of nn_mapped_solve (lifted right-hand side).
        report%residual = norm2(r)/max(tiny(1.0_dp), scale, norm2(merge(problem%rhs, &
            0.0_dp, .not. problem%fixed)), norm2(merge(sx, 0.0_dp, .not. problem%fixed)))
        report%energy = (dot_product(x, sx) - 2*dot_product(x, problem%rhs) &
            + energy_h)/(2*field%mu0)
        bh = csc_matvec(problem%c, x)
        n2 = field%fs%sp(2)%ndof
        field%b(:, 0) = cmplx(bh(:n2), 0.0_dp, dp)
        do n = 1, problem%modes
            field%b(:, n) = cmplx(bh((2*n - 1)*n2 + 1:2*n*n2), &
                -bh(2*n*n2 + 1:(2*n + 1)*n2), dp)/2
        end do
        call nn_mr_free(problem)
        call nn_divergence(field, report%divergence, st)
        mr%total_s = clock() - t0
        if (st /= 0) return
        status = 0
    end subroutine

    subroutine nn_mr_prepare(field, metric, nphi, options, problem, energy_h, mr, &
                             status)
        ! Forcing, fixed traces, the Hessian action and the preconditioner.
        type(free_b_field_t), intent(in), target :: field
        type(map_metric_t), intent(in), target :: metric
        integer, intent(in) :: nphi
        type(nn_mr_options_t), intent(in) :: options
        type(mapped_vacuum_t), intent(out) :: problem
        real(dp), intent(out) :: energy_h
        type(nn_mr_report_t), intent(inout) :: mr
        integer, intent(out) :: status
        type(fortsparse_status_t) :: fst
        real(dp), allocatable :: g(:)
        real(dp) :: volume, t0
        integer :: b, st
        status = 1
        problem%modes = ubound(field%b, 2)
        if (nphi < 2*problem%modes + 1) return
        problem%nb = 2*problem%modes + 1
        problem%n1 = field%fs%sp(1)%ndof; problem%nn = field%fs%nn
        problem%kind = options%preconditioner
        problem%matrix_free = options%matrix_free
        if (problem%kind < NN_MR_PC_NONE .or. problem%kind > NN_MR_PC_LUMPED) return
        call meridian_fourier_curl(field%fs, problem%modes, problem%c, st)
        if (st /= 0) return
        call csc_transpose(problem%c, problem%ct, fst)
        if (.not. status_ok(fst)) return
        call nn_harmonic_forms(field, metric, nphi, g, energy_h, volume, st)
        if (st /= 0) return
        g = field%f_vac*g
        energy_h = field%f_vac**2*energy_h
        problem%rhs = -csc_matvec(problem%ct, g)
        allocate (problem%fixed(problem%nb*problem%n1))
        do b = 0, problem%nb - 1
            problem%fixed(b*problem%n1 + 1:(b + 1)*problem%n1) = field%fs%sp(1)%bnd
            if (b > 0) problem%fixed(b*problem%n1 + problem%nn + 1:(b + 1)*problem%n1) &
                = .true.
        end do
        t0 = clock()
        if (problem%matrix_free) then
            call meridian_fourier_op_init(field%fs, 2, problem%modes, nphi, metric, &
                problem%op, st)
            mr%operator_bytes = meridian_fourier_op_bytes(problem%op)
        else
            call meridian_fourier_curl_curl(field%fs, problem%modes, nphi, metric, &
                problem%s, st)
            mr%operator_bytes = 12_int64*problem%s%nnz
        end if
        mr%operator_s = clock() - t0
        if (st /= 0) return
        t0 = clock()
        if (problem%kind /= NN_MR_PC_NONE) then
            call lumped_blocks(field, metric, merge(options%nphi_lump, nphi, &
                options%nphi_lump > 0), problem, st)
            if (st /= 0) return
        end if
        mr%preconditioner_s = clock() - t0
        status = 0
    end subroutine

    subroutine lumped_blocks(field, metric, nphi, problem, status)
        ! Per-mode operators of the phi-averaged metric Kbar. With M0 the
        ! meridional V2 mass of Kbar, the Fourier V2 mass of Kbar is
        ! blockdiag(M0, M0/2, M0/2, ...) exactly, and with the curl
        ! C_n = [[a, n b], [-n b, a]] on (cos, sin) the lumped operator of
        ! mode n is C_n^T (M0/2 (+) M0/2) C_n; n = 0 is a^T M0 a + eps M_1
        ! on the Nedelec block (gauge regularisation of grad V0).
        type(free_b_field_t), intent(in), target :: field
        type(map_metric_t), intent(in), target :: metric
        integer, intent(in) :: nphi
        type(mapped_vacuum_t), intent(inout) :: problem
        integer, intent(out) :: status
        type(lumped_metric_t) :: lumped
        type(csc_t) :: m0, cn, mn, mc, cnt, sn, sreg
        type(fortsparse_status_t) :: fst
        logical, allocatable :: fixed(:)
        integer :: n, n1, k, st
        status = 1
        n1 = problem%n1
        lumped%base => metric; lumped%nphi = nphi
        call meridian_fourier_mass(field%fs, 2, 0, 1, lumped, m0, st)
        if (st /= 0) return
        allocate (problem%block(0:problem%modes))
        if (problem%kind == NN_MR_PC_JACOBI) allocate (problem%diag(problem%nb*n1))
        do n = 0, problem%modes
            call mode_curl(field%fs, n, cn, st)
            if (st /= 0) return
            call mode_mass(m0, n, mn, st)
            if (st /= 0) return
            call csc_matmul(mn, cn, mc, fst)
            if (.not. status_ok(fst)) return
            call csc_transpose(cn, cnt, fst)
            if (.not. status_ok(fst)) return
            call csc_matmul(cnt, mc, sn, fst)
            if (.not. status_ok(fst)) return
            if (n == 0) then
                call gauge_regularised(sn, field%fs%mass(1), problem%nn, sreg, st)
                if (st /= 0) return
                sn = sreg
                fixed = field%fs%sp(1)%bnd
            else
                fixed = [field%fs%sp(1)%bnd, field%fs%sp(1)%bnd]
                fixed(problem%nn + 1:n1) = .true.
                fixed(n1 + problem%nn + 1:2*n1) = .true.
            end if
            if (problem%kind == NN_MR_PC_JACOBI) then
                k = merge(0, (2*n - 1)*n1, n == 0)
                problem%diag(k + 1:k + size(fixed)) = csc_diagonal(sn)
                cycle
            end if
            call nn_factor(sn, fixed, problem%block(n), st)
            if (st /= 0) return
        end do
        status = 0
    end subroutine

    subroutine mode_curl(fs, n, c, status)
        ! Real curl of mode n: a (n = 0) or [[a, n b], [-n b, a]].
        type(meridian_feec_t), intent(in) :: fs
        integer, intent(in) :: n
        type(csc_t), intent(out) :: c
        integer, intent(out) :: status
        integer, allocatable :: rows(:), cols(:)
        real(dp), allocatable :: vals(:)
        type(fortsparse_status_t) :: fst
        integer :: pos, nr, nc
        status = 1
        associate (a => fs%curl%a, b => fs%curl%b)
            nr = a%nrow; nc = a%ncol
            if (n == 0) then
                c = a; status = 0; return
            end if
            allocate (rows(2*a%nnz + 2*b%nnz), cols(2*a%nnz + 2*b%nnz), &
                vals(2*a%nnz + 2*b%nnz))
            pos = 0
            call put(a, 0, 0, 1.0_dp)
            call put(a, nr, nc, 1.0_dp)
            call put(b, 0, nc, real(n, dp))
            call put(b, nr, 0, -real(n, dp))
            call csc_from_triplet(2*nr, 2*nc, rows(:pos), cols(:pos), vals(:pos), c, fst)
        end associate
        if (status_ok(fst)) status = 0
    contains
        subroutine put(m, ro, co, s)
            type(csc_t), intent(in) :: m
            integer, intent(in) :: ro, co
            real(dp), intent(in) :: s
            integer :: j, p
            do j = 1, m%ncol
                do p = m%col_ptr(j), m%col_ptr(j + 1) - 1
                    pos = pos + 1
                    rows(pos) = m%row_idx(p) + ro; cols(pos) = j + co
                    vals(pos) = s*m%val(p)
                end do
            end do
        end subroutine
    end subroutine

    subroutine mode_mass(m0, n, m, status)
        ! Fourier V2 mass of a phi-independent tensor for mode n: M0 or
        ! blockdiag(M0/2, M0/2) (exact orthogonality, no roundoff couplings).
        type(csc_t), intent(in) :: m0
        integer, intent(in) :: n
        type(csc_t), intent(out) :: m
        integer, intent(out) :: status
        integer, allocatable :: rows(:), cols(:)
        real(dp), allocatable :: vals(:)
        type(fortsparse_status_t) :: fst
        integer :: j, p, pos, b
        status = 1
        if (n == 0) then
            m = m0; status = 0; return
        end if
        allocate (rows(2*m0%nnz), cols(2*m0%nnz), vals(2*m0%nnz))
        pos = 0
        do b = 0, 1
            do j = 1, m0%ncol
                do p = m0%col_ptr(j), m0%col_ptr(j + 1) - 1
                    pos = pos + 1
                    rows(pos) = m0%row_idx(p) + b*m0%nrow; cols(pos) = j + b*m0%ncol
                    vals(pos) = 0.5_dp*m0%val(p)
                end do
            end do
        end do
        call csc_from_triplet(2*m0%nrow, 2*m0%ncol, rows, cols, vals, m, fst)
        if (status_ok(fst)) status = 0
    end subroutine

    subroutine gauge_regularised(a, m, nn, out, status)
        ! a + eps M_1 on the leading nn x nn (Nedelec) block, eps = 1e-8
        ! trace(a)/trace(M_1) over that block (as nn_mapped_solve).
        type(csc_t), intent(in) :: a, m
        integer, intent(in) :: nn
        type(csc_t), intent(out) :: out
        integer, intent(out) :: status
        type(fortsparse_status_t) :: fst
        integer, allocatable :: row(:), col(:)
        real(dp), allocatable :: val(:), da(:), dm(:)
        real(dp) :: eps
        integer :: j, q, cnt
        da = csc_diagonal(a); dm = csc_diagonal(m)
        eps = 1.0e-8_dp*sum(da(:nn))/max(sum(dm(:nn)), tiny(1.0_dp))
        allocate (row(a%nnz + m%nnz), col(a%nnz + m%nnz), val(a%nnz + m%nnz))
        cnt = 0
        do j = 1, a%ncol
            do q = a%col_ptr(j), a%col_ptr(j + 1) - 1
                cnt = cnt + 1
                row(cnt) = a%row_idx(q); col(cnt) = j; val(cnt) = a%val(q)
            end do
        end do
        do j = 1, nn
            do q = m%col_ptr(j), m%col_ptr(j + 1) - 1
                if (m%row_idx(q) > nn) cycle
                cnt = cnt + 1
                row(cnt) = m%row_idx(q); col(cnt) = j; val(cnt) = eps*m%val(q)
            end do
        end do
        call csc_from_triplet(a%nrow, a%ncol, row(:cnt), col(:cnt), val(:cnt), out, fst)
        status = merge(0, 1, status_ok(fst))
    end subroutine

    function csc_diagonal(a) result(d)
        type(csc_t), intent(in) :: a
        real(dp), allocatable :: d(:)
        integer :: j, q
        allocate (d(a%ncol))
        d = 0
        do j = 1, a%ncol
            do q = a%col_ptr(j), a%col_ptr(j + 1) - 1
                if (a%row_idx(q) == j) d(j) = d(j) + a%val(q)
            end do
        end do
    end function

    subroutine full_apply(self, x, y)
        ! y = S x on all dofs (fixed data included).
        class(mapped_vacuum_t), intent(inout) :: self
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: y(:)
        if (self%matrix_free) then
            call meridian_fourier_curl_curl_apply(self%op, self%c, self%ct, x, y)
        else
            y = csc_matvec(self%s, x)
        end if
    end subroutine

    subroutine vacuum_apply(self, x, y)
        ! Hessian on the free dofs: P S P with P the free-dof projection.
        class(mapped_vacuum_t), intent(inout) :: self
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: y(:)
        call full_apply(self, merge(0.0_dp, x, self%fixed), y)
        where (self%fixed) y = 0
    end subroutine

    subroutine vacuum_evaluate(self, x, f, g, status)
        class(mapped_vacuum_t), intent(inout) :: self
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: f, g(:)
        integer, intent(out) :: status
        call full_apply(self, x, g)
        f = 0.5_dp*dot_product(x, g) - dot_product(x, self%rhs)
        g = g - self%rhs
        where (self%fixed) g = 0
        status = merge(0, 1, ieee_is_finite(f))
    end subroutine

    subroutine vacuum_precondition(self, r, z, status)
        class(mapped_vacuum_t), intent(inout) :: self
        real(dp), intent(in) :: r(:)
        real(dp), intent(out) :: z(:)
        integer, intent(out) :: status
        integer :: n, k, m
        status = 0
        select case (self%kind)
        case (NN_MR_PC_NONE)
            z = r
        case (NN_MR_PC_JACOBI)
            z = r/merge(self%diag, 1.0_dp, self%diag > 0)
        case default
            do n = 0, self%modes
                k = merge(0, (2*n - 1)*self%n1, n == 0)
                m = merge(self%n1, 2*self%n1, n == 0)
                call nn_factor_solve_real(self%block(n), r(k + 1:k + m), &
                    z(k + 1:k + m), status)
                if (status /= 0) return
            end do
        end select
        where (self%fixed) z = 0
    end subroutine

    subroutine nn_mr_free(problem)
        ! Release the per-mode preconditioner factors.
        type(mapped_vacuum_t), intent(inout) :: problem
        integer :: n
        if (allocated(problem%block)) then
            do n = lbound(problem%block, 1), ubound(problem%block, 1)
                call nn_factor_free(problem%block(n))
            end do
            deallocate (problem%block)
        end if
    end subroutine

    subroutine lumped_tensor(self, t, x, k, status)
        ! Kbar(R, Z) = mean over nphi angles of the mapped K at (R, Z); the
        ! phi of the query points is ignored.
        class(lumped_metric_t), intent(in) :: self
        integer, intent(in) :: t
        real(dp), intent(in) :: x(:, :)
        real(dp), intent(out) :: k(:, :, :)
        integer, intent(out) :: status
        real(dp), allocatable :: pts(:, :), kt(:, :, :)
        integer :: p, j, np
        status = 1
        k = 0
        if (.not. associated(self%base) .or. self%nphi < 1) return
        np = size(x, 2)
        allocate (pts(3, np*self%nphi), kt(3, 3, np*self%nphi))
        do p = 1, np
            do j = 1, self%nphi
                pts(:, (p - 1)*self%nphi + j) = [x(1:2, p), 2*pi*(j - 1)/self%nphi]
            end do
        end do
        call self%base%tensor(t, pts, kt, status)
        if (status /= 0) return
        do p = 1, np
            k(:, :, p) = sum(kt(:, :, (p - 1)*self%nphi + 1:p*self%nphi), dim=3) &
                /self%nphi
        end do
    end subroutine

    real(dp) function clock()
        integer(int64) :: count, rate
        call system_clock(count, rate)
        clock = real(count, dp)/real(rate, dp)
    end function
end module kin6d_nonnested_newton_mr
