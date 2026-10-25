module kin6d_nonnested_coupled
    ! NN2: the coupled resistive equilibrium on the NN1 FEEC complex times a
    ! real Fourier series in phi with modes 0..N (N = ubound(field%b, 2)).
    ! Every V0/V1/V2 unknown is stored in FortFEM's real block layout
    ! (fortfem_meridian_fourier_curl): block 0 is mode 0, blocks 2n - 1 and
    ! 2n the cos(n phi) and sin(n phi) amplitudes. Unknowns, in this order:
    !     x = [A (V1 blocks: A_pol Nedelec, a = R A_phi) | lambda (V0 blocks) |
    !          p (V0 blocks) | Phi (V0 blocks) | V],
    ! with B = curl A + F_vac/R e_phi. The residuals are
    !     Ampere   <curl dA, B>/mu0 + <dA, grad lambda> - <dA, J> = 0,
    !     gauge    <A, grad chi> = 0,
    !     pressure <grad v, K grad p> - <v, S> = 0,
    !     current  <grad chi, J> = 0,
    !     loop     phi-mean of int J_phi dR dZ - I_phi = 0,
    ! for all test functions with zero trace, where
    !     J = B x grad p/|B|^2 + g_par b b.E - (g_par b b + g_perp (I - b b))
    !         grad Phi,  E = V/(2 pi R) e_phi,
    !     K = chi_par b b + chi_perp (I - b b).
    ! Inner products are R dR dZ times the phi mean on nphi uniform points
    ! (FortFEM's measure). The Coulomb gauge is used in every block: the
    ! current equation is Ampere tested with grad chi, so lambda vanishes
    ! whenever it holds (discrete solvability of weak Ampere).
    ! curl^T M2 curl, M1 grad and the block curl are FortFEM's Fourier forms
    ! with the identity tensor. J and K grad p and their JVPs are the
    ! generated leaves coupled_current and anisotropic_flux; their phi
    ! components carry the toroidal derivatives. A b that varies in phi
    ! couples the modes through these pointwise terms only.
    ! A floating conductor (problem%floating, V0 boundary nodes) ties the
    ! mode-0 lambda and Phi on those nodes to one unknown each. This removes
    ! the harmonic A of a shell domain and imposes zero net current into
    ! the conductor. See docs/equilibrium-nonnested-nn2.md.
    ! Data hooks: source / g_par_shape take (R, Z) and are axisymmetric;
    ! source_3d / g_par_shape_3d (nc_scalar3_i) take (R, Z, phi) and are
    ! evaluated at every phi point of the element quadrature (exactly one of
    ! source and source_3d, at most one g_par shape). The 3D source is loaded
    ! into every Fourier block with the residual's phi rule. Landreman family
    ! 1 uses them: S = -chi_perp lap p = -chi_perp div(J x B) (B.grad p = 0),
    ! Phi = 0 (a nonzero g_perp forbids grad Phi perpendicular to b), and
    ! g_par = g0 w/w_ref with w = R J.B/B_phi so that J.B = g_par V B_phi/
    ! (2 pi R) with V = 2 pi w_ref/g0 (the Pfirsch-Schlueter current included).
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, &
        ieee_quiet_nan
    use fortfem_meridian_feec_n, only: meridian_qpoint, meridian_qgrad, &
        meridian_weight_i
    use fortfem_meridian_fourier_curl, only: meridian_fourier_curl, &
        meridian_fourier_identity_forms
    use fortsparse, only: csc_t, csc_matvec, csc_from_triplet, csc_matmul, &
        csc_transpose, sparse_solver_t, sparse_factor, sparse_solve, sparse_destroy, &
        fortsparse_status_t, status_ok, FORTSPARSE_BACKEND_NATIVE, FORTSPARSE_ORDER_AMD, &
        FORTSPARSE_SYM_GENERAL, FORTSPARSE_ORDER_ND, FORTSPARSE_MEMORY_LOW, &
        native_ordering_available, FORTSPARSE_LU_SUPERNODAL
    use fortnum_krylov, only: real_fgmres_operator, KRYLOV_OK
    use fortnum_fft, only: fortnum_fft_plan_t, fortnum_fft_workspace_t, &
        fft_c2c_plan_init, fft_workspace_init, fft_c2c
    use kin6d_island_closure, only: island_closure_t, island_closure_check, &
        ISLAND_INVALID, ISLAND_NULL_FIELD, ISLAND_BAD_DATA
    use kin6d_gs_direct, only: gs_direct_available, gs_direct_mpi_ready
    use kin6d_equilibrium_nonnested, only: free_b_field_t, nn_divergence
    use kin6d_nonnested_leaves, only: coupled_current, coupled_current_jvp, &
        anisotropic_flux, anisotropic_flux_jvp, anisotropic_flux_vjp, &
        coupled_current_vjp
    implicit none
    private
    public :: nc_problem_t, nc_layout_t, nc_report_t, nc_layout, nc_residual
    public :: nc_jacobian, nc_solve, nc_fixed, nc_store, nc_reduce
    public :: nc_scalar3_i, nc_continue, nc_validate, nc_jvp, nc_vjp

    abstract interface
        function nc_scalar3_i(x, phi) result(w)
            ! A declared datum at (R, Z) = x and toroidal angle phi.
            import :: dp
            real(dp), intent(in) :: x(2), phi
            real(dp) :: w
        end function
    end interface

    type :: nc_problem_t
        type(island_closure_t) :: closure
        real(dp) :: voltage = 0 ! full-torus loop voltage when drive = V
        real(dp) :: chi_par = 1, chi_perp = 1, g_par = 1, g_perp = 1
        real(dp) :: current = 0 ! net toroidal current, phi-mean int J_phi dR dZ
        integer :: nphi = 0     ! phi points; 0 selects 4 N + 1
        procedure(meridian_weight_i), pointer, nopass :: source => null()
        procedure(meridian_weight_i), pointer, nopass :: g_par_shape => null()
        ! Phi-dependent alternatives; at most one of each pair is associated.
        procedure(nc_scalar3_i), pointer, nopass :: source_3d => null()
        procedure(nc_scalar3_i), pointer, nopass :: g_par_shape_3d => null()
        logical, allocatable :: floating(:) ! V0 boundary nodes of a conductor
    end type

    ! Global index of local entry j of block b of a space is offset +
    ! b*(per-mode size) + j. na, nv are whole-space sizes (all blocks).
    type :: nc_layout_t
        integer :: modes = 0, nb = 1, na1 = 0, n0 = 0, nb2 = 0
        integer :: na = 0, nv = 0, n = 0
        integer :: ia = 0, il = 0, ip = 0, iphi = 0, iv = 0
    end type

    type :: nc_report_t
        integer :: newton_steps = 0
        real(dp) :: residual(5) = 0     ! blockwise backward errors at exit
        real(dp) :: history(0:40) = -1  ! max blockwise backward error per step
        real(dp) :: voltage = 0, multiplier = 0, divergence = 0
        integer :: stages = 0, total_steps = 0 ! nc_continue only
        character(160) :: message = ''
        character(32) :: claim = 'candidate'
    end type

    type :: nc_ops_t
        type(csc_t) :: c, s, mg, mgt
        integer :: nphi = 1
    end type

    integer, parameter :: chunk = 64
contains
    subroutine nc_validate(problem, status, message)
        type(nc_problem_t), intent(in) :: problem
        integer, intent(out) :: status
        character(*), intent(out) :: message
        call island_closure_check(problem%closure, status, message)
        if (status /= 0) return
        status = ISLAND_INVALID
        message = 'Require finite drive and parallel >= perpendicular > 0'
        if (.not. all(ieee_is_finite([problem%chi_par, problem%chi_perp, &
            problem%g_par, problem%g_perp, problem%current, problem%voltage]))) return
        if (problem%chi_perp <= 0 .or. problem%chi_par < problem%chi_perp) return
        if (problem%g_perp <= 0 .or. problem%g_par < problem%g_perp) return
        message = 'Declare exactly one pressure source hook (zero source is explicit)'
        if (associated(problem%source) .eqv. associated(problem%source_3d)) return
        message = 'Declare at most one parallel conductivity shape hook'
        if (associated(problem%g_par_shape) .and. associated(problem%g_par_shape_3d)) &
            return
        status = 0
        message = 'A+D regularised candidate; Dirichlet data supplied in the seed'
    end subroutine

    function nc_layout(field) result(lay)
        type(free_b_field_t), intent(in) :: field
        type(nc_layout_t) :: lay
        lay%modes = ubound(field%b, 2); lay%nb = 2*lay%modes + 1
        lay%na1 = field%fs%sp(1)%ndof; lay%n0 = field%fs%n0
        lay%nb2 = field%fs%sp(2)%ndof
        lay%na = lay%na1*lay%nb; lay%nv = lay%n0*lay%nb
        lay%ia = 0; lay%il = lay%na; lay%ip = lay%na + lay%nv
        lay%iphi = lay%na + 2*lay%nv; lay%iv = lay%na + 3*lay%nv + 1
        lay%n = lay%iv
    end function

    function nc_fixed(field, problem) result(mask)
        ! Dirichlet dofs: A trace (normal B and toroidal flux), lambda,
        ! p = p_b and Phi = 0 on the boundary, in every block. V is free.
        ! Floating-conductor nodes of the mode-0 lambda and Phi are free
        ! (tied, see nc_reduce).
        type(free_b_field_t), intent(in) :: field
        type(nc_problem_t), intent(in), optional :: problem
        logical, allocatable :: mask(:)
        type(nc_layout_t) :: lay
        integer :: b
        lay = nc_layout(field)
        allocate (mask(lay%n))
        do b = 0, lay%nb - 1
            mask(lay%ia + b*lay%na1 + 1:lay%ia + (b + 1)*lay%na1) = &
                field%fs%sp(1)%bnd
            mask(lay%il + b*lay%n0 + 1:lay%il + (b + 1)*lay%n0) = field%fs%sp(0)%bnd
            mask(lay%ip + b*lay%n0 + 1:lay%ip + (b + 1)*lay%n0) = field%fs%sp(0)%bnd
            mask(lay%iphi + b*lay%n0 + 1:lay%iphi + (b + 1)*lay%n0) = &
                field%fs%sp(0)%bnd
        end do
        mask(lay%iv) = .false.
        if (.not. present(problem)) return
        if (problem%closure%drive == 'V') mask(lay%iv) = .true.
        if (.not. allocated(problem%floating)) return
        where (problem%floating) mask(lay%il + 1:lay%il + lay%n0) = .false.
        where (problem%floating) mask(lay%iphi + 1:lay%iphi + lay%n0) = .false.
    end function

    function representative(field, problem) result(rep)
        ! rep(i) = i, except tied floating dofs, which map to the first one.
        type(free_b_field_t), intent(in) :: field
        type(nc_problem_t), intent(in) :: problem
        integer, allocatable :: rep(:)
        type(nc_layout_t) :: lay
        integer :: i, first
        lay = nc_layout(field)
        rep = [(i, i=1, lay%n)]
        if (.not. allocated(problem%floating)) return
        first = findloc(problem%floating, .true., 1)
        if (first == 0) return
        do i = 1, lay%n0
            if (.not. problem%floating(i)) cycle
            rep(lay%il + i) = lay%il + first
            rep(lay%iphi + i) = lay%iphi + first
        end do
    end function

    subroutine nc_reduce(field, problem, r)
        ! Rows of tied dofs summed into their representative; the others
        ! become zero. Fixed rows are zero already.
        type(free_b_field_t), intent(in) :: field
        type(nc_problem_t), intent(in) :: problem
        real(dp), intent(inout) :: r(:)
        integer, allocatable :: rep(:)
        integer :: i
        rep = representative(field, problem)
        do i = 1, size(r)
            if (rep(i) == i) cycle
            r(rep(i)) = r(rep(i)) + r(i)
            r(i) = 0
        end do
    end subroutine

    subroutine operators(field, problem, ops, status)
        type(free_b_field_t), intent(in) :: field
        type(nc_problem_t), intent(in) :: problem
        type(nc_ops_t), intent(out) :: ops
        integer, intent(out) :: status
        type(fortsparse_status_t) :: st
        integer :: t, nln, modes
        character(160) :: message
        call nc_validate(problem, status, message)
        if (status /= 0) return
        status = 1
        if (.not. allocated(field%b)) return
        if (field%fs%cyl_r0 > 0) return
        modes = ubound(field%b, 2)
        ops%nphi = 4*modes + 1
        if (problem%nphi > 0) ops%nphi = problem%nphi
        if (ops%nphi < 2*modes + 1) return
        if (allocated(problem%floating)) then
            if (size(problem%floating) /= field%fs%n0) return
            if (any(problem%floating .and. .not. field%fs%sp(0)%bnd)) return
        end if
        nln = field%fs%nln
        do t = 1, field%fs%mm%topology%n_triangles
            if (any(field%fs%sp(1)%dof(nln + 1:, t) &
                /= field%fs%nn + field%fs%sp(0)%dof(:, t))) return
        end do
        ! Exact constant forms avoid mode-pair quadrature storage. Both the
        ! assembled reference and the action retain the same FEEC operators.
        call meridian_fourier_identity_forms(field%fs, modes, &
            ops%c, ops%s, ops%mg, status)
        if (status /= 0) return
        status = 1
        call csc_transpose(ops%mg, ops%mgt, st)
        if (.not. status_ok(st)) return
        status = 0
    end subroutine

    subroutine nc_residual(field, problem, x, r, status, boundary_rows)
        ! Default: fixed rows zero and ties reduced. boundary_rows returns
        ! all weak rows before those reductions and the current drive target.
        type(free_b_field_t), intent(in) :: field
        type(nc_problem_t), intent(in) :: problem
        real(dp), intent(in) :: x(:)
        real(dp), allocatable, intent(out) :: r(:)
        integer, intent(out) :: status
        logical, intent(in), optional :: boundary_rows
        type(nc_ops_t) :: ops
        call operators(field, problem, ops, status)
        if (status /= 0) return
        call residual_ops(field, problem, ops, x, r, status, boundary_rows)
    end subroutine

    subroutine residual_ops(field, problem, ops, x, r, status, boundary_rows)
        type(free_b_field_t), intent(in) :: field
        type(nc_problem_t), intent(in) :: problem
        type(nc_ops_t), intent(in) :: ops
        real(dp), intent(in) :: x(:)
        real(dp), allocatable, intent(out) :: r(:)
        integer, intent(out) :: status
        logical, intent(in), optional :: boundary_rows
        type(nc_layout_t) :: lay
        real(dp), allocatable :: b(:), rl(:, :), jb(:, :), jd(:, :)
        integer :: t, nt, st, failure
        status = 1
        lay = nc_layout(field)
        if (size(x) /= lay%n) return
        if (.not. all(ieee_is_finite(x))) return
        b = csc_matvec(ops%c, x(lay%ia + 1:lay%ia + lay%na))
        allocate (r(lay%n))
        r = 0
        r(lay%ia + 1:lay%ia + lay%na) = &
            csc_matvec(ops%s, x(lay%ia + 1:lay%ia + lay%na))/field%mu0 &
            + csc_matvec(ops%mg, x(lay%il + 1:lay%il + lay%nv))
        r(lay%il + 1:lay%il + lay%nv) = &
            csc_matvec(ops%mgt, x(lay%ia + 1:lay%ia + lay%na))
        nt = field%fs%mm%topology%n_triangles
        allocate (rl(local_rows(field, lay), nt), jb(1, 1), jd(1, 1))
        failure = 0
        !$omp parallel do default(shared) private(t, st) firstprivate(jb, jd) &
        !$omp& schedule(static) reduction(max:failure)
        do t = 1, nt
            call element(field, problem, lay, ops%nphi, b, x, t, .false., &
                rl(:, t), jb, jd, st)
            failure = max(failure, st)
        end do
        !$omp end parallel do
        if (failure /= 0) then
            status = failure
            return
        end if
        do t = 1, nt
            call scatter_residual(field, lay, t, rl(:, t), r)
        end do
        if (present(boundary_rows)) then
            if (boundary_rows) then
                status = 0
                if (.not. all(ieee_is_finite(r))) status = 1
                return
            end if
        end if
        r(lay%iv) = r(lay%iv) - problem%current
        where (nc_fixed(field, problem)) r = 0
        call nc_reduce(field, problem, r)
        if (.not. all(ieee_is_finite(r))) return
        status = 0
    end subroutine

    subroutine nc_jacobian(field, problem, x, jac, status)
        ! Exact Jacobian of the unreduced residual, before the Dirichlet
        ! and tie reductions.
        type(free_b_field_t), intent(in) :: field
        type(nc_problem_t), intent(in) :: problem
        real(dp), intent(in) :: x(:)
        type(csc_t), intent(out) :: jac
        integer, intent(out) :: status
        type(nc_ops_t) :: ops
        call operators(field, problem, ops, status)
        if (status /= 0) return
        call jacobian_ops(field, problem, ops, x, jac, status)
    end subroutine

    subroutine jacobian_ops(field, problem, ops, x, jac, status)
        type(free_b_field_t), intent(in) :: field
        type(nc_problem_t), intent(in) :: problem
        type(nc_ops_t), intent(in) :: ops
        real(dp), intent(in) :: x(:)
        type(csc_t), intent(out) :: jac
        integer, intent(out) :: status
        type(nc_layout_t) :: lay
        type(csc_t) :: xb, xbc
        type(fortsparse_status_t) :: st
        real(dp), allocatable :: b(:), rl(:, :), jb(:, :, :), jd(:, :, :)
        real(dp), allocatable :: vals(:), bvals(:)
        integer, allocatable :: rows(:), cols(:), brows(:), bcols(:)
        integer :: t, t0, t1, nt, nr, n2, nd, ierr, failure, count, bcount
        status = 1
        lay = nc_layout(field)
        if (size(x) /= lay%n) return
        if (.not. all(ieee_is_finite(x))) return
        b = csc_matvec(ops%c, x(lay%ia + 1:lay%ia + lay%na))
        nt = field%fs%mm%topology%n_triangles
        nr = local_rows(field, lay); n2 = field%fs%sp(2)%nloc*lay%nb
        nd = 2*field%fs%nl0*lay%nb + 1
        allocate (rl(nr, chunk), jb(nr, n2, chunk), jd(nr, nd, chunk))
        allocate (rows(16), cols(16), vals(16), brows(16), bcols(16), bvals(16))
        count = 0; bcount = 0
        do t0 = 1, nt, chunk
            t1 = min(nt, t0 + chunk - 1)
            failure = 0
            !$omp parallel do default(shared) private(t, ierr) schedule(static) &
            !$omp& reduction(max:failure)
            do t = t0, t1
                call element(field, problem, lay, ops%nphi, b, x, t, .true., &
                    rl(:, t - t0 + 1), jb(:, :, t - t0 + 1), jd(:, :, t - t0 + 1), &
                    ierr)
                failure = max(failure, ierr)
            end do
            !$omp end parallel do
            if (failure /= 0) then
                status = failure
                return
            end if
            do t = t0, t1
                call scatter_jacobian(field, lay, t, jb(:, :, t - t0 + 1), &
                    jd(:, :, t - t0 + 1), rows, cols, vals, count, brows, bcols, &
                    bvals, bcount)
            end do
        end do
        call csc_from_triplet(lay%n, lay%nb2*lay%nb, brows(:bcount), &
            bcols(:bcount), bvals(:bcount), xb, st)
        if (.not. status_ok(st)) return
        deallocate (brows, bcols, bvals)
        call csc_matmul(xb, ops%c, xbc, st)
        if (.not. status_ok(st)) return
        call append(xbc, 0, lay%ia, 1.0_dp, rows, cols, vals, count)
        call append(ops%s, lay%ia, lay%ia, 1/field%mu0, rows, cols, vals, count)
        call append(ops%mg, lay%ia, lay%il, 1.0_dp, rows, cols, vals, count)
        call append(ops%mgt, lay%il, lay%ia, 1.0_dp, rows, cols, vals, count)
        call csc_from_triplet(lay%n, lay%n, rows(:count), cols(:count), &
            vals(:count), jac, st)
        if (.not. status_ok(st)) return
        if (.not. all(ieee_is_finite(jac%val))) return
        status = 0
    end subroutine

    subroutine nc_solve(field, problem, x, report, status, tol, max_steps)
        ! Newton with the exact Jacobian and residual backtracking. On entry
        ! x is the seed; its fixed entries are the boundary data and tied
        ! entries must agree. On exit the field holds B = curl A (all
        ! modes), psi = a of mode 0 and the transported p blocks.
        ! Tangents and adjoints of the converged state are solves with the
        ! same Jacobian (equation level); the iteration is not differentiated.
        type(free_b_field_t), intent(inout) :: field
        type(nc_problem_t), intent(in) :: problem
        real(dp), intent(inout) :: x(:)
        type(nc_report_t), intent(out) :: report
        integer, intent(out) :: status
        real(dp), intent(in), optional :: tol
        integer, intent(in), optional :: max_steps
        type(nc_layout_t) :: lay
        type(nc_ops_t) :: ops
        type(csc_t) :: jac
        type(sparse_solver_t) :: solver
        real(dp), allocatable :: r(:), dx(:), trial(:), scale(:), rt(:), load(:)
        logical, allocatable :: mask(:)
        integer, allocatable :: rep(:)
        real(dp) :: target, merit, merit_trial, alpha, w(5)
        integer :: it, steps, st, ls
        logical :: factored
        status = 1
        factored = .false.
        target = 1.0e-11_dp; steps = 30
        if (present(tol)) target = tol
        if (present(max_steps)) steps = min(max_steps, 40)
        call nc_validate(problem, status, report%message)
        if (status /= 0) return
        status = 1
        call operators(field, problem, ops, st)
        if (st /= 0) return
        lay = nc_layout(field)
        if (size(x) /= lay%n) return
        if (problem%closure%drive == 'V') x(lay%iv) = problem%voltage
        mask = nc_fixed(field, problem)
        rep = representative(field, problem)
        if (any(abs(x - x(rep)) > 0)) return
        call source_load(field, problem, lay, ops%nphi, load, st)
        if (st /= 0) return
        w = huge(1.0_dp)
        do it = 0, steps
            call residual_ops(field, problem, ops, x, r, st)
            if (st /= 0) exit
            call jacobian_ops(field, problem, ops, x, jac, st)
            if (st /= 0) exit
            ! Normwise backward error |r|/(|J||x| + |load|): the load keeps
            ! the pressure rows scaled when the seed has p = 0.
            scale = abs_matvec(jac, magnitude(x, lay)) + load
            call nc_reduce(field, problem, scale)
            call backward(r, scale, lay, w)
            report%history(it) = maxval(w); report%residual = w
            report%newton_steps = it
            if (maxval(w) <= target) exit
            if (it == steps) exit
            call newton_step(jac, mask, rep, -r, solver, factored, dx, st)
            if (st /= 0) exit
            ! The line search scales by |J|(|x| + |dx|): a block that is zero
            ! in x (p or Phi of a vacuum seed) is scaled by its step.
            scale = abs_matvec(jac, magnitude(abs(x) + abs(dx), lay)) + load
            call nc_reduce(field, problem, scale)
            call backward(r, scale, lay, w)
            merit = sum(w**2)
            alpha = 1
            do ls = 1, 12
                trial = x + alpha*dx
                call residual_ops(field, problem, ops, trial, rt, st)
                if (st == 0) then
                    call backward(rt, scale, lay, w)
                    merit_trial = sum(w**2)
                    if (merit_trial < merit) exit
                end if
                alpha = alpha/2
            end do
            if (ls > 12) then
                if (st == 0) st = 1
                exit
            end if
            x = trial
        end do
        if (factored) call sparse_destroy(solver)
        if (st /= 0) then
            status = st
            select case (st)
            case (ISLAND_NULL_FIELD)
                report%message = &
                    'Null field: A+D requires |B| > 0, with no floor'
            case (ISLAND_BAD_DATA)
                report%message = &
                    'Nonfinite source/field or nonpositive conductivity shape'
            case default
                report%message = 'NN2 assembly or nonlinear solve failed'
            end select
            return
        end if
        if (maxval(w) > target) return
        call nc_store(field, x, st)
        if (st /= 0) return
        report%voltage = x(lay%iv)
        report%multiplier = maxval(abs(x(lay%il + 1:lay%il + lay%nv)))
        call nn_divergence(field, report%divergence, st)
        if (st /= 0) return
        status = 0
    end subroutine

    subroutine nc_continue(field, problem, x, report, status, tol, decades, &
            max_steps)
        ! Natural-parameter continuation in the anisotropy: chi_par/chi_perp
        ! = rc^s and g_par/g_perp = rg^s (rc, rg the targets of problem; the
        ! perpendicular chi_perp and parallel g_par are kept) for s from 0 to
        ! 1. Each stage is nc_solve to tol seeded by the previous converged
        ! state, so the end point is the discrete solution at the target and
        ! does not depend on the path. The Newton basin shrinks like
        ! 1/ratio (the pressure rows carry chi_par b(A).grad p, so a b error
        ! enters times chi_par), while the solution itself moves little with
        ! the ratio. Steps are at most `decades` of the larger ratio (default
        ! 1); a stage that fails or needs more than max_steps (default 12)
        ! Newton steps is retried from the last converged state with half
        ! the step, down to 1/64 decade. report is the last stage's, with
        ! the stage count and the total Newton steps.
        type(free_b_field_t), intent(inout) :: field
        type(nc_problem_t), intent(in) :: problem
        real(dp), intent(inout) :: x(:)
        type(nc_report_t), intent(out) :: report
        integer, intent(out) :: status
        real(dp), intent(in), optional :: tol, decades
        integer, intent(in), optional :: max_steps
        type(nc_problem_t) :: stage
        real(dp), allocatable :: xs(:)
        real(dp) :: span, s, s_ok, ds, ds_max, target
        integer :: steps, stages, total, st
        status = 1
        target = 1.0e-11_dp; ds_max = 1; steps = 12
        if (present(tol)) target = tol
        if (present(decades)) ds_max = decades
        if (present(max_steps)) steps = max_steps
        call nc_validate(problem, status, report%message)
        if (status /= 0) return
        status = 1
        if (ds_max <= 0) return
        if (problem%chi_perp <= 0 .or. problem%chi_par < problem%chi_perp) return
        if (problem%g_perp <= 0 .or. problem%g_par < problem%g_perp) return
        span = max(log10(problem%chi_par/problem%chi_perp), &
            log10(problem%g_par/problem%g_perp))
        ds_max = merge(1.0_dp, min(1.0_dp, ds_max/max(span, tiny(1.0_dp))), &
            span <= 0)
        stage = problem
        s = merge(1.0_dp, 0.0_dp, span <= 0) ! isotropic: one plain solve
        s_ok = 0; ds = ds_max; stages = 0; total = 0
        do
            xs = x
            call set_stage(s)
            call nc_solve(field, stage, xs, report, st, tol=target, max_steps=steps)
            stages = stages + 1; total = total + report%newton_steps
            if (st == 0) then
                x = xs; s_ok = s
                if (s >= 1) exit
                ds = min(ds_max, 2*ds)
            else
                if (s <= 0 .or. ds*span < 1.0_dp/64) then
                    report%stages = stages; report%total_steps = total
                    return
                end if
                ds = ds/2
            end if
            s = min(1.0_dp, s_ok + ds)
        end do
        report%stages = stages; report%total_steps = total
        status = 0
    contains
        subroutine set_stage(t)
            real(dp), intent(in) :: t
            stage%chi_par = problem%chi_perp*(problem%chi_par/problem%chi_perp)**t
            stage%g_perp = problem%g_par/(problem%g_par/problem%g_perp)**t
            if (t >= 1) then
                stage%chi_par = problem%chi_par; stage%g_perp = problem%g_perp
            end if
        end subroutine
    end subroutine

    subroutine source_load(field, problem, lay, nphi, load, status)
        ! |<v, S>| on the pressure rows of every block (the phi rule of the
        ! residual) and |I_phi| on the loop row; zero elsewhere.
        type(free_b_field_t), intent(in) :: field
        type(nc_problem_t), intent(in) :: problem
        type(nc_layout_t), intent(in) :: lay
        integer, intent(in) :: nphi
        real(dp), allocatable, intent(out) :: load(:)
        integer, intent(out) :: status
        real(dp) :: val(field%fs%nl0), grad(2, field%fs%nl0), xq(2), wq
        real(dp) :: tr(lay%nb), dtr(lay%nb), phi
        integer :: t, q, j, bk
        allocate (load(lay%n))
        load = 0
        do t = 1, field%fs%mm%topology%n_triangles
            do q = 1, field%fs%nq
                call meridian_qgrad(field%fs, t, q, xq, wq, val, grad, status)
                if (status /= 0) return
                do j = 1, nphi
                    phi = 2*acos(-1.0_dp)*(j - 1)/nphi
                    call trig(lay%modes, phi, tr, dtr)
                    do bk = 1, lay%nb
                        associate (g => lay%ip + (bk - 1)*lay%n0 &
                                   + field%fs%sp(0)%dof(:, t))
                            load(g) = load(g) + wq/nphi*tr(bk)*val &
                                *source_value(problem, xq, phi)
                        end associate
                    end do
                end do
            end do
        end do
        load = abs(load)
        load(lay%iv) = abs(problem%current)
        where (nc_fixed(field, problem)) load = 0
        status = 0
    end subroutine

    subroutine nc_store(field, x, status)
        ! B modes in NN1's convention B = B_0 + 2 Re(B_n e^{i n phi}), so
        ! B_n = (c - i s)/2 for the cosine/sine blocks c, s.
        type(free_b_field_t), intent(inout) :: field
        real(dp), intent(in) :: x(:)
        integer, intent(out) :: status
        type(nc_layout_t) :: lay
        type(csc_t) :: c
        real(dp), allocatable :: bv(:)
        integer :: n, m
        lay = nc_layout(field)
        status = 1
        if (size(x) /= lay%n) return
        call meridian_fourier_curl(field%fs, lay%modes, c, status)
        if (status /= 0) return
        bv = csc_matvec(c, x(lay%ia + 1:lay%ia + lay%na))
        m = lay%nb2
        field%b(:, 0) = cmplx(bv(:m), 0.0_dp, dp)
        do n = 1, lay%modes
            field%b(:, n) = cmplx(bv((2*n - 1)*m + 1:2*n*m), &
                -bv(2*n*m + 1:(2*n + 1)*m), dp)/2
        end do
        field%psi = x(lay%ia + field%fs%nn + 1:lay%ia + lay%na1)
        field%p = x(lay%ip + 1:lay%ip + lay%nv)
        status = 0
    end subroutine

    real(dp) function source_value(problem, x, phi) result(s)
        type(nc_problem_t), intent(in) :: problem
        real(dp), intent(in) :: x(2), phi
        if (associated(problem%source_3d)) then
            s = problem%source_3d(x, phi)
        else
            s = problem%source(x)
        end if
    end function

    subroutine point_data(problem, xq, phi, bq, src, gpar, status)
        ! The action and reference share the residual's declared-data domain.
        type(nc_problem_t), intent(in) :: problem
        real(dp), intent(in) :: xq(2), phi, bq(3)
        real(dp), intent(out) :: src, gpar
        integer, intent(out) :: status
        gpar = problem%g_par
        if (associated(problem%g_par_shape)) then
            gpar = gpar*problem%g_par_shape(xq)
        else if (associated(problem%g_par_shape_3d)) then
            gpar = gpar*problem%g_par_shape_3d(xq, phi)
        end if
        status = ISLAND_BAD_DATA
        if (.not. ieee_is_finite(gpar)) return
        if (gpar <= 0) return
        src = source_value(problem, xq, phi)
        if (.not. ieee_is_finite(src)) return
        if (.not. all(ieee_is_finite(bq))) return
        status = ISLAND_NULL_FIELD
        if (dot_product(bq, bq) <= tiny(1.0_dp)) return
        status = 0
    end subroutine

    pure integer function local_rows(field, lay) result(nr)
        type(free_b_field_t), intent(in) :: field
        type(nc_layout_t), intent(in) :: lay
        nr = (field%fs%sp(1)%nloc + 2*field%fs%nl0)*lay%nb + 1
    end function

    pure subroutine trig(modes, phi, tr, dtr)
        ! tr = 1, cos phi, sin phi, cos 2 phi, ...; dtr = d tr / d phi.
        integer, intent(in) :: modes
        real(dp), intent(in) :: phi
        real(dp), intent(out) :: tr(:), dtr(:)
        integer :: n
        tr(1) = 1; dtr(1) = 0
        do n = 1, modes
            tr(2*n) = cos(n*phi); tr(2*n + 1) = sin(n*phi)
            dtr(2*n) = -n*tr(2*n + 1); dtr(2*n + 1) = n*tr(2*n)
        end do
    end subroutine

    subroutine element(field, problem, lay, nphi, b, x, t, want, rl, jb, jd, &
                       status)
        ! Local rows [A (n1 per block) | p (n0 per block) | Phi (n0 per
        ! block) | V], block-major. jb columns are local V2 coefficients per
        ! block (chained with the block curl globally), jd columns local p,
        ! Phi per block and V. Orientation signs are applied when scattering.
        type(free_b_field_t), intent(in) :: field
        type(nc_problem_t), intent(in) :: problem
        type(nc_layout_t), intent(in) :: lay
        integer, intent(in) :: nphi
        real(dp), intent(in) :: b(:), x(:)
        integer, intent(in) :: t
        logical, intent(in) :: want
        real(dp), intent(out) :: rl(:), jb(:, :), jd(:, :)
        integer, intent(out) :: status
        real(dp) :: v2(3, field%fs%sp(2)%nloc), v1(3, field%fs%sp(1)%nloc)
        real(dp) :: val(field%fs%nl0), grad(2, field%fs%nl0), xq(2), wq, rr, w
        real(dp) :: bl(field%fs%sp(2)%nloc, lay%nb), pl(field%fs%nl0, lay%nb)
        real(dp) :: fl(field%fs%nl0, lay%nb), bm(3, lay%nb)
        real(dp) :: pv(lay%nb), pg(2, lay%nb), fv(lay%nb), fg(2, lay%nb)
        real(dp) :: tr(lay%nb), dtr(lay%nb), phi
        real(dp) :: bq(3), gp(3), gf(3), ephi, gpar, src, cj(3), fq(3)
        real(dp) :: dj(3, 10), dq(3, 6), e(12), two_pi
        real(dp), allocatable :: bop(:, :), dop(:, :), w1(:, :), jc(:, :), qc(:, :)
        integer :: q, j, k, bk, n1, n0, n2, nb, ra, rp, rf, rv, cb, cp, cf
        status = 1
        rl = 0
        n1 = field%fs%sp(1)%nloc; n0 = field%fs%nl0; n2 = field%fs%sp(2)%nloc
        nb = lay%nb
        ! Row and column offsets of the local blocks.
        ra = 0; rp = n1*nb; rf = rp + n0*nb; rv = rf + n0*nb + 1
        cb = n2*nb; cp = 0; cf = n0*nb
        if (want) then
            jb = 0; jd = 0
        end if
        allocate (bop(3, cb), dop(3, n0*nb), w1(3, n1*nb), &
            jc(3, cb + 2*n0*nb + 1), qc(3, cb + 2*n0*nb + 1))
        two_pi = 2*acos(-1.0_dp)
        do bk = 1, nb
            bl(:, bk) = b((bk - 1)*lay%nb2 + field%fs%sp(2)%dof(:, t)) &
                *field%fs%sp(2)%sgn(:, t)
            pl(:, bk) = x(lay%ip + (bk - 1)*lay%n0 + field%fs%sp(0)%dof(:, t))
            fl(:, bk) = x(lay%iphi + (bk - 1)*lay%n0 + field%fs%sp(0)%dof(:, t))
        end do
        do q = 1, field%fs%nq
            call meridian_qpoint(field%fs, 2, t, q, xq, wq, v2, status)
            if (status /= 0) return
            call meridian_qpoint(field%fs, 1, t, q, xq, wq, v1, status)
            if (status /= 0) return
            call meridian_qgrad(field%fs, t, q, xq, wq, val, grad, status)
            if (status /= 0) return
            status = 1
            rr = xq(1)
            bm = matmul(v2, bl)
            pv = matmul(val, pl); pg = matmul(grad, pl)
            fv = matmul(val, fl); fg = matmul(grad, fl)
            ephi = x(lay%iv)/(two_pi*rr)
            w = wq/nphi
            do j = 1, nphi
                phi = two_pi*(j - 1)/nphi
                call trig(lay%modes, phi, tr, dtr)
                bq = matmul(bm, tr)
                bq(2) = bq(2) + field%f_vac/rr
                call point_data(problem, xq, phi, bq, src, gpar, status)
                if (status /= 0) return
                status = 1
                gp = [dot_product(pg(1, :), tr), dot_product(pv, dtr)/rr, &
                    dot_product(pg(2, :), tr)]
                gf = [dot_product(fg(1, :), tr), dot_product(fv, dtr)/rr, &
                    dot_product(fg(2, :), tr)]
                call coupled_current(bq(1), bq(2), bq(3), gp(1), gp(2), gp(3), &
                    gf(1), gf(2), gf(3), ephi, gpar, problem%g_perp, cj(1), &
                    cj(2), cj(3))
                call anisotropic_flux(bq(1), bq(2), bq(3), gp(1), gp(2), gp(3), &
                    problem%chi_par, problem%chi_perp, fq(1), fq(2), fq(3))
                do bk = 1, nb
                    rl(ra + (bk - 1)*n1 + 1:ra + bk*n1) = &
                        rl(ra + (bk - 1)*n1 + 1:ra + bk*n1) - w*tr(bk)*matmul(cj, v1)
                    rl(rp + (bk - 1)*n0 + 1:rp + bk*n0) = &
                        rl(rp + (bk - 1)*n0 + 1:rp + bk*n0) + w*(tr(bk) &
                        *(grad(1, :)*fq(1) + grad(2, :)*fq(3) - val*src) &
                        + dtr(bk)*val*fq(2)/rr)
                    rl(rf + (bk - 1)*n0 + 1:rf + bk*n0) = &
                        rl(rf + (bk - 1)*n0 + 1:rf + bk*n0) + w*(tr(bk) &
                        *(grad(1, :)*cj(1) + grad(2, :)*cj(3)) &
                        + dtr(bk)*val*cj(2)/rr)
                end do
                rl(rv) = rl(rv) + w*cj(2)/rr
                if (.not. want) cycle
                ! Unit tangents: B (3), grad p (3), grad Phi (3), ephi.
                do k = 1, 10
                    e = 0
                    e(k) = 1
                    call coupled_current_jvp(bq(1), bq(2), bq(3), gp(1), gp(2), &
                        gp(3), gf(1), gf(2), gf(3), ephi, gpar, problem%g_perp, &
                        e(1), e(2), e(3), e(4), e(5), e(6), e(7), e(8), e(9), &
                        e(10), e(11), e(12), dj(1, k), dj(2, k), dj(3, k))
                end do
                do k = 1, 6
                    e = 0
                    e(k) = 1
                    call anisotropic_flux_jvp(bq(1), bq(2), bq(3), gp(1), gp(2), &
                        gp(3), problem%chi_par, problem%chi_perp, e(1), e(2), &
                        e(3), e(4), e(5), e(6), 0.0_dp, 0.0_dp, dq(1, k), &
                        dq(2, k), dq(3, k))
                end do
                ! Trial and test operators: B per V2 basis, the 3D gradient
                ! of a V0 basis, and the V1 test field, each times its trig.
                do bk = 1, nb
                    bop(:, (bk - 1)*n2 + 1:bk*n2) = tr(bk)*v2
                    w1(:, (bk - 1)*n1 + 1:bk*n1) = tr(bk)*v1
                    dop(1, (bk - 1)*n0 + 1:bk*n0) = tr(bk)*grad(1, :)
                    dop(2, (bk - 1)*n0 + 1:bk*n0) = dtr(bk)*val/rr
                    dop(3, (bk - 1)*n0 + 1:bk*n0) = tr(bk)*grad(2, :)
                end do
                jc(:, :cb) = matmul(dj(:, 1:3), bop)
                jc(:, cb + 1:cb + n0*nb) = matmul(dj(:, 4:6), dop)
                jc(:, cb + n0*nb + 1:cb + 2*n0*nb) = matmul(dj(:, 7:9), dop)
                jc(:, cb + 2*n0*nb + 1) = dj(:, 10)/(two_pi*rr)
                qc = 0
                qc(:, :cb) = matmul(dq(:, 1:3), bop)
                qc(:, cb + 1:cb + n0*nb) = matmul(dq(:, 4:6), dop)
                call add_block(ra + 1, ra + n1*nb, -w*transpose(w1), jc)
                call add_block(rp + 1, rp + n0*nb, w*transpose(dop), qc)
                call add_block(rf + 1, rf + n0*nb, w*transpose(dop), jc)
                jb(rv, :) = jb(rv, :) + w*jc(2, :cb)/rr
                jd(rv, :) = jd(rv, :) + w*jc(2, cb + 1:)/rr
            end do
        end do
        status = 0
    contains
        subroutine add_block(r0, r1, test, cols)
            integer, intent(in) :: r0, r1
            real(dp), intent(in) :: test(:, :), cols(:, :)
            jb(r0:r1, :) = jb(r0:r1, :) + matmul(test, cols(:, :cb))
            jd(r0:r1, :) = jd(r0:r1, :) + matmul(test, cols(:, cb + 1:))
        end subroutine
    end subroutine

    subroutine local_maps(field, lay, t, grow, srow, gcol)
        ! Global rows (with signs) and p/Phi/V columns of triangle t.
        type(free_b_field_t), intent(in) :: field
        type(nc_layout_t), intent(in) :: lay
        integer, intent(in) :: t
        integer, intent(out) :: grow(:), gcol(:)
        real(dp), intent(out) :: srow(:)
        integer :: bk, n1, n0, nb, o
        n1 = field%fs%sp(1)%nloc; n0 = field%fs%nl0; nb = lay%nb
        srow = 1
        do bk = 0, nb - 1
            o = bk*n1
            grow(o + 1:o + n1) = lay%ia + bk*lay%na1 + field%fs%sp(1)%dof(:, t)
            srow(o + 1:o + n1) = field%fs%sp(1)%sgn(:, t)
            o = n1*nb + bk*n0
            grow(o + 1:o + n0) = lay%ip + bk*lay%n0 + field%fs%sp(0)%dof(:, t)
            gcol(bk*n0 + 1:(bk + 1)*n0) = grow(o + 1:o + n0)
            o = n1*nb + n0*nb + bk*n0
            grow(o + 1:o + n0) = lay%iphi + bk*lay%n0 + field%fs%sp(0)%dof(:, t)
            gcol(n0*nb + bk*n0 + 1:n0*nb + (bk + 1)*n0) = grow(o + 1:o + n0)
        end do
        grow(size(grow)) = lay%iv
        gcol(size(gcol)) = lay%iv
    end subroutine

    subroutine scatter_residual(field, lay, t, rl, r)
        type(free_b_field_t), intent(in) :: field
        type(nc_layout_t), intent(in) :: lay
        integer, intent(in) :: t
        real(dp), intent(in) :: rl(:)
        real(dp), intent(inout) :: r(:)
        integer :: grow(size(rl)), gcol(2*field%fs%nl0*lay%nb + 1), i
        real(dp) :: srow(size(rl))
        call local_maps(field, lay, t, grow, srow, gcol)
        do i = 1, size(rl)
            r(grow(i)) = r(grow(i)) + srow(i)*rl(i)
        end do
    end subroutine

    subroutine scatter_jacobian(field, lay, t, jb, jd, rows, cols, vals, count, &
                                brows, bcols, bvals, bcount)
        type(free_b_field_t), intent(in) :: field
        type(nc_layout_t), intent(in) :: lay
        integer, intent(in) :: t
        real(dp), intent(in) :: jb(:, :), jd(:, :)
        integer, allocatable, intent(inout) :: rows(:), cols(:), brows(:), bcols(:)
        real(dp), allocatable, intent(inout) :: vals(:), bvals(:)
        integer, intent(inout) :: count, bcount
        integer :: grow(size(jb, 1)), gcol(size(jd, 2)), bcol(size(jb, 2))
        real(dp) :: srow(size(jb, 1)), scol(size(jb, 2))
        integer :: i, j, n2, bk
        call local_maps(field, lay, t, grow, srow, gcol)
        n2 = field%fs%sp(2)%nloc
        do bk = 0, lay%nb - 1
            bcol(bk*n2 + 1:(bk + 1)*n2) = bk*lay%nb2 + field%fs%sp(2)%dof(:, t)
            scol(bk*n2 + 1:(bk + 1)*n2) = field%fs%sp(2)%sgn(:, t)
        end do
        call reserve(brows, bcols, bvals, bcount + size(jb))
        do j = 1, size(jb, 2)
            do i = 1, size(jb, 1)
                if (.not. abs(jb(i, j)) > 0) cycle
                bcount = bcount + 1
                brows(bcount) = grow(i); bcols(bcount) = bcol(j)
                bvals(bcount) = srow(i)*scol(j)*jb(i, j)
            end do
        end do
        call reserve(rows, cols, vals, count + size(jd))
        do j = 1, size(jd, 2)
            do i = 1, size(jd, 1)
                if (.not. abs(jd(i, j)) > 0) cycle
                count = count + 1
                rows(count) = grow(i); cols(count) = gcol(j)
                vals(count) = srow(i)*jd(i, j)
            end do
        end do
    end subroutine

    subroutine append(a, row0, col0, factor, rows, cols, vals, count)
        type(csc_t), intent(in) :: a
        integer, intent(in) :: row0, col0
        real(dp), intent(in) :: factor
        integer, allocatable, intent(inout) :: rows(:), cols(:)
        real(dp), allocatable, intent(inout) :: vals(:)
        integer, intent(inout) :: count
        integer :: j, p
        call reserve(rows, cols, vals, count + a%nnz)
        do j = 1, a%ncol
            do p = a%col_ptr(j), a%col_ptr(j + 1) - 1
                count = count + 1
                rows(count) = row0 + a%row_idx(p); cols(count) = col0 + j
                vals(count) = factor*a%val(p)
            end do
        end do
    end subroutine

    subroutine reserve(rows, cols, vals, need)
        integer, allocatable, intent(inout) :: rows(:), cols(:)
        real(dp), allocatable, intent(inout) :: vals(:)
        integer, intent(in) :: need
        integer, allocatable :: ri(:), ci(:)
        real(dp), allocatable :: vi(:)
        integer :: n
        if (need <= size(rows)) return
        n = max(need, 2*size(rows))
        allocate (ri(n), ci(n), vi(n))
        ri(:size(rows)) = rows; ci(:size(cols)) = cols; vi(:size(vals)) = vals
        call move_alloc(ri, rows); call move_alloc(ci, cols)
        call move_alloc(vi, vals)
    end subroutine

    function abs_matvec(a, x) result(y)
        type(csc_t), intent(in) :: a
        real(dp), intent(in) :: x(:)
        real(dp) :: y(a%nrow)
        integer :: j, p
        y = 0
        do j = 1, a%ncol
            do p = a%col_ptr(j), a%col_ptr(j + 1) - 1
                y(a%row_idx(p)) = y(a%row_idx(p)) + abs(a%val(p)*x(j))
            end do
        end do
    end function

    function magnitude(x, lay) result(m)
        ! |x| plus the block maximum, so that a block that vanishes (the
        ! seed's A_pol or Phi) still sets the scale of the rows it couples to.
        real(dp), intent(in) :: x(:)
        type(nc_layout_t), intent(in) :: lay
        real(dp) :: m(size(x))
        integer :: lo(5), hi(5), k
        lo = [lay%ia, lay%il, lay%ip, lay%iphi, lay%iv - 1] + 1
        hi = [lay%il, lay%ip, lay%iphi, lay%iv - 1, lay%iv]
        do k = 1, 5
            m(lo(k):hi(k)) = abs(x(lo(k):hi(k))) + maxval(abs(x(lo(k):hi(k))))
        end do
    end function

    subroutine backward(r, scale, lay, w)
        ! Blockwise normwise backward error |r_blk| / | |J| |x| |_blk.
        real(dp), intent(in) :: r(:), scale(:)
        type(nc_layout_t), intent(in) :: lay
        real(dp), intent(out) :: w(5)
        integer :: lo(5), hi(5), k
        lo = [lay%ia, lay%il, lay%ip, lay%iphi, lay%iv - 1] + 1
        hi = [lay%il, lay%ip, lay%iphi, lay%iv - 1, lay%iv]
        do k = 1, 5
            w(k) = norm2(r(lo(k):hi(k)))/max(tiny(1.0_dp), &
                norm2(scale(lo(k):hi(k))))
        end do
    end subroutine

    subroutine newton_step(a, fixed, rep, rhs, solver, factored, x, status)
        ! Newton increment from the reduced system: fixed dofs have zero
        ! increment, tied dofs share the increment of their representative
        ! (rows summed, columns summed). The sparse LU of an earlier Jacobian
        ! of the same solve is kept as a right preconditioner of GMRES on the
        ! current exact Jacobian; the step is accepted only when the
        ! residual of the current reduced system falls below gmres_tol
        ! relative to the right-hand side, otherwise the current Jacobian is
        ! factored and the step is the direct solve, as without reuse.
        type(csc_t), intent(in) :: a
        logical, intent(in) :: fixed(:)
        integer, intent(in) :: rep(:)
        real(dp), intent(in) :: rhs(:)
        type(sparse_solver_t), intent(inout) :: solver
        logical, intent(inout) :: factored
        real(dp), allocatable, intent(out) :: x(:)
        integer, intent(out) :: status
        real(dp), parameter :: gmres_tol = 1.0e-10_dp
        integer, parameter :: gmres_max = 20
        type(csc_t) :: reduced
        type(fortsparse_status_t) :: st
        real(dp), allocatable :: b(:), y(:)
        real(dp) :: bnorm, rnorm
        integer :: info, iterations
        status = 1
        call reduced_system(a, fixed, rep, reduced, status)
        if (status /= 0) return
        status = 1
        ! rhs is already reduced (tied rows summed, the others zero).
        b = merge(0.0_dp, rhs, fixed)
        allocate (y(size(b)))
        y = 0
        bnorm = norm2(b)
        if (factored .and. bnorm > 0) then
            call real_fgmres_operator(apply, precondition, b/bnorm, y, gmres_tol, &
                gmres_max, gmres_max, info, iterations, rnorm)
            if (info == KRYLOV_OK .and. all(ieee_is_finite(y))) then
                x = bnorm*y(rep)
                status = 0
                return
            end if
        end if
        if (factored) call sparse_destroy(solver)
        factored = .false.
        call factor_reduced(reduced, solver, status)
        if (status /= 0) return
        status = 1
        factored = .true.
        call sparse_solve(solver, b, y, st)
        if (.not. status_ok(st)) return
        x = y(rep)
        if (.not. all(ieee_is_finite(x))) return
        status = 0
    contains
        subroutine apply(input, output)
            real(dp), intent(in) :: input(:)
            real(dp), intent(out) :: output(:)
            output = csc_matvec(reduced, input)
        end subroutine

        subroutine precondition(input, output)
            real(dp), intent(in) :: input(:)
            real(dp), intent(out) :: output(:)
            type(fortsparse_status_t) :: pst
            call sparse_solve(solver, input, output, pst)
            if (.not. status_ok(pst)) output = ieee_value(1.0_dp, ieee_quiet_nan)
        end subroutine
    end subroutine

    subroutine reduced_system(a, fixed, rep, reduced, status)
        ! The Jacobian with the Dirichlet and tie reductions applied: unit
        ! rows and columns for fixed and non-representative dofs.
        type(csc_t), intent(in) :: a
        logical, intent(in) :: fixed(:)
        integer, intent(in) :: rep(:)
        type(csc_t), intent(out) :: reduced
        integer, intent(out) :: status
        type(fortsparse_status_t) :: st
        integer, allocatable :: row(:), col(:)
        real(dp), allocatable :: val(:)
        integer :: n, j, p, count
        n = size(fixed); status = 1
        allocate (row(a%nnz + n), col(a%nnz + n), val(a%nnz + n))
        count = 0
        do j = 1, n
            if (fixed(j) .or. rep(j) /= j) then
                count = count + 1
                row(count) = j; col(count) = j; val(count) = 1
            end if
            if (fixed(j)) cycle
            do p = a%col_ptr(j), a%col_ptr(j + 1) - 1
                if (fixed(a%row_idx(p))) cycle
                count = count + 1
                row(count) = rep(a%row_idx(p)); col(count) = rep(j)
                val(count) = a%val(p)
            end do
        end do
        call csc_from_triplet(n, n, row(:count), col(:count), val(:count), &
            reduced, st)
        if (.not. status_ok(st)) return
        status = 0
    end subroutine

    subroutine factor_reduced(reduced, solver, status)
!$      use omp_lib, only: omp_get_max_threads
        ! General sparse LU (the Jacobian is not symmetric).
        type(csc_t), intent(in) :: reduced
        type(sparse_solver_t), intent(inout) :: solver
        integer, intent(out) :: status
        type(fortsparse_status_t) :: st
        status = 1
        if (gs_direct_available()) then
            ! The assembled reduced Jacobian requires native unsymmetric LU.
            if (.not. gs_direct_mpi_ready()) return
            solver%backend_id = FORTSPARSE_BACKEND_NATIVE
            solver%options%symmetry = FORTSPARSE_SYM_GENERAL
            solver%options%lu_method = FORTSPARSE_LU_SUPERNODAL
            solver%options%ordering = FORTSPARSE_ORDER_AMD
            solver%options%structural_matching = .true.
            ! ND on the original graph avoids the observed saddle fill growth.
            ! Keep bounded AMD when METIS or its memory policy is unavailable.
            if (native_ordering_available(FORTSPARSE_ORDER_ND) .and. &
                solver%options%host_memory_limit_bytes == 0 .and. &
                solver%options%memory_policy /= FORTSPARSE_MEMORY_LOW) then
                solver%options%ordering = FORTSPARSE_ORDER_ND
                solver%options%structural_matching = .false.
            end if
            solver%options%pivot_tolerance = 0.1_dp
!$          solver%options%nthreads = omp_get_max_threads()
        else
            solver%backend_id = FORTSPARSE_BACKEND_NATIVE
        end if
        call sparse_factor(solver, reduced, st)
        if (.not. status_ok(st)) return
        status = 0
    end subroutine
    ! I3 Fourier actions retained from islands/i34mf 953a01e. Keep them in
    ! the host so native Release links preserve PRIVATE assembly helpers.
    subroutine nc_jvp(field, problem, x, direction, output, status)
        ! Exact directional derivative, with the same phi quadrature as the
        ! dense reference. Fixed rows are zero and floating rows are summed.
        type(free_b_field_t), intent(in) :: field
        type(nc_problem_t), intent(in) :: problem
        real(dp), intent(in) :: x(:), direction(:)
        real(dp), allocatable, intent(out) :: output(:)
        integer, intent(out) :: status
        type(nc_ops_t) :: ops
        call operators(field, problem, ops, status)
        if (status /= 0) return
        call jvp_ops(field, problem, ops, x, direction, output, status)
    end subroutine

    subroutine nc_vjp(field, problem, x, weight, output, status)
        ! Transpose of nc_jvp in the Euclidean coefficient pairing. Reverse
        ! the row reduction before testing and retain boundary input columns.
        type(free_b_field_t), intent(in) :: field
        type(nc_problem_t), intent(in) :: problem
        real(dp), intent(in) :: x(:), weight(:)
        real(dp), allocatable, intent(out) :: output(:)
        integer, intent(out) :: status
        type(nc_ops_t) :: ops
        type(nc_layout_t) :: lay
        type(csc_t) :: ct, st
        type(fortsparse_status_t) :: sparse_status
        type(fortnum_fft_plan_t) :: plan
        type(fortnum_fft_workspace_t) :: work
        real(dp), allocatable :: b(:), w(:), rb(:, :), rs(:, :), dual_b(:)
        integer, allocatable :: rep(:)
        integer :: t, nt, bk, j, failure, local_status, n0, n2
        status = 1
        if (.not. allocated(field%b)) return
        lay = nc_layout(field)
        if (size(x) /= lay%n .or. size(weight) /= lay%n) return
        if (.not. all(ieee_is_finite(x))) return
        if (.not. all(ieee_is_finite(weight))) return
        call operators(field, problem, ops, status)
        if (status /= 0) return
        status = 1
        call csc_transpose(ops%c, ct, sparse_status)
        if (.not. status_ok(sparse_status)) return
        call csc_transpose(ops%s, st, sparse_status)
        if (.not. status_ok(sparse_status)) return
        rep = representative(field, problem)
        w = weight(rep)
        where (nc_fixed(field, problem)) w = 0
        b = csc_matvec(ops%c, x(:lay%na))
        call fft_c2c_plan_init(plan, ops%nphi)
        call fft_workspace_init(work, plan)
        nt = field%fs%mm%topology%n_triangles
        n0 = field%fs%nl0; n2 = field%fs%sp(2)%nloc
        allocate (rb(n2*lay%nb, nt), rs(2*n0*lay%nb + 1, nt))
        failure = 0
        !$omp parallel do default(shared) private(t, local_status) &
        !$omp& firstprivate(work) schedule(static) reduction(max:failure)
        do t = 1, nt
            call element_vjp(field, problem, lay, ops%nphi, b, x, w, t, &
                rb(:, t), rs(:, t), plan, work, local_status)
            failure = max(failure, local_status)
        end do
        !$omp end parallel do
        status = failure
        if (status /= 0) return
        status = 1
        allocate (output(lay%n), dual_b(lay%nb2*lay%nb))
        output = 0; dual_b = 0
        output(:lay%na) = csc_matvec(st, w(:lay%na))/field%mu0 &
            + csc_matvec(ops%mg, w(lay%il + 1:lay%ip))
        output(lay%il + 1:lay%ip) = csc_matvec(ops%mgt, w(:lay%na))
        do t = 1, nt
            do bk = 0, lay%nb - 1
                do j = 1, n2
                    local_status = bk*lay%nb2 + field%fs%sp(2)%dof(j, t)
                    dual_b(local_status) = dual_b(local_status) &
                        + field%fs%sp(2)%sgn(j, t)*rb(bk*n2 + j, t)
                end do
                do j = 1, n0
                    local_status = bk*lay%n0 + field%fs%sp(0)%dof(j, t)
                    output(lay%ip + local_status) = output(lay%ip + local_status) &
                        + rs(bk*n0 + j, t)
                    output(lay%iphi + local_status) = &
                        output(lay%iphi + local_status) + rs((lay%nb + bk)*n0 + j, t)
                end do
            end do
            output(lay%iv) = output(lay%iv) + rs(size(rs, 1), t)
        end do
        output(:lay%na) = output(:lay%na) + csc_matvec(ct, dual_b)
        if (.not. all(ieee_is_finite(output))) return
        status = 0
    end subroutine

    subroutine element_vjp(field, problem, lay, nphi, b, x, weight, t, rb, rs, &
                           plan, work, status)
        type(free_b_field_t), intent(in) :: field
        type(nc_problem_t), intent(in) :: problem
        type(nc_layout_t), intent(in) :: lay
        integer, intent(in) :: nphi, t
        real(dp), intent(in) :: b(:), x(:), weight(:)
        real(dp), intent(out) :: rb(:), rs(:)
        type(fortnum_fft_plan_t), intent(in) :: plan
        type(fortnum_fft_workspace_t), intent(inout) :: work
        integer, intent(out) :: status
        real(dp) :: v1(3, field%fs%sp(1)%nloc), v2(3, field%fs%sp(2)%nloc)
        real(dp) :: val(field%fs%nl0), grad(2, field%fs%nl0), xq(2), wq
        real(dp) :: u(10, nphi), w(10, nphi), gu(10, nphi), gq(8)
        real(dp) :: am(3, lay%nb), wa(3, nphi), gm(9, lay%nb)
        real(dp) :: al(field%fs%sp(1)%nloc, lay%nb), dummy(2), phi, gpar, src
        real(dp) :: dpmean(lay%nb), dfmean(lay%nb)
        integer :: q, j, k, bk, n, n0, n2
        n0 = field%fs%nl0
        n2 = field%fs%sp(2)%nloc
        rb = 0; rs = 0
        do bk = 1, lay%nb
            al(:, bk) = weight((bk - 1)*lay%na1 + field%fs%sp(1)%dof(:, t)) &
                *field%fs%sp(1)%sgn(:, t)
        end do
        do q = 1, field%fs%nq
            call point_series(field, lay, t, q, b, x, u, xq, wq, &
                v1, v2, val, grad, plan, work, status)
            if (status /= 0) return
            call point_series(field, lay, t, q, b, weight, w, xq, wq, &
                v1, v2, val, grad, plan, work, status)
            if (status /= 0) return
            am = matmul(v1, al)
            do k = 1, 3
                call synthesize(am(k, :), wa(k, :), plan, work)
            end do
            u(2, :) = u(2, :) + field%f_vac/xq(1)
            w(7:9, :) = w(7:9, :) - wa
            w(8, :) = w(8, :) + weight(lay%iv)/xq(1)
            do j = 1, nphi
                phi = 2*acos(-1.0_dp)*(j - 1)/nphi
                call point_data(problem, xq, phi, u(:3, j), src, gpar, status)
                if (status /= 0) return
                call coupled_current_vjp(u(1, j), u(2, j), u(3, j), u(4, j), &
                    u(5, j), u(6, j), u(7, j), u(8, j), u(9, j), u(10, j), &
                    gpar, problem%g_perp, w(7, j), w(8, j), w(9, j), &
                    gu(1, j), gu(2, j), gu(3, j), gu(4, j), gu(5, j), gu(6, j), &
                    gu(7, j), gu(8, j), gu(9, j), gu(10, j), dummy(1), dummy(2))
                call anisotropic_flux_vjp(u(1, j), u(2, j), u(3, j), u(4, j), &
                    u(5, j), u(6, j), problem%chi_par, problem%chi_perp, &
                    w(4, j), w(5, j), w(6, j), gq(1), gq(2), gq(3), &
                    gq(4), gq(5), gq(6), gq(7), gq(8))
                gu(:6, j) = gu(:6, j) + gq(:6)
            end do
            do k = 1, 9
                call project(gu(k, :), gm(k, :), plan, work)
            end do
            dpmean = 0; dfmean = 0
            do n = 1, lay%modes
                dpmean(2*n) = -n*gm(5, 2*n + 1)
                dpmean(2*n + 1) = n*gm(5, 2*n)
                dfmean(2*n) = -n*gm(8, 2*n + 1)
                dfmean(2*n + 1) = n*gm(8, 2*n)
            end do
            do bk = 1, lay%nb
                rb((bk - 1)*n2 + 1:bk*n2) = rb((bk - 1)*n2 + 1:bk*n2) &
                    + wq*matmul(gm(:3, bk), v2)
                rs((bk - 1)*n0 + 1:bk*n0) = rs((bk - 1)*n0 + 1:bk*n0) &
                    + wq*(grad(1, :)*gm(4, bk) + grad(2, :)*gm(6, bk) &
                          + val*dpmean(bk)/xq(1))
                rs((lay%nb + bk - 1)*n0 + 1:(lay%nb + bk)*n0) = &
                    rs((lay%nb + bk - 1)*n0 + 1:(lay%nb + bk)*n0) &
                    + wq*(grad(1, :)*gm(7, bk) + grad(2, :)*gm(9, bk) &
                          + val*dfmean(bk)/xq(1))
            end do
            rs(size(rs)) = rs(size(rs)) &
                + wq*sum(gu(10, :))/(nphi*2*acos(-1.0_dp)*xq(1))
        end do
        status = 0
    end subroutine

    subroutine synthesize(coeff, values, plan, work, derivative)
        real(dp), intent(in) :: coeff(:)
        real(dp), intent(out) :: values(:)
        type(fortnum_fft_plan_t), intent(in) :: plan
        type(fortnum_fft_workspace_t), intent(inout) :: work
        logical, intent(in), optional :: derivative
        complex(dp) :: z(size(values)), a
        integer :: n
        z = 0
        z(1) = coeff(1)
        if (present(derivative)) then
            if (derivative) z(1) = 0
        end if
        do n = 1, (size(coeff) - 1)/2
            a = cmplx(coeff(2*n), -coeff(2*n + 1), dp)/2
            if (present(derivative)) then
                if (derivative) a = cmplx(0.0_dp, real(n, dp), dp)*a
            end if
            z(n + 1) = a
            z(size(z) - n + 1) = conjg(a)
        end do
        call fft_c2c(z, 1, plan, work)
        values = real(z, dp)
    end subroutine

    subroutine project(values, coeff, plan, work)
        ! Test coefficients are phi means, not Fourier amplitudes: no factor 2.
        real(dp), intent(in) :: values(:)
        real(dp), intent(out) :: coeff(:)
        type(fortnum_fft_plan_t), intent(in) :: plan
        type(fortnum_fft_workspace_t), intent(inout) :: work
        complex(dp) :: z(size(values))
        integer :: n
        z = cmplx(values, 0.0_dp, dp)
        call fft_c2c(z, -1, plan, work)
        z = z/size(values)
        coeff = 0
        coeff(1) = real(z(1), dp)
        do n = 1, (size(coeff) - 1)/2
            coeff(2*n) = real(z(n + 1), dp)
            coeff(2*n + 1) = -aimag(z(n + 1))
        end do
    end subroutine

    subroutine point_series(field, lay, t, q, b, x, values, xq, wq, &
                            v1, v2, val, grad, plan, work, status)
        ! Geometry/basis evaluation belongs to FortFEM. Here we only gather
        ! coefficients and schedule Fourier transforms of B, grad p, grad Phi.
        type(free_b_field_t), intent(in) :: field
        type(nc_layout_t), intent(in) :: lay
        integer, intent(in) :: t, q
        real(dp), intent(in) :: b(:), x(:)
        real(dp), intent(out) :: values(:, :), xq(2), wq
        real(dp), intent(out) :: v1(:, :), v2(:, :), val(:), grad(:, :)
        type(fortnum_fft_plan_t), intent(in) :: plan
        type(fortnum_fft_workspace_t), intent(inout) :: work
        integer, intent(out) :: status
        real(dp) :: bl(field%fs%sp(2)%nloc, lay%nb)
        real(dp) :: pl(field%fs%nl0, lay%nb), fl(field%fs%nl0, lay%nb)
        real(dp) :: bm(3, lay%nb), pg(2, lay%nb), fg(2, lay%nb)
        integer :: bk, k
        call meridian_qpoint(field%fs, 2, t, q, xq, wq, v2, status)
        if (status /= 0) return
        call meridian_qpoint(field%fs, 1, t, q, xq, wq, v1, status)
        if (status /= 0) return
        call meridian_qgrad(field%fs, t, q, xq, wq, val, grad, status)
        if (status /= 0) return
        do bk = 1, lay%nb
            bl(:, bk) = b((bk - 1)*lay%nb2 + field%fs%sp(2)%dof(:, t)) &
                        *field%fs%sp(2)%sgn(:, t)
            pl(:, bk) = x(lay%ip + (bk - 1)*lay%n0 + field%fs%sp(0)%dof(:, t))
            fl(:, bk) = x(lay%iphi + (bk - 1)*lay%n0 + field%fs%sp(0)%dof(:, t))
        end do
        bm = matmul(v2, bl); pg = matmul(grad, pl); fg = matmul(grad, fl)
        do k = 1, 3
            call synthesize(bm(k, :), values(k, :), plan, work)
        end do
        call synthesize(pg(1, :), values(4, :), plan, work)
        call synthesize(matmul(val, pl), values(5, :), plan, work, .true.)
        values(5, :) = values(5, :)/xq(1)
        call synthesize(pg(2, :), values(6, :), plan, work)
        call synthesize(fg(1, :), values(7, :), plan, work)
        call synthesize(matmul(val, fl), values(8, :), plan, work, .true.)
        values(8, :) = values(8, :)/xq(1)
        call synthesize(fg(2, :), values(9, :), plan, work)
        values(10, :) = x(lay%iv)/(2*acos(-1.0_dp)*xq(1))
    end subroutine

    subroutine jvp_ops(field, problem, ops, x, direction, output, status)
        type(free_b_field_t), intent(in) :: field
        type(nc_problem_t), intent(in) :: problem
        type(nc_ops_t), intent(in) :: ops
        real(dp), intent(in) :: x(:), direction(:)
        real(dp), allocatable, intent(out) :: output(:)
        integer, intent(out) :: status
        type(nc_layout_t) :: lay
        type(fortnum_fft_plan_t) :: plan
        type(fortnum_fft_workspace_t) :: work
        real(dp), allocatable :: b(:), db(:), rl(:, :)
        integer :: t, nt, st, failure
        status = 1
        if (.not. allocated(field%b)) return
        lay = nc_layout(field)
        if (size(x) /= lay%n .or. size(direction) /= lay%n) return
        if (.not. all(ieee_is_finite(x))) return
        if (.not. all(ieee_is_finite(direction))) return
        call fft_c2c_plan_init(plan, ops%nphi)
        call fft_workspace_init(work, plan)
        b = csc_matvec(ops%c, x(:lay%na))
        db = csc_matvec(ops%c, direction(:lay%na))
        allocate (output(lay%n)); output = 0
        output(:lay%na) = csc_matvec(ops%s, direction(:lay%na))/field%mu0 &
                          + csc_matvec(ops%mg, direction(lay%il + 1:lay%ip))
        output(lay%il + 1:lay%ip) = csc_matvec(ops%mgt, direction(:lay%na))
        nt = field%fs%mm%topology%n_triangles
        allocate (rl(local_rows(field, lay), nt)); failure = 0
        !$omp parallel do default(shared) private(t, st) firstprivate(work) &
        !$omp& schedule(static) reduction(max:failure)
        do t = 1, nt
            call element_jvp(field, problem, lay, ops%nphi, b, db, x, direction, &
                             t, rl(:, t), plan, work, st)
            failure = max(failure, st)
        end do
        !$omp end parallel do
        status = failure
        if (status /= 0) return
        status = 1
        do t = 1, nt
            call scatter_residual(field, lay, t, rl(:, t), output)
        end do
        where (nc_fixed(field, problem)) output = 0
        call nc_reduce(field, problem, output)
        if (.not. all(ieee_is_finite(output))) return
        status = 0
    end subroutine

    subroutine element_jvp(field, problem, lay, nphi, b, db, x, dx, t, rl, &
                           plan, work, status)
        type(free_b_field_t), intent(in) :: field
        type(nc_problem_t), intent(in) :: problem
        type(nc_layout_t), intent(in) :: lay
        integer, intent(in) :: nphi, t
        real(dp), intent(in) :: b(:), db(:), x(:), dx(:)
        real(dp), intent(out) :: rl(:)
        type(fortnum_fft_plan_t), intent(in) :: plan
        type(fortnum_fft_workspace_t), intent(inout) :: work
        integer, intent(out) :: status
        real(dp) :: v1(3, field%fs%sp(1)%nloc), v2(3, field%fs%sp(2)%nloc)
        real(dp) :: val(field%fs%nl0), grad(2, field%fs%nl0), xq(2), wq
        real(dp) :: u(10, nphi), du(10, nphi), cj(3, nphi), fq(3, nphi)
        real(dp) :: cm(3, lay%nb), fm(3, lay%nb), dc(lay%nb), df(lay%nb)
        real(dp) :: rr, phi, gpar, twopi, src
        integer :: q, j, k, bk, n1, n0, rp, rf, n
        n1 = field%fs%sp(1)%nloc; n0 = field%fs%nl0
        rp = n1*lay%nb; rf = rp + n0*lay%nb
        twopi = 2*acos(-1.0_dp); rl = 0
        do q = 1, field%fs%nq
            call point_series(field, lay, t, q, b, x, u, xq, wq, &
                              v1, v2, val, grad, plan, work, status)
            if (status /= 0) return
            call point_series(field, lay, t, q, db, dx, du, xq, wq, &
                              v1, v2, val, grad, plan, work, status)
            if (status /= 0) return
            rr = xq(1); u(2, :) = u(2, :) + field%f_vac/rr
            do j = 1, nphi
                phi = twopi*(j - 1)/nphi
                call point_data(problem, xq, phi, u(:3, j), src, gpar, status)
                if (status /= 0) return
                call coupled_current_jvp(u(1, j), u(2, j), u(3, j), u(4, j), u(5, j), &
                    u(6, j), u(7, j), u(8, j), u(9, j), u(10, j), &
                    gpar, problem%g_perp, du(1, j), du(2, j), du(3, j), &
                    du(4, j), du(5, j), du(6, j), du(7, j), &
                    du(8, j), du(9, j), du(10, j), 0.0_dp, 0.0_dp, &
                    cj(1, j), cj(2, j), cj(3, j))
                call anisotropic_flux_jvp(u(1, j), u(2, j), u(3, j), u(4, j), u(5, j), &
                    u(6, j), problem%chi_par, problem%chi_perp, du(1, j), du(2, j), &
                    du(3, j), du(4, j), du(5, j), du(6, j), 0.0_dp, 0.0_dp, &
                    fq(1, j), fq(2, j), fq(3, j))
            end do
            do k = 1, 3
                call project(cj(k, :), cm(k, :), plan, work)
                call project(fq(k, :), fm(k, :), plan, work)
            end do
            dc = 0; df = 0
            do n = 1, lay%modes
                dc(2*n) = -n*cm(2, 2*n + 1); dc(2*n + 1) = n*cm(2, 2*n)
                df(2*n) = -n*fm(2, 2*n + 1); df(2*n + 1) = n*fm(2, 2*n)
            end do
            do bk = 1, lay%nb
                rl((bk - 1)*n1 + 1:bk*n1) = rl((bk - 1)*n1 + 1:bk*n1) &
                                            - wq*matmul(cm(:, bk), v1)
                rl(rp + (bk - 1)*n0 + 1:rp + bk*n0) = &
                    rl(rp + (bk - 1)*n0 + 1:rp + bk*n0) &
                    + wq*(grad(1, :)*fm(1, bk) + grad(2, :)*fm(3, bk) &
                          + val*df(bk)/rr)
                rl(rf + (bk - 1)*n0 + 1:rf + bk*n0) = &
                    rl(rf + (bk - 1)*n0 + 1:rf + bk*n0) &
                    + wq*(grad(1, :)*cm(1, bk) + grad(2, :)*cm(3, bk) &
                          + val*dc(bk)/rr)
            end do
            rl(size(rl)) = rl(size(rl)) + wq*cm(2, 1)/rr
        end do
        status = 0
    end subroutine

end module kin6d_nonnested_coupled
