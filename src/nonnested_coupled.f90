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
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_meridian_feec_n, only: meridian_qpoint, meridian_qgrad, &
        meridian_weight_i
    use fortfem_meridian_fourier_scalar, only: meridian_fourier_tensor_t
    use fortfem_meridian_fourier_curl, only: meridian_fourier_curl, &
        meridian_fourier_curl_curl, meridian_fourier_gradient_coupling
    use fortsparse, only: csc_t, csc_matvec, csc_from_triplet, csc_matmul, &
        csc_transpose, sparse_solver_t, sparse_factor, sparse_solve, &
        sparse_destroy, fortsparse_status_t, status_ok, FORTSPARSE_BACKEND_SUPERLU, &
        FORTSPARSE_BACKEND_MUMPS, FORTSPARSE_SYM_GENERAL, FORTSPARSE_ORDER_ND
    use kin6d_gs_direct, only: gs_direct_available
    use kin6d_equilibrium_nonnested, only: free_b_field_t, nn_divergence
    use kin6d_nonnested_leaves, only: coupled_current, coupled_current_jvp, &
        anisotropic_flux, anisotropic_flux_jvp
    implicit none
    private
    public :: nc_problem_t, nc_layout_t, nc_report_t, nc_layout, nc_residual
    public :: nc_jacobian, nc_solve, nc_fixed, nc_store, nc_reduce
    public :: nc_scalar3_i, nc_continue

    abstract interface
        function nc_scalar3_i(x, phi) result(w)
            ! A declared datum at (R, Z) = x and toroidal angle phi.
            import :: dp
            real(dp), intent(in) :: x(2), phi
            real(dp) :: w
        end function
    end interface

    type :: nc_problem_t
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
        character(32) :: claim = 'candidate'
    end type

    type :: nc_ops_t
        type(csc_t) :: c, s, mg, mgt
        integer :: nphi = 1
    end type

    type, extends(meridian_fourier_tensor_t) :: identity_t
    contains
        procedure :: tensor => identity_tensor
    end type

    integer, parameter :: chunk = 64
contains
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
        type(identity_t) :: id
        integer :: t, nln, modes
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
        call meridian_fourier_curl(field%fs, modes, ops%c, status)
        if (status /= 0) return
        call meridian_fourier_curl_curl(field%fs, modes, ops%nphi, id, ops%s, status)
        if (status /= 0) return
        call meridian_fourier_gradient_coupling(field%fs, modes, ops%nphi, id, &
            ops%mg, status)
        if (status /= 0) return
        status = 1
        call csc_transpose(ops%mg, ops%mgt, st)
        if (.not. status_ok(st)) return
        status = 0
    end subroutine

    subroutine identity_tensor(self, t, x, k, status)
        class(identity_t), intent(in) :: self
        integer, intent(in) :: t
        real(dp), intent(in) :: x(:, :)
        real(dp), intent(out) :: k(:, :, :)
        integer, intent(out) :: status
        integer :: m
        status = 1
        if (t < 1 .or. size(k, 3) /= size(x, 2)) return
        k = 0
        do m = 1, size(x, 2)
            k(1, 1, m) = 1; k(2, 2, m) = 1; k(3, 3, m) = 1
        end do
        status = 0
    end subroutine

    subroutine nc_residual(field, problem, x, r, status)
        ! Residual with zero rows at fixed dofs and tied rows reduced.
        type(free_b_field_t), intent(in) :: field
        type(nc_problem_t), intent(in) :: problem
        real(dp), intent(in) :: x(:)
        real(dp), allocatable, intent(out) :: r(:)
        integer, intent(out) :: status
        type(nc_ops_t) :: ops
        call operators(field, problem, ops, status)
        if (status /= 0) return
        call residual_ops(field, problem, ops, x, r, status)
    end subroutine

    subroutine residual_ops(field, problem, ops, x, r, status)
        type(free_b_field_t), intent(in) :: field
        type(nc_problem_t), intent(in) :: problem
        type(nc_ops_t), intent(in) :: ops
        real(dp), intent(in) :: x(:)
        real(dp), allocatable, intent(out) :: r(:)
        integer, intent(out) :: status
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
        if (failure /= 0) return
        do t = 1, nt
            call scatter_residual(field, lay, t, rl(:, t), r)
        end do
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
            if (failure /= 0) return
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
        real(dp), allocatable :: r(:), dx(:), trial(:), scale(:), rt(:), load(:)
        logical, allocatable :: mask(:)
        integer, allocatable :: rep(:)
        real(dp) :: target, merit, merit_trial, alpha, w(5)
        integer :: it, steps, st, ls
        status = 1
        target = 1.0e-11_dp; steps = 30
        if (present(tol)) target = tol
        if (present(max_steps)) steps = min(max_steps, 40)
        if (problem%chi_perp <= 0 .or. problem%chi_par < problem%chi_perp) return
        if (problem%g_perp <= 0 .or. problem%g_par < problem%g_perp) return
        if (associated(problem%source) .eqv. associated(problem%source_3d)) return
        if (associated(problem%g_par_shape) .and. associated(problem%g_par_shape_3d)) &
            return
        call operators(field, problem, ops, st)
        if (st /= 0) return
        lay = nc_layout(field)
        if (size(x) /= lay%n) return
        mask = nc_fixed(field, problem)
        rep = representative(field, problem)
        if (any(abs(x - x(rep)) > 0)) return
        call source_load(field, problem, lay, ops%nphi, load, st)
        if (st /= 0) return
        do it = 0, steps
            call residual_ops(field, problem, ops, x, r, st)
            if (st /= 0) return
            call jacobian_ops(field, problem, ops, x, jac, st)
            if (st /= 0) return
            ! Normwise backward error |r|/(|J||x| + |load|): the load keeps
            ! the pressure rows scaled when the seed has p = 0.
            scale = abs_matvec(jac, magnitude(x, lay)) + load
            call nc_reduce(field, problem, scale)
            call backward(r, scale, lay, w)
            report%history(it) = maxval(w); report%residual = w
            report%newton_steps = it
            if (maxval(w) <= target) exit
            if (it == steps) return
            call lu_solve(jac, mask, rep, -r, dx, st)
            if (st /= 0) return
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
            if (ls > 12) return
            x = trial
        end do
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
                gpar = problem%g_par
                if (associated(problem%g_par_shape)) then
                    gpar = gpar*problem%g_par_shape(xq)
                else if (associated(problem%g_par_shape_3d)) then
                    gpar = gpar*problem%g_par_shape_3d(xq, phi)
                end if
                src = source_value(problem, xq, phi)
                call trig(lay%modes, phi, tr, dtr)
                bq = matmul(bm, tr)
                bq(2) = bq(2) + field%f_vac/rr
                ! |B| = 0 at a quadrature point is a reported failure.
                if (.not. all(ieee_is_finite(bq))) return
                if (dot_product(bq, bq) <= tiny(1.0_dp)) return
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

    subroutine lu_solve(a, fixed, rep, rhs, x, status)
!$      use omp_lib, only: omp_get_max_threads
        ! General sparse LU (the Jacobian is not symmetric) of the reduced
        ! system: fixed dofs have zero increment, tied dofs share the
        ! increment of their representative (rows summed, columns summed).
        type(csc_t), intent(in) :: a
        logical, intent(in) :: fixed(:)
        integer, intent(in) :: rep(:)
        real(dp), intent(in) :: rhs(:)
        real(dp), allocatable, intent(out) :: x(:)
        integer, intent(out) :: status
        type(csc_t) :: reduced
        type(sparse_solver_t) :: solver
        type(fortsparse_status_t) :: st
        integer, allocatable :: row(:), col(:)
        real(dp), allocatable :: val(:), b(:), y(:)
        integer :: n, j, p, count
        n = size(rhs); status = 1
        allocate (row(a%nnz + n), col(a%nnz + n), val(a%nnz + n), y(n))
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
        if (gs_direct_available()) then
            ! Threaded MUMPS unsymmetric LU, nested dissection.
            solver%backend_id = FORTSPARSE_BACKEND_MUMPS
            solver%options%symmetry = FORTSPARSE_SYM_GENERAL
            solver%options%ordering = FORTSPARSE_ORDER_ND
!$          solver%options%nthreads = omp_get_max_threads()
        else
            solver%backend_id = FORTSPARSE_BACKEND_SUPERLU
        end if
        call sparse_factor(solver, reduced, st)
        if (.not. status_ok(st)) return
        ! rhs is already reduced (tied rows summed, the others zero).
        b = merge(0.0_dp, rhs, fixed)
        call sparse_solve(solver, b, y, st)
        call sparse_destroy(solver)
        if (.not. status_ok(st)) return
        x = y(rep)
        if (.not. all(ieee_is_finite(x))) return
        status = 0
    end subroutine
end module kin6d_nonnested_coupled
