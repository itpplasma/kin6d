module kin6d_nonnested_coupled
    ! NN2 second increment: the axisymmetric (n = 0) coupled resistive
    ! equilibrium on the NN1 FEEC complex. Unknowns, in this order:
    !     x = [A (V1: A_pol Nedelec, a = R A_phi) | lambda (V0) | p (V0) |
    !          Phi (V0) | V],
    ! with B = curl A + F_vac/R e_phi. The residuals are
    !     Ampere   <curl dA, B>/mu0 + <dA, grad lambda> - <dA, J> = 0,
    !     gauge    <A, grad chi> = 0,
    !     pressure <grad v, K grad p> - <v, S> = 0,
    !     current  <grad chi, J> = 0,
    !     loop     int J_phi dR dZ - I_phi = 0,
    ! for all test functions with zero trace, where
    !     J = B x grad p/|B|^2 + g_par b b.E - (g_par b b + g_perp (I - b b))
    !         grad Phi,  E = V/(2 pi R) e_phi,
    !     K = chi_par b b + chi_perp (I - b b).
    ! The current equation is Ampere tested with grad chi, so lambda vanishes
    ! whenever it holds: it is the discrete solvability of weak Ampere.
    ! J and K grad p and their JVPs are the generated leaves coupled_current
    ! and anisotropic_flux. Newton uses that exact Jacobian; the D3 stress
    ! residual (nn_residual) is evaluated afterwards as the independent check.
    ! See docs/equilibrium-nonnested-nn2.md.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_meridian_feec_n, only: meridian_qpoint, meridian_qgrad, &
        meridian_weight_i
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
    public :: nc_jacobian, nc_solve, nc_fixed

    type :: nc_problem_t
        real(dp) :: chi_par = 1, chi_perp = 1, g_par = 1, g_perp = 1
        real(dp) :: current = 0 ! net toroidal current, int J_phi dR dZ
        procedure(meridian_weight_i), pointer, nopass :: source => null()
        procedure(meridian_weight_i), pointer, nopass :: g_par_shape => null()
    end type

    ! Global index of local entry j of a block is offset + j.
    type :: nc_layout_t
        integer :: na = 0, n0 = 0, nb = 0, n = 0
        integer :: ia = 0, il = 0, ip = 0, iphi = 0, iv = 0
    end type

    type :: nc_report_t
        integer :: newton_steps = 0
        real(dp) :: residual(5) = 0     ! blockwise backward errors at exit
        real(dp) :: history(0:40) = -1  ! max blockwise backward error per step
        real(dp) :: voltage = 0, multiplier = 0, divergence = 0
        character(32) :: claim = 'candidate'
    end type

    type :: nc_ops_t
        type(csc_t) :: c, ct, m1, m2, g, gt
    end type

    integer, parameter :: chunk = 256
contains
    function nc_layout(field) result(lay)
        type(free_b_field_t), intent(in) :: field
        type(nc_layout_t) :: lay
        lay%na = field%fs%sp(1)%ndof; lay%n0 = field%fs%n0
        lay%nb = field%fs%sp(2)%ndof
        lay%ia = 0; lay%il = lay%na; lay%ip = lay%na + lay%n0
        lay%iphi = lay%na + 2*lay%n0; lay%iv = lay%na + 3*lay%n0 + 1
        lay%n = lay%iv
    end function

    function nc_fixed(field) result(mask)
        ! Dirichlet dofs: A trace (normal B and toroidal flux), lambda,
        ! p = p_b and Phi = 0 on the boundary. V is free.
        type(free_b_field_t), intent(in) :: field
        logical, allocatable :: mask(:)
        type(nc_layout_t) :: lay
        lay = nc_layout(field)
        allocate (mask(lay%n))
        mask(lay%ia + 1:lay%ia + lay%na) = field%fs%sp(1)%bnd
        mask(lay%il + 1:lay%il + lay%n0) = field%fs%sp(0)%bnd
        mask(lay%ip + 1:lay%ip + lay%n0) = field%fs%sp(0)%bnd
        mask(lay%iphi + 1:lay%iphi + lay%n0) = field%fs%sp(0)%bnd
        mask(lay%iv) = .false.
    end function

    subroutine operators(field, ops, status)
        type(free_b_field_t), intent(in) :: field
        type(nc_ops_t), intent(out) :: ops
        integer, intent(out) :: status
        type(fortsparse_status_t) :: st
        integer :: t, nln
        status = 1
        if (.not. allocated(field%b)) return
        if (ubound(field%b, 2) /= 0) return
        if (field%fs%cyl_r0 > 0) return
        nln = field%fs%nln
        do t = 1, field%fs%mm%topology%n_triangles
            if (any(field%fs%sp(1)%dof(nln + 1:, t) &
                /= field%fs%nn + field%fs%sp(0)%dof(:, t))) return
        end do
        ops%c = field%fs%curl%a; ops%g = field%fs%grad%a
        ops%m1 = field%fs%mass(1); ops%m2 = field%fs%mass(2)
        call csc_transpose(ops%c, ops%ct, st)
        if (.not. status_ok(st)) return
        call csc_transpose(ops%g, ops%gt, st)
        if (.not. status_ok(st)) return
        status = 0
    end subroutine

    subroutine nc_residual(field, problem, x, r, status)
        ! Residual with zero rows at fixed dofs.
        type(free_b_field_t), intent(in) :: field
        type(nc_problem_t), intent(in) :: problem
        real(dp), intent(in) :: x(:)
        real(dp), allocatable, intent(out) :: r(:)
        integer, intent(out) :: status
        type(nc_ops_t) :: ops
        type(nc_layout_t) :: lay
        real(dp), allocatable :: b(:), rl(:, :), jb(:, :), jd(:, :)
        integer :: t, nt, st, failure
        status = 1
        call operators(field, ops, st)
        if (st /= 0) return
        lay = nc_layout(field)
        if (size(x) /= lay%n) return
        if (.not. all(ieee_is_finite(x))) return
        b = csc_matvec(ops%c, x(lay%ia + 1:lay%ia + lay%na))
        allocate (r(lay%n))
        r = 0
        r(lay%ia + 1:lay%ia + lay%na) = &
            csc_matvec(ops%ct, csc_matvec(ops%m2, b))/field%mu0 &
            + csc_matvec(ops%m1, csc_matvec(ops%g, x(lay%il + 1:lay%il + lay%n0)))
        r(lay%il + 1:lay%il + lay%n0) = &
            csc_matvec(ops%gt, csc_matvec(ops%m1, x(lay%ia + 1:lay%ia + lay%na)))
        nt = field%fs%mm%topology%n_triangles
        allocate (rl(local_rows(field), nt), jb(1, 1), jd(1, 1))
        failure = 0
        !$omp parallel do default(shared) private(t, st) firstprivate(jb, jd) &
        !$omp& schedule(static) reduction(max:failure)
        do t = 1, nt
            call element(field, problem, lay, b, x, t, .false., rl(:, t), jb, &
                jd, st)
            failure = max(failure, st)
        end do
        !$omp end parallel do
        if (failure /= 0) return
        do t = 1, nt
            call scatter_residual(field, lay, t, rl(:, t), r)
        end do
        r(lay%iv) = r(lay%iv) - problem%current
        where (nc_fixed(field)) r = 0
        if (.not. all(ieee_is_finite(r))) return
        status = 0
    end subroutine

    subroutine nc_jacobian(field, problem, x, jac, status)
        ! Exact Jacobian of nc_residual, before the Dirichlet reduction.
        type(free_b_field_t), intent(in) :: field
        type(nc_problem_t), intent(in) :: problem
        real(dp), intent(in) :: x(:)
        type(csc_t), intent(out) :: jac
        integer, intent(out) :: status
        type(nc_ops_t) :: ops
        type(nc_layout_t) :: lay
        type(csc_t) :: xb, xbc, kcc, tmp, mg, mgt
        type(fortsparse_status_t) :: st
        real(dp), allocatable :: b(:), rl(:, :), jb(:, :, :), jd(:, :, :)
        real(dp), allocatable :: vals(:), bvals(:)
        integer, allocatable :: rows(:), cols(:), brows(:), bcols(:)
        integer :: t, t0, t1, nt, nr, n2, nd, ierr, failure, count, bcount
        status = 1
        call operators(field, ops, ierr)
        if (ierr /= 0) return
        lay = nc_layout(field)
        if (size(x) /= lay%n) return
        if (.not. all(ieee_is_finite(x))) return
        b = csc_matvec(ops%c, x(lay%ia + 1:lay%ia + lay%na))
        nt = field%fs%mm%topology%n_triangles
        nr = local_rows(field); n2 = field%fs%sp(2)%nloc
        nd = 2*field%fs%nl0 + 1
        allocate (rl(nr, chunk), jb(nr, n2, chunk), jd(nr, nd, chunk))
        allocate (rows(16), cols(16), vals(16), brows(16), bcols(16), bvals(16))
        count = 0; bcount = 0
        do t0 = 1, nt, chunk
            t1 = min(nt, t0 + chunk - 1)
            failure = 0
            !$omp parallel do default(shared) private(t, ierr) schedule(static) &
            !$omp& reduction(max:failure)
            do t = t0, t1
                call element(field, problem, lay, b, x, t, .true., &
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
        call csc_from_triplet(lay%n, lay%nb, brows(:bcount), bcols(:bcount), &
            bvals(:bcount), xb, st)
        if (.not. status_ok(st)) return
        call csc_matmul(xb, ops%c, xbc, st)
        if (.not. status_ok(st)) return
        call append(xbc, 0, lay%ia, 1.0_dp, rows, cols, vals, count)
        call csc_matmul(ops%m2, ops%c, tmp, st)
        if (.not. status_ok(st)) return
        call csc_matmul(ops%ct, tmp, kcc, st)
        if (.not. status_ok(st)) return
        call append(kcc, lay%ia, lay%ia, 1/field%mu0, rows, cols, vals, count)
        call csc_matmul(ops%m1, ops%g, mg, st)
        if (.not. status_ok(st)) return
        call csc_transpose(mg, mgt, st)
        if (.not. status_ok(st)) return
        call append(mg, lay%ia, lay%il, 1.0_dp, rows, cols, vals, count)
        call append(mgt, lay%il, lay%ia, 1.0_dp, rows, cols, vals, count)
        call csc_from_triplet(lay%n, lay%n, rows(:count), cols(:count), &
            vals(:count), jac, st)
        if (.not. status_ok(st)) return
        if (.not. all(ieee_is_finite(jac%val))) return
        status = 0
    end subroutine

    subroutine nc_solve(field, problem, x, report, status, tol, max_steps)
        ! Newton with the exact Jacobian and residual backtracking. On entry
        ! x is the seed; its fixed entries are the boundary data. On exit the
        ! field holds B = curl A (mode 0), psi = a and the transported p.
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
        type(csc_t) :: jac
        real(dp), allocatable :: r(:), dx(:), trial(:), scale(:), rt(:)
        logical, allocatable :: mask(:)
        real(dp) :: target, merit, merit_trial, alpha, w(5)
        integer :: it, steps, st, ls
        status = 1
        target = 1.0e-11_dp; steps = 30
        if (present(tol)) target = tol
        if (present(max_steps)) steps = min(max_steps, 40)
        if (problem%chi_perp <= 0 .or. problem%chi_par < problem%chi_perp) return
        if (problem%g_perp <= 0 .or. problem%g_par < problem%g_perp) return
        if (.not. associated(problem%source)) return
        lay = nc_layout(field)
        if (size(x) /= lay%n) return
        mask = nc_fixed(field)
        do it = 0, steps
            call nc_residual(field, problem, x, r, st)
            if (st /= 0) return
            call nc_jacobian(field, problem, x, jac, st)
            if (st /= 0) return
            scale = abs_matvec(jac, magnitude(x, lay))
            call backward(r, scale, lay, w)
            report%history(it) = maxval(w); report%residual = w
            report%newton_steps = it
            if (maxval(w) <= target) exit
            if (it == steps) return
            call lu_solve(jac, mask, -r, dx, st)
            if (st /= 0) return
            merit = sum(w**2)
            alpha = 1
            do ls = 1, 12
                trial = x + alpha*dx
                call nc_residual(field, problem, trial, rt, st)
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
        report%multiplier = maxval(abs(x(lay%il + 1:lay%il + lay%n0)))
        call nn_divergence(field, report%divergence, st)
        if (st /= 0) return
        status = 0
    end subroutine

    subroutine nc_store(field, x, status)
        type(free_b_field_t), intent(inout) :: field
        real(dp), intent(in) :: x(:)
        integer, intent(out) :: status
        type(nc_layout_t) :: lay
        lay = nc_layout(field)
        status = 1
        if (size(x) /= lay%n) return
        field%b(:, 0) = cmplx(csc_matvec(field%fs%curl%a, &
            x(lay%ia + 1:lay%ia + lay%na)), 0.0_dp, dp)
        field%psi = x(lay%ia + field%fs%nn + 1:lay%ia + lay%na)
        field%p = x(lay%ip + 1:lay%ip + lay%n0)
        status = 0
    end subroutine

    pure integer function local_rows(field) result(nr)
        type(free_b_field_t), intent(in) :: field
        nr = field%fs%sp(1)%nloc + 2*field%fs%nl0 + 1
    end function

    subroutine element(field, problem, lay, b, x, t, want, rl, jb, jd, status)
        ! Local rows [A (n1) | p (n0) | Phi (n0) | V]; jb columns are local V2
        ! coefficients (chained with curl globally), jd columns local p, Phi
        ! and V. Orientation signs are applied when scattering.
        type(free_b_field_t), intent(in) :: field
        type(nc_problem_t), intent(in) :: problem
        type(nc_layout_t), intent(in) :: lay
        real(dp), intent(in) :: b(:), x(:)
        integer, intent(in) :: t
        logical, intent(in) :: want
        real(dp), intent(out) :: rl(:), jb(:, :), jd(:, :)
        integer, intent(out) :: status
        real(dp) :: v2(3, field%fs%sp(2)%nloc), v1(3, field%fs%sp(1)%nloc)
        real(dp) :: val(field%fs%nl0), grad(2, field%fs%nl0), xq(2), wq, rr
        real(dp) :: bl(field%fs%sp(2)%nloc), pl(field%fs%nl0), fl(field%fs%nl0)
        real(dp) :: bq(3), gp(2), gf(2), ephi, gpar, src, cj(3), fq(3)
        real(dp) :: dj(3, 8), dq(3, 5), dcj(3), dfq(3), e(12), two_pi
        integer, parameter :: slot(8) = [1, 2, 3, 4, 6, 7, 9, 10]
        integer :: q, j, k, n1, n0, n2
        status = 1
        rl = 0
        if (want) then
            jb = 0; jd = 0
        end if
        n1 = field%fs%sp(1)%nloc; n0 = field%fs%nl0; n2 = field%fs%sp(2)%nloc
        two_pi = 2*acos(-1.0_dp)
        bl = b(field%fs%sp(2)%dof(:, t))*field%fs%sp(2)%sgn(:, t)
        pl = x(lay%ip + field%fs%sp(0)%dof(:, t))
        fl = x(lay%iphi + field%fs%sp(0)%dof(:, t))
        do q = 1, field%fs%nq
            call meridian_qpoint(field%fs, 2, t, q, xq, wq, v2, status)
            if (status /= 0) return
            call meridian_qpoint(field%fs, 1, t, q, xq, wq, v1, status)
            if (status /= 0) return
            call meridian_qgrad(field%fs, t, q, xq, wq, val, grad, status)
            if (status /= 0) return
            status = 1
            rr = xq(1)
            bq = matmul(v2, bl)
            bq(2) = bq(2) + field%f_vac/rr
            ! |B| = 0 at a quadrature point is a reported failure.
            if (.not. all(ieee_is_finite(bq))) return
            if (dot_product(bq, bq) <= tiny(1.0_dp)) return
            gp = matmul(grad, pl); gf = matmul(grad, fl)
            ephi = x(lay%iv)/(two_pi*rr)
            gpar = problem%g_par
            if (associated(problem%g_par_shape)) gpar = gpar*problem%g_par_shape(xq)
            src = problem%source(xq)
            call coupled_current(bq(1), bq(2), bq(3), gp(1), 0.0_dp, gp(2), &
                gf(1), 0.0_dp, gf(2), ephi, gpar, problem%g_perp, cj(1), cj(2), &
                cj(3))
            call anisotropic_flux(bq(1), bq(2), bq(3), gp(1), 0.0_dp, gp(2), &
                problem%chi_par, problem%chi_perp, fq(1), fq(2), fq(3))
            call add_rows(rl, cj, fq, -val*src)
            if (.not. want) cycle
            ! Unit tangents: B (3), grad p (R, Z), grad Phi (R, Z), ephi.
            do k = 1, 8
                e = 0
                e(slot(k)) = 1
                call coupled_current_jvp(bq(1), bq(2), bq(3), gp(1), 0.0_dp, &
                    gp(2), gf(1), 0.0_dp, gf(2), ephi, gpar, problem%g_perp, &
                    e(1), e(2), e(3), e(4), e(5), e(6), e(7), e(8), e(9), e(10), &
                    e(11), e(12), dj(1, k), dj(2, k), dj(3, k))
            end do
            do k = 1, 5
                e = 0
                e(slot(k)) = 1
                call anisotropic_flux_jvp(bq(1), bq(2), bq(3), gp(1), 0.0_dp, &
                    gp(2), problem%chi_par, problem%chi_perp, e(1), e(2), e(3), &
                    e(4), e(5), e(6), 0.0_dp, 0.0_dp, dq(1, k), dq(2, k), dq(3, k))
            end do
            do j = 1, n2
                dcj = matmul(dj(:, 1:3), v2(:, j)); dfq = matmul(dq(:, 1:3), v2(:, j))
                call add_rows(jb(:, j), dcj, dfq)
            end do
            do j = 1, n0
                dcj = dj(:, 4)*grad(1, j) + dj(:, 5)*grad(2, j)
                dfq = dq(:, 4)*grad(1, j) + dq(:, 5)*grad(2, j)
                call add_rows(jd(:, j), dcj, dfq)
                dcj = dj(:, 6)*grad(1, j) + dj(:, 7)*grad(2, j)
                call add_rows(jd(:, n0 + j), dcj, [0.0_dp, 0.0_dp, 0.0_dp])
            end do
            call add_rows(jd(:, 2*n0 + 1), dj(:, 8)/(two_pi*rr), &
                [0.0_dp, 0.0_dp, 0.0_dp])
        end do
        status = 0
    contains
        subroutine add_rows(col, cur, flux, extra)
            real(dp), intent(inout) :: col(:)
            real(dp), intent(in) :: cur(3), flux(3)
            real(dp), intent(in), optional :: extra(:)
            col(1:n1) = col(1:n1) - wq*matmul(cur, v1)
            col(n1 + 1:n1 + n0) = col(n1 + 1:n1 + n0) &
                + wq*(grad(1, :)*flux(1) + grad(2, :)*flux(3))
            if (present(extra)) col(n1 + 1:n1 + n0) = col(n1 + 1:n1 + n0) + wq*extra
            col(n1 + n0 + 1:n1 + 2*n0) = col(n1 + n0 + 1:n1 + 2*n0) &
                + wq*(grad(1, :)*cur(1) + grad(2, :)*cur(3))
            col(n1 + 2*n0 + 1) = col(n1 + 2*n0 + 1) + wq*cur(2)/rr
        end subroutine
    end subroutine

    subroutine scatter_residual(field, lay, t, rl, r)
        type(free_b_field_t), intent(in) :: field
        type(nc_layout_t), intent(in) :: lay
        integer, intent(in) :: t
        real(dp), intent(in) :: rl(:)
        real(dp), intent(inout) :: r(:)
        integer :: i, n1, n0
        n1 = field%fs%sp(1)%nloc; n0 = field%fs%nl0
        do i = 1, n1
            associate (g => lay%ia + field%fs%sp(1)%dof(i, t))
                r(g) = r(g) + field%fs%sp(1)%sgn(i, t)*rl(i)
            end associate
        end do
        do i = 1, n0
            associate (g => field%fs%sp(0)%dof(i, t))
                r(lay%ip + g) = r(lay%ip + g) + rl(n1 + i)
                r(lay%iphi + g) = r(lay%iphi + g) + rl(n1 + n0 + i)
            end associate
        end do
        r(lay%iv) = r(lay%iv) + rl(n1 + 2*n0 + 1)
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
        integer :: grow(size(jb, 1)), gcol(size(jd, 2))
        real(dp) :: srow(size(jb, 1))
        integer :: i, j, n1, n0, nr
        n1 = field%fs%sp(1)%nloc; n0 = field%fs%nl0; nr = size(jb, 1)
        grow(1:n1) = lay%ia + field%fs%sp(1)%dof(:, t)
        srow(1:n1) = field%fs%sp(1)%sgn(:, t)
        grow(n1 + 1:n1 + n0) = lay%ip + field%fs%sp(0)%dof(:, t)
        grow(n1 + n0 + 1:n1 + 2*n0) = lay%iphi + field%fs%sp(0)%dof(:, t)
        grow(nr) = lay%iv
        srow(n1 + 1:) = 1
        gcol(1:n0) = lay%ip + field%fs%sp(0)%dof(:, t)
        gcol(n0 + 1:2*n0) = lay%iphi + field%fs%sp(0)%dof(:, t)
        gcol(2*n0 + 1) = lay%iv
        call reserve(brows, bcols, bvals, bcount + nr*size(jb, 2))
        do j = 1, size(jb, 2)
            do i = 1, nr
                if (.not. abs(jb(i, j)) > 0) cycle
                bcount = bcount + 1
                brows(bcount) = grow(i); bcols(bcount) = field%fs%sp(2)%dof(j, t)
                bvals(bcount) = srow(i)*field%fs%sp(2)%sgn(j, t)*jb(i, j)
            end do
        end do
        call reserve(rows, cols, vals, count + nr*size(jd, 2))
        do j = 1, size(jd, 2)
            do i = 1, nr
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

    subroutine lu_solve(a, fixed, rhs, x, status)
!$      use omp_lib, only: omp_get_max_threads
        ! General sparse LU (the Jacobian is not symmetric); the increment
        ! vanishes on fixed dofs.
        type(csc_t), intent(in) :: a
        logical, intent(in) :: fixed(:)
        real(dp), intent(in) :: rhs(:)
        real(dp), allocatable, intent(out) :: x(:)
        integer, intent(out) :: status
        type(csc_t) :: reduced
        type(sparse_solver_t) :: solver
        type(fortsparse_status_t) :: st
        integer, allocatable :: row(:), col(:)
        real(dp), allocatable :: val(:), b(:)
        integer :: n, j, p, count
        n = size(rhs); status = 1
        allocate (row(a%nnz + n), col(a%nnz + n), val(a%nnz + n), x(n))
        count = 0
        do j = 1, n
            if (fixed(j)) then
                count = count + 1
                row(count) = j; col(count) = j; val(count) = 1
                cycle
            end if
            do p = a%col_ptr(j), a%col_ptr(j + 1) - 1
                if (fixed(a%row_idx(p))) cycle
                count = count + 1
                row(count) = a%row_idx(p); col(count) = j; val(count) = a%val(p)
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
        b = merge(0.0_dp, rhs, fixed)
        call sparse_solve(solver, b, x, st)
        call sparse_destroy(solver)
        if (.not. status_ok(st)) return
        if (.not. all(ieee_is_finite(x))) return
        status = 0
    end subroutine
end module kin6d_nonnested_coupled
