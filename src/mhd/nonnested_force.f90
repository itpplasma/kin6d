module kin6d_nonnested_force
    ! Projected force on the NN1 real-Fourier FEEC complex
    ! (docs/relaxation-comparators.md, items 3 and 5; docs/nonnested-force.md).
    ! The force of a state is P(J x B), the L2 projection of J x B onto the
    ! velocities u in V2_0 (RT x DG with zero normal trace, all Fourier
    ! blocks) with div u = 0. Since P removes every gradient, P(J x B - grad p)
    ! = P(J x B): the pressure is the multiplier and need not be given.
    ! Two discrete forms of the same projection are provided:
    !   Leray (mixed):  [M2 D^T 0; D 0 l; 0 l^T 0] (u, phi, mu) = (F, 0, 0),
    !       D the coefficient divergence V2 -> V3 and l the n = 0 integral
    !       functional that removes the constant from ker D^T on V2_0;
    !   potential velocity:  u = C psi + sum_j <F, h_j> h_j,
    !       K psi = C^T F on V1_0 (K = C^T M2 C, curl (curl curl)^-1 curl),
    !       deflated by the gauge multiplier lambda in V0_0
    !       [K M1 G; G^T M1 0] (psi, lambda) = (C^T F, 0)
    !   with F = <V2, J x B> the dual load and h_j the M2-orthonormal discrete
    !   harmonic 2-forms of V2_0 (ker D, orthogonal to C V1_0).
    ! On an exact complex ker D in V2_0 = C V1_0 + H, so both forms are the
    ! same M2-orthogonal projection and agree to solver roundoff; D u = 0
    ! holds as a matrix identity for the potential form.
    ! Harmonic forms are computed per geometry (nf_prepare): the toroidal
    ! class (n = 0 DG field, nonzero toroidal flux) always, and for a shell
    ! (inner boundary nodes given) the poloidal class C a_in with a = R A_phi
    ! = 1 on the inner boundary. Each candidate g becomes h = g - C psi(g),
    ! the discrete Hodge projection. A shell also has one relative harmonic
    ! 1-form G phi_in in ker K on V1_0; it is deflated by tying the n = 0
    ! gauge multipliers on the inner boundary to one unknown (the NN2
    ! floating-conductor tie), so that psi is M1-orthogonal to all of ker K.
    ! The Lorentz load uses J = M1^-1 (C^T M2 B + f_vac C^T <V2, e_phi/R>)
    ! in V1_0 (the FEEC curl* of B, as kin6d_nonnested_relax) and B = C A +
    ! f_vac/R e_phi at the quadrature points, nphi-point phi mean.
    ! All results are candidate diagnostics (no enclosure).
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_meridian_feec_n, only: meridian_qpoint, meridian_load
    use fortfem_meridian_fourier_curl, only: meridian_fourier_curl, &
        meridian_fourier_gradient, meridian_fourier_diagonal_mass, &
        meridian_fourier_pairing
    use fortsparse, only: csc_t, csc_matvec, csc_transpose, csc_matmul, &
        csc_from_triplet, sparse_solver_t, sparse_factor, sparse_solve, &
        sparse_destroy, fortsparse_status_t, status_ok, &
        FORTSPARSE_BACKEND_NATIVE, FORTSPARSE_SYM_GENERAL, &
        FORTSPARSE_LU_SUPERNODAL, FORTSPARSE_ORDER_AMD
    use kin6d_equilibrium_nonnested, only: free_b_field_t, nn_factor_t, &
        nn_factor, nn_factor_solve_real, nn_factor_free
    use kin6d_nonnested_mapped, only: map_metric_t, METRIC_IDENTITY, nn_map_t, &
        nn_harmonic_forms, nn_vector_field_i
    implicit none
    private
    public :: nf_system_t, nf_report_t, nf_prepare, nf_free, nf_load
    public :: nf_lorentz_load, nf_project, nf_potential_solve, nf_deflate
    public :: nf_l2_error, nf_field_coefficients

    real(dp), parameter :: pi = acos(-1.0_dp)

    type :: nf_system_t
        integer :: modes = -1, nb = 0, n0 = 0, n1 = 0, n2 = 0, n3 = 0, nn = 0
        integer :: nt = 0, nq = 0, nphi = 0, nh = 0
        type(csc_t) :: c, ct, m1, m2, d, k, m1g
        ! Dirichlet masks (boundary traces) of the Fourier V0/V1/V2 blocks
        ! and the gauge-multiplier representative (shell tie).
        logical, allocatable :: fixed0(:), fixed1(:), fixed2(:)
        integer, allocatable :: rep0(:)
        type(sparse_solver_t) :: pot, ler
        logical :: pot_ready = .false., ler_ready = .false.
        type(nn_factor_t) :: mass1
        real(dp), allocatable :: weight(:), ell(:), harm(:, :), ck(:)
        real(dp), allocatable :: basis1(:, :, :, :), basis2(:, :, :, :)
        real(dp), allocatable :: wq(:, :), xq(:, :, :), trig(:, :)
        integer, allocatable :: dof1(:, :), sgn1(:, :), dof2(:, :), sgn2(:, :)
        real(dp) :: f_vac = 0, length = 1, volume = 0
        ! Wall seconds of the two factorisations (nf_prepare).
        real(dp) :: potential_seconds = 0, leray_seconds = 0
        logical :: ready = .false.
    end type

    type :: nf_report_t
        ! Projected force |u|_M2 of the two forms, their relative difference,
        ! max |D u| relative to |D| |u| (both), harmonic orthogonality
        ! max_j |<h_j, C psi>_M2| / |C psi|_M2, the harmonic part |H^T F|,
        ! the kernel part absorbed by the gauge multiplier relative to
        ! |C^T F| (deflation of the right-hand side), and the scales
        ! <|J|^2 |B|^2>^1/2 and <|B|^2> of a Lorentz load (0 otherwise).
        real(dp) :: force_potential = 0, force_leray = 0, difference = 0
        real(dp) :: divergence_potential = 0, divergence_leray = 0
        real(dp) :: harmonic_orthogonality = 0, harmonic_part = 0
        real(dp) :: kernel_part = 0, jb_scale = 0, b2 = 0
        real(dp) :: relative = 0 ! force_potential / jb_scale (Lorentz load)
        integer :: status = 1
    end type
contains
    subroutine nf_prepare(field, sys, status, nphi, inner, leray)
        ! inner(n0): V0 boundary nodes of the inner wall of a shell domain.
        ! leray (default true): also factor the mixed Leray system (the
        ! comparison path; its saddle-point LU costs far more than the
        ! potential form, see docs/nonnested-force.md).
        type(free_b_field_t), intent(in), target :: field
        type(nf_system_t), intent(out) :: sys
        integer, intent(out) :: status
        integer, intent(in), optional :: nphi
        logical, intent(in), optional :: inner(:), leray
        type(nn_map_t), target :: identity
        type(map_metric_t) :: metric
        type(fortsparse_status_t) :: fst
        type(csc_t) :: g, tmp
        real(dp), allocatable :: gh(:), cand(:, :), a(:), r1(:), psi(:), h(:)
        real(dp) :: henergy, hvolume, x(2), nrm
        integer :: b, t, q, st, ip, n, j, i, ncand, first
        integer(int64) :: c0, c1, rate
        logical :: shell, want_leray
        status = 1
        if (.not. allocated(field%b)) return
        want_leray = .true.
        if (present(leray)) want_leray = leray
        shell = present(inner)
        if (shell) then
            if (size(inner) /= field%fs%n0) return
            if (any(inner .and. .not. field%fs%sp(0)%bnd)) return
            if (.not. any(inner)) return
        end if
        sys%modes = ubound(field%b, 2); sys%nb = 2*sys%modes + 1
        sys%n0 = field%fs%n0; sys%n1 = field%fs%sp(1)%ndof
        sys%n2 = field%fs%sp(2)%ndof; sys%n3 = field%fs%sp(3)%ndof
        sys%nn = field%fs%nn
        sys%nt = field%fs%mm%topology%n_triangles; sys%nq = field%fs%nq
        sys%f_vac = field%f_vac; sys%length = field%radial_scale
        sys%nphi = 4*sys%modes + 1
        if (present(nphi)) sys%nphi = nphi
        if (sys%nphi < 2*sys%modes + 1) return
        allocate (sys%weight(0:sys%nb - 1)); sys%weight = 0.5_dp; sys%weight(0) = 1
        call meridian_fourier_curl(field%fs, sys%modes, sys%c, st)
        if (st /= 0) return
        call csc_transpose(sys%c, sys%ct, fst)
        if (.not. status_ok(fst)) return
        call meridian_fourier_diagonal_mass(field%fs, 1, sys%modes, sys%m1, st)
        if (st /= 0) return
        call meridian_fourier_diagonal_mass(field%fs, 2, sys%modes, sys%m2, st)
        if (st /= 0) return
        call meridian_fourier_gradient(field%fs, sys%modes, g, st)
        if (st /= 0) return
        call csc_matmul(sys%m1, g, sys%m1g, fst)
        if (.not. status_ok(fst)) return
        call csc_matmul(sys%m2, sys%c, tmp, fst)
        if (.not. status_ok(fst)) return
        call csc_matmul(sys%ct, tmp, sys%k, fst)
        if (.not. status_ok(fst)) return
        call divergence_blocks(field, sys%modes, sys%d, st)
        if (st /= 0) return
        allocate (sys%fixed0(sys%nb*sys%n0), sys%fixed1(sys%nb*sys%n1), &
            sys%fixed2(sys%nb*sys%n2), sys%rep0(sys%nb*sys%n0))
        do b = 0, sys%nb - 1
            sys%fixed0(b*sys%n0 + 1:(b + 1)*sys%n0) = field%fs%sp(0)%bnd
            sys%fixed1(b*sys%n1 + 1:(b + 1)*sys%n1) = field%fs%sp(1)%bnd
            sys%fixed2(b*sys%n2 + 1:(b + 1)*sys%n2) = field%fs%sp(2)%bnd
        end do
        sys%rep0 = [(i, i=1, sys%nb*sys%n0)]
        if (shell) then
            ! n = 0 inner-wall multipliers share one free unknown.
            first = findloc(inner, .true., dim=1)
            do i = 1, sys%n0
                if (.not. inner(i)) cycle
                sys%fixed0(i) = .false.; sys%rep0(i) = first
            end do
        end if
        ! Quadrature caches (physical V1/V2 bases, weights, points, trig).
        allocate (sys%basis1(3, field%fs%sp(1)%nloc, sys%nq, sys%nt), &
            sys%basis2(3, field%fs%sp(2)%nloc, sys%nq, sys%nt), &
            sys%wq(sys%nq, sys%nt), sys%xq(2, sys%nq, sys%nt))
        sys%dof1 = field%fs%sp(1)%dof; sys%sgn1 = field%fs%sp(1)%sgn
        sys%dof2 = field%fs%sp(2)%dof; sys%sgn2 = field%fs%sp(2)%sgn
        st = 0
        !$omp parallel do default(shared) private(t, q, x, n) reduction(max:st) &
        !$omp& schedule(static)
        do t = 1, sys%nt
            do q = 1, sys%nq
                call meridian_qpoint(field%fs, 1, t, q, x, sys%wq(q, t), &
                    sys%basis1(:, :, q, t), n)
                st = max(st, n)
                call meridian_qpoint(field%fs, 2, t, q, sys%xq(:, q, t), &
                    sys%wq(q, t), sys%basis2(:, :, q, t), n)
                st = max(st, n)
            end do
        end do
        !$omp end parallel do
        if (st /= 0) return
        sys%volume = sum(sys%wq)
        allocate (sys%trig(0:sys%nb - 1, sys%nphi))
        do ip = 1, sys%nphi
            sys%trig(0, ip) = 1
            do n = 1, sys%modes
                sys%trig(2*n - 1, ip) = cos(n*2*pi*(ip - 1)/sys%nphi)
                sys%trig(2*n, ip) = sin(n*2*pi*(ip - 1)/sys%nphi)
            end do
        end do
        call integral_functional(field, sys, st)
        if (st /= 0) return
        ! Curl* of the analytic harmonic f_vac/R e_phi (identity metric).
        metric%fs => field%fs
        identity%modes = 0
        allocate (identity%d(field%fs%n0, 2, 0:0)); identity%d = 0
        metric%map => identity; metric%kind = METRIC_IDENTITY
        call nn_harmonic_forms(field, metric, max(sys%nphi, 1), gh, henergy, &
            hvolume, st)
        if (st /= 0) return
        sys%ck = field%f_vac*csc_matvec(sys%ct, gh)
        call nn_factor(field%fs%mass(1), field%fs%sp(1)%bnd, sys%mass1, st)
        if (st /= 0) return
        call system_clock(c0, rate)
        call factor_potential(sys, st)
        if (st /= 0) return
        call system_clock(c1)
        sys%potential_seconds = real(c1 - c0, dp)/rate
        if (want_leray) then
            call factor_leray(sys, st)
            if (st /= 0) return
            call system_clock(c0)
            sys%leray_seconds = real(c0 - c1, dp)/rate
        end if
        ! Discrete harmonic 2-forms: Hodge projection of the class candidates.
        ncand = merge(2, 1, shell)
        allocate (cand(sys%nb*sys%n2, ncand)); cand = 0
        cand(sys%nn + 1:sys%n2, 1) = 1 ! toroidal: n = 0 DG field
        if (shell) then
            allocate (a(sys%nb*sys%n1)); a = 0
            where (inner) a(sys%nn + 1:sys%n1) = 1
            cand(:, 2) = csc_matvec(sys%c, a)
            if (maxval(abs(pack(cand(:, 2), sys%fixed2))) > &
                1.0e-12_dp*maxval(abs(cand(:, 2)))) return
            where (sys%fixed2) cand(:, 2) = 0
        end if
        allocate (sys%harm(sys%nb*sys%n2, ncand)); sys%nh = 0
        do j = 1, ncand
            r1 = csc_matvec(sys%ct, csc_matvec(sys%m2, cand(:, j)))
            call nf_potential_solve(sys, r1, psi, st)
            if (st /= 0) return
            h = cand(:, j) - csc_matvec(sys%c, psi)
            do i = 1, sys%nh
                h = h - dot_product(sys%harm(:, i), csc_matvec(sys%m2, h))*sys%harm(:, i)
            end do
            nrm = sqrt(max(dot_product(h, csc_matvec(sys%m2, h)), 0.0_dp))
            if (.not. nrm > 1.0e-8_dp*sqrt(dot_product(cand(:, j), &
                csc_matvec(sys%m2, cand(:, j))))) return
            sys%nh = sys%nh + 1
            sys%harm(:, sys%nh) = h/nrm
        end do
        if (.not. all(ieee_is_finite(sys%harm))) return
        sys%ready = .true.
        status = 0
    end subroutine

    subroutine nf_free(sys)
        type(nf_system_t), intent(inout) :: sys
        if (sys%pot_ready) call sparse_destroy(sys%pot)
        if (sys%ler_ready) call sparse_destroy(sys%ler)
        call nn_factor_free(sys%mass1)
        sys%pot_ready = .false.; sys%ler_ready = .false.; sys%ready = .false.
    end subroutine

    subroutine divergence_blocks(field, modes, d, status)
        ! Real-block divergence V2 -> V3 by the block rule of FortFEM's
        ! curl/gradient (n = 0: a; n > 0: [[a, n b], [-n b, a]]) applied to
        ! div_n = a + i n b. FortFEM exposes no Fourier-block divergence;
        ! this is the only operator built here (candidate for upstream).
        type(free_b_field_t), intent(in) :: field
        integer, intent(in) :: modes
        type(csc_t), intent(out) :: d
        integer, intent(out) :: status
        integer, allocatable :: rows(:), cols(:)
        real(dp), allocatable :: vals(:)
        type(fortsparse_status_t) :: fst
        integer :: nr, nc, pos, n
        nr = field%fs%div%a%nrow; nc = field%fs%div%a%ncol
        allocate (rows(field%fs%div%a%nnz*(2*modes + 1) + 2*modes*field%fs%div%b%nnz))
        allocate (cols(size(rows)), vals(size(rows)))
        pos = 0
        call put(field%fs%div%a, 0, 0, 1.0_dp)
        do n = 1, modes
            call put(field%fs%div%a, (2*n - 1)*nr, (2*n - 1)*nc, 1.0_dp)
            call put(field%fs%div%a, 2*n*nr, 2*n*nc, 1.0_dp)
            call put(field%fs%div%b, (2*n - 1)*nr, 2*n*nc, real(n, dp))
            call put(field%fs%div%b, 2*n*nr, (2*n - 1)*nc, -real(n, dp))
        end do
        call csc_from_triplet((2*modes + 1)*nr, (2*modes + 1)*nc, rows(:pos), &
            cols(:pos), vals(:pos), d, fst)
        status = merge(0, 1, status_ok(fst))
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

    subroutine integral_functional(field, sys, status)
        ! l_i = int g_i dR dZ of the V3 coefficient field g = R div B (the
        ! metric-free integral, = sum_q w_q phi_i(x_q) since the V3 Piola is
        ! td/(det R) and w_q carries det R). D V2_0 = ker l at n = 0.
        type(free_b_field_t), intent(in) :: field
        type(nf_system_t), intent(inout) :: sys
        integer, intent(out) :: status
        real(dp) :: x(2), w, phi(3, field%fs%sp(3)%nloc)
        integer :: t, q
        allocate (sys%ell(sys%n3)); sys%ell = 0
        do t = 1, sys%nt
            do q = 1, sys%nq
                call meridian_qpoint(field%fs, 3, t, q, x, w, phi, status)
                if (status /= 0) return
                sys%ell(field%fs%sp(3)%dof(:, t)) = sys%ell(field%fs%sp(3)%dof(:, t)) &
                    + w*phi(1, :)*field%fs%sp(3)%sgn(:, t)
            end do
        end do
        status = merge(0, 1, all(ieee_is_finite(sys%ell)) .and. maxval(abs(sys%ell)) > 0)
    end subroutine

    subroutine factor_potential(sys, status)
        ! [K M1 G; G^T M1 0] on (psi in V1_0, lambda in V0_0 + tie).
        type(nf_system_t), intent(inout) :: sys
        integer, intent(out) :: status
        integer, allocatable :: row(:), col(:), rep(:)
        real(dp), allocatable :: val(:)
        logical, allocatable :: fixed(:)
        integer :: n1b, n0b, count
        n1b = sys%nb*sys%n1; n0b = sys%nb*sys%n0
        fixed = [sys%fixed1, sys%fixed0]
        rep = [[(count, count=1, n1b)], n1b + sys%rep0]
        allocate (row(sys%k%nnz + 2*sys%m1g%nnz))
        allocate (col(size(row)), val(size(row)))
        count = 0
        call append(sys%k, 0, 0)
        call append(sys%m1g, 0, n1b)
        call append_transpose(sys%m1g, n1b, 0)
        call reduce_and_factor(n1b + n0b, row(:count), col(:count), val(:count), &
            fixed, rep, sys%pot, status)
        sys%pot_ready = status == 0
    contains
        subroutine append(m, ro, co)
            type(csc_t), intent(in) :: m
            integer, intent(in) :: ro, co
            integer :: j, p
            do j = 1, m%ncol
                do p = m%col_ptr(j), m%col_ptr(j + 1) - 1
                    count = count + 1
                    row(count) = m%row_idx(p) + ro; col(count) = j + co
                    val(count) = m%val(p)
                end do
            end do
        end subroutine
        subroutine append_transpose(m, ro, co)
            type(csc_t), intent(in) :: m
            integer, intent(in) :: ro, co
            integer :: j, p
            do j = 1, m%ncol
                do p = m%col_ptr(j), m%col_ptr(j + 1) - 1
                    count = count + 1
                    row(count) = j + ro; col(count) = m%row_idx(p) + co
                    val(count) = m%val(p)
                end do
            end do
        end subroutine
    end subroutine

    subroutine factor_leray(sys, status)
        ! [M2 D^T 0; D 0 l; 0 l^T 0] on (u in V2_0, phi in V3, mu).
        type(nf_system_t), intent(inout) :: sys
        integer, intent(out) :: status
        integer, allocatable :: row(:), col(:), rep(:)
        real(dp), allocatable :: val(:)
        logical, allocatable :: fixed(:)
        integer :: n2b, n3b, count, j, p, n
        n2b = sys%nb*sys%n2; n3b = sys%nb*sys%n3
        n = n2b + n3b + 1
        fixed = [sys%fixed2, spread(.false., 1, n3b + 1)]
        rep = [(j, j=1, n)]
        allocate (row(sys%m2%nnz + 2*sys%d%nnz + 2*sys%n3))
        allocate (col(size(row)), val(size(row)))
        count = 0
        do j = 1, sys%m2%ncol
            do p = sys%m2%col_ptr(j), sys%m2%col_ptr(j + 1) - 1
                count = count + 1
                row(count) = sys%m2%row_idx(p); col(count) = j
                val(count) = sys%m2%val(p)
            end do
        end do
        do j = 1, sys%d%ncol
            do p = sys%d%col_ptr(j), sys%d%col_ptr(j + 1) - 1
                count = count + 1
                row(count) = n2b + sys%d%row_idx(p); col(count) = j
                val(count) = sys%d%val(p)
                count = count + 1
                row(count) = j; col(count) = n2b + sys%d%row_idx(p)
                val(count) = sys%d%val(p)
            end do
        end do
        do j = 1, sys%n3
            count = count + 1
            row(count) = n2b + j; col(count) = n; val(count) = sys%ell(j)
            count = count + 1
            row(count) = n; col(count) = n2b + j; val(count) = sys%ell(j)
        end do
        call reduce_and_factor(n, row(:count), col(:count), val(:count), fixed, &
            rep, sys%ler, status)
        sys%ler_ready = status == 0
    end subroutine

    subroutine reduce_and_factor(n, row, col, val, fixed, rep, solver, status)
        ! Identity rows/columns at fixed and non-representative dofs, tied
        ! rows and columns summed into their representative; general LU.
        integer, intent(in) :: n, row(:), col(:), rep(:)
        real(dp), intent(in) :: val(:)
        logical, intent(in) :: fixed(:)
        type(sparse_solver_t), intent(inout) :: solver
        integer, intent(out) :: status
        integer, allocatable :: r(:), c(:)
        real(dp), allocatable :: v(:)
        type(csc_t) :: a
        type(fortsparse_status_t) :: st
        integer :: i, count
        allocate (r(size(row) + n), c(size(row) + n), v(size(row) + n))
        count = 0
        do i = 1, n
            if (fixed(i) .or. rep(i) /= i) then
                count = count + 1
                r(count) = i; c(count) = i; v(count) = 1
            end if
        end do
        do i = 1, size(row)
            if (fixed(row(i)) .or. fixed(col(i))) cycle
            count = count + 1
            r(count) = rep(row(i)); c(count) = rep(col(i)); v(count) = val(i)
        end do
        call csc_from_triplet(n, n, r(:count), c(:count), v(:count), a, st)
        status = 1
        if (.not. status_ok(st)) return
        solver%backend_id = FORTSPARSE_BACKEND_NATIVE
        solver%options%symmetry = FORTSPARSE_SYM_GENERAL
        solver%options%lu_method = FORTSPARSE_LU_SUPERNODAL
        solver%options%ordering = FORTSPARSE_ORDER_AMD
        solver%options%structural_matching = .false.
        solver%options%pivot_tolerance = 0.1_dp
        call sparse_factor(solver, a, st)
        if (status_ok(st)) status = 0
    end subroutine

    subroutine nf_potential_solve(sys, rhs, psi, status, kernel_part)
        ! psi in V1_0, M1-orthogonal to ker K (gradients and the shell's
        ! relative harmonic 1-form), with K psi = rhs - (kernel part of rhs).
        ! kernel_part: |M1 G lambda| / |rhs| (zero for a consistent rhs).
        type(nf_system_t), intent(inout) :: sys
        real(dp), intent(in) :: rhs(:)
        real(dp), allocatable, intent(out) :: psi(:)
        integer, intent(out) :: status
        real(dp), intent(out), optional :: kernel_part
        real(dp), allocatable :: b(:), y(:), lam(:)
        type(fortsparse_status_t) :: st
        integer :: n1b, n0b, i
        status = 1
        if (.not. sys%pot_ready) return
        n1b = sys%nb*sys%n1; n0b = sys%nb*sys%n0
        if (size(rhs) /= n1b) return
        allocate (b(n1b + n0b), y(n1b + n0b)); b = 0
        b(:n1b) = merge(0.0_dp, rhs, sys%fixed1)
        call sparse_solve(sys%pot, b, y, st)
        if (.not. status_ok(st)) return
        psi = merge(0.0_dp, y(:n1b), sys%fixed1)
        if (present(kernel_part)) then
            allocate (lam(n0b))
            do i = 1, n0b
                lam(i) = y(n1b + sys%rep0(i))
                if (sys%fixed0(i)) lam(i) = 0
            end do
            kernel_part = norm2(merge(0.0_dp, csc_matvec(sys%m1g, lam), sys%fixed1)) &
                /max(norm2(b(:n1b)), tiny(1.0_dp))
        end if
        if (all(ieee_is_finite(psi))) status = 0
    end subroutine

    function nf_deflate(sys, v) result(w)
        ! v minus its M2-orthogonal harmonic part (a V2 coefficient vector).
        type(nf_system_t), intent(in) :: sys
        real(dp), intent(in) :: v(:)
        real(dp), allocatable :: w(:)
        real(dp), allocatable :: mv(:)
        integer :: j
        w = v
        mv = csc_matvec(sys%m2, v)
        do j = 1, sys%nh
            w = w - dot_product(sys%harm(:, j), mv)*sys%harm(:, j)
        end do
    end function

    subroutine nf_project(sys, load, u_potential, u_leray, report)
        ! Both forms of P for the dual V2 load F (nb*n2, boundary entries
        ! ignored), with the diagnostics of nf_report_t.
        type(nf_system_t), intent(inout) :: sys
        real(dp), intent(in) :: load(:)
        real(dp), allocatable, intent(out) :: u_potential(:), u_leray(:)
        type(nf_report_t), intent(inout) :: report
        real(dp), allocatable :: f(:), r1(:), psi(:), cpsi(:), b(:), y(:), hf(:)
        type(fortsparse_status_t) :: st
        real(dp) :: dnorm, kp
        integer :: n2b, j
        report%status = 1
        if (.not. sys%ready) return
        n2b = sys%nb*sys%n2
        if (size(load) /= n2b) return
        f = merge(0.0_dp, load, sys%fixed2)
        ! Potential velocity: curl (curl curl)^-1 curl plus harmonic part.
        r1 = csc_matvec(sys%ct, f)
        call nf_potential_solve(sys, r1, psi, report%status, kp)
        if (report%status /= 0) return
        report%status = 1
        cpsi = csc_matvec(sys%c, psi)
        allocate (hf(sys%nh))
        u_potential = cpsi
        do j = 1, sys%nh
            hf(j) = dot_product(sys%harm(:, j), f)
            u_potential = u_potential + hf(j)*sys%harm(:, j)
        end do
        report%force_potential = m2norm(sys, u_potential)
        ! Leray: mixed RT/DG saddle point (when factored; otherwise u_leray
        ! is the potential form and difference = -1).
        if (sys%ler_ready) then
            allocate (b(n2b + sys%nb*sys%n3 + 1)); b = 0
            allocate (y(size(b)))
            b(:n2b) = f
            call sparse_solve(sys%ler, b, y, st)
            if (.not. status_ok(st)) return
            u_leray = merge(0.0_dp, y(:n2b), sys%fixed2)
            report%force_leray = m2norm(sys, u_leray)
            report%difference = m2norm(sys, u_potential - u_leray) &
                /max(report%force_potential, report%force_leray, tiny(1.0_dp))
        else
            u_leray = u_potential
            report%force_leray = report%force_potential
            report%difference = -1
        end if
        dnorm = d_scale(sys)
        report%divergence_potential = maxval(abs(csc_matvec(sys%d, u_potential))) &
            /max(dnorm*maxval(abs(u_potential)), tiny(1.0_dp))
        report%divergence_leray = maxval(abs(csc_matvec(sys%d, u_leray))) &
            /max(dnorm*maxval(abs(u_leray)), tiny(1.0_dp))
        report%harmonic_orthogonality = 0
        do j = 1, sys%nh
            report%harmonic_orthogonality = max(report%harmonic_orthogonality, &
                abs(dot_product(sys%harm(:, j), csc_matvec(sys%m2, cpsi))) &
                /max(m2norm(sys, cpsi), tiny(1.0_dp)))
        end do
        report%harmonic_part = norm2(hf)
        report%kernel_part = kp
        if (report%jb_scale > 0) report%relative = report%force_potential/report%jb_scale
        if (all(ieee_is_finite([u_potential, u_leray]))) report%status = 0
    end subroutine

    real(dp) function m2norm(sys, v)
        type(nf_system_t), intent(in) :: sys
        real(dp), intent(in) :: v(:)
        m2norm = sqrt(max(dot_product(v, csc_matvec(sys%m2, v)), 0.0_dp))
    end function

    real(dp) function d_scale(sys)
        ! max row sum of |D|.
        type(nf_system_t), intent(in) :: sys
        real(dp), allocatable :: s(:)
        integer :: j, p
        allocate (s(sys%d%nrow)); s = 0
        do j = 1, sys%d%ncol
            do p = sys%d%col_ptr(j), sys%d%col_ptr(j + 1) - 1
                s(sys%d%row_idx(p)) = s(sys%d%row_idx(p)) + abs(sys%d%val(p))
            end do
        end do
        d_scale = maxval(s)
    end function

    subroutine nf_load(sys, fun, load)
        ! Dual V2 load <psi_i, f> of a field f(R, Z, phi) with physical
        ! components (R, phi, Z), nphi-point phi mean.
        type(nf_system_t), intent(in) :: sys
        procedure(nn_vector_field_i) :: fun
        real(dp), allocatable, intent(out) :: load(:)
        real(dp), allocatable :: local(:, :, :)
        real(dp) :: f(3), w
        integer :: t, q, ip, b, nl
        nl = size(sys%basis2, 2)
        allocate (local(nl, 0:sys%nb - 1, sys%nt))
        !$omp parallel do default(shared) private(t, q, ip, b, f, w) schedule(static)
        do t = 1, sys%nt
            local(:, :, t) = 0
            do q = 1, sys%nq
                w = sys%wq(q, t)/sys%nphi
                do ip = 1, sys%nphi
                    call fun([sys%xq(:, q, t), 2*pi*(ip - 1)/sys%nphi], f)
                    do b = 0, sys%nb - 1
                        local(:, b, t) = local(:, b, t) + (w*sys%trig(b, ip)) &
                            *matmul(f, sys%basis2(:, :, q, t))
                    end do
                end do
            end do
        end do
        !$omp end parallel do
        call scatter2(sys, local, load)
    end subroutine

    subroutine scatter2(sys, local, load)
        type(nf_system_t), intent(in) :: sys
        real(dp), intent(in) :: local(:, 0:, :)
        real(dp), allocatable, intent(out) :: load(:)
        integer :: t, b
        allocate (load(sys%nb*sys%n2)); load = 0
        do t = 1, sys%nt
            do b = 0, sys%nb - 1
                load(b*sys%n2 + sys%dof2(:, t)) = load(b*sys%n2 + sys%dof2(:, t)) &
                    + local(:, b, t)*sys%sgn2(:, t)
            end do
        end do
    end subroutine

    function nf_field_coefficients(sys, field) result(b2)
        ! Real Fourier-block V2 coefficients of the NN1 complex amplitudes
        ! (b_n = (c - i s)/2, inverse of relax_store).
        type(nf_system_t), intent(in) :: sys
        type(free_b_field_t), intent(in) :: field
        real(dp), allocatable :: b2(:)
        integer :: n
        allocate (b2(sys%nb*sys%n2))
        b2(:sys%n2) = real(field%b(:, 0), dp)
        do n = 1, sys%modes
            b2((2*n - 1)*sys%n2 + 1:2*n*sys%n2) = 2*real(field%b(:, n), dp)
            b2(2*n*sys%n2 + 1:(2*n + 1)*sys%n2) = -2*aimag(field%b(:, n))
        end do
    end function

    subroutine nf_lorentz_load(sys, field, load, report, status, primal)
        ! F = <V2, J x B>, B = field + harmonic. Default (dual) current:
        ! J = M1^-1 (C^T M2 B + ck) in V1_0, the FEEC curl* (as
        ! kin6d_nonnested_relax); its zero tangential trace leaves an O(1)
        ! wall layer of width h when n x B /= 0 at the wall, so P(J x B)
        ! decreases like h^1/2 there. primal: J = C B1 in V2 with B1 the
        ! full-V1 L2 projection of B (no boundary condition), the strong
        ! curl of the 1-form proxy, without that layer.
        ! report%jb_scale = <|J|^2 |B|^2>^1/2, report%b2 = <|B|^2>.
        type(nf_system_t), intent(inout) :: sys
        type(free_b_field_t), intent(in) :: field
        real(dp), allocatable, intent(out) :: load(:)
        type(nf_report_t), intent(inout) :: report
        integer, intent(out) :: status
        logical, intent(in), optional :: primal
        type(csc_t) :: pair
        type(nn_factor_t) :: full
        complex(dp), allocatable :: hload(:)
        real(dp), allocatable :: b2(:), jd(:), j1(:), local(:, :, :), j2(:)
        logical, allocatable :: none(:)
        real(dp) :: jb(3, 0:sys%nb - 1), bb(3, 0:sys%nb - 1), jx(3), bx(3), f(3)
        real(dp) :: w, radius, s_jb, s_b
        integer :: t, q, ip, b, nl1, nl2
        logical :: strong
        status = 1
        if (.not. sys%ready) return
        strong = .false.
        if (present(primal)) strong = primal
        b2 = nf_field_coefficients(sys, field)
        allocate (j1(sys%nb*sys%n1)); j1 = 0
        if (strong) then
            call meridian_fourier_pairing(field%fs, sys%modes, pair, status)
            if (status /= 0) return
            call meridian_load(field%fs, 1, harmonic_unit, hload, status)
            if (status /= 0) return
            jd = csc_matvec(pair, b2)
            jd(:sys%n1) = jd(:sys%n1) + field%f_vac*real(hload, dp)
            allocate (none(sys%n1)); none = .false.
            call nn_factor(field%fs%mass(1), none, full, status)
            if (status /= 0) return
        else
            jd = csc_matvec(sys%ct, csc_matvec(sys%m2, b2)) + sys%ck
        end if
        do b = 0, sys%nb - 1
            if (strong) then
                call nn_factor_solve_real(full, jd(b*sys%n1 + 1:(b + 1)*sys%n1), &
                    j1(b*sys%n1 + 1:(b + 1)*sys%n1), status)
            else
                call nn_factor_solve_real(sys%mass1, jd(b*sys%n1 + 1:(b + 1)*sys%n1), &
                    j1(b*sys%n1 + 1:(b + 1)*sys%n1), status)
            end if
            if (status /= 0) return
            j1(b*sys%n1 + 1:(b + 1)*sys%n1) = j1(b*sys%n1 + 1:(b + 1)*sys%n1) &
                /sys%weight(b)
        end do
        if (strong) then
            call nn_factor_free(full)
            j2 = csc_matvec(sys%c, j1)
        end if
        status = 1
        nl1 = size(sys%basis1, 2); nl2 = size(sys%basis2, 2)
        allocate (local(nl2, 0:sys%nb - 1, sys%nt))
        s_jb = 0; s_b = 0
        !$omp parallel do default(shared) private(t, q, ip, b, jb, bb, jx, bx, f, &
        !$omp& w, radius) reduction(+:s_jb, s_b) schedule(static)
        do t = 1, sys%nt
            local(:, :, t) = 0
            do q = 1, sys%nq
                do b = 0, sys%nb - 1
                    if (strong) then
                        jb(:, b) = matmul(sys%basis2(:, :, q, t), &
                            sys%sgn2(:, t)*j2(b*sys%n2 + sys%dof2(:, t)))
                    else
                        jb(:, b) = matmul(sys%basis1(:, :, q, t), &
                            sys%sgn1(:, t)*j1(b*sys%n1 + sys%dof1(:, t)))
                    end if
                    bb(:, b) = matmul(sys%basis2(:, :, q, t), &
                        sys%sgn2(:, t)*b2(b*sys%n2 + sys%dof2(:, t)))
                end do
                radius = sys%xq(1, q, t)
                if (field%fs%cyl_r0 > 0) radius = field%fs%cyl_r0
                bb(2, 0) = bb(2, 0) + sys%f_vac/radius
                w = sys%wq(q, t)/sys%nphi
                do ip = 1, sys%nphi
                    jx = matmul(jb, sys%trig(:, ip)); bx = matmul(bb, sys%trig(:, ip))
                    ! (R, phi, Z) is right-handed.
                    f = [jx(2)*bx(3) - jx(3)*bx(2), jx(3)*bx(1) - jx(1)*bx(3), &
                        jx(1)*bx(2) - jx(2)*bx(1)]
                    s_jb = s_jb + w*sum(jx**2)*sum(bx**2)
                    s_b = s_b + w*sum(bx**2)
                    do b = 0, sys%nb - 1
                        local(:, b, t) = local(:, b, t) + (w*sys%trig(b, ip)) &
                            *matmul(f, sys%basis2(:, :, q, t))
                    end do
                end do
            end do
        end do
        !$omp end parallel do
        call scatter2(sys, local, load)
        report%jb_scale = sqrt(s_jb); report%b2 = s_b
        if (all(ieee_is_finite(load))) status = 0
    contains
        subroutine harmonic_unit(x, v)
            real(dp), intent(in) :: x(2)
            complex(dp), intent(out) :: v(3)
            v = 0; v(2) = 1/merge(field%fs%cyl_r0, x(1), field%fs%cyl_r0 > 0)
        end subroutine
    end subroutine

    subroutine nf_l2_error(sys, u, fun, error, norm)
        ! |u_h - u|_L2 and |u|_L2 of a V2 coefficient vector against f.
        type(nf_system_t), intent(in) :: sys
        real(dp), intent(in) :: u(:)
        procedure(nn_vector_field_i) :: fun
        real(dp), intent(out) :: error, norm
        real(dp) :: ub(3, 0:sys%nb - 1), ux(3), f(3), w
        integer :: t, q, ip, b
        error = 0; norm = 0
        !$omp parallel do default(shared) private(t, q, ip, b, ub, ux, f, w) &
        !$omp& reduction(+:error, norm) schedule(static)
        do t = 1, sys%nt
            do q = 1, sys%nq
                do b = 0, sys%nb - 1
                    ub(:, b) = matmul(sys%basis2(:, :, q, t), &
                        sys%sgn2(:, t)*u(b*sys%n2 + sys%dof2(:, t)))
                end do
                w = sys%wq(q, t)/sys%nphi
                do ip = 1, sys%nphi
                    ux = matmul(ub, sys%trig(:, ip))
                    call fun([sys%xq(:, q, t), 2*pi*(ip - 1)/sys%nphi], f)
                    error = error + w*sum((ux - f)**2)
                    norm = norm + w*sum(f**2)
                end do
            end do
        end do
        !$omp end parallel do
        error = sqrt(error); norm = sqrt(norm)
    end subroutine
end module kin6d_nonnested_force
