module kin6d_gs_free_boundary
    !! Free-boundary Grad-Shafranov Newton solver on the FEM-BEM coupled
    !! operator (DESIGN sections 1.2 and 2).
    !!
    !! Unknown: nodal flux psi inside the artificial curve Gamma_c. Residual
    !!
    !!     F(psi; I) = A psi - b(psi) - sum_k I_k load_k - extra,
    !!
    !! A = K + S_h the coupled operator of kin6d_gs_coupling, load_k the coil
    !! loads, extra a psi-independent nodal load, and b the plasma source
    !!
    !!     b_i = int f(s, R) 1{plasma} v_i dA,  s = (psi - psi_a)/(psi_b - psi_a),
    !!
    !! with the normalised profile law of kin6d_gs_profile_law (s = 0 on the
    !! axis, s = 1 on the last closed surface, f = mu0 R p' + F F'/R inside
    !! and zero outside). The plasma is the connected set of s < 1 around
    !! the axis. The cut at s = 1 is one more level of the exact cubic
    !! level-cut quadrature of kin6d_gs_nonlinear.
    !!
    !! Flux levels. psi_a is the extremum of the represented flux near the
    !! axis (a maximum when axis_sign = +1). psi_b is a linear functional of
    !! point values of psi, whose weight vector g_b is the nodal vector of
    !! the point evaluations (the stationary points do not move to first
    !! order, so the envelope argument gives their derivative exactly):
    !!   GS_FB_XPOINT       psi_b = psi_a + c (psi_X - psi_a), psi_X the
    !!                      saddle near xpoint_guess; c = cut_fraction, 1 is
    !!                      the separatrix, c < 1 the FB-cut definition;
    !!   GS_FB_LIMITER      psi_b = the largest sigma*psi over the limiter
    !!                      points (sigma = axis_sign), the plasma touching
    !!                      the limiter point of highest flux;
    !!   GS_FB_TARGET_MEAN  psi_b = the mean of psi over the target points
    !!                      (FreeGS parity).
    !!
    !! Jacobian. With phi = psi - psi_b and a = psi_a - psi_b,
    !!
    !!     d b = M (d psi - 1 d psi_b) + L (d psi_a - d psi_b),
    !!     J   = A - M - u_a g_a^T - u_b g_b^T,
    !!
    !! M the cell mass matrix of the source derivative plus the interface
    !! term of the moving cut, sigma int_{Gamma_p} f_edge v_i v_j/|grad psi| dl
    !! (f_edge the law at s = 1; the cut moves with the normal velocity
    !! sigma (d psi - d psi_b)/|grad psi|), u_a = L the derivative of the load
    !! with respect to the axis flux, u_b = -M 1 - L. J is solved with the
    !! sparse direct factorisation of J0 = A - M (indefinite at fixed coil
    !! currents, so no CG) and a rank-two Woodbury update.
    !!
    !! The interface integral is evaluated along the Hamiltonian flow of
    !! phi in reference coordinates, int g det dt, which has no singularity
    !! at tangencies of the cut with the cell. It diverges logarithmically
    !! when the cut passes through a saddle: in GS_FB_XPOINT mode with
    !! cut_fraction = 1 an arc that stalls at the X point is dropped, so a
    !! nonzero edge source then makes the Jacobian inexact in that cell
    !! (DESIGN FB1, deferred); a law that vanishes at s = 1 is exact.
    use, intrinsic :: iso_fortran_env, only: dp => real64, i64 => int64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_axisymmetric_gs_kernel, only: gs_mu0
    use fortfem_sparse_direct, only: sparse_direct_factor_t, &
                         sparse_direct_factor_csc, sparse_direct_factor_transpose_csc, &
                                     sparse_direct_solve_factored, sparse_direct_free
    use fortnum_quadrature, only: gauss_legendre
    use kin6d_gs_direct, only: gs_direct_t, gs_direct_available, gs_direct_factor, &
                               gs_direct_solve, gs_direct_solve_block, gs_direct_free, &
                               GS_DIRECT_INDEFINITE, GS_DIRECT_ORDER_AMD, &
                               GS_DIRECT_ORDER_ND
    use kin6d_sort_keys, only: sort_key_order
    use kin6d_gs_timers, only: gs_phase_tick, gs_profile_value
    use kin6d_grad_shafranov, only: gs_state_t, gs_cell_data, initialize_gs
    use kin6d_gs_coils, only: gs_circuit_set_t
    use kin6d_gs_coupling, only: gs_coupling_t, gs_coupling_build, gs_coupled_apply, &
                                 gs_coupling_coil_load, GS_COUPLING_OK
    use kin6d_gs_nonlinear, only: gs_cubic_load, gs_profile_jacobian_t
    use kin6d_gs_profile_law, only: gs_profile_law_t, gs_law_sources
    use kin6d_triangle_basis, only: triangle_basis_jet, cubic_segment_roots, &
                                    cubic_segment_stationary, &
                                    triangle_bernstein_controls
    use kin6d_curved_triangle, only: triangle_locate, boundary_curve
    implicit none
    private

    integer, parameter, public :: GS_FB_XPOINT = 1
    integer, parameter, public :: GS_FB_LIMITER = 2
    integer, parameter, public :: GS_FB_TARGET_MEAN = 3

    integer, parameter, public :: GS_FB_OK = 0
    integer, parameter, public :: GS_FB_BAD_INPUT = 1
    integer, parameter, public :: GS_FB_NO_AXIS = 2
    integer, parameter, public :: GS_FB_NO_XPOINT = 3
    integer, parameter, public :: GS_FB_BAD_LCFS = 4
    integer, parameter, public :: GS_FB_SOURCE_FAILED = 5
    integer, parameter, public :: GS_FB_FACTOR_FAILED = 6
    integer, parameter, public :: GS_FB_NOT_CONVERGED = 7
    integer, parameter, public :: GS_FB_LINE_SEARCH = 8
    integer, parameter, public :: GS_FB_REFIT_FAILED = 9
    integer, parameter, public :: GS_FB_REFIT_NOT_CONVERGED = 10

    type, public :: gs_fb_config_t
        integer :: mode = GS_FB_XPOINT
        integer :: axis_sign = 1   !< +1: the axis is the flux maximum, -1 minimum
        real(dp) :: axis_guess(2) = 0.0_dp      !< first guess of the magnetic axis
        real(dp) :: xpoint_guess(2) = 0.0_dp    !< first guess of the X point
        real(dp) :: search_radius = 0.5_dp      !< window of the stationary-point search
        real(dp) :: cut_fraction = 1.0_dp       !< psi_b = psi_a + c (psi_X - psi_a)
        real(dp), allocatable :: points(:, :)   !< limiter / target points (2, m)
        real(dp) :: quadrature_tolerance = 1.0e-9_dp
        real(dp) :: jacobian_tolerance = 1.0e-9_dp
        real(dp) :: tolerance = 1.0e-10_dp      !< Newton: ||F|| / ||b + loads||
        !> A cell face whose nodes all satisfy |psi - psi_b| <= face_tolerance
        !! |psi_a - psi_b| lies on the plasma boundary (fitted mesh). 0 disables.
        real(dp) :: face_tolerance = 1.0e-10_dp
        integer :: max_iterations = 30
    end type gs_fb_config_t

    type, public :: gs_fb_system_t
        type(gs_coupling_t) :: coupling
        type(gs_profile_law_t) :: law          !< normalised law, axis set per call
        type(gs_fb_config_t) :: cfg
        real(dp) :: mu0 = gs_mu0
        real(dp), allocatable :: coil_load(:, :)   !< (n, circuits)
        real(dp), allocatable :: extra_load(:)     !< psi-independent nodal load
        ! CSC pattern of J0 = A - M and the maps into it.
        integer, allocatable :: colptr(:), rowind(:)
        integer, allocatable :: slot_k(:), slot_trace(:, :), slot_cell(:, :, :)
        real(dp) :: axis_hint(2) = 0.0_dp, xpoint_hint(2) = 0.0_dp
        integer :: n = 0
        logical :: ready = .false.
    end type gs_fb_system_t

    type, public :: gs_fb_eval_t
        real(dp), allocatable :: residual(:), source(:)
        real(dp), allocatable :: ga(:), gb(:), ua(:), ub(:)
        real(dp), allocatable :: mass(:, :, :)   !< M per cell, interface included
        real(dp), allocatable :: interface_mass(:, :, :)  !< the interface part of M
        logical, allocatable :: plasma_cell(:)
        real(dp) :: psi_axis = 0.0_dp, psi_bound = 0.0_dp
        real(dp) :: axis(2) = 0.0_dp, xpoint(2) = 0.0_dp
        real(dp) :: plasma_current = 0.0_dp      !< int mu0 j dA / mu0 (A)
        real(dp) :: interface_dropped = 0.0_dp   !< arcs that stalled at a saddle
        integer :: face_cells = 0                !< plasma cell faces on the boundary
        logical :: has_jacobian = .false.
    end type gs_fb_eval_t

    type, public :: gs_fb_factor_t
        type(sparse_direct_factor_t) :: lu, lu_t
        !> Threaded MUMPS LDL^T of the symmetric J0 (replaces lu and lu_t).
        type(gs_direct_t) :: direct
        logical :: use_direct = .false.
        real(dp), allocatable :: values(:), z(:, :), zt(:, :)
        real(dp) :: schur(2, 2) = 0.0_dp, schur_t(2, 2) = 0.0_dp
        logical :: ready = .false., have_transpose = .false.
    end type gs_fb_factor_t

    type, public :: gs_fb_history_t
        integer :: iterations = 0
        real(dp), allocatable :: residual(:)     !< ||F|| / scale per iterate
        real(dp), allocatable :: step(:)         !< ||step||_inf per iteration
        real(dp), allocatable :: damping(:)
        real(dp), allocatable :: ip(:)
    end type gs_fb_history_t

    type, public :: gs_fb_refit_config_t
        integer :: max_refits = 4
        real(dp) :: tolerance = 1.0e-3_dp   !< stop when Gamma_p moves < this * h
        integer :: dense_factor = 8         !< rays per chain point (Gamma_p trace)
        integer :: project_iterations = 4
        !> A refit that leaves the distance above contraction * (previous distance)
        !! ends the loop at the discretisation floor (reported as stagnated).
        real(dp) :: contraction = 0.7_dp
    end type gs_fb_refit_config_t

    type, public :: gs_fb_refit_report_t
        integer :: refits = 0               !< mesh rebuilds performed
        logical :: converged = .false.      !< chain within tolerance * h of boundary
        logical :: stagnated = .false.      !< movement stopped contracting (floor)
        real(dp) :: h = 0.0_dp              !< chain spacing of the final mesh
        real(dp) :: psi_bound = 0.0_dp
        !> Distance of the chain from the computed boundary, in units of h, on the
        !! mesh of pass 0, 1, ... (the last entry used is the final mesh).
        real(dp), allocatable :: movement(:)
        integer, allocatable :: newton_iterations(:)
    end type gs_fb_refit_report_t

    type :: gs_fb_grid_t
        integer :: nx = 0, ny = 0
        real(dp) :: origin(2) = 0.0_dp, step(2) = 1.0_dp
        integer, allocatable :: start(:), list(:)
    end type gs_fb_grid_t

    abstract interface
        subroutine gs_fb_extra_cb(state, sys)
            import :: gs_state_t, gs_fb_system_t
            type(gs_state_t), intent(inout) :: state
            type(gs_fb_system_t), intent(inout) :: sys
        end subroutine gs_fb_extra_cb
    end interface

    public :: gs_fb_refit_solve, gs_fb_trace_lcfs, gs_fb_refit_mesh
    public :: gs_fb_chain_movement
    public :: gs_fb_setup, gs_fb_evaluate, gs_fb_jacobian_apply
    public :: gs_fb_factor, gs_fb_solve, gs_fb_solve_block, gs_fb_release, gs_fb_newton
    public :: gs_fb_interface_mass, gs_fb_lcfs_levels
    public :: gs_fb_locate

contains

    ! ------------------------------------------------------------------
    ! Setup: coupling, coil loads, sparsity pattern of J0.
    ! ------------------------------------------------------------------

    subroutine gs_fb_setup(state, sys, law, cfg, set, status, panels_per_edge, &
                           flux_degree)
        type(gs_state_t), intent(in) :: state
        type(gs_fb_system_t), intent(out) :: sys
        type(gs_profile_law_t), intent(in) :: law
        type(gs_fb_config_t), intent(in) :: cfg
        type(gs_circuit_set_t), intent(in), optional :: set
        integer, intent(out) :: status
        integer, intent(in), optional :: panels_per_edge, flux_degree
        integer :: stat
        real(dp) :: t0
        t0 = -1; call gs_phase_tick('', t0)

        status = GS_FB_BAD_INPUT
        if (.not. law%normalized) return
        if (abs(cfg%axis_sign) /= 1) return
        if (cfg%mode /= GS_FB_XPOINT) then
            if (.not. allocated(cfg%points)) return
            if (size(cfg%points, 1) /= 2 .or. size(cfg%points, 2) < 1) return
        end if
        sys%law = law
        sys%cfg = cfg
        sys%axis_hint = cfg%axis_guess
        sys%xpoint_hint = cfg%xpoint_guess
        call gs_coupling_build(state, sys%coupling, stat, panels_per_edge, flux_degree)
        if (stat /= GS_COUPLING_OK) then
            status = GS_FB_BAD_INPUT + 100*stat
            return
        end if
        call gs_phase_tick('setup.coupling_build', t0)
        sys%n = sys%coupling%n_nodes
        allocate (sys%extra_load(sys%n))
        sys%extra_load = 0.0_dp
        if (present(set)) then
            call gs_coupling_coil_load(state, sys%coupling, set, sys%coil_load, stat)
            if (stat /= GS_COUPLING_OK) then
                status = GS_FB_BAD_INPUT + 200*stat
                return
            end if
        else
            allocate (sys%coil_load(sys%n, 0))
        end if
        call gs_phase_tick('setup.coil_load', t0)
        call build_pattern(state, sys)
        call gs_phase_tick('setup.pattern', t0)
        sys%ready = .true.
        status = GS_FB_OK
    end subroutine gs_fb_setup

    integer function find_slot(keys, key) result(slot)
        integer(i64), intent(in) :: keys(:), key
        integer :: lo, hi, mid

        lo = 1
        hi = size(keys)
        slot = 0
        do while (lo <= hi)
            mid = (lo + hi)/2
            if (keys(mid) == key) then
                slot = mid
                return
            else if (keys(mid) < key) then
                lo = mid + 1
            else
                hi = mid - 1
            end if
        end do
    end function find_slot

    subroutine build_pattern(state, sys)
        !! CSC pattern: the CSR pattern of K united with the dense Gamma_c
        !! block, and the slots of every contribution.
        type(gs_state_t), intent(in) :: state
        type(gs_fb_system_t), intent(inout) :: sys
        integer(i64), allocatable :: keys(:), unique(:)
        integer, allocatable :: order(:)
        integer :: n, nk, nt, cnt, i, j, p, nu, e, k, l, ndof, ne, col
        integer(i64) :: previous

        n = sys%n
        nt = sys%coupling%n_trace
        nk = size(sys%coupling%col)
        allocate (keys(nk + nt*nt), order(nk + nt*nt))
        cnt = nk + nt*nt
        !$omp parallel do default(shared) private(i, p) schedule(static)
        do i = 1, n
            do p = sys%coupling%row_ptr(i), sys%coupling%row_ptr(i + 1) - 1
                keys(p) = int(sys%coupling%col(p) - 1, i64)*n + (i - 1)
            end do
        end do
        !$omp end parallel do
        !$omp parallel do default(shared) private(i, j) schedule(static)
        do j = 1, nt
            do i = 1, nt
                keys(nk + (j - 1)*nt + i) = int(sys%coupling%trace_node(j) - 1, i64)*n &
                            + (sys%coupling%trace_node(i) - 1)
            end do
        end do
        !$omp end parallel do
        call sort_key_order(keys, order)
        allocate (unique(cnt))
        nu = 0
        previous = -1
        do i = 1, cnt
            if (keys(order(i)) /= previous) then
                nu = nu + 1
                unique(nu) = keys(order(i))
                previous = unique(nu)
            end if
        end do
        allocate (sys%colptr(n + 1), sys%rowind(nu))
        sys%colptr = 0
        do i = 1, nu
            col = int(unique(i)/n) + 1
            sys%rowind(i) = int(mod(unique(i), int(n, i64))) + 1
            sys%colptr(col + 1) = sys%colptr(col + 1) + 1
        end do
        sys%colptr(1) = 1
        do i = 1, n
            sys%colptr(i + 1) = sys%colptr(i + 1) + sys%colptr(i)
        end do
        allocate (sys%slot_k(nk), sys%slot_trace(nt, nt))
        !$omp parallel do default(shared) private(i, p) schedule(static)
        do i = 1, n
            do p = sys%coupling%row_ptr(i), sys%coupling%row_ptr(i + 1) - 1
                sys%slot_k(p) = find_slot(unique(:nu), &
                                          int(sys%coupling%col(p) - 1, i64)*n + (i - 1))
            end do
        end do
        !$omp end parallel do
        !$omp parallel do default(shared) private(i, j) schedule(static)
        do j = 1, nt
            do i = 1, nt
                sys%slot_trace(i, j) = find_slot(unique(:nu), &
                                            int(sys%coupling%trace_node(j) - 1, i64)*n &
                                                 + (sys%coupling%trace_node(i) - 1))
            end do
        end do
        !$omp end parallel do
        ndof = size(state%triangles, 1)
        ne = size(state%triangles, 2)
        allocate (sys%slot_cell(ndof, ndof, ne))
        !$omp parallel do default(shared) private(e, k, l) schedule(static)
        do e = 1, ne
            do l = 1, ndof
                do k = 1, ndof
                    sys%slot_cell(k, l, e) = find_slot(unique(:nu), &
                                                 int(state%triangles(l, e) - 1, i64)*n &
                                                       + (state%triangles(k, e) - 1))
                end do
            end do
        end do
        !$omp end parallel do
    end subroutine build_pattern

    ! ------------------------------------------------------------------
    ! Point location and reference-space extrema.
    ! ------------------------------------------------------------------

    subroutine gs_fb_locate(state, r, z, element, xi, found, hint)
        !! Reference coordinates of the physical point (r, z): the first cell
        !! that contains it (the search starts at `hint`).
        type(gs_state_t), intent(in) :: state
        real(dp), intent(in) :: r, z
        integer, intent(out) :: element
        real(dp), intent(out) :: xi(2)
        logical, intent(out) :: found
        integer, intent(in), optional :: hint
        integer :: ne, count, e, status
        real(dp) :: local(2)

        found = .false.
        element = 0
        xi = 0.0_dp
        ne = size(state%triangles, 2)
        do count = 1, ne
            e = count
            if (present(hint)) then
                if (hint >= 1 .and. hint <= ne) e = modulo(hint + count - 2, ne) + 1
            end if
            local = matmul(transpose(state%gradient(:, 2:3, e)), &
                           [r, z] - state%points(:, state%triangles(1, e)))
            if (minval([1.0_dp - sum(local), local]) < -0.5_dp) cycle
            if (state%mapped .and. state%curved_cell(e)) then
                call triangle_locate(state%points(:, state%triangles(:, e)), &
                                     [r, z], local, status)
                if (status /= 0) cycle
            end if
            if (minval([1.0_dp - sum(local), local]) < -1.0e-10_dp) cycle
            element = e
            xi = local
            found = .true.
            return
        end do
    end subroutine gs_fb_locate

    subroutine point_shape(state, e, xi, shape, point, status)
        type(gs_state_t), intent(in) :: state
        integer, intent(in) :: e
        real(dp), intent(in) :: xi(2)
        real(dp), intent(out) :: shape(:), point(2)
        integer, intent(out) :: status
        real(dp) :: grads(2, 10), hess(2, 2, 10), det
        integer :: ndof

        ndof = size(state%triangles, 1)
        call gs_cell_data(state, e, xi, point, det, shape(:ndof), grads(:, :ndof), &
                          hess(:, :, :ndof), status)
    end subroutine point_shape

    subroutine cell_maximum(values, sigma, xi, best, found)
        !! Maximum of sigma * psi_hat over the closed reference triangle for
        !! the nodal values of one cell (vertices, edge stationary points,
        !! interior stationary points).
        real(dp), intent(in) :: values(:), sigma
        real(dp), intent(out) :: xi(2), best
        logical, intent(out) :: found
        real(dp), parameter :: verts(2, 3) = reshape([0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, &
                                                      0.0_dp, 1.0_dp], [2, 3])
        real(dp) :: shape(10), grad(2, 10), hess(2, 2, 10), samples(4), crit(2)
        real(dp) :: p(2), q(2), value, h(2, 2), g(2), det, step(2), start(2)
        integer :: ndof, k, j, count, iter, s

        ndof = size(values)
        best = -huge(1.0_dp)
        found = .false.
        xi = 0.0_dp
        do k = 1, 3
            value = sigma*values(k)
            if (value > best) then
                best = value
                xi = verts(:, k)
                found = .true.
            end if
            p = verts(:, k)
            q = verts(:, mod(k, 3) + 1)
            do j = 1, 4
                call triangle_basis_jet(p + (q - p)*real(j - 1, dp)/3.0_dp, &
                                        shape(:ndof), grad(:, :ndof))
                samples(j) = sigma*dot_product(shape(:ndof), values)
            end do
            call cubic_segment_stationary(samples, crit, count)
            do j = 1, count
                call triangle_basis_jet(p + (q - p)*crit(j), shape(:ndof), &
                                        grad(:, :ndof))
                value = sigma*dot_product(shape(:ndof), values)
                if (value > best) then
                    best = value
                    xi = p + (q - p)*crit(j)
                end if
            end do
        end do
        do s = 1, 4
            select case (s)
            case (1)
                start = [1.0_dp/3, 1.0_dp/3]
            case (2)
                start = [0.6_dp, 0.2_dp]
            case (3)
                start = [0.2_dp, 0.6_dp]
            case default
                start = [0.2_dp, 0.2_dp]
            end select
            p = start
            do iter = 1, 40
             call triangle_basis_jet(p, shape(:ndof), grad(:, :ndof), hess(:, :, :ndof))
                g = matmul(grad(:, :ndof), values)
                h = 0.0_dp
                do j = 1, ndof
                    h = h + hess(:, :, j)*values(j)
                end do
                det = h(1, 1)*h(2, 2) - h(1, 2)*h(2, 1)
                if (abs(det) <= tiny(1.0_dp)) exit
                step = -[h(2, 2)*g(1) - h(1, 2)*g(2), &
                         -h(2, 1)*g(1) + h(1, 1)*g(2)]/det
                p = p + step
                if (maxval(abs(p)) > 4.0_dp) exit
                if (maxval(abs(step)) < 1.0e-15_dp) exit
            end do
            if (maxval(abs(p)) > 4.0_dp) cycle
            if (min(1.0_dp - sum(p), p(1), p(2)) < -1.0e-12_dp) cycle
            call triangle_basis_jet(p, shape(:ndof), grad(:, :ndof), hess(:, :, :ndof))
            g = matmul(grad(:, :ndof), values)
            if (maxval(abs(g)) > 1.0e-9_dp*max(1.0_dp, maxval(abs(values)))) cycle
            h = 0.0_dp
            do j = 1, ndof
                h = h + hess(:, :, j)*values(j)
            end do
            if (sigma*(h(1, 1) + h(2, 2)) >= 0.0_dp) cycle
            if (sigma*sigma*(h(1, 1)*h(2, 2) - h(1, 2)*h(2, 1)) <= 0.0_dp) cycle
            value = sigma*dot_product(shape(:ndof), values)
            if (value > best) then
                best = value
                xi = p
            end if
        end do
    end subroutine cell_maximum

    subroutine cell_saddle(values, xi, found)
        !! Interior saddle (negative Hessian determinant) of the cell
        !! polynomial nearest to the centroid.
        real(dp), intent(in) :: values(:)
        real(dp), intent(out) :: xi(2)
        logical, intent(out) :: found
        real(dp) :: shape(10), grad(2, 10), hess(2, 2, 10), p(2), g(2), h(2, 2)
        real(dp) :: det, step(2), start(2)
        integer :: ndof, iter, j, s

        ndof = size(values)
        found = .false.
        xi = 0.0_dp
        do s = 1, 4
            select case (s)
            case (1)
                start = [1.0_dp/3, 1.0_dp/3]
            case (2)
                start = [0.6_dp, 0.2_dp]
            case (3)
                start = [0.2_dp, 0.6_dp]
            case default
                start = [0.2_dp, 0.2_dp]
            end select
            p = start
            do iter = 1, 40
             call triangle_basis_jet(p, shape(:ndof), grad(:, :ndof), hess(:, :, :ndof))
                g = matmul(grad(:, :ndof), values)
                h = 0.0_dp
                do j = 1, ndof
                    h = h + hess(:, :, j)*values(j)
                end do
                det = h(1, 1)*h(2, 2) - h(1, 2)*h(2, 1)
                if (abs(det) <= tiny(1.0_dp)) exit
                step = -[h(2, 2)*g(1) - h(1, 2)*g(2), &
                         -h(2, 1)*g(1) + h(1, 1)*g(2)]/det
                p = p + step
                if (maxval(abs(p)) > 4.0_dp) exit
                if (maxval(abs(step)) < 1.0e-15_dp) exit
            end do
            if (maxval(abs(p)) > 4.0_dp) cycle
            if (min(1.0_dp - sum(p), p(1), p(2)) < -1.0e-12_dp) cycle
            call triangle_basis_jet(p, shape(:ndof), grad(:, :ndof), hess(:, :, :ndof))
            g = matmul(grad(:, :ndof), values)
            if (maxval(abs(g)) > 1.0e-9_dp*max(1.0_dp, maxval(abs(values)))) cycle
            h = 0.0_dp
            do j = 1, ndof
                h = h + hess(:, :, j)*values(j)
            end do
            if (h(1, 1)*h(2, 2) - h(1, 2)*h(2, 1) >= 0.0_dp) cycle
            xi = p
            found = .true.
            return
        end do
    end subroutine cell_saddle

    ! ------------------------------------------------------------------
    ! Boundary flux levels and their nodal gradients.
    ! ------------------------------------------------------------------

    subroutine gs_fb_lcfs_levels(state, psi, cfg, axis_hint, xpoint_hint, &
                                 psi_axis, psi_bound, axis, xpoint, axis_cell, ga, gb, &
                                 status)
        !! psi_a, psi_b, the stationary points and the nodal vectors g_a,
        !! g_b with d psi_a = g_a . d psi and d psi_b = g_b . d psi.
        type(gs_state_t), intent(inout) :: state
        real(dp), intent(in) :: psi(:)
        type(gs_fb_config_t), intent(in) :: cfg
        real(dp), intent(in) :: axis_hint(2), xpoint_hint(2)
        real(dp), intent(out) :: psi_axis, psi_bound, axis(2), xpoint(2)
        integer, intent(out) :: axis_cell
        real(dp), allocatable, intent(out) :: ga(:), gb(:)
        integer, intent(out) :: status
        real(dp) :: shape(10), values(10), xi(2), value, best, point(2), dist, best_dist
        real(dp) :: psi_x, level, local_best, local_axis(2), local_dist, local_point(2)
        real(dp) :: local_shape(10), local_psi
        integer :: local_cell
        real(dp), allocatable :: shape_x(:)
        integer :: e, ndof, n, k, cell_best, stat, c, element
        logical :: found, any_found, near

        status = GS_FB_NO_AXIS
        ndof = size(state%triangles, 1)
        n = size(psi)
        allocate (ga(n), gb(n))
        ga = 0.0_dp
        gb = 0.0_dp
        axis = axis_hint
        xpoint = xpoint_hint
        psi_axis = 0.0_dp
        psi_bound = 0.0_dp
        axis_cell = 0
        best = -huge(1.0_dp)
        cell_best = 0
        ! Threaded scan; the lowest cell index wins a tie, as in a serial scan.
        !$omp parallel default(shared) private(e, k, near, values, xi, value, found, &
        !$omp& local_best, local_cell, local_axis)
        local_best = -huge(1.0_dp)
        local_cell = 0
        local_axis = 0.0_dp
        !$omp do schedule(static)
        do e = 1, size(state%triangles, 2)
            near = .false.
            do k = 1, ndof
                if (norm2(state%points(:, state%triangles(k, e)) - axis_hint) &
                    <= cfg%search_radius) near = .true.
            end do
            if (.not. near) cycle
            values(:ndof) = psi(state%triangles(:, e))
            call cell_maximum(values(:ndof), real(cfg%axis_sign, dp), xi, value, found)
            if (.not. found) cycle
            if (value > local_best) then
                local_best = value
                local_cell = e
                local_axis = xi
            end if
        end do
        !$omp end do nowait
        !$omp critical(gs_fb_axis_scan)
        if (local_cell /= 0) then
            if (local_best > best .or. (local_best == best .and. local_cell < cell_best)) &
                then
                best = local_best
                cell_best = local_cell
                axis = local_axis
            end if
        end if
        !$omp end critical(gs_fb_axis_scan)
        !$omp end parallel
        if (cell_best == 0) return
        call point_shape(state, cell_best, axis, shape, point, stat)
        if (stat /= 0) return
        axis_cell = cell_best
        axis = point
        psi_axis = real(cfg%axis_sign, dp)*best
        ga(state%triangles(:, cell_best)) = shape(:ndof)

        select case (cfg%mode)
        case (GS_FB_XPOINT)
            status = GS_FB_NO_XPOINT
            best_dist = huge(1.0_dp)
            any_found = .false.
            element = 0
            allocate (shape_x(ndof))
            !$omp parallel default(shared) private(e, k, near, values, xi, found, &
            !$omp& shape, point, stat, dist, local_dist, local_cell, local_point, &
            !$omp& local_shape, local_psi)
            local_dist = huge(1.0_dp)
            local_cell = 0
            local_point = 0.0_dp
            local_shape = 0.0_dp
            local_psi = 0.0_dp
            !$omp do schedule(static)
            do e = 1, size(state%triangles, 2)
                near = .false.
                do k = 1, ndof
                    if (norm2(state%points(:, state%triangles(k, e)) - xpoint_hint) &
                        <= cfg%search_radius) near = .true.
                end do
                if (.not. near) cycle
                values(:ndof) = psi(state%triangles(:, e))
                call cell_saddle(values(:ndof), xi, found)
                if (.not. found) cycle
                call point_shape(state, e, xi, shape, point, stat)
                if (stat /= 0) cycle
                dist = norm2(point - xpoint_hint)
                if (dist < local_dist) then
                    local_dist = dist
                    local_cell = e
                    local_point = point
                    local_shape(:ndof) = shape(:ndof)
                    local_psi = dot_product(shape(:ndof), values(:ndof))
                end if
            end do
            !$omp end do nowait
            !$omp critical(gs_fb_saddle_scan)
            if (local_cell /= 0) then
                if (local_dist < best_dist .or. &
                    (local_dist == best_dist .and. local_cell < element)) then
                    best_dist = local_dist
                    any_found = .true.
                    element = local_cell
                    xpoint = local_point
                    shape_x = local_shape(:ndof)
                    psi_x = local_psi
                end if
            end if
            !$omp end critical(gs_fb_saddle_scan)
            !$omp end parallel
            if (.not. any_found) return
            level = cfg%cut_fraction
            psi_bound = psi_axis + level*(psi_x - psi_axis)
            gb = (1.0_dp - level)*ga
            gb(state%triangles(:, element)) = gb(state%triangles(:, element)) &
                                              + level*shape_x
        case (GS_FB_LIMITER)
            status = GS_FB_BAD_LCFS
            best = -huge(1.0_dp)
            do c = 1, size(cfg%points, 2)
                call gs_fb_locate(state, cfg%points(1, c), cfg%points(2, c), &
                                  element, xi, found)
                if (.not. found) return
                call point_shape(state, element, xi, shape, point, stat)
                if (stat /= 0) return
                value = real(cfg%axis_sign, dp)*dot_product(shape(:ndof), &
                                                       psi(state%triangles(:, element)))
                if (value > best) then
                    best = value
                    cell_best = element
                    gb = 0.0_dp
                    gb(state%triangles(:, element)) = shape(:ndof)
                    xpoint = cfg%points(:, c)
                end if
            end do
            psi_bound = real(cfg%axis_sign, dp)*best
        case (GS_FB_TARGET_MEAN)
            status = GS_FB_BAD_LCFS
            do c = 1, size(cfg%points, 2)
                call gs_fb_locate(state, cfg%points(1, c), cfg%points(2, c), &
                                  element, xi, found)
                if (.not. found) return
                call point_shape(state, element, xi, shape, point, stat)
                if (stat /= 0) return
                gb(state%triangles(:, element)) = gb(state%triangles(:, element)) &
                                                  + shape(:ndof)
            end do
            gb = gb/real(size(cfg%points, 2), dp)
            psi_bound = dot_product(gb, psi)
        case default
            status = GS_FB_BAD_INPUT
            return
        end select
        if (real(cfg%axis_sign, dp)*(psi_axis - psi_bound) <= 0.0_dp) then
            status = GS_FB_BAD_LCFS
            return
        end if
        status = GS_FB_OK
    end subroutine gs_fb_lcfs_levels

  subroutine plasma_cells(state, psi, sigma, psi_bound, axis_cell, mode, cut_fraction, &
                            axis, xpoint, flags)
        !! Cells reached from the axis cell through cells that contain flux
        !! beyond the boundary level on the axis side. For the separatrix the
        !! half plane behind the X point (seen from the axis) is excluded:
        !! it holds the private flux region.
        type(gs_state_t), intent(in) :: state
      real(dp), intent(in) :: psi(:), sigma, psi_bound, cut_fraction, axis(2), xpoint(2)
        integer, intent(in) :: axis_cell, mode
        logical, allocatable, intent(out) :: flags(:)
        integer, allocatable :: queue(:)
        real(dp) :: controls(10), centroid(2)
        integer :: head, tail, e, k, nb, ndof, ne
        logical, allocatable :: seen(:)

        ndof = size(state%triangles, 1)
        ne = size(state%triangles, 2)
        allocate (flags(ne), seen(ne), queue(ne))
        flags = .false.
        seen = .false.
        head = 1
        tail = 1
        queue(1) = axis_cell
        seen(axis_cell) = .true.
        do while (head <= tail)
            e = queue(head)
            head = head + 1
            flags(e) = .true.
            do k = 1, 3
                nb = state%neighbors(k, e)
                if (nb <= 0) cycle
                if (seen(nb)) cycle
                call triangle_bernstein_controls(psi(state%triangles(:, nb)), &
                                                 controls(:ndof))
                if (maxval(sigma*controls(:ndof)) <= sigma*psi_bound) cycle
                if (mode == GS_FB_XPOINT .and. cut_fraction >= 1.0_dp) then
                    centroid = sum(state%points(:, state%triangles(1:3, nb)), 2)/3.0_dp
                    if (dot_product(centroid - xpoint, axis - xpoint) < 0.0_dp) cycle
                end if
                seen(nb) = .true.
                tail = tail + 1
                queue(tail) = nb
            end do
        end do
    end subroutine plasma_cells

    ! ------------------------------------------------------------------
    ! Interface term of the moving cut.
    ! ------------------------------------------------------------------

    subroutine gs_fb_interface_mass(cell_points, affine, phi, mu0, pp_edge, ff_edge, &
                                    sigma, mass, dropped, status)
        !! sigma int_{phi = 0 in the cell} f_edge v_i v_j / |grad psi| dl, with
        !! f_edge(R) = ff_edge/R + mu0 pp_edge R, by the flow of phi in
        !! reference coordinates: int g det dt. `cell_points` are the nodal
        !! positions (isoparametric map; `affine` selects the vertex-only
        !! straight map), `phi` the nodal flux measured from the boundary level.
        real(dp), intent(in) :: cell_points(:, :), phi(:), mu0, pp_edge, ff_edge, sigma
        logical, intent(in) :: affine
        real(dp), intent(out) :: mass(:, :), dropped
        integer, intent(out) :: status
        real(dp), parameter :: verts(2, 3) = reshape([0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, &
                                                      0.0_dp, 1.0_dp], [2, 3])
      real(dp), parameter :: inward(2, 3) = reshape([0.0_dp, 1.0_dp, -1.0_dp, -1.0_dp, &
                                                       1.0_dp, 0.0_dp], [2, 3])
        real(dp) :: samples(4), roots(3), p(2), q(2), t, v(2), yacc(102)
        real(dp) :: shape(10), grad(2, 10), vmax
        integer :: ndof, k, j, nr, nstate
        logical :: stalled

        ndof = size(phi)
        nstate = 2 + ndof*ndof
        mass = 0.0_dp
        dropped = 0.0_dp
        status = 0
        call grad_scale(vmax)
        do k = 1, 3
            p = verts(:, k)
            q = verts(:, mod(k, 3) + 1)
            do j = 1, 4
                call triangle_basis_jet(p + (q - p)*real(j - 1, dp)/3.0_dp, &
                                        shape(:ndof), grad(:, :ndof))
                samples(j) = dot_product(shape(:ndof), phi)
            end do
            call cubic_segment_roots(samples, 0.0_dp, roots, nr)
            do j = 1, nr
                t = roots(j)
                if (t < 0.0_dp .or. t >= 1.0_dp) cycle
                y_start: block
                    real(dp) :: y0(102)
                    y0 = 0.0_dp
                    y0(1:2) = p + (q - p)*t
                    call triangle_basis_jet(y0(1:2), shape(:ndof), grad(:, :ndof))
                    v = matmul(grad(:, :ndof), phi)
                    v = [-v(2), v(1)]
                    if (dot_product(v, inward(:, k)) <= 0.0_dp) exit y_start
                    call trace_arc(y0, yacc, stalled)
                    if (stalled) then
                        dropped = dropped + 1.0_dp
                    else
                        mass = mass + reshape(yacc(3:nstate), [ndof, ndof])
                    end if
                end block y_start
            end do
        end do
        if (.not. all(ieee_is_finite(mass))) status = 1
        mass = sigma*mass
    contains
        subroutine grad_scale(gmax)
            real(dp), intent(out) :: gmax
            real(dp) :: gg(2)
            integer :: a
            gmax = 0.0_dp
            do a = 1, 3
                call triangle_basis_jet(verts(:, a), shape(:ndof), grad(:, :ndof))
                gg = matmul(grad(:, :ndof), phi)
                gmax = max(gmax, norm2(gg))
            end do
            call triangle_basis_jet([1.0_dp/3, 1.0_dp/3], shape(:ndof), grad(:, :ndof))
            gg = matmul(grad(:, :ndof), phi)
            gmax = max(gmax, norm2(gg), tiny(1.0_dp))
        end subroutine grad_scale

        subroutine rhs(y, dy, speed, ok)
            real(dp), intent(in) :: y(:)
            real(dp), intent(out) :: dy(:)
            real(dp), intent(out) :: speed
            logical, intent(out) :: ok
            real(dp) :: n(10), gr(2, 10), gv(2), r_phys, det, jac(2, 2), weight
            integer :: a, b

            call triangle_basis_jet(y(1:2), n(:ndof), gr(:, :ndof))
            gv = matmul(gr(:, :ndof), phi)
            speed = norm2(gv)
            dy = 0.0_dp
            dy(1) = -gv(2)
            dy(2) = gv(1)
            ok = .false.
            if (affine) then
                r_phys = dot_product([1.0_dp - y(1) - y(2), y(1), y(2)], &
                                     cell_points(1, 1:3))
                jac = reshape([cell_points(1, 2) - cell_points(1, 1), &
                               cell_points(2, 2) - cell_points(2, 1), &
                               cell_points(1, 3) - cell_points(1, 1), &
                               cell_points(2, 3) - cell_points(2, 1)], [2, 2])
            else
                r_phys = dot_product(n(:ndof), cell_points(1, :ndof))
                jac = matmul(cell_points(:, :ndof), transpose(gr(:, :ndof)))
            end if
            det = jac(1, 1)*jac(2, 2) - jac(1, 2)*jac(2, 1)
            if (r_phys <= 0.0_dp .or. det <= 0.0_dp) return
            weight = ff_edge/r_phys + mu0*pp_edge*r_phys
            weight = weight*det
            do b = 1, ndof
                do a = 1, ndof
                    dy(2 + (b - 1)*ndof + a) = weight*n(a)*n(b)
                end do
            end do
            ok = .true.
        end subroutine rhs

        subroutine trace_arc(y0, y_end, stalled_out)
            !! Dormand-Prince 5(4) along the flow from an entry point to the
            !! exit point of the reference triangle.
            real(dp), intent(in) :: y0(:)
            real(dp), intent(out) :: y_end(:)
            logical, intent(out) :: stalled_out
            real(dp), parameter :: rtol = 1.0e-12_dp
            real(dp) :: y(102), ynew(102), err(102), kk(102, 7), h, hlo, hhi, hm, sp
            real(dp) :: atol(102), enorm, edge_min, ymid(102), speed_min, qscale
            integer :: steps, i, bis
            logical :: ok

            stalled_out = .false.
            y = y0
            call rhs(y, kk(:, 1), sp, ok)
            if (.not. ok) then
                stalled_out = .true.
                y_end = 0.0_dp
                return
            end if
            speed_min = 1.0e-7_dp*vmax
            if (sp < speed_min) then
                stalled_out = .true.
                y_end = 0.0_dp
                return
            end if
            qscale = maxval(abs(kk(3:nstate, 1)))/sp
            atol = 1.0e-14_dp*max(qscale, tiny(1.0_dp))
            atol(1:2) = 1.0e-14_dp
            h = 0.02_dp/sp
            do steps = 1, 20000
                call dp_step(y, h, ynew, err, kk, ok)
                if (.not. ok) then
                    h = 0.5_dp*h
                    if (h < 1.0e-14_dp/vmax) then
                        stalled_out = .true.
                        exit
                    end if
                    cycle
                end if
                enorm = 0.0_dp
                do i = 1, nstate
                    enorm = max(enorm, abs(err(i))/(atol(i) + rtol*abs(ynew(i))))
                end do
                if (enorm > 1.0_dp) then
                    h = h*max(0.2_dp, 0.9_dp*enorm**(-0.2_dp))
                    cycle
                end if
                edge_min = min(ynew(1), ynew(2), 1.0_dp - ynew(1) - ynew(2))
                if (edge_min < 0.0_dp) then
                    ! The step left the triangle: find the exit by bisection
                    ! on the step length, then finish with the partial step.
                    hlo = 0.0_dp
                    hhi = h
                    do bis = 1, 70
                        hm = 0.5_dp*(hlo + hhi)
                        call dp_step(y, hm, ymid, err, kk, ok)
                        if (ok .and. min(ymid(1), ymid(2), 1.0_dp - ymid(1) - ymid(2)) &
                            >= 0.0_dp) then
                            hlo = hm
                        else
                            hhi = hm
                        end if
                    end do
                    call dp_step(y, hlo, ymid, err, kk, ok)
                    y_end = ymid
                    return
                end if
                y = ynew
                call rhs(y, kk(:, 1), sp, ok)
                if (.not. ok .or. sp < speed_min) then
                    stalled_out = .true.
                    exit
                end if
              h = h*min(5.0_dp, max(0.2_dp, 0.9_dp*(max(enorm, 1.0e-10_dp))**(-0.2_dp)))
                h = min(h, 0.2_dp/sp)
            end do
            stalled_out = .true.
            y_end = y
        end subroutine trace_arc

        subroutine dp_step(y, h, ynew, err, kk, ok)
            real(dp), intent(in) :: y(:), h
            real(dp), intent(out) :: ynew(:), err(:)
            real(dp), intent(inout) :: kk(:, :)
            logical, intent(out) :: ok
            real(dp) :: yt(102), sp
            logical :: ok_stage
            real(dp), parameter :: a21 = 1.0_dp/5
            real(dp), parameter :: a31 = 3.0_dp/40, a32 = 9.0_dp/40
            real(dp), parameter :: a41 = 44.0_dp/45, a42 = -56.0_dp/15, a43 = 32.0_dp/9
            real(dp), parameter :: a51 = 19372.0_dp/6561, a52 = -25360.0_dp/2187, &
                                   a53 = 64448.0_dp/6561, a54 = -212.0_dp/729
            real(dp), parameter :: a61 = 9017.0_dp/3168, a62 = -355.0_dp/33, &
                        a63 = 46732.0_dp/5247, a64 = 49.0_dp/176, a65 = -5103.0_dp/18656
            real(dp), parameter :: a71 = 35.0_dp/384, a73 = 500.0_dp/1113, &
                             a74 = 125.0_dp/192, a75 = -2187.0_dp/6784, a76 = 11.0_dp/84
            real(dp), parameter :: e1 = 71.0_dp/57600, e3 = -71.0_dp/16695, &
                         e4 = 71.0_dp/1920, e5 = -17253.0_dp/339200, e6 = 22.0_dp/525, &
                                   e7 = -1.0_dp/40
            integer :: ns

            ns = nstate
            ok = .false.
            ! Stage 1 is the derivative at y, supplied or recomputed here.
            call rhs(y, kk(:, 1), sp, ok_stage)
            if (.not. ok_stage) return
            yt(:ns) = y(:ns) + h*a21*kk(:ns, 1)
            call rhs(yt, kk(:, 2), sp, ok_stage)
            if (.not. ok_stage) return
            yt(:ns) = y(:ns) + h*(a31*kk(:ns, 1) + a32*kk(:ns, 2))
            call rhs(yt, kk(:, 3), sp, ok_stage)
            if (.not. ok_stage) return
            yt(:ns) = y(:ns) + h*(a41*kk(:ns, 1) + a42*kk(:ns, 2) + a43*kk(:ns, 3))
            call rhs(yt, kk(:, 4), sp, ok_stage)
            if (.not. ok_stage) return
            yt(:ns) = y(:ns) + h*(a51*kk(:ns, 1) + a52*kk(:ns, 2) &
                                  + a53*kk(:ns, 3) + a54*kk(:ns, 4))
            call rhs(yt, kk(:, 5), sp, ok_stage)
            if (.not. ok_stage) return
            yt(:ns) = y(:ns) + h*(a61*kk(:ns, 1) + a62*kk(:ns, 2) &
                                  + a63*kk(:ns, 3) + a64*kk(:ns, 4) &
                                  + a65*kk(:ns, 5))
            call rhs(yt, kk(:, 6), sp, ok_stage)
            if (.not. ok_stage) return
            ynew(:ns) = y(:ns) + h*(a71*kk(:ns, 1) + a73*kk(:ns, 3) &
                                    + a74*kk(:ns, 4) + a75*kk(:ns, 5) &
                                    + a76*kk(:ns, 6))
            call rhs(ynew, kk(:, 7), sp, ok_stage)
            if (.not. ok_stage) return
            err(:ns) = h*(e1*kk(:ns, 1) + e3*kk(:ns, 3) + e4*kk(:ns, 4) &
                          + e5*kk(:ns, 5) + e6*kk(:ns, 6) + e7*kk(:ns, 7))
            ok = .true.
        end subroutine dp_step
    end subroutine gs_fb_interface_mass

    ! ------------------------------------------------------------------
    ! Residual and Jacobian pieces.
    ! ------------------------------------------------------------------

    ! ------------------------------------------------------------------
    ! Faces that coincide with the plasma boundary (fitted mesh).
    ! ------------------------------------------------------------------

    subroutine edge_nodes(order, k, nodes)
        !! Local node numbers of edge k (vertex k, vertex k+1, interior nodes).
        integer, intent(in) :: order, k
        integer, intent(out) :: nodes(order + 1)
        integer :: j

        nodes(1) = k
        nodes(2) = mod(k, 3) + 1
        do j = 1, order - 1
            nodes(2 + j) = 3 + (k - 1)*(order - 1) + j
        end do
    end subroutine edge_nodes

    subroutine snap_faces(state, psi, psi_bound, psi_axis, sigma, tol, plasma, &
                          cell_phi, face_edge)
        !! cell_phi(:, e) = psi - psi_b at the nodes of cell e. A face of a
        !! plasma-flagged cell whose nodes all lie within tol |psi_a - psi_b| of
        !! the boundary level is snapped to exactly zero, which removes the
        !! sliver cuts that rounding noise would create on both sides of it, and
        !! the cell is assigned by the sign of the flux at its centroid: a cell
        !! on the plasma side keeps the face as its boundary (face_edge), a cell
        !! on the vacuum side is removed from the plasma. A cell that is zero
        !! at its centroid too counts as vacuum.
        type(gs_state_t), intent(in) :: state
        real(dp), intent(in) :: psi(:), psi_bound, psi_axis, sigma, tol
        logical, intent(inout) :: plasma(:)
        real(dp), intent(out) :: cell_phi(:, :)
        logical, intent(out) :: face_edge(:, :)
        real(dp) :: shp(10), grd(2, 10), phic, limit
        integer :: e, k, ndof, order, nodes(4)
        logical :: zero_edge(3)

        ndof = size(state%triangles, 1)
        order = merge(2, 3, ndof == 6)
        face_edge = .false.
        limit = tol*abs(psi_axis - psi_bound)
        call triangle_basis_jet([1.0_dp/3.0_dp, 1.0_dp/3.0_dp], shp(:ndof), &
                                grd(:, :ndof))
        do e = 1, size(state%triangles, 2)
            cell_phi(:, e) = psi(state%triangles(:, e)) - psi_bound
            if (tol <= 0.0_dp) cycle
            if (.not. plasma(e)) cycle
            zero_edge = .false.
            do k = 1, 3
                call edge_nodes(order, k, nodes(:order + 1))
                zero_edge(k) = all(abs(cell_phi(nodes(:order + 1), e)) <= limit)
            end do
            if (.not. any(zero_edge)) cycle
            phic = sigma*dot_product(shp(:ndof), cell_phi(:, e))
            do k = 1, 3
                if (.not. zero_edge(k)) cycle
                call edge_nodes(order, k, nodes(:order + 1))
                cell_phi(nodes(:order + 1), e) = 0.0_dp
            end do
            if (phic > limit) then
                face_edge(:, e) = zero_edge
            else
                plasma(e) = .false.
            end if
        end do
    end subroutine snap_faces

    subroutine face_interface_mass(state, e, k, psi, mu0, pp_edge, ff_edge, sigma, &
                                   mass, status)
        !! sigma int_face f_edge v_i v_j / |grad psi| dl over local edge k of
        !! cell e, with the gradient of this cell: the interface term of a cut
        !! that lies on the face (same expression as gs_fb_interface_mass).
        type(gs_state_t), intent(in) :: state
        integer, intent(in) :: e, k
        real(dp), intent(in) :: psi(:), mu0, pp_edge, ff_edge, sigma
        real(dp), intent(out) :: mass(:, :)
        integer, intent(out) :: status
        real(dp), parameter :: verts(2, 3) = reshape([0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, &
                                                      0.0_dp, 1.0_dp], [2, 3])
        real(dp), parameter :: dt = 1.0e-5_dp
        real(dp) :: x(8), w(8), t, xi(2), point(2), pa(2), pb(2), det, dl, gnorm
        real(dp) :: shape(10), grads(2, 10), hess(2, 2, 10), g(2), fe, local(10)
        integer :: q, ndof, st, a, b

        status = 1
        ndof = size(state%triangles, 1)
        mass = 0.0_dp
        local(:ndof) = psi(state%triangles(:, e))
        a = k
        b = mod(k, 3) + 1
        call gauss_legendre(8, x, w)
        do q = 1, 8
            t = 0.5_dp*(x(q) + 1.0_dp)
            xi = (1.0_dp - t)*verts(:, a) + t*verts(:, b)
            call gs_cell_data(state, e, xi, point, det, shape(:ndof), &
                              grads(:, :ndof), hess(:, :, :ndof), st)
            if (st /= 0) return
            call gs_cell_data(state, e, xi + dt*(verts(:, b) - verts(:, a)), pb, det, &
                              shape(:ndof), grads(:, :ndof), hess(:, :, :ndof), st)
            if (st /= 0) return
            call gs_cell_data(state, e, xi - dt*(verts(:, b) - verts(:, a)), pa, det, &
                              shape(:ndof), grads(:, :ndof), hess(:, :, :ndof), st)
            if (st /= 0) return
            dl = norm2(pb - pa)/(2.0_dp*dt)
            call gs_cell_data(state, e, xi, point, det, shape(:ndof), &
                              grads(:, :ndof), hess(:, :, :ndof), st)
            if (st /= 0) return
            g = matmul(grads(:, :ndof), local(:ndof))
            gnorm = norm2(g)
            if (gnorm <= tiny(1.0_dp)) return
            fe = ff_edge/point(1) + mu0*pp_edge*point(1)
            mass(:ndof, :ndof) = mass(:ndof, :ndof) + 0.5_dp*w(q)*sigma*fe*dl/gnorm &
                                 *spread(shape(:ndof), 2, ndof) &
                                 *spread(shape(:ndof), 1, ndof)
        end do
        status = 0
    end subroutine face_interface_mass

    subroutine gs_fb_evaluate(state, sys, psi, currents, ev, status, jacobian)
        type(gs_state_t), intent(inout) :: state
        type(gs_fb_system_t), intent(in) :: sys
        real(dp), intent(in) :: psi(:), currents(:)
        type(gs_fb_eval_t), intent(inout) :: ev
        integer, intent(out) :: status
        logical, intent(in), optional :: jacobian
        type(gs_profile_law_t) :: law
        type(gs_profile_jacobian_t) :: jac
        real(dp), allocatable :: saved(:), action(:), load(:)
    real(dp) :: pp_edge, ff_edge, controls(10), phi(10), dropped_cell, cell_mass(10, 10)
        real(dp) :: sigma, integral
        real(dp), allocatable :: cell_phi(:, :)
        logical, allocatable :: face_edge(:, :)
        integer :: e, ndof, n, ne, stat, k, i
        integer :: dropped_total
        logical :: need_jacobian
        real(dp) :: t0
        t0 = -1; call gs_phase_tick('', t0)

        status = GS_FB_BAD_INPUT
        need_jacobian = .true.
        if (present(jacobian)) need_jacobian = jacobian
        call reset_eval(ev)
        n = sys%n
        if (size(psi) /= n .or. size(currents) /= size(sys%coil_load, 2)) return
        ndof = size(state%triangles, 1)
        ne = size(state%triangles, 2)
        sigma = real(sys%cfg%axis_sign, dp)
        call gs_fb_lcfs_levels(state, psi, sys%cfg, sys%axis_hint, sys%xpoint_hint, &
                               ev%psi_axis, ev%psi_bound, ev%axis, ev%xpoint, &
                               k, ev%ga, ev%gb, stat)
        if (stat /= GS_FB_OK) then
            status = stat
            return
        end if
        call gs_phase_tick('eval.lcfs_levels', t0)
        call plasma_cells(state, psi, sigma, ev%psi_bound, k, sys%cfg%mode, &
                          sys%cfg%cut_fraction, ev%axis, ev%xpoint, ev%plasma_cell)
        law = sys%law
        law%axis = ev%psi_axis - ev%psi_bound
        ! Per-cell boundary-level flux; faces that coincide with the boundary
        ! are snapped onto it and their cells assigned by the centroid sign.
        allocate (cell_phi(ndof, ne), face_edge(3, ne))
        call snap_faces(state, psi, ev%psi_bound, ev%psi_axis, sigma, &
                        sys%cfg%face_tolerance, ev%plasma_cell, cell_phi, face_edge)
        ev%face_cells = count(face_edge)
        allocate (saved(n), ev%source(n))
        call gs_phase_tick('eval.plasma_snap', t0)
        saved = state%psi
        state%psi = psi - ev%psi_bound
        call gs_cubic_load(state, sys%mu0, law, ev%source, jac, stat, &
                           sys%cfg%quadrature_tolerance, &
                           jacobian_tolerance=sys%cfg%jacobian_tolerance, &
                           load_only=.not. need_jacobian, source_integral=integral, &
                           plasma_cell=ev%plasma_cell, cell_psi=cell_phi)
        state%psi = saved
        call gs_phase_tick('eval.cubic_load', t0)
        if (stat /= 0) then
            status = GS_FB_SOURCE_FAILED
            return
        end if
        ev%plasma_current = integral/sys%mu0
        allocate (action(n), load(n))
        call gs_coupled_apply(sys%coupling, psi, action)
        load = sys%extra_load + ev%source
        if (size(currents) > 0) load = load + matmul(sys%coil_load, currents)
        if (allocated(ev%residual)) deallocate (ev%residual)
        allocate (ev%residual(n))
        ev%residual = action - load
        ev%has_jacobian = need_jacobian
        ev%interface_dropped = 0.0_dp
        if (need_jacobian) then
            call move_alloc(jac%mass, ev%mass)
            allocate (ev%interface_mass(ndof, ndof, ne))
            ev%interface_mass = 0.0_dp
            call gs_law_sources(law, 0.0_dp, pp_edge, ff_edge)
            dropped_total = 0
            if (abs(pp_edge) + abs(ff_edge) > 0.0_dp) then
                !$omp parallel do default(shared) &
                !$omp private(e, controls, phi, cell_mass, dropped_cell, stat, i) &
                !$omp reduction(+:dropped_total) schedule(dynamic, 4)
                do e = 1, ne
                    if (.not. ev%plasma_cell(e)) cycle
                    phi(:ndof) = cell_phi(:ndof, e)
                    call triangle_bernstein_controls(phi(:ndof), controls(:ndof))
                    if (minval(controls(:ndof)) >= 0.0_dp) cycle
                    if (maxval(controls(:ndof)) <= 0.0_dp) cycle
                    call gs_fb_interface_mass(state%points(:, state%triangles(:, e)), &
                             .not. state%curved_cell(e), phi(:ndof), sys%mu0, pp_edge, &
                            ff_edge, sigma, cell_mass(:ndof, :ndof), dropped_cell, stat)
                    if (stat /= 0) cycle
                    ev%mass(:, :, e) = ev%mass(:, :, e) + cell_mass(:ndof, :ndof)
                    ev%interface_mass(:, :, e) = cell_mass(:ndof, :ndof)
                    dropped_total = dropped_total + int(dropped_cell)
                end do
                !$omp end parallel do
                ! Faces on the boundary: the cut is the face itself.
                do e = 1, ne
                    if (.not. any(face_edge(:, e))) cycle
                    do k = 1, 3
                        if (.not. face_edge(k, e)) cycle
                        call face_interface_mass(state, e, k, psi, sys%mu0, pp_edge, &
                                                 ff_edge, sigma, &
                                                 cell_mass(:ndof, :ndof), stat)
                        if (stat /= 0) cycle
                        ev%mass(:, :, e) = ev%mass(:, :, e) + cell_mass(:ndof, :ndof)
                        ev%interface_mass(:, :, e) = ev%interface_mass(:, :, e) &
                                                     + cell_mass(:ndof, :ndof)
                    end do
                end do
            end if
            call gs_phase_tick('eval.interface', t0)
            ev%interface_dropped = real(dropped_total, dp)
            ev%ua = jac%axis_load
            allocate (ev%ub(n))
            ev%ub = -ev%ua
            do e = 1, ne
                do i = 1, ndof
                    ev%ub(state%triangles(i, e)) = ev%ub(state%triangles(i, e)) &
                                                   - sum(ev%mass(i, :, e))
                end do
            end do
        end if
        status = GS_FB_OK
    end subroutine gs_fb_evaluate

    subroutine reset_eval(ev)
        type(gs_fb_eval_t), intent(inout) :: ev

        if (allocated(ev%residual)) deallocate (ev%residual)
        if (allocated(ev%source)) deallocate (ev%source)
        if (allocated(ev%ga)) deallocate (ev%ga)
        if (allocated(ev%gb)) deallocate (ev%gb)
        if (allocated(ev%ua)) deallocate (ev%ua)
        if (allocated(ev%ub)) deallocate (ev%ub)
        if (allocated(ev%mass)) deallocate (ev%mass)
        if (allocated(ev%interface_mass)) deallocate (ev%interface_mass)
        if (allocated(ev%plasma_cell)) deallocate (ev%plasma_cell)
        ev%has_jacobian = .false.
    end subroutine reset_eval

    subroutine gs_fb_jacobian_apply(state, sys, ev, input, output, transpose)
        !! output = J input, or J^T input.
        type(gs_state_t), intent(in) :: state
        type(gs_fb_system_t), intent(in) :: sys
        type(gs_fb_eval_t), intent(in) :: ev
        real(dp), intent(in) :: input(:)
        real(dp), intent(out) :: output(:)
        logical, intent(in), optional :: transpose
        real(dp) :: local(10), product(10)
        integer :: e, i, ndof
        logical :: adjoint

        adjoint = .false.
        if (present(transpose)) adjoint = transpose
        ndof = size(state%triangles, 1)
        ! A and the cell mass matrices are symmetric.
        call gs_coupled_apply(sys%coupling, input, output)
        do e = 1, size(state%triangles, 2)
            local(:ndof) = input(state%triangles(:, e))
            product(:ndof) = matmul(ev%mass(:ndof, :ndof, e), local(:ndof))
            do i = 1, ndof
                output(state%triangles(i, e)) = output(state%triangles(i, e)) &
                                                - product(i)
            end do
        end do
        if (adjoint) then
            output = output - ev%ga*dot_product(ev%ua, input) &
                     - ev%gb*dot_product(ev%ub, input)
        else
            output = output - ev%ua*dot_product(ev%ga, input) &
                     - ev%ub*dot_product(ev%gb, input)
        end if
    end subroutine gs_fb_jacobian_apply

    subroutine gs_fb_factor(state, sys, ev, fac, status)
        !! Sparse factor of J0 = A - M and the two Woodbury columns. J0 is
        !! symmetric (stiffness, coupling, source and interface mass), so the
        !! threaded MUMPS LDL^T serves J and J^T; an asymmetric J0, or a build
        !! without MUMPS, uses the serial SuperLU LU and a transposed LU.
        type(gs_state_t), intent(in) :: state
        type(gs_fb_system_t), intent(in) :: sys
        type(gs_fb_eval_t), intent(in) :: ev
        type(gs_fb_factor_t), intent(inout) :: fac
        integer, intent(out) :: status
        real(dp), allocatable :: rhs2(:), col(:)
        real(dp) :: t0, asym, vmax
        integer :: n, nt, ndof, e, k, l, i, j, p, stat
        logical :: use_mumps, ok
        t0 = -1; call gs_phase_tick('', t0)

        status = GS_FB_BAD_INPUT
        if (.not. ev%has_jacobian) return
        call gs_fb_release(fac)
        n = sys%n
        nt = sys%coupling%n_trace
        ndof = size(state%triangles, 1)
        allocate (fac%values(size(sys%rowind)))
        fac%values = 0.0_dp
        !$omp parallel do default(shared) private(p) schedule(static)
        do p = 1, size(sys%slot_k)
            fac%values(sys%slot_k(p)) = sys%coupling%val(p)
        end do
        !$omp end parallel do
        !$omp parallel do default(shared) private(i, j) schedule(static)
        do j = 1, nt
            do i = 1, nt
                !$omp atomic update
                fac%values(sys%slot_trace(i, j)) = &
                    fac%values(sys%slot_trace(i, j)) + sys%coupling%s_trace(i, j)
            end do
        end do
        !$omp end parallel do
        !$omp parallel do default(shared) private(e, k, l) schedule(static)
        do e = 1, size(state%triangles, 2)
            do l = 1, ndof
                do k = 1, ndof
                    !$omp atomic update
                    fac%values(sys%slot_cell(k, l, e)) = &
                        fac%values(sys%slot_cell(k, l, e)) - ev%mass(k, l, e)
                end do
            end do
        end do
        !$omp end parallel do
        call gs_phase_tick('factor.assemble', t0)
        use_mumps = mumps_enabled()
        if (use_mumps) then
            call asymmetry(n, sys%colptr, sys%rowind, fac%values, asym, vmax)
            use_mumps = asym <= 1.0e-12_dp*vmax
        end if
        allocate (fac%z(n, 2), fac%zt(n, 2), rhs2(n), col(n))
        if (use_mumps) then
            call gs_direct_factor(fac%direct, n, sys%colptr, sys%rowind, fac%values, &
                                  ok, GS_DIRECT_INDEFINITE, direct_ordering())
            call gs_phase_tick('factor.mumps', t0)
            call gs_profile_value('mumps.analysis', fac%direct%analysis_s)
            call gs_profile_value('mumps.factor', fac%direct%factor_s)
            if (.not. ok) then
                call gs_direct_free(fac%direct)
                status = GS_FB_FACTOR_FAILED
                return
            end if
            fac%use_direct = .true.
            block
                real(dp), allocatable :: two(:, :)
                allocate (two(n, 2))
                two(:, 1) = ev%ua
                two(:, 2) = ev%ub
                call gs_direct_solve_block(fac%direct, two, fac%z, ok)
            end block
            if (.not. ok) then
                status = GS_FB_FACTOR_FAILED
                return
            end if
            call gs_phase_tick('factor.z', t0)
            fac%schur(1, 1) = 1.0_dp - dot_product(ev%ga, fac%z(:, 1))
            fac%schur(1, 2) = -dot_product(ev%ga, fac%z(:, 2))
            fac%schur(2, 1) = -dot_product(ev%gb, fac%z(:, 1))
            fac%schur(2, 2) = 1.0_dp - dot_product(ev%gb, fac%z(:, 2))
            fac%ready = .true.
            status = GS_FB_OK
            return
        end if
        call sparse_direct_factor_csc(fac%lu, n, sys%colptr, sys%rowind, &
                                      fac%values, stat)
        call gs_phase_tick('factor.lu', t0)
        if (stat /= 0) then
            status = GS_FB_FACTOR_FAILED
            return
        end if
        rhs2 = ev%ua
        call sparse_direct_solve_factored(fac%lu, rhs2, col, stat)
        if (stat /= 0) then
            status = GS_FB_FACTOR_FAILED
            return
        end if
        fac%z(:, 1) = col
        rhs2 = ev%ub
        call sparse_direct_solve_factored(fac%lu, rhs2, col, stat)
        if (stat /= 0) then
            status = GS_FB_FACTOR_FAILED
            return
        end if
        fac%z(:, 2) = col
        call gs_phase_tick('factor.z', t0)
        fac%schur(1, 1) = 1.0_dp - dot_product(ev%ga, fac%z(:, 1))
        fac%schur(1, 2) = -dot_product(ev%ga, fac%z(:, 2))
        fac%schur(2, 1) = -dot_product(ev%gb, fac%z(:, 1))
        fac%schur(2, 2) = 1.0_dp - dot_product(ev%gb, fac%z(:, 2))
        fac%ready = .true.
        status = GS_FB_OK
    end subroutine gs_fb_factor

    integer function direct_ordering() result(ordering)
        !! Fill-reducing ordering of the MUMPS analysis. The factor is redone
        !! every Newton iteration and the analysis is serial, so approximate
        !! minimum degree (cheap analysis, fill slightly above nested
        !! dissection) is faster overall; KIN6D_FB_ORDERING = nd selects nested
        !! dissection.
        character(len=8) :: choice
        integer :: length, stat

        ordering = GS_DIRECT_ORDER_AMD
        call get_environment_variable('KIN6D_FB_ORDERING', choice, length, stat)
        if (stat == 0 .and. length > 0) then
            if (trim(choice) == 'nd') ordering = GS_DIRECT_ORDER_ND
        end if
    end function direct_ordering

    logical function mumps_enabled() result(enabled)
        !! Threaded MUMPS when built in; KIN6D_FB_SOLVER=superlu forces the
        !! serial SuperLU path (A/B comparisons).
        character(len=16) :: choice
        integer :: length, stat

        enabled = gs_direct_available()
        if (.not. enabled) return
        call get_environment_variable('KIN6D_FB_SOLVER', choice, length, stat)
        if (stat == 0 .and. length > 0) enabled = trim(choice) /= 'superlu'
    end function mumps_enabled

    subroutine asymmetry(n, colptr, rowind, values, asym, vmax)
        !! max |A_ij - A_ji| and max |A_ij| of a CSC matrix with sorted rows.
        integer, intent(in) :: n, colptr(:), rowind(:)
        real(dp), intent(in) :: values(:)
        real(dp), intent(out) :: asym, vmax
        integer :: c, p, r, lo, hi, mid

        asym = 0.0_dp
        vmax = 0.0_dp
        !$omp parallel do default(shared) private(c, p, r, lo, hi, mid) &
        !$omp reduction(max:asym, vmax) schedule(static)
        do c = 1, n
            do p = colptr(c), colptr(c + 1) - 1
                r = rowind(p)
                vmax = max(vmax, abs(values(p)))
                ! Entry (c, r) lives in column r.
                lo = colptr(r)
                hi = colptr(r + 1) - 1
                do while (lo < hi)
                    mid = (lo + hi)/2
                    if (rowind(mid) < c) then
                        lo = mid + 1
                    else
                        hi = mid
                    end if
                end do
                if (lo <= colptr(r + 1) - 1) then
                    if (rowind(lo) == c) then
                        asym = max(asym, abs(values(p) - values(lo)))
                    else
                        asym = max(asym, abs(values(p)))
                    end if
                else
                    asym = max(asym, abs(values(p)))
                end if
            end do
        end do
        !$omp end parallel do
    end subroutine asymmetry

    subroutine gs_fb_solve(sys, ev, fac, rhs, solution, status, transpose)
        !! solution = J^{-1} rhs (or J^{-T} rhs) with the factor of J0 and the
        !! rank-two Woodbury update.
        type(gs_fb_system_t), intent(in) :: sys
        type(gs_fb_eval_t), intent(in) :: ev
        type(gs_fb_factor_t), intent(inout) :: fac
        real(dp), intent(in) :: rhs(:)
        real(dp), intent(out) :: solution(:)
        integer, intent(out) :: status
        logical, intent(in), optional :: transpose
        real(dp), allocatable :: y(:), col(:), rhs_v(:)
        real(dp) :: w(2), c(2, 2), det, d(2)
        integer :: n, stat, j
        logical :: adjoint, ok

        status = GS_FB_BAD_INPUT
        if (.not. fac%ready) return
        adjoint = .false.
        if (present(transpose)) adjoint = transpose
        n = sys%n
        allocate (y(n), col(n), rhs_v(n))
        if (.not. adjoint) then
            if (fac%use_direct) then
                call gs_direct_solve(fac%direct, rhs, y, ok)
                stat = merge(0, 1, ok)
            else
                call sparse_direct_solve_factored(fac%lu, rhs, y, stat)
            end if
            if (stat /= 0) then
                status = GS_FB_FACTOR_FAILED
                return
            end if
            w = [dot_product(ev%ga, y), dot_product(ev%gb, y)]
            c = fac%schur
        else
            if (fac%use_direct .and. .not. fac%have_transpose) then
                ! J0 is symmetric: the same factor solves J0^T.
                block
                    real(dp), allocatable :: two(:, :)
                    allocate (two(n, 2))
                    two(:, 1) = ev%ga
                    two(:, 2) = ev%gb
                    call gs_direct_solve_block(fac%direct, two, fac%zt, ok)
                end block
                if (.not. ok) then
                    status = GS_FB_FACTOR_FAILED
                    return
                end if
                fac%schur_t(1, 1) = 1.0_dp - dot_product(ev%ua, fac%zt(:, 1))
                fac%schur_t(1, 2) = -dot_product(ev%ua, fac%zt(:, 2))
                fac%schur_t(2, 1) = -dot_product(ev%ub, fac%zt(:, 1))
                fac%schur_t(2, 2) = 1.0_dp - dot_product(ev%ub, fac%zt(:, 2))
                fac%have_transpose = .true.
            end if
            if (.not. fac%have_transpose) then
                call sparse_direct_factor_transpose_csc(fac%lu_t, n, sys%colptr, &
                                                        sys%rowind, fac%values, stat)
                if (stat /= 0) then
                    status = GS_FB_FACTOR_FAILED
                    return
                end if
                do j = 1, 2
                    rhs_v = ev%ga
                    if (j == 2) rhs_v = ev%gb
                    call sparse_direct_solve_factored(fac%lu_t, rhs_v, col, stat)
                    if (stat /= 0) then
                        status = GS_FB_FACTOR_FAILED
                        return
                    end if
                    fac%zt(:, j) = col
                end do
                ! J^T = J0^T - V U^T: Schur matrix I - U^T J0^{-T} V.
                fac%schur_t(1, 1) = 1.0_dp - dot_product(ev%ua, fac%zt(:, 1))
                fac%schur_t(1, 2) = -dot_product(ev%ua, fac%zt(:, 2))
                fac%schur_t(2, 1) = -dot_product(ev%ub, fac%zt(:, 1))
                fac%schur_t(2, 2) = 1.0_dp - dot_product(ev%ub, fac%zt(:, 2))
                fac%have_transpose = .true.
            end if
            if (fac%use_direct) then
                call gs_direct_solve(fac%direct, rhs, y, ok)
                stat = merge(0, 1, ok)
            else
                call sparse_direct_solve_factored(fac%lu_t, rhs, y, stat)
            end if
            if (stat /= 0) then
                status = GS_FB_FACTOR_FAILED
                return
            end if
            w = [dot_product(ev%ua, y), dot_product(ev%ub, y)]
            c = fac%schur_t
        end if
        det = c(1, 1)*c(2, 2) - c(1, 2)*c(2, 1)
        if (abs(det) <= 1.0e-14_dp*max(1.0_dp, abs(c(1, 1)*c(2, 2)))) then
            status = GS_FB_FACTOR_FAILED
            return
        end if
        d = [c(2, 2)*w(1) - c(1, 2)*w(2), -c(2, 1)*w(1) + c(1, 1)*w(2)]/det
        if (.not. adjoint) then
            solution = y + fac%z(:, 1)*d(1) + fac%z(:, 2)*d(2)
        else
            solution = y + fac%zt(:, 1)*d(1) + fac%zt(:, 2)*d(2)
        end if
        status = GS_FB_OK
    end subroutine gs_fb_solve

    subroutine gs_fb_solve_block(sys, ev, fac, rhs, solution, status)
        !! solution(:, k) = J^{-1} rhs(:, k) for every column. A MUMPS factor
        !! is solved with one threaded block solve (its handle is not safe to
        !! share between threads); SuperLU solves are threaded over columns.
        type(gs_fb_system_t), intent(in) :: sys
        type(gs_fb_eval_t), intent(in) :: ev
        type(gs_fb_factor_t), intent(inout) :: fac
        real(dp), intent(in) :: rhs(:, :)
        real(dp), intent(out) :: solution(:, :)
        integer, intent(out) :: status
        real(dp), allocatable :: y(:, :)
        real(dp) :: c(2, 2), det, w(2), d(2)
        integer, allocatable :: stats(:)
        integer :: k, nrhs
        logical :: ok

        status = GS_FB_BAD_INPUT
        if (.not. fac%ready) return
        nrhs = size(rhs, 2)
        if (.not. fac%use_direct) then
            allocate (stats(nrhs))
            stats = GS_FB_OK
            !$omp parallel do default(shared) private(k) schedule(dynamic, 1)
            do k = 1, nrhs
                call gs_fb_solve(sys, ev, fac, rhs(:, k), solution(:, k), stats(k))
            end do
            !$omp end parallel do
            status = GS_FB_OK
            if (any(stats /= GS_FB_OK)) status = GS_FB_FACTOR_FAILED
            return
        end if
        allocate (y(sys%n, nrhs))
        call gs_direct_solve_block(fac%direct, rhs, y, ok)
        if (.not. ok) then
            status = GS_FB_FACTOR_FAILED
            return
        end if
        c = fac%schur
        det = c(1, 1)*c(2, 2) - c(1, 2)*c(2, 1)
        if (abs(det) <= 1.0e-14_dp*max(1.0_dp, abs(c(1, 1)*c(2, 2)))) then
            status = GS_FB_FACTOR_FAILED
            return
        end if
        !$omp parallel do default(shared) private(k, w, d) schedule(static)
        do k = 1, nrhs
            w = [dot_product(ev%ga, y(:, k)), dot_product(ev%gb, y(:, k))]
            d = [c(2, 2)*w(1) - c(1, 2)*w(2), -c(2, 1)*w(1) + c(1, 1)*w(2)]/det
            solution(:, k) = y(:, k) + fac%z(:, 1)*d(1) + fac%z(:, 2)*d(2)
        end do
        !$omp end parallel do
        status = GS_FB_OK
    end subroutine gs_fb_solve_block

    subroutine gs_fb_release(fac)
        type(gs_fb_factor_t), intent(inout) :: fac

        if (fac%use_direct) then
            call gs_direct_free(fac%direct)
        else
            if (fac%ready) call sparse_direct_free(fac%lu)
            if (fac%have_transpose) call sparse_direct_free(fac%lu_t)
        end if
        fac%use_direct = .false.
        fac%ready = .false.
        fac%have_transpose = .false.
        if (allocated(fac%values)) deallocate (fac%values)
        if (allocated(fac%z)) deallocate (fac%z)
        if (allocated(fac%zt)) deallocate (fac%zt)
    end subroutine gs_fb_release

    ! ------------------------------------------------------------------
    ! Damped Newton.
    ! ------------------------------------------------------------------

    subroutine gs_fb_newton(state, sys, currents, status, history)
        !! Newton iteration from the flux in state%psi; the converged flux is
        !! returned in state%psi. Armijo backtracking on ||F||.
        type(gs_state_t), intent(inout) :: state
        type(gs_fb_system_t), intent(inout) :: sys
        real(dp), intent(in) :: currents(:)
        integer, intent(out) :: status
        type(gs_fb_history_t), intent(out), optional :: history
        type(gs_fb_eval_t) :: ev, trial
        type(gs_fb_factor_t) :: fac
        real(dp), allocatable :: psi(:), step(:), candidate(:), scale_load(:)
        real(dp), allocatable :: res(:), stp(:), damp(:), ipv(:)
        real(dp) :: old_norm, new_norm, damping, scale, axis_save(2), x_save(2)
        integer :: iteration, backtrack, stat, n

        status = GS_FB_BAD_INPUT
        if (.not. sys%ready) return
        n = sys%n
        allocate (psi(n), step(n), candidate(n), res(0:sys%cfg%max_iterations), &
                  stp(0:sys%cfg%max_iterations), damp(0:sys%cfg%max_iterations), &
                  ipv(0:sys%cfg%max_iterations))
        res = 0.0_dp
        stp = 0.0_dp
        damp = 0.0_dp
        ipv = 0.0_dp
        psi = state%psi
        call gs_fb_evaluate(state, sys, psi, currents, ev, stat)
        if (stat /= GS_FB_OK) then
            status = stat
            return
        end if
        scale = 0.0_dp
        do iteration = 0, sys%cfg%max_iterations
            scale_load = ev%source + sys%extra_load
       if (size(currents) > 0) scale_load = scale_load + matmul(sys%coil_load, currents)
            scale = max(norm2(scale_load), tiny(1.0_dp))
            old_norm = norm2(ev%residual)
            res(iteration) = old_norm/scale
            ipv(iteration) = ev%plasma_current
            sys%axis_hint = ev%axis
            if (cfg_x(sys)) sys%xpoint_hint = ev%xpoint
            if (res(iteration) <= sys%cfg%tolerance) then
                status = GS_FB_OK
                exit
            end if
            if (iteration == sys%cfg%max_iterations) then
                status = GS_FB_NOT_CONVERGED
                exit
            end if
            call gs_fb_factor(state, sys, ev, fac, stat)
            if (stat /= GS_FB_OK) then
                status = stat
                exit
            end if
            call gs_fb_solve(sys, ev, fac, -ev%residual, step, stat)
            if (stat /= GS_FB_OK) then
                status = stat
                exit
            end if
            stp(iteration) = maxval(abs(step))
            damping = 1.0_dp
            axis_save = sys%axis_hint
            x_save = sys%xpoint_hint
            do backtrack = 0, 20
                candidate = psi + damping*step
                call gs_fb_evaluate(state, sys, candidate, currents, trial, stat)
                if (stat == GS_FB_OK) then
                    new_norm = norm2(trial%residual)
                    if (new_norm <= (1.0_dp - 1.0e-4_dp*damping)*old_norm) exit
                end if
                damping = 0.5_dp*damping
                sys%axis_hint = axis_save
                sys%xpoint_hint = x_save
            end do
            if (backtrack > 20) then
                status = GS_FB_LINE_SEARCH
                exit
            end if
            damp(iteration) = damping
            psi = candidate
            call move_ev(trial, ev)
        end do
        state%psi = psi
        call gs_fb_release(fac)
        if (present(history)) then
            history%iterations = min(iteration, sys%cfg%max_iterations)
            allocate (history%residual(0:history%iterations), &
                      history%step(0:history%iterations), &
                      history%damping(0:history%iterations), &
                      history%ip(0:history%iterations))
            history%residual = res(0:history%iterations)
            history%step = stp(0:history%iterations)
            history%damping = damp(0:history%iterations)
            history%ip = ipv(0:history%iterations)
        end if
    contains
        logical function cfg_x(s) result(flag)
            type(gs_fb_system_t), intent(in) :: s
            flag = s%cfg%mode == GS_FB_XPOINT
        end function cfg_x

        subroutine move_ev(from, to)
            type(gs_fb_eval_t), intent(inout) :: from, to
            to%psi_axis = from%psi_axis
            to%psi_bound = from%psi_bound
            to%axis = from%axis
            to%xpoint = from%xpoint
            to%plasma_current = from%plasma_current
            to%interface_dropped = from%interface_dropped
            to%face_cells = from%face_cells
            to%has_jacobian = from%has_jacobian
            call move_alloc(from%residual, to%residual)
            call move_alloc(from%source, to%source)
            call move_alloc(from%ga, to%ga)
            call move_alloc(from%gb, to%gb)
            call move_alloc(from%ua, to%ua)
            call move_alloc(from%ub, to%ub)
            call move_alloc(from%mass, to%mass)
            call move_alloc(from%interface_mass, to%interface_mass)
            call move_alloc(from%plasma_cell, to%plasma_cell)
        end subroutine move_ev
    end subroutine gs_fb_newton

    ! ------------------------------------------------------------------
    ! P8: outer refit loop (DESIGN 1.2, step 2).
    ! ------------------------------------------------------------------

    subroutine cell_contains(state, e, r, z, local, ok)
        !! Reference coordinates of (r, z) in cell e, ok when inside it.
        type(gs_state_t), intent(in) :: state
        integer, intent(in) :: e
        real(dp), intent(in) :: r, z
        real(dp), intent(out) :: local(2)
        logical, intent(out) :: ok
        integer :: status

        ok = .false.
        local = matmul(transpose(state%gradient(:, 2:3, e)), &
                       [r, z] - state%points(:, state%triangles(1, e)))
        if (minval([1.0_dp - sum(local), local]) < -0.5_dp) return
        if (state%mapped .and. state%curved_cell(e)) then
            call triangle_locate(state%points(:, state%triangles(:, e)), [r, z], &
                                 local, status)
            if (status /= 0) return
        end if
        ok = minval([1.0_dp - sum(local), local]) >= -1.0e-10_dp
    end subroutine cell_contains

    subroutine build_grid(state, grid)
        !! Uniform bucket grid over the bounding boxes of the cells (CSR lists).
        type(gs_state_t), intent(in) :: state
        type(gs_fb_grid_t), intent(out) :: grid
        real(dp) :: lo(2), hi(2)
        integer, allocatable :: fill(:)
        integer :: e, ne, i, j, i0, i1, j0, j1, nb, pass, b

        ne = size(state%triangles, 2)
        lo = minval(state%points, dim=2)
        hi = maxval(state%points, dim=2)
        grid%nx = max(1, nint(sqrt(real(ne, dp)*(hi(1) - lo(1))/(hi(2) - lo(2)))))
        grid%ny = max(1, nint(sqrt(real(ne, dp)*(hi(2) - lo(2))/(hi(1) - lo(1)))))
        grid%origin = lo
        grid%step = (hi - lo)/real([grid%nx, grid%ny], dp)
        nb = grid%nx*grid%ny
        allocate (grid%start(nb + 1), fill(nb))
        grid%start = 0
        fill = 0
        do pass = 1, 2
            do e = 1, ne
                call cell_range(e, i0, i1, j0, j1)
                do j = j0, j1
                    do i = i0, i1
                        b = (j - 1)*grid%nx + i
                        if (pass == 1) then
                            grid%start(b + 1) = grid%start(b + 1) + 1
                        else
                            grid%list(grid%start(b) + fill(b)) = e
                            fill(b) = fill(b) + 1
                        end if
                    end do
                end do
            end do
            if (pass == 1) then
                grid%start(1) = 1
                do b = 1, nb
                    grid%start(b + 1) = grid%start(b) + grid%start(b + 1)
                end do
                allocate (grid%list(grid%start(nb + 1) - 1))
            end if
        end do
    contains
        subroutine cell_range(cell, ia, ib, ja, jb)
            integer, intent(in) :: cell
            integer, intent(out) :: ia, ib, ja, jb
            real(dp) :: low(2), high(2), pad

            low = minval(state%points(:, state%triangles(:, cell)), dim=2)
            high = maxval(state%points(:, state%triangles(:, cell)), dim=2)
            pad = 0.05_dp*maxval(high - low)
            low = low - pad
            high = high + pad
            ia = bucket(low(1), 1, grid%nx)
            ib = bucket(high(1), 1, grid%nx)
            ja = bucket(low(2), 2, grid%ny)
            jb = bucket(high(2), 2, grid%ny)
        end subroutine cell_range

        integer function bucket(x, dim, n) result(index)
            real(dp), intent(in) :: x
            integer, intent(in) :: dim, n
            index = min(n, max(1, 1 + int((x - grid%origin(dim))/grid%step(dim))))
        end function bucket
    end subroutine build_grid

    subroutine grid_locate(state, grid, r, z, element, xi, found)
        type(gs_state_t), intent(in) :: state
        type(gs_fb_grid_t), intent(in) :: grid
        real(dp), intent(in) :: r, z
        integer, intent(out) :: element
        real(dp), intent(out) :: xi(2)
        logical, intent(out) :: found
        integer :: i, j, b, k, e
        logical :: ok

        found = .false.
        element = 0
        xi = 0.0_dp
        i = 1 + int((r - grid%origin(1))/grid%step(1))
        j = 1 + int((z - grid%origin(2))/grid%step(2))
        if (i < 1 .or. i > grid%nx .or. j < 1 .or. j > grid%ny) return
        b = (j - 1)*grid%nx + i
        do k = grid%start(b), grid%start(b + 1) - 1
            e = grid%list(k)
            call cell_contains(state, e, r, z, xi, ok)
            if (ok) then
                element = e
                found = .true.
                return
            end if
        end do
    end subroutine grid_locate

    subroutine grid_flux(state, grid, psi, r, z, value, gradient, found)
        !! Value and physical gradient of the finite element flux at (r, z).
        type(gs_state_t), intent(in) :: state
        type(gs_fb_grid_t), intent(in) :: grid
        real(dp), intent(in) :: psi(:), r, z
        real(dp), intent(out) :: value, gradient(2)
        logical, intent(out) :: found
        real(dp) :: xi(2), point(2), det, shape(10), grads(2, 10), hess(2, 2, 10)
        integer :: e, ndof, status

        value = 0.0_dp
        gradient = 0.0_dp
        call grid_locate(state, grid, r, z, e, xi, found)
        if (.not. found) return
        ndof = size(state%triangles, 1)
        call gs_cell_data(state, e, xi, point, det, shape(:ndof), grads(:, :ndof), &
                          hess(:, :, :ndof), status)
        found = status == 0
        if (.not. found) return
        value = dot_product(shape(:ndof), psi(state%triangles(:, e)))
        gradient = matmul(grads(:, :ndof), psi(state%triangles(:, e)))
    end subroutine grid_flux

    subroutine gs_fb_trace_lcfs(state, psi, axis, psi_bound, sigma, n_points, chain, &
                                status, dense_factor, project_iterations)
        !! n_points points on the connected curve {psi = psi_b} around the axis,
        !! equally spaced in arc length of the traced polyline and projected
        !! onto the curve (Newton along grad psi), in counter-clockwise order.
        !! The curve is found along rays from the axis (star-shaped plasma).
        type(gs_state_t), intent(in) :: state
        real(dp), intent(in) :: psi(:), axis(2), psi_bound, sigma
        integer, intent(in) :: n_points
        real(dp), allocatable, intent(out) :: chain(:, :)
        integer, intent(out) :: status
        integer, intent(in), optional :: dense_factor, project_iterations
        type(gs_fb_grid_t) :: grid
        real(dp), allocatable :: dense(:, :), arc(:)
        real(dp) :: step, theta, tau, tau_prev, v, v_prev, g(2), x(2), total, target, f
        real(dp) :: lo, hi, mid
        integer :: m, nd, k, it, project, kk, bis, factor
        logical :: found, failed

        status = GS_FB_BAD_LCFS
        factor = 8
        if (present(dense_factor)) factor = dense_factor
        project = 4
        if (present(project_iterations)) project = project_iterations
        nd = factor*n_points
        call build_grid(state, grid)
        step = 0.5_dp*sqrt(sum(state%area)/size(state%area))
        allocate (dense(2, nd), arc(nd + 1))
        failed = .false.
        !$omp parallel do default(shared) &
        !$omp private(m, theta, tau, tau_prev, v, v_prev, g, found, lo, hi, mid) &
        !$omp private(bis, x) &
        !$omp schedule(static)
        do m = 1, nd
            theta = 2.0_dp*acos(-1.0_dp)*real(m - 1, dp)/nd
            tau_prev = 0.0_dp
            tau = 0.0_dp
            v_prev = 1.0_dp
            found = .false.
            do
                tau = tau + step
                x = axis + tau*[cos(theta), sin(theta)]
                call grid_flux(state, grid, psi, x(1), x(2), v, g, found)
                if (.not. found) then
                    failed = .true.
                    exit
                end if
                v = sigma*(v - psi_bound)
                if (v <= 0.0_dp) exit
                tau_prev = tau
                v_prev = v
            end do
            if (.not. found) cycle
            lo = tau_prev
            hi = tau
            do bis = 1, 60
                mid = 0.5_dp*(lo + hi)
                x = axis + mid*[cos(theta), sin(theta)]
                call grid_flux(state, grid, psi, x(1), x(2), v, g, found)
                if (.not. found) exit
                if (sigma*(v - psi_bound) > 0.0_dp) then
                    lo = mid
                else
                    hi = mid
                end if
            end do
            dense(:, m) = axis + 0.5_dp*(lo + hi)*[cos(theta), sin(theta)]
        end do
        !$omp end parallel do
        if (failed) return
        arc(1) = 0.0_dp
        do m = 1, nd
            k = mod(m, nd) + 1
            arc(m + 1) = arc(m) + norm2(dense(:, k) - dense(:, m))
        end do
        total = arc(nd + 1)
        allocate (chain(2, n_points))
        kk = 1
        do it = 1, n_points
            target = total*real(it - 1, dp)/n_points
            do while (arc(kk + 1) < target .and. kk < nd)
                kk = kk + 1
            end do
            f = (target - arc(kk))/max(arc(kk + 1) - arc(kk), tiny(1.0_dp))
            chain(:, it) = (1.0_dp - f)*dense(:, kk) + f*dense(:, mod(kk, nd) + 1)
        end do
        !$omp parallel do default(shared) private(it, k, v, g, found) schedule(static)
        do it = 1, n_points
            do k = 1, project
                call grid_flux(state, grid, psi, chain(1, it), chain(2, it), v, g, &
                               found)
                if (.not. found .or. dot_product(g, g) <= tiny(1.0_dp)) exit
                chain(:, it) = chain(:, it) - (v - psi_bound)*g/dot_product(g, g)
            end do
        end do
        !$omp end parallel do
        status = GS_FB_OK
    end subroutine gs_fb_trace_lcfs

    function chain_distance(state, psi, psi_bound, chain) result(distance)
        !! Largest first-order distance |psi - psi_b| / |grad psi| of the chain
        !! points from the computed boundary.
        type(gs_state_t), intent(in) :: state
        real(dp), intent(in) :: psi(:), psi_bound, chain(:, :)
        real(dp) :: distance
        type(gs_fb_grid_t) :: grid
        real(dp) :: v, g(2)
        integer :: j
        logical :: found

        distance = huge(1.0_dp)
        call build_grid(state, grid)
        distance = 0.0_dp
        do j = 1, size(chain, 2)
            call grid_flux(state, grid, psi, chain(1, j), chain(2, j), v, g, found)
            if (.not. found .or. norm2(g) <= tiny(1.0_dp)) then
                distance = huge(1.0_dp)
                return
            end if
            distance = max(distance, abs(v - psi_bound)/norm2(g))
        end do
    end function chain_distance

    subroutine transfer_flux(old, old_psi, fresh, new_psi)
        !! Finite element interpolation of the flux of mesh `old` at the nodes of
        !! mesh `fresh`. A node outside every old cell (curved outer boundary
        !! discretised differently) takes the value of the nearest cell, with the
        !! reference coordinates clamped to the triangle.
        type(gs_state_t), intent(in) :: old, fresh
        real(dp), intent(in) :: old_psi(:)
        real(dp), allocatable, intent(out) :: new_psi(:)
        type(gs_fb_grid_t) :: grid
        real(dp), allocatable :: centroid(:, :)
        real(dp) :: xi(2), point(2), det, shape(10), grads(2, 10), hess(2, 2, 10)
        real(dp) :: best, d, x(2)
        integer :: i, e, ndof, status, ne, k
        logical :: found

        ndof = size(old%triangles, 1)
        ne = size(old%triangles, 2)
        allocate (new_psi(size(fresh%points, 2)), centroid(2, ne))
        do e = 1, ne
            centroid(:, e) = sum(old%points(:, old%triangles(1:3, e)), dim=2)/3.0_dp
        end do
        call build_grid(old, grid)
        !$omp parallel do default(shared) &
        !$omp private(i, e, xi, point, det, shape, grads, hess, status, best, d) &
        !$omp private(k, x, found) &
        !$omp schedule(static)
        do i = 1, size(fresh%points, 2)
            x = fresh%points(:, i)
            call grid_locate(old, grid, x(1), x(2), e, xi, found)
            if (.not. found) then
                best = huge(1.0_dp)
                e = 1
                do k = 1, ne
                    d = sum((centroid(:, k) - x)**2)
                    if (d < best) then
                        best = d
                        e = k
                    end if
                end do
                xi = matmul(transpose(old%gradient(:, 2:3, e)), &
                            x - old%points(:, old%triangles(1, e)))
                xi = max(xi, 0.0_dp)
                if (sum(xi) > 1.0_dp) xi = xi/sum(xi)
            end if
            call gs_cell_data(old, e, xi, point, det, shape(:ndof), grads(:, :ndof), &
                              hess(:, :, :ndof), status)
            new_psi(i) = dot_product(shape(:ndof), old_psi(old%triangles(:, e)))
        end do
        !$omp end parallel do
    end subroutine transfer_flux

    function gs_fb_chain_movement(state, psi, psi_bound, chain) result(movement)
        !! Largest first-order distance |psi - psi_b| / |grad psi| of the chain
        !! points from the computed boundary, in units of the chain spacing.
        type(gs_state_t), intent(in) :: state
        real(dp), intent(in) :: psi(:), psi_bound, chain(:, :)
        real(dp) :: movement

        movement = chain_distance(state, psi, psi_bound, chain) &
                   /(sum(norm2(chain - cshift(chain, 1, dim=2), dim=1))/size(chain, 2))
    end function gs_fb_chain_movement

    subroutine gs_fb_refit_mesh(state, axis, psi_bound, sigma, contour, max_area, &
                                chain, rcfg, status, curve)
        !! One refit: trace the boundary of the flux in state%psi, rebuild the mesh
        !! (same Gamma_c, element order and curvature) with the Gamma_p chain on
        !! it, and interpolate the flux onto the new mesh. state and chain are
        !! replaced; on failure they are unchanged.
        type(gs_state_t), intent(inout) :: state
        real(dp), intent(in) :: axis(2), psi_bound, sigma, contour(:, :), max_area
        real(dp), intent(inout) :: chain(:, :)
        type(gs_fb_refit_config_t), intent(in) :: rcfg
        integer, intent(out) :: status
        procedure(boundary_curve), optional :: curve
        type(gs_state_t) :: old
        real(dp), allocatable :: psi_new(:), new_chain(:, :)
        integer :: stat
        real(dp) :: t0
        t0 = -1; call gs_phase_tick('', t0)

        call gs_fb_trace_lcfs(state, state%psi, axis, psi_bound, sigma, &
                              size(chain, 2), new_chain, status, &
                              rcfg%dense_factor, rcfg%project_iterations)
        call gs_phase_tick('refit.trace', t0)
        if (status /= GS_FB_OK) return
        old = state
        if (present(curve)) then
            call initialize_gs(state, contour, max_area, stat, old%element_order, &
                               curved=old%mapped, curve=curve, interior_chain=new_chain)
        else
            call initialize_gs(state, contour, max_area, stat, old%element_order, &
                               curved=old%mapped, interior_chain=new_chain)
        end if
        if (stat /= 0) then
            state = old
            status = GS_FB_REFIT_FAILED
            return
        end if
        call gs_phase_tick('refit.initialize_gs', t0)
        chain = new_chain
        call transfer_flux(old, old%psi, state, psi_new)
        call gs_phase_tick('refit.transfer', t0)
        state%psi = psi_new
        status = GS_FB_OK
    end subroutine gs_fb_refit_mesh

    subroutine gs_fb_refit_solve(state, sys, law, cfg, currents, contour, max_area, &
                                 chain, rcfg, report, status, set, curve, extra_cb, &
                                 history)
        !! Newton on the fixed mesh of `state` (psi seed in state%psi), then up to
        !! rcfg%max_refits times: trace the computed boundary, rebuild the mesh with
        !! the Gamma_p chain on it (same Gamma_c, same element order), transfer the
        !! flux by interpolation and re-solve. Stops when the chain lies within
        !! rcfg%tolerance * h of the boundary of the solution computed on its own
        !! mesh (h = chain length / chain points); that mesh is then the final one.
        !! state, sys and chain are returned for the final mesh. extra_cb, if
        !! present, rebuilds the psi-independent load sys%extra_load after every
        !! mesh change (it is mesh dependent).
        type(gs_state_t), intent(inout) :: state
        type(gs_fb_system_t), intent(out) :: sys
        type(gs_profile_law_t), intent(in) :: law
        type(gs_fb_config_t), intent(in) :: cfg
        real(dp), intent(in) :: currents(:), contour(:, :), max_area
        real(dp), intent(inout) :: chain(:, :)
        type(gs_fb_refit_config_t), intent(in) :: rcfg
        type(gs_fb_refit_report_t), intent(out) :: report
        integer, intent(out) :: status
        type(gs_circuit_set_t), intent(in), optional :: set
        procedure(boundary_curve), optional :: curve
        procedure(gs_fb_extra_cb), optional :: extra_cb
        type(gs_fb_history_t), intent(out), optional :: history
        type(gs_fb_config_t) :: current
        type(gs_fb_eval_t) :: ev
        type(gs_fb_history_t) :: last
        integer :: pass, stat, n_chain

        n_chain = size(chain, 2)
        current = cfg
        allocate (report%movement(0:rcfg%max_refits), &
                  report%newton_iterations(0:rcfg%max_refits))
        report%movement = 0.0_dp
        report%newton_iterations = 0
        call solve_stage(state, sys, law, current, currents, set, extra_cb, last, &
                         status)
        report%newton_iterations(0) = last%iterations
        if (status /= GS_FB_OK) return
        do pass = 0, rcfg%max_refits
            call gs_fb_evaluate(state, sys, state%psi, currents, ev, stat, &
                                jacobian=.false.)
            if (stat /= GS_FB_OK) then
                status = stat
                return
            end if
            report%h = sum(norm2(chain - cshift(chain, 1, dim=2), dim=1))/n_chain
            report%movement(pass) = gs_fb_chain_movement(state, state%psi, &
                                                         ev%psi_bound, chain)
            report%psi_bound = ev%psi_bound
            if (report%movement(pass) <= rcfg%tolerance) then
                report%converged = .true.
                exit
            end if
            ! The chain cannot follow the boundary closer than the discretisation
            ! error of the flux: stop once a refit no longer contracts the distance.
            if (pass >= 1) then
                if (report%movement(pass) >= &
                    rcfg%contraction*report%movement(pass - 1)) then
                    report%stagnated = .true.
                    exit
                end if
            end if
            if (pass == rcfg%max_refits) exit
            call gs_fb_refit_mesh(state, ev%axis, ev%psi_bound, &
                                  real(sys%cfg%axis_sign, dp), contour, max_area, &
                                  chain, rcfg, stat, curve)
            if (stat /= GS_FB_OK) then
                status = stat
                return
            end if
            current = sys%cfg
            current%axis_guess = sys%axis_hint
            current%xpoint_guess = sys%xpoint_hint
            call solve_stage(state, sys, law, current, currents, set, extra_cb, last, &
                             status)
            report%newton_iterations(pass + 1) = last%iterations
            if (status /= GS_FB_OK) return
            report%refits = pass + 1
        end do
        if (present(history)) history = last
        status = GS_FB_OK
        if (.not. (report%converged .or. report%stagnated) .and. rcfg%max_refits > 0) &
            status = GS_FB_REFIT_NOT_CONVERGED

    end subroutine gs_fb_refit_solve

    subroutine solve_stage(state, sys, law, cfg, currents, set, extra_cb, history, &
                           status)
        !! Setup, optional rebuild of the nodal extra load, Newton.
        type(gs_state_t), intent(inout) :: state
        type(gs_fb_system_t), intent(out) :: sys
        type(gs_profile_law_t), intent(in) :: law
        type(gs_fb_config_t), intent(in) :: cfg
        real(dp), intent(in) :: currents(:)
        type(gs_circuit_set_t), intent(in), optional :: set
        procedure(gs_fb_extra_cb), optional :: extra_cb
        type(gs_fb_history_t), intent(out) :: history
        integer, intent(out) :: status

        call gs_fb_setup(state, sys, law, cfg, set, status)
        if (status /= GS_FB_OK) return
        if (present(extra_cb)) call extra_cb(state, sys)
        call gs_fb_newton(state, sys, currents, status, history)
    end subroutine solve_stage

end module kin6d_gs_free_boundary
