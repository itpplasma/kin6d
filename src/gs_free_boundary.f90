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
    use kin6d_sort_keys, only: sort_key_order
    use kin6d_grad_shafranov, only: gs_state_t, gs_cell_data
    use kin6d_gs_coils, only: gs_circuit_set_t
    use kin6d_gs_coupling, only: gs_coupling_t, gs_coupling_build, gs_coupled_apply, &
                                 gs_coupling_coil_load, GS_COUPLING_OK
    use kin6d_gs_nonlinear, only: gs_cubic_load, gs_profile_jacobian_t
    use kin6d_gs_profile_law, only: gs_profile_law_t, gs_law_sources
    use kin6d_triangle_basis, only: triangle_basis_jet, cubic_segment_roots, &
                                    cubic_segment_stationary, &
                                    triangle_bernstein_controls
    use kin6d_curved_triangle, only: triangle_locate
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
        logical :: has_jacobian = .false.
    end type gs_fb_eval_t

    type, public :: gs_fb_factor_t
        type(sparse_direct_factor_t) :: lu, lu_t
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

    public :: gs_fb_setup, gs_fb_evaluate, gs_fb_jacobian_apply
    public :: gs_fb_factor, gs_fb_solve, gs_fb_release, gs_fb_newton
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
        call build_pattern(state, sys)
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
        cnt = 0
        do i = 1, n
            do p = sys%coupling%row_ptr(i), sys%coupling%row_ptr(i + 1) - 1
                cnt = cnt + 1
                keys(cnt) = int(sys%coupling%col(p) - 1, i64)*n + (i - 1)
            end do
        end do
        do j = 1, nt
            do i = 1, nt
                cnt = cnt + 1
                keys(cnt) = int(sys%coupling%trace_node(j) - 1, i64)*n &
                            + (sys%coupling%trace_node(i) - 1)
            end do
        end do
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
        do i = 1, n
            do p = sys%coupling%row_ptr(i), sys%coupling%row_ptr(i + 1) - 1
                sys%slot_k(p) = find_slot(unique(:nu), &
                                          int(sys%coupling%col(p) - 1, i64)*n + (i - 1))
            end do
        end do
        do j = 1, nt
            do i = 1, nt
                sys%slot_trace(i, j) = find_slot(unique(:nu), &
                                            int(sys%coupling%trace_node(j) - 1, i64)*n &
                                                 + (sys%coupling%trace_node(i) - 1))
            end do
        end do
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
        real(dp) :: psi_x, level
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
            if (value > best) then
                best = value
                cell_best = e
                axis = xi
            end if
        end do
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
            allocate (shape_x(ndof))
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
                if (dist < best_dist) then
                    best_dist = dist
                    any_found = .true.
                    element = e
                    xpoint = point
                    shape_x = shape(:ndof)
                    psi_x = dot_product(shape(:ndof), values(:ndof))
                end if
            end do
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
        integer :: e, ndof, n, ne, stat, k, i
        integer :: dropped_total
        logical :: need_jacobian

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
        call plasma_cells(state, psi, sigma, ev%psi_bound, k, sys%cfg%mode, &
                          sys%cfg%cut_fraction, ev%axis, ev%xpoint, ev%plasma_cell)
        law = sys%law
        law%axis = ev%psi_axis - ev%psi_bound
        allocate (saved(n), ev%source(n))
        saved = state%psi
        state%psi = psi - ev%psi_bound
        call gs_cubic_load(state, sys%mu0, law, ev%source, jac, stat, &
                           sys%cfg%quadrature_tolerance, &
                           jacobian_tolerance=sys%cfg%jacobian_tolerance, &
                           load_only=.not. need_jacobian, source_integral=integral, &
                           plasma_cell=ev%plasma_cell)
        state%psi = saved
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
                    phi(:ndof) = psi(state%triangles(:, e)) - ev%psi_bound
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
            end if
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
        !! Sparse LU of J0 = A - M and the two Woodbury columns.
        type(gs_state_t), intent(in) :: state
        type(gs_fb_system_t), intent(in) :: sys
        type(gs_fb_eval_t), intent(in) :: ev
        type(gs_fb_factor_t), intent(inout) :: fac
        integer, intent(out) :: status
        real(dp), allocatable :: rhs2(:), col(:)
        integer :: n, nt, ndof, e, k, l, i, j, p, stat

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
        call sparse_direct_factor_csc(fac%lu, n, sys%colptr, sys%rowind, &
                                      fac%values, stat)
        if (stat /= 0) then
            status = GS_FB_FACTOR_FAILED
            return
        end if
        allocate (fac%z(n, 2), fac%zt(n, 2), rhs2(n), col(n))
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
        fac%schur(1, 1) = 1.0_dp - dot_product(ev%ga, fac%z(:, 1))
        fac%schur(1, 2) = -dot_product(ev%ga, fac%z(:, 2))
        fac%schur(2, 1) = -dot_product(ev%gb, fac%z(:, 1))
        fac%schur(2, 2) = 1.0_dp - dot_product(ev%gb, fac%z(:, 2))
        fac%ready = .true.
        status = GS_FB_OK
    end subroutine gs_fb_factor

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
        logical :: adjoint

        status = GS_FB_BAD_INPUT
        if (.not. fac%ready) return
        adjoint = .false.
        if (present(transpose)) adjoint = transpose
        n = sys%n
        allocate (y(n), col(n), rhs_v(n))
        if (.not. adjoint) then
            call sparse_direct_solve_factored(fac%lu, rhs, y, stat)
            if (stat /= 0) then
                status = GS_FB_FACTOR_FAILED
                return
            end if
            w = [dot_product(ev%ga, y), dot_product(ev%gb, y)]
            c = fac%schur
        else
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
            call sparse_direct_solve_factored(fac%lu_t, rhs, y, stat)
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

    subroutine gs_fb_release(fac)
        type(gs_fb_factor_t), intent(inout) :: fac

        if (fac%ready) call sparse_direct_free(fac%lu)
        if (fac%have_transpose) call sparse_direct_free(fac%lu_t)
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

end module kin6d_gs_free_boundary
