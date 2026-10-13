module kin6d_gs_coupling
    !! FEM-BEM coupling of the free-boundary Grad-Shafranov problem on a fixed
    !! artificial curve Gamma_c (DESIGN section 1.1).
    !!
    !! Inside Gamma_c the weak form a(psi, v) = int (1/R) grad psi . grad v is
    !! discretised by the curved P1/P2/P3 elements of kin6d_grad_shafranov
    !! (no Dirichlet rows: the nodes flagged `boundary` are ordinary
    !! unknowns). Outside, psi is the coil flux plus an exterior-decaying
    !! harmonic field, and the exterior enters through the symmetric
    !! (Costabel) Steklov-Poincare operator S_h of the fortfem module
    !! fortfem_axisymmetric_gs_boundary_operators_2d:
    !!
    !!     a(psi, v) + <S psi, v>_Gamma_c = int mu0 j v + sum_k I_k <g_k, v>,
    !!     g_k = (S + d_nu) G_k,   d_nu = (1/R) d_n,  n pointing outward.
    !!
    !! The operator K + R^T S_h R is symmetric positive definite (K carries the
    !! constants in its kernel, S_h does not); its quadratic form is the
    !! interior energy int (1/R)|grad psi|^2 plus the energy of the exterior
    !! decaying extension. psi = 0 on the symmetry axis fixes the gauge.
    !!
    !! Trace space. Gamma_c is the union of the boundary edges of the mesh in
    !! counter-clockwise order. The BEM uses continuous piecewise quadratic
    !! traces on panels, `panels_per_edge` panels per mesh edge, with nodes
    !! taken from the curved mesh edge. For P2 elements the restriction R from
    !! the FEM trace to the BEM trace is exact for every panel count; for P3
    !! (and P1) it is interpolation, whose error falls like panels_per_edge**3.
    !! The dense block on the FEM trace dofs is S_trace = R^T S_h R.
    !!
    !! Flux convention. The source of the weak form is f = mu0 j_phi (A/m^2 in,
    !! Wb/rad out), the coil flux is that of kin6d_gs_coils.
    !!
    !! Core boundary layer. Gamma_c, its panels and restriction, the cached
    !! Steklov-Poincare map and the coil flux come from the D6 layer
    !! kin6d_equilibrium_boundary (docs/equilibrium-core.md, M4); this module
    !! is its psi-lane (n = 0, kind GAMMA_C_FLUX) adapter and keeps the
    !! interior stiffness, the coupled operator and the loads.
    use, intrinsic :: iso_fortran_env, only: dp => real64, i64 => int64
    use fortfem_axisymmetric_gs_kernel, only: gs_mu0
    use fortfem_axisymmetric_gs_boundary_operators_2d, only: gs_bem_curve_t
    use fortnum_quadrature, only: gauss_legendre
    use kin6d_equilibrium_boundary, only: gamma_c_t, gamma_c_chain, gamma_c_build, &
        gamma_c_edge_weights, gamma_c_edge_point, gamma_c_quadrature, &
        gamma_c_assemble, gamma_c_steklov, gamma_c_steklov_cache_clear, &
        gamma_c_exterior, axisymmetric_coils_t, GAMMA_C_FLUX, GAMMA_C_OK
    use fortnum_krylov, only: real_conjugate_gradient_operator, KRYLOV_OK
    use kin6d_sort_keys, only: sort_key_order
    use kin6d_gs_timers, only: gs_phase_tick
    use kin6d_grad_shafranov, only: gs_state_t, gs_cell_data
    use kin6d_gs_coils, only: gs_circuit_set_t
    implicit none
    private
    public :: gs_coupling_steklov, gs_coupling_steklov_cache_clear

    integer, parameter, public :: GS_COUPLING_OK = 0
    integer, parameter, public :: GS_COUPLING_BAD_STATE = 1
    integer, parameter, public :: GS_COUPLING_BAD_BOUNDARY = 2
    integer, parameter, public :: GS_COUPLING_BEM_FAILED = 3
    integer, parameter, public :: GS_COUPLING_BAD_INPUT = 4

    !> Pointwise toroidal current density j_phi(R, Z), independent of psi.
    abstract interface
        function gs_density_fn(r, z) result(j)
            import :: dp
            real(dp), intent(in) :: r, z
            real(dp) :: j
        end function gs_density_fn
    end interface
    public :: gs_density_fn

    type, public :: gs_coupling_t
        integer :: order = 0, n_nodes = 0, n_edges = 0, n_trace = 0
        integer :: panels_per_edge = 1, flux_degree = 2
        !> Boundary edge i is local edge edge_side(i) of triangle edge_element(i),
        !> counter-clockwise; its trace dofs are (i-1)*order + 1 .. i*order.
        integer, allocatable :: edge_element(:), edge_side(:)
        integer, allocatable :: trace_node(:)       !< state node of trace dof
        real(dp), allocatable :: s_trace(:, :)      !< R^T S_h R, symmetric
        integer, allocatable :: row_ptr(:), col(:)  !< CSR stiffness
        real(dp), allocatable :: val(:), diag(:)    !< diag: full operator
        type(gamma_c_t) :: gamma_c                  !< the core curve
    end type gs_coupling_t

    public :: gs_coupling_build, gs_coupled_apply, gs_coupled_matrix
    public :: gs_coupling_coil_load, gs_coupling_source_load, gs_coupled_solve
    public :: gs_coupling_energies, gs_coupling_edge_points
    public :: gs_coupling_bem_nodes, gs_coupling_edge_weights


contains

    subroutine gs_coupling_edge_points(state, coupling, edge, t, x, tangent)
        !! Position and (non-normalised) tangent dx/dt of boundary edge `edge`
        !! at t in [0, 1].
        type(gs_state_t), intent(in) :: state
        type(gs_coupling_t), intent(in) :: coupling
        integer, intent(in) :: edge
        real(dp), intent(in) :: t
        real(dp), intent(out) :: x(2), tangent(2)
        if (.not. state % ready) then
            x = 0.0_dp; tangent = 0.0_dp
            return
        end if
        call gamma_c_edge_point(coupling % gamma_c, edge, t, x, tangent)
    end subroutine gs_coupling_edge_points


    subroutine gs_coupling_build(state, coupling, status, panels_per_edge, flux_degree)
        !! Assemble the sparse stiffness and the dense Gamma_c block.
        type(gs_state_t), intent(in) :: state
        type(gs_coupling_t), intent(out) :: coupling
        integer, intent(out) :: status
        integer, intent(in), optional :: panels_per_edge, flux_degree
        integer :: m, stat
        real(dp) :: t0
        t0 = -1; call gs_phase_tick('', t0)

        status = GS_COUPLING_BAD_STATE
        if (.not. state % ready) return
        m = 1
        if (present(panels_per_edge)) m = panels_per_edge
        coupling % panels_per_edge = m
        if (present(flux_degree)) coupling % flux_degree = flux_degree
        if (m < 1) return
        coupling % order = state % element_order
        coupling % n_nodes = size(state % points, 2)

        ! Boundary edges in counter-clockwise order, chained by their vertices;
        ! Gamma_c is a single loop.
        status = GS_COUPLING_BAD_BOUNDARY
        call gamma_c_chain(state % points, state % triangles, state % neighbors, &
                           coupling % order, coupling % edge_element, &
                           coupling % edge_side, coupling % trace_node, stat)
        if (stat /= GAMMA_C_OK) return
        coupling % n_edges = size(coupling % edge_element)
        coupling % n_trace = size(coupling % trace_node)
        call gamma_c_build(state % points(:, coupling % trace_node), coupling % order, &
                           m, coupling % gamma_c, stat)
        if (stat /= GAMMA_C_OK) return
        call gs_phase_tick('cb.chain', t0)
        status = GS_COUPLING_BEM_FAILED
        call gamma_c_exterior(coupling % gamma_c, GAMMA_C_FLUX, 0, &
                              coupling % flux_degree, coupling % s_trace, stat)
        if (stat /= GAMMA_C_OK) return
        t0 = -1; call gs_phase_tick('', t0)
        call build_stiffness(state, coupling)
        call gs_phase_tick('cb.stiffness', t0)
        status = GS_COUPLING_OK
    end subroutine gs_coupling_build

    subroutine gs_coupling_steklov(nodes, flux_degree, curve, s_bem, stat)
        !! Curve and symmetric Steklov-Poincare matrix S_h of the BEM trace
        !! nodes: the n = 0 flux map of the core layer, cached there.
        real(dp), intent(in) :: nodes(:, :)
        integer, intent(in) :: flux_degree
        type(gs_bem_curve_t), intent(out) :: curve
        real(dp), allocatable, intent(out) :: s_bem(:, :)
        integer, intent(out) :: stat
        call gamma_c_steklov(nodes, GAMMA_C_FLUX, 0, flux_degree, curve, s_bem, stat)
    end subroutine gs_coupling_steklov

    subroutine gs_coupling_steklov_cache_clear()
        call gamma_c_steklov_cache_clear()
    end subroutine gs_coupling_steklov_cache_clear

    subroutine gs_coupling_bem_nodes(state, coupling, m, nodes, restriction)
        !! BEM trace nodes of m panels per mesh edge of Gamma_c (node 2p-1 at
        !! the start of panel p, 2p at its middle) and the restriction from
        !! the FEM trace dofs to them (exact for P2 elements).
        type(gs_state_t), intent(in) :: state
        type(gs_coupling_t), intent(in) :: coupling
        integer, intent(in) :: m
        real(dp), allocatable, intent(out) :: nodes(:, :), restriction(:, :)
        type(gamma_c_t) :: panels
        integer :: stat
        call gamma_c_build(state % points(:, coupling % trace_node), coupling % order, &
                           m, panels, stat)
        call move_alloc(panels % nodes, nodes)
        call move_alloc(panels % restriction, restriction)
    end subroutine gs_coupling_bem_nodes

    pure subroutine gs_coupling_edge_weights(coupling, edge, t, cols, w, dw)
        !! Trace dofs cols(0:order) of boundary edge `edge` and the Lagrange
        !! weights (and t-derivatives) of a trace function at t in [0, 1].
        type(gs_coupling_t), intent(in) :: coupling
        integer, intent(in) :: edge
        real(dp), intent(in) :: t
        integer, intent(out) :: cols(0:)
        real(dp), intent(out) :: w(0:)
        real(dp), intent(out), optional :: dw(0:)
        call gamma_c_edge_weights(coupling % gamma_c, edge, t, cols, w, dw)
    end subroutine gs_coupling_edge_weights

    subroutine build_stiffness(state, coupling)
        !! CSR copy of the element stiffness matrices. Duplicates are summed
        !! in element order, so symmetric entries carry identical sums.
        type(gs_state_t), intent(in) :: state
        type(gs_coupling_t), intent(inout) :: coupling
        integer(i64), allocatable :: keys(:)
        real(dp), allocatable :: raw(:)
        integer, allocatable :: order_idx(:), rows(:)
        integer :: e, k, l, ndof, n, cnt, i, nnz, idx
        integer(i64) :: previous, key

        n = coupling % n_nodes
        ndof = size(state % triangles, 1)
        cnt = size(state % triangles, 2) * ndof * ndof
        allocate (keys(cnt), raw(cnt), order_idx(cnt))
        cnt = 0
        do e = 1, size(state % triangles, 2)
            do l = 1, ndof
                do k = 1, ndof
                    cnt = cnt + 1
                    keys(cnt) = int(state % triangles(k, e) - 1, i64) * n &
                                + state % triangles(l, e) - 1
                    raw(cnt) = state % stiffness(k, l, e)
                end do
            end do
        end do
        call sort_key_order(keys, order_idx)
        allocate (coupling % col(cnt), coupling % val(cnt), rows(cnt))
        nnz = 0
        previous = -1
        do i = 1, cnt
            idx = order_idx(i)
            key = keys(idx)
            if (key /= previous) then
                nnz = nnz + 1
                rows(nnz) = int(key / n) + 1
                coupling % col(nnz) = int(mod(key, int(n, i64))) + 1
                coupling % val(nnz) = 0.0_dp
                previous = key
            end if
            coupling % val(nnz) = coupling % val(nnz) + raw(idx)
        end do
        coupling % col = coupling % col(:nnz)
        coupling % val = coupling % val(:nnz)
        allocate (coupling % row_ptr(n + 1), coupling % diag(n))
        coupling % row_ptr = 0
        coupling % diag = 0.0_dp
        do i = 1, nnz
            coupling % row_ptr(rows(i) + 1) = coupling % row_ptr(rows(i) + 1) + 1
          if (coupling % col(i) == rows(i)) coupling % diag(rows(i)) = coupling % val(i)
        end do
        coupling % row_ptr(1) = 1
        do i = 1, n
           coupling % row_ptr(i + 1) = coupling % row_ptr(i + 1) + coupling % row_ptr(i)
        end do
        do i = 1, coupling % n_trace
            idx = coupling % trace_node(i)
            coupling % diag(idx) = coupling % diag(idx) + coupling % s_trace(i, i)
        end do
    end subroutine build_stiffness

    subroutine gs_coupled_apply(coupling, input, output)
        !! output = (K + scatter(S_trace gather)) input.
        type(gs_coupling_t), intent(in) :: coupling
        real(dp), intent(in) :: input(:)
        real(dp), intent(out) :: output(:)
        real(dp), allocatable :: trace(:), gathered(:)
        real(dp) :: acc
        integer :: i, p

        allocate (trace(coupling % n_trace), gathered(coupling % n_trace))
        gathered = input(coupling % trace_node)
        !$omp parallel default(shared) private(i, p, acc)
        !$omp do schedule(static)
        do i = 1, coupling % n_nodes
            acc = 0.0_dp
            do p = coupling % row_ptr(i), coupling % row_ptr(i + 1) - 1
                acc = acc + coupling % val(p) * input(coupling % col(p))
            end do
            output(i) = acc
        end do
        !$omp end do
        !$omp do schedule(static)
        do i = 1, coupling % n_trace
            ! S_trace is symmetric: the column is the contiguous row.
            trace(i) = dot_product(coupling % s_trace(:, i), gathered)
        end do
        !$omp end do
        !$omp end parallel
        ! trace_node entries are distinct, so this scatter has no conflicts.
        output(coupling % trace_node) = output(coupling % trace_node) + trace
    end subroutine gs_coupled_apply

    subroutine gs_coupled_matrix(coupling, a)
        !! Dense coupled matrix; for tests on small meshes only.
        type(gs_coupling_t), intent(in) :: coupling
        real(dp), allocatable, intent(out) :: a(:, :)
        integer :: i, p

        allocate (a(coupling % n_nodes, coupling % n_nodes))
        a = 0.0_dp
        do i = 1, coupling % n_nodes
            do p = coupling % row_ptr(i), coupling % row_ptr(i + 1) - 1
                a(i, coupling % col(p)) = coupling % val(p)
            end do
        end do
        a(coupling % trace_node, coupling % trace_node) = &
            a(coupling % trace_node, coupling % trace_node) + coupling % s_trace
    end subroutine gs_coupled_matrix

    subroutine gs_coupling_energies(coupling, psi, interior, exterior)
        !! Interior energy psi^T K psi and exterior energy t^T S_trace t.
        type(gs_coupling_t), intent(in) :: coupling
        real(dp), intent(in) :: psi(:)
        real(dp), intent(out) :: interior, exterior
        real(dp), allocatable :: action(:), trace(:)
        integer :: i, p

        allocate (action(coupling % n_nodes))
        do i = 1, coupling % n_nodes
            action(i) = 0.0_dp
            do p = coupling % row_ptr(i), coupling % row_ptr(i + 1) - 1
                action(i) = action(i) + coupling % val(p) * psi(coupling % col(p))
            end do
        end do
        interior = dot_product(psi, action)
        trace = psi(coupling % trace_node)
        exterior = dot_product(trace, matmul(coupling % s_trace, trace))
    end subroutine gs_coupling_energies

    subroutine gs_coupling_coil_load(state, coupling, set, load, status, n_gauss)
        !! load(:, k) = scatter of (S_trace gamma_k + c_k) for each circuit k:
        !! gamma_k the trace values of the circuit flux g_k and
        !! c_k(i) = int_Gamma_c d_nu g_k phi_i ds by Gauss quadrature on the
        !! curved edges. The coil flux psi = sum_k I_k g_k then enters the right
        !! hand side as sum_k I_k load(:, k). The flux comes through the core
        !! coil interface (axisymmetric circuits).
        type(gs_state_t), intent(in) :: state
        type(gs_coupling_t), intent(in) :: coupling
        type(gs_circuit_set_t), intent(in) :: set
        real(dp), allocatable, intent(out) :: load(:, :)
        integer, intent(out) :: status
        integer, intent(in), optional :: n_gauss
        type(axisymmetric_coils_t) :: coils
        real(dp), allocatable :: gq(:, :), gr(:, :), gz(:, :), gt(:, :), dnu(:, :)
        real(dp), allocatable :: xq(:, :), nq(:, :), wq(:), tq(:), c(:, :)
        integer :: ng, k, nk, point, stat
        status = GS_COUPLING_BAD_INPUT
        if (.not. allocated(set % incidence)) return
        nk = size(set % incidence, 2)
        ng = 10
        if (present(n_gauss)) ng = n_gauss
        if (ng < 1) return
        coils % set = set
        call gamma_c_quadrature(coupling % gamma_c, ng, xq, nq, wq, tq)
        allocate (gq(size(xq, 2), nk), gr(size(xq, 2), nk), gz(size(xq, 2), nk), &
                  gt(coupling % n_trace, nk), c(coupling % n_trace, nk), &
                  dnu(size(xq, 2), nk))
        call coils % flux(xq, gq, gr, gz, stat)
        if (stat /= GAMMA_C_OK) return
        call coils % flux(state % points(:, coupling % trace_node), gt, status=stat)
        if (stat /= GAMMA_C_OK) return
        do point = 1, size(xq, 2)
            do k = 1, nk
                dnu(point, k) = (nq(1, point) * gr(point, k) &
                                 + nq(2, point) * gz(point, k)) / xq(1, point)
            end do
        end do
        call gamma_c_assemble(coupling % gamma_c, tq, wq, dnu, c)
        c = c + matmul(coupling % s_trace, gt)
        allocate (load(coupling % n_nodes, nk))
        load = 0.0_dp
        do k = 1, nk
            load(coupling % trace_node, k) = load(coupling % trace_node, k) + c(:, k)
        end do
        status = GS_COUPLING_OK
    end subroutine gs_coupling_coil_load

    subroutine gs_coupling_source_load(state, density, load, status, n_gauss)
        !! load(i) = int mu0 j_phi(R, Z) phi_i dA, collapsed tensor Gauss rule of
        !! n_gauss x n_gauss points per curved cell. j_phi may be supported on
        !! the interior of a fitted Gamma_p chain and need only be smooth on
        !! each cell.
        type(gs_state_t), intent(in) :: state
        procedure(gs_density_fn) :: density
        real(dp), allocatable, intent(out) :: load(:)
        integer, intent(out) :: status
        integer, intent(in), optional :: n_gauss
        real(dp), allocatable :: gx(:), gw(:)
        real(dp) :: shape(10), grads(2, 10), hess(2, 2, 10), point(2), det, xi(2)
        real(dp) :: weight
        integer :: ng, e, a, b, l, ndof, stat

        status = GS_COUPLING_BAD_INPUT
        if (.not. state % ready) return
        ng = 8
        if (present(n_gauss)) ng = n_gauss
        if (ng < 1) return
        allocate (gx(ng), gw(ng))
        call gauss_legendre(ng, gx, gw)
        gx = 0.5_dp * (gx + 1.0_dp)
        gw = 0.5_dp * gw
        ndof = size(state % triangles, 1)
        allocate (load(size(state % points, 2)))
        load = 0.0_dp
        do e = 1, size(state % triangles, 2)
            do a = 1, ng
                do b = 1, ng
                    xi = [gx(a), (1.0_dp - gx(a)) * gx(b)]
                    call gs_cell_data(state, e, xi, point, det, shape(:ndof), &
                                      grads(:, :ndof), hess(:, :, :ndof), stat)
                    if (stat /= 0) return
                    weight = gw(a) * gw(b) * (1.0_dp - gx(a)) * det * gs_mu0 * &
                             density(point(1), point(2))
                    do l = 1, ndof
                        load(state % triangles(l, e)) = load(state % triangles(l, e)) &
                                                        + weight * shape(l)
                    end do
                end do
            end do
        end do
        status = GS_COUPLING_OK
    end subroutine gs_coupling_source_load

    subroutine gs_coupled_solve(coupling, rhs, psi, status, relative_tolerance, &
                                iterations)
        !! Jacobi-preconditioned CG on the symmetric positive definite operator.
        type(gs_coupling_t), intent(in) :: coupling
        real(dp), intent(in) :: rhs(:)
        real(dp), allocatable, intent(out) :: psi(:)
        integer, intent(out) :: status
        real(dp), intent(in), optional :: relative_tolerance
        integer, intent(out), optional :: iterations
        real(dp) :: tolerance, residual
        integer :: iters

        status = GS_COUPLING_BAD_INPUT
        if (size(rhs) /= coupling % n_nodes) return
        allocate (psi(coupling % n_nodes))
        psi = 0.0_dp
        tolerance = 1.0e-12_dp
        if (present(relative_tolerance)) tolerance = relative_tolerance
        tolerance = tolerance * min(1.0_dp, norm2(rhs))
        if (tolerance <= 0.0_dp) tolerance = tiny(1.0_dp)
        call real_conjugate_gradient_operator(matvec, rhs, psi, tolerance, &
                  max(2000, 20 * ceiling(sqrt(real(coupling % n_nodes, dp)))), status, &
                                              iters, residual, preconditioner)
        if (present(iterations)) iterations = iters
        if (status == KRYLOV_OK) status = GS_COUPLING_OK
    contains
        subroutine matvec(input, output)
            real(dp), intent(in) :: input(:)
            real(dp), intent(out) :: output(:)
            call gs_coupled_apply(coupling, input, output)
        end subroutine matvec
        subroutine preconditioner(input, output)
            real(dp), intent(in) :: input(:)
            real(dp), intent(out) :: output(:)
            output = input / coupling % diag
        end subroutine preconditioner
    end subroutine gs_coupled_solve

end module kin6d_gs_coupling
