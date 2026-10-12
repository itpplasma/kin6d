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
    use, intrinsic :: iso_fortran_env, only: dp => real64, i64 => int64
    use fortfem_axisymmetric_gs_kernel, only: gs_mu0
    use fortfem_axisymmetric_gs_boundary_operators_2d, only: gs_bem_curve_t, &
                         gs_bem_operators_t, gs_bem_curve_from_nodes, gs_bem_assemble, &
                                                             gs_bem_steklov_symmetric
    use fortnum_quadrature, only: gauss_legendre
    use fortnum_krylov, only: real_conjugate_gradient_operator, KRYLOV_OK
    use kin6d_sort_keys, only: sort_key_order
    use kin6d_grad_shafranov, only: gs_state_t, gs_cell_data
    use kin6d_gs_coils, only: gs_circuit_set_t, gs_circuit_flux_vectors
    implicit none
    private

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
    end type gs_coupling_t

    public :: gs_coupling_build, gs_coupled_apply, gs_coupled_matrix
    public :: gs_coupling_coil_load, gs_coupling_source_load, gs_coupled_solve
    public :: gs_coupling_energies, gs_coupling_edge_points

contains

    pure subroutine edge_lagrange(order, t, w, dw)
        !! Lagrange basis (and derivative) on the equispaced nodes j/order.
        integer, intent(in) :: order
        real(dp), intent(in) :: t
        real(dp), intent(out) :: w(0:order), dw(0:order)
        real(dp) :: prod
        integer :: j, l, m

        do j = 0, order
            w(j) = 1.0_dp
            dw(j) = 0.0_dp
            do l = 0, order
                if (l == j) cycle
                w(j) = w(j) * (t - real(l, dp) / order) / (real(j - l, dp) / order)
            end do
            do m = 0, order
                if (m == j) cycle
                prod = 1.0_dp / (real(j - m, dp) / order)
                do l = 0, order
                    if (l == j .or. l == m) cycle
                    prod = prod * (t - real(l, dp) / order) / (real(j - l, dp) / order)
                end do
                dw(j) = dw(j) + prod
            end do
        end do
    end subroutine edge_lagrange

    pure subroutine edge_node_ids(state, e, k, ids)
        !! State nodes of local edge k of triangle e, from vertex k to k+1.
        type(gs_state_t), intent(in) :: state
        integer, intent(in) :: e, k
        integer, intent(out) :: ids(0:)
        integer :: order, j

        order = state % element_order
        ids(0) = state % triangles(k, e)
        do j = 1, order - 1
            ids(j) = state % triangles(3 + (k - 1) * (order - 1) + j, e)
        end do
        ids(order) = state % triangles(mod(k, 3) + 1, e)
    end subroutine edge_node_ids

    subroutine gs_coupling_edge_points(state, coupling, edge, t, x, tangent)
        !! Position and (non-normalised) tangent dx/dt of boundary edge `edge`
        !! at t in [0, 1].
        type(gs_state_t), intent(in) :: state
        type(gs_coupling_t), intent(in) :: coupling
        integer, intent(in) :: edge
        real(dp), intent(in) :: t
        real(dp), intent(out) :: x(2), tangent(2)
        integer :: ids(0:3), order, j
        real(dp) :: w(0:3), dw(0:3)

        order = coupling % order
        call edge_node_ids(state, coupling % edge_element(edge), &
                           coupling % edge_side(edge), ids(0:order))
        call edge_lagrange(order, t, w(0:order), dw(0:order))
        x = 0.0_dp
        tangent = 0.0_dp
        do j = 0, order
            x = x + w(j) * state % points(:, ids(j))
            tangent = tangent + dw(j) * state % points(:, ids(j))
        end do
    end subroutine gs_coupling_edge_points

    subroutine gs_coupling_build(state, coupling, status, panels_per_edge, flux_degree)
        !! Assemble the sparse stiffness and the dense Gamma_c block.
        type(gs_state_t), intent(in) :: state
        type(gs_coupling_t), intent(out) :: coupling
        integer, intent(out) :: status
        integer, intent(in), optional :: panels_per_edge, flux_degree
        type(gs_bem_curve_t) :: curve
        type(gs_bem_operators_t) :: ops
        real(dp), allocatable :: s_bem(:, :), bem_nodes(:, :), restriction(:, :)
        integer, allocatable :: first_edge(:)
        integer :: e, k, i, j, p, q, n, order, ne, np, m, ids(0:3), row, col_j
        integer :: next_edge, stat
        real(dp) :: w(0:3), dw(0:3), t, area

        status = GS_COUPLING_BAD_STATE
        if (.not. state % ready) return
        order = state % element_order
        m = 1
        if (present(panels_per_edge)) m = panels_per_edge
        coupling % panels_per_edge = m
        if (present(flux_degree)) coupling % flux_degree = flux_degree
        if (m < 1) return
        coupling % order = order
        n = size(state % points, 2)
        ne = size(state % triangles, 2)
        coupling % n_nodes = n

        ! Boundary edges in counter-clockwise order, chained by their vertices.
        status = GS_COUPLING_BAD_BOUNDARY
        allocate (first_edge(n))
        first_edge = 0
        coupling % n_edges = 0
        do e = 1, ne
            do k = 1, 3
                if (state % neighbors(k, e) > 0) cycle
                if (first_edge(state % triangles(k, e)) /= 0) return
                first_edge(state % triangles(k, e)) = e * 3 + k - 3
                coupling % n_edges = coupling % n_edges + 1
            end do
        end do
        if (coupling % n_edges * m < 4) return
        allocate (coupling % edge_element(coupling % n_edges), &
                  coupling % edge_side(coupling % n_edges))
        coupling % edge_element = 0
        coupling % edge_side = 0
        do e = 1, ne
            do k = 1, 3
                if (state % neighbors(k, e) <= 0) then
                    coupling % edge_element(1) = e
                    coupling % edge_side(1) = k
                    exit
                end if
            end do
            if (coupling % edge_element(1) /= 0) exit
        end do
        do i = 2, coupling % n_edges
            e = coupling % edge_element(i - 1)
            k = coupling % edge_side(i - 1)
            next_edge = first_edge(state % triangles(mod(k, 3) + 1, e))
            if (next_edge == 0) return
            coupling % edge_element(i) = (next_edge - 1) / 3 + 1
            coupling % edge_side(i) = next_edge - 3 * (coupling % edge_element(i) - 1)
        end do
        e = coupling % edge_element(coupling % n_edges)
        k = coupling % edge_side(coupling % n_edges)
        ! The chain must close on its first edge, visiting every boundary
        ! edge once: Gamma_c is a single loop.
        if (first_edge(state % triangles(mod(k, 3) + 1, e)) /= &
            coupling % edge_element(1) * 3 + coupling % edge_side(1) - 3) return
        first_edge = 0
        do i = 1, coupling % n_edges
            e = coupling % edge_element(i)
            if (first_edge(state % triangles(coupling % edge_side(i), e)) /= 0) return
            first_edge(state % triangles(coupling % edge_side(i), e)) = 1
        end do

        coupling % n_trace = coupling % n_edges * order
        allocate (coupling % trace_node(coupling % n_trace))
        do i = 1, coupling % n_edges
            call edge_node_ids(state, coupling % edge_element(i), &
                               coupling % edge_side(i), ids(0:order))
            coupling % trace_node((i - 1) * order + 1:i * order) = ids(0:order - 1)
        end do

        ! BEM panels: node 2p-1 at the panel start, 2p at its middle.
        np = coupling % n_edges * m
        allocate (bem_nodes(2, 2 * np), restriction(2 * np, coupling % n_trace))
        restriction = 0.0_dp
        do i = 1, coupling % n_edges
            call edge_node_ids(state, coupling % edge_element(i), &
                               coupling % edge_side(i), ids(0:order))
            do q = 1, m
                p = (i - 1) * m + q
                do j = 1, 2
                    t = (real(q - 1, dp) + 0.5_dp * real(j - 1, dp)) / m
                    call edge_lagrange(order, t, w(0:order), dw(0:order))
                    row = 2 * p - 2 + j
                    bem_nodes(:, row) = 0.0_dp
                    do k = 0, order
                        bem_nodes(:, row) = bem_nodes(:, row) &
                                            + w(k) * state % points(:, ids(k))
                        col_j = (i - 1) * order + k + 1
                        if (k == order) col_j = mod(i, coupling % n_edges) * order + 1
                        restriction(row, col_j) = restriction(row, col_j) + w(k)
                    end do
                end do
            end do
        end do
        ! Counter-clockwise orientation (shoelace) is required by the BEM.
        area = 0.0_dp
        do j = 1, 2 * np
            k = mod(j, 2 * np) + 1
            area = area + bem_nodes(1, j) * bem_nodes(2, k) &
                   - bem_nodes(1, k) * bem_nodes(2, j)
        end do
        if (area <= 0.0_dp) return

        status = GS_COUPLING_BEM_FAILED
        call gs_bem_curve_from_nodes(bem_nodes, curve, stat)
        if (stat /= 0) return
        call gs_bem_assemble(curve, coupling % flux_degree, ops, stat)
        if (stat /= 0) return
        call gs_bem_steklov_symmetric(ops, s_bem, stat)
        if (stat /= 0) return
        coupling % s_trace = matmul(transpose(restriction), matmul(s_bem, restriction))
      coupling % s_trace = 0.5_dp * (coupling % s_trace + transpose(coupling % s_trace))

        call build_stiffness(state, coupling)
        status = GS_COUPLING_OK
    end subroutine gs_coupling_build

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
        real(dp), allocatable :: trace(:)
        integer :: i, p

        do i = 1, coupling % n_nodes
            output(i) = 0.0_dp
            do p = coupling % row_ptr(i), coupling % row_ptr(i + 1) - 1
                output(i) = output(i) + coupling % val(p) * input(coupling % col(p))
            end do
        end do
        trace = matmul(coupling % s_trace, input(coupling % trace_node))
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
        !! hand side as sum_k I_k load(:, k).
        type(gs_state_t), intent(in) :: state
        type(gs_coupling_t), intent(in) :: coupling
        type(gs_circuit_set_t), intent(in) :: set
        real(dp), allocatable, intent(out) :: load(:, :)
        integer, intent(out) :: status
        integer, intent(in), optional :: n_gauss
        real(dp), allocatable :: gq(:, :), gr(:, :), gz(:, :), gt(:, :)
        real(dp), allocatable :: xq(:, :), nq(:, :), wq(:), tq(:), c(:, :)
        real(dp) :: x(2), tan(2), speed, w(0:3), dw(0:3), dnu
        real(dp), allocatable :: gx(:), gw(:)
        integer :: ng, order, i, q, j, k, nk, col_j, point

        status = GS_COUPLING_BAD_INPUT
        if (.not. allocated(set % incidence)) return
        nk = size(set % incidence, 2)
        ng = 10
        if (present(n_gauss)) ng = n_gauss
        if (ng < 1) return
        order = coupling % order
        allocate (gx(ng), gw(ng))
        call gauss_legendre(ng, gx, gw)
        gx = 0.5_dp * (gx + 1.0_dp)
        gw = 0.5_dp * gw
        allocate (xq(2, coupling % n_edges * ng), nq(2, coupling % n_edges * ng), &
                  wq(coupling % n_edges * ng), tq(coupling % n_edges * ng))
        point = 0
        do i = 1, coupling % n_edges
            do q = 1, ng
                point = point + 1
                call gs_coupling_edge_points(state, coupling, i, gx(q), x, tan)
                speed = norm2(tan)
                xq(:, point) = x
                nq(:, point) = [tan(2), -tan(1)] / speed
                wq(point) = gw(q) * speed
                tq(point) = gx(q)
            end do
        end do
        allocate (gq(size(xq, 2), nk), gr(size(xq, 2), nk), gz(size(xq, 2), nk), &
                  gt(coupling % n_trace, nk), c(coupling % n_trace, nk))
        call gs_circuit_flux_vectors(set, xq(1, :), xq(2, :), gq, gr, gz)
        call gs_circuit_flux_vectors(set, state % points(1, coupling % trace_node), &
                                     state % points(2, coupling % trace_node), gt)
        c = 0.0_dp
        point = 0
        do i = 1, coupling % n_edges
            do q = 1, ng
                point = point + 1
                call edge_lagrange(order, tq(point), w(0:order), dw(0:order))
                do k = 1, nk
                    dnu = (nq(1, point) * gr(point, k) + nq(2, point) * gz(point, k)) &
                          / xq(1, point)
                    do j = 0, order
                        col_j = (i - 1) * order + j + 1
                        if (j == order) col_j = mod(i, coupling % n_edges) * order + 1
                        c(col_j, k) = c(col_j, k) + wq(point) * dnu * w(j)
                    end do
                end do
            end do
        end do
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
