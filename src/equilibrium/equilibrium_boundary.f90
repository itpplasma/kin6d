module kin6d_equilibrium_boundary
    ! One boundary layer for every lane (docs/equilibrium-core.md, D6,
    ! migration M4). A fixed axisymmetric curve Gamma_c in vacuum encloses the
    ! plasma, so the exterior separates in the toroidal mode n. The layer owns
    !   the curve     Gamma_c as a closed counter-clockwise chain of curved
    !                 Lagrange edges of order 1..3 (the mesh boundary of a
    !                 lane), its P2 BEM panels (panels_per_edge per edge) and
    !                 the restriction R from the FEM trace to the panel trace,
    !                 exact for order 2;
    !   the exterior  per mode n the symmetric (Costabel) Steklov-Poincare map
    !                 S_n of the decaying exterior field, on the panel trace
    !                 (cached per curve, kind, n and flux degree) and on the
    !                 FEM trace, S_trace = R^T S_n R;
    !   the coils     one coil interface (axisymmetric circuits of
    !                 kin6d_gs_coils, 3D polylines of kin6d_coils_3d) that
    !                 returns per unit circuit current the toroidal harmonics
    !                 of B and, where the set has one, the poloidal flux;
    !   the traces    the per-n traces of a coil field on Gamma_c.
    ! Two exterior kinds:
    !   GAMMA_C_SCALAR  the scalar potential Phi_n of B = grad Phi outside
    !                   Gamma_c, L_n Phi_n = -div(R grad Phi_n) + n^2 Phi_n/R
    !                   = 0 (fortfem_toroidal_mode_boundary_operators_2d); the
    !                   trace is Phi_n, the conormal flux R d_n Phi_n = R B.n,
    !                   and (S_n u, v) = -int R B.n v ds;
    !   GAMMA_C_FLUX    the poloidal flux psi = R A_phi of the n = 0 field,
    !                   weight 1/R (fortfem_axisymmetric_gs_boundary_
    !                   operators_2d); this is the n = 1 vector mode of the
    !                   toroidal kernel and the psi lane's coupling.
    ! Units are SI. Fourier harmonics are two-sided, X = sum_n X_n e^{i n phi}
    ! with X_{-n} = conjg(X_n), as in kin6d_coils_3d.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortfem_axisymmetric_gs_boundary_operators_2d, only: gs_bem_curve_t, &
        gs_bem_operators_t, gs_bem_curve_from_nodes, gs_bem_assemble, &
        gs_bem_steklov_symmetric
    use fortfem_toroidal_mode_boundary_operators_2d, only: tm_bem_assemble, &
        tm_bem_steklov_symmetric
    use fortnum_quadrature, only: gauss_legendre
    use kin6d_gs_timers, only: gs_phase_tick
    use kin6d_gs_coils, only: gs_circuit_set_t, gs_circuit_flux_vectors
    use kin6d_coils_3d, only: coil_set_t, coil_set_modes
    use kin6d_equilibrium_state, only: EQ_UNAVAILABLE
    use kin6d_real_identity, only: same_real_value
    implicit none
    private
    public :: gamma_c_t, gamma_c_chain, gamma_c_build, gamma_c_edge_weights
    public :: gamma_c_edge_point, gamma_c_quadrature, gamma_c_assemble
    public :: gamma_c_steklov, gamma_c_steklov_cache_clear, gamma_c_exterior
    public :: boundary_coils_t, axisymmetric_coils_t, polyline_coils_t
    public :: gamma_c_coil_traces

    integer, parameter, public :: GAMMA_C_SCALAR = 1, GAMMA_C_FLUX = 2
    integer, parameter, public :: GAMMA_C_OK = 0, GAMMA_C_BAD_INPUT = 1, &
        GAMMA_C_BAD_BOUNDARY = 2, GAMMA_C_BEM_FAILED = 3

    type :: gamma_c_t
        integer :: order = 0, n_edges = 0, n_trace = 0, panels_per_edge = 1
        ! Trace dof positions: dofs (i-1)*order + 1 .. i*order are the nodes
        ! 0..order-1 of edge i; node `order` is the first dof of edge i+1.
        real(dp), allocatable :: trace_points(:, :)
        ! BEM nodes (2, 2 panels): node 2p-1 at the start of panel p, 2p at
        ! its middle; restriction (2 panels, n_trace) from the FEM trace.
        real(dp), allocatable :: nodes(:, :), restriction(:, :)
    end type

    ! One cached Steklov-Poincare operator of a fixed curve, keyed by the BEM
    ! nodes, the kind, the mode and the flux degree. Gamma_c does not change
    ! between refits, Newton iterations, inverse passes or modes of a sweep,
    ! so the dense BEM assembly is done once per (kind, n).
    type :: steklov_entry_t
        integer :: kind = 0, n = 0, flux_degree = -1
        real(dp), allocatable :: nodes(:, :), s(:, :)
        type(gs_bem_curve_t) :: curve
    end type
    integer, parameter :: steklov_slots = 8
    type(steklov_entry_t), save :: steklov_cache(steklov_slots)
    integer, save :: steklov_next = 1

    ! The coil interface: per unit current of each circuit, the toroidal
    ! harmonics b(c, m, i, k) of the cylindrical field components c = (R,
    ! phi, Z) for modes(m) at points(:, i) = (R, Z), and the n = 0 poloidal
    ! flux g(i, k) (Wb/rad) with its (R, Z) gradient when the set has one.
    type, abstract :: boundary_coils_t
    contains
        procedure(coils_count), deferred :: circuits
        procedure(coils_modes), deferred :: field_modes
        procedure(coils_flux), deferred :: flux
    end type
    abstract interface
        integer function coils_count(coils)
            import :: boundary_coils_t
            class(boundary_coils_t), intent(in) :: coils
        end function
        subroutine coils_modes(coils, points, modes, b, status)
            import :: boundary_coils_t, dp
            class(boundary_coils_t), intent(in) :: coils
            real(dp), intent(in) :: points(:, :)
            integer, intent(in) :: modes(:)
            complex(dp), intent(out) :: b(:, :, :, :)
            integer, intent(out) :: status
        end subroutine
        subroutine coils_flux(coils, points, g, dg_dr, dg_dz, status)
            import :: boundary_coils_t, dp
            class(boundary_coils_t), intent(in) :: coils
            real(dp), intent(in) :: points(:, :)
            real(dp), intent(out) :: g(:, :)
            real(dp), intent(out), optional :: dg_dr(:, :), dg_dz(:, :)
            integer, intent(out) :: status
        end subroutine
    end interface

    ! Axisymmetric circuits (PF/CS coils): n = 0 only, from the exact flux.
    type, extends(boundary_coils_t) :: axisymmetric_coils_t
        type(gs_circuit_set_t) :: set
    contains
        procedure :: circuits => axisymmetric_circuits
        procedure :: field_modes => axisymmetric_modes
        procedure :: flux => axisymmetric_flux
    end type

    ! 3D polylines; circuit k carries incidence(:, k) times its current in
    ! the coils (identity when unallocated). The harmonics come from an
    ! nphi-point DFT of the Biot-Savart field; there is no flux.
    type, extends(boundary_coils_t) :: polyline_coils_t
        type(coil_set_t) :: set
        real(dp), allocatable :: incidence(:, :)
        integer :: nphi = 128
    contains
        procedure :: circuits => polyline_circuits
        procedure :: field_modes => polyline_modes
        procedure :: flux => polyline_flux
    end type
contains
    ! ------------------------------------------------------------------
    ! The curve
    ! ------------------------------------------------------------------

    subroutine gamma_c_chain(points, triangles, neighbors, order, edge_element, &
                             edge_side, trace_node, status)
        ! Counter-clockwise chain of the boundary edges of a mesh (local edge
        ! k of triangle e runs from vertex k to k+1; neighbors(k, e) <= 0 on
        ! the boundary), started at the first boundary edge in element order.
        ! Gamma_c must be a single loop. trace_node: mesh node of each trace
        ! dof. Higher-order edge nodes of edge k follow the vertices in
        ! triangles(4:, e), order - 1 per edge.
        real(dp), intent(in) :: points(:, :)
        integer, intent(in) :: triangles(:, :), neighbors(:, :), order
        integer, allocatable, intent(out) :: edge_element(:), edge_side(:), &
            trace_node(:)
        integer, intent(out) :: status
        integer, allocatable :: first_edge(:)
        integer :: e, k, i, j, n, ne, n_edges, next_edge, ids(0:3)
        status = GAMMA_C_BAD_BOUNDARY
        n = size(points, 2)
        ne = size(triangles, 2)
        allocate (first_edge(n))
        first_edge = 0
        n_edges = 0
        do e = 1, ne
            do k = 1, 3
                if (neighbors(k, e) > 0) cycle
                if (first_edge(triangles(k, e)) /= 0) return
                first_edge(triangles(k, e)) = e*3 + k - 3
                n_edges = n_edges + 1
            end do
        end do
        if (n_edges < 1) return
        allocate (edge_element(n_edges), edge_side(n_edges))
        edge_element = 0
        edge_side = 0
        do e = 1, ne
            do k = 1, 3
                if (neighbors(k, e) <= 0) then
                    edge_element(1) = e
                    edge_side(1) = k
                    exit
                end if
            end do
            if (edge_element(1) /= 0) exit
        end do
        do i = 2, n_edges
            e = edge_element(i - 1)
            k = edge_side(i - 1)
            next_edge = first_edge(triangles(mod(k, 3) + 1, e))
            if (next_edge == 0) return
            edge_element(i) = (next_edge - 1)/3 + 1
            edge_side(i) = next_edge - 3*(edge_element(i) - 1)
        end do
        e = edge_element(n_edges)
        k = edge_side(n_edges)
        if (first_edge(triangles(mod(k, 3) + 1, e)) /= &
            edge_element(1)*3 + edge_side(1) - 3) return
        first_edge = 0
        do i = 1, n_edges
            e = edge_element(i)
            if (first_edge(triangles(edge_side(i), e)) /= 0) return
            first_edge(triangles(edge_side(i), e)) = 1
        end do
        allocate (trace_node(n_edges*order))
        do i = 1, n_edges
            e = edge_element(i)
            k = edge_side(i)
            ids(0) = triangles(k, e)
            do j = 1, order - 1
                ids(j) = triangles(3 + (k - 1)*(order - 1) + j, e)
            end do
            trace_node((i - 1)*order + 1:i*order) = ids(0:order - 1)
        end do
        status = GAMMA_C_OK
    end subroutine

    subroutine gamma_c_build(trace_points, order, panels_per_edge, gc, status)
        ! The curve from its trace dof positions (see gamma_c_t), its BEM
        ! nodes and the FEM -> BEM trace restriction. The panel chain must be
        ! counter-clockwise and have at least 4 panels.
        real(dp), intent(in) :: trace_points(:, :)
        integer, intent(in) :: order, panels_per_edge
        type(gamma_c_t), intent(out) :: gc
        integer, intent(out) :: status
        real(dp) :: t, w(0:3), area
        integer :: i, q, p, j, k, row, np, cols(0:3)
        status = GAMMA_C_BAD_INPUT
        if (order < 1 .or. order > 3 .or. panels_per_edge < 1) return
        if (size(trace_points, 1) /= 2) return
        if (mod(size(trace_points, 2), order) /= 0) return
        gc%order = order
        gc%panels_per_edge = panels_per_edge
        gc%n_trace = size(trace_points, 2)
        gc%n_edges = gc%n_trace/order
        gc%trace_points = trace_points
        status = GAMMA_C_BAD_BOUNDARY
        if (gc%n_edges*panels_per_edge < 4) return
        np = gc%n_edges*panels_per_edge
        allocate (gc%nodes(2, 2*np), gc%restriction(2*np, gc%n_trace))
        gc%restriction = 0.0_dp
        do i = 1, gc%n_edges
            do q = 1, panels_per_edge
                p = (i - 1)*panels_per_edge + q
                do j = 1, 2
                    t = (real(q - 1, dp) + 0.5_dp*real(j - 1, dp))/panels_per_edge
                    call gamma_c_edge_weights(gc, i, t, cols, w)
                    row = 2*p - 2 + j
                    gc%nodes(:, row) = 0.0_dp
                    do k = 0, order
                        gc%nodes(:, row) = gc%nodes(:, row) + w(k) &
                            *trace_points(:, cols(k))
                        gc%restriction(row, cols(k)) = gc%restriction(row, cols(k)) &
                            + w(k)
                    end do
                end do
            end do
        end do
        ! Counter-clockwise orientation (shoelace) is required by the BEM.
        area = 0.0_dp
        do j = 1, 2*np
            k = mod(j, 2*np) + 1
            area = area + gc%nodes(1, j)*gc%nodes(2, k) - gc%nodes(1, k)*gc%nodes(2, j)
        end do
        if (area <= 0.0_dp) return
        status = GAMMA_C_OK
    end subroutine

    pure subroutine edge_lagrange(order, t, w, dw)
        ! Lagrange basis (and derivative) on the equispaced nodes j/order.
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
                w(j) = w(j)*(t - real(l, dp)/order)/(real(j - l, dp)/order)
            end do
            do m = 0, order
                if (m == j) cycle
                prod = 1.0_dp/(real(j - m, dp)/order)
                do l = 0, order
                    if (l == j .or. l == m) cycle
                    prod = prod*(t - real(l, dp)/order)/(real(j - l, dp)/order)
                end do
                dw(j) = dw(j) + prod
            end do
        end do
    end subroutine

    pure subroutine gamma_c_edge_weights(gc, edge, t, cols, w, dw)
        ! Trace dofs cols(0:order) of edge `edge` and the Lagrange weights
        ! (and t-derivatives) of a trace function at t in [0, 1].
        type(gamma_c_t), intent(in) :: gc
        integer, intent(in) :: edge
        real(dp), intent(in) :: t
        integer, intent(out) :: cols(0:)
        real(dp), intent(out) :: w(0:)
        real(dp), intent(out), optional :: dw(0:)
        real(dp) :: d(0:3)
        integer :: j, order
        order = gc%order
        call edge_lagrange(order, t, w(0:order), d(0:order))
        if (present(dw)) dw(0:order) = d(0:order)
        do j = 0, order
            cols(j) = (edge - 1)*order + j + 1
        end do
        cols(order) = mod(edge, gc%n_edges)*order + 1
    end subroutine

    pure subroutine gamma_c_edge_point(gc, edge, t, x, tangent)
        ! Position and (non-normalised) tangent dx/dt of edge `edge` at t.
        type(gamma_c_t), intent(in) :: gc
        integer, intent(in) :: edge
        real(dp), intent(in) :: t
        real(dp), intent(out) :: x(2), tangent(2)
        real(dp) :: w(0:3), dw(0:3)
        integer :: cols(0:3), j
        call gamma_c_edge_weights(gc, edge, t, cols, w, dw)
        x = 0.0_dp
        tangent = 0.0_dp
        do j = 0, gc%order
            x = x + w(j)*gc%trace_points(:, cols(j))
            tangent = tangent + dw(j)*gc%trace_points(:, cols(j))
        end do
    end subroutine

    subroutine gamma_c_quadrature(gc, ng, x, normal, weight, t)
        ! ng-point Gauss rule on every edge, edge-major: points x, outward
        ! unit normals, arc-length weights and edge parameters t.
        type(gamma_c_t), intent(in) :: gc
        integer, intent(in) :: ng
        real(dp), allocatable, intent(out) :: x(:, :), normal(:, :), weight(:), t(:)
        real(dp) :: gx(ng), gw(ng), tan(2), speed
        integer :: i, q, point
        call gauss_legendre(ng, gx, gw)
        gx = 0.5_dp*(gx + 1.0_dp)
        gw = 0.5_dp*gw
        allocate (x(2, gc%n_edges*ng), normal(2, gc%n_edges*ng), &
                  weight(gc%n_edges*ng), t(gc%n_edges*ng))
        point = 0
        do i = 1, gc%n_edges
            do q = 1, ng
                point = point + 1
                call gamma_c_edge_point(gc, i, gx(q), x(:, point), tan)
                speed = norm2(tan)
                normal(:, point) = [tan(2), -tan(1)]/speed
                weight(point) = gw(q)*speed
                t(point) = gx(q)
            end do
        end do
    end subroutine

    subroutine gamma_c_assemble(gc, t, weight, values, c)
        ! c(j, k) = sum_points weight values(point, k) N_j for the trace basis
        ! N_j, at the edge-major points of gamma_c_quadrature (t, weight).
        type(gamma_c_t), intent(in) :: gc
        real(dp), intent(in) :: t(:), weight(:), values(:, :)
        real(dp), intent(out) :: c(:, :)
        real(dp) :: w(0:3)
        integer :: ng, i, q, k, j, point, cols(0:3)
        ng = size(t)/gc%n_edges
        c = 0.0_dp
        point = 0
        do i = 1, gc%n_edges
            do q = 1, ng
                point = point + 1
                call gamma_c_edge_weights(gc, i, t(point), cols, w)
                do k = 1, size(values, 2)
                    do j = 0, gc%order
                        c(cols(j), k) = c(cols(j), k) + weight(point)*values(point, k) &
                            *w(j)
                    end do
                end do
            end do
        end do
    end subroutine

    ! ------------------------------------------------------------------
    ! The exterior
    ! ------------------------------------------------------------------

    subroutine gamma_c_steklov(nodes, kind, n, flux_degree, curve, s_bem, stat)
        ! Curve and symmetric Steklov-Poincare matrix S_n of the BEM trace
        ! nodes, from the cache when the same nodes were assembled before for
        ! the same (kind, n, flux_degree), bitwise equal coordinates, also
        ! when the closed curve is entered at another panel. Not for use
        ! inside a parallel region; the assembly itself is threaded.
        real(dp), intent(in) :: nodes(:, :)
        integer, intent(in) :: kind, n, flux_degree
        type(gs_bem_curve_t), intent(out) :: curve
        real(dp), allocatable, intent(out) :: s_bem(:, :)
        integer, intent(out) :: stat
        type(gs_bem_operators_t) :: ops
        integer :: slot, shift, nn, i, j
        real(dp) :: t0
        t0 = -1
        stat = GAMMA_C_BAD_INPUT
        if (kind /= GAMMA_C_SCALAR .and. kind /= GAMMA_C_FLUX) return
        if (kind == GAMMA_C_FLUX .and. n /= 0) return
        do slot = 1, steklov_slots
            if (steklov_cache(slot)%flux_degree /= flux_degree) cycle
            if (steklov_cache(slot)%kind /= kind) cycle
            if (steklov_cache(slot)%n /= n) cycle
            if (.not. allocated(steklov_cache(slot)%nodes)) cycle
            if (any(shape(steklov_cache(slot)%nodes) /= shape(nodes))) cycle
            shift = steklov_shift(steklov_cache(slot)%nodes, nodes)
            if (shift < 0) cycle
            nn = size(nodes, 2)
            if (shift == 0) then
                curve = steklov_cache(slot)%curve
                s_bem = steklov_cache(slot)%s
            else
                ! The same closed curve entered at another panel: relabel.
                call gs_bem_curve_from_nodes(nodes, curve, stat)
                if (stat /= 0) return
                allocate (s_bem(nn, nn))
                !$omp parallel do default(shared) private(j, i) schedule(static)
                do j = 1, nn
                    do i = 1, nn
                        s_bem(i, j) = steklov_cache(slot)%s( &
                            mod(i - 1 + shift, nn) + 1, mod(j - 1 + shift, nn) + 1)
                    end do
                end do
                !$omp end parallel do
            end if
            stat = 0
            return
        end do
        call gs_bem_curve_from_nodes(nodes, curve, stat)
        if (stat /= 0) return
        call gs_phase_tick('', t0)
        if (kind == GAMMA_C_FLUX) then
            call gs_bem_assemble(curve, flux_degree, ops, stat)
        else
            call tm_bem_assemble(curve, n, flux_degree, ops, stat)
        end if
        if (stat /= 0) return
        call gs_phase_tick('bem.assemble', t0)
        if (kind == GAMMA_C_FLUX) then
            call gs_bem_steklov_symmetric(ops, s_bem, stat)
        else
            call tm_bem_steklov_symmetric(ops, s_bem, stat)
        end if
        call gs_phase_tick('bem.steklov_dense', t0)
        if (stat /= 0) return
        slot = steklov_next
        steklov_next = mod(steklov_next, steklov_slots) + 1
        steklov_cache(slot)%kind = kind
        steklov_cache(slot)%n = n
        steklov_cache(slot)%flux_degree = flux_degree
        steklov_cache(slot)%nodes = nodes
        steklov_cache(slot)%curve = curve
        steklov_cache(slot)%s = s_bem
    end subroutine

    pure integer function steklov_shift(cached, nodes) result(shift)
        ! Even s >= 0 with nodes(:, j) == cached(:, j + s) cyclically for every
        ! j (panels hold two nodes), or -1 when the point lists differ.
        real(dp), intent(in) :: cached(:, :), nodes(:, :)
        integer :: nn, j, s0
        nn = size(nodes, 2)
        shift = -1
        do s0 = 0, nn - 2, 2
            if (all(same_real_value(cached(:, s0 + 1), nodes(:, 1)))) then
                do j = 1, nn
                    if (any(.not. same_real_value( &
                        cached(:, mod(j - 1 + s0, nn) + 1), nodes(:, j)))) exit
                    if (j == nn) shift = s0
                end do
                if (shift >= 0) return
            end if
        end do
    end function

    subroutine gamma_c_steklov_cache_clear()
        integer :: slot
        do slot = 1, steklov_slots
            steklov_cache(slot)%flux_degree = -1
            steklov_cache(slot)%kind = 0
            steklov_cache(slot)%n = 0
            if (allocated(steklov_cache(slot)%nodes)) &
                deallocate (steklov_cache(slot)%nodes)
            if (allocated(steklov_cache(slot)%s)) deallocate (steklov_cache(slot)%s)
        end do
        steklov_next = 1
    end subroutine

    subroutine gamma_c_exterior(gc, kind, n, flux_degree, s_trace, status, curve, &
                                s_bem)
        ! S_trace = R^T S_n R on the FEM trace of Gamma_c, symmetrised; the
        ! panel curve and the panel-level S_n on request.
        type(gamma_c_t), intent(in) :: gc
        integer, intent(in) :: kind, n, flux_degree
        real(dp), allocatable, intent(out) :: s_trace(:, :)
        integer, intent(out) :: status
        type(gs_bem_curve_t), intent(out), optional :: curve
        real(dp), allocatable, intent(out), optional :: s_bem(:, :)
        type(gs_bem_curve_t) :: panels
        real(dp), allocatable :: s(:, :), work(:, :)
        integer :: j, stat
        real(dp) :: t0
        t0 = -1; call gs_phase_tick('', t0)
        status = GAMMA_C_BAD_INPUT
        if (.not. allocated(gc%nodes)) return
        status = GAMMA_C_BEM_FAILED
        call gamma_c_steklov(gc%nodes, kind, n, flux_degree, panels, s, stat)
        if (stat /= 0) return
        call gs_phase_tick('cb.steklov', t0)
        ! One column per thread task.
        allocate (s_trace(gc%n_trace, gc%n_trace))
        allocate (work(size(gc%nodes, 2), gc%n_trace))
        !$omp parallel do default(shared) private(j) schedule(static)
        do j = 1, gc%n_trace
            work(:, j) = matmul(s, gc%restriction(:, j))
        end do
        !$omp end parallel do
        !$omp parallel do default(shared) private(j) schedule(static)
        do j = 1, gc%n_trace
            s_trace(:, j) = matmul(transpose(gc%restriction), work(:, j))
        end do
        !$omp end parallel do
        s_trace = 0.5_dp*(s_trace + transpose(s_trace))
        call gs_phase_tick('cb.strace', t0)
        if (present(curve)) curve = panels
        if (present(s_bem)) call move_alloc(s, s_bem)
        status = GAMMA_C_OK
    end subroutine

    ! ------------------------------------------------------------------
    ! The coils
    ! ------------------------------------------------------------------

    integer function axisymmetric_circuits(coils) result(k)
        class(axisymmetric_coils_t), intent(in) :: coils
        k = 0
        if (allocated(coils%set%incidence)) k = size(coils%set%incidence, 2)
    end function

    subroutine axisymmetric_flux(coils, points, g, dg_dr, dg_dz, status)
        class(axisymmetric_coils_t), intent(in) :: coils
        real(dp), intent(in) :: points(:, :)
        real(dp), intent(out) :: g(:, :)
        real(dp), intent(out), optional :: dg_dr(:, :), dg_dz(:, :)
        integer, intent(out) :: status
        status = GAMMA_C_BAD_INPUT
        if (coils%circuits() == 0) return
        call gs_circuit_flux_vectors(coils%set, points(1, :), points(2, :), g, &
                                     dg_dr, dg_dz)
        status = GAMMA_C_OK
    end subroutine

    subroutine axisymmetric_modes(coils, points, modes, b, status)
        ! n = 0: B_R = -(1/R) d psi/dZ, B_Z = (1/R) d psi/dR, B_phi = 0.
        class(axisymmetric_coils_t), intent(in) :: coils
        real(dp), intent(in) :: points(:, :)
        integer, intent(in) :: modes(:)
        complex(dp), intent(out) :: b(:, :, :, :)
        integer, intent(out) :: status
        real(dp), allocatable :: g(:, :), gr(:, :), gz(:, :)
        integer :: m, i, nk
        nk = coils%circuits()
        status = GAMMA_C_BAD_INPUT
        if (nk == 0) return
        if (any(shape(b) /= [3, size(modes), size(points, 2), nk])) return
        allocate (g(size(points, 2), nk), gr(size(points, 2), nk), &
                  gz(size(points, 2), nk))
        call gs_circuit_flux_vectors(coils%set, points(1, :), points(2, :), g, gr, gz)
        b = (0.0_dp, 0.0_dp)
        do m = 1, size(modes)
            if (modes(m) /= 0) cycle
            do i = 1, size(points, 2)
                b(1, m, i, :) = cmplx(-gz(i, :)/points(1, i), 0.0_dp, dp)
                b(3, m, i, :) = cmplx(gr(i, :)/points(1, i), 0.0_dp, dp)
            end do
        end do
        status = GAMMA_C_OK
    end subroutine

    integer function polyline_circuits(coils) result(k)
        class(polyline_coils_t), intent(in) :: coils
        k = coils%set%ncoil
        if (allocated(coils%incidence)) k = size(coils%incidence, 2)
    end function

    subroutine polyline_flux(coils, points, g, dg_dr, dg_dz, status)
        class(polyline_coils_t), intent(in) :: coils
        real(dp), intent(in) :: points(:, :)
        real(dp), intent(out) :: g(:, :)
        real(dp), intent(out), optional :: dg_dr(:, :), dg_dz(:, :)
        integer, intent(out) :: status
        status = EQ_UNAVAILABLE
        if (size(points, 2) > size(g, 1)) return
        g = 0
        if (present(dg_dr)) dg_dr = 0
        if (present(dg_dz)) dg_dz = 0
    end subroutine

    subroutine polyline_modes(coils, points, modes, b, status)
        class(polyline_coils_t), intent(in) :: coils
        real(dp), intent(in) :: points(:, :)
        integer, intent(in) :: modes(:)
        complex(dp), intent(out) :: b(:, :, :, :)
        integer, intent(out) :: status
        real(dp), allocatable :: weights(:)
        integer :: k, nk
        nk = coils%circuits()
        status = GAMMA_C_BAD_INPUT
        if (nk == 0) return
        if (any(shape(b) /= [3, size(modes), size(points, 2), nk])) return
        allocate (weights(coils%set%ncoil))
        do k = 1, nk
            if (allocated(coils%incidence)) then
                weights = coils%incidence(:, k)
            else
                weights = 0
                weights(k) = 1
            end if
            call coil_set_modes(coils%set, weights, points, coils%nphi, modes, &
                                b(:, :, :, k), status)
            if (status /= 0) return
        end do
        status = GAMMA_C_OK
    end subroutine

    subroutine gamma_c_coil_traces(gc, coils, modes, ng, potential, normal_load, &
                                   status, closure)
        ! Per-n traces on Gamma_c of the field of the coils outside it, per
        ! unit circuit current (arrays (n_trace, size(modes), circuits)):
        !   potential    Phi_n at the trace dofs, B = grad Phi outside the
        !                coils: Phi_n = R B_phi,n/(i n) for n /= 0; for n = 0
        !                the line integral of B.t from trace dof 1 (single
        !                valued: the coils do not link Gamma_c);
        !   normal_load  int R B_n.n N_j ds, the conormal flux R d_n Phi_n
        !                tested with the trace basis (ng Gauss points per
        !                edge).
        ! closure: the n = 0 mismatch of the line integral after one turn,
        ! relative to the largest |Phi_0| (0 without an n = 0 mode).
        type(gamma_c_t), intent(in) :: gc
        class(boundary_coils_t), intent(in) :: coils
        integer, intent(in) :: modes(:), ng
        complex(dp), allocatable, intent(out) :: potential(:, :, :), normal_load(:, :, :)
        integer, intent(out) :: status
        real(dp), intent(out), optional :: closure
        real(dp), allocatable :: xq(:, :), nq(:, :), wq(:), tq(:), re(:, :), im(:, :)
        real(dp), allocatable :: cr(:, :), ci(:, :), sub(:, :), xs(:, :)
        real(dp) :: gx(ng), gw(ng), x(2), tangent(2), scale
        complex(dp), allocatable :: bq(:, :, :, :), bt(:, :, :, :), bs(:, :, :, :)
        complex(dp), allocatable :: run(:)
        integer :: nm, nk, m, i, j, q, p, point, nq_total
        status = GAMMA_C_BAD_INPUT
        if (present(closure)) closure = 0
        if (.not. allocated(gc%nodes) .or. ng < 1) return
        nm = size(modes)
        nk = coils%circuits()
        if (nk == 0) return
        allocate (potential(gc%n_trace, nm, nk), normal_load(gc%n_trace, nm, nk))
        ! Conormal loads.
        call gamma_c_quadrature(gc, ng, xq, nq, wq, tq)
        nq_total = size(wq)
        allocate (bq(3, nm, nq_total, nk))
        call coils%field_modes(xq, modes, bq, status)
        if (status /= GAMMA_C_OK) return
        allocate (re(nq_total, nk), im(nq_total, nk), cr(gc%n_trace, nk), &
                  ci(gc%n_trace, nk))
        do m = 1, nm
            do point = 1, nq_total
                re(point, :) = xq(1, point)*real(nq(1, point)*bq(1, m, point, :) &
                    + nq(2, point)*bq(3, m, point, :))
                im(point, :) = xq(1, point)*aimag(nq(1, point)*bq(1, m, point, :) &
                    + nq(2, point)*bq(3, m, point, :))
            end do
            call gamma_c_assemble(gc, tq, wq, re, cr)
            call gamma_c_assemble(gc, tq, wq, im, ci)
            normal_load(:, m, :) = cmplx(cr, ci, dp)
        end do
        ! Potential traces, n /= 0 at the trace dofs.
        allocate (bt(3, nm, gc%n_trace, nk))
        call coils%field_modes(gc%trace_points, modes, bt, status)
        if (status /= GAMMA_C_OK) return
        do m = 1, nm
            if (modes(m) == 0) cycle
            do j = 1, gc%n_trace
                potential(j, m, :) = gc%trace_points(1, j)*bt(2, m, j, :) &
                    /cmplx(0.0_dp, real(modes(m), dp), dp)
            end do
        end do
        ! n = 0: integrate B.t along each edge from its first dof to every
        ! later node j/order (ng Gauss points per sub-interval [0, j/order]).
        if (any(modes == 0)) then
            call gauss_legendre(ng, gx, gw)
            gx = 0.5_dp*(gx + 1.0_dp)
            gw = 0.5_dp*gw
            allocate (xs(2, gc%n_edges*gc%order*ng), sub(gc%n_edges*gc%order*ng, 3))
            p = 0
            do i = 1, gc%n_edges
                do j = 1, gc%order
                    do q = 1, ng
                        p = p + 1
                        call gamma_c_edge_point(gc, i, gx(q)*j/gc%order, x, tangent)
                        xs(:, p) = x
                        sub(p, 1:2) = tangent
                        sub(p, 3) = gw(q)*real(j, dp)/gc%order
                    end do
                end do
            end do
            allocate (bs(3, nm, size(xs, 2), nk), run(nk))
            call coils%field_modes(xs, modes, bs, status)
            if (status /= GAMMA_C_OK) return
            scale = 0
            do m = 1, nm
                if (modes(m) /= 0) cycle
                run = 0
                p = 0
                do i = 1, gc%n_edges
                    potential((i - 1)*gc%order + 1, m, :) = run
                    do j = 1, gc%order
                        block
                            complex(dp) :: acc(nk)
                            acc = 0
                            do q = 1, ng
                                p = p + 1
                                acc = acc + sub(p, 3)*(sub(p, 1)*bs(1, m, p, :) &
                                    + sub(p, 2)*bs(3, m, p, :))
                            end do
                            if (j < gc%order) then
                                potential((i - 1)*gc%order + j + 1, m, :) = run + acc
                            else
                                run = run + acc
                            end if
                        end block
                    end do
                end do
                scale = max(scale, maxval(abs(potential(:, m, :))))
                if (present(closure) .and. scale > 0) &
                    closure = max(closure, maxval(abs(run))/scale)
            end do
        end if
        status = GAMMA_C_OK
    end subroutine
end module kin6d_equilibrium_boundary
