module kin6d_ideal_free_boundary
    ! 3D free boundary of the nested map lane (DESIGN 7.3, core note D6, N4).
    ! One reference mesh of the region inside a fixed axisymmetric coupling
    ! curve Gamma_c, with the plasma boundary Gamma_p as an internal chain:
    !   plasma cells: the frozen-in map of kin6d_ideal_map (unchanged), with the
    !     Gamma_p nodes free in their transversal (normal) direction;
    !   shell cells (Gamma_p -> Gamma_c): a vacuum field B = grad Phi,
    !       Phi = F_vac q3 + c theta/(2 pi) + Phi~,
    !     q3 the physical toroidal angle of the shell map, c the poloidal
    !     circulation (mu0 times the plasma current), theta the reference
    !     poloidal angle about the reference axis and Phi~ single valued
    !     (Lagrange P_k x real Fourier in phi-hat, nfp*Z). The shell map is the
    !     harmonic extension of the Gamma_p displacement (mesh motion only);
    !   outside Gamma_c: per toroidal mode n the decaying exterior field enters
    !     through the exact Steklov-Poincare map S_n of FortFEM's toroidal-mode
    !     BEM, and the coils (PF/CS at n = 0, RMP at n in nfp*Z) through their
    !     traces on the fixed Gamma_c.
    ! Equilibrium is a stationary point of
    !   L = W_p[X] - (1/mu0) E_v[X, Phi~] + lambda . C(X),
    !   E_v = int_shell |grad Phi|^2/2 + sum_m c_m [ <S_n u_m, u_m>/2
    !         - <h_m, Phi~_m>_Gamma_c ],  u_m = Phi~_m - a_m on Gamma_c,
    ! with W_p the plasma energy at frozen fluxes, -E_v/mu0 the vacuum
    ! co-energy at fixed coil currents and fixed c, a_m the single-valued trace
    ! of the external potential (coils plus the unit-circulation part), h_m the
    ! external normal field and C the shell-map constraints (harmonic
    ! extension, Gamma_p tie, Gamma_c fixed) with multipliers lambda. The Phi~
    ! equations are vacuum magnetostatics with B.n = 0 on Gamma_p (natural)
    ! and continuity of B.n and B_t on Gamma_c (flux matching through S_n);
    ! the Gamma_p map equations give [[p + B^2/2mu0]] = 0 as the natural
    ! traction. The plasma position is not a [[B^2]] relaxation: it follows
    ! from the coupled Newton system, which is a saddle (vacuum block negative
    ! definite; with fixed coil currents the vertical mode is a negative
    ! direction of the reduced plasma Hessian). Newton therefore solves the
    ! KKT system with an indefinite LDL^T and globalises on a scaled residual
    ! merit, never on the energy.
    ! Everything is per field period; c_m = 2 pi/nfp for the n = 0 mode and
    ! pi/nfp for cos/sin modes.
    use, intrinsic :: iso_fortran_env, only: dp => real64, i64 => int64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use kin6d_grad_shafranov, only: gs_state_t
    use kin6d_ideal_map, only: ideal_map_t, ideal_chart_t, ideal_frozen_t, &
        ideal_map_init, ideal_map_freeze, ideal_map_energy, ideal_map_hessian, &
        ideal_map_reduce, ideal_map_expand, ideal_map_rdof, IDEAL_MAP_OK
    use kin6d_vacuum_shell_leaves, only: vacuum_energy_density, &
        vacuum_energy_gradient, vacuum_energy_hessian
    use kin6d_triangle_basis, only: triangle_basis_jet
    use kin6d_sort_keys, only: sort_key_order
    use kin6d_gs_direct, only: gs_direct_available
    use fortfem_element_to_csr, only: csr_matrix_t, csr_build_pattern, &
        csr_assemble_blocks
    use fortfem_axisymmetric_gs_boundary_operators_2d, only: gs_bem_curve_t, &
        gs_bem_operators_t, gs_bem_curve_from_nodes
    use fortfem_toroidal_mode_boundary_operators_2d, only: tm_bem_assemble, &
        tm_bem_steklov_symmetric
    use fortnum_quadrature, only: gauss_legendre
    use fortsparse, only: sparse_solver_t, sparse_factor, sparse_solve, &
        sparse_destroy, fortsparse_status_t, status_ok, FORTSPARSE_BACKEND_SUPERLU, &
        FORTSPARSE_BACKEND_MUMPS, FORTSPARSE_SYM_GENERAL, FORTSPARSE_SYM_INDEFINITE
    implicit none
    private
    public :: ideal_fb_t, ideal_fb_report_t, ideal_fb_init, ideal_fb_set_drive
    public :: ideal_fb_gradient, ideal_fb_hessian, ideal_fb_step, ideal_fb_solve
    public :: ideal_fb_shell_field, ideal_fb_sampler, ideal_fb_mode_weight
    public :: ideal_fb_solve_vacuum, ideal_fb_solve_multipliers
    integer, parameter, public :: IDEAL_FB_OK = 0, IDEAL_FB_BAD_INPUT = 1, &
        IDEAL_FB_MAP_FAILED = 2, IDEAL_FB_BEM_FAILED = 3, IDEAL_FB_FOLDED = 4, &
        IDEAL_FB_SOLVE_FAILED = 5, IDEAL_FB_NOT_CONVERGED = 6
    real(dp), parameter :: pi = acos(-1.0_dp)

    ! External field of real Fourier harmonic k (n = k*nfp) at (R, Z):
    ! B = bcos cos(n phi) + bsin sin(n phi), components (B_R, B_phi, B_Z).
    abstract interface
        subroutine ideal_fb_sampler(r, z, k, bcos, bsin)
            import :: dp
            real(dp), intent(in) :: r, z
            integer, intent(in) :: k
            real(dp), intent(out) :: bcos(3), bsin(3)
        end subroutine
    end interface

    type :: ideal_fb_t
        type(ideal_map_t) :: plasma
        integer :: nfp = 1, ntor = 0, nmode = 1, nloc = 0, order = 0
        real(dp) :: mu0 = 1, f_vac = 0, circulation = 0, theta_centre(2) = 0
        ! Shell: ns nodes, ncell cells, local connectivity, reference nodes,
        ! node role (1 Gamma_p, 2 interior, 3 Gamma_c), Gamma_p node -> plasma
        ! node, and the global mesh numbering of both parts.
        integer :: ns = 0, ncell = 0
        integer, allocatable :: conn(:, :), role(:), plasma_node(:)
        integer, allocatable :: shell_global(:), plasma_global(:)
        real(dp), allocatable :: reference(:, :)
        ! Master derivatives of theta/(2 pi) and the reference det per point.
        real(dp), allocatable :: theta_s(:, :, :), det(:, :)
        ! Reference Laplacian of the shell (harmonic extension), CSR.
        integer, allocatable :: lap_ptr(:), lap_col(:)
        real(dp), allocatable :: lap_val(:)
        ! Gamma_c: nt trace dofs (shell nodes), ccw edges (0:order, nedge),
        ! Galerkin S per harmonic k (nt, nt, 0:ntor), external trace a and
        ! normal-field load h per real mode (nt, nmode).
        integer :: nt = 0, nedge = 0
        integer, allocatable :: trace(:), edge(:, :)
        real(dp), allocatable :: s_trace(:, :, :), drive(:, :), load(:, :)
        ! State: shell map coefficients (3, nmode, ns), potential (nmode, ns),
        ! multipliers (3, nmode, ns).
        real(dp), allocatable :: xs(:, :, :), phi(:, :), lam(:, :, :)
        ! Unknown layout: plasma reduced, shell map, multipliers, potential.
        integer :: np = 0, nx = 0, nf = 0, n = 0
        ! Energy scale sigma*V (J per field period) for the reported residual.
        real(dp) :: energy_scale = 1
    end type

    type :: ideal_fb_report_t
        integer :: iterations = 0
        real(dp) :: merit0 = 0, merit = 0, residual = huge(1.0_dp)
        real(dp) :: step = huge(1.0_dp), boundary_step = huge(1.0_dp)
        real(dp) :: constraint = 0, min_jacobian = 0, damping = 0
        real(dp) :: assemble_s = 0, factor_s = 0
    end type
contains
    pure real(dp) function ideal_fb_mode_weight(fb, m) result(c)
        ! Field-period integral of the squared real Fourier mode m.
        type(ideal_fb_t), intent(in) :: fb
        integer, intent(in) :: m
        c = merge(2*pi, pi, m == 1)/fb%nfp
    end function

    pure integer function harmonic(m)
        integer, intent(in) :: m
        harmonic = m/2
    end function

    pure integer function ix(fb, c, m, s)
        type(ideal_fb_t), intent(in) :: fb
        integer, intent(in) :: c, m, s
        ix = fb%np + c + 3*(m - 1) + 3*fb%nmode*(s - 1)
    end function

    pure integer function il(fb, c, m, s)
        type(ideal_fb_t), intent(in) :: fb
        integer, intent(in) :: c, m, s
        il = fb%np + fb%nx + c + 3*(m - 1) + 3*fb%nmode*(s - 1)
    end function

    pure integer function iphi(fb, m, s)
        type(ideal_fb_t), intent(in) :: fb
        integer, intent(in) :: m, s
        iphi = fb%np + 2*fb%nx + m + fb%nmode*(s - 1)
    end function

    subroutine ideal_fb_init(fb, mesh, chain, chart, field, nfp, ntor, f_vac, &
                             circulation, status, mu0, nphi, panels_per_edge, &
                             quadrature_order, omega0_free)
        ! mesh: curved P2/P3 mesh of the Gamma_c region with Gamma_p as the
        ! constrained internal chain (initialize_gs with interior_chain);
        ! chain: the Gamma_p vertices used for it; field: frozen plasma data.
        type(ideal_fb_t), intent(out) :: fb
        type(gs_state_t), intent(in) :: mesh
        real(dp), intent(in) :: chain(:, :), f_vac, circulation
        class(ideal_chart_t), intent(in) :: chart
        class(ideal_frozen_t), intent(in) :: field
        integer, intent(in) :: nfp, ntor
        integer, intent(out) :: status
        real(dp), intent(in), optional :: mu0
        integer, intent(in), optional :: nphi, panels_per_edge, quadrature_order
        logical, intent(in), optional :: omega0_free
        type(gs_state_t) :: sub
        logical, allocatable :: inside(:), in_plasma(:), in_shell(:)
        integer :: e, a, m, nodes, cells, nplasma, k, st, ppe
        real(dp) :: centre(2)
        status = IDEAL_FB_BAD_INPUT
        fb%nloc = size(mesh%triangles, 1)
        if (fb%nloc /= 6 .and. fb%nloc /= 10) return
        fb%order = merge(2, 3, fb%nloc == 6)
        if (.not. allocated(mesh%neighbors)) return
        nodes = size(mesh%points, 2); cells = size(mesh%triangles, 2)
        fb%nfp = nfp; fb%ntor = ntor; fb%nmode = 2*ntor + 1
        fb%f_vac = f_vac; fb%circulation = circulation
        if (present(mu0)) fb%mu0 = mu0
        ppe = merge(1, 2, fb%order == 2)
        if (present(panels_per_edge)) ppe = panels_per_edge
        ! Plasma cells: centroid inside the Gamma_p polygon.
        allocate (inside(cells), in_plasma(nodes), in_shell(nodes))
        do e = 1, cells
            centre = sum(mesh%points(:, mesh%triangles(1:3, e)), dim=2)/3
            inside(e) = point_in_polygon(chain, centre)
        end do
        in_plasma = .false.; in_shell = .false.
        do e = 1, cells
            if (inside(e)) then
                in_plasma(mesh%triangles(:, e)) = .true.
            else
                in_shell(mesh%triangles(:, e)) = .true.
            end if
        end do
        if (count(inside) == 0 .or. count(.not. inside) == 0) return
        if (any(mesh%boundary .and. in_plasma)) return
        ! Plasma sub-mesh in the form ideal_map_init reads.
        nplasma = count(in_plasma)
        allocate (fb%plasma_global(nplasma), sub%points(2, nplasma), &
            sub%boundary(nplasma), sub%triangles(fb%nloc, count(inside)))
        block
            integer, allocatable :: number(:)
            allocate (number(nodes))
            number = 0; k = 0
            do a = 1, nodes
                if (.not. in_plasma(a)) cycle
                k = k + 1; number(a) = k; fb%plasma_global(k) = a
                sub%points(:, k) = mesh%points(:, a)
                sub%boundary(k) = in_shell(a)
            end do
            k = 0
            do e = 1, cells
                if (.not. inside(e)) cycle
                k = k + 1
                sub%triangles(:, k) = number(mesh%triangles(:, e))
            end do
            ! Shell numbering.
            fb%ns = count(in_shell); fb%ncell = cells - count(inside)
            allocate (fb%shell_global(fb%ns), fb%reference(2, fb%ns), &
                fb%role(fb%ns), fb%plasma_node(fb%ns), fb%conn(fb%nloc, fb%ncell))
            k = 0
            do a = 1, nodes
                if (.not. in_shell(a)) cycle
                k = k + 1; fb%shell_global(k) = a
                fb%reference(:, k) = mesh%points(:, a)
                fb%role(k) = 2
                fb%plasma_node(k) = number(a)
                if (in_plasma(a)) fb%role(k) = 1
                if (mesh%boundary(a)) fb%role(k) = 3
            end do
            number = 0
            do k = 1, fb%ns
                number(fb%shell_global(k)) = k
            end do
            k = 0
            do e = 1, cells
                if (inside(e)) cycle
                k = k + 1
                fb%conn(:, k) = number(mesh%triangles(:, e))
            end do
            call gamma_c_chain(fb, mesh, number, st)
            if (st /= 0) return
        end block
        status = IDEAL_FB_MAP_FAILED
        call ideal_map_init(fb%plasma, sub, chart, nfp, ntor, st, nphi, &
            quadrature_order, mu0)
        if (st /= IDEAL_MAP_OK) return
        ! Free boundary: Gamma_p nodes move along their transversal.
        do a = 1, nplasma
            if (.not. sub%boundary(a)) cycle
            do m = 1, fb%nmode
                fb%plasma%fixed(ideal_map_rdof(fb%plasma, 1, m, a)) = .false.
            end do
        end do
        ! omega_0 (the n = 0 toroidal displacement) is held at zero. Its
        ! surface-function part omega_0 = f(chi) is the continuum relabelling
        ! null of the map lane; with frozen data sampled from a finite-element
        ! flux (gradients only C0) the weak N1 gauge term does not hold it and
        ! Newton drifts along it, dragging the surfaces (measured: 1-5 cm
        ! interior surface defect on the gs_free fixture, fixed or free
        ! boundary). Its non-f(chi) part twists field lines and is zero at
        ! any axisymmetric equilibrium. The cost is the O(eps^2) n = 0
        ! toroidal shear of a 3D state, outside the linear and NL-3a gates.
        if (present(omega0_free)) then
            if (.not. omega0_free) call pin_omega0()
        else
            call pin_omega0()
        end if
        call ideal_map_freeze(fb%plasma, field, st)
        if (st /= IDEAL_MAP_OK) return
        fb%mu0 = fb%plasma%mu0
        ! Centre of the poloidal-angle lift: the reference axis node, or the
        ! interior plasma node nearest the plasma centroid when the frozen
        ! labels put the extremum on Gamma_p.
        fb%theta_centre = fb%plasma%reference(:, fb%plasma%axis_node)
        if (sub%boundary(fb%plasma%axis_node)) then
            centre = sum(sub%points, dim=2)/nplasma
            k = minloc(norm2(sub%points - spread(centre, 2, nplasma), dim=1), dim=1, &
                mask=.not. sub%boundary)
            fb%theta_centre = sub%points(:, k)
        end if
        call shell_geometry(fb)
        call shell_laplacian(fb)
        status = IDEAL_FB_BEM_FAILED
        call gamma_c_operators(fb, ppe, st)
        if (st /= 0) return
        fb%np = fb%plasma%nred; fb%nx = 3*fb%nmode*fb%ns; fb%nf = fb%nmode*fb%ns
        fb%n = fb%np + 2*fb%nx + fb%nf
        allocate (fb%xs(3, fb%nmode, fb%ns), fb%phi(fb%nmode, fb%ns), &
            fb%lam(3, fb%nmode, fb%ns), fb%drive(fb%nt, fb%nmode), &
            fb%load(fb%nt, fb%nmode))
        fb%xs = 0; fb%phi = 0; fb%lam = 0; fb%drive = 0; fb%load = 0
        fb%xs(1:2, 1, :) = fb%reference
        do k = 1, fb%ns
            if (fb%role(k) == 1) fb%xs(:, :, k) = fb%plasma%x(:, :, fb%plasma_node(k))
        end do
        status = IDEAL_FB_OK
    contains
        subroutine pin_omega0()
            integer :: node
            do node = 1, fb%plasma%nodes
                fb%plasma%fixed(ideal_map_rdof(fb%plasma, 2, 1, node)) = .true.
            end do
        end subroutine
    end subroutine

    pure logical function point_in_polygon(poly, p) result(inside)
        real(dp), intent(in) :: poly(:, :), p(2)
        integer :: i, j, n
        n = size(poly, 2)
        inside = .false.
        j = n
        do i = 1, n
            if ((poly(2, i) > p(2)) .neqv. (poly(2, j) > p(2))) then
                if (p(1) < (poly(1, j) - poly(1, i))*(p(2) - poly(2, i)) &
                    /(poly(2, j) - poly(2, i)) + poly(1, i)) inside = .not. inside
            end if
            j = i
        end do
    end function

    subroutine gamma_c_chain(fb, mesh, number, status)
        ! Counter-clockwise chain of the outer boundary edges; trace dofs are
        ! the edge nodes 0..order-1 of each edge, as shell nodes.
        type(ideal_fb_t), intent(inout) :: fb
        type(gs_state_t), intent(in) :: mesh
        integer, intent(in) :: number(:)
        integer, intent(out) :: status
        integer, allocatable :: first(:), ee(:), ek(:)
        integer :: e, k, i, nb, nxt, j
        status = 1
        allocate (first(size(mesh%points, 2)))
        first = 0; nb = 0
        do e = 1, size(mesh%triangles, 2)
            do k = 1, 3
                if (mesh%neighbors(k, e) > 0) cycle
                if (first(mesh%triangles(k, e)) /= 0) return
                first(mesh%triangles(k, e)) = 3*(e - 1) + k
                nb = nb + 1
            end do
        end do
        if (nb < 4) return
        allocate (ee(nb), ek(nb))
        do i = 1, size(first)
            if (first(i) /= 0) exit
        end do
        ee(1) = (first(i) - 1)/3 + 1; ek(1) = first(i) - 3*(ee(1) - 1)
        do i = 2, nb
            nxt = first(mesh%triangles(mod(ek(i - 1), 3) + 1, ee(i - 1)))
            if (nxt == 0) return
            ee(i) = (nxt - 1)/3 + 1; ek(i) = nxt - 3*(ee(i) - 1)
        end do
        if (first(mesh%triangles(mod(ek(nb), 3) + 1, ee(nb))) &
            /= 3*(ee(1) - 1) + ek(1)) return
        fb%nedge = nb; fb%nt = nb*fb%order
        allocate (fb%edge(0:fb%order, nb), fb%trace(fb%nt))
        do i = 1, nb
            fb%edge(0, i) = number(mesh%triangles(ek(i), ee(i)))
            do j = 1, fb%order - 1
                fb%edge(j, i) = number(mesh%triangles(3 + (ek(i) - 1)*(fb%order - 1) &
                    + j, ee(i)))
            end do
            fb%edge(fb%order, i) = number(mesh%triangles(mod(ek(i), 3) + 1, ee(i)))
            fb%trace((i - 1)*fb%order + 1:i*fb%order) = fb%edge(0:fb%order - 1, i)
        end do
        if (any(fb%trace == 0)) return
        status = 0
    end subroutine

    pure subroutine edge_lagrange(order, t, w, dw)
        integer, intent(in) :: order
        real(dp), intent(in) :: t
        real(dp), intent(out) :: w(0:order), dw(0:order)
        real(dp) :: prod
        integer :: j, l, m
        do j = 0, order
            w(j) = 1; dw(j) = 0
            do l = 0, order
                if (l == j) cycle
                w(j) = w(j)*(t - real(l, dp)/order)/(real(j - l, dp)/order)
            end do
            do m = 0, order
                if (m == j) cycle
                prod = 1/(real(j - m, dp)/order)
                do l = 0, order
                    if (l == j .or. l == m) cycle
                    prod = prod*(t - real(l, dp)/order)/(real(j - l, dp)/order)
                end do
                dw(j) = dw(j) + prod
            end do
        end do
    end subroutine

    subroutine edge_point(fb, i, t, x, tangent, w)
        ! Position, tangent dx/dt and edge basis of Gamma_c edge i at t.
        type(ideal_fb_t), intent(in) :: fb
        integer, intent(in) :: i
        real(dp), intent(in) :: t
        real(dp), intent(out) :: x(2), tangent(2), w(0:fb%order)
        real(dp) :: dw(0:3)
        integer :: j
        call edge_lagrange(fb%order, t, w, dw(0:fb%order))
        x = 0; tangent = 0
        do j = 0, fb%order
            x = x + w(j)*fb%reference(:, fb%edge(j, i))
            tangent = tangent + dw(j)*fb%reference(:, fb%edge(j, i))
        end do
    end subroutine

    subroutine gamma_c_operators(fb, ppe, status)
        ! S_trace(:, :, k) = R^T S_h(n = k nfp) R with R the restriction of
        ! the P_k FEM trace to the P2 panel trace (ppe panels per edge).
        type(ideal_fb_t), intent(inout) :: fb
        integer, intent(in) :: ppe
        integer, intent(out) :: status
        type(gs_bem_curve_t) :: curve
        type(gs_bem_operators_t) :: ops
        real(dp), allocatable :: nodes(:, :), restriction(:, :), s_bem(:, :)
        real(dp) :: w(0:3), tangent(2), t
        integer :: i, q, j, k, p, row, col, npan, st
        status = 1
        npan = fb%nedge*ppe
        allocate (nodes(2, 2*npan), restriction(2*npan, fb%nt))
        restriction = 0
        do i = 1, fb%nedge
            do q = 1, ppe
                p = (i - 1)*ppe + q
                do j = 1, 2
                    t = (real(q - 1, dp) + 0.5_dp*(j - 1))/ppe
                    row = 2*p - 2 + j
                    call edge_point(fb, i, t, nodes(:, row), tangent, w(0:fb%order))
                    do k = 0, fb%order
                        col = (i - 1)*fb%order + k + 1
                        if (k == fb%order) col = mod(i, fb%nedge)*fb%order + 1
                        restriction(row, col) = restriction(row, col) + w(k)
                    end do
                end do
            end do
        end do
        call gs_bem_curve_from_nodes(nodes, curve, st)
        if (st /= 0) return
        allocate (fb%s_trace(fb%nt, fb%nt, 0:fb%ntor))
        do k = 0, fb%ntor
            call tm_bem_assemble(curve, k*fb%nfp, 2, ops, st)
            if (st /= 0) return
            call tm_bem_steklov_symmetric(ops, s_bem, st)
            if (st /= 0) return
            fb%s_trace(:, :, k) = matmul(transpose(restriction), &
                matmul(s_bem, restriction))
            fb%s_trace(:, :, k) = 0.5_dp*(fb%s_trace(:, :, k) &
                + transpose(fb%s_trace(:, :, k)))
        end do
        status = 0
    end subroutine

    subroutine ideal_fb_set_drive(fb, sample, status, closure)
        ! External traces from the field of everything outside Gamma_c plus
        ! the unit-circulation exterior part: drive a_m (single-valued
        ! external potential on the trace nodes) and load h_m = int (B.n)
        ! N_i R ds. n = 0: a = int (B - c grad theta/2pi).t ds from trace
        ! node 1; closure returns its mismatch after one turn (T m). n /= 0:
        ! a from B_phi = (1/R) dPhi/dphi; closure also returns the largest
        ! mismatch between d a/ds and B.t (curl-free check, T).
        type(ideal_fb_t), intent(inout) :: fb
        procedure(ideal_fb_sampler) :: sample
        integer, intent(out) :: status
        real(dp), intent(out), optional :: closure
        integer, parameter :: ng = 10
        real(dp) :: gx(ng), gw(ng), x(2), tangent(2), w(0:3), bc(3), bs(3)
        real(dp) :: grad_t(2), acc, run, speed, worst, nn, at(2), dmax
        integer :: i, q, j, k, m, col
        status = IDEAL_FB_OK
        call gauss_legendre(ng, gx, gw)
        gx = 0.5_dp*(gx + 1); gw = 0.5_dp*gw
        fb%drive = 0; fb%load = 0; worst = 0
        ! Loads.
        do i = 1, fb%nedge
            do q = 1, ng
                call edge_point(fb, i, gx(q), x, tangent, w(0:fb%order))
                do k = 0, fb%ntor
                    call sample(x(1), x(2), k, bc, bs)
                    do j = 0, fb%order
                        col = (i - 1)*fb%order + j + 1
                        if (j == fb%order) col = mod(i, fb%nedge)*fb%order + 1
                        if (k == 0) then
                            fb%load(col, 1) = fb%load(col, 1) + gw(q)*w(j)*x(1) &
                                *(bc(1)*tangent(2) - bc(3)*tangent(1))
                        else
                            fb%load(col, 2*k) = fb%load(col, 2*k) + gw(q)*w(j)*x(1) &
                                *(bc(1)*tangent(2) - bc(3)*tangent(1))
                            fb%load(col, 2*k + 1) = fb%load(col, 2*k + 1) &
                                + gw(q)*w(j)*x(1)*(bs(1)*tangent(2) - bs(3)*tangent(1))
                        end if
                    end do
                end do
            end do
        end do
        ! n = 0 trace by integration along Gamma_c.
        run = 0
        do i = 1, fb%nedge
            fb%drive((i - 1)*fb%order + 1, 1) = run
            do j = 1, fb%order
                acc = 0
                do q = 1, ng
                    call edge_point(fb, i, gx(q)*j/fb%order, x, tangent, w(0:fb%order))
                    call sample(x(1), x(2), 0, bc, bs)
                    grad_t = theta_gradient(fb, x)
                    acc = acc + gw(q)*real(j, dp)/fb%order*(dot_product(bc([1, 3]) &
                        - fb%circulation*grad_t, tangent))
                end do
                if (j < fb%order) then
                    fb%drive((i - 1)*fb%order + j + 1, 1) = run + acc
                else
                    run = run + acc
                end if
            end do
        end do
        worst = abs(run)
        ! n /= 0 traces from B_phi at the trace nodes, and the curl-free check.
        do k = 1, fb%ntor
            nn = k*fb%nfp
            do i = 1, fb%nedge
                do j = 0, fb%order - 1
                    col = (i - 1)*fb%order + j + 1
                    x = fb%reference(:, fb%trace(col))
                    call sample(x(1), x(2), k, bc, bs)
                    fb%drive(col, 2*k) = -x(1)*bs(2)/nn
                    fb%drive(col, 2*k + 1) = x(1)*bc(2)/nn
                end do
            end do
            dmax = 0
            do i = 1, fb%nedge
                call edge_point(fb, i, 0.5_dp, x, tangent, w(0:fb%order))
                call sample(x(1), x(2), k, bc, bs)
                speed = norm2(tangent)
                block
                    real(dp) :: dw(0:3), wv(0:3)
                    call edge_lagrange(fb%order, 0.5_dp, wv(0:fb%order), dw(0:fb%order))
                    at = 0
                    do j = 0, fb%order
                        col = (i - 1)*fb%order + j + 1
                        if (j == fb%order) col = mod(i, fb%nedge)*fb%order + 1
                        at = at + dw(j)*fb%drive(col, 2*k:2*k + 1)
                    end do
                end block
                dmax = max(dmax, abs(at(1)/speed - dot_product(bc([1, 3]), &
                    tangent)/speed), abs(at(2)/speed - dot_product(bs([1, 3]), &
                    tangent)/speed))
            end do
            worst = max(worst, dmax)
        end do
        do m = 1, fb%nmode
            fb%load(:, m) = fb%load(:, m)
        end do
        if (present(closure)) closure = worst
    end subroutine

    pure function theta_gradient(fb, y) result(g)
        ! grad (theta/2 pi) about the reference axis.
        type(ideal_fb_t), intent(in) :: fb
        real(dp), intent(in) :: y(2)
        real(dp) :: g(2), d(2)
        d = y - fb%theta_centre
        g = [-d(2), d(1)]/(2*pi*sum(d**2))
    end function

    subroutine shell_geometry(fb)
        type(ideal_fb_t), intent(inout) :: fb
        real(dp) :: dg(2, 2), y(2)
        integer :: e, g, nq
        nq = fb%plasma%nq
        allocate (fb%theta_s(2, nq, fb%ncell), fb%det(nq, fb%ncell))
        do e = 1, fb%ncell
            do g = 1, nq
                y = matmul(fb%reference(:, fb%conn(:, e)), fb%plasma%shape(:, g))
                dg = matmul(fb%reference(:, fb%conn(:, e)), &
                    transpose(fb%plasma%dshape(:, :, g)))
                fb%det(g, e) = dg(1, 1)*dg(2, 2) - dg(1, 2)*dg(2, 1)
                fb%theta_s(:, g, e) = matmul(transpose(dg), theta_gradient(fb, y))
            end do
        end do
    end subroutine

    subroutine shell_laplacian(fb)
        ! Reference stiffness int grad N_i . grad N_j dy over the shell.
        type(ideal_fb_t), intent(inout) :: fb
        integer(i64), allocatable :: keys(:)
        real(dp), allocatable :: raw(:)
        integer, allocatable :: order(:), rows(:)
        real(dp) :: dg(2, 2), ginv(2, 2), gr(2, 10), k(10, 10)
        integer :: e, g, a, b, cnt, i, nnz, nl
        integer(i64) :: key, previous
        nl = fb%nloc
        allocate (keys(fb%ncell*nl*nl), raw(fb%ncell*nl*nl), order(fb%ncell*nl*nl))
        cnt = 0
        do e = 1, fb%ncell
            k = 0
            do g = 1, fb%plasma%nq
                dg = matmul(fb%reference(:, fb%conn(:, e)), &
                    transpose(fb%plasma%dshape(:, :, g)))
                ginv = reshape([dg(2, 2), -dg(2, 1), -dg(1, 2), dg(1, 1)], [2, 2]) &
                    /fb%det(g, e)
                gr(:, :nl) = matmul(transpose(ginv), fb%plasma%dshape(:, :, g))
                k(:nl, :nl) = k(:nl, :nl) + fb%plasma%weight(g)*fb%det(g, e) &
                    *matmul(transpose(gr(:, :nl)), gr(:, :nl))
            end do
            do b = 1, nl
                do a = 1, nl
                    cnt = cnt + 1
                    keys(cnt) = int(fb%conn(a, e) - 1, i64)*fb%ns + fb%conn(b, e) - 1
                    raw(cnt) = k(a, b)
                end do
            end do
        end do
        call sort_key_order(keys, order)
        allocate (fb%lap_col(cnt), fb%lap_val(cnt), rows(cnt))
        nnz = 0; previous = -1
        do i = 1, cnt
            key = keys(order(i))
            if (key /= previous) then
                nnz = nnz + 1
                rows(nnz) = int(key/fb%ns) + 1
                fb%lap_col(nnz) = int(mod(key, int(fb%ns, i64))) + 1
                fb%lap_val(nnz) = 0
                previous = key
            end if
            fb%lap_val(nnz) = fb%lap_val(nnz) + raw(order(i))
        end do
        fb%lap_col = fb%lap_col(:nnz); fb%lap_val = fb%lap_val(:nnz)
        allocate (fb%lap_ptr(fb%ns + 1))
        fb%lap_ptr = 0
        do i = 1, nnz
            fb%lap_ptr(rows(i) + 1) = fb%lap_ptr(rows(i) + 1) + 1
        end do
        fb%lap_ptr(1) = 1
        do i = 1, fb%ns
            fb%lap_ptr(i + 1) = fb%lap_ptr(i + 1) + fb%lap_ptr(i)
        end do
    end subroutine

    subroutine shell_jet(fb, xs, phi, e, g, j, q, qd, gt)
        ! Chart coordinates q, Q(k,i) = dq_k/ds_i and dPhi~/ds at a shell point.
        type(ideal_fb_t), intent(in) :: fb
        real(dp), intent(in) :: xs(:, :, :), phi(:, :)
        integer, intent(in) :: e, g, j
        real(dp), intent(out) :: q(3), qd(3, 3), gt(3)
        real(dp) :: t, dt
        integer :: a, m, node
        q = 0; qd = 0; gt = 0
        do a = 1, fb%nloc
            node = fb%conn(a, e)
            do m = 1, fb%nmode
                t = fb%plasma%trig(m, j); dt = fb%plasma%dtrig(m, j)
                q = q + xs(:, m, node)*fb%plasma%shape(a, g)*t
                qd(:, 1) = qd(:, 1) + xs(:, m, node)*fb%plasma%dshape(1, a, g)*t
                qd(:, 2) = qd(:, 2) + xs(:, m, node)*fb%plasma%dshape(2, a, g)*t
                qd(:, 3) = qd(:, 3) + xs(:, m, node)*fb%plasma%shape(a, g)*dt
                gt(1) = gt(1) + phi(m, node)*fb%plasma%dshape(1, a, g)*t
                gt(2) = gt(2) + phi(m, node)*fb%plasma%dshape(2, a, g)*t
                gt(3) = gt(3) + phi(m, node)*fb%plasma%shape(a, g)*dt
            end do
        end do
        q(3) = q(3) + fb%plasma%phi(j)
        qd(3, 3) = qd(3, 3) + 1
    end subroutine

    subroutine shell_point(fb, xs, phi, e, g, j, w, jac, gv, hv, want_hessian)
        ! Vacuum energy density and its derivatives with respect to
        ! v = (q(1:3), Q(k,i) at 3 + k + 3(i-1), dPhi~/ds_i at 12 + i).
        type(ideal_fb_t), intent(in) :: fb
        real(dp), intent(in) :: xs(:, :, :), phi(:, :)
        integer, intent(in) :: e, g, j
        real(dp), intent(out) :: w, jac, gv(15), hv(15, 15)
        logical, intent(in) :: want_hessian
        real(dp) :: q(3), qd(3, 3), gt(3), gg(3), pos(3), d1(3, 3), d2(3, 3, 3)
        real(dp) :: d3(3, 3, 3, 3), dx(3, 3), gd(12), hd(12, 12), dv(12, 15)
        real(dp) :: hflat(78), gm(3, 3)
        integer :: a, i, k, m, n, r, c
        call shell_jet(fb, xs, phi, e, g, j, q, qd, gt)
        call fb%plasma%chart%jet(q, pos, d1, d2, d3)
        dx = matmul(d1, qd)
        gg = gt + fb%f_vac*qd(3, :) + fb%circulation*[fb%theta_s(:, g, e), 0.0_dp]
        call vacuum_energy_density(dx(1, 1), dx(2, 1), dx(3, 1), dx(1, 2), &
            dx(2, 2), dx(3, 2), dx(1, 3), dx(2, 3), dx(3, 3), gg(1), gg(2), gg(3), &
            fb%mu0, fb%plasma%sigma, w, jac)
        call vacuum_energy_gradient(dx(1, 1), dx(2, 1), dx(3, 1), dx(1, 2), &
            dx(2, 2), dx(3, 2), dx(1, 3), dx(2, 3), dx(3, 3), gg(1), gg(2), gg(3), &
            fb%mu0, fb%plasma%sigma, gd(1), gd(2), gd(3), gd(4), gd(5), gd(6), &
            gd(7), gd(8), gd(9), gd(10), gd(11), gd(12))
        dv = 0
        do i = 1, 3
            do a = 1, 3
                r = a + 3*(i - 1)
                do m = 1, 3
                    dv(r, m) = dot_product(d2(a, :, m), qd(:, i))
                end do
                do k = 1, 3
                    dv(r, 3 + k + 3*(i - 1)) = d1(a, k)
                end do
            end do
            dv(9 + i, 6 + 3*(i - 1)) = fb%f_vac
            dv(9 + i, 12 + i) = 1
        end do
        gv = matmul(gd, dv)
        hv = 0
        if (.not. want_hessian) return
        call hessian_leaf(hflat)
        n = 0
        do c = 1, 12
            do r = 1, c
                n = n + 1
                hd(r, c) = hflat(n); hd(c, r) = hflat(n)
            end do
        end do
        hv = matmul(transpose(dv), matmul(hd, dv))
        gm = reshape(gd(1:9), [3, 3])
        do m = 1, 3
            do n = 1, 3
                do i = 1, 3
                    do a = 1, 3
                        hv(m, n) = hv(m, n) + gm(a, i)*dot_product(d3(a, :, m, n), &
                            qd(:, i))
                    end do
                end do
            end do
            do i = 1, 3
                do k = 1, 3
                    r = 3 + k + 3*(i - 1)
                    hv(m, r) = hv(m, r) + dot_product(gm(:, i), d2(:, k, m))
                    hv(r, m) = hv(m, r)
                end do
            end do
        end do
    contains
        subroutine hessian_leaf(t)
            real(dp), intent(out) :: t(78)
            call vacuum_energy_hessian(dx(1, 1), dx(2, 1), dx(3, 1), dx(1, 2), &
                dx(2, 2), dx(3, 2), dx(1, 3), dx(2, 3), dx(3, 3), gg(1), gg(2), &
                gg(3), fb%mu0, fb%plasma%sigma, t(1), t(2), t(3), t(4), t(5), t(6), &
                t(7), t(8), t(9), t(10), t(11), t(12), t(13), t(14), t(15), t(16), &
                t(17), t(18), t(19), t(20), t(21), t(22), t(23), t(24), t(25), &
                t(26), t(27), t(28), t(29), t(30), t(31), t(32), t(33), t(34), &
                t(35), t(36), t(37), t(38), t(39), t(40), t(41), t(42), t(43), &
                t(44), t(45), t(46), t(47), t(48), t(49), t(50), t(51), t(52), &
                t(53), t(54), t(55), t(56), t(57), t(58), t(59), t(60), t(61), &
                t(62), t(63), t(64), t(65), t(66), t(67), t(68), t(69), t(70), &
                t(71), t(72), t(73), t(74), t(75), t(76), t(77), t(78))
        end subroutine
    end subroutine

    subroutine local_rows(fb, g, j, rows)
        ! rows(v, l): dv/dU_l for local dof l = (component c = 1..4, mode m,
        ! node a) at 4*(m-1) + c + 4*nmode*(a-1); c = 4 is Phi~.
        type(ideal_fb_t), intent(in) :: fb
        integer, intent(in) :: g, j
        real(dp), intent(out) :: rows(:, :)
        real(dp) :: f(0:3)
        integer :: a, m, c, l
        rows = 0
        do a = 1, fb%nloc
            do m = 1, fb%nmode
                f(0) = fb%plasma%shape(a, g)*fb%plasma%trig(m, j)
                f(1) = fb%plasma%dshape(1, a, g)*fb%plasma%trig(m, j)
                f(2) = fb%plasma%dshape(2, a, g)*fb%plasma%trig(m, j)
                f(3) = fb%plasma%shape(a, g)*fb%plasma%dtrig(m, j)
                l = 4*(m - 1) + 4*fb%nmode*(a - 1)
                do c = 1, 3
                    rows(c, l + c) = f(0)
                    rows(3 + c, l + c) = f(1)
                    rows(6 + c, l + c) = f(2)
                    rows(9 + c, l + c) = f(3)
                end do
                rows(13:15, l + 4) = f(1:3)
            end do
        end do
    end subroutine

    subroutine shell_energy(fb, xs, phi, energy, gx, gf, min_jacobian, blocks)
        ! E_v interior part (with 1/mu0), gradients with respect to xs and
        ! Phi~, and optionally the local Hessian blocks (shell-local order).
        type(ideal_fb_t), intent(in) :: fb
        real(dp), intent(in) :: xs(:, :, :), phi(:, :)
        real(dp), intent(out) :: energy, gx(:, :, :), gf(:, :), min_jacobian
        real(dp), intent(out), optional :: blocks(:, :, :)
        real(dp), allocatable :: local(:, :), rows(:, :), tmp(:, :)
        real(dp) :: w, jac, gv(15), hv(15, 15), weight, dphi, minjac
        integer :: e, g, j, nl, a, m, l
        logical :: want
        want = present(blocks)
        nl = 4*fb%nmode*fb%nloc
        allocate (local(nl, fb%ncell))
        dphi = 2*pi/(fb%nfp*fb%plasma%nphi)
        energy = 0; minjac = huge(1.0_dp)
        !$omp parallel default(shared) private(e, g, j, w, jac, gv, hv, weight, rows, tmp)
        allocate (rows(15, nl), tmp(15, nl))
        !$omp do reduction(+: energy) reduction(min: minjac) schedule(static)
        do e = 1, fb%ncell
            local(:, e) = 0
            if (want) blocks(:, :, e) = 0
            do j = 1, fb%plasma%nphi
                do g = 1, fb%plasma%nq
                    call shell_point(fb, xs, phi, e, g, j, w, jac, gv, hv, want)
                    minjac = min(minjac, jac/fb%det(g, e))
                    weight = fb%plasma%weight(g)*dphi
                    energy = energy + weight*w
                    call local_rows(fb, g, j, rows)
                    local(:, e) = local(:, e) + weight*matmul(gv, rows)
                    if (want) then
                        tmp = matmul(hv, rows)
                        blocks(:, :, e) = blocks(:, :, e) + weight &
                            *matmul(transpose(rows), tmp)
                    end if
                end do
            end do
        end do
        !$omp end do
        deallocate (rows, tmp)
        !$omp end parallel
        min_jacobian = minjac
        gx = 0; gf = 0
        do e = 1, fb%ncell
            do a = 1, fb%nloc
                do m = 1, fb%nmode
                    l = 4*(m - 1) + 4*fb%nmode*(a - 1)
                    gx(:, m, fb%conn(a, e)) = gx(:, m, fb%conn(a, e)) &
                        + local(l + 1:l + 3, e)
                    gf(m, fb%conn(a, e)) = gf(m, fb%conn(a, e)) + local(l + 4, e)
                end do
            end do
        end do
    end subroutine

    subroutine constraint_residual(fb, res)
        type(ideal_fb_t), intent(in) :: fb
        real(dp), intent(out) :: res(:, :, :)
        real(dp), allocatable :: u(:, :, :)
        integer :: s, p
        u = fb%xs
        u(1:2, 1, :) = u(1:2, 1, :) - fb%reference
        do s = 1, fb%ns
            select case (fb%role(s))
            case (1)
                res(:, :, s) = fb%xs(:, :, s) - fb%plasma%x(:, :, fb%plasma_node(s))
            case (3)
                res(:, :, s) = u(:, :, s)
            case default
                res(:, :, s) = 0
                do p = fb%lap_ptr(s), fb%lap_ptr(s + 1) - 1
                    res(:, :, s) = res(:, :, s) + fb%lap_val(p)*u(:, :, fb%lap_col(p))
                end do
            end select
        end do
    end subroutine

    subroutine ideal_fb_gradient(fb, grad, status, lagrangian, min_jacobian)
        ! Gradient of L with respect to (plasma reduced, xs, lambda, Phi~);
        ! fixed plasma dofs have zero rows.
        type(ideal_fb_t), intent(in) :: fb
        real(dp), intent(out) :: grad(:)
        integer, intent(out) :: status
        real(dp), intent(out), optional :: lagrangian, min_jacobian
        real(dp), allocatable :: gfull(:), gx(:, :, :), gf(:, :), res(:, :, :), u(:)
        real(dp), allocatable :: cx(:, :, :)
        real(dp) :: wp, ev, ec, minp, mins, cm, d(2)
        integer :: st, s, m, c, p, a, k
        allocate (gfull(fb%plasma%ndof), gx(3, fb%nmode, fb%ns), gf(fb%nmode, fb%ns), &
            res(3, fb%nmode, fb%ns), cx(3, fb%nmode, fb%ns))
        status = IDEAL_FB_OK
        call ideal_map_energy(fb%plasma, fb%plasma%x, wp, st, gfull, minp)
        if (st /= IDEAL_MAP_OK) status = IDEAL_FB_FOLDED
        call ideal_map_reduce(fb%plasma, gfull, grad(1:fb%np))
        call shell_energy(fb, fb%xs, fb%phi, ev, gx, gf, mins)
        if (.not. mins > 0 .or. .not. ieee_is_finite(ev)) status = IDEAL_FB_FOLDED
        ! Gamma_c terms.
        ec = 0
        do m = 1, fb%nmode
            cm = ideal_fb_mode_weight(fb, m)/fb%mu0
            k = harmonic(m)
            u = fb%phi(m, fb%trace) - fb%drive(:, m)
            ec = ec + cm*(0.5_dp*dot_product(u, matmul(fb%s_trace(:, :, k), u)) &
                - dot_product(fb%load(:, m), fb%phi(m, fb%trace)))
            gf(m, fb%trace) = gf(m, fb%trace) + cm*(matmul(fb%s_trace(:, :, k), u) &
                - fb%load(:, m))
        end do
        call constraint_residual(fb, res)
        ! C^T lambda on xs and on the plasma Gamma_p dofs.
        cx = 0
        do s = 1, fb%ns
            select case (fb%role(s))
            case (1)
                cx(:, :, s) = cx(:, :, s) + fb%lam(:, :, s)
                a = fb%plasma_node(s)
                d = fb%plasma%direction(:, a)
                do m = 1, fb%nmode
                    p = ideal_map_rdof(fb%plasma, 1, m, a)
                    grad(p) = grad(p) - dot_product(d, fb%lam(1:2, m, s))
                    p = ideal_map_rdof(fb%plasma, 2, m, a)
                    grad(p) = grad(p) - fb%lam(3, m, s)
                end do
            case (3)
                cx(:, :, s) = cx(:, :, s) + fb%lam(:, :, s)
            case default
                do p = fb%lap_ptr(s), fb%lap_ptr(s + 1) - 1
                    cx(:, :, fb%lap_col(p)) = cx(:, :, fb%lap_col(p)) &
                        + fb%lap_val(p)*fb%lam(:, :, s)
                end do
            end select
        end do
        do s = 1, fb%ns
            do m = 1, fb%nmode
                do c = 1, 3
                    grad(ix(fb, c, m, s)) = -gx(c, m, s) + cx(c, m, s)
                    grad(il(fb, c, m, s)) = res(c, m, s)
                end do
                grad(iphi(fb, m, s)) = -gf(m, s)
            end do
        end do
        where (fb%plasma%fixed) grad(1:fb%np) = 0
        if (present(lagrangian)) lagrangian = wp - ev - ec + sum(fb%lam*res)
        if (present(min_jacobian)) min_jacobian = min(minp, mins)
    end subroutine

    subroutine ideal_fb_hessian(fb, rowptr, colidx, val)
        ! Assembled KKT matrix of L (symmetric, indefinite), CSR with sorted
        ! columns; fixed plasma dofs get identity rows and columns.
        type(ideal_fb_t), intent(in) :: fb
        integer, allocatable, intent(out) :: rowptr(:), colidx(:)
        real(dp), allocatable, intent(out) :: val(:)
        type(csr_matrix_t) :: hp, hs
        real(dp), allocatable :: blocks(:, :, :), gx(:, :, :), gf(:, :)
        integer, allocatable :: sconn(:, :), ci(:), cj(:)
        real(dp), allocatable :: cv(:)
        real(dp) :: ev, mins, cm, d(2)
        integer :: e, a, m, c, l, nl, cnt, i, k, p, s, t1, t2, nsb, cap
        nl = 4*fb%nmode*fb%nloc
        call ideal_map_hessian(fb%plasma, fb%plasma%x, hp)
        allocate (blocks(nl, nl, fb%ncell), gx(3, fb%nmode, fb%ns), gf(fb%nmode, fb%ns))
        call shell_energy(fb, fb%xs, fb%phi, ev, gx, gf, mins, blocks)
        ! Shell block numbering: xs as 1..nx, Phi~ as nx + 1..nx + nf.
        allocate (sconn(nl, fb%ncell))
        do e = 1, fb%ncell
            do a = 1, fb%nloc
                do m = 1, fb%nmode
                    l = 4*(m - 1) + 4*fb%nmode*(a - 1)
                    do c = 1, 3
                        sconn(l + c, e) = ix(fb, c, m, fb%conn(a, e)) - fb%np
                    end do
                    sconn(l + 4, e) = fb%nx + m + fb%nmode*(fb%conn(a, e) - 1)
                end do
            end do
        end do
        nsb = fb%nx + fb%nf
        call csr_build_pattern(nsb, sconn, .false., hs)
        call csr_assemble_blocks(hs, sconn, blocks)
        deallocate (blocks)
        cap = size(hp%val) + size(hs%val) + fb%nmode*fb%nt**2 &
            + 2*3*fb%nmode*(size(fb%lap_val) + 2*fb%ns) + fb%np
        allocate (ci(cap), cj(cap), cv(cap))
        cnt = 0
        do i = 1, fb%np
            if (fb%plasma%fixed(i)) then
                call push(i, i, 1.0_dp)
                cycle
            end if
            do p = hp%rowptr(i), hp%rowptr(i + 1) - 1
                if (fb%plasma%fixed(hp%colidx(p))) cycle
                call push(i, hp%colidx(p), hp%val(p))
            end do
        end do
        do i = 1, nsb
            do p = hs%rowptr(i), hs%rowptr(i + 1) - 1
                call push(shell_global(i), shell_global(hs%colidx(p)), -hs%val(p))
            end do
        end do
        do m = 1, fb%nmode
            cm = ideal_fb_mode_weight(fb, m)/fb%mu0
            k = harmonic(m)
            do t2 = 1, fb%nt
                do t1 = 1, fb%nt
                    call push(iphi(fb, m, fb%trace(t1)), iphi(fb, m, fb%trace(t2)), &
                        -cm*fb%s_trace(t1, t2, k))
                end do
            end do
        end do
        do s = 1, fb%ns
            do m = 1, fb%nmode
                do c = 1, 3
                    select case (fb%role(s))
                    case (1)
                        call pair(il(fb, c, m, s), ix(fb, c, m, s), 1.0_dp)
                        a = fb%plasma_node(s)
                        d = fb%plasma%direction(:, a)
                        if (c < 3) then
                            call pair(il(fb, c, m, s), ideal_map_rdof(fb%plasma, 1, m, &
                                a), -d(min(c, 2)))
                        else
                            call pair(il(fb, c, m, s), ideal_map_rdof(fb%plasma, 2, m, &
                                a), -1.0_dp)
                        end if
                    case (3)
                        call pair(il(fb, c, m, s), ix(fb, c, m, s), 1.0_dp)
                    case default
                        do p = fb%lap_ptr(s), fb%lap_ptr(s + 1) - 1
                            call pair(il(fb, c, m, s), ix(fb, c, m, fb%lap_col(p)), &
                                fb%lap_val(p))
                        end do
                    end select
                end do
            end do
        end do
        call coo_to_csr(fb%n, ci(:cnt), cj(:cnt), cv(:cnt), rowptr, colidx, val)
    contains
        integer function shell_global(q)
            integer, intent(in) :: q
            shell_global = fb%np + q
            if (q > fb%nx) shell_global = fb%np + fb%nx + q
        end function

        subroutine push(r, col, v)
            integer, intent(in) :: r, col
            real(dp), intent(in) :: v
            cnt = cnt + 1
            ci(cnt) = r; cj(cnt) = col; cv(cnt) = v
        end subroutine

        subroutine pair(r, col, v)
            ! Constraint entry and its transpose; plasma fixed columns skipped.
            integer, intent(in) :: r, col
            real(dp), intent(in) :: v
            if (col <= fb%np) then
                if (fb%plasma%fixed(col)) return
            end if
            call push(r, col, v)
            call push(col, r, v)
        end subroutine
    end subroutine

    subroutine coo_to_csr(n, ci, cj, cv, rowptr, colidx, val)
        integer, intent(in) :: n, ci(:), cj(:)
        real(dp), intent(in) :: cv(:)
        integer, allocatable, intent(out) :: rowptr(:), colidx(:)
        real(dp), allocatable, intent(out) :: val(:)
        integer(i64), allocatable :: keys(:)
        integer, allocatable :: order(:), rows(:)
        integer(i64) :: previous
        integer :: i, nnz
        allocate (keys(size(ci)), order(size(ci)))
        keys = int(ci - 1, i64)*n + cj - 1
        call sort_key_order(keys, order)
        allocate (colidx(size(ci)), val(size(ci)), rows(size(ci)))
        nnz = 0; previous = -1
        do i = 1, size(ci)
            if (keys(order(i)) /= previous) then
                nnz = nnz + 1
                previous = keys(order(i))
                rows(nnz) = ci(order(i)); colidx(nnz) = cj(order(i)); val(nnz) = 0
            end if
            val(nnz) = val(nnz) + cv(order(i))
        end do
        colidx = colidx(:nnz); val = val(:nnz)
        allocate (rowptr(n + 1))
        rowptr = 0
        do i = 1, nnz
            rowptr(rows(i) + 1) = rowptr(rows(i) + 1) + 1
        end do
        rowptr(1) = 1
        do i = 1, n
            rowptr(i + 1) = rowptr(i + 1) + rowptr(i)
        end do
    end subroutine

    subroutine ideal_fb_step(fb, du, alpha)
        ! U <- U + alpha du in the global layout.
        type(ideal_fb_t), intent(inout) :: fb
        real(dp), intent(in) :: du(:), alpha
        real(dp), allocatable :: full(:), dr(:)
        integer :: s, m, c
        allocate (full(fb%plasma%ndof))
        dr = alpha*du(1:fb%np)
        where (fb%plasma%fixed) dr = 0
        call ideal_map_expand(fb%plasma, dr, full)
        fb%plasma%x = fb%plasma%x + reshape(full, shape(fb%plasma%x))
        do s = 1, fb%ns
            do m = 1, fb%nmode
                do c = 1, 3
                    fb%xs(c, m, s) = fb%xs(c, m, s) + alpha*du(ix(fb, c, m, s))
                    fb%lam(c, m, s) = fb%lam(c, m, s) + alpha*du(il(fb, c, m, s))
                end do
                fb%phi(m, s) = fb%phi(m, s) + alpha*du(iphi(fb, m, s))
            end do
        end do
    end subroutine

    subroutine ideal_fb_solve(fb, tolerance, max_iterations, report, status, verbose, &
                              step_tolerance)
        ! Newton on the KKT system of L (symmetric, indefinite: negative
        ! vacuum block, multipliers, and with fixed coil currents a negative
        ! vertical direction of the reduced plasma Hessian). Globalisation:
        !   - merit = sum g_i^2/|H_ii| over primal rows, g_i^2 |H_xx| over
        !     constraint rows; a step is accepted when the merit decreases and
        !     every Jacobian stays positive (never on the energy, which is a
        !     saddle);
        !   - Levenberg damping mu |H_ii| on the plasma map rows only. It
        !     regularises the continuum gauge directions of the map (rigid
        !     rotation of surfaces, omega_0 = f(chi)) that are near-null at
        !     n = 0; it is not applied to the vacuum or multiplier blocks and is
        !     reduced by 100 after every full step, so Newton is recovered;
        !   - trust radii: displacement <= 2 % of the Gamma_c extent, omega
        !     <= 0.5 rad; a capped direction raises mu.
        ! Converged when sqrt(merit/merit0) <= tolerance or the last
        ! displacement increment is below step_tolerance (m).
        type(ideal_fb_t), intent(inout) :: fb
        real(dp), intent(in) :: tolerance
        integer, intent(in) :: max_iterations
        type(ideal_fb_report_t), intent(out) :: report
        integer, intent(out) :: status
        logical, intent(in), optional :: verbose
        real(dp), intent(in), optional :: step_tolerance
        integer, allocatable :: rowptr(:), colidx(:), dpos(:)
        real(dp), allocatable :: val(:), dval(:), grad(:), du(:), diag(:), wts(:)
        real(dp), allocatable :: rhs(:), x0(:, :, :), xs0(:, :, :), phi0(:, :)
        real(dp), allocatable :: lam0(:, :, :)
        type(sparse_solver_t) :: solver
        type(fortsparse_status_t) :: sst
        real(dp) :: merit, tmerit, alpha, cap_x, cap_w, disp, rot, minjac, t0
        real(dp) :: step_tol, mu, scale_p
        integer :: it, i, p, st, tries
        logical :: talk, accepted
        talk = .false.
        if (present(verbose)) talk = verbose
        step_tol = 1.0e-10_dp
        if (present(step_tolerance)) step_tol = step_tolerance
        cap_x = 0.02_dp*maxval(maxval(fb%reference, dim=2) - minval(fb%reference, dim=2))
        cap_w = 0.5_dp
        mu = 1.0e-8_dp
        allocate (grad(fb%n), du(fb%n), diag(fb%n), wts(fb%n), rhs(fb%n), dpos(fb%np))
        call ideal_fb_gradient(fb, grad, st)
        if (st /= IDEAL_FB_OK) then
            status = IDEAL_FB_FOLDED; return
        end if
        status = IDEAL_FB_NOT_CONVERGED
        do it = 1, max_iterations
            t0 = wall()
            call ideal_fb_hessian(fb, rowptr, colidx, val)
            report%assemble_s = report%assemble_s + wall() - t0
            diag = 0
            do i = 1, fb%n
                do p = rowptr(i), rowptr(i + 1) - 1
                    if (colidx(p) == i) then
                        diag(i) = abs(val(p))
                        if (i <= fb%np) dpos(i) = p
                    end if
                end do
            end do
            scale_p = maxval(diag(1:fb%np), mask=.not. fb%plasma%fixed)
            call merit_weights(fb, diag, wts)
            merit = sum(wts*grad**2)
            if (it == 1) report%merit0 = merit
            report%merit = merit
            report%residual = sqrt(merit/fb%energy_scale)
            if (talk) print '(a,i3,a,es10.3,a,es10.3,a,es10.3,a,es9.2)', 'fb newton', &
                it - 1, ' merit/m0 ', sqrt(merit/report%merit0), ' res/sigmaV ', &
                report%residual, ' step ', report%step, ' mu ', mu
            if (.not. sqrt(merit) > tolerance*sqrt(report%merit0)) then
                status = IDEAL_FB_OK; exit
            end if
            if (report%step < step_tol .and. mu < 1.0e-10_dp) then
                status = IDEAL_FB_OK; exit
            end if
            accepted = .false.
            x0 = fb%plasma%x; xs0 = fb%xs; phi0 = fb%phi; lam0 = fb%lam
            do tries = 1, 40
                dval = val
                do i = 1, fb%np
                    if (fb%plasma%fixed(i)) cycle
                    dval(dpos(i)) = dval(dpos(i)) + mu*(diag(i) + 1.0e-12_dp*scale_p)
                end do
                t0 = wall()
                call configure(solver)
                call sparse_factor(solver, fb%n, fb%n, size(dval), rowptr, colidx, &
                    dval, sst)
                if (.not. status_ok(sst)) then
                    call sparse_destroy(solver)
                    mu = 10*mu; cycle
                end if
                rhs = -grad
                call sparse_solve(solver, rhs, du, sst)
                call sparse_destroy(solver)
                report%factor_s = report%factor_s + wall() - t0
                if (.not. status_ok(sst)) then
                    status = IDEAL_FB_SOLVE_FAILED; return
                end if
                call map_steps(fb, du, disp, rot)
                alpha = min(1.0_dp, cap_x/max(disp, tiny(1.0_dp)), &
                    cap_w/max(rot, tiny(1.0_dp)))
                if (talk) print '(a,i3,4es10.2)', '   try mu disp rot alpha', tries, &
                    mu, disp, rot, alpha
                if (alpha < 1 .and. mu < 1.0e8_dp) then
                    mu = 10*mu; cycle
                end if
                call ideal_fb_step(fb, du, alpha)
                call ideal_fb_gradient(fb, rhs, st, min_jacobian=minjac)
                if (st == IDEAL_FB_OK) then
                    tmerit = sum(wts*rhs**2)
                    if (tmerit < (1 - 1.0e-4_dp*alpha)*merit) then
                        accepted = .true.; exit
                    end if
                end if
                fb%plasma%x = x0; fb%xs = xs0; fb%phi = phi0; fb%lam = lam0
                mu = 10*mu
            end do
            if (.not. accepted) exit
            grad = rhs
            report%iterations = it
            report%step = alpha*disp
            report%min_jacobian = minjac
            report%damping = mu
            if (.not. alpha < 1) mu = max(mu/100, 1.0e-14_dp)
        end do
        call constraint_norm(fb, report%constraint)
    contains
        real(dp) function wall()
            integer(i64) :: count, rate
            call system_clock(count, rate)
            wall = real(count, dp)/rate
        end function
    end subroutine

    subroutine configure(solver)
        ! Threaded MUMPS LDL^T when available, serial SuperLU LU otherwise.
        type(sparse_solver_t), intent(inout) :: solver
        if (gs_direct_available()) then
            solver%backend_id = FORTSPARSE_BACKEND_MUMPS
            solver%options%symmetry = FORTSPARSE_SYM_INDEFINITE
            solver%options%nthreads = 0
        else
            solver%backend_id = FORTSPARSE_BACKEND_SUPERLU
            solver%options%symmetry = FORTSPARSE_SYM_GENERAL
        end if
    end subroutine

    subroutine map_steps(fb, du, disp, rot)
        ! Largest displacement (m) and omega (rad) increment of a direction.
        type(ideal_fb_t), intent(in) :: fb
        real(dp), intent(in) :: du(:)
        real(dp), intent(out) :: disp, rot
        real(dp), allocatable :: full(:), dr(:)
        integer :: s, m
        allocate (full(fb%plasma%ndof))
        dr = du(1:fb%np)
        where (fb%plasma%fixed) dr = 0
        call ideal_map_expand(fb%plasma, dr, full)
        disp = max(maxval(abs(full(1::3))), maxval(abs(full(2::3))))
        rot = maxval(abs(full(3::3)))
        do s = 1, fb%ns
            do m = 1, fb%nmode
                disp = max(disp, abs(du(ix(fb, 1, m, s))), abs(du(ix(fb, 2, m, s))))
                rot = max(rot, abs(du(ix(fb, 3, m, s))))
            end do
        end do
    end subroutine

    subroutine merit_weights(fb, diag, w)
        type(ideal_fb_t), intent(in) :: fb
        real(dp), intent(in) :: diag(:)
        real(dp), intent(out) :: w(:)
        real(dp) :: floor
        integer :: s, m, c
        floor = 1.0e-12_dp*maxval(diag)
        w = 1/max(diag, floor)
        do s = 1, fb%ns
            do m = 1, fb%nmode
                do c = 1, 3
                    w(il(fb, c, m, s)) = max(diag(ix(fb, c, m, s)), floor)
                end do
            end do
        end do
        where (fb%plasma%fixed) w(1:fb%np) = 0
    end subroutine

    subroutine constraint_norm(fb, value)
        type(ideal_fb_t), intent(in) :: fb
        real(dp), intent(out) :: value
        real(dp), allocatable :: res(:, :, :)
        allocate (res(3, fb%nmode, fb%ns))
        call constraint_residual(fb, res)
        value = maxval(abs(res))
    end subroutine

    subroutine ideal_fb_solve_vacuum(fb, status)
        ! Phi~ at fixed geometry: the vacuum block of the KKT system (linear in
        ! Phi~, negative definite) solved exactly. Used to start Newton from a
        ! consistent vacuum field and by the vacuum oracles.
        type(ideal_fb_t), intent(inout) :: fb
        integer, intent(out) :: status
        real(dp), allocatable :: du(:)
        integer :: off
        off = fb%np + 2*fb%nx
        call block_solve(fb, off, off, fb%nf, du, status)
        if (status == IDEAL_FB_OK) fb%phi = fb%phi + reshape(du, shape(fb%phi))
    end subroutine

    subroutine ideal_fb_solve_multipliers(fb, status)
        ! Multipliers at fixed geometry and potential: the shell-map rows
        ! -dE_v/dx + C^T lambda = 0 (C square and invertible), so that Newton
        ! starts with the shell-map residual removed.
        type(ideal_fb_t), intent(inout) :: fb
        integer, intent(out) :: status
        real(dp), allocatable :: du(:)
        call block_solve(fb, fb%np, fb%np + fb%nx, fb%nx, du, status)
        if (status == IDEAL_FB_OK) fb%lam = fb%lam + reshape(du, shape(fb%lam))
    end subroutine

    subroutine block_solve(fb, row0, col0, nb, du, status)
        ! Solve A(rows row0+1.., cols col0+1..) du = -g(rows) of the KKT system.
        type(ideal_fb_t), intent(in) :: fb
        integer, intent(in) :: row0, col0, nb
        real(dp), allocatable, intent(out) :: du(:)
        integer, intent(out) :: status
        integer, allocatable :: rowptr(:), colidx(:), rp(:), ci(:)
        real(dp), allocatable :: val(:), v(:), grad(:), rhs(:)
        type(sparse_solver_t) :: solver
        type(fortsparse_status_t) :: sst
        integer :: i, p, nnz, st
        status = IDEAL_FB_SOLVE_FAILED
        allocate (grad(fb%n), du(nb))
        call ideal_fb_gradient(fb, grad, st)
        call ideal_fb_hessian(fb, rowptr, colidx, val)
        allocate (rp(nb + 1), ci(size(val)), v(size(val)))
        nnz = 0; rp(1) = 1
        do i = 1, nb
            do p = rowptr(row0 + i), rowptr(row0 + i + 1) - 1
                if (colidx(p) <= col0 .or. colidx(p) > col0 + nb) cycle
                nnz = nnz + 1
                ci(nnz) = colidx(p) - col0; v(nnz) = val(p)
            end do
            rp(i + 1) = nnz + 1
        end do
        if (row0 == col0) then
            call configure(solver)
        else
            solver%backend_id = FORTSPARSE_BACKEND_SUPERLU
            if (gs_direct_available()) solver%backend_id = FORTSPARSE_BACKEND_MUMPS
            solver%options%symmetry = FORTSPARSE_SYM_GENERAL
        end if
        call sparse_factor(solver, nb, nb, nnz, rp, ci(:nnz), v(:nnz), sst)
        if (.not. status_ok(sst)) return
        rhs = -grad(row0 + 1:row0 + nb)
        call sparse_solve(solver, rhs, du, sst)
        call sparse_destroy(solver)
        if (.not. status_ok(sst)) return
        status = IDEAL_FB_OK
    end subroutine

    subroutine ideal_fb_shell_field(fb, e, xi, phi_hat, x, b, status)
        ! Physical point x = (R, Z, phi) and cylindrical field (B_R, B_phi,
        ! B_Z) of the vacuum potential at master point xi of shell cell e.
        type(ideal_fb_t), intent(in) :: fb
        integer, intent(in) :: e
        real(dp), intent(in) :: xi(2), phi_hat
        real(dp), intent(out) :: x(3), b(3)
        integer, intent(out) :: status
        real(dp) :: shape(10), grad(2, 10), q(3), qd(3, 3), gt(3), t, dt, nn
        real(dp) :: dg(2, 2), y(2), gg(3), cart(3), d1(3, 3), d2(3, 3, 3), dx(3, 3)
        real(dp) :: d3(3, 3, 3, 3), inv(3, 3), bc(3), det
        integer :: a, m, node, k
        status = 1
        call triangle_basis_jet(xi, shape(:fb%nloc), grad(:, :fb%nloc))
        q = 0; qd = 0; gt = 0
        do a = 1, fb%nloc
            node = fb%conn(a, e)
            do m = 1, fb%nmode
                k = harmonic(m); nn = k*fb%nfp
                if (m == 1) then
                    t = 1; dt = 0
                else if (mod(m, 2) == 0) then
                    t = cos(nn*phi_hat); dt = -nn*sin(nn*phi_hat)
                else
                    t = sin(nn*phi_hat); dt = nn*cos(nn*phi_hat)
                end if
                q = q + fb%xs(:, m, node)*shape(a)*t
                qd(:, 1) = qd(:, 1) + fb%xs(:, m, node)*grad(1, a)*t
                qd(:, 2) = qd(:, 2) + fb%xs(:, m, node)*grad(2, a)*t
                qd(:, 3) = qd(:, 3) + fb%xs(:, m, node)*shape(a)*dt
                gt(1) = gt(1) + fb%phi(m, node)*grad(1, a)*t
                gt(2) = gt(2) + fb%phi(m, node)*grad(2, a)*t
                gt(3) = gt(3) + fb%phi(m, node)*shape(a)*dt
            end do
        end do
        q(3) = q(3) + phi_hat; qd(3, 3) = qd(3, 3) + 1
        y = matmul(fb%reference(:, fb%conn(:, e)), shape(:fb%nloc))
        dg = matmul(fb%reference(:, fb%conn(:, e)), transpose(grad(:, :fb%nloc)))
        gg = gt + fb%f_vac*qd(3, :) + fb%circulation*[matmul(transpose(dg), &
            theta_gradient(fb, y)), 0.0_dp]
        call fb%plasma%chart%jet(q, cart, d1, d2, d3)
        dx = matmul(d1, qd)
        call inverse3(dx, inv, det)
        if (.not. abs(det) > 0) return
        ! Cartesian B = DX^{-T} g, then cylindrical components at q.
        bc = matmul(transpose(inv), gg)
        b(1) = bc(1)*cos(q(3)) + bc(2)*sin(q(3))
        b(2) = -bc(1)*sin(q(3)) + bc(2)*cos(q(3))
        b(3) = bc(3)
        x = q([1, 2, 3])
        status = 0
    end subroutine

    pure subroutine inverse3(a, inv, det)
        real(dp), intent(in) :: a(3, 3)
        real(dp), intent(out) :: inv(3, 3), det
        inv(1, 1) = a(2, 2)*a(3, 3) - a(2, 3)*a(3, 2)
        inv(1, 2) = a(1, 3)*a(3, 2) - a(1, 2)*a(3, 3)
        inv(1, 3) = a(1, 2)*a(2, 3) - a(1, 3)*a(2, 2)
        inv(2, 1) = a(2, 3)*a(3, 1) - a(2, 1)*a(3, 3)
        inv(2, 2) = a(1, 1)*a(3, 3) - a(1, 3)*a(3, 1)
        inv(2, 3) = a(1, 3)*a(2, 1) - a(1, 1)*a(2, 3)
        inv(3, 1) = a(2, 1)*a(3, 2) - a(2, 2)*a(3, 1)
        inv(3, 2) = a(1, 2)*a(3, 1) - a(1, 1)*a(3, 2)
        inv(3, 3) = a(1, 1)*a(2, 2) - a(1, 2)*a(2, 1)
        det = a(1, 1)*inv(1, 1) + a(1, 2)*inv(2, 1) + a(1, 3)*inv(3, 1)
        if (abs(det) > 0) inv = inv/det
    end subroutine
end module kin6d_ideal_free_boundary
