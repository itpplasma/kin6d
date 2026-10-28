module kin6d_ideal_free_boundary
    ! 3D free boundary of the nested map lane (DESIGN 7.3, core note D6, N4).
    ! One reference mesh of the region inside a fixed axisymmetric coupling
    ! curve Gamma_c, with the plasma boundary Gamma_p as an internal chain:
    !   plasma cells: the frozen-in energy of kin6d_ideal_map (unchanged), in
    !     the position gauge: every node (Gamma_p included) carries its full
    !     poloidal position (q1, q2) per Fourier mode and the toroidal
    !     displacement omega is zero. Relabelling along B needs a toroidal
    !     component, so this gauge has no relabelling null; it represents
    !     rigid motions exactly. (N1's transversal + omega gauge cannot slide
    !     material along surfaces without omega_0, whose surface-function part
    !     is a continuum null; pinning omega_0 there over-constrains the map
    !     and stalls an offset start 7 mm short of the equilibrium; leaving it
    !     free makes Newton crawl along the null. Measured on the gs_free
    !     fixture.)
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
    !     traces on the fixed Gamma_c. The per-n map needs an axisymmetric
    !     Gamma_c, i.e. exactly the cylindrical chart. Any other chart (curve
    !     chart, nested-map chart) gives the 3D Gamma_c = H(gamma x [0, 2 pi))
    !     and takes the full-surface exterior (R2) through
    !     ideal_fb_attach_exterior_3d / ideal_fb_use_bem3d and its traces
    !     through ideal_fb_drive_traces_3d; the residual is unchanged.
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
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64, i64 => int64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, &
        ieee_quiet_nan
    use kin6d_grad_shafranov, only: gs_state_t
    use kin6d_equilibrium_chart, only: ideal_chart_t, cylindrical_chart_t
    use kin6d_ideal_map, only: ideal_map_t, ideal_frozen_t, &
        ideal_map_init, ideal_map_freeze, ideal_map_energy, ideal_map_dof, &
        IDEAL_MAP_OK
    use kin6d_ideal_map_leaves, only: ideal_energy_gradient, ideal_energy_hessian
    use kin6d_vacuum_shell_leaves, only: vacuum_energy_density, &
        vacuum_energy_gradient, vacuum_energy_hessian
    use kin6d_triangle_basis, only: triangle_basis_jet
    use kin6d_sort_keys, only: sort_key_order
    use kin6d_gs_direct, only: gs_direct_available, gs_direct_mpi_ready
    use fortfem_element_to_csr, only: csr_matrix_t, csr_build_pattern, &
        csr_assemble_blocks
    use fortfem_axisymmetric_gs_boundary_operators_2d, only: gs_bem_curve_t, &
        gs_bem_operators_t, gs_bem_curve_from_nodes
    use fortfem_toroidal_mode_boundary_operators_2d, only: tm_bem_assemble, &
        tm_bem_steklov_symmetric
    use fortnum_quadrature, only: gauss_legendre
    use fortnum_krylov, only: real_fgmres_operator, KRYLOV_OK
    use fortnum_eigen, only: lanczos, eigen_options_t, eigen_stats_t, &
        EIGEN_OK, EIGEN_NEAREST
    use fortsparse, only: sparse_solver_t, sparse_factor, sparse_solve, &
        sparse_refactor, sparse_destroy, fortsparse_status_t, status_ok, &
        FORTSPARSE_BACKEND_NATIVE, FORTSPARSE_SYM_GENERAL, &
        csc_normalize_symmetric_raw, &
        FORTSPARSE_ORDER_AMD, FORTSPARSE_LU_SUPERNODAL, &
        FORTSPARSE_SYM_INDEFINITE, sparse_inertia, FORTSPARSE_INERTIA_RELIABLE
    use kin6d_equilibrium_boundary_3d, only: gamma_c_3d_t, gamma_c_3d_energy, &
        assignment(=)
    implicit none
    private
    public :: ideal_fb_t, ideal_fb_report_t, ideal_fb_init, ideal_fb_set_drive
    public :: ideal_fb_gradient, ideal_fb_hessian, ideal_fb_step, ideal_fb_solve
    public :: ideal_fb_shell_field, ideal_fb_sampler, ideal_fb_mode_weight
    public :: ideal_fb_solve_vacuum, ideal_fb_solve_multipliers
    public :: ideal_fb_solve_shell_map, ideal_fb_rigid_stiffness
    public :: ideal_fb_mode_stability, ideal_fb_mode_energy, ideal_fb_stability_t
    public :: ideal_fb_mode_inertia, ideal_fb_mode_lanczos
    public :: ideal_fb_newton_step
    public :: ideal_fb_slide_system, ideal_fb_hypotheses_t, ideal_fb_hypotheses
    public :: ideal_fb_hypotheses_print
    public :: ideal_fb_attach_exterior_3d, ideal_fb_exterior_ready
    public :: ideal_fb_field_basis_t, ideal_fb_drive_traces, ideal_fb_full_turn_traces
    public :: ideal_fb_set_traces, ideal_fb_traces_jvp, ideal_fb_traces_vjp
    public :: ideal_fb_drive_traces_3d, ideal_fb_set_drive_3d
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

    ! nfield independent external fields (one sampler, or the unit-current
    ! fields b_k of circuits) at the meridional points rz(:, q), all real
    ! harmonics k = 0..ntor at once:
    ! B_f = bcos(:, k, f, q) cos(n phi) + bsin(:, k, f, q) sin(n phi), n = k nfp.
    ! Cartesian contract (sample_points): B_f at Cartesian points x(:, q) in
    ! Cartesian components b(:, f, q), for coupling surfaces that are not
    ! axisymmetric. The default sums the harmonics k = 0..ntor of sample at
    ! (R, Z, phi) of each point (exact for fields without other harmonics);
    ! a basis with a direct evaluation (coils) overrides it.
    type, abstract :: ideal_fb_field_basis_t
    contains
        procedure(basis_nfield), deferred :: nfield
        procedure(basis_sample), deferred :: sample
        procedure :: sample_points => basis_sample_points
    end type
    abstract interface
        integer function basis_nfield(self)
            import :: ideal_fb_field_basis_t
            class(ideal_fb_field_basis_t), intent(in) :: self
        end function
        subroutine basis_sample(self, rz, nfp, ntor, bcos, bsin, status)
            import :: dp, ideal_fb_field_basis_t
            class(ideal_fb_field_basis_t), intent(in) :: self
            real(dp), intent(in) :: rz(:, :)
            integer, intent(in) :: nfp, ntor
            real(dp), intent(out) :: bcos(:, 0:, :, :), bsin(:, 0:, :, :)
            integer, intent(out) :: status
        end subroutine
    end interface

    ! One ideal_fb_sampler as a single-field basis (ideal_fb_set_drive).
    type, extends(ideal_fb_field_basis_t) :: sampler_basis_t
        procedure(ideal_fb_sampler), pointer, nopass :: f => null()
    contains
        procedure :: nfield => sampler_nfield
        procedure :: sample => sampler_sample
    end type

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
        ! Optional full-surface R2 operator. Flattened traces are node-fast,
        ! then real Fourier mode. The surface integrals cover a full turn.
        type(gamma_c_3d_t), allocatable :: exterior_3d
        real(dp), allocatable :: drive_3d(:), load_3d(:)
        ! State: shell map coefficients (3, nmode, ns), potential (nmode, ns),
        ! multipliers (3, nmode, ns).
        real(dp), allocatable :: xs(:, :, :), phi(:, :), lam(:, :, :)
        ! Unknown layout: plasma reduced, shell map, multipliers, potential.
        integer :: np = 0, nx = 0, nf = 0, n = 0
        ! Energy scale sigma*V (J per field period) for the reported residual.
        real(dp) :: energy_scale = 1
        ! Tangential damping (Newton matrix only) on the plasma nodes, relative
        ! to their diagonal: bounds steps along the n = 0 continuum null of the
        ! position gauge (surface rotation v ~ grad f(chi) x grad phi / B_phi,
        ! which keeps the surfaces, the poloidal flux and the toroidal-flux
        ! 2-form). Must stay below the vertical eigenvalue (~1e-7 in this
        ! metric) or the damped step turns against the vertical mode.
        real(dp) :: tangential_damping = 0
        ! After every accepted Newton step, try the rigid n = 0 Z correction
        ! dZ = -F_Z/K_Z (K_Z from ideal_fb_rigid_stiffness), kept only when the
        ! merit decreases.
        logical :: rigid_vertical = .false.
        ! Evaluate ideal_fb_hypotheses at every Newton state and print it
        ! (verbose solves only; a diagnostic, never a stopping criterion).
        logical :: monitor_hypotheses = .false.
        ! Relative asymmetry ||A - A^T||_F/||A||_F of the assembled Galerkin
        ! hypersingular W (1) and single-layer V (2) operators per harmonic k,
        ! before FortFEM symmetrises them into S_k (FB-A2).
        real(dp), allocatable :: bem_asymmetry(:, :)
        ! Slide-moment rows of the n = 0 relabelling gauge (2, plasma node,
        ! layer), unit Euclidean norm; see slide_directions.
        real(dp), allocatable :: slide(:, :, :)
        ! Newton-homotopy offset: when allocated, the gradient returned is
        ! g(U) - shift (shift = (1 - t) g(U0) solves g = (1 - t) g(U0)).
        real(dp), allocatable :: shift(:)
    end type

    type :: ideal_fb_report_t
        integer :: iterations = 0
        real(dp) :: merit0 = 0, merit = 0, residual = huge(1.0_dp)
        real(dp) :: step = huge(1.0_dp), boundary_step = huge(1.0_dp)
        real(dp) :: constraint = 0, min_jacobian = 0, damping = 0
        ! Stopped at the roundoff floor of the merit.
        logical :: floor = .false.
        real(dp) :: assemble_s = 0, factor_s = 0, gradient_s = 0, line_s = 0
        real(dp) :: assembly_roundoff = 0
        integer :: factorisations = 0, reused = 0
    end type

    type :: ideal_fb_stability_t
        ! Raw second-variation ledger, per field period. These are half of
        ! u^T H u in the fixed-current co-energy formulation, not a kinetic
        ! norm or DCON's boundary-normalized response eigenvalues.
        real(dp) :: plasma = 0, vacuum = 0, exterior = 0, constraint = 0
        real(dp) :: total = 0, metric_norm = 0, relative_residual = huge(1.0_dp)
        real(dp) :: residual_radius = huge(1.0_dp), backward_error = huge(1.0_dp)
    end type

    type :: ideal_fb_hypotheses_t
        ! Per-state checks of the open discrete hypotheses of
        ! docs/free-boundary-theory.md (FB-A2, FB-C5, FB-D2, FB-E1). These are
        ! diagnostics, never stopping criteria. status = IDEAL_FB_OK when every
        ! entry was evaluated; the dense parts are limited to fixture sizes.
        integer :: status = IDEAL_FB_NOT_CONVERGED
        ! (V), FB-A2: largest relative asymmetry of the Galerkin W and V
        ! before FortFEM symmetrises them (S_k is symmetric by construction),
        ! smallest eigenvalue of S_k over k; extreme eigenvalues and relative
        ! asymmetry of -mu0 H_PhiPhi, the assembled A_h + S_h block (times
        ! the mode weights c_m).
        real(dp) :: w_asymmetry = 0, v_asymmetry = 0, s_min = 0
        real(dp) :: vacuum_min = 0, vacuum_max = 0, vacuum_asymmetry = 0
        ! (R): R = [grad C; G] on the map unknowns (G the unit slide rows):
        ! rows, extreme singular values (inverse/power iteration) and full row
        ! rank, decided as sigma_min > rank_tolerance sigma_max.
        integer :: constraint_rows = 0
        real(dp) :: constraint_sigma_min = 0, constraint_sigma_max = 0
        real(dp) :: rank_tolerance = 1.0e-10_dp
        logical :: full_rank = .false.
        ! S_red = Z^T S_map Z, Z = [I; E] Z_G: E the shell-map extension of
        ! the plasma positions, Z_G an orthonormal basis of ker G in the
        ! Euclidean metric of the plasma position coefficients (m), so the
        ! eigenvalues are J/m^2 per field period. Inertia (n+, n-, n0) with
        ! n0 counting |lambda| <= zero_tolerance = 1e-10 max |lambda|;
        ! vertical_overlap = |<z_min, unit rigid Z shift>| of the lowest
        ! eigenvector.
        integer :: reduced_dim = 0, n_positive = 0, n_negative = 0, n_zero = 0
        real(dp) :: zero_tolerance = 0, reduced_min = 0, reduced_min_abs = 0
        real(dp) :: reduced_max = 0, vertical_overlap = 0
        ! Inertia of the bordered KKT matrix implied by FB-D2 under (V) and (R):
        ! (n+ + r, n- + dim Phi~ + r, n0), r = constraint_rows.
        integer :: kkt_inertia(3) = 0
        ! FB-C5: largest slide-gauge multiplier (J/m per field period, unit
        ! slide rows) and the residual it removes, in sigma V units.
        real(dp) :: slide_multiplier = 0, slide_residual = 0
    end type
contains
    pure real(dp) function ideal_fb_mode_weight(fb, m) result(c)
        ! Field-period integral of the squared real Fourier mode m.
        type(ideal_fb_t), intent(in) :: fb
        integer, intent(in) :: m
        c = merge(2*pi, pi, m == 1)/fb%nfp
    end function

    pure logical function ideal_fb_exterior_ready(fb) result(ready)
        ! An exterior is attached: the per-n map (cylindrical chart) or R2.
        type(ideal_fb_t), intent(in) :: fb
        ready = allocated(fb%exterior_3d) .or. allocated(fb%s_trace)
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

    pure integer function pidx(fb, c, m, a)
        ! Plasma unknown: position component c (q1, q2) of mode m at node a.
        type(ideal_fb_t), intent(in) :: fb
        integer, intent(in) :: c, m, a
        pidx = c + 2*(m - 1) + 2*fb%nmode*(a - 1)
    end function

    pure integer function iphi(fb, m, s)
        type(ideal_fb_t), intent(in) :: fb
        integer, intent(in) :: m, s
        iphi = fb%np + 2*fb%nx + m + fb%nmode*(s - 1)
    end function

    subroutine ideal_fb_init(fb, mesh, chain, chart, field, nfp, ntor, f_vac, &
                             circulation, status, mu0, nphi, panels_per_edge, &
                             quadrature_order)
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
        type(gs_state_t) :: sub
        logical, allocatable :: inside(:), in_plasma(:), in_shell(:)
        integer :: e, a, nodes, cells, nplasma, k, st, ppe
        real(dp) :: centre(2)
        logical :: axisymmetric
        status = IDEAL_FB_BAD_INPUT
        ! Only the cylindrical chart itself (not an extension that overrides
        ! its jet) has the axisymmetric Gamma_c of the per-n exterior. Every
        ! other chart leaves the exterior to the R2 attachment.
        axisymmetric = .false.
        select type (chart)
        type is (cylindrical_chart_t)
            axisymmetric = .true.
        end select
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
        call ideal_map_freeze(fb%plasma, field, st)
        if (st /= IDEAL_MAP_OK) return
        call slide_directions(fb)
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
        if (axisymmetric) then
            call gamma_c_operators(fb, ppe, st)
            if (st /= 0) return
        end if
        fb%np = 2*fb%nmode*fb%plasma%nodes; fb%nx = 3*fb%nmode*fb%ns; fb%nf = fb%nmode*fb%ns
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
        allocate (fb%s_trace(fb%nt, fb%nt, 0:fb%ntor), fb%bem_asymmetry(2, 0:fb%ntor))
        do k = 0, fb%ntor
            call tm_bem_assemble(curve, k*fb%nfp, 2, ops, st)
            if (st /= 0) return
            fb%bem_asymmetry(:, k) = [norm2(ops%w - transpose(ops%w))/norm2(ops%w), &
                norm2(ops%v - transpose(ops%v))/norm2(ops%v)]
            call tm_bem_steklov_symmetric(ops, s_bem, st)
            if (st /= 0) return
            fb%s_trace(:, :, k) = matmul(transpose(restriction), &
                matmul(s_bem, restriction))
            fb%s_trace(:, :, k) = 0.5_dp*(fb%s_trace(:, :, k) &
                + transpose(fb%s_trace(:, :, k)))
        end do
        status = 0
    end subroutine

    subroutine ideal_fb_attach_exterior_3d(fb, boundary, drive, load, status)
        type(ideal_fb_t), intent(inout) :: fb
        type(gamma_c_3d_t), intent(in) :: boundary
        real(dp), intent(in) :: drive(:), load(:)
        integer, intent(out) :: status
        integer :: n
        status = IDEAL_FB_BAD_INPUT
        n = fb%nt*fb%nmode
        if (boundary%genus /= 1) return
        if (.not. allocated(boundary%s)) return
        if (any(shape(boundary%s) /= [n,n])) return
        if (size(drive) /= n .or. size(load) /= n) return
        if (.not. all(ieee_is_finite(boundary%s))) return
        if (.not. all(ieee_is_finite(drive))) return
        if (.not. all(ieee_is_finite(load))) return
        if (allocated(fb%exterior_3d)) deallocate(fb%exterior_3d)
        allocate(fb%exterior_3d)
        fb%exterior_3d = boundary
        fb%drive_3d = drive
        fb%load_3d = load
        status = IDEAL_FB_OK
    end subroutine

    subroutine ideal_fb_full_turn_traces(fb, drive, load, drive3, load3)
        ! Per-period Fourier traces (nt, nmode) -> the flattened full-turn
        ! pairing of the R2 exterior (node fast, then mode; load times
        ! nfp c_m), the convention of ideal_fb_use_bem3d.
        type(ideal_fb_t), intent(in) :: fb
        real(dp), intent(in) :: drive(:, :), load(:, :)
        real(dp), intent(out) :: drive3(:), load3(:)
        integer :: m
        drive3 = reshape(drive, [fb%nt*fb%nmode])
        do m = 1, fb%nmode
            load3((m - 1)*fb%nt + 1:m*fb%nt) = fb%nfp*ideal_fb_mode_weight(fb, m) &
                *load(:, m)
        end do
    end subroutine

    subroutine ideal_fb_set_traces(fb, drive, load, f_vac, status)
        ! Install external traces and the toroidal period F_vac (for example
        ! sum_k I_k of circuit traces). With an attached R2 exterior its
        ! full-turn drive and load are replaced consistently.
        type(ideal_fb_t), intent(inout) :: fb
        real(dp), intent(in) :: drive(:, :), load(:, :), f_vac
        integer, intent(out) :: status
        status = IDEAL_FB_BAD_INPUT
        if (any(shape(drive) /= [fb%nt, fb%nmode])) return
        if (any(shape(load) /= [fb%nt, fb%nmode])) return
        if (.not. all(ieee_is_finite(drive))) return
        if (.not. all(ieee_is_finite(load))) return
        if (.not. ieee_is_finite(f_vac)) return
        fb%drive = drive; fb%load = load; fb%f_vac = f_vac
        if (allocated(fb%exterior_3d)) then
            if (.not. allocated(fb%drive_3d)) allocate (fb%drive_3d(fb%nt*fb%nmode))
            if (.not. allocated(fb%load_3d)) allocate (fb%load_3d(fb%nt*fb%nmode))
            call ideal_fb_full_turn_traces(fb, drive, load, fb%drive_3d, fb%load_3d)
        end if
        status = IDEAL_FB_OK
    end subroutine

    subroutine ideal_fb_traces_jvp(fb, ddrive, dload, df, jv, status)
        ! Directional derivative of the residual (ideal_fb_gradient) at the
        ! current state along the external data (drive, load, F_vac). The
        ! residual is affine in drive and load (Gamma_c rows of Phi~ only);
        ! its F_vac part is exact through the generated shell Hessian leaf.
        type(ideal_fb_t), intent(in) :: fb
        real(dp), intent(in) :: ddrive(:, :), dload(:, :), df
        real(dp), intent(out) :: jv(:)
        integer, intent(out) :: status
        real(dp), allocatable :: t(:, :), d3(:), l3(:), gx(:, :, :), gf(:, :)
        real(dp), allocatable :: fgx(:, :, :), fgf(:, :)
        real(dp) :: ev, mins
        integer :: m, s, c
        status = IDEAL_FB_BAD_INPUT
        if (.not. ideal_fb_exterior_ready(fb)) return
        if (size(jv) /= fb%n) return
        if (any(shape(ddrive) /= [fb%nt, fb%nmode])) return
        if (any(shape(dload) /= [fb%nt, fb%nmode])) return
        jv = 0
        allocate (t(fb%nt, fb%nmode))
        if (allocated(fb%exterior_3d)) then
            allocate (d3(fb%nt*fb%nmode), l3(fb%nt*fb%nmode))
            call ideal_fb_full_turn_traces(fb, ddrive, dload, d3, l3)
            t = reshape((matmul(fb%exterior_3d%s, d3) + l3)/(fb%mu0*fb%nfp), &
                [fb%nt, fb%nmode])
        else
            do m = 1, fb%nmode
                t(:, m) = ideal_fb_mode_weight(fb, m)/fb%mu0 &
                    *(matmul(fb%s_trace(:, :, harmonic(m)), ddrive(:, m)) + dload(:, m))
            end do
        end if
        do m = 1, fb%nmode
            do s = 1, fb%nt
                jv(iphi(fb, m, fb%trace(s))) = t(s, m)
            end do
        end do
        if (abs(df) > 0) then
            allocate (gx(3, fb%nmode, fb%ns), gf(fb%nmode, fb%ns), &
                fgx(3, fb%nmode, fb%ns), fgf(fb%nmode, fb%ns))
            call shell_energy(fb, fb%xs, fb%phi, ev, gx, gf, mins, f_gx=fgx, f_gf=fgf)
            do s = 1, fb%ns
                do m = 1, fb%nmode
                    do c = 1, 3
                        jv(ix(fb, c, m, s)) = jv(ix(fb, c, m, s)) - df*fgx(c, m, s)
                    end do
                    jv(iphi(fb, m, s)) = jv(iphi(fb, m, s)) - df*fgf(m, s)
                end do
            end do
        end if
        status = IDEAL_FB_OK
    end subroutine

    subroutine ideal_fb_traces_vjp(fb, v, drive_bar, load_bar, f_bar, status)
        ! Transpose of ideal_fb_traces_jvp: drive_bar, load_bar and f_bar
        ! with <v, J (dd, dl, df)> = <drive_bar, dd> + <load_bar, dl> + f_bar df.
        type(ideal_fb_t), intent(in) :: fb
        real(dp), intent(in) :: v(:)
        real(dp), intent(out) :: drive_bar(:, :), load_bar(:, :), f_bar
        integer, intent(out) :: status
        real(dp), allocatable :: vt(:, :), v3(:), gx(:, :, :), gf(:, :)
        real(dp), allocatable :: fgx(:, :, :), fgf(:, :)
        real(dp) :: ev, mins, cm
        integer :: m, s, c
        status = IDEAL_FB_BAD_INPUT
        if (.not. ideal_fb_exterior_ready(fb)) return
        if (size(v) /= fb%n) return
        if (any(shape(drive_bar) /= [fb%nt, fb%nmode])) return
        if (any(shape(load_bar) /= [fb%nt, fb%nmode])) return
        allocate (vt(fb%nt, fb%nmode))
        do m = 1, fb%nmode
            do s = 1, fb%nt
                vt(s, m) = v(iphi(fb, m, fb%trace(s)))
            end do
        end do
        if (allocated(fb%exterior_3d)) then
            v3 = reshape(vt, [fb%nt*fb%nmode])/(fb%mu0*fb%nfp)
            drive_bar = reshape(matmul(transpose(fb%exterior_3d%s), v3), &
                [fb%nt, fb%nmode])
            do m = 1, fb%nmode
                load_bar(:, m) = fb%nfp*ideal_fb_mode_weight(fb, m) &
                    *v3((m - 1)*fb%nt + 1:m*fb%nt)
            end do
        else
            do m = 1, fb%nmode
                cm = ideal_fb_mode_weight(fb, m)/fb%mu0
                drive_bar(:, m) = cm*matmul(transpose(fb%s_trace(:, :, harmonic(m))), &
                    vt(:, m))
                load_bar(:, m) = cm*vt(:, m)
            end do
        end if
        allocate (gx(3, fb%nmode, fb%ns), gf(fb%nmode, fb%ns), &
            fgx(3, fb%nmode, fb%ns), fgf(fb%nmode, fb%ns))
        call shell_energy(fb, fb%xs, fb%phi, ev, gx, gf, mins, f_gx=fgx, f_gf=fgf)
        f_bar = 0
        do s = 1, fb%ns
            do m = 1, fb%nmode
                do c = 1, 3
                    f_bar = f_bar - v(ix(fb, c, m, s))*fgx(c, m, s)
                end do
                f_bar = f_bar - v(iphi(fb, m, s))*fgf(m, s)
            end do
        end do
        status = IDEAL_FB_OK
    end subroutine

    subroutine ideal_fb_set_drive(fb, sample, status, closure)
        ! External traces from the field of everything outside Gamma_c plus
        ! the unit-circulation exterior part: drive a_m (single-valued
        ! external potential on the trace nodes) and load h_m = int (B.n)
        ! N_i R ds. n = 0: a = int (B - c grad theta/2pi).t ds from trace
        ! node 1; closure returns its mismatch after one turn (T m). n /= 0:
        ! a from B_phi = (1/R) dPhi/dphi; closure also returns the largest
        ! mismatch between d a/ds and B.t (curl-free check, T).
        ! The traces are those of ideal_fb_drive_traces for this one field.
        type(ideal_fb_t), intent(inout) :: fb
        procedure(ideal_fb_sampler) :: sample
        integer, intent(out) :: status
        real(dp), intent(out), optional :: closure
        type(sampler_basis_t) :: one
        real(dp), allocatable :: drive(:, :, :), load(:, :, :)
        real(dp) :: worst(1), f(1), spread_f(1)
        status = IDEAL_FB_BAD_INPUT
        if (present(closure)) closure = huge(1.0_dp)
        if (allocated(fb%exterior_3d)) return
        one%f => sample
        allocate (drive(fb%nt, fb%nmode, 1), load(fb%nt, fb%nmode, 1))
        call ideal_fb_drive_traces(fb, one, [fb%circulation], drive, load, worst, &
            f, spread_f, status)
        if (status /= IDEAL_FB_OK) return
        fb%drive = drive(:, :, 1); fb%load = load(:, :, 1)
        if (present(closure)) closure = worst(1)
    end subroutine

    integer function sampler_nfield(self)
        class(sampler_basis_t), intent(in) :: self
        sampler_nfield = 1
        if (.not. associated(self%f)) sampler_nfield = 0
    end function

    subroutine sampler_sample(self, rz, nfp, ntor, bcos, bsin, status)
        ! Serial: an ideal_fb_sampler may cache between calls.
        class(sampler_basis_t), intent(in) :: self
        real(dp), intent(in) :: rz(:, :)
        integer, intent(in) :: nfp, ntor
        real(dp), intent(out) :: bcos(:, 0:, :, :), bsin(:, 0:, :, :)
        integer, intent(out) :: status
        integer :: q, k
        status = IDEAL_FB_BAD_INPUT
        if (.not. associated(self%f) .or. nfp < 1) return
        do q = 1, size(rz, 2)
            do k = 0, ntor
                call self%f(rz(1, q), rz(2, q), k, bcos(:, k, 1, q), bsin(:, k, 1, q))
            end do
        end do
        status = IDEAL_FB_OK
    end subroutine

    subroutine basis_sample_points(self, x, nfp, ntor, b, status)
        ! Default Cartesian evaluation from the harmonic contract, serial
        ! (an ideal_fb_sampler may cache between calls).
        class(ideal_fb_field_basis_t), intent(in) :: self
        real(dp), intent(in) :: x(:, :)
        integer, intent(in) :: nfp, ntor
        real(dp), intent(out) :: b(:, :, :)
        integer, intent(out) :: status
        real(dp), allocatable :: bc(:, :, :, :), bs(:, :, :, :)
        real(dp) :: rz(2, 1), phi, bcyl(3)
        integer :: q, f, k, nf
        status = IDEAL_FB_BAD_INPUT
        nf = self%nfield()
        if (size(x, 1) /= 3 .or. nfp < 1 .or. ntor < 0 .or. nf < 1) return
        if (any(shape(b) /= [3, nf, size(x, 2)])) return
        allocate (bc(3, 0:ntor, nf, 1), bs(3, 0:ntor, nf, 1))
        do q = 1, size(x, 2)
            rz(:, 1) = [hypot(x(1, q), x(2, q)), x(3, q)]
            phi = atan2(x(2, q), x(1, q))
            call self%sample(rz, nfp, ntor, bc, bs, status)
            if (status /= IDEAL_FB_OK) return
            do f = 1, nf
                bcyl = bc(:, 0, f, 1)
                do k = 1, ntor
                    bcyl = bcyl + bc(:, k, f, 1)*cos(k*nfp*phi) &
                        + bs(:, k, f, 1)*sin(k*nfp*phi)
                end do
                b(:, f, q) = [bcyl(1)*cos(phi) - bcyl(2)*sin(phi), &
                    bcyl(1)*sin(phi) + bcyl(2)*cos(phi), bcyl(3)]
            end do
        end do
        status = IDEAL_FB_OK
    end subroutine

    subroutine ideal_fb_set_drive_3d(fb, sample, nphi, status, closure)
        ! ideal_fb_set_drive on the coupling surface of the plasma chart
        ! (ideal_fb_drive_traces_3d). With an attached R2 exterior its
        ! full-turn drive and load are replaced consistently.
        type(ideal_fb_t), intent(inout) :: fb
        procedure(ideal_fb_sampler) :: sample
        integer, intent(in) :: nphi
        integer, intent(out) :: status
        real(dp), intent(out), optional :: closure
        type(sampler_basis_t) :: one
        real(dp), allocatable :: drive(:, :, :), load(:, :, :)
        real(dp) :: worst(1), f(1), spread_f(1)
        status = IDEAL_FB_BAD_INPUT
        if (present(closure)) closure = huge(1.0_dp)
        one%f => sample
        allocate (drive(fb%nt, fb%nmode, 1), load(fb%nt, fb%nmode, 1))
        call ideal_fb_drive_traces_3d(fb, one, [fb%circulation], nphi, drive, load, &
            worst, f, spread_f, status)
        if (status /= IDEAL_FB_OK) return
        call ideal_fb_set_traces(fb, drive(:, :, 1), load(:, :, 1), fb%f_vac, status)
        if (present(closure)) closure = worst(1)
    end subroutine

    subroutine ideal_fb_drive_traces_3d(fb, basis, circulation, nphi, drive, load, &
                                        closure, f_vac, f_spread, status)
        ! The functionals of ideal_fb_drive_traces on the coupling surface
        ! Gamma_c = H(gamma x [0, 2 pi)) of the plasma chart H, gamma the
        ! reference boundary (the geometry of ideal_fb_surface), evaluated at
        ! Cartesian points X = H(gamma(t), phi) through the Cartesian basis
        ! contract. With x_t = dH (gamma', 0) and x_phi = dH e_3:
        !   load   nfp c_m h_{i,m} = int B.(x_phi x x_t) N_i(t) trig_m dt dphi
        !   n = 0  a(t) = int_0^t (<B.x_t> - c dtheta(gamma')/2 pi) dt
        !   n /= 0 a from B.x_phi = dPhi/dphi at the trace nodes
        !   F_vac  mean over the trace nodes of <B.x_phi>, f_spread its spread
        ! <.> and the harmonics use the nphi-point trapezoid rule over the
        ! full turn; t uses the Gauss rule of ideal_fb_drive_traces. With the
        ! cylindrical chart these are its meridional formulas at the same
        ! points. Outputs are per-period arrays as in ideal_fb_drive_traces
        ! (ideal_fb_full_turn_traces gives the R2 pairing).
        type(ideal_fb_t), intent(in) :: fb
        class(ideal_fb_field_basis_t), intent(in) :: basis
        real(dp), intent(in) :: circulation(:)
        integer, intent(in) :: nphi
        real(dp), intent(out) :: drive(:, :, :), load(:, :, :), closure(:)
        real(dp), intent(out) :: f_vac(:), f_spread(:)
        integer, intent(out) :: status
        integer, parameter :: ng = 10
        real(dp), allocatable :: q2(:, :), tq(:, :), wq(:, :), xc(:, :), xt(:, :)
        real(dp), allocatable :: xp(:, :), bf(:, :, :), ck(:, :), sk(:, :)
        real(dp), allocatable :: fnode(:), g(:)
        real(dp) :: gx(ng), gw(ng), w(0:3), dw(0:3), d1(3, 3), d2(3, 3, 3)
        real(dp) :: d3(3, 3, 3, 3), x3(3), coef(2), acc, run, worst, dmax, nn, speed
        real(dp) :: at(2), gc(0:fb%ntor), gs(0:fb%ntor)
        integer :: nf, i, j, q, k, col, b, p, nl, nd, nm, nb, np, l, nk, j0
        status = IDEAL_FB_BAD_INPUT
        nf = basis%nfield()
        np = nphi
        if (nf < 1 .or. size(circulation) /= nf) return
        if (any(shape(drive) /= [fb%nt, fb%nmode, nf])) return
        if (any(shape(load) /= [fb%nt, fb%nmode, nf])) return
        if (size(closure) /= nf .or. size(f_vac) /= nf .or. size(f_spread) /= nf) return
        if (np < 2 .or. np <= 2*fb%ntor*fb%nfp) return
        if (.not. allocated(fb%plasma%chart)) return
        call gauss_legendre(ng, gx, gw)
        gx = 0.5_dp*(gx + 1); gw = 0.5_dp*gw
        ! Base points on gamma: load quadrature, n = 0 path quadrature,
        ! trace nodes, edge midpoints; each is swept over nphi angles.
        nl = fb%nedge*ng; nd = fb%nedge*fb%order*ng; nm = fb%nedge
        nb = nl + nd + fb%nt + nm
        allocate (q2(2, nb), tq(2, nb), wq(0:fb%order, nl))
        do i = 1, fb%nedge
            do q = 1, ng
                p = (i - 1)*ng + q
                call edge_point(fb, i, gx(q), q2(:, p), tq(:, p), wq(:, p))
            end do
            do j = 1, fb%order
                do q = 1, ng
                    p = nl + ((i - 1)*fb%order + j - 1)*ng + q
                    call edge_point(fb, i, gx(q)*j/fb%order, q2(:, p), tq(:, p), &
                        w(0:fb%order))
                end do
            end do
            p = nl + nd + fb%nt + i
            call edge_point(fb, i, 0.5_dp, q2(:, p), tq(:, p), w(0:fb%order))
        end do
        do col = 1, fb%nt
            q2(:, nl + nd + col) = fb%reference(:, fb%trace(col))
            tq(:, nl + nd + col) = 0
        end do
        ! Trapezoid weights of the mean and the real harmonics (the DFT of
        ! coil_circuit_basis_t%sample, so both paths share their rounding).
        allocate (ck(np, 0:fb%ntor), sk(np, 0:fb%ntor))
        do j = 1, np
            do k = 0, fb%ntor
                nk = mod(k*fb%nfp*(j - 1), np)
                ck(j, k) = merge(1.0_dp, 2.0_dp, k == 0)*cos(2*pi*nk/np)/np
                sk(j, k) = 2*sin(2*pi*nk/np)/np
            end do
        end do
        allocate (xc(3, nb*np), xt(3, nb*np), xp(3, nb*np), bf(3, nf, nb*np))
        !$omp parallel do default(none) schedule(static) &
        !$omp shared(fb, q2, tq, nb, np, xc, xt, xp) private(p, j, l, x3, d1, d2, d3)
        do l = 1, nb*np
            p = (l - 1)/np + 1; j = l - (p - 1)*np
            call fb%plasma%chart%jet([q2(:, p), 2*pi*(j - 1)/np], x3, d1, d2, d3)
            xc(:, l) = x3
            xt(:, l) = matmul(d1(:, 1:2), tq(:, p))
            xp(:, l) = d1(:, 3)
        end do
        !$omp end parallel do
        if (.not. all(ieee_is_finite(xc)) .or. .not. all(ieee_is_finite(xt))) return
        if (.not. all(ieee_is_finite(xp))) return
        call basis%sample_points(xc, fb%nfp, fb%ntor, bf, status)
        if (status /= IDEAL_FB_OK) return
        status = IDEAL_FB_BAD_INPUT
        if (.not. all(ieee_is_finite(bf))) return
        status = IDEAL_FB_OK
        drive = 0; load = 0
        allocate (fnode(fb%nt), g(np))
        do b = 1, nf
            ! Loads: harmonics of B.(x_phi x x_t) at each Gauss point.
            do i = 1, fb%nedge
                do q = 1, ng
                    p = (i - 1)*ng + q; j0 = (p - 1)*np
                    do j = 1, np
                        g(j) = dot_product(bf(:, b, j0 + j), &
                            cross3(xp(:, j0 + j), xt(:, j0 + j)))
                    end do
                    call harmonics(g, gc, gs)
                    w(0:fb%order) = wq(:, p)
                    do j = 0, fb%order
                        col = (i - 1)*fb%order + j + 1
                        if (j == fb%order) col = mod(i, fb%nedge)*fb%order + 1
                        load(col, 1, b) = load(col, 1, b) + gw(q)*w(j)*gc(0)
                        do k = 1, fb%ntor
                            load(col, 2*k, b) = load(col, 2*k, b) + gw(q)*w(j)*gc(k)
                            load(col, 2*k + 1, b) = load(col, 2*k + 1, b) &
                                + gw(q)*w(j)*gs(k)
                        end do
                    end do
                end do
            end do
            ! n = 0 trace by integration of <B.x_t> along gamma.
            run = 0
            do i = 1, fb%nedge
                drive((i - 1)*fb%order + 1, 1, b) = run
                do j = 1, fb%order
                    acc = 0
                    do q = 1, ng
                        p = nl + ((i - 1)*fb%order + j - 1)*ng + q; j0 = (p - 1)*np
                        do l = 1, np
                            g(l) = dot_product(bf(:, b, j0 + l), xt(:, j0 + l))
                        end do
                        acc = acc + gw(q)*real(j, dp)/fb%order*(sum(ck(:, 0)*g) &
                            - circulation(b)*dot_product(theta_gradient(fb, q2(:, p)), &
                            tq(:, p)))
                    end do
                    if (j < fb%order) then
                        drive((i - 1)*fb%order + j + 1, 1, b) = run + acc
                    else
                        run = run + acc
                    end if
                end do
            end do
            worst = abs(run)
            ! n /= 0 traces and F_vac from B.x_phi at the trace nodes.
            do col = 1, fb%nt
                p = nl + nd + col; j0 = (p - 1)*np
                do l = 1, np
                    g(l) = dot_product(bf(:, b, j0 + l), xp(:, j0 + l))
                end do
                call harmonics(g, gc, gs)
                fnode(col) = gc(0)
                do k = 1, fb%ntor
                    nn = k*fb%nfp
                    drive(col, 2*k, b) = -gs(k)/nn
                    drive(col, 2*k + 1, b) = gc(k)/nn
                end do
            end do
            ! Curl-free check: d a/dt against <B.x_t> harmonics at midpoints.
            do i = 1, fb%nedge
                p = nl + nd + fb%nt + i; j0 = (p - 1)*np
                do l = 1, np
                    g(l) = dot_product(bf(:, b, j0 + l), xt(:, j0 + l))
                end do
                call harmonics(g, gc, gs)
                speed = norm2(tq(:, p))
                call edge_lagrange(fb%order, 0.5_dp, w(0:fb%order), dw(0:fb%order))
                dmax = 0
                do k = 1, fb%ntor
                    at = 0
                    do j = 0, fb%order
                        col = (i - 1)*fb%order + j + 1
                        if (j == fb%order) col = mod(i, fb%nedge)*fb%order + 1
                        at = at + dw(j)*drive(col, 2*k:2*k + 1, b)
                    end do
                    coef = [gc(k), gs(k)]
                    dmax = max(dmax, maxval(abs(at - coef))/speed)
                end do
                worst = max(worst, dmax)
            end do
            closure(b) = worst
            f_vac(b) = sum(fnode)/fb%nt
            f_spread(b) = maxval(abs(fnode - f_vac(b)))
        end do
    contains
        subroutine harmonics(values, c, s)
            real(dp), intent(in) :: values(:)
            real(dp), intent(out) :: c(0:), s(0:)
            integer :: kk
            do kk = 0, fb%ntor
                c(kk) = sum(ck(:, kk)*values); s(kk) = sum(sk(:, kk)*values)
            end do
            s(0) = 0
        end subroutine
    end subroutine

    pure function cross3(u, v) result(w)
        real(dp), intent(in) :: u(3), v(3)
        real(dp) :: w(3)
        w = [u(2)*v(3) - u(3)*v(2), u(3)*v(1) - u(1)*v(3), u(1)*v(2) - u(2)*v(1)]
    end function

    subroutine ideal_fb_drive_traces(fb, basis, circulation, drive, load, closure, &
                                     f_vac, f_spread, status)
        ! Traces of every field f of the basis on Gamma_c, linear in the
        ! field (formulas of ideal_fb_set_drive): drive(:, m, f), load(:, m, f)
        ! and the closure/curl-free mismatch closure(f). circulation(f) is the
        ! poloidal circulation removed from the n = 0 trace of field f (the
        ! unit-circulation filament part; zero for coils outside Gamma_c).
        ! f_vac(f) is the mean of R B_phi^(n=0) over the trace nodes (the
        ! toroidal period / 2 pi, T m), f_spread(f) its largest deviation
        ! there (zero for a field that is curl free across the trace).
        type(ideal_fb_t), intent(in) :: fb
        class(ideal_fb_field_basis_t), intent(in) :: basis
        real(dp), intent(in) :: circulation(:)
        real(dp), intent(out) :: drive(:, :, :), load(:, :, :), closure(:)
        real(dp), intent(out) :: f_vac(:), f_spread(:)
        integer, intent(out) :: status
        integer, parameter :: ng = 10
        real(dp), allocatable :: pl(:, :), tl(:, :), wl(:, :), pd(:, :), td(:, :)
        real(dp), allocatable :: pt(:, :), pm(:, :), tm(:, :)
        real(dp), allocatable :: cl(:, :, :, :), sl(:, :, :, :), cd(:, :, :, :)
        real(dp), allocatable :: sd(:, :, :, :), ct(:, :, :, :), st_(:, :, :, :)
        real(dp), allocatable :: cmid(:, :, :, :), smid(:, :, :, :)
        real(dp) :: gx(ng), gw(ng), x(2), tangent(2), w(0:3), bc(3), bs(3)
        real(dp) :: grad_t(2), acc, run, speed, worst, nn, at(2), dmax
        real(dp) :: dw(0:3), wv(0:3)
        integer :: nf, i, q, j, k, col, b, p, nl, nd
        status = IDEAL_FB_BAD_INPUT
        ! Meridional (R, Z) sampling: axisymmetric Gamma_c only. Other charts
        ! use ideal_fb_drive_traces_3d.
        if (.not. allocated(fb%s_trace)) return
        nf = basis%nfield()
        if (nf < 1 .or. size(circulation) /= nf) return
        if (any(shape(drive) /= [fb%nt, fb%nmode, nf])) return
        if (any(shape(load) /= [fb%nt, fb%nmode, nf])) return
        if (size(closure) /= nf .or. size(f_vac) /= nf .or. size(f_spread) /= nf) return
        call gauss_legendre(ng, gx, gw)
        gx = 0.5_dp*(gx + 1); gw = 0.5_dp*gw
        ! Sample points: load quadrature, n = 0 path quadrature, trace nodes
        ! and edge midpoints (curl-free check).
        nl = fb%nedge*ng; nd = fb%nedge*fb%order*ng
        allocate (pl(2, nl), tl(2, nl), wl(0:fb%order, nl), pd(2, nd), td(2, nd), &
            pt(2, fb%nt), pm(2, fb%nedge), tm(2, fb%nedge))
        do i = 1, fb%nedge
            do q = 1, ng
                p = (i - 1)*ng + q
                call edge_point(fb, i, gx(q), pl(:, p), tl(:, p), wl(:, p))
            end do
            do j = 1, fb%order
                do q = 1, ng
                    p = ((i - 1)*fb%order + j - 1)*ng + q
                    call edge_point(fb, i, gx(q)*j/fb%order, pd(:, p), td(:, p), &
                        w(0:fb%order))
                end do
            end do
            call edge_point(fb, i, 0.5_dp, pm(:, i), tm(:, i), w(0:fb%order))
        end do
        do col = 1, fb%nt
            pt(:, col) = fb%reference(:, fb%trace(col))
        end do
        allocate (cl(3, 0:fb%ntor, nf, nl), sl(3, 0:fb%ntor, nf, nl), &
            cd(3, 0:0, nf, nd), sd(3, 0:0, nf, nd), &
            ct(3, 0:fb%ntor, nf, fb%nt), st_(3, 0:fb%ntor, nf, fb%nt), &
            cmid(3, 0:fb%ntor, nf, fb%nedge), smid(3, 0:fb%ntor, nf, fb%nedge))
        call basis%sample(pl, fb%nfp, fb%ntor, cl, sl, status)
        if (status /= IDEAL_FB_OK) return
        call basis%sample(pd, fb%nfp, 0, cd, sd, status)
        if (status /= IDEAL_FB_OK) return
        call basis%sample(pt, fb%nfp, fb%ntor, ct, st_, status)
        if (status /= IDEAL_FB_OK) return
        if (fb%ntor > 0) then
            call basis%sample(pm, fb%nfp, fb%ntor, cmid, smid, status)
            if (status /= IDEAL_FB_OK) return
        end if
        status = IDEAL_FB_BAD_INPUT
        if (.not. all(ieee_is_finite(cl)) .or. .not. all(ieee_is_finite(sl))) return
        if (.not. all(ieee_is_finite(cd)) .or. .not. all(ieee_is_finite(ct))) return
        if (.not. all(ieee_is_finite(st_))) return
        status = IDEAL_FB_OK
        drive = 0; load = 0
        do b = 1, nf
            ! Loads.
            do i = 1, fb%nedge
                do q = 1, ng
                    p = (i - 1)*ng + q
                    x = pl(:, p); tangent = tl(:, p); w(0:fb%order) = wl(:, p)
                    do k = 0, fb%ntor
                        bc = cl(:, k, b, p); bs = sl(:, k, b, p)
                        do j = 0, fb%order
                            col = (i - 1)*fb%order + j + 1
                            if (j == fb%order) col = mod(i, fb%nedge)*fb%order + 1
                            if (k == 0) then
                                load(col, 1, b) = load(col, 1, b) + gw(q)*w(j)*x(1) &
                                    *(bc(1)*tangent(2) - bc(3)*tangent(1))
                            else
                                load(col, 2*k, b) = load(col, 2*k, b) + gw(q)*w(j) &
                                    *x(1)*(bc(1)*tangent(2) - bc(3)*tangent(1))
                                load(col, 2*k + 1, b) = load(col, 2*k + 1, b) &
                                    + gw(q)*w(j)*x(1)*(bs(1)*tangent(2) &
                                    - bs(3)*tangent(1))
                            end if
                        end do
                    end do
                end do
            end do
            ! n = 0 trace by integration along Gamma_c.
            run = 0
            do i = 1, fb%nedge
                drive((i - 1)*fb%order + 1, 1, b) = run
                do j = 1, fb%order
                    acc = 0
                    do q = 1, ng
                        p = ((i - 1)*fb%order + j - 1)*ng + q
                        x = pd(:, p); tangent = td(:, p); bc = cd(:, 0, b, p)
                        grad_t = theta_gradient(fb, x)
                        acc = acc + gw(q)*real(j, dp)/fb%order*(dot_product(bc([1, 3]) &
                            - circulation(b)*grad_t, tangent))
                    end do
                    if (j < fb%order) then
                        drive((i - 1)*fb%order + j + 1, 1, b) = run + acc
                    else
                        run = run + acc
                    end if
                end do
            end do
            worst = abs(run)
            ! n /= 0 traces from B_phi at the trace nodes, and the curl-free check.
            do k = 1, fb%ntor
                nn = k*fb%nfp
                do col = 1, fb%nt
                    x = pt(:, col)
                    drive(col, 2*k, b) = -x(1)*st_(2, k, b, col)/nn
                    drive(col, 2*k + 1, b) = x(1)*ct(2, k, b, col)/nn
                end do
                dmax = 0
                do i = 1, fb%nedge
                    tangent = tm(:, i); bc = cmid(:, k, b, i); bs = smid(:, k, b, i)
                    speed = norm2(tangent)
                    call edge_lagrange(fb%order, 0.5_dp, wv(0:fb%order), dw(0:fb%order))
                    at = 0
                    do j = 0, fb%order
                        col = (i - 1)*fb%order + j + 1
                        if (j == fb%order) col = mod(i, fb%nedge)*fb%order + 1
                        at = at + dw(j)*drive(col, 2*k:2*k + 1, b)
                    end do
                    dmax = max(dmax, abs(at(1)/speed - dot_product(bc([1, 3]), &
                        tangent)/speed), abs(at(2)/speed - dot_product(bs([1, 3]), &
                        tangent)/speed))
                end do
                worst = max(worst, dmax)
            end do
            closure(b) = worst
            f_vac(b) = sum(pt(1, :)*ct(2, 0, b, :))/fb%nt
            f_spread(b) = maxval(abs(pt(1, :)*ct(2, 0, b, :) - f_vac(b)))
        end do
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

    subroutine shell_point(fb, xs, phi, e, g, j, w, jac, gv, hv, want_hessian, fdot)
        ! Vacuum energy density and its derivatives with respect to
        ! v = (q(1:3), Q(k,i) at 3 + k + 3(i-1), dPhi~/ds_i at 12 + i).
        ! fdot: d gv / d F_vac (exact, from the generated Hessian leaf).
        type(ideal_fb_t), intent(in) :: fb
        real(dp), intent(in) :: xs(:, :, :), phi(:, :)
        integer, intent(in) :: e, g, j
        real(dp), intent(out) :: w, jac, gv(15), hv(15, 15)
        logical, intent(in) :: want_hessian
        real(dp), intent(out), optional :: fdot(15)
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
        if (.not. (want_hessian .or. present(fdot))) return
        call hessian_leaf(hflat)
        n = 0
        do c = 1, 12
            do r = 1, c
                n = n + 1
                hd(r, c) = hflat(n); hd(c, r) = hflat(n)
            end do
        end do
        if (present(fdot)) then
            ! gg depends on F_vac through F_vac*Q(3,:), and dv(9+i, 6+3(i-1))
            ! = F_vac: d gv/dF = (hd d(dx,gg)/dF) dv + gd d(dv)/dF.
            fdot = matmul(matmul(hd(:, 10:12), qd(3, :)), dv)
            do i = 1, 3
                fdot(6 + 3*(i - 1)) = fdot(6 + 3*(i - 1)) + gd(9 + i)
            end do
        end if
        if (.not. want_hessian) return
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

    subroutine shell_energy(fb, xs, phi, energy, gx, gf, min_jacobian, blocks, &
                            f_gx, f_gf)
        ! E_v interior part (with 1/mu0), gradients with respect to xs and
        ! Phi~, optionally the local Hessian blocks (shell-local order) and
        ! the exact derivatives f_gx, f_gf of gx, gf with respect to F_vac.
        type(ideal_fb_t), intent(in) :: fb
        real(dp), intent(in) :: xs(:, :, :), phi(:, :)
        real(dp), intent(out) :: energy, gx(:, :, :), gf(:, :), min_jacobian
        real(dp), intent(out), optional :: blocks(:, :, :)
        real(dp), intent(out), optional :: f_gx(:, :, :), f_gf(:, :)
        real(dp), allocatable :: local(:, :), rows(:, :), tmp(:, :), flocal(:, :)
        real(dp) :: w, jac, gv(15), hv(15, 15), weight, dphi, minjac, fd(15)
        integer :: e, g, j, nl, a, m, l
        logical :: want, want_f
        want = present(blocks)
        want_f = present(f_gx) .and. present(f_gf)
        nl = 4*fb%nmode*fb%nloc
        allocate (local(nl, fb%ncell))
        if (want_f) then
            allocate (flocal(nl, fb%ncell))
        else
            allocate (flocal(1, 1))
        end if
        dphi = 2*pi/(fb%nfp*fb%plasma%nphi)
        energy = 0; minjac = huge(1.0_dp)
        !$omp parallel default(shared) &
        !$omp private(e, g, j, w, jac, gv, hv, weight, rows, tmp, fd)
        allocate (rows(15, nl), tmp(15, nl))
        !$omp do reduction(+: energy) reduction(min: minjac) schedule(static)
        do e = 1, fb%ncell
            local(:, e) = 0
            if (want) blocks(:, :, e) = 0
            if (want_f) flocal(:, e) = 0
            do j = 1, fb%plasma%nphi
                do g = 1, fb%plasma%nq
                    if (want_f) then
                        call shell_point(fb, xs, phi, e, g, j, w, jac, gv, hv, want, fd)
                    else
                        call shell_point(fb, xs, phi, e, g, j, w, jac, gv, hv, want)
                    end if
                    minjac = min(minjac, jac/fb%det(g, e))
                    weight = fb%plasma%weight(g)*dphi
                    energy = energy + weight*w
                    call local_rows(fb, g, j, rows)
                    local(:, e) = local(:, e) + weight*matmul(gv, rows)
                    if (want_f) flocal(:, e) = flocal(:, e) + weight*matmul(fd, rows)
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
        if (want_f) then
            f_gx = 0; f_gf = 0
        end if
        do e = 1, fb%ncell
            do a = 1, fb%nloc
                do m = 1, fb%nmode
                    l = 4*(m - 1) + 4*fb%nmode*(a - 1)
                    gx(:, m, fb%conn(a, e)) = gx(:, m, fb%conn(a, e)) &
                        + local(l + 1:l + 3, e)
                    gf(m, fb%conn(a, e)) = gf(m, fb%conn(a, e)) + local(l + 4, e)
                    if (want_f) then
                        f_gx(:, m, fb%conn(a, e)) = f_gx(:, m, fb%conn(a, e)) &
                            + flocal(l + 1:l + 3, e)
                        f_gf(m, fb%conn(a, e)) = f_gf(m, fb%conn(a, e)) &
                            + flocal(l + 4, e)
                    end if
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
        ! Gradient of L with respect to (plasma positions, xs, lambda, Phi~).
        type(ideal_fb_t), intent(in) :: fb
        real(dp), intent(out) :: grad(:)
        integer, intent(out) :: status
        real(dp), intent(out), optional :: lagrangian, min_jacobian
        real(dp), allocatable :: gfull(:), gx(:, :, :), gf(:, :), res(:, :, :), u(:)
        real(dp), allocatable :: cx(:, :, :), exterior_load(:)
        real(dp) :: wp, ev, ec, minp, mins, cm
        integer :: st, s, m, c, p, a, k
        status = IDEAL_FB_BAD_INPUT
        if (.not. ideal_fb_exterior_ready(fb)) return
        allocate (gfull(fb%plasma%ndof), gx(3, fb%nmode, fb%ns), gf(fb%nmode, fb%ns), &
            res(3, fb%nmode, fb%ns), cx(3, fb%nmode, fb%ns))
        status = IDEAL_FB_OK
        call ideal_map_energy(fb%plasma, fb%plasma%x, wp, st, gfull, minp)
        if (st /= IDEAL_MAP_OK) status = IDEAL_FB_FOLDED
        do a = 1, fb%plasma%nodes
            do m = 1, fb%nmode
                do c = 1, 2
                    grad(pidx(fb, c, m, a)) = gfull(ideal_map_dof(fb%plasma, c, m, a))
                end do
            end do
        end do
        call shell_energy(fb, fb%xs, fb%phi, ev, gx, gf, mins)
        if (.not. mins > 0 .or. .not. ieee_is_finite(ev)) status = IDEAL_FB_FOLDED
        ! Gamma_c terms.
        ec = 0
        if (allocated(fb%exterior_3d)) then
            u = reshape(transpose(fb%phi(:,fb%trace)),[fb%nt*fb%nmode])
            allocate(exterior_load(size(u)))
            call gamma_c_3d_energy(fb%exterior_3d,u,fb%drive_3d,fb%load_3d, &
                fb%mu0,ec,exterior_load,st)
            if (st /= 0) then
                status = IDEAL_FB_BEM_FAILED
                return
            end if
            ec = ec/fb%nfp
            gf(:,fb%trace) = gf(:,fb%trace) + &
                transpose(reshape(exterior_load,[fb%nt,fb%nmode]))/fb%nfp
        else
            do m = 1, fb%nmode
                cm = ideal_fb_mode_weight(fb, m)/fb%mu0
                k = harmonic(m)
                u = fb%phi(m, fb%trace) - fb%drive(:, m)
                ec = ec + cm*(0.5_dp*dot_product(u, matmul(fb%s_trace(:, :, k), u)) &
                    - dot_product(fb%load(:, m), fb%phi(m, fb%trace)))
                gf(m, fb%trace) = gf(m, fb%trace) + cm*(matmul(fb%s_trace(:, :, k), u) &
                    - fb%load(:, m))
            end do
        end if
        call constraint_residual(fb, res)
        ! C^T lambda on xs and on the plasma Gamma_p dofs.
        cx = 0
        do s = 1, fb%ns
            select case (fb%role(s))
            case (1)
                cx(:, :, s) = cx(:, :, s) + fb%lam(:, :, s)
                a = fb%plasma_node(s)
                do m = 1, fb%nmode
                    do c = 1, 2
                        p = pidx(fb, c, m, a)
                        grad(p) = grad(p) - fb%lam(c, m, s)
                    end do
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
        if (allocated(fb%shift)) grad = grad - fb%shift
        if (present(lagrangian)) lagrangian = wp - ev - ec + sum(fb%lam*res)
        if (present(min_jacobian)) min_jacobian = min(minp, mins)
    end subroutine

    subroutine ideal_fb_hessian(fb, rowptr, colidx, val)
        ! Assembled KKT matrix of L (symmetric, indefinite), CSR with sorted
        ! columns (equal to its CSC form, which is what fortsparse reads).
        type(ideal_fb_t), intent(in) :: fb
        integer, allocatable, intent(out) :: rowptr(:), colidx(:)
        real(dp), allocatable, intent(out) :: val(:)
        type(csr_matrix_t) :: hp, hs
        real(dp), allocatable :: blocks(:, :, :), gx(:, :, :), gf(:, :)
        integer, allocatable :: sconn(:, :), ci(:), cj(:)
        real(dp), allocatable :: cv(:)
        real(dp) :: ev, mins, cm
        integer :: e, a, m, c, l, nl, cnt, i, k, p, s, t1, t2, nsb, cap, m2
        nl = 4*fb%nmode*fb%nloc
        call plasma_hessian(fb, hp)
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
        cap = fb%nmode*fb%nt**2
        if (allocated(fb%exterior_3d)) cap = (fb%nmode*fb%nt)**2
        cap = size(hp%val) + size(hs%val) + cap &
            + 2*3*fb%nmode*(size(fb%lap_val) + 2*fb%ns) + fb%np
        allocate (ci(cap), cj(cap), cv(cap))
        cnt = 0
        do i = 1, fb%np
            do p = hp%rowptr(i), hp%rowptr(i + 1) - 1
                call push(i, hp%colidx(p), hp%val(p))
            end do
        end do
        do i = 1, nsb
            do p = hs%rowptr(i), hs%rowptr(i + 1) - 1
                call push(shell_global(i), shell_global(hs%colidx(p)), -hs%val(p))
            end do
        end do
        if (allocated(fb%exterior_3d)) then
            do m2 = 1, fb%nmode
                do m = 1, fb%nmode
                    do t2 = 1, fb%nt
                        do t1 = 1, fb%nt
                            call push(iphi(fb,m,fb%trace(t1)), &
                                iphi(fb,m2,fb%trace(t2)), &
                                -fb%exterior_3d%s(t1+(m-1)*fb%nt, &
                                t2+(m2-1)*fb%nt)/(fb%mu0*fb%nfp))
                        end do
                    end do
                end do
            end do
        else
            do m = 1, fb%nmode
                cm = ideal_fb_mode_weight(fb, m)/fb%mu0
                k = harmonic(m)
                do t2 = 1, fb%nt
                    do t1 = 1, fb%nt
                        call push(iphi(fb,m,fb%trace(t1)), &
                            iphi(fb,m,fb%trace(t2)),-cm*fb%s_trace(t1,t2,k))
                    end do
                end do
            end do
        end if
        do s = 1, fb%ns
            do m = 1, fb%nmode
                do c = 1, 3
                    select case (fb%role(s))
                    case (1)
                        ! Gamma_p tie: shell (q1, q2) = plasma (q1, q2), shell
                        ! omega = 0 (position gauge of the plasma).
                        call pair(il(fb, c, m, s), ix(fb, c, m, s), 1.0_dp)
                        if (c < 3) call pair(il(fb, c, m, s), &
                            pidx(fb, c, m, fb%plasma_node(s)), -1.0_dp)
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
            ! Constraint entry and its transpose.
            integer, intent(in) :: r, col
            real(dp), intent(in) :: v
            call push(r, col, v)
            call push(col, r, v)
        end subroutine
    end subroutine

    subroutine plasma_hessian(fb, matrix)
        ! Hessian of the plasma energy W_p in the position gauge (q1, q2 per
        ! node and mode, omega = 0): the second variation of N1's frozen-in
        ! density (same generated leaves, chart chain rule and quadrature as
        ! ideal_map_hessian), kept in full coefficients instead of N1's
        ! (transversal, omega) reduction.
        type(ideal_fb_t), intent(in) :: fb
        type(csr_matrix_t), intent(out) :: matrix
        real(dp), allocatable :: blocks(:, :, :), rows(:, :), tmp(:, :)
        integer, allocatable :: conn(:, :)
        real(dp) :: hv(12, 12), weight, dphi, f(0:3)
        integer :: e, g, j, a, m, c, l, nl
        nl = 2*fb%nmode*fb%nloc
        allocate (blocks(nl, nl, fb%plasma%cells), conn(nl, fb%plasma%cells))
        dphi = 2*pi/(fb%nfp*fb%plasma%nphi)
        !$omp parallel default(shared) private(e, g, j, a, m, c, l, hv, weight, f, &
        !$omp rows, tmp)
        allocate (rows(12, nl), tmp(12, nl))
        !$omp do schedule(static)
        do e = 1, fb%plasma%cells
            blocks(:, :, e) = 0
            do j = 1, fb%plasma%nphi
                do g = 1, fb%plasma%nq
                    call plasma_point_hessian(fb%plasma, e, g, j, hv)
                    weight = fb%plasma%weight(g)*dphi
                    rows = 0
                    do a = 1, fb%nloc
                        do m = 1, fb%nmode
                            f(0) = fb%plasma%shape(a, g)*fb%plasma%trig(m, j)
                            f(1) = fb%plasma%dshape(1, a, g)*fb%plasma%trig(m, j)
                            f(2) = fb%plasma%dshape(2, a, g)*fb%plasma%trig(m, j)
                            f(3) = fb%plasma%shape(a, g)*fb%plasma%dtrig(m, j)
                            do c = 1, 2
                                l = c + 2*(m - 1) + 2*fb%nmode*(a - 1)
                                rows(c, l) = f(0)
                                rows(3 + c, l) = f(1)
                                rows(6 + c, l) = f(2)
                                rows(9 + c, l) = f(3)
                            end do
                        end do
                    end do
                    tmp = matmul(hv, rows)
                    blocks(:, :, e) = blocks(:, :, e) + weight*matmul(transpose(rows), tmp)
                end do
            end do
        end do
        !$omp end do
        deallocate (rows, tmp)
        !$omp end parallel
        do e = 1, fb%plasma%cells
            do a = 1, fb%nloc
                do m = 1, fb%nmode
                    do c = 1, 2
                        conn(c + 2*(m - 1) + 2*fb%nmode*(a - 1), e) = &
                            pidx(fb, c, m, fb%plasma%conn(a, e))
                    end do
                end do
            end do
        end do
        call csr_build_pattern(fb%np, conn, .false., matrix)
        call csr_assemble_blocks(matrix, conn, blocks)
    end subroutine

    subroutine plasma_point_hessian(map, e, g, j, hv)
        ! Second derivatives of N1's energy density w(DX) with respect to
        ! v = (q(1:3), Q(k,i) at 3 + k + 3(i-1)), DX = dH(q) Q (chart chain
        ! rule as in ideal_map's local_derivatives).
        type(ideal_map_t), intent(in) :: map
        integer, intent(in) :: e, g, j
        real(dp), intent(out) :: hv(12, 12)
        real(dp) :: q(3), qd(3, 3), pos(3), d1(3, 3), d2(3, 3, 3), d3(3, 3, 3, 3)
        real(dp) :: dx(3, 3), gd(3, 3), hd(9, 9), dvx(9, 12), t(45), tt, dt
        integer :: a, i, k, m, n, r, c, node
        q = 0; qd = 0
        do a = 1, map%nloc
            node = map%conn(a, e)
            do m = 1, map%nmode
                tt = map%trig(m, j); dt = map%dtrig(m, j)
                q = q + map%x(:, m, node)*map%shape(a, g)*tt
                qd(:, 1) = qd(:, 1) + map%x(:, m, node)*map%dshape(1, a, g)*tt
                qd(:, 2) = qd(:, 2) + map%x(:, m, node)*map%dshape(2, a, g)*tt
                qd(:, 3) = qd(:, 3) + map%x(:, m, node)*map%shape(a, g)*dt
            end do
        end do
        q(3) = q(3) + map%phi(j); qd(3, 3) = qd(3, 3) + 1
        call map%chart%jet(q, pos, d1, d2, d3)
        dx = matmul(d1, qd)
        call ideal_energy_gradient(dx(1, 1), dx(2, 1), dx(3, 1), dx(1, 2), &
            dx(2, 2), dx(3, 2), dx(1, 3), dx(2, 3), dx(3, 3), map%b(1, g, e), &
            map%b(2, g, e), map%b(3, g, e), map%p(g, e), map%mu0, map%sigma, &
            gd(1, 1), gd(2, 1), gd(3, 1), gd(1, 2), gd(2, 2), gd(3, 2), gd(1, 3), &
            gd(2, 3), gd(3, 3))
        call ideal_energy_hessian(dx(1, 1), dx(2, 1), dx(3, 1), dx(1, 2), &
            dx(2, 2), dx(3, 2), dx(1, 3), dx(2, 3), dx(3, 3), map%b(1, g, e), &
            map%b(2, g, e), map%b(3, g, e), map%p(g, e), map%mu0, map%sigma, &
            t(1), t(2), t(3), t(4), t(5), t(6), t(7), t(8), t(9), t(10), t(11), &
            t(12), t(13), t(14), t(15), t(16), t(17), t(18), t(19), t(20), t(21), &
            t(22), t(23), t(24), t(25), t(26), t(27), t(28), t(29), t(30), t(31), &
            t(32), t(33), t(34), t(35), t(36), t(37), t(38), t(39), t(40), t(41), &
            t(42), t(43), t(44), t(45))
        dvx = 0
        do i = 1, 3
            do a = 1, 3
                r = a + 3*(i - 1)
                do m = 1, 3
                    dvx(r, m) = dot_product(d2(a, :, m), qd(:, i))
                end do
                do k = 1, 3
                    dvx(r, 3 + k + 3*(i - 1)) = d1(a, k)
                end do
            end do
        end do
        n = 0
        do c = 1, 9
            do r = 1, c
                n = n + 1
                hd(r, c) = t(n); hd(c, r) = t(n)
            end do
        end do
        hv = matmul(transpose(dvx), matmul(hd, dvx))
        do m = 1, 3
            do n = 1, 3
                do i = 1, 3
                    do a = 1, 3
                        hv(m, n) = hv(m, n) + gd(a, i)*dot_product(d3(a, :, m, n), &
                            qd(:, i))
                    end do
                end do
            end do
            do i = 1, 3
                do k = 1, 3
                    r = 3 + k + 3*(i - 1)
                    hv(m, r) = hv(m, r) + dot_product(gd(:, i), d2(:, k, m))
                    hv(r, m) = hv(m, r)
                end do
            end do
        end do
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
        integer :: s, m, c, a
        do a = 1, fb%plasma%nodes
            do m = 1, fb%nmode
                do c = 1, 2
                    fb%plasma%x(c, m, a) = fb%plasma%x(c, m, a) &
                        + alpha*du(pidx(fb, c, m, a))
                end do
            end do
        end do
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

    subroutine slide_directions(fb)
        ! Gauge of the n = 0 relabelling null. A map of the reference onto
        ! itself that keeps chi and the toroidal-flux 2-form b3 dR dZ moves
        ! points with v = f(chi) (-d chi/dZ, d chi/dR)/b3 for any f
        ! (div(b3 v) = 0, v.grad chi = 0): fields, pressure and surfaces are
        ! unchanged, an exact continuum null that is near-null discretely.
        ! It is fixed by the b3-weighted slide moments of the displacement u,
        !   M_l(u) = int phi_l(rho) b3 v[f = 1] . u dR dZ
        !          = int phi_l(rho) (DG b_s) . u ds   (b_s the frozen poloidal
        !            flux density, master coordinates, reference quadrature),
        ! with the C1 partition of unity phi_l = cos^2(pi/2 (L rho - l)),
        ! l = 0..L, rho = (1 - T)^1/2 and one layer per element layer, as in
        ! the nested omega gauge of kin6d_ideal_map. Integrating by parts,
        ! M_l vanishes for every curl-free u (translations, elongation), and
        ! a smooth f restores the moments for any continuum equilibrium: the
        ! gauge selects the labelling only (to discretisation accuracy) and
        ! leaves the geometry and the fields free. The f(chi) left free
        ! oscillate within a layer and are stiff.
        type(ideal_fb_t), intent(inout) :: fb
        real(dp) :: dg(2, 2), flux(2), rho, x, w
        integer :: e, g, l, nl, k
        nl = max(1, nint(sqrt(fb%plasma%cells*sqrt(3.0_dp)/(4*pi))))
        allocate (fb%slide(2, fb%plasma%nodes, nl + 1))
        fb%slide = 0
        associate (map => fb%plasma)
            do e = 1, map%cells
                do g = 1, map%nq
                    rho = nl*sqrt(max(0.0_dp, &
                        1 - dot_product(map%axis_t(map%conn(:, e)), map%shape(:, g))))
                    dg = matmul(map%reference(:, map%conn(:, e)), &
                        transpose(map%dshape(:, :, g)))
                    flux = matmul(dg, map%b(1:2, g, e))*map%weight(g)
                    do l = max(0, floor(rho)), min(nl, ceiling(rho))
                        x = rho - l
                        if (abs(x) >= 1) cycle
                        w = cos(pi*x/2)**2
                        do k = 1, map%nloc
                            fb%slide(:, map%conn(k, e), l + 1) = &
                                fb%slide(:, map%conn(k, e), l + 1) &
                                + w*map%shape(k, g)*flux
                        end do
                    end do
                end do
            end do
        end associate
        do l = 1, nl + 1
            fb%slide(:, :, l) = fb%slide(:, :, l)/norm2(fb%slide(:, :, l))
        end do
    end subroutine

    subroutine ideal_fb_slide_system(fb, rowptr, colidx, val, scale, rows, cols, &
                                     matrix)
        ! Newton matrix bordered by the slide gauge: [H, s G^T; s G, 0], with
        ! G the slide directions on the n = 0 plasma rows (fb%n + size(slide, 3)
        ! rows; the gauge rows follow the KKT unknowns).
        type(ideal_fb_t), intent(in) :: fb
        integer, intent(in) :: rowptr(:), colidx(:)
        real(dp), intent(in) :: val(:), scale
        integer, allocatable, intent(out) :: rows(:), cols(:)
        real(dp), allocatable, intent(out) :: matrix(:)
        integer, allocatable :: ci(:), cj(:)
        real(dp), allocatable :: cv(:)
        integer :: i, p, k, a, c, m, nl
        nl = size(fb%slide, 3)
        m = size(val) + 2*count(abs(fb%slide) > 0)
        allocate (ci(m), cj(m), cv(m))
        m = 0
        do i = 1, fb%n
            do p = rowptr(i), rowptr(i + 1) - 1
                m = m + 1
                ci(m) = i; cj(m) = colidx(p); cv(m) = val(p)
            end do
        end do
        do k = 1, nl
            do a = 1, fb%plasma%nodes
                do c = 1, 2
                    if (.not. abs(fb%slide(c, a, k)) > 0) cycle
                    ci(m + 1:m + 2) = [fb%n + k, pidx(fb, c, 1, a)]
                    cj(m + 1:m + 2) = [pidx(fb, c, 1, a), fb%n + k]
                    cv(m + 1:m + 2) = scale*fb%slide(c, a, k)
                    m = m + 2
                end do
            end do
        end do
        call coo_to_csr(fb%n + nl, ci(:m), cj(:m), cv(:m), rows, cols, matrix)
    end subroutine

    subroutine project_slide(fb, w, g, multiplier)
        ! Remove the slide component of the plasma gradient in the merit
        ! metric: g <- g - G^T mu, mu = argmin sum w (g + G^T mu)^2. This is
        ! the gauge multiplier of the bordered Newton system, so the merit
        ! measures only the physical (gauge-invariant) residual. multiplier
        ! returns G^T-coefficients b with g_removed = G^T b (FB-C5: b -> 0 at
        ! an equilibrium in the gauge).
        type(ideal_fb_t), intent(in) :: fb
        real(dp), intent(in) :: w(:)
        real(dp), intent(inout) :: g(:)
        real(dp), intent(out), optional :: multiplier(:)
        real(dp) :: a(size(fb%slide, 3), size(fb%slide, 3)), b(size(fb%slide, 3))
        integer :: i, j, k, n, c, nl
        nl = size(fb%slide, 3)
        a = 0; b = 0
        do n = 1, fb%plasma%nodes
            do c = 1, 2
                i = pidx(fb, c, 1, n)
                b = b + w(i)*fb%slide(c, n, :)*g(i)
                do k = 1, nl
                    a(:, k) = a(:, k) + w(i)*fb%slide(c, n, :)*fb%slide(c, n, k)
                end do
            end do
        end do
        ! Cholesky solve a mu = b (a is SPD: the slide directions are independent).
        do j = 1, nl
            a(j, j) = sqrt(a(j, j) - sum(a(j, :j - 1)**2))
            do i = j + 1, nl
                a(i, j) = (a(i, j) - sum(a(i, :j - 1)*a(j, :j - 1)))/a(j, j)
            end do
        end do
        do i = 1, nl
            b(i) = (b(i) - sum(a(i, :i - 1)*b(:i - 1)))/a(i, i)
        end do
        do i = nl, 1, -1
            b(i) = (b(i) - sum(a(i + 1:, i)*b(i + 1:)))/a(i, i)
        end do
        do n = 1, fb%plasma%nodes
            do c = 1, 2
                i = pidx(fb, c, 1, n)
                g(i) = g(i) - dot_product(fb%slide(c, n, :), b)
            end do
        end do
        if (present(multiplier)) multiplier = b
    end subroutine

    subroutine ideal_fb_solve(fb, tolerance, max_iterations, report, status, verbose, &
                              step_tolerance, residual_tolerance)
        ! Newton on the slide-gauged KKT system of L (symmetric, indefinite:
        ! negative vacuum block, multipliers, and with fixed coil currents a
        ! negative vertical direction of the reduced plasma Hessian). The n = 0
        ! relabelling null is fixed by the bordered slide moments
        ! (slide_directions); the merit is taken with their multiplier removed.
        ! Globalisation:
        !   - merit = sum g_i^2/|H_ii| over primal rows, g_i^2 |H_xx| over
        !     constraint rows; a step is accepted when the merit decreases and
        !     every Jacobian stays positive (never on the energy, which is a
        !     saddle);
        !   - backtracking (alpha = 1, 1/2, .., 1/128) along the Newton
        !     direction; only when that fails, Levenberg damping mu |H_ii| on
        !     the plasma rows (never on the vacuum or multiplier blocks),
        !     reduced by 100 after every full step;
        !   - trust radius: displacement <= 2 % of the Gamma_c extent.
        ! With a positive residual_tolerance, convergence requires the residual
        ! sqrt(merit/(sigma V)) <= residual_tolerance (sigma V = energy_scale,
        ! the poloidal-field and pressure energy: toroidal field removed from
        ! the normalisation). Otherwise use relative merit or the legacy floor:
        ! the merit stagnates
        ! (less than a factor 2 over 5 accepted steps) with displacement
        ! increments below step_tolerance (m) or no acceptable step remains
        ! below residual 1e-6.
        ! max_iterations caps accepted steps, not residual evaluations: zero
        ! checks the supplied state, and the final accepted step is checked too.
        type(ideal_fb_t), intent(inout) :: fb
        real(dp), intent(in) :: tolerance
        integer, intent(in) :: max_iterations
        type(ideal_fb_report_t), intent(out) :: report
        integer, intent(out) :: status
        logical, intent(in), optional :: verbose
        real(dp), intent(in), optional :: step_tolerance, residual_tolerance
        integer, allocatable :: rowptr(:), colidx(:), dpos(:)
        real(dp), allocatable :: val(:), dval(:), grad(:), du(:), diag(:), wts(:)
        real(dp), allocatable :: rhs(:), x0(:, :, :), xs0(:, :, :), phi0(:, :)
        real(dp), allocatable :: lam0(:, :, :)
        type(sparse_solver_t) :: solver
        type(fortsparse_status_t) :: sst
        real(dp) :: merit, tmerit, alpha, cap_x, cap_w, disp, rot, minjac, t0
        real(dp) :: step_tol, mu, scale_p, res_tol, history(0:5), tangent_damping
        integer :: it, i, p, st, tries, back, nback
        integer, allocatable :: tpos(:, :)
        real(dp), allocatable :: trial_du(:)
        type(sparse_solver_t) :: csolver
        type(ideal_fb_hypotheses_t) :: hyp
        logical :: cready
        logical :: talk, accepted, factored, kept, reused_factor
        real(dp) :: kz, fz, fzr, ftry
        integer, allocatable :: srows(:), scols(:)
        real(dp), allocatable :: smat(:), sol(:)
        real(dp) :: sscale
        integer :: nslide
        report = ideal_fb_report_t()
        report%step = 0
        status = IDEAL_FB_BAD_INPUT
        if (max_iterations < 0) return
        if (.not. ideal_fb_exterior_ready(fb)) return
        kz = 0
        nslide = size(fb%slide, 3)
        allocate (sol(fb%n + nslide))
        talk = .false.
        if (present(verbose)) talk = verbose
        step_tol = 1.0e-5_dp
        if (present(step_tolerance)) step_tol = step_tolerance
        res_tol = 0
        if (present(residual_tolerance)) res_tol = residual_tolerance
        history = huge(1.0_dp)
        cap_x = 0.02_dp*maxval(maxval(fb%reference, dim=2) - minval(fb%reference, dim=2))
        cap_w = 0.5_dp
        mu = 1.0e-8_dp
        tangent_damping = fb%tangential_damping
        allocate (grad(fb%n), du(fb%n), diag(fb%n), wts(fb%n), rhs(fb%n), dpos(fb%np))
        factored = .false.; cready = .false.
        call ideal_fb_gradient(fb, grad, st, min_jacobian=minjac)
        if (st /= IDEAL_FB_OK) then
            status = IDEAL_FB_FOLDED; return
        end if
        status = IDEAL_FB_NOT_CONVERGED
        report%min_jacobian = minjac
        if (fb%rigid_vertical) then
            call ideal_fb_rigid_stiffness(fb, 2, 1.0e-3_dp, fz, kz, st)
            if (st /= IDEAL_FB_OK .or. .not. kz < 0) kz = 0
            if (talk) print '(a,es12.4)', '   rigid vertical stiffness K_Z ', kz
        end if
        do it = 0, max_iterations
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
            scale_p = maxval(diag(1:fb%np))
            sscale = sum(diag(1:fb%np))/fb%np
            if (.not. cready) then
                ! C^T (shell-map rows x multiplier columns) depends on the
                ! reference only: factorise once; after every trial step the
                ! multipliers are re-solved exactly, so the mesh-motion rows
                ! (toroidal-field-scale forces on shell nodes, no physics)
                ! never enter the merit.
                call factor_constraint(fb, rowptr, colidx, val, csolver, cready)
                if (cready) then
                    call project_multipliers(fb, csolver, grad)
                    call merit_weights(fb, diag, wts)
                end if
            end if
            if (fb%tangential_damping > 0) &
                call tangential_positions(fb, rowptr, colidx, tpos)
            call merit_weights(fb, diag, wts)
            if (it == 0) fzr = rigid_force(grad)
            call project_slide(fb, wts, grad)
            merit = sum(wts*grad**2)
            if (it == 0) report%merit0 = merit
            report%merit = merit
            report%residual = sqrt(merit/fb%energy_scale)
            if (talk) print '(a,i3,a,es10.3,a,es10.3,a,es10.3,a,es9.2)', 'fb newton', &
                it, ' merit/m0 ', sqrt(merit/max(report%merit0, tiny(1.0_dp))), &
                ' res/sigmaV ', &
                report%residual, ' step ', report%step, ' mu ', mu
            if (talk) call boundary_report(fb)
            if (talk) print '(a,4es10.2)', '   merit/sigmaV by block (plasma, shell map, &
                &constraint, potential)', sum(wts(1:fb%np)*grad(1:fb%np)**2) &
                /fb%energy_scale, sum(wts(fb%np + 1:fb%np + fb%nx)*grad(fb%np + 1:fb%np &
                + fb%nx)**2)/fb%energy_scale, sum(wts(fb%np + fb%nx + 1:fb%np + 2*fb%nx) &
                *grad(fb%np + fb%nx + 1:fb%np + 2*fb%nx)**2)/fb%energy_scale, &
                sum(wts(fb%np + 2*fb%nx + 1:)*grad(fb%np + 2*fb%nx + 1:)**2)/fb%energy_scale
            if (talk .and. fb%monitor_hypotheses) then
                call ideal_fb_hypotheses(fb, hyp, st)
                call ideal_fb_hypotheses_print(hyp, it)
            end if
            if (res_tol > 0) then
                if (report%residual <= res_tol) then
                    status = IDEAL_FB_OK; exit
                end if
            else if (sqrt(merit) <= tolerance*sqrt(report%merit0)) then
                status = IDEAL_FB_OK; exit
            end if
            if (it == max_iterations) exit
            history = eoshift(history, 1, merit)
            if (history(0) < huge(1.0_dp)) then
                if (history(5) > 0.5_dp*history(0) .and. report%step < step_tol &
                    .and. res_tol <= 0) then
                    report%floor = .true.; status = IDEAL_FB_OK; exit
                end if
            end if
            accepted = .false.
            x0 = fb%plasma%x; xs0 = fb%xs; phi0 = fb%phi; lam0 = fb%lam
            do tries = 1, 16
                dval = val
                do i = 1, fb%np
                    dval(dpos(i)) = dval(dpos(i)) + mu*(diag(i) + 1.0e-12_dp*scale_p)
                end do
                if (fb%tangential_damping > 0) &
                    call add_tangential(fb, diag, tpos, dval, tangent_damping)
                call ideal_fb_slide_system(fb, rowptr, colidx, dval, sscale, &
                                           srows, scols, smat)
                t0 = wall()
                rhs = -grad
                ! The first try of a Newton step reuses the factor of the
                ! previous step as the preconditioner of FGMRES on the current
                ! KKT matrix (accepted only at full solve accuracy); damped
                ! retries and Krylov failures refactorise as before.
                reused_factor = .false.
                if (factored .and. tries == 1) &
                    call preconditioned_solve(solver, srows, scols, smat, &
                        [rhs, (0.0_dp, i = 1, nslide)], sol, reused_factor)
                if (.not. reused_factor) then
                    ! One handle for the whole solve: the KKT pattern is fixed,
                    ! so the native backend may reuse its symbolic analysis.
                    if (.not. factored) then
                        if (.not. configure(solver)) exit
                    end if
                    call refactor_kkt_raw(solver, fb%n + nslide, fb%n + nslide, &
                        size(smat), srows, scols, smat, sst, kept, &
                        report%assembly_roundoff)
                    report%factorisations = report%factorisations + 1
                    if (kept) report%reused = report%reused + 1
                    if (.not. status_ok(sst)) then
                        factored = .false.
                        mu = 10*mu; cycle
                    end if
                    factored = .true.
                    call sparse_solve(solver, [rhs, (0.0_dp, i = 1, nslide)], sol, sst)
                    if (.not. status_ok(sst)) then
                        status = IDEAL_FB_SOLVE_FAILED; return
                    end if
                end if
                du = sol(:fb%n)
                report%factor_s = report%factor_s + wall() - t0
                call map_steps(fb, du, disp, rot)
                alpha = min(1.0_dp, cap_x/max(disp, tiny(1.0_dp)), &
                    cap_w/max(rot, tiny(1.0_dp)))
                if (talk) print '(a,i3,4es10.2)', '   try mu disp rot alpha', tries, &
                    mu, disp, rot, alpha
                if (talk .and. tries == 1) call direction_report(fb, du)
                ! Backtrack along this direction first: -H^-1 g descends the
                ! weighted residual merit for any SPD weight, whatever the
                ! inertia of H. Raising mu instead would turn the step against
                ! the vertical mode (negative curvature with fixed currents)
                ! once mu |H_ii| exceeds its eigenvalue.
                ! With n /= 0 modes, first damp only the n = 0 part of the step
                ! (its near-null floor components spoil the merit while the
                ! n /= 0 part is nearly linear), then the whole step.
                nback = 8
                if (fb%nmode > 1 .and. alpha >= 1) nback = 14
                do back = 1, nback
                    t0 = wall()
                    if (nback == 14 .and. back <= 6) then
                        trial_du = du
                        call scale_mode(fb, trial_du, 1, 0.5_dp**(back - 1))
                        call ideal_fb_step(fb, trial_du, 1.0_dp)
                    else
                        call ideal_fb_step(fb, du, alpha)
                    end if
                    call ideal_fb_gradient(fb, rhs, st, min_jacobian=minjac)
                    if (cready .and. st == IDEAL_FB_OK) call project_multipliers(fb, &
                        csolver, rhs)
                    if (st == IDEAL_FB_OK) then
                        ftry = rigid_force(rhs)
                        call project_slide(fb, wts, rhs)
                    end if
                    report%line_s = report%line_s + wall() - t0
                    if (st == IDEAL_FB_OK) then
                        tmerit = sum(wts*rhs**2)
                        if (tmerit < (1 - 1.0e-4_dp*alpha)*merit) then
                            accepted = .true.; exit
                        end if
                    end if
                    fb%plasma%x = x0; fb%xs = xs0; fb%phi = phi0; fb%lam = lam0
                    if (.not. (nback == 14 .and. back <= 6)) alpha = alpha/2
                end do
                if (accepted) exit
                mu = 10*mu
            end do
            if (.not. accepted) then
                report%floor = report%residual < 1.0e-6_dp
                if (report%floor .and. res_tol <= 0) then
                    status = IDEAL_FB_OK
                end if
                exit
            end if
            grad = rhs; fzr = ftry
            if (kz < 0 .and. cready) call rigid_vertical_correct()
            report%iterations = it + 1
            report%step = alpha*disp
            report%min_jacobian = minjac
            report%damping = mu
            if (.not. alpha < 1) then
                mu = max(mu/100, 1.0e-14_dp)
                ! This is step regularisation, not part of the equations.
                ! A fixed tangential penalty stalls near-null n /= 0 modes
                ! on fine meshes. Recover Newton as full steps succeed,
                ! retaining the slide gauge and the residual acceptance test.
                tangent_damping = tangent_damping/100
            end if
        end do
        if (factored) call sparse_destroy(solver)
        if (cready) call sparse_destroy(csolver)
        call constraint_norm(fb, report%constraint)
        if (talk) print '(a,4f10.1,a,i0,a,i0)', '   seconds: assemble factor+solve &
            &line-search', report%assemble_s, report%factor_s, report%line_s, &
            report%gradient_s, '  factorisations ', report%factorisations, &
            ' analysis kept ', report%reused
    contains
        real(dp) function wall()
            integer(i64) :: count, rate
            call system_clock(count, rate)
            wall = real(count, dp)/rate
        end function

        real(dp) function rigid_force(g)
            ! Net Z force on the plasma nodes (before the gauge projection).
            real(dp), intent(in) :: g(:)
            integer :: a
            rigid_force = 0
            do a = 1, fb%plasma%nodes
                rigid_force = rigid_force + g(pidx(fb, 2, 1, a))
            end do
        end function

        subroutine rigid_vertical_correct()
            ! dZ = -F_Z/K_Z on the rigid n = 0 Z mode, kept only when the
            ! weighted (slide-projected) merit decreases.
            real(dp), allocatable :: g2(:)
            real(dp) :: dz, m2
            integer :: s3
            dz = -fzr/kz
            ftry = fzr
            x0 = fb%plasma%x; xs0 = fb%xs; phi0 = fb%phi; lam0 = fb%lam
            fb%plasma%x(2, 1, :) = fb%plasma%x(2, 1, :) + dz
            call ideal_fb_solve_shell_map(fb, s3)
            if (s3 == IDEAL_FB_OK) call ideal_fb_solve_vacuum(fb, s3)
            if (s3 == IDEAL_FB_OK) call ideal_fb_solve_multipliers(fb, s3)
            allocate (g2(fb%n))
            if (s3 == IDEAL_FB_OK) call ideal_fb_gradient(fb, g2, s3, &
                min_jacobian=minjac)
            if (s3 == IDEAL_FB_OK) then
                call project_multipliers(fb, csolver, g2)
                ftry = rigid_force(g2)
                call project_slide(fb, wts, g2)
                m2 = sum(wts*g2**2)
                if (m2 < sum(wts*grad**2)) then
                    if (talk) print '(a,es10.2,a,es10.3,a,es10.3,a,es10.3,a,es10.3)', &
                        '   rigid dZ ', dz, ' merit ', sum(wts*grad**2), ' ->', m2, &
                        ' F_Z ', fzr, ' ->', ftry
                    grad = g2; fzr = ftry
                    return
                end if
            end if
            fb%plasma%x = x0; fb%xs = xs0; fb%phi = phi0; fb%lam = lam0
            if (talk) print '(a,es10.2,a,es10.3,a,es10.3)', '   rigid dZ ', dz, &
                ' rejected F_Z ', fzr, ' ->', ftry
        end subroutine
    end subroutine

    subroutine boundary_report(fb)
        ! Mean (n = 0) displacement of Gamma_p from its reference.
        type(ideal_fb_t), intent(in) :: fb
        real(dp) :: mean(2), rms
        integer :: s, nb
        mean = 0; nb = 0; rms = 0
        do s = 1, fb%ns
            if (fb%role(s) /= 1) cycle
            mean = mean + fb%xs(1:2, 1, s) - fb%reference(:, s); nb = nb + 1
            if (fb%nmode > 1) rms = rms + sum(fb%xs(1:2, 2:fb%nmode, s)**2)
        end do
        print '(a,2es12.4)', '   Gamma_p mean shift', mean/nb
        if (fb%nmode > 1) print '(a,es12.4)', '   Gamma_p n/=0 rms (m)', sqrt(rms/nb)
    end subroutine

    subroutine direction_report(fb, du)
        ! Rigid content of a Newton direction on Gamma_p (n = 0 mode): mean
        ! (dR, dZ) and the rms remainder.
        type(ideal_fb_t), intent(in) :: fb
        real(dp), intent(in) :: du(:)
        real(dp) :: mean(2), rms
        integer :: s, nb
        mean = 0; nb = 0; rms = 0
        do s = 1, fb%ns
            if (fb%role(s) /= 1) cycle
            mean = mean + [du(ix(fb, 1, 1, s)), du(ix(fb, 2, 1, s))]; nb = nb + 1
        end do
        mean = mean/nb
        do s = 1, fb%ns
            if (fb%role(s) /= 1) cycle
            rms = rms + sum(([du(ix(fb, 1, 1, s)), du(ix(fb, 2, 1, s))] - mean)**2)
        end do
        print '(a,2es10.2,a,es10.2)', '   direction on Gamma_p: mean dR dZ', mean, &
            '  rms remainder', sqrt(rms/nb)
        block
            integer :: k, a
            k = maxloc(abs(du(1:fb%np)), dim=1)
            a = (k - 1)/(2*fb%nmode) + 1
            print '(a,es10.2,a,2f9.4,a,l2)', '   largest plasma increment', du(k), &
                ' at', fb%plasma%reference(:, a), ' boundary', fb%plasma%boundary(a)
        end block
    end subroutine

    subroutine factor_constraint(fb, rowptr, colidx, val, csolver, ready)
        ! Factor the C^T block (rows: shell map, columns: multipliers); the
        ! CSR rows of the transposed block are its CSC columns.
        type(ideal_fb_t), intent(in) :: fb
        integer, intent(in) :: rowptr(:), colidx(:)
        real(dp), intent(in) :: val(:)
        type(sparse_solver_t), intent(inout) :: csolver
        logical, intent(out) :: ready
        integer, allocatable :: rp(:), ci(:)
        real(dp), allocatable :: v(:)
        type(fortsparse_status_t) :: sst
        integer :: i, p, nnz, row0, col0
        row0 = fb%np; col0 = fb%np + fb%nx
        allocate (rp(fb%nx + 1), ci(size(val)), v(size(val)))
        nnz = 0; rp(1) = 1
        do i = 1, fb%nx
            do p = rowptr(col0 + i), rowptr(col0 + i + 1) - 1
                if (colidx(p) <= row0 .or. colidx(p) > row0 + fb%nx) cycle
                nnz = nnz + 1
                ci(nnz) = colidx(p) - row0; v(nnz) = val(p)
            end do
            rp(i + 1) = nnz + 1
        end do
        csolver%backend_id = FORTSPARSE_BACKEND_NATIVE
        if (gs_direct_available()) then
            ready = gs_direct_mpi_ready()
            if (.not. ready) return
            csolver%backend_id = FORTSPARSE_BACKEND_NATIVE
        end if
        csolver%options%symmetry = FORTSPARSE_SYM_GENERAL
        csolver%options%lu_method = FORTSPARSE_LU_SUPERNODAL
        call factor_kkt_raw(csolver, fb%nx, fb%nx, nnz, rp, ci(:nnz), &
                            v(:nnz), sst, symmetric=.false.)
        ready = status_ok(sst)
    end subroutine

    subroutine project_multipliers(fb, csolver, grad)
        ! lambda += delta with C^T delta = -g_x, then the gradient updated:
        ! the shell-map rows become zero and the Gamma_p plasma rows receive
        ! -delta (their tie coefficient is -1).
        type(ideal_fb_t), intent(inout) :: fb
        type(sparse_solver_t), intent(inout) :: csolver
        real(dp), intent(inout) :: grad(:)
        real(dp), allocatable :: rhs(:), d(:)
        type(fortsparse_status_t) :: sst
        integer :: s2, m, c, k
        rhs = -grad(fb%np + 1:fb%np + fb%nx)
        allocate (d(fb%nx))
        call sparse_solve(csolver, rhs, d, sst)
        if (.not. status_ok(sst)) return
        fb%lam = fb%lam + reshape(d, shape(fb%lam))
        grad(fb%np + 1:fb%np + fb%nx) = 0
        do s2 = 1, fb%ns
            if (fb%role(s2) /= 1) cycle
            do m = 1, fb%nmode
                do c = 1, 2
                    k = c + 3*(m - 1) + 3*fb%nmode*(s2 - 1)
                    grad(pidx(fb, c, m, fb%plasma_node(s2))) = &
                        grad(pidx(fb, c, m, fb%plasma_node(s2))) - d(k)
                end do
            end do
        end do
    end subroutine

    subroutine scale_mode(fb, du, m, factor)
        ! Scale the components of real Fourier mode m of a global increment.
        type(ideal_fb_t), intent(in) :: fb
        real(dp), intent(inout) :: du(:)
        integer, intent(in) :: m
        real(dp), intent(in) :: factor
        integer :: a, s2, c
        do a = 1, fb%plasma%nodes
            do c = 1, 2
                du(pidx(fb, c, m, a)) = factor*du(pidx(fb, c, m, a))
            end do
        end do
        do s2 = 1, fb%ns
            do c = 1, 3
                du(ix(fb, c, m, s2)) = factor*du(ix(fb, c, m, s2))
                du(il(fb, c, m, s2)) = factor*du(il(fb, c, m, s2))
            end do
            du(iphi(fb, m, s2)) = factor*du(iphi(fb, m, s2))
        end do
    end subroutine

    subroutine tangential_positions(fb, rowptr, colidx, tpos)
        ! CSR positions of the (q1, q2) x (q1, q2) block of every plasma node
        ! and mode: tpos(1:4, k) = (11, 12, 21, 22), k = (mode, node).
        type(ideal_fb_t), intent(in) :: fb
        integer, intent(in) :: rowptr(:), colidx(:)
        integer, allocatable, intent(inout) :: tpos(:, :)
        integer :: a, m, k, c, r, p, col
        if (.not. allocated(tpos)) allocate (tpos(4, fb%nmode*fb%plasma%nodes))
        do a = 1, fb%plasma%nodes
            do m = 1, fb%nmode
                k = m + fb%nmode*(a - 1)
                do c = 1, 2
                    r = pidx(fb, c, m, a)
                    do p = rowptr(r), rowptr(r + 1) - 1
                        col = colidx(p)
                        if (col == pidx(fb, 1, m, a)) tpos(2*c - 1, k) = p
                        if (col == pidx(fb, 2, m, a)) tpos(2*c, k) = p
                    end do
                end do
            end do
        end do
    end subroutine

    subroutine add_tangential(fb, diag, tpos, val, damping)
        ! val += tangential_damping s_a t t^T, t the unit tangent of the
        ! reference surface at node a, s_a the mean (q1, q2) diagonal.
        type(ideal_fb_t), intent(in) :: fb
        real(dp), intent(in) :: diag(:)
        integer, intent(in) :: tpos(:, :)
        real(dp), intent(inout) :: val(:)
        real(dp), intent(in) :: damping
        real(dp) :: t(2), sa
        integer :: a, m, k
        do a = 1, fb%plasma%nodes
            t = [-fb%plasma%direction(2, a), fb%plasma%direction(1, a)]
            do m = 1, fb%nmode
                k = m + fb%nmode*(a - 1)
                sa = damping*0.5_dp*(diag(pidx(fb, 1, m, a)) &
                                                   + diag(pidx(fb, 2, m, a)))
                val(tpos(1, k)) = val(tpos(1, k)) + sa*t(1)*t(1)
                val(tpos(2, k)) = val(tpos(2, k)) + sa*t(1)*t(2)
                val(tpos(3, k)) = val(tpos(3, k)) + sa*t(2)*t(1)
                val(tpos(4, k)) = val(tpos(4, k)) + sa*t(2)*t(2)
            end do
        end do
    end subroutine

    subroutine preconditioned_solve(solver, col_ptr, row_idx, val, b, x, ok)
        ! x solves A x = b, A given in the compressed layout the sparse
        ! factor reads (column pointers, row indices), by FGMRES right-
        ! preconditioned with the factor of an earlier matrix of the same
        ! pattern; ok only when the relative residual reaches krylov_tol.
        type(sparse_solver_t), intent(inout) :: solver
        integer, intent(in) :: col_ptr(:), row_idx(:)
        real(dp), intent(in) :: val(:), b(:)
        real(dp), intent(out) :: x(:)
        logical, intent(out) :: ok
        real(dp), parameter :: krylov_tol = 1.0e-10_dp
        integer, parameter :: krylov_max = 20
        real(dp) :: bnorm, rnorm
        integer :: info, iterations
        ok = .false.
        x = 0
        bnorm = norm2(b)
        if (.not. bnorm > 0) return
        call real_fgmres_operator(apply, precondition, b/bnorm, x, krylov_tol, &
            krylov_max, krylov_max, info, iterations, rnorm)
        ok = info == KRYLOV_OK .and. all(ieee_is_finite(x))
        x = bnorm*x
    contains
        subroutine apply(input, output)
            real(dp), intent(in) :: input(:)
            real(dp), intent(out) :: output(:)
            integer :: j, p
            output = 0
            do j = 1, size(input)
                do p = col_ptr(j), col_ptr(j + 1) - 1
                    output(row_idx(p)) = output(row_idx(p)) + val(p)*input(j)
                end do
            end do
        end subroutine

        subroutine precondition(input, output)
            real(dp), intent(in) :: input(:)
            real(dp), intent(out) :: output(:)
            type(fortsparse_status_t) :: pst
            call sparse_solve(solver, input, output, pst)
            if (.not. status_ok(pst)) output = ieee_value(1.0_dp, ieee_quiet_nan)
        end subroutine
    end subroutine

    subroutine prepare_kkt_values(solver, n, ptr, idx, values, correction, &
                                  status, symmetric)
        ! CSR is also CSC after this checked, explicit symmetry completion.
        ! The underlying scalar Lagrangian defines the symmetric matrix;
        ! independent accumulation may differ within the roundoff budget.
        type(sparse_solver_t), intent(inout) :: solver
        integer, intent(in) :: n, ptr(:), idx(:)
        real(dp), intent(inout) :: values(:)
        real(dp), intent(out) :: correction
        type(fortsparse_status_t), intent(out) :: status
        logical, optional, intent(in) :: symmetric
        correction = 0
        status = fortsparse_status_t()
        if (present(symmetric)) then
            if (.not. symmetric) return
        end if
        call csc_normalize_symmetric_raw(n, ptr, idx, values, 1.0e-12_dp, &
            correction, status)
        if (.not. status_ok(status)) solver%factored = .false.
    end subroutine

    subroutine factor_kkt_raw(solver, nr, nc, nnz, ptr, idx, values, status, &
                              symmetric)
        type(sparse_solver_t), intent(inout) :: solver
        integer, intent(in) :: nr, nc, nnz, ptr(:), idx(:)
        real(dp), intent(inout) :: values(:)
        type(fortsparse_status_t), intent(out) :: status
        logical, optional, intent(in) :: symmetric
        real(dp) :: correction
        call prepare_kkt_values(solver, nc, ptr, idx, values, correction, &
                                status, symmetric)
        if (.not. status_ok(status)) return
        call sparse_factor(solver, nr, nc, nnz, ptr, idx, values, status)
    end subroutine

    subroutine refactor_kkt_raw(solver, nr, nc, nnz, ptr, idx, values, &
                                status, reused, maximum_correction)
        type(sparse_solver_t), intent(inout) :: solver
        integer, intent(in) :: nr, nc, nnz, ptr(:), idx(:)
        real(dp), intent(inout) :: values(:)
        type(fortsparse_status_t), intent(out) :: status
        logical, intent(out) :: reused
        real(dp), optional, intent(inout) :: maximum_correction
        real(dp) :: correction
        reused = .false.
        call prepare_kkt_values(solver, nc, ptr, idx, values, correction, status)
        if (.not. status_ok(status)) return
        if (present(maximum_correction)) &
            maximum_correction = max(maximum_correction, correction)
        call sparse_refactor(solver, nr, nc, nnz, ptr, idx, values, status, reused)
    end subroutine

    logical function configure(solver) result(ready)
        ! Native pivoted LU for the KKT system. Symmetry completion remains
        ! an assembly contract, independent of this factorization choice.
        ! Native LDL^T has not met the runtime budget for this system.
        type(sparse_solver_t), intent(inout) :: solver
        if (gs_direct_available()) then
            ready = gs_direct_mpi_ready()
            if (.not. ready) return
            solver%backend_id = FORTSPARSE_BACKEND_NATIVE
            solver%options%symmetry = FORTSPARSE_SYM_GENERAL
            solver%options%lu_method = FORTSPARSE_LU_SUPERNODAL
            solver%options%ordering = FORTSPARSE_ORDER_AMD
            solver%options%nthreads = 0
        else
            solver%backend_id = FORTSPARSE_BACKEND_NATIVE
            solver%options%symmetry = FORTSPARSE_SYM_GENERAL
            solver%options%lu_method = FORTSPARSE_LU_SUPERNODAL
        end if
        ready = .true.
    end function

    subroutine map_steps(fb, du, disp, rot)
        ! Largest displacement (m) and omega (rad) increment of a direction.
        type(ideal_fb_t), intent(in) :: fb
        real(dp), intent(in) :: du(:)
        real(dp), intent(out) :: disp, rot
        integer :: s, m
        disp = maxval(abs(du(1:fb%np)))
        rot = 0
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

    subroutine ideal_fb_solve_shell_map(fb, status)
        ! Shell map from the constraints at fixed plasma map: Gamma_p copies,
        ! harmonic extension, Gamma_c fixed.
        type(ideal_fb_t), intent(inout) :: fb
        integer, intent(out) :: status
        real(dp), allocatable :: du(:)
        call block_solve(fb, fb%np + fb%nx, fb%np, fb%nx, du, status)
        if (status == IDEAL_FB_OK) fb%xs = fb%xs + reshape(du, shape(fb%xs))
    end subroutine

    subroutine ideal_fb_rigid_stiffness(fb, c, delta, force, stiffness, status)
        ! Generalised force sum_a dL/dq_c(a) (n = 0) on a rigid plasma shift in
        ! R (c = 1) or Z (c = 2), with shell map, Phi~ and lambda re-solved at
        ! each shifted geometry, and its central-difference stiffness
        ! dforce/dshift (J/m^2 per field period). stiffness < 0 in Z is the
        ! vertical instability of the fixed-current saddle; the state is
        ! restored on return.
        type(ideal_fb_t), intent(inout) :: fb
        integer, intent(in) :: c
        real(dp), intent(in) :: delta
        real(dp), intent(out) :: force, stiffness
        integer, intent(out) :: status
        real(dp), allocatable :: x0(:, :, :), xs0(:, :, :), phi0(:, :), lam0(:, :, :)
        real(dp) :: f(-1:1)
        integer :: k
        x0 = fb%plasma%x; xs0 = fb%xs; phi0 = fb%phi; lam0 = fb%lam
        do k = -1, 1
            fb%plasma%x(c, 1, :) = x0(c, 1, :) + k*delta
            call settle(f(k))
            if (status /= IDEAL_FB_OK) exit
        end do
        fb%plasma%x = x0; fb%xs = xs0; fb%phi = phi0; fb%lam = lam0
        force = f(0)
        stiffness = (f(1) - f(-1))/(2*delta)
    contains
        subroutine settle(value)
            real(dp), intent(out) :: value
            real(dp), allocatable :: g(:)
            integer :: a, st
            call ideal_fb_solve_shell_map(fb, status)
            if (status /= IDEAL_FB_OK) return
            call ideal_fb_solve_vacuum(fb, status)
            if (status /= IDEAL_FB_OK) return
            call ideal_fb_solve_multipliers(fb, status)
            if (status /= IDEAL_FB_OK) return
            allocate (g(fb%n))
            call ideal_fb_gradient(fb, g, st)
            value = 0
            do a = 1, fb%plasma%nodes
                value = value + g(pidx(fb, c, 1, a))
            end do
        end subroutine
    end subroutine

    subroutine ideal_fb_newton_step(fb, status, nonaxisymmetric_only)
        ! One undamped Newton step on the slide-gauged KKT system (no line
        ! search): [H, s G^T; s G, 0] [du; mu] = [-g; 0].
        type(ideal_fb_t), intent(inout) :: fb
        integer, intent(out) :: status
        logical, intent(in), optional :: nonaxisymmetric_only
        integer, allocatable :: rowptr(:), colidx(:), rows(:), cols(:)
        real(dp), allocatable :: val(:), grad(:), du(:), matrix(:), sol(:)
        type(sparse_solver_t) :: solver
        type(fortsparse_status_t) :: sst
        integer :: st, nslide
        status = IDEAL_FB_SOLVE_FAILED
        nslide = size(fb%slide, 3)
        allocate (grad(fb%n), du(fb%n), sol(fb%n + nslide))
        if (.not. ideal_fb_exterior_ready(fb)) return
        call ideal_fb_gradient(fb, grad, st)
        call ideal_fb_hessian(fb, rowptr, colidx, val)
        call ideal_fb_slide_system(fb, rowptr, colidx, val, plasma_diagonal_mean(fb, &
            rowptr, colidx, val), rows, cols, matrix)
        if (.not. configure(solver)) return
        call factor_kkt_raw(solver, fb%n + nslide, fb%n + nslide, size(matrix), &
                           rows, cols, matrix, sst)
        if (.not. status_ok(sst)) return
        call sparse_solve(solver, [-grad, (0.0_dp, st = 1, nslide)], sol, sst)
        call sparse_destroy(solver)
        if (.not. status_ok(sst)) return
        du = sol(:fb%n)
        if (present(nonaxisymmetric_only)) then
            ! At an axisymmetric state the modes decouple: keep the n /= 0
            ! part (the linear response to the n /= 0 drive) only.
            if (nonaxisymmetric_only) call zero_mode(1)
        end if
        call ideal_fb_step(fb, du, 1.0_dp)
        status = IDEAL_FB_OK
    contains
        subroutine zero_mode(m)
            integer, intent(in) :: m
            integer :: a, s2, c
            do a = 1, fb%plasma%nodes
                do c = 1, 2
                    du(pidx(fb, c, m, a)) = 0
                end do
            end do
            do s2 = 1, fb%ns
                do c = 1, 3
                    du(ix(fb, c, m, s2)) = 0; du(il(fb, c, m, s2)) = 0
                end do
                du(iphi(fb, m, s2)) = 0
            end do
        end subroutine
    end subroutine

    real(dp) function plasma_diagonal_mean(fb, rowptr, colidx, val) result(mean)
        ! Mean |H_ii| over the plasma rows: scale of the slide-gauge border.
        type(ideal_fb_t), intent(in) :: fb
        integer, intent(in) :: rowptr(:), colidx(:)
        real(dp), intent(in) :: val(:)
        integer :: i, p
        mean = 0
        do i = 1, fb%np
            do p = rowptr(i), rowptr(i + 1) - 1
                if (colidx(p) == i) mean = mean + abs(val(p))
            end do
        end do
        mean = mean/fb%np
    end function

    subroutine mode_matrix(fb, m, rp, ci, v, mass, map, np_m, status)
        ! An axisymmetric state decouples by harmonic, NOT by real phase:
        ! B_pol.grad and B_phi.d/dphi couple cosine and sine coefficients.
        type(ideal_fb_t), intent(in) :: fb
        integer, intent(in) :: m
        integer, allocatable, intent(out) :: rp(:), ci(:), map(:)
        real(dp), allocatable, intent(out) :: v(:), mass(:)
        integer, intent(out) :: np_m, status
        integer, allocatable :: rowptr(:), colidx(:), inverse(:)
        real(dp), allocatable :: val(:)
        real(dp) :: correction
        type(fortsparse_status_t) :: sst
        integer :: first, last, a, c, phase, i, p, n, nnz, s
        status = IDEAL_FB_BAD_INPUT
        np_m = 0
        if (allocated(fb%exterior_3d)) return
        if (.not. allocated(fb%s_trace)) return
        if (m < 1 .or. m > fb%nmode) return
        if (fb%nmode > 1) then
            if (maxval(abs(fb%plasma%x(:, 2:, :))) > 1.0e-10_dp) return
            if (maxval(abs(fb%xs(:, 2:, :))) > 1.0e-10_dp) return
            if (maxval(abs(fb%phi(2:, :))) > 1.0e-10_dp) return
        end if
        first = 1; last = 1
        if (m > 1) then
            first = 2*(m/2); last = first + 1
        end if
        call ideal_fb_hessian(fb, rowptr, colidx, val)
        ! Apply the same checked symmetry completion as the Newton solve,
        ! before sub-block extraction changes the matrix's scale.
        call csc_normalize_symmetric_raw(fb%n, rowptr, colidx, val, &
            1.0e-12_dp, correction, sst)
        if (.not. status_ok(sst)) then
            status = IDEAL_FB_SOLVE_FAILED
            return
        end if
        allocate (map(fb%n))
        map = 0; n = 0
        do a = 1, fb%plasma%nodes
            do phase = first, last
                do c = 1, 2
                    n = n + 1; map(pidx(fb, c, phase, a)) = n
                end do
            end do
        end do
        np_m = n
        do s = 1, fb%ns
            do phase = first, last
                do c = 1, 3
                    n = n + 1; map(ix(fb, c, phase, s)) = n
                    n = n + 1; map(il(fb, c, phase, s)) = n
                end do
                n = n + 1; map(iphi(fb, phase, s)) = n
            end do
        end do
        allocate (rp(n + 1), ci(size(val)), v(size(val)), mass(np_m), inverse(n))
        do i = 1, fb%n
            if (map(i) > 0) inverse(map(i)) = i
        end do
        nnz = 0; rp(1) = 1; mass = 0
        do i = 1, n
            do p = rowptr(inverse(i)), rowptr(inverse(i) + 1) - 1
                if (map(colidx(p)) == 0) cycle
                nnz = nnz + 1
                ci(nnz) = map(colidx(p)); v(nnz) = val(p)
                if (ci(nnz) == i .and. i <= np_m) mass(i) = abs(val(p))
            end do
            rp(i + 1) = nnz + 1
            call sort_mode_row(ci(rp(i):nnz), v(rp(i):nnz))
        end do
        if (any(mass <= 0)) return
        ci = ci(:nnz); v = v(:nnz)
        status = IDEAL_FB_OK
    end subroutine

    subroutine sort_mode_row(col, value)
        integer, intent(inout) :: col(:)
        real(dp), intent(inout) :: value(:)
        integer :: j, k, key
        real(dp) :: tv
        do j = 2, size(col)
            key = col(j); tv = value(j); k = j - 1
            do while (k >= 1)
                if (col(k) <= key) exit
                col(k + 1) = col(k); value(k + 1) = value(k); k = k - 1
            end do
            col(k + 1) = key; value(k + 1) = tv
        end do
    end subroutine

    subroutine mode_ledger(fb, rp, ci, v, map, np_m, z, mass, report, action)
        type(ideal_fb_t), intent(in) :: fb
        integer, intent(in) :: rp(:), ci(:), map(:), np_m
        real(dp), intent(in) :: v(:), z(:), mass(:)
        type(ideal_fb_stability_t), intent(out) :: report
        real(dp), intent(out) :: action(:)
        real(dp) :: u(fb%n), cm, trace(fb%nt)
        integer :: i, p, j, m, a
        report = ideal_fb_stability_t()
        action = 0; u = 0
        do i = 1, fb%n
            if (map(i) > 0) u(i) = z(map(i))
        end do
        do i = 1, size(rp) - 1
            do p = rp(i), rp(i + 1) - 1
                action(i) = action(i) + v(p)*z(ci(p))
                if (i <= np_m .and. ci(p) <= np_m) &
                    report%plasma = report%plasma + 0.5_dp*z(i)*v(p)*z(ci(p))
            end do
        end do
        report%total = 0.5_dp*dot_product(z, action)
        do i = fb%np + fb%nx + 1, fb%np + 2*fb%nx
            if (map(i) > 0) report%constraint = report%constraint + u(i)*action(map(i))
        end do
        do m = 1, fb%nmode
            cm = ideal_fb_mode_weight(fb, m)/fb%mu0
            do a = 1, fb%nt
                trace(a) = u(iphi(fb, m, fb%trace(a)))
            end do
            do j = 1, fb%nt
                report%exterior = report%exterior &
                    - 0.5_dp*cm*trace(j) &
                    *dot_product(trace, fb%s_trace(:, j, harmonic(m)))
            end do
        end do
        report%vacuum = report%total - report%plasma - report%constraint
        report%metric_norm = sum(mass*z(:np_m)**2)
    end subroutine

    subroutine ideal_fb_mode_energy(fb, m, displacement, report, status)
        ! Schur energy of a prescribed plasma displacement. For n>0 the
        ! layout is (Rcos,Zcos,Rsin,Zsin) per node; n=0 is (R,Z) per node.
        ! Shell map, potential and multipliers are solved at fixed plasma.
        type(ideal_fb_t), intent(in) :: fb
        integer, intent(in) :: m
        real(dp), intent(in) :: displacement(:)
        type(ideal_fb_stability_t), intent(out) :: report
        integer, intent(out) :: status
        integer, allocatable :: rp(:), ci(:), map(:), ar(:), ac(:)
        real(dp), allocatable :: v(:), mass(:), av(:), rhs(:), z(:), action(:)
        type(sparse_solver_t) :: solver
        type(fortsparse_status_t) :: sst
        integer :: np_m, n, i, p, count
        report = ideal_fb_stability_t()
        call mode_matrix(fb, m, rp, ci, v, mass, map, np_m, status)
        if (status /= IDEAL_FB_OK) return
        status = IDEAL_FB_BAD_INPUT
        if (size(displacement) /= np_m) return
        if (.not. all(ieee_is_finite(displacement))) return
        n = size(rp) - 1
        allocate (ar(n - np_m + 1), ac(size(v)), av(size(v)), &
            rhs(n - np_m), z(n), action(n))
        z(:np_m) = displacement; rhs = 0; count = 0; ar(1) = 1
        do i = np_m + 1, n
            do p = rp(i), rp(i + 1) - 1
                if (ci(p) <= np_m) then
                    rhs(i - np_m) = rhs(i - np_m) - v(p)*z(ci(p))
                else
                    count = count + 1
                    ac(count) = ci(p) - np_m; av(count) = v(p)
                end if
            end do
            ar(i - np_m + 1) = count + 1
        end do
        status = IDEAL_FB_SOLVE_FAILED
        if (.not. configure(solver)) return
        call factor_kkt_raw(solver, n - np_m, n - np_m, count, ar, ac(:count), &
            av(:count), sst)
        if (.not. status_ok(sst)) return
        call sparse_solve(solver, rhs, z(np_m + 1:), sst)
        call sparse_destroy(solver)
        if (.not. status_ok(sst)) return
        call mode_ledger(fb, rp, ci, v, map, np_m, z, mass, report, action)
        status = IDEAL_FB_OK
    end subroutine

    subroutine mode_rest_block(rp, ci, v, np, ar, ac, av)
        ! Extract the eliminated shell-map/potential/multiplier block C.
        integer, intent(in) :: rp(:), ci(:), np
        real(dp), intent(in) :: v(:)
        integer, allocatable, intent(out) :: ar(:), ac(:)
        real(dp), allocatable, intent(out) :: av(:)
        integer :: i, p, count, n
        n = size(rp) - 1
        allocate (ar(n - np + 1), ac(size(v)), av(size(v)))
        count = 0; ar(1) = 1
        do i = np + 1, n
            do p = rp(i), rp(i + 1) - 1
                if (ci(p) <= np) cycle
                count = count + 1
                ac(count) = ci(p) - np; av(count) = v(p)
            end do
            ar(i - np + 1) = count + 1
        end do
        ac = ac(:count); av = av(:count)
    end subroutine

    subroutine ideal_fb_mode_inertia(fb, m, shifts, counts, status)
        ! Haynsworth: In(S - sigma M) = In(K - diag(sigma M,0)) - In(C).
        ! No dense Schur matrix. Counts are (negative,zero,positive) computed
        ! pivot signs, NOT rigorous enclosures of the exact-real operator.
        ! C must be nonsingular; uncertain pivot signs reject the diagnostic.
        type(ideal_fb_t), intent(in) :: fb
        integer, intent(in) :: m
        real(dp), intent(in) :: shifts(:)
        integer, intent(out) :: counts(3, size(shifts)), status
        integer, allocatable :: rp(:), ci(:), map(:), ar(:), ac(:)
        real(dp), allocatable :: v(:), mass(:), av(:), shifted(:)
        integer :: np, n, j, i, p, base(3), full(3), st
        counts = -1
        call mode_matrix(fb, m, rp, ci, v, mass, map, np, status)
        if (status /= IDEAL_FB_OK) return
        status = IDEAL_FB_BAD_INPUT
        if (.not. all(ieee_is_finite(shifts))) return
        n = size(rp) - 1
        call mode_rest_block(rp, ci, v, np, ar, ac, av)
        call count_factor(ar, ac, av, base, st)
        status = IDEAL_FB_SOLVE_FAILED
        if (st /= IDEAL_FB_OK .or. base(2) /= 0) return
        allocate (shifted(size(v)))
        do j = 1, size(shifts)
            shifted = v
            do i = 1, np
                do p = rp(i), rp(i + 1) - 1
                    if (ci(p) == i) shifted(p) = v(p) - shifts(j)*mass(i)
                end do
            end do
            call count_factor(rp, ci, shifted, full, st)
            if (st /= IDEAL_FB_OK) return
            full = full - base
            if (any(full < 0) .or. sum(full) /= np) return
            counts(:, j) = full
        end do
        status = IDEAL_FB_OK
    contains
        subroutine count_factor(rows, cols, values, inertia, result)
            integer, intent(in) :: rows(:), cols(:)
            real(dp), intent(in) :: values(:)
            integer, intent(out) :: inertia(3), result
            type(sparse_solver_t) :: solver
            type(fortsparse_status_t) :: sst
            real(dp), allocatable :: scale(:), scaled(:)
            integer :: nr, i, p, reliability
            integer(int64) :: signs(3)
            result = IDEAL_FB_SOLVE_FAILED; inertia = -1
            nr = size(rows) - 1
            allocate (scale(nr), scaled(size(values)))
            ! Positive diagonal congruence controls mixed physical units;
            ! it cannot change inertia. A zero row is a rejected singularity.
            do i = 1, nr
                scale(i) = maxval(abs(values(rows(i):rows(i + 1) - 1)))
            end do
            if (any(scale <= 0)) return
            scale = 1/sqrt(scale)
            !$omp parallel do default(shared) private(i,p) schedule(static)
            do i = 1, nr
                do p = rows(i), rows(i + 1) - 1
                    scaled(p) = values(p)*(scale(i)*scale(cols(p)))
                end do
            end do
            !$omp end parallel do
            solver%backend_id = FORTSPARSE_BACKEND_NATIVE
            solver%options%symmetry = FORTSPARSE_SYM_INDEFINITE
            solver%options%ordering = FORTSPARSE_ORDER_AMD
            solver%options%ldlt_blocked = .true.
            call sparse_factor(solver, nr, nr, size(values), rows, cols, scaled, sst)
            if (status_ok(sst)) then
                call sparse_inertia(solver, signs(1), signs(2), signs(3), &
                    reliability, sst)
                if (status_ok(sst)) then
                    if (reliability == FORTSPARSE_INERTIA_RELIABLE) then
                        inertia = int(signs)
                        result = IDEAL_FB_OK
                    end if
                end if
            end if
            if (.not. status_ok(sst)) print '(a,a)', 'mode inertia: ', sst%msg
            call sparse_destroy(solver)
        end subroutine
    end subroutine

    subroutine ideal_fb_mode_lanczos(fb, m, lambda, report, status, shift)
        ! Opt-in near-shift Schur eigenpair. FortNum #103 checks the original
        ! residual with an observed LOWER operator-norm estimate (stricter
        ! than the exact norm) and exact ||M||_2 = max(mass). No unit floor.
        ! This gate does not establish extremal ordering or physical accuracy.
        type(ideal_fb_t), intent(in) :: fb
        integer, intent(in) :: m
        real(dp), intent(out) :: lambda
        type(ideal_fb_stability_t), intent(out) :: report
        integer, intent(out) :: status
        real(dp), intent(in), optional :: shift
        integer, allocatable :: rp(:), ci(:), map(:), ar(:), ac(:)
        real(dp), allocatable :: v(:), mass(:), av(:), shifted(:), start(:)
        real(dp), allocatable :: rhs(:), z(:), crhs(:), cz(:), images(:)
        real(dp), allocatable :: vectors(:, :)
        real(dp) :: values(1), residuals(1), sigma, relative_residual, observed_scale
        real(dp) :: radius, backward
        type(sparse_solver_t) :: full_solver, rest_solver
        type(fortsparse_status_t) :: sst
        type(eigen_options_t) :: options
        type(eigen_stats_t) :: stats
        integer :: np, n, i, p, info
        logical :: action_failed
        lambda = ieee_value(0.0_dp, ieee_quiet_nan)
        report = ideal_fb_stability_t()
        call mode_matrix(fb, m, rp, ci, v, mass, map, np, status)
        if (status /= IDEAL_FB_OK) return
        sigma = 0
        if (present(shift)) sigma = shift
        status = IDEAL_FB_BAD_INPUT
        if (.not. ieee_is_finite(sigma)) return
        n = size(rp) - 1
        call mode_rest_block(rp, ci, v, np, ar, ac, av)
        shifted = v
        do i = 1, np
            do p = rp(i), rp(i + 1) - 1
                if (ci(p) == i) shifted(p) = v(p) - sigma*mass(i)
            end do
        end do
        status = IDEAL_FB_SOLVE_FAILED
        if (.not. configure(full_solver)) return
        if (.not. configure(rest_solver)) return
        call factor_kkt_raw(rest_solver, n - np, n - np, size(av), ar, ac, av, sst)
        if (.not. status_ok(sst)) then
            call sparse_destroy(rest_solver)
            return
        end if
        call factor_kkt_raw(full_solver, n, n, size(v), rp, ci, shifted, sst)
        if (.not. status_ok(sst)) then
            call sparse_destroy(rest_solver)
            call sparse_destroy(full_solver)
            return
        end if
        allocate (start(np), vectors(np, 1), rhs(n), z(n), &
            crhs(n - np), cz(n - np), images(n))
        do i = 1, np
            start(i) = 1 + 0.1_dp*sin(1.7_dp*i)
        end do
        options%max_basis = min(32, np)
        options%retain = min(16, np - 1)
        options%max_restarts = 25
        options%which = EIGEN_NEAREST
        options%mass_scale = maxval(mass)
        action_failed = .false.; observed_scale = 0
        call lanczos(apply_schur, start, values, vectors, residuals, &
            1.0e-10_dp, options, info, stats, mass=apply_mass, &
            solve_shift=solve_shift, sigma=sigma)
        status = IDEAL_FB_NOT_CONVERGED
        if (info == EIGEN_OK .and. .not. action_failed) then
            ! Reconstruct the SAME physical Schur vector, then its ledger.
            call apply_schur(vectors(:, 1), start)
            if (.not. action_failed) then
                z(:np) = vectors(:, 1); z(np + 1:) = -cz
                lambda = values(1)
                radius = norm2((start - lambda*mass*z(:np))/sqrt(mass))
                relative_residual = radius/max(abs(lambda), 1.0e-12_dp)
                backward = norm2(start - lambda*mass*z(:np)) &
                    /((observed_scale + abs(lambda)*maxval(mass))*norm2(z(:np)))
                call mode_ledger(fb, rp, ci, v, map, np, z, mass, report, images)
                report%relative_residual = relative_residual
                report%residual_radius = radius
                report%backward_error = backward
                status = IDEAL_FB_OK
            end if
        end if
        if (action_failed) status = IDEAL_FB_SOLVE_FAILED
        call sparse_destroy(full_solver)
        call sparse_destroy(rest_solver)
    contains
        subroutine apply_mass(x, y)
            real(dp), intent(in) :: x(:)
            real(dp), intent(out) :: y(:)
            y = mass*x
        end subroutine
        subroutine solve_shift(x, y, info)
            real(dp), intent(in) :: x(:)
            real(dp), intent(out) :: y(:)
            integer, intent(out) :: info
            type(fortsparse_status_t) :: st
            rhs = 0; rhs(:np) = x
            call sparse_solve(full_solver, rhs, z, st)
            info = 1
            if (status_ok(st)) info = 0
            y = z(:np)
        end subroutine
        subroutine apply_schur(x, y)
            real(dp), intent(in) :: x(:)
            real(dp), intent(out) :: y(:)
            type(fortsparse_status_t) :: st
            integer :: i, p
            crhs = 0
            !$omp parallel do default(shared) private(i,p) schedule(static)
            do i = np + 1, n
                do p = rp(i), rp(i + 1) - 1
                    if (ci(p) <= np) crhs(i - np) = &
                        crhs(i - np) + v(p)*x(ci(p))
                end do
            end do
            !$omp end parallel do
            call sparse_solve(rest_solver, crhs, cz, st)
            if (.not. status_ok(st)) action_failed = .true.
            y = ieee_value(0.0_dp, ieee_quiet_nan)
            if (action_failed) return
            y = 0
            !$omp parallel do default(shared) private(i,p) schedule(static)
            do i = 1, np
                do p = rp(i), rp(i + 1) - 1
                    if (ci(p) <= np) then
                        y(i) = y(i) + v(p)*x(ci(p))
                    else
                        y(i) = y(i) - v(p)*cz(ci(p) - np)
                    end if
                end do
            end do
            !$omp end parallel do
            if (norm2(x) > 0) observed_scale = max(observed_scale, norm2(y)/norm2(x))
        end subroutine
    end subroutine

    subroutine ideal_fb_mode_stability(fb, m, lambda, status, iterations, shift, report)
        ! Shifted inverse iteration on the full real harmonic pair, with
        ! shell map, multipliers and potential eliminated by the KKT solve.
        ! A negative converged Rayleigh quotient witnesses instability of
        ! this discrete gamma=0 model. A positive quotient is NOT an inertia
        ! count or proof of stability: only the eigenvalue near shift is sought.
        type(ideal_fb_t), intent(in) :: fb
        integer, intent(in) :: m
        real(dp), intent(out) :: lambda
        integer, intent(out) :: status
        integer, intent(in), optional :: iterations
        real(dp), intent(in), optional :: shift
        type(ideal_fb_stability_t), intent(out), optional :: report
        integer, allocatable :: rp(:), ci(:), map(:)
        real(dp), allocatable :: v(:), original(:), mass(:), x(:), z(:), rhs(:), hz(:)
        type(sparse_solver_t) :: solver
        type(fortsparse_status_t) :: sst
        type(ideal_fb_stability_t) :: ledger
        integer :: np_m, n, i, p, it, nit, nnz
        real(dp) :: sigma, norm
        lambda = ieee_value(0.0_dp, ieee_quiet_nan)
        if (present(report)) report = ideal_fb_stability_t()
        call mode_matrix(fb, m, rp, ci, v, mass, map, np_m, status)
        if (status /= IDEAL_FB_OK) return
        sigma = 0
        if (present(shift)) sigma = shift
        ! Fine meshes can need more than 80 steps near a clustered spectrum.
        ! This is a work budget; the unshifted residual gate stays unchanged.
        nit = 800
        if (present(iterations)) nit = iterations
        status = IDEAL_FB_BAD_INPUT
        if (nit < 1 .or. .not. ieee_is_finite(sigma)) return
        n = size(rp) - 1; nnz = rp(n + 1) - 1
        original = v
        do i = 1, np_m
            do p = rp(i), rp(i + 1) - 1
                if (ci(p) == i) v(p) = v(p) - sigma*mass(i)
            end do
        end do
        status = IDEAL_FB_SOLVE_FAILED
        if (.not. configure(solver)) return
        call factor_kkt_raw(solver, n, n, nnz, rp, ci(:nnz), v(:nnz), sst)
        if (.not. status_ok(sst)) return
        allocate (x(np_m), z(n), rhs(n), hz(n))
        do i = 1, np_m
            x(i) = 1 + 0.1_dp*sin(1.7_dp*i)
        end do
        x = x/sqrt(sum(mass*x**2))
        status = IDEAL_FB_NOT_CONVERGED
        do it = 1, nit
            rhs = 0; rhs(1:np_m) = mass*x
            call sparse_solve(solver, rhs, z, sst)
            if (.not. status_ok(sst)) then
                status = IDEAL_FB_SOLVE_FAILED
                exit
            end if
            norm = sqrt(sum(mass*z(:np_m)**2))
            if (.not. ieee_is_finite(norm) .or. norm <= 0) exit
            z = z/norm
            ! The unshifted action gives an actual Rayleigh quotient and a
            ! metric-scaled eigen-residual, even before iteration converges.
            !$omp parallel do default(shared) private(i, p) schedule(static)
            do i = 1, n
                hz(i) = 0
                do p = rp(i), rp(i + 1) - 1
                    hz(i) = hz(i) + original(p)*z(ci(p))
                end do
            end do
            !$omp end parallel do
            lambda = dot_product(z, hz)
            ledger%relative_residual = norm2((hz(:np_m) - lambda*mass*z(:np_m)) &
                /sqrt(mass))/max(abs(lambda), 1.0e-12_dp)
            x = z(:np_m)
            if (ledger%relative_residual < 1.0e-6_dp) then
                status = IDEAL_FB_OK
                exit
            end if
        end do
        call sparse_destroy(solver)
        if (present(report)) then
            norm = ledger%relative_residual
            call mode_ledger(fb, rp, ci, original, map, np_m, z, mass, ledger, hz)
            ledger%relative_residual = norm
            report = ledger
        end if
    end subroutine

    subroutine ideal_fb_hypotheses(fb, h, status, dense_limit)
        ! Checks of the open discrete hypotheses at the current state
        ! (docs/free-boundary-theory.md FB-A2, FB-C5, FB-D2, FB-E1; derivation
        ! theory/derivations/free-boundary-gauge-vertical.md section 5):
        !   (V) A_h + S_h > 0 (H_PhiPhi < 0): dense eigenvalues of every S_k
        !       and of the assembled -mu0 H_PhiPhi;
        !   (R) full row rank of R = [grad C; G] on the map unknowns:
        !       lambda_min(R R^T) by inverse iteration through the augmented
        !       factor [I R^T; R 0] (nonsingular iff R has full row rank; the
        !       Rayleigh value is an upper bound), sigma_max by power iteration;
        !   S_red: the Schur complement S_q of the exact KKT matrix onto the
        !       plasma positions q (shell map, Phi~ and lambda eliminated).
        !       grad C fixes the shell map as xs = E q, so S_q = [I; E]^T S_map
        !       [I; E], and S_red = Z_G^T S_q Z_G with Z_G a QR basis of ker G.
        !       Only the Gamma_p tie rows couple q to the other unknowns, so
        !       S_q = H_qq - K_qw K_ww^-1 K_wq needs one solve per coupled
        !       plasma unknown. Dense eigenvalues give its inertia;
        !   FB-C5: the slide multiplier of the merit projection.
        ! Under (V) and (R), FB-D2 (Haynsworth inertia additivity and the
        ! bordering Lemma B2) gives In(K_full) = In(S_red) + (r, r + dim Phi~,
        ! 0) for the slide-bordered KKT matrix; this is reported as kkt_inertia
        ! and checked against a dense eigen-solve in test_ideal_free_boundary.
        ! FortSparse (pinned) exposes no MUMPS inertia, hence the dense path:
        ! nf and np must not exceed dense_limit (default 8000).
        type(ideal_fb_t), intent(in) :: fb
        type(ideal_fb_hypotheses_t), intent(out) :: h
        integer, intent(out) :: status
        integer, intent(in), optional :: dense_limit
        interface
            subroutine dgeqrf(m, n, a, lda, tau, work, lwork, info)
                import :: dp
                integer, intent(in) :: m, n, lda, lwork
                real(dp), intent(inout) :: a(lda, *)
                real(dp), intent(out) :: tau(*), work(*)
                integer, intent(out) :: info
            end subroutine
            subroutine dormqr(side, trans, m, n, k, a, lda, tau, c, ldc, work, &
                              lwork, info)
                import :: dp
                character, intent(in) :: side, trans
                integer, intent(in) :: m, n, k, lda, ldc, lwork
                real(dp), intent(in) :: a(lda, *), tau(*)
                real(dp), intent(inout) :: c(ldc, *)
                real(dp), intent(out) :: work(*)
                integer, intent(out) :: info
            end subroutine
        end interface
        integer, allocatable :: rowptr(:), colidx(:), ci(:), cj(:), arow(:), acol(:)
        integer, allocatable :: coupled(:)
        real(dp), allocatable :: val(:), a(:, :), w(:), cv(:), aval(:), y(:), z(:)
        real(dp), allocatable :: rhs(:), sol(:), grad(:), diag(:), wts(:), gs(:)
        real(dp), allocatable :: mult(:), gt(:, :), tau(:), work(:), x(:), v(:)
        integer, allocatable :: ri(:), rj(:)
        real(dp), allocatable :: rv(:)
        type(sparse_solver_t) :: solver
        type(fortsparse_status_t) :: sst
        integer :: limit, i, j, p, k, q, it, info, nl, i0, nmap, nr, nred, nw, nc, cnt
        integer :: c, aa, st
        real(dp) :: est, num, den, tol, query(1)
        logical :: ok
        limit = 8000
        if (present(dense_limit)) limit = dense_limit
        status = IDEAL_FB_BAD_INPUT
        h%status = status
        if (allocated(fb%exterior_3d)) return
        if (.not. allocated(fb%s_trace)) return
        if (fb%nf > limit .or. fb%np > limit) return
        status = IDEAL_FB_SOLVE_FAILED
        nl = size(fb%slide, 3)
        call ideal_fb_hessian(fb, rowptr, colidx, val)

        ! (V): S_k and -mu0 H_PhiPhi.
        h%w_asymmetry = maxval(fb%bem_asymmetry(1, :))
        h%v_asymmetry = maxval(fb%bem_asymmetry(2, :))
        h%s_min = huge(1.0_dp)
        do k = 0, fb%ntor
            a = fb%s_trace(:, :, k)
            call symmetric_eigen(a, w, .false., info)
            if (info /= 0) return
            h%s_min = min(h%s_min, w(1))
        end do
        i0 = fb%np + 2*fb%nx
        deallocate (a)
        allocate (a(fb%nf, fb%nf))
        a = 0
        do i = 1, fb%nf
            do p = rowptr(i0 + i), rowptr(i0 + i + 1) - 1
                if (colidx(p) > i0) a(i, colidx(p) - i0) = -fb%mu0*val(p)
            end do
        end do
        num = 0; den = 0
        do j = 1, fb%nf
            do i = 1, fb%nf
                num = num + (a(i, j) - a(j, i))**2; den = den + a(i, j)**2
            end do
        end do
        h%vacuum_asymmetry = sqrt(num/den)
        call symmetric_eigen(a, w, .false., info)
        if (info /= 0) return
        h%vacuum_min = w(1); h%vacuum_max = w(fb%nf)
        deallocate (a)

        ! (R): R = [grad C; G] on the map unknowns (plasma, shell map).
        nmap = fb%np + fb%nx; nr = fb%nx + nl
        h%constraint_rows = nr
        cnt = count(abs(fb%slide) > 0)
        do i = fb%np + fb%nx + 1, fb%np + 2*fb%nx
            cnt = cnt + count(colidx(rowptr(i):rowptr(i + 1) - 1) <= nmap)
        end do
        allocate (ri(cnt), rj(cnt), rv(cnt))
        cnt = 0
        do i = fb%np + fb%nx + 1, fb%np + 2*fb%nx
            do p = rowptr(i), rowptr(i + 1) - 1
                if (colidx(p) > nmap) cycle
                cnt = cnt + 1
                ri(cnt) = i - fb%np - fb%nx; rj(cnt) = colidx(p); rv(cnt) = val(p)
            end do
        end do
        do k = 1, nl
            do aa = 1, fb%plasma%nodes
                do c = 1, 2
                    if (.not. abs(fb%slide(c, aa, k)) > 0) cycle
                    cnt = cnt + 1
                    ri(cnt) = fb%nx + k; rj(cnt) = pidx(fb, c, 1, aa)
                    rv(cnt) = fb%slide(c, aa, k)
                end do
            end do
        end do
        ! Augmented [I R^T; R 0]: its (2, 2) inverse block is -(R R^T)^-1.
        allocate (ci(nmap + 2*cnt), cj(nmap + 2*cnt), cv(nmap + 2*cnt))
        ci(:nmap) = [(i, i=1, nmap)]; cj(:nmap) = ci(:nmap); cv(:nmap) = 1
        ci(nmap + 1:nmap + cnt) = nmap + ri; cj(nmap + 1:nmap + cnt) = rj
        ci(nmap + cnt + 1:) = rj; cj(nmap + cnt + 1:) = nmap + ri
        cv(nmap + 1:nmap + cnt) = rv; cv(nmap + cnt + 1:) = rv
        call coo_to_csr(nmap + nr, ci, cj, cv, arow, acol, aval)
        deallocate (ci, cj, cv)
        allocate (y(nr), z(nr), rhs(nmap + nr), sol(nmap + nr), x(nmap))
        ok = configure(solver)
        if (ok) then
            call sparse_factor(solver, nmap + nr, nmap + nr, size(aval), arow, acol, &
                               aval, sst)
            ok = status_ok(sst)
        end if
        h%constraint_sigma_min = 0
        est = 0
        if (ok) then
            y = [(1 + 0.1_dp*sin(1.3_dp*i), i=1, nr)]
            y = y/norm2(y)
            do it = 1, 40
                rhs = 0; rhs(nmap + 1:) = y
                call sparse_solve(solver, rhs, sol, sst)
                if (.not. status_ok(sst)) exit
                z = -sol(nmap + 1:)
                est = dot_product(y, z)
                y = z/norm2(z)
            end do
            if (status_ok(sst) .and. est > 0) h%constraint_sigma_min = sqrt(1/est)
            call sparse_destroy(solver)
        end if
        y = [(1 + 0.1_dp*cos(0.7_dp*i), i=1, nr)]
        y = y/norm2(y)
        do it = 1, 40
            x = 0
            do p = 1, cnt
                x(rj(p)) = x(rj(p)) + rv(p)*y(ri(p))
            end do
            z = 0
            do p = 1, cnt
                z(ri(p)) = z(ri(p)) + rv(p)*x(rj(p))
            end do
            est = norm2(z)
            y = z/est
        end do
        h%constraint_sigma_max = sqrt(est)
        h%full_rank = h%constraint_sigma_min > h%rank_tolerance*h%constraint_sigma_max
        deallocate (arow, acol, aval, y, z, rhs, sol, x)

        ! S_q: K_ww (all unknowns but the plasma positions), factorised once.
        nw = fb%n - fb%np
        cnt = 0
        do i = fb%np + 1, fb%n
            cnt = cnt + count(colidx(rowptr(i):rowptr(i + 1) - 1) > fb%np)
        end do
        allocate (ci(cnt), cj(cnt), cv(cnt))
        cnt = 0
        do i = fb%np + 1, fb%n
            do p = rowptr(i), rowptr(i + 1) - 1
                if (colidx(p) <= fb%np) cycle
                cnt = cnt + 1
                ci(cnt) = i - fb%np; cj(cnt) = colidx(p) - fb%np; cv(cnt) = val(p)
            end do
        end do
        call coo_to_csr(nw, ci, cj, cv, arow, acol, aval)
        deallocate (ci, cj, cv)
        if (.not. configure(solver)) return
        call sparse_factor(solver, nw, nw, size(aval), arow, acol, aval, sst)
        if (.not. status_ok(sst)) return
        allocate (a(fb%np, fb%np), coupled(fb%np), rhs(nw), sol(nw))
        a = 0; nc = 0
        do i = 1, fb%np
            do p = rowptr(i), rowptr(i + 1) - 1
                if (colidx(p) <= fb%np) then
                    a(i, colidx(p)) = val(p)
                else if (nc == 0) then
                    nc = nc + 1; coupled(nc) = i
                else if (coupled(nc) /= i) then
                    nc = nc + 1; coupled(nc) = i
                end if
            end do
        end do
        do q = 1, nc
            j = coupled(q)
            rhs = 0
            do p = rowptr(j), rowptr(j + 1) - 1
                if (colidx(p) > fb%np) rhs(colidx(p) - fb%np) = val(p)
            end do
            call sparse_solve(solver, rhs, sol, sst)
            if (.not. status_ok(sst)) then
                call sparse_destroy(solver); return
            end if
            do k = 1, nc
                i = coupled(k)
                do p = rowptr(i), rowptr(i + 1) - 1
                    if (colidx(p) > fb%np) a(i, j) = a(i, j) &
                        - val(p)*sol(colidx(p) - fb%np)
                end do
            end do
        end do
        call sparse_destroy(solver)
        do j = 1, fb%np
            do i = j + 1, fb%np
                a(i, j) = 0.5_dp*(a(i, j) + a(j, i)); a(j, i) = a(i, j)
            end do
        end do
        ! Z_G: the trailing columns of Q in G^T = Q R.
        allocate (gt(fb%np, nl), tau(nl), work(1))
        gt = 0
        do k = 1, nl
            do aa = 1, fb%plasma%nodes
                do c = 1, 2
                    gt(pidx(fb, c, 1, aa), k) = fb%slide(c, aa, k)
                end do
            end do
        end do
        call dgeqrf(fb%np, nl, gt, fb%np, tau, query, -1, info)
        if (info /= 0) return
        deallocate (work); allocate (work(max(1, nint(query(1)))))
        call dgeqrf(fb%np, nl, gt, fb%np, tau, work, size(work), info)
        if (info /= 0) return
        deallocate (work); allocate (work(max(1, 64*fb%np)))
        call dormqr('L', 'T', fb%np, fb%np, nl, gt, fb%np, tau, a, fb%np, work, &
                    size(work), info)
        if (info /= 0) return
        call dormqr('R', 'N', fb%np, fb%np, nl, gt, fb%np, tau, a, fb%np, work, &
                    size(work), info)
        if (info /= 0) return
        nred = fb%np - nl
        h%reduced_dim = nred
        block
            real(dp), allocatable :: b(:, :)
            b = a(nl + 1:, nl + 1:)
            deallocate (a)
            call symmetric_eigen(b, w, .true., info)
            if (info /= 0) return
            tol = 1.0e-10_dp*maxval(abs(w))
            h%zero_tolerance = tol
            h%n_positive = count(w > tol); h%n_negative = count(w < -tol)
            h%n_zero = nred - h%n_positive - h%n_negative
            h%reduced_min = w(1); h%reduced_max = w(nred)
            h%reduced_min_abs = minval(abs(w))
            allocate (x(fb%np), v(fb%np))
            x = 0; x(nl + 1:) = b(:, 1)
            call dormqr('L', 'N', fb%np, 1, nl, gt, fb%np, tau, x, fb%np, work, &
                        size(work), info)
            if (info /= 0) return
            v = 0
            do aa = 1, fb%plasma%nodes
                v(pidx(fb, 2, 1, aa)) = 1
            end do
            h%vertical_overlap = abs(dot_product(x, v))/(norm2(x)*norm2(v))
        end block
        h%kkt_inertia = [h%n_positive + nr, h%n_negative + nr + fb%nf, h%n_zero]

        ! FB-C5: slide multiplier of the merit projection.
        allocate (grad(fb%n), diag(fb%n), wts(fb%n), mult(nl))
        call ideal_fb_gradient(fb, grad, st)
        if (st /= IDEAL_FB_OK) return
        diag = 0
        do i = 1, fb%n
            do p = rowptr(i), rowptr(i + 1) - 1
                if (colidx(p) == i) diag(i) = abs(val(p))
            end do
        end do
        call merit_weights(fb, diag, wts)
        gs = grad
        call project_slide(fb, wts, gs, mult)
        h%slide_multiplier = maxval(abs(mult))
        h%slide_residual = sqrt(sum(wts*(grad - gs)**2)/fb%energy_scale)
        status = IDEAL_FB_OK
        h%status = status
    end subroutine

    subroutine ideal_fb_hypotheses_print(h, iteration)
        ! One line of the per-state hypothesis report.
        type(ideal_fb_hypotheses_t), intent(in) :: h
        integer, intent(in), optional :: iteration
        integer :: it
        it = -1
        if (present(iteration)) it = iteration
        print '(a,i4,a,i0,a,2es9.2,a,es10.3,a,es10.3,a,l1,a,3(1x,i0),a,es11.3, &
            &a,es10.3,a,es10.3)', '   hypotheses', it, ' status ', h%status, &
            ' W,V asym ', h%w_asymmetry, h%v_asymmetry, ' min eig(A+S) ', &
            h%vacuum_min, ' sigma_min(R) ', h%constraint_sigma_min, ' full rank ', &
            h%full_rank, ' In(S_red)', h%n_positive, h%n_negative, h%n_zero, &
            ' min ', h%reduced_min, ' min|.| ', h%reduced_min_abs, ' slide mu ', &
            h%slide_multiplier
    end subroutine

    subroutine symmetric_eigen(a, w, vectors, info)
        ! Eigenvalues (ascending) and optionally eigenvectors (in a) of a dense
        ! symmetric matrix, LAPACK dsyev with a workspace query.
        real(dp), intent(inout) :: a(:, :)
        real(dp), allocatable, intent(out) :: w(:)
        logical, intent(in) :: vectors
        integer, intent(out) :: info
        interface
            subroutine dsyev(jobz, uplo, n, a, lda, w, work, lwork, info)
                import :: dp
                character, intent(in) :: jobz, uplo
                integer, intent(in) :: n, lda, lwork
                real(dp), intent(inout) :: a(lda, *)
                real(dp), intent(out) :: w(*), work(*)
                integer, intent(out) :: info
            end subroutine
        end interface
        real(dp) :: query(1)
        real(dp), allocatable :: work(:)
        character :: job
        integer :: n
        n = size(a, 1)
        job = merge('V', 'N', vectors)
        allocate (w(n))
        call dsyev(job, 'U', n, a, n, w, query, -1, info)
        if (info /= 0) return
        allocate (work(max(1, nint(query(1)))))
        call dsyev(job, 'U', n, a, n, w, work, size(work), info)
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
        if (.not. ideal_fb_exterior_ready(fb)) return
        call ideal_fb_gradient(fb, grad, st)
        call ideal_fb_hessian(fb, rowptr, colidx, val)
        ! fortsparse takes CSC: the CSR rows of the transposed block (rows
        ! col0.., columns row0.. of the symmetric KKT matrix) are the CSC
        ! columns of the wanted block.
        allocate (rp(nb + 1), ci(size(val)), v(size(val)))
        nnz = 0; rp(1) = 1
        do i = 1, nb
            do p = rowptr(col0 + i), rowptr(col0 + i + 1) - 1
                if (colidx(p) <= row0 .or. colidx(p) > row0 + nb) cycle
                nnz = nnz + 1
                ci(nnz) = colidx(p) - row0; v(nnz) = val(p)
            end do
            rp(i + 1) = nnz + 1
        end do
        if (row0 == col0) then
            if (.not. configure(solver)) return
        else
            solver%backend_id = FORTSPARSE_BACKEND_NATIVE
            if (gs_direct_available()) then
                if (.not. gs_direct_mpi_ready()) return
                solver%backend_id = FORTSPARSE_BACKEND_NATIVE
            end if
            solver%options%symmetry = FORTSPARSE_SYM_GENERAL
            solver%options%lu_method = FORTSPARSE_LU_SUPERNODAL
        end if
        call factor_kkt_raw(solver, nb, nb, nnz, rp, ci(:nnz), v(:nnz), &
                            sst, symmetric=(row0 == col0))
        if (.not. status_ok(sst)) return
        rhs = -grad(row0 + 1:row0 + nb)
        call sparse_solve(solver, rhs, du, sst)
        call sparse_destroy(solver)
        if (.not. status_ok(sst)) return
        status = IDEAL_FB_OK
    end subroutine

    subroutine ideal_fb_shell_field(fb, e, xi, phi_hat, x, b, status)
        ! Physical point x = (R, Z, phi) and cylindrical field (B_R, B_phi,
        ! B_Z) of the vacuum potential at master point xi of shell cell e,
        ! for any chart (cylindrical: x = q exactly, as before).
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
        select type (chart => fb%plasma%chart)
        type is (cylindrical_chart_t)
            b(1) = bc(1)*cos(q(3)) + bc(2)*sin(q(3))
            b(2) = -bc(1)*sin(q(3)) + bc(2)*cos(q(3))
            b(3) = bc(3)
            x = q([1, 2, 3])
        class default
            ! Any other chart: cylindrical coordinates of the Cartesian point,
            ! the angle unwrapped to the branch nearest the chart angle q3.
            x(3) = atan2(cart(2), cart(1))
            x(3) = x(3) + 2*pi*nint((q(3) - x(3))/(2*pi))
            x(1) = hypot(cart(1), cart(2)); x(2) = cart(3)
            b(1) = bc(1)*cos(x(3)) + bc(2)*sin(x(3))
            b(2) = -bc(1)*sin(x(3)) + bc(2)*cos(x(3))
            b(3) = bc(3)
        end select
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
