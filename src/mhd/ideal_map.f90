module kin6d_ideal_map
    ! Nonlinear 3D ideal equilibrium as frozen-in energy minimisation over maps
    ! (DESIGN 7.1, choice A). The unknown is a map from abstract coordinates
    ! s = (s1, s2, s3): (s1, s2) the master coordinates of each triangle of a
    ! 2D reference mesh and s3 = phi-hat, into Cartesian space,
    !     X(s) = H(q(s)),  q = (q1, q2, phi-hat + omega),
    ! where H is a chart and q1, q2, omega are Lagrange P_k in
    ! the plane times a real Fourier series in phi-hat with n in nfp*Z, |n| <=
    ! nfp*ntor. The frozen flux density b (a closed 2-form, contravariant density
    ! in s) and pressure p (gamma = 0) are carried by the reference; the field is
    ! B = DX b/J with J = sigma*det DX, and
    !     W[X] = int ( |DX b|^2/(2 mu0 J) - p J ) ds   (one field period).
    ! The energy density, its DX-gradient and DX-Hessian are generated leaves;
    ! the energy, quadrature and solver see the chart only through H(q) and its
    ! q-derivatives, so a curve-following frame is a further ideal_chart_t.
    ! Relabelling gauge (DESIGN 7.1): each node moves in the poloidal plane only
    ! along a fixed transversal d_a = grad chi/|grad chi| of the reference label,
    !     (q1, q2)(a) = y_a + nu_a(phi-hat) d_a,
    ! so the unknowns are (nu, omega) per node and mode (a VMEC-like fixed-angle
    ! gauge; omega plays the role of lambda). Surfaces still move freely; only
    ! the tangential reparametrisation along them is fixed. A fixed direction
    ! cannot translate surfaces smaller than the axis shift, and any transversal
    ! field is singular at the O-point (a cone-type map that loses order). A
    ! global axis expansion (VMEC's raxis/zaxis generalised) is therefore added:
    !     (q1, q2)(a) = y_a + nu_a d_a + T_a sum_j sigma_j P_j(y_a),
    ! with T = (chi - chi_edge)/(chi_axis - chi_edge) in [0, 1] and P_j the
    ! monomials of degree <= k-1 about the O-point node, so the remaining
    ! transversal part is smooth enough for the P_k interpolation order.
    ! Fixed boundary: nu on boundary nodes is held (all modes); omega's n = 0
    ! value at one node fixes the rigid toroidal rotation. At n = 0 any
    ! omega = f(chi) is a continuum null direction (rigid rotation of each
    ! surface); discretely it is near-null and Newton wanders along it. A gauge
    ! term fixes it: (gauge_weight/2) int R^2 (det G grad chi . grad omega_0)^2
    ! / det G over the reference, the perpendicular counterpart of the physical
    ! R^2 (b_pol . grad omega)^2 / J. The physics never sees the surface
    ! function it selects (continuum); the discrete bias is measured by the
    ! tests at two weights.
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use kin6d_triangle_basis, only: triangle_basis_jet
    use kin6d_curved_triangle, only: triangle_locate
    use kin6d_grad_shafranov, only: gs_state_t
    use kin6d_ideal_map_leaves, only: ideal_point_gradient_jets, &
        ideal_point_hessian_jets
    use kin6d_equilibrium_chart, only: ideal_chart_t, cylindrical_chart_t, curve_chart_t
    use kin6d_curve_chart_leaves, only: curve_omega_metric
    use fortfem_gauss_quadrature_2d, only: gauss_quadrature_triangle_t
    use fortfem_triangle_duffy_quadrature, only: triangle_duffy_quadrature
    use fortfem_element_to_csr, only: csr_matrix_t, csr_build_pattern, &
        csr_assemble_blocks
    use fortnum_krylov, only: real_conjugate_gradient_operator, KRYLOV_OK
    use fortsparse, only: sparse_solver_t, sparse_factor, sparse_solve, &
        sparse_destroy, fortsparse_status_t, status_ok, &
        FORTSPARSE_BACKEND_NATIVE, csc_normalize_symmetric_raw, &
        FORTSPARSE_ORDER_AMD, &
        FORTSPARSE_SYM_GENERAL, FORTSPARSE_LU_SUPERNODAL, &
        sparse_factor_stats, sparse_factor_stats_t, &
        sparse_refactor, FORTSPARSE_SYM_SPD, FORTSPARSE_LU_SCALAR, &
        FORTSPARSE_NOT_POSITIVE_DEFINITE, csc_t
    use kin6d_gs_direct, only: gs_direct_available, gs_direct_mpi_ready
    implicit none
    private
    public :: ideal_map_t, ideal_frozen_t
    public :: ideal_map_init, ideal_map_freeze, ideal_map_energy, ideal_map_hessian
    public :: ideal_map_solve, ideal_map_point, ideal_map_dof, ideal_map_mode_norms
    public :: ideal_map_reduce, ideal_map_expand, ideal_map_rdof, ideal_map_axis_dof
    public :: ideal_map_prolong, ideal_map_residual, ideal_map_position_gauge
    public :: ideal_map_solve_report_t, ideal_map_eliminate_omega
    public :: ideal_map_gauge_project
    ! PF1 profiling instrumentation (islands/pf1prof, islands/sp2meas only, not
    ! for main). pf1_t/pf1_n: 1 energy, 2 energy+gradient, 3 Hessian cell loop
    ! (wall), 4 CSR pattern+assembly, 5 leaf (thread sum), 6 rank updates
    ! (thread sum), 7 P^T H P (thread sum), 8 factor_masked cache update and
    ! diagonal, 9 sparse_factor, 10 sparse_refactor, 11 sparse_solve,
    ! 12 refactor reported reused.
    real(dp), public :: pf1_t(12) = 0
    integer, public :: pf1_n(12) = 0
    type(sparse_factor_stats_t), public :: pf1_stats
    ! Damping above which factor_masked no longer tries the SPD Cholesky.
    real(dp), parameter :: spd_mu_limit = 1.0e2_dp
    integer, parameter, public :: IDEAL_MAP_OK = 0, IDEAL_MAP_FOLDED = 1, &
        IDEAL_MAP_NOT_CONVERGED = 2, IDEAL_MAP_SOLVE_FAILED = 3, &
        IDEAL_MAP_BAD_INPUT = 4

    ! Frozen reference data at a reference point y = (q1, q2): poloidal flux
    ! label chi and its gradient, toroidal flux density b3 (contravariant density
    ! in (q1, q2, phi-hat)) and pressure. b_pol = (-d chi/dq2, d chi/dq1).
    type, abstract :: ideal_frozen_t
    contains
        procedure(frozen_evaluate), deferred :: evaluate
    end type
    abstract interface
        subroutine frozen_evaluate(field, point, chi, chi_gradient, b3, pressure)
            import :: ideal_frozen_t, dp
            class(ideal_frozen_t), intent(in) :: field
            real(dp), intent(in) :: point(2)
            real(dp), intent(out) :: chi, chi_gradient(2), b3, pressure
        end subroutine
    end interface

    type :: ideal_map_t
        class(ideal_chart_t), allocatable :: chart
        integer :: nfp = 1, ntor = 0, nmode = 1, nphi = 1
        integer :: nodes = 0, cells = 0, nloc = 0, nq = 0, ndof = 0, gauge = 0
        integer :: nred = 0
        real(dp) :: mu0 = 1, sigma = 1
        integer, allocatable :: conn(:, :), dofconn(:, :), redconn(:, :)
        ! Transversal directions d(2, node); fixed(nred) flags reduced dofs.
        ! axis_basis(j, node) = T P_j at the node; npoly monomials per mode.
        ! basis(:, j, node): the two reduced increment directions in
        ! (dq1, dq2, domega). Transversal gauge: (d, 0) and (0, 0, 1).
        real(dp), allocatable :: direction(:, :), axis_basis(:, :), basis(:, :, :)
        ! Nested solve: vector axis basis axis_vec(:, p, node) of axis dof
        ! p = i + 2 (j - 1) (tangential fields, see nested_axis_basis); when
        ! allocated it replaces e_i axis_basis(j, node). axis_removed(1, :):
        ! tangential fraction^2 of the dropped near-normal directions;
        ! axis_removed(2, :): cos^2 of the dropped directions with the
        ! field-line slides.
        real(dp), allocatable :: axis_vec(:, :, :), axis_removed(:, :)
        ! Axis expansion: polynomial degree (-1: k - 1) and power of T (set
        ! before ideal_map_init / ideal_map_freeze respectively).
        integer :: axis_degree = -1, axis_power = 1
        ! Native factor thread budget. Newton follows near-null (relabelling)
        ! directions, so roundoff differences between factorisations move the
        ! discrete solution along them (measured: 3.9e-7 against 5.4e-7 on the
        ! finest A3 P3 seed mesh). One thread keeps the factorisation, and with
        ! the serial cell sums the whole solve, reproducible; 0 uses all
        ! threads (faster numeric factorisation, not reproducible).
        ! 0 follows the caller's OpenMP budget; positive values cap the factor.
        integer :: factor_threads = 0
        ! Debug: verify every sparse factorisation and cross-check the cached
        ! symmetric factor input (off in production, see factor_masked).
        logical :: verify_factor = .false.
        ! Position gauge instead of transversals + axis expansion + omega.
        logical :: position_gauge = .false.
        real(dp) :: initial_damping = 1.0e-6_dp
        ! Non-nested Newton step by preconditioned CG on the assembled damped
        ! Hessian with the last factorisation as preconditioner (inexact
        ! Newton, Eisenstat-Walker forcing). The factorisation is redone when
        ! CG needs more than krylov_refactor iterations or does not converge
        ! within krylov_max. False (default): factor and solve at every trial.
        ! Off by default because the P3 seed ladder (stop at 1e-11) crawls
        ! along near-null directions that inexact steps do not resolve.
        logical :: krylov_step = .false.
        integer :: krylov_refactor = 8, krylov_max = 20
        integer :: npoly = 0
        ! Gauge quadratic form on the n = 0 omega, per cell (nloc, nloc, cell).
        real(dp) :: gauge_weight = 1.0e-4_dp
        ! True: the gauge form regularises the Hessian only (energy and
        ! gradient are the physical ones, so the end point is the stationary
        ! point of W and independent of gauge_weight); False: it is part of
        ! the energy and biases the end point.
        logical :: gauge_hessian_only = .false.
        ! True: nested solve. Newton on the transversal dofs of the reduced
        ! functional W_red(nu) = min_omega W(nu, omega); omega is eliminated
        ! exactly after every step (W is quadratic in omega at fixed nu since J
        ! does not depend on omega at ntor = 0; ntor > 0 is rejected) and
        ! the outer step is the Schur-complement Newton step, damped on the
        ! transversal block only. No gauge form: only the exact null (the
        ! omega pin) is fixed. Set before ideal_map_freeze.
        logical :: omega_elimination = .false.
        ! Nested solve: surface-label gauge of omega (nested_omega_gauge).
        ! omega_layers bumps in rho = (1 - T)^1/2 (0: none; -1 before freeze:
        ! one per element layer); omega_c(l, node): moment weights, omega_cgi:
        ! inverse Gram matrix of their free-omega columns.
        integer :: omega_layers = -1
        real(dp), allocatable :: omega_c(:, :), omega_cgi(:, :)
        real(dp), allocatable :: omega_gauge(:, :, :)
        ! Nodal label, axis-expansion weight T and the axis node.
        real(dp), allocatable :: chi_node(:), axis_t(:)
        integer :: axis_node = 0
        real(dp), allocatable :: reference(:, :)
        logical, allocatable :: boundary(:), fixed(:)
        real(dp), allocatable :: weight(:), shape(:, :), dshape(:, :, :)
        real(dp), allocatable :: phi(:), trig(:, :), dtrig(:, :)
        ! Fourier-line products of the Hessian contraction over the nphi planes,
        ! pphi(j, m + nmode (m' - 1), t) = T_t(m, j) T'_t(m', j) for the pairs
        ! (trig, trig), (trig, dtrig), (dtrig, dtrig), t = 1, 2, 3.
        real(dp), allocatable :: pphi(:, :, :)
        ! Frozen master-coordinate data b(3, q, e), p(q, e), reference det(q, e).
        real(dp), allocatable :: b(:, :, :), p(:, :), det(:, :)
        ! Coefficients x(3, mode, node): q1, q2 and omega.
        real(dp), allocatable :: x(:, :, :)
    end type

    type :: factor_cache_t
        ! Per-pattern data of factor_masked. base is the masked, symmetrised
        ! matrix with the diagonal unshifted, so that a trial only rewrites the
        ! diagonal of val. keep/gauge/pattern identify the entries, raw the
        ! Hessian values that base was built from.
        integer :: nc = -1
        logical :: pattern_ok = .false., base_ok = .false.
        real(dp) :: scale = 0, correction = 0
        logical, allocatable :: keep(:)
        integer, allocatable :: in_rowptr(:), in_colidx(:)
        ! Factor input pattern, source index in matrix%val (0: gauge entry),
        ! diagonal positions, masked entries and reflections of the upper
        ! triangle onto the lower one.
        integer, allocatable :: rowptr(:), colidx(:), src(:), dpos(:), refl(:)
        integer, allocatable :: mpos(:), gpos(:)
        real(dp), allocatable :: mval(:), gval(:)
        real(dp), allocatable :: raw(:), base(:), bdiag(:)
        ! The factor input itself (a%val: base with the trial diagonal). It
        ! goes through FortSparse's csc_t entry points, which keep the symbolic
        ! analysis after a failed numeric factor, unlike the raw-array ones.
        type(csc_t) :: a
    end type

    type :: ideal_map_solve_report_t
        integer :: iterations = 0
        real(dp) :: energy = 0, gradient = huge(1.0_dp), step = huge(1.0_dp)
        real(dp) :: min_jacobian = 0, damping = 0, initial_gradient = 0
        ! Wall seconds: energy/gradient, Hessian assembly, factor + solve.
        real(dp) :: energy_s = 0, hessian_s = 0, factor_s = 0
        ! Backend split of factor_s where reported: symbolic analysis
        ! and numerical factorisation.
        real(dp) :: analysis_s = 0, numeric_s = 0
        real(dp) :: assembly_roundoff = 0
        ! Factorisation attempts, and symbolic analyses (fresh factors and
        ! refactors that report no reuse).
        integer :: analyses = 0, factorizations = 0
        ! Non-nested Newton: step trials, CG solves and CG iterations.
        integer :: trials = 0
        integer :: krylov_solves = 0, krylov_iterations = 0
        type(factor_cache_t), private :: cache
    end type
contains
    pure integer function ideal_map_dof(map, component, mode, node) result(dof)
        type(ideal_map_t), intent(in) :: map
        integer, intent(in) :: component, mode, node
        dof = component + 3*(mode - 1) + 3*map%nmode*(node - 1)
    end function

    pure integer function ideal_map_rdof(map, j, mode, node) result(dof)
        ! Reduced dof: j = 1 transversal displacement nu, j = 2 omega.
        type(ideal_map_t), intent(in) :: map
        integer, intent(in) :: j, mode, node
        dof = j + 2*(mode - 1) + 2*map%nmode*(node - 1)
    end function

    pure integer function ideal_map_axis_dof(map, i, j, mode) result(dof)
        ! Reduced dof of axis-expansion monomial j, component i (q1 or q2).
        type(ideal_map_t), intent(in) :: map
        integer, intent(in) :: i, j, mode
        dof = 2*map%nmode*map%nodes + i + 2*(j - 1) + 2*map%npoly*(mode - 1)
    end function

    subroutine ideal_map_reduce(map, g, gr)
        ! Coefficient gradient (or any dual vector) to the reduced dofs.
        type(ideal_map_t), intent(in) :: map
        real(dp), intent(in) :: g(:)
        real(dp), intent(out) :: gr(:)
        integer :: a, m, k, j, l
        do a = 1, map%nodes
            do m = 1, map%nmode
                k = ideal_map_dof(map, 1, m, a)
                gr(ideal_map_rdof(map, 1, m, a)) = dot_product(map%basis(:, 1, a), &
                    g(k:k + 2))
                gr(ideal_map_rdof(map, 2, m, a)) = dot_product(map%basis(:, 2, a), &
                    g(k:k + 2))
            end do
        end do
        do m = 1, map%nmode
            do j = 1, map%npoly
                do k = 1, 2
                    gr(ideal_map_axis_dof(map, k, j, m)) = 0
                    if (allocated(map%axis_vec)) then
                        do a = 1, map%nodes
                            l = ideal_map_dof(map, 1, m, a)
                            gr(ideal_map_axis_dof(map, k, j, m)) = &
                                gr(ideal_map_axis_dof(map, k, j, m)) &
                                + dot_product(map%axis_vec(:, k + 2*(j - 1), a), &
                                g(l:l + 1))
                        end do
                        cycle
                    end if
                    do a = 1, map%nodes
                        gr(ideal_map_axis_dof(map, k, j, m)) = &
                            gr(ideal_map_axis_dof(map, k, j, m)) &
                            + map%axis_basis(j, a)*g(ideal_map_dof(map, k, m, a))
                    end do
                end do
            end do
        end do
    end subroutine

    subroutine ideal_map_expand(map, vr, v)
        ! Reduced increment to coefficient increment.
        type(ideal_map_t), intent(in) :: map
        real(dp), intent(in) :: vr(:)
        real(dp), intent(out) :: v(:)
        integer :: a, m, k, j
        do a = 1, map%nodes
            do m = 1, map%nmode
                k = ideal_map_dof(map, 1, m, a)
                v(k:k + 2) = vr(ideal_map_rdof(map, 1, m, a))*map%basis(:, 1, a) &
                    + vr(ideal_map_rdof(map, 2, m, a))*map%basis(:, 2, a)
                if (allocated(map%axis_vec)) then
                    do j = 1, 2*map%npoly
                        v(k:k + 1) = v(k:k + 1) + map%axis_vec(:, j, a) &
                            *vr(ideal_map_axis_dof(map, 1, 1, m) + j - 1)
                    end do
                    cycle
                end if
                do j = 1, map%npoly
                    v(k:k + 1) = v(k:k + 1) + map%axis_basis(j, a) &
                        *vr(ideal_map_axis_dof(map, 1, j, m) + [0, 1])
                end do
            end do
        end do
    end subroutine

    subroutine ideal_map_init(map, mesh, chart, nfp, ntor, status, nphi, &
                              quadrature_order, mu0, axis_degree)
        ! Reference mesh, element order and boundary flags from a KIN6D GS mesh.
        ! The reference chart coordinates of the nodes are mesh%points.
        type(ideal_map_t), intent(out) :: map
        type(gs_state_t), intent(in) :: mesh
        class(ideal_chart_t), intent(in) :: chart
        integer, intent(in) :: nfp, ntor
        integer, intent(out) :: status
        integer, intent(in), optional :: nphi, quadrature_order
        real(dp), intent(in), optional :: mu0
        integer, intent(in), optional :: axis_degree
        type(gauss_quadrature_triangle_t) :: rule
        real(dp) :: h(2, 2, 10)
        real(dp), allocatable :: qxi(:), qeta(:), qw(:)
        integer :: g, j, m, k, a, c, e, l, order
        status = IDEAL_MAP_BAD_INPUT
        if (nfp < 1 .or. ntor < 0) return
        select type (chart)
        class is (curve_chart_t)
            ! Geometry couples modes even when the unknown has ntor=0.
            ! Require an explicit periodic quadrature and the admitted period.
            if (.not. present(nphi)) return
            if (nfp /= chart%field_periods()) return
        end select
        map%nloc = size(mesh%triangles, 1)
        if (map%nloc /= 6 .and. map%nloc /= 10) return
        allocate (map%chart, source=chart)
        map%nfp = nfp; map%ntor = ntor; map%nmode = 2*ntor + 1
        map%nphi = 3*ntor + 1
        if (present(nphi)) map%nphi = nphi
        if (map%nphi < 2*ntor + 1) return
        if (present(mu0)) map%mu0 = mu0
        map%nodes = size(mesh%points, 2); map%cells = size(mesh%triangles, 2)
        map%conn = mesh%triangles
        map%reference = mesh%points
        map%boundary = mesh%boundary
        order = 7
        if (present(quadrature_order)) order = quadrature_order
        ! Symmetric Gauss rules up to degree 7; collapsed (Duffy) Gauss-Legendre
        ! rules of the requested polynomial degree above that.
        if (order <= 7) then
            call rule%init(order)
            qxi = rule%xi; qeta = rule%eta; qw = rule%weights
        else
            call triangle_duffy_quadrature(order, qxi, qeta, qw, g)
            if (g /= 0) return
        end if
        map%nq = size(qw)
        allocate (map%weight(map%nq), map%shape(map%nloc, map%nq), &
            map%dshape(2, map%nloc, map%nq))
        do g = 1, map%nq
            map%weight(g) = qw(g)
            call triangle_basis_jet([qxi(g), qeta(g)], map%shape(:, g), &
                map%dshape(:, :, g), h(:, :, :map%nloc))
        end do
        allocate (map%phi(map%nphi), map%trig(map%nmode, map%nphi), &
            map%dtrig(map%nmode, map%nphi))
        do j = 1, map%nphi
            map%phi(j) = 2*acos(-1.0_dp)*(j - 1)/(nfp*map%nphi)
            map%trig(1, j) = 1; map%dtrig(1, j) = 0
            do k = 1, ntor
                map%trig(2*k, j) = cos(k*nfp*map%phi(j))
                map%trig(2*k + 1, j) = sin(k*nfp*map%phi(j))
                map%dtrig(2*k, j) = -k*nfp*sin(k*nfp*map%phi(j))
                map%dtrig(2*k + 1, j) = k*nfp*cos(k*nfp*map%phi(j))
            end do
        end do
        allocate (map%pphi(map%nphi, map%nmode**2, 3))
        do k = 1, map%nmode
            do m = 1, map%nmode
                map%pphi(:, m + map%nmode*(k - 1), 1) = map%trig(m, :)*map%trig(k, :)
                map%pphi(:, m + map%nmode*(k - 1), 2) = map%trig(m, :)*map%dtrig(k, :)
                map%pphi(:, m + map%nmode*(k - 1), 3) = &
                    map%dtrig(m, :)*map%dtrig(k, :)
            end do
        end do
        map%ndof = 3*map%nmode*map%nodes
        k = merge(2, 3, map%nloc == 6)
        if (present(axis_degree)) map%axis_degree = axis_degree
        if (map%axis_degree >= 0) k = map%axis_degree + 1
        map%npoly = k*(k + 1)/2
        map%nred = 2*map%nmode*map%nodes + 2*map%npoly*map%nmode
        allocate (map%dofconn(3*map%nmode*map%nloc, map%cells))
        allocate (map%redconn(2*map%nmode*map%nloc + 2*map%npoly*map%nmode, &
            map%cells))
        do m = 1, map%nmode
            do j = 1, map%npoly
                do c = 1, 2
                    map%redconn(2*map%nmode*map%nloc + c + 2*(j - 1) &
                        + 2*map%npoly*(m - 1), :) = ideal_map_axis_dof(map, c, j, m)
                end do
            end do
        end do
        do e = 1, map%cells
            l = 0
            do a = 1, map%nloc
                do m = 1, map%nmode
                    do c = 1, 3
                        l = l + 1
                        map%dofconn(l, e) = ideal_map_dof(map, c, m, map%conn(a, e))
                        if (c < 3) map%redconn(2*((l - 1)/3) + c, e) = &
                            ideal_map_rdof(map, c, m, map%conn(a, e))
                    end do
                end do
            end do
        end do
        allocate (map%x(3, map%nmode, map%nodes))
        map%x = 0
        map%x(1:2, 1, :) = map%reference
        block
            real(dp) :: pos(3), d1(3, 3), d2(3, 3, 3), d3(3, 3, 3, 3), det
            call map%chart%jet([map%reference(:, 1), 0.0_dp], pos, d1, d2, d3)
            det = d1(1, 1)*(d1(2, 2)*d1(3, 3) - d1(2, 3)*d1(3, 2)) &
                - d1(1, 2)*(d1(2, 1)*d1(3, 3) - d1(2, 3)*d1(3, 1)) &
                + d1(1, 3)*(d1(2, 1)*d1(3, 2) - d1(2, 2)*d1(3, 1))
            if (.not. abs(det) > 0) return
            map%sigma = sign(1.0_dp, det)
        end block
        ! Gauge node: the interior node nearest the reference centroid.
        map%gauge = 0
        block
            real(dp) :: centre(2), best, distance
            integer :: i
            centre = sum(map%reference, dim=2)/map%nodes
            best = huge(1.0_dp)
            do i = 1, map%nodes
                if (map%boundary(i)) cycle
                distance = norm2(map%reference(:, i) - centre)
                if (distance < best) then
                    best = distance; map%gauge = i
                end if
            end do
        end block
        if (map%gauge == 0) return
        allocate (map%fixed(map%nred))
        map%fixed = .false.
        do a = 1, map%nodes
            if (.not. map%boundary(a)) cycle
            do m = 1, map%nmode
                map%fixed(ideal_map_rdof(map, 1, m, a)) = .true.
            end do
        end do
        map%fixed(ideal_map_rdof(map, 2, 1, map%gauge)) = .true.
        allocate (map%direction(2, map%nodes), map%axis_basis(map%npoly, map%nodes))
        map%direction(1, :) = 1; map%direction(2, :) = 0
        allocate (map%basis(3, 2, map%nodes))
        map%basis = 0
        map%basis(1, 1, :) = 1; map%basis(3, 2, :) = 1
        map%axis_basis = 0
        allocate (map%b(3, map%nq, map%cells), map%p(map%nq, map%cells), &
            map%det(map%nq, map%cells))
        map%b = 0; map%p = 0; map%det = 0
        status = IDEAL_MAP_OK
    end subroutine

    subroutine ideal_map_freeze(map, field, status)
        ! Pull the frozen reference data back to master coordinates:
        ! b^s = (-d chi/ds2, d chi/ds1), b^3 = b3*det(DG), with G the isoparametric
        ! reference geometry. Exactly closed for any label chi.
        type(ideal_map_t), intent(inout) :: map
        class(ideal_frozen_t), intent(in) :: field
        integer, intent(out) :: status
        real(dp) :: dg(2, 2), y(2), grad(2), ds(2), b3, pressure, chi, lo, hi, t
        real(dp), allocatable :: labels(:, :), slope(:)
        integer :: e, g, a, k
        status = IDEAL_MAP_OK
        allocate (labels(map%nodes, 1), map%axis_t(map%nodes), slope(map%nodes))
        do a = 1, map%nodes
            call field%evaluate(map%reference(:, a), chi, grad, b3, pressure)
            slope(a) = norm2(grad)
            if (norm2(grad) > 0) map%direction(:, a) = grad/norm2(grad)
            map%basis(1:2, 1, a) = map%direction(:, a)
            labels(a, 1) = chi
        end do
        ! Axis expansion: T in [0, 1], monomials about the extreme-label node,
        ! scaled by the reference radius.
        lo = sum(labels(:, 1), mask=map%boundary)/count(map%boundary)
        k = maxloc(abs(labels(:, 1) - lo), dim=1)
        hi = labels(k, 1)
        map%axis_node = k
        map%chi_node = labels(:, 1)
        do a = 1, map%nodes
            t = 0
            if (.not. map%boundary(a)) &
                t = min(1.0_dp, max(0.0_dp, (labels(a, 1) - lo)/(hi - lo)))
            map%axis_t(a) = t
        end do
        call axis_expansion(map, map%reference)
        do e = 1, map%cells
            do g = 1, map%nq
                y = matmul(map%reference(:, map%conn(:, e)), map%shape(:, g))
                dg = matmul(map%reference(:, map%conn(:, e)), &
                    transpose(map%dshape(:, :, g)))
                map%det(g, e) = dg(1, 1)*dg(2, 2) - dg(1, 2)*dg(2, 1)
                if (.not. map%det(g, e) > 0) status = IDEAL_MAP_FOLDED
                call field%evaluate(y, chi, grad, b3, pressure)
                ds = matmul(transpose(dg), grad)
                map%b(:, g, e) = [-ds(2), ds(1), b3*map%det(g, e)]
                map%p(g, e) = pressure
            end do
        end do
        call build_gauge(map)
        if (map%position_gauge) then
            map%omega_layers = 0
            call ideal_map_position_gauge(map)
        else if (map%omega_elimination) then
            call nested_axis_basis(map, slope, k)
            if (k /= IDEAL_MAP_OK) status = k
            call nested_omega_gauge(map)
        else
            map%omega_layers = 0
        end if
    end subroutine

    subroutine nested_omega_gauge(map)
        ! At n = 0 every omega = f(chi) rotates each surface rigidly in phi:
        ! b_pol . grad f(chi) = 0, so W is exactly invariant in the continuum
        ! for any f (a gauge of the same kind as the field-line slides of the
        ! axis basis; the gauge probe measures its discrete curvature vanishing
        ! with h). Discretely the smooth f are near-nulls of H_oo; omega* then
        ! jitters along them with the roundoff of g_omega and the nu gradient
        ! inherits it. They are removed from the unknowns by constraints, the
        ! L2 moments int phi_l(rho) omega dA = 0 (reference quadrature) of
        ! the C1 partition of unity phi_l = cos^2(pi/2 (L rho - l)),
        ! l = 0..L, rho = (1 - T)^1/2, L = omega_layers (one per element layer).
        ! The moments are smooth functionals, so a smooth f restores them for
        ! any continuum minimiser omega* + f: no physics is fixed and the
        ! order is kept; the f(chi) left free oscillate within a layer and are
        ! stiff. They enter the omega solves as Lagrange rows.
        type(ideal_map_t), intent(inout) :: map
        real(dp), allocatable :: gram(:, :), lam(:)
        logical, allocatable :: free(:)
        real(dp) :: rho, w, x, pi
        integer :: e, g, l, nc, k, info
        pi = acos(-1.0_dp)
        if (map%omega_layers < 0) map%omega_layers = max(1, &
            nint(sqrt(map%cells*sqrt(3.0_dp)/(4*pi))))
        if (map%omega_layers == 0) return
        nc = map%omega_layers + 1
        ! The moments exclude the constant (sum phi_l = 1): the omega pin, a
        ! point condition that would force layer-scale f, is released.
        map%fixed(ideal_map_rdof(map, 2, 1, map%gauge)) = .false.
        allocate (map%omega_c(nc, map%nodes))
        map%omega_c = 0
        do e = 1, map%cells
            do g = 1, map%nq
                rho = map%omega_layers*sqrt(max(0.0_dp, &
                    1 - dot_product(map%axis_t(map%conn(:, e)), map%shape(:, g))))
                w = map%weight(g)*map%det(g, e)
                do l = max(0, floor(rho)), min(nc - 1, ceiling(rho))
                    x = rho - l
                    if (abs(x) >= 1) cycle
                    map%omega_c(l + 1, map%conn(:, e)) = map%omega_c(l + 1, &
                        map%conn(:, e)) + w*cos(pi*x/2)**2*map%shape(:, g)
                end do
            end do
        end do
        ! Scale: unit largest weight.
        map%omega_c = map%omega_c/maxval(abs(map%omega_c))
        allocate (free(map%nodes))
        do k = 1, map%nodes
            free(k) = .not. map%fixed(ideal_map_rdof(map, 2, 1, k))
        end do
        gram = matmul(map%omega_c, transpose(map%omega_c) &
            *spread(merge(1.0_dp, 0.0_dp, free), 2, nc))
        call symmetric_eigen(gram, lam, info)
        map%omega_cgi = matmul(gram, transpose(gram)/spread(lam, 2, nc))
    end subroutine

    subroutine ideal_map_gauge_project(map, gred)
        ! Remove the Lagrange part C lambda from the free omega gradient
        ! (orthogonal projection onto ker C^T): the gradient of W on the
        ! constrained omega space.
        type(ideal_map_t), intent(in) :: map
        real(dp), intent(inout) :: gred(:)
        real(dp), allocatable :: go(:)
        integer :: a
        if (.not. allocated(map%omega_c)) return
        allocate (go(map%nodes))
        do a = 1, map%nodes
            go(a) = gred(ideal_map_rdof(map, 2, 1, a))
            if (map%fixed(ideal_map_rdof(map, 2, 1, a))) go(a) = 0
        end do
        go = go - matmul(matmul(matmul(map%omega_cgi, map%omega_c), go), map%omega_c)
        do a = 1, map%nodes
            if (map%fixed(ideal_map_rdof(map, 2, 1, a))) cycle
            gred(ideal_map_rdof(map, 2, 1, a)) = go(a)
        end do
    end subroutine

    subroutine nested_axis_basis(map, slope, status)
        ! Axis basis of the nested solve, free of relabelling directions.
        ! The nodal span {nu_a d_a} + span{e_i T P_j} equals {nu_a d_a} + the
        ! tangential parts t_a (t_a . e_i) T P_j (t_a = d_a rotated by +90 deg),
        ! so only tangential content h_a is parametrised. Of the 2 npoly
        ! polynomial fields of degree <= K - 1 two kinds of directions are not
        ! unknowns (ng = npoly - K of each, leaving 2K):
        ! (1) near-normal fields, duplicated by nu: the ng directions with the
        !     smallest tangential fraction |h|/|F| (circular surfaces: the
        !     dilation u and the radial fields P u, deg P <= K - 2);
        ! (2) field-line slides: any v = g b_pol, b_pol = |grad chi| t, with
        !     omega* following preserves the frozen flux and pressure, so W_red
        !     is exactly invariant in the continuum. Smooth slides of the axis
        !     degree, g = T P (deg P <= K - 2: the T-weighted rotation, its
        !     u-linear 1/R part and the r^2-weighted translation), are removed:
        !     the ng remaining directions with the largest cosine to them.
        ! Kept: the 2K harmonics r^d e^{i(d+1) theta} for circular surfaces
        ! (shift, elongation, triangularity). Contents use the nodal sum.
        type(ideal_map_t), intent(inout) :: map
        real(dp), intent(in) :: slope(:)
        integer, intent(out) :: status
        real(dp), allocatable :: f(:, :), h(:, :), mm(:, :), w(:, :), lam(:), y(:, :)
        real(dp), allocatable :: c(:, :), q(:, :), mu(:)
        real(dp) :: t(2)
        integer :: n, np2, kd, ng, nk, a, p, i, j
        n = map%nodes; np2 = 2*map%npoly
        kd = nint((sqrt(8.0_dp*map%npoly + 1) - 1)/2)
        ng = map%npoly - kd; nk = np2 - 2*ng
        allocate (f(2*n, np2), h(n, np2))
        do a = 1, n
            t = [-map%direction(2, a), map%direction(1, a)]
            do j = 1, map%npoly
                do i = 1, 2
                    p = i + 2*(j - 1)
                    f(2*a - 1:2*a, p) = 0
                    f(2*a - 2 + i, p) = map%axis_basis(j, a)
                    h(a, p) = t(i)*map%axis_basis(j, a)
                end do
            end do
        end do
        ! Whiten with the field Gram matrix: |h y| / |y| is the tangential fraction.
        mm = matmul(transpose(f), f)
        call symmetric_eigen(mm, lam, status)
        if (status /= IDEAL_MAP_OK) return
        if (.not. minval(lam) > 0) then
            status = IDEAL_MAP_BAD_INPUT; return
        end if
        w = matmul(mm, matmul(diag_matrix(1/sqrt(lam)), transpose(mm)))
        h = matmul(h, w)
        y = matmul(transpose(h), h)
        call symmetric_eigen(y, lam, status)
        if (status /= IDEAL_MAP_OK) return
        allocate (map%axis_removed(2, ng))
        map%axis_removed(1, :) = lam(1:ng)
        ! Orthonormal contents of the kept near-tangential directions.
        c = matmul(h, matmul(y(:, ng + 1:), diag_matrix(1/sqrt(lam(ng + 1:)))))
        if (ng > 0) then
            ! Orthonormal slide contents T P_k |grad chi|, k <= ng.
            allocate (q(n, ng))
            do j = 1, ng
                q(:, j) = map%axis_basis(j, :)*slope
            end do
            mm = matmul(transpose(q), q)
            call symmetric_eigen(mm, lam, status)
            if (status /= IDEAL_MAP_OK) return
            q = matmul(q, matmul(mm, diag_matrix(1/sqrt(lam))))
            mm = matmul(transpose(c), q)
            mm = matmul(mm, transpose(mm))
            call symmetric_eigen(mm, mu, status)
            if (status /= IDEAL_MAP_OK) return
            map%axis_removed(2, :) = mu(nk + 1:)
            c = matmul(c, mm(:, :nk))
        end if
        allocate (map%axis_vec(2, np2, n))
        map%axis_vec = 0
        do p = 1, nk
            c(:, p) = c(:, p)/maxval(abs(c(:, p)))
            do a = 1, n
                map%axis_vec(:, p, a) = [-map%direction(2, a), map%direction(1, a)] &
                    *c(a, p)
            end do
        end do
        do p = 1, map%nmode
            do j = nk + 1, np2
                map%fixed(ideal_map_axis_dof(map, 1, 1, p) + j - 1) = .true.
            end do
        end do
    contains
        pure function diag_matrix(v) result(d)
            real(dp), intent(in) :: v(:)
            real(dp) :: d(size(v), size(v))
            integer :: k
            d = 0
            do k = 1, size(v)
                d(k, k) = v(k)
            end do
        end function
    end subroutine

    subroutine symmetric_eigen(a, w, status)
        ! Eigenvalues ascending in w, orthonormal eigenvectors overwrite a.
        real(dp), intent(inout) :: a(:, :)
        real(dp), allocatable, intent(out) :: w(:)
        integer, intent(out) :: status
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
        real(dp), allocatable :: work(:)
        integer :: n, info
        n = size(a, 1)
        allocate (w(n), work(max(1, 66*n)))
        call dsyev('V', 'U', n, a, n, w, work, size(work), info)
        status = IDEAL_MAP_OK
        if (info /= 0) status = IDEAL_MAP_BAD_INPUT
    end subroutine

    subroutine ideal_map_position_gauge(map)
        ! Position gauge (as N4): full (q1, q2) per interior node and mode, no
        ! toroidal displacement (omega = 0, toroidal labels follow phi-hat). The
        ! poloidal slide along surfaces is then physical (it redistributes the
        ! toroidal flux) instead of a relabelling null, rigid motions are
        ! exact, and neither the axis expansion nor an omega pin is needed.
        type(ideal_map_t), intent(inout) :: map
        integer :: a, m, j
        if (allocated(map%axis_vec)) deallocate (map%axis_vec)
        map%axis_basis = 0
        do a = 1, map%nodes
            map%basis(:, 1, a) = [1.0_dp, 0.0_dp, 0.0_dp]
            map%basis(:, 2, a) = [0.0_dp, 1.0_dp, 0.0_dp]
            do m = 1, map%nmode
                map%fixed(ideal_map_rdof(map, 1, m, a)) = map%boundary(a)
                map%fixed(ideal_map_rdof(map, 2, m, a)) = map%boundary(a)
            end do
        end do
        do m = 1, map%nmode
            do j = 1, map%npoly
                map%fixed(ideal_map_axis_dof(map, 1, j, m)) = .true.
                map%fixed(ideal_map_axis_dof(map, 2, j, m)) = .true.
            end do
        end do
        map%omega_gauge = 0
    end subroutine

    subroutine axis_expansion(map, position)
        ! axis_basis(j, a) = T_a P_j((y_a - y_axis)/radius) on node positions.
        type(ideal_map_t), intent(inout) :: map
        real(dp), intent(in) :: position(:, :)
        real(dp) :: centre(2), radius, u(2)
        integer :: a, j, d, i
        centre = position(:, map%axis_node)
        radius = maxval(norm2(position - spread(centre, 2, map%nodes), dim=1))
        do a = 1, map%nodes
            u = (position(:, a) - centre)/radius
            j = 0
            do d = 0, nint((sqrt(8.0_dp*map%npoly + 1) - 1)/2) - 1
                do i = 0, d
                    j = j + 1
                    map%axis_basis(j, a) = map%axis_t(a)**map%axis_power &
                        *u(1)**(d - i)*u(2)**i
                end do
            end do
        end do
    end subroutine

    subroutine ideal_map_prolong(fine, coarse, status)
        ! Grid sequencing (VMEC's ns ladder): initialise the coefficients of a
        ! map on a finer reference mesh of the same domain and labels from a
        ! solved coarse map, evaluated at the fine reference nodes. Boundary
        ! nodes keep their reference position (held dofs). Re-gauging the
        ! transversals to the coarse surfaces (cell-averaged F) was tried and
        ! did not change the observed order.
        type(ideal_map_t), intent(inout) :: fine
        type(ideal_map_t), intent(in) :: coarse
        integer, intent(out) :: status
        real(dp) :: xi(2), shape(10), grad(2, 10), y(2)
        integer :: a, e, k, owner, located, j
        status = IDEAL_MAP_BAD_INPUT
        if (fine%nmode /= coarse%nmode .or. fine%nfp /= coarse%nfp) return
        owner = 1
        do a = 1, fine%nodes
            y = fine%reference(:, a)
            located = 0
            do k = 0, coarse%cells - 1
                e = 1 + mod(owner - 1 + merge(k/2, coarse%cells - (k + 1)/2, &
                    mod(k, 2) == 0), coarse%cells)
                call triangle_locate(coarse%reference(:, coarse%conn(:, e)), y, xi, j)
                if (j /= 0) cycle
                if (minval([xi, 1 - sum(xi)]) < -1.0e-10_dp) cycle
                located = e; exit
            end do
            if (located == 0) then
                if (.not. fine%boundary(a)) return
                cycle
            end if
            owner = located
            call triangle_basis_jet(xi, shape(:coarse%nloc), grad(:, :coarse%nloc))
            do k = 1, fine%nmode
                fine%x(:, k, a) = matmul(coarse%x(:, k, coarse%conn(:, located)), &
                    shape(:coarse%nloc))
            end do
        end do
        do a = 1, fine%nodes
            if (.not. fine%boundary(a)) cycle
            fine%x(1:2, 1, a) = fine%reference(:, a)
            fine%x(1:2, 2:, a) = 0
        end do
        fine%x(3, 1, :) = fine%x(3, 1, :) - fine%x(3, 1, fine%gauge)
        status = IDEAL_MAP_OK
    end subroutine

    subroutine build_gauge(map)
        ! Reference-metric form of the omega gauge term; weight relative to the
        ! largest physical omega stiffness R^2 |b_pol|^2/det G of the reference.
        type(ideal_map_t), intent(inout) :: map
        real(dp) :: dg(2, 2), ginv(2, 2), y(2), chi_s(2), v(10), scale, w, metric, h33
        real(dp) :: pos(3), d1(3, 3), d2(3, 3, 3), d3(3, 3, 3, 3)
        integer :: e, g, n, j
        n = map%nloc
        allocate (map%omega_gauge(n, n, map%cells))
        map%omega_gauge = 0
        ! The phi-fixed increment space has no omega unknown. Its gauge
        ! matrix is identically zero; avoid assembling a form discarded by
        ! ideal_map_position_gauge immediately after this call.
        if (map%position_gauge) return
        scale = 0
        !$omp parallel do default(shared) schedule(static) &
        !$omp private(e,g,dg,ginv,y,chi_s,v,w,metric,h33,pos,d1,d2,d3,j) &
        !$omp reduction(max:scale)
        do e = 1, map%cells
            do g = 1, map%nq
                dg = matmul(map%reference(:, map%conn(:, e)), &
                    transpose(map%dshape(:, :, g)))
                y = matmul(map%reference(:, map%conn(:, e)), map%shape(:, g))
                ! Inverse metric of s -> y; chi_s = (b2, -b1) in master coordinates.
                ginv = matmul(transpose(dg), dg)
                ginv = reshape([ginv(2, 2), -ginv(2, 1), -ginv(1, 2), ginv(1, 1)], &
                    [2, 2])/(ginv(1, 1)*ginv(2, 2) - ginv(1, 2)*ginv(2, 1))
                chi_s = [map%b(2, g, e), -map%b(1, g, e)]
                v(:n) = map%det(g, e)*matmul(matmul(chi_s, ginv), map%dshape(:, :n, g))
                metric = 0
                do j = 1, map%nphi
                    call map%chart%jet([y, map%phi(j)], pos, d1, d2, d3)
                    call curve_omega_metric(d1(1, 3), d1(2, 3), d1(3, 3), h33)
                    metric = metric + h33/map%nphi
                end do
                w = map%weight(g)*metric/map%det(g, e)
                scale = max(scale, metric*sum(map%b(1:2, g, e)**2)/map%det(g, e)**2)
                map%omega_gauge(:, :, e) = map%omega_gauge(:, :, e) &
                    + w*spread(v(:n), 2, n)*spread(v(:n), 1, n)
            end do
        end do
        !$omp end parallel do
        map%omega_gauge = map%gauge_weight*map%omega_gauge/max(scale, tiny(1.0_dp))
        if (map%omega_elimination) map%omega_gauge = 0
    end subroutine

    subroutine point_jet(map, e, g, j, x, q, qd)
        ! q(3) and Q(k, i) = dq_k/ds_i at quadrature point g, plane j.
        type(ideal_map_t), intent(in) :: map
        real(dp), intent(in) :: x(:, :, :)
        integer, intent(in) :: e, g, j
        real(dp), intent(out) :: q(3), qd(3, 3)
        real(dp) :: v
        integer :: a, m, node
        q = 0; qd = 0
        do a = 1, map%nloc
            node = map%conn(a, e)
            do m = 1, map%nmode
                v = map%trig(m, j)
                q = q + x(:, m, node)*map%shape(a, g)*v
                qd(:, 1) = qd(:, 1) + x(:, m, node)*map%dshape(1, a, g)*v
                qd(:, 2) = qd(:, 2) + x(:, m, node)*map%dshape(2, a, g)*v
                qd(:, 3) = qd(:, 3) + x(:, m, node)*map%shape(a, g)*map%dtrig(m, j)
            end do
        end do
        q(3) = q(3) + map%phi(j)
        qd(3, 3) = qd(3, 3) + 1
    end subroutine

    subroutine ideal_map_point(map, e, xi, phi, q, y)
        ! Chart coordinates q of the map at master point xi of cell e and angle
        ! phi-hat, and the reference point y.
        type(ideal_map_t), intent(in) :: map
        integer, intent(in) :: e
        real(dp), intent(in) :: xi(2), phi
        real(dp), intent(out) :: q(3), y(2)
        real(dp) :: shape(10), grad(2, 10), t
        integer :: a, k
        call triangle_basis_jet(xi, shape(:map%nloc), grad(:, :map%nloc))
        q = 0
        do a = 1, map%nloc
            q = q + map%x(:, 1, map%conn(a, e))*shape(a)
            do k = 1, map%ntor
                t = k*map%nfp*phi
                q = q + (map%x(:, 2*k, map%conn(a, e))*cos(t) &
                    + map%x(:, 2*k + 1, map%conn(a, e))*sin(t))*shape(a)
            end do
        end do
        q(3) = q(3) + phi
        y = matmul(map%reference(:, map%conn(:, e)), shape(:map%nloc))
    end subroutine

    subroutine point_gradient(map, e, g, j, x, w, jac, gv)
        ! Energy density w, Jacobian and the gradient with respect to
        ! v = (q(1:3), Q(1:3, 1:3)) (index 3 + k + 3*(i-1) for Q(k,i)): the
        ! generated chain-rule leaf (gen/gen_ideal_map.f90) on the chart jet.
        type(ideal_map_t), intent(in) :: map
        integer, intent(in) :: e, g, j
        real(dp), intent(in) :: x(:, :, :)
        real(dp), intent(out) :: w, jac, gv(12)
        real(dp) :: q(3), qd(3, 3), pos(3), d1(3, 3), d2(3, 3, 3), d3(3, 3, 3, 3)
        call point_jet(map, e, g, j, x, q, qd)
        call map%chart%jet(q, pos, d1, d2, d3)
        call ideal_point_gradient_jets(d1, d2, d3, qd, map%b(:, g, e), map%p(g, e), &
            map%mu0, map%sigma, w, jac, gv)
    end subroutine

    subroutine point_hessian(map, e, g, j, x, w, jac, gv, hp)
        ! As point_gradient with the Hessian with respect to v, packed in the
        ! group order of gen/gen_ideal_map.f90: hp(1:45) the upper triangle of
        ! the block of v(1:9), hp(46:72) the rows v(1:9) of the columns
        ! v(10:12), hp(73:78) the upper triangle of the block of v(10:12).
        type(ideal_map_t), intent(in) :: map
        integer, intent(in) :: e, g, j
        real(dp), intent(in) :: x(:, :, :)
        real(dp), intent(out) :: w, jac, gv(12), hp(78)
        real(dp) :: q(3), qd(3, 3), pos(3), d1(3, 3), d2(3, 3, 3), d3(3, 3, 3, 3)
        call point_jet(map, e, g, j, x, q, qd)
        call map%chart%jet(q, pos, d1, d2, d3)
        call ideal_point_hessian_jets(d1, d2, d3, qd, map%b(:, g, e), map%p(g, e), &
            map%mu0, map%sigma, w, jac, gv, hp)
    end subroutine

    pure subroutine test_functions(map, g, j, a, m, f)
        ! f(0:3) = (N T, dN/ds1 T, dN/ds2 T, N dT/ds3) of node a, mode m.
        type(ideal_map_t), intent(in) :: map
        integer, intent(in) :: g, j, a, m
        real(dp), intent(out) :: f(0:3)
        f(0) = map%shape(a, g)*map%trig(m, j)
        f(1) = map%dshape(1, a, g)*map%trig(m, j)
        f(2) = map%dshape(2, a, g)*map%trig(m, j)
        f(3) = map%shape(a, g)*map%dtrig(m, j)
    end subroutine

    subroutine ideal_map_energy(map, x, energy, status, gradient, min_jacobian)
        ! Energy per field period and its coefficient gradient. status
        ! IDEAL_MAP_FOLDED when J <= 0 at any quadrature point.
        type(ideal_map_t), intent(in) :: map
        real(dp), intent(in) :: x(:, :, :)
        real(dp), intent(out) :: energy
        integer, intent(out) :: status
        real(dp), intent(out), optional :: gradient(:), min_jacobian
        real(dp), allocatable :: local(:, :), cell_energy(:), cell_jac(:)
        real(dp) :: w, jac, gv(12), f(0:3), dphi, weight, minjac
        integer :: e, g, j, a, m, c, l
        logical :: want
        real(dp) :: pf1_t0
        pf1_t0 = wall_seconds()
        want = present(gradient)
        allocate (local(3*map%nmode*map%nloc, map%cells), cell_energy(map%cells), &
            cell_jac(map%cells))
        dphi = 2*acos(-1.0_dp)/(map%nfp*map%nphi)
        ! Cell sums are combined serially below, so the energy (which decides
        ! step acceptance near roundoff) does not depend on the thread count
        ! or on the order of an OpenMP reduction.
        !$omp parallel do default(none) shared(map, x, local, want, dphi, &
        !$omp cell_energy, cell_jac) &
        !$omp private(e, g, j, a, m, c, l, w, jac, gv, f, weight) schedule(static)
        do e = 1, map%cells
            local(:, e) = 0
            cell_energy(e) = 0; cell_jac(e) = huge(1.0_dp)
            do j = 1, map%nphi
                do g = 1, map%nq
                    call point_gradient(map, e, g, j, x, w, jac, gv)
                    cell_jac(e) = min(cell_jac(e), jac/map%det(g, e))
                    weight = map%weight(g)*dphi
                    cell_energy(e) = cell_energy(e) + weight*w
                    if (.not. want) cycle
                    l = 0
                    do a = 1, map%nloc
                        do m = 1, map%nmode
                            call test_functions(map, g, j, a, m, f)
                            do c = 1, 3
                                l = l + 1
                                local(l, e) = local(l, e) + weight*(gv(c)*f(0) &
                                    + gv(3 + c)*f(1) + gv(6 + c)*f(2) + gv(9 + c)*f(3))
                            end do
                        end do
                    end do
                end do
            end do
        end do
        !$omp end parallel do
        energy = 0
        do e = 1, map%cells
            energy = energy + cell_energy(e)
        end do
        minjac = minval(cell_jac)
        if (.not. map%gauge_hessian_only) then
            do e = 1, map%cells
                energy = energy + 0.5_dp*dot_product(x(3, 1, map%conn(:, e)), &
                    matmul(map%omega_gauge(:, :, e), x(3, 1, map%conn(:, e))))
            end do
        end if
        status = IDEAL_MAP_OK
        if (.not. (minjac > 0) .or. .not. ieee_is_finite(energy)) &
            status = IDEAL_MAP_FOLDED
        if (present(min_jacobian)) min_jacobian = minjac
        if (.not. want) then
            pf1_t(1) = pf1_t(1) + wall_seconds() - pf1_t0; pf1_n(1) = pf1_n(1) + 1
            return
        end if
        gradient = 0
        do e = 1, map%cells
            gradient(map%dofconn(:, e)) = gradient(map%dofconn(:, e)) + local(:, e)
            if (map%gauge_hessian_only) cycle
            do a = 1, map%nloc
                c = ideal_map_dof(map, 3, 1, map%conn(a, e))
                gradient(c) = gradient(c) + dot_product(map%omega_gauge(a, :, e), &
                    x(3, 1, map%conn(:, e)))
            end do
        end do
        pf1_t(2) = pf1_t(2) + wall_seconds() - pf1_t0; pf1_n(2) = pf1_n(2) + 1
    end subroutine

    subroutine ideal_map_hessian(map, x, matrix)
        ! Assembled coefficient Hessian (all Fourier modes coupled) as CSR. The
        ! sparsity pattern of a matrix already built for this map is reused.
        ! Per cell the weighted point Hessians hv(12, 12) of v = (q, Q) are
        ! contracted in stages instead of one rank update per point:
        !   K[(r,c),(r',c')](m,m',g) = sum_j T_r(m,j) T_r'(m',j) hv_j      (planes)
        !   full[(a,m,c),(a',m',c')] = sum_g,s,s' phi_s(a,g) phi_s'(a',g) K_ss'
        !                                                                   (space)
        !   blocks = P^T full P, reflected                              (reduction)
        ! with r = 0, 1, 2 for (N T, dN/ds1 T, dN/ds2 T) and r = 3 for N dT/ds3.
        ! The plane stage is one DGEMM per trig/dtrig pairing over the map-level
        ! products pphi. The space stage merges r = 0 and 3 (both carry N) into
        ! s = 1, 2, 3 and runs two DGEMM families. Workspaces are thread-private
        ! and allocated once per parallel region.
        type(ideal_map_t), intent(in) :: map
        real(dp), intent(in) :: x(:, :, :)
        type(csr_matrix_t), intent(inout) :: matrix
        interface
            subroutine dgemm(ta, tb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
                import :: dp
                character, intent(in) :: ta, tb
                integer, intent(in) :: m, n, k, lda, ldb, ldc
                real(dp), intent(in) :: alpha, beta, a(lda, *), b(ldb, *)
                real(dp), intent(inout) :: c(ldc, *)
            end subroutine
        end interface
        real(dp), allocatable :: blocks(:, :, :), phig(:, :, :), phik(:, :)
        real(dp), allocatable :: hwa(:, :), hwb(:, :), hwc(:, :), ka(:, :), kb(:, :)
        real(dp), allocatable :: kc(:, :), ks(:, :), yw(:, :), ft(:, :)
        real(dp), allocatable :: full(:, :), half(:, :)
        real(dp) :: w, jac, gv(12), hp(78), dphi, weight
        integer :: e, g, j, a, c, l, k, r, nl, nr, nlm, nmm, nmc, nloc, nq, nphi
        integer :: nmode
        real(dp) :: pf1_t0
        pf1_t0 = wall_seconds()
        nloc = map%nloc; nq = map%nq; nphi = map%nphi; nmode = map%nmode
        nlm = nloc*nmode
        nl = 3*nlm
        nr = 2*nlm + 2*map%npoly*nmode
        nmm = nmode**2
        nmc = 9*nmm
        allocate (blocks(nr, nr, map%cells), phig(3, nloc, nq), phik(3*nq, nloc))
        dphi = 2*acos(-1.0_dp)/(map%nfp*map%nphi)
        ! phi_s(a, g) for s = 1, 2, 3: N, dN/ds1, dN/ds2.
        do g = 1, nq
            do a = 1, nloc
                phig(1, a, g) = map%shape(a, g)
                phig(2, a, g) = map%dshape(1, a, g)
                phig(3, a, g) = map%dshape(2, a, g)
                do k = 1, 3
                    phik(k + 3*(g - 1), a) = phig(k, a, g)
                end do
            end do
        end do
        !$omp parallel default(none) shared(map, x, blocks, phig, phik, dphi, nl, nr, &
        !$omp nlm, nmm, nmc, nloc, nq, nphi, nmode) &
        !$omp private(e, g, j, a, c, l, k, r, w, jac, gv, hp, weight, hwa, hwb, hwc, &
        !$omp ka, kb, kc, ks, yw, ft, full, half)
        allocate (hwa(nphi, 45*nq), hwb(nphi, 27*nq), hwc(nphi, 6*nq), &
            ka(nmm, 45*nq), kb(nmm, 27*nq), kc(nmm, 6*nq), ks(3*nmc, 3*nq), &
            yw(nmc, 3*nq*nloc), ft(nmc, nloc*nloc), full(nl, nl), half(nl, nr))
        !$omp do schedule(static)
        do e = 1, map%cells
            ! Weighted point Hessians; column (pair, g) of the plane stage input.
            do g = 1, nq
                weight = map%weight(g)*dphi
                do j = 1, nphi
                    call point_hessian(map, e, g, j, x, w, jac, gv, hp)
                    do k = 1, 45
                        hwa(j, k + 45*(g - 1)) = weight*hp(k)
                    end do
                    do k = 1, 27
                        hwb(j, k + 27*(g - 1)) = weight*hp(45 + k)
                    end do
                    do k = 1, 6
                        hwc(j, k + 6*(g - 1)) = weight*hp(72 + k)
                    end do
                end do
            end do
            ! Plane stage: K(m + nmode (m' - 1), pair, g) per pairing.
            call dgemm('T', 'N', nmm, 45*nq, nphi, 1.0_dp, map%pphi(1, 1, 1), nphi, &
                hwa, nphi, 0.0_dp, ka, nmm)
            call dgemm('T', 'N', nmm, 27*nq, nphi, 1.0_dp, map%pphi(1, 1, 2), nphi, &
                hwb, nphi, 0.0_dp, kb, nmm)
            call dgemm('T', 'N', nmm, 6*nq, nphi, 1.0_dp, map%pphi(1, 1, 3), nphi, &
                hwc, nphi, 0.0_dp, kc, nmm)
            call hessian_merge_pairs(nmode, nq, ka, kb, kc, ks)
            ! Space stage: contract s' with phi_s'(a', g), then (s, g) with
            ! phi_s(a, g); yw(mmcc, s + 3 (g - 1) + 3 nq (a' - 1)).
            do g = 1, nq
                call dgemm('N', 'N', 3*nmc, nloc, 3, 1.0_dp, ks(1, 3*(g - 1) + 1), &
                    3*nmc, phig(1, 1, g), 3, 0.0_dp, yw(1, 3*(g - 1) + 1), &
                    3*nmc*nq)
            end do
            do a = 1, nloc
                call dgemm('N', 'N', nmc, nloc, 3*nq, 1.0_dp, &
                    yw(1, 3*nq*(a - 1) + 1), nmc, phik, 3*nq, 0.0_dp, &
                    ft(1, nloc*(a - 1) + 1), nmc)
            end do
            call hessian_scatter(nmode, nloc, ft, full)
            do a = 1, nloc
                do c = 1, nloc
                    l = 3*(1 + nmode*(a - 1))
                    r = 3*(1 + nmode*(c - 1))
                    full(l, r) = full(l, r) + map%omega_gauge(a, c, e)
                end do
            end do
            call hessian_reduce(map, e, full, half, blocks(:, :, e))
        end do
        !$omp end do
        deallocate (hwa, hwb, hwc, ka, kb, kc, ks, yw, ft, full, half)
        !$omp end parallel
        pf1_t(3) = pf1_t(3) + wall_seconds() - pf1_t0; pf1_n(3) = pf1_n(3) + 1
        pf1_t0 = wall_seconds()
        if (matrix%n /= map%nred .or. matrix%nloc /= nr) &
            call csr_build_pattern(map%nred, map%redconn, .false., matrix)
        call csr_assemble_blocks(matrix, map%redconn, blocks)
        pf1_t(4) = pf1_t(4) + wall_seconds() - pf1_t0; pf1_n(4) = pf1_n(4) + 1
    end subroutine

    subroutine hessian_merge_pairs(nmode, nq, ka, kb, kc, ks)
        ! Packed pair blocks K_uu'(m, m') (u, u' = 1..12 with u = c + 3 r) to
        ! ks(mmcc + nmc (s - 1), s' + 3 (g - 1)), mmcc = m + nmode (m' - 1)
        ! + nmm (c - 1 + 3 (c' - 1)). The lower triangle is the transposed
        ! upper one, K_u'u(m', m) = K_uu'(m, m'), and r = 0 and r = 3 both
        ! carry N, so they share s = 1.
        integer, intent(in) :: nmode, nq
        real(dp), intent(in) :: ka(:, :), kb(:, :), kc(:, :)
        real(dp), intent(out) :: ks(:, :)
        integer :: g, u, v, nmm, nmc
        nmm = nmode**2
        nmc = 9*nmm
        ks = 0
        do g = 1, nq
            do v = 1, 9
                do u = 1, v
                    call add_pair(u, v, ka(:, u + v*(v - 1)/2 + 45*(g - 1)))
                end do
            end do
            do v = 10, 12
                do u = 1, 9
                    call add_pair(u, v, kb(:, u + 9*(v - 10) + 27*(g - 1)))
                end do
                do u = 10, v
                    call add_pair(u, v, kc(:, u - 9 + (v - 10)*(v - 9)/2 + 6*(g - 1)))
                end do
            end do
        end do
    contains
        subroutine add_pair(u, v, kv)
            integer, intent(in) :: u, v
            real(dp), intent(in) :: kv(:)
            integer :: m, mp, cu, cv, su, sv, row
            cu = mod(u - 1, 3) + 1
            cv = mod(v - 1, 3) + 1
            su = merge(1, (u - 1)/3 + 1, (u - 1)/3 == 3)
            sv = merge(1, (v - 1)/3 + 1, (v - 1)/3 == 3)
            do mp = 1, nmode
                do m = 1, nmode
                    row = m + nmode*(mp - 1) + nmm*(cu - 1 + 3*(cv - 1)) + nmc*(su - 1)
                    ks(row, sv + 3*(g - 1)) = ks(row, sv + 3*(g - 1)) &
                        + kv(m + nmode*(mp - 1))
                end do
            end do
            if (u == v) return
            do mp = 1, nmode
                do m = 1, nmode
                    row = mp + nmode*(m - 1) + nmm*(cv - 1 + 3*(cu - 1)) + nmc*(sv - 1)
                    ks(row, su + 3*(g - 1)) = ks(row, su + 3*(g - 1)) &
                        + kv(m + nmode*(mp - 1))
                end do
            end do
        end subroutine
    end subroutine

    subroutine hessian_scatter(nmode, nloc, ft, full)
        ! ft(mmcc, a + nloc (a' - 1)) to full(c + 3 (m - 1) + 3 nmode (a - 1),
        ! c' + 3 (m' - 1) + 3 nmode (a' - 1)).
        integer, intent(in) :: nmode, nloc
        real(dp), intent(in) :: ft(:, :)
        real(dp), intent(out) :: full(:, :)
        integer :: a, ap, c, cp, m, mp, nmm
        nmm = nmode**2
        do ap = 1, nloc
            do a = 1, nloc
                do cp = 1, 3
                    do c = 1, 3
                        do mp = 1, nmode
                            do m = 1, nmode
                                full(c + 3*(m - 1) + 3*nmode*(a - 1), &
                                    cp + 3*(mp - 1) + 3*nmode*(ap - 1)) = &
                                    ft(m + nmode*(mp - 1) + nmm*(c - 1 + 3*(cp - 1)), &
                                    a + nloc*(ap - 1))
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end subroutine

    subroutine hessian_reduce(map, e, full, half, blk)
        ! Cell block P^T full P in the reduced (nu, omega, axis) dofs, reflected
        ! from the lower triangle. P maps (nu, omega) of local node-mode l to
        ! (d1 nu, d2 nu, omega) through map%basis and the axis modes to the
        ! components of the node vectors.
        type(ideal_map_t), intent(in) :: map
        integer, intent(in) :: e
        real(dp), intent(in) :: full(:, :)
        real(dp), intent(out) :: half(:, :), blk(:, :)
        real(dp) :: b1(3), b2(3), av(2), t, s
        integer :: i, l, a, node, m, jp, c, r, nlm, nl, nr, nmode, nloc
        logical :: nested
        nmode = map%nmode; nloc = map%nloc
        nlm = nloc*nmode
        nl = 3*nlm
        nr = 2*nlm + 2*map%npoly*nmode
        nested = allocated(map%axis_vec)
        do l = 1, nlm
            node = map%conn(1 + (l - 1)/nmode, e)
            b1 = map%basis(:, 1, node)
            b2 = map%basis(:, 2, node)
            do i = 1, nl
                half(i, 2*l - 1) = b1(1)*full(i, 3*l - 2) + b1(2)*full(i, 3*l - 1) &
                    + b1(3)*full(i, 3*l)
                half(i, 2*l) = b2(1)*full(i, 3*l - 2) + b2(2)*full(i, 3*l - 1) &
                    + b2(3)*full(i, 3*l)
            end do
        end do
        do m = 1, nmode
            do jp = 1, map%npoly
                do c = 1, 2
                    r = 2*nlm + c + 2*(jp - 1) + 2*map%npoly*(m - 1)
                    half(:, r) = 0
                    do a = 1, nloc
                        l = m + nmode*(a - 1)
                        node = map%conn(a, e)
                        if (nested) then
                            av = map%axis_vec(:, c + 2*(jp - 1), node)
                            do i = 1, nl
                                half(i, r) = half(i, r) + av(1)*full(i, 3*l - 2) &
                                    + av(2)*full(i, 3*l - 1)
                            end do
                        else
                            t = map%axis_basis(jp, node)
                            do i = 1, nl
                                half(i, r) = half(i, r) + t*full(i, 3*l - 3 + c)
                            end do
                        end if
                    end do
                end do
            end do
        end do
        do c = 1, nr
            do l = 1, nlm
                node = map%conn(1 + (l - 1)/nmode, e)
                b1 = map%basis(:, 1, node)
                b2 = map%basis(:, 2, node)
                blk(2*l - 1, c) = b1(1)*half(3*l - 2, c) + b1(2)*half(3*l - 1, c) &
                    + b1(3)*half(3*l, c)
                blk(2*l, c) = b2(1)*half(3*l - 2, c) + b2(2)*half(3*l - 1, c) &
                    + b2(3)*half(3*l, c)
            end do
            do m = 1, nmode
                do jp = 1, map%npoly
                    do i = 1, 2
                        r = 2*nlm + i + 2*(jp - 1) + 2*map%npoly*(m - 1)
                        s = 0
                        do a = 1, nloc
                            l = m + nmode*(a - 1)
                            node = map%conn(a, e)
                            if (nested) then
                                av = map%axis_vec(:, i + 2*(jp - 1), node)
                                s = s + av(1)*half(3*l - 2, c) + av(2)*half(3*l - 1, c)
                            else
                                s = s + map%axis_basis(jp, node)*half(3*l - 3 + i, c)
                            end if
                        end do
                        blk(r, c) = s
                    end do
                end do
            end do
        end do
        ! The scalar-energy Hessian has one authoritative triangle.
        ! Reflect after P^T H P: cancellation in independently evaluated
        ! contractions must not create an antisymmetric Newton operator.
        do c = 1, nr
            do r = 1, c - 1
                blk(r, c) = blk(c, r)
            end do
        end do
    end subroutine

    subroutine ideal_map_solve(map, tolerance, max_iterations, report, status, verbose, &
                               reference_gradient)
        ! Damped Newton on the reduced (nu, omega) coefficients: Levenberg-Marquardt
        ! with the Hessian diagonal as metric, shrinking the damping after every
        ! accepted step. Steps keep J > 0 and do not raise W beyond roundoff.
        ! Converged when the free reduced-gradient norm is below tolerance times
        ! its initial value or the step stalls at roundoff.
        ! With map%omega_elimination the iteration is Newton on W_red(nu) =
        ! min_omega W(nu, omega): omega is eliminated exactly at the start and
        ! at every trial point, so its gradient vanishes and the full-system
        ! solve with right-hand side -g is the Schur-complement step
        ! (H_tt - H_to H_oo^-1 H_ot) d nu = -g_t; damping acts on the
        ! transversal block only and the trust region on nu only.
        type(ideal_map_t), intent(inout) :: map
        real(dp), intent(in) :: tolerance
        integer, intent(in) :: max_iterations
        type(ideal_map_solve_report_t), intent(out) :: report
        integer, intent(out) :: status
        logical, intent(in), optional :: verbose
        ! Absolute gradient the tolerance refers to (default: the initial one),
        ! e.g. the first-phase gradient for a polishing solve.
        real(dp), intent(in), optional :: reference_gradient
        type(csr_matrix_t) :: matrix, omatrix
        type(sparse_solver_t) :: solver, osolver
        type(fortsparse_status_t) :: sstatus
        real(dp), allocatable :: grad(:), gred(:), step(:), full(:), diag(:)
        real(dp), allocatable :: trial(:, :, :), rhs(:), shift(:)
        logical, allocatable :: omask(:), damped(:)
        real(dp) :: energy, trial_energy, g0, gnorm, mu, minjac, scale, alpha
        real(dp) :: cap_x, cap_w, t0, trial_minjac, pf1_t0
        character(8) :: int5_env
        real(dp) :: eta, eta_prev, gnorm_prev
        real(dp), allocatable :: krylov_rhs(:), krylov_mask(:)
        real(dp), parameter :: eta_max = 0.1_dp, eta_min = 1.0e-4_dp
        real(dp), parameter :: eta_gamma = 0.9_dp
        integer :: it, i, k, tries, estatus, kinfo, kiters
        logical :: talk, accepted, factored, ofactored, nested, refresh, pc_failed
        logical :: retained
        talk = .false.
        ! INT5 measurement switch (never for a PR): INT5_KRYLOV=1 enables the CG step.
        call get_environment_variable('INT5_KRYLOV', int5_env)
        if (trim(int5_env) == '1') map%krylov_step = .true.
        factored = .false.; ofactored = .false.; refresh = .false.
        ! factored: solver holds the symbolic analysis, which factor_masked
        ! keeps after a pivot that is not positive; retained: it also holds a
        ! valid numerical factor, the CG preconditioner.
        retained = .false.
        pc_failed = .false.
        eta = eta_max; eta_prev = eta_max; gnorm_prev = 0
        if (present(verbose)) talk = verbose
        nested = map%omega_elimination .and. .not. map%position_gauge
        allocate (grad(map%ndof), full(map%ndof), gred(map%nred), step(map%nred), &
            diag(map%nred), rhs(map%nred), shift(map%nred), krylov_rhs(map%nred), &
            krylov_mask(map%nred))
        call omega_mask(map, nested, omask)
        damped = .not. (map%fixed .or. omask)
        t0 = wall_seconds()
        call ideal_map_energy(map, map%x, energy, estatus, grad, minjac)
        report%energy_s = report%energy_s + wall_seconds() - t0
        if (estatus /= IDEAL_MAP_OK) then
            status = IDEAL_MAP_FOLDED; return
        end if
        call ideal_map_reduce(map, grad, gred)
        where (map%fixed) gred = 0
        g0 = norm2(gred); gnorm = g0
        report%initial_gradient = g0
        if (present(reference_gradient)) g0 = reference_gradient
        if (nested) then
            call eliminate_omega(map, map%x, omask, omatrix, osolver, ofactored, &
                energy, gred, minjac, report, estatus)
            if (estatus /= IDEAL_MAP_OK) then
                status = estatus; call cleanup(); return
            end if
            gnorm = norm2(gred)
        end if
        ! A failed trial multiplies mu by 10, so mu stays positive (the floor of
        ! the accept rule below).
        mu = max(map%initial_damping, 1.0e-14_dp)
        cap_x = 0.05_dp*maxval(maxval(map%reference, dim=2) &
            - minval(map%reference, dim=2))
        cap_w = 1
        if (nested) cap_w = huge(1.0_dp)
        alpha = 1
        status = IDEAL_MAP_NOT_CONVERGED
        report%energy = energy; report%gradient = gnorm
        do it = 1, max_iterations
            if (.not. gnorm > tolerance*g0) then
                status = IDEAL_MAP_OK; exit
            end if
            report%iterations = it
            t0 = wall_seconds()
            call ideal_map_hessian(map, map%x, matrix)
            report%hessian_s = report%hessian_s + wall_seconds() - t0
            do i = 1, map%nred
                do k = matrix%rowptr(i), matrix%rowptr(i + 1) - 1
                    if (matrix%colidx(k) == i) diag(i) = matrix%val(k)
                end do
            end do
            scale = maxval(abs(diag), mask=.not. map%fixed)
            accepted = .false.
            if (nested) then
                call nested_line_search()
                if (status == IDEAL_MAP_SOLVE_FAILED) then
                    call cleanup(); return
                end if
            end if
            if (.not. nested) then
            ! Eisenstat-Walker forcing (choice 2, safeguarded) for the CG step.
            if (it > 1) then
                eta = eta_gamma*(gnorm/gnorm_prev)**2
                if (eta_gamma*eta_prev**2 > 0.1_dp) &
                    eta = max(eta, eta_gamma*eta_prev**2)
            end if
            eta = min(max(eta, eta_min), eta_max)
            do tries = 1, 30
                t0 = wall_seconds()
                shift = merge(mu*(abs(diag) + 1.0e-12_dp*scale), 0.0_dp, damped)
                report%trials = report%trials + 1
                rhs = -gred
                kinfo = -1
                if (map%krylov_step .and. retained .and. .not. refresh) then
                    call krylov_solve()
                    if (kiters > map%krylov_refactor) refresh = .true.
                end if
                if (kinfo /= KRYLOV_OK) then
                    call factor_masked(map, matrix, .not. map%fixed, shift, solver, &
                        factored, estatus, report)
                    retained = estatus == 0
                    if (estatus /= 0) then
                        mu = mu*10; cycle
                    end if
                    refresh = .false.
                    pf1_t0 = wall_seconds()
                    call sparse_solve(solver, rhs, step, sstatus)
                    pf1_t(11) = pf1_t(11) + wall_seconds() - pf1_t0; pf1_n(11) = pf1_n(11) + 1
                    if (.not. status_ok(sstatus)) then
                        status = IDEAL_MAP_SOLVE_FAILED; call cleanup(); return
                    end if
                end if
                report%factor_s = report%factor_s + wall_seconds() - t0
                where (map%fixed) step = 0
                call ideal_map_expand(map, step, full)
                ! Trust region: near-null (relabelling) directions can carry
                ! large Newton components; cap the displacement and rotation.
                call trust_region()
                full = alpha*full; step = alpha*step
                trial = map%x + reshape(full, shape(map%x))
                t0 = wall_seconds()
                call ideal_map_energy(map, trial, trial_energy, estatus)
                report%energy_s = report%energy_s + wall_seconds() - t0
                if (estatus == IDEAL_MAP_OK .and. trial_energy < energy &
                    - 64*epsilon(1.0_dp)*abs(energy)) then
                    accepted = .true.; exit
                end if
                ! Within the roundoff band of W the energy cannot rank steps.
                ! Rank them by the trapezoidal energy change
                ! (g(x) + g(x + s)) . s / 2 (error O(|s|^3), no cancellation
                ! against W); near-null directions need full Newton steps that
                ! do not decrease |g| monotonically.
                if (estatus == IDEAL_MAP_OK .and. trial_energy <= energy &
                    + 1024*epsilon(1.0_dp)*abs(energy)) then
                    call ideal_map_energy(map, trial, trial_energy, estatus, grad)
                    call ideal_map_reduce(map, grad, rhs)
                    where (map%fixed) rhs = 0
                    if (dot_product(gred + rhs, step) <= 0) then
                        accepted = .true.; exit
                    end if
                end if
                mu = mu*10
            end do
            end if
            if (.not. accepted) exit
            report%step = maxval(abs(trial - map%x))
            map%x = trial
            if (nested) then
                energy = trial_energy; gred = rhs; minjac = trial_minjac
            else
                t0 = wall_seconds()
                call ideal_map_energy(map, map%x, energy, estatus, grad, minjac)
                report%energy_s = report%energy_s + wall_seconds() - t0
                call ideal_map_reduce(map, grad, gred)
                where (map%fixed) gred = 0
            end if
            gnorm_prev = gnorm; eta_prev = eta
            gnorm = norm2(gred)
            report%energy = energy; report%gradient = gnorm/g0
            report%min_jacobian = minjac; report%damping = mu
            if (talk) print '(a,i3,5es12.4,i3,3es10.2)', 'ideal map newton', it, &
                energy, gnorm/g0, maxval(abs(step(1::2))), maxval(abs(step(2::2))), mu, &
                tries, alpha, norm2(gred(2:2*map%nodes:2))/gnorm, &
                norm2(step(2:2*map%nodes:2))/norm2(step)
            if (.not. alpha < 1) mu = max(mu/100, 1.0e-14_dp)
            if (report%step <= 1.0e-13_dp*maxval(abs(map%reference))) then
                status = IDEAL_MAP_OK; exit
            end if
        end do
        if (.not. gnorm > tolerance*g0) status = IDEAL_MAP_OK
        call cleanup()
    contains
        subroutine trust_region()
            ! alpha = largest fraction of full within the displacement caps.
            alpha = 1
            do i = 1, map%nodes
                k = ideal_map_dof(map, 1, 1, i)
                alpha = min(alpha, cap_x/max(tiny(1.0_dp), &
                    maxval(abs(full(k:k + 3*map%nmode - 1:3))), &
                    maxval(abs(full(k + 1:k + 3*map%nmode - 1:3)))), &
                    cap_w/max(tiny(1.0_dp), maxval(abs(full(k + 2:k + 3*map%nmode &
                    - 1:3)))))
            end do
        end subroutine

        subroutine krylov_solve()
            ! step = -(masked H + shift)^-1 gred by CG to the relative residual
            ! eta, preconditioned with the retained factorisation (identity on
            ! the fixed rows, as in factor_masked). kinfo = KRYLOV_OK only for
            ! a converged solve; otherwise the caller refactors. The system is
            ! scaled to a unit right-hand side: the fortnum tolerance is
            ! absolute below unit norm.
            real(dp) :: bnorm, residual
            kiters = 0; pc_failed = .false.
            bnorm = norm2(rhs)
            if (.not. bnorm > 0) return
            krylov_rhs = rhs/bnorm
            step = 0
            call real_conjugate_gradient_operator(apply_damped, krylov_rhs, step, &
                eta, map%krylov_max, kinfo, kiters, residual, apply_factor)
            report%krylov_solves = report%krylov_solves + 1
            report%krylov_iterations = report%krylov_iterations + kiters
            if (pc_failed) kinfo = -1
            if (kinfo == KRYLOV_OK) step = bnorm*step
        end subroutine

        subroutine apply_damped(input, output)
            real(dp), intent(in) :: input(:)
            real(dp), intent(out) :: output(:)
            real(dp) :: s
            integer :: row, kk
            krylov_mask = merge(0.0_dp, input, map%fixed)
            !$omp parallel do private(row, kk, s) schedule(static)
            do row = 1, map%nred
                if (map%fixed(row)) then
                    output(row) = input(row)
                    cycle
                end if
                s = shift(row)*input(row)
                do kk = matrix%rowptr(row), matrix%rowptr(row + 1) - 1
                    s = s + matrix%val(kk)*krylov_mask(matrix%colidx(kk))
                end do
                output(row) = s
            end do
            !$omp end parallel do
        end subroutine

        subroutine apply_factor(input, output)
            real(dp), intent(in) :: input(:)
            real(dp), intent(out) :: output(:)
            type(fortsparse_status_t) :: pstatus
            call sparse_solve(solver, input, output, pstatus)
            if (.not. status_ok(pstatus)) then
                pc_failed = .true.; output = input
            end if
        end subroutine

        subroutine nested_line_search()
            ! Nested step: the damped Schur-complement Newton direction d of
            ! W_red, then a line search along it (W_red is strongly nonlinear
            ! along near-null shaping directions, where full steps zigzag).
            ! d must be a descent direction; g.d >= 0 means negative
            ! curvature of the damped Schur complement along d and raises the
            ! damping. The step length t (from the trust-region cap) is
            ! accepted when W_red(t) - W_red(0) <= 1e-4 t g.d (Armijo) and
            ! g(t).d <= |g.d|/2 (no overshoot past the line minimum); otherwise
            ! t moves to the minimiser of the cubic through W_red and g.d at 0
            ! and t, kept in [0.1 t, 0.9 t]. Within the roundoff band of W the
            ! change is the trapezoidal (g.d + g(t).d) t/2.
            real(dp) :: gd, g1d, t, dw, d1, d2
            integer :: ls
            status = IDEAL_MAP_NOT_CONVERGED
            do tries = 1, 30
                t0 = wall_seconds()
                shift = merge(mu*(abs(diag) + 1.0e-12_dp*scale), 0.0_dp, damped)
                call factor_masked(map, matrix, .not. map%fixed, shift, solver, &
                    factored, estatus, report, gauge=.true.)
                if (estatus /= 0) then
                    mu = mu*10; cycle
                end if
                rhs = -gred
                call gauge_solve(map, solver, rhs, 0*gauge_moments(map, map%x), step, &
                    sstatus)
                report%factor_s = report%factor_s + wall_seconds() - t0
                if (.not. status_ok(sstatus)) then
                    status = IDEAL_MAP_SOLVE_FAILED; return
                end if
                ! omega is re-eliminated at every trial point.
                where (map%fixed .or. omask) step = 0
                gd = dot_product(gred, step)
                if (.not. gd < 0) then
                    mu = mu*10; cycle
                end if
                call ideal_map_expand(map, step, full)
                call trust_region()
                t = alpha
                do ls = 1, 30
                    trial = map%x + t*reshape(full, shape(map%x))
                    call eliminate_omega(map, trial, omask, omatrix, osolver, &
                        ofactored, trial_energy, rhs, trial_minjac, report, estatus)
                    if (estatus == IDEAL_MAP_SOLVE_FAILED) then
                        status = estatus; return
                    end if
                    if (estatus /= IDEAL_MAP_OK) then
                        t = t/4; cycle
                    end if
                    g1d = dot_product(rhs, step)
                    dw = trial_energy - energy
                    if (abs(dw) <= 64*epsilon(1.0_dp)*abs(energy)) dw = (gd + g1d)*t/2
                    if (dw <= 1.0e-4_dp*t*gd .and. g1d <= abs(gd)/2) then
                        accepted = .true.; alpha = t
                        step = t*step
                        return
                    end if
                    ! Cubic in s in [0, t]: f(0) = 0, f'(0) = gd, f(t) = dw,
                    ! f'(t) = g1d.
                    d1 = gd + g1d - 3*dw/t
                    d2 = d1*d1 - gd*g1d
                    if (d2 > 0) then
                        d2 = sqrt(d2)
                        d1 = t - t*(g1d + d2 - d1)/(g1d - gd + 2*d2)
                    else
                        d1 = t/2
                    end if
                    if (.not. ieee_is_finite(d1)) d1 = t/2
                    t = min(0.9_dp*t, max(0.1_dp*t, d1))
                end do
                mu = mu*10
            end do
        end subroutine

        subroutine cleanup()
            if (factored) call sparse_destroy(solver)
            if (ofactored) call sparse_destroy(osolver)
            factored = .false.; ofactored = .false.
        end subroutine
    end subroutine

    subroutine omega_mask(map, nested, omask)
        ! Free omega dofs (all modes) when nested; none otherwise.
        type(ideal_map_t), intent(in) :: map
        logical, intent(in) :: nested
        logical, allocatable, intent(out) :: omask(:)
        integer :: a, m
        allocate (omask(map%nred))
        omask = .false.
        if (.not. nested) return
        do a = 1, map%nodes
            do m = 1, map%nmode
                omask(ideal_map_rdof(map, 2, m, a)) = .true.
            end do
        end do
        omask = omask .and. .not. map%fixed
    end subroutine

    subroutine ideal_map_eliminate_omega(map, x, energy, gradient, status)
        ! Reduced functional W_red(nu) = min_omega W(nu, omega) at the
        ! transversal dofs of x: x(3, mode, node) is returned with omega
        ! replaced by the minimiser omega*(nu) (fixed dofs kept), energy =
        ! W_red, gradient(nred) the reduced gradient there (fixed dofs zero;
        ! its omega part vanishes, the rest is grad W_red by the envelope
        ! property). Transversal gauge only; the gauge form is not used when
        ! map%omega_elimination is set before ideal_map_freeze.
        type(ideal_map_t), intent(in) :: map
        real(dp), intent(inout) :: x(:, :, :)
        real(dp), intent(out) :: energy, gradient(:)
        integer, intent(out) :: status
        type(csr_matrix_t) :: matrix
        type(sparse_solver_t) :: slv
        type(ideal_map_solve_report_t) :: report
        logical, allocatable :: omask(:)
        logical :: done
        real(dp) :: minjac
        done = .false.
        if (map%position_gauge) then
            status = IDEAL_MAP_BAD_INPUT; return
        end if
        call omega_mask(map, .true., omask)
        call eliminate_omega(map, x, omask, matrix, slv, done, energy, gradient, &
            minjac, report, status)
        if (done) call sparse_destroy(slv)
    end subroutine

    subroutine eliminate_omega(map, x, omask, matrix, slv, done, energy, gred, &
                               minjac, report, info)
        ! Minimise W over the free omega dofs (omask) at fixed transversal
        ! dofs: one undamped solve H_oo d omega = -g_o (the pin removes the
        ! exact null). At ntor = 0 with the cylindrical chart J does not depend
        ! on omega and W is exactly quadratic in it, so the solve is exact.
        ! ntor > 0 is rejected (IDEAL_MAP_BAD_INPUT): there J depends on
        ! omega, and an undamped inner Newton probe (P2, 16 segments, nfp 3)
        ! converged only linearly (|g_o| halving per pass).
        type(ideal_map_t), intent(in) :: map
        real(dp), intent(inout) :: x(:, :, :)
        logical, intent(in) :: omask(:)
        type(csr_matrix_t), intent(inout) :: matrix
        type(sparse_solver_t), intent(inout) :: slv
        logical, intent(inout) :: done
        real(dp), intent(out) :: energy, gred(:), minjac
        type(ideal_map_solve_report_t), intent(inout) :: report
        integer, intent(out) :: info
        type(fortsparse_status_t) :: sstatus
        real(dp), allocatable :: grad(:), rhs(:), dw(:), full(:), shift(:)
        real(dp) :: t0
        if (map%ntor /= 0) then
            info = IDEAL_MAP_BAD_INPUT; return
        end if
        ! The one-solve quadratic elimination is specific to an axisymmetric
        ! cylindrical chart, not to a zero Fourier mode in an arbitrary chart.
        select type (chart => map%chart)
        class is (cylindrical_chart_t)
        class default
            info = IDEAL_MAP_BAD_INPUT; return
        end select
        allocate (grad(map%ndof), full(map%ndof), rhs(map%nred), dw(map%nred), &
            shift(map%nred))
        shift = 0
        call evaluate()
        if (info /= IDEAL_MAP_OK) return
        t0 = wall_seconds()
        call ideal_map_hessian(map, x, matrix)
        report%hessian_s = report%hessian_s + wall_seconds() - t0
        t0 = wall_seconds()
        call factor_masked(map, matrix, omask, shift, slv, done, info, report, &
            gauge=.true.)
        if (info /= 0) then
            info = IDEAL_MAP_SOLVE_FAILED; return
        end if
        rhs = merge(-gred, 0.0_dp, omask)
        ! The step also restores the gauge moments to zero (after prolong).
        call gauge_solve(map, slv, rhs, -gauge_moments(map, x), dw, sstatus)
        report%factor_s = report%factor_s + wall_seconds() - t0
        if (.not. status_ok(sstatus)) then
            info = IDEAL_MAP_SOLVE_FAILED; return
        end if
        dw = merge(dw, 0.0_dp, omask)
        call ideal_map_expand(map, dw, full)
        x = x + reshape(full, shape(x))
        call evaluate()
    contains
        subroutine evaluate()
            t0 = wall_seconds()
            call ideal_map_energy(map, x, energy, info, grad, minjac)
            report%energy_s = report%energy_s + wall_seconds() - t0
            call ideal_map_reduce(map, grad, gred)
            where (map%fixed) gred = 0
            call ideal_map_gauge_project(map, gred)
        end subroutine
    end subroutine

    pure integer function gauge_rows(map) result(nc)
        ! Number of omega gauge constraints (Lagrange rows after the nred).
        type(ideal_map_t), intent(in) :: map
        nc = 0
        if (allocated(map%omega_c)) nc = size(map%omega_c, 1)
    end function

    subroutine gauge_solve(map, slv, b, rc, x, sstatus)
        ! Solve with a factorisation from factor_masked(gauge = .true.):
        ! right-hand side b (nred) and constraint values rc, x(nred) returned.
        type(ideal_map_t), intent(in) :: map
        type(sparse_solver_t), intent(inout) :: slv
        real(dp), intent(in) :: b(:), rc(:)
        real(dp), intent(out) :: x(:)
        type(fortsparse_status_t), intent(out) :: sstatus
        real(dp), allocatable :: bb(:), xx(:)
        allocate (bb(map%nred + gauge_rows(map)), xx(map%nred + gauge_rows(map)))
        bb(:map%nred) = b; bb(map%nred + 1:) = rc
        call sparse_solve(slv, bb, xx, sstatus)
        x = xx(:map%nred)
    end subroutine

    function gauge_moments(map, x) result(rc)
        ! Constraint values sum_a phi_l(rho_a) omega_a of the coefficients x.
        type(ideal_map_t), intent(in) :: map
        real(dp), intent(in) :: x(:, :, :)
        real(dp), allocatable :: rc(:)
        allocate (rc(gauge_rows(map)))
        if (size(rc) == 0) return
        rc = matmul(map%omega_c, x(3, 1, :))
    end function

    subroutine factor_masked(map, matrix, keep, shift, slv, done, info, report, gauge)
        ! Factorise the reduced Hessian restricted to the keep rows/columns
        ! (identity elsewhere), plus shift on the diagonal. done: slv holds a
        ! factorisation of this pattern (only the numerical one is redone).
        ! gauge: append the omega gauge constraints (nested_omega_gauge) on
        ! the free n = 0 omega dofs as Lagrange rows/columns (saddle point;
        ! solve with gauge_solve). The masked, symmetrised matrix is cached in
        ! report per pattern; a trial only rewrites the diagonal.
        type(ideal_map_t), intent(in) :: map
        type(csr_matrix_t), intent(in) :: matrix
        logical, intent(in) :: keep(:)
        real(dp), intent(in) :: shift(:)
        type(sparse_solver_t), intent(inout) :: slv
        logical, intent(inout) :: done
        integer, intent(out) :: info
        type(ideal_map_solve_report_t), intent(inout) :: report
        logical, intent(in), optional :: gauge
        type(sparse_factor_stats_t) :: stats
        type(fortsparse_status_t) :: sstatus
        logical :: reused, spd
        integer :: row, n, nc
        real(dp) :: dscale, mu, pf1_t0
        pf1_t0 = wall_seconds()
        nc = 0
        if (present(gauge)) then
            if (gauge) nc = gauge_rows(map)
        end if
        n = map%nred + nc
        call factor_cache_update(map, matrix, keep, nc, report%cache, info)
        if (info /= 0) then
            if (done) call sparse_destroy(slv)
            done = .false.
            return
        end if
        associate (c => report%cache)
            report%assembly_roundoff = max(report%assembly_roundoff, c%correction)
            do row = 1, map%nred
                if (keep(row) .and. c%dpos(row) > 0) &
                    c%a%val(c%dpos(row)) = c%bdiag(row) + shift(row)
            end do
            ! The damped Hessian is symmetric and positive definite near the
            ! solution, so try the SPD Cholesky first: a trial that is not PD
            ! fails cheaply and the caller raises the damping. The relative
            ! damping mu is recovered from the shift mu*(|d| + 1e-12 max|d|).
            ! The gauge saddle point is indefinite, and at a damping above
            ! spd_mu_limit a failed Cholesky is no longer informative, so those
            ! use the native supernodal pivoted LU action (also used for the
            ! free-boundary KKT): scalar LDLT Schur-row updates are too costly
            ! on the map ladders. This changes only the factor, not the Hessian.
            spd = nc == 0
            if (spd) then
                dscale = maxval(abs(c%bdiag), mask=keep(:map%nred))
                mu = 0
                do row = 1, map%nred
                    if (keep(row) .and. shift(row) > 0) mu = max(mu, shift(row)/ &
                        (abs(c%bdiag(row)) + 1.0e-12_dp*dscale))
                end do
                spd = .not. mu > spd_mu_limit
            end if
            if (spd) then
                slv%options%symmetry = FORTSPARSE_SYM_SPD
                slv%options%lu_method = FORTSPARSE_LU_SCALAR
            else
                slv%options%symmetry = FORTSPARSE_SYM_GENERAL
                slv%options%lu_method = FORTSPARSE_LU_SUPERNODAL
            end if
            slv%options%ordering = FORTSPARSE_ORDER_AMD
            slv%options%verify_factor = map%verify_factor
            pf1_t(8) = pf1_t(8) + wall_seconds() - pf1_t0; pf1_n(8) = pf1_n(8) + 1
            pf1_t0 = wall_seconds()
            if (gs_direct_available()) then
                if (.not. gs_direct_mpi_ready()) then
                    info = 1
                    if (done) call sparse_destroy(slv)
                    done = .false.
                    return
                end if
                slv%backend_id = FORTSPARSE_BACKEND_NATIVE
                slv%options%nthreads = map%factor_threads
            else
                slv%backend_id = FORTSPARSE_BACKEND_NATIVE
            end if
            ! The pattern is fixed per map: after the first factorisation only
            ! the numerical factorisation is redone when the backend can reuse
            ! it.
            reused = .false.
            if (done) then
                call sparse_refactor(slv, c%a, sstatus, reused)
                pf1_t(10) = pf1_t(10) + wall_seconds() - pf1_t0
                pf1_n(10) = pf1_n(10) + 1
                if (reused) pf1_n(12) = pf1_n(12) + 1
            else
                call sparse_factor(slv, c%a, sstatus)
                pf1_t(9) = pf1_t(9) + wall_seconds() - pf1_t0; pf1_n(9) = pf1_n(9) + 1
            end if
        end associate
        info = merge(0, 1, status_ok(sstatus))
        report%factorizations = report%factorizations + 1
        if (info /= 0) then
            ! A pivot that is not positive keeps the symbolic analysis: the
            ! next trial refactors on it without a new analysis.
            if (spd .and. sstatus%code == FORTSPARSE_NOT_POSITIVE_DEFINITE) then
                if (.not. done) report%analyses = report%analyses + 1
                done = .true.
                return
            end if
            if (done) call sparse_destroy(slv)
            done = .false.
            return
        end if
        done = .true.
        if (.not. reused) report%analyses = report%analyses + 1
        call sparse_factor_stats(slv, stats, sstatus)
        pf1_stats = stats
        if (status_ok(sstatus)) then
            if (.not. reused) report%analysis_s = report%analysis_s &
                + max(0.0_dp, stats%analysis_s)
            report%numeric_s = report%numeric_s + max(0.0_dp, stats%factor_s)
        end if
    end subroutine

    subroutine factor_cache_update(map, matrix, keep, nc, c, info)
        ! Bring the cache up to date with the Hessian: rebuild the pattern data
        ! when the pattern, the mask or the constraints changed, and the base
        ! matrix when the Hessian values changed (once per Newton iteration,
        ! not per trial). info /= 0: invalid pattern or asymmetry above 1e-12.
        type(ideal_map_t), intent(in) :: map
        type(csr_matrix_t), intent(in) :: matrix
        logical, intent(in) :: keep(:)
        integer, intent(in) :: nc
        type(factor_cache_t), intent(inout) :: c
        integer, intent(out) :: info
        logical :: same
        same = c%nc == nc .and. allocated(c%in_colidx)
        if (same) same = size(c%in_colidx) == size(matrix%colidx) .and. &
            size(c%in_rowptr) == size(matrix%rowptr) .and. &
            all(c%keep .eqv. keep(:map%nred))
        if (same) same = all(c%in_rowptr == matrix%rowptr) .and. &
            all(c%in_colidx == matrix%colidx)
        if (.not. same) then
            call factor_cache_build(map, matrix, keep, nc, c)
            c%base_ok = .false.
        end if
        info = 1
        if (.not. c%pattern_ok) return
        if (c%base_ok) c%base_ok = all(abs(c%raw - matrix%val) <= 0.0_dp)
        if (.not. c%base_ok) call factor_cache_refresh(map, matrix, c, info)
        if (c%base_ok) info = 0
    end subroutine

    subroutine factor_cache_build(map, matrix, keep, nc, c)
        ! Factor input pattern: the Hessian pattern, extended by the gauge
        ! constraint rows/columns (nc > 0), with the data of the per-trial and
        ! per-iteration updates. Every nonzero moment is kept, including tiny
        ! signed P_k weights: this is sparse support, not an approximation.
        type(ideal_map_t), intent(in) :: map
        type(csr_matrix_t), intent(in) :: matrix
        logical, intent(in) :: keep(:)
        integer, intent(in) :: nc
        type(factor_cache_t), intent(inout) :: c
        logical, allocatable :: omega_active(:, :)
        integer, allocatable :: node(:)
        integer :: row, kk, a, l, n, k, nz, nm, ng, previous, col, p
        c%nc = nc
        c%pattern_ok = .false.
        c%keep = keep(:map%nred)
        c%in_rowptr = matrix%rowptr
        c%in_colidx = matrix%colidx
        n = map%nred + nc
        nz = size(matrix%val)
        if (nc == 0) then
            c%rowptr = matrix%rowptr
            c%colidx = matrix%colidx
            c%src = [(k, k=1, nz)]
            c%gpos = [integer ::]
            c%gval = [real(dp) ::]
        else
            ! Constrained dofs: free n = 0 omega (kept in every gauge use).
            allocate (node(map%nred))
            node = 0
            do a = 1, map%nodes
                row = ideal_map_rdof(map, 2, 1, a)
                if (.not. map%fixed(row)) node(row) = a
            end do
            omega_active = .not. (abs(map%omega_c) <= 0.0_dp)
            ng = count(omega_active)
            c%rowptr = [(0, k=1, n + 1)]
            c%colidx = [(0, k=1, nz + 2*ng + nc)]
            c%src = c%colidx
            c%gpos = [(0, k=1, 2*ng + nc)]
            c%gval = [(0.0_dp, k=1, 2*ng + nc)]
            k = 0; nm = 0; c%rowptr(1) = 1
            do row = 1, map%nred
                do kk = matrix%rowptr(row), matrix%rowptr(row + 1) - 1
                    k = k + 1; c%colidx(k) = matrix%colidx(kk); c%src(k) = kk
                end do
                a = node(row)
                if (a > 0) then
                    do l = 1, nc
                        if (.not. omega_active(l, a)) cycle
                        k = k + 1; c%colidx(k) = map%nred + l; c%src(k) = 0
                        nm = nm + 1; c%gpos(nm) = k; c%gval(nm) = map%omega_c(l, a)
                    end do
                end if
                c%rowptr(row + 1) = k + 1
            end do
            do l = 1, nc
                do row = 1, map%nred
                    a = node(row)
                    if (a == 0) cycle
                    if (.not. omega_active(l, a)) cycle
                    k = k + 1; c%colidx(k) = row; c%src(k) = 0
                    nm = nm + 1; c%gpos(nm) = k; c%gval(nm) = map%omega_c(l, a)
                end do
                k = k + 1; c%colidx(k) = map%nred + l; c%src(k) = 0
                nm = nm + 1; c%gpos(nm) = k; c%gval(nm) = 0
                c%rowptr(map%nred + l + 1) = k + 1
            end do
            c%colidx = c%colidx(:k); c%src = c%src(:k)
            c%gpos = c%gpos(:nm); c%gval = c%gval(:nm)
        end if
        nz = size(c%colidx)
        ! Diagonal positions and the masked Hessian entries.
        c%dpos = [(0, k=1, map%nred)]
        nm = 0
        do row = 1, map%nred
            do k = c%rowptr(row), c%rowptr(row + 1) - 1
                if (c%src(k) == 0) cycle
                if (c%colidx(k) == row) c%dpos(row) = k
                if (.not. keep(row) .or. .not. keep(c%colidx(k))) nm = nm + 1
            end do
        end do
        c%mpos = [(0, k=1, nm)]
        c%mval = [(0.0_dp, k=1, nm)]
        nm = 0
        do row = 1, map%nred
            do k = c%rowptr(row), c%rowptr(row + 1) - 1
                if (c%src(k) == 0) cycle
                col = c%colidx(k)
                if (keep(row) .and. keep(col)) cycle
                nm = nm + 1; c%mpos(nm) = k
                c%mval(nm) = merge(1.0_dp, 0.0_dp, col == row .and. .not. keep(row))
            end do
        end do
        ! A sorted, unique and structurally symmetric full pattern is required
        ! by the reflection of the upper triangle onto the lower one.
        c%refl = [(0, k=1, nz)]
        do row = 1, n
            previous = 0
            do k = c%rowptr(row), c%rowptr(row + 1) - 1
                col = c%colidx(k)
                if (col <= previous .or. col > n) return
                previous = col
                if (col == row) cycle
                p = find_entry(c, col, row)
                if (p == 0) return
                if (col > row) c%refl(k) = p
            end do
        end do
        c%raw = [(0.0_dp, k=1, size(matrix%val))]
        c%base = [(0.0_dp, k=1, nz)]
        c%bdiag = [(0.0_dp, k=1, map%nred)]
        c%a%nrow = n; c%a%ncol = n; c%a%nnz = nz
        c%a%col_ptr = c%rowptr; c%a%row_idx = c%colidx; c%a%val = c%base
        c%pattern_ok = .true.
    end subroutine

    integer function find_entry(c, row, col) result(found)
        ! Position of (row, col) in the sorted factor pattern, 0 if absent.
        type(factor_cache_t), intent(in) :: c
        integer, intent(in) :: row, col
        integer :: lo, hi, mid
        found = 0
        lo = c%rowptr(row); hi = c%rowptr(row + 1) - 1
        do while (lo <= hi)
            mid = lo + (hi - lo)/2
            if (c%colidx(mid) == col) then
                found = mid
                return
            else if (c%colidx(mid) < col) then
                lo = mid + 1
            else
                hi = mid - 1
            end if
        end do
    end function

    subroutine factor_cache_refresh(map, matrix, c, info)
        ! base = masked Hessian (+ gauge rows) on one canonical triangle: the
        ! strict upper entries are copied onto the lower ones, which preserves
        ! the scalar Hessian's exact symmetry after independently rounded
        ! cell/constraint accumulation. The asymmetry removed (outward one-ulp
        ! bound) must not exceed 1e-12 of the largest entry.
        type(ideal_map_t), intent(in) :: map
        type(csr_matrix_t), intent(in) :: matrix
        type(factor_cache_t), intent(inout) :: c
        integer, intent(out) :: info
        real(dp), parameter :: rtol = 1.0e-12_dp
        real(dp), allocatable :: reference(:)
        type(fortsparse_status_t) :: sstatus
        integer :: k, p, row
        real(dp) :: delta
        info = 1
        c%correction = 0
        if (c%nc == 0) then
            c%base = matrix%val
        else
            do k = 1, size(c%base)
                if (c%src(k) > 0) c%base(k) = matrix%val(c%src(k))
            end do
        end if
        c%base(c%gpos) = c%gval
        c%base(c%mpos) = c%mval
        if (.not. all(ieee_is_finite(c%base))) return
        c%scale = maxval(abs(c%base))
        do k = 1, size(c%base)
            p = c%refl(k)
            if (p == 0) cycle
            delta = abs(c%base(k) - c%base(p))
            if (delta <= 0.0_dp) cycle
            c%correction = max(c%correction, nearest(delta, 1.0_dp))
        end do
        if (c%correction > rtol*c%scale) return
        do k = 1, size(c%base)
            if (c%refl(k) > 0) c%base(c%refl(k)) = c%base(k)
        end do
        c%a%val = c%base
        c%raw = matrix%val
        do row = 1, map%nred
            if (c%dpos(row) > 0) c%bdiag(row) = c%base(c%dpos(row))
        end do
        if (map%verify_factor) then
            ! Debug: FortSparse's canonical projection must find nothing to do.
            reference = c%base
            call csc_normalize_symmetric_raw(size(c%rowptr) - 1, c%rowptr, c%colidx, &
                reference, rtol, delta, sstatus)
            if (.not. status_ok(sstatus) .or. delta > 0.0_dp .or. &
                any(abs(reference - c%base) > 0.0_dp)) &
                error stop 'ideal_map: cached factor input is not canonical'
        end if
        c%base_ok = .true.
        info = 0
    end subroutine

    real(dp) function wall_seconds()
        integer(int64) :: count, rate
        call system_clock(count, rate)
        wall_seconds = real(count, dp)/rate
    end function

    subroutine ideal_map_residual(map, f_vac, residual, sigma, status)
        ! Weak Maxwell-stress residual (equilibrium-core D3) of the current map,
        ! normalised to the poloidal-field scale. For the map lane R[xi] = -dW[xi],
        ! so the coefficient gradient g holds R tested with every basis
        ! displacement xi = (dR, dZ, R d omega) (boundary dR, dZ held). The dual
        ! norm uses the consistent L2(V) mass of xi: f = M^-1 g is the force
        ! density and
        !     residual = ||f||_L2 L / (sigma |V|^1/2),  L = (A_ref/pi)^1/2,
        !     sigma = <|B - B_0|^2>/mu0 + <p>,  B_0 = F_vac/R e_phi,
        ! averages over the plasma volume, so a tokamak is measured against
        ! B_p^2 and pressure, never against B_phi^2.
        type(ideal_map_t), intent(in) :: map
        real(dp), intent(in) :: f_vac
        real(dp), intent(out) :: residual, sigma
        integer, intent(out) :: status
        real(dp), allocatable :: g(:), f(:), mass(:, :, :, :)
        real(dp) :: energy, q(3), qd(3, 3), pos(3), d1(3, 3), d2(3, 3, 3), d3(3, 3, 3, 3)
        real(dp) :: dx(3, 3), bx(3), jac, r, ephi(3), weight, volume, b2, pr, info_res
        real(dp) :: dphi, y(2), wv
        integer :: e, gq, j, a, b, m, info, iterations
        logical, allocatable :: held(:)
        ! This diagnostic's mass and background field are cylindrical.
        residual = huge(1.0_dp); sigma = 0
        status = IDEAL_MAP_BAD_INPUT
        select type (chart => map%chart)
        class is (cylindrical_chart_t)
        class default
            return
        end select
        allocate (g(map%ndof), f(map%ndof), held(map%ndof))
        call ideal_map_energy(map, map%x, energy, status, g)
        if (status /= IDEAL_MAP_OK) return
        held = .false.
        do a = 1, map%nodes
            if (.not. map%boundary(a)) cycle
            do m = 1, map%nmode
                held(ideal_map_dof(map, 1, m, a)) = .true.
                held(ideal_map_dof(map, 2, m, a)) = .true.
            end do
        end do
        where (held) g = 0
        ! Reference-cell mass blocks mass(a, b, component, cell): R and R^3
        ! weights (n = 0 geometry), times the mode norms 2 pi/nfp or pi/nfp.
        allocate (mass(map%nloc, map%nloc, 3, map%cells))
        mass = 0
        do e = 1, map%cells
            do gq = 1, map%nq
                y = matmul(map%x(1:2, 1, map%conn(:, e)), map%shape(:, gq))
                wv = map%weight(gq)*map%det(gq, e)*y(1)
                do b = 1, map%nloc
                    do a = 1, map%nloc
                        mass(a, b, 1:2, e) = mass(a, b, 1:2, e) &
                            + wv*map%shape(a, gq)*map%shape(b, gq)
                        mass(a, b, 3, e) = mass(a, b, 3, e) &
                            + wv*y(1)**2*map%shape(a, gq)*map%shape(b, gq)
                    end do
                end do
            end do
        end do
        f = 0
        call real_conjugate_gradient_operator(apply_mass, g, f, 1.0e-12_dp, &
            10*map%ndof, info, iterations, info_res)
        status = IDEAL_MAP_OK
        if (info /= KRYLOV_OK) status = IDEAL_MAP_SOLVE_FAILED
        ! sigma and the volume from the current map.
        dphi = 2*acos(-1.0_dp)/(map%nfp*map%nphi)
        volume = 0; b2 = 0; pr = 0
        do e = 1, map%cells
            do j = 1, map%nphi
                do gq = 1, map%nq
                    call point_jet(map, e, gq, j, map%x, q, qd)
                    call map%chart%jet(q, pos, d1, d2, d3)
                    dx = matmul(d1, qd)
                    jac = map%sigma*(dx(1, 1)*(dx(2, 2)*dx(3, 3) - dx(2, 3)*dx(3, 2)) &
                        - dx(1, 2)*(dx(2, 1)*dx(3, 3) - dx(2, 3)*dx(3, 1)) &
                        + dx(1, 3)*(dx(2, 1)*dx(3, 2) - dx(2, 2)*dx(3, 1)))
                    bx = matmul(dx, map%b(:, gq, e))/jac
                    r = norm2(pos(1:2))
                    ephi = [-pos(2), pos(1), 0.0_dp]/r
                    weight = map%weight(gq)*dphi*jac
                    volume = volume + weight
                    b2 = b2 + weight*sum((bx - f_vac/r*ephi)**2)
                    pr = pr + weight*map%p(gq, e)
                end do
            end do
        end do
        sigma = (b2/map%mu0 + pr)/volume
        residual = sqrt(max(0.0_dp, dot_product(f, g))) &
            *sqrt(area()/acos(-1.0_dp)) &
            /(sigma*sqrt(volume))
    contains
        real(dp) function area()
            integer :: ee
            area = 0
            do ee = 1, map%cells
                area = area + sum(map%weight*map%det(:, ee))
            end do
        end function

        subroutine apply_mass(input, output)
            real(dp), intent(in) :: input(:)
            real(dp), intent(out) :: output(:)
            integer :: ee, aa, bb, cc, mm, ka, kb
            real(dp) :: modal
            output = 0
            do ee = 1, map%cells
                do mm = 1, map%nmode
                    modal = merge(2.0_dp, 1.0_dp, mm == 1)*acos(-1.0_dp)/map%nfp
                    do cc = 1, 3
                        do bb = 1, map%nloc
                            kb = ideal_map_dof(map, cc, mm, map%conn(bb, ee))
                            if (held(kb)) cycle
                            do aa = 1, map%nloc
                                ka = ideal_map_dof(map, cc, mm, map%conn(aa, ee))
                                if (held(ka)) cycle
                                output(ka) = output(ka) &
                                    + modal*mass(aa, bb, cc, ee)*input(kb)
                            end do
                        end do
                    end do
                end do
            end do
            where (held) output = input
        end subroutine
    end subroutine

    subroutine ideal_map_mode_norms(map, norms)
        ! Max-norm of the coefficients per Fourier mode, excluding the n = 0
        ! reference part of q1, q2: norms(component, mode).
        type(ideal_map_t), intent(in) :: map
        real(dp), intent(out) :: norms(3, map%nmode)
        integer :: c, m
        do m = 1, map%nmode
            do c = 1, 3
                norms(c, m) = maxval(abs(map%x(c, m, :)))
            end do
        end do
    end subroutine
end module kin6d_ideal_map
