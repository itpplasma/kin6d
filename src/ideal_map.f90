module kin6d_ideal_map
    ! Nonlinear 3D ideal equilibrium as frozen-in energy minimisation over maps
    ! (DESIGN 7.1, choice A). The unknown is a map from abstract coordinates
    ! s = (s1, s2, s3): (s1, s2) the master coordinates of each triangle of a
    ! 2D reference mesh and s3 = phi-hat, into Cartesian space,
    !     X(s) = H(q(s)),  q = (q1, q2, phi-hat + omega),
    ! where H is a chart (cylindrical here) and q1, q2, omega are Lagrange P_k in
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
    use kin6d_ideal_map_leaves, only: ideal_energy_density, ideal_energy_gradient, &
        ideal_energy_hessian
    use kin6d_equilibrium_chart, only: ideal_chart_t
    use fortfem_gauss_quadrature_2d, only: gauss_quadrature_triangle_t
    use fortfem_triangle_duffy_quadrature, only: triangle_duffy_quadrature
    use fortfem_element_to_csr, only: csr_matrix_t, csr_build_pattern, &
        csr_assemble_blocks
    use fortnum_krylov, only: real_conjugate_gradient_operator, KRYLOV_OK
    use fortsparse, only: sparse_solver_t, sparse_factor, sparse_solve, &
        sparse_destroy, fortsparse_status_t, status_ok, FORTSPARSE_BACKEND_SUPERLU, &
        FORTSPARSE_BACKEND_MUMPS, FORTSPARSE_SYM_GENERAL, FORTSPARSE_SYM_INDEFINITE, &
        sparse_factor_stats, sparse_factor_stats_t, sparse_refactor
    use kin6d_gs_direct, only: gs_direct_available
    implicit none
    private
    public :: ideal_map_t, ideal_frozen_t
    public :: ideal_map_init, ideal_map_freeze, ideal_map_energy, ideal_map_hessian
    public :: ideal_map_solve, ideal_map_point, ideal_map_dof, ideal_map_mode_norms
    public :: ideal_map_reduce, ideal_map_expand, ideal_map_rdof, ideal_map_axis_dof
    public :: ideal_map_prolong, ideal_map_residual, ideal_map_position_gauge
    public :: ideal_map_solve_report_t, ideal_map_eliminate_omega
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
        ! Axis expansion: polynomial degree (-1: k - 1) and power of T (set
        ! before ideal_map_init / ideal_map_freeze respectively).
        integer :: axis_degree = -1, axis_power = 1
        ! MUMPS threads. The Newton path follows near-null (relabelling)
        ! directions, so roundoff differences between factorisations move the
        ! discrete solution along them (measured: 3.9e-7 against 5.4e-7 on the
        ! finest A3 P3 seed mesh). One thread keeps the factorisation, and with
        ! the serial cell sums the whole solve, reproducible; 0 uses all
        ! threads (faster numeric factorisation, not reproducible).
        integer :: factor_threads = 1
        ! Position gauge instead of transversals + axis expansion + omega.
        logical :: position_gauge = .false.
        real(dp) :: initial_damping = 1.0e-6_dp
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
        real(dp), allocatable :: omega_gauge(:, :, :)
        ! Nodal label, axis-expansion weight T and the axis node.
        real(dp), allocatable :: chi_node(:), axis_t(:)
        integer :: axis_node = 0
        real(dp), allocatable :: reference(:, :)
        logical, allocatable :: boundary(:), fixed(:)
        real(dp), allocatable :: weight(:), shape(:, :), dshape(:, :, :)
        real(dp), allocatable :: phi(:), trig(:, :), dtrig(:, :)
        ! Frozen master-coordinate data b(3, q, e), p(q, e), reference det(q, e).
        real(dp), allocatable :: b(:, :, :), p(:, :), det(:, :)
        ! Coefficients x(3, mode, node): q1, q2 and omega.
        real(dp), allocatable :: x(:, :, :)
    end type

    type :: ideal_map_solve_report_t
        integer :: iterations = 0
        real(dp) :: energy = 0, gradient = huge(1.0_dp), step = huge(1.0_dp)
        real(dp) :: min_jacobian = 0, damping = 0, initial_gradient = 0
        ! Wall seconds: energy/gradient, Hessian assembly, factor + solve.
        real(dp) :: energy_s = 0, hessian_s = 0, factor_s = 0
        ! Backend split of factor_s where reported (MUMPS): symbolic analysis
        ! and numerical factorisation.
        real(dp) :: analysis_s = 0, numeric_s = 0
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
        integer :: a, m, k, j
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
        real(dp), allocatable :: labels(:, :)
        integer :: e, g, a, k
        status = IDEAL_MAP_OK
        allocate (labels(map%nodes, 1), map%axis_t(map%nodes))
        do a = 1, map%nodes
            call field%evaluate(map%reference(:, a), chi, grad, b3, pressure)
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
        if (map%position_gauge) call ideal_map_position_gauge(map)
    end subroutine

    subroutine ideal_map_position_gauge(map)
        ! Position gauge (as N4): full (q1, q2) per interior node and mode, no
        ! toroidal displacement (omega = 0, toroidal labels follow phi-hat). The
        ! poloidal slide along surfaces is then physical (it redistributes the
        ! toroidal flux) instead of a relabelling null, rigid motions are
        ! exact, and neither the axis expansion nor an omega pin is needed.
        type(ideal_map_t), intent(inout) :: map
        integer :: a, m, j
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
        real(dp) :: dg(2, 2), ginv(2, 2), y(2), chi_s(2), v(10), scale, w
        integer :: e, g, n
        n = map%nloc
        allocate (map%omega_gauge(n, n, map%cells))
        map%omega_gauge = 0
        scale = 0
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
                w = map%weight(g)*y(1)**2/map%det(g, e)
                scale = max(scale, y(1)**2*sum(map%b(1:2, g, e)**2)/map%det(g, e)**2)
                map%omega_gauge(:, :, e) = map%omega_gauge(:, :, e) &
                    + w*spread(v(:n), 2, n)*spread(v(:n), 1, n)
            end do
        end do
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

    subroutine local_derivatives(map, e, g, j, x, w, jac, gv, hv, want_hessian)
        ! Energy density w, Jacobian and its first/second derivatives with
        ! respect to v = (q(1:3), Q(1:3, 1:3)) (index 3 + k + 3*(i-1) for Q(k,i)).
        type(ideal_map_t), intent(in) :: map
        integer, intent(in) :: e, g, j
        real(dp), intent(in) :: x(:, :, :)
        real(dp), intent(out) :: w, jac, gv(12), hv(12, 12)
        logical, intent(in) :: want_hessian
        real(dp) :: q(3), qd(3, 3), pos(3), d1(3, 3), d2(3, 3, 3), d3(3, 3, 3, 3)
        real(dp) :: dx(3, 3), gd(3, 3), hd(9, 9), dvx(9, 12), hflat(45)
        integer :: a, i, k, m, n, r, c
        call point_jet(map, e, g, j, x, q, qd)
        call map%chart%jet(q, pos, d1, d2, d3)
        dx = matmul(d1, qd)
        call ideal_energy_density(dx(1, 1), dx(2, 1), dx(3, 1), dx(1, 2), &
            dx(2, 2), dx(3, 2), dx(1, 3), dx(2, 3), dx(3, 3), map%b(1, g, e), &
            map%b(2, g, e), map%b(3, g, e), map%p(g, e), map%mu0, &
            map%sigma, w, jac)
        call ideal_energy_gradient(dx(1, 1), dx(2, 1), dx(3, 1), dx(1, 2), &
            dx(2, 2), dx(3, 2), dx(1, 3), dx(2, 3), dx(3, 3), map%b(1, g, e), &
            map%b(2, g, e), map%b(3, g, e), map%p(g, e), map%mu0, &
            map%sigma, gd(1, 1), gd(2, 1), gd(3, 1), gd(1, 2), &
            gd(2, 2), gd(3, 2), gd(1, 3), gd(2, 3), gd(3, 3))
        ! dDX(a,i)/dv: q_m -> sum_k d2(a,k,m) Q(k,i); Q(k,i') -> d1(a,k) delta.
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
        gv = matmul(reshape(gd, [9]), dvx)
        hv = 0
        if (.not. want_hessian) return
        call hessian_leaf(hflat)
        n = 0
        do c = 1, 9
            do r = 1, c
                n = n + 1
                hd(r, c) = hflat(n); hd(c, r) = hflat(n)
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
    contains
        subroutine hessian_leaf(t)
            real(dp), intent(out) :: t(45)
            call ideal_energy_hessian(dx(1, 1), dx(2, 1), dx(3, 1), dx(1, 2), &
                dx(2, 2), dx(3, 2), dx(1, 3), dx(2, 3), dx(3, 3), map%b(1, g, e), &
                map%b(2, g, e), map%b(3, g, e), map%p(g, e), map%mu0, &
                map%sigma, t(1), t(2), t(3), t(4), t(5), t(6), t(7), &
                t(8), t(9), t(10), t(11), t(12), t(13), t(14), t(15), t(16), t(17), &
                t(18), t(19), t(20), t(21), t(22), t(23), t(24), t(25), t(26), &
                t(27), t(28), t(29), t(30), t(31), t(32), t(33), t(34), t(35), &
                t(36), t(37), t(38), t(39), t(40), t(41), t(42), t(43), t(44), t(45))
        end subroutine
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
        real(dp) :: w, jac, gv(12), hv(12, 12), f(0:3), dphi, weight, minjac
        integer :: e, g, j, a, m, c, l
        logical :: want
        want = present(gradient)
        allocate (local(3*map%nmode*map%nloc, map%cells), cell_energy(map%cells), &
            cell_jac(map%cells))
        dphi = 2*acos(-1.0_dp)/(map%nfp*map%nphi)
        ! Cell sums are combined serially below, so the energy (which decides
        ! step acceptance near roundoff) does not depend on the thread count
        ! or on the order of an OpenMP reduction.
        !$omp parallel do default(none) shared(map, x, local, want, dphi, &
        !$omp cell_energy, cell_jac) &
        !$omp private(e, g, j, a, m, c, l, w, jac, gv, hv, f, weight) schedule(static)
        do e = 1, map%cells
            local(:, e) = 0
            cell_energy(e) = 0; cell_jac(e) = huge(1.0_dp)
            do j = 1, map%nphi
                do g = 1, map%nq
                    call local_derivatives(map, e, g, j, x, w, jac, gv, hv, .false.)
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
        if (.not. want) return
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
    end subroutine

    subroutine ideal_map_hessian(map, x, matrix)
        ! Assembled coefficient Hessian (all Fourier modes coupled) as CSR. The
        ! sparsity pattern of a matrix already built for this map is reused.
        type(ideal_map_t), intent(in) :: map
        real(dp), intent(in) :: x(:, :, :)
        type(csr_matrix_t), intent(inout) :: matrix
        real(dp), allocatable :: blocks(:, :, :), fa(:, :), fb(:, :), full(:, :)
        real(dp), allocatable :: half(:, :)
        real(dp) :: w, jac, gv(12), hv(12, 12), dphi, weight, s
        integer :: e, g, j, a, m, c, r, l, lb, nl, c2, r2, nr
        nl = 3*map%nmode*map%nloc
        nr = 2*map%nmode*map%nloc + 2*map%npoly*map%nmode
        allocate (blocks(nr, nr, map%cells))
        dphi = 2*acos(-1.0_dp)/(map%nfp*map%nphi)
        !$omp parallel do default(none) shared(map, x, blocks, dphi, nl, nr) &
        !$omp private(e, g, j, a, m, c, r, l, lb, w, jac, gv, hv, weight, fa, fb, &
        !$omp s, c2, r2, full, half) schedule(static)
        do e = 1, map%cells
            allocate (fa(0:3, map%nloc*map%nmode), fb(nl, 12), full(nl, nl), &
                half(nl, nr))
            full = 0
            do j = 1, map%nphi
                do g = 1, map%nq
                    call local_derivatives(map, e, g, j, x, w, jac, gv, hv, .true.)
                    weight = map%weight(g)*dphi
                    l = 0
                    do a = 1, map%nloc
                        do m = 1, map%nmode
                            l = l + 1
                            call test_functions(map, g, j, a, m, fa(:, l))
                        end do
                    end do
                    ! fb(row, v) = sum_r f_r(row) hv(idx(c,r), v)
                    do l = 1, map%nloc*map%nmode
                        do c = 1, 3
                            lb = 3*(l - 1) + c
                            fb(lb, :) = 0
                            do r = 0, 3
                                fb(lb, :) = fb(lb, :) + fa(r, l)*hv(c + 3*r, :)
                            end do
                        end do
                    end do
                    do l = 1, map%nloc*map%nmode
                        do c2 = 1, 3
                            lb = 3*(l - 1) + c2
                            do r2 = 0, 3
                                s = weight*fa(r2, l)
                                full(:, lb) = full(:, lb) + s*fb(:, c2 + 3*r2)
                            end do
                        end do
                    end do
                end do
            end do
            do a = 1, map%nloc
                do c = 1, map%nloc
                    l = 3*(1 + map%nmode*(a - 1))
                    r = 3*(1 + map%nmode*(c - 1))
                    full(l, r) = full(l, r) + map%omega_gauge(a, c, e)
                end do
            end do
            ! Reduce to (nu, omega): P maps (nu, omega) of local node-mode l to
            ! (d1 nu, d2 nu, omega).
            do l = 1, map%nloc*map%nmode
                a = map%conn(1 + (l - 1)/map%nmode, e)
                half(:, 2*l - 1) = matmul(full(:, 3*l - 2:3*l), map%basis(:, 1, a))
                half(:, 2*l) = matmul(full(:, 3*l - 2:3*l), map%basis(:, 2, a))
            end do
            do m = 1, map%nmode
                do j = 1, map%npoly
                    do c = 1, 2
                        r = 2*map%nmode*map%nloc + c + 2*(j - 1) + 2*map%npoly*(m - 1)
                        half(:, r) = 0
                        do a = 1, map%nloc
                            l = m + map%nmode*(a - 1)
                            half(:, r) = half(:, r) &
                                + map%axis_basis(j, map%conn(a, e))*full(:, 3*l - 3 + c)
                        end do
                    end do
                end do
            end do
            do l = 1, map%nloc*map%nmode
                a = map%conn(1 + (l - 1)/map%nmode, e)
                blocks(2*l - 1, :, e) = matmul(map%basis(:, 1, a), half(3*l - 2:3*l, :))
                blocks(2*l, :, e) = matmul(map%basis(:, 2, a), half(3*l - 2:3*l, :))
            end do
            do m = 1, map%nmode
                do j = 1, map%npoly
                    do c = 1, 2
                        r = 2*map%nmode*map%nloc + c + 2*(j - 1) + 2*map%npoly*(m - 1)
                        blocks(r, :, e) = 0
                        do a = 1, map%nloc
                            l = m + map%nmode*(a - 1)
                            blocks(r, :, e) = blocks(r, :, e) &
                                + map%axis_basis(j, map%conn(a, e))*half(3*l - 3 + c, :)
                        end do
                    end do
                end do
            end do
            deallocate (fa, fb, full, half)
        end do
        !$omp end parallel do
        if (matrix%n /= map%nred .or. matrix%nloc /= nr) &
            call csr_build_pattern(map%nred, map%redconn, .false., matrix)
        call csr_assemble_blocks(matrix, map%redconn, blocks)
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
        real(dp) :: cap_x, cap_w, t0, trial_minjac
        integer :: it, i, k, tries, estatus
        logical :: talk, accepted, factored, ofactored, nested
        talk = .false.
        factored = .false.; ofactored = .false.
        if (present(verbose)) talk = verbose
        nested = map%omega_elimination .and. .not. map%position_gauge
        allocate (grad(map%ndof), full(map%ndof), gred(map%nred), step(map%nred), &
            diag(map%nred), rhs(map%nred), shift(map%nred))
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
        mu = map%initial_damping
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
            do tries = 1, 30
                t0 = wall_seconds()
                shift = merge(mu*(abs(diag) + 1.0e-12_dp*scale), 0.0_dp, damped)
                call factor_masked(map, matrix, .not. map%fixed, shift, solver, &
                    factored, estatus, report)
                if (estatus /= 0) then
                    mu = mu*10; cycle
                end if
                rhs = -gred
                call sparse_solve(solver, rhs, step, sstatus)
                report%factor_s = report%factor_s + wall_seconds() - t0
                if (.not. status_ok(sstatus)) then
                    status = IDEAL_MAP_SOLVE_FAILED; call cleanup(); return
                end if
                where (map%fixed) step = 0
                call ideal_map_expand(map, step, full)
                ! Trust region: near-null (relabelling) directions can carry
                ! large Newton components; cap the displacement and rotation
                ! (nested: omega is re-eliminated, only nu is capped).
                alpha = 1
                do i = 1, map%nodes
                    k = ideal_map_dof(map, 1, 1, i)
                    alpha = min(alpha, cap_x/max(tiny(1.0_dp), &
                        maxval(abs(full(k:k + 3*map%nmode - 1:3))), &
                        maxval(abs(full(k + 1:k + 3*map%nmode - 1:3)))), &
                        cap_w/max(tiny(1.0_dp), maxval(abs(full(k + 2:k + 3*map%nmode &
                        - 1:3)))))
                end do
                full = alpha*full; step = alpha*step
                trial = map%x + reshape(full, shape(map%x))
                if (nested) then
                    call eliminate_omega(map, trial, omask, omatrix, osolver, ofactored, &
                        trial_energy, rhs, trial_minjac, report, estatus)
                    if (estatus == IDEAL_MAP_SOLVE_FAILED) then
                        status = estatus; call cleanup(); return
                    end if
                    ! Both omega gradients vanish: rank by the nu part.
                    where (omask) step = 0
                else
                    t0 = wall_seconds()
                    call ideal_map_energy(map, trial, trial_energy, estatus)
                    report%energy_s = report%energy_s + wall_seconds() - t0
                end if
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
                    if (.not. nested) then
                        call ideal_map_energy(map, trial, trial_energy, estatus, grad)
                        call ideal_map_reduce(map, grad, rhs)
                        where (map%fixed) rhs = 0
                    end if
                    if (dot_product(gred + rhs, step) <= 0) then
                        accepted = .true.; exit
                    end if
                end if
                mu = mu*10
            end do
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
            gnorm = norm2(gred)
            report%energy = energy; report%gradient = gnorm/g0
            report%min_jacobian = minjac; report%damping = mu
            if (talk) print '(a,i3,5es12.4,i3,3es10.2)', 'ideal map newton', it, energy, &
                gnorm/g0, maxval(abs(step(1::2))), maxval(abs(step(2::2))), mu, &
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
        allocate (grad(map%ndof), full(map%ndof), rhs(map%nred), dw(map%nred), &
            shift(map%nred))
        shift = 0
        call evaluate()
        if (info /= IDEAL_MAP_OK) return
        t0 = wall_seconds()
        call ideal_map_hessian(map, x, matrix)
        report%hessian_s = report%hessian_s + wall_seconds() - t0
        t0 = wall_seconds()
        call factor_masked(map, matrix, omask, shift, slv, done, info, report)
        if (info /= 0) then
            info = IDEAL_MAP_SOLVE_FAILED; return
        end if
        rhs = merge(-gred, 0.0_dp, omask)
        call sparse_solve(slv, rhs, dw, sstatus)
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
        end subroutine
    end subroutine

    subroutine factor_masked(map, matrix, keep, shift, slv, done, info, report)
        ! Factorise the reduced Hessian restricted to the keep rows/columns
        ! (identity elsewhere), plus shift on the diagonal. done: slv holds a
        ! factorisation of this pattern (only the numerical one is redone).
        type(ideal_map_t), intent(in) :: map
        type(csr_matrix_t), intent(in) :: matrix
        logical, intent(in) :: keep(:)
        real(dp), intent(in) :: shift(:)
        type(sparse_solver_t), intent(inout) :: slv
        logical, intent(inout) :: done
        integer, intent(out) :: info
        type(ideal_map_solve_report_t), intent(inout) :: report
        type(sparse_factor_stats_t) :: stats
        type(fortsparse_status_t) :: sstatus
        logical :: reused
        real(dp), allocatable :: val(:)
        integer :: row, kk
        val = matrix%val
        do row = 1, map%nred
            do kk = matrix%rowptr(row), matrix%rowptr(row + 1) - 1
                if (.not. keep(row)) then
                    val(kk) = merge(1.0_dp, 0.0_dp, matrix%colidx(kk) == row)
                else if (.not. keep(matrix%colidx(kk))) then
                    val(kk) = 0
                else if (matrix%colidx(kk) == row) then
                    val(kk) = val(kk) + shift(row)
                end if
            end do
        end do
        ! MUMPS LDL^T (nested dissection) when KIN6D is built with MUMPS;
        ! serial SuperLU LU otherwise. The damped system is symmetric and
        ! may be indefinite.
        if (gs_direct_available()) then
            slv%backend_id = FORTSPARSE_BACKEND_MUMPS
            slv%options%symmetry = FORTSPARSE_SYM_INDEFINITE
            slv%options%nthreads = map%factor_threads
        else
            slv%backend_id = FORTSPARSE_BACKEND_SUPERLU
            slv%options%symmetry = FORTSPARSE_SYM_GENERAL
        end if
        ! The pattern is fixed per map: after the first factorisation only
        ! the numerical factorisation is redone (MUMPS keeps its analysis).
        reused = .false.
        if (done) then
            call sparse_refactor(slv, map%nred, map%nred, size(val), &
                matrix%rowptr, matrix%colidx, val, sstatus, reused)
        else
            call sparse_factor(slv, map%nred, map%nred, size(val), &
                matrix%rowptr, matrix%colidx, val, sstatus)
        end if
        info = merge(0, 1, status_ok(sstatus))
        if (info /= 0) then
            if (done) call sparse_destroy(slv)
            done = .false.
            return
        end if
        done = .true.
        call sparse_factor_stats(slv, stats, sstatus)
        if (status_ok(sstatus)) then
            if (.not. reused) report%analysis_s = report%analysis_s &
                + max(0.0_dp, stats%analysis_s)
            report%numeric_s = report%numeric_s + max(0.0_dp, stats%factor_s)
        end if
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
