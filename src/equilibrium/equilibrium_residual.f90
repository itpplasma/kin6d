module kin6d_equilibrium_residual
    ! One residual for every lane (docs/equilibrium-core.md, D3): the weak
    ! Maxwell-stress form
    !   R[xi] = int (B.grad xi.B/mu0 - (|B|^2/(2 mu0) + p) div xi) dV,
    ! equal to -<J x B - grad p, xi> for smooth fields. It needs no
    ! derivatives of B or p and is evaluated on any state through its point
    ! evaluator. The curl-free background B0 = F_vac/R e_phi is removed
    ! analytically (T(B) - T(B0)); the scale is
    !   sigma = <|b|^2>/mu0 + <|p|>,  b = B - B0  (volume means),
    ! never <B^2>. Displacements: axisymmetric xi = (xi_R, 0, xi_Z) in vector
    ! P_k on a KIN6D triangulation, zero on its boundary nodes, with the norm
    ! ||xi||^2 = int |grad xi_R|^2 + |grad xi_Z|^2 dV. The reported residual is
    ! the dual norm sup R[xi]/||xi|| (Riesz representer, one sparse SPD
    ! factorisation) divided by sigma sqrt(V), V the volume of the mesh.
    ! The edge consistency |F_edge - F_vac|/F_vac is a separate gate of the
    ! caller.
    !
    ! The dual norm over the finite test space is a lower bound of the dual
    ! norm over H1_0. With error_estimate the Riesz problem of each displacement
    ! component, int A grad r . grad v = int f v + g . grad v with A = 2 pi R,
    ! f = 2 pi c1 and g = 2 pi (c2, c3) (c the pointwise stress load, see
    ! equilibrium_residual_load), is bounded by FortFEM's equilibrated-flux
    ! estimate with the discrete flux A grad r_h - g: |||r - r_h||| <= M_j.
    ! Galerkin orthogonality gives ||R||_*^2 = ||r_h||^2 + sum_j |||r_j - r_h,j|||^2,
    ! so dual <= dual_upper = sqrt(dual^2 + M^2). The estimator integrals use the
    ! quadrature rule of degree 2k + 6 on the sampled state and are not
    ! enclosed: the claim is a practical upper bound, never a certificate.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use kin6d_grad_shafranov, only: gs_state_t
    use kin6d_curved_triangle, only: triangle_geometry
    use kin6d_triangle_basis, only: triangle_basis_jet, triangle_bernstein_controls
    use kin6d_gs_csr, only: csr_matrix_t, gs_csr_build, gs_csr_spd_matrix
    use kin6d_gs_direct, only: gs_direct_t, gs_direct_available, gs_direct_factor, &
        gs_direct_solve, gs_direct_free
    use kin6d_equilibrium_state, only: equilibrium_state_t, EQ_BR, EQ_BZ, &
        EQ_BPHI, EQ_P, EQ_VALUES
    use kin6d_equilibrium_core_leaves, only: stress_residual_density
    use kin6d_equilibrium_core_leaves, only: riesz_source, riesz_coefficient, &
        riesz_flux, riesz_friedrichs
    use fortfem_triangle_duffy_quadrature, only: triangle_duffy_quadrature
    use fortfem_equilibrated_flux_2d, only: equilibrated_flux_estimate, &
        equilibrated_flux_estimate_t, EQUILIBRATION_OK
    use fortsparse, only: sparse_solver_t, sparse_factor, sparse_solve, &
        sparse_destroy, fortsparse_status_t, status_ok, FORTSPARSE_BACKEND_NATIVE, &
        FORTSPARSE_SYM_GENERAL
    implicit none
    private
    public :: equilibrium_residual_t, equilibrium_residual, equilibrium_residual_load

    type :: equilibrium_residual_t
        ! dual = ||R||_* [Pa m^(3/2)], normalized = dual/(sigma sqrt(volume)),
        ! sigma [Pa] and the means over the mesh volume [m^3].
        real(dp) :: dual = 0, sigma = 0, volume = 0, normalized = 0
        real(dp) :: b2_mean = 0, p_mean = 0
        ! error_estimate only: dual_upper >= ||R||_* over H1_0 (practical, the
        ! estimator quadrature is not enclosed), riesz_error = M above,
        ! normalized_upper = dual_upper/(sigma sqrt(volume)); bound_status 0
        ! when available, -1 when not requested, > 0 when the estimator failed.
        real(dp) :: dual_upper = huge(1.0_dp), riesz_error = huge(1.0_dp)
        real(dp) :: normalized_upper = huge(1.0_dp)
        integer :: bound_status = -1
        character(len=32) :: claim = 'unavailable-not-requested'
        character(len=32) :: quadrature_status = 'approximate'
    end type
contains
    pure subroutine stress_load(eq, x, values, c)
        ! Pointwise stress load c(5) per dR dZ of R[xi] = int 2 pi (c1 xi_R +
        ! c2 d_R xi_R + c3 d_Z xi_R + c4 d_R xi_Z + c5 d_Z xi_Z) dR dZ.
        type(equilibrium_state_t), intent(in) :: eq
        real(dp), intent(in) :: x(2), values(:)
        real(dp), intent(out) :: c(5)
        call stress_residual_density(x(1), values(EQ_BR), values(EQ_BZ), &
            values(EQ_BPHI), eq%f_vac/x(1), values(EQ_P), eq%mu0, c(1), c(2), &
            c(3), c(4), c(5))
    end subroutine

    subroutine equilibrium_residual_load(eq, points, load, status)
        ! load(:, k): the five coefficients c of the stress load at points(:, k),
        ! the integrand of R[xi] per dR dZ before the factor 2 pi (background
        ! B0 removed). The same leaf feeds the Riesz right-hand side.
        type(equilibrium_state_t), intent(in) :: eq
        real(dp), intent(in) :: points(:, :)
        real(dp), intent(out) :: load(:, :)
        integer, intent(out) :: status
        real(dp), allocatable :: values(:, :)
        integer, allocatable :: statuses(:)
        integer :: k
        status = 1
        if (.not. allocated(eq%field)) return
        if (size(points, 1) /= 2 .or. size(load, 1) /= 5) return
        if (size(load, 2) /= size(points, 2)) return
        if (any(points(1, :) <= 0)) return
        allocate (values(EQ_VALUES, size(points, 2)), statuses(size(points, 2)))
        call eq%field%sample(points, values, statuses)
        if (any(statuses /= 0)) return
        if (.not. all(ieee_is_finite(values))) return
        !$omp parallel do default(shared) private(k) schedule(static)
        do k = 1, size(points, 2)
            call stress_load(eq, points(:, k), values(:, k), load(:, k))
        end do
        !$omp end parallel do
        status = 0
    end subroutine

    subroutine equilibrium_residual(eq, mesh, report, status, residual, &
                                    error_estimate)
        ! residual(2, node) returns R[phi_node e_R], R[phi_node e_Z] (zero on
        ! boundary nodes) when present. error_estimate (default .false.) adds
        ! the equilibrated upper bound of the dual norm over H1_0.
        type(equilibrium_state_t), intent(in) :: eq
        type(gs_state_t), intent(in) :: mesh
        type(equilibrium_residual_t), intent(out) :: report
        integer, intent(out) :: status
        real(dp), intent(out), optional :: residual(:, :)
        logical, intent(in), optional :: error_estimate
        real(dp), allocatable :: xi(:), eta(:), weight(:), shapes(:, :)
        real(dp), allocatable :: dshape(:, :, :), points(:, :), values(:, :)
        real(dp), allocatable :: measure(:, :), gradient(:, :, :, :)
        real(dp), allocatable :: blocks(:, :, :), rhs(:, :), riesz(:, :)
        real(dp), allocatable :: cell_rhs(:, :, :)
        integer, allocatable :: statuses(:)
        real(dp) :: hessians(2, 2, 10), jac(2, 2), curvature(2, 2, 2), det
        real(dp) :: inverse(2, 2), c(5), b0, two_pi, wq, b2, p_sum, b2_sum, volume
        integer :: nloc, nq, ncell, nnode, e, g, a, b, k, local
        status = 1
        if (.not. allocated(eq%field) .or. .not. mesh%ready) return
        nloc = size(mesh%triangles, 1)
        if (nloc /= 6 .and. nloc /= 10) return
        ncell = size(mesh%triangles, 2)
        nnode = size(mesh%points, 2)
        two_pi = 2*acos(-1.0_dp)
        call triangle_duffy_quadrature(2*mesh%element_order + 2, xi, eta, weight, &
            local)
        if (local /= 0) return
        nq = size(weight)
        allocate (shapes(nloc, nq), dshape(2, nloc, nq))
        do g = 1, nq
            call triangle_basis_jet([xi(g), eta(g)], shapes(:, g), dshape(:, :, g), &
                hessians(:, :, :nloc))
        end do
        ! Physical quadrature points, measure 2 pi R det w and basis gradients.
        allocate (points(2, nq*ncell), measure(nq, ncell), &
            gradient(2, nloc, nq, ncell))
        local = 0
        !$omp parallel do default(shared) schedule(static) &
        !$omp& private(e, g, a, jac, curvature, det, inverse, k) reduction(max: local)
        do e = 1, ncell
            do g = 1, nq
                k = g + nq*(e - 1)
                call triangle_geometry(mesh%points(:, mesh%triangles(:, e)), &
                    [xi(g), eta(g)], points(:, k), jac, curvature, det, inverse, a)
                if (a /= 0) local = 1
                measure(g, e) = two_pi*points(1, k)*abs(det)*weight(g)
                do a = 1, nloc
                    gradient(:, a, g, e) = matmul(transpose(inverse), dshape(:, a, g))
                end do
            end do
        end do
        !$omp end parallel do
        if (local /= 0) return
        allocate (values(EQ_VALUES, size(points, 2)), statuses(size(points, 2)))
        call eq%field%sample(points, values, statuses)
        if (any(statuses /= 0)) return
        if (.not. all(ieee_is_finite(values))) return
        ! Cell contributions in parallel, scattered in cell order afterwards.
        allocate (rhs(nnode, 2), blocks(nloc, nloc, ncell), cell_rhs(nloc, 2, ncell))
        b2_sum = 0; p_sum = 0; volume = 0
        !$omp parallel do default(shared) schedule(static) &
        !$omp& private(e, g, k, a, b, b0, wq, b2, c) &
        !$omp& reduction(+: b2_sum, p_sum, volume)
        do e = 1, ncell
            blocks(:, :, e) = 0
            cell_rhs(:, :, e) = 0
            do g = 1, nq
                k = g + nq*(e - 1)
                b0 = eq%f_vac/points(1, k)
                wq = measure(g, e)
                b2 = values(EQ_BR, k)**2 + values(EQ_BZ, k)**2 &
                    + (values(EQ_BPHI, k) - b0)**2
                b2_sum = b2_sum + wq*b2
                p_sum = p_sum + wq*abs(values(EQ_P, k))
                volume = volume + wq
                call stress_load(eq, points(:, k), values(:, k), c)
                ! c is per dR dZ; the measure carries 2 pi R.
                c = c*wq/points(1, k)
                do a = 1, nloc
                    cell_rhs(a, 1, e) = cell_rhs(a, 1, e) + c(1)*shapes(a, g) &
                        + c(2)*gradient(1, a, g, e) + c(3)*gradient(2, a, g, e)
                    cell_rhs(a, 2, e) = cell_rhs(a, 2, e) &
                        + c(4)*gradient(1, a, g, e) + c(5)*gradient(2, a, g, e)
                    do b = 1, nloc
                        blocks(a, b, e) = blocks(a, b, e) + wq &
                            *dot_product(gradient(:, a, g, e), gradient(:, b, g, e))
                    end do
                end do
            end do
        end do
        !$omp end parallel do
        report%volume = volume
        rhs = 0
        do e = 1, ncell
            rhs(mesh%triangles(:, e), :) = rhs(mesh%triangles(:, e), :) &
                + cell_rhs(:, :, e)
        end do
        where (spread(mesh%boundary, 2, 2)) rhs = 0
        report%b2_mean = b2_sum/report%volume
        report%p_mean = p_sum/report%volume
        report%sigma = report%b2_mean/eq%mu0 + report%p_mean
        allocate (riesz(nnode, 2))
        call solve_riesz(mesh, blocks, rhs, riesz, local)
        if (local /= 0) return
        report%dual = sqrt(max(0.0_dp, sum(rhs*riesz)))
        if (report%sigma <= 0) return
        report%normalized = report%dual/(report%sigma*sqrt(report%volume))
        if (present(residual)) then
            if (size(residual, 1) /= 2 .or. size(residual, 2) /= nnode) return
            residual = transpose(rhs)
        end if
        report%claim = 'unavailable-not-requested'
        if (present(error_estimate)) then
            if (error_estimate) call riesz_error_bound(eq, mesh, riesz, report)
        end if
        status = 0
    end subroutine

    subroutine riesz_error_bound(eq, mesh, riesz, report)
        ! Equilibrated-flux bound of both Riesz components (see the header).
        ! A failure leaves the lower bound intact and reports bound_status > 0.
        type(equilibrium_state_t), intent(in) :: eq
        type(gs_state_t), intent(in) :: mesh
        real(dp), intent(in) :: riesz(:, :)
        type(equilibrium_residual_t), intent(inout) :: report
        type(equilibrated_flux_estimate_t) :: estimate
        real(dp), allocatable :: xi(:), eta(:), w(:), shapes(:, :), dshape(:, :, :)
        real(dp), allocatable :: x(:, :), load(:, :), coefficient(:, :), f(:, :)
        real(dp), allocatable :: flux(:, :, :), lower(:), geometry(:, :, :)
        real(dp), allocatable :: gradient(:, :, :, :)
        real(dp) :: hessians(2, 2, 10), jac(2, 2), curvature(2, 2, 2), det
        real(dp) :: inverse(2, 2), controls(2, 10), lo(2), hi(2), cf, slope(2)
        real(dp) :: squared
        integer :: nloc, nq, ncell, k, e, g, a, axis, st, failed, comp, q

        report%bound_status = 1
        report%claim = 'unavailable'
        nloc = size(mesh%triangles, 1)
        ncell = size(mesh%triangles, 2)
        k = mesh%element_order
        call triangle_duffy_quadrature(2*k + 6, xi, eta, w, st)
        if (st /= 0) return
        nq = size(w)
        allocate (shapes(nloc, nq), dshape(2, nloc, nq))
        do g = 1, nq
            call triangle_basis_jet([xi(g), eta(g)], shapes(:, g), dshape(:, :, g), &
                hessians(:, :, :nloc))
        end do
        allocate (geometry(2, nloc, ncell), lower(ncell))
        lo = huge(1.0_dp); hi = -huge(1.0_dp)
        failed = 0
        !$omp parallel do default(shared) schedule(static) &
        !$omp& private(e, axis, controls) reduction(max:failed) &
        !$omp& reduction(min:lo) reduction(max:hi)
        do e = 1, ncell
            geometry(:, :, e) = mesh%points(:, mesh%triangles(:, e))
            do axis = 1, 2
                call triangle_bernstein_controls(geometry(axis, :, e), &
                    controls(axis, :nloc))
            end do
            if (minval(controls(1, :nloc)) <= 0) failed = 1
            call riesz_coefficient(minval(controls(1, :nloc)), lower(e))
            lo = min(lo, minval(controls(:, :nloc), dim=2))
            hi = max(hi, maxval(controls(:, :nloc), dim=2))
        end do
        !$omp end parallel do
        if (failed /= 0) return
        if (any(hi <= lo)) return
        call riesz_friedrichs(lo(1), hi(1) - lo(1), hi(2) - lo(2), cf)
        allocate (x(2, nq*ncell), gradient(2, nloc, nq, ncell))
        failed = 0
        !$omp parallel do default(shared) schedule(static) &
        !$omp& private(e, g, a, q, jac, curvature, det, inverse, st) &
        !$omp& reduction(max: failed)
        do e = 1, ncell
            do g = 1, nq
                q = g + nq*(e - 1)
                call triangle_geometry(geometry(:, :, e), [xi(g), eta(g)], x(:, q), &
                    jac, curvature, det, inverse, st)
                if (st /= 0) failed = 1
                do a = 1, nloc
                    gradient(:, a, g, e) = matmul(transpose(inverse), dshape(:, a, g))
                end do
            end do
        end do
        !$omp end parallel do
        if (failed /= 0) return
        allocate (load(5, nq*ncell), coefficient(nq, ncell), f(nq, ncell), &
            flux(2, nq, ncell))
        call equilibrium_residual_load(eq, x, load, st)
        if (st /= 0) return
        squared = 0
        do comp = 1, 2
            !$omp parallel do default(shared) schedule(static) &
            !$omp& private(e, g, a, q, slope)
            do e = 1, ncell
                do g = 1, nq
                    q = g + nq*(e - 1)
                    slope = 0
                    do a = 2, nloc
                        slope = slope + (riesz(mesh%triangles(a, e), comp) &
                            - riesz(mesh%triangles(1, e), comp))*gradient(:, a, g, e)
                    end do
                    call riesz_coefficient(x(1, q), coefficient(g, e))
                    f(g, e) = 0
                    if (comp == 1) call riesz_source(load(1, q), f(g, e))
                    call riesz_flux(x(1, q), slope(1), slope(2), &
                        load(2*comp, q), load(2*comp + 1, q), flux(1, g, e), &
                        flux(2, g, e))
                end do
            end do
            !$omp end parallel do
            call equilibrated_flux_estimate(mesh%points, mesh%triangles(1:3, :), k, &
                xi, eta, w, coefficient, f, flux, lower, cf, estimate, st, &
                geometry_nodes=geometry)
            report%bound_status = 2
            if (st /= EQUILIBRATION_OK) return
            squared = squared + estimate%bound**2
        end do
        report%riesz_error = sqrt(squared)
        report%dual_upper = sqrt(report%dual**2 + squared)
        if (report%sigma > 0) report%normalized_upper = &
            report%dual_upper/(report%sigma*sqrt(report%volume))
        report%bound_status = 0
        report%claim = 'practical-upper-bound'
    end subroutine

    subroutine solve_riesz(mesh, blocks, rhs, riesz, status)
        ! A riesz = rhs with A the weighted H1_0 stiffness (Dirichlet rows
        ! identity), solved with native sparse Cholesky.
        type(gs_state_t), intent(in) :: mesh
        real(dp), intent(in) :: blocks(:, :, :), rhs(:, :)
        real(dp), intent(out) :: riesz(:, :)
        integer, intent(out) :: status
        type(csr_matrix_t) :: op
        type(gs_direct_t) :: direct
        type(sparse_solver_t) :: solver
        type(fortsparse_status_t) :: sstatus
        integer, allocatable :: rowptr(:), colidx(:)
        real(dp), allocatable :: val(:)
        logical :: ok
        integer :: j, n
        status = 1
        n = size(mesh%points, 2)
        call gs_csr_build(op, mesh%triangles, blocks, n)
        call gs_csr_spd_matrix(op, mesh%boundary, rowptr, colidx, val)
        if (gs_direct_available()) then
            call gs_direct_factor(direct, n, rowptr, colidx, val, ok)
            if (.not. ok) return
            do j = 1, 2
                call gs_direct_solve(direct, rhs(:, j), riesz(:, j), ok)
                if (.not. ok) exit
            end do
            call gs_direct_free(direct)
            if (.not. ok) return
        else
            solver%backend_id = FORTSPARSE_BACKEND_NATIVE
            solver%options%symmetry = FORTSPARSE_SYM_GENERAL
            call sparse_factor(solver, n, n, size(val), rowptr, colidx, val, sstatus)
            if (.not. status_ok(sstatus)) return
            do j = 1, 2
                call sparse_solve(solver, rhs(:, j), riesz(:, j), sstatus)
                if (.not. status_ok(sstatus)) exit
            end do
            call sparse_destroy(solver)
            if (.not. status_ok(sstatus)) return
        end if
        status = 0
    end subroutine
end module kin6d_equilibrium_residual
