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
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use kin6d_grad_shafranov, only: gs_state_t
    use kin6d_curved_triangle, only: triangle_geometry
    use kin6d_triangle_basis, only: triangle_basis_jet
    use kin6d_gs_csr, only: csr_matrix_t, gs_csr_build, gs_csr_spd_matrix
    use kin6d_gs_direct, only: gs_direct_t, gs_direct_available, gs_direct_factor, &
        gs_direct_solve, gs_direct_free
    use kin6d_equilibrium_state, only: equilibrium_state_t, EQ_BR, EQ_BZ, &
        EQ_BPHI, EQ_P, EQ_VALUES
    use kin6d_equilibrium_core_leaves, only: stress_residual_density
    use fortfem_triangle_duffy_quadrature, only: triangle_duffy_quadrature
    use fortsparse, only: sparse_solver_t, sparse_factor, sparse_solve, &
        sparse_destroy, fortsparse_status_t, status_ok, &
        FORTSPARSE_BACKEND_SUPERLU, FORTSPARSE_SYM_GENERAL
    implicit none
    private
    public :: equilibrium_residual_t, equilibrium_residual

    type :: equilibrium_residual_t
        ! dual = ||R||_* [Pa m^(3/2)], normalized = dual/(sigma sqrt(volume)),
        ! sigma [Pa] and the means over the mesh volume [m^3].
        real(dp) :: dual = 0, sigma = 0, volume = 0, normalized = 0
        real(dp) :: b2_mean = 0, p_mean = 0
    end type
contains
    subroutine equilibrium_residual(eq, mesh, report, status, residual)
        ! residual(2, node) returns R[phi_node e_R], R[phi_node e_Z] (zero on
        ! boundary nodes) when present.
        type(equilibrium_state_t), intent(in) :: eq
        type(gs_state_t), intent(in) :: mesh
        type(equilibrium_residual_t), intent(out) :: report
        integer, intent(out) :: status
        real(dp), intent(out), optional :: residual(:, :)
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
                call stress_residual_density(points(1, k), values(EQ_BR, k), &
                    values(EQ_BZ, k), values(EQ_BPHI, k), b0, values(EQ_P, k), &
                    eq%mu0, c(1), c(2), c(3), c(4), c(5))
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
        status = 0
    end subroutine

    subroutine solve_riesz(mesh, blocks, rhs, riesz, status)
        ! A riesz = rhs with A the weighted H1_0 stiffness (Dirichlet rows
        ! identity): MUMPS Cholesky when built, serial SuperLU otherwise.
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
            solver%backend_id = FORTSPARSE_BACKEND_SUPERLU
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
