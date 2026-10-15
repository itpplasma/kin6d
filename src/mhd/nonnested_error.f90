module kin6d_nonnested_error
    ! Equilibrated-flux (Prager--Synge) energy bound for the N = 0 restriction
    ! of the NN1 free-B lane, nn_solve(field, 0, ...): on the meridian V0 space
    !   -div(grad psi / R) = mu0 p' R   (measure dR dZ),
    ! with the Dirichlet trace of the discrete solution. FortFEM's patchwise RT_k
    ! equilibration (equilibrated_flux_2d, isoparametric cells) bounds
    !   |||psi - psi_h|||^2 = int |grad (psi - psi_h)|^2 / R dR dZ
    ! for the continuous solution psi that shares that trace. A = 1/R has the
    ! lower bound 1/max_T R from the Bernstein hull of each cell, and C_F that
    ! of the geometry bounding box, as for fixed-boundary GS. The integrals use
    ! the meridian complex's own quadrature and are not enclosed, so the claim
    ! is practical. The trace error of an approximated boundary and the n > 0
    ! modes are outside this bound.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use kin6d_equilibrium_nonnested, only: free_b_field_t
    use kin6d_triangle_basis, only: triangle_bernstein_controls
    use kin6d_gs_error_leaves, only: gs_error_coefficient, gs_error_flux, &
        gs_error_friedrichs
    use kin6d_grad_shafranov_leaves, only: gs_source
    use fortfem_meridian_feec_n, only: meridian_qgrad
    use fortfem_triangle_duffy_quadrature, only: triangle_duffy_quadrature
    use fortfem_triangle_lagrange_arbitrary_order, only: triangle_lagrange_basis_t, &
        initialize_triangle_lagrange_basis, triangle_lagrange_nodes
    use fortfem_equilibrated_flux_2d, only: equilibrated_flux_estimate, &
        equilibrated_flux_estimate_t, EQUILIBRATION_OK
    implicit none
    private
    public :: nn_error_t, nn_axisymmetric_error

    integer, parameter, public :: NN_ERROR_OK = 0, NN_ERROR_BAD_INPUT = 1
    integer, parameter, public :: NN_ERROR_BAD_CELL = 2
    integer, parameter, public :: NN_ERROR_EQUILIBRATION_FAILED = 3

    type :: nn_error_t
        ! bound >= |||psi - psi_h||| [flux units / m^(1/2)] when available.
        real(dp) :: bound = huge(1.0_dp)
        real(dp) :: flux_defect = 0, oscillation = 0, divergence_defect = 0
        character(len=40) :: claim = 'unavailable'
        character(len=40) :: quadrature_status = 'approximate'
    end type

contains

    subroutine nn_axisymmetric_error(field, solution, error, status, pprime)
        ! solution: the V0 coefficients returned by nn_solve(field, 0, ...)
        ! (real); pprime as passed to nn_solve (default field%pprime).
        type(free_b_field_t), intent(in) :: field
        complex(dp), intent(in) :: solution(:)
        type(nn_error_t), intent(out) :: error
        integer, intent(out) :: status
        real(dp), intent(in), optional :: pprime
        type(triangle_lagrange_basis_t) :: geometry_basis
        type(equilibrated_flux_estimate_t) :: estimate
        real(dp), allocatable :: xi(:), eta(:), w(:), lattice(:, :), vertices(:, :)
        real(dp), allocatable :: geometry(:, :, :), coefficient(:, :), f(:, :)
        real(dp), allocatable :: flux(:, :, :), lower(:), psi(:)
        real(dp) :: lo(2), hi(2), cf, pp, x(2), wq, gradient(2), controls(2, 10)
        real(dp) :: val(28), grad(2, 28)
        integer, allocatable :: triangles(:, :), order(:)
        integer :: nt, ngeo, nl, t, q, axis, st, failed

        status = NN_ERROR_BAD_INPUT
        pp = field%pprime
        if (present(pprime)) pp = pprime
        if (.not. ieee_is_finite(pp)) return
        if (field%fs%cyl_r0 > 0 .or. field%mu0 <= 0) return
        if (size(solution) /= field%fs%n0) return
        if (.not. all(ieee_is_finite(real(solution, dp)))) return
        if (any(abs(aimag(solution)) > 0)) return
        nt = field%fs%mm%topology%n_triangles
        ngeo = size(field%fs%mm%nodes, 2)
        nl = field%fs%nl0
        if (ngeo /= 3 .and. ngeo /= 6 .and. ngeo /= 10) return
        if (field%fs%k < 1 .or. field%fs%k > 6) return
        ! The meridian rule (degree 2k + 2g + 2) samples the data.
        call triangle_duffy_quadrature(2*field%fs%k + 2*field%fs%mm%geo_order + 2, &
            xi, eta, w, st)
        if (st /= 0 .or. size(w) /= field%fs%nq) return
        if (maxval(abs(w - field%fs%qw)) > 0) return
        call initialize_triangle_lagrange_basis(field%fs%mm%geo_order, &
            geometry_basis, st)
        if (st /= 0) return
        call triangle_lagrange_nodes(geometry_basis, lattice)
        call geometry_order(lattice, order, st)
        if (st /= 0) return
        psi = real(solution, dp)
        allocate (triangles(3, nt), vertices(2, field%fs%mm%topology%n_vertices))
        allocate (geometry(2, ngeo, nt), lower(nt), coefficient(size(w), nt))
        allocate (f(size(w), nt), flux(2, size(w), nt))
        status = NN_ERROR_BAD_CELL
        triangles = field%fs%mm%topology%triangles(1:3, :)
        lo = huge(1.0_dp); hi = -huge(1.0_dp); failed = 0
        ! Each cell maps its own copy of a shared vertex; the copies agree to
        ! rounding and are snapped to one position, as the estimator requires.
        do t = nt, 1, -1
            geometry(:, :, t) = field%fs%mm%nodes(:, order, t)
            vertices(:, triangles(:, t)) = geometry(:, 1:3, t)
        end do
        do t = 1, nt
            if (any(abs(geometry(:, 1:3, t) - vertices(:, triangles(:, t))) > &
                1.0e-12_dp*maxval(abs(geometry(:, :, t))))) failed = 1
            geometry(:, 1:3, t) = vertices(:, triangles(:, t))
            do axis = 1, 2
                call triangle_bernstein_controls(geometry(axis, :, t), &
                    controls(axis, :ngeo))
            end do
            if (minval(controls(1, :ngeo)) <= 0) failed = 1
            call gs_error_coefficient(maxval(controls(1, :ngeo)), lower(t))
            lo = min(lo, minval(controls(:, :ngeo), dim=2))
            hi = max(hi, maxval(controls(:, :ngeo), dim=2))
        end do
        if (failed /= 0 .or. any(hi <= lo)) return
        !$omp parallel do schedule(static) default(shared) &
        !$omp private(t, q, x, wq, val, grad, gradient, st) reduction(max:failed)
        do t = 1, nt
            do q = 1, size(w)
                call meridian_qgrad(field%fs, t, q, x, wq, val(:nl), grad(:, :nl), st)
                if (st /= 0) then
                    failed = 1
                    cycle
                end if
                gradient = matmul(grad(:, :nl), psi(field%fs%sp(0)%dof(:, t)))
                call gs_error_coefficient(x(1), coefficient(q, t))
                call gs_source(x(1), field%mu0, pp, 0.0_dp, f(q, t))
                call gs_error_flux(x(1), gradient(1), gradient(2), flux(1, q, t), &
                    flux(2, q, t))
            end do
        end do
        !$omp end parallel do
        if (failed /= 0) return
        call gs_error_friedrichs(hi(1), hi(1) - lo(1), hi(2) - lo(2), cf)
        call equilibrated_flux_estimate(vertices, triangles, field%fs%k, xi, eta, w, &
            coefficient, f, flux, lower, cf, estimate, st, geometry_nodes=geometry)
        status = NN_ERROR_EQUILIBRATION_FAILED
        if (st /= EQUILIBRATION_OK) return
        error%bound = estimate%bound
        error%flux_defect = estimate%flux_defect
        error%oscillation = estimate%oscillation
        error%divergence_defect = estimate%divergence_defect
        error%claim = 'practical-fixed-dirichlet-trace'
        status = NN_ERROR_OK
    end subroutine

    subroutine geometry_order(lattice, order, status)
        ! FortFEM's Lagrange nodes are lexicographic in the barycentric powers;
        ! the estimator's geometry nodes are the vertices, then per edge e (e to
        ! e + 1) its interior nodes, then the centroid (P3).
        real(dp), intent(in) :: lattice(:, :)
        integer, allocatable, intent(out) :: order(:)
        integer, intent(out) :: status
        real(dp) :: wanted(2, 10), corner(2, 3), s
        integer :: n, g, e, j, i, m

        status = 1
        n = size(lattice, 2)
        g = nint((sqrt(8.0_dp*n + 1) - 3)/2)
        if ((g + 1)*(g + 2)/2 /= n .or. g < 1 .or. g > 3) return
        corner = reshape([0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2, 3])
        wanted(:, 1:3) = corner
        m = 3
        do e = 1, 3
            do j = 1, g - 1
                s = real(j, dp)/g
                m = m + 1
                wanted(:, m) = (1 - s)*corner(:, e) + s*corner(:, mod(e, 3) + 1)
            end do
        end do
        if (g == 3) wanted(:, 10) = 1.0_dp/3
        allocate (order(n))
        order = 0
        do j = 1, n
            do i = 1, n
                if (all(abs(lattice(:, i) - wanted(:, j)) < 1.0e-12_dp)) order(j) = i
            end do
        end do
        if (any(order == 0)) return
        status = 0
    end subroutine
end module kin6d_nonnested_error
