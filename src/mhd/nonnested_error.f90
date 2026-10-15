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
    !
    ! nt_axisymmetric_error does the same for the NN2 transport solve
    ! nt_solve: -div(K grad p) = f0 - div f with the Fourier-V0 measure R dR dZ
    ! times the phi mean. For K independent of phi (checked at the nphi
    ! sample angles) the N = 0 block of p_h is the Galerkin solution of the
    ! planar problem -div(A grad p) = R f0 - div(R f_pol) with the SPD tensor
    ! A = R K_pol, decoupled from the N > 0 modes; FortFEM's tensor
    ! equilibration bounds int R grad e . K_pol grad e dR dZ for the
    ! continuous solution with the discrete Dirichlet trace. A_T,min =
    ! R_min,T chi_lower (Cauchy interlacing: K >= chi_lower I gives K_pol >=
    ! chi_lower I). A phi-dependent K couples the modes and no bound is
    ! claimed for modes > 0; with modes = 0 the bound is for the planar model
    ! with the phi-sampled mean coefficient (scope string).
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use kin6d_equilibrium_nonnested, only: free_b_field_t
    use kin6d_triangle_basis, only: triangle_bernstein_controls
    use kin6d_gs_error_leaves, only: gs_error_coefficient, gs_error_flux, &
        gs_error_friedrichs
    use kin6d_grad_shafranov_leaves, only: gs_source
    use kin6d_nonnested_leaves, only: transport_error_data, &
        transport_error_lower, transport_error_friedrichs
    use fortfem_meridian_feec_n, only: meridian_feec_t
    use fortfem_meridian_fourier_scalar, only: meridian_fourier_tensor_t, &
        meridian_fourier_source_t
    use fortfem_meridian_feec_n, only: meridian_qgrad
    use fortfem_triangle_duffy_quadrature, only: triangle_duffy_quadrature
    use fortfem_triangle_lagrange_arbitrary_order, only: triangle_lagrange_basis_t, &
        initialize_triangle_lagrange_basis, triangle_lagrange_nodes
    use fortfem_equilibrated_flux_2d, only: equilibrated_flux_estimate, &
        equilibrated_flux_estimate_tensor, equilibrated_flux_estimate_t, &
        EQUILIBRATION_OK
    implicit none
    private
    public :: nn_error_t, nn_axisymmetric_error, nt_axisymmetric_error

    integer, parameter, public :: NN_ERROR_OK = 0, NN_ERROR_BAD_INPUT = 1
    integer, parameter, public :: NN_ERROR_BAD_CELL = 2
    integer, parameter, public :: NN_ERROR_EQUILIBRATION_FAILED = 3
    integer, parameter, public :: NN_ERROR_MODE_COUPLING = 4

    type :: nn_error_t
        ! bound >= |||psi - psi_h||| [flux units / m^(1/2)] when available.
        real(dp) :: bound = huge(1.0_dp)
        real(dp) :: flux_defect = 0, oscillation = 0, divergence_defect = 0
        character(len=40) :: claim = 'unavailable'
        character(len=40) :: quadrature_status = 'approximate'
        character(len=40) :: scope = 'unavailable'
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
        type(equilibrated_flux_estimate_t) :: estimate
        real(dp), allocatable :: xi(:), eta(:), w(:), vertices(:, :)
        real(dp), allocatable :: geometry(:, :, :), coefficient(:, :), f(:, :)
        real(dp), allocatable :: flux(:, :, :), lower(:), psi(:), rlo(:), rhi(:)
        real(dp) :: lo(2), hi(2), cf, pp, x(2), wq, gradient(2)
        real(dp) :: val(28), grad(2, 28)
        integer, allocatable :: triangles(:, :)
        integer :: nt, nl, t, q, st, failed

        status = NN_ERROR_BAD_INPUT
        pp = field%pprime
        if (present(pprime)) pp = pprime
        if (.not. ieee_is_finite(pp)) return
        if (field%fs%cyl_r0 > 0 .or. field%mu0 <= 0) return
        if (size(solution) /= field%fs%n0) return
        if (.not. all(ieee_is_finite(real(solution, dp)))) return
        if (any(abs(aimag(solution)) > 0)) return
        nt = field%fs%mm%topology%n_triangles
        nl = field%fs%nl0
        call meridian_rule(field%fs, xi, eta, w, st)
        if (st /= 0) return
        psi = real(solution, dp)
        allocate (coefficient(size(w), nt), f(size(w), nt), flux(2, size(w), nt))
        status = NN_ERROR_BAD_CELL
        call meridian_cells(field%fs, vertices, triangles, geometry, rlo, rhi, lo, &
            hi, st)
        if (st /= 0) return
        allocate (lower(nt))
        do t = 1, nt
            call gs_error_coefficient(rhi(t), lower(t))
        end do
        failed = 0
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
        error%scope = 'n0-restriction'
        status = NN_ERROR_OK
    end subroutine

    subroutine nt_axisymmetric_error(fs, modes, nphi, conductivity, chi_lower, &
            source, p, error, status)
        ! p: the Fourier-V0 solution of nt_solve(fs, modes, nphi, conductivity,
        ! source, ...); chi_lower > 0 a lower bound of the eigenvalues of K
        ! (chi_perp for nt_conductivity_t). Bounds the N = 0 block.
        type(meridian_feec_t), intent(in) :: fs
        integer, intent(in) :: modes, nphi
        class(meridian_fourier_tensor_t), intent(in) :: conductivity
        real(dp), intent(in) :: chi_lower
        class(meridian_fourier_source_t), intent(in) :: source
        real(dp), intent(in) :: p(:)
        type(nn_error_t), intent(out) :: error
        integer, intent(out) :: status
        type(equilibrated_flux_estimate_t) :: estimate
        real(dp), allocatable :: xi(:), eta(:), w(:), vertices(:, :)
        real(dp), allocatable :: geometry(:, :, :), tensor(:, :, :, :), f(:, :)
        real(dp), allocatable :: flux(:, :, :), lower(:), rlo(:), rhi(:)
        real(dp), allocatable :: xq(:, :), k(:, :, :), f0(:), fl(:, :)
        real(dp) :: lo(2), hi(2), cf, x(2), wq, gradient(2), val(28), grad(2, 28)
        real(dp) :: km(3, 3), fm0, fm(3), rr, scale, a11, a12, a22
        integer, allocatable :: triangles(:, :)
        integer :: nt, nl, t, q, j, pt, st, failed, coupled, bad
        real(dp), parameter :: pi = acos(-1.0_dp)

        status = NN_ERROR_BAD_INPUT
        if (modes < 0 .or. nphi < 2*modes + 1) return
        if (.not. (chi_lower > 0) .or. .not. ieee_is_finite(chi_lower)) return
        if (size(p) /= fs%n0*(2*modes + 1)) return
        if (.not. all(ieee_is_finite(p))) return
        nt = fs%mm%topology%n_triangles
        nl = fs%nl0
        call meridian_rule(fs, xi, eta, w, st)
        if (st /= 0) return
        status = NN_ERROR_BAD_CELL
        call meridian_cells(fs, vertices, triangles, geometry, rlo, rhi, lo, hi, st)
        if (st /= 0) return
        allocate (tensor(2, 2, size(w), nt), f(size(w), nt), flux(2, size(w), nt))
        allocate (lower(nt))
        failed = 0
        coupled = 0
        ! The reduction copies start at -huge; bad is the per-cell flag.
        !$omp parallel default(shared) private(t, q, j, pt, x, wq, val, grad, &
        !$omp gradient, st, xq, k, f0, fl, km, fm0, fm, rr, scale, a11, a12, a22, &
        !$omp bad) reduction(max:failed, coupled)
        allocate (xq(3, size(w)*nphi), k(3, 3, size(w)*nphi), f0(size(w)*nphi), &
            fl(3, size(w)*nphi))
        !$omp do schedule(static)
        do t = 1, nt
            call transport_error_lower(rlo(t), chi_lower, lower(t))
            bad = 0
            do q = 1, size(w)
                call meridian_qgrad(fs, t, q, x, wq, val(:nl), grad(:, :nl), st)
                if (st /= 0) bad = 1
                do j = 1, nphi
                    xq(:, (q - 1)*nphi + j) = [x, 2*pi*(j - 1)/nphi]
                end do
            end do
            if (bad == 0) call conductivity%tensor(t, xq, k, st)
            if (st /= 0) bad = 1
            if (bad == 0) call source%source(t, xq, f0, fl, st)
            if (st /= 0) bad = 1
            failed = max(failed, bad)
            if (bad /= 0) cycle
            do q = 1, size(w)
                call meridian_qgrad(fs, t, q, x, wq, val(:nl), grad(:, :nl), st)
                rr = x(1)
                if (fs%cyl_r0 > 0) rr = fs%cyl_r0
                ! phi means of K, f0 and f at the sample angles.
                pt = (q - 1)*nphi
                km = sum(k(:, :, pt + 1:pt + nphi), dim=3)/nphi
                fm0 = sum(f0(pt + 1:pt + nphi))/nphi
                fm = sum(fl(:, pt + 1:pt + nphi), dim=2)/nphi
                scale = maxval(abs(km))
                do j = 1, nphi
                    if (any(abs(k(:, :, pt + j) - km) > 1.0e-12_dp*scale)) &
                        coupled = 1
                end do
                gradient = matmul(grad(:, :nl), p(fs%sp(0)%dof(:, t)))
                ! (R, phi, Z) physical order: the poloidal block is (1, 3).
                call transport_error_data(rr, km(1, 1), km(1, 3), km(3, 3), fm0, &
                    fm(1), fm(3), gradient(1), gradient(2), a11, a12, a22, &
                    f(q, t), flux(1, q, t), flux(2, q, t))
                tensor(:, :, q, t) = reshape([a11, a12, a12, a22], [2, 2])
            end do
        end do
        !$omp end do
        !$omp end parallel
        if (failed /= 0) return
        if (coupled /= 0 .and. modes > 0) then
            status = NN_ERROR_MODE_COUPLING
            error%claim = 'unavailable-mode-coupling'
            return
        end if
        call transport_error_friedrichs(minval(lower), hi(1) - lo(1), hi(2) - lo(2), &
            cf)
        call equilibrated_flux_estimate_tensor(vertices, triangles, fs%k, xi, eta, &
            w, tensor, f, flux, lower, cf, estimate, st, geometry_nodes=geometry)
        status = NN_ERROR_EQUILIBRATION_FAILED
        if (st /= EQUILIBRATION_OK) return
        error%bound = estimate%bound
        error%flux_defect = estimate%flux_defect
        error%oscillation = estimate%oscillation
        error%divergence_defect = estimate%divergence_defect
        error%claim = 'practical-fixed-dirichlet-trace'
        if (coupled == 0 .and. nphi > 1) then
            error%scope = 'n0-block-of-phi-independent-k'
        else
            error%scope = 'phi-sampled-mean-planar-model'
        end if
        status = NN_ERROR_OK
    end subroutine

    subroutine meridian_rule(fs, xi, eta, w, status)
        ! The meridian rule (degree 2k + 2g + 2), which samples the data.
        type(meridian_feec_t), intent(in) :: fs
        real(dp), allocatable, intent(out) :: xi(:), eta(:), w(:)
        integer, intent(out) :: status
        integer :: ngeo, st

        status = 1
        ngeo = size(fs%mm%nodes, 2)
        if (ngeo /= 3 .and. ngeo /= 6 .and. ngeo /= 10) return
        if (fs%k < 1 .or. fs%k > 6) return
        call triangle_duffy_quadrature(2*fs%k + 2*fs%mm%geo_order + 2, xi, eta, w, &
            st)
        if (st /= 0 .or. size(w) /= fs%nq) return
        if (maxval(abs(w - fs%qw)) > 0) return
        status = 0
    end subroutine

    subroutine meridian_cells(fs, vertices, triangles, geometry, rlo, rhi, lo, hi, &
            status)
        ! Estimator geometry of the meridian mesh: snapped vertices, cells,
        ! geometry nodes in the estimator order, per-cell R range of the
        ! Bernstein hull (the metric radius R0 on the cylinder chart) and the
        ! bounding box lo, hi.
        type(meridian_feec_t), intent(in) :: fs
        real(dp), allocatable, intent(out) :: vertices(:, :), geometry(:, :, :)
        integer, allocatable, intent(out) :: triangles(:, :)
        real(dp), allocatable, intent(out) :: rlo(:), rhi(:)
        real(dp), intent(out) :: lo(2), hi(2)
        integer, intent(out) :: status
        type(triangle_lagrange_basis_t) :: geometry_basis
        real(dp), allocatable :: lattice(:, :)
        real(dp) :: controls(2, 10)
        integer, allocatable :: order(:)
        integer :: nt, ngeo, t, axis, st, failed

        status = 1
        nt = fs%mm%topology%n_triangles
        ngeo = size(fs%mm%nodes, 2)
        call initialize_triangle_lagrange_basis(fs%mm%geo_order, geometry_basis, st)
        if (st /= 0) return
        call triangle_lagrange_nodes(geometry_basis, lattice)
        call geometry_order(lattice, order, st)
        if (st /= 0) return
        allocate (triangles(3, nt), vertices(2, fs%mm%topology%n_vertices))
        allocate (geometry(2, ngeo, nt), rlo(nt), rhi(nt))
        triangles = fs%mm%topology%triangles(1:3, :)
        lo = huge(1.0_dp); hi = -huge(1.0_dp); failed = 0
        ! Each cell maps its own copy of a shared vertex; the copies agree to
        ! rounding and are snapped to one position, as the estimator requires.
        do t = nt, 1, -1
            geometry(:, :, t) = fs%mm%nodes(:, order, t)
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
            rlo(t) = minval(controls(1, :ngeo))
            rhi(t) = maxval(controls(1, :ngeo))
            if (fs%cyl_r0 > 0) then
                rlo(t) = fs%cyl_r0
                rhi(t) = fs%cyl_r0
            end if
            if (rlo(t) <= 0) failed = 1
            lo = min(lo, minval(controls(:, :ngeo), dim=2))
            hi = max(hi, maxval(controls(:, :ngeo), dim=2))
        end do
        if (failed /= 0 .or. any(hi <= lo)) return
        status = 0
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
