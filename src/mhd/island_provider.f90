module kin6d_island_provider
    ! Closure A+D forward provider for inference (TIAGO DESIGN.md seam):
    ! named parameters, a warm-started stationary solve, B and p sampled in
    ! real space, equation-level JVP/VJP and provenance/status.
    !
    ! Parameters theta = [sigma, epsilon, drive, chi_par]:
    !   sigma    scale of the declared heat source, S = sigma S_0(R, Z);
    !   epsilon  amplitude of a declared boundary-data direction, the fixed
    !            dofs being base + epsilon perturbation (a resonant coil
    !            trace for an island; epsilon = 0 is the unperturbed case);
    !   drive    I_phi, or the loop voltage V when the closure drive is V;
    !   chi_par  parallel heat diffusivity (chi_perp, g_par, g_perp fixed).
    ! The residual is affine in sigma, I_phi and chi_par, so their residual
    ! derivatives are exact affine increments; epsilon and V enter through
    ! fixed dofs and their directions are columns of the exact Jacobian.
    ! dx = -J_red^{-1} dR/dtheta (+ the fixed-dof direction); the adjoint
    ! uses the transposed solve with the same factor (nc_tangent_*). The
    ! Newton iteration is never differentiated. The map is smooth wherever
    ! the reduced Jacobian is nonsingular, including across epsilon = 0 at
    ! finite chi_par/chi_perp: no topology detection enters the residual.
    !
    ! Samples are Cartesian (x, y, z) points with values [B_x, B_y, B_z, p].
    ! One provider at a time per process binds the source scale (module
    ! state read inside the residual's OpenMP regions).
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_meridian_feec_n, only: meridian_weight_i, meridian_point_basis
    use fortfem_meridian_fourier_curl, only: meridian_fourier_curl
    use fortsparse, only: csc_t, csc_matvec, csc_transpose, &
        fortsparse_status_t, status_ok
    use kin6d_mhd_mesh, only: flux_mesh_locate
    use kin6d_equilibrium_nonnested, only: free_b_field_t, nn_point
    use kin6d_equilibrium_state, only: EQ_VALUES, EQ_BR, EQ_BZ, EQ_BPHI, EQ_P
    use kin6d_nonnested_coupled, only: nc_problem_t, nc_layout_t, nc_layout, &
        nc_report_t, nc_solve, nc_continue, nc_residual, nc_store, nc_fixed, &
        nc_reduce, nc_tangent_t, nc_tangent_factor, nc_tangent_solve, &
        nc_tangent_free
    use kin6d_island_report, only: island_report_t, island_residual_report
    implicit none
    private
    public :: island_provider_t, island_provider_init, island_provider_solve
    public :: island_provider_sample, island_provider_jvp, island_provider_vjp
    public :: island_provider_free, island_provider_provenance
    public :: IP_SOURCE, IP_EPSILON, IP_DRIVE, IP_CHI_PAR, IP_NTHETA
    public :: IP_OK, IP_NOT_SOLVED, IP_OUTSIDE, IP_BAD_INPUT, IP_SOLVE_FAILED, &
        IP_DERIVATIVE_UNAVAILABLE

    integer, parameter :: IP_SOURCE = 1, IP_EPSILON = 2, IP_DRIVE = 3
    integer, parameter :: IP_CHI_PAR = 4, IP_NTHETA = 4
    integer, parameter :: IP_OK = 0, IP_NOT_SOLVED = 21, IP_OUTSIDE = 22
    integer, parameter :: IP_BAD_INPUT = 23, IP_SOLVE_FAILED = 24
    integer, parameter :: IP_DERIVATIVE_UNAVAILABLE = 25

    type :: island_provider_t
        type(free_b_field_t) :: field
        type(nc_problem_t) :: problem
        real(dp), allocatable :: x(:), base(:), perturbation(:)
        logical, allocatable :: fixed(:)
        real(dp) :: theta(IP_NTHETA) = [1, 0, 0, 1]
        real(dp) :: length_scale = 1, pressure_scale = 1, tol = 1.0e-11_dp
        type(nc_report_t) :: solve_report
        type(island_report_t) :: report
        type(nc_tangent_t) :: tangent
        type(csc_t) :: curl_t
        procedure(meridian_weight_i), pointer, nopass :: shape => null()
        integer :: status = IP_NOT_SOLVED, report_status = -1
        integer :: continued_stages = 0
        logical :: solved = .false.
        character(128) :: case_id = ''
        character(160) :: message = 'not solved'
    end type

    procedure(meridian_weight_i), pointer, save :: active_shape => null()
    real(dp), save :: active_scale = 1
contains
    function scaled_source(x) result(s)
        real(dp), intent(in) :: x(2)
        real(dp) :: s
        s = active_scale*active_shape(x)
    end function

    subroutine bind(prov)
        type(island_provider_t), intent(inout) :: prov
        active_shape => prov%shape
        active_scale = prov%theta(IP_SOURCE)
        prov%problem%source => scaled_source
        prov%problem%chi_par = prov%theta(IP_CHI_PAR)
        if (prov%problem%closure%drive == 'V') then
            prov%problem%voltage = prov%theta(IP_DRIVE)
        else
            prov%problem%current = prov%theta(IP_DRIVE)
        end if
    end subroutine

    subroutine island_provider_init(prov, field, problem, base, perturbation, &
        shape, theta, length_scale, pressure_scale, case_id, status)
        ! base: seed state whose fixed entries are the epsilon = 0 boundary
        ! data; perturbation: a state direction whose fixed entries are the
        ! declared boundary-data direction (its free entries only improve the
        ! seed). problem%source is replaced by sigma*shape.
        type(island_provider_t), intent(inout) :: prov
        type(free_b_field_t), intent(in) :: field
        type(nc_problem_t), intent(in) :: problem
        real(dp), intent(in) :: base(:), perturbation(:), theta(IP_NTHETA)
        procedure(meridian_weight_i) :: shape
        real(dp), intent(in) :: length_scale, pressure_scale
        character(*), intent(in) :: case_id
        integer, intent(out) :: status
        type(nc_layout_t) :: lay
        call island_provider_free(prov)
        status = IP_BAD_INPUT
        lay = nc_layout(field)
        if (size(base) /= lay%n .or. size(perturbation) /= lay%n) return
        if (.not. all(ieee_is_finite(base))) return
        if (.not. all(ieee_is_finite(perturbation))) return
        if (.not. all(ieee_is_finite(theta))) return
        if (length_scale <= 0 .or. pressure_scale <= 0) return
        if (associated(problem%source_3d) .or. allocated(problem%fitted)) return
        prov%field = field; prov%problem = problem
        prov%base = base; prov%perturbation = perturbation
        prov%fixed = nc_fixed(field, problem)
        prov%shape => shape
        prov%theta = theta
        prov%length_scale = length_scale; prov%pressure_scale = pressure_scale
        prov%case_id = case_id
        prov%x = base + theta(IP_EPSILON)*perturbation
        status = IP_OK
        prov%status = IP_NOT_SOLVED
    end subroutine

    subroutine island_provider_solve(prov, theta, status, tol)
        ! Warm start from the previous solution with the new boundary data;
        ! on failure, anisotropy continuation from the declared seed.
        type(island_provider_t), intent(inout) :: prov
        real(dp), intent(in) :: theta(IP_NTHETA)
        integer, intent(out) :: status
        real(dp), intent(in), optional :: tol
        real(dp), allocatable :: x(:), seed(:)
        integer :: st
        status = IP_BAD_INPUT
        if (.not. allocated(prov%base)) return
        if (.not. all(ieee_is_finite(theta))) return
        if (present(tol)) prov%tol = tol
        call nc_tangent_free(prov%tangent)
        prov%theta = theta
        call bind(prov)
        seed = prov%base + theta(IP_EPSILON)*prov%perturbation
        x = seed
        if (prov%solved) x = merge(seed, prov%x, prov%fixed)
        prov%solved = .false.; prov%continued_stages = 0
        call nc_solve(prov%field, prov%problem, x, prov%solve_report, st, &
            tol=prov%tol)
        if (st /= 0) then
            x = seed
            call nc_continue(prov%field, prov%problem, x, prov%solve_report, st, &
                tol=prov%tol)
            prov%continued_stages = prov%solve_report%stages
        end if
        if (st /= 0) then
            prov%status = IP_SOLVE_FAILED
            prov%message = prov%solve_report%message
            if (len_trim(prov%message) == 0) prov%message = 'A+D solve failed'
            status = prov%status
            return
        end if
        prov%x = x
        prov%solved = .true.
        ! Diagnostics are a separate pass: pressure-transport and force (D3)
        ! residuals are reported independently of Newton convergence.
        call island_residual_report(prov%field, prov%problem, prov%x, &
            prov%length_scale, prov%pressure_scale, prov%report, &
            prov%report_status)
        call nc_store(prov%field, prov%x, st)
        if (st /= 0) then
            prov%solved = .false.; prov%status = IP_SOLVE_FAILED
            status = prov%status
            return
        end if
        prov%status = IP_OK
        prov%message = 'A+D regularised candidate at finite chi_par/chi_perp'
        status = IP_OK
    end subroutine

    subroutine island_provider_sample(prov, xyz, values, status)
        type(island_provider_t), intent(in) :: prov
        real(dp), intent(in) :: xyz(:, :)
        real(dp), intent(out) :: values(:, :)
        integer, intent(out) :: status
        status = IP_NOT_SOLVED
        values = 0
        if (.not. prov%solved) return
        call evaluate(prov%field, xyz, values, status)
    end subroutine

    subroutine evaluate(field, xyz, values, status)
        type(free_b_field_t), intent(in) :: field
        real(dp), intent(in) :: xyz(:, :)
        real(dp), intent(out) :: values(:, :)
        integer, intent(out) :: status
        real(dp) :: val(EQ_VALUES), rzp(3)
        integer :: j, st, failure
        status = IP_BAD_INPUT
        if (size(xyz, 1) /= 3 .or. size(values, 1) /= 4) return
        if (size(values, 2) /= size(xyz, 2)) return
        failure = 0
        !$omp parallel do default(shared) private(j, st, val, rzp) &
        !$omp& schedule(static) reduction(max:failure)
        do j = 1, size(xyz, 2)
            rzp = cylindrical(xyz(:, j))
            call nn_point(field, rzp, val, st)
            if (st /= 0) then
                failure = 1
                values(:, j) = 0
            else
                values(:, j) = [val(EQ_BR)*cos(rzp(3)) - val(EQ_BPHI)*sin(rzp(3)), &
                    val(EQ_BR)*sin(rzp(3)) + val(EQ_BPHI)*cos(rzp(3)), &
                    val(EQ_BZ), val(EQ_P)]
            end if
        end do
        !$omp end parallel do
        status = IP_OK
        if (failure /= 0) status = IP_OUTSIDE
    end subroutine

    pure function cylindrical(p) result(rzp)
        real(dp), intent(in) :: p(3)
        real(dp) :: rzp(3)
        rzp = [hypot(p(1), p(2)), p(3), atan2(p(2), p(1))]
    end function

    subroutine linearise(prov, g, d, status)
        ! g(:, k): reduced residual derivative dR/dtheta_k at the solution;
        ! d(:, k): fixed-dof direction of theta_k (zero for affine data).
        type(island_provider_t), intent(inout) :: prov
        real(dp), allocatable, intent(out) :: g(:, :), d(:, :)
        integer, intent(out) :: status
        type(nc_layout_t) :: lay
        real(dp), allocatable :: r0(:), r1(:)
        real(dp) :: theta0(IP_NTHETA), delta
        integer :: k
        status = IP_NOT_SOLVED
        if (.not. prov%solved) return
        status = IP_DERIVATIVE_UNAVAILABLE
        lay = nc_layout(prov%field)
        if (.not. prov%tangent%factored) then
            call bind(prov)
            call nc_tangent_factor(prov%field, prov%problem, prov%x, &
                prov%tangent, status)
            if (status /= 0) then
                status = IP_DERIVATIVE_UNAVAILABLE
                return
            end if
        end if
        allocate (g(lay%n, IP_NTHETA), d(lay%n, IP_NTHETA))
        g = 0; d = 0
        theta0 = prov%theta
        call bind(prov)
        call nc_residual(prov%field, prov%problem, prov%x, r0, status)
        if (status /= 0) return
        do k = 1, IP_NTHETA
            select case (k)
            case (IP_EPSILON)
                d(:, k) = merge(prov%perturbation, 0.0_dp, prov%fixed)
            case (IP_DRIVE)
                if (prov%problem%closure%drive == 'V') d(lay%iv, k) = 1
            end select
            if (k == IP_EPSILON .or. (k == IP_DRIVE .and. &
                prov%problem%closure%drive == 'V')) then
                g(:, k) = csc_matvec(prov%tangent%jac, d(:, k))
                call nc_reduce(prov%field, prov%problem, g(:, k))
                g(:, k) = merge(0.0_dp, g(:, k), prov%fixed)
                cycle
            end if
            ! Affine data: an exact increment (chi_par doubled, the others
            ! shifted by max(1, |theta_k|)), not a solution difference.
            delta = max(1.0_dp, abs(theta0(k)))
            if (k == IP_CHI_PAR) delta = theta0(k)
            prov%theta = theta0; prov%theta(k) = theta0(k) + delta
            call bind(prov)
            call nc_residual(prov%field, prov%problem, prov%x, r1, status)
            prov%theta = theta0
            call bind(prov)
            if (status /= 0) return
            g(:, k) = (r1 - r0)/delta
        end do
        status = IP_OK
    end subroutine

    subroutine island_provider_jvp(prov, xyz, dtheta, dvalues, status, dx)
        ! dvalues = d[B_x, B_y, B_z, p](xyz)/dtheta . dtheta.
        type(island_provider_t), intent(inout) :: prov
        real(dp), intent(in) :: xyz(:, :), dtheta(IP_NTHETA)
        real(dp), intent(out) :: dvalues(:, :)
        integer, intent(out) :: status
        real(dp), allocatable, intent(out), optional :: dx(:)
        type(free_b_field_t) :: tangent_field
        real(dp), allocatable :: g(:, :), d(:, :), y(:)
        integer :: st
        dvalues = 0
        call linearise(prov, g, d, status)
        if (status /= IP_OK) return
        call nc_tangent_solve(prov%tangent, -matmul(g, dtheta), y, st)
        status = IP_DERIVATIVE_UNAVAILABLE
        if (st /= 0) return
        y = y + matmul(d, dtheta)
        ! B and p are linear in the state: evaluate the tangent state with
        ! no vacuum toroidal field.
        tangent_field = prov%field
        tangent_field%f_vac = 0
        call nc_store(tangent_field, y, st)
        if (st /= 0) return
        call evaluate(tangent_field, xyz, dvalues, status)
        if (present(dx)) call move_alloc(y, dx)
    end subroutine

    subroutine island_provider_vjp(prov, xyz, weight, theta_bar, status)
        ! theta_bar = (d[B_x, B_y, B_z, p](xyz)/dtheta)^T weight, by one
        ! transposed solve: c = H^T weight, J_red^T mu = E^T c,
        ! theta_bar_k = -mu . dR/dtheta_k + c . d_k.
        type(island_provider_t), intent(inout) :: prov
        real(dp), intent(in) :: xyz(:, :), weight(:, :)
        real(dp), intent(out) :: theta_bar(IP_NTHETA)
        integer, intent(out) :: status
        real(dp), allocatable :: g(:, :), d(:, :), c(:), mu(:)
        integer :: st, k
        theta_bar = 0
        status = IP_BAD_INPUT
        if (any(shape(weight) /= [4, size(xyz, 2)])) return
        call linearise(prov, g, d, status)
        if (status /= IP_OK) return
        call observation_transpose(prov, xyz, weight, c, status)
        if (status /= IP_OK) return
        call nc_tangent_solve(prov%tangent, c, mu, st, adjoint=.true.)
        status = IP_DERIVATIVE_UNAVAILABLE
        if (st /= 0) return
        do k = 1, IP_NTHETA
            theta_bar(k) = -dot_product(mu, g(:, k)) + dot_product(c, d(:, k))
        end do
        status = IP_OK
    end subroutine

    subroutine observation_transpose(prov, xyz, weight, c, status)
        ! c = H^T weight, with H the linear map from the state to Cartesian
        ! [B, p] samples (B = curl A in FortFEM's real Fourier blocks).
        type(island_provider_t), intent(inout) :: prov
        real(dp), intent(in) :: xyz(:, :), weight(:, :)
        real(dp), allocatable, intent(out) :: c(:)
        integer, intent(out) :: status
        type(nc_layout_t) :: lay
        type(csc_t) :: curl
        type(fortsparse_status_t) :: fst
        real(dp), allocatable :: bbar(:)
        integer :: j, failure, st
        lay = nc_layout(prov%field)
        status = IP_DERIVATIVE_UNAVAILABLE
        if (prov%curl_t%nrow == 0) then
            call meridian_fourier_curl(prov%field%fs, lay%modes, curl, st)
            if (st /= 0) return
            call csc_transpose(curl, prov%curl_t, fst)
            if (.not. status_ok(fst)) return
        end if
        allocate (c(lay%n), bbar(lay%nb2*lay%nb))
        c = 0; bbar = 0
        failure = 0
        !$omp parallel do default(shared) private(j, st) schedule(static) &
        !$omp& reduction(max:failure)
        do j = 1, size(xyz, 2)
            call scatter_point(prov%field, lay, xyz(:, j), weight(:, j), c, bbar, st)
            failure = max(failure, st)
        end do
        !$omp end parallel do
        status = IP_OUTSIDE
        if (failure /= 0) return
        c(lay%ia + 1:lay%ia + lay%na) = c(lay%ia + 1:lay%ia + lay%na) &
            + csc_matvec(prov%curl_t, bbar)
        status = IP_OK
    end subroutine

    subroutine scatter_point(field, lay, p, w, c, bbar, status)
        ! Transpose of nn_point's B and p evaluation at one Cartesian point.
        type(free_b_field_t), intent(in) :: field
        type(nc_layout_t), intent(in) :: lay
        real(dp), intent(in) :: p(3), w(4)
        real(dp), intent(inout) :: c(:), bbar(:)
        integer, intent(out) :: status
        real(dp) :: rzp(3), u, theta, xi, eta, x(2), wc(3), tr(lay%nb)
        real(dp) :: v2(3, field%fs%sp(2)%nloc), v0(3, field%fs%nl0), contrib
        integer :: t, n, blk, j, i
        status = 1
        rzp = cylindrical(p)
        call field%coordinates(field%centre, field%radial_scale, rzp(:2), u, theta)
        if (field%mesh%tlast(1) == 0) then
            if (u <= field%mesh%u(1)) return
        end if
        call flux_mesh_locate(field%mesh, rzp(:2), u, theta, t, xi, eta, status)
        if (status /= 0) return
        call meridian_point_basis(field%fs, field%pb, 2, t, xi, eta, x, v2, status)
        if (status /= 0) return
        call meridian_point_basis(field%fs, field%pb, 0, t, xi, eta, x, v0, status)
        if (status /= 0) return
        ! Cylindrical cotangent (R, phi, Z) of the Cartesian B weight.
        wc = [w(1)*cos(rzp(3)) + w(2)*sin(rzp(3)), &
            -w(1)*sin(rzp(3)) + w(2)*cos(rzp(3)), w(3)]
        tr(1) = 1
        do n = 1, lay%modes
            tr(2*n) = cos(n*rzp(3)); tr(2*n + 1) = sin(n*rzp(3))
        end do
        do blk = 1, lay%nb
            do j = 1, field%fs%sp(2)%nloc
                i = (blk - 1)*lay%nb2 + field%fs%sp(2)%dof(j, t)
                contrib = tr(blk)*field%fs%sp(2)%sgn(j, t)*dot_product(v2(:, j), wc)
                !$omp atomic update
                bbar(i) = bbar(i) + contrib
            end do
            do j = 1, field%fs%nl0
                i = lay%ip + (blk - 1)*lay%n0 + field%fs%sp(0)%dof(j, t)
                contrib = tr(blk)*v0(1, j)*w(4)
                !$omp atomic update
                c(i) = c(i) + contrib
            end do
        end do
        status = 0
    end subroutine

    function island_provider_provenance(prov) result(text)
        ! key=value record for the consumer's provenance ledger. The KIN6D
        ! revision is not compiled in; the consumer records its binary.
        type(island_provider_t), intent(in) :: prov
        character(:), allocatable :: text
        character(1024) :: buf
        type(nc_layout_t) :: lay
        lay = nc_layout(prov%field)
        write (buf, '(9a,i0,a,i0,a,i0,a,i0,a,4(es23.15e3,1x),a,4(es12.4,1x),a,i0,a,i0,a,es10.2,a,i0,a,i0)') &
            'provider=kin6d_island_provider case=', trim(prov%case_id), &
            ' closure=', trim(prov%problem%closure%pressure), '+', &
            trim(prov%problem%closure%current), ' drive=', &
            trim(prov%problem%closure%drive), ' order=', prov%field%fs%k, &
            ' triangles=', prov%field%fs%mm%topology%n_triangles, &
            ' modes=', lay%modes, ' unknowns=', lay%n, ' theta=', prov%theta, &
            ' chi_perp,g_par,g_perp,nphi=', prov%problem%chi_perp, &
            prov%problem%g_par, prov%problem%g_perp, real(prov%problem%nphi, dp), &
            ' status=', prov%status, ' newton_steps=', &
            prov%solve_report%newton_steps, ' tol=', prov%tol, &
            ' continuation_stages=', prov%continued_stages, &
            ' report_status=', prov%report_status
        text = trim(buf)//' claim='//trim(prov%report%claim)// &
            ' force_status='//trim(prov%report%force_status)// &
            ' certificate='//trim(prov%report%certificate_status)// &
            ' derivatives=equation-level implicit tangent/adjoint, fixed branch'
    end function

    subroutine island_provider_free(prov)
        type(island_provider_t), intent(inout) :: prov
        call nc_tangent_free(prov%tangent)
        if (allocated(prov%x)) deallocate (prov%x)
        if (allocated(prov%base)) deallocate (prov%base)
        if (allocated(prov%perturbation)) deallocate (prov%perturbation)
        if (allocated(prov%fixed)) deallocate (prov%fixed)
        prov%solved = .false.; prov%status = IP_NOT_SOLVED
        prov%curl_t = csc_t()
    end subroutine
end module
