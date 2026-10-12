module kin6d_gs_readback
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
    use kin6d_grad_shafranov, only: gs_state_t, gs_axis_candidate, &
        gs_profile_integral
    use kin6d_gs_profile_law, only: gs_profile_law_t, gs_law_primitives
    use kin6d_gs_coarea, only: gs_level_coarea
    use kin6d_gs_axis_recovery, only: gs_magnetic_axis
    use fortfem_triangle_duffy_quadrature, only: triangle_duffy_quadrature
    use kin6d_curved_triangle, only: triangle_geometry
    use kin6d_triangle_basis, only: triangle_basis_jet
    use fortnum_polynomial, only: lagrange_weights
    use fortnum_quadrature, only: gauss_legendre_ab
    implicit none
    private
    public :: gs_readback_t, read_gs_equilibrium, gs_profile_f, gs_readback_q
    public :: gs_readback_resolution
    type :: gs_readback_t
        real(dp) :: axis(2) = 0, psi_axis = 0, volume = 0, phi_edge = 0
        real(dp) :: flux_consistency = 0, contour_s = 0, integral_s = 0
        integer, allocatable :: contour_work(:)
        real(dp), allocatable :: psi(:), s_tor(:), q(:)
    end type
contains
    subroutine gs_readback_resolution(state, radial_points, angular_points)
        type(gs_state_t), intent(in) :: state
        integer, intent(inout) :: radial_points, angular_points
        ! Resolve radial changes on the element scale. Cubic profile
        ! interpolation/integration is at least the P2/P3 field order;
        ! each smooth cell arc also shrinks with the mesh.
        if (radial_points == 0) radial_points = max(17, &
            ceiling(sqrt(real(size(state%triangles, 2), dp))) + 1)
        if (angular_points == 0) angular_points = 16*(state%element_order + 2)
    end subroutine

    real(dp) function gs_readback_q(result, s) result(q)
        type(gs_readback_t), intent(in) :: result
        real(dp), intent(in) :: s
        real(dp) :: weights(4)
        integer :: i, first
        q = ieee_value(0.0_dp, ieee_quiet_nan)
        if (.not. ieee_is_finite(s)) return
        if (s < 0 .or. s > 1) return
        i = min(size(result%q) - 1, max(1, count(result%s_tor < s)))
        first = max(1, min(i - 1, size(result%q) - 3))
        call lagrange_weights(4, s, result%s_tor(first:first + 3), weights)
        q = dot_product(weights, result%q(first:first + 3))
    end function

    real(dp) function gs_profile_f(psi, f_edge, ffprime, law) result(f)
        real(dp), intent(in) :: psi, f_edge, ffprime
        type(gs_profile_law_t), intent(in), optional :: law
        real(dp) :: squared, pressure
        integer :: status
        if (present(law)) then
            call gs_law_primitives(law, psi, pressure, f, status)
            if (status /= 0) f = 0
            return
        end if
        squared = f_edge**2 + 2*gs_profile_integral( &
            [-1.0_dp, 0.0_dp], [ffprime, ffprime], psi)
        f = 0
        if (squared <= 0) return
        f = sign(sqrt(squared), f_edge)
    end function

    subroutine read_gs_equilibrium(state, f_edge, ffprime, radial_points, &
                                   angular_points, result, status, law, rho_grid)
        ! rho_grid, when present, replaces the uniform sqrt(s_pol) levels: it
        ! must increase strictly from 0 to 1 and radial_points is ignored.
        type(gs_state_t), intent(in) :: state
        real(dp), intent(in) :: f_edge, ffprime
        integer, intent(in) :: radial_points, angular_points
        type(gs_profile_law_t), intent(in), optional :: law
        real(dp), intent(in), optional :: rho_grid(:)
        type(gs_state_t) :: view
        type(gs_readback_t), intent(out) :: result
        integer, intent(out) :: status
        real(dp), allocatable :: coarea(:), x(:), y(:), w(:), density(:), rho(:)
        real(dp), allocatable :: shapes(:, :)
        real(dp) :: hessian(2, 2), gradients(2, 10), hessians(2, 2, 10)
        real(dp) :: point(2), det, flux, f, two_pi, integral, gx(2), gw(2), weights(4)
        real(dp) :: jac(2, 2), curvature(2, 2, 2), inverse(2, 2), r, orientation
        integer(int64) :: tick, tock, clock_rate
        real(dp) :: volume_sum, flux_sum
        integer :: axis_element, i, e, node, ndof, nr, na, first, k
        integer :: failed_element, failed_status, local_status, snapshot
        status = 1
        view = state
        orientation = 1
        if (maxval(view%psi) > -minval(view%psi)) orientation = -1
        view%psi = orientation*view%psi
        if (.not. view%ready .or. view%element_order < 2) return
        nr = radial_points; na = angular_points
        call gs_readback_resolution(view, nr, na)
        if (present(rho_grid)) then
            nr = size(rho_grid)
            if (nr < 5 .or. nr > 4097) return
            if (.not. all(ieee_is_finite(rho_grid))) return
            if (rho_grid(1) > 0 .or. rho_grid(1) < 0) return
            if (rho_grid(nr) > 1 .or. rho_grid(nr) < 1) return
            if (any(rho_grid(2:) <= rho_grid(:nr - 1))) return
        end if
        if (nr < 5 .or. nr > 4097) return
        if (.not. all(ieee_is_finite([f_edge, ffprime]))) return
        if (abs(f_edge) <= tiny(1.0_dp)) return
        if (any(abs(pack(view%psi, view%boundary)) > 1.0e-12_dp)) return
        call gs_axis_candidate(view, result%axis, result%psi_axis, hessian, &
            axis_element, status)
        if (status /= 0) return
        status = 1
        if (result%psi_axis >= 0) return
        if (abs(gs_profile_f(orientation*result%psi_axis, f_edge, &
            ffprime, law)) <= tiny(1.0_dp)) return
        ndof = size(view%triangles, 1)
        two_pi = 2*acos(-1.0_dp)
        allocate (result%psi(nr), result%s_tor(nr), result%q(nr), &
            coarea(nr), density(nr), rho(nr), result%contour_work(nr))
        do i = 1, nr
            rho(i) = real(i - 1, dp)/(nr - 1)
            if (present(rho_grid)) rho(i) = rho_grid(i)
            result%psi(i) = result%psi_axis*(1 - rho(i)**2)
        end do
        call system_clock(tick, clock_rate)
        call gs_level_coarea(view, result%psi, na, coarea, &
            result%axis, status, work=result%contour_work)
        call system_clock(tock)
        result%contour_s = real(tock - tick, dp)/clock_rate
        if (status /= 0) return
        do i = 1, nr
            f = gs_profile_f(orientation*result%psi(i), f_edge, ffprime, law)
            density(i) = f*coarea(i)
        end do
        result%q = orientation*density/two_pi
        result%s_tor(1) = 0
        ! Integrate the cubic interpolant of dPhi/drho. Working in sqrt(s_pol)
        ! resolves the axis cell without giving its Hessian a full O(h) weight.
        density = -2*result%psi_axis*rho*density
        call gauss_legendre_ab(2, 0.0_dp, 1.0_dp, gx, gw)
        do i = 2, nr
            result%s_tor(i) = result%s_tor(i - 1)
            first = max(1, min(i - 2, nr - 3))
            do k = 1, 2
                r = rho(i - 1) + (rho(i) - rho(i - 1))*gx(k)
                call lagrange_weights(4, r, rho(first:first + 3), weights)
                result%s_tor(i) = result%s_tor(i) + (rho(i) - rho(i - 1)) &
                    *gw(k)*dot_product(weights, density(first:first + 3))
            end do
        end do
        integral = result%s_tor(nr)
        if (abs(integral) <= tiny(1.0_dp)) then
            status = 1
            return
        end if
        result%s_tor = result%s_tor/integral
        if (.not. all(ieee_is_finite(result%s_tor))) then
            status = 1
            return
        end if
        if (any(result%s_tor(2:) <= result%s_tor(:nr - 1))) then
            status = 1
            return
        end if
        call system_clock(tick)
        call triangle_duffy_quadrature(view%element_order + 2, x, y, w, status)
        if (status /= 0) return
        allocate (shapes(ndof, size(w)))
        do node = 1, size(w)
            call triangle_basis_jet([x(node), y(node)], shapes(:, node), &
                gradients(:, :ndof), hessians(:, :, :ndof))
        end do
        ! Element-parallel volume and toroidal-flux sums: static partition per
        ! thread, reduction over the thread partial sums.
        volume_sum = 0; flux_sum = 0
        failed_element = huge(1); failed_status = 1
        !$omp parallel do default(shared) schedule(static) &
        !$omp& private(node, point, jac, curvature, det, inverse, flux, f) &
        !$omp& private(local_status, snapshot) &
        !$omp& reduction(+: volume_sum, flux_sum)
        do e = 1, size(view%triangles, 2)
            !$omp atomic read
            snapshot = failed_element
            if (e > snapshot) cycle
            do node = 1, size(w)
                if (view%curved_cell(e)) then
                    call triangle_geometry(view%points(:, view%triangles(:, e)), &
                        [x(node), y(node)], point, jac, curvature, det, inverse, &
                        local_status)
                    if (local_status /= 0) then
                        call record_failure(e, local_status)
                        exit
                    end if
                else
                    point = matmul(view%points(:, view%triangles(1:3, e)), &
                        [1 - x(node) - y(node), x(node), y(node)])
                    det = 2*view%area(e)
                end if
                flux = dot_product(shapes(:, node), view%psi(view%triangles(:, e)))
                f = gs_profile_f(orientation*flux, f_edge, ffprime, law)
                if (abs(f) <= tiny(1.0_dp)) then
                    call record_failure(e, 1)
                    exit
                end if
                volume_sum = volume_sum + two_pi*w(node)*det*point(1)
                flux_sum = flux_sum + w(node)*det*f/point(1)
            end do
        end do
        !$omp end parallel do
        if (failed_element /= huge(1)) then
            status = failed_status
            return
        end if
        result%volume = volume_sum
        result%phi_edge = flux_sum
        call system_clock(tock)
        result%integral_s = real(tock - tick, dp)/clock_rate
        result%flux_consistency = abs(integral/result%phi_edge - 1)
        call gs_magnetic_axis(view, result%axis, status)
        result%psi = orientation*result%psi
        result%psi_axis = orientation*result%psi_axis
    contains
        subroutine record_failure(element, code)
            ! The lowest failing element decides the status, as in a serial pass.
            integer, intent(in) :: element, code
            !$omp critical(gs_readback_failure)
            if (element < failed_element) then
                failed_element = element
                failed_status = code
            end if
            !$omp end critical(gs_readback_failure)
        end subroutine
    end subroutine
end module
