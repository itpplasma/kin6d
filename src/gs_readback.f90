module kin6d_gs_readback
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
    use kin6d_grad_shafranov, only: gs_state_t, gs_axis_candidate, &
        gs_profile_integral
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

    real(dp) function gs_profile_f(psi, f_edge, ffprime) result(f)
        real(dp), intent(in) :: psi, f_edge, ffprime
        real(dp) :: squared
        squared = f_edge**2 + 2*gs_profile_integral( &
            [-1.0_dp, 0.0_dp], [ffprime, ffprime], psi)
        f = 0
        if (squared <= 0) return
        f = sign(sqrt(squared), f_edge)
    end function

    subroutine read_gs_equilibrium(state, f_edge, ffprime, radial_points, &
                                   angular_points, result, status)
        type(gs_state_t), intent(in) :: state
        real(dp), intent(in) :: f_edge, ffprime
        integer, intent(in) :: radial_points, angular_points
        type(gs_readback_t), intent(out) :: result
        integer, intent(out) :: status
        real(dp), allocatable :: coarea(:), x(:), y(:), w(:), density(:), rho(:)
        real(dp), allocatable :: shapes(:, :)
        real(dp) :: hessian(2, 2), gradients(2, 10), hessians(2, 2, 10)
        real(dp) :: point(2), det, flux, f, two_pi, integral, gx(2), gw(2), weights(4)
        real(dp) :: jac(2, 2), curvature(2, 2, 2), inverse(2, 2), r
        integer(int64) :: tick, tock, clock_rate
        integer :: axis_element, i, e, node, ndof, nr, na, first, k
        status = 1
        if (.not. state%ready .or. state%element_order < 2) return
        nr = radial_points; na = angular_points
        call gs_readback_resolution(state, nr, na)
        if (nr < 5 .or. nr > 4097) return
        if (.not. all(ieee_is_finite([f_edge, ffprime]))) return
        if (abs(f_edge) <= tiny(1.0_dp)) return
        if (any(abs(pack(state%psi, state%boundary)) > 1.0e-12_dp)) return
        call gs_axis_candidate(state, result%axis, result%psi_axis, hessian, &
            axis_element, status)
        if (status /= 0) return
        status = 1
        if (result%psi_axis >= 0) return
        if (f_edge**2 + 2*ffprime*result%psi_axis <= 0) return
        ndof = size(state%triangles, 1)
        two_pi = 2*acos(-1.0_dp)
        allocate (result%psi(nr), result%s_tor(nr), result%q(nr), &
            coarea(nr), density(nr), rho(nr), result%contour_work(nr))
        do i = 1, nr
            rho(i) = real(i - 1, dp)/(nr - 1)
            result%psi(i) = result%psi_axis*(1 - rho(i)**2)
        end do
        call system_clock(tick, clock_rate)
        call gs_level_coarea(state, result%psi, na, coarea, &
            result%axis, status, work=result%contour_work)
        call system_clock(tock)
        result%contour_s = real(tock - tick, dp)/clock_rate
        if (status /= 0) return
        do i = 1, nr
            f = gs_profile_f(result%psi(i), f_edge, ffprime)
            density(i) = f*coarea(i)
        end do
        result%q = density/two_pi
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
        call triangle_duffy_quadrature(state%element_order + 2, x, y, w, status)
        if (status /= 0) return
        allocate (shapes(ndof, size(w)))
        do node = 1, size(w)
            call triangle_basis_jet([x(node), y(node)], shapes(:, node), &
                gradients(:, :ndof), hessians(:, :, :ndof))
        end do
        do e = 1, size(state%triangles, 2)
            do node = 1, size(w)
                if (state%curved_cell(e)) then
                    call triangle_geometry(state%points(:, state%triangles(:, e)), &
                        [x(node), y(node)], point, jac, curvature, det, inverse, status)
                    if (status /= 0) return
                else
                    point = matmul(state%points(:, state%triangles(1:3, e)), &
                        [1 - x(node) - y(node), x(node), y(node)])
                    det = 2*state%area(e)
                end if
                flux = dot_product(shapes(:, node), state%psi(state%triangles(:, e)))
                f = gs_profile_f(flux, f_edge, ffprime)
                if (abs(f) <= tiny(1.0_dp)) then
                    status = 1
                    return
                end if
                result%volume = result%volume + two_pi*w(node)*det*point(1)
                result%phi_edge = result%phi_edge + w(node)*det*f/point(1)
            end do
        end do
        call system_clock(tock)
        result%integral_s = real(tock - tick, dp)/clock_rate
        result%flux_consistency = abs(integral/result%phi_edge - 1)
        call gs_magnetic_axis(state, result%axis, status)
    end subroutine
end module
