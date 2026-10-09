module kin6d_gs_readback
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use kin6d_grad_shafranov, only: gs_state_t, gs_cell_data, gs_axis_candidate, &
        gs_profile_integral
    use kin6d_gs_coarea, only: gs_level_coarea
    use kin6d_gs_axis_recovery, only: gs_magnetic_axis
    use fortfem_triangle_duffy_quadrature, only: triangle_duffy_quadrature
    implicit none
    private
    public :: gs_readback_t, read_gs_equilibrium, gs_profile_f
    type :: gs_readback_t
        real(dp) :: axis(2) = 0, psi_axis = 0, volume = 0, phi_edge = 0
        real(dp) :: flux_consistency = 0
        real(dp), allocatable :: psi(:), s_tor(:), q(:)
    end type
contains
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
        real(dp), allocatable :: coarea(:), x(:), y(:), w(:), density(:)
        real(dp) :: hessian(2, 2), shape(6), gradients(2, 6), hessians(2, 2, 6)
        real(dp) :: point(2), det, flux, f, two_pi, integral
        integer :: axis_element, i, e, node
        status = 1
        if (.not. state%ready .or. state%element_order /= 2) return
        if (radial_points < 5 .or. radial_points > 4097) return
        if (.not. all(ieee_is_finite([f_edge, ffprime]))) return
        if (abs(f_edge) <= tiny(1.0_dp)) return
        if (any(abs(pack(state%psi, state%boundary)) > 1.0e-12_dp)) return
        call gs_axis_candidate(state, result%axis, result%psi_axis, hessian, &
            axis_element, status)
        if (status /= 0) return
        status = 1
        if (result%psi_axis >= 0) return
        if (f_edge**2 + 2*ffprime*result%psi_axis <= 0) return
        two_pi = 2*acos(-1.0_dp)
        allocate (result%psi(radial_points), result%s_tor(radial_points), &
            result%q(radial_points), coarea(radial_points), density(radial_points))
        do i = 1, radial_points
            result%psi(i) = result%psi_axis &
                *(1 - (real(i - 1, dp)/(radial_points - 1))**2)
        end do
        call gs_level_coarea(state, result%psi, angular_points, coarea, &
            result%axis, status)
        if (status /= 0) return
        do i = 1, radial_points
            f = gs_profile_f(result%psi(i), f_edge, ffprime)
            density(i) = f*coarea(i)
        end do
        result%q = density/two_pi
        result%s_tor(1) = 0
        do i = 2, radial_points
            result%s_tor(i) = result%s_tor(i - 1) &
                + (result%psi(i) - result%psi(i - 1))*(density(i) + density(i - 1))/2
        end do
        integral = result%s_tor(radial_points)
        if (abs(integral) <= tiny(1.0_dp)) then
            status = 1
            return
        end if
        result%s_tor = result%s_tor/integral
        call triangle_duffy_quadrature(10, x, y, w, status)
        if (status /= 0) return
        do e = 1, size(state%triangles, 2)
            do node = 1, size(w)
                call gs_cell_data(state, e, [x(node), y(node)], point, det, shape, &
                    gradients, hessians, status)
                if (status /= 0) return
                flux = dot_product(shape, state%psi(state%triangles(:, e)))
                f = gs_profile_f(flux, f_edge, ffprime)
                if (abs(f) <= tiny(1.0_dp)) then
                    status = 1
                    return
                end if
                result%volume = result%volume + two_pi*w(node)*det*point(1)
                result%phi_edge = result%phi_edge + w(node)*det*f/point(1)
            end do
        end do
        result%flux_consistency = abs(integral/result%phi_edge - 1)
        call gs_magnetic_axis(state, result%axis, status)
    end subroutine
end module
