! Frozen NN1/NN2 Solov'ev oracle (docs/equilibrium-nonnested.md).
module nn_solovev_fixture
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    real(dp), parameter :: r0 = 6.2_dp, minor = 0.62_dp, ftor = 32.86_dp
    real(dp), parameter :: c = ftor/(12*r0**3), pp = -16*c
contains
    subroutine mapping(u, theta, x)
        real(dp), intent(in) :: u, theta
        real(dp), intent(out) :: x(2)
        x(1) = sqrt(r0*r0 + 2*r0*minor*u*cos(theta))
        x(2) = r0*minor*u*sin(theta)/x(1)
    end subroutine
    subroutine coordinates(centre, scale, x, u, theta)
        real(dp), intent(in) :: centre(2), scale, x(2)
        real(dp), intent(out) :: u, theta
        real(dp) :: v(2)
        v = [(x(1)**2 - centre(1)**2)/2, x(1)*x(2)]/(centre(1)*scale)
        u = norm2(v); theta = atan2(v(2), v(1))
    end subroutine
    subroutine boundary(t, x)
        real(dp), intent(in) :: t
        real(dp), intent(out) :: x(2)
        call mapping(1.0_dp, 2*acos(-1.0_dp)*t, x)
    end subroutine
    subroutine exact(x, psi, b)
        real(dp), intent(in) :: x(2)
        real(dp), intent(out) :: psi, b(2)
        real(dp) :: r, z
        r = x(1); z = x(2)
        psi = c*((r*r - r0*r0)**2 + 4*r*r*z*z - 4*r0*r0*minor*minor)
        b = [-8*c*r*z, 4*c*(r*r - r0*r0 + 2*z*z)]
    end subroutine
    pure real(dp) function laplacian(x)
        ! Cylindrical Laplacian of psi: Delta* psi + 2 psi_R / R.
        real(dp), intent(in) :: x(2)
        laplacian = c*(24*x(1)**2 - 8*r0*r0 + 16*x(2)**2)
    end function
end module

