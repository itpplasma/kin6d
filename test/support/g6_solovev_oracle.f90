! Independent flux integrals for the frozen Solov'ev fixture, mu0 = 1.
! With X=(R^2-R0^2)/2, Y=R Z, the disk map has
! dR dZ = R0^2 a^2 u/R^2 du dtheta and psi=4 c R0^2 a^2 (u^2-1).
module g6_solovev_oracle
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_quadrature, only: gauss_legendre_ab
    use nn_solovev_fixture, only: r0, minor, ftor, c
    implicit none
    private
    public :: exact_flux, exact_iota
contains
    function exact_flux(radius) result(flux)
        real(dp), intent(in) :: radius
        real(dp) :: flux, u(32), w(32), theta, r
        integer :: i, j
        call gauss_legendre_ab(32, 0.0_dp, radius, u, w)
        flux = 0
        do j = 1, 256
            theta = 2*acos(-1.0_dp)*(j - 0.5_dp)/256
            do i = 1, size(u)
                r = sqrt(r0*r0 + 2*r0*minor*u(i)*cos(theta))
                flux = flux + w(i)*u(i)/r**3
            end do
        end do
        flux = flux*(2*acos(-1.0_dp)/256)*ftor*r0*r0*minor*minor
    end function

    function exact_iota(radius) result(iota)
        real(dp), intent(in) :: radius
        real(dp) :: iota, theta, r, q
        integer :: j
        q = 0
        do j = 1, 256
            theta = 2*acos(-1.0_dp)*(j - 0.5_dp)/256
            r = sqrt(r0*r0 + 2*r0*minor*radius*cos(theta))
            q = q + ftor/(8*c*r**3)/256
        end do
        iota = 1/q
    end function
end module
