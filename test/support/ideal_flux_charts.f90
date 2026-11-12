module ideal_flux_charts
    ! Analytic test equilibria for milca_ideal_flux.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use milca_ideal_flux, only: flux_equilibrium_t, flux_point_t
    implicit none
    private
    public :: cylinder_t
    !> Straight periodic cylinder x = (a rho cos t, a rho sin t, R0 zeta).
    type, extends(flux_equilibrium_t) :: cylinder_t
        real(dp) :: a = 1, r0 = 1
        integer :: profile = 1
    contains
        procedure :: evaluate => cylinder
    end type cylinder_t
contains
    subroutine cylinder(self, rho, theta, zeta, point)
        class(cylinder_t), intent(in) :: self
        real(dp), intent(in) :: rho, theta, zeta
        type(flux_point_t), intent(out) :: point
        real(dp) :: r, bz, bzp, p, pp, c, s, cyl_a, cyl_r0
        cyl_a = self%a; cyl_r0 = self%r0
        r = cyl_a*rho; c = cos(theta + 0*zeta); s = sin(theta)
        if (self%profile == 1) then
            bz = sqrt(0.25_dp + r*r); bzp = r/bz
            p = 1.5_dp*(cyl_a**2 - r*r); pp = -3*r
        else
            bz = 1; bzp = 0; p = 0.3_dp; pp = 0
        end if
        point%xu(:, 1) = [cyl_a*c, cyl_a*s, 0.0_dp]
        point%xu(:, 2) = [-r*s, r*c, 0.0_dp]
        point%xu(:, 3) = [0.0_dp, 0.0_dp, cyl_r0]
        point%xuu = 0
        point%xuu(:, 1, 2) = [-cyl_a*s, cyl_a*c, 0.0_dp]
        point%xuu(:, 2, 1) = point%xuu(:, 1, 2)
        point%xuu(:, 2, 2) = [-r*c, -r*s, 0.0_dp]
        ! psi' = a r Bz, iota = R0 Btheta/(r Bz) in rho; zeta = z/R0.
        point%psi1 = cyl_a*r*bz
        point%psi2 = cyl_a**2*(bz + r*bzp)
        if (self%profile == 1) then
            point%iota = cyl_r0/bz
            point%iota1 = -cyl_a*cyl_r0*bzp/bz**2
        else
            point%iota = 0; point%iota1 = 0
        end if
        point%lam = 0
        point%p = p; point%pp = cyl_a*pp; point%density = 1
    end subroutine cylinder
end module ideal_flux_charts
