!> Independent linear-response oracle for Q7 part 2: dielectric function of
!> drifting Maxwellian electron beams under the Wigner (hbar > 0) or Vlasov
!> (hbar = 0) equation with a neutralizing background, in the units of
!> kin6d_wigner_leaves (omega_p = 1, m = 1, v/v_ref). The plasma dispersion
!> function is integrated from its ODE Z' = -2(1 + zeta Z), Z(0) = i sqrt(pi),
!> which yields the Landau continuation for Im(zeta) < 0 automatically.
!> This module shares no formula with the split-step solver.
module kin6d_wigner_dispersion
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    private
    public :: beam_t, plasma_z, dielectric, dispersion_root

    real(dp), parameter :: pi = acos(-1.0_dp)
    complex(dp), parameter :: iu = (0.0_dp, 1.0_dp)

    type :: beam_t
        real(dp) :: weight = 1.0_dp, drift = 0.0_dp, sigma = 1.0_dp
    end type

contains

    !> Z(zeta) and Z'(zeta) by RK4 along the straight path from 0.
    pure subroutine plasma_z(zeta, z, dz)
        complex(dp), intent(in) :: zeta
        complex(dp), intent(out) :: z, dz
        complex(dp) :: k1, k2, k3, k4
        real(dp) :: h, t
        integer :: steps, j
        steps = 400 + ceiling(4000*abs(zeta))
        h = 1.0_dp/steps
        z = iu*sqrt(pi)
        do j = 0, steps - 1
            t = j*h
            k1 = rhs(t, z)
            k2 = rhs(t + h/2, z + h*k1/2)
            k3 = rhs(t + h/2, z + h*k2/2)
            k4 = rhs(t + h, z + h*k3)
            z = z + h*(k1 + 2*k2 + 2*k3 + k4)/6
        end do
        dz = -2*(1 + zeta*z)
    contains
        pure complex(dp) function rhs(s, y)
            real(dp), intent(in) :: s
            complex(dp), intent(in) :: y
            rhs = -2*zeta*(1 + s*zeta*y)
        end function
    end subroutine

    !> eps(k, omega) and its omega- and hbar-derivatives.
    pure subroutine dielectric(k, omega, hbar, beams, eps, deps_domega, deps_dhbar)
        real(dp), intent(in) :: k, hbar
        complex(dp), intent(in) :: omega
        type(beam_t), intent(in) :: beams(:)
        complex(dp), intent(out) :: eps, deps_domega, deps_dhbar
        complex(dp) :: c, xp, xm, zp, zm, dzp, dzm, x0, z0, dz0, d2z0
        real(dp) :: a, scale, delta
        integer :: b
        c = omega/k
        eps = 1; deps_domega = 0; deps_dhbar = 0
        do b = 1, size(beams)
            scale = sqrt(2.0_dp)*beams(b)%sigma
            if (hbar > 1e-6_dp) then
                delta = hbar*k/2
                xp = (c + delta - beams(b)%drift)/scale
                xm = (c - delta - beams(b)%drift)/scale
                call plasma_z(xp, zp, dzp)
                call plasma_z(xm, zm, dzm)
                a = beams(b)%weight/(k**3*scale)
                eps = eps - a*(zp - zm)/hbar
                deps_domega = deps_domega - a*(dzp - dzm)/(hbar*k*scale)
                deps_dhbar = deps_dhbar + a*(zp - zm)/hbar**2 &
                    - a*(dzp + dzm)*k/(2*scale*hbar)
            else
                x0 = (c - beams(b)%drift)/scale
                call plasma_z(x0, z0, dz0)
                d2z0 = -2*(z0 + x0*dz0)
                a = beams(b)%weight/(2*k**2*beams(b)%sigma**2)
                eps = eps - a*dz0
                deps_domega = deps_domega - a*d2z0/(k*scale)
            end if
        end do
    end subroutine

    !> Newton root of eps(k, omega) = 0 from a guess; domega_dhbar is the
    !> equation-level implicit derivative -eps_hbar/eps_omega at the root.
    pure subroutine dispersion_root(k, hbar, beams, guess, omega, domega_dhbar, ok)
        real(dp), intent(in) :: k, hbar
        type(beam_t), intent(in) :: beams(:)
        complex(dp), intent(in) :: guess
        complex(dp), intent(out) :: omega, domega_dhbar
        logical, intent(out) :: ok
        complex(dp) :: eps, de, dh, step
        integer :: iteration
        omega = guess
        ok = .false.
        do iteration = 1, 60
            call dielectric(k, omega, hbar, beams, eps, de, dh)
            step = eps/de
            omega = omega - step
            if (abs(step) < 1e-13_dp*max(1.0_dp, abs(omega))) then
                ok = .true.
                exit
            end if
        end do
        call dielectric(k, omega, hbar, beams, eps, de, dh)
        domega_dhbar = -dh/de
    end subroutine
end module kin6d_wigner_dispersion
