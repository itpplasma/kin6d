!> Local linear roots of the fluid drift-wave models of
!> docs/drift-turbulence.md for one wave vector, exp[i(kx x + ky y - omega t)].
!> The formulas are the generated leaves of gen/gen_drift_fluid.f90; this
!> module selects branches, solves the drift-Alfven cubic and applies the
!> implicit rule for root tangents.
module milca_drift_dispersion
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use milca_drift_fluid_leaves, only: hw_coupling, hw_coupling_jvp, &
        drift_hyperviscosity, drift_hyperviscosity_jvp, hw_roots, hw_roots_jvp, &
        hw_density_response, drift_alfven_cubic, drift_alfven_cubic_jvp
    implicit none
    private
    public :: hw_parameters_t, hw_linear_roots, hw_linear_roots_jvp
    public :: drift_alfven_parameters_t, drift_alfven_roots, &
        drift_alfven_root_tangent
    public :: drift_root_ok, drift_root_decoupled, drift_root_invalid

    integer, parameter :: drift_root_ok = 0
    !> No resistive coupling acts on this mode: both roots are -i nu and the
    !> density response is undefined.
    integer, parameter :: drift_root_decoupled = 1
    integer, parameter :: drift_root_invalid = 2

    !> Hasegawa-Wakatani parameters: adiabaticity, kappa = rho_s/L_n,
    !> hyperviscosity d (-laplace)^order, and the modified model in which the
    !> coupling acts on ky /= 0 only.
    type :: hw_parameters_t
        real(dp) :: alpha = 1, kappa = 1, d = 0
        integer :: order = 2
        logical :: modified = .false.
    end type

    !> beta = (beta_e/2)(L_par/L_n)^2, mu = (m_e/m_i)(L_par/L_n)^2,
    !> c = 0.51 (nu_e L_n/c_s) mu.
    type :: drift_alfven_parameters_t
        real(dp) :: beta = 0, mu = 0, c = 0
    end type

    interface
        subroutine zgeev(jobvl, jobvr, n, a, lda, w, vl, ldvl, vr, ldvr, work, &
                         lwork, rwork, info)
            import :: dp
            character, intent(in) :: jobvl, jobvr
            integer, intent(in) :: n, lda, ldvl, ldvr, lwork
            complex(dp), intent(inout) :: a(lda, *)
            complex(dp), intent(out) :: w(*), vl(ldvl, *), vr(ldvr, *), work(*)
            real(dp), intent(out) :: rwork(*)
            integer, intent(out) :: info
        end subroutine
    end interface
contains
    !> omega(1) is the root of larger growth rate; response = n/phi of each.
    subroutine hw_linear_roots(par, kx, ky, omega, response, status)
        type(hw_parameters_t), intent(in) :: par
        real(dp), intent(in) :: kx, ky
        complex(dp), intent(out) :: omega(2), response(2)
        integer, intent(out) :: status
        real(dp) :: b, wstar, nu, g(2), w(2), nr, ni
        integer :: j

        omega = 0; response = 0
        status = drift_root_invalid
        if (kx == 0 .and. ky == 0) return
        if (par%alpha < 0 .or. par%d < 0) return
        call drift_hyperviscosity(par%d, real(par%order, dp), kx, ky, nu)
        if (decoupled(par, ky)) then
            omega = cmplx(0, -nu, dp)
            status = drift_root_decoupled
            return
        end if
        call hw_coupling(par%alpha, par%kappa, kx, ky, b, wstar)
        call hw_roots(b, wstar, nu, g(1), w(1), g(2), w(2))
        do j = 1, 2
            call hw_density_response(par%alpha, kx, ky, nu, g(j), w(j), nr, ni)
            omega(j) = cmplx(w(j), g(j), dp)
            response(j) = cmplx(nr, ni, dp)
        end do
        status = drift_root_ok
    end subroutine

    !> Tangent of both roots along (dalpha, dkappa, dd, dkx, dky).
    subroutine hw_linear_roots_jvp(par, kx, ky, dalpha, dkappa, dd, dkx, dky, &
                                   domega, status)
        type(hw_parameters_t), intent(in) :: par
        real(dp), intent(in) :: kx, ky, dalpha, dkappa, dd, dkx, dky
        complex(dp), intent(out) :: domega(2)
        integer, intent(out) :: status
        real(dp) :: b, wstar, nu, db, dwstar, dnu, dg(2), dw(2)

        domega = 0
        status = drift_root_invalid
        if (kx == 0 .and. ky == 0) return
        if (par%alpha < 0 .or. par%d < 0) return
        call drift_hyperviscosity(par%d, real(par%order, dp), kx, ky, nu)
        call drift_hyperviscosity_jvp(par%d, real(par%order, dp), kx, ky, dd, &
                                      0.0_dp, dkx, dky, dnu)
        if (decoupled(par, ky)) then
            domega = cmplx(0, -dnu, dp)
            status = drift_root_decoupled
            return
        end if
        call hw_coupling(par%alpha, par%kappa, kx, ky, b, wstar)
        call hw_coupling_jvp(par%alpha, par%kappa, kx, ky, dalpha, dkappa, dkx, &
                             dky, db, dwstar)
        call hw_roots_jvp(b, wstar, nu, db, dwstar, dnu, dg(1), dw(1), dg(2), dw(2))
        domega = cmplx(dw, dg, dp)
        status = drift_root_ok
    end subroutine

    pure logical function decoupled(par, ky)
        type(hw_parameters_t), intent(in) :: par
        real(dp), intent(in) :: ky
        decoupled = par%alpha == 0 .or. (par%modified .and. ky == 0)
    end function

    !> The three roots ordered by decreasing growth rate. residual is
    !> max |P(omega)| over the roots divided by the sum of |a_j omega^j|.
    subroutine drift_alfven_roots(par, kpar, kx, ky, omega, residual, status)
        type(drift_alfven_parameters_t), intent(in) :: par
        real(dp), intent(in) :: kpar, kx, ky
        complex(dp), intent(out) :: omega(3)
        real(dp), intent(out) :: residual
        integer, intent(out) :: status
        complex(dp) :: a(0:3), companion(3, 3), vl(1, 1), vr(1, 1), work(12), z
        real(dp) :: rwork(6), scale
        integer :: info, j, k

        omega = 0; residual = huge(1.0_dp)
        status = drift_root_invalid
        call coefficients(par, kpar, kx, ky, a)
        if (abs(a(3)) == 0) return
        companion = 0
        companion(2, 1) = 1; companion(3, 2) = 1
        companion(:, 3) = -a(0:2)/a(3)
        call zgeev('N', 'N', 3, companion, 3, omega, vl, 1, vr, 1, work, &
                   size(work), rwork, info)
        if (info /= 0) return
        do j = 1, 3
            call polish(a, omega(j))
        end do
        do j = 2, 3
            z = omega(j)
            k = j - 1
            do while (k >= 1)
                if (aimag(omega(k)) >= aimag(z)) exit
                omega(k + 1) = omega(k)
                k = k - 1
            end do
            omega(k + 1) = z
        end do
        residual = 0
        do j = 1, 3
            z = omega(j)
            scale = abs(a(3)*z**3) + abs(a(2)*z**2) + abs(a(1)*z) + abs(a(0))
            if (scale > 0) residual = max(residual, abs(polynomial(a, z))/scale)
        end do
        status = drift_root_ok
    end subroutine

    !> Implicit rule d omega = -dP/P'(omega) at a simple root, along
    !> (dbeta, dmu, dc, dkpar, dkx, dky).
    subroutine drift_alfven_root_tangent(par, kpar, kx, ky, omega, dpar, dkpar, &
                                         dkx, dky, domega, status)
        type(drift_alfven_parameters_t), intent(in) :: par, dpar
        real(dp), intent(in) :: kpar, kx, ky, dkpar, dkx, dky
        complex(dp), intent(in) :: omega
        complex(dp), intent(out) :: domega
        integer, intent(out) :: status
        complex(dp) :: a(0:3), da(0:3), slope
        real(dp) :: d3, d2r, d2i, d1, d0

        call coefficients(par, kpar, kx, ky, a)
        call drift_alfven_cubic_jvp(par%beta, par%mu, par%c, kpar, kx, ky, &
            dpar%beta, dpar%mu, dpar%c, dkpar, dkx, dky, d3, d2r, d2i, d1, d0)
        da = [cmplx(d0, 0, dp), cmplx(d1, 0, dp), cmplx(d2r, d2i, dp), &
              cmplx(d3, 0, dp)]
        slope = 3*a(3)*omega**2 + 2*a(2)*omega + a(1)
        domega = 0
        status = drift_root_invalid
        if (abs(slope) == 0) return
        domega = -polynomial(da, omega)/slope
        status = drift_root_ok
    end subroutine

    subroutine coefficients(par, kpar, kx, ky, a)
        type(drift_alfven_parameters_t), intent(in) :: par
        real(dp), intent(in) :: kpar, kx, ky
        complex(dp), intent(out) :: a(0:3)
        real(dp) :: a3, a2r, a2i, a1, a0
        call drift_alfven_cubic(par%beta, par%mu, par%c, kpar, kx, ky, a3, a2r, &
                                a2i, a1, a0)
        a = [cmplx(a0, 0, dp), cmplx(a1, 0, dp), cmplx(a2r, a2i, dp), &
             cmplx(a3, 0, dp)]
    end subroutine

    pure complex(dp) function polynomial(a, z)
        complex(dp), intent(in) :: a(0:3), z
        polynomial = ((a(3)*z + a(2))*z + a(1))*z + a(0)
    end function

    !> Two Newton steps on the eigenvalue of the companion matrix.
    pure subroutine polish(a, z)
        complex(dp), intent(in) :: a(0:3)
        complex(dp), intent(inout) :: z
        complex(dp) :: slope
        integer :: it
        do it = 1, 2
            slope = (3*a(3)*z + 2*a(2))*z + a(1)
            if (abs(slope) == 0) return
            z = z - polynomial(a, z)/slope
        end do
    end subroutine
end module milca_drift_dispersion
