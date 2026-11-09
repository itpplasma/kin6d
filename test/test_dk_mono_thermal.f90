!> Thermal coefficients from monoenergetic ones with the momentum correction
!! of Sugama and Nishimura, Phys. Plasmas 9, 4637 (2002) (mono_thermal_
!! coefficients), electrons and deuterium at 2 keV, against relations the
!! paper derives independently of the numerical solution:
!!   symmetric   for B = B(theta - N zeta) the energy integrals satisfy
!!               N_aj/M_aj = L_aj/N_aj = (G + N I)/(Z_a e (N q - 1) <B^2>)
!!               for j = 1, 2, 3 and both species (Eqs. (D4), (D5) with the
!!               poloidal flux per radian as radial label, chi' = 1, psi' = q);
!!               axisymmetric (N = 0) and quasi-helical (N = 1);
!!   ambipolar   in the symmetric case the particle fluxes are intrinsically
!!               ambipolar for arbitrary forces, Gamma_i = Gamma_e/Z_i
!!               (Appendix D);
!!   onsager     L_jk^ab = L_kj^ba and L_jE = -L_Ej (Eq. (C14)) on a field
!!               without symmetry.
program test_dk_mono_thermal
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use milca_dk_mono_thermal, only: mono_thermal_t, mono_thermal_coefficients, &
        THERMAL_OK
    use milca_dk_collop, only: collop_species_t, collop_grid_t, collop_make_grid, &
        collop_nuhat, collop_vth, collop_qe, COLLOP_OK
    implicit none
    real(dp), parameter :: pi = 4*atan(1.0_dp), kev = 1.602176634e-16_dp
    real(dp), parameter :: q = 2.2_dp, ig = 0.05_dp, b0 = 2.0_dp, rmajor = 1.7_dp
    character(64) :: check
    logical :: ok

    call get_command_argument(1, check)
    select case (trim(check))
    case ('symmetric')
        ok = symmetric(0)
        ok = symmetric(1) .and. ok
    case ('ambipolar'); ok = ambipolar()
    case ('onsager'); ok = onsager()
    case default
        print '(a)', 'unknown check '//trim(check)
        ok = .false.
    end select
    if (.not. ok) error stop 1

contains

    subroutine solve(b, nper, res, status)
        real(dp), intent(in) :: b(:, :)
        integer, intent(in) :: nper
        type(mono_thermal_t), intent(out) :: res
        integer, intent(out) :: status
        type(collop_species_t) :: sp(2)
        type(collop_grid_t) :: grid
        real(dp) :: lnl(2, 2)
        integer :: a
        sp(1) = collop_species_t(9.1093837015e-31_dp, -1.0_dp, 5e19_dp, 2*kev)
        sp(2) = collop_species_t(3.3435837768e-27_dp, 1.0_dp, 5e19_dp, 2*kev)
        ! NEO-2 collisionality 0.01/m for both species.
        do a = 1, 2
            lnl(a, :) = 3*sqrt(pi)/8*collop_vth(sp(a))*1e-2_dp/collop_nuhat(sp(a), sp(a), &
                1.0_dp)
        end do
        call collop_make_grid(5.0_dp, 4, 5, grid, status)
        if (status /= COLLOP_OK) return
        call mono_thermal_coefficients(b, nper, q, ig, 1/(q*rmajor), b0, sp, lnl, grid, &
            40, res, status)
    end subroutine solve

    !> b(theta - n zeta) on nt x nz points, zeta over one period (nper = max(n, 1)).
    function field(n, nt, nz) result(b)
        integer, intent(in) :: n, nt, nz
        real(dp) :: b(nt, nz)
        integer :: i, j
        do j = 1, nz
            do i = 1, nt
                b(i, j) = 1 - 0.1_dp*cos(2*pi*(i - 1)/nt - n*2*pi*(j - 1)/nz) &
                    + 0.02_dp*cos(2*(2*pi*(i - 1)/nt - n*2*pi*(j - 1)/nz))
            end do
        end do
    end function field

    logical function symmetric(n) result(ok)
        integer, intent(in) :: n
        type(mono_thermal_t) :: res
        real(dp) :: g, want, z(2)
        integer :: st, a, j
        call solve(field(n, 15, merge(1, 15, n == 0)), max(n, 1), res, st)
        ok = st == THERMAL_OK
        if (.not. ok) return
        ! kappa = iota B0/(G + iota I) with kappa = 1/(q R).
        g = b0/((1/(q*rmajor))*q*(1 + ig/q))
        z = [-1.0_dp, 1.0_dp]
        do a = 1, 2
            want = (g + n*ig*g)/(z(a)*collop_qe*(n*q - 1)*b0*b0*res%b2)
            do j = 1, 3
                ok = close('N/M', res%mnl(j, 2, a)/res%mnl(j, 1, a), want, 1e-8_dp) .and. ok
                ok = close('L/N', res%mnl(j, 3, a)/res%mnl(j, 2, a), want, 1e-8_dp) .and. ok
            end do
        end do
    end function symmetric

    logical function ambipolar() result(ok)
        type(mono_thermal_t) :: res
        integer :: st, j
        call solve(field(1, 15, 15), 1, res, st)
        ok = st == THERMAL_OK
        if (.not. ok) return
        do j = 1, 5
            ! The symmetric identities hold to 1e-9 (test symmetric); entries
            ! formed from their differences carry that error amplified.
            ok = close('Z Gamma_i = Gamma_e', res%l(3, j), res%l(1, j), 1e-7_dp) .and. ok
        end do
    end function ambipolar

    logical function onsager() result(ok)
        type(mono_thermal_t) :: res
        real(dp) :: b(15, 9), scale
        integer :: st, i, j
        do j = 1, 9
            do i = 1, 15
                b(i, j) = 1 - 0.1_dp*cos(2*pi*(i - 1)/15) &
                    - 0.04_dp*cos(2*pi*(i - 1)/15 - 2*pi*(j - 1)/9)
            end do
        end do
        call solve(b, 1, res, st)
        ok = st == THERMAL_OK
        if (.not. ok) return
        do i = 1, 5
            do j = 1, 5
                scale = sqrt(abs(res%l(i, i)*res%l(j, j)))
                if (i < 5 .and. j < 5) then
                    ok = ok .and. abs(res%l(i, j) - res%l(j, i)) <= 1e-12_dp*scale
                else if (i /= j) then
                    ok = ok .and. abs(res%l(i, j) + res%l(j, i)) <= 1e-12_dp*scale
                end if
            end do
        end do
        print '(a,l2)', 'Onsager symmetry of the transport matrix:', ok
    end function onsager

    logical function close(label, got, want, rtol)
        character(*), intent(in) :: label
        real(dp), intent(in) :: got, want, rtol
        close = abs(got - want) <= rtol*abs(want)
        print '(a,2es20.11,a,es9.2,a,l2)', label//' ', got, want, ' rel ', &
            abs(got - want)/abs(want), ' pass', close
    end function close
end program test_dk_mono_thermal
