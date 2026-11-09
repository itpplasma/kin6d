!> Axisymmetric transport matrix of the energy-coupled kernel
!! (ntv_axisymmetric_coefficients) on the circular surface b = 1 - eps cos(theta)
!! in Boozer angles, b.grad theta = k b:
!!   onsager        d(1, 1, 1, 3) = -d(1, 1, 3, 1) and d(1, 1, 2, 3) =
!!                  -d(1, 1, 3, 2), with pitch-angle scattering and with the
!!                  full linearised operator;
!!   momentum       a single species with momentum-conserving self-collisions:
!!                  A_1 drives only the rigid rotation, so the flow response
!!                  is -(T/(Z e))/(2 k q) exactly (D_31 = G T/(Z e));
!!   local          the local flow response, whose Boozer surface average is
!!                  the flow coefficient, is a + c b^2 over the surface
!!                  (u_par B = omega G + K B^2, from particle conservation);
!!   collisionless  with pitch-angle scattering every speed node solves the
!!                  problem of dk_mono, driven by (T/(Z e)) x r times its radial
!!                  drift, r = -2/(k q). As nu -> 0 the parallel-flow response
!!                  to A_1 therefore tends to (2 T r/(Z e sqrt(pi))) D31_mono
!!                  int x^4 e^-x^2 dx = (3/4) (T/(Z e)) r D31_mono at the same
!!                  theta and Legendre resolution, and the response to A_2 to
!!                  int x^6 e^-x^2 dx/int x^4 e^-x^2 dx = 5/2 times that. The
!!                  sqrt(nu) corrections are removed by two collisionalities.
!!                  (D31_mono -> f_t/3 needs far more Legendre modes; dk_mono
!!                  checks that limit.)
program test_ntv_axisymmetric
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use milca_dk_ntv, only: ntv_surface_t, ntv_axisymmetric_coefficients, NTV_OK, &
        NTV_WEIGHT_FULL
    use milca_dk_mono, only: mono_coefficients_t, mono_coefficients, MONO_OK
    use milca_dk_collop, only: collop_species_t, collop_grid_t, collop_make_grid, &
        collop_nuhat, collop_vth, COLLOP_OK, COLLOP_LORENTZ, COLLOP_FULL
    implicit none
    real(dp), parameter :: pi = 4*atan(1.0_dp), kev = 1.602176634e-16_dp
    real(dp), parameter :: charge = 1.602176634e-19_dp
    real(dp), parameter :: eps = 0.1_dp, q = 2.2_dp, rmajor = 6.2_dp
    integer, parameter :: n_theta = 33, n_leg = 40
    character(64) :: check
    logical :: ok

    call get_command_argument(1, check)
    select case (trim(check))
    case ('onsager'); ok = onsager()
    case ('momentum'); ok = momentum()
    case ('local'); ok = local_flow()
    case ('collisionless'); ok = collisionless()
    case default
        print '(a)', 'unknown check '//trim(check)
        ok = .false.
    end select
    if (.not. ok) error stop 1

contains

    function surface() result(sf)
        type(ntv_surface_t) :: sf
        sf%eps = eps
        sf%q = q
        sf%k = 1/(q*rmajor)
        ! One harmonic: the field setup needs it; the axisymmetric solve
        ! drops the perturbation.
        sf%n = 3
        sf%m = [6]
        sf%c = [(1e-3_dp, 0.0_dp)]
        sf%weight = NTV_WEIGHT_FULL
        sf%boozer_stream = .true.
    end function surface

    !> Deuterium at 2 keV; kappa is the NEO-2 collisionality 2/(v tau) [1/m].
    subroutine solve(kappa, model, d, status)
        real(dp), intent(in) :: kappa
        integer, intent(in) :: model
        real(dp), intent(out) :: d(1, 1, 3, 3)
        integer, intent(out) :: status
        type(collop_species_t) :: sp(1)
        type(collop_grid_t) :: grid
        real(dp) :: lnl(1, 1)
        integer :: its
        sp(1) = collop_species_t(3.3435837768e-27_dp, 1.0_dp, 3e19_dp, 2*kev)
        lnl = 3*sqrt(pi)/8*collop_vth(sp(1))*kappa/collop_nuhat(sp(1), sp(1), 1.0_dp)
        call collop_make_grid(5.0_dp, 6, 5, grid, status)
        if (status /= COLLOP_OK) return
        call ntv_axisymmetric_coefficients(surface(), sp, lnl, grid, model, n_theta, n_leg, &
            d, its, status)
    end subroutine solve

    logical function onsager() result(ok)
        real(dp) :: d(1, 1, 3, 3)
        integer :: status, m
        integer, parameter :: models(2) = [COLLOP_LORENTZ, COLLOP_FULL]
        ok = .true.
        do m = 1, 2
            call solve(1e-3_dp, models(m), d, status)
            if (status /= NTV_OK) then
                print '(a,i0)', 'solve failed with status ', status
                ok = .false.
                return
            end if
            ok = ok .and. close('d13 = -d31', d(1, 1, 1, 3), -d(1, 1, 3, 1), 1e-6_dp) &
                .and. close('d23 = -d32', d(1, 1, 2, 3), -d(1, 1, 3, 2), 1e-6_dp)
        end do
    end function onsager

    !> Fit a + c b^2 to the local A_1 and A_2 flow responses with the full
    !> operator; the residual and the average must close. Particle
    !> conservation of the discrete operator holds to the speed quadrature
    !> (2e-7 for A_2), hence the tolerance.
    logical function local_flow() result(ok)
        type(collop_species_t) :: sp(1)
        type(collop_grid_t) :: grid
        real(dp) :: lnl(1, 1), d(1, 1, 3, 3), loc(n_theta, 1, 1, 3), theta(n_theta), &
            b(n_theta), w(n_theta), m(n_theta, 2), coef(2), resid
        integer :: status, its, i, j
        sp(1) = collop_species_t(3.3435837768e-27_dp, 1.0_dp, 3e19_dp, 2*kev)
        lnl = 3*sqrt(pi)/8*collop_vth(sp(1))*1e-3_dp/collop_nuhat(sp(1), sp(1), 1.0_dp)
        call collop_make_grid(5.0_dp, 6, 5, grid, status)
        call ntv_axisymmetric_coefficients(surface(), sp, lnl, grid, COLLOP_FULL, n_theta, &
            n_leg, d, its, status, local=loc)
        ok = status == NTV_OK
        if (.not. ok) return
        theta = [(2*pi*(i - 1)/n_theta, i=1, n_theta)]
        b = 1 - eps*cos(theta)
        w = 1/b**2
        w = w/sum(w)
        do j = 1, 2
            m(:, 1) = 1
            m(:, 2) = b**2
            coef = least_squares(m, loc(:, 1, 1, j))
            resid = maxval(abs(matmul(m, coef) - loc(:, 1, 1, j)))/maxval(abs(loc(:, 1, 1, j)))
            ok = ok .and. close('local flow average', sum(w*loc(:, 1, 1, j)), d(1, 1, 3, j), &
                1e-12_dp) .and. resid < 1e-6_dp
            print '(a,es10.2)', 'residual of a + c b^2', resid
        end do
    end function local_flow

    function least_squares(m, y) result(c)
        real(dp), intent(in) :: m(:, :), y(:)
        real(dp) :: c(2), a(2, 2), r(2), det
        a = matmul(transpose(m), m)
        r = matmul(transpose(m), y)
        det = a(1, 1)*a(2, 2) - a(1, 2)*a(2, 1)
        c = [a(2, 2)*r(1) - a(1, 2)*r(2), a(1, 1)*r(2) - a(2, 1)*r(1)]/det
    end function least_squares

    logical function momentum() result(ok)
        real(dp) :: d(1, 1, 3, 3)
        integer :: status
        ok = .false.
        call solve(1e-3_dp, COLLOP_FULL, d, status)
        if (status /= NTV_OK) return
        ok = close('A_1 flow, full operator', d(1, 1, 3, 1), &
            -(2*kev/charge)/(2*(1/(q*rmajor))*q), 1e-5_dp)
    end function momentum

    logical function collisionless() result(ok)
        real(dp), parameter :: nu(2) = [6.25e-6_dp, 1.5625e-6_dp], mu(2) = [1e-6_dp, 2.5e-7_dp]
        type(mono_coefficients_t) :: c(2)
        real(dp) :: d(1, 1, 3, 3, 2), theta(n_theta), b(n_theta), k, m31, d31, d32
        integer :: status, i
        theta = [(2*pi*(i - 1)/n_theta, i=1, n_theta)]
        b = 1 - eps*cos(theta)
        k = 1/(q*rmajor)
        ok = .false.
        do i = 1, 2
            call mono_coefficients(b, k, mu(i), n_leg, c(i), status)
            if (status /= MONO_OK) return
            call solve(nu(i), COLLOP_LORENTZ, d(:, :, :, :, i), status)
            if (status /= NTV_OK) return
        end do
        m31 = richardson(c(1)%d31, mu(1), c(2)%d31, mu(2))
        d31 = richardson(d(1, 1, 3, 1, 1), nu(1), d(1, 1, 3, 1, 2), nu(2))
        d32 = richardson(d(1, 1, 3, 2, 1), nu(1), d(1, 1, 3, 2, 2), nu(2))
        ok = close('A_1 flow against dk_mono', d31, &
            0.75_dp*(2*kev/charge)*(-2/(k*q))*m31, 1e-3_dp) &
            .and. close('A_2 flow / A_1 flow', d32/d31, 2.5_dp, 1e-3_dp)
    end function collisionless

    !> Richardson extrapolation to nu -> 0 of a sqrt(nu) correction.
    real(dp) function richardson(a1, nu1, a2, nu2)
        real(dp), intent(in) :: a1, nu1, a2, nu2
        richardson = (a2*sqrt(nu1) - a1*sqrt(nu2))/(sqrt(nu1) - sqrt(nu2))
    end function richardson

    logical function close(label, got, want, rtol)
        character(*), intent(in) :: label
        real(dp), intent(in) :: got, want, rtol
        close = abs(got - want) <= rtol*abs(want)
        print '(a,2es20.11,a,es9.2,a,l2)', label//' ', got, want, ' rel ', &
            abs(got - want)/abs(want), ' pass', close
    end function close
end program test_ntv_axisymmetric
