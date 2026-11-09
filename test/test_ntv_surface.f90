!> Transport matrix of the energy-coupled operator on a 3D surface in Boozer
!! angles (ntv_surface_coefficients), deuterium at 2 keV, against oracles that
!! do not share its 3D assembly:
!!   axisymmetric   b(theta) on a theta x zeta grid reproduces the axisymmetric
!!                  matrix (ntv_axisymmetric_coefficients, analytic leaf) with
!!                  the full linearised operator;
!!   quasisymmetric b = b(theta - N zeta) is the axisymmetric problem in
!!                  chi = theta - N zeta: streaming kappa (1 - q N) b d_chi
!!                  and radial drive (1 + N I/G) b'/q, so radial rows and
!!                  columns scale with 1 + N I/G (full operator);
!!   monoenergetic  with pitch-angle scattering on a non-symmetric field every
!!                  speed node x_k solves mono_coefficients_3d at
!!                  nu = nu_D(x_k)/(v_th x_k), with nu_D of Chandrasekhar, so
!!                  d11 = 8 T^2/(Z^2 e^2 sqrt(pi) v_th kappa^2 q^2)
!!                        sum_k w4_k^2/w3_k D11(nu_k),
!!                  d31 = -4 T/(Z e sqrt(pi) kappa q) sum_k w4_k D31(nu_k),
!!                  d33 = 2 v_th/sqrt(pi) sum_k w3_k D33(nu_k),
!!                  w_p = w_k x_k^p exp(-x_k^2) the GLL weights;
!!   onsager        d13 = -d31, d23 = -d32 on the non-symmetric field with
!!                  the full operator, at a collisionality the grid resolves
!!                  (the discretisation defect falls spectrally with the
!!                  theta x zeta grid), to the GMRES tolerance;
!!   ambipolar      intrinsic ambipolarity of a quasi-helical field: with
!!                  momentum-conserving self-collisions the particle flux
!!                  vanishes, d11 and d12 -> 0 as the speed grid is refined
!!                  (Helander and Simakov, PRL 101, 145003, 2008), here
!!                  against the pitch-angle d11 on the same grid;
!!   iterative      GMRES right-preconditioned by the collision-mode
!!                  factorisations, without (pairs = 0) and with (pairs = 1)
!!                  the exact l <= 1 block, reproduces the direct solve;
!!   exact          with pitch-angle scattering every preconditioner is the
!!                  inverse, so each GMRES solve converges in one iteration:
!!                  the modal chain with 0 and 1 exact pair blocks here, and
!!                  the per-node chain of the NTV kernel
!!                  (ntv_energy_coefficients, two exact pair blocks).
program test_ntv_surface
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use milca_dk_ntv, only: ntv_surface_t, ntv_axisymmetric_coefficients, &
                            ntv_surface_coefficients, ntv_energy_coefficients, NTV_OK, &
                            NTV_WEIGHT_FULL
    use milca_dk_mono, only: mono_coefficients_t, mono_coefficients_3d, MONO_OK
    use milca_dk_collop, only: collop_species_t, collop_grid_t, collop_make_grid, &
                        collop_nuhat, collop_vth, COLLOP_OK, COLLOP_LORENTZ, COLLOP_FULL
    implicit none
    real(dp), parameter :: pi = 4*atan(1.0_dp), kev = 1.602176634e-16_dp
    real(dp), parameter :: charge = 1.602176634e-19_dp
    real(dp), parameter :: eps = 0.1_dp, q = 2.2_dp, rmajor = 6.2_dp, ig = 0.05_dp
    integer, parameter :: n_leg = 20
    character(64) :: check
    logical :: ok

    call get_command_argument(1, check)
    select case (trim(check))
    case ('axisymmetric'); ok = axisymmetric()
    case ('quasisymmetric'); ok = quasisymmetric()
    case ('monoenergetic'); ok = monoenergetic()
    case ('onsager'); ok = onsager()
    case ('ambipolar'); ok = ambipolar()
    case ('iterative'); ok = iterative()
    case ('exact'); ok = exact()
    case default
        print '(a)', 'unknown check '//trim(check)
        ok = .false.
    end select
    if (.not. ok) error stop 1

contains

    subroutine plasma(kappa_coll, sp, lnl, grid, nelem, order)
        real(dp), intent(in) :: kappa_coll
        type(collop_species_t), intent(out) :: sp(1)
        real(dp), intent(out) :: lnl(1, 1)
        type(collop_grid_t), intent(out) :: grid
        integer, intent(in), optional :: nelem, order
        integer :: st, ne, no
        ne = 3
        no = 4
        if (present(nelem)) ne = nelem
        if (present(order)) no = order
        sp(1) = collop_species_t(3.3435837768e-27_dp, 1.0_dp, 3e19_dp, 2*kev)
        lnl = 3*sqrt(pi)/8*collop_vth(sp(1))*kappa_coll/collop_nuhat(sp(1), sp(1), 1.0_dp)
        call collop_make_grid(5.0_dp, ne, no, grid, st)
        if (st /= COLLOP_OK) error stop 'grid'
    end subroutine plasma

    !> A tokamak with a helical 1/1 core ripple: not quasisymmetric.
    function helical(nt, nz) result(b)
        integer, intent(in) :: nt, nz
        real(dp) :: b(nt, nz), th, ze
        integer :: i, j
        do j = 1, nz
            do i = 1, nt
                th = 2*pi*(i - 1)/nt
                ze = 2*pi*(j - 1)/nz
                b(i, j) = 1 - eps*cos(th) - 0.04_dp*cos(th - ze) + 0.02_dp*cos(2*th - ze)
            end do
        end do
    end function helical

    logical function axisymmetric() result(ok)
        integer, parameter :: nt = 17, nz = 3
        type(collop_species_t) :: sp(1)
        type(collop_grid_t) :: grid
        type(ntv_surface_t) :: sf
        real(dp) :: lnl(1, 1), d2(1, 1, 3, 3), d3(1, 1, 3, 3), b(nt, nz)
        integer :: its, st, i
        call plasma(1e-3_dp, sp, lnl, grid)
        sf%eps = eps
        sf%q = q
        sf%k = 1/(q*rmajor)
        sf%n = 3
        sf%m = [6]
        sf%c = [(1e-3_dp, 0.0_dp)]
        sf%weight = NTV_WEIGHT_FULL
        sf%boozer_stream = .true.
        call ntv_axisymmetric_coefficients(sf, sp, lnl, grid, COLLOP_FULL, nt, n_leg, d2, &
                                              its, st)
        ok = st == NTV_OK
        b = spread([(1 - eps*cos(2*pi*(i - 1)/nt), i=1, nt)], 2, nz)
        call ntv_surface_coefficients(b, 1, q, ig, sf%k, sp, lnl, grid, COLLOP_FULL, n_leg, &
                                           d3, its, st)
        ok = ok .and. st == NTV_OK
        if (.not. ok) return
        ok = matrix_close('axisymmetric field on the 3D grid', d3(1, 1, :, :), &
                          d2(1, 1, :, :), 1e-10_dp)
    end function axisymmetric

    logical function quasisymmetric() result(ok)
        integer, parameter :: nt = 11, nz = 11, nper = 1
        real(dp), parameter :: iota = 1.25_dp, kappa = 0.3_dp
        type(collop_species_t) :: sp(1)
        type(collop_grid_t) :: grid
        type(ntv_surface_t) :: sf
        real(dp) :: lnl(1, 1), d2(1, 1, 3, 3), d3(1, 1, 3, 3), b(nt, nz), want(3, 3), &
                    scale(3)
        integer :: its, st, i, j
        call plasma(1e-3_dp, sp, lnl, grid)
        do j = 1, nz
            do i = 1, nt
                b(i, j) = 1 - eps*cos(2*pi*(i - 1)/nt - nper*2*pi*(j - 1)/(nz*nper))
            end do
        end do
        call ntv_surface_coefficients(b, nper, 1/iota, ig, kappa, sp, lnl, grid, &
                                      COLLOP_FULL, n_leg, d3, its, st)
        ok = st == NTV_OK
        sf%eps = eps
        sf%q = 1/iota
        sf%k = kappa*(1 - nper/iota)
        sf%n = 3
        sf%m = [6]
        sf%c = [(1e-3_dp, 0.0_dp)]
        sf%weight = NTV_WEIGHT_FULL
        sf%boozer_stream = .true.
        call ntv_axisymmetric_coefficients(sf, sp, lnl, grid, COLLOP_FULL, nt, n_leg, d2, &
                                              its, st)
        ok = ok .and. st == NTV_OK
        if (.not. ok) return
        scale = [1 + nper*ig, 1 + nper*ig, 1.0_dp]
        do j = 1, 3
            do i = 1, 3
                want(i, j) = scale(i)*scale(j)*d2(1, 1, i, j)
            end do
        end do
        ok = matrix_close('quasi-helical against the isomorphic tokamak', d3(1, 1, :, :), &
                             want, 1e-10_dp)
    end function quasisymmetric

    logical function monoenergetic() result(ok)
        integer, parameter :: nt = 11, nz = 7
        real(dp), parameter :: kappa = 1/(q*rmajor)
        type(collop_species_t) :: sp(1)
        type(collop_grid_t) :: grid
        type(mono_coefficients_t) :: c
        real(dp) :: lnl(1, 1), d(1, 1, 3, 3), b(nt, nz), want(3), x, w, nu, tz, vt, nuhat
        integer :: its, st, k
        call plasma(1e-2_dp, sp, lnl, grid)
        b = helical(nt, nz)
        call ntv_surface_coefficients(b, 1, q, ig, kappa, sp, lnl, grid, COLLOP_LORENTZ, &
                                        n_leg, d, its, st)
        ok = st == NTV_OK
        if (.not. ok) return
        tz = sp(1)%t/(sp(1)%z*charge)
        vt = collop_vth(sp(1))
        nuhat = collop_nuhat(sp(1), sp(1), lnl(1, 1))
        want = 0
        do k = 2, grid%n
            x = grid%x(k)
            w = grid%w(k)*exp(-x**2)
            nu = nuhat*(erf(x) - chandrasekhar(x))/x**3/(vt*x)
            call mono_coefficients_3d(b, 1, q, ig, kappa, nu, n_leg, c, st)
            if (st /= MONO_OK) then
                ok = .false.
                return
            end if
            want(1) = want(1) + (w*x**4)**2/(w*x**3)*c%d11
            want(2) = want(2) + w*x**4*c%d31
            want(3) = want(3) + w*x**3*c%d33
        end do
        want(1) = 8*tz**2/(sqrt(pi)*vt*kappa**2*q**2)*want(1)
        want(2) = -4*tz/(sqrt(pi)*kappa*q)*want(2)
        want(3) = 2*vt/sqrt(pi)*want(3)
        ok = close ('d11 against mono_coefficients_3d', d(1, 1, 1, 1), want(1), 1e-9_dp) &
            .and. close ('d31 against mono_coefficients_3d', d(1, 1, 3, 1), want(2), 1e-9_dp) &
            .and. close ('d33 against mono_coefficients_3d', d(1, 1, 3, 3), want(3), 1e-9_dp)
    end function monoenergetic

    logical function onsager() result(ok)
        integer, parameter :: nt = 15, nz = 9
        type(collop_species_t) :: sp(1)
        type(collop_grid_t) :: grid
        real(dp) :: lnl(1, 1), d(1, 1, 3, 3)
        integer :: its, st
        call plasma(1.0_dp, sp, lnl, grid)
        call ntv_surface_coefficients(helical(nt, nz), 1, q, ig, 1/(q*rmajor), sp, lnl, &
                                       grid, COLLOP_FULL, n_leg, d, its, st, pairs=0)
        ok = st == NTV_OK
        if (.not. ok) return
        ! GMRES stops at relative residual 1e-7; the discretisation defect at
        ! this collisionality is 4e-12 (pitch-angle scattering, direct).
        ok = close ('d13 = -d31', d(1, 1, 1, 3), -d(1, 1, 3, 1), 1e-7_dp) &
             .and. close ('d23 = -d32', d(1, 1, 2, 3), -d(1, 1, 3, 2), 1e-7_dp)
    end function onsager

    logical function ambipolar() result(ok)
        integer, parameter :: nt = 7, nz = 7
        type(collop_species_t) :: sp(1)
        type(collop_grid_t) :: grid
        real(dp) :: lnl(1, 1), d(1, 1, 3, 3), dl(1, 1, 3, 3), b(nt, nz)
        integer :: its, st, i, j
        call plasma(1e-1_dp, sp, lnl, grid, 6, 6)
        do j = 1, nz
            do i = 1, nt
                b(i, j) = 1 - eps*cos(2*pi*(i - 1)/nt - 2*pi*(j - 1)/nz)
            end do
        end do
        call ntv_surface_coefficients(b, 1, 0.8_dp, ig, 0.3_dp, sp, lnl, grid, COLLOP_FULL, &
                                           12, d, its, st, pairs=0)
        ok = st == NTV_OK
        call ntv_surface_coefficients(b, 1, 0.8_dp, ig, 0.3_dp, sp, lnl, grid, &
                                      COLLOP_LORENTZ, 12, dl, its, st, pairs=0)
        ok = ok .and. st == NTV_OK
        if (.not. ok) return
        print '(a,2es10.2)', 'd11/d11_lorentz, d12/d11_lorentz:', &
            d(1, 1, 1, 1)/dl(1, 1, 1, 1), &
            d(1, 1, 1, 2)/dl(1, 1, 1, 1)
        ok = abs(d(1, 1, 1, 1)) < 1e-5_dp*dl(1, 1, 1, 1) .and. &
             abs(d(1, 1, 1, 2)) < 1e-5_dp*dl(1, 1, 1, 1)
    end function ambipolar

    logical function iterative() result(ok)
        integer, parameter :: nt = 11, nz = 7
        type(collop_species_t) :: sp(1)
        type(collop_grid_t) :: grid
        real(dp) :: lnl(1, 1), d(1, 1, 3, 3), dg(1, 1, 3, 3)
        integer :: its, st
        call plasma(1e-2_dp, sp, lnl, grid)
        call ntv_surface_coefficients(helical(nt, nz), 1, q, ig, 1/(q*rmajor), sp, lnl, &
                                       grid, COLLOP_FULL, n_leg, d, its, st)
        ok = st == NTV_OK
        call ntv_surface_coefficients(helical(nt, nz), 1, q, ig, 1/(q*rmajor), sp, lnl, &
                                       grid, COLLOP_FULL, n_leg, dg, its, st, pairs=0)
        ok = ok .and. st == NTV_OK
        if (.not. ok) return
        print '(a,i0)', 'GMRES iterations: ', its
        ok = matrix_close('GMRES on collision-mode factorisations', dg(1, 1, :, :), &
                          d(1, 1, :, :), 1e-6_dp)
        call ntv_surface_coefficients(helical(nt, nz), 1, q, ig, 1/(q*rmajor), sp, lnl, &
                                      grid, COLLOP_FULL, n_leg, dg, its, st, pairs=1)
        ok = ok .and. st == NTV_OK
        if (.not. ok) return
        print '(a,i0)', 'GMRES iterations with the exact l <= 1 block: ', its
        ok = matrix_close('GMRES with the exact l <= 1 block', dg(1, 1, :, :), &
                          d(1, 1, :, :), 1e-6_dp) .and. ok
    end function iterative

    logical function exact() result(ok)
        integer, parameter :: nt = 11, nz = 7
        type(collop_species_t) :: sp(1)
        type(collop_grid_t) :: grid
        type(ntv_surface_t) :: sf
        real(dp) :: lnl(1, 1), d(1, 1, 3, 3), k(1, 1, 2)
        integer :: its, st, p
        call plasma(1e-2_dp, sp, lnl, grid)
        ok = .true.
        do p = 0, 1
            call ntv_surface_coefficients(helical(nt, nz), 1, q, ig, 1/(q*rmajor), sp, &
                                          lnl, grid, COLLOP_LORENTZ, n_leg, d, its, st, &
                                          pairs=p)
            print '(a,i0,a,i0)', 'surface, exact pair blocks ', p, ': iterations ', its
            ok = ok .and. st == NTV_OK .and. its == 3
        end do
        sf%eps = eps
        sf%q = q
        sf%k = 1/(q*rmajor)
        sf%n = 3
        sf%m = [6]
        sf%c = [(1e-3_dp, 0.0_dp)]
        sf%weight = NTV_WEIGHT_FULL
        sf%boozer_stream = .true.
        call ntv_energy_coefficients(sf, sp, lnl, grid, COLLOP_LORENTZ, 2e3_dp, 17, n_leg, &
                                     k, its, st)
        print '(a,i0)', 'NTV kernel, per-node chain: iterations ', its
        ok = ok .and. st == NTV_OK .and. its == 2
    end function exact

    !> G(x) = (erf(x) - 2 x exp(-x^2)/sqrt(pi))/(2 x^2).
    real(dp) function chandrasekhar(x)
        real(dp), intent(in) :: x
        chandrasekhar = (erf(x) - 2*x*exp(-x*x)/sqrt(pi))/(2*x*x)
    end function chandrasekhar

    !> Entrywise against the largest entry of each row and column pair.
    logical function matrix_close(label, got, want, rtol) result(ok)
        character(*), intent(in) :: label
        real(dp), intent(in) :: got(3, 3), want(3, 3), rtol
        real(dp) :: err, s
        integer :: i, j
        err = 0
        do j = 1, 3
            do i = 1, 3
                s = sqrt(abs(want(i, i)*want(j, j)))
                err = max(err, abs(got(i, j) - want(i, j))/s)
            end do
        end do
        ok = err <= rtol
        print '(a,es10.2,a,l2)', label//': largest error ', err, ' pass', ok
    end function matrix_close

    logical function close (label, got, want, rtol)
        character(*), intent(in) :: label
        real(dp), intent(in) :: got, want, rtol
        close = abs(got - want) <= rtol*abs(want)
        print '(a,2es20.11,a,es9.2,a,l2)', label//' ', got, want, ' rel ', &
            abs(got - want)/abs(want), ' pass', close
    end function close
end program test_ntv_surface
