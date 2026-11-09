!> Finite-width guiding-centre orbits on a Boozer file (milca_boozer_orbit),
!! on a synthetic circular surface family, B = b0 (1 - 0.1 sqrt(s) cos(theta)):
!!   deeply_trapped  a deeply trapped thin orbit precesses at dH/dP_zeta for
!!                   H = mu B(s*, 0) + Z e Phi(s*), i.e. at
!!                   omega_e - mu dB/ds(s*, 0)/(Z e iota psi'), with and
!!                   without E x B rotation;
!!   thin_bounce     the bounce frequency tends to the thin-orbit quadrature
!!                   2 pi/oint G dtheta/(v_par iota B) at s* as the mass goes
!!                   to zero at fixed speed; the correction is quadratic in the
!!                   mass (orbit width) because the field is up-down symmetric;
!!   invariants      every orbit point carries the energy of the orbit;
!!   shear           the precession of a trapped orbit with magnetic shear is the
!!                   thin bounce average of the binormal drift with the
!!                   current and shear term of NEO-2-QL,
!!                   (m v^2/(2 e)) [(1 + xi^2) d + xi^2 d2], d = -(1/iota)
!!                   (dB/ds)/(psi' B), d2 = 2 (I' + (q G)')/(psi' (G + iota I));
!!                   without d2 it misses the shear contribution;
!!   pitch_map       the pitch map returns the magnetic moment of the orbit at
!!                   every point of it (mu is the stream function of the
!!                   orbits at fixed s* and H), and mu B(s*, theta)/E tends to
!!                   1 - xi^2 for thin orbits.
program test_boozer_orbit
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use milca_boozer_file, only: boozer_file_t, load_boozer_file, boozer_surface_t, &
        boozer_surface, BOOZER_OK
    use milca_boozer_orbit, only: orbit_particle_t, orbit_point_t, orbit_point, &
        orbit_frequencies, orbit_state, orbit_pitch_point, ORBIT_OK
    use fortnum_quadrature, only: gauss_legendre
    implicit none
    real(dp), parameter :: pi = 4*atan(1.0_dp), e = 1.602176634e-19_dp
    real(dp), parameter :: md = 3.3435837768e-27_dp, kev = 1e3_dp*e
    real(dp), parameter :: s_star = 0.45_dp
    character(*), parameter :: file = 'boozer_orbit.bc'
    type(boozer_file_t) :: bf
    character(64) :: check
    logical :: ok
    integer :: status

    call write_file(file)
    call load_boozer_file(file, bf, status)
    if (status /= BOOZER_OK) error stop 'cannot read the test file'
    call get_command_argument(1, check)
    select case (trim(check))
    case ('deeply_trapped'); ok = deeply_trapped()
    case ('thin_bounce'); ok = thin_bounce()
    case ('invariants'); ok = invariants()
    case ('pitch_map'); ok = pitch_map()
    case ('shear'); ok = shear()
    case default
        print '(a)', 'unknown check '//trim(check)
        ok = .false.
    end select
    if (.not. ok) error stop 1

contains

    !> Deuteron-like particle of the given mass factor at 10 keV kinetic energy
    !> with kinetic pitch lambda = mu B(s*, 0)/E_kin.
    function particle(mass_factor, lambda, omega_e) result(p)
        real(dp), intent(in) :: mass_factor, lambda, omega_e
        type(orbit_particle_t) :: p
        real(dp) :: b, b_s, g, i, g_s, i_s, iota, psi_prime
        integer :: st
        call orbit_state(bf, s_star, 0.0_dp, b, b_s, g, i, g_s, i_s, iota, psi_prime, st)
        p%mass = md*mass_factor
        p%charge = e
        p%energy = 10*kev*mass_factor
        p%mu = lambda*p%energy/b
        p%omega_e = omega_e
        p%s_ref = s_star
    end function particle

    logical function deeply_trapped() result(ok)
        real(dp), parameter :: factor = 1e-4_dp
        real(dp) :: ob, oz, want, b, b_s, g, i, g_s, i_s, iota, psi_prime, omega_e
        type(orbit_particle_t) :: p
        logical :: trapped
        integer :: st, k
        ok = .true.
        call orbit_state(bf, s_star, 0.0_dp, b, b_s, g, i, g_s, i_s, iota, psi_prime, st)
        do k = 1, 2
            omega_e = merge(0.0_dp, 2e3_dp, k == 1)
            ! lambda just below 1: bounce amplitude of a few degrees.
            p = particle(factor, 0.9999_dp, omega_e)
            call orbit_frequencies(bf, p, s_star, 1.0_dp, 64, ob, oz, trapped, st)
            want = omega_e - p%mu*b_s/(p%charge*iota*psi_prime)
            ok = ok .and. st == ORBIT_OK .and. trapped .and. close('precession', oz, want, 2e-3_dp)
        end do
    end function deeply_trapped

    logical function thin_bounce() result(ok)
        real(dp) :: ob(3), oz, thin, f(3)
        type(orbit_particle_t) :: p
        logical :: trapped
        integer :: st, k
        f = [1e-2_dp, 5e-3_dp, 2.5e-3_dp]
        do k = 1, 3
            p = particle(f(k), 0.97_dp, 0.0_dp)
            call orbit_frequencies(bf, p, s_star, 1.0_dp, 64, ob(k), oz, trapped, st)
            if (st /= ORBIT_OK .or. .not. trapped) then
                ok = .false.
                return
            end if
        end do
        thin = thin_bounce_frequency(particle(1.0_dp, 0.97_dp, 0.0_dp))
        ! Quadratic in the mass: in the up-down symmetric field the first-order
        ! shifts of the two legs cancel, so the differences quarter.
        ok = close('bounce frequency, mass -> 0', (4*ob(3) - ob(2))/3, thin, 1e-6_dp) &
            .and. close('second-order convergence', (ob(1) - ob(2))/(ob(2) - ob(3)), 4.0_dp, &
            5e-2_dp)
    end function thin_bounce

    !> 2 pi/oint G dtheta/(v_par iota B) at s*, both legs, by the substitution
    !> theta = Theta sin(u) between the turning points +-Theta (B even in theta).
    real(dp) function thin_bounce_frequency(p) result(omega)
        type(orbit_particle_t), intent(in) :: p
        real(dp) :: x(64), w(64), b, b_s, g, i, g_s, i_s, iota, psi_prime, tmax, lo, hi, &
            th, v_par, tau
        integer :: k, st
        call gauss_legendre(64, x, w)
        lo = 0
        hi = pi
        do k = 1, 200
            tmax = 0.5_dp*(lo + hi)
            call orbit_state(bf, s_star, tmax, b, b_s, g, i, g_s, i_s, iota, psi_prime, st)
            if (p%mu*b > p%energy) then
                hi = tmax
            else
                lo = tmax
            end if
        end do
        tau = 0
        do k = 1, 64
            th = tmax*sin(0.5_dp*pi*x(k))
            call orbit_state(bf, s_star, th, b, b_s, g, i, g_s, i_s, iota, psi_prime, st)
            v_par = sqrt(2*(p%energy - p%mu*b)/p%mass)
            tau = tau + 2*w(k)*tmax*cos(0.5_dp*pi*x(k))*0.5_dp*pi*abs((g + iota*i)/(iota*b))/v_par
        end do
        omega = 2*pi/tau
    end function thin_bounce_frequency

    logical function invariants() result(ok)
        type(orbit_particle_t) :: p
        type(orbit_point_t) :: pt
        real(dp) :: th, worst
        integer :: k, st
        p = particle(1.0_dp, 0.97_dp, 3e3_dp)
        worst = 0
        do k = -4, 4
            th = 0.1_dp*k
            call orbit_point(bf, p, s_star, th, merge(1.0_dp, -1.0_dp, k >= 0), pt, st)
            if (st /= ORBIT_OK) then
                ok = .false.
                return
            end if
            worst = max(worst, abs(p%mass*pt%v_par**2/2 + p%mu*pt%b &
                - p%charge*p%omega_e*psi_p(pt%s) - p%energy)/p%energy)
        end do
        ok = worst < 1e-10_dp
        print '(a,es10.2)', 'energy error', worst
    end function invariants

    logical function shear() result(ok)
        real(dp), parameter :: lams(3) = [0.999_dp, 0.98_dp, 0.95_dp]
        type(orbit_particle_t) :: p
        real(dp) :: x(96), w(96), ob, oz, b, b_s, g, i, g_s, i_s, iota, psi_prime, iota_s, &
            tmax, lo, hi, th, v_par, dt, tau, num, d, d2, xi2, jac
        type(boozer_surface_t) :: surf
        logical :: trapped
        integer :: st, k, j, it
        call gauss_legendre(96, x, w)
        call boozer_surface(bf, s_star, surf, st)
        iota_s = surf%diota_ds
        ok = .true.
        do j = 1, 3
            p = particle(1e-3_dp, 1.0_dp, 0.0_dp)
            p%mu = lams(j)*p%energy/2.3_dp
            call orbit_frequencies(bf, p, s_star, 1.0_dp, 96, ob, oz, trapped, st)
            lo = 0
            hi = pi
            do it = 1, 200
                tmax = 0.5_dp*(lo + hi)
                call orbit_state(bf, s_star, tmax, b, b_s, g, i, g_s, i_s, iota, psi_prime, st)
                if (p%mu*b > p%energy) then
                    hi = tmax
                else
                    lo = tmax
                end if
            end do
            tau = 0
            num = 0
            do k = 1, 96
                th = tmax*sin(0.5_dp*pi*x(k))
                jac = tmax*cos(0.5_dp*pi*x(k))*0.5_dp*pi
                call orbit_state(bf, s_star, th, b, b_s, g, i, g_s, i_s, iota, psi_prime, st)
                v_par = sqrt(2*(p%energy - p%mu*b)/p%mass)
                xi2 = p%mass*v_par**2/(2*p%energy)
                dt = w(k)*jac*abs((g + iota*i)/(iota*b))/v_par
                d = -(b_s/(psi_prime*b))/iota
                d2 = 2*(i_s + g_s/iota - g*iota_s/(iota*iota))/(psi_prime*(g + iota*i))
                tau = tau + dt
                num = num + dt*(p%energy/p%charge)*((1 + xi2)*d + xi2*d2)
            end do
            ok = ok .and. st == ORBIT_OK .and. trapped .and. &
                close('precession with shear', oz, num/tau, 1e-4_dp)
        end do
    end function shear

    logical function pitch_map() result(ok)
        type(orbit_particle_t) :: p, thin
        type(orbit_point_t) :: pt, back
        real(dp) :: th, worst, xi, v, b, b_s, g, i, g_s, i_s, iota, psi_prime, e_kin
        integer :: k, st
        p = particle(1.0_dp, 0.97_dp, 3e3_dp)
        worst = 0
        do k = -4, 4
            th = 0.1_dp*k
            call orbit_point(bf, p, s_star, th, merge(1.0_dp, -1.0_dp, k >= 0), pt, st)
            if (st /= ORBIT_OK) exit
            ! Kinetic energy E - Z e Phi with Phi = -omega_e (psi_p - psi_p(s*)).
            v = sqrt(2*(p%energy + p%charge*p%omega_e*psi_p(pt%s))/p%mass)
            xi = pt%v_par/v
            call orbit_pitch_point(bf, p, s_star, th, xi, back, st)
            if (st /= ORBIT_OK) exit
            worst = max(worst, abs(back%mu/p%mu - 1), abs(back%s - pt%s))
        end do
        ok = st == ORBIT_OK .and. worst < 1e-10_dp
        print '(a,es10.2)', 'mu and s along the orbit, worst deviation', worst
        thin = particle(1e-6_dp, 0.97_dp, 0.0_dp)
        xi = 0.3_dp
        call orbit_pitch_point(bf, thin, s_star, 1.0_dp, xi, back, st)
        call orbit_state(bf, s_star, 1.0_dp, b, b_s, g, i, g_s, i_s, iota, psi_prime, st)
        e_kin = thin%energy
        ok = ok .and. close('thin mu B/E', back%mu*b/e_kin, 1 - xi**2, 1e-5_dp)
    end function pitch_map

    !> psi_p(s) - psi_p(s*) = psi' int iota ds, exact by the trapezoid for the
    !> linear iota of the file.
    real(dp) function psi_p(s)
        real(dp), intent(in) :: s
        type(boozer_surface_t) :: a, b
        integer :: st
        call boozer_surface(bf, s_star, a, st)
        call boozer_surface(bf, s, b, st)
        psi_p = a%psi_prime*0.5_dp*(a%iota + b%iota)*(s - s_star)
    end function psi_p

    logical function close(label, got, want, rtol)
        character(*), intent(in) :: label
        real(dp), intent(in) :: got, want, rtol
        close = abs(got - want) <= rtol*abs(want)
        print '(a,2es20.11,a,es9.2,a,l2)', label//' ', got, want, ' rel ', &
            abs(got - want)/abs(want), ' pass', close
    end function close

    !> Concentric circles r = 0.5 sqrt(s) around R = 1.7 m, B = 2.3 (1 - 0.1
    !> sqrt(s) cos(theta)) T, iota = 0.9 - 0.3 s, constant poloidal current.
    subroutine write_file(name)
        character(*), intent(in) :: name
        integer :: unit, k
        real(dp) :: s
        open (newunit=unit, file=name, status='replace', action='write')
        write (unit, '(a)') 'CC synthetic test file'
        write (unit, '(a)') ' m0b n0b nsurf nper flux/[Tm^2] a/[m] R/[m]'
        write (unit, '(a)') ' 1 0 19 1 -1.3 0.5 1.7'
        do k = 1, 19
            s = 0.05_dp*k
            write (unit, '(a)') ' s iota curr_pol/nper curr_tor pprime sqrt g(0,0)'
            write (unit, '(a)') ' [A] [A] dp/ds,[Pa] (dV/ds)/nper'
            write (unit, '(6es24.16)') s, 0.9_dp - 0.3_dp*s, -1.2e7_dp, 0.0_dp, 0.0_dp, 0.0_dp
            write (unit, '(a)') ' m n r/[m] z/[m] (phib-phi)*nper/twopi bmn/[T]'
            write (unit, '(2i5,4es24.16)') 0, 0, 1.7_dp, 0.0_dp, 0.0_dp, 2.3_dp
            write (unit, '(2i5,4es24.16)') 1, 0, 0.5_dp*sqrt(s), 0.5_dp*sqrt(s), 0.0_dp, &
                -0.23_dp*sqrt(s)
        end do
        close (unit)
    end subroutine write_file
end program test_boozer_orbit
