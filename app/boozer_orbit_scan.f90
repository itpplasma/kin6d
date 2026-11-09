!> Poloidal and toroidal frequencies of finite-width guiding-centre orbits
!! started on the outboard midplane (Z = 0) of a Boozer file, the scan of
!! POTATO's itest_type = 5 (milca_boozer_orbit).
!! Usage: boozer_orbit_scan FILE A Z E_keV PITCH R_MIN R_MAX N [OMEGA_E]
!! with the mass number A, charge number Z, kinetic energy, v_par/v at the
!! start point, the start radii [m] and the E x B rotation [rad/s].
!! Prints per start radius: R [m], s, s*, trapped, omega_theta and omega_zeta
!! [rad/s] (the poloidal transit or bounce frequency and the toroidal advance
!! per poloidal period over that period).
program boozer_orbit_scan
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use milca_boozer_file, only: boozer_file_t, load_boozer_file, boozer_surface_t, &
        boozer_surface, boozer_point, boozer_point_t, BOOZER_OK
    use milca_boozer_orbit, only: orbit_particle_t, orbit_frequencies, orbit_state, &
        orbit_flux_p, ORBIT_OK
    implicit none
    real(dp), parameter :: e = 1.602176634e-19_dp, amu = 1.66053906660e-27_dp
    type(boozer_file_t) :: bf
    type(orbit_particle_t) :: p
    character(1024) :: arg
    real(dp) :: a_mass, z, e_kev, pitch, r_min, r_max, omega_e, r, s, s_star, v, b, b_s, g, &
        i_cov, g_s, i_s, iota, psi_prime, rho, ob, oz
    integer :: n, k, st
    logical :: trapped

    if (command_argument_count() < 8) &
        error stop 'usage: boozer_orbit_scan FILE A Z E_keV PITCH R_MIN R_MAX N [OMEGA_E]'
    call get_command_argument(1, arg)
    call load_boozer_file(trim(arg), bf, st)
    if (st /= BOOZER_OK) error stop 'cannot read the Boozer file'
    a_mass = real_arg(2)
    z = real_arg(3)
    e_kev = real_arg(4)
    pitch = real_arg(5)
    r_min = real_arg(6)
    r_max = real_arg(7)
    n = nint(real_arg(8))
    omega_e = 0
    if (command_argument_count() >= 9) omega_e = real_arg(9)
    p%mass = a_mass*amu
    p%charge = z*e
    p%omega_e = omega_e
    v = sqrt(2*e_kev*1e3_dp*e/p%mass)
    do k = 0, n - 1
        r = r_min + (r_max - r_min)*k/max(n - 1, 1)
        s = surface_at(r, st)
        if (st /= 0) cycle
        call orbit_state(bf, s, midplane_theta(s), b, b_s, g, i_cov, g_s, i_s, iota, &
            psi_prime, st)
        p%mu = p%mass*v**2*(1 - pitch**2)/(2*b)
        p%s_ref = s
        p%energy = p%mass*v**2/2
        rho = p%mass*pitch*v/(p%charge*b)
        ! P_zeta = rho G - psi_p(s) = -psi_p(s*).
        s_star = orbit_surface(s, rho*g, st)
        if (st /= 0) cycle
        call orbit_frequencies(bf, p, s_star, sign(1.0_dp, pitch), 96, ob, oz, trapped, st)
        if (st /= ORBIT_OK) then
            write (error_unit, '(a,f8.4,a,i0)') 'R = ', r, ': orbit status ', st
            cycle
        end if
        print '(f10.5,2f12.8,l3,2es16.7)', r, s, s_star, trapped, ob, oz
    end do

contains

    real(dp) function real_arg(i)
        integer, intent(in) :: i
        character(64) :: buf
        call get_command_argument(i, buf)
        read (buf, *) real_arg
    end function real_arg

    !> s with R = r on the outboard midplane by bisection (R grows outwards).
    real(dp) function surface_at(r, st)
        real(dp), intent(in) :: r
        integer, intent(out) :: st
        real(dp) :: lo, hi
        integer :: it
        lo = bf%d%s(1)
        hi = bf%d%s(size(bf%d%s))
        st = 1
        if ((radius(lo) - r)*(radius(hi) - r) > 0) return
        do it = 1, 100
            surface_at = 0.5_dp*(lo + hi)
            if ((radius(surface_at) - r)*(radius(lo) - r) > 0) then
                lo = surface_at
            else
                hi = surface_at
            end if
        end do
        st = 0
    end function surface_at

    real(dp) function radius(s)
        real(dp), intent(in) :: s
        type(boozer_surface_t) :: surf
        type(boozer_point_t) :: pt
        integer :: st
        call boozer_surface(bf, s, surf, st)
        pt = boozer_point(surf, midplane_theta(s), 0.0_dp)
        radius = pt%r
    end function radius

    !> Boozer theta of the outboard point with Z = 0 on surface s, by
    !> bisection in [-pi/2, pi/2], where Z increases with theta.
    real(dp) function midplane_theta(s) result(theta)
        real(dp), intent(in) :: s
        type(boozer_surface_t) :: surf
        type(boozer_point_t) :: pt
        real(dp) :: lo, hi
        integer :: it, st
        call boozer_surface(bf, s, surf, st)
        lo = -2*atan(1.0_dp)
        hi = -lo
        do it = 1, 80
            theta = 0.5_dp*(lo + hi)
            pt = boozer_point(surf, theta, 0.0_dp)
            if (pt%z > 0) then
                hi = theta
            else
                lo = theta
            end if
        end do
    end function midplane_theta

    !> s* with psi_p(s*) = psi_p(s) - rho G.
    real(dp) function orbit_surface(s, shift, st)
        real(dp), intent(in) :: s, shift
        integer, intent(out) :: st
        real(dp) :: lo, hi, f_lo, f
        integer :: it
        lo = bf%d%s(1)
        hi = bf%d%s(size(bf%d%s))
        st = 1
        f_lo = orbit_flux_p(bf, s, lo) + shift
        if (f_lo*(orbit_flux_p(bf, s, hi) + shift) > 0) return
        do it = 1, 100
            orbit_surface = 0.5_dp*(lo + hi)
            f = orbit_flux_p(bf, s, orbit_surface) + shift
            if (f*f_lo > 0) then
                lo = orbit_surface
                f_lo = f
            else
                hi = orbit_surface
            end if
        end do
        st = 0
    end function orbit_surface
end program boozer_orbit_scan
