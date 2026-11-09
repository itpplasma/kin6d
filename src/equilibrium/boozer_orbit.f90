!> Unperturbed guiding-centre orbits at finite width on the axisymmetric
!! equilibrium of a Boozer file, from the canonical Lagrangian of White and
!! Chance (Phys. Fluids 27, 2455 (1984)) per unit charge,
!!   L = P_theta theta_dot + P_zeta zeta_dot - H,
!!   P_theta = psi + rho I,  P_zeta = rho G - psi_p,
!!   H = (Z e/m) rho^2 B^2/2 + (mu/(Z e)) B + Phi(psi_p),
!! with psi = psi_tor/(2 pi) = psi' s, psi_p the poloidal flux per radian,
!! dpsi_p/dpsi = iota, rho = m v_par/(Z e B), G and I the covariant Boozer
!! components [T m] and H in volts. For fixed (P_zeta, H, mu) the orbit is a
!! curve psi(theta); its turning points (rho = 0) lie on psi_p = -P_zeta,
!! the surface s* of psi* = p_phi/(Z e). With D = G + iota I + rho (G I' - I G'),
!! ' = d/dpsi, Hamilton's equations give
!!   theta_dot = (H_psi G - H_rho (rho G' - iota))/D,
!!   zeta_dot = (H_rho (1 + rho I') - H_psi I)/D,
!! exact at any orbit width within the model of the file (axisymmetric,
!! no radial covariant component). The electrostatic potential is linear,
!! dPhi/dpsi_p = -omega_e, so that omega_e is the E x B rotation frequency,
!! unless the particle carries a potential profile Phi(s) [V] (orbit_set_potential),
!! which then enters the energy and the radial derivative of H exactly
!! (E_r shear and orbit squeezing).
!! Lengths in m, fluxes in T m^2, energies in J (mu in J/T), times in s.
module milca_boozer_orbit
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_quadrature, only: gauss_legendre
    use fortnum_cubic_spline, only: cubic_spline_t, cubic_spline_init, cubic_spline_eval
    use milca_boozer_file, only: boozer_file_t, boozer_surface_t, boozer_surface, BOOZER_OK
    use milca_boozer_orbit_leaves, only: orbit_pitch_residual, orbit_pitch_residual_jvp, &
        orbit_moment, orbit_moment_jvp, orbit_energy_residual, orbit_energy_residual_jvp
    implicit none
    private
    public :: orbit_particle_t, orbit_point_t, orbit_state, orbit_point, orbit_frequencies, &
        orbit_flux_p, orbit_pitch_point, orbit_set_potential, orbit_potential
    integer, parameter, public :: ORBIT_OK = 0, ORBIT_BAD_INPUT = 1, ORBIT_NO_ROOT = 2

    real(dp), parameter :: pi = 4*atan(1.0_dp)

    !> One particle class: mass [kg], charge Z e [C], magnetic moment [J/T],
    !> total energy [J] including Z e Phi, dPhi/dpsi_p = -omega_e [rad/s], and
    !> the reference surface s_ref where Phi = 0. The orbit itself is fixed by
    !> its surface s* (P_zeta = -psi_p(s*)).
    type :: orbit_particle_t
        real(dp) :: mass = 0, charge = 0, mu = 0, energy = 0
        real(dp) :: omega_e = 0, s_ref = 0
        logical :: has_potential = .false.
        type(cubic_spline_t) :: potential
    end type

    !> The orbit at one theta: s, rho [m], v_par [m/s], B [T], theta_dot and
    !> zeta_dot [rad/s], and the magnetic moment mu [J/T] of the orbit through
    !> the point (set by orbit_pitch_point).
    type :: orbit_point_t
        real(dp) :: s = 0, rho = 0, v_par = 0, b = 0, theta_dot = 0, zeta_dot = 0, mu = 0
        ! dmu/dxi along the ray and ds/ds* at fixed (theta, xi, H).
        real(dp) :: mu_xi = 0, s_star_s = 0
    end type

contains

    !> Potential profile Phi(s) [V] on increasing s for particle p (cubic
    !> spline); Phi then replaces the linear potential of omega_e.
    subroutine orbit_set_potential(p, s, phi, status)
        type(orbit_particle_t), intent(inout) :: p
        real(dp), intent(in) :: s(:), phi(:)
        integer, intent(out) :: status
        call cubic_spline_init(s, reshape(phi, [size(phi), 1]), p%potential, status)
        p%has_potential = status == 0
    end subroutine orbit_set_potential

    !> Phi(s) - Phi(s_ref) [V] and dPhi/ds of the particle's potential.
    subroutine orbit_potential(bf, p, s, phi, phi_s)
        type(boozer_file_t), intent(in) :: bf
        type(orbit_particle_t), intent(in) :: p
        real(dp), intent(in) :: s
        real(dp), intent(out) :: phi, phi_s
        real(dp) :: v(1), dv(1), ref(1), ref_s(1)
        type(boozer_surface_t) :: surf
        integer :: st
        if (p%has_potential) then
            call cubic_spline_eval(p%potential, s, v, dv)
            call cubic_spline_eval(p%potential, p%s_ref, ref, ref_s)
            phi = v(1) - ref(1)
            phi_s = dv(1)
        else
            call boozer_surface(bf, s, surf, st)
            phi = -p%omega_e*orbit_flux_p(bf, p%s_ref, s)
            phi_s = -p%omega_e*surf%iota*surf%psi_prime
        end if
    end subroutine orbit_potential

    !> psi_p(s) - psi_p(s0) = psi' int_s0^s iota ds by 8-point Gauss-Legendre.
    function orbit_flux_p(bf, s0, s) result(dpsi_p)
        type(boozer_file_t), intent(in) :: bf
        real(dp), intent(in) :: s0, s
        real(dp) :: dpsi_p, x(8), w(8)
        type(boozer_surface_t) :: surf
        integer :: k, st
        call gauss_legendre(8, x, w)
        dpsi_p = 0
        do k = 1, 8
            call boozer_surface(bf, 0.5_dp*(s0 + s) + 0.5_dp*(s - s0)*x(k), surf, st)
            dpsi_p = dpsi_p + w(k)*surf%iota
        end do
        dpsi_p = 0.5_dp*(s - s0)*dpsi_p*surf%psi_prime
    end function orbit_flux_p

    !> Field quantities at (s, theta): B, dB/ds, G, I, dG/ds, dI/ds, iota, psi'.
    subroutine orbit_state(bf, s, theta, b, b_s, g, i, g_s, i_s, iota, psi_prime, status, &
            b_t, b_tt)
        type(boozer_file_t), intent(in) :: bf
        real(dp), intent(in) :: s, theta
        real(dp), intent(out) :: b, b_s, g, i, g_s, i_s, iota, psi_prime
        integer, intent(out) :: status
        ! dB/dtheta and d2B/dtheta2 of the Fourier series.
        real(dp), intent(out), optional :: b_t, b_tt
        type(boozer_surface_t) :: surf
        integer :: k
        call boozer_surface(bf, s, surf, status)
        if (status /= BOOZER_OK) then
            status = ORBIT_BAD_INPUT
            return
        end if
        b = 0
        b_s = 0
        if (present(b_t)) b_t = 0
        if (present(b_tt)) b_tt = 0
        do k = 1, size(surf%m)
            if (surf%n(k) /= 0) cycle
            b = b + surf%bmn(k)*cos(surf%m(k)*theta) + surf%bma(k)*sin(surf%m(k)*theta)
            b_s = b_s + surf%bmn_s(k)*cos(surf%m(k)*theta) + surf%bma_s(k)*sin(surf%m(k)*theta)
            if (present(b_t)) b_t = b_t + surf%m(k)*(surf%bma(k)*cos(surf%m(k)*theta) &
                - surf%bmn(k)*sin(surf%m(k)*theta))
            if (present(b_tt)) b_tt = b_tt - surf%m(k)**2*(surf%bmn(k)*cos(surf%m(k)*theta) &
                + surf%bma(k)*sin(surf%m(k)*theta))
        end do
        g = surf%g_cov
        i = surf%i_cov
        g_s = surf%g_cov_s
        i_s = surf%i_cov_s
        iota = surf%iota
        psi_prime = surf%psi_prime
        status = ORBIT_OK
    end subroutine orbit_state

    !> The orbit point at theta on the leg sigma (sign of v_par). rho follows
    !> from P_zeta, and s is the root of H(s) on the side of s* that sigma
    !> selects, bracketed outwards from s* and refined by Newton's method on
    !> the generated energy residual.
    subroutine orbit_point(bf, p, s_star, theta, sigma, pt, status)
        type(boozer_file_t), intent(in) :: bf
        type(orbit_particle_t), intent(in) :: p
        real(dp), intent(in) :: s_star, theta, sigma
        type(orbit_point_t), intent(out) :: pt
        integer, intent(out) :: status
        real(dp) :: lo, hi, mid, f_lo, f_hi, f_mid, step, dir, df
        integer :: it
        status = ORBIT_BAD_INPUT
        if (p%mass <= 0 .or. p%charge == 0) return
        dir = sigma*sign(1.0_dp, leg_sign())
        lo = s_star
        f_lo = residual(lo)
        ! At a turning point f_lo vanishes; rounding must not make it a miss.
        if (f_lo > 0 .and. f_lo <= 1e-12_dp*abs(p%energy)) then
            call evaluate(s_star, pt, status)
            return
        end if
        if (f_lo > 0) then
            ! The orbit does not reach this theta.
            status = ORBIT_NO_ROOT
            return
        end if
        step = 1e-3_dp
        do it = 1, 60
            hi = s_star + dir*step
            if (hi <= 0 .or. hi >= 1) then
                status = ORBIT_NO_ROOT
                return
            end if
            f_hi = residual(hi)
            if (f_hi > 0) exit
            lo = hi
            f_lo = f_hi
            step = 2*step
        end do
        if (f_hi <= 0) then
            status = ORBIT_NO_ROOT
            return
        end if
        ! Newton on the generated energy residual, bisecting a step that
        ! leaves the bracket [lo, hi] (residual <= 0 at lo, > 0 at hi).
        mid = 0.5_dp*(lo + hi)
        do it = 1, 100
            call residual_ds(mid, f_mid, df)
            if (f_mid == 0) exit
            if (f_mid > 0) then
                hi = mid
            else
                lo = mid
            end if
            step = f_mid/df
            if (df == 0 .or. (mid - step - lo)*(mid - step - hi) > 0) step = mid - 0.5_dp*(lo + hi)
            mid = mid - step
            if (abs(step) <= 4*epsilon(1.0_dp)*max(1.0_dp, abs(mid))) exit
        end do
        call evaluate(mid, pt, status)
    contains
        !> sign of rho on the side s > s*: rho G = P_zeta + psi_p(s) grows like
        !> iota psi' (s - s*).
        real(dp) function leg_sign()
            real(dp) :: b, b_s, g, i, g_s, i_s, iota, psi_prime
            integer :: st
            call orbit_state(bf, s_star, theta, b, b_s, g, i, g_s, i_s, iota, psi_prime, st)
            leg_sign = iota*psi_prime*g
        end function leg_sign

        real(dp) function residual(s) result(r)
            real(dp), intent(in) :: s
            real(dp) :: r_s
            call residual_ds(s, r, r_s)
        end function residual

        !> The energy residual and its s derivative from the generated leaf.
        subroutine residual_ds(s, r, r_s)
            real(dp), intent(in) :: s
            real(dp), intent(out) :: r, r_s
            real(dp) :: bb, bs, gg, ii, gs, is, io, pp, ph, ph_s, dpsi_p
            integer :: st
            call orbit_state(bf, s, theta, bb, bs, gg, ii, gs, is, io, pp, st)
            ! Outside the file counts as beyond the orbit.
            r = huge(1.0_dp)
            r_s = 0
            if (st /= ORBIT_OK) return
            call orbit_potential(bf, p, s, ph, ph_s)
            dpsi_p = orbit_flux_p(bf, s_star, s)
            call orbit_energy_residual(dpsi_p, gg, bb, p%mass, p%charge, p%mu, ph, p%energy, r)
            call orbit_energy_residual_jvp(dpsi_p, gg, bb, p%mass, p%charge, p%mu, ph, &
                p%energy, io*pp, gs, bs, 0.0_dp, 0.0_dp, 0.0_dp, ph_s, 0.0_dp, r_s)
        end subroutine residual_ds

        subroutine evaluate(s, q, st)
            real(dp), intent(in) :: s
            type(orbit_point_t), intent(out) :: q
            integer, intent(out) :: st
            call orbit_evaluate(bf, p, s_star, s, theta, q, st)
        end subroutine evaluate
    end subroutine orbit_point

    !> The orbit of surface s_star at (s, theta): rho from P_zeta, v_par, B and
    !> Hamilton's theta_dot, zeta_dot (module header) with the moment of p.
    subroutine orbit_evaluate(bf, p, s_star, s, theta, q, st)
        type(boozer_file_t), intent(in) :: bf
        type(orbit_particle_t), intent(in) :: p
        real(dp), intent(in) :: s_star, s, theta
        type(orbit_point_t), intent(out) :: q
        integer, intent(out) :: st
        real(dp) :: b, b_s, g, i, g_s, i_s, iota, psi_prime, h_psi, h_rho, d, rho, phi, phi_s
        call orbit_state(bf, s, theta, b, b_s, g, i, g_s, i_s, iota, psi_prime, st)
        if (st /= ORBIT_OK) return
        rho = orbit_flux_p(bf, s_star, s)/g
        q%s = s
        q%rho = rho
        q%b = b
        q%v_par = p%charge*rho*b/p%mass
        q%mu = p%mu
        ! Derivatives per unit psi (toroidal flux per radian) and per volt.
        call orbit_potential(bf, p, s, phi, phi_s)
        h_psi = (p%charge*rho**2*b/p%mass + p%mu/p%charge)*b_s/psi_prime + phi_s/psi_prime
        h_rho = p%charge*rho*b**2/p%mass
        d = g + iota*i + rho*(g*i_s - i*g_s)/psi_prime
        q%theta_dot = (h_psi*g - h_rho*(rho*g_s/psi_prime - iota))/d
        q%zeta_dot = (h_rho*(1 + rho*i_s/psi_prime) - h_psi*i)/d
    end subroutine orbit_evaluate

    !> The point at theta with local pitch xi = v_par/v of the orbit family
    !> with surface s_star and total energy p%energy (p%mu is ignored): s is
    !> the root of the generated residual R = rho G - (psi_p(s) - psi_p(s*)),
    !> rho = m xi v/(Z e B), v the speed from the energy at s, by Newton's
    !> method on its exact tangent in s, kept inside a bracket that grows
    !> outwards from s* on the side of the sign of rho G/iota. pt%mu is the
    !> magnetic moment of the orbit through the point, so that mu over
    !> (theta, xi) is the stream function of the orbits at fixed (s*, H);
    !> pt%mu_xi = dmu/dxi and pt%s_star_s = ds/ds* at fixed (theta, H) follow
    !> from the tangents by the implicit function theorem.
    subroutine orbit_pitch_point(bf, p, s_star, theta, xi, pt, status)
        type(boozer_file_t), intent(in) :: bf
        type(orbit_particle_t), intent(in) :: p
        real(dp), intent(in) :: s_star, theta, xi
        type(orbit_point_t), intent(out) :: pt
        integer, intent(out) :: status
        type(orbit_particle_t) :: q
        real(dp) :: lo, hi, f_lo, f_hi, s, f, f_s, step, dir, b, b_s, g, i, g_s, i_s, &
            iota, psi_prime, iota_star, psi_prime_star, phi, phi_s, r_xi, mu, mu_s, mu_xi
        integer :: it
        status = ORBIT_BAD_INPUT
        if (p%mass <= 0 .or. p%charge == 0 .or. abs(xi) > 1) return
        call orbit_state(bf, s_star, theta, b, b_s, g, i, g_s, i_s, iota_star, &
            psi_prime_star, status)
        if (status /= ORBIT_OK) return
        ! Bracket outwards from s*, where R = rho G has the sign of the root's
        ! displacement times iota psi'.
        call residual(s_star, f_lo, f_s)
        if (status /= ORBIT_OK) return
        s = s_star
        if (f_lo /= 0) then
            dir = sign(1.0_dp, f_lo*iota_star*psi_prime_star)
            lo = s_star
            step = 1e-4_dp
            do it = 1, 60
                hi = s_star + dir*step
                if (hi <= bf%d%s(1) .or. hi >= bf%d%s(size(bf%d%s))) then
                    status = ORBIT_NO_ROOT
                    return
                end if
                call residual(hi, f_hi, f_s)
                if (f_hi*f_lo <= 0) exit
                lo = hi
                f_lo = f_hi
                step = 2*step
            end do
            if (f_hi*f_lo > 0) then
                status = ORBIT_NO_ROOT
                return
            end if
            ! Newton from the first-order root, bisecting a step that leaves
            ! the bracket.
            s = s_star + f_lo/(iota_star*psi_prime_star)
            if ((s - lo)*(s - hi) > 0) s = 0.5_dp*(lo + hi)
            do it = 1, 100
                call residual(s, f, f_s)
                if (status /= ORBIT_OK) return
                if (f == 0) exit
                if (f*f_lo > 0) then
                    lo = s
                    f_lo = f
                else
                    hi = s
                end if
                step = f/f_s
                if (f_s == 0 .or. (s - step - lo)*(s - step - hi) > 0) then
                    step = s - 0.5_dp*(lo + hi)
                end if
                s = s - step
                if (abs(step) <= 4*epsilon(1.0_dp)*max(1.0_dp, abs(s))) exit
            end do
        end if
        call orbit_state(bf, s, theta, b, b_s, g, i, g_s, i_s, iota, psi_prime, status)
        if (status /= ORBIT_OK) return
        call orbit_potential(bf, p, s, phi, phi_s)
        call residual(s, f, f_s)
        call orbit_pitch_residual_jvp(xi, p%mass, p%charge, p%energy, phi, b, g, 0.0_dp, &
            1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, r_xi)
        call orbit_moment(xi, p%mass, p%charge, p%energy, phi, b, mu)
        call orbit_moment_jvp(xi, p%mass, p%charge, p%energy, phi, b, 0.0_dp, 0.0_dp, 0.0_dp, &
            0.0_dp, phi_s, b_s, mu_s)
        call orbit_moment_jvp(xi, p%mass, p%charge, p%energy, phi, b, 1.0_dp, 0.0_dp, 0.0_dp, &
            0.0_dp, 0.0_dp, 0.0_dp, mu_xi)
        q = p
        q%mu = mu
        ! Evaluated at the root itself: near a turning point the leg may pass
        ! theta at s /= s* where s* itself is out of reach.
        call orbit_evaluate(bf, q, s_star, s, theta, pt, status)
        ! ds/dxi = -R_xi/R_s; R depends on s* through -psi_p(s*) alone.
        pt%mu_xi = mu_xi - mu_s*r_xi/f_s
        pt%s_star_s = -iota_star*psi_prime_star/f_s
    contains
        !> R(s) and dR/ds from the generated leaf and its tangent along s.
        subroutine residual(s, r, r_s)
            real(dp), intent(in) :: s
            real(dp), intent(out) :: r, r_s
            real(dp) :: bb, bs, gg, ii, gs, is, io, pp, ph, ph_s
            call orbit_state(bf, s, theta, bb, bs, gg, ii, gs, is, io, pp, status)
            if (status /= ORBIT_OK) return
            call orbit_potential(bf, p, s, ph, ph_s)
            call orbit_pitch_residual(xi, p%mass, p%charge, p%energy, ph, bb, gg, &
                orbit_flux_p(bf, s_star, s), r)
            call orbit_pitch_residual_jvp(xi, p%mass, p%charge, p%energy, ph, bb, gg, 0.0_dp, &
                0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, ph_s, bs, gs, io*pp, r_s)
        end subroutine residual
    end subroutine orbit_pitch_point

    !> Poloidal and toroidal frequencies of the orbit with surface s_star and
    !> the energy and magnetic moment of p. A trapped
    !> orbit bounces between the turning points theta_1 < theta_2 where
    !> mu B(s*, theta) + Z e Phi(s*) = H; omega_b = 2 pi/tau_b and the
    !> precession is the toroidal advance per bounce over tau_b. A passing
    !> orbit (sigma given) transits theta over 2 pi. Gauss-Legendre in u with
    !> theta = theta_m + Delta sin(u) removes the inverse square root at the
    !> turning points. trapped is set from the field at s*.
    subroutine orbit_frequencies(bf, p, s_star, sigma, n_quad, omega_b, omega_zeta, trapped, &
            status)
        type(boozer_file_t), intent(in) :: bf
        type(orbit_particle_t), intent(in) :: p
        real(dp), intent(in) :: s_star, sigma
        integer, intent(in) :: n_quad
        real(dp), intent(out) :: omega_b, omega_zeta
        logical, intent(out) :: trapped
        integer, intent(out) :: status
        real(dp) :: x(n_quad), w(n_quad), t1, t2, tau, dzeta, u, theta, jac, sg
        type(orbit_point_t) :: pt
        integer :: k, leg
        call gauss_legendre(n_quad, x, w)
        call turning_points(t1, t2, trapped, status)
        if (status /= ORBIT_OK) return
        tau = 0
        dzeta = 0
        if (trapped) then
            do leg = 1, 2
                sg = merge(1.0_dp, -1.0_dp, leg == 1)
                do k = 1, n_quad
                    u = 0.5_dp*pi*x(k)
                    theta = 0.5_dp*(t1 + t2) + 0.5_dp*(t2 - t1)*sin(u)
                    jac = 0.5_dp*(t2 - t1)*cos(u)*0.5_dp*pi
                    call orbit_point(bf, p, s_star, theta, sg, pt, status)
                    if (status /= ORBIT_OK) return
                    tau = tau + w(k)*jac/abs(pt%theta_dot)
                    dzeta = dzeta + w(k)*jac*pt%zeta_dot/abs(pt%theta_dot)
                end do
            end do
        else
            do k = 1, n_quad
                theta = pi*(x(k) + 1)
                call orbit_point(bf, p, s_star, theta, sigma, pt, status)
                if (status /= ORBIT_OK) return
                tau = tau + w(k)*pi/abs(pt%theta_dot)
                dzeta = dzeta + w(k)*pi*pt%zeta_dot/abs(pt%theta_dot)
            end do
        end if
        omega_b = 2*pi/tau
        omega_zeta = dzeta/tau
        status = ORBIT_OK
    contains
        !> Roots of mu B(s*, theta) = H - Z e Phi(s*) around the minimum of B
        !> bracketed on a 512-point scan and refined by Newton on mu dB/dtheta;
        !> none means passing.
        subroutine turning_points(a, b, trap, st)
            real(dp), intent(out) :: a, b
            logical, intent(out) :: trap
            integer, intent(out) :: st
            integer, parameter :: m = 512
            real(dp) :: th(0:m), f(0:m)
            integer :: j, jmin
            do j = 0, m
                th(j) = 2*pi*j/m
                f(j) = level(th(j), st)
                if (st /= ORBIT_OK) return
            end do
            trap = any(f > 0)
            a = 0
            b = 0
            if (.not. trap) return
            jmin = minloc(f, 1) - 1
            ! Walk from the minimum of B to both sides until the level is crossed.
            a = root(th(jmin), th(jmin) - 2*pi, st)
            if (st /= ORBIT_OK) return
            b = root(th(jmin), th(jmin) + 2*pi, st)
            if (a > b) a = a - 2*pi
        end subroutine turning_points

        real(dp) function level(th, st, level_t)
            real(dp), intent(in) :: th
            integer, intent(out) :: st
            real(dp), intent(out), optional :: level_t
            real(dp) :: bb, b_s, g, i, g_s, i_s, iota, psi_prime, phi, phi_s, b_t
            call orbit_state(bf, s_star, th, bb, b_s, g, i, g_s, i_s, iota, psi_prime, st, &
                b_t=b_t)
            call orbit_potential(bf, p, s_star, phi, phi_s)
            level = p%mu*bb + p%charge*phi - p%energy
            if (present(level_t)) level_t = p%mu*b_t
        end function level

        !> First crossing of level = 0 going from inside (level < 0) at th_in
        !> towards th_out.
        real(dp) function root(th_in, th_out, st)
            real(dp), intent(in) :: th_in, th_out
            integer, intent(out) :: st
            real(dp) :: lo, hi, mid, step, f, f_t
            integer :: j
            lo = th_in
            step = (th_out - th_in)/512
            hi = lo
            do j = 1, 512
                hi = th_in + j*step
                if (level(hi, st) > 0) exit
                lo = hi
            end do
            ! Newton on level' = mu dB/dtheta inside the bracket [lo, hi].
            mid = 0.5_dp*(lo + hi)
            do j = 1, 100
                f = level(mid, st, f_t)
                if (f == 0) exit
                if (f > 0) then
                    hi = mid
                else
                    lo = mid
                end if
                step = f/f_t
                if (f_t == 0 .or. (mid - step - lo)*(mid - step - hi) > 0) &
                    step = mid - 0.5_dp*(lo + hi)
                mid = mid - step
                if (abs(step) <= 4*epsilon(1.0_dp)*max(1.0_dp, abs(mid))) exit
            end do
            root = mid
        end function root
    end subroutine orbit_frequencies

end module milca_boozer_orbit
