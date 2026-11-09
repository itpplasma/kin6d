!> Finite orbit width on the orbit-class spectral elements (milca_dk_hp): the
!! orbit geometry of one surface psi* and total energy H of a Boozer file
!! (milca_boozer_orbit), in the conventions of the thin-orbit field so that
!! it reduces to it for vanishing orbit width (plan step 3).
!!
!! At (theta, xi), xi = v_par/v at the orbit point, the pitch map gives the
!! orbit's magnetic moment mu and its theta_dot, zeta_dot. With the kinetic
!! energy E* and speed v* at s* and B_ref,
!!   stream    -mu_hat/2, mu_hat = mu B_ref/E* (= lambda for thin orbits),
!!   measure   the phase-space Jacobian of (psi*, theta, xi): 1/(b k) at the
!!             orbit point times dpsi/dpsi* (= 1/(b k) for thin orbits); by
!!             Liouville -(d mu_hat/d xi)/(2 theta_dot/v*), a ratio that turns
!!             0/0 where the drift stops the poloidal motion,
!! so that the bracket of the stream function is the streaming term in this
!! measure (Liouville at fixed psi*, H). The toroidal phase of one harmonic
!! in the frame exp(-i n zeta) is -i n zeta_dot/v* times the measure; it is
!! passed as the thin field's phase m zeta_dot/(q v*), so the operator must
!! carry omega_E = pre = 0 (the orbit's zeta_dot holds both). The drive and
!! torque weights are those of the full drift weight (wsign = +1) with b at
!! the orbit point; the perturbation harmonics are taken at s* (first stage).
!! Pitch-angle scattering acts at fixed psi* (first stage of collisions).
module milca_ntv_orbit_geometry
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use milca_boozer_file, only: boozer_file_t
    use milca_boozer_orbit, only: orbit_particle_t, orbit_point_t, orbit_pitch_point, &
        orbit_state, orbit_potential, ORBIT_OK
    use milca_dk_hp, only: hp_geometry_t, hp_operator_t
    implicit none
    private
    public :: hp_orbit_geometry_t, hp_orbit_geometry

    real(dp), parameter :: pi = 4*atan(1.0_dp)

    type, extends(hp_geometry_t) :: hp_orbit_geometry_t
        type(boozer_file_t), pointer :: bf => null()
        type(orbit_particle_t) :: p
        real(dp) :: s_star = 0, b_ref = 0, e_star = 0, v_star = 0
    contains
        procedure :: stream => orbit_stream
        procedure :: level => orbit_level
        procedure :: point => orbit_point_coefficients
        procedure :: mu_hat
        procedure :: measure
    end type

contains

    !> Geometry of the orbits with surface s_star and the total energy of p
    !> (Phi(s_ref) = 0) on bf, which must outlive g; b_ref the reference field
    !> of mu_hat. theta_x and lambda_c are those of B(s*, theta): the turning
    !> points of all orbits lie on s*.
    subroutine hp_orbit_geometry(bf, p, s_star, b_ref, g, status)
        type(boozer_file_t), target, intent(in) :: bf
        type(orbit_particle_t), intent(in) :: p
        real(dp), intent(in) :: s_star, b_ref
        type(hp_orbit_geometry_t), intent(out) :: g
        integer, intent(out) :: status
        integer, parameter :: n = 2048
        real(dp) :: b, b_s, gg, i, g_s, i_s, iota, psi_prime, th, b_max, phi, phi_s
        integer :: k
        g%bf => bf
        g%p = p
        g%s_star = s_star
        g%b_ref = b_ref
        b_max = 0
        do k = 0, n - 1
            th = 2*pi*k/n
            call orbit_state(bf, s_star, th, b, b_s, gg, i, g_s, i_s, iota, psi_prime, status)
            if (status /= ORBIT_OK) return
            if (b > b_max) then
                b_max = b
                g%theta_x = th
            end if
        end do
        ! Golden-section refinement of the maximum: the separatrix closes at
        ! theta_x, and the mesh columns end there.
        block
            real(dp) :: a, c, x1, x2, f1, f2, r
            integer :: it
            r = (sqrt(5.0_dp) - 1)/2
            a = g%theta_x - 2*pi/n
            c = g%theta_x + 2*pi/n
            x1 = c - r*(c - a)
            x2 = a + r*(c - a)
            f1 = field(x1)
            f2 = field(x2)
            do it = 1, 80
                if (f1 > f2) then
                    c = x2; x2 = x1; f2 = f1
                    x1 = c - r*(c - a); f1 = field(x1)
                else
                    a = x1; x1 = x2; f1 = f2
                    x2 = a + r*(c - a); f2 = field(x2)
                end if
            end do
            g%theta_x = modulo(0.5_dp*(a + c), 2*pi)
            b_max = field(g%theta_x)
        end block
        g%lambda_c = b_ref/b_max
        g%kappa = iota*b_ref/(gg + iota*i)
        ! Kinetic energy at s*: H - Z e Phi(s*).
        call orbit_potential(bf, p, s_star, phi, phi_s)
        g%e_star = p%energy - p%charge*phi
        g%v_star = sqrt(2*g%e_star/p%mass)
        status = ORBIT_OK
    contains
        real(dp) function field(th)
            real(dp), intent(in) :: th
            real(dp) :: bb, bs, g1, i1, gs, is, io, pp
            integer :: st
            call orbit_state(bf, s_star, th, bb, bs, g1, i1, gs, is, io, pp, st)
            field = bb
        end function field
    end subroutine hp_orbit_geometry

    real(dp) function mu_hat(g, th, xi, pt, status)
        class(hp_orbit_geometry_t), intent(in) :: g
        real(dp), intent(in) :: th, xi
        type(orbit_point_t), intent(out) :: pt
        integer, intent(out) :: status
        call orbit_pitch_point(g%bf, g%p, g%s_star, th, xi, pt, status)
        mu_hat = pt%mu*g%b_ref/g%e_star
    end function mu_hat

    !> The measure at (theta, xi) as the phase-space Jacobian of (psi*,
    !> theta, xi) at fixed H: the thin field's 1/(b k) at the orbit point
    !> times dpsi/dpsi* at fixed (theta, xi) and v/v*. By Liouville it equals
    !> -(d mu_hat/d xi)/(2 theta_dot/v*), but stays regular where theta_dot
    !> and d mu_hat/d xi vanish together.
    real(dp) function measure(g, th, xi, op) result(m)
        class(hp_orbit_geometry_t), intent(in) :: g
        real(dp), intent(in) :: th, xi
        type(hp_operator_t), intent(in) :: op
        type(orbit_point_t) :: pt, up, down
        type(hp_orbit_geometry_t) :: shifted
        real(dp), parameter :: h = 1e-6_dp
        real(dp) :: b, b_s, gg, i, g_s, i_s, iota, psi_prime, iota_star, k, v, ds
        integer :: st
        b = g%mu_hat(th, xi, pt, st)
        shifted = g
        shifted%s_star = g%s_star + h
        b = shifted%mu_hat(th, xi, up, st)
        shifted%s_star = g%s_star - h
        b = shifted%mu_hat(th, xi, down, st)
        ds = (up%s - down%s)/(2*h)
        call orbit_state(g%bf, g%s_star, th, b, b_s, gg, i, g_s, i_s, iota_star, psi_prime, &
            st)
        call orbit_state(g%bf, pt%s, th, b, b_s, gg, i, g_s, i_s, iota, psi_prime, st)
        b = b/g%b_ref
        k = op%k*b**op%k_power*(iota*(g%kappa/op%k)/iota_star)
        v = abs(pt%v_par)/max(abs(xi), tiny(1.0_dp))
        if (abs(xi) < 1e-12_dp) v = g%v_star
        ! dpsi/dpsi* with psi = int iota psi' ds: iota ds/(iota* ds*).
        m = 1/(b*k)*ds*iota/iota_star*v/g%v_star
    end function measure

    real(dp) function orbit_stream(g, th, xi) result(psi)
        class(hp_orbit_geometry_t), intent(in) :: g
        real(dp), intent(in) :: th, xi
        type(orbit_point_t) :: pt
        integer :: st
        psi = -g%mu_hat(th, max(-1.0_dp, min(1.0_dp, xi)), pt, st)/2
    end function orbit_stream

    !> Passing level mu_hat = val on the side kind = +-1, or the scaled
    !> separatrix xi = val xi_c(theta) on the side of val, xi_c the root of
    !> mu_hat = lambda_c; by bisection in xi, where mu_hat falls from its value
    !> on s* at xi = 0 to zero at xi = +-1.
    subroutine orbit_level(g, kind, val, th, xi, dxi)
        class(hp_orbit_geometry_t), intent(in) :: g
        integer, intent(in) :: kind
        real(dp), intent(in) :: val, th
        real(dp), intent(out) :: xi, dxi
        dxi = 0
        if (kind /= 0) then
            xi = root(real(kind, dp), val)
        else if (val == 0) then
            xi = 0
        else
            xi = abs(val)*root(sign(1.0_dp, val), g%lambda_c)
            xi = sign(xi, val)
        end if
    contains
        !> |xi| on the side sigma where mu_hat = level; zero where mu_hat on s*
        !> is already below the level (beyond the separatrix X point).
        real(dp) function root(sigma, level) result(x)
            real(dp), intent(in) :: sigma, level
            type(orbit_point_t) :: pt
            real(dp) :: lo, hi, mid
            integer :: it, st
            x = 0
            if (g%mu_hat(th, 0.0_dp, pt, st) <= level) return
            lo = 0
            hi = 1
            do it = 1, 60
                mid = 0.5_dp*(lo + hi)
                if (g%mu_hat(th, sigma*mid, pt, st) > level) then
                    lo = mid
                else
                    hi = mid
                end if
            end do
            x = sigma*0.5_dp*(lo + hi)
        end function root
    end subroutine orbit_level

    subroutine orbit_point_coefficients(g, th, xi, op, collide, drive_r, drive_e, density, &
            phase, precess, drive, torque)
        class(hp_orbit_geometry_t), intent(in) :: g
        real(dp), intent(in) :: th, xi
        type(hp_operator_t), intent(in) :: op
        real(dp), intent(out) :: collide, drive_r, drive_e, density, phase, precess, drive, &
            torque
        type(orbit_point_t) :: pt
        real(dp) :: m, b, k
        integer :: st
        b = g%mu_hat(th, xi, pt, st)
        m = g%measure(th, xi, op)
        b = pt%b/g%b_ref
        k = op%k*b**op%k_power
        collide = m*op%nu*(1 - xi*xi)/2
        drive_r = 0
        drive_e = 0
        density = m
        ! The axisymmetric operator carries no q; the phase enters only for n /= 0.
        phase = 0
        if (op%q /= 0) phase = m*pt%zeta_dot/(op%q*g%v_star)
        precess = 0
        drive = m*(1 + xi*xi)/b
        torque = 2*k*drive
    end subroutine orbit_point_coefficients

end module milca_ntv_orbit_geometry
