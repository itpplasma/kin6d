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
    use milca_boozer_file, only: boozer_file_t, boozer_field_t, boozer_field, boozer_field_eval, &
        boozer_harmonics, BOOZER_OK
    use milca_boozer_orbit, only: orbit_particle_t, orbit_point_t, orbit_pitch_point, &
        orbit_point, orbit_state, orbit_potential, ORBIT_OK, ORBIT_NO_ROOT
    use milca_dk_hp, only: hp_geometry_t, hp_operator_t
    implicit none
    private
    public :: hp_orbit_geometry_t, hp_orbit_geometry, hp_orbit_perturbation, hp_orbit_confined

    real(dp), parameter :: pi = 4*atan(1.0_dp)

    type, extends(hp_geometry_t) :: hp_orbit_geometry_t
        type(boozer_file_t), pointer :: bf => null()
        ! Perturbation file of the toroidal harmonic n_tor; when associated, the
        ! drive is taken at the orbit point instead of from the operator.
        type(boozer_file_t), pointer :: pert => null()
        integer :: n_tor = 0
        real(dp) :: pert_ref = 0
        ! The harmonics of n_tor alone, cut from the perturbation file.
        type(boozer_field_t) :: pert_field
        type(orbit_particle_t) :: p
        real(dp) :: s_star = 0, b_ref = 0, e_star = 0, v_star = 0
    contains
        procedure :: stream => orbit_stream
        procedure :: level => orbit_level
        procedure :: point => orbit_point_coefficients
        procedure :: mu_hat
        procedure :: measure
        procedure :: perturbation => orbit_perturbation
        procedure :: point_source => orbit_point_source
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
        ! Newton on dB/dtheta = 0 from the sampled maximum: the separatrix
        ! closes at theta_x, and the mesh columns end there.
        block
            real(dp) :: b_t, b_tt, step
            integer :: it
            do it = 1, 50
                call orbit_state(bf, s_star, g%theta_x, b, b_s, gg, i, g_s, i_s, iota, &
                    psi_prime, status, b_t=b_t, b_tt=b_tt)
                if (b_tt >= 0) exit
                step = b_t/b_tt
                g%theta_x = g%theta_x - max(-pi/n, min(pi/n, step))
                if (abs(step) <= 4*epsilon(1.0_dp)) exit
            end do
            g%theta_x = modulo(g%theta_x, 2*pi)
            call orbit_state(bf, s_star, g%theta_x, b_max, b_s, gg, i, g_s, i_s, iota, &
                psi_prime, status)
        end block
        g%lambda_c = b_ref/b_max
        g%kappa = iota*b_ref/(gg + iota*i)
        ! Kinetic energy at s*: H - Z e Phi(s*).
        call orbit_potential(bf, p, s_star, phi, phi_s)
        g%e_star = p%energy - p%charge*phi
        g%v_star = sqrt(2*g%e_star/p%mass)
        status = ORBIT_OK
    end subroutine hp_orbit_geometry

    !> Takes the drive of harmonic n_tor from pert at the orbit point; b_ref
    !> normalises the amplitudes as the thin-orbit device does (c_m = (bmnc +
    !> i bmns)/B_ref of mode -m).
    subroutine hp_orbit_perturbation(g, pert, n_tor, b_ref)
        type(hp_orbit_geometry_t), intent(inout) :: g
        type(boozer_file_t), target, intent(in) :: pert
        integer, intent(in) :: n_tor
        real(dp), intent(in) :: b_ref
        g%pert => pert
        g%n_tor = n_tor
        g%pert_ref = b_ref
        call boozer_field(pert, n_tor, g%pert_field)
    end subroutine hp_orbit_perturbation

    !> Whether every orbit of g through an n_theta x n_xi sample of (theta, xi)
    !> has its point inside the files' radial range. Orbits that leave it
    !> (near the axis or at high energy) are outside the model; the kernel
    !> must not be solved for them.
    logical function hp_orbit_confined(g, n_theta, n_xi) result(ok)
        type(hp_orbit_geometry_t), intent(in) :: g
        integer, intent(in) :: n_theta, n_xi
        type(orbit_point_t) :: pt
        real(dp) :: mu, th, xi
        integer :: j, k, st
        ok = .false.
        do j = 0, n_theta - 1
            th = g%theta_x - 2*pi + 2*pi*(j + 0.5_dp)/n_theta
            do k = 0, n_xi - 1
                xi = -1 + 2*(k + 0.5_dp)/n_xi
                mu = g%mu_hat(th, xi, pt, st)
                if (st /= ORBIT_OK) return
                if (associated(g%pert)) then
                    if (pt%s < g%pert%d%s(1) .or. pt%s > g%pert%d%s(g%pert%d%nsurf)) return
                end if
            end do
        end do
        ok = .true.
    end function hp_orbit_confined

    complex(dp) function orbit_perturbation(g, th, xi, op) result(src)
        class(hp_orbit_geometry_t), intent(in) :: g
        real(dp), intent(in) :: th, xi
        type(hp_operator_t), intent(in) :: op
        type(orbit_point_t) :: pt
        real(dp) :: mu
        integer :: st
        if (.not. associated(g%pert)) then
            src = op_drive(op, th)
            return
        end if
        src = 0
        mu = g%mu_hat(th, xi, pt, st)
        if (st /= ORBIT_OK) return
        src = perturbation_at(g, th, op, pt)
    contains
        complex(dp) function op_drive(o, t)
            type(hp_operator_t), intent(in) :: o
            real(dp), intent(in) :: t
            integer :: j
            op_drive = 0
            if (.not. allocated(o%m)) return
            do j = 1, size(o%m)
                op_drive = op_drive + (0, 1)*o%n*o%c(j)*exp((0, 1)*o%m(j)*t)
            end do
        end function op_drive
    end function orbit_perturbation

    !> The drive of n_tor at the pitch point pt of theta (orbit_perturbation).
    complex(dp) function perturbation_at(g, th, op, pt) result(src)
        class(hp_orbit_geometry_t), intent(in) :: g
        real(dp), intent(in) :: th
        type(hp_operator_t), intent(in) :: op
        type(orbit_point_t), intent(in) :: pt
        real(dp), dimension(size(g%pert_field%m)) :: bmn, bma, bmn_s, bma_s, cm, sm
        real(dp) :: iota, gc, ic, gc_s, ic_s
        integer :: k, st
        src = 0
        call boozer_field_eval(g%pert_field, pt%s, iota, gc, ic, gc_s, ic_s, bmn, bma, bmn_s, &
            bma_s, st)
        if (st /= BOOZER_OK) return
        call boozer_harmonics(g%pert_field%m, th, cm, sm)
        do k = 1, size(g%pert_field%m)
            src = src + (0, 1)*op%n*cmplx(bmn(k), bma(k), dp)/g%pert_ref &
                *cmplx(cm(k), -sm(k), dp)
        end do
    end function perturbation_at

    !> orbit_perturbation and orbit_point_coefficients on one pitch point.
    subroutine orbit_point_source(g, th, xi, op, src, collide, drive_r, drive_e, density, &
            phase, precess, drive, torque)
        class(hp_orbit_geometry_t), intent(in) :: g
        real(dp), intent(in) :: th, xi
        type(hp_operator_t), intent(in) :: op
        complex(dp), intent(out) :: src
        real(dp), intent(out) :: collide, drive_r, drive_e, density, phase, precess, drive, &
            torque
        type(orbit_point_t) :: pt
        real(dp) :: mu
        integer :: st
        mu = g%mu_hat(th, xi, pt, st)
        if (.not. associated(g%pert)) then
            src = g%perturbation(th, xi, op)
        else if (st /= ORBIT_OK) then
            src = 0
        else
            src = perturbation_at(g, th, op, pt)
        end if
        call point_at(g, th, xi, op, pt, collide, drive_r, drive_e, density, phase, precess, &
            drive, torque)
    end subroutine orbit_point_source

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
    !> times dpsi/dpsi* at fixed (theta, xi) and v/v*, with ds/ds* from the
    !> implicit function theorem on the generated pitch residual. By Liouville it equals
    !> -(d mu_hat/d xi)/(2 theta_dot/v*), but stays regular where theta_dot
    !> and d mu_hat/d xi vanish together.
    real(dp) function measure(g, th, xi, op) result(m)
        class(hp_orbit_geometry_t), intent(in) :: g
        real(dp), intent(in) :: th, xi
        type(hp_operator_t), intent(in) :: op
        type(orbit_point_t) :: pt
        real(dp) :: b
        integer :: st
        b = g%mu_hat(th, xi, pt, st)
        m = measure_at(g, th, xi, op, pt)
    end function measure

    !> The measure at the pitch point pt of (theta, xi).
    real(dp) function measure_at(g, th, xi, op, pt) result(m)
        class(hp_orbit_geometry_t), intent(in) :: g
        real(dp), intent(in) :: th, xi
        type(hp_operator_t), intent(in) :: op
        type(orbit_point_t), intent(in) :: pt
        real(dp) :: b, b_s, gg, i, g_s, i_s, iota, psi_prime, iota_star, k, v, ds
        integer :: st
        ds = pt%s_star_s
        call orbit_state(g%bf, g%s_star, th, b, b_s, gg, i, g_s, i_s, iota_star, psi_prime, &
            st)
        call orbit_state(g%bf, pt%s, th, b, b_s, gg, i, g_s, i_s, iota, psi_prime, st)
        b = b/g%b_ref
        k = op%k*b**op%k_power*(iota*(g%kappa/op%k)/iota_star)
        v = abs(pt%v_par)/max(abs(xi), tiny(1.0_dp))
        if (abs(xi) < 1e-12_dp) v = g%v_star
        ! dpsi/dpsi* with psi = int iota psi' ds: iota ds/(iota* ds*).
        m = 1/(b*k)*ds*iota/iota_star*v/g%v_star
    end function measure_at

    real(dp) function orbit_stream(g, th, xi) result(psi)
        class(hp_orbit_geometry_t), intent(in) :: g
        real(dp), intent(in) :: th, xi
        type(orbit_point_t) :: pt
        integer :: st
        psi = -g%mu_hat(th, max(-1.0_dp, min(1.0_dp, xi)), pt, st)/2
    end function orbit_stream

    !> Passing level mu_hat = val on the side kind = +-1, or the scaled
    !> separatrix xi = val xi_c(theta) on the side of val, xi_c the root of
    !> mu_hat = lambda_c; by Newton in xi, where mu_hat falls from its value
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
        !> xi on the side sigma where mu_hat = level. At fixed mu, H and P_zeta
        !> the orbit point at theta is one root in s (orbit_point), and xi is
        !> v_par/v there; it does not reach theta beyond the X point (xi = 0).
        !> The Newton iteration in xi on nested pitch points is the fallback.
        real(dp) function root(sigma, level) result(x)
            real(dp), intent(in) :: sigma, level
            type(orbit_particle_t) :: q
            type(orbit_point_t) :: pt
            real(dp) :: phi, phi_s
            integer :: st
            x = sigma
            if (level <= 0) return
            q = g%p
            q%mu = level*g%e_star/g%b_ref
            call orbit_point(g%bf, q, g%s_star, th, sigma, pt, st)
            if (st == ORBIT_OK) then
                call orbit_potential(g%bf, q, pt%s, phi, phi_s)
                x = pt%v_par/sqrt(2*(q%energy - q%charge*phi)/q%mass)
            else if (st == ORBIT_NO_ROOT .and. g%mu_hat(th, 0.0_dp, pt, st) <= level) then
                x = 0
            else
                x = root_newton(sigma, level)
            end if
        end function root

        !> |xi| on the side sigma where mu_hat = level; zero where mu_hat on s*
        !> is already below the level (beyond the separatrix X point). mu_hat
        !> falls monotonically to zero at |xi| = 1, so [0, 1] brackets the root;
        !> Newton on the exact dmu_hat/dxi of the pitch point, bisecting a step
        !> that leaves the bracket.
        real(dp) function root_newton(sigma, level) result(x)
            real(dp), intent(in) :: sigma, level
            type(orbit_point_t) :: pt
            real(dp) :: lo, hi, f, f_x, step
            integer :: it, st
            x = 0
            if (g%mu_hat(th, 0.0_dp, pt, st) <= level) return
            x = sigma
            if (level <= 0) return
            lo = 0
            hi = 1
            ! First guess from mu_hat ~ mu_hat(0) (1 - xi^2).
            x = sqrt(max(0.0_dp, 1 - level/g%mu_hat(th, 0.0_dp, pt, st)))
            do it = 1, 100
                f = g%mu_hat(th, sigma*x, pt, st) - level
                f_x = sigma*pt%mu_xi*g%b_ref/g%e_star
                if (f == 0) exit
                if (f > 0) then
                    lo = x
                else
                    hi = x
                end if
                step = f/f_x
                if (f_x == 0 .or. (x - step - lo)*(x - step - hi) > 0) step = x - 0.5_dp*(lo + hi)
                x = x - step
                ! mu_hat carries the rounding of its own root, so steps of a
                ! few epsilon persist near the X point, where f_x is small.
                if (abs(step) <= 64*epsilon(1.0_dp) .or. hi - lo <= 64*epsilon(1.0_dp)) exit
            end do
            x = sigma*x
        end function root_newton
    end subroutine orbit_level

    subroutine orbit_point_coefficients(g, th, xi, op, collide, drive_r, drive_e, density, &
            phase, precess, drive, torque)
        class(hp_orbit_geometry_t), intent(in) :: g
        real(dp), intent(in) :: th, xi
        type(hp_operator_t), intent(in) :: op
        real(dp), intent(out) :: collide, drive_r, drive_e, density, phase, precess, drive, &
            torque
        type(orbit_point_t) :: pt
        real(dp) :: b
        integer :: st
        b = g%mu_hat(th, xi, pt, st)
        call point_at(g, th, xi, op, pt, collide, drive_r, drive_e, density, phase, precess, &
            drive, torque)
    end subroutine orbit_point_coefficients

    !> The coefficients of orbit_point_coefficients at the pitch point pt.
    subroutine point_at(g, th, xi, op, pt, collide, drive_r, drive_e, density, phase, &
            precess, drive, torque)
        class(hp_orbit_geometry_t), intent(in) :: g
        real(dp), intent(in) :: th, xi
        type(hp_operator_t), intent(in) :: op
        type(orbit_point_t), intent(in) :: pt
        real(dp), intent(out) :: collide, drive_r, drive_e, density, phase, precess, drive, &
            torque
        real(dp) :: m, b, k
        m = measure_at(g, th, xi, op, pt)
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
        ! As the thin field: the flux weight with the surface measure 1/b**k_power.
        torque = 2*k*drive/b**op%k_power
    end subroutine point_at

end module milca_ntv_orbit_geometry
