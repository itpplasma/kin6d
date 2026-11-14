!> One-dimensional tokamak discharge for teaching: operating limits and the
!> time history of a disruption.
!>
!> The plasma is a straight cylinder of area-equivalent radius a_e = a sqrt(kappa);
!> shaping and toroidicity enter q through the ITER Physics Basis q95 factor.
!> State: electron = ion temperature T(x,t) and enclosed current I(x,t), x = r/a_e.
!> - heat: 3n dT/dt = div(2 n chi grad T) + eta j^2 + p_aux - n^2 L(T), with chi
!>   scaled to the IPB98(y,2) confinement time;
!> - current: resistive diffusion dI/dt = (x/(mu0 a_e^2)) d/dx((eta/x) dI/dx) with
!>   Spitzer eta(T); the edge holds I_p (feedback) until the thermal quench and
!>   then follows the external inductance;
!> - islands: Rutherford growth from the cylindrical Delta' of the 2/1 and 3/2
!>   surfaces, which flatten T across their width;
!> - ideal limits: cylindrical external kink (m = 2..4) and local s-alpha
!>   ballooning; their growth is an estimate from the energy principle;
!> - disruption: thermal quench by stochastic transport, current flattening with
!>   poloidal-flux conservation (current spike), current quench with cold
!>   resistivity, and a rigid vertical displacement on the wall time.
!> Surrogates are marked where they appear; nothing here is a validated
!> disruption predictor.
module milca_tokamak_discharge
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: iso_c_binding, only: c_double, c_int
    use milca_cylinder_stability, only: cyl_current_from_q, cyl_kink_energy, &
        cyl_tearing_delta_prime, CYL_OK, CYL_RESONANT_INSIDE
    use milca_ballooning_s_alpha, only: ballooning_critical_alpha, &
        ballooning_first_zero, BALLOON_OK
    implicit none
    private

    integer, parameter, public :: TD_NX = 65, TD_NPAR = 20, TD_NOUT = 40, &
        TD_NSER = 28, TD_NPROF = 5, TD_NRAMP = 8
    ! Phases of a discharge
    integer, parameter, public :: TD_FLATTOP = 0, TD_TEARING = 1, TD_IDEAL = 2, &
        TD_THERMAL_QUENCH = 3, TD_CURRENT_QUENCH = 4, TD_TERMINATED = 5, &
        TD_SURVIVED = 6
    ! Causes of a thermal quench
    integer, parameter, public :: TD_CAUSE_NONE = 0, TD_CAUSE_ISLANDS = 1, &
        TD_CAUSE_KINK = 2, TD_CAUSE_BALLOONING = 3, TD_CAUSE_RADIATIVE = 4

    public :: td_limits, td_run, milca_td_limits, milca_td_run, milca_td_sizes

    real(dp), parameter :: pi = acos(-1.0_dp), mu0 = 4e-7_dp*pi, &
        kev = 1.602176634e-16_dp, mp = 1.67262192e-27_dp, lnlam = 17
    integer, parameter :: NBAL = 30
    real(dp), save :: bal_s(NBAL), bal_alpha(NBAL)
    logical, save :: bal_ready = .false.

    !> Machine and plasma: everything one discharge needs.
    type :: plasma_t
        real(dp) :: r0, a, kappa, delta, b0, ip, q0, nbar, npeak, t0, tpeak
        real(dp) :: zeff, fimp, wall, tauw, h98, aion, tedge
        real(dp) :: ae, cq, lext, chi0, chi_ref_tau, paux
        real(dp) :: x(TD_NX), t(TD_NX), cur(TD_NX), n(TD_NX), heat(TD_NX)
    end type

contains

    ! ------------------------------------------------------------------ setup

    subroutine plasma_from_par(par, p, status)
        real(dp), intent(in) :: par(TD_NPAR)
        type(plasma_t), intent(out) :: p
        integer, intent(out) :: status
        real(dp) :: eps, nu, q95, s
        integer :: i
        status = 1
        p%r0 = par(1); p%a = par(2); p%kappa = par(3); p%delta = par(4)
        p%b0 = par(5); p%ip = par(6)*1e6_dp; p%q0 = par(7)
        p%nbar = par(8)*1e20_dp; p%npeak = par(9); p%t0 = par(10)
        p%tpeak = par(11); p%zeff = par(12); p%fimp = par(13); p%wall = par(14)
        p%tauw = par(15); p%h98 = par(16); p%aion = par(17); p%tedge = par(18)
        if (p%r0 <= 0 .or. p%a <= 0 .or. p%a >= p%r0 .or. p%kappa < 1 .or. &
            p%b0 <= 0 .or. p%ip <= 0 .or. p%q0 <= 0 .or. p%nbar <= 0 .or. &
            p%t0 <= p%tedge .or. p%tedge <= 0 .or. p%zeff < 1 .or. p%fimp < 0 &
            .or. p%wall <= 1 .or. p%tauw <= 0 .or. p%h98 <= 0 .or. p%aion <= 0) return
        eps = p%a/p%r0
        p%ae = p%a*sqrt(p%kappa)
        s = (1 + p%kappa**2*(1 + 2*p%delta**2 - 1.2_dp*p%delta**3))/2 &
            *(1.17_dp - 0.65_dp*eps)/(1 - eps*eps)**2
        p%cq = 2*pi*p%a**2*p%b0/(mu0*p%r0)*s
        p%lext = mu0*p%r0*max(0.1_dp, log(8*p%r0/p%ae) - 2)
        q95 = p%cq/p%ip
        nu = max(0.05_dp, q95/p%q0 - 1)
        do i = 1, TD_NX
            p%x(i) = real(i-1, dp)/real(TD_NX-1, dp)
            p%cur(i) = p%ip*(1 - (1 - p%x(i)**2)**(nu + 1))
            p%t(i) = p%tedge + (p%t0 - p%tedge)*(1 - p%x(i)**2)**p%tpeak
            p%heat(i) = exp(-((p%x(i) - 0.25_dp)/0.3_dp)**2)
        end do
        call set_density(p, p%nbar)
        status = 0
    end subroutine

    !> n = n_edge + (n_0 - n_edge)(1 - x^2)^peaking with n_edge = 0.3 nbar
    !> rescaled to the line average nbar.
    subroutine set_density(p, nbar)
        type(plasma_t), intent(inout) :: p
        real(dp), intent(in) :: nbar
        real(dp) :: avg
        integer :: i
        do i = 1, TD_NX
            p%n(i) = 0.3_dp + 0.7_dp*(1 - p%x(i)**2)**p%npeak
        end do
        avg = line_average(p%n)
        p%n = p%n*nbar/avg
        p%nbar = nbar
    end subroutine

    pure real(dp) function line_average(f) result(v)
        real(dp), intent(in) :: f(TD_NX)
        v = (sum(f) - (f(1) + f(TD_NX))/2)/real(TD_NX-1, dp)
    end function

    !> Volume of the cell around node i (2 pi R0 . 2 pi a_e^2 x dx).
    pure real(dp) function cell_volume(p, i) result(v)
        type(plasma_t), intent(in) :: p
        integer, intent(in) :: i
        real(dp) :: h
        h = 1.0_dp/real(TD_NX-1, dp)
        if (i == 1) then
            v = h*h/8
        else if (i == TD_NX) then
            v = (1 - h/4)*h/2
        else
            v = p%x(i)*h
        end if
        v = v*4*pi*pi*p%r0*p%ae**2
    end function

    pure real(dp) function spitzer(p, t) result(eta)
        type(plasma_t), intent(in) :: p
        real(dp), intent(in) :: t
        real(dp) :: fz
        fz = p%zeff*(0.29_dp + 0.457_dp/(1.077_dp + p%zeff)) &
            /(0.29_dp + 0.457_dp/2.077_dp)
        eta = 1.65e-9_dp*lnlam*fz/max(t, 1e-3_dp)**1.5_dp
    end function

    !> Radiated power density / n_e^2 [W m^3]: bremsstrahlung and a carbon-like
    !> impurity cooling curve (surrogate: peak 1.5e-31 W m^3 at 10 eV, T^-1
    !> above and T^3 below the peak).
    pure subroutine radiation(p, t, l, dl)
        type(plasma_t), intent(in) :: p
        real(dp), intent(in) :: t
        real(dp), intent(out) :: l, dl
        real(dp) :: u, den, dden
        u = t/0.010_dp
        den = u + u**(-3)
        dden = (1 - 3*u**(-4))/0.010_dp
        l = p%fimp*2*1.5e-31_dp/den + 5.35e-37_dp*p%zeff*sqrt(t)
        dl = -p%fimp*2*1.5e-31_dp*dden/den**2 + 0.5_dp*5.35e-37_dp*p%zeff/sqrt(t)
    end subroutine

    pure subroutine q_profile(p, q)
        type(plasma_t), intent(in) :: p
        real(dp), intent(out) :: q(TD_NX)
        integer :: i
        do i = 2, TD_NX
            q(i) = p%cq*p%x(i)**2/max(p%cur(i), 1e-12_dp*p%ip)
        end do
        q(1) = 2*q(2) - q(3)
        if (q(1) <= 0) q(1) = q(2)
    end subroutine

    !> Current density [A/m^2] on the nodes.
    pure subroutine current_density(p, j)
        type(plasma_t), intent(in) :: p
        real(dp), intent(out) :: j(TD_NX)
        real(dp) :: h
        integer :: i
        h = 1.0_dp/real(TD_NX-1, dp)
        do i = 2, TD_NX-1
            j(i) = (p%cur(i+1) - p%cur(i-1))/(2*h*p%x(i))
        end do
        j(1) = 2*p%cur(2)/(p%x(2)**2)/2
        j(TD_NX) = (p%cur(TD_NX) - p%cur(TD_NX-1))/h
        j = j/(2*pi*p%ae**2)
    end subroutine

    !> Internal inductance li = 2 int (I(x)/x)^2 x dx / I_a^2.
    pure real(dp) function internal_inductance(p) result(li)
        type(plasma_t), intent(in) :: p
        real(dp) :: h, f
        integer :: i
        h = 1.0_dp/real(TD_NX-1, dp)
        li = 0
        do i = 2, TD_NX
            f = p%cur(i)**2/p%x(i)
            if (i == TD_NX) f = f/2
            li = li + 2*f*h
        end do
        li = li/max(p%cur(TD_NX), 1e-30_dp)**2
    end function

    pure real(dp) function stored_energy(p) result(w)
        type(plasma_t), intent(in) :: p
        integer :: i
        w = 0
        do i = 1, TD_NX
            w = w + 3*p%n(i)*p%t(i)*kev*cell_volume(p, i)
        end do
    end function

    pure real(dp) function plasma_volume(p) result(v)
        type(plasma_t), intent(in) :: p
        v = 2*pi*pi*p%r0*p%ae**2
    end function

    !> IPB98(y,2) energy confinement time [s] at heating power p_mw [MW].
    pure real(dp) function tau98(p, ip, p_mw) result(tau)
        type(plasma_t), intent(in) :: p
        real(dp), intent(in) :: ip, p_mw
        tau = 0.0562_dp*p%h98*(ip/1e6_dp)**0.93_dp*p%b0**0.15_dp &
            *max(p_mw, 0.01_dp)**(-0.69_dp)*(p%nbar/1e19_dp)**0.41_dp &
            *p%aion**0.19_dp*p%r0**1.97_dp*(p%a/p%r0)**0.58_dp*p%kappa**0.78_dp
    end function

    !> Ohmic, auxiliary and radiated power [W] and loop voltage [V].
    subroutine powers(p, pohm, paux, prad, vloop)
        type(plasma_t), intent(in) :: p
        real(dp), intent(out) :: pohm, paux, prad, vloop
        real(dp) :: j(TD_NX), l, dl, norm, v
        integer :: i
        call current_density(p, j)
        pohm = 0; prad = 0; norm = 0
        do i = 1, TD_NX
            v = cell_volume(p, i)
            call radiation(p, p%t(i), l, dl)
            pohm = pohm + spitzer(p, p%t(i))*j(i)**2*v
            prad = prad + p%n(i)**2*l*v
            norm = norm + p%heat(i)*v
        end do
        paux = p%paux
        vloop = 2*pi*p%r0*spitzer(p, p%t(TD_NX))*j(TD_NX)
    end subroutine

    !> Calibrates the auxiliary power and the transport level so that the
    !> initial profiles are close to a steady state with tau_E = H tau98:
    !> P = (W/C)^(1/0.31) from tau98 = C P^-0.69, and chi at x = 0.75 carries
    !> the net power deposited inside.
    subroutine calibrate(p)
        type(plasma_t), intent(inout) :: p
        real(dp) :: w, c, pneed, pohm, paux, prad, vloop, pin, j(TD_NX), l, dl
        real(dp) :: v, norm, grad, h
        integer :: i, ic
        w = stored_energy(p)
        p%paux = 0
        call powers(p, pohm, paux, prad, vloop)
        c = tau98(p, p%ip, 1.0_dp)
        pneed = (w/1e6_dp/c)**(1/0.31_dp)*1e6_dp
        p%paux = max(0.0_dp, pneed - pohm)
        p%chi_ref_tau = tau98(p, p%ip, (pohm + p%paux)/1e6_dp)
        call current_density(p, j)
        norm = 0
        do i = 1, TD_NX
            norm = norm + p%heat(i)*cell_volume(p, i)
        end do
        h = 1.0_dp/real(TD_NX-1, dp)
        ic = nint(0.75_dp*real(TD_NX-1, dp)) + 1
        pin = 0
        do i = 1, ic
            v = cell_volume(p, i)
            call radiation(p, p%t(i), l, dl)
            pin = pin + (spitzer(p, p%t(i))*j(i)**2 + p%paux*p%heat(i)/norm &
                - p%n(i)**2*l)*v
        end do
        grad = abs(p%t(ic+1) - p%t(ic-1))/(2*h)*kev/p%ae
        p%chi0 = 0.05_dp
        if (pin > 0 .and. grad > 0) p%chi0 = pin/(4*pi*pi*p%r0*p%ae*p%x(ic) &
            *2*p%n(ic)*grad*(1 + 3*p%x(ic)**2))
        p%chi0 = min(max(p%chi0, 0.02_dp), 50.0_dp)
    end subroutine

    ! -------------------------------------------------------------- stability

    !> Table of the first s-alpha ballooning boundary alpha_c(s), s = 0.1..3.
    subroutine ballooning_table()
        integer :: k, rc
        real(dp) :: ac
        if (bal_ready) return
        do k = 1, NBAL
            bal_s(k) = 0.1_dp*real(k, dp)
            call ballooning_critical_alpha(bal_s(k), 4.0_dp, 12*pi, 768, ac, rc)
            if (rc /= BALLOON_OK) ac = 4
            bal_alpha(k) = ac
        end do
        bal_ready = .true.
    end subroutine

    pure real(dp) function alpha_crit(s) result(ac)
        real(dp), intent(in) :: s
        real(dp) :: u, w
        integer :: k
        if (s <= bal_s(1)) then
            ac = bal_alpha(1)*max(s, 0.0_dp)/bal_s(1)
            return
        end if
        u = min(s, bal_s(NBAL))/0.1_dp
        k = min(int(u), NBAL-1)
        w = u - real(k, dp)
        ac = (1-w)*bal_alpha(k) + w*bal_alpha(k+1)
    end function

    !> alpha/alpha_c on 0.2 <= x <= 0.95 (zero elsewhere); alpha = -2 mu0 R0 q^2 p'/B0^2.
    !> The axis is left out: there s -> 0 and the local s-alpha model does not apply.
    subroutine ballooning_ratio(p, q, ratio, shear, drive)
        type(plasma_t), intent(in) :: p
        real(dp), intent(in) :: q(TD_NX)
        real(dp), intent(out) :: ratio(TD_NX)
        real(dp), intent(out), optional :: shear(TD_NX), drive(TD_NX)
        real(dp) :: h, dp_dr, s, alpha
        integer :: i
        h = 1.0_dp/real(TD_NX-1, dp)
        ratio = 0
        if (present(shear)) shear = 0
        if (present(drive)) drive = 0
        do i = 2, TD_NX-1
            if (p%x(i) < 0.2_dp .or. p%x(i) > 0.95_dp) cycle
            dp_dr = 2*kev*(p%n(i+1)*p%t(i+1) - p%n(i-1)*p%t(i-1))/(2*h*p%ae)
            ! Elongation raises the circular s-alpha boundary; surrogate factor
            ! (1 + kappa^2)/2 on alpha.
            alpha = -2*mu0*p%r0*q(i)**2*dp_dr/p%b0**2/((1 + p%kappa**2)/2)
            s = p%x(i)*(q(i+1) - q(i-1))/(2*h*q(i))
            if (alpha > 0) ratio(i) = alpha/max(alpha_crit(s), 1e-3_dp)
            if (present(shear)) shear(i) = s
            if (present(drive)) drive(i) = alpha
        end do
    end subroutine

    pure real(dp) function alfven_time(p) result(tau)
        type(plasma_t), intent(in) :: p
        tau = p%r0*sqrt(mu0*line_average(p%n)*p%aion*mp)/p%b0
    end function

    ! ----------------------------------------------------------- static limits

    !> Operating point and stability of the parameter set (see TD_NOUT layout
    !> in docs/tokamak-limits-course.md); prof(:,k) = T, n, q, j, alpha/alpha_c.
    subroutine td_limits(par, out, prof, status)
        real(dp), intent(in) :: par(TD_NPAR)
        real(dp), intent(out) :: out(TD_NOUT), prof(TD_NX, TD_NPROF)
        integer, intent(out) :: status
        type(plasma_t) :: p
        real(dp) :: q(TD_NX), cur(TD_NX), dcur(TD_NX), ratio(TD_NX), j(TD_NX)
        real(dp) :: w, betat, pavg, ipma, pohm, paux, prad, vloop, dw, g, xs, rd
        integer :: m, rc, k
        out = 0
        prof = 0
        call plasma_from_par(par, p, status)
        if (status /= 0) return
        call ballooning_table()
        call calibrate(p)
        call q_profile(p, q)
        ipma = p%ip/1e6_dp
        w = stored_energy(p)
        pavg = w/(1.5_dp*plasma_volume(p))
        betat = 2*mu0*pavg/p%b0**2
        call powers(p, pohm, paux, prad, vloop)
        out(1) = q(1)
        out(2) = q(TD_NX)
        out(3) = internal_inductance(p)
        out(4) = 100*betat
        out(5) = 100*betat*p%a*p%b0/ipma
        out(6) = 2*mu0*pavg/(mu0*p%ip/(2*pi*p%a*sqrt((1 + p%kappa**2)/2)))**2
        out(7) = ipma/(pi*p%a**2)
        out(8) = p%nbar/1e20_dp/out(7)
        out(9) = w/1e6_dp
        out(10) = tau98(p, p%ip, (pohm + paux)/1e6_dp)
        out(11) = (pohm + paux)/1e6_dp
        out(12) = pohm/1e6_dp
        out(13) = paux/1e6_dp
        out(14) = prad/1e6_dp
        out(15) = vloop
        out(16) = alfven_time(p)
        do m = 1, 4
            call cyl_kink_energy(TD_NX, q, m, 1, p%wall, dw, g, rc)
            out(16+m) = min(dw, 1e3_dp)
            out(20+m) = real(rc, dp)
            out(24+m) = g/out(16)
        end do
        call cyl_current_from_q(TD_NX, q, cur, dcur)
        call cyl_tearing_delta_prime(TD_NX, q, dcur, 2, 1, p%wall, 2e-3_dp, xs, rd, rc)
        out(29) = xs
        out(30) = merge(rd, 0.0_dp, rc == CYL_OK)
        call cyl_tearing_delta_prime(TD_NX, q, dcur, 3, 2, p%wall, 2e-3_dp, xs, rd, rc)
        out(31) = xs
        out(32) = merge(rd, 0.0_dp, rc == CYL_OK)
        call ballooning_ratio(p, q, ratio)
        k = maxloc(ratio, 1)
        out(33) = ratio(k)
        out(34) = p%x(k)
        if (ratio(k) > 0) then
            out(35) = out(4)/ratio(k)
            out(36) = out(5)/ratio(k)
        end if
        out(37) = 2.8_dp
        out(38) = 4*out(3)
        out(39) = mu0*p%ae**2/spitzer(p, p%t(1))
        out(40) = p%chi0
        call current_density(p, j)
        prof(:, 1) = p%t
        prof(:, 2) = p%n/1e20_dp
        prof(:, 3) = q
        prof(:, 4) = j/1e6_dp
        prof(:, 5) = ratio
    end subroutine

    ! ---------------------------------------------------------------- dynamics

    !> Thomas algorithm for a tridiagonal system (lower l, diagonal d, upper u).
    pure subroutine tridiag(nn, l, d, u, r, sol)
        integer, intent(in) :: nn
        real(dp), intent(in) :: l(nn), d(nn), u(nn), r(nn)
        real(dp), intent(out) :: sol(nn)
        real(dp) :: c(nn), g(nn), den
        integer :: i
        c(1) = u(1)/d(1)
        g(1) = r(1)/d(1)
        do i = 2, nn
            den = d(i) - l(i)*c(i-1)
            c(i) = u(i)/den
            g(i) = (r(i) - l(i)*g(i-1))/den
        end do
        sol(nn) = g(nn)
        do i = nn-1, 1, -1
            sol(i) = g(i) - c(i)*sol(i+1)
        end do
    end subroutine

    !> One implicit heat step; chi(i) is the diffusivity on the nodes.
    subroutine heat_step(p, chi, dt)
        type(plasma_t), intent(inout) :: p
        real(dp), intent(in) :: chi(TD_NX), dt
        real(dp) :: l(TD_NX), d(TD_NX), u(TD_NX), r(TD_NX), j(TD_NX), sol(TD_NX)
        real(dp) :: h, v, cr, cl, xf, kf, src, dsrc, rl, drl, eta, norm, tmin
        integer :: i
        h = 1.0_dp/real(TD_NX-1, dp)
        call current_density(p, j)
        norm = 0
        do i = 1, TD_NX
            norm = norm + p%heat(i)*cell_volume(p, i)
        end do
        tmin = 0.001_dp
        l = 0; u = 0
        do i = 1, TD_NX-1
            v = cell_volume(p, i)/(4*pi*pi*p%r0*p%ae**2)
            xf = (real(i, dp) - 0.5_dp)*h
            kf = (p%n(i)*chi(i) + p%n(i+1)*chi(i+1))
            cr = xf*kf/(p%ae**2*h)
            cl = 0
            if (i > 1) then
                xf = (real(i, dp) - 1.5_dp)*h
                kf = (p%n(i)*chi(i) + p%n(i-1)*chi(i-1))
                cl = xf*kf/(p%ae**2*h)
            end if
            eta = spitzer(p, p%t(i))
            call radiation(p, p%t(i), rl, drl)
            ! Heating explicit; radiation as the Picard sink -(n^2 L/T) T_new,
            ! which keeps T positive through the radiative collapse.
            src = (eta*j(i)**2 + p%paux*p%heat(i)/norm)/kev
            dsrc = p%n(i)**2*rl/(p%t(i)*kev)
            d(i) = 3*p%n(i)*v/dt + cr + cl + v*dsrc
            u(i) = -cr
            l(i) = -cl
            r(i) = 3*p%n(i)*v/dt*p%t(i) + v*src
        end do
        d(TD_NX) = 1; l(TD_NX) = 0; u(TD_NX) = 0; r(TD_NX) = p%tedge
        call tridiag(TD_NX, l, d, u, r, sol)
        p%t = max(sol, tmin)
    end subroutine

    !> One implicit current-diffusion step. With feedback the edge holds
    !> ip_target; without, the edge obeys L_ext dI/dt = -2 pi R0 eta j.
    subroutine current_step(p, dt, feedback, ip_target)
        type(plasma_t), intent(inout) :: p
        real(dp), intent(in) :: dt, ip_target
        logical, intent(in) :: feedback
        real(dp) :: l(TD_NX), d(TD_NX), u(TD_NX), r(TD_NX), sol(TD_NX)
        real(dp) :: h, el, er, xf, c
        integer :: i
        h = 1.0_dp/real(TD_NX-1, dp)
        l = 0; u = 0
        d(1) = 1; r(1) = 0
        do i = 2, TD_NX-1
            xf = (real(i, dp) - 0.5_dp)*h
            er = (spitzer(p, p%t(i)) + spitzer(p, p%t(i+1)))/(2*xf)
            xf = (real(i, dp) - 1.5_dp)*h
            el = (spitzer(p, p%t(i)) + spitzer(p, p%t(i-1)))/(2*xf)
            c = dt*p%x(i)/(mu0*p%ae**2*h*h)
            d(i) = 1 + c*(er + el)
            u(i) = -c*er
            l(i) = -c*el
            r(i) = p%cur(i)
        end do
        if (feedback) then
            d(TD_NX) = 1; r(TD_NX) = ip_target
        else
            c = dt*p%r0*spitzer(p, p%t(TD_NX))/(p%lext*p%ae**2*h)
            d(TD_NX) = 1 + c
            l(TD_NX) = -c
            r(TD_NX) = p%cur(TD_NX)
        end if
        call tridiag(TD_NX, l, d, u, r, sol)
        p%cur = sol
        p%cur(1) = 0
    end subroutine

    !> Thermal-quench current flattening: J moves a third of the way to uniform
    !> (surrogate for the observed 5-10 % spikes), and the
    !> total current follows poloidal-flux conservation (L_ext + L_i) I = const.
    subroutine flatten_current(p)
        type(plasma_t), intent(inout) :: p
        real(dp) :: li0, li1, ip
        integer :: i
        li0 = internal_inductance(p)
        ip = p%cur(TD_NX)
        do i = 1, TD_NX
            p%cur(i) = (2*p%cur(i) + ip*p%x(i)**2)/3
        end do
        li1 = internal_inductance(p)
        p%cur = p%cur*(p%lext + mu0*p%r0*li0/2)/(p%lext + mu0*p%r0*li1/2)
    end subroutine

    !> Sawtooth crash (Kadomtsev-like surrogate): T and q flattened inside 1.4 x_1,
    !> blended into the old profiles over a layer of 0.1 in x.
    subroutine sawtooth(p, q)
        type(plasma_t), intent(inout) :: p
        real(dp), intent(in) :: q(TD_NX)
        real(dp) :: w, tm, vtot, v, xmix, b
        integer :: i, imix
        imix = 1
        do i = 2, TD_NX
            if (q(i) < 1) imix = i
        end do
        xmix = min(0.9_dp, 1.4_dp*p%x(imix))
        imix = min(TD_NX-2, nint(xmix*real(TD_NX-1, dp)) + 1)
        w = 0; vtot = 0
        do i = 1, imix
            v = cell_volume(p, i)
            w = w + p%n(i)*p%t(i)*v
            vtot = vtot + p%n(i)*v
        end do
        tm = w/vtot
        do i = 1, TD_NX-1
            b = min(max((p%x(i) - xmix)/0.1_dp, 0.0_dp), 1.0_dp)
            if (b >= 1) exit
            p%t(i) = (1 - b)*tm + b*p%t(i)
            p%cur(i) = (1 - b)*min(p%cq*p%x(i)**2/1.01_dp, p%cur(i)) + b*p%cur(i)
        end do
    end subroutine

    pure real(dp) function ramp_value(t, t0, dur, a, b) result(v)
        real(dp), intent(in) :: t, t0, dur, a, b
        real(dp) :: s
        s = min(max((t - t0)/max(dur, 1e-9_dp), 0.0_dp), 1.0_dp)
        v = a + (b - a)*s
    end function

    !> Runs a discharge: after `settle` seconds the line density, auxiliary
    !> power and plasma current ramp linearly over `duration` to their final
    !> values; the run ends at t_end or when the current is gone.
    !> ramp = [duration, nbar_final(1e20), Paux_final(MW, <0 keeps), Ip_final(MA),
    !>         t_end, settle, -, -]. series(:,k) holds TD_NSER quantities per
    !> record; profiles(:,:,k) T, n, q, j, alpha/alpha_c.
    subroutine td_run(par, ramp, maxrec, nrec, series, profiles, status)
        real(dp), intent(in) :: par(TD_NPAR), ramp(TD_NRAMP)
        integer, intent(in) :: maxrec
        integer, intent(out) :: nrec, status
        real(dp), intent(out) :: series(TD_NSER, maxrec), &
            profiles(TD_NX, TD_NPROF, maxrec)
        type(plasma_t) :: p
        real(dp) :: q(TD_NX), cur(TD_NX), dcur(TD_NX), chi(TD_NX), ratio(TD_NX)
        real(dp) :: j(TD_NX), shear(TD_NX), drive(TD_NX)
        logical :: crossed
        real(dp) :: t, dt, dt0, tend, settle, dur, n_i, n_f, pa_i, pa_f, ip_i, ip_f
        real(dp) :: ip_now, w21, w32, x21, x32, rd21, rd32, wseed, chi_mult, taua
        real(dp) :: xi, gam, z, gamv, gap, xnow, ip0, t0prev, ipprev, wprev, rel
        real(dp) :: w, pohm, paux, prad, vloop, tau, wsat, rate, dw, g, gmax, br
        real(dp) :: t_last_saw, dtcap, wtq
        integer :: phase, cause, m, rc, rc21, rc32, step, mgrow, k
        logical :: feedback
        nrec = 0
        series = 0
        profiles = 0
        call plasma_from_par(par, p, status)
        if (status /= 0) return
        if (maxrec < 2) then
            status = 1
            return
        end if
        call ballooning_table()
        call calibrate(p)
        dur = max(ramp(1), 1e-3_dp)
        n_i = p%nbar; n_f = ramp(2)*1e20_dp
        if (n_f <= 0) n_f = n_i
        pa_i = p%paux; pa_f = ramp(3)*1e6_dp
        if (ramp(3) < 0) pa_f = pa_i
        ip_i = p%ip; ip_f = ramp(4)*1e6_dp
        if (ip_f <= 0) ip_f = ip_i
        tend = max(ramp(5), 1e-3_dp)
        settle = max(ramp(6), 0.0_dp)
        dt0 = min(5e-3_dp, max(2e-5_dp, (settle + dur)/1500))
        dt = dt0/20
        t = 0
        phase = TD_FLATTOP
        cause = TD_CAUSE_NONE
        feedback = .true.
        wseed = 2e-3_dp*p%ae
        w21 = wseed; w32 = wseed
        x21 = 0; x32 = 0; rd21 = 0; rd32 = 0
        chi_mult = 1
        xi = 0; gam = 0; mgrow = 0
        z = 0
        gamv = (p%kappa - 1)/(max(0.1_dp, 2.3_dp - p%kappa)*p%tauw)
        gap = (p%wall - 1)*p%a*p%kappa
        xnow = 1
        ip0 = p%ip
        t_last_saw = -1
        wtq = 1
        do step = 1, 10*maxrec
            ip_now = ramp_value(t, settle, dur, ip_i, ip_f)
            if (feedback) then
                call set_density(p, ramp_value(t, settle, dur, n_i, n_f))
                p%paux = ramp_value(t, settle, dur, pa_i, pa_f)
            else
                p%paux = 0
            end if
            call q_profile(p, q)
            ! stability of the present profiles
            call cyl_current_from_q(TD_NX, q, cur, dcur)
            call cyl_tearing_delta_prime(TD_NX, q, dcur, 2, 1, p%wall, 2e-3_dp, &
                x21, rd21, rc21)
            call cyl_tearing_delta_prime(TD_NX, q, dcur, 3, 2, p%wall, 2e-3_dp, &
                x32, rd32, rc32)
            call ballooning_ratio(p, q, ratio, shear, drive)
            ! Beyond the first boundary the point may lie past the unstable
            ! tongue (second stability at low shear): confirm by the direct
            ! Newcomb test at the local (s, alpha).
            br = 0
            do k = 2, TD_NX-1
                ! The local high-n model is pessimistic at low shear (finite
                ! Larmor radius, diamagnetic and sawtooth effects): s >= 0.4 only.
                if (p%x(k) > 0.7_dp .or. shear(k) < 0.4_dp) cycle
                g = (ratio(k-1) + ratio(k) + ratio(k+1))/3
                if (g > 1) then
                    call ballooning_first_zero(max(shear(k), 1e-3_dp), drive(k), &
                        12*pi, 768, crossed, rc)
                    if (rc /= BALLOON_OK .or. .not. crossed) g = min(g, 0.99_dp)
                end if
                br = max(br, g)
            end do
            taua = alfven_time(p)
            ! Radiative collapse: the cold radiating front has reached x = 0.8.
            if (phase <= TD_TEARING .and. interp(p%t, 0.8_dp) < 0.3_dp) &
                call start_quench(TD_CAUSE_RADIATIVE)
            if (phase <= TD_TEARING) then
                gmax = 0
                ! Only m = 2 disrupts: toroidal coupling closes the cylindrical
                ! m >= 3 windows, so the low-n limit is q95 = 2. Marginal kinks
                ! (|dW| below a tenth of the m = 2 window depth) saturate.
                do m = 2, 2
                    call cyl_kink_energy(TD_NX, q, m, 1, p%wall, dw, g, rc)
                    if (rc == CYL_OK .and. dw < -2e-3_dp .and. g > gmax) then
                        gmax = g
                        mgrow = m
                    end if
                end do
                if (gmax > 0) then
                    phase = TD_IDEAL
                    cause = TD_CAUSE_KINK
                    gam = gmax/taua
                    xi = 1e-4_dp*p%a
                else if (br > 1) then
                    phase = TD_IDEAL
                    cause = TD_CAUSE_BALLOONING
                    mgrow = 10
                    gam = sqrt(br - 1)/(q((TD_NX+1)/2)*taua)
                    xi = 1e-4_dp*p%a
                end if
            end if
            ! islands (Rutherford with a saturation surrogate)
            if (phase <= TD_IDEAL) then
                w21 = island(w21, rc21, x21, rd21)
                w32 = island(w32, rc32, x32, rd32)
                if (phase == TD_FLATTOP .and. max(w21, w32) > 5*wseed) &
                    phase = TD_TEARING
            end if
            if (phase <= TD_TEARING .and. rc21 == CYL_OK) then
                if ((rc32 == CYL_OK .and. (w21 + w32)/(2*p%ae) >= abs(x21 - x32)) &
                    .or. w21/(2*p%ae) >= 1 - x21) call start_quench(TD_CAUSE_ISLANDS)
            end if
            if (phase == TD_IDEAL) then
                xi = min(xi*exp(min(gam*dt, 50.0_dp)), p%a)
                if (xi > 0.1_dp*p%a) call start_quench(cause)
            end if
            ! sawteeth
            if (phase <= TD_TEARING .and. q(1) < 0.97_dp .and. &
                t - t_last_saw > 0.05_dp*max(p%chi_ref_tau, 1e-3_dp)) then
                call sawtooth(p, q)
                t_last_saw = t
            end if
            ! transport
            tau = tau98(p, p%cur(TD_NX), (p%paux + ohmic_power(p))/1e6_dp)
            do k = 1, TD_NX
                chi(k) = p%chi0*p%chi_ref_tau/tau*(1 + 3*p%x(k)**2)*chi_mult
                if (phase <= TD_IDEAL) then
                    if (rc21 == CYL_OK .and. abs(p%x(k) - x21) < w21/(2*p%ae)) &
                        chi(k) = 30*chi(k)
                    if (rc32 == CYL_OK .and. abs(p%x(k) - x32) < w32/(2*p%ae)) &
                        chi(k) = 30*chi(k)
                end if
            end do
            t0prev = p%t(1); ipprev = p%cur(TD_NX); wprev = stored_energy(p)
            call heat_step(p, chi, dt)
            call current_step(p, dt, feedback, ip_now)
            if (phase == TD_THERMAL_QUENCH .and. stored_energy(p) < 0.05_dp*wtq) then
                phase = TD_CURRENT_QUENCH
                chi_mult = 300
            end if
            ! vertical displacement once position control is lost
            if (phase >= TD_THERMAL_QUENCH) z = z*exp(gamv*dt)
            xnow = 1
            if (abs(z) > gap) xnow = max(0.0_dp, 1 - (abs(z) - gap)/(p%a*p%kappa))
            t = t + dt
            call record()
            if (nrec >= maxrec) exit
            if (phase >= TD_THERMAL_QUENCH .and. (xnow < 0.25_dp .or. &
                p%cur(TD_NX) < 0.01_dp*ip0)) then
                phase = TD_TERMINATED
                series(21, nrec) = phase
                exit
            end if
            if (t >= tend) then
                if (phase < TD_THERMAL_QUENCH) series(21, nrec) = TD_SURVIVED
                exit
            end if
            ! time step from the relative change of T0, I_p and W
            w = stored_energy(p)
            rel = max(abs(p%t(1) - t0prev)/t0prev, abs(p%cur(TD_NX) - ipprev) &
                /max(ipprev, 1.0_dp), abs(w - wprev)/max(wprev, 1.0_dp), 1e-12_dp)
            dt = dt*min(1.5_dp, max(0.3_dp, 0.03_dp/rel))
            dtcap = dt0
            if (phase == TD_TEARING) dtcap = min(dtcap, island_time())
            if (phase == TD_IDEAL) dtcap = min(dtcap, 0.1_dp/gam)
            if (phase >= TD_THERMAL_QUENCH) dtcap = min(dtcap, 0.05_dp/gamv)
            dt = min(max(dt, 1e-9_dp), dtcap)
        end do
        status = 0
    contains
        real(dp) function island(wi, rci, xs, rsd) result(wn)
            real(dp), intent(in) :: wi, xs, rsd
            integer, intent(in) :: rci
            real(dp) :: rs
            wn = wseed
            if (rci /= CYL_OK) return
            rs = xs*p%ae
            if (rsd > 0) then
                wsat = rs*rsd/(rsd + 12)
                rate = 1.22_dp*spitzer(p, interp(p%t, xs))/mu0*(rsd/rs)*(1 - wi/wsat)
            else
                rate = 1.22_dp*spitzer(p, interp(p%t, xs))/mu0*rsd/rs
            end if
            wn = max(wseed, wi + dt*rate)
        end function
        real(dp) function island_time() result(tc)
            tc = 1
            if (rate > 0) tc = 0.05_dp*max(w21, w32)/rate
            tc = max(tc, 1e-7_dp)
        end function
        subroutine start_quench(why)
            integer, intent(in) :: why
            if (phase >= TD_THERMAL_QUENCH) return
            cause = why
            phase = TD_THERMAL_QUENCH
            feedback = .false.
            chi_mult = 3000
            ! Surrogate for the impurity influx from the wall during the quench;
            ! the edge falls to wall temperature.
            p%fimp = max(p%fimp, 0.03_dp)
            p%tedge = min(p%tedge, 0.003_dp)
            wtq = stored_energy(p)
            call flatten_current(p)
            z = 0.02_dp*p%a
            dt = min(dt, 1e-6_dp)
        end subroutine
        subroutine record()
            real(dp) :: ipx, pavg, betat
            nrec = nrec + 1
            call q_profile(p, q)
            call powers(p, pohm, paux, prad, vloop)
            call current_density(p, j)
            w = stored_energy(p)
            ipx = interp(p%cur, xnow)
            pavg = w/(1.5_dp*plasma_volume(p))
            betat = 2*mu0*pavg/p%b0**2
            series(1, nrec) = t
            series(2, nrec) = ipx/1e6_dp
            series(3, nrec) = p%t(1)
            series(4, nrec) = w/1e6_dp
            series(5, nrec) = p%nbar/1e20_dp
            series(6, nrec) = p%nbar/1e20_dp/(ipx/1e6_dp/(pi*p%a**2) + 1e-30_dp)
            series(7, nrec) = 100*betat*p%a*p%b0/max(ipx/1e6_dp, 1e-6_dp)
            series(8, nrec) = q(1)
            series(9, nrec) = interp(q, xnow)
            series(10, nrec) = internal_inductance(p)
            series(11, nrec) = z
            series(12, nrec) = xnow*p%a
            series(13, nrec) = merge(w21, 0.0_dp, rc21 == CYL_OK)
            series(14, nrec) = merge(w32, 0.0_dp, rc32 == CYL_OK)
            series(15, nrec) = merge(x21, -1.0_dp, rc21 == CYL_OK)
            series(16, nrec) = merge(x32, -1.0_dp, rc32 == CYL_OK)
            series(17, nrec) = prad/1e6_dp
            series(18, nrec) = pohm/1e6_dp
            series(19, nrec) = p%paux/1e6_dp
            series(20, nrec) = vloop
            series(21, nrec) = phase
            series(22, nrec) = cause
            series(23, nrec) = maxval(ratio)
            series(24, nrec) = xi
            series(25, nrec) = mgrow
            series(26, nrec) = merge(rd21, 0.0_dp, rc21 == CYL_OK)
            series(27, nrec) = merge(rd32, 0.0_dp, rc32 == CYL_OK)
            series(28, nrec) = dt
            profiles(:, 1, nrec) = p%t
            profiles(:, 2, nrec) = p%n/1e20_dp
            profiles(:, 3, nrec) = q
            profiles(:, 4, nrec) = j/1e6_dp
            profiles(:, 5, nrec) = ratio
        end subroutine
    end subroutine

    real(dp) function ohmic_power(p) result(pw)
        type(plasma_t), intent(in) :: p
        real(dp) :: pohm, paux, prad, vloop
        call powers(p, pohm, paux, prad, vloop)
        pw = pohm
    end function

    pure real(dp) function interp(f, x) result(v)
        real(dp), intent(in) :: f(TD_NX), x
        real(dp) :: s, w
        integer :: i
        s = min(max(x, 0.0_dp), 1.0_dp)*real(TD_NX-1, dp)
        i = min(int(s) + 1, TD_NX-1)
        w = s - real(i-1, dp)
        v = (1-w)*f(i) + w*f(i+1)
    end function

    ! --------------------------------------------------------------- C/WASM ABI

    subroutine milca_td_sizes(nx, npar, nout, nser, nprof, nramp) &
            bind(C, name='milca_td_sizes')
        integer(c_int), intent(out) :: nx, npar, nout, nser, nprof, nramp
        nx = TD_NX; npar = TD_NPAR; nout = TD_NOUT
        nser = TD_NSER; nprof = TD_NPROF; nramp = TD_NRAMP
    end subroutine

    subroutine milca_td_limits(par, out, prof, status) bind(C, name='milca_td_limits')
        real(c_double), intent(in) :: par(TD_NPAR)
        real(c_double), intent(out) :: out(TD_NOUT), prof(TD_NX, TD_NPROF)
        integer(c_int), intent(out) :: status
        integer :: rc
        call td_limits(par, out, prof, rc)
        status = rc
    end subroutine

    subroutine milca_td_run(par, ramp, maxrec, nrec, series, profiles, status) &
            bind(C, name='milca_td_run')
        real(c_double), intent(in) :: par(TD_NPAR), ramp(TD_NRAMP)
        integer(c_int), value, intent(in) :: maxrec
        integer(c_int), intent(out) :: nrec, status
        real(c_double), intent(out) :: series(TD_NSER, maxrec), &
            profiles(TD_NX, TD_NPROF, maxrec)
        integer :: n, rc
        call td_run(par, ramp, maxrec, n, series, profiles, rc)
        nrec = n
        status = rc
    end subroutine

end module milca_tokamak_discharge
