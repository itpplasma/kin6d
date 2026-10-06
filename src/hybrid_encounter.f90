!> QG6: local quantum treatment of binary close encounters embedded in a
!> classical test-particle plasma.
!>
!> Physical model (declared): test charges in a fixed Maxwellian background of
!> field charges; binary encounters with the bare Coulomb pair potential and a
!> sharp impact-parameter cut at b_max = lambda_D (Spitzer screening model);
!> encounters are instantaneous on the plasma clock (binary-collision
!> approximation). Relative motion uses the shared nbody KDK step with
!> encounter-resolved step control (reduced mass, fixed source; no softening).
!>
!> Model switch (interface, not a flag):
!> - indicator: lifted asymptotic eta = g/(hbar v_inf) < eta_switch and lifted
!>   impact parameter b < b_s, detected when the classical integrator first
!>   enters the approach sphere r_det = 2 (b90 + sqrt(b90^2 + b_s^2)), which
!>   every trajectory with b < b_s enters;
!> - lifting: the current classical state is mapped to its exact incoming
!>   Kepler asymptote (v_inf, unit u_in, angular momentum h) through the
!>   conserved energy, angular momentum and Runge-Lenz vector;
!> - quantum kernel: partial-wave ring l with l(l+1) <= (k b)^2 < (l+1)(l+2)
!>   (ring areas 2 pi (l+1)/k^2 match the weights of the exact l-local
!>   momentum-transfer sum); outgoing asymptotic direction rotated by
!>   chi_l = 2 (sigma_{l+1} - sigma_l) = 2 atan(eta/(l+1)) in the orbit plane;
!> - projection: the outgoing classical state on the sphere r_det is the point
!>   of the repulsive hyperbola with asymptote u_out and impact parameter
!>   b' = b90 cot(chi/2) = (l+1)/k; the classical integrator continues from it.
!>   Energy (|v| at r_det) and pair momentum are conserved by construction; the
!>   angular momentum changes by |b' - b| mu v_inf <= hbar (ring width).
!> The quantum event replaces the whole classical hyperbola (both legs are
!> represented by the lifted/projected asymptotes), so no near force is
!> counted twice; the far zone b > b_s stays classical. The naive diagnostic
!> (rotate v at r_det) double counts the outer-leg deflection and is reported.
module kin6d_hybrid_encounter
    use, intrinsic :: iso_fortran_env, only: dp => real64, i64 => int64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortnum_status, only: fortnum_status_t, FORTNUM_OK
    use fortnum_rng, only: rng_t, rng_seed, rng_uniform, rng_normal
    use fortnum_quadrature, only: gauss_legendre_ab
    use kin6d_species, only: units_t, species_t, pair_scales_t, pair_scales
    use kin6d_nbody, only: nbody_system_t, nbody_isolated, nbody_add_source, &
                           nbody_step, nbody_encounter_step
    use kin6d_quantum_collision_kernel, only: quantum_cut_log, classical_cut_log
    use kin6d_quantum_collision_leaves
    implicit none
    private
    public :: hybrid_config_t, hybrid_model_t, encounter_record_t
    public :: hybrid_init, hybrid_encounter, hybrid_zone, hybrid_transfer_log
    public :: lift_asymptote, orbit_point, hybrid_rate_samples, kernel_slowing_rate
    public :: hybrid_beam_step
    public :: MODEL_CLASSICAL, MODEL_QUANTUM, MODEL_HYBRID

    integer, parameter :: MODEL_CLASSICAL = 0, MODEL_QUANTUM = 1, MODEL_HYBRID = 2
    real(dp), parameter :: pi = 3.14159265358979323846_dp

    type :: hybrid_config_t
        real(dp) :: eta_switch = 1        ! promote when lifted eta < eta_switch
        real(dp) :: zone = 4              ! near rings L = ceiling(zone max(eta, 1))
        real(dp) :: step_eta = 0.02_dp    ! encounter-resolved KDK step factor
        real(dp) :: start_factor = 3      ! start/exit sphere R0 / max(b, r_det)
        integer :: max_steps = 400000
        logical :: promote = .true.       ! .false.: all-classical reference
        logical :: naive = .false.        ! diagnostic: rotate v at r_det only
    end type

    type :: hybrid_model_t
        type(hybrid_config_t) :: cfg
        real(dp) :: g = 0, mu = 0, hbar = 1, mt = 0, mf = 0, bmax = 0
        type(nbody_system_t) :: relative
    end type

    type :: encounter_record_t
        real(dp) :: b = 0, eta = 0, theta = 0, energy_defect = 0
        real(dp) :: angular_change = 0     ! |h' - h| mu / hbar at the event
        real(dp) :: chi = 0                ! quantum ring deflection (promoted)
        integer :: steps = 0, ring = -1
        logical :: promoted = .false.
    end type
contains
    !> Pair of a test species and a field species with impact cut bmax.
    subroutine hybrid_init(model, units, test, field, bmax, cfg, ok)
        type(hybrid_model_t), intent(out) :: model
        type(units_t), intent(in) :: units
        type(species_t), intent(in) :: test, field
        real(dp), intent(in) :: bmax
        type(hybrid_config_t), intent(in) :: cfg
        logical, intent(out) :: ok
        type(pair_scales_t) :: scales
        ok = .false.
        scales = pair_scales(units, test, field)
        if (.not. (scales%coupling > 0 .and. bmax > 0 .and. cfg%step_eta > 0 .and. &
                   cfg%start_factor > 1 .and. cfg%zone > 0)) return
        model%cfg = cfg; model%g = scales%coupling; model%mu = scales%reduced_mass
        model%hbar = units%hbar; model%mt = test%mass; model%mf = field%mass
        model%bmax = bmax
        call nbody_isolated(model%relative, units, &
                            [species_t('relative', test%charge, model%mu)], ok)
        if (.not. ok) return
        call nbody_add_source(model%relative, [0.0_dp, 0.0_dp, 0.0_dp], field%charge, ok)
    end subroutine hybrid_init

    !> Near zone at asymptotic relative speed v: rings l < nring, radius bs.
    !> bs = sqrt(L(L+1))/k with L = ceiling(zone max(eta,1)), capped at bmax
    !> (then the ring kernel covers the whole disk, last ring partial).
    pure subroutine hybrid_zone(model, v, eta, k, b90, nring, bs)
        type(hybrid_model_t), intent(in) :: model
        real(dp), intent(in) :: v
        real(dp), intent(out) :: eta, k, b90, bs
        integer, intent(out) :: nring
        real(dp) :: lambdabar
        call qc_pair_scales(model%g, model%mu, model%hbar, v, eta, k, lambdabar, b90)
        nring = ceiling(model%cfg%zone*max(eta, 1.0_dp))
        bs = sqrt(real(nring, dp)*(nring + 1))/k
        if (bs >= model%bmax) then
            bs = model%bmax
            nring = floor(k*model%bmax) + 1
        end if
    end subroutine hybrid_zone

    pure function cross(a, b) result(c)
        real(dp), intent(in) :: a(3), b(3)
        real(dp) :: c(3)
        c = [a(2)*b(3) - a(3)*b(2), a(3)*b(1) - a(1)*b(3), a(1)*b(2) - a(2)*b(1)]
    end function cross

    !> Exact repulsive Kepler asymptotes of the relative state (x, v):
    !> v_inf, angular momentum h = x cross v, Runge-Lenz A = v cross h + (g/mu) x/|x|,
    !> incoming/outgoing asymptotic unit velocities.
    pure subroutine lift_asymptote(model, x, v, vinf, h, a, uin, uout)
        type(hybrid_model_t), intent(in) :: model
        real(dp), intent(in) :: x(3), v(3)
        real(dp), intent(out) :: vinf, h(3), a(3), uin(3), uout(3)
        real(dp) :: kg, hm, den, ha(3)
        kg = model%g/model%mu
        call qc_asymptotic_speed(sum(v*v), norm2(x), model%g, model%mu, vinf)
        h = cross(x, v); hm = norm2(h)
        a = cross(v, h) + kg*x/norm2(x)
        den = kg**2 + (vinf*hm)**2
        ha = cross(h/hm, a)
        uout = (kg*a + vinf*hm*ha)/den
        uin = (-kg*a + vinf*hm*ha)/den
    end subroutine lift_asymptote

    !> State at radius r on the hyperbola (A, h, v_inf); branch -1 incoming,
    !> +1 outgoing. The in-plane sign is the one reproducing A.
    pure subroutine orbit_point(model, a, h, vinf, r, branch, x, v, ok)
        type(hybrid_model_t), intent(in) :: model
        real(dp), intent(in) :: a(3), h(3), vinf, r
        integer, intent(in) :: branch
        real(dp), intent(out) :: x(3), v(3)
        logical, intent(out) :: ok
        real(dp) :: hm, am, c, s, vr, kg, xr(3), vel(3), err, best, ah(3), hh(3)
        integer :: sgn
        kg = model%g/model%mu
        hm = norm2(h); am = norm2(a); hh = h/hm; ah = a/am
        call qc_orbit_cosine(hm**2, r, model%g, model%mu, am, c)
        ok = c <= 1
        if (.not. ok) return
        s = sqrt(max(0.0_dp, 1 - c*c))
        call qc_radial_speed(vinf, r, model%g, model%mu, hm, vr)
        vr = branch*vr
        best = huge(1.0_dp)
        do sgn = -1, 1, 2
            xr = c*ah + sgn*s*cross(hh, ah)
            vel = vr*xr + (hm/r)*cross(hh, xr)
            err = norm2(cross(vel, h) + kg*xr - a)
            if (err < best) then
                best = err; x = r*xr; v = vel
            end if
        end do
        ok = ieee_is_finite(best)
    end subroutine orbit_point

    !> One binary encounter with asymptotic relative velocity gvec and impact
    !> vector bvec (perpendicular to gvec). Returns the outgoing asymptotic
    !> relative velocity gout (|gout| = |gvec| by construction; the integrator's
    !> energy error is recorded, not propagated) and the encounter record.
    subroutine hybrid_encounter(model, gvec, bvec, gout, rec, ok)
        type(hybrid_model_t), intent(inout) :: model
        real(dp), intent(in) :: gvec(3), bvec(3)
        real(dp), intent(out) :: gout(3)
        type(encounter_record_t), intent(out) :: rec
        logical, intent(out) :: ok
        real(dp) :: vinf, u(3), h(3), a(3), kg, eta, k, b90, bs, rdet, r0, y(6), dt
        real(dp) :: vl, hl(3), al(3), uin(3), uout(3), bl, etal, kl, b90l, bsl, chi
        real(dp) :: bhat(3), bp, hp(3), e0, e1, ul(3)
        integer :: nring, nr, l
        type(fortnum_status_t) :: status
        logical :: inside, detected
        ok = .false.
        kg = model%g/model%mu
        vinf = norm2(gvec); u = gvec/vinf
        rec%b = norm2(bvec)
        h = cross(bvec, gvec)
        a = vinf*norm2(h)*cross(u, h/norm2(h)) - kg*u
        call hybrid_zone(model, vinf, eta, k, b90, nring, bs)
        rec%eta = eta
        rdet = 2*(b90 + sqrt(b90**2 + bs**2))
        r0 = model%cfg%start_factor*max(rec%b, rdet)
        call orbit_point(model, a, h, vinf, r0, -1, y(1:3), y(4:6), ok)
        if (.not. ok) return
        e0 = sum(y(4:6)**2)/2 + kg/norm2(y(1:3))
        detected = .false.; inside = .false.
        do while (rec%steps < model%cfg%max_steps)
            call nbody_encounter_step(model%relative, y, model%cfg%step_eta, dt, status)
            if (status%code /= FORTNUM_OK) return
            call nbody_step(model%relative, y, dt, status)
            if (status%code /= FORTNUM_OK) return
            rec%steps = rec%steps + 1
            if (.not. detected .and. norm2(y(1:3)) < rdet) then
                detected = .true.
                ! Indicator from the lifted asymptotic state (approach detection).
                call lift_asymptote(model, y(1:3), y(4:6), vl, hl, al, uin, uout)
                bl = norm2(hl)/vl
                call hybrid_zone(model, vl, etal, kl, b90l, nr, bsl)
                if (model%cfg%promote .and. etal < model%cfg%eta_switch .and. bl < bsl) then
                    rec%promoted = .true.
                    l = ring_index(kl*bl)
                    rec%ring = l
                    call qc_ring_deflection(etal, real(l + 1, dp), chi)
                    rec%chi = chi
                    if (model%cfg%naive) then
                        y(4:6) = rotate(y(4:6), hl/norm2(hl), chi)
                    else
                        bhat = cross(uin, hl/norm2(hl))
                        ul = cos(chi)*uin + sin(chi)*bhat
                        call qc_rutherford_impact(chi, b90l, bp)
                        hp = bp*vl*hl/norm2(hl)
                        al = vl*norm2(hp)*cross(ul, hp/norm2(hp)) + kg*ul
                        call orbit_point(model, al, hp, vl, norm2(y(1:3)), 1, y(1:3), &
                                         y(4:6), ok)
                        if (.not. ok) return
                        rec%angular_change = abs(bp - bl)*vl*model%mu/model%hbar
                    end if
                end if
            end if
            if (norm2(y(1:3)) > r0 .and. dot_product(y(1:3), y(4:6)) > 0) then
                inside = .true.; exit
            end if
        end do
        if (.not. inside) return
        e1 = sum(y(4:6)**2)/2 + kg/norm2(y(1:3))
        rec%energy_defect = (e1 - e0)/e0
        call lift_asymptote(model, y(1:3), y(4:6), vl, hl, al, uin, uout)
        gout = vinf*uout/norm2(uout)
        rec%theta = atan2(norm2(cross(u, uout)), dot_product(u, uout))
        ok = all(ieee_is_finite(gout))
    contains
        pure integer function ring_index(kb) result(l)
            real(dp), intent(in) :: kb
            l = max(0, floor((sqrt(1 + 4*kb**2) - 1)/2))
            do while (real(l + 1, dp)*(l + 2) <= kb**2)
                l = l + 1
            end do
            do while (l > 0 .and. real(l, dp)*(l + 1) > kb**2)
                l = l - 1
            end do
        end function ring_index
        pure function rotate(w, axis, angle) result(wr)
            real(dp), intent(in) :: w(3), axis(3), angle
            real(dp) :: wr(3)
            ! Repulsive sense: rotation by -angle about h (same as the orbit).
            wr = w*cos(-angle) + cross(axis, w)*sin(-angle) + &
                 axis*dot_product(axis, w)*(1 - cos(-angle))
        end function rotate
    end subroutine hybrid_encounter

    !> Normalized momentum-transfer log sigma_tr/(4 pi b90^2) at speed v for the
    !> classical cut model, the all-quantum ring (partial-wave cut) model, or the
    !> hybrid model (quantum rings l < L, classical Rutherford b_s < b < b_max).
    !> dlnl = d lnL/d eta at fixed kappa and L.
    subroutine hybrid_transfer_log(model, v, which, lnl, dlnl)
        type(hybrid_model_t), intent(in) :: model
        real(dp), intent(in) :: v
        integer, intent(in) :: which
        real(dp), intent(out) :: lnl, dlnl
        real(dp) :: eta, k, b90, bs, kappa, t, dt, far, dfar
        integer :: nring, n
        call hybrid_zone(model, v, eta, k, b90, nring, bs)
        kappa = k*model%bmax
        select case (which)
        case (MODEL_QUANTUM)
            call quantum_cut_log(eta, kappa, lnl, dlnl)
        case (MODEL_HYBRID)
            if (.not. (model%cfg%promote .and. eta < model%cfg%eta_switch)) then
                call classical_cut_log(eta, kappa, lnl, dlnl)
            else if (bs >= model%bmax) then
                call quantum_cut_log(eta, kappa, lnl, dlnl)
            else
                lnl = 0; dlnl = 0
                do n = nring, 1, -1
                    call qc_ring_transfer(eta, real(n, dp), t)
                    call qc_ring_transfer_jvp(eta, real(n, dp), 1.0_dp, 0.0_dp, dt)
                    lnl = lnl + t; dlnl = dlnl + dt
                end do
                call qc_far_classical_log(eta, kappa, real(nring, dp)*(nring + 1), far)
                call qc_far_classical_log_jvp(eta, kappa, real(nring, dp)*(nring + 1), &
                                              1.0_dp, 0.0_dp, 0.0_dp, dfar)
                lnl = lnl + far; dlnl = dlnl + dfar
            end if
        case default
            call classical_cut_log(eta, kappa, lnl, dlnl)
        end select
    end subroutine hybrid_transfer_log

    !> Kernel quadrature of the test-particle slowing-down rate
    !> nu = -<dv_z/dt>/v_b = (n mu/(m_t v_b)) int W(g) J(a) sigma_tr(g) dg for a test
    !> particle at v_b z in a field Maxwellian at temperature T; sigma_tr from
    !> hybrid_transfer_log (which >= 0) or constant lnL = lnl_const (which < 0).
    function kernel_slowing_rate(model, vb, temperature, density, which, lnl_const) &
        result(nu)
        type(hybrid_model_t), intent(in) :: model
        real(dp), intent(in) :: vb, temperature, density
        integer, intent(in) :: which
        real(dp), intent(in), optional :: lnl_const
        real(dp) :: nu, vth, gmax, x(16), w(16), wt, j, aa, lnl, dlnl, eta, k, b90, bs
        integer :: ip, i, nring
        integer, parameter :: npanel = 96
        vth = sqrt(temperature/model%mf)
        gmax = vb + 14*vth
        nu = 0
        do ip = 1, npanel
            call gauss_legendre_ab(16, gmax*(ip - 1)/npanel, gmax*ip/npanel, x, w)
            do i = 1, 16
                call qc_maxwell_weight(x(i), vb, temperature, model%mf, wt)
                aa = model%mf*x(i)*vb/temperature
                if (aa < 0.05_dp) then
                    call qc_angle_moment_series(aa, j)
                else
                    call qc_angle_moment(aa, j)
                end if
                call hybrid_zone(model, x(i), eta, k, b90, nring, bs)
                if (which < 0) then
                    lnl = lnl_const
                else
                    call hybrid_transfer_log(model, x(i), which, lnl, dlnl)
                end if
                nu = nu + w(i)*wt*j*4*pi*b90**2*lnl
            end do
        end do
        nu = nu*density*model%mu/(model%mt*vb)
    end function kernel_slowing_rate

    !> Monte Carlo samples of the slowing-down rate nu = -<dv_z/dt>/v_b for a
    !> test particle at v_b z: field velocity from the Maxwellian, impact
    !> parameter by importance sampling (uniform disk below b_lo, log-uniform up
    !> to b_max, b_lo = 1e-3 min(b90, lambdabar)) with the exact uniform-disk
    !> weight, Rao-Blackwellized azimuth (E_phi[dg] = -(1 - cos theta) g).
    !> xs(i) are independent identically distributed estimates of nu; the same
    !> seed gives paired samples for promote on/off. Statistical law: CLT.
    subroutine hybrid_rate_samples(model, vb, temperature, density, seed, xs, nprom, &
                                   max_energy_defect, max_angular_change, ok)
        type(hybrid_model_t), intent(inout) :: model
        real(dp), intent(in) :: vb, temperature, density
        integer(i64), intent(in) :: seed
        real(dp), intent(out) :: xs(:)
        integer, intent(out) :: nprom
        real(dp), intent(out) :: max_energy_defect, max_angular_change
        logical, intent(out) :: ok
        type(rng_t) :: rng
        type(fortnum_status_t) :: status
        type(encounter_record_t) :: rec
        real(dp) :: vf(3), gvec(3), gout(3), bvec(3), vinf, eta, k, b90, bs, blo, u1, u2
        real(dp) :: b, p, w, vth
        real(dp), parameter :: eps = 0.05_dp
        integer :: i, nring, c
        ok = .false.; nprom = 0; max_energy_defect = 0; max_angular_change = 0
        call rng_seed(rng, seed, status)
        if (status%code /= FORTNUM_OK) return
        vth = sqrt(temperature/model%mf)
        do i = 1, size(xs)
            do c = 1, 3
                call rng_normal(rng, vf(c))
            end do
            gvec = [0.0_dp, 0.0_dp, vb] - vth*vf
            vinf = norm2(gvec)
            call hybrid_zone(model, vinf, eta, k, b90, nring, bs)
            blo = 1.0e-3_dp*min(b90, 1/k)
            call rng_uniform(rng, u1); call rng_uniform(rng, u2)
            if (u1 < eps) then
                b = blo*sqrt(u2)
                p = eps*2*b/blo**2
            else
                b = blo*(model%bmax/blo)**u2
                p = (1 - eps)/(b*log(model%bmax/blo))
            end if
            w = (2*b/model%bmax**2)/p
            call random_perpendicular(rng, gvec/vinf, bvec)
            call hybrid_encounter(model, gvec, b*bvec, gout, rec, ok)
            if (.not. ok) return
            if (rec%promoted) nprom = nprom + 1
            max_energy_defect = max(max_energy_defect, abs(rec%energy_defect))
            max_angular_change = max(max_angular_change, rec%angular_change)
            xs(i) = density*pi*model%bmax**2*vinf*w*(model%mu/model%mt)* &
                    2*sin(rec%theta/2)**2*gvec(3)/vb
        end do
        ok = .true.
    end subroutine hybrid_rate_samples

    !> Uniformly random unit vector perpendicular to unit u.
    subroutine random_perpendicular(rng, u, e)
        type(rng_t), intent(inout) :: rng
        real(dp), intent(in) :: u(3)
        real(dp), intent(out) :: e(3)
        real(dp) :: e1(3), e2(3), phi
        if (abs(u(1)) < 0.9_dp) then
            e1 = cross(u, [1.0_dp, 0.0_dp, 0.0_dp])
        else
            e1 = cross(u, [0.0_dp, 1.0_dp, 0.0_dp])
        end if
        e1 = e1/norm2(e1); e2 = cross(u, e1)
        call rng_uniform(rng, phi)
        phi = 2*pi*phi
        e = cos(phi)*e1 + sin(phi)*e2
    end subroutine random_perpendicular

    !> Advance test velocities v(3,:) by one plasma time step dt with physical
    !> (uniform-disk) encounter sampling: majorant rate n pi b_max^2 (|v| + 8 v_th)
    !> with thinning by |g|. promoted(i) flags a promoted encounter in this step.
    subroutine hybrid_beam_step(model, v, temperature, density, dt, rng, promoted, &
                                nenc, ok)
        type(hybrid_model_t), intent(inout) :: model
        real(dp), intent(inout) :: v(:, :)
        real(dp), intent(in) :: temperature, density, dt
        type(rng_t), intent(inout) :: rng
        logical, intent(out) :: promoted(:)
        integer, intent(inout) :: nenc
        logical, intent(out) :: ok
        type(encounter_record_t) :: rec
        real(dp) :: vth, gm, rate, t, u, vf(3), gvec(3), gout(3), bvec(3), b
        integer :: i, c
        ok = .true.; promoted = .false.
        vth = sqrt(temperature/model%mf)
        do i = 1, size(v, 2)
            gm = norm2(v(:, i)) + 8*vth
            rate = density*pi*model%bmax**2*gm
            t = 0
            do
                call rng_uniform(rng, u)
                t = t - log(max(u, tiny(1.0_dp)))/rate
                if (t > dt) exit
                do c = 1, 3
                    call rng_normal(rng, vf(c))
                end do
                gvec = v(:, i) - vth*vf
                call rng_uniform(rng, u)
                if (u*gm > norm2(gvec)) cycle
                call rng_uniform(rng, u)
                b = model%bmax*sqrt(u)
                call random_perpendicular(rng, gvec/norm2(gvec), bvec)
                call hybrid_encounter(model, gvec, b*bvec, gout, rec, ok)
                if (.not. ok) return
                nenc = nenc + 1
                if (rec%promoted) promoted(i) = .true.
                v(:, i) = v(:, i) + (model%mu/model%mt)*(gout - gvec)
                gm = max(gm, norm2(v(:, i)) + 8*vth); rate = density*pi*model%bmax**2*gm
            end do
        end do
    end subroutine hybrid_beam_step
end module kin6d_hybrid_encounter
