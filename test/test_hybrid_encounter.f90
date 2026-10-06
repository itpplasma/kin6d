! Independent behavioural oracles for the QG6 hybrid classical/quantum encounter.
! Oracles: exact Kepler invariants and the Rutherford deflection 2 atan(b90/b),
! lab-frame momentum/energy balances, the Chandrasekhar/Landau slowing rate for
! a constant Coulomb logarithm, kernel quadratures of the cut models (Monte
! Carlo agreement within CLT intervals), the eta >> 1 classical limit with
! forced promotion (paired samples) and central finite differences.
program test_hybrid_encounter
    use, intrinsic :: iso_fortran_env, only: dp => real64, i64 => int64
    use kin6d_species, only: atomic_units, electron
    use kin6d_hybrid_encounter
    use kin6d_quantum_collision_leaves, only: qc_chandrasekhar, qc_slowing_rate
    implicit none
    integer :: failures = 0
    real(dp), parameter :: pi = acos(-1.0_dp)
    call test_kepler_maps()
    call test_classical_encounter()
    call test_promoted_encounter()
    call test_kernel_rate()
    call test_monte_carlo()
    call test_classical_limit()
    if (failures > 0) then
        print '(i0,a)', failures, ' hybrid encounter checks failed'
        error stop 1
    end if
    print '(a)', 'all hybrid encounter checks passed'
contains
    subroutine check(name, condition, value)
        character(*), intent(in) :: name
        logical, intent(in) :: condition
        real(dp), intent(in) :: value
        if (condition) then
            print '(a,a,es11.3)', 'PASS ', name//' = ', value
        else
            print '(a,a,es11.3)', 'FAIL ', name//' = ', value
            failures = failures + 1
        end if
    end subroutine check

    subroutine model_at(model, bmax, cfg)
        type(hybrid_model_t), intent(out) :: model
        real(dp), intent(in) :: bmax
        type(hybrid_config_t), intent(in) :: cfg
        logical :: ok
        call hybrid_init(model, atomic_units(), electron(), electron(), bmax, cfg, ok)
        if (.not. ok) error stop 'hybrid_init failed'
    end subroutine model_at

    ! Lift a hyperbola state, place it at another radius on both branches and
    ! lift again; the asymptotic deflection equals 2 atan(b90/b), b90 = (g/mu)/v^2.
    subroutine test_kepler_maps()
        type(hybrid_model_t) :: model
        type(hybrid_config_t) :: cfg
        real(dp) :: x(3), v(3), vinf, h(3), a(3), uin(3), uout(3), x2(3), v2(3)
        real(dp) :: vinf2, h2(3), a2(3), uin2(3), uout2(3), err, defl, b, b90
        integer :: br
        logical :: ok
        call model_at(model, 1.0e3_dp, cfg)
        x = [0.7_dp, -0.4_dp, 2.1_dp]; v = [0.3_dp, 1.1_dp, -0.6_dp]
        call lift_asymptote(model, x, v, vinf, h, a, uin, uout)
        err = 0
        do br = -1, 1, 2
            call orbit_point(model, a, h, vinf, 9.0_dp, br, x2, v2, ok)
            call lift_asymptote(model, x2, v2, vinf2, h2, a2, uin2, uout2)
            err = max(err, abs(vinf2 - vinf), maxval(abs(h2 - h)), maxval(abs(a2 - a)), &
                      maxval(abs(uin2 - uin)), maxval(abs(uout2 - uout)))
            err = max(err, abs(norm2(x2) - 9))
            if (.not. (ok .and. br*dot_product(x2, v2) > 0)) err = huge(1.0_dp)
        end do
        call check('lift/project invariants on both branches', err < 1.0e-12_dp, err)
        b = norm2(h)/vinf; b90 = (model%g/model%mu)/vinf**2
        defl = atan2(norm2([uin(2)*uout(3) - uin(3)*uout(2), uin(3)*uout(1) - &
                            uin(1)*uout(3), uin(1)*uout(2) - uin(2)*uout(1)]), &
                     dot_product(uin, uout))
        call check('lifted asymptotes: Rutherford deflection', &
                   abs(defl - 2*atan(b90/b)) < 1.0e-12_dp, defl - 2*atan(b90/b))
    end subroutine test_kepler_maps

    ! Classical branch (KDK, encounter-resolved steps) against Rutherford.
    subroutine test_classical_encounter()
        type(hybrid_model_t) :: model
        type(hybrid_config_t) :: cfg
        type(encounter_record_t) :: rec
        real(dp) :: gout(3), b, err, ed
        integer :: i
        logical :: ok
        cfg%promote = .false.
        call model_at(model, 1.0e4_dp, cfg)
        err = 0; ed = 0
        do i = 0, 12
            b = 0.5_dp*10.0_dp**(-2 + 0.5_dp*i)    ! b90 = 0.5 at |g| = 2
            call hybrid_encounter(model, [0.0_dp, 0.0_dp, 2.0_dp], [b, 0.0_dp, 0.0_dp], &
                                  gout, rec, ok)
            if (.not. ok) error stop 'classical encounter failed'
            err = max(err, abs(rec%theta/(2*atan(0.5_dp/b)) - 1))
            ed = max(ed, abs(rec%energy_defect))
        end do
        call check('classical encounter vs Rutherford, b/b90 in [1e-2,1e4]', &
                   err < 2.0e-4_dp, err)
        call check('classical encounter energy defect', ed < 1.0e-5_dp, ed)
    end subroutine test_classical_encounter

    ! Promoted encounter: outgoing deflection equals the ring deflection, |g| and
    ! the lab-frame pair momentum/energy are conserved, the angular-momentum
    ! change is at most hbar; the naive velocity rotation double counts.
    subroutine test_promoted_encounter()
        type(hybrid_model_t) :: model
        type(hybrid_config_t) :: cfg
        type(encounter_record_t) :: rec, naive
        real(dp) :: g(3), gout(3), vt(3), vf(3), vt2(3), vf2(3), p0(3), p1(3), e0, e1
        real(dp) :: mt, mf, mu
        logical :: ok
        cfg%eta_switch = 10
        call model_at(model, 1.0e3_dp, cfg)
        vt = [0.1_dp, -0.2_dp, 2.3_dp]; vf = [0.05_dp, 0.1_dp, 0.3_dp]; g = vt - vf
        call hybrid_encounter(model, g, [0.3_dp, -0.1_dp, -0.0225_dp], &
                              gout, rec, ok)
        call check('promoted: deflection = ring deflection chi_l', &
                   ok .and. rec%promoted .and. abs(rec%theta/rec%chi - 1) < 1.0e-3_dp, &
                   rec%theta/rec%chi - 1)
        mt = model%mt; mf = model%mf; mu = model%mu
        vt2 = vt + (mu/mt)*(gout - g); vf2 = vf - (mu/mf)*(gout - g)
        p0 = mt*vt + mf*vf; p1 = mt*vt2 + mf*vf2
        e0 = (mt*sum(vt**2) + mf*sum(vf**2))/2; e1 = (mt*sum(vt2**2) + mf*sum(vf2**2))/2
        call check('promoted: lab momentum conserved', maxval(abs(p1 - p0)) < 1.0e-14_dp, &
                   maxval(abs(p1 - p0)))
        call check('promoted: lab energy conserved', abs(e1/e0 - 1) < 1.0e-14_dp, e1/e0 - 1)
        call check('promoted: angular momentum change <= hbar', &
                   rec%angular_change <= 1 + 1.0e-12_dp, rec%angular_change)
        model%cfg%naive = .true.
        call hybrid_encounter(model, g, [0.3_dp, -0.1_dp, -0.0225_dp], &
                              gout, naive, ok)
        print '(a,f9.5)', '  naive rotation at r_det, theta/chi - 1 = ', naive%theta/naive%chi - 1
        call check('naive rotation double counts outer legs (> 2%)', &
                   naive%theta/naive%chi - 1 > 0.02_dp, naive%theta/naive%chi - 1)
    end subroutine test_promoted_encounter

    ! Kernel quadrature with a constant Coulomb logarithm reproduces the
    ! Chandrasekhar/Landau rate (1 + m/mf) 4 pi n g^2 lnL Psi(x)/(m^2 v^3).
    subroutine test_kernel_rate()
        type(hybrid_model_t) :: model
        type(hybrid_config_t) :: cfg
        real(dp) :: nu, psi, nus, vb, t, n
        t = 0.7_dp; n = 1.0e-6_dp; vb = 1.9_dp
        call model_at(model, 50.0_dp, cfg)
        nu = kernel_slowing_rate(model, vb, t, n, -1, 7.5_dp)
        call qc_chandrasekhar(vb/sqrt(2*t), psi)
        call qc_slowing_rate(n, model%g, model%mt, model%mf, 7.5_dp, vb, psi, nus)
        call check('kernel quadrature = Chandrasekhar rate (constant lnL)', &
                   abs(nu/nus - 1) < 1.0e-10_dp, nu/nus - 1)
    end subroutine test_kernel_rate

    subroutine stats(xs, mean, se)
        real(dp), intent(in) :: xs(:)
        real(dp), intent(out) :: mean, se
        mean = sum(xs)/size(xs)
        se = sqrt(sum((xs - mean)**2)/(size(xs) - 1)/size(xs))
    end subroutine stats

    ! Monte Carlo slowing rates (classical, hybrid) against the kernel
    ! quadrature of the same model; |difference| < 3 standard errors (CLT).
    ! Hybrid transfer: interface defect against the all-quantum ring model and
    ! d lnL/d eta against differences in the coupling g (kappa, L fixed).
    subroutine test_monte_carlo()
        type(hybrid_model_t) :: model
        type(hybrid_config_t) :: cfg
        real(dp) :: xs(4000), t, n, vb, lam, mc, se, ex, ed, ac, l1, l2, d1, dum, h, g0
        integer :: np
        logical :: ok
        t = 1; n = 1.0e-6_dp; vb = 2*sqrt(2*t); lam = sqrt(t/(4*pi*n))
        cfg%promote = .false.
        call model_at(model, lam, cfg)
        call hybrid_rate_samples(model, vb, t, n, 11_i64, xs, np, ed, ac, ok)
        call stats(xs, mc, se)
        ex = kernel_slowing_rate(model, vb, t, n, MODEL_CLASSICAL)
        call check('MC classical vs kernel (|d|/se < 3, statistical)', &
                   ok .and. abs(mc - ex) < 3*se, abs(mc - ex)/se)
        model%cfg%promote = .true.
        call hybrid_rate_samples(model, vb, t, n, 12_i64, xs, np, ed, ac, ok)
        call stats(xs, mc, se)
        ex = kernel_slowing_rate(model, vb, t, n, MODEL_HYBRID)
        call check('MC hybrid vs kernel (|d|/se < 3, statistical)', &
                   ok .and. np > 0 .and. abs(mc - ex) < 3*se, abs(mc - ex)/se)
        call hybrid_transfer_log(model, 1.7_dp, MODEL_HYBRID, l1, d1)
        call hybrid_transfer_log(model, 1.7_dp, MODEL_QUANTUM, l2, dum)
        print '(a,es10.3)', '  hybrid - all-quantum lnL at v=1.7: ', l1 - l2
        call check('hybrid vs all-quantum transfer (interface defect)', &
                   abs(l1 - l2) < 2.0e-2_dp, l1 - l2)
        h = 1.0e-6_dp; g0 = model%g
        model%g = g0*(1 + h)
        call hybrid_transfer_log(model, 1.7_dp, MODEL_HYBRID, l2, dum)
        model%g = g0*(1 - h)
        call hybrid_transfer_log(model, 1.7_dp, MODEL_HYBRID, ex, dum)
        model%g = g0
        ! eta = g/(hbar v): d lnL = (d lnL/d eta) eta (dg/g).
        call check('hybrid d lnL/d eta vs FD in g', &
                   abs((l2 - ex)/(2*h*g0/(model%hbar*1.7_dp)) - d1) < 1.0e-6_dp*abs(d1), &
                   (l2 - ex)/(2*h*g0/(model%hbar*1.7_dp)) - d1)
    end subroutine test_monte_carlo

    ! eta >> 1: no promotion by default (bitwise classical); forced promotion
    ! agrees with the classical rate within paired statistics and the kernel
    ! formulas agree to O(1/eta^2).
    subroutine test_classical_limit()
        type(hybrid_model_t) :: model
        type(hybrid_config_t) :: cfg
        real(dp) :: xc(3000), xh(3000), xd(3000), t, n, vb, lam, mc, se, md, sd
        real(dp) :: ed, ac, kc, kh
        integer :: np
        logical :: ok
        t = 1.0e-3_dp; n = 1.0e-12_dp; vb = 2*sqrt(2*t); lam = sqrt(t/(4*pi*n))
        cfg%promote = .false.
        call model_at(model, lam, cfg)
        call hybrid_rate_samples(model, vb, t, n, 21_i64, xc, np, ed, ac, ok)
        model%cfg%promote = .true.
        call hybrid_rate_samples(model, vb, t, n, 21_i64, xh, np, ed, ac, ok)
        call check('eta>>1, default switch: no promotion, identical samples', &
                   np == 0 .and. maxval(abs(xh - xc)) == 0, real(np, dp))
        model%cfg%eta_switch = 1.0e3_dp
        call hybrid_rate_samples(model, vb, t, n, 21_i64, xh, np, ed, ac, ok)
        xd = xh - xc
        call stats(xd, md, sd)
        call stats(xc, mc, se)
        print '(a,i0,a,es10.3,a,es10.3)', '  forced promotions ', np, ', paired rel diff ', &
            md/mc, ' +- ', sd/mc
        call check('eta>>1 forced hybrid == classical (paired, |d|/se < 3)', &
                   np > 0 .and. abs(md) < 3*sd, abs(md)/sd)
        kc = kernel_slowing_rate(model, vb, t, n, MODEL_CLASSICAL)
        kh = kernel_slowing_rate(model, vb, t, n, MODEL_HYBRID)
        call check('eta>>1 kernel hybrid vs classical', abs(kh/kc - 1) < 5.0e-3_dp, kh/kc - 1)
    end subroutine test_classical_limit
end program test_hybrid_encounter
