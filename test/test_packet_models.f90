! Independent behavioural oracles for the Q2 packet models.
! Oracles: binary128 Gauss-Legendre quadrature of the complex Boys integrals;
! direct 3D quadrature of overlap, moment, kinetic (gradient form) and Coulomb
! pair elements; analytic free-Gaussian spreading; exact harmonic-oscillator
! Gaussian motion (TDVP is exact for quadratic H, also for a parity-tied sum);
! the Q1 partial-wave reference (bound validity, K enrichment); energy/norm
! conservation; finite differences of the bin JVP.
program test_packet_models
    use, intrinsic :: iso_fortran_env, only: dp => real64, qp => real128, i64 => int64
    use kin6d_gaussian_pair
    use kin6d_packet_models
    use kin6d_coulomb_encounter
    use kin6d_coulomb_encounter_leaves, only: ce_free_width
    use fortnum_quadrature, only: gauss_legendre
    implicit none
    integer :: failures = 0
    real(dp), parameter :: pi = acos(-1.0_dp)
    type(encounter_t) :: enc
    real(dp) :: tfin
    call test_boys()
    call test_pair_elements()
    call test_free_gaussian()
    call test_harmonic()
    call setup(1.0_dp, 3.0_dp, 0.0_dp, enc, tfin)
    call test_engine_cross_checks()
    call test_reference_comparison()
    call test_bin_jvp()
    if (failures > 0) then
        print '(i0,a)', failures, ' packet model checks failed'
        error stop 1
    end if
    print '(a)', 'all packet model checks passed'
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

    ! Q1 setup in atomic-like units (hbar = m = g = 1): eta, sigma/lambdabar, b/lambdabar.
    subroutine setup(eta, sig_l, b_l, e, t)
        real(dp), intent(in) :: eta, sig_l, b_l
        type(encounter_t), intent(out) :: e
        real(dp), intent(out) :: t
        type(encounter_config_t) :: cfg
        character(:), allocatable :: msg
        logical :: ok
        cfg%g = 1; cfg%m1 = 1; cfg%m2 = 1; cfg%k = 0.5_dp/eta
        cfg%sigma = sig_l/cfg%k; cfg%z0 = 6*cfg%sigma; cfg%b = b_l/cfg%k
        t = 2*cfg%z0/(cfg%hbar*cfg%k/0.5_dp); cfg%t_final = t
        call encounter_init(e, cfg, ok, msg)
        if (.not. ok) error stop 'encounter setup failed: '//msg
    end subroutine setup

    ! Binary128 Gauss-Legendre nodes on [-1,1] by Newton iteration.
    subroutine gl_qp(n, x, w)
        integer, intent(in) :: n
        real(qp), intent(out) :: x(n), w(n)
        real(qp) :: z, p0, p1, p2, dpz
        integer :: i, j, it
        do i = 1, n
            z = cos(acos(-1.0_qp)*(i - 0.25_qp)/(n + 0.5_qp))
            do it = 1, 100
                p0 = 1; p1 = z
                do j = 2, n
                    p2 = ((2*j - 1)*z*p1 - (j - 1)*p0)/j; p0 = p1; p1 = p2
                end do
                dpz = n*(z*p1 - p0)/(z*z - 1)
                z = z - p1/dpz
                if (abs(p1/dpz) < 1.0e-32_qp) exit
            end do
            x(i) = z; w(i) = 2/((1 - z*z)*dpz**2)
        end do
    end subroutine gl_qp

    ! E_n(T) = int_0^1 u^2n exp(T(1-u^2)) du vs composite binary128 GL.
    subroutine test_boys()
        type(pair_rule_t) :: rule
        real(qp) :: xq(32), wq(32), u, a, b
        complex(qp) :: ref(0:2), tq
        complex(dp) :: t, tau, f(0:2)
        real(dp), parameter :: mags(9) = [1.0e-3_dp, 0.5_dp, 5.0_dp, 20.0_dp, 35.9_dp, &
                                          36.1_dp, 60.0_dp, 200.0_dp, 900.0_dp]
        real(dp), parameter :: args(7) = [0.0_dp, 0.25_dp, 0.5_dp, 0.75_dp, 0.99_dp, &
                                          -0.5_dp, -0.9_dp]
        real(dp) :: err, scale
        integer :: im, ia, p, q, n
        call pair_rule_init(rule, 96)
        call gl_qp(32, xq, wq)
        err = 0
        do im = 1, size(mags)
            do ia = 1, size(args)
                t = mags(im)*exp(cmplx(0, pi*args(ia), dp)); tq = t
                ref = 0
                do p = 1, 400
                    a = real(p - 1, qp)/400; b = real(p, qp)/400
                    do q = 1, 32
                        u = (a + b)/2 + (b - a)/2*xq(q)
                        do n = 0, 2
                            ref(n) = ref(n) + (b - a)/2*wq(q)*u**(2*n)*exp(tq*(1 - u*u))
                        end do
                    end do
                end do
                call boys_scaled(rule, t, tau, f)
                ! Compare E_n exp(-tau) on the bounded scale.
                ref = ref*exp(-cmplx(tau, kind=qp))
                scale = maxval(abs(real(ref(0:0), dp)))
                do n = 0, 2
                    err = max(err, real(abs(f(n) - ref(n)), dp)/max(scale, 1.0e-300_dp))
                end do
            end do
        end do
        call check('complex Boys E_0..2 vs binary128 quadrature (63 T)', err < 1.0e-12_dp, err)
    end subroutine test_boys

    ! Pair elements vs direct spherical quadrature about the Coulomb centre
    ! (kinetic term in the independent gradient form (hbar^2/2m) int grad e* . grad g).
    subroutine test_pair_elements()
        type(gaussian_t) :: gk, gl
        integer :: icase
        real(dp) :: err_s, err_t, err_v
        err_s = 0; err_t = 0; err_v = 0
        do icase = 1, 2
            gk%q = [0.7_dp, -0.4_dp]; gk%kappa = [0.3_dp, 1.1_dp]
            gk%alpha = (0.6_dp, -0.25_dp)
            gl%q = [-0.2_dp, 0.5_dp]; gl%kappa = [-0.5_dp, 0.4_dp]
            gl%alpha = (0.9_dp, 0.3_dp)
            if (icase == 2) then     ! far from the origin: |T| > 36 branch
                gk%q = [5.0_dp, -6.0_dp]; gl%q = [5.5_dp, -5.0_dp]
            end if
            call one_case(gk, gl, err_s, err_t, err_v)
        end do
        call check('pair overlap and moments vs 3D quadrature', err_s < 1.0e-10_dp, err_s)
        call check('pair kinetic+harmonic vs gradient-form quadrature', err_t < 1.0e-10_dp, &
                   err_t)
        call check('pair Coulomb (erf/Boys) vs 3D quadrature', err_v < 1.0e-10_dp, err_v)
    end subroutine test_pair_elements

    subroutine one_case(gk, gl, err_s, err_t, err_v)
        type(gaussian_t), intent(in) :: gk, gl
        real(dp), intent(inout) :: err_s, err_t, err_v
        integer, parameter :: nr = 200, nt = 80, np = 96
        real(dp), parameter :: hbar = 1.3_dp, mass = 0.7_dp, kap = 0.45_dp
        type(pair_rule_t) :: rule
        real(dp) :: xr(nr), wr(nr), xt(nt), wt(nt), rmax, r, st, ph, x, y, z, w, uk(3), ul(3)
        complex(dp) :: s, sig(4, 4), hel(4), hk(4), hv(4), ssig(4, 4), tref(4), vref(4)
        complex(dp) :: e(4), el(4), ge(3, 4), gg(3), vk, vl, sref(4, 4)
        integer :: i, j, q, a, b
        call pair_rule_init(rule, 96)
        call pair_elements(rule, gk, gl, hbar, mass, 0.0_dp, kap, s, sig, hk)
        call pair_elements(rule, gk, gl, hbar, 0.0_dp, 1.0_dp, 0.0_dp, s, sig, hv)
        ssig = s*sig
        rmax = max(norm2(gk%q), norm2(gl%q)) + 9
        call gauss_legendre(nr, xr, wr); xr = (xr + 1)*rmax/2; wr = wr*rmax/2
        call gauss_legendre(nt, xt, wt)
        sref = 0; tref = 0; vref = 0
        do i = 1, nr
            r = xr(i)
            do j = 1, nt
                st = sqrt(1 - xt(j)**2)
                do q = 1, np
                    ph = 2*pi*(q - 1)/np
                    x = r*st*cos(ph); y = r*st*sin(ph); z = r*xt(j)
                    w = wr(i)*wt(j)*2*pi/np*r*r
                    uk = [x - gk%q(1), y, z - gk%q(2)]; ul = [x - gl%q(1), y, z - gl%q(2)]
                    vk = gaussian_value(gk, x, y, z); vl = gaussian_value(gl, x, y, z)
                    e = [vk, uk(1)*vk, uk(3)*vk, sum(uk**2)*vk]
                    ! gradients of the bra elementary functions and of g_l
                    gg = vk*(-2*gk%alpha*uk + cmplx(0, 1, dp)*[gk%kappa(1), 0.0_dp, gk%kappa(2)])
                    ge(:, 1) = gg
                    ge(:, 2) = uk(1)*gg + [vk, (0.0_dp, 0.0_dp), (0.0_dp, 0.0_dp)]
                    ge(:, 3) = uk(3)*gg + [(0.0_dp, 0.0_dp), (0.0_dp, 0.0_dp), vk]
                    ge(:, 4) = sum(uk**2)*gg + 2*uk*vk
                    gg = vl*(-2*gl%alpha*ul + cmplx(0, 1, dp)*[gl%kappa(1), 0.0_dp, gl%kappa(2)])
                    el = [vl, ul(1)*vl, ul(3)*vl, sum(ul**2)*vl]
                    do a = 1, 4
                        do b = 1, 4
                            sref(a, b) = sref(a, b) + w*conjg(e(a))*el(b)
                        end do
                        tref(a) = tref(a) + w*(hbar**2/(2*mass))*sum(conjg(ge(:, a))*gg) &
                                  + w*conjg(e(a))*(kap/2)*r*r*vl
                        vref(a) = vref(a) + w*conjg(e(a))*vl/r
                    end do
                end do
            end do
        end do
        err_s = max(err_s, maxval(abs(ssig - sref))/maxval(abs(sref)))
        err_t = max(err_t, maxval(abs(hk - tref))/maxval(abs(tref)))
        err_v = max(err_v, maxval(abs(hv - vref))/maxval(abs(vref)))
    end subroutine one_case

    subroutine run(e, spec, times, smp, tr)
        type(encounter_t), intent(in) :: e
        type(packet_spec_t), intent(in) :: spec
        real(dp), intent(in) :: times(:)
        type(packet_sample_t), allocatable, intent(out) :: smp(:)
        type(packet_trace_t), intent(out), optional :: tr
        character(:), allocatable :: msg
        logical :: ok
        integer(i64) :: nfev
        call packet_samples(e, spec, times, smp, ok, msg, nfev)
        if (.not. ok) error stop 'packet model failed: '//msg
        if (present(tr)) call packet_trace(e, spec, smp, tr)
    end subroutine run

    function linspace(a, b, n) result(x)
        real(dp), intent(in) :: a, b
        integer, intent(in) :: n
        real(dp) :: x(n)
        integer :: i
        x = [(a + (b - a)*(i - 1)/real(n - 1, dp), i=1, n)]
    end function linspace

    ! Free relative motion (g = 0): TDVP Gaussian = exact free Gaussian.
    subroutine test_free_gaussian()
        type(encounter_config_t) :: cfg
        type(encounter_t) :: e
        type(packet_spec_t) :: spec
        type(packet_sample_t), allocatable :: smp(:)
        type(packet_trace_t) :: tr
        character(:), allocatable :: msg
        real(dp) :: t, w, err
        logical :: ok
        integer :: i
        cfg%g = 0; cfg%k = 0.5_dp; cfg%sigma = 6; cfg%z0 = 30; cfg%b = 4
        t = 2*cfg%z0/(cfg%k/0.5_dp); cfg%t_final = t
        call encounter_init(e, cfg, ok, msg)
        if (.not. ok) error stop msg
        spec%model = MODEL_M3
        call run(e, spec, linspace(0.0_dp, t, 9), smp, tr)
        err = 0
        do i = 1, 9
            call ce_free_width(cfg%sigma, e%mu, cfg%hbar, smp(i)%t, w)
            err = max(err, abs(1/(4*real(smp(i)%g(1)%alpha)) - w**2)/w**2)
        end do
        call check('free TDVP width vs analytic spreading', err < 1.0e-9_dp, err)
        call check('free Gaussian vs Q1 g=0 reference (discrete error)', &
                   maxval(tr%error) < 1.0e-5_dp, maxval(tr%error))
        call check('free Gaussian full-H residual bound', tr%bound(9) < 1.0e-4_dp, tr%bound(9))
    end subroutine test_free_gaussian

    ! Isotropic oscillator kappa r^2/2 (g = 0): a thawed Gaussian, and a
    ! parity-tied Gaussian pair, evolve exactly (analytic centroid and width).
    subroutine test_harmonic()
        type(encounter_config_t) :: cfg
        type(encounter_t) :: e
        type(packet_spec_t) :: spec
        type(packet_sample_t), allocatable :: smp(:)
        character(:), allocatable :: msg
        real(dp) :: om, times(13), qa(2), sa, errq, errs, erre, errn, c, sn
        logical :: ok
        integer :: i, im
        integer, parameter :: models(4) = [MODEL_M3, MODEL_M3, MODEL_M3S, MODEL_M3T]
        cfg%g = 0; cfg%k = 0.5_dp; cfg%sigma = 3; cfg%z0 = 10; cfg%b = 4; cfg%t_final = 1
        call encounter_init(e, cfg, ok, msg)
        if (.not. ok) error stop msg
        spec%kappa = 0.02_dp
        om = sqrt(spec%kappa/e%mu)
        times = linspace(0.0_dp, 2*pi/om, 13)
        errq = 0; errs = 0; erre = 0; errn = 0
        do im = 1, 4
            spec%model = models(im); spec%metric_engine = im == 2
            call run(e, spec, times, smp)
            do i = 1, size(times)
                c = cos(om*times(i)); sn = sin(om*times(i))
                qa = [cfg%b, -cfg%z0]*c + [0.0_dp, cfg%k]/(e%mu*om)*sn
                sa = cfg%sigma**2*c**2 + (1/(2*e%mu*om))**2/cfg%sigma**2*sn**2
                errq = max(errq, norm2(smp(i)%g(1)%q - qa)/cfg%sigma)
                errs = max(errs, abs(1/(4*real(smp(i)%g(1)%alpha)) - sa)/sa)
                erre = max(erre, abs(smp(i)%energy - smp(1)%energy)/abs(smp(1)%energy))
                errn = max(errn, abs(smp(i)%norm - 1))
            end do
        end do
        call check('harmonic TDVP centroid exact (M3, engine, M3S, M3T)', errq < 1.0e-8_dp, errq)
        call check('harmonic TDVP width exact', errs < 1.0e-8_dp, errs)
        call check('harmonic energy conservation', erre < 1.0e-9_dp, erre)
        call check('harmonic norm conservation', errn < 1.0e-8_dp, errn)
    end subroutine test_harmonic

    ! Canonical (FortSym) and metric-engine TDVP are two implementations of the
    ! same equations; conservation of energy and norm along Coulomb runs.
    subroutine test_engine_cross_checks()
        type(packet_spec_t) :: spec
        type(packet_sample_t), allocatable :: a(:), b(:)
        real(dp) :: times(9), d, erre, errn
        integer :: im, i
        integer, parameter :: models(4) = [MODEL_M4, MODEL_M2S, MODEL_M2T, MODEL_M2]
        times = linspace(0.0_dp, tfin, 9)
        spec%model = MODEL_M3
        call run(enc, spec, times, a)
        spec%metric_engine = .true.
        call run(enc, spec, times, b)
        d = max(norm2(a(9)%g(1)%q - b(9)%g(1)%q)/enc%cfg%sigma, &
                abs(a(9)%g(1)%alpha - b(9)%g(1)%alpha)/abs(a(9)%g(1)%alpha), &
                abs(a(9)%g(1)%d - b(9)%g(1)%d)/abs(a(9)%g(1)%d))
        call check('M3 canonical vs metric engine (eta=1)', d < 1.0e-7_dp, d)
        spec%model = MODEL_M2; spec%metric_engine = .false.
        spec%rtol = 1.0e-12_dp; spec%atol = 1.0e-14_dp
        call run(enc, spec, times, a)
        spec%metric_engine = .true.
        call run(enc, spec, times, b)
        d = max(norm2(a(9)%g(1)%q - b(9)%g(1)%q)/enc%cfg%sigma, &
                abs(a(9)%g(1)%alpha - b(9)%g(1)%alpha)/abs(a(9)%g(1)%alpha), &
                abs(a(9)%g(1)%d - b(9)%g(1)%d)/abs(a(9)%g(1)%d), &
                abs(a(9)%cm_overlap - b(9)%cm_overlap))
        call check('M2 canonical vs metric engine with CM tie', d < 1.0e-7_dp, d)
        call check('M2 symmetric preparation keeps s1 = s2', a(9)%width_asymmetry < 1.0e-12_dp, &
                   a(9)%width_asymmetry)
        spec%metric_engine = .false.; spec%k = 4
        spec%rtol = 1.0e-10_dp; spec%atol = 1.0e-12_dp
        erre = 0; errn = 0
        do im = 1, 4
            spec%model = models(im)
            call run(enc, spec, times, a)
            do i = 1, 9
                erre = max(erre, abs(a(i)%energy - a(1)%energy)/abs(a(1)%energy))
                errn = max(errn, abs(a(i)%norm - 1))
            end do
        end do
        call check('TDVP energy conservation (M4 K=4, M2S, M2T, M2)', erre < 1.0e-7_dp, erre)
        call check('TDVP norm conservation', errn < 1.0e-8_dp, errn)
    end subroutine test_engine_cross_checks

    ! Against the Q1 reference: identical preparation, Duhamel bound >= actual
    ! discrete error at every sample, enrichment K=1 -> 4 reduces the error.
    subroutine test_reference_comparison()
        type(packet_spec_t) :: spec
        type(packet_sample_t), allocatable :: smp(:)
        type(packet_trace_t) :: tr
        real(dp) :: times(41), e0, slack, ek(3)
        integer :: im
        integer, parameter :: models(6) = [MODEL_M3, MODEL_M4, MODEL_M4, MODEL_M2S, &
                                           MODEL_M2T, MODEL_M2]
        integer, parameter :: ks(6) = [1, 2, 4, 1, 1, 1]
        times = linspace(0.0_dp, tfin, 41)
        e0 = 0; slack = huge(1.0_dp)
        do im = 1, 6
            spec%model = models(im); spec%k = ks(im)
            call run(enc, spec, times, smp, tr)
            e0 = max(e0, tr%initial_error)
            slack = min(slack, minval(tr%bound - tr%error), minval(tr%total_bound - tr%total_error))
            if (im <= 3) ek(im) = tr%error(41)
            print '(a,a,5es11.3)', '  ', model_name(spec)//' error, bound, phase-free, slacks: ', &
                tr%error(41), tr%bound(41), tr%phase_error(41), minval(tr%bound - tr%error), &
                minval(tr%total_bound - tr%total_error)
        end do
        print '(a,es11.3)', '  Q1 spectral-window e0 = ', enc%ledger%e0
        ! Sector normalization (norm ~ 1/2) amplifies the window defect by <= sqrt(2).
        call check('model preparation equals Q1 state (<= sqrt2 window e0 + 1e-9)', &
                   e0 <= sqrt(2.0_dp)*enc%ledger%e0 + 1.0e-9_dp, e0)
        call check('Duhamel full-H bound >= actual error (all samples)', slack > -1.0e-9_dp, &
                   slack)
        call check('enrichment K=1 -> 2 -> 4 reduces the state error', &
                   ek(2) <= ek(1) .and. ek(3) < ek(1), ek(3)/ek(1))
    end subroutine test_reference_comparison

    ! JVP of an outgoing bin probability (M3) vs central differences, frozen basis.
    subroutine test_bin_jvp()
        type(packet_spec_t) :: spec
        type(encounter_t) :: ep, eb
        real(dp) :: p, dpv(2), fd(2), pp, pm, dum, h(2), err, tb
        logical :: ok
        integer :: j
        call setup(1.0_dp, 3.0_dp, 2.0_dp, eb, tb)
        spec%model = MODEL_M3; spec%rtol = 1.0e-12_dp; spec%atol = 1.0e-14_dp
        h = [1.0e-3_dp*eb%cfg%sigma, 1.0e-4_dp*eb%cfg%sigma]
        err = 0
        call packet_bin_jvp(eb, spec, tb, SECTOR_PRODUCT, eb%cfg%z0/4, pi/6, pi/3, &
                            [1.0_dp, 0.0_dp], p, dpv(1), ok)
        call packet_bin_jvp(eb, spec, tb, SECTOR_PRODUCT, eb%cfg%z0/4, pi/6, pi/3, &
                            [0.0_dp, 1.0_dp], p, dpv(2), ok)
        do j = 1, 2
            ep = eb
            if (j == 1) ep%cfg%b = eb%cfg%b + h(1)
            if (j == 2) ep%cfg%sigma = eb%cfg%sigma + h(2)
            call packet_bin_jvp(ep, spec, tb, SECTOR_PRODUCT, eb%cfg%z0/4, pi/6, pi/3, &
                                [0.0_dp, 0.0_dp], pp, dum, ok)
            ep = eb
            if (j == 1) ep%cfg%b = eb%cfg%b - h(1)
            if (j == 2) ep%cfg%sigma = eb%cfg%sigma - h(2)
            call packet_bin_jvp(ep, spec, tb, SECTOR_PRODUCT, eb%cfg%z0/4, pi/6, pi/3, &
                                [0.0_dp, 0.0_dp], pm, dum, ok)
            fd(j) = (pp - pm)/(2*h(j))
            err = max(err, abs(dpv(j) - fd(j))/max(abs(fd(j)), 1.0e-12_dp))
        end do
        print '(a,es12.4,a,2es12.4,a,2es12.4)', '  P(30-60 deg) = ', p, ' JVP ', dpv, ' FD ', fd
        call check('M3 bin JVP (b, sigma) vs central differences', err < 1.0e-5_dp, err)
    end subroutine test_bin_jvp
end program test_packet_models
