program response_strip_ladder
    !! R-B1 S-ladder of the resolved resistive strip (periodic cylinder,
    !! sharp current step) with the pre-declared decision rules:
    !!   R1 constant-psi (2/1, q0 1.2, r1 0.5, r_s Delta' = 5), S = 1e6..1e10:
    !!      |gamma/gamma_RI - 1| <= 0.02 at every S; fitted exponent of gamma
    !!      within 0.03 of -3/5; of the current HWHM within 0.03 of -2/5.
    !!   R2 non-constant-psi (1/1, q0 0.7, Delta' = infinity), S = 1e6..1e10:
    !!      |gamma/(eta F'^2)^(1/3) - 1| <= 0.02; exponents of gamma and of the
    !!      HWHM within 0.03 of -1/3.
    !!   R3 ideal limit eta = 0, gamma = 1e-6: layer current within 1e-3 of
    !!      2m psi_ext r_s^m.
    !!   R4 steady gamma = 1e-9, S = 1e6: psi_s/psi_vac within 1 % of
    !!      -2m/(r_s Delta').
    !! Reported only: mesh (hcore, ratio) and X_s spread at S = 1e8, 1e10;
    !! eta -> 0 ladder at gamma = 1e-3 against the eta = 0 strip; nodes per
    !! strip against ln S.  Cases run in parallel (OpenMP over cases).
    !$ use omp_lib, only: omp_get_max_threads
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_response_strip, only: cyl_step_t, strip_options_t, &
        strip_solution_t, step_rational, step_fprime, step_delta_prime, &
        strip_growth_rate, strip_mesh, strip_solve, layer_width, &
        ri_dispersion_gamma
    implicit none

    integer, parameter :: ns = 9, nconv = 9
    real(dp) :: slist(ns), gnum(ns, 2), gref(ns, 2), wid(ns, 2), wref(ns, 2)
    real(dp) :: conv(nconv, 2, 2)
    integer :: nodes(ns, 2), its(ns, 2), sts(ns, 2)
    type(cyl_step_t) :: eqs(2)
    type(strip_options_t) :: copt(nconv)
    logical :: ok(4)
    real(dp) :: sl(2), slw(2), dev
    integer :: k, c, idx, nthreads

    nthreads = 1
    !$ nthreads = omp_get_max_threads()
    print '(a,i0)', "response_strip_ladder threads ", nthreads
    eqs(1) = cyl_step_t(m=2, n=1, q0=1.2_dp, r1=0.5_dp)
    eqs(2) = cyl_step_t(m=1, n=1, q0=0.7_dp, r1=0.5_dp)
    do k = 1, ns
        slist(k) = 10.0_dp**(6.0_dp + 0.5_dp*(k - 1))
    end do
    copt(1) = strip_options_t()
    copt(2) = strip_options_t(hcore=0.1_dp)
    copt(3) = strip_options_t(hcore=0.025_dp)
    copt(4) = strip_options_t(ratio=1.1_dp)
    copt(5) = strip_options_t(ratio=1.025_dp)
    copt(6) = strip_options_t(xs=0.04_dp)
    copt(7) = strip_options_t(xs=0.06_dp)
    copt(8) = strip_options_t(hmax=2.0e-3_dp)
    copt(9) = strip_options_t(hmax=5.0e-4_dp)

    !$omp parallel do schedule(dynamic) private(k, c)
    do idx = 1, 2*ns
        c = 1 + (idx - 1)/ns
        k = 1 + mod(idx - 1, ns)
        call one(eqs(c), slist(k), strip_options_t(), gnum(k, c), gref(k, c), &
            wid(k, c), wref(k, c), nodes(k, c), its(k, c), sts(k, c))
    end do
    !$omp end parallel do
    !$omp parallel do schedule(dynamic) private(k, c, dev)
    do idx = 1, 4*nconv
        c = 1 + mod(idx - 1, 2)
        k = 1 + mod((idx - 1)/2, 2)
        call one(eqs(c), merge(1.0e8_dp, 1.0e10_dp, k == 1), copt(1 + (idx - 1)/4), &
            conv(1 + (idx - 1)/4, k, c), dev, dev, dev, nodes(1, 1), its(1, 1), &
            sts(1, 1))
    end do
    !$omp end parallel do

    do c = 1, 2
        print '(/,a,i0,a,i0,a,f5.2,a,es10.3)', "case m/n = ", eqs(c)%m, "/", eqs(c)%n, &
            " q0 = ", eqs(c)%q0, " r_s Delta' = ", step_delta_prime(eqs(c))
        print '(a)', "        S        gamma        gamma_ref    ratio-1    " // &
            "HWHM       width_ref  nodes it st"
        do k = 1, ns
            print '(es10.2,2es13.5,es11.2,2es11.3,i6,i3,i3)', slist(k), gnum(k, c), &
                gref(k, c), gnum(k, c)/gref(k, c) - 1.0_dp, wid(k, c), wref(k, c), &
                nodes(k, c), its(k, c), sts(k, c)
        end do
        sl(c) = slope(log10(slist), log10(gnum(:, c)))
        slw(c) = slope(log10(slist), log10(wid(:, c)))
        print '(a,f9.5,a,f9.5)', "fitted exponent gamma ", sl(c), "  HWHM ", slw(c)
        do k = 1, 2
            print '(a,es8.1,a,/,9es11.3)', "convergence S = ", merge(1.0e8_dp, 1.0e10_dp, &
                k == 1), " gamma/base-1 (hcore .1 .025 ratio 1.1 1.025 xs .04 .06 hmax 2e-3 5e-4):", &
                conv(2:, k, c)/conv(1, k, c) - 1.0_dp
        end do
    end do

    ok(1) = all(sts(:, 1) == 0) .and. maxval(abs(gnum(:, 1)/gref(:, 1) - 1.0_dp)) <= &
        0.02_dp .and. abs(sl(1) + 0.6_dp) <= 0.03_dp .and. abs(slw(1) + 0.4_dp) <= 0.03_dp
    ok(2) = all(sts(:, 2) == 0) .and. maxval(abs(gnum(:, 2)/gref(:, 2) - 1.0_dp)) <= &
        0.02_dp .and. abs(sl(2) + 1.0_dp/3) <= 0.03_dp .and. &
        abs(slw(2) + 1.0_dp/3) <= 0.03_dp
    call ideal_and_steady(ok(3), ok(4))
    print '(/,a,l1)', "R1 constant-psi FKR (gamma, S^-3/5, width S^-2/5): PASS=", ok(1)
    print '(a,l1)', "R2 non-constant-psi (gamma0, S^-1/3): PASS=", ok(2)
    print '(a,l1)', "R3 ideal limit eta = 0 within 1e-3: PASS=", ok(3)
    print '(a,l1)', "R4 steady P = -S/(r_s Delta') within 1 %: PASS=", ok(4)
    print '(a,a)', "VERDICT ", merge("PASS", "FAIL", all(ok))

contains

    subroutine one(eq, s, opt, g, gr, w, wr, nn, it, st)
        type(cyl_step_t), intent(in) :: eq
        real(dp), intent(in) :: s
        type(strip_options_t), intent(in) :: opt
        real(dp), intent(out) :: g, gr, w, wr
        integer, intent(out) :: nn, it, st
        type(strip_solution_t) :: sol
        complex(dp) :: gc
        real(dp) :: eta, rs, dpr, fp
        integer :: st2
        eta = 1.0_dp/s
        rs = step_rational(eq)
        fp = step_fprime(eq)
        dpr = step_delta_prime(eq)
        if (dpr < 0.5_dp*huge(1.0_dp)) dpr = dpr/rs
        call ri_dispersion_gamma(eta, fp, dpr, gr, st2)
        call strip_growth_rate(eq, eta, cmplx(1.3_dp*gr, 0.0_dp, dp), opt, gc, sol, &
            it, st)
        if (st2 /= 0 .or. sol%info /= 0) st = st + 10
        if (abs(aimag(gc)) > 1.0e-8_dp*abs(gc)) st = st + 100
        g = real(gc)
        w = sol%width
        nn = size(sol%r)
        if (eq%m == 1) then
            wr = (eta**2/fp**2)**(1.0_dp/6.0_dp)
        else
            wr = (gr*eta/fp**2)**0.25_dp
        end if
    end subroutine one

    pure real(dp) function slope(x, y) result(b)
        real(dp), intent(in) :: x(:), y(:)
        real(dp) :: xm, ym
        xm = sum(x)/size(x)
        ym = sum(y)/size(y)
        b = sum((x - xm)*(y - ym))/sum((x - xm)**2)
    end function slope

    subroutine ideal_and_steady(ok3, ok4)
        logical, intent(out) :: ok3, ok4
        type(cyl_step_t) :: eq
        type(strip_options_t) :: opt
        type(strip_solution_t) :: sol, sol0
        real(dp), allocatable :: r(:)
        complex(dp) :: g, p
        real(dp) :: rs, fp, ideal, eta, pref
        integer :: k
        eq = eqs(1)
        rs = step_rational(eq)
        fp = step_fprime(eq)
        ideal = 2*eq%m*real(eq%psi_ext)*rs**eq%m
        print '(/,a)', "ideal limit, eta = 0: gamma, I/I_ideal - 1, |psi_s|/psi_vac, nodes"
        do k = 4, 7
            g = cmplx(10.0_dp**(-k), 0.0_dp, dp)
            call strip_mesh(rs, layer_width(0.0_dp, g, fp), opt, r)
            call strip_solve(eq, 0.0_dp, g, r, sol)
            print '(es10.2,2es13.4,i6)', real(g), abs(sol%current/ideal - 1.0_dp), &
                abs(sol%psi_s)/rs**eq%m, size(r)
            if (k == 6) ok3 = sol%info == 0 .and. &
                abs(sol%current/ideal - 1.0_dp) <= 1.0e-3_dp
        end do
        g = (1.0e-3_dp, 0.0_dp)
        call strip_mesh(rs, layer_width(0.0_dp, g, fp), opt, r)
        call strip_solve(eq, 0.0_dp, g, r, sol0)
        print '(a,es13.5)', "eta -> 0 at gamma = 1e-3: I(eta=0)/I_ideal = ", &
            real(sol0%current/ideal)
        do k = 8, 18, 2
            eta = 10.0_dp**(-k)
            call strip_mesh(rs, layer_width(eta, g, fp), opt, r)
            call strip_solve(eq, eta, g, r, sol)
            print '(a,es9.1,a,es12.4,a,i6)', "  S ", 1.0_dp/eta, &
                "  |I(eta)/I(0) - 1| ", abs(sol%current/sol0%current - 1.0_dp), &
                " nodes ", size(r)
        end do
        eta = 1.0e-6_dp
        g = (1.0e-9_dp, 0.0_dp)
        call strip_mesh(rs, layer_width(eta, g, fp), opt, r)
        call strip_solve(eq, eta, g, r, sol)
        p = sol%psi_s/(eq%psi_ext*rs**eq%m)
        pref = -2.0_dp*eq%m/step_delta_prime(eq)
        print '(a,2es14.6,a,es14.6)', "steady P ", p, " ref ", pref
        ok4 = sol%info == 0 .and. abs(p/pref - 1.0_dp) <= 0.01_dp
    end subroutine ideal_and_steady

end program response_strip_ladder
