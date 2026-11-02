program test_rm_oracles
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use rm_mhd_oracles
    use kin6d_response_strip, only: cyl_step_t, step_delta_prime, &
        ri_dispersion_gamma
    use fortnum_roots, only: root_brent
    use fortnum_status, only: fortnum_status_t
    implicit none
    character(32) :: mode
    integer :: failures = 0
    call get_command_argument(1, mode)
    select case (trim(mode))
    case ('waves')
        call waves()
    case ('rates')
        call rates()
    case ('outer')
        call step_profile()
    case ('newcomb')
        call force_free()
    case ('convergence')
        call convergence()
    case ('saturation')
        call saturation()
    case ('table')
        call table()
    case default
        error stop 'expected waves|rates|outer|newcomb|convergence|saturation|table'
    end select
    if (failures /= 0) error stop 'RM oracle checks failed'
contains
    subroutine check(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description
        if (condition) then
            print '(2a)', 'PASS ', description
        else
            print '(2a)', 'FAIL ', description
            failures = failures + 1
        end if
    end subroutine

    subroutine waves()
        real(dp) :: k(3), b(3), roots(3), matrix(3, 3), expected(3), work(24)
        real(dp) :: rho, pressure, adiabatic, kb, b2, cs2, identity
        integer :: i, j, sample, status, info
        external :: dsyev
        rho = 2.7_dp; pressure = 0.9_dp; adiabatic = 5.0_dp/3
        cs2 = adiabatic*pressure/rho
        ! These O(1 s^-2) spectra allow LAPACK roundoff in the absolute budget.
        do sample = 1, 12
            k = [0.1_dp*sample, -0.3_dp, 0.07_dp*(sample - 3)]
            b = sqrt(oracle_mu0)*[0.2_dp, 0.11_dp*sample, -0.7_dp]
            call uniform_waves(k, b, rho, pressure, adiabatic, roots, status)
            kb = dot_product(k, b); b2 = dot_product(b, b)
            ! Independent Cartesian linearized momentum/induction operator.
            do j = 1, 3
                do i = 1, 3
                    identity = 0
                    if (i == j) identity = 1
                    matrix(i, j) = cs2*k(i)*k(j) + (kb**2*identity &
                        - kb*(b(i)*k(j) + k(i)*b(j)) + b2*k(i)*k(j)) &
                        /(oracle_mu0*rho)
                end do
            end do
            call dsyev('N', 'U', 3, matrix, 3, expected, work, size(work), info)
            call check(status == 0 .and. info == 0, 'wave operator admissible')
            call check(maxval(abs(roots([3, 1, 2]) - expected)) < 2e-13_dp, &
                'oblique waves agree with independent 3x3 force eigenvalues')
        end do
        b = [0.0_dp, 0.0_dp, sqrt(oracle_mu0*rho)]
        call uniform_waves([0.0_dp, 0.0_dp, 2.0_dp], b, rho, pressure, &
            adiabatic, roots, status)
        call check(maxval(abs(roots - [4.0_dp, 4.0_dp, 4*cs2])) < 2e-14_dp, &
            'parallel Alfven and acoustic limits')
        call uniform_waves([2.0_dp, 0.0_dp, 0.0_dp], b, rho, pressure, &
            adiabatic, roots, status)
        call check(maxval(abs(roots - [0.0_dp, 4*(1 + cs2), 0.0_dp])) < 2e-14_dp, &
            'perpendicular fast wave and two zero roots')
        call uniform_waves([1.0_dp, 2.0_dp, 3.0_dp], b, rho, 0.0_dp, &
            adiabatic, roots, status)
        call check(maxval(abs(roots - [9.0_dp, 14.0_dp, 0.0_dp])) < 2e-13_dp, &
            'cold plasma limit')
        call uniform_waves([1.0_dp, 0.0_dp, 0.0_dp], [0.0_dp, 0.0_dp, 0.0_dp], &
            rho, pressure, adiabatic, roots, status)
        call check(maxval(abs(roots - [0.0_dp, cs2, 0.0_dp])) < 2e-14_dp, &
            'unmagnetized acoustic limit')
        call uniform_waves([0.0_dp, 0.0_dp, 0.0_dp], b, rho, 0.0_dp, &
            adiabatic, roots, status)
        call check(status == 0 .and. maxval(abs(roots)) <= tiny(1.0_dp), &
            'zero wave vector')
        b = sqrt(oracle_mu0*rho)*[1.0_dp, 1.0_dp, 1.0_dp]
        k = [1.0_dp, 1.0_dp, 1.0_dp]
        call uniform_waves(k, b, rho, 3*rho/adiabatic, adiabatic, roots, status)
        call check(status == 0 .and. maxval(abs(roots - 9)) < 2e-13_dp, &
            'coincident Alfven and acoustic speeds remain finite')
        call uniform_waves(k, b, -rho, pressure, adiabatic, roots, status)
        call check(status == oracle_domain, 'negative density rejected')
    end subroutine

    subroutine rates()
        real(dp) :: eta, rho, fp, delta, fkr, kink, ri, h, width, second, speed
        integer :: status
        eta = oracle_mu0*1e-10_dp; rho = 1/oracle_mu0
        fp = 0.7_dp; delta = 0.2_dp
        call tearing_rates(eta, rho, fp, delta, fkr, kink, status)
        call ri_dispersion_gamma(eta/oracle_mu0, fp, delta, ri, status)
        ! The full layer retains a finite-Q correction at D=1e-10 m^2/s.
        call check(status == 0 .and. abs(fkr/ri - 1) < 3e-5_dp, &
            'FKR agrees with full Ara dispersion in constant-psi limit')
        h = (eta/oracle_mu0)**(1.0_dp/3)*fp**(2.0_dp/3)
        call ri_dispersion_gamma(eta/oracle_mu0, fp, huge(1.0_dp), ri, status)
        call check(abs(kink/h - 1) < 2e-14_dp .and. abs(kink/ri - 1) < 2e-14_dp, &
            'marginal-ideal internal kink agrees with Q=1 layer root')
        call tearing_rates(32*eta, 32*rho, 2*fp, delta, ri, second, status)
        call check(abs(ri/fkr - 64.0_dp**0.4_dp) < 2e-12_dp, &
            'SI FKR eta/rho/Fprime scaling')
        call check(abs(second/kink - 4.0_dp**(1.0_dp/3)) < 2e-13_dp, &
            'SI internal kink eta/rho/Fprime scaling')
        speed = rutherford_speed(oracle_mu0*0.25_dp, -2.0_dp)
        call check(abs(speed + 0.61_dp) < 1e-15_dp, &
            'Rutherford full-width decay retains sign and SI diffusivity')
        ! Harris B_y=tanh(x/a): J''(0)/J(0)=-2/a^2, a=1 m.
        call symmetric_width(0.0164_dp, -2.0_dp, width, status)
        call check(status == 0 .and. abs(width - 0.02_dp) < 1e-15_dp, &
            'Militello-Porcelli/Escande-Ottaviani symmetric Harris saturation')
        call cylinder_saturated_width(0.0164_dp, 1.0_dp, 0.0_dp, -2.0_dp, &
            0.0_dp, 0.5_dp, [1e-5_dp, 0.1_dp], second, status)
        call check(status == 0 .and. abs(second/width - 1) < 5e-4_dp, &
            'Arcis cylindrical limit agrees within rounding of published 0.41')
        call tearing_rates(eta, rho, fp, -delta, fkr, kink, status)
        call check(status == oracle_domain, 'unstable FKR branch rejects Delta<0')
    end subroutine

    subroutine step_profile()
        type(cyl_step_t) :: model
        type(outer_result_t) :: value
        real(dp) :: q0, rs, g, exact, wall_exact, error
        integer :: i
        do i = 1, 3
            q0 = 0.6_dp + 0.2_dp*i
            model = cyl_step_t(m=2, n=1, q0=q0, r1=0.3_dp)
            rs = model%r1*sqrt(2/q0)
            g = 2/(1 - q0/2)
            ! The analytic axis branch r^m crosses the exact current jump.
            call shoot_outer(vacuum_coefficients, 2, 1.0_dp, rs, 0.0_dp, &
                1e-8_dp, 1e-11_dp, value, axis_start=model%r1, &
                axis_logderivative=(2 - g)/model%r1, wall_logderivative=-2.0_dp)
            exact = 16*(q0 - 1)/(q0 - 2)**2/rs
            error = abs(value%delta - exact)
            print '(a,f6.3,a,es12.4)', 'top-hat q0=', q0, ' Delta error [1/m]=', error
            ! The 1e-8 m cut leaves matching error beyond the ODE tolerance.
            call check(value%status == 0 .and. error < 3e-5_dp, &
                'top-hat cylinder matches independently derived power solutions')
            call check(abs(step_delta_prime(model)/rs - exact) < 2e-12_dp, &
                'existing response strip has the same no-wall convention')
            if (i /= 2) call check(value%delta*exact > 0, &
                'top-hat family changes sign across q0=1')
            call shoot_outer(vacuum_coefficients, 2, 1.0_dp, rs, 0.0_dp, &
                1e-8_dp, 1e-11_dp, value, axis_start=model%r1, &
                axis_logderivative=(2 - g)/model%r1)
            wall_exact = exact - 4*rs**3/(1 - rs**4)
            call check(value%status == 0 .and. &
                abs(value%delta - wall_exact) < 3e-5_dp, &
                'conducting wall has the analytic stabilizing correction')
        end do
    end subroutine

    subroutine vacuum_coefficients(r, p, u)
        real(dp), intent(in) :: r
        real(dp), intent(out) :: p, u
        p = 1/r; u = 4/r**2
    end subroutine

    subroutine force_free()
        type(fortnum_status_t) :: stat
        type(outer_result_t) :: result
        real(dp) :: rs, pole, exact, jm, djm, ym, dym, jw, djw, yw, dyw
        call root_brent(bessel_resonance, 0.01_dp, 1.0_dp, rs, stat, xtol=1e-14_dp)
        call check(stat%code == 0, 'force-free Bessel field has one resonance')
        ! Constant lambda gives no current-gradient logarithm.
        pole = 0
        call shoot_outer(bessel_coefficients, 1, 1.0_dp, rs, pole, &
            1e-8_dp, 1e-11_dp, result)
        call helical_bessel(rs, jm, djm, ym, dym)
        call helical_bessel(1.0_dp, jw, djw, yw, dyw)
        exact = (djm*yw - dym*jw)/(jm*yw - ym*jw) - djm/jm
        print '(a,2es18.9)', 'Newcomb Delta numerical/Bessel [1/m]: ', &
            result%delta, exact
        ! Allow finite-cut and axis-start errors as in the power-solution check.
        call check(result%status == 0 .and. abs(result%delta - exact) < 3e-5_dp, &
            'full zero-pressure Newcomb agrees with force-free Bessel solution')
    end subroutine

    pure real(dp) function bessel_resonance(r) result(f)
        real(dp), intent(in) :: r
        f = bessel_j1(1.8_dp*r)/r - bessel_j0(1.8_dp*r)
    end function

    subroutine bessel_coefficients(r, p, u)
        real(dp), intent(in) :: r
        real(dp), intent(out) :: p, u
        real(dp) :: bz, bt, dbt, dbz, f, fp, fpp, residue
        bz = bessel_j0(1.8_dp*r); bt = bessel_j1(1.8_dp*r)
        dbz = -1.8_dp*bt; dbt = 1.8_dp*bz - bt/r
        f = bt/r - bz
        fp = dbt/r - bt/r**2 - dbz
        fpp = (-1.8_dp**2*bt - dbt/r + bt/r**2)/r &
            - 2*dbt/r**2 + 2*bt/r**3 + 1.8_dp*dbt
        call newcomb_coefficients(r, 1.0_dp, -1.0_dp, bt, bz, f, fp, fpp, &
            p, u, residue)
    end subroutine

    subroutine helical_bessel(r, j, dj, y, dy)
        real(dp), intent(in) :: r
        real(dp), intent(out) :: j, dj, y, dy
        real(dp) :: kap, x, bj, dbj, by, dby
        kap = sqrt(1.8_dp**2 - 1); x = kap*r
        bj = bessel_j1(x); dbj = (bessel_j0(x) - bessel_jn(2, x))/2
        by = bessel_y1(x); dby = (bessel_y0(x) - bessel_yn(2, x))/2
        j = 1.8_dp*bj - r*kap*dbj
        y = 1.8_dp*by - r*kap*dby
        dj = 1.8_dp*kap*dbj + r*(kap**2 - 1/r**2)*bj
        dy = 1.8_dp*kap*dby + r*(kap**2 - 1/r**2)*by
    end subroutine

    subroutine convergence()
        type(wesson_t) :: model
        type(outer_result_t) :: coarse, medium, fine, cut1, cut2
        model = wesson_t(q0=1.1_dp)
        call wesson_reference(model, 1e-6_dp, 1e-5_dp, coarse)
        call wesson_reference(model, 1e-6_dp, 1e-8_dp, medium)
        call wesson_reference(model, 1e-6_dp, 1e-11_dp, fine)
        print '(a,3es18.9)', 'ODE tolerance ladder Delta [1/m]: ', &
            coarse%delta, medium%delta, fine%delta
        call check(coarse%status == 0 .and. medium%status == 0 .and. fine%status == 0, &
            'ODE tolerance ladder succeeds')
        ! 1e-5/m is below a 2e-6 relative budget for this order-one index.
        call check(abs(medium%delta - fine%delta) < 1e-5_dp .and. &
            abs(medium%delta - fine%delta) < abs(coarse%delta - fine%delta), &
            'Delta converges under adaptive ODE tolerance refinement')
        call wesson_reference(model, 1e-4_dp, 1e-11_dp, cut1)
        call wesson_reference(model, 1e-7_dp, 1e-11_dp, cut2)
        print '(a,3es18.9)', 'Frobenius cut ladder Delta [1/m]: ', &
            cut1%delta, fine%delta, cut2%delta
        call check(cut1%status == 0 .and. cut2%status == 0, 'cut ladder succeeds')
        call check(abs(fine%delta - cut2%delta) < 1e-5_dp .and. &
            abs(fine%delta - cut2%delta) < abs(cut1%delta - cut2%delta), &
            'Frobenius finite parts converge independently of ODE tolerance')
        model%m = 1
        call wesson_reference(model, 1e-6_dp, 1e-10_dp, fine)
        call check(fine%status == oracle_no_resonance, 'q0=1.1 has no 1/1 resonance')
        model%q0 = 0.9_dp
        call wesson_reference(model, 1e-7_dp, 1e-11_dp, fine)
        call check(fine%status == oracle_singular, &
            'reduced 1/1 kink has zero outer flux, no finite Delta')
        call logarithmic_solution()
    end subroutine

    subroutine logarithmic_solution()
        type(outer_result_t) :: coarse, fine
        real(dp) :: x, large, derivative, exact
        x = 0.1_dp - 0.4_dp
        large = 1 + 0.7_dp*x*log(abs(x))
        derivative = 0.7_dp*(log(abs(x)) + 1)
        ! Exact independent branches: 1+A x log|x| and x, A=0.7/m.
        ! The right combination vanishes at wall=1 m; left has zero finite slope.
        x = 1.0_dp - 0.4_dp
        exact = -(1 + 0.7_dp*x*log(x))/x
        call shoot_outer(logarithmic_coefficients, 1, 1.0_dp, 0.4_dp, 0.7_dp, &
            1e-3_dp, 1e-11_dp, coarse, axis_start=0.1_dp, &
            axis_logderivative=derivative/large)
        call shoot_outer(logarithmic_coefficients, 1, 1.0_dp, 0.4_dp, 0.7_dp, &
            1e-5_dp, 1e-11_dp, fine, axis_start=0.1_dp, &
            axis_logderivative=derivative/large)
        call check(coarse%status == 0 .and. fine%status == 0, &
            'manufactured logarithmic branches integrate from both sides')
        call check(abs(fine%delta - exact) < 2e-7_dp, &
            'Frobenius finite parts match independent exact logarithmic solution')
    end subroutine

    subroutine logarithmic_coefficients(r, p, u)
        real(dp), intent(in) :: r
        real(dp), intent(out) :: p, u
        p = 0.7_dp/(1 - 0.7_dp*(r - 0.4_dp))
        u = p/(r - 0.4_dp)
    end subroutine

    subroutine saturation()
        type(wesson_t) :: model
        type(outer_result_t) :: value
        real(dp) :: rs, t, f, fp, fpp, jp, jpp, js, alpha, beta, asymmetry
        real(dp) :: width, scaled, positive_delta, delta
        integer :: i, status
        positive_delta = 0
        ! Hastie et al. (IAEA 2004 PD/1-1, section 5): nu=2 marginal near rs/a=0.88.
        do i = 1, 2
            rs = 0.865_dp + 0.01_dp*i
            t = rs**2
            model%q0 = 2*(1 - t + t*t/3); model%qa = 3*model%q0
            call wesson_reference(model, 1e-5_dp, 1e-11_dp, value)
            call check(value%status == 0, 'published nu=2 family resolves')
            if (i == 1) positive_delta = value%delta
        end do
        call check(positive_delta > 0 .and. value%delta < 0, &
            'published nu=2 marginal at rs/a about 0.88 is bracketed')
        rs = 0.8788_dp; t = rs**2
        model%q0 = 2*(1 - t + t*t/3); model%qa = 3*model%q0
        call wesson_reference(model, 1e-5_dp, 1e-11_dp, value)
        call check(value%status == 0, 'near-marginal nu=2 profile resolves')
        call model%evaluate(value%rs, f, fp, fpp, jp, jpp, js)
        delta = value%delta; alpha = jp/js; beta = jpp/js
        asymmetry = value%slope_minus + value%slope_plus
        call cylinder_saturated_width(delta, 1.0_dp, alpha, beta, asymmetry, &
            rs, [1e-7_dp, 0.02_dp], width, status)
        print '(a,3es16.7)', 'nu=2 cylinder rs[m], Delta[1/m], Wsat[m]: ', &
            value%rs, delta, width
        call check(status == 0 .and. width > 0 .and. width < 0.02_dp*(1 - rs), &
            'saturation reference stays in the small-island regime')
        call cylinder_saturated_width(delta/2, 2.0_dp, alpha/2, beta/4, &
            asymmetry/2, 2*rs, [2e-7_dp, 0.04_dp], scaled, status)
        call check(status == 0 .and. abs(scaled/width - 2) < 1e-8_dp, &
            'cylindrical saturation is invariant under a change of length unit')
    end subroutine

    subroutine table()
        type(wesson_t) :: model
        type(outer_result_t) :: result
        real(dp) :: f, fp, fpp, jp, jpp, jstar, fkr, kink
        integer :: iq, m, status, expected_status
        print '(a)', 'Wesson reduced cylinder: a=1 m, Bz=1 T, R0=3 m, qa=3'
        print '(a)', 'rho=1e-6 kg/m^3, eta=1e-8 Ohm m; wall at a; phase=m theta-n z/R0'
        print '(a)', 'q0 nu m n rs[m] Delta[1/m] status; rates in 1/s, dW/dt in m/s'
        do iq = 1, 2
            model%q0 = 0.7_dp + 0.2_dp*iq
            do m = 2, 1, -1
                model%m = m
                call wesson_reference(model, 1e-5_dp, 1e-11_dp, result)
                print '(2f9.5,2i3,2es18.9,i5)', model%q0, model%qa/model%q0 - 1, m, &
                    model%n, result%rs, result%delta, result%status
                expected_status = oracle_ok
                if (m == 1 .and. iq == 1) expected_status = oracle_singular
                if (m == 1 .and. iq == 2) expected_status = oracle_no_resonance
                call check(result%status == expected_status, &
                    'reference profile resolved')
                if (result%status == oracle_no_resonance) cycle
                call model%evaluate(result%rs, f, fp, fpp, jp, jpp, jstar)
                call tearing_rates(1e-8_dp, 1e-6_dp, fp/3, max(0.0_dp, result%delta), &
                    fkr, kink, status)
                if (m == 1) then
                    print '(a,es18.9)', 'marginal-ideal kink asymptote only: ', kink
                else if (result%delta > 0) then
                    print '(a,es18.9)', 'constant-psi FKR asymptote: ', fkr
                else
                    print '(a)', 'Delta<=0: no unstable FKR branch'
                end if
                if (result%status == oracle_ok) print '(a,es18.9)', &
                    'Rutherford full-width speed: ', &
                    rutherford_speed(1e-8_dp, result%delta)
            end do
        end do
        print '(a)', 'status: 0 finite Delta; -2 no resonance; '// &
            '-3 singular kink outer flux'
        print '(a)', 'Candidate references; ODE and Frobenius errors '// &
            'are not certified bounds.'
        call saturation()
    end subroutine
end program test_rm_oracles
