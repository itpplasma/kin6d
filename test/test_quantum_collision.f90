! Independent behavioural oracles for the QG6 binary Coulomb collision kernel.
! Oracles: digamma series, closed-form Rutherford/Mott angular integrals,
! Legendre Q_l by quad-precision Miller recurrence (first Born Yukawa phases),
! direct angular quadrature of the Born amplitude (transfer/viscosity sums),
! Coulomb phase identity sigma_{l+1}-sigma_l = atan(eta/(l+1)) at weak
! screening, classical Yukawa orbits by a test-local RK4 integration, the
! classical limit eta >> 1 and central finite differences.
program test_quantum_collision
    use, intrinsic :: iso_fortran_env, only: dp => real64, qp => real128
    use kin6d_quantum_collision_kernel
    use fortnum_quadrature, only: gauss_legendre_ab
    implicit none
    integer :: failures = 0
    real(dp), parameter :: pi = acos(-1.0_dp)
    call test_digamma_cut_log()
    call test_mott_integrals()
    call test_free_and_born()
    call test_coulomb_limit()
    call test_classical_limit_and_derivative()
    call test_classical_yukawa_orbit()
    if (failures > 0) then
        print '(i0,a)', failures, ' quantum collision kernel checks failed'
        error stop 1
    end if
    print '(a)', 'all quantum collision kernel checks passed'
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

    ! Re psi(1+i eta) = -gamma + eta^2 sum 1/(n(n^2+eta^2)); the long-sum branch
    ! of quantum_cut_log against the direct sum; d lnL/d eta by differences.
    subroutine test_digamma_cut_log()
        complex(dp) :: p, dp1
        real(qp) :: s, eta
        real(qp), parameter :: etas(4) = [0.05_qp, 0.7_qp, 3.0_qp, 25.0_qp]
        real(dp) :: err, lnl, dl, lp, lm, dum, direct, kappa
        integer :: n, i
        err = 0
        do i = 1, 4
            eta = etas(i)
            s = 0
            do n = 2000000, 1, -1  ! real(n, qp)**2: n**2 overflows default integers
                s = s + 1/(real(n, qp)*(real(n, qp)**2 + eta**2))
            end do
            s = -0.57721566490153286060651209_qp + eta**2*(s + 1/(2*real(2000000, qp)**2))
            call complex_digamma(cmplx(1.0_dp, real(eta, dp), dp), p, dp1)
            err = max(err, abs(real(p) - real(s, dp)))
        end do
        call check('Re psi(1+i eta) vs series', err < 1.0e-12_dp, err)
        ! kappa with 6000 full rings (digamma branch) vs the direct sum.
        kappa = sqrt(6000.0_dp*6001.0_dp) + 1.0e-9_dp
        call quantum_cut_log(0.8_dp, kappa, lnl, dl)
        direct = 0
        do n = 6000, 1, -1
            direct = direct + n/(n**2 + 0.64_dp)
        end do
        call check('ring sum, digamma branch vs direct', abs(lnl - direct) < 1.0e-11_dp, &
                   abs(lnl - direct))
        call quantum_cut_log(0.8_dp + 1.0e-5_dp, kappa, lp, dum)
        call quantum_cut_log(0.8_dp - 1.0e-5_dp, kappa, lm, dum)
        call check('d lnL_q/d eta (trigamma) vs FD', &
                   abs(dl - (lp - lm)/2.0e-5_dp) < 1.0e-8_dp*abs(dl), &
                   abs(dl - (lp - lm)/2.0e-5_dp))
        call quantum_cut_log(0.8_dp, 37.3_dp, lnl, dl)
        call quantum_cut_log(0.8_dp + 1.0e-5_dp, 37.3_dp, lp, dum)
        call quantum_cut_log(0.8_dp - 1.0e-5_dp, 37.3_dp, lm, dum)
        call check('d lnL_q/d eta (direct) vs FD', &
                   abs(dl - (lp - lm)/2.0e-5_dp) < 1.0e-8_dp*abs(dl), &
                   abs(dl - (lp - lm)/2.0e-5_dp))
        ! Limits: eta -> 0 gives ln kappa + gamma; eta >> 1 the classical log.
        call quantum_cut_log(1.0e-6_dp, 1.0e6_dp, lnl, dl)
        call check('quantum cut log, eta->0: ln kappa + gamma', &
                   abs(lnl - log(1.0e6_dp) - 0.5772156649_dp) < 1.0e-6_dp, &
                   lnl - log(1.0e6_dp) - 0.5772156649_dp)
        call quantum_cut_log(50.0_dp, 1.0e6_dp, lnl, dl)
        call classical_cut_log(50.0_dp, 1.0e6_dp, lp, dum)
        call check('quantum cut log, eta=50 vs classical', abs(lnl - lp) < 1.0e-4_dp, &
                   lnl - lp)
    end subroutine test_digamma_cut_log

    ! Closed-form angular integrals: distinguishable Rutherford with the
    ! classical angle cut; identical classical folding (1-2 ln2)/2; Mott
    ! unpolarized eta -> 0 adds -ln2/2; eta >> 1 the interference averages out.
    subroutine test_mott_integrals()
        real(dp) :: lnl, lc, d
        call classical_cut_log(2.0_dp, 1.0e3_dp, lc, d)
        lnl = mott_transfer_log(2.0_dp, 1.0e3_dp, SPIN_DISTINGUISHABLE)
        call check('Rutherford quadrature vs cut log', abs(lnl - lc) < 1.0e-10_dp, lnl - lc)
        call classical_cut_log(0.01_dp, 1.0e3_dp, lc, d)
        lnl = mott_transfer_log(0.01_dp, 1.0e3_dp, SPIN_IDENTICAL_CLASSICAL) - lc
        call check('identical classical folding (1-2ln2)/2', &
                   abs(lnl - (1 - 2*log(2.0_dp))/2) < 1.0e-10_dp, lnl - (1 - 2*log(2.0_dp))/2)
        call classical_cut_log(1.0e-4_dp, 1.0e3_dp, lc, d)
        lnl = mott_transfer_log(1.0e-4_dp, 1.0e3_dp, SPIN_UNPOLARIZED) - lc
        call check('Mott unpolarized eta->0: 1/2 - 3/2 ln2', &
                   abs(lnl - (0.5_dp - 1.5_dp*log(2.0_dp))) < 1.0e-5_dp, &
                   lnl - (0.5_dp - 1.5_dp*log(2.0_dp)))
        lnl = mott_transfer_log(1.0e-4_dp, 1.0e3_dp, SPIN_SINGLET) + &
              3*mott_transfer_log(1.0e-4_dp, 1.0e3_dp, SPIN_TRIPLET) - &
              4*mott_transfer_log(1.0e-4_dp, 1.0e3_dp, SPIN_UNPOLARIZED)
        call check('singlet + 3 triplet = 4 unpolarized', abs(lnl) < 1.0e-9_dp, lnl)
        lnl = mott_transfer_log(20.0_dp, 1.0e4_dp, SPIN_UNPOLARIZED) - &
              mott_transfer_log(20.0_dp, 1.0e4_dp, SPIN_IDENTICAL_CLASSICAL)
        call check('Mott eta=20: exchange term averages out', abs(lnl) < 2.0e-3_dp, lnl)
    end subroutine test_mott_integrals

    ! Legendre Q_l(z), z > 1, l = 0..lmax: quad-precision Miller recurrence
    ! normalized by Q_0 = artanh(1/z) (test-local oracle).
    subroutine legendre_q_miller(lmax, z, q)
        integer, intent(in) :: lmax
        real(qp), intent(in) :: z
        real(qp), intent(out) :: q(0:lmax)
        real(qp) :: a, b, c, scale
        integer :: l, lstart
        lstart = lmax + 60 + ceiling(80/acosh(z))
        a = 0; b = 1.0e-30_qp
        do l = lstart, 1, -1
            c = ((2*l + 1)*z*b - (l + 1)*a)/l    ! Q_{l-1}
            a = b; b = c
            if (l - 1 <= lmax) q(l - 1) = c
            if (abs(b) > 1.0e300_qp) then
                scale = 1/abs(b); a = a*scale; b = b*scale
                q(max(0, l - 1):lmax) = q(max(0, l - 1):lmax)*scale
            end if
        end do
        q = q*(log((z + 1)/(z - 1))/2)/q(0)
    end subroutine legendre_q_miller

    ! Free phases vanish; first Born phases -eta Q_l(1 + 1/(2 kappa^2)) and the
    ! closed-form Born transfer log; transfer and viscosity sums of Born phases
    ! against direct angular quadrature of the Born cross section.
    subroutine test_free_and_born()
        type(yukawa_plan_t) :: plan
        type(yukawa_config_t) :: cfg
        real(dp), allocatable :: d(:), dd(:), qd(:)
        real(qp), allocatable :: q(:)
        real(dp) :: kappa, eta, err, lnl, dl, ang, vis, x(16), w(16), th, qq, f2, eps
        integer :: l, ip
        logical :: ok
        kappa = 10
        call yukawa_plan(plan, 1.0_dp, kappa, cfg, ok)
        call yukawa_phases(plan, 0.0_dp, d, dd, ok)
        call check('free phases (DVR dispersion indicator)', &
                   ok .and. maxval(abs(sin(d))) < 1.0e-7_dp, maxval(abs(sin(d))))
        eta = 1.0e-4_dp
        call yukawa_plan(plan, eta, kappa, cfg, ok)
        call yukawa_phases(plan, eta, d, dd, ok)
        allocate (q(0:plan%lmax), qd(0:plan%lmax))
        call legendre_q_miller(plan%lmax, 1 + 1/(2*real(kappa, qp)**2), q)
        qd = real(q, dp)
        err = 0
        do l = 0, 20
            err = max(err, abs(-d(l)/(eta*qd(l)) - 1))
        end do
        call check('Born phases -eta Q_l(1+1/(2 kappa^2)), l<=20', ok .and. err < 1.0e-3_dp, err)
        call yukawa_transfer_log(plan, eta, lnl, dl, ok)
        call check('Yukawa lnL vs closed-form Born, eta=1e-4', &
                   abs(lnl - born_yukawa_log(kappa)) < 1.0e-4_dp, lnl - born_yukawa_log(kappa))
        eps = 1.0e-6_dp
        lnl = transfer_from_phases(eps*qd)/eps**2
        call check('Born phase sum = closed-form Born log', &
                   abs(lnl/born_yukawa_log(kappa) - 1) < 1.0e-9_dp, &
                   lnl/born_yukawa_log(kappa) - 1)
        ! Normalized sigma^(2) = 2 int sin^2 sin dtheta/(q^2 + 1/kappa^2)^2, q = 2 sin(theta/2).
        ang = 0
        do ip = 1, 400
            call gauss_legendre_ab(16, pi*(ip - 1)/400, pi*ip/400, x, w)
            do l = 1, 16
                th = x(l); qq = 2*sin(th/2)
                f2 = 1/(qq**2 + 1/kappa**2)**2
                ang = ang + w(l)*2*sin(th)**3*f2
            end do
        end do
        vis = viscosity_from_phases(eps*qd)/eps**2
        call check('Born viscosity sum vs angular quadrature', abs(vis/ang - 1) < 1.0e-9_dp, &
                   vis/ang - 1)
    end subroutine test_free_and_born

    ! Weak screening: low-l Yukawa phase differences approach the Coulomb
    ! identity sigma_{l+1} - sigma_l = atan(eta/(l+1)) (deviation O(eta^2/kappa)).
    subroutine test_coulomb_limit()
        type(yukawa_plan_t) :: plan
        type(yukawa_config_t) :: cfg
        real(dp), allocatable :: d(:), dd(:)
        real(dp) :: eta, kappa, err
        integer :: l
        logical :: ok
        eta = 0.3_dp; kappa = 100
        call yukawa_plan(plan, eta, kappa, cfg, ok)
        call yukawa_phases(plan, eta, d, dd, ok)
        err = 0
        do l = 0, 4
            err = max(err, abs(sin(d(l + 1) - d(l) - atan(eta/(l + 1)))))
        end do
        print '(a,es10.3)', '  Coulomb-limit deviation x kappa/eta^2: ', err*kappa/eta**2
        call check('Yukawa -> Coulomb phase differences, l<=4', ok .and. err < 2*eta**2/kappa, &
                   err)
    end subroutine test_coulomb_limit

    ! eta >> 1: quantum Yukawa transfer approaches the classical Yukawa orbit
    ! integral; implicit d lnL/d eta against central differences on the frozen
    ! plan; resolution indicator (10 vs 14 nodes per wavelength).
    subroutine test_classical_limit_and_derivative()
        type(yukawa_plan_t) :: plan
        type(yukawa_config_t) :: cfg, fine
        real(dp) :: lnl, dl, lp, lm, dum, lcl, lf
        logical :: ok
        call yukawa_plan(plan, 5.0_dp, 30.0_dp, cfg, ok)
        call yukawa_transfer_log(plan, 5.0_dp, lnl, dl, ok)
        lcl = classical_yukawa_log(6.0_dp, 200, 16)
        call check('eta=5: quantum Yukawa lnL vs classical orbit integral', &
                   ok .and. abs(lnl/lcl - 1) < 1.0e-2_dp, lnl/lcl - 1)
        call yukawa_plan(plan, 0.7_dp, 30.0_dp, cfg, ok)
        call yukawa_transfer_log(plan, 0.7_dp, lnl, dl, ok)
        call yukawa_transfer_log(plan, 0.7_dp + 1.0e-5_dp, lp, dum, ok)
        call yukawa_transfer_log(plan, 0.7_dp - 1.0e-5_dp, lm, dum, ok)
        call check('implicit d lnL_Y/d eta vs FD', &
                   abs(dl - (lp - lm)/2.0e-5_dp) < 1.0e-6_dp*abs(dl), &
                   abs(dl - (lp - lm)/2.0e-5_dp)/abs(dl))
        fine = cfg; fine%nodes_per_wavelength = 14; fine%tail = 34; fine%barrier = 44
        call yukawa_plan(plan, 0.7_dp, 30.0_dp, fine, ok)
        call yukawa_transfer_log(plan, 0.7_dp, lf, dum, ok)
        call check('Yukawa lnL resolution indicator', abs(lf - lnl) < 1.0e-7_dp, lf - lnl)
    end subroutine test_classical_limit_and_derivative

    ! Classical Yukawa deflection quadrature against a test-local RK4 orbit
    ! (b90 = 1 units: mu = v_inf = g = 1) and the Rutherford limit.
    subroutine test_classical_yukawa_orbit()
        real(dp) :: lam, b, y(6), k1(6), k2(6), k3(6), k4(6), dt, th, err, r
        real(dp), parameter :: bs(3) = [0.5_dp, 2.0_dp, 10.0_dp]
        integer :: i, it
        lam = 20; err = 0
        do i = 1, 3
            b = bs(i)
            y = [b, 0.0_dp, -40*lam, 0.0_dp, 0.0_dp, 1.0_dp]
            do it = 1, 10000000
                r = norm2(y(1:3))
                dt = 2.0e-3_dp*max(r, 0.5_dp)
                k1 = rhs(y, lam); k2 = rhs(y + dt/2*k1, lam)
                k3 = rhs(y + dt/2*k2, lam); k4 = rhs(y + dt*k3, lam)
                y = y + dt/6*(k1 + 2*k2 + 2*k3 + k4)
                if (norm2(y(1:3)) > 40*lam .and. dot_product(y(1:3), y(4:6)) > 0) exit
            end do
            th = atan2(norm2(y(4:5)), y(6))
            err = max(err, abs(classical_yukawa_deflection(b, lam, 200)/th - 1))
        end do
        call check('classical Yukawa deflection vs RK4 orbit', err < 1.0e-6_dp, err)
        th = classical_yukawa_deflection(3.0_dp, 1.0e9_dp, 200)
        call check('classical Yukawa deflection, Rutherford limit', &
                   abs(th/(2*atan(1/3.0_dp)) - 1) < 1.0e-8_dp, th/(2*atan(1/3.0_dp)) - 1)
    end subroutine test_classical_yukawa_orbit

    function rhs(s, lam) result(f)
        real(dp), intent(in) :: s(6), lam
        real(dp) :: f(6), rr
        rr = norm2(s(1:3))
        f(1:3) = s(4:6)
        f(4:6) = exp(-rr/lam)*(1 + rr/lam)/rr**3*s(1:3)
    end function rhs
end program test_quantum_collision
