! Independent behavioural oracles for the Q1/Q1b partial-wave Coulomb encounter.
! Oracles: GLL/Legendre exactness, exact hydrogen levels, Coulomb phases
! arg Gamma(l+1+i eta) (quad-precision Stirling) from box quantization,
! Coulomb F_l power series (quad precision), analytic free Gaussian, exchange
! identities, Rutherford/Mott plane-wave formulas, the classical Wigner
! ensemble, finite differences and the dot-product identity.
program test_coulomb_encounter
    use, intrinsic :: iso_fortran_env, only: dp => real64, qp => real128, i64 => int64
    use kin6d_radial_dvr
    use kin6d_coulomb_encounter
    use kin6d_coulomb_kepler
    use kin6d_coulomb_encounter_leaves, only: ce_mott, ce_rutherford, ce_free_width
    use fortnum_quadrature, only: gauss_legendre, gauss_legendre_ab, gauss_lobatto_legendre
    use fortnum_special_legendre, only: legendre_p_normalized_table
    implicit none
    integer :: failures = 0
    real(dp), parameter :: pi = acos(-1.0_dp)
    call test_gll_legendre()
    call test_hydrogen_spectrum()
    call test_coulomb_phases()
    call test_coulomb_function()
    call test_free_gaussian()
    call test_exchange_sectors()
    call test_mott_rutherford()
    call test_classical_ensemble()
    call test_derivatives()
    call test_kepler()
    if (failures > 0) then
        print '(i0,a)', failures, ' coulomb encounter checks failed'
        error stop 1
    end if
    print '(a)', 'all coulomb encounter checks passed'
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

    subroutine test_gll_legendre()
        real(dp) :: x(10), w(10), xq(40), wq(40), p(0:20, 0:5, 40), gram, err
        integer :: k, l, l2, m, q
        call gauss_lobatto_legendre(10, x, w)
        err = 0
        do k = 0, 2*10 - 3
            err = max(err, abs(sum(w*x**k) - merge(2.0_dp/(k + 1), 0.0_dp, mod(k, 2) == 0)))
        end do
        call check('GLL exactness degree 2n-3', err < 1.0e-14_dp, err)
        call gauss_legendre(40, xq, wq)
        do q = 1, 40
            call legendre_p_normalized_table(20, 5, xq(q), p(:, :, q))
        end do
        err = 0
        do m = 0, 5
            do l = m, 20
                do l2 = m, 20
                    gram = sum(wq*p(l, m, :)*p(l2, m, :))
                    if (l == l2) gram = gram - 1
                    err = max(err, abs(gram))
                end do
            end do
        end do
        call check('normalized Legendre orthonormality', err < 1.0e-13_dp, err)
    end subroutine test_gll_legendre

    ! Exact Coulomb bound levels E_n = -mu g^2/(2 hbar^2 n^2), mu = hbar = 1, g = -1.
    subroutine test_hydrogen_spectrum()
        type(radial_mesh_t) :: mesh
        real(dp), allocatable :: e(:), v(:, :), res(:)
        real(dp) :: orth, err, exact, worst_res
        integer :: l, j
        logical :: ok
        call build_radial_mesh(mesh, 1500.0_dp, 11, &
                               mesh_spacing(sqrt(2/(1500.0_dp/8)), 11, 6.0_dp), &
                               0.5_dp, 1.08_dp, ok)
        call check('graded mesh built', ok, real(mesh%n, dp))
        err = 0; worst_res = 0
        do l = 0, 3
            call radial_eigen_window(mesh, l, -1.0_dp, 1.0_dp, 1.0_dp, -1.0_dp/(2*24.0_dp**2), &
                                     e, v, res, orth, ok)
            if (.not. ok) err = huge(1.0_dp)
            do j = 1, 20 - l
                exact = -0.5_dp/real(j + l, dp)**2
                err = max(err, abs(e(j) - exact)/abs(exact))
            end do
            worst_res = max(worst_res, maxval(res), orth)
        end do
        call check('hydrogen levels n<=20, l<=3 relative error', err < 1.0e-10_dp, err)
        call check('eigen residual/orthogonality', worst_res < 1.0e-11_dp, worst_res)
    end subroutine test_hydrogen_spectrum

    ! sigma_L = arg Gamma(L+1+i eta): shifted Stirling series in binary128.
    real(qp) function coulomb_sigma(l, eta) result(s)
        integer, intent(in) :: l
        real(qp), intent(in) :: eta
        real(qp), parameter :: b(7) = [1/6.0_qp, -1/30.0_qp, 1/42.0_qp, -1/30.0_qp, &
                                       5/66.0_qp, -691/2730.0_qp, 7/6.0_qp]
        complex(qp) :: z, w, lg
        integer :: j, k
        z = cmplx(1.0_qp, eta, qp)
        w = z + 40
        lg = (w - 0.5_qp)*log(w) - w + 0.5_qp*log(2*acos(-1.0_qp))
        do k = 1, 7
            lg = lg + b(k)/(real(2*k*(2*k - 1), qp)*w**(2*k - 1))
        end do
        s = aimag(lg)
        do j = 0, 39
            s = s - aimag(log(z + j))
        end do
        do j = 1, l
            s = s + atan(eta/real(j, qp))
        end do
    end function coulomb_sigma

    ! Box quantization u(R)=0 <=> F_L(eta,kR)=0. With the asymptotic expansion
    ! F_L = g cos(theta) + f sin(theta) (A&S 14.5), theta = rho - eta ln(2 rho)
    ! - L pi/2 + sigma_L, each DVR eigenvalue yields sigma_L modulo pi.
    real(qp) function box_sigma(l, eta, rho) result(s)
        integer, intent(in) :: l
        real(qp), intent(in) :: eta, rho
        real(qp) :: f, g, fk, gk, fn, gn, ak, bk, theta
        integer :: k
        f = 1; g = 0; fk = 1; gk = 0
        do k = 0, 200
            ak = real(2*k + 1, qp)*eta/(real(2*k + 2, qp)*rho)
            bk = (real(l*(l + 1) - k*(k + 1), qp) + eta**2)/(real(2*k + 2, qp)*rho)
            fn = ak*fk - bk*gk; gn = ak*gk + bk*fk
            if (abs(fn) + abs(gn) > abs(fk) + abs(gk) .and. k > 2) exit
            fk = fn; gk = gn; f = f + fk; g = g + gk
            if (abs(fk) + abs(gk) < 1.0e-32_qp) exit
        end do
        theta = -atan(g/f)
        s = theta - (rho - eta*log(2*rho) - real(l, qp)*acos(-1.0_qp)/2)
    end function box_sigma

    subroutine test_coulomb_phases()
        type(radial_mesh_t) :: mesh
        real(dp), allocatable :: e(:), v(:, :), res(:)
        real(dp) :: orth, g, err, h0, growth
        real(qp) :: k, eta, diff
        integer :: l, j, is, count
        logical :: ok
        real(dp), parameter :: couplings(3) = [1.0_dp, -1.0_dp, 0.03_dp]
        character(48), parameter :: names(3) = [character(48) :: &
            'Coulomb phase sigma_l, repulsive', 'Coulomb phase sigma_l, attractive', &
            'Coulomb phase sigma_l, repulsive eta ~ 0.03']
        ! Repulsive and attractive at eta ~ 1, and weak repulsive coupling at
        ! eta ~ 0.03 (the hot-electron encounter movie).
        do is = 1, 3
            g = couplings(is)
            h0 = mesh_spacing(1.6_dp, 12, 8.0_dp); growth = 1
            if (g < 0) then
                h0 = 0.5_dp; growth = 1.15_dp
            end if
            call build_radial_mesh(mesh, 60.0_dp, 12, mesh_spacing(1.6_dp, 12, 8.0_dp), &
                                   h0, growth, ok)
            err = 0; count = 0
            do l = 0, 4
                call radial_eigen_window(mesh, l, g, 1.0_dp, 1.0_dp, 0.8_dp, e, v, res, &
                                         orth, ok)
                do j = 1, size(e)
                    if (e(j) <= 0) cycle
                    k = sqrt(2*real(e(j), qp))
                    if (k < 0.8_qp .or. k > 1.2_qp) cycle
                    eta = real(g, qp)/k
                    diff = box_sigma(l, eta, k*60) - coulomb_sigma(l, eta)
                    diff = diff - acos(-1.0_qp)*nint(diff/acos(-1.0_qp))
                    err = max(err, real(abs(diff), dp)); count = count + 1
                end do
            end do
            call check(trim(names(is)), err < 1.0e-8_dp .and. count > 20, err)
        end do
    end subroutine test_coulomb_phases

    ! Regular Coulomb function shape: binary128 power series (A&S 14.1.3-14.1.7)
    ! against the DVR box eigenvector on r <= R/2, after a least-squares amplitude.
    subroutine test_coulomb_function()
        type(radial_mesh_t) :: mesh
        real(dp), allocatable :: e(:), v(:, :), res(:)
        real(dp) :: orth, err, g
        real(qp) :: k, eta, f(1000), u(1000), amp
        integer :: l, j, i, is, n
        logical :: ok
        err = 0
        do is = 1, 2
            g = merge(1.0_dp, -1.0_dp, is == 1)
            call build_radial_mesh(mesh, 60.0_dp, 12, mesh_spacing(1.6_dp, 12, 8.0_dp), &
                                   merge(mesh_spacing(1.6_dp, 12, 8.0_dp), 0.5_dp, is == 1), &
                                   merge(1.0_dp, 1.15_dp, is == 1), ok)
            do l = 0, 3
                call radial_eigen_window(mesh, l, g, 1.0_dp, 1.0_dp, 0.6_dp, e, v, res, orth, ok)
                j = size(e)
                k = sqrt(2*real(e(j), qp)); eta = real(g, qp)/k
                n = 0
                do i = 1, mesh%n
                    if (mesh%r(i) > 30) exit
                    n = n + 1
                    u(n) = real(v(i, j)/sqrt(mesh%w(i)), qp)
                    f(n) = coulomb_series(l, eta, k*real(mesh%r(i), qp))
                end do
                amp = sum(u(1:n)*f(1:n))/sum(f(1:n)**2)
                err = max(err, real(sqrt(sum((u(1:n) - amp*f(1:n))**2)/sum(u(1:n)**2)), dp))
            end do
        end do
        call check('Coulomb F_l shape vs binary128 series', err < 1.0e-8_dp, err)
    end subroutine test_coulomb_function

    real(qp) function coulomb_series(l, eta, rho) result(f)
        integer, intent(in) :: l
        real(qp), intent(in) :: eta, rho
        real(qp) :: am2, am1, ak, term, big
        integer :: kk
        am2 = 1; am1 = eta/real(l + 1, qp)
        f = am2 + am1*rho; big = abs(f)
        term = 1
        do kk = l + 3, 600
            ak = (2*eta*am1 - am2)/real((kk + l)*(kk - l - 1), qp)
            term = ak*rho**(kk - l - 1)
            f = f + term; big = max(big, abs(term))
            am2 = am1; am1 = ak
            if (abs(term) < 1.0e-34_qp*big .and. kk > l + 3 + int(2*rho)) exit
        end do
        f = f*rho**(l + 1)
    end function coulomb_series

    ! g = 0, b /= 0: analytic free 3D Gaussian (all l and m participate).
    subroutine test_free_gaussian()
        type(encounter_config_t) :: cfg
        type(encounter_t) :: enc
        character(:), allocatable :: msg
        complex(dp), allocatable :: a(:, :, :)
        integer, parameter :: np = 150
        real(dp) :: x(np), y(np), z(np), c(3), t, tau, sig, err, n0, e0, dn, de, q
        complex(dp) :: psi(np), ref(np)
        integer :: p, it
        logical :: ok
        cfg%g = 0; cfg%b = 4; cfg%sigma = 3; cfg%k = 1; cfg%z0 = 12; cfg%t_final = 8
        call encounter_init(enc, cfg, ok, msg)
        call check('free encounter setup', ok, real(enc%lmax, dp))
        call check('free: m sectors present (b /= 0)', enc%mmax >= 3, real(enc%mmax, dp))
        call check('free: initial reconstruction rms', enc%ledger%reconstruction_rms < 2.0e-6_dp, &
                   enc%ledger%reconstruction_rms)
        call check('free: projection norm defect', abs(enc%ledger%projection_defect) < 1.0e-9_dp, &
                   enc%ledger%projection_defect)
        call check('free: l and m edge weights', max(enc%ledger%l_edge2, enc%ledger%m_edge2) &
                   < 1.0e-20_dp, max(enc%ledger%l_edge2, enc%ledger%m_edge2))
        allocate (a(enc%mesh%n, 0:enc%mmax, 0:enc%lmax))
        err = 0; dn = 0; de = 0
        do it = 0, 2
            t = 4.0_dp*it
            call encounter_state(enc, t, a)
            if (it == 0) then
                n0 = encounter_norm(enc, a, SECTOR_PRODUCT)
                e0 = encounter_energy(enc, a, SECTOR_PRODUCT)
            end if
            dn = max(dn, abs(encounter_norm(enc, a, SECTOR_PRODUCT) - n0))
            de = max(de, abs(encounter_energy(enc, a, SECTOR_PRODUCT) - e0)/e0)
            tau = cfg%hbar*t/(2*enc%mu*cfg%sigma**2)
            c = [cfg%b, 0.0_dp, -cfg%z0 + cfg%hbar*cfg%k*t/enc%mu]
            sig = cfg%sigma*sqrt(1 + tau**2)
            do p = 1, np
                q = real(p, dp)
                x(p) = c(1) + 2*sig*(2*modulo(q*0.7548776662_dp, 1.0_dp) - 1)
                y(p) = c(2) + 2*sig*(2*modulo(q*0.5698402910_dp, 1.0_dp) - 1)
                z(p) = c(3) + 2*sig*(2*modulo(q*0.4142135624_dp, 1.0_dp) - 1)
                ref(p) = (2*pi*cfg%sigma**2)**(-0.75_dp)*(1 + (0, 1)*tau)**(-1.5_dp)* &
                         exp(-((x(p) - c(1))**2 + (y(p) - c(2))**2 + (z(p) - c(3))**2)/ &
                             (4*cfg%sigma**2*(1 + (0, 1)*tau)) + (0, 1)*cfg%k*z(p) - &
                             (0, 1)*cfg%hbar*cfg%k**2*t/(2*enc%mu))
            end do
            call encounter_points(enc, a, x, y, z, psi)
            err = max(err, sqrt(sum(abs(psi - ref)**2)/sum(abs(ref)**2)))
        end do
        call check('free Gaussian vs analytic (t <= 8)', err < 2.0e-6_dp, err)
        call check('norm conservation', dn < 1.0e-12_dp, dn)
        call check('energy conservation (independent banded action)', de < 1.0e-12_dp, de)
    end subroutine test_free_gaussian

    ! Exchange: Y_lm(-r) = (-1)^l Y_lm(r); singlet = even l, triplet = odd l.
    subroutine test_exchange_sectors()
        type(encounter_config_t) :: cfg
        type(encounter_t) :: enc
        character(:), allocatable :: msg
        complex(dp), allocatable :: a(:, :, :), ae(:, :, :), ao(:, :, :)
        integer, parameter :: nb = 12
        real(dp) :: edges(0:nb), pp(nb), ps(nb), pt(nb), ns, nt, n, err, x(3), y(3), z(3)
        real(dp) :: redges(0:nb), pr(nb)
        complex(dp) :: psp(3), psm(3), ptp(3), ptm(3)
        integer :: ib, l
        logical :: ok
        cfg%g = 1; cfg%k = 0.5_dp; cfg%sigma = 8; cfg%z0 = 40; cfg%b = 3; cfg%t_final = 100
        call encounter_init(enc, cfg, ok, msg)
        call check('Coulomb encounter setup (b /= 0)', ok, real(enc%lmax, dp))
        allocate (a(enc%mesh%n, 0:enc%mmax, 0:enc%lmax))
        call encounter_state(enc, 90.0_dp, a)
        ns = encounter_norm(enc, a, SECTOR_SINGLET); nt = encounter_norm(enc, a, SECTOR_TRIPLET)
        n = encounter_norm(enc, a, SECTOR_PRODUCT)
        call check('sector norms add', abs(ns + nt - n) < 1.0e-14_dp, abs(ns + nt - n))
        call check('even/odd l weights ~ 1/2 (separated mirror)', abs(ns - 0.5_dp) < 1.0e-3_dp, ns)
        do ib = 0, nb
            edges(ib) = pi*ib/nb
        end do
        call encounter_angular_bins(enc, a, SECTOR_PRODUCT, 0.0_dp, edges, pp)
        call encounter_angular_bins(enc, a, SECTOR_SINGLET, 0.0_dp, edges, ps)
        call encounter_angular_bins(enc, a, SECTOR_TRIPLET, 0.0_dp, edges, pt)
        err = 0
        do ib = 1, nb
            err = max(err, abs((pp(ib) + pp(nb + 1 - ib))/2*n - (ns*ps(ib) + nt*pt(ib))))
        end do
        call check('spin-unresolved product = weighted singlet+triplet', err < 1.0e-13_dp, err)
        call check('angular bins sum to 1', abs(sum(pp) - 1) < 1.0e-12_dp, abs(sum(pp) - 1))
        do ib = 0, nb
            redges(ib) = enc%mesh%rmax*ib/nb
        end do
        call encounter_separation_density(enc, a, SECTOR_TRIPLET, redges, pr)
        call check('separation density (interpolated) vs DVR norm', abs(sum(pr) - 1) &
                   < 1.0e-6_dp, abs(sum(pr) - 1))
        ae = a; ao = a
        do l = 0, enc%lmax
            if (mod(l, 2) == 1) ae(:, :, l) = 0
            if (mod(l, 2) == 0) ao(:, :, l) = 0
        end do
        x = [5.0_dp, -7.0_dp, 2.0_dp]; y = [1.0_dp, 3.0_dp, -4.0_dp]; z = [-9.0_dp, 6.0_dp, 11.0_dp]
        call encounter_points(enc, ae, x, y, z, psp); call encounter_points(enc, ae, -x, -y, -z, psm)
        call encounter_points(enc, ao, x, y, z, ptp); call encounter_points(enc, ao, -x, -y, -z, ptm)
        err = max(maxval(abs(psp - psm)), maxval(abs(ptp + ptm)))
        call check('singlet even / triplet odd under r -> -r', err < 1.0e-13_dp, err)
    end subroutine test_exchange_sectors

    ! Late-time outgoing angular bins of a near-monochromatic packet (eta = 1,
    ! sigma = 10 lambdabar, b = 0) against plane-wave Rutherford/Mott. Finite-time
    ! bins are not asymptotic cross sections: comparisons are restricted to angles
    ! with impact parameter b(theta) << sigma, and the Rutherford fluence uses the
    ! transverse free-Gaussian density at the collision time z0/v.
    subroutine test_mott_rutherford()
        type(encounter_config_t) :: cfg
        type(encounter_t) :: enc
        character(:), allocatable :: msg
        complex(dp), allocatable :: a(:, :, :)
        integer, parameter :: nb = 18
        real(dp) :: edges(0:nb), pp(nb), ps(nb), pt(nb), t, xq(16), wq(16), d, ex, cr, ru
        real(dp) :: mott_s, mott_d, ruth, width, nperp, err_m, err_r, eb, rad, wall, psym
        integer :: ib, q
        logical :: ok
        cfg%g = 1; cfg%k = 0.5_dp; cfg%sigma = 20; cfg%z0 = 120; cfg%b = 0
        t = 2*cfg%z0/(cfg%hbar*cfg%k/0.5_dp); cfg%t_final = t
        call encounter_init(enc, cfg, ok, msg)
        call check('Mott case setup (eta = 1)', ok .and. abs(enc%eta - 1) < 1.0e-14_dp, enc%eta)
        call check('head-on turning radius r0 = 2 eta lambdabar', &
                   abs(enc%r0 - 2*enc%eta*enc%lambdabar) < 1.0e-12_dp, enc%r0)
        allocate (a(enc%mesh%n, 0:enc%mmax, 0:enc%lmax))
        call encounter_state(enc, t, a)
        wall = encounter_wall_probability(enc, a, SECTOR_PRODUCT, 2*cfg%sigma)
        call check('wall-layer probability at final time [indicator]', wall < 1.0e-4_dp, wall)
        do ib = 0, nb
            edges(ib) = pi*ib/nb
        end do
        call encounter_angular_bins(enc, a, SECTOR_PRODUCT, cfg%z0/4, edges, pp)
        call encounter_angular_bins(enc, a, SECTOR_SINGLET, cfg%z0/4, edges, ps)
        call encounter_angular_bins(enc, a, SECTOR_TRIPLET, cfg%z0/4, edges, pt)
        call ce_free_width(cfg%sigma, enc%mu, cfg%hbar, cfg%z0/enc%speed, width)
        nperp = 1/(2*pi*width**2)
        err_m = 0; err_r = 0
        do ib = 1, nb
            call gauss_legendre_ab(16, cos(edges(ib)), cos(edges(ib - 1)), xq, wq)
            mott_s = 0; mott_d = 0; ruth = 0
            do q = 1, 16
                call ce_mott(acos(xq(q)), enc%eta, cfg%k, d, ex, cr)
                call ce_rutherford(acos(xq(q)), enc%eta, cfg%k, ru)
                mott_s = mott_s + 2*pi*wq(q)*(d + ex + cr)
                mott_d = mott_d + 2*pi*wq(q)*(d + ex)
                ruth = ruth + 2*pi*wq(q)*ru
            end do
            if (edges(ib - 1) >= pi/3 - 1.0e-12_dp .and. edges(ib) <= 2*pi/3 + 1.0e-12_dp) &
                err_m = max(err_m, abs(ps(ib)/((pp(ib) + pp(nb + 1 - ib))/2) - mott_s/mott_d))
            if (edges(ib - 1) >= 2*pi/3 - 1.0e-12_dp) &
                err_r = max(err_r, abs(pp(ib)/(nperp*ruth) - 1))
        end do
        call check('singlet Mott ratio 60-120 deg (abs)', err_m < 0.03_dp, err_m)
        call check('product vs Rutherford fluence 120-180 deg (rel)', err_r < 0.02_dp, err_r)
        psym = (pp(nb/2) + pp(nb/2 + 1))/2
        call check('triplet suppressed at 80-90 deg', pt(nb/2)/psym < 0.05_dp, pt(nb/2)/psym)
        call check('singlet enhanced at 80-90 deg', ps(nb/2)/psym > 1.9_dp, ps(nb/2)/psym)
        call check('spectral window discarded weight [bound]', &
                   enc%ledger%spectral_discarded2 <= cfg%tolerance, enc%ledger%spectral_discarded2)
        call encounter_state_error(enc, t, 1.0_dp, pp(nb/2), eb, rad)
        call check('finite-space state bound e(t) [bound]', eb < 1.0e-4_dp, eb)
        call check('bin radius from state bound [bound]', rad < 1.0e-3_dp, rad)
        call check('initial-domain tail [bound]', enc%ledger%radial_tail2 < 1.0e-12_dp, &
                   enc%ledger%radial_tail2)
    end subroutine test_mott_rutherford

    ! eta = 10, sigma = r0/2: the quantum distinguishable-particle bins approach
    ! the classical Wigner-sampled point-pair ensemble (same finite-time bins).
    subroutine test_classical_ensemble()
        type(encounter_config_t) :: cfg
        type(encounter_t) :: enc
        character(:), allocatable :: msg
        complex(dp), allocatable :: a(:, :, :)
        integer, parameter :: nb = 12, ns = 4000
        real(dp) :: edges(0:nb), pq(nb), pc(nb), t, th, drift, worst, sd, rel(3, ns, 1)
        real(dp) :: cm(3, ns, 1)
        integer :: ib, i
        logical :: ok
        cfg%g = 1; cfg%k = 0.05_dp; cfg%sigma = 10/cfg%k; cfg%z0 = 60/cfg%k
        t = 2*cfg%z0/(cfg%hbar*cfg%k/0.5_dp); cfg%t_final = t
        call encounter_init(enc, cfg, ok, msg)
        allocate (a(enc%mesh%n, 0:enc%mmax, 0:enc%lmax))
        call encounter_state(enc, t, a)
        do ib = 0, nb
            edges(ib) = pi*ib/nb
        end do
        call encounter_angular_bins(enc, a, SECTOR_PRODUCT, cfg%z0/4, edges, pq)
        call encounter_classical_ensemble(enc, ns, 2026_i64, [t], rel, cm, drift, ok)
        call check('classical ensemble integration', ok, drift)
        call check('classical relative energy drift', drift < 1.0e-4_dp, drift)
        pc = 0
        do i = 1, ns
            if (norm2(rel(:, i, 1)) < cfg%z0/4) cycle
            th = acos(rel(3, i, 1)/norm2(rel(:, i, 1)))
            ib = min(nb, int(th/(pi/nb)) + 1)
            pc(ib) = pc(ib) + 1.0_dp/ns
        end do
        worst = 0
        do ib = 1, nb
            sd = sqrt(max(pc(ib)*(1 - pc(ib)), 1.0_dp/ns)/ns)
            worst = max(worst, abs(pq(ib) - pc(ib))/(4*sd + 0.05_dp*pq(ib)))
        end do
        call check('eta=10 quantum vs classical bins (4 sd + 5%)', worst < 1, worst)
    end subroutine test_classical_ensemble

    ! JVP/VJP of an outgoing angular-bin probability w.r.t. (b, sigma, g):
    ! central finite differences on a frozen discretization and the dot-product
    ! identity <grad, v> = JVP(v).
    subroutine test_derivatives()
        type(encounter_config_t) :: cfg, cp, cm
        type(encounter_t) :: enc, ep, em
        character(:), allocatable :: msg
        real(dp) :: t, th1, th2, rcut, p, grad(3), dj, pp, pm, h(3), dir(3), fd, err
        integer :: d
        logical :: ok
        cfg%g = 1; cfg%k = 0.5_dp; cfg%sigma = 8; cfg%z0 = 40; cfg%b = 4
        t = 2*cfg%z0/(cfg%hbar*cfg%k/0.5_dp); cfg%t_final = t
        cfg%rmax = 130
        call encounter_init(enc, cfg, ok, msg)
        th1 = 70*pi/180; th2 = 110*pi/180; rcut = 10
        call encounter_bin_vjp(enc, t, SECTOR_PRODUCT, rcut, th1, th2, p, grad)
        call check('derivative case bin probability', p > 1.0e-3_dp, p)
        h = [1.0e-3_dp, 1.0e-3_dp, 1.0e-4_dp]
        err = 0
        do d = 1, 3
            dir = 0; dir(d) = 1
            call encounter_bin_jvp(enc, t, SECTOR_PRODUCT, rcut, th1, th2, dir, p, dj)
            cp = cfg; cm = cfg
            cp%hmax = enc%hmax; cm%hmax = enc%hmax
            cp%lmax = enc%lmax; cp%mmax = enc%mmax; cm%lmax = enc%lmax; cm%mmax = enc%mmax
            select case (d)
            case (1)
                cp%b = cfg%b + h(1); cm%b = cfg%b - h(1)
            case (2)
                cp%sigma = cfg%sigma + h(2); cm%sigma = cfg%sigma - h(2)
            case default
                cp%g = cfg%g + h(3); cm%g = cfg%g - h(3)
            end select
            call encounter_init(ep, cp, ok, msg)
            call encounter_init(em, cm, ok, msg)
            call encounter_bin_full(ep, t, SECTOR_PRODUCT, rcut, th1, th2, pp)
            call encounter_bin_full(em, t, SECTOR_PRODUCT, rcut, th1, th2, pm)
            fd = (pp - pm)/(2*h(d))
            err = max(err, abs(dj - fd)/abs(fd), abs(grad(d) - fd)/abs(fd))
        end do
        call check('JVP/VJP vs central differences (b, sigma, g)', err < 1.0e-6_dp, err)
        dir = [0.3_dp, -0.7_dp, 0.2_dp]
        call encounter_bin_jvp(enc, t, SECTOR_SINGLET, rcut, th1, th2, dir, p, dj)
        call encounter_bin_vjp(enc, t, SECTOR_SINGLET, rcut, th1, th2, p, grad)
        call check('dot-product identity (singlet sector)', &
                   abs(dj - sum(grad*dir)) < 1.0e-12_dp*max(1.0_dp, abs(dj)), abs(dj - sum(grad*dir)))
    end subroutine test_derivatives

    ! Q1b: hydrogen-like radial Kepler packet, nbar = 40, l = 1, rms Delta n ~ 1.5.
    ! Breathing amplitude per Kepler period: collapse, half revival at T_rev/2,
    ! (cubic-dispersion degraded) full revival near T_rev = (2 nbar/3) T_K.
    subroutine test_kepler()
        type(kepler_t) :: kep
        logical :: ok
        integer, parameter :: nt = 3000
        real(dp) :: t, dmax, tret, best, a2, amp0, ampc, amph, ampf, amp, eb, rm
        complex(dp) :: an, ao
        integer :: i, k, nper
        call kepler_init(kep, -1.0_dp, 1.0_dp, 1.0_dp, 1, 40.0_dp, 1.5_dp, 12, ok)
        call check('Kepler packet setup', ok, real(kep%nbound, dp))
        call check('identified levels vs exact hydrogen (weighted states)', &
                   kep%max_level_error < 1.0e-10_dp, kep%max_level_error)
        call check('continuum + window weight', kep%continuum_weight + kep%window_discarded2 &
                   < 1.0e-10_dp, kep%continuum_weight + kep%window_discarded2)
        call check('mean level near nbar', abs(kep%mean_n - 40) < 0.3_dp, kep%mean_n)
        dmax = 0; best = 0; tret = 0
        do i = 0, nt
            t = 1.1_dp*kep%revival*i/nt
            call kepler_autocorrelation(kep, t, an, ao)
            a2 = abs(an)**2
            dmax = max(dmax, abs(an - ao))
            if (t > 0.5_dp*kep%period .and. t < 1.5_dp*kep%period .and. a2 > best) then
                best = a2; tret = t
            end if
        end do
        call check('autocorrelation vs exact-level oracle (t <= 1.1 T_rev)', dmax < 1.0e-8_dp, dmax)
        call check('first return near Kepler period T_K', abs(tret/kep%period - 1) < 0.05_dp, &
                   tret/kep%period)
        nper = ceiling(1.1_dp*kep%revival/kep%period)
        amp0 = kepler_breathing(kep, 0.0_dp, kep%period, 80)
        ampc = huge(1.0_dp); amph = 0; ampf = 0
        do k = 1, nper - 1
            t = (k + 0.5_dp)*kep%period/kep%revival
            amp = kepler_breathing(kep, k*kep%period, (k + 1)*kep%period, 80)
            if (t > 0.1_dp .and. t < 0.35_dp) ampc = min(ampc, amp)
            if (t > 0.4_dp .and. t < 0.6_dp) amph = max(amph, amp)
            if (t > 0.85_dp .and. t < 1.1_dp) ampf = max(ampf, amp)
        end do
        call check('collapse: breathing amplitude ratio', ampc < 0.3_dp*amp0, ampc/amp0)
        call check('half revival near T_rev/2', amph > 0.6_dp*amp0, amph/amp0)
        call check('full revival near T_rev (above collapse)', ampf > 3*ampc, ampf/amp0)
        call kepler_mean_radius(kep, 0.5_dp*kep%period, rm, eb)
        call check('Kepler finite-space state bound at T_K/2 [bound]', eb < 1.0e-5_dp, eb)
    end subroutine test_kepler
end program test_coulomb_encounter
