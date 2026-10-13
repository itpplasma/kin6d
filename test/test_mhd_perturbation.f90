program test_mhd_perturbation
    !! L7 oracles on the analytic circular equilibrium psi = c ((R-R0)^2 + Z^2),
    !! F = const (kinematic: no force balance is needed):
    !!  1. a known complex displacement xi e^{i n phi} -> dB = curl(xi x B0) by
    !!     central differences of the cross product (independent of the
    !!     reconstruction formulas) -> recovered normal and binormal parts;
    !!  2. dB_L = dB_E + (xi.grad) B0, dBmag_L = b.dB_E + xi.grad|B| and the
    !!     curvature against finite displacements and field-line steps of the
    !!     Cartesian fields (rotation terms checked independently);
    !!  3. SFL chart maps (PEST, Hamada) and the refusal at a rational surface;
    !!  4. the export file hash (SHA-256 vectors, sidecar, tamper detection).
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_vacuum_surfaces
    use kin6d_mhd_ideal_n, only: ideal_options_t
    use kin6d_mhd_perturbation_n
    implicit none
    real(dp), parameter :: pi = 3.14159265358979323846_dp
    real(dp), parameter :: r0 = 6.5_dp, c = 0.8_dp, amax = 1.8_dp, fconst = 33.0_dp
    real(dp), parameter :: ip = 1.5e7_dp
    integer, parameter :: n = 3
    integer :: failures = 0
    type(psi_grid_t) :: g
    type(ideal_options_t) :: opt

    call make_grid(g)
    call test_gradient()
    call test_round_trip_and_identities()
    call test_charts_and_resonance()
    call test_hash()
    if (failures > 0) then
        print '(a,i0,a)', 'test_mhd_perturbation: ', failures, ' failure(s)'
        error stop 1
    end if
    print '(a)', 'test_mhd_perturbation: all passed'

contains

    subroutine check(name, value, limit)
        character(len=*), intent(in) :: name
        real(dp), intent(in) :: value, limit
        logical :: ok
        ok = value <= limit
        print '(a,a,es12.4,a,es10.2,a)', merge('ok   ', 'FAIL ', ok), name, value, &
            ' (limit ', limit, ')'
        if (.not. ok) failures = failures + 1
    end subroutine

    subroutine make_grid(gg)
        type(psi_grid_t), intent(out) :: gg
        real(dp), allocatable :: psirz(:, :), f(:)
        integer, parameter :: nr = 201, nz = 201
        integer :: i, j, st
        allocate (psirz(nr, nz), f(nr))
        do j = 1, nz
            do i = 1, nr
                psirz(i, j) = c*((4.5_dp + 4.0_dp*(i - 1)/(nr - 1) - r0)**2 &
                                 + (-2.0_dp + 4.0_dp*(j - 1)/(nz - 1))**2)
            end do
        end do
        f = fconst
        call psi_grid_init(gg, psirz, 4.5_dp, 4.0_dp, 0.0_dp, 4.0_dp, 6.45_dp, &
            0.05_dp, c*amax**2, ip, f, st)
        if (st /= 0) error stop 'grid'
    end subroutine

    function b0_exact(r, z) result(b)
        real(dp), intent(in) :: r, z
        real(dp) :: b(3)
        b = [g%s_pol*2*c*z/r, fconst/r, -g%s_pol*2*c*(r - r0)/r]
    end function

    subroutine gradb_exact(r, z, gb)
        real(dp), intent(in) :: r, z
        real(dp), intent(out) :: gb(3, 2)
        gb(1, 1) = -g%s_pol*2*c*z/r**2
        gb(1, 2) = g%s_pol*2*c/r
        gb(2, 1) = -fconst/r**2
        gb(2, 2) = 0.0_dp
        gb(3, 1) = -g%s_pol*2*c*(1.0_dp/r - (r - r0)/r**2)
        gb(3, 2) = 0.0_dp
    end subroutine

    function xi_fun(r, z) result(v)
        !! known displacement amplitude (R, phi, Z), all three components
        real(dp), intent(in) :: r, z
        complex(dp) :: v(3)
        v(1) = cmplx(0.010_dp, 0.004_dp, dp)*cos(0.9_dp*z)*(1 + 0.08_dp*(r - r0))
        v(2) = cmplx(-0.006_dp, 0.003_dp, dp)*sin(0.7_dp*(r - r0))*cos(0.5_dp*z)
        v(3) = cmplx(0.008_dp, -0.002_dp, dp)*(1 + 0.1_dp*z)*sin(0.6_dp*(r - r0) + 0.3_dp)
    end function

    function e_field(r, z) result(e)
        !! xi x B0 (R, phi, Z)
        real(dp), intent(in) :: r, z
        complex(dp) :: e(3), v(3)
        real(dp) :: b(3)
        v = xi_fun(r, z)
        b = b0_exact(r, z)
        e(1) = v(2)*b(3) - v(3)*b(2)
        e(2) = v(3)*b(1) - v(1)*b(3)
        e(3) = v(1)*b(2) - v(2)*b(1)
    end function

    function curl_xb(r, z) result(db)
        !! curl(xi x B0) e^{i n phi} by central differences in (R, Z)
        real(dp), intent(in) :: r, z
        complex(dp) :: db(3), ep(3), em(3)
        real(dp), parameter :: h = 1.0e-4_dp
        complex(dp), parameter :: ii = (0.0_dp, 1.0_dp)
        complex(dp) :: e0(3), dezdr, dephidz, dercdz, dephidr
        e0 = e_field(r, z)
        ep = e_field(r, z + h); em = e_field(r, z - h)
        dephidz = (ep(2) - em(2))/(2*h)
        dercdz = (ep(1) - em(1))/(2*h)
        db(1) = ii*n*e0(3)/r - dephidz
        ep = e_field(r + h, z); em = e_field(r - h, z)
        dezdr = (ep(3) - em(3))/(2*h)
        dephidr = ((r + h)*ep(2) - (r - h)*em(2))/(2*h)
        db(2) = dercdz - dezdr
        db(3) = (dephidr - ii*n*e0(1))/r
    end function

    subroutine test_gradient()
        real(dp) :: b(3), gp(2), gb(3, 2), gx(3, 2), e
        integer :: i
        e = 0.0_dp
        do i = 1, 5
            call equilibrium_gradient(g, opt, r0 + 0.3_dp*i, 0.2_dp*i - 0.5_dp, b, gp, gb)
            call gradb_exact(r0 + 0.3_dp*i, 0.2_dp*i - 0.5_dp, gx)
            e = max(e, maxval(abs(gb - gx))/maxval(abs(gx)))
        end do
        call check('grad B0 (spline, central difference) vs exact', e, 1.0e-7_dp)
    end subroutine

    subroutine test_round_trip_and_identities()
        type(surface_t) :: s
        type(perturbation_t) :: res
        complex(dp), allocatable :: dbe(:, :)
        real(dp) :: e_n, e_b, e_l, e_m, e_k, xmax, lmax, kmax, nh(3), eb(3)
        real(dp) :: sc
        complex(dp) :: v(3), xt
        integer :: k, st, np
        np = 128
        call surface_build(g, 0.5_dp, np, s, st)
        allocate (dbe(3, np))
        do k = 1, np
            dbe(:, k) = curl_xb(s%r(k), s%z(k))
        end do
        call perturbation_on_surface(g, opt, s, n, dbe, res, st)
        if (st /= 0) error stop 'perturbation_on_surface'
        call check('q_eff vs closed form', abs(abs(res%q_eff) - fconst/ &
            (2*c*sqrt(r0**2 - (amax**2*0.5_dp)))), 1.0e-9_dp)
        e_n = 0; e_b = 0; xmax = 0
        do k = 1, np
            v = xi_fun(s%r(k), s%z(k))
            nh = res%nhat(:, k); eb = res%ebin(:, k)
            xt = sum(nh*v)
            e_n = max(e_n, abs(res%xi_n(k) - xt))
            e_b = max(e_b, abs(res%xi_b(k) - sum(eb*v)))
            xmax = max(xmax, abs(xt), abs(sum(eb*v)))
        end do
        call check('xi_n round trip (rel. to |xi|)', e_n/xmax, 2.0e-6_dp)
        call check('xi_b round trip (rel. to |xi|)', e_b/xmax, 2.0e-6_dp)

        ! Lagrangian identities against finite displacements (Cartesian, phi=0)
        sc = 1.0e-3_dp
        e_l = 0; e_m = 0; e_k = 0; lmax = 0; kmax = 0
        do k = 1, np
            block
                real(dp) :: xi_re(3), bplus(3), bminus(3), lag(3), dl(3)
                real(dp) :: dmag, mp, mm, kc(3), bh(3)
                xi_re = real(res%xi(:, k))
                bplus = cart_total(s%r(k), s%z(k), xi_re, sc)
                bminus = cart_total(s%r(k), s%z(k), xi_re, -sc)
                lag = (bplus - bminus)/(2*sc)
                dl = real(res%dbl(:, k))
                e_l = max(e_l, maxval(abs(lag - dl)))
                lmax = max(lmax, maxval(abs(dl)))
                mp = norm2(bplus)
                mm = norm2(bminus)
                dmag = (mp - mm)/(2*sc)
                e_m = max(e_m, abs(dmag - real(res%dbmag_l(k))))
                ! curvature = d b/dl along the field line (Cartesian)
                bh = b0_exact(s%r(k), s%z(k))
                bh = bh/norm2(bh)
                kc = curvature_fd(s%r(k), s%z(k))
                e_k = max(e_k, maxval(abs(kc - res%kappa(:, k))))
                kmax = max(kmax, maxval(abs(kc)))
            end block
        end do
        call check('dB_L = dB_E + xi.grad B0 (vector, finite displacement)', &
            e_l/lmax, 1.0e-6_dp)
        call check('dBmag_L = b.dB_E + xi.grad|B| (finite displacement)', &
            e_m/lmax, 1.0e-6_dp)
        call check('kappa = (b.grad)b (field-line step)', e_k/kmax, 1.0e-6_dp)
        ! kappa.xi
        e_k = 0
        do k = 1, np
            e_k = max(e_k, abs(res%kappa_xi(k) - sum(res%kappa(:, k)*res%xi(:, k))))
        end do
        call check('kappa.xi consistency', e_k, 1.0e-14_dp)
        ! kappa_n, kappa_b components reproduce kappa.xi_perp
        e_k = 0
        do k = 1, np
            e_k = max(e_k, abs(res%kappa_xi(k) - (sum(res%kappa(:, k)*res%nhat(:, k))* &
                res%xi_n(k) + sum(res%kappa(:, k)*res%ebin(:, k))*res%xi_b(k))))
        end do
        call check('kappa.xi = kappa_n xi_n + kappa_b xi_b', e_k, 1.0e-14_dp)

        ! export and hash
        call perturbation_add_spectra(s, res, -30, 30)
        ! phi_n: one-sided flux-density harmonics = half the peak-amplitude phi_mn
        ! of surface_harmonics (independent implementation), both charts
        block
            complex(dp) :: bmn(-30:30), pmn(-30:30)
            integer :: ic
            do ic = angle_pest, angle_hamada
                call surface_harmonics(s, n, dbe([1, 3], :), ic, -30, 30, bmn, pmn)
                call check('phi_n spectrum = phi_mn/2 of surface_harmonics', &
                    maxval(abs(res%spec(:, nq_spec, ic) - 0.5_dp*pmn))/ &
                    maxval(abs(pmn)), 1.0e-12_dp)
            end do
        end block
        call perturbation_export('l7_test_export.dat', [res], 'test_mhd_perturbation', &
            3, g%s_pol, 1.0_dp, st)
        call check('export written', real(st, dp), 0.0_dp)
        block
            logical :: ok
            call sha256_verify('l7_test_export.dat', ok)
            call check('sidecar hash verifies', merge(0.0_dp, 1.0_dp, ok), 0.0_dp)
        end block
    end subroutine

    function cart_total(r, z, xi, sc) result(bt)
        !! Cartesian B_tot(x + sc xi) with B_tot = B0 + sc Re(dB e^{i n phi}),
        !! x = (R, 0, Z); xi is the real (R, phi, Z) displacement at x.
        real(dp), intent(in) :: r, z, xi(3), sc
        real(dp) :: bt(3), rp, php, zp, b0(3), dbc(3), cyl(3), ca, sa
        complex(dp) :: dbh(3)
        rp = sqrt((r + sc*xi(1))**2 + (sc*xi(2))**2)
        php = atan2(sc*xi(2), r + sc*xi(1))
        zp = z + sc*xi(3)
        b0 = b0_exact(rp, zp)
        dbh = curl_xb(rp, zp)
        dbc = real(dbh*exp(cmplx(0.0_dp, n*php, dp)))
        cyl = b0 + sc*dbc
        ca = cos(php); sa = sin(php)
        bt = [cyl(1)*ca - cyl(2)*sa, cyl(1)*sa + cyl(2)*ca, cyl(3)]
    end function

    function curvature_fd(r, z) result(k)
        !! (b.grad) b at phi = 0 by a Cartesian step along b; (R, phi, Z).
        real(dp), intent(in) :: r, z
        real(dp) :: k(3), b(3), bn(3), hp(3), hm(3), h
        real(dp), parameter :: step = 1.0e-5_dp
        b = b0_exact(r, z)
        bn = b/norm2(b)
        h = step
        hp = bhat_cart(r + h*bn(1), h*bn(2), z + h*bn(3))
        hm = bhat_cart(r - h*bn(1), -h*bn(2), z - h*bn(3))
        k = (hp - hm)/(2*h)
    end function

    function bhat_cart(x, y, z) result(bh)
        real(dp), intent(in) :: x, y, z
        real(dp) :: bh(3), r, ph, b(3)
        r = hypot(x, y)
        ph = atan2(y, x)
        b = b0_exact(r, z)
        b = b/norm2(b)
        bh = [b(1)*cos(ph) - b(2)*sin(ph), b(1)*sin(ph) + b(2)*cos(ph), b(3)]
    end function

    subroutine test_charts_and_resonance()
        type(surface_t) :: s
        type(perturbation_t) :: res
        complex(dp), allocatable :: f(:), f2(:), fm_p(:), fm_h(:), fm_h2(:), fm_p2(:), dbe(:, :)
        complex(dp), allocatable :: x(:)
        integer, parameter :: np = 256, mm = 40
        real(dp) :: e1, e2, e3, e4, th, psir, qex
        integer :: j, m, ic, st, nres, mres
        allocate (f(np), f2(np), fm_p(-mm:mm), fm_h(-mm:mm), fm_h2(-mm:mm), fm_p2(-mm:mm))
        call surface_build(g, 0.5_dp, np, s, st)
        ! single harmonic in each chart
        do ic = 1, 2
            allocate (x(np))
            call perturbation_synthesis(s, ic, n, [(merge((1.0_dp, 0.0_dp), &
                (0.0_dp, 0.0_dp), m == 5), m=-mm, mm)], -mm, mm, x)
            call perturbation_spectrum(s, ic, n, x, -mm, mm, fm_p)
            fm_p(5) = fm_p(5) - 1.0_dp
            e1 = maxval(abs(fm_p))
            call check('single harmonic recovered in chart (1 PEST, 2 Hamada)', e1, &
                1.0e-10_dp)
            deallocate (x)
        end do
        ! smooth field: PEST -> Hamada -> PEST
        do j = 1, np
            th = s%theta_g(j)
            f(j) = cmplx(1 + 0.3_dp*cos(th) + 0.2_dp*sin(2*th), 0.5_dp - 0.1_dp*cos(3*th), dp)
        end do
        call perturbation_spectrum(s, angle_pest, n, f, -mm, mm, fm_p)
        call perturbation_spectrum(s, angle_hamada, n, f, -mm, mm, fm_h)
        call perturbation_chart_map(s, n, angle_pest, angle_hamada, fm_p, -mm, mm, fm_h2)
        call perturbation_chart_map(s, n, angle_hamada, angle_pest, fm_h2, -mm, mm, fm_p2)
        e2 = maxval(abs(fm_h2 - fm_h))/maxval(abs(fm_h))
        e3 = maxval(abs(fm_p2 - fm_p))/maxval(abs(fm_p))
        call check('chart map PEST -> Hamada vs direct Hamada spectrum', e2, 1.0e-8_dp)
        call check('chart map round trip PEST -> Hamada -> PEST', e3, 1.0e-8_dp)
        ! synthesis reproduces the field
        call perturbation_synthesis(s, angle_pest, n, fm_p, -mm, mm, f2)
        call check('synthesis of PEST spectrum reproduces the field', &
            maxval(abs(f2 - f))/maxval(abs(f)), 1.0e-9_dp)
        ! resonance refusal: n = 7 at q = m/7
        nres = 7
        qex = 23.0_dp/7.0_dp
        psir = (r0**2 - (fconst/(2*c*qex))**2)/amax**2
        call surface_build(g, psir, 128, s, st)
        allocate (dbe(3, 128))
        dbe = (1.0e-3_dp, 0.0_dp)
        call perturbation_on_surface(g, opt, s, nres, dbe, res, st)
        mres = nint(nres*abs(res%q_eff))
        e4 = abs(nres*abs(res%q_eff) - mres)
        call check('rational surface located (n q - m)', e4, 1.0e-4_dp)
        call check('resonant surface refused (status 2)', real(abs(st - 2), dp), 0.0_dp)
    end subroutine

    subroutine test_hash()
        character(len=64) :: hex
        integer :: u, st, i
        character(len=1), allocatable :: big(:)
        logical :: ok
        open (newunit=u, file='l7_hash_abc.txt', access='stream', form='unformatted', &
            status='replace')
        write (u) 'abc'
        close (u)
        call sha256_file('l7_hash_abc.txt', hex, st)
        call check('sha256("abc")', merge(0.0_dp, 1.0_dp, hex == &
            'ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad'), 0.0_dp)
        allocate (big(1000000))
        big = 'a'
        open (newunit=u, file='l7_hash_big.txt', access='stream', form='unformatted', &
            status='replace')
        write (u) big
        close (u)
        call sha256_file('l7_hash_big.txt', hex, st)
        call check('sha256(10^6 x "a")', merge(0.0_dp, 1.0_dp, hex == &
            'cdc76e5c9914fb9281a1c7e284d73e67f1809a48a497200e046d39ccc7112cd0'), 0.0_dp)
        ! tamper detection on the export file
        open (newunit=u, file='l7_test_export.dat', access='stream', form='unformatted', &
            status='old', action='readwrite')
        write (u, pos=3) 'X'
        close (u)
        call sha256_verify('l7_test_export.dat', ok)
        i = merge(0, 1, .not. ok)
        call check('tampered export fails the sidecar hash', real(i, dp), 0.0_dp)
    end subroutine

end program test_mhd_perturbation
