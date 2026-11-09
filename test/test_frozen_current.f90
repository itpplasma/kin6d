program test_frozen_current
    ! Oracles of the frozen-current field (milca_frozen_current).
    !   ring k nr      smooth axisymmetric toroidal + poloidal current in a
    !                  disk; B_pol against the exact ring field integrated
    !                  over the cross-section (AGM elliptic integrals, polar
    !                  Gauss quadrature, all in this file), B_phi = mu0 F/R
    !                  and the toroidal flux of the poloidal current, which
    !                  fixes the mode-0 circulation.
    !   helical k      manufactured A = (0, w, U(phi)) with modes 0..2 and
    !                  exact boundary data; B errors at nr = 2, 4, 8 and the
    !                  observed orders.
    !   exterior k nr nmax tol
    !                  localized current with all toroidal harmonics; the
    !                  exterior field from the FE boundary traces (Kirchhoff
    !                  surface integral) against direct Biot-Savart of J at
    !                  points outside the domain, for N = 2, 4, .., nmax.
    !   source k nr tol
    !                  two-periodic current J = curl F/mu0: exterior data
    !                  from samples of F (fc_casing) against direct
    !                  Biot-Savart of J, the fields solved with either, and
    !                  the period-2 solve against the period-1 solve with
    !                  all modes up to 4.
    ! The removed non-solenoidal fraction is printed at every resolution.
    ! ncore = 16 keeps the polar core cells close to isotropic.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_mhd_mesh, only: flux_mesh_t, flux_mesh_build
    use milca_equilibrium_nonnested, only: free_b_field_t, nn_build, nn_point
    use milca_equilibrium_state, only: EQ_VALUES, EQ_BR, EQ_BZ, EQ_BPHI
    use milca_frozen_current, only: fc_exterior_t, fc_report_t, fc_exterior, &
        fc_solve, fc_source_t, fc_casing
    use fortfem_meridian_feec_n, only: meridian_boundary_trace_points
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp), mu0 = 4.0e-7_dp*pi
    real(dp), parameter :: r0 = 3, a = 0.4_dp, bdom = 0.6_dp
    real(dp), parameter :: j0 = 1.0e6_dp, f0 = 1.0e5_dp, h0 = 2.0e5_dp
    character(32) :: arg
    integer :: k, nr, ncore = 16, source_period = 1
    real(dp) :: mu = mu0

    call get_command_argument(1, arg)
    k = int_arg(2, 3); nr = int_arg(3, 8)
    select case (trim(arg))
    case ('ring')
        call ring_case(k, nr)
    case ('helical')
        call helical_case(k)
    case ('exterior')
        call exterior_case(k, nr, int_arg(4, 4), real_arg(5, 1.5e-3_dp))
    case ('source')
        call source_case(k, nr, real_arg(4, 1.0e-3_dp))
    case default
        error stop 'usage: test_frozen_current ring|helical|exterior|source k [nr]'
    end select
contains
    integer function int_arg(i, default) result(v)
        integer, intent(in) :: i, default
        character(16) :: s
        v = default
        call get_command_argument(i, s)
        if (len_trim(s) > 0) read (s, *) v
    end function

    real(dp) function real_arg(i, default) result(v)
        integer, intent(in) :: i
        real(dp), intent(in) :: default
        character(32) :: s
        v = default
        call get_command_argument(i, s)
        if (len_trim(s) > 0) read (s, *) v
    end function

    subroutine disk_field(k, nr, modes, field)
        integer, intent(in) :: k, nr, modes
        type(free_b_field_t), intent(out) :: field
        type(flux_mesh_t) :: fm
        real(dp) :: rings(nr)
        integer :: i, st
        do i = 1, nr
            rings(i) = bdom*i/nr
        end do
        call flux_mesh_build(rings, ncore, 4096, min(k, 3), disk_map, fm, st)
        if (st /= 0) error stop 'mesh'
        call nn_build(fm, k, modes, [r0, 0.0_dp], 1.0_dp, 0.0_dp, mu, field, st)
        if (st /= 0) error stop 'build'
    end subroutine

    subroutine disk_map(u, theta, x)
        real(dp), intent(in) :: u, theta
        real(dp), intent(out) :: x(2)
        x = [r0 + u*cos(theta), u*sin(theta)]
    end subroutine

    function fe_b(field, x, phi) result(b)
        type(free_b_field_t), intent(in) :: field
        real(dp), intent(in) :: x(2), phi
        real(dp) :: b(3), v(EQ_VALUES)
        integer :: st
        call nn_point(field, [x, phi], v, st)
        if (st /= 0) error stop 'nn_point'
        b = [v(EQ_BR), v(EQ_BPHI), v(EQ_BZ)]
    end function

    ! ---------------------------------------------------------------- ring
    subroutine ring_current(x, phi, j)
        real(dp), intent(in) :: x(2), phi
        real(dp), intent(out) :: j(3)
        real(dp) :: q, fr, fz
        q = 1 - ((x(1) - r0)**2 + x(2)**2)/a**2
        j = 0
        if (q <= 0) return
        ! J_pol = grad F x grad phi with F = f0 q^4.
        fr = f0*4*q**3*(-2*(x(1) - r0)/a**2)
        fz = f0*4*q**3*(-2*x(2)/a**2)
        j = [-fz/x(1), j0*q**3, fr/x(1)]
        if (.false.) j(1) = phi
    end subroutine

    subroutine ring_case(k, nr)
        integer, intent(in) :: k, nr
        type(free_b_field_t) :: field
        type(fc_exterior_t) :: ext
        type(fc_report_t) :: report
        real(dp), allocatable :: x(:)
        real(dp), parameter :: outer(3) = [0.45_dp, 0.5_dp, 0.57_dp]
        real(dp), parameter :: inner(4) = [0.05_dp, 0.15_dp, 0.25_dp, 0.35_dp]
        real(dp) :: p(2), b(3), exact(2), errp, errt, scale, flux, q, rho, th
        real(dp) :: xg(48), wg(48)
        integer :: st, i, l
        call disk_field(k, nr, 0, field)
        call fc_exterior(field, ring_current, 16, ext, st)
        if (st /= 0) error stop 'exterior data'
        call fc_solve(field, ring_current, ext, x, report, st)
        if (st /= 0) error stop 'solve'
        errp = 0; errt = 0; scale = 0
        do i = 1, 3
            rho = outer(i)
            do l = 0, 7
                th = 2*pi*l/8 + 0.1_dp
                p = [r0 + rho*cos(th), rho*sin(th)]
                b = fe_b(field, p, 0.0_dp)
                exact = ring_integral(p)
                errp = max(errp, norm2(b([1, 3]) - exact))
                scale = max(scale, norm2(exact))
                errt = max(errt, abs(b(2)))
            end do
        end do
        do i = 1, 4
            rho = inner(i)
            do l = 0, 5
                th = 2*pi*l/6 + 0.3_dp
                p = [r0 + rho*cos(th), rho*sin(th)]
                b = fe_b(field, p, 0.0_dp)
                q = 1 - rho**2/a**2
                errt = max(errt, abs(b(2) - mu*f0*q**4/p(1)))
            end do
        end do
        ! Exact toroidal flux mu0 int F/R dR dZ of the poloidal current.
        flux = 0
        call gauss(48, xg, wg)
        do i = 1, 48
            rho = a*(xg(i) + 1)/2
            do l = 1, 128
                th = 2*pi*(l - 1)/128
                q = 1 - rho**2/a**2
                flux = flux + mu*f0*q**4/(r0 + rho*cos(th))*rho*(a/2)*wg(i) &
                    *(2*pi/128)
            end do
        end do
        print '(a,i0,a,i0)', 'ring P', k, ' nr ', nr
        print '(a,es12.4)', '  B_pol relative error outside support ', errp/scale
        print '(a,es12.4)', '  B_phi error / max B_phi             ', &
            errt/(mu*f0/(r0 - a))
        print '(a,2es16.8)', '  toroidal flux FE, exact [Wb]         ', &
            -report%circulation, flux
        print '(a,es12.4)', '  removed non-solenoidal fraction     ', report%projection
        print '(a,es12.4)', '  linear residual                     ', report%residual
        if (errp/scale > 2.0e-3_dp) error stop 'B_pol against the ring integral'
        if (errt/(mu*f0/(r0 - a)) > 2.0e-2_dp) error stop 'B_phi against mu0 F/R'
        if (abs(-report%circulation/flux - 1) > 1.0e-5_dp) error stop 'toroidal flux'
        if (report%projection > 1.0e-6_dp) error stop 'projection of a solenoidal J'
    end subroutine

    function ring_integral(p) result(b)
        ! B_pol of J_phi integrated over the disk rho' < a: Gauss-Legendre in
        ! rho' (the integrand is polynomial in rho'^2 times a smooth kernel)
        ! and the trapezoid rule in theta'.
        real(dp), intent(in) :: p(2)
        real(dp) :: b(2), xg(48), wg(48), rho, th, q, src(2), br, bz
        integer :: i, l
        call gauss(48, xg, wg)
        b = 0
        do i = 1, 48
            rho = a*(xg(i) + 1)/2
            q = 1 - rho**2/a**2
            do l = 1, 256
                th = 2*pi*(l - 1)/256
                src = [r0 + rho*cos(th), rho*sin(th)]
                call loop_field(src(1), p(1), p(2) - src(2), br, bz)
                b = b + [br, bz]*j0*q**3*rho*(a/2)*wg(i)*(2*pi/256)
            end do
        end do
    end function

    subroutine loop_field(rl, r, z, br, bz)
        ! Field of a unit-current circular loop of radius rl at height 0.
        real(dp), intent(in) :: rl, r, z
        real(dp), intent(out) :: br, bz
        real(dp) :: alpha2, beta, m, kk, ee, c
        alpha2 = rl**2 + r**2 + z**2 - 2*rl*r
        beta = sqrt(rl**2 + r**2 + z**2 + 2*rl*r)
        m = 1 - alpha2/beta**2
        call elliptic(m, kk, ee)
        c = mu/pi
        br = c*z/(2*alpha2*beta*r)*((rl**2 + r**2 + z**2)*ee - alpha2*kk)
        bz = c/(2*alpha2*beta)*((rl**2 - r**2 - z**2)*ee + alpha2*kk)
    end subroutine

    subroutine elliptic(m, kk, ee)
        ! Complete elliptic integrals K(m), E(m) by the AGM.
        real(dp), intent(in) :: m
        real(dp), intent(out) :: kk, ee
        real(dp) :: aa, bb, cc, an, sum, pw
        aa = 1; bb = sqrt(1 - m); sum = m/2; pw = 0.5_dp
        do while (abs(aa - bb) > 1.0e-15_dp*aa)
            an = (aa + bb)/2
            cc = (aa - bb)/2
            bb = sqrt(aa*bb); aa = an
            pw = 2*pw
            sum = sum + pw*cc**2
        end do
        kk = pi/(2*aa)
        ee = kk*(1 - sum)
    end subroutine

    subroutine gauss(n, x, w)
        ! Gauss-Legendre nodes and weights on [-1, 1] by Newton.
        integer, intent(in) :: n
        real(dp), intent(out) :: x(n), w(n)
        real(dp) :: z, p0, p1, p2, dp1
        integer :: i, j, it
        do i = 1, n
            z = cos(pi*(i - 0.25_dp)/(n + 0.5_dp))
            do it = 1, 100
                p0 = 1; p1 = z
                do j = 2, n
                    p2 = ((2*j - 1)*z*p1 - (j - 1)*p0)/j
                    p0 = p1; p1 = p2
                end do
                dp1 = n*(z*p1 - p0)/(z*z - 1)
                if (abs(p1/dp1) < 1.0e-16_dp) exit
                z = z - p1/dp1
            end do
            x(i) = z; w(i) = 2/((1 - z*z)*dp1**2)
        end do
    end subroutine

    ! -------------------------------------------------------------- source
    ! F = (-psi_Z/R, H/R, psi_R/R) + grad chi on the torus rho < a, with
    ! psi = c1 q (1 + e sin P phi), H = c2 q and the harmonic
    ! chi = c3 Z R^P cos P phi: div F = 0, J = curl F/mu0 (fourth-order
    ! differences of F) jumps at rho = a, and modes 0 and P have a field
    ! outside.
    subroutine source_field(x, phi, f)
        real(dp), intent(in) :: x(2), phi
        real(dp), intent(out) :: f(3)
        real(dp), parameter :: c1 = 2.0e-3_dp, c2 = 4.0e-2_dp, c3 = 5.0e-3_dp
        real(dp), parameter :: e = 0.4_dp
        real(dp) :: q, qr, qz, r, z, s
        integer :: p
        p = source_period
        r = x(1); z = x(2)
        q = 1 - ((r - r0)**2 + z**2)/a**2
        f = c3*[p*z*r**(p - 1)*cos(p*phi), -p*z*r**(p - 1)*sin(p*phi), &
            r**p*cos(p*phi)]
        if (q <= 0) return
        qr = -2*(r - r0)/a**2; qz = -2*z/a**2
        s = c1*(1 + e*sin(p*phi))
        f = f + [-s*qz/r, c2*q/r, s*qr/r]
    end subroutine

    subroutine source_current(x, phi, j)
        real(dp), intent(in) :: x(2), phi
        real(dp), intent(out) :: j(3)
        real(dp) :: dr(3), dz(3), dq(3), f(3)
        j = 0
        if ((x(1) - r0)**2 + x(2)**2 >= a**2) return
        dr = source_diff(x, phi, [1.0_dp, 0.0_dp], 0.0_dp)
        dz = source_diff(x, phi, [0.0_dp, 1.0_dp], 0.0_dp)
        dq = source_diff(x, phi, [0.0_dp, 0.0_dp], 1.0_dp)
        call source_field(x, phi, f)
        j = [dq(3)/x(1) - dz(2), dz(1) - dr(3), &
            f(2)/x(1) + dr(2) - dq(1)/x(1)]/mu
    end subroutine

    function source_diff(x, phi, dx, dphi) result(d)
        ! Fourth-order central difference of F along (dx, dphi).
        real(dp), intent(in) :: x(2), phi, dx(2), dphi
        real(dp), parameter :: h = 1.0e-4_dp
        real(dp) :: d(3), f1(3), f2(3), f3(3), f4(3)
        call source_field(x + 2*h*dx, phi + 2*h*dphi, f1)
        call source_field(x + h*dx, phi + h*dphi, f2)
        call source_field(x - h*dx, phi - h*dphi, f3)
        call source_field(x - 2*h*dx, phi - 2*h*dphi, f4)
        d = (-f1 + 8*f2 - 8*f3 + f4)/(12*h)
    end function

    subroutine source_samples(nrho, ntheta, nphi, src)
        ! Gauss in rho, uniform theta and phi over one period (weights x P).
        integer, intent(in) :: nrho, ntheta, nphi
        type(fc_source_t), intent(out) :: src
        real(dp) :: xg(nrho), wg(nrho), rho, th, phi, x(2)
        integer :: i, l, m, k, p
        p = source_period
        src%period = p
        call gauss(nrho, xg, wg)
        allocate (src%xv(3, nrho*ntheta*nphi), src%wv(nrho*ntheta*nphi), &
            src%fv(3, nrho*ntheta*nphi), src%xs(3, ntheta*nphi), &
            src%ns(3, ntheta*nphi), src%fs(3, ntheta*nphi))
        k = 0
        do m = 1, nphi
            phi = 2*pi*(m - 1)/(nphi*p)
            do l = 1, ntheta
                th = 2*pi*(l - 1)/ntheta
                do i = 1, nrho
                    k = k + 1
                    rho = a*(xg(i) + 1)/2
                    x = [r0 + rho*cos(th), rho*sin(th)]
                    src%xv(:, k) = [x, phi]
                    src%wv(k) = x(1)*rho*(a/2)*wg(i)*(2*pi/ntheta)*(2*pi/nphi)
                    call source_field(x, phi, src%fv(:, k))
                end do
                x = [r0 + a*cos(th), a*sin(th)]
                i = (m - 1)*ntheta + l
                src%xs(:, i) = [x, phi]
                src%ns(:, i) = [cos(th), 0.0_dp, sin(th)]*x(1)*a &
                    *(2*pi/ntheta)*(2*pi/nphi)
                ! Outer trace of F inside V (q -> 0+ keeps only grad chi).
                call source_field([r0 + a*(1 - 1.0e-12_dp)*cos(th), &
                    a*(1 - 1.0e-12_dp)*sin(th)], phi, src%fs(:, i))
            end do
        end do
    end subroutine

    subroutine source_case(k, nr, tol)
        ! (1) fc_casing from samples of F against direct Biot-Savart of
        !     J = curl F/mu0 by the same polar quadrature of V, at every
        !     fourth boundary trace point (psi up to a constant);
        ! (2) the period-P solve against the period-1 solve with modes
        !     0..2P on the shared modes 0, P, 2P.
        integer, intent(in) :: k, nr
        real(dp), intent(in) :: tol
        integer, parameter :: nb = 9
        type(free_b_field_t) :: field, full
        type(fc_source_t) :: src
        type(fc_exterior_t) :: ext, cas
        type(fc_report_t) :: report, rfull
        real(dp), allocatable :: x(:), psi(:), bc(:, :), bs(:, :)
        real(dp) :: b0(3), b2(3), p(2), ep, scale, eext, dpsi
        integer :: st, i, l, n
        source_period = 2
        call disk_field(k, nr, 2, field)
        field%period = source_period
        call source_samples(16, 128, 96, src)
        call fc_casing(field, src, nb, cas, st)
        if (st /= 0) error stop 'casing'
        n = (size(cas%xy, 2) + 3)/4
        allocate (psi(n), bc(2, n), bs(2, n))
        !$omp parallel do schedule(dynamic)
        do i = 1, n
            call direct_exterior(cas%xy(:, 4*i - 3), nb, psi(i), bc(:, i), bs(:, i))
        end do
        !$omp end parallel do
        dpsi = sum(cas%psi(1::4) - psi)/n
        eext = max(maxval(abs(cas%bc(:, 1::4, 1) - bc)), &
            maxval(abs(cas%bs(:, 1::4, 1) - bs))) &
            /max(maxval(abs(bc)), maxval(abs(bs)))
        eext = max(eext, maxval(abs(cas%psi(1::4) - psi - dpsi))/maxval(abs(psi)))
        call fc_exterior(field, source_current, 64, ext, st)
        if (st /= 0) error stop 'exterior data'
        call fc_solve(field, source_current, ext, x, report, st, 64)
        if (st /= 0) error stop 'solve, period 2'
        call disk_field(k, nr, 2*source_period, full)
        call fc_exterior(full, source_current, 64, ext, st)
        if (st /= 0) error stop 'exterior data, period 1'
        call fc_solve(full, source_current, ext, x, rfull, st, 64)
        if (st /= 0) error stop 'solve, period 1'
        ep = 0; scale = 0
        do i = 1, 5
            do l = 0, 7
                p = [r0 + 0.11_dp*i*cos(2*pi*l/8 + 0.2_dp), &
                    0.11_dp*i*sin(2*pi*l/8 + 0.2_dp)]
                b0 = fe_b(field, p, 0.3_dp*l)
                b2 = fe_b(full, p, 0.3_dp*l)
                ep = max(ep, norm2(b2 - b0)); scale = max(scale, norm2(b0))
            end do
        end do
        print '(a,i0,a,i0)', 'source P', k, ' nr ', nr
        print '(a,es12.4)', '  casing vs direct Biot-Savart of J  ', eext
        print '(a,es12.4)', '  period 2 vs period 1               ', ep/scale
        print '(a,es12.4)', '  removed non-solenoidal fraction    ', report%projection
        if (eext > tol) error stop 'casing exterior data'
        if (ep/scale > 1.0e-10_dp) error stop 'period solve'
    end subroutine

    subroutine direct_exterior(xy, nb, psi, bc, bs)
        ! R A_phi (phi mean) and the cos/sin amplitudes of (B_R, B_Z) of
        ! mode P at (R, Z) = xy from mu0/(4 pi) sum W J/|d| and
        ! W J x d/|d|^3 over the polar samples of V on the whole torus.
        real(dp), intent(in) :: xy(2)
        integer, intent(in) :: nb
        real(dp), intent(out) :: psi, bc(2), bs(2)
        integer, parameter :: nrho = 16, ntheta = 128, nphi = 192
        real(dp) :: xg(nrho), wg(nrho), rho, th, ph, w, xs(2), j(3), jc(3)
        real(dp) :: d(3), x(3), b(3), aa(3), c, s, ri, ang
        complex(dp) :: bn(2)
        integer :: i, l, m, t, pp
        pp = source_period
        call gauss(nrho, xg, wg)
        psi = 0; bn = 0
        do t = 0, nb - 1
            ang = 2*pi*t/(nb*pp)
            x = [xy(1)*cos(ang), xy(1)*sin(ang), xy(2)]
            b = 0; aa = 0
            do m = 1, nphi
                ph = 2*pi*(m - 1)/nphi
                c = cos(ph); s = sin(ph)
                do l = 1, ntheta
                    th = 2*pi*(l - 1)/ntheta
                    do i = 1, nrho
                        rho = a*(xg(i) + 1)/2
                        xs = [r0 + rho*cos(th), rho*sin(th)]
                        w = xs(1)*rho*(a/2)*wg(i)*(2*pi/ntheta)*(2*pi/nphi)
                        call source_current(xs, ph, j)
                        jc = [j(1)*c - j(2)*s, j(1)*s + j(2)*c, j(3)]*w
                        d = x - [xs(1)*c, xs(1)*s, xs(2)]
                        ri = 1/norm2(d)
                        aa = aa + jc*ri
                        b = b + [jc(2)*d(3) - jc(3)*d(2), jc(3)*d(1) - jc(1)*d(3), &
                            jc(1)*d(2) - jc(2)*d(1)]*ri**3
                    end do
                end do
            end do
            b = mu*b/(4*pi); aa = mu*aa/(4*pi)
            psi = psi + xy(1)*(-aa(1)*sin(ang) + aa(2)*cos(ang))/nb
            bn = bn + [b(1)*cos(ang) + b(2)*sin(ang), b(3)] &
                *cmplx(cos(2*pi*t/nb), -sin(2*pi*t/nb), dp)/nb
        end do
        bc = 2*real(bn, dp); bs = -2*aimag(bn)
    end subroutine

    ! ------------------------------------------------------------- helical
    subroutine manufactured(x, phi, bf, jf, w)
        ! A = (0, w, U), U = s + t Theta(phi), mu0 = 1.
        real(dp), intent(in) :: x(2), phi
        real(dp), intent(out) :: bf(3), jf(3), w
        real(dp), parameter :: e1 = 0.3_dp, e2 = 0.2_dp
        real(dp) :: r, xx, z, d, sx, sz, sxx, sxz, t, th, th1, th2
        real(dp) :: ur, urr, urz, uphi, uphiphi, uphiz, wr, wz, wrr, wzz
        r = x(1); xx = r - r0; z = x(2)
        d = bdom**2 - xx**2 - z**2
        sx = -4*xx*d; sz = -4*z*d; sxx = -4*d + 8*xx**2; sxz = 8*xx*z
        t = 1 + 0.5_dp*xx + 0.3_dp*z
        th = e1*cos(phi) + e2*sin(2*phi)
        th1 = -e1*sin(phi) + 2*e2*cos(2*phi)
        th2 = -e1*cos(phi) - 4*e2*sin(2*phi)
        ur = sx + 0.5_dp*th; urr = sxx; urz = sxz
        uphi = t*th1; uphiphi = t*th2; uphiz = 0.3_dp*th1
        w = 0.5_dp + d*(1 + 0.2_dp*xx)
        wr = -2*xx*(1 + 0.2_dp*xx) + 0.2_dp*d
        wz = -2*z*(1 + 0.2_dp*xx)
        wrr = -2 - 1.2_dp*xx
        wzz = -2*(1 + 0.2_dp*xx)
        bf = [uphi/r - wz, -ur, wr + w/r]
        jf = [urz, uphiz/r - wzz - wrr - wr/r + w/r**2, -ur/r - urr - uphiphi/r**2]
        if (.false.) jf(1) = sz
    end subroutine

    subroutine manufactured_current(x, phi, j)
        real(dp), intent(in) :: x(2), phi
        real(dp), intent(out) :: j(3)
        real(dp) :: bf(3), w
        call manufactured(x, phi, bf, j, w)
    end subroutine

    subroutine helical_case(k)
        integer, intent(in) :: k
        type(free_b_field_t) :: field
        type(fc_exterior_t) :: ext
        type(fc_report_t) :: report
        real(dp), allocatable :: x(:)
        real(dp) :: err(3), bf(3), jf(3), w, p(2), b(3), phi, scale, rho, th
        complex(dp) :: modes(2, 2)
        integer :: level, st, i, l, m, n, jp
        mu = 1
        do level = 1, 3
            call disk_field(k, 2**(level + 2), 2, field)
            call meridian_boundary_trace_points(field%fs, ext%xy)
            allocate (ext%psi(size(ext%xy, 2)), ext%bc(2, size(ext%xy, 2), 2), &
                ext%bs(2, size(ext%xy, 2), 2))
            do i = 1, size(ext%xy, 2)
                modes = 0
                do jp = 0, 15
                    phi = 2*pi*jp/16
                    call manufactured(ext%xy(:, i), phi, bf, jf, w)
                    do n = 1, 2
                        modes(:, n) = modes(:, n) + bf([1, 3]) &
                            *cmplx(cos(n*phi), sin(n*phi), dp)/8
                    end do
                end do
                ext%psi(i) = ext%xy(1, i)*w
                ext%bc(:, i, :) = real(modes, dp); ext%bs(:, i, :) = aimag(modes)
            end do
            call fc_solve(field, manufactured_current, ext, x, report, st)
            if (st /= 0) error stop 'manufactured solve'
            err(level) = 0; scale = 0
            do i = 1, 5
                rho = 0.11_dp*i
                do l = 0, 6
                    th = 2*pi*l/7 + 0.2_dp
                    p = [r0 + rho*cos(th), rho*sin(th)]
                    do m = 0, 4
                        phi = 2*pi*m/5 + 0.1_dp
                        call manufactured(p, phi, bf, jf, w)
                        b = fe_b(field, p, phi)
                        err(level) = max(err(level), norm2(b - bf))
                        scale = max(scale, norm2(bf))
                    end do
                end do
            end do
            err(level) = err(level)/scale
            print '(a,i0,a,i0,a,es11.3,a,es11.3,a,es11.3)', 'helical P', k, &
                ' nr ', 2**(level + 2), ' max|B-B_ex|/max|B| ', err(level), &
                ' removed fraction ', report%projection, ' circulation ', &
                report%circulation
            deallocate (ext%xy, ext%psi, ext%bc, ext%bs)
        end do
        print '(a,2f8.3)', '  observed orders ', log(err(1:2)/err(2:3))/log(2.0_dp)
        if (log(err(2)/err(3))/log(2.0_dp) < k - 0.5_dp) error stop 'order'
        if (report%projection > 1.0e-6_dp) error stop 'projection'
    end subroutine

    ! ------------------------------------------------------------ exterior
    subroutine localized_current(x, phi, j)
        ! The ring current plus curl(gamma e_Z),
        ! gamma = h0 q^4 exp(0.5 cos phi) (1 + sin 2 phi)/2.
        real(dp), intent(in) :: x(2), phi
        real(dp), intent(out) :: j(3)
        real(dp) :: q, g, gphi, gr, jring(3)
        q = 1 - ((x(1) - r0)**2 + x(2)**2)/a**2
        j = 0
        if (q <= 0) return
        g = exp(0.5_dp*cos(phi))*(1 + sin(2*phi))/2
        gphi = exp(0.5_dp*cos(phi))*(-0.5_dp*sin(phi)*(1 + sin(2*phi)) &
            + 2*cos(2*phi))/2
        gr = h0*4*q**3*(-2*(x(1) - r0)/a**2)*g
        j = [h0*q**4*gphi/x(1), -gr, 0.0_dp]
        call ring_current(x, phi, jring)
        j = j + jring
    end subroutine

    subroutine exterior_case(k, nr, nmax, tol)
        integer, intent(in) :: k, nr, nmax
        real(dp), intent(in) :: tol
        type(free_b_field_t) :: field
        type(fc_exterior_t) :: ext
        type(fc_report_t) :: report
        real(dp), allocatable :: x(:), surf(:, :, :, :)
        real(dp), parameter :: dist(2) = [1.0_dp, 1.3_dp]
        real(dp) :: targets(3, 6), bfe(3, 6), bbs(3, 6), err, scale
        integer :: modes, st, i
        do i = 1, 6
            targets(:, i) = [r0 + dist(mod(i, 2) + 1)*cos(1.1_dp*i), &
                dist(mod(i, 2) + 1)*sin(1.1_dp*i), 0.7_dp*i]
        end do
        call direct_biot_savart(targets, bbs)
        print '(a,f6.3,a)', 'exterior: domain boundary ', bdom - a, &
            ' m outside the current support; targets 0.4 and 0.7 m outside the domain'
        do modes = 2, nmax, 2
            call disk_field(k, nr, modes, field)
            call fc_exterior(field, localized_current, 256, ext, st)
            if (st /= 0) error stop 'exterior data'
            call fc_solve(field, localized_current, ext, x, report, st, &
                nphi=32)
            if (st /= 0) error stop 'solve'
            call kirchhoff(field, targets, bfe)
            err = 0; scale = 0
            do i = 1, 6
                err = max(err, norm2(bfe(:, i) - bbs(:, i)))
                scale = max(scale, norm2(bbs(:, i)))
            end do
            print '(a,i0,a,i0,a,i0,a,es11.3,a,es11.3)', 'P', k, ' nr ', nr, &
                ' N ', modes, ' exterior max|dB|/max|B| ', err/scale, &
                ' removed fraction ', report%projection
        end do
        if (err/scale > tol) error stop 'exterior field'
        if (allocated(surf)) deallocate (surf)
    end subroutine

    subroutine kirchhoff(field, targets, bout)
        ! B(x) = (1/4 pi) oint [(n x B') x d + (n . B') d] / |d|^3 dS',
        ! d = x - x', for x outside the surface rho = bdom (inner trace B').
        type(free_b_field_t), intent(in) :: field
        real(dp), intent(in) :: targets(:, :)
        real(dp), intent(out) :: bout(:, :)
        integer, parameter :: nth = 256, nph = 256
        real(dp) :: th, ph, xs(2), b(3), nrm(3), bc(3), nc(3), xc(3), d(3)
        real(dp) :: ds, acc(3, size(targets, 2)), tc(3)
        integer :: i, l, it
        acc = 0
        !$omp parallel do default(shared) private(i, l, it, th, ph, xs, b, nrm, &
        !$omp& bc, nc, xc, d, ds, tc) reduction(+:acc) collapse(2)
        do i = 1, nth
            do l = 1, nph
                th = 2*pi*(i - 1)/nth; ph = 2*pi*(l - 1)/nph
                xs = [r0 + bdom*(1 - 1.0e-10_dp)*cos(th), &
                    bdom*(1 - 1.0e-10_dp)*sin(th)]
                b = fe_b(field, xs, ph)
                nrm = [cos(th), 0.0_dp, sin(th)]
                bc = [b(1)*cos(ph) - b(2)*sin(ph), b(1)*sin(ph) + b(2)*cos(ph), b(3)]
                nc = [nrm(1)*cos(ph), nrm(1)*sin(ph), nrm(3)]
                xc = [xs(1)*cos(ph), xs(1)*sin(ph), xs(2)]
                ds = xs(1)*bdom*(2*pi/nth)*(2*pi/nph)
                do it = 1, size(targets, 2)
                    tc = [targets(1, it)*cos(targets(3, it)), &
                        targets(1, it)*sin(targets(3, it)), targets(2, it)]
                    d = tc - xc
                    acc(:, it) = acc(:, it) + ds/(4*pi)*(cross(cross(nc, bc), d) &
                        + dot_product(nc, bc)*d)/norm2(d)**3
                end do
            end do
        end do
        !$omp end parallel do
        do it = 1, size(targets, 2)
            call to_cylindrical(acc(:, it), targets(3, it), bout(:, it))
        end do
    end subroutine

    subroutine direct_biot_savart(targets, bout)
        ! Independent 3D quadrature of J over the support rho' < a.
        real(dp), intent(in) :: targets(:, :)
        real(dp), intent(out) :: bout(:, :)
        integer, parameter :: nrho = 32, nth = 96, nph = 256
        real(dp) :: xg(nrho), wg(nrho), rho, th, ph, p(2), j(3), jc(3), xc(3)
        real(dp) :: d(3), tc(3), dv, acc(3, size(targets, 2))
        integer :: i, l, m, it
        call gauss(nrho, xg, wg)
        acc = 0
        !$omp parallel do default(shared) private(i, l, m, it, rho, th, ph, p, j, &
        !$omp& jc, xc, d, tc, dv) reduction(+:acc) collapse(2)
        do i = 1, nrho
            do l = 1, nth
                rho = a*(xg(i) + 1)/2
                th = 2*pi*(l - 1)/nth
                p = [r0 + rho*cos(th), rho*sin(th)]
                do m = 1, nph
                    ph = 2*pi*(m - 1)/nph
                    call localized_current(p, ph, j)
                    jc = [j(1)*cos(ph) - j(2)*sin(ph), j(1)*sin(ph) + j(2)*cos(ph), &
                        j(3)]
                    xc = [p(1)*cos(ph), p(1)*sin(ph), p(2)]
                    dv = rho*(a/2)*wg(i)*(2*pi/nth)*p(1)*(2*pi/nph)
                    do it = 1, size(targets, 2)
                        tc = [targets(1, it)*cos(targets(3, it)), &
                            targets(1, it)*sin(targets(3, it)), targets(2, it)]
                        d = tc - xc
                        acc(:, it) = acc(:, it) + mu/(4*pi)*dv*cross(jc, d) &
                            /norm2(d)**3
                    end do
                end do
            end do
        end do
        !$omp end parallel do
        do it = 1, size(targets, 2)
            call to_cylindrical(acc(:, it), targets(3, it), bout(:, it))
        end do
    end subroutine

    pure function cross(u, v) result(w)
        real(dp), intent(in) :: u(3), v(3)
        real(dp) :: w(3)
        w = [u(2)*v(3) - u(3)*v(2), u(3)*v(1) - u(1)*v(3), u(1)*v(2) - u(2)*v(1)]
    end function

    pure subroutine to_cylindrical(bc, ph, b)
        real(dp), intent(in) :: bc(3), ph
        real(dp), intent(out) :: b(3)
        b = [bc(1)*cos(ph) + bc(2)*sin(ph), -bc(1)*sin(ph) + bc(2)*cos(ph), bc(3)]
    end subroutine
end program
