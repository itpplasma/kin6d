program test_coils_3d
    !! K-C1 oracles: segment Biot-Savart against quadrature, polygon ring against
    !! the on-axis ring formula, div/curl of the field, toroidal selection rules,
    !! Fourier synthesis against the direct 3D field and thread invariance.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_coils_3d
!$  use omp_lib
    implicit none
    real(dp), parameter :: pi = 3.14159265358979323846_dp
    integer :: failures = 0

    call test_segment_quadrature()
    call test_ring_on_axis()
    call test_div_curl()
    call test_selection_and_synthesis()
    if (failures > 0) then
        print '(a,i0,a)', 'test_coils_3d: ', failures, ' failure(s)'
        error stop 1
    end if
    print '(a)', 'test_coils_3d: all passed'

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

    subroutine gauss_legendre(n, x, w)
        integer, intent(in) :: n
        real(dp), intent(out) :: x(n), w(n)
        integer :: i, it
        real(dp) :: z, p1, p2, p3, pp
        integer :: j
        do i = 1, n
            z = cos(pi*(i - 0.25_dp)/(n + 0.5_dp))
            do it = 1, 100
                p1 = 1; p2 = 0
                do j = 1, n
                    p3 = p2; p2 = p1
                    p1 = ((2*j - 1)*z*p2 - (j - 1)*p3)/j
                end do
                pp = n*(z*p1 - p2)/(z*z - 1)
                z = z - p1/pp
                if (abs(p1/pp) < 1.0e-16_dp) exit
            end do
            x(i) = z
            w(i) = 2/((1 - z*z)*pp*pp)
        end do
    end subroutine

    subroutine test_segment_quadrature()
        real(dp) :: x(64), w(64), a(3), l(3), p(3), b(3), bq(3), r(3), err, scale
        integer :: i, trial
        call gauss_legendre(64, x, w)
        err = 0
        scale = 0
        do trial = 1, 5
            a = [1.0_dp, -0.4_dp, 0.2_dp*trial]
            l = [0.3_dp, 0.8_dp, -0.5_dp + 0.1_dp*trial]
            p = [0.2_dp*trial, 0.9_dp, 1.5_dp - 0.3_dp*trial]
            bq = 0
            do i = 1, 64
                r = p - (a + l*(x(i) + 1)/2)
                bq = bq + w(i)/2*cross3(l, r)/norm2(r)**3
            end do
            bq = bq*coil_mu0/(4*pi)
            b = segment_b(a, l, p)
            err = max(err, norm2(b - bq))
            scale = max(scale, norm2(bq))
        end do
       call check('segment vs 64-pt Gauss quadrature (relative)', err/scale, 1.0e-13_dp)
    end subroutine

    pure function cross3(u, v) result(w)
        real(dp), intent(in) :: u(3), v(3)
        real(dp) :: w(3)
        w = [u(2)*v(3) - u(3)*v(2), u(3)*v(1) - u(1)*v(3), u(1)*v(2) - u(2)*v(1)]
    end function

    subroutine ring_set(set, n, rho, z0)
        type(coil_set_t), intent(out) :: set
        integer, intent(in) :: n
        real(dp), intent(in) :: rho, z0
        real(dp) :: pts(3, n)
        integer :: i, st
        do i = 1, n
            pts(:, i) = [rho*cos(2*pi*(i - 1)/n), rho*sin(2*pi*(i - 1)/n), z0]
        end do
        call coil_set_from_polylines(set, [n], pts, status=st)
        if (st /= 0) error stop 'ring'
    end subroutine

    subroutine test_ring_on_axis()
        type(coil_set_t) :: set
        real(dp) :: b(3), exact, e1, e2, p(3)
        real(dp), parameter :: rho = 2.0_dp, z = 1.3_dp
        p = [0.0_dp, 0.0_dp, 0.0_dp]
        exact = coil_mu0*rho**2/(2*(rho**2 + z**2)**1.5_dp)
        call ring_set(set, 64, rho, z)
        call coil_set_b_cart(set, [1.0_dp], p, b)
        e1 = abs(b(3) - exact)/exact
        call ring_set(set, 128, rho, z)
        call coil_set_b_cart(set, [1.0_dp], p, b)
        e2 = abs(b(3) - exact)/exact
        call check('polygon ring on axis, N=128', e2, 2.0e-4_dp)
        ! inscribed N-gon: error is O(N^-2)
        call check('polygon ring N^-2 convergence (4 - e1/e2 deviation)', &
                   abs(e1/e2 - 4), 0.05_dp)
    end subroutine

    subroutine tilted_coil_set(set, ncopy)
        !! ncopy exact 40-degree rotated copies of one nonplanar loop near R = 8.
        type(coil_set_t), intent(out) :: set
        integer, intent(in) :: ncopy
        real(dp) :: base(3, 12), pts(3, 12*ncopy), c, s, ang
        integer :: i, k, st, npts(ncopy)
        do i = 1, 12
            ang = 2*pi*(i - 1)/12
            base(:, i) = [8.0_dp + 0.45_dp*cos(ang) + 0.05_dp*sin(2*ang), &
                          0.35_dp*sin(ang), 2.5_dp + 0.4_dp*sin(ang) + 0.1_dp*cos(ang)]
        end do
        do k = 1, ncopy
            ang = 2*pi*(k - 1)/ncopy
            c = cos(ang); s = sin(ang)
            do i = 1, 12
                pts(:, 12*(k - 1) + i) = [c*base(1, i) - s*base(2, i), &
                                          s*base(1, i) + c*base(2, i), base(3, i)]
            end do
            npts(k) = 12
        end do
        call coil_set_from_polylines(set, npts, pts, status=st)
        if (st /= 0) error stop 'tilted'
    end subroutine

    subroutine test_div_curl()
        type(coil_set_t) :: set
        real(dp) :: w(9), p(3), h, g(3, 3), bp(3), bm(3), bmag, divb, curlb(3)
        real(dp), parameter :: e(3, 3) = reshape([1, 0, 0, 0, 1, 0, 0, 0, 1], [3, 3])
        integer :: i
        call tilted_coil_set(set, 9)
        w = 1.0e4_dp
        p = [6.1_dp, 1.7_dp, 0.8_dp]
        h = 1.0e-3_dp
        do i = 1, 3
            call coil_set_b_cart(set, w, p + h*e(:, i), bp)
            call coil_set_b_cart(set, w, p - h*e(:, i), bm)
            g(:, i) = (bp - bm)/(2*h)         ! g(c, i) = dB_c / dx_i
        end do
        call coil_set_b_cart(set, w, p, bp)
        bmag = norm2(bp)
        divb = g(1, 1) + g(2, 2) + g(3, 3)
        curlb = [g(3, 2) - g(2, 3), g(1, 3) - g(3, 1), g(2, 1) - g(1, 2)]
        call check('div B / (|B| / 1 m)', abs(divb)/bmag, 1.0e-6_dp)
        call check('|curl B| / (|B| / 1 m)', norm2(curlb)/bmag, 1.0e-6_dp)
    end subroutine

    subroutine test_selection_and_synthesis()
        integer, parameter :: nphi = 256, nmode = 128
        type(coil_set_t) :: set
       real(dp) :: w(9), rz(2, 3), amp(0:nmode - 1), off, on, phi, b3(3), bsyn(3), bd(3)
        real(dp) :: err, bmag, cyl(3)
        complex(dp) :: bn(3, nmode, 3), bn1(3, nmode, 3)
        integer :: modes(nmode), k, ip, im, st, i
        real(dp) :: i0, th(9)
        call tilted_coil_set(set, 9)
        i0 = 9.0e4_dp
        do k = 1, 9
            th(k) = 2*pi*(k - 1)/9
            w(k) = i0*cos(3*th(k) + 0.3_dp)
        end do
        modes = [(i, i=0, nmode - 1)]
        rz = reshape([6.0_dp, 0.5_dp, 6.4_dp, -0.3_dp, 5.2_dp, 1.0_dp], [2, 3])
        call coil_set_modes(set, w, rz, nphi, modes, bn, st)
        if (st /= 0) error stop 'modes'
        ! selection rule: only n = 3, 6 (mod 9) survive for a cos(3 phi) pattern
        off = 0; on = 0
        do ip = 1, 3
            do im = 1, nmode
                amp(im - 1) = maxval(abs(bn(:, im, ip)))
            end do
            do im = 0, nmode - 1
                if (mod(im, 9) == 3 .or. mod(im, 9) == 6) then
                    on = max(on, amp(im)/maxval(amp))
                else
                    off = max(off, amp(im)/maxval(amp))
                end if
            end do
        end do
        call check('non-selected n amplitude / max (roundoff)', off, 1.0e-13_dp)
        call check('selected modes are populated (inverted: 1 if max<0.5)', &
                   merge(0.0_dp, 1.0_dp, on > 0.5_dp), 0.5_dp)
        call check('n=3,6,12 nonzero (inverted)', merge(0.0_dp, 1.0_dp, &
                           minval([maxval(abs(bn(:, 4, 1))), maxval(abs(bn(:, 7, 1))), &
                                       maxval(abs(bn(:, 13, 1)))]) > 1.0e-9_dp), 0.5_dp)
        ! synthesis against the direct field at off-grid angles
        err = 0
        do ip = 1, 3
            do k = 1, 7
                phi = 0.1_dp + 0.9_dp*k
                bsyn = coil_synthesis(bn(:, :, ip), modes, phi)
                call coil_set_b_cart(set, w, [rz(1, ip)*cos(phi), rz(1, ip)*sin(phi), &
                                              rz(2, ip)], b3)
        cyl = [b3(1)*cos(phi) + b3(2)*sin(phi), -b3(1)*sin(phi) + b3(2)*cos(phi), b3(3)]
                bd = cyl
                bmag = norm2(bd)
                err = max(err, norm2(bsyn - bd)/bmag)
            end do
        end do
        call check('Fourier synthesis vs direct Biot-Savart (relative)', err, 1.0e-8_dp)
        ! thread invariance
!$      call omp_set_num_threads(1)
        call coil_set_modes(set, w, rz, nphi, modes, bn, st)
!$      call omp_set_num_threads(4)
        call coil_set_modes(set, w, rz, nphi, modes, bn1, st)
        call check('thread invariance max|diff| (1 vs 4 threads)', &
                   maxval(abs(bn - bn1)), 0.0_dp)
    end subroutine

end program
