program test_vacuum_surfaces
    !! L3 oracles on an analytic circular equilibrium psi = c ((R-R0)^2 + Z^2):
    !! exact q, Hamada and PEST angles, single-harmonic recovery in both angles,
    !! G-EQDSK round trip, and resolution convergence / thread invariance of the
    !! vacuum b_mn of synthetic 3D coils.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_coils_3d
    use kin6d_vacuum_surfaces
!$  use omp_lib
    implicit none
    real(dp), parameter :: pi = 3.14159265358979323846_dp
    real(dp), parameter :: r0 = 6.5_dp, c = 0.8_dp, amax = 1.8_dp, fconst = 33.0_dp
    real(dp), parameter :: ip = 1.5e7_dp
    integer :: failures = 0
    type(psi_grid_t) :: g

    call make_grid(g)
    call test_geometry()
    call test_single_harmonic()
    call test_eqdsk_round_trip()
    call test_resolution_and_threads()
    if (failures > 0) then
        print '(a,i0,a)', 'test_vacuum_surfaces: ', failures, ' failure(s)'
        error stop 1
    end if
    print '(a)', 'test_vacuum_surfaces: all passed'

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

    subroutine grid_arrays(psirz, f)
        real(dp), allocatable, intent(out) :: psirz(:, :), f(:)
        integer, parameter :: nr = 201, nz = 201
        integer :: i, j
        allocate (psirz(nr, nz), f(nr))
        do j = 1, nz
            do i = 1, nr
                psirz(i, j) = c*((4.5_dp + 4.0_dp*(i - 1)/(nr - 1) - r0)**2 &
                                 + (-2.0_dp + 4.0_dp*(j - 1)/(nz - 1))**2)
            end do
        end do
        f = fconst
    end subroutine

    subroutine make_grid(gg)
        type(psi_grid_t), intent(out) :: gg
        real(dp), allocatable :: psirz(:, :), f(:)
        integer :: st
        call grid_arrays(psirz, f)
       call psi_grid_init(gg, psirz, 4.5_dp, 4.0_dp, 0.0_dp, 4.0_dp, 6.45_dp, 0.05_dp, &
                           c*amax**2, ip, f, st)
        if (st /= 0) error stop 'grid'
    end subroutine

    subroutine test_geometry()
        type(surface_t) :: s
        real(dp) :: r, kappa, qex, e_q, e_h, e_p, th, hp, pp
        integer :: j, st
       call check('axis position error', abs(g%r_axis - r0) + abs(g%z_axis), 1.0e-12_dp)
        call surface_build(g, 0.5_dp, 96, s, st)
        r = amax*sqrt(0.5_dp)
        qex = -fconst/(2*c*sqrt(r0**2 - r**2))
        e_q = abs(s%q - qex)/abs(qex)
        call check('signed q vs closed form', e_q, 1.0e-10_dp)
        kappa = sqrt((r0 - r)/(r0 + r))
        e_h = 0; e_p = 0
        do j = 1, s%n
            th = s%theta_g(j)
            hp = (r0*th + r*sin(th))/r0
            pp = 2*atan2(kappa*sin(th/2), cos(th/2))
            e_h = max(e_h, abs(s%theta_s(j, angle_hamada) - hp))
            e_p = max(e_p, abs(s%theta_s(j, angle_pest) - pp))
        end do
        call check('Hamada angle vs closed form', e_h, 1.0e-10_dp)
        call check('PEST angle vs closed form', e_p, 1.0e-10_dp)
    end subroutine

    subroutine test_single_harmonic()
        type(surface_t) :: s
        complex(dp), allocatable :: b(:, :)
        complex(dp) :: bmn(-20:20), phimn(-20:20)
        real(dp) :: err
        integer :: st, j, angle, m0, n, ia
        n = 3
        m0 = -10
        call surface_build(g, 0.8_dp, 128, s, st)
        allocate (b(2, s%n))
        err = 0
        do ia = 1, 2
            angle = ia
            do j = 1, s%n
                ! normal component 1/2 exp(-i m0 theta_s) exp(-i n (nu - q theta_s))
                b(1, j) = 0.5_dp*exp(cmplx(0.0_dp, -m0*s%theta_s(j, angle) &
                                   - n*(s%nu(j) - s%q*s%theta_s(j, angle)), dp))*s%nr(j)
                b(2, j) = 0.5_dp*exp(cmplx(0.0_dp, -m0*s%theta_s(j, angle) &
                                   - n*(s%nu(j) - s%q*s%theta_s(j, angle)), dp))*s%nz(j)
            end do
            call surface_harmonics(s, n, b, angle, -20, 20, bmn, phimn)
            err = max(err, abs(bmn(m0) - 1))
            bmn(m0) = 0
            err = max(err, maxval(abs(bmn)))
        end do
       call check('single harmonic recovered in PEST and Hamada angle', err, 1.0e-12_dp)
    end subroutine

    subroutine test_eqdsk_round_trip()
        real(dp), allocatable :: psirz(:, :), f(:)
        type(psi_grid_t) :: h
        integer :: unit, st, i, j
        real(dp) :: p1, p2, d1(2), d2(2), z1, z2
        call grid_arrays(psirz, f)
        open (newunit=unit, file='test_vacuum_surfaces.eqdsk', status='replace')
        write (unit, '(a,3i5)') ' TEST  00/00/00  #000000,00000ms20050328', 3, &
            size(psirz, 1), size(psirz, 2)
        write (unit, '(5es16.9)') 4.0_dp, 4.0_dp, r0, 4.5_dp, 0.0_dp
        write (unit, '(5es16.9)') r0, 0.0_dp, 0.0_dp, c*amax**2, 5.3_dp
        write (unit, '(5es16.9)') ip, 0.0_dp, 0.0_dp, r0, 0.0_dp
        write (unit, '(5es16.9)') 0.0_dp, 0.0_dp, c*amax**2, 0.0_dp, 0.0_dp
        write (unit, '(5es16.9)') (f(i), i=1, size(f))
        do j = 1, 3
            write (unit, '(5es16.9)') (0.0_dp, i=1, size(f))
        end do
     write (unit, '(5es16.9)') ((psirz(i, j), i=1, size(psirz, 1)), j=1, size(psirz, 2))
        close (unit, status='keep')
        call eqdsk_read('test_vacuum_surfaces.eqdsk', h, st)
        call check('EQDSK read status', real(st, dp), 0.0_dp)
        call psi_eval(g, 7.3_dp, 0.4_dp, p1, d1(1), d1(2))
        call psi_eval(h, 7.3_dp, 0.4_dp, p2, d2(1), d2(2))
        z1 = abs(p1 - p2) + norm2(d1 - d2)
        call check('EQDSK round trip psi and gradient', z1, 1.0e-6_dp)
        z2 = abs(p1 - c*((7.3_dp - r0)**2 + 0.4_dp**2))
        call check('spline psi vs analytic', z2, 1.0e-12_dp)
    end subroutine

    subroutine synthetic_coils(set, w)
        type(coil_set_t), intent(out) :: set
        real(dp), intent(out) :: w(9)
        real(dp) :: base(3, 16), pts(3, 16*9), ang, cc, ss
        integer :: i, k, st, npts(9)
        do i = 1, 16
            ang = 2*pi*(i - 1)/16
            base(:, i) = [r0 + 2.7_dp + 0.4_dp*cos(ang), 0.35_dp*sin(ang), &
                          1.2_dp*sin(ang)*0.5_dp + 0.3_dp*cos(ang)]
        end do
        do k = 1, 9
            ang = 2*pi*(k - 1)/9
            cc = cos(ang); ss = sin(ang)
            do i = 1, 16
                pts(:, 16*(k - 1) + i) = [cc*base(1, i) - ss*base(2, i), &
                                          ss*base(1, i) + cc*base(2, i), base(3, i)]
            end do
            npts(k) = 16
            w(k) = 9.0e4_dp*cos(3*ang + 0.4_dp)
        end do
        call coil_set_from_polylines(set, npts, pts, status=st)
        if (st /= 0) error stop 'coils'
    end subroutine

    subroutine vacuum_bmn(set, w, psin, nth, nphi, n, bmn, angle)
        type(coil_set_t), intent(in) :: set
        real(dp), intent(in) :: w(:), psin
        integer, intent(in) :: nth, nphi, n, angle
        complex(dp), intent(out) :: bmn(-24:24)
        type(surface_t) :: s
        complex(dp) :: bn(3, 1, nth), bb(2, nth), phimn(-24:24)
        real(dp), allocatable :: rz(:, :)
        integer :: st, j
        call surface_build(g, psin, nth, s, st)
        allocate (rz(2, nth))
        do j = 1, nth
            rz(:, j) = [s%r(j), s%z(j)]
        end do
        call coil_set_modes(set, w, rz, nphi, [n], bn, st)
        do j = 1, nth
            bb(1, j) = bn(1, 1, j)
            bb(2, j) = bn(3, 1, j)
        end do
        call surface_harmonics(s, n, bb, angle, -24, 24, bmn, phimn)
    end subroutine

    subroutine test_resolution_and_threads()
        type(coil_set_t) :: set
        real(dp) :: w(9), scale, e1, e2, e3
        complex(dp) :: ref(-24:24), a(-24:24), b(-24:24), cc(-24:24), t1(-24:24)
        complex(dp) :: t4(-24:24)
        call synthetic_coils(set, w)
        call vacuum_bmn(set, w, 0.8_dp, 128, 384, 3, ref, angle_pest)
        scale = maxval(abs(ref))
        call vacuum_bmn(set, w, 0.8_dp, 32, 384, 3, a, angle_pest)
        call vacuum_bmn(set, w, 0.8_dp, 48, 384, 3, b, angle_pest)
        call vacuum_bmn(set, w, 0.8_dp, 64, 384, 3, cc, angle_pest)
        e1 = maxval(abs(a - ref))/scale
        e2 = maxval(abs(b - ref))/scale
        e3 = maxval(abs(cc - ref))/scale
      print '(a,3es12.3)', 'b_mn error vs ntheta=128, ntheta = 32, 48, 64: ', e1, e2, e3
        call check('b_mn error ntheta=48 below ntheta=32 (ratio)', e2/e1, 1.0e-3_dp)
        call check('b_mn converged at ntheta = 64 (relative)', e3, 1.0e-8_dp)
        call vacuum_bmn(set, w, 0.8_dp, 128, 96, 3, a, angle_pest)
        call vacuum_bmn(set, w, 0.8_dp, 128, 192, 3, b, angle_pest)
        e1 = maxval(abs(a - ref))/scale
        e2 = maxval(abs(b - ref))/scale
        print '(a,2es12.3)', 'b_mn error vs nphi=384, nphi = 96, 192: ', e1, e2
        call check('b_mn error nphi=192 below nphi=96 (ratio)', e2/e1, 0.01_dp)
        call check('b_mn converged at nphi = 192 (relative)', e2, 1.0e-6_dp)
        call vacuum_bmn(set, w, 0.8_dp, 128, 192, 3, b, angle_hamada)
        call vacuum_bmn(set, w, 0.8_dp, 256, 192, 3, cc, angle_hamada)
        call check('Hamada b_mn converged ntheta 128 vs 256 (relative)', &
                   maxval(abs(b - cc))/maxval(abs(cc)), 1.0e-6_dp)
!$      call omp_set_num_threads(1)
        call vacuum_bmn(set, w, 0.8_dp, 128, 192, 3, t1, angle_pest)
!$      call omp_set_num_threads(4)
        call vacuum_bmn(set, w, 0.8_dp, 128, 192, 3, t4, angle_pest)
        call check('b_mn thread invariance max|diff|', maxval(abs(t1 - t4)), 0.0_dp)
    end subroutine

end program
