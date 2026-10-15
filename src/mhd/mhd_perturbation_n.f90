module kin6d_mhd_perturbation_n
    !! L7 core: displacement, Lagrangian field and curvature coupling of a per-n
    !! linear response on flux surfaces of an axisymmetric equilibrium.
    !!
    !! Input per surface: the Eulerian field dB_E of the response (physical
    !! field Re[dB_E e^{i n phi}] times two, i.e. dB_E is the coefficient of
    !! e^{i n phi}; cylindrical components (R, phi, Z)) at np points that are
    !! uniform in the periodic geometric angle theta_g of the surface, and the
    !! axisymmetric equilibrium B0, grad psi and the gradient tensor of B0 there.
    !!
    !! Ideal reconstruction (xi_par = 0 gauge).  With xi = xi_n nhat + xi_b ebin,
    !! nhat = sgn grad psi/|grad psi| (outward), ebin = b x nhat,
    !! dB = curl(xi x B0) = (B0.grad) xi - (xi.grad) B0 - B0 div xi, and
    !!   dB.grad psi        = B0.grad (xi.grad psi)                   (normal)
    !!   ebin.dB            = ebin.[(B0.grad) xi - (xi.grad) B0]      (binormal)
    !! (div xi multiplies B0 and drops out of both).  Both are first-order
    !! ODEs along the field line, B0.grad f = Bth d_theta f + i n (B_phi/R) f,
    !! Bth = B0.grad theta_g, solved with a periodic Fourier collocation matrix;
    !! they are singular only at n q = integer, which is reported as
    !! dist_rational and refused below a tolerance.  (B0.grad) and (xi.grad)
    !! act on cylindrical components with the basis-rotation terms +-v_phi/R.
    !!
    !! Lagrangian field: dB_L = dB_E + (xi.grad) B0, |B|_L = b.dB_E + xi.grad|B|
    !! = -B (kappa.xi + div xi_perp).  The parallel displacement is gauge:
    !! xi_par shifts |B|_L by xi_par b.grad|B|; we set xi_par = 0.
    !!
    !! Charts: spectra use the SFL angles of kin6d_vacuum_surfaces (PEST,
    !! Hamada) in the convention of `surface_harmonics`: X e^{i n phi} =
    !! e^{i n zeta} sum_m X_m e^{-i m theta_s}, zeta = phi - (nu - q theta_s).
    !! Boozer angles are not provided by the core yet (no D7 Boozer map).
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64, int8
    use kin6d_vacuum_surfaces, only: psi_grid_t, surface_t, angle_pest, &
        angle_hamada
    use kin6d_mhd_ideal_n, only: equilibrium_point, ideal_options_t
    implicit none
    private

    real(dp), parameter :: pi = 3.14159265358979323846_dp
    integer, parameter, public :: nq_spec = 6
    character(len=10), parameter, public :: spec_names(nq_spec) = &
        [character(len=10) :: 'xi_n', 'xi_b', 'dB_n_E', 'dB_par_E', 'dBmag_L', &
         'kappa_xi']

    type, public :: perturbation_t
        integer :: n = 0, np = 0
        real(dp) :: psin = 0.0_dp, q_eff = 0.0_dp, dist_rational = 0.0_dp
        real(dp) :: sgn_dpsi = 1.0_dp
        real(dp), allocatable :: theta_g(:), x(:, :)       !! (2, np) R, Z
        real(dp), allocatable :: b0(:, :), bmag(:)         !! (3, np) R,phi,Z
        real(dp), allocatable :: gradb(:, :)               !! (2, np) d_R, d_Z |B|
        real(dp), allocatable :: nhat(:, :), ebin(:, :)    !! (3, np)
        real(dp), allocatable :: kappa(:, :)               !! (3, np) (b.grad) b
        complex(dp), allocatable :: dbe(:, :), dbl(:, :)   !! (3, np)
        complex(dp), allocatable :: xi(:, :)               !! (3, np) cylindrical
        complex(dp), allocatable :: xi_psi(:), xi_n(:), xi_b(:)
        complex(dp), allocatable :: dbpar_e(:), dbmag_l(:), kappa_xi(:)
        !! spectra(m, quantity, chart), chart 1 = PEST, 2 = Hamada
        integer :: mmin = 0, mmax = -1
        complex(dp), allocatable :: spec(:, :, :)
    end type perturbation_t

    public :: perturbation_solve, perturbation_on_surface
    public :: equilibrium_gradient, perturbation_spectrum
    public :: perturbation_synthesis, perturbation_chart_map
    public :: perturbation_add_spectra, perturbation_export
    public :: sha256_file, sha256_verify

contains

    ! ------------------------------------------------------------ equilibrium

    subroutine equilibrium_gradient(g, opt, r, z, b0, gpsi, gb)
        !! B0 (R, phi, Z), grad psi (R, Z) and gb(i, j) = d_j B0_i (j = R, Z)
        !! by central differences of the exact point evaluation.
        type(psi_grid_t), intent(in) :: g
        type(ideal_options_t), intent(in) :: opt
        real(dp), intent(in) :: r, z
        real(dp), intent(out) :: b0(3), gpsi(2), gb(3, 2)
        real(dp) :: bp(3), bm(3), jj(3), gg(2), pp, h

        call equilibrium_point(g, opt, r, z, b0, jj, gpsi, pp)
        h = 1.0e-5_dp*r
        call equilibrium_point(g, opt, r + h, z, bp, jj, gg, pp)
        call equilibrium_point(g, opt, r - h, z, bm, jj, gg, pp)
        gb(:, 1) = (bp - bm)/(2*h)
        call equilibrium_point(g, opt, r, z + h, bp, jj, gg, pp)
        call equilibrium_point(g, opt, r, z - h, bm, jj, gg, pp)
        gb(:, 2) = (bp - bm)/(2*h)
    end subroutine equilibrium_gradient

    subroutine perturbation_on_surface(g, opt, s, n, dbe, res, status, tol_res)
        !! Reconstruct the perturbation on the surface s (surface_build) from
        !! dbe(3, s%n) = total Eulerian dB_E at the surface points.
        type(psi_grid_t), intent(in) :: g
        type(ideal_options_t), intent(in) :: opt
        type(surface_t), intent(in) :: s
        integer, intent(in) :: n
        complex(dp), intent(in) :: dbe(:, :)
        type(perturbation_t), intent(out) :: res
        integer, intent(out) :: status
        real(dp), intent(in), optional :: tol_res
        real(dp), allocatable :: x(:, :), b0(:, :), gpsi(:, :), gb(:, :, :)
        integer :: k, np

        np = s%n
        allocate (x(2, np), b0(3, np), gpsi(2, np), gb(3, 2, np))
        do k = 1, np
            x(:, k) = [s%r(k), s%z(k)]
            call equilibrium_gradient(g, opt, s%r(k), s%z(k), b0(:, k), &
                gpsi(:, k), gb(:, :, k))
        end do
        call perturbation_solve(n, g%sgn_dpsi, x, b0, gpsi, gb, dbe, res, status, &
            tol_res)
        res%psin = s%psin
        res%theta_g = s%theta_g
    end subroutine perturbation_on_surface

    ! ------------------------------------------------------------------- core

    subroutine fourier_diff_matrix(np, d)
        !! Periodic spectral first-derivative matrix on np uniform nodes of
        !! [0, 2 pi).
        integer, intent(in) :: np
        real(dp), intent(out) :: d(np, np)
        integer :: i, j, k
        real(dp) :: h2

        h2 = pi/np
        d = 0.0_dp
        do j = 1, np
            do i = 1, np
                if (i == j) cycle
                k = i - j
                if (mod(np, 2) == 0) then
                    d(i, j) = 0.5_dp*(-1)**k/tan(k*h2)
                else
                    d(i, j) = 0.5_dp*(-1)**k/sin(k*h2)
                end if
            end do
        end do
    end subroutine fourier_diff_matrix

    subroutine perturbation_solve(n, sgn, x, b0, gpsi, gb, dbe, res, status, tol_res)
        !! Core reconstruction; see the module header.  Status 0 ok, 1 bad
        !! input, 2 within tol_res (default 1e-3) of a rational n q = m (no
        !! displacement is computed), 3 singular collocation system.
        integer, intent(in) :: n
        real(dp), intent(in) :: sgn, x(:, :), b0(:, :), gpsi(:, :), gb(:, :, :)
        complex(dp), intent(in) :: dbe(:, :)
        type(perturbation_t), intent(inout) :: res
        integer, intent(out) :: status
        real(dp), intent(in), optional :: tol_res
        real(dp), allocatable :: d(:, :), dr(:), dz(:), bth(:), nh(:, :), &
            eb(:, :), bn(:, :), cc(:), aa(:), dnr(:), dnz(:), bb(:), bhat(:, :)
        complex(dp), allocatable :: amat(:, :), rhs(:), rhs2(:)
        complex(dp) :: v(3), w(3), vb(3)
        real(dp) :: tol, gn, bv(3), kk(3), t(3), pr
        integer, allocatable :: ipiv(:)
        integer :: np, k, j, info

        status = 1
        np = size(x, 2)
        if (np < 8 .or. size(b0, 2) /= np .or. size(dbe, 2) /= np .or. n == 0) return
        tol = 1.0e-3_dp
        if (present(tol_res)) tol = tol_res
        res%n = n
        res%np = np
        res%sgn_dpsi = sgn
        allocate (d(np, np), dr(np), dz(np), bth(np), nh(3, np), eb(3, np), &
            bb(np), bhat(3, np), aa(np), cc(np), dnr(np), dnz(np), bn(3, np))
        call fourier_diff_matrix(np, d)
        dr = matmul(d, x(1, :))
        dz = matmul(d, x(2, :))
        bth = (b0(1, :)*dr + b0(3, :)*dz)/(dr**2 + dz**2)
        res%q_eff = sum(b0(2, :)/x(1, :)/bth)/np
        res%dist_rational = abs(n*res%q_eff - nint(n*res%q_eff))
        res%x = x
        res%b0 = b0
        res%dbe = dbe
        if (.not. allocated(res%theta_g)) then
            allocate (res%theta_g(np))
            do k = 1, np
                res%theta_g(k) = 2*pi*(k - 1)/np
            end do
        end if
        do k = 1, np
            gn = norm2(gpsi(:, k))
            nh(:, k) = [sgn*gpsi(1, k)/gn, 0.0_dp, sgn*gpsi(2, k)/gn]
            bb(k) = norm2(b0(:, k))
            bhat(:, k) = b0(:, k)/bb(k)
            eb(:, k) = [bhat(2, k)*nh(3, k), &
                bhat(3, k)*nh(1, k) - bhat(1, k)*nh(3, k), -bhat(2, k)*nh(1, k)]
        end do
        res%bmag = bb
        res%nhat = nh
        res%ebin = eb
        if (res%dist_rational < tol) then
            status = 2
            return
        end if
        dnr = matmul(d, nh(1, :))
        dnz = matmul(d, nh(3, :))

        ! collocation matrix of B0.grad on the e^{i n phi} sector
        allocate (amat(np, np), rhs(np), rhs2(np), ipiv(np))
        do j = 1, np
            do k = 1, np
                amat(k, j) = cmplx(bth(k)*d(k, j), 0.0_dp, dp)
            end do
        end do
        do k = 1, np
            amat(k, k) = amat(k, k) + cmplx(0.0_dp, n*b0(2, k)/x(1, k), dp)
        end do
        ! normal part: B0.grad (xi.grad psi) = dB.grad psi
        do k = 1, np
            rhs(k) = dbe(1, k)*gpsi(1, k) + dbe(3, k)*gpsi(2, k)
        end do
        block
            complex(dp), allocatable :: a1(:, :)
            a1 = amat
            call zgesv(np, 1, a1, np, ipiv, rhs, np, info)
        end block
        if (info /= 0) then
            status = 3
            return
        end if
        allocate (res%xi_psi(np), res%xi_n(np), res%xi_b(np))
        res%xi_psi = rhs
        do k = 1, np
            res%xi_n(k) = sgn*rhs(k)/norm2(gpsi(:, k))
        end do
        ! binormal part: B0.grad xi_b - a xi_b = ebin.dB - xi_n c
        do k = 1, np
            v = cmplx(eb(:, k), 0.0_dp, dp)
            w = adv(k, v)
            aa(k) = real(dot_product(eb(:, k), real(w, dp)), dp)
            bn(:, k) = [bth(k)*dnr(k), b0(2, k)*nh(1, k)/x(1, k), bth(k)*dnz(k)]
            v = cmplx(nh(:, k), 0.0_dp, dp)
            w = adv(k, v)
            cc(k) = dot_product(eb(:, k), bn(:, k)) - dot_product(eb(:, k), &
                real(w, dp))
            rhs2(k) = dbe(1, k)*eb(1, k) + dbe(2, k)*eb(2, k) + dbe(3, k)*eb(3, k) &
                - res%xi_n(k)*cc(k)
            amat(k, k) = amat(k, k) - aa(k)
        end do
        call zgesv(np, 1, amat, np, ipiv, rhs2, np, info)
        if (info /= 0) then
            status = 3
            return
        end if
        res%xi_b = rhs2
        ! fields
        allocate (res%xi(3, np), res%dbl(3, np), res%dbpar_e(np), &
            res%dbmag_l(np), res%kappa_xi(np), res%kappa(3, np), res%gradb(2, np))
        do k = 1, np
            res%xi(:, k) = res%xi_n(k)*nh(:, k) + res%xi_b(k)*eb(:, k)
            res%dbl(:, k) = dbe(:, k) + adv(k, res%xi(:, k))
            do j = 1, 2
                res%gradb(j, k) = dot_product(bhat(:, k), gb(:, j, k))
            end do
            res%dbpar_e(k) = dot_product(bhat(:, k), dbe(:, k))
            res%dbmag_l(k) = res%dbpar_e(k) + res%xi(1, k)*res%gradb(1, k) &
                + res%xi(3, k)*res%gradb(2, k)
            vb = adv(k, cmplx(b0(:, k), 0.0_dp, dp))
            bv = real(vb, dp)
            pr = dot_product(bhat(:, k), bv)
            kk = (bv - bhat(:, k)*pr)/bb(k)**2
            res%kappa(:, k) = kk
            t = kk
            res%kappa_xi(k) = t(1)*res%xi(1, k) + t(2)*res%xi(2, k) + &
                t(3)*res%xi(3, k)
        end do
        status = 0

    contains

        function adv(k, v) result(o)
            !! (v.grad) B0 in cylindrical components at point k.
            integer, intent(in) :: k
            complex(dp), intent(in) :: v(3)
            complex(dp) :: o(3)
            o(1) = v(1)*gb(1, 1, k) + v(3)*gb(1, 2, k) - v(2)*b0(2, k)/x(1, k)
            o(2) = v(1)*gb(2, 1, k) + v(3)*gb(2, 2, k) + v(2)*b0(1, k)/x(1, k)
            o(3) = v(1)*gb(3, 1, k) + v(3)*gb(3, 2, k)
        end function adv

    end subroutine perturbation_solve

    ! ----------------------------------------------------------------- spectra

    subroutine perturbation_spectrum(s, angle, n, f, mmin, mmax, fm)
        !! Harmonics fm(m) (coefficient of e^{i n zeta - i m theta_s}) of the
        !! field f e^{i n phi} sampled at the surface points of s.
        type(surface_t), intent(in) :: s
        integer, intent(in) :: angle, n, mmin, mmax
        complex(dp), intent(in) :: f(:)
        complex(dp), intent(out) :: fm(mmin:mmax)
        complex(dp), allocatable :: c(:)
        integer :: j, m

        allocate (c(s%n))
        do j = 1, s%n
            c(j) = f(j)*exp(cmplx(0.0_dp, n*(s%nu(j) - s%q*s%theta_s(j, angle)), dp))
        end do
        do m = mmin, mmax
            fm(m) = (0.0_dp, 0.0_dp)
            do j = 1, s%n
                fm(m) = fm(m) + c(j)*s%dtheta_s(j, angle)* &
                    exp(cmplx(0.0_dp, m*s%theta_s(j, angle), dp))
            end do
            fm(m) = fm(m)/s%n
        end do
    end subroutine perturbation_spectrum

    subroutine perturbation_synthesis(s, angle, n, fm, mmin, mmax, f)
        !! Inverse of perturbation_spectrum on the window mmin..mmax.
        type(surface_t), intent(in) :: s
        integer, intent(in) :: angle, n, mmin, mmax
        complex(dp), intent(in) :: fm(mmin:mmax)
        complex(dp), intent(out) :: f(:)
        integer :: j, m

        do j = 1, s%n
            f(j) = (0.0_dp, 0.0_dp)
            do m = mmin, mmax
                f(j) = f(j) + fm(m)*exp(cmplx(0.0_dp, -m*s%theta_s(j, angle), dp))
            end do
            f(j) = f(j)*exp(cmplx(0.0_dp, -n*(s%nu(j) - s%q*s%theta_s(j, angle)), dp))
        end do
    end subroutine perturbation_synthesis

    subroutine perturbation_chart_map(s, n, from, to, fm_from, mmin, mmax, fm_to)
        !! Harmonics in chart `to` of the field whose harmonics in chart `from`
        !! are fm_from (exact up to the truncation of the window).
        type(surface_t), intent(in) :: s
        integer, intent(in) :: n, from, to, mmin, mmax
        complex(dp), intent(in) :: fm_from(mmin:mmax)
        complex(dp), intent(out) :: fm_to(mmin:mmax)
        complex(dp), allocatable :: f(:)

        allocate (f(s%n))
        call perturbation_synthesis(s, from, n, fm_from, mmin, mmax, f)
        call perturbation_spectrum(s, to, n, f, mmin, mmax, fm_to)
    end subroutine perturbation_chart_map

    subroutine perturbation_add_spectra(s, res, mmin, mmax)
        !! Fill res%spec for the PEST and Hamada charts of surface s.
        type(surface_t), intent(in) :: s
        type(perturbation_t), intent(inout) :: res
        integer, intent(in) :: mmin, mmax
        complex(dp), allocatable :: f(:, :)
        integer :: iq, ic, k

        allocate (f(res%np, nq_spec))
        f(:, 1) = res%xi_n
        f(:, 2) = res%xi_b
        do k = 1, res%np
            f(k, 3) = sum(res%nhat(:, k)*res%dbe(:, k))
        end do
        f(:, 4) = res%dbpar_e
        f(:, 5) = res%dbmag_l
        f(:, 6) = res%kappa_xi
        res%mmin = mmin
        res%mmax = mmax
        if (allocated(res%spec)) deallocate (res%spec)
        allocate (res%spec(mmin:mmax, nq_spec, 2))
        do ic = 1, 2
            do iq = 1, nq_spec
                call perturbation_spectrum(s, ic, res%n, f(:, iq), mmin, mmax, &
                    res%spec(:, iq, ic))
            end do
        end do
    end subroutine perturbation_add_spectra

    ! ------------------------------------------------------------------ export

    subroutine perturbation_export(path, res, source, cocos_declared, s_pol, ip_sign, &
            status)
        !! One documented ASCII file `path` (format kin6d-l7-v1) with all
        !! surfaces of res(:), and the sidecar path//'.sha256' in sha256sum
        !! format.  status 0 on success.
        character(len=*), intent(in) :: path, source
        type(perturbation_t), intent(in) :: res(:)
        integer, intent(in) :: cocos_declared
        real(dp), intent(in) :: s_pol, ip_sign
        integer, intent(out) :: status
        character(len=64) :: hex
        integer :: u, i, k, m, iq, ic, ios, slash

        status = 1
        open (newunit=u, file=path, status='replace', action='write', iostat=ios)
        if (ios /= 0) return
        write (u, '(a)') '# KIN6D L7 perturbation export'
        write (u, '(a)') 'format = kin6d-l7-v1'
        write (u, '(a)') 'source = '//trim(source)
        write (u, '(a,i0)') 'n = ', res(1)%n
        write (u, '(a)') 'n_sign = physical field X_phys = '// &
            'X exp(+i n phi) + c.c.; X stored is the complex coefficient of '// &
            'exp(+i n phi); n > 0 is exp(+i n phi) with phi the right-handed '// &
            '(R,phi,Z) toroidal angle'
        write (u, '(a,i0)') 'cocos_declared = ', cocos_declared
        write (u, '(a,f6.1,a,f6.1)') 'sgn_dpsi_outward = ', res(1)%sgn_dpsi, &
            '  s_pol = ', s_pol
        write (u, '(a,f6.1)') 'ip_sign = ', ip_sign
        write (u, '(a)') 'units = SI (T, m); xi in m, kappa in 1/m'
        write (u, '(a)') 'chart_points = cylindrical (R, phi, Z) components; '// &
            'points on flux surfaces psin = const, uniform in the geometric '// &
            'poloidal angle theta_g about the magnetic axis'
        write (u, '(a)') 'chart_spectra = harmonics X_m of X exp(i n phi) = '// &
            'exp(i n zeta) sum_m X_m exp(-i m theta_s), zeta = phi - (nu - '// &
            'q theta_s), theta_s = PEST (chart 1) or Hamada (chart 2); '// &
            'Boozer not provided'
        write (u, '(a)') 'displacement = ideal xi = xi_n nhat + xi_b ebin, '// &
            'xi_par = 0, nhat = sgn_dpsi grad psi/|grad psi|, ebin = b x nhat'
        write (u, '(a)') 'dB_L = dB_E + (xi.grad) B0 (vector), dBmag_L = '// &
            'b.dB_E + xi.grad|B|, kappa = (b.grad) b'
        write (u, '(a,i0)') 'surfaces = ', size(res)
        do i = 1, size(res)
            write (u, '(a,i0)') '[surface] ', i
            write (u, '(a,es25.16e3)') 'psin = ', res(i)%psin
            write (u, '(a,es25.16e3)') 'q_eff = ', res(i)%q_eff
            write (u, '(a,es25.16e3)') 'dist_rational = ', res(i)%dist_rational
            write (u, '(a,i0)') 'np = ', res(i)%np
            write (u, '(a)') '# columns: k theta_g R Z B0_R B0_phi B0_Z |B| '// &
                '(Re Im) xi_n xi_b dBE_R dBE_phi dBE_Z dBL_R dBL_phi dBL_Z '// &
                'dB_par_E dBmag_L kappa_xi; then real kappa_R kappa_phi kappa_Z'
            do k = 1, res(i)%np
                write (u, '(i6,38es25.16e3)') k, res(i)%theta_g(k), res(i)%x(:, k), &
                    res(i)%b0(:, k), res(i)%bmag(k), &
                    cx(res(i)%xi_n(k)), cx(res(i)%xi_b(k)), cx3(res(i)%dbe(:, k)), &
                    cx3(res(i)%dbl(:, k)), cx(res(i)%dbpar_e(k)), &
                    cx(res(i)%dbmag_l(k)), cx(res(i)%kappa_xi(k)), res(i)%kappa(:, k)
            end do
            if (allocated(res(i)%spec)) then
                write (u, '(a,2i6)') 'spectra mmin mmax = ', res(i)%mmin, res(i)%mmax
                write (u, '(a)') '# columns: chart quantity m Re Im'
                do ic = 1, 2
                    do iq = 1, nq_spec
                        do m = res(i)%mmin, res(i)%mmax
                            write (u, '(i2,1x,a10,i6,2es25.16e3)') ic, spec_names(iq), &
                                m, res(i)%spec(m, iq, ic)
                        end do
                    end do
                end do
            end if
        end do
        close (u)
        call sha256_file(path, hex, ios)
        if (ios /= 0) return
        slash = index(path, '/', back=.true.)
        open (newunit=u, file=path//'.sha256', status='replace', action='write', &
            iostat=ios)
        if (ios /= 0) return
        write (u, '(a)') trim(hex)//'  '//path(slash + 1:)
        close (u)
        status = 0

    contains

        function cx(z) result(a)
            complex(dp), intent(in) :: z
            real(dp) :: a(2)
            a = [real(z), aimag(z)]
        end function cx

        function cx3(z) result(a)
            complex(dp), intent(in) :: z(3)
            real(dp) :: a(6)
            a = [real(z(1)), aimag(z(1)), real(z(2)), aimag(z(2)), real(z(3)), &
                aimag(z(3))]
        end function cx3

    end subroutine perturbation_export

    ! ------------------------------------------------------------------ SHA-256

    subroutine sha256_file(path, hex, status)
        !! SHA-256 (FIPS 180-4) of the bytes of a file, lowercase hex.
        character(len=*), intent(in) :: path
        character(len=64), intent(out) :: hex
        integer, intent(out) :: status
        integer(int8), allocatable :: buf(:)
        integer(int64) :: sz
        integer :: u, ios

        hex = ''
        status = 1
        open (newunit=u, file=path, access='stream', form='unformatted', &
            status='old', action='read', iostat=ios)
        if (ios /= 0) return
        inquire (unit=u, size=sz)
        allocate (buf(sz))
        if (sz > 0) read (u, iostat=ios) buf
        close (u)
        if (ios /= 0) return
        call sha256_bytes(buf, hex)
        status = 0
    end subroutine sha256_file

    subroutine sha256_verify(path, ok)
        !! True when path//'.sha256' holds the SHA-256 of path.
        character(len=*), intent(in) :: path
        logical, intent(out) :: ok
        character(len=64) :: hex
        character(len=200) :: line
        integer :: u, ios

        ok = .false.
        call sha256_file(path, hex, ios)
        if (ios /= 0) return
        open (newunit=u, file=path//'.sha256', status='old', action='read', &
            iostat=ios)
        if (ios /= 0) return
        read (u, '(a)', iostat=ios) line
        close (u)
        if (ios /= 0) return
        ok = line(1:64) == hex
    end subroutine sha256_verify

    pure function rotr(x, k) result(y)
        integer(int64), intent(in) :: x
        integer, intent(in) :: k
        integer(int64) :: y
        y = iand(ior(shiftr(x, k), shiftl(x, 32 - k)), int(z'FFFFFFFF', int64))
    end function rotr

    subroutine sha256_bytes(buf, hex)
        integer(int8), intent(in) :: buf(:)
        character(len=64), intent(out) :: hex
        integer(int64), parameter :: mask = int(z'FFFFFFFF', int64)
        integer(int64) :: kc(64)
        data kc /z'428a2f98', z'71374491', z'b5c0fbcf', z'e9b5dba5', &
            z'3956c25b', z'59f111f1', z'923f82a4', z'ab1c5ed5', &
            z'd807aa98', z'12835b01', z'243185be', z'550c7dc3', &
            z'72be5d74', z'80deb1fe', z'9bdc06a7', z'c19bf174', &
            z'e49b69c1', z'efbe4786', z'0fc19dc6', z'240ca1cc', &
            z'2de92c6f', z'4a7484aa', z'5cb0a9dc', z'76f988da', &
            z'983e5152', z'a831c66d', z'b00327c8', z'bf597fc7', &
            z'c6e00bf3', z'd5a79147', z'06ca6351', z'14292967', &
            z'27b70a85', z'2e1b2138', z'4d2c6dfc', z'53380d13', &
            z'650a7354', z'766a0abb', z'81c2c92e', z'92722c85', &
            z'a2bfe8a1', z'a81a664b', z'c24b8b70', z'c76c51a3', &
            z'd192e819', z'd6990624', z'f40e3585', z'106aa070', &
            z'19a4c116', z'1e376c08', z'2748774c', z'34b0bcb5', &
            z'391c0cb3', z'4ed8aa4a', z'5b9cca4f', z'682e6ff3', &
            z'748f82ee', z'78a5636f', z'84c87814', z'8cc70208', &
            z'90befffa', z'a4506ceb', z'bef9a3f7', z'c67178f2'/
        integer(int64) :: h0(8), h(8), wv(64), a, b, c, d, e, f, gg, hh, t1, t2, s0, s1
        data h0 /z'6a09e667', z'bb67ae85', z'3c6ef372', z'a54ff53a', z'510e527f', &
            z'9b05688c', z'1f83d9ab', z'5be0cd19'/
        integer(int8), allocatable :: m(:)
        integer(int64) :: nbytes, nblk, bitlen
        integer :: i, j, blk
        character(len=16), parameter :: digits = '0123456789abcdef'

        nbytes = size(buf)
        nblk = (nbytes + 9 + 63)/64
        allocate (m(nblk*64))
        m = 0_int8
        if (nbytes > 0) m(1:nbytes) = buf
        m(nbytes + 1) = int(z'80', int8)
        bitlen = nbytes*8
        do i = 0, 7
            m(nblk*64 - i) = int(iand(shiftr(bitlen, 8*i), 255_int64), int8)
        end do
        h = h0
        do blk = 0, int(nblk) - 1
            do j = 1, 16
                wv(j) = 0
                do i = 0, 3
                    wv(j) = ior(shiftl(wv(j), 8), &
                        iand(int(m(blk*64 + (j - 1)*4 + i + 1), int64), 255_int64))
                end do
            end do
            do j = 17, 64
                s0 = ieor(ieor(rotr(wv(j - 15), 7), rotr(wv(j - 15), 18)), &
                    shiftr(wv(j - 15), 3))
                s1 = ieor(ieor(rotr(wv(j - 2), 17), rotr(wv(j - 2), 19)), &
                    shiftr(wv(j - 2), 10))
                wv(j) = iand(wv(j - 16) + s0 + wv(j - 7) + s1, mask)
            end do
            a = h(1); b = h(2); c = h(3); d = h(4)
            e = h(5); f = h(6); gg = h(7); hh = h(8)
            do j = 1, 64
                s1 = ieor(ieor(rotr(e, 6), rotr(e, 11)), rotr(e, 25))
                t1 = iand(hh + s1 + ieor(iand(e, f), iand(ieor(e, mask), gg)) &
                    + kc(j) + wv(j), mask)
                s0 = ieor(ieor(rotr(a, 2), rotr(a, 13)), rotr(a, 22))
                t2 = iand(s0 + ieor(ieor(iand(a, b), iand(a, c)), iand(b, c)), mask)
                hh = gg; gg = f; f = e
                e = iand(d + t1, mask)
                d = c; c = b; b = a
                a = iand(t1 + t2, mask)
            end do
            h = iand(h + [a, b, c, d, e, f, gg, hh], mask)
        end do
        do i = 1, 8
            do j = 0, 7
                hex((i - 1)*8 + j + 1:(i - 1)*8 + j + 1) = &
                    digits(int(iand(shiftr(h(i), 4*(7 - j)), 15_int64)) + 1: &
                    int(iand(shiftr(h(i), 4*(7 - j)), 15_int64)) + 1)
            end do
        end do
    end subroutine sha256_bytes

end module kin6d_mhd_perturbation_n
