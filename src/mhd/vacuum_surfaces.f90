module kin6d_vacuum_surfaces
    !! L3: flux-surface geometry from a poloidal-flux grid psi(R, Z) (G-EQDSK)
    !! and the projection of a per-n field onto straight-field-line harmonics.
    !!
    !! Conventions (all stated once here).
    !!  - Cylindrical (R, phi, Z) is right handed; phi is the coil phi.
    !!  - psi is the poloidal flux per radian (Wb/rad), psin = (psi - psi_axis) /
    !!    (psi_edge - psi_axis). The outward unit normal is n = sgn(dpsi) grad psi /
    !!    |grad psi|.
    !!  - Poloidal field B_pol = (s/R)(psi_Z, -psi_R) in (R, Z) with
    !!    s = sgn(Ip) sgn(psi_edge - psi_axis) (Ampere-consistent, any COCOS), and
    !!    B_phi = F/R.
    !!  - theta_g is the geometric angle about the magnetic axis, increasing
    !!    counter-clockwise in the (R, Z) plane. nu(theta_g) is the toroidal
    !!    advance of a field line, d nu/d theta_g = F |x'| / (R^2 B_t), with
    !!    B_t = B_pol . t the field component along the contour tangent t. The
    !!    signed safety factor is qs = (1/2pi) nu(2 pi); qs < 0 when the field
    !!    lines run clockwise in the poloidal plane while phi increases.
    !!  - Straight-field-line angle theta_s (theta_s = theta_pest or theta_hamada)
    !!    is a monotone function of theta_g; the straight toroidal angle is
    !!    zeta = phi - nu(theta_g) + qs theta_s, so field lines are zeta = const +
    !!    qs theta_s. PEST: d theta_s ~ d nu; Hamada: d theta_s ~ R |x'|/|grad psi|
    !!    (Jacobian constant).
    !!  - A per-n field B_n e^{i n phi} (physical field Re[2 B_n e^{i n phi}], n > 0)
    !!    has the harmonic expansion on a surface
    !!        2 B_n.n (theta_s, zeta) = sum_m b_mn e^{i (n zeta - m theta_s)}
    !!    and b_mn = (2/2pi) int B_n.n e^{i n (nu - qs theta_s)} e^{+i m theta_s}
    !!    d theta_s. It is a peak amplitude in T. The field line is resonant with
    !!    the harmonic m = n qs. The flux-density coefficient phi_mn uses
    !!    f = (B.n) R |x'| (T m^2 per rad^2, = J B.grad psi) in place of B.n and
    !!    no Jacobian factor d theta_s/d theta_g.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_real_identity, only: same_real_value
    implicit none
    private

    real(dp), parameter :: pi = 3.14159265358979323846_dp
    integer, parameter, public :: angle_pest = 1, angle_hamada = 2

    type, public :: psi_grid_t
        integer :: nr = 0, nz = 0, nf = 0
        real(dp) :: rleft = 0, dr = 0, zbot = 0, dz = 0
        real(dp) :: psi_axis = 0, psi_edge = 0, r_axis = 0, z_axis = 0
        real(dp) :: ip = 0, sgn_dpsi = 1, s_pol = 1
        real(dp), allocatable :: coef(:, :)    ! B-spline coefficients (0:nr+1, 0:nz+1)
        real(dp), allocatable :: fcoef(:)      ! F(psin) B-spline coefficients, (0:nf+1)
        real(dp), allocatable :: pcoef(:)      ! p(psin) (Pa), set by eqdsk_read
    end type

    type, public :: surface_t
        integer :: n = 0
        real(dp) :: psin = 0, q = 0
        real(dp), allocatable :: theta_g(:), r(:), z(:)
        real(dp), allocatable :: nr(:), nz(:)       ! outward unit normal
        real(dp), allocatable :: dl(:)              ! |x'| (m per rad of theta_g)
        real(dp), allocatable :: nu(:)              ! toroidal advance
        real(dp), allocatable :: theta_s(:, :)      ! (n, 2) pest, hamada
        real(dp), allocatable :: dtheta_s(:, :)     ! d theta_s / d theta_g
    end type

    public :: psi_grid_init, eqdsk_read, psi_eval, psin_of, f_of_psin, profiles_of_psin
    public :: surface_build, surface_harmonics, find_surface_q, flux_ray_point

contains

    ! ---------------------------------------------------------------- splines

    subroutine bspline_coefficients_1d(f, c)
        !! Cubic B-spline interpolation coefficients on a uniform grid with
        !! natural end conditions. c(0:n+1) such that
        !! f(i) = (c(i) + 4 c(i+1) + c(i+2))/6 for the point index i = 0..n-1
        !! (c(k) holds the coefficient of knot k-1).
        real(dp), intent(in) :: f(:)
        real(dp), intent(out) :: c(0:)
        real(dp), allocatable :: d(:), w(:)
        integer :: n, i
        real(dp) :: m
        n = size(f)
        allocate (d(n), w(n))
        ! unknowns u(i) = c(i+1), i = 1..n; interior rows u(i-1)+4u(i)+u(i+1) = 6 f(i)
        ! natural ends: u(1) = f(1), u(n) = f(n)
        d = 6*f
        w = 1                                 ! scratch for the super-diagonal
        c(1) = f(1)
        c(n) = f(n)
        ! Thomas algorithm for i = 2..n-1 with u(1), u(n) known
        d(2) = d(2) - f(1)
        d(n - 1) = d(n - 1) - f(n)
        w(2) = 4
        do i = 3, n - 1
            m = 1/w(i - 1)
            w(i) = 4 - m
            d(i) = d(i) - m*d(i - 1)
        end do
        c(n - 1) = d(n - 1)/w(n - 1)
        do i = n - 2, 2, -1
            c(i) = (d(i) - c(i + 1))/w(i)
        end do
        c(0) = 2*c(1) - c(2)
        c(n + 1) = 2*c(n) - c(n - 1)
    end subroutine

    pure subroutine basis(t, b, db, d2b)
        real(dp), intent(in) :: t
        real(dp), intent(out) :: b(4), db(4), d2b(4)
     b = [(1 - t)**3/6, (3*t**3 - 6*t**2 + 4)/6, (-3*t**3 + 3*t**2 + 3*t + 1)/6, t**3/6]
        db = [-(1 - t)**2/2, (3*t**2 - 4*t)/2, (-3*t**2 + 2*t + 1)/2, t**2/2]
        d2b = [1 - t, 3*t - 2, -3*t + 1, t]
    end subroutine

    subroutine psi_grid_init(g, psirz, rleft, rdim, zmid, zdim, &
                             r_axis_guess, z_axis_guess, psi_edge, ip, fpol, status)
        !! psirz(nr, nz) on the uniform grid R = rleft + (i-1) dR, Z = zbot + (j-1) dZ
        type(psi_grid_t), intent(out) :: g
        real(dp), intent(in) :: psirz(:, :), rleft, rdim, zmid, zdim
        real(dp), intent(in) :: r_axis_guess, z_axis_guess, psi_edge, ip
        real(dp), intent(in) :: fpol(:)   ! F on uniform psin = 0..1
        integer, intent(out) :: status
        real(dp), allocatable :: tmp(:, :)
        real(dp) :: r, z, psi, pr, pz, prr, prz, pzz, det, dr, dz
        integer :: i, j, it

        status = 1
        g%nr = size(psirz, 1); g%nz = size(psirz, 2); g%nf = size(fpol)
        if (g%nr < 8 .or. g%nz < 8 .or. g%nf < 8) return
        g%rleft = rleft; g%dr = rdim/(g%nr - 1)
        g%zbot = zmid - zdim/2; g%dz = zdim/(g%nz - 1)
        g%ip = ip; g%psi_edge = psi_edge
        allocate (g%coef(0:g%nr + 1, 0:g%nz + 1), tmp(0:g%nr + 1, g%nz), &
                  g%fcoef(0:g%nf + 1))
        do j = 1, g%nz
            call bspline_coefficients_1d(psirz(:, j), tmp(:, j))
        end do
        do i = 0, g%nr + 1
            call bspline_coefficients_1d(tmp(i, :), g%coef(i, :))
        end do
        call bspline_coefficients_1d(fpol, g%fcoef)
        ! magnetic axis: Newton on grad psi = 0
        r = r_axis_guess; z = z_axis_guess
        do it = 1, 50
            call psi_eval(g, r, z, psi, pr, pz, prr, prz, pzz)
            det = prr*pzz - prz*prz
            if (abs(det) <= 0.0_dp) return
            dr = -(pzz*pr - prz*pz)/det
            dz = -(-prz*pr + prr*pz)/det
            r = r + dr; z = z + dz
            if (abs(dr) + abs(dz) < 1.0e-14_dp) exit
        end do
        if (it > 50) return
        call psi_eval(g, r, z, psi, pr, pz)
        g%r_axis = r; g%z_axis = z; g%psi_axis = psi
        if (same_real_value(psi_edge, psi)) return
        g%sgn_dpsi = sign(1.0_dp, psi_edge - psi)
        g%s_pol = sign(1.0_dp, ip)*g%sgn_dpsi
        status = 0
    end subroutine

    subroutine psi_eval(g, r, z, psi, psi_r, psi_z, psi_rr, psi_rz, psi_zz)
        type(psi_grid_t), intent(in) :: g
        real(dp), intent(in) :: r, z
        real(dp), intent(out) :: psi, psi_r, psi_z
        real(dp), intent(out), optional :: psi_rr, psi_rz, psi_zz
        real(dp) :: u, v, tr, tz, br(4), dbr(4), d2br(4), bz(4), dbz(4), d2bz(4)
        real(dp) :: blk(4, 4)
        integer :: i0, j0
        u = (r - g%rleft)/g%dr + 1       ! data point i sits at u = i
        v = (z - g%zbot)/g%dz + 1
        i0 = max(1, min(g%nr - 1, floor(u)))
        j0 = max(1, min(g%nz - 1, floor(v)))
        tr = u - i0
        tz = v - j0
        call basis(tr, br, dbr, d2br)
        call basis(tz, bz, dbz, d2bz)
        blk = g%coef(i0 - 1:i0 + 2, j0 - 1:j0 + 2)
        psi = dot_product(br, matmul(blk, bz))
        psi_r = dot_product(dbr, matmul(blk, bz))/g%dr
        psi_z = dot_product(br, matmul(blk, dbz))/g%dz
        if (present(psi_rr)) then
            psi_rr = dot_product(d2br, matmul(blk, bz))/g%dr**2
            psi_rz = dot_product(dbr, matmul(blk, dbz))/(g%dr*g%dz)
            psi_zz = dot_product(br, matmul(blk, d2bz))/g%dz**2
        end if
    end subroutine

    pure function psin_of(g, psi) result(x)
        type(psi_grid_t), intent(in) :: g
        real(dp), intent(in) :: psi
        real(dp) :: x
        x = (psi - g%psi_axis)/(g%psi_edge - g%psi_axis)
    end function

    function f_of_psin(g, x) result(f)
        type(psi_grid_t), intent(in) :: g
        real(dp), intent(in) :: x
        real(dp) :: f, u, t, b(4), db(4), d2b(4)
        integer :: i0
        u = max(0.0_dp, min(1.0_dp, x))*(g%nf - 1) + 1
        i0 = max(1, min(g%nf - 1, floor(u)))
        t = u - i0
        call basis(t, b, db, d2b)
        f = dot_product(b, g%fcoef(i0 - 1:i0 + 2))
    end function

    subroutine profiles_of_psin(g, x, f, dfdx, p, dpdx)
        !! F and p with their derivatives in psin (cubic B-splines on the
        !! uniform psin grid of the G-EQDSK profiles; p = 0 without pcoef).
        type(psi_grid_t), intent(in) :: g
        real(dp), intent(in) :: x
        real(dp), intent(out) :: f, dfdx, p, dpdx
        real(dp) :: u, t, b(4), db(4), d2b(4)
        integer :: i0
        u = max(0.0_dp, min(1.0_dp, x))*(g%nf - 1) + 1
        i0 = max(1, min(g%nf - 1, floor(u)))
        t = u - i0
        call basis(t, b, db, d2b)
        f = dot_product(b, g%fcoef(i0 - 1:i0 + 2))
        dfdx = dot_product(db, g%fcoef(i0 - 1:i0 + 2))*(g%nf - 1)
        p = 0.0_dp
        dpdx = 0.0_dp
        if (allocated(g%pcoef)) then
            p = dot_product(b, g%pcoef(i0 - 1:i0 + 2))
            dpdx = dot_product(db, g%pcoef(i0 - 1:i0 + 2))*(g%nf - 1)
        end if
    end subroutine profiles_of_psin

    ! ----------------------------------------------------------------- EQDSK

    subroutine eqdsk_read(path, g, status)
        type(psi_grid_t), intent(out) :: g
        character(len=*), intent(in) :: path
        integer, intent(out) :: status
        integer :: unit, ios, nw, nh, k
        character(len=80) :: head
        real(dp) :: rdim, zdim, rcentr, rleft, zmid, rmaxis, zmaxis, simag, sibry
        real(dp) :: bcentr, current, xdum
        real(dp), allocatable :: fpol(:), work(:), psirz(:, :), pres(:)

        status = 1
        open (newunit=unit, file=path, status='old', action='read', iostat=ios)
        if (ios /= 0) return
        read (unit, '(a)', iostat=ios) head
        if (ios /= 0) return
        k = index(trim(head), ' ', back=.true.)
        read (head(k + 1:), *, iostat=ios) nh
        if (ios /= 0) return
        head = head(1:k - 1)
        k = index(trim(head), ' ', back=.true.)
        read (head(k + 1:), *, iostat=ios) nw
        if (ios /= 0) return
        read (unit, '(5e16.9)', iostat=ios) rdim, zdim, rcentr, rleft, zmid
        read (unit, '(5e16.9)', iostat=ios) rmaxis, zmaxis, simag, sibry, bcentr
        read (unit, '(5e16.9)', iostat=ios) current, simag, xdum, rmaxis, xdum
        read (unit, '(5e16.9)', iostat=ios) zmaxis, xdum, sibry, xdum, xdum
        if (ios /= 0) return
        allocate (fpol(nw), work(nw), psirz(nw, nh), pres(nw))
        read (unit, '(5e16.9)', iostat=ios) fpol
        read (unit, '(5e16.9)', iostat=ios) pres
        read (unit, '(5e16.9)', iostat=ios) work       ! ffprim
        read (unit, '(5e16.9)', iostat=ios) work       ! pprime
        read (unit, '(5e16.9)', iostat=ios) psirz
        close (unit)
        if (ios /= 0) return
        call psi_grid_init(g, psirz, rleft, rdim, zmid, zdim, rmaxis, zmaxis, &
                           sibry, current, fpol, status)
        if (status /= 0) return
        allocate (g%pcoef(0:nw + 1))
        call bspline_coefficients_1d(pres, g%pcoef)
    end subroutine

    ! -------------------------------------------------------------- surfaces

    subroutine periodic_antiderivative(gfun, big_g, mean)
        !! Spectral antiderivative on theta_j = 2 pi (j-1)/n of the periodic
        !! function gfun: G(theta) = mean*theta + periodic, G(0) = 0.
        real(dp), intent(in) :: gfun(:)
        real(dp), intent(out) :: big_g(:)
        real(dp), intent(out) :: mean
        complex(dp), allocatable :: ghat(:), tab(:)
        integer :: n, j, k, kk
        complex(dp) :: s
        n = size(gfun)
        allocate (ghat(-n/2:n/2), tab(0:n - 1))
        ghat = 0
        do k = 0, n - 1
            tab(k) = cmplx(cos(2*pi*k/n), sin(2*pi*k/n), dp)
        end do
        do kk = -n/2 + 1, n/2 - 1
            s = 0
            do j = 1, n
                s = s + gfun(j)*conjg(tab(modulo(kk*(j - 1), n)))
            end do
            ghat(kk) = s/n
        end do
        mean = real(ghat(0), dp)
        do j = 1, n
            s = 0
            do kk = -n/2 + 1, n/2 - 1
                if (kk == 0) cycle
                s = s + ghat(kk)*(tab(modulo(kk*(j - 1), n)) - 1)/cmplx(0.0_dp, kk, dp)
            end do
            big_g(j) = mean*2*pi*(j - 1)/n + real(s, dp)
        end do
    end subroutine

    subroutine flux_ray_point(g, psin, th, rho, ok)
        !! Distance rho from the magnetic axis along the ray of geometric angle
        !! th to the first crossing of the level psin (march, then safeguarded
        !! Newton).  ok = .false. if the ray leaves the grid first.
        type(psi_grid_t), intent(in) :: g
        real(dp), intent(in) :: psin, th
        real(dp), intent(out) :: rho
        logical, intent(out) :: ok
        real(dp) :: er(2), lo, hi, psi, pr, pz, target, f, rp, ds, hi_grid
        integer :: it

        ok = .false.
        target = g%psi_axis + psin*(g%psi_edge - g%psi_axis)
        er = [cos(th), sin(th)]
        ! longest ray that stays inside the grid (one cell margin)
        hi_grid = huge(1.0_dp)
        if (er(1) > 1.0e-12_dp) hi_grid = &
            min(hi_grid, (g%rleft + g%dr*(g%nr - 2) - g%r_axis)/er(1))
        if (er(1) < -1.0e-12_dp) hi_grid = min(hi_grid, (g%rleft + g%dr - g%r_axis)/er(1))
        if (er(2) > 1.0e-12_dp) hi_grid = &
            min(hi_grid, (g%zbot + g%dz*(g%nz - 2) - g%z_axis)/er(2))
        if (er(2) < -1.0e-12_dp) hi_grid = min(hi_grid, (g%zbot + g%dz - g%z_axis)/er(2))
        ! march outward to the first crossing of the target level
        lo = 0
        ds = 0.02_dp
        hi = -1
        rho = 0
        do while (rho < hi_grid)
            rho = min(rho + ds, hi_grid)
            call psi_eval(g, g%r_axis + rho*er(1), g%z_axis + rho*er(2), psi, pr, pz)
            if ((psi - target)*g%sgn_dpsi > 0) then
                hi = rho
                exit
            end if
            lo = rho
        end do
        if (hi < 0) return
        rho = 0.5_dp*(lo + hi)
        do it = 1, 200
            call psi_eval(g, g%r_axis + rho*er(1), g%z_axis + rho*er(2), psi, pr, pz)
            f = psi - target
            if (abs(f) < 1.0e-12_dp*abs(g%psi_edge - g%psi_axis)) then
                ok = .true.
                return
            end if
            if (f*g%sgn_dpsi > 0) then
                hi = rho
            else
                lo = rho
            end if
            rp = pr*er(1) + pz*er(2)
            if (.not. (abs(rp) <= 0.0_dp)) then
                rho = rho - f/rp
            else
                rho = 0.5_dp*(lo + hi)
            end if
            if (rho <= lo .or. rho >= hi) rho = 0.5_dp*(lo + hi)
        end do
    end subroutine flux_ray_point

    subroutine surface_build(g, psin, n, s, status)
        !! Flux surface psin sampled at n equally spaced geometric angles.
        type(psi_grid_t), intent(in) :: g
        real(dp), intent(in) :: psin
        integer, intent(in) :: n
        type(surface_t), intent(out) :: s
        integer, intent(out) :: status
        real(dp), allocatable :: dnu(:), wj(:)
        real(dp) :: th, er(2), et(2), rho, psi, pr, pz, fval
        real(dp) :: rp, t(2), bt, numean
        integer :: j
        logical :: ok

        status = 1
        s%n = n; s%psin = psin
        allocate (s%theta_g(n), s%r(n), s%z(n), s%nr(n), s%nz(n), s%dl(n), s%nu(n), &
                  s%theta_s(n, 2), s%dtheta_s(n, 2), dnu(n), wj(n))
        fval = f_of_psin(g, psin)
        do j = 1, n
            th = 2*pi*(j - 1)/n
            er = [cos(th), sin(th)]
            et = [-sin(th), cos(th)]
            call flux_ray_point(g, psin, th, rho, ok)
            if (.not. ok) return
            s%theta_g(j) = th
            s%r(j) = g%r_axis + rho*er(1)
            s%z(j) = g%z_axis + rho*er(2)
            call psi_eval(g, s%r(j), s%z(j), psi, pr, pz)
            ! contour tangent dx/dtheta = rho' e_r + rho e_theta with
            ! rho' = -rho (grad psi.e_t)/(grad psi.e_r)
            rp = -rho*(pr*et(1) + pz*et(2))/(pr*er(1) + pz*er(2))
            t = rp*er + rho*et
            s%dl(j) = norm2(t)
            s%nr(j) = g%sgn_dpsi*pr/norm2([pr, pz])
            s%nz(j) = g%sgn_dpsi*pz/norm2([pr, pz])
            ! B_pol.t = (s/R)(psi_Z t_R - psi_R t_Z); dnu/dtheta = F |x'| / (R^2 B_t)
            bt = g%s_pol/s%r(j)*(pz*t(1) - pr*t(2))/s%dl(j)
            dnu(j) = fval*s%dl(j)/(s%r(j)**2*bt)
            wj(j) = s%r(j)*s%dl(j)/norm2([pr, pz])
        end do
        call periodic_antiderivative(dnu, s%nu, numean)
        s%q = numean       ! (1/2pi) * int dnu = mean of dnu
        s%dtheta_s(:, angle_pest) = dnu/numean
        s%theta_s(:, angle_pest) = s%nu/numean
        call hamada_angle(wj, s)
        status = 0
    end subroutine

    subroutine hamada_angle(wj, s)
        real(dp), intent(in) :: wj(:)
        type(surface_t), intent(inout) :: s
        real(dp), allocatable :: big_g(:)
        real(dp) :: wmean
        allocate (big_g(size(wj)))
        call periodic_antiderivative(wj, big_g, wmean)
        s%theta_s(:, angle_hamada) = big_g/wmean
        s%dtheta_s(:, angle_hamada) = wj/wmean
    end subroutine

    function find_surface_q(g, qtarget, lo, hi, n, tol, psin, status) result(found)
        !! Bisection on psin in [lo, hi] for qs = qtarget (signed q monotone).
        type(psi_grid_t), intent(in) :: g
        real(dp), intent(in) :: qtarget, lo, hi, tol
        integer, intent(in) :: n
        real(dp), intent(out) :: psin
        integer, intent(out) :: status
        logical :: found
        type(surface_t) :: s
        real(dp) :: a, b, fa, fb, fm, m
        integer :: it, st
        found = .false.; status = 1
        a = lo; b = hi
        call surface_build(g, a, n, s, st); if (st /= 0) return
        fa = s%q - qtarget
        call surface_build(g, b, n, s, st); if (st /= 0) return
        fb = s%q - qtarget
        if (fa*fb > 0) return
        do it = 1, 80
            m = 0.5_dp*(a + b)
            call surface_build(g, m, n, s, st); if (st /= 0) return
            fm = s%q - qtarget
            if (fa*fm <= 0) then
                b = m; fb = fm
            else
                a = m; fa = fm
            end if
            if (b - a < tol) exit
        end do
        psin = 0.5_dp*(a + b)
        found = .true.; status = 0
    end function

    subroutine surface_harmonics(s, n, bn_rz, angle, mmin, mmax, bmn, phimn)
        !! Peak-amplitude harmonics b_mn (T) and flux-density coefficients phi_mn
        !! (T m^2) of the per-n field with meridional components bn_rz(1:2, j)
        !! = (B_R, B_Z) of B_n at the surface points. m runs over mmin..mmax.
        type(surface_t), intent(in) :: s
        integer, intent(in) :: n, angle, mmin, mmax
        complex(dp), intent(in) :: bn_rz(:, :)
        complex(dp), intent(out) :: bmn(mmin:mmax), phimn(mmin:mmax)
        complex(dp), allocatable :: c(:), ph(:)
        integer :: j, m, nn, kk
        real(dp) :: ts
        nn = s%n
        allocate (c(nn), ph(nn))
        do j = 1, nn
            ts = s%theta_s(j, angle)
            c(j) = (bn_rz(1, j)*s%nr(j) + bn_rz(2, j)*s%nz(j))* &
                   exp(cmplx(0.0_dp, n*(s%nu(j) - s%q*ts), dp))
        end do
        do m = mmin, mmax
            bmn(m) = 0; phimn(m) = 0
            do j = 1, nn
                ph(j) = exp(cmplx(0.0_dp, m*s%theta_s(j, angle), dp))
                bmn(m) = bmn(m) + c(j)*s%dtheta_s(j, angle)*ph(j)
                phimn(m) = phimn(m) + c(j)*s%r(j)*s%dl(j)*ph(j)
            end do
        end do
        kk = 2
        bmn = kk*bmn/nn
        phimn = kk*phimn/nn
    end subroutine

end module
