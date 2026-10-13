module kin6d_equilibrium_geometry
    ! Kinetic geometry view of the D7 evaluator (docs/equilibrium-core.md, D7).
    ! For an axisymmetric nested state it reads one flux surface psi = psi_0
    ! through the batched evaluator only (B, p, psi, grad psi); it carries no
    ! flux-coordinate assumption beyond: the surface is star-shaped about a
    ! given interior axis point and psi is monotone along rays from it.
    !
    ! Conventions (cylindrical basis (R, Z, phi) as in the evaluator rows):
    !   theta_g   geometric angle about the axis, counter-clockwise, 0 outboard;
    !   t         unit tangent of the surface in the (R, Z) plane along +theta_g;
    !   q         signed (1/2 pi) oint B_phi dl /(R B_pol.t) = (1/2 pi) oint
    !             B_phi dl/|grad psi| (positive when a field line advances in
    !             +phi along +theta_g);
    !   nu        nu(theta_g) = int_0^theta_g B_phi/(R B_pol.t) dl, nu(2 pi) = 2 pi q;
    !   alpha     field-line label alpha = phi - nu(psi, theta_g), so B.grad alpha
    !             = 0 and alpha = phi on the outboard midplane;
    !   PEST      theta* = nu/q (straight field lines, b.grad theta* = B_phi/(q R B));
    !   Boozer    theta_B = int h/ oint h with h = R B^2 dl/|grad psi| (Jacobian
    !             ~ 1/B^2), zeta_B = phi + q (theta_B - theta*); the covariant
    !             B_phi R = F and the axisymmetric poloidal current I are not
    !             stored (F is R B_phi);
    !   V'        dV/dpsi = sign(psi_0 - psi_axis) 2 pi oint R dl/|grad psi|;
    !   <f>       oint f R dl/|grad psi| / oint R dl/|grad psi|.
    ! All contour integrals are uniform sums in theta_g (spectrally accurate for
    ! smooth periodic integrands). nu is held as a Fourier series, so any angle
    ! (unwrapped for field lines) is evaluated without interpolation; d nu/d psi
    ! (the shear term of alpha) and q' follow from five-point differences of the
    ! series of five surfaces at fixed theta_g. grad|B| is a five-point stencil
    ! of the evaluator in (R, Z). The drift factors use grad-B drift, i.e. low-beta
    ! curvature = grad_perp B/B; they are v_d = (m (v_par^2 + v_perp^2/2)/(Z e))
    ! b x grad B/B^2 per unit factor. Curvature at finite beta is not included.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use kin6d_equilibrium_state, only: equilibrium_field_t, equilibrium_state_t, &
        labelled_field_t, EQ_BR, EQ_BZ, EQ_BPHI, EQ_PSI, EQ_PSI_R, EQ_PSI_Z, &
        EQ_VALUES
    implicit none
    private
    public :: kinetic_surface_t, kinetic_samples_t, build_kinetic_surface
    public :: build_kinetic_surface_from_state, sample_kinetic_geometry
    public :: surface_average, surface_average_function
    public :: KG_ANGLE_GEOMETRIC, KG_ANGLE_PEST, KG_ANGLE_BOOZER

    integer, parameter :: KG_ANGLE_GEOMETRIC = 1, KG_ANGLE_PEST = 2, &
        KG_ANGLE_BOOZER = 3
    real(dp), parameter :: two_pi = 2*acos(-1.0_dp)

    ! f(theta) = mean + sum_k a_k cos(k theta) + b_k sin(k theta).
    type :: series_t
        real(dp) :: mean = 0
        real(dp), allocatable :: a(:), b(:)
    end type

    ! Geometry at a list of angles; arrays are (components, points).
    type :: kinetic_samples_t
        integer :: kind = KG_ANGLE_GEOMETRIC
        real(dp), allocatable :: theta(:)      ! requested angle (of `kind`)
        real(dp), allocatable :: theta_g(:), theta_pest(:), theta_boozer(:)
        real(dp), allocatable :: nu(:)         ! unwrapped; alpha = phi - nu
        real(dp), allocatable :: boozer_shift(:) ! zeta_B - phi
        real(dp), allocatable :: r(:), z(:)
        real(dp), allocatable :: b(:, :)       ! (R, Z, phi) components
        real(dp), allocatable :: bmag(:)
        real(dp), allocatable :: grad_b(:, :)  ! grad|B| (R, Z, phi), phi = 0
        real(dp), allocatable :: grad_psi(:, :) ! (R, Z)
        real(dp), allocatable :: grad_psi_mag(:)
        real(dp), allocatable :: grad_alpha(:, :) ! (R, Z, phi)
        real(dp), allocatable :: dalpha_dpsi(:) ! at fixed theta_g
        real(dp), allocatable :: b_grad_theta(:) ! b.grad of the angle of `kind`
        real(dp), allocatable :: bxgradpsi_gradtheta(:) ! (b x grad psi).grad theta
        real(dp), allocatable :: drift_psi(:)   ! (b x grad|B|).grad psi/|B|^2
        real(dp), allocatable :: drift_alpha(:) ! (b x grad|B|).grad alpha/|B|^2
    end type

    type :: kinetic_surface_t
        integer :: n_theta = 0
        real(dp) :: axis(2) = 0, psi = 0, psi_axis = 0
        real(dp) :: q = 0, dq_dpsi = 0, dv_dpsi = 0, fd_step = 0
        real(dp), allocatable :: weight(:)   ! flux-surface average weights, sum 1
        type(kinetic_samples_t) :: node      ! uniform theta_g nodes
        type(series_t) :: nu_series, nu_psi_series, boozer_series
    end type

    type :: contour_t
        real(dp), allocatable :: r(:), x(:, :), val(:, :), tangent(:, :)
        real(dp), allocatable :: dl(:), g(:), h(:), jac(:)
    end type
contains
    subroutine build_kinetic_surface_from_state(eq, rho, surf, status, n_theta)
        !! Surface rho = sqrt((psi - psi_b)/(psi_a - psi_b)), the readback levels.
        type(equilibrium_state_t), intent(in) :: eq
        real(dp), intent(in) :: rho
        type(kinetic_surface_t), intent(out) :: surf
        integer, intent(out) :: status
        integer, intent(in), optional :: n_theta
        real(dp) :: axis(2), psi_axis
        status = 1
        if (.not. allocated(eq%field)) return
        select type (field => eq%field)
        class is (labelled_field_t)
            call field%label_axis(axis, psi_axis, status)
            if (status /= 0) then
                status = 1
                return
            end if
            call build_kinetic_surface(field, axis, &
                eq%psi_bound + (psi_axis - eq%psi_bound)*(1 - rho**2), surf, &
                status, n_theta)
        end select
    end subroutine

    subroutine build_kinetic_surface(field, axis, psi, surf, status, n_theta, &
                                     dpsi)
        !! Surface psi = psi of `field`, star-shaped about `axis` (O-point of psi).
        !! dpsi: step of the five-point psi differences (default 4e-3 of the
        !! flux depth).
        class(equilibrium_field_t), intent(in) :: field
        real(dp), intent(in) :: axis(2), psi
        type(kinetic_surface_t), intent(out) :: surf
        integer, intent(out) :: status
        integer, intent(in), optional :: n_theta
        real(dp), intent(in), optional :: dpsi
        type(contour_t) :: c(-2:2)
        type(kinetic_samples_t) :: node
        real(dp), allocatable :: theta(:), gk(:, :)
        real(dp) :: values(EQ_VALUES, 1), depth, delta, outward, point(2, 1)
        integer :: n, k, j, st(1)
        status = 1
        n = 128
        if (present(n_theta)) n = n_theta
        if (n < 16) return
        point(:, 1) = axis
        call field%sample(point, values, st)
        if (st(1) /= 0) return
        surf%axis = axis; surf%psi = psi; surf%psi_axis = values(EQ_PSI, 1)
        depth = psi - surf%psi_axis
        if (abs(depth) <= tiny(1.0_dp) .or. .not. ieee_is_finite(depth)) return
        outward = sign(1.0_dp, depth)
        delta = 4.0e-3_dp*depth
        if (present(dpsi)) delta = dpsi
        if (abs(delta) <= 0 .or. abs(2*delta) >= abs(depth)) return
        allocate (theta(n))
        theta = [(two_pi*(j - 1)/n, j=1, n)]
        allocate (gk(n, -2:2))
        do k = -2, 2
            call trace_contour(field, axis, psi + k*delta, surf%psi_axis, theta, &
                c(k), status)
            if (status /= 0) return
            gk(:, k) = c(k)%g
        end do
        surf%n_theta = n
        surf%fd_step = 2.0e-3_dp*sum(c(0)%r)/n
        surf%nu_series = fourier(c(0)%g)
        surf%q = surf%nu_series%mean
        if (abs(surf%q) <= tiny(1.0_dp)) then
            status = 1
            return
        end if
        surf%nu_psi_series = fourier((-gk(:, 2) + 8*gk(:, 1) - 8*gk(:, -1) &
            + gk(:, -2))/(12*delta))
        surf%dq_dpsi = surf%nu_psi_series%mean
        surf%boozer_series = fourier(c(0)%h)
        surf%dv_dpsi = outward*two_pi*sum(c(0)%jac)*two_pi/n
        allocate (surf%weight(n))
        surf%weight = c(0)%jac/sum(c(0)%jac)
        call sample_kinetic_geometry(field, surf, KG_ANGLE_GEOMETRIC, theta, &
            node, status)
        if (status == 0) surf%node = node
    end subroutine

    subroutine sample_kinetic_geometry(field, surf, kind, theta, s, status)
        !! Geometry at the angles `theta` of the given kind (geometric, PEST or
        !! Boozer). Angles need not lie in [0, 2 pi): nu and the angles are
        !! unwrapped (field lines), the geometry is periodic.
        class(equilibrium_field_t), intent(in) :: field
        type(kinetic_surface_t), intent(in) :: surf
        integer, intent(in) :: kind
        real(dp), intent(in) :: theta(:)
        type(kinetic_samples_t), intent(out) :: s
        integer, intent(out) :: status
        type(contour_t) :: c
        real(dp), allocatable :: tg(:), pts(:, :), vals(:, :), bm(:, :)
        real(dp) :: h, cs, sn, rr, dth, bvec(3), gb(3), ga(3), gpsi(2), nu_psi
        real(dp) :: bx(3), gt(2), cr(3)
        integer, allocatable :: sts(:)
        integer :: m, i, p
        status = 1
        m = size(theta)
        if (kind < KG_ANGLE_GEOMETRIC .or. kind > KG_ANGLE_BOOZER) return
        if (surf%n_theta == 0 .or. m == 0) return
        allocate (tg(m))
        do i = 1, m
            call to_geometric(surf, kind, theta(i), tg(i), status)
            if (status /= 0) return
        end do
        status = 1
        call trace_contour(field, surf%axis, surf%psi, surf%psi_axis, tg, c, status)
        if (status /= 0) return
        status = 1
        h = surf%fd_step
        ! Five-point stencil of |B| in R and Z: points 1..4 are R +-h, +-2h.
        allocate (pts(2, 8*m), vals(EQ_VALUES, 8*m), sts(8*m), bm(8, m))
        do i = 1, m
            do p = 1, 4
                pts(:, 8*(i - 1) + p) = c%x(:, i) + [stencil(p)*h, 0.0_dp]
                pts(:, 8*(i - 1) + 4 + p) = c%x(:, i) + [0.0_dp, stencil(p)*h]
            end do
        end do
        call field%sample(pts, vals, sts)
        if (any(sts /= 0)) return
        do i = 1, m
            do p = 1, 8
                bm(p, i) = norm2(vals([EQ_BR, EQ_BZ, EQ_BPHI], 8*(i - 1) + p))
            end do
        end do
        s%kind = kind
        allocate (s%theta(m), s%theta_g(m), s%theta_pest(m), s%theta_boozer(m), &
            s%nu(m), s%boozer_shift(m), s%r(m), s%z(m), s%b(3, m), s%bmag(m), &
            s%grad_b(3, m), s%grad_psi(2, m), s%grad_psi_mag(m), &
            s%grad_alpha(3, m), s%dalpha_dpsi(m), s%b_grad_theta(m), &
            s%bxgradpsi_gradtheta(m), s%drift_psi(m), s%drift_alpha(m))
        do i = 1, m
            cs = cos(tg(i)); sn = sin(tg(i))
            rr = c%r(i)
            s%theta(i) = theta(i); s%theta_g(i) = tg(i)
            s%r(i) = c%x(1, i); s%z(i) = c%x(2, i)
            bvec = c%val([EQ_BR, EQ_BZ, EQ_BPHI], i)
            s%b(:, i) = bvec; s%bmag(i) = norm2(bvec)
            gpsi = c%val([EQ_PSI_R, EQ_PSI_Z], i)
            s%grad_psi(:, i) = gpsi; s%grad_psi_mag(i) = norm2(gpsi)
            gb = 0
            gb(1) = dot_product(bm(1:4, i), [-1.0_dp, 8.0_dp, -8.0_dp, 1.0_dp]) &
                /(12*h)
            gb(2) = dot_product(bm(5:8, i), [-1.0_dp, 8.0_dp, -8.0_dp, 1.0_dp]) &
                /(12*h)
            s%grad_b(:, i) = gb
            s%nu(i) = integral(surf%nu_series, tg(i))
            s%theta_pest(i) = s%nu(i)/surf%q
            s%theta_boozer(i) = integral(surf%boozer_series, tg(i)) &
                /surf%boozer_series%mean
            s%boozer_shift(i) = surf%q*(s%theta_boozer(i) - s%theta_pest(i))
            ! grad theta_g = (-(Z - Z_a), R - R_a)/r^2 in (R, Z).
            gt = [-sn, cs]/rr
            nu_psi = integral(surf%nu_psi_series, tg(i))
            s%dalpha_dpsi(i) = -nu_psi
            dth = value(surf%nu_series, tg(i))
            ga(1:2) = -(nu_psi*gpsi + dth*gt)
            ga(3) = 1/s%r(i)
            s%grad_alpha(:, i) = ga
            select case (kind)
            case (KG_ANGLE_PEST)
                dth = dth/surf%q
            case (KG_ANGLE_BOOZER)
                dth = value(surf%boozer_series, tg(i))/surf%boozer_series%mean
            case default
                dth = 1
            end select
            s%b_grad_theta(i) = dth*dot_product(bvec(1:2), gt)/s%bmag(i)
            ! (b x grad psi).grad theta: only the in-plane theta_g part survives.
            cr = cross(bvec/s%bmag(i), [gpsi, 0.0_dp])
            s%bxgradpsi_gradtheta(i) = dth*dot_product(cr(1:2), gt)
            bx = cross(bvec/s%bmag(i), gb)
            s%drift_psi(i) = dot_product(bx(1:2), gpsi)/s%bmag(i)**2
            s%drift_alpha(i) = dot_product(bx, ga)/s%bmag(i)**2
        end do
        status = 0
    end subroutine

    pure function stencil(p) result(x)
        integer, intent(in) :: p
        real(dp), parameter :: offsets(4) = [2.0_dp, 1.0_dp, -1.0_dp, -2.0_dp]
        real(dp) :: x
        x = offsets(p)
    end function

    pure function cross(a, b) result(c)
        ! a x b for (R, Z, phi) ordered components; (R, phi, Z) is right-handed.
        real(dp), intent(in) :: a(3), b(3)
        real(dp) :: c(3), u(3), v(3), w(3)
        u = [a(1), a(3), a(2)]; v = [b(1), b(3), b(2)]
        w = [u(2)*v(3) - u(3)*v(2), u(3)*v(1) - u(1)*v(3), u(1)*v(2) - u(2)*v(1)]
        c = [w(1), w(3), w(2)]
    end function

    real(dp) function surface_average(surf, values) result(avg)
        !! <f> from f at the nodes surf%node (uniform theta_g).
        type(kinetic_surface_t), intent(in) :: surf
        real(dp), intent(in) :: values(:)
        avg = ieee_nan()
        if (size(values) /= surf%n_theta) return
        avg = dot_product(surf%weight, values)
    end function

    real(dp) function surface_average_function(surf, f) result(avg)
        !! <f(R, Z)> for a function of position.
        type(kinetic_surface_t), intent(in) :: surf
        interface
            real(dp) function f(r, z)
                import :: dp
                real(dp), intent(in) :: r, z
            end function
        end interface
        real(dp), allocatable :: v(:)
        integer :: i
        allocate (v(surf%n_theta))
        do i = 1, surf%n_theta
            v(i) = f(surf%node%r(i), surf%node%z(i))
        end do
        avg = surface_average(surf, v)
    end function

    real(dp) function ieee_nan() result(x)
        use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
        x = ieee_value(0.0_dp, ieee_quiet_nan)
    end function

    ! ---- Fourier series -------------------------------------------------

    function fourier(f) result(s)
        real(dp), intent(in) :: f(:)
        type(series_t) :: s
        real(dp), allocatable :: ck(:), sk(:)
        integer :: n, nk, k, j
        n = size(f); nk = (n - 1)/2
        allocate (s%a(nk), s%b(nk), ck(0:n - 1), sk(0:n - 1))
        ck = [(cos(two_pi*j/n), j=0, n - 1)]
        sk = [(sin(two_pi*j/n), j=0, n - 1)]
        s%mean = sum(f)/n
        do k = 1, nk
            s%a(k) = 0; s%b(k) = 0
            do j = 0, n - 1
                s%a(k) = s%a(k) + f(j + 1)*ck(mod(j*k, n))
                s%b(k) = s%b(k) + f(j + 1)*sk(mod(j*k, n))
            end do
            s%a(k) = 2*s%a(k)/n; s%b(k) = 2*s%b(k)/n
        end do
    end function

    pure real(dp) function value(s, theta)
        type(series_t), intent(in) :: s
        real(dp), intent(in) :: theta
        integer :: k
        value = s%mean
        do k = 1, size(s%a)
            value = value + s%a(k)*cos(k*theta) + s%b(k)*sin(k*theta)
        end do
    end function

    pure real(dp) function integral(s, theta)
        ! int_0^theta of the series; unwrapped.
        type(series_t), intent(in) :: s
        real(dp), intent(in) :: theta
        integer :: k
        integral = s%mean*theta
        do k = 1, size(s%a)
            integral = integral + (s%a(k)*sin(k*theta) &
                + s%b(k)*(1 - cos(k*theta)))/k
        end do
    end function

    subroutine to_geometric(surf, kind, theta, tg, status)
        ! Invert theta = Theta(theta_g) (monotone, Theta(x + 2 pi) = Theta + 2 pi,
        ! Theta(2 pi k) = 2 pi k): the root lies in the cell holding theta.
        type(kinetic_surface_t), intent(in) :: surf
        integer, intent(in) :: kind
        real(dp), intent(in) :: theta
        real(dp), intent(out) :: tg
        integer, intent(out) :: status
        real(dp) :: lo, hi, f, d, step
        integer :: it
        status = 1
        tg = theta
        if (.not. ieee_is_finite(theta)) return
        if (kind == KG_ANGLE_GEOMETRIC) then
            status = 0
            return
        end if
        lo = two_pi*floor(theta/two_pi); hi = lo + two_pi
        tg = min(hi, max(lo, theta))
        do it = 1, 100
            if (kind == KG_ANGLE_PEST) then
                f = integral(surf%nu_series, tg)/surf%q - theta
                d = value(surf%nu_series, tg)/surf%q
            else
                f = integral(surf%boozer_series, tg)/surf%boozer_series%mean - theta
                d = value(surf%boozer_series, tg)/surf%boozer_series%mean
            end if
            if (f > 0) then
                hi = tg
            else
                lo = tg
            end if
            if (d > 0) then
                step = -f/d
            else
                step = hi - lo
            end if
            if (tg + step <= lo .or. tg + step >= hi) then
                step = (lo + hi)/2 - tg
            end if
            tg = tg + step
            if (abs(step) <= 4*epsilon(1.0_dp)*max(1.0_dp, abs(tg))) then
                status = 0
                return
            end if
        end do
    end subroutine

    ! ---- Contour and rays -------------------------------------------------

    subroutine trace_contour(field, axis, psi, psi_axis, theta, c, status)
        ! Surface points at the rays theta and the series integrands
        !   g = B_phi |dx/dtheta|/(R B_pol.t),   h = R B^2 |dx/dtheta|/|grad psi|,
        !   jac = R |dx/dtheta|/|grad psi|
        ! with dr/dtheta = -psi_theta/psi_r from the evaluator gradient.
        class(equilibrium_field_t), intent(in) :: field
        real(dp), intent(in) :: axis(2), psi, psi_axis, theta(:)
        type(contour_t), intent(out) :: c
        integer, intent(out) :: status
        real(dp) :: cs, sn, psi_r, psi_t, rt, dx(2), rr, bpt, b2, gp
        integer, allocatable :: sts(:)
        integer :: m, i
        m = size(theta)
        allocate (c%r(m), c%x(2, m), c%val(EQ_VALUES, m), c%tangent(2, m), &
            c%dl(m), c%g(m), c%h(m), c%jac(m), sts(m))
        call trace_rays(field, axis, psi, psi_axis, theta, c%r, status)
        if (status /= 0) return
        status = 1
        do i = 1, m
            c%x(:, i) = axis + c%r(i)*[cos(theta(i)), sin(theta(i))]
        end do
        call field%sample(c%x, c%val, sts)
        if (any(sts /= 0)) return
        do i = 1, m
            cs = cos(theta(i)); sn = sin(theta(i)); rr = c%r(i)
            psi_r = cs*c%val(EQ_PSI_R, i) + sn*c%val(EQ_PSI_Z, i)
            psi_t = rr*(-sn*c%val(EQ_PSI_R, i) + cs*c%val(EQ_PSI_Z, i))
            if (abs(psi_r) <= tiny(1.0_dp)) return
            rt = -psi_t/psi_r
            dx = [rt*cs - rr*sn, rt*sn + rr*cs]
            c%dl(i) = norm2(dx)
            c%tangent(:, i) = dx/c%dl(i)
            bpt = c%val(EQ_BR, i)*c%tangent(1, i) + c%val(EQ_BZ, i)*c%tangent(2, i)
            b2 = sum(c%val([EQ_BR, EQ_BZ, EQ_BPHI], i)**2)
            gp = norm2(c%val([EQ_PSI_R, EQ_PSI_Z], i))
            if (abs(bpt) <= tiny(1.0_dp) .or. gp <= tiny(1.0_dp) .or. &
                c%x(1, i) <= 0) return
            c%g(i) = c%val(EQ_BPHI, i)*c%dl(i)/(c%x(1, i)*bpt)
            c%jac(i) = c%x(1, i)*c%dl(i)/gp
            c%h(i) = b2*c%jac(i)
        end do
        if (.not. all(ieee_is_finite(c%g)) .or. .not. all(ieee_is_finite(c%h))) return
        status = 0
    end subroutine

    subroutine trace_rays(field, axis, psi_t, psi_axis, theta, r, status)
        ! First crossing of psi = psi_t along each ray from the axis: domain edge
        ! by doubling and bisection, a table to bracket, then bracketed Newton.
        class(equilibrium_field_t), intent(in) :: field
        real(dp), intent(in) :: axis(2), psi_t, psi_axis, theta(:)
        real(dp), intent(out) :: r(:)
        integer, intent(out) :: status
        integer, parameter :: table = 256, max_iter = 100
        real(dp), allocatable :: pts(:, :), vals(:, :), inside(:), outside(:)
        real(dp), allocatable :: lo(:), hi(:), cs(:), sn(:)
        integer, allocatable :: sts(:)
        logical, allocatable :: done(:), capped(:), active(:)
        real(dp) :: s0, cap, tol, d, step, f, fprime, tm
        integer :: m, j, i, it, k
        status = 1
        m = size(theta)
        s0 = sign(1.0_dp, psi_axis - psi_t)
        cap = 8*abs(axis(1))
        tol = 8*epsilon(1.0_dp)*max(abs(psi_t), abs(psi_axis - psi_t))
        allocate (cs(m), sn(m), inside(m), outside(m), done(m), capped(m), &
            active(m), lo(m), hi(m), pts(2, max(m, m*(table + 1))), &
            vals(EQ_VALUES, max(m, m*(table + 1))), sts(max(m, m*(table + 1))))
        cs = cos(theta); sn = sin(theta)
        inside = 0; outside = 1.0e-3_dp*abs(axis(1)); done = .false.
        capped = .false.
        do it = 1, 64
            if (all(done)) exit
            do j = 1, m
                pts(:, j) = axis + outside(j)*[cs(j), sn(j)]
            end do
            call field%sample(pts(:, :m), vals(:, :m), sts(:m))
            do j = 1, m
                if (done(j)) cycle
                if (sts(j) == 0 .and. pts(1, j) > 0) then
                    inside(j) = outside(j)
                    if (outside(j) >= cap) then
                        done(j) = .true.; capped(j) = .true.
                    else
                        outside(j) = min(2*outside(j), cap)
                    end if
                else
                    done(j) = .true.
                end if
            end do
        end do
        if (.not. all(done) .or. any(inside <= 0)) return
        do it = 1, 40
            do j = 1, m
                pts(:, j) = axis + (inside(j) + outside(j))/2*[cs(j), sn(j)]
            end do
            call field%sample(pts(:, :m), vals(:, :m), sts(:m))
            do j = 1, m
                if (capped(j)) cycle
                if (sts(j) == 0 .and. pts(1, j) > 0) then
                    inside(j) = (inside(j) + outside(j))/2
                else
                    outside(j) = (inside(j) + outside(j))/2
                end if
            end do
        end do
        do j = 1, m
            do i = 0, table
                pts(:, (j - 1)*(table + 1) + i + 1) = axis &
                    + inside(j)*i/table*[cs(j), sn(j)]
            end do
        end do
        call field%sample(pts(:, :m*(table + 1)), vals(:, :m*(table + 1)), &
            sts(:m*(table + 1)))
        if (any(sts(:m*(table + 1)) /= 0)) return
        active = .false.
        do j = 1, m
            do i = 1, table
                k = (j - 1)*(table + 1) + i + 1
                if (s0*(vals(EQ_PSI, k) - psi_t) <= 0) then
                    lo(j) = inside(j)*(i - 1)/table; hi(j) = inside(j)*i/table
                    active(j) = .true.
                    exit
                end if
            end do
            if (.not. active(j)) return
            r(j) = (lo(j) + hi(j))/2
        end do
        do it = 1, max_iter
            if (.not. any(active)) exit
            do j = 1, m
                pts(:, j) = axis + r(j)*[cs(j), sn(j)]
            end do
            call field%sample(pts(:, :m), vals(:, :m), sts(:m))
            do j = 1, m
                if (.not. active(j)) cycle
                if (sts(j) /= 0) return
                f = vals(EQ_PSI, j) - psi_t
                if (abs(f) <= tol .or. hi(j) - lo(j) <= &
                    4*epsilon(1.0_dp)*(abs(axis(1)) + inside(j))) then
                    active(j) = .false.
                    cycle
                end if
                if (s0*f > 0) then
                    lo(j) = r(j)
                else
                    hi(j) = r(j)
                end if
                fprime = cs(j)*vals(EQ_PSI_R, j) + sn(j)*vals(EQ_PSI_Z, j)
                step = 0; tm = (lo(j) + hi(j))/2
                if (abs(fprime) > tiny(1.0_dp)) step = -f/fprime
                d = r(j) + step
                if (abs(fprime) > tiny(1.0_dp) .and. d > lo(j) .and. d < hi(j)) then
                    r(j) = d
                else
                    r(j) = tm
                end if
            end do
        end do
        if (any(active)) return
        status = 0
    end subroutine
end module kin6d_equilibrium_geometry
