!> Kinetic geometry on one surface of a Boozer-coordinate `.bc` file, the
!! input format of NEO-2, for comparisons on the same field.
!!
!! The file is read by libneo's `bc_file`, in the six-column
!! (stellarator-symmetric) or ten-column layout. On the Boozer angles
!! (theta, phi_B), with n the file's mode number without the field period N_p,
!!   R = sum rmnc cos(a) + rmns sin(a), Z = sum zmns sin(a) + zmnc cos(a),
!!   B = sum bmnc cos(a) + bmns sin(a),
!!   phi = phi_B - (2 pi/N_p) sum (vmns sin(a) + vmnc cos(a)),
!!   a = m theta - n N_p phi_B,
!! as in NEO-2 (`neo_sub.f90`). The radial factor s^p uses
!! p = min(|m|, 12 + modulo(|m|, 2))/2, matching NEO-2's bounded regularity
!! normalization. The cap limits roundoff amplification near the axis while
!! preserving even/odd parity. Normalized coefficients and iota use not-a-knot
!! cubic splines; surfaces outside the tabulated range are refused.
!! Lengths in m, |B| in T, |grad s| in 1/m.
module milca_boozer_file
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use bc_file, only: bc_data_t, read_bc_file
    use fortnum_cubic_spline, only: cubic_spline_t, cubic_spline_init, &
        cubic_spline_eval
    use milca_boozer_file_leaves, only: boozer_axis_factor, boozer_metric
    implicit none
    private
    public :: boozer_file_t, boozer_surface_t, boozer_point_t
    public :: load_boozer_file, boozer_surface, boozer_point, boozer_average
    public :: boozer_jacobian_ratio
    public :: boozer_field_t, boozer_field, boozer_field_eval, boozer_iota, &
        boozer_iota_integral, boozer_harmonics
    integer, parameter, public :: BOOZER_OK = 0, BOOZER_OUTSIDE = 1

    real(dp), parameter :: twopi = 8*atan(1.0_dp)
    real(dp), parameter :: mu0 = 2*twopi*1e-7_dp

    !> iota, G, I and the B harmonics of one toroidal mode number, cut from the
    !> file's spline: the columns keep their nodes and slopes, so the values
    !> are those of boozer_surface, without evaluating the other columns.
    type :: boozer_field_t
        integer :: n_tor = 0
        real(dp) :: psi_prime = 0
        ! File indices of the modes with n nper = n_tor, in file order.
        integer, allocatable :: mode(:), m(:)
        ! Radial power of each mode as a power of sqrt(s), and the largest.
        integer, allocatable :: half_power(:)
        integer :: max_half_power = 0
        ! Hermite data of the columns node by node: values (column, node) and
        ! slopes, so that one evaluation reads contiguous memory.
        real(dp), allocatable :: y(:, :), slope(:, :)
        ! int_{s_1}^{s_k} iota ds of the iota spline at its nodes.
        real(dp), allocatable :: iota_int(:)
        ! Columns: iota, G, I, bmnc and bmns of the modes divided by s^power.
        type(cubic_spline_t) :: spline
        type(cubic_spline_t) :: iota
    end type

    type :: boozer_file_t
        type(bc_data_t) :: d
        ! Columns: iota, then rmnc, zmns, vmns, bmnc, rmns, zmnc, vmnc, bmns
        ! divided by s^radial_power(m), then the poloidal and toroidal currents.
        type(cubic_spline_t) :: spline
        ! The axisymmetric part (n = 0), for guiding-centre orbits.
        type(boozer_field_t) :: axi
    end type

    type :: boozer_surface_t
        real(dp) :: s = 0, iota = 0, diota_ds = 0
        ! Covariant Boozer components G, I [T m], their s derivatives and
        ! dpsi_tor/ds [T m^2].
        real(dp) :: g_cov = 0, i_cov = 0, g_cov_s = 0, i_cov_s = 0, psi_prime = 0
        integer :: nper = 1
        integer, allocatable :: m(:), n(:)
        ! Symmetric parts (rmnc, zmns, vmns, bmnc) and their s derivatives,
        ! then the asymmetric partners (rmns, zmnc, vmnc, bmns), zero for a
        ! stellarator-symmetric file.
        real(dp), allocatable :: rmn(:), zmn(:), vmn(:), bmn(:)
        real(dp), allocatable :: rmn_s(:), zmn_s(:), vmn_s(:), bmn_s(:)
        real(dp), allocatable :: rma(:), zma(:), vma(:), bma(:)
        real(dp), allocatable :: rma_s(:), zma_s(:), vma_s(:), bma_s(:)
    end type

    type :: boozer_point_t
        real(dp) :: r = 0, z = 0, phi = 0, bmod = 0
        real(dp) :: jacobian = 0 ! X_s . (X_theta x X_phiB), m^3
        real(dp) :: grad_s = 0 ! |grad s|, 1/m
    end type

contains

    pure real(dp) function radial_power(m) result(power)
        integer, intent(in) :: m
        ! Limit amplification of roundoff near the axis, preserving parity,
        ! as in NEO-2's radial spline normalization.
        integer, parameter :: max_m = 12
        power = 0.5_dp*min(abs(m), max_m + modulo(abs(m), 2))
    end function radial_power

    subroutine load_boozer_file(filename, bf, status)
        character(*), intent(in) :: filename
        type(boozer_file_t), intent(out) :: bf
        integer, intent(out) :: status
        real(dp), allocatable :: y(:, :)
        real(dp) :: f, f_s
        integer :: k, i, nm

        call read_bc_file(filename, bf%d)
        nm = bf%d%nmode
        allocate (y(bf%d%nsurf, 3 + 8*nm))
        do k = 1, bf%d%nsurf
            y(k, 1) = bf%d%iota(k)
            do i = 1, nm
                call boozer_axis_factor(bf%d%s(k), radial_power(bf%d%m(i)), f, f_s)
                y(k, 1 + i) = bf%d%rmnc(k, i)/f
                y(k, 1 + nm + i) = bf%d%zmns(k, i)/f
                y(k, 1 + 2*nm + i) = bf%d%vmns(k, i)/f
                y(k, 1 + 3*nm + i) = bf%d%bmnc(k, i)/f
                y(k, 1 + 4*nm + i) = bf%d%rmns(k, i)/f
                y(k, 1 + 5*nm + i) = bf%d%zmnc(k, i)/f
                y(k, 1 + 6*nm + i) = bf%d%vmnc(k, i)/f
                y(k, 1 + 7*nm + i) = bf%d%bmns(k, i)/f
            end do
            y(k, 2 + 8*nm) = bf%d%Jpol_over_nper(k)*bf%d%nper
            y(k, 3 + 8*nm) = bf%d%Itor(k)
        end do
        call cubic_spline_init(bf%d%s, y, bf%spline, status)
        if (status /= 0) then
            status = BOOZER_OUTSIDE
            return
        end if
        call boozer_field(bf, 0, bf%axi)
    end subroutine load_boozer_file

    !> The field of the modes with n nper = n_tor of a loaded file.
    subroutine boozer_field(bf, n_tor, fld)
        type(boozer_file_t), intent(in) :: bf
        integer, intent(in) :: n_tor
        type(boozer_field_t), intent(out) :: fld
        integer, allocatable :: col(:)
        real(dp) :: pw
        integer :: nm, nk
        nm = bf%d%nmode
        fld%n_tor = n_tor
        fld%psi_prime = bf%d%flux/twopi
        fld%mode = pack([(nk, nk=1, nm)], bf%d%n*bf%d%nper == n_tor)
        fld%m = bf%d%m(fld%mode)
        allocate (fld%half_power(size(fld%m)))
        do nk = 1, size(fld%m)
            fld%half_power(nk) = nint(2*radial_power(fld%m(nk)))
        end do
        fld%max_half_power = max(0, maxval(fld%half_power))
        col = [1, 2 + 8*nm, 3 + 8*nm, 1 + 3*nm + fld%mode, 1 + 7*nm + fld%mode]
        fld%spline%x = bf%spline%x
        allocate (fld%spline%y(size(bf%spline%x), size(col)), &
            fld%spline%slope(size(bf%spline%x), size(col)))
        do nk = 1, size(col)
            fld%spline%y(:, nk) = bf%spline%y(:, col(nk))
            fld%spline%slope(:, nk) = bf%spline%slope(:, col(nk))
        end do
        fld%y = transpose(fld%spline%y)
        fld%slope = transpose(fld%spline%slope)
        fld%iota%x = bf%spline%x
        fld%iota%y = bf%spline%y(:, 1:1)
        fld%iota%slope = bf%spline%slope(:, 1:1)
        allocate (fld%iota_int(size(fld%iota%x)))
        fld%iota_int(1) = 0
        do nk = 2, size(fld%iota%x)
            pw = fld%iota%x(nk) - fld%iota%x(nk - 1)
            fld%iota_int(nk) = fld%iota_int(nk - 1) + pw*(fld%iota%y(nk - 1, 1) &
                + fld%iota%y(nk, 1))/2 + pw**2*(fld%iota%slope(nk - 1, 1) &
                - fld%iota%slope(nk, 1))/12
        end do
    end subroutine boozer_field

    !> iota, G, I [T m], dG/ds, dI/ds and the harmonics bmnc, bmns [T] of fld
    !> with their s derivatives at s, as boozer_surface gives them.
    subroutine boozer_field_eval(fld, s, iota, g_cov, i_cov, g_cov_s, i_cov_s, bmn, bma, &
            bmn_s, bma_s, status)
        type(boozer_field_t), intent(in) :: fld
        real(dp), intent(in) :: s
        real(dp), intent(out) :: iota, g_cov, i_cov, g_cov_s, i_cov_s
        real(dp), intent(out) :: bmn(:), bma(:), bmn_s(:), bma_s(:)
        integer, intent(out) :: status
        real(dp) :: g(3 + 2*size(fld%m)), g_s(3 + 2*size(fld%m)), f, f_s
        real(dp) :: rs(0:fld%max_half_power), h, t, h00, h10, h01, h11, d00, d10, d01, d11
        integer :: k, nm, n, lo, hi, mid
        status = BOOZER_OUTSIDE
        if (.not. allocated(fld%y)) return
        n = size(fld%spline%x)
        if (s < fld%spline%x(1) .or. s > fld%spline%x(n)) return
        nm = size(fld%m)
        ! The Hermite form of cubic_spline_eval on the node-major data.
        lo = 1
        hi = n - 1
        do while (hi > lo)
            mid = (lo + hi + 1)/2
            if (s >= fld%spline%x(mid)) then
                lo = mid
            else
                hi = mid - 1
            end if
        end do
        h = fld%spline%x(lo + 1) - fld%spline%x(lo)
        t = (s - fld%spline%x(lo))/h
        h00 = (1 + 2*t)*(1 - t)**2
        h10 = h*t*(1 - t)**2
        h01 = t*t*(3 - 2*t)
        h11 = h*t*t*(t - 1)
        d00 = 6*t*(t - 1)/h
        d10 = (1 - t)*(1 - 3*t)
        d01 = -d00
        d11 = t*(3*t - 2)
        g = h00*fld%y(:, lo) + h10*fld%slope(:, lo) + h01*fld%y(:, lo + 1) &
            + h11*fld%slope(:, lo + 1)
        g_s = d00*fld%y(:, lo) + d10*fld%slope(:, lo) + d01*fld%y(:, lo + 1) &
            + d11*fld%slope(:, lo + 1)
        iota = g(1)
        g_cov = mu0*g(2)/twopi
        i_cov = mu0*g(3)/twopi
        g_cov_s = mu0*g_s(2)/twopi
        i_cov_s = mu0*g_s(3)/twopi
        ! s^(j/2) of radial_power by products of sqrt(s), d/ds = (j/2) s^(j/2)/s.
        rs(0) = 1
        if (fld%max_half_power > 0) rs(1) = sqrt(s)
        do k = 2, fld%max_half_power
            rs(k) = rs(k - 1)*rs(1)
        end do
        do k = 1, nm
            f = rs(fld%half_power(k))
            f_s = 0.5_dp*fld%half_power(k)*f/s
            bmn(k) = f*g(3 + k)
            bma(k) = f*g(3 + nm + k)
            bmn_s(k) = f_s*g(3 + k) + f*g_s(3 + k)
            bma_s(k) = f_s*g(3 + nm + k) + f*g_s(3 + nm + k)
        end do
        status = BOOZER_OK
    end subroutine boozer_field_eval

    !> int_{s0}^s iota ds of the iota spline, exactly: the closed-form
    !> integrals of its cubic pieces, node sums and partial pieces taken
    !> apart so that short ranges keep their digits; zero outside the file.
    pure real(dp) function boozer_iota_integral(fld, s0, s) result(v)
        type(boozer_field_t), intent(in) :: fld
        real(dp), intent(in) :: s0, s
        integer :: n, k0, k1
        v = 0
        if (.not. allocated(fld%iota_int)) return
        n = size(fld%iota%x)
        if (min(s0, s) < fld%iota%x(1) .or. max(s0, s) > fld%iota%x(n)) return
        k0 = interval(s0)
        k1 = interval(s)
        v = (fld%iota_int(k1) - fld%iota_int(k0)) + (piece(k1, s) - piece(k0, s0))
    contains
        pure integer function interval(x) result(lo)
            real(dp), intent(in) :: x
            integer :: hi, mid
            lo = 1
            hi = n - 1
            do while (hi > lo)
                mid = (lo + hi + 1)/2
                if (x >= fld%iota%x(mid)) then
                    lo = mid
                else
                    hi = mid - 1
                end if
            end do
        end function interval

        !> int from node lo to x of the Hermite piece (basis h00, h10, h01, h11).
        pure real(dp) function piece(lo, x) result(q)
            integer, intent(in) :: lo
            real(dp), intent(in) :: x
            real(dp) :: h, t
            h = fld%iota%x(lo + 1) - fld%iota%x(lo)
            t = (x - fld%iota%x(lo))/h
            q = h*(fld%iota%y(lo, 1)*(t**4/2 - t**3 + t) &
                + h*fld%iota%slope(lo, 1)*(t**4/4 - 2*t**3/3 + t**2/2) &
                + fld%iota%y(lo + 1, 1)*(t**3 - t**4/2) &
                + h*fld%iota%slope(lo + 1, 1)*(t**4/4 - t**3/3))
        end function piece
    end function boozer_iota_integral

    !> cos(m theta) and sin(m theta) by angle addition from one sincos (rounding
    !> grows like |m| epsilon).
    pure subroutine boozer_harmonics(m, theta, cm, sm)
        integer, intent(in) :: m(:)
        real(dp), intent(in) :: theta
        real(dp), intent(out) :: cm(:), sm(:)
        complex(dp) :: z(0:max(0, maxval(abs(m)))), e
        integer :: k
        e = cmplx(cos(theta), sin(theta), dp)
        z(0) = 1
        do k = 1, ubound(z, 1)
            z(k) = z(k - 1)*e
        end do
        cm = real(z(abs(m)), dp)
        sm = sign(1, m)*aimag(z(abs(m)))
    end subroutine boozer_harmonics

    !> iota at s, zero outside the file as from boozer_surface.
    real(dp) function boozer_iota(fld, s) result(iota)
        type(boozer_field_t), intent(in) :: fld
        real(dp), intent(in) :: s
        real(dp) :: v(1), v_s(1)
        iota = 0
        if (.not. allocated(fld%iota%x)) return
        if (s < fld%iota%x(1) .or. s > fld%iota%x(size(fld%iota%x))) return
        call cubic_spline_eval(fld%iota, s, v, v_s)
        iota = v(1)
    end function boozer_iota

    subroutine boozer_surface(bf, s, surf, status)
        type(boozer_file_t), intent(in) :: bf
        real(dp), intent(in) :: s
        type(boozer_surface_t), intent(out) :: surf
        integer, intent(out) :: status
        real(dp), allocatable :: g(:), g_s(:)
        real(dp) :: f, f_s
        integer :: nm, i

        status = BOOZER_OUTSIDE
        if (.not. allocated(bf%spline%x)) return
        if (s < bf%d%s(1) .or. s > bf%d%s(bf%d%nsurf)) return
        nm = bf%d%nmode
        allocate (g(3 + 8*nm), g_s(3 + 8*nm))
        call cubic_spline_eval(bf%spline, s, g, g_s)
        surf%s = s
        surf%nper = bf%d%nper
        surf%m = bf%d%m
        surf%n = bf%d%n
        surf%iota = g(1)
        surf%diota_ds = g_s(1)
        surf%g_cov = mu0*g(2 + 8*nm)/twopi
        surf%i_cov = mu0*g(3 + 8*nm)/twopi
        surf%g_cov_s = mu0*g_s(2 + 8*nm)/twopi
        surf%i_cov_s = mu0*g_s(3 + 8*nm)/twopi
        surf%psi_prime = bf%d%flux/twopi
        allocate (surf%rmn(nm), surf%zmn(nm), surf%vmn(nm), surf%bmn(nm), &
            surf%rmn_s(nm), surf%zmn_s(nm), surf%vmn_s(nm), surf%bmn_s(nm), &
            surf%rma(nm), surf%zma(nm), surf%vma(nm), surf%bma(nm), &
            surf%rma_s(nm), surf%zma_s(nm), surf%vma_s(nm), surf%bma_s(nm))
        do i = 1, nm
            call boozer_axis_factor(s, radial_power(bf%d%m(i)), f, f_s)
            surf%rmn(i) = f*g(1 + i)
            surf%zmn(i) = f*g(1 + nm + i)
            surf%vmn(i) = f*g(1 + 2*nm + i)
            surf%bmn(i) = f*g(1 + 3*nm + i)
            surf%rmn_s(i) = f_s*g(1 + i) + f*g_s(1 + i)
            surf%zmn_s(i) = f_s*g(1 + nm + i) + f*g_s(1 + nm + i)
            surf%vmn_s(i) = f_s*g(1 + 2*nm + i) + f*g_s(1 + 2*nm + i)
            surf%bmn_s(i) = f_s*g(1 + 3*nm + i) + f*g_s(1 + 3*nm + i)
            surf%rma(i) = f*g(1 + 4*nm + i)
            surf%zma(i) = f*g(1 + 5*nm + i)
            surf%vma(i) = f*g(1 + 6*nm + i)
            surf%bma(i) = f*g(1 + 7*nm + i)
            surf%rma_s(i) = f_s*g(1 + 4*nm + i) + f*g_s(1 + 4*nm + i)
            surf%zma_s(i) = f_s*g(1 + 5*nm + i) + f*g_s(1 + 5*nm + i)
            surf%vma_s(i) = f_s*g(1 + 6*nm + i) + f*g_s(1 + 6*nm + i)
            surf%bma_s(i) = f_s*g(1 + 7*nm + i) + f*g_s(1 + 7*nm + i)
        end do
        status = BOOZER_OK
    end subroutine boozer_surface

    pure function boozer_point(surf, theta, phi_b) result(p)
        type(boozer_surface_t), intent(in) :: surf
        real(dp), intent(in) :: theta, phi_b
        type(boozer_point_t) :: p
        real(dp) :: a, c, sn, mm, nn, v, r_s, r_t, r_p, z_s, z_t, z_p
        real(dp) :: v_s, v_t, v_p, shift
        integer :: i

        p%r = 0; p%z = 0; p%bmod = 0; v = 0
        r_s = 0; r_t = 0; r_p = 0; z_s = 0; z_t = 0; z_p = 0
        v_s = 0; v_t = 0; v_p = 0
        do i = 1, size(surf%m)
            mm = surf%m(i)
            nn = surf%n(i)*surf%nper
            a = mm*theta - nn*phi_b
            c = cos(a); sn = sin(a)
            p%r = p%r + surf%rmn(i)*c + surf%rma(i)*sn
            p%z = p%z + surf%zmn(i)*sn + surf%zma(i)*c
            p%bmod = p%bmod + surf%bmn(i)*c + surf%bma(i)*sn
            v = v + surf%vmn(i)*sn + surf%vma(i)*c
            r_s = r_s + surf%rmn_s(i)*c + surf%rma_s(i)*sn
            r_t = r_t + mm*(surf%rma(i)*c - surf%rmn(i)*sn)
            r_p = r_p + nn*(surf%rmn(i)*sn - surf%rma(i)*c)
            z_s = z_s + surf%zmn_s(i)*sn + surf%zma_s(i)*c
            z_t = z_t + mm*(surf%zmn(i)*c - surf%zma(i)*sn)
            z_p = z_p + nn*(surf%zma(i)*sn - surf%zmn(i)*c)
            v_s = v_s + surf%vmn_s(i)*sn + surf%vma_s(i)*c
            v_t = v_t + mm*(surf%vmn(i)*c - surf%vma(i)*sn)
            v_p = v_p + nn*(surf%vma(i)*sn - surf%vmn(i)*c)
        end do
        shift = twopi/surf%nper
        p%phi = phi_b - shift*v
        call boozer_metric(p%r, r_s, r_t, r_p, z_s, z_t, z_p, &
            -shift*v_s, -shift*v_t, 1 - shift*v_p, &
            p%jacobian, p%grad_s)
    end function boozer_point

    !> Ratio of the surface integrals of the geometric Jacobian
    !! X_s.(X_theta x X_phiB) and of the Boozer flux-function Jacobian
    !! psi'(G + iota I)/B^2. It differs from 1 where the file's geometry and
    !! currents disagree. NEO-2 averages with the flux-function form, so its
    !! <|grad s|> is this ratio times the geometric one.
    function boozer_jacobian_ratio(surf, n_theta) result(ratio)
        type(boozer_surface_t), intent(in) :: surf
        integer, intent(in) :: n_theta
        real(dp) :: ratio, geometric, flux
        type(boozer_point_t) :: p
        integer :: it
        geometric = 0; flux = 0
        do it = 0, n_theta - 1
            p = boozer_point(surf, twopi*it/n_theta, 0.0_dp)
            geometric = geometric + abs(p%jacobian)
            flux = flux + abs(surf%psi_prime*(surf%g_cov + surf%iota*surf%i_cov)) &
                /p%bmod**2
        end do
        ratio = geometric/flux
    end function boozer_jacobian_ratio

    !> Flux-surface averages <|grad s|> and <B^2> with the Jacobian weight on
    !! a uniform (theta, phi_B) grid over one field period.
    subroutine boozer_average(surf, n_theta, n_phi, grad_s, b2)
        type(boozer_surface_t), intent(in) :: surf
        integer, intent(in) :: n_theta, n_phi
        real(dp), intent(out) :: grad_s, b2
        type(boozer_point_t) :: p
        real(dp) :: weight, total
        integer :: it, ip

        grad_s = 0; b2 = 0; total = 0
        do ip = 0, n_phi - 1
            do it = 0, n_theta - 1
                p = boozer_point(surf, twopi*it/n_theta, &
                    twopi*ip/(n_phi*surf%nper))
                weight = abs(p%jacobian)
                total = total + weight
                grad_s = grad_s + weight*p%grad_s
                b2 = b2 + weight*p%bmod**2
            end do
        end do
        grad_s = grad_s/total
        b2 = b2/total
    end subroutine boozer_average

end module milca_boozer_file
