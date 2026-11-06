!> Kinetic geometry on one surface of a Boozer-coordinate `.bc` file, the
!! input format of NEO-2, for comparisons on the same field.
!!
!! The file is read by libneo's `bc_file` (stellarator-symmetric files). On
!! the Boozer angles (theta, phi_B), with n the file's mode number without
!! the field period N_p,
!!   R = sum rmnc cos(a), Z = sum zmns sin(a), B = sum bmnc cos(a),
!!   phi = phi_B - (2 pi/N_p) sum vmns sin(a),  a = m theta - n N_p phi_B,
!! as in NEO-2 (`neo_sub.f90`). Coefficients divided by s^(m/2) (m capped
!! at 12 or 13) and iota
!! are interpolated in s with not-a-knot cubic splines, as NEO-2 splines
!! c_mn/s^(m/2); a surface outside the tabulated range is refused. Lengths in m, |B| in T, |grad s| in 1/m.
module kin6d_boozer_file
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use bc_file, only: bc_data_t, read_bc_file
    use fortnum_cubic_spline, only: cubic_spline_t, cubic_spline_init, &
        cubic_spline_eval
    use kin6d_boozer_file_leaves, only: boozer_axis_factor, boozer_metric
    implicit none
    private
    public :: boozer_file_t, boozer_surface_t, boozer_point_t
    public :: load_boozer_file, boozer_surface, boozer_point, boozer_average
    public :: boozer_jacobian_ratio
    integer, parameter, public :: BOOZER_OK = 0, BOOZER_OUTSIDE = 1

    real(dp), parameter :: twopi = 8*atan(1.0_dp)
    real(dp), parameter :: mu0 = 2*twopi*1e-7_dp

    type :: boozer_file_t
        type(bc_data_t) :: d
        ! Columns: iota, then rmnc, zmns, vmns, bmnc each divided by
        ! s^(m/2), then the poloidal and toroidal currents.
        type(cubic_spline_t) :: spline
    end type

    type :: boozer_surface_t
        real(dp) :: s = 0, iota = 0, diota_ds = 0
        ! Covariant Boozer components G, I [T m] and dpsi_tor/ds [T m^2].
        real(dp) :: g_cov = 0, i_cov = 0, psi_prime = 0
        integer :: nper = 1
        integer, allocatable :: m(:), n(:)
        real(dp), allocatable :: rmn(:), zmn(:), vmn(:), bmn(:)
        real(dp), allocatable :: rmn_s(:), zmn_s(:), vmn_s(:)
    end type

    type :: boozer_point_t
        real(dp) :: r = 0, z = 0, phi = 0, bmod = 0
        real(dp) :: jacobian = 0 ! X_s . (X_theta x X_phiB), m^3
        real(dp) :: grad_s = 0 ! |grad s|, 1/m
    end type

contains

    subroutine load_boozer_file(filename, bf, status)
        character(*), intent(in) :: filename
        type(boozer_file_t), intent(out) :: bf
        integer, intent(out) :: status
        real(dp), allocatable :: y(:, :)
        real(dp) :: f, f_s
        integer :: k, i, nm

        call read_bc_file(filename, bf%d)
        nm = bf%d%nmode
        allocate (y(bf%d%nsurf, 3 + 4*nm))
        do k = 1, bf%d%nsurf
            y(k, 1) = bf%d%iota(k)
            do i = 1, nm
                call boozer_axis_factor(bf%d%s(k), half_m(bf%d%m(i)), f, f_s)
                y(k, 1 + i) = bf%d%rmnc(k, i)/f
                y(k, 1 + nm + i) = bf%d%zmns(k, i)/f
                y(k, 1 + 2*nm + i) = bf%d%vmns(k, i)/f
                y(k, 1 + 3*nm + i) = bf%d%bmnc(k, i)/f
            end do
            y(k, 2 + 4*nm) = bf%d%Jpol_over_nper(k)*bf%d%nper
            y(k, 3 + 4*nm) = bf%d%Itor(k)
        end do
        call cubic_spline_init(bf%d%s, y, bf%spline, status)
        if (status /= 0) status = BOOZER_OUTSIDE
    end subroutine load_boozer_file

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
        allocate (g(3 + 4*nm), g_s(3 + 4*nm))
        call cubic_spline_eval(bf%spline, s, g, g_s)
        surf%s = s
        surf%nper = bf%d%nper
        surf%m = bf%d%m
        surf%n = bf%d%n
        surf%iota = g(1)
        surf%diota_ds = g_s(1)
        surf%g_cov = mu0*g(2 + 4*nm)/twopi
        surf%i_cov = mu0*g(3 + 4*nm)/twopi
        surf%psi_prime = bf%d%flux/twopi
        allocate (surf%rmn(nm), surf%zmn(nm), surf%vmn(nm), surf%bmn(nm), &
            surf%rmn_s(nm), surf%zmn_s(nm), surf%vmn_s(nm))
        do i = 1, nm
            call boozer_axis_factor(s, half_m(bf%d%m(i)), f, f_s)
            surf%rmn(i) = f*g(1 + i)
            surf%zmn(i) = f*g(1 + nm + i)
            surf%vmn(i) = f*g(1 + 2*nm + i)
            surf%bmn(i) = f*g(1 + 3*nm + i)
            surf%rmn_s(i) = f_s*g(1 + i) + f*g_s(1 + i)
            surf%zmn_s(i) = f_s*g(1 + nm + i) + f*g_s(1 + nm + i)
            surf%vmn_s(i) = f_s*g(1 + 2*nm + i) + f*g_s(1 + 2*nm + i)
        end do
        status = BOOZER_OK
    end subroutine boozer_surface

    ! Axis exponent m/2, capped at 6 or 6.5 (keeping the parity of m) for
    ! m > 12 as in NEO-2 (`neo_sub.f90`, m_max_sp), so high harmonics are not
    ! divided by large powers of s.
    pure real(dp) function half_m(m)
        integer, intent(in) :: m
        if (m <= 12) then
            half_m = 0.5_dp*m
        else
            half_m = 0.5_dp*(12 + modulo(m, 2))
        end if
    end function half_m

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
            p%r = p%r + surf%rmn(i)*c
            p%z = p%z + surf%zmn(i)*sn
            p%bmod = p%bmod + surf%bmn(i)*c
            v = v + surf%vmn(i)*sn
            r_s = r_s + surf%rmn_s(i)*c
            r_t = r_t - mm*surf%rmn(i)*sn
            r_p = r_p + nn*surf%rmn(i)*sn
            z_s = z_s + surf%zmn_s(i)*sn
            z_t = z_t + mm*surf%zmn(i)*c
            z_p = z_p - nn*surf%zmn(i)*c
            v_s = v_s + surf%vmn_s(i)*sn
            v_t = v_t + mm*surf%vmn(i)*c
            v_p = v_p - nn*surf%vmn(i)*c
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

end module kin6d_boozer_file
