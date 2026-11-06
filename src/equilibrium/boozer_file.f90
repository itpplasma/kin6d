!> Kinetic geometry on one surface of a Boozer-coordinate `.bc` file, the
!! input format of NEO-2, for comparisons on the same field.
!!
!! The file is read by libneo's `bc_file` (stellarator-symmetric files). On
!! the Boozer angles (theta, phi_B), with n the file's mode number without
!! the field period N_p,
!!   R = sum rmnc cos(a), Z = sum zmns sin(a), B = sum bmnc cos(a),
!!   phi = phi_B - (2 pi/N_p) sum vmns sin(a),  a = m theta - n N_p phi_B,
!! as in NEO-2 (`neo_sub.f90`). Coefficients and iota are interpolated in s
!! with four-point Lagrange weights; a surface outside the tabulated range
!! is refused. Lengths in m, |B| in T, |grad s| in 1/m.
module kin6d_boozer_file
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use bc_file, only: bc_data_t, read_bc_file
    use kin6d_boozer_file_leaves, only: boozer_lagrange4, boozer_metric
    implicit none
    private
    public :: boozer_file_t, boozer_surface_t, boozer_point_t
    public :: load_boozer_file, boozer_surface, boozer_point, boozer_average
    integer, parameter, public :: BOOZER_OK = 0, BOOZER_OUTSIDE = 1

    real(dp), parameter :: twopi = 8*atan(1.0_dp)

    type :: boozer_file_t
        type(bc_data_t) :: d
    end type

    type :: boozer_surface_t
        real(dp) :: s = 0, iota = 0, diota_ds = 0
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

    subroutine load_boozer_file(filename, bf)
        character(*), intent(in) :: filename
        type(boozer_file_t), intent(out) :: bf
        call read_bc_file(filename, bf%d)
    end subroutine load_boozer_file

    subroutine boozer_surface(bf, s, surf, status)
        type(boozer_file_t), intent(in) :: bf
        real(dp), intent(in) :: s
        type(boozer_surface_t), intent(out) :: surf
        integer, intent(out) :: status
        integer :: ns, k, i
        real(dp) :: w(4), dw(4)

        ns = bf%d%nsurf
        status = BOOZER_OUTSIDE
        if (ns < 4) return
        if (s < bf%d%s(1) .or. s > bf%d%s(ns)) return
        k = 1
        do while (k < ns - 3)
            if (s < bf%d%s(k + 2)) exit
            k = k + 1
        end do
        call boozer_lagrange4(s, bf%d%s(k), bf%d%s(k + 1), bf%d%s(k + 2), &
            bf%d%s(k + 3), w(1), w(2), w(3), w(4), &
            dw(1), dw(2), dw(3), dw(4))
        surf%s = s
        surf%nper = bf%d%nper
        surf%m = bf%d%m
        surf%n = bf%d%n
        surf%iota = dot_product(w, bf%d%iota(k:k + 3))
        surf%diota_ds = dot_product(dw, bf%d%iota(k:k + 3))
        allocate (surf%rmn(bf%d%nmode), surf%zmn(bf%d%nmode), &
            surf%vmn(bf%d%nmode), surf%bmn(bf%d%nmode), &
            surf%rmn_s(bf%d%nmode), surf%zmn_s(bf%d%nmode), &
            surf%vmn_s(bf%d%nmode))
        do i = 1, bf%d%nmode
            surf%rmn(i) = dot_product(w, bf%d%rmnc(k:k + 3, i))
            surf%zmn(i) = dot_product(w, bf%d%zmns(k:k + 3, i))
            surf%vmn(i) = dot_product(w, bf%d%vmns(k:k + 3, i))
            surf%bmn(i) = dot_product(w, bf%d%bmnc(k:k + 3, i))
            surf%rmn_s(i) = dot_product(dw, bf%d%rmnc(k:k + 3, i))
            surf%zmn_s(i) = dot_product(dw, bf%d%zmns(k:k + 3, i))
            surf%vmn_s(i) = dot_product(dw, bf%d%vmns(k:k + 3, i))
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
