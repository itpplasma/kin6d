module kin6d_coils_3d
    !! K-C1: closed 3D polyline coils, closed-form straight-segment Biot-Savart
    !! field, and projection of the field onto toroidal Fourier modes in the
    !! meridional plane.
    !!
    !! Conventions. Lengths in m, currents in A, field in T. Coil points are
    !! Cartesian (x, y, z) with z up and phi = atan2(y, x), so (R, phi, Z) is
    !! right handed. Each polyline is a closed loop: the last point connects
    !! back to the first (a duplicated closing point is dropped on input).
    !! The current flows from point i to point i+1.
    !!
    !! Fourier convention (two-sided):
    !!     B(R, phi, Z) = sum_n bn(n; R, Z) exp(i n phi),
    !!     bn(n)        = (1/Nphi) sum_j B(R, phi_j, Z) exp(-i n phi_j),
    !! with B_{-n} = conjg(B_n). The physical n-th harmonic has peak
    !! amplitude 2 |B_n| for n > 0. Components are cylindrical (R, phi, Z).
    !! The DFT aliases n -> n + k Nphi; the error is exp(-Nphi d / R) with d
    !! the distance to the nearest coil.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    private

    real(dp), parameter, public :: coil_mu0 = 4.0e-7_dp*3.14159265358979323846_dp
    real(dp), parameter :: pi = 3.14159265358979323846_dp
    integer, parameter, public :: coil_name_len = 32

    type, public :: coil_set_t
        integer :: ncoil = 0
        integer :: nseg = 0
        integer, allocatable :: first(:)  ! segments of coil k: first(k):first(k+1)-1
        real(dp), allocatable :: a(:, :)        ! segment start (3, nseg)
        real(dp), allocatable :: l(:, :)        ! segment vector end - start (3, nseg)
        character(len=coil_name_len), allocatable :: name(:)
    end type

    public :: coil_set_from_polylines, coil_set_read, coil_set_rotated_copy
    public :: coil_read_currents
    public :: coil_set_b_cart, coil_set_modes, coil_synthesis, segment_b

contains

    pure function cross(u, v) result(w)
        real(dp), intent(in) :: u(3), v(3)
        real(dp) :: w(3)
        w = [u(2)*v(3) - u(3)*v(2), u(3)*v(1) - u(1)*v(3), u(1)*v(2) - u(2)*v(1)]
    end function

    pure function segment_b(a, l, p) result(b)
        !! Field in T per ampere of the straight segment a -> a + l at p:
        !! B = mu0/(4 pi) (lh x r1)/|lh x r1|^2 (lh.r1/|r1| - lh.r2/|r2|).
        real(dp), intent(in) :: a(3), l(3), p(3)
        real(dp) :: b(3)
        real(dp) :: r1(3), r2(3), lh(3), c(3), c2, n1, n2, len

        r1 = p - a
        r2 = r1 - l
        len = norm2(l)
        lh = l/len
        c = cross(lh, r1)
        c2 = dot_product(c, c)
        n1 = norm2(r1)
        n2 = norm2(r2)
        if (c2 <= 1.0e-24_dp*len*len .or. n1 == 0 .or. n2 == 0) then
            b = 0     ! on the supporting line: field of the segment is zero or singular
            return
        end if
        b = coil_mu0/(4*pi)*c/c2*(dot_product(lh, r1)/n1 - dot_product(lh, r2)/n2)
    end function

    subroutine coil_set_from_polylines(set, npts, pts, names, status)
        !! pts(3, sum(npts)) holds the coils one after another.
        type(coil_set_t), intent(out) :: set
        integer, intent(in) :: npts(:)
        real(dp), intent(in) :: pts(:, :)
        character(len=*), intent(in), optional :: names(:)
        integer, intent(out) :: status
        integer :: k, i, n, off, iseg, nclosed

        status = 1
        if (size(pts, 1) /= 3 .or. sum(npts) /= size(pts, 2) .or. any(npts < 3)) return
        set%ncoil = size(npts)
        allocate (set%first(set%ncoil + 1), set%name(set%ncoil))
        set%name = ''
        if (present(names)) then
            if (size(names) /= set%ncoil) return
            set%name = names
        end if
        ! drop a duplicated closing point; count segments
        nclosed = 0
        off = 0
        do k = 1, set%ncoil
            n = npts(k)
            if (norm2(pts(:, off + 1) - pts(:, off + n)) == 0) n = n - 1
            nclosed = nclosed + n
            off = off + npts(k)
        end do
        set%nseg = nclosed
        allocate (set%a(3, nclosed), set%l(3, nclosed))
        iseg = 0
        off = 0
        do k = 1, set%ncoil
            n = npts(k)
            if (norm2(pts(:, off + 1) - pts(:, off + n)) == 0) n = n - 1
            set%first(k) = iseg + 1
            do i = 1, n
                iseg = iseg + 1
                set%a(:, iseg) = pts(:, off + i)
                if (i < n) then
                    set%l(:, iseg) = pts(:, off + i + 1) - pts(:, off + i)
                else
                    set%l(:, iseg) = pts(:, off + 1) - pts(:, off + n)
                end if
                if (norm2(set%l(:, iseg)) == 0) return
            end do
            off = off + npts(k)
        end do
        set%first(set%ncoil + 1) = iseg + 1
        status = 0
    end subroutine

    subroutine coil_set_read(set, path, status)
        !! Neutral coil table: "ncoil"; then per coil "name npts" followed by
        !! npts lines "x y z" in m. Lines starting with '#' are comments.
        type(coil_set_t), intent(out) :: set
        character(len=*), intent(in) :: path
        integer, intent(out) :: status
        integer :: unit, ios, k, i, ncoil, total
        integer, allocatable :: npts(:)
        real(dp), allocatable :: pts(:, :)
        character(len=coil_name_len), allocatable :: names(:)
        character(len=256) :: line

        status = 1
        open (newunit=unit, file=path, status='old', action='read', iostat=ios)
        if (ios /= 0) return
        call next_data_line(unit, line, ios)
        if (ios /= 0) return
        read (line, *, iostat=ios) ncoil
        if (ios /= 0 .or. ncoil < 1) return
        allocate (npts(ncoil), names(ncoil), pts(3, 0))
        total = 0
        do k = 1, ncoil
            call next_data_line(unit, line, ios)
            if (ios /= 0) return
            read (line, *, iostat=ios) names(k), npts(k)
            if (ios /= 0 .or. npts(k) < 3) return
            call grow(pts, total + npts(k))
            do i = 1, npts(k)
                call next_data_line(unit, line, ios)
                if (ios /= 0) return
                read (line, *, iostat=ios) pts(:, total + i)
                if (ios /= 0) return
            end do
            total = total + npts(k)
        end do
        close (unit)
        call coil_set_from_polylines(set, npts, pts(:, 1:total), names, status)
    contains
        subroutine grow(x, n)
            real(dp), allocatable, intent(inout) :: x(:, :)
            integer, intent(in) :: n
            real(dp), allocatable :: tmp(:, :)
            allocate (tmp(3, n))
            tmp(:, 1:size(x, 2)) = x
            call move_alloc(tmp, x)
        end subroutine
    end subroutine

    subroutine next_data_line(unit, line, ios)
        integer, intent(in) :: unit
        character(len=*), intent(out) :: line
        integer, intent(out) :: ios
        do
            read (unit, '(a)', iostat=ios) line
            if (ios /= 0) return
            if (len_trim(line) > 0 .and. line(1:1) /= '#') return
        end do
    end subroutine

    subroutine coil_read_currents(path, ncoil, currents, status)
        !! One current in A per coil, in coil-table order ('#' comments allowed).
        character(len=*), intent(in) :: path
        integer, intent(in) :: ncoil
        real(dp), intent(out) :: currents(ncoil)
        integer, intent(out) :: status
        integer :: unit, ios
        character(len=4096) :: line
        status = 1
        open (newunit=unit, file=path, status='old', action='read', iostat=ios)
        if (ios /= 0) return
        call next_data_line(unit, line, ios)
        if (ios == 0) read (line, *, iostat=ios) currents
        close (unit)
        if (ios == 0) status = 0
    end subroutine

    subroutine coil_set_rotated_copy(src, k, angle, dst, status)
        !! Coil k of src rotated about the z axis by angle (rad), as a 1-coil set.
        type(coil_set_t), intent(in) :: src
        integer, intent(in) :: k
        real(dp), intent(in) :: angle
        type(coil_set_t), intent(out) :: dst
        integer, intent(out) :: status
        real(dp), allocatable :: pts(:, :)
        integer :: n, i, s
        n = src%first(k + 1) - src%first(k)
        allocate (pts(3, n))
        do i = 1, n
            s = src%first(k) + i - 1
            pts(:, i) = [cos(angle)*src%a(1, s) - sin(angle)*src%a(2, s), &
                         sin(angle)*src%a(1, s) + cos(angle)*src%a(2, s), src%a(3, s)]
        end do
        call coil_set_from_polylines(dst, [n], pts, status=status)
    end subroutine

    subroutine coil_set_b_cart(set, weights, p, b)
        !! Direct 3D Biot-Savart field (T) of the weighted coils (A) at p (m).
        type(coil_set_t), intent(in) :: set
        real(dp), intent(in) :: weights(:), p(3)
        real(dp), intent(out) :: b(3)
        integer :: k, s
        b = 0
        do k = 1, set%ncoil
            if (weights(k) == 0) cycle
            do s = set%first(k), set%first(k + 1) - 1
                b = b + weights(k)*segment_b(set%a(:, s), set%l(:, s), p)
            end do
        end do
    end subroutine

    subroutine coil_set_modes(set, weights, rz, nphi, modes, bn, status)
        !! Toroidal Fourier components bn(c, im, ip) of the cylindrical field
        !! components c = (R, phi, Z) at the meridional points rz(:, ip) = (R, Z).
        !! Every point is independent, so the result does not depend on the
        !! number of OpenMP threads.
        type(coil_set_t), intent(in) :: set
        real(dp), intent(in) :: weights(:), rz(:, :)
        integer, intent(in) :: nphi, modes(:)
        complex(dp), intent(out) :: bn(:, :, :)
        integer, intent(out) :: status
        real(dp), allocatable :: cp(:), sp(:)
        complex(dp), allocatable :: ph(:, :)
        integer, allocatable :: act(:)
        integer :: j, im, ip, nm, nact, k
        real(dp) :: bc(3), p(3), r, z, bcyl(3), arg
        complex(dp) :: acc(3, size(modes))

        status = 1
        nm = size(modes)
        if (size(weights) /= set%ncoil .or. size(rz, 1) /= 2 .or. nphi < 2) return
        if (any(shape(bn) /= [3, nm, size(rz, 2)])) return
        allocate (cp(nphi), sp(nphi), ph(nphi, nm))
        do j = 1, nphi
            cp(j) = cos(2*pi*(j - 1)/nphi)
            sp(j) = sin(2*pi*(j - 1)/nphi)
        end do
        do im = 1, nm
            do j = 1, nphi
            arg = 2*pi*real(mod(int(modes(im), 8)*int(j - 1, 8), int(nphi, 8)), dp)/nphi
                ph(j, im) = cmplx(cos(arg), -sin(arg), dp)
            end do
        end do
        nact = count(weights /= 0)
        allocate (act(nact))
        act = pack([(k, k=1, set%ncoil)], weights /= 0)
        !$omp parallel do default(none) schedule(dynamic, 1) &
        !$omp shared(set, weights, rz, nphi, nm, cp, sp, ph, bn, act, nact) &
        !$omp private(ip, j, im, k, p, r, z, bc, bcyl, acc)
        do ip = 1, size(rz, 2)
            r = rz(1, ip)
            z = rz(2, ip)
            acc = 0
            do j = 1, nphi
                p = [r*cp(j), r*sp(j), z]
                bc = 0
                do k = 1, nact
                    bc = bc + weights(act(k))*coil_b_one(set, act(k), p)
                end do
                bcyl = [bc(1)*cp(j) + bc(2)*sp(j), -bc(1)*sp(j) + bc(2)*cp(j), bc(3)]
                do im = 1, nm
                    acc(:, im) = acc(:, im) + bcyl*ph(j, im)
                end do
            end do
            bn(:, :, ip) = acc/nphi
        end do
        !$omp end parallel do
        status = 0
    end subroutine

    pure function coil_b_one(set, k, p) result(b)
        type(coil_set_t), intent(in) :: set
        integer, intent(in) :: k
        real(dp), intent(in) :: p(3)
        real(dp) :: b(3)
        integer :: s
        b = 0
        do s = set%first(k), set%first(k + 1) - 1
            b = b + segment_b(set%a(:, s), set%l(:, s), p)
        end do
    end function

    pure function coil_synthesis(bn, modes, phi) result(b)
        !! Real field at toroidal angle phi from the components of the modes
        !! n >= 0 (B_{-n} = conjg(B_n)): Re B_0 + 2 Re sum_{n>0} B_n exp(i n phi).
        complex(dp), intent(in) :: bn(:, :)
        integer, intent(in) :: modes(:)
        real(dp), intent(in) :: phi
        real(dp) :: b(3)
        integer :: im
        b = 0
        do im = 1, size(modes)
            if (modes(im) < 0) cycle
            b = b + merge(1.0_dp, 2.0_dp, modes(im) == 0)* &
                real(bn(:, im)*cmplx(cos(modes(im)*phi), sin(modes(im)*phi), dp), dp)
        end do
    end function

end module
