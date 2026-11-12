module coil_box_geometry
    !! Star-shaped wall given as a Fourier surface, R = sum rc cos(m theta -
    !! n P phi), Z = sum zs sin(m theta - n P phi), and the map of the
    !! reference domain (the phi-mean cross-section) onto the cross-section
    !! at angle phi: equal star coordinates (u, theta) about the area
    !! centroids, y = c + u (s w + (1 - s) A w_ref) with w, w_ref relative to
    !! the centroids, A the least-squares linear map w_ref -> w and s(u) a
    !! cubic ramp from 0 at u = 0.1 to 1 at u = 0.6. The star map alone
    !! (s = 1) has a direction-dependent slope at the centre, which no
    !! polynomial element there represents; the affine core removes it. The
    !! map cannot fold while every cross-section is star-shaped about its
    !! centroid in theta (star_quality) and A stays close to the shape.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use milca_coils_3d, only: coil_set_t, coil_set_a_cart, coil_set_b_cart
    implicit none
    private
    public :: wall_read, wall_point, wall_centroid, star_quality, ref_polar, ref_centre
    public :: box_map, box_potential, box_field, set_map_planes, box_coordinates
    public :: coils, weights, nfp
    public :: box_source, axis_rc, axis_zs, source_peak, source_width

    real(dp), parameter :: pi = acos(-1.0_dp)
    integer, parameter :: nsample = 1024
    real(dp), allocatable :: wn(:), wm(:), wr(:), wz(:)
    integer :: nfp = 1
    type(coil_set_t) :: coils
    real(dp), allocatable :: weights(:)
    real(dp) :: cref(2)
    real(dp), allocatable :: plane_phi(:), plane_c(:, :), plane_a(:, :, :)
    ! Heating source of the finite-pressure box: a Gaussian of width w about
    ! the axis R = sum rc cos(n P phi), Z = sum zs sin(n P phi), in W/m^3
    ! per unit chi_perp (the steady p equation is -div(K grad p) = S).
    real(dp) :: axis_rc(0:16) = 0, axis_zs(0:16) = 0, source_peak = 0, &
        source_width = 0.1_dp
    real(dp), parameter :: ramp0 = 0.1_dp, ramp1 = 0.6_dp

contains

    subroutine wall_read(path, period, status)
        character(len=*), intent(in) :: path
        integer, intent(in) :: period
        integer, intent(out) :: status
        integer :: unit, ios, n
        real(dp) :: row(4)
        status = 1
        open (newunit=unit, file=path, status='old', action='read', iostat=ios)
        if (ios /= 0) return
        n = 0
        do
            read (unit, *, iostat=ios) row
            if (ios /= 0) exit
            n = n + 1
        end do
        if (n < 1) return
        rewind (unit)
        allocate (wn(n), wm(n), wr(n), wz(n))
        do n = 1, size(wn)
            read (unit, *) wn(n), wm(n), wr(n), wz(n)
        end do
        close (unit)
        nfp = period
        call wall_centroid(0.0_dp, cref, mean=.true.)
        status = 0
    end subroutine

    pure function wall_point(theta, phi, mean) result(w)
        ! mean: the phi-mean cross-section, the n = 0 harmonics.
        real(dp), intent(in) :: theta, phi
        logical, intent(in) :: mean
        real(dp) :: w(2), a(size(wn))
        a = wm*theta - wn*nfp*phi
        if (mean) then
            w = [sum(merge(wr, 0.0_dp, wn == 0)*cos(wm*theta)), &
                 sum(merge(wz, 0.0_dp, wn == 0)*sin(wm*theta))]
        else
            w = [sum(wr*cos(a)), sum(wz*sin(a))]
        end if
    end function

    pure function wall_tangent(theta, phi, mean) result(t)
        real(dp), intent(in) :: theta, phi
        logical, intent(in) :: mean
        real(dp) :: t(2), a(size(wn))
        a = wm*theta - wn*nfp*phi
        if (mean) then
            t = [-sum(merge(wm*wr, 0.0_dp, wn == 0)*sin(wm*theta)), &
                 sum(merge(wm*wz, 0.0_dp, wn == 0)*cos(wm*theta))]
        else
            t = [-sum(wm*wr*sin(a)), sum(wm*wz*cos(a))]
        end if
    end function

    pure subroutine wall_centroid(phi, c, mean)
        real(dp), intent(in) :: phi
        real(dp), intent(out) :: c(2)
        logical, intent(in) :: mean
        real(dp) :: p(2), q(2), cr, area
        integer :: k
        area = 0; c = 0
        do k = 0, nsample - 1
            p = wall_point(2*pi*k/nsample, phi, mean)
            q = wall_point(2*pi*(k + 1)/nsample, phi, mean)
            cr = p(1)*q(2) - q(1)*p(2)
            area = area + cr/2
            c = c + (p + q)*cr/6
        end do
        c = c/area
    end subroutine

    real(dp) function star_quality(nphi) result(q)
        ! Minimum over theta and nphi planes of sin(angle between the ray
        ! from the centroid and the wall tangent); positive iff star-shaped.
        integer, intent(in) :: nphi
        real(dp) :: c(2), p(2), t(2), phi
        integer :: j, k
        q = huge(1.0_dp)
        do j = 0, nphi
            phi = 2*pi*j/(nphi*nfp)
            if (j < nphi) then
                call wall_centroid(phi, c, .false.)
            else
                c = cref
            end if
            do k = 0, nsample - 1
                p = wall_point(2*pi*k/nsample, phi, j == nphi) - c
                t = wall_tangent(2*pi*k/nsample, phi, j == nphi)
                q = min(q, (p(1)*t(2) - p(2)*t(1))/(norm2(p)*norm2(t)))
            end do
        end do
    end function

    pure function ref_centre() result(c)
        real(dp) :: c(2)
        c = cref
    end function

    subroutine ref_polar(u, theta, x)
        ! Reference mesh: star scaling of the phi-mean cross-section.
        real(dp), intent(in) :: u, theta
        real(dp), intent(out) :: x(2)
        x = cref + u*(wall_point(theta, 0.0_dp, .true.) - cref)
    end subroutine

    subroutine set_map_planes(nphi)
        ! Centroids of the planes on which nn_map_build samples the map.
        integer, intent(in) :: nphi
        integer :: j, k
        real(dp) :: g(2, 2), h(2, 2), wr(2), w(2), det
        if (allocated(plane_phi)) deallocate (plane_phi, plane_c, plane_a)
        allocate (plane_phi(nphi), plane_c(2, nphi), plane_a(2, 2, nphi))
        do j = 1, nphi
            plane_phi(j) = 2*pi*(j - 1)/(nphi*nfp)
            call wall_centroid(plane_phi(j), plane_c(:, j), .false.)
            ! A = H G^-1 with G = sum wr wr^T, H = sum w wr^T
            g = 0; h = 0
            do k = 0, nsample - 1
                wr = wall_point(2*pi*k/nsample, 0.0_dp, .true.) - cref
                w = wall_point(2*pi*k/nsample, plane_phi(j), .false.) - plane_c(:, j)
                g = g + spread(wr, 2, 2)*spread(wr, 1, 2)
                h = h + spread(w, 2, 2)*spread(wr, 1, 2)
            end do
            det = g(1, 1)*g(2, 2) - g(1, 2)*g(2, 1)
            plane_a(:, :, j) = matmul(h, reshape([g(2, 2), -g(2, 1), -g(1, 2), g(1, 1)], &
                [2, 2])/det)
        end do
    end subroutine

    subroutine star_coordinates(d, u, theta)
        ! Star coordinates of the offset d from cref on the mean wall:
        ! d = u (w_ref(theta) - cref), u >= 0.
        real(dp), intent(in) :: d(2)
        real(dp), intent(out) :: u, theta
        real(dp) :: f, df, w(2), t(2)
        integer :: k, it
        theta = 0; f = huge(1.0_dp)
        do k = 0, 255
            w = wall_point(2*pi*k/256, 0.0_dp, .true.) - cref
            if (dot_product(w, d) > 0 .and. abs(w(1)*d(2) - w(2)*d(1))/norm2(w) < f) then
                f = abs(w(1)*d(2) - w(2)*d(1))/norm2(w); theta = 2*pi*k/256
            end if
        end do
        do it = 1, 50
            w = wall_point(theta, 0.0_dp, .true.) - cref
            t = wall_tangent(theta, 0.0_dp, .true.)
            f = w(1)*d(2) - w(2)*d(1)
            df = t(1)*d(2) - t(2)*d(1)
            theta = theta - f/df
            if (abs(f) <= 1.0e-15_dp*norm2(w)*norm2(d)) exit
        end do
        w = wall_point(theta, 0.0_dp, .true.) - cref
        if (dot_product(w, d) <= 0) error stop 'star inversion on the wrong ray'
        u = norm2(d)/norm2(w)
        theta = modulo(theta, 2*pi)
    end subroutine

    subroutine box_coordinates(centre, scale, x, u, theta)
        ! Locator coordinates of the reference mesh built from ref_polar.
        real(dp), intent(in) :: centre(2), scale, x(2)
        real(dp), intent(out) :: u, theta
        if (norm2(x - cref) <= 1.0e-14_dp) then
            u = 0; theta = 0
            return
        end if
        call star_coordinates(x - cref, u, theta)
    end subroutine

    subroutine box_map(x, phi, y)
        ! Star coordinates of x about cref on the mean wall, then the same
        ! coordinates on the wall at phi.
        real(dp), intent(in) :: x(2), phi
        real(dp), intent(out) :: y(2)
        real(dp) :: d(2), theta, u, c(2), w(2), s, z
        integer :: j
        j = minloc(abs(plane_phi - phi), 1)
        if (abs(plane_phi(j) - phi) > 1.0e-12_dp) error stop 'box_map: phi not a sampling plane'
        c = plane_c(:, j)
        d = x - cref
        if (norm2(d) <= 1.0e-14_dp) then
            y = c
            return
        end if
        if (norm2(d) < 0.5_dp*ramp0*minval(wall_radii())) then
            y = c + matmul(plane_a(:, :, j), d)
            return
        end if
        call star_coordinates(d, u, theta)
        w = wall_point(theta, 0.0_dp, .true.) - cref
        z = min(max((u - ramp0)/(ramp1 - ramp0), 0.0_dp), 1.0_dp)
        s = z*z*(3 - 2*z)
        y = c + u*(s*(wall_point(theta, phi, .false.) - c) &
            + (1 - s)*matmul(plane_a(:, :, j), w))
    end subroutine

    function wall_radii() result(r)
        ! Distances from cref to the mean wall at 64 angles.
        real(dp) :: r(64)
        integer :: k
        do k = 1, 64
            r(k) = norm2(wall_point(2*pi*(k - 1)/64, 0.0_dp, .true.) - cref)
        end do
    end function

    function box_source(x, phi) result(s)
        real(dp), intent(in) :: x(2), phi
        real(dp) :: s, ax(2)
        integer :: n
        ax = 0
        do n = 0, 16
            ax = ax + [axis_rc(n)*cos(n*nfp*phi), axis_zs(n)*sin(n*nfp*phi)]
        end do
        s = source_peak*exp(-sum((x - ax)**2)/source_width**2)
    end function

    subroutine box_potential(x, a)
        ! Cartesian vector potential of the weighted coils (T m).
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: a(3)
        call coil_set_a_cart(coils, weights, x, a)
    end subroutine

    subroutine box_field(x, b)
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: b(3)
        call coil_set_b_cart(coils, weights, x, b)
    end subroutine
end module coil_box_geometry
