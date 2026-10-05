module kin6d_debye_sphere_reference
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortnum_interval, only: interval_t, interval, hull, is_empty, iadd, &
        linear_interpolation_remainder
    use kin6d_debye_sphere_leaves, only: sphere_inside, sphere_outside, &
        sphere_inside_tangent, sphere_outside_tangent, sphere_inside_packet, &
        sphere_outside_packet, sphere_bare_inside, sphere_bare_outside, &
        sphere_bare_inside_interval, sphere_bare_outside_interval, &
        sphere_bare_curvature
    implicit none
    private
    public :: debye_sphere_values, debye_sphere_enclosure, debye_sphere_bands

contains
    ! Analytical matched-equilibrium reference. Fixed R, infinity gauge.
    ! No numerical BVP solve, no kinetic transient, no model-error estimate.
    subroutine debye_sphere_values(r, k, phi, bare, tangent, ok)
        real(dp), intent(in) :: r, k
        real(dp), intent(out) :: phi, bare, tangent
        logical, intent(out) :: ok

        ok = .false.
        if (.not. ieee_is_finite(r)) return
        if (.not. ieee_is_finite(k)) return
        if (r < 0.0_dp .or. r > 6.0_dp) return
        if (k < 0.5_dp .or. k > 2.0_dp) return
        if (r <= 1.0_dp) then
            call sphere_inside(r, k, phi)
            call sphere_inside_tangent(r, k, tangent)
            call sphere_bare_inside(r, bare)
        else
            call sphere_outside(r, k, phi)
            call sphere_outside_tangent(r, k, tangent)
            call sphere_bare_outside(r, bare)
        end if
        ok = all(ieee_is_finite([phi, bare, tangent]))
    end subroutine debye_sphere_values

    recursive subroutine debye_sphere_enclosure(r, k, phi, bare, tangent, &
                                                m, mb, mt, ok)
        type(interval_t), intent(in) :: r, k
        type(interval_t), intent(out) :: phi, bare, tangent, m, mb, mt
        logical, intent(out) :: ok
        type(interval_t) :: f, g, ef, eg
        type(interval_t) :: f2, b2, g2, m2, mb2, mt2
        logical :: right_ok

        ok = .false.
        if (.not. valid_box(r)) return
        if (.not. valid_box(k)) return
        if (r%lo < 0.0_dp .or. r%hi > 6.0_dp) return
        if (k%lo < 0.5_dp .or. k%hi > 2.0_dp) return
        if (r%lo < 1.0_dp .and. r%hi > 1.0_dp) then
            call debye_sphere_enclosure(interval(r%lo, 1.0_dp), k, &
                phi, bare, tangent, m, mb, mt, right_ok)
            if (.not. right_ok) return
            call debye_sphere_enclosure(interval(1.0_dp, r%hi), k, &
                f2, b2, g2, m2, mb2, mt2, right_ok)
            if (.not. right_ok) return
            phi = hull(phi, f2); bare = hull(bare, b2)
            tangent = hull(tangent, g2); m = hull(m, m2)
            mb = hull(mb, mb2); mt = hull(mt, mt2)
        else if (r%hi <= 1.0_dp) then
            call sphere_inside_packet(r, k, f, g, m, mt, ef, eg)
            if (.not. valid_box(ef)) return
            if (.not. valid_box(eg)) return
            phi = widen(f, max(0.0_dp, ef%hi))
            tangent = widen(g, max(0.0_dp, eg%hi))
            call sphere_bare_inside_interval(r, bare)
            mb = interval(1.0_dp)
        else
            call sphere_outside_packet(r, k, phi, tangent, m, mt)
            call sphere_bare_outside_interval(r, bare)
            call sphere_bare_curvature(r, mb)
        end if
        if (.not. valid_box(m)) return
        if (.not. valid_box(mt)) return
        if (.not. valid_box(mb)) return
        m = interval(0.0_dp, max(0.0_dp, m%hi))
        mt = interval(0.0_dp, max(0.0_dp, mt%hi))
        mb = interval(0.0_dp, max(0.0_dp, mb%hi))
        ok = valid_box(phi) .and. valid_box(bare) .and. valid_box(tangent)
    end subroutine debye_sphere_enclosure

    ! Endpoint strict intervals plus max adjacent C2 chord errors cover every
    ! radial segment. Crossing r=1 is allowed: phi and d_k phi are C1 there,
    ! with bounded piecewise second derivative; the W^{2,infinity} chord
    ! theorem therefore applies across the interface as well.
    ! All buffers are caller-preallocated. No per-call allocation.
    subroutine debye_sphere_bands(r, k, phi, bare, tangent, lower, upper, work, ok)
        real(dp), intent(in) :: r(:), k
        real(dp), intent(out) :: phi(:), bare(:), tangent(:)
        real(dp), intent(out) :: lower(:, :), upper(:, :), work(:, :)
        logical, intent(out) :: ok
        type(interval_t) :: f, b, g, m, mb, mt, pad, lo
        real(dp) :: padding
        logical :: local_ok
        integer :: i, q, n

        ok = .false.
        n = size(r)
        if (n < 2) return
        if (size(phi) /= n .or. size(bare) /= n .or. size(tangent) /= n) return
        if (size(lower, 1) /= 3 .or. size(upper, 1) /= 3) return
        if (size(lower, 2) /= n .or. size(upper, 2) /= n) return
        if (size(work, 1) /= 3 .or. size(work, 2) /= n - 1) return
        if (.not. all(ieee_is_finite(r))) return
        do i = 2, n
            if (r(i) <= r(i - 1)) return
        end do
        do i = 1, n
            call debye_sphere_values(r(i), k, phi(i), bare(i), tangent(i), local_ok)
            if (.not. local_ok) return
            call debye_sphere_enclosure(interval(r(i)), interval(k), &
                f, b, g, m, mb, mt, local_ok)
            if (.not. local_ok) return
            lower(:, i) = [f%lo, b%lo, g%lo]
            upper(:, i) = [f%hi, b%hi, g%hi]
        end do
        do i = 1, n - 1
            call debye_sphere_enclosure(interval(r(i), r(i + 1)), interval(k), &
                f, b, g, m, mb, mt, local_ok)
            if (.not. local_ok) return
            pad = linear_interpolation_remainder(interval(r(i)), interval(r(i + 1)), m)
            if (.not. valid_box(pad)) return
            work(1, i) = pad%hi
            pad = linear_interpolation_remainder(interval(r(i)), interval(r(i + 1)), mb)
            if (.not. valid_box(pad)) return
            work(2, i) = pad%hi
            pad = linear_interpolation_remainder(interval(r(i)), interval(r(i + 1)), mt)
            if (.not. valid_box(pad)) return
            work(3, i) = pad%hi
        end do
        do i = 1, n
            do q = 1, 3
                padding = work(q, max(1, i - 1))
                if (i < n) padding = max(padding, work(q, i))
                lo = widen(interval(lower(q, i), upper(q, i)), padding)
                if (.not. valid_box(lo)) return
                lower(q, i) = lo%lo; upper(q, i) = lo%hi
            end do
        end do
        ok = .true.
    end subroutine debye_sphere_bands

    pure logical function valid_box(box)
        type(interval_t), intent(in) :: box

        valid_box = .false.
        if (is_empty(box)) return
        if (.not. ieee_is_finite(box%lo)) return
        if (.not. ieee_is_finite(box%hi)) return
        valid_box = .true.
    end function valid_box
    !> Outward-rounded symmetric widening [lo - radius, hi + radius].
    elemental function widen(value, radius) result(w)
        type(interval_t), intent(in) :: value
        real(dp), intent(in) :: radius
        type(interval_t) :: w
        w = iadd(value, interval(-radius, radius))
    end function
end module kin6d_debye_sphere_reference
