module kin6d_solovev_geometry
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortnum_roots, only: root_bisect
    use fortnum_linalg, only: det2
    use fortnum_status, only: fortnum_status_t, status_ok
    use kin6d_solovev_geometry_leaves, only: solovev_geometry_level, &
        solovev_geometry_rectangle, solovev_boundary_jet, solovev_inverse_scale, &
        solovev_inverse_discriminant, solovev_forward_scale, &
        solovev_h_first, solovev_h_second, &
        solovev_facet_first, solovev_facet_second, solovev_c_first, &
        solovev_c_second, solovev_map_first, solovev_map_second
    implicit none
    private
    integer, parameter, public :: geometry_ok = 0, geometry_invalid = 1
    integer, parameter, public :: geometry_interface = 2, geometry_arithmetic = 3
    type, public :: solovev_geometry_t
        private
        real(dp) :: r0 = 0, level = 0, tau = 0, rectangle(3) = 0
        real(dp), allocatable :: normals(:, :), beta(:)
        logical :: ready = .false.
    end type
    type, public :: geometry_sample_t
        ! Tensor order: physical component, first coordinate, second coordinate.
        real(dp) :: value(2) = 0, jacobian(2, 2) = 0, hessian(2, 2, 2) = 0
        real(dp) :: root_residual = 0
        integer :: facet = 0, status = geometry_invalid
        logical :: identity = .false., inner_interface = .false.
        logical :: sector_interface = .false., certified = .false.
    end type
    type, public :: geometry_position_t
        real(dp) :: value(2) = 0, root_residual = 0
        integer :: facet = 0, status = geometry_invalid
        logical :: identity = .false., inner_interface = .false.
        logical :: sector_interface = .false., certified = .false.
    end type
    public :: initialize_solovev_geometry, sample_solovev_geometry
    public :: position_solovev_geometry, locate_solovev_geometry
    public :: solovev_geometry_jvp, solovev_geometry_vjp
contains
    subroutine initialize_solovev_geometry(map, r0, level, tau, polygon, status)
        type(solovev_geometry_t), intent(out) :: map
        real(dp), intent(in) :: r0, level, tau, polygon(:, :)
        integer, intent(out) :: status
        real(dp) :: edge(2), normal(2), offset, value, tolerance
        integer :: i, j, count
        status = geometry_invalid
        if (.not. all(ieee_is_finite([r0, level, tau]))) return
        if (r0 <= 0 .or. level <= 0 .or. tau <= 0 .or. tau >= 1) return
        if (sqrt(level) >= r0**2) return
        if (size(polygon, 1) /= 2) return
        count = size(polygon, 2)
        if (count < 3) return
        if (.not. all(ieee_is_finite(polygon))) return
        if (any(polygon(1, :) <= 0)) return
        allocate (map%normals(2, count), map%beta(count))
        do i = 1, count
            edge = polygon(:, mod(i, count)+1)-polygon(:, i)
            normal = [edge(2), -edge(1)]
            if (norm2(normal) <= 0) return
            offset = dot_product(normal, polygon(:, i)-[r0, 0.0_dp])
            if (offset <= 0) return
            tolerance = 64*epsilon(1.0_dp)*max(tiny(1.0_dp), abs(offset))
            do j = 1, count
                value = dot_product(normal, polygon(:, j)-[r0, 0.0_dp])
                if (value > offset+tolerance) return
            end do
            call solovev_geometry_level(polygon(1, i), polygon(2, i), r0, value)
            if (.not. ieee_is_finite(value)) return
            if (value > level+64*epsilon(1.0_dp)*max(tiny(1.0_dp), level)) return
            map%normals(:, i) = normal
            map%beta(i) = offset
        end do
        map%r0 = r0; map%level = level; map%tau = tau
        call solovev_geometry_rectangle(r0, level, map%rectangle(1), &
            map%rectangle(2), map%rectangle(3))
        if (.not. all(ieee_is_finite(map%rectangle))) return
        ! Finite checks above do not certify convex-Omega, q*Omega inclusion,
        ! positive radial derivative or floating root/Jacobian enclosures.
        map%ready = .true.
        status = geometry_ok
    end subroutine

    subroutine sample_solovev_geometry(map, point, sample)
        type(solovev_geometry_t), intent(in) :: map
        real(dp), intent(in) :: point(2)
        type(geometry_sample_t), intent(out) :: sample
        call evaluate_solovev_geometry(map, point, sample, .false.)
    end subroutine

    subroutine position_solovev_geometry(map, point, position)
        type(solovev_geometry_t), intent(in) :: map
        real(dp), intent(in) :: point(2)
        type(geometry_position_t), intent(out) :: position
        type(geometry_sample_t) :: sample
        call evaluate_solovev_geometry(map, point, sample, .true.)
        position = geometry_position_t()
        position%status = sample%status
        position%facet = sample%facet
        position%identity = sample%identity
        position%inner_interface = sample%inner_interface
        position%sector_interface = sample%sector_interface
        if (sample%status /= geometry_ok) return
        position%value = sample%value
        position%root_residual = sample%root_residual
    end subroutine

    subroutine evaluate_solovev_geometry(map, point, sample, position_only)
        type(solovev_geometry_t), intent(in) :: map
        real(dp), intent(in) :: point(2)
        type(geometry_sample_t), intent(out) :: sample
        logical, intent(in) :: position_only
        real(dp) :: d(2), inner(2), q, value, score, best, runner, tolerance
        real(dp) :: h, p, c, discriminant, v, boundary(2), g(2)
        real(dp) :: k(2, 2), hi(2), hij(2, 2), pi(2), pij(2, 2), ci(2), cij
        real(dp) :: bj(2, 2), normal(2), beta, ei(2), ej(2), determinant
        real(dp) :: mapped(2), jacobian(2, 2), hessian(2, 2, 2), residual
        integer :: i, j, a, facet, ray_status
        sample = geometry_sample_t()
        if (.not. map%ready) return
        if (.not. all(ieee_is_finite(point))) return
        d = point-[map%r0, 0.0_dp]
        best = -huge(1.0_dp); runner = -huge(1.0_dp); facet = 0
        do i = 1, size(map%beta)
            score = dot_product(map%normals(:, i), d)/map%beta(i)
            if (.not. ieee_is_finite(score)) return
            if (score > best) then
                runner = best; best = score; facet = i
            else
                runner = max(runner, score)
            end if
        end do
        tolerance = 128*epsilon(1.0_dp)*max(1.0_dp, abs(best))
        if (best > 1+tolerance) return
        q = 1-map%tau
        inner = [map%r0, 0.0_dp]+d/q
        if (.not. all(ieee_is_finite(inner))) return
        if (inner(1) > 0) then
            call solovev_geometry_level(inner(1), inner(2), map%r0, value)
            if (.not. ieee_is_finite(value)) return
            tolerance = 128*epsilon(1.0_dp)*max(tiny(1.0_dp), map%level, abs(value))
            if (abs(value-map%level) <= tolerance) then
                sample%inner_interface = .true.
                if (.not. position_only) then
                    sample%status = geometry_interface
                    return
                end if
            end if
            if (value < map%level) then
                sample%value = point
                sample%jacobian(1, 1) = 1; sample%jacobian(2, 2) = 1
                sample%identity = .true.
                sample%status = geometry_ok
                return
            end if
        end if
        sample%facet = facet
        tolerance = 128*epsilon(1.0_dp)*max(1.0_dp, abs(best))
        if (best-runner <= tolerance) then
            sample%sector_interface = .true.
            if (.not. position_only) then
                sample%status = geometry_interface
                return
            end if
        end if
        if (best <= 0) return
        normal = map%normals(:, facet); beta = map%beta(facet)
        p = 1/best
        sample%status = geometry_arithmetic
        call solovev_ray_scale(map, d, h, residual, ray_status)
        if (ray_status /= geometry_ok) return
        if (h < p) return
        if (.not. position_only) then
            if (q*h >= 1) return
        end if
        call solovev_inverse_discriminant(map%tau, h, p, discriminant)
        if (.not. ieee_is_finite(discriminant)) return
        if (discriminant <= 0) return
        call solovev_inverse_scale(map%tau, h, p, c)
        if (.not. ieee_is_finite(c)) return
        ! Range guards allow local floating evaluation noise without clipping
        ! the root or claiming an error enclosure. Exact predicates stay external.
        tolerance = 64*epsilon(1.0_dp)*max(1.0_dp, h, abs(c))
        if (c < 1-tolerance .or. c > h+tolerance) return
        if (position_only) then
            mapped(1) = map%r0+c*d(1)
            mapped(2) = 0.0_dp+c*d(2)
            if (.not. all(ieee_is_finite(mapped))) return
            sample%value = mapped
            sample%root_residual = residual
            sample%status = geometry_ok
            return
        end if
        boundary = [map%r0, 0.0_dp]+h*d
        call solovev_boundary_jet(boundary(1), boundary(2), map%r0, &
            g(1), g(2), k(1, 1), k(1, 2), k(2, 2))
        k(2, 1) = k(1, 2)
        v = dot_product(g, d)
        if (.not. ieee_is_finite(v)) return
        if (v <= 0) return
        do i = 1, 2
            call solovev_h_first(h, g(i), v, hi(i))
            call solovev_facet_first(p, beta, normal(i), pi(i))
            bj(:, i) = d*hi(i)
            bj(i, i) = bj(i, i)+h
            call solovev_c_first(map%tau, h, p, c, hi(i), pi(i), ci(i))
        end do
        do j = 1, 2
            do i = 1, 2
                call solovev_h_second(g(i), g(j), v, bj(1, i), bj(2, i), &
                    bj(1, j), bj(2, j), k(1, 1), k(1, 2), k(2, 2), &
                    hi(i), hi(j), hij(i, j))
                call solovev_facet_second(p, beta, normal(i), normal(j), pij(i, j))
                call solovev_c_second(map%tau, h, p, c, hi(i), hi(j), hij(i, j), &
                    pi(i), pi(j), pij(i, j), ci(i), ci(j), cij)
                do a = 1, 2
                    ei = 0; ej = 0; ei(i) = 1; ej(j) = 1
                    call solovev_map_second(d(a), ci(i), ci(j), cij, ei(a), &
                        ej(a), hessian(a, i, j))
                end do
            end do
        end do
        do i = 1, 2
            do a = 1, 2
                ei = 0; ei(i) = 1
                call solovev_map_first(c, d(a), ci(i), ei(a), jacobian(a, i))
            end do
        end do
        mapped = [map%r0, 0.0_dp]+c*d
        if (.not. ieee_is_finite(residual)) return
        if (.not. all(ieee_is_finite(mapped))) return
        if (.not. all(ieee_is_finite(jacobian))) return
        if (.not. all(ieee_is_finite(hessian))) return
        determinant = det2(jacobian)
        if (.not. ieee_is_finite(determinant)) return
        if (determinant <= 0) return
        sample%value = mapped
        sample%jacobian = jacobian
        sample%hessian = hessian
        sample%root_residual = residual
        sample%status = geometry_ok
    end subroutine

    subroutine solovev_ray_scale(map, d, h, residual, status)
        type(solovev_geometry_t), intent(in) :: map
        real(dp), intent(in) :: d(2)
        real(dp), intent(out) :: h, residual
        integer, intent(out) :: status
        type(fortnum_status_t) :: root_status
        real(dp) :: upper, boundary(2), candidate
        h = 0; residual = 0; status = geometry_arithmetic
        upper = huge(1.0_dp)
        if (d(1) > 0) upper = (map%rectangle(2)-map%r0)/d(1)
        if (d(1) < 0) upper = (map%rectangle(1)-map%r0)/d(1)
        if (abs(d(2)) > 0) upper = min(upper, map%rectangle(3)/abs(d(2)))
        if (.not. ieee_is_finite(upper)) return
        if (upper <= 0) return
        ! Floating bracket recovery, not a certified root enclosure.
        upper = upper*(1+64*epsilon(1.0_dp))
        boundary = [map%r0, 0.0_dp]+upper*d
        if (.not. all(ieee_is_finite(boundary))) return
        if (boundary(1) <= 0) return
        call root_bisect(ray_residual, 0.0_dp, upper, candidate, root_status, &
            xtol=8*epsilon(1.0_dp)*max(1.0_dp, upper), ftol=0.0_dp)
        if (.not. status_ok(root_status)) return
        if (.not. ieee_is_finite(candidate)) return
        if (candidate <= 0) return
        residual = ray_residual(candidate)
        if (.not. ieee_is_finite(residual)) then
            residual = 0
            return
        end if
        h = candidate; status = geometry_ok
    contains
        pure real(dp) function ray_residual(scale) result(value)
            real(dp), intent(in) :: scale
            real(dp) :: point(2)
            point = [map%r0, 0.0_dp]+scale*d
            call solovev_geometry_level(point(1), point(2), map%r0, value)
            value = value-map%level
        end function
    end subroutine

    subroutine locate_solovev_geometry(map, physical, position)
        type(solovev_geometry_t), intent(in) :: map
        real(dp), intent(in) :: physical(2)
        type(geometry_position_t), intent(out) :: position
        real(dp) :: d(2), inner(2), value, tolerance, q, best, runner, score
        real(dp) :: h, p, scale, residual, point(2)
        integer :: i, facet, ray_status
        position = geometry_position_t()
        if (.not. map%ready) return
        if (.not. all(ieee_is_finite(physical))) return
        if (physical(1) <= 0) return
        call solovev_geometry_level(physical(1), physical(2), map%r0, value)
        if (.not. ieee_is_finite(value)) return
        tolerance = 128*epsilon(1.0_dp)*max(tiny(1.0_dp), map%level, abs(value))
        if (value > map%level+tolerance) return
        d(1) = physical(1)-map%r0
        d(2) = physical(2)-0.0_dp
        q = 1-map%tau
        inner(1) = map%r0+d(1)/q
        inner(2) = 0.0_dp+d(2)/q
        if (.not. all(ieee_is_finite(inner))) return
        if (inner(1) > 0) then
            call solovev_geometry_level(inner(1), inner(2), map%r0, value)
            if (.not. ieee_is_finite(value)) return
            tolerance = 128*epsilon(1.0_dp)*max(tiny(1.0_dp), map%level, abs(value))
            position%inner_interface = abs(value-map%level) <= tolerance
            if (value < map%level) then
                position%value = physical
                position%identity = .true.
                position%status = geometry_ok
                return
            end if
        end if
        best = -huge(1.0_dp); runner = best; facet = 0
        do i = 1, size(map%beta)
            score = dot_product(map%normals(:, i), d)/map%beta(i)
            if (.not. ieee_is_finite(score)) return
            if (score > best) then
                runner = best; best = score; facet = i
            else
                runner = max(runner, score)
            end if
        end do
        if (best <= 0) return
        position%facet = facet
        tolerance = 128*epsilon(1.0_dp)*max(1.0_dp, abs(best))
        position%sector_interface = best-runner <= tolerance
        position%status = geometry_arithmetic
        call solovev_ray_scale(map, d, h, residual, ray_status)
        if (ray_status /= geometry_ok) return
        p = 1/best
        call solovev_forward_scale(map%tau, h, p, scale)
        if (.not. ieee_is_finite(scale)) return
        if (scale <= 0) return
        point(1) = map%r0+scale*d(1)
        point(2) = 0.0_dp+scale*d(2)
        if (.not. all(ieee_is_finite(point))) return
        position%value = point
        position%root_residual = residual
        position%status = geometry_ok
    end subroutine

    subroutine solovev_geometry_jvp(map, point, direction, product, status)
        type(solovev_geometry_t), intent(in) :: map
        real(dp), intent(in) :: point(2), direction(2)
        real(dp), intent(out) :: product(2)
        integer, intent(out) :: status
        type(geometry_sample_t) :: sample
        product = 0; status = geometry_invalid
        if (.not. all(ieee_is_finite(direction))) return
        call sample_solovev_geometry(map, point, sample)
        status = sample%status
        if (status /= geometry_ok) return
        product = matmul(sample%jacobian, direction)
        if (.not. all(ieee_is_finite(product))) then
            product = 0
            status = geometry_arithmetic
        end if
    end subroutine

    subroutine solovev_geometry_vjp(map, point, weight, product, status)
        type(solovev_geometry_t), intent(in) :: map
        real(dp), intent(in) :: point(2), weight(2)
        real(dp), intent(out) :: product(2)
        integer, intent(out) :: status
        type(geometry_sample_t) :: sample
        product = 0; status = geometry_invalid
        if (.not. all(ieee_is_finite(weight))) return
        call sample_solovev_geometry(map, point, sample)
        status = sample%status
        if (status /= geometry_ok) return
        product = matmul(transpose(sample%jacobian), weight)
        if (.not. all(ieee_is_finite(product))) then
            product = 0
            status = geometry_arithmetic
        end if
    end subroutine

end module kin6d_solovev_geometry
