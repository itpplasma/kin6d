module kin6d_gs_coarea
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use kin6d_grad_shafranov, only: gs_state_t, &
                                    gs_element_jet, gs_axis_candidate, &
                                    gs_axis_safety_factor
    use kin6d_gs_coarea_leaves, only: gs_contour_weight, &
                                      gs_contour_weight_jvp, &
                                      gs_ray_trace, &
                                      gs_contour_root_response, &
            gs_moving_jet_response, &
                                      gs_edge_trace_coefficients, &
                                      gs_quadratic_translate, gs_polar_angle, &
                                      gs_polar_angle_jvp, &
                                      gs_boundary_weight, gs_boundary_weight_jvp
    use kin6d_grad_shafranov_leaves, &
        only: gs_smooth_axis_q_jvp
    use fortnum_linalg, only: dense_solve
    use fortnum_quadrature, only: gauss_legendre_ab
    use fortnum_multiroot, only: argsort
    use fortnum_roots, only: root_brent
    use fortnum_status, only: fortnum_status_t, FORTNUM_OK
    implicit none
    private
    public :: gs_level_coarea
    type :: event_t
        real(dp) :: angle = 0, tangent = 0
        integer :: left = 0, right = 0
    end type
contains
    subroutine gs_level_coarea(state, levels, &
                               angular_points, coarea, axis, status, &
                               psi_tangent, levels_tangent, coarea_tangent)
        type(gs_state_t), intent(in) :: state
        real(dp), intent(in) :: levels(:)
        integer, intent(in) :: angular_points
        real(dp), intent(out) :: coarea(:), axis(2)
        integer, intent(out) :: status
        real(dp), optional, intent(in) :: psi_tangent(:), &
                                          levels_tangent(:)
        real(dp), optional, intent(out) :: coarea_tangent(:)
        type(gs_state_t) :: tangent_state
        type(event_t), allocatable :: events(:), events_copy(:)
        integer, allocatable :: permutation(:)
        real(dp), allocatable :: x(:), w(:), jet(:, :), &
                                 hess(:, :, :)
        real(dp) :: axis_flux, axis_h(2, 2), axis_delta(2), &
                    axis_flux_delta, gh(2, 2), gg(2)
        real(dp) :: q_geometry, dq, value, grad(2), h(2, 2), &
                    center(2), shifted_grad(2)
        real(dp) :: point(2), vec(2), coeff_b, coeff_c, &
                    roots(2)
        real(dp) :: lo, hi, dlo, dhi, angle
        real(dp) :: t, dt, tv, tg(2), th(2, 2), two_pi, &
                    tolerance, level_delta
        integer :: e, k, j, m, stat, axis_element, n, count, &
                   cell, first_cell, next_cell
        integer :: order, nr, a, b, mid, index, nb
        logical :: with_tangent, is_boundary
        coarea = 0
        axis = 0
        status = 1
        with_tangent = present(psi_tangent) .and. &
                       present(levels_tangent) .and. &
                       present(coarea_tangent)
        if ((present(psi_tangent) .or. present(levels_tangent) .or. &
             present(coarea_tangent)) .and. &
            .not. with_tangent) return
        if (size(levels) < 2 .or. &
            size(coarea) /= size(levels)) return
        if (state%element_order /= 2) return
        if (.not. all(ieee_is_finite(levels))) return
        if (any(levels(2:) <= levels(:size(levels) - &
                                     1))) return
        if (angular_points < 16 .or. &
            angular_points > 16384) return
        if (with_tangent) then
            if (size(psi_tangent) /= size(state%psi)) return
            if (size(levels_tangent) /= size(levels) .or. &
                size(coarea_tangent) /= size(levels)) return
            if (.not. all(ieee_is_finite(psi_tangent))) return
            if (.not. all(ieee_is_finite(levels_tangent))) return
            tangent_state = state
            tangent_state%psi = psi_tangent
            coarea_tangent = 0
        end if
        call gs_axis_candidate(state, axis, axis_flux, &
                               axis_h, axis_element, status)
        if (status /= 0) return
        status = 1
        two_pi = 2*acos(-1.0_dp)
        tolerance = 1024*epsilon(1.0_dp)*max(1.0_dp, &
                                             maxval(abs(state%psi)))
        if (abs(levels(1) - axis_flux) > tolerance) then
            status = 13
            return
        end if
        call gs_axis_safety_factor(axis(1), 1.0_dp, axis_h, &
                                   q_geometry, stat)
        if (stat /= 0) then
            status = 14
            return
        end if
        coarea(1) = two_pi*q_geometry
        axis_delta = 0
        if (with_tangent) then
            call gs_element_jet(tangent_state, axis_element, &
                                axis(1), axis(2), &
                                axis_flux_delta, gg, gh, stat)
            if (stat /= 0) return
            if (abs(levels_tangent(1) - &
                    axis_flux_delta) > tolerance) then
                status = 13
                return
            end if
            call dense_solve(axis_h, -gg, axis_delta, stat)
            if (stat /= 0) return
            call gs_smooth_axis_q_jvp(axis(1), 1.0_dp, axis_h(1, &
                                                              1), axis_h(1, 2), &
            axis_h(2, 2), &
                                      axis_delta(1), 0.0_dp, gh(1, 1), gh(1, 2), &
            gh(2, 2), &
                                      dq)
            coarea_tangent(1) = two_pi*dq
        end if
        n = size(state%triangles, 2)
        allocate (jet(3, n), hess(2, 2, n), events(6*n), events_copy(6*n), permutation(6*n))
        ! Exact quadratic extension of each cell to the axis; roots remain
        ! qualified against the original cell, never against the extension.
        do e = 1, n
            center = sum(state%points(:, state%triangles(1:3, &
                                                         e)), dim=2)/3
            call gs_element_jet(state, e, center(1), center(2), &
                                value, grad, h, stat)
            if (stat /= 0) return
            call gs_quadratic_translate(value, grad(1), grad(2), &
                                        h(1, 1), h(1, 2), h(2, 2), &
                                        axis(1) - center(1), axis(2) - center(2), &
            jet(1, e), &
                                        shifted_grad(1), shifted_grad(2))
            jet(2:3, e) = shifted_grad
            hess(:, :, e) = h
        end do
        ! angular_points is a quadrature budget: integrate every smooth arc
        ! at order budget/16, rather than sample across jumps. Doubling
        ! the budget doubles the order; there is no hidden saturation.
        order = angular_points/16
        allocate (x(order), w(order))
        call gauss_legendre_ab(order, 0.0_dp, 1.0_dp, &
                               x, w)
        do j = 2, size(levels)
            level_delta = 0
            if (with_tangent) level_delta = levels_tangent(j)
            is_boundary = .true.
            nb = 0
            do e = 1, n
                do k = 1, 3
                    if (state%neighbors(k, e) > 0) cycle
                    nb = nb + 1
                    a = state%triangles(k, &
                                        e)
                    b = state%triangles(mod(k, 3) + 1, e)
                    mid = state%triangles(k + 3, e)
                    if (maxval(abs(state%psi([a, mid, b]) - &
                                   levels(j))) > tolerance) is_boundary = .false.
                end do
            end do
            if (nb == 0) is_boundary = .false.
            if (is_boundary) then
                call integrate_boundary(j, level_delta, stat)
                if (stat /= 0) then
                    status = stat
                    return
                end if
                cycle
            end if
            count = 0
            do e = 1, n
                do k = 1, 3
                    next_cell = state%neighbors(k, e)
                    if (next_cell > 0) then
                        if (next_cell < e) cycle
                    end if
                    a = state%triangles(k, e)
                    b = state%triangles(mod(k, 3) + 1, &
                                        e)
                    mid = state%triangles(k + 3, e)
                    call gs_edge_trace_coefficients(state%psi(a), &
                                                    state%psi(mid), &
                                                    state%psi(b), coeff_b, coeff_c)
                    call bracket_roots(coeff_c, coeff_b, &
                                       state%psi(a) - levels(j), 0.0_dp, 1.0_dp, &
            roots, nr)
                    do m = 1, nr
                        t = roots(m)
                        if (t < 0 .or. t > 1) cycle
                        ! A contour through a mesh vertex has no unique face
                        ! event derivative. Refuse that branch explicitly.
                        if (min(t, 1 - t) < 1.0e-10_dp) then
                            status = 18
                            return
                        end if
                        if (next_cell < 1) then
                            status = 7
                            return
                        end if
                        vec = state%points(:, b) - state%points(:, a)
                        point = state%points(:, a) + t*vec
                        call gs_polar_angle(point(1) - axis(1), point(2) - &
                                            axis(2), angle)
                        count = count + 1
                        events(count)%angle = modulo(angle, two_pi)
                        events(count)%left = e
                        events(count)%right = next_cell
                        events(count)%tangent = 0
                        if (with_tangent) then
                            call gs_element_jet(state, e, point(1), point(2), &
                                                value, grad, h, stat)
                            if (stat /= 0) return
                            call gs_element_jet(tangent_state, e, point(1), &
                                                point(2), tv, tg, th, stat)
                            if (stat /= 0) return
                            if (abs(dot_product(grad, vec)) < tolerance) then
                                status = 18
                                return
                            end if
                            call gs_contour_root_response(dot_product(grad, vec), &
                                                          level_delta, tv, dt)
                            call gs_polar_angle_jvp(point(1) - axis(1), point(2) - &
                                                    axis(2), &
                                                    vec(1)*dt - axis_delta(1), &
            vec(2)*dt - axis_delta(2), &
                                                    events(count)%tangent)
                        end if
                    end do
                end do
            end do
            if (count == 0) then
                ! A contour wholly inside the axis cell has no face breaks.
                ! Sixteen regular sectors keep full-circle Gauss quadrature
                ! independent of the size/anisotropy of the local ellipse.
                cell = axis_element
                do index = 1, 16
                    lo = (index - 1)*two_pi/16
                    hi = index*two_pi/16
                    call integrate_arc(cell, lo, hi, 0.0_dp, 0.0_dp, j, &
                                       level_delta, stat)
                    if (stat /= 0) then
                        status = stat
                        return
                    end if
                end do
                cycle
            end if
            call argsort(events(:count)%angle, permutation(:count))
            events_copy(:count) = events(:count)
            do index = 1, count
                events(index) = events_copy(permutation(index))
            end do
            if (count > 1) then
                if (minval(events(2:count)%angle - events(1:count - &
                                                          1)%angle) < 1.0e-10_dp) then
                    status = 18
                    return
                end if
            end if
            lo = events(1)%angle
            hi = lo + two_pi
            if (count > 1) hi = events(2)%angle
            cell = 0
            do e = 1, n
                call point_on_cell(e, 0.5_dp*(lo + hi), levels(j), &
                                   point, grad, h, stat)
                if (stat /= 0) cycle
                cell = e
                exit
            end do
            if (cell == 0) then
                status = 7
                return
            end if
            first_cell = cell
            do index = 1, count
                lo = events(index)%angle
                dlo = events(index)%tangent
                k = mod(index, count) + 1
                hi = events(k)%angle
                dhi = events(k)%tangent
                if (k == 1) hi = hi + two_pi
                call integrate_arc(cell, lo, hi, dlo, dhi, j, &
                                   level_delta, stat)
                if (stat /= 0) then
                    status = stat
                    return
                end if
                if (events(k)%left == cell) then
                    cell = events(k)%right
                else if (events(k)%right == cell) then
                    cell = events(k)%left
                else
                    status = 8
                    return
                end if
            end do
            if (cell /= first_cell) then
                status = 8
                return
            end if
        end do
        if (.not. all(ieee_is_finite(coarea)) .or. &
            minval(coarea) <= 0) then
            status = 15
            return
        end if
        if (with_tangent) then
            if (.not. all(ieee_is_finite(coarea_tangent))) then
                status = 17
                return
            end if
        end if
        status = 0
    contains
        subroutine point_on_cell(e, theta, level, point, &
                                 grad, h, status)
            integer, intent(in)::e
            real(dp), intent(in)::theta, level
            real(dp), intent(out)::point(2), grad(2), h(2, 2)
            integer, intent(out)::status
            real(dp)::dir(2), bary(3), slope(3), lo, hi, aa, bb, &
                       cc, root(2), r, value
            integer::k, nr, found, stat
            status = 7
            point = 0
            grad = 0
            h = 0
            dir = [cos(theta), sin(theta)]
            point = axis - state%points(:, state%triangles(1, e))
            bary(2:3) = matmul(transpose(state%gradient(:, 2:3, e)), point)
            bary(1) = 1 - sum(bary(2:3))
            slope = matmul(transpose(state%gradient(:, &
                                                    :, e)), dir)
            lo = 0
            hi = huge(1.0_dp)
            do k = 1, 3
                if (slope(k) > 0) then
                    lo = max(lo, -bary(k)/slope(k))
                else if (slope(k) < 0) then
                    hi = min(hi, -bary(k)/slope(k))
                else
                    if (bary(k) < -1.0e-12_dp) return
                end if
            end do
            if (hi <= lo .or. hi >= huge(1.0_dp)) return
            aa = 0.5_dp*dot_product(dir, matmul(hess(:, :, e), &
                                                dir))
            bb = dot_product(jet(2:3, e), dir)
            cc = jet(1, e) - level
            ! Monotone outward trace is part of this native coarea contract.
            if (min(bb + 2*aa*lo, bb + 2*aa*hi) < -tolerance) then
                status = 6
                return
            end if
            call bracket_roots(aa, bb, cc, lo, hi, root, nr)
            if (nr == 0) then
                if (abs(cc + lo*(bb + aa*lo)) <= tolerance) then
                    root(1) = lo
                    nr = 1
                else if (abs(cc + hi*(bb + aa*hi)) <= tolerance) then
                    root(1) = hi
                    nr = 1
                end if
            end if
            found = 0
            do k = 1, nr
                r = root(k)
                if (r <= 0) cycle
                if (bb + 2*aa*r <= 0) cycle
                if (found /= 0) return
                point = axis + r*dir
                found = 1
            end do
            if (found == 0) return
            call gs_element_jet(state, e, point(1), point(2), &
                                value, grad, h, stat)
            if (stat /= 0) return
            if (abs(value - level) > tolerance) then
                status = 12
                return
            end if
            status = 0
        end subroutine point_on_cell

        subroutine arc_weight(e, theta, level, level_delta, &
                              weight, delta_weight, status)
            integer, intent(in)::e
            real(dp), intent(in)::theta, level, level_delta
            real(dp), intent(out)::weight, delta_weight
            integer, intent(out)::status
            real(dp)::point(2), grad(2), h(2, 2), dir(2), r, rg, droot, &
                       tv, tg(2), th(2, 2), dv, dg(2), ignored
            weight = 0
            delta_weight = 0
            call point_on_cell(e, theta, level, point, grad, h, &
                               status)
            if (status /= 0) return
            dir = [cos(theta), &
                   sin(theta)]
            r = norm2(point - axis)
            rg = dot_product(grad, dir)
            call gs_contour_weight(r, point(1), rg, weight)
            if (.not. with_tangent) return
            call gs_element_jet(tangent_state, e, point(1), &
                                point(2), tv, tg, th, status)
            if (status /= 0) return
            call gs_moving_jet_response(grad(1), grad(2), h(1, &
                                                            1), h(1, 2), h(2, 2), &
                                        tv, tg(1), tg(2), axis_delta(1), &
            axis_delta(2), dv, &
                                        dg(1), dg(2))
            call gs_contour_root_response(rg, level_delta, dv, &
                                          droot)
            call gs_moving_jet_response(grad(1), grad(2), h(1, 1), &
                                        h(1, 2), h(2, 2), tv, tg(1), tg(2), &
                                        axis_delta(1) + droot*dir(1), axis_delta(2) + &
                                        droot*dir(2), ignored, dg(1), dg(2))
            call gs_contour_weight_jvp(r, point(1), rg, droot, &
                                       axis_delta(1) + droot*dir(1), &
                                       dot_product(dg, dir), delta_weight)
        end subroutine arc_weight

        subroutine integrate_arc(e, lo, hi, dlo, dhi, j, &
                                 level_delta, status)
            integer, intent(in)::e, j
            real(dp), intent(in)::lo, hi, dlo, dhi, level_delta
            integer, intent(out)::status
            real(dp)::weight, delta_weight, angle
            integer::q
            do q = 1, order
                angle = lo + (hi - lo)*x(q)
                call arc_weight(e, angle, levels(j), level_delta, &
                                weight, delta_weight, status)
                if (status /= 0) return
                coarea(j) = coarea(j) + (hi - lo)*w(q)*weight
                if (with_tangent) coarea_tangent(j) = coarea_tangent(j) + (hi - &
            lo)*w(q)*delta_weight
            end do
            if (with_tangent) then
                ! Leibniz endpoint terms supply moving gradient-jump products;
                ! adjacent arcs cancel them only for a C1 native polynomial.
                call arc_weight(e, lo, levels(j), level_delta, &
                                weight, delta_weight, status)
                if (status /= 0) return
                coarea_tangent(j) = coarea_tangent(j) - weight*dlo
                call arc_weight(e, hi, levels(j), level_delta, &
                                weight, delta_weight, status)
                if (status /= 0) return
                coarea_tangent(j) = coarea_tangent(j) + weight*dhi
            end if
            status = 0
        end subroutine integrate_arc

        subroutine integrate_boundary(j, level_delta, status)
            integer, intent(in)::j
            real(dp), intent(in)::level_delta
            integer, intent(out)::status
            real(dp)::point(2), vec(2), value, grad(2), h(2, 2), tv, &
                       tg(2), th(2, 2), weight, dweight
            integer::e, k, a, b, mid, q, stat
            status = 1
            do e = 1, n
                do k = 1, 3
                    if (state%neighbors(k, e) > 0) cycle
                    a = state%triangles(k, e)
                    b = state%triangles(mod(k, 3) + 1, &
                                        e)
                    mid = state%triangles(k + 3, e)
                    vec = state%points(:, b) - state%points(:, a)
                    if (with_tangent) then
                        if (maxval(abs(psi_tangent([a, mid, b]) - &
                                       level_delta)) > tolerance) then
                            status = 18
                            return
                        end if
                    end if
                    do q = 0, order + 1
                        if (q == 0) then
                            point = state%points(:, a)
                        else if (q == order + 1) then
                            point = state%points(:, b)
                        else
                            point = state%points(:, a) + x(q)*vec
                        end if
                        call gs_element_jet(state, e, point(1), point(2), &
                                            value, grad, h, stat)
                        if (stat /= 0) return
                        if (norm2(grad) <= tolerance) then
                            status = 19
                            return
                        end if
                        if (q == 0 .or. q == order + 1) cycle
                        call gs_boundary_weight(norm2(vec), point(1), &
                                                grad(1), grad(2), weight)
                        coarea(j) = coarea(j) + w(q)*weight
                        if (with_tangent) then
                            call gs_element_jet(tangent_state, e, point(1), &
                                                point(2), tv, tg, th, stat)
                            if (stat /= 0) return
                            call gs_boundary_weight_jvp(norm2(vec), point(1), &
                                                        grad(1), grad(2), &
                                                        0.0_dp, 0.0_dp, tg(1), &
            tg(2), dweight)
                            coarea_tangent(j) = coarea_tangent(j) + w(q)*dweight
                        end if
                    end do
                end do
            end do
            status = 0
        end subroutine integrate_boundary
    end subroutine gs_level_coarea

    subroutine bracket_roots(a, b, c, lo, hi, roots, n)
        real(dp), intent(in)::a, b, c, lo, hi
        real(dp), intent(out)::roots(2)
        integer, intent(out)::n
        real(dp)::knots(3), fl, fr, ftol, xtol, x, critical
        integer::nk, k
        type(fortnum_status_t)::stat
        roots = 0
        n = 0
        knots(1) = lo
        knots(2) = hi
        nk = 2
        if (abs(a) > tiny(1.0_dp)) then
            critical = -b/(2*a)
            if (critical > lo .and. critical < hi) then
                nk = 3
                knots(2) = critical
                knots(3) = hi
            end if
        end if
        ftol = 16*epsilon(1.0_dp)*max(abs(a)*hi**2, &
                                      abs(b)*hi, abs(c), tiny(1.0_dp))
        xtol = 16*epsilon(1.0_dp)*max(1.0_dp, abs(lo), abs(hi))
        do k = 1, nk - 1
            fl = polynomial(knots(k))
            fr = polynomial(knots(k + 1))
            if (abs(fl) <= ftol) then
                x = knots(k)
            else if (abs(fr) <= ftol) then
                x = knots(k + 1)
            else
                if ((fl > 0 .and. fr > 0) .or. (fl < 0 .and. &
                                                fr < 0)) cycle
                call root_brent(polynomial, knots(k), knots(k + 1), &
                                x, stat, xtol, ftol, 64)
                if (stat%code /= FORTNUM_OK) cycle
            end if
            if (n > 0) then
                if (abs(x - roots(n)) <= xtol) cycle
            end if
            n = n + 1
            roots(n) = x
        end do
    contains
        pure function polynomial(x) result(y)
            real(dp), intent(in)::x
            real(dp)::y, derivative
            call gs_ray_trace(c, b, 2*a, x, y, derivative)
        end function polynomial
    end subroutine bracket_roots
end module kin6d_gs_coarea
