module kin6d_gs_coarea
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use kin6d_grad_shafranov, only: gs_state_t, &
                                    gs_element_jet, gs_axis_candidate, &
                                    gs_axis_safety_factor, GS_UNSUPPORTED_GEOMETRY
    use kin6d_curved_triangle, only: triangle_geometry, triangle_locate
    use kin6d_triangle_basis, only: cubic_segment_roots, cubic_segment_stationary, &
        triangle_basis_jet
    use kin6d_gs_coarea_leaves, only: gs_contour_weight, gs_cubic_flux, &
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
    use fortnum_polynomial, only: lagrange_weights
    use fortnum_linalg, only: dense_solve
    use fortnum_quadrature, only: gauss_legendre_ab
    use fortnum_multiroot, only: argsort
    use fortnum_roots, only: root_brent
    use fortnum_status, only: fortnum_status_t, FORTNUM_OK
    implicit none
    private
    public :: gs_level_coarea
    type :: event_t
        ! No default initialisation: the per-thread buffers are large and only
        ! a small prefix is ever written.
        real(dp) :: angle, tangent, point(2)
        integer :: left, right
    end type
contains
    subroutine gs_level_coarea(state, levels, &
                               angular_points, coarea, axis, status, &
                               psi_tangent, levels_tangent, coarea_tangent, work, &
                               reference_cubic)
        type(gs_state_t), intent(in) :: state
        real(dp), intent(in) :: levels(:)
        integer, intent(in) :: angular_points
        real(dp), intent(out) :: coarea(:), axis(2)
        integer, intent(out) :: status
        real(dp), optional, intent(in) :: psi_tangent(:), &
                                          levels_tangent(:)
        real(dp), optional, intent(out) :: coarea_tangent(:)
        integer, optional, intent(out) :: work(:)
        logical, optional, intent(in) :: reference_cubic
        type(gs_state_t) :: tangent_state
        integer, allocatable :: edge_cells(:, :), edge_nodes(:, :)
        real(dp), allocatable :: x(:), w(:), jet(:, :), &
                                 hess(:, :, :), edge_range(:, :), cubic(:, :)
        real(dp) :: physical_axis(2), local_axis(2)
        real(dp) :: edge_values(4), critical(2), weights(4), extrema(4), boundary_level
        integer :: ne, nc, first_failure
        integer, allocatable :: level_status(:)
        logical :: constant_boundary, use_reference
        real(dp) :: transform(10, 10), shape(10), basis_g(2, 10), basis_h(2, 2, 10)
        real(dp) :: axis_flux, axis_h(2, 2), axis_delta(2), &
                    axis_flux_delta, gh(2, 2), gg(2)
        real(dp) :: q_geometry, dq, value, grad(2), h(2, 2), &
                    center(2), shifted_grad(2)
        real(dp) :: coeff_b, coeff_c
        real(dp) :: t, two_pi, tolerance
        integer :: e, k, m, stat, axis_element, n, next_cell
        integer :: order, nr, nb, edge_indices(4)
        logical :: with_tangent
        if (present(work)) then
            if (size(work) /= size(levels)) then
                status = 1
                return
            end if
            work = 0
        end if
        coarea = 0
        axis = 0
        status = 1
        if (state%mapped .or. state%element_order == 3) then
            if (present(psi_tangent) .or. present(levels_tangent) .or. &
                present(coarea_tangent)) then
                if (present(coarea_tangent)) coarea_tangent = 0
                status = GS_UNSUPPORTED_GEOMETRY
                return
            end if
        end if
        use_reference = .false.
        if (present(reference_cubic)) use_reference = reference_cubic
        with_tangent = present(psi_tangent) .and. &
                       present(levels_tangent) .and. &
                       present(coarea_tangent)
        if ((present(psi_tangent) .or. present(levels_tangent) .or. &
             present(coarea_tangent)) .and. &
            .not. with_tangent) return
        if (size(levels) < 2 .or. &
            size(coarea) /= size(levels)) return
        if (state%element_order < 2 .or. state%element_order > 3) return
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
        physical_axis = axis
        if (state%mapped) then
            if (state%curved_cell(axis_element)) then
                call triangle_locate( &
                    state%points(:, state%triangles(:, axis_element)), &
                    axis, local_axis, status)
                if (status /= 0) return
                axis = matmul(state%points(:, state%triangles(1:3, axis_element)), &
                    [1 - sum(local_axis), local_axis])
            end if
        end if
        status = 1
        two_pi = 2*acos(-1.0_dp)
        tolerance = 1024*epsilon(1.0_dp)*max(1.0_dp, &
                                             maxval(abs(state%psi)))
        if (abs(levels(1) - axis_flux) > tolerance) then
            status = 13
            return
        end if
        call gs_axis_safety_factor(physical_axis(1), 1.0_dp, axis_h, &
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
                                axis_flux_delta, gg, gh, stat, affine_chart=.true.)
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
        if (state%element_order == 3) then
            ! Contract the authoritative FortFEM basis once, then reuse its
            ! cubic coefficients for every ray. The reference schedule below
            ! contracts that same basis at each evaluation instead.
            allocate (cubic(10, n))
            call triangle_basis_jet([0.0_dp, 0.0_dp], shape, basis_g, basis_h)
            transform(1, :) = shape
            transform(2:3, :) = basis_g
            transform(4, :) = basis_h(1, 1, :)/2
            transform(5, :) = basis_h(1, 2, :)
            transform(6, :) = basis_h(2, 2, :)/2
            call triangle_basis_jet([1.0_dp, 0.0_dp], shape, basis_g, basis_h)
            transform(7, :) = (basis_h(1, 1, :) - 2*transform(4, :))/6
            transform(8, :) = (basis_h(1, 2, :) - transform(5, :))/2
            call triangle_basis_jet([0.0_dp, 1.0_dp], shape, basis_g, basis_h)
            transform(9, :) = (basis_h(1, 2, :) - transform(5, :))/2
            transform(10, :) = (basis_h(2, 2, :) - 2*transform(6, :))/6
            do e = 1, n
                cubic(:, e) = matmul(transform, state%psi(state%triangles(:, e)) &
                    - state%psi(state%triangles(1, e)))
            end do
        else
            allocate (jet(3, n), hess(2, 2, n))
        end if
        ! Exact quadratic extension of each cell to the axis; roots remain
        ! qualified against the original cell, never against the extension.
        do e = 1, n
            if (state%element_order == 3) exit
            center = sum(state%points(:, state%triangles(1:3, &
                                                         e)), dim=2)/3
            call gs_element_jet(state, e, center(1), center(2), &
                                value, grad, h, stat, affine_chart=.true.)
            if (stat /= 0) return
            call gs_quadratic_translate(value, grad(1), grad(2), &
                                        h(1, 1), h(1, 2), h(2, 2), &
                                        axis(1) - center(1), axis(2) - center(2), &
            jet(1, e), &
                                        shifted_grad(1), shifted_grad(2))
            jet(2:3, e) = shifted_grad
            hess(:, :, e) = h
        end do
        ! Precompute edge ranges once. Reject nonintersecting edges before
        ! solving their polynomial; nodal min/max alone is not a valid range.
        allocate (edge_cells(2, 3*n), edge_nodes(4, 3*n), edge_range(2, 3*n))
        ne = 0; nb = 0; constant_boundary = .true.; boundary_level = 0
        do e = 1, n
            do k = 1, 3
                next_cell = state%neighbors(k, e)
                if (next_cell > 0) then
                    if (next_cell < e) cycle
                end if
                call edge_dofs(e, k, edge_indices)
                ne = ne + 1
                edge_nodes(:, ne) = edge_indices
                edge_cells(:, ne) = [e, next_cell]
                edge_values = 0
                edge_values(:state%element_order + 1) = &
                    state%psi(edge_indices(:state%element_order + 1))
                if (next_cell < 1) then
                    if (nb == 0) boundary_level = edge_values(1)
                    nb = nb + 1
                    if (maxval(abs(edge_values(:state%element_order + 1) &
                        - boundary_level)) > tolerance) constant_boundary = .false.
                end if
                if (state%element_order == 2) then
                    call gs_edge_trace_coefficients(edge_values(1), &
                        edge_values(2), edge_values(3), coeff_b, coeff_c)
                    extrema(1:2) = edge_values([1, 3]); nc = 2
                    if (abs(coeff_c) > tiny(1.0_dp)) then
                        t = -coeff_b/(2*coeff_c)
                        if (t > 0 .and. t < 1) then
                            nc = 3
                            extrema(3) = edge_values(1) + t*(coeff_b + t*coeff_c)
                        end if
                    end if
                else
                    call cubic_segment_stationary(edge_values, critical, nr)
                    extrema(1:2) = edge_values([1, 4]); nc = nr + 2
                    do m = 1, nr
                        call lagrange_weights(4, critical(m), &
                            [0.0_dp, 1.0_dp/3, 2.0_dp/3, 1.0_dp], weights)
                        extrema(m + 2) = dot_product(weights, edge_values)
                    end do
                end if
                edge_range(:, ne) = [minval(extrema(:nc)), maxval(extrema(:nc))]
            end do
        end do
        if (nb == 0) constant_boundary = .false.
        ! angular_points is a quadrature budget: integrate every smooth arc
        ! at order budget/16, rather than sample across jumps. Doubling
        ! the budget doubles the order; there is no hidden saturation.
        order = angular_points/16
        allocate (x(order), w(order))
        call gauss_legendre_ab(order, 0.0_dp, 1.0_dp, &
                               x, w)
        allocate (level_status(size(levels)))
        level_status = 0
        first_failure = huge(1)
        !$omp parallel default(shared)
        block
            type(event_t), allocatable :: level_events(:), level_events_copy(:)
            integer, allocatable :: level_permutation(:)
            integer :: jl, snapshot
            allocate (level_events(6*n), level_events_copy(6*n), &
                      level_permutation(6*n))
            !$omp do schedule(dynamic, 1)
            do jl = 2, size(levels)
                !$omp atomic read
                snapshot = first_failure
                if (jl > snapshot) cycle
                call integrate_level(jl, level_events, level_events_copy, &
                                     level_permutation, level_status(jl))
                if (level_status(jl) /= 0) then
                    !$omp critical(gs_coarea_failure)
                    first_failure = min(first_failure, jl)
                    !$omp end critical(gs_coarea_failure)
                end if
            end do
            !$omp end do
        end block
        !$omp end parallel
        if (first_failure /= huge(1)) then
            status = level_status(first_failure)
            return
        end if
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
        axis = physical_axis
        status = 0
    contains
        subroutine integrate_level(j, events, events_copy, permutation, level_status)
            ! One contour level. Concurrent calls for distinct j are independent:
            ! the arrays are per-thread buffers and only coarea(j) is written.
            integer, intent(in) :: j
            type(event_t), intent(inout) :: events(:), events_copy(:)
            integer, intent(inout) :: permutation(:)
            integer, intent(out) :: level_status
            real(dp) :: level_delta, t, dt, tv, tg(2), th(2, 2), value, grad(2)
            real(dp) :: h(2, 2), point(2), vec(2), coeff_b, coeff_c, roots(3)
            real(dp) :: lo, hi, dlo, dhi, angle
            integer :: contour_work, count, edge, e, next_cell, edge_indices(4)
            integer :: a, b, mid, nr, m, stat, cell, first_cell, index, k
            logical :: is_boundary
            level_status = 1
            level_delta = 0
            if (with_tangent) level_delta = levels_tangent(j)
            contour_work = 0
            is_boundary = constant_boundary .and. &
                abs(levels(j) - boundary_level) <= tolerance
            if (is_boundary) then
                call integrate_boundary(j, level_delta, stat)
                if (stat /= 0) then
                    level_status = stat
                    return
                end if
                if (present(work)) work(j) = contour_work
                level_status = 0
                return
            end if
            count = 0
            do edge = 1, ne
                if (levels(j) < edge_range(1, edge) - tolerance .or. &
                    levels(j) > edge_range(2, edge) + tolerance) cycle
                contour_work = contour_work + 1
                e = edge_cells(1, edge); next_cell = edge_cells(2, edge)
                edge_indices = edge_nodes(:, edge)
                a = edge_indices(1); b = edge_indices(state%element_order + 1)
                if (state%element_order == 2) then
                    mid = edge_indices(2)
                    call gs_edge_trace_coefficients(state%psi(a), &
                        state%psi(mid), state%psi(b), coeff_b, coeff_c)
                    call bracket_roots(coeff_c, coeff_b, &
                        state%psi(a) - levels(j), 0.0_dp, 1.0_dp, roots(:2), nr)
                else
                    call cubic_segment_roots(state%psi(edge_indices), &
                        levels(j), roots, nr)
                end if
                do m = 1, nr
                    t = roots(m)
                    if (t < 0 .or. t > 1) cycle
                    ! A contour through a mesh vertex has no unique face
                    ! event derivative. Refuse that branch explicitly.
                    if (min(t, 1 - t) < 1.0e-10_dp) then
                        level_status = 18
                        return
                    end if
                    if (next_cell < 1) then
                        level_status = 7
                        return
                    end if
                    vec = state%points(:, b) - state%points(:, a)
                    point = state%points(:, a) + t*vec
                    call gs_polar_angle(point(1) - axis(1), point(2) - &
                                        axis(2), angle)
                    count = count + 1
                    events(count)%angle = modulo(angle, two_pi)
                    events(count)%point = point
                    events(count)%left = e
                    events(count)%right = next_cell
                    events(count)%tangent = 0
                    if (with_tangent) then
                        call gs_element_jet(state, e, point(1), point(2), &
                                            value, grad, h, stat, &
                                            affine_chart=.true.)
                        if (stat /= 0) return
                        call gs_element_jet(tangent_state, e, point(1), &
                                            point(2), tv, tg, th, stat, &
                                            affine_chart=.true.)
                        if (stat /= 0) return
                        if (abs(dot_product(grad, vec)) < tolerance) then
                            level_status = 18
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
            if (count == 0) then
                ! A contour wholly inside the axis cell has no face breaks.
                ! Sixteen regular sectors keep full-circle Gauss quadrature
                ! independent of the size/anisotropy of the local ellipse.
                cell = axis_element
                do index = 1, 16
                    lo = (index - 1)*two_pi/16
                    hi = index*two_pi/16
                    call integrate_arc(cell, lo, hi, 0.0_dp, 0.0_dp, j, &
                                       level_delta, stat, contour_work)
                    if (stat /= 0) then
                        level_status = stat
                        return
                    end if
                end do
                if (present(work)) work(j) = contour_work
                level_status = 0
                return
            end if
            call argsort(events(:count)%angle, permutation(:count))
            events_copy(:count) = events(:count)
            do index = 1, count
                events(index) = events_copy(permutation(index))
            end do
            ! Nearly equal event angles are a regular contour passing near a
            ! mesh vertex: the face events stay distinct and ordered, and the
            ! cell walk below rejects any inconsistent order (status 8).
            lo = events(1)%angle
            hi = lo + two_pi
            if (count > 1) hi = events(2)%angle
            cell = 0
            ! The first arc belongs to one of the first event's two cells.
            do k = 1, 2
                e = events(1)%left
                if (k == 2) e = events(1)%right
                call point_on_cell(e, 0.5_dp*(lo + hi), levels(j), &
                                   point, grad, h, stat, contour_work)
                if (stat /= 0) cycle
                cell = e
                exit
            end do
            if (cell == 0) then
                level_status = 7
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
                                   level_delta, stat, contour_work, &
                                   events(index)%point, events(k)%point)
                if (stat /= 0) then
                    level_status = stat
                    return
                end if
                if (events(k)%left == cell) then
                    cell = events(k)%right
                else if (events(k)%right == cell) then
                    cell = events(k)%left
                else
                    level_status = 8
                    return
                end if
            end do
            if (present(work)) work(j) = contour_work
            if (cell /= first_cell) then
                level_status = 8
                return
            end if
            level_status = 0
        end subroutine integrate_level

        subroutine edge_dofs(e, edge, indices)
            integer, intent(in) :: e, edge
            integer, intent(out) :: indices(4)
            integer :: k, degree
            degree = state%element_order
            indices = 0
            indices(1) = state%triangles(edge, e)
            indices(degree + 1) = state%triangles(mod(edge, 3) + 1, e)
            do k = 1, degree - 1
                indices(k + 1) = state%triangles(3 + (edge - 1)*(degree - 1) + k, e)
            end do
        end subroutine

        subroutine cubic_jet(e, point, value, grad, h, status)
            integer, intent(in) :: e
            real(dp), intent(in) :: point(2)
            real(dp), intent(out) :: value, grad(2), h(2, 2)
            integer, intent(out) :: status
            real(dp) :: local(2), g(2)
            if (use_reference) then
                call gs_element_jet(state, e, point(1), point(2), value, grad, h, &
                    status, affine_chart=.true.)
                return
            end if
            value = 0; grad = 0; h = 0; status = 7
            local = matmul(transpose(state%gradient(:, 2:3, e)), &
                point - state%points(:, state%triangles(1, e)))
            if (minval([1 - sum(local), local]) < -1.0e-10_dp) return
            call gs_cubic_flux(local(1), local(2), cubic(1, e), cubic(2, e), &
                cubic(3, e), cubic(4, e), cubic(5, e), cubic(6, e), cubic(7, e), &
                cubic(8, e), cubic(9, e), cubic(10, e), value, g(1), g(2))
            value = value + state%psi(state%triangles(1, e))
            grad = matmul(state%gradient(:, 2:3, e), g)
            ! P3 contour sensitivities are unavailable; the primal arc uses
            ! only value and gradient, so no Hessian is reconstructed here.
            status = 0
        end subroutine

        subroutine cubic_ray_root(e, lo, hi, direction, level, point, grad, h, status)
            integer, intent(in) :: e
            real(dp), intent(in) :: lo, hi, direction(2), level
            real(dp), intent(out) :: point(2), grad(2), h(2, 2)
            integer, intent(out) :: status
            real(dp) :: values(4), roots(3), value, radius, selected(2)
            integer :: i, count, found, stat
            do i = 1, 4
                radius = lo + (hi - lo)*real(i - 1, dp)/3
                point = axis + radius*direction
                call cubic_jet(e, point, values(i), &
                    grad, h, status)
                if (status /= 0) return
            end do
            call cubic_segment_roots(values, level, roots, count)
            found = 0
            do i = 1, count
                radius = lo + (hi - lo)*roots(i)
                if (radius <= 0) cycle
                point = axis + radius*direction
                call cubic_jet(e, point, value, &
                    grad, h, stat)
                if (stat /= 0) cycle
                if (dot_product(grad, direction) <= 0) cycle
                if (abs(value - level) > tolerance) cycle
                found = found + 1
                selected = point
            end do
            status = 7
            if (found /= 1) return
            point = selected
            call cubic_jet(e, point, value, &
                grad, h, status)
        end subroutine

        subroutine measure_factor(e, chart_point, factor, status)
            integer, intent(in) :: e
            real(dp), intent(in) :: chart_point(2)
            real(dp), intent(out) :: factor
            integer, intent(out) :: status
            real(dp) :: xi(2), point(2), jac(2, 2), curvature(2, 2, 2)
            real(dp) :: det, inverse(2, 2)
            factor = 1; status = 0
            if (.not. state%mapped) return
            if (.not. state%curved_cell(e)) return
            xi = matmul(transpose(state%gradient(:, 2:3, e)), &
                chart_point - state%points(:, state%triangles(1, e)))
            call triangle_geometry(state%points(:, state%triangles(:, e)), xi, &
                point, jac, curvature, det, inverse, status)
            if (status /= 0) return
            factor = det/(2*state%area(e))*chart_point(1)/point(1)
        end subroutine

        subroutine point_on_cell(e, theta, level, point, &
                                 grad, h, status, cw)
            integer, intent(in)::e
            integer, intent(inout)::cw
            real(dp), intent(in)::theta, level
            real(dp), intent(out)::point(2), grad(2), h(2, 2)
            integer, intent(out)::status
            real(dp)::dir(2), bary(3), slope(3), lo, hi, aa, bb, &
                       cc, root(2), r, value
            integer::k, nr, found, stat
            cw = cw + 1
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
            if (state%element_order == 3) then
                call cubic_ray_root(e, lo, hi, dir, level, point, grad, h, status)
                return
            end if
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
                                value, grad, h, stat, affine_chart=.true.)
            if (stat /= 0) return
            if (abs(value - level) > tolerance) then
                status = 12
                return
            end if
            status = 0
        end subroutine point_on_cell

        subroutine arc_weight(e, theta, level, level_delta, &
                              weight, delta_weight, status, cw)
            integer, intent(in)::e
            integer, intent(inout)::cw
            real(dp), intent(in)::theta, level, level_delta
            real(dp), intent(out)::weight, delta_weight
            integer, intent(out)::status
            real(dp)::point(2), grad(2), h(2, 2), dir(2), r, rg, droot, &
                       tv, tg(2), th(2, 2), dv, dg(2), ignored, factor
            weight = 0
            delta_weight = 0
            call point_on_cell(e, theta, level, point, grad, h, &
                               status, cw)
            if (status /= 0) return
            dir = [cos(theta), &
                   sin(theta)]
            r = norm2(point - axis)
            rg = dot_product(grad, dir)
            call gs_contour_weight(r, point(1), rg, weight)
            call measure_factor(e, point, factor, status)
            if (status /= 0) return
            weight = weight*factor
            if (.not. with_tangent) return
            call gs_element_jet(tangent_state, e, point(1), &
                                point(2), tv, tg, th, status, affine_chart=.true.)
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
                                 level_delta, status, cw, first_point, last_point)
            integer, intent(in)::e, j
            integer, intent(inout)::cw
            real(dp), intent(in)::lo, hi, dlo, dhi, level_delta
            real(dp), optional, intent(in)::first_point(2), last_point(2)
            integer, intent(out)::status
            real(dp)::weight, delta_weight, angle
            integer::q
            do q = 1, order
                angle = lo + (hi - lo)*x(q)
                call arc_weight(e, angle, levels(j), level_delta, &
                                weight, delta_weight, status, cw)
                if (status /= 0) return
                coarea(j) = coarea(j) + (hi - lo)*w(q)*weight
                if (with_tangent) coarea_tangent(j) = coarea_tangent(j) + (hi - &
            lo)*w(q)*delta_weight
            end do
            if (with_tangent) then
                ! Leibniz endpoint terms supply moving gradient-jump products;
                ! adjacent arcs cancel them only for a C1 native polynomial.
                if (present(first_point)) then
                    call face_weight(e, first_point, weight, status)
                else
                    call arc_weight(e, lo, levels(j), level_delta, &
                                    weight, delta_weight, status, cw)
                end if
                if (status /= 0) return
                coarea_tangent(j) = coarea_tangent(j) - weight*dlo
                if (present(last_point)) then
                    call face_weight(e, last_point, weight, status)
                else
                    call arc_weight(e, hi, levels(j), level_delta, &
                                    weight, delta_weight, status, cw)
                end if
                if (status /= 0) return
                coarea_tangent(j) = coarea_tangent(j) + weight*dhi
            end if
            status = 0
        end subroutine integrate_arc

        subroutine face_weight(e, face_point, weight, status)
            ! The event already owns the native face intersection. Re-solving
            ! its radial root can place a rounded endpoint outside this cell.
            integer, intent(in)::e
            real(dp), intent(in)::face_point(2)
            real(dp), intent(out)::weight
            integer, intent(out)::status
            real(dp)::value, gradient(2), hessian(2, 2), direction(2), radius, radial_gradient
            weight = 0
            call gs_element_jet(state, e, face_point(1), face_point(2), &
                                value, gradient, hessian, status, affine_chart=.true.)
            if (status /= 0) return
            direction = face_point - axis
            radius = norm2(direction)
            if (radius <= 0) then
                status = 7
                return
            end if
            direction = direction/radius
            radial_gradient = dot_product(gradient, direction)
            if (radial_gradient <= 0) then
                status = 6
                return
            end if
            call gs_contour_weight(radius, face_point(1), radial_gradient, weight)
        end subroutine face_weight

        subroutine integrate_boundary(j, level_delta, status)
            integer, intent(in)::j
            real(dp), intent(in)::level_delta
            integer, intent(out)::status
            real(dp)::point(2), vec(2), value, grad(2), h(2, 2), tv, &
                       tg(2), th(2, 2), weight, dweight, factor
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
                                            value, grad, h, stat, affine_chart=.true.)
                        if (stat /= 0) return
                        if (norm2(grad) <= tolerance) then
                            status = 19
                            return
                        end if
                        if (q == 0 .or. q == order + 1) cycle
                        call gs_boundary_weight(norm2(vec), point(1), &
                                                grad(1), grad(2), weight)
                        call measure_factor(e, point, factor, stat)
                        if (stat /= 0) return
                        coarea(j) = coarea(j) + w(q)*weight*factor
                        if (with_tangent) then
                            call gs_element_jet(tangent_state, e, point(1), &
                                                point(2), tv, tg, th, stat, &
                                                affine_chart=.true.)
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
