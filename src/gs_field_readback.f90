module kin6d_gs_field_readback
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use kin6d_grad_shafranov, only: gs_state_t, gs_element_jet, sample_gs
    use kin6d_grad_shafranov_leaves, only: gs_field
    use delaunay_types, only: mesh_t, point_t
    use geometric_predicates, only: orientation, ORIENTATION_CCW, ORIENTATION_CW
    use point_location, only: locate_point, LOCATION_INSIDE, LOCATION_NOT_FOUND
    implicit none
    private
    public :: sample_gs_many
contains
    subroutine sample_gs_many(state, points, f, values, statuses)
        type(gs_state_t), intent(in) :: state
        real(dp), intent(in) :: points(:, :), f
        real(dp), intent(out) :: values(:, :)
        integer, intent(out) :: statuses(:)
        type(mesh_t) :: locator
        integer, allocatable :: owner(:), owner_edge(:), boundary_edges(:, :)
        type(point_t) :: query
        real(dp) :: grad(2), hessian(2, 2), local(2)
        real(dp) :: xmin, xmax, ymin, ymax, band, dx, dy, edge_len2
        real(dp) :: fraction, rx, ry
        integer :: e, k, key, other, i, cell, location, seed, degree
        integer :: nboundary, crossings, side, a, b
        logical :: boundary_valid, ambiguous, upward, downward

        ! Adapt the native FE connectivity to FortFEM's walking locator once
        ! per batch. Curved geometry and field values remain in the GS state.
        values = 0; statuses = 1
        if (size(points, 1) /= 2 .or. size(values, 1) /= 4) return
        if (size(values, 2) /= size(points, 2)) return
        if (size(statuses) /= size(points, 2)) return
        if (.not. ieee_is_finite(f)) return
        degree = state%element_order
        if (degree < 2) then
            do i = 1, size(points, 2)
                call sample_gs(state, points(1, i), points(2, i), f, &
                    values(1, i), values(2, i), values(3, i), values(4, i), statuses(i))
            end do
            return
        end if
        locator%npoints = size(state%points, 2)
        locator%ntriangles = size(state%triangles, 2)
        allocate (locator%points(locator%npoints), &
            locator%triangles(locator%ntriangles), owner(locator%npoints), &
            owner_edge(locator%npoints))
        owner = 0; owner_edge = 0
        locator%points%x = state%points(1, :)
        locator%points%y = state%points(2, :)
        do e = 1, locator%ntriangles
            locator%triangles(e)%vertices = state%triangles(1:3, e)
            locator%triangles(e)%neighbors = 0
        end do
        do e = 1, locator%ntriangles
            do k = 1, 3
                ! Shared edge DOFs provide an O(N) adjacency map, also after
                ! mesh-file readback (which deliberately has no solve operator).
                key = minval(state%triangles(4 + (k - 1)*(degree - 1): &
                    3 + k*(degree - 1), e))
                other = owner(key)
                if (other == 0) then
                    owner(key) = e; owner_edge(key) = k
                else
                    locator%triangles(e)%neighbors(k) = other
                    locator%triangles(other)%neighbors(owner_edge(key)) = e
                end if
            end do
        end do
        ! Cover the scalar evaluator's 1e-10 barycentric edge allowance.
        band = 2.0e-10_dp*norm2([maxval(state%points(1, :)) &
            - minval(state%points(1, :)), maxval(state%points(2, :)) &
            - minval(state%points(2, :))])
        nboundary = 0
        do e = 1, locator%ntriangles
            nboundary = nboundary + count(locator%triangles(e)%neighbors == 0)
        end do
        allocate (boundary_edges(2, nboundary))
        boundary_valid = nboundary > 0
        nboundary = 0
        do e = 1, locator%ntriangles
            do k = 1, 3
                if (locator%triangles(e)%neighbors(k) /= 0) cycle
                nboundary = nboundary + 1
                a = locator%triangles(e)%vertices(k)
                b = locator%triangles(e)%vertices(mod(k, 3) + 1)
                boundary_edges(:, nboundary) = [a, b]
                if (.not. state%boundary(a) .or. .not. state%boundary(b)) &
                    boundary_valid = .false.
            end do
        end do
        seed = 1
        do i = 1, size(points, 2)
            if (.not. all(ieee_is_finite(points(:, i)))) cycle
            call locate_point(locator, points(1, i), points(2, i), cell, location, seed)
            if (location == LOCATION_INSIDE) then
                local = matmul(transpose(state%gradient(:, 2:3, cell)), &
                    points(:, i) - state%points(:, state%triangles(1, cell)))
                if (minval([1 - sum(local), local]) <= 1.0e-10_dp) location = 0
            end if
            if (location == LOCATION_INSIDE) then
                call gs_element_jet(state, cell, points(1, i), points(2, i), &
                    values(1, i), grad, hessian, statuses(i))
                if (statuses(i) == 0) then
                    call gs_field(points(1, i), grad(1), grad(2), f, &
                        values(2, i), values(3, i), values(4, i))
                    seed = cell
                    cycle
                end if
            end if
            ambiguous = .true.
            crossings = 0
            if (location == LOCATION_NOT_FOUND .and. boundary_valid) then
                ! A walk can exit and reenter a concave or disconnected mesh.
                ! Even-odd boundary crossings prove affine absence globally.
                ambiguous = .false.
                query%x = points(1, i); query%y = points(2, i)
                do k = 1, nboundary
                    a = boundary_edges(1, k); b = boundary_edges(2, k)
                    xmin = min(locator%points(a)%x, locator%points(b)%x)
                    xmax = max(locator%points(a)%x, locator%points(b)%x)
                    ymin = min(locator%points(a)%y, locator%points(b)%y)
                    ymax = max(locator%points(a)%y, locator%points(b)%y)
                    if (query%x >= xmin - band .and. query%x <= xmax + band .and. &
                        query%y >= ymin - band .and. query%y <= ymax + band) then
                        dx = locator%points(b)%x - locator%points(a)%x
                        dy = locator%points(b)%y - locator%points(a)%y
                        edge_len2 = dx*dx + dy*dy
                        if (edge_len2 <= tiny(1.0_dp)) then
                            ambiguous = .true.
                            exit
                        end if
                        fraction = max(0.0_dp, min(1.0_dp, &
                            ((query%x - locator%points(a)%x)*dx &
                            + (query%y - locator%points(a)%y)*dy)/edge_len2))
                        rx = query%x - locator%points(a)%x - fraction*dx
                        ry = query%y - locator%points(a)%y - fraction*dy
                        if (rx*rx + ry*ry <= band*band) then
                            ambiguous = .true.
                            exit
                        end if
                    end if
                    if (query%y < ymin .or. query%y > ymax) cycle
                    upward = locator%points(a)%y <= query%y .and. &
                        query%y < locator%points(b)%y
                    downward = locator%points(b)%y <= query%y .and. &
                        query%y < locator%points(a)%y
                    if (.not. upward .and. .not. downward) cycle
                    side = orientation(locator%points(a), locator%points(b), query)
                    if (upward .and. side == ORIENTATION_CCW) crossings = crossings + 1
                    if (downward .and. side == ORIENTATION_CW) crossings = crossings + 1
                end do
            end if
            if (location == LOCATION_NOT_FOUND .and. .not. ambiguous .and. &
                mod(crossings, 2) == 0) then
                ! Only curved cells can extend beyond the affine triangulation.
                ! Keep their physical inverse map, without scanning every cell
                ! for the exterior points of a rectangular export grid.
                do e = 1, locator%ntriangles
                    if (.not. state%curved_cell(e)) cycle
                    call gs_element_jet(state, e, points(1, i), points(2, i), &
                        values(1, i), grad, hessian, statuses(i))
                    if (statuses(i) /= 0) cycle
                    call gs_field(points(1, i), grad(1), grad(2), f, &
                        values(2, i), values(3, i), values(4, i))
                    exit
                end do
                if (statuses(i) /= 0) values(:, i) = 0
                cycle
            end if
            ! Preserve first-containing-cell semantics on faces and retain
            ! the physical inverse map for ambiguous affine-interior points.
            call sample_gs(state, points(1, i), points(2, i), f, &
                values(1, i), values(2, i), values(3, i), values(4, i), statuses(i))
        end do
    end subroutine
end module
