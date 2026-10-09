module kin6d_gs_field_readback
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use kin6d_grad_shafranov, only: gs_state_t, gs_element_jet, sample_gs
    use kin6d_grad_shafranov_leaves, only: gs_field
    use delaunay_types, only: mesh_t
    use point_location, only: locate_point, LOCATION_INSIDE
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
        integer, allocatable :: owner(:), owner_edge(:)
        real(dp) :: grad(2), hessian(2, 2), local(2)
        integer :: e, k, key, other, i, cell, location, seed, degree

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
            ! Preserve first-containing-cell semantics on faces and retain
            ! the physical inverse map outside the affine boundary polygon.
            call sample_gs(state, points(1, i), points(2, i), f, &
                values(1, i), values(2, i), values(3, i), values(4, i), statuses(i))
        end do
    end subroutine
end module
