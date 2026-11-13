module milca_nested_ring_mesh
    ! Flux-aligned reference mesh of a region with nested surfaces: rings of
    ! curved P2/P3 triangles between consecutive surfaces and a fan at the
    ! axis. Every node of a ring cell lies on the image of a regular lattice
    ! in (rho, theta), so the radial spacing can be far smaller than the
    ! poloidal one (thin layers at rational surfaces) without the minimum
    ! angle refinement of the unstructured mesher.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use milca_grad_shafranov, only: gs_state_t
    implicit none
    private
    public :: nested_ring_mesh, nested_surface

    abstract interface
        subroutine nested_surface(rho, theta, point)
            !! Point of the surface with label rho in (0, 1] at the poloidal
            !! parameter theta in [0, 1), counter-clockwise; rho = 0 is the axis.
            import :: dp
            real(dp), intent(in) :: rho, theta
            real(dp), intent(out) :: point(2)
        end subroutine
    end interface
contains
    subroutine nested_ring_mesh(surface, radii, ntheta, order, state, status)
        ! radii(1) < ... < radii(nr) = outermost label: the surfaces that carry
        ! mesh edges. Only the outermost surface is flagged as boundary.
        procedure(nested_surface) :: surface
        real(dp), intent(in) :: radii(:)
        integer, intent(in) :: ntheta, order
        type(gs_state_t), intent(out) :: state
        integer, intent(out) :: status
        real(dp) :: axis(2), outer(2)
        integer :: nr, nloc, cells, k, j, e, per, i, jj, a
        integer :: lattice(2, 10, 2), tri(10)
        status = 1
        nr = size(radii)
        if (nr < 1 .or. ntheta < 3) return
        if (order /= 2 .and. order /= 3) return
        if (.not. all(ieee_is_finite(radii))) return
        if (radii(1) <= 0) return
        if (any(radii(2:) <= radii(:nr - 1))) return
        nloc = (order + 1)*(order + 2)/2
        per = order*ntheta
        cells = ntheta*(2*nr - 1)
        ! Node (i, j): i = 0 the axis, i = order*k on surface k; j counts
        ! theta in steps of 1/per. Fan centroids follow.
        allocate (state%points(2, 1 + order*nr*per + merge(ntheta, 0, order == 3)), &
            state%triangles(nloc, cells))
        allocate (state%boundary(size(state%points, 2)))
        state%element_order = order
        state%boundary = .false.
        call surface(0.0_dp, 0.0_dp, axis)
        state%points(:, 1) = axis
        do i = 1, order*nr
            k = (i - 1)/order + 1
            do jj = 0, per - 1
                a = node(i, jj)
                if (k == 1) then
                    ! Straight radial edges of the fan, curved outer edge.
                    call surface(radii(1), real(jj, dp)/per, outer)
                    if (i < order) then
                        if (mod(jj, order) /= 0) cycle
                        state%points(:, a) = axis + real(i, dp)/order*(outer - axis)
                    else
                        state%points(:, a) = outer
                    end if
                else
                    call surface(radii(k - 1) + real(i - order*(k - 1), dp)/order &
                        *(radii(k) - radii(k - 1)), real(jj, dp)/per, &
                        state%points(:, a))
                end if
                state%boundary(a) = i == order*nr
            end do
        end do
        ! Local lattice offsets of the two triangles of a ring cell, in the
        ! node order of higher_order_mesh (vertices, edges, centroid).
        call triangle_lattice(order, [0, 0], [order, 0], [order, order], &
            lattice(:, :nloc, 1))
        call triangle_lattice(order, [0, 0], [order, order], [0, order], &
            lattice(:, :nloc, 2))
        e = 0
        do j = 0, ntheta - 1
            ! Fan: axis, surface 1 at theta_j and theta_(j+1).
            e = e + 1
            tri(1) = 1; tri(2) = node(order, order*j)
            tri(3) = node(order, order*(j + 1))
            do i = 1, order - 1
                tri(3 + i) = node(i, order*j)
                tri(3 + (order - 1) + i) = node(order, order*j + i)
                tri(3 + 2*(order - 1) + i) = node(order - i, order*(j + 1))
            end do
            if (order == 3) then
                tri(10) = 1 + order*nr*per + j + 1
                state%points(:, tri(10)) = (state%points(:, tri(1)) &
                    + state%points(:, tri(2)) + state%points(:, tri(3)))/3
            end if
            state%triangles(:, e) = tri(:nloc)
            do k = 2, nr
                do i = 1, 2
                    e = e + 1
                    do a = 1, nloc
                        state%triangles(a, e) = node(order*(k - 1) + lattice(1, a, i), &
                            order*j + lattice(2, a, i))
                    end do
                end do
            end do
        end do
        ! Unused fan lattice slots (i < order, j off the radial edges) would
        ! leave gaps in the numbering: compress.
        call compress(state)
        status = 0
    contains
        integer function node(i, jj)
            integer, intent(in) :: i, jj
            node = 2 + (i - 1)*per + modulo(jj, per)
        end function
    end subroutine nested_ring_mesh

    pure subroutine triangle_lattice(order, v1, v2, v3, offsets)
        ! Lattice points of a triangle with vertices on a lattice of step
        ! 1/order: vertices, then order - 1 nodes per edge (1->2, 2->3, 3->1),
        ! then the centroid of the cubic.
        integer, intent(in) :: order, v1(2), v2(2), v3(2)
        integer, intent(out) :: offsets(:, :)
        integer :: v(2, 4), k, i
        v(:, 1) = v1; v(:, 2) = v2; v(:, 3) = v3; v(:, 4) = v1
        offsets(:, 1:3) = v(:, 1:3)
        do k = 1, 3
            do i = 1, order - 1
                offsets(:, 3 + (k - 1)*(order - 1) + i) = v(:, k) &
                    + i*(v(:, k + 1) - v(:, k))/order
            end do
        end do
        if (order == 3) offsets(:, 10) = (v1 + v2 + v3)/3
    end subroutine triangle_lattice

    subroutine compress(state)
        type(gs_state_t), intent(inout) :: state
        integer, allocatable :: number(:)
        real(dp), allocatable :: points(:, :)
        logical, allocatable :: boundary(:)
        integer :: a, e, count
        allocate (number(size(state%points, 2)))
        number = 0
        do e = 1, size(state%triangles, 2)
            number(state%triangles(:, e)) = 1
        end do
        count = 0
        do a = 1, size(number)
            if (number(a) == 0) cycle
            count = count + 1; number(a) = count
        end do
        allocate (points(2, count), boundary(count))
        do a = 1, size(number)
            if (number(a) == 0) cycle
            points(:, number(a)) = state%points(:, a)
            boundary(number(a)) = state%boundary(a)
        end do
        do e = 1, size(state%triangles, 2)
            state%triangles(:, e) = number(state%triangles(:, e))
        end do
        call move_alloc(points, state%points)
        call move_alloc(boundary, state%boundary)
    end subroutine compress
end module milca_nested_ring_mesh
