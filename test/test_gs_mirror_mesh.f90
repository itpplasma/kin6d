program test_gs_mirror_mesh
    ! initialize_gs with mirror_x: on the FB3D synthetic normal-plane loops
    ! (circle Gamma_c, chain Gamma_p, both mirror symmetric about x = 0.5)
    ! every node of the curved P2 mesh, loop nodes included, has an image
    ! under x -> 1 - x, every cell has an image cell, the loop nodes keep
    ! their input numbers, and the area is that of the curved disc. Loops
    ! that are not symmetric, or that meet the line between nodes, fail.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_grad_shafranov, only: gs_state_t, initialize_gs
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp), rc = 0.2_dp, xm = 0.5_dp
    integer, parameter :: levels(2, 3) = reshape([12, 16, 16, 20, 20, 24], [2, 3])
    type(gs_state_t) :: mesh
    real(dp), allocatable :: contour(:, :), chain(:, :)
    integer :: level, st
    do level = 1, 3
        call loops(levels(1, level), levels(2, level), 0.0_dp, contour, chain)
        call initialize_gs(mesh, contour, 0.5_dp*(2*pi*rc/size(contour, 2))**2, st, &
            2, curve=gamma_c, interior_chain=chain, mirror_x=xm)
        if (st /= 0) error stop 'mirror mesh status'
        call check(mesh, contour, chain)
    end do
    ! Rotated chain: not mirror symmetric.
    call loops(16, 20, 0.05_dp, contour, chain)
    call initialize_gs(mesh, contour, 0.5_dp*(2*pi*rc/16)**2, st, 2, curve=gamma_c, &
        interior_chain=chain, mirror_x=xm)
    if (st == 0) error stop 'asymmetric chain accepted'
    ! 10 contour nodes: the line crosses two contour segments.
    call loops(10, 20, 0.0_dp, contour, chain)
    call initialize_gs(mesh, contour, 0.5_dp*(2*pi*rc/10)**2, st, 2, curve=gamma_c, &
        interior_chain=chain, mirror_x=xm)
    if (st == 0) error stop 'line crossing a segment accepted'
    print '(a)', 'mirror mesh ok'
contains
    subroutine gamma_c(t, point)
        real(dp), intent(in) :: t
        real(dp), intent(out) :: point(2)
        point = [xm + rc*cos(2*pi*t), rc*sin(2*pi*t)]
    end subroutine

    subroutine loops(n, nc, turn, contour, chain)
        integer, intent(in) :: n, nc
        real(dp), intent(in) :: turn
        real(dp), allocatable, intent(out) :: contour(:, :), chain(:, :)
        integer :: j
        allocate (contour(2, n), chain(2, nc))
        do j = 1, n
            call gamma_c(real(j - 1, dp)/n, contour(:, j))
        end do
        do j = 1, nc
            chain(:, j) = [xm + 0.1_dp*cos(2*pi*(j - 1)/nc + turn), &
                0.1_dp*sin(2*pi*(j - 1)/nc + turn)]
        end do
    end subroutine

    subroutine check(mesh, contour, chain)
        type(gs_state_t), intent(in) :: mesh
        real(dp), intent(in) :: contour(:, :), chain(:, :)
        real(dp), parameter :: tol = 1.0e-12_dp
        real(dp) :: target(2), area
        integer, allocatable :: image(:), cell(:)
        integer :: k, j, e, f, n, nc, np, nt
        n = size(contour, 2); nc = size(chain, 2)
        np = size(mesh%points, 2); nt = size(mesh%triangles, 2)
        if (maxval(abs(mesh%points(:, :n) - contour)) > 0) error stop 'contour numbering'
        if (maxval(abs(mesh%points(:, n + 1:n + nc) - chain)) > 0) &
            error stop 'chain numbering'
        allocate (image(np), cell(nt))
        image = 0
        do k = 1, np
            target = [2*xm - mesh%points(1, k), mesh%points(2, k)]
            do j = 1, np
                if (maxval(abs(mesh%points(:, j) - target)) <= tol) image(k) = j
            end do
        end do
        if (any(image == 0)) error stop 'node without mirror image'
        cell = 0
        do e = 1, nt
            do f = 1, nt
                if (all([(any(mesh%triangles(:, f) == image(mesh%triangles(j, e))), &
                    j=1, size(mesh%triangles, 1))])) cell(e) = f
            end do
        end do
        if (any(cell == 0)) error stop 'cell without mirror image'
        if (any(cell(cell) /= [(e, e=1, nt)])) error stop 'cell image not an involution'
        if (any(mesh%boundary .neqv. mesh%boundary(image))) error stop 'boundary flags'
        area = sum(mesh%area)
        ! Straight-sided area: the inscribed n-gon of Gamma_c.
        if (abs(area - 0.5_dp*n*rc**2*sin(2*pi/n)) > 1.0e-12_dp) error stop 'area'
        print '(a,i0,a,i0,a,i0,a,i0)', 'nedge ', n, ' nchain ', nc, ' nodes ', np, &
            ' cells ', nt
    end subroutine
end program test_gs_mirror_mesh
