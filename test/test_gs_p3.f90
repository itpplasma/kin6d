program test_gs_p3
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_triangle_basis, only: triangle_basis_jet, cubic_segment_roots
    use kin6d_curved_triangle, only: triangle_geometry, triangle_shapes, triangle_locate
    use kin6d_grad_shafranov, only: gs_state_t, initialize_gs, solve_gs, &
        gs_linear_response, write_gs_mesh, read_gs_mesh, sample_gs
    implicit none
    real(dp), parameter :: reference(2, 10) = reshape([ &
        0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, &
        1.0_dp/3, 0.0_dp, 2.0_dp/3, 0.0_dp, 2.0_dp/3, 1.0_dp/3, &
        1.0_dp/3, 2.0_dp/3, 0.0_dp, 2.0_dp/3, 0.0_dp, 1.0_dp/3, &
        1.0_dp/3, 1.0_dp/3], [2, 10])
    real(dp), parameter :: seeds(2, 3) = reshape([ &
        6.13_dp, 0.071_dp, 6.45_dp, -0.125_dp, 6.57_dp, 0.2_dp], [2, 3])
    type(gs_state_t) :: state, restored
    real(dp) :: nodes(2, 10), shape(10), gradients(2, 10), hessians(2, 2, 10)
    real(dp) :: nodal(10), xi(2), point(2), expected(2), jac(2, 2), inverse(2, 2)
    real(dp) :: curvature(2, 2, 2), det, found(2), gradient(2), hessian(2, 2)
    real(dp) :: values(4), roots(3), t, contour(2, 32), actual(4), reloaded(4)
    real(dp), allocatable :: source(:), goal(:), tangent(:), adjoint(:)
    integer :: i, j, px, py, status, count
    do i = 1, 10
        call triangle_basis_jet(reference(:, i), shape, gradients, hessians)
        shape(i) = shape(i) - 1
        if (maxval(abs(shape)) > 2.0e-14_dp) error stop 'P3 cardinality'
        call map(reference(:, i), nodes(:, i))
    end do
    xi = [0.21_dp, 0.32_dp]
    call triangle_basis_jet(xi, shape, gradients, hessians)
    do px = 0, 3
        do py = 0, 3 - px
            nodal = reference(1, :)**px*reference(2, :)**py
            if (abs(dot_product(shape, nodal) - xi(1)**px*xi(2)**py) &
                > 2.0e-14_dp) error stop 'P3 polynomial value'
            gradient = 0; hessian = 0
            if (px > 0) gradient(1) = px*xi(1)**(px - 1)*xi(2)**py
            if (py > 0) gradient(2) = py*xi(1)**px*xi(2)**(py - 1)
            if (px > 1) hessian(1, 1) = px*(px - 1)*xi(1)**(px - 2)*xi(2)**py
            if (py > 1) hessian(2, 2) = py*(py - 1)*xi(1)**px*xi(2)**(py - 2)
            if (px > 0 .and. py > 0) hessian(1, 2) = &
                px*py*xi(1)**(px - 1)*xi(2)**(py - 1)
            hessian(2, 1) = hessian(1, 2)
            if (norm2(matmul(gradients, nodal) - gradient) > 3.0e-14_dp) &
                error stop 'P3 polynomial gradient'
            do j = 1, 10
                hessian = hessian - nodal(j)*hessians(:, :, j)
            end do
            if (maxval(abs(hessian)) > 8.0e-14_dp) error stop 'P3 Hessian'
        end do
    end do
    call map(xi, expected)
    call triangle_geometry(nodes, xi, point, jac, curvature, det, inverse, status)
    if (status /= 0 .or. norm2(point - expected) > 2.0e-14_dp) error stop 'cubic map'
    call triangle_locate(nodes, expected, found, status)
    if (status /= 0 .or. norm2(found - xi) > 2.0e-13_dp) error stop 'cubic inverse'
    call triangle_shapes(nodes, xi, point, det, shape, gradients, hessians, status)
    if (status /= 0) error stop 'cubic physical shapes'
    nodal = 2*nodes(1, :) - 3*nodes(2, :)
    if (norm2(matmul(gradients, nodal - nodal(1)) - [2.0_dp, -3.0_dp]) &
        > 2.0e-13_dp) error stop 'cubic physical gradient'
    hessian = 0
    do j = 1, 10
        hessian = hessian + (nodal(j) - nodal(1))*hessians(:, :, j)
    end do
    if (maxval(abs(hessian)) > 3.0e-13_dp) error stop 'cubic connection'
    do i = 1, 4
        t = real(i - 1, dp)/3
        values(i) = (t - 0.2_dp)*(t - 0.5_dp)*(t - 0.8_dp)
    end do
    call cubic_segment_roots(values, 0.0_dp, roots, count)
    if (count /= 3) error stop 'all cubic intersections'
    if (maxval(abs(roots - [0.2_dp, 0.5_dp, 0.8_dp])) > 1.0e-13_dp) &
        error stop 'cubic intersection values'
    do i = 1, 32
        t = 2*acos(-1.0_dp)*(i - 1)/32
        contour(:, i) = [6.2_dp + cos(t), sin(t)]
    end do
    call initialize_gs(state, contour, 0.08_dp, status, 3, interior_points=seeds)
    if (status /= 0) error stop 'affine P3 mesh'
    do i = 1, size(seeds, 2)
        j = minloc(sum((state%points - spread(seeds(:, i), 2, &
            size(state%points, 2)))**2, dim=1), dim=1)
        if (norm2(state%points(:, j) - seeds(:, i)) > 1.0e-14_dp) &
            error stop 'interior seed missing'
        if (state%boundary(j)) error stop 'interior seed constrained as boundary'
    end do
    where (state%boundary) state%psi = state%points(2, :)**2
    call solve_gs(state, 1.0_dp, 0.0_dp, -2.0_dp, status)
    if (status /= 0) error stop 'P3 weak source'
    if (maxval(abs(state%psi - state%points(2, :)**2)) > 2.0e-8_dp) &
        error stop 'independent P3 quadratic weak solution'
    allocate (source(size(state%psi)), goal(size(state%psi)))
    do i = 1, size(source)
        source(i) = sin(real(i, dp)); goal(i) = cos(real(i, dp))
    end do
    call gs_linear_response(state, source, tangent, status)
    if (status /= 0) error stop 'P3 equation tangent'
    call gs_linear_response(state, goal, adjoint, status)
    if (status /= 0) error stop 'P3 equation adjoint'
    if (abs(dot_product(goal, tangent) - dot_product(adjoint, source)) &
        > 1.0e-8_dp) error stop 'P3 equation duality'
    call write_gs_mesh(state, 'p3_native_mesh.dat', status)
    if (status /= 0) error stop 'P3 export'
    call read_gs_mesh(restored, 'p3_native_mesh.dat', status)
    if (status /= 0) error stop 'P3 readback'
    call sample_gs(state, 6.21_dp, 0.17_dp, -30.0_dp, actual(1), actual(2), &
        actual(3), actual(4), status)
    if (status /= 0) error stop 'P3 sample'
    call sample_gs(restored, 6.21_dp, 0.17_dp, -30.0_dp, reloaded(1), reloaded(2), &
        reloaded(3), reloaded(4), status)
    if (status /= 0) error stop 'P3 reloaded sample'
    if (maxval(abs(actual - reloaded)) > 1.0e-13_dp) error stop 'P3 round trip'
contains
    pure subroutine map(x, y)
        real(dp), intent(in) :: x(2)
        real(dp), intent(out) :: y(2)
        y = [5 + 0.6_dp*x(1) + 0.1_dp*x(2) + 0.05_dp*x(1)*x(2)*(1 - sum(x)), &
            -0.1_dp + 0.2_dp*x(1) + 0.5_dp*x(2) + 0.02_dp*x(1)*x(2)**2]
    end subroutine
end program
