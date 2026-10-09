program test_curved_triangle
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_curved_triangle, only: p2_geometry, p2_shapes, p2_locate, &
        periodic_curve_t, initialize_periodic_curve, curve_value
    implicit none
    real(dp), parameter :: reference(2, 6) = reshape([ &
        0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, &
        0.5_dp, 0.0_dp, 0.5_dp, 0.5_dp, 0.0_dp, 0.5_dp], [2, 6])
    type(periodic_curve_t) :: curve
    real(dp) :: nodes(2, 6), xi(2), found(2), point(2), expected(2), jac(2, 2)
    real(dp) :: curvature(2, 2, 2), inverse(2, 2), det, shape(6), gradients(2, 6)
    real(dp) :: hessians(2, 2, 6), flux(6), hessian(2, 2), samples(2, 32)
    real(dp) :: errors(2), t
    integer :: i, j, n, status
    do j = 1, 6
        call map(reference(:, j), nodes(:, j))
    end do
    xi = [0.21_dp, 0.32_dp]
    call map(xi, expected)
    call p2_geometry(nodes, xi, point, jac, curvature, det, inverse, status)
    if (status /= 0) error stop 'valid map rejected'
    if (norm2(point - expected) > 1.0e-14_dp) error stop 'map position'
    if (maxval(abs(jac - reshape([2 + 0.2_dp*xi(2), 0.2_dp*xi(1), &
        0.2_dp*xi(1), 0.7_dp], [2, 2]))) > 1.0e-14_dp) error stop 'map Jacobian'
    call p2_locate(nodes, expected, found, status)
    if (status /= 0 .or. norm2(found - xi) > 1.0e-13_dp) error stop 'map inverse'
    call p2_shapes(nodes, xi, point, det, shape, gradients, hessians, status)
    flux = nodes(1, :) + 2*nodes(2, :)
    if (abs(dot_product(shape, flux) - point(1) - 2*point(2)) > 1.0e-13_dp) &
        error stop 'linear physical value'
    if (norm2(matmul(gradients, flux - flux(1)) - [1.0_dp, 2.0_dp]) > &
        1.0e-13_dp) error stop 'signed physical gradient'
    hessian = 0
    do j = 1, 6
        hessian = hessian + (flux(j) - flux(1))*hessians(:, :, j)
    end do
    if (maxval(abs(hessian)) > 1.0e-13_dp) error stop 'curvature connection'
    do i = 1, 2
        n = 8*2**i
        do j = 1, n
            t = 2*acos(-1.0_dp)*(j - 1)/n
            samples(:, j) = [5 + cos(t), sin(t)]
        end do
        call initialize_periodic_curve(curve, samples(:, :n), status)
        if (status /= 0) error stop 'periodic spline initialization'
        call curve_value(curve, 0.5_dp/n, point)
        t = acos(-1.0_dp)/n
        errors(i) = norm2(point - [5 + cos(t), sin(t)])
        call curve_value(curve, -0.5_dp/n, expected)
        if (abs(expected(1) - point(1)) > 1.0e-14_dp) error stop 'periodic wrap'
        if (abs(expected(2) + point(2)) > 1.0e-14_dp) error stop 'periodic parity'
    end do
    if (errors(1)/errors(2) < 14) error stop 'cubic spline fourth order'
    print *, 'PASS: generic P2 geometry, inverse, signed jets and periodic spline'
contains
    pure subroutine map(x, y)
        real(dp), intent(in) :: x(2)
        real(dp), intent(out) :: y(2)
        y = [5 + 2*x(1) + 0.2_dp*x(1)*x(2), &
            -0.4_dp + 0.7_dp*x(2) + 0.1_dp*x(1)**2]
    end subroutine
end program
