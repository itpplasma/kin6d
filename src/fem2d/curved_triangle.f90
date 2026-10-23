module kin6d_curved_triangle
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use kin6d_triangle_basis, only: triangle_basis_jet, TRIANGLE_MAX_NODES
    use fortfem_geometry_jet_2d, only: transform_scalar_jet_2d
    implicit none
    private
    public :: periodic_curve_t, initialize_periodic_curve, curve_value
    public :: p2_geometry, p2_shapes, p2_locate, boundary_curve
    public :: triangle_geometry, triangle_shapes, triangle_locate
    interface p2_geometry
        module procedure triangle_geometry
    end interface
    interface p2_shapes
        module procedure triangle_shapes
    end interface
    interface p2_locate
        module procedure triangle_locate
    end interface

    ! Uniform periodic parameter t in [0,1); samples omit the closing duplicate.
    type :: periodic_curve_t
        real(dp), allocatable :: points(:, :), slopes(:, :)
    end type
    abstract interface
        subroutine boundary_curve(t, point)
            import :: dp
            real(dp), intent(in) :: t
            real(dp), intent(out) :: point(2)
        end subroutine
    end interface
contains
    subroutine initialize_periodic_curve(curve, points, status)
        type(periodic_curve_t), intent(out) :: curve
        real(dp), intent(in) :: points(:, :)
        integer, intent(out) :: status
        real(dp), allocatable :: rhs(:, :), old(:, :)
        integer :: n, i, previous, next, iteration
        status = 1
        if (size(points, 1) /= 2 .or. size(points, 2) < 4) return
        if (.not. all(ieee_is_finite(points))) return
        n = size(points, 2)
        curve%points = points
        allocate (curve%slopes(2, n), rhs(2, n), old(2, n))
        do i = 1, n
            previous = modulo(i - 2, n) + 1
            next = modulo(i, n) + 1
            if (norm2(points(:, next) - points(:, i)) <= tiny(1.0_dp)) return
            rhs(:, i) = 3*(points(:, next) - points(:, previous))
        end do
        ! The cyclic cubic-spline system is strictly diagonally dominant.
        ! Jacobi contracts by 1/2, independently of the number of samples.
        curve%slopes = rhs/6
        do iteration = 1, 64
            old = curve%slopes
            do i = 1, n
                previous = modulo(i - 2, n) + 1
                next = modulo(i, n) + 1
                curve%slopes(:, i) = &
                    (rhs(:, i) - old(:, previous) - old(:, next))/4
            end do
        end do
        status = 0
    end subroutine

    pure subroutine curve_value(curve, t, point)
        type(periodic_curve_t), intent(in) :: curve
        real(dp), intent(in) :: t
        real(dp), intent(out) :: point(2)
        real(dp) :: s, u
        integer :: i, j, n
        n = size(curve%points, 2)
        s = modulo(t, 1.0_dp)*n
        i = min(n, int(s) + 1)
        j = modulo(i, n) + 1
        u = s - (i - 1)
        point = (2*u**3 - 3*u**2 + 1)*curve%points(:, i) &
                + (u**3 - 2*u**2 + u)*curve%slopes(:, i) &
                + (-2*u**3 + 3*u**2)*curve%points(:, j) &
                + (u**3 - u**2)*curve%slopes(:, j)
    end subroutine

    pure subroutine triangle_geometry(nodes, xi, point, jac, curvature, det, &
                                          inverse, status)
        real(dp), intent(in) :: nodes(:, :), xi(2)
        real(dp), intent(out) :: point(2), jac(2, 2), curvature(2, 2, 2)
        real(dp), intent(out) :: det, inverse(2, 2)
        integer, intent(out) :: status
        real(dp) :: values(TRIANGLE_MAX_NODES), gradients(2, TRIANGLE_MAX_NODES)
        real(dp) :: hessians(2, 2, TRIANGLE_MAX_NODES)
        real(dp) :: offset(2), g(2), h(2, 2), unused_g(2), unused_h(2, 2)
        integer :: j, a, n
        n = size(nodes, 2)
        call triangle_basis_jet(xi, values(:n), gradients(:, :n), hessians(:, :, :n))
        point = nodes(:, 1); jac = 0; curvature = 0
        do j = 2, n
            offset = nodes(:, j) - nodes(:, 1)
            point = point + offset*values(j)
            g = gradients(:, j)
            h = hessians(:, :, j)
            do a = 1, 2
                jac(a, :) = jac(a, :) + offset(a)*g
                curvature(:, :, a) = curvature(:, :, a) + offset(a)*h
            end do
        end do
        call transform_scalar_jet_2d(point, jac, curvature, [0.0_dp, 0.0_dp], &
            reshape([0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], [2, 2]), &
            det, inverse, unused_g, unused_h, status)
        if (point(1) <= 0) status = 1
    end subroutine

    pure subroutine triangle_shapes(nodes, xi, point, det, shape, gradient, &
                                        hessian, status)
        real(dp), intent(in) :: nodes(:, :), xi(2)
        real(dp), intent(out) :: point(2), det, shape(:), gradient(:, :)
        real(dp), intent(out) :: hessian(:, :, :)
        integer, intent(out) :: status
        real(dp) :: values(TRIANGLE_MAX_NODES), gradients(2, TRIANGLE_MAX_NODES)
        real(dp) :: hessians(2, 2, TRIANGLE_MAX_NODES)
        real(dp) :: jac(2, 2), curvature(2, 2, 2), inverse(2, 2), corrected(2, 2)
        integer :: j, n
        n = size(nodes, 2)
        call triangle_basis_jet(xi, values(:n), gradients(:, :n), hessians(:, :, :n))
        shape = 0; gradient = 0; hessian = 0
        call triangle_geometry(nodes, xi, point, jac, curvature, det, inverse, status)
        if (status /= 0) return
        do j = 1, n
            shape(j) = values(j)
            gradient(:, j) = matmul(transpose(inverse), gradients(:, j))
            corrected = hessians(:, :, j) &
                        - gradient(1, j)*curvature(:, :, 1) &
                        - gradient(2, j)*curvature(:, :, 2)
            hessian(:, :, j) = matmul(transpose(inverse), matmul(corrected, inverse))
        end do
    end subroutine

    pure subroutine triangle_locate(nodes, point, xi, status)
        real(dp), intent(in) :: nodes(:, :), point(2)
        real(dp), intent(out) :: xi(2)
        integer, intent(out) :: status
        real(dp) :: mapped(2), jac(2, 2), curvature(2, 2, 2), inverse(2, 2)
        real(dp) :: det, step(2), scale
        integer :: iteration
        xi = [1.0_dp/3, 1.0_dp/3]
        scale = max(1.0_dp, maxval(abs(nodes)))
        do iteration = 1, 12
            call triangle_geometry(nodes, xi, mapped, jac, curvature, det, &
                inverse, status)
            if (status /= 0) return
            step = matmul(inverse, mapped - point)
            xi = xi - step
            if (maxval(abs(mapped - point)) < 8*epsilon(1.0_dp)*scale) exit
        end do
        status = 1
        if (iteration > 12) return
        if (minval([xi, 1 - sum(xi)]) < -1.0e-10_dp) return
        status = 0
    end subroutine
end module
