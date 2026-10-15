module kin6d_triangle_basis
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortfem_generated_p1_basis, only: generated_p1_jet
    use fortfem_generated_p2_basis, only: generated_p2_jet
    use fortfem_generated_cardinal_product_jet1_order2, only: &
        generated_cardinal_product_jet1_order2
    use fortfem_generated_factor_product_jet2_order1, only: &
        generated_factor_product_jet2_order1
    use fortfem_generated_factor_product_jet3_order2, only: &
        generated_factor_product_jet3_order2
    use fortnum_roots, only: root_brent
    use fortnum_status, only: fortnum_status_t, FORTNUM_OK
    implicit none
    private
    public :: triangle_basis_jet, cubic_segment_roots, cubic_segment_stationary
    public :: triangle_bernstein_controls
    real(dp), parameter :: segment_nodes(4) = [0.0_dp, 1.0_dp/3, 2.0_dp/3, 1.0_dp]
    ! Derivatives of the cubic Lagrange basis on segment_nodes at t = 0, 1/2,
    ! 1 (columns), and the barycentric constants 1/prod_(k/=i)(x_i - x_k);
    ! all exact binary fractions.
    real(dp), parameter :: segment_slopes(4, 3) = reshape([ &
        -5.5_dp, 9.0_dp, -4.5_dp, 1.0_dp, &
        0.125_dp, -3.375_dp, 3.375_dp, -0.125_dp, &
        -1.0_dp, 4.5_dp, -9.0_dp, 5.5_dp], [4, 3])
    real(dp), parameter :: segment_barycentric(4) = &
        [-4.5_dp, 13.5_dp, -13.5_dp, 4.5_dp]
contains
    pure subroutine triangle_bernstein_controls(values, controls)
        real(dp), intent(in) :: values(:)
        real(dp), intent(out) :: controls(:)
        integer :: j, next, edge
        controls = values
        if (size(values) == 6) then
            do j = 1, 3
                next = mod(j, 3) + 1
                controls(j + 3) = 2*values(j + 3) - (values(j) + values(next))/2
            end do
        else if (size(values) == 10) then
            do j = 1, 3
                next = mod(j, 3) + 1; edge = 2*j + 2
                controls(edge) = (-5*values(j) + 18*values(edge) &
                    - 9*values(edge + 1) + 2*values(next))/6
                controls(edge + 1) = (2*values(j) - 9*values(edge) &
                    + 18*values(edge + 1) - 5*values(next))/6
            end do
            controls(10) = (27*values(10) - sum(values(:3)) &
                - 3*sum(controls(4:9)))/6
        end if
    end subroutine
    pure subroutine triangle_basis_jet(xi, value, gradient, hessian)
        real(dp), intent(in) :: xi(2)
        real(dp), intent(out) :: value(:), gradient(:, :)
        real(dp), intent(out), optional :: hessian(:, :, :)
        integer, parameter :: powers(3, 10) = reshape([ &
            3, 0, 0, 0, 3, 0, 0, 0, 3, &
            2, 1, 0, 1, 2, 0, 0, 2, 1, 0, 1, 2, 1, 0, 2, 2, 0, 1, &
            1, 1, 1], [3, 10])
        real(dp) :: v(0:3, 3), g(0:3, 3), h(0:3, 3), lambda(3)
        real(dp) :: next_v, next_g, next_h, unused, unused_z
        real(dp) :: local_hessian(2, 2)
        integer :: i, j, k, a, b, c
        if (size(value) == 3) then
            do i = 1, 3
                call generated_p1_jet(i, xi(1), xi(2), value(i), &
                    gradient(:, i), local_hessian)
                if (present(hessian)) hessian(:, :, i) = local_hessian
            end do
        else if (size(value) == 6) then
            do i = 1, 6
                call generated_p2_jet(i, xi(1), xi(2), value(i), &
                    gradient(:, i), local_hessian)
                if (present(hessian)) hessian(:, :, i) = local_hessian
            end do
        else
            ! Use FortFEM's generated cardinal products. Its arbitrary-order
            ! triangle API lacks Hessians and allocates on each evaluation.
            ! Restrict the generated tetrahedral Hessian action to z=0 with
            ! a constant fourth factor, keeping the same symbolic definition.
            lambda = [1 - sum(xi), xi]
            v(0, :) = 1; g(0, :) = 0; h(0, :) = 0
            do j = 1, 3
                do k = 1, 3
                    call generated_cardinal_product_jet1_order2(v(k - 1, j), &
                        g(k - 1, j), h(k - 1, j), 3.0_dp, real(k - 1, dp), &
                        lambda(j), next_v, next_g, next_h)
                    v(k, j) = next_v/k; g(k, j) = next_g/k; h(k, j) = next_h/k
                end do
            end do
            do i = 1, 10
                a = powers(1, i); b = powers(2, i); c = powers(3, i)
                call generated_factor_product_jet2_order1(v(a, 1), g(a, 1), &
                    v(b, 2), g(b, 2), v(c, 3), g(c, 3), value(i), &
                    gradient(1, i), gradient(2, i))
                if (.not. present(hessian)) cycle
                do j = 1, 2
                    call generated_factor_product_jet3_order2( &
                        v(a, 1), g(a, 1), h(a, 1), v(b, 2), g(b, 2), h(b, 2), &
                        v(c, 3), g(c, 3), h(c, 3), 1.0_dp, 0.0_dp, 0.0_dp, &
                        real(2 - j, dp), real(j - 1, dp), 0.0_dp, unused, &
                        hessian(1, j, i), hessian(2, j, i), unused_z)
                end do
            end do
        end if
    end subroutine

    subroutine cubic_segment_stationary(values, roots, count)
        real(dp), intent(in) :: values(4)
        real(dp), intent(out) :: roots(2)
        integer, intent(out) :: count
        real(dp) :: g(3), a, b, c, discriminant, q, trial(2)
        integer :: j
        do j = 1, 3
            g(j) = dot_product(segment_slopes(:, j), values - values(1))
        end do
        a = 2*(g(1) + g(3) - 2*g(2)); b = g(3) - g(1) - a; c = g(1)
        count = 0; roots = 0
        if (abs(a) <= 32*epsilon(1.0_dp)*max(abs(b), abs(c))) then
            if (abs(b) <= tiny(1.0_dp)) return
            trial = -c/b
        else
            discriminant = b*b - 4*a*c
            if (discriminant < 0) return
            q = -(b + sign(sqrt(discriminant), b))/2
            if (abs(q) <= tiny(1.0_dp)) return
            trial = [q/a, c/q]
        end if
        if (trial(2) < trial(1)) trial = trial(2:1:-1)
        do j = 1, 2
            if (trial(j) <= 0 .or. trial(j) >= 1) cycle
            if (count > 0) then
                if (abs(trial(j) - roots(count)) < 32*epsilon(1.0_dp)) cycle
            end if
            count = count + 1
            roots(count) = trial(j)
        end do
    end subroutine

    subroutine cubic_segment_roots(values, level, roots, count)
        real(dp), intent(in) :: values(4), level
        real(dp), intent(out) :: roots(3)
        integer, intent(out) :: count
        real(dp) :: critical(2), knots(4), left, right, root, tolerance
        integer :: nc, j
        type(fortnum_status_t) :: status
        call cubic_segment_stationary(values, critical, nc)
        knots(1) = 0; knots(nc + 2) = 1
        knots(2:nc + 1) = critical(:nc)
        tolerance = 16*epsilon(1.0_dp)*max(maxval(abs(values - level)), tiny(1.0_dp))
        count = 0; roots = 0
        do j = 1, nc + 1
            left = polynomial(knots(j)); right = polynomial(knots(j + 1))
            if (abs(left) <= tolerance) then
                root = knots(j)
            else if (abs(right) <= tolerance) then
                root = knots(j + 1)
            else
                if ((left > 0) .eqv. (right > 0)) cycle
                call root_brent(polynomial, knots(j), knots(j + 1), root, status, &
                    16*epsilon(1.0_dp), tolerance, 64)
                if (status%code /= FORTNUM_OK) cycle
            end if
            if (count > 0) then
                if (abs(root - roots(count)) < 32*epsilon(1.0_dp)) cycle
            end if
            count = count + 1
            roots(count) = root
        end do
    contains
        pure function polynomial(t) result(value)
            real(dp), intent(in) :: t
            real(dp) :: value, d(4)
            d = t - segment_nodes
            value = segment_barycentric(1)*d(2)*d(3)*d(4)*(values(1) - level) &
                + segment_barycentric(2)*d(1)*d(3)*d(4)*(values(2) - level) &
                + segment_barycentric(3)*d(1)*d(2)*d(4)*(values(3) - level) &
                + segment_barycentric(4)*d(1)*d(2)*d(3)*(values(4) - level)
        end function
    end subroutine
end module
