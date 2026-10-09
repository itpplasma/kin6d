module kin6d_solovev_cell_weak
    !! Candidate finite-node forms on T = G o A; no integration certificate.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use basis_p2_2d_module, only: basis_p2_2d_t
    use kin6d_grad_shafranov_leaves, only: gs_source
    use kin6d_solovev_geometry, only: solovev_geometry_t, geometry_sample_t, &
        geometry_ok, geometry_invalid, geometry_arithmetic
    use kin6d_solovev_cell_jet, only: solovev_cell_scalar_jet
    implicit none
    private
    public :: solovev_cell_p2_forms
contains
    subroutine solovev_cell_p2_forms(map, vertices, nodes, weights, mu0, &
            pprime, ffprime, stiffness, load, status, certified)
        !! Weights integrate dxi deta on Delta (area 1/2), NOT a unit-area rule.
        !! Strictly positive weights preserve the ideal finite-sum PSD structure.
        !! No sum/degree check, integral-error bound, or PDE admission is implied.
        !! At fixed inputs the coefficient residual is K*u-b, JVP K*du, VJP K^T*z
        !! in Euclidean coefficient duals; geometry/parameter products unavailable.
        type(solovev_geometry_t), intent(in) :: map
        real(dp), intent(in) :: vertices(2, 3), nodes(:, :), weights(:)
        real(dp), intent(in) :: mu0, pprime, ffprime
        real(dp), intent(out) :: stiffness(6, 6), load(6)
        integer, intent(out) :: status
        logical, intent(out) :: certified
        type(basis_p2_2d_t) :: basis
        type(geometry_sample_t) :: geometry
        real(dp) :: matrix(6, 6), rhs(6), shape(6), reference(2, 6)
        real(dp) :: gradient(2, 6), inverse(2, 2), determinant, radius, source
        real(dp) :: unused_gradient(2), unused_hessian(2, 2), zero_hessian(2, 2)
        real(dp) :: area_weight, energy_weight, xi(2)
        integer :: q, i, j

        stiffness = 0; load = 0; certified = .false.
        status = geometry_invalid
        if (size(nodes, 1) /= 2) return
        if (size(weights) < 1) return
        if (size(nodes, 2) /= size(weights)) return
        if (.not. all(ieee_is_finite(vertices))) return
        if (.not. all(ieee_is_finite(nodes))) return
        if (.not. all(ieee_is_finite(weights))) return
        if (any(weights <= 0)) return
        if (.not. ieee_is_finite(mu0)) return
        if (mu0 <= 0) return
        if (.not. ieee_is_finite(pprime)) return
        if (.not. ieee_is_finite(ffprime)) return
        matrix = 0; rhs = 0; zero_hessian = 0
        do q = 1, size(weights)
            ! One common geometry sample/inverse for all six trial/test functions.
            xi(1) = nodes(1, q); xi(2) = nodes(2, q)
            call solovev_cell_scalar_jet(map, vertices, xi, &
                [0.0_dp, 0.0_dp], zero_hessian, geometry, determinant, inverse, &
                unused_gradient, unused_hessian, status)
            if (status /= geometry_ok) return
            radius = geometry%value(1)
            status = geometry_arithmetic
            if (radius <= 0) return
            do i = 1, 6
                shape(i) = basis%eval(i, nodes(1, q), nodes(2, q))
                reference(:, i) = basis%grad(i, nodes(1, q), nodes(2, q))
            end do
            gradient = matmul(transpose(inverse), reference)
            ! gs_source already contains the physical weak source's 1/R factor.
            call gs_source(radius, mu0, pprime, ffprime, source)
            area_weight = weights(q)*determinant
            energy_weight = area_weight/radius
            if (.not. all(ieee_is_finite(gradient))) return
            if (.not. ieee_is_finite(source)) return
            if (.not. ieee_is_finite(area_weight)) return
            if (.not. ieee_is_finite(energy_weight)) return
            if (area_weight <= 0 .or. energy_weight <= 0) return
            do j = 1, 6
                rhs(j) = rhs(j) + area_weight*source*shape(j)
                do i = 1, 6
                    matrix(i, j) = matrix(i, j) &
                        + energy_weight*dot_product(gradient(:, i), gradient(:, j))
                end do
            end do
            if (.not. all(ieee_is_finite(matrix))) return
            if (.not. all(ieee_is_finite(rhs))) return
        end do
        stiffness = matrix; load = rhs; status = geometry_ok
    end subroutine
end module
