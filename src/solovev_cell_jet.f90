module kin6d_solovev_cell_jet
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use basis_p1_2d_module, only: basis_p1_2d_t
    use fortfem_geometry_jet_2d, only: transform_scalar_jet_2d, &
        geometry_jet_success, geometry_jet_invalid_input
    use kin6d_solovev_geometry, only: solovev_geometry_t, geometry_sample_t, &
        sample_solovev_geometry, geometry_ok, geometry_invalid, geometry_arithmetic
    implicit none
    private
    public :: solovev_cell_scalar_jet
contains
    subroutine solovev_cell_scalar_jet(map, vertices, xi, reference_gradient, &
            reference_hessian, geometry, determinant, inverse_jacobian, &
            gradient, hessian, status)
        type(solovev_geometry_t), intent(in) :: map
        real(dp), intent(in) :: vertices(2, 3), xi(2)
        real(dp), intent(in) :: reference_gradient(2), reference_hessian(2, 2)
        type(geometry_sample_t), intent(out) :: geometry
        real(dp), intent(out) :: determinant, inverse_jacobian(2, 2)
        real(dp), intent(out) :: gradient(2), hessian(2, 2)
        integer, intent(out) :: status
        type(basis_p1_2d_t) :: affine_basis
        type(geometry_sample_t) :: sample
        real(dp) :: point(2), affine(2, 2), affine_det
        real(dp) :: composed_jacobian(2, 2), composed_hessian(2, 2, 2)
        real(dp) :: fortfem_hessian(2, 2, 2), det, inverse(2, 2)
        real(dp) :: physical_gradient(2), physical_hessian(2, 2)
        integer :: a, i, j, jet_status

        geometry = geometry_sample_t()
        determinant = 0; inverse_jacobian = 0; gradient = 0; hessian = 0
        status = geometry_invalid
        if (.not. all(ieee_is_finite(vertices))) return
        if (.not. all(ieee_is_finite(xi))) return
        if (.not. all(ieee_is_finite(reference_gradient))) return
        if (.not. all(ieee_is_finite(reference_hessian))) return
        if (any(xi < 0) .or. sum(xi) > 1) return
        call affine_basis%compute_jacobian(vertices, affine, affine_det)
        if (.not. all(ieee_is_finite(affine))) return
        if (.not. ieee_is_finite(affine_det)) return
        if (affine_det <= 0) return
        call affine_basis%transform_to_physical(xi(1), xi(2), vertices, &
            point(1), point(2))
        if (.not. all(ieee_is_finite(point))) return

        call sample_solovev_geometry(map, point, sample)
        ! Retain diagnostic metadata, without publishing partial numerical jets.
        geometry%facet = sample%facet
        geometry%identity = sample%identity
        geometry%inner_interface = sample%inner_interface
        geometry%sector_interface = sample%sector_interface
        status = sample%status; geometry%status = status
        if (status /= geometry_ok) return

        ! T = G o A: KIN uses H(a,i,j), FortFEM uses H(i,j,a).
        composed_jacobian = matmul(sample%jacobian, affine)
        do a = 1, 2
            composed_hessian(a, :, :) = matmul(transpose(affine), &
                matmul(sample%hessian(a, :, :), affine))
            do j = 1, 2
                do i = 1, 2
                    fortfem_hessian(i, j, a) = composed_hessian(a, i, j)
                end do
            end do
        end do
        status = geometry_arithmetic; geometry%status = status
        if (.not. all(ieee_is_finite(composed_jacobian))) return
        if (.not. all(ieee_is_finite(composed_hessian))) return
        call transform_scalar_jet_2d(sample%value, composed_jacobian, fortfem_hessian, &
            reference_gradient, reference_hessian, det, inverse, &
            physical_gradient, physical_hessian, jet_status)
        if (jet_status /= geometry_jet_success) then
            ! FortFEM's nonpositive-J status2 is not KIN's interface status2.
            if (jet_status == geometry_jet_invalid_input) status = geometry_invalid
            geometry%status = status
            return
        end if
        geometry = sample
        geometry%jacobian = composed_jacobian
        geometry%hessian = composed_hessian
        determinant = det; inverse_jacobian = inverse
        gradient = physical_gradient; hessian = physical_hessian
        status = geometry_ok
    end subroutine
end module
