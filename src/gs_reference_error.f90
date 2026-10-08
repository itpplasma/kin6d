! Independent known-solution validation, not a general GS certifier.
module kin6d_gs_reference_error
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortnum_interval, only: interval_t, interval, hull, operator(+), &
        operator(-), operator(*), operator(/), abs
    use kin6d_gs_reference_error_leaves, only: gs_reference_error_density_interval
    use fortfem_reference_scalar_intervals, only: generated_affine_p2_scalar_interval
    implicit none
    private
    public :: gs_reference_parameters_t, enclose_gs_reference_cell
    type :: gs_reference_parameters_t
        real(dp) :: r0, c, edge, f, mu0, pp
    end type
contains
    subroutine enclose_gs_reference_cell(vertices, nodal_psi, parameters, depth, &
            integral, area, volume, max_field2, ok)
        real(dp), intent(in) :: vertices(2, 3), nodal_psi(:)
        type(gs_reference_parameters_t), intent(in) :: parameters
        integer, intent(in) :: depth
        type(interval_t), intent(out) :: integral(3), area, volume, max_field2
        logical, intent(out) :: ok
        type(interval_t) :: points(2, 3), det, gr, gz, origin(2), ph0
        integer :: i, j
        ok = .false.
        integral = interval(0); area = interval(0); volume = interval(0)
        max_field2 = interval(0)
        if (size(nodal_psi) /= 3 .and. size(nodal_psi) /= 6) return
        if (.not. all(ieee_is_finite(vertices))) return
        if (.not. all(ieee_is_finite(nodal_psi))) return
        if (minval(vertices(1, :)) <= 0) return
        if (depth < 0 .or. depth > 8) return
        if (.not. all(ieee_is_finite([parameters%r0, parameters%c, &
            parameters%edge, parameters%f, parameters%mu0, parameters%pp]))) return
        if (parameters%mu0 <= 0) return
        do j = 1, 3
            do i = 1, 2
                points(i, j) = interval(vertices(i, j))
            end do
        end do
        det = (points(1,2)-points(1,1))*(points(2,3)-points(2,1)) &
            - (points(1,3)-points(1,1))*(points(2,2)-points(2,1))
        if (det%lo <= 0 .and. det%hi >= 0) return
        if (size(nodal_psi) == 3) then
            gr = ((interval(nodal_psi(2))-nodal_psi(1))*(points(2,3)-points(2,1)) &
                -(interval(nodal_psi(3))-nodal_psi(1))*(points(2,2)-points(2,1)))/det
            gz = ((points(1,2)-points(1,1))*(interval(nodal_psi(3))-nodal_psi(1)) &
                -(points(1,3)-points(1,1))*(interval(nodal_psi(2))-nodal_psi(1)))/det
        end if
        area = abs(det)/2
        volume = area*(points(1,1)+points(1,2)+points(1,3))/3
        origin = points(:, 1); ph0 = interval(nodal_psi(1))
        call cover(points, area, depth)
        ok = all(ieee_is_finite(integral%lo)) &
            .and. all(ieee_is_finite(integral%hi))
    contains
        recursive subroutine cover(p, weight, remaining)
            type(interval_t), intent(in) :: p(2,3), weight
            integer, intent(in) :: remaining
            type(interval_t) :: mid(2,3), child(2,3), r, z, ph, density(3)
            integer :: k
            if (remaining > 0) then
                mid(:,1) = (p(:,1)+p(:,2))/2
                mid(:,2) = (p(:,2)+p(:,3))/2
                mid(:,3) = (p(:,3)+p(:,1))/2
                child(:,1)=p(:,1); child(:,2)=mid(:,1); child(:,3)=mid(:,3)
                call cover(child, weight/4, remaining-1)
                child(:,1)=mid(:,1); child(:,2)=p(:,2); child(:,3)=mid(:,2)
                call cover(child, weight/4, remaining-1)
                child(:,1)=mid(:,3); child(:,2)=mid(:,2); child(:,3)=p(:,3)
                call cover(child, weight/4, remaining-1)
                call cover(mid, weight/4, remaining-1)
                return
            end if
            r = hull(hull(p(1,1), p(1,2)), p(1,3))
            z = hull(hull(p(2,1), p(2,2)), p(2,3))
            if (size(nodal_psi) == 3) then
                ph = ph0 + gr*(r-origin(1)) + gz*(z-origin(2))
            else
                call generated_affine_p2_scalar_interval(r, z, &
                    points(1,1), points(2,1), points(1,2), points(2,2), &
                    points(1,3), points(2,3), interval(nodal_psi(1)), &
                    interval(nodal_psi(2)), interval(nodal_psi(3)), &
                    interval(nodal_psi(4)), interval(nodal_psi(5)), &
                    interval(nodal_psi(6)), ph, gr, gz, det)
            end if
            call gs_reference_error_density_interval(r, z, interval(parameters%r0), &
                interval(parameters%c), interval(parameters%edge), &
                interval(parameters%f), interval(parameters%mu0), &
                interval(parameters%pp), ph, gr, gz, &
                density(1), density(2), density(3))
            do k = 1, 3
                ! The density is a sum of exact squares. Natural interval
                ! arithmetic may yield a negative lower endpoint through CSE.
                density(k)%lo = max(0.0_dp, density(k)%lo)
            end do
            integral = integral + weight*density
            max_field2%hi = max(max_field2%hi, density(1)%hi)
        end subroutine cover
    end subroutine enclose_gs_reference_cell
end module kin6d_gs_reference_error
