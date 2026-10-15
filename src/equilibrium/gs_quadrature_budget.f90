module kin6d_gs_quadrature_budget
    !! Outward-rounded continuous energy audit for a concrete production source:
    !! -div(grad psi/R)=mu0*R*pprime+ffprime/R on an affine polygon, psi=0 on
    !! its boundary. The candidate can be affine P1 or P2. Its sampled RTp
    !! reconstruction (p = element order) only supplies candidate
    !! coefficients; the independent FortFEM checker defines an exactly H(div)
    !! degree-(p+1) field with shared edge normal components and encloses its
    !! full residual and flux defect, including quadrature and arithmetic.
    !! The source is passed exactly as (mu0*pprime*R^2+ffprime)/R.
    !! Arbitrary callbacks, nonlinear profiles, curved cells, and pointwise
    !! error bands are outside this claim. This bound says nothing about error
    !! between the polygon and an intended curved physical domain.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortnum_interval, only: interval_t, interval, operator(-)
    use fortfem_equilibrated_flux_2d, only: equilibrated_flux_estimate_t
    use fortfem_verified_bernstein_majorant, only: verified_bernstein_majorant, &
        VERIFIED_BERNSTEIN_OK
    use kin6d_grad_shafranov, only: gs_state_t
    use kin6d_gs_error_leaves, only: gs_error_friedrichs_interval, &
                                     gs_error_source_numerator_interval
    implicit none
    private
    public :: gs_checked_equilibrated_majorant
    integer, parameter, public :: GS_CHECKED_OK = 0, GS_CHECKED_UNAVAILABLE = 1
contains
    subroutine gs_checked_equilibrated_majorant(state, psi, estimate, &
                                  mu0, pprime, ffprime, energy_bound, status)
        type(gs_state_t), intent(in) :: state
        real(dp), intent(in) :: psi(:), mu0, pprime, ffprime
        type(equilibrated_flux_estimate_t), intent(in) :: estimate
        type(interval_t), intent(out) :: energy_bound
        integer, intent(out) :: status
        type(interval_t) :: lower(2), upper(2), widths(2), cf, flux, residual
        integer :: nd, checked

        status = GS_CHECKED_UNAVAILABLE
        energy_bound = interval(0.0_dp, huge(1.0_dp))
        if (.not. state%ready) return
        if (.not. allocated(state%points) .or. .not. allocated(state%triangles)) return
        if (.not. allocated(state%boundary)) return
        if (size(psi) /= size(state%points, 2)) return
        if (size(state%boundary) /= size(psi)) return
        if (state%element_order < 1 .or. state%element_order > 2) return
        if (allocated(state%curved_cell)) then
            if (any(state%curved_cell)) return
        end if
        if (any(psi /= 0.0_dp .and. state%boundary)) return
        if (.not. all(ieee_is_finite(psi))) return
        if (.not. all(ieee_is_finite([mu0, pprime, ffprime]))) return
        if (mu0 <= 0) return
        if (.not. allocated(estimate%flux_coefficients)) return
        if (size(estimate%flux_coefficients, 2) /= size(state%triangles, 2)) return
        nd = size(state%triangles, 1)
        if (nd /= merge(3, 6, state%element_order == 1)) return
        lower = interval(minval(state%points, dim=2))
        upper = interval(maxval(state%points, dim=2))
        if (lower(1)%lo <= 0) return
        widths = upper - lower
        if (any(widths%lo <= 0)) return
        call gs_error_friedrichs_interval(upper(1), widths(1), widths(2), cf)
        ! A = 1/R and f = s/R with s = mu0*pprime*R^2 + ffprime of degree two.
        call verified_bernstein_majorant(state%points, state%triangles(1:3, :), &
            state%element_order, estimate%flux_coefficients, state%triangles, psi, &
            state%points(1, :), numerator, 2, 1, cf, energy_bound, flux, residual, &
            checked)
        if (checked /= VERIFIED_BERNSTEIN_OK) return
        status = GS_CHECKED_OK
    contains
        subroutine numerator(cell, coefficients, ok)
            integer, intent(in) :: cell
            type(interval_t), intent(out) :: coefficients(:)
            logical, intent(out) :: ok
            type(interval_t) :: r(3)
            integer :: j
            ok = .false.
            if (size(coefficients) /= 6) return
            do j = 1, 3
                r(j) = interval(state%points(1, state%triangles(j, cell)))
            end do
            call gs_error_source_numerator_interval(interval(mu0), interval(pprime), &
                interval(ffprime), r(1), r(2), r(3), coefficients(1), &
                coefficients(2), coefficients(3), coefficients(4), coefficients(5), &
                coefficients(6))
            ok = .true.
        end subroutine numerator
    end subroutine gs_checked_equilibrated_majorant
end module kin6d_gs_quadrature_budget
