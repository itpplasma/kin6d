module kin6d_gs_quadrature_budget
    !! Outward-rounded continuous energy audit for a concrete production source:
    !! -div(grad psi/R)=mu0*R*pprime+ffprime/R on an affine polygon, psi=0 on
    !! its boundary. The candidate can be affine P1, P2 or P3. Its sampled RTp
    !! reconstruction (p = element order) only supplies candidate
    !! coefficients; the independent FortFEM checker defines an exactly H(div)
    !! degree-(p+1) field with shared edge normal components and encloses its
    !! full residual and flux defect, including quadrature and arithmetic.
    !! The source is passed exactly as (mu0*pprime*R^2+ffprime)/R.
    !! Unchecked callbacks, nonlinear solution error, curved cells, and pointwise
    !! error bands are outside this claim. This bound says nothing about error
    !! between the polygon and an intended curved physical domain.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortnum_interval, only: interval_t, interval, operator(-)
    use fortfem_equilibrated_flux_2d, only: equilibrated_flux_estimate_t
    use fortfem_verified_bernstein_majorant, only: verified_bernstein_majorant, &
        VERIFIED_BERNSTEIN_OK, verified_bernstein_source, &
        verified_bernstein_source_remainder
    use kin6d_grad_shafranov, only: gs_state_t
    use kin6d_gs_error_leaves, only: gs_error_friedrichs_interval, &
                                     gs_error_source_numerator_interval
    use fortfem_verified_bernstein_majorant, only: bernstein_source_t, &
        bernstein_remainder_t
    use fortfem_verified_bernstein_majorant, only: bernstein_source_t
    implicit none
    private
    public :: gs_checked_equilibrated_majorant, gs_checked_fixed_source_majorant
    integer, parameter, public :: GS_CHECKED_OK = 0, GS_CHECKED_UNAVAILABLE = 1
    type, extends(bernstein_source_t) :: bernstein_source_adapter_t
        procedure(verified_bernstein_source), pointer, nopass :: callback => null()
    contains
        procedure :: evaluate => bernstein_source_adapter_t_evaluate
    end type
    type, extends(bernstein_remainder_t) :: bernstein_remainder_adapter_t
        procedure(verified_bernstein_source_remainder), pointer, nopass :: callback &
            => null()
    contains
        procedure :: evaluate => bernstein_remainder_adapter_t_evaluate
    end type
    interface gs_checked_fixed_source_majorant
        module procedure gs_checked_fixed_source_majorant_object, &
            gs_checked_fixed_source_majorant_procedure
    end interface
    type, extends(bernstein_source_t) :: gs_checked_equilibrated_majorant_numerator_t
        real(dp), pointer :: ffprime => null()
        real(dp), pointer :: mu0 => null()
        real(dp), pointer :: pprime => null()
        type(gs_state_t), pointer :: state => null()
    contains
        procedure :: evaluate => gs_checked_equilibrated_majorant_numerator_evaluate
    end type
contains
    subroutine gs_checked_equilibrated_majorant(state, psi, estimate, &
                                  mu0, pprime, ffprime, energy_bound, status)
        type(gs_state_t), intent(in), target :: state
        real(dp), intent(in), target :: psi(:), mu0, pprime, ffprime
        type(equilibrated_flux_estimate_t), intent(in) :: estimate
        type(interval_t), intent(out) :: energy_bound
        integer, intent(out) :: status
        type(gs_checked_equilibrated_majorant_numerator_t) :: numerator_context
        status = GS_CHECKED_UNAVAILABLE
        energy_bound = interval(0.0_dp, huge(1.0_dp))
        if (.not. all(ieee_is_finite([mu0, pprime, ffprime]))) return
        if (mu0 <= 0) return
        numerator_context%ffprime => ffprime
        numerator_context%mu0 => mu0
        numerator_context%pprime => pprime
        numerator_context%state => state
        call gs_checked_fixed_source_majorant(state, psi, estimate, numerator_context, &
            2, 1, energy_bound, status)
    end subroutine gs_checked_equilibrated_majorant
    subroutine gs_checked_fixed_source_majorant_object(state, psi, estimate, &
        numerator, &
            source_degree, weight_power, energy_bound, status, source_remainder)
        !! Certifies the frozen equation -div(grad psi/R) = (p + delta)/R^m.
        !! The caller supplies outward Bernstein coefficients of p and, for a
        !! non-polynomial source, a proved whole-cell |delta| upper bound.
        !! This does not certify a nonlinear profile-law solution error.
        type(gs_state_t), intent(in) :: state
        real(dp), intent(in) :: psi(:)
        type(equilibrated_flux_estimate_t), intent(in) :: estimate
        class(bernstein_source_t), target, intent(in) :: numerator
        integer, intent(in) :: source_degree, weight_power
        type(interval_t), intent(out) :: energy_bound
        integer, intent(out) :: status
        class(bernstein_remainder_t), target, intent(in), optional :: source_remainder
        type(interval_t) :: lower(2), upper(2), widths(2), cf, flux, residual
        integer :: nd, checked
        status = GS_CHECKED_UNAVAILABLE
        energy_bound = interval(0.0_dp, huge(1.0_dp))
        if (.not. state%ready) return
        if (.not. allocated(state%points)) return
        if (.not. allocated(state%triangles)) return
        if (.not. allocated(state%boundary)) return
        if (size(psi) /= size(state%points, 2)) return
        if (size(state%boundary) /= size(psi)) return
        if (state%element_order < 1 .or. state%element_order > 3) return
        if (allocated(state%curved_cell)) then
            if (any(state%curved_cell)) return
        end if
        if (any(psi /= 0.0_dp .and. state%boundary)) return
        if (.not. all(ieee_is_finite(psi))) return
        if (.not. allocated(estimate%flux_coefficients)) return
        if (size(estimate%flux_coefficients, 2) /= size(state%triangles, 2)) return
        nd = size(state%triangles, 1)
        if (nd /= (state%element_order + 1)*(state%element_order + 2)/2) return
        lower = interval(minval(state%points, dim=2))
        upper = interval(maxval(state%points, dim=2))
        if (lower(1)%lo <= 0) return
        widths = upper - lower
        if (any(widths%lo <= 0)) return
        call gs_error_friedrichs_interval(upper(1), widths(1), widths(2), cf)
        call verified_bernstein_majorant(state%points, state%triangles(1:3, :), &
            state%element_order, estimate%flux_coefficients, state%triangles, psi, &
            state%points(1, :), numerator, source_degree, weight_power, cf, &
            energy_bound, flux, residual, checked, source_remainder)
        if (checked /= VERIFIED_BERNSTEIN_OK) return
        status = GS_CHECKED_OK
    end subroutine gs_checked_fixed_source_majorant_object
    subroutine gs_checked_fixed_source_majorant_procedure(state, psi, estimate, &
        numerator, &
            source_degree, weight_power, energy_bound, status, source_remainder)
        !! Certifies the frozen equation -div(grad psi/R) = (p + delta)/R^m.
        !! The caller supplies outward Bernstein coefficients of p and, for a
        !! non-polynomial source, a proved whole-cell |delta| upper bound.
        !! This does not certify a nonlinear profile-law solution error.
        type(gs_state_t), intent(in) :: state
        real(dp), intent(in) :: psi(:)
        type(equilibrated_flux_estimate_t), intent(in) :: estimate
        procedure(verified_bernstein_source) :: numerator
        integer, intent(in) :: source_degree, weight_power
        type(interval_t), intent(out) :: energy_bound
        integer, intent(out) :: status
        procedure(verified_bernstein_source_remainder), optional :: source_remainder
        type(bernstein_source_adapter_t), target :: numerator_adapter
        type(bernstein_remainder_adapter_t), target :: source_remainder_adapter
        class(bernstein_remainder_t), pointer :: source_remainder_object
        nullify(source_remainder_object)
        numerator_adapter%callback => numerator
        if (present(source_remainder)) then
            source_remainder_adapter%callback => source_remainder
            source_remainder_object => source_remainder_adapter
        end if
        call gs_checked_fixed_source_majorant_object(state, &
            psi, &
            estimate, &
            numerator_adapter, &
            source_degree, &
            weight_power, &
            energy_bound, &
            status, &
            source_remainder_object)
    end subroutine
    subroutine bernstein_source_adapter_t_evaluate(self, cell, coefficients, ok)
        class(bernstein_source_adapter_t), intent(in) :: self
            !! Outward Bernstein coefficients of the polynomial numerator (or
            !! its approximation when source_remainder is supplied), using
            !! the vertex order of `triangles`. Called
            !! concurrently from OpenMP threads.
            integer, intent(in) :: cell
            type(interval_t), intent(out) :: coefficients(:)
            logical, intent(out) :: ok
        call self%callback(cell, coefficients, ok)
    end subroutine
    subroutine bernstein_remainder_adapter_t_evaluate(self, cell, radius, ok)
        class(bernstein_remainder_adapter_t), intent(in) :: self
            !! Nonnegative outward upper bound radius%hi for the absolute
            !! numerator remainder on the WHOLE cell, not sampled points.
            !! The source polynomial is supplied by verified_bernstein_source.
            !! The caller proves this contract. Called concurrently.
            integer, intent(in) :: cell
            type(interval_t), intent(out) :: radius
            logical, intent(out) :: ok
        call self%callback(cell, radius, ok)
    end subroutine
    subroutine gs_checked_equilibrated_majorant_numerator_evaluate(self, cell, &
        coefficients, ok)
        class(gs_checked_equilibrated_majorant_numerator_t), intent(in) :: self
            integer, intent(in) :: cell
            type(interval_t), intent(out) :: coefficients(:)
            logical, intent(out) :: ok
            type(interval_t) :: r(3)
            integer :: j
            ok = .false.
            if (size(coefficients) /= 6) return
            do j = 1, 3
                r(j) = interval(self%state%points(1, self%state%triangles(j, cell)))
            end do
            call gs_error_source_numerator_interval(interval(self%mu0), &
                interval(self%pprime), &
                interval(self%ffprime), r(1), r(2), r(3), coefficients(1), &
                coefficients(2), coefficients(3), coefficients(4), coefficients(5), &
                coefficients(6))
            ok = .true.
    end subroutine gs_checked_equilibrated_majorant_numerator_evaluate
end module kin6d_gs_quadrature_budget
