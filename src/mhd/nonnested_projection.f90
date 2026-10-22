module kin6d_nonnested_projection
    ! I5b diagnostic only. Component membership and positive quadrature of
    ! invariant measure are declared inputs, never inferred inside Newton.
    ! Closed lines need dl/|B|, mixing volumes need volume weights. Finite
    ! traces do not prove ergodicity. Open lines require heat boundary data
    ! and are deliberately unsupported here.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    private
    public :: nn_projection_t, nn_projection_report_t, nn_projection_build
    integer, parameter, public :: NN_CLOSED_LINE = 1, NN_INVARIANT_SURFACE = 2
    integer, parameter, public :: NN_MIXING_VOLUME = 3

    type :: nn_projection_report_t
        real(dp) :: target_mismatch = 0, mean_defect = 0
        real(dp) :: trace_change = -1
        character(32) :: claim = 'declared-component-diagnostic'
        character(32) :: derivative_status = 'fixed-components-only'
        character(32) :: bound_status = 'unavailable'
    end type
    type :: nn_projection_t
        private
        integer, allocatable :: component(:), kind(:)
        real(dp), allocatable :: weight(:), mass(:)
    contains
        procedure :: apply => projection_apply
    end type
contains
    subroutine nn_projection_build(component, kind, weight, projection, status)
        integer, intent(in) :: component(:), kind(:)
        real(dp), intent(in) :: weight(:)
        type(nn_projection_t), intent(out) :: projection
        integer, intent(out) :: status
        integer :: i
        status = 1
        if (size(component) == 0 .or. size(kind) == 0) return
        if (size(component) /= size(weight)) return
        if (any(component < 1) .or. any(component > size(kind))) return
        if (any(kind < NN_CLOSED_LINE) .or. any(kind > NN_MIXING_VOLUME)) return
        if (.not. all(ieee_is_finite(weight))) return
        if (any(weight <= 0)) return
        projection%component = component; projection%kind = kind
        projection%weight = weight
        allocate (projection%mass(size(kind)))
        projection%mass = 0
        do i = 1, size(weight)
            projection%mass(component(i)) = projection%mass(component(i)) + weight(i)
        end do
        if (any(projection%mass <= 0)) return
        if (.not. all(ieee_is_finite(projection%mass))) return
        status = 0
    end subroutine

    subroutine projection_apply(self, target, projected, report, status, previous)
        class(nn_projection_t), intent(in) :: self
        real(dp), intent(in) :: target(:)
        real(dp), intent(out) :: projected(:)
        type(nn_projection_report_t), intent(out) :: report
        integer, intent(out) :: status
        real(dp), intent(in), optional :: previous(:)
        real(dp), allocatable :: means(:), after(:)
        integer :: i, c
        status = 1; projected = 0
        if (.not. allocated(self%mass)) return
        if (size(target) /= size(self%component)) return
        if (size(projected) /= size(target)) return
        if (.not. all(ieee_is_finite(target))) return
        if (present(previous)) then
            if (size(previous) /= size(target)) return
            if (.not. all(ieee_is_finite(previous))) return
        end if
        allocate (means(size(self%mass)), after(size(self%mass)))
        means = 0; after = 0
        !$omp parallel do default(shared) private(i,c) reduction(+:means) &
        !$omp& schedule(static)
        do i = 1, size(target)
            c = self%component(i)
            means(c) = means(c) + self%weight(i)*target(i)
        end do
        !$omp end parallel do
        means = means/self%mass
        !$omp parallel do default(shared) private(i,c) reduction(+:after) &
        !$omp& schedule(static)
        do i = 1, size(target)
            c = self%component(i)
            projected(i) = means(c)
            after(c) = after(c) + self%weight(i)*projected(i)
        end do
        !$omp end parallel do
        report%mean_defect = maxval(abs(after/self%mass - means))
        report%target_mismatch = sqrt(sum(self%weight*(projected - target)**2) &
                                      /sum(self%mass))
        if (present(previous)) report%trace_change = &
            sqrt(sum(self%weight*(projected - previous)**2)/sum(self%mass))
        if (.not. all(ieee_is_finite(projected))) return
        status = 0
    end subroutine
end module kin6d_nonnested_projection
