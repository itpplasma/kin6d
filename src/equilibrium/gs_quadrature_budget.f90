module kin6d_gs_quadrature_budget
    !! Outward-rounded continuous energy audit for a concrete production source:
    !! -div(grad psi/R)=mu0*R*pprime+ffprime/R on an affine polygon, psi=0 on
    !! its boundary. The candidate can be affine P1 or P2. Its sampled RTp
    !! reconstruction only supplies candidate edge moments; the independent
    !! FortFEM checker constructs an exactly H(div) RT0 field and encloses its
    !! full residual and flux defect, including quadrature and arithmetic.
    !! Arbitrary callbacks, nonlinear profiles, curved cells, and pointwise
    !! error bands are outside this claim. This bound says nothing about error
    !! between the polygon and an intended curved physical domain.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortnum_interval, only: interval_t, interval, operator(+), operator(-), &
                                operator(*), operator(/)
    use fortfem_reference_scalar_intervals, only: generated_affine_p2_scalar_interval
    use fortfem_equilibrated_flux_2d, only: equilibrated_flux_estimate_t
    use fortfem_verified_rt0_majorant, only: verified_rt0_majorant
    use kin6d_grad_shafranov, only: gs_state_t
    use kin6d_grad_shafranov_leaves, only: gs_source_interval
    use kin6d_gs_error_leaves, only: gs_error_coefficient_interval, &
                                     gs_error_friedrichs_interval
    implicit none
    private
    public :: gs_checked_equilibrated_majorant
    integer, parameter, public :: GS_CHECKED_OK = 0, GS_CHECKED_UNAVAILABLE = 1
contains
    subroutine gs_checked_equilibrated_majorant(state, psi, estimate, &
                          mu0, pprime, ffprime, energy_bound, status, subdivision_level)
        type(gs_state_t), intent(in) :: state
        real(dp), intent(in) :: psi(:), mu0, pprime, ffprime
        type(equilibrated_flux_estimate_t), intent(in) :: estimate
        type(interval_t), intent(out) :: energy_bound
        integer, intent(out) :: status
        integer, intent(in), optional :: subdivision_level
        type(interval_t) :: lower(2), upper(2), widths(2), cf, flux, residual
        real(dp), allocatable :: moments(:, :)
        integer :: level, nd, e, edge, stride

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
        stride = state%element_order + 1
        if (size(estimate%flux_coefficients, 1) < 3*stride) return
        level = 1
        if (present(subdivision_level)) level = subdivision_level
        if (level < 0 .or. level > 5) return
        nd = size(state%triangles, 1)
        if (nd /= merge(3, 6, state%element_order == 1)) return
        lower = interval(minval(state%points, dim=2))
        upper = interval(maxval(state%points, dim=2))
        if (lower(1)%lo <= 0) return
        widths = upper - lower
        if (any(widths%lo <= 0)) return
        call gs_error_friedrichs_interval(upper(1), widths(1), widths(2), cf)
        allocate (moments(3, size(state%triangles, 2)))
        do e = 1, size(moments, 2)
            do edge = 1, 3
                moments(edge, e) = estimate%flux_coefficients((edge - 1)*stride + 1, e)
            end do
        end do
        call verified_rt0_majorant(state%points, state%triangles(1:3, :), moments, &
                               cf, problem, energy_bound, flux, residual, status, level)
    contains
        subroutine problem(cell, r, z, coefficient, source, gradient, ok)
            integer, intent(in) :: cell
            type(interval_t), intent(in) :: r, z
            type(interval_t), intent(out) :: coefficient, source, gradient(2)
            logical, intent(out) :: ok
            type(interval_t) :: vertices(2, 3), nodes(6), value, det
            integer :: j, next
            ok = .false.
            if (r%lo <= 0) return
            do j = 1, 3
                vertices(:, j) = interval(state%points(:, state%triangles(j, cell)))
                nodes(j) = interval(psi(state%triangles(j, cell)))
            end do
            if (state%element_order == 1) then
                ! Exact P1 embedded in P2. Interval midpoint values enclose
                ! their exact dyadic averages without changing the FE function.
                do j = 1, 3
                    next = mod(j, 3) + 1
                    nodes(3 + j) = (nodes(j) + nodes(next))*0.5_dp
                end do
            else
                nodes(4:6) = interval(psi(state%triangles(4:6, cell)))
            end if
            call generated_affine_p2_scalar_interval(r, z, &
                       vertices(1, 1), vertices(2, 1), vertices(1, 2), vertices(2, 2), &
                         vertices(1, 3), vertices(2, 3), nodes(1), nodes(2), nodes(3), &
                     nodes(4), nodes(5), nodes(6), value, gradient(1), gradient(2), det)
            if (det%lo <= 0) return
            call gs_error_coefficient_interval(r, coefficient)
            call gs_source_interval(r, interval(mu0), interval(pprime), &
                                    interval(ffprime), source)
            ok = .true.
        end subroutine problem
    end subroutine gs_checked_equilibrated_majorant
end module kin6d_gs_quadrature_budget
