module kin6d_fp_laboratory
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortnum_interval, only: interval_t, interval, mag
    use kin6d_fp_leaves
    implicit none
    private
    public :: fp_state_t, initialize_fp, solve_fp, sample_fp, fp_certificate_radius, &
              fp_certificate_accepts, fp_sample_band, fp_reference_band, fp_mass_interval, &
              fp_normalization_interval
    type :: fp_state_t
        real(dp), allocatable :: v(:), f(:)
        real(dp) :: theta = 1.0_dp
        real(dp) :: mass = 1.0_dp
    end type
contains
    subroutine initialize_fp(state, cells, theta, mass)
        type(fp_state_t), intent(out) :: state
        integer, intent(in) :: cells
        real(dp), intent(in), optional :: theta, mass
        integer :: j
        if (cells < 2) error stop 'FP grid needs at least two cells'
        if (present(theta)) state%theta = theta
        if (present(mass)) state%mass = mass
        allocate (state%v(cells + 1), state%f(cells + 1))
        do j = 1, cells + 1
            call fp_coordinate(-2.0_dp, 2.0_dp, real(j - 1, dp), &
                               real(cells, dp), state%v(j))
        end do
        state%f = 1.0_dp
        if (.not. valid_fp(state)) error stop 'invalid FP heat bath or mass'
    end subroutine

    subroutine solve_fp(state, residual)
        ! Direct zero-flux recurrence, then the prescribed discrete mass.
        ! O(cells) work and two state vectors; no dense matrix or trajectory.
        type(fp_state_t), intent(inout) :: state
        real(dp), intent(out) :: residual
        real(dp) :: midpoint, width, area, total, updated, flux
        integer :: j
        if (.not. valid_fp(state)) error stop 'invalid FP state'
        state%f(1) = 1.0_dp
        total = 0.0_dp
        do j = 1, size(state%v) - 1
            call fp_cell(state%v(j), state%v(j + 1), midpoint, width)
            call fp_step(midpoint, state%theta, width, state%f(j), state%f(j + 1))
            if (.not. ieee_is_finite(state%f(j + 1))) error stop 'invalid FP recurrence'
            if (state%f(j + 1) <= 0.0_dp) error stop 'FP positivity needs finer grid'
            call fp_area(width, state%f(j), state%f(j + 1), area)
            call fp_add(total, area, updated)
            total = updated
        end do
        do j = 1, size(state%f)
            call fp_scale(state%f(j), state%mass, total, updated)
            state%f(j) = updated
        end do
        residual = 0.0_dp
        do j = 1, size(state%v) - 1
            call fp_cell(state%v(j), state%v(j + 1), midpoint, width)
            call fp_flux(midpoint, state%theta, width, state%f(j), state%f(j + 1), flux)
            residual = max(residual, abs(flux))
        end do
    end subroutine

    logical function valid_fp(state)
        type(fp_state_t), intent(in) :: state
        integer :: j
        valid_fp = .false.
        if (.not. allocated(state%v)) return
        if (.not. allocated(state%f)) return
        if (size(state%v) < 3) return
        if (size(state%f) /= size(state%v)) return
        if (.not. ieee_is_finite(state%theta)) return
        if (.not. ieee_is_finite(state%mass)) return
        if (state%theta <= 0.0_dp) return
        if (state%mass <= 0.0_dp) return
        if (.not. all(ieee_is_finite(state%v))) return
        if (.not. all(ieee_is_finite(state%f))) return
        if (any(state%f <= 0.0_dp)) return
        if (state%v(1) < -2.0_dp .or. state%v(1) > -2.0_dp) return
        if (state%v(size(state%v)) < 2.0_dp .or. &
            state%v(size(state%v)) > 2.0_dp) return
        do j = 1, size(state%v) - 1
            if (state%v(j + 1) <= state%v(j)) return
        end do
        valid_fp = .true.
    end function

    subroutine normalizer(state, z)
        ! Midpoint integral with rigorous cellwise second-derivative remainder.
        ! Z is the exact finite-domain Gaussian integral, not a trapezoidal value.
        type(fp_state_t), intent(in) :: state
        type(interval_t), intent(out) :: z
        type(interval_t) :: midpoint, width, lower, upper, term, updated
        real(dp) :: vm
        integer :: j
        z = interval(0.0_dp)
        do j = 1, size(state%v) - 1
            call fp_interval_cell(interval(state%v(j)), interval(state%v(j + 1)), &
                                  midpoint, width)
            vm = max(abs(state%v(j)), abs(state%v(j + 1)))
            call fp_normalizer_cell(midpoint, interval(state%theta), width, &
                                    interval(vm), lower, upper)
            term = interval(lower%lo, upper%hi)
            call fp_interval_add(z, term, updated)
            z = updated
        end do
    end subroutine

    function fp_certificate_radius(state) result(radius)
        ! L-infinity enclosure of the piecewise-linear candidate on [-2,2].
        ! Exact ODE plus mass proves unique normalized Gaussian equilibrium.
        ! Node evaluation, exact Z enclosure, and interpolation remainder are
        ! separate from the discrete flux diagnostic. This is not a Lean proof.
        type(fp_state_t), intent(in) :: state
        real(dp) :: radius, node_bound, cell_bound, vm
        type(interval_t) :: z, exact, error, midpoint, width, remainder, sum_bound
        integer :: j
        radius = huge(0.0_dp)
        if (.not. valid_fp(state)) return
        call normalizer(state, z)
        if (z%lo <= 0.0_dp) return
        if (.not. ieee_is_finite(z%lo)) return
        if (.not. ieee_is_finite(z%hi)) return
        node_bound = 0.0_dp; cell_bound = 0.0_dp
        do j = 1, size(state%v)
            call fp_reference(interval(state%v(j)), interval(state%theta), &
                              interval(state%mass), z, exact)
            call fp_difference(interval(state%f(j)), exact, error)
            if (.not. ieee_is_finite(error%lo)) return
            if (.not. ieee_is_finite(error%hi)) return
            node_bound = max(node_bound, mag(error))
            if (j == size(state%v)) cycle
            call fp_interval_cell(interval(state%v(j)), interval(state%v(j + 1)), &
                                  midpoint, width)
            vm = max(abs(state%v(j)), abs(state%v(j + 1)))
            call fp_interpolation_error(interval(state%theta), interval(state%mass), &
                                        z, width, interval(vm), remainder)
            if (.not. ieee_is_finite(remainder%hi)) return
            cell_bound = max(cell_bound, remainder%hi)
        end do
        call fp_interval_add(interval(node_bound), interval(cell_bound), sum_bound)
        radius = sum_bound%hi
    end function

    logical function fp_certificate_accepts(state, radius)
        type(fp_state_t), intent(in) :: state
        real(dp), intent(in) :: radius
        real(dp) :: required
        fp_certificate_accepts = .false.
        if (.not. valid_fp(state)) return
        if (.not. ieee_is_finite(radius)) return
        if (radius < 0.0_dp) return
        required = fp_certificate_radius(state)
        if (.not. ieee_is_finite(required)) return
        if (required >= huge(0.0_dp)) return
        fp_certificate_accepts = radius >= required
    end function

    integer function locate(state, velocity) result(j)
        type(fp_state_t), intent(in) :: state
        real(dp), intent(in) :: velocity
        if (.not. valid_fp(state)) error stop 'invalid FP sampling state'
        if (.not. ieee_is_finite(velocity)) error stop 'invalid FP sample velocity'
        if (velocity < -2.0_dp .or. velocity > 2.0_dp) error stop 'FP sample outside domain'
        j = 1
        do while (j < size(state%v) - 1)
            if (velocity <= state%v(j + 1)) exit
            j = j + 1
        end do
    end function

    subroutine sample_fp(state, velocity, value)
        type(fp_state_t), intent(in) :: state
        real(dp), intent(in) :: velocity
        real(dp), intent(out) :: value
        real(dp) :: weight
        integer :: j
        j = locate(state, velocity)
        call fp_weight(state%v(j), state%v(j + 1), velocity, weight)
        call fp_interpolate(state%f(j), state%f(j + 1), weight, value)
    end subroutine

    subroutine fp_sample_band(state, velocity, lower, upper)
        type(fp_state_t), intent(in) :: state
        real(dp), intent(in) :: velocity
        real(dp), intent(out) :: lower, upper
        type(interval_t) :: weight, value, lo, hi
        real(dp) :: radius
        integer :: j
        j = locate(state, velocity)
        radius = fp_certificate_radius(state)
        if (.not. fp_certificate_accepts(state, radius)) error stop 'FP enclosure unavailable'
        call fp_interval_weight(interval(state%v(j)), interval(state%v(j + 1)), &
                                interval(velocity), weight)
        call fp_interval_interpolate(interval(state%f(j)), interval(state%f(j + 1)), &
                                     weight, value)
        call fp_band(value, interval(radius), lo, hi)
        lower = lo%lo; upper = hi%hi
    end subroutine

    subroutine fp_reference_band(state, velocity, lower, upper)
        type(fp_state_t), intent(in) :: state
        real(dp), intent(in) :: velocity
        real(dp), intent(out) :: lower, upper
        type(interval_t) :: z, exact
        integer :: j
        j = locate(state, velocity)
        call normalizer(state, z)
        if (z%lo <= 0.0_dp) error stop 'FP reference normalization unavailable'
        if (.not. ieee_is_finite(z%lo)) error stop 'FP reference normalization overflow'
        if (.not. ieee_is_finite(z%hi)) error stop 'FP reference normalization overflow'
        call fp_reference(interval(velocity), interval(state%theta), &
                          interval(state%mass), z, exact)
        lower = exact%lo; upper = exact%hi
    end subroutine

    subroutine fp_mass_interval(state, lower, upper)
        type(fp_state_t), intent(in) :: state
        real(dp), intent(out) :: lower, upper
        type(interval_t) :: midpoint, width, area, total, updated
        integer :: j
        if (.not. valid_fp(state)) error stop 'invalid FP mass state'
        total = interval(0.0_dp)
        do j = 1, size(state%v) - 1
            call fp_interval_cell(interval(state%v(j)), interval(state%v(j + 1)), &
                                  midpoint, width)
            call fp_interval_area(width, interval(state%f(j)), interval(state%f(j + 1)), area)
            call fp_interval_add(total, area, updated)
            total = updated
        end do
        lower = total%lo; upper = total%hi
    end subroutine

    subroutine fp_normalization_interval(state, lower, upper)
        type(fp_state_t), intent(in) :: state
        real(dp), intent(out) :: lower, upper
        type(interval_t) :: z
        if (.not. valid_fp(state)) error stop 'invalid FP reference state'
        call normalizer(state, z)
        if (z%lo <= 0.0_dp) error stop 'FP normalization enclosure unavailable'
        if (.not. ieee_is_finite(z%lo)) error stop 'FP normalization overflow'
        if (.not. ieee_is_finite(z%hi)) error stop 'FP normalization overflow'
        lower = z%lo; upper = z%hi
    end subroutine
end module
