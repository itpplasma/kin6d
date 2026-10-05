module kin6d_shielding_dynamics
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortnum_interval, only: interval_t, interval, is_empty, &
                                linear_interpolation_remainder
    use kin6d_shielding_dynamics_leaves, only: shielding_step, &
                        shielding_interval_length, shielding_interval_length_interval, &
     shielding_slope_interval, shielding_residual_interval, shielding_radius_interval, &
                    shielding_mode, shielding_mode_interval, shielding_error_interval, &
       shielding_add_interval, shielding_band_interval, shielding_amplitudes_interval, &
                                               shielding_curvature_interval
    implicit none
    private
    public :: shielding_state_t, init_shielding, advance_shielding, sample_shielding
    public :: shielding_coefficients, shielding_plot_padding
    type :: shielding_state_t
        private
        real(dp), allocatable :: source(:), a(:)
        type(interval_t), allocatable :: error(:)
        real(dp) :: t = 0.0_dp
        logical :: ready = .false., certified = .false.
    end type
contains
    subroutine init_shielding(state, source, ok)
        type(shielding_state_t), intent(out) :: state
        real(dp), intent(in) :: source(:)
        logical, intent(out) :: ok
        ok = .false.
        if (size(source) < 1) return
        if (.not. all(ieee_is_finite(source))) return
        allocate (state%source(size(source)), state%a(size(source)), &
                  state%error(size(source)))
        state%source = source
        state%a = 0.0_dp
        state%error = interval(0.0_dp)
        state%ready = .true.
        state%certified = .true.
        ok = .true.
    end subroutine

    subroutine advance_shielding(state, t, max_step, ok, certify)
        type(shielding_state_t), intent(inout) :: state
        real(dp), intent(in) :: t, max_step
        logical, intent(out) :: ok
        logical, intent(in), optional :: certify
        real(dp) :: next_time, h, b, k
        type(interval_t) :: hi, slope, rl, rr, residual, radius
        integer :: j
        logical :: checking
        ok = .false.
        if (.not. state%ready) return
        if (.not. ieee_is_finite(t)) return
        if (.not. ieee_is_finite(max_step)) return
        if (max_step <= 0.0_dp .or. t < state%t) return
        checking = .true.
        if (present(certify)) checking = certify
        if (checking .and. .not. state%certified) return
        state%certified = .false.
        do while (state%t < t)
            next_time = min(t, state%t + max_step)
            if (next_time <= state%t) return
            call shielding_interval_length(state%t, next_time, h)
            call shielding_interval_length_interval(interval(state%t), &
                                                    interval(next_time), hi)
            if (hi%lo <= 0.0_dp) return
            do j = 1, size(state%a)
                k = real(j, dp)
                call shielding_step(state%a(j), k, state%source(j), h, b)
                if (.not. ieee_is_finite(b)) return
                if (checking) then
                    ! Segment slope uses the exact represented endpoint times;
                    ! all floating-point candidate defects enter this residual.
                    call shielding_slope_interval(interval(state%a(j)), interval(b), &
                                                  hi, slope)
                   call shielding_residual_interval(interval(state%a(j)), interval(k), &
                                                   interval(state%source(j)), slope, rl)
                    call shielding_residual_interval(interval(b), interval(k), &
                                                   interval(state%source(j)), slope, rr)
                    residual = interval(max(abs(rl%lo), abs(rl%hi), &
                                            abs(rr%lo), abs(rr%hi)))
                    call shielding_radius_interval(state%error(j), residual, &
                                                   interval(k), hi, radius)
                    if (is_empty(radius)) return
                    if (.not. ieee_is_finite(radius%hi)) return
                    state%error(j) = interval(0.0_dp, radius%hi)
                end if
                state%a(j) = b
            end do
            state%t = next_time
        end do
        state%certified = checking
        ok = .true.
    end subroutine

    subroutine shielding_coefficients(state, source, candidate, radii, ok)
        type(shielding_state_t), intent(in) :: state
        real(dp), intent(out) :: source(:), candidate(:), radii(:)
        logical, intent(out) :: ok
        ok = .false.
        if (.not. state%ready) return
        if (size(source) /= size(state%source)) return
        if (size(candidate) /= size(state%a) .or. size(radii) /= size(state%a)) return
        source = state%source
        candidate = state%a
        radii = huge(0.0_dp)
        if (state%certified) radii = state%error%hi
        ok = .true.
    end subroutine

    subroutine sample_shielding(state, x, values, lower, upper, ok)
        type(shielding_state_t), intent(in) :: state
        real(dp), intent(in) :: x
        real(dp), intent(out) :: values(3), lower(3), upper(3)
        logical, intent(out) :: ok
        type(interval_t) :: modes(3), centers(3), errors(3), next(3), lo, hi
        real(dp) :: values_mode(3), k
        integer :: j, q
        ok = .false.
        if (.not. state%ready .or. .not. state%certified) return
        if (.not. ieee_is_finite(x)) return
        centers = interval(0.0_dp)
        errors = interval(0.0_dp)
        values = 0.0_dp
        do j = 1, size(state%a)
            k = real(j, dp)
            call shielding_mode(state%a(j), k, state%source(j), x, &
                                values_mode(1), values_mode(2), values_mode(3))
            values = values + values_mode
            call shielding_mode_interval(interval(state%a(j)), interval(k), &
                   interval(state%source(j)), interval(x), modes(1), modes(2), modes(3))
            call shielding_error_interval(interval(0.0_dp), state%error(j), &
                                          interval(k), next(1), next(2), next(3))
            do q = 1, 3
                call shielding_add_interval(centers(q), modes(q), lo)
                centers(q) = lo
                call shielding_add_interval(errors(q), next(q), lo)
                errors(q) = lo
            end do
        end do
        do q = 1, 3
            call shielding_band_interval(centers(q), errors(q), lo, hi)
            lower(q) = min(lo%lo, values(q))
            upper(q) = max(hi%hi, values(q))
            if (.not. ieee_is_finite(lower(q))) return
            if (.not. ieee_is_finite(upper(q))) return
        end do
        ok = .true.
    end subroutine

    subroutine shielding_plot_padding(state, left, right, padding, ok)
        type(shielding_state_t), intent(in) :: state
        real(dp), intent(in) :: left, right
        real(dp), intent(out) :: padding(3)
        logical, intent(out) :: ok
        type(interval_t) :: bound(3), density, potential, field, next
        integer :: j, q
        real(dp) :: k
        type(interval_t) :: amplitudes(3)
        ok = .false.
        if (.not. state%ready .or. .not. state%certified) return
        if (.not. ieee_is_finite(left) .or. .not. ieee_is_finite(right)) return
        if (right <= left) return
        bound = interval(0.0_dp)
        do j = 1, size(state%a)
            k = real(j, dp)
            call shielding_amplitudes_interval(interval(state%a(j)), &
                      interval(state%source(j)), interval(k), density, potential, field)
            amplitudes = [density, potential, field]
            do q = 1, 3
                next = amplitudes(q)
                call shielding_curvature_interval(bound(q), next, interval(k), density)
                bound(q) = density
            end do
        end do
        do q = 1, 3
            next = linear_interpolation_remainder(interval(left), interval(right), &
                                             interval(0.0_dp, max(0.0_dp, bound(q)%hi)))
            if (is_empty(next)) return
            padding(q) = next%hi
            if (.not. ieee_is_finite(padding(q))) return
        end do
        ok = .true.
    end subroutine
end module kin6d_shielding_dynamics
