module kin6d_laboratory
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use kin6d_leaves, only: coordinate, initial_coefficients, &
                            stream_factors, stream_cached, &
                            interpolate, interpolation_weight, reconstruct, &
                            steady_source, &
                            steady_solve, steady_residual, &
                            exact_coefficients, interpolation_radius, interval_band, &
                            radius_sum, clock_step
    use fortnum_interval, only: interval_t, interval, mag
    implicit none
    private
    public :: laboratory_t, initialize, advance, sample, solve_steady, &
              certificate_radius, certificate_accepts
    type :: laboratory_t
        real(dp), allocatable :: v(:), m(:), a(:), b(:), work_a(:), work_b(:)
        real(dp), allocatable :: cv(:), sv(:)
        real(dp) :: cached_dt = 0.0_dp
        logical :: cached = .false.
        real(dp) :: amplitude = 0.5_dp
        real(dp) :: time = 0.0_dp
    end type
contains
    subroutine initialize(state, cells)
        type(laboratory_t), intent(out) :: state
        integer, intent(in) :: cells
        integer :: j
        if (cells < 2) error stop 'velocity grid requires at least two cells'
        allocate (state%v(cells + 1), state%m(cells + 1), state%a(cells + 1), &
                  state%b(cells + 1), state%work_a(cells + 1), &
                  state%work_b(cells + 1), &
                  state%cv(cells + 1), state%sv(cells + 1))
        do j = 1, cells + 1
            call coordinate(-4.0_dp, 4.0_dp, real(j - 1, dp), real(cells, dp), &
                            state%v(j))
        end do
        call initial_coefficients(state%v, state%amplitude, state%m, state%a, state%b)
    end subroutine

    subroutine advance(state, dt)
        type(laboratory_t), intent(inout) :: state
        real(dp), intent(in) :: dt
        real(dp) :: next_time
        if (.not. ieee_is_finite(dt)) error stop 'nonfinite timestep'
        if (.not. state%cached .or. &
            transfer(dt, 0_int64) /= transfer(state%cached_dt, 0_int64)) then
            call stream_factors(state%v, dt, state%cv, state%sv)
            state%cached_dt = dt
            state%cached = .true.
        end if
        call stream_cached(state%a, state%b, state%cv, state%sv, &
                           state%work_a, state%work_b)
        state%a = state%work_a
        state%b = state%work_b
        call clock_step(state%time, dt, next_time)
        state%time = next_time
    end subroutine

    subroutine sample(state, x, velocity, value)
        type(laboratory_t), intent(in) :: state
        real(dp), intent(in) :: x, velocity
        real(dp), intent(out) :: value
        real(dp) :: w, m, a, b
        integer :: j
        if (velocity < state%v(1)) error stop 'sample below velocity domain'
        if (velocity > state%v(size(state%v))) error stop 'sample above velocity domain'
        j = 1
        do while (j < size(state%v) - 1)
            if (velocity <= state%v(j + 1)) exit
            j = j + 1
        end do
        call interpolation_weight(state%v(j), state%v(j + 1), velocity, w)
        call interpolate(state%m(j), state%m(j + 1), w, m)
        call interpolate(state%a(j), state%a(j + 1), w, a)
        call interpolate(state%b(j), state%b(j + 1), w, b)
        call reconstruct(m, a, b, x, value)
    end subroutine

    subroutine solve_steady(state, nu, residual)
        type(laboratory_t), intent(inout) :: state
        real(dp), intent(in) :: nu
        real(dp), intent(out) :: residual
        real(dp) :: sc, ss, rc, rs
        integer :: j
        if (.not. ieee_is_finite(nu)) error stop 'nonfinite reservoir rate'
        if (nu <= 0.0_dp) error stop 'steady reservoir rate must be positive'
        state%time = 0.0_dp
        residual = 0.0_dp
        do j = 1, size(state%v)
            call steady_source(state%v(j), nu, state%amplitude, sc, ss)
            call steady_solve(state%v(j), nu, sc, ss, state%a(j), state%b(j))
            call steady_residual(state%a(j), state%b(j), state%v(j), nu, sc, ss, rc, rs)
            residual = max(residual, abs(rc), abs(rs))
        end do
    end subroutine

    function certificate_radius(state) result(radius)
        ! Snapshot enclosure of the whole continuous x,v reconstruction.
        ! Independently enclose exact node coefficients and the between-node
        ! interpolation remainder. No extrapolation beyond [-4,4], no Lean claim.
        type(laboratory_t), intent(in) :: state
        real(dp) :: radius, node_bound, cell_bound, vm
        type(interval_t) :: em, ea, eb, dm, da, db, unused, r, s
        integer :: j
        node_bound = 0.0_dp
        cell_bound = 0.0_dp
        do j = 1, size(state%v)
            call exact_coefficients(interval(state%v(j)), interval(state%time), &
                                    interval(state%amplitude), em, ea, eb)
            call interval_band(interval(state%m(j)), em, dm, unused)
            call interval_band(interval(state%a(j)), ea, da, unused)
            call interval_band(interval(state%b(j)), eb, db, unused)
            call radius_sum(interval(mag(dm)), interval(mag(da)), r)
            call radius_sum(r, interval(mag(db)), s)
            node_bound = max(node_bound, s%hi)
            if (j == size(state%v)) cycle
            call exact_coefficients(interval(state%v(j), state%v(j + 1)), &
                                    interval(0.0_dp), interval(state%amplitude), &
                                    em, ea, eb)
            vm = max(abs(state%v(j)), abs(state%v(j + 1)))
            call interpolation_radius(interval(vm), interval(em%hi), &
                                      interval(abs(state%time)), &
                                      interval(state%amplitude), &
                                      interval(state%v(j + 1)), interval(state%v(j)), r)
            cell_bound = max(cell_bound, r%hi)
        end do
        call radius_sum(interval(node_bound), interval(cell_bound), r)
        radius = r%hi
    end function

    logical function certificate_accepts(state, claimed_radius)
        type(laboratory_t), intent(in) :: state
        real(dp), intent(in) :: claimed_radius
        certificate_accepts = .false.
        if (.not. ieee_is_finite(claimed_radius)) return
        if (claimed_radius < 0.0_dp) return
        if (.not. ieee_is_finite(state%time)) return
        if (.not. ieee_is_finite(state%amplitude)) return
        if (state%amplitude < 0.0_dp .or. state%amplitude > 1.0_dp) return
        if (.not. allocated(state%v)) return
        if (.not. allocated(state%m)) return
        if (.not. allocated(state%a)) return
        if (.not. allocated(state%b)) return
        if (size(state%v) < 3) return
        if (any(.not. ieee_is_finite(state%v))) return
        if (size(state%m) /= size(state%v)) return
        if (size(state%a) /= size(state%v)) return
        if (size(state%b) /= size(state%v)) return
        if (any(state%v(2:) <= state%v(:size(state%v) - 1))) return
        if (any(.not. ieee_is_finite(state%m))) return
        if (any(.not. ieee_is_finite(state%a))) return
        if (any(.not. ieee_is_finite(state%b))) return
        certificate_accepts = claimed_radius >= certificate_radius(state)
    end function
end module
