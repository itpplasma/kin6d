module kin6d_transport_laboratory
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortnum_interval, only: interval_t, interval, hull
    use kin6d_transport_leaves, only: transport_grid_coordinate, transport_node, &
        transport_node_jvp, transport_node_vjp, transport_scalar_add, &
        transport_interpolate, transport_mode_accumulate, &
        transport_interpolate_interval, transport_mode_accumulate_interval
    implicit none
    private
    public :: transport_state_t, initialize_transport, transport_valid, solve_transport
    public :: solve_transport_tangent, solve_transport_vjp
    public :: sample_transport, sample_transport_arrays
    public :: sample_transport_coefficients_interval, sample_transport_input_interval

    ! s(x,v) is supplied independently as a finite Fourier-x/P1-v function.
    ! v*f_x+nu*f=s on periodic x, v in [-1,1]; no velocity boundary condition
    ! is needed because this prescribed E=0 operator has no v derivative.
    type :: transport_state_t
        real(dp) :: nu = 1.0_dp
        real(dp), allocatable :: v(:)
        real(dp), allocatable :: source_cos(:, :), source_sin(:, :)
        real(dp), allocatable :: cosine(:, :), sine(:, :)
    end type transport_state_t

contains

    subroutine initialize_transport(state, cells, modes, nu, ok)
        type(transport_state_t), intent(out) :: state
        integer, intent(in) :: cells, modes
        real(dp), intent(in) :: nu
        logical, intent(out), optional :: ok
        integer :: j

        if (present(ok)) ok = .false.
        if (cells < 1 .or. modes < 0) return
        if (.not. ieee_is_finite(nu)) return
        if (nu <= 0.0_dp) return
        state%nu = nu
        allocate (state%v(cells + 1))
        allocate (state%source_cos(0:modes, cells + 1), &
                  state%source_sin(0:modes, cells + 1), &
                  state%cosine(0:modes, cells + 1), state%sine(0:modes, cells + 1))
        do j = 1, cells + 1
            call transport_grid_coordinate(real(j - 1, dp), real(cells, dp), state%v(j))
        end do
        state%source_cos = 0.0_dp
        state%source_sin = 0.0_dp
        state%cosine = 0.0_dp
        state%sine = 0.0_dp
        if (present(ok)) ok = transport_valid(state)
    end subroutine initialize_transport

    logical function transport_valid(state) result(ok)
        type(transport_state_t), intent(in) :: state
        integer :: j

        ok = .false.
        if (.not. allocated(state%v)) return
        if (.not. allocated(state%source_cos)) return
        if (.not. allocated(state%source_sin)) return
        if (.not. allocated(state%cosine)) return
        if (.not. allocated(state%sine)) return
        if (size(state%v) < 2) return
        if (lbound(state%v, 1) /= 1) return
        if (lbound(state%source_cos, 2) /= 1 .or. &
            lbound(state%source_sin, 2) /= 1) return
        if (lbound(state%cosine, 2) /= 1 .or. lbound(state%sine, 2) /= 1) return
        if (lbound(state%source_cos, 1) /= 0 .or. &
            lbound(state%source_sin, 1) /= 0) return
        if (lbound(state%cosine, 1) /= 0 .or. lbound(state%sine, 1) /= 0) return
        if (size(state%source_cos, 1) < 1) return
        if (any(shape(state%source_cos) /= shape(state%source_sin))) return
        if (any(shape(state%source_cos) /= shape(state%cosine))) return
        if (any(shape(state%source_cos) /= shape(state%sine))) return
        if (size(state%source_cos, 2) /= size(state%v)) return
        if (.not. ieee_is_finite(state%nu)) return
        if (state%nu <= 0.0_dp) return
        if (.not. all(ieee_is_finite(state%v))) return
        if (.not. all(ieee_is_finite(state%source_cos))) return
        if (.not. all(ieee_is_finite(state%source_sin))) return
        if (.not. all(ieee_is_finite(state%cosine))) return
        if (.not. all(ieee_is_finite(state%sine))) return
        if (any(abs(state%source_sin(0, :)) > 0.0_dp)) return
        if (any(abs(state%sine(0, :)) > 0.0_dp)) return
        if (state%v(1) < -1.0_dp .or. state%v(1) > -1.0_dp) return
        if (state%v(size(state%v)) < 1.0_dp .or. state%v(size(state%v)) > 1.0_dp) return
        do j = 1, size(state%v) - 1
            if (state%v(j + 1) <= state%v(j)) return
        end do
        ok = .true.
    end function transport_valid

    subroutine solve_transport(state, ok)
        type(transport_state_t), intent(inout) :: state
        logical, intent(out) :: ok
        integer :: j, k

        ok = .false.
        if (.not. transport_valid(state)) return
        do j = 1, size(state%v)
            do k = 0, size(state%cosine, 1) - 1
                call transport_node(state%v(j), state%nu, real(k, dp), &
                                    state%source_cos(k, j), state%source_sin(k, j), &
                                    state%cosine(k, j), state%sine(k, j))
            end do
        end do
        ok = transport_valid(state)
    end subroutine solve_transport

    subroutine solve_transport_tangent(state, dnu, dsource_cos, dsource_sin, &
                                       dcos, dsin, ok)
        ! The generated inverse JVP solves the differentiated node equation.
        ! Buffers are caller-owned; no solver trajectory or Jacobian is stored.
        type(transport_state_t), intent(in) :: state
        real(dp), intent(in) :: dnu, dsource_cos(0:, :), dsource_sin(0:, :)
        real(dp), intent(out) :: dcos(0:, :), dsin(0:, :)
        logical, intent(out) :: ok
        integer :: j, k

        ok = .false.
        if (.not. transport_valid(state)) return
        if (.not. ieee_is_finite(dnu)) return
        if (any(shape(dsource_cos) /= shape(state%cosine))) return
        if (any(shape(dsource_sin) /= shape(state%cosine))) return
        if (any(shape(dcos) /= shape(state%cosine))) return
        if (any(shape(dsin) /= shape(state%cosine))) return
        if (.not. all(ieee_is_finite(dsource_cos))) return
        if (.not. all(ieee_is_finite(dsource_sin))) return
        if (any(abs(dsource_sin(0, :)) > 0.0_dp)) return
        do j = 1, size(state%v)
            do k = 0, size(state%cosine, 1) - 1
                call transport_node_jvp(state%v(j), state%nu, real(k, dp), &
                                       state%source_cos(k, j), state%source_sin(k, j), &
                            0.0_dp, dnu, 0.0_dp, dsource_cos(k, j), dsource_sin(k, j), &
                                        dcos(k, j), dsin(k, j))
            end do
        end do
        if (.not. all(ieee_is_finite(dcos))) return
        if (.not. all(ieee_is_finite(dsin))) return
        ok = .true.
    end subroutine solve_transport_tangent

    subroutine solve_transport_vjp(state, bar_cos, bar_sin, bar_nu, &
                                   bar_source_cos, bar_source_sin, ok)
        ! Raw Euclidean nodal coefficient pairing plus one global nu scalar.
        ! This is an implicit node-solve transpose, not a continuous adjoint.
        type(transport_state_t), intent(in) :: state
        real(dp), intent(in) :: bar_cos(0:, :), bar_sin(0:, :)
        real(dp), intent(out) :: bar_nu, bar_source_cos(0:, :), bar_source_sin(0:, :)
        logical, intent(out) :: ok
        real(dp) :: gv, gnu, gk, updated
        integer :: j, k

        ok = .false.
        bar_nu = 0.0_dp
        if (.not. transport_valid(state)) return
        if (any(shape(bar_cos) /= shape(state%cosine))) return
        if (any(shape(bar_sin) /= shape(state%cosine))) return
        if (any(shape(bar_source_cos) /= shape(state%cosine))) return
        if (any(shape(bar_source_sin) /= shape(state%cosine))) return
        if (.not. all(ieee_is_finite(bar_cos))) return
        if (.not. all(ieee_is_finite(bar_sin))) return
        do j = 1, size(state%v)
            do k = 0, size(state%cosine, 1) - 1
                call transport_node_vjp(state%v(j), state%nu, real(k, dp), &
                                       state%source_cos(k, j), state%source_sin(k, j), &
                                        bar_cos(k, j), bar_sin(k, j), gv, gnu, gk, &
                                        bar_source_cos(k, j), bar_source_sin(k, j))
                call transport_scalar_add(bar_nu, gnu, updated)
                bar_nu = updated
            end do
        end do
        ! Canonical dual on the physical mean-sine constraint subspace.
        bar_source_sin(0, :) = 0.0_dp
        if (.not. ieee_is_finite(bar_nu)) return
        if (.not. all(ieee_is_finite(bar_source_cos))) return
        if (.not. all(ieee_is_finite(bar_source_sin))) return
        ok = .true.
    end subroutine solve_transport_vjp

    subroutine sample_transport(state, x, velocity, value, ok)
        type(transport_state_t), intent(in) :: state
        real(dp), intent(in) :: x, velocity
        real(dp), intent(out) :: value
        logical, intent(out) :: ok

        ok = .false.
        value = 0.0_dp
        if (.not. allocated(state%v)) return
        if (.not. allocated(state%cosine)) return
        if (.not. allocated(state%sine)) return
        call sample_transport_arrays(state%v, state%cosine, state%sine, &
                                     x, velocity, value, ok)
    end subroutine sample_transport

    subroutine sample_transport_arrays(grid, cosine, sine, x, velocity, value, ok)
        ! Fast evaluation requires a globally validated grid/coefficient packet.
        ! ok checks this evaluation, not the validity of unvisited cells.
        real(dp), intent(in) :: grid(:), cosine(0:, :), sine(0:, :), x, velocity
        real(dp), intent(out) :: value
        logical, intent(out) :: ok
        real(dp) :: a, b, updated
        integer :: j, k

        value = 0.0_dp
        ok = .false.
        if (.not. ieee_is_finite(x)) return
        if (size(cosine, 1) < 1) return
        if (size(cosine, 2) /= size(grid)) return
        if (any(shape(cosine) /= shape(sine))) return
        call locate_cell(grid, velocity, j, ok)
        if (.not. ok) return
        ok = .false.
        if (.not. all(ieee_is_finite(cosine(:, j:j + 1)))) return
        if (.not. all(ieee_is_finite(sine(:, j:j + 1)))) return
        if (any(abs(sine(0, j:j + 1)) > 0.0_dp)) return
        do k = 0, size(cosine, 1) - 1
            call transport_interpolate(velocity, grid(j), grid(j + 1), &
                                       cosine(k, j), cosine(k, j + 1), a)
            call transport_interpolate(velocity, grid(j), grid(j + 1), &
                                       sine(k, j), sine(k, j + 1), b)
            call transport_mode_accumulate(x, real(k, dp), a, b, value, updated)
            value = updated
        end do
        ok = ieee_is_finite(value)
    end subroutine sample_transport_arrays

    subroutine locate_cell(grid, velocity, j, ok)
        real(dp), intent(in) :: grid(:), velocity
        integer, intent(out) :: j
        logical, intent(out) :: ok
        integer :: left, right, middle

        ok = .false.
        j = 1
        if (size(grid) < 2) return
        if (.not. ieee_is_finite(velocity)) return
        if (velocity < -1.0_dp .or. velocity > 1.0_dp) return
        if (grid(1) < -1.0_dp .or. grid(1) > -1.0_dp) return
        if (grid(size(grid)) < 1.0_dp .or. grid(size(grid)) > 1.0_dp) return
        left = 1
        right = size(grid)
        do while (right - left > 1)
            middle = left + (right - left)/2
            if (velocity < grid(middle)) then
                right = middle
            else
                left = middle
            end if
        end do
        j = left
        if (.not. ieee_is_finite(grid(j))) return
        if (.not. ieee_is_finite(grid(j + 1))) return
        if (grid(j + 1) <= grid(j)) return
        if (velocity < grid(j) .or. velocity > grid(j + 1)) return
        ok = .true.
    end subroutine locate_cell

    subroutine sample_transport_input_interval(state, x, velocity, value, ok)
        type(transport_state_t), intent(in) :: state
        type(interval_t), intent(in) :: x, velocity
        type(interval_t), intent(out) :: value
        logical, intent(out) :: ok

        ok = .false.
        value = interval(0.0_dp)
        if (.not. allocated(state%source_cos)) return
        if (.not. allocated(state%source_sin)) return
        call sample_transport_coefficients_interval(state, state%source_cos, &
                                               state%source_sin, x, velocity, value, ok)
    end subroutine sample_transport_input_interval

    subroutine sample_transport_coefficients_interval(state, cosine, sine, &
                                                      x, velocity, value, ok)
        ! Local interval reconstruction, including every touched velocity cell.
        ! No per-call allocations or global packet rescans. Certified wrappers
        ! prepare/validate the complete immutable packet separately.
        type(transport_state_t), intent(in) :: state
        real(dp), intent(in) :: cosine(0:, :), sine(0:, :)
        type(interval_t), intent(in) :: x, velocity
        type(interval_t), intent(out) :: value
        logical, intent(out) :: ok
        type(interval_t) :: local_v, a, b, accumulator, updated
        integer :: first, last, j, k
        logical :: located

        ok = .false.
        value = interval(0.0_dp)
        if (.not. allocated(state%v)) return
        if (size(cosine, 1) < 1) return
        if (size(cosine, 2) /= size(state%v)) return
        if (any(shape(cosine) /= shape(sine))) return
        if (.not. ieee_is_finite(x%lo) .or. .not. ieee_is_finite(x%hi)) return
        if (.not. ieee_is_finite(velocity%lo) .or. &
            .not. ieee_is_finite(velocity%hi)) return
        if (x%lo > x%hi .or. velocity%lo > velocity%hi) return
        call locate_cell(state%v, velocity%lo, first, located)
        if (.not. located) return
        call locate_cell(state%v, velocity%hi, last, located)
        if (.not. located) return
        if (last < first) return
        do j = first, last
            if (.not. ieee_is_finite(state%v(j))) return
            if (.not. ieee_is_finite(state%v(j + 1))) return
            if (state%v(j + 1) <= state%v(j)) return
            if (.not. all(ieee_is_finite(cosine(:, j:j + 1)))) return
            if (.not. all(ieee_is_finite(sine(:, j:j + 1)))) return
            if (any(abs(sine(0, j:j + 1)) > 0.0_dp)) return
            local_v = interval(max(velocity%lo, state%v(j)), &
                               min(velocity%hi, state%v(j + 1)))
            accumulator = interval(0.0_dp)
            do k = 0, size(cosine, 1) - 1
                call transport_interpolate_interval(local_v, interval(state%v(j)), &
                                     interval(state%v(j + 1)), interval(cosine(k, j)), &
                                                    interval(cosine(k, j + 1)), a)
                call transport_interpolate_interval(local_v, interval(state%v(j)), &
                                       interval(state%v(j + 1)), interval(sine(k, j)), &
                                                    interval(sine(k, j + 1)), b)
                call transport_mode_accumulate_interval(x, interval(real(k, dp)), &
                                                        a, b, accumulator, updated)
                accumulator = updated
            end do
            if (j == first) then
                value = accumulator
            else
                value = hull(value, accumulator)
            end if
        end do
        ok = ieee_is_finite(value%lo) .and. ieee_is_finite(value%hi)
    end subroutine sample_transport_coefficients_interval

end module kin6d_transport_laboratory
