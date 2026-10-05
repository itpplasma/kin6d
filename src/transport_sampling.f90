module kin6d_transport_sampling
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortnum_interval, only: interval_t, interval, hull
    use kin6d_transport_laboratory, only: transport_state_t, locate_cell
    use kin6d_transport_certificate, only: transport_certificate_t, &
        transport_observable_band
    use kin6d_transport_leaves, only: transport_basis_interval, &
        transport_weighted_mode_accumulate_interval, transport_interpolate_interval, &
        transport_scalar_add_interval
    implicit none
    private
    public :: transport_basis_t, transport_sampling_workspace_t
    public :: initialize_transport_sampling, prepare_transport_basis
    public :: sample_transport_cached, sample_transport_input_cached
    public :: sample_transport_coefficients_cached, certify_transport_cached
    public :: certify_transport_tangent_cached, transport_sampling_bytes

    ! Only prepare_transport_basis can create the captured rigorous factors.
    ! State/certificate/tangent packets remain immutable while being sampled.
    type :: transport_basis_t
        private
        integer :: modes = -1
        logical :: prepared = .false.
        type(interval_t), allocatable :: x(:), cosine(:, :), sine(:, :)
    end type
    type :: transport_sampling_workspace_t
        private
        type(interval_t), allocatable :: cell_range(:), candidate(:)
    end type
contains
    subroutine initialize_transport_sampling(basis, work, modes, cells, ok)
        type(transport_basis_t), intent(out) :: basis
        type(transport_sampling_workspace_t), intent(out) :: work
        integer, intent(in) :: modes, cells
        logical, intent(out) :: ok

        ok = .false.
        if (modes < 0 .or. cells < 1) return
        allocate (basis%x(cells), basis%cosine(cells, 0:modes), &
                  basis%sine(cells, 0:modes))
        allocate (work%cell_range(cells), work%candidate(cells))
        basis%modes = modes
        ok = .true.
    end subroutine

    subroutine prepare_transport_basis(basis, x_cells, ok)
        type(transport_basis_t), intent(inout) :: basis
        type(interval_t), intent(in) :: x_cells(:)
        logical, intent(out) :: ok
        integer :: i, k

        ok = .false.
        basis%prepared = .false.
        if (.not. allocated(basis%x)) return
        if (size(x_cells) /= size(basis%x)) return
        do i = 1, size(x_cells)
            if (.not. finite_interval(x_cells(i))) return
        end do
        basis%x = x_cells
        basis%cosine(:, 0) = interval(1.0_dp)
        basis%sine(:, 0) = interval(0.0_dp)
        do k = 1, basis%modes
            do i = 1, size(x_cells)
                call transport_basis_interval(x_cells(i), interval(real(k, dp)), &
                    basis%sine(i, k), basis%cosine(i, k))
                if (.not. finite_interval(basis%sine(i, k))) return
                if (.not. finite_interval(basis%cosine(i, k))) return
            end do
        end do
        basis%prepared = .true.
        ok = .true.
    end subroutine

    logical function finite_interval(value) result(ok)
        type(interval_t), intent(in) :: value

        ok = .false.
        if (.not. ieee_is_finite(value%lo)) return
        if (.not. ieee_is_finite(value%hi)) return
        ok = value%lo <= value%hi
    end function

    logical function matching_geometry(basis, work, x_cells) result(ok)
        type(transport_basis_t), intent(in) :: basis
        type(transport_sampling_workspace_t), intent(in) :: work
        type(interval_t), intent(in) :: x_cells(:)
        integer :: i

        ok = .false.
        if (.not. basis%prepared) return
        if (.not. allocated(basis%x)) return
        if (.not. allocated(work%candidate)) return
        if (.not. allocated(work%cell_range)) return
        if (size(x_cells) /= size(basis%x)) return
        if (size(work%candidate) /= size(x_cells)) return
        if (size(work%cell_range) /= size(x_cells)) return
        do i = 1, size(x_cells)
            if (transfer(x_cells(i)%lo, 0_int64) /= &
                transfer(basis%x(i)%lo, 0_int64)) return
            if (transfer(x_cells(i)%hi, 0_int64) /= &
                transfer(basis%x(i)%hi, 0_int64)) return
        end do
        ok = .true.
    end function

    subroutine reconstruct(state, cosine, sine, basis, x_cells, velocity, work, ok)
        type(transport_state_t), intent(in) :: state
        real(dp), intent(in) :: cosine(0:, :), sine(0:, :)
        type(transport_basis_t), intent(in) :: basis
        type(interval_t), intent(in) :: x_cells(:), velocity
        type(transport_sampling_workspace_t), intent(inout) :: work
        logical, intent(out) :: ok
        type(interval_t) :: local_v, a, b, updated
        integer :: first, last, j, k, i
        logical :: located

        ok = .false.
        if (.not. matching_geometry(basis, work, x_cells)) return
        if (.not. allocated(state%v)) return
        if (size(cosine, 1) /= basis%modes + 1) return
        if (size(cosine, 2) /= size(state%v)) return
        if (any(shape(cosine) /= shape(sine))) return
        if (.not. finite_interval(velocity)) return
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
            call transport_interpolate_interval(local_v, interval(state%v(j)), &
                interval(state%v(j + 1)), interval(cosine(0, j)), &
                interval(cosine(0, j + 1)), a)
            ! Mode zero has no trigonometry and no sine contribution.
            call transport_scalar_add_interval(interval(0.0_dp), a, updated)
            work%cell_range = updated
            do k = 1, basis%modes
                call transport_interpolate_interval(local_v, interval(state%v(j)), &
                    interval(state%v(j + 1)), interval(cosine(k, j)), &
                    interval(cosine(k, j + 1)), a)
                call transport_interpolate_interval(local_v, interval(state%v(j)), &
                    interval(state%v(j + 1)), interval(sine(k, j)), &
                    interval(sine(k, j + 1)), b)
                do i = 1, size(x_cells)
                    call transport_weighted_mode_accumulate_interval( &
                        basis%sine(i, k), basis%cosine(i, k), a, b, &
                        work%cell_range(i), &
                        updated)
                    work%cell_range(i) = updated
                end do
            end do
            do i = 1, size(x_cells)
                if (.not. finite_interval(work%cell_range(i))) return
                if (j == first) then
                    work%candidate(i) = work%cell_range(i)
                else
                    work%candidate(i) = hull(work%candidate(i), work%cell_range(i))
                end if
            end do
        end do
        ok = .true.
    end subroutine

    subroutine sample_transport_coefficients_cached(state, cosine, sine, basis, &
        x_cells, velocity, work, values, ok)
        type(transport_state_t), intent(in) :: state
        real(dp), intent(in) :: cosine(0:, :), sine(0:, :)
        type(transport_basis_t), intent(in) :: basis
        type(interval_t), intent(in) :: x_cells(:), velocity
        type(transport_sampling_workspace_t), intent(inout) :: work
        type(interval_t), intent(out) :: values(:)
        logical, intent(out) :: ok

        values = interval(0.0_dp)
        ok = .false.
        if (size(values) /= size(x_cells)) return
        call reconstruct(state, cosine, sine, basis, x_cells, velocity, work, ok)
        if (ok) values = work%candidate
    end subroutine

    subroutine sample_transport_cached(state, basis, x_cells, velocity, work, &
        values, ok)
        type(transport_state_t), intent(in) :: state
        type(transport_basis_t), intent(in) :: basis
        type(interval_t), intent(in) :: x_cells(:), velocity
        type(transport_sampling_workspace_t), intent(inout) :: work
        type(interval_t), intent(out) :: values(:)
        logical, intent(out) :: ok

        ok = .false.
        values = interval(0.0_dp)
        if (.not. allocated(state%cosine)) return
        if (.not. allocated(state%sine)) return
        call sample_transport_coefficients_cached(state, state%cosine, state%sine, &
            basis, x_cells, velocity, work, values, ok)
    end subroutine

    subroutine sample_transport_input_cached(state, basis, x_cells, velocity, work, &
        values, ok)
        type(transport_state_t), intent(in) :: state
        type(transport_basis_t), intent(in) :: basis
        type(interval_t), intent(in) :: x_cells(:), velocity
        type(transport_sampling_workspace_t), intent(inout) :: work
        type(interval_t), intent(out) :: values(:)
        logical, intent(out) :: ok

        ok = .false.
        values = interval(0.0_dp)
        if (.not. allocated(state%source_cos)) return
        if (.not. allocated(state%source_sin)) return
        call sample_transport_coefficients_cached(state, state%source_cos, &
            state%source_sin, basis, x_cells, velocity, work, values, ok)
    end subroutine

    logical function prepared_layout(state, certificate) result(ok)
        type(transport_state_t), intent(in) :: state
        type(transport_certificate_t), intent(in) :: certificate

        ok = .false.
        if (.not. certificate%valid) return
        if (.not. allocated(state%v)) return
        if (.not. allocated(state%cosine)) return
        if (.not. allocated(state%sine)) return
        if (size(state%v) /= certificate%nodes) return
        if (size(state%cosine, 1) /= certificate%modes + 1) return
        if (size(state%cosine, 2) /= certificate%nodes) return
        if (any(shape(state%cosine) /= shape(state%sine))) return
        if (.not. ieee_is_finite(state%nu)) return
        if (state%nu < certificate%nu .or. state%nu > certificate%nu) return
        ok = .true.
    end function

    subroutine widen(work, radius, lower, upper, ok)
        type(transport_sampling_workspace_t), intent(in) :: work
        real(dp), intent(in) :: radius
        real(dp), intent(out) :: lower(:), upper(:)
        logical, intent(out) :: ok
        integer :: i

        ok = .false.
        if (size(lower) /= size(work%candidate)) return
        if (size(upper) /= size(work%candidate)) return
        do i = 1, size(lower)
            call transport_observable_band(work%candidate(i), interval(1.0_dp), &
                radius, lower(i), upper(i), ok)
            if (.not. ok) return
        end do
        ok = .true.
    end subroutine

    subroutine certify_transport_cached(state, certificate, basis, x_cells, &
        velocity, work, lower, upper, ok)
        type(transport_state_t), intent(in) :: state
        type(transport_certificate_t), intent(in) :: certificate
        type(transport_basis_t), intent(in) :: basis
        type(interval_t), intent(in) :: x_cells(:), velocity
        type(transport_sampling_workspace_t), intent(inout) :: work
        real(dp), intent(out) :: lower(:), upper(:)
        logical, intent(out) :: ok

        ok = .false.
        lower = -huge(0.0_dp); upper = huge(0.0_dp)
        if (.not. prepared_layout(state, certificate)) return
        call reconstruct(state, state%cosine, state%sine, basis, x_cells, &
            velocity, work, ok)
        if (.not. ok) return
        call widen(work, certificate%state_radius, lower, upper, ok)
    end subroutine

    subroutine certify_transport_tangent_cached(state, certificate, dcos, dsin, &
        basis, x_cells, velocity, work, lower, upper, ok)
        type(transport_state_t), intent(in) :: state
        type(transport_certificate_t), intent(in) :: certificate
        real(dp), intent(in) :: dcos(0:, :), dsin(0:, :)
        type(transport_basis_t), intent(in) :: basis
        type(interval_t), intent(in) :: x_cells(:), velocity
        type(transport_sampling_workspace_t), intent(inout) :: work
        real(dp), intent(out) :: lower(:), upper(:)
        logical, intent(out) :: ok

        ok = .false.
        lower = -huge(0.0_dp); upper = huge(0.0_dp)
        if (.not. certificate%tangent_valid) return
        if (.not. prepared_layout(state, certificate)) return
        call reconstruct(state, dcos, dsin, basis, x_cells, velocity, work, ok)
        if (.not. ok) return
        call widen(work, certificate%tangent_radius, lower, upper, ok)
    end subroutine

    integer(int64) function transport_sampling_bytes(basis, work) result(bytes)
        type(transport_basis_t), intent(in) :: basis
        type(transport_sampling_workspace_t), intent(in) :: work
        type(interval_t) :: cell
        integer(int64) :: count

        count = 0_int64
        if (allocated(basis%x)) count = count + size(basis%x, kind=int64)
        if (allocated(basis%cosine)) count = count + size(basis%cosine, kind=int64)
        if (allocated(basis%sine)) count = count + size(basis%sine, kind=int64)
        if (allocated(work%candidate)) count = count + size(work%candidate, kind=int64)
        if (allocated(work%cell_range)) &
            count = count + size(work%cell_range, kind=int64)
        bytes = count*int(storage_size(cell)/8, int64)
    end function
end module kin6d_transport_sampling
