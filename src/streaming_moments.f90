module kin6d_streaming_moments
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortnum_interval, only: interval_t, interval, hull
    use kin6d_streaming_comparison, only: streaming_state_t, init_streaming, &
        advance_streaming, density_coefficients_streaming, streaming_radius, &
        evaluate_streaming
    use kin6d_streaming_leaves, only: streaming_rotation, streaming_rotation_interval, &
        streaming_mode_accumulate, streaming_mode_accumulate_interval, &
        streaming_enlarge_interval
    use kin6d_streaming_moments_leaves, only: streaming_moments_current_cell, &
        streaming_moments_current_cell_interval, &
        streaming_moments_current_radius_interval, &
        streaming_moments_bulk, streaming_moments_bulk_interval
    implicit none
    private
    public :: streaming_moments_t, initialize_streaming_moments
    public :: advance_streaming_moments
    public :: sample_streaming_moments, sample_streaming_moments_box
    public :: sample_streaming_moments_distribution

    type :: streaming_moments_t
        private
        type(streaming_state_t) :: distribution
        real(dp), allocatable :: v(:), a0(:, :), b0(:, :), a(:, :), b(:, :)
        real(dp), allocatable :: nc(:), ns(:), jc(:), js(:)
        real(dp), allocatable :: ncl(:), nch(:), nsl(:), nsh(:)
        type(interval_t), allocatable :: ai(:, :), bi(:, :), ni(:), nsi(:)
        type(interval_t), allocatable :: ji(:), jsi(:)
        type(interval_t) :: nr, jr
        logical :: initialized = .false., ready = .false.
    end type
contains
    subroutine initialize_streaming_moments(packet, velocity, cosine0, sine0, ok)
        type(streaming_moments_t), intent(out) :: packet
        real(dp), intent(in) :: velocity(:), cosine0(0:, :), sine0(0:, :)
        logical, intent(out) :: ok
        integer :: modes, nodes

        call init_streaming(packet%distribution, velocity, cosine0, sine0, ok)
        if (.not. ok) return
        modes = size(cosine0, 1) - 1
        nodes = size(velocity)
        allocate (packet%v(nodes), packet%a0(0:modes, nodes), &
            packet%b0(0:modes, nodes), packet%a(0:modes, nodes), &
            packet%b(0:modes, nodes), packet%ai(0:modes, nodes), &
            packet%bi(0:modes, nodes), packet%nc(0:modes), packet%ns(0:modes), &
            packet%jc(0:modes), packet%js(0:modes), packet%ncl(0:modes), &
            packet%nch(0:modes), packet%nsl(0:modes), packet%nsh(0:modes), &
            packet%ni(0:modes), packet%nsi(0:modes), packet%ji(0:modes), &
            packet%jsi(0:modes))
        packet%v = velocity
        packet%a0 = cosine0
        packet%b0 = sine0
        packet%initialized = .true.
        call advance_streaming_moments(packet, 0.0_dp, ok)
    end subroutine

    subroutine advance_streaming_moments(packet, t, ok)
        type(streaming_moments_t), intent(inout) :: packet
        real(dp), intent(in) :: t
        logical, intent(out) :: ok
        type(interval_t) :: next, ai, bi
        real(dp) :: fnext, radius
        integer :: k, j

        ok = .false.
        packet%ready = .false.
        if (.not. packet%initialized) return
        call advance_streaming(packet%distribution, t, ok)
        if (.not. ok) return
        ok = .false.
        do j = 1, size(packet%v)
            do k = 0, size(packet%a0, 1) - 1
                ! Reuse the authoritative propagator; no second physics formula.
                call streaming_rotation(t, packet%v(j), real(k, dp), &
                    packet%a0(k, j), packet%b0(k, j), packet%a(k, j), packet%b(k, j))
                call streaming_rotation_interval(interval(t), interval(packet%v(j)), &
                    interval(real(k, dp)), interval(packet%a0(k, j)), &
                    interval(packet%b0(k, j)), ai, bi)
                packet%ai(k, j) = hull(ai, interval(packet%a(k, j)))
                packet%bi(k, j) = hull(bi, interval(packet%b(k, j)))
            end do
        end do
        packet%jc = 0.0_dp; packet%js = 0.0_dp
        packet%ji = interval(0.0_dp); packet%jsi = interval(0.0_dp)
        do k = 0, size(packet%jc) - 1
            do j = 1, size(packet%v) - 1
                call streaming_moments_current_cell(packet%v(j), packet%v(j + 1), &
                    packet%a(k, j), packet%a(k, j + 1), packet%jc(k), fnext)
                packet%jc(k) = fnext
                call streaming_moments_current_cell(packet%v(j), packet%v(j + 1), &
                    packet%b(k, j), packet%b(k, j + 1), packet%js(k), fnext)
                packet%js(k) = fnext
                call streaming_moments_current_cell_interval(interval(packet%v(j)), &
                    interval(packet%v(j + 1)), packet%ai(k, j), packet%ai(k, j + 1), &
                    packet%ji(k), next)
                packet%ji(k) = next
                call streaming_moments_current_cell_interval(interval(packet%v(j)), &
                    interval(packet%v(j + 1)), packet%bi(k, j), packet%bi(k, j + 1), &
                    packet%jsi(k), next)
                packet%jsi(k) = next
            end do
        end do
        if (.not. all(finite_interval(packet%ji))) return
        if (.not. all(finite_interval(packet%jsi))) return
        if (.not. all(ieee_is_finite(packet%jc))) return
        if (.not. all(ieee_is_finite(packet%js))) return
        call density_coefficients_streaming(packet%distribution, packet%nc, packet%ns, &
            packet%ncl, packet%nch, packet%nsl, packet%nsh, radius, ok)
        if (.not. ok) return
        ok = .false.
        packet%ni = interval(packet%ncl, packet%nch)
        packet%nsi = interval(packet%nsl, packet%nsh)
        packet%nr = interval(0.0_dp, radius)
        call streaming_radius(packet%distribution, radius, ok)
        if (.not. ok) return
        ok = .false.
        call streaming_moments_current_radius_interval(interval(0.0_dp, radius), &
            packet%jr)
        if (.not. finite_interval(packet%nr)) return
        if (.not. finite_interval(packet%jr)) return
        packet%ready = .true.
        ok = .true.
    end subroutine

    subroutine sample_streaming_moments_box(packet, x, density, current, bulk, ok)
        type(streaming_moments_t), intent(in) :: packet
        type(interval_t), intent(in) :: x
        type(interval_t), intent(out) :: density, current, bulk
        logical, intent(out) :: ok
        type(interval_t) :: n, j, next, lo, hi
        integer :: k

        ok = .false.
        density = interval(0.0_dp); current = density; bulk = density
        if (.not. packet%ready) return
        if (.not. finite_interval(x)) return
        n = interval(0.0_dp); j = n
        do k = 0, size(packet%nc) - 1
            call streaming_mode_accumulate_interval(x, interval(real(k, dp)), &
                packet%ni(k), packet%nsi(k), n, next)
            n = next
            call streaming_mode_accumulate_interval(x, interval(real(k, dp)), &
                packet%ji(k), packet%jsi(k), j, next)
            j = next
        end do
        call streaming_enlarge_interval(n, packet%nr, lo, hi)
        density = hull(lo, hi)
        call streaming_enlarge_interval(j, packet%jr, lo, hi)
        current = hull(lo, hi)
        if (.not. finite_interval(density) .or. .not. finite_interval(current)) return
        ! Bulk velocity is unsupported where a strictly positive mass is unproved.
        if (density%lo <= 0.0_dp) return
        call streaming_moments_bulk_interval(density, current, bulk)
        ok = finite_interval(bulk)
    end subroutine

    subroutine sample_streaming_moments(packet, x, n, j, u, nband, jband, uband, ok)
        type(streaming_moments_t), intent(in) :: packet
        real(dp), intent(in) :: x
        real(dp), intent(out) :: n, j, u
        type(interval_t), intent(out) :: nband, jband, uband
        logical, intent(out) :: ok
        real(dp) :: next
        integer :: k

        n = 0.0_dp; j = 0.0_dp; u = 0.0_dp
        call sample_streaming_moments_box(packet, interval(x), nband, jband, uband, ok)
        if (.not. ok) return
        ok = .false.
        do k = 0, size(packet%nc) - 1
            call streaming_mode_accumulate(x, real(k, dp), packet%nc(k), &
                packet%ns(k), n, next)
            n = next
            call streaming_mode_accumulate(x, real(k, dp), packet%jc(k), &
                packet%js(k), j, next)
            j = next
        end do
        if (.not. ieee_is_finite(n) .or. .not. ieee_is_finite(j)) return
        if (n <= 0.0_dp) return
        call streaming_moments_bulk(n, j, u)
        if (.not. ieee_is_finite(u)) return
        nband = hull(nband, interval(n))
        jband = hull(jband, interval(j))
        uband = hull(uband, interval(u))
        ok = .true.
    end subroutine

    subroutine sample_streaming_moments_distribution(packet, x, v, value, lo, hi, ok)
        type(streaming_moments_t), intent(in) :: packet
        real(dp), intent(in) :: x, v
        real(dp), intent(out) :: value, lo, hi
        logical, intent(out) :: ok

        ok = .false.
        if (.not. packet%ready) return
        call evaluate_streaming(packet%distribution, x, v, value, lo, hi, ok)
    end subroutine

    pure elemental logical function finite_interval(value) result(ok)
        type(interval_t), intent(in) :: value

        ok = ieee_is_finite(value%lo) .and. ieee_is_finite(value%hi) .and. &
            value%lo <= value%hi
    end function
end module kin6d_streaming_moments
