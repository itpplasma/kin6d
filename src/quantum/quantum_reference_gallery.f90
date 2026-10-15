module kin6d_quantum_reference_gallery
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortnum_interval, only: interval_t, interval, is_empty, mag, iadd, isub, &
        linear_interpolation_remainder
    use kin6d_quantum_reference_gallery_leaves, only: qgallery_accumulate, &
        qgallery_accumulate_interval, &
        qgallery_norm_add_interval, &
        qgallery_error_interval, qgallery_curvature_interval, &
        qgallery_band_interval, &
        qgallery_normalize_interval, &
        qgallery_hamiltonian_defect_interval, qgallery_voxel_interval, &
        qgallery_outer_interval, &
        qgallery_sector_error_interval, &
        qgallery_phase, qgallery_phase_interval, &
        qgallery_contract, &
        qgallery_contract_interval
    use kin6d_quantum_correlations_leaves, only: quantum_cumulant_pair_interval
    implicit none
    private
    public :: quantum_observable_packet_t, prepare_quantum_observables
    public :: sample_quantum_observables, quantum_occupation_bands
    public :: prepare_quantum_pair_sector
    public :: quantum_initial_norm, quantum_matrix_defect, quantum_physical_radius
    integer, parameter, public :: quantum_observable_count = 6
    type :: quantum_observable_packet_t
        private
        integer :: m = 0, modes = -1
        real(dp) :: voxel_width(3) = 0.0_dp
        type(interval_t), allocatable :: cr(:, :), ci(:, :), gr(:, :), gi(:, :)
        type(interval_t), allocatable :: occupations(:)
        type(interval_t) :: error(6), padding(6)
        logical :: ready = .false., rdm_scratch = .false.
    end type
contains
    subroutine quantum_initial_norm(re, im, norm2, ok)
        real(dp), intent(in) :: re(:), im(:)
        type(interval_t), intent(out) :: norm2
        logical, intent(out) :: ok
        type(interval_t) :: next
        integer :: i
        ok = .false.; norm2 = interval(0.0_dp)
        if (size(re) /= size(im)) return
        if (size(re) < 1) return
        if (.not. all(ieee_is_finite(re))) return
        if (.not. all(ieee_is_finite(im))) return
        do i = 1, size(re)
            call qgallery_norm_add_interval(norm2, interval(re(i)), interval(im(i)), &
                next)
            norm2 = next
        end do
        if (.not. valid_interval(norm2)) return
        ok = norm2%lo > 0.0_dp
    end subroutine

    subroutine quantum_matrix_defect(hr, hi, hri, hii, defect, ok)
        real(dp), intent(in) :: hr(:, :), hi(:, :)
        type(interval_t), intent(in) :: hri(:, :), hii(:, :)
        real(dp), intent(out) :: defect
        logical, intent(out) :: ok
        type(interval_t) :: acc, next, dr, di
        integer :: i, j
        ok = .false.; defect = huge(1.0_dp)
        if (any(shape(hr) /= shape(hi))) return
        if (any(shape(hr) /= shape(hri))) return
        if (any(shape(hr) /= shape(hii))) return
        acc = interval(0.0_dp)
        do j = 1, size(hr, 2)
            do i = 1, size(hr, 1)
                if (.not. valid_interval(hri(i, j))) return
                if (.not. valid_interval(hii(i, j))) return
                dr = isub(hri(i, j), interval(hr(i, j)))
                di = isub(hii(i, j), interval(hi(i, j)))
                call qgallery_norm_add_interval(acc, interval(mag(dr)), &
                    interval(mag(di)), next)
                acc = next
            end do
        end do
        ! The generic interval square root encloses the Frobenius norm.
        call norm_root(acc, next)
        if (.not. valid_interval(next)) return
        defect = next%hi; ok = .true.
    end subroutine

    subroutine norm_root(value, result)
        use fortnum_interval, only: sqrt
        type(interval_t), intent(in) :: value
        type(interval_t), intent(out) :: result
        result = sqrt(value)
    end subroutine

    subroutine quantum_physical_radius(norm0, radius, defect, time, hbar, &
        physical_radius, ok)
        type(interval_t), intent(in) :: norm0
        real(dp), intent(in) :: radius, defect, time, hbar
        real(dp), intent(out) :: physical_radius
        logical, intent(out) :: ok
        type(interval_t) :: result
        ok = .false.; physical_radius = huge(1.0_dp)
        if (.not. valid_interval(norm0)) return
        if (norm0%lo <= 0.0_dp) return
        if (.not. all(ieee_is_finite([radius, defect, time, hbar]))) return
        if (min(radius, defect, time) < 0.0_dp) return
        if (hbar <= 0.0_dp) return
        call qgallery_hamiltonian_defect_interval(norm0, interval(radius), &
            interval(defect), interval(time), interval(hbar), &
            result)
        if (.not. valid_interval(result)) return
        physical_radius = result%hi; ok = .true.
    end subroutine

    subroutine prepare_quantum_observables(packet, momentum, spin, pi, pj, &
        gamma_re, gamma_im, pair_re, pair_im, candidate_norm2, norm0, &
        radius, &
        max_gap, ok, voxel_width)
        type(quantum_observable_packet_t), intent(inout) :: packet
        integer, intent(in) :: momentum(:, :), spin(:), pi(:), pj(:)
        type(interval_t), intent(in) :: gamma_re(:, :), gamma_im(:, :)
        type(interval_t), intent(in) :: pair_re(:, :), pair_im(:, :)
        type(interval_t), intent(in) :: candidate_norm2, norm0
        real(dp), intent(in) :: radius, max_gap
        real(dp), intent(in), optional :: voxel_width(3)
        logical, intent(out) :: ok
        type(interval_t) :: pr, pimag, wr, wi, kr, ki, occupation_error
        type(interval_t) :: pair_error, connected_error, curvature, next
        integer :: m, p, modes, i, j, k, l, a, b, channel, s, count1, count2
        integer :: pair_index, allocation_status
        logical, allocatable :: seen_pairs(:)
        ok = .false.; packet%ready = .false.
        packet%voxel_width = 0.0_dp
        if (present(voxel_width)) then
            if (.not. all(ieee_is_finite(voxel_width))) return
            if (minval(voxel_width) < 0.0_dp .or. maxval(voxel_width) > 1.0_dp) return
            packet%voxel_width = voxel_width
        end if
        m = size(spin); p = size(pi)
        if (m < 1 .or. m > 10000) return
        if (size(momentum, 1) /= 3 .or. size(momentum, 2) /= m) return
        if (any(abs(spin) /= 1)) return
        if (any(momentum > 4096) .or. any(momentum < -4096)) return
        if (size(pj) /= p) return
        if (p /= m*(m - 1)/2) return
        if (any(shape(gamma_re) /= [m, m])) return
        if (any(shape(gamma_im) /= [m, m])) return
        if (any(shape(pair_re) /= [p, p])) return
        if (any(shape(pair_im) /= [p, p])) return
        if (any(pi < 1) .or. any(pj > m)) return
        if (any(pj <= pi)) return
        do j = 2, m
            do i = 1, j - 1
                if (spin(i) == spin(j)) then
                    if (all(momentum(:, i) == momentum(:, j))) return
                end if
            end do
        end do
        allocate (seen_pairs(p), stat=allocation_status)
        if (allocation_status /= 0) return
        seen_pairs = .false.
        do a = 1, p
            i = pi(a); j = pj(a)
            pair_index = (i - 1)*(2*m - i)/2 + j - i
            if (seen_pairs(pair_index)) return
            seen_pairs(pair_index) = .true.
        end do
        if (.not. valid_interval(candidate_norm2)) return
        if (.not. valid_interval(norm0)) return
        if (candidate_norm2%lo < 0.0_dp .or. norm0%lo <= 0.0_dp) return
        if (.not. all(ieee_is_finite([radius, max_gap]))) return
        if (min(radius, max_gap) < 0.0_dp) return
        modes = maxval(momentum(1, :)) - minval(momentum(1, :))
        if (packet%m /= m .or. packet%modes /= modes .or. .not. packet%rdm_scratch) then
            if (allocated(packet%cr)) deallocate (packet%cr, packet%ci, &
                packet%gr, packet%gi, &
                packet%occupations)
            allocate (packet%cr(-modes:modes, 6), packet%ci(-modes:modes, 6), &
                packet%gr(m, m), packet%gi(m, m), packet%occupations(m))
            packet%m = m; packet%modes = modes; packet%rdm_scratch = .true.
        end if
        packet%cr = interval(0.0_dp); packet%ci = interval(0.0_dp)
        do j = 1, m
            do i = 1, m
                if (.not. valid_interval(gamma_re(i, j))) return
                if (.not. valid_interval(gamma_im(i, j))) return
                call qgallery_normalize_interval(gamma_re(i, j), gamma_im(i, j), &
                    norm0, packet%gr(i, j), packet%gi(i, j))
            end do
        end do
        do i = 1, m
            packet%occupations(i) = packet%gr(i, i)
        end do
        do b = 1, p
            k = pi(b); l = pj(b)
            do a = 1, p
                i = pi(a); j = pj(a)
                if (.not. valid_interval(pair_re(a, b))) return
                if (.not. valid_interval(pair_im(a, b))) return
                call qgallery_normalize_interval(pair_re(a, b), pair_im(a, b), &
                    norm0, pr, pimag)
                call quantum_cumulant_pair_interval(pr, pimag, &
                    packet%gr(i, k), packet%gi(i, k), packet%gr(j, l), &
                    packet%gi(j, l), packet%gr(i, l), packet%gi(i, l), &
                    packet%gr(j, k), packet%gi(j, k), &
                    kr, ki)
                call add_pair(packet, momentum, spin, i, j, k, l, 1, -1, 1, pr, pimag)
                call add_pair(packet, momentum, spin, i, j, k, l, 1, -1, 2, kr, ki)
                call add_pair(packet, momentum, spin, i, j, k, l, 1, 1, 3, pr, pimag)
                call add_pair(packet, momentum, spin, i, j, k, l, 1, 1, 4, kr, ki)
            end do
        end do
        do j = 1, m
            do i = 1, m
                if (spin(i) /= spin(j)) cycle
                channel = 5
                if (spin(i) == -1) channel = 6
                call add_voxel_coefficient(packet, channel, momentum(:, &
                    i) - momentum(:, j), &
                    1, packet%gr(i, j), packet%gi(i, j))
            end do
        end do
        count1 = count(spin == 1); count2 = count(spin == -1)
        call qgallery_error_interval(candidate_norm2, interval(radius), &
            interval(real(count1, dp)), interval(real(count2, &
            dp)), &
            interval(0.0_dp), norm0, occupation_error, pair_error, &
            connected_error)
        packet%error(1) = pair_error; packet%error(2) = connected_error
        call qgallery_error_interval(candidate_norm2, interval(radius), &
            interval(real(count1, dp)), interval(real(count1, &
            dp)), &
            interval(1.0_dp), norm0, occupation_error, pair_error, &
            connected_error)
        packet%error(3) = pair_error; packet%error(4) = connected_error
        call qgallery_accumulate_interval(interval(0.0_dp), interval(0.0_dp), &
            occupation_error, interval(0.0_dp), interval(real(count1, &
            dp)), &
            packet%error(5), next)
        call qgallery_accumulate_interval(interval(0.0_dp), interval(0.0_dp), &
            occupation_error, interval(0.0_dp), interval(real(count2, &
            dp)), &
            packet%error(6), next)
        do i = 1, m
            call qgallery_band_interval(packet%occupations(i), occupation_error, wr, wi)
            packet%occupations(i) = interval(wr%lo, wi%hi)
            if (.not. valid_interval(packet%occupations(i))) return
        end do
        do channel = 1, 6
            curvature = interval(0.0_dp)
            do s = -modes, modes
                if (.not. valid_interval(packet%cr(s, channel))) return
                if (.not. valid_interval(packet%ci(s, channel))) return
                call qgallery_curvature_interval(curvature, packet%cr(s, channel), &
                    packet%ci(s, channel), interval(real(s, &
                    dp)), next)
                curvature = next
            end do
            ! The exact sum of magnitudes is nonnegative; discard a negative
            ! lower endpoint introduced by outward rounding of a zero sum.
            curvature = interval(0.0_dp, curvature%hi)
            packet%padding(channel) = linear_interpolation_remainder( &
                interval(0.0_dp), interval(max_gap), curvature)
            if (.not. valid_interval(packet%padding(channel))) return
            if (.not. valid_interval(packet%error(channel))) return
        end do
        packet%ready = .true.; ok = .true.
    end subroutine

    subroutine prepare_quantum_pair_sector(packet, momentum, re, im, norm0, &
        radius, max_gap, voxel_width, ok)
        type(quantum_observable_packet_t), intent(inout) :: packet
        integer, intent(in) :: momentum(:, :)
        real(dp), intent(in) :: re(:), im(:), radius, max_gap, voxel_width(3)
        type(interval_t), intent(in) :: norm0
        logical, intent(out) :: ok
        type(interval_t) :: norm2, probability, pr, pi, nr, ni, occupation_error
        type(interval_t) :: pair_error, exchange_error, next, curvature, unused
        integer :: count_modes, modes, a, b, channel, n, blocks, block_size
        ok = .false.; packet%ready = .false.
        count_modes = size(re)
        if (count_modes < 1 .or. count_modes > 10000) return
        if (size(im) /= count_modes) return
        if (any(shape(momentum) /= [3, count_modes])) return
        if (any(momentum > 4096) .or. any(momentum < -4096)) return
        if (.not. valid_interval(norm0)) return
        if (norm0%lo <= 0.0_dp) return
        if (.not. all(ieee_is_finite([radius, max_gap, voxel_width]))) return
        if (min(radius, max_gap, minval(voxel_width)) < 0.0_dp) return
        if (maxval(voxel_width) > 1.0_dp) return
        call quantum_initial_norm(re, im, norm2, ok)
        if (.not. ok) return
        ok = .false.
        do a = 1, count_modes
            do b = 1, a - 1
                if (all(momentum(:, a) == momentum(:, b))) return
            end do
        end do
        modes = maxval(momentum(1, :)) - minval(momentum(1, :))
    if (packet%m /= count_modes .or. packet%modes /= modes .or. packet%rdm_scratch) then
            if (allocated(packet%cr)) deallocate (packet%cr, packet%ci, &
                packet%gr, packet%gi, &
                packet%occupations)
            allocate (packet%cr(-modes:modes, 6), packet%ci(-modes:modes, 6), &
                packet%gr(0, 0), &
                packet%gi(0, 0), &
                packet%occupations(count_modes))
            packet%m = count_modes; packet%modes = modes; packet%rdm_scratch = .false.
        end if
        packet%voxel_width = voxel_width
        packet%cr = interval(0.0_dp); packet%ci = interval(0.0_dp)
        do a = 1, count_modes
            call qgallery_norm_add_interval(interval(0.0_dp), interval(re(a)), &
                interval(im(a)), probability)
            call qgallery_normalize_interval(probability, interval(0.0_dp), &
                norm0, nr, unused)
            packet%occupations(a) = nr
        end do
        do b = 1, count_modes
            do a = 1, count_modes
                if (voxel_width(2) >= 1.0_dp) then
                    if (momentum(2, a) /= momentum(2, b)) cycle
                end if
                if (voxel_width(3) >= 1.0_dp) then
                    if (momentum(3, a) /= momentum(3, b)) cycle
                end if
                call qgallery_outer_interval(interval(re(a)), interval(im(a)), &
                    interval(re(b)), interval(im(b)), &
                    pr, pi)
                call qgallery_normalize_interval(pr, pi, norm0, nr, ni)
                call add_voxel_coefficient(packet, 1, &
                    momentum(:, a) - momentum(:, b), 1, nr, ni)
                call qgallery_accumulate_interval(interval(0.0_dp), interval(0.0_dp), &
                    packet%occupations(a), interval(0.0_dp), &
                    packet%occupations(b), &
                    nr, ni)
                call add_voxel_coefficient(packet, 4, &
                    momentum(:, a) - momentum(:, b), 1, nr, ni)
            end do
        end do
        packet%cr(:, 2) = packet%cr(:, 1); packet%ci(:, 2) = packet%ci(:, 1)
        packet%cr(0, 2) = isub(packet%cr(0, 2), interval(1.0_dp))
        packet%cr(0, 4) = isub(packet%cr(0, 4), interval(1.0_dp))
        packet%cr(0, 5:6) = interval(1.0_dp)
        ! Averaging complete transverse periods removes nonzero ky/kz transfers.
        ! On each surviving block, every pair-observable matrix entry has |a|<=1.
        blocks = count_modes
        if (voxel_width(2) >= 1.0_dp .and. voxel_width(3) >= 1.0_dp) then
            blocks = 0
            do a = 1, count_modes
                block_size = 0
                do b = 1, count_modes
                    if (all(momentum(2:3, a) == momentum(2:3, b))) &
                        block_size = block_size + 1
                end do
                blocks = max(blocks, block_size)
            end do
        end if
        call qgallery_sector_error_interval(norm2, interval(radius), &
            interval(real(blocks, dp)), norm0, occupation_error, pair_error, &
                exchange_error)
        packet%error = interval(0.0_dp)
        packet%error(1:2) = pair_error; packet%error(4) = exchange_error
        do a = 1, count_modes
            pr = isub(packet%occupations(a), occupation_error)
            pi = iadd(packet%occupations(a), occupation_error)
            packet%occupations(a) = interval(pr%lo, pi%hi)
        end do
        do channel = 1, 6
            curvature = interval(0.0_dp)
            do n = -modes, modes
                call qgallery_curvature_interval(curvature, packet%cr(n, channel), &
                    packet%ci(n, channel), interval(real(n, &
                    dp)), next)
                curvature = next
            end do
            if (.not. valid_interval(curvature)) return
            curvature = interval(0.0_dp, curvature%hi)
            packet%padding(channel) = linear_interpolation_remainder( &
                interval(0.0_dp), interval(max_gap), curvature)
            if (.not. valid_interval(packet%padding(channel))) return
            if (.not. valid_interval(packet%error(channel))) return
        end do
        packet%ready = .true.; ok = .true.
    end subroutine

    subroutine add_pair(packet, momentum, spin, i, j, k, l, sigma, tau, &
        channel, re, im)
        type(quantum_observable_packet_t), intent(inout) :: packet
        integer, intent(in) :: momentum(:, :), spin(:), i, j, k, l, sigma, tau
        integer, intent(in) :: channel
        type(interval_t), intent(in) :: re, im
        integer :: a, b, left(2), right(2), lp(2), rp(2), difference(3)
        left = [i, j]; right = [k, l]; lp = 0; rp = 0
        if (spin(i) == sigma .and. spin(j) == tau) lp(1) = 1
        if (spin(j) == sigma .and. spin(i) == tau) lp(2) = -1
        if (spin(k) == sigma .and. spin(l) == tau) rp(1) = 1
        if (spin(l) == sigma .and. spin(k) == tau) rp(2) = -1
        do b = 1, 2
            if (rp(b) == 0) cycle
            do a = 1, 2
                if (lp(a) == 0) cycle
                difference = momentum(:, left(a)) - momentum(:, right(b))
                call add_voxel_coefficient(packet, channel, difference, &
                    lp(a)*rp(b), re, im)
            end do
        end do
    end subroutine

    subroutine add_voxel_coefficient(packet, channel, difference, sign, re, im)
        type(quantum_observable_packet_t), intent(inout) :: packet
        integer, intent(in) :: channel, difference(3), sign
        type(interval_t), intent(in) :: re, im
        type(interval_t) :: weight, factor, next, ignored, rr, ri
        integer :: axis
        weight = interval(real(sign, dp))
        do axis = 1, 3
            if (difference(axis) == 0) cycle
            if (packet%voxel_width(axis) <= 0.0_dp) cycle
            if (packet%voxel_width(axis) >= 1.0_dp) then
                weight = interval(0.0_dp)
                exit
            end if
            call qgallery_voxel_interval(interval(real(difference(axis), dp)), &
                interval(packet%voxel_width(axis)), factor)
            call qgallery_accumulate_interval(interval(0.0_dp), interval(0.0_dp), &
                weight, interval(0.0_dp), factor, next, ignored)
            weight = next
        end do
        call qgallery_accumulate_interval(interval(0.0_dp), interval(0.0_dp), &
            re, im, weight, rr, ri)
        call add_coefficient(packet, channel, difference(1), 1, rr, ri)
    end subroutine

    subroutine add_coefficient(packet, channel, mode, weight, re, im)
        type(quantum_observable_packet_t), intent(inout) :: packet
        integer, intent(in) :: channel, mode, weight
        type(interval_t), intent(in) :: re, im
        type(interval_t) :: rr, ri
        call qgallery_accumulate_interval(packet%cr(mode, channel), &
            packet%ci(mode, channel), re, im, interval(real(weight, dp)), &
            rr, ri)
        packet%cr(mode, channel) = rr; packet%ci(mode, channel) = ri
    end subroutine

    subroutine sample_quantum_observables(packet, x, value, lower, upper, ok)
        use fortnum_interval, only: mid
        type(quantum_observable_packet_t), intent(in) :: packet
        real(dp), intent(in) :: x
        real(dp), intent(out) :: value(:), lower(:), upper(:)
        logical, intent(out) :: ok
        type(interval_t) :: acc(6), term, budget, lo, hi, cosine, sine
        real(dp) :: candidate, c, s, next_value, unused_imag
        integer :: channel, mode
        ok = .false.; value = 0.0_dp; lower = 0.0_dp; upper = 0.0_dp
        if (.not. packet%ready) return
        if (.not. ieee_is_finite(x)) return
        if (any([size(value), size(lower), size(upper)] /= 6)) return
        acc = interval(0.0_dp)
        do mode = -packet%modes, packet%modes
            call qgallery_phase(real(mode, dp), x, c, s)
            call qgallery_phase_interval(interval(real(mode, dp)), interval(x), &
                cosine, sine)
            do channel = 1, 6
                call qgallery_contract(mid(packet%cr(mode, channel)), &
                    mid(packet%ci(mode, channel)), c, s, candidate)
                call qgallery_accumulate(value(channel), 0.0_dp, candidate, &
                    0.0_dp, 1.0_dp, next_value, unused_imag)
                value(channel) = next_value
                call qgallery_contract_interval(packet%cr(mode, channel), &
                    packet%ci(mode, channel), cosine, &
                    sine, term)
                acc(channel) = iadd(acc(channel), term)
            end do
        end do
        do channel = 1, 6
            budget = iadd(packet%error(channel), packet%padding(channel))
            call qgallery_band_interval(acc(channel), budget, lo, hi)
            if (.not. valid_interval(lo)) return
            if (.not. valid_interval(hi)) return
            lower(channel) = lo%lo; upper(channel) = hi%hi
        end do
        ok = all(ieee_is_finite(value))
    end subroutine

    subroutine quantum_occupation_bands(packet, value, lower, upper, ok)
        use fortnum_interval, only: mid
        type(quantum_observable_packet_t), intent(in) :: packet
        real(dp), intent(out) :: value(:), lower(:), upper(:)
        logical, intent(out) :: ok
        integer :: i
        ok = .false.; value = 0.0_dp; lower = 0.0_dp; upper = 0.0_dp
        if (.not. packet%ready) return
        if (any([size(value), size(lower), size(upper)] /= packet%m)) return
        do i = 1, packet%m
            value(i) = mid(packet%occupations(i))
            lower(i) = packet%occupations(i)%lo; upper(i) = packet%occupations(i)%hi
        end do
        ok = .true.
    end subroutine

    pure logical function valid_interval(value)
        type(interval_t), intent(in) :: value
        valid_interval = .false.
        if (is_empty(value)) return
        if (.not. ieee_is_finite(value%lo)) return
        if (.not. ieee_is_finite(value%hi)) return
        valid_interval = value%lo <= value%hi
    end function
end module kin6d_quantum_reference_gallery
