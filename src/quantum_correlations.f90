module kin6d_quantum_correlations
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortnum_interval, only: interval_t, interval
    use kin6d_quantum_correlations_leaves, only: quantum_wick_pair, &
        quantum_wick_pair_jvp, quantum_wick_pair_vjp, quantum_wick_pair_interval, &
        quantum_cumulant_pair, quantum_cumulant_pair_jvp, quantum_cumulant_pair_vjp, &
        quantum_cumulant_pair_interval, quantum_add_interval, quantum_negate, &
        quantum_magnitude_bound_interval, quantum_gershgorin_interval
    implicit none
    private
    public :: quantum_wick_pair, quantum_wick_pair_jvp, quantum_wick_pair_vjp
    public :: quantum_wick_pair_interval, quantum_cumulant_pair
    public :: quantum_cumulant_pair_jvp, quantum_cumulant_pair_vjp
    public :: quantum_cumulant_pair_interval, quantum_pauli
    public :: pauli_certified, pauli_violated, pauli_undetermined
    integer, parameter :: pauli_certified = 1, pauli_violated = -1
    integer, parameter :: pauli_undetermined = 0
contains
    ! Sufficient Pauli test of the EXACT supplied binary64 matrix, not of an
    ! uncertain physical matrix. Failure of Gershgorin sufficiency is undetermined.
    ! Non-Hermitian data or outside-[0,1] diagonal data violate necessary conditions.
    subroutine quantum_pauli(gr, gi, status, ok)
        real(dp), intent(in) :: gr(:, :), gi(:, :)
        integer, intent(out) :: status
        logical, intent(out) :: ok
        integer :: n, i, j
        real(dp) :: negative
        type(interval_t) :: radius, magnitude, next, lower, upper
        logical :: sufficient, offdiagonal, hermitian
        ok = .false.; status = pauli_undetermined
        n = size(gr,1)
        if (n < 1) return
        if (.not. all(shape(gr) == [n,n])) return
        if (.not. all(shape(gi) == [n,n])) return
        if (.not. all(ieee_is_finite(gr))) return
        if (.not. all(ieee_is_finite(gi))) return
        hermitian = .true.
        do j = 1, n
            do i = 1, n
                call quantum_negate(gi(j,i),negative)
                if (gr(i,j) < gr(j,i) .or. gr(i,j) > gr(j,i)) hermitian = .false.
                if (gi(i,j) < negative .or. gi(i,j) > negative) hermitian = .false.
            end do
        end do
        ok = .true.; status = pauli_violated
        if (.not. hermitian) return
        do i = 1, n
            if (gr(i,i) < 0 .or. gr(i,i) > 1) return
        end do
        status = pauli_undetermined; sufficient = .true.
        do i = 1, n
            radius = interval(0.0_dp); offdiagonal = .false.
            do j = 1, n
                if (j == i) cycle
                if (abs(gr(i,j)) <= 0 .and. abs(gi(i,j)) <= 0) cycle
                offdiagonal = .true.
                call quantum_magnitude_bound_interval(interval(gr(i,j)), &
                    interval(gi(i,j)),magnitude)
                call quantum_add_interval(radius,magnitude,next)
                radius = next
                if (.not. finite_interval(radius)) then
                    ok = .false.
                    return
                end if
            end do
            ! Exact diagonal rows include eigenvalue endpoints 0 and 1.
            if (.not. offdiagonal) cycle
            call quantum_gershgorin_interval(interval(gr(i,i)),radius,lower,upper)
            if (.not. finite_interval(lower)) then
                ok = .false.
                return
            end if
            if (.not. finite_interval(upper)) then
                ok = .false.
                return
            end if
            if (lower%lo < 0 .or. upper%hi > 1) sufficient = .false.
        end do
        if (sufficient) status = pauli_certified
    end subroutine

    logical function finite_interval(value) result(ok)
        type(interval_t), intent(in) :: value
        ok = ieee_is_finite(value%lo) .and. ieee_is_finite(value%hi)
        if (.not. ok) return
        ok = value%lo <= value%hi
    end function
end module kin6d_quantum_correlations
