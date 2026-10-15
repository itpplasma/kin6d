module kin6d_quantum_correlations
    use, intrinsic :: iso_fortran_env, only: dp => real64, i64 => int64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use kin6d_sort_keys, only: sort_key_order
    use fortnum_interval, only: interval_t, interval
    use kin6d_quantum_correlations_leaves, only: quantum_wick_pair, &
        quantum_wick_pair_jvp, quantum_wick_pair_vjp, quantum_wick_pair_interval, &
        quantum_cumulant_pair, quantum_cumulant_pair_jvp, quantum_cumulant_pair_vjp, &
        quantum_cumulant_pair_interval, quantum_add_interval, quantum_negate, &
        quantum_magnitude_bound_interval, quantum_gershgorin_interval, &
        quantum_wick_triple, quantum_triple_product, quantum_square_add
    implicit none
    private
    public :: quantum_wick_pair, quantum_wick_pair_jvp, quantum_wick_pair_vjp
    public :: quantum_wick_pair_interval, quantum_cumulant_pair
    public :: quantum_cumulant_pair_jvp, quantum_cumulant_pair_vjp
    public :: quantum_cumulant_pair_interval, quantum_pauli
    public :: pauli_certified, pauli_violated, pauli_undetermined
    public :: quantum_cumulant_norms
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
    !> Connected two- and three-body cumulants of real RDMs (real states or
    !> ensembles). gamma(m,m) with gamma_ij=<a_j^dagger a_i>; gamma2(p,p) in
    !> ordered-pair storage; Gamma3 as sparse reduced coordinates with triple
    !> keys ((t1-1)m+t2-1)m+t3-1 (as returned by fock_rdm3_entries).
    !>   lambda2 = Gamma2 - W,  W = gamma ^ gamma (pair Wick term);
    !>   lambda3 = Gamma3 - det3[gamma(t_a,u_b)]
    !>             - sum_ab (-1)^(a+b) gamma(t_a,u_b) lambda2(t\a,u\b).
    !> Norms are Frobenius norms over ordered pairs i<j / triples i<j<k.
    !> lambda3 is assembled on the union of the supports of its terms.
    subroutine quantum_cumulant_norms(m, gamma, gamma2, tkey, ukey, g3, &
                                      norm_gamma2, norm_lambda2, norm_gamma3, &
                                      norm_lambda3, ok)
        integer, intent(in) :: m
        real(dp), intent(in) :: gamma(:, :), gamma2(:, :), g3(:)
        integer(i64), intent(in) :: tkey(:), ukey(:)
        real(dp), intent(out) :: norm_gamma2, norm_lambda2, norm_gamma3, norm_lambda3
        logical, intent(out) :: ok
        real(dp), allocatable :: lambda2(:, :), val(:)
        integer(i64), allocatable :: key(:)
        integer, allocatable :: perm(:), pi(:), pj(:)
        integer :: p, q, i, j, k, x, y, a, b, nnz, cursor, np
        integer :: t(3), u(3), rows(3)
        integer(i64) :: m3
        real(dp) :: rr, ri, next, total, term, gm(3, 3)

        ok = .false.
        norm_gamma2 = 0; norm_lambda2 = 0; norm_gamma3 = 0; norm_lambda3 = 0
        np = m*(m - 1)/2
        if (m < 1 .or. m > 1000) return
        if (any(shape(gamma) /= [m, m]) .or. any(shape(gamma2) /= [np, np])) return
        if (size(tkey) /= size(ukey) .or. size(tkey) /= size(g3)) return
        if (.not. all(ieee_is_finite(gamma)) .or. .not. all(ieee_is_finite(gamma2))) return
        if (.not. all(ieee_is_finite(g3))) return
        m3 = int(m, i64)**3
        allocate (pi(np), pj(np), lambda2(np, np))
        p = 0
        do i = 1, m - 1
            do j = i + 1, m
                p = p + 1; pi(p) = i; pj(p) = j
            end do
        end do
        do q = 1, np
            do p = 1, np
                call quantum_cumulant_pair(gamma2(p, q), 0.0_dp, &
                    gamma(pi(p), pi(q)), 0.0_dp, gamma(pj(p), pj(q)), 0.0_dp, &
                    gamma(pi(p), pj(q)), 0.0_dp, gamma(pj(p), pi(q)), 0.0_dp, rr, ri)
                lambda2(p, q) = rr
                call quantum_square_add(norm_gamma2, gamma2(p, q), next)
                norm_gamma2 = next
                call quantum_square_add(norm_lambda2, rr, next)
                norm_lambda2 = next
            end do
        end do
        do k = 1, size(g3)
            call quantum_square_add(norm_gamma3, g3(k), next)
            norm_gamma3 = next
        end do
        ! Contribution count: Gamma3 + Wick candidates + gamma x lambda2.
        nnz = size(g3)
        do i = 1, m
            do j = 1, m
                if (abs(gamma(i, j)) > 0) nnz = nnz + count(abs(lambda2) > 0)
            end do
        end do
        nnz = nnz + wick_candidates(.false.)
        allocate (key(nnz), val(nnz))
        cursor = 0
        do k = 1, size(g3)
            cursor = cursor + 1; key(cursor) = tkey(k)*m3 + ukey(k); val(cursor) = g3(k)
        end do
        cursor = cursor + wick_candidates(.true.)
        do q = 1, np
            do p = 1, np
                if (.not. abs(lambda2(p, q)) > 0) cycle
                do x = 1, m
                    if (x == pi(p) .or. x == pj(p)) cycle
                    do y = 1, m
                        if (.not. abs(gamma(x, y)) > 0) cycle
                        if (y == pi(q) .or. y == pj(q)) cycle
                        call insert(pi(p), pj(p), x, t, a)
                        call insert(pi(q), pj(q), y, u, b)
                        call quantum_triple_product(real(-(-1)**(a + b), dp), &
                                                    gamma(x, y), lambda2(p, q), term)
                        cursor = cursor + 1
                        key(cursor) = triple_key(t)*m3 + triple_key(u)
                        val(cursor) = term
                    end do
                end do
            end do
        end do
        allocate (perm(cursor))
        call sort_key_order(key(1:cursor), perm)
        x = 1
        do while (x <= cursor)
            total = 0
            y = x
            do while (y <= cursor)
                if (key(perm(y)) /= key(perm(x))) exit
                total = total + val(perm(y))
                y = y + 1
            end do
            call quantum_square_add(norm_lambda3, total, next)
            norm_lambda3 = next
            x = y
        end do
        norm_gamma2 = sqrt(norm_gamma2); norm_lambda2 = sqrt(norm_lambda2)
        norm_gamma3 = sqrt(norm_gamma3); norm_lambda3 = sqrt(norm_lambda3)
        ok = ieee_is_finite(norm_lambda3) .and. ieee_is_finite(norm_lambda2)
    contains
        integer(i64) function triple_key(s)
            integer, intent(in) :: s(3)
            triple_key = (int(s(1) - 1, i64)*m + (s(2) - 1))*m + (s(3) - 1)
        end function
        subroutine insert(c1, c2, z, s, position)
            integer, intent(in) :: c1, c2, z
            integer, intent(out) :: s(3), position
            if (z < c1) then
                s = [z, c1, c2]; position = 1
            else if (z < c2) then
                s = [c1, z, c2]; position = 2
            else
                s = [c1, c2, z]; position = 3
            end if
        end subroutine
        !> Wick det3 contributions -det[gamma(t_a,u_b)] for every triple t and
        !> every sorted triple u of columns reachable from the nonzero rows.
        integer function wick_candidates(store) result(added)
            logical, intent(in) :: store
            integer :: c1, c2, c3, s1, s2, s3, e
            real(dp) :: det
            integer(i64), allocatable :: seen(:)
            integer :: nseen
            added = 0
            allocate (seen(m**3))
            do s1 = 1, m - 2
                do s2 = s1 + 1, m - 1
                    do s3 = s2 + 1, m
                        rows = [s1, s2, s3]
                        nseen = 0
                        do c1 = 1, m
                            if (.not. abs(gamma(s1, c1)) > 0) cycle
                            do c2 = 1, m
                                if (c2 == c1 .or. .not. abs(gamma(s2, c2)) > 0) cycle
                                do c3 = 1, m
                                    if (c3 == c1 .or. c3 == c2) cycle
                                    if (.not. abs(gamma(s3, c3)) > 0) cycle
                                    u = [c1, c2, c3]
                                    call sort3(u)
                                    if (any(seen(1:nseen) == triple_key(u))) cycle
                                    nseen = nseen + 1; seen(nseen) = triple_key(u)
                                    added = added + 1
                                    if (.not. store) cycle
                                    do e = 1, 3
                                        gm(e, :) = gamma(rows(e), u)
                                    end do
                                    call quantum_wick_triple(gm(1, 1), gm(1, 2), &
                                        gm(1, 3), gm(2, 1), gm(2, 2), gm(2, 3), &
                                        gm(3, 1), gm(3, 2), gm(3, 3), det)
                                    key(cursor + added) = triple_key(rows)*m3 + &
                                                          triple_key(u)
                                    val(cursor + added) = -det
                                end do
                            end do
                        end do
                    end do
                end do
            end do
        end function
        subroutine sort3(s)
            integer, intent(inout) :: s(3)
            integer :: tmp
            if (s(1) > s(2)) then
                tmp = s(1); s(1) = s(2); s(2) = tmp
            end if
            if (s(2) > s(3)) then
                tmp = s(2); s(2) = s(3); s(3) = tmp
            end if
            if (s(1) > s(2)) then
                tmp = s(1); s(1) = s(2); s(2) = tmp
            end if
        end subroutine
    end subroutine
end module kin6d_quantum_correlations
