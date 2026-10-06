!> H1b quantum side for N = 3, 4 and the common difference measures.
!>
!> - jellium_full_spectrum: every eigenpair of every (S_z, K) sector of the
!>   Q4 cell by dense assembly of the matrix-free coupled action and LAPACK
!>   dsyev. Eigenvectors are kept up to a window above each sector minimum;
!>   the ordering is exact (dense), so the Q4 omitted-weight estimate
!>   (jellium_thermal_weights) is a bound up to the dense-eigensolver accuracy.
!> - jellium_diagonal_moments: the diagonal pair/triplet Fourier moments
!>   f2(Q) = <sum_{i/=j} exp(iQ.(x_i-x_j))>/(N(N-1)) and
!>   f3(Q1,Q2) = <sum_{i,j,k distinct} exp(i(Q1.x_i+Q2.x_j+Q3.x_k))>/(N(N-1)(N-2)),
!>   Q3 = -Q1-Q2, from the ensemble 2-RDM and the per-state 3-RDM entries. At
!>   nonzero Q these are the Fourier coefficients of the Ursell (connected)
!>   diagonal cumulants of a homogeneous state: the common classical/quantum
!>   measure (kin6d_classical_cell computes the same moments).
!> - ideal_cell_pair: spin-summed/-resolved shell-averaged pair functions of
!>   ideal fermions in the same box and canonical ensemble over all S_z,
!>   assembled from the Q3 fixed-N_sigma states (kin6d_fermion_ensemble).
!> - difference norms used by the regime map.
module kin6d_quantum_classical_map
    use, intrinsic :: iso_fortran_env, only: dp => real64, i64 => int64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortnum_quadrature, only: gauss_legendre_ab
    use kin6d_jellium_cell, only: jellium_cell_t, jellium_sector_t, &
        jellium_spectrum_t, jellium_sector, jellium_sector_momenta, &
        jellium_thermal_rdms
    use kin6d_quantum_fock, only: fock_hamiltonian_t, fock_matvec, fock_rdm3_entries
    use kin6d_fermion_ensemble, only: fermion_spectrum_t, canonical_state_t, &
        box_spectrum, box_cutoff, box_parameters, solve_canonical, same_spin_pair
    implicit none
    private
    public :: jellium_full_spectrum, jellium_diagonal_moments, ideal_cell_pair
    public :: shell_rms, pair_channels_classical
    real(dp), parameter :: pi = acos(-1.0_dp)

    interface
        subroutine dsyev(jobz, uplo, n, a, lda, w, work, lwork, info)
            import :: dp
            character, intent(in) :: jobz, uplo
            integer, intent(in) :: n, lda, lwork
            real(dp), intent(inout) :: a(lda, *)
            real(dp), intent(out) :: w(*), work(*)
            integer, intent(out) :: info
        end subroutine dsyev
    end interface
contains
    !> Full sector spectra (all eigenvalues; vectors with E - E_min(sector)
    !> <= window). kept(s)/top(s) record the stored set for the Q4 weights.
    subroutine jellium_full_spectrum(cell, model, window, spectrum, ok)
        type(jellium_cell_t), intent(in) :: cell
        type(fock_hamiltonian_t), intent(in) :: model
        real(dp), intent(in) :: window
        type(jellium_spectrum_t), intent(out) :: spectrum
        logical, intent(out) :: ok
        integer, allocatable :: ks(:, :)
        real(dp), allocatable :: h(:, :), e(:), x(:), zero(:), y(:), yi(:), work(:)
        integer :: nup, l, s, d, j, k, total, info, count
        integer(i64) :: used, need

        ok = .false.
        total = 0
        do nup = 0, cell%n
            if (nup > cell%nmodes .or. cell%n - nup > cell%nmodes) cycle
            call jellium_sector_momenta(cell, nup, ks)
            total = total + size(ks, 2)
        end do
        allocate (spectrum%sectors(total), spectrum%kept(total), spectrum%top(total))
        allocate (spectrum%sector_of(0), spectrum%offset(0), spectrum%energy(0), &
                  spectrum%residual(0), spectrum%flat(0))
        s = 0; used = 0; count = 0
        do nup = 0, cell%n
            if (nup > cell%nmodes .or. cell%n - nup > cell%nmodes) cycle
            call jellium_sector_momenta(cell, nup, ks)
            do l = 1, size(ks, 2)
                s = s + 1
                call jellium_sector(cell, nup, ks(:, l:l), spectrum%sectors(s), ok)
                if (.not. ok) return
                d = spectrum%sectors(s)%d
                allocate (h(d, d), e(d), x(d), zero(d), y(d), yi(d), work(8*max(1, 3*d)))
                zero = 0
                do j = 1, d
                    x = 0; x(j) = 1
                    call fock_matvec(model, spectrum%sectors(s)%basis, x, zero, y, yi, ok)
                    if (.not. ok) return
                    h(:, j) = y
                end do
                call dsyev('V', 'U', d, h, d, e, work, size(work), info)
                ok = info == 0
                if (.not. ok) return
                k = count_window(e, window)
                spectrum%kept(s) = k; spectrum%top(s) = e(k)
                need = used + int(d, i64)*k
                call grow(spectrum, count + k, need)
                do j = 1, k
                    count = count + 1
                    spectrum%sector_of(count) = s
                    spectrum%offset(count) = used
                    spectrum%energy(count) = e(j)
                    spectrum%residual(count) = 0
                    spectrum%flat(used + 1:used + d) = h(:, j)
                    used = used + d
                end do
                deallocate (h, e, x, zero, y, yi, work)
            end do
        end do
        spectrum%count = count; spectrum%nsectors = s
        ok = .true.
    end subroutine

    integer function count_window(e, window) result(k)
        real(dp), intent(in) :: e(:), window
        k = 1
        do while (k < size(e))
            if (e(k + 1) - e(1) > window) exit
            k = k + 1
        end do
    end function

    subroutine grow(spectrum, states, flat)
        type(jellium_spectrum_t), intent(inout) :: spectrum
        integer, intent(in) :: states
        integer(i64), intent(in) :: flat
        integer, allocatable :: gi(:)
        integer(i64), allocatable :: go(:)
        real(dp), allocatable :: gr(:)
        integer :: c
        integer(i64) :: f
        c = size(spectrum%energy)
        if (states > c) then
            c = max(states, 2*c)
            allocate (gi(c)); gi(1:size(spectrum%sector_of)) = spectrum%sector_of
            call move_alloc(gi, spectrum%sector_of)
            allocate (go(c)); go(1:size(spectrum%offset)) = spectrum%offset
            call move_alloc(go, spectrum%offset)
            allocate (gr(c)); gr(1:size(spectrum%energy)) = spectrum%energy
            call move_alloc(gr, spectrum%energy)
            allocate (gr(c)); gr(1:size(spectrum%residual)) = spectrum%residual
            call move_alloc(gr, spectrum%residual)
        end if
        f = size(spectrum%flat, kind=i64)
        if (flat > f) then
            f = max(flat, 2*f)
            allocate (gr(f)); gr(1:size(spectrum%flat)) = spectrum%flat
            call move_alloc(gr, spectrum%flat)
        end if
    end subroutine

    !> Diagonal Fourier moments f2 on q(:, k) and f3 on triple(:, k) (indices
    !> into q, triple(3, k) the index of -(Q1+Q2)) from the ensemble 2-RDM
    !> gamma2 (pair storage of kin6d_quantum_fock) and the 3-RDM entries of
    !> every state with weight above threshold*max(weight).
    subroutine jellium_diagonal_moments(cell, spectrum, weights, gamma2, q, triple, &
                                        f2, f3, ok, threshold)
        type(jellium_cell_t), intent(in) :: cell
        type(jellium_spectrum_t), intent(in) :: spectrum
        real(dp), intent(in) :: weights(:), gamma2(:, :)
        integer, intent(in) :: q(:, :), triple(:, :)
        real(dp), intent(out) :: f2(:), f3(:)
        logical, intent(out) :: ok
        real(dp), intent(in), optional :: threshold
        integer, parameter :: perm(3, 6) = reshape([1, 2, 3, 2, 3, 1, 3, 1, 2, &
                                                    2, 1, 3, 1, 3, 2, 3, 2, 1], [3, 6])
        real(dp), parameter :: psign(6) = [1, 1, 1, -1, -1, -1]
        real(dp), allocatable :: t2(:, :, :), t3(:, :, :, :, :, :), zero(:), vr(:), vi(:)
        integer(i64), allocatable :: tk(:), uk(:)
        integer, allocatable :: pi_(:), pj_(:)
        integer :: rr, np, a, b, i, j, s, d, count, e, ps, ss, t(3), u(3), w(3, 3), m
        integer(i64) :: o, m64
        real(dp) :: cut, val
        logical :: spin_ok

        ok = .false.
        rr = 2*maxval(abs(cell%modes))
        m = cell%m; np = m*(m - 1)/2
        allocate (t2(-rr:rr, -rr:rr, -rr:rr), pi_(np), pj_(np))
        t2 = 0
        a = 0
        do i = 1, m - 1
            do j = i + 1, m
                a = a + 1; pi_(a) = i; pj_(a) = j
            end do
        end do
        ! F2(Q) = sum <a+_a a+_b a_d a_c>, Q = k_a - k_c, spins s_a=s_c, s_b=s_d.
        do b = 1, np
            do a = 1, np
                val = gamma2(b, a)
                if (.not. abs(val) > 0) cycle
                call add2(pi_(a), pj_(a), pi_(b), pj_(b), val)
                call add2(pj_(a), pi_(a), pi_(b), pj_(b), -val)
                call add2(pi_(a), pj_(a), pj_(b), pi_(b), -val)
                call add2(pj_(a), pi_(a), pj_(b), pi_(b), val)
            end do
        end do
        do a = 1, size(q, 2)
            f2(a) = t2(q(1, a), q(2, a), q(3, a))/(cell%n*(cell%n - 1.0_dp))
        end do
        f3 = 0
        if (cell%n < 3) then
            ok = .true.; return
        end if
        allocate (t3(-rr:rr, -rr:rr, -rr:rr, -rr:rr, -rr:rr, -rr:rr))
        t3 = 0
        cut = 1.0e-12_dp
        if (present(threshold)) cut = threshold
        cut = cut*maxval(weights(1:spectrum%count))
        m64 = int(m, i64)
        do e = 1, spectrum%count
            if (weights(e) <= cut) cycle
            s = spectrum%sector_of(e); d = spectrum%sectors(s)%d; o = spectrum%offset(e)
            if (allocated(zero)) deallocate (zero)
            allocate (zero(d)); zero = 0
            call fock_rdm3_entries(spectrum%sectors(s)%basis, &
                spectrum%flat(o + 1:o + d), zero, tk, uk, vr, vi, count, ok)
            if (.not. ok) return
            do i = 1, count
                t = decode(tk(i)); u = decode(uk(i))
                do ps = 1, 6
                    do ss = 1, 6
                        spin_ok = .true.
                        do j = 1, 3
                            if (cell%spin(u(perm(j, ps))) /= cell%spin(t(perm(j, ss)))) &
                                spin_ok = .false.
                            w(:, j) = cell%momentum(:, u(perm(j, ps))) - &
                                      cell%momentum(:, t(perm(j, ss)))
                        end do
                        if (.not. spin_ok) cycle
                        t3(w(1, 1), w(2, 1), w(3, 1), w(1, 2), w(2, 2), w(3, 2)) = &
                            t3(w(1, 1), w(2, 1), w(3, 1), w(1, 2), w(2, 2), w(3, 2)) + &
                            weights(e)*psign(ps)*psign(ss)*vr(i)
                    end do
                end do
            end do
        end do
        do a = 1, size(triple, 2)
            i = triple(1, a); j = triple(2, a)
            f3(a) = t3(q(1, i), q(2, i), q(3, i), q(1, j), q(2, j), q(3, j))/ &
                    (cell%n*(cell%n - 1.0_dp)*(cell%n - 2.0_dp))
        end do
        ok = all(ieee_is_finite(f2)) .and. all(ieee_is_finite(f3))
    contains
        subroutine add2(ca, cb, an, bn, v)
            integer, intent(in) :: ca, cb, an, bn
            real(dp), intent(in) :: v
            integer :: qq(3)
            if (cell%spin(ca) /= cell%spin(an) .or. cell%spin(cb) /= cell%spin(bn)) return
            qq = cell%momentum(:, ca) - cell%momentum(:, an)
            t2(qq(1), qq(2), qq(3)) = t2(qq(1), qq(2), qq(3)) + v
        end subroutine
        function decode(key) result(idx)
            integer(i64), intent(in) :: key
            integer :: idx(3)
            idx(3) = int(mod(key, m64)) + 1
            idx(2) = int(mod(key/m64, m64)) + 1
            idx(1) = int(key/(m64*m64)) + 1
        end function
    end subroutine

    !> Ideal fermions (no interaction) in the periodic box of N electrons at
    !> Theta, canonical over all S_z: weights Z_{N_up} Z_{N-N_up} from the Q3
    !> fixed-N_sigma states; shell averages (Gauss--Legendre, 8 points per
    !> shell) of the Q3 sphere-averaged same-spin pair function. Normalization
    !> as the Q4 cell: g_same, g_opp by 2 (N/2)^2, g_sum by N(N-1).
    subroutine ideal_cell_pair(n, theta, length, edges, g_sum, g_same, g_opp, ok)
        integer, intent(in) :: n
        real(dp), intent(in) :: theta, length, edges(:)
        real(dp), intent(out) :: g_sum(:), g_same(:), g_opp(:)
        logical, intent(out) :: ok
        type(fermion_spectrum_t) :: spec
        type(canonical_state_t) :: up, dn
        real(dp), allocatable :: lz(:), gu(:), gd(:), xr(:), wr(:)
        real(dp) :: ef, b, wsum, w, trunc, same, opp, shell, half
        integer :: nup, nb, i, k
        ok = .false.
        nb = size(edges) - 1
        call box_parameters(n, theta, ef, b)
        call box_spectrum(box_cutoff(b, n), spec, ok)
        if (.not. ok) return
        allocate (lz(0:n), xr(8), wr(8), gu(8), gd(8))
        do nup = 0, n
            lz(nup) = 0
            if (nup > 0) then
                call solve_canonical(spec, nup, b, up, ok, .false.); if (.not. ok) return
                lz(nup) = up%log_z
            end if
        end do
        g_sum = 0; g_same = 0; g_opp = 0; wsum = 0
        half = n*n/2.0_dp
        do nup = 0, n
            w = exp(lz(nup) + lz(n - nup) - maxval([(lz(k) + lz(n - k), k=0, n)]))
            wsum = wsum + w
            if (nup > 0) call solve_canonical(spec, nup, b, up, ok, .false.)
            if (n - nup > 0) call solve_canonical(spec, n - nup, b, dn, ok, .false.)
            do i = 1, nb
                call gauss_legendre_ab(8, edges(i), edges(i + 1), xr, wr)
                gu = 0; gd = 0
                if (nup > 1) call same_spin_pair(spec, up, xr/length, gu, ok, 1.0e-16_dp, trunc)
                if (n - nup > 1) call same_spin_pair(spec, dn, xr/length, gd, ok, &
                                                     1.0e-16_dp, trunc)
                shell = sum(wr*xr**2)
                same = sum(wr*xr**2*(nup**2*gu + (n - nup)**2*gd))/shell
                opp = 2.0_dp*nup*(n - nup)
                g_same(i) = g_same(i) + w*same/half
                g_opp(i) = g_opp(i) + w*opp/half
                g_sum(i) = g_sum(i) + w*(same + opp)/(n*(n - 1.0_dp))
            end do
        end do
        g_sum = g_sum/wsum; g_same = g_same/wsum; g_opp = g_opp/wsum
        ok = all(ieee_is_finite(g_sum))
    end subroutine

    !> Volume-weighted RMS over the shells (|r| < L/2): sqrt(sum V_b a_b^2/sum V_b).
    real(dp) function shell_rms(edges, a) result(r)
        real(dp), intent(in) :: edges(:), a(:)
        real(dp) :: v, total
        integer :: i
        r = 0; total = 0
        do i = 1, size(a)
            v = edges(i + 1)**3 - edges(i)**3
            r = r + v*a(i)**2; total = total + v
        end do
        r = sqrt(r/total)
    end function

    !> Spin-blind classical channels in the Q4 normalization 2 (N/2)^2:
    !> independent uniform spin labels give g_same = g_opp = g (N-1)/N.
    subroutine pair_channels_classical(n, g, g_same, g_opp)
        integer, intent(in) :: n
        real(dp), intent(in) :: g(:)
        real(dp), intent(out) :: g_same(:), g_opp(:)
        g_same = g*(n - 1.0_dp)/n; g_opp = g_same
    end subroutine
end module kin6d_quantum_classical_map
