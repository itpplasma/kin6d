module kin6d_quantum_fock
    use, intrinsic :: iso_fortran_env, only: dp => real64, i64 => int64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortnum_interval, only: interval_t, interval, ipoint, entire
    use kin6d_quantum_fock_leaves, only: &
        qf_accumulate, &
        qf_accumulate_interval, &
        qf_accumulate_jvp, &
        qf_accumulate_vjp, &
        qf_add, qf_add_interval, &
        qf_average, &
        qf_coulomb, &
        qf_coulomb_interval, &
        qf_density_contribution, &
        qf_enclose_interval, &
        qf_expect, &
        qf_expect_interval, &
        qf_kinetic, &
        qf_kinetic_interval, &
        qf_norm_accumulate_interval, &
        qf_observable_radius_interval, &
        qf_pair_action_accumulate, &
        qf_pair_action_accumulate_interval, &
        qf_pair_integral, &
        qf_pair_integral_interval
    use kin6d_sort_keys, only: sort_key_order
    implicit none
    private
    public :: fock_basis_t, fock_hamiltonian_t, init_fock_basis, &
              fock_basis_info, &
              get_fock_states, fock_pair_indices, fermion_apply, &
              fock_apply_word, &
              prepare_fock_hamiltonian, prepare_jellium_hamiltonian, &
              assemble_fock_hamiltonian, &
              assemble_fock_hamiltonian_interval, &
              fock_matvec, fock_matvec_interval, &
              fock_matvec_jvp, fock_matvec_vjp, &
              fock_density_rhs, fock_operator_expectation, &
              fock_rdms, fock_rdms_interval, &
              get_fock_integrals, fock_model_info, &
              fock_psi_operator_expectation, &
              fock_pair_action, fock_pair_action_interval, &
              assemble_fock_hamiltonian_block, &
              assemble_fock_hamiltonian_block_interval, &
              assemble_jellium_pair_sector, &
              init_fock_sector_basis, fock_pair_index, add_fock_one_body, &
              fock_matvec_between, fock_matvec_between_interval, &
              fock_transition_rdms, fock_rdm3_entries, fock_rdm2_entries
    type :: fock_basis_t
        private
        integer :: m = 0, n = 0, d = 0, p = 0
        integer(i64) :: bytes = 0
        integer(i64), allocatable :: states(:)
        integer, allocatable :: pi(:), pj(:)
        logical :: valid = .false., sector = .false.
    end type
    type :: fock_hamiltonian_t
        private
        integer :: m = 0, n = 0, d = 0, p = 0
        real(dp), allocatable :: h(:, :), v(:, :)
        type(interval_t), allocatable :: hb(:, :), vb(:, :)
        integer, allocatable :: hptr(:), hcol(:), vptr(:), vcol(:)
        logical :: valid = .false., bounded = .false.
    end type
contains
    subroutine init_fock_basis(basis, m, n, max_states, &
                               max_bytes, ok, status)
        type(fock_basis_t), intent(out) :: basis
        integer, intent(in) :: m, n
        integer(i64), intent(in) :: max_states, max_bytes
        logical, intent(out) :: ok
        character(*), intent(out) :: status
        integer(i64) :: count, cap, mask, low, ripple, bytes
        integer :: k, i, j, p, allocation
        ok = .false.; status = 'invalid fixed-N basis or budget'
        if (m < 1 .or. m > 62 .or. n < 0 .or. n > m) return
        if (max_states < 1 .or. max_bytes < 1) return
        cap = min(max_states, int(huge(1), i64)); count = 1
        do k = 1, min(n, m - n)
            ! count never exceeds default-integer capacity before multiplication;
            ! m<=62 therefore keeps this int64 intermediate below 2^63.
            count = count*int(m - k + 1, i64)/int(k, i64)
            if (count > cap) then
                status = 'Hilbert dimension exceeds explicit sta'// &
                         'te budget'; return
            end if
        end do
        p = m*(m - 1)/2
        bytes = count*8_i64 + int(p, i64)*2*int(storage_size(1)/8, i64)
        if (bytes > max_bytes) then
            status = 'basis payload exceeds explicit byte budget'; return
        end if
        allocate (basis%states(int(count)), basis%pi(p), &
                  basis%pj(p), stat=allocation)
        if (allocation /= 0) then
            status = 'basis allocation failed'; return
        end if
        mask = ishft(1_i64, n) - 1
        do k = 1, int(count)
            basis%states(k) = mask
            if (k == int(count)) exit
            low = iand(mask, -mask); ripple = mask + low
            mask = ior(ripple, ishft(ieor(mask, ripple)/low, -2))
        end do
        k = 0
        do i = 1, m - 1
            do j = i + 1, m
                k = k + 1; basis%pi(k) = i; basis%pj(k) = j
            end do
        end do
        basis%m = m; basis%n = n; basis%d = int(count); basis%p = p; basis%bytes = bytes
        basis%valid = .true.
        ok = .true.
        status = 'exact CAR basis; payload budget met'
    end subroutine
    subroutine fock_basis_info(basis, m, n, dimension, pairs, bytes, ok)
        type(fock_basis_t), intent(in) :: basis
        integer, intent(out) :: m, n, dimension, pairs
        integer(i64), intent(out) :: bytes
        logical, intent(out) :: ok
        m = basis%m; n = basis%n; dimension = basis%d; pairs = basis%p
        bytes = basis%bytes; ok = basis%valid
    end subroutine
    subroutine get_fock_states(basis, states, ok)
        type(fock_basis_t), intent(in) :: basis
        integer(i64), intent(out) :: states(:)
        logical, intent(out) :: ok
        ok = .false.; states = 0
        if (.not. basis%valid) return
        if (size(states) /= basis%d) return
        states = basis%states; ok = .true.
    end subroutine
    subroutine fock_pair_indices(basis, i, j, ok)
        type(fock_basis_t), intent(in) :: basis
        integer, intent(out) :: i(:), j(:)
        logical, intent(out) :: ok
        ok = .false.; i = 0; j = 0
        if (.not. basis%valid) return
        if (size(i) /= basis%p .or. size(j) /= basis%p) return
        i = basis%pi; j = basis%pj; ok = .true.
    end subroutine
    pure subroutine fermion_apply(mask, orbital, creation, &
                                  result, phase, nonzero)
        integer(i64), intent(in) :: mask
        integer, intent(in) :: orbital
        logical, intent(in) :: creation
        integer(i64), intent(out) :: result
        integer, intent(out) :: phase
        logical, intent(out) :: nonzero
        result = mask; phase = 0; nonzero = .false.
        if (mask < 0 .or. mask >= ishft(1_i64, &
                                        62) .or. orbital < 1 .or. orbital > 62) return
        if (btest(mask, orbital - 1) .eqv. creation) return
        phase = 1
        if (mod(popcnt(iand(mask, ishft(1_i64, &
                                        orbital - 1) - 1)), 2) == 1) phase = -1
        result = ieor(mask, ishft(1_i64, orbital - 1)); nonzero = .true.
    end subroutine
    subroutine fock_apply_word(basis, ket, orbitals, creation, &
                               bra, phase, nonzero)
        type(fock_basis_t), intent(in) :: basis
        integer, intent(in) :: ket, orbitals(:)
        logical, intent(in) :: creation(:)
        integer, intent(out) :: bra, phase
        logical, intent(out) :: nonzero
        integer(i64) :: mask, next
        integer :: k, local_phase, lo, hi, mid
        logical :: active
        bra = 0; phase = 0; nonzero = .false.
        if (.not. basis%valid) return
        if (ket < 1 .or. ket > basis%d) return
        if (size(orbitals) /= size(creation)) return
        if (any(orbitals < 1) .or. any(orbitals > basis%m)) return
        mask = basis%states(ket); phase = 1
        do k = size(orbitals), 1, -1
            call fermion_apply(mask, orbitals(k), creation(k), &
                               next, local_phase, active)
            if (.not. active) then
                phase = 0; return
            end if
            mask = next; phase = phase*local_phase
        end do
        lo = 1; hi = basis%d
        do while (lo <= hi)
            mid = lo + (hi - lo)/2
            if (basis%states(mid) == mask) then
                bra = mid; nonzero = .true.; return
            else if (basis%states(mid) < mask) then
                lo = mid + 1
            else
                hi = mid - 1
            end if
        end do
        phase = 0
    end subroutine
    subroutine allocate_model(model, basis, max_bytes, ok, status)
        type(fock_hamiltonian_t), intent(out) :: model
        type(fock_basis_t), intent(in) :: basis
        integer(i64), intent(in) :: max_bytes
        logical, intent(out) :: ok
        character(*), intent(out) :: status
        integer(i64) :: bytes
        integer :: allocation
        ok = .false.; status = 'invalid basis'
        if (.not. basis%valid) return
        bytes = (int(basis%m, i64)**2 + int(basis%p, i64)**2)*24_i64
        if (bytes > max_bytes) then
            status = 'mode-pair Hamiltonian payload exceeds byte'// &
                     ' budget'; return
        end if
        allocate (model%h(basis%m, basis%m), model%v(basis%p, basis%p), &
                  model%hb(basis%m, basis%m), model%vb(basis%p, &
                                                       basis%p), stat=allocation)
        if (allocation /= 0) then
            status = 'Hamiltonian allocation failed'; return
        end if
        model%m = basis%m; model%n = basis%n; model%d = basis%d; model%p = basis%p
        model%h = 0; model%v = 0; model%hb = ipoint(0._dp); model%vb = ipoint(0._dp)
        ok = .true.; status = 'mode-pair payload budget met'
    end subroutine
    subroutine prepare_fock_hamiltonian(model, basis, h, v, &
                                        max_bytes, ok, status)
        type(fock_hamiltonian_t), intent(out) :: model
        type(fock_basis_t), intent(in) :: basis
        real(dp), intent(in) :: h(:, :), v(:, :)
        integer(i64), intent(in) :: max_bytes
        logical, intent(out) :: ok
        character(*), intent(out) :: status
        ok = .false.; status = 'invalid real symmetric one/pair Hamiltonian'
        if (.not. basis%valid) return
        if (any(shape(h) /= [basis%m, basis%m])) return
        if (any(shape(v) /= [basis%p, basis%p])) return
        if (.not. all(ieee_is_finite(h)) .or. .not. all(ieee_is_finite(v))) return
        if (any(h > transpose(h)) .or. any(h < transpose(h))) return
        if (any(v > transpose(v)) .or. any(v < transpose(v))) return
        call allocate_model(model, basis, max_bytes, ok, status)
        if (.not. ok) return
        model%h = h; model%v = v; model%hb = ipoint(h); model%vb = ipoint(v)
        call index_couplings(model, max_bytes, ok, status)
        if (.not. ok) return
        model%valid = .true.; model%bounded = .true.
        status = 'certified-discrete matrix data: exact binary64'// &
                 ' h and pair v'
    end subroutine
    subroutine prepare_jellium_hamiltonian(model, basis, &
                                           momentum, spin, length, &
                                           hbar, mass, coupling, &
                                           max_bytes, ok, status)
        type(fock_hamiltonian_t), intent(out) :: model
        type(fock_basis_t), intent(in) :: basis
        integer, intent(in) :: momentum(:, :), spin(:)
        real(dp), intent(in) :: length, hbar, mass, coupling
        integer(i64), intent(in) :: max_bytes
        logical, intent(out) :: ok
        character(*), intent(out) :: status
        integer :: i, j, k, l, p, q
        integer(i64) :: norm2
        real(dp) :: direct, exchange
        logical :: direct_active, exchange_active
        type(interval_t) :: db, xb
        ok = .false.; status = 'invalid periodic 3D jellium para'// &
 'meters/modes'
        if (.not. basis%valid) return
        if (any(shape(momentum) /= [3, &
                                    basis%m]) .or. size(spin) /= basis%m) return
        if (any(momentum < -1000000) .or. any(momentum > 1000000)) return
        if (any(abs(spin) /= 1)) return
        if (.not. all(ieee_is_finite([length, hbar, mass, &
                                      coupling]))) return
        if (length <= 0 .or. hbar <= 0 .or. mass <= 0 .or. coupling < 0) return
        do i = 1, basis%m
            do j = i + 1, basis%m
                if (all(momentum(:, i) == momentum(:, &
                                                   j)) .and. spin(i) == spin(j)) return
            end do
        end do
        call allocate_model(model, basis, max_bytes, ok, status)
        if (.not. ok) return
        do i = 1, basis%m
            norm2 = sum(int(momentum(:, i), i64)**2)
            call qf_kinetic(length, hbar, mass, real(norm2, dp), &
                            model%h(i, i))
            call qf_kinetic_interval(ipoint(length), &
                                     ipoint(hbar), ipoint(mass), &
                                     ipoint(real(norm2, dp)), &
                                     model%hb(i, i))
        end do
        do p = 1, basis%p
            i = basis%pi(p); j = basis%pj(p)
            do q = p, basis%p
                k = basis%pi(q); l = basis%pj(q)
                if (any(momentum(:, i) + momentum(:, j) /= momentum(:, &
                                                             k) + momentum(:, l))) cycle
                direct = 0; exchange = 0; db = ipoint(0._dp); xb = ipoint(0._dp)
                direct_active = .false.; exchange_active = .false.
                if (spin(i) == spin(k) .and. spin(j) == spin(l)) then
                    norm2 = sum((int(momentum(:, i), &
                                     i64) - int(momentum(:, k), i64))**2)
                    if (norm2 > 0) then
                        direct_active = .true.
                        call qf_coulomb(length, coupling, &
                                        real(norm2, dp), direct)
                        call qf_coulomb_interval(ipoint(length), &
                                                 ipoint(coupling), &
                                                 ipoint(real(norm2, dp)), db)
                    end if
                end if
                if (spin(i) == spin(l) .and. spin(j) == spin(k)) then
                    norm2 = sum((int(momentum(:, i), &
                                     i64) - int(momentum(:, l), i64))**2)
                    if (norm2 > 0) then
                        exchange_active = .true.
                        call qf_coulomb(length, coupling, &
                                        real(norm2, dp), exchange)
                        call qf_coulomb_interval(ipoint(length), &
                                                 ipoint(coupling), &
                                                 ipoint(real(norm2, dp)), xb)
                    end if
                end if
                if (.not. direct_active .and. .not. exchange_active) cycle
                if (coupling <= 0) cycle
                call qf_pair_integral(direct, exchange, model%v(p, q))
                call qf_pair_integral_interval(db, xb, model%vb(p, q))
                model%v(q, p) = model%v(p, q); model%vb(q, &
 p) = model%vb(p, q)
            end do
        end do
        if (.not. all(ieee_is_finite(model%h)) .or. &
            .not. all(ieee_is_finite(model%v))) then
            ok = .false.; status = 'nonfinite candidate Hamiltonian'; return
        end if
        call index_couplings(model, max_bytes, ok, status)
        if (.not. ok) return
        model%valid = .true.
        model%bounded = valid_matrix(model%hb) .and. valid_matrix(model%vb)
        status = 'certified-discrete jellium matrix elements; co'// &
                 'ntinuum cuto'// &
                 'ff unavailable'
        if (.not. model%bounded) status = 'candidate jellium mat'// &
                                          'rix; interval '// &
                                          'construction failed'
    end subroutine
    subroutine assemble_jellium_pair_sector(momentum, length, hbar, mass, &
                                            coupling, max_bytes, hr, hi, hrb, hib, &
                                            ok, status, certify)
        integer, intent(in) :: momentum(:, :)
        real(dp), intent(in) :: length, hbar, mass, coupling
        integer(i64), intent(in) :: max_bytes
        real(dp), intent(out) :: hr(:, :), hi(:, :)
        type(interval_t), intent(out) :: hrb(:, :), hib(:, :)
        logical, intent(out) :: ok
        character(*), intent(out) :: status
        integer :: k, i, j
        integer(i64) :: count, norm2
        real(dp) :: kinetic
        type(interval_t) :: kinetic_b
        logical, intent(in), optional :: certify
        logical :: strict
        strict = .true.
        if (present(certify)) strict = certify
        ok = .false.; status = 'invalid physical opposite-spin pair sector'
        k = size(momentum, 2); count = int(k, i64)
        if (size(momentum, 1) /= 3 .or. k < 1) return
        if (any(shape(hr) /= [k, k]) .or. any(shape(hi) /= [k, k])) return
        if (any(shape(hrb) /= [k, k]) .or. any(shape(hib) /= [k, k])) return
        if (max_bytes < 1_i64) return
        if (count > max_bytes/(48_i64*count)) then
            status = 'pair-sector output payload exceeds byte budget'; return
        end if
        if (any(momentum < -1000000) .or. any(momentum > 1000000)) return
        if (.not. all(ieee_is_finite([length, hbar, mass, coupling]))) return
        if (min(length, hbar, mass) <= 0._dp .or. coupling < 0._dp) return
        do i = 1, k
            do j = 1, i - 1
                if (all(momentum(:, i) == momentum(:, j))) return
            end do
        end do
        hr = 0._dp; hi = 0._dp; hrb = entire(); hib = entire()
        if (strict) then
            hrb = ipoint(0._dp); hib = ipoint(0._dp)
        end if
        do i = 1, k
            norm2 = sum(int(momentum(:, i), i64)**2)
            call qf_kinetic(length, hbar, mass, real(norm2, dp), kinetic)
            call qf_add(kinetic, kinetic, hr(i, i))
            if (strict) then
                call qf_kinetic_interval(ipoint(length), ipoint(hbar), ipoint(mass), &
                                         ipoint(real(norm2, dp)), kinetic_b)
                call qf_add_interval(kinetic_b, kinetic_b, hrb(i, i))
            end if
            if (coupling <= 0._dp) cycle
            do j = 1, i - 1
                norm2 = sum((int(momentum(:, i), i64) - &
                             int(momentum(:, j), i64))**2)
                call qf_coulomb(length, coupling, real(norm2, dp), hr(i, j))
                if (strict) then
                    call qf_coulomb_interval(ipoint(length), ipoint(coupling), &
                                             ipoint(real(norm2, dp)), hrb(i, j))
                end if
                hr(j, i) = hr(i, j); hrb(j, i) = hrb(i, j)
            end do
        end do
        if (.not. all(ieee_is_finite(hr))) then
            status = 'nonfinite pair-sector candidate'; return
        end if
        if (strict .and. .not. valid_matrix(hrb)) then
            status = 'pair-sector strict physical packet failed'; return
        end if
        ok = .true.
        status = 'candidate exact spin/momentum sector; strict packet unavailable'
        if (strict) status = &
            'certified-discrete exact spin/momentum sector; cutoff unavailable'
    end subroutine
    logical function compatible(model, basis)
        type(fock_hamiltonian_t), intent(in) :: model
        type(fock_basis_t), intent(in) :: basis
        compatible = .false.
        if (.not. model%valid .or. .not. basis%valid) return
        ! Full-scan schedules silently drop out-of-basis targets; they are
        ! restricted to complete fixed-N bases.
        if (basis%sector) return
        compatible = model%m == basis%m .and. model%n == basis%n .and. &
                     model%d == basis%d
    end function
    subroutine next_transition(model, basis, cursor, ket, bra, &
                               phase, sector, ci, cj, more)
        type(fock_hamiltonian_t), intent(in) :: model
        type(fock_basis_t), intent(in) :: basis
        integer, intent(inout) :: cursor(4)
        integer, intent(out) :: ket, bra, phase, sector, ci, cj
        logical, intent(out) :: more
        integer :: i, j, p, q
        integer(i64) :: mask
        logical :: active
        more = .false.; ket = 0; bra = 0; phase = 0; sector = 0; ci = 0; cj = 0
        do while (cursor(1) <= model%d)
            ket = cursor(1); mask = basis%states(ket)
            if (cursor(2) == 1) then
                do while (cursor(3) <= model%m)
                    j = cursor(3)
                    if (.not. btest(mask, j - 1)) then
                        cursor(3) = cursor(3) + 1; cursor(4) = 1; cycle
                    end if
                    do while (cursor(4) <= model%m)
                        i = cursor(4); cursor(4) = cursor(4) + 1
                        call fock_apply_word(basis, ket, [i, j], &
                                             [.true., .false.], &
                                             bra, phase, active)
                        if (active) then
                            sector = 1; ci = i; cj = j; more = .true.; return
                        end if
                    end do
                    cursor(3) = cursor(3) + 1; cursor(4) = 1
                end do
                cursor(2) = 2; cursor(3) = 1; cursor(4) = 1
            end if
            do while (cursor(3) <= model%p)
                q = cursor(3)
                if (.not. btest(mask, basis%pi(q) - 1) .or. &
                    .not. btest(mask, basis%pj(q) - 1)) then
                    cursor(3) = cursor(3) + 1; cursor(4) = 1; cycle
                end if
                do while (cursor(4) <= model%p)
                    p = cursor(4); cursor(4) = cursor(4) + 1
                    call fock_apply_word(basis, ket, &
                                         [basis%pi(p), basis%pj(p), &
                                          basis%pj(q), basis%pi(q)], &
                                         [.true., .true., .false., &
                                          .false.], bra, phase, active)
                    if (active) then
                        sector = 2; ci = p; cj = q; more = .true.; return
                    end if
                end do
                cursor(3) = cursor(3) + 1; cursor(4) = 1
            end do
            cursor = [cursor(1) + 1, 1, 1, 1]
        end do
    end subroutine
    subroutine transition_coefficient(model, sector, i, j, coef, bound)
        type(fock_hamiltonian_t), intent(in) :: model
        integer, intent(in) :: sector, i, j
        real(dp), intent(out) :: coef
        type(interval_t), intent(out) :: bound
        if (sector == 1) then
            coef = model%h(i, j); bound = model%hb(i, j)
        else
            coef = model%v(i, j); bound = model%vb(i, j)
        end if
    end subroutine
    subroutine assemble_fock_hamiltonian(model, basis, hre, him, ok)
        type(fock_hamiltonian_t), intent(in) :: model
        type(fock_basis_t), intent(in) :: basis
        real(dp), intent(out) :: hre(:, :), him(:, :)
        logical, intent(out) :: ok
        integer :: cursor(4), ket, bra, phase, sector, i, j
        real(dp) :: coef, next
        type(interval_t) :: bound
        logical :: more
        ok = .false.; hre = 0; him = 0
        if (.not. compatible(model, basis)) return
        if (model%d > 46340) return
        if (any(shape(hre) /= [model%d, model%d])) return
        if (any(shape(him) /= [model%d, model%d])) return
        cursor = [1, 1, 1, 1]
        do
            call next_transition(model, basis, cursor, ket, bra, &
                                 phase, sector, i, j, more)
            if (.not. more) exit
            call transition_coefficient(model, sector, i, j, coef, bound)
            call qf_accumulate(hre(bra, ket), coef, 1._dp, &
                               real(phase, dp), next)
            hre(bra, ket) = next
        end do
        do i = 1, model%d
            do j = i + 1, model%d
                call qf_average(hre(i, j), hre(j, i), next)
                hre(i, j) = next; hre(j, i) = next
            end do
        end do
        ok = all(ieee_is_finite(hre))
    end subroutine
    subroutine assemble_fock_hamiltonian_interval(model, basis, &
                                                  hre, him, ok)
        type(fock_hamiltonian_t), intent(in) :: model
        type(fock_basis_t), intent(in) :: basis
        type(interval_t), intent(out) :: hre(:, :), him(:, :)
        logical, intent(out) :: ok
        integer :: cursor(4), ket, bra, phase, sector, i, j
        real(dp) :: coef
        type(interval_t) :: bound, next
        logical :: more
        ok = .false.; hre = ipoint(0._dp); him = ipoint(0._dp)
        if (.not. compatible(model, basis)) return
        if (.not. model%bounded) return
        if (model%d > 46340) return
        if (any(shape(hre) /= [model%d, model%d])) return
        if (any(shape(him) /= [model%d, model%d])) return
        cursor = [1, 1, 1, 1]
        do
            call next_transition(model, basis, cursor, ket, bra, &
                                 phase, sector, i, j, more)
            if (.not. more) exit
            call transition_coefficient(model, sector, i, j, coef, bound)
            call qf_accumulate_interval(hre(bra, ket), bound, &
                                        ipoint(1._dp), &
                                        ipoint(real(phase, dp)), next)
            hre(bra, ket) = next
        end do
        ok = valid_matrix(hre)
    end subroutine
    logical function valid_selection(model, basis, indices)
        type(fock_hamiltonian_t), intent(in) :: model
        type(fock_basis_t), intent(in) :: basis
        integer, intent(in) :: indices(:)
        integer :: i, j
        valid_selection = .false.
        if (.not. compatible(model, basis)) return
        if (size(indices) < 1 .or. size(indices) > 46340) return
        if (any(indices < 1) .or. any(indices > model%d)) return
        do i = 1, size(indices)
            do j = i + 1, size(indices)
                if (indices(i) == indices(j)) return
            end do
        end do
        valid_selection = .true.
    end function
    integer function selected_index(indices, bra) result(index)
        integer, intent(in) :: indices(:), bra
        integer :: i
        index = 0
        do i = 1, size(indices)
            if (indices(i) == bra) then
                index = i; return
            end if
        end do
    end function
    subroutine assemble_fock_hamiltonian_block(model, basis, &
                                               indices, hre, him, &
                                               ok, invariant)
        type(fock_hamiltonian_t), intent(in) :: model
        type(fock_basis_t), intent(in) :: basis
        integer, intent(in) :: indices(:)
        real(dp), intent(out) :: hre(:, :), him(:, :)
        logical, intent(out) :: ok, invariant
        integer :: cursor(4), ket, bra, phase, sector, i, j, col, row, k
        real(dp) :: coef, next
        type(interval_t) :: bound
        logical :: more
        ok = .false.; invariant = .false.; hre = 0; him = 0
        if (.not. valid_selection(model, basis, indices)) return
        k = size(indices)
        if (any(shape(hre) /= [k, k]) .or. any(shape(him) /= [k, k])) return
        invariant = model%bounded
        do col = 1, k
            cursor = [indices(col), 1, 1, 1]
            do
                call next_transition(model, basis, cursor, ket, &
                                     bra, phase, &
                                     sector, i, j, more)
                if (.not. more) exit
                if (ket /= indices(col)) exit
                call transition_coefficient(model, sector, i, j, &
                                            coef, bound)
                row = selected_index(indices, bra)
                if (row < 1) then
                    if (bound%lo < 0 .or. bound%hi > 0) invariant = .false.
                    cycle
                end if
                call qf_accumulate(hre(row, col), coef, 1._dp, &
                                   real(phase, dp), next)
                hre(row, col) = next
            end do
        end do
        do i = 1, k
            do j = i + 1, k
                call qf_average(hre(i, j), hre(j, i), next)
                hre(i, j) = next; hre(j, i) = next
            end do
        end do
        ok = all(ieee_is_finite(hre))
    end subroutine
    subroutine assemble_fock_hamiltonian_block_interval(model, &
                                                        basis, indices, &
                                                        hre, him, ok, invariant)
        type(fock_hamiltonian_t), intent(in) :: model
        type(fock_basis_t), intent(in) :: basis
        integer, intent(in) :: indices(:)
        type(interval_t), intent(out) :: hre(:, :), him(:, :)
        logical, intent(out) :: ok, invariant
        integer :: cursor(4), ket, bra, phase, sector, i, j, col, row, k
        real(dp) :: coef
        type(interval_t) :: bound, next
        logical :: more
        ok = .false.; invariant = .false.; hre = ipoint(0._dp); him = ipoint(0._dp)
        if (.not. valid_selection(model, basis, indices)) return
        if (.not. model%bounded) return
        k = size(indices)
        if (any(shape(hre) /= [k, k]) .or. any(shape(him) /= [k, k])) return
        invariant = .true.
        do col = 1, k
            cursor = [indices(col), 1, 1, 1]
            do
                call next_transition(model, basis, cursor, ket, &
                                     bra, phase, &
                                     sector, i, j, more)
                if (.not. more) exit
                if (ket /= indices(col)) exit
                call transition_coefficient(model, sector, i, j, &
                                            coef, bound)
                row = selected_index(indices, bra)
                if (row < 1) then
                    if (bound%lo < 0 .or. bound%hi > 0) invariant = .false.
                    cycle
                end if
                call qf_accumulate_interval(hre(row, col), &
                                            bound, ipoint(1._dp), &
                                            ipoint(real(phase, dp)), next)
                hre(row, col) = next
            end do
        end do
        ok = valid_matrix(hre)
    end subroutine
    subroutine fock_matvec(model, basis, re, im, yr, yi, ok)
        type(fock_hamiltonian_t), intent(in) :: model
        type(fock_basis_t), intent(in) :: basis
        real(dp), intent(in) :: re(:), im(:)
        real(dp), intent(out) :: yr(:), yi(:)
        logical, intent(out) :: ok
        ok = .false.; yr = 0; yi = 0
        if (.not. model%valid .or. .not. basis%valid) return
        if (model%d /= basis%d .and. .not. basis%sector) return
        call fock_matvec_between(model, basis, basis, re, im, yr, yi, ok)
    end subroutine
    subroutine fock_matvec_interval(model, basis, re, im, yr, yi, ok)
        type(fock_hamiltonian_t), intent(in) :: model
        type(fock_basis_t), intent(in) :: basis
        type(interval_t), intent(in) :: re(:), im(:)
        type(interval_t), intent(out) :: yr(:), yi(:)
        logical, intent(out) :: ok
        ok = .false.; yr = ipoint(0._dp); yi = ipoint(0._dp)
        if (.not. model%valid .or. .not. basis%valid) return
        if (model%d /= basis%d .and. .not. basis%sector) return
        call fock_matvec_between_interval(model, basis, basis, re, im, yr, yi, ok)
    end subroutine
    subroutine fock_matvec_jvp(model, basis, re, im, dh, dv, &
                               dre, dim, dyr, dyi, ok)
        type(fock_hamiltonian_t), intent(in) :: model
        type(fock_basis_t), intent(in) :: basis
        real(dp), intent(in) :: re(:), im(:), dh(:, :), dv(:, &
                                                           :), dre(:), dim(:)
        real(dp), intent(out) :: dyr(:), dyi(:)
        logical, intent(out) :: ok
        integer :: cursor(4), ket, bra, phase, sector, i, j
        real(dp) :: coef, dcoef, next
        type(interval_t) :: bound
        logical :: more
        ok = .false.; dyr = 0; dyi = 0
        if (.not. compatible(model, basis)) return
        if (any([size(re), size(im), size(dre), size(dim), size(dyr), &
                 size(dyi)] /= model%d)) return
        if (any(shape(dh) /= [model%m, &
                              model%m]) .or. any(shape(dv) /= [model%p, &
                                                               model%p])) return
        if (.not. all(ieee_is_finite(re)) .or. .not. all(ieee_is_finite(im))) return
        if (.not. all(ieee_is_finite(dre)) .or. .not. all(ieee_is_finite(dim))) return
        if (.not. all(ieee_is_finite(dh)) .or. .not. all(ieee_is_finite(dv))) return
        if (any(dh > transpose(dh)) .or. any(dh < transpose(dh))) return
        if (any(dv > transpose(dv)) .or. any(dv < transpose(dv))) return
        cursor = [1, 1, 1, 1]
        do
            call next_transition(model, basis, cursor, ket, bra, &
                                 phase, sector, i, j, more)
            if (.not. more) exit
            call transition_coefficient(model, sector, i, j, coef, bound)
            if (sector == 1) then
                dcoef = dh(i, j)
            else
                dcoef = dv(i, j)
            end if
            call qf_accumulate_jvp(0._dp, coef, re(ket), real(phase, dp), &
                                   dyr(bra), dcoef, dre(ket), 0._dp, next)
            dyr(bra) = next
            call qf_accumulate_jvp(0._dp, coef, im(ket), real(phase, dp), &
                                   dyi(bra), dcoef, dim(ket), 0._dp, next)
            dyi(bra) = next
        end do
        ok = all(ieee_is_finite(dyr)) .and. all(ieee_is_finite(dyi))
    end subroutine
    subroutine fock_matvec_vjp(model, basis, re, im, byr, byi, &
                               gh, gv, gre, gim, ok)
        type(fock_hamiltonian_t), intent(in) :: model
        type(fock_basis_t), intent(in) :: basis
        real(dp), intent(in) :: re(:), im(:), byr(:), byi(:)
        real(dp), intent(out) :: gh(:, :), gv(:, :), gre(:), gim(:)
        logical, intent(out) :: ok
        integer :: cursor(4), ket, bra, phase, sector, i, j
        real(dp) :: coef, gacc, gcoef, gvalue, gphase, total, next, average
        type(interval_t) :: bound
        logical :: more
        ok = .false.; gh = 0; gv = 0; gre = 0; gim = 0
        if (.not. compatible(model, basis)) return
        if (any([size(re), size(im), size(byr), size(byi), size(gre), &
                 size(gim)] /= model%d)) return
        if (any(shape(gh) /= [model%m, &
                              model%m]) .or. any(shape(gv) /= [model%p, &
                                                               model%p])) return
        if (.not. all(ieee_is_finite(re)) .or. .not. all(ieee_is_finite(im))) return
        if (.not. all(ieee_is_finite(byr)) .or. .not. all(ieee_is_finite(byi))) return
        cursor = [1, 1, 1, 1]
        do
            call next_transition(model, basis, cursor, ket, bra, &
                                 phase, sector, i, j, more)
            if (.not. more) exit
            call transition_coefficient(model, sector, i, j, coef, bound)
            call qf_accumulate_vjp(0._dp, coef, re(ket), &
                                   real(phase, dp), byr(bra), &
                                   gacc, gcoef, gvalue, gphase)
            total = gcoef; call qf_add(gre(ket), gvalue, &
 next); gre(ket) = next
            call qf_accumulate_vjp(0._dp, coef, im(ket), &
                                   real(phase, dp), byi(bra), &
                                   gacc, gcoef, gvalue, gphase)
            call qf_add(total, gcoef, next); total = next
            call qf_add(gim(ket), gvalue, next); gim(ket) = next
            if (sector == 1) then
                call qf_add(gh(i, j), total, next); gh(i, j) = next
            else
                call qf_add(gv(i, j), total, next); gv(i, j) = next
            end if
        end do
        ! Orthogonal projection gives the Frobenius Riesz gradient on symmetric
        ! real input matrices. No continuous physical adjoint is inferred.
        do i = 1, model%m
            do j = i + 1, model%m
                call qf_average(gh(i, j), gh(j, i), average); gh(i, &
 j) = average; gh(j, i) = average
            end do
        end do
        do i = 1, model%p
            do j = i + 1, model%p
                call qf_average(gv(i, j), gv(j, i), average); gv(i, &
 j) = average; gv(j, i) = average
            end do
        end do
        ok = all(ieee_is_finite(gh)) .and. all(ieee_is_finite(gv)) .and. &
             all(ieee_is_finite(gre)) .and. all(ieee_is_finite(gim))
    end subroutine
    subroutine fock_density_rhs(model, basis, re, im, hbar, dr, di, ok)
        type(fock_hamiltonian_t), intent(in) :: model
        type(fock_basis_t), intent(in) :: basis
        real(dp), intent(in) :: re(:, :), im(:, :), hbar
        real(dp), intent(out) :: dr(:, :), di(:, :)
        logical, intent(out) :: ok
        integer :: cursor(4), ket, bra, phase, sector, i, j, k
        real(dp) :: coef, nr, ni
        type(interval_t) :: bound
        logical :: more
        ok = .false.; dr = 0; di = 0
        if (.not. compatible(model, basis)) return
        if (.not. ieee_is_finite(hbar)) return
        if (hbar <= 0) return
        if (any(shape(re) /= [model%d, &
                              model%d]) .or. any(shape(im) /= [model%d, &
                                                               model%d])) return
        if (any(shape(dr) /= [model%d, &
                              model%d]) .or. any(shape(di) /= [model%d, &
                                                               model%d])) return
        if (.not. all(ieee_is_finite(re)) .or. .not. all(ieee_is_finite(im))) return
        cursor = [1, 1, 1, 1]
        do
            call next_transition(model, basis, cursor, ket, bra, &
                                 phase, sector, i, j, more)
            if (.not. more) exit
            call transition_coefficient(model, sector, i, j, coef, bound)
            do k = 1, model%d
                call qf_density_contribution(dr(bra, k), di(bra, k), coef, &
                                             real(phase, dp), &
                                             re(ket, k), im(ket, &
                                                            k), hbar, nr, ni)
                dr(bra, k) = nr; di(bra, k) = ni
                call qf_density_contribution(dr(k, ket), di(k, ket), coef, &
                                             real(-phase, dp), &
                                             re(k, bra), im(k, &
                                                            bra), hbar, nr, ni)
                dr(k, ket) = nr; di(k, ket) = ni
            end do
        end do
        ok = all(ieee_is_finite(dr)) .and. all(ieee_is_finite(di))
    end subroutine
    subroutine fock_operator_expectation(basis, re, im, &
                                         orbitals, creation, vr, vi, ok)
        type(fock_basis_t), intent(in) :: basis
        real(dp), intent(in) :: re(:, :), im(:, :)
        integer, intent(in) :: orbitals(:)
        logical, intent(in) :: creation(:)
        real(dp), intent(out) :: vr, vi
        logical, intent(out) :: ok
        integer :: ket, bra, phase
        real(dp) :: next
        logical :: active
        ok = .false.; vr = 0; vi = 0
        if (.not. basis%valid) return
        if (any(shape(re) /= [basis%d, &
                              basis%d]) .or. any(shape(im) /= [basis%d, &
                                                               basis%d])) return
        if (size(orbitals) /= size(creation)) return
        if (any(orbitals < 1) .or. any(orbitals > basis%m)) return
        if (.not. all(ieee_is_finite(re)) .or. .not. all(ieee_is_finite(im))) return
        do ket = 1, basis%d
            call fock_apply_word(basis, ket, orbitals, creation, &
                                 bra, phase, active)
            if (.not. active) cycle
            call qf_accumulate(vr, re(ket, bra), 1._dp, &
                               real(phase, dp), next)
            vr = next
            call qf_accumulate(vi, im(ket, bra), 1._dp, &
                               real(phase, dp), next)
            vi = next
        end do
        ok = all(ieee_is_finite([vr, vi]))
    end subroutine
    subroutine psi_expectation(basis, re, im, orbitals, creation, vr, vi)
        type(fock_basis_t), intent(in) :: basis
        real(dp), intent(in) :: re(:), im(:)
        integer, intent(in) :: orbitals(:)
        logical, intent(in) :: creation(:)
        real(dp), intent(out) :: vr, vi
        integer :: ket, bra, phase
        real(dp) :: nr, ni
        logical :: active
        vr = 0; vi = 0
        do ket = 1, basis%d
            call fock_apply_word(basis, ket, orbitals, creation, &
                                 bra, phase, active)
            if (.not. active) cycle
            call qf_expect(vr, vi, re(bra), im(bra), re(ket), &
                           im(ket), real(phase, &
                                         dp), nr, ni)
            vr = nr; vi = ni
        end do
    end subroutine
    subroutine fock_rdms(basis, re, im, gr, gi, g2r, g2i, ok)
        type(fock_basis_t), intent(in) :: basis
        real(dp), intent(in) :: re(:), im(:)
        real(dp), intent(out) :: gr(:, :), gi(:, :), g2r(:, :), g2i(:, :)
        logical, intent(out) :: ok
        call fock_transition_rdms(basis, re, im, re, im, gr, gi, g2r, g2i, ok)
    end subroutine
    subroutine psi_expectation_interval(basis, re, im, orbitals, &
                                        creation, budget, vr, vi)
        type(fock_basis_t), intent(in) :: basis
        real(dp), intent(in) :: re(:), im(:)
        integer, intent(in) :: orbitals(:)
        logical, intent(in) :: creation(:)
        type(interval_t), intent(in) :: budget
        type(interval_t), intent(out) :: vr, vi
        integer :: ket, bra, phase
        type(interval_t) :: nr, ni, lower, upper
        logical :: active
        vr = ipoint(0._dp); vi = ipoint(0._dp)
        do ket = 1, basis%d
            call fock_apply_word(basis, ket, orbitals, creation, &
                                 bra, phase, active)
            if (.not. active) cycle
            call qf_expect_interval(vr, vi, ipoint(re(bra)), &
                                    ipoint(im(bra)), &
                                    ipoint(re(ket)), ipoint(im(ket)), &
                                    ipoint(real(phase, dp)), nr, ni)
            vr = nr; vi = ni
        end do
        call qf_enclose_interval(vr, budget, lower, &
                                 upper); vr = interval(lower%lo, upper%hi)
        call qf_enclose_interval(vi, budget, lower, &
                                 upper); vi = interval(lower%lo, upper%hi)
    end subroutine
    subroutine fock_rdms_interval(basis, re, im, state_radius, &
                                  gr, gi, g2r, g2i, ok)
        type(fock_basis_t), intent(in) :: basis
        real(dp), intent(in) :: re(:), im(:), state_radius
        type(interval_t), intent(out) :: gr(:, :), gi(:, :), &
                                         g2r(:, :), g2i(:, :)
        logical, intent(out) :: ok
        integer :: i, j, p, q, k
        type(interval_t) :: norm2, next, budget
        ok = .false.; gr = ipoint(0._dp); gi = ipoint(0._dp)
        g2r = ipoint(0._dp); g2i = ipoint(0._dp)
        if (.not. basis%valid) return
        if (size(re) /= basis%d .or. size(im) /= basis%d) return
        if (any(shape(gr) /= [basis%m, &
                              basis%m]) .or. any(shape(gi) /= [basis%m, &
                                                               basis%m])) return
        if (any(shape(g2r) /= [basis%p, &
                               basis%p]) .or. any(shape(g2i) /= [basis%p, &
                                                                 basis%p])) return
        if (.not. all(ieee_is_finite(re)) .or. .not. all(ieee_is_finite(im))) return
        if (.not. ieee_is_finite(state_radius)) return
        if (state_radius < 0) return
        norm2 = ipoint(0._dp)
        do k = 1, basis%d
            call qf_norm_accumulate_interval(norm2, &
                                             ipoint(re(k)), ipoint(im(k)), next)
            norm2 = next
        end do
        norm2%lo = max(0._dp, norm2%lo)
        call qf_observable_radius_interval(norm2, &
                                           ipoint(state_radius), budget)
        do i = 1, basis%m
            do j = 1, basis%m
                call psi_expectation_interval(basis, re, im, [j, &
                                                              i], [.true., .false.], &
                                              budget, gr(i, j), gi(i, j))
            end do
        end do
        do p = 1, basis%p
            do q = 1, basis%p
                call psi_expectation_interval(basis, re, im, &
                                              [basis%pi(q), basis%pj(q), &
                                               basis%pj(p), basis%pi(p)], &
                                              [.true., .true., .false., .false.], &
                                              budget, g2r(p, q), g2i(p, q))
            end do
        end do
        ok = valid_matrix(gr) .and. valid_matrix(gi) .and. &
             valid_matrix(g2r) .and. valid_matrix(g2i)
    end subroutine
    subroutine get_fock_integrals(model, h, v, ok, hb, vb)
        type(fock_hamiltonian_t), intent(in) :: model
        real(dp), intent(out) :: h(:, :), v(:, :)
        logical, intent(out) :: ok
        type(interval_t), intent(out), optional :: hb(:, :), vb(:, :)
        ok = .false.; h = 0; v = 0
        if (.not. model%valid) return
        if (any(shape(h) /= [model%m, &
                             model%m]) .or. any(shape(v) /= [model%p, &
                                                             model%p])) return
        if (present(hb)) then
            if (any(shape(hb) /= [model%m, model%m])) return
            if (.not. model%bounded) return
            hb = model%hb
        end if
        if (present(vb)) then
            if (any(shape(vb) /= [model%p, model%p])) return
            if (.not. model%bounded) return
            vb = model%vb
        end if
        h = model%h; v = model%v; ok = .true.
    end subroutine
    subroutine fock_model_info(model, bytes, bounded, ok)
        type(fock_hamiltonian_t), intent(in) :: model
        integer(i64), intent(out) :: bytes
        logical, intent(out) :: bounded, ok
        bytes = (int(model%m, i64)**2 + int(model%p, i64)**2)*24_i64
        bounded = model%bounded; ok = model%valid
    end subroutine
    subroutine fock_psi_operator_expectation(basis, re, im, &
                                             orbitals, creation, vr, vi, ok)
        type(fock_basis_t), intent(in) :: basis
        real(dp), intent(in) :: re(:), im(:)
        integer, intent(in) :: orbitals(:)
        logical, intent(in) :: creation(:)
        real(dp), intent(out) :: vr, vi
        logical, intent(out) :: ok
        ok = .false.; vr = 0; vi = 0
        if (.not. basis%valid) return
        if (size(re) /= basis%d .or. size(im) /= basis%d) return
        if (size(orbitals) /= size(creation)) return
        if (any(orbitals < 1) .or. any(orbitals > basis%m)) return
        if (.not. all(ieee_is_finite(re)) .or. .not. all(ieee_is_finite(im))) return
        call psi_expectation(basis, re, im, orbitals, creation, vr, vi)
        ok = all(ieee_is_finite([vr, vi]))
    end subroutine
    integer function state_index(basis, mask) result(index)
        type(fock_basis_t), intent(in) :: basis
        integer(i64), intent(in) :: mask
        integer :: lo, hi, mid
        index = 0; lo = 1; hi = basis%d
        do while (lo <= hi)
            mid = lo + (hi - lo)/2
            if (basis%states(mid) == mask) then
                index = mid; return
            else if (basis%states(mid) < mask) then
                lo = mid + 1
            else
                hi = mid - 1
            end if
        end do
    end function
    subroutine fock_pair_action(basis, remaining, fr, fi, sr, si, re, &
                                im, out_re, out_im, ok)
        type(fock_basis_t), intent(in) :: basis, remaining
        real(dp), intent(in) :: fr(:), fi(:), sr(:), si(:), re(:), im(:)
        real(dp), intent(out) :: out_re(:), out_im(:)
        logical, intent(out) :: ok
        integer :: ket, i, j, phase1, phase2, bra
        integer(i64) :: first, second
        real(dp) :: nr, ni
        logical :: active
        ok = .false.; out_re = 0; out_im = 0
        if (.not. basis%valid .or. .not. remaining%valid) return
        if (remaining%m /= basis%m .or. remaining%n /= basis%n - 2) return
        if (any([size(fr), size(fi), size(sr), size(si)] /= basis%m)) return
        if (size(re) /= basis%d .or. size(im) /= basis%d) return
        if (size(out_re) /= remaining%d .or. size(out_im) /= remaining%d) return
        if (.not. all(ieee_is_finite(fr)) .or. .not. all(ieee_is_finite(fi))) return
        if (.not. all(ieee_is_finite(sr)) .or. .not. all(ieee_is_finite(si))) return
        if (.not. all(ieee_is_finite(re)) .or. .not. all(ieee_is_finite(im))) return
        do ket = 1, basis%d
            if (abs(re(ket)) <= 0 .and. abs(im(ket)) <= 0) cycle
            do i = 1, basis%m
                if (.not. btest(basis%states(ket), i - 1)) cycle
                call fermion_apply(basis%states(ket), i, &
                                   .false., first, phase1, active)
                do j = 1, basis%m
                    if (.not. btest(first, j - 1)) cycle
                    call fermion_apply(first, j, .false., &
                                       second, phase2, active)
                    bra = state_index(remaining, second)
                    if (bra < 1) return
                    call qf_pair_action_accumulate(out_re(bra), &
                                                   out_im(bra), &
                                                   fr(i), fi(i), sr(j), &
                                                   si(j), re(ket), im(ket), &
                                                   real(phase1*phase2, dp), nr, ni)
                    out_re(bra) = nr; out_im(bra) = ni
                end do
            end do
        end do
        ok = all(ieee_is_finite(out_re)) .and. all(ieee_is_finite(out_im))
    end subroutine
    subroutine fock_pair_action_interval(basis, remaining, fr, &
                                         fi, sr, si, re, im, &
                                         out_re, out_im, ok)
        type(fock_basis_t), intent(in) :: basis, remaining
        type(interval_t), intent(in) :: fr(:), fi(:), sr(:), &
                                        si(:), re(:), im(:)
        type(interval_t), intent(out) :: out_re(:), out_im(:)
        logical, intent(out) :: ok
        integer :: ket, i, j, phase1, phase2, bra
        integer(i64) :: first, second
        type(interval_t) :: nr, ni
        logical :: active
        ok = .false.; out_re = ipoint(0._dp); out_im = ipoint(0._dp)
        if (.not. basis%valid .or. .not. remaining%valid) return
        if (remaining%m /= basis%m .or. remaining%n /= basis%n - 2) return
        if (any([size(fr), size(fi), size(sr), size(si)] /= basis%m)) return
        if (size(re) /= basis%d .or. size(im) /= basis%d) return
        if (size(out_re) /= remaining%d .or. size(out_im) /= remaining%d) return
        if (.not. valid_intervals(fr) .or. .not. valid_intervals(fi)) return
        if (.not. valid_intervals(sr) .or. .not. valid_intervals(si)) return
        if (.not. valid_intervals(re) .or. .not. valid_intervals(im)) return
        do ket = 1, basis%d
            if (re(ket)%lo >= 0 .and. re(ket)%hi <= 0 .and. &
                im(ket)%lo >= 0 .and. im(ket)%hi <= 0) cycle
            do i = 1, basis%m
                if (.not. btest(basis%states(ket), i - 1)) cycle
                call fermion_apply(basis%states(ket), i, &
                                   .false., first, phase1, active)
                do j = 1, basis%m
                    if (.not. btest(first, j - 1)) cycle
                    call fermion_apply(first, j, .false., &
                                       second, phase2, active)
                    bra = state_index(remaining, second)
                    if (bra < 1) return
                    call qf_pair_action_accumulate_interval(out_re(bra), &
                                                            out_im(bra), &
                                                            fr(i), fi(i), sr(j), &
                                                            si(j), re(ket), im(ket), &
                                                ipoint(real(phase1*phase2, dp)), nr, ni)
                    out_re(bra) = nr; out_im(bra) = ni
                end do
            end do
        end do
        ok = valid_intervals(out_re) .and. valid_intervals(out_im)
    end subroutine
    pure logical function valid_intervals(values)
        type(interval_t), intent(in) :: values(:)
        valid_intervals = all(ieee_is_finite(values%lo)) .and. &
                        all(ieee_is_finite(values%hi)) .and. all(values%lo <= values%hi)
    end function
    pure logical function valid_matrix(values)
        type(interval_t), intent(in) :: values(:, :)
        integer :: i, j
        valid_matrix = .false.
        do j = 1, size(values, 2)
            do i = 1, size(values, 1)
                if (.not. ieee_is_finite(values(i, j)%lo)) return
                if (.not. ieee_is_finite(values(i, j)%hi)) return
                if (values(i, j)%lo > values(i, j)%hi) return
            end do
        end do
        valid_matrix = .true.
    end function
    ! ------------------------------------------------------------------
    ! Symmetry-sector bases, coupled-transition actions and factorized RDMs.
    ! ------------------------------------------------------------------

    !> Basis from a caller-supplied strictly increasing list of N-particle
    !> masks (a symmetry sector or a union of sectors). The coupled actions
    !> below fail closed if a nonzero Hamiltonian coefficient leaves the
    !> basis; the legacy full-scan schedules reject sector bases.
    subroutine init_fock_sector_basis(basis, m, n, states, max_bytes, ok, status)
        type(fock_basis_t), intent(out) :: basis
        integer, intent(in) :: m, n
        integer(i64), intent(in) :: states(:), max_bytes
        logical, intent(out) :: ok
        character(*), intent(out) :: status
        integer(i64) :: bytes
        integer :: k, i, j, p, allocation

        ok = .false.; status = 'invalid fixed-N sector basis or budget'
        if (m < 1 .or. m > 62 .or. n < 0 .or. n > m) return
        if (size(states) < 1 .or. max_bytes < 1) return
        if (any(states < 0) .or. any(states >= ishft(1_i64, m))) return
        if (any(popcnt(states) /= n)) return
        do k = 2, size(states)
            if (states(k) <= states(k - 1)) return
        end do
        p = m*(m - 1)/2
        bytes = size(states, kind=i64)*8_i64 + int(p, i64)*2*int(storage_size(1)/8, i64)
        if (bytes > max_bytes) then
            status = 'basis payload exceeds explicit byte budget'; return
        end if
        allocate (basis%states(size(states)), basis%pi(p), basis%pj(p), stat=allocation)
        if (allocation /= 0) then
            status = 'basis allocation failed'; return
        end if
        basis%states = states
        k = 0
        do i = 1, m - 1
            do j = i + 1, m
                k = k + 1; basis%pi(k) = i; basis%pj(k) = j
            end do
        end do
        basis%m = m; basis%n = n; basis%d = size(states); basis%p = p
        basis%bytes = bytes; basis%valid = .true.; basis%sector = .true.
        ok = .true.
        status = 'exact CAR sector basis; payload budget met'
    end subroutine

    !> Lexicographic index of the ordered pair i<j among m orbitals.
    pure integer function fock_pair_index(m, i, j) result(index)
        integer, intent(in) :: m, i, j
        index = (i - 1)*m - ((i - 1)*i)/2 + (j - i)
    end function

    !> Sparse lists of nonzero one-body (per annihilated orbital) and pair
    !> (per annihilated pair) coefficients: candidate or interval packet not
    !> an exact zero. Skipping exact-zero terms leaves every sum unchanged.
    subroutine index_couplings(model, max_bytes, ok, status)
        type(fock_hamiltonian_t), intent(inout) :: model
        integer(i64), intent(in) :: max_bytes
        logical, intent(out) :: ok
        character(*), intent(out) :: status
        integer(i64) :: bytes
        integer :: i, j, nh, nv, allocation

        ok = .false.
        nh = 0; nv = 0
        do j = 1, model%m
            do i = 1, model%m
                if (active_coefficient(model%h(i, j), model%hb(i, j))) nh = nh + 1
            end do
        end do
        do j = 1, model%p
            do i = 1, model%p
                if (active_coefficient(model%v(i, j), model%vb(i, j))) nv = nv + 1
            end do
        end do
        bytes = (int(model%m, i64)**2 + int(model%p, i64)**2)*24_i64 + &
                4_i64*(int(model%m + model%p + 2, i64) + nh + nv)
        if (bytes > max_bytes) then
            status = 'coupling index exceeds byte budget'; return
        end if
        if (allocated(model%hptr)) deallocate (model%hptr, model%hcol, &
                                               model%vptr, model%vcol)
        allocate (model%hptr(0:model%m), model%hcol(nh), model%vptr(0:model%p), &
                  model%vcol(nv), stat=allocation)
        if (allocation /= 0) then
            status = 'coupling index allocation failed'; return
        end if
        model%hptr(0) = 0; nh = 0
        do j = 1, model%m
            do i = 1, model%m
                if (.not. active_coefficient(model%h(i, j), model%hb(i, j))) cycle
                nh = nh + 1; model%hcol(nh) = i
            end do
            model%hptr(j) = nh
        end do
        model%vptr(0) = 0; nv = 0
        do j = 1, model%p
            do i = 1, model%p
                if (.not. active_coefficient(model%v(i, j), model%vb(i, j))) cycle
                nv = nv + 1; model%vcol(nv) = i
            end do
            model%vptr(j) = nv
        end do
        ok = .true.
    end subroutine

    pure logical function active_coefficient(coef, bound)
        real(dp), intent(in) :: coef
        type(interval_t), intent(in) :: bound
        active_coefficient = abs(coef) > 0 .or. abs(bound%lo) > 0 .or. &
                             abs(bound%hi) > 0
    end function

    !> Add a real symmetric one-body matrix dh (exact binary64 data, e.g. an
    !> external potential) to h and, outward, to its interval packets.
    subroutine add_fock_one_body(model, dh, max_bytes, ok, status)
        type(fock_hamiltonian_t), intent(inout) :: model
        real(dp), intent(in) :: dh(:, :)
        integer(i64), intent(in) :: max_bytes
        logical, intent(out) :: ok
        character(*), intent(out) :: status
        integer :: i, j
        real(dp) :: next
        type(interval_t) :: bound

        ok = .false.; status = 'invalid one-body perturbation'
        if (.not. model%valid) return
        if (any(shape(dh) /= [model%m, model%m])) return
        if (.not. all(ieee_is_finite(dh))) return
        if (any(dh > transpose(dh)) .or. any(dh < transpose(dh))) return
        do j = 1, model%m
            do i = 1, model%m
                call qf_add(model%h(i, j), dh(i, j), next); model%h(i, j) = next
                call qf_add_interval(model%hb(i, j), ipoint(dh(i, j)), bound)
                model%hb(i, j) = bound
            end do
        end do
        model%bounded = model%bounded .and. valid_matrix(model%hb)
        call index_couplings(model, max_bytes, ok, status)
        if (ok) status = 'one-body perturbation added'
    end subroutine

    logical function coupled_compatible(model, from, to)
        type(fock_hamiltonian_t), intent(in) :: model
        type(fock_basis_t), intent(in) :: from, to
        coupled_compatible = .false.
        if (.not. model%valid .or. .not. from%valid .or. .not. to%valid) return
        if (.not. allocated(model%hptr)) return
        coupled_compatible = model%m == from%m .and. model%m == to%m .and. &
                             model%n == from%n .and. model%n == to%n
    end function

    !> Next nonzero-coefficient transition out of determinant `ket` of `from`.
    !> cursor = [stage, outer, inner]; returns bra = 0 for a Pauli-allowed
    !> target outside `to` (leakage), which callers treat as failure.
    subroutine next_coupled(model, from, to, ket, cursor, bra, phase, sector, &
                            ci, cj, more)
        type(fock_hamiltonian_t), intent(in) :: model
        type(fock_basis_t), intent(in) :: from, to
        integer, intent(in) :: ket
        integer, intent(inout) :: cursor(3)
        integer, intent(out) :: bra, phase, sector, ci, cj
        logical, intent(out) :: more
        integer(i64) :: mask, next, step1, step2, step3
        integer :: j, k, q, p, local
        logical :: active

        more = .false.; bra = 0; phase = 0; sector = 0; ci = 0; cj = 0
        mask = from%states(ket)
        if (cursor(1) == 1) then
            do while (cursor(2) <= model%m)
                j = cursor(2)
                if (btest(mask, j - 1)) then
                    do while (cursor(3) <= model%hptr(j))
                        k = model%hcol(cursor(3)); cursor(3) = cursor(3) + 1
                        call fermion_apply(mask, j, .false., step1, phase, active)
                        call fermion_apply(step1, k, .true., next, local, active)
                        if (.not. active) cycle
                        phase = phase*local; sector = 1; ci = k; cj = j
                        bra = basis_index(to, next); more = .true.; return
                    end do
                end if
                cursor(2) = cursor(2) + 1
                if (cursor(2) <= model%m) cursor(3) = model%hptr(cursor(2) - 1) + 1
            end do
            cursor = [2, 1, 0]
        end if
        do while (cursor(2) <= model%p)
            q = cursor(2)
            if (btest(mask, from%pi(q) - 1) .and. btest(mask, from%pj(q) - 1)) then
                if (cursor(3) == 0) cursor(3) = model%vptr(q - 1) + 1
                do while (cursor(3) <= model%vptr(q))
                    p = model%vcol(cursor(3)); cursor(3) = cursor(3) + 1
                    call fermion_apply(mask, from%pi(q), .false., step1, phase, active)
                    call fermion_apply(step1, from%pj(q), .false., step2, local, active)
                    phase = phase*local
                    call fermion_apply(step2, from%pj(p), .true., step3, local, active)
                    if (.not. active) cycle
                    phase = phase*local
                    call fermion_apply(step3, from%pi(p), .true., next, local, active)
                    if (.not. active) cycle
                    phase = phase*local; sector = 2; ci = p; cj = q
                    bra = basis_index(to, next); more = .true.; return
                end do
            end if
            cursor(2) = cursor(2) + 1; cursor(3) = 0
        end do
    end subroutine

    integer function basis_index(basis, mask) result(index)
        type(fock_basis_t), intent(in) :: basis
        integer(i64), intent(in) :: mask
        integer :: lo, hi, mid
        index = 0; lo = 1; hi = basis%d
        do while (lo <= hi)
            mid = lo + (hi - lo)/2
            if (basis%states(mid) == mask) then
                index = mid; return
            else if (basis%states(mid) < mask) then
                lo = mid + 1
            else
                hi = mid - 1
            end if
        end do
    end function

    !> y = H x from basis `from` to basis `to` using only nonzero couplings.
    !> Fails closed (ok=.false.) if a nonzero term leaves `to`.
    subroutine fock_matvec_between(model, from, to, re, im, yr, yi, ok)
        type(fock_hamiltonian_t), intent(in) :: model
        type(fock_basis_t), intent(in) :: from, to
        real(dp), intent(in) :: re(:), im(:)
        real(dp), intent(out) :: yr(:), yi(:)
        logical, intent(out) :: ok
        integer :: cursor(3), ket, bra, phase, sector, i, j
        real(dp) :: coef, next
        type(interval_t) :: bound
        logical :: more

        ok = .false.
        if (size(yr) == to%d .and. size(yi) == to%d) then
            yr = 0; yi = 0
        end if
        if (.not. coupled_compatible(model, from, to)) return
        if (size(re) /= from%d .or. size(im) /= from%d) return
        if (size(yr) /= to%d .or. size(yi) /= to%d) return
        if (.not. all(ieee_is_finite(re)) .or. .not. all(ieee_is_finite(im))) return
        do ket = 1, from%d
            cursor = [1, 1, model%hptr(0) + 1]
            do
                call next_coupled(model, from, to, ket, cursor, bra, phase, &
                                  sector, i, j, more)
                if (.not. more) exit
                if (bra == 0) return
                call transition_coefficient(model, sector, i, j, coef, bound)
                call qf_accumulate(yr(bra), coef, re(ket), real(phase, dp), next)
                yr(bra) = next
                call qf_accumulate(yi(bra), coef, im(ket), real(phase, dp), next)
                yi(bra) = next
            end do
        end do
        ok = all(ieee_is_finite(yr)) .and. all(ieee_is_finite(yi))
    end subroutine

    !> Strict interval twin of fock_matvec_between.
    subroutine fock_matvec_between_interval(model, from, to, re, im, yr, yi, ok)
        type(fock_hamiltonian_t), intent(in) :: model
        type(fock_basis_t), intent(in) :: from, to
        type(interval_t), intent(in) :: re(:), im(:)
        type(interval_t), intent(out) :: yr(:), yi(:)
        logical, intent(out) :: ok
        integer :: cursor(3), ket, bra, phase, sector, i, j
        real(dp) :: coef
        type(interval_t) :: bound, next
        logical :: more

        ok = .false.
        if (size(yr) == to%d .and. size(yi) == to%d) then
            yr = ipoint(0._dp); yi = ipoint(0._dp)
        end if
        if (.not. coupled_compatible(model, from, to)) return
        if (.not. model%bounded) return
        if (size(re) /= from%d .or. size(im) /= from%d) return
        if (size(yr) /= to%d .or. size(yi) /= to%d) return
        if (.not. valid_intervals(re) .or. .not. valid_intervals(im)) return
        do ket = 1, from%d
            cursor = [1, 1, model%hptr(0) + 1]
            do
                call next_coupled(model, from, to, ket, cursor, bra, phase, &
                                  sector, i, j, more)
                if (.not. more) exit
                if (bra == 0) return
                call transition_coefficient(model, sector, i, j, coef, bound)
                call qf_accumulate_interval(yr(bra), bound, re(ket), &
                                            ipoint(real(phase, dp)), next)
                yr(bra) = next
                call qf_accumulate_interval(yi(bra), bound, im(ket), &
                                            ipoint(real(phase, dp)), next)
                yi(bra) = next
            end do
        end do
        ok = valid_intervals(yr) .and. valid_intervals(yi)
    end subroutine

    !> Residual amplitudes <r| a_{s_k} ... a_{s_1} |psi> for every determinant
    !> and every occupied ordered k-subset s (a_{s_1} applied first). Labels:
    !> orbital (k=1), pair index (k=2), triple key (k=3). Exact-zero
    !> amplitudes are skipped.
    subroutine residual_amplitudes(basis, re, im, order, residual, label, ar, ai, &
                                   count)
        type(fock_basis_t), intent(in) :: basis
        real(dp), intent(in) :: re(:), im(:)
        integer, intent(in) :: order
        integer(i64), allocatable, intent(out) :: residual(:)
        integer(i64), allocatable, intent(out) :: label(:)
        real(dp), allocatable, intent(out) :: ar(:), ai(:)
        integer, intent(out) :: count
        integer :: occupied(62), ket, a, b, c, n, phase, local, local3, total, cap
        integer(i64) :: mask, next, mask_b, mask_c, m64
        logical :: active

        n = basis%n; m64 = basis%m
        select case (order)
        case (1)
            cap = n
        case (2)
            cap = n*(n - 1)/2
        case default
            cap = n*(n - 1)*(n - 2)/6
        end select
        total = 0
        do ket = 1, basis%d
            if (abs(re(ket)) > 0 .or. abs(im(ket)) > 0) total = total + cap
        end do
        allocate (residual(total), label(total), ar(total), ai(total))
        count = 0
        do ket = 1, basis%d
            if (.not. (abs(re(ket)) > 0 .or. abs(im(ket)) > 0)) cycle
            mask = basis%states(ket)
            c = 0
            do a = 1, basis%m
                if (btest(mask, a - 1)) then
                    c = c + 1; occupied(c) = a
                end if
            end do
            do a = 1, n
                call fermion_apply(mask, occupied(a), .false., next, phase, active)
                if (order == 1) then
                    call push(next, int(occupied(a), i64), phase)
                    cycle
                end if
                do b = a + 1, n
                    call fermion_apply(next, occupied(b), .false., mask_b, local, active)
                    if (order == 2) then
                        call push(mask_b, int(fock_pair_index(basis%m, occupied(a), &
                                  occupied(b)), i64), phase*local)
                        cycle
                    end if
                    do c = b + 1, n
                        call fermion_apply(mask_b, occupied(c), .false., mask_c, &
                                           local3, active)
                        call push(mask_c, (int(occupied(a) - 1, i64)*m64 + &
                                  (occupied(b) - 1))*m64 + (occupied(c) - 1), &
                                  phase*local*local3)
                    end do
                end do
            end do
        end do
    contains
        subroutine push(r, l, sign)
            integer(i64), intent(in) :: r, l
            integer, intent(in) :: sign
            count = count + 1
            residual(count) = r; label(count) = l
            ar(count) = real(sign, dp)*re(ket); ai(count) = real(sign, dp)*im(ket)
        end subroutine
    end subroutine

    !> Factorized k-body transition density: for every residual r, collect
    !> ket amplitudes A_s(r) and bra amplitudes B_u(r); each product
    !> conj(B_u(r)) A_s(r) is passed to `emit(s, u, pr, pi)`.
    subroutine factorized_products(basis, bra_re, bra_im, ket_re, ket_im, order, &
                                   emit)
        type(fock_basis_t), intent(in) :: basis
        real(dp), intent(in) :: bra_re(:), bra_im(:), ket_re(:), ket_im(:)
        integer, intent(in) :: order
        interface
            subroutine emit(s, u, pr, pi)
                import :: i64, dp
                integer(i64), intent(in) :: s, u
                real(dp), intent(in) :: pr, pi
            end subroutine
        end interface
        integer(i64), allocatable :: rk(:), lk(:), rb(:), lb(:), keys(:)
        real(dp), allocatable :: akr(:), aki(:), abr(:), abi(:)
        integer, allocatable :: perm(:)
        integer :: nk, nb, first, last, x, y, ix, iy
        real(dp) :: pr, pi

        call residual_amplitudes(basis, ket_re, ket_im, order, rk, lk, akr, aki, nk)
        call residual_amplitudes(basis, bra_re, bra_im, order, rb, lb, abr, abi, nb)
        allocate (keys(nk + nb), perm(nk + nb))
        keys(1:nk) = rk(1:nk); keys(nk + 1:nk + nb) = rb(1:nb)
        call sort_key_order(keys, perm)
        first = 1
        do while (first <= nk + nb)
            last = first
            do while (last < nk + nb)
                if (keys(perm(last + 1)) /= keys(perm(first))) exit
                last = last + 1
            end do
            do x = first, last
                ix = perm(x)
                if (ix > nk) cycle
                do y = first, last
                    iy = perm(y) - nk
                    if (iy < 1) cycle
                    pr = abr(iy)*akr(ix) + abi(iy)*aki(ix)
                    pi = abr(iy)*aki(ix) - abi(iy)*akr(ix)
                    call emit(lk(ix), lb(iy), pr, pi)
                end do
            end do
            first = last + 1
        end do
    end subroutine

    !> Transition RDMs gamma_ij=<bra|a_j^dagger a_i|ket> and
    !> Gamma_ij,kl=<bra|a_k^dagger a_l^dagger a_j a_i|ket> (ordered-pair
    !> storage) by residual factorization: cost sum_r s_r^2 over the
    !> (N-1)/(N-2)-particle residuals instead of a scan of every pair word.
    subroutine fock_transition_rdms(basis, bra_re, bra_im, ket_re, ket_im, &
                                    gr, gi, g2r, g2i, ok)
        type(fock_basis_t), intent(in) :: basis
        real(dp), intent(in) :: bra_re(:), bra_im(:), ket_re(:), ket_im(:)
        real(dp), intent(out) :: gr(:, :), gi(:, :), g2r(:, :), g2i(:, :)
        logical, intent(out) :: ok

        ok = .false.
        gr = 0; gi = 0; g2r = 0; g2i = 0
        if (.not. basis%valid) return
        if (any([size(bra_re), size(bra_im), size(ket_re), size(ket_im)] /= basis%d)) return
        if (any(shape(gr) /= [basis%m, basis%m]) .or. &
            any(shape(gi) /= [basis%m, basis%m])) return
        if (any(shape(g2r) /= [basis%p, basis%p]) .or. &
            any(shape(g2i) /= [basis%p, basis%p])) return
        if (.not. all(ieee_is_finite(bra_re)) .or. .not. all(ieee_is_finite(bra_im))) return
        if (.not. all(ieee_is_finite(ket_re)) .or. .not. all(ieee_is_finite(ket_im))) return
        if (basis%n >= 1) call factorized_products(basis, bra_re, bra_im, ket_re, &
                                                   ket_im, 1, add_one)
        if (basis%n >= 2) call factorized_products(basis, bra_re, bra_im, ket_re, &
                                                   ket_im, 2, add_two)
        ok = all(ieee_is_finite(gr)) .and. all(ieee_is_finite(gi)) .and. &
             all(ieee_is_finite(g2r)) .and. all(ieee_is_finite(g2i))
    contains
        subroutine add_one(s, u, pr, pi)
            integer(i64), intent(in) :: s, u
            real(dp), intent(in) :: pr, pi
            real(dp) :: next
            call qf_add(gr(s, u), pr, next); gr(s, u) = next
            call qf_add(gi(s, u), pi, next); gi(s, u) = next
        end subroutine
        subroutine add_two(s, u, pr, pi)
            integer(i64), intent(in) :: s, u
            real(dp), intent(in) :: pr, pi
            real(dp) :: next
            call qf_add(g2r(s, u), pr, next); g2r(s, u) = next
            call qf_add(g2i(s, u), pi, next); g2i(s, u) = next
        end subroutine
    end subroutine

    !> Three-body RDM Gamma3_(t),(u) = <psi|a_u1^dagger a_u2^dagger a_u3^dagger
    !> a_t3 a_t2 a_t1|psi> for t1<t2<t3, u1<u2<u3, as reduced sparse
    !> coordinates with keys key(t)=((t1-1)m+t2-1)m+t3-1 (zero entries absent).
    subroutine fock_rdm3_entries(basis, re, im, tkey, ukey, vr, vi, count, ok)
        type(fock_basis_t), intent(in) :: basis
        real(dp), intent(in) :: re(:), im(:)
        integer(i64), allocatable, intent(out) :: tkey(:), ukey(:)
        real(dp), allocatable, intent(out) :: vr(:), vi(:)
        integer, intent(out) :: count
        logical, intent(out) :: ok
        call sparse_rdm_entries(basis, re, im, 3, int(basis%m, i64)**3, tkey, ukey, &
                                vr, vi, count, ok)
    end subroutine

    !> Two-body RDM Gamma_(s),(u) = <psi|a_u1^dagger a_u2^dagger a_s2 a_s1|psi>
    !> (pair indices s, u of fock_pair_indices; the fock_rdms convention) as
    !> reduced sparse coordinates (zero entries absent). Avoids the dense
    !> p x p matrix per state when only contractions are needed.
    subroutine fock_rdm2_entries(basis, re, im, skey, ukey, vr, vi, count, ok)
        type(fock_basis_t), intent(in) :: basis
        real(dp), intent(in) :: re(:), im(:)
        integer(i64), allocatable, intent(out) :: skey(:), ukey(:)
        real(dp), allocatable, intent(out) :: vr(:), vi(:)
        integer, intent(out) :: count
        logical, intent(out) :: ok
        call sparse_rdm_entries(basis, re, im, 2, int(basis%p, i64) + 1, skey, ukey, &
                                vr, vi, count, ok)
    end subroutine

    !> Reduced sparse k-body RDM entries (k = 2, 3) with label base m3.
    subroutine sparse_rdm_entries(basis, re, im, order, m3, tkey, ukey, vr, vi, &
                                  count, ok)
        type(fock_basis_t), intent(in) :: basis
        real(dp), intent(in) :: re(:), im(:)
        integer, intent(in) :: order
        integer(i64), intent(in) :: m3
        integer(i64), allocatable, intent(out) :: tkey(:), ukey(:)
        real(dp), allocatable, intent(out) :: vr(:), vi(:)
        integer, intent(out) :: count
        logical, intent(out) :: ok
        integer(i64), allocatable :: keys(:)
        real(dp), allocatable :: pr_(:), pi_(:)
        integer, allocatable :: perm(:)
        integer :: raw, x, y

        ok = .false.; count = 0
        allocate (tkey(0), ukey(0), vr(0), vi(0))
        if (.not. basis%valid .or. basis%n < order) return
        if (size(re) /= basis%d .or. size(im) /= basis%d) return
        if (.not. all(ieee_is_finite(re)) .or. .not. all(ieee_is_finite(im))) return
        allocate (keys(1024), pr_(1024), pi_(1024))
        raw = 0
        call factorized_products(basis, re, im, re, im, order, store)
        allocate (perm(raw))
        call sort_key_order(keys(1:raw), perm)
        deallocate (tkey, ukey, vr, vi)
        allocate (tkey(raw), ukey(raw), vr(raw), vi(raw))
        x = 1
        do while (x <= raw)
            y = x
            count = count + 1
            tkey(count) = keys(perm(x))/m3; ukey(count) = mod(keys(perm(x)), m3)
            vr(count) = 0; vi(count) = 0
            do while (y <= raw)
                if (keys(perm(y)) /= keys(perm(x))) exit
                vr(count) = vr(count) + pr_(perm(y))
                vi(count) = vi(count) + pi_(perm(y))
                y = y + 1
            end do
            x = y
        end do
        ok = all(ieee_is_finite(vr(1:count))) .and. all(ieee_is_finite(vi(1:count)))
    contains
        subroutine store(s, u, pr, pi)
            integer(i64), intent(in) :: s, u
            real(dp), intent(in) :: pr, pi
            integer(i64), allocatable :: kt(:)
            real(dp), allocatable :: rt(:)
            if (raw == size(keys)) then
                allocate (kt(2*raw)); kt(1:raw) = keys; call move_alloc(kt, keys)
                allocate (rt(2*raw)); rt(1:raw) = pr_; call move_alloc(rt, pr_)
                allocate (rt(2*raw)); rt(1:raw) = pi_; call move_alloc(rt, pi_)
            end if
            raw = raw + 1
            keys(raw) = s*m3 + u; pr_(raw) = pr; pi_(raw) = pi
        end subroutine
    end subroutine
end module kin6d_quantum_fock
