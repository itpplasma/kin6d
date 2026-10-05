program test_quantum_fock
    use, intrinsic :: iso_fortran_env, only: dp => real64, &
                                                             qp => real128, i64 => int64
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_is_finite
    use fortnum_interval, only: interval_t, ipoint
    use kin6d_quantum_fock, only: fock_basis_t, fock_hamiltonian_t, init_fock_basis, &
        fock_basis_info, get_fock_states, fock_pair_indices, fermion_apply, &
        prepare_fock_hamiltonian, prepare_jellium_hamiltonian, &
        assemble_fock_hamiltonian, assemble_fock_hamiltonian_interval, fock_matvec, &
        fock_matvec_interval, fock_matvec_jvp, fock_matvec_vjp, fock_density_rhs, &
        fock_operator_expectation, fock_rdms, fock_rdms_interval, &
        get_fock_integrals, fock_model_info, fock_pair_action, &
        fock_pair_action_interval, assemble_fock_hamiltonian_block, &
        assemble_fock_hamiltonian_block_interval, assemble_jellium_pair_sector
    implicit none
    integer, parameter :: m = 6, n = 3, d = 20, p = 15
    integer(i64), parameter :: budget = 100000000_i64
    type(fock_basis_t) :: basis, bad_basis
    type(fock_hamiltonian_t) :: model, plus, minus, jellium
    integer(i64) :: masks(d), bytes
    integer :: pi(p), pj(p), i, j, k, dim, mm, nn, pp
    integer :: momentum(3, m), spin(m)
    real(dp) :: h(m, m), v(p, p), hr(d, d), hi(d, d), re(d), &
                im(d), yr(d), yi(d)
    real(dp) :: dh(m, m), dv(p, p), dre(d), dimag(d), jyr(d), &
                jyi(d), byr(d), byi(d)
    real(dp) :: gh(m, m), gv(p, p), gre(d), gim(d), pr(d), &
                qi(d), mr(d), mi(d)
    real(dp) :: gr(m, m), gi(m, m), g2r(p, p), g2i(p, p), &
                rho_r(d, d), rho_i(d, d)
    real(dp) :: rr(d, d), ri(d, d), value_r, value_i, norm, eps, left, right
    type(interval_t) :: hb(d, d), hib(d, d), ybr(d), ybi(d)
    type(interval_t) :: grb(m, m), gib(m, m), g2rb(p, p), g2ib(p, p)
    complex(qp) :: ref_h(d, d), psi(d), y(d), ref_gamma(m, m), &
                   ref_two(m, m, m, m)
    complex(qp) :: rho(d, d), rho_dot(d, d), perturbed(d)
    complex(dp) :: orbitals(m, n), slater(d), wick
    integer :: occupied(m), count
    logical :: ok, bounded
    character(160) :: status
    call car_oracle()
    call pair_sector_constructor_oracle()
    call init_fock_basis(basis, m, n, 1000_i64, budget, ok, status)
    call check(ok, 'arbitrary fixed-N basis')
    call fock_basis_info(basis, mm, nn, dim, pp, bytes, ok)
    call check(ok .and. mm == m .and. nn == n .and. dim == d .and. pp == p, &
               'combinatorial dimensions')
    call get_fock_states(basis, masks, &
                         ok); call fock_pair_indices(basis, pi, pj, ok)
    call check(all(masks(2:) > masks(:d - 1)), &
               'ascending determinant order')
    do i = 1, m
        do j = i, m
            h(i, j) = sin(real(i + j, dp))*.2_dp; h(j, i) = h(i, j)
            dh(i, j) = cos(real(i + j, dp))*.1_dp; dh(j, i) = dh(i, j)
        end do
    end do
    do i = 1, p
        do j = i, p
            v(i, j) = cos(real(i + 2*j, dp))*.07_dp; v(j, i) = v(i, j)
            dv(i, j) = sin(real(2*i + j, dp))*.03_dp; dv(j, i) = dv(i, j)
        end do
    end do
    call prepare_fock_hamiltonian(model, basis, h, v, budget, ok, status)
    call check(ok, 'general signed interacting Hermitian input')
    call fock_model_info(model, bytes, bounded, ok)
    call check(ok .and. bounded, 'exact binary64 coefficient packet')
    call assemble_fock_hamiltonian(model, basis, hr, hi, ok)
    call check(ok .and. maxval(abs(hr - transpose(hr))) < 1.e-15_dp, &
               'Hermitian assembled candidate')
    call assemble_fock_hamiltonian_interval(model, basis, hb, hib, ok)
    call check(ok, 'strict assembled finite-basis matrix')
    call reference_hamiltonian(masks, pi, pj, h, v, ref_h)
    call check(maxval(abs(real(hr, qp) - real(ref_h, qp))) < 2.e-15_qp, &
               'independent full four-index 1/4 Hamiltonian oracle')
    call check(all(real(hb%lo, qp) <= real(ref_h, qp)) .and. &
               all(real(ref_h, qp) <= real(hb%hi, qp)), &
               'real128 matrix element containment')
    do i = 1, d
        re(i) = cos(real(i, dp)); im(i) = sin(real(2*i, dp))
        dre(i) = sin(real(i, &
                          dp))*.03_dp; dimag(i) = cos(real(3*i, dp))*.02_dp
        byr(i) = cos(real(4*i, dp)); byi(i) = sin(real(3*i, dp))
    end do
    norm = sqrt(sum(re**2 + im**2)); re = re/norm; im = im/norm
    psi = cmplx(real(re, qp), real(im, qp), qp); y = matmul(ref_h, psi)
    call fock_matvec(model, basis, re, im, yr, yi, ok)
    call check(ok .and. maxval(abs(cmplx(real(yr, qp), real(yi, qp), &
                                         qp) - y)) < 2.e-15_qp, &
               'matrix-free action independent oracle')
    call fock_matvec_interval(model, basis, ipoint(re), &
                              ipoint(im), ybr, ybi, ok)
    call check(ok .and. all(real(ybr%lo, qp) <= real(y, qp)) .and. &
               all(real(y, qp) <= real(ybr%hi, qp)) .and. all(real(ybi%lo, &
                                                                qp) <= aimag(y)) .and. &
               all(aimag(y) <= real(ybi%hi, qp)), &
               'strict matrix-free action containment')
    call fock_matvec_jvp(model, basis, re, im, dh, dv, dre, &
                         dimag, jyr, jyi, ok)
    call check(ok, 'composed matrix-free JVP')
    eps = 1.e-5_dp
    call prepare_fock_hamiltonian(plus, basis, h + eps*dh, &
                                  v + eps*dv, budget, ok, status)
    call prepare_fock_hamiltonian(minus, basis, h - eps*dh, &
                                  v - eps*dv, budget, ok, status)
    call fock_matvec(plus, basis, re + eps*dre, im + eps*dimag, pr, qi, ok)
    call fock_matvec(minus, basis, re - eps*dre, im - eps*dimag, mr, mi, ok)
    call check(maxval(abs((pr - mr)/(2*eps) - jyr)) < 2.e-10_dp .and. &
               maxval(abs((qi - mi)/(2*eps) - jyi)) < 2.e-10_dp, &
               'independent JVP finite difference')
    call fock_matvec_vjp(model, basis, re, im, byr, byi, gh, gv, &
                         gre, gim, ok)
    left = sum(byr*jyr + byi*jyi)
    right = sum(gh*dh) + sum(gv*dv) + sum(gre*dre + gim*dimag)
    call check(ok .and. abs(left - right) < 3.e-14_dp, &
               'symmetric-Frobenius/complex-real VJP pairing')
    call fock_rdms(basis, re, im, gr, gi, g2r, g2i, ok)
    call check(ok, 'arbitrary correlated full-state RDM query')
    call reference_rdms(masks, psi, ref_gamma, ref_two)
    call check(maxval(abs(cmplx(real(gr, qp), real(gi, qp), &
                                qp) - ref_gamma)) < 2.e-15_qp, &
               'one-body RDM exact ordering')
    call check(abs(sum([(gr(i, i), i=1, &
                         m)]) - n*sum(re**2 + im**2)) < 2.e-14_dp, &
               'Tr gamma=N norm2')
    call check(abs(2*sum([(g2r(i, i), i=1, &
                           p)]) - n*(n - 1)*sum(re**2 + im**2)) < 4.e-14_dp, &
               'ordered-pair RDM trace normalization')
    do i = 1, m
        do j = 1, m
            call check(abs(sum([(ref_two(i, k, j, k), k=1, &
                                 m)]) - (n - 1)*ref_gamma(i, &
                                                          j)) < 1.e-30_qp, &
                       'fixed-N contraction identity')
        end do
    end do
    call fock_rdms_interval(basis, re, im, 0._dp, grb, gib, g2rb, g2ib, ok)
    call check(ok, 'candidate-only strict RDM arithmetic')
    call check(all(real(grb%lo, qp) <= real(ref_gamma, qp)) .and. &
               all(real(ref_gamma, qp) <= real(grb%hi, qp)), &
               'one-body interval containment')
    do i = 1, p
        do j = 1, p
            call check(real(g2rb(i, j)%lo, &
                            qp) <= real(ref_two(pi(i), pj(i), &
                                                pi(j), pj(j)), qp) .and. &
                       real(ref_two(pi(i), pj(i), pi(j), pj(j)), &
                            qp) <= real(g2rb(i, j)%hi, qp), &
                       'pair RDM interval containment')
        end do
    end do
    perturbed = psi + cmplx(real(dre, qp), real(dimag, qp), qp)*.05_qp
    call check(sqrt(sum(abs(perturbed - psi)**2)) < .03_qp, &
               'independent state error premise')
    call reference_rdms(masks, perturbed, ref_gamma, ref_two)
    call fock_rdms_interval(basis, re, im, .03_dp, grb, gib, g2rb, g2ib, ok)
    call check(ok .and. all(real(grb%lo, qp) <= real(ref_gamma, qp)) .and. &
               all(real(ref_gamma, qp) <= real(grb%hi, qp)) .and. &
               all(real(gib%lo, &
             qp) <= aimag(ref_gamma)) .and. all(aimag(ref_gamma) <= real(gib%hi, qp)), &
               'state-radius observable coupling containment')
    rho = spread(psi, 2, d)*spread(conjg(psi), 1, d)
    rho_dot = cmplx(0._qp, -1._qp, qp)*(matmul(ref_h, &
                                               rho) - matmul(rho, ref_h))
    rho_r = real(rho, dp); rho_i = real(aimag(rho), dp)
    call fock_density_rhs(model, basis, rho_r, rho_i, 1._dp, rr, ri, ok)
    call check(ok .and. maxval(abs(cmplx(real(rr, qp), real(ri, qp), &
                                         qp) - rho_dot)) < 3.e-15_qp, &
               'full density commutator preserves all correlations')
    call check(abs(sum([(rr(i, i), i=1, d)])) < 3.e-15_dp .and. &
               abs(sum([(ri(i, i), i=1, d)])) < 3.e-15_dp, &
               'commutator trace conservation')
    call fock_operator_expectation(basis, rho_r, rho_i, [2, 1], &
                                   [.true., .false.], &
                                   value_r, value_i, ok)
    call check(ok .and. abs(cmplx(value_r, value_i, &
                                  dp) - cmplx(gr(1, 2), gi(1, 2), &
                                              dp)) < 2.e-15_dp, &
               'rho versus psi observable consistency')
    orbitals = cmplx(0._dp, 0._dp, dp)
    do i = 1, n
        orbitals(i, i) = cos(.17_dp*i)*exp(cmplx(0._dp, .2_dp*i, dp))
        orbitals(i + n, i) = sin(.17_dp*i)*exp(cmplx(0._dp, &
                                                     .2_dp*(i + n), dp))
    end do
    do k = 1, d
        call occupation_list(masks(k), occupied, count)
        slater(k) = det3(orbitals(occupied(:n), :))
    end do
    call fock_rdms(basis, real(slater, dp), aimag(slater), gr, &
                   gi, g2r, g2i, ok)
    do i = 1, p
        do j = 1, p
            wick = cmplx(gr(pi(i), pi(j)), gi(pi(i), pi(j)), dp)* &
                   cmplx(gr(pj(i), pj(j)), gi(pj(i), pj(j)), dp) - &
                   cmplx(gr(pi(i), pj(j)), gi(pi(i), pj(j)), dp)* &
                   cmplx(gr(pj(i), pi(j)), gi(pj(i), pi(j)), dp)
            call check(abs(cmplx(g2r(i, j), g2i(i, j), &
                                 dp) - wick) < 3.e-15_dp, &
                       'nontrivial complex Slater Wick oracle')
        end do
    end do
    call jellium_oracle()
    call pair_action_oracle()
    call block_oracle()
    call noninteracting_sectors()
    call init_fock_basis(bad_basis, 62, 31, 1000_i64, budget, ok, status)
    call check(.not. ok, 'combinatorial budget rejected before enumeration')
    call init_fock_basis(bad_basis, 6, 3, 1000_i64, 1_i64, ok, status)
    call check(.not. ok, 'basis memory budget failure')
    h(1, 2) = h(1, 2) + .1_dp
    call prepare_fock_hamiltonian(plus, basis, h, v, budget, ok, status)
    call check(.not. ok, 'non-Hermitian input rejected')
    re(1) = ieee_value(0._dp, ieee_quiet_nan)
    call fock_matvec(model, basis, re, im, yr, yi, ok)
    call check(.not. ok, 'nonfinite state rejected')
    call fock_rdms_interval(basis, real(slater, dp), &
                            aimag(slater), -1._dp, grb, &
                            gib, g2rb, g2ib, ok)
    call check(.not. ok, 'negative state radius rejected')
    print '(a)', &
        'PASS: general fixed-N CAR, full four-index oracle, JVP/'// &
        'VJP, RDM/W'// &
        'ick, jellium and budgets'
contains
    subroutine pair_sector_constructor_oracle()
        integer :: modes(3, 3)
        real(dp) :: real_h(3, 3), imag_h(3, 3)
        type(interval_t) :: real_b(3, 3), imag_b(3, 3)
        real(qp) :: ell, wave, piq, reference(3, 3)
        integer :: i, j
        modes(:, 1) = [0, 0, 0]; modes(:, 2) = [1, -2, 3]
        modes(:, 3) = [-2, 1, -1]
        call assemble_jellium_pair_sector(modes, 5.2_dp, 0.9_dp, 1.3_dp, &
                             0.8_dp, budget, real_h, imag_h, real_b, imag_b, ok, status)
        call check(ok, 'direct arbitrary physical 3D pair sector')
        ell = real(5.2_dp, qp); piq = acos(-1._qp); wave = 2*piq/ell
        reference = 0
        do i = 1, 3
            reference(i, i) = real(0.9_dp, qp)**2*wave**2* &
                              sum(modes(:, i)**2)/real(1.3_dp, qp)
            do j = 1, i - 1
                reference(i, j) = 4*piq*real(0.8_dp, qp)/ &
                                  (ell**3*wave**2*sum((modes(:, i) - modes(:, j))**2))
                reference(j, i) = reference(i, j)
            end do
        end do
        call check(all(real(real_b%lo, qp) <= reference) .and. &
                all(real(real_b%hi, qp) >= reference), 'direct sector physical real128')
        call check(.not. any(real_h > transpose(real_h)) .and. &
                 .not. any(real_h < transpose(real_h)), 'direct sector exact Hermitian')
        call check(maxval(abs(imag_h)) <= 0._dp .and. &
                maxval(abs(imag_b%lo)) <= 0._dp .and. maxval(abs(imag_b%hi)) <= 0._dp, &
                   'real physical Hamiltonian')
        call assemble_jellium_pair_sector(modes, 5.2_dp, 0.9_dp, 1.3_dp, &
                            0.8_dp, 431_i64, real_h, imag_h, real_b, imag_b, ok, status)
        call check(.not. ok, 'direct sector output budget before work')
        call assemble_jellium_pair_sector(modes, 5.2_dp, 0._dp, 1.3_dp, &
                             0.8_dp, budget, real_h, imag_h, real_b, imag_b, ok, status)
        call check(.not. ok, 'direct sector positive hbar contract')
        modes(:, 3) = modes(:, 2)
        call assemble_jellium_pair_sector(modes, 5.2_dp, 0.9_dp, 1.3_dp, &
                             0.8_dp, budget, real_h, imag_h, real_b, imag_b, ok, status)
        call check(.not. ok, 'direct sector duplicate momentum rejected')
        modes(:, 3) = [-2, 1, -1]
        call assemble_jellium_pair_sector(modes, 5.2_dp, 0.9_dp, 1.3_dp, &
                           0.8_dp, budget, real_h, imag_h, real_b, imag_b, ok, status, &
                                          certify=.false.)
        call check(ok .and. maxval(abs(real(real_h, qp) - reference)) < 1.e-14_qp, &
                   'candidate-only constructor retains physical candidate')
        call check(.not. any(ieee_is_finite(real_b%lo)) .and. &
                   .not. any(ieee_is_finite(real_b%hi)), &
                   'candidate-only constructor strict packet explicitly unavailable')

    end subroutine
    subroutine check(condition, message)
        logical, intent(in) :: condition
        character(*), intent(in) :: message
        if (.not. condition) then
            print '(a)', 'FAIL: '//message; error stop 1
        end if
    end subroutine
    subroutine car_oracle()
        integer, parameter :: cm = 5, cd = 32
        integer :: a(cd, cd, cm), eye(cd, cd), i, j, k, phase, &
                   expected(cd, cd)
        integer(i64) :: target
        logical :: active
        a = 0; eye = 0
        do k = 1, cd
            eye(k, k) = 1
            do i = 1, cm
                call fermion_apply(int(k - 1, i64), i, .false., &
                                   target, phase, active)
                if (active) a(int(target) + 1, k, i) = phase
            end do
        end do
        do i = 1, cm
            do j = 1, cm
                expected = 0
                if (i == j) expected = eye
                call check(all(matmul(a(:, :, i), transpose(a(:, :, j))) + &
                               matmul(transpose(a(:, :, j)), &
                                      a(:, :, i)) == expected), &
                           'exact mixed CAR')
                call check(all(matmul(a(:, :, i), a(:, :, &
                                                    j)) + matmul(a(:, :, j), &
                                                                 a(:, :, i)) == 0), &
                           'exact annihilator CAR')
                call check(all(matmul(transpose(a(:, :, i)), &
                                      transpose(a(:, :, j))) + &
                               matmul(transpose(a(:, :, j)), &
                                      transpose(a(:, :, i))) == 0), &
                           'exact creator CAR')
            end do
        end do
    end subroutine
    subroutine occupation_list(mask, list, count)
        integer(i64), intent(in) :: mask
        integer, intent(out) :: list(:), count
        integer(i64) :: work
        integer :: i
        work = mask; count = 0; list = 0
        do i = 1, size(list)
            if (mod(work, 2_i64) == 1) then
                count = count + 1; list(count) = i
            end if
            work = work/2
        end do
    end subroutine
    subroutine reference_word(mask, word, create, target, phase, active)
        ! Independent ordered-list insertion/deletion oracle; no bit parity,
        ! production CAR calls, lookup routine, or numerical leaves.
        integer(i64), intent(in) :: mask
        integer, intent(in) :: word(:)
        logical, intent(in) :: create(:)
        integer(i64), intent(out) :: target
        integer, intent(out) :: phase
        logical, intent(out) :: active
        integer :: list(m), count, k, pos, i
        logical :: occupied
        call occupation_list(mask, list, count)
        phase = 1; active = .false.; target = 0
        do k = size(word), 1, -1
            pos = 1
            do while (pos <= count)
                if (list(pos) >= word(k)) exit
                pos = pos + 1
            end do
            occupied = .false.
            if (pos <= count) occupied = list(pos) == word(k)
            if (occupied .eqv. create(k)) then
                phase = 0; return
            end if
            phase = phase*(-1)**(pos - 1)
            if (create(k)) then
                do i = count, pos, -1
                    list(i + 1) = list(i)
                end do
                list(pos) = word(k); count = count + 1
            else
                do i = pos, count - 1
                    list(i) = list(i + 1)
                end do
                count = count - 1
            end if
        end do
        do i = 1, count
            target = target + 2_i64**(list(i) - 1)
        end do
        active = .true.
    end subroutine
    integer function reference_index(states, target) result(index)
        integer(i64), intent(in) :: states(:), target
        integer :: i
        index = 0
        do i = 1, size(states)
            if (states(i) == target) then
                index = i; return
            end if
        end do
    end function
    real(qp) function full_pair(values, pi, pj, i, j, k, l) result(value)
        real(dp), intent(in) :: values(:, :)
        integer, intent(in) :: pi(:), pj(:), i, j, k, l
        integer :: a, b, sign, index
        value = 0
        if (i == j .or. k == l) return
        sign = 1
        if (i > j) sign = -sign
        if (k > l) sign = -sign
        a = 0; b = 0
        do index = 1, size(pi)
            if (pi(index) == min(i, j) .and. pj(index) == max(i, &
                                                              j)) a = index
            if (pi(index) == min(k, l) .and. pj(index) == max(k, &
                                                              l)) b = index
        end do
        value = real(sign, qp)*real(values(a, b), qp)
    end function
    subroutine reference_hamiltonian(states, pi, pj, h, v, matrix)
        integer(i64), intent(in) :: states(:)
        integer, intent(in) :: pi(:), pj(:)
        real(dp), intent(in) :: h(:, :), v(:, :)
        complex(qp), intent(out) :: matrix(:, :)
        integer :: ket, bra, i, j, k, l, phase
        integer(i64) :: target
        logical :: active
        matrix = cmplx(0._qp, 0._qp, qp)
        do ket = 1, size(states)
            do i = 1, m
                do j = 1, m
                    call reference_word(states(ket), [i, j], &
                                        [.true., .false.], &
                                        target, phase, active)
                    if (active) then
                        bra = reference_index(states, target)
                        matrix(bra, ket) = matrix(bra, &
                                                  ket) + phase*real(h(i, j), qp)
                    end if
                    do k = 1, m
                        do l = 1, m
                            call reference_word(states(ket), [i, j, l, k], &
                                                [.true., .true., .false., &
                                                 .false.], target, phase, active)
                            if (.not. active) cycle
                            bra = reference_index(states, target)
                            matrix(bra, ket) = matrix(bra, ket) + &
                                               real(phase, &
                                                    qp)*full_pair(v, pi, &
                                                                  pj, i, j, k, l)/4
                        end do
                    end do
                end do
            end do
        end do
    end subroutine
    subroutine reference_rdms(states, psi, gamma, two)
        integer(i64), intent(in) :: states(:)
        complex(qp), intent(in) :: psi(:)
        complex(qp), intent(out) :: gamma(:, :), two(:, :, :, :)
        integer :: ket, bra, i, j, k, l, phase
        integer(i64) :: target
        logical :: active
        gamma = 0; two = 0
        do i = 1, m
            do j = 1, m
                do ket = 1, size(states)
                    call reference_word(states(ket), [j, i], &
                                        [.true., .false.], &
                                        target, phase, active)
                    if (.not. active) cycle
                    bra = reference_index(states, target)
                    gamma(i, j) = gamma(i, &
                                        j) + phase*conjg(psi(bra))*psi(ket)
                end do
                do k = 1, m
                    do l = 1, m
                        do ket = 1, size(states)
                            call reference_word(states(ket), [k, l, j, i], &
                                                [.true., .true., .false., &
                                                 .false.], target, phase, active)
                            if (.not. active) cycle
                            bra = reference_index(states, target)
                            two(i, j, k, l) = two(i, j, k, &
                                                  l) + phase*conjg(psi(bra))*psi(ket)
                        end do
                    end do
                end do
            end do
        end do
    end subroutine
    pure complex(dp) function det3(a) result(value)
        complex(dp), intent(in) :: a(3, 3)
        value = a(1, 1)*a(2, 2)*a(3, 3) + a(1, 2)*a(2, 3)*a(3, &
                                                        1) + a(1, 3)*a(2, 1)*a(3, 2) - &
                a(1, 3)*a(2, 2)*a(3, 1) - a(1, 2)*a(2, 1)*a(3, &
                                                            3) - a(1, 1)*a(2, 3)*a(3, 2)
    end function
    subroutine jellium_oracle()
        real(dp) :: one(m, m), pair(p, p)
        type(interval_t) :: ob(m, m), pb(p, p)
        real(qp) :: exact_one, exact_pair, piq, direct, exchange, qn2
        integer :: i, j, k, l, a, b
        momentum(:, 1) = [0, 0, 0]; momentum(:, 2) = [0, 0, 0]
        momentum(:, 3) = [1, 0, 0]; momentum(:, 4) = [-1, 0, 0]
        momentum(:, 5) = [0, 1, 0]; momentum(:, 6) = [0, 0, 1]
        spin = [1, -1, 1, -1, 1, -1]
        call prepare_jellium_hamiltonian(jellium, basis, &
                                         momentum, spin, 5.2_dp, &
                                         1.1_dp, &
                                         .75_dp, .67_dp, budget, ok, status)
        call check(ok, 'configurable three-dimensional jellium basis')
        call get_fock_integrals(jellium, one, pair, ok, ob, pb)
        call check(ok, 'physical coefficient interval packet')
        piq = acos(-1._qp)
        do i = 1, m
            exact_one = real(1.1_dp, qp)**2*(2*piq/real(5.2_dp, qp))**2* &
                        sum(real(momentum(:, i), &
                                 qp)**2)/(2*real(.75_dp, qp))
            call check(real(ob(i, i)%lo, &
                            qp) <= exact_one .and. exact_one <= real(ob(i, i)%hi, qp), &
                       'physical kinetic energy real128 oracle containment')
        end do
        do a = 1, p
            i = pi(a); j = pj(a)
            do b = 1, p
                k = pi(b); l = pj(b); direct = 0; exchange = 0
                if (all(momentum(:, i) + momentum(:, j) == momentum(:, &
                                                              k) + momentum(:, l))) then
                    if (spin(i) == spin(k) .and. spin(j) == spin(l)) then
                        qn2 = sum((real(momentum(:, i), &
                                        qp) - real(momentum(:, k), qp))**2)
                        if (qn2 > 0) direct = real(.67_dp, &
                                                   qp)/(piq*real(5.2_dp, qp)*qn2)
                    end if
                    if (spin(i) == spin(l) .and. spin(j) == spin(k)) then
                        qn2 = sum((real(momentum(:, i), &
                                        qp) - real(momentum(:, l), qp))**2)
                        if (qn2 > 0) exchange = real(.67_dp, &
                                                     qp)/(piq*real(5.2_dp, qp)*qn2)
                    end if
                end if
                exact_pair = direct - exchange
                call check(abs(real(pair(a, b), &
                                    qp) - exact_pair) < 3.e-16_qp, &
                           'momentum/spin selection and unscreen'// &
                           'ed Coulomb oracle')
                call check(real(pb(a, b)%lo, &
                          qp) <= exact_pair .and. exact_pair <= real(pb(a, b)%hi, qp), &
                           'physical pair matrix real128 containment')
            end do
        end do
        call check(abs(pair(1, 1)) < 1.e-30_dp, &
                   'zero-transfer/background mode removed')
        call check(maxval(abs(pair)) > .01_dp, &
                   'nonzero genuine interaction retained')
        spin(2) = 1
        call prepare_jellium_hamiltonian(plus, basis, momentum, &
                                         spin, 5.2_dp, &
                                         1._dp, 1._dp, &
                                         1._dp, budget, ok, status)
        call check(.not. ok, 'duplicate spin-orbital rejected')
        spin(2) = -1; momentum(1, 1) = -huge(1) - 1
        call prepare_jellium_hamiltonian(plus, basis, momentum, &
                                         spin, 5.2_dp, &
                                         1._dp, 1._dp, &
                                         1._dp, budget, ok, status)
        call check(.not. ok, &
                   'integer mode guard rejects extreme input wit'// &
                   'hout overflow')
    end subroutine
    subroutine pair_action_oracle()
        type(fock_basis_t) :: remaining
        integer(i64) :: residue(m), target
        real(dp) :: fr(m), fi(m), sr(m), si(m), ar(m), ai(m)
        type(interval_t) :: ab_r(m), ab_i(m)
        complex(qp) :: exact(m), fc(m), sc(m), state(d)
        integer :: i, j, ket, bra, phase
        logical :: active
        call init_fock_basis(remaining, m, n - 2, 1000_i64, &
                             budget, ok, status)
        call check(ok, 'residual N-2 basis')
        call get_fock_states(remaining, residue, ok)
        do i = 1, m
            fr(i) = cos(.2_dp*i); fi(i) = sin(.3_dp*i)
            sr(i) = sin(.17_dp*i); si(i) = cos(.13_dp*i)
        end do
        fc = cmplx(real(fr, qp), real(fi, qp), qp)
        sc = cmplx(real(sr, qp), real(si, qp), qp)
        state = cmplx(real(slater, qp), real(aimag(slater), qp), &
                      qp); exact = 0
        do ket = 1, d
            do i = 1, m
                do j = 1, m
                    call reference_word(masks(ket), [j, i], &
                                        [.false., .false.], &
                                        target, phase, active)
                    if (.not. active) cycle
                    bra = reference_index(residue, target)
                    exact(bra) = exact(bra) + real(phase, &
                                                   qp)*fc(i)*sc(j)*state(ket)
                end do
            end do
        end do
        call fock_pair_action(basis, remaining, fr, fi, sr, si, &
                              real(slater, dp), &
                              aimag(slater), ar, ai, ok)
        call check(ok .and. maxval(abs(cmplx(real(ar, qp), real(ai, qp), &
                                             qp) - exact)) < 3.e-15_qp, &
                   'generic weighted pair annihilation independe'// &
                   'nt ordered-list oracle')
        call fock_pair_action_interval(basis, remaining, ipoint(fr), &
                                       ipoint(fi), ipoint(sr), &
                                       ipoint(si), ipoint(real(slater, dp)), &
                                       ipoint(aimag(slater)), ab_r, ab_i, ok)
        call check(ok .and. all(real(ab_r%lo, qp) <= real(exact, &
                                                          qp)) .and. &
                   all(real(exact, qp) <= real(ab_r%hi, &
                                               qp)) .and. all(real(ab_i%lo, &
                                                            qp) <= aimag(exact)) .and. &
                   all(aimag(exact) <= real(ab_i%hi, qp)), &
                   'strict pair annihilation containment')
        call fock_pair_action(basis, remaining, fr, fi, fr, fi, &
                              real(slater, dp), &
                              aimag(slater), ar, ai, ok)
        call check(ok .and. maxval(abs(ar) + abs(ai)) < 4.e-15_dp, &
                   'Pauli equal-field pair exclusion')
    end subroutine
    subroutine block_oracle()
        type(fock_basis_t) :: b
        type(fock_hamiltonian_t) :: a
        integer(i64) :: states(15)
        integer :: wave(3, m), spins(m), selection(3), list(m), &
                   count, i, j, k, q
        integer :: plus_wave(3, 3), orientation(3)
        real(dp) :: block(3, 3), imaginary(3, 3), single(1, 1), &
                    single_i(1, 1)
        type(interval_t) :: block_b(3, 3), imaginary_b(3, 3)
        real(qp) :: expected, piq, transfer
        logical :: invariant
        wave(:, 1) = [0, 0, 0]; wave(:, 2) = [0, 0, 0]
        wave(:, 3) = [1, 0, 0]; wave(:, 4) = [-1, 0, 0]
        wave(:, 5) = [0, 1, 1]; wave(:, 6) = [0, -1, -1]
        spins = [1, -1, -1, 1, 1, -1]
        call init_fock_basis(b, m, 2, 1000_i64, budget, ok, status)
        call get_fock_states(b, states, ok)
        q = 0
        do k = 1, 15
            call occupation_list(states(k), list, count)
            if (spins(list(1)) + spins(list(2)) /= 0) cycle
            if (any(wave(:, list(1)) + wave(:, list(2)) /= 0)) cycle
            q = q + 1; selection(q) = k
            orientation(q) = 1; plus_wave(:, q) = wave(:, list(1))
            if (spins(list(1)) < 0) then
                orientation(q) = -1; plus_wave(:, q) = wave(:, list(2))
            end if
        end do
        call check(q == 3, 'exact integer momentum/spin sector selection')
        call prepare_jellium_hamiltonian(a, b, wave, spins, &
                                         5.2_dp, 1._dp, 1._dp, &
                                         .67_dp, budget, ok, status)
        call assemble_fock_hamiltonian_block(a, b, selection, &
                                             block, imaginary, &
                                             ok, invariant)
        call check(ok .and. invariant, &
                   'structural zero leakage proof for exact sector')
        call assemble_fock_hamiltonian_block_interval(a, b, &
                                                      selection, block_b, &
                                                      imaginary_b, ok, invariant)
        call check(ok .and. invariant, 'strict structural sector proof')
        piq = acos(-1._qp)
        do i = 1, 3
            do j = 1, 3
                expected = 0
                if (i == j) then
                    expected = (2*piq/real(5.2_dp, &
                                           qp))**2*sum(real(plus_wave(:, i), qp)**2)
                else
                    transfer = sum((real(plus_wave(:, i), &
                                         qp) - real(plus_wave(:, &
                                                              j), qp))**2)
                    expected = orientation(i)*orientation(j)* &
                               real(.67_dp, &
                                    qp)/(piq*real(5.2_dp, qp)*transfer)
                end if
                call check(abs(real(block(i, j), &
                                    qp) - expected) < 2.e-15_qp, &
                           'selected-block independent plane-wav'// &
                           'e physical oracle')
                call check(real(block_b(i, j)%lo, qp) <= expected .and. &
                           expected <= real(block_b(i, j)%hi, qp), &
                           'physical selected-block containment')
            end do
        end do
        call assemble_fock_hamiltonian_block(a, b, selection(:1), single, &
                                             single_i, ok, invariant)
        call check(ok .and. .not. invariant, &
                   'arbitrary projection cannot claim closed dynamics')
        selection(2) = selection(1)
        call assemble_fock_hamiltonian_block(a, b, selection, &
                                             block, imaginary, &
                                             ok, invariant)
        call check(.not. ok, 'duplicate block indices rejected')
    end subroutine
    subroutine noninteracting_sectors()
        type(fock_basis_t) :: b
        type(fock_hamiltonian_t) :: a
        integer :: particles, dim, pairs, mm, nn, i, j, occupied(m), count
        integer, parameter :: sectors(4) = [0, 1, 3, 6]
        integer(i64), allocatable :: states(:)
        real(dp), allocatable :: one(:, :), two(:, :), mr(:, :), mi(:, :)
        real(dp) :: expected
        do particles = 1, size(sectors)
            call init_fock_basis(b, m, sectors(particles), &
                                 1000_i64, budget, ok, status)
            call fock_basis_info(b, mm, nn, dim, pairs, bytes, ok)
            allocate (states(dim), one(m, m), two(pairs, pairs), &
                      mr(dim, dim), mi(dim, dim))
            call get_fock_states(b, states, ok)
            one = 0; two = 0
            do i = 1, m
                one(i, i) = real(i, dp)/2
            end do
            call prepare_fock_hamiltonian(a, b, one, two, &
                                          budget, ok, status)
            call assemble_fock_hamiltonian(a, b, mr, mi, ok)
            do i = 1, dim
                call occupation_list(states(i), occupied, count)
                expected = sum(real(occupied(:count), dp))/2
                call check(abs(mr(i, i) - expected) < 1.e-15_dp, &
                           'N-independent free energies')
                do j = 1, dim
                    if (i /= j) call check(abs(mr(i, j)) < 1.e-30_dp, &
                                           'free diagonal Hamiltonian')
                end do
            end do
            deallocate (states, one, two, mr, mi)
        end do
    end subroutine
end program test_quantum_fock
