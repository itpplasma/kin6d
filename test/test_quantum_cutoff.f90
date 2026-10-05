program test_quantum_cutoff
    use, intrinsic :: iso_fortran_env, only: dp => real64, qp => real128, i64 => int64
    use fortnum_interval, only: interval_t
    use kin6d_quantum_fock, only: fock_basis_t, fock_hamiltonian_t, &
                                  init_fock_basis, fock_basis_info, get_fock_states, &
                         prepare_jellium_hamiltonian, assemble_fock_hamiltonian_block, &
                  assemble_fock_hamiltonian_block_interval, assemble_jellium_pair_sector
    use kin6d_quantum_evolution, only: quantum_evolution_t, &
         initialize_quantum_evolution, solve_quantum_evolution, sample_quantum_evolution
    use kin6d_quantum_reference_gallery, only: quantum_observable_packet_t, &
                 quantum_initial_norm, quantum_matrix_defect, quantum_physical_radius, &
                                 prepare_quantum_pair_sector, sample_quantum_observables
    implicit none
    integer, parameter :: expected(3) = [7, 19, 27], times(4) = [0, 1, 4, 16]
    real(dp), parameter :: length = 6.2831853071795864769_dp, h = 1._dp/64
    real(dp) :: bins(16, 4, 3), uncertainty(4, 3), differences(4, 3)
    real(qp) :: truth_bins(16, 4, 3)
    integer :: cutoff, time_index, unit
    do cutoff = 1, 3
        call run_cutoff(cutoff, bins(:, :, cutoff), truth_bins(:, :, cutoff), &
                        uncertainty(:, cutoff))
    end do
    call require(maxval(abs(bins(:, 1, :) - 1._dp)) < 1.e-12_dp, &
                 'all cutoffs start from the same uniform physical pair state')
    differences(:, 1) = maxval(abs(bins(:, :, 1) - bins(:, :, 2)), dim=1)
    differences(:, 2) = maxval(abs(bins(:, :, 2) - bins(:, :, 3)), dim=1)
    differences(:, 3) = maxval(abs(bins(:, :, 1) - bins(:, :, 3)), dim=1)
    call require(maxval(differences(2:, :)) > 1.e-3_dp, &
                 'actual interaction quench produces resolved cutoff differences')
    open (newunit=unit, file='quantum-cutoff-evidence.csv', status='replace')
 write (unit, '(a)') '# finite-basis model differences; continuum remainder unavailable'
    write (unit, '(a)') 'time,max_bin_7v19,max_bin_19v27,max_bin_7v27,'// &
        'max_radius_7,max_radius_19,max_radius_27'
    do time_index = 1, 4
        write (unit, '(i0,6(",",es24.16e3))') times(time_index), &
            differences(time_index, :), uncertainty(time_index, :)
    end do
    close (unit)
    print *, 'PASS: physical quench 7/19/27, real128 evolution and 16-bin voxel oracle'
    do time_index = 1, 4
        print '(a,i2,a,3es12.4)', 'time=', times(time_index), ' cutoff differences=', &
            differences(time_index, :)
    end do
contains
    subroutine run_cutoff(cutoff, values, truths, radii)
        integer, intent(in) :: cutoff
        real(dp), intent(out) :: values(16, 4), radii(4)
        real(qp), intent(out) :: truths(16, 4)
        type(fock_basis_t) :: basis
        type(fock_hamiltonian_t) :: model
        type(quantum_evolution_t) :: state
        type(quantum_observable_packet_t) :: packet
        type(interval_t) :: norm0
        type(interval_t), allocatable :: hrb(:, :), hib(:, :)
        integer, allocatable :: modes(:, :), momentum(:, :), spin(:), indices(:)
        integer(i64), allocatable :: masks(:)
        real(dp), allocatable :: hr(:, :), hi(:, :), re(:), im(:), sign(:)
        real(dp), allocatable :: direct_h(:, :), direct_i(:, :)
        type(interval_t), allocatable :: direct_b(:, :), direct_ib(:, :)
        real(qp), allocatable :: physical_h(:, :)
        complex(qp), allocatable :: oracle(:)
        integer :: kx, ky, kz, shell, k, j, count, m, n, d, p, up, down, cell, b, t
        integer :: previous_time
        integer(i64) :: bytes, mask
        real(dp) :: radius, physical_radius, defect, x, result(6), lo(6), upper(6)
        real(qp) :: truth(6), oracle_error, difference
        logical :: ok, invariant, bounded
        character(160) :: status
        count = expected(cutoff)
        allocate (modes(3, count), momentum(3, 2*count), spin(2*count))
        k = 0
        do shell = 0, cutoff
            do kz = -1, 1
                do ky = -1, 1
                    do kx = -1, 1
                        if (kx*kx + ky*ky + kz*kz /= shell) cycle
                        k = k + 1; modes(:, k) = [kx, ky, kz]
                        momentum(:, 2*k - 1) = modes(:, k)
                        momentum(:, 2*k) = modes(:, k)
                        spin(2*k - 1:2*k) = [1, -1]
                    end do
                end do
            end do
        end do
        call require(k == count, 'genuine three-dimensional shell cardinality')
        call init_fock_basis(basis, 2*count, 2, 10000_i64, 100000000_i64, ok, status)
        call require(ok, 'parent basis construction')
        call fock_basis_info(basis, m, n, d, p, bytes, ok)
        allocate (masks(d), indices(count), sign(count))
        call get_fock_states(basis, masks, ok)
        call require(ok, 'parent determinant masks')
        do k = 1, count
            up = 2*k - 1; down = 0
            do j = 1, count
                if (all(modes(:, j) == -modes(:, k))) down = 2*j
            end do
            call require(down > 0, 'inversion symmetric three-dimensional modes')
            mask = ior(shiftl(1_i64, up - 1), shiftl(1_i64, down - 1))
            indices(k) = findloc(masks, mask, dim=1)
            call require(indices(k) > 0, 'selected exact parent determinant')
            sign(k) = 1._dp
            if (up > down) sign(k) = -1._dp
        end do
        call prepare_jellium_hamiltonian(model, basis, momentum, spin, length, &
                                         1._dp, 1._dp, 2._dp, &
                                         100000000_i64, ok, status)
        call require(ok, 'physical unscreened jellium Hamiltonian')
        allocate (hr(count, count), hi(count, count), hrb(count, count), &
                  hib(count, count), re(count), im(count), physical_h(count, count), &
                  oracle(count))
        call assemble_fock_hamiltonian_block(model, basis, indices, hr, hi, &
                                             ok, invariant)
        call require(ok .and. invariant, 'exact parent spin/momentum sector closure')
        call assemble_fock_hamiltonian_block_interval(model, basis, indices, hrb, hib, &
                                                      ok, invariant)
        call require(ok .and. invariant, 'strict physical sector packet')
        allocate (direct_h(count, count), direct_i(count, count), &
                  direct_b(count, count), direct_ib(count, count))
        call assemble_jellium_pair_sector(modes, length, 1._dp, 1._dp, 2._dp, &
                     100000000_i64, direct_h, direct_i, direct_b, direct_ib, ok, status)
        call require(ok, 'direct physical symmetry constructor')
        do k = 1, count
            do j = 1, count
                call require(abs(direct_h(k, j) - sign(k)*sign(j)*hr(k, j)) < &
                        5.e-15_dp, 'direct constructor matches complete fixed-N parent')
                call require(direct_b(k, j)%lo <= direct_h(k, j) .and. &
                          direct_b(k, j)%hi >= direct_h(k, j), 'direct physical packet')
            end do
        end do
        call quantum_matrix_defect(hr, hi, hrb, hib, defect, ok)
        call require(ok, 'physical coefficient defect')
        re = 0; im = 0; re(1) = 1
        call quantum_initial_norm(re, im, norm0, ok)
        call require(ok, 'identical initial uniform pair state')
        call initialize_quantum_evolution(state, hr, hi, re, im, 1024, h, 1._dp, ok)
        call require(ok, 'shared finite-basis IVP initialization')
        call solve_quantum_evolution(state, ok)
        call require(ok, 'shared finite-basis IVP candidate')
        call independent_hamiltonian(modes, physical_h)
        oracle = cmplx(0._qp, 0._qp, qp); oracle(1) = cmplx(1._qp, 0._qp, qp)
        previous_time = 0
        do t = 1, 4
            call reference_evolve(physical_h, oracle, times(t) - previous_time)
            previous_time = times(t)
            cell = max(1, times(t)*64)
            call sample_quantum_evolution(state, cell, merge(0._dp, 1._dp, t == 1), &
                                          re, im, radius, ok, bounded)
            call require(ok .and. bounded, 'finite-basis continuous trajectory radius')
            call quantum_physical_radius(norm0, radius, defect, real(times(t), dp), &
                                         1._dp, physical_radius, ok)
            call require(ok, 'physical coefficient rounding added once')
            re = sign*re; im = sign*im
            oracle_error = sqrt(sum(abs(cmplx(real(re, qp), real(im, qp), qp) - &
                                        oracle)**2))
            call require(oracle_error <= real(physical_radius, qp), &
                         'independent real128 physical evolution Hilbert containment')
            call prepare_quantum_pair_sector(packet, modes, re, im, norm0, &
                                   physical_radius, 0._dp, [1._dp/16, 1._dp, 1._dp], ok)
            call require(ok, 'same 16-bin spatial observable packet')
            radii(t) = 0
            do b = 1, 16
                x = -0.5_dp + (real(b, dp) - 0.5_dp)/16
                call sample_quantum_observables(packet, x, result, lo, upper, ok)
                call require(ok, 'physical bin sampling')
                call reference_voxel(modes, oracle, real(x, qp), truth)
                call require(all(real(lo, qp) <= truth) .and. &
                             all(real(upper, qp) >= truth), &
                             'independent exact Fourier slab integral containment')
                values(b, t) = result(1); truths(b, t) = truth(1)
                radii(t) = max(radii(t), result(1) - lo(1), upper(1) - result(1))
                difference = abs(real(result(1), qp) - truth(1))
                call require(difference <= real(radii(t), qp), 'reported bin radius')
            end do
            call require(abs(sum(values(:, t))/16 - 1._dp) < 1.e-10_dp, &
                         'uniform partition integrates to normalized pair probability')
        end do
    end subroutine
    subroutine independent_hamiltonian(modes, matrix)
        integer, intent(in) :: modes(:, :)
        real(qp), intent(out) :: matrix(:, :)
        real(qp) :: pi, ell, wave, coupling
        integer :: a, b, q2
        pi = acos(-1._qp); ell = real(length, qp); wave = 2*pi/ell
        coupling = 4*pi*2/(ell**3*wave**2)
        matrix = 0
        do a = 1, size(modes, 2)
            matrix(a, a) = wave**2*sum(modes(:, a)**2)
            do b = 1, a - 1
                q2 = sum((modes(:, a) - modes(:, b))**2)
                call require(q2 > 0, 'distinct reference momenta')
                matrix(a, b) = coupling/real(q2, qp)
                matrix(b, a) = matrix(a, b)
            end do
        end do
    end subroutine
    subroutine reference_evolve(matrix, state, elapsed)
        real(qp), intent(in) :: matrix(:, :)
        complex(qp), intent(inout) :: state(:)
        integer, intent(in) :: elapsed
        complex(qp) :: term(size(state)), total(size(state))
        complex(qp), parameter :: factor = cmplx(0._qp, -1._qp/32, qp)
        integer :: step, order
        ! Independent high-precision Taylor exponential action, not solver input.
        do step = 1, 32*elapsed
            term = state; total = state
            do order = 1, 32
                term = factor*matmul(matrix, term)/real(order, qp)
                total = total + term
            end do
            call require(maxval(abs(term)) < 1.e-32_qp, &
                         'high-precision exponential oracle convergence')
            state = total
        end do
    end subroutine
    subroutine reference_voxel(modes, state, x, truth)
        integer, intent(in) :: modes(:, :)
        complex(qp), intent(in) :: state(:)
        real(qp), intent(in) :: x
        real(qp), intent(out) :: truth(6)
        real(qp) :: pi, sinc, phase, probability(size(state)), norm, g, exchange, z
        integer :: a, b, delta(3)
        pi = acos(-1._qp); probability = abs(state)**2
        norm = sum(probability); g = 0; exchange = 0
        do a = 1, size(state)
            do b = 1, size(state)
                delta = modes(:, a) - modes(:, b)
                ! Full transverse integration kills all nonzero transverse modes.
                if (any(delta(2:3) /= 0)) cycle
                z = real(delta(1), qp)*pi/16
                sinc = 1
                if (delta(1) /= 0) sinc = sin(z)/z
                phase = 2*pi*real(delta(1), qp)*x
                g = g + real(state(a)*conjg(state(b))* &
                             cmplx(cos(phase), sin(phase), qp), qp)*sinc
                exchange = exchange + probability(a)*probability(b)*cos(phase)*sinc
            end do
        end do
        truth(1) = g; truth(2) = g - norm**2
        truth(3) = 0; truth(4) = exchange - norm**2
        truth(5:6) = norm
    end subroutine
    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(*), intent(in) :: message
        if (.not. condition) then
            print *, 'FAIL: '//message
            error stop 1
        end if
    end subroutine
end program
