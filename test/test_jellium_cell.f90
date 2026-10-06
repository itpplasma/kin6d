program test_jellium_cell
    !> Independent oracles for the interacting jellium cell (Q4):
    !> full parent-basis dense diagonalization (legacy full-scan assembly plus
    !> FortNum dense eigen) versus sector bases with the coupled action and
    !> block Lanczos; the existing Q0 opposite-spin pair-sector constructor;
    !> ideal-gas limits (analytic T=0 response, canonical determinant sums);
    !> Pauli/structural pair-density zeros; Slater and product-state cumulants;
    !> finite-difference coupling derivatives and VJP pairing; a dense
    !> spectral oracle for the certified kick response.
    use, intrinsic :: iso_fortran_env, only: dp => real64, i64 => int64
    use fortnum_symmetric_eigen, only: symmetric_eigen
    use fortnum_interval, only: interval_t
    use kin6d_quantum_fock
    use kin6d_quantum_correlations, only: quantum_cumulant_norms
    use kin6d_jellium_cell
    implicit none
    integer, allocatable :: modes(:, :)
    type(jellium_cell_t) :: cell
    logical :: ok

    call jellium_modes('spherical', 1, modes)
    call check(size(modes, 2) == 7, 'seven |n|^2<=1 modes')
    call init_jellium_cell(cell, 2, modes, 2.0_dp, ok)
    call check(ok .and. abs(cell%rs - 2.0_dp) < 1.0e-13_dp, 'cell from shared rs')
    call full_spectrum_oracle(cell, 'N=2')
    call init_jellium_cell(cell, 3, modes, 2.0_dp, ok)
    call full_spectrum_oracle(cell, 'N=3')
    call pair_sector_oracle()
    call ideal_response_oracles()
    call pair_density_and_cumulant_oracles()
    call derivative_oracles()
    call kick_oracle()
    print '(a)', 'PASS: jellium cell sector, spectrum, response, cumulant, '// &
        'derivative and kick oracles'

contains

    subroutine full_spectrum_oracle(cell, label)
        type(jellium_cell_t), intent(in) :: cell
        character(*), intent(in) :: label
        type(fock_basis_t) :: full
        type(fock_hamiltonian_t) :: model
        type(jellium_spectrum_t) :: spectrum
        real(dp), allocatable :: hr(:, :), hi(:, :), values(:), vectors(:, :), mine(:)
        character(120) :: status
        integer :: m, n, d, p, info, i, j
        integer(i64) :: bytes
        real(dp) :: tmp

        call init_fock_basis(full, cell%m, cell%n, 100000_i64, 2_i64**32, ok, status)
        call check(ok, label//' full parent basis')
        call fock_basis_info(full, m, n, d, p, bytes, ok)
        call jellium_hamiltonian(cell, full, model, ok)
        call check(ok, label//' jellium model')
        allocate (hr(d, d), hi(d, d), values(d), vectors(d, d))
        call assemble_fock_hamiltonian(model, full, hr, hi, ok)
        call check(ok, label//' legacy dense assembly')
        call symmetric_eigen(hr, values, vectors, info)
        call check(info == 0, label//' dense eigen')
        call jellium_build_spectrum(cell, model, 100000, spectrum, ok)
        call check(ok .and. spectrum%count == d, label//' sector spectrum covers parent')
        mine = spectrum%energy(1:spectrum%count)
        do i = 2, d
            tmp = mine(i); j = i - 1
            do while (j >= 1)
                if (mine(j) <= tmp) exit
                mine(j + 1) = mine(j); j = j - 1
            end do
            mine(j + 1) = tmp
        end do
        call check(maxval(abs(mine - values)) < 1.0e-9_dp*maxval(abs(values)), &
                   label//' sector Lanczos spectrum equals parent dense spectrum')
    end subroutine

    subroutine pair_sector_oracle()
        type(jellium_cell_t) :: c2
        type(jellium_sector_t) :: sector
        type(fock_hamiltonian_t) :: model
        integer, allocatable :: m2(:, :)
        real(dp), allocatable :: hr(:, :), hi(:, :), values(:), vectors(:, :), lv(:, :)
        type(interval_t), allocatable :: hrb(:, :), hib(:, :)
        real(dp) :: low(4), res(4)
        integer :: info, k, iterations
        character(120) :: status

        call jellium_modes('spherical', 2, m2)
        call init_jellium_cell(c2, 2, m2, 1.5_dp, ok)
        call jellium_sector(c2, 1, reshape([0, 0, 0], [3, 1]), sector, ok)
        k = size(m2, 2)
        call check(ok .and. sector%d == k, 'N=2 K=0 S_z=0 sector equals pair list')
        call jellium_hamiltonian(c2, sector%basis, model, ok)
        allocate (lv(k, 4))
        call jellium_lowest(model, sector%basis, 4, low, lv, res, ok, iterations)
        call check(ok, 'pair-sector Lanczos')
        allocate (hr(k, k), hi(k, k), hrb(k, k), hib(k, k), values(k), vectors(k, k))
        call assemble_jellium_pair_sector(m2, c2%length, c2%hbar, c2%mass, &
            c2%coupling, 2_i64**30, hr, hi, hrb, hib, ok, status)
        call check(ok, 'Q0 pair-sector constructor')
        call symmetric_eigen(hr, values, vectors, info)
        call check(maxval(abs(low - values(1:4))) < 1.0e-10_dp, &
                   'N=2 sector agrees with Q0 pair-sector constructor')
    end subroutine

    subroutine ideal_response_oracles()
        type(jellium_cell_t) :: c
        type(fock_hamiltonian_t) :: model
        type(jellium_spectrum_t) :: spectrum
        real(dp), allocatable :: w(:), chis(:)
        real(dp) :: chi, can, gc, rcan, rgc, truncation, ground, eq
        integer, parameter :: q(3) = [1, 0, 0]

        call init_jellium_cell(c, 2, modes, 2.0_dp, ok, coupling=0.0_dp)
        call jellium_hamiltonian(c, dummy_basis(c), model, ok)
        call jellium_build_spectrum(c, model, 100000, spectrum, ok)
        call check(ok, 'ideal N=2 spectrum')
        allocate (w(spectrum%count), chis(spectrum%count))
        call jellium_thermal_weights(c, spectrum, 0.0_dp, w, truncation, ground)
        call jellium_static_response(c, model, spectrum, w, q, chi, chis, ok)
        call check(ok, 'ideal T=0 static response solve')
        eq = 0.5_dp*(2*acos(-1.0_dp)/c%length)**2
        call jellium_ideal_response(c, 0.0_dp, q, can, gc, rcan, rgc, ok)
        call check(abs(chi + 4/(c%volume*eq)) < 1.0e-10_dp*abs(chi) .and. &
                   abs(can + 4/(c%volume*eq)) < 1.0e-12_dp*abs(can), &
                   'closed-shell T=0 Lindhard -4/(V e_q)')
        call check(abs(rcan - can) < 1.0e-15_dp, 'RPA equals Lindhard at g=0')
        deallocate (w, chis)
        call init_jellium_cell(c, 3, modes, 2.0_dp, ok, coupling=1.0e-8_dp)
        call jellium_hamiltonian(c, dummy_basis(c), model, ok)
        call jellium_build_spectrum(c, model, 100000, spectrum, ok)
        allocate (w(spectrum%count), chis(spectrum%count))
        call jellium_thermal_weights(c, spectrum, 0.5_dp, w, truncation, ground)
        call check(truncation <= 0, 'complete spectrum has no truncation')
        call jellium_static_response(c, model, spectrum, w, q, chi, chis, ok, 0.0_dp)
        call check(ok, 'thermal N=3 response')
        call jellium_ideal_response(c, 0.5_dp, q, can, gc, rcan, rgc, ok)
        call check(abs(chi - can) < 1.0e-5_dp*abs(can), &
                   'g->0 canonical response equals ideal determinant sum')
        call check(abs(gc - can) > 0, 'grand-canonical Lindhard differs at finite N')
    end subroutine

    function dummy_basis(c) result(basis)
        type(jellium_cell_t), intent(in) :: c
        type(fock_basis_t) :: basis
        type(jellium_sector_t) :: sector
        call jellium_sector(c, c%n/2, reshape([0, 0, 0], [3, 1]), sector, ok)
        if (.not. ok) call jellium_sector(c, c%n, reshape(sum(c%modes(:, 1:c%n), 2), &
                                          [3, 1]), sector, ok)
        call check(ok, 'model basis')
        basis = sector%basis
    end function

    subroutine pair_density_and_cumulant_oracles()
        type(jellium_cell_t) :: c
        type(fock_hamiltonian_t) :: model
        type(jellium_spectrum_t) :: spectrum
        real(dp), allocatable :: w(:), gamma(:, :), gamma2(:, :)
        real(dp) :: edges(6), gs(5), go(5), truncation, ground, tiny(2), t1(1), t2(1)
        real(dp) :: norms(4)
        integer :: b, i

        call init_jellium_cell(c, 2, modes, 2.0_dp, ok, coupling=0.0_dp)
        call jellium_hamiltonian(c, dummy_basis(c), model, ok)
        call jellium_build_spectrum(c, model, 4, spectrum, ok)
        allocate (w(spectrum%count))
        call jellium_thermal_weights(c, spectrum, 0.0_dp, w, truncation, ground)
        allocate (gamma(c%m, c%m), gamma2(c%m*(c%m - 1)/2, c%m*(c%m - 1)/2))
        call jellium_thermal_rdms(c, spectrum, w, gamma, gamma2, ok)
        edges = [(0.1_dp*b*c%length, b=0, 5)]
        call jellium_pair_bins(c, gamma2, edges, gs, go)
        call check(maxval(abs(go - 1)) < 1.0e-12_dp .and. maxval(abs(gs)) < 1.0e-14_dp, &
                   'ideal k=0 singlet: uniform g_opp=1, g_same=0')
        deallocate (w)
        call init_jellium_cell(c, 3, modes, 2.0_dp, ok)
        call jellium_hamiltonian(c, dummy_basis(c), model, ok)
        call jellium_build_spectrum(c, model, 100000, spectrum, ok)
        allocate (w(spectrum%count))
        call jellium_thermal_weights(c, spectrum, 0.5_dp, w, truncation, ground)
        call jellium_thermal_rdms(c, spectrum, w, gamma, gamma2, ok, 0.0_dp)
        call check(abs(sum([(gamma(i, i), i=1, c%m)]) - 3) < 1.0e-12_dp .and. &
                   abs(sum([(gamma2(i, i), i=1, size(gamma2, 1))]) - 3) < 1.0e-12_dp, &
                   'thermal RDM traces N and N(N-1)/2')
        tiny = [0.0_dp, 1.0e-7_dp*c%length]
        call jellium_pair_bins(c, gamma2, tiny, t1, t2)
        call check(abs(t1(1)) < 1.0e-9_dp, 'Pauli exchange hole g_same(0)=0')
        call jellium_pair_bins(c, gamma2, edges, gs, go)
        call check(t2(1) < go(5), 'Coulomb hole: g_opp(0) below its far value')
        call jellium_cumulants(c, spectrum, w, norms, ok, 0.0_dp)
        call check(ok .and. norms(4) > 0 .and. norms(2) > 0, 'interacting cumulants')
        call slater_and_product_cumulants()
    end subroutine

    subroutine slater_and_product_cumulants()
        type(fock_basis_t) :: basis
        integer(i64), allocatable :: states(:)
        real(dp), allocatable :: psi(:), zero(:), gr(:, :), gi(:, :), g2r(:, :), &
                                 g2i(:, :)
        real(dp) :: orb(6, 3), sub(3, 3), norms(4), pair(6), pair_norms(4)
        type(fock_basis_t) :: pair_basis
        integer(i64), allocatable :: pair_states(:)
        real(dp), allocatable :: pair_psi(:)
        integer :: d2
        integer :: d, k, a, cnt, occ(3), i, j, mm, nn, pp
        integer(i64) :: bytes
        character(120) :: status

        call init_fock_basis(basis, 6, 3, 1000_i64, 2_i64**24, ok, status)
        call fock_basis_info(basis, mm, nn, d, pp, bytes, ok)
        allocate (states(d), psi(d), zero(d), gr(6, 6), gi(6, 6), g2r(15, 15), &
                  g2i(15, 15))
        call get_fock_states(basis, states, ok)
        do j = 1, 3
            do i = 1, 6
                orb(i, j) = cos(0.9_dp*i*j + 0.2_dp*j*j)
            end do
            do k = 1, j - 1
                orb(:, j) = orb(:, j) - dot_product(orb(:, k), orb(:, j))*orb(:, k)
            end do
            orb(:, j) = orb(:, j)/norm2(orb(:, j))
        end do
        do k = 1, d
            cnt = 0
            do a = 1, 6
                if (btest(states(k), a - 1)) then
                    cnt = cnt + 1; occ(cnt) = a
                end if
            end do
            sub = orb(occ, :)
            psi(k) = sub(1, 1)*(sub(2, 2)*sub(3, 3) - sub(2, 3)*sub(3, 2)) - &
                     sub(1, 2)*(sub(2, 1)*sub(3, 3) - sub(2, 3)*sub(3, 1)) + &
                     sub(1, 3)*(sub(2, 1)*sub(3, 2) - sub(2, 2)*sub(3, 1))
        end do
        zero = 0
        call cumulant_of(basis, psi, norms)
        call check(norms(2) < 1.0e-13_dp .and. norms(4) < 1.0e-13_dp .and. &
                   norms(3) > 0.1_dp, 'real Slater determinant: zero cumulants')
        pair = [0.6_dp, -0.3_dp, 0.25_dp, 0.5_dp, -0.35_dp, 0.2_dp]
        pair = pair/norm2(pair)
        psi = 0; cnt = 0
        do i = 1, 3
            do j = i + 1, 4
                cnt = cnt + 1
                do k = 1, d
                    if (states(k) == ibset(ibset(ibset(0_i64, i - 1), j - 1), 4)) &
                        psi(k) = pair(cnt)
                end do
            end do
        end do
        call cumulant_of(basis, psi, norms)
        ! Cumulants are additive over disjoint independent subsystems: the
        ! spectator determinant adds nothing, so lambda2 and lambda3 equal
        ! those of the fixed-N=2 pair state alone (whose Gamma3 vanishes but
        ! whose lambda3 does not: a fixed particle number is itself a
        ! correlation).
        call init_fock_basis(pair_basis, 6, 2, 1000_i64, 2_i64**24, ok, status)
        call fock_basis_info(pair_basis, mm, nn, d2, pp, bytes, ok)
        allocate (pair_states(d2), pair_psi(d2))
        call get_fock_states(pair_basis, pair_states, ok)
        pair_psi = 0; cnt = 0
        do i = 1, 3
            do j = i + 1, 4
                cnt = cnt + 1
                pair_psi(findloc(pair_states, ibset(ibset(0_i64, i - 1), j - 1), 1)) = &
                    pair(cnt)
            end do
        end do
        call cumulant_of(pair_basis, pair_psi, pair_norms)
        call check(norms(2) > 0.1_dp .and. norms(4) > 0.01_dp .and. &
                   abs(norms(2) - pair_norms(2)) < 1.0e-13_dp .and. &
                   abs(norms(4) - pair_norms(4)) < 1.0e-13_dp, &
                   'cumulant additivity: correlated pair times spectator')
    end subroutine

    subroutine cumulant_of(basis, psi, result)
        type(fock_basis_t), intent(in) :: basis
        real(dp), intent(in) :: psi(:)
        real(dp), intent(out) :: result(4)
        integer(i64), allocatable :: tk(:), uk(:)
        real(dp), allocatable :: zero(:), gr(:, :), gi(:, :), g2r(:, :), g2i(:, :), &
                                 vr(:), vi(:)
        integer :: cnt
        allocate (zero(size(psi)), gr(6, 6), gi(6, 6), g2r(15, 15), g2i(15, 15))
        zero = 0
        call fock_rdms(basis, psi, zero, gr, gi, g2r, g2i, ok)
        call fock_rdm3_entries(basis, psi, zero, tk, uk, vr, vi, cnt, ok)
        call quantum_cumulant_norms(6, gr, g2r, tk(1:cnt), uk(1:cnt), vr(1:cnt), &
            result(1), result(2), result(3), result(4), ok)
        call check(ok, 'cumulant norms')
    end subroutine

    subroutine ground_pair(c, g, sector, model, psi, energy)
        type(jellium_cell_t), intent(in) :: c
        real(dp), intent(in) :: g
        type(jellium_sector_t), intent(out) :: sector
        type(fock_hamiltonian_t), intent(out) :: model
        real(dp), allocatable, intent(out) :: psi(:)
        real(dp), intent(out) :: energy
        real(dp) :: val(2), res(2)
        real(dp), allocatable :: vec(:, :)
        integer :: iterations
        call jellium_sector(c, 1, reshape([0, 0, 0], [3, 1]), sector, ok)
        call jellium_hamiltonian(c, sector%basis, model, ok, g)
        allocate (vec(sector%d, 2))
        call jellium_lowest(model, sector%basis, 2, val, vec, res, ok, iterations)
        call check(ok .and. val(2) - val(1) > 1.0e-3_dp, 'nondegenerate N=2 ground state')
        psi = vec(:, 1); energy = val(1)
    end subroutine

    subroutine derivative_oracles()
        type(jellium_cell_t) :: c
        type(jellium_sector_t) :: sector
        type(fock_hamiltonian_t) :: model, interaction, mp, mm
        real(dp), allocatable :: psi(:), pp(:), pm(:), zero(:), gr(:, :), gi(:, :), &
                                 g2r(:, :), g2i(:, :)
        real(dp) :: energy, ep, em, de, edges(4), gs(3), go(3), dgs(3), dgo(3)
        real(dp) :: gsp(3), gop(3), gsm(3), gom(3), h, gbar, bs(3), bo(3)
        integer :: iterations, p, b

        call init_jellium_cell(c, 2, modes, 2.0_dp, ok)
        call ground_pair(c, c%coupling, sector, model, psi, energy)
        call jellium_interaction_model(c, sector%basis, interaction, ok)
        edges = [(0.15_dp*b*c%length, b=0, 3)]
        call jellium_coupling_derivative(c, sector%basis, model, interaction, psi, &
            energy, edges, de, gs, go, dgs, dgo, ok, iterations)
        call check(ok, 'equation-level eigen derivative')
        h = 1.0e-4_dp
        call ground_pair(c, c%coupling + h, sector, mp, pp, ep)
        call ground_pair(c, c%coupling - h, sector, mm, pm, em)
        p = c%m*(c%m - 1)/2
        allocate (zero(sector%d), gr(c%m, c%m), gi(c%m, c%m), g2r(p, p), g2i(p, p))
        zero = 0
        call fock_rdms(sector%basis, pp, zero, gr, gi, g2r, g2i, ok)
        call jellium_pair_bins(c, g2r, edges, gsp, gop)
        call fock_rdms(sector%basis, pm, zero, gr, gi, g2r, g2i, ok)
        call jellium_pair_bins(c, g2r, edges, gsm, gom)
        call check(abs(de - (ep - em)/(2*h)) < 1.0e-7_dp*abs(de), &
                   'Hellmann-Feynman dE/dg vs central difference')
        call check(maxval(abs(dgo - (gop - gom)/(2*h))) < 1.0e-6_dp*maxval(abs(dgo)), &
                   'eigen-derivative dg_opp/dg vs central difference')
        call check(maxval(abs(dgo)) > 1.0e-4_dp, 'nontrivial pair derivative')
        bs = [0.3_dp, -0.2_dp, 0.1_dp]; bo = [0.7_dp, 0.4_dp, -0.5_dp]
        call jellium_coupling_vjp(c, sector%basis, model, interaction, psi, energy, &
                                  edges, 0.25_dp, bs, bo, gbar, ok)
        call check(ok .and. abs(gbar - (0.25_dp*de + sum(bs*dgs) + sum(bo*dgo))) < &
                   1.0e-8_dp*abs(gbar), 'adjoint VJP pairs with tangent JVP')
    end subroutine

    subroutine kick_oracle()
        type(jellium_cell_t) :: c
        type(jellium_sector_t) :: sector, union
        type(fock_hamiltonian_t) :: model, full_model, density
        type(fock_basis_t) :: full
        real(dp), allocatable :: psi(:), times(:), response(:), band(:)
        real(dp), allocatable :: hr(:, :), hi(:, :), ar(:, :), ai(:, :), sub(:, :)
        real(dp), allocatable :: values(:), vectors(:, :), full_psi(:), phi(:), amp(:)
        integer(i64), allocatable :: full_states(:), sector_states(:), union_states(:)
        real(dp) :: energy, tc, tz, cr, ci, exact
        integer, parameter :: steps = 200, q(3) = [1, 0, 0]
        integer :: du, k, i, d, mm, nn, pp, info, ks(3, 2)
        integer(i64) :: bytes
        character(120) :: status

        call init_jellium_cell(c, 2, modes, 2.0_dp, ok)
        call ground_pair(c, c%coupling, sector, model, psi, energy)
        allocate (times(0:steps), response(0:steps), band(0:steps))
        call jellium_kick_evolution(c, model, sector, psi, energy, q, steps, 0.01_dp, &
                                    times, response, band, ok, du, tc, tz)
        call check(ok, 'certified kick evolution')
        call check(band(steps) < 1.0e-6_dp, 'finite-basis response band below 1e-6')
        call init_fock_basis(full, c%m, c%n, 100000_i64, 2_i64**32, ok, status)
        call fock_basis_info(full, mm, nn, d, pp, bytes, ok)
        allocate (hr(d, d), hi(d, d), ar(d, d), ai(d, d), full_states(d), full_psi(d))
        call jellium_hamiltonian(c, full, full_model, ok)
        call assemble_fock_hamiltonian(full_model, full, hr, hi, ok)
        call jellium_density_model(c, full, q, density, ok)
        call assemble_fock_hamiltonian(density, full, ar, ai, ok)
        call check(ok, 'legacy dense density operator')
        call get_fock_states(full, full_states, ok)
        allocate (sector_states(sector%d))
        call get_fock_states(sector%basis, sector_states, ok)
        full_psi = 0
        do k = 1, sector%d
            full_psi(findloc(full_states, sector_states(k), 1)) = psi(k)
        end do
        ks(:, 1) = q; ks(:, 2) = -q
        call jellium_sector(c, 1, ks, union, ok)
        allocate (union_states(union%d), sub(union%d, union%d), values(union%d), &
                  vectors(union%d, union%d), phi(union%d), amp(union%d))
        call get_fock_states(union%basis, union_states, ok)
        call check(union%d == du, 'union dimension')
        do k = 1, union%d
            i = findloc(full_states, union_states(k), 1)
            phi(k) = dot_product(ar(i, :), full_psi)
            do du = 1, union%d
                sub(k, du) = hr(i, findloc(full_states, union_states(du), 1))
            end do
        end do
        call symmetric_eigen(sub, values, vectors, info)
        amp = matmul(transpose(vectors), phi)
        do k = 0, steps
            cr = sum(amp**2*cos(values*times(k))); ci = -sum(amp**2*sin(values*times(k)))
            exact = 2*(sin(energy*times(k))*cr + cos(energy*times(k))*ci)
            call check(abs(exact - response(k)) <= band(k) + 1.0e-12_dp, &
                       'dense spectral oracle inside certified response band')
        end do
        call check(maxval(abs(response)) > 1.0e-2_dp, 'nontrivial density-mode response')
    end subroutine

    subroutine check(condition, message)
        logical, intent(in) :: condition
        character(*), intent(in) :: message
        if (.not. condition) then
            print '(a)', 'FAIL: '//message
            error stop 1
        end if
    end subroutine
end program test_jellium_cell
