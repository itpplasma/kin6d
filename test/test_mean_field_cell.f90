program test_mean_field_cell
    !> Independent oracles for the Q5 mean-field models of the jellium cell:
    !> the mean-field energy functional against <Phi|H|Phi> of the exact
    !> Fock-space Hamiltonian (matrix-free engine action) for converged and
    !> arbitrary symmetry-broken determinants; N=2 opposite-spin UHF against a
    !> brute-force angle minimization of <Phi|H|Phi> in a two-mode cell;
    !> Koopmans and total-energy identities; the g -> 0 ideal limit (energy,
    !> canonical Lindhard response); TDHF without exchange against the
    !> same-cell RPA; variational and Temple ordering against exact sector
    !> spectra; Hellmann-Feynman against finite differences; vanishing Slater
    !> cumulants; Mermin limits (ideal grand-canonical gas, Theta -> 0).
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_quantum_fock, only: fock_hamiltonian_t, fock_matvec
    use kin6d_jellium_cell
    use kin6d_mean_field_cell
    implicit none
    integer, allocatable :: modes(:, :)

    call jellium_modes('spherical', 1, modes)
    call functional_oracles()
    call two_mode_minimization_oracle()
    call identity_oracles()
    call ideal_limit_oracles()
    call variational_and_temple_oracles()
    call derivative_oracle()
    call cumulant_oracle()
    call mermin_oracles()
    print '(a)', 'PASS: mean-field functional, minimization, Koopmans, ideal/RPA '// &
        'limits, variational/Temple, Hellmann-Feynman, cumulant and Mermin oracles'

contains

    subroutine check(condition, label)
        logical, intent(in) :: condition
        character(*), intent(in) :: label
        if (.not. condition) then
            print '(2a)', 'FAIL: ', label
            error stop 1
        end if
    end subroutine

    real(dp) function rel(a, b)
        real(dp), intent(in) :: a, b
        rel = abs(a - b)/max(1.0e-300_dp, abs(b))
    end function

    !> <Phi|H|Phi> of a Slater state in the union basis of its spin filling.
    real(dp) function parent_energy(cell, state) result(e)
        type(jellium_cell_t), intent(in) :: cell
        type(mf_state_t), intent(in) :: state
        type(jellium_sector_t) :: sector
        type(fock_hamiltonian_t) :: model
        integer, allocatable :: ks(:, :)
        real(dp), allocatable :: psi(:), hpsi(:), zero(:), scratch(:)
        real(dp) :: norm
        logical :: ok
        call jellium_sector_momenta(cell, state%nup, ks)
        call jellium_sector(cell, state%nup, ks, sector, ok)
        call check(ok, 'union sector')
        call jellium_hamiltonian(cell, sector%basis, model, ok)
        call mf_slater_vector(cell, state, sector%basis, psi, norm, ok)
        call check(ok .and. abs(norm - 1) < 1.0e-12_dp, 'Slater vector normalized')
        allocate (hpsi(size(psi)), zero(size(psi)), scratch(size(psi)))
        zero = 0
        call fock_matvec(model, sector%basis, psi, zero, hpsi, scratch, ok)
        call check(ok, 'parent action')
        e = dot_product(psi, hpsi)
    end function

    subroutine functional_oracles()
        type(jellium_cell_t) :: cell
        type(mf_interaction_t) :: inter
        type(mf_state_t) :: state, broken
        real(dp) :: ek, eh, ex
        integer :: ns, nc, n
        logical :: ok
        do n = 2, 3
            call init_jellium_cell(cell, n, modes, 2.0_dp, ok)
            call mf_interaction(cell, inter)
            call mf_ground(cell, inter, .false., .true., 4, state, ok, ns, nc)
            call check(ok .and. state%converged, 'UHF converged')
            call check(rel(state%energy, parent_energy(cell, state)) < 1.0e-12_dp, &
                       'HF energy = <Phi|H|Phi> (converged)')
            ! An arbitrary symmetry-broken determinant (not stationary).
            broken = state
            call mf_start_orbitals(inter, 3, 2*cell%fermi, broken%cup)
            call mf_start_orbitals(inter, 7, 2*cell%fermi, broken%cdn)
            call mf_energy(inter, broken, ek, eh, ex)
            call check(rel(ek + eh + ex, parent_energy(cell, broken)) < 1.0e-12_dp, &
                       'functional = <Phi|H|Phi> (broken symmetry)')
            call check(abs(eh) > 1.0e-6_dp, 'broken state has Hartree energy')
        end do
    end subroutine

    !> N = 2, S_z = 0 in modes {0, (100)}: orbitals cos t|0> + sin t|q> per
    !> spin. Grid plus alternating golden-section refinement of <Phi|H|Phi>.
    subroutine two_mode_minimization_oracle()
        type(jellium_cell_t) :: cell
        type(mf_interaction_t) :: inter
        type(mf_state_t) :: state, probe, trial
        real(dp) :: cu(2, 2), cd(2, 2)
        integer :: two(3, 2), i, j, nc, pass, c
        real(dp) :: best, t(2), e, couplings(2), step
        logical :: ok
        two = reshape([0, 0, 0, 1, 0, 0], [3, 2])
        couplings = [1.0_dp, 30.0_dp]
        step = acos(-1.0_dp)/180
        do c = 1, 2
            call init_jellium_cell(cell, 2, two, 2.0_dp, ok, couplings(c))
            call mf_interaction(cell, inter)
            ! Opposite-spin filling only (the S_z = 1 triplet is a different sector).
            nc = 0
            do i = 0, 7
                call mf_start_orbitals(inter, i, cell%fermi, cu)
                call mf_start_orbitals(inter, merge(0, i + 101, i == 0), cell%fermi, cd)
                call mf_scf(cell, inter, 1, .false., .true., cu, cd, trial, ok)
                if (.not. ok) cycle
                nc = nc + 1
                if (nc == 1) state = trial
                if (trial%energy < state%energy) state = trial
            end do
            call check(nc > 0, 'two-mode UHF')
            probe = state
            probe%fup = [1, 0]; probe%fdn = [1, 0]; probe%nup = 1; probe%ndn = 1
            best = huge(1.0_dp); t = 0
            do i = 0, 179
                do j = 0, 179
                    e = angle_energy(cell, probe, [i*step, j*step])
                    if (e < best) then
                        best = e; t = [i*step, j*step]
                    end if
                end do
            end do
            do pass = 1, 20
                call golden(cell, probe, t, 1); call golden(cell, probe, t, 2)
            end do
            best = angle_energy(cell, probe, t)
            call check(abs(state%energy - best) < 1.0e-9_dp*max(abs(best), cell%fermi), &
                       'UHF = brute-force minimum of <Phi|H|Phi>')
            ! Strong coupling: spin-density-wave UHF minimum (uniform total
            ! density, so zero Hartree energy; modulated spin density).
            if (c == 2) call check(abs(state%cup(2, 1)) > 1.0e-3_dp .and. &
                abs(state%hartree) < 1.0e-12_dp, 'strong coupling: broken-symmetry UHF')
        end do
    end subroutine

    real(dp) function angle_energy(cell, probe, t) result(val)
        type(jellium_cell_t), intent(in) :: cell
        type(mf_state_t), intent(inout) :: probe
        real(dp), intent(in) :: t(2)
        probe%cup = reshape([cos(t(1)), sin(t(1)), -sin(t(1)), cos(t(1))], [2, 2])
        probe%cdn = reshape([cos(t(2)), sin(t(2)), -sin(t(2)), cos(t(2))], [2, 2])
        val = parent_energy(cell, probe)
    end function

    subroutine golden(cell, probe, t, which)
        type(jellium_cell_t), intent(in) :: cell
        type(mf_state_t), intent(inout) :: probe
        real(dp), intent(inout) :: t(2)
        integer, intent(in) :: which
        real(dp) :: lo, hi, x1, x2, f1, f2, g
        integer :: k
        g = 0.5_dp*(sqrt(5.0_dp) - 1)
        lo = t(which) - 0.02_dp; hi = lo + 0.04_dp
        x1 = hi - g*(hi - lo); x2 = lo + g*(hi - lo)
        f1 = along(cell, probe, t, which, x1); f2 = along(cell, probe, t, which, x2)
        do k = 1, 60
            if (f1 < f2) then
                hi = x2; x2 = x1; f2 = f1; x1 = hi - g*(hi - lo); f1 = along(cell, probe, t, which, x1)
            else
                lo = x1; x1 = x2; f1 = f2; x2 = lo + g*(hi - lo); f2 = along(cell, probe, t, which, x2)
            end if
        end do
        t(which) = 0.5_dp*(lo + hi)
    end subroutine

    real(dp) function along(cell, probe, t, which, x)
        type(jellium_cell_t), intent(in) :: cell
        type(mf_state_t), intent(inout) :: probe
        real(dp), intent(in) :: t(2), x
        integer, intent(in) :: which
        real(dp) :: s(2)
        s = t; s(which) = x
        along = angle_energy(cell, probe, s)
    end function

    !> Koopmans (frozen-orbital N-1 energies) and E = (1/2) sum_occ (eps + t).
    subroutine identity_oracles()
        type(jellium_cell_t) :: cell
        type(mf_interaction_t) :: inter
        type(mf_state_t) :: state
        real(dp) :: half
        integer :: ns, nc, j
        logical :: ok
        call init_jellium_cell(cell, 3, modes, 4.0_dp, ok)
        call mf_interaction(cell, inter)
        call mf_ground(cell, inter, .false., .true., 4, state, ok, ns, nc)
        call check(ok, 'N=3 UHF')
        half = 0
        do j = 1, state%nup
            call check(abs(state%energy - mf_frozen_energy(inter, state, j, 1) - &
                state%eup(j)) < 1.0e-12_dp*abs(state%energy), 'Koopmans up')
            half = half + state%eup(j) + sum(inter%kinetic*state%cup(:, j)**2)
        end do
        do j = 1, state%ndn
            call check(abs(state%energy - mf_frozen_energy(inter, state, j, -1) - &
                state%edn(j)) < 1.0e-12_dp*abs(state%energy), 'Koopmans down')
            half = half + state%edn(j) + sum(inter%kinetic*state%cdn(:, j)**2)
        end do
        call check(rel(0.5_dp*half, state%energy) < 1.0e-12_dp, 'E = (eps + t)/2 sum')
    end subroutine

    !> g -> 0: HF equals the ideal ground energy and TDHF the canonical
    !> same-cell Lindhard response; exchange off: TDH equals same-cell RPA.
    subroutine ideal_limit_oracles()
        type(jellium_cell_t) :: cell
        type(mf_interaction_t) :: inter
        type(mf_state_t) :: state, hartree
        real(dp) :: chi, hmin, chi_can, chi_gc, rpa_can, rpa_gc, ideal
        integer :: ns, nc, nexc, n
        logical :: ok
        do n = 2, 3
            call init_jellium_cell(cell, n, modes, 2.0_dp, ok, 1.0e-10_dp)
            call mf_interaction(cell, inter)
            call mf_ground(cell, inter, .false., .true., 3, state, ok, ns, nc)
            ideal = inter%kinetic(1)*2 + merge(inter%kinetic(2), 0.0_dp, n == 3)
            call check(ok .and. abs(state%energy - ideal) < 1.0e-9_dp*cell%fermi, &
                       'g->0 HF = ideal ground')
        end do
        call init_jellium_cell(cell, 2, modes, 2.0_dp, ok, 1.0e-10_dp)
        call mf_interaction(cell, inter)
        call mf_ground(cell, inter, .false., .true., 1, state, ok, ns, nc)
        call mf_static_response(cell, inter, state, [1, 0, 0], chi, hmin, ok, nexc)
        call jellium_ideal_response(cell, 0.0_dp, [1, 0, 0], chi_can, chi_gc, &
                                    rpa_can, rpa_gc, ok)
        call check(ok .and. rel(chi, chi_can) < 1.0e-8_dp, 'g->0 TDHF = Lindhard')
        call init_jellium_cell(cell, 3, modes, 2.0_dp, ok)
        call mf_interaction(cell, inter)
        call mf_ground(cell, inter, .false., .false., 1, hartree, ok, ns, nc)
        call check(ok .and. abs(hartree%hartree) < 1.0e-14_dp .and. &
                   abs(hartree%energy - hartree%kinetic) < 1.0e-14_dp, &
                   'homogeneous Hartree potential vanishes (N=3)')
        call init_jellium_cell(cell, 2, modes, 2.0_dp, ok)
        call mf_interaction(cell, inter)
        call mf_ground(cell, inter, .false., .false., 1, hartree, ok, ns, nc)
        call mf_static_response(cell, inter, hartree, [1, 0, 0], chi, hmin, ok, nexc)
        call jellium_ideal_response(cell, 0.0_dp, [1, 0, 0], chi_can, chi_gc, &
                                    rpa_can, rpa_gc, ok)
        call check(ok .and. rel(chi, rpa_can) < 1.0e-12_dp, 'TDH (no exchange) = RPA')
        call mf_ground(cell, inter, .false., .true., 1, state, ok, ns, nc)
        call mf_static_response(cell, inter, state, [1, 0, 0], chi, hmin, ok, nexc)
        call check(ok .and. hmin > 0 .and. abs(chi) > abs(rpa_can), &
                   'TDHF stable, exchange enhances |chi| over RPA')
    end subroutine

    !> E_HF >= E0 (exact, all sectors); Temple lower <= sector E0 <= rho.
    subroutine variational_and_temple_oracles()
        type(jellium_cell_t) :: cell
        type(mf_interaction_t) :: inter
        type(mf_state_t) :: state, rhf
        type(fock_hamiltonian_t) :: model
        type(jellium_spectrum_t) :: spectrum
        type(jellium_sector_t) :: sector
        type(mf_estimate_t) :: est
        real(dp), allocatable :: psi(:)
        real(dp) :: rs(2), e0, norm
        integer :: n, r, ns, nc, nk
        logical :: ok
        rs = [1.0_dp, 4.0_dp]
        do n = 2, 3
            do r = 1, 2
                call init_jellium_cell(cell, n, modes, rs(r), ok)
                call mf_interaction(cell, inter)
                call mf_ground(cell, inter, .false., .true., 4, state, ok, ns, nc)
                call mf_slater_sector(cell, state, sector, psi, norm, nk, ok)
                call check(ok .and. abs(norm - 1) < 1.0e-12_dp, 'Slater sector')
                call jellium_hamiltonian(cell, sector%basis, model, ok)
                call jellium_build_spectrum(cell, model, 4, spectrum, ok)
                e0 = minval(spectrum%energy(1:spectrum%count))
                call check(state%energy >= e0, 'variational E_HF >= E0')
                call mf_residual(model, sector, psi, est, ok)
                call check(ok .and. rel(est%rho, state%energy) < 1.0e-12_dp, 'rho = E_HF')
                call check(est%e0 >= e0 - 1.0e-10_dp .and. est%e0 <= est%rho, &
                           'sector ground between E0 and rho')
                if (est%temple_valid) call check(est%temple <= est%e0 + est%r0, &
                                                 'Temple lower bound')
                if (n == 2 .and. r == 1) call check(est%temple_valid, 'Temple applies')
                if (n == 2) then
                    call mf_ground(cell, inter, .true., .true., 2, rhf, ok, ns, nc)
                    call check(ok .and. state%energy <= rhf%energy + 1.0e-12_dp, &
                               'UHF <= RHF')
                end if
            end do
        end do
    end subroutine

    subroutine derivative_oracle()
        type(jellium_cell_t) :: cell
        type(mf_interaction_t) :: inter
        type(mf_state_t) :: state, sp, sm
        real(dp) :: g, h, de, fd
        integer :: ns, nc
        logical :: ok
        call init_jellium_cell(cell, 3, modes, 2.0_dp, ok)
        g = cell%coupling; h = 1.0e-4_dp
        call mf_interaction(cell, inter)
        call mf_ground(cell, inter, .false., .true., 1, state, ok, ns, nc)
        call check(ok .and. abs(state%exchange_energy) > 1.0e-3_dp, 'N=3 exchange')
        de = mf_coupling_derivative(cell, state)
        call mf_interaction(cell, inter, g + h)
        call mf_scf(cell, inter, state%nup, .false., .true., state%cup, state%cdn, sp, ok)
        call mf_interaction(cell, inter, g - h)
        call mf_scf(cell, inter, state%nup, .false., .true., state%cup, state%cdn, sm, ok)
        fd = (sp%energy - sm%energy)/(2*h)
        call check(rel(de, fd) < 1.0e-7_dp, 'Hellmann-Feynman dE_HF/dg = FD')
    end subroutine

    !> All connected cumulants of an arbitrary Slater state vanish.
    subroutine cumulant_oracle()
        type(jellium_cell_t) :: cell
        type(mf_interaction_t) :: inter
        type(mf_state_t) :: state
        type(jellium_sector_t) :: sector
        type(jellium_spectrum_t) :: spectrum
        real(dp), allocatable :: psi(:)
        real(dp) :: norms(4), norm
        integer :: ns, nc, nk
        logical :: ok
        call init_jellium_cell(cell, 3, modes, 2.0_dp, ok)
        call mf_interaction(cell, inter)
        call mf_ground(cell, inter, .false., .true., 1, state, ok, ns, nc)
        call mf_start_orbitals(inter, 5, cell%fermi, state%cup)
        call mf_start_orbitals(inter, 9, cell%fermi, state%cdn)
        call mf_slater_sector(cell, state, sector, psi, norm, nk, ok)
        call check(ok .and. nk > 1, 'broken determinant spans several K')
        call mf_as_spectrum(sector, psi, 0.0_dp, spectrum)
        call jellium_cumulants(cell, spectrum, [1.0_dp], norms, ok, 1.0e-12_dp)
        call check(ok .and. norms(2) < 1.0e-12_dp*norms(1) .and. &
                   norms(4) < 1.0e-12_dp*max(norms(3), 1.0_dp), 'Slater cumulants vanish')
    end subroutine

    !> Mermin: ideal grand-canonical gas at g = 0; Theta -> 0 HF limit.
    subroutine mermin_oracles()
        type(jellium_cell_t) :: cell
        type(mf_interaction_t) :: inter
        type(mf_state_t) :: thermal, ground, up, down
        real(dp) :: beta, lo, hi, mu, f(64), e, h
        integer :: it, ns, nc
        logical :: ok
        call init_jellium_cell(cell, 2, modes, 2.0_dp, ok, 0.0_dp)
        call mf_interaction(cell, inter)
        call mf_mermin(cell, inter, 0.5_dp, .true., thermal, ok)
        call check(ok, 'Mermin converged (g=0)')
        beta = 1/(0.5_dp*cell%fermi); lo = -50*cell%fermi; hi = 50*cell%fermi
        do it = 1, 200
            mu = 0.5_dp*(lo + hi)
            f(1:inter%nmodes) = 1/(1 + exp(beta*(inter%kinetic - mu)))
            if (2*sum(f(1:inter%nmodes)) > 2) then
                hi = mu
            else
                lo = mu
            end if
        end do
        e = 2*sum(f(1:inter%nmodes)*inter%kinetic)
        call check(rel(thermal%energy, e) < 1.0e-10_dp, 'Mermin g=0 = ideal GC energy')
        call init_jellium_cell(cell, 2, modes, 2.0_dp, ok)
        call mf_interaction(cell, inter)
        call mf_ground(cell, inter, .true., .true., 1, ground, ok, ns, nc)
        call mf_mermin(cell, inter, 0.01_dp, .true., thermal, ok)
        call check(ok .and. abs(thermal%energy - ground%energy) < 1.0e-8_dp*cell%fermi, &
                   'Mermin Theta->0 = RHF')
        ! Envelope identity at a stationary Mermin state: dF/dT = -S (fixed N).
        call mf_mermin(cell, inter, 0.5_dp, .true., thermal, ok)
        h = 1.0e-4_dp
        call mf_mermin(cell, inter, 0.5_dp + h, .true., up, ok)
        call mf_mermin(cell, inter, 0.5_dp - h, .true., down, ok)
        call check(ok .and. rel((up%free - down%free)/(2*h*cell%fermi), -thermal%entropy) &
                   < 1.0e-6_dp, 'Mermin dF/dT = -S')
    end subroutine
end program test_mean_field_cell
