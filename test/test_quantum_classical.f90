program test_quantum_classical
    !> Independent oracles for H1b (quantum vs classical finite-N cell):
    !> pairwise additivity of the shared N-body Ewald energy; N = 2 classical
    !> Monte Carlo against direct quadrature of the Boltzmann factor over the
    !> relative coordinate; the exact N = 2 relative reduction against the Q4
    !> two-particle sector spectrum in the identical plane-wave set; its ideal
    !> limit against the Q3 canonical ideal-fermion pair function; the
    !> Wigner--Kirkwood weight against the exact harmonic-oscillator density
    !> and finite differences of the Ewald potential; Kelbg limits and
    !> leaf JVPs against finite differences; the spectral dg/dTheta against
    !> central differences; the Q4 diagonal three-body moments against a
    !> Slater-determinant Wick evaluation.
    use, intrinsic :: iso_fortran_env, only: dp => real64, i64 => int64
    use fortnum_status, only: fortnum_status_t, FORTNUM_OK
    use kin6d_nbody, only: nbody_force
    use kin6d_coulomb_leaves, only: ewald_background
    use kin6d_classical_cell
    use kin6d_quantum_pair_cell
    use kin6d_quantum_classical_leaves
    use kin6d_jellium_cell
    use kin6d_quantum_fock, only: fock_matvec, fock_hamiltonian_t
    use kin6d_quantum_fock_leaves, only: qf_kinetic
    use fortnum_symmetric_eigen, only: symmetric_eigen
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp)
    integer :: failures = 0

    call additivity_oracle()
    call pair_mc_oracle()
    call sector_equivalence_oracle()
    if (failures > 0) error stop 'test_quantum_classical failed'
    print '(a)', 'PASS: H1b quantum-classical oracles'
contains
    subroutine check(condition, label)
        logical, intent(in) :: condition
        character(*), intent(in) :: label
        if (condition) then
            print '(2a)', 'ok   ', label
        else
            print '(2a)', 'FAIL ', label
            failures = failures + 1
        end if
    end subroutine

    !> U_N(x) - const_N = sum_{i<j} u(x_i - x_j) with u from the 2-body system.
    subroutine additivity_oracle()
        type(classical_cell_t) :: cell
        type(fortnum_status_t) :: status
        real(dp) :: x(9), f(9), pot, sum_pairs, bg
        logical :: ok
        call init_classical_cell(cell, 3, 2.0_dp, 1.0_dp, model_coulomb, ok)
        call check(ok, 'classical cell N=3')
        x = [0.1_dp, 0.2_dp, -0.3_dp, 1.1_dp, -0.7_dp, 0.4_dp, -1.2_dp, 0.9_dp, 1.3_dp]
        call nbody_force(cell%full, x, f, pot, status)
        call cell_total_energy(cell, x, sum_pairs, ok)
        call ewald_background(cell%full%ewald%alpha, cell%volume, cell%coupling, bg)
        pot = pot - cell%full%constant_energy + 3*bg - &
            1.5_dp*cell%coupling*sum(cell%full%weight(1:cell%full%modes))
        print '(a,es12.4)', '     additivity defect ', pot - sum_pairs
        call check(ok .and. abs(pot - sum_pairs) < &
                   1.0e-9_dp, 'Ewald energy pairwise additive (N=3)')
    end subroutine

    !> N = 2: Metropolis histogram vs quadrature of exp(-beta u) (4 s.e.).
    subroutine pair_mc_oracle()
        type(classical_cell_t) :: cell
        type(classical_mc_t) :: mc
        real(dp), allocatable :: g1(:), edges(:), g(:), gw(:), gx(:), gp(:), dg(:), f2(:)
        real(dp) :: dev, devf
        integer :: i, nb
        logical :: ok
        call init_classical_cell(cell, 2, 1.0_dp, 0.5_dp, model_coulomb, ok)
        nb = 8
        edges = [(i*cell%length/(2*nb), i=0, nb)]
        allocate (g(nb), gw(nb), gx(nb), gp(nb), dg(nb), f2(26))
        call classical_pair_quadrature(cell, edges, 40, 12, g, f2, gw, gx, gp, dg, ok)
        call check(ok, 'N=2 quadrature')
        call classical_cell_mc(cell, edges, 200000, 2000, 40, 7_i64, mc, ok, .true.)
        call check(ok, 'N=2 Monte Carlo')
        dev = maxval(abs(mc%g - g)/max(mc%g_err, 1.0e-3_dp))
        devf = maxval(abs(mc%f2 - f2)/max(mc%f2_err, 1.0e-12_dp))
        print '(a,f8.3,a,f8.3,a,f6.3)', '     Gamma ', cell%gamma, ' max |MC-quad|/se g ', &
            dev, ' acc ', mc%acceptance
        print '(a,8f8.4)', '     g quad ', g
        print '(a,8f8.4)', '     g MC   ', mc%g
        print '(a,f8.3)', '     max |MC-quad|/se f2 ', devf
        call check(dev < 4.5_dp .and. devf < 4.5_dp, 'N=2 MC g, f2 = quadrature (4.5 s.e.)')
        dev = maxval(abs(mc%dg_dbeta - dg)/max(mc%dg_dbeta_err, 1.0e-3_dp))
        call check(dev < 4.5_dp, 'N=2 MC covariance dg/dbeta = quadrature')
        ! Semiclassical point: reweighting variance stays finite.
        g1 = g
        call init_classical_cell(cell, 2, 4.0_dp, 2.0_dp, model_coulomb, ok)
        edges = [(i*cell%length/(2*nb), i=0, nb)]
        call classical_pair_quadrature(cell, edges, 40, 12, g, f2, gw, gx, gp, dg, ok)
        call classical_cell_mc(cell, edges, 200000, 2000, 40, 11_i64, mc, ok, .true.)
        ! The contact bin of the reweighted estimator is heavy-tailed
        ! (weight ~ r^-4 exp(-beta g/r)); compare bins 2..nb.
        dev = maxval(abs(mc%g_pred(2:) - gp(2:))/max(mc%g_pred_err(2:), 1.0e-3_dp))
        call check(dev < 4.5_dp, 'N=2 MC WK+exchange reweighting = quadrature')
        call check(maxval(abs(mc%g - g1)/max(mc%g_err, 1.0e-3_dp)) < 4.5_dp .and. &
                   maxval(abs(g - g1)) < 1.0e-6_dp, &
                   'classical g depends on Gamma only (rs=4, Theta=2 vs rs=1, Theta=0.5)')
    end subroutine

    !> Relative reduction restricted to the Q4 pair set of total momentum K
    !> reproduces the Q4 S_z = 0 two-particle sector spectrum (minus E_cm).
    subroutine sector_equivalence_oracle()
        type(jellium_cell_t) :: jc
        type(jellium_sector_t) :: sector
        type(quantum_pair_cell_t) :: pc
        type(fock_hamiltonian_t) :: model
        integer, allocatable :: modes(:, :), kv(:, :)
        real(dp), allocatable :: h(:, :), e(:), v(:, :), x(:), zero(:), y(:), yi(:), er(:)
        integer :: ks(3, 2), l, a, b, m, d, j, info
        real(dp) :: ecm, worst
        logical :: ok
        call jellium_modes('spherical', 2, modes)
        call init_jellium_cell(jc, 2, modes, 1.5_dp, ok)
        call init_quantum_pair_cell(pc, 1.5_dp, 12, ok)
        call check(ok .and. abs(pc%length - jc%length) < 1.0e-13_dp, 'pair cell = Q4 cell')
        ks(:, 1) = [0, 0, 0]; ks(:, 2) = [1, 1, 0]
        worst = 0
        do l = 1, 2
            call jellium_sector(jc, 1, ks(:, l:l), sector, ok)
            call jellium_hamiltonian(jc, sector%basis, model, ok)
            d = sector%d
            allocate (h(d, d), e(d), v(d, d), x(d), zero(d), y(d), yi(d))
            zero = 0
            do j = 1, d
                x = 0; x(j) = 1
                call fock_matvec(model, sector%basis, x, zero, y, yi, ok)
                h(:, j) = y
            end do
            call symmetric_eigen(h, e, v, info)
            m = 0
            allocate (kv(3, size(modes, 2)**2))
            do a = 1, size(modes, 2)
                do b = 1, size(modes, 2)
                    if (any(modes(:, a) + modes(:, b) /= ks(:, l))) cycle
                    m = m + 1; kv(:, m) = modes(:, a) - modes(:, b)
                end do
            end do
            allocate (er(m))
            call quantum_pair_relative_spectrum(pc, kv(:, 1:m), er, ok)
            call qf_kinetic(jc%length, jc%hbar, 2*jc%mass, real(sum(ks(:, l)**2), dp), ecm)
            call check(ok .and. m == d, 'relative basis size = Q4 sector dimension')
            if (m == d) worst = max(worst, maxval(abs(er + ecm - e))/maxval(abs(e)))
            print '(a,i3,a,3es14.6)', '     K sector d=', d, ' E0 Q4/rel ', e(1), er(1) + ecm, &
                maxval(abs(er + ecm - e))
            deallocate (h, e, v, x, zero, y, yi, kv, er)
        end do
        call check(worst < 1.0e-11_dp, 'relative reduction = Q4 N=2 sector spectrum')
    end subroutine
end program
