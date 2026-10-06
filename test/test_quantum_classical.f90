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
    use kin6d_quantum_classical_map
    use kin6d_quantum_fock, only: fock_matvec, fock_hamiltonian_t
    use kin6d_quantum_fock_leaves, only: qf_kinetic
    use fortnum_symmetric_eigen, only: symmetric_eigen
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp)
    integer :: failures = 0

    call additivity_oracle()
    call pair_mc_oracle()
    call sector_equivalence_oracle()
    call ideal_limit_oracle()
    call wigner_kirkwood_oracle()
    call kelbg_oracle()
    call pair_derivative_oracle()
    call moments_oracle()
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

    !> Ideal limit of the relative reduction (g = 0) = Q3 canonical ideal
    !> fermions over all S_z (independent finite-N recursion).
    subroutine ideal_limit_oracle()
        type(quantum_pair_cell_t) :: pc
        real(dp), allocatable :: edges(:), gs(:), gss(:), go(:), is(:), iss(:), io(:)
        real(dp) :: f2(26), theta(2), worst
        integer, allocatable :: q(:, :)
        integer :: i, j, nb
        logical :: ok
        theta = [0.5_dp, 2.0_dp]
        call cell_q_set(q)
        worst = 0
        do j = 1, 2
            call init_quantum_pair_cell(pc, 1.0_dp, quantum_pair_cutoff(1.0_dp, &
                                        theta(j), 40.0_dp, 16), ok, 0.0_dp)
            nb = 10
            edges = [(i*pc%length/(2*nb), i=0, nb)]
            allocate (gs(nb), gss(nb), go(nb), is(nb), iss(nb), io(nb))
            call quantum_pair_thermal(pc, theta(j), edges, q, gs, gss, go, f2, ok)
            call ideal_cell_pair(2, theta(j), pc%length, edges, is, iss, io, ok)
            worst = max(worst, maxval(abs(gs - is)), maxval(abs(gss - iss)), &
                        maxval(abs(go - io)))
            print '(a,f5.2,a,3f9.5,a,es10.2)', '     ideal Theta', theta(j), ' g(0..) ', &
                gs(1:3), ' max dev Q3 ', worst
            deallocate (gs, gss, go, is, iss, io)
        end do
        call check(worst < 1.0e-9_dp, 'N=2 ideal relative reduction = Q3 canonical (all S_z)')
    end subroutine

    !> Leading Wigner--Kirkwood weight vs the exact 1D harmonic-oscillator
    !> diagonal density (hbar = m = 1): error O((beta omega)^4).
    subroutine wigner_kirkwood_oracle()
        type(classical_cell_t) :: cell
        real(dp) :: y(2), err(2), beta, x, w, exact, wk0, ex0, d(3), u0, up, um, &
            f(3), lap, h, fd(3)
        integer :: k, i, a
        logical :: ok
        y = [0.2_dp, 0.1_dp]
        do k = 1, 2
            beta = y(k)       ! omega = 1
            err(k) = 0
            call qcm_wigner_kirkwood(beta, 1.0_dp, 1.0_dp, 1.0_dp, 0.0_dp, wk0)
            ex0 = sqrt(beta/sinh(beta))
            do i = -20, 20
                x = 3*i/(20*sqrt(beta))   ! three thermal widths
                call qcm_wigner_kirkwood(beta, 1.0_dp, 1.0_dp, 1.0_dp, x**2, w)
                exact = ex0*exp(-x**2*(tanh(beta/2) - beta/2))
                ! WK is first order in hbar^2: compare against the exact
                ! expansion to the same order, i.e. the exact ratio itself.
                err(k) = max(err(k), abs(w - exact))
            end do
        end do
        print '(a,2es11.3,a,f6.2)', '     WK harmonic errors ', err, ' ratio ', err(1)/err(2)
        call check(err(1)/err(2) > 12 .and. err(1)/err(2) < 20, &
                   'Wigner-Kirkwood weight = exact harmonic density to O(beta^4)')
        ! Ewald gradient and Laplacian by central differences (the WK inputs).
        call init_classical_cell(cell, 2, 2.0_dp, 1.0_dp, model_coulomb, ok)
        d = [0.31_dp, -0.17_dp, 0.52_dp]*cell%length
        h = 1.0e-3_dp*cell%length
        call cell_pair_energy(cell, d, u0, ok, f)
        lap = 0
        do a = 1, 3
            call cell_pair_energy(cell, d + h*unit(a), up, ok)
            call cell_pair_energy(cell, d - h*unit(a), um, ok)
            fd(a) = -(up - um)/(2*h)
            lap = lap + (up - 2*u0 + um)/h**2
        end do
        print '(a,es11.3,a,2es12.4)', '     |F + grad u|_FD ', maxval(abs(f - fd)), &
            ' lap u FD vs 4 pi g/V ', lap, 4*pi*cell%coupling/cell%volume
        call check(maxval(abs(f - fd)) < 1.0e-5_dp*maxval(abs(f)) .and. &
                   abs(lap - 4*pi*cell%coupling/cell%volume) < &
                   1.0e-3_dp*4*pi*cell%coupling/cell%volume, &
                   'Ewald force and Laplacian (4 pi g/V off contact) vs differences')
    end subroutine

    function unit(a) result(e)
        integer, intent(in) :: a
        real(dp) :: e(3)
        e = 0; e(a) = 1
    end function

    !> Kelbg: Coulomb for r >> lambda, sqrt(pi) g/lambda at contact, series =
    !> closed form at the switch; JVP leaves vs central differences.
    subroutine kelbg_oracle()
        real(dp) :: lam, g, u, us, v, du, up, um, h, dfd, w, dw
        lam = 0.7_dp; g = 1.0_dp
        call qcm_kelbg(8*lam, lam, g, u)
        call check(abs(u - g/(8*lam)) < 1.0e-14_dp*g/lam, 'Kelbg -> g/r for r = 8 lambda')
        call qcm_kelbg_series(0.0_dp, lam, g, us)
        call check(abs(us - sqrt(pi)*g/lam) < 1.0e-15_dp, 'Kelbg contact value sqrt(pi) g/lambda')
        call qcm_kelbg(0.1_dp*lam, lam, g, u)
        call qcm_kelbg_series(0.1_dp*lam, lam, g, us)
        call check(abs(u - us) < 1.0e-13_dp*abs(u), 'Kelbg series = closed form at x = 0.1')
        call qcm_kelbg_excess(0.9_dp, lam, g, v)
        call qcm_kelbg(0.9_dp, lam, g, u)
        call check(abs(v - (u - g/0.9_dp)) < 1.0e-14_dp, 'Kelbg excess = Kelbg - Coulomb')
        h = 1.0e-5_dp
        call qcm_kelbg_jvp(0.9_dp, lam, g, 1.0_dp, 0.0_dp, 0.0_dp, du)
        call qcm_kelbg(0.9_dp + h, lam, g, up); call qcm_kelbg(0.9_dp - h, lam, g, um)
        dfd = (up - um)/(2*h)
        call qcm_pauli_jvp(0.3_dp, 0.8_dp, 0.5_dp, 1.0_dp, 0.0_dp, 0.0_dp, dw)
        call qcm_pauli(0.3_dp + h, 0.8_dp, 0.5_dp, up); call qcm_pauli(0.3_dp - h, 0.8_dp, 0.5_dp, um)
        w = (up - um)/(2*h)
        call check(abs(du - dfd) < 1.0e-8_dp*abs(du) .and. abs(dw - w) < 1.0e-8_dp*abs(w), &
                   'Kelbg/Pauli JVP leaves = central differences')
    end subroutine

    !> Exact thermal-covariance dg/dTheta of the N = 2 quantum cell vs central
    !> differences of the full pipeline.
    subroutine pair_derivative_oracle()
        type(quantum_pair_cell_t) :: pc
        real(dp), allocatable :: edges(:), g(:), gs(:), go(:), dg(:), gp(:), gm(:)
        real(dp) :: f2(26), h, theta, rel
        integer, allocatable :: q(:, :)
        integer :: i, nb
        logical :: ok
        call cell_q_set(q)
        call init_quantum_pair_cell(pc, 2.0_dp, 60, ok)
        nb = 8; theta = 1.0_dp; h = 1.0e-4_dp
        edges = [(i*pc%length/(2*nb), i=0, nb)]
        allocate (g(nb), gs(nb), go(nb), dg(nb), gp(nb), gm(nb))
        call quantum_pair_thermal(pc, theta, edges, q, g, gs, go, f2, ok, dg)
        call quantum_pair_thermal(pc, theta + h, edges, q, gp, gs, go, f2, ok)
        call quantum_pair_thermal(pc, theta - h, edges, q, gm, gs, go, f2, ok)
        rel = maxval(abs(dg - (gp - gm)/(2*h)))/maxval(abs(dg))
        print '(a,es11.3)', '     dg/dTheta spectral vs FD rel ', rel
        call check(ok .and. rel < 1.0e-6_dp, 'N=2 spectral dg/dTheta = central differences')
    end subroutine

    !> Q4 diagonal moments of a single Slater determinant (one occupation
    !> vector) = brute-force Wick sum over the occupied plane-wave orbitals;
    !> dense full-sector ground energy = block-Lanczos ground energy.
    subroutine moments_oracle()
        type(jellium_cell_t) :: jc
        type(jellium_sector_t) :: sector
        type(jellium_spectrum_t) :: one, full, lanczos
        type(fock_hamiltonian_t) :: model
        integer, allocatable :: modes(:, :), q(:, :), triple(:, :)
        real(dp), allocatable :: gamma(:, :), gamma2(:, :), f2(:), f3(:), b2(:), b3(:), &
            w(:)
        integer :: occ(3), k, a, b, c, p, s
        integer, parameter :: perm(3, 6) = reshape([1, 2, 3, 2, 3, 1, 3, 1, 2, &
                                                    2, 1, 3, 1, 3, 2, 3, 2, 1], [3, 6])
        real(dp), parameter :: psign(6) = [1, 1, 1, -1, -1, -1]
        integer :: o(3), ks(3, 1)
        real(dp) :: tr, ground_full, ground_lanczos, trunc, ground, edges(7)
        real(dp), allocatable :: gs1(:), go1(:), gs2(:), go2(:)
        logical :: ok
        call jellium_modes('spherical', 1, modes)
        call init_jellium_cell(jc, 3, modes, 1.0_dp, ok)
        ks(:, 1) = [1, 1, 0]
        ! Three same-spin electrons {0, e_x, e_y}: exchange 3-cycles with all
        ! transfers in S give f3 /= 0.
        call jellium_sector(jc, 3, ks, sector, ok)
        call check(ok, 'N=3 sector K=(110), all up')
        ! The single determinant: first basis state of the sector.
        one%count = 1; one%nsectors = 1
        allocate (one%sectors(1), one%kept(1), one%top(1), one%sector_of(1), &
                  one%offset(1), one%energy(1), one%residual(1), one%flat(sector%d))
        one%sectors(1) = sector; one%kept = 1; one%top = 0; one%sector_of = 1
        one%offset = 0; one%energy = 0; one%residual = 0
        one%flat = 0; one%flat(1) = 1
        allocate (gamma(jc%m, jc%m), gamma2(jc%m*(jc%m - 1)/2, jc%m*(jc%m - 1)/2))
        call jellium_thermal_rdms(jc, one, [1.0_dp], gamma, gamma2, ok)
        call cell_q_set(q); call cell_triple_set(q, triple)
        allocate (f2(size(q, 2)), f3(size(triple, 2)), b2(size(q, 2)), b3(size(triple, 2)))
        allocate (gs1(4), go1(4))
        call jellium_diagonal_moments(jc, one, [1.0_dp], [(a*jc%length/8, a=0, 4)], q, &
                                      triple, gs1, go1, f2, f3, ok)
        call check(ok, 'diagonal moments of a determinant')
        ! Occupied spin orbitals from the 1-RDM diagonal.
        k = 0
        do a = 1, jc%m
            if (gamma(a, a) > 0.5_dp) then
                k = k + 1; occ(k) = a
            end if
        end do
        b2 = 0; b3 = 0
        do s = 1, size(q, 2)
            do a = 1, 3
                do b = 1, 3
                    if (a == b) cycle
                    ! direct minus exchange
                    b2(s) = b2(s) + me(jc, occ(a), occ(a), q(:, s))*me(jc, occ(b), occ(b), -q(:, s)) - &
                            me(jc, occ(a), occ(b), q(:, s))*me(jc, occ(b), occ(a), -q(:, s))
                end do
            end do
        end do
        do s = 1, size(triple, 2)
            do a = 1, 3
                do b = 1, 3
                    do c = 1, 3
                        if (a == b .or. b == c .or. a == c) cycle
                        o = [occ(a), occ(b), occ(c)]
                        do p = 1, 6
                            b3(s) = b3(s) + psign(p)*me(jc, o(1), o(perm(1, p)), q(:, triple(1, s)))* &
                                    me(jc, o(2), o(perm(2, p)), q(:, triple(2, s)))* &
                                    me(jc, o(3), o(perm(3, p)), q(:, triple(3, s)))
                        end do
                    end do
                end do
            end do
        end do
        b2 = b2/6; b3 = b3/6     ! N(N-1) = N(N-1)(N-2) = 6
        tr = maxval(abs(f2 - b2)) + maxval(abs(f3 - b3))
        print '(a,es11.3,a,2es11.3)', '     moments vs Wick ', tr, ' max|f2|,|f3| ', &
            maxval(abs(b2)), maxval(abs(b3))
        call check(tr < 1.0e-13_dp .and. maxval(abs(b3)) > 0.01_dp, &
                   'Q4 diagonal f2/f3 = Slater Wick sum')
        ! Dense full-sector spectrum vs block Lanczos (lowest per sector).
        call jellium_hamiltonian(jc, sector%basis, model, ok)
        call jellium_full_spectrum(jc, model, 1.0e30_dp, full, ok)
        call jellium_build_spectrum(jc, model, 4, lanczos, ok)
        ground_full = minval(full%energy(1:full%count))
        ground_lanczos = minval(lanczos%energy(1:lanczos%count))
        call check(ok .and. abs(ground_full - ground_lanczos) < 1.0e-10_dp*abs(ground_full) &
                   .and. full%count == 364, 'dense full spectrum (364 states) = Lanczos ground')
        ! Sparse-entry shell bins of a thermal ensemble = dense 2-RDM bins.
        allocate (w(full%count), gs2(6), go2(6))
        deallocate (gs1, go1); allocate (gs1(6), go1(6))
        edges = [(a*jc%length/12, a=0, 6)]
        call jellium_thermal_weights(jc, full, 1.0_dp, w, trunc, ground)
        call jellium_thermal_rdms(jc, full, w, gamma, gamma2, ok)
        call jellium_pair_bins(jc, gamma2, edges, gs2, go2)
        call jellium_diagonal_moments(jc, full, w, edges, q, triple, gs1, go1, f2, f3, ok)
        tr = maxval(abs(gs1 - gs2)) + maxval(abs(go1 - go2))
        print '(a,es11.3,a,es10.2)', '     sparse vs dense bins ', tr, ' trunc ', trunc
        call check(ok .and. tr < 1.0e-12_dp .and. trunc == 0, &
                   'sparse 2-RDM shell bins = dense jellium_pair_bins (Theta = 1, N = 3)')
    end subroutine

    real(dp) function me(jc, x, y, qq) result(v)
        type(jellium_cell_t), intent(in) :: jc
        integer, intent(in) :: x, y, qq(3)
        v = 0
        if (jc%spin(x) == jc%spin(y) .and. all(jc%momentum(:, x) - &
                                               jc%momentum(:, y) == qq)) v = 1
    end function
end program
