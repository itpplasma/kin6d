!> H1b classical parent: N point electrons (N = 2..4) in the Q4 periodic cube
!> with neutralizing background, canonical at (n, T). Spin is irrelevant to
!> the classical Boltzmann weight; spin labels are spectators.
!>
!> The periodic energy is pairwise additive, U = const + sum_{i<j} u(x_i-x_j),
!> with u the zero-mean Ewald pair potential of the shared kin6d_nbody/kin6d_ewald path
!> (a two-particle periodic system). Optional semiclassical pair models add
!> the periodic-image sum of the Kelbg diffraction excess and the
!> spin-averaged Pauli potential (kin6d_quantum_classical_leaves).
!>
!> Estimators (common with the quantum cell, docs/quantum-classical-
!> correlations.md): shell-averaged spin-summed g on equal shells in |r|<L/2,
!> diagonal pair/triplet Fourier moments f2(Q), f3(Q1,Q2) on the set
!> S = {1 <= |n|^2 <= 3}, Wigner--Kirkwood and leading-exchange reweighting,
!> and the beta derivative of g by the energy covariance. Monte Carlo
!> intervals are batch-means statistical intervals (labelled), not bounds.
module kin6d_classical_cell
    use, intrinsic :: iso_fortran_env, only: dp => real64, i64 => int64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortnum_status, only: fortnum_status_t, FORTNUM_OK
    use fortnum_quadrature, only: gauss_legendre, gauss_legendre_ab
    use kin6d_species, only: units_t, species_t, plasma_scales_t, atomic_units, &
        electron, plasma_scales, pair_coupling
    use kin6d_ewald, only: ewald_parameters_t, ewald_select
    use kin6d_coulomb_leaves, only: ewald_background
    use kin6d_nbody, only: nbody_system_t, nbody_periodic, nbody_force, nbody_wrap
    use kin6d_quantum_classical_leaves, only: qcm_pair_length, qcm_kelbg, &
        qcm_kelbg_series, qcm_kelbg_excess, qcm_exchange, qcm_pauli, &
        qcm_wigner_kirkwood
    implicit none
    private
    public :: classical_cell_t, classical_mc_t, init_classical_cell
    public :: cell_pair_energy, cell_exchange_sum, cell_total_energy
    public :: cell_q_set, cell_triple_set, classical_cell_mc
    public :: classical_pair_quadrature, cell_shell_volume
    integer, parameter, public :: model_coulomb = 1, model_kelbg = 2, &
        model_kelbg_pauli = 3
    real(dp), parameter :: pi = acos(-1.0_dp)

    type :: classical_cell_t
        integer :: n = 0, model = model_coulomb, images = 0
        real(dp) :: length = 0, volume = 0, coupling = 1, mass = 1, hbar = 1
        real(dp) :: rs = 0, theta = 0, gamma = 0, fermi = 0, kt = 0, beta = 0
        !> Pair thermal length hbar/sqrt(2 mu kT), mu = m/2 (Kelbg/exchange).
        real(dp) :: lambda = 0
        !> Same-spin probability of a pair for spin-blind labels (1/2).
        real(dp) :: same_spin = 0.5_dp
        real(dp) :: pair_constant = 0
        type(nbody_system_t) :: pair, full
    end type

    !> Monte Carlo result; *_err are batch-means standard errors.
    type :: classical_mc_t
        real(dp), allocatable :: g(:), g_err(:), g_wk(:), g_x(:), g_pred(:)
        real(dp), allocatable :: g_pred_err(:), dg_dbeta(:), dg_dbeta_err(:)
        real(dp), allocatable :: f2(:), f2_err(:), f3(:), f3_err(:)
        real(dp) :: c2 = 0, c2_err = 0, c3 = 0, c3_err = 0
        real(dp) :: acceptance = 0, energy = 0, energy_err = 0, step = 0
        integer :: samples = 0
    end type
contains
    !> Cell of n electrons at density parameter rs (L from the shared rs map,
    !> as kin6d_jellium_cell) and Theta = kT/E_F; model selects the pair law.
    subroutine init_classical_cell(cell, n, rs, theta, model, ok, tolerance)
        type(classical_cell_t), intent(out) :: cell
        integer, intent(in) :: n, model
        real(dp), intent(in) :: rs, theta
        logical, intent(out) :: ok
        real(dp), intent(in), optional :: tolerance
        type(units_t) :: units
        type(species_t) :: e
        type(species_t), allocatable :: all(:)
        type(plasma_scales_t) :: unit_cell, scales
        type(ewald_parameters_t) :: p
        real(dp) :: tol, bg
        ok = .false.
        if (n < 2 .or. .not. (rs > 0) .or. .not. (theta > 0)) return
        if (model < model_coulomb .or. model > model_kelbg_pauli) return
        tol = 1.0e-10_dp
        if (present(tolerance)) tol = tolerance
        units = atomic_units(); e = electron()
        unit_cell = plasma_scales(units, e, real(n, dp), 1.0_dp)
        cell%length = rs/unit_cell%rs; cell%volume = cell%length**3
        scales = plasma_scales(units, e, real(n, dp)/cell%volume, 1.0_dp)
        cell%fermi = scales%fermi; cell%kt = theta*scales%fermi
        scales = plasma_scales(units, e, real(n, dp)/cell%volume, cell%kt)
        cell%rs = scales%rs; cell%theta = scales%theta; cell%gamma = scales%gamma
        cell%beta = 1/cell%kt; cell%n = n; cell%model = model
        cell%mass = e%mass; cell%hbar = units%hbar
        cell%coupling = pair_coupling(units, e, e)
        call qcm_pair_length(cell%hbar, cell%mass/2, cell%kt, cell%lambda)
        ! Images |d + nL| >= (images - 1/2) L carry Kelbg/exchange terms below
        ! exp(-36) relative to their contact values.
        cell%images = ceiling((6*cell%lambda + cell%length/2)/cell%length)
        call ewald_select(cell%length, tol, p, ok, 2)
        if (.not. ok) return
        call nbody_periodic(cell%pair, units, [e, e], p, ok)
        if (.not. ok) return
        ! Zero-mean pair potential (Fourier Q = 0 removed, as in the quantum
        ! cell): remove the two reciprocal self terms (inside |S(k)|^2) and add
        ! back the pair's share g of the background constant.
        call ewald_background(p%alpha, cell%volume, cell%coupling, bg)
        cell%pair_constant = cell%pair%constant_energy - bg + &
            cell%coupling*sum(cell%pair%weight(1:cell%pair%modes))
        allocate (all(n)); all = e
        call ewald_select(cell%length, tol, p, ok, n)
        if (.not. ok) return
        call nbody_periodic(cell%full, units, all, p, ok)
    end subroutine

    !> Pair energy u(d) of the cell model and the force -grad u on the first
    !> particle (Coulomb part only; the Kelbg/Pauli terms enter u only).
    subroutine cell_pair_energy(cell, d, u, ok, force)
        type(classical_cell_t), intent(inout) :: cell
        real(dp), intent(in) :: d(3)
        real(dp), intent(out) :: u
        logical, intent(out) :: ok
        real(dp), intent(out), optional :: force(3)
        type(fortnum_status_t) :: status
        real(dp) :: x(6), f(6), r, e(3), v, xs
        integer :: a, b, c
        x = 0; x(1:3) = nbody_wrap(d, cell%length)
        call nbody_force(cell%pair, x, f, u, status)
        ok = status%code == FORTNUM_OK
        if (.not. ok) return
        u = u - cell%pair_constant
        if (present(force)) force = f(1:3)
        if (cell%model == model_coulomb) return
        r = norm2(x(1:3))
        ! Replace the bare minimum-image g/r by the full Kelbg value (finite
        ! at contact) and add the Kelbg excess of every other image.
        if (r < 0.1_dp*cell%lambda) then
            call qcm_kelbg_series(r, cell%lambda, cell%coupling, v)
        else
            call qcm_kelbg(r, cell%lambda, cell%coupling, v)
        end if
        u = u - cell%coupling/r + v
        xs = 0
        do a = -cell%images, cell%images
            do b = -cell%images, cell%images
                do c = -cell%images, cell%images
                    e = x(1:3) + cell%length*[a, b, c]
                    r = norm2(e)
                    if (a /= 0 .or. b /= 0 .or. c /= 0) then
                        call qcm_kelbg_excess(r, cell%lambda, cell%coupling, v)
                        u = u + v
                    end if
                    call qcm_exchange(r, cell%lambda, v)
                    xs = xs + v
                end do
            end do
        end do
        if (cell%model == model_kelbg_pauli) then
            call qcm_pauli(xs, cell%kt, cell%same_spin, v)
            u = u + v
        end if
        ok = ieee_is_finite(u)
    end subroutine

    !> Periodic leading-order exchange factor sum_n exp(-|d + nL|^2/lambda^2).
    real(dp) function cell_exchange_sum(cell, d) result(xs)
        type(classical_cell_t), intent(in) :: cell
        real(dp), intent(in) :: d(3)
        real(dp) :: e(3), v
        integer :: a, b, c
        xs = 0
        do a = -cell%images, cell%images
            do b = -cell%images, cell%images
                do c = -cell%images, cell%images
                    e = nbody_wrap(d, cell%length) + cell%length*[a, b, c]
                    call qcm_exchange(norm2(e), cell%lambda, v)
                    xs = xs + v
                end do
            end do
        end do
    end function

    !> Configuration energy sum_{i<j} u(x_i - x_j) (constant omitted).
    subroutine cell_total_energy(cell, x, energy, ok)
        type(classical_cell_t), intent(inout) :: cell
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: energy
        logical, intent(out) :: ok
        real(dp) :: u
        integer :: i, j
        energy = 0; ok = .true.
        do i = 1, cell%n - 1
            do j = i + 1, cell%n
                call cell_pair_energy(cell, x(3*i - 2:3*i) - x(3*j - 2:3*j), u, ok)
                if (.not. ok) return
                energy = energy + u
            end do
        end do
    end subroutine

    !> Common Fourier set S = {n in Z^3: 1 <= |n|^2 <= 3}, Q = 2 pi n/L.
    subroutine cell_q_set(q)
        integer, allocatable, intent(out) :: q(:, :)
        integer :: a, b, c, m, buffer(3, 27)
        m = 0
        do a = -1, 1
            do b = -1, 1
                do c = -1, 1
                    if (a*a + b*b + c*c == 0) cycle
                    m = m + 1; buffer(:, m) = [a, b, c]
                end do
            end do
        end do
        q = buffer(:, 1:m)
    end subroutine

    !> Ordered pairs (a, b) of S with Q3 = -(Q_a + Q_b) also in S;
    !> triple(3, :) is the index of Q3.
    subroutine cell_triple_set(q, triple)
        integer, intent(in) :: q(:, :)
        integer, allocatable, intent(out) :: triple(:, :)
        integer :: a, b, c, m, buffer(3, 27*27)
        m = 0
        do a = 1, size(q, 2)
            do b = 1, size(q, 2)
                do c = 1, size(q, 2)
                    if (all(q(:, a) + q(:, b) + q(:, c) == 0)) then
                        m = m + 1; buffer(:, m) = [a, b, c]
                    end if
                end do
            end do
        end do
        triple = buffer(:, 1:m)
    end subroutine

    real(dp) function cell_shell_volume(ra, rb) result(v)
        real(dp), intent(in) :: ra, rb
        v = 4*pi*(rb**3 - ra**3)/3
    end function

    !> Canonical Metropolis Monte Carlo of the cell (single-particle uniform
    !> cube moves, step adapted to acceptance 1/2 during equilibration only).
    !> One measurement per sweep of n moves; batch means over `batches`.
    !> reweight (Coulomb model) adds the Wigner--Kirkwood and leading-order
    !> exchange reweighted g and the covariance derivative dg/dbeta.
    subroutine classical_cell_mc(cell, edges, sweeps, equilibration, batches, &
                                 seed, mc, ok, reweight, sample_seconds, summary_seconds)
        type(classical_cell_t), intent(inout) :: cell
        real(dp), intent(in) :: edges(:)
        integer, intent(in) :: sweeps, equilibration, batches
        integer(i64), intent(in) :: seed
        type(classical_mc_t), intent(out) :: mc
        logical, intent(out) :: ok
        logical, intent(in), optional :: reweight
        real(dp), intent(out), optional :: sample_seconds, summary_seconds
        integer, allocatable :: q(:, :), triple(:, :)
        real(dp), allocatable :: x(:), pe(:, :), bg(:, :), bwk(:, :), bx(:, :), &
            bp(:, :), bd(:, :), bf2(:, :), bf3(:, :), be(:, :), cnt(:), f(:)
        real(dp), allocatable :: s_c(:), s_cwk(:), s_cx(:), s_cp(:), s_ce(:)
        real(dp), allocatable :: s_f2(:), s_f3(:), tmp(:), tmp_err(:)
        complex(dp), allocatable :: rho(:)
        real(dp) :: clock_start, clock_summary, clock_end
        real(dp) :: step, trial(3), du, u, unew(64), e_total, s_wk, s_x, s_p, s_e
        real(dp) :: w_wk, w_x, lap, grad2, pot, norm, r, accepted, attempts
        integer(i64) :: state(2)
        integer :: nb, nq, nt, i, j, k, s, b, per, sweep, n, bin
        type(fortnum_status_t) :: status
        logical :: rw, move_ok

        ok = .false.
        if (present(sample_seconds)) sample_seconds = 0
        if (present(summary_seconds)) summary_seconds = 0
        call cpu_time(clock_start)
        rw = .false.
        if (present(reweight)) rw = reweight .and. cell%model == model_coulomb
        n = cell%n; nb = size(edges) - 1
        if (nb < 1 .or. batches < 2 .or. sweeps < batches .or. n > 64) return
        if (edges(nb + 1) > cell%length/2*(1 + 1.0e-12_dp)) return
        call cell_q_set(q); call cell_triple_set(q, triple)
        nq = size(q, 2); nt = size(triple, 2)
        allocate (x(3*n), pe(n, n), f(3*n), rho(nq), cnt(nb))
        allocate (bg(nb, batches), bwk(nb, batches), bx(nb, batches), &
                  bp(nb, batches), bd(nb, batches), bf2(nq, batches), &
                  bf3(nt, batches), be(1, batches))
        allocate (s_c(nb), s_cwk(nb), s_cx(nb), s_cp(nb), s_ce(nb), s_f2(nq), &
                  s_f3(nt))
        state = [1 + modulo(seed, 2147483562_i64), 1 + modulo(7*seed + 12345, 2147483398_i64)]
        do i = 1, 3*n
            x(i) = (uniform(state) - 0.5_dp)*cell%length
        end do
        pe = 0
        do i = 1, n - 1
            do j = i + 1, n
                call cell_pair_energy(cell, x(3*i - 2:3*i) - x(3*j - 2:3*j), u, ok)
                if (.not. ok) return
                pe(i, j) = u; pe(j, i) = u
            end do
        end do
        ok = .false.
        step = 0.3_dp*cell%length
        accepted = 0; attempts = 0
        per = sweeps/batches
        bg = 0; bwk = 0; bx = 0; bp = 0; bd = 0; bf2 = 0; bf3 = 0; be = 0
        do sweep = 1, equilibration + per*batches
            do k = 1, n
                i = min(n, 1 + int(uniform(state)*n))
                do j = 1, 3
                    trial(j) = x(3*i - 3 + j) + step*(2*uniform(state) - 1)
                end do
                trial = nbody_wrap(trial, cell%length)
                du = 0; move_ok = .true.
                do j = 1, n
                    if (j == i) cycle
                    call cell_pair_energy(cell, trial - x(3*j - 2:3*j), unew(j), move_ok)
                    if (.not. move_ok) exit
                    du = du + unew(j) - pe(i, j)
                end do
                attempts = attempts + 1
                if (.not. move_ok) cycle
                if (cell%beta*du <= 0 .or. uniform(state) < exp(-cell%beta*du)) then
                    accepted = accepted + 1
                    x(3*i - 2:3*i) = trial
                    do j = 1, n
                        if (j == i) cycle
                        pe(i, j) = unew(j); pe(j, i) = unew(j)
                    end do
                end if
            end do
            if (sweep <= equilibration) then
                if (mod(sweep, 50) == 0) then
                    step = min(cell%length/2, step*max(0.5_dp, min(2.0_dp, &
                                                                   2*accepted/attempts)))
                    accepted = 0; attempts = 0
                end if
                cycle
            end if
            if (sweep == equilibration + 1) then
                accepted = 0; attempts = 0
            end if
            s = sweep - equilibration; b = (s - 1)/per + 1
            if (mod(s - 1, per) == 0) then
                s_c = 0; s_cwk = 0; s_cx = 0; s_cp = 0; s_ce = 0; s_f2 = 0; s_f3 = 0
                s_wk = 0; s_x = 0; s_p = 0; s_e = 0
            end if
            cnt = 0; w_x = 1
            do i = 1, n - 1
                do j = i + 1, n
                    trial = nbody_wrap(x(3*i - 2:3*i) - x(3*j - 2:3*j), cell%length)
                    r = norm2(trial)
                    if (rw) w_x = w_x - cell%same_spin*cell_exchange_sum(cell, trial)
                    if (r >= edges(nb + 1) .or. r < edges(1)) cycle
                    do bin = 1, nb - 1
                        if (r < edges(bin + 1)) exit
                    end do
                    cnt(bin) = cnt(bin) + 1
                end do
            end do
            e_total = 0
            do i = 1, n - 1
                e_total = e_total + sum(pe(i, i + 1:n))
            end do
            w_wk = 1
            if (rw) then
                call nbody_force(cell%full, x, f, pot, status)
                if (status%code /= FORTNUM_OK) return
                grad2 = sum(f**2)
                lap = n*(n - 1)*4*pi*cell%coupling/cell%volume
                call qcm_wigner_kirkwood(cell%beta, cell%hbar, cell%mass, lap, &
                                         grad2, w_wk)
            end if
            do k = 1, nq
                rho(k) = 0
                do i = 1, n
                    rho(k) = rho(k) + exp(cmplx(0.0_dp, 2*pi*dot_product( &
                                   real(q(:, k), dp), x(3*i - 2:3*i))/cell%length, dp))
                end do
                s_f2(k) = s_f2(k) + abs(rho(k))**2 - n
            end do
            do k = 1, nt
                s_f3(k) = s_f3(k) + real(rho(triple(1, k))*rho(triple(2, k))* &
                          rho(triple(3, k))) - abs(rho(triple(1, k)))**2 - &
                          abs(rho(triple(2, k)))**2 - abs(rho(triple(3, k)))**2 + 2*n
            end do
            s_c = s_c + cnt; s_cwk = s_cwk + w_wk*cnt; s_cx = s_cx + w_x*cnt
            s_cp = s_cp + w_wk*w_x*cnt; s_ce = s_ce + e_total*cnt
            s_wk = s_wk + w_wk; s_x = s_x + w_x; s_p = s_p + w_wk*w_x
            s_e = s_e + e_total
            if (mod(s, per) == 0) then
                do bin = 1, nb
                    norm = per*n*(n - 1)/2.0_dp* &
                           cell_shell_volume(edges(bin), edges(bin + 1))/cell%volume
                    bg(bin, b) = s_c(bin)/norm
                    bwk(bin, b) = s_cwk(bin)/norm*per/s_wk
                    bx(bin, b) = s_cx(bin)/norm*per/s_x
                    bp(bin, b) = s_cp(bin)/norm*per/s_p
                    bd(bin, b) = -(s_ce(bin)/norm - bg(bin, b)*s_e/per)
                end do
                bf2(:, b) = s_f2/(per*real(n*(n - 1), dp))
                if (n >= 3) bf3(:, b) = s_f3/(per*real(n*(n - 1)*(n - 2), dp))
                be(1, b) = s_e/per
            end if
        end do
        call cpu_time(clock_summary)
        mc%acceptance = accepted/max(1.0_dp, attempts); mc%step = step
        mc%samples = per*batches
        call batch_stats(bg, mc%g, mc%g_err)
        call batch_stats(bwk, mc%g_wk, tmp_err)
        call batch_stats(bx, mc%g_x, tmp_err)
        call batch_stats(bp, mc%g_pred, mc%g_pred_err)
        call batch_stats(bd, mc%dg_dbeta, mc%dg_dbeta_err)
        call batch_stats(bf2, mc%f2, mc%f2_err)
        call batch_stats(bf3, mc%f3, mc%f3_err)
        call batch_stats(be, tmp, tmp_err)
        mc%energy = tmp(1); mc%energy_err = tmp_err(1)
        call jackknife_norm(bf2, mc%c2, mc%c2_err)
        mc%c3 = 0; mc%c3_err = 0
        if (n >= 3) call jackknife_norm(bf3, mc%c3, mc%c3_err)
        call cpu_time(clock_end)
        if (present(sample_seconds)) sample_seconds = clock_summary - clock_start
        if (present(summary_seconds)) summary_seconds = clock_end - clock_summary
        ok = all(ieee_is_finite(mc%g)) .and. ieee_is_finite(mc%c2)
    end subroutine

    !> Mean and standard error over batch columns.
    subroutine batch_stats(a, mean, err)
        real(dp), intent(in) :: a(:, :)
        real(dp), allocatable, intent(out) :: mean(:), err(:)
        integer :: nb
        nb = size(a, 2)
        mean = sum(a, 2)/nb
        err = sqrt(max(0.0_dp, sum((a - spread(mean, 2, nb))**2, 2))/ &
                   (nb*(nb - 1.0_dp)))
    end subroutine

    !> Root-mean-square norm sqrt(mean_k abar_k^2) with the jackknife bias
    !> correction (removes the O(1/batches) noise-squared bias) and error.
    subroutine jackknife_norm(a, c, err)
        real(dp), intent(in) :: a(:, :)
        real(dp), intent(out) :: c, err
        real(dp), allocatable :: total(:), leave(:)
        integer :: nb, j
        nb = size(a, 2); c = 0; err = 0
        if (size(a, 1) == 0) return
        total = sum(a, 2)
        allocate (leave(nb))
        do j = 1, nb
            leave(j) = sqrt(sum(((total - a(:, j))/(nb - 1))**2)/size(a, 1))
        end do
        c = sqrt(sum((total/nb)**2)/size(a, 1))
        c = nb*c - (nb - 1)*sum(leave)/nb
        err = sqrt((nb - 1.0_dp)/nb*sum((leave - sum(leave)/nb)**2))
    end subroutine

    !> L'Ecuyer (1988) combined multiplicative generator, period ~2.3e18,
    !> exact in 64-bit integers; deviate in (0, 1).
    real(dp) function uniform(state) result(u)
        integer(i64), intent(inout) :: state(2)
        integer(i64) :: z
        state(1) = mod(40014_i64*state(1), 2147483563_i64)
        state(2) = mod(40692_i64*state(2), 2147483399_i64)
        z = state(1) - state(2)
        if (z < 1) z = z + 2147483562_i64
        u = real(z, dp)*4.656613057391769e-10_dp
    end function

    !> Deterministic N = 2 reference: the relative coordinate of two particles
    !> in the cell has density exp(-beta u(d))/Z on the torus (the centre of
    !> mass integrates out). Z and f2 by the periodic midpoint rule on an
    !> ngrid^3 torus grid (spectrally accurate for the smooth periodic
    !> integrand); shell averages by Gauss--Legendre in r and cos(theta) and the
    !> midpoint rule in phi with `nshell` points per direction scale. Also the
    !> Wigner--Kirkwood, leading-exchange and combined predictions and, for the
    !> Coulomb model, dg/dbeta. Outputs are on the common estimator contract.
    subroutine classical_pair_quadrature(cell, edges, ngrid, nshell, g, f2, g_wk, &
                                         g_x, g_pred, dg_dbeta, ok)
        type(classical_cell_t), intent(inout) :: cell
        real(dp), intent(in) :: edges(:)
        integer, intent(in) :: ngrid, nshell
        real(dp), intent(out) :: g(:), f2(:), g_wk(:), g_x(:), g_pred(:), dg_dbeta(:)
        logical, intent(out) :: ok
        integer, allocatable :: q(:, :)
        real(dp), allocatable :: u(:, :, :), w(:, :, :), wx(:, :, :), xr(:), wr(:), &
            xt(:), wt(:)
        real(dp) :: h, d(3), f(3), umin, b, z, zw, zx, zp, zu, lap, val, s(6), &
            st, ph, vs
        integer :: i, j, k, m, nb, a, c, e
        ok = .false.
        nb = size(edges) - 1
        if (cell%n /= 2 .or. ngrid < 4 .or. nshell < 2) return
        if (size(g) /= nb .or. size(f2) /= 26) return
        call cell_q_set(q)
        h = cell%length/ngrid
        allocate (u(ngrid, ngrid, ngrid), w(ngrid, ngrid, ngrid), wx(ngrid, ngrid, ngrid))
        lap = 2*4*pi*cell%coupling/cell%volume
        do k = 1, ngrid
            do j = 1, ngrid
                do i = 1, ngrid
                    d = -cell%length/2 + h*([i, j, k] - 0.5_dp)
                    call point(d, u(i, j, k), w(i, j, k), wx(i, j, k))
                    if (.not. ok) return
                end do
            end do
        end do
        umin = minval(u)
        z = 0; zw = 0; zx = 0; zp = 0; zu = 0; f2 = 0
        do k = 1, ngrid
            do j = 1, ngrid
                do i = 1, ngrid
                    d = -cell%length/2 + h*([i, j, k] - 0.5_dp)
                    b = exp(-cell%beta*(u(i, j, k) - umin))
                    z = z + b; zw = zw + b*w(i, j, k); zx = zx + b*wx(i, j, k)
                    zp = zp + b*w(i, j, k)*wx(i, j, k); zu = zu + b*u(i, j, k)
                    do m = 1, size(q, 2)
                        f2(m) = f2(m) + b*cos(2*pi*dot_product(real(q(:, m), dp), d)/ &
                                              cell%length)
                    end do
                end do
            end do
        end do
        f2 = f2/z
        ! Shell integrals: radial and polar Gauss--Legendre, midpoint in phi.
        allocate (xr(nshell), wr(nshell), xt(2*nshell), wt(2*nshell))
        call gauss_legendre(2*nshell, xt, wt)
        do m = 1, nb
            call gauss_legendre_ab(nshell, edges(m), edges(m + 1), xr, wr)
            s = 0
            do a = 1, nshell
                do c = 1, 2*nshell
                    st = sqrt(max(0.0_dp, 1 - xt(c)**2))
                    do e = 1, 4*nshell
                        ph = 2*pi*(e - 0.5_dp)/(4*nshell)
                        d = xr(a)*[st*cos(ph), st*sin(ph), xt(c)]
                        call point(d, val, f(1), f(2))
                        if (.not. ok) return
                        b = exp(-cell%beta*(val - umin))*wr(a)*xr(a)**2*wt(c)* &
                            2*pi/(4*nshell)
                        s(1) = s(1) + b; s(2) = s(2) + b*f(1); s(3) = s(3) + b*f(2)
                        s(4) = s(4) + b*f(1)*f(2); s(5) = s(5) + b*val
                    end do
                end do
            end do
            vs = cell_shell_volume(edges(m), edges(m + 1))
            ! Grid sums carry h^3: normalized densities V p = V s/(Z h^3).
            g(m) = cell%volume*s(1)/vs/(z*h**3)
            g_wk(m) = cell%volume*s(2)/vs/(zw*h**3)
            g_x(m) = cell%volume*s(3)/vs/(zx*h**3)
            g_pred(m) = cell%volume*s(4)/vs/(zp*h**3)
            dg_dbeta(m) = 0
            if (cell%model == model_coulomb) dg_dbeta(m) = &
                -(cell%volume*s(5)/vs/(z*h**3) - g(m)*zu/z)
        end do
        ok = all(ieee_is_finite(g)) .and. all(ieee_is_finite(g_pred))
    contains
        subroutine point(dd, uu, ww, xx)
            real(dp), intent(in) :: dd(3)
            real(dp), intent(out) :: uu, ww, xx
            real(dp) :: ff(3)
            call cell_pair_energy(cell, dd, uu, ok, ff)
            if (.not. ok) return
            ! Two particles: sum_i |grad_i U|^2 = 2 |grad u|^2.
            call qcm_wigner_kirkwood(cell%beta, cell%hbar, cell%mass, lap, &
                                     2*sum(ff**2), ww)
            xx = 1 - cell%same_spin*cell_exchange_sum(cell, dd)
        end subroutine
    end subroutine
end module kin6d_classical_cell
