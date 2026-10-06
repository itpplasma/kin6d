program test_correlation_hierarchy
    ! Independent oracles of the H1a correlation-hierarchy map:
    ! - MC Ewald energy against the tiled direct Ewald sum of kin6d_nbody with
    !   different (tighter) parameters; incremental energy against a full
    !   recompute; MC forces against central differences of the MC energy;
    ! - the Metropolis chain for N = 2 against <U> by direct quadrature of
    !   the periodic pair Boltzmann factor (kin6d_nbody energy);
    ! - triplet bin measure against the analytic two-ball overlap integral;
    !   ideal-gas g3 = 1 and c3 = 0;
    ! - YBG first-order closure force against the analytic Gaussian-cloud
    !   potential;
    ! - HNC against the Caillol-Gilles (2010) HNC energies and the Debye
    !   limit; DST against a direct sine sum;
    ! - autocorrelation time of an AR(1) series; the fit-derivative leaf
    !   against finite differences.
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use fortnum_status, only: fortnum_status_t, FORTNUM_OK
    use fortnum_rng, only: rng_t, rng_seed, rng_uniform, rng_normal
    use kin6d_species, only: atomic_units, electron, species_t
    use kin6d_ewald, only: ewald_parameters_t, ewald_select
    use kin6d_nbody, only: nbody_system_t, nbody_periodic, nbody_energy
    use kin6d_ocp_monte_carlo
    use kin6d_correlation_hierarchy
    use kin6d_ocp_hnc
    use kin6d_series_statistics
    use kin6d_correlation_hierarchy_leaves, only: ch_ocp_energy_fit, &
        ch_ocp_energy_fit_jvp, ch_debye_hueckel_energy, ch_debye_wavenumber
    implicit none
    real(dp), parameter :: pi = 3.14159265358979323846_dp

    call monte_carlo_energy()
    call two_particle_chain()
    call triplet_ideal()
    call ybg_gaussian()
    call hnc_oracles()
    call statistics_and_fit()
    print '(a)', 'PASS: MC Ewald/forces/N=2 chain, g3 measure/ideal gas, YBG, HNC, stats'
contains
    subroutine check(condition, message)
        logical, intent(in) :: condition
        character(*), intent(in) :: message
        if (.not. condition) then
            print '(a)', 'FAIL: '//message
            error stop 1
        end if
    end subroutine

    !> kin6d_nbody direct Ewald system of n electrons (tight, real_cutoff 2).
    subroutine direct_system(n, length, system)
        integer, intent(in) :: n
        real(dp), intent(in) :: length
        type(nbody_system_t), intent(out) :: system
        type(ewald_parameters_t) :: p
        type(species_t) :: particles(n)
        logical :: ok
        call ewald_select(length, 1.0e-13_dp, p, ok, particles=n, real_cutoff=2)
        call check(ok, 'direct Ewald parameters')
        particles = electron()
        call nbody_periodic(system, atomic_units(), particles, p, ok)
        call check(ok, 'nbody periodic setup')
    end subroutine

    !> Reduced energy of positions x(3, n) by the direct Ewald sum.
    real(dp) function direct_energy(system, x, n) result(u)
        type(nbody_system_t), intent(inout) :: system
        real(dp), intent(in) :: x(:, :)
        integer, intent(in) :: n
        type(fortnum_status_t) :: status
        real(dp) :: state(6*n)
        state = 0; state(1:3*n) = reshape(x(:, 1:n), [3*n])
        call nbody_energy(system, state, u, status)
        call check(status%code == FORTNUM_OK, 'nbody energy')
    end function

    subroutine monte_carlo_energy()
        type(ocp_mc_t) :: mc
        type(nbody_system_t) :: system
        logical :: ok
        integer :: i, j, a
        real(dp) :: u_direct, u_full, f(3, 64), fd, h, up, um, worst, scale
        call mc_setup(mc, 64, 2.0_dp, 1.0e-8_dp, 11_int64, ok)
        call check(ok, 'MC setup')
        call check(mc%ewald%real_cutoff == 0, 'minimum-image real space')
        call direct_system(64, mc%length, system)
        u_direct = direct_energy(system, mc%x, 64)
        ! Truncation bound per pair in 1/L, N(N+1)/2 pair/self terms.
        print '(a,es10.3,a,es10.3)', 'MC vs direct Ewald energy: ', &
            abs(mc%energy - u_direct), ' bound ', &
            mc%ewald%potential_bound*64*65/2/mc%length
        call check(abs(mc%energy - u_direct) <= &
                   mc%ewald%potential_bound*64*65/2/mc%length, 'MC Ewald energy')
        do i = 1, 30
            call mc_sweep(mc, i <= 15)
        end do
        u_full = mc_energy(mc)
        call check(abs(u_full - mc%energy) <= 1.0e-10_dp*abs(u_full), &
                   'incremental energy drift')
        call check(abs(u_full - direct_energy(system, mc%x, 64)) <= &
                   mc%ewald%potential_bound*64*65/2/mc%length, 'MC energy after sweeps')
        call mc_forces(mc, f)
        h = 1.0e-5_dp; worst = 0; scale = maxval(abs(f))
        do i = 1, 64, 9
            do a = 1, 3
                mc%x(a, i) = mc%x(a, i) + h; up = mc_energy(mc)
                mc%x(a, i) = mc%x(a, i) - 2*h; um = mc_energy(mc)
                mc%x(a, i) = mc%x(a, i) + h
                fd = -(up - um)/(2*h)
                worst = max(worst, abs(fd - f(a, i))/scale)
            end do
        end do
        j = 0
        print '(a,es10.3)', 'MC force vs central difference (rel. max): ', worst
        call check(worst < 1.0e-7_dp, 'MC forces')
        call check(abs(sum(f(1, :))) + abs(sum(f(2, :))) + abs(sum(f(3, :))) < &
                   1.0e-10_dp*scale*64, 'momentum conservation of forces')
        u_full = mc_energy(mc)
        call check(j == 0 .and. abs(u_full - mc%energy) <= 1.0e-10_dp*abs(u_full), &
                   'energy restored')
    end subroutine

    !> N = 2: the relative minimum-image vector d is distributed with
    !> exp(-Gamma phi(d)) on the cube; <U> by a midpoint grid of d.
    subroutine two_particle_chain()
        type(ocp_mc_t) :: mc
        type(nbody_system_t) :: system
        logical :: ok
        integer, parameter :: m = 32, sweeps = 100000
        real(dp) :: x(3, 2), length, w, z, uz, phi, mean, err, tau, ess
        real(dp), allocatable :: series(:)
        integer :: i, j, k, s
        real(dp), parameter :: gamma = 2.0_dp
        call mc_setup(mc, 2, gamma, 1.0e-6_dp, 5_int64, ok)
        call check(ok, 'N=2 setup')
        length = mc%length
        call direct_system(2, length, system)
        z = 0; uz = 0
        do k = 1, m
            do j = 1, m
                do i = 1, m
                    x(:, 1) = 0
                    x(:, 2) = length*([i, j, k] - 0.5_dp)/m
                    phi = direct_energy(system, x, 2)
                    w = exp(-gamma*phi)
                    z = z + w; uz = uz + w*phi
                end do
            end do
        end do
        allocate (series(sweeps))
        mc%step = length/2
        do s = 1, sweeps
            call mc_sweep(mc, .false.)
            series(s) = mc%energy
        end do
        call series_autocorrelation(series, mean, err, tau, ess)
        print '(a,f10.6,a,f10.6,a,f9.6,a,f6.2)', 'N=2 <U>: MC ', mean, ' quadrature ', &
            uz/z, ' stat. err ', err, ' tau ', tau
        ! Midpoint grid of d: m = 24 and m = 32 agree to 1e-6.
        call check(abs(mean - uz/z) < 4*err + 2.0e-4_dp, 'N=2 Metropolis <U>')
    end subroutine

    !> Triplet measure: sum over bins = int_{|r2|,|r3|,|r2-r3| < R} d3r2 d3r3
    !> = int_0^R 4 pi d^2 V(d) dd, V(d) = pi (4R + d)(2R - d)^2/12 (two-ball
    !> overlap at distance d); ideal gas g3 = 1.
    subroutine triplet_ideal()
        type(triplet_histogram_t) :: t, fine, ref
        real(dp) :: worst
        type(rng_t) :: g
        type(fortnum_status_t) :: status
        real(dp) :: exact, d, rmax, x(3, 300), length, g3(8, 8, 8), ratio, sigma
        real(dp) :: per(200), prev
        real(dp) :: rg(40), one(40), kirk(8, 8, 8), lin(8, 8, 8), cen(8, 8, 8)
        integer :: i, s, a
        rmax = 2.0_dp
        call triplet_init(t, 8, rmax, 8)
        exact = 0
        do i = 1, 20000
            d = (i - 0.5_dp)*rmax/20000
            exact = exact + 4*pi*d**2*pi*(4*rmax + d)*(2*rmax - d)**2/12*rmax/20000
        end do
        print '(a,es10.3,a,es10.3)', 'triplet measure rel. error ', &
            abs(sum(t%volume)/exact - 1), ' indicator ', &
            abs(sum(t%volume_coarse)/exact - 1)
        call check(abs(sum(t%volume)/exact - 1) < 2.0e-3_dp, 'triplet bin measure')
        ! Per-bin measure (normalizes g3 directly): sub = 32 against sub = 96.
        call triplet_init(fine, 8, rmax, 32)
        call triplet_init(ref, 8, rmax, 96)
        worst = maxval(abs(fine%volume - ref%volume)/max(ref%volume, 1.0e-300_dp), &
                       mask=ref%volume > 1.0e-12_dp)
        print '(a,es10.3)', 'per-bin triplet measure rel. error (sub 32) ', worst
        call check(worst < 2.0e-3_dp, 'per-bin triplet measure')
        call rng_seed(g, 3_int64, status)
        length = mc_wigner_seitz_length(300)
        do s = 1, 200
            prev = sum(t%counts)
            do i = 1, 300
                do a = 1, 3
                    call rng_uniform(g, x(a, i))
                end do
            end do
            x = x*length
            call triplet_accumulate(t, x, 300, length)
            per(s) = (sum(t%counts) - prev)/(300.0_dp*299*298/length**6*sum(t%volume))
        end do
        call triplet_distribution(t, t%counts, t%snapshots, 300, length, g3)
        ratio = sum(t%counts)/(t%snapshots*300.0_dp*299*298/length**6*sum(t%volume))
        sigma = sqrt(sum((per - ratio)**2)/(200*199))
        print '(a,f9.5,a,es9.2)', 'ideal-gas triplet count ratio ', ratio, &
            ' +- ', sigma
        call check(abs(ratio - 1) < 4*sigma, 'ideal-gas g3 normalization')
        do i = 1, 40
            rg(i) = (i - 0.5_dp)*0.1_dp
        end do
        one = 1
        call triplet_references(t, rg, one, kirk, lin, cen)
        call check(maxval(abs(merge(kirk - 1, 0.0_dp, t%volume > 0))) < 1.0e-12_dp .and. &
                   maxval(abs(merge(lin - 1, 0.0_dp, t%volume > 0))) < 1.0e-12_dp, &
                   'ideal references K = P = 1')
    end subroutine

    !> First-order YBG closure with h(t) = A exp(-t^2): n int F(r13) h(r23)
    !> = -d/dR [n Gamma A pi^(3/2) erf(R)/R] (potential of a Gaussian cloud).
    subroutine ybg_gaussian()
        integer, parameter :: m = 800
        real(dp) :: rg(m), g(m), force(m), exact, worst, r, n
        real(dp), parameter :: gamma = 0.7_dp, amp = 0.3_dp
        integer :: k
        n = 3/(4*pi)
        do k = 1, m
            rg(k) = (k - 0.5_dp)*0.01_dp
            g(k) = 1 + amp*exp(-rg(k)**2)
        end do
        call ybg_closure_force(rg, g, gamma, CLOSURE_LINEAR, 0.005_dp, 48, force)
        worst = 0
        do k = 25, 400, 25
            r = rg(k)
            exact = -n*gamma*amp*pi**1.5_dp*(2/sqrt(pi)*exp(-r*r)/r - erf(r)/r**2)
            worst = max(worst, abs(force(k) - exact))
        end do
        print '(a,es10.3)', 'YBG first-order closure vs Gaussian-cloud force: ', worst
        call check(worst < 2.0e-4_dp*gamma*amp, 'YBG closure quadrature')
    end subroutine

    subroutine hnc_oracles()
        type(hnc_t) :: h
        logical :: ok
        real(dp) :: a(63), b(63), direct, worst, kappa
        integer :: i, j
        call hnc_init(h, 64, 10.0_dp, 1.0_dp)
        do i = 1, 63
            a(i) = cos(0.37_dp*i) + 0.1_dp*i
        end do
        call hnc_dst(h, a, b)
        worst = 0
        do j = 1, 63
            direct = 0
            do i = 1, 63
                direct = direct + a(i)*sin(pi*i*j/64.0_dp)
            end do
            worst = max(worst, abs(direct - b(j)))
        end do
        call check(worst < 1.0e-11_dp, 'DST-I by FFT')
        ! Caillol & Gilles, J. Stat. Phys. 100, 933 (2000)/2010 HNC values of
        ! -u/Gamma as tabulated by Demyanov & Levashov (arXiv:2509.02390).
        ! Gamma = 0.01: the core exp(-Gamma/r) needs dr < Gamma; reference
        ! is the Ortner (1999) cluster expansion 0.086193 (same table).
        call hnc_init(h, 32768, 100.0_dp, 1.0_dp)
        call hnc_solve(h, 0.01_dp, 1.0e-12_dp, 2000, ok)
        call ch_debye_wavenumber(0.01_dp, kappa)
        print '(a,es14.6,a,i0,a)', 'HNC u(0.01) ', h%energy, ' (', h%iterations, ' it)'
        call check(ok .and. abs(h%energy + 0.086193e-2_dp) < 1.0e-7_dp, 'HNC Gamma=0.01')
        call hnc_init(h, 16384, 60.0_dp, 1.0_dp)
        call hnc_solve(h, 0.1_dp, 1.0e-12_dp, 2000, ok)
        print '(a,f12.7)', 'HNC u(0.1) ', h%energy
        call check(ok .and. abs(h%energy + 0.0256885_dp) < 2.0e-6_dp, 'HNC Gamma=0.1')
        call hnc_solve(h, 1.0_dp, 1.0e-12_dp, 2000, ok)
        print '(a,f12.7)', 'HNC u(1) ', h%energy
        call check(ok .and. abs(h%energy + 0.570455_dp) < 2.0e-5_dp, 'HNC Gamma=1')
        call check(abs(h%s(1) - h%k(1)**2/(h%k(1)**2 + 3)) < 0.05_dp*h%s(1), &
                   'HNC S(k->0) Debye perfect screening')
    end subroutine

    subroutine statistics_and_fit()
        type(rng_t) :: g
        type(fortnum_status_t) :: status
        integer, parameter :: n = 200000
        real(dp), allocatable :: x(:)
        real(dp) :: e, mean, err, tau, ess, bm, berr, u, du, up, um, kappa
        real(dp), parameter :: rho = 0.8_dp
        integer :: i
        allocate (x(n))
        call rng_seed(g, 9_int64, status)
        x(1) = 0
        do i = 2, n
            call rng_normal(g, e)
            x(i) = rho*x(i - 1) + e
        end do
        call series_autocorrelation(x, mean, err, tau, ess)
        call series_batch_means(x, 40, bm, berr)
        print '(a,f7.3,a,f7.3,a,es10.3,a,es10.3)', 'AR(1) tau ', tau, ' exact ', &
            (1 + rho)/(2*(1 - rho)), ' err ', err, ' batch ', berr
        call check(abs(tau/4.5_dp - 1) < 0.1_dp, 'AR(1) autocorrelation time')
        call check(abs(err/sqrt(2*4.5_dp/(1 - rho**2)/n) - 1) < 0.1_dp, 'AR(1) error')
        call check(abs(berr/err - 1) < 0.35_dp, 'batch means vs autocorrelation error')
        call ocp_energy_literature(1.0_dp, u, du)
        call check(abs(u + 0.571414_dp) < 1.0e-3_dp, 'PC2000 fit at Gamma = 1')
        call ocp_energy_literature(3.0_dp + 1.0e-5_dp, up, e)
        call ocp_energy_literature(3.0_dp - 1.0e-5_dp, um, e)
        call ocp_energy_literature(3.0_dp, u, du)
        call check(abs(du - (up - um)/2.0e-5_dp) < 1.0e-8_dp, 'fit JVP leaf vs FD')
        call ch_debye_hueckel_energy(0.2_dp, u)
        call ch_debye_wavenumber(0.2_dp, kappa)
        call check(abs(u + sqrt(3.0_dp)/2*0.2_dp**1.5_dp) < 1.0e-15_dp .and. &
                   abs(kappa - sqrt(0.6_dp)) < 1.0e-15_dp, 'Debye leaves')
    end subroutine
end program test_correlation_hierarchy
