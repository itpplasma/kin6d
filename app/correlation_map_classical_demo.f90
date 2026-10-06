!> H1a: classical correlation-hierarchy map of the 3D OCP versus Gamma
!> (docs/correlation-hierarchy.md). Canonical Metropolis MC (incremental
!> Ewald, kin6d_ocp_monte_carlo) at each Gamma and two N; one MD cross-check
!> (kin6d_ocp, PME + encounter splitting, thermostat only in equilibration).
!> Estimators: g2, S(k), g3, cumulants c3/U3, Kirkwood defect, conditional
!> mean force; YBG closures (exact force, Kirkwood, first order); HNC and
!> Debye-Hueckel reduced models; fluctuation derivative du/dGamma.
!> Intervals are statistical (batch means / jackknife over batches).
!> Writes output/correlation-map-*.csv, output/correlation_map_classical_demo.nc
!> and output/correlation-map-*.png.
program correlation_map_classical_demo
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use kin6d_ocp_monte_carlo
    use kin6d_correlation_hierarchy
    use kin6d_ocp_hnc
    use kin6d_series_statistics
    use kin6d_pair_correlation, only: pair_counts, pair_distribution, structure_factor
    use kin6d_correlation_hierarchy_leaves, only: ch_debye_wavenumber, &
        ch_debye_hueckel_g, ch_debye_hueckel_linear_g, ch_debye_hueckel_energy, &
        ch_energy_fluctuation_derivative
    use kin6d_result_io, only: run_t, result_file_t
    use correlation_map_md, only: md_result_t, md_run
    use correlation_map_plots, only: plot_all
    use correlation_map_types
    implicit none
    type(run_t) :: run
    type(result_file_t) :: result
    type(chain_t), allocatable :: chains(:)
    type(md_result_t) :: md
    type(hnc_t), allocatable :: hnc(:)
    integer :: config_unit, ios, i, ng, task, status
    character(256) :: msg
    logical :: ok
    real(dp) :: t0, t1
    ! Configuration (namelist /correlation_map/).
    real(dp) :: gammas(6) = [0.1_dp, 0.3_dp, 1.0_dp, 3.0_dp, 10.0_dp, 30.0_dp]
    integer :: particles = 1000, particles_small = 500
    integer :: equilibration = 2000, production = 10000, batches = 20
    integer :: sample_every = 5, force_every = 10
    real(dp) :: tolerance = 1.0e-5_dp, triplet_rmax = 3.0_dp, pair_width = 0.05_dp
    integer :: triplet_bins = 15, triplet_sub = 32
    integer(int64) :: seed = 20261006_int64
    real(dp) :: md_gamma = 1.0_dp, md_equilibration = 50.0_dp, md_time = 1000.0_dp
    integer :: md_particles = 500
    !> Re-analyse the raw accumulators of a previous run instead of sampling.
    logical :: reuse_samples = .false.
    namelist /correlation_map/ gammas, particles, particles_small, equilibration, &
        production, batches, sample_every, force_every, tolerance, triplet_rmax, &
        pair_width, triplet_bins, triplet_sub, seed, md_gamma, md_equilibration, &
        md_time, md_particles, reuse_samples

    call run%start('correlation_map_classical_demo')
    if (run%argument_count() > 0) error stop 'usage: correlation_map_classical_demo [--config FILE]'
    ios = 0; msg = ''
    call run%open_config(config_unit)
    if (run%has_config()) read (config_unit, nml=correlation_map, iostat=ios, iomsg=msg)
    call run%close_config(config_unit, 'correlation_map', ios, msg)
    if (.not. (all(gammas > 0) .and. particles >= 50 .and. particles_small >= 50 .and. &
               production >= batches*sample_every .and. batches >= 4 .and. &
               mod(production, batches) == 0 .and. tolerance > 0 .and. &
               triplet_rmax > 0 .and. md_time > 0)) &
        error stop 'correlation_map: invalid configuration'
    call run%open_resolved(config_unit)
    write (config_unit, nml=correlation_map)
    call run%close_resolved(config_unit)
    call execute_command_line('mkdir -p output', exitstat=status)

    ng = size(gammas)
    allocate (chains(2*ng), hnc(ng))
    do i = 1, ng
        chains(i)%gamma = gammas(i); chains(i)%n = particles
        chains(ng + i)%gamma = gammas(i); chains(ng + i)%n = particles_small
    end do
    call cpu_time(t0)
    if (reuse_samples) then
        call raw_io(chains, md, .false.)
    else
        ! Independent chains (and the MD member) on OpenMP threads; the
        ! library is reentrant (no module state).
        !$omp parallel do schedule(dynamic, 1)
        do task = 1, 2*ng + 1
            if (task == 2*ng + 1) then
                call md_run(md_gamma, md_particles, md_equilibration, md_time, &
                            pair_width, seed + 999, md)
            else
                call run_chain(chains(task), seed + task)
            end if
        end do
        !$omp end parallel do
        call raw_io(chains, md, .true.)
    end if
    call cpu_time(t1)
    print '(a,f10.1,a)', 'sampling cpu ', t1 - t0, ' s'
    ! HNC reduced model by continuation in Gamma (ascending order assumed).
    do i = 1, ng
        if (i == 1) call hnc_init(hnc(i), 16384, 60.0_dp, 1.0_dp)
        if (i > 1) hnc(i) = hnc(i - 1)
        call hnc_solve(hnc(i), gammas(i), 1.0e-12_dp, 5000, ok)
        if (.not. ok) error stop 'HNC did not converge'
    end do
    !$omp parallel do schedule(dynamic, 1)
    do task = 1, 2*ng
        call analyze(chains(task), hnc(1 + mod(task - 1, ng)))
    end do
    !$omp end parallel do
    call write_csv(chains, ng, md)
    call write_result(chains, ng, md)
    call plot_all(chains, ng, hnc, md)
    print '(a)', 'wrote output/correlation-map-*.csv/.png and correlation_map_classical_demo.nc'
contains
    !> Raw accumulators to/from output/correlation-map-raw.bin (stream).
    subroutine raw_io(chains, md, save)
        type(chain_t), intent(inout) :: chains(:)
        type(md_result_t), intent(inout) :: md
        logical, intent(in) :: save
        integer :: unit, i, nb, m, nt, ne, np, nr
        if (save) then
            open (newunit=unit, file='output/correlation-map-raw.bin', access='stream', &
                  form='unformatted', status='replace')
            do i = 1, size(chains)
                associate (c => chains(i))
                    write (unit) c%gamma, c%n, c%bins, size(c%equil), size(c%energy), &
                        size(c%snaps), c%tri%bins, c%length, c%alpha, c%nk, c%bound, &
                        c%rmax_g, c%acceptance, c%step, c%cpu, c%energy_drift
                    write (unit) c%equil, c%energy, c%pair, c%pair_force, c%force_sum, &
                        c%triplet, c%sk, c%sk_vectors, c%snaps, c%fsnaps
                end associate
            end do
            write (unit) md%gamma, md%dt, md%temperature, md%temperature_err, md%u, &
                md%u_err, md%tau, md%drift, md%wall, md%a, md%n, md%steps, size(md%r), &
                md%failed_step, md%failure
            write (unit) md%r, md%g
            close (unit)
            return
        end if
        open (newunit=unit, file='output/correlation-map-raw.bin', access='stream', &
              form='unformatted', status='old', action='read')
        do i = 1, size(chains)
            associate (c => chains(i))
                read (unit) c%gamma, c%n, m, ne, np, nb, nt, c%length, c%alpha, c%nk, &
                    c%bound, c%rmax_g, c%acceptance, c%step, c%cpu, c%energy_drift
                if (nb /= batches .or. nt /= triplet_bins .or. np /= production) &
                    error stop 'reuse_samples: raw file does not match configuration'
                c%bins = m
                allocate (c%equil(ne), c%energy(np), c%pair(m, nb), c%pair_force(m, nb), &
                          c%force_sum(m, nb), c%triplet(nt, nt, nt, nb), &
                          c%sk(c%shells, nb), c%sk_vectors(c%shells), c%snaps(nb), &
                          c%fsnaps(nb))
                read (unit) c%equil, c%energy, c%pair, c%pair_force, c%force_sum, &
                    c%triplet, c%sk, c%sk_vectors, c%snaps, c%fsnaps
                call triplet_init(c%tri, nt, triplet_rmax, triplet_sub)
            end associate
        end do
        read (unit) md%gamma, md%dt, md%temperature, md%temperature_err, md%u, &
            md%u_err, md%tau, md%drift, md%wall, md%a, md%n, md%steps, nr, &
            md%failed_step, md%failure
        allocate (md%r(nr), md%g(nr))
        read (unit) md%r, md%g
        close (unit)
    end subroutine

    !> Equilibration (step tuning in its first half only), then production
    !> with per-batch accumulation of all estimators.
    subroutine run_chain(c, chain_seed)
        type(chain_t), intent(inout) :: c
        integer(int64), intent(in) :: chain_seed
        type(ocp_mc_t) :: mc
        logical :: ok
        integer :: s, b, len, nt
        integer(int64) :: c0, c1, rate
        real(dp), allocatable :: f(:, :), flat(:), sbuf(:)
        integer, allocatable :: vbuf(:)
        call system_clock(c0, rate)
        call mc_setup(mc, c%n, c%gamma, tolerance, chain_seed, ok)
        if (.not. ok) error stop 'MC setup failed'
        c%length = mc%length; c%alpha = mc%ewald%alpha; c%nk = mc%nk
        c%bound = mc%ewald%potential_bound
        c%bins = int(c%length/2/pair_width); c%rmax_g = c%bins*pair_width
        nt = triplet_bins
        allocate (c%equil(equilibration), c%energy(production))
        allocate (c%pair(c%bins, batches), c%pair_force(c%bins, batches), &
                  c%force_sum(c%bins, batches), c%triplet(nt, nt, nt, batches))
        allocate (c%sk(c%shells, batches), c%sk_vectors(c%shells), c%snaps(batches), &
                  c%fsnaps(batches))
        c%pair = 0; c%pair_force = 0; c%force_sum = 0; c%triplet = 0; c%sk = 0
        c%snaps = 0; c%fsnaps = 0
        call triplet_init(c%tri, nt, triplet_rmax, triplet_sub)
        allocate (f(3, c%n), flat(3*c%n), sbuf(c%shells), vbuf(c%shells))
        mc%step = 0.5_dp
        do s = 1, equilibration
            call mc_sweep(mc, s <= equilibration/2)
            c%equil(s) = c%gamma*mc%energy/c%n
        end do
        len = production/batches
        do s = 1, production
            b = 1 + (s - 1)/len
            call mc_sweep(mc, .false.)
            c%energy(s) = c%gamma*mc%energy/c%n
            if (mod(s, sample_every) == 0) then
                flat = reshape(mc%x, [3*c%n])
                call pair_counts(flat, c%n, c%length, c%rmax_g, c%pair(:, b))
                c%tri%counts = 0
                call triplet_accumulate(c%tri, mc%x, c%n, c%length)
                c%triplet(:, :, :, b) = c%triplet(:, :, :, b) + c%tri%counts
                call structure_factor(flat, c%n, c%length, 4, sbuf, vbuf)
                c%sk(:, b) = c%sk(:, b) + sbuf; c%sk_vectors = vbuf
                c%snaps(b) = c%snaps(b) + 1
            end if
            if (mod(s, force_every) == 0) then
                flat = reshape(mc%x, [3*c%n])
                call mc_forces(mc, f, ok)
                if (.not. ok) error stop 'MC forces'
                call mean_force_accumulate(mc%x, f, c%n, c%length, c%rmax_g, &
                                           c%force_sum(:, b))
                call pair_counts(flat, c%n, c%length, c%rmax_g, c%pair_force(:, b))
                c%fsnaps(b) = c%fsnaps(b) + 1
            end if
        end do
        c%energy_drift = abs(mc_energy(mc) - mc%energy)*c%gamma/c%n
        c%acceptance = mc_acceptance(mc); c%step = mc%step
        call system_clock(c1)
        c%cpu = real(c1 - c0, dp)/rate
        print '(a,f6.2,a,i0,a,f8.1,a,f6.3)', 'chain Gamma=', c%gamma, ' N=', c%n, &
            ' wall ', c%cpu, ' s, acceptance ', c%acceptance
    end subroutine

    !> Energy statistics, S(k), and jackknife over batches of every
    !> configuration-derived output (evaluate).
    subroutine analyze(c, h)
        type(chain_t), intent(inout) :: c
        type(hnc_t), intent(in) :: h
        real(dp) :: um, v(batches), outs(nouts, batches), mean(nouts), kappa, dr
        real(dp) :: gb(c%bins), rb(c%bins), gsum(c%bins), gsq(c%bins), ideal, k
        integer :: b, len, a, m
        logical :: use(batches)
        len = production/batches
        call series_autocorrelation(c%energy, c%u, c%u_err, c%tau, c%ess)
        call series_batch_means(c%energy, batches, um, c%u_batch_err)
        do b = 1, batches
            associate (x => c%energy((b - 1)*len + 1:b*len))
                v(b) = c%n*sum((x - sum(x)/len)**2)/len
            end associate
        end do
        c%var = sum(v)/batches
        c%var_err = sqrt(sum((v - c%var)**2)/(batches*(batches - 1)))
        call ch_energy_fluctuation_derivative(c%u, c%var, c%gamma, c%dudg)
        c%dudg_err = sqrt(max(c%u_err, c%u_batch_err)**2 + c%var_err**2)/c%gamma
        call ch_debye_wavenumber(c%gamma, kappa)
        allocate (c%r(c%bins), c%g(c%bins), c%g_err(c%bins), c%g_exact(c%bins), &
                  c%g_kirk(c%bins), c%g_lin(c%bins), c%g_hnc(c%bins), c%g_dh(c%bins), &
                  c%g_dhl(c%bins))
        gsum = 0; gsq = 0
        do b = 1, batches
            call pair_distribution(c%pair(:, b), c%n, c%length, c%rmax_g, c%snaps(b), &
                                   rb, gb)
            gsum = gsum + gb; gsq = gsq + gb**2
        end do
        c%g_err = sqrt(max(0.0_dp, gsq/batches - (gsum/batches)**2)/(batches - 1))
        c%r = rb
        call ch_debye_hueckel_g(c%r, c%gamma, kappa, c%g_dh)
        call ch_debye_hueckel_linear_g(c%r, c%gamma, kappa, c%g_dhl)
        do a = 1, c%bins
            c%g_hnc(a) = pair_value(h%r, h%g, c%r(a))
        end do
        dr = h%r(2) - h%r(1)
        c%u_tail = 1.5_dp*c%gamma*sum(h%r*(h%g - 1), mask=h%r > c%rmax_g)*dr
        c%u_hnc = h%energy
        call ch_debye_hueckel_energy(c%gamma, c%u_dh)
        ! Jackknife over batches.
        use = .true.
        call evaluate(c, h, use, c%out, .true.)
        do b = 1, batches
            use = .true.; use(b) = .false.
            call evaluate(c, h, use, outs(:, b), .false.)
        end do
        mean = sum(outs, dim=2)/batches
        c%out_err = sqrt((batches - 1)*sum((outs - spread(mean, 2, batches))**2, &
                                           dim=2)/batches)
        c%out = batches*c%out - (batches - 1)*mean
        ! Equilateral g3 batch errors.
        m = c%tri%bins
        allocate (c%eq_g3_err(m))
        do a = 1, m
            ideal = real(c%n, dp)*(c%n - 1)*real(c%n - 2, dp)/c%length**6* &
                    c%tri%volume(a, a, a)
            v = c%triplet(a, a, a, :)/(c%snaps*ideal)
            c%eq_g3_err(a) = sqrt(sum((v - sum(v)/batches)**2)/(batches*(batches - 1)))
        end do
        ! Structure factor per lattice shell.
        allocate (c%s_k(0), c%s_mc(0), c%s_err(0), c%s_dh(0), c%s_hnc(0))
        do a = 1, c%shells
            if (c%sk_vectors(a) == 0) cycle
            k = 2*acos(-1.0_dp)*sqrt(real(a, dp))/c%length
            v = c%sk(a, :)/c%snaps
            c%s_k = [c%s_k, k]
            c%s_mc = [c%s_mc, sum(v)/batches]
            c%s_err = [c%s_err, sqrt(sum((v - sum(v)/batches)**2)/(batches*(batches - 1)))]
            c%s_dh = [c%s_dh, k**2/(k**2 + kappa**2)]
            c%s_hnc = [c%s_hnc, pair_value(h%k, h%s, k)]
        end do
    end subroutine

    !> All configuration-derived outputs from the batches in `use`; with
    !> `store`, the tables (g2 and closure g2, equilateral slices) are kept.
    subroutine evaluate(c, h, use, out, store)
        type(chain_t), intent(inout) :: c
        type(hnc_t), intent(in) :: h
        logical, intent(in) :: use(:), store
        real(dp), intent(out) :: out(:)
        integer :: m, b, k, kmin, snaps, a
        real(dp) :: counts(c%bins), r(c%bins), g(c%bins), fc(c%bins), fs(c%bins)
        real(dp) :: force(c%bins), w(c%bins), gex(c%bins), gk(c%bins), gl(c%bins)
        real(dp) :: gn(c%bins)
        real(dp), allocatable :: tc(:, :, :), g3(:, :, :), kp(:, :, :), lp(:, :, :)
        real(dp), allocatable :: cen(:, :, :)
        real(dp) :: dr, pi, n, wt, m0, m2, hh
        pi = acos(-1.0_dp); n = 3/(4*pi)
        m = c%tri%bins; dr = pair_width
        allocate (tc(m, m, m), g3(m, m, m), kp(m, m, m), lp(m, m, m), cen(m, m, m))
        counts = 0; fc = 0; fs = 0; tc = 0; snaps = 0
        do b = 1, size(use)
            if (.not. use(b)) cycle
            counts = counts + c%pair(:, b); fc = fc + c%pair_force(:, b)
            fs = fs + c%force_sum(:, b); tc = tc + c%triplet(:, :, :, b)
            snaps = snaps + c%snaps(b)
        end do
        call pair_distribution(counts, c%n, c%length, c%rmax_g, snaps, r, g)
        call triplet_distribution(c%tri, tc, snaps, c%n, c%length, g3)
        ! Canonical finite-N normalization: the histograms are normalized by
        ! N(N-1)/V^2 and N(N-1)(N-2)/V^3, the thermodynamic-limit densities
        ! are n^2 and n^3; with S(0) = 0 (OCP) this is the leading 1/N term.
        g = g*(c%n - 1)/real(c%n, dp)
        g3 = g3*(c%n - 1)*real(c%n - 2, dp)/real(c%n, dp)**2
        call triplet_references(c%tri, r, g, kp, lp, cen)
        ! Closure input: h shifted on r > rmax/2 so that the table satisfies
        ! perfect screening n int h d3r = -1 (the far field of a non-neutral
        ! hole is a spurious long-range force; the shift is reported as
        ! neutral_mc + 1).
        call pair_moments(r, g, c%gamma, m0, m2)
        gn = g
        where (r > c%rmax_g/2) gn = g + (-1 - m0)/(3*sum(r**2, mask=r > c%rmax_g/2)*dr)
        out = 0
        ! Sizes over the triplet domain (all sides < triplet_rmax).
        out(O_H2) = n*sum(4*pi*r**2*(g - 1)**2, mask=r < c%tri%rmax)*dr
        where (c%tri%volume <= 0) g3 = kp
        out(O_C3) = triplet_norm(c%tri, g3 - kp)
        out(O_U3) = triplet_norm(c%tri, g3 - lp)
        out(O_RED) = triplet_norm(c%tri, kp - lp)
        ! Core-excluded sizes (all sides >= 0.4a): separates the near-contact
        ! region, where the g2 table (width pair_width) resolves exp(-Gamma/r)
        ! poorly at weak coupling.
        where (c%tri%volume > 0) cen = 1
        cen(1:2, :, :) = 0; cen(:, 1:2, :) = 0; cen(:, :, 1:2) = 0
        out(O_C3OUT) = triplet_norm(c%tri, (g3 - kp)*cen)
        out(O_U3OUT) = triplet_norm(c%tri, (g3 - lp)*cen)
        out(O_REDOUT) = triplet_norm(c%tri, (kp - lp)*cen)
        ! Resolution indicator: ||c3|| from the same counts merged into
        ! 3x3x3 blocks (bin width 3 rmax/bins).
        out(O_C3COARSE) = coarse_c3(c, tc, snaps, r, g)
        call triplet_references(c%tri, r, g, kp, lp, cen)
        wt = sum(kp*c%tri%volume, mask=kp > 1.0e-3_dp)
        out(O_DK) = sqrt(sum((g3 - kp)**2/kp*c%tri%volume, mask=kp > 1.0e-3_dp)/wt)
        ! Exact YBG: conditional mean force where >= 20 pair samples.
        force = 0
        where (fc >= 20) force = c%gamma*fs/max(fc, 1.0_dp)
        kmin = c%bins
        do k = c%bins, 1, -1
            if (fc(k) < 20) exit
            kmin = k
        end do
        gex = 0
        call ybg_potential(r(kmin:), force(kmin:), c%gamma, .false., w(kmin:))
        gex(kmin:) = exp(-w(kmin:))
        call ybg_closure_force(r, gn, c%gamma, CLOSURE_KIRKWOOD, dr/2, 32, force)
        call ybg_potential(r, force, c%gamma, .true., w)
        gk = exp(-w)
        call ybg_closure_force(r, gn, c%gamma, CLOSURE_LINEAR, dr/2, 32, force)
        call ybg_potential(r, force, c%gamma, .true., w)
        gl = exp(-w)
        out(O_UG) = pair_excess_energy(r, g, c%gamma) + c%u_tail
        out(O_UEX) = pair_excess_energy(r, gex, c%gamma) + c%u_tail
        out(O_UK) = pair_excess_energy(r, gk, c%gamma) + c%u_tail
        out(O_UL) = pair_excess_energy(r, gl, c%gamma) + c%u_tail
        out(O_EEX) = l2(r, gex, g, dr); out(O_EK) = l2(r, gk, g, dr); out(O_EL) = l2(r, gl, g, dr)
        out(O_EHNC) = l2(r, c%g_hnc, g, dr); out(O_EDH) = l2(r, c%g_dh, g, dr)
        out(O_XEX) = maxval(abs(gex - g)); out(O_XK) = maxval(abs(gk - g))
        out(O_XL) = maxval(abs(gl - g)); out(O_XHNC) = maxval(abs(c%g_hnc - g))
        out(O_XDH) = maxval(abs(c%g_dh - g))
        call pair_moments(r, g, c%gamma, m0, m2); out(O_M0MC) = m0; out(O_M2MC) = m2
        call pair_moments(r, gk, c%gamma, m0, m2); out(O_M2K) = m2
        call pair_moments(r, gl, c%gamma, m0, m2); out(O_M2L) = m2
        ! First maximum of g beyond r = 1 (parabolic refinement).
        k = maxloc(g, dim=1, mask=r > 1)
        k = max(2, min(c%bins - 1, k))
        hh = g(k - 1) - 2*g(k) + g(k + 1)
        out(O_RPEAK) = r(k)
        if (hh < 0) out(O_RPEAK) = r(k) - dr*(g(k + 1) - g(k - 1))/(2*hh)
        out(O_GPEAK) = g(k)
        if (.not. store) return
        c%g = g; c%g_exact = gex; c%g_kirk = gk; c%g_lin = gl
        allocate (c%eq_r(m), c%eq_g3(m), c%eq_kirk(m), c%eq_lin(m))
        do a = 1, m
            c%eq_r(a) = (a - 0.5_dp)*c%tri%rmax/m
            c%eq_g3(a) = g3(a, a, a); c%eq_kirk(a) = kp(a, a, a); c%eq_lin(a) = lp(a, a, a)
        end do
        c%bin_indicator = triplet_norm(c%tri, kp - cen)
    end subroutine

    real(dp) function coarse_c3(c, tc, snaps, r, g) result(norm)
        type(chain_t), intent(inout) :: c
        real(dp), intent(in) :: tc(:, :, :), r(:), g(:)
        integer, intent(in) :: snaps
        real(dp), allocatable :: cc(:, :, :), g3(:, :, :), kp(:, :, :), lp(:, :, :), &
            cen(:, :, :)
        integer :: m, a, b, d
        m = c%tri%bins/3
        if (.not. allocated(c%coarse%volume)) call triplet_init(c%coarse, m, &
            3*m*c%tri%rmax/c%tri%bins, triplet_sub)
        allocate (cc(m, m, m), g3(m, m, m), kp(m, m, m), lp(m, m, m), cen(m, m, m))
        do d = 1, m
            do b = 1, m
                do a = 1, m
                    cc(a, b, d) = sum(tc(3*a - 2:3*a, 3*b - 2:3*b, 3*d - 2:3*d))
                end do
            end do
        end do
        call triplet_distribution(c%coarse, cc, snaps, c%n, c%length, g3)
        g3 = g3*(c%n - 1)*real(c%n - 2, dp)/real(c%n, dp)**2
        call triplet_references(c%coarse, r, g, kp, lp, cen)
        norm = triplet_norm(c%coarse, g3 - kp)
    end function

    !> sqrt(n int (x - g)^2 d3r) on the pair table.
    real(dp) function l2(r, x, g, dr)
        real(dp), intent(in) :: r(:), x(:), g(:), dr
        l2 = sqrt(3*sum(r**2*(x - g)**2)*dr)
    end function

    subroutine write_csv(chains, ng, md)
        type(chain_t), intent(in) :: chains(:)
        integer, intent(in) :: ng
        type(md_result_t), intent(in) :: md
        integer :: unit, i, k, j
        real(dp) :: ul, dl
        open (newunit=unit, file='output/correlation-map-summary.csv', status='replace')
        write (unit, '(a)') '# H1a classical OCP correlation map; statistical intervals '// &
            '(1 s.e.; autocorrelation, batch means, jackknife over batches).'
        write (unit, '(a)', advance='no') 'gamma,N,L,alpha,ewald_pair_bound,nk,acceptance,'// &
            'step,wall_s,u,u_err_tau,u_err_batch,tau,ess,var_betaU_per_N,var_err,'// &
            'dudg_fluct,dudg_err,u_lit_pc2000,dudg_lit,u_hnc,u_dh,u_tail_hnc,'// &
            'energy_drift,kirkwood_centre_minus_bin_sq'
        do j = 1, nouts
            write (unit, '(a)', advance='no') ','//trim(output_names(j))//','// &
                trim(output_names(j))//'_err'
        end do
        write (unit, '(a)') ''
        do i = 1, size(chains)
            associate (c => chains(i))
                call ocp_energy_literature(c%gamma, ul, dl)
                write (unit, '(f0.3,",",i0,",",3(es14.6,","),i0,19(",",es14.6))', &
                       advance='no') c%gamma, c%n, c%length, c%alpha, c%bound, c%nk, &
                    c%acceptance, c%step, c%cpu, c%u, c%u_err, c%u_batch_err, c%tau, &
                    c%ess, c%var, c%var_err, c%dudg, c%dudg_err, ul, dl, c%u_hnc, &
                    c%u_dh, c%u_tail, c%energy_drift, c%bin_indicator
                do j = 1, nouts
                    write (unit, '(2(",",es14.6))', advance='no') c%out(j), c%out_err(j)
                end do
                write (unit, '(a)') ''
            end associate
        end do
        close (unit)
        open (newunit=unit, file='output/correlation-map-g2.csv', status='replace')
        write (unit, '(a)') 'gamma,N,r,g_mc,g_err,g_ybg_exact,g_ybg_kirkwood,'// &
            'g_ybg_linear,g_hnc,g_dh,g_dh_linear'
        do i = 1, size(chains)
            associate (c => chains(i))
                do k = 1, c%bins
                    write (unit, '(f0.3,",",i0,9(",",es14.6))') c%gamma, c%n, c%r(k), &
                        c%g(k), c%g_err(k), c%g_exact(k), c%g_kirk(k), c%g_lin(k), &
                        c%g_hnc(k), c%g_dh(k), c%g_dhl(k)
                end do
            end associate
        end do
        close (unit)
        open (newunit=unit, file='output/correlation-map-triplet.csv', status='replace')
        write (unit, '(a)') 'gamma,N,r_equilateral,g3,g3_err,kirkwood,cumulant_truncation'
        do i = 1, size(chains)
            associate (c => chains(i))
                do k = 1, size(c%eq_r)
                    write (unit, '(f0.3,",",i0,5(",",es14.6))') c%gamma, c%n, c%eq_r(k), &
                        c%eq_g3(k), c%eq_g3_err(k), c%eq_kirk(k), c%eq_lin(k)
                end do
            end associate
        end do
        close (unit)
        open (newunit=unit, file='output/correlation-map-sk.csv', status='replace')
        write (unit, '(a)') 'gamma,N,k,s_mc,s_err,s_dh,s_hnc'
        do i = 1, size(chains)
            associate (c => chains(i))
                do k = 1, size(c%s_k)
                    write (unit, '(f0.3,",",i0,5(",",es14.6))') c%gamma, c%n, c%s_k(k), &
                        c%s_mc(k), c%s_err(k), c%s_dh(k), c%s_hnc(k)
                end do
            end associate
        end do
        close (unit)
        open (newunit=unit, file='output/correlation-map-equilibration.csv', status='replace')
        write (unit, '(a)') 'gamma,N,sweep,beta_U_per_N'
        do i = 1, size(chains)
            do k = 10, size(chains(i)%equil), 10
                write (unit, '(f0.3,",",i0,",",i0,",",es14.6)') chains(i)%gamma, &
                    chains(i)%n, k, chains(i)%equil(k)
            end do
        end do
        close (unit)
        open (newunit=unit, file='output/correlation-map-md.csv', status='replace')
        write (unit, '(a)') '# MD (PME, encounter splitting, NVE after thermostat) vs MC'
        write (unit, '(a,f0.3,a,i0,a,f0.5,a,i0,a,es12.4,a,f8.1)') '# gamma ', md%gamma, &
            ' N ', md%n, ' dt ', md%dt, ' steps ', md%steps, ' energy_drift_kT_per_particle ', &
            md%drift, ' wall_s ', md%wall
        write (unit, '(a,i0,a)') '# failed_step ', md%failed_step, ' '//trim(md%failure)
        write (unit, '(a,f10.6,a,es10.3,a,f10.6,a,es10.3,a,f8.2)') '# T/T_ref ', &
            md%temperature, ' +- ', md%temperature_err, ' u ', md%u, ' +- ', md%u_err, &
            ' tau_samples ', md%tau
        write (unit, '(a)') 'r,g_md'
        do k = 1, size(md%r)
            write (unit, '(es14.6,",",es14.6)') md%r(k), md%g(k)
        end do
        close (unit)
        ! Weak-coupling power laws of the sizes (N1 chains): ||X|| = A Gamma^p
        ! over Gamma <= 1 (Gamma <= 3 for c3, which is unresolved at 0.1).
        open (newunit=unit, file='output/correlation-map-scaling.csv', status='replace')
        write (unit, '(a)') 'quantity,gamma_max,exponent,exponent_err,amplitude,points'
        call fit_row(unit, 'h2', O_H2, 1.0_dp)
        call fit_row(unit, 'reducible', O_RED, 1.0_dp)
        call fit_row(unit, 'ursell_u3', O_U3, 1.0_dp)
        call fit_row(unit, 'connected_c3', O_C3, 3.0_dp)
        call fit_row(unit, 'connected_c3_strong', O_C3, 30.0_dp)
        close (unit)
        ! Finite-size check: two N and the linear 1/N extrapolation of u.
        open (newunit=unit, file='output/correlation-map-finite-size.csv', status='replace')
        write (unit, '(a)') 'gamma,u_N1,err_N1,u_N2,err_N2,u_extrap_1_over_N,err_extrap,'// &
            'u_lit_pc2000,c3_sq_N1,c3_sq_err_N1,c3_sq_N2,c3_sq_err_N2'
        do i = 1, ng
            associate (a => chains(i), b => chains(ng + i))
                call ocp_energy_literature(a%gamma, ul, dl)
                write (unit, '(f0.3,11(",",es14.6))') a%gamma, a%u, &
                    max(a%u_err, a%u_batch_err), b%u, max(b%u_err, b%u_batch_err), &
                    (a%n*a%u - b%n*b%u)/(a%n - b%n), &
                    sqrt((a%n*max(a%u_err, a%u_batch_err))**2 + &
                         (b%n*max(b%u_err, b%u_batch_err))**2)/abs(a%n - b%n), ul, &
                    a%out(O_C3), a%out_err(O_C3), b%out(O_C3), b%out_err(O_C3)
            end associate
        end do
        close (unit)
    end subroutine

    subroutine fit_row(unit, name, idx, gmax)
        integer, intent(in) :: unit
        character(*), intent(in) :: name
        integer, intent(in) :: idx
        real(dp), intent(in) :: gmax
        real(dp) :: y(ng), e(ng), p, a, pe
        integer :: i, used, m
        m = count(chains(1:ng)%gamma <= gmax)
        do i = 1, ng
            y(i) = sqrt(max(chains(i)%out(idx), 0.0_dp))
            e(i) = chains(i)%out_err(idx)/(2*max(y(i), 1.0e-300_dp))
        end do
        call power_fit(chains(1:m)%gamma, y(1:m), e(1:m), p, a, pe, used)
        write (unit, '(a,",",f0.2,3(",",es14.6),",",i0)') name, gmax, p, pe, a, used
    end subroutine

    subroutine write_result(chains, ng, md)
        type(chain_t), intent(in) :: chains(:)
        integer, intent(in) :: ng
        type(md_result_t), intent(in) :: md
        real(dp) :: table(nouts, size(chains)), errs(nouts, size(chains))
        integer :: i
        do i = 1, size(chains)
            table(:, i) = chains(i)%out; errs(:, i) = chains(i)%out_err
        end do
        call result%create(run)
        call result%write('gamma', chains(1:ng)%gamma, '1', 'coupling parameter', &
                          dims='gamma')
        call result%write('particles', real([chains(1)%n, chains(ng + 1)%n], dp), '1', &
                          'N of the large and small chains', dims='size')
        call result%write('u_excess', reshape(chains%u, [ng, 2]), '1', &
                          'beta U_ex/N (MC energy route)', dims='gamma,size', &
                          status='statistical')
        call result%write('u_excess_err', reshape(max(chains%u_err, chains%u_batch_err), &
                          [ng, 2]), '1', 'one standard error', dims='gamma,size')
        call result%write('dudgamma_fluctuation', reshape(chains%dudg, [ng, 2]), '1', &
                          'd u/d Gamma = (u - Var(beta U)/N)/Gamma', dims='gamma,size', &
                          status='statistical')
        call result%write('dudgamma_fluctuation_err', reshape(chains%dudg_err, [ng, 2]), &
                          '1', 'one standard error', dims='gamma,size')
        call result%write('hierarchy_outputs', reshape(table, [nouts, ng, 2]), '1', &
                          'jackknife outputs in output_names order (summary CSV header)', &
                          dims='output,gamma,size', status='statistical')
        call result%write('hierarchy_outputs_err', reshape(errs, [nouts, ng, 2]), '1', &
                          'jackknife standard errors', dims='output,gamma,size')
        call result%write('md_u_excess', md%u, '1', 'MD beta U_ex/N at T_ref', &
                          status='statistical')
        call result%write('md_temperature', md%temperature, '1', 'MD mean T/T_ref')
        call result%close()
    end subroutine
end program correlation_map_classical_demo
