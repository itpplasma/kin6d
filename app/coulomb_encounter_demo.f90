! Q1/Q1b evidence: Sommerfeld-eta and impact-parameter scans of the localized
! two-body quantum Coulomb encounter (product, singlet, triplet, classical),
! head-on and off-axis encounter movies of the relative density (singlet,
! triplet, spin-unresolved, classical ensemble) with a turning-point cut, and
! the hydrogen Kepler packet. Writes CSV and figures to output/
! and the evidence file output/coulomb_encounter_demo.nc (kin6d_result_io).
! Usage: coulomb_encounter_demo [--config FILE] [--no-movie | --movie-only]
! Namelist group &coulomb_encounter: etas (four Sommerfeld parameters), nclass
! (classical pairs per b-scan ensemble), nclass_scan (pairs per eta-scan
! ensemble), seed/movie_seed (classical ensembles), movie.
program coulomb_encounter_demo
    use, intrinsic :: iso_fortran_env, only: dp => real64, i64 => int64
    use kin6d_coulomb_encounter
    use kin6d_coulomb_kepler
    use kin6d_coulomb_encounter_leaves, only: ce_mott, ce_free_width
    use kin6d_result_io, only: run_t, result_file_t
    use fortplot, only: figure, get_global_figure, figure_t, subplot, plot, title, &
        xlabel, ylabel, xlim, ylim, legend, savefig, pcolormesh, tight_layout, &
        set_yscale, suptitle, colorbar, text, fill_between
    use fortplot_animation, only: animation_t, animation_func => FuncAnimation, &
        save_animation
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp)
    integer, parameter :: nb = 36
    real(dp), parameter :: bs(4) = [0.0_dp, 5.0_dp, 10.0_dp, 20.0_dp]
    real(dp) :: etas(4) = [0.3_dp, 1.0_dp, 3.0_dp, 10.0_dp]
    integer :: nclass = 20000, nclass_scan = 100000
    integer(i64) :: seed = 1234_i64, movie_seed = 77_i64
    logical :: movie = .true.
    namelist /coulomb_encounter/ etas, nclass, nclass_scan, seed, movie_seed, movie
    type(run_t) :: run
    type(result_file_t) :: result
    character(160) :: arg
    character(256) :: msg
    integer :: status, i, unit, ios
    integer, parameter :: nbf = 180
    real(dp) :: ratio(nbf, 4, 4), cratio(nb, 2, 4), theta(nb), theta_fine(nbf)
    ! Evidence for the result file, filled by the scans.
    real(dp) :: eta_summary(6, 4), bscan_quantum(nb, 4), bscan_classical(nb, 4)
    real(dp) :: kepler_errors(2), bin_vjp(4), cost_seconds(6)
    real(dp), allocatable :: kepler_series(:, :)
    ! Movie state shared with the animation callback (host association).
    integer, parameter :: nmove = 264, nhold = 24, npx = 240, ncx = 120
    real(dp), parameter :: mv_slab = 1.0_dp   ! classical slab thickness k h
    real(dp), allocatable :: mq(:, :, :, :), mc(:, :, :)
    real(dp) :: mtimes(nmove), mslow(nmove), xq(0:npx), xc(0:ncx), mv_lambdabar, &
        mv_speed, mv_z0, mv_b, mv_r0, mv_tc, mv_lmin, mv_lmax
    logical :: movie_only = .false.
    type(figure_t), pointer :: canvas => null()
    call run%start('coulomb_encounter_demo')
    call run%open_config(unit)
    ios = 0; msg = ''
    if (run%has_config()) read (unit, nml=coulomb_encounter, iostat=ios, iomsg=msg)
    call run%close_config(unit, 'coulomb_encounter', ios, msg)
    do i = 1, run%argument_count()
        call run%argument(i, arg)
        if (trim(arg) == '--no-movie') movie = .false.
        if (trim(arg) == '--movie-only') movie_only = .true.
    end do
    if (any(etas <= 0) .or. min(nclass, nclass_scan) < 1) &
        error stop 'etas > 0, nclass >= 1 and nclass_scan >= 1 required'
    call run%open_resolved(unit)
    write (unit, nml=coulomb_encounter)
    call run%close_resolved(unit)
    call execute_command_line('mkdir -p output', exitstat=status)
    if (status /= 0) error stop 'cannot create output/'
    if (movie_only) then
        call encounter_movies()
        stop
    end if
    call eta_scan()
    call b_scan()
    call kepler_case()
    call cost_ledger()
    if (movie) call encounter_movies()
    call write_result()
contains
    subroutine write_result()
        call result%create(run)
        call result%attribute('model', 'localized two-body quantum Coulomb '// &
            'encounter (partial waves, exact retained-spectrum evolution) and '// &
            'hydrogen radial Kepler packet; atomic units')
        call result%attribute('status', 'candidate finite-time wavepacket bins '// &
            'versus plane-wave references; state_bound is the retained-state '// &
            'error bound, no pointwise bin enclosure')
        call result%write('theta', theta, 'deg', 'CM scattering angle bin centre', &
            dims='bin')
        call result%write('theta_fine', theta_fine, 'deg', 'CM scattering angle fine '// &
            'bin centre', dims='fine_bin')
        call result%write('eta', etas, '1', 'Sommerfeld parameter', dims='eta')
        call result%write('eta_bin_ratio', ratio, '1', 'bin probability / plane-wave '// &
            'Rutherford: symmetrized product, singlet, triplet, Mott singlet', &
            dims='fine_bin,channel,eta', status='candidate')
        call result%write('eta_classical_ratio', cratio, '1', 'symmetrized classical '// &
            'bin probability / plane-wave Rutherford and its standard error', &
            dims='bin,estimate,eta', status='statistical (Monte Carlo, fixed seed)')
        call result%write('eta_summary', eta_summary, 'mixed', 'k, r0/sigma, '// &
            'state_bound, projection_defect, reconstruction_rms, wall_probability', &
            dims='summary,eta')
        call result%write('b_over_lambdabar', bs, '1', 'impact parameter (eta = 1)', &
            dims='impact')
        call result%write('b_scan_quantum', bscan_quantum, '1', &
            'product-state bin probability', dims='bin,impact', status='candidate')
        call result%write('b_scan_classical', bscan_classical, '1', &
            'classical-ensemble bin probability', dims='bin,impact', &
            status='statistical (Monte Carlo, fixed seed)')
        call result%write('kepler_series', kepler_series, 'mixed', 't/T_rev, '// &
            '<r>/r_c, |A|^2 numeric, |A|^2 exact levels', dims='kepler_time,series')
        call result%write('kepler_errors', kepler_errors, '1', &
            'max level relative error, max |A - A_exact|', dims='kepler_error')
        call result%write('bin_vjp', bin_vjp, 'mixed', 'bin probability (pi/3..2pi/3) '// &
            'and its VJP gradient d/db, d/dsigma, d/dg', dims='vjp')
        call result%write('cost_seconds', cost_seconds, 's', 'state, bins, state '// &
            'bound, xz projection, JVP, VJP wall time', dims='cost')
        call result%close()
    end subroutine write_result

    ! Encounter for one eta in reduced units hbar = m_e = e^2/(4 pi eps0) = 1
    ! (atomic units); lengths scale with lambdabar = 1/k.
    subroutine setup_case(eta, sig_l, z0_l, b_l, cfg, enc, t)
        real(dp), intent(in) :: eta, sig_l, z0_l, b_l
        type(encounter_config_t), intent(out) :: cfg
        type(encounter_t), intent(out) :: enc
        real(dp), intent(out) :: t
        character(:), allocatable :: msg
        logical :: ok
        cfg%g = 1; cfg%m1 = 1; cfg%m2 = 1
        cfg%k = 0.5_dp/eta
        cfg%sigma = sig_l/cfg%k; cfg%z0 = z0_l/cfg%k; cfg%b = b_l/cfg%k
        t = 2*cfg%z0/(cfg%hbar*cfg%k/0.5_dp)
        cfg%t_final = t
        call encounter_init(enc, cfg, ok, msg)
        if (.not. ok) error stop 'encounter setup failed: '//msg
    end subroutine setup_case

    subroutine classical_bins(enc, t, rcut, edges, nsamp, pc)
        type(encounter_t), intent(in) :: enc
        integer, intent(in) :: nsamp
        real(dp), intent(in) :: t, rcut, edges(0:)
        real(dp), intent(out) :: pc(:)
        real(dp), allocatable :: rel(:, :, :), cm(:, :, :)
        real(dp) :: drift, th
        integer :: n, ib
        logical :: ok
        allocate (rel(3, nsamp, 1), cm(3, nsamp, 1))
        call encounter_classical_ensemble(enc, nsamp, seed, [t], rel, cm, drift, ok)
        if (.not. ok) error stop 'classical ensemble failed'
        pc = 0
        do n = 1, nsamp
            if (norm2(rel(:, n, 1)) < rcut) cycle
            th = acos(rel(3, n, 1)/norm2(rel(:, n, 1)))
            ib = min(size(pc), int(th/(edges(1) - edges(0))) + 1)
            pc(ib) = pc(ib) + 1.0_dp/nsamp
        end do
    end subroutine classical_bins

    ! Sommerfeld-eta scan: finite-time angular bins versus plane-wave Mott and
    ! Rutherford references (fluence from the free transverse width at z0/v).
    subroutine eta_scan()
        type(encounter_config_t) :: cfg
        type(encounter_t) :: enc
        complex(dp), allocatable :: a(:, :, :)
        real(dp) :: edges(0:nb), fedges(0:nbf), pp(nbf), ps(nbf), pt(nbf), pc(nb), psym(nbf), &
            t, refc(nb), sig_l, width, nperp, ref, ms, eb, khi, uniform, dof, wall, &
            ymax, lo(nb), hi(nb)
        integer :: ie, ib, u, us, uc, l, i1, i2, c1, c2
        open (newunit=u, file='output/coulomb-eta-bins.csv', status='replace')
        write (u, '(a)') 'eta,theta_lo_deg,theta_hi_deg,P_product,P_product_sym,'// &
            'P_singlet,P_triplet,R_sym_fluence,mott_singlet_ratio,mott_triplet_ratio'
        open (newunit=uc, file='output/coulomb-eta-classical-bins.csv', status='replace')
        write (uc, '(a)') 'eta,theta_lo_deg,theta_hi_deg,P_classical_sym,'// &
            'standard_error,R_sym_fluence'
        open (newunit=us, file='output/coulomb-eta-summary.csv', status='replace')
        write (us, '(a)') 'eta,k,lambdabar,r0,sigma_over_lambdabar,r0_over_sigma,z0,'// &
            't_final,rmax,n_radial,lmax,mmax,kept_states,partial_wave_dof,'// &
            'uniform_grid_points,setup_s,eigen_s,projection_s,stored_MB,e0,rho,'// &
            'state_bound_t,projection_defect,reconstruction_rms,wall_prob,l_edge2,'// &
            'spectral_discarded2'
        do ib = 0, nb
            edges(ib) = pi*ib/nb
        end do
        do ib = 0, nbf
            fedges(ib) = pi*ib/nbf
        end do
        theta = (edges(1:) + edges(:nb - 1))/2*180/pi
        theta_fine = (fedges(1:) + fedges(:nbf - 1))/2*180/pi
        do ie = 1, size(etas)
            sig_l = max(10.0_dp, 4*etas(ie))
            call setup_case(etas(ie), sig_l, 6*sig_l, 0.0_dp, cfg, enc, t)
            if (allocated(a)) deallocate (a)
            allocate (a(enc%mesh%n, 0:enc%mmax, 0:enc%lmax))
            call encounter_state(enc, t, a)
            call encounter_angular_bins(enc, a, SECTOR_PRODUCT, cfg%z0/4, fedges, pp)
            call encounter_angular_bins(enc, a, SECTOR_SINGLET, cfg%z0/4, fedges, ps)
            call encounter_angular_bins(enc, a, SECTOR_TRIPLET, cfg%z0/4, fedges, pt)
            call classical_bins(enc, t, cfg%z0/4, edges, nclass_scan, pc)
            call ce_free_width(cfg%sigma, enc%mu, cfg%hbar, cfg%z0/enc%speed, width)
            nperp = 1/(2*pi*width**2)
            psym = (pp + pp(nbf:1:-1))/2
            do ib = 1, nbf
                call reference_bin(fedges(ib - 1), fedges(ib), enc%eta, cfg%k, nperp, ref, ms)
                write (u, '(f6.2,2(",",f8.3),7(",",es14.6))') etas(ie), &
                    fedges(ib - 1)*180/pi, fedges(ib)*180/pi, pp(ib), psym(ib), ps(ib), &
                    pt(ib), ref, ms/ref, 2 - ms/ref
                ratio(ib, :, ie) = [psym(ib)/ref, ps(ib)/ref, pt(ib)/ref, ms/ref]
            end do
            ! Classical: symmetrized counts in coarse bins, binomial standard error.
            do ib = 1, nb
                call reference_bin(edges(ib - 1), edges(ib), enc%eta, cfg%k, nperp, refc(ib), ms)
                cratio(ib, 1, ie) = (pc(ib) + pc(nb + 1 - ib))/2/refc(ib)
                cratio(ib, 2, ie) = sqrt((pc(ib) + pc(nb + 1 - ib))/nclass_scan)/2/refc(ib)
                write (uc, '(f6.2,2(",",f8.3),3(",",es14.6))') etas(ie), edges(ib - 1)*180/pi, &
                    edges(ib)*180/pi, cratio(ib, 1, ie)*refc(ib), cratio(ib, 2, ie)*refc(ib), &
                    refc(ib)
            end do
            khi = cfg%k + 8/(2*cfg%sigma)
            uniform = (2*enc%mesh%rmax*khi*cfg%nodes_per_wavelength/(2*pi))**3
            dof = 0
            do l = 0, enc%lmax
                dof = dof + real(enc%mesh%n*(min(l, enc%mmax) + 1), dp)
            end do
            call encounter_state_error(enc, t, 1.0_dp, 0.5_dp, eb, wall)
            wall = encounter_wall_probability(enc, a, SECTOR_PRODUCT, 2*cfg%sigma)
            eta_summary(:, ie) = [cfg%k, enc%r0/cfg%sigma, eb, &
                enc%ledger%projection_defect, enc%ledger%reconstruction_rms, wall]
            write (us, '(f6.2,26(",",es13.5))') etas(ie), cfg%k, enc%lambdabar, enc%r0, &
                sig_l, enc%r0/cfg%sigma, cfg%z0, t, enc%mesh%rmax, real(enc%mesh%n, dp), &
                real(enc%lmax, dp), real(enc%mmax, dp), real(enc%ledger%kept_states, dp), &
                dof, uniform, enc%ledger%setup_seconds, enc%ledger%eigen_seconds, &
                enc%ledger%projection_seconds, real(enc%ledger%stored_bytes, dp)/1.0e6_dp, &
                enc%ledger%e0, enc%ledger%residual_rate, eb, enc%ledger%projection_defect, &
                enc%ledger%reconstruction_rms, wall, enc%ledger%l_edge2, &
                enc%ledger%spectral_discarded2
            print '(a,f5.1,a,f6.2,a,i0,a,i0,a,f7.2,a,es9.2)', 'eta=', etas(ie), &
                ' r0/sigma=', enc%r0/cfg%sigma, ' n=', enc%mesh%n, ' lmax=', enc%lmax, &
                ' setup s=', enc%ledger%setup_seconds, ' e(t)=', eb
        end do
        close (u); close (uc); close (us)
        ! Plot window 50-130 deg (finite-time bins diverge from Rutherford at small
        ! angles): fine bins i1..i2, coarse classical bins c1..c2.
        i1 = 5*nbf/18 + 1; i2 = nbf - 5*nbf/18; c1 = 5*nb/18 + 1; c2 = nb - 5*nb/18
        call figure(figsize=[13.2_dp, 8.8_dp])
        do ie = 1, size(etas)
            call subplot(2, 2, ie)
            ymax = 1.08_dp*max(maxval(ratio(i1:i2, :, ie)), maxval(2 - ratio(i1:i2, 4, ie)), &
                               maxval(cratio(c1:c2, 1, ie)))
            if (ie == 1) ymax = 1.75_dp*ymax
            write (arg, '(a,f4.1,a,f5.2,a)') '$\eta$ = ', etas(ie), ',  $r_0/\sigma$ = ', &
                2*etas(ie)/max(10.0_dp, 4*etas(ie)), ',  1° bins'
            call title(trim(arg))
            call xlabel('CM scattering angle $\theta$ (deg), at $t = 2 z_0/v$')
            call ylabel('bin probability / plane-wave Rutherford')
            call xlim(50.0_dp, 130.0_dp); call ylim(0.0_dp, ymax)
            lo = max(0.0_dp, min(ymax, cratio(:, 1, ie) - cratio(:, 2, ie)))
            hi = max(0.0_dp, min(ymax, cratio(:, 1, ie) + cratio(:, 2, ie)))
            call fill_between(theta(c1:c2), lo(c1:c2), hi(c1:c2), &
                              color=[0.0_dp, 0.62_dp, 0.45_dp], alpha=0.25_dp)
            call plot(theta_fine(i1:i2), ratio(i1:i2, 4, ie), label='Mott singlet (plane wave)', &
                      linestyle='-', color=[0.6_dp, 0.6_dp, 0.6_dp])
            call plot(theta_fine(i1:i2), 2 - ratio(i1:i2, 4, ie), &
                      label='Mott triplet (plane wave)', linestyle='--', &
                      color=[0.6_dp, 0.6_dp, 0.6_dp])
            call plot(theta_fine(i1:i2), ratio(i1:i2, 2, ie), label='quantum singlet', &
                      linestyle='-', color=[0.0_dp, 0.447_dp, 0.698_dp])
            call plot(theta_fine(i1:i2), ratio(i1:i2, 3, ie), label='quantum triplet', &
                      linestyle='-', color=[0.835_dp, 0.369_dp, 0.0_dp])
            call plot(theta_fine(i1:i2), ratio(i1:i2, 1, ie), &
                      label='quantum product, spin-unresolved', linestyle='-', &
                      color=[0.0_dp, 0.0_dp, 0.0_dp])
            call plot(theta(c1:c2), min(ymax, cratio(c1:c2, 1, ie)), &
                      label='classical ensemble, 5° bins $\pm 1$ s.e.', &
                      linestyle=':', color=[0.0_dp, 0.62_dp, 0.45_dp])
            if (ie == 1) call legend(loc='upper center')
        end do
        call tight_layout()
        call savefig('output/coulomb-mott-eta-scan.png')
    end subroutine eta_scan

    ! Plane-wave Rutherford (spin-averaged, no cross term) and Mott singlet bin
    ! fluences for a CM angle bin, 16-point Gauss-Legendre in cos theta.
    subroutine reference_bin(th1, th2, eta, k, nperp, ref, ms)
        real(dp), intent(in) :: th1, th2, eta, k, nperp
        real(dp), intent(out) :: ref, ms
        real(dp) :: xq(16), wq(16), d, ex, cr
        integer :: q
        call gauss_nodes(th1, th2, xq, wq)
        ref = 0; ms = 0
        do q = 1, 16
            call ce_mott(acos(xq(q)), eta, k, d, ex, cr)
            ref = ref + 2*pi*wq(q)*(d + ex)/2*nperp
            ms = ms + 2*pi*wq(q)*(d + ex + cr)/2*nperp
        end do
    end subroutine reference_bin

    subroutine gauss_nodes(th1, th2, xq, wq)
        use fortnum_quadrature, only: gauss_legendre_ab
        real(dp), intent(in) :: th1, th2
        real(dp), intent(out) :: xq(:), wq(:)
        call gauss_legendre_ab(size(xq), cos(th2), cos(th1), xq, wq)
    end subroutine gauss_nodes

    ! Impact-parameter scan at eta = 1 (distinguishable product state).
    subroutine b_scan()
        type(encounter_config_t) :: cfg
        type(encounter_t) :: enc
        complex(dp), allocatable :: a(:, :, :)
        real(dp) :: edges(0:nb), pq(nb, 4), pc(nb, 4), t
        integer :: ib, ibs, u
        do ib = 0, nb
            edges(ib) = pi*ib/nb
        end do
        open (newunit=u, file='output/coulomb-b-scan.csv', status='replace')
        write (u, '(a)') 'b_over_lambdabar,b_over_sigma,theta_lo_deg,theta_hi_deg,'// &
            'P_quantum_product,P_classical,lmax,mmax'
        do ibs = 1, size(bs)
            call setup_case(1.0_dp, 10.0_dp, 60.0_dp, bs(ibs), cfg, enc, t)
            if (allocated(a)) deallocate (a)
            allocate (a(enc%mesh%n, 0:enc%mmax, 0:enc%lmax))
            call encounter_state(enc, t, a)
            call encounter_angular_bins(enc, a, SECTOR_PRODUCT, cfg%z0/4, edges, pq(:, ibs))
            call classical_bins(enc, t, cfg%z0/4, edges, nclass, pc(:, ibs))
            do ib = 1, nb
                write (u, '(f6.2,",",f6.3,2(",",f8.3),2(",",es14.6),2(",",i0))') bs(ibs), &
                    bs(ibs)/10, edges(ib - 1)*180/pi, edges(ib)*180/pi, pq(ib, ibs), &
                    pc(ib, ibs), enc%lmax, enc%mmax
            end do
            print '(a,f5.1,a,i0,a,i0)', 'b/lambdabar=', bs(ibs), ' lmax=', enc%lmax, &
                ' mmax=', enc%mmax
        end do
        close (u)
        bscan_quantum = pq; bscan_classical = pc
        call figure(figsize=[12.8_dp, 6.4_dp])
        call title('Impact-parameter scan, $\eta = 1$, $\sigma = 10/k$ (distinguishable pair)')
        call xlabel('CM scattering angle $\theta$ (deg), at $t = 2 z_0/v$')
        call ylabel('bin probability (5° bins, $|r| > z_0/4$)')
        call set_yscale('log')
        call xlim(0.0_dp, 180.0_dp); call ylim(1.0e-6_dp, 1.0_dp)
        do ibs = 1, size(bs)
            write (arg, '(a,i0)') 'quantum, $kb$ = ', nint(bs(ibs))
            call plot(theta, max(pq(:, ibs), 1.0e-6_dp), label=trim(arg), linestyle='-', &
                      color=palette(ibs))
            write (arg, '(a,i0,a,i0,a)') 'classical, $kb$ = ', nint(bs(ibs)), ' (', nclass, &
                ' pairs)'
            call plot(pack(theta, pc(:, ibs) > 0), max(pack(pc(:, ibs), pc(:, ibs) > 0), &
                      1.0e-6_dp), label=trim(arg), linestyle='none', marker='o', &
                      markersize=4.0_dp, color=palette(ibs))
        end do
        call legend(loc='upper right')
        call savefig('output/coulomb-b-scan.png')
    end subroutine b_scan

    function palette(i) result(c)
        integer, intent(in) :: i
        real(dp) :: c(3)
        real(dp), parameter :: p(3, 4) = reshape([0.0_dp, 0.447_dp, 0.698_dp, &
            0.835_dp, 0.369_dp, 0.0_dp, 0.0_dp, 0.62_dp, 0.45_dp, 0.8_dp, 0.475_dp, &
            0.655_dp], [3, 4])
        c = p(:, 1 + mod(i - 1, 4))
    end function palette

    ! Q1b Rydberg radial packet: breathing, collapse, revivals; exact-level oracle.
    subroutine kepler_case()
        type(kepler_t) :: kep
        logical :: ok
        integer, parameter :: ns = 1800, nr = 160
        real(dp) :: t(ns), rmean(ns), a2(ns), o2(ns), eb, edges(0:nr), dens(nr, ns)
        real(dp) :: img(nr, ns - 1), tedge(ns), dmax
        complex(dp) :: an, ao
        integer :: i, u
        call kepler_init(kep, -1.0_dp, 1.0_dp, 1.0_dp, 1, 40.0_dp, 1.5_dp, 12, ok)
        if (.not. ok) error stop 'kepler setup failed'
        do i = 0, nr
            edges(i) = 1.6_dp*kep%rc*i/nr
        end do
        open (newunit=u, file='output/coulomb-kepler.csv', status='replace')
        write (u, '(a)') 't,t_over_TK,t_over_Trev,mean_r_over_rc,A2_numeric,A2_exact_levels,'// &
            'state_bound'
        dmax = 0
        do i = 1, ns
            t(i) = 1.1_dp*kep%revival*real(i - 1, dp)/real(ns - 1, dp)
            call kepler_mean_radius(kep, t(i), rmean(i), eb)
            call kepler_autocorrelation(kep, t(i), an, ao)
            a2(i) = abs(an)**2; o2(i) = abs(ao)**2; dmax = max(dmax, abs(an - ao))
            call kepler_radial_density(kep, t(i), edges, dens(:, i))
            write (u, '(es14.6,6(",",es14.6))') t(i), t(i)/kep%period, t(i)/kep%revival, &
                rmean(i)/kep%rc, a2(i), o2(i), eb
        end do
        close (u)
        open (newunit=u, file='output/coulomb-kepler-summary.csv', status='replace')
        write (u, '(a)') 'nbar,l,mean_n,rms_n,T_K,T_rev,rc,width,n_radial,kept,nbound,'// &
            'max_level_rel_error,continuum_weight,window_discarded2,e0,rho,max_A_minus_oracle,'// &
            'setup_s'
        write (u, '(f6.1,",",i0,16(",",es13.5))') kep%nbar, kep%l, kep%mean_n, kep%rms_n, &
            kep%period, kep%revival, kep%rc, kep%s, real(kep%mesh%n, dp), &
            real(size(kep%e), dp), real(kep%nbound, dp), kep%max_level_error, &
            kep%continuum_weight, kep%window_discarded2, kep%e0, kep%rho, dmax, &
            kep%setup_seconds
        close (u)
        print '(a,f8.3,a,es9.2)', 'Kepler: <n>=', kep%mean_n, ' max|A-A_exact|=', dmax
        kepler_series = reshape([t/kep%revival, rmean/kep%rc, a2, o2], [ns, 4])
        kepler_errors = [kep%max_level_error, dmax]
        call figure(figsize=[12.8_dp, 11.2_dp])
        call suptitle('Hydrogen radial Rydberg packet ($n_0 = 40$, $l = 1$, rms $\delta n = 1.5$):'// &
                      ' breathing, collapse and revivals')
        call subplot(3, 1, 1)
        call title('mean radius')
        call xlabel('$t / T_{rev}$,   $T_{rev} = (2 n_0/3) T_K$')
        call ylabel('mean radius / $r_{out}$')
        call xlim(0.0_dp, 1.1_dp); call ylim(0.3_dp, 1.3_dp)
        call plot(t/kep%revival, rmean/kep%rc, linestyle='-', color=[0.0_dp, 0.447_dp, 0.698_dp])
        call time_marks(kep%period/kep%revival, 0.3_dp, 1.3_dp, [0.3_dp, 0.3_dp, 0.3_dp], .true.)
        call legend(loc='upper center')
        call subplot(3, 1, 2)
        call title('autocorrelation $A(t)$ = overlap of $\psi(t)$ with $\psi(0)$: '// &
                   'represented spectrum versus exact levels')
        call xlabel('$t / T_{rev}$'); call ylabel('$|A(t)|^2$')
        call xlim(0.0_dp, 1.1_dp); call ylim(0.0_dp, 1.4_dp)
        call plot(t/kep%revival, a2, label='FEM-DVR spectrum', linestyle='-', &
                  color=[0.0_dp, 0.447_dp, 0.698_dp])
        call plot(t/kep%revival, o2, label='exact levels $E_n = -1/(2n^2)$', linestyle='--', &
                  color=[0.0_dp, 0.0_dp, 0.0_dp])
        call time_marks(kep%period/kep%revival, 0.0_dp, 1.4_dp, [0.3_dp, 0.3_dp, 0.3_dp], .false.)
        call legend(loc='upper center')
        call subplot(3, 1, 3)
        call title('radial probability density')
        call xlabel('$t / T_{rev}$'); call ylabel('$r / r_{out}$')
        tedge = t/kep%revival
        do i = 1, ns - 1
            img(:, i) = dens(:, i)*kep%rc
        end do
        call pcolormesh(tedge, edges/kep%rc, img, cmap='viridis', vmin=0.0_dp, &
                        vmax=maxval(img)*0.6_dp)
        call xlim(0.0_dp, 1.1_dp); call ylim(0.0_dp, 1.6_dp)
        call colorbar(label='$P(r,t) r_{out}$ (saturated above)')
        call time_marks(kep%period/kep%revival, 0.0_dp, 1.6_dp, [1.0_dp, 1.0_dp, 1.0_dp], .false.)
        call tight_layout()
        call savefig('output/coulomb-kepler.png')
    end subroutine kepler_case

    ! Vertical markers at T_K, T_rev/2 and T_rev (x in units of T_rev), told
    ! apart by line style; labelled through the legend when `labelled`.
    subroutine time_marks(tk, y0, y1, c, labelled)
        real(dp), intent(in) :: tk, y0, y1, c(3)
        logical, intent(in) :: labelled
        character(4), parameter :: ls(3) = [character(4) :: ':', '-.', '--']
        character(40), parameter :: lab(3) = [character(40) :: &
            'Kepler period $T_K$', 'half revival $T_{rev}/2$', 'full revival $T_{rev}$']
        real(dp) :: x
        integer :: i
        do i = 1, 3
            x = merge(tk, 0.5_dp*(i - 1), i == 1)
            if (labelled) then
                call plot([x, x], [y0, y1], linestyle=trim(ls(i)), color=c, label=trim(lab(i)))
            else
                call plot([x, x], [y0, y1], linestyle=trim(ls(i)), color=c)
            end if
        end do
    end subroutine time_marks

    ! Solve-only, check and derivative times with retained memory (warm, one run).
    subroutine cost_ledger()
        type(encounter_config_t) :: cfg
        type(encounter_t) :: enc
        complex(dp), allocatable :: a(:, :, :)
        real(dp) :: t, edges(0:nb), p(nb), eb, rad, pb, dpj, grad(3), img(96, 96), xe(0:96)
        real(dp) :: ts, tb, tc, tj, tv, tp
        integer(i64) :: c0, c1, rate
        integer :: u, i
        call setup_case(3.0_dp, 3.0_dp, 20.0_dp, 1.0_dp, cfg, enc, t)
        allocate (a(enc%mesh%n, 0:enc%mmax, 0:enc%lmax))
        do i = 0, nb
            edges(i) = pi*i/nb
        end do
        do i = 0, 96
            xe(i) = -enc%mesh%rmax*0.7_dp + 1.4_dp*enc%mesh%rmax*i/96
        end do
        call system_clock(c0, rate)
        do i = 1, 5
            call encounter_state(enc, t*i/5, a)
        end do
        call system_clock(c1); ts = real(c1 - c0, dp)/rate/5
        call system_clock(c0)
        call encounter_angular_bins(enc, a, SECTOR_PRODUCT, cfg%z0/4, edges, p)
        call system_clock(c1); tb = real(c1 - c0, dp)/rate
        call system_clock(c0)
        call encounter_state_error(enc, t, 1.0_dp, p(nb/2), eb, rad)
        call system_clock(c1); tc = real(c1 - c0, dp)/rate
        call system_clock(c0)
        call encounter_relative_projection(enc, a, SECTOR_PRODUCT, .false., xe, xe, 96, &
                                           enc%lambdabar/4, img)
        call system_clock(c1); tp = real(c1 - c0, dp)/rate
        call system_clock(c0)
        call encounter_bin_jvp(enc, t, SECTOR_PRODUCT, cfg%z0/4, pi/3, 2*pi/3, &
                               [1.0_dp, 0.0_dp, 0.0_dp], pb, dpj)
        call system_clock(c1); tj = real(c1 - c0, dp)/rate
        call system_clock(c0)
        call encounter_bin_vjp(enc, t, SECTOR_PRODUCT, cfg%z0/4, pi/3, 2*pi/3, pb, grad)
        call system_clock(c1); tv = real(c1 - c0, dp)/rate
        open (newunit=u, file='output/coulomb-cost-ledger.csv', status='replace')
        write (u, '(a)') 'item,seconds,bytes,note'
        write (u, '(a,es12.4,a,i0,a)') 'setup_total,', enc%ledger%setup_seconds, ',', &
            enc%ledger%stored_bytes, ',eta=3 sigma=3 lambdabar z0=20 lambdabar b=1 lambdabar'
        write (u, '(a,es12.4,a)') 'setup_projection,', enc%ledger%projection_seconds, &
            ',0,angular/radial projection incl. l/m estimate'
        write (u, '(a,es12.4,a)') 'setup_eigen_windows,', enc%ledger%eigen_seconds, &
            ',0,banded dsbevx per l plus residual norms'
        write (u, '(a,es12.4,a,i0,a)') 'solve_state_at_t,', ts, ',', &
            int(16*size(a), i64), ',exact evolution of retained coefficients'
        write (u, '(a,es12.4,a)') 'observable_36_bins,', tb, ',0,angular bins one sector'
        write (u, '(a,es12.4,a)') 'check_state_bound,', tc, ',0,interval e0+t rho/hbar and radius'
        write (u, '(a,es12.4,a)') 'projection_xz_96x96,', tp, ',0,visualization path'
        write (u, '(a,es12.4,a,i0,a)') 'derivative_jvp,', tj, ',', &
            int(32*enc%mesh%n**2, i64), ',full eigenbasis per l recomputed; O(n^2) scratch'
        write (u, '(a,es12.4,a,i0,a)') 'derivative_vjp_3_params,', tv, ',', &
            int(32*enc%mesh%n**2, i64), ',one backward pass gives d/db d/dsigma d/dg'
        close (u)
        bin_vjp = [pb, grad]
        cost_seconds = [ts, tb, tc, tp, tj, tv]
        print '(a,es9.2,a,es9.2,a,es9.2,a,es9.2)', 'cost: state ', ts, ' s, bins ', tb, &
            ' s, jvp ', tj, ' s, vjp ', tv
    end subroutine cost_ledger

    subroutine encounter_movies()
        allocate (mq(npx, npx, 3, nmove), mc(ncx, ncx, nmove))
        call encounter_movie(0.0_dp, 'headon')
        call encounter_movie(2.0_dp, 'offaxis')
    end subroutine encounter_movies

    ! Encounter movie in the relative coordinate r = x1 - x2 (eta = 3,
    ! sigma = 3 lambdabar = r0/2, z0 = 20 lambdabar, impact parameter b_l in
    ! lambdabar): slice y = 0 of the sector-normalized relative density for the
    ! singlet, triplet and spin-unresolved product beside the classical Wigner
    ! ensemble (slab |y| < h/2), log colour scale, slow motion around closest
    ! approach. Visualization path: interpolated synthesis, no pointwise bound.
    subroutine encounter_movie(b_l, tag)
        real(dp), intent(in) :: b_l
        character(*), intent(in) :: tag
        integer, parameter :: nscan = 81, nsb = 240, nbatch = 40, nbs = 20000
        real(dp), parameter :: slow = 4
        type(encounter_config_t) :: cfg
        type(encounter_t) :: enc
        type(animation_t) :: anim
        complex(dp), allocatable :: a(:, :, :), psi(:)
        real(dp), allocatable :: px(:), py(:), pz(:), rel(:, :, :), cm(:, :, :), tfine(:), &
            wcum(:)
        real(dp) :: tf, sedges(0:nsb), sp(nsb), rmean(nscan), ts(nscan), tau, tw, l, h, &
            ns, nt, np, w, x, z, drift, rho2, cell, dmax
        complex(dp) :: pp, pm
        integer :: f, i, j, ib, n, u, s
        logical :: ok
        call setup_case(3.0_dp, 3.0_dp, 20.0_dp, b_l, cfg, enc, tf)
        mv_lambdabar = enc%lambdabar; mv_speed = enc%speed; mv_z0 = cfg%z0
        mv_b = b_l; mv_r0 = enc%r0
        l = min(enc%mesh%rmax, 34*enc%lambdabar)
        h = mv_slab*enc%lambdabar
        do i = 0, npx
            xq(i) = -l + 2*l*i/npx
        end do
        do i = 0, ncx
            xc(i) = -l + 2*l*i/ncx
        end do
        ! Closest approach from the quantum mean separation.
        allocate (a(enc%mesh%n, 0:enc%mmax, 0:enc%lmax))
        do i = 0, nsb
            sedges(i) = enc%mesh%rmax*i/nsb
        end do
        do i = 1, nscan
            ts(i) = tf*(i - 1)/real(nscan - 1, dp)
            call encounter_state(enc, ts(i), a)
            call encounter_separation_density(enc, a, SECTOR_PRODUCT, sedges, sp)
            rmean(i) = sum((sedges(1:) + sedges(:nsb - 1))/2*sp)/sum(sp)
        end do
        mv_tc = ts(minloc(rmean, 1))
        ! Frame clock: playback rate 1/slow inside a Gaussian window of width tau
        ! centred after closest approach (reflection and outgoing interference).
        tau = 0.15_dp*tf; tw = mv_tc + 0.15_dp*tf
        allocate (tfine(0:4000), wcum(0:4000))
        wcum(0) = 0
        do i = 0, 4000
            tfine(i) = tf*i/4000.0_dp
            if (i > 0) wcum(i) = wcum(i - 1) + (tfine(i) - tfine(i - 1))* &
                (1 + (slow - 1)*exp(-((tfine(i) - tw)/tau)**2))
        end do
        j = 0
        do f = 1, nmove
            w = wcum(4000)*(f - 1)/real(nmove - 1, dp)
            do while (j < 3999 .and. wcum(j + 1) < w)
                j = j + 1
            end do
            x = (w - wcum(j))/max(wcum(j + 1) - wcum(j), tiny(1.0_dp))
            mtimes(f) = tfine(j) + min(1.0_dp, x)*(tfine(j + 1) - tfine(j))
            mslow(f) = 1 + (slow - 1)*exp(-((mtimes(f) - tw)/tau)**2)
        end do
        ! Quantum slice y = 0 on pixel centres; psi(-r) is the flipped image.
        allocate (px(npx*npx), py(npx*npx), pz(npx*npx), psi(npx*npx))
        n = 0
        do i = 1, npx
            do j = 1, npx
                n = n + 1
                px(n) = (xq(i - 1) + xq(i))/2; pz(n) = (xq(j - 1) + xq(j))/2
            end do
        end do
        py = 0
        open (newunit=u, file='output/coulomb-movie-'//tag//'.csv', status='replace')
        write (u, '(a)') 't,t_v_over_z0,slow_motion,norm_product,norm_singlet,'// &
            'norm_triplet,energy,wall_prob'
        do f = 1, nmove
            call encounter_state(enc, mtimes(f), a)
            call encounter_points(enc, a, px, py, pz, psi)
            np = encounter_norm(enc, a, SECTOR_PRODUCT)
            ns = encounter_norm(enc, a, SECTOR_SINGLET)
            nt = encounter_norm(enc, a, SECTOR_TRIPLET)
            do i = 1, npx
                do j = 1, npx
                    pp = psi((i - 1)*npx + j); pm = psi((npx - i)*npx + npx + 1 - j)
                    mq(j, i, 1, f) = abs(pp + pm)**2/(4*ns)
                    mq(j, i, 2, f) = abs(pp - pm)**2/(4*nt)
                    mq(j, i, 3, f) = (abs(pp)**2 + abs(pm)**2)/(2*np)
                end do
            end do
            mq(:, :, :, f) = mq(:, :, :, f)*enc%lambdabar**3
            write (u, '(es13.5,7(",",es13.5))') mtimes(f), mtimes(f)*enc%speed/cfg%z0, &
                mslow(f), np, ns, nt, encounter_energy(enc, a, SECTOR_PRODUCT), &
                encounter_wall_probability(enc, a, SECTOR_PRODUCT, 2*cfg%sigma)
        end do
        close (u)
        ! Classical ensemble, spin-unresolved (r and -r with weight 1/2), slab |y| < h/2.
        allocate (rel(3, nbs, nmove), cm(3, nbs, nmove))
        mc = 0
        cell = (xc(1) - xc(0))**2*h*nbatch*nbs/enc%lambdabar**3
        do ib = 1, nbatch
            call encounter_classical_ensemble(enc, nbs, movie_seed + ib, mtimes, rel, cm, &
                                              drift, ok)
            if (.not. ok) error stop 'classical movie ensemble failed'
            do f = 1, nmove
                do n = 1, nbs
                    if (abs(rel(2, n, f)) >= h/2) cycle
                    do s = -1, 1, 2
                        x = s*rel(1, n, f); z = s*rel(3, n, f)
                        i = floor((x - xc(0))/(xc(1) - xc(0))) + 1
                        j = floor((z - xc(0))/(xc(1) - xc(0))) + 1
                        if (i >= 1 .and. i <= ncx .and. j >= 1 .and. j <= ncx) &
                            mc(j, i, f) = mc(j, i, f) + 0.5_dp/cell
                    end do
                end do
            end do
        end do
        ! Fixed log colour range: four decades below the global quantum maximum.
        dmax = maxval(mq)
        mv_lmax = ceiling(2*log10(dmax))/2.0_dp; mv_lmin = mv_lmax - 4
        mq = log10(max(mq, 10.0_dp**(mv_lmin - 1)))
        mc = log10(max(mc, 10.0_dp**(mv_lmin - 1)))
        print '(3a,f7.3,a,es9.2,a,es9.2)', 'movie ', tag, ': t_c v/z0=', &
            mv_tc*enc%speed/cfg%z0, ' max density=', dmax, ' classical drift=', drift
        if (b_l == 0) call turning_cut(enc, cfg, a)
        call figure(figsize=[20.0_dp, 4.4_dp])
        canvas => get_global_figure()
        do i = 1, 4
            select case (i)
            case (1); f = 1
            case (2); f = minloc(abs(mtimes - mv_tc), 1)
            case (3); f = minloc(abs(mtimes - 0.72_dp*tf), 1)
            case default; f = nmove
            end select
            call draw(f)
            write (arg, '(4a,i0,a)') 'output/coulomb-encounter-', tag, '-frame', '-', f, '.png'
            call savefig(trim(arg))
            call canvas%clear()
        end do
        anim = animation_func(update, frames=nmove + nhold, interval=42, fig=canvas)
        call save_animation(anim, 'output/coulomb-encounter-'//tag//'.mp4', fps=24, &
                            status=status)
        if (status /= 0) print '(a,i0)', 'movie encoding status ', status
        call canvas%clear()
    end subroutine encounter_movie

    ! Head-on cuts at closest approach t_c along the collision axis z.
    ! (a) x = 0, density averaged over the cylinder rho < 1/k: distinguishable
    ! quantum pair versus the classical ensemble (tunnelling tail inside r0).
    ! (b) x = 2 r0, y = 0: singlet/triplet exchange fringes (spacing -> lambda/2)
    ! versus the spin-unresolved product and the classical ensemble (rod of
    ! cross-section (1/k)^2). Classical bands are +-1 binomial standard error.
    subroutine turning_cut(enc, cfg, a)
        use fortnum_quadrature, only: gauss_legendre_ab
        type(encounter_t), intent(in) :: enc
        type(encounter_config_t), intent(in) :: cfg
        complex(dp), intent(inout) :: a(:, 0:, 0:)
        integer, parameter :: nz = 600, nzc = 96, nq = 6, nbt = 10, nbc = 40000
        real(dp) :: zq(nz), dq(nz), zc(0:nzc), zm(nzc), cnt(nzc, 2), dc(nzc, 2), ec(nzc, 2), &
            sq(nq), wq(nq), xx(nz), yy(nz), lo, hi, drift, lb, vol(2), ytop, ylo, xo, &
            ds(nz), dt(nz), du(nz), peak(nz), ymark, nrm(3), r(3)
        real(dp), allocatable :: rel(:, :, :), cm(:, :, :)
        complex(dp) :: psi(nz), psim(nz)
        integer :: i, q, ib, n, np, u, s
        logical :: ok
        lb = enc%lambdabar; xo = 2*enc%r0
        lo = -24*lb; hi = 6*lb
        call encounter_state(enc, mv_tc, a)
        nrm = [encounter_norm(enc, a, SECTOR_PRODUCT), encounter_norm(enc, a, SECTOR_SINGLET), &
               encounter_norm(enc, a, SECTOR_TRIPLET)]
        call gauss_legendre_ab(nq, 0.0_dp, 1.0_dp, sq, wq)
        do i = 1, nz
            zq(i) = lo + (hi - lo)*(i - 0.5_dp)/nz
        end do
        dq = 0; yy = 0
        do q = 1, nq
            xx = lb*sqrt(sq(q))
            call encounter_points(enc, a, xx, yy, zq, psi)
            dq = dq + wq(q)*abs(psi)**2
        end do
        dq = dq/nrm(1)*lb**3
        xx = xo
        call encounter_points(enc, a, xx, yy, zq, psi)
        call encounter_points(enc, a, -xx, yy, -zq, psim)
        ds = abs(psi + psim)**2/(4*nrm(2))*lb**3
        dt = abs(psi - psim)**2/(4*nrm(3))*lb**3
        du = (abs(psi)**2 + abs(psim)**2)/(2*nrm(1))*lb**3
        do i = 0, nzc
            zc(i) = lo + (hi - lo)*i/nzc
        end do
        zm = (zc(1:) + zc(:nzc - 1))/2
        allocate (rel(3, nbc, 1), cm(3, nbc, 1))
        cnt = 0
        do ib = 1, nbt
            call encounter_classical_ensemble(enc, nbc, movie_seed + 100 + ib, [mv_tc], rel, &
                                              cm, drift, ok)
            if (.not. ok) error stop 'classical cut ensemble failed'
            do n = 1, nbc
                r = rel(:, n, 1)
                i = floor((r(3) - lo)/(zc(1) - zc(0))) + 1
                if (r(1)**2 + r(2)**2 < lb**2 .and. i >= 1 .and. i <= nzc) &
                    cnt(i, 1) = cnt(i, 1) + 1
                do s = -1, 1, 2     ! spin-unresolved: r and -r with weight 1/2
                    i = floor((s*r(3) - lo)/(zc(1) - zc(0))) + 1
                    if (abs(s*r(1) - xo) < lb/2 .and. abs(r(2)) < lb/2 .and. i >= 1 .and. &
                        i <= nzc) cnt(i, 2) = cnt(i, 2) + 0.5_dp
                end do
            end do
        end do
        vol = [acos(-1.0_dp), 1.0_dp]*lb**2*(zc(1) - zc(0))*nbt*nbc/lb**3
        do s = 1, 2
            dc(:, s) = cnt(:, s)/vol(s); ec(:, s) = sqrt(cnt(:, s)*s)/s/vol(s)
        end do
        open (newunit=u, file='output/coulomb-turning-cut.csv', status='replace')
        write (u, '(a)') 'kz,axis_quantum_product,offaxis_singlet,offaxis_triplet,'// &
            'offaxis_spin_unresolved   (densities times lambdabar^3)'
        do i = 1, nz
            write (u, '(es13.5,4(",",es13.5))') zq(i)/lb, dq(i), ds(i), dt(i), du(i)
        end do
        write (u, '(a)') 'kz_centre,axis_classical,axis_se,offaxis_classical,offaxis_se'
        do i = 1, nzc
            write (u, '(es13.5,4(",",es13.5))') zm(i)/lb, dc(i, 1), ec(i, 1), dc(i, 2), ec(i, 2)
        end do
        close (u)
        call figure(figsize=[12.8_dp, 6.4_dp])
        write (arg, '(a,f5.2,a)') 'Head-on encounter at closest approach, $t v/z_0$ = ', &
            mv_tc*enc%speed/cfg%z0, ',  $\eta = 3$, $\sigma = r_0/2$'
        call suptitle(trim(arg))
        call subplot(1, 2, 1)
        call title('on the axis: no fringes, tunnelling into $|z| < r_0$')
        call xlabel('$k z$ (relative coordinate, $x = y = 0$)')
        call ylabel('density $\times k^{-3}$ (cylinder $\rho < 1/k$)')
        call set_yscale('log')
        ytop = 10*max(maxval(dq), maxval(dc(:, 1) + ec(:, 1))); ylo = ytop*1.0e-5_dp
        call xlim(lo/lb, hi/lb); call ylim(ylo, ytop)
        call classical_curve(zm/lb, dc(:, 1), ec(:, 1), cnt(:, 1), ylo)
        call plot(pack(zq/lb, dq > ylo), pack(dq, dq > ylo), label='quantum, distinguishable pair', &
                  linestyle='-', color=[0.0_dp, 0.447_dp, 0.698_dp])
        call plot([-enc%r0/lb, -enc%r0/lb], [ylo, ytop], linestyle=':', &
                  color=[0.0_dp, 0.0_dp, 0.0_dp], label='classical turning point $-r_0$')
        call legend(loc='upper left')
        call subplot(1, 2, 2)
        call title('beside it ($k x = 2 k r_0 = 12$): exchange fringes')
        call xlabel('$k z$ (relative coordinate, $y = 0$)')
        call ylabel('density $\times k^{-3}$')
        call set_yscale('log')
        ytop = 6*max(maxval(ds), maxval(dt)); ylo = ytop*1.0e-3_dp
        call xlim(lo/lb, hi/lb); call ylim(ylo, ytop)
        call classical_curve(zm/lb, dc(:, 2), ec(:, 2), cnt(:, 2), ylo)
        call plot(pack(zq/lb, du > ylo), pack(du, du > ylo), label='spin-unresolved product', &
                  linestyle='-', &
                  color=[0.0_dp, 0.0_dp, 0.0_dp])
        call plot(pack(zq/lb, du > ylo), max(pack(ds, du > ylo), ylo), label='singlet', &
                  linestyle='-', &
                  color=[0.0_dp, 0.447_dp, 0.698_dp])
        call plot(pack(zq/lb, du > ylo), max(pack(dt, du > ylo), ylo), label='triplet', &
                  linestyle='--', &
                  color=[0.835_dp, 0.369_dp, 0.0_dp])
        ! lambda/2 bracket between the second and third singlet maxima (z < 0).
        np = 0
        do i = nz - 1, 2, -1
            if (zq(i) < 0 .and. ds(i) > ds(i - 1) .and. ds(i) >= ds(i + 1) .and. &
                ds(i) > 1.0e-3_dp*maxval(ds)) then
                np = np + 1; peak(np) = zq(i)/lb
            end if
        end do
        if (np >= 3) then
            ymark = 1.4_dp*maxval(ds)
            call plot([peak(2) - acos(-1.0_dp), peak(2)], [ymark, ymark], linestyle='-', &
                      color=[0.0_dp, 0.0_dp, 0.0_dp], label='$\lambda/2 = \pi/k$ (asymptotic)')
            do i = 0, 1
                call plot([peak(2) - i*acos(-1.0_dp), peak(2) - i*acos(-1.0_dp)], &
                          [ymark/1.25_dp, ymark*1.25_dp], linestyle='-', &
                          color=[0.0_dp, 0.0_dp, 0.0_dp])
            end do
            print '(a,5f8.3)', 'turning cut: singlet maxima kz =', peak(1:min(5, np))
        end if
        call legend(loc='upper left')
        call tight_layout()
        call savefig('output/coulomb-turning-cut.png')
    end subroutine turning_cut

    subroutine classical_curve(z, d, e, cnt, ylo)
        real(dp), intent(in) :: z(:), d(:), e(:), cnt(:), ylo
        call fill_between(z, max(d - e, ylo), max(d + e, ylo), &
                          color=[0.0_dp, 0.62_dp, 0.45_dp], alpha=0.25_dp)
        call plot(pack(z, d > ylo), pack(d, d > ylo), &
                  label='classical ensemble ($\pm 1$ s.e.)', linestyle=':', &
                  color=[0.0_dp, 0.62_dp, 0.45_dp])
    end subroutine classical_curve


    subroutine update(number)
        integer, intent(in) :: number
        call canvas%clear()
        call draw(min(number, nmove))
        call canvas%set_rendered(.false.)
    end subroutine update

    subroutine draw(f)
        integer, intent(in) :: f
        integer :: p, nt
        real(dp) :: lb, tk(8)
        character(16) :: tl(8)
        character(48), parameter :: names(4) = [character(48) :: 'singlet (even $l$)', &
            'triplet (odd $l$)', 'spin-unresolved product', 'classical ensemble']
        lb = mv_lambdabar
        nt = 0
        do p = ceiling(mv_lmin), floor(mv_lmax)
            nt = nt + 1; tk(nt) = p
            write (tl(nt), '(a,i0,a)') '$10^{', p, '}$'
        end do
        if (mslow(f) > 1.05_dp) then
            write (arg, '(a,f3.1,a,f5.2,a,f4.1,a)') 'Relative density in the plane $y = 0$, '// &
                '$\eta = 3$, $\sigma = r_0/2 = 3/k$, $kb = ', mv_b*1.0_dp, '$:   $t v/z_0$ = ', &
                mtimes(f)*mv_speed/mv_z0, '   (slow motion ', mslow(f), '$\times$)'
        else
            write (arg, '(a,f3.1,a,f5.2)') 'Relative density in the plane $y = 0$, '// &
                '$\eta = 3$, $\sigma = r_0/2 = 3/k$, $kb = ', mv_b*1.0_dp, '$:   $t v/z_0$ = ', &
                mtimes(f)*mv_speed/mv_z0
        end if
        call suptitle(trim(arg))
        do p = 1, 4
            call subplot(1, 4, p)
            call title(trim(names(p)))
            call xlabel('$k x$')
            if (p == 1) call ylabel('$k z$')
            if (p < 4) then
                call pcolormesh(xq/lb, xq/lb, mq(:, :, p, f), cmap='viridis', &
                                vmin=mv_lmin, vmax=mv_lmax)
            else
                call pcolormesh(xc/lb, xc/lb, mc(:, :, f), cmap='viridis', &
                                vmin=mv_lmin, vmax=mv_lmax)
            end if
            call xlim(xq(0)/lb, xq(npx)/lb); call ylim(xq(0)/lb, xq(npx)/lb)
            if (p == 4) then
                call colorbar(label='$|ψ(x,0,z)|^2 k^{-3}$', ticks=tk(:nt), ticklabels=tl(:nt))
            else
                call colorbar(ticks=tk(:nt), ticklabels=tl(:nt))
            end if
        end do
        call tight_layout()
    end subroutine draw
end program coulomb_encounter_demo
