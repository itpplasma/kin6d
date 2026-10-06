!> QG6 evidence: quantum vs classical Coulomb logarithm as a function of the
!> Sommerfeld parameter (Yukawa partial waves, cut models, Born, Mott), test
!> electron slowing-down rate vs temperature (all-classical / hybrid Monte Carlo,
!> kernel quadratures, Landau/Spitzer with classical or quantum logarithm), cost
!> at matched statistical accuracy, and a beam-relaxation movie with promoted
!> encounters highlighted. Atomic units (hbar = m_e = e^2/(4 pi eps0) = 1).
!> Usage: hybrid_encounter_demo [--config FILE] (namelist /hybrid_encounter/).
program hybrid_encounter_demo
    use, intrinsic :: iso_fortran_env, only: dp => real64, i64 => int64
    use fortnum_rng, only: rng_t, rng_seed, rng_normal
    use fortnum_status, only: fortnum_status_t
    use fortplot, only: figure, get_global_figure, figure_t, subplot, plot, scatter, &
        errorbar, title, suptitle, xlabel, ylabel, xlim, ylim, legend, tight_layout, &
        savefig, fill_between
    use fortplot_animation, only: animation_t, animation_func => FuncAnimation, &
        save_animation
    use kin6d_species, only: atomic_units, electron
    use kin6d_quantum_collision_kernel
    use kin6d_hybrid_encounter, single_encounter => hybrid_encounter
    use kin6d_result_io, only: run_t, result_file_t
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp), parameter :: blue(3) = [0.000_dp, 0.447_dp, 0.698_dp]
    real(dp), parameter :: verm(3) = [0.835_dp, 0.369_dp, 0.000_dp]
    real(dp), parameter :: green(3) = [0.000_dp, 0.620_dp, 0.451_dp]
    real(dp), parameter :: orange(3) = [0.902_dp, 0.624_dp, 0.000_dp]
    real(dp), parameter :: purple(3) = [0.800_dp, 0.475_dp, 0.655_dp]
    real(dp), parameter :: sky(3) = [0.337_dp, 0.706_dp, 0.914_dp]
    real(dp), parameter :: black(3) = [0.0_dp, 0.0_dp, 0.0_dp]
    real(dp), parameter :: grey(3) = [0.6_dp, 0.6_dp, 0.6_dp]
    ! Configuration (namelist /hybrid_encounter/).
    real(dp) :: kappa_log = 100, density = 2.0e-12_dp, tmin = 1.0e-3_dp, tmax = 100
    real(dp) :: eta_switch = 1, zone = 4, step_eta = 0.02_dp, beam_x = 2
    real(dp) :: movie_temperature = 1, movie_density = 2.5e-4_dp
    integer :: neta = 25, ntemp = 11, nsample = 20000, movie_particles = 300
    integer :: movie_frames = 200
    integer(i64) :: seed = 20261006_i64
    logical :: movie = .true.
    namelist /hybrid_encounter/ kappa_log, density, tmin, tmax, eta_switch, zone, &
        step_eta, beam_x, movie_temperature, movie_density, neta, ntemp, nsample, &
        movie_particles, movie_frames, seed, movie
    type(run_t) :: run
    type(result_file_t) :: result
    type(hybrid_config_t) :: hcfg
    character(256) :: msg
    integer :: ios, cu, status
    ! Section 1 results.
    real(dp), allocatable :: etas(:), ly(:), dly(:), lyc(:), lcq(:), lcc(:), lborn(:)
    real(dp), allocatable :: lmott(:, :), visr(:), vise(:), ycost(:), yfit(:), ytail(:)
    ! Section 2 results.
    real(dp), allocatable :: temps(:), etau(:), leff(:, :), lse(:, :), lspitzer(:, :)
    real(dp), allocatable :: cost(:, :), defect(:), nprom(:), edef(:), angch(:), lkern(:, :)
    real(dp) :: fd_check(3)
    ! Section 3 (movie) state shared with the frame callback.
    real(dp), allocatable :: vel(:, :, :, :), tf(:)
    logical, allocatable :: flag(:, :, :)
    real(dp) :: nu0, vth, lims(4)
    type(figure_t), pointer :: canvas => null()

    call run%start('hybrid_encounter_demo')
    if (run%argument_count() > 0) error stop 'usage: hybrid_encounter_demo [--config FILE]'
    ios = 0; msg = ''
    call run%open_config(cu)
    if (run%has_config()) read (cu, nml=hybrid_encounter, iostat=ios, iomsg=msg)
    call run%close_config(cu, 'hybrid_encounter', ios, msg)
    if (.not. (kappa_log > 1 .and. density > 0 .and. tmax > tmin .and. tmin > 0 .and. &
               neta >= 3 .and. ntemp >= 2 .and. nsample >= 100 .and. zone > 0 .and. &
               step_eta > 0 .and. beam_x > 0 .and. movie_particles > 0 .and. &
               movie_frames >= 2)) error stop 'hybrid_encounter: invalid configuration'
    call run%open_resolved(cu)
    write (cu, nml=hybrid_encounter)
    call run%close_resolved(cu)
    call execute_command_line('mkdir -p output', exitstat=status)
    hcfg%eta_switch = eta_switch; hcfg%zone = zone; hcfg%step_eta = step_eta

    call coulomb_logarithm_scan()
    call rate_scan()
    call write_results()
    call coulomb_log_figure()
    call rate_figure()
    if (movie) call beam_movie()
contains
    ! 1. Normalized momentum-transfer logarithm vs eta at kappa = k lambda_D.
    subroutine coulomb_logarithm_scan()
        type(yukawa_plan_t) :: plan
        type(yukawa_config_t) :: ycfg
        real(dp) :: d, lp, lm, t0, t1
        integer :: i, u, nrings
        logical :: ok
        allocate (etas(neta), ly(neta), dly(neta), lyc(neta), lcq(neta), lcc(neta), &
                  lborn(neta), lmott(neta, 4), visr(neta), vise(neta), ycost(neta), &
                  yfit(neta), ytail(neta))
        nrings = floor(kappa_log)
        open (newunit=u, file='output/hybrid-coulomb-log.csv', status='replace')
        write (u, '(a,f0.3)') '# QG6 lnL = sigma_tr/(4 pi b90^2) vs eta, kappa = k lambda_D = ', &
            kappa_log
        write (u, '(a)') '# eta,quantum_yukawa,dlnL_deta,classical_yukawa,quantum_cut,'// &
            'classical_cut,born_yukawa,mott_unpolarized,mott_singlet,mott_triplet,'// &
            'identical_classical,visc_ring,visc_exact,yukawa_cpu_s,fit_residual,tail_estimate'
        do i = 1, neta
            etas(i) = 10.0_dp**(-2 + 3.5_dp*(i - 1)/real(neta - 1, dp))
            call cpu_time(t0)
            call yukawa_plan(plan, etas(i), kappa_log, ycfg, ok)
            if (.not. ok) error stop 'yukawa plan failed'
            call yukawa_transfer_log(plan, etas(i), ly(i), dly(i), ok)
            if (.not. ok) error stop 'yukawa solve failed'
            call cpu_time(t1); ycost(i) = t1 - t0
            yfit(i) = plan%ledger%fit_residual; ytail(i) = plan%ledger%tail_estimate
            lyc(i) = classical_yukawa_log(kappa_log/etas(i), 200, 16)
            call quantum_cut_log(etas(i), kappa_log, lcq(i), d)
            call classical_cut_log(etas(i), kappa_log, lcc(i), d)
            lborn(i) = born_yukawa_log(kappa_log)
            lmott(i, 1) = mott_transfer_log(etas(i), kappa_log, SPIN_UNPOLARIZED)
            lmott(i, 2) = mott_transfer_log(etas(i), kappa_log, SPIN_SINGLET)
            lmott(i, 3) = mott_transfer_log(etas(i), kappa_log, SPIN_TRIPLET)
            lmott(i, 4) = mott_transfer_log(etas(i), kappa_log, SPIN_IDENTICAL_CLASSICAL)
            visr(i) = coulomb_ring_viscosity(etas(i), nrings)
            vise(i) = coulomb_viscosity_exact(etas(i), nrings)
            write (u, '(es13.6,15(",",es15.7e3))') etas(i), ly(i), dly(i), lyc(i), lcq(i), &
                lcc(i), lborn(i), lmott(i, :), visr(i), vise(i), ycost(i), yfit(i), ytail(i)
        end do
        close (u)
        ! Equation-level derivative against central differences on a frozen plan.
        call yukawa_plan(plan, 0.7_dp, kappa_log, ycfg, ok)
        call yukawa_transfer_log(plan, 0.7_dp, fd_check(1), fd_check(2), ok)
        call yukawa_transfer_log(plan, 0.7_dp + 1.0e-5_dp, lp, d, ok)
        call yukawa_transfer_log(plan, 0.7_dp - 1.0e-5_dp, lm, d, ok)
        fd_check(3) = (lp - lm)/2.0e-5_dp
        print '(a,f8.4,a,es10.3)', 'Yukawa lnL(eta=0.7) = ', fd_check(1), &
            ', implicit vs FD d/deta rel. diff = ', abs(fd_check(2)/fd_check(3) - 1)
    end subroutine coulomb_logarithm_scan

    subroutine mean_se(xs, mean, se)
        real(dp), intent(in) :: xs(:)
        real(dp), intent(out) :: mean, se
        mean = sum(xs)/size(xs)
        se = sqrt(sum((xs - mean)**2)/(size(xs) - 1)/size(xs))
    end subroutine mean_se

    ! 2. Slowing-down rate of a test electron at v_b = beam_x sqrt(2T) in a
    ! Maxwellian electron background at density n, as an effective logarithm
    ! L_eff = nu/nu_Landau(lnL = 1). Paired samples (same seed) for the
    ! classical and hybrid Monte Carlo.
    subroutine rate_scan()
        type(hybrid_model_t) :: model
        real(dp), allocatable :: xc(:), xh(:)
        real(dp) :: t, lam, vb, uu, nu1, mc, se, mh, sh, md, sd, ed, ac
        integer(i64) :: c0, c1, rate
        integer :: i, np, u
        logical :: ok
        allocate (temps(ntemp), etau(ntemp), leff(ntemp, 3), lse(ntemp, 3), &
                  lspitzer(ntemp, 2), cost(ntemp, 2), defect(ntemp), nprom(ntemp), &
                  edef(ntemp), angch(ntemp), lkern(ntemp, 3), xc(nsample), xh(nsample))
        open (newunit=u, file='output/hybrid-relaxation-rate.csv', status='replace')
        write (u, '(a,es10.3,a,f5.2,a,i0)') '# QG6 test-electron slowing rate, n = ', density, &
            ' a.u., v_b = ', beam_x, ' sqrt(2T), samples = ', nsample
        write (u, '(a)') '# T,eta_u,Leff_classical_mc,se,Leff_hybrid_mc,se,paired_diff,se,'// &
            'Leff_kernel_classical,Leff_kernel_hybrid,Leff_kernel_quantum,ln(lD/b90),'// &
            'ln(lD/lambdabar),cost1pct_classical_s,cost1pct_hybrid_s,interface_defect,'// &
            'promoted_fraction,max_energy_defect,max_angular_change_hbar'
        do i = 1, ntemp
            t = tmin*(tmax/tmin)**((i - 1)/real(ntemp - 1, dp))
            temps(i) = t
            lam = sqrt(t/(4*pi*density)); vb = beam_x*sqrt(2*t); uu = sqrt(vb**2 + 3*t)
            hcfg%promote = .false.
            call hybrid_init(model, atomic_units(), electron(), electron(), lam, hcfg, ok)
            if (.not. ok) error stop 'hybrid_init failed'
            etau(i) = model%g/(model%hbar*uu)
            nu1 = kernel_slowing_rate(model, vb, t, density, -1, 1.0_dp)
            call system_clock(c0, rate)
            call hybrid_rate_samples(model, vb, t, density, seed + i, xc, np, ed, ac, ok)
            if (.not. ok) error stop 'classical samples failed'
            call system_clock(c1)
            call mean_se(xc, mc, se)
            cost(i, 1) = real(c1 - c0, dp)/rate/nsample*(se*sqrt(real(nsample, dp))/mc)**2/1.0e-4_dp
            model%cfg%promote = .true.
            call system_clock(c0)
            call hybrid_rate_samples(model, vb, t, density, seed + i, xh, np, ed, ac, ok)
            if (.not. ok) error stop 'hybrid samples failed'
            call system_clock(c1)
            call mean_se(xh, mh, sh)
            cost(i, 2) = real(c1 - c0, dp)/rate/nsample*(sh*sqrt(real(nsample, dp))/mh)**2/1.0e-4_dp
            call mean_se(xh - xc, md, sd)
            nprom(i) = real(np, dp)/nsample; edef(i) = ed; angch(i) = ac
            leff(i, :) = [mc, mh, md]/nu1; lse(i, :) = [se, sh, sd]/nu1
            lkern(i, 1) = kernel_slowing_rate(model, vb, t, density, MODEL_CLASSICAL)/nu1
            lkern(i, 2) = kernel_slowing_rate(model, vb, t, density, MODEL_HYBRID)/nu1
            lkern(i, 3) = kernel_slowing_rate(model, vb, t, density, MODEL_QUANTUM)/nu1
            lspitzer(i, 1) = log(lam*model%mu*uu**2/model%g)
            lspitzer(i, 2) = log(lam*model%mu*uu/model%hbar)
            defect(i) = lkern(i, 2)/lkern(i, 3) - 1
            write (u, '(es13.6,18(",",es15.7e3))') t, etau(i), leff(i, 1), lse(i, 1), &
                leff(i, 2), lse(i, 2), leff(i, 3), lse(i, 3), lkern(i, :), lspitzer(i, :), &
                cost(i, :), defect(i), nprom(i), edef(i), angch(i)
            print '(a,es9.2,a,f7.3,a,3f8.3,a,3f8.3,a,f6.3)', 'T=', t, ' eta_u=', etau(i), &
                ' MC c/h/diff', leff(i, :), ' kernel c/h/q', lkern(i, :), ' prom', nprom(i)
        end do
        close (u)
    end subroutine rate_scan

    subroutine write_results()
        call result%create(run)
        call result%write('eta', etas, '1', 'Sommerfeld parameter of the lnL scan', dims='eta')
        call result%write('kappa', kappa_log, '1', 'k lambda_D of the lnL scan')
        call result%write('lnl_quantum_yukawa', ly, '1', &
            'sigma_tr/(4 pi b90^2), Yukawa partial waves on the radial DVR', dims='eta', &
            status='candidate; Born-tail estimate and fit residual indicators')
        call result%write('dlnl_deta_quantum_yukawa', dly, '1', &
            'implicit equation-level derivative d lnL/d eta', dims='eta', status='candidate')
        call result%write('lnl_classical_yukawa', lyc, '1', 'classical Yukawa orbit integral', &
            dims='eta', status='quadrature candidate')
        call result%write('lnl_quantum_cut', lcq, '1', 'pure Coulomb, partial-wave cut', &
            dims='eta', status='closed form')
        call result%write('lnl_classical_cut', lcc, '1', 'Rutherford, impact cut', dims='eta', &
            status='closed form')
        call result%write('lnl_born_yukawa', lborn, '1', 'first Born Yukawa', dims='eta')
        call result%write('lnl_mott', lmott, '1', &
            'tagged identical-particle transfer: unpolarized, singlet, triplet, classical', &
            dims='eta,spin', status='quadrature candidate')
        call result%write('viscosity_ring', visr, '1', 'ring-kernel sigma^(2)/(4 pi b90^2)', &
            dims='eta')
        call result%write('viscosity_exact', vise, '1', 'partial-wave sigma^(2)/(4 pi b90^2)', &
            dims='eta')
        call result%write('yukawa_fd_check', fd_check, '1', &
            'lnL, implicit d/deta, central FD at eta=0.7', dims='fd')
        call result%write('temperature', temps, 'E_h', 'background temperature', dims='T')
        call result%write('eta_u', etau, '1', 'eta at u = sqrt(v_b^2 + 3T)', dims='T')
        call result%write('leff_mc', leff, '1', &
            'nu/nu_Landau(lnL=1): classical MC, hybrid MC, paired hybrid-classical', &
            dims='T,mc_method', status='statistical (CLT), standard error in leff_mc_se')
        call result%write('leff_mc_se', lse, '1', 'standard errors of leff_mc', &
            dims='T,mc_method', status='statistical')
        call result%write('leff_kernel', lkern, '1', &
            'kernel quadrature: classical cut, hybrid, all-quantum ring cut', &
            dims='T,kernel_model', &
            status='quadrature candidate')
        call result%write('leff_spitzer', lspitzer, '1', &
            'ln(lambda_D/b90(u)), ln(lambda_D/lambdabar(u))', dims='T,landau_log', &
            status='analytic approximation')
        call result%write('cost_1pct', cost, 's', &
            'CPU seconds to 1% relative standard error: classical, hybrid', &
            dims='T,cost_method')
        call result%write('interface_defect', defect, '1', &
            'hybrid/all-quantum kernel rate - 1', dims='T', status='model/interface')
        call result%write('promoted_fraction', nprom, '1', 'promoted encounter fraction', &
            dims='T')
        call result%write('max_energy_defect', edef, '1', 'max relative KDK energy defect', &
            dims='T', status='diagnostic')
        call result%write('max_angular_change', angch, 'hbar', &
            'max angular-momentum change at promotion', dims='T', status='diagnostic')
        call result%close()
    end subroutine write_results

    subroutine coulomb_log_figure()
        real(dp) :: le(neta), ref(neta), x(2)
        complex(dp) :: p, dpsi
        integer :: i
        character(64) :: heading
        le = log10(etas)
        call figure(figsize=[13.2_dp, 8.8_dp])
        call subplot(2, 2, 1)
        write (heading, '(a,i0)') 'Coulomb logarithm, $k\lambda_D$ = ', nint(kappa_log)
        call title(trim(heading))
        call plot(le, lyc, color=verm, label='classical Yukawa orbit')
        call plot(le, lcc, color=verm, linestyle='--', label='classical, impact cut')
        call plot(le, ly, color=blue, marker='o', markersize=4.0_dp, &
                  label='quantum Yukawa (partial waves)')
        call plot(le, lcq, color=blue, linestyle='--', label='quantum, partial-wave cut')
        call plot(le, lborn, color=black, linestyle=':', label='first Born Yukawa')
        call plot(le, log(kappa_log/etas), color=black, linestyle='--', &
                  label='$\ln(\lambda_D/b_{90})$')
        x = [le(1), le(neta)]
        call plot(x, [1.0_dp, 1.0_dp]*log(kappa_log), color=black, linestyle='-.', &
                  label='$\ln(k\lambda_D)$')
        call xlim(le(1), le(neta)); call ylim(0.0_dp, 1.15_dp*maxval(log(kappa_log/etas)))
        call xlabel('$\log_{10}\eta$'); call ylabel('$\sigma_{tr}/(4\pi b_{90}^2)$')
        call legend(loc='upper right')
        call subplot(2, 2, 2)
        call title('Quantum minus classical logarithm')
        do i = 1, neta
            call complex_digamma(cmplx(1.0_dp, etas(i), dp), p, dpsi)
            ref(i) = log(etas(i)) - real(p)
        end do
        call plot(le, ly - lyc, color=blue, marker='o', markersize=4.0_dp, label='Yukawa')
        call plot(le, lcq - lcc, color=verm, label='cut models')
        call plot(le, ref, color=black, linestyle='--', &
                  label='$\ln\eta$ - Re $\psi(1+i\eta)$')
        call xlim(le(1), le(neta))
        call xlabel('$\log_{10}\eta$'); call ylabel('$\Delta \ln\Lambda$')
        call legend(loc='lower right')
        call subplot(2, 2, 3)
        call title('Exchange (Mott) correction, identical electrons')
        call plot(le, lmott(:, 1) - lmott(:, 4), color=blue, label='unpolarized')
        call plot(le, lmott(:, 2) - lmott(:, 4), color=verm, linestyle='--', label='singlet')
        call plot(le, lmott(:, 3) - lmott(:, 4), color=green, linestyle='-.', label='triplet')
        call plot(x, [1.0_dp, 1.0_dp]*(-log(2.0_dp)/2), color=black, linestyle=':', &
                  label='$-\ln 2/2$')
        call xlim(le(1), le(neta))
        call xlabel('$\log_{10}\eta$'); call ylabel('$\Delta \ln\Lambda$ vs classical identical')
        call legend(loc='upper right')
        call subplot(2, 2, 4)
        call title('Ring kernel: second-moment defect')
        call plot(le, visr/vise - 1, color=blue, label='$\sigma^{(2)}$ ring / partial waves - 1')
        call plot(x, [0.0_dp, 0.0_dp], color=black, linestyle='--', &
                  label='$\sigma_{tr}$: exact by construction')
        call xlim(le(1), le(neta))
        call xlabel('$\log_{10}\eta$'); call ylabel('relative defect')
        call legend(loc='upper right')
        call suptitle('QG6 binary Coulomb kernel; candidates with quadrature/tail indicators')
        call tight_layout()
        call savefig('output/hybrid-coulomb-log.png')
    end subroutine coulomb_log_figure

    subroutine rate_figure()
        real(dp) :: lt(ntemp), x(2), rc(ntemp), rh(ntemp)
        lt = log10(temps); x = [lt(1), lt(ntemp)]
        call figure(figsize=[13.2_dp, 8.8_dp])
        call subplot(2, 2, 1)
        call title('Slowing-down rate as effective logarithm')
        call plot(lt, lkern(:, 1), color=verm, label='classical kernel')
        call plot(lt, lkern(:, 2), color=blue, label='hybrid kernel')
        call plot(lt, lkern(:, 3), color=green, linestyle=':', label='all-quantum kernel')
        call errorbar(lt, leff(:, 1), yerr=2*lse(:, 1), fmt='s', color=verm, &
                      capsize=3.0_dp, label='classical MC, 95% CLT')
        call errorbar(lt, leff(:, 2), yerr=2*lse(:, 2), fmt='o', color=blue, &
                      capsize=3.0_dp, label='hybrid MC, 95% CLT')
        call plot(lt, lspitzer(:, 1), color=black, linestyle='--', &
                  label='Landau, $\ln(\lambda_D/b_{90,u})$')
        call plot(lt, lspitzer(:, 2), color=black, linestyle='-.', &
                  label='Landau, $\ln(k_u\lambda_D)$')
        call xlim(x(1), x(2))
        call xlabel('$\log_{10}(T/E_h)$'); call ylabel('$\nu_s/\nu_{Landau}(\ln\Lambda=1)$')
        call legend(loc='upper left')
        call subplot(2, 2, 2)
        call title('Monte Carlo minus kernel (statistical)')
        rc = leff(:, 1)/lkern(:, 1) - 1; rh = leff(:, 2)/lkern(:, 2) - 1
        call fill_between(lt, -2*lse(:, 1)/lkern(:, 1), 2*lse(:, 1)/lkern(:, 1), &
                          color=verm, alpha=0.25_dp)
        call fill_between(lt, -2*lse(:, 2)/lkern(:, 2), 2*lse(:, 2)/lkern(:, 2), &
                          color=blue, alpha=0.25_dp)
        call plot(lt, rc, color=verm, marker='s', markersize=5.0_dp, &
                  label='classical MC/kernel - 1')
        call plot(lt, rh, color=blue, marker='o', markersize=5.0_dp, &
                  label='hybrid MC/kernel - 1')
        call plot(x, [0.0_dp, 0.0_dp], color=black, linestyle='--', &
                  label='kernel; bands 2 s.e. (paired, overlapping)')
        call xlim(x(1), x(2))
        call xlabel('$\log_{10}(T/E_h)$'); call ylabel('relative difference')
        call legend(loc='upper left')
        call subplot(2, 2, 3)
        call title('Model and interface error vs all-quantum kernel')
        call plot(lt, log10(abs(lkern(:, 1)/lkern(:, 3) - 1) + 1.0e-16_dp), color=verm, &
                  marker='s', markersize=5.0_dp, label='classical (model error)')
        call plot(lt, log10(abs(defect) + 1.0e-16_dp), color=blue, marker='o', &
                  markersize=5.0_dp, label='hybrid (interface defect)')
        call xlim(x(1), x(2))
        call xlabel('$\log_{10}(T/E_h)$'); call ylabel('$\log_{10}$ |rate/quantum - 1|')
        call legend(loc='lower right')
        call subplot(2, 2, 4)
        call title('CPU time to 1% statistical error')
        call plot(lt, log10(cost(:, 1)), color=verm, marker='s', markersize=5.0_dp, &
                  label='all-classical')
        call plot(lt, log10(cost(:, 2)), color=blue, marker='o', markersize=5.0_dp, &
                  label='hybrid')
        call xlim(x(1), x(2))
        call xlabel('$\log_{10}(T/E_h)$'); call ylabel('$\log_{10}$(CPU s)')
        call legend(loc='upper left')
        call suptitle('QG6 test electron in a Maxwellian electron background, cut screening')
        call tight_layout()
        call savefig('output/hybrid-relaxation-rate.png')
    end subroutine rate_figure

    ! 3. Beam relaxation movie: the same test electrons and random stream,
    ! all-classical (left) and hybrid (right); markers flag frames in which an
    ! electron had an encounter inside the near zone (classical: treated
    ! classically; hybrid: promoted to the quantum ring kernel).
    subroutine beam_movie()
        type(hybrid_model_t) :: mc, mh
        type(rng_t) :: rc, rh
        type(fortnum_status_t) :: st
        type(animation_t) :: anim
        real(dp) :: lam, vb, dt
        integer :: f, k, u, nenc(2), nall
        logical :: ok
        character(64) :: fname
        lam = sqrt(movie_temperature/(4*pi*movie_density))
        vb = beam_x*sqrt(2*movie_temperature); vth = sqrt(2*movie_temperature)
        hcfg%promote = .false.
        call hybrid_init(mc, atomic_units(), electron(), electron(), lam, hcfg, ok)
        hcfg%promote = .true.
        call hybrid_init(mh, atomic_units(), electron(), electron(), lam, hcfg, ok)
        if (.not. ok) error stop 'movie model failed'
        nu0 = kernel_slowing_rate(mh, vb, movie_temperature, movie_density, MODEL_HYBRID)
        dt = 1.5_dp/nu0/(movie_frames - 1)
        allocate (vel(3, movie_particles, movie_frames, 2), flag(movie_particles, &
                  movie_frames, 2), tf(movie_frames))
        call rng_seed(rc, seed + 1000, st); call rng_seed(rh, seed + 1000, st)
        vel = 0; vel(3, :, 1, :) = vb; flag = .false.; nenc = 0; tf(1) = 0
        do f = 2, movie_frames
            vel(:, :, f, 1) = vel(:, :, f - 1, 1); vel(:, :, f, 2) = vel(:, :, f - 1, 2)
            call hybrid_beam_step(mc, vel(:, :, f, 1), movie_temperature, movie_density, &
                                  dt, rc, flag(:, f, 1), nenc(1), ok)
            if (.not. ok) error stop 'classical beam step failed'
            call hybrid_beam_step(mh, vel(:, :, f, 2), movie_temperature, movie_density, &
                                  dt, rh, flag(:, f, 2), nenc(2), ok)
            if (.not. ok) error stop 'hybrid beam step failed'
            tf(f) = tf(f - 1) + dt
        end do
        open (newunit=u, file='output/hybrid-beam.csv', status='replace')
        write (u, '(a,es10.3,a,es10.3,a,i0,a,2i9)') '# beam movie: T = ', movie_temperature, &
            ', n = ', movie_density, ', N = ', movie_particles, ', encounters c/h ', nenc
        write (u, '(a)') '# t*nu0,mean_vz_classical/vb,mean_vz_hybrid/vb,near_classical,promoted'
        do f = 1, movie_frames
            write (u, '(f10.5,2(",",es15.7),2(",",i6))') tf(f)*nu0, &
                sum(vel(3, :, f, 1))/movie_particles/vb, sum(vel(3, :, f, 2))/movie_particles/vb, &
                count(flag(:, f, 1)), count(flag(:, f, 2))
        end do
        close (u)
        print '(a,2i9,a,f7.3,a,f7.3)', 'beam movie encounters c/h ', nenc, &
            '; final <v_z>/v_b classical ', sum(vel(3, :, movie_frames, 1))/movie_particles/vb, &
            ' hybrid ', sum(vel(3, :, movie_frames, 2))/movie_particles/vb
        lims(1) = minval(vel(3, :, :, :))/vth; lims(2) = maxval(vel(3, :, :, :))/vth
        lims(3) = 0; lims(4) = maxval(sqrt(vel(1, :, :, :)**2 + vel(2, :, :, :)**2))/vth
        lims(1:2) = lims(1:2) + [-0.05_dp, 0.05_dp]*(lims(2) - lims(1))
        lims(4) = 1.05_dp*lims(4)
        nall = movie_frames + 24
        call figure(figsize=[12.8_dp, 6.4_dp])
        canvas => get_global_figure()
        do k = 1, 3
            f = max(1, (k - 1)*movie_frames/2)
            call draw(f)
            write (fname, '(a,i3.3,a)') 'output/hybrid-beam-frame', f, '.png'
            call savefig(trim(fname))
        end do
        anim = animation_func(draw, frames=nall, interval=42, fig=canvas)
        call save_animation(anim, 'output/hybrid-beam.mp4', fps=24, status=status)
        if (status /= 0) print '(a,i0)', 'movie encoding status ', status
    end subroutine beam_movie

        subroutine draw(frame)
            integer, intent(in) :: frame
            integer :: ff, side, nf
            real(dp), allocatable :: xs(:), ys(:)
            character(64) :: heading
            ff = min(frame, movie_frames)
            call canvas%clear()
            do side = 1, 2
                call subplot(1, 2, side)
                xs = vel(3, :, ff, side)/vth
                ys = sqrt(vel(1, :, ff, side)**2 + vel(2, :, ff, side)**2)/vth
                nf = count(flag(:, ff, side))
                if (nf > 0) then
                    call scatter(pack(xs, flag(:, ff, side)), pack(ys, flag(:, ff, side)), &
                        color=black, marker='s', s=30.0_dp, &
                        label=merge('near-zone encounter (classical)', &
                                    'promoted (quantum ring kernel) ', side == 1))
                end if
                call scatter(xs, ys, color=merge(verm, blue, side == 1), marker='o', &
                             s=9.0_dp, label='test electrons')
                write (heading, '(a,f5.2)') trim(merge('all-classical,', 'hybrid,       ', &
                    side == 1))//' $t\nu_0$ = ', tf(ff)*nu0
                call title(trim(heading))
                call xlim(lims(1), lims(2)); call ylim(lims(3), lims(4))
                call xlabel('$v_z/v_{th}$'); call ylabel('transverse speed / $v_{th}$')
                call legend(loc='upper left')
            end do
            call suptitle('QG6 beam relaxation, same electrons and random stream')
            call tight_layout()
        end subroutine draw
end program hybrid_encounter_demo
