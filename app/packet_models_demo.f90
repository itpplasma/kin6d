! Q2 evidence: reduced wave-packet models of the Q1 two-electron Coulomb
! encounter (classical Wigner ensemble M1, product WPMD M2, symmetrized and
! antisymmetrized WPMD M2S/M2T, relative TDVP Gaussian M3 and M3S/M3T,
! coherent Gaussian sums M4 with K = 2, 4, 8) against the Q1 partial-wave
! reference at the same physical parameters, observables and bins, with the
! full-Hamiltonian Duhamel estimator. Writes CSV, figures and a movie to
! output/ and output/packet_models_demo.nc.
! Usage: packet_models_demo [--config FILE]   (namelist &packet_models)
program packet_models_demo
    use, intrinsic :: iso_fortran_env, only: dp => real64, i64 => int64
    use kin6d_coulomb_encounter
    use kin6d_packet_models
    use kin6d_coulomb_encounter_leaves, only: ce_observable_radius
    use kin6d_result_io, only: run_t, result_file_t
    use fortplot, only: figure, get_global_figure, figure_t, subplot, plot, title, &
        xlabel, ylabel, xlim, ylim, legend, savefig, pcolormesh, tight_layout, &
        set_xscale, set_yscale, suptitle, colorbar
    use fortplot_animation, only: animation_t, animation_func => FuncAnimation, &
        save_animation
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp)
    integer, parameter :: nb = 36, nmod = 9, nsep = 60
    ! Wave-packet models in table order; M1 (classical) is column nmod + 1.
    integer, parameter :: models(nmod) = [MODEL_M2, MODEL_M2S, MODEL_M2T, MODEL_M3, &
                                          MODEL_M3S, MODEL_M3T, MODEL_M4, MODEL_M4, MODEL_M4]
    integer, parameter :: kvals(nmod) = [1, 1, 1, 1, 1, 1, 2, 4, 8]
    ! Configuration (&packet_models).
    integer :: samples = 61, nclass = 20000, movie_frames = 192
    real(dp) :: etas(4) = [0.3_dp, 1.0_dp, 3.0_dp, 10.0_dp]
    real(dp) :: b_values(4) = [0.0_dp, 5.0_dp, 10.0_dp, 20.0_dp]
    logical :: movie = .true., b_scan = .true.
    namelist /packet_models/ samples, nclass, movie_frames, etas, b_values, movie, b_scan
    type(run_t) :: run
    type(result_file_t) :: result
    ! Results kept for the figures and the NetCDF file.
    real(dp) :: err(nmod + 1, 4), bnd(nmod + 1, 4), tv(nmod + 1, 4), terr(nmod + 1, 4)
    real(dp) :: tbnd(nmod + 1, 4), eff(nmod + 1, 4), secs(nmod + 1, 4)
    real(dp) :: bins_ref(nb, 3, 4), bins_mod(nb, nmod + 1, 4), theta(nb)
    real(dp) :: berr(3, 4), bbnd(3, 4), btv(3, 4)
    ! Movie state shared with the animation callback (host association).
    integer, parameter :: npx = 112
    real(dp), allocatable :: img(:, :, :, :), mtime(:)
    real(dp) :: xe(0:npx), mv_lambdabar, mv_speed, mv_z0, vlo, vhi
    type(figure_t), pointer :: canvas => null()
    integer :: unit, ios, status, i
    character(256) :: msg
    character(200) :: arg
    call run%start('packet_models_demo')
    ios = 0; msg = ''
    call run%open_config(unit)
    if (run%has_config()) read (unit, nml=packet_models, iostat=ios, iomsg=msg)
    call run%close_config(unit, 'packet_models', ios, msg)
    call run%open_resolved(unit)
    write (unit, nml=packet_models)
    call run%close_resolved(unit)
    if (samples < 3 .or. nclass < 100 .or. movie_frames < 2) error stop &
        'packet_models: need samples >= 3, nclass >= 100, movie_frames >= 2'
    call execute_command_line('mkdir -p output', exitstat=status)
    if (status /= 0) error stop 'cannot create output/'
    do i = 1, nb
        theta(i) = (i - 0.5_dp)*180.0_dp/nb
    end do
    call eta_scan()
    call eta_figures()
    if (b_scan) call impact_scan()
    if (movie) call density_movie()
    call result%create(run)
    call result%attribute('estimator_scope', 'Duhamel bound in the Q1 represented '// &
        'space (binary64, trapezoid time integral of sampled ||r_h||); continuum '// &
        'adds Q1 discretization error and model projection defect (indicators)')
    call result%write('eta', etas, '1', 'Sommerfeld parameter', dims='eta')
    call result%write('theta_deg', theta, 'deg', 'CM-frame bin centre', dims='bin')
    call result%write('state_error', err, '1', 'model-reference distance at t_final '// &
        '(relative state, sector-normalized)', dims='model,eta')
    call result%write('state_bound', bnd, '1', 'Duhamel full-H estimator at t_final', &
        dims='model,eta')
    call result%write('total_error', terr, '1', 'two-body (CM x relative) distance', &
        dims='model,eta')
    call result%write('total_bound', tbnd, '1', 'two-body Duhamel estimator', &
        dims='model,eta')
    call result%write('bin_tv_error', tv, '1', 'sum_b |P_model - P_ref| (36 bins)', &
        dims='model,eta')
    call result%write('efficiency_index', eff, '1', 'estimate/actual state error', &
        dims='model,eta')
    call result%write('model_seconds', secs, 's', 'model integration time', &
        dims='model,eta')
    call result%write('bins_reference', bins_ref, '1', 'Q1 bins (product, singlet, '// &
        'triplet)', dims='bin,sector,eta')
    call result%write('bins_model', bins_mod, '1', 'model bins (own sector)', &
        dims='bin,model,eta')
    call result%close()
contains
    ! Q1 case for one eta in reduced units hbar = m_e = e^2/(4 pi eps0) = 1,
    ! identical to coulomb_encounter_demo (sigma = max(10, 4 eta) lambdabar,
    ! z0 = 6 sigma, t_final = 2 z0/v).
    subroutine setup_case(eta, b_l, enc, t)
        real(dp), intent(in) :: eta, b_l
        type(encounter_t), intent(out) :: enc
        real(dp), intent(out) :: t
        type(encounter_config_t) :: cfg
        character(:), allocatable :: m
        logical :: ok
        real(dp) :: sig_l
        sig_l = max(10.0_dp, 4*eta)
        cfg%g = 1; cfg%m1 = 1; cfg%m2 = 1; cfg%k = 0.5_dp/eta
        cfg%sigma = sig_l/cfg%k; cfg%z0 = 6*cfg%sigma; cfg%b = b_l/cfg%k
        t = 2*cfg%z0/(cfg%hbar*cfg%k/0.5_dp); cfg%t_final = t
        call encounter_init(enc, cfg, ok, m)
        if (.not. ok) error stop 'encounter setup failed: '//m
    end subroutine setup_case

    subroutine run_model(enc, spec, t, tr, seconds, nfev, ok)
        type(encounter_t), intent(in) :: enc
        type(packet_spec_t), intent(in) :: spec
        real(dp), intent(in) :: t
        type(packet_trace_t), intent(out) :: tr
        real(dp), intent(out) :: seconds
        integer(i64), intent(out) :: nfev
        logical, intent(out) :: ok
        type(packet_sample_t), allocatable :: smp(:)
        character(:), allocatable :: m
        integer(i64) :: c0, c1, rate
        integer :: j
        call system_clock(c0, rate)
        call packet_samples(enc, spec, [(t*(j - 1)/real(samples - 1, dp), j=1, samples)], &
                            smp, ok, m, nfev)
        call system_clock(c1)
        if (.not. ok) then
            print '(a)', model_name(spec)//': '//m
            return
        end if
        seconds = real(c1 - c0, dp)/real(rate, dp)
        call packet_trace(enc, spec, smp, tr)
    end subroutine run_model

    subroutine classical_bins(enc, t, rcut, edges, pc)
        type(encounter_t), intent(in) :: enc
        real(dp), intent(in) :: t, rcut, edges(0:)
        real(dp), intent(out) :: pc(:)
        real(dp), allocatable :: rel(:, :, :), cm(:, :, :)
        real(dp) :: drift, th
        integer :: n, ib
        logical :: ok
        allocate (rel(3, nclass, 1), cm(3, nclass, 1))
        call encounter_classical_ensemble(enc, nclass, 1234_i64, [t], rel, cm, drift, ok)
        if (.not. ok) error stop 'classical ensemble failed'
        pc = 0
        do n = 1, nclass
            if (norm2(rel(:, n, 1)) < rcut) cycle
            th = acos(rel(3, n, 1)/norm2(rel(:, n, 1)))
            ib = min(size(pc), int(th/(edges(1) - edges(0))) + 1)
            pc(ib) = pc(ib) + 1.0_dp/nclass
        end do
    end subroutine classical_bins

    pure function bin_edges() result(edges)
        real(dp) :: edges(0:nb)
        integer :: ib
        edges = [(pi*ib/nb, ib=0, nb)]
    end function bin_edges

    ! Sommerfeld-eta scan at b = 0: every model against the Q1 reference.
    subroutine eta_scan()
        type(encounter_t) :: enc
        type(packet_spec_t) :: spec
        type(packet_trace_t) :: tr
        complex(dp), allocatable :: a(:, :, :)
        real(dp) :: t, edges(0:nb), pm(nb), rad, sd
        integer(i64) :: nfev
        integer :: ie, im, ue, ub, uc, ib, sec, sectors(3)
        logical :: ok
        sectors = [SECTOR_PRODUCT, SECTOR_SINGLET, SECTOR_TRIPLET]
        edges = bin_edges()
        open (newunit=ue, file='output/packet-eta-errors.csv', status='replace')
        write (ue, '(a)') 'eta,model,sector,state_error,phase_free_error,state_bound,'// &
            'efficiency,total_error,total_bound,bin_tv_error,bin_max_error,'// &
            'bin_radius_from_bound,energy_drift,norm_drift,projection_defect,'// &
            'max_overlap_condition,max_regularized,status'
        open (newunit=ub, file='output/packet-eta-bins.csv', status='replace')
        write (ub, '(a)') 'eta,theta_lo_deg,theta_hi_deg,ref_product,ref_singlet,'// &
            'ref_triplet,M2,M2S,M2T,M3,M3S,M3T,M4K2,M4K4,M4K8,M1_classical'
        open (newunit=uc, file='output/packet-cost-ledger.csv', status='replace')
        write (uc, '(a)') 'eta,model,real_parameters,model_seconds,rhs_evaluations,'// &
            'diagnostic_seconds,q1_setup_seconds,q1_kept_states,q1_stored_bytes'
        do ie = 1, size(etas)
            call setup_case(etas(ie), 0.0_dp, enc, t)
            allocate (a(enc%mesh%n, 0:enc%mmax, 0:enc%lmax))
            call encounter_state(enc, t, a)
            do sec = 1, 3
                call encounter_angular_bins(enc, a, sectors(sec), enc%cfg%z0/4, edges, &
                                            bins_ref(:, sec, ie))
            end do
            deallocate (a)
            do im = 1, nmod
                spec = packet_spec_t()
                spec%model = models(im); spec%k = kvals(im)
                call run_model(enc, spec, t, tr, secs(im, ie), nfev, ok)
                if (.not. ok) error stop 'eta scan model refused'
                sec = model_sector(spec)
                call encounter_angular_bins(enc, tr%a, sec, enc%cfg%z0/4, edges, pm)
                bins_mod(:, im, ie) = pm
                sec = merge(1, merge(2, 3, sec == SECTOR_SINGLET), sec == SECTOR_PRODUCT)
                err(im, ie) = tr%error(samples); bnd(im, ie) = tr%bound(samples)
                terr(im, ie) = tr%total_error(samples); tbnd(im, ie) = tr%total_bound(samples)
                tv(im, ie) = sum(abs(pm - bins_ref(:, sec, ie)))
                eff(im, ie) = bnd(im, ie)/max(err(im, ie), tiny(1.0_dp))
                call ce_observable_radius(bnd(im, ie), rad)
                write (ue, '(f6.2,",",a,",",i0,13(",",es12.4),",",i0,",",a)') etas(ie), &
                    model_name(spec), model_sector(spec), err(im, ie), &
                    tr%phase_error(samples), bnd(im, ie), eff(im, ie), terr(im, ie), &
                    tbnd(im, ie), tv(im, ie), maxval(abs(pm - bins_ref(:, sec, ie))), rad, &
                    maxval(abs(tr%energy - tr%energy(1)))/abs(tr%energy(1)), &
                    maxval(abs(tr%norm - 1)), 1 - tr%projected_norm(samples)/tr%norm(samples), &
                    tr%max_overlap_condition, tr%max_regularized, tr%status
                write (uc, '(f6.2,",",a,",",i0,",",es12.4,",",i0,2(",",es12.4),2(",",i0))') &
                    etas(ie), model_name(spec), nparam(spec), secs(im, ie), nfev, &
                    tr%seconds_projection, enc%ledger%setup_seconds, enc%ledger%kept_states, &
                    enc%ledger%stored_bytes
                print '(a,f5.1,1x,a6,a,es10.3,a,es10.3,a,f7.2,a,es10.3,a,f8.2,a)', 'eta=', &
                    etas(ie), model_name(spec), ' e=', err(im, ie), ' bound=', bnd(im, ie), &
                    ' I=', eff(im, ie), ' TV=', tv(im, ie), ' t=', secs(im, ie) + &
                    tr%seconds_projection, ' s'
            end do
            ! M1: classical Wigner ensemble against the product reference.
            call classical_bins(enc, t, enc%cfg%z0/4, edges, pm)
            bins_mod(:, nmod + 1, ie) = pm
            tv(nmod + 1, ie) = sum(abs(pm - bins_ref(:, 1, ie)))
            err(nmod + 1, ie) = -1; bnd(nmod + 1, ie) = -1; terr(nmod + 1, ie) = -1
            tbnd(nmod + 1, ie) = -1; eff(nmod + 1, ie) = -1; secs(nmod + 1, ie) = 0
            sd = sum(sqrt(pm*(1 - pm)/nclass))
            write (ue, '(f6.2,",M1,0,n/a,n/a,n/a,n/a,n/a,n/a,",es12.4,",",es12.4,'// &
                '",n/a,n/a,n/a,n/a,n/a,n/a,binomial sum of sd ",es10.3)') etas(ie), &
                tv(nmod + 1, ie), maxval(abs(pm - bins_ref(:, 1, ie))), sd
            print '(a,f5.1,a,es10.3,a,es10.3)', 'eta=', etas(ie), ' M1     TV=', &
                tv(nmod + 1, ie), ' (sum sd ', sd
            do ib = 1, nb
                write (ub, '(f6.2,2(",",f8.3),16(",",es13.5))') etas(ie), edges(ib - 1)*180/pi, &
                    edges(ib)*180/pi, bins_ref(ib, :, ie), bins_mod(ib, :, ie)
            end do
        end do
        close (ue); close (ub); close (uc)
    end subroutine eta_scan

    pure integer function nparam(spec)
        type(packet_spec_t), intent(in) :: spec
        select case (spec%model)
        case (MODEL_M2); nparam = 12
        case (MODEL_M3); nparam = 6
        case (MODEL_M4); nparam = 8*spec%k
        case default; nparam = 8
        end select
    end function nparam

    function okabe(i) result(c)
        integer, intent(in) :: i
        real(dp) :: c(3)
        real(dp), parameter :: pal(3, 7) = reshape([0.000_dp, 0.447_dp, 0.698_dp, &
            0.835_dp, 0.369_dp, 0.000_dp, 0.000_dp, 0.620_dp, 0.451_dp, 0.902_dp, &
            0.624_dp, 0.000_dp, 0.800_dp, 0.475_dp, 0.655_dp, 0.337_dp, 0.706_dp, &
            0.914_dp, 0.0_dp, 0.0_dp, 0.0_dp], [3, 7])
        c = pal(:, i)
    end function okabe

    ! (i) state error and estimator vs eta; bin total-variation error vs eta.
    ! (ii) angular distributions at eta = 1 and 10, models vs reference.
    subroutine eta_figures()
        integer, parameter :: shown(6) = [1, 2, 4, 7, 8, 9]   ! M2 M2S M3 M4K2 M4K4 M4K8
        integer, parameter :: col(6) = [4, 2, 1, 6, 5, 3]
        integer :: j, ie, ip, idx(2)
        real(dp) :: ymin
        call figure(figsize=[12.8_dp, 6.4_dp])
        call subplot(1, 2, 1)
        call title('state error at $t = 2 z_0/v$ (b = 0)')
        call xlabel('Sommerfeld $\eta$'); call ylabel('$\|\psi_{model} - \psi_{Q1}\|$')
        do j = 1, size(shown)
            call plot(etas, err(shown(j), :), label=model_label(shown(j))//' actual', &
                      linestyle='-', marker='o', color=okabe(col(j)))
            call plot(etas, bnd(shown(j), :), label=model_label(shown(j))//' estimator', &
                      linestyle='--', marker='x', color=okabe(col(j)))
        end do
        call set_xscale('log'); call set_yscale('log')
        call legend(loc='lower left')
        call subplot(1, 2, 2)
        call title('angular-bin error $\sum_b |P_b - P_b^{Q1}|$')
        call xlabel('Sommerfeld $\eta$'); call ylabel('total-variation distance (36 bins)')
        do j = 1, size(shown)
            call plot(etas, max(tv(shown(j), :), 1.0e-6_dp), label=model_label(shown(j)), &
                      linestyle='-', marker='o', color=okabe(col(j)))
        end do
        call plot(etas, max(tv(nmod + 1, :), 1.0e-6_dp), label='M1 classical', &
                  linestyle=':', marker='s', color=okabe(7))
        call set_xscale('log'); call set_yscale('log')
        call legend(loc='lower left')
        call tight_layout()
        call savefig('output/packet-error-vs-eta.png')
        idx = [minloc(abs(etas - 1), 1), minloc(abs(etas - 10), 1)]
        call figure(figsize=[13.2_dp, 8.8_dp])
        do ip = 1, 2
            ie = idx(ip)
            ymin = 1.0e-7_dp
            call subplot(2, 2, 2*ip - 1)
            write (arg, '(a,f4.1,a)') 'distinguishable (product), $\eta$ = ', etas(ie)
            call title(trim(arg))
            call xlabel('$\theta$ (deg), CM frame, $t = 2 z_0/v$')
            call ylabel('bin probability (5 deg, $|r| > z_0/4$)')
            call plot(theta, max(bins_mod(:, nmod + 1, ie), ymin), label='M1 classical', &
                      linestyle=':', color=okabe(5))
            call plot(theta, max(bins_mod(:, 1, ie), ymin), label='M2 WPMD product', &
                      linestyle='-', color=okabe(4))
            call plot(theta, max(bins_mod(:, 4, ie), ymin), label='M3 one Gaussian', &
                      linestyle='-', color=okabe(1))
            call plot(theta, max(bins_mod(:, 9, ie), ymin), label='M4 K = 8', &
                      linestyle='-', color=okabe(3))
            call plot(theta, max(bins_ref(:, 1, ie), ymin), label='Q1 reference', &
                      linestyle='--', color=okabe(7))
            call set_yscale('log'); call ylim(ymin, 1.5_dp)
            call legend(loc='upper right')
            call subplot(2, 2, 2*ip)
            write (arg, '(a,f4.1,a)') 'singlet / triplet, $\eta$ = ', etas(ie)
            call title(trim(arg))
            call xlabel('$\theta$ (deg), CM frame, $t = 2 z_0/v$')
            call ylabel('sector-normalized bin probability')
            call plot(theta, max(bins_mod(:, 2, ie), ymin), label='M2S symmetrized WPMD', &
                      linestyle='-', color=okabe(2))
            call plot(theta, max(bins_mod(:, 3, ie), ymin), label='M2T antisymmetrized WPMD', &
                      linestyle='-', color=okabe(6))
            call plot(theta, max(bins_ref(:, 2, ie), ymin), label='Q1 singlet', &
                      linestyle='--', color=okabe(7))
            call plot(theta, max(bins_ref(:, 3, ie), ymin), label='Q1 triplet', &
                      linestyle=':', color=okabe(7))
            call set_yscale('log'); call ylim(ymin, 1.5_dp)
            call legend(loc='upper right')
        end do
        call suptitle('Q2 packet models vs the Q1 partial-wave reference (b = 0)')
        call tight_layout()
        call savefig('output/packet-angular-eta1-eta10.png')
    end subroutine eta_figures

    function model_label(im) result(s)
        integer, intent(in) :: im
        character(:), allocatable :: s
        type(packet_spec_t) :: spec
        spec%model = models(im); spec%k = kvals(im)
        s = model_name(spec)
    end function model_label

    ! Impact-parameter scan at eta = 1 (product sector): M2, M3, M4 K = 8.
    subroutine impact_scan()
        integer, parameter :: bm(3) = [MODEL_M2, MODEL_M3, MODEL_M4], bk(3) = [1, 1, 8]
        type(encounter_t) :: enc
        type(packet_spec_t) :: spec
        type(packet_trace_t) :: tr
        complex(dp), allocatable :: a(:, :, :)
        real(dp) :: t, edges(0:nb), pr(nb), pm(nb), sx
        integer(i64) :: nfev
        integer :: ib, j, u
        logical :: ok
        edges = bin_edges()
        open (newunit=u, file='output/packet-b-scan.csv', status='replace')
        write (u, '(a)') 'b_over_lambdabar,model,state_error,state_bound,efficiency,'// &
            'bin_tv_error,model_seconds,max_overlap_condition'
        do ib = 1, size(b_values)
            call setup_case(1.0_dp, b_values(ib), enc, t)
            allocate (a(enc%mesh%n, 0:enc%mmax, 0:enc%lmax))
            call encounter_state(enc, t, a)
            call encounter_angular_bins(enc, a, SECTOR_PRODUCT, enc%cfg%z0/4, edges, pr)
            deallocate (a)
            do j = 1, 3
                spec = packet_spec_t(); spec%model = bm(j); spec%k = bk(j)
                call run_model(enc, spec, t, tr, sx, nfev, ok)
                if (.not. ok) then
                    write (u, '(f6.2,",",a,a)') b_values(ib), model_name(spec), &
                        ',refused,refused,refused,refused,refused,refused'
                    berr(j, ib) = -1; bbnd(j, ib) = -1; btv(j, ib) = -1
                    cycle
                end if
                call encounter_angular_bins(enc, tr%a, SECTOR_PRODUCT, enc%cfg%z0/4, edges, pm)
                berr(j, ib) = tr%error(samples); bbnd(j, ib) = tr%bound(samples)
                btv(j, ib) = sum(abs(pm - pr))
                write (u, '(f6.2,",",a,6(",",es12.4))') b_values(ib), model_name(spec), &
                    berr(j, ib), bbnd(j, ib), bbnd(j, ib)/berr(j, ib), btv(j, ib), sx, &
                    tr%max_overlap_condition
                print '(a,f5.1,1x,a6,a,es10.3,a,es10.3,a,es10.3)', 'b/lambdabar=', &
                    b_values(ib), model_name(spec), ' e=', berr(j, ib), ' bound=', &
                    bbnd(j, ib), ' TV=', btv(j, ib)
            end do
        end do
        close (u)
    end subroutine impact_scan

    ! (iii) relative (x,z) density, y-integrated, at eta = 1: exact Q1 vs one
    ! Gaussian (M3) vs K = 8 sum, same axes, clock and log colour range.
    subroutine density_movie()
        type(encounter_t) :: enc
        type(packet_spec_t) :: spec
        type(packet_sample_t), allocatable :: s3(:), s8(:)
        complex(dp), allocatable :: a(:, :, :), ad(:, :, :)
        real(dp) :: t, lr, sedges(0:nsep), sep(nsep, 3)
        character(:), allocatable :: m
        integer(i64) :: nfev
        integer :: f, ix, u
        logical :: ok
        type(animation_t) :: anim
        call setup_case(1.0_dp, 0.0_dp, enc, t)
        mv_lambdabar = enc%lambdabar; mv_speed = enc%speed; mv_z0 = enc%cfg%z0
        lr = min(enc%mesh%rmax, 1.6_dp*enc%cfg%z0)
        do ix = 0, npx
            xe(ix) = -lr + 2*lr*ix/npx
        end do
        sedges = [(lr*ix/nsep, ix=0, nsep)]
        allocate (mtime(movie_frames), img(npx, npx, 3, movie_frames))
        mtime = [(t*(f - 1)/real(movie_frames - 1, dp), f=1, movie_frames)]
        spec%model = MODEL_M3
        call packet_samples(enc, spec, mtime, s3, ok, m, nfev)
        if (.not. ok) error stop 'movie M3: '//m
        spec%model = MODEL_M4; spec%k = 8
        call packet_samples(enc, spec, mtime, s8, ok, m, nfev)
        if (.not. ok) error stop 'movie M4: '//m
        allocate (a(enc%mesh%n, 0:enc%mmax, 0:enc%lmax))
        open (newunit=u, file='output/packet-movie-separation.csv', status='replace')
        write (u, '(a)') 't,t_v_over_z0,r_over_lambdabar,P_q1,P_m3,P_m4k8'
        do f = 1, movie_frames
            call encounter_state(enc, mtime(f), a)
            call encounter_relative_projection(enc, a, SECTOR_PRODUCT, .false., xe, xe, npx, &
                                               enc%lambdabar/4, img(:, :, 1, f))
            call encounter_separation_density(enc, a, SECTOR_PRODUCT, sedges, sep(:, 1))
            call packet_project(enc, s3(f)%g, s3(f)%rate, a, ad)
            call encounter_relative_projection(enc, a, SECTOR_PRODUCT, .false., xe, xe, npx, &
                                               enc%lambdabar/4, img(:, :, 2, f))
            call encounter_separation_density(enc, a, SECTOR_PRODUCT, sedges, sep(:, 2))
            call packet_project(enc, s8(f)%g, s8(f)%rate, a, ad)
            call encounter_relative_projection(enc, a, SECTOR_PRODUCT, .false., xe, xe, npx, &
                                               enc%lambdabar/4, img(:, :, 3, f))
            call encounter_separation_density(enc, a, SECTOR_PRODUCT, sedges, sep(:, 3))
            do ix = 1, nsep
                write (u, '(es13.5,5(",",es13.5))') mtime(f), mtime(f)*enc%speed/enc%cfg%z0, &
                    (sedges(ix - 1) + sedges(ix))/(2*enc%lambdabar), sep(ix, :)
            end do
        end do
        close (u)
        vhi = log10(maxval(img)); vlo = vhi - 4
        img = log10(max(img, 10**vlo))
        call figure(figsize=[19.2_dp, 6.4_dp])
        canvas => get_global_figure()
        call draw(movie_frames)
        call savefig('output/packet-density-final.png')
        call canvas%clear()
        call draw(movie_frames/2)
        call savefig('output/packet-density-mid.png')
        anim = animation_func(update, frames=movie_frames + 24, interval=42, fig=canvas)
        call save_animation(anim, 'output/packet-density.mp4', fps=24, status=status)
        if (status /= 0) print '(a,i0)', 'movie encoding status ', status
    end subroutine density_movie

    subroutine update(number)
        integer, intent(in) :: number
        call canvas%clear()
        call draw(min(number, movie_frames))     ! hold the last frame for 1 s
        call canvas%set_rendered(.false.)
    end subroutine update

    subroutine draw(f)
        integer, intent(in) :: f
        character(40), parameter :: names(3) = [character(40) :: 'exact (Q1 partial waves)', &
            'one Gaussian (M3)', 'coherent sum K = 8 (M4)']
        integer :: p
        write (arg, '(a,f5.2)') 'relative density, $\eta = 1$, b = 0:  $t v/z_0 = $', &
            mtime(f)*mv_speed/mv_z0
        call suptitle(trim(arg))
        do p = 1, 3
            call subplot(1, 3, p)
            call title(trim(names(p)))
            call xlabel('$x/\bar\lambda$'); call ylabel('$z/\bar\lambda$')
            call pcolormesh(xe/mv_lambdabar, xe/mv_lambdabar, img(:, :, p, f), cmap='viridis', &
                            vmin=vlo, vmax=vhi)
            call xlim(xe(0)/mv_lambdabar, xe(npx)/mv_lambdabar)
            call ylim(xe(0)/mv_lambdabar, xe(npx)/mv_lambdabar)
            call colorbar(label='log10 y-integrated density')
        end do
        call tight_layout()
    end subroutine draw
end program packet_models_demo
