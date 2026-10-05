! Q1/Q1b evidence: Sommerfeld-eta and impact-parameter scans of the localized
! two-body quantum Coulomb encounter (product, singlet, triplet, classical),
! an encounter movie with relative and lab-frame densities beside the classical
! ensemble, and the hydrogen Kepler packet. Writes CSV and figures to output/.
! Usage: coulomb_encounter_demo [--no-movie]
program coulomb_encounter_demo
    use, intrinsic :: iso_fortran_env, only: dp => real64, i64 => int64
    use kin6d_coulomb_encounter
    use kin6d_coulomb_kepler
    use kin6d_coulomb_encounter_leaves, only: ce_mott, ce_free_width
    use fortplot, only: figure, get_global_figure, figure_t, subplot, plot, title, &
        xlabel, ylabel, xlim, ylim, legend, savefig, pcolormesh, tight_layout, &
        set_yscale, suptitle
    use fortplot_animation, only: animation_t, animation_func => FuncAnimation, &
        save_animation
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp), parameter :: etas(4) = [0.3_dp, 1.0_dp, 3.0_dp, 10.0_dp]
    integer, parameter :: nb = 36, nclass = 20000
    logical :: movie
    character(160) :: arg
    integer :: status, i
    real(dp) :: ratio(nb, 5, 4), theta(nb)
    ! Movie state shared with the animation callback (host association).
    integer, parameter :: nf = 192, npx = 128, nsep = 60
    real(dp), allocatable :: qr(:, :, :), ql(:, :, :), tl(:, :, :), cr(:, :, :), &
        cl(:, :, :), sep(:, :, :)
    real(dp) :: times(nf), xr(0:npx), xl(0:npx), scent(nsep), lr, mv_lambdabar, &
        mv_speed, mv_z0, mv_r0_sigma
    type(figure_t), pointer :: canvas => null()
    movie = .true.
    do i = 1, command_argument_count()
        call get_command_argument(i, arg)
        if (trim(arg) == '--no-movie') movie = .false.
    end do
    call execute_command_line('mkdir -p output', exitstat=status)
    if (status /= 0) error stop 'cannot create output/'
    call eta_scan()
    call b_scan()
    call kepler_case()
    call cost_ledger()
    if (movie) call encounter_movie()
contains
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

    ! Sommerfeld-eta scan: finite-time angular bins versus plane-wave Mott and
    ! Rutherford references (fluence from the free transverse width at z0/v).
    subroutine eta_scan()
        type(encounter_config_t) :: cfg
        type(encounter_t) :: enc
        complex(dp), allocatable :: a(:, :, :)
        real(dp) :: edges(0:nb), pp(nb), ps(nb), pt(nb), pc(nb), psym(nb), pcs(nb), t
        real(dp) :: sig_l, width, nperp, ref, ms, xq(16), wq(16), d, ex, cr, eb, rad
        real(dp) :: khi, uniform, dof, wall
        integer :: ie, ib, q, u, us, l
        open (newunit=u, file='output/coulomb-eta-bins.csv', status='replace')
        write (u, '(a)') 'eta,theta_lo_deg,theta_hi_deg,P_product,P_product_sym,'// &
            'P_singlet,P_triplet,P_classical_sym,R_sym_fluence,mott_singlet_ratio,'// &
            'mott_triplet_ratio'
        open (newunit=us, file='output/coulomb-eta-summary.csv', status='replace')
        write (us, '(a)') 'eta,k,lambdabar,r0,sigma_over_lambdabar,r0_over_sigma,z0,'// &
            't_final,rmax,n_radial,lmax,mmax,kept_states,partial_wave_dof,'// &
            'uniform_grid_points,setup_s,eigen_s,projection_s,stored_MB,e0,rho,'// &
            'state_bound_t,projection_defect,reconstruction_rms,wall_prob,l_edge2,'// &
            'spectral_discarded2'
        do ib = 0, nb
            edges(ib) = pi*ib/nb
        end do
        theta = (edges(1:) + edges(:nb - 1))/2*180/pi
        do ie = 1, size(etas)
            sig_l = max(10.0_dp, 4*etas(ie))
            call setup_case(etas(ie), sig_l, 6*sig_l, 0.0_dp, cfg, enc, t)
            if (allocated(a)) deallocate (a)
            allocate (a(enc%mesh%n, 0:enc%mmax, 0:enc%lmax))
            call encounter_state(enc, t, a)
            call encounter_angular_bins(enc, a, SECTOR_PRODUCT, cfg%z0/4, edges, pp)
            call encounter_angular_bins(enc, a, SECTOR_SINGLET, cfg%z0/4, edges, ps)
            call encounter_angular_bins(enc, a, SECTOR_TRIPLET, cfg%z0/4, edges, pt)
            call classical_bins(enc, t, cfg%z0/4, edges, pc)
            call ce_free_width(cfg%sigma, enc%mu, cfg%hbar, cfg%z0/enc%speed, width)
            nperp = 1/(2*pi*width**2)
            do ib = 1, nb
                psym(ib) = (pp(ib) + pp(nb + 1 - ib))/2
                pcs(ib) = (pc(ib) + pc(nb + 1 - ib))/2
            end do
            do ib = 1, nb
                call gauss_nodes(edges(ib - 1), edges(ib), xq, wq)
                ref = 0; ms = 0
                do q = 1, 16
                    call ce_mott(acos(xq(q)), enc%eta, cfg%k, d, ex, cr)
                    ref = ref + 2*pi*wq(q)*(d + ex)/2*nperp
                    ms = ms + 2*pi*wq(q)*(d + ex + cr)/2*nperp
                end do
                write (u, '(f6.2,2(",",f8.3),9(",",es14.6))') etas(ie), edges(ib - 1)*180/pi, &
                    edges(ib)*180/pi, pp(ib), psym(ib), ps(ib), pt(ib), pcs(ib), ref, &
                    ms/ref, 2 - ms/ref
                ratio(ib, :, ie) = [psym(ib)/ref, ps(ib)/ref, pt(ib)/ref, pcs(ib)/ref, ms/ref]
            end do
            khi = cfg%k + 8/(2*cfg%sigma)
            uniform = (2*enc%mesh%rmax*khi*cfg%nodes_per_wavelength/(2*pi))**3
            dof = 0
            do l = 0, enc%lmax
                dof = dof + real(enc%mesh%n*(min(l, enc%mmax) + 1), dp)
            end do
            call encounter_state_error(enc, t, 1.0_dp, 0.5_dp, eb, rad)
            wall = encounter_wall_probability(enc, a, SECTOR_PRODUCT, 2*cfg%sigma)
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
        close (u); close (us)
        call figure(figsize=[12.0_dp, 8.5_dp])
        do ie = 1, size(etas)
            call subplot(2, 2, ie)
            write (arg, '(a,f4.1,a,f5.2)') '$\eta$ = ', etas(ie), ',  $r_0/\sigma$ = ', &
                2*etas(ie)/max(10.0_dp, 4*etas(ie))
            call title(trim(arg))
            call xlabel('$\theta$ (deg), CM frame, bins at $t = 2 z_0/v$')
            call ylabel('bin probability / plane-wave Rutherford (sym.)')
            call xlim(30.0_dp, 150.0_dp); call ylim(0.0_dp, 2.6_dp)
            call plot(theta, ratio(:, 5, ie), label='Mott singlet (plane wave)', &
                      linestyle='-', color=[0.6_dp, 0.6_dp, 0.6_dp])
            call plot(theta, 2 - ratio(:, 5, ie), label='Mott triplet (plane wave)', &
                      linestyle='--', color=[0.6_dp, 0.6_dp, 0.6_dp])
            call plot(theta, ratio(:, 1, ie), label='quantum product, spin-unresolved', &
                      linestyle='-', color=[0.0_dp, 0.0_dp, 0.0_dp])
            call plot(theta, ratio(:, 2, ie), label='quantum singlet', linestyle='-', &
                      color=[0.0_dp, 0.447_dp, 0.698_dp])
            call plot(theta, ratio(:, 3, ie), label='quantum triplet', linestyle='-', &
                      color=[0.835_dp, 0.369_dp, 0.0_dp])
            call plot(theta, ratio(:, 4, ie), label='classical ensemble', linestyle=':', &
                      color=[0.0_dp, 0.62_dp, 0.45_dp])
            if (ie == 1) call legend(loc='upper left')
        end do
        call tight_layout()
        call savefig('output/coulomb-mott-eta-scan.png')
    end subroutine eta_scan

    subroutine gauss_nodes(th1, th2, xq, wq)
        use fortnum_quadrature, only: gauss_legendre_ab
        real(dp), intent(in) :: th1, th2
        real(dp), intent(out) :: xq(:), wq(:)
        call gauss_legendre_ab(size(xq), cos(th2), cos(th1), xq, wq)
    end subroutine gauss_nodes

    ! Impact-parameter scan at eta = 1 (distinguishable product state).
    subroutine b_scan()
        real(dp), parameter :: bs(4) = [0.0_dp, 5.0_dp, 10.0_dp, 20.0_dp]
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
            call classical_bins(enc, t, cfg%z0/4, edges, pc(:, ibs))
            do ib = 1, nb
                write (u, '(f6.2,",",f6.3,2(",",f8.3),2(",",es14.6),2(",",i0))') bs(ibs), &
                    bs(ibs)/10, edges(ib - 1)*180/pi, edges(ib)*180/pi, pq(ib, ibs), &
                    pc(ib, ibs), enc%lmax, enc%mmax
            end do
            print '(a,f5.1,a,i0,a,i0)', 'b/lambdabar=', bs(ibs), ' lmax=', enc%lmax, &
                ' mmax=', enc%mmax
        end do
        close (u)
        call figure(figsize=[9.0_dp, 6.0_dp])
        call title('eta = 1, sigma = 10 lambdabar: impact-parameter scan (product state)')
        call xlabel('$\theta$ (deg), CM frame, bins at $t = 2 z_0/v$')
        call ylabel('bin probability (5 deg bins, |r| > z0/4)')
        call set_yscale('log')
        call xlim(0.0_dp, 180.0_dp); call ylim(1.0e-6_dp, 1.0_dp)
        do ibs = 1, size(bs)
            write (arg, '(a,f4.1,a)') 'quantum b = ', bs(ibs), ' lambdabar'
            call plot(theta, max(pq(:, ibs), 1.0e-7_dp), label=trim(arg), linestyle='-', &
                      color=palette(ibs))
            write (arg, '(a,f4.1)') 'classical b = ', bs(ibs)
            call plot(theta, max(pc(:, ibs), 1.0e-7_dp), label=trim(arg), linestyle=':', &
                      color=palette(ibs))
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
        call figure(figsize=[11.0_dp, 10.0_dp])
        call subplot(3, 1, 1)
        call title('Q1b: hydrogen radial packet, nbar = 40, l = 1, rms dn = 1.5')
        call xlabel('t / T_rev   (T_rev = (2 nbar/3) T_K)')
        call ylabel('<r> / r_out')
        call xlim(0.0_dp, 1.1_dp)
        call plot(t/kep%revival, rmean/kep%rc, linestyle='-', color=[0.0_dp, 0.447_dp, 0.698_dp])
        call subplot(3, 1, 2)
        call title('autocorrelation |A(t)|^2: represented spectrum vs exact hydrogen levels')
        call xlabel('t / T_rev'); call ylabel('|<psi(0)|psi(t)>|^2')
        call xlim(0.0_dp, 1.1_dp); call ylim(0.0_dp, 1.05_dp)
        call plot(t/kep%revival, a2, label='FEM-DVR spectrum', linestyle='-', &
                  color=[0.0_dp, 0.0_dp, 0.0_dp])
        call plot(t/kep%revival, o2, label='exact E_n oracle', linestyle=':', &
                  color=[0.835_dp, 0.369_dp, 0.0_dp])
        call legend(loc='upper right')
        call subplot(3, 1, 3)
        call title('radial probability density P(r, t)')
        call xlabel('t / T_rev'); call ylabel('r / r_out')
        tedge = t/kep%revival
        do i = 1, ns - 1
            img(:, i) = dens(:, i)*kep%rc
        end do
        call pcolormesh(tedge, edges/kep%rc, img, cmap='viridis', vmin=0.0_dp, &
                        vmax=maxval(img)*0.6_dp)
        call tight_layout()
        call savefig('output/coulomb-kepler.png')
    end subroutine kepler_case

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
        print '(a,es9.2,a,es9.2,a,es9.2,a,es9.2)', 'cost: state ', ts, ' s, bins ', tb, &
            ' s, jvp ', tj, ' s, vjp ', tv
    end subroutine cost_ledger

    ! Encounter movie (eta = 3, sigma = 3 lambdabar, b = 1 lambdabar): quantum
    ! relative and lab-frame one-body densities beside the classical ensemble on
    ! the same axes and clock, plus relative-separation densities.
    subroutine encounter_movie()
        type(encounter_config_t) :: cfg
        type(encounter_t) :: enc
        type(animation_t) :: anim
        complex(dp), allocatable :: a(:, :, :)
        real(dp), allocatable :: rel(:, :, :), cm(:, :, :), tri(:, :), tmp(:, :)
        real(dp) :: t, tf, sedges(0:nsep), ll, wcm, drift, x1(3), x2(3)
        integer :: f, n, ix, u, ib
        logical :: ok
        call setup_case(3.0_dp, 3.0_dp, 20.0_dp, 1.0_dp, cfg, enc, tf)
        mv_lambdabar = enc%lambdabar; mv_speed = enc%speed; mv_z0 = cfg%z0
        mv_r0_sigma = enc%r0/cfg%sigma
        call ce_free_width(cfg%sigma, enc%mu, cfg%hbar, tf, wcm)
        lr = min(enc%mesh%rmax, cfg%z0 + 3*wcm)
        call ce_free_width(enc%sigma_cm, enc%mtotal, cfg%hbar, tf, wcm)
        ll = lr/2 + 2*wcm
        do ix = 0, npx
            xr(ix) = -lr + 2*lr*ix/npx; xl(ix) = -ll + 2*ll*ix/npx
        end do
        do ib = 0, nsep
            sedges(ib) = lr*ib/nsep
        end do
        scent = (sedges(1:) + sedges(:nsep - 1))/2
        do f = 1, nf
            times(f) = tf*real(f - 1, dp)/real(nf - 1, dp)
        end do
        allocate (a(enc%mesh%n, 0:enc%mmax, 0:enc%lmax), rel(3, nclass, nf), cm(3, nclass, nf))
        allocate (qr(npx, npx, nf), ql(npx, npx, nf), tl(npx, npx, nf), cr(npx, npx, nf), &
                  cl(npx, npx, nf), sep(nsep, 4, nf), tri(npx, npx), tmp(npx, npx))
        call encounter_classical_ensemble(enc, nclass, 77_i64, times, rel, cm, drift, ok)
        if (.not. ok) error stop 'classical movie ensemble failed'
        open (newunit=u, file='output/coulomb-movie-separation.csv', status='replace')
        write (u, '(a)') 't,t_v_over_z0,r_over_lambdabar,P_product,P_singlet,P_triplet,'// &
            'P_classical,norm,energy,wall_prob'
        do f = 1, nf
            t = times(f)
            call encounter_state(enc, t, a)
            call encounter_relative_projection(enc, a, SECTOR_PRODUCT, .false., xr, xr, npx, &
                                               enc%lambdabar/4, qr(:, :, f))
            call encounter_relative_projection(enc, a, SECTOR_TRIPLET, .false., xr, xr, npx, &
                                               enc%lambdabar/4, tri)
            call encounter_lab_projection(enc, qr(:, :, f), xr, xr, t, 1, xl, xl, ql(:, :, f))
            call encounter_lab_projection(enc, qr(:, :, f), xr, xr, t, 2, xl, xl, tmp)
            ql(:, :, f) = ql(:, :, f) + tmp
            call encounter_lab_projection(enc, tri, xr, xr, t, 1, xl, xl, tl(:, :, f))
            call encounter_lab_projection(enc, tri, xr, xr, t, 2, xl, xl, tmp)
            tl(:, :, f) = tl(:, :, f) + tmp
            call encounter_separation_density(enc, a, SECTOR_PRODUCT, sedges, sep(:, 1, f))
            call encounter_separation_density(enc, a, SECTOR_SINGLET, sedges, sep(:, 2, f))
            call encounter_separation_density(enc, a, SECTOR_TRIPLET, sedges, sep(:, 3, f))
            cr(:, :, f) = 0; cl(:, :, f) = 0; sep(:, 4, f) = 0
            do n = 1, nclass
                call deposit(cr(:, :, f), xr, rel(1, n, f), rel(3, n, f))
                x1 = cm(:, n, f) + rel(:, n, f)/2; x2 = cm(:, n, f) - rel(:, n, f)/2
                call deposit(cl(:, :, f), xl, x1(1), x1(3))
                call deposit(cl(:, :, f), xl, x2(1), x2(3))
                ib = int(norm2(rel(:, n, f))/(sedges(1) - sedges(0))) + 1
                if (ib <= nsep) sep(ib, 4, f) = sep(ib, 4, f) + 1.0_dp/nclass
            end do
            cr(:, :, f) = cr(:, :, f)/(nclass*(xr(1) - xr(0))**2)
            cl(:, :, f) = cl(:, :, f)/(nclass*(xl(1) - xl(0))**2)
            do ib = 1, nsep
                write (u, '(es13.5,9(",",es13.5))') t, t*enc%speed/cfg%z0, &
                    scent(ib)/enc%lambdabar, sep(ib, :, f), &
                    encounter_norm(enc, a, SECTOR_PRODUCT), &
                    encounter_energy(enc, a, SECTOR_PRODUCT), &
                    encounter_wall_probability(enc, a, SECTOR_PRODUCT, 2*cfg%sigma)
            end do
        end do
        close (u)
        print '(a,i0,a,es9.2)', 'movie frames=', nf, ' classical drift=', drift
        call figure(figsize=[15.0_dp, 9.0_dp])
        canvas => get_global_figure()
        do f = 1, nf, nf/3
            call draw(f)
            write (arg, '(a,i3.3,a)') 'output/coulomb-encounter-frame', f, '.png'
            call savefig(trim(arg))
            call canvas%clear()
        end do
        call draw(nf)
        call savefig('output/coulomb-encounter-final.png')
        anim = animation_func(update, frames=nf, interval=42, fig=canvas)
        call save_animation(anim, 'output/coulomb-encounter.mp4', fps=24, status=status)
        if (status /= 0) print '(a,i0)', 'movie encoding status ', status
    end subroutine encounter_movie
    subroutine deposit(img, edges, x, z)
        real(dp), intent(inout) :: img(:, :)
        real(dp), intent(in) :: edges(0:), x, z
        integer :: i, j
        i = floor((x - edges(0))/(edges(1) - edges(0))) + 1
        j = floor((z - edges(0))/(edges(1) - edges(0))) + 1
        if (i >= 1 .and. i <= size(img, 2) .and. j >= 1 .and. j <= size(img, 1)) &
            img(j, i) = img(j, i) + 1
    end subroutine deposit

    subroutine update(number)
        integer, intent(in) :: number
        call canvas%clear()
        call draw(number)
        call canvas%set_rendered(.false.)
    end subroutine update

    subroutine draw(f)
        integer, intent(in) :: f
        real(dp) :: vr, vl, lb
        lb = mv_lambdabar
        ! Fixed colour range for the whole movie (quantum peak).
        vr = maxval(qr); vl = max(maxval(ql), maxval(tl))
        write (arg, '(a,f5.2,a,f5.2,a)') 'Q1 two-electron encounter, $\eta = 3$, '// &
            '$r_0/\sigma = $', mv_r0_sigma, ',  $t v/z_0 = $', times(f)*mv_speed/mv_z0, &
            '; fixed colour range'
        call suptitle(trim(arg))
        call subplot(2, 3, 1)
        call title('quantum relative density (x,z), y-integrated')
        call xlabel('$k x$'); call ylabel('$k z$')
        call pcolormesh(xr/lb, xr/lb, qr(:, :, f), cmap='viridis', vmin=0.0_dp, vmax=vr)
        call subplot(2, 3, 2)
        call title('quantum one-body density, both electrons (up,down)')
        call xlabel('$k x$'); call ylabel('$k z$')
        call pcolormesh(xl/lb, xl/lb, ql(:, :, f), cmap='viridis', vmin=0.0_dp, vmax=vl)
        call subplot(2, 3, 3)
        call title('quantum one-body density, same spin (triplet)')
        call xlabel('$k x$'); call ylabel('$k z$')
        call pcolormesh(xl/lb, xl/lb, tl(:, :, f), cmap='viridis', vmin=0.0_dp, vmax=vl)
        call subplot(2, 3, 4)
        write (arg, '(a,i0,a)') 'classical relative ensemble (', nclass, ' pairs)'
        call title(trim(arg))
        call xlabel('$k x$'); call ylabel('$k z$')
        call pcolormesh(xr/lb, xr/lb, cr(:, :, f), cmap='viridis', vmin=0.0_dp, vmax=vr)
        call subplot(2, 3, 5)
        call title('classical one-body density, both electrons')
        call xlabel('$k x$'); call ylabel('$k z$')
        call pcolormesh(xl/lb, xl/lb, cl(:, :, f), cmap='viridis', vmin=0.0_dp, vmax=vl)
        call subplot(2, 3, 6)
        call title('relative separation |r|: probability per bin')
        call xlabel('$k |r|$'); call ylabel('probability per bin')
        call xlim(0.0_dp, lr/lb); call ylim(0.0_dp, 1.2_dp*maxval(sep(:, :, :)))
        call plot(scent/lb, sep(:, 1, f), label='quantum product', linestyle='-', &
                  color=[0.0_dp, 0.0_dp, 0.0_dp])
        call plot(scent/lb, sep(:, 2, f), label='quantum singlet', linestyle='-', &
                  color=[0.0_dp, 0.447_dp, 0.698_dp])
        call plot(scent/lb, sep(:, 3, f), label='quantum triplet', linestyle='--', &
                  color=[0.835_dp, 0.369_dp, 0.0_dp])
        call plot(scent/lb, sep(:, 4, f), label='classical ensemble', linestyle=':', &
                  color=[0.0_dp, 0.62_dp, 0.45_dp])
        call legend(loc='upper right')
    end subroutine draw
end program coulomb_encounter_demo
