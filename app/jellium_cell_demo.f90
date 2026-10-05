!> Q4 interacting jellium cell: finite-N electrons in a periodic cube with a
!> neutralizing background and bare Coulomb interaction, solved in exact
!> (S_z, K) sectors with the matrix-free Fock engine. Writes CSV evidence,
!> output/jellium_cell_demo.nc and four figures (pair correlation, static
!> response, cumulant size, kicked density-mode history).
program jellium_cell_demo
    use, intrinsic :: iso_fortran_env, only: dp => real64, i64 => int64
    use kin6d_quantum_fock, only: fock_hamiltonian_t, fock_model_info, &
        fock_rdms
    use kin6d_jellium_cell
    use kin6d_result_io, only: run_t, result_file_t
    use fortplot, only: figure, subplot, plot, fill_between, title, xlabel, &
        ylabel, xlim, ylim, legend, tight_layout, savefig, xscale, yscale, suptitle
    implicit none
    integer, parameter :: max_rs = 8, max_theta = 4, max_bins = 32
    real(dp), parameter :: blue(3) = [0.000_dp, 0.447_dp, 0.698_dp], &
        vermillion(3) = [0.835_dp, 0.369_dp, 0.000_dp], &
        green(3) = [0.000_dp, 0.620_dp, 0.451_dp], &
        orange(3) = [0.902_dp, 0.624_dp, 0.000_dp], &
        purple(3) = [0.800_dp, 0.475_dp, 0.655_dp], black(3) = 0.0_dp
    ! Configuration (&jellium_cell).
    integer :: cutoff = 2, nev = 8, nrs = 5, ntheta = 2, bins = 12
    integer :: scan_cutoff_max(4) = [0, 3, 3, 2], kick_steps = 800
    real(dp) :: rs_values(max_rs) = [0.5_dp, 1.0_dp, 2.0_dp, 4.0_dp, 8.0_dp, &
                                     0.0_dp, 0.0_dp, 0.0_dp]
    real(dp) :: thetas(max_theta) = [0.0_dp, 0.5_dp, 0.0_dp, 0.0_dp]
    real(dp) :: kick_dt = 0.01_dp, kick_rs = 2.0_dp, response_threshold = 1.0e-6_dp
    namelist /jellium_cell/ cutoff, nev, nrs, ntheta, bins, scan_cutoff_max, &
        kick_steps, rs_values, thetas, kick_dt, kick_rs, response_threshold
    integer, parameter :: q_list(3, 2) = reshape([1, 0, 0, 1, 1, 1], [3, 2])
    type(run_t) :: run
    type(result_file_t) :: result
    integer :: unit, ios
    character(256) :: msg
    ! Results gathered for figures and NetCDF.
    real(dp), allocatable :: g_same(:, :, :, :), g_opp(:, :, :, :), centers(:)
    real(dp), allocatable :: energy(:, :, :), cumulant(:, :, :, :), chi(:, :, :, :)
    real(dp), allocatable :: chi0(:, :, :, :), rpa(:, :, :, :), trunc(:, :, :)
    real(dp), allocatable :: times(:), response(:), band(:)
    real(dp) :: derivative_data(8)
    integer :: n, r, t

    call run%start('jellium_cell_demo')
    call run%open_config(unit)
    if (run%has_config()) read (unit, nml=jellium_cell, iostat=ios, iomsg=msg)
    call run%close_config(unit, 'jellium_cell', ios, msg)
    if (nrs < 1 .or. nrs > max_rs .or. ntheta < 1 .or. ntheta > max_theta) &
        error stop 'rs/theta list length'
    if (bins < 2 .or. bins > max_bins .or. cutoff < 1 .or. cutoff > 3) &
        error stop 'bins or cutoff out of range'
    if (any(rs_values(1:nrs) <= 0) .or. any(thetas(1:ntheta) < 0)) &
        error stop 'rs > 0 and theta >= 0 required'
    if (kick_steps < 2 .or. .not. (kick_dt > 0)) error stop 'kick time grid'
    call run%open_resolved(unit)
    write (unit, nml=jellium_cell)
    call run%close_resolved(unit)
    call execute_command_line('mkdir -p output', exitstat=ios)
    if (ios /= 0) error stop 'cannot create output'

    allocate (g_same(bins, ntheta, nrs, 2:4), g_opp(bins, ntheta, nrs, 2:4))
    allocate (energy(ntheta, nrs, 2:4), cumulant(4, ntheta, nrs, 2:4))
    allocate (chi(2, ntheta, nrs, 2:4), chi0(2, ntheta, nrs, 2:4), &
              rpa(2, ntheta, nrs, 2:4), trunc(ntheta, nrs, 2:4), centers(bins))
    allocate (times(0:kick_steps), response(0:kick_steps), band(0:kick_steps))
    g_same = 0; g_opp = 0; energy = 0; cumulant = 0; chi = 0; chi0 = 0; rpa = 0
    trunc = 0
    open (newunit=unit, file='output/jellium-cost-ledger.csv', status='replace')
    write (unit, '(a)') 'stage,N,cutoff,M,total_dim,max_sector_dim,rs,seconds,'// &
        'model_bytes,vector_bytes,matvecs_or_steps'
    call cutoff_scan(unit)
    do n = 2, 4
        do r = 1, nrs
            call rs_case(n, rs_values(r), r, unit)
        end do
    end do
    call kick_case(unit)
    call derivative_case(unit)
    close (unit)
    call write_tables()
    call write_netcdf()
    call make_figures()
    print '(a)', 'jellium_cell_demo: wrote output/jellium-*.csv, '// &
        'output/jellium_cell_demo.nc and output/jellium-*.png'
    if (.false.) print *, t

contains

    subroutine setup(n, c, rs, cell, model, spectrum, seconds, ok, g)
        integer, intent(in) :: n, c
        real(dp), intent(in) :: rs
        type(jellium_cell_t), intent(out) :: cell
        type(fock_hamiltonian_t), intent(out) :: model
        type(jellium_spectrum_t), intent(out) :: spectrum
        real(dp), intent(out) :: seconds
        logical, intent(out) :: ok
        real(dp), intent(in), optional :: g
        type(jellium_sector_t) :: sector
        integer, allocatable :: modes(:, :)
        real(dp) :: t0, t1
        call jellium_modes('spherical', c, modes)
        call init_jellium_cell(cell, n, modes, rs, ok, g)
        if (.not. ok) error stop 'cell'
        call jellium_sector(cell, n/2, reshape([0, 0, 0], [3, 1]), sector, ok)
        call jellium_hamiltonian(cell, sector%basis, model, ok)
        if (.not. ok) error stop 'model'
        call cpu_time(t0)
        call jellium_build_spectrum(cell, model, nev, spectrum, ok)
        call cpu_time(t1)
        if (.not. ok) error stop 'spectrum'
        seconds = t1 - t0
    end subroutine

    subroutine ledger(lu, stage, cell, spectrum, model, rs, seconds, count)
        integer, intent(in) :: lu, count
        character(*), intent(in) :: stage
        type(jellium_cell_t), intent(in) :: cell
        type(jellium_spectrum_t), intent(in) :: spectrum
        type(fock_hamiltonian_t), intent(in) :: model
        real(dp), intent(in) :: rs, seconds
        integer(i64) :: bytes
        logical :: bounded, ok
        call fock_model_info(model, bytes, bounded, ok)
        write (lu, '(a,",",i0,",",i0,",",i0,",",i0,",",i0,",",f6.2,",",es11.4,'// &
            '",",i0,",",i0,",",i0)') stage, cell%n, cutoff_of(cell), cell%m, &
            sum(spectrum%sectors(1:spectrum%nsectors)%d), &
            maxval(spectrum%sectors(1:spectrum%nsectors)%d), rs, seconds, bytes, &
            8_i64*size(spectrum%flat, kind=i64), count
    end subroutine

    integer function cutoff_of(cell) result(c)
        type(jellium_cell_t), intent(in) :: cell
        c = maxval(sum(cell%modes**2, 1))
    end function

    subroutine cutoff_scan(lu)
        integer, intent(in) :: lu
        type(jellium_cell_t) :: cell
        type(fock_hamiltonian_t) :: model
        type(jellium_spectrum_t) :: spectrum
        real(dp) :: seconds, ground, ideal, w(1)
        integer :: n, c, out
        logical :: ok
        open (newunit=out, file='output/jellium-cutoff-scan.csv', status='replace')
        write (out, '(a)') 'N,cutoff_n2,M,total_dim,rs,E0_hartree,E0_per_N,'// &
            'ideal_E0_per_N,correlation_indicator_per_N'
        do n = 2, 4
            do c = 1, scan_cutoff_max(n)
                call setup(n, c, 2.0_dp, cell, model, spectrum, seconds, ok)
                ground = minval(spectrum%energy(1:spectrum%count))
                ideal = ideal_ground(cell)
                write (out, '(i0,",",i0,",",i0,",",i0,",",f5.2,4(",",es16.8))') n, c, &
                    cell%m, sum(spectrum%sectors(1:spectrum%nsectors)%d), 2.0_dp, &
                    ground, ground/n, ideal/n, (ground - ideal)/n
                call ledger(lu, 'cutoff_spectrum', cell, spectrum, model, 2.0_dp, &
                            seconds, spectrum%count)
            end do
        end do
        close (out)
        if (.false.) w = 0
    end subroutine

    real(dp) function ideal_ground(cell) result(e)
        use kin6d_quantum_fock_leaves, only: qf_kinetic
        type(jellium_cell_t), intent(in) :: cell
        real(dp) :: eps(62), tmp
        integer :: i, j
        do i = 1, cell%m
            call qf_kinetic(cell%length, cell%hbar, cell%mass, &
                            real(sum(cell%momentum(:, i)**2), dp), eps(i))
        end do
        do i = 2, cell%m
            tmp = eps(i); j = i - 1
            do while (j >= 1)
                if (eps(j) <= tmp) exit
                eps(j + 1) = eps(j); j = j - 1
            end do
            eps(j + 1) = tmp
        end do
        e = sum(eps(1:cell%n))
    end function

    subroutine rs_case(n, rs, r, lu)
        integer, intent(in) :: n, r, lu
        real(dp), intent(in) :: rs
        type(jellium_cell_t) :: cell
        type(fock_hamiltonian_t) :: model
        type(jellium_spectrum_t) :: spectrum
        real(dp), allocatable :: w(:), chis(:), gamma(:, :), gamma2(:, :)
        real(dp) :: seconds, ground, edges(max_bins + 1), t0, t1, gc, rgc, norms(4)
        integer :: t, b, qi
        logical :: ok

        call setup(n, cutoff, rs, cell, model, spectrum, seconds, ok)
        call ledger(lu, 'spectrum', cell, spectrum, model, rs, seconds, spectrum%count)
        allocate (w(spectrum%count), chis(spectrum%count), gamma(cell%m, cell%m), &
                  gamma2(cell%m*(cell%m - 1)/2, cell%m*(cell%m - 1)/2))
        edges(1:bins + 1) = [(0.5_dp*cell%length*b/bins, b=0, bins)]
        centers = 0.5_dp*(edges(1:bins) + edges(2:bins + 1))/cell%length
        do t = 1, ntheta
            call jellium_thermal_weights(cell, spectrum, thetas(t), w, trunc(t, r, n), &
                                         ground)
            energy(t, r, n) = sum(w*spectrum%energy(1:spectrum%count))
            call cpu_time(t0)
            call jellium_thermal_rdms(cell, spectrum, w, gamma, gamma2, ok)
            call jellium_pair_bins(cell, gamma2, edges(1:bins + 1), &
                                   g_same(:, t, r, n), g_opp(:, t, r, n))
            call cpu_time(t1)
            call ledger(lu, 'rdm_pair_bins', cell, spectrum, model, rs, t1 - t0, 0)
            if (n >= 3) then
                call cpu_time(t0)
                call jellium_cumulants(cell, spectrum, w, norms, ok, 1.0e-8_dp)
                call cpu_time(t1)
                if (.not. ok) error stop 'cumulants'
                cumulant(:, t, r, n) = norms
                call ledger(lu, 'cumulants', cell, spectrum, model, rs, t1 - t0, 0)
            end if
            do qi = 1, 2
                call cpu_time(t0)
                call jellium_static_response(cell, model, spectrum, w, q_list(:, qi), &
                    chi(qi, t, r, n), chis, ok, response_threshold)
                call cpu_time(t1)
                if (.not. ok) error stop 'static response'
                call ledger(lu, 'static_response', cell, spectrum, model, rs, t1 - t0, &
                            count(chis /= 0))
                call jellium_ideal_response(cell, thetas(t), q_list(:, qi), &
                    chi0(qi, t, r, n), gc, rpa(qi, t, r, n), rgc, ok)
            end do
        end do
    end subroutine

    subroutine ground_state(cell, model, sector, psi, e0, gap)
        type(jellium_cell_t), intent(in) :: cell
        type(fock_hamiltonian_t), intent(in) :: model
        type(jellium_sector_t), intent(out) :: sector
        real(dp), allocatable, intent(out) :: psi(:)
        real(dp), intent(out) :: e0, gap
        real(dp) :: val(2), res(2)
        real(dp), allocatable :: vec(:, :)
        integer :: iterations
        logical :: ok
        call jellium_sector(cell, 1, reshape([0, 0, 0], [3, 1]), sector, ok)
        allocate (vec(sector%d, 2))
        call jellium_lowest(model, sector%basis, 2, val, vec, res, ok, iterations)
        if (.not. ok) error stop 'N=2 ground state'
        psi = vec(:, 1); e0 = val(1); gap = val(2) - val(1)
    end subroutine

    subroutine kick_case(lu)
        integer, intent(in) :: lu
        type(jellium_cell_t) :: cell
        type(fock_hamiltonian_t) :: model
        type(jellium_spectrum_t) :: spectrum
        type(jellium_sector_t) :: sector
        real(dp), allocatable :: psi(:)
        real(dp) :: seconds, e0, gap, tc, tz
        integer :: du, k, out
        logical :: ok
        call setup(2, cutoff, kick_rs, cell, model, spectrum, seconds, ok)
        call ground_state(cell, model, sector, psi, e0, gap)
        call jellium_kick_evolution(cell, model, sector, psi, e0, q_list(:, 1), &
            kick_steps, kick_dt, times, response, band, ok, du, tc, tz)
        if (.not. ok) error stop 'kick evolution certificate'
        call ledger(lu, 'kick_candidate_rk54', cell, spectrum, model, kick_rs, tc, &
                    kick_steps)
        call ledger(lu, 'kick_certificate', cell, spectrum, model, kick_rs, tz, du)
        open (newunit=out, file='output/jellium-mode-history.csv', status='replace')
        write (out, '(a)') 'time_hartree_inverse,response_per_kappa,certified_band'
        do k = 0, kick_steps
            write (out, '(es16.8,",",es18.10,",",es12.4)') times(k), response(k), band(k)
        end do
        close (out)
    end subroutine

    subroutine derivative_case(lu)
        integer, intent(in) :: lu
        type(jellium_cell_t) :: cell
        type(fock_hamiltonian_t) :: model, interaction, mp, mm
        type(jellium_spectrum_t) :: spectrum
        type(jellium_sector_t) :: sector
        real(dp), allocatable :: psi(:), pp(:), pm(:), zero(:), gr(:, :), gi(:, :), &
                                 g2r(:, :), g2i(:, :)
        real(dp) :: seconds, e0, gap, ep, em, de, edges(4), gs(3), go(3), dgs(3), dgo(3)
        real(dp) :: gsp(3), gop(3), gsm(3), gom(3), h, gbar, t0, t1, dummy
        integer :: iterations, p, b, out
        logical :: ok
        call setup(2, cutoff, kick_rs, cell, model, spectrum, seconds, ok)
        call ground_state(cell, model, sector, psi, e0, gap)
        call jellium_interaction_model(cell, sector%basis, interaction, ok)
        edges = [(0.5_dp*cell%length*b/3, b=0, 3)]
        call cpu_time(t0)
        call jellium_coupling_derivative(cell, sector%basis, model, interaction, psi, &
            e0, edges, de, gs, go, dgs, dgo, ok, iterations)
        call cpu_time(t1)
        if (.not. ok) error stop 'coupling derivative'
        call ledger(lu, 'derivative_jvp', cell, spectrum, model, kick_rs, t1 - t0, &
                    iterations)
        call cpu_time(t0)
        call jellium_coupling_vjp(cell, sector%basis, model, interaction, psi, e0, &
            edges, 1.0_dp, [0.0_dp, 0.0_dp, 0.0_dp], [1.0_dp, 0.0_dp, 0.0_dp], gbar, ok)
        call cpu_time(t1)
        call ledger(lu, 'derivative_vjp', cell, spectrum, model, kick_rs, t1 - t0, 0)
        h = 1.0e-4_dp
        call jellium_hamiltonian(cell, sector%basis, mp, ok, cell%coupling + h)
        call ground_state(cell, mp, sector, pp, ep, dummy)
        call jellium_hamiltonian(cell, sector%basis, mm, ok, cell%coupling - h)
        call ground_state(cell, mm, sector, pm, em, dummy)
        p = cell%m*(cell%m - 1)/2
        allocate (zero(sector%d), gr(cell%m, cell%m), gi(cell%m, cell%m), g2r(p, p), &
                  g2i(p, p))
        zero = 0
        call fock_rdms(sector%basis, pp, zero, gr, gi, g2r, g2i, ok)
        call jellium_pair_bins(cell, g2r, edges, gsp, gop)
        call fock_rdms(sector%basis, pm, zero, gr, gi, g2r, g2i, ok)
        call jellium_pair_bins(cell, g2r, edges, gsm, gom)
        derivative_data = [e0, de, (ep - em)/(2*h), go(1), dgo(1), &
                           (gop(1) - gom(1))/(2*h), gbar, de + dgo(1)]
        open (newunit=out, file='output/jellium-derivative.csv', status='replace')
        write (out, '(a)') 'quantity,value'
        write (out, '(a,es20.12)') 'E0,', e0
        write (out, '(a,es20.12)') 'dE_dg_hellmann_feynman,', de
        write (out, '(a,es20.12)') 'dE_dg_central_difference,', (ep - em)/(2*h)
        write (out, '(a,es20.12)') 'g_opp_bin1,', go(1)
        write (out, '(a,es20.12)') 'dg_opp_bin1_dg_eigen_derivative,', dgo(1)
        write (out, '(a,es20.12)') 'dg_opp_bin1_dg_central_difference,', &
            (gop(1) - gom(1))/(2*h)
        write (out, '(a,es20.12)') 'vjp_gbar_for_bar_E1_bar_bin1,', gbar
        write (out, '(a,es20.12)') 'jvp_pairing_dE_plus_dg_bin1,', de + dgo(1)
        write (out, '(a,es20.12)') 'spectral_gap,', gap
        close (out)
        if (.false.) print *, gs, dgs, gsp, gsm
    end subroutine

    subroutine write_tables()
        integer :: out, n, r, t, b, qi
        open (newunit=out, file='output/jellium-pair-correlation.csv', status='replace')
        write (out, '(a)') 'N,rs,theta,r_over_L,g_same,g_opp'
        do n = 2, 4
            do r = 1, nrs
                do t = 1, ntheta
                    do b = 1, bins
                        write (out, '(i0,",",f6.2,",",f5.2,",",f8.5,2(",",es15.7))') &
                            n, rs_values(r), thetas(t), centers(b), g_same(b, t, r, n), &
                            g_opp(b, t, r, n)
                    end do
                end do
            end do
        end do
        close (out)
        open (newunit=out, file='output/jellium-response.csv', status='replace')
        write (out, '(a)') 'N,rs,theta,q_n,chi,chi0_canonical_lindhard,'// &
            'chi_rpa_canonical,chi_over_chi0,energy_hartree,truncation_weight'
        do n = 2, 4
            do r = 1, nrs
                do t = 1, ntheta
                    do qi = 1, 2
                        write (out, '(i0,",",f6.2,",",f5.2,",",3i1,6(",",es15.7))') &
                            n, rs_values(r), thetas(t), q_list(:, qi), chi(qi, t, r, n), &
                            chi0(qi, t, r, n), rpa(qi, t, r, n), &
                            chi(qi, t, r, n)/chi0(qi, t, r, n), energy(t, r, n), &
                            trunc(t, r, n)
                    end do
                end do
            end do
        end do
        close (out)
        open (newunit=out, file='output/jellium-cumulants.csv', status='replace')
        write (out, '(a)') 'N,rs,theta,norm_Gamma2,norm_lambda2,norm_Gamma3,'// &
            'norm_lambda3,lambda2_over_Gamma2,lambda3_over_Gamma3'
        do n = 3, 4
            do r = 1, nrs
                do t = 1, ntheta
                    write (out, '(i0,",",f6.2,",",f5.2,6(",",es15.7))') n, &
                        rs_values(r), thetas(t), cumulant(:, t, r, n), &
                        cumulant(2, t, r, n)/cumulant(1, t, r, n), &
                        cumulant(4, t, r, n)/cumulant(3, t, r, n)
                end do
            end do
        end do
        close (out)
    end subroutine

    subroutine write_netcdf()
        character(8) :: tag
        integer :: n
        call result%create(run)
        call result%attribute('model', 'periodic jellium cell, bare Coulomb, '// &
            'neutralizing background, spherical plane-wave cutoff')
        call result%attribute('status', 'finite-basis candidate eigenpairs with '// &
            'recomputed residuals; kick band certified-discrete; cutoff, '// &
            'finite-size and thermodynamic-limit errors are indicators')
        call result%write('rs', rs_values(1:nrs), 'a0', 'Wigner-Seitz radius', dims='rs')
        call result%write('theta', thetas(1:ntheta), '1', 'T/T_F (0: ground multiplet)', &
                          dims='theta')
        call result%write('r_over_L', centers, '1', 'radial shell center / L', dims='bin')
        do n = 2, 4
            write (tag, '(a,i0)') '_N', n
            call result%write('g_same'//trim(tag), g_same(:, :, :, n), '1', &
                'same-spin pair correlation', dims='bin,theta,rs', status='candidate')
            call result%write('g_opp'//trim(tag), g_opp(:, :, :, n), '1', &
                'opposite-spin pair correlation', dims='bin,theta,rs', status='candidate')
            call result%write('chi'//trim(tag), chi(:, :, :, n), 'a0^-3 Ha^-1', &
                'static density response', dims='q,theta,rs', status='candidate')
            call result%write('chi0'//trim(tag), chi0(:, :, :, n), 'a0^-3 Ha^-1', &
                'canonical finite-cell Lindhard', dims='q,theta,rs')
            call result%write('chi_rpa'//trim(tag), rpa(:, :, :, n), 'a0^-3 Ha^-1', &
                'RPA from canonical Lindhard', dims='q,theta,rs')
            call result%write('cumulant_norms'//trim(tag), cumulant(:, :, :, n), '1', &
                '|Gamma2|,|lambda2|,|Gamma3|,|lambda3|', dims='norm,theta,rs')
        end do
        call result%write('kick_time', times, 'hbar/Ha', 'time after kick', dims='frame')
        call result%write('kick_response', response, 'per kappa', &
            'delta<A>(t)/kappa', dims='frame', status='candidate center')
        call result%write('kick_band', band, 'per kappa', 'finite-basis response bound', &
            dims='frame', status='certified-discrete')
        call result%write('derivatives', derivative_data, 'mixed', &
            'E0,dE/dg,FD,g1,dg1/dg,FD,vjp,jvp pairing', dims='derivative')
        call result%close()
    end subroutine

    subroutine make_figures()
        real(dp), parameter :: palette(3, 5) = reshape([0.000_dp, 0.447_dp, 0.698_dp, &
            0.835_dp, 0.369_dp, 0.000_dp, 0.000_dp, 0.620_dp, 0.451_dp, &
            0.902_dp, 0.624_dp, 0.000_dp, 0.800_dp, 0.475_dp, 0.655_dp], [3, 5])
        character(32) :: label
        integer :: r, t, qi, n, k
        real(dp), allocatable :: ratio(:)
        t = min(2, ntheta)
        ! Pair correlation, N=4 at the second temperature.
        call figure(figsize=[12.8_dp, 6.4_dp])
        do k = 1, 2
            call subplot(1, 2, k)
            if (k == 1) call title('Same-spin pair correlation, N = 4')
            if (k == 2) call title('Opposite-spin pair correlation, N = 4')
            call xlabel('r / L'); call ylabel('g(r)')
            call xlim(0.0_dp, 0.5_dp); call ylim(0.0_dp, 1.4_dp)
            do r = 1, nrs
                write (label, '(a,f4.1)') '$r_s$ = ', rs_values(r)
                if (k == 1) call plot(centers, g_same(:, t, r, 4), &
                    color=palette(:, 1 + mod(r - 1, 5)), marker='o', label=trim(label))
                if (k == 2) call plot(centers, g_opp(:, t, r, 4), &
                    color=palette(:, 1 + mod(r - 1, 5)), marker='o', label=trim(label))
            end do
            call plot([0.0_dp, 0.5_dp], [1.0_dp, 1.0_dp], color=black, linestyle='--', &
                      label='uncorrelated')
            call legend(loc='lower right')
        end do
        write (label, '(a,f4.2)') 'Θ = ', thetas(t)
        call suptitle('Exchange and Coulomb holes in the jellium cell, '//trim(label))
        call tight_layout()
        call savefig('output/jellium-pair-correlation.png')
        ! Static response relative to the same-cell canonical Lindhard function.
        call figure(figsize=[12.8_dp, 6.4_dp])
        allocate (ratio(nrs))
        do qi = 1, 2
            call subplot(1, 2, qi)
            write (label, '(a,3i1,a)') 'q L/2π = (', q_list(:, qi), ')'
            call title('Static response, '//trim(label))
            call xlabel('$r_s$'); call ylabel('$\chi / \chi_0$')
            call xscale('log')
            do n = 2, 4, 2
                k = merge(1, t, n == 2)
                ratio = chi(qi, k, 1:nrs, n)/chi0(qi, k, 1:nrs, n)
                write (label, '(a,i0)') 'interacting N = ', n
                call plot(rs_values(1:nrs), ratio, color=blue, marker='o', &
                    linestyle=merge('- ', ': ', n == 2), label=trim(label))
                ratio = rpa(qi, k, 1:nrs, n)/chi0(qi, k, 1:nrs, n)
                write (label, '(a,i0)') 'RPA N = ', n
                call plot(rs_values(1:nrs), ratio, color=vermillion, marker='s', &
                    linestyle=merge('- ', ': ', n == 2), label=trim(label))
            end do
            call plot([rs_values(1), rs_values(nrs)], [1.0_dp, 1.0_dp], color=black, &
                      linestyle='--', label='Lindhard (same cell)')
            call legend(loc='lower left')
        end do
        write (label, '(a,f4.2)') 'N = 4 at Θ = ', thetas(t)
        call suptitle('N = 2 ground state, '//trim(label))
        call tight_layout()
        call savefig('output/jellium-response.png')
        ! Cumulant size versus coupling.
        call figure(figsize=[12.8_dp, 6.4_dp])
        do k = 1, 2
            call subplot(1, 2, k)
            if (k == 1) call title('Pair cumulant $|\lambda_2| / |\Gamma_2|$')
            if (k == 2) call title('Three-body cumulant $|\lambda_3| / |\Gamma_3|$')
            call xlabel('$r_s$'); call ylabel('relative Frobenius norm')
            call xscale('log')
            do n = 3, 4
                do r = 1, ntheta
                    ratio = cumulant(2*k, r, 1:nrs, n)/cumulant(2*k - 1, r, 1:nrs, n)
                    write (label, '(a,i0,a,f4.2)') 'N = ', n, ', Θ = ', thetas(r)
                    call plot(rs_values(1:nrs), ratio, color=palette(:, 2*(n - 3) + r), &
                        marker='o', linestyle=merge('- ', ': ', r == 1), label=trim(label))
                end do
            end do
            call legend(loc='upper left')
        end do
        call suptitle('Connected correlations of the finite jellium cell')
        call tight_layout()
        call savefig('output/jellium-cumulants.png')
        ! Kicked density-mode history with its certified finite-basis band.
        call figure(figsize=[12.8_dp, 6.4_dp])
        call subplot(1, 2, 1)
        call title('Density mode after a weak kick')
        call xlabel('t [ħ/Ha]'); call ylabel('δ<A>(t) / κ')
        call xlim(0.0_dp, times(kick_steps))
        call fill_between(times, response - band, response + band, color=blue, &
                          alpha=0.25_dp)
        call plot(times, response, color=blue, label='numerical ± certified bound')
        call legend(loc='upper right')
        call subplot(1, 2, 2)
        call title('Certified finite-basis response bound')
        call xlabel('t [ħ/Ha]'); call ylabel('bound per κ')
        call xlim(0.0_dp, times(kick_steps))
        call yscale('log')
        call plot(times(1:kick_steps), band(1:kick_steps), color=blue, label='2|φ0| e(t)')
        call legend(loc='lower right')
        write (label, '(a,f4.1)') '$r_s$ = ', kick_rs
        call suptitle('N = 2, q L/2π = (100), '//trim(label)//'; recurrences, not damping')
        call tight_layout()
        call savefig('output/jellium-mode-history.png')
        if (.false.) print *, green, orange, purple
    end subroutine
end program jellium_cell_demo
