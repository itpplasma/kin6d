!> Q5 mean-field models of the Q4 jellium cell: Hartree, RHF/UHF, Mermin
!> finite-T HF and static TDHF response in the same plane-wave basis, N, spin
!> sector, r_s and cutoff as the exact sector reference, compared one
!> assumption at a time. Writes output/mean-field-*.csv,
!> output/mean_field_cell_demo.nc and four figures.
program mean_field_cell_demo
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_quantum_fock, only: fock_hamiltonian_t
    use kin6d_jellium_cell
    use kin6d_mean_field_cell
    use kin6d_mean_field_cell_leaves, only: mf_pz81_high, mf_pz81_low
    use kin6d_result_io, only: run_t, result_file_t
    use fortplot, only: figure, subplot, plot, title, xlabel, ylabel, xlim, ylim, &
        legend, tight_layout, savefig, xscale, yscale, suptitle
    implicit none
    integer, parameter :: max_rs = 8, max_bins = 32, nq = 1
    integer, parameter :: q(3) = [1, 0, 0]
    real(dp), parameter :: blue(3) = [0.000_dp, 0.447_dp, 0.698_dp], &
        vermillion(3) = [0.835_dp, 0.369_dp, 0.000_dp], &
        green(3) = [0.000_dp, 0.620_dp, 0.451_dp], &
        orange(3) = [0.902_dp, 0.624_dp, 0.000_dp], &
        purple(3) = [0.800_dp, 0.475_dp, 0.655_dp], &
        sky(3) = [0.337_dp, 0.706_dp, 0.914_dp], black(3) = 0.0_dp
    ! Configuration (&mean_field_cell).
    integer :: cutoff = 2, nev = 8, nrs = 5, nstarts = 6, bins = 12
    real(dp) :: rs_values(max_rs) = [0.5_dp, 1.0_dp, 2.0_dp, 4.0_dp, 8.0_dp, &
                                     0.0_dp, 0.0_dp, 0.0_dp]
    real(dp) :: theta_thermal = 0.5_dp, response_threshold = 1.0e-6_dp
    namelist /mean_field_cell/ cutoff, nev, nrs, nstarts, bins, rs_values, &
        theta_thermal, response_threshold
    type(run_t) :: run
    type(result_file_t) :: result
    integer :: unit, ios, n, r
    character(256) :: msg
    ! energy(:, r, n): exact, Hartree, RHF, UHF, UHF solutions, UHF converged
    ! runs, UHF K-sector count, TDHF Hessian minimum, UHF nup.
    real(dp), allocatable :: energy(:, :, :), estimate(:, :, :), response(:, :, :)
    real(dp), allocatable :: thermal(:, :, :), cumulant(:, :, :), derivative(:, :, :)
    real(dp), allocatable :: pair(:, :, :, :), centers(:)

    call run%start('mean_field_cell_demo')
    call run%open_config(unit)
    if (run%has_config()) read (unit, nml=mean_field_cell, iostat=ios, iomsg=msg)
    call run%close_config(unit, 'mean_field_cell', ios, msg)
    if (nrs < 1 .or. nrs > max_rs .or. any(rs_values(1:nrs) <= 0)) &
        error stop 'rs list'
    if (bins < 2 .or. bins > max_bins .or. cutoff < 1 .or. cutoff > 3) &
        error stop 'bins or cutoff out of range'
    if (nstarts < 1 .or. .not. (theta_thermal > 0)) error stop 'nstarts, theta'
    call run%open_resolved(unit)
    write (unit, nml=mean_field_cell)
    call run%close_resolved(unit)
    call execute_command_line('mkdir -p output', exitstat=ios)
    if (ios /= 0) error stop 'cannot create output'
    allocate (energy(9, nrs, 2:4), estimate(10, nrs, 2:4), response(5, nrs, 2:4), &
              thermal(9, nrs, 2:4), cumulant(4, nrs, 2:4), derivative(4, nrs, 2:4), &
              pair(bins, 4, nrs, 2:4), centers(bins))
    energy = 0; estimate = 0; response = 0; thermal = 0; cumulant = 0
    derivative = 0; pair = 0
    open (newunit=unit, file='output/mean-field-cost.csv', status='replace')
    write (unit, '(a)') 'stage,N,rs,seconds'
    do n = 2, 4
        do r = 1, nrs
            call rs_case(n, r, unit)
        end do
    end do
    close (unit)
    call write_tables()
    call write_netcdf()
    call make_figures()
    print '(a)', 'mean_field_cell_demo: wrote output/mean-field-*.csv, '// &
        'output/mean_field_cell_demo.nc and output/mean-field-*.png'

contains

    subroutine cost(lu, stage, n, rs, t0)
        integer, intent(in) :: lu, n
        character(*), intent(in) :: stage
        real(dp), intent(in) :: rs
        real(dp), intent(inout) :: t0
        real(dp) :: t1
        call cpu_time(t1)
        write (lu, '(a,",",i0,",",f6.2,",",es11.4)') stage, n, rs, t1 - t0
        t0 = t1
    end subroutine

    subroutine rs_case(n, r, lu)
        integer, intent(in) :: n, r, lu
        type(jellium_cell_t) :: cell, ideal_cell
        type(fock_hamiltonian_t) :: model, sector_model
        type(jellium_spectrum_t) :: spectrum, slater
        type(jellium_sector_t) :: sector
        type(mf_interaction_t) :: inter, ideal_inter
        type(mf_state_t) :: hartree, rhf, uhf, mermin, mermin_h, ideal, sp, sm
        type(mf_estimate_t) :: est
        integer, allocatable :: modes(:, :)
        real(dp), allocatable :: w(:), chis(:), gamma(:, :), gamma2(:, :), psi(:)
        real(dp) :: rs, t0, edges(max_bins + 1), ground, trunc, chi0, gc, rpa, rgc
        real(dp) :: chi, hmin, norm, norms(4), h
        integer :: ns, nc, nk, nexc, b, ns2, nc2
        logical :: ok

        rs = rs_values(r)
        call cpu_time(t0)
        call jellium_modes('spherical', cutoff, modes)
        call init_jellium_cell(cell, n, modes, rs, ok)
        if (.not. ok) error stop 'cell'
        call jellium_sector(cell, n/2, reshape([0, 0, 0], [3, 1]), sector, ok)
        call jellium_hamiltonian(cell, sector%basis, model, ok)
        if (.not. ok) error stop 'model'
        ! Exact sector reference (Q4 queries, unchanged).
        call jellium_build_spectrum(cell, model, nev, spectrum, ok)
        if (.not. ok) error stop 'spectrum'
        call cost(lu, 'exact_spectrum', n, rs, t0)
        allocate (w(spectrum%count), chis(spectrum%count), gamma(cell%m, cell%m), &
                  gamma2(cell%m*(cell%m - 1)/2, cell%m*(cell%m - 1)/2))
        edges(1:bins + 1) = [(0.5_dp*cell%length*b/bins, b=0, bins)]
        centers = 0.5_dp*(edges(1:bins) + edges(2:bins + 1))/cell%length
        call jellium_thermal_weights(cell, spectrum, 0.0_dp, w, trunc, ground)
        energy(1, r, n) = sum(w*spectrum%energy(1:spectrum%count))
        call jellium_thermal_rdms(cell, spectrum, w, gamma, gamma2, ok)
        call jellium_pair_bins(cell, gamma2, edges(1:bins + 1), pair(:, 1, r, n), &
                               pair(:, 2, r, n))
        call jellium_cumulants(cell, spectrum, w, cumulant(:, r, n), ok, 1.0e-8_dp)
        if (.not. ok) error stop 'exact cumulants'
        call jellium_static_response(cell, model, spectrum, w, q, response(1, r, n), &
                                     chis, ok, response_threshold)
        if (.not. ok) error stop 'exact response'
        call jellium_ideal_response(cell, 0.0_dp, q, chi0, gc, rpa, rgc, ok)
        response(4, r, n) = rpa; response(5, r, n) = chi0
        call jellium_thermal_weights(cell, spectrum, theta_thermal, w, trunc, ground)
        thermal(1, r, n) = sum(w*spectrum%energy(1:spectrum%count))
        thermal(9, r, n) = trunc
        call cost(lu, 'exact_observables', n, rs, t0)
        ! Mean-field models, one assumption at a time.
        call mf_interaction(cell, inter)
        call mf_ground(cell, inter, .false., .false., nstarts, hartree, ok, ns, nc)
        if (.not. ok) error stop 'Hartree'
        energy(2, r, n) = hartree%energy
        energy(3, r, n) = huge(1.0_dp)
        if (mod(n, 2) == 0) then
            call mf_ground(cell, inter, .true., .true., nstarts, rhf, ok, ns2, nc2)
            if (.not. ok) error stop 'RHF'
            energy(3, r, n) = rhf%energy
        end if
        call mf_ground(cell, inter, .false., .true., nstarts, uhf, ok, ns, nc)
        if (.not. ok) error stop 'UHF'
        call cost(lu, 'scf', n, rs, t0)
        energy(4, r, n) = uhf%energy; energy(5, r, n) = ns; energy(6, r, n) = nc
        energy(9, r, n) = uhf%nup
        call mf_static_response(cell, inter, uhf, q, chi, hmin, ok, nexc)
        if (.not. ok) chi = 0
        response(2, r, n) = chi; energy(8, r, n) = hmin
        call mf_static_response(cell, inter, hartree, q, chi, hmin, ok, nexc)
        if (.not. ok) chi = 0
        response(3, r, n) = chi
        call cost(lu, 'tdhf_response', n, rs, t0)
        ! Full-parent residual estimator and HF pair correlation.
        call mf_slater_sector(cell, uhf, sector, psi, norm, nk, ok)
        if (.not. ok) error stop 'Slater sector'
        energy(7, r, n) = nk
        call jellium_hamiltonian(cell, sector%basis, sector_model, ok)
        call mf_residual(sector_model, sector, psi, est, ok)
        if (.not. ok) error stop 'residual'
        estimate(:, r, n) = [est%rho, est%sigma, est%e0, est%r0, est%e1, est%beta, &
            est%temple, merge(1.0_dp, 0.0_dp, est%temple_valid), &
            real(est%dimension, dp), norm]
        call mf_as_spectrum(sector, psi, uhf%energy, slater)
        call jellium_thermal_rdms(cell, slater, [1.0_dp], gamma, gamma2, ok)
        call jellium_pair_bins(cell, gamma2, edges(1:bins + 1), pair(:, 3, r, n), &
                               pair(:, 4, r, n))
        call cost(lu, 'residual_estimator', n, rs, t0)
        ! Hellmann-Feynman derivative, finite-difference check (warm starts).
        derivative(1, r, n) = mf_coupling_derivative(cell, uhf)
        h = 1.0e-4_dp*cell%coupling
        call mf_interaction(cell, ideal_inter, cell%coupling + h)
        call mf_scf(cell, ideal_inter, uhf%nup, .false., .true., uhf%cup, uhf%cdn, sp, ok)
        call mf_interaction(cell, ideal_inter, cell%coupling - h)
        call mf_scf(cell, ideal_inter, uhf%nup, .false., .true., uhf%cup, uhf%cdn, sm, ok)
        derivative(2, r, n) = (sp%energy - sm%energy)/(2*h)
        if (n == 2) call exact_derivative(cell, model, derivative(3, r, n))
        if (n == 2) derivative(4, r, n) = derivative(3, r, n) - derivative(1, r, n)
        ! Mermin finite-T HF / Hartree and the ideal grand-canonical gas.
        call mf_mermin(cell, inter, theta_thermal, .true., mermin, ok)
        if (.not. ok) error stop 'Mermin HF'
        call mf_mermin(cell, inter, theta_thermal, .false., mermin_h, ok)
        if (.not. ok) error stop 'Mermin Hartree'
        call init_jellium_cell(ideal_cell, n, modes, rs, ok, 0.0_dp)
        call mf_interaction(ideal_cell, ideal_inter)
        call mf_mermin(ideal_cell, ideal_inter, theta_thermal, .true., ideal, ok)
        if (.not. ok) error stop 'ideal GC'
        thermal(2:8, r, n) = [mermin%energy, mermin%free, mermin%entropy, mermin%mu, &
                              mermin_h%energy, ideal%energy, real(mermin%iterations, dp)]
        call cost(lu, 'mermin', n, rs, t0)
        print '(a,i0,a,f5.2,a,es13.5,a,es13.5,a,i0,a,i0)', 'N=', n, ' rs=', rs, &
            ' E_exact=', energy(1, r, n), ' E_UHF=', uhf%energy, ' solutions=', ns, &
            ' K=', nk
    end subroutine

    !> Exact dE/dg of the nondegenerate N = 2 ground state (Q4 query).
    subroutine exact_derivative(cell, model, de)
        type(jellium_cell_t), intent(in) :: cell
        type(fock_hamiltonian_t), intent(in) :: model
        real(dp), intent(out) :: de
        type(jellium_sector_t) :: sector
        type(fock_hamiltonian_t) :: interaction
        real(dp) :: val(1), res(1), edges(2), gs(1), go(1), dgs(1), dgo(1)
        real(dp), allocatable :: vec(:, :)
        integer :: iterations
        logical :: ok
        call jellium_sector(cell, 1, reshape([0, 0, 0], [3, 1]), sector, ok)
        allocate (vec(sector%d, 1))
        call jellium_lowest(model, sector%basis, 1, val, vec, res, ok, iterations)
        if (.not. ok) error stop 'N=2 ground'
        call jellium_interaction_model(cell, sector%basis, interaction, ok)
        edges = [0.0_dp, 0.5_dp*cell%length]
        call jellium_coupling_derivative(cell, sector%basis, model, interaction, &
            vec(:, 1), val(1), edges, de, gs, go, dgs, dgo, ok, iterations)
        if (.not. ok) error stop 'exact derivative'
    end subroutine

    real(dp) function pz81(rs) result(ec)
        real(dp), intent(in) :: rs
        if (rs < 1) then
            call mf_pz81_high(rs, ec)
        else
            call mf_pz81_low(rs, ec)
        end if
    end function

    subroutine write_tables()
        integer :: out, n, r, b
        real(dp) :: rhf
        open (newunit=out, file='output/mean-field-energy.csv', status='replace')
        write (out, '(a)') 'N,rs,E_exact,E_hartree,E_RHF,E_UHF,Ec_UHF,Ec_RHF,'// &
            'Ec_UHF_per_N,pz81_ec_per_N_TL,uhf_nup,uhf_solutions,uhf_converged_runs,'// &
            'uhf_K_sectors,tdhf_hessian_min'
        do n = 2, 4
            do r = 1, nrs
                rhf = energy(3, r, n)
                if (mod(n, 2) /= 0) rhf = 0
                write (out, '(i0,",",f6.2,8(",",es16.8),4(",",i0),",",es16.8)') n, &
                    rs_values(r), &
                    energy(1:2, r, n), rhf, energy(4, r, n), &
                    energy(1, r, n) - energy(4, r, n), &
                    merge(energy(1, r, n) - rhf, 0.0_dp, mod(n, 2) == 0), &
                    (energy(1, r, n) - energy(4, r, n))/n, pz81(rs_values(r)), &
                    nint(energy(9, r, n)), nint(energy(5, r, n)), nint(energy(6, r, n)), &
                    nint(energy(7, r, n)), energy(8, r, n)
            end do
        end do
        close (out)
        open (newunit=out, file='output/mean-field-estimator.csv', status='replace')
        write (out, '(a)') 'N,rs,rho_E_UHF,sigma_residual,E0_sector,E0_ritz_residual,'// &
            'E1_sector,beta_gap_input,temple_lower,temple_valid,true_error,'// &
            'temple_error_bound,sector_dimension,slater_norm'
        do n = 2, 4
            do r = 1, nrs
                write (out, '(i0,",",f6.2,7(",",es16.8),",",i0,2(",",es16.8),",",i0,'// &
                    '",",es16.8)') n, rs_values(r), estimate(1:7, r, n), &
                    nint(estimate(8, r, n)), estimate(1, r, n) - estimate(3, r, n), &
                    merge(estimate(1, r, n) - estimate(7, r, n), -1.0_dp, &
                          estimate(8, r, n) > 0.5_dp), nint(estimate(9, r, n)), &
                    estimate(10, r, n)
            end do
        end do
        close (out)
        open (newunit=out, file='output/mean-field-response.csv', status='replace')
        write (out, '(a)') 'N,rs,q_n,chi_exact,chi_TDHF,chi_TDH_RPA,'// &
            'chi_RPA_canonical,chi0_canonical_lindhard,tdhf_hessian_min'
        do n = 2, 4
            do r = 1, nrs
                write (out, '(i0,",",f6.2,",",3i1,6(",",es16.8))') n, rs_values(r), q, &
                    response(:, r, n), energy(8, r, n)
            end do
        end do
        close (out)
        open (newunit=out, file='output/mean-field-thermal.csv', status='replace')
        write (out, '(a)') 'N,rs,theta,U_exact_canonical,U_mermin_HF,F_mermin_HF,'// &
            'S_mermin_HF,mu_mermin_HF,U_mermin_hartree,U_ideal_grand_canonical,'// &
            'mermin_iterations,exact_truncation_weight'
        do n = 2, 4
            do r = 1, nrs
                write (out, '(i0,",",f6.2,",",f5.2,9(",",es16.8))') n, rs_values(r), &
                    theta_thermal, thermal(:, r, n)
            end do
        end do
        close (out)
        open (newunit=out, file='output/mean-field-cumulants.csv', status='replace')
        write (out, '(a)') 'N,rs,norm_Gamma2,norm_lambda2_exact,norm_Gamma3,'// &
            'norm_lambda3_exact,lambda2_HF,Ec_UHF'
        do n = 2, 4
            do r = 1, nrs
                write (out, '(i0,",",f6.2,6(",",es16.8))') n, rs_values(r), &
                    cumulant(:, r, n), 0.0_dp, energy(1, r, n) - energy(4, r, n)
            end do
        end do
        close (out)
        open (newunit=out, file='output/mean-field-derivative.csv', status='replace')
        write (out, '(a)') 'N,rs,dE_UHF_dg_hellmann_feynman,dE_UHF_dg_central_FD,'// &
            'dE_exact_dg_N2_only,dEc_dg_N2_only'
        do n = 2, 4
            do r = 1, nrs
                write (out, '(i0,",",f6.2,4(",",es18.10))') n, rs_values(r), &
                    derivative(:, r, n)
            end do
        end do
        close (out)
        open (newunit=out, file='output/mean-field-pair-correlation.csv', status='replace')
        write (out, '(a)') 'N,rs,r_over_L,g_same_exact,g_opp_exact,g_same_UHF,g_opp_UHF'
        do n = 2, 4
            do r = 1, nrs
                do b = 1, bins
                    write (out, '(i0,",",f6.2,",",f8.5,4(",",es15.7))') n, &
                        rs_values(r), centers(b), pair(b, :, r, n)
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
            'neutralizing background, spherical plane-wave cutoff; Hartree, RHF, '// &
            'UHF, Mermin HF and static TDHF in the same spin-orbital basis')
        call result%attribute('status', 'finite-basis candidates: exact = Q4 '// &
            'sector Lanczos (ground multiplet average at T=0); mean-field = lowest '// &
            'converged SCF of several starts; Temple bound uses a candidate gap')
        call result%write('rs', rs_values(1:nrs), 'a0', 'Wigner-Seitz radius', dims='rs')
        call result%write('r_over_L', centers, '1', 'radial shell center / L', dims='bin')
        do n = 2, 4
            write (tag, '(a,i0)') '_N', n
            call result%write('energy'//trim(tag), energy(:, :, n), 'Ha', &
                'exact,Hartree,RHF,UHF,solutions,runs,K sectors,Hessian min,nup', &
                dims='energy_item,rs', status='candidate')
            call result%write('estimator'//trim(tag), estimate(:, :, n), 'Ha', &
                'rho,sigma,E0,r0,E1,beta,Temple,valid,dimension,norm', &
                dims='estimator_item,rs', status='Temple: candidate gap input')
            call result%write('chi'//trim(tag), response(:, :, n), 'a0^-3 Ha^-1', &
                'exact,TDHF,TDH(RPA),RPA canonical,Lindhard canonical', &
                dims='response_item,rs', status='candidate')
            call result%write('thermal'//trim(tag), thermal(:, :, n), 'Ha', &
                'U exact,U MHF,F MHF,S,mu,U M-Hartree,U ideal GC,iterations,'// &
                'truncation', dims='thermal_item,rs', status='candidate')
            call result%write('cumulant_norms'//trim(tag), cumulant(:, :, n), '1', &
                '|Gamma2|,|lambda2|,|Gamma3|,|lambda3| (exact; HF = 0)', &
                dims='norm,rs')
            call result%write('coupling_derivative'//trim(tag), derivative(:, :, n), &
                'Ha', 'dE_HF/dg HF, FD, dE_exact/dg, dEc/dg', dims='derivative,rs')
            call result%write('pair'//trim(tag), pair(:, :, :, n), '1', &
                'g_same exact, g_opp exact, g_same UHF, g_opp UHF', &
                dims='bin,pair_item,rs', status='candidate')
        end do
        call result%close()
    end subroutine

    subroutine make_figures()
        real(dp) :: palette(3, 3), x(nrs), y(nrs), z(nrs), tl(200), ptl(200)
        character(48) :: label
        integer :: n, r, k
        palette(:, 1) = blue; palette(:, 2) = vermillion; palette(:, 3) = green
        x = rs_values(1:nrs)
        do k = 1, 200
            tl(k) = exp(log(x(1)) + (log(x(nrs)) - log(x(1)))*(k - 1)/199.0_dp)
            ptl(k) = -pz81(tl(k))
        end do
        ! Correlation energy and the missing cumulant.
        call figure(figsize=[12.8_dp, 6.4_dp])
        call subplot(1, 2, 1)
        call title('Correlation energy $E_c = E_{exact} - E_{UHF}$')
        call xlabel('$r_s$'); call ylabel('$-E_c / N$ [Ha]')
        call xscale('log'); call yscale('log')
        do n = 2, 4
            y = -(energy(1, 1:nrs, n) - energy(4, 1:nrs, n))/n
            write (label, '(a,i0,a)') 'cell, N = ', n, ' (UHF reference)'
            call plot(x, y, color=palette(:, n - 1), marker='o', label=trim(label))
        end do
        y = -(energy(1, 1:nrs, 4) - energy(3, 1:nrs, 4))/4
        call plot(x, y, color=palette(:, 3), marker='s', linestyle=': ', &
                  label='cell, N = 4 (RHF reference)')
        call plot(tl, ptl, color=black, linestyle='--', &
                  label='Perdew-Zunger fit to Ceperley-Alder, N = ∞')
        call ylim(2.0e-3_dp, 0.5_dp)
        call legend(loc='upper right')
        call subplot(1, 2, 2)
        call title('Exact pair cumulant (zero for every Slater state)')
        call xlabel('$r_s$'); call ylabel('$|\lambda_2| / |\Gamma_2|$')
        call xscale('log')
        do n = 2, 4
            y = cumulant(2, 1:nrs, n)/cumulant(1, 1:nrs, n)
            write (label, '(a,i0)') 'exact, N = ', n
            call plot(x, y, color=palette(:, n - 1), marker='o', label=trim(label))
        end do
        call plot([x(1), x(nrs)], [0.0_dp, 0.0_dp], color=black, linestyle='--', &
                  label='HF / Hartree')
        call legend(loc='upper left')
        call suptitle('Mean field versus the exact jellium cell (same basis); '// &
                      'N = 3, 4 exact: ground-multiplet average')
        call tight_layout()
        call savefig('output/mean-field-correlation.png')
        ! Pair correlation: Coulomb hole missing in HF.
        call figure(figsize=[12.8_dp, 6.4_dp])
        do k = 1, 2
            call subplot(1, 2, k)
            n = merge(2, 4, k == 1)
            if (k == 1) call title('Opposite-spin pair correlation, N = 2')
            if (k == 2) call title('Same-spin pair correlation, N = 4')
            call xlabel('r / L'); call ylabel('g(r)')
            call xlim(0.0_dp, 0.5_dp); call ylim(0.0_dp, 1.5_dp)
            do r = 1, nrs, max(1, (nrs - 1)/2)
                write (label, '(a,f4.1)') 'exact, $r_s$ = ', rs_values(r)
                call plot(centers, pair(:, 3 - k, r, n), color=palette(:, 1 + mod(r - 1, 3)), &
                          marker='o', label=trim(label))
                write (label, '(a,f4.1)') 'UHF, $r_s$ = ', rs_values(r)
                call plot(centers, pair(:, 5 - k, r, n), &
                          color=palette(:, 1 + mod(r - 1, 3)), linestyle='--', &
                          label=trim(label))
            end do
            call legend(loc='lower right')
        end do
        call suptitle('HF keeps the exchange hole but has no Coulomb hole')
        call tight_layout()
        call savefig('output/mean-field-pair-correlation.png')
        ! Static response hierarchy relative to the same-cell Lindhard function.
        call figure(figsize=[12.8_dp, 6.4_dp])
        do k = 1, 2
            call subplot(1, 2, k)
            n = merge(2, 4, k == 1)
            write (label, '(a,i0)') 'Static response, q L/2π = (100), N = ', n
            call title(trim(label))
            call xlabel('$r_s$'); call ylabel('$\chi / \chi_0$')
            call xscale('log')
            y = response(1, 1:nrs, n)/response(5, 1:nrs, n)
            call plot(x, y, color=blue, marker='o', label='exact')
            y = response(2, 1:nrs, n)/response(5, 1:nrs, n)
            call plot(x, y, color=vermillion, marker='s', label='TDHF (RPA + exchange)')
            y = response(3, 1:nrs, n)/response(5, 1:nrs, n)
            call plot(x, y, color=green, marker='^', label='TDH = RPA (mean-field state)')
            y = response(4, 1:nrs, n)/response(5, 1:nrs, n)
            call plot(x, y, color=orange, linestyle=': ', label='RPA (canonical Lindhard)')
            call plot([x(1), x(nrs)], [1.0_dp, 1.0_dp], color=black, linestyle='--', &
                      label='Lindhard (same cell)')
            call legend(loc='lower left')
        end do
        call suptitle('Exchange in the kernel and orbital energies versus correlation')
        call tight_layout()
        call savefig('output/mean-field-response.png')
        ! Full-parent residual estimator versus the true HF energy error.
        call figure(figsize=[12.8_dp, 6.4_dp])
        call subplot(1, 2, 1)
        call title('HF energy error in its symmetry sector')
        call xlabel('$r_s$'); call ylabel('error [Ha]')
        call xscale('log'); call yscale('log')
        do n = 2, 4
            y = estimate(1, 1:nrs, n) - estimate(3, 1:nrs, n)
            write (label, '(a,i0)') 'true $E_{HF} - E_0$, N = ', n
            call plot(x, y, color=palette(:, n - 1), marker='o', label=trim(label))
            z = merge(estimate(1, 1:nrs, n) - estimate(7, 1:nrs, n), y, &
                      estimate(8, 1:nrs, n) > 0.5_dp)
            write (label, '(a,i0)') 'Temple bound, N = ', n
            call plot(x, z, color=palette(:, n - 1), linestyle='--', marker='s', &
                      label=trim(label))
        end do
        call legend(loc='lower left')
        call subplot(1, 2, 2)
        call title('Estimator efficiency')
        call xlabel('$r_s$'); call ylabel('bound / true error')
        call xscale('log')
        do n = 2, 4
            y = (estimate(1, 1:nrs, n) - estimate(7, 1:nrs, n))/ &
                (estimate(1, 1:nrs, n) - estimate(3, 1:nrs, n))
            write (label, '(a,i0)') 'Temple, N = ', n
            call plot(x, y, color=palette(:, n - 1), marker='o', label=trim(label))
        end do
        call plot([x(1), x(nrs)], [1.0_dp, 1.0_dp], color=black, linestyle='--', &
                  label='exact error')
        call legend(loc='upper left')
        call suptitle('Residual $||(H - E_{HF})\Phi_{HF}||$ with Temple; gap input '// &
                      'from Lanczos (candidate)')
        call tight_layout()
        call savefig('output/mean-field-estimator.png')
        if (.false.) print *, purple, sky
    end subroutine
end program mean_field_cell_demo
