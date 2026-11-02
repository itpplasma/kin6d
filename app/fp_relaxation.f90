program fp_relaxation_main
    ! B04 homogeneous relaxation of the drifting W3+ test population. Two
    ! resolutions (cells and time step both doubled) give the error estimate.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_fp_relaxation, only: fp_relax_t, fp_relax_setup, fp_relax_initial, &
        fp_relax_run, fp_relax_cells, fp_relax_nmom
    use kin6d_result_io, only: run_t, result_file_t, runtime_git
    implicit none
    integer, parameter :: ntimes = 12
    real(dp), parameter :: amu = 1.66053906660e-27_dp, e_charge = 1.602176634e-19_dp
    real(dp), parameter :: eps0 = 8.8541878128e-12_dp, mp = 1.67262192369e-27_dp
    type(run_t) :: run
    type(result_file_t) :: result
    type(fp_relax_t) :: coarse, fine
    integer :: cells_per_unit = 64, ncoarse, nfine, unit, ios, it, k
    real(dp) :: dt_tau = 0.01_dp, tail = 1.0e-14_dp
    real(dp) :: ma_amu = 183.84_dp, tb_ev = 50.0_dp, tpar_ev = 200.0_dp, tperp_ev = 50.0_dp
    real(dp) :: drift = 7.2e3_dp, dens = 1.0_dp, nb = 1.0e20_dp, lnl = 15.0_dp
    real(dp) :: output_times(ntimes)
    real(dp) :: mom_c(fp_relax_nmom, ntimes), mom_f(fp_relax_nmom, ntimes)
    real(dp) :: err(fp_relax_nmom, ntimes), row(11), tau_s
    character(len=256) :: case_file = '', case_sha256 = 'not supplied'
    character(len=:), allocatable :: commit, tree
    character(len=256) :: message
    character(len=*), parameter :: names(fp_relax_nmom) = &
        [character(len=8) :: 'n', 'u_z', 'p_par', 'p_perp', 'v4']
    character(len=*), parameter :: units(fp_relax_nmom) = &
        [character(len=8) :: 'm^-3', 'm/s', 'Pa', 'Pa', 'm^4/s^4']
    namelist /fp_relaxation/ cells_per_unit, dt_tau, tail, ma_amu, tb_ev, tpar_ev, &
        tperp_ev, drift, dens, nb, lnl, output_times, case_file, case_sha256
    output_times = [0.0_dp, 0.01_dp, 0.03_dp, 0.1_dp, 0.3_dp, 1.0_dp, 3.0_dp, &
        10.0_dp, 30.0_dp, 100.0_dp, 300.0_dp, 1000.0_dp]

    call run%start('fp_relaxation')
    ios = 0; message = ''
    call run%open_config(unit)
    if (run%has_config()) read (unit, nml=fp_relaxation, iostat=ios, iomsg=message)
    call run%close_config(unit, 'fp_relaxation', ios, message)
    if (run%argument_count() /= 0) error stop 'usage: fp_relaxation [--config FILE]'
    if (cells_per_unit < 8) error stop 'fp_relaxation needs at least eight cells per unit'
    if (.not. (dt_tau > 0.0_dp .and. tail > 0.0_dp)) error stop 'invalid time step or tail'
    if (any(output_times(2:) <= output_times(:ntimes - 1))) error stop 'output times must increase'

    call run%open_resolved(unit)
    write (unit, nml=fp_relaxation)
    call run%close_resolved(unit)

    ncoarse = fp_relax_cells(cells_per_unit, ma_amu*amu, tb_ev*e_charge, drift, &
        tpar_ev*e_charge, tail)
    nfine = fp_relax_cells(2*cells_per_unit, ma_amu*amu, tb_ev*e_charge, drift, &
        tpar_ev*e_charge, tail)
    call solve(cells_per_unit, dt_tau, mom_c)
    call solve(2*cells_per_unit, dt_tau/2.0_dp, mom_f)
    err = abs(mom_c - mom_f)
    call tau_of_b04(tau_s)

    call execute_command_line('mkdir -p output/homma')
    call runtime_git(commit, tree, 'fp_relaxation')
    open (newunit=unit, file='output/homma/B04.csv', status='replace', action='write')
    write (unit, '(a)') '# code: kin6d (fp_relaxation, Eulerian Legendre-mode FP)'
    write (unit, '(a)') '# code_revision: '//commit
    write (unit, '(a)') '# case_file: '//trim(case_file)
    write (unit, '(a)') '# case_sha256: '//trim(case_sha256)
    write (unit, '(a)') '# convention_deviations: n normalised to 1 m^-3 at t=0 (equation is linear in f); '// &
        'tau_s uses v_a = sqrt(3 T_b/m_a) as the case file states (conventions.md default is T_a); '// &
        'results written to output/homma/B04.csv (not results/<code>/); '// &
        'no electron or field effect (B = 0 in the moments); heat flux zero'
    write (unit, '(a)') '# precision: double (IEEE binary64)'
    write (unit, '(a,i0,a,i0,a,es10.3,a,es10.3)') '# resolution: cells per u_a=', cells_per_unit, &
        ' (fine: ', 2*cells_per_unit, ') dt_tau=', dt_tau, ' (fine: ', dt_tau/2.0_dp
    write (unit, '(a,es25.16e3)') '# tau_s_seconds: ', tau_s
    write (unit, '(a)') '# units: t_over_tau_s 1; n m^-3; u_z m/s; p_par Pa; p_perp Pa; v4 m^4/s^4; '// &
        'err_* in the same units as the moment, |coarse - fine| (two-resolution estimate)'
    write (unit, '(a)') '# status: candidate'
    write (unit, '(a)') 't_over_tau_s,n,u_z,p_par,p_perp,v4,err_n,err_u_z,err_p_par,err_p_perp,err_v4,status'
    do it = 1, ntimes
        row(1) = output_times(it)
        row(2:6) = mom_f(:, it)
        row(7:11) = err(:, it)
        write (unit, '(11(es25.16e3,","),a)') row, 'candidate'
    end do
    close (unit)

    call execute_command_line('mkdir -p output')
    call result%create(run)
    call result%write('time_over_tau_s', output_times, '1', 'time in tau_s', dims='time', &
        status='deterministic')
    do k = 1, fp_relax_nmom
        call result%write(trim(names(k)), mom_f(k, :), units(k), 'fine-resolution moment', &
            dims='time', status='candidate')
        call result%write('err_'//trim(names(k)), err(k, :), units(k), &
            'two-resolution estimate |coarse-fine|', dims='time', status='candidate')
    end do
    call result%close()

    print '(a,i0,a,i0)', 'cells coarse=', ncoarse, ' fine=', nfine
    do it = 1, ntimes
        print '(f9.3,5es14.6e3)', output_times(it), mom_f(:, it)
    end do
contains

    subroutine solve(cpu, dt, mom)
        integer, intent(in) :: cpu
        real(dp), intent(in) :: dt
        real(dp), intent(out) :: mom(fp_relax_nmom, ntimes)
        type(fp_relax_t) :: s
        integer :: n

        n = fp_relax_cells(cpu, ma_amu*amu, tb_ev*e_charge, drift, tpar_ev*e_charge, tail)
        call fp_relax_setup(s, n, ma_amu*amu, mp, 3.0_dp*e_charge, e_charge, nb, &
            tb_ev*e_charge, lnl, eps0, tpar_ev*e_charge, tperp_ev*e_charge, &
            drift, dens, tail)
        call fp_relax_initial(s)
        call fp_relax_run(s, output_times, dt, mom)
    end subroutine solve

    subroutine tau_of_b04(tau)
        ! tau_s at v_a = sqrt(3 T_b/m_a) (case file), from the leaves.
        use kin6d_fp_relaxation_leaves, only: fp_rel_gamma, fp_rel_speed, fp_rel_xs, &
            fp_rel_mu, fp_rel_tau_s
        real(dp), intent(out) :: tau
        real(dp) :: gam, vt, va, x, muf, mr

        mr = ma_amu*amu/mp
        call fp_rel_gamma(3.0_dp*e_charge, e_charge, lnl, eps0, ma_amu*amu, gam)
        call fp_rel_speed(tb_ev*e_charge, mp, 2.0_dp, vt)
        call fp_rel_speed(tb_ev*e_charge, ma_amu*amu, 3.0_dp, va)
        call fp_rel_xs(va, vt, x)
        call fp_rel_mu(x, muf)
        call fp_rel_tau_s(va, muf, gam, nb, mr, tau)
    end subroutine tau_of_b04

end program fp_relaxation_main
