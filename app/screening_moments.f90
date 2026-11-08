program screening_moments_main
    ! Thermal-force benchmarks B05-B07 from the first-order spatial-moment equations
    ! (milca_screening_moments). Two resolutions (cells, oracle panels and window
    ! steps all doubled) give the error estimate |coarse - fine|.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use milca_screening_moments, only: scr_t, scr_setup, scr_steady, scr_rate0, scr_transient
    use milca_fp_relaxation, only: fp_relax_cells
    use milca_rosenbluth_reference, only: homma_pair_t, homma_closure
    use milca_result_io, only: run_t, result_file_t, runtime_git
    implicit none
    real(dp), parameter :: amu = 1.66053906660e-27_dp, e_charge = 1.602176634e-19_dp
    real(dp), parameter :: mp = 1.67262192369e-27_dp
    integer, parameter :: nb07 = 4, nwin = 5, nsample = 4
    character(len=*), parameter :: code = 'milca screening_moments (Eulerian, first order in grad T)'
    type(run_t) :: run
    type(result_file_t) :: result
    integer :: cells_per_unit = 64, window_steps = 100, unit, ios, k
    real(dp) :: ppu = 4.0_dp, tail = 1.0e-14_dp, ma_amu = 183.84_dp, t0_ev = 50.0_dp
    real(dp) :: nb = 1.0e20_dp, lnl = 15.0_dp, b_tesla = 1.0_dp
    integer :: window_charge = 3, b07_charge(nb07) = [3, 1, 10, 3]
    real(dp) :: b07_grad_ev_m(nb07) = [100.0_dp, 100.0_dp, 100.0_dp, 300.0_dp]
    real(dp) :: window_grad_ev_m(3, nwin)
    character(len=3) :: window_case(nwin) = ['B05', 'B05', 'B05', 'B06', 'B06']
    character(len=256) :: case_b05 = '', case_b06 = '', case_b07 = ''
    character(len=256) :: sha_b05 = 'not supplied', sha_b06 = 'not supplied', &
        sha_b07 = 'not supplied'
    character(len=256) :: message
    character(len=:), allocatable :: commit, tree
    ! B07: (variant full/gc, component X/Y, row); window: rate and samples per component.
    real(dp) :: drift(2, 2, nb07), drift_err(2, 2, nb07)
    real(dp) :: rate(3, nwin), rate_err(3, nwin)
    real(dp) :: sample(3, nsample, nwin), sample_err(3, nsample, nwin)
    namelist /screening_moments/ cells_per_unit, ppu, window_steps, tail, ma_amu, t0_ev, nb, &
        lnl, b_tesla, window_charge, b07_charge, b07_grad_ev_m, window_grad_ev_m, &
        case_b05, sha_b05, case_b06, sha_b06, case_b07, sha_b07
    window_grad_ev_m = reshape([0.0_dp, 0.0_dp, 5.0_dp, 0.0_dp, 0.0_dp, -5.0_dp, &
        300.0_dp, 0.0_dp, 5.0_dp, 100.0_dp, 0.0_dp, 0.0_dp, 300.0_dp, 0.0_dp, 0.0_dp], [3, nwin])

    call run%start('screening_moments')
    ios = 0; message = ''
    call run%open_config(unit)
    if (run%has_config()) read (unit, nml=screening_moments, iostat=ios, iomsg=message)
    call run%close_config(unit, 'screening_moments', ios, message)
    if (run%argument_count() /= 0) error stop 'usage: screening_moments [--config FILE]'
    if (cells_per_unit < 8) error stop 'screening_moments needs at least eight cells per unit'
    if (.not. (ppu > 0.0_dp .and. tail > 0.0_dp .and. b_tesla > 0.0_dp)) &
        error stop 'invalid ppu, tail or field'
    if (window_steps < nsample .or. mod(window_steps, nsample) /= 0) &
        error stop 'window_steps must be a positive multiple of four'
    call run%open_resolved(unit)
    write (unit, nml=screening_moments)
    call run%close_resolved(unit)

    call b07_drifts(cells_per_unit, ppu, drift)
    call b07_drifts(2*cells_per_unit, 2.0_dp*ppu, drift_err)
    drift_err = abs(drift_err - drift)
    call b07_drifts(2*cells_per_unit, 2.0_dp*ppu, drift)
    call window(cells_per_unit, ppu, window_steps, rate, sample)
    call window(2*cells_per_unit, 2.0_dp*ppu, 2*window_steps, rate_err, sample_err)
    rate_err = abs(rate_err - rate)
    sample_err = abs(sample_err - sample)
    call window(2*cells_per_unit, 2.0_dp*ppu, 2*window_steps, rate, sample)

    call execute_command_line('mkdir -p output/homma')
    call runtime_git(commit, tree, 'screening_moments')
    call write_b07()
    call write_window('B05', case_b05, sha_b05)
    call write_window('B06', case_b06, sha_b06)

    call result%create(run)
    call result%write('b07_charge_number', real(b07_charge, dp), '1', 'test-ion charge number', &
        dims='b07_case', status='deterministic')
    call result%write('b07_grad_t', b07_grad_ev_m, 'eV/m', 'grad_x T_b', dims='b07_case', &
        status='deterministic')
    call result%write('b07_drift', drift, 'm/s', &
        'steady d<X,Y>/dt (variant full, gc; component X, Y; case)', &
        dims='variant,component,b07_case', status='candidate')
    call result%write('b07_drift_err', drift_err, 'm/s', 'two-resolution estimate |coarse-fine|', &
        dims='variant,component,b07_case', status='candidate')
    call result%write('window_grad_t', window_grad_ev_m, 'eV/m', 'grad T_b of B05/B06 cases', &
        dims='xyz,window_case', status='deterministic')
    call result%write('window_rate0', rate, 'm/s^2', 'd<v>/dt at t = 0', &
        dims='xyz,window_case', status='candidate')
    call result%write('window_rate0_err', rate_err, 'm/s^2', 'two-resolution estimate', &
        dims='xyz,window_case', status='candidate')
    call result%write('window_mean_v', sample, 'm/s', '<v> at Omega_a t = 0.0025 k', &
        dims='xyz,sample,window_case', status='candidate')
    call result%write('window_mean_v_err', sample_err, 'm/s', 'two-resolution estimate', &
        dims='xyz,sample,window_case', status='candidate')
    call result%close()

    do k = 1, nb07
        print '(a,i0,a,f6.1,a,2es14.6,a,2es14.6)', 'B07 Z=', b07_charge(k), ' grad=', &
            b07_grad_ev_m(k), ' full X,Y=', drift(1, :, k), ' gc X,Y=', drift(2, :, k)
    end do
    do k = 1, nwin
        print '(a,1x,3f8.1,a,3es14.6)', window_case(k), window_grad_ev_m(:, k), &
            ' d<v>/dt(0)=', rate(:, k)
    end do

contains

    subroutine species(z, s, cpu, ppu_in, omega)
        integer, intent(in) :: z, cpu
        real(dp), intent(in) :: ppu_in
        type(scr_t), intent(out) :: s
        real(dp), intent(out) :: omega
        real(dp) :: ma, t0
        integer :: n

        ma = ma_amu*amu
        t0 = t0_ev*e_charge
        n = fp_relax_cells(cpu, ma, t0, 0.0_dp, t0, tail)
        call scr_setup(s, n, ma, mp, z*e_charge, e_charge, nb, t0, lnl, tail, ppu_in)
        omega = z*e_charge*b_tesla/ma
    end subroutine species

    ! homma2014 closure at T_0 with b = e_z.
    function heat_flux(grad_ev_m) result(q)
        real(dp), intent(in) :: grad_ev_m(3)
        real(dp) :: q(3)
        type(homma_pair_t) :: pair

        pair%n_b = nb; pair%t_b = t0_ev*e_charge; pair%m_b = mp; pair%q_b = e_charge
        pair%lnlambda = lnl
        call homma_closure(pair, [0.0_dp, 0.0_dp, 1.0_dp], b_tesla, grad_ev_m*e_charge, q)
    end function heat_flux

    subroutine b07_drifts(cpu, ppu_in, d)
        integer, intent(in) :: cpu
        real(dp), intent(in) :: ppu_in
        real(dp), intent(out) :: d(2, 2, nb07)
        type(scr_t) :: s
        real(dp) :: omega, g(3), v(3)
        integer :: k, variant

        do k = 1, nb07
            call species(b07_charge(k), s, cpu, ppu_in, omega)
            g = [b07_grad_ev_m(k), 0.0_dp, 0.0_dp]
            do variant = 1, 2
                call scr_steady(s, omega, g*e_charge, heat_flux(g), variant == 2, v)
                d(variant, :, k) = v(1:2)
            end do
        end do
    end subroutine b07_drifts

    ! Full variant (particles start at the origin, collisions at the particle).
    subroutine window(cpu, ppu_in, nsteps, r, smp)
        integer, intent(in) :: cpu, nsteps
        real(dp), intent(in) :: ppu_in
        real(dp), intent(out) :: r(3, nwin), smp(3, nsample, nwin)
        type(scr_t) :: s
        real(dp) :: omega, q(3), g(3)
        real(dp), allocatable :: vel(:, :)
        integer :: k, j

        call species(window_charge, s, cpu, ppu_in, omega)
        allocate (vel(3, 0:nsteps))
        do k = 1, nwin
            g = window_grad_ev_m(:, k)
            q = heat_flux(g)
            call scr_rate0(s, omega, g*e_charge, q, .false., r(:, k))
            call scr_transient(s, omega, g*e_charge, q, .false., 0.01_dp/omega, nsteps, vel)
            do j = 1, nsample
                smp(:, j, k) = vel(:, j*nsteps/nsample)
            end do
        end do
    end subroutine window

    subroutine header(case_id, case_file, sha)
        character(len=*), intent(in) :: case_id, case_file, sha
        write (unit, '(a)') '# code: '//code
        write (unit, '(a)') '# code_revision: '//commit//trim(merge(' (dirty)', '        ', tree == 'dirty'))
        write (unit, '(a)') '# case_file: '//trim(case_file)
        write (unit, '(a)') '# case_file_sha256: '//trim(sha)
        write (unit, '(a)') '# case_id: '//case_id
        write (unit, '(a)') '# precision: double (IEEE binary64)'
        write (unit, '(a,i0,a,i0,a,f0.1,a,f0.1,a,f0.1,a)') '# resolution: radial cells per '// &
            'sqrt(2 T_0/m_a) = ', cells_per_unit, ' (fine ', 2*cells_per_unit, &
            '), Rosenbluth-oracle panels per v_t,b = ', ppu, ' (fine ', 2.0_dp*ppu, &
            '), oracle range 8 v_t,b, test-Maxwellian tail ', log10(tail), ' decades'
    end subroutine header

    subroutine write_b07()
        character(len=11), parameter :: mode(3) = ['particle   ', 'gc         ', 'particle-gc']
        real(dp) :: val(2), err(2)
        integer :: k, variant

        open (newunit=unit, file='output/homma/B07.csv', status='replace', action='write')
        call header('B07', case_b07, sha_b07)
        write (unit, '(a)') '# quantity: steady d<X>/dt, d<Y>/dt of the first-order moment '// &
            'equations (the t -> infinity slope; d<X>/dt = d<x>/dt there); mode particle = '// &
            'collisions at T_b(x), gc = collisions at T_b(X), X = x - b x v/Omega_a, '// &
            'particle-gc = their difference'
        write (unit, '(a)') '# convention_deviations: first order in grad T (no O(grad T^2) '// &
            'cloud spreading, so the 300 eV/m rows are exactly three times the 100 eV/m rows); '// &
            'steady slope instead of a fit over 100 <= Omega_a t <= 700; deterministic, so se '// &
            'columns are 0 and dt_coll_over_tau_s is the label converged; '// &
            'results written to output/homma/ (copy to results/milca/)'
        write (unit, '(a)') '# error_estimate: err_dt_* = |coarse - fine| over the two '// &
            'resolutions above; values are the fine resolution'
        write (unit, '(a)') '# units: drifts m/s along e_x and e_y; grad_x eV/m'
        write (unit, '(a)') '# status: candidate'
        write (unit, '(a)') 'charge_number,grad_x_eV_per_m,mode,estimate,dt_coll_over_tau_s,'// &
            'd_X_dt_m_s,se_d_X_dt,err_dt_d_X_dt,d_Y_dt_m_s,se_d_Y_dt,err_dt_d_Y_dt,status'
        do k = 1, nb07
            do variant = 1, 3
                if (variant < 3) then
                    val = drift(variant, :, k)
                    err = drift_err(variant, :, k)
                else
                    val = drift(1, :, k) - drift(2, :, k)
                    err = hypot(drift_err(1, :, k), drift_err(2, :, k))
                end if
                write (unit, '(i0,",",a,",",a,",fit,converged,",es24.16e3,",0,",es10.3e3,'// &
                    '",",es24.16e3,",0,",es10.3e3,",candidate")') b07_charge(k), &
                    decimal(b07_grad_ev_m(k)), trim(mode(variant)), val(1), err(1), val(2), err(2)
            end do
        end do
        close (unit)
    end subroutine write_b07

    subroutine write_window(case_id, case_file, sha)
        character(len=*), intent(in) :: case_id, case_file, sha
        character(len=1), parameter :: comp(3) = ['x', 'y', 'z']
        integer :: k, c

        open (newunit=unit, file='output/homma/'//case_id//'.csv', status='replace', action='write')
        call header(case_id, case_file, sha)
        write (unit, '(a,i0,a)') '# quantity: d<v>/dt at t = 0 of the semi-discrete moment '// &
            'equations (no window fit, no control run needed: the equations are linear in '// &
            'grad T) and <v>(t) at Omega_a t = 0.0025, 0.005, 0.0075, 0.01 by TR-BDF2 with ', &
            window_steps, ' steps over the window (fine: twice as many); W3+, particles start '// &
            'at the origin, collisions at the particle position'
        write (unit, '(a)') '# convention_deviations: first order in grad T; deterministic, '// &
            'so se_window is 0; d_mean_v_dt_window_m_s2 holds the t = 0 rate; '// &
            'results written to output/homma/ (copy to results/milca/)'
        write (unit, '(a)') '# error_estimate: err_* = |coarse - fine| over the two resolutions; '// &
            'values are the fine resolution'
        write (unit, '(a)') '# units: rates m/s^2; velocities m/s; grad eV/m'
        write (unit, '(a)') '# status: candidate'
        write (unit, '(a)') 'case,grad_eV_per_m,component,d_mean_v_dt_window_m_s2,se_window,'// &
            'err_d_mean_v_dt,v_at_Omega_t_0p0025,v_at_Omega_t_0p005,v_at_Omega_t_0p0075,'// &
            'v_at_Omega_t_0p01,err_v_at_Omega_t_0p01,status'
        do k = 1, nwin
            if (window_case(k) /= case_id) cycle
            do c = 1, 3
                write (unit, '(a,",""(",a,", ",a,", ",a,")"",",a,",",es24.16e3,'// &
                    '",0,",es10.3e3,4(",",es24.16e3),",",es10.3e3,",candidate")') case_id, &
                    decimal(window_grad_ev_m(1, k)), decimal(window_grad_ev_m(2, k)), &
                    decimal(window_grad_ev_m(3, k)), comp(c), rate(c, k), rate_err(c, k), &
                    sample(c, :, k), sample_err(c, nsample, k)
            end do
        end do
        close (unit)
    end subroutine write_window

    ! One decimal with a leading zero, as the reference files print gradients.
    function decimal(x) result(text)
        real(dp), intent(in) :: x
        character(len=:), allocatable :: text
        character(len=32) :: buffer
        write (buffer, '(f0.1)') x
        text = trim(buffer)
        if (text(1:1) == '.') text = '0'//text
        if (text(1:2) == '-.') text = '-0'//text(2:)
    end function decimal

end program screening_moments_main
