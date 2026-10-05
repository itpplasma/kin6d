program benchmark_fermi_gas
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_fermi_gas, only: fermi_workspace_t, fermi_state_t, fermi_si_state_t, &
        initialize_fermi_workspace, solve_fermi_gas, solve_fermi_gas_si
    use kin6d_fermi_gas_leaves, only: fermi_electron_density, &
                                    fermi_relativity_indicator
    implicit none
    integer, parameter :: repetitions = 9, solves_per_sample = 20
    real(dp), parameter :: temperatures(3) = [0.05_dp, 1.0_dp, 10.0_dp]
    type(fermi_workspace_t) :: work
    type(fermi_state_t) :: state
    type(fermi_si_state_t) :: physical
    real(dp) :: start, finish, solve_time, complete_time, density, relativity
    integer :: sample, iteration, case_index, unit, si_unit
    logical :: ok
    character(256) :: path

    path = 'fermi-gas-timings.csv'
    if (command_argument_count() > 0) call get_command_argument(1, path)
    open (newunit=unit, file=trim(path), status='replace', action='write')
    write (unit, '(a)') 'theta,sample,solve_seconds,complete_seconds,factor,'// &
                       'candidate_nodes,strict_cells,mu_radius,pressure_radius'
    do case_index = 1, size(temperatures)
        call initialize_fermi_workspace(work, 16384, ok)
        if (.not. ok) error stop 'workspace failed'
        ! All numerical validation and adaptive setup are outside timed sections.
        call solve_fermi_gas(work, temperatures(case_index), state, ok)
        if (.not. ok .or. .not. state%certified) error stop 'setup certificate failed'
        do sample = 1, repetitions
            call cpu_time(start)
            do iteration = 1, solves_per_sample
                call solve_fermi_gas(work, temperatures(case_index), state, ok, &
                                     certify=.false.)
            end do
            call cpu_time(finish)
            solve_time = (finish - start)/real(solves_per_sample, dp)
            if (.not. ok .or. .not. state%valid) error stop 'candidate failed'
            call cpu_time(start)
            call solve_fermi_gas(work, temperatures(case_index), state, ok)
            call cpu_time(finish)
            complete_time = finish - start
            if (.not. ok .or. .not. state%certified) error stop 'certificate failed'
            write (unit, '(*(g0,:,","))') state%theta, sample, solve_time, &
                complete_time, complete_time/solve_time, state%candidate_nodes, &
                state%cells, (state%mu_bound%hi - state%mu_bound%lo)/2.0_dp, &
                (state%pressure_bound%hi - state%pressure_bound%lo)/2.0_dp
        end do
    end do
    close (unit)
    ! Illustrative envelope inputs, not a stellar structure or relativistic EOS.
    ! Fixed supplied CODATA-2022 central values; their measurement uncertainty
    ! is excluded from this numerical certificate and remains a caller ledger.
    call fermi_electron_density(1.0e6_dp, 0.5_dp, 1.66053906892e-27_dp, density)
    call initialize_fermi_workspace(work, 16384, ok)
    call solve_fermi_gas_si(work, density, 1.0e6_dp, 9.1093837139e-31_dp, &
        1.054571817e-34_dp, 1.380649e-23_dp, 1.602176634e-19_dp, &
        8.8541878188e-12_dp, physical, ok)
    if (.not. ok) error stop 'illustrative SI pressure failed'
    call fermi_relativity_indicator(physical%fermi_momentum, 9.1093837139e-31_dp, &
                                    299792458.0_dp, relativity)
    open (newunit=si_unit, file='fermi-gas-si.csv', status='replace', action='write')
    write (si_unit, '(a)') 'scope,n_m3,T_K,theta,pF_over_mc,P_Pa,P_lower_Pa,P_upper_Pa'
    write (si_unit, '(*(g0,:,","))') 'ideal-nonrel-electron-envelope-only', &
        density, 1.0e6_dp, physical%normalized%theta, relativity, &
        physical%pressure, physical%pressure_bound%lo, physical%pressure_bound%hi
    close (si_unit)
    print '(a)', 'PASS: raw matched-input timing rows and nonrelativistic SI pressure'
end program benchmark_fermi_gas
