program landreman_prof
    ! Measurement-only driver (perf/landreman-w3-20261025): L/I cold start through
    ! the PF2 boundary/profile adapter, phase timers and solver counters.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_equilibrium_harness, only: equilibrium_resolution_t, equilibrium_case_t, &
        equilibrium_start_t, equilibrium_run_report_t, run_equilibrium_case, &
        write_equilibrium_report
    use kin6d_equilibrium_state, only: equilibrium_state_t
    use kin6d_equilibrium_3d_driver, only: equilibrium_3d_driver_t, &
        nested_integer_driver, run_boundary_profile_equilibrium, prof_x_seed
    use kin6d_map_spectrum_probe, only: map_spectrum_probe
    use kin6d_nested_oracles, only: nested_oracle_t
    use omp_lib, only: omp_get_wtime, omp_get_max_threads
    implicit none
    type(equilibrium_3d_driver_t) :: driver
    type(nested_oracle_t) :: oracle
    type(equilibrium_state_t) :: equilibrium
    type(equilibrium_run_report_t) :: report
    type(equilibrium_resolution_t) :: resolution
    character(64) :: arg, name
    integer :: mpol, ntor
    integer :: status
    real(dp) :: started, finished
    started = omp_get_wtime()
    oracle%eps = 0.5_dp; oracle%delta = 1.0_dp/64
    call get_command_argument(1, arg); read (arg, *) oracle%tau
    call get_command_argument(2, arg); read (arg, *) mpol
    call get_command_argument(3, arg); read (arg, *) ntor
    call get_command_argument(4, arg); read (arg, *) resolution%boundary_nodes
    if (command_argument_count() >= 5) then
        call get_command_argument(5, arg); read (arg, *) resolution%degree
    end if
    name = 'landreman_3d_nested'
    if (abs(oracle%tau) > tiny(1.0_dp)) name = 'issan_3d_nonsymmetric'
    call nested_integer_driver(oracle, mpol, ntor, driver, status)
    if (status /= 0) error stop 'integer-family boundary/profile data'
    call run_boundary_profile_equilibrium(equilibrium_case_t(trim(name), 'map', &
        'integer'), resolution, equilibrium_start_t(), driver, report, equilibrium)
    finished = omp_get_wtime()
    call write_equilibrium_report(6, report)
    print '(a,i0)', 'free_dofs=', count(.not. driver%map%fixed)
    print '(a,i0)', 'coefficient_dofs=', driver%map%ndof
    print '(a,i0)', 'reduced_dofs=', driver%map%nred
    print '(a,i0)', 'nodes=', driver%map%nodes
    print '(a,i0)', 'cells=', driver%map%cells
    print '(a,i0)', 'nmode=', driver%map%nmode
    print '(a,i0)', 'nphi=', driver%map%nphi
    print '(a,es24.16)', 'B_sample_rms_relative=', driver%magnetic_error
    print '(a,es24.16)', 'p_sample_rms_relative=', driver%pressure_error
    print '(a,es24.16)', 'surface_label_max_relative=', driver%label_error
    print '(a,es24.16)', 'workflow_wall_s=', finished - started
    print '(a,es24.16)', 'energy_gradient_wall_s=', driver%solve%energy_s
    print '(a,es24.16)', 'hessian_wall_s=', driver%solve%hessian_s
    print '(a,es24.16)', 'factor_solve_wall_s=', driver%solve%factor_s
    print '(a,es24.16)', 'analysis_s=', driver%solve%analysis_s
    print '(a,es24.16)', 'numeric_s=', driver%solve%numeric_s
    print '(a,i0)', 'newton_iterations=', driver%solve%iterations
    print '(a,i0)', 'factorizations=', driver%solve%factorizations
    print '(a,i0)', 'analyses=', driver%solve%analyses
#ifdef KIN6D_PROF_KRYLOV
    print '(a,3(1x,i0))', 'trials_ksolves_kits=', driver%solve%trials, &
        driver%solve%krylov_solves, driver%solve%krylov_iterations
#endif
    print '(a,es24.16)', 'final_residual=', report%residual
    print '(a,i0)', 'cores=', omp_get_max_threads()
    block
        character(8) :: flag
        call get_environment_variable('KIN6D_PROF_SPECTRUM', flag)
        if (flag == '1') call map_spectrum_probe(driver%map, prof_x_seed)
    end block
    if (report%status /= 0) error stop 'native D8 3D equilibrium'
    if (report%residual > driver%tolerance) error stop 'Newton tolerance not attained'
end program
