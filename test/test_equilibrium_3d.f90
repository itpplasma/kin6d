program test_equilibrium_3d
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_equilibrium_harness, only: equilibrium_resolution_t, &
        equilibrium_run_report_t, run_equilibrium_case, write_equilibrium_report
    use kin6d_equilibrium_state, only: equilibrium_state_t
    use kin6d_equilibrium_3d_driver, only: equilibrium_3d_driver_t, &
        nested_integer_driver, run_boundary_profile_equilibrium
    use kin6d_nested_oracles, only: nested_oracle_t
    use omp_lib, only: omp_get_wtime, omp_get_max_threads
    implicit none
    type(equilibrium_3d_driver_t) :: driver
    type(nested_oracle_t) :: oracle
    type(equilibrium_state_t) :: equilibrium
    type(equilibrium_run_report_t) :: report
    type(equilibrium_resolution_t) :: resolution
    character(64) :: arg, name
    integer :: status, mpol, ntor
    real(dp) :: started, finished
    oracle%eps=0.5_dp;oracle%delta=1.0_dp/64
    call get_command_argument(1,arg);read(arg,*) oracle%tau
    call get_command_argument(2,arg);read(arg,*) mpol
    call get_command_argument(3,arg);read(arg,*) ntor
    call get_command_argument(4,arg);read(arg,*) resolution%boundary_nodes
    if(command_argument_count()>=5) then
        call get_command_argument(5,arg);read(arg,*) resolution%degree
    end if
    name='landreman_3d_nested'
    if(abs(oracle%tau)>tiny(1.0_dp)) name='issan_3d_nonsymmetric'
    started=omp_get_wtime()
    call nested_integer_driver(oracle,mpol,ntor,driver,status)
    if(status/=0) error stop 'integer-family boundary/profile data'
    call run_equilibrium_case(trim(name),'map',resolution,'legacy', &
        run_boundary_profile_equilibrium,driver,report,equilibrium)
    finished=omp_get_wtime()
    call write_equilibrium_report(6,report)
    print '(a,es24.16)', 'B_sample_rms_relative=',driver%magnetic_error
    print '(a,es24.16)', 'p_sample_rms_relative=',driver%pressure_error
    print '(a,es24.16)', 'surface_label_max_relative=',driver%label_error
    print '(a,es24.16)', 'workflow_wall_s=',finished-started
    print '(a,es24.16)', 'energy_gradient_wall_s=',driver%solve%energy_s
    print '(a,es24.16)', 'hessian_wall_s=',driver%solve%hessian_s
    print '(a,es24.16)', 'factor_solve_wall_s=',driver%solve%factor_s
    print '(a,i0)', 'cores=',omp_get_max_threads()
    print '(a)', '[V] 1728 equal-weight candidate-coordinate samples; no volume norm'
    print '(a)', '[V] force and a posteriori estimator effectivity unavailable'
    if(report%status/=0) error stop 'native D8 3D equilibrium'
    if(report%residual>driver%tolerance) error stop 'Newton tolerance not attained'
    ! Coarse smoke: reject a wrong field or unchanged cold start; these are
    ! not the physical-accuracy gates of the benchmark ladder.
    if(report%oracle_error>0.2_dp) error stop 'coarse B oracle sanity'
    if(report%oracle_error>=report%seed_error) error stop 'seed did not improve'
end program
