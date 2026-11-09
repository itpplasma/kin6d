program test_equilibrium_3d
    ! test_equilibrium_3d TAU MPOL NTOR NODES DEGREE GATE_B GATE_P GATE_LABEL
    !                     [COARSE_MPOL COARSE_NTOR COARSE_NODES COARSE_DEGREE]
    ! Gates are about twice the observed sampled errors (docs/ideal-seed-3d.md).
    ! Without a coarse level, a reversed-iota control must be clearly worse.
    ! With a coarse level, every sampled error must decrease under refinement.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use milca_equilibrium_harness, only: equilibrium_resolution_t, &
        equilibrium_run_report_t, run_equilibrium_case, write_equilibrium_report
    use milca_equilibrium_state, only: equilibrium_state_t
    use kin6d_equilibrium_3d_driver, only: equilibrium_3d_driver_t, &
        nested_integer_driver, run_boundary_profile_equilibrium
    use kin6d_nested_oracles, only: nested_oracle_t
    use omp_lib, only: omp_get_wtime, omp_get_max_threads
    implicit none
    type(nested_oracle_t) :: oracle
    type(equilibrium_resolution_t) :: resolution, coarse
    character(64) :: arg
    integer :: mpol, ntor, cmpol, cntor
    real(dp) :: gates(3), fine(3), base(3), reversed(3)
    oracle%eps=0.5_dp;oracle%delta=1.0_dp/64
    if(command_argument_count()/=8 .and. command_argument_count()/=12) &
        error stop 'usage: TAU MPOL NTOR NODES DEGREE GATE_B GATE_P GATE_LABEL [coarse]'
    call get_command_argument(1,arg);read(arg,*) oracle%tau
    call get_command_argument(2,arg);read(arg,*) mpol
    call get_command_argument(3,arg);read(arg,*) ntor
    call get_command_argument(4,arg);read(arg,*) resolution%boundary_nodes
    call get_command_argument(5,arg);read(arg,*) resolution%degree
    call get_command_argument(6,arg);read(arg,*) gates(1)
    call get_command_argument(7,arg);read(arg,*) gates(2)
    call get_command_argument(8,arg);read(arg,*) gates(3)
    call solve(mpol,ntor,resolution,1.0_dp,.true.,fine)
    if(fine(1)>gates(1)) error stop 'sampled B error above recorded gate'
    if(fine(2)>gates(2)) error stop 'sampled pressure error above recorded gate'
    if(fine(3)>gates(3)) error stop 'sampled surface-label error above recorded gate'
    if(command_argument_count()==12) then
        call get_command_argument(9,arg);read(arg,*) cmpol
        call get_command_argument(10,arg);read(arg,*) cntor
        call get_command_argument(11,arg);read(arg,*) coarse%boundary_nodes
        call get_command_argument(12,arg);read(arg,*) coarse%degree
        call solve(cmpol,cntor,coarse,1.0_dp,.false.,base)
        print '(a,3es12.4)', 'refinement_ratio_B_p_label=',base/fine
        if(any(fine>=base)) error stop 'sampled errors did not decrease under refinement'
    else
        ! The prescribed iota sign is physical: the reversed transform gives
        ! the same pressure and labels but must not reproduce the oracle B.
        call solve(mpol,ntor,resolution,-1.0_dp,.false.,reversed)
        print '(a,es12.4)', 'reversed_iota_B_sample_rms_relative=',reversed(1)
        if(reversed(1)<5*fine(1)) error stop 'B is insensitive to the iota sign'
    end if
contains
    subroutine solve(mpol,ntor,resolution,iota_sign,gate,errors)
        integer, intent(in) :: mpol, ntor
        type(equilibrium_resolution_t), intent(in) :: resolution
        real(dp), intent(in) :: iota_sign
        logical, intent(in) :: gate
        real(dp), intent(out) :: errors(3)
        type(equilibrium_3d_driver_t) :: driver
        type(equilibrium_state_t) :: equilibrium
        type(equilibrium_run_report_t) :: report
        character(64) :: name
        integer :: status
        real(dp) :: started, finished
        name='landreman_3d_nested'
        if(abs(oracle%tau)>tiny(1.0_dp)) name='issan_3d_nonsymmetric'
        started=omp_get_wtime()
        call nested_integer_driver(oracle,mpol,ntor,driver,status)
        if(status/=0) error stop 'integer-family boundary/profile data'
        if(driver%profile%iota(1)<=0) error stop 'integer-family iota must be positive'
        driver%profile%iota=iota_sign*driver%profile%iota
        call run_equilibrium_case(trim(name),'map',resolution,'legacy', &
            run_boundary_profile_equilibrium,driver,report,equilibrium)
        finished=omp_get_wtime()
        errors=[driver%magnetic_error,driver%pressure_error,driver%label_error]
        print '(a,2i3,2i4,f6.1)', 'case mpol ntor nodes degree iota_sign=', &
            mpol,ntor,resolution%boundary_nodes,resolution%degree,iota_sign
        call write_equilibrium_report(6,report)
        print '(a,es24.16)', 'B_sample_rms_relative=',errors(1)
        print '(a,es24.16)', 'p_sample_rms_relative=',errors(2)
        print '(a,es24.16)', 'surface_label_max_relative=',errors(3)
        print '(a,es24.16)', 'workflow_wall_s=',finished-started
        print '(a,es24.16)', 'energy_gradient_wall_s=',driver%solve%energy_s
        print '(a,es24.16)', 'hessian_wall_s=',driver%solve%hessian_s
        print '(a,es24.16)', 'factor_solve_wall_s=',driver%solve%factor_s
        print '(a,i0)', 'cores=',omp_get_max_threads()
        if(.not. gate) return
        print '(a)', '[V] 1728 equal-weight candidate-coordinate samples; no volume norm'
        print '(a)', '[V] force and a posteriori estimator effectivity unavailable'
        if(report%status/=0) error stop 'native D8 3D equilibrium'
        if(report%residual>driver%tolerance) error stop 'Newton tolerance not attained'
        if(report%oracle_error>=report%seed_error) error stop 'seed did not improve'
    end subroutine
end program
