program benchmark_ocp_app
    ! Cost of the scalable periodic OCP action against N (C3/C4): PME
    ! parameter selection, setup, one force evaluation, one encounter-
    ! resolved step and its tangent, and the attributed scratch memory, at
    ! fixed Lambda and requested force tolerance. For N <= direct_limit the
    ! direct Ewald sum (oracle) is timed on the same prepared configuration
    ! and the measured rms force error is reported against the estimate.
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use fortnum_status, only: fortnum_status_t, FORTNUM_OK
    use kin6d_species, only: electron, species_t
    use kin6d_ewald, only: ewald_parameters_t, ewald_select
    use kin6d_nbody, only: nbody_system_t, nbody_periodic, nbody_force, nbody_step, &
                           nbody_step_jvp
    use kin6d_ocp, only: ocp_model_t, ocp_system, ocp_prepare, ocp_units
    use kin6d_result_io, only: run_t, result_file_t
    implicit none
    integer, parameter :: cases = 9
    real(dp) :: debye_volume(cases) = [100.0_dp, 100.0_dp, 100.0_dp, 100.0_dp, 300.0_dp, &
                                       100.0_dp, 1000.0_dp, 100.0_dp, 100.0_dp]
    real(dp) :: box(cases) = [2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, 4.0_dp, 6.0_dp, 4.0_dp, &
                              8.0_dp, 10.0_dp]
    real(dp) :: tolerance = 1.0e-4_dp
    integer :: steps = 10, direct_limit = 6400
    type(run_t) :: run
    type(result_file_t) :: result
    type(ocp_model_t) :: model
    type(nbody_system_t) :: system, direct
    type(fortnum_status_t) :: status
    type(ewald_parameters_t) :: ep
    type(species_t), allocatable :: particles(:)
    real(dp), allocatable :: state(:), tangent(:), force(:), reference(:)
    real(dp) :: table(cases, 13), t0, t1, u, a, memory
    integer :: c, k, unit, ios, rejected
    logical :: ok
    character(256) :: msg
    namelist /benchmark_ocp/ debye_volume, box, tolerance, steps, direct_limit
    call run%start('benchmark_ocp')
    call run%open_config(unit)
    if (run%has_config()) read (unit, nml=benchmark_ocp, iostat=ios, iomsg=msg)
    call run%close_config(unit, 'benchmark_ocp', ios, msg)
    call run%open_resolved(unit)
    write (unit, nml=benchmark_ocp)
    call run%close_resolved(unit)
    table = -1
    open (newunit=unit, file='output/benchmark_ocp.csv', status='replace', action='write')
    write (unit, '(a)') 'lambda,box,n,cutoff_over_a,mesh,order,select_setup_s,reserved,'// &
        'force_s,step_s,tangent_step_s,memory_mb,estimate_over_f0,direct_s,'// &
        'measured_over_f0'
    do c = 1, cases
        model = ocp_model_t()
        model%debye_volume = debye_volume(c); model%box = box(c)
        model%tolerance = tolerance
        call cpu_time(t0)
        call ocp_system(model, system, ok)
        call cpu_time(t1)
        if (.not. ok) error stop 'OCP system setup failed'
        allocate (state(6*model%particles), tangent(6*model%particles))
        allocate (force(3*model%particles), reference(3*model%particles))
        call ocp_prepare(model, 17_int64, state, rejected, status)
        table(c, 1:3) = [model%debye_volume, model%box, real(model%electrons, dp)]
        a = (3/(4*3.14159265358979_dp*model%debye_volume))**(1.0_dp/3)
        table(c, 4:6) = [model%pme%cutoff/a, real(model%pme%mesh, dp), &
                         real(model%pme%order, dp)]
        table(c, 7) = t1 - t0
        call cpu_time(t0)
        call nbody_force(system, state(:3*model%particles), force, u, status)
        call cpu_time(t1)
        table(c, 9) = t1 - t0
        call cpu_time(t0)
        do k = 1, steps
            call nbody_step(system, state, model%dt, status)
            if (status%code /= FORTNUM_OK) error stop 'step failed'
        end do
        call cpu_time(t1)
        table(c, 10) = (t1 - t0)/steps
        tangent = 1.0e-3_dp
        call cpu_time(t0)
        do k = 1, max(1, steps/2)
            call nbody_step_jvp(system, state, tangent, model%dt, status)
        end do
        call cpu_time(t1)
        table(c, 11) = (t1 - t0)/max(1, steps/2)
        ! Mesh (two complex K^3 grids, influence) and per-particle scratch.
        memory = (3*16 + 0)*real(model%pme%mesh, dp)**3/2 + &
                 model%particles*(8*(3*3*model%pme%order + 6*6 + 3*12) + 4*8)
        table(c, 12) = memory/1.0e6_dp
        table(c, 13) = sqrt(model%pme%real_estimate**2 + &
                            model%pme%reciprocal_estimate**2)/model%pme%force_scale
        if (model%electrons <= direct_limit) then
            allocate (particles(model%particles)); particles = electron()
            call ewald_select(model%box, 1.0e-10_dp, ep, ok, particles=model%particles)
            call nbody_periodic(direct, ocp_units(model), particles, ep, ok)
            call cpu_time(t0)
            call nbody_force(direct, state(:3*model%particles), reference, u, status)
            call cpu_time(t1)
            call nbody_force(system, state(:3*model%particles), force, u, status)
            write (unit, '(*(g0.6,:,","))') table(c, :), t1 - t0, &
                sqrt(sum((force - reference)**2)/model%electrons)/model%pme%force_scale
            deallocate (particles)
        else
            write (unit, '(*(g0.6,:,","))') table(c, :), -1.0_dp, -1.0_dp
        end if
        flush (unit)
        print '(a,i7,a,f8.4,a,f8.4,a,i4)', 'N=', model%electrons, ' force s=', &
            table(c, 9), ' step s=', table(c, 10), ' K=', model%pme%mesh
        deallocate (state, tangent, force, reference)
    end do
    close (unit)
    call result%create(run)
    call result%write('benchmark', table, '1', &
        'lambda, box, N, rc/a, K, p, select/setup/force/step/tangent s, MB, est', &
        dims='case,column')
    call result%close()
end program benchmark_ocp_app
