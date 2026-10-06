program benchmark_error_microscopic
    use, intrinsic :: iso_fortran_env, only: dp => real64, i64 => int64
    use fortnum_status, only: fortnum_status_t, FORTNUM_OK
    use fortnum_rng, only: rng_t, rng_seed, rng_uniform
    use kin6d_species, only: electron, atomic_units, species_t
    use kin6d_pme, only: pme_parameters_t, pme_select, pme_estimate
    use kin6d_nbody, only: nbody_system_t, nbody_periodic_mesh, nbody_force
    use kin6d_ocp_monte_carlo, only: ocp_mc_t, mc_setup, mc_sweep, mc_energy
    use kin6d_series_statistics, only: series_batch_means
    use kin6d_classical_cell, only: classical_cell_t, classical_mc_t, &
        init_classical_cell, classical_cell_mc, model_coulomb
    use error_benchmark_support
    implicit none
    integer, parameter :: n=160
    type(species_t) :: particles(n)
    type(nbody_system_t) :: system
    type(pme_parameters_t) :: parameters
    type(rng_t) :: random
    type(fortnum_status_t) :: status
    type(ocp_mc_t) :: chain
    type(classical_cell_t) :: cell
    type(classical_mc_t) :: mc
    real(dp) :: y(6*n), force(3*n), potential, a,b,c,solve,estimate,radius
    real(dp) :: energy(1024),mean,error,edges(7), scale
    integer :: trial,i,which
    logical :: ok
    call header()
    particles=electron()
    call rng_seed(random,17631_i64,status)
    do i=1,3*n
        call rng_uniform(random,y(i))
    end do
    y(:3*n)=5*y(:3*n); y(3*n+1:)=0
    do trial=1,5
        do which=1,2
            call pme_select(5._dp,n,1._dp,merge(1.e-4_dp,1.e-5_dp,which==1), &
                parameters,ok)
            call require(ok)
            call nbody_periodic_mesh(system,atomic_units(),particles,parameters,ok)
            call require(ok)
            call nbody_force(system,y,force,potential,status)
            call require(status%code==FORTNUM_OK)
            call cpu_time(a)
            do i=1,10
                call nbody_force(system,y,force,potential,status)
            end do
            call cpu_time(b)
            call require(status%code==FORTNUM_OK); solve=(b-a)/10
            call cpu_time(a)
            call pme_estimate(parameters,n,1._dp)
            call cpu_time(b)
            radius=2*sqrt(parameters%real_estimate**2+parameters%reciprocal_estimate**2)
            call row(merge('C3-PME','C4-PME',which==1),'E',trial,solve,b-a,0._dp, &
                radius,parameters%force_scale,'indicator;uncorrelated-position-law', &
                merge('N=160;L=5;tol=1e-4','N=160;L=5;tol=1e-5',which==1), &
                'direct;selection-and-setup-excluded')
        end do
        call mc_setup(chain,16,1._dp,1.e-8_dp,27_i64,ok)
        call require(ok)
        call cpu_time(a)
        do i=1,256
            call mc_sweep(chain,tune=.true.)
        end do
        do i=1,size(energy)
            call mc_sweep(chain,tune=.false.)
            energy(i)=mc_energy(chain)/16
        end do
        call cpu_time(b)
        call series_batch_means(energy,16,mean,error)
        call cpu_time(c)
        call row('H1a-statistical','statistical',trial,b-a,c-b,0._dp, &
            1.96_dp*error,1._dp,'statistical;asymptotic-95-percent', &
            'N=16;Gamma=1;sweeps=1024;equil=256;batches=16;seed=27','direct')
        call init_classical_cell(cell,3,2._dp,1._dp,model_coulomb,ok)
        call require(ok)
        edges=[(real(i,dp)*cell%length/12,i=0,6)]
        call classical_cell_mc(cell,edges,1024,256,16,27_i64,mc,ok, &
            sample_seconds=solve,summary_seconds=estimate)
        call require(ok)
        scale=max(1._dp,abs(mc%energy))
        call row('H1b-statistical','statistical',trial,solve,estimate,0._dp, &
            1.96_dp*mc%energy_err,scale,'statistical;asymptotic-95-percent', &
            'N=3;rs=2;Theta=1;sweeps=1024;equil=256;batches=16;seed=27','direct')
    end do
end program
