program benchmark_error_quantum
    use, intrinsic :: iso_fortran_env, only: dp => real64, i64 => int64
    use fortnum_interval, only: interval_t
    use kin6d_quantum_fock, only: assemble_jellium_pair_sector, fock_hamiltonian_t
    use kin6d_quantum_evolution
    use kin6d_coulomb_encounter
    use kin6d_packet_models
    use kin6d_fermion_ensemble
    use kin6d_fermi_gas
    use kin6d_jellium_cell
    use kin6d_mean_field_cell
    use kin6d_screening_hierarchy
    use error_benchmark_support
    implicit none
    integer :: trial
    call header()
    do trial = 1, 5
        call pair(trial)
        call encounter(trial)
        call fermions(trial)
        call fermi(trial, 0._dp)
        call fermi(trial, 1._dp)
        call jellium(trial)
        call screening(trial)
    end do
contains
    subroutine pair(trial)
        integer, intent(in) :: trial
        integer, parameter :: n=125, steps=128
        integer :: modes(3,n), ix, iy, iz, k, shell
        real(dp) :: hr(n,n), hi(n,n), re(n), im(n), a, b, solve, radius
        type(interval_t) :: hrb(n,n), hib(n,n)
        type(quantum_evolution_t) :: state
        logical :: ok, bounded
        character(160) :: status
        k=0
        do shell=0,12
            do iz=-2,2
                do iy=-2,2
                    do ix=-2,2
                        if (ix*ix+iy*iy+iz*iz/=shell) cycle
                        k=k+1; modes(:,k)=[ix,iy,iz]
                    end do
                end do
            end do
        end do
        call assemble_jellium_pair_sector(modes, 2*acos(-1._dp), 1._dp, &
            1._dp, 2._dp, 100000000_i64, hr, hi, hrb, hib, ok, status)
        call require(ok)
        re=0; im=0; re(1)=1
        call initialize_quantum_evolution(state,hr,hi,re,im,steps,1._dp/128,1._dp,ok)
        call require(ok)
        call cpu_time(a)
        call solve_quantum_evolution(state,ok,certify=.false.)
        call cpu_time(b)
        call require(ok); solve=b-a
        call cpu_time(a)
        call solve_quantum_evolution(state,ok,error_mode=QE_STRICT)
        call cpu_time(b)
        call require(ok)
        call sample_quantum_evolution(state,steps,1._dp,re,im,radius,ok,bounded)
        call require(ok.and.bounded)
        call row('Q0-pair','S',trial,solve,0._dp,max(0._dp,b-a-solve), &
            radius,1._dp,'certified-discrete','n=125;steps=128;t=1;g=2', &
            'paired-difference;matrix-setup-excluded')
    end subroutine

    subroutine encounter(trial)
        integer, intent(in) :: trial
        type(encounter_config_t) :: cfg
        type(encounter_t) :: enc
        type(packet_spec_t) :: spec
        type(packet_sample_t), allocatable :: samples(:)
        type(packet_trace_t) :: trace
        complex(dp), allocatable :: state(:,:,:)
        character(:), allocatable :: message
        real(dp) :: a,b,c,d,solve,p,e,radius,setup
        integer(i64) :: nfev
        logical :: ok
        cfg%g=.3_dp; cfg%k=1; cfg%sigma=2; cfg%z0=6; cfg%t_final=1
        cfg%rmax=18; cfg%hmax=1; cfg%order=6; cfg%lmax=8; cfg%mmax=0
        call cpu_time(a)
        call encounter_init(enc,cfg,ok,message)
        call cpu_time(b)
        call require(ok)
        setup=b-a
        allocate(state(enc%mesh%n,0:enc%mmax,0:enc%lmax))
        call cpu_time(c)
        call encounter_state(enc,1._dp,state)
        call cpu_time(d)
        solve=setup+d-c-enc%ledger%estimator_seconds-enc%ledger%indicator_seconds
        call cpu_time(c)
        call encounter_state_error(enc,1._dp,1._dp,1._dp,e,radius)
        call cpu_time(a)
        call row('Q1-state','E',trial,solve, &
            max(0._dp,a-c)+enc%ledger%estimator_seconds+enc%ledger%indicator_seconds, &
            0._dp,radius,1._dp,'indicator;finite-space-formula', &
            'g=.3;sigma=2;z0=6;lmax=8;rmax=18;t=1', &
            'instrumented-embedded-residual-and-reconstruction')
        spec%model=MODEL_M3
        call cpu_time(a)
        call packet_samples(enc,spec,[0._dp,.25_dp,.5_dp,.75_dp,1._dp], &
            samples,ok,message,nfev)
        call cpu_time(b)
        call require(ok)
        call packet_trace(enc,spec,samples,trace)
        call cpu_time(c)
        call row('Q2-residual','E',trial,b-a,c-b,0._dp,trace%bound(5), &
            1._dp,'indicator;sampled-Duhamel','M3;Q1-mesh;samples=5;t=1', 'direct')
    end subroutine

    subroutine fermions(trial)
        integer, intent(in) :: trial
        type(fermion_spectrum_t) :: spectrum
        type(canonical_state_t) :: state
        real(dp) :: a,b,solve,radius
        logical :: ok
        call box_spectrum(60,spectrum,ok)
        call require(ok)
        call cpu_time(a)
        call solve_canonical(spectrum,7,.4_dp,state,ok,certify=.false.)
        call cpu_time(b)
        call require(ok); solve=b-a
        call cpu_time(a)
        call solve_canonical(spectrum,7,.4_dp,state,ok)
        call cpu_time(b)
        call require(ok.and.state%certified)
        radius=max(state%energy-state%energy_bound%lo, &
            state%energy_bound%hi-state%energy)
        call row('Q3-fermions','S',trial,solve,0._dp,max(0._dp,b-a-solve), &
            radius,7._dp,'certified-discrete;orbital-tail-included', &
            'cutoff=60;Nspin=7;b=.4','paired-difference')
    end subroutine

    subroutine fermi(trial,relativity)
        integer, intent(in) :: trial
        real(dp), intent(in) :: relativity
        type(fermi_workspace_t) :: work
        type(fermi_state_t) :: state
        real(dp) :: a,b,solve,radius
        logical :: ok
        call initialize_fermi_workspace(work,16384,ok)
        call solve_fermi_gas(work,1._dp,state,ok,relativity=relativity)
        call require(ok.and.state%certified)
        call cpu_time(a)
        call solve_fermi_gas(work,1._dp,state,ok,certify=.false.,relativity=relativity)
        call cpu_time(b)
        call require(ok); solve=b-a
        call cpu_time(a)
        call solve_fermi_gas(work,1._dp,state,ok,relativity=relativity)
        call cpu_time(b)
        call require(ok.and.state%certified)
        radius=max(state%pressure-state%pressure_bound%lo, &
            state%pressure_bound%hi-state%pressure)
        call row(merge('Q8-relativistic','Q8-Fermi-gas   ',relativity>0), &
            'D',trial,solve,0._dp,max(0._dp,b-a-solve),radius, &
            max(1._dp,abs(state%pressure)),'certified-continuous', &
            merge('theta=1;xF=1','theta=1;xF=0',relativity>0),'paired-difference;warm')
    end subroutine

    subroutine jellium(trial)
        integer, intent(in) :: trial
        type(jellium_cell_t) :: cell
        type(jellium_sector_t) :: sector, mfsector
        type(fock_hamiltonian_t) :: model, mfmodel
        type(mf_interaction_t) :: interaction
        type(mf_state_t) :: mf
        type(mf_estimate_t) :: estimate
        integer, allocatable :: modes(:,:)
        real(dp), allocatable :: vectors(:,:), psi(:)
        real(dp) :: values(2), residuals(2), times(0:128), response(0:128), band(0:128)
        real(dp) :: a,b,c,solve,check,norm
        integer :: iterations,dim,solutions,converged,nk
        logical :: ok
        call jellium_modes('spherical',1,modes)
        call init_jellium_cell(cell,2,modes,2._dp,ok)
        call jellium_sector(cell,1,reshape([0,0,0],[3,1]),sector,ok)
        call jellium_hamiltonian(cell,sector%basis,model,ok)
        allocate(vectors(sector%d,2))
        call jellium_lowest(model,sector%basis,2,values,vectors,residuals,ok,iterations)
        call require(ok)
        call jellium_kick_evolution(cell,model,sector,vectors(:,1),values(1),[1,0,0], &
            128,.01_dp,times,response,band,ok,dim,solve,check)
        call require(ok)
        call row('Q4-kick','S',trial,solve,0._dp,check,band(128),2._dp, &
            'certified-discrete','N=2;modes=7;steps=128;dt=.01','direct')
        call mf_interaction(cell,interaction)
        call cpu_time(a)
        call mf_ground(cell,interaction,.true.,.true.,1,mf,ok,solutions,converged)
        call cpu_time(b)
        call require(ok); solve=b-a
        call mf_slater_sector(cell,mf,mfsector,psi,norm,nk,ok)
        call jellium_hamiltonian(cell,mfsector%basis,mfmodel,ok)
        call cpu_time(c)
        call mf_residual(mfmodel,mfsector,psi,estimate,ok)
        call cpu_time(a)
        call require(ok)
        call row('Q5-Temple','E',trial,solve,a-c,0._dp, &
            abs(estimate%rho-estimate%temple),max(1._dp,abs(estimate%rho)), &
            'indicator;candidate-gap','N=2;modes=7;rs=2;RHF','direct')
    end subroutine

    subroutine screening(trial)
        integer, intent(in) :: trial
        type(screening_state_t) :: state
        type(screening_q_table_t) :: table
        real(dp) :: a,b,c,ratio(1),tail,indicator
        logical :: ok
        call screening_setup_normalized(1._dp,1._dp,state,ok)
        call require(ok)
        call cpu_time(a)
        call screening_q_table(state,32._dp,.5_dp,8._dp,table)
        call cpu_time(b)
        call screened_ratio_from_table(state,table,[1._dp],ratio,tail,indicator)
        call cpu_time(c)
        call row('Q7-screening','E',trial,b-a,c-b,0._dp,tail+indicator,1._dp, &
            'rigorous-tail;quadrature-indicator','theta=1;kappa2=1;qmax=32;r=1', &
            'embedded-tail-and-indicator;query-only')
    end subroutine
end program
