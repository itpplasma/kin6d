program benchmark_error_modes
    use, intrinsic :: iso_fortran_env, only: dp => real64, i64 => int64
    use fortnum_interval, only: interval_t
    use fortnum_rng, only: rng_t, rng_seed
    use fortnum_status, only: fortnum_status_t
    use kin6d_quantum_fock, only: assemble_jellium_pair_sector
    use kin6d_quantum_evolution
    use error_benchmark_support
    implicit none
    integer :: trial, extent
    call header()
    do trial=1,5
        do extent=2,3
            call measure_pair(extent,trial)
        end do
    end do
contains
    subroutine measure_pair(extent,trial)
        integer, intent(in) :: extent,trial
        integer :: n,steps,ix,iy,iz,k,shell,mode,order(3),which
        integer, allocatable :: modes(:,:)
        real(dp), allocatable :: hr(:,:),hi(:,:),re(:),im(:)
        type(interval_t), allocatable :: hrb(:,:),hib(:,:)
        type(quantum_evolution_t) :: state
        type(rng_t) :: draw
        type(fortnum_status_t) :: rng_status
        real(dp) :: a,b,c,d,solve,estimate,radius,epsilon_b,allowance,indicator
        logical :: ok
        character(160) :: status
        character(80) :: inputs
        character(1) :: tier
        n=(2*extent+1)**3; steps=128
        if (extent==3) steps=64
        allocate(modes(3,n),hr(n,n),hi(n,n),re(n),im(n),hrb(n,n),hib(n,n))
        k=0
        do shell=0,3*extent**2
            do iz=-extent,extent
                do iy=-extent,extent
                    do ix=-extent,extent
                        if (ix*ix+iy*iy+iz*iz/=shell) cycle
                        k=k+1; modes(:,k)=[ix,iy,iz]
                    end do
                end do
            end do
        end do
        call assemble_jellium_pair_sector(modes,2*acos(-1._dp),1._dp,1._dp, &
            2._dp,100000000_i64,hr,hi,hrb,hib,ok,status)
        call require(ok)
        re=0; im=0; re(1)=1
        call initialize_quantum_evolution(state,hr,hi,re,im,steps, &
            1._dp/steps,1._dp,ok)
        call require(ok)
        call solve_quantum_evolution(state,ok,certify=.false.) ! warmup
        call require(ok)
        call cpu_time(a)
        call solve_quantum_evolution(state,ok,certify=.false.)
        call cpu_time(b)
        call require(ok); solve=b-a
        write(inputs,'(a,i0,a,i0,a)') 'n=',n,';steps=',steps,';t=1;g=2'
        call cpu_time(a)
        do k=1,1000
            call quantum_evolution_indicator(state,indicator,status)
        end do
        call cpu_time(b)
        call row('Q0-pair','E',trial,solve,(b-a)/1000,0._dp,indicator,1._dp, &
            'indicator',trim(inputs),'direct;query-only;embedded-in-solve')
        ! Alternate audit/workhorse order; the candidate bytes stay fixed.
        order=[QE_WORKHORSE,QE_STRICT,QE_APPROXIMATE]
        if (mod(trial,2)==0) order=[QE_STRICT,QE_WORKHORSE,QE_APPROXIMATE]
        do which=1,3
            mode=order(which)
            tier='D'
            if (mode==QE_STRICT) tier='S'
            if (mode==QE_APPROXIMATE) tier='A'
            call cpu_time(a)
            call certify_quantum_evolution(state,ok,mode)
            call cpu_time(b)
            call require(ok)
            call cpu_time(c)
            call quantum_evolution_error(state,steps,radius,status,epsilon_b,ok, &
                evaluation_error=allowance)
            call cpu_time(d)
            call require(ok)
            if (mode==QE_APPROXIMATE) then
                call row('Q0-pair',tier,trial,solve,d-c+b-a,0._dp,radius,1._dp, &
                    trim(status),trim(inputs),'direct;approximate-evaluation')
            else
                call row('Q0-pair',tier,trial,solve,d-c,b-a,radius,1._dp,trim(status), &
                    trim(inputs),'direct;deferred-all-cell-audit')
            end if
        end do
        call rng_seed(draw,27182_i64,rng_status)
        call cpu_time(a)
        call quantum_evolution_sampled_residual(state,steps,.5_dp,24,.05_dp,draw, &
            estimate,radius,ok)
        call cpu_time(b)
        call require(ok)
        call row('Q0-residual-norm','P',trial,solve,b-a,0._dp,radius,1._dp, &
            'probabilistically-certified residual norm; confidence=.95; ideal-iid', &
            trim(inputs)//';samples=24;seed=27182', &
            'direct;one-time-norm;state-error-unavailable')
    end subroutine
end program
