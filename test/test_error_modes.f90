program test_error_modes
    use, intrinsic :: iso_fortran_env, only: dp => real64, qp => real128, i64 => int64
    use kin6d_quantum_evolution
    use kin6d_slab_oscillation, only: slab_oscillation_state_t, &
        initialize_slab_oscillation, slab_oscillation_indicator
    use kin6d_fp_laboratory, only: fp_state_t, initialize_fp, solve_fp, fp_indicator
    use fortnum_rng, only: rng_t, rng_seed
    use fortnum_status, only: fortnum_status_t
    implicit none
    type(quantum_evolution_t) :: state
    type(slab_oscillation_state_t) :: slab
    type(fp_state_t) :: fp
    type(rng_t) :: draw
    type(fortnum_status_t) :: rng_status
    real(dp) :: hr(2,2),hi(2,2),initial(2),zero(2),re(2),im(2),h,radius
    real(dp) :: epsilon_b, audit(3), eta(3), allowance, indicator, residual
    real(dp) :: sampled, upper, actual, worst_under, loosen, worst_audit
    complex(qp) :: exact(2)
    integer :: mode,level,point,failures,seed
    logical :: ok,bounded
    character(128) :: status
    hr=0; hr(1,1)=1; hr(2,2)=-.75_dp; hi=0
    initial=[.6_dp,.8_dp]; zero=0
    worst_under=0; loosen=0; worst_audit=0
    do level=1,3
        h=.05_dp*2**(level-1)
        do mode=QE_WORKHORSE,QE_APPROXIMATE
            call initialize_quantum_evolution(state,hr,hi,initial,zero,1,h,1._dp,ok)
            call require(ok,'initialize')
            call solve_quantum_evolution(state,ok,certify=.false.)
            call require(ok,'candidate')
            call certify_quantum_evolution(state,ok,mode)
            call require(ok,'deferred error evaluation')
            do point=1,7
                call sample_quantum_evolution(state,1,real(point,dp)/8,re,im, &
                    radius,ok,bounded)
                call require(ok,'candidate sample')
                call require(bounded .eqv. (mode/=QE_APPROXIMATE),'status separation')
                call quantum_evolution_error(state,1,radius,status,epsilon_b,ok, &
                    fraction=real(point,dp)/8,evaluation_error=allowance)
                call require(ok,'error query')
                ! Exact diagonal-unitary evolution is independent of RK and
                ! interval residual reconstruction, evaluated in binary128.
                exact(1)=real(initial(1),qp)*exp(cmplx(0._qp, &
                    -real(point,qp)*real(h,qp)/8,qp))
                exact(2)=real(initial(2),qp)*exp(cmplx(0._qp, &
                    .75_qp*real(point,qp)*real(h,qp)/8,qp))
                actual=real(sqrt(sum(abs( &
                    cmplx(real(re,qp),real(im,qp),qp)-exact)**2)),dp)
                if (mode==QE_APPROXIMATE) then
                    call require(index(status,'approximate-bound')==1,'A label')
                    call require(actual<=(1+epsilon_b)*radius,'A true-error oracle')
                    call require(allowance<=epsilon_b*radius,'A omitted-term gate')
                    worst_under=max(worst_under,actual/radius)
                else
                    call require(actual<=radius,'D/S containment oracle')
                end if
                if (point==4) then
                    eta(mode)=radius
                    if (mode==QE_WORKHORSE) loosen=max(loosen,radius/actual)
                end if
            end do
            call quantum_evolution_indicator(state,indicator,status)
            call require(status=='indicator'.and.indicator>=0,'embedded indicator')
        end do
        worst_audit=max(worst_audit,eta(QE_STRICT)/eta(QE_APPROXIMATE)-1)
        call require(eta(QE_STRICT)<=1.05_dp*eta(QE_APPROXIMATE),'strict audit of A')
    end do
    print '(a,es12.4)', 'A worst true-error/eta=',worst_under
    print '(a,es12.4)', 'A worst strict-radius/eta minus 1=',worst_audit
    print '(a,es12.4)', 'D max radius/true-error at midpoint=',loosen
    ! A relative allowance below arithmetic resolution must fail explicitly.
    call initialize_quantum_evolution(state,hr,hi,initial,zero,1,1.e-12_dp,1._dp,ok)
    call solve_quantum_evolution(state,ok,certify=.false.)
    call certify_quantum_evolution(state,ok,QE_APPROXIMATE,epsilon_b=1.e-20_dp)
    call require(.not.ok,'A gate must reject unresolved relative estimator error')
    call quantum_evolution_error(state,1,radius,status,epsilon_b,ok)
    call require(.not.ok.and.status=='certification-failed','A failure status')
    ! P estimates a residual NORM at one time. Independent exact residual:
    ! diagonal H acts on an exact Hermite polynomial of the candidate endpoints.
    call initialize_quantum_evolution(state,hr,hi,initial,zero,1,.2_dp,1._dp,ok)
    call solve_quantum_evolution(state,ok,certify=.false.)
    failures=0
    do seed=1,200
        call rng_seed(draw,int(seed,i64),rng_status)
        call quantum_evolution_sampled_residual(state,1,.5_dp,32,.05_dp,draw, &
            sampled,upper,ok)
        call require(ok,'P estimator')
        ! Analytic midpoint residual of cubic Hermite based on RK endpoint,
        ! in quad precision; no production interval/action routine is used.
        call sample_quantum_evolution(state,1,1._dp,re,im,radius,ok,bounded)
        call midpoint_residual(cmplx(real(re,qp),real(im,qp),qp),residual)
        if (residual>upper) failures=failures+1
    end do
    call require(failures<=20,'P empirical coverage with binomial margin')
    print '(a,i0,a)', 'P Q0 residual coverage failures=',failures,'/200'
    call initialize_slab_oscillation(slab,32,1._dp,1._dp,0._dp,ok)
    call slab_oscillation_indicator(slab,indicator,status)
    call require(status=='indicator'.and.indicator>0,'slab embedded indicator')
    call initialize_fp(fp,32)
    call solve_fp(fp,residual)
    call fp_indicator(fp,indicator,status)
    call require(status=='indicator'.and.indicator==residual,'FP indicator')
contains
    subroutine require(condition,message)
        logical, intent(in) :: condition
        character(*), intent(in) :: message
        if (.not.condition) then
            print *, 'FAIL: ',message
            error stop 1
        end if
    end subroutine
    subroutine midpoint_residual(endpoint,result)
        complex(qp), intent(in) :: endpoint(2)
        real(dp), intent(out) :: result
        complex(qp) :: y0(2),f0(2),f1(2),c1(2),c2(2),c3(2),y(2),rate(2),r(2)
        real(qp) :: energies(2),dt
        energies=[1._qp,-.75_qp]; dt=real(.2_dp,qp)
        y0=cmplx(real(initial,qp),0._qp,qp)
        f0=cmplx(0._qp,-energies,qp)*y0
        f1=cmplx(0._qp,-energies,qp)*endpoint
        c1=dt*f0; c2=3*(endpoint-y0)-dt*(2*f0+f1)
        c3=2*(y0-endpoint)+dt*(f0+f1)
        y=y0+.5_qp*(c1+.5_qp*(c2+.5_qp*c3))
        rate=(c1+c2+.75_qp*c3)/dt
        r=rate-cmplx(0._qp,-energies,qp)*y
        result=real(sqrt(sum(abs(r)**2)),dp)
    end subroutine
end program
