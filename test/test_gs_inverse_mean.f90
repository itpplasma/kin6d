program test_gs_inverse_mean
    use, intrinsic :: iso_fortran_env, only: dp=>real64
    use kin6d_gs_inverse_leaves, only: gs_profile_mean_last, &
        gs_profile_mean_last_jvp, gs_profile_mean_last_vjp
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use kin6d_gs_inverse, only: gs_inverse_result_t, gs_inverse_f, gs_inverse_profile_status
    implicit none
    type(gs_inverse_result_t)::profile
    real(dp)::f
    integer::status
    real(dp)::last, scaled, delta, gmean, gpartial, gweight, pairing
    real(dp)::weights(3), coefficients(3), mean, partial, amplitude, tangent(3)
    integer::i
    ! Independent trapezoidal integral on the nonuniform knots [0,.25,1].
    weights=[.125_dp,.5_dp,.375_dp]
    coefficients(:2)=[-2.0_dp,-.5_dp]
    mean=-1;partial=dot_product(weights(:2),coefficients(:2))
    call gs_profile_mean_last(mean,partial,weights(3),last)
    coefficients(3)=last
    if(abs(dot_product(weights,coefficients)-mean)>2.0e-15_dp) &
        error stop 'complete linear-profile mean constraint'
    if(abs(last+4.0_dp/3)>2.0e-15_dp)error stop 'independent nonuniform integral oracle'
    tangent=[.1_dp,-.05_dp,.02_dp]
    call gs_profile_mean_last_jvp(mean,partial,weights(3),tangent(1), &
        tangent(2),tangent(3),delta)
    ! Differentiated defining constraint: w*dlast + last*dw + dpartial = dmean.
    if(abs(weights(3)*delta+last*tangent(3)+tangent(2)-tangent(1))>2.0e-15_dp) &
        error stop 'mean constraint equation-level JVP'
    call gs_profile_mean_last_vjp(mean,partial,weights(3),.7_dp,gmean,gpartial,gweight)
    pairing=dot_product([gmean,gpartial,gweight],tangent)
    if(abs(pairing-.7_dp*delta)>2.0e-15_dp)error stop 'mean Euclidean JVP/VJP pairing'
    do i=1,3
        amplitude=10.0_dp**(i-2)
        call gs_profile_mean_last(amplitude*mean,amplitude*partial,weights(3),scaled)
        if(abs(scaled/amplitude-last)>2.0e-15_dp) &
            error stop 'positive zero-beta profile homogeneity'
    end do
    profile%flux=[-1.0_dp,0.0_dp];profile%ffprime=[0.0_dp,0.0_dp]
    profile%f_squared=[1.0_dp,1.0_dp];profile%q_target=-1
    call gs_inverse_profile_status(profile,status)
    if(status/=0)error stop 'positive constant F profile rejected'
    call gs_inverse_f(profile,-.5_dp,f,status)
    if(status/=0.or.f/=-1)error stop 'signed positive F branch'
    profile%ffprime=[-4.0_dp,4.0_dp]
    ! Both endpoints remain +1, while exact integration gives F^2(-.5)=-1.
    call gs_inverse_profile_status(profile,status)
    if(status==0)error stop 'positive knots hid a negative interior F squared'
    call gs_inverse_f(profile,-.5_dp,f,status)
    if(status==0)error stop 'negative F squared query accepted'
    profile%ffprime=0;profile%q_target=ieee_value(1.0_dp,ieee_quiet_nan)
    call gs_inverse_f(profile,-.5_dp,f,status)
    if(status==0)error stop 'NaN signed q branch accepted'
    print *, 'PASS: exact integral mean, equation derivative, adjoint, positive homogeneity'
end program test_gs_inverse_mean
