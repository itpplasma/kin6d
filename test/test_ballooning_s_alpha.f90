program test_ballooning_s_alpha
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use milca_ballooning_s_alpha
    implicit none
    real(dp) :: a, b, c, beta, err, beta2, err2
    logical :: crossed
    integer :: rc
    ! The zero-drive marginal equation has F=1: cannot have a zero.
    call ballooning_first_zero(1.0_dp,0.0_dp,24*acos(-1.0_dp),1536,crossed,rc)
    if (rc /= BALLOON_OK .or. crossed) error stop 'zero-drive analytically stable'
    ! A strong pressure gradient is unstable (a first zero of the even mode).
    call ballooning_first_zero(1.0_dp,1.0_dp,24*acos(-1.0_dp),1536,crossed,rc)
    if (rc /= BALLOON_OK .or. .not. crossed) error stop 'driven mode must cross'
    call ballooning_critical_alpha(1.0_dp,2.0_dp,16*acos(-1.0_dp),1024,a,rc)
    if (rc /= BALLOON_OK) error stop 'first stability boundary'
    call ballooning_critical_alpha(1.0_dp,2.0_dp,16*acos(-1.0_dp),2048,b,rc)
    if (rc /= BALLOON_OK) error stop 'step refinement'
    call ballooning_critical_alpha(1.0_dp,2.0_dp,32*acos(-1.0_dp),2048,c,rc)
    if (rc /= BALLOON_OK) error stop 'domain refinement'
    ! Independent published s-alpha first-stability diagram (s=1, alpha about 0.6).
    if (c < 0.58_dp .or. c > 0.65_dp) error stop 's-alpha first boundary'
    if (abs(a-b) > 0.003_dp) error stop 'RK4 convergence'
    if (abs(b-c) > 0.01_dp) error stop 'angle horizon sensitivity'
    call ballooning_beta_limit(1.0_dp,2.0_dp,0.8_dp,1.0_dp,3.0_dp, &
        beta,err,rc)
    if (rc /= BALLOON_OK) error stop 'local beta calculation'
    if (abs(beta-c/(4*3.0_dp*4*0.8_dp)) > 2e-10_dp) &
        error stop 'beta pressure profile factor'
    if (err <= 0 .or. err > 0.02_dp*beta) error stop 'error estimate'
    call ballooning_beta_limit(1.0_dp,4.0_dp,0.8_dp,1.0_dp,3.0_dp, &
        beta2,err2,rc)
    if (rc /= BALLOON_OK) error stop 'second q run'
    if (abs(4*beta2/beta-1) > 1e-10_dp) error stop 'beta scales like q^-2'
    call ballooning_beta_limit(0.0_dp,2.0_dp,0.8_dp,1.0_dp,3.0_dp, &
        beta,err,rc)
    if (rc /= BALLOON_INVALID) error stop 'invalid shear accepted'
    print '(a,3f12.7)', 'alpha boundary; coarse, refined, long ', a,b,c
    print '(a)', 'PASS s-alpha ballooning beta model'
end program
