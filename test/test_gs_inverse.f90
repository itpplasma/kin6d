program test_gs_inverse
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_grad_shafranov, only: gs_state_t, initialize_gs, gs_axis_candidate, gs_profile_implicit_response
    use kin6d_gs_inverse, only: gs_inverse_result_t, gs_toroidal_flux, &
        solve_gs_constant_q
    implicit none
    type(gs_state_t) :: state, circular, reversed
    type(gs_inverse_result_t) :: profiles, inverse, negative
    real(dp) :: rectangle(2, 4), boundary(2, 63), phi, expected, theta
    real(dp) :: axis(2), axis_flux, hessian(2, 2)
    real(dp),allocatable::response(:)
    integer :: status, i, element
    character(16)::mode
    integer::knots
    real(dp)::area
    call get_command_argument(1,mode)
    if(len_trim(mode)==0)mode="coarse"
    knots=5; area=.0025_dp
    if(trim(mode)=="fine")then
        knots=9; area=.00125_dp
    end if
    rectangle = reshape([5.8_dp, -0.4_dp, 6.6_dp, -0.4_dp, &
        6.6_dp, 0.4_dp, 5.8_dp, 0.4_dp], [2, 4])
    call initialize_gs(state, rectangle, 0.008_dp, status, 2)
    if (status /= 0) error stop 'inverse flux-oracle mesh failed'
    profiles%q_target = -1
    profiles%flux = [-1.0_dp, 1.0_dp]
    profiles%f_squared = [900.0_dp, 900.0_dp]
    profiles%ffprime = [0.0_dp,0.0_dp]
    call gs_toroidal_flux(state, profiles, phi, status, 12)
    if (status /= 0) error stop 'native physical toroidal flux integral failed'
    expected = -30.0_dp*0.8_dp*log(6.6_dp/5.8_dp)
    if (abs(phi/expected-1) > 2.0e-12_dp) &
        error stop 'independent signed rectangular toroidal flux oracle'
    do i = 1, size(boundary, 2)
        theta = 2*acos(-1.0_dp)*(i-1)/size(boundary, 2)
        boundary(:, i) = [4.0_dp+0.4_dp*cos(theta), 0.4_dp*sin(theta)]
    end do
    call initialize_gs(circular, boundary, area, status, 2)
    if (status /= 0) error stop 'constant-q unit mesh failed'
    call solve_gs_constant_q(circular, 1.5_dp, 1.0_dp, inverse, status, &
        flux_points=knots, angular_points=128, &
        relative_tolerance=1.0e-6_dp,q_tolerance=.005_dp, &
        verify_jacobian=(trim(mode)/="fine"))
    print *, 'constant-q status/iterations/GS residual:', status, inverse%iterations, &
        inverse%coupled_residual
    print *, 'constant-q Jacobian independent check:', inverse%jacobian_check_error
    print *, 'constant-q angular/profile indicators:', inverse%angular_indicator, &
        inverse%profile_indicator
    call gs_axis_candidate(circular, axis, axis_flux, hessian, element, i)
    print *, 'native axis status/position/flux:', i, axis, axis_flux
    print *, 'native Hessian/asymmetry:', hessian, abs(hessian(1,2)-hessian(2,1)), &
        64*epsilon(1.0_dp)*maxval(abs(hessian))
    if (status /= 0) error stop 'constant-q coupled GS unit problem failed'
    if (inverse%coupled_residual > 1.0e-6_dp) error stop 'coupled weak GS residual failed'
    if (inverse%profile_indicator > 0.01_dp) error stop 'off-node q profile failed'
    call gs_toroidal_flux(circular, inverse, phi, status, 12)
    if (status /= 0) error stop 'normalized independent Phi check failed'
    if (abs(phi-1) > 2.0e-5_dp) error stop 'fixed signed toroidal flux failed'
    call gs_profile_implicit_response(circular,inverse%flux,inverse%ffprime, &
        inverse%ffprime,response,status)
    if(status/=0)error stop 'native implicit homogeneous response failed'
    if(norm2(response-circular%psi)/norm2(circular%psi)>2.0e-7_dp) &
        error stop 'independent zero-beta amplitude response oracle'
    if(trim(mode)/='fine')then
        reversed=circular
        call solve_gs_constant_q(reversed,-1.5_dp,-1.0_dp,negative,status, &
            flux_points=knots,angular_points=128,relative_tolerance=1.0e-6_dp,q_tolerance=.005_dp)
        if(status/=0)error stop 'signed q/Phi reversal control failed'
        ! Known limitation: the nonlinear solve is converged only to about 1e-7 in a
        ! -march=x86-64-v3 build (5 instead of 3 Newton steps; tolerances 1e-7 to
        ! 1e-6 give identical output and 1e-8 is not attainable), so there the
        ! psi difference is 5.2e-8 and the ffprime difference 7.4e-6 of max|ffprime|,
        ! above these gates. The registered CTest uses default flags: 1.6e-11 and 4e-10.
        if(maxval(abs(reversed%psi-circular%psi))>2.0e-8_dp) &
            error stop 'F-only reversal changed poloidal flux'
        if(maxval(abs(negative%ffprime-inverse%ffprime))>2.0e-5_dp) &
            error stop 'F-only reversal changed GS toroidal-current source'
        call gs_toroidal_flux(reversed,negative,phi,status,12)
        if(status/=0.or.abs(phi+1)>2.0e-6_dp)error stop 'negative physical Phi control'
        call solve_gs_constant_q(reversed,-1.5_dp,1.0_dp,negative,status)
        if(status==0)error stop 'inconsistent signed q/Phi accepted'
        call solve_gs_constant_q(reversed,1.5_dp,1.0_dp,negative,status, &
            flux_points=5,angular_points=128,relative_tolerance=1.0e-6_dp,q_tolerance=1.0e-4_dp)
        if(status/=9)error stop 'underresolved off-knot q accepted as converged'
        call gs_toroidal_flux(reversed,negative,phi,status,12)
        if(status/=0.or.abs(phi-1)>2.0e-6_dp) &
            error stop 'underresolved candidate lost the declared physical Phi'
    end if
    print *, 'PASS: native coupled constant-q/zero-beta/Phi unit controls'
end program test_gs_inverse
