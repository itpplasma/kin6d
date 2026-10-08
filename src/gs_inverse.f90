module kin6d_gs_inverse
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use kin6d_grad_shafranov, only: gs_state_t, solve_gs, gs_axis_candidate, &
        gs_element_jet, gs_profile_response, gs_profile_implicit_response, &
        gs_profile_response_cache_t, prepare_gs_profile_response_cache, &
        gs_profile_value, gs_profile_integral
    use kin6d_gs_coarea, only: gs_level_coarea
    use kin6d_gs_inverse_leaves, only: gs_f_from_coarea, gs_square_slope, gs_f_from_square, &
        gs_f_from_square_jvp, gs_profile_square_response, gs_q_log_residual, gs_q_log_residual_jvp, &
        gs_profile_mean_last, gs_profile_mean_last_jvp, gs_ff_zero
    use fortnum_linalg, only: dense_solve
    use fortfem_triangle_duffy_quadrature, only: triangle_duffy_quadrature
    implicit none
    private
    public :: gs_inverse_result_t, solve_gs_constant_q, gs_inverse_f, gs_toroidal_flux, gs_inverse_profile_status
    type :: gs_inverse_result_t
        real(dp), allocatable :: flux(:), f_squared(:), ffprime(:)
        real(dp) :: q_target = 0, phi_target = 0, phi_native = 0
        real(dp) :: positive_normalization = 1, coupled_residual = huge(1.0_dp)
        real(dp) :: angular_indicator = huge(1.0_dp), profile_indicator = huge(1.0_dp)
        real(dp) :: flux_indicator = huge(1.0_dp), knot_indicator = huge(1.0_dp)
        real(dp) :: jacobian_check_error = huge(1.0_dp)
        real(dp) :: half_indicator = huge(1.0_dp), validation_indicator = huge(1.0_dp)
        integer :: iterations = 0, angular_points = 0
        character(32) :: claim = 'candidate'
    end type gs_inverse_result_t
contains
    subroutine solve_gs_constant_q(state, q_target, phi_target, result, status, &
            flux_points, angular_points, relative_tolerance, q_tolerance, flux_tolerance, verify_jacobian)
        type(gs_state_t), intent(inout) :: state
        real(dp), intent(in) :: q_target, phi_target
        type(gs_inverse_result_t), intent(out) :: result
        integer, intent(out) :: status
        integer, intent(in), optional :: flux_points, angular_points
        real(dp), intent(in), optional :: relative_tolerance, q_tolerance, flux_tolerance
        logical,intent(in),optional::verify_jacobian
        type(gs_state_t) :: trial_state
        type(gs_inverse_result_t) :: trial_result
        real(dp), allocatable :: x(:), xp(:), residual(:), rp(:), rm(:), jac(:,:), step(:)
        real(dp), allocatable :: coarea(:), refined(:), solution(:), half_levels(:)
        real(dp)::validation_levels(128),validation_coarse(129),validation_fine(129),ignored_indicator
        real(dp) :: axis(2), axis_flux, h(2,2), f, two_pi, tolerance, factor, phi
        real(dp) :: length, weak_residual, objective, qtolerance, ptolerance, phi_coarse, check_error, delta
        integer :: n, angles, stat, element, i, iteration, backtrack
        check_error=huge(1.0_dp)
        status=1; n=9; angles=256; tolerance=1.0e-8_dp; qtolerance=1.0e-3_dp; ptolerance=1.0e-6_dp
        if(present(flux_points)) n=flux_points
        if(present(angular_points)) angles=angular_points
        if(present(relative_tolerance)) tolerance=relative_tolerance
        if(present(q_tolerance))qtolerance=q_tolerance
        if(present(flux_tolerance))ptolerance=flux_tolerance
        if(.not.all(ieee_is_finite([q_target,phi_target,tolerance,qtolerance,ptolerance])))return
        if(abs(q_target)<=tiny(q_target).or.abs(phi_target)<=tiny(phi_target))return
        if(phi_target/q_target<=0.or.tolerance<=0.or.tolerance>=1)return
        if(qtolerance<=0.or.qtolerance>=1.or.ptolerance<=0.or.ptolerance>=1)return
        if(state%element_order/=2.or.n<5.or.n>65.or.angles<16.or.2*angles>16384)return
        if(any(abs(pack(state%psi,state%boundary))>0))return
        two_pi=2*acos(-1.0_dp)
        if(maxval(abs(state%psi))<=tiny(1.0_dp))then
            call solve_gs(state,1.0_dp,0.0_dp,-1.0_dp,stat)
            if(stat/=0)return
        end if
        call gs_axis_candidate(state,axis,axis_flux,h,element,stat)
        if(stat/=0)then
            status=stat
            return
        end if
        allocate(x(n),xp(n),residual(n),rp(n),rm(n),jac(n,n),step(n),coarea(n),refined(n))
        allocate(result%flux(n),result%f_squared(n),result%ffprime(n))
        do i=1,n
            result%flux(i)=axis_flux*(1-real(i-1,dp)/(n-1))
        end do
        call gs_level_coarea(state,result%flux,angles,coarea,axis,stat)
        if(stat/=0)then
            status=100+stat
            return
        end if
        call gs_f_from_coarea(q_target,two_pi,coarea(n),f)
        x(:n-1)=-1; x(n)=log(f*f)
        ! Fix the normalized-flux mean FF' to -1 to remove positive zero-beta
        ! homogeneity over the complete profile. Phi selects the physical
        ! amplitude. The final coefficient follows the exact linear-profile mean.
        call profile_candidate(state,x,q_target,angles,tolerance,result,residual,stat)
        if(stat/=0)then
            status=200+stat
            return
        end if
        do iteration=1,24
            objective=maxval(abs(residual))
            result%iterations=iteration
            if(objective<=min(1.0e-6_dp,10*tolerance))exit
            call profile_jacobian(state,result,angles,jac,stat)
            if(stat/=0)then
                status=300+stat
                return
            end if
            if(present(verify_jacobian))then
                if(verify_jacobian.and.iteration==1)then
                    check_error=0
                    do i=1,n
                        delta=1.0e-5_dp*max(1.0_dp,abs(x(i)))
                        xp=x; xp(i)=x(i)+delta; trial_state=state
                        call profile_candidate(trial_state,xp,q_target,angles,tolerance,trial_result,rp,stat)
                        if(stat/=0)return
                        xp=x; xp(i)=x(i)-delta; trial_state=state
                        call profile_candidate(trial_state,xp,q_target,angles,tolerance,trial_result,rm,stat)
                        if(stat/=0)return
                        check_error=max(check_error,maxval(abs((rp-rm)/(2*delta)-jac(:,i))) &
                            /max(maxval(abs(jac(:,i))),1.0e-8_dp))
                    end do
                    if(check_error>1.0e-3_dp)then
                        status=11
                        return
                    end if
                end if
            end if
            call dense_solve(jac,-residual,step,stat)
            if(stat/=0)then
                status=3
                return
            end if
            length=1
            do backtrack=1,12
                xp=x+length*step; trial_state=state
                call profile_candidate(trial_state,xp,q_target,angles,tolerance, &
                    trial_result,rp,stat)
                if(stat==0)then
                    if(maxval(abs(rp))<objective)exit
                end if
                length=length/2
            end do
            if(backtrack>12)then
                status=4
                return
            end if
            x=xp; state=trial_state; result=trial_result; residual=rp
            result%iterations=iteration
        end do
        if(maxval(abs(residual))>min(1.0e-6_dp,10*tolerance))then
            status=2
            return
        end if
        result%jacobian_check_error=check_error
        result%knot_indicator=maxval(abs(exp(residual)-1))
        call gs_level_coarea(state,result%flux,angles,coarea,axis,stat)
        if(stat/=0)return
        call gs_level_coarea(state,result%flux,2*angles,refined,axis,stat)
        if(stat/=0)return
        result%angular_indicator=maxval(abs(refined/coarea-1))
        half_levels=.5_dp*(result%flux(:n-1)+result%flux(2:))
        call check_half_nodes(state,result,half_levels,2*angles,result%half_indicator,stat)
        if(stat/=0)return
        ! A fixed rho_pol grid is independent of the fitted flux knots and
        ! includes small near-axis surfaces and the near-boundary interval.
        do i=1,128
            validation_levels(i)=result%flux(1)*(1-(real(i,dp)/128)**2)
        end do
        call check_half_nodes(state,result,validation_levels,angles,ignored_indicator, &
            stat,validation_coarse)
        if(stat/=0)return
        call check_half_nodes(state,result,validation_levels,2*angles,result%validation_indicator, &
            stat,validation_fine)
        if(stat/=0)return
        result%profile_indicator=max(result%half_indicator,result%validation_indicator)
        result%angular_indicator=max(result%angular_indicator, &
            maxval(abs(validation_fine/validation_coarse-1)))
        call gs_toroidal_flux(state,result,phi_coarse,stat,8)
        if(stat/=0)return
        call gs_toroidal_flux(state,result,phi,stat,12)
        if(stat/=0)return
        if(phi/q_target<=0)return
        result%flux_indicator=abs(phi_coarse/phi-1)
        if(result%flux_indicator>ptolerance)then
            status=10
            return
        end if
        factor=phi_target/phi
        state%psi=state%psi*factor; result%flux=result%flux*factor
        result%f_squared=result%f_squared*factor*factor; result%ffprime=result%ffprime*factor
        result%positive_normalization=factor; result%phi_native=phi; result%phi_target=phi_target
        call gs_profile_response(state,result%flux,result%ffprime,solution,weak_residual,stat)
        if(stat/=0)return
        result%coupled_residual=weak_residual
        if(weak_residual>tolerance)then
            status=2
            return
        end if
        state%residual=weak_residual; state%iterations=result%iterations
        ! Complete the declared zero-pressure amplitude constraint before
        ! returning an underresolved q candidate. Failure remains status9;
        ! diagnostics then refer to the requested physical Phi, not a seed.
        if(result%knot_indicator+result%profile_indicator+result%angular_indicator>qtolerance)then
            status=9
            return
        end if
        status=0
    end subroutine solve_gs_constant_q

    subroutine profile_candidate(state,x,q_target,angles,tolerance,result,residual,status)
        type(gs_state_t),intent(inout)::state
        real(dp),intent(in)::x(:),q_target,tolerance
        integer,intent(in)::angles
        type(gs_inverse_result_t),intent(out)::result
        real(dp),intent(out)::residual(:)
        integer,intent(out)::status
        real(dp),allocatable::solution(:),coarea(:)
        real(dp)::axis(2),axis_flux,h(2,2),weak_residual,f,square
        integer::n,i,stat,element,iteration
        status=1; residual=huge(1.0_dp); n=size(x)
        if(.not.all(ieee_is_finite(x)).or.abs(x(n))>100)return
        allocate(result%flux(n),result%ffprime(n),result%f_squared(n),coarea(n))
        result%flux=0; result%ffprime=0; result%f_squared=0
        result%q_target=q_target; result%angular_points=angles
        result%ffprime(:n-1)=x(:n-1)
        call gs_profile_mean_last(-1.0_dp, &
            (.5_dp*x(1)+sum(x(2:n-1)))/real(n-1,dp), &
            .5_dp/real(n-1,dp),result%ffprime(n))
        do iteration=1,150
            call gs_axis_candidate(state,axis,axis_flux,h,element,stat)
            if(stat/=0)then
                status=10+stat
                return
            end if
            do i=1,n
                result%flux(i)=axis_flux*(1-real(i-1,dp)/(n-1))
            end do
            call gs_profile_response(state,result%flux,result%ffprime,solution,weak_residual,stat)
            if(stat/=0)return
            if(weak_residual<=min(1.0e-10_dp,.01_dp*tolerance))exit
            state%psi=state%psi+.6_dp*(solution-state%psi)
        end do
        if(iteration>150)then
            status=2
            return
        end if
        result%coupled_residual=weak_residual
        do i=1,n
            square=exp(x(n))+2*gs_profile_integral(result%flux,result%ffprime,result%flux(i))
            if(square<=0)return
            result%f_squared(i)=square
        end do
        call gs_inverse_profile_status(result,stat)
        if(stat/=0)then
            status=15
            return
        end if
        call gs_level_coarea(state,result%flux,angles,coarea,axis,stat)
        if(stat/=0)then
            status=30+stat
            return
        end if
        do i=1,n
            call gs_inverse_f(result,result%flux(i),f,stat)
            if(stat/=0)return
            call gs_q_log_residual(f,coarea(i),2*acos(-1.0_dp),q_target,residual(i))
        end do
        status=0
    end subroutine profile_candidate

    subroutine profile_jacobian(state,result,angles,jac,status)
        type(gs_state_t),intent(in)::state
        type(gs_inverse_result_t),intent(in)::result
        integer,intent(in)::angles
        real(dp),intent(out)::jac(:,:)
        integer,intent(out)::status
        type(gs_state_t)::tangent
        type(gs_profile_response_cache_t) :: cache
        real(dp),allocatable::dff(:),dpsi(:),levels_tangent(:),coarea(:),coarea_tangent(:)
        real(dp)::axis(2),axis_flux,h(2,2),value_tangent,grad(2),ht(2,2),square_tangent,f,df
        real(dp)::integral_tangent,edge_tangent,partial_mean,last_weight,partial_tangent
        integer::n,i,j,element,stat
        status=1; n=size(result%flux)
        if(any(shape(jac)/=[n,n]))return
        call gs_axis_candidate(state,axis,axis_flux,h,element,stat)
        if(stat/=0)return
        allocate(dff(n),levels_tangent(n),coarea(n),coarea_tangent(n))
        tangent=state
        partial_mean=(.5_dp*result%ffprime(1)+sum(result%ffprime(2:n-1)))/real(n-1,dp)
        last_weight=.5_dp/real(n-1,dp)
        call prepare_gs_profile_response_cache(state, result%flux, &
            result%ffprime, cache, stat)
        if(stat/=0)return
        do j=1,n
            dff=0
            if(j<n)then
                dff(j)=1
                partial_tangent=merge(.5_dp,1.0_dp,j==1)/real(n-1,dp)
                call gs_profile_mean_last_jvp(-1.0_dp,partial_mean,last_weight, &
                    0.0_dp,partial_tangent,0.0_dp,dff(n))
                call gs_profile_implicit_response(state, result%flux, result%ffprime, &
                    dff, dpsi, stat, cache)
                if(stat/=0)return
            else
                if(allocated(dpsi))deallocate(dpsi)
                allocate(dpsi(size(state%psi))); dpsi=0
            end if
            tangent%psi=dpsi
            call gs_element_jet(tangent,element,axis(1),axis(2),value_tangent,grad,ht,stat)
            if(stat/=0)return
            levels_tangent=result%flux*value_tangent/axis_flux
            call gs_level_coarea(state,result%flux,angles,coarea,axis,stat, &
                dpsi,levels_tangent,coarea_tangent)
            if(stat/=0)then
                status=30+stat
                return
            end if
            edge_tangent=0
            if(j==n)edge_tangent=result%f_squared(n)
            do i=1,n
                integral_tangent=gs_profile_integral(result%flux,dff,result%flux(i))
                call gs_profile_square_response(result%f_squared(i),result%f_squared(n), &
                    axis_flux,value_tangent,integral_tangent,edge_tangent,square_tangent)
                call gs_f_from_square(result%f_squared(i),sign(1.0_dp,result%q_target),f)
                call gs_f_from_square_jvp(result%f_squared(i),sign(1.0_dp,result%q_target), &
                    square_tangent,0.0_dp,df)
                call gs_q_log_residual_jvp(f,coarea(i),2*acos(-1.0_dp),result%q_target, &
                    df,coarea_tangent(i),0.0_dp,0.0_dp,jac(i,j))
            end do
        end do
        status=0
    end subroutine profile_jacobian

    subroutine gs_inverse_profile_status(result, status)
        ! Positivity over the complete closed axis-to-edge interval. With
        ! piecewise-linear FF', F^2 is quadratic; test its interior minima.
        type(gs_inverse_result_t), intent(in)::result
        integer, intent(out)::status
        real(dp)::critical, square
        integer::n, i
        status=1
        if (.not.allocated(result%flux)) return
        if (.not.allocated(result%ffprime)) return
        if (.not.allocated(result%f_squared)) return
        if (.not.ieee_is_finite(result%q_target)) return
        if (abs(result%q_target)<=tiny(1.0_dp)) return
        n=size(result%flux)
        if (n<2) return
        if (size(result%ffprime)/=n.or.size(result%f_squared)/=n) return
        if (.not.all(ieee_is_finite(result%flux))) return
        if (.not.all(ieee_is_finite(result%ffprime))) return
        if (.not.all(ieee_is_finite(result%f_squared))) return
        if (any(result%flux(2:)<=result%flux(:n-1))) return
        do i=1,n
            square=result%f_squared(n)+2*gs_profile_integral(result%flux,result%ffprime,result%flux(i))
            if (.not.ieee_is_finite(square)) return
            if (square<=0) return
        end do
        do i=1,n-1
            if (result%ffprime(i)>=0.or.result%ffprime(i+1)<=0) cycle
            call gs_ff_zero(result%flux(i),result%flux(i+1),result%ffprime(i), &
                result%ffprime(i+1),critical)
            square=result%f_squared(n)+2*gs_profile_integral(result%flux,result%ffprime,critical)
            if (.not.ieee_is_finite(square)) return
            if (square<=0) return
        end do
        status=0
    end subroutine gs_inverse_profile_status

    subroutine gs_inverse_f(result, psi, f, status)
        type(gs_inverse_result_t), intent(in) :: result
        real(dp), intent(in) :: psi
        real(dp), intent(out) :: f
        integer, intent(out) :: status
        real(dp) :: square
        integer :: n
        f = 0; status = 1
        if (.not. allocated(result%flux)) return
        if (.not. allocated(result%f_squared)) return
        if (.not. allocated(result%ffprime)) return
        if (.not. ieee_is_finite(psi)) return
        if (.not.ieee_is_finite(result%q_target)) return
        if (abs(result%q_target) <= tiny(1.0_dp)) return
        n = size(result%flux)
        if (n < 2) return
        if (size(result%f_squared) /= n .or. size(result%ffprime) /= n) return
        square=result%f_squared(n)+2*gs_profile_integral(result%flux,result%ffprime,psi)
        if (.not. ieee_is_finite(square)) return
        if (square <= 0) return
        call gs_f_from_square(square, sign(1.0_dp, result%q_target), f)
        status = 0
    end subroutine gs_inverse_f

    subroutine gs_toroidal_flux(state, result, phi, status, polynomial_degree)
        type(gs_state_t), intent(in) :: state
        type(gs_inverse_result_t), intent(in) :: result
        real(dp), intent(out) :: phi
        integer, intent(out) :: status
        integer, intent(in), optional :: polynomial_degree
        real(dp), allocatable :: xi(:), eta(:), weights(:)
        real(dp) :: point(2), value, grad(2), hessian(2, 2), f
        integer :: e, i, stat, degree
        phi = 0; status = 1; degree = 8
        if (present(polynomial_degree)) degree = polynomial_degree
        if (degree < 2 .or. degree > 30) return
        call triangle_duffy_quadrature(degree, xi, eta, weights, stat)
        if (stat /= 0) return
        do e = 1, size(state%triangles, 2)
            do i = 1, size(xi)
                point = matmul(state%points(:, state%triangles(1:3, e)), &
                    [1-xi(i)-eta(i), xi(i), eta(i)])
                call gs_element_jet(state, e, point(1), point(2), &
                    value, grad, hessian, stat)
                if (stat /= 0) return
                call gs_inverse_f(result, value, f, stat)
                if (stat /= 0) return
                phi = phi+2*state%area(e)*weights(i)*f/point(1)
            end do
        end do
        if (.not. ieee_is_finite(phi)) return
        status = 0
    end subroutine gs_toroidal_flux

    subroutine check_half_nodes(state, result, half_levels, angles, indicator, status, coarea_out)
        type(gs_state_t), intent(in) :: state
        type(gs_inverse_result_t), intent(in) :: result
        real(dp), intent(in) :: half_levels(:)
        integer, intent(in) :: angles
        real(dp), intent(out) :: indicator
        integer, intent(out) :: status
        real(dp),optional,intent(out)::coarea_out(:)
        real(dp), allocatable :: levels(:), coarea(:)
        real(dp) :: axis(2), f, q
        integer :: i, stat
        indicator = huge(1.0_dp); status = 1
        allocate (levels(size(half_levels)+1), coarea(size(half_levels)+1))
        levels(1) = result%flux(1); levels(2:) = half_levels
        call gs_level_coarea(state, levels, angles, coarea, axis, stat)
        if (stat /= 0) return
        if(present(coarea_out))then
            if(size(coarea_out)/=size(coarea))return
            coarea_out=coarea
        end if
        indicator = 0
        do i = 1, size(half_levels)
            call gs_inverse_f(result, half_levels(i), f, stat)
            if (stat /= 0) return
            q = f*coarea(i+1)/(2*acos(-1.0_dp))
            indicator = max(indicator, abs(q/result%q_target-1))
        end do
        status = 0
    end subroutine check_half_nodes
end module kin6d_gs_inverse
