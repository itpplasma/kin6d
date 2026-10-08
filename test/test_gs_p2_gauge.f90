program test_gs_p2_gauge
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_grad_shafranov, only: gs_state_t, initialize_gs, gs_element_jet, &
        gs_classical_current, gs_axis_candidate, gs_axis_safety_factor
    use kin6d_grad_shafranov_leaves, only: gs_smooth_axis_q_jvp, gs_current_density_jvp
    implicit none
    type(gs_state_t) :: state
    real(dp) :: boundary(2,4), point(2), center(2), h(2,2), grad(2), value, current, q
    real(dp) :: wanted_h(2,2), displacement(2), gauge, rounding, jet_bound, scale, derivative
    real(dp), parameter :: gauges(3)=[0.0_dp,1.0_dp,2.0_dp**20]
    integer :: stat, i, e, kind, element, failures, count
    failures=0; count=0
    boundary=reshape([5.8_dp,-.4_dp,6.6_dp,-.4_dp,6.6_dp,.4_dp,5.8_dp,.4_dp],[2,4])
    call initialize_gs(state,boundary,.012_dp,stat,2)
    if(stat/=0) error stop 'gauge oracle mesh failed'
    ! An exactly represented constant has exactly zero physical derivatives.
    state%psi=2.0_dp**30
    do e=1,size(state%triangles,2)
        point=sum(state%points(:,state%triangles(1:3,e)),dim=2)/3
        call gs_element_jet(state,e,point(1),point(2),value,grad,h,stat)
        call require(stat==0,'constant jet status')
        call require(all(grad==0).and.all(h==0),'exact constant derivatives')
        call gs_classical_current(state,e,point(1),point(2),1.0_dp,current,stat)
        call require(stat==0.and.current==0,'exact constant classical current')
        ! The physical jet is linear in nodal flux. A constant nodal tangent
        ! must therefore induce zero current and axis-q through generated JVPs.
        call gs_current_density_jvp(point(1),1.0_dp,0.0_dp,2.0_dp,3.0_dp, &
             0.0_dp,0.0_dp,grad(1),h(1,1),h(2,2),derivative)
        call require(derivative==0,'constant-gauge physical current JVP')
        call gs_smooth_axis_q_jvp(point(1),-30.0_dp,2.0_dp,.5_dp,3.0_dp, &
             0.0_dp,0.0_dp,h(1,1),h(1,2),h(2,2),derivative)
        call require(derivative==0,'constant-gauge smooth axis-q JVP')
    end do
    center=sum(state%points(:,state%triangles(1:3,1)),dim=2)/3
    wanted_h=reshape([2.0_dp,.5_dp,.5_dp,3.0_dp],[2,2])
    do kind=1,3
        gauge=gauges(kind)
        do i=1,size(state%psi)
            displacement=state%points(:,i)-center
            state%psi(i)=gauge+dot_product(displacement,matmul(wanted_h,displacement))/2
        end do
        do e=1,size(state%triangles,2)
            point=sum(state%points(:,state%triangles(1:3,e)),dim=2)/3
            ! A conservative physical-jet amplification bound on binary64 nodal
            ! input rounding; no claim of exact arbitrary-gauge representability.
            rounding=spacing(max(1.0_dp,maxval(abs(state%psi(state%triangles(:,e))))))
            scale=maxval(abs(state%gradient(:,:,e)))
            jet_bound=128*rounding*max(1.0_dp,scale**2)
            call gs_element_jet(state,e,point(1),point(2),value,grad,h,stat)
            call require(stat==0,'quadratic gauge jet status')
            call require(maxval(abs(grad-matmul(wanted_h,point-center)))<=jet_bound, &
                         'quadratic physical gradient with input rounding')
            call require(maxval(abs(h-wanted_h))<=jet_bound, &
                         'quadratic Hessian with input rounding')
            call gs_classical_current(state,e,point(1),point(2),2.0_dp,current,stat)
            displacement=matmul(wanted_h,point-center)
            call require(abs(current+(5-displacement(1)/point(1))/(2*point(1)))<=jet_bound, &
                         'quadratic physical current with input rounding')
        end do
        call gs_axis_candidate(state,point,value,h,element,stat)
        call require(stat==0,'quadratic interior axis survives gauge')
        if(stat/=0) cycle
        call gs_axis_safety_factor(point(1),-30.0_dp,h,q,stat)
        call require(stat==0,'quadratic axis-q guard survives gauge')
        call require(abs(q+30/(center(1)*sqrt(5.75_dp)))<=1.0e-5_dp, &
                     'independent elliptic coarea axis-q with input rounding')
    end do
    print *, 'gauge controls/failures:',count,failures
    if(failures/=0) error stop 'physical P2 gauge stability failed'
contains
    subroutine require(condition,label)
        logical,intent(in)::condition
        character(*),intent(in)::label
        count=count+1
        if(condition)return
        failures=failures+1
        if(failures<=5) print *, 'FAIL: ',label
    end subroutine
end program
