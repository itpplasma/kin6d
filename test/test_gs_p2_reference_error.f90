program test_gs_p2_reference_error
    use, intrinsic :: iso_fortran_env, only: dp=>real64
    use fortnum_interval, only: interval_t, interval, log, operator(+), &
        operator(-), operator(*), operator(/)
    use kin6d_gs_reference_error, only: gs_reference_parameters_t, &
        enclose_gs_reference_cell
    implicit none
    type(gs_reference_parameters_t) :: parameters
    type(interval_t) :: integral(3), coarse(3), area, volume, peak, exact_b, exact_e
    real(dp) :: vertices(2,3), nodes(6)
    integer :: checks
    logical :: ok
    checks=0
    vertices=reshape([1.0_dp,0.0_dp,2.0_dp,0.0_dp,1.0_dp,1.0_dp],[2,3])
    parameters=gs_reference_parameters_t(1.0_dp,0.0_dp,0.0_dp,0.0_dp,1.0_dp,0.0_dp)
    ! Genuine quadratic R^2: area field integral 2, weighted energy 8/3,
    ! weighted Maxwell-stress difference 8. No production shape in oracle.
    nodes=[1.0_dp,4.0_dp,1.0_dp,2.25_dp,2.25_dp,1.0_dp]
    call enclose_gs_reference_cell(vertices,nodes,parameters,0,coarse,area,volume,peak,ok)
    call require(ok,'coarse quadratic supported')
    call enclose_gs_reference_cell(vertices,nodes,parameters,3,integral,area,volume,peak,ok)
    call require(ok,'refined quadratic supported')
    call enclosed(integral(1),interval(2))
    call enclosed(integral(2),interval(8)/3)
    call enclosed(integral(3),interval(8))
    call require(integral(2)%hi-integral(2)%lo < coarse(2)%hi-coarse(2)%lo, &
        'refinement tightens whole-cell cover')
    call enclosed(area,interval(1)/2)
    call enclosed(volume,interval(2)/3)
    call enclose_gs_reference_cell(vertices(:,[1,3,2]),nodes([1,3,2,6,5,4]), &
        parameters,3,integral,area,volume,peak,ok)
    call require(ok,'reversed affine P2 winding supported')
    call enclosed(integral(2),interval(8)/3)
    call enclose_gs_reference_cell(vertices,nodes+1048576.0_dp,parameters,3, &
        integral,area,volume,peak,ok)
    call require(ok,'large representable gauge supported')
    call enclosed(integral(1),interval(2))
    call enclosed(integral(2),interval(8)/3)
    ! Same physical quadratic on a translated positive-R triangle.
    vertices(1,:)=vertices(1,:)+32
    nodes=[1089.0_dp,1156.0_dp,1089.0_dp,1122.25_dp,1122.25_dp,1089.0_dp]
    call enclose_gs_reference_cell(vertices,nodes,parameters,3, &
        integral,area,volume,peak,ok)
    call require(ok,'translated physical quadratic supported')
    call enclosed(integral(1),interval(2))
    call enclosed(integral(2),interval(200)/3)
    call enclosed(integral(3),interval(200))
    vertices(1,:)=vertices(1,:)-32
    ! RZ has nonconstant physical gradient. Integrate over R=1..2,Z=0..2-R:
    ! IB=10/3-4log(2); IE=(8/3)log(2)-10/9.
    nodes=[0.0_dp,0.0_dp,1.0_dp,0.0_dp,0.75_dp,0.5_dp]
    exact_b=interval(10)/3-4*log(interval(2))
    exact_e=(interval(8)/3)*log(interval(2))-interval(10)/9
    call enclose_gs_reference_cell(vertices,nodes,parameters,4,integral,area,volume,peak,ok)
    call require(ok,'mixed quadratic supported')
    call enclosed(integral(1),exact_b)
    call enclosed(integral(2),exact_e)
    call enclose_gs_reference_cell(vertices,nodes(:4),parameters,1, &
        integral,area,volume,peak,ok)
    call require(.not. ok,'unsupported DOF count rejected')
    vertices(1,1)=0
    call enclose_gs_reference_cell(vertices,nodes,parameters,1, &
        integral,area,volume,peak,ok)
    call require(.not. ok,'nonpositive R rejected')
    print '(a,i0,a)', 'PASS: ',checks,' independent quadratic/log P2 enclosure controls'
contains
    subroutine require(condition,message)
        logical, intent(in) :: condition
        character(*), intent(in) :: message
        if (.not. condition) then
            print '(a)',message
            error stop 'P2 enclosure independent control failed'
        end if
        checks=checks+1
    end subroutine require
    subroutine enclosed(actual,expected)
        type(interval_t), intent(in) :: actual, expected
        call require(actual%lo<=expected%lo .and. actual%hi>=expected%hi, &
            'exact physical polynomial/log integral outside whole-cell interval')
    end subroutine enclosed
end program test_gs_p2_reference_error
