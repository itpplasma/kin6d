program test_gs_reference_error
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_interval, only: interval_t, interval, log, contains, operator(-), &
        operator(*), operator(/)
    use kin6d_gs_reference_error, only: gs_reference_parameters_t, &
        enclose_gs_reference_cell
    implicit none
    type(gs_reference_parameters_t) :: parameters
    type(interval_t) :: integral(3), fine(3), area, volume, maximum, oracle(3)
    real(dp) :: p(2,3), nodal(3)
    logical :: ok
    integer :: permutation(3), i
    parameters = gs_reference_parameters_t(1.0_dp, 0.0_dp, 0.0_dp, &
                                           0.0_dp, 1.0_dp, 0.0_dp)
    p(:,1)=[1.0_dp,0.0_dp]; p(:,2)=[2.0_dp,0.0_dp]; p(:,3)=[1.0_dp,1.0_dp]
    nodal = [0.0_dp,1.0_dp,0.0_dp]
    ! Independently integrate psi_h=R-1 on Z in [0,2-R]. The exact
    ! field, energy and R-weighted Maxwell-stress-square integrals are:
    oracle = [1-log(interval(2)), 2*log(interval(2))-1, interval(3)/16]
    call enclose_gs_reference_cell(p, nodal, parameters, 0, integral, &
                                    area, volume, maximum, ok)
    if (.not. ok) error stop 'coarse interval failed'
    do i=1,3
        if (.not. contains(integral(i),oracle(i))) error stop 'oracle not enclosed'
    end do
    call enclose_gs_reference_cell(p, nodal, parameters, 5, fine, &
                                    area, volume, maximum, ok)
    if (.not. ok) error stop 'fine interval failed'
    do i=1,3
        if (.not. contains(fine(i),oracle(i))) error stop 'fine oracle not enclosed'
        if (fine(i)%hi >= integral(i)%hi) error stop 'upper bound did not tighten'
        if (fine(i)%lo <= integral(i)%lo) error stop 'lower bound did not tighten'
    end do
    if (.not. contains(area,0.5_dp)) error stop 'area oracle'
    if (.not. contains(volume,2.0_dp/3)) error stop 'volume oracle'
    permutation=[1,3,2]
    call enclose_gs_reference_cell(p(:,permutation), nodal(permutation), &
        parameters, 5, integral, area, volume, maximum, ok)
    if (.not. ok) error stop 'reordered triangle failed'
    do i=1,3
        if (.not. contains(integral(i),oracle(i))) error stop 'winding sign'
    end do
    p(1,1)=0
    call enclose_gs_reference_cell(p,nodal,parameters,0,integral, &
                                    area,volume,maximum,ok)
    if (ok) error stop 'R zero must reject'
    print '(a)', 'PASS: exact field/current-energy/stress integrals; ' &
        //'refinement; winding; R guard'
end program test_gs_reference_error
