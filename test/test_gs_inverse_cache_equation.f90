program inverse_cache_parent_polynomial
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_grad_shafranov, only: gs_state_t, gs_profile_implicit_response, &
        gs_element_jet
    implicit none
    type(gs_state_t)::s, derivative
    real(dp), allocatable::response(:)
    real(dp) :: flux(3), ff(3), dff(3), axis(2), daxis, grad(2)
    real(dp) :: hess(2, 2), area, residual
    integer::i, status
    allocate (s%points(2, 6), s%triangles(6, 1), s%psi(6), s%boundary(6), &
        s%stiffness(6, 6, 1), s%gradient(2, 3, 1), s%area(1), s%diagonal(6))
    s%points = reshape([2.0_dp, 0.0_dp, 3.0_dp, 0.0_dp, 2.0_dp, 1.0_dp, &
        2.5_dp, 0.0_dp, 2.5_dp, 0.5_dp, 2.0_dp, 0.5_dp], [2, 6])
    s%triangles(:, 1) = [1, 2, 3, 4, 5, 6]; s%element_order = 2
    s%gradient(:, :, 1) = reshape([-1.0_dp, -1.0_dp, 1.0_dp, 0.0_dp, &
        0.0_dp, 1.0_dp], [2, 3])
    s%boundary = .false.; s%area = .5_dp; s%diagonal = 1; s%stiffness = 0
    do i = 1, 6
        s%stiffness(i, i, 1) = 1
    end do
    axis = [7.0_dp/3, 1.0_dp/3]
    s%psi = (s%points(1, :) - axis(1))**2 + (s%points(2, :) - axis(2))**2 - .01_dp
    flux = [-.01_dp, -.005_dp, 1.0_dp]; ff = [0.0_dp, .005_dp, .005_dp]; dff = 1
    call gs_profile_implicit_response(s, flux, ff, dff, response, status)
    print *, 'implicit finite-dimensional response status', status
    if (status /= 0) error stop 'response unavailable'
    derivative = s; derivative%psi = response
    call gs_element_jet(derivative, 1, axis(1), axis(2), daxis, grad, hess, status)
    if (status /= 0) error stop 'response quadratic unavailable'
    ! Independent disk moments; the active slope disk has radius^2=.005.
    ! Integral of quadratic u is area*(u_axis+r^2*trace(Hu)/8).
    ! Integral psi/psi_axis over this disk is .75*area.
    ! Testing the weak source against represented R cancels its 1/R exactly.
    area = acos(-1.0_dp)*.005_dp
    residual = dot_product(s%points(1, :), response) - &
        area*(.25_dp*daxis + .005_dp*(hess(1, 1) + hess(2, 2))/8) - .5_dp
    print *, 'independent exact disk moment equation residual', residual
    print *, 'scope polynomial weak response with identity operator'
    print *, 'no equilibrium solve'
    if (abs(residual) > 1e-9_dp) error stop 'profile cut mass/axis response alias'
end program inverse_cache_parent_polynomial
