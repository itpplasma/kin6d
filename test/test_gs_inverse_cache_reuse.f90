program test_gs_inverse_cache_reuse
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_grad_shafranov, only: gs_state_t, gs_profile_implicit_response, &
        gs_profile_response_cache_t, &
        prepare_gs_profile_response_cache, gs_element_jet
    implicit none
    type(gs_state_t)::s, derivative, changed
    type(gs_profile_response_cache_t)::cache, renewed, empty
    real(dp), allocatable::response(:), original(:), repeated(:)
    real(dp) :: flux(3), ff(3), dff(3), axis(2), daxis, grad(2)
    real(dp) :: hess(2, 2), area, residual
    integer::i, status, passed
    passed = 0
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
    call require(abs(residual) < 1e-9_dp, 'independent exact disk moment equation')
    original = response
    call prepare_gs_profile_response_cache(s, flux, ff, cache, status)
    call require(status == 0, 'prepare one immutable response context')
    call gs_profile_implicit_response(s, flux, ff, dff, response, status, cache)
    call require(status == 0, 'reuse prepared context')
    call require(maxval(abs(response - original)) < 1e-13_dp, &
        'cached and independently prepared response agree')
    call gs_profile_implicit_response(s, flux, ff, 2*dff, repeated, status, cache)
    call require(status == 0, 'second RHS shares same context')
    call require(maxval(abs(repeated - 2*response)) < 1e-12_dp, &
        'independent RHS homogeneity')
    call gs_profile_implicit_response(s, flux, ff, dff, repeated, status, empty)
    call require(status /= 0, 'unprepared context fails closed')
    changed = s; changed%psi = 1.1_dp*s%psi
    call gs_profile_implicit_response(changed, flux, ff, dff, repeated, status, cache)
    call require(status /= 0, 'changed flux state invalidates context')
    changed = s; changed%points(1, 1) = changed%points(1, 1) + .1_dp
    call gs_profile_implicit_response(changed, flux, ff, dff, repeated, status, cache)
    call require(status /= 0, 'changed geometry invalidates context')
    changed = s; changed%gradient(1, 1, 1) = 2
    call gs_profile_implicit_response(changed, flux, ff, dff, repeated, status, cache)
    call require(status /= 0, 'changed reference map invalidates context')
    changed = s; changed%area(1) = 2
    call gs_profile_implicit_response(changed, flux, ff, dff, repeated, status, cache)
    call require(status /= 0, 'changed integration area invalidates context')
    changed = s; changed%triangles(2:3, 1) = [3, 2]
    call gs_profile_implicit_response(changed, flux, ff, dff, repeated, status, cache)
    call require(status /= 0, 'changed connectivity invalidates context')
    changed = s; changed%element_order = 1
    call gs_profile_implicit_response(changed, flux, ff, dff, repeated, status, cache)
    call require(status /= 0, 'changed element order invalidates context')
    changed = s; deallocate (changed%psi)
    call gs_profile_implicit_response(changed, flux, ff, dff, repeated, status, cache)
    call require(status /= 0, 'incomplete numerical state fails closed')
    changed = s; changed%stiffness(1, 1, 1) = 2
    call gs_profile_implicit_response(changed, flux, ff, dff, repeated, status, cache)
    call require(status /= 0, 'changed weak operator invalidates context')
    changed = s; changed%diagonal(1) = 2
    call gs_profile_implicit_response(changed, flux, ff, dff, repeated, status, cache)
    call require(status /= 0, 'changed preconditioner invalidates context')
    changed = s; changed%boundary(1) = .true.
    call gs_profile_implicit_response(changed, flux, ff, dff, repeated, status, cache)
    call require(status /= 0, 'changed boundary domain invalidates context')
    call gs_profile_implicit_response(s, flux, ff + .01_dp, dff, &
        repeated, status, cache)
    call require(status /= 0, 'changed FF law invalidates context')
    call gs_profile_implicit_response(s, 1.1_dp*flux, ff, dff, repeated, status, cache)
    call require(status /= 0, 'changed profile knots invalidate context')
    changed = s; changed%psi = 1.1_dp*s%psi
    call prepare_gs_profile_response_cache(changed, 1.1_dp*flux, ff, renewed, status)
    call require(status == 0, 'new state prepares a distinct valid context')
    call gs_profile_implicit_response(changed, 1.1_dp*flux, ff, dff, &
        repeated, status, renewed)
    call require(status == 0, 'new state reuses its own context')
    derivative = changed; derivative%psi = repeated
    call gs_element_jet(derivative, 1, axis(1), axis(2), daxis, grad, hess, status)
    call require(status == 0, 'new-state physical response jet')
    residual = dot_product(s%points(1, :), repeated) - &
        area/1.1_dp*(.25_dp*daxis + .005_dp*(hess(1, 1) + hess(2, 2))/8) - .5_dp
    call require(abs(residual) < 1e-9_dp, 'new-state exact disk moment equation')
    call prepare_gs_profile_response_cache(s, flux, ff(:2), renewed, status)
    call require(status /= 0, 'failed preparation cannot retain old ready state')
    call gs_profile_implicit_response(s, flux, ff, dff, repeated, status, renewed)
    call require(status /= 0, 'failed preparation leaves unavailable context')
    print *, 'PASS: independent response-context checks ', passed
contains
    subroutine require(condition, message)
        logical, intent(in)::condition
        character(*), intent(in)::message
        if (.not. condition) then
            print *, message
            error stop 2
        end if
        passed = passed + 1
    end subroutine require
end program test_gs_inverse_cache_reuse
