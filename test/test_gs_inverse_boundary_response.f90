program test_gs_inverse_boundary_response
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_grad_shafranov, only: gs_state_t, gs_profile_response, &
        gs_profile_implicit_response, gs_profile_response_cache_t, &
        prepare_gs_profile_response_cache, gs_profile_load
    implicit none
    type(gs_state_t) :: state
    type(gs_profile_response_cache_t) :: cache
    real(dp), allocatable :: response(:), uncached(:), twice(:), combined(:)
    real(dp), allocatable :: primal(:), positive(:), negative(:)
    real(dp) :: flux(3), ff(3), first(3), second(3), zero(3), rhs(6), base(6)
    real(dp) :: expected(6), axis(2), residual, h
    integer :: i, status, passed
    passed = 0
    allocate (state%points(2, 6), state%triangles(6, 1), state%psi(6), &
        state%boundary(6), state%stiffness(6, 6, 1), state%gradient(2, 3, 1), &
        state%area(1), state%diagonal(6))
    state%points = reshape([2.0_dp, 0.0_dp, 3.0_dp, 0.0_dp, 2.0_dp, 1.0_dp, &
        2.5_dp, 0.0_dp, 2.5_dp, 0.5_dp, 2.0_dp, 0.5_dp], [2, 6])
    state%triangles(:, 1) = [1, 2, 3, 4, 5, 6]; state%element_order = 2
    state%gradient(:, :, 1) = reshape([-1.0_dp, -1.0_dp, 1.0_dp, 0.0_dp, &
        0.0_dp, 1.0_dp], [2, 3])
    state%boundary = .false.; state%boundary(1) = .true.
    state%area = .5_dp; state%diagonal = 1; state%stiffness = 0
    do i = 1, 6
        state%stiffness(i, i, 1) = 1
    end do
    state%stiffness(2, 1, 1) = .1_dp; state%stiffness(1, 2, 1) = .1_dp
    axis = [7.0_dp/3, 1.0_dp/3]
    state%psi = (state%points(1, :) - axis(1))**2 + &
        (state%points(2, :) - axis(2))**2 - .01_dp
    flux = [-.01_dp, -.005_dp, 1.0_dp]; ff = 0; zero = 0
    ! The eliminated free operator is identity. Primal nonzero Dirichlet data
    ! contributes -.1*psi_boundary at free node2, and has zero coefficient JVP.
    expected = 0; expected(1) = state%psi(1); expected(2) = -.1_dp*state%psi(1)
    call gs_profile_response(state, flux, ff, primal, residual, status)
    call require(status == 0, 'primal Picard map preserves nonzero Dirichlet data')
    call close (maxval(abs(primal - expected)), 1e-12_dp, 'exact affine primal map')
    call prepare_gs_profile_response_cache(state, flux, ff, cache, status)
    call require(status == 0, 'nonzero fixed Dirichlet response context is supported')
    call gs_profile_implicit_response(state, flux, ff, zero, response, status, cache)
    call require(status == 0, 'cached zero coefficient direction')
    call close (norm2(response), 1e-12_dp, 'zero coefficient JVP has zero response')
    call gs_profile_implicit_response(state, flux, ff, zero, uncached, status)
    call require(status == 0, 'independently prepared zero coefficient direction')
    call close (norm2(uncached), 1e-12_dp, 'uncached zero coefficient JVP')
    first = [.2_dp, -.3_dp, .7_dp]; second = [-.1_dp, .4_dp, .2_dp]
    call gs_profile_implicit_response(state, flux, ff, first, response, status, cache)
    call require(status == 0, 'first nonzero coefficient response')
    call gs_profile_load(state, 1.0_dp, flux, zero, first, rhs, status)
    call require(status == 0, 'independent primal source load')
    call gs_profile_load(state, 1.0_dp, flux, zero, zero, base, status)
    call require(status == 0, 'independent fixed boundary affine offset')
    call close (maxval(abs(response - (rhs - base))), 1e-11_dp, &
        'exact identity free-operator coefficient response')
    call gs_profile_implicit_response(state, flux, ff, 2*first, twice, status, cache)
    call require(status == 0, 'second scaled coefficient response')
    call close (maxval(abs(twice - 2*response)), 1e-11_dp, 'response homogeneity')
    call gs_profile_implicit_response(state, flux, ff, second, twice, status, cache)
    call require(status == 0, 'independent second direction')
    call gs_profile_implicit_response(state, flux, ff, first + second, &
        combined, status, cache)
    call require(status == 0, 'combined coefficient direction')
    call close (maxval(abs(combined - response - twice)), 1e-11_dp, 'response additivity')
    h = 1e-3_dp
    call gs_profile_response(state, flux, h*first, positive, residual, status)
    call require(status == 0, 'positive independently solved affine primal map')
    call gs_profile_response(state, flux, -h*first, negative, residual, status)
    call require(status == 0, 'negative independently solved affine primal map')
    call close (maxval(abs((positive - negative)/(2*h) - response)), &
        1e-10_dp, 'nonzero Dirichlet primal coefficient finite difference')
    print *, 'PASS: independent fixed-Dirichlet response checks ', passed
contains
    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(*), intent(in) :: message
        if (.not. condition) then
            print *, message
            error stop 2
        end if
        passed = passed + 1
    end subroutine require
    subroutine close (error, tolerance, message)
        real(dp), intent(in) :: error, tolerance
        character(*), intent(in) :: message
        call require(error <= tolerance, message)
    end subroutine close
end program test_gs_inverse_boundary_response
