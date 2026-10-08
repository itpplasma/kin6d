program test_gs_inverse_rhs_batch
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use kin6d_grad_shafranov, only: gs_state_t, gs_profile_ff_load_basis, &
        gs_profile_load, gs_profile_load_jvp
    implicit none
    type(gs_state_t) :: state, affine
    real(dp), allocatable :: basis(:, :), errors(:, :)
    real(dp), allocatable :: refined(:, :), refined_error(:, :)
    real(dp) :: flux(3), pp(3), ff(3), load(6), tangent(6), expected(3)
    real(dp) :: direction(3), zero_state(6), zero_profile(3), shifted(3), exact(6)
    integer :: status, column, passed
    passed = 0
    allocate (state%points(2, 6), state%triangles(6, 1), state%psi(6), &
        state%boundary(6), state%stiffness(6, 6, 1), state%area(1))
    state%points = reshape([2.0_dp, 0.0_dp, 3.0_dp, 0.0_dp, 2.0_dp, 1.0_dp, &
        2.5_dp, 0.0_dp, 2.5_dp, 0.5_dp, 2.0_dp, 0.5_dp], [2, 6])
    state%triangles(:, 1) = [1, 2, 3, 4, 5, 6]; state%element_order = 2
    state%area = .5_dp; state%boundary = .false.; state%stiffness = 0
    state%psi = state%points(2, :)**2
    flux = [0.0_dp, .25_dp, 1.0_dp]; pp = 0; zero_state = 0; zero_profile = 0
    call gs_profile_ff_load_basis(state, flux, basis, errors, status)
    call require(status == 0, 'quadratic native-cut profile basis')
    ! Independent integrals int_0^1 (1-eta)*hat_j(eta^2) deta.
    ! The represented affine test function R cancels the physical source 1/R.
    expected = [13.0_dp/48, 7.0_dp/36, 5.0_dp/144]
    do column = 1, 3
        call close (dot_product(state%points(1, :), basis(:, column)), &
            expected(column), 1e-11_dp, 'exact native-cut hat moment')
        ff = 0; ff(column) = 1
        call gs_profile_load(state, 1.0_dp, flux, pp, ff, load, status)
        call require(status == 0, 'independent complete column load')
        call close (maxval(abs(load - basis(:, column))), 0.0_dp, 1e-11_dp, &
            'each batch column against original native API')
        call require(norm2(errors(:, column)) <= 1e-11_dp*norm2(basis(:, column)), &
            'each column retains native global estimator budget')
    end do
    call gs_profile_ff_load_basis(state, flux, refined, refined_error, status, &
        relative_tolerance=1e-13_dp)
    call require(status == 0, 'independent tighter vector quadrature')
    call close (maxval(abs(refined - basis)), 0.0_dp, 1e-11_dp, 'column refinement')
    direction = [-.3_dp, .7_dp, -.2_dp]
    call gs_profile_load(state, 1.0_dp, flux, pp, direction, load, status)
    call require(status == 0, 'independent signed combined load')
    call close (maxval(abs(matmul(basis, direction) - load)), 0.0_dp, 1e-11_dp, &
        'profile coefficient linearity')
    ff = [.1_dp, .2_dp, .3_dp]
    call gs_profile_load_jvp(state, 1.0_dp, flux, pp, ff, zero_state, 0.0_dp, &
        zero_profile, zero_profile, direction, tangent, status)
    call require(status == 0, 'independent full coefficient JVP')
    call close (maxval(abs(tangent - matmul(basis, direction))), 0.0_dp, 1e-11_dp, &
        'batch contraction agrees with original coefficient JVP')
    shifted = flux + 1048576.0_dp; state%psi = state%psi + 1048576.0_dp
    call gs_profile_ff_load_basis(state, shifted, refined, refined_error, status)
    call require(status == 0, 'large simultaneous exactly representable gauge')
    call close (maxval(abs(basis - refined)), 0.0_dp, 1e-11_dp, 'gauge-stable columns')
    ! New state: FFprime=psi=R gives source1 exactly, so every P2 row has
    ! an independent polynomial moment: zero vertices and 1/6 at each edge.
    state%psi = state%points(1, :); flux = [2.0_dp, 2.5_dp, 3.0_dp]
    call gs_profile_ff_load_basis(state, flux, refined, refined_error, status)
    call require(status == 0, 'new state computes new basis')
    exact = [0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp/6, 1.0_dp/6, 1.0_dp/6]
    call close (maxval(abs(matmul(refined, flux) - exact)), 0.0_dp, 1e-11_dp, &
        'six exact P2 source polynomial moments')
    do column = 1, 3
        ff = 0; ff(column) = 1
        call gs_profile_load(state, 1.0_dp, flux, pp, ff, load, status)
        call require(status == 0, 'new-state native column')
        call close (maxval(abs(refined(:, column) - load)), 0.0_dp, 1e-11_dp, &
            'each new-state column against independent native load')
    end do
    call require(maxval(abs(refined - basis)) > 1e-3_dp, &
        'different state cannot retain old numerical columns')
    call gs_profile_ff_load_basis(state, flux, basis, errors, status, max_evaluations=1)
    call require(status == 3, 'same unresolved quadrature-budget status')
    call require(all(basis == 0), 'budget failure does not expose partial loads')
    shifted = flux; shifted(2) = ieee_value(0.0_dp, ieee_quiet_nan)
    call gs_profile_ff_load_basis(state, shifted, basis, errors, status)
    call require(status == 1, 'nonfinite profile domain fails closed')
    call require(all(basis == 0), 'invalid input does not expose old columns')
    call gs_profile_ff_load_basis(state, flux(3:1:-1), basis, errors, status)
    call require(status == 1, 'unordered profile domain fails closed')
    allocate (affine%points(2, 3), affine%triangles(3, 1), affine%psi(3), &
        affine%boundary(3), affine%stiffness(3, 3, 1), affine%area(1))
    affine%points = state%points(:, :3); affine%triangles(:, 1) = [1, 2, 3]
    affine%psi = affine%points(1, :); affine%area = .5_dp
    affine%boundary = .false.; affine%stiffness = 0; affine%element_order = 1
    call gs_profile_ff_load_basis(affine, flux, basis, errors, status)
    call require(status == 0, 'P1 vector load basis')
    call close (maxval(abs(matmul(basis, flux) - 1.0_dp/6)), &
        0.0_dp, 1e-11_dp, 'three exact P1 polynomial moments')
    do column = 1, 3
        ff = 0; ff(column) = 1
        call gs_profile_load(affine, 1.0_dp, flux, pp, ff, load(:3), status)
        call require(status == 0, 'P1 original column API')
        call close (maxval(abs(basis(:, column) - load(:3))), 0.0_dp, &
            1e-11_dp, 'each P1 column against original native load')
    end do
    affine%boundary(1) = .true.
    call gs_profile_ff_load_basis(affine, flux, basis, errors, status)
    call require(status == 0, 'Dirichlet response-source basis')
    call require(all(basis(1, :) == 0), 'Dirichlet source rows are excluded')
    affine%points(1, 1) = 0
    call gs_profile_ff_load_basis(affine, flux, basis, errors, status)
    call require(status == 1, 'nonpositive physical R fails closed')
    call require(all(basis == 0), 'failed geometry clears partial columns')
    print *, 'PASS: independent batched profile-source checks ', passed
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
    subroutine close (actual, wanted, tolerance, message)
        real(dp), intent(in) :: actual, wanted, tolerance
        character(*), intent(in) :: message
        call require(abs(actual - wanted) <= tolerance, message)
    end subroutine close
end program test_gs_inverse_rhs_batch
