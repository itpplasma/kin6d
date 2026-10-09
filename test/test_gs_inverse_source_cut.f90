program test_gs_inverse_source_cut
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use kin6d_grad_shafranov, only: gs_state_t, gs_profile_tangent_cache, &
        gs_profile_load_jvp
    implicit none
    type(gs_state_t) :: s
    real(dp), allocatable :: mass(:, :, :), axis_response(:)
    real(dp) :: flux(3), ff(3), pp(3), zero(3), direction(6), load_jvp(6)
    real(dp) :: axis_flux, eps, lhs, original(6)
    integer :: status, passed
    passed = 0
    allocate (s%points(2, 6), s%triangles(6, 1), s%psi(6), s%boundary(6), &
        s%stiffness(6, 6, 1), s%gradient(2, 3, 1), s%area(1))
    s%points = reshape([2.0_dp, 0.0_dp, 3.0_dp, 0.0_dp, 2.0_dp, 1.0_dp, &
        2.5_dp, 0.0_dp, 2.5_dp, 0.5_dp, 2.0_dp, 0.5_dp], [2, 6])
    s%triangles(:, 1) = [1, 2, 3, 4, 5, 6]; s%element_order = 2
    s%area = .5_dp; s%boundary = .false.; s%stiffness = 0
    s%psi = s%points(2, :)**2; original = s%psi
    flux = [0.0_dp, .25_dp, 1.0_dp]; ff = [0.0_dp, 0.0_dp, .75_dp]; pp = 0; zero = 0
    ! This tests the algebraic axis dependence independently of finding an axis.
    axis_flux = -.25_dp
    call gs_profile_tangent_cache(s, flux, ff, axis_flux, mass, axis_response, status)
    call require(status == 0, 'cut-conforming tangent cache status')
    call near(maxval(abs(mass(:, :, 1) - transpose(mass(:, :, 1)))), 0.0_dp, &
        1.0e-12_dp, 'symmetric weak source mass')
    direction = 1
    lhs = dot_product(s%points(1, :), matmul(mass(:, :, 1), direction))
    call near(lhs, 1.0_dp/8, 1.0e-11_dp, 'exact step area .5<eta<1')
    lhs = dot_product(s%points(1, :), matmul(mass(:, :, 1), s%points(2, :)))
    call near(lhs, 1.0_dp/12, 1.0e-11_dp, 'exact first step moment')
    lhs = dot_product(s%points(1, :), matmul(mass(:, :, 1), s%psi))
    call near(lhs, 11.0_dp/192, 1.0e-11_dp, 'exact quadratic step moment')
    lhs = dot_product(s%points(1, :), axis_response)
    call near(lhs, -11.0_dp/48, 1.0e-11_dp, 'signed axis-response moment')
    call near(maxval(abs(matmul(mass(:, :, 1), s%psi) - &
        axis_flux*axis_response)), 0.0_dp, 1.0e-10_dp, &
        'exact normalized-flux homogeneity')
    direction = [.13_dp, -.2_dp, .4_dp, .1_dp, -.3_dp, .6_dp]
    call gs_profile_load_jvp(s, 1.0_dp, flux, pp, ff, direction, 0.0_dp, &
        zero, zero, zero, load_jvp, status)
    call require(status == 0, 'independently assembled psi-source tangent')
    call near(maxval(abs(matmul(mass(:, :, 1), direction) - load_jvp)), &
        0.0_dp, 1.0e-10_dp, 'cached mass matches all-active load derivative')
    direction = 0
    call gs_profile_load_jvp(s, 1.0_dp, flux, pp, ff, direction, 0.0_dp, &
        flux/axis_flux, zero, zero, load_jvp, status)
    call require(status == 0, 'independently assembled moving profile knots')
    call near(maxval(abs(axis_response + load_jvp)), 0.0_dp, 1.0e-10_dp, &
        'axis rank-one sign and knot dependence')
    ! Reversing every FFprime value reverses both weak derivative products.
    call gs_profile_tangent_cache(s, flux, -ff, axis_flux, mass, axis_response, status)
    call require(status == 0, 'signed negative source derivative')
    call near(dot_product(s%points(1, :), axis_response), 11.0_dp/48, &
        1.0e-11_dp, 'negative derivative axis moment')
    call gs_profile_tangent_cache(s, flux, ff, axis_flux, mass, axis_response, status, &
        max_evaluations=1)
    call require(status /= 0, 'bounded cache integration fails closed')
    call near(maxval(abs(mass)), 0.0_dp, 0.0_dp, 'failed cache matrix is zero')
    call near(maxval(abs(axis_response)), 0.0_dp, 0.0_dp, 'failed axis vector is zero')
    call gs_profile_tangent_cache(s, flux, ff, 0.0_dp, mass, axis_response, status)
    call require(status /= 0, 'undefined zero axis denominator rejected')
    eps = ieee_value(1.0_dp, ieee_quiet_nan)
    call gs_profile_tangent_cache(s, flux, ff, eps, mass, axis_response, status)
    call require(status /= 0, 'nonfinite axis denominator rejected')
    s%psi = .25_dp
    call gs_profile_tangent_cache(s, flux, ff, axis_flux, mass, axis_response, status)
    call require(status == 4, 'whole-cell kink derivative unavailable')
    s%psi = original; s%boundary(2) = .true.
    call gs_profile_tangent_cache(s, flux, ff, axis_flux, mass, axis_response, status)
    call require(status == 0, 'fixed Dirichlet axis cache')
    call near(axis_response(2), 0.0_dp, 0.0_dp, 'boundary axis response excluded')
    print *, 'PASS: independent cut-cache checks ', passed
contains
    subroutine require(condition, message)
        logical, intent(in)::condition
        character(*), intent(in)::message
        if (.not. condition) then
            print *, message
            error stop 1
        end if
        passed = passed + 1
    end subroutine require
    subroutine near(value, expected, tolerance, message)
        real(dp), intent(in)::value, expected, tolerance
        character(*), intent(in)::message
        if (abs(value - expected) > tolerance) then
            print *, message, value, expected
            error stop 2
        end if
        passed = passed + 1
    end subroutine near
end program test_gs_inverse_source_cut
