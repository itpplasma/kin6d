program test_tiago_island
    use, intrinsic :: iso_fortran_env, only: dp => real64, int32
    use island_a_fixture, only: island_a_setup, shell_source, island_a_sample_points
    use nn_vacuum_fixture, only: field, coils, current, ftor
    use milca_coils_3d, only: coil_set_b_cart
    use milca_nonnested_coupled, only: nc_problem_t
    use milca_island_provider, only: island_provider_t, island_provider_init, &
        island_provider_solve, island_provider_sample, island_provider_vjp, &
        island_provider_free, island_provider_provenance, IP_OK, IP_NTHETA, IP_EPSILON
    use milca_tiago_island_mse, only: island_mse_posterior_t, island_mse_predict
    use tiago_diagnostic_types, only: mse_channel_t
    use tiago_mse, only: mse_angles
    use tiago_inference_engines, only: map_fit
    implicit none
    type(island_mse_posterior_t) :: model
    type(island_provider_t) :: reference, cold
    type(mse_channel_t), allocatable :: channels(:)
    real(dp), allocatable :: xyz(:, :), truth(:), prediction(:), derivative(:)
    real(dp), allocatable :: plus(:), minus(:), fd(:), grad(:, :), angles(:)
    real(dp), allocatable :: weights(:, :), values(:, :), angle_weights(:)
    integer, allocatable :: train(:), held(:)
    real(dp), parameter :: h = 1.0e-4_dp
    real(dp), allocatable :: epsilon(:)
    real(dp) :: source_scale, theta(1), theta_bar(IP_NTHETA), relative, lhs, rhs
    real(dp) :: radius, phi, bcart(3), parameter_tolerance, held_tolerance
    integer :: i, e, steps, status, n, nr, order, modes
    integer(int32) :: ierr
    character(16) :: case_name
    character(32) :: argument
    character(:), allocatable :: message
    logical :: ok, converged

    call get_command_argument(1, case_name)
    nr = 4; order = 2; modes = 1
    epsilon = [0.0_dp, 0.5_dp]
    parameter_tolerance = 1.0e-3_dp; held_tolerance = 1.0e-2_dp
    select case (trim(case_name))
    case ('vacuum')
        source_scale = 0
    case ('pressure')
        source_scale = 1
    case ('direct')
        source_scale = 0
        epsilon = [1.0_dp]
        parameter_tolerance = 1.0e-2_dp; held_tolerance = 1.0_dp
    case default
        error stop 'expected vacuum, pressure or direct'
    end select
    if (command_argument_count() >= 2) then
        call get_command_argument(2, argument); read(argument, *) nr
    end if
    if (command_argument_count() >= 3) then
        call get_command_argument(3, argument); read(argument, *) order
    end if
    if (command_argument_count() >= 4) then
        call get_command_argument(4, argument); read(argument, *) modes
    end if
    if (nr < 1 .or. order < 1 .or. modes < 1) error stop 'invalid resolution'
    print '(a,3i6)', 'radial subdivisions / element order / Fourier modes ', &
        nr, order, modes
    xyz = island_a_sample_points()
    n = size(xyz, 2)
    allocate (channels(n), values(4, n), weights(4, n))
    train = pack([(i, i=1, n)], [(mod(i, 3) /= 0, i=1, n)])
    held = pack([(i, i=1, n)], [(mod(i, 3) == 0, i=1, n)])
    do i = 1, n
        channels(i)%beam_point = xyz(:, i)
        channels(i)%sight_point = xyz(:, i)
        channels(i)%beam_direction = [1.0_dp, 0.3_dp, 0.2_dp]
        channels(i)%sight_direction = [0.2_dp, 1.0_dp, -0.1_dp]
        channels(i)%reference = [0.0_dp, 0.0_dp, 1.0_dp]
        channels(i)%max_gap = 1.0e-8_dp
    end do
    model%channels = channels(train)
    angle_weights = sin(real([(i, i=1, n)], dp))

    do e = 1, size(epsilon)
        call initialize(reference, epsilon(e))
        call island_mse_predict(reference, channels, epsilon(e), truth, ok, &
            derivative)
        if (.not. ok) error stop 'truth prediction'
        print '(a)', island_provider_provenance(reference)

        ! Both differences start from a fresh state, away from a warm branch.
        call initialize(cold, epsilon(e) + h)
        call island_mse_predict(cold, channels, epsilon(e) + h, plus, ok)
        if (.not. ok) error stop 'cold plus'
        call initialize(cold, epsilon(e) - h)
        call island_mse_predict(cold, channels, epsilon(e) - h, minus, ok)
        if (.not. ok) error stop 'cold minus'
        fd = (plus - minus)/(2*h)
        relative = maxval(abs(fd - derivative))/maxval(abs(fd))
        print '(a,f6.3,a,es14.6)', 'epsilon ', epsilon(e), &
            ' cold composed derivative relative error ', relative
        if (relative > 2.0e-5_dp) error stop 'composed cold derivative'

        call island_provider_sample(reference, xyz, values, status)
        if (status /= IP_OK) error stop 'adjoint sample'
        call mse_angles(channels, values(1:3, :), angles, ierr, message, grad)
        if (ierr /= 0) error stop 'adjoint MSE geometry'
        weights = 0
        weights(1:3, :) = grad*spread(angle_weights, 1, 3)
        call island_provider_vjp(reference, xyz, weights, theta_bar, status)
        if (status /= IP_OK) error stop 'composed adjoint'
        lhs = dot_product(angle_weights, derivative)
        rhs = theta_bar(IP_EPSILON)
        if (abs(lhs - rhs) > 1.0e-9_dp*max(abs(lhs), abs(rhs))) &
            error stop 'composed adjoint pairing'

        if (trim(case_name) == 'direct') then
            ! Independent observations retain all coil harmonics and no FE solve.
            do i = 1, n
                call coil_set_b_cart(coils, [current], xyz(:, i), bcart)
                radius = norm2(xyz(:2, i)); phi = atan2(xyz(2, i), xyz(1, i))
                values(1:3, i) = bcart + &
                    [-ftor*sin(phi)/radius, ftor*cos(phi)/radius, 0.0_dp]
            end do
            call mse_angles(channels, values(1:3, :), truth, ierr, message)
            if (ierr /= 0) error stop 'direct-coil truth geometry'
        end if
        model%observed = truth(train)
        call initialize(model%provider, 0.2_dp)
        theta = 0.2_dp
        call map_fit(model, theta, 12, steps, message, converged, 1.0e-8_dp)
        print '(a)', message
        print '(a,2f12.7,a,i0)', 'truth / reconstructed epsilon ', epsilon(e), &
            theta(1), ' accepted steps ', steps
        if (.not. converged) error stop 'TIAGO MAP did not converge'
        call island_mse_predict(model%provider, channels(held), &
            theta(1), prediction, ok)
        if (.not. ok) error stop 'held prediction'
        relative = maxval(abs(prediction - truth(held)))/model%sigma
        print '(a,es14.6)', 'held angle maximum error / assigned sigma ', relative
        if (abs(theta(1) - epsilon(e)) > parameter_tolerance) &
            error stop 'parameter recovery'
        if (relative > held_tolerance) error stop 'held prediction recovery'
    end do
    call island_provider_free(reference)
    call island_provider_free(cold)
    call island_provider_free(model%provider)
    print '(a)', 'PASS: Cartesian KIN6D provider with existing TIAGO MSE and MAP'
contains
    subroutine initialize(provider, eps)
        type(island_provider_t), intent(inout) :: provider
        real(dp), intent(in) :: eps
        type(nc_problem_t) :: problem
        real(dp), allocatable :: base(:), perturbation(:)
        integer :: st
        call island_a_setup(nr, order, modes, base, perturbation, problem)
        call island_provider_init(provider, field, problem, base, perturbation, &
            shell_source, [source_scale, eps, 0.0_dp, 100.0_dp], 0.9_dp, &
            1.0e5_dp, 'synthetic-helical-coil-tiago', st)
        if (st /= IP_OK) error stop 'provider initialization'
        call island_provider_solve(provider, provider%theta, st, tol=1.0e-13_dp)
        if (st /= IP_OK) then
            print '(a)', trim(provider%message)
            print '(a,5es16.7)', 'failed solve residuals ', provider%solve_report%residual
            print '(a,4es16.7)', 'failed solve theta ', provider%theta
            error stop 'provider initial solve'
        end if
    end subroutine
end program
