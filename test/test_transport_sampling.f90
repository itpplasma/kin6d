program test_transport_sampling
    use, intrinsic :: iso_fortran_env, only: dp => real64, qp => real128
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use fortnum_interval, only: interval_t, interval
    use kin6d_transport_laboratory, only: transport_state_t, initialize_transport, &
        solve_transport, sample_transport_coefficients_interval
    use kin6d_transport_certificate, only: transport_certificate_t, &
        certify_transport_tangent, certify_transport_sample, &
        certify_transport_tangent_sample
    use kin6d_transport_sampling, only: transport_basis_t, &
        transport_sampling_workspace_t, initialize_transport_sampling, &
        prepare_transport_basis, sample_transport_cached, &
        sample_transport_input_cached, &
        sample_transport_coefficients_cached, certify_transport_cached, &
        certify_transport_tangent_cached, transport_sampling_bytes
    implicit none
    integer, parameter :: nx = 3, nodes = 5, modes = 2
    type(transport_state_t) :: state, empty
    type(transport_certificate_t) :: certificate
    type(transport_basis_t) :: basis, unprepared, zero_basis
    type(transport_sampling_workspace_t) :: work, empty_work, zero_work
    type(interval_t) :: geometry(nx), bad_geometry(nx), range(nx), direct, velocity
    real(dp) :: dcos(0:modes, nodes), dsin(0:modes, nodes)
    real(dp) :: dsource_cos(0:modes, nodes), dsource_sin(0:modes, nodes)
    real(dp) :: lower(nx), upper(nx), dlower(nx), dupper(nx), lo, hi
    real(dp) :: x, v, nan, saved
    real(qp) :: candidate, exact, derivative
    logical :: ok, strict
    integer :: j, i, r, evaluator

    nan = ieee_value(0.0_dp, ieee_quiet_nan)
    geometry(1) = interval(0.17_dp, 0.33_dp)
    geometry(2) = interval(1.21_dp, 1.43_dp)
    geometry(3) = interval(4.65_dp, 4.94_dp)
    call initialize_transport(state, 4, modes, 1.0_dp, ok)
    if (.not. ok) error stop 'sampling fixture initialization'
    ! Arbitrary unsolved coefficients, unrelated to nodal inversion.
    state%cosine(0, :) = [1.0_dp, 0.25_dp, 1.5_dp, 0.75_dp, 1.25_dp]
    state%cosine(1, :) = [0.25_dp, -0.5_dp, 0.75_dp, 0.5_dp, -0.25_dp]
    state%cosine(2, :) = [-0.125_dp, 0.375_dp, 0.25_dp, -0.5_dp, 0.125_dp]
    state%sine(1, :) = [-0.25_dp, 0.5_dp, 0.125_dp, -0.75_dp, 0.375_dp]
    state%sine(2, :) = [0.5_dp, -0.25_dp, 0.75_dp, 0.125_dp, -0.5_dp]
    state%source_cos = state%cosine
    state%source_sin = state%sine
    call initialize_transport_sampling(basis, work, modes, nx, ok)
    if (.not. ok) error stop 'basis initialization'
    if (transport_sampling_bytes(basis, work) /= 432) &
        error stop 'owned interval buffer accounting'
    call prepare_transport_basis(basis, geometry, ok)
    if (.not. ok) error stop 'basis preparation'
    velocity = interval(-0.77_dp, 0.66_dp)
    call sample_transport_cached(state, basis, geometry, velocity, work, range, ok)
    if (.not. ok) error stop 'cross-velocity-cell cached reconstruction'
    do i = 1, nx
        call sample_transport_coefficients_interval(state, state%cosine, state%sine, &
            geometry(i), velocity, direct, ok)
        if (.not. ok) error stop 'scalar reconstruction fallback'
        ! Both interval schedules must contain the same independent interior values.
        do r = 1, 11
            x = geometry(i)%lo + real(r, dp)*(geometry(i)%hi - geometry(i)%lo)/12
            v = velocity%lo + real(r, dp)*(velocity%hi - velocity%lo)/12
            candidate = independent_candidate(state%v, state%cosine, state%sine, &
                real(x, qp), real(v, qp))
            call contains(range(i), candidate)
            call contains(direct, candidate)
        end do
        ! The split schedule uses identical rigorous factors and arithmetic order;
        ! skipping neutral mode-zero terms can only remove rounding inflation.
        if (range(i)%lo < direct%lo .or. range(i)%hi > direct%hi) &
            error stop 'cached range unexpectedly widens scalar enclosure'
    end do
    state%source_cos(0, :) = -0.75_dp
    call sample_transport_input_cached(state, basis, geometry, velocity, work, &
        range, ok)
    if (.not. ok) error stop 'cached supplied-source reconstruction'
    do i = 1, nx
        candidate = independent_candidate(state%v, state%source_cos, &
            state%source_sin, real(geometry(i)%hi, qp), 0.17_qp)
        call contains(range(i), candidate)
    end do
    do j = 1, nodes
        call sample_transport_cached(state, basis, geometry, interval(state%v(j)), &
            work, range, ok)
        if (.not. ok) error stop 'exact velocity node or domain boundary rejected'
        do i = 1, nx
            candidate = independent_candidate(state%v, state%cosine, state%sine, &
                real(geometry(i)%lo, qp), real(state%v(j), qp))
            call contains(range(i), candidate)
        end do
    end do

    ! Independently supplied affine source; continuous resolvent is TEST ONLY.
    state%source_cos(0, :) = [0.875_dp, 0.9375_dp, 1.0_dp, 1.0625_dp, 1.125_dp]
    state%source_cos(1, :) = [0.5_dp, 0.375_dp, 0.25_dp, 0.125_dp, 0.0_dp]
    state%source_cos(2, :) = 0.0_dp
    state%source_sin = 0.0_dp
    state%source_sin(2, :) = [-0.125_dp, -0.0625_dp, 0.0_dp, 0.0625_dp, 0.125_dp]
    call solve_transport(state, ok)
    if (.not. ok) error stop 'solved certificate fixture'
    ! An arbitrary tangent candidate is also bounded by its computed residual.
    dcos = state%cosine
    dsin = state%sine
    dsource_cos = 0.0_dp; dsource_sin = 0.0_dp
    do evaluator = 1, 2
        strict = evaluator == 2
        call certify_transport_tangent(state, 1.0_dp, dsource_cos, dsource_sin, &
            dcos, dsin, certificate, ok, strict=strict)
        if (.not. ok) error stop 'certificate preparation for sampling'
        call certify_transport_cached(state, certificate, basis, geometry, &
            velocity, work, lower, upper, ok)
        if (.not. ok) error stop 'cached continuous state bands'
        call certify_transport_tangent_cached(state, certificate, dcos, dsin, &
            basis, geometry, velocity, work, dlower, dupper, ok)
        if (.not. ok) error stop 'cached continuous tangent bands'
        do i = 1, nx
            call certify_transport_sample(state, certificate, geometry(i), &
                velocity, lo, hi, ok)
            if (.not. ok) error stop 'direct continuous state bands'
            if (lower(i) < lo .or. upper(i) > hi) error stop 'state band schedule'
            call certify_transport_tangent_sample(state, certificate, dcos, dsin, &
                geometry(i), velocity, lo, hi, ok)
            if (.not. ok) error stop 'direct continuous tangent bands'
            if (dlower(i) < lo .or. dupper(i) > hi) error stop 'tangent band schedule'
            do r = 1, 9
                x = geometry(i)%lo + real(r, dp)*(geometry(i)%hi - geometry(i)%lo)/10
                v = velocity%lo + real(r, dp)*(velocity%hi - velocity%lo)/10
                call independent_resolvent(real(x, qp), real(v, qp), exact, derivative)
                call contains(interval(lower(i), upper(i)), exact)
                call contains(interval(dlower(i), dupper(i)), derivative)
            end do
        end do
    end do

    bad_geometry = geometry
    bad_geometry(2)%hi = 1.44_dp
    call sample_transport_cached(state, basis, bad_geometry, velocity, work, range, ok)
    if (ok) error stop 'stale captured geometry accepted'
    bad_geometry = geometry
    bad_geometry(1)%lo = nan
    call prepare_transport_basis(basis, bad_geometry, ok)
    if (ok) error stop 'NaN geometry accepted'
    call sample_transport_cached(state, basis, geometry, velocity, work, range, ok)
    if (ok) error stop 'failed preparation retained a valid cache'
    call prepare_transport_basis(basis, geometry, ok)
    if (.not. ok) error stop 'valid cache reprepare'
    bad_geometry = geometry
    bad_geometry(1)%lo = 1.0_dp
    bad_geometry(1)%hi = 0.0_dp
    call prepare_transport_basis(basis, bad_geometry, ok)
    if (ok) error stop 'reversed interval geometry accepted'
    call prepare_transport_basis(basis, geometry, ok)
    call sample_transport_cached(state, unprepared, geometry, velocity, work, range, ok)
    if (ok) error stop 'uninitialized basis accepted'
    call sample_transport_cached(state, basis, geometry, velocity, &
        empty_work, range, ok)
    if (ok) error stop 'uninitialized workspace accepted'
    call sample_transport_cached(empty, basis, geometry, velocity, work, range, ok)
    if (ok) error stop 'uninitialized state accepted'
    call sample_transport_cached(state, basis, geometry, velocity, work, range(:2), ok)
    if (ok) error stop 'wrong result shape accepted'
    call sample_transport_cached(state, basis, geometry(:2), velocity, work, range, ok)
    if (ok) error stop 'wrong geometry shape accepted'
    call sample_transport_coefficients_cached(state, state%cosine(:1, :), &
        state%sine(:1, :), basis, geometry, velocity, work, range, ok)
    if (ok) error stop 'wrong mode shape accepted'
    call sample_transport_coefficients_cached(state, state%cosine(:, :4), &
        state%sine(:, :4), basis, geometry, velocity, work, range, ok)
    if (ok) error stop 'wrong node shape accepted'
    call certify_transport_tangent_cached(state, certificate, dcos(:1, :), &
        dsin(:1, :), basis, geometry, velocity, work, lower, upper, ok)
    if (ok) error stop 'wrong certified tangent shape accepted'
    call sample_transport_cached(state, basis, geometry, interval(nan), work, range, ok)
    if (ok) error stop 'NaN velocity accepted'
    call sample_transport_cached(state, basis, geometry, interval(-1.01_dp), &
        work, range, ok)
    if (ok) error stop 'outside velocity domain accepted'
    saved = state%cosine(1, 2)
    state%cosine(1, 2) = nan
    call sample_transport_cached(state, basis, geometry, velocity, work, range, ok)
    if (ok) error stop 'NaN touched coefficient accepted'
    state%cosine(1, 2) = saved
    state%nu = 2.0_dp
    call certify_transport_cached(state, certificate, basis, geometry, velocity, &
        work, lower, upper, ok)
    if (ok) error stop 'stale certificate metadata accepted'
    state%nu = 1.0_dp
    ! Exact mode-zero scheduling has no interval trig call even for huge x.
    call initialize_transport_sampling(zero_basis, zero_work, 0, 1, ok)
    call prepare_transport_basis(zero_basis, [interval(huge(1.0_dp))], ok)
    if (.not. ok) error stop 'mode-zero-only huge finite geometry'
    call sample_transport_coefficients_cached(state, state%cosine(:0, :), &
        state%sine(:0, :), zero_basis, [interval(huge(1.0_dp))], interval(0.0_dp), &
        zero_work, range(:1), ok)
    if (.not. ok) error stop 'mode-zero-only reconstruction'
    call contains(range(1), 1.0_qp)
    print '(a)', 'cached transport interior, boundary and guard oracles PASS'
contains
    subroutine contains(bounds, value)
        type(interval_t), intent(in) :: bounds
        real(qp), intent(in) :: value
        if (value < real(bounds%lo, qp) .or. value > real(bounds%hi, qp)) &
            error stop 'interval misses independent binary128 value'
    end subroutine

    function independent_candidate(grid, a, b, xq, vq) result(value)
        ! Deliberately independent TEST ONLY binary128 P1/Fourier evaluation.
        real(dp), intent(in) :: grid(:), a(0:, :), b(0:, :)
        real(qp), intent(in) :: xq, vq
        real(qp) :: value, t, aq, bq
        integer :: cell, k

        cell = 1
        do while (cell < size(grid) - 1)
            if (vq < real(grid(cell + 1), qp)) exit
            cell = cell + 1
        end do
        t = (vq - real(grid(cell), qp))/ &
            (real(grid(cell + 1), qp) - real(grid(cell), qp))
        value = 0.0_qp
        do k = 0, size(a, 1) - 1
            aq = (1 - t)*real(a(k, cell), qp) + t*real(a(k, cell + 1), qp)
            bq = (1 - t)*real(b(k, cell), qp) + t*real(b(k, cell + 1), qp)
            value = value + aq*cos(real(k, qp)*xq) + bq*sin(real(k, qp)*xq)
        end do
    end function

    subroutine independent_resolvent(xq, vq, value, derivative)
        real(qp), intent(in) :: xq, vq
        real(qp), intent(out) :: value, derivative
        complex(qp) :: r1, r2

        r1 = cmplx((1 - vq)/4, 0.0_qp, qp)/cmplx(1.0_qp, vq, qp)
        r2 = cmplx(0.0_qp, -vq/8, qp)/cmplx(1.0_qp, 2*vq, qp)
        value = 1 + vq/8 + real(r1*exp(cmplx(0.0_qp, xq, qp)), qp) + &
            real(r2*exp(cmplx(0.0_qp, 2*xq, qp)), qp)
        derivative = -1 - vq/8 - &
            real(r1/cmplx(1.0_qp, vq, qp)*exp(cmplx(0.0_qp, xq, qp)), qp) - &
            real(r2/cmplx(1.0_qp, 2*vq, qp)*exp(cmplx(0.0_qp, 2*xq, qp)), qp)
    end subroutine
end program test_transport_sampling
