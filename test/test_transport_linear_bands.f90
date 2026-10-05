program test_transport_linear_bands
    use, intrinsic :: iso_fortran_env, only: dp => real64, qp => real128
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use fortnum_interval, only: interval_t, interval
    use kin6d_transport_laboratory, only: transport_state_t, initialize_transport, &
        solve_transport, solve_transport_tangent
    use kin6d_transport_certificate, only: transport_certificate_t, &
        certify_transport_tangent
    use kin6d_transport_sampling, only: transport_basis_t, &
        transport_sampling_workspace_t, initialize_transport_sampling, &
        prepare_transport_basis, transport_linear_enclosure, &
        certify_transport_linear, certify_transport_tangent_linear
    implicit none
    integer, parameter :: modes = 3, nodes = 5, points = 4
    type(transport_state_t) :: state, empty
    type(transport_certificate_t) :: certificate
    type(transport_basis_t) :: basis, peak_basis, zero_basis, missing
    type(transport_sampling_workspace_t) :: work, peak_work, zero_work, empty_work
    type(interval_t) :: geometry(points), bad(points), peak(2), velocity
    real(dp) :: lower(points), upper(points), tangent_lo(points), tangent_hi(points)
    real(dp) :: peak_lo(2), peak_hi(2), zero_lo(2), zero_hi(2)
    real(dp) :: dsource_cos(0:modes, nodes), dsource_sin(0:modes, nodes)
    real(dp) :: dcos(0:modes, nodes), dsin(0:modes, nodes)
    real(dp) :: nan, saved
    real(qp) :: xq, vq, value, tangent, t
    integer :: i, r, j, evaluator
    logical :: ok

    nan = ieee_value(0.0_dp, ieee_quiet_nan)
    geometry = [interval(-1.0_dp), interval(-0.125_dp), &
                interval(0.5_dp), interval(2.75_dp)]
    velocity = interval(-0.77_dp, 0.66_dp)
    call initialize_transport(state, 4, modes, 1.0_dp, ok)
    if (.not. ok) error stop 'linear fixture initialization'
    call initialize_transport_sampling(basis, work, modes, points, ok)
    if (.not. ok) error stop 'linear basis initialization'
    call prepare_transport_basis(basis, geometry, ok)
    if (.not. ok) error stop 'linear basis preparation'
    ! Coarse, nonuniform x cells and arbitrary unsolved nodal coefficients.
    state%cosine(0, :) = [1.0_dp, 0.25_dp, 1.5_dp, 0.75_dp, 1.25_dp]
    state%cosine(1, :) = [0.25_dp, -0.5_dp, 0.75_dp, 0.5_dp, -0.25_dp]
    state%cosine(2, :) = [-0.125_dp, 0.375_dp, 0.25_dp, -0.5_dp, 0.125_dp]
    state%cosine(3, :) = [0.375_dp, 0.125_dp, -0.25_dp, 0.5_dp, 0.25_dp]
    state%sine(1, :) = [-0.25_dp, 0.5_dp, 0.125_dp, -0.75_dp, 0.375_dp]
    state%sine(2, :) = [0.5_dp, -0.25_dp, 0.75_dp, 0.125_dp, -0.5_dp]
    state%sine(3, :) = [-0.5_dp, 0.25_dp, 0.125_dp, -0.125_dp, 0.375_dp]
    call transport_linear_enclosure(state, state%cosine, state%sine, basis, &
        geometry, velocity, work, 0.0_dp, lower, upper, ok)
    if (.not. ok) error stop 'arbitrary linear candidate enclosure'
    call check_candidate(state%cosine, state%sine, geometry, lower, upper, velocity)
    ! Source enclosure has radius zero; distinct coefficients detect wrong selection.
    state%source_cos = state%cosine
    state%source_sin = state%sine
    state%source_cos(0, :) = -0.75_dp
    call transport_linear_enclosure(state, state%source_cos, state%source_sin, &
        basis, geometry, velocity, work, 0.0_dp, lower, upper, ok)
    if (.not. ok) error stop 'source-only linear enclosure'
    call check_candidate(state%source_cos, state%source_sin, geometry, &
        lower, upper, velocity)

    ! Endpoint-only intervals for cos(x) miss its interior maximum at x=0.
    peak = [interval(-1.0_dp), interval(1.0_dp)]
    call initialize_transport_sampling(peak_basis, peak_work, modes, 2, ok)
    call prepare_transport_basis(peak_basis, peak, ok)
    if (.not. ok) error stop 'peak point basis'
    state%cosine = 0.0_dp; state%sine = 0.0_dp
    state%cosine(1, :) = 1.0_dp
    call transport_linear_enclosure(state, state%cosine, state%sine, peak_basis, &
        peak, interval(0.0_dp), peak_work, 0.0_dp, peak_lo, peak_hi, ok)
    if (.not. ok) error stop 'cosine peak linear enclosure'
    if (cos(-1.0_qp) >= 1.0_qp) error stop 'endpoint-only failure fixture invalid'
    call linear_contains(peak, peak_lo, peak_hi, 1, 0.5_qp, 1.0_qp)
    do r = 0, 20
        t = real(r, qp)/20
        xq = -1 + 2*t
        call linear_contains(peak, peak_lo, peak_hi, 1, t, cos(xq))
    end do
    ! Both vertices of a coarse cell border tiny cells: adjacent MAX is required.
    bad = [interval(-1.03125_dp), interval(-1.0_dp), &
           interval(1.0_dp), interval(1.03125_dp)]
    call prepare_transport_basis(basis, bad, ok)
    if (.not. ok) error stop 'adjacent curvature fixture basis'
    call transport_linear_enclosure(state, state%cosine, state%sine, basis, bad, &
        interval(0.0_dp), work, 0.0_dp, lower, upper, ok)
    if (.not. ok) error stop 'adjacent curvature fixture enclosure'
    call linear_contains(bad, lower, upper, 2, 0.5_qp, 1.0_qp)
    call prepare_transport_basis(basis, geometry, ok)
    if (.not. ok) error stop 'restore nonuniform geometry'
    ! Mode-zero curvature must vanish; large coarse spacing adds no chord width.
    call initialize_transport_sampling(zero_basis, zero_work, 0, 2, ok)
    call prepare_transport_basis(zero_basis, peak, ok)
    state%cosine(0, :) = 0.75_dp
    call transport_linear_enclosure(state, state%cosine(:0, :), state%sine(:0, :), &
        zero_basis, peak, interval(0.0_dp), zero_work, 0.0_dp, zero_lo, zero_hi, ok)
    if (.not. ok) error stop 'zero-mode linear enclosure'
    if (maxval(zero_hi - zero_lo) > 1.0e-12_dp) &
        error stop 'zero-mode curvature unexpectedly adds width'
    call linear_contains(peak, zero_lo, zero_hi, 1, 0.37_qp, 0.75_qp)

    ! Continuous resolvent oracle is independent and used only by this TEST.
    state%source_cos = 0.0_dp; state%source_sin = 0.0_dp
    state%source_cos(0, :) = [0.875_dp, 0.9375_dp, 1.0_dp, 1.0625_dp, 1.125_dp]
    state%source_cos(1, :) = [0.5_dp, 0.375_dp, 0.25_dp, 0.125_dp, 0.0_dp]
    state%source_sin(2, :) = [-0.125_dp, -0.0625_dp, 0.0_dp, 0.0625_dp, 0.125_dp]
    call solve_transport(state, ok)
    if (.not. ok) error stop 'linear solved fixture'
    dsource_cos = 0.0_dp; dsource_sin = 0.0_dp
    call solve_transport_tangent(state, 1.0_dp, dsource_cos, dsource_sin, &
        dcos, dsin, ok)
    if (.not. ok) error stop 'linear tangent fixture'
    do evaluator = 1, 2
        call certify_transport_tangent(state, 1.0_dp, dsource_cos, dsource_sin, &
            dcos, dsin, certificate, ok, strict=evaluator == 2)
        if (.not. ok) error stop 'linear certificate preparation'
        call certify_transport_linear(state, certificate, basis, geometry, &
            velocity, work, lower, upper, ok)
        if (.not. ok) error stop 'continuous linear state band'
        call certify_transport_tangent_linear(state, certificate, dcos, dsin, &
            basis, geometry, velocity, work, tangent_lo, tangent_hi, ok)
        if (.not. ok) error stop 'continuous linear tangent band'
        do i = 1, points - 1
            do r = 0, 18
                t = real(r, qp)/18
                xq = (1 - t)*real(geometry(i)%lo, qp) + &
                    t*real(geometry(i + 1)%lo, qp)
                do j = 0, 12
                    vq = (1 - real(j, qp)/12)*real(velocity%lo, qp) + &
                        real(j, qp)/12*real(velocity%hi, qp)
                    call independent_resolvent(xq, vq, value, tangent)
                    call linear_contains(geometry, lower, upper, i, t, value)
                    call linear_contains(geometry, tangent_lo, tangent_hi, &
                        i, t, tangent)
                end do
            end do
        end do
    end do

    bad = geometry; bad(2)%lo = -0.1_dp; bad(2)%hi = -0.1_dp
    call transport_linear_enclosure(state, state%cosine, state%sine, basis, bad, &
        velocity, work, 0.0_dp, lower, upper, ok)
    if (ok) error stop 'stale point geometry accepted'
    bad = geometry; bad(2)%hi = 0.0_dp
    call prepare_transport_basis(basis, bad, ok)
    if (.not. ok) error stop 'finite interval basis fixture'
    call transport_linear_enclosure(state, state%cosine, state%sine, basis, bad, &
        velocity, work, 0.0_dp, lower, upper, ok)
    if (ok) error stop 'nonpoint linear geometry accepted'
    bad = geometry; bad(2) = bad(1)
    call prepare_transport_basis(basis, bad, ok)
    call transport_linear_enclosure(state, state%cosine, state%sine, basis, bad, &
        velocity, work, 0.0_dp, lower, upper, ok)
    if (ok) error stop 'duplicate x nodes accepted'
    bad = geometry; bad(2) = interval(-2.0_dp)
    call prepare_transport_basis(basis, bad, ok)
    call transport_linear_enclosure(state, state%cosine, state%sine, basis, bad, &
        velocity, work, 0.0_dp, lower, upper, ok)
    if (ok) error stop 'decreasing x nodes accepted'
    bad = geometry; bad(1)%lo = nan; bad(1)%hi = nan
    call prepare_transport_basis(basis, bad, ok)
    if (ok) error stop 'NaN point basis accepted'
    call transport_linear_enclosure(state, state%cosine, state%sine, basis, &
        geometry, velocity, work, 0.0_dp, lower, upper, ok)
    if (ok) error stop 'failed preparation retained linear validity'
    call prepare_transport_basis(basis, geometry, ok)
    if (.not. ok) error stop 'restore valid linear basis'
    call transport_linear_enclosure(state, state%cosine, state%sine, basis, &
        geometry, velocity, work, -1.0_dp, lower, upper, ok)
    if (ok) error stop 'negative supplied radius accepted'
    call transport_linear_enclosure(state, state%cosine, state%sine, basis, &
        geometry, velocity, work, nan, lower, upper, ok)
    if (ok) error stop 'NaN supplied radius accepted'
    call transport_linear_enclosure(state, state%cosine, state%sine, basis, &
        geometry, interval(nan), work, 0.0_dp, lower, upper, ok)
    if (ok) error stop 'NaN velocity accepted'
    call transport_linear_enclosure(state, state%cosine, state%sine, basis, &
        geometry, interval(-1.01_dp), work, 0.0_dp, lower, upper, ok)
    if (ok) error stop 'outside velocity domain accepted'
    call transport_linear_enclosure(state, state%cosine, state%sine, basis, &
        geometry, velocity, work, 0.0_dp, lower(:3), upper, ok)
    if (ok) error stop 'wrong lower shape accepted'
    call transport_linear_enclosure(state, state%cosine(:2, :), state%sine(:2, :), &
        basis, geometry, velocity, work, 0.0_dp, lower, upper, ok)
    if (ok) error stop 'wrong coefficient modes accepted'
    call transport_linear_enclosure(state, state%cosine, state%sine, missing, &
        geometry, velocity, work, 0.0_dp, lower, upper, ok)
    if (ok) error stop 'uninitialized basis accepted'
    call transport_linear_enclosure(state, state%cosine, state%sine, basis, &
        geometry, velocity, empty_work, 0.0_dp, lower, upper, ok)
    if (ok) error stop 'uninitialized workspace accepted'
    call transport_linear_enclosure(empty, state%cosine, state%sine, basis, &
        geometry, velocity, work, 0.0_dp, lower, upper, ok)
    if (ok) error stop 'uninitialized state accepted'
    saved = state%cosine(1, 2); state%cosine(1, 2) = nan
    call transport_linear_enclosure(state, state%cosine, state%sine, basis, &
        geometry, velocity, work, 0.0_dp, lower, upper, ok)
    if (ok) error stop 'NaN touched coefficient accepted'
    state%cosine(1, 2) = saved
    state%nu = 2.0_dp
    call certify_transport_linear(state, certificate, basis, geometry, velocity, &
        work, lower, upper, ok)
    if (ok) error stop 'stale certificate nu accepted'
    state%nu = 1.0_dp
    certificate%tangent_valid = .false.
    call certify_transport_tangent_linear(state, certificate, dcos, dsin, basis, &
        geometry, velocity, work, tangent_lo, tangent_hi, ok)
    if (ok) error stop 'missing tangent certificate accepted'
    ! Finite inputs can still overflow the chord estimate; never return a band.
    bad = [interval(0.0_dp), interval(1.0e155_dp), &
           interval(2.0e155_dp), interval(3.0e155_dp)]
    call prepare_transport_basis(basis, bad, ok)
    if (.not. ok) error stop 'finite large-phase basis fixture'
    call transport_linear_enclosure(state, state%cosine, state%sine, basis, bad, &
        interval(0.0_dp), work, 0.0_dp, lower, upper, ok)
    if (ok) error stop 'overflowing chord estimate accepted'
    ! The first small cells permit finite widening. Only the final large cell
    ! adds enough remainder to overflow radius .99*huge at a later vertex.
    bad = [interval(0.0_dp), interval(1.0_dp), &
           interval(2.0_dp), interval(1.0e154_dp)]
    call prepare_transport_basis(basis, bad, ok)
    if (.not. ok) error stop 'late overflow basis fixture'
    call transport_linear_enclosure(state, state%cosine, state%sine, basis, bad, &
        interval(0.0_dp), work, 0.99_dp*huge(0.0_dp), lower, upper, ok)
    if (ok) error stop 'later-node interval widening overflow accepted'
    print '(a)', 'continuous linear band and curvature interior oracles PASS'
contains
    subroutine linear_contains(x, lo, hi, cell, weight, exact)
        type(interval_t), intent(in) :: x(:)
        real(dp), intent(in) :: lo(:), hi(:)
        integer, intent(in) :: cell
        real(qp), intent(in) :: weight, exact
        real(qp) :: bottom, top

        if (cell < 1 .or. cell >= size(x)) error stop 'oracle cell invalid'
        ! Independent binary128 interpolation of the bounds actually plotted.
        bottom = (1 - weight)*real(lo(cell), qp) + weight*real(lo(cell + 1), qp)
        top = (1 - weight)*real(hi(cell), qp) + weight*real(hi(cell + 1), qp)
        if (exact < bottom .or. exact > top) &
            error stop 'linear plot band misses independent interior value'
    end subroutine

    subroutine check_candidate(a, b, x, lo, hi, vrange)
        real(dp), intent(in) :: a(0:, :), b(0:, :), lo(:), hi(:)
        type(interval_t), intent(in) :: x(:), vrange
        real(qp) :: weight, xq, vq, exact
        integer :: cell, p, q

        do cell = 1, size(x) - 1
            do p = 0, 18
                weight = real(p, qp)/18
                xq = (1 - weight)*real(x(cell)%lo, qp) + &
                    weight*real(x(cell + 1)%lo, qp)
                do q = 0, 20
                    vq = (1 - real(q, qp)/20)*real(vrange%lo, qp) + &
                        real(q, qp)/20*real(vrange%hi, qp)
                    exact = independent_candidate(state%v, a, b, xq, vq)
                    call linear_contains(x, lo, hi, cell, weight, exact)
                end do
            end do
        end do
    end subroutine

    function independent_candidate(grid, a, b, x, v) result(value)
        ! TEST ONLY binary128 interpolation and Fourier evaluation.
        real(dp), intent(in) :: grid(:), a(0:, :), b(0:, :)
        real(qp), intent(in) :: x, v
        real(qp) :: value, weight, aq, bq
        integer :: cell, k

        cell = 1
        do while (cell < size(grid) - 1)
            if (v < real(grid(cell + 1), qp)) exit
            cell = cell + 1
        end do
        weight = (v - real(grid(cell), qp))/ &
            (real(grid(cell + 1), qp) - real(grid(cell), qp))
        value = 0.0_qp
        do k = 0, size(a, 1) - 1
            aq = (1 - weight)*real(a(k, cell), qp) + weight*real(a(k, cell + 1), qp)
            bq = (1 - weight)*real(b(k, cell), qp) + weight*real(b(k, cell + 1), qp)
            value = value + aq*cos(real(k, qp)*x) + bq*sin(real(k, qp)*x)
        end do
    end function

    subroutine independent_resolvent(x, v, value, tangent)
        ! Independent complex Fourier resolvent, including its nu derivative.
        real(qp), intent(in) :: x, v
        real(qp), intent(out) :: value, tangent
        complex(qp) :: r1, r2

        r1 = cmplx((1 - v)/4, 0.0_qp, qp)/cmplx(1.0_qp, v, qp)
        r2 = cmplx(0.0_qp, -v/8, qp)/cmplx(1.0_qp, 2*v, qp)
        value = 1 + v/8 + real(r1*exp(cmplx(0.0_qp, x, qp)), qp) + &
            real(r2*exp(cmplx(0.0_qp, 2*x, qp)), qp)
        tangent = -1 - v/8 - &
            real(r1/cmplx(1.0_qp, v, qp)*exp(cmplx(0.0_qp, x, qp)), qp) - &
            real(r2/cmplx(1.0_qp, 2*v, qp)*exp(cmplx(0.0_qp, 2*x, qp)), qp)
    end subroutine
end program test_transport_linear_bands
