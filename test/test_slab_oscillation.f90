program test_slab_oscillation
    use, intrinsic :: iso_fortran_env, only: dp => real64, qp => real128
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use fortnum_interval, only: interval_t, interval, contains, operator(+)
    use kin6d_slab_oscillation
    use kin6d_slab_oscillation_leaves, only: slab_oscillation_rhs_jvp, &
                                      slab_oscillation_rhs_vjp, &
                                      slab_oscillation_residual_interval, &
                                      slab_oscillation_hermite_interval
    implicit none
    type(slab_oscillation_state_t) :: state, coarse, fine, tangent, zero
    type(interval_t) :: qi, pi, qpp, residual_coefficients(4), midpoint_residual
    real(dp) :: q, p, qlo, qhi, plo, phi, tau, radius, rc, rf, padding
    real(dp) :: time(65), qtrace(65), ptrace(65), radii(65), residual(64)
    real(dp) :: dq, dm, gt, gq, gp, tangent_q, tangent_p, qlower, qupper
    real(qp) :: reference_q, reference_p, x
    integer :: i, j
    logical :: ok

    ! Arbitrary initial displacement/velocity, independently evaluated in binary128.
    call initialize_slab_oscillation(state, 64, 7.0_dp, 0.375_dp, -0.625_dp, ok)
    if (.not. ok) error stop 'general initial state rejected'
    do i = 0, 200
        tau = 7.0_dp*real(i, dp)/200.0_dp
        call sample_slab_oscillation(state, tau, q, p, qlo, qhi, plo, phi, ok)
        if (.not. ok) error stop 'point certificate failed'
        reference_q = 0.375_qp*cos(real(tau, qp)) - 0.625_qp*sin(real(tau, qp))
        reference_p = -0.625_qp*cos(real(tau, qp)) - 0.375_qp*sin(real(tau, qp))
        call enclosed(qlo, qhi, reference_q)
        call enclosed(plo, phi, reference_p)
    end do
    call sample_slab_oscillation_interval(state, interval(1.125_dp, 1.75_dp), &
                                           qi, pi, ok)
    if (.not. ok) error stop 'whole time interval failed'
    do j = 0, 64
        x = 1.125_qp + 0.625_qp*real(j, qp)/64.0_qp
        call enclosed(qi%lo, qi%hi, 0.375_qp*cos(x) - 0.625_qp*sin(x))
        call enclosed(pi%lo, pi%hi, -0.625_qp*cos(x) - 0.375_qp*sin(x))
    end do
    call initialize_slab_oscillation(coarse, 40, 8.0_dp, 1.0_dp, 0.0_dp, ok)
    if (.not. ok) error stop 'coarse trace'
    call initialize_slab_oscillation(fine, 80, 8.0_dp, 1.0_dp, 0.0_dp, ok)
    if (.not. ok) error stop 'fine trace'
    call slab_oscillation_radius(coarse, rc, ok)
    call slab_oscillation_radius(fine, rf, ok)
    if (rf >= 0.35_dp*rc) error stop 'residual radius lacks h-squared refinement'
    call initialize_slab_oscillation(zero, 16, 2.0_dp, 0.0_dp, 0.0_dp, ok)
    if (.not. ok) error stop 'zero state'
    call sample_slab_oscillation(zero, 1.0_dp, q, p, qlo, qhi, plo, phi, ok)
    if (.not. ok) error stop 'zero certificate'
    if (abs(q) > 0.0_dp .or. abs(p) > 0.0_dp) error stop 'zero solution moved'
    call enclosed(qlo, qhi, 0.0_qp)
    call enclosed(plo, phi, 0.0_qp)

    ! The unsolved candidate q=t^3,p=3t^2 has r=6t+t^3. Its exact cubic
    ! Bernstein residual coefficients on [0,1] are [0,2,4,7], midpoint25/8.
    call slab_oscillation_residual_interval(interval(0.0_dp), interval(1.0_dp), &
               interval(0.0_dp), interval(1.0_dp), interval(0.0_dp), interval(3.0_dp), &
                                   residual_coefficients(1), residual_coefficients(2), &
                                     residual_coefficients(3), residual_coefficients(4))
    do i = 1, 4
        select case (i)
        case (1)
            radius = 0.0_dp
        case (2)
            radius = 2.0_dp
        case (3)
            radius = 4.0_dp
        case (4)
            radius = 7.0_dp
        end select
        if (.not. contains(residual_coefficients(i), radius)) &
            error stop 'independent Bernstein polynomial excluded'
    end do

    call slab_oscillation_hermite_interval(interval(0.5_dp), interval(0.0_dp), &
                                          interval(1.0_dp), interval(0.0_dp), &
                                          interval(1.0_dp), interval(0.0_dp), &
                                          interval(3.0_dp), qi, pi, qpp)
    midpoint_residual = qi + qpp
    if (.not. contains(midpoint_residual, 3.125_dp)) &
        error stop 'exact cubic interior residual excluded'

    call get_slab_oscillation_trace(state, time, qtrace, ptrace, radii, ok)
    if (.not. ok) error stop 'trace getter'
    call certify_slab_oscillation_trace(time, qtrace, ptrace, 0.375_dp, -0.625_dp, &
                                        radii, residual, ok)
    if (.not. ok) error stop 'raw certificate'
    rc = radii(65)
    qtrace(31) = qtrace(31) + 0.5_dp
    call certify_slab_oscillation_trace(time, qtrace, ptrace, 0.375_dp, -0.625_dp, &
                                        radii, residual, ok)
    if (.not. ok) error stop 'finite corrupted candidate needs a valid broad bound'
    if (radii(65) <= 1.0_dp .or. radii(65) <= rc) error stop 'corruption undetected'
    time(32) = time(31)
    call certify_slab_oscillation_trace(time, qtrace, ptrace, 0.375_dp, -0.625_dp, &
                                        radii, residual, ok)
    if (ok) error stop 'duplicate time accepted'

    ! Equation-level tangent and adjoint candidates have separate continuous
    ! accuracy tests, beyond the generated local transpose identity.
    call slab_oscillation_rhs_jvp(0.0_dp, 0.375_dp, -0.625_dp, &
                                  1.0_dp, 0.25_dp, -0.5_dp, dq, dm)
    if (abs(dq + 0.5_dp) > 0.0_dp .or. abs(dm + 0.25_dp) > 0.0_dp) &
        error stop 'RHS tangent rational oracle'
    call slab_oscillation_rhs_vjp(0.0_dp, 0.375_dp, -0.625_dp, &
                                  0.75_dp, -0.25_dp, gt, gq, gp)
    if (abs(gt) > 0.0_dp .or. abs(gq - 0.25_dp) > 0.0_dp .or. &
        abs(gp - 0.75_dp) > 0.0_dp) error stop 'RHS adjoint rational oracle'
    call initialize_slab_oscillation_tangent(state, 0.25_dp, -0.5_dp, tangent, ok)
    if (.not. ok) error stop 'variational IVP'
    call sample_slab_oscillation(tangent, 7.0_dp, tangent_q, tangent_p, &
                                 qlo, qhi, plo, phi, ok)
    call enclosed(qlo, qhi, 0.25_qp*cos(7.0_qp) - 0.5_qp*sin(7.0_qp))
    call enclosed(plo, phi, -0.5_qp*cos(7.0_qp) - 0.25_qp*sin(7.0_qp))
    call slab_oscillation_terminal_vjp(state, 0.75_dp, -0.25_dp, gq, gp, ok)
    if (.not. ok) error stop 'adjoint IVP candidate'
    reference_q = 0.75_qp*cos(7.0_qp) + 0.25_qp*sin(7.0_qp)
    reference_p = 0.75_qp*sin(7.0_qp) - 0.25_qp*cos(7.0_qp)
    if (abs(real(gq, qp) - reference_q) > 1.0e-7_qp .or. &
        abs(real(gp, qp) - reference_p) > 1.0e-7_qp) error stop 'adjoint accuracy'
    if (abs(0.25_dp*gq - 0.5_dp*gp - &
            (0.75_dp*tangent_q - 0.25_dp*tangent_p)) > 1.0e-12_dp) &
        error stop 'frozen linear RK transpose pairing'

    ! Connected sample bands use a global residual radius plus a W2,infinity
    ! chord margin for the candidate, including across its C1 reconstruction knots.
    call slab_oscillation_radius(state, radius, ok)
    call slab_oscillation_chord(state, 1.0_dp, 2.0_dp, padding, ok)
    if (.not. ok) error stop 'chord query'
    call sample_slab_oscillation_candidate_interval(state, interval(1.0_dp), qi, pi, ok)
    qlower = qi%lo - radius - padding
    qupper = qi%hi + radius + padding
    call sample_slab_oscillation_candidate_interval(state, interval(2.0_dp), qi, pi, ok)
    do j = 0, 64
        x = real(j, qp)/64.0_qp
        reference_q = 0.375_qp*cos(1.0_qp + x) - 0.625_qp*sin(1.0_qp + x)
        call enclosed(real((1 - x)*real(qlower, qp) + &
                           x*(real(qi%lo, qp) - radius - padding), dp), &
                      real((1 - x)*real(qupper, qp) + &
                           x*(real(qi%hi, qp) + radius + padding), dp), reference_q)
    end do
    call solve_slab_oscillation(state, ok, certify=.false.)
    call sample_slab_oscillation(state, 1.0_dp, q, p, qlo, qhi, plo, phi, ok)
    if (ok) error stop 'candidate-only state advertised as certified'
    call initialize_slab_oscillation(state, 0, 1.0_dp, 1.0_dp, 0.0_dp, ok)
    if (ok) error stop 'zero steps accepted'
    call initialize_slab_oscillation(state, 16, -1.0_dp, 1.0_dp, 0.0_dp, ok)
    if (ok) error stop 'negative duration accepted'
    call initialize_slab_oscillation(state, 16, 1.0_dp, &
                                     ieee_value(1.0_dp, ieee_quiet_nan), 0.0_dp, ok)
    if (ok) error stop 'NaN initial condition accepted'
    print '(a)', 'Cold planar slab numerical/residual/tangent oracles PASS'

contains

    subroutine enclosed(lower, upper, reference)
        real(dp), intent(in) :: lower, upper
        real(qp), intent(in) :: reference

        if (real(lower, qp) > reference .or. real(upper, qp) < reference) &
            error stop 'independent harmonic oscillator excluded'
    end subroutine enclosed

end program test_slab_oscillation
