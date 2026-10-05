program test_laboratory
    use, intrinsic :: iso_fortran_env, only: dp => real64, qp => real128
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use kin6d_laboratory, only: laboratory_t, initialize, advance, sample, &
                                solve_steady, certificate_radius, certificate_accepts
    use kin6d_leaves, only: characteristic_reference_qp, stream_jvp, stream_vjp, &
                            stream_step, squared_norm, dot_four, steady_solve_jvp, &
                            steady_solve_vjp, steady_residual_jvp
    implicit none
    type(laboratory_t) :: coarse, fine
    real(dp) :: value, bound, coarse_bound, residual, measured, velocity, x
    real(dp) :: ja, jb, ga, gb, gv, gdt, left, right, an, bn, norm0, norm1
    real(qp) :: exact
    integer :: i, j
    call initialize(coarse, 64)
    call initialize(fine, 128)
    do i = 1, 64
        call advance(coarse, 0.0625_dp)
        call advance(fine, 0.0625_dp)
    end do
    coarse_bound = certificate_radius(coarse)
    bound = certificate_radius(fine)
    if (bound >= coarse_bound*0.26_dp) error stop 'velocity refinement not quadratic'
    measured = 0.0_dp
    do j = 1, 1025
        velocity = -4.0_dp + real(j - 1, dp)/128.0_dp
        do i = 0, 12
            x = real(i, dp)/2.0_dp
            call sample(fine, x, velocity, value)
            call characteristic_reference_qp(real(x, qp), real(velocity, qp), &
                                             real(fine%time, qp), &
                                             real(fine%amplitude, qp), exact)
            measured = max(measured, real(abs(real(value, qp) - exact), dp))
        end do
    end do
    if (measured > bound) error stop 'continuous snapshot bound missed exact flow'
    if (.not. certificate_accepts(fine, bound)) error stop 'valid bound rejected'
    if (certificate_accepts(fine, 0.0_dp)) error stop 'zero-width bound accepted'
    if (certificate_accepts(fine, ieee_value(0.0_dp, ieee_quiet_nan))) &
        error stop 'NaN bound accepted'
    fine%time = ieee_value(0.0_dp, ieee_quiet_nan)
    if (certificate_accepts(fine, bound)) error stop 'NaN time accepted'
    fine%time = 4.0_dp
    fine%amplitude = -0.5_dp
    if (certificate_accepts(fine, bound)) error stop 'invalid amplitude accepted'
    fine%amplitude = 0.5_dp
    fine%a(2) = 3.0_dp
    if (certificate_accepts(fine, bound)) error stop 'corrupted candidate accepted'
    ! Changed timestep and reverse flow invalidate the cached phase factors.
    call initialize(fine, 128)
    call advance(fine, 0.0625_dp)
    call advance(fine, -0.125_dp)
    coarse_bound = certificate_radius(fine)
    call sample(fine, 0.75_dp, 0.5_dp, value)
    call characteristic_reference_qp(0.75_qp, 0.5_qp, -0.0625_qp, 0.5_qp, exact)
    if (abs(real(value, qp) - exact) > real(coarse_bound, qp)) &
        error stop 'cached timestep change or reverse flow failed'
    call initialize(fine, 128)
    call solve_steady(fine, 0.75_dp, residual)
    if (residual > 2.0e-15_dp) error stop 'direct stationary solve residual'
    do i = 1, size(fine%v)
        if (abs(fine%a(i) - 0.5_dp*fine%m(i)) > 2.0e-15_dp) &
            error stop 'stationary manufactured Fourier amplitude'
        if (abs(fine%b(i)) > 2.0e-15_dp) error stop 'stationary sine mode'
    end do
    ! Closed streaming rotates coefficients; known input has squared norm 5.
    call stream_step(1.0_dp, 2.0_dp, 0.75_dp, 0.125_dp, an, bn)
    call squared_norm(1.0_dp, 2.0_dp, norm0)
    call squared_norm(an, bn, norm1)
    if (abs(norm1 - norm0) > 2.0e-15_dp) error stop 'streaming norm invariant'
    call stream_jvp(1.0_dp, 2.0_dp, 0.75_dp, 0.125_dp, &
                    0.3_dp, -0.5_dp, 0.2_dp, 0.1_dp, ja, jb)
    call stream_vjp(1.0_dp, 2.0_dp, 0.75_dp, 0.125_dp, 0.4_dp, 0.7_dp, &
                    ga, gb, gv, gdt)
    call dot_four(ja, jb, 0.0_dp, 0.0_dp, 0.4_dp, 0.7_dp, 0.0_dp, 0.0_dp, left)
    call dot_four(ga, gb, gv, gdt, 0.3_dp, -0.5_dp, 0.2_dp, 0.1_dp, right)
    if (abs(left - right) > 1.0e-14_dp) error stop 'JVP/VJP pairing'
    ! A separate centered parameter perturbation checks the derivative, too.
    call stream_step(1.0_dp, 2.0_dp, 0.750001_dp, 0.125_dp, an, bn)
    call stream_step(1.0_dp, 2.0_dp, 0.749999_dp, 0.125_dp, ga, gb)
    call stream_jvp(1.0_dp, 2.0_dp, 0.75_dp, 0.125_dp, &
                    0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, ja, jb)
    if (abs((an - ga)/0.000002_dp - ja) > 1.0e-9_dp) &
        error stop 'independent JVP perturbation'
    if (abs((bn - gb)/0.000002_dp - jb) > 1.0e-9_dp) &
        error stop 'independent JVP perturbation'
    ! Implicit tangent: A du = ds for a unit cosine-source perturbation.
    ! At nu=v=3/4 the independently known solution is du=(2/3,2/3).
    call steady_solve_jvp(0.75_dp, 0.75_dp, 0.3_dp, -0.4_dp, &
                          0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, ja, jb)
    if (abs(ja - 2.0_dp/3.0_dp) > 1.0e-14_dp) error stop 'steady tangent amplitude'
    if (abs(jb - 2.0_dp/3.0_dp) > 1.0e-14_dp) error stop 'steady tangent sine'
    call steady_residual_jvp(0.2_dp, -0.1_dp, 0.75_dp, 0.75_dp, &
                             0.3_dp, -0.4_dp, ja, jb, 0.0_dp, 0.0_dp, &
                             1.0_dp, 0.0_dp, an, bn)
    if (max(abs(an), abs(bn)) > 1.0e-14_dp) error stop 'implicit residual tangent'
    call steady_solve_jvp(0.75_dp, 0.75_dp, 0.3_dp, -0.4_dp, &
                          0.25_dp, 0.4_dp, -0.2_dp, 0.1_dp, ja, jb)
    call steady_solve_vjp(0.75_dp, 0.75_dp, 0.3_dp, -0.4_dp, &
                          0.4_dp, 0.7_dp, ga, gb, gv, gdt)
    call dot_four(ja, jb, 0.0_dp, 0.0_dp, 0.4_dp, 0.7_dp, 0.0_dp, 0.0_dp, left)
    call dot_four(ga, gb, gv, gdt, 0.25_dp, 0.4_dp, -0.2_dp, 0.1_dp, right)
    if (abs(left - right) > 1.0e-14_dp) error stop 'steady JVP/VJP pairing'

    print '(a,es12.4)', 'PASS continuous snapshot bound: ', bound
    print '(a,es12.4)', 'PASS independent real128 maximum sampled error: ', measured
    print '(a,es12.4)', 'PASS direct steady residual: ', residual
    print '(a)', 'PASS corrupted/zero/NaN rejection, refinement, invariant, derivatives'
end program
