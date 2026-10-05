program test_poisson
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_is_finite
    use fortnum_interval, only: interval_t, interval
    use poisson_laboratory, only: poisson_mode_t, coupled_steady_t, &
                                  solve_poisson_mode, poisson_sample, &
                                  solve_coupled_steady, coupled_sample, &
                                  coupled_reference, coupled_certificate_radius, &
                                  coupled_certificate_accepts
    use kin6d_poisson_leaves, only: poisson_coefficients_jvp, &
                                   poisson_coefficients_vjp, &
                                   poisson_reconstruct_interval, &
                                   coupled_distribution_interval, &
                                   coupled_amplitude_jvp, coupled_amplitude_vjp, &
                                   coupled_density, coupled_collision, &
                                   coupled_collision_jvp, coupled_collision_vjp, &
                                   coupled_source_jvp
    use kin6d_leaves, only: coordinate
    implicit none
    real(dp), parameter :: pi = 3.141592653589793238462643383279502884_dp
    real(dp), parameter :: tolerance = 2.0e-14_dp
    type(poisson_mode_t) :: field
    type(coupled_steady_t) :: state, krook, changed, pure_fp
    type(interval_t) :: energy, cell_f, cell_phi, cell_e
    real(dp) :: f, e, p, source, residual, fr, pr, er, ref_f, ref_e, ref_p
    real(dp) :: pc, ps, ec, es, gk, geps, grc, grs, da, gd, gn, gb, gt
    real(dp) :: density_mean, density_cos, collision, nan, x, v
    logical :: ok
    integer :: i, j
    nan = ieee_value(0.0_dp, ieee_quiet_nan)
    ! Independent Fourier sign/scale fixture: k=2, epsilon0=3, rho=(12,6).
    call solve_poisson_mode(2, 3.0_dp, 0.0_dp, 12.0_dp, 6.0_dp, field, ok)
    if (.not. ok) error stop 'compatible Poisson mode rejected'
    call expect(field%phi_mean, 0.0_dp, 'zero potential gauge')
    call expect(field%phi_cos, 1.0_dp, 'Poisson cosine scale')
    call expect(field%phi_sin, 0.5_dp, 'Poisson sine scale')
    call expect(field%e_cos, -1.0_dp, 'negative potential derivative cosine')
    call expect(field%e_sin, 2.0_dp, 'negative potential derivative sine')
    call poisson_sample(field, 0.0_dp, p, e, ok)
    if (.not. ok) error stop 'valid field evaluation rejected'
    call expect(p, 1.0_dp, 'potential at zero')
    call expect(e, -1.0_dp, 'field at zero')
    call poisson_sample(field, pi/4.0_dp, p, e, ok)
    call expect(p, 0.5_dp, 'quarter-mode potential')
    call expect(e, 2.0_dp, 'quarter-mode field')
    call solve_poisson_mode(2, 3.0_dp, 1.0e-300_dp, 12.0_dp, 6.0_dp, field, ok)
    if (ok) error stop 'nonzero net charge silently removed'
    call solve_poisson_mode(0, 3.0_dp, 0.0_dp, 12.0_dp, 6.0_dp, field, ok)
    if (ok) error stop 'zero mode inverted'
    call solve_poisson_mode(-2, 3.0_dp, 0.0_dp, 12.0_dp, 6.0_dp, field, ok)
    if (ok) error stop 'noncanonical negative mode accepted'
    call solve_poisson_mode(2, 0.0_dp, 0.0_dp, 12.0_dp, 6.0_dp, field, ok)
    if (ok) error stop 'nonpositive permittivity accepted'
    call solve_poisson_mode(2, 3.0_dp, 0.0_dp, nan, 6.0_dp, field, ok)
    if (ok) error stop 'NaN charge accepted'
    ! Analytic differential of the mode resolvent: unit cosine input tangent.
    call poisson_coefficients_jvp(2.0_dp, 3.0_dp, 12.0_dp, 6.0_dp, &
                                   0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, pc, ps, ec, es)
    call expect(pc, 1.0_dp/12.0_dp, 'independent Poisson tangent potential')
    call expect(ps, 0.0_dp, 'Poisson tangent sine potential')
    call expect(ec, 0.0_dp, 'Poisson tangent cosine field')
    call expect(es, 1.0_dp/6.0_dp, 'independent Poisson tangent field')
    call poisson_coefficients_vjp(2.0_dp, 3.0_dp, 12.0_dp, 6.0_dp, &
                                   1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, gk, geps, grc, grs)
    call expect(gk, -1.0_dp, 'independent Poisson transpose wavenumber')
    call expect(geps, -1.0_dp/3.0_dp, 'independent Poisson transpose epsilon')
    call expect(grc, 1.0_dp/12.0_dp, 'independent Poisson transpose charge')
    call expect(grs, 0.0_dp, 'independent Poisson transpose sine charge')
    ! Direct coupled FP/Krook fixture. Moment source drive=7/32 gives a=1/4.
    call solve_coupled_steady(7.0_dp/32.0_dp, 0.5_dp, 2, state, ok, &
                              fp_rate=0.25_dp, temperature=0.5_dp)
    if (.not. ok) error stop 'manufactured coupled VPFP solve rejected'
    call expect(state%amplitude, 0.25_dp, 'source-moment amplitude solve')
    call expect(state%field%phi_cos, 1.0_dp/16.0_dp, 'self-consistent potential')
    call expect(state%field%e_sin, 1.0_dp/8.0_dp, 'self-consistent electric field')
    call coupled_density(state%amplitude, density_mean, density_cos)
    call expect(density_mean, 1.0_dp, 'exact finite-domain polynomial density')
    call expect(density_cos, 0.25_dp, 'exact active-population charge')
    call coupled_sample(state, 0.0_dp, 0.5_dp, f, e, p, source, residual, ok)
    if (.not. ok) error stop 'coupled sample rejected'
    call expect(f, 45.0_dp/64.0_dp, 'independent rational distribution fixture')
    call expect(source, 63.0_dp/256.0_dp, 'independent rational FP source fixture')
    call coupled_collision(0.25_dp, 2.0_dp, 0.0_dp, 0.5_dp, 0.25_dp, 0.5_dp, &
                            collision)
    call expect(collision, -45.0_dp/256.0_dp, 'independent OU collision fixture')
    call coupled_sample(state, pi/4.0_dp, 0.5_dp, f, e, p, source, residual, ok)
    call expect(f, 9.0_dp/16.0_dp, 'quarter-mode distribution')
    call expect(source, -3.0_dp/32.0_dp, 'independent electric/streaming source')
    call coupled_sample(state, pi/4.0_dp, -0.5_dp, f, e, p, source, residual, ok)
    if (.not. ok) error stop 'negative-velocity polynomial domain lost'
    call expect(source, 3.0_dp/8.0_dp, 'independent negative-velocity FP source')
    call coupled_sample(state, 0.0_dp, 0.0_dp, f, e, p, source, residual, ok)
    if (.not. ok) error stop 'zero-velocity polynomial domain lost'
    call expect(source, 3.0_dp/32.0_dp, 'independent zero-velocity FP source')
    call coupled_collision_jvp(0.25_dp, 2.0_dp, 0.0_dp, 0.0_dp, &
                                0.25_dp, 0.5_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                                1.0_dp, 0.0_dp, 0.0_dp, da)
    call expect(da, 0.0_dp, 'polynomial collision tangent at zero velocity')
    call coupled_collision_jvp(0.25_dp, 2.0_dp, 0.0_dp, -0.5_dp, &
                                0.25_dp, 0.5_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                                1.0_dp, 0.0_dp, 0.0_dp, da)
    call expect(da, 45.0_dp/64.0_dp, 'independent negative-velocity collision tangent')
    call coupled_collision_vjp(0.25_dp, 2.0_dp, 0.0_dp, -0.5_dp, &
                                0.25_dp, 0.5_dp, 1.0_dp, gd, gk, gn, geps, gb, gt)
    call expect(gd, -9.0_dp/64.0_dp, 'independent collision amplitude adjoint')
    call expect(gk, 0.0_dp, 'independent collision mode adjoint at zero coordinate')
    call expect(gn, 0.0_dp, 'independent collision coordinate adjoint')
    call expect(geps, 45.0_dp/64.0_dp, 'independent collision velocity adjoint')
    call expect(gb, -45.0_dp/64.0_dp, 'independent collision rate adjoint')
    call expect(gt, -15.0_dp/32.0_dp, 'independent collision temperature adjoint')
    call coupled_source_jvp(7.0_dp/32.0_dp, 0.5_dp, 2.0_dp, pi/4.0_dp, &
                             0.0_dp, 0.25_dp, 0.5_dp, 0.0_dp, 0.0_dp, &
                             0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, da)
    call expect(da, -9.0_dp/16.0_dp, 'independent source tangent at zero velocity')
    ! At k*x=pi/4 the quadratic sin(2*k*x) term is nonzero.
    call coupled_sample(state, pi/8.0_dp, 0.5_dp, f, e, p, source, residual, ok)
    call expect(source, 33.0_dp/256.0_dp*(1.0_dp - &
                0.707106781186547524400844362104849039_dp), &
                'independent nonlinear second-harmonic source')
    call coupled_certificate_radius(state, fr, pr, er, energy, ok)
    if (.not. ok) error stop 'manufactured interval enclosure rejected'
    if (energy%lo > 1.0_dp/256.0_dp) error stop 'field energy lower bound missed'
    if (energy%hi < 1.0_dp/256.0_dp) error stop 'field energy upper bound missed'
    if (.not. coupled_certificate_accepts(state, fr, pr, er)) &
        error stop 'valid pointwise radii rejected'
    if (coupled_certificate_accepts(state, 0.0_dp, 0.0_dp, 0.0_dp)) &
        error stop 'zero-width certificate accepted'
    if (coupled_certificate_accepts(state, nan, pr, er)) &
        error stop 'NaN certificate accepted'
    changed = state; changed%field%phi_cos = 0.2_dp
    if (coupled_certificate_accepts(changed, fr, pr, er)) &
        error stop 'corrupted potential escaped enclosure check'
    changed = state; changed%field%e_sin = 0.2_dp
    if (coupled_certificate_accepts(changed, fr, pr, er)) &
        error stop 'corrupted electric field escaped enclosure check'
    changed = state; changed%amplitude = 0.5_dp
    if (coupled_certificate_accepts(changed, fr, pr, er)) &
        error stop 'corrupted kinetic state escaped charge compatibility'
    changed = state; changed%field%phi_mean = 0.1_dp
    if (coupled_certificate_accepts(changed, fr, pr, er)) &
        error stop 'nonzero gauge accepted'
    call coupled_distribution_interval(interval(state%amplitude), interval(2.0_dp), &
                                        interval(0.0_dp, pi), interval(0.0_dp), cell_f)
    if (any(.not. ieee_is_finite([cell_f%lo, cell_f%hi]))) &
        error stop 'nonfinite whole-cell distribution range'
    if (cell_f%lo > 9.0_dp/16.0_dp .or. cell_f%hi < 15.0_dp/16.0_dp) &
        error stop 'whole-cell distribution range missed extrema'
    call poisson_reconstruct_interval(interval(2.0_dp), interval(1.0_dp/16.0_dp), &
                                       interval(0.0_dp), interval(0.0_dp), &
                                       interval(1.0_dp/8.0_dp), interval(0.0_dp, pi), &
                                       cell_phi, cell_e)
    if (any(.not. ieee_is_finite([cell_phi%lo, cell_phi%hi, cell_e%lo, cell_e%hi]))) &
        error stop 'nonfinite whole-cell field range'
    if (cell_phi%lo > -1.0_dp/16.0_dp .or. cell_phi%hi < 1.0_dp/16.0_dp) &
        error stop 'whole-cell potential missed extrema'
    if (cell_e%lo > -1.0_dp/8.0_dp .or. cell_e%hi < 1.0_dp/8.0_dp) &
        error stop 'whole-cell electric field missed extrema'
    do j = 0, 16
        call coordinate(-1.0_dp, 1.0_dp, real(j, dp), 16.0_dp, v)
        do i = 0, 16
            call coordinate(0.0_dp, 2.0_dp*pi, real(i, dp), 16.0_dp, x)
            call coupled_sample(state, x, v, f, e, p, source, residual, ok)
            if (.not. ok) error stop 'finite-domain coupled grid sample failed'
            if (abs(residual) > tolerance) error stop 'full VPFP source balance'
            if (f < -tolerance) error stop 'manufactured distribution is negative'
            call coupled_reference(state, x, v, ref_f, ref_e, ref_p, ok)
            if (.not. ok) error stop 'manufactured reference rejected'
            if (abs(f - ref_f) > fr + tolerance) error stop 'distribution missed'
            if (abs(e - ref_e) > er + tolerance) error stop 'field missed'
            if (abs(p - ref_p) > pr + tolerance) error stop 'potential missed'
        end do
    end do
    call coupled_sample(state, 0.5_dp, -1.0_dp, f, e, p, source, residual, ok)
    call expect(f, 0.0_dp, 'left Dirichlet velocity trace')
    call coupled_sample(state, 0.5_dp, 1.0_dp, f, e, p, source, residual, ok)
    call expect(f, 0.0_dp, 'right Dirichlet velocity trace')
    call coupled_sample(state, 0.0_dp, 1.0001_dp, f, e, p, source, residual, ok)
    if (ok) error stop 'finite velocity domain silently extrapolated'
    ! Independent implicit derivative of (nu+3*beta*theta)*a=drive.
    call coupled_amplitude_jvp(7.0_dp/32.0_dp, 0.5_dp, 0.25_dp, 0.5_dp, &
                                 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, da)
    call expect(da, 8.0_dp/7.0_dp, 'independent coupled source tangent')
    call coupled_amplitude_vjp(7.0_dp/32.0_dp, 0.5_dp, 0.25_dp, 0.5_dp, &
                                 1.0_dp, gd, gn, gb, gt)
    call expect(gd, 8.0_dp/7.0_dp, 'independent coupled source adjoint')
    call expect(gn, -2.0_dp/7.0_dp, 'independent coupled Krook adjoint')
    call expect(gb, -3.0_dp/7.0_dp, 'independent coupled FP-rate adjoint')
    call expect(gt, -3.0_dp/14.0_dp, &
                'independent coupled bath-temperature adjoint')
    ! Retained beta=0 Krook limit and pure FP source-balanced state.
    call solve_coupled_steady(1.0_dp/8.0_dp, 0.5_dp, 2, krook, ok)
    if (.not. ok) error stop 'retained coupled Krook limit rejected'
    call expect(krook%amplitude, 0.25_dp, 'Krook limit source moment')
    call expect(krook%field%e_sin, state%field%e_sin, 'Krook retained field')
    call solve_coupled_steady(3.0_dp/32.0_dp, 0.0_dp, 2, pure_fp, ok, &
                              fp_rate=0.25_dp, temperature=0.5_dp)
    if (.not. ok) error stop 'pure FP source-balanced state rejected'
    call expect(pure_fp%amplitude, 0.25_dp, 'pure FP source moment')
    call solve_coupled_steady(0.0_dp, 0.0_dp, 2, changed, ok)
    if (ok) error stop 'singular undamped family silently inverted'
    call solve_coupled_steady(1.0_dp, 0.5_dp, 2, changed, ok)
    if (ok) error stop 'nonpositive distribution family accepted'
    call solve_coupled_steady(0.1_dp, 0.5_dp, 2, changed, ok, fp_rate=-0.25_dp)
    if (ok) error stop 'negative FP rate accepted'
    call solve_coupled_steady(0.1_dp, 0.5_dp, 2, changed, ok, temperature=0.0_dp)
    if (ok) error stop 'zero bath temperature accepted'
    call solve_coupled_steady(nan, 0.5_dp, 2, changed, ok)
    if (ok) error stop 'NaN drive accepted'
    call solve_coupled_steady(0.1_dp, 0.5_dp, 0, changed, ok)
    if (ok) error stop 'zero coupled mode accepted'
    print '(a)', 'PASS independent Poisson sign/gauge/charge and derivative fixtures'
    print '(a)', 'PASS direct coupled source-balanced FP/Krook and pure FP fixtures'
    print '(a,3es12.4)', 'PASS whole-domain pointwise f/phi/E radii: ', fr, pr, er
    print '(a)', 'PASS full-cell ranges, velocity traces, corrupted-input rejection'
contains
    subroutine expect(actual, expected, label)
        real(dp), intent(in) :: actual, expected
        character(*), intent(in) :: label
        if (.not. ieee_is_finite(actual)) then
            print '(a)', label
            error stop 'nonfinite value in independent fixture'
        end if
        if (abs(actual - expected) > tolerance) then
            print '(a,2es22.12)', label, actual, expected
            error stop 'independent analytical fixture failed'
        end if
    end subroutine
end program
