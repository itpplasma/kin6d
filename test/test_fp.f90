program test_fp
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use kin6d_fp_laboratory, only: fp_state_t, initialize_fp, solve_fp, sample_fp, &
                                 fp_certificate_radius, fp_certificate_accepts, &
                                 fp_sample_band, fp_reference_band, fp_mass_interval, &
                                 fp_normalization_interval
    use kin6d_fp_leaves, only: fp_flux, fp_flux_jvp, fp_flux_vjp, fp_step_jvp, &
                              fp_mass_error_band
    use fortnum_interval, only: interval_t, interval
    implicit none
    type(fp_state_t) :: coarse, fine, warm, invalid
    real(dp) :: rcoarse, rfine, residual, value, lo, hi, rlo, rhi, tangent
    real(dp) :: gv, gt, gh, ga, gb
    integer :: j
    type(interval_t) :: mlower, mupper
    ! Exact fixtures from independently evaluated erf normalization:
    ! Z=sqrt(2*pi*theta)*erf(2/sqrt(2*theta)), f(v)=mass*M(v)/Z.
    ! Values retained to 34 digits; independent of the midpoint normalizer.
    real(dp), parameter :: exact0 = 0.41795955023513457256435545185990697_dp
    real(dp), parameter :: exact1 = 0.25350528173731172246311001951610784_dp
    real(dp), parameter :: exact2 = 0.056564674112519188205021581580777325_dp
    real(dp), parameter :: warm0 = 0.58436856725681664457100810290418899_dp
    real(dp), parameter :: exactz = 2.3925760266452164058628475409486855_dp
    real(dp), parameter :: exactmid = 0.41775551824004887499433640664745709_dp
    call initialize_fp(coarse, 32)
    call initialize_fp(fine, 64)
    call solve_fp(coarse, residual)
    call solve_fp(fine, residual)
    rcoarse = fp_certificate_radius(coarse)
    rfine = fp_certificate_radius(fine)
    if (rfine >= 0.27_dp*rcoarse) error stop 'FP enclosure refinement not quadratic'
    call sample_fp(fine, 0.0_dp, value)
    if (abs(value - exact0) > rfine) error stop 'FP equilibrium at zero missed exact fixture'
    if (abs(value - exact0) > 0.0001_dp) error stop 'FP direct candidate is inaccurate'
    call sample_fp(fine, 1.0_dp, value)
    if (abs(value - exact1) > rfine) error stop 'FP equilibrium at one missed exact fixture'
    call sample_fp(fine, 2.0_dp, value)
    if (abs(value - exact2) > rfine) error stop 'FP equilibrium at endpoint missed exact fixture'
    call fp_sample_band(fine, 0.0_dp, lo, hi)
    if (exact0 < lo .or. exact0 > hi) error stop 'FP pointwise band missed exact fixture'
    call fp_sample_band(fine, 0.03125_dp, lo, hi)
    if (exactmid < lo .or. exactmid > hi) error stop 'FP between-node band missed exact fixture'
    call fp_reference_band(fine, 0.0_dp, rlo, rhi)
    if (exact0 < rlo .or. exact0 > rhi) error stop 'FP rigorous normalization missed erf oracle'
    call fp_reference_band(fine, 1.0_dp, rlo, rhi)
    if (exact1 < rlo .or. exact1 > rhi) error stop 'FP reference at one missed erf oracle'
    call fp_mass_interval(fine, lo, hi)
    if (lo > 1.000000000001_dp .or. hi < 0.999999999999_dp) &
        error stop 'FP prescribed particle mass'
    call fp_mass_error_band(interval(lo), interval(hi), interval(rfine), mlower, mupper)
    if (mlower%lo > 1.0_dp .or. mupper%hi < 1.0_dp) &
        error stop 'FP continuous particle mass enclosure'
    call fp_normalization_interval(fine, lo, hi)
    if (exactz < lo .or. exactz > hi) error stop 'FP normalizer missed independent erf fixture'
    ! Even heat bath and symmetric reflecting domain fix an even equilibrium.
    do j = 1, size(fine%f)
        if (abs(fine%f(j) - fine%f(size(fine%f) + 1 - j)) > 1.0e-14_dp) &
            error stop 'FP equilibrium symmetry'
    end do
    ! Independent exact flux fixture: 3*(5-4)+2*(4+5)/2=12.
    call fp_flux(2.0_dp, 3.0_dp, 1.0_dp, 4.0_dp, 5.0_dp, value)
    if (abs(value - 12.0_dp) > 1.0e-14_dp) error stop 'FP flux exact fixture'
    call fp_flux_jvp(2.0_dp, 3.0_dp, 1.0_dp, 4.0_dp, 5.0_dp, &
                     1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, tangent)
    if (abs(tangent - 4.5_dp) > 1.0e-14_dp) error stop 'FP flux analytical tangent'
    call fp_flux_vjp(2.0_dp, 3.0_dp, 1.0_dp, 4.0_dp, 5.0_dp, &
                     2.0_dp, gv, gt, gh, ga, gb)
    if (max(abs(gv - 9.0_dp), abs(gt - 2.0_dp), abs(gh + 6.0_dp), &
            abs(ga + 4.0_dp), abs(gb - 8.0_dp)) > 1.0e-14_dp) &
        error stop 'FP flux analytical adjoint'
    ! d[a(2 theta-h v)/(2 theta+h v)]/dv=-3/4 for these exact inputs.
    call fp_step_jvp(2.0_dp, 3.0_dp, 1.0_dp, 4.0_dp, &
                     1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, tangent)
    if (abs(tangent + 0.75_dp) > 1.0e-14_dp) error stop 'FP recurrence analytical tangent'
    call initialize_fp(warm, 64, theta=4.0_dp, mass=2.0_dp)
    call solve_fp(warm, residual)
    call fp_sample_band(warm, 0.0_dp, lo, hi)
    if (warm0 < lo .or. warm0 > hi) error stop 'FP heat-bath temperature or mass scaling'
    if (.not. fp_certificate_accepts(fine, rfine)) error stop 'FP valid radius rejected'
    if (fp_certificate_accepts(fine, 0.0_dp)) error stop 'FP zero radius accepted'
    if (fp_certificate_accepts(fine, ieee_value(0.0_dp, ieee_quiet_nan))) &
        error stop 'FP NaN radius accepted'
    fine%f(33) = 3.0_dp
    if (fp_certificate_accepts(fine, rfine)) error stop 'FP corrupt candidate accepted'
    fine%theta = 0.0_dp
    if (fp_certificate_accepts(fine, rfine)) error stop 'FP invalid bath accepted'
    if (fp_certificate_accepts(invalid, rfine)) error stop 'FP unallocated candidate accepted'
    invalid = coarse
    invalid%v(2) = invalid%v(1)
    if (fp_certificate_accepts(invalid, rcoarse)) error stop 'FP collapsed cell accepted'
    invalid = coarse
    invalid%v(1) = -3.0_dp
    if (fp_certificate_accepts(invalid, rcoarse)) error stop 'FP changed boundary accepted'
    invalid = coarse
    invalid%f(1) = ieee_value(0.0_dp, ieee_quiet_nan)
    if (fp_certificate_accepts(invalid, rcoarse)) error stop 'FP NaN candidate accepted'
    write (*, '(a,2es15.6)') 'FP continuous Linf radii (32/64 cells): ', rcoarse, rfine
end program
