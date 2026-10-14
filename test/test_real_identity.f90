program test_real_identity
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, &
        ieee_positive_inf, ieee_negative_inf
    use kin6d_real_identity, only: same_real_value
    implicit none
    real(dp) :: nan, pinf, ninf, subnormal

    nan = ieee_value(0.0_dp, ieee_quiet_nan)
    pinf = ieee_value(0.0_dp, ieee_positive_inf)
    ninf = ieee_value(0.0_dp, ieee_negative_inf)
    subnormal = nearest(0.0_dp, 1.0_dp)
    if (.not. same_real_value(0.0_dp, -0.0_dp)) error stop 'signed zeros'
    if (.not. same_real_value(1.0_dp, 1.0_dp)) error stop 'identical input'
    if (same_real_value(1.0_dp, nearest(1.0_dp, 1.0_dp))) &
        error stop 'one-ULP input change must invalidate cache'
    if (same_real_value(0.0_dp, subnormal)) error stop 'retain subnormal input'
    if (same_real_value(subnormal, -subnormal)) error stop 'retain sign'
    if (.not. same_real_value(pinf, pinf)) error stop 'positive infinity'
    if (.not. same_real_value(ninf, ninf)) error stop 'negative infinity'
    if (same_real_value(pinf, ninf)) error stop 'infinity sign'
    if (same_real_value(nan, nan)) error stop 'NaN cannot be a cache hit'
    if (same_real_value(nan, 1.0_dp)) error stop 'NaN input'
    if (same_real_value(1.0_dp, nan)) error stop 'NaN cached value'
    if (.not. all(same_real_value([0.0_dp, 2.0_dp], [-0.0_dp, 2.0_dp]))) &
        error stop 'array identity'
end program test_real_identity
