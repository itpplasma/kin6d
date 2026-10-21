program benchmark_nested_reference
    ! Independent diagnostic: external samples, never solver forcing.
    ! Rows: s, x(3), B(3), mu0*p, mu0*J(3), mu0*grad(p)(3).
    ! s is normalized toroidal flux; only the integer family is admitted.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use kin6d_nested_oracles, only: nested_oracle_t
    implicit none
    type(nested_oracle_t) :: oracle
    real(dp) :: row(14), exact(20), force(3), sums(6), label_error
    integer :: unit, status, count
    character(1024) :: path, argument
    if (command_argument_count() /= 4) &
        error stop 'usage: benchmark_nested_reference samples eps tau delta'
    call get_command_argument(1, path)
    call get_command_argument(2, argument); read (argument, *) oracle%eps
    call get_command_argument(3, argument); read (argument, *) oracle%tau
    call get_command_argument(4, argument); read (argument, *) oracle%delta
    if (.not. oracle%valid()) error stop 'invalid oracle domain'
    open (newunit=unit, file=trim(path), status='old', action='read')
    sums = 0; count = 0; label_error = 0
    do
        read (unit, *, iostat=status) row
        if (status < 0) exit
        if (status /= 0) error stop 'malformed sample'
        if (.not. all(ieee_is_finite(row))) error stop 'nonfinite sample'
        if (row(1) <= 0 .or. row(1) > 1) error stop 'invalid flux coordinate'
        call oracle%point(row(2:4), exact, status)
        if (status /= 0) error stop 'sample outside analytic field domain'
        force = [row(10)*row(7) - row(11)*row(6), &
                 row(11)*row(5) - row(9)*row(7), &
                 row(9)*row(6) - row(10)*row(5)] - row(12:14)
        sums = sums + [sum((row(5:7) - exact(3:5))**2), sum(exact(3:5)**2), &
            (row(8) - exact(2))**2, exact(2)**2, &
            sum(force**2), sum(exact(9:11)**2)]
        label_error = max(label_error, abs(exact(1)/oracle%delta - row(1)))
        count = count + 1
    end do
    close (unit)
    if (count == 0 .or. minval(sums([2, 4, 6])) <= 0) error stop 'empty norm scale'
    print '(a,i0)', 'samples=', count
    print '(a,es24.16)', 'B_sample_rms_relative=', sqrt(sums(1)/sums(2))
    print '(a,es24.16)', 'p_sample_rms_relative=', sqrt(sums(3)/sums(4))
    print '(a,es24.16)', 'surface_label_max_relative=', label_error
    print '(a,es24.16)', 'force_sample_rms_relative=', sqrt(sums(5)/sums(6))
    print '(a)', 'state_bound=unavailable; sampled diagnostics; no volume-norm claim'
end program
