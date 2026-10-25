module g6_external_comparison
    ! Equal-weight physical sample diagnostics, never a volume norm or bound.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, &
        ieee_quiet_nan
    use g6_sampled_force, only: g6_force_accumulator_t, g6_force_add, &
        g6_force_result
    implicit none
    private
    public :: g6_comparison_t, compare_external_samples
    type :: g6_comparison_t
        integer :: samples = 0
        real(dp) :: b_error = huge(1.0_dp), p_error = huge(1.0_dp)
        real(dp) :: p_tolerance = 0
        real(dp) :: force_rms = huge(1.0_dp)
        integer :: force_status = 1
        character(160) :: force_message = 'external sampled force unavailable'
        character(16) :: pressure_kind = 'unavailable'
        character(160) :: message = 'unavailable'
    end type
    type :: norm_t
        real(dp) :: scale = 0, ssq = 1
    end type
    abstract interface
        subroutine physical_sampler(xyz, b_cartesian, mu0_p, status)
            import dp
            real(dp), intent(in) :: xyz(3)
            real(dp), intent(out) :: b_cartesian(3), mu0_p
            integer, intent(out) :: status
        end subroutine
    end interface
contains
    subroutine compare_external_samples(path, sampler, report, status, &
                                        zero_pressure_atol)
        character(*), intent(in) :: path
        procedure(physical_sampler) :: sampler
        type(g6_comparison_t), intent(out) :: report
        integer, intent(out) :: status
        ! Absolute RMS tolerance for mu0*p, in T^2, only for p_ref == 0.
        real(dp), intent(in), optional :: zero_pressure_atol
        type(norm_t) :: bdiff, pdiff, bref, pref
        type(g6_force_accumulator_t) :: force_samples
        real(dp) :: row(14), extra, b(3), pressure, delta
        character(4096) :: line
        integer :: unit, ios, st, j
        logical :: ok
        report = g6_comparison_t()
        report%force_rms = ieee_value(0.0_dp, ieee_quiet_nan)
        status = 1
        open (newunit=unit, file=path, status='old', action='read', iostat=ios)
        if (ios /= 0) then
            report%message = 'G6 cannot open external samples'
            return
        end if
        report%message = ''
        do
            read (unit, '(a)', iostat=ios) line
            if (ios < 0) exit
            report%message = 'G6 malformed external sample: expected 14 values'
            if (ios /= 0) exit
            if (len_trim(line) == len(line)) exit
            row = ieee_value(0.0_dp, ieee_quiet_nan)
            read (line, *, iostat=ios) row
            if (ios /= 0) exit
            if (.not. all(ieee_is_finite(row))) then
                report%message = 'G6 nonfinite external sample'
                exit
            end if
            ! Reject extra columns as well as truncated/missing values.
            read (line, *, iostat=ios) row, extra
            if (ios >= 0) exit
            call sampler(row(2:4), b, pressure, st)
            if (st /= 0) then
                report%message = 'G6 external point outside domain'
                exit
            end if
            if (.not. all(ieee_is_finite(b)) .or. &
                .not. ieee_is_finite(pressure)) then
                report%message = 'G6 nonfinite sampled B or pressure'
                exit
            end if
            report%message = 'G6 invalid norm accumulation'
            ok = .true.
            do j = 1, 3
                call difference(b(j), row(4 + j), delta, ok)
                if (.not. ok) exit
                call norm_add(bdiff, delta)
                call norm_add(bref, row(4 + j))
            end do
            if (.not. ok) exit
            call difference(pressure, row(8), delta, ok)
            if (.not. ok) exit
            call norm_add(pdiff, delta)
            call norm_add(pref, row(8))
            call g6_force_add(force_samples, row(5:7), row(9:11), row(12:14), st)
            if (report%samples == huge(report%samples)) exit
            report%samples = report%samples + 1
            report%message = ''
        end do
        close (unit)
        if (len_trim(report%message) /= 0) return
        if (report%samples >= 3) then
            call g6_force_result(force_samples, report%force_rms, &
                                 report%force_status)
            if (report%force_status == 0) then
                report%force_message = 'external sampled force available; no force gate'
            else
                report%force_message = 'external sampled force arithmetic unavailable'
            end if
        end if
        report%message = 'G6 insufficient samples or zero magnetic scale'
        if (report%samples < 3 .or. bref%scale <= 0) return
        call norm_ratio(bdiff, bref, report%b_error, ok)
        report%message = 'G6 invalid norm accumulation'
        if (.not. ok) return
        if (pref%scale > 0) then
            report%pressure_kind = 'relative RMS'
            report%p_tolerance = 1.0e-3_dp
            call norm_ratio(pdiff, pref, report%p_error, ok)
            if (.not. ok) return
        else
            report%pressure_kind = 'absolute RMS'
            report%message = 'G6 zero pressure requires positive finite '// &
                'absolute mu0*p RMS tolerance (T^2)'
            if (.not. present(zero_pressure_atol)) return
            if (.not. ieee_is_finite(zero_pressure_atol)) return
            if (zero_pressure_atol <= 0) return
            report%p_tolerance = zero_pressure_atol
            report%p_error = pdiff%scale &
                *sqrt(pdiff%ssq/real(report%samples, dp))
        end if
        report%message = 'G6 invalid norm accumulation'
        if (.not. ieee_is_finite(report%b_error)) return
        if (.not. ieee_is_finite(report%p_error)) return
        report%message = 'G6 external B/p agreement failed'
        if (report%b_error > 1.0e-3_dp) return
        if (report%p_error > report%p_tolerance) return
        report%message = 'PASS: G6 sampled B/p agreement only'
        status = 0
    end subroutine

    subroutine norm_add(acc, value)
        ! Scaled sum of squares avoids overflow/underflow in raw squares.
        type(norm_t), intent(inout) :: acc
        real(dp), intent(in) :: value
        real(dp) :: a
        a = abs(value)
        if (a == 0) return
        if (acc%scale < a) then
            acc%ssq = 1 + acc%ssq*(acc%scale/a)**2
            acc%scale = a
        else
            acc%ssq = acc%ssq + (a/acc%scale)**2
        end if
    end subroutine

    subroutine difference(a, b, delta, ok)
        real(dp), intent(in) :: a, b
        real(dp), intent(out) :: delta
        logical, intent(out) :: ok
        delta = 0
        ok = .false.
        if (b < 0 .and. a > 0) then
            if (a > huge(a) + b) return
        else if (b > 0 .and. a < 0) then
            if (a < -huge(a) + b) return
        end if
        delta = a - b
        ok = ieee_is_finite(delta)
    end subroutine

    subroutine norm_ratio(numerator, denominator, value, ok)
        type(norm_t), intent(in) :: numerator, denominator
        real(dp), intent(out) :: value
        logical, intent(out) :: ok
        real(dp) :: ratio, factor
        value = huge(1.0_dp)
        ok = .false.
        if (denominator%scale <= 0) return
        if (numerator%scale == 0) then
            value = 0
            ok = .true.
            return
        end if
        if (denominator%scale < 1) then
            if (numerator%scale > huge(value)*denominator%scale) return
        end if
        ratio = numerator%scale/denominator%scale
        factor = sqrt(numerator%ssq/denominator%ssq)
        if (factor > 1) then
            if (ratio > huge(value)/factor) return
        end if
        value = ratio*factor
        ok = ieee_is_finite(value)
    end subroutine
end module
