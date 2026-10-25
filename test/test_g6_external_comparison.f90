module g6_external_test_fixture
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, &
        ieee_positive_inf
    implicit none
    private
    public :: set_fixture, sample_fixture
    real(dp) :: fixture_b(3) = [2.0_dp, 0.0_dp, 0.0_dp]
    real(dp) :: fixture_p(3) = 0
    logical :: invalid_b = .false., invalid_p = .false.
contains
    subroutine set_fixture(b, pressure, nonfinite_b, nonfinite_p)
        real(dp), intent(in) :: b(3), pressure(3)
        logical, intent(in), optional :: nonfinite_b, nonfinite_p
        fixture_b = b
        fixture_p = pressure
        invalid_b = .false.
        invalid_p = .false.
        if (present(nonfinite_b)) invalid_b = nonfinite_b
        if (present(nonfinite_p)) invalid_p = nonfinite_p
    end subroutine set_fixture

    subroutine sample_fixture(xyz, b_cartesian, mu0_p, status)
        real(dp), intent(in) :: xyz(3)
        real(dp), intent(out) :: b_cartesian(3), mu0_p
        integer, intent(out) :: status
        integer :: row
        b_cartesian = fixture_b
        mu0_p = 0
        status = 1
        if (xyz(1) < 1 .or. xyz(1) > 3) return
        row = nint(xyz(1))
        mu0_p = fixture_p(row)
        if (invalid_b) b_cartesian(1) = ieee_value(0.0_dp, ieee_quiet_nan)
        if (invalid_p) mu0_p = ieee_value(0.0_dp, ieee_positive_inf)
        status = 0
    end subroutine sample_fixture
end module g6_external_test_fixture

program test_g6_external_comparison
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, &
        ieee_positive_inf, ieee_is_nan
    use g6_external_comparison, only: g6_comparison_t, compare_external_samples
    use g6_external_test_fixture, only: set_fixture, sample_fixture
    implicit none
    character(*), parameter :: path = 'g6_external_reader_fixture.dat'
    real(dp), parameter :: b(3) = [2.0_dp, 0.0_dp, 0.0_dp]
    real(dp), parameter :: zero(3) = 0
    type(g6_comparison_t) :: report
    real(dp) :: refs(3), sampled(3), delta(3), expected, atol, invalid(4)
    real(dp) :: vector_b(3), extreme_b(3), scaled_refs(3), pressure_scales(2)
    integer :: status, i, checks
    checks = 0

    call write_fixture(zero, b)
    call set_fixture(b, zero)
    call compare_external_samples(path, sample_fixture, report, status, 1.0e-6_dp)
    call passed('zero pressure with explicit tolerance')
    call require(report%samples == 3, 'three samples retained')
    call require(trim(report%pressure_kind) == 'absolute RMS', &
                 'zero-pressure comparison kind')
    call require(report%p_error == 0, 'zero-pressure exact agreement')
    call require(report%p_tolerance == 1.0e-6_dp, 'absolute tolerance retained')

    call set_fixture(b, [2.0e-6_dp, 2.0e-6_dp, 2.0e-6_dp])
    call compare_external_samples(path, sample_fixture, report, status, 1.0e-6_dp)
    call rejected('zero-pressure perturbation exceeds tolerance')
    call set_fixture(b, zero)
    call compare_external_samples(path, sample_fixture, report, status)
    call rejected('zero pressure requires an explicit tolerance')
    invalid = [0.0_dp, -1.0_dp, ieee_value(0.0_dp, ieee_quiet_nan), &
               ieee_value(0.0_dp, ieee_positive_inf)]
    do i = 1, size(invalid)
        call compare_external_samples(path, sample_fixture, report, status, &
                                      invalid(i))
        call rejected('invalid zero-pressure tolerance')
    end do

    atol = 0.125_dp
    call set_fixture(b, [atol, atol, atol])
    call compare_external_samples(path, sample_fixture, report, status, atol)
    call passed('absolute RMS at the exact tolerance boundary')
    call require(report%p_error == atol, 'independent absolute RMS boundary')
    call set_fixture(b, [atol, atol, atol]*(1 + 1.0e-6_dp))
    call compare_external_samples(path, sample_fixture, report, status, atol)
    call rejected('absolute RMS above tolerance boundary')

    refs = [2.0_dp, 4.0_dp, 4.0_dp]
    delta = [0.001_dp, 0.002_dp, 0.002_dp]
    sampled = refs + delta
    expected = sqrt(sum(delta**2)/sum(refs**2))
    call write_fixture(refs, b)
    call set_fixture(b, sampled)
    call compare_external_samples(path, sample_fixture, report, status)
    call passed('unchanged nonzero relative-pressure convention')
    call require(trim(report%pressure_kind) == 'relative RMS', &
                 'nonzero-pressure comparison kind')
    call require(abs(report%p_error - expected) < 1.0e-15_dp, &
                 'independent nonuniform relative-pressure RMS')
    call require(abs(report%p_error - 5.0e-4_dp) < 1.0e-15_dp, &
                 'nonzero reference-pressure denominator')
    call require(report%p_tolerance == 1.0e-3_dp, &
                 'relative tolerance remains unchanged')
    call compare_external_samples(path, sample_fixture, report, status, -1.0_dp)
    call passed('irrelevant absolute tolerance ignored for nonzero pressure')
    call set_fixture(b, refs*(1 + 2.0e-3_dp))
    call compare_external_samples(path, sample_fixture, report, status, 1.0e9_dp)
    call rejected('absolute tolerance cannot relax nonzero relative gate')

    pressure_scales = [1.0e-200_dp, 1.0e200_dp]
    do i = 1, size(pressure_scales)
        scaled_refs = refs*pressure_scales(i)
        call write_fixture(scaled_refs, b)
        call set_fixture(b, scaled_refs + delta*pressure_scales(i))
        call compare_external_samples(path, sample_fixture, report, status)
        call passed('relative pressure survives extreme finite common scaling')
        call require(trim(report%pressure_kind) == 'relative RMS', &
                     'small nonzero pressure is not treated as zero')
        call require(abs(report%p_error - expected) < 1.0e-15_dp, &
                     'relative RMS is homogeneous under finite scaling')
    end do

    vector_b = [2.0_dp, 3.0_dp, 6.0_dp]
    call write_fixture(refs, vector_b)
    call set_fixture(vector_b + [0.006_dp, 0.0_dp, 0.0_dp], refs)
    call compare_external_samples(path, sample_fixture, report, status)
    call passed('total-vector magnetic denominator')
    expected = 0.006_dp/sqrt(sum(vector_b**2))
    call require(abs(report%b_error - expected) < 1.0e-15_dp, &
                 'independent total-vector magnetic RMS')
    call set_fixture(vector_b + [0.008_dp, 0.0_dp, 0.0_dp], refs)
    call compare_external_samples(path, sample_fixture, report, status)
    call rejected('magnetic relative RMS exceeds unchanged tolerance')

    ! Parsed physical columns carry an independent external-only diagnostic.
    vector_b = [2.0_dp, -3.0_dp, 4.0_dp]
    call set_fixture(vector_b, refs)
    call write_fixture(refs, vector_b, 'force')
    call compare_external_samples(path, sample_fixture, report, status)
    call passed('nonzero external force does not invent an admission gate')
    call require(report%force_status == 0, 'parsed physical force available')
    call require(abs(report%force_rms - 3) < 1.0e-14_dp, &
                 'parsed current/field/gradient order and units')
    call write_fixture(refs, vector_b, 'force_overflow')
    call compare_external_samples(path, sample_fixture, report, status)
    call passed('force arithmetic unavailability preserves B/p admission')
    call require(report%force_status /= 0, 'parsed force overflow unavailable')
    call require(ieee_is_nan(report%force_rms), 'unavailable force remains NaN')

    call set_fixture(b, refs)
    call write_fixture(refs, b, 'malformed')
    call compare_external_samples(path, sample_fixture, report, status)
    call rejected('malformed sample row')
    call write_fixture(refs, b, 'nonfinite')
    call compare_external_samples(path, sample_fixture, report, status)
    call rejected('nonfinite reference sample')
    call write_fixture(refs, b, 'outdomain')
    call compare_external_samples(path, sample_fixture, report, status)
    call rejected('sampler rejects an out-of-domain point')
    call write_fixture(refs, b, 'few')
    call compare_external_samples(path, sample_fixture, report, status)
    call rejected('fewer than three samples')
    call write_fixture(refs, b, 'extra')
    call compare_external_samples(path, sample_fixture, report, status)
    call rejected('extra numeric sample column')
    call write_fixture(refs, b, 'short')
    call compare_external_samples(path, sample_fixture, report, status)
    call rejected('missing numeric sample column')

    call write_fixture(refs, b)
    call set_fixture(b, refs, nonfinite_b=.true.)
    call compare_external_samples(path, sample_fixture, report, status)
    call rejected('nonfinite sampler field')
    call set_fixture(b, refs, nonfinite_p=.true.)
    call compare_external_samples(path, sample_fixture, report, status)
    call rejected('nonfinite sampler pressure')
    call write_fixture(refs, zero)
    call set_fixture(zero, refs)
    call compare_external_samples(path, sample_fixture, report, status)
    call rejected('zero magnetic reference norm')
    extreme_b = [huge(1.0_dp), 0.0_dp, 0.0_dp]
    call write_fixture(refs, extreme_b)
    call set_fixture(-extreme_b, refs)
    call compare_external_samples(path, sample_fixture, report, status)
    call rejected('overflowing finite magnetic difference is rejected')

    call delete_fixture()
    print '(a,i0,a)', 'PASS: ', checks, ' G6 external-reader checks'
contains
    subroutine require(condition, label)
        logical, intent(in) :: condition
        character(*), intent(in) :: label
        checks = checks + 1
        if (condition) return
        write (error_unit, '(2a)') 'FAIL: ', label
        call delete_fixture()
        error stop 1
    end subroutine require

    subroutine passed(label)
        character(*), intent(in) :: label
        call require(status == 0, label)
    end subroutine passed

    subroutine rejected(label)
        character(*), intent(in) :: label
        call require(status /= 0, label)
        call require(len_trim(report%message) > 0, 'explicit failure message')
    end subroutine rejected

    subroutine write_fixture(pressure, field, variant)
        real(dp), intent(in) :: pressure(3), field(3)
        character(*), intent(in), optional :: variant
        real(dp) :: row(14)
        character(16) :: kind
        integer :: unit, j, rows
        kind = ''
        if (present(variant)) kind = variant
        rows = 3
        if (kind == 'few') rows = 2
        open (newunit=unit, file=path, status='replace', action='write')
        do j = 1, rows
            row = 0
            row(1) = real(j, dp)/4
            row(2) = real(j, dp)
            row(5:7) = field
            row(8) = pressure(j)
            if (kind == 'force') then
                row(9:11) = [5.0_dp, 6.0_dp, -7.0_dp]
                row(12:14) = [2.0_dp, -32.0_dp, -29.0_dp]
            else if (kind == 'force_overflow') then
                row(10) = huge(1.0_dp)
            end if
            if (j == 2) then
                select case (kind)
                case ('malformed')
                    write (unit, '(a)') 'not a fourteen-column numeric row'
                    cycle
                case ('nonfinite')
                    row(8) = ieee_value(0.0_dp, ieee_quiet_nan)
                case ('outdomain')
                    row(2) = 99
                case ('extra')
                    write (unit, '(*(es25.16e3,1x))') row, 1.0_dp
                    cycle
                case ('short')
                    write (unit, '(*(es25.16e3,1x))') row(:13)
                    cycle
                end select
            end if
            write (unit, '(*(es25.16e3,1x))') row
        end do
        close (unit)
    end subroutine write_fixture

    subroutine delete_fixture()
        integer :: unit, ios
        open (newunit=unit, file=path, status='old', iostat=ios)
        if (ios == 0) close (unit, status='delete')
    end subroutine delete_fixture
end program test_g6_external_comparison
