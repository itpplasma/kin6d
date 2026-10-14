program test_response_qualification
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use kin6d_response_qualification, only: response_record_t, qualification_t, &
                                            qualify_response, ladder_gate
    implicit none
    type(response_record_t) :: rows(3, 8)
    type(qualification_t) :: result
    logical :: ok
    real(dp) :: d(2)
    integer :: i, axis
    do axis = 1, 8
        do i = 1, 3
            rows(i, axis)%available = .true.
            rows(i, axis)%coupled = .true.
            rows(i, axis)%calibrated = .true.
            rows(i, axis)%full_pressure = .true.
            rows(i, axis)%all_rationals = .true.
            rows(i, axis)%residual = 1.0e-9_dp
            rows(i, axis)%balance = 1.0e-7_dp
            rows(i, axis)%extent = [real(i, dp), real(i, dp)]
            rows(i, axis)%b_vac = (1.0_dp, 0.0_dp)
            rows(i, axis)%b = cmplx(1 + 0.01_dp/2**(i - 1), 0.2_dp, dp)
        end do
    end do
    call qualify_response(rows, rows(3, 1)%b, 2.0_dp, .true., result)
    call check(result%passed .and. result%verdict == 1, 'exclusive GPEC verdict')
    call qualify_response(rows, rows(3, 1)%b, abs(rows(3, 1)%b), .true., result)
    call check(result%passed .and. result%verdict == 0, 'both codes no verdict')
    call qualify_response(rows, 2*rows(3, 1)%b, 2.0_dp, .true., result)
    call check(result%passed .and. result%verdict == 0, 'neither code no verdict')
    call qualify_response(rows, 2*rows(3, 1)%b, abs(rows(3, 1)%b), .true., result)
    call check(result%passed .and. result%verdict == 2, 'exclusive MARS verdict')
    call qualify_response(rows, rows(3, 1)%b, 2.0_dp, .false., result)
    call check(result%passed .and. result%verdict == 0, 'unadmitted MARS no verdict')
    rows(3, 1)%b = -rows(3, 1)%b
    call ladder_gate(rows(:, 1), 1, ok, d)
    call check(.not. ok, 'complex phase reversal fails despite same amplitude')
    rows(3, 1)%b = -rows(3, 1)%b
    rows(2, 1)%coupled = .false.
    call qualify_response(rows, rows(3, 1)%b, 2.0_dp, .true., result)
    call check(.not. result%passed .and. result%verdict == 0, 'single ring no verdict')
    rows(2, 1)%coupled = .true.
    rows(2, 1)%residual = 1.1e-8_dp
    call ladder_gate(rows(:, 1), 1, ok, d)
    call check(.not. ok, 'true residual threshold')
    rows(2, 1)%residual = -1.0_dp
    call ladder_gate(rows(:, 1), 1, ok, d)
    call check(.not. ok, 'negative residual rejected')
    rows(2, 1)%residual = 1.0e-9_dp
    rows(:, 6)%extent(1) = 1.0_dp
    rows(:, 6)%extent(2) = 1.0_dp
    call ladder_gate(rows(:, 6), 6, ok, d)
    call check(.not. ok, 'identical capped extents fail')
    rows(:, 6)%whole_edge = .true.
    call ladder_gate(rows(:, 6), 6, ok, d)
    call check(ok, 'whole edge strip accepts capped extents')
    rows(3, 8)%b = (1.05_dp, 0.2_dp)
    call ladder_gate(rows(:, 8), 8, ok, d)
    call check(.not. ok, 'ideal-limit final change requires one percent')
    rows(3, 8)%b = cmplx(ieee_value(0.0_dp, ieee_quiet_nan), 0.0_dp, dp)
    call ladder_gate(rows(:, 8), 8, ok, d)
    call check(.not. ok, 'nonfinite result rejected')
    rows(2, 2)%available = .false.
    call ladder_gate(rows(:, 2), 2, ok, d)
    call check(.not. ok, 'missing ladder member rejected')
    ! Explicit near-zero floor; small absolute noise must not divide by b=0.
    rows(:, 3)%b = (0.0_dp, 0.0_dp)
    rows(3, 3)%b = (4.0e-6_dp, 0.0_dp)
    call ladder_gate(rows(:, 3), 3, ok, d)
    call check(ok .and. abs(d(2) - 0.004_dp) < 1.0e-14_dp, 'vacuum scale floor')
    print '(a)', 'frozen TC24 qualification rules: PASS'
contains
    subroutine check(condition, message)
        logical, intent(in) :: condition
        character(*), intent(in) :: message
        if (.not. condition) then
            print '(a)', 'FAIL: '//message
            error stop 1
        end if
    end subroutine check
end program test_response_qualification
