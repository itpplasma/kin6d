!> Title numbers: hand-written expected strings are the oracle.
program test_demo_style
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_demo_style, only: num
    implicit none
    integer :: failures = 0
    call check(num(2.0e-323_dp), '0')
    call check(num(0.0_dp), '0')
    call check(num(1.234e-3_dp), '1.23e-03')
    call check(num(-4.5e-12_dp), '-4.50e-12')
    call check(num(6.02e23_dp, 3), '6.020e23')
    call check(num(0.5_dp), '0.500')
    call check(num(0.0123_dp), '0.0123')
    call check(num(3.14159_dp, 3), '3.142')
    call check(num(1234.0_dp), '1234')
    call check(num(12345.0_dp), '1.23e04')
    call check(num(1.0e-14_dp, 0), '1e-14')
    if (failures > 0) error stop 'demo number formatting failed'
    print '(a)', 'demo number formatting passed'
contains
    subroutine check(got, expected)
        character(*), intent(in) :: got, expected
        if (got /= expected) then
            print '(5a)', 'got "', got, '" expected "', expected, '"'
            failures = failures + 1
        end if
    end subroutine check
end program test_demo_style
