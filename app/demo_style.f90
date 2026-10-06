!> Presentation helpers shared by the gallery demos (docs/figure-style.md).
!> Display only: no physics, numerics or certificate logic lives here.
module kin6d_demo_style
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
    implicit none
    private
    public :: blue, vermillion, green, orange, purple, sky, black, grey
    public :: movie_fps, hold_frames, num
    ! Okabe-Ito categorical palette; one colour per role across all demos.
    real(dp), parameter :: blue(3) = [0.000_dp, 0.447_dp, 0.698_dp]
    real(dp), parameter :: vermillion(3) = [0.835_dp, 0.369_dp, 0.000_dp]
    real(dp), parameter :: green(3) = [0.000_dp, 0.620_dp, 0.451_dp]
    real(dp), parameter :: orange(3) = [0.902_dp, 0.624_dp, 0.000_dp]
    real(dp), parameter :: purple(3) = [0.800_dp, 0.475_dp, 0.655_dp]
    real(dp), parameter :: sky(3) = [0.337_dp, 0.706_dp, 0.914_dp]
    real(dp), parameter :: black(3) = [0.0_dp, 0.0_dp, 0.0_dp]
    real(dp), parameter :: grey(3) = [0.55_dp, 0.55_dp, 0.55_dp]
    !> Movie frame rate and the closing hold (about one second).
    integer, parameter :: movie_fps = 24, hold_frames = 24
contains
    !> Fixed-format number for titles: `1.23e-03`, never a denormal or a
    !> Fortran `E+00` exponent; `digits + 1` significant digits. Magnitudes
    !> below 1e-300 print as `0`.
    pure function num(x, digits) result(text)
        real(dp), intent(in) :: x
        integer, intent(in), optional :: digits
        character(:), allocatable :: text
        character(32) :: buffer, form
        integer :: d, e, mantissa_end
        d = 2
        if (present(digits)) d = max(0, min(digits, 12))
        if (ieee_is_nan(x)) then
            text = 'NaN'
            return
        end if
        if (abs(x) < 1.0e-300_dp) then
            text = '0'
            return
        end if
        write (form, '(a,i0,a,i0,a)') '(es', d + 9, '.', d, 'e3)'
        write (buffer, form) x
        buffer = adjustl(buffer)
        mantissa_end = index(buffer, 'E') - 1
        read (buffer(mantissa_end + 2:), *) e
        if (e >= -2 .and. e <= 3) then
            if (d - e <= 0) then
                write (buffer, '(i0)') nint(x)
            else
                write (form, '(a,i0,a,i0,a)') '(f', d + 8, '.', d - e, ')'
                write (buffer, form) x
            end if
            text = trim(adjustl(buffer))
            return
        end if
        write (form, '(a,i2.2)') 'e', abs(e)
        if (e < 0) form = 'e-'//form(2:)
        text = buffer(1:mantissa_end)//trim(form)
    end function num
end module kin6d_demo_style
