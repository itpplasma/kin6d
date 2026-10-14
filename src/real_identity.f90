module kin6d_real_identity
    ! Value identity for cache keys and exact model/representation choices.
    ! Approximate agreement must not reuse a cache for different inputs or
    ! declare a piecewise-linear profile differentiable at a genuine kink.
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
    implicit none
    private
    public :: same_real_value
contains
    pure elemental logical function same_real_value(a, b) result(same)
        real(dp), intent(in) :: a, b
        integer(int64) :: ia, ib

        same = .false.
        if (ieee_is_nan(a) .or. ieee_is_nan(b)) return
        ia = transfer(a, 0_int64)
        ib = transfer(b, 0_int64)
        ! IEEE signed zeros denote the same value; all other non-NaN values
        ! have a unique binary64 representation (including infinities).
        if (iand(ia, huge(0_int64)) == 0_int64) ia = 0_int64
        if (iand(ib, huge(0_int64)) == 0_int64) ib = 0_int64
        same = ia == ib
    end function same_real_value
end module kin6d_real_identity
