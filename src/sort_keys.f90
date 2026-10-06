!> Stable merge-sort permutation for 64-bit integer keys.
!> Used to group sparse many-body contributions by residual determinant or
!> by reduced-density-matrix coordinate. No floating-point data are touched.
module kin6d_sort_keys
    use, intrinsic :: iso_fortran_env, only: i64 => int64
    implicit none
    private
    public :: sort_key_order

contains

    !> order is the stable permutation with keys(order) nondecreasing.
    subroutine sort_key_order(keys, order)
        integer(i64), intent(in) :: keys(:)
        integer, intent(out) :: order(:)
        integer, allocatable :: scratch(:)
        integer :: n, width, lo, mid, hi, i

        n = size(keys)
        do i = 1, n
            order(i) = i
        end do
        if (n < 2) return
        allocate (scratch(n))
        width = 1
        do while (width < n)
            lo = 1
            do while (lo <= n - width)
                mid = lo + width - 1
                hi = min(lo + 2*width - 1, n)
                call merge_runs(keys, order, scratch, lo, mid, hi)
                lo = lo + 2*width
            end do
            width = 2*width
        end do
    end subroutine sort_key_order

    subroutine merge_runs(keys, order, scratch, lo, mid, hi)
        integer(i64), intent(in) :: keys(:)
        integer, intent(inout) :: order(:), scratch(:)
        integer, intent(in) :: lo, mid, hi
        integer :: i, j, k

        i = lo; j = mid + 1; k = lo
        do while (i <= mid .and. j <= hi)
            if (keys(order(j)) < keys(order(i))) then
                scratch(k) = order(j); j = j + 1
            else
                scratch(k) = order(i); i = i + 1
            end if
            k = k + 1
        end do
        do while (i <= mid)
            scratch(k) = order(i); i = i + 1; k = k + 1
        end do
        do while (j <= hi)
            scratch(k) = order(j); j = j + 1; k = k + 1
        end do
        order(lo:hi) = scratch(lo:hi)
    end subroutine merge_runs
end module kin6d_sort_keys
