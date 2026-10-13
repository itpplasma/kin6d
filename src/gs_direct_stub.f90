module kin6d_gs_direct
    ! Build without MUMPS (KIN6D_ENABLE_MUMPS=OFF): no direct factor is available
    ! and every caller keeps the iterative path.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    private
    public :: gs_direct_t, gs_direct_available, gs_direct_factor, gs_direct_solve
    public :: gs_direct_free, gs_direct_solve_block, GS_DIRECT_SPD, GS_DIRECT_INDEFINITE, &
        GS_DIRECT_ORDER_AMD, GS_DIRECT_ORDER_ND
    integer, parameter :: GS_DIRECT_ORDER_AMD = 1, GS_DIRECT_ORDER_ND = 2
    integer, parameter :: GS_DIRECT_SPD = 1, GS_DIRECT_INDEFINITE = 2
    type :: gs_direct_t
        logical :: factored = .false.
        logical :: failed = .false.
        ! Dirichlet mask the factor was built for; a change refactors.
        logical, allocatable :: mask(:)
        integer :: n = 0
        integer :: factor_nnz = 0
        real(dp) :: factor_s = 0, analysis_s = 0
    end type
contains
    logical function gs_direct_available() result(available)
        available = .false.
    end function

    subroutine gs_direct_factor(direct, n, rowptr, colidx, val, ok, symmetry, ordering)
        type(gs_direct_t), intent(inout) :: direct
        integer, intent(in) :: n, rowptr(:), colidx(:)
        real(dp), intent(in) :: val(:)
        logical, intent(out) :: ok
        integer, intent(in), optional :: symmetry, ordering
        ok = .false.
        direct%failed = .true.
    end subroutine

    subroutine gs_direct_solve(direct, rhs, x, ok)
        type(gs_direct_t), intent(inout) :: direct
        real(dp), intent(in) :: rhs(:)
        real(dp), intent(out) :: x(:)
        logical, intent(out) :: ok
        ok = .false.
        x = 0
    end subroutine

    subroutine gs_direct_solve_block(direct, rhs, x, ok)
        type(gs_direct_t), intent(inout) :: direct
        real(dp), intent(in) :: rhs(:, :)
        real(dp), intent(out) :: x(:, :)
        logical, intent(out) :: ok
        ok = .false.
        x = 0
    end subroutine

    subroutine gs_direct_free(direct)
        type(gs_direct_t), intent(inout) :: direct
        direct%factored = .false.
    end subroutine
end module
