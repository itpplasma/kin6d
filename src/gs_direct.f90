module kin6d_gs_direct
    ! Cached sparse Cholesky factor of the Grad--Shafranov stiffness matrix
    ! (threaded MUMPS, METIS nested dissection, through fortsparse). One factor
    ! serves every solve on the mesh: fixed-source solves, the A^-1 Newton
    ! preconditioner and the majorant recovery.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortsparse, only: sparse_solver_t, sparse_factor, sparse_solve, &
        sparse_destroy, sparse_factor_stats, sparse_factor_stats_t, &
        fortsparse_status_t, status_ok, FORTSPARSE_BACKEND_MUMPS, FORTSPARSE_SYM_SPD
    implicit none
    private
    public :: gs_direct_t, gs_direct_available, gs_direct_factor, gs_direct_solve
    public :: gs_direct_free
    type :: gs_direct_t
        type(sparse_solver_t) :: solver
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
        available = .true.
    end function

    subroutine gs_direct_factor(direct, n, rowptr, colidx, val, ok)
        ! Factor the full symmetric positive definite matrix given as CSR
        ! (identical to CSC by symmetry). A failure (for example a matrix that is
        ! not positive definite) is recorded and not retried.
        type(gs_direct_t), intent(inout) :: direct
        integer, intent(in) :: n, rowptr(:), colidx(:)
        real(dp), intent(in) :: val(:)
        logical, intent(out) :: ok
        type(fortsparse_status_t) :: status
        type(sparse_factor_stats_t) :: stats
        call gs_direct_free(direct)
        direct%solver%backend_id = FORTSPARSE_BACKEND_MUMPS
        direct%solver%refine = .false.
        direct%solver%options%symmetry = FORTSPARSE_SYM_SPD
        direct%solver%options%nthreads = 0
        call sparse_factor(direct%solver, n, n, rowptr(n + 1) - 1, rowptr, &
                           colidx, val, status)
        ok = status_ok(status)
        direct%failed = .not. ok
        direct%factored = ok
        direct%n = n
        if (.not. ok) return
        call sparse_factor_stats(direct%solver, stats, status)
        if (status_ok(status)) then
            direct%factor_s = max(0.0_dp, stats%factor_s)
            direct%analysis_s = max(0.0_dp, stats%analysis_s)
            direct%factor_nnz = int(min(int(huge(1), 8), stats%factor_nnz))
        end if
    end subroutine

    subroutine gs_direct_solve(direct, rhs, x, ok)
        type(gs_direct_t), intent(inout) :: direct
        real(dp), intent(in) :: rhs(:)
        real(dp), intent(out) :: x(:)
        logical, intent(out) :: ok
        type(fortsparse_status_t) :: status
        ok = .false.
        x = 0
        if (.not. direct%factored) return
        call sparse_solve(direct%solver, rhs, x, status)
        ok = status_ok(status)
    end subroutine

    subroutine gs_direct_free(direct)
        type(gs_direct_t), intent(inout) :: direct
        if (direct%factored .or. direct%failed) call sparse_destroy(direct%solver)
        direct%factored = .false.
        direct%failed = .false.
    end subroutine
end module
