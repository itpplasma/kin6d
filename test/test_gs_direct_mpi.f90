program test_gs_direct_mpi
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use mpi, only: MPI_Initialized, MPI_Finalized, MPI_Init_thread, MPI_Finalize, &
        MPI_THREAD_SERIALIZED, MPI_SUCCESS
    use kin6d_gs_direct, only: gs_direct_t, gs_direct_factor, gs_direct_solve, &
        gs_direct_free
    implicit none
    type(gs_direct_t) :: solver
    character(16) :: mode
    integer :: ierr, provided, pass
    logical :: ok, initialized, finalized
    real(dp) :: x(2)

    call get_command_argument(1, mode)
    if (mode /= 'automatic') then
        call MPI_Init_thread(MPI_THREAD_SERIALIZED, provided, ierr)
        if (ierr /= MPI_SUCCESS) error stop 'test MPI initialization'
        if (provided < MPI_THREAD_SERIALIZED) error stop 'test MPI thread support'
    end if
    if (mode == 'finalized') then
        call MPI_Finalize(ierr)
        if (ierr /= MPI_SUCCESS) error stop 'test MPI finalization'
        call factor()
        if (ok .or. .not. solver%failed) error stop 'accepted finalized MPI runtime'
        call gs_direct_free(solver)
        stop
    end if

    ! Repeated factors preserve the runtime and the independent exact solution.
    do pass = 1, 2
        call factor()
        if (.not. ok) error stop 'MUMPS factor failed'
        call MPI_Initialized(initialized, ierr)
        if (ierr /= MPI_SUCCESS) error stop 'MPI initialization query'
        if (.not. initialized) error stop 'runtime not initialized'
        call gs_direct_solve(solver, [0.0_dp, 3.0_dp], x, ok)
        if (.not. ok) error stop 'MUMPS solve failed'
        if (maxval(abs(x - [1.0_dp, 2.0_dp])) > 1.0e-12_dp) &
            error stop 'MUMPS solution differs from exact solution'
        call gs_direct_free(solver)
        call MPI_Finalized(finalized, ierr)
        if (ierr /= MPI_SUCCESS) error stop 'MPI finalization query'
        if (finalized) error stop 'solver finalized application runtime'
    end do
    call MPI_Finalize(ierr)
    if (ierr /= MPI_SUCCESS) error stop 'application MPI finalization'
contains
    subroutine factor()
        call gs_direct_factor(solver, 2, [1, 3, 5], [1, 2, 1, 2], &
            [2.0_dp, -1.0_dp, -1.0_dp, 2.0_dp], ok)
    end subroutine
end program
