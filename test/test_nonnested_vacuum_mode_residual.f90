program test_nonnested_vacuum_mode_residual
    ! The n = 1 helical-coil vacuum solve on the refined annulus (P3,
    ! nr = 8 and 16 by default, or the rings named on the command line,
    ! e.g. 32). Its right-hand side is zero and the boundary
    ! trace is the only forcing, so the reported residual must be relative
    ! to the lifted right-hand side. Before that normalisation the absolute
    ! residual of the converged direct solve (1.2e-9 at nr = 16, 3.84e-9 at
    ! nr = 32 in the Newton loop) failed the 1e-9 gate.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use nn_vacuum_fixture, only: r0, ftor, field, active_mode, coil_mu0, &
        make_coil, gauge_potential, annulus_mesh
    use kin6d_mhd_mesh, only: flux_mesh_t
    use kin6d_equilibrium_nonnested, only: nn_report_t, nn_build, nn_vacuum_mode
    implicit none
    integer, parameter :: k = 3
    type(flux_mesh_t) :: fm
    type(nn_report_t) :: report
    integer, allocatable :: rings(:)
    character(16) :: arg
    integer :: st, nr, i
    if (command_argument_count() == 0) then
        rings = [8, 16]
    else
        allocate (rings(command_argument_count()))
        do i = 1, size(rings)
            call get_command_argument(i, arg)
            read (arg, *) rings(i)
        end do
    end if
    call make_coil()
    do i = 1, size(rings)
        nr = rings(i)
        call annulus_mesh(nr, k, fm)
        call nn_build(fm, k, 1, [r0, 0.0_dp], 1.0_dp, ftor, coil_mu0, field, st)
        if (st /= 0) error stop 'build'
        active_mode = 1
        call nn_vacuum_mode(field, 1, gauge_potential, report, st)
        print '(a,i0,a,i0,a,es10.3)', 'nr ', nr, ' status ', st, &
            ' residual ', report%residual
        if (st /= 0) error stop 'vacuum mode'
        if (report%residual > 1.0e-12_dp) error stop 'relative residual'
    end do
end program
