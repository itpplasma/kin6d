program test_nonnested_vacuum_mode_residual
    ! The n = 1 helical-coil vacuum solve on the refined annulus (P3,
    ! nr = 16). Its right-hand side is zero and the boundary trace is the
    ! only forcing, so the reported residual must be relative to the lifted
    ! right-hand side. Before that normalisation the absolute residual of
    ! the converged direct solve (1.2e-9) failed the 1e-9 gate.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use nn_vacuum_fixture, only: r0, ftor, field, active_mode, coil_mu0, &
        make_coil, gauge_potential, annulus_mesh
    use kin6d_mhd_mesh, only: flux_mesh_t
    use kin6d_equilibrium_nonnested, only: nn_report_t, nn_build, nn_vacuum_mode
    implicit none
    integer, parameter :: k = 3
    type(flux_mesh_t) :: fm
    type(nn_report_t) :: report
    integer :: st, nr
    call make_coil()
    do nr = 8, 16, 8
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
