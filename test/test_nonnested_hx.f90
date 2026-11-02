program test_nonnested_hx
    ! The n > 0 NN1 vacuum solve on the helical-coil annulus by
    ! Hiptmair-Xu preconditioned CG (field%linear_solver = NN_SOLVER_HX)
    ! against the default direct factor. Both must pass nn_vacuum_mode's
    ! unchanged 1e-9 equation-residual gate and give the same field; the CG
    ! iteration count must stay mesh-independent (usage: ORDER NR...).
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_mhd_mesh, only: flux_mesh_t
    use kin6d_equilibrium_nonnested, only: nn_report_t, nn_build, nn_vacuum_mode, &
        NN_SOLVER_DIRECT, NN_SOLVER_HX
    use kin6d_coils_3d, only: coil_mu0
    use nn_vacuum_fixture, only: field, annulus_mesh, make_coil, gauge_potential, &
        active_mode, nmodes, r0, ftor
    implicit none
    type(flux_mesh_t) :: fm
    type(nn_report_t) :: rd, rh
    character(16) :: arg
    integer :: order, nr, n, st, i
    real(dp) :: diff
    complex(dp), allocatable :: bd(:)
    if (command_argument_count() < 2) error stop 'usage: ORDER NR...'
    call get_command_argument(1, arg); read (arg, *) order
    call make_coil()
    do i = 2, command_argument_count()
        call get_command_argument(i, arg); read (arg, *) nr
        call annulus_mesh(nr, 2, fm)
        call nn_build(fm, order, nmodes, [r0, 0.0_dp], 1.0_dp, ftor, coil_mu0, &
                      field, st)
        if (st /= 0) error stop 'nn_build'
        if (field%linear_solver /= NN_SOLVER_DIRECT) error stop 'default solver'
        do n = 1, nmodes
            active_mode = n
            field%linear_solver = NN_SOLVER_DIRECT
            call nn_vacuum_mode(field, n, gauge_potential, rd, st)
            if (st /= 0) error stop 'direct vacuum mode'
            bd = field%b(:, n)
            field%linear_solver = NN_SOLVER_HX
            call nn_vacuum_mode(field, n, gauge_potential, rh, st)
            if (st /= 0) error stop 'HX vacuum mode'
            diff = norm2(abs(field%b(:, n) - bd))/norm2(abs(bd))
            print '(a,i0,a,i0,a,i0,a,i0,a,i0,a,3es10.2)', 'P', order, ' nr ', nr, &
                ' n ', n, ' HX its ', rh%linear_iterations, ' refinements ', &
                rh%newton_steps - 1, ' rel diff/residual direct/HX ', diff, &
                rd%residual, rh%residual
            if (rh%linear_iterations < 1) error stop 'HX path not taken'
            if (diff > 1.0e-8_dp) error stop 'HX field differs from direct'
            ! Two real solves per step, 15 to 20 iterations each at P2-P4.
            if (rh%linear_iterations > 60*rh%newton_steps) &
                error stop 'HX iteration count not mesh independent'
        end do
    end do
    print '(a)', 'PASS: HX-preconditioned CG matches the direct NN1 vacuum solve'
end program
