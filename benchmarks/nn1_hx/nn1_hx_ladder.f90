program nn1_hx_ladder
    ! NN1 per-mode vacuum solves (nn_vacuum_mode, the helical-coil annulus
    ! of test_nonnested_vacuum) with the native/MUMPS Cholesky and with
    ! HX-preconditioned CG (FortFEM fortfem_hcurl_auxiliary_space).
    ! usage: nn1_hx_ladder ORDER NR [GEO]
    ! CSV per mode: order,nr,n,threads,ndof,direct_s,hx_s,hx_its,
    !   hx_refinements,rel_diff_B,direct_residual,hx_residual,
    !   direct_linear_s,hx_linear_s (the *_s columns include the load)
    use, intrinsic :: iso_fortran_env, only: dp => real64
    !$ use omp_lib, only: omp_get_max_threads
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
    integer :: order, nr, geo, n, st, threads
    integer(8) :: c0, c1, rate
    real(dp) :: td, th
    complex(dp), allocatable :: bd(:)
    character(*), parameter :: fmt = &
        '(i0,",",i0,",",i0,",",i0,",",i0,2(",",es11.4),2(",",i0),5(",",es11.4))'

    call get_command_argument(1, arg); read (arg, *) order
    call get_command_argument(2, arg); read (arg, *) nr
    geo = 2
    if (command_argument_count() >= 3) then
        call get_command_argument(3, arg); read (arg, *) geo
    end if
    threads = 1
    !$ threads = omp_get_max_threads()
    call make_coil()
    call annulus_mesh(nr, geo, fm)
    call nn_build(fm, order, nmodes, [r0, 0.0_dp], 1.0_dp, ftor, coil_mu0, field, st)
    if (st /= 0) error stop 'nn_build'
    do n = 1, nmodes
        active_mode = n
        field%linear_solver = NN_SOLVER_DIRECT
        call system_clock(c0, rate)
        call nn_vacuum_mode(field, n, gauge_potential, rd, st)
        call system_clock(c1)
        if (st /= 0) print '(a,i0,a,es11.4)', '# direct mode solve failed, n=', n, &
            ' residual ', rd%residual
        td = real(c1 - c0, dp)/real(rate, dp)
        bd = field%b(:, n)
        field%linear_solver = NN_SOLVER_HX
        call system_clock(c0)
        call nn_vacuum_mode(field, n, gauge_potential, rh, st)
        call system_clock(c1)
        if (st /= 0) print '(a,i0,a,es11.4)', '# hx mode solve failed, n=', n, &
            ' residual ', rh%residual
        th = real(c1 - c0, dp)/real(rate, dp)
        write (*, fmt) &
            order, nr, n, threads, field%fs%nn, td, th, rh%linear_iterations, &
            rh%newton_steps, norm2(abs(field%b(:, n) - bd))/norm2(abs(bd)), &
            rd%residual, rh%residual, rd%linear_seconds, rh%linear_seconds
    end do
end program nn1_hx_ladder
