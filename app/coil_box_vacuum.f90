program coil_box_vacuum
    !! Vacuum field of MAKEGRID coils in a domain bounded by a star-shaped
    !! Fourier wall with P field periods, from the tangential coil vector
    !! potential on the wall (mapped NN1), checked against direct
    !! Biot-Savart in the L2 norm over the domain.
    !! Usage: coil_box_vacuum input.nml   (namelist /box/)
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use coil_box_geometry
    use coil_box_setup, only: box_build, poincare
    use milca_coils_3d, only: coil_mu0
    use milca_equilibrium_nonnested, only: free_b_field_t
    use milca_nonnested_mapped, only: nn_map_t, nn_mapped_report_t, &
        nn_mapped_cartesian_data, nn_mapped_solve, nn_mapped_cartesian_error
    implicit none
    character(len=512) :: coils_file = '', wall_file = '', arg
    real(dp) :: extcur(64) = 0, tolerance = 1.0e-10_dp
    integer :: period = 1, nring = 8, ncore = 8, nthmax = 256, order = 2, modes = 4
    integer :: map_modes = 8, nphi = 0, nphi_error = 0
    integer :: trace_lines = 0, trace_periods = 100
    real(dp) :: trace_u0 = 0.5_dp, trace_u1 = 0.95_dp
    character(len=512) :: trace_file = 'poincare.txt'
    ! A test gate: stop with an error if the L2 field error exceeds it.
    real(dp) :: max_error = huge(1.0_dp)
    namelist /box/ coils_file, wall_file, extcur, period, nring, ncore, nthmax, &
        order, modes, map_modes, nphi, nphi_error, tolerance, trace_lines, &
        trace_periods, trace_u0, trace_u1, trace_file, max_error
    type(free_b_field_t), target :: field
    type(nn_map_t), target :: map
    type(nn_mapped_report_t) :: report
    real(dp), allocatable :: data(:)
    real(dp) :: err, eh, ee
    integer :: unit, st, c0, c1, rate

    call get_command_argument(1, arg)
    open (newunit=unit, file=trim(arg), status='old', action='read')
    read (unit, nml=box)
    close (unit)
    if (nphi == 0) nphi = 4*modes + 4
    if (nphi_error == 0) nphi_error = 2*nphi + 1
    call system_clock(c0, rate)

    call box_build(trim(coils_file), trim(wall_file), extcur, period, nring, ncore, &
        nthmax, order, modes, map_modes, nphi, field, map)
    call nn_mapped_cartesian_data(field, map, nphi, box_potential, data, st)
    if (st /= 0) error stop 'data (map Jacobian not positive?)'
    call nn_mapped_solve(field, map, nphi, data, tolerance, report, st)
    if (st /= 0) error stop 'solve'
    call nn_mapped_cartesian_error(field, map, nphi_error, box_field, err, eh, ee, st)
    if (st /= 0) error stop 'error norm'
    call system_clock(c1)
    write (*, '(a,i0,a,i0,a,i0,a,i0,a,i0)') 'period=', period, ' rings=', nring, &
        ' order=', order, ' modes=', modes, ' ndof=', report%ndof
    write (*, '(a,es10.3,a,es10.3,a,i0,a,es9.2,a,es9.2,a,f9.2)') &
        'field error / |B| (L2)= ', err/sqrt(2*coil_mu0*ee), &
        ' energy rel= ', abs(eh - ee)/ee, ' fgmres= ', report%iterations, &
        ' residual= ', report%residual, ' div= ', report%divergence, &
        ' s= ', real(c1 - c0, dp)/rate
    if (err/sqrt(2*coil_mu0*ee) > max_error) error stop 'field error above max_error'
    if (trace_lines > 0) call poincare(field, period, trace_lines, trace_periods, &
        trace_u0, trace_u1, trim(trace_file))
end program coil_box_vacuum
