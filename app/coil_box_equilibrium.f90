program coil_box_equilibrium
    !! Finite-pressure resistive equilibrium (NN2, closure A+D) of MAKEGRID
    !! coils in the mapped coil box: heating source about a given axis,
    !! p = 0 on the wall, zero net toroidal current, continuation in
    !! chi_par/chi_perp and g_par/g_perp. Islands are not imposed; they
    !! follow from the coils and the plasma response.
    !! Usage: coil_box_equilibrium input.nml   (namelists /box/, /plasma/)
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use coil_box_geometry
    use coil_box_setup, only: box_build, poincare
    use milca_coils_3d, only: coil_mu0
    use milca_equilibrium_nonnested, only: free_b_field_t, nn_fourier_value
    use milca_equilibrium_state, only: EQ_VALUES, EQ_BR, EQ_BZ, EQ_BPHI
    use milca_nonnested_mapped, only: nn_map_t, nn_mapped_cartesian_data, nn_map_point
    use milca_nonnested_coupled, only: nc_problem_t, nc_layout_t, nc_report_t, &
        nc_layout, nc_continue, nc_solve
    use fortfem_meridian_feec_n, only: meridian_qgrad
    implicit none
    character(len=512) :: coils_file = '', wall_file = '', arg
    real(dp) :: extcur(64) = 0, tolerance = 1.0e-10_dp
    integer :: period = 1, nring = 8, ncore = 8, nthmax = 256, order = 2, modes = 4
    integer :: map_modes = 8, nphi = 0, nphi_error = 0
    integer :: trace_lines = 0, trace_periods = 100
    real(dp) :: trace_u0 = 0.5_dp, trace_u1 = 0.95_dp, max_error = huge(1.0_dp)
    character(len=512) :: trace_file = 'poincare.txt'
    namelist /box/ coils_file, wall_file, extcur, period, nring, ncore, nthmax, &
        order, modes, map_modes, nphi, nphi_error, tolerance, trace_lines, &
        trace_periods, trace_u0, trace_u1, trace_file, max_error
    real(dp) :: chi_perp = 1, chi_ratio = 1.0e4_dp, g_par = 1, g_ratio = 1.0e4_dp
    real(dp) :: decades = 1, newton_tolerance = 1.0e-10_dp
    ! A test gate on <beta>; the default accepts any value.
    real(dp) :: beta_range(2) = [-huge(1.0_dp), huge(1.0_dp)]
    ! Continuation in the source: the anisotropy continuation runs at
    ! source_peak/source_steps, then each step adds that amount and is
    ! solved by Newton from the previous state at the target anisotropy.
    integer :: source_steps = 1
    namelist /plasma/ axis_rc, axis_zs, source_peak, source_width, chi_perp, &
        chi_ratio, g_par, g_ratio, decades, newton_tolerance, beta_range, &
        source_steps
    type(free_b_field_t), target :: field
    type(nn_map_t), target :: map
    type(nc_problem_t) :: problem
    type(nc_layout_t) :: lay
    type(nc_report_t) :: report
    real(dp), allocatable :: data(:), x(:)
    real(dp) :: pmax, pmean, b2mean, volume, peak
    integer :: unit, st, c0, c1, rate, step

    call get_command_argument(1, arg)
    open (newunit=unit, file=trim(arg), status='old', action='read')
    read (unit, nml=box)
    read (unit, nml=plasma)
    close (unit)
    if (nphi == 0) nphi = 4*modes + 4
    call system_clock(c0, rate)
    call box_build(trim(coils_file), trim(wall_file), extcur, period, nring, ncore, &
        nthmax, order, modes, map_modes, nphi, field, map)
    ! The projected coil potential seeds A and fixes its wall trace.
    call nn_mapped_cartesian_data(field, map, nphi, box_potential, data, st)
    if (st /= 0) error stop 'data'
    lay = nc_layout(field)
    allocate (x(lay%n))
    x = 0
    x(lay%ia + 1:lay%ia + lay%na) = data
    problem%chi_perp = chi_perp; problem%chi_par = chi_ratio*chi_perp
    problem%g_par = g_par; problem%g_perp = g_par/g_ratio
    problem%current = 0
    problem%nphi = nphi
    problem%source_3d => box_source
    problem%map => map
    if (source_steps < 1) error stop 'source_steps'
    peak = source_peak
    source_peak = peak/source_steps
    call nc_continue(field, problem, x, report, st, tol=newton_tolerance, &
        decades=decades)
    write (*, '(a,i0,a,i0,a,i0,a,es9.2,a,a)') 'source step 1 Newton stages ', &
        report%stages, ' steps ', report%total_steps, ' status ', st, &
        ' divergence ', report%divergence, ' ', trim(report%message)
    write (*, '(a,*(es9.1))') 'last stage backward errors ', &
        report%history(:report%newton_steps)
    if (st /= 0) error stop 'coupled solve'
    do step = 2, source_steps
        source_peak = peak*step/source_steps
        call nc_solve(field, problem, x, report, st, tol=newton_tolerance)
        write (*, '(a,i0,a,i0,a,i0,a,es9.2)') 'source step ', step, &
            ' Newton steps ', report%newton_steps, ' status ', st, &
            ' divergence ', report%divergence
        if (st /= 0) error stop 'coupled solve (source step)'
    end do
    call moments(field, map, nphi, pmax, pmean, b2mean, volume)
    call system_clock(c1)
    write (*, '(a,es11.3,a,es11.3,a,es11.3,a,es11.3,a,es11.3,a,f9.1)') &
        'p_max (Pa)= ', pmax, ' <p>= ', pmean, ' <beta>= ', &
        2*coil_mu0*pmean/b2mean, ' volume (m^3)= ', volume, ' V= ', report%voltage, &
        ' s= ', real(c1 - c0, dp)/rate
    if (report%divergence > 1.0e-12_dp) error stop 'divergence'
    if (2*coil_mu0*pmean/b2mean < beta_range(1) .or. &
        2*coil_mu0*pmean/b2mean > beta_range(2)) error stop '<beta> outside beta_range'
    if (trace_lines > 0) call poincare(field, period, trace_lines, trace_periods, &
        trace_u0, trace_u1, trim(trace_file))
contains
    subroutine moments(field, map, nphi, pmax, pmean, b2mean, volume)
        ! Physical-volume means of p and |B|^2 and the largest p at the
        ! quadrature points of nphi planes of one period.
        type(free_b_field_t), intent(in) :: field
        type(nn_map_t), intent(in) :: map
        integer, intent(in) :: nphi
        real(dp), intent(out) :: pmax, pmean, b2mean, volume
        real(dp) :: xq(2), wq, val(field%fs%nl0), grad(2, field%fs%nl0), pos(3)
        real(dp) :: p1(3, 3), p2(3, 3), jw, ph, p, values(EQ_VALUES, 1), pts(3, 1), bref(3)
        integer :: t, q, j, st, sts(1)
        pmax = 0; pmean = 0; b2mean = 0; volume = 0
        do t = 1, field%fs%mm%topology%n_triangles
            do q = 1, field%fs%nq
                call meridian_qgrad(field%fs, t, q, xq, wq, val, grad, st)
                do j = 1, nphi
                    ph = 2*acos(-1.0_dp)*(j - 1)/(nphi*field%period)
                    call nn_map_point(map, field%fs%sp(0)%dof(:, t), val, grad, xq, &
                        ph, pos, p1, p2, jw, st)
                    p = nn_fourier_value(field%p, field%fs%n0, &
                        field%fs%sp(0)%dof(:, t), val, field%period*ph)
                    pts(:, 1) = [xq, ph]
                    call field%sample(pts, values, sts)
                    bref = [values(EQ_BR, 1), values(EQ_BPHI, 1), values(EQ_BZ, 1)]
                    pmax = max(pmax, p)
                    pmean = pmean + wq*jw*p/nphi
                    b2mean = b2mean + wq*jw*sum(matmul(p2, bref)**2)/nphi
                    volume = volume + wq*jw/nphi
                end do
            end do
        end do
        pmean = pmean/volume; b2mean = b2mean/volume
        volume = 2*acos(-1.0_dp)*volume
    end subroutine
end program coil_box_equilibrium
