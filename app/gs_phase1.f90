program gs_phase1_main
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64, input_unit, output_unit
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use kin6d_grad_shafranov, only: gs_state_t, initialize_gs, solve_gs, &
        write_gs_mesh, read_gs_mesh, sample_gs
    use kin6d_gs_field_readback, only: sample_gs_many
    use kin6d_curved_triangle, only: periodic_curve_t, initialize_periodic_curve, &
        curve_value
    use kin6d_gs_readback, only: gs_readback_t, read_gs_equilibrium, gs_profile_f, &
        gs_readback_resolution, gs_readback_q
    use kin6d_gs_majorant, only: gs_functional_majorant
    use kin6d_gs_profile_law, only: gs_profile_law_t, read_gs_profile_law, &
        write_gs_profile_law, gs_law_sources, gs_law_primitives
    use kin6d_gs_nonlinear, only: solve_gs_cubic_profiles, gs_normalize_current
    use kin6d_result_io, only: run_t, result_file_t, runtime_git
    use fortplot_system_runtime, only: create_directory_runtime
    implicit none
    type(gs_state_t) :: state
    type(gs_profile_law_t), allocatable :: law
    type(periodic_curve_t) :: curve
    type(gs_readback_t) :: readback
    type(run_t) :: run
    type(result_file_t) :: result
    real(dp), allocatable :: samples(:, :), contour(:, :)
    real(dp) :: max_area = 0, pprime = 0, ffprime = 0, f_edge = 0
    real(dp) :: mu0 = 1.2566370614359173e-6_dp, wall_s, perimeter
    real(dp) :: nonlinear_tolerance = 1.0e-9_dp, pressure, f, pp, ff
    real(dp) :: quadrature_tolerance = 1.0e-8_dp
    real(dp) :: total_current = 0, measured_current = 0
    real(dp) :: majorant, flux_defect, balance_defect, estimator_s
    integer :: degree = 2
    integer :: boundary_nodes = 96, radial_points = 0, angular_points = 0
    integer :: unit, status, n, i, next, majorant_status, recovery_iterations
    integer(int64) :: start, finish, rate
    logical :: curved = .true., directory_ok, recovery_converged
    character(512) :: boundary_file = '', mesh_file = 'mesh.dat', profile_file = ''
    character(256) :: message
    character(512) :: argument
    character(:), allocatable :: commit, dirty
    namelist /gs_phase1/ boundary_file, mesh_file, boundary_nodes, max_area, &
        pprime, ffprime, f_edge, mu0, curved, radial_points, angular_points, degree, &
        profile_file, nonlinear_tolerance, quadrature_tolerance, total_current
    call run%start('gs_phase1')
    if (run%argument_count() > 0) then
        call run%argument(1, argument)
        select case (trim(argument))
        case ('--evaluate', '--evaluate-stream')
            call evaluate_file(trim(argument) == '--evaluate-stream')
        case ('--q')
            call evaluate_q_file()
        case ('--source')
            call evaluate_source_file()
        case default
            error stop 'expected --config, --evaluate or --q'
        end select
        stop
    end if
    status = 0; message = ''
    call run%open_config(unit)
    if (run%has_config()) read (unit, nml=gs_phase1, iostat=status, iomsg=message)
    call run%close_config(unit, 'gs_phase1', status, message)
    if (len_trim(boundary_file) == 0) error stop 'boundary_file required'
    if (degree /= 2 .and. degree /= 3) error stop 'degree must be 2 or 3'
    if (boundary_nodes < 8) error stop 'at least eight boundary nodes required'
    call system_clock(start, rate)
    open (newunit=unit, file=trim(boundary_file), status='old', action='read')
    read (unit, *) n
    if (n < 4) error stop 'boundary requires four samples without closing duplicate'
    allocate (samples(2, n), contour(2, boundary_nodes))
    do i = 1, n
        read (unit, *) samples(:, i)
    end do
    close (unit)
    call initialize_periodic_curve(curve, samples, status)
    if (status /= 0) error stop 'invalid closed curve'
    do i = 1, boundary_nodes
        call boundary(real(i - 1, dp)/boundary_nodes, contour(:, i))
    end do
    if (max_area <= 0) then
        perimeter = 0
        do i = 1, boundary_nodes
            next = mod(i, boundary_nodes) + 1
            perimeter = perimeter + norm2(contour(:, next) - contour(:, i))
        end do
        max_area = (perimeter/boundary_nodes)**2/2
    end if
    if (curved) then
        call initialize_gs(state, contour, max_area, status, degree, curve=boundary)
    else
        call initialize_gs(state, contour, max_area, status, degree)
    end if
    if (status /= 0) error stop 'mesh/assembly failed'
    if (len_trim(profile_file) > 0) then
        allocate (law)
        call read_gs_profile_law(trim(profile_file), law, status)
        if (status /= 0) error stop 'invalid cubic source profile file'
        f_edge = law%f_edge
        call solve_gs_cubic_profiles(state, mu0, law, status, nonlinear_tolerance, &
            quadrature_tolerance=quadrature_tolerance)
    else
        call solve_gs(state, mu0, pprime, ffprime, status)
    end if
    if (status /= 0) then
        print *, 'solve status, iteration, residual:', status, state%iterations, &
            state%residual
        error stop 'equilibrium solve failed'
    end if
    if (total_current /= 0) then
        if (.not. allocated(law)) error stop 'current constraint requires normalized law'
        call gs_normalize_current(state, mu0, law, total_current, measured_current, &
            status, quadrature_tolerance)
        if (status /= 0) error stop 'invalid normalized current constraint'
    end if
    call system_clock(finish)
    wall_s = real(finish - start, dp)/rate
    if (wall_s > 300) error stop 'solve exceeded 300 seconds'
    call system_clock(start)
    call gs_functional_majorant(state, mu0, pprime, ffprime, majorant, &
        flux_defect, balance_defect, majorant_status, law=law, &
        recovery_converged=recovery_converged, recovery_iterations=recovery_iterations)
    if (majorant_status /= 0 .and. majorant_status /= 4) &
        error stop 'recovered-flux recovery failed'
    call system_clock(finish)
    estimator_s = real(finish - start, dp)/rate
    call gs_readback_resolution(state, radial_points, angular_points)
    call read_gs_equilibrium(state, f_edge, ffprime, radial_points, angular_points, &
        readback, status, law)
    if (status /= 0) then
        print *, 'readback status:', status
        error stop 'equilibrium contour readback failed'
    end if
    call write_gs_mesh(state, trim(mesh_file), status)
    if (status /= 0) error stop 'mesh export failed'
    if (allocated(law)) call write_gs_profile_law('profile-law.dat', law)
    call run%open_resolved(unit)
    write (unit, nml=gs_phase1)
    call run%close_resolved(unit)
    call create_directory_runtime('output', directory_ok)
    if (.not. directory_ok) error stop 'cannot create output directory'
    call result%create(run)
    call result%attribute('convention', 'COCOS 3; SI; psi in Wb/rad; psi_edge=0')
    call result%write('boundary_samples', samples, 'm', 'supplied closed R Z curve')
    call result%write('points', state%points, 'm', &
        'physical isoparametric geometry nodes')
    call result%write('psi', state%psi, 'Wb/rad', 'native FE coefficients')
    call result%write('triangles', real(state%triangles, dp), '1', &
        'isoparametric connectivity')
    if (allocated(law)) then
        call result%write('source_knots', law%knots, 'profile coordinate', &
            'executed cubic-source knots; coordinate defined by source_normalized')
        call result%write('source_pprime', law%pp, 'Pa/(Wb/rad)', &
            'cubic coefficients of dp/dpsi in ascending coordinate powers')
        call result%write('source_ffprime', law%ff, 'T2 m2/(Wb/rad)', &
            'cubic coefficients of FFprime in ascending coordinate powers')
        call result%write('source_normalized', real(merge(1, 0, law%normalized), dp), &
            '1', '0 absolute psi; 1 normalized flux s=1-psi/psi_axis')
        call result%write('source_axis', law%axis, 'Wb/rad', &
            'normalized axis or absolute-law initialization seed')
        call result%write('source_p_edge', law%p_edge, 'Pa', 'edge pressure')
        call result%write('source_f_edge', law%f_edge, 'T m', 'signed edge F')
    end if
    call result%write('axis', readback%axis, 'm', 'physical magnetic axis')
    call result%write('volume', readback%volume, 'm3', 'revolved plasma volume')
    call result%write('phi_edge', readback%phi_edge, 'Wb', 'signed full toroidal flux')
    call result%write('s_tor', readback%s_tor, '1', 'normalized toroidal flux')
    call result%write('q', readback%q, '1', 'signed safety factor on FE contours')
    call result%write('profile_psi', readback%psi, 'Wb/rad', 'native contour levels')
    call result%write('contour_s', readback%contour_s, 's', 'contour readback time')
    call result%write('integral_s', readback%integral_s, 's', 'area integrals time')
    call result%write('wall_s', wall_s, 's', 'mesh assembly and solve wall time')
    call result%write('dof', real(count(.not. state%boundary), dp), '1', &
        'free unknowns')
    call result%write('iterations', real(state%iterations, dp), '1', &
        'linear CG iterations or nonlinear Newton steps')
    call result%write('residual', state%residual, '1', 'relative weak residual')
    call result%write('flux_consistency', readback%flux_consistency, '1', &
        'relative contour versus area flux quadrature difference')
    call result%write('majorant_status', real(majorant_status, dp), '1', &
        '0 energy estimate; 4 frozen-source defect: nonlinear stability unavailable')
    call result%write('majorant', majorant, 'Wb/m^0.5', 'recovered-flux energy error estimate')
    call result%write('majorant_flux', flux_defect, 'Wb/m^0.5', 'weighted flux defect')
    call result%write('majorant_balance', balance_defect, 'Wb/m^0.5', 'weighted balance defect')
    call result%write('majorant_recovery_converged', merge(1.0_dp, 0.0_dp, recovery_converged), &
        '1', 'recovery minimization met its algebraic tolerance')
    call result%write('majorant_iterations', real(recovery_iterations, dp), &
        '1', 'recovery CG iterations')
    call result%write('estimator_s', estimator_s, 's', 'recovered-flux estimator wall time')
    call result%close()
    open (newunit=unit, file='profiles.dat', status='replace')
    write (unit, '(a)') '# psi[Wb/rad] s_tor q p[Pa] F[Tm] pprime FFprime; COCOS 3'
    do i = 1, radial_points
        pressure = pprime*readback%psi(i); pp = pprime; ff = ffprime
        f = gs_profile_f(readback%psi(i), f_edge, ffprime)
        if (allocated(law)) then
            call gs_law_primitives(law, readback%psi(i), pressure, f, status)
            if (status /= 0) error stop 'nonphysical profile primitive'
            call gs_law_sources(law, readback%psi(i), pp, ff)
        end if
        write (unit, '(7es25.16e3)') readback%psi(i), readback%s_tor(i), &
            readback%q(i), pressure, f, pp, ff
    end do
    close (unit)
    call runtime_git(commit, dirty, 'gs_phase1')
    open (newunit=unit, file='summary.dat', status='replace')
    write (unit, '(a,1x,es25.16e3)') 'total_current', measured_current
    write (unit, '(a,1x,i0)') 'degree', state%element_order
    write (unit, '(a,1x,i0)') 'iterations', state%iterations
    write (unit, '(a,1x,i0)') 'dof', count(.not. state%boundary)
    write (unit, '(a,1x,es25.16e3)') 'wall_s', wall_s
    write (unit, '(a,1x,es25.16e3)') 'contour_s', readback%contour_s
    write (unit, '(a,1x,es25.16e3)') 'integral_s', readback%integral_s
    write (unit, '(a,1x,i0)') 'radial_points', radial_points
    write (unit, '(a,1x,i0)') 'angular_points', angular_points
    write (unit, '(a,1x,es25.16e3)') 'psi_axis', readback%psi_axis
    write (unit, '(a,1x,i0)') 'majorant_status', majorant_status
    write (unit, '(a,1x,es25.16e3)') 'axis_R', readback%axis(1)
    write (unit, '(a,1x,es25.16e3)') 'axis_Z', readback%axis(2)
    write (unit, '(a,1x,es25.16e3)') 'volume', readback%volume
    write (unit, '(a,1x,es25.16e3)') 'phi_edge', readback%phi_edge
    write (unit, '(a,1x,es25.16e3)') 'residual', state%residual
    write (unit, '(a,1x,es25.16e3)') 'flux_consistency', readback%flux_consistency
    write (unit, '(a,1x,es25.16e3)') 'majorant', majorant
    write (unit, '(a,1x,es25.16e3)') 'majorant_flux', flux_defect
    write (unit, '(a,1x,es25.16e3)') 'majorant_balance', balance_defect
    write (unit, '(a,1x,i0)') 'majorant_recovery_converged', merge(1, 0, recovery_converged)
    write (unit, '(a,1x,i0)') 'majorant_iterations', recovery_iterations
    write (unit, '(a,1x,es25.16e3)') 'estimator_s', estimator_s
    write (unit, '(a,1x,a)') 'code_commit', commit
    write (unit, '(a,1x,a)') 'code_tree', dirty
    close (unit)
    print '(a,i0,a,es13.5)', 'DOF=', count(.not. state%boundary), ' wall_s=', wall_s
contains
    subroutine evaluate_source_file()
        character(512) :: path
        real(dp) :: psi
        integer :: count, i, stat
        if (run%argument_count() /= 2) error stop '--source profile-law.dat'
        call run%argument(2, path)
        allocate (law)
        call read_gs_profile_law(trim(path), law, stat)
        if (stat /= 0) error stop 'invalid executed profile law'
        read (input_unit, *) count
        do i = 1, count
            read (input_unit, *) psi
            call gs_law_primitives(law, psi, pressure, f, stat)
            if (stat /= 0) error stop 'nonphysical profile primitive'
            call gs_law_sources(law, psi, pp, ff)
            write (*, '(4es25.16e3)') pressure, f, pp, ff
        end do
    end subroutine

    subroutine boundary(t, point)
        real(dp), intent(in) :: t
        real(dp), intent(out) :: point(2)
        call curve_value(curve, t, point)
    end subroutine

    subroutine evaluate_q_file()
        character(512) :: path
        real(dp) :: row(3), s
        integer :: n, i, unit, stat
        if (run%argument_count() /= 2) &
            error stop '--q profiles.dat; stdin: count, s_tor'
        call run%argument(2, path)
        open (newunit=unit, file=trim(path), status='old', action='read')
        read (unit, *)
        n = 0
        do
            read (unit, *, iostat=stat) row
            if (stat /= 0) exit
            n = n + 1
        end do
        if (n < 4) error stop 'incomplete q profile'
        rewind (unit)
        read (unit, *)
        allocate (readback%psi(n), readback%s_tor(n), readback%q(n))
        do i = 1, n
            read (unit, *) readback%psi(i), readback%s_tor(i), readback%q(i)
        end do
        close (unit)
        read (input_unit, *) n
        do i = 1, n
            read (input_unit, *) s
            write (*, '(es25.16e3)') gs_readback_q(readback, s)
        end do
    end subroutine

    subroutine evaluate_file(streaming)
        logical, intent(in) :: streaming
        character(512) :: path, text
        real(dp) :: f, nan
        real(dp), allocatable :: points(:, :), values(:, :)
        integer, allocatable :: statuses(:)
        integer :: count, i, stat
        if (run%argument_count() /= 4 .and. run%argument_count() /= 5) &
            error stop '--evaluate mesh F_edge FFprime [profile-law.dat]'
        if (run%argument_count() == 5) then
            allocate (law)
            call run%argument(5, path)
            call read_gs_profile_law(trim(path), law, stat)
            if (stat /= 0) error stop 'invalid executed profile law'
        end if
        call run%argument(2, path)
        call run%argument(3, text)
        read (text, *) f_edge
        call run%argument(4, text)
        read (text, *) ffprime
        call read_gs_mesh(state, path, stat)
        if (stat /= 0) error stop 'invalid native field export'
        nan = ieee_value(0.0_dp, ieee_quiet_nan)
        do
            read (input_unit, *, iostat=stat) count
            if (stat < 0) exit
            if (stat /= 0 .or. count < 0) error stop 'invalid point count'
            if (allocated(points)) then
                if (size(points, 2) /= count) deallocate (points, values, statuses)
            end if
            if (.not. allocated(points)) &
                allocate (points(2, count), values(4, count), statuses(count))
            do i = 1, count
                read (input_unit, *) points(:, i)
            end do
            call sample_gs_many(state, points, f_edge, values, statuses)
            do i = 1, count
                if (statuses(i) == 0) then
                    f = gs_profile_f(values(1, i), f_edge, ffprime, law)
                    if (abs(f) <= tiny(1.0_dp)) statuses(i) = 1
                    values(4, i) = f/points(1, i)
                end if
                if (statuses(i) /= 0) values(:, i) = nan
                write (*, '(4es25.16e3)') values(:, i)
            end do
            flush (output_unit)
            if (.not. streaming) exit
        end do
    end subroutine
end program
