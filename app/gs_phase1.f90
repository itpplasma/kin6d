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
    use kin6d_gs_certificate, only: gs_fixed_certifier_t
    use kin6d_gs_error_control, only: gs_linear_goal_t, gs_fixed_linear_goal
    use kin6d_grad_shafranov_leaves, only: gs_source
    use kin6d_equilibrium_certificate, only: equilibrium_certificate_t, &
        certify, claim_name
    use kin6d_gs_profile_law, only: gs_profile_law_t, read_gs_profile_law, &
        write_gs_profile_law, gs_law_sources, gs_law_primitives
    use kin6d_gs_nonlinear, only: solve_gs_cubic_profiles, gs_normalize_current
    use kin6d_gs_prescribed_q, only: solve_gs_prescribed_q
    use kin6d_result_io, only: run_t, result_file_t, runtime_git
    use kin6d_gs_timers, only: gs_wall_s, gs_timer_total, gs_timer_calls, &
        gs_timer_iterations, gs_timer_event_count, gs_timer_event_get, &
        gs_timer_reset, gs_progress_set, gs_progress, GS_T_ASSEMBLY, GS_T_FIXED_CG, &
        GS_T_JAC_CG, GS_T_JAC_GMRES, GS_T_PRECOND_CG, GS_T_FACTOR, GS_T_DIRECT, GS_T_LOAD_CELLS
    use kin6d_gs_options, only: gs_options_text
    use fortplot_system_runtime, only: create_directory_runtime
    !$ use omp_lib, only: omp_get_max_threads
    implicit none
    type(gs_state_t), target :: state
    type(gs_fixed_certifier_t) :: certifier
    type(equilibrium_certificate_t) :: certificate
    type(gs_linear_goal_t) :: flux_goal
    type(gs_profile_law_t), allocatable :: law
    type(gs_profile_law_t) :: pressure_law
    type(periodic_curve_t) :: curve
    type(gs_readback_t) :: readback
    type(run_t) :: run
    type(result_file_t) :: result
    real(dp), allocatable :: samples(:, :), contour(:, :), interior_points(:, :)
    real(dp), allocatable :: q_s(:), q_values(:), inverse_knots(:)
    real(dp) :: max_area = 0, pprime = 0, ffprime = 0, f_edge = 0
    real(dp) :: mu0 = 1.2566370614359173e-6_dp, wall_s, perimeter
    real(dp) :: nonlinear_tolerance = 1.0e-9_dp, pressure, f, pp, ff
    real(dp) :: quadrature_tolerance = 1.0e-8_dp
    real(dp) :: solve_wall_limit_s = 300.0_dp
    real(dp) :: total_current = 0, measured_current = 0
    real(dp) :: majorant, flux_defect, balance_defect, estimator_s, stability_estimate
    real(dp) :: phi_target = 0, inverse_q_error = 0, q_tolerance = 1.0e-3_dp
    real(dp) :: inverse_lambda = 0, inverse_lambda_used = 0, inverse_q_accuracy = 0
    integer :: inverse_intervals = 2, inverse_evaluations = 0, inverse_newton = 0
    integer :: degree = 2
    integer :: boundary_nodes = 96, radial_points = 0, angular_points = 0
    integer :: unit, status, n, i, next, majorant_status, recovery_iterations
    integer(int64) :: start, finish, rate
    real(dp) :: mesh_s = 0, newton_s = 0, readback_s = 0, io_s = 0, io_total_s = 0
    real(dp) :: clock, program_start_s
    integer :: threads = 1
    character(16) :: determinism = 'fast'
    logical :: curved = .true., directory_ok, recovery_converged
    logical :: inverse_pressure_knots = .false., inverse_verbose = .false.
    character(512) :: boundary_file = '', mesh_file = 'mesh.dat', profile_file = ''
    character(512) :: q_file = '', interior_points_file = '', pprime_law_file = ''
    character(512) :: inverse_knots_file = ''
    character(16) :: inverse_jacobian = 'exact', inverse_regularization = 'none'
    character(256) :: message
    character(16) :: phase
    character(512) :: argument
    character(:), allocatable :: commit, dirty
    namelist /gs_phase1/ boundary_file, mesh_file, boundary_nodes, max_area, &
        pprime, ffprime, f_edge, mu0, curved, radial_points, angular_points, degree, &
        profile_file, nonlinear_tolerance, quadrature_tolerance, total_current, q_file, &
        phi_target, inverse_intervals, q_tolerance, solve_wall_limit_s, &
        interior_points_file, pprime_law_file, inverse_knots_file, &
        inverse_pressure_knots, inverse_jacobian, inverse_verbose, &
        inverse_regularization, inverse_lambda, inverse_q_accuracy
    call run%start('gs_phase1')
    program_start_s = gs_wall_s()
    call gs_timer_reset()
    call gs_progress_set(.true.)
    !$ threads = omp_get_max_threads()
    call get_environment_variable('KIN6D_DETERMINISM', determinism)
    if (len_trim(determinism) == 0) determinism = 'fast'
    if (trim(determinism) /= 'deterministic' .and. trim(determinism) /= 'fast') &
        error stop 'KIN6D_DETERMINISM must be fast (default) or deterministic (debug)'
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
    clock = gs_wall_s()
    write (output_unit, '(a,i0,a,a)') 'gs_phase1 start: omp_threads=', threads, &
        ' determinism=', trim(determinism)
    flush (output_unit)
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
    allocate (interior_points(2, 0))
    if (len_trim(interior_points_file) > 0) then
        open (newunit=unit, file=trim(interior_points_file), status='old', action='read')
        read (unit, *) n
        if (n < 0) error stop 'negative interior point count'
        deallocate (interior_points); allocate (interior_points(2, n))
        do i = 1, n
            read (unit, *) interior_points(:, i)
        end do
        close (unit)
    end if
    io_s = gs_wall_s() - clock
    clock = gs_wall_s()
    if (curved) then
        call initialize_gs(state, contour, max_area, status, degree, curve=boundary, &
            interior_points=interior_points)
    else
        call initialize_gs(state, contour, max_area, status, degree, &
            interior_points=interior_points)
    end if
    if (status /= 0) error stop 'mesh/assembly failed'
    mesh_s = gs_wall_s() - clock
    write (output_unit, '(a,i0,a,f10.3)') 'mesh done: dof=', &
        count(.not. state%boundary), ' mesh_s=', mesh_s
    flush (output_unit)
    if (len_trim(pprime_law_file) > 0 .and. len_trim(q_file) == 0) &
        error stop 'pprime_law_file requires q_file'
    if (len_trim(q_file) > 0) then
        if (len_trim(profile_file) > 0) error stop 'q_file and profile_file conflict'
        open (newunit=unit, file=trim(q_file), status='old', action='read')
        read (unit, *) n
        if (n < 4) error stop 'q table requires at least four rows'
        allocate (q_s(n), q_values(n), law)
        do i = 1, n
            read (unit, *) q_s(i), q_values(i)
        end do
        close (unit)
        call gs_readback_resolution(state, radial_points, angular_points)
        if (len_trim(inverse_knots_file) > 0 .and. inverse_pressure_knots) &
            error stop 'inverse_knots_file and inverse_pressure_knots conflict'
        if (len_trim(inverse_knots_file) > 0) then
            ! Row count followed by ascending FFprime spline break points 0..1.
            open (newunit=unit, file=trim(inverse_knots_file), status='old', &
                  action='read')
            read (unit, *) n
            if (n < 2) error stop 'knot file requires at least two rows'
            allocate (inverse_knots(n))
            read (unit, *) inverse_knots
            close (unit)
        end if
        clock = gs_wall_s()
        if (len_trim(pprime_law_file) > 0) then
            ! Variable p'(s): the pp rows of a normalized law file replace pprime.
            if (pprime /= 0) error stop 'pprime and pprime_law_file conflict'
            call read_gs_profile_law(trim(pprime_law_file), pressure_law, status)
            if (status /= 0 .or. .not. pressure_law%normalized) &
                error stop 'invalid normalized pprime_law_file'
            if (inverse_pressure_knots) inverse_knots = pressure_law%knots
        else if (inverse_pressure_knots) then
            error stop 'inverse_pressure_knots requires pprime_law_file'
        end if
        if (len_trim(pprime_law_file) > 0) then
            call solve_gs_prescribed_q(state, mu0, pprime, q_s, q_values, phi_target, &
                law, status, intervals=inverse_intervals, radial_points=radial_points, &
                angular_points=angular_points, tolerance=nonlinear_tolerance, &
                q_error=inverse_q_error, evaluations=inverse_evaluations, &
                q_tolerance=q_tolerance, pressure_law=pressure_law, &
                knots=inverse_knots, jacobian=trim(inverse_jacobian), &
                newton_iterations=inverse_newton, verbose=inverse_verbose, &
                regularization=trim(inverse_regularization), &
                regularization_weight=inverse_lambda, &
                regularization_used=inverse_lambda_used, q_accuracy=inverse_q_accuracy)
        else
            call solve_gs_prescribed_q(state, mu0, pprime, q_s, q_values, phi_target, &
                law, status, intervals=inverse_intervals, radial_points=radial_points, &
                angular_points=angular_points, tolerance=nonlinear_tolerance, &
                q_error=inverse_q_error, evaluations=inverse_evaluations, &
                q_tolerance=q_tolerance, knots=inverse_knots, &
                jacobian=trim(inverse_jacobian), newton_iterations=inverse_newton, &
                verbose=inverse_verbose, regularization=trim(inverse_regularization), &
                regularization_weight=inverse_lambda, &
                regularization_used=inverse_lambda_used, q_accuracy=inverse_q_accuracy)
        end if
        f_edge = law%f_edge
        newton_s = gs_wall_s() - clock
    else if (len_trim(profile_file) > 0) then
        allocate (law)
        call read_gs_profile_law(trim(profile_file), law, status)
        if (status /= 0) error stop 'invalid cubic source profile file'
        f_edge = law%f_edge
        clock = gs_wall_s()
        call solve_gs_cubic_profiles(state, mu0, law, status, nonlinear_tolerance, &
            quadrature_tolerance=quadrature_tolerance)
        newton_s = gs_wall_s() - clock
    else
        clock = gs_wall_s()
        call solve_gs(state, mu0, pprime, ffprime, status)
        newton_s = gs_wall_s() - clock
    end if
    if (len_trim(q_file) > 0) print '(a,i0,a,es10.3,a,es10.3,a,f10.1)', &
        'inverse status ', status, ' q_error ', inverse_q_error, ' lambda ', &
        inverse_lambda_used, &
        ' newton_s ', newton_s
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
    ! Report the converged state before the wall limit is enforced: a solve
    ! that converged slowly still has to leave its numbers behind.
    print '(a,i0,a,es13.5,a,es13.5)', 'DOF=', count(.not. state%boundary), &
        ' wall_s=', wall_s, ' residual=', state%residual
    if (wall_s > solve_wall_limit_s) error stop 'solve exceeded solve_wall_limit_s'
    flush (output_unit)
    clock = gs_wall_s()
    call write_gs_mesh(state, trim(mesh_file), status)
    if (status /= 0) error stop 'mesh export failed'
    if (allocated(law)) call write_gs_profile_law('profile-law.dat', law)
    io_s = io_s + gs_wall_s() - clock
    call system_clock(start)
    certifier%state => state
    certifier%mu0 = mu0; certifier%pprime = pprime; certifier%ffprime = ffprime
    if (allocated(law)) certifier%law = law
    call certify(certifier, certificate, status)
    if (status /= 0) error stop 'GS error estimator failed'
    majorant = certificate%bound
    majorant_status = certificate%lane_status
    flux_defect = certifier%flux_defect; balance_defect = certifier%balance_defect
    stability_estimate = certificate%stability
    recovery_converged = certifier%recovery_converged
    recovery_iterations = certifier%recovery_iterations
    if (.not. allocated(law)) then
        call gs_fixed_linear_goal(state, state%psi, fixed_source, unit_weight, &
                                  flux_goal, status)
        if (status /= 0) flux_goal%available = .false.
    end if
    call system_clock(finish)
    estimator_s = real(finish - start, dp)/rate
    write (output_unit, '(a,f10.3,a,i0)') 'majorant done: estimator_s=', estimator_s, &
        ' recovery_iterations=', recovery_iterations
    flush (output_unit)
    clock = gs_wall_s()
    call gs_readback_resolution(state, radial_points, angular_points)
    call read_gs_equilibrium(state, f_edge, ffprime, radial_points, angular_points, &
        readback, status, law)
    readback_s = gs_wall_s() - clock
    write (output_unit, '(a,f10.3)') 'readback done: readback_s=', readback_s
    flush (output_unit)
    if (status /= 0) then
        print *, 'readback status:', status
        error stop 'equilibrium contour readback failed'
    end if
    clock = gs_wall_s()
    call run%open_resolved(unit)
    write (unit, nml=gs_phase1)
    call run%close_resolved(unit)
    call create_directory_runtime('output', directory_ok)
    if (.not. directory_ok) error stop 'cannot create output directory'
    call result%create(run)
    call result%attribute('convention', 'COCOS 3; SI; psi in Wb/rad; psi_edge=0')
    call result%attribute('omp_threads', int_text(threads))
    call result%attribute('determinism', trim(determinism))
    call result%attribute('gs_linear_algebra', gs_options_text())
    call result%write('boundary_samples', samples, 'm', 'supplied closed R Z curve')
    call result%write('points', state%points, 'm', &
        'physical isoparametric geometry nodes')
    call result%write('psi', state%psi, 'Wb/rad', 'native FE coefficients')
    call result%write('triangles', real(state%triangles, dp), '1', &
        'isoparametric connectivity')
    if (allocated(q_s)) then
        call result%write('prescribed_q_s', q_s, '1', 'input normalized poloidal flux')
        call result%write('prescribed_q', q_values, '1', 'input signed COCOS 3 q')
        call result%write('phi_target', phi_target, 'Wb', 'prescribed toroidal flux')
        call result%write('inverse_q_error', inverse_q_error, '1', &
            'maximum sampled relative q constraint defect')
        call result%write('inverse_lambda', inverse_lambda_used, '1', &
            'FFprime smoothness regularization weight used (0: exact moments)')
    end if
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
    call result%write('omp_threads', real(threads, dp), '1', 'OpenMP maximum threads')
    call result%write('mesh_s', mesh_s, 's', 'mesh generation and element setup')
    call result%write('assembly_s', gs_timer_total(GS_T_ASSEMBLY), 's', &
        'source and Jacobian assembly (gs_cubic_load)')
    call result%write('assembly_cells_s', gs_timer_total(GS_T_LOAD_CELLS), 's', &
        'element loop inside gs_cubic_load (cell integration only)')
    call result%write('assembly_calls', real(gs_timer_calls(GS_T_ASSEMBLY), dp), '1', &
        'gs_cubic_load calls')
    call result%write('fixed_cg_s', gs_timer_total(GS_T_FIXED_CG), 's', &
        'fixed-source linear CG solves')
    call result%write('fixed_cg_iterations', &
        real(gs_timer_iterations(GS_T_FIXED_CG), dp), '1', 'fixed-source CG iterations')
    call result%write('jacobian_cg_s', gs_timer_total(GS_T_JAC_CG), 's', &
        'Newton Jacobian symmetric CG solves')
    call result%write('jacobian_cg_iterations', &
        real(gs_timer_iterations(GS_T_JAC_CG), dp), '1', 'Jacobian CG iterations')
    call result%write('jacobian_gmres_s', gs_timer_total(GS_T_JAC_GMRES), 's', &
        'Newton Jacobian GMRES fallback, including nested preconditioner CG')
    call result%write('jacobian_gmres_iterations', &
        real(gs_timer_iterations(GS_T_JAC_GMRES), dp), '1', 'GMRES iterations')
    call result%write('precond_cg_s', gs_timer_total(GS_T_PRECOND_CG), 's', &
        'nested preconditioner CG solves (inside jacobian_gmres_s)')
    call result%write('precond_cg_iterations', &
        real(gs_timer_iterations(GS_T_PRECOND_CG), dp), '1', 'nested CG iterations')
    call result%write('jacobian_cg_calls', real(gs_timer_calls(GS_T_JAC_CG), dp), '1', &
        'Jacobian CG solves')
    call result%write('jacobian_gmres_calls', &
        real(gs_timer_calls(GS_T_JAC_GMRES), dp), '1', 'Jacobian GMRES solves')
    call result%write('factor_s', gs_timer_total(GS_T_FACTOR), 's', &
        'Cholesky factorisation of the stiffness matrix')
    call result%write('direct_solves', real(gs_timer_calls(GS_T_DIRECT), dp), '1', &
        'back-substitutions with the cached factor')
    call result%write('direct_s', gs_timer_total(GS_T_DIRECT), 's', &
        'time in back-substitutions with the cached factor')
    call result%write('newton_s', newton_s, 's', &
        'whole nonlinear or linear solve call, including assembly and linear solves')
    call result%write('readback_s', readback_s, 's', 'equilibrium readback')
    call result%write('io_s', io_s, 's', &
        'input reading, mesh export and resolved config; summary.dat io_total_s adds '// &
        'NetCDF and profiles.dat')
    call result%write('dof', real(count(.not. state%boundary), dp), '1', &
        'free unknowns')
    call result%write('iterations', real(state%iterations, dp), '1', &
        'linear CG iterations or nonlinear Newton steps')
    call result%write('residual', state%residual, '1', 'relative weak residual')
    call result%write('flux_consistency', readback%flux_consistency, '1', &
        'relative contour versus area flux quadrature difference')
    call result%attribute('integrated_psi_goal_status', &
                          merge('approximate', 'unavailable', flux_goal%available))
    if (flux_goal%available) then
        call result%write('integrated_psi', flux_goal%value, 'Wb m2', &
                          'area integral of poloidal flux on represented domain')
        call result%write('integrated_psi_correction', flux_goal%correction, 'Wb m2', &
                          'adjoint weighted residual correction')
        call result%write('integrated_psi_radius', flux_goal%radius, 'Wb m2', &
                          'primal-adjoint energy product; approximate quadrature')
        call result%write('integrated_psi_lower', flux_goal%lower, 'Wb m2', &
                          'practical lower endpoint; not a pointwise flux band')
        call result%write('integrated_psi_upper', flux_goal%upper, 'Wb m2', &
                          'practical upper endpoint; not a pointwise flux band')
    end if
    call result%attribute('majorant_method', trim(certifier%method))
    call result%attribute('majorant_claim', trim(claim_name(certificate%claim)))
    call result%attribute('majorant_quadrature_status', &
        trim(certifier%quadrature_status))
    call result%attribute('majorant_stability_status', &
                         merge('available  ', 'unavailable', certificate%stability > 0))
    call result%attribute('majorant_scope', &
                     merge('solution energy error     ', 'frozen-source energy error', &
                                certificate%stability > 0))
    call result%write('majorant_status', real(majorant_status, dp), '1', &
                      '4 frozen-source residual only; 5 practical solution estimate')
    call result%write('majorant', majorant, 'Wb/m^0.5', &
        'energy error estimate; see majorant_scope and majorant_claim')
    call result%write('stability_estimate', stability_estimate, '1', &
                   'stability divisor; residual quadrature remains separately declared')
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
    io_total_s = io_s + gs_wall_s() - clock
    call runtime_git(commit, dirty, 'gs_phase1')
    open (newunit=unit, file='summary.dat', status='replace')
    write (unit, '(a,1x,es25.16e3)') 'total_current', measured_current
    write (unit, '(a,1x,i0)') 'degree', state%element_order
    write (unit, '(a,1x,i0)') 'inverse_evaluations', inverse_evaluations
    write (unit, '(a,1x,i0)') 'inverse_newton_iterations', inverse_newton
    write (unit, '(a,1x,es25.16e3)') 'inverse_q_error', inverse_q_error
    write (unit, '(a,1x,es25.16e3)') 'inverse_lambda', inverse_lambda_used
    write (unit, '(a,1x,i0)') 'iterations', state%iterations
    write (unit, '(a,1x,i0)') 'dof', count(.not. state%boundary)
    write (unit, '(a,1x,es25.16e3)') 'wall_s', wall_s
    write (unit, '(a,1x,es25.16e3)') 'contour_s', readback%contour_s
    write (unit, '(a,1x,es25.16e3)') 'integral_s', readback%integral_s
    write (unit, '(a,1x,i0)') 'radial_points', radial_points
    write (unit, '(a,1x,i0)') 'angular_points', angular_points
    write (unit, '(a,1x,es25.16e3)') 'psi_axis', readback%psi_axis
    write (unit, '(a,1x,a)') 'integrated_psi_goal_status', &
        merge('approximate', 'unavailable', flux_goal%available)
    if (flux_goal%available) then
        write (unit, '(a,1x,es25.16e3)') 'integrated_psi', flux_goal%value
        write (unit, '(a,1x,es25.16e3)') &
            'integrated_psi_correction', flux_goal%correction
        write (unit, '(a,1x,es25.16e3)') 'integrated_psi_radius', flux_goal%radius
        write (unit, '(a,1x,es25.16e3)') 'integrated_psi_lower', flux_goal%lower
        write (unit, '(a,1x,es25.16e3)') 'integrated_psi_upper', flux_goal%upper
    end if
    write (unit, '(a,1x,a)') 'majorant_method', trim(certifier%method)
    write (unit, '(a,1x,a)') 'majorant_claim', trim(claim_name(certificate%claim))
    write (unit, '(a,1x,a)') 'majorant_quadrature_status', &
        trim(certifier%quadrature_status)
    write (unit, '(a,1x,a)') 'majorant_scope', &
        merge('solution-energy-error     ', 'frozen-source-energy-error', &
              certificate%stability > 0)
    write (unit, '(a,1x,i0)') 'majorant_status', majorant_status
    write (unit, '(a,1x,es25.16e3)') 'stability_estimate', stability_estimate
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
    write (unit, '(a,1x,i0)') 'omp_threads', threads
    write (unit, '(a,1x,a)') 'determinism', trim(determinism)
    write (unit, '(a,1x,a)') 'gs_linear_algebra', gs_options_text()
    write (unit, '(a,1x,es25.16e3)') 'mesh_s', mesh_s
    write (unit, '(a,1x,es25.16e3)') 'assembly_s', gs_timer_total(GS_T_ASSEMBLY)
    write (unit, '(a,1x,i0)') 'assembly_calls', gs_timer_calls(GS_T_ASSEMBLY)
    write (unit, '(a,1x,es25.16e3)') 'fixed_cg_s', gs_timer_total(GS_T_FIXED_CG)
    write (unit, '(a,1x,i0)') 'fixed_cg_iterations', gs_timer_iterations(GS_T_FIXED_CG)
    write (unit, '(a,1x,es25.16e3)') 'jacobian_cg_s', gs_timer_total(GS_T_JAC_CG)
    write (unit, '(a,1x,i0)') 'jacobian_cg_iterations', gs_timer_iterations(GS_T_JAC_CG)
    write (unit, '(a,1x,es25.16e3)') 'jacobian_gmres_s', gs_timer_total(GS_T_JAC_GMRES)
    write (unit, '(a,1x,i0)') 'jacobian_gmres_iterations', &
        gs_timer_iterations(GS_T_JAC_GMRES)
    write (unit, '(a,1x,es25.16e3)') 'precond_cg_s', gs_timer_total(GS_T_PRECOND_CG)
    write (unit, '(a,1x,i0)') 'precond_cg_iterations', &
        gs_timer_iterations(GS_T_PRECOND_CG)
    write (unit, '(a,1x,i0)') 'jacobian_cg_calls', gs_timer_calls(GS_T_JAC_CG)
    write (unit, '(a,1x,i0)') 'jacobian_gmres_calls', gs_timer_calls(GS_T_JAC_GMRES)
    write (unit, '(a,1x,es25.16e3)') 'factor_s', gs_timer_total(GS_T_FACTOR)
    write (unit, '(a,1x,i0)') 'direct_solves', gs_timer_calls(GS_T_DIRECT)
    write (unit, '(a,1x,es25.16e3)') 'direct_s', gs_timer_total(GS_T_DIRECT)
    write (unit, '(a,1x,es25.16e3)') 'assembly_cells_s', gs_timer_total(GS_T_LOAD_CELLS)
    write (unit, '(a,1x,es25.16e3)') 'newton_s', newton_s
    write (unit, '(a,1x,es25.16e3)') 'readback_s', readback_s
    write (unit, '(a,1x,es25.16e3)') 'io_total_s', io_total_s
    write (unit, '(a,1x,es25.16e3)') 'total_s', gs_wall_s() - program_start_s
    write (unit, '(a,1x,a)') 'code_commit', commit
    write (unit, '(a,1x,a)') 'code_tree', dirty
    close (unit)
    open (newunit=unit, file='timing.dat', status='replace')
    write (unit, '(a)') '# per-call log: phase seconds iterations (nested '// &
        'preconditioner CG is aggregated in summary.dat only)'
    do i = 1, gs_timer_event_count()
        call gs_timer_event_get(i, phase, clock, next)
        write (unit, '(a16,1x,es25.16e3,1x,i0)') phase, clock, next
    end do
    close (unit)
contains
    function fixed_source(cell, r, z, psi) result(source)
        integer, intent(in) :: cell
        real(dp), intent(in) :: r, z, psi
        real(dp) :: source
        call gs_source(r, mu0, pprime, ffprime, source)
        source = source + 0.0_dp*(cell + z + psi)
    end function fixed_source

    function unit_weight(r, z) result(weight)
        real(dp), intent(in) :: r, z
        real(dp) :: weight
        weight = 1.0_dp + 0.0_dp*(r + z)
    end function unit_weight

    function int_text(value) result(text)
        integer, intent(in) :: value
        character(16) :: text
        write (text, '(i0)') value
    end function

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
