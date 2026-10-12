program gs_free_main
    !! Free-boundary Grad-Shafranov application (DESIGN P11): FEM-BEM coupled
    !! solve on a mesh whose plasma-boundary chain is refitted onto the computed
    !! boundary, with an optional coil-current inverse. See
    !! docs/gs-free-boundary.md for the namelist and the output contract.
    use, intrinsic :: iso_fortran_env, only: dp => real64, input_unit, output_unit
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use kin6d_grad_shafranov, only: gs_state_t, initialize_gs, write_gs_mesh, &
                                    read_gs_mesh
    use kin6d_gs_field_readback, only: sample_gs_many
    use kin6d_curved_triangle, only: periodic_curve_t, initialize_periodic_curve, &
                                     curve_value
    use kin6d_gs_readback, only: gs_profile_f
    use kin6d_gs_profile_law, only: gs_profile_law_t, read_gs_profile_law, &
                                    write_gs_profile_law, gs_law_primitives, &
                                    gs_law_sources
    use kin6d_gs_coils, only: gs_circuit_set_t, gs_read_coil_table, &
                              gs_circuit_flux_vectors
    use kin6d_gs_free_boundary, only: gs_fb_config_t, gs_fb_system_t, gs_fb_eval_t, &
        gs_fb_refit_config_t, gs_fb_refit_report_t, gs_fb_setup, gs_fb_newton, &
        gs_fb_evaluate, gs_fb_refit_solve, gs_fb_refit_mesh, gs_fb_chain_movement, &
        gs_fb_trace_lcfs, GS_FB_OK, GS_FB_XPOINT, GS_FB_LIMITER, GS_FB_TARGET_MEAN
    use kin6d_gs_coil_inverse, only: gs_ci_row_t, gs_ci_model_t, gs_ci_config_t, &
        gs_ci_result_t, gs_ci_build_model, gs_ci_fit, GS_CI_ISOFLUX_PAIR, &
        GS_CI_XPOINT_BR, GS_CI_XPOINT_BZ, GS_CI_PSI_POINT, GS_CI_ISOFLUX_LCFS, &
        GS_CI_PSI_BOUND, GS_CI_UNCONSTRAINED, GS_CI_BOUNDED, GS_CI_FIXED_AT_LIMIT, &
        GS_CI_OK, GS_CI_NOT_CONVERGED
    use kin6d_result_io, only: run_t, result_file_t, runtime_git
    use fortplot_system_runtime, only: create_directory_runtime
    !$ use omp_lib, only: omp_get_max_threads
    implicit none

    integer, parameter :: max_circuits = 64, max_ties = 16
    type(run_t) :: run
    type(result_file_t) :: result
    type(gs_state_t) :: state
    type(gs_profile_law_t) :: law
    type(gs_circuit_set_t) :: set
    type(gs_fb_config_t) :: cfg
    type(gs_fb_system_t) :: sys
    type(gs_fb_eval_t) :: ev
    type(gs_fb_refit_config_t) :: rcfg
    type(gs_fb_refit_report_t) :: report
    type(periodic_curve_t) :: curve
    real(dp), allocatable :: prof(:, :)
    type(gs_ci_model_t) :: model
    type(gs_ci_config_t) :: ccfg
    type(gs_ci_result_t) :: fit_result
    integer :: fit_stat = 0
    type(gs_ci_row_t), allocatable :: rows(:)
    real(dp), allocatable :: samples(:, :), contour(:, :), chain(:, :), points(:, :)
    real(dp), allocatable :: cur(:), lcfs(:, :), seed(:), vac(:, :)
    real(dp) :: max_area = 0.0_dp, mu0 = 1.2566370614359173e-6_dp
    real(dp) :: plasma_ellipse(3) = 0.0_dp, axis_guess(2) = 0.0_dp
    real(dp) :: xpoint_guess(2) = 0.0_dp, search_radius = 0.4_dp, cut_fraction = 1.0_dp
    real(dp) :: refit_tolerance = 1.0e-3_dp, newton_tolerance = 1.0e-10_dp
    real(dp) :: quadrature_tolerance = 1.0e-9_dp, seed_amplitude = 1.0_dp
    real(dp) :: currents(max_circuits) = 0.0_dp
    real(dp) :: reference_currents(max_circuits) = 0.0_dp
    real(dp) :: lambda = 0.0_dp, relaxation = 1.0_dp, fit_step_tolerance = 1.0e-8_dp
    real(dp) :: volume = 0.0_dp, perimeter, h, wall_start, wall_s, psi_a, psi_b
    integer :: boundary_nodes = 96, degree = 2, chain_nodes = 48, axis_sign = 1
    integer :: max_refits = 4, max_newton_iterations = 30, max_fit_iterations = 30
    integer :: n_reference = 0, profile_points = 65
    integer :: tie_a(max_ties) = 0, tie_b(max_ties) = 0
    integer :: unit, status, i, j, n, nc, pass, variant, stat, threads = 1
    logical :: curved = .true., emulate_freegs = .false., directory_ok, inverse
    logical :: refit_converged = .false.
    character(512) :: boundary_file = '', coil_file = '', profile_file = ''
    character(512) :: constraints_file = '', points_file = '', mesh_file = 'mesh.dat'
    character(512) :: argument
    character(32) :: lcfs_mode = 'xpoint', fit = 'none'
    character(256) :: message
    character(:), allocatable :: commit, dirty
    namelist /gs_free/ boundary_file, boundary_nodes, max_area, degree, curved, &
        coil_file, profile_file, mu0, plasma_ellipse, chain_nodes, lcfs_mode, &
        axis_sign, axis_guess, xpoint_guess, search_radius, cut_fraction, points_file, &
        max_refits, refit_tolerance, newton_tolerance, quadrature_tolerance, &
        max_newton_iterations, seed_amplitude, currents, fit, constraints_file, &
        lambda, reference_currents, n_reference, tie_a, tie_b, emulate_freegs, &
        relaxation, max_fit_iterations, fit_step_tolerance, mesh_file, profile_points

    call run%start('gs_free')
    wall_start = wall_clock()
    !$ threads = omp_get_max_threads()
    if (run%argument_count() > 0) then
        call run%argument(1, argument)
        if (trim(argument) /= '--evaluate' .and. &
            trim(argument) /= '--evaluate-stream') &
            error stop 'expected --config or --evaluate'
        call evaluate_file(trim(argument) == '--evaluate-stream')
        stop
    end if
    status = 0
    message = ''
    call run%open_config(unit)
    if (run%has_config()) read (unit, nml=gs_free, iostat=status, iomsg=message)
    call run%close_config(unit, 'gs_free', status, message)
    call check_inputs()

    ! Inputs: coils, profile law, boundary curve, LCFS points.
    call gs_read_coil_table(trim(coil_file), set, stat, message)
    if (stat /= 0) error stop 'invalid coil table'
    nc = size(set%incidence, 2)
    if (nc > max_circuits) error stop 'too many circuits'
    call read_gs_profile_law(trim(profile_file), law, stat)
    if (stat /= 0 .or. .not. law%normalized) &
        error stop 'invalid normalized profile_file'
    call read_curve(trim(boundary_file), samples)
    call initialize_periodic_curve(curve, samples, stat)
    if (stat /= 0) error stop 'invalid closed curve'
    allocate (contour(2, boundary_nodes))
    do i = 1, boundary_nodes
        call boundary(real(i - 1, dp)/boundary_nodes, contour(:, i))
    end do
    perimeter = sum(norm2(contour - cshift(contour, 1, dim=2), dim=1))
    if (max_area <= 0.0_dp) max_area = 0.4_dp*(perimeter/boundary_nodes)**2
    allocate (chain(2, chain_nodes))
    do j = 1, chain_nodes
        chain(:, j) = [plasma_ellipse(1) + plasma_ellipse(2)* &
                       cos(2.0_dp*acos(-1.0_dp)*(j - 1)/chain_nodes), &
                       plasma_ellipse(3)*sin(2.0_dp*acos(-1.0_dp)*(j - 1)/chain_nodes)]
    end do
    cfg%axis_sign = axis_sign
    cfg%axis_guess = axis_guess
    cfg%xpoint_guess = xpoint_guess
    cfg%search_radius = search_radius
    cfg%cut_fraction = cut_fraction
    cfg%quadrature_tolerance = quadrature_tolerance
    cfg%jacobian_tolerance = quadrature_tolerance
    cfg%tolerance = newton_tolerance
    cfg%max_iterations = max_newton_iterations
    select case (trim(lcfs_mode))
    case ('xpoint')
        cfg%mode = GS_FB_XPOINT
    case ('limiter')
        cfg%mode = GS_FB_LIMITER
    case ('target_mean')
        cfg%mode = GS_FB_TARGET_MEAN
    case default
        error stop 'lcfs_mode must be xpoint, limiter or target_mean'
    end select
    if (cfg%mode /= GS_FB_XPOINT) then
        if (len_trim(points_file) == 0) error stop 'points_file required'
        call read_curve(trim(points_file), points)
        call move_alloc(points, cfg%points)
    end if
    rcfg%max_refits = max_refits
    rcfg%tolerance = refit_tolerance
    cur = currents(:nc)
    inverse = trim(fit) /= 'none'
    h = perimeter/boundary_nodes

    call initialize_gs(state, contour, max_area, stat, degree, curved=curved, &
                       curve=boundary, interior_chain=chain)
    if (stat /= 0) error stop 'mesh/assembly failed'
    allocate (seed(size(state%points, 2)), vac(size(state%points, 2), nc))
    call gs_circuit_flux_vectors(set, state%points(1, :), state%points(2, :), vac)
    do i = 1, size(seed)
        seed(i) = dot_product(vac(i, :), cur) + &
                  real(axis_sign, dp)*seed_amplitude* &
                  exp(-(((state%points(1, i) - plasma_ellipse(1)) &
                         /plasma_ellipse(2))**2 &
                        + (state%points(2, i)/plasma_ellipse(3))**2))
    end do
    state%psi = seed
    write (output_unit, '(a,i0,a,i0,a,i0)') 'gs_free start: omp_threads=', threads, &
        ' nodes=', size(state%points, 2), ' circuits=', nc
    flush (output_unit)

    if (.not. inverse) then
        call gs_fb_refit_solve(state, sys, law, cfg, cur, contour, max_area, chain, &
                               rcfg, report, status, set=set, curve=boundary)
        if (status /= GS_FB_OK .and. max_refits > 0 .and. &
            .not. (report%converged .or. report%stagnated)) &
            write (output_unit, '(a,i0)') &
                'refit loop ended unconverged, status ', status
        if (allocated(report%movement)) refit_converged = report%converged
        if (.not. allocated(report%movement) .or. &
            (status /= GS_FB_OK .and. status /= 10)) &
            error stop 'free-boundary solve failed'
    else
        call read_rows(trim(constraints_file), rows)
        ccfg%lambda = lambda
        ccfg%emulate_freegs = emulate_freegs
        ccfg%relaxation = relaxation
        ccfg%max_iterations = max_fit_iterations
        ccfg%step_tolerance = fit_step_tolerance
        if (n_reference > 0) ccfg%i_ref = reference_currents(:nc)
        call build_ties(ccfg)
        select case (trim(fit))
        case ('unconstrained')
            variant = GS_CI_UNCONSTRAINED
        case ('bounded')
            variant = GS_CI_BOUNDED
        case ('fixed_at_limit')
            variant = GS_CI_FIXED_AT_LIMIT
        case default
            error stop 'fit must be none, unconstrained, bounded or fixed_at_limit'
        end select
        allocate (report%movement(0:max_refits), report%newton_iterations(0:max_refits))
        report%movement = 0.0_dp
        report%newton_iterations = 0
        do pass = 0, max_refits
            call gs_fb_setup(state, sys, law, cfg, set, status)
            if (status /= GS_FB_OK) error stop 'free-boundary setup failed'
            call gs_ci_build_model(state, sys, set, rows, model, stat)
            if (stat /= GS_CI_OK) error stop 'constraint rows could not be located'
            call gs_fb_newton(state, sys, cur, status)
            if (status /= GS_FB_OK) error stop 'free-boundary Newton failed'
            call gs_ci_fit(state, sys, model, ccfg, variant, cur, state%psi, &
                           fit_result, fit_stat)
            if (fit_stat /= GS_CI_OK .and. fit_stat /= GS_CI_NOT_CONVERGED) &
                error stop 'coil-current fit failed'
            cur = fit_result%currents
            state%psi = fit_result%psi
            call gs_fb_evaluate(state, sys, state%psi, cur, ev, status, &
                                jacobian=.false.)
            if (status /= GS_FB_OK) error stop 'evaluation failed'
            report%movement(pass) = gs_fb_chain_movement(state, state%psi, &
                                                         ev%psi_bound, chain)
            write (output_unit, '(a,i0,a,es10.3,a,es10.3)') 'inverse pass ', pass, &
                ' phi=', fit_result%phi, ' chain distance/h=', report%movement(pass)
            flush (output_unit)
            report%refits = pass
            if (report%movement(pass) <= refit_tolerance) then
                refit_converged = .true.
                exit
            end if
            if (pass == max_refits) exit
            call gs_fb_refit_mesh(state, ev%axis, ev%psi_bound, real(axis_sign, dp), &
                                  contour, max_area, chain, rcfg, stat, boundary)
            if (stat /= GS_FB_OK) error stop 'mesh refit failed'
            cfg%axis_guess = ev%axis
            cfg%xpoint_guess = ev%xpoint
        end do
    end if

    ! Final evaluation, boundary curve, readback.
    call gs_fb_evaluate(state, sys, state%psi, cur, ev, status, jacobian=.false.)
    if (status /= GS_FB_OK) error stop 'final evaluation failed'
    psi_a = ev%psi_axis
    psi_b = ev%psi_bound
    call gs_fb_trace_lcfs(state, state%psi, ev%axis, psi_b, real(axis_sign, dp), 256, &
                          lcfs, stat)
    if (stat /= GS_FB_OK) error stop 'plasma boundary trace failed'
    call do_readback()
    call plasma_volume()
    wall_s = wall_clock() - wall_start
    call write_outputs()

contains

    real(dp) function wall_clock() result(t)
        integer(8) :: count, rate

        call system_clock(count, rate)
        t = real(count, dp)/real(rate, dp)
    end function wall_clock

    subroutine check_inputs()
        if (status /= 0) error stop 'invalid gs_free namelist'
        if (len_trim(boundary_file) == 0) error stop 'boundary_file required'
        if (len_trim(coil_file) == 0) error stop 'coil_file required'
        if (len_trim(profile_file) == 0) error stop 'profile_file required'
        if (degree /= 2 .and. degree /= 3) error stop 'degree must be 2 or 3'
        if (boundary_nodes < 8) error stop 'at least eight boundary nodes required'
        if (chain_nodes < 8) error stop 'at least eight chain nodes required'
        if (abs(axis_sign) /= 1) error stop 'axis_sign must be +1 or -1'
        if (plasma_ellipse(2) <= 0.0_dp .or. plasma_ellipse(3) <= 0.0_dp .or. &
            plasma_ellipse(1) <= plasma_ellipse(2)) &
            error stop 'plasma_ellipse = R0, a, b with R0 > a > 0 and b > 0 required'
        if (max_refits < 0) error stop 'max_refits must not be negative'
    end subroutine check_inputs

    subroutine boundary(t, point)
        real(dp), intent(in) :: t
        real(dp), intent(out) :: point(2)

        call curve_value(curve, t, point)
    end subroutine boundary

    subroutine read_curve(path, array)
        !! "count" then count rows "R Z" (the gs_phase1 boundary format).
        character(*), intent(in) :: path
        real(dp), allocatable, intent(out) :: array(:, :)
        integer :: u, count, k

        open (newunit=u, file=path, status='old', action='read')
        read (u, *) count
        if (count < 1) error stop 'point file needs at least one row'
        allocate (array(2, count))
        do k = 1, count
            read (u, *) array(:, k)
        end do
        close (u)
    end subroutine read_curve

    subroutine read_rows(path, list)
        !! Constraint rows, one per line: KIND R Z R2 Z2 TARGET SCALE with KIND
        !! isoflux_pair, xpoint_br, xpoint_bz, psi_point, isoflux_lcfs or psi_bound.
        !! '#' starts a comment. For single-point kinds R2 Z2 are ignored.
        character(*), intent(in) :: path
        type(gs_ci_row_t), allocatable, intent(out) :: list(:)
        type(gs_ci_row_t), allocatable :: work(:)
        character(512) :: line
        character(32) :: kind
        real(dp) :: v(6)
        integer :: u, ios, count, hash

        if (len_trim(path) == 0) error stop 'constraints_file required for a fit'
        open (newunit=u, file=path, status='old', action='read')
        allocate (work(1024))
        count = 0
        do
            read (u, '(a)', iostat=ios) line
            if (ios /= 0) exit
            hash = index(line, '#')
            if (hash > 0) line = line(:hash - 1)
            if (len_trim(line) == 0) cycle
            read (line, *, iostat=ios) kind, v
            if (ios /= 0) error stop 'invalid constraint row'
            count = count + 1
            if (count > size(work)) error stop 'too many constraint rows'
            select case (trim(kind))
            case ('isoflux_pair')
                work(count)%kind = GS_CI_ISOFLUX_PAIR
            case ('xpoint_br')
                work(count)%kind = GS_CI_XPOINT_BR
            case ('xpoint_bz')
                work(count)%kind = GS_CI_XPOINT_BZ
            case ('psi_point')
                work(count)%kind = GS_CI_PSI_POINT
            case ('isoflux_lcfs')
                work(count)%kind = GS_CI_ISOFLUX_LCFS
            case ('psi_bound')
                work(count)%kind = GS_CI_PSI_BOUND
            case default
                error stop 'unknown constraint kind'
            end select
            work(count)%p = v(1:2)
            work(count)%q = v(3:4)
            work(count)%target = v(5)
            work(count)%scale = v(6)
            if (work(count)%scale <= 0.0_dp) &
                error stop 'constraint scale must be positive'
        end do
        close (u)
        if (count < 1) error stop 'no constraint rows'
        list = work(:count)
    end subroutine read_rows

    subroutine build_ties(c)
        !! tie_a(k) = tie_b(k) /= 0 ties the currents of two circuits equal.
        type(gs_ci_config_t), intent(inout) :: c
        integer :: k, count

        count = count_nonzero(tie_a)
        if (count == 0) return
        allocate (c%tie_c(count, nc), c%tie_d(count))
        c%tie_c = 0.0_dp
        c%tie_d = 0.0_dp
        do k = 1, count
            if (tie_a(k) < 1 .or. tie_a(k) > nc .or. tie_b(k) < 1 .or. tie_b(k) > nc) &
                error stop 'tie circuit index out of range'
            c%tie_c(k, tie_a(k)) = 1.0_dp
            c%tie_c(k, tie_b(k)) = -1.0_dp
        end do
    end subroutine build_ties

    integer function count_nonzero(a) result(count)
        integer, intent(in) :: a(:)

        count = 0
        do while (count < size(a))
            if (a(count + 1) == 0) exit
            count = count + 1
        end do
    end function count_nonzero

    subroutine plasma_volume()
        !! V = 2 pi int_plasma R dA: the plasma current of a law with dp/dpsi = 1
        !! and F F' = 0 on the same cut cells.
        type(gs_fb_system_t) :: probe
        type(gs_fb_eval_t) :: pev
        integer :: st

        probe = sys
        probe%law%pp = 0.0_dp
        probe%law%pp(1, :) = 1.0_dp
        probe%law%ff = 0.0_dp
        call gs_fb_evaluate(state, probe, state%psi, cur, pev, st, jacobian=.false.)
        if (st /= GS_FB_OK) error stop 'volume evaluation failed'
        volume = 2.0_dp*acos(-1.0_dp)*pev%plasma_current
    end subroutine plasma_volume

    subroutine do_readback()
        !! Profile readback of the executed law on the plasma flux range: rows
        !! psi, s, p, F, dp/dpsi, FF' at profile_points normalised fluxes.
        type(gs_profile_law_t) :: executed
        real(dp) :: s, phi, pressure, f, pp, ff
        integer :: k, st

        executed = law
        executed%axis = psi_a - psi_b
        allocate (prof(6, profile_points))
        do k = 1, profile_points
            s = real(k - 1, dp)/(profile_points - 1)
            phi = executed%axis*(1.0_dp - s)
            call gs_law_primitives(executed, phi, pressure, f, st)
            if (st /= 0) error stop 'nonphysical profile primitive'
            call gs_law_sources(executed, phi, pp, ff)
            prof(:, k) = [psi_b + phi, s, pressure, f, pp, ff]
        end do
    end subroutine do_readback

    subroutine write_outputs()
        type(gs_state_t) :: shifted
        type(gs_profile_law_t) :: executed
        real(dp), allocatable :: coil_i(:)
        integer :: u, k

        ! Mesh and law for --evaluate: flux measured from the boundary.
        shifted = state
        shifted%psi = state%psi - psi_b
        call write_gs_mesh(shifted, trim(mesh_file), stat)
        if (stat /= 0) error stop 'mesh export failed'
        executed = law
        executed%axis = psi_a - psi_b
        call write_gs_profile_law('profile-law.dat', executed)
        call run%open_resolved(u)
        write (u, nml=gs_free)
        call run%close_resolved(u)
        call create_directory_runtime('output', directory_ok)
        if (.not. directory_ok) error stop 'cannot create output directory'
        call result%create(run)
        call result%attribute('convention', &
            'SI; psi in Wb/rad absolute (gauge psi=0 on the axis of symmetry); '// &
            'profile tables of the executed law on the plasma flux range')
        call result%attribute('lcfs_mode', trim(lcfs_mode))
        call result%attribute('fit', trim(fit))
        call result%attribute('coil_file', trim(coil_file))
        call result%attribute('profile_file', trim(profile_file))
        call result%attribute('omp_threads', int_text(threads))
        call runtime_git(commit, dirty, 'gs_free')
        call result%attribute('code_commit', commit)
        call result%attribute('code_tree', dirty)
        call result%write('points', state%points, 'm', 'isoparametric geometry nodes')
        call result%write('triangles', real(state%triangles, dp), '1', &
                          'isoparametric connectivity')
        call result%write('psi', state%psi, 'Wb/rad', 'poloidal flux at the nodes')
        call result%write('psi_axis', psi_a, 'Wb/rad', 'flux on the magnetic axis')
        call result%write('psi_bound', psi_b, 'Wb/rad', 'flux on the plasma boundary')
        call result%write('axis', ev%axis, 'm', 'magnetic axis R Z')
        call result%write('xpoint', ev%xpoint, 'm', 'X point (xpoint mode) or limiter')
        call result%write('plasma_current', ev%plasma_current, 'A', 'toroidal current')
        call result%write('lcfs', lcfs, 'm', 'plasma boundary, 256 points R Z')
        call result%write('chain', chain, 'm', 'Gamma_p chain of the final mesh')
        call result%write('circuit_currents', cur, 'A', &
                          'circuit currents in table order')
        coil_i = matmul(set%incidence, cur)
        call result%write('coil_currents', coil_i, 'A', 'per-coil conductor currents')
        call result%write('circuit_i_min', set%i_min, 'A', 'circuit lower limits')
        call result%write('circuit_i_max', set%i_max, 'A', 'circuit upper limits')
        call result%write('refits', real(report%refits, dp), '1', &
                          'mesh refits performed')
        call result%write('refit_converged', &
                          merge(1.0_dp, 0.0_dp, refit_converged), '1', &
                          '1 when the chain is within refit_tolerance*h of the '// &
                          'boundary')
        call result%write('refit_movement', report%movement(0:report%refits), '1', &
                          'chain distance from the boundary / h before each solve')
        call result%write('boundary_faces', real(ev%face_cells, dp), '1', &
                          'plasma cell faces on the boundary')
        call result%write('profile_psi', prof(1, :), 'Wb/rad', 'absolute flux levels')
        call result%write('profile_s', prof(2, :), '1', &
                          'normalised flux s = 1 - (psi - psi_b)/(psi_a - psi_b)')
        call result%write('profile_pressure', prof(3, :), 'Pa', &
                          'p(psi), edge pressure p_edge')
        call result%write('profile_f', prof(4, :), 'T m', 'F(psi), signed')
        call result%write('profile_pprime', prof(5, :), 'Pa/(Wb/rad)', 'dp/dpsi')
        call result%write('profile_ffprime', prof(6, :), 'T2 m2/(Wb/rad)', 'F F-prime')
        call result%write('volume', volume, 'm3', 'plasma volume, 2 pi int R dA')
        call result%write('wall_s', wall_s, 's', 'solve and readback wall time')
        if (inverse) then
            call result%write('fit_status', real(fit_stat, dp), '1', &
                              '0 converged, 5 iteration limit')
            call result%write('fit_iterations', real(fit_result%iterations, dp), '1', &
                              'Gauss-Newton iterations of the last fit')
            call result%write('fit_phi', fit_result%phi, '1', 'objective Phi')
            call result%write('fit_phi_data', fit_result%phi_data, '1', &
                              'data part of Phi')
            call result%write('fit_phi_reg', fit_result%phi_reg, '1', &
                              'regulariser part of Phi')
            call result%write('fit_residual', fit_result%residual, '1', &
                              'scaled constraint residuals')
            call result%write('fit_gradient', fit_result%gradient, '1', 'dPhi/dI')
            if (allocated(fit_result%pinned)) &
                call result%write('fit_pinned', &
                                  merge(1.0_dp, 0.0_dp, fit_result%pinned), &
                                  '1', 'circuits pinned at a limit')
            call result%write('gsvd_values', fit_result%gsvd%gsv, '1', &
                              'generalised singular values of (D, L), ascending')
            call result%write('gsvd_modes', fit_result%gsvd%modes, 'A', &
                              'coil combination of each mode (columns), |L v| = 1')
            call result%write('gsvd_factor_x', fit_result%gsvd%xmat, '1', &
                              'GSVD factor X')
            call result%write('gsvd_null_mode', real(fit_result%gsvd%null_mode, dp), &
                              '1', 'index of the smallest generalised singular value')
            call result%write('gsvd_condition', fit_result%gsvd%condition, '1', &
                              'ratio of largest to smallest generalised singular value')
        end if
        call result%close()
        open (newunit=u, file='summary.dat', status='replace')
        write (u, '(a,1x,es25.16e3)') 'psi_axis', psi_a
        write (u, '(a,1x,es25.16e3)') 'psi_bound', psi_b
        write (u, '(a,1x,es25.16e3)') 'plasma_current', ev%plasma_current
        write (u, '(a,1x,es25.16e3)') 'axis_R', ev%axis(1)
        write (u, '(a,1x,es25.16e3)') 'axis_Z', ev%axis(2)
        write (u, '(a,1x,i0)') 'refits', report%refits
        write (u, '(a,1x,i0)') 'refit_converged', merge(1, 0, refit_converged)
        write (u, '(a,1x,es25.16e3)') 'wall_s', wall_s
        write (u, '(a,1x,i0)') 'omp_threads', threads
        do k = 1, nc
            write (u, '(a,i0,1x,es25.16e3)') 'circuit_current_', k, cur(k)
        end do
        write (u, '(a,1x,a)') 'code_commit', commit
        write (u, '(a,1x,a)') 'code_tree', dirty
        close (u)
    end subroutine write_outputs

    function int_text(value) result(text)
        integer, intent(in) :: value
        character(16) :: text

        write (text, '(i0)') value
    end function int_text

    subroutine evaluate_file(streaming)
        !! gs_free --evaluate mesh.dat F_edge FFprime [profile-law.dat [psi_bound]]
        !! with the points of stdin (count, then R Z rows). The mesh carries the
        !! flux measured from the boundary; the answer is psi (absolute when
        !! psi_bound is given, else the measured flux), B_R, B_Z, B_T = F/R. F is
        !! the executed law inside the plasma and F_edge outside it.
        logical, intent(in) :: streaming
        type(gs_profile_law_t), allocatable :: elaw
        type(gs_state_t) :: emesh
        character(512) :: path, text
        real(dp) :: f, nan, f_edge, ffprime, offset
        real(dp), allocatable :: pts(:, :), values(:, :)
        integer, allocatable :: statuses(:)
        integer :: count, k, st

        if (run%argument_count() < 4 .or. run%argument_count() > 6) &
            error stop '--evaluate mesh F_edge FFprime [profile-law.dat [psi_bound]]'
        offset = 0.0_dp
        if (run%argument_count() >= 5) then
            allocate (elaw)
            call run%argument(5, path)
            call read_gs_profile_law(trim(path), elaw, st)
            if (st /= 0) error stop 'invalid executed profile law'
        end if
        if (run%argument_count() == 6) then
            call run%argument(6, text)
            read (text, *) offset
        end if
        call run%argument(2, path)
        call run%argument(3, text)
        read (text, *) f_edge
        call run%argument(4, text)
        read (text, *) ffprime
        call read_gs_mesh(emesh, path, st)
        if (st /= 0) error stop 'invalid native field export'
        nan = ieee_value(0.0_dp, ieee_quiet_nan)
        do
            read (input_unit, *, iostat=st) count
            if (st < 0) exit
            if (st /= 0 .or. count < 0) error stop 'invalid point count'
            if (allocated(pts)) then
                if (size(pts, 2) /= count) deallocate (pts, values, statuses)
            end if
            if (.not. allocated(pts)) allocate (pts(2, count), values(4, count), &
                                                statuses(count))
            do k = 1, count
                read (input_unit, *) pts(:, k)
            end do
            call sample_gs_many(emesh, pts, f_edge, values, statuses)
            do k = 1, count
                if (statuses(k) == 0) then
                    if (allocated(elaw)) then
                        f = gs_profile_f(values(1, k), f_edge, ffprime, elaw)
                        if (abs(f) <= tiny(1.0_dp)) f = f_edge
                    else
                        f = f_edge
                    end if
                    values(4, k) = f/pts(1, k)
                    values(1, k) = values(1, k) + offset
                end if
                if (statuses(k) /= 0) values(:, k) = nan
                write (*, '(4es25.16e3)') values(:, k)
            end do
            flush (output_unit)
            if (.not. streaming) exit
        end do
    end subroutine evaluate_file
end program gs_free_main
