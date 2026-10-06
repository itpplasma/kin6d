program transport_demo
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_transport_laboratory, only: transport_state_t, initialize_transport, &
        solve_transport, solve_transport_tangent, sample_transport, &
        sample_transport_arrays
    use kin6d_transport_certificate, only: transport_certificate_t, &
        certify_transport, certify_transport_tangent, transport_certificate_accepts, &
        transport_observable_band
    use kin6d_transport_sampling, only: transport_basis_t, &
        transport_sampling_workspace_t, initialize_transport_sampling, &
        prepare_transport_basis, sample_transport_input_cached, &
        certify_transport_cached, certify_transport_tangent_cached, &
        transport_linear_enclosure, certify_transport_linear, &
        certify_transport_tangent_linear
    use kin6d_transport_leaves, only: transport_source_affine, transport_plot_x, &
        transport_cosine_moment_cell_interval, transport_scalar_add_interval, &
        transport_accuracy_interval
    use kin6d_poisson_leaves, only: coupled_distribution_interval
    use fortnum_interval, only: interval_t, interval, rad
    use fortplot, only: figure, subplot, plot, fill_between, title, xlabel, &
        ylabel, legend, tight_layout, savefig, ylim, xlim, suptitle
    use kin6d_demo_style, only: blue, vermillion, black, grey, num
    use fortplot_system_runtime, only: create_directory_runtime
    use kin6d_result_io, only: run_t, result_file_t
    implicit none
    ! Linear plotted bands include a rigorous chord interpolation remainder.
    integer, parameter :: modes = 2
    integer :: nx = 96, coarse_cells = 16, fine_cells = 64, limit_cells = 16
    real(dp) :: velocity = 0.3125_dp, rtol = 0.01_dp, atol = 1.0e-8_dp
    real(dp) :: goal_scale = 0.25_dp, derivative_scale = 0.25_dp
    namelist /transport/ nx, coarse_cells, fine_cells, limit_cells, velocity, rtol, &
        atol, goal_scale, derivative_scale
    type(run_t) :: run
    type(result_file_t) :: result
    type :: case_t
        type(transport_state_t) :: state
        type(transport_certificate_t) :: certificate
        real(dp), allocatable :: dsource_cos(:, :), dsource_sin(:, :)
        real(dp), allocatable :: dcos(:, :), dsin(:, :)
    end type
    type(case_t) :: coarse, fine
    type(transport_state_t) :: limiting
    type(transport_certificate_t) :: limit_certificate, audit
    type(transport_basis_t) :: basis, point_basis, limit_basis
    type(transport_sampling_workspace_t) :: work, point_work, limit_work
    type(interval_t) :: reference
    type(interval_t), allocatable :: x_cells(:), sources(:), x_points(:)
    real(dp), allocatable :: band_lo(:), band_hi(:)
    real(dp), allocatable :: plot_source_lo(:), plot_source_hi(:)
    real(dp), allocatable :: plot_coarse_lo(:), plot_coarse_hi(:)
    real(dp), allocatable :: plot_fine_lo(:), plot_fine_hi(:)
    real(dp), allocatable :: plot_coarse_dlo(:), plot_coarse_dhi(:)
    real(dp), allocatable :: plot_fine_dlo(:), plot_fine_dhi(:)
    real(dp), allocatable :: x(:), source_line(:)
    real(dp), allocatable :: coarse_line(:), fine_line(:)
    real(dp), allocatable :: coarse_tangent(:), fine_tangent(:)
    real(dp), allocatable :: source_lo(:), source_hi(:)
    real(dp), allocatable :: coarse_lo(:), coarse_hi(:), fine_lo(:), fine_hi(:)
    real(dp), allocatable :: coarse_dlo(:), coarse_dhi(:)
    real(dp), allocatable :: fine_dlo(:), fine_dhi(:)
    real(dp), allocatable :: limit_line(:), limit_lo(:), limit_hi(:)
    real(dp), allocatable :: reference_lo(:), reference_hi(:), reference_line(:)
    real(dp) :: scale, band_max
    integer :: i, output, unit, ios
    logical :: ok
    character(256) :: msg
    character(64) :: label, velocity_text

    call run%start('transport_demo')
    if (run%argument_count() > 0) error stop 'usage: transport_demo [--config FILE]'
    ios = 0; msg = ''
    call run%open_config(unit)
    if (run%has_config()) read (unit, nml=transport, iostat=ios, iomsg=msg)
    call run%close_config(unit, 'transport', ios, msg)
    if (nx < 1 .or. min(coarse_cells, fine_cells, limit_cells) < 2) &
        error stop 'nx >= 1 and cell counts >= 2 required'
    if (.not. (abs(velocity) <= 1.0_dp)) error stop 'velocity outside v_domain [-1,1]'
    if (.not. (rtol >= 0.0_dp .and. atol >= 0.0_dp .and. goal_scale > 0.0_dp .and. &
        derivative_scale > 0.0_dp)) error stop 'invalid accuracy request'
    call run%open_resolved(unit)
    write (unit, nml=transport)
    call run%close_resolved(unit)
    allocate (x_cells(nx), sources(nx), x_points(nx + 1), band_lo(nx), band_hi(nx), &
        plot_source_lo(nx + 1), plot_source_hi(nx + 1), plot_coarse_lo(nx + 1), &
        plot_coarse_hi(nx + 1), plot_fine_lo(nx + 1), plot_fine_hi(nx + 1), &
        plot_coarse_dlo(nx + 1), plot_coarse_dhi(nx + 1), plot_fine_dlo(nx + 1), &
        plot_fine_dhi(nx + 1), x(nx + 1), source_line(nx + 1), coarse_line(nx + 1), &
        fine_line(nx + 1), coarse_tangent(nx + 1), fine_tangent(nx + 1), &
        source_lo(2*nx), source_hi(2*nx), coarse_lo(2*nx), coarse_hi(2*nx), &
        fine_lo(2*nx), fine_hi(2*nx), coarse_dlo(2*nx), coarse_dhi(2*nx), &
        fine_dlo(2*nx), fine_dhi(2*nx), limit_line(nx + 1), limit_lo(nx + 1), &
        limit_hi(nx + 1), reference_lo(nx + 1), reference_hi(nx + 1), &
        reference_line(nx + 1))
    velocity_text = '5/16'
    if (velocity /= 0.3125_dp) write (velocity_text, '(g0)') velocity
    call create_directory_runtime('output', ok)
    if (.not. ok) error stop 'cannot create output directory'
    call prepare_case(coarse, coarse_cells)
    call prepare_case(fine, fine_cells)
    call certify_transport(fine%state, audit, ok, strict=.true.)
    if (.not. ok) error stop 'strict-box transport audit failed'
    print '(a,es24.16)', 'strict_box_audit_state_radius=', audit%state_radius
    print '(a)', 'certificate_kind=continuous_sup_residual_inverse'
    print '(a)', 'source=exact_supplied_Fourier_x_P1_v; nu=1; E=0; v_domain=[-1,1]'
    print '(a)', 'nu_tangent=certified_continuous; adjoint_goal_gradient=unavailable'

    do i = 1, nx + 1
        call transport_plot_x(real(i - 1, dp), real(nx + 1, dp), x(i))
        x_points(i) = interval(x(i))
        call sample_transport(coarse%state, x(i), velocity, coarse_line(i), ok)
        if (.not. ok) error stop 'coarse sample failed'
        call sample_transport(fine%state, x(i), velocity, fine_line(i), ok)
        if (.not. ok) error stop 'fine sample failed'
        call sample_transport_arrays(fine%state%v, fine%state%source_cos, &
            fine%state%source_sin, x(i), velocity, source_line(i), ok)
        if (.not. ok) error stop 'source sample failed'
        call sample_transport_arrays(coarse%state%v, coarse%dcos, coarse%dsin, &
            x(i), velocity, coarse_tangent(i), ok)
        if (.not. ok) error stop 'coarse tangent sample failed'
        call sample_transport_arrays(fine%state%v, fine%dcos, fine%dsin, &
            x(i), velocity, fine_tangent(i), ok)
        if (.not. ok) error stop 'fine tangent sample failed'
    end do
    do i = 1, nx
        x_cells(i) = interval(x(i), x(i + 1))
    end do
    call initialize_transport_sampling(basis, work, modes, nx, ok)
    if (.not. ok) error stop 'plot basis initialization failed'
    call prepare_transport_basis(basis, x_cells, ok)
    if (.not. ok) error stop 'plot basis preparation failed'
    call sample_transport_input_cached(fine%state, basis, x_cells, &
        interval(velocity), work, sources, ok)
    if (.not. ok) error stop 'source batch interval failed'
    do i = 1, nx
        source_lo(2*i - 1:2*i) = sources(i)%lo
        source_hi(2*i - 1:2*i) = sources(i)%hi
    end do
    call certify_transport_cached(coarse%state, coarse%certificate, basis, x_cells, &
        interval(velocity), work, band_lo, band_hi, ok)
    if (.not. ok) error stop 'coarse batch certificate failed'
    call expand_band(band_lo, band_hi, coarse_lo, coarse_hi)
    call certify_transport_cached(fine%state, fine%certificate, basis, x_cells, &
        interval(velocity), work, band_lo, band_hi, ok)
    if (.not. ok) error stop 'fine batch certificate failed'
    call expand_band(band_lo, band_hi, fine_lo, fine_hi)
    call certify_transport_tangent_cached(coarse%state, coarse%certificate, &
        coarse%dcos, coarse%dsin, basis, x_cells, interval(velocity), &
        work, band_lo, band_hi, ok)
    if (.not. ok) error stop 'coarse tangent batch certificate failed'
    call expand_band(band_lo, band_hi, coarse_dlo, coarse_dhi)
    call certify_transport_tangent_cached(fine%state, fine%certificate, &
        fine%dcos, fine%dsin, basis, x_cells, interval(velocity), &
        work, band_lo, band_hi, ok)
    if (.not. ok) error stop 'fine tangent batch certificate failed'
    call expand_band(band_lo, band_hi, fine_dlo, fine_dhi)
    open (newunit=output, file='output/transport-cell-bands.csv', status='replace')
    write (output, '(a)') 'x_left,x_right,source_lower,source_upper,'// &
        'coarse_lower,coarse_upper,fine_lower,fine_upper,'// &
        'coarse_dnu_lower,coarse_dnu_upper,fine_dnu_lower,fine_dnu_upper'
    do i = 1, nx
        write (output, '(es24.16,11(",",es24.16))') x(i), x(i + 1), &
            sources(i)%lo, sources(i)%hi, coarse_lo(2*i), coarse_hi(2*i), &
            fine_lo(2*i), fine_hi(2*i), coarse_dlo(2*i), coarse_dhi(2*i), &
            fine_dlo(2*i), fine_dhi(2*i)
    end do
    close (output)

    call initialize_transport_sampling(point_basis, point_work, modes, nx + 1, ok)
    if (.not. ok) error stop 'point basis initialization failed'
    call prepare_transport_basis(point_basis, x_points, ok)
    if (.not. ok) error stop 'point basis preparation failed'
    call transport_linear_enclosure(fine%state, fine%state%source_cos, &
        fine%state%source_sin, point_basis, x_points, interval(velocity), point_work, &
        0.0_dp, plot_source_lo, plot_source_hi, ok)
    if (.not. ok) error stop 'linear source enclosure failed'
    call certify_transport_linear(coarse%state, coarse%certificate, point_basis, &
        x_points, interval(velocity), point_work, plot_coarse_lo, plot_coarse_hi, ok)
    if (.not. ok) error stop 'linear coarse certificate failed'
    call certify_transport_linear(fine%state, fine%certificate, point_basis, &
        x_points, interval(velocity), point_work, plot_fine_lo, plot_fine_hi, ok)
    if (.not. ok) error stop 'linear fine certificate failed'
    call certify_transport_tangent_linear(coarse%state, coarse%certificate, &
        coarse%dcos, coarse%dsin, point_basis, x_points, interval(velocity), &
        point_work, plot_coarse_dlo, plot_coarse_dhi, ok)
    if (.not. ok) error stop 'linear coarse tangent certificate failed'
    call certify_transport_tangent_linear(fine%state, fine%certificate, &
        fine%dcos, fine%dsin, point_basis, x_points, interval(velocity), &
        point_work, plot_fine_dlo, plot_fine_dhi, ok)
    if (.not. ok) error stop 'linear fine tangent certificate failed'
    open (newunit=output, file='output/transport-linear-bands.csv', status='replace')
    write (output, '(a)') 'x,source_lower,source_upper,'// &
        'coarse_lower,coarse_upper,fine_lower,fine_upper,'// &
        'coarse_dnu_lower,coarse_dnu_upper,fine_dnu_lower,fine_dnu_upper'
    do i = 1, nx + 1
        write (output, '(es24.16,10(",",es24.16))') x(i), &
            plot_source_lo(i), plot_source_hi(i), plot_coarse_lo(i), &
            plot_coarse_hi(i), plot_fine_lo(i), plot_fine_hi(i), &
            plot_coarse_dlo(i), plot_coarse_dhi(i), plot_fine_dlo(i), plot_fine_dhi(i)
    end do
    close (output)

    call figure(figsize=[12.8_dp, 6.4_dp])
    call suptitle('Damped transport with a supplied source: refining 16 to 64 cells '// &
        'tightens the certified bands of f and of its rate derivative')
    call subplot(1, 2, 1)
    call title('Distribution and source at v = '//trim(velocity_text))
    call xlabel('x  (period 2π)')
    call ylabel('f, g')
    call xlim(x(1), x(nx + 1))
    call envelope(x, plot_source_lo, plot_source_hi, grey)
    call envelope(x, plot_coarse_lo, plot_coarse_hi, vermillion)
    call envelope(x, plot_fine_lo, plot_fine_hi, blue)
    call plot(x, source_line, label='source g ± enclosure', linestyle=':', &
        color=grey, linewidth=1.5_dp)
    write (label, '(i0,a)') coarse_cells, ' cells ± bound'
    call plot(x, coarse_line, label=trim(label), linestyle='--', color=vermillion, &
        linewidth=1.5_dp)
    write (label, '(i0,a)') fine_cells, ' cells ± bound'
    call plot(x, fine_line, label=trim(label), color=blue, linewidth=2.0_dp)
    call legend(loc='lower right')
    call subplot(1, 2, 2)
    call title('Implicit derivative of f in the collision rate (source fixed)')
    call xlabel('x  (period 2π)')
    call ylabel('$\partial f/\partial$(rate) at v = '//trim(velocity_text))
    call xlim(x(1), x(nx + 1))
    call envelope(x, plot_coarse_dlo, plot_coarse_dhi, vermillion)
    call envelope(x, plot_fine_dlo, plot_fine_dhi, blue)
    write (label, '(i0,a)') coarse_cells, ' cells ± tangent bound'
    call plot(x, coarse_tangent, label=trim(label), linestyle='--', color=vermillion, &
        linewidth=1.5_dp)
    write (label, '(i0,a)') fine_cells, ' cells ± tangent bound'
    call plot(x, fine_tangent, label=trim(label), color=blue, linewidth=2.0_dp)
    call legend(loc='upper left')
    call tight_layout()
    call savefig('output/general-transport.png')

    ! The earlier reference is used only in this separately declared slice.
    ! Its full velocity source differs; no reference enters either checker.
    call initialize_transport(limiting, limit_cells, 0, 1.0_dp, ok)
    if (.not. ok) error stop 'limiting initialization failed'
    limiting%source_cos(0, :) = 0.75_dp
    call solve_transport(limiting, ok)
    if (.not. ok) error stop 'limiting solve failed'
    call certify_transport(limiting, limit_certificate, ok)
    if (.not. ok) error stop 'limiting certificate failed'
    do i = 1, nx + 1
        call sample_transport(limiting, x(i), 0.0_dp, limit_line(i), ok)
        if (.not. ok) error stop 'limiting sample failed'
        call coupled_distribution_interval(interval(0.0_dp), interval(1.0_dp), &
            interval(x(i)), interval(0.0_dp), reference)
        reference_line(i) = reference%lo
    end do
    call initialize_transport_sampling(limit_basis, limit_work, 0, nx + 1, ok)
    if (.not. ok) error stop 'limiting point basis initialization failed'
    call prepare_transport_basis(limit_basis, x_points, ok)
    if (.not. ok) error stop 'limiting point basis preparation failed'
    call certify_transport_linear(limiting, limit_certificate, limit_basis, &
        x_points, interval(0.0_dp), limit_work, limit_lo, limit_hi, ok)
    if (.not. ok) error stop 'limiting linear certificate failed'
    call coupled_distribution_interval(interval(0.0_dp), interval(1.0_dp), &
        interval(x(1), x(nx + 1)), interval(0.0_dp), reference)
    reference_lo = reference%lo; reference_hi = reference%hi
    ! Both solutions equal 3/4 to rounding: plot offsets from the reference
    ! centre on a scale that resolves the slice bands.
    scale = max(maxval(abs(limit_lo - reference_line)), &
        maxval(abs(limit_hi - reference_line)), maxval(abs(reference_lo - &
        reference_line)), maxval(abs(reference_hi - reference_line)), tiny(1.0_dp))
    band_max = scale
    scale = 10.0_dp**floor(log10(scale))
    call figure(figsize=[9.6_dp, 6.0_dp])
    call subplot(1, 1, 1)
    call suptitle('Field-free zero-mode limit: the general solver reproduces '// &
        'the Krook value 3/4 within both bands')
    call title('Offset from the Krook reference at v = 0  (g = 3/4, collision rate 1)')
    call xlabel('x  (period 2π)')
    call ylabel('(f - 3/4) / '//num(scale, 0))
    call xlim(x(1), x(nx + 1)); call ylim(-1.8_dp*band_max/scale, 1.8_dp*band_max/scale)
    call fill_between(x, (reference_lo - reference_line)/scale, &
        (reference_hi - reference_line)/scale, color=grey, alpha=0.25_dp)
    call fill_between(x, (limit_lo - reference_line)/scale, &
        (limit_hi - reference_line)/scale, color=blue, alpha=0.25_dp)
    call plot(x, (limit_line - reference_line)/scale, label='general solver ± slice band', &
        color=blue, linewidth=2.0_dp)
    call plot(x, 0.0_dp*x, label='Krook reference ± band', linestyle='--', &
        color=black, linewidth=1.2_dp)
    call legend(loc='lower right')
    call tight_layout()
    call savefig('output/general-transport-krook-limit.png')
    call write_result()
contains
    subroutine write_result()
        character(len=*), parameter :: bound = 'certified-continuous linear plot '// &
            'band (sup residual inverse + chord remainder)'
        character(len=*), parameter :: tangent = 'certified-continuous nu-tangent '// &
            'band; adjoint goal gradient unavailable'
        call result%create(run)
        call result%write('x', x, '1', 'plot points', dims='x')
        call result%write('source_lower', plot_source_lo, '1', 'supplied source g lower', &
            dims='x', status=bound)
        call result%write('source_upper', plot_source_hi, '1', 'supplied source g upper', &
            dims='x', status=bound)
        call result%write('coarse_f', coarse_line, '1', 'coarse f candidate', dims='x', &
            status='candidate center')
        call result%write('coarse_lower', plot_coarse_lo, '1', 'coarse f lower', &
            dims='x', status=bound)
        call result%write('coarse_upper', plot_coarse_hi, '1', 'coarse f upper', &
            dims='x', status=bound)
        call result%write('fine_f', fine_line, '1', 'fine f candidate', dims='x', &
            status='candidate center')
        call result%write('fine_lower', plot_fine_lo, '1', 'fine f lower', dims='x', &
            status=bound)
        call result%write('fine_upper', plot_fine_hi, '1', 'fine f upper', dims='x', &
            status=bound)
        call result%write('coarse_dnu_lower', plot_coarse_dlo, '1', &
            'coarse df/dnu lower', dims='x', status=tangent)
        call result%write('coarse_dnu_upper', plot_coarse_dhi, '1', &
            'coarse df/dnu upper', dims='x', status=tangent)
        call result%write('fine_dnu_lower', plot_fine_dlo, '1', 'fine df/dnu lower', &
            dims='x', status=tangent)
        call result%write('fine_dnu_upper', plot_fine_dhi, '1', 'fine df/dnu upper', &
            dims='x', status=tangent)
        call result%write('krook_limit_lower', limit_lo, '1', &
            'Krook limiting slice f(x, v=0) lower', dims='x', status=bound)
        call result%write('krook_limit_upper', limit_hi, '1', &
            'Krook limiting slice f(x, v=0) upper', dims='x', status=bound)
        call result%write('fine_state_sup_radius', fine%certificate%state_radius, '1', &
            'fine state sup radius', status='certified-continuous')
        call result%write('fine_tangent_sup_radius', fine%certificate%tangent_radius, &
            '1', 'fine nu-tangent sup radius', status='certified-continuous')
        call result%close()
    end subroutine

    subroutine expand_band(lower, upper, step_lower, step_upper)
        real(dp), intent(in) :: lower(:), upper(:)
        real(dp), intent(out) :: step_lower(:), step_upper(:)
        integer :: j
        do j = 1, size(lower)
            step_lower(2*j - 1:2*j) = lower(j)
            step_upper(2*j - 1:2*j) = upper(j)
        end do
    end subroutine
    subroutine prepare_case(problem, cells)
        type(case_t), intent(out) :: problem
        integer, intent(in) :: cells
        integer :: j
        logical :: prepared

        call initialize_transport(problem%state, cells, modes, 1.0_dp, prepared)
        if (.not. prepared) error stop 'transport initialization failed'
        allocate (problem%dsource_cos(0:modes, cells + 1), &
                  problem%dsource_sin(0:modes, cells + 1), &
                  problem%dcos(0:modes, cells + 1), problem%dsin(0:modes, cells + 1))
        problem%dsource_cos = 0.0_dp; problem%dsource_sin = 0.0_dp
        do j = 1, cells + 1
            call transport_source_affine(problem%state%v(j), 1.0_dp, 0.125_dp, &
                problem%state%source_cos(0, j))
            call transport_source_affine(problem%state%v(j), 0.25_dp, -0.25_dp, &
                problem%state%source_cos(1, j))
            call transport_source_affine(problem%state%v(j), 0.0_dp, 0.125_dp, &
                problem%state%source_sin(2, j))
        end do
        call solve_transport(problem%state, prepared)
        if (.not. prepared) error stop 'transport solve failed'
        call solve_transport_tangent(problem%state, 1.0_dp, problem%dsource_cos, &
            problem%dsource_sin, problem%dcos, problem%dsin, prepared)
        if (.not. prepared) error stop 'transport nu tangent solve failed'
        call certify_transport_tangent(problem%state, 1.0_dp, problem%dsource_cos, &
            problem%dsource_sin, problem%dcos, problem%dsin, &
            problem%certificate, prepared)
        if (.not. prepared) error stop 'transport state/tangent certification failed'
        if (.not. transport_certificate_accepts(problem%state, &
            problem%certificate%state_radius)) &
            error stop 'transport certificate rejected'
        print '(a,i0,a,es24.16,a,es24.16)', 'cells=', cells, &
            ' state_sup_radius=', problem%certificate%state_radius, &
            ' nu_tangent_sup_radius=', problem%certificate%tangent_radius
        call report_goal(problem%state, problem%state%cosine, &
            problem%certificate%state_radius, goal_scale, 'Q normalized cosine-density')
        call report_goal(problem%state, problem%dcos, &
            problem%certificate%tangent_radius, derivative_scale, 'dQ/dnu fixed source')
    end subroutine

    subroutine report_goal(state, coefficients, radius, scale, label)
        type(transport_state_t), intent(in) :: state
        real(dp), intent(in) :: coefficients(0:, :), radius, scale
        character(len=*), intent(in) :: label
        type(interval_t) :: total, local, updated, relative, target, goal_radius
        real(dp) :: goal_lower, goal_upper
        integer :: j
        logical :: enclosed

        total = interval(0.0_dp)
        do j = 1, size(state%v) - 1
            call transport_cosine_moment_cell_interval(interval(state%v(j)), &
                interval(state%v(j + 1)), interval(coefficients(1, j)), &
                interval(coefficients(1, j + 1)), local)
            call transport_scalar_add_interval(total, local, updated)
            total = updated
        end do
        call transport_observable_band(total, interval(2.0_dp), radius, &
            goal_lower, goal_upper, enclosed)
        if (.not. enclosed) error stop 'cosine moment enclosure failed'
        goal_radius = interval(rad(interval(goal_lower, goal_upper)))
        call transport_accuracy_interval(goal_radius, interval(scale), &
            interval(rtol), interval(atol), relative, target)
        print '(a,a,2es24.16)', trim(label), ' interval=', goal_lower, goal_upper
        print '(a,es24.16,a,es24.16,a,l1)', 'declared_quantity_scale=', scale, &
            ' target_lower=', target%lo, &
            ' requested_radius_met=', goal_radius%hi <= target%lo
    end subroutine

    subroutine envelope(horizontal, lower, upper, color)
        real(dp), intent(in) :: horizontal(:), lower(:), upper(:), color(3)

        call fill_between(horizontal, lower, upper, color=color, alpha=0.25_dp)
    end subroutine
end program transport_demo
