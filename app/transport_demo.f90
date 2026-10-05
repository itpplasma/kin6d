program transport_demo
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_transport_laboratory, only: transport_state_t, initialize_transport, &
        solve_transport, solve_transport_tangent, sample_transport, &
        sample_transport_arrays
    use kin6d_transport_certificate, only: transport_certificate_t, &
        certify_transport, certify_transport_tangent, transport_certificate_accepts, &
        certify_transport_sample, certify_transport_tangent_sample, &
        transport_observable_band
    use kin6d_transport_sampling, only: transport_basis_t, &
        transport_sampling_workspace_t, initialize_transport_sampling, &
        prepare_transport_basis, sample_transport_input_cached, &
        certify_transport_cached, certify_transport_tangent_cached
    use kin6d_transport_leaves, only: transport_source_affine, transport_plot_x, &
        transport_cosine_moment_cell_interval, transport_scalar_add_interval, &
        transport_accuracy_interval
    use kin6d_poisson_leaves, only: coupled_distribution_interval
    use fortnum_interval, only: interval_t, interval, rad
    use fortplot, only: figure, subplot, plot, fill_between, title, xlabel, &
        ylabel, legend, tight_layout, savefig, ylim
    use fortplot_system_runtime, only: create_directory_runtime
    implicit none
    ! Subpixel plotting cells retain full enclosures without visible stairs.
    ! This display resolution is independent of the solver and benchmark grids.
    integer, parameter :: modes = 2, nx = 1536
    real(dp), parameter :: velocity = 0.3125_dp, rtol = 0.01_dp, atol = 1.0e-8_dp
    real(dp), parameter :: goal_scale = 0.25_dp, derivative_scale = 0.25_dp
    type :: case_t
        type(transport_state_t) :: state
        type(transport_certificate_t) :: certificate
        real(dp), allocatable :: dsource_cos(:, :), dsource_sin(:, :)
        real(dp), allocatable :: dcos(:, :), dsin(:, :)
    end type
    type(case_t) :: coarse, fine
    type(transport_state_t) :: limiting
    type(transport_certificate_t) :: limit_certificate, audit
    type(transport_basis_t) :: basis
    type(transport_sampling_workspace_t) :: work
    type(interval_t) :: cell, reference, x_cells(nx), sources(nx)
    real(dp) :: band_lo(nx), band_hi(nx)
    real(dp) :: x(nx + 1), edges(2*nx), source_line(nx + 1)
    real(dp) :: coarse_line(nx + 1), fine_line(nx + 1)
    real(dp) :: coarse_tangent(nx + 1), fine_tangent(nx + 1)
    real(dp) :: source_lo(2*nx), source_hi(2*nx)
    real(dp) :: coarse_lo(2*nx), coarse_hi(2*nx), fine_lo(2*nx), fine_hi(2*nx)
    real(dp) :: coarse_dlo(2*nx), coarse_dhi(2*nx)
    real(dp) :: fine_dlo(2*nx), fine_dhi(2*nx)
    real(dp) :: limit_line(nx + 1), limit_lo(2*nx), limit_hi(2*nx)
    real(dp) :: reference_lo(2*nx), reference_hi(2*nx), reference_line(nx + 1)
    real(dp) :: lower, upper
    integer :: i, output
    logical :: ok

    call create_directory_runtime('output', ok)
    if (.not. ok) error stop 'cannot create output directory'
    call prepare_case(coarse, 16)
    call prepare_case(fine, 64)
    call certify_transport(fine%state, audit, ok, strict=.true.)
    if (.not. ok) error stop 'strict-box transport audit failed'
    print '(a,es24.16)', 'strict_box_audit_state_radius=', audit%state_radius
    print '(a)', 'certificate_kind=continuous_sup_residual_inverse'
    print '(a)', 'source=exact_supplied_Fourier_x_P1_v; nu=1; E=0; v_domain=[-1,1]'
    print '(a)', 'nu_tangent=certified_continuous; adjoint_goal_gradient=unavailable'

    do i = 1, nx + 1
        call transport_plot_x(real(i - 1, dp), real(nx, dp), x(i))
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
        edges(2*i - 1) = x(i); edges(2*i) = x(i + 1)
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

    call figure(figsize=[9.6_dp, 7.2_dp])
    call subplot(2, 1, 1)
    call title('Unknown-solution damped transport | full x-cell enclosures')
    call xlabel('x / dimensionless periodic coordinate (period 2 pi)')
    call ylabel('f and supplied source g at v = 5/16')
    call envelope(edges, source_lo, source_hi, 'gainsboro', [0.4_dp, 0.4_dp, 0.4_dp])
    call envelope(edges, coarse_lo, coarse_hi, 'moccasin', &
        [0.8353_dp, 0.3686_dp, 0.0_dp])
    call envelope(edges, fine_lo, fine_hi, 'lightblue', &
        [0.0_dp, 0.4471_dp, 0.6980_dp])
    call plot(x, source_line, label='supplied source + cell range', linestyle=':', &
        color=[0.4_dp, 0.4_dp, 0.4_dp])
    call plot(x, coarse_line, label='16 cells + full band', linestyle='--', &
        color=[0.8353_dp, 0.3686_dp, 0.0_dp])
    call plot(x, fine_line, label='64 cells + full band', &
        color=[0.0_dp, 0.4471_dp, 0.6980_dp])
    call legend(loc='lower right')
    call subplot(2, 1, 2)
    call title('Implicit nu tangent | fixed independently supplied source')
    call xlabel('x / dimensionless periodic coordinate (period 2 pi)')
    call ylabel('d f / d nu at v = 5/16')
    call envelope(edges, coarse_dlo, coarse_dhi, 'moccasin', &
        [0.8353_dp, 0.3686_dp, 0.0_dp])
    call envelope(edges, fine_dlo, fine_dhi, 'lightblue', &
        [0.0_dp, 0.4471_dp, 0.6980_dp])
    call plot(x, coarse_tangent, label='16 cells + full tangent band', linestyle='--', &
        color=[0.8353_dp, 0.3686_dp, 0.0_dp])
    call plot(x, fine_tangent, label='64 cells + full tangent band', &
        color=[0.0_dp, 0.4471_dp, 0.6980_dp])
    call legend(loc='lower right')
    call tight_layout()
    call savefig('output/general-transport.png')

    ! The earlier reference is used only in this separately declared slice.
    ! Its full velocity source differs; no reference enters either checker.
    call initialize_transport(limiting, 16, 0, 1.0_dp, ok)
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
    do i = 1, nx
        cell = interval(x(i), x(i + 1))
        call certify_transport_sample(limiting, limit_certificate, cell, &
            interval(0.0_dp), lower, upper, ok)
        if (.not. ok) error stop 'limiting cell certificate failed'
        limit_lo(2*i - 1:2*i) = lower; limit_hi(2*i - 1:2*i) = upper
        call coupled_distribution_interval(interval(0.0_dp), interval(1.0_dp), &
            cell, interval(0.0_dp), reference)
        reference_lo(2*i - 1:2*i) = reference%lo
        reference_hi(2*i - 1:2*i) = reference%hi
    end do
    call figure(figsize=[9.6_dp, 4.8_dp])
    call title('Field-free zero-mode slice | same g = 3/4, nu = 1 at v = 0')
    call xlabel('x / dimensionless periodic coordinate (period 2 pi)')
    call ylabel('f(x, v=0) / dimensionless distribution')
    call ylim(0.70_dp, 0.80_dp)
    call envelope(edges, limit_lo, limit_hi, 'lightblue', [0.0_dp, 0.45_dp, 0.7_dp])
    call envelope(edges, reference_lo, reference_hi, 'gainsboro', &
        [0.0_dp, 0.0_dp, 0.0_dp])
    call plot(x, limit_line, label='general solver + full slice band', &
        color=[0.0_dp, 0.45_dp, 0.7_dp])
    call plot(x, reference_line, label='earlier Krook reference + full slice band', &
        linestyle='--', color=[0.0_dp, 0.0_dp, 0.0_dp])
    call legend(loc='upper right')
    call tight_layout()
    call savefig('output/general-transport-krook-limit.png')
contains
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

    subroutine envelope(horizontal, lower, upper, fill, color)
        real(dp), intent(in) :: horizontal(:), lower(:), upper(:), color(3)
        character(len=*), intent(in) :: fill

        call fill_between(horizontal, lower, upper, color=fill, alpha=0.55_dp)
        call plot(horizontal, lower, linestyle=':', color=color)
        call plot(horizontal, upper, linestyle=':', color=color)
    end subroutine
end program transport_demo
