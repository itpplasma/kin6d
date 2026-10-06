program shielding_dynamics_demo
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_interval, only: interval_t, interval, linear_interpolation_remainder
    use kin6d_shielding_dynamics, only: shielding_state_t, init_shielding, &
                             advance_shielding, sample_shielding, shielding_plot_padding
    use kin6d_shielding_dynamics_leaves, only: shielding_coordinate, shielding_time, &
                       shielding_mode, shielding_equilibrium, shielding_band_interval, &
        shielding_time_curvature_interval, shielding_radius_width_interval, &
        shielding_interval_length_interval
    use fortplot, only: figure, get_global_figure, figure_t, subplot, plot, &
          fill_between, title, xlabel, ylabel, xlim, ylim, legend, tight_layout, savefig, &
          suptitle, set_yscale
    use kin6d_demo_style, only: blue, black, grey, movie_fps, hold_frames, num
    use fortplot_animation, only: animation_t, animation_func => FuncAnimation, &
                                  save_animation
    use kin6d_result_io, only: run_t, result_file_t
    implicit none
    integer, parameter :: modes = 6
    ! 193 frames span t = 0..4 in 8 s at 24 fps; a one-second hold follows.
    integer :: points = 193, frames = 193
    real(dp) :: source(modes) = &
                           [0.04_dp, 0.025_dp, 0.0125_dp, 0.006_dp, 0.003_dp, 0.0015_dp]
    real(dp) :: tmax = 4.0_dp, step = 0.0025_dp
    namelist /shielding_dynamics/ points, frames, source, tmax, step
    type(run_t) :: run
    type(result_file_t) :: result
    integer :: ios, center
    character(len=256) :: msg
    character(len=*), parameter :: bound = 'certified-continuous finite-Fourier '// &
        'residual/contractivity band; physical/model error unavailable'
    type(shielding_state_t) :: state, baseline
    type(figure_t), pointer :: canvas
    type(animation_t) :: animation
    real(dp), allocatable :: x(:), times(:), values(:, :, :)
    real(dp), allocatable :: lower(:, :, :), upper(:, :, :)
    real(dp), allocatable :: equilibrium(:, :), bare(:, :)
    real(dp) :: a, d, p, e, padding(3)
    real(dp) :: first, last, solve_time, certified_time, render_first, render_last
    real(dp), allocatable :: candidate_frame_s(:), certified_frame_s(:)
    real(dp) :: frame_first, frame_last, max_x_gap, max_t_gap
    type(interval_t) :: padded_lo, padded_hi, temporal_curvature, next_bound
    real(dp), allocatable :: history_lower(:), history_upper(:)
    real(dp) :: max_radius(3)
    integer :: f, j, k, q, unit, status, snapshot
    logical :: ok, snapshot_only
    character(128) :: argument
    call run%start('shielding_dynamics_demo')
    ios = 0; msg = ''
    call run%open_config(unit)
    if (run%has_config()) read (unit, nml=shielding_dynamics, iostat=ios, iomsg=msg)
    call run%close_config(unit, 'shielding_dynamics', ios, msg)
    if (points < 3 .or. mod(points, 2) == 0 .or. frames < 2) &
        error stop 'shielding_dynamics needs odd points >= 3 and frames >= 2'
    if (.not. (tmax > 0.0_dp .and. step > 0.0_dp)) &
        error stop 'shielding_dynamics tmax and step must be positive'
    call run%open_resolved(unit)
    write (unit, nml=shielding_dynamics)
    call run%close_resolved(unit)
    center = (points + 1)/2
    allocate (x(points), times(frames), values(3, points, frames), &
        lower(3, points, frames), upper(3, points, frames), equilibrium(3, points), &
        bare(3, points), candidate_frame_s(frames), certified_frame_s(frames), &
        history_lower(frames), history_upper(frames))
    snapshot_only = .false.
    snapshot = frames
    if (run%argument_count() > 2) error stop 'use [--config FILE] [--snapshot [frame]]'
    if (run%argument_count() > 0) then
        call run%argument(1, argument)
        if (trim(argument) /= '--snapshot') error stop 'use --snapshot [frame]'
        snapshot_only = .true.
        if (run%argument_count() > 1) then
            call run%argument(2, argument)
            read (argument, *, iostat=status) snapshot
            if (status /= 0) error stop 'invalid snapshot frame'
            if (snapshot < 1 .or. snapshot > frames) error stop 'snapshot frame range'
        end if
    end if
    call execute_command_line('mkdir -p output', exitstat=status)
    if (status /= 0) error stop 'output directory failed'
    do j = 1, points
        call shielding_coordinate(real(j - 1, dp), real(points, dp), x(j))
    end do
    do f = 1, frames
        call shielding_time(real(f - 1, dp), real(frames, dp), tmax, times(f))
    end do
    max_x_gap = 0.0_dp
    max_t_gap = 0.0_dp
    do j = 1, points - 1
        call shielding_interval_length_interval(interval(x(j)), interval(x(j + 1)), &
            next_bound)
        max_x_gap = max(max_x_gap, next_bound%hi)
    end do
    do f = 1, frames - 1
        call shielding_interval_length_interval(interval(times(f)), &
            interval(times(f + 1)), next_bound)
        max_t_gap = max(max_t_gap, next_bound%hi)
    end do
    bare = 0.0_dp
    equilibrium = 0.0_dp
    do j = 1, points
        do k = 1, modes
            call shielding_mode(0.0_dp, real(k, dp), source(k), x(j), d, p, e)
            bare(:, j) = bare(:, j) + [d, p, e]
            call shielding_equilibrium(source(k), real(k, dp), a)
            call shielding_mode(a, real(k, dp), source(k), x(j), d, p, e)
            equilibrium(:, j) = equilibrium(:, j) + [d, p, e]
        end do
    end do
    call cpu_time(first)
    call init_shielding(baseline, source, ok)
    if (.not. ok) error stop 'candidate initialization failed'
    do f = 1, frames
        call cpu_time(frame_first)
        call advance_shielding(baseline, times(f), step, ok, certify=.false.)
        call cpu_time(frame_last)
        candidate_frame_s(f) = frame_last - frame_first
        if (.not. ok) error stop 'candidate advance failed'
    end do
    call cpu_time(last)
    solve_time = last - first
    call cpu_time(first)
    call init_shielding(state, source, ok)
    if (.not. ok) error stop 'checked initialization failed'
    do f = 1, frames
        call cpu_time(frame_first)
        call advance_shielding(state, times(f), step, ok)
        call cpu_time(frame_last)
        certified_frame_s(f) = frame_last - frame_first
        if (.not. ok) error stop 'checked advance failed'
    end do
    call cpu_time(last)
    certified_time = last - first
    open (newunit=unit, file='output/shielding-performance.csv', status='replace')
write (unit, '(a)') '# matched local CPU single trial; step solves only; setup excluded'
  write (unit, '(a)') 'frame,t,candidate_s,certified_s,extra_s,certified_over_candidate'
    do f = 1, frames
        if (candidate_frame_s(f) > 0.0_dp) then
            write (unit, '(i0,5(",",es24.16e3))') f, times(f), candidate_frame_s(f), &
                certified_frame_s(f), certified_frame_s(f) - candidate_frame_s(f), &
                certified_frame_s(f)/candidate_frame_s(f)
        else
            write (unit, '(i0,4(",",es24.16e3),a)') f, times(f), candidate_frame_s(f), &
       certified_frame_s(f), certified_frame_s(f) - candidate_frame_s(f), ',unavailable'
        end if
    end do
    close (unit)
    call init_shielding(state, source, ok)
    open (newunit=unit, file='output/shielding-dynamics.csv', status='replace')
 write (unit, '(a)') '# periodic planar linear drift-diffusion-Poisson, x/lambda_D,'// &
        ' t/(lambda_D**2/D); finite supplied cosine source; eta(t=0)=0'
 write(unit,'(a)') '# continuous numerical bands; physical/model error unavailable;'// &
        ' smooth leaf JVP/VJP generated; continuous derivative bounds unavailable'
    write (unit, '(a,es24.16e3)') '# candidate_init_solve_seconds=', solve_time
    write (unit, '(a,es24.16e3)') '# certified_init_solve_seconds=', certified_time
    write (unit, '(a,es24.16e3)') '# extra_seconds=', certified_time - solve_time
    if (solve_time > 0.0_dp) write (unit, '(a,es24.16e3)') &
        '# certified_over_candidate_factor=', certified_time/solve_time
    write (unit, '(a)') '# single local CPU trial; same frame endpoints and setup;'// &
        ' sampling/rendering excluded from both solve costs'
    write (unit, '(a)') 'frame,t,x,eta,eta_lower,eta_upper,phi,phi_lower,phi_upper,'// &
        'E,E_lower,E_upper'
    max_radius = 0.0_dp
    call cpu_time(render_first)
    do f = 1, frames
        call advance_shielding(state, times(f), step, ok)
        if (.not. ok) error stop 'trajectory advance failed'
        call shielding_plot_padding(state, 0.0_dp, max_x_gap, padding, ok)
        if (.not. ok) error stop 'display chord bound failed'
        do j = 1, points
            call sample_shielding(state, x(j), values(:, j, f), lower(:, j, f), &
                                  upper(:, j, f), ok)
            if (.not. ok) error stop 'pointwise certificate failed'
            do q = 1, 3
                call shielding_band_interval(interval(lower(q, j, f)), &
                                             interval(padding(q)), padded_lo, padded_hi)
                lower(q, j, f) = padded_lo%lo
                call shielding_band_interval(interval(upper(q, j, f)), &
                                             interval(padding(q)), padded_lo, padded_hi)
                upper(q, j, f) = padded_hi%hi
                call shielding_radius_width_interval(interval(lower(q, j, f)), &
                                                   interval(upper(q, j, f)), next_bound)
                max_radius(q) = max(max_radius(q), next_bound%hi)
            end do
            write (unit, '(i0,11(",",es24.16e3))') f, times(f), x(j), &
                (values(q, j, f), lower(q, j, f), upper(q, j, f), q=1, 3)
        end do
    end do
    temporal_curvature = interval(0.0_dp)
    do k = 1, modes
        call shielding_time_curvature_interval(temporal_curvature, &
                                 interval(source(k)), interval(real(k, dp)), next_bound)
        temporal_curvature = next_bound
    end do
    ! The derivative itself obeys a contractive modal ODE, hence its second
    ! derivative is bounded without using the exact time solution as checker.
    next_bound = linear_interpolation_remainder(interval(0.0_dp), &
        interval(max_t_gap), temporal_curvature)
    do f = 1, frames
        call shielding_band_interval(interval(lower(2, center, f)), next_bound, &
                                     padded_lo, padded_hi)
        history_lower(f) = padded_lo%lo
        call shielding_band_interval(interval(upper(2, center, f)), next_bound, &
                                     padded_lo, padded_hi)
        history_upper(f) = padded_hi%hi
    end do
write(unit,'(a)') '# observable_scales=0.05,0.05,0.05; requested_rtol=0.01; floor=1e-12'
    write (unit, '(a,3(es24.16e3,1x))') '# max_spatial_band_radii=', max_radius
    if (any(max_radius > 0.000500000001_dp)) &
        error stop 'requested one-percent band radius not achieved'
    call cpu_time(render_last)
    write (unit, '(a,es24.16e3)') '# checked_solve_sampling_seconds=', &
        render_last - render_first
    close (unit)
    call result%create(run)
    call result%write('t', times, 'lambda_D**2/D', 'time', dims='frame')
    call result%write('x', x, 'lambda_D', 'position', dims='x')
    call write_field(1, 'eta', 'electron excess')
    call write_field(2, 'phi', 'potential')
    call write_field(3, 'E', 'electric field')
    call result%write('phi_center_history_lower', history_lower, '1', &
        'center potential history lower (time-interpolation padded)', dims='frame', &
        status=bound)
    call result%write('phi_center_history_upper', history_upper, '1', &
        'center potential history upper (time-interpolation padded)', dims='frame', &
        status=bound)
    call result%write('max_band_radius', max_radius, '1', &
        'max spatial band radius of eta, phi, E', dims='field', status=bound)
    call result%close()
    call figure(figsize=[13.2_dp, 8.8_dp])
    canvas => get_global_figure()
    call draw_frame(snapshot)
    call savefig('output/debye-shielding.png')
    if (.not. snapshot_only) then
        animation = animation_func(update_frame, frames=frames + hold_frames, &
            interval=42, fig=canvas)
        call save_animation(animation, 'output/debye-shielding.mp4', fps=movie_fps, &
                            status=status)
        if (status /= 0) error stop 'movie encoding failed'
    end if
contains
    subroutine write_field(q, name, long_name)
        integer, intent(in) :: q
        character(len=*), intent(in) :: name, long_name
        call result%write(name, values(q, :, :), '1', long_name//' center', &
            dims='x,frame', status='candidate center')
        call result%write(name//'_lower', lower(q, :, :), '1', long_name//' lower', &
            dims='x,frame', status=bound)
        call result%write(name//'_upper', upper(q, :, :), '1', long_name//' upper', &
            dims='x,frame', status=bound)
    end subroutine
    subroutine update_frame(number)
        integer, intent(in) :: number
        call canvas%clear()
        call draw_frame(min(number, frames))
        call canvas%set_rendered(.false.)
    end subroutine
    subroutine draw_frame(number)
        integer, intent(in) :: number
        real(dp), parameter :: floor = 1.0e-6_dp  ! bottom of the log axis
        integer :: panel
        call suptitle('Electrons pile up on a fixed source and screen its potential, '// &
            'relaxing to the linear Debye equilibrium')
        do panel = 1, 3
            call subplot(2, 2, panel)
            call xlabel('$x/\lambda_D$  (periodic planar model)')
            call xlim(x(1), x(points))
            select case (panel)
            case (1)
                call title('Electron excess builds up  |  t = '//num(times(number), 2))
                call ylabel('$\eta = \delta n_e/n_0$')
                call ylim(-0.034_dp, 0.05_dp)
            case (2)
                call title('Potential drops as the excess screens the source')
                call ylabel('$e\phi/k_BT$')
                call ylim(-0.04_dp, 0.08_dp)
                call plot(x, bare(2, :), linestyle=':', color=grey, linewidth=1.5_dp, &
                          label='bare source (t = 0)')
            case (3)
                call title('Electric field is reduced, not replaced by a Yukawa force')
                call ylabel('$e\lambda_D E/k_BT$')
                call ylim(-0.058_dp, 0.082_dp)
                call plot(x, bare(3, :), linestyle=':', color=grey, linewidth=1.5_dp, &
                          label='bare source (t = 0)')
            end select
            call fill_between(x, lower(panel, :, number), upper(panel, :, number), &
                              color=blue, alpha=0.25_dp)
            call plot(x, values(panel, :, number), color=blue, linewidth=2.0_dp, &
                      label='numerical ± certified band')
            call plot(x, equilibrium(panel, :), linestyle='--', color=black, &
                      linewidth=1.2_dp, label='Debye equilibrium')
            call legend(loc='upper right')
        end do
        call subplot(2, 2, 4)
        call title('Exponential approach to equilibrium at x = 0 (band resolved)')
        call xlabel('$tD/\lambda_D^2$')
        call ylabel('$e|\phi(0,t) - \phi_{eq}(0)|/k_BT$')
        call xlim(0.0_dp, tmax)
        call ylim(1.0e-6_dp, 0.1_dp)
        call set_yscale('log')
        ! Distance of the certified band from the equilibrium value; a band that
        ! contains the equilibrium reaches the axis floor.
        if (number >= 2) then
            call fill_between(times(:number), &
                merge(floor, max(floor, min(abs(history_lower(:number) - &
                equilibrium(2, center)), abs(history_upper(:number) - &
                equilibrium(2, center)))), history_lower(:number) <= equilibrium(2, center) &
                .and. history_upper(:number) >= equilibrium(2, center)), &
                max(floor, max(abs(history_lower(:number) - equilibrium(2, center)), &
                abs(history_upper(:number) - equilibrium(2, center)))), &
                color=blue, alpha=0.25_dp)
            call plot(times(:number), max(floor, abs(values(2, center, :number) - &
                equilibrium(2, center))), color=blue, linewidth=2.0_dp, &
                label='numerical ± certified time band')
        end if
        call plot([times(number)], [max(floor, abs(values(2, center, number) - &
            equilibrium(2, center)))], color=blue, linestyle='none', marker='o', &
            markersize=6.0_dp)
        if (number >= 2) call legend(loc='upper right')
        call tight_layout()
    end subroutine
end program shielding_dynamics_demo
