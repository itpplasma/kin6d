program streaming_comparison_demo
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_interval, only: interval_t, interval
    use kin6d_streaming_comparison, only: streaming_state_t, streaming_particles_t, &
        init_streaming, advance_streaming, evaluate_streaming, streaming_radius, &
        evaluate_streaming_candidate_interval, &
        init_streaming_particles, get_streaming_particles, density_streaming, &
        density_streaming_chord
    use kin6d_streaming_particle_density, only: particle_density_certificate_t, &
        particle_density_workspace_t, prepare_particle_density_workspace, &
        project_particle_density, certify_particle_density, &
        sample_particle_density_band, particle_density_chord
    use kin6d_streaming_leaves, only: streaming_plot_x, streaming_velocity_coordinate, &
        streaming_frame_time, streaming_characteristic, streaming_display_shift, &
        streaming_particle_midpoint, streaming_particle_glyph_size, &
        streaming_mode_accumulate, streaming_enlarge_interval
    use fortplot, only: figure, get_global_figure, figure_t, subplot, scatter, &
        pcolormesh, plot, fill_between, title, xlabel, ylabel, xlim, ylim, &
        legend, tight_layout, savefig
    use fortplot_animation, only: animation_t, animation_func => FuncAnimation, &
        save_animation
    use fortplot_system_runtime, only: create_directory_runtime
    use kin6d_result_io, only: run_t, result_file_t
    implicit none
    integer, parameter :: modes = 1, fps = 24, cohort_size = 3
    integer :: nx = 96, nv = 128, particle_nx = 32, particle_nv = 48
    integer :: frames = 121, repeats = 20, np
    real(dp) :: duration = 12.0_dp
    namelist /streaming_comparison/ nx, nv, particle_nx, particle_nv, frames, &
        repeats, duration
    type(run_t) :: run
    type(result_file_t) :: result
    type(streaming_state_t) :: eulerian, initial, particle_initial
    type(streaming_particles_t) :: particles
    type(particle_density_certificate_t) :: particle_certificate
    type(particle_density_workspace_t) :: particle_work
    type(figure_t), pointer :: canvas
    type(animation_t) :: animation
    real(dp), allocatable :: v(:), x(:), cosine0(:, :), sine0(:, :)
    real(dp) :: initial_cosine(0:modes) = [1.0_dp, 0.5_dp]
    real(dp) :: particle_cosine0(0:modes, 2), particle_sine0(0:modes, 2)
    real(dp) :: particle_v0(2) = [-1.0_dp, 1.0_dp]
    real(dp), allocatable :: x0(:), particle_v(:), weight(:), unwrapped(:)
    real(dp), allocatable :: particle_x(:), glyph(:)
    real(dp) :: cohort_x(cohort_size), cohort_v(cohort_size)
    integer :: cohorts(cohort_size, 3)
    type(interval_t), allocatable :: x0_bounds(:), velocity_bounds(:), weight_bounds(:)
    real(dp), allocatable :: image(:, :), slice(:), slice_lo(:), slice_hi(:)
    real(dp), allocatable :: initial_slice(:), initial_lo(:), initial_hi(:)
    real(dp), allocatable :: density(:), density_lo(:), density_hi(:)
    real(dp), allocatable :: particle_density(:), particle_lo(:), particle_hi(:)
    real(dp), allocatable :: density_padding(:), particle_padding(:)
    ! Per-frame evidence retained for the NetCDF result.
    real(dp), allocatable :: frame_time(:), frame_radius(:), store_density(:, :, :)
    real(dp), allocatable :: store_particle(:, :, :), store_slice(:, :, :)
    real(dp) :: pcos(0:modes), psin(0:modes), time, radius, period, measure
    real(dp) :: candidate_cpu, complete_cpu, particle_push_cpu, particle_project_cpu
    real(dp) :: particle_complete_cpu, postprocess_cpu, ignored_x, ignored_v
    integer :: i, frame, status, evidence, bands, slices, particle_data, snapshot_frame
    integer :: unit
    logical :: ok, snapshot_only
    character(len=64) :: argument
    character(len=256) :: msg

    call run%start('streaming_comparison_demo')
    status = 0; msg = ''
    call run%open_config(unit)
    if (run%has_config()) read (unit, nml=streaming_comparison, iostat=status, &
        iomsg=msg)
    call run%close_config(unit, 'streaming_comparison', status, msg)
    if (nx < 1 .or. nv < 1 .or. particle_nx < 1 .or. particle_nv < 1) &
        error stop 'grid sizes must be positive'
    if (frames < 2 .or. repeats < 1 .or. .not. (duration > 0.0_dp)) &
        error stop 'frames >= 2, repeats >= 1 and duration > 0 required'
    call run%open_resolved(unit)
    write (unit, nml=streaming_comparison)
    call run%close_resolved(unit)
    np = particle_nx*particle_nv
    allocate (v(nv + 1), x(nx + 1), cosine0(0:modes, nv + 1), sine0(0:modes, nv + 1), &
        x0(np), particle_v(np), weight(np), unwrapped(np), particle_x(np), glyph(np), &
        x0_bounds(np), velocity_bounds(np), weight_bounds(np), &
        image(nv + 1, nx + 1), slice(nv + 1), slice_lo(nv + 1), slice_hi(nv + 1), &
        initial_slice(nv + 1), initial_lo(nv + 1), initial_hi(nv + 1), &
        density(nx + 1), density_lo(nx + 1), density_hi(nx + 1), &
        particle_density(nx + 1), particle_lo(nx + 1), particle_hi(nx + 1), &
        density_padding(nx), particle_padding(nx), frame_time(frames), &
        frame_radius(frames), store_density(nx + 1, 2, frames), &
        store_particle(nx + 1, 2, frames), store_slice(nv + 1, 2, frames))
    snapshot_only = .false.; snapshot_frame = (frames + 1)/2
    if (run%argument_count() > 2) error stop &
        'usage: demo [--config FILE] [--snapshot [frame]]'
    if (run%argument_count() >= 1) then
        call run%argument(1, argument)
        if (trim(argument) /= '--snapshot') error stop 'expected --snapshot'
        snapshot_only = .true.
    end if
    if (run%argument_count() == 2) then
        call run%argument(2, argument)
        read (argument, *, iostat=status) snapshot_frame
        if (status /= 0) error stop 'snapshot frame must be an integer'
    end if
    if (snapshot_frame < 1 .or. snapshot_frame > frames) &
        error stop 'snapshot frame outside 1..frames'
    call create_directory_runtime('output', ok)
    if (.not. ok) error stop 'cannot create output directory'
    do i = 1, nv + 1
        call streaming_velocity_coordinate(real(i - 1, dp), real(nv, dp), v(i))
        cosine0(:, i) = initial_cosine
    end do
    sine0 = 0.0_dp
    do i = 1, nx + 1
        call streaming_plot_x(real(i - 1, dp), real(nx + 1, dp), x(i))
    end do
    period = x(nx + 1)
    particle_cosine0(:, 1) = initial_cosine
    particle_cosine0(:, 2) = initial_cosine
    particle_sine0 = 0.0_dp
    call init_streaming(eulerian, v, cosine0, sine0, ok)
    if (.not. ok) error stop 'Eulerian initialization failed'
    call init_streaming(initial, v, cosine0, sine0, ok)
    if (.not. ok) error stop 'initial slice initialization failed'
    call init_streaming(particle_initial, particle_v0, particle_cosine0, &
        particle_sine0, ok)
    if (.not. ok) error stop 'particle source initialization failed'
    call init_streaming_particles(particle_initial, particle_nx, particle_nv, &
        particles, ok)
    if (.not. ok) error stop 'weighted particle initialization failed'
    call prepare_particle_density_workspace(particle_initial, particles, &
        particle_work, ok)
    if (.not. ok) error stop 'particle workspace preparation failed'
    call certify_particle_density(particle_initial, particles, 0.0_dp, &
        particle_certificate, ok, particle_work)
    if (.not. ok) error stop 'particle certificate initialization failed'
    call get_streaming_particles(particles, x0, particle_v, weight, ok, &
        x0_bounds, velocity_bounds, weight_bounds)
    if (.not. ok) error stop 'particle geometry retrieval failed'
    call streaming_particle_midpoint(0.0_dp, real(particle_nx, dp), 0.0_dp, &
        real(particle_nv, dp), ignored_x, ignored_v, measure)
    do i = 1, np
        call streaming_particle_glyph_size(weight(i), measure, glyph(i))
    end do
    ! Fixed identifiers selected by initial phase and three separated velocities.
    call select_cohorts()
    call write_particle_data()
    print '(a)', 'model=free_streaming; E=0; x_period=2*pi; velocity_domain=[-1,1]'
    print '(a)', 'initial=f0=1+0.5*cos(x); particle_weights=f0*dx*dv'
    print '(a)', 'Eulerian=Fourier_x_P1_v'
    print '(a)', 'particles=weighted_characteristic_midpoint_carrier'
    print '(a)', 'density=normalized_0.5_integral_f_dv; particle_curve=Fourier_PIF'
    print '(a)', 'particle_pointwise_Dirac_bound=unavailable; field_gather=absent'
    print '(a)', 'weak_density_sensitivity=unavailable; chart_wrap=display_only'
    print '(a)', 'science_timers=CPU_seconds_20_repeats; rendering_excluded'
    print '(a)', 'glyph_area=12*f0; tracked_cohorts=9_fixed_particle_IDs'
    print '(a)', 'timer_factor=-1_means_baseline_below_timer_resolution'
    open (newunit=evidence, file='output/streaming-comparison.csv', status='replace')
    write (evidence, '(a)') 'frame,t,velocity_interpolation_radius,'// &
        'Eulerian_candidate_s,Eulerian_complete_s,Eulerian_factor,'// &
        'particle_advance_s,particle_projection_s,particle_complete_s,'// &
        'particle_factor,particle_observable_factor,particle_observable_extra_s,'// &
        'Eulerian_extra_s,'// &
        'bound_postprocess_s,density_x0_lower,density_x0_upper,'// &
        'particle_density_x0_lower,particle_density_x0_upper'
    open (newunit=bands, file='output/streaming-density-bands.csv', status='replace')
    write (bands, '(a)') 'frame,t,x,Eulerian_density_lower,Eulerian_density_upper,'// &
        'particle_density_lower,particle_density_upper'
    open (newunit=slices, file='output/streaming-slice-bands.csv', status='replace')
    write (slices, '(a)') 'frame,t,v,lower,upper,initial_lower,initial_upper'
    do frame = 1, frames
        call compute_frame(frame, timed=.true.)
        call write_evidence(frame)
    end do
    close (evidence); close (bands); close (slices)
    call write_result()
    call figure(figsize=[13.2_dp, 8.8_dp])
    canvas => get_global_figure()
    call compute_frame(snapshot_frame, timed=.false.)
    call draw_frame()
    call savefig('output/particles-vlasov.png')
    if (.not. snapshot_only) then
        animation = animation_func(update_frame, frames=frames, interval=42, fig=canvas)
        call save_animation(animation, 'output/particles-vlasov.mp4', fps=fps, &
            status=status)
        if (status /= 0) error stop 'animation rendering or encoding failed'
    end if
contains
    subroutine advance_particles()
        integer :: p
        do p = 1, np
            call streaming_characteristic(x0(p), particle_v(p), time, unwrapped(p))
        end do
    end subroutine

    subroutine compute_frame(number, timed)
        integer, intent(in) :: number
        logical, intent(in) :: timed
        real(dp) :: start, finish, lower, upper, padding, next
        type(interval_t) :: enlarged_lo, enlarged_hi, point_lo, point_hi, chord_padding
        type(interval_t) :: candidate_packet
        integer :: repeat, ix, iv, p, turns

        call streaming_frame_time(real(number - 1, dp), real(frames, dp), &
            duration, time)
        if (timed) then
            call cpu_time(start)
            do repeat = 1, repeats
                call advance_streaming(eulerian, time, ok, certify=.false.)
            end do
            call cpu_time(finish)
            candidate_cpu = (finish - start)/repeats
            if (.not. ok) error stop 'uncertified Eulerian advance failed'
            call cpu_time(start)
            do repeat = 1, repeats
                call advance_streaming(eulerian, time, ok)
            end do
            call cpu_time(finish)
            complete_cpu = (finish - start)/repeats
            if (.not. ok) error stop 'certified Eulerian advance failed'
            call cpu_time(start)
            do repeat = 1, repeats
                call advance_particles()
            end do
            call cpu_time(finish)
            particle_push_cpu = (finish - start)/repeats
            call cpu_time(start)
            do repeat = 1, repeats
                call project_particle_density(particles, time, pcos, psin, &
                ok, particle_work)
            end do
            call cpu_time(finish)
            particle_project_cpu = (finish - start)/repeats
            if (.not. ok) error stop 'particle Fourier projection failed'
            call cpu_time(start)
            do repeat = 1, repeats
                call advance_particles()
                call project_particle_density(particles, time, pcos, psin, &
                ok, particle_work)
                call certify_particle_density(particle_initial, particles, time, &
                    particle_certificate, ok, particle_work)
            end do
            call cpu_time(finish)
            particle_complete_cpu = (finish - start)/repeats
            if (.not. ok) error stop 'particle weak-density checker failed'
        else
            call advance_streaming(eulerian, time, ok)
            if (.not. ok) error stop 'Eulerian animation advance failed'
            call advance_particles()
            call project_particle_density(particles, time, pcos, psin, &
                ok, particle_work)
            if (.not. ok) error stop 'particle animation projection failed'
            call certify_particle_density(particle_initial, particles, time, &
                particle_certificate, ok, particle_work)
            if (.not. ok) error stop 'particle animation checker failed'
        end if
        ! Render scheduling never changes the unwrapped scientific coordinates.
        do p = 1, np
            turns = floor(unwrapped(p)/period)
            call streaming_display_shift(unwrapped(p), real(turns, dp), particle_x(p))
        end do
        call streaming_radius(eulerian, radius, ok)
        if (.not. ok) error stop 'Eulerian velocity interpolation radius unavailable'
        call cpu_time(start)
        do ix = 1, nx + 1
            do iv = 1, nv + 1
                call evaluate_streaming(eulerian, x(ix), v(iv), image(iv, ix), &
                    lower, upper, ok)
                if (.not. ok) error stop 'Eulerian phase-space sampling failed'
            end do
        end do
        do iv = 1, nv + 1
            call evaluate_streaming(eulerian, 0.0_dp, v(iv), slice(iv), &
                slice_lo(iv), slice_hi(iv), ok)
            if (.not. ok) error stop 'Eulerian slice sampling failed'
            call evaluate_streaming_candidate_interval(eulerian, interval(0.0_dp), &
                interval(v(iv)), candidate_packet, ok)
            if (.not. ok) error stop 'strict slice candidate packet unavailable'
            call streaming_enlarge_interval(candidate_packet, interval(radius), &
                enlarged_lo, enlarged_hi)
            slice_lo(iv) = enlarged_lo%lo; slice_hi(iv) = enlarged_hi%hi
            call evaluate_streaming(initial, 0.0_dp, v(iv), initial_slice(iv), &
                initial_lo(iv), initial_hi(iv), ok)
            if (.not. ok) error stop 'initial slice sampling failed'
        end do
        do ix = 1, nx + 1
            call density_streaming(eulerian, x(ix), density(ix), density_lo(ix), &
                density_hi(ix), ok)
            if (.not. ok) error stop 'Eulerian density sampling failed'
            particle_density(ix) = 0.0_dp
            do p = 0, modes
                call streaming_mode_accumulate(x(ix), real(p, dp), pcos(p), &
                    psin(p), particle_density(ix), next)
                particle_density(ix) = next
            end do
            call sample_particle_density_band(particle_initial, particles, time, &
                particle_certificate, interval(x(ix)), point_lo, point_hi, ok)
            if (.not. ok) error stop 'particle density sampling failed'
            particle_lo(ix) = point_lo%lo; particle_hi(ix) = point_hi%hi
        end do
        do ix = 1, nx
            call density_streaming_chord(eulerian, x(ix), x(ix + 1), &
                density_padding(ix), ok)
            if (.not. ok) error stop 'Eulerian density chord bound failed'
            call particle_density_chord(particle_initial, particles, time, &
                particle_certificate, interval(x(ix)), interval(x(ix + 1)), &
                chord_padding, ok)
            if (.not. ok) error stop 'particle density chord bound failed'
            particle_padding(ix) = chord_padding%hi
        end do
        do ix = 1, nx + 1
            padding = adjacent_padding(density_padding, ix)
            call streaming_enlarge_interval(interval(density_lo(ix), density_hi(ix)), &
                interval(padding), enlarged_lo, enlarged_hi)
            density_lo(ix) = enlarged_lo%lo; density_hi(ix) = enlarged_hi%hi
            padding = adjacent_padding(particle_padding, ix)
            call streaming_enlarge_interval( &
                interval(particle_lo(ix), particle_hi(ix)), &
                interval(padding), enlarged_lo, enlarged_hi)
            particle_lo(ix) = enlarged_lo%lo; particle_hi(ix) = enlarged_hi%hi
        end do
        call cpu_time(finish)
        postprocess_cpu = finish - start
    end subroutine

    real(dp) function adjacent_padding(padding, node) result(value)
        real(dp), intent(in) :: padding(:)
        integer, intent(in) :: node
        value = 0.0_dp
        if (node > 1) value = padding(node - 1)
        if (node <= size(padding)) value = max(value, padding(node))
    end function

    real(dp) function timer_factor(total, baseline) result(value)
        real(dp), intent(in) :: total, baseline
        value = -1.0_dp
        if (baseline > 0.0_dp) value = total/baseline
    end function

    subroutine write_evidence(number)
        integer, intent(in) :: number
        integer :: j
        write (evidence, '(i0,17(",",es24.16e3))') number, time, radius, &
            candidate_cpu, complete_cpu, timer_factor(complete_cpu, candidate_cpu), &
            particle_push_cpu, particle_project_cpu, particle_complete_cpu, &
            timer_factor(particle_complete_cpu, particle_push_cpu), &
            timer_factor(particle_complete_cpu, &
                particle_push_cpu + particle_project_cpu), &
            particle_complete_cpu - particle_push_cpu - particle_project_cpu, &
            complete_cpu - candidate_cpu, postprocess_cpu, &
            density_lo(1), density_hi(1), particle_lo(1), particle_hi(1)
        frame_time(number) = time; frame_radius(number) = radius
        store_density(:, 1, number) = density_lo; store_density(:, 2, number) = density_hi
        store_particle(:, 1, number) = particle_lo
        store_particle(:, 2, number) = particle_hi
        store_slice(:, 1, number) = slice_lo; store_slice(:, 2, number) = slice_hi
        do j = 1, nx + 1
            write (bands, '(i0,6(",",es24.16e3))') number, time, x(j), &
                density_lo(j), density_hi(j), particle_lo(j), particle_hi(j)
        end do
        do j = 1, nv + 1
            write (slices, '(i0,6(",",es24.16e3))') number, time, v(j), &
                slice_lo(j), slice_hi(j), initial_lo(j), initial_hi(j)
        end do
    end subroutine

    subroutine write_result()
        character(len=*), parameter :: bound = 'certified-continuous finite '// &
            'Fourier/P1 enclosure (lower, upper)'
        call result%create(run)
        call result%attribute('bound_index', '1=lower, 2=upper')
        call result%write('time', frame_time, '1', 'frame time', dims='frame')
        call result%write('x', x, '1', 'x (normalized; period 2 pi)', dims='x')
        call result%write('v', v, '1', 'velocity nodes', dims='v')
        call result%write('velocity_interpolation_radius', frame_radius, '1', &
            'Eulerian velocity interpolation radius', dims='frame', &
            status='certified-continuous')
        call result%write('eulerian_density_bounds', store_density, '1', &
            'normalized Eulerian density enclosure', dims='x,bound,frame', status=bound)
        call result%write('particle_density_bounds', store_particle, '1', &
            'normalized particle Fourier-PIF density enclosure', dims='x,bound,frame', &
            status='certified deterministic particle Fourier quadrature bound; '// &
            'pointwise Dirac bound unavailable')
        call result%write('slice_bounds', store_slice, '1', &
            'distribution slice enclosure at x=0', dims='v,bound,frame', status=bound)
        call result%write('initial_slice_bounds', reshape([initial_lo, initial_hi], &
            [nv + 1, 2]), '1', 'initial distribution slice enclosure', dims='v,bound', &
            status=bound)
        call result%close()
    end subroutine

    subroutine select_cohorts()
        integer, parameter :: phase_ids(3) = [0, 11, 21]
        integer, parameter :: velocity_ids(3) = [7, 23, 39]
        real(dp) :: target_x, target_v, unused_measure
        integer :: phase, row, p
        cohorts = 0
        do phase = 1, 3
            do row = 1, cohort_size
                call streaming_particle_midpoint(real(phase_ids(phase), dp), &
                    real(particle_nx, dp), real(velocity_ids(row), dp), &
                    real(particle_nv, dp), target_x, target_v, unused_measure)
                do p = 1, np
                    if (abs(x0(p) - target_x) > 1.0e-12_dp) cycle
                    if (abs(particle_v(p) - target_v) > 1.0e-12_dp) cycle
                    cohorts(row, phase) = p
                    exit
                end do
            end do
        end do
        if (any(cohorts == 0)) error stop 'fixed particle cohort not found'
    end subroutine

    subroutine write_particle_data()
        integer :: p
        open (newunit=particle_data, file='output/streaming-particles.csv', &
            status='replace')
        write (particle_data, '(a)') 'id,x0,v,weight,x0_lower,x0_upper,'// &
            'v_lower,v_upper,weight_lower,weight_upper,glyph_area'
        do p = 1, np
            write (particle_data, '(i0,10(",",es24.16e3))') p, x0(p), particle_v(p), &
                weight(p), x0_bounds(p)%lo, x0_bounds(p)%hi, &
                velocity_bounds(p)%lo, velocity_bounds(p)%hi, weight_bounds(p)%lo, &
                weight_bounds(p)%hi, glyph(p)
        end do
        close (particle_data)
    end subroutine

    subroutine update_frame(number)
        integer, intent(in) :: number
        call compute_frame(number, timed=.false.)
        call canvas%clear()
        call draw_frame()
        call canvas%set_rendered(.false.)
    end subroutine

    subroutine draw_frame()
        character(len=96) :: heading
        call subplot(2, 2, 1)
        write (heading, '(a,f5.2,a)') 'Characteristic particles | t = ', time, &
            ' | area ~ weight'
        call title(trim(heading))
        call xlabel('x (normalized; period 2 pi)')
        call ylabel('v (normalized)')
        call xlim(0.0_dp, period); call ylim(-1.0_dp, 1.0_dp)
        call scatter(particle_x, particle_v, s=glyph, &
            color=[0.5_dp, 0.5_dp, 0.5_dp], alpha=0.65_dp)
        call cohort(1, 'o', [0.0_dp, 0.4471_dp, 0.6980_dp], 'cohort A')
        call cohort(2, '^', [0.8353_dp, 0.3686_dp, 0.0_dp], 'cohort B')
        call cohort(3, 's', [0.8_dp, 0.475_dp, 0.655_dp], 'cohort C')
        call legend(loc='upper right')
        call subplot(2, 2, 2)
        write (heading, '(a,es10.1e3)') &
            'Eulerian | f color [0.5,1.5] | v-radius = ', radius
        call title(trim(heading))
        call xlabel('x (normalized; period 2 pi)')
        call ylabel('v (normalized)')
        call xlim(0.0_dp, period); call ylim(-1.0_dp, 1.0_dp)
        call pcolormesh(x, v, image, cmap='viridis', vmin=0.5_dp, vmax=1.5_dp)
        call subplot(2, 2, 3)
        call title('Density | certified Eulerian and particle Fourier projection')
        call xlabel('x (normalized; period 2 pi)')
        call ylabel('n = (1/2) integral f dv')
        call xlim(0.0_dp, period); call ylim(0.25_dp, 1.55_dp)
        call envelope(x, particle_lo, particle_hi, 'moccasin', &
            [0.8353_dp, 0.3686_dp, 0.0_dp])
        call envelope(x, density_lo, density_hi, 'lightblue', &
            [0.0_dp, 0.4471_dp, 0.6980_dp])
        call plot(x, particle_density, linestyle='--', &
            label='particle PIF + weak bound', &
            color=[0.8353_dp, 0.3686_dp, 0.0_dp])
        call plot(x, density, label='Eulerian + continuous bound', &
            color=[0.0_dp, 0.4471_dp, 0.6980_dp])
        call legend(loc='lower right')
        call subplot(2, 2, 4)
        call title('Slice x = 0 | Eulerian P1 and full continuous bound')
        call xlabel('v (normalized)')
        call ylabel('f(x=0,v) (normalized)')
        call xlim(-1.0_dp, 1.0_dp); call ylim(0.45_dp, 1.80_dp)
        call envelope(v, slice_lo, slice_hi, 'lightblue', &
            [0.0_dp, 0.4471_dp, 0.6980_dp])
        call envelope(v, initial_lo, initial_hi, 'gainsboro', &
            [0.4_dp, 0.4_dp, 0.4_dp])
        call plot(v, initial_slice, linestyle='--', label='initial profile + bound', &
            color=[0.4_dp, 0.4_dp, 0.4_dp])
        call plot(v, slice, label='Eulerian + continuous bound', &
            color=[0.0_dp, 0.4471_dp, 0.6980_dp])
        call legend(loc='best')
        call tight_layout()
    end subroutine

    subroutine cohort(number, marker, color, label)
        integer, intent(in) :: number
        character(len=*), intent(in) :: marker, label
        real(dp), intent(in) :: color(3)
        integer :: j, p
        do j = 1, cohort_size
            p = cohorts(j, number)
            cohort_x(j) = particle_x(p); cohort_v(j) = particle_v(p)
        end do
        call scatter(cohort_x, cohort_v, s=30.0_dp, marker=marker, color=color, &
            label=label)
    end subroutine

    subroutine envelope(coordinate, lower, upper, fill, color)
        real(dp), intent(in) :: coordinate(:), lower(:), upper(:), color(3)
        character(len=*), intent(in) :: fill
        call fill_between(coordinate, lower, upper, color=fill, alpha=0.5_dp)
        call plot(coordinate, lower, linestyle=':', color=color)
        call plot(coordinate, upper, linestyle=':', color=color)
    end subroutine
end program streaming_comparison_demo
