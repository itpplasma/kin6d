program slab_oscillation_demo
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_interval, only: interval_t, interval, hull, mag
    use kin6d_slab_oscillation, only: slab_oscillation_state_t, &
                                 initialize_slab_oscillation, sample_slab_oscillation, &
                  sample_slab_oscillation_candidate_interval, slab_oscillation_radius, &
                                      slab_oscillation_chord
    use kin6d_slab_oscillation_leaves
    use fortplot, only: figure, get_global_figure, figure_t, subplot, scatter, &
                        plot, fill_between, title, xlabel, ylabel, xlim, ylim, legend, &
                        tight_layout, savefig, suptitle
    use kin6d_demo_style, only: vermillion, black, grey, hold_frames, num
    use fortplot_animation, only: animation_t, animation_func => FuncAnimation, &
                                  save_animation
    use fortplot_system_runtime, only: create_directory_runtime
    use kin6d_result_io, only: run_t, result_file_t
    use kin6d_animation_export, only: animation_export_t
    implicit none
    integer, parameter :: markers = 17, fps = 24
    ! 193 frames span two periods in 8 s at 24 fps; a one-second hold follows.
    integer :: steps = 128, frames = 193
    !> Duration in plasma periods and requested observable radius (1% of q0).
    real(dp) :: periods = 2.0_dp, requested_radius = 0.01_dp
    namelist /slab_oscillation/ steps, frames, periods, requested_radius
    type(run_t) :: run
    type(result_file_t) :: result
    integer :: unit, ios
    character(len=256) :: msg
    real(dp), parameter :: blue(3) = [0.0_dp, 0.4471_dp, 0.6980_dp]
    real(dp), parameter :: ratio = 0.1_dp
    type(slab_oscillation_state_t) :: state
    type(figure_t), pointer :: canvas
    type(animation_t) :: animation
    type(interval_t), allocatable :: qpacket(:), reference_packet(:)
    type(interval_t) :: unused, lo, hi
    type(interval_t) :: padding, margin, error_packet, shift_i, left_i, right_i
    real(dp), allocatable :: time(:), q(:), p(:), reference(:), zero(:)
    real(dp), allocatable :: qlo(:), qhi(:), elo(:), ehi(:), error(:)
    real(dp), allocatable :: chord(:), error_chord(:)
    real(dp) :: initial_x(markers), ion_y(markers), electron_y(markers)
    real(dp) :: electron_x(markers), duration, radius, ql, qu, pl, pu
    real(dp) :: shift, left, right, exact_left, exact_right, exact_shift, ignored
    real(dp) :: field, force, arrow_field, arrow_force
    integer :: i, frame, status, snapshot_frame, evidence
    logical :: ok, snapshot_only
    character(len=256) :: argument, export_dir

    call run%start('slab_oscillation_demo')
    ios = 0; msg = ''
    call run%open_config(unit)
    if (run%has_config()) read (unit, nml=slab_oscillation, iostat=ios, iomsg=msg)
    call run%close_config(unit, 'slab_oscillation', ios, msg)
    if (steps < 1 .or. frames < 2) error stop 'slab_oscillation needs steps >= 1, frames >= 2'
    if (.not. (periods > 0.0_dp .and. requested_radius > 0.0_dp)) &
        error stop 'slab_oscillation periods and requested_radius must be positive'
    call run%open_resolved(unit)
    write (unit, nml=slab_oscillation)
    call run%close_resolved(unit)
    allocate (qpacket(frames), reference_packet(frames), time(frames), q(frames), &
        p(frames), reference(frames), zero(frames), qlo(frames), qhi(frames), &
        elo(frames), ehi(frames), error(frames), chord(frames - 1), &
        error_chord(frames - 1))
    snapshot_only = .false.; snapshot_frame = (frames + 1)/2; export_dir = ''
    i = 1
    do while (i <= run%argument_count())
        call run%argument(i, argument)
        select case (trim(argument))
        case ('--snapshot')
            snapshot_only = .true.
            if (i < run%argument_count()) then
                call run%argument(i + 1, argument)
                read (argument, *, iostat=status) snapshot_frame
                if (status == 0) i = i + 1
            end if
        case ('--export')
            if (i == run%argument_count()) error stop '--export needs a directory'
            call run%argument(i + 1, export_dir); i = i + 1
        case default
            error stop 'usage: demo [--config FILE] [--snapshot [frame]] [--export DIR]'
        end select
        i = i + 1
    end do
    if (snapshot_frame < 1 .or. snapshot_frame > frames) &
        error stop 'frame outside 1..frames'
    call create_directory_runtime('output', ok)
    if (.not. ok) error stop 'output directory unavailable'
    call slab_oscillation_period(periods, duration)
    call initialize_slab_oscillation(state, steps, duration, 1.0_dp, 0.0_dp, ok)
    if (.not. ok) error stop 'numerical slab solve or certificate failed'
    call slab_oscillation_radius(state, radius, ok)
    if (.not. ok .or. radius > requested_radius) error stop 'one-percent certificate failed'
    zero = 0.0_dp
    ion_y = 0.3_dp; electron_y = -0.3_dp
    do i = 1, markers
     call slab_oscillation_marker_grid(real(i - 1, dp), real(markers, dp), initial_x(i))
    end do
    do i = 1, frames
        call slab_oscillation_time(real(i - 1, dp), real(frames - 1, dp), &
                                   duration, time(i))
        call sample_slab_oscillation(state, time(i), q(i), p(i), ql, qu, pl, pu, ok)
        if (.not. ok) error stop 'time sample unavailable'
        call sample_slab_oscillation_candidate_interval(state, interval(time(i)), &
                                                        qpacket(i), unused, ok)
        if (.not. ok) error stop 'candidate packet unavailable'
        call slab_oscillation_reference(time(i), 1.0_dp, 0.0_dp, reference(i), ignored)
        call slab_oscillation_reference_interval(interval(time(i)), interval(1.0_dp), &
                                          interval(0.0_dp), reference_packet(i), unused)
        call slab_oscillation_difference(q(i), reference(i), error(i))
    end do
    do i = 1, frames - 1
        call slab_oscillation_chord(state, time(i), time(i + 1), chord(i), ok)
        if (.not. ok) error stop 'continuous trace chord unavailable'
        call slab_oscillation_chord(state, time(i), time(i + 1), error_chord(i), &
                                    ok, error_curve=.true.)
        if (.not. ok) error stop 'error trace chord unavailable'
    end do
    do i = 1, frames
        padding = interval(0.0_dp)
        if (i > 1) padding = hull(padding, interval(chord(i - 1)))
        if (i < frames) padding = hull(padding, interval(chord(i)))
        call slab_oscillation_margin_interval(interval(radius), &
                                              interval(padding%hi), margin)
        call slab_oscillation_enlarge_interval(qpacket(i), margin, lo, hi)
        qlo(i) = lo%lo; qhi(i) = hi%hi
        call slab_oscillation_difference_interval(interval(qlo(i), qhi(i)), &
                                                    interval(q(i)), error_packet)
        if (mag(error_packet) > requested_radius) error stop 'one-percent full band failed'
        padding = interval(0.0_dp)
        if (i > 1) padding = hull(padding, interval(error_chord(i - 1)))
        if (i < frames) padding = hull(padding, interval(error_chord(i)))
        call slab_oscillation_margin_interval(interval(radius), &
                                              interval(padding%hi), margin)
        call slab_oscillation_difference_interval(qpacket(i), reference_packet(i), &
                                                  error_packet)
        call slab_oscillation_enlarge_interval(error_packet, margin, lo, hi)
        elo(i) = lo%lo; ehi(i) = hi%hi
    end do
    open (newunit=evidence, file='output/plasma-oscillation.csv', status='replace')
    write (evidence, '(a)') 'frame,tau,q_candidate,p_candidate,q_analytical,'// &
        'q_lower,q_upper,error_lower,error_upper,global_residual_radius'
    do i = 1, frames
        write (evidence, '(i0,9(a,es24.16))') i, ',', time(i), ',', q(i), ',', p(i), &
      ',', reference(i), ',', qlo(i), ',', qhi(i), ',', elo(i), ',', ehi(i), ',', radius
    end do
    close (evidence)
    call result%create(run)
    call result%write('tau', time, '1/omega_pe', 'tau = omega_pe t', dims='frame')
    call result%write('q_candidate', q, 'xi0', 'candidate displacement q = xi/xi0', &
        dims='frame', status='candidate')
    call result%write('p_candidate', p, 'xi0 omega_pe', 'candidate velocity', &
        dims='frame', status='candidate')
    call result%write('q_analytical', reference, 'xi0', 'analytical cos oracle', &
        dims='frame', status='oracle')
    call result%write('q_lower', qlo, 'xi0', 'q trace band lower', dims='frame', &
        status='certified-continuous residual radius plus candidate chord')
    call result%write('q_upper', qhi, 'xi0', 'q trace band upper', dims='frame', &
        status='certified-continuous residual radius plus candidate chord')
    call result%write('error_lower', elo, 'xi0', 'q - analytical lower', dims='frame', &
        status='certified-continuous residual radius plus candidate chord')
    call result%write('error_upper', ehi, 'xi0', 'q - analytical upper', dims='frame', &
        status='certified-continuous residual radius plus candidate chord')
    call result%write('global_residual_radius', radius, 'xi0', &
        'global Hermite/Bernstein residual radius for q and p', &
        status='certified-continuous')
    call result%close()
    print '(a)', 'model=cold_planar_electron_slab; ions=fixed; tau=omega_pe*t'
    print '(a)', 'q=xi/xi0; equation=qprimeprime=-q; solver=FortNum_Cash_Karp_RK5(4)'
    print '(a)', 'checker=Hermite_residual_Bernstein_Duhamel; '// &
        'analytical_oracle_not_required'
    print '(a)', 'markers=collective_slab_teaching_glyphs; not_3D_periodic_k0_OCP'
    print '(a)', 'Eulerian=translated_box_density; edge_uncertainty=[0,n0]'
    print '(a)', 'position_x_over_a; display_xi0_over_a=0.1; normalized_outputs_unit1'
    print '(a)', 'time_trace_band=global_residual_radius_plus_candidate_chord'
    print '(a)', 'continuous_VJP_bound=unavailable; '// &
        'fixed_mesh_tangent_IVP_bound=available'
    print '(a,es24.16)', 'global_q_and_p_residual_radius=', radius
    if (len_trim(export_dir) > 0) call export_animation()
    call figure(figsize=[13.2_dp, 8.8_dp])
    canvas => get_global_figure()
    frame = snapshot_frame
    call draw_frame()
    call savefig('output/plasma-oscillation.png')
    if (.not. snapshot_only) then
        animation = animation_func(update_frame, frames=frames + hold_frames, &
            interval=42, fig=canvas)
        call save_animation(animation, 'output/plasma-oscillation.mp4', &
                              fps=fps, status=status)
        if (status /= 0) error stop 'animation encoding failed'
    end if

contains

    subroutine update_frame(number)
        integer, intent(in) :: number

        frame = min(number, frames)
        call canvas%clear()
        call draw_frame()
        call canvas%set_rendered(.false.)
    end subroutine update_frame

    subroutine draw_frame()

        call slab_oscillation_display_displacement(q(frame), ratio, shift)
        call slab_oscillation_display_displacement(reference(frame), ratio, exact_shift)
        call slab_oscillation_edges(shift, 1.0_dp, left, right)
        call slab_oscillation_edges(exact_shift, 1.0_dp, exact_left, exact_right)
        call slab_oscillation_enlarge_interval(qpacket(frame), interval(radius), lo, hi)
        call slab_oscillation_display_displacement_interval(hull(lo, hi), &
                                                            interval(ratio), shift_i)
        call slab_oscillation_edges_interval(shift_i, interval(1.0_dp), left_i, right_i)
        do i = 1, markers
            call slab_oscillation_marker(initial_x(i), shift, electron_x(i))
        end do
        call suptitle('A displaced electron slab feels a linear restoring field and '// &
            'oscillates at the plasma frequency')
        call subplot(2, 2, 1)
        call title('Electrons (o) over fixed ions (+)  |  $\tau$ = '//num(time(frame), 2))
        call xlabel('x / a   (display amplification: $\xi_0/a$ = 0.1)')
        call ylabel('species lanes (display only)')
        call xlim(-1.3_dp, 1.3_dp); call ylim(-0.7_dp, 1.1_dp)
        call scatter(initial_x, ion_y, marker='+', s=45.0_dp, color=blue)
        call scatter(electron_x, electron_y, marker='o', s=35.0_dp, color=vermillion)
        call slab_oscillation_field_force(q(frame), field, force)
        call slab_oscillation_display_displacement(field, 0.3_dp, arrow_field)
        call slab_oscillation_display_displacement(force, 0.3_dp, arrow_force)
        call plot([0.0_dp, arrow_field], [0.65_dp, 0.65_dp], color=blue, linewidth=2.0_dp, &
                  label='field E / E0 (from x = 0)')
        call plot([0.0_dp, arrow_force], [0.9_dp, 0.9_dp], color=vermillion, &
                  linestyle='--', linewidth=2.0_dp, label='force on electrons F / F0')
        call legend(loc='lower center')
        call subplot(2, 2, 2)
        call title('Eulerian slab densities (bands: uncertain edges)')
        call xlabel('x / a'); call ylabel('$n/n_0$')
        call xlim(-1.3_dp, 1.3_dp); call ylim(-0.05_dp, 1.15_dp)
        call edge_band(left_i); call edge_band(right_i)
        call plot([-1.3_dp, -1.0_dp, -1.0_dp, 1.0_dp, 1.0_dp, 1.3_dp], &
                  [0.0_dp, 0.0_dp, 1.0_dp, 1.0_dp, 0.0_dp, 0.0_dp], &
                  color=grey, linestyle=':', linewidth=1.5_dp, label='fixed ions')
        call plot([-1.3_dp, left, left, right, right, 1.3_dp], &
                  [0.0_dp, 0.0_dp, 1.0_dp, 1.0_dp, 0.0_dp, 0.0_dp], &
                  color=vermillion, linewidth=2.0_dp, label='numerical electrons')
        call plot([-1.3_dp, exact_left, exact_left, exact_right, exact_right, 1.3_dp], &
                  [0.0_dp, 0.0_dp, 1.0_dp, 1.0_dp, 0.0_dp, 0.0_dp], &
                  color=black, linestyle='--', linewidth=1.2_dp, label='analytical electrons')
        call legend(loc='lower center')
        call subplot(2, 2, 3)
        call title('Displacement: numerical trajectory and analytical cos τ')
        call xlabel('$\tau = \omega_{pe} t$'); call ylabel('$\xi/\xi_0$')
        call xlim(0.0_dp, duration); call ylim(-1.25_dp, 1.5_dp)
        call fill_between(time, qlo, qhi, color=blue, alpha=0.25_dp)
        call plot(time, q, color=blue, linewidth=2.0_dp, label='numerical RK5 ± bound')
        call plot(time, reference, color=black, linestyle='--', linewidth=1.2_dp, &
                  label='analytical cos τ')
        call scatter([time(frame)], [q(frame)], color=blue, marker='o', s=40.0_dp)
        call legend(loc='upper center')
        call subplot(2, 2, 4)
        call title('Numerical minus analytical, with the certified bound')
        call xlabel('$\tau = \omega_{pe} t$'); call ylabel('$\Delta\xi/\xi_0$')
        call xlim(0.0_dp, duration); call ylim(-0.003_dp, 0.003_dp)
        call fill_between(time, elo, ehi, color=blue, alpha=0.25_dp)
        call plot(time, zero, color=black, linestyle='--', linewidth=1.2_dp)
        call plot(time, error, color=blue, linewidth=2.0_dp, &
                  label='numerical − analytical ± bound (incl. chord margin)')
        call scatter([time(frame)], [error(frame)], color=blue, marker='o', s=35.0_dp)
        call legend(loc='upper center')
        call tight_layout()
    end subroutine draw_frame

    !> Animation data for the lecture scene (slug plasma-oscillation): the
    !> per-frame candidate trajectory, its connected certified band, the
    !> analytical oracle, slab edges at the declared display ratio xi0/a and
    !> the normalized field/force, all from the quantities drawn above.
    subroutine export_animation()
        type(animation_export_t) :: export
        real(dp) :: edges(6, frames), loads(2, frames)
        do i = 1, frames
            call slab_oscillation_display_displacement(q(i), ratio, shift)
            call slab_oscillation_edges(shift, 1.0_dp, left, right)
            call slab_oscillation_enlarge_interval(qpacket(i), interval(radius), lo, hi)
            call slab_oscillation_display_displacement_interval(hull(lo, hi), &
                                                                interval(ratio), shift_i)
            call slab_oscillation_edges_interval(shift_i, interval(1.0_dp), left_i, &
                                                 right_i)
            edges(:, i) = [left, right, left_i%lo, left_i%hi, right_i%lo, right_i%hi]
            call slab_oscillation_field_force(q(i), field, force)
            loads(:, i) = [field, force]
        end do
        call export%create(run, export_dir, 'plasma-oscillation')
        call export%note('model', 'cold planar electron slab, fixed ions, '// &
            'tau = omega_pe t, q = xi/xi0, q'''' = -q, q(0) = 1, p(0) = 0')
        call export%note('solver', 'FortNum Cash-Karp RK5(4), fixed mesh (candidate)')
        call export%note('band', 'q_lower/q_upper: global Hermite/Bernstein '// &
            'residual radius (certified-continuous) plus candidate chord margin')
        call export%note('oracle', 'q_analytical = cos tau, independent overlay')
        call export%parameter('steps', real(steps, dp), '1', 'RK steps')
        call export%parameter('periods', periods, '1', 'duration in plasma periods')
        call export%parameter('display_ratio', ratio, '1', 'xi0/a, display only')
        call export%parameter('residual_radius', radius, 'xi0', &
            'global certified residual radius of q and p')
        call export%real32('tau', time, [frames], '1/omega_pe', 'omega_pe t per frame')
        call export%real32('q', q, [frames], 'xi0', 'candidate displacement xi/xi0')
        call export%real32('p', p, [frames], 'xi0 omega_pe', 'candidate velocity')
        call export%real32('q_analytical', reference, [frames], 'xi0', 'cos tau')
        call export%real32('q_lower', qlo, [frames], 'xi0', 'connected band lower')
        call export%real32('q_upper', qhi, [frames], 'xi0', 'connected band upper')
        call export%real32('edges', reshape(edges, [6*frames]), [6, frames], 'a', &
            'x/a of left, right edge and enclosures left lo/hi, right lo/hi')
        call export%real32('field_force', reshape(loads, [2*frames]), [2, frames], &
            'E0, F0', 'E/E0 and F/F0 inside the slab (candidate)')
        call export%close()
        print '(2a)', 'animation data written to ', trim(export_dir)//'/plasma-oscillation'
    end subroutine export_animation

    subroutine edge_band(edge)
        type(interval_t), intent(in) :: edge

        call fill_between([edge%lo, edge%hi], [0.0_dp, 0.0_dp], [1.0_dp, 1.0_dp], &
                          color=vermillion, alpha=0.25_dp)
    end subroutine edge_band

end program slab_oscillation_demo
