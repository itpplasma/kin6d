program interaction_demo
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_interval, only: interval_t, ipoint
    use kin6d_interaction_laboratory
    use kin6d_interaction_leaves, only: interaction_energy, &
                                        interaction_particle_error_interval, &
                                        interaction_frame_time, &
                                        interaction_frame_time_interval, &
                                        interaction_particle_coordinate, &
                                        interaction_particle_coordinate_interval
    use kin6d_coulomb_leaves, only: coulomb_pair
    use kin6d_species, only: units_t, atomic_units, electron, pair_coupling
    use fortplot, only: figure, figure_t, get_global_figure, subplot, plot, &
                        scatter, quiver, &
                        title, xlabel, ylabel, xlim, ylim, legend, &
                        tight_layout, savefig
    use fortplot_animation, only: animation_t, &
                                  animation_func => FuncAnimation, save_animation
    use fortplot_system_runtime, only: create_directory_runtime
    use kin6d_result_io, only: run_t, result_file_t
    implicit none
    integer, parameter :: fps = 24
    integer :: frames = 121, substeps = 100
    real(dp) :: h = .0005_dp, duration
    !> Neutral elastic contact radius and requested particle-coordinate radius.
    real(dp) :: contact_radius = 1._dp, radius_target = .01_dp
    real(dp), parameter :: two_pi = 6.2831853071795864769_dp
    type(units_t) :: units
    real(dp) :: coupling
    real(dp), parameter :: blue(3) = [.12_dp, .38_dp, .7_dp], &
                           orange(3) = [.85_dp, .28_dp, .12_dp]
    type(interaction_state_t) :: state, fast
    type(interval_t) :: rb(3), ub(3), error_radius(3), particle_box(3)
    type(interval_t), allocatable :: time_bounds(:)
    type(figure_t), pointer :: canvas
    type(animation_t) :: animation
    real(dp) :: r0(3) = [-6._dp, .6_dp, 0._dp], u0(3) = [2._dp, 0._dp, 0._dp]
    namelist /interaction/ frames, substeps, h, r0, u0, contact_radius, radius_target
    type(run_t) :: run
    type(result_file_t) :: result
    integer :: unit, ios
    character(len=256) :: msg
    real(dp), allocatable :: charged(:, :), neutral(:, :), times(:), energy(:)
    real(dp), allocatable :: radii(:), boxlo(:, :, :), boxhi(:, :, :)
    real(dp) :: r(3), u(3), en0
    real(dp), allocatable :: radius_plot(:), energy_plot(:)
    real(dp), allocatable :: charged_position(:, :, :), neutral_position(:, :, :)
    real(dp) :: particle_sign
    real(dp) :: start, finish, solve_seconds, check_seconds, factor, f(3), potential
    integer :: frame, k, status, out, snapshot_frame, last_bounded, particle
    logical :: ok, bounded, snapshot_only
    logical, allocatable :: certified(:), budget_met(:)
    character(64) :: argument
    character(150) :: heading
    call run%start('interaction_demo')
    ios = 0; msg = ''
    call run%open_config(unit)
    if (run%has_config()) read (unit, nml=interaction, iostat=ios, iomsg=msg)
    call run%close_config(unit, 'interaction', ios, msg)
    if (frames < 2 .or. substeps < 1) error stop 'interaction needs frames >= 2, substeps >= 1'
    if (.not. (h > 0.0_dp .and. contact_radius > 0.0_dp .and. radius_target > 0.0_dp)) &
        error stop 'interaction h, contact_radius and radius_target must be positive'
    call run%open_resolved(unit)
    write (unit, nml=interaction)
    call run%close_resolved(unit)
    duration = real((frames - 1)*substeps, dp)*h
    allocate (charged(3, frames), neutral(3, frames), times(frames), energy(frames), &
        radii(frames), boxlo(3, 2, frames), boxhi(3, 2, frames), radius_plot(frames), &
        energy_plot(frames), charged_position(3, 2, frames), &
        neutral_position(3, 2, frames), time_bounds(frames), certified(frames), &
        budget_met(frames))
    snapshot_only = .false.; snapshot_frame = (frames + 1)/2
    if (run%argument_count() > 2) error stop 'usage: interaction_demo ['// &
        '--config FILE] [--snapshot [frame]]'
    if (run%argument_count() >= 1) then
        call run%argument(1, argument)
        if (trim(argument) /= '--snapshot') error stop 'expected --snapshot'
        snapshot_only = .true.
    end if
    if (run%argument_count() == 2) then
        call run%argument(2, argument)
        read (argument, *, iostat=status) snapshot_frame
        if (status /= 0) error stop 'invalid snapshot frame'
    end if
  if (snapshot_frame < 1 .or. snapshot_frame > frames) error stop 'frame outside 1..frames'
    call create_directory_runtime('output', ok)
    if (.not. ok) error stop 'cannot create output'
    ! Two electrons in Hartree atomic units: g = 1, mu = 1/2, coupling g/mu = 2.
    units = atomic_units()
    call init_interaction(fast, r0, u0, units, electron(), electron(), ok)
    if (.not. ok) error stop 'interaction setup failed'
    coupling = interaction_coupling(fast)
    call cpu_time(start)
    do k = 1, (frames - 1)*substeps
        call step_interaction(fast, h, ok, certify=.false.)
        if (.not. ok) error stop 'candidate integration failed'
    end do
    call cpu_time(finish); solve_seconds = finish - start
    call init_interaction(state, r0, u0, units, electron(), electron(), ok)
    call cpu_time(start)
    do k = 1, (frames - 1)*substeps
        call step_interaction(state, h, ok)
        if (.not. ok) error stop 'checked candidate integration failed'
    end do
    call cpu_time(finish); check_seconds = finish - start
    call init_interaction(state, r0, u0, units, electron(), electron(), ok)
    call interaction_energy(r0(1), r0(2), r0(3), u0(1), u0(2), u0(3), coupling, en0)
    do frame = 1, frames
        if (frame > 1) then
            do k = 1, substeps
                call step_interaction(state, h, ok)
                if (.not. ok) error stop 'candidate integration failed'
            end do
        end if
        call interaction_frame_time(real(frame - 1, dp), real(substeps, dp), h, &
                                    times(frame))
        call interaction_frame_time_interval(ipoint(real(frame - 1, dp)), &
                                             ipoint(real(substeps, dp)), ipoint(h), &
                                             time_bounds(frame))
        call sample_interaction(state, r, u, rb, ub, ok, bounded)
        charged(:, frame) = r; certified(frame) = bounded
        do particle = 1, 2
            particle_sign = real(3 - 2*particle, dp)
            call interaction_particle_coordinate(r, particle_sign, &
                                                 charged_position(:, particle, frame))
            call interaction_particle_coordinate_interval(rb, ipoint(particle_sign), &
                                                          particle_box)
            boxlo(:, particle, frame) = particle_box%lo
            boxhi(:, particle, frame) = particle_box%hi
        end do
        radii(frame) = -1
        if (bounded) then
            call interaction_particle_error_interval(rb, ipoint(r), error_radius)
            radii(frame) = maxval(error_radius%hi)
        end if
        budget_met(frame) = bounded
        if (bounded) budget_met(frame) = radii(frame) <= radius_target
        call interaction_energy(r(1), r(2), r(3), u(1), u(2), u(3), coupling, &
                                energy(frame))
        call neutral_contact_trajectory(r0, u0, contact_radius, times(frame), r, u, ok)
        if (.not. ok) error stop 'neutral contact failed'
        neutral(:, frame) = r
        do particle = 1, 2
            call interaction_particle_coordinate(r, real(3 - 2*particle, dp), &
                                                 neutral_position(:, particle, frame))
        end do
    end do
    factor = -1
    if (solve_seconds > 0) factor = check_seconds/solve_seconds
    last_bounded = 0
    do frame = 1, frames
        if (certified(frame)) last_bounded = frame
    end do
    radius_plot = 1.e5_dp*radii
    energy_plot = (energy - en0)*1.e8_dp/en0
    open (newunit=out, file='output/interaction-evidence.csv', status='replace')
    write (out, &
           '(a)') 'frame,t,t_lower,t_upper,endpoint_enclosure_available,one_pe'// &
        'rcent_budget_met,'// &
        'particle_coordinate_radius,relative_specific_energy_diagnostic'
    do frame = 1, frames
        write (out, '(i0,3(a,es25.16e3),2(a,l1),2(a,es25.16e3))') frame, ',', &
            times(frame), &
            ',', time_bounds(frame)%lo, ',', time_bounds(frame)%hi, ',', &
            certified(frame), &
            ',', budget_met(frame), ',', radii(frame), ',', energy(frame)
    end do
    close (out)
    call result%create(run)
    call result%attribute('model', 'two electrons, Hartree atomic units, '// &
        'unscreened Coulomb vs neutral elastic contact')
    call result%write('t', times, 'hbar/E_h', 'frame time', dims='frame')
    call result%write('t_lower', time_bounds%lo, 'hbar/E_h', 'frame time lower', &
        dims='frame', status='certified rounding enclosure')
    call result%write('t_upper', time_bounds%hi, 'hbar/E_h', 'frame time upper', &
        dims='frame', status='certified rounding enclosure')
    call result%write('relative_position', charged, 'a_0', &
        'Coulomb relative position r1 - r2', dims='xyz,frame', status='candidate')
    call result%write('particle_lower', boxlo, 'a_0', &
        'Coulomb particle coordinate enclosure lower', dims='xyz,particle,frame', &
        status='certified-continuous where endpoint_enclosure_available')
    call result%write('particle_upper', boxhi, 'a_0', &
        'Coulomb particle coordinate enclosure upper', dims='xyz,particle,frame', &
        status='certified-continuous where endpoint_enclosure_available')
    call result%write('endpoint_enclosure_available', merge(1.0_dp, 0.0_dp, certified), &
        '1', '1 if the interval Picard tube reached this frame', dims='frame')
    call result%write('particle_coordinate_radius', radii, 'a_0', &
        'max particle coordinate radius (-1 = unavailable)', dims='frame', &
        status='certified-continuous where endpoint_enclosure_available')
    call result%write('relative_specific_energy', energy, 'atomic units (per mass)', &
        'specific relative energy', dims='frame', status='diagnostic')
    call result%write('neutral_relative_position', neutral, 'a_0', &
        'neutral elastic contact relative position', dims='xyz,frame', &
        status='analytical reference; event sensitivity unavailable')
    call result%close()
    print '(a)', &
        'model=two_equal_mass_particles; Cartesian_3D; C=1; reduced_mass=0.5'
    print '(a)', &
        'force=unscreened_repulsive_Coulomb; no_softening; no_prescribed_screening'
    print '(a)', &
        'neutral=ideal_elastic_contact_radius_1; event_sensitivity=unavailable'
    print '(a)', &
        'trajectory=interval_Picard_tube_plus_integral_Taylor; failures_explicit'
    print '(a)', &
        'energy=diagnostic_only; relativistic_and_physical_model_remainders'// &
        '=unavailable'
    print '(a)', &
        'timers=CPU_seconds; rendering_and_sampling_excluded; checker_attem'// &
        'pts_until_failure'
    print '(a)', &
        'radius_target=0.01; physical_coordinate_scale=1; absolute_floor=1e-12'
    print '(a)', &
        'colored_traces=candidate_polylines; interval_boxes=observation_end'// &
        'points_only'
    print '(a)', &
        'energy_observable_and_force_arrow_bounds=unavailable; arrows=candi'// &
        'date_forces'
    print '(a,es12.4)', 'candidate_s=', solve_seconds
    print '(a,es12.4)', 'candidate_plus_checker_s=', check_seconds
    print '(a,es12.4)', 'checker_extra_s=', max(0._dp, &
                                                check_seconds - solve_seconds)
    print '(a,es12.4)', 'complete_factor=', factor
    print '(a,i0,a,i0)', 'bounded_frames=', last_bounded, '/', frames
    call figure(figsize=[13._dp, 9._dp]); canvas => get_global_figure()
    frame = snapshot_frame; call draw_frame()
    call savefig('output/collective-response.png')
    if (.not. snapshot_only) then
        animation = animation_func(update_frame, frames=frames, interval=42, &
                                   fig=canvas)
        call save_animation(animation, 'output/collective-response.mp4', fps=fps, &
                            status=status)
        if (status /= 0) error stop 'animation failed'
    end if
contains
    subroutine update_frame(index)
        integer, intent(in)::index
        frame = index; call draw_frame()
    end subroutine
    subroutine draw_frame()
        integer :: particle
        integer :: outline_index
        real(dp) :: sign, bx(5), by(5), angles(65), outline_x(65), outline_y(65)
        do outline_index = 1, 65
            angles(outline_index) = two_pi*real(outline_index - 1, dp)/64
        end do
        call canvas%clear()
        call subplot(2, 2, 1)
        write (heading, '(a,f5.2,a)') 'Neutral contact | t=', times(frame), &
            ' | event derivative unavailable'
        call title(trim(heading))
        do particle = 1, 2
            sign = real(3 - 2*particle, dp)
            call plot(neutral_position(1, particle, 1:frame), &
                      neutral_position(2, particle, 1:frame), &
                      color=blue, linewidth=1.7_dp)
            call scatter([neutral_position(1, particle, frame)], &
                         [neutral_position(2, particle, frame)], &
                         color=blue, s=95._dp, marker='o')
            ! Contact sphere of radius 1/2 around each neutral particle.
            outline_x = neutral_position(1, particle, frame) + cos(angles)/2
            outline_y = neutral_position(2, particle, frame) + sin(angles)/2
            call plot(outline_x, outline_y, color=blue, linewidth=1.5_dp)
        end do
        call xlabel('x [1] | Cartesian 3D orbit, z=0 by symmetry')
        call ylabel('y [1] | contact distance=1; event bound unavailable')
        call xlim(-3.6_dp, 3.6_dp); call ylim(-4.3_dp, 4.3_dp)
        call subplot(2, 2, 2)
        if (budget_met(frame)) then
            write (heading, '(a,f5.2,a)') 'Coulomb | t=', times(frame), &
                ' | endpoint boxes; 1% budget met'
        else if (certified(frame)) then
            write (heading, '(a,f5.2,a)') 'Coulomb | t=', times(frame), &
                ' | endpoint boxes; 1% budget missed'
        else
            write (heading, '(a,f5.2,a)') 'Unscreened repulsive Coulomb | t=', &
                times(frame), ' | candidate: bound failed'
        end if
        call title(trim(heading))
        call coulomb_pair(charged(1, frame), charged(2, frame), charged(3, frame), &
                          pair_coupling(units, electron(), electron()), &
                          f(1), f(2), f(3), potential)
        do particle = 1, 2
            sign = real(3 - 2*particle, dp)
            call plot(charged_position(1, particle, 1:frame), &
                      charged_position(2, particle, 1:frame), &
                      color=orange, linewidth=1.7_dp)
            call scatter([charged_position(1, particle, frame)], &
                         [charged_position(2, particle, frame)], &
                         color=orange, s=95._dp, marker='o')
            call quiver([charged_position(1, particle, frame)], &
                        [charged_position(2, particle, frame)], &
                        [sign*.35_dp*f(1)], [sign*.35_dp*f(2)], color=orange, &
                        scale=1._dp, angles='xy', scale_units='xy', units='xy', &
                        width=.025_dp)
            if (certified(frame)) then
                bx = [boxlo(1, particle, frame), boxhi(1, particle, frame), boxhi(1, &
                                                                     particle, frame), &
                      boxlo(1, particle, frame), boxlo(1, particle, frame)]
                by = [boxlo(2, particle, frame), boxlo(2, particle, frame), boxhi(2, &
                                                                     particle, frame), &
                      boxhi(2, particle, frame), boxlo(2, particle, frame)]
                call plot(bx, by, color=orange, linewidth=1._dp)
            end if
        end do
        call xlabel('x [1] | candidate arrows: 0.35 F; equal opposite forces')
        call ylabel('y [1] | no screening or confining force')
        call xlim(-3.5_dp, 3.5_dp); call ylim(-3._dp, 3._dp)
        call subplot(2, 2, 3)
        call title('Particle coordinate radius | observation endpoints only')
        if (last_bounded > 0) call scatter(times(1:last_bounded), &
                                     radius_plot(1:last_bounded), color=orange, s=9._dp)
        if (certified(frame)) call scatter([times(frame)], [radius_plot(frame)], &
                                           color=orange, s=35._dp)
        call xlabel('t [1] | endpoint boxes only')
        call ylabel('max coordinate radius / 1e-5 [1]')
        call xlim(0._dp, duration)
        call subplot(2, 2, 4)
        call title('Energy drift / 1e-8 | diagnostic, trajectory bound separate')
        call plot(times, energy_plot, color=orange, linewidth=1.7_dp)
        call scatter([times(frame)], [energy_plot(frame)], &
                     color=orange, s=35._dp)
        call xlabel('t [1] | physical-model uncertainty unavailable')
        call ylabel('relative specific-energy drift / 1e-8 [1]')
        call xlim(0._dp, duration)
        call tight_layout()
    end subroutine
end program
