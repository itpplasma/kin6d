program microscopic_cloud_demo
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_status, only: fortnum_status_t, FORTNUM_OK
    use kin6d_microscopic_cloud, only: cloud_model_t, cloud_system, cloud_initialize
    use kin6d_nbody, only: nbody_system_t, nbody_step, nbody_energy, &
                           nbody_encounter_step
    use fortplot, only: figure, figure_t, get_global_figure, subplot, scatter, &
                  plot, title, xlabel, ylabel, xlim, ylim, legend, tight_layout, savefig, &
                  suptitle
    use kin6d_demo_style, only: vermillion, black, grey, hold_frames, num
    use fortplot_animation, only: animation_t, animation_func => FuncAnimation, &
                                  save_animation
    use fortplot_system_runtime, only: create_directory_runtime
    use kin6d_result_io, only: run_t, result_file_t
    implicit none
    integer, parameter :: fps = 24
    ! Namelist &microscopic_cloud: electrons, output frames, substeps of dt per
    ! frame, encounter control eta, fixed central impurity charge of run 2 and
    ! the neutralizing sphere (defaults from kin6d_microscopic_cloud).
    ! 241 frames of 50 substeps: the same 12000 steps as before, shown in 10 s.
    integer :: particles = 32, frames = 241, substeps = 50
    real(dp) :: dt = 0.001_dp, encounter_eta = 0.05_dp, impurity_charge = -2.0_dp
    real(dp) :: sphere_charge, sphere_radius
    namelist /microscopic_cloud/ particles, frames, substeps, dt, encounter_eta, &
        impurity_charge, sphere_charge, sphere_radius
    type(run_t) :: app
    type(result_file_t) :: result
    integer :: ios
    character(len=256) :: msg
    real(dp), parameter :: two_pi = 6.2831853071795864769_dp
    real(dp), parameter :: blue(3) = [0.0_dp, 0.4471_dp, 0.6980_dp]
    ! Width/height of the cloud axes boxes in pixels (measured), for equal aspect.
    real(dp), parameter :: box_aspect = 1.594_dp
    type(cloud_model_t) :: model(2)
    type(nbody_system_t) :: system
    type(fortnum_status_t) :: status
    type(figure_t), pointer :: canvas
    type(animation_t) :: animation
    real(dp), allocatable :: state(:), trajectory(:, :, :)
    real(dp), allocatable :: time(:), defect(:, :), display_defect(:, :)
    real(dp) :: initial_energy, energy
    real(dp), allocatable :: candidate_cpu(:, :), diagnostic_cpu(:, :)
    real(dp) :: start, finish
    real(dp), allocatable :: xp(:), yp(:), radii(:), fraction(:), cdf_x(:), cdf_y(:)
    real(dp) :: circle_x(101), circle_y(101), phase, scratch, max_view, energy_view
    real(dp) :: limit
    integer :: run, frame, step, p, j, k, unit, snapshot_frame, encode_status
    integer :: substep, pieces, refined(2)
    logical :: ok, snapshot_only
    character(len=64) :: argument
    call app%start('microscopic_cloud_demo')
    sphere_charge = model(1)%sphere_charge; sphere_radius = model(1)%sphere_radius
    ios = 0; msg = ''
    call app%open_config(unit)
    if (app%has_config()) read (unit, nml=microscopic_cloud, iostat=ios, iomsg=msg)
    call app%close_config(unit, 'microscopic_cloud', ios, msg)
    if (particles < 1 .or. frames < 2 .or. substeps < 1) error stop &
        'particles, frames and substeps must be positive (frames >= 2)'
    if (.not. (dt > 0.0_dp .and. encounter_eta > 0.0_dp .and. sphere_radius > 0.0_dp &
        .and. sphere_charge >= 0.0_dp .and. impurity_charge <= 0.0_dp)) error stop &
        'invalid cloud step control or charges (impurity must be <= 0)'
    call app%open_resolved(unit)
    write (unit, nml=microscopic_cloud)
    call app%close_resolved(unit)
    allocate (state(6*particles), trajectory(6*particles, frames, 2), time(frames), &
        defect(frames, 2), display_defect(frames, 2), candidate_cpu(frames, 2), &
        diagnostic_cpu(frames, 2), xp(particles), yp(particles), radii(particles), &
        fraction(particles), cdf_x(2*particles + 2), cdf_y(2*particles + 2))
    snapshot_only = .false.; snapshot_frame = (frames + 1)/2
    if (app%argument_count() > 2) error stop &
        'usage: microscopic_cloud_demo [--config FILE] [--snapshot [frame]]'
    if (app%argument_count() >= 1) then
        call app%argument(1, argument)
        if (trim(argument) /= '--snapshot') error stop 'expected --snapshot'
        snapshot_only = .true.
    end if
    if (app%argument_count() == 2) then
        call app%argument(2, argument)
        read (argument, *, iostat=encode_status) snapshot_frame
        if (encode_status /= 0) error stop 'snapshot frame must be an integer'
    end if
    if (snapshot_frame < 1 .or. snapshot_frame > frames) &
        error stop 'snapshot frame outside 1..frames'
    call create_directory_runtime('output', ok)
    if (.not. ok) error stop 'cannot create output directory'
    ! Fixed negative impurity at the origin (repulsive for electrons).
    model%sphere_charge = sphere_charge; model%sphere_radius = sphere_radius
    model(2)%impurity_charge = impurity_charge
    candidate_cpu = 0.0_dp; diagnostic_cpu = 0.0_dp; refined = 0
    do run = 1, 2
        call cloud_system(model(run), particles, system, ok)
        if (.not. ok) error stop 'cloud system setup failed'
        call cloud_initialize(state)
        call nbody_energy(system, state, initial_energy, status)
        if (status%code /= FORTNUM_OK) error stop 'cloud energy initialization failed'
        trajectory(:, 1, run) = state; defect(1, run) = 0.0_dp
        time(1) = 0.0_dp
        do frame = 2, frames
            call cpu_time(start)
            do step = 1, substeps
                ! Encounter-resolved control: split a substep only when the
                ! Kepler/crossing time of the closest encounter requires it.
                call nbody_encounter_step(system, state, encounter_eta, limit, status)
                if (status%code /= FORTNUM_OK) error stop status%msg
                pieces = max(1, ceiling(dt/limit))
                if (pieces > 1) refined(run) = refined(run) + 1
                do substep = 1, pieces
                    call nbody_step(system, state, dt/pieces, status)
                    if (status%code /= FORTNUM_OK) error stop status%msg
                end do
            end do
            call cpu_time(finish)
            candidate_cpu(frame, run) = finish - start
            trajectory(:, frame, run) = state
            call cpu_time(start)
            call nbody_energy(system, state, energy, status)
            if (status%code /= FORTNUM_OK) error stop 'cloud energy diagnostic failed'
            defect(frame, run) = (energy - initial_energy)/abs(initial_energy)
            call cpu_time(finish)
            diagnostic_cpu(frame, run) = finish - start
            time(frame) = real(frame - 1, dp)*real(substeps, dp)*dt
        end do
    end do
    open (newunit=unit, file='output/microscopic-cloud.csv', status='replace')
    write (unit, '(a)') '# candidate; isolated 3D3V Coulomb; no screened forces'
    write (unit, '(a)') '# trajectory/numerical/model/statistical bounds unavailable'
    write (unit, '(a)') '# q=-1,m=1,Coulomb prefactor=1; sphere Q='// &
        trim(number(sphere_charge))//',R='//trim(number(sphere_radius))
   write (unit, '(a)') '# no thermostat, wall, cutoff, self force or collision operator'
    write (unit, '(a)') 'run,frame,t,solve_cpu_s,diagnostic_extra_cpu_s,'// &
        'diagnostic_multiplier,relative_energy_defect'
    do run = 1, 2
        do frame = 2, frames
            write (unit, '(i0,a,i0,5(a,es24.16e3))') run, ',', frame, ',', &
                time(frame), ',', candidate_cpu(frame, run), ',', &
                diagnostic_cpu(frame, run), ',', &
                1.0_dp + diagnostic_cpu(frame, run)/candidate_cpu(frame, run), ',', &
                defect(frame, run)
        end do
    end do
    close (unit)
    open (newunit=unit, file='output/microscopic-cloud-particles.csv', status='replace')
    write (unit, '(a)') 'run,frame,particle,t,x,y,z,vx,vy,vz,status'
    do run = 1, 2
        do frame = 1, frames
            do p = 1, particles
                j = 3*p - 2
                write (unit, '(i0,2(a,i0),7(a,es24.16e3),a)') run, ',', frame, &
                    ',', p, ',', time(frame), ',', trajectory(j, frame, run), ',', &
                    trajectory(j + 1, frame, run), ',', trajectory(j + 2, frame, run), &
                    ',', trajectory(3*particles + j, frame, run), ',', &
                    trajectory(3*particles + j + 1, frame, run), ',', &
                    trajectory(3*particles + j + 2, frame, run), ',candidate'
            end do
        end do
    end do
    close (unit)
    do p = 1, 101
        phase = two_pi*real(p - 1, dp)/100
        circle_x(p) = model(1)%sphere_radius*cos(phase)
        circle_y(p) = model(1)%sphere_radius*sin(phase)
    end do
    max_view = 4.0_dp
    do run = 1, 2
        do frame = 1, frames
            do p = 1, particles
                j = 3*p - 2
                call projection(trajectory(j:j + 2, frame, run), scratch, phase)
                max_view = max(max_view, abs(scratch)*1.05_dp, abs(phase)*1.05_dp)
            end do
        end do
    end do
    display_defect = 1.0e6_dp*defect
    energy_view = max(1.0e-6_dp, 1.15_dp*maxval(abs(display_defect)))
    call figure(figsize=[13.2_dp, 8.8_dp])
    canvas => get_global_figure()
    frame = snapshot_frame
    call draw_frame()
    call savefig('output/microscopic-cloud.png')
    if (.not. snapshot_only) then
        animation = animation_func(update_frame, frames=frames + hold_frames, &
            interval=42, fig=canvas)
        call save_animation(animation, 'output/microscopic-cloud.mp4', fps=fps, &
                            status=encode_status)
        if (encode_status /= 0) error stop 'cloud animation encoding failed'
    end if
    call write_result()
    print '(a,es12.4)', 'Cloud max relative energy defect (diagnostic): ', &
        maxval(abs(defect))
    print '(a,2(1x,i0))', 'Encounter-refined substeps (control, impurity):', refined
    print '(a)', 'Candidate only: trajectory enclosure unavailable; no Debye fit.'
contains
    !> Integral values as i0 (header "Q=32,R=3"), others as g0.
    function number(value) result(text)
        real(dp), intent(in) :: value
        character(len=32) :: text
        if (value == anint(value) .and. abs(value) < 1.0e9_dp) then
            write (text, '(i0)') nint(value)
        else
            write (text, '(g0)') value
        end if
    end function

    subroutine write_result()
        character(len=*), parameter :: none = 'candidate; trajectory/numerical/'// &
            'model/statistical bounds unavailable'
        call result%create(app)
        call result%attribute('run_index', '1=control, 2=fixed central impurity')
        call result%write('time', time, '1', 'frame time (reduced units)', dims='frame')
        call result%write('relative_energy_defect', defect, '1', &
            '(H(t)-H(0))/|H(0)|', dims='frame,run', status='diagnostic')
        call result%write('positions', trajectory(:3*particles, :, :), '1', &
            'positions [x1,y1,z1,x2,...]', dims='position_component,frame,run', &
            status=none)
        call result%write('velocities', trajectory(3*particles + 1:, :, :), '1', &
            'velocities [vx1,vy1,vz1,vx2,...]', dims='velocity_component,frame,run', &
            status=none)
        call result%write('solve_cpu', candidate_cpu, 's', 'candidate solve CPU time', &
            dims='frame,run')
        call result%write('encounter_refined_substeps', real(refined, dp), '1', &
            'substeps split by encounter control', dims='run')
        call result%close()
    end subroutine
    !> Display chart X = x + z/4, Y = y + z/6 of a 3D position.
    pure subroutine projection(r, x, y)
        real(dp), intent(in) :: r(3)
        real(dp), intent(out) :: x, y
        x = r(1) + r(3)/4; y = r(2) + r(3)/6
    end subroutine

    subroutine update_frame(index)
        integer, intent(in) :: index
        if (index < 1 .or. index > frames + hold_frames) error stop 'invalid animation frame'
        frame = min(index, frames)
        call canvas%clear()
        call draw_frame()
        call canvas%set_rendered(.false.)
    end subroutine
    subroutine draw_frame()
        call suptitle('A fixed negative impurity pushes the nearest electrons outward; '// &
            'the radial distribution shifts (exploratory, N = '// &
            trim(number(real(particles, dp)))//')')
        do run = 1, 2
            call subplot(2, 2, run)
            do p = 1, particles
                j = 3*p - 2
                call projection(trajectory(j:j + 2, frame, run), xp(p), yp(p))
            end do
            ! Dotted circle: sphere equator, a neutralizing background but no wall.
            call plot(circle_x, circle_y, linestyle=':', color=grey, linewidth=1.5_dp)
            if (run == 1) then
                call title('Coulomb electrons, control  |  t = '//num(time(frame), 2))
                call scatter(xp, yp, s=36.0_dp, color=blue, marker='o')
            else
                call title('Same start, fixed impurity q = '// &
                    trim(number(impurity_charge))//' (x) at the centre')
                call scatter(xp, yp, s=36.0_dp, color=vermillion, marker='o')
                call scatter([0.0_dp], [0.0_dp], s=90.0_dp, color=black, marker='x')
            end if
            call xlabel('projection X = x + z/4')
            call ylabel('projection Y = y + z/6')
            call xlim(-box_aspect*max_view, box_aspect*max_view)
            call ylim(-max_view, max_view)
        end do
        call subplot(2, 2, 3)
        do run = 1, 2
            do p = 1, particles
                j = 3*p - 2
                radii(p) = norm2(trajectory(j:j + 2, frame, run))
                fraction(p) = real(p, dp)/real(particles, dp)
            end do
            ! Sort radii: the empirical CDF is an observable, not a smooth density.
            do p = 2, particles
                scratch = radii(p); k = p - 1
                do while (k >= 1)
                    if (radii(k) <= scratch) exit
                    radii(k + 1) = radii(k); k = k - 1
                end do
                radii(k + 1) = scratch
            end do
            ! Post-step staircase drawn as a polyline (fortplot step() does not
            ! route to subplots at the pinned revision).
            cdf_x(1) = 0.0_dp; cdf_y(1) = 0.0_dp
            do p = 1, particles
                cdf_x(2*p) = radii(p); cdf_y(2*p) = cdf_y(2*p - 1)
                cdf_x(2*p + 1) = radii(p); cdf_y(2*p + 1) = fraction(p)
            end do
            cdf_x(2*particles + 2) = 5.0_dp; cdf_y(2*particles + 2) = 1.0_dp
            if (run == 1) then
                call plot(cdf_x, cdf_y, color=blue, linewidth=2.0_dp, label='control')
            else
                call plot(cdf_x, cdf_y, color=vermillion, linestyle='--', &
                    linewidth=2.0_dp, label='negative impurity')
            end if
        end do
        call title('Empirical radial CDF (no statistical band)')
        call xlabel('3D radius r')
        call ylabel('fraction of electrons inside r')
        call xlim(0.0_dp, 5.0_dp); call ylim(0.0_dp, 1.05_dp)
        call legend(loc='lower right')
        call subplot(2, 2, 4)
        call plot(time(:frame), display_defect(:frame, 1), color=blue, &
                  linewidth=2.0_dp, label='control')
        call plot(time(:frame), display_defect(:frame, 2), color=vermillion, &
                  linestyle='--', linewidth=2.0_dp, label='negative impurity')
        call title('Energy drift diagnostic (no trajectory bound)')
        call xlabel('t'); call ylabel('$10^6 (H(t) - H(0))/|H(0)|$')
        call xlim(0.0_dp, time(frames)); call ylim(-energy_view, energy_view)
        call legend(loc='upper left')
        call tight_layout()
    end subroutine
end program
