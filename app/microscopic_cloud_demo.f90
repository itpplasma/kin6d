program microscopic_cloud_demo
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_status, only: fortnum_status_t, FORTNUM_OK
    use kin6d_microscopic_cloud, only: cloud_model_t, cloud_system, cloud_initialize
    use kin6d_nbody, only: nbody_system_t, nbody_step, nbody_energy, &
                           nbody_encounter_step
    use fortplot, only: figure, figure_t, get_global_figure, subplot, scatter, &
                  plot, title, xlabel, ylabel, xlim, ylim, legend, tight_layout, savefig
    use fortplot_animation, only: animation_t, animation_func => FuncAnimation, &
                                  save_animation
    use fortplot_system_runtime, only: create_directory_runtime
    use kin6d_result_io, only: run_t, result_file_t
    implicit none
    integer, parameter :: fps = 24
    ! Namelist &microscopic_cloud: electrons, output frames, substeps of dt per
    ! frame, encounter control eta, fixed central impurity charge of run 2 and
    ! the neutralizing sphere (defaults from kin6d_microscopic_cloud).
    integer :: particles = 32, frames = 121, substeps = 100
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
    real(dp), parameter :: orange(3) = [0.8353_dp, 0.3686_dp, 0.0_dp]
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
    real(dp), allocatable :: xp(:), yp(:), radii(:), fraction(:)
    real(dp) :: circle_x(101), circle_y(101), phase, scratch, max_view, energy_view
    real(dp) :: limit
    integer :: run, frame, step, p, j, k, unit, snapshot_frame, encode_status
    integer :: substep, pieces, refined(2)
    logical :: ok, snapshot_only
    character(len=64) :: argument
    character(len=150) :: heading
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
        fraction(particles))
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
    energy_view = max(0.1_dp, 1.4_dp*maxval(abs(display_defect)))
    call figure(figsize=[13.2_dp, 8.8_dp])
    canvas => get_global_figure()
    frame = snapshot_frame
    call draw_frame()
    call savefig('output/microscopic-cloud.png')
    if (.not. snapshot_only) then
        animation = animation_func(update_frame, frames=frames, interval=42, fig=canvas)
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
        if (index < 1 .or. index > frames) error stop 'invalid animation frame'
        frame = index
        call canvas%clear()
        call draw_frame()
        call canvas%set_rendered(.false.)
    end subroutine
    subroutine draw_frame()
        do run = 1, 2
            call subplot(2, 2, run)
            do p = 1, particles
                j = 3*p - 2
                call projection(trajectory(j:j + 2, frame, run), xp(p), yp(p))
            end do
            if (run == 1) then
                write (heading, '(i0,a,f5.2)') particles, &
                    ' Coulomb electrons | control | t = ', time(frame)
                call scatter(xp, yp, s=36.0_dp, color=blue, marker='o')
            else
                write (heading, '(a,f5.2)') 'Same start + negative impurity | t = ', &
                    time(frame)
                call scatter(xp, yp, s=36.0_dp, color=orange, marker='o')
                call scatter([0.0_dp], [0.0_dp], s=90.0_dp, &
                         color=[0.0_dp, 0.0_dp, 0.0_dp], marker='x', &
                         label='fixed q='//trim(number(impurity_charge)))
                call legend(loc='upper right')
            end if
            call plot(circle_x, circle_y, linestyle=':', &
                    color=[0.55_dp, 0.55_dp, 0.55_dp], label='sphere equator (no wall)')
            call legend(loc='upper right')
            call title(trim(heading))
            call xlabel('projection X = x + z/4 [1]')
            call ylabel('projection Y = y + z/6 [1]')
            call xlim(-max_view, max_view); call ylim(-max_view, max_view)
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
            if (run == 1) then
                call scatter(radii, fraction, s=18.0_dp, color=blue, &
                             marker='o', label='control')
            else
                call scatter(radii, fraction, s=18.0_dp, color=orange, &
                             marker='s', label='negative impurity')
            end if
        end do
        call title('Empirical radial CDF | candidate; no statistical band')
        call xlabel('3D radius r [1]')
        call ylabel('fraction of electrons inside r [1]')
        call xlim(0.0_dp, 5.0_dp); call ylim(0.0_dp, 1.15_dp)
        call legend(loc='upper left')
        call subplot(2, 2, 4)
        call plot(time(:frame), display_defect(:frame, 1), color=blue, label='control')
        call plot(time(:frame), display_defect(:frame, 2), color=orange, &
                  linestyle='--', &
                  label='negative impurity')
        call title('Energy diagnostic | trajectory bounds unavailable')
        call xlabel('t [1]'); call ylabel('10^6 (H(t)-H(0))/|H(0)| [1]')
        call xlim(0.0_dp, time(frames)); call ylim(-energy_view, energy_view)
        call legend(loc='upper right')
        call tight_layout()
    end subroutine
end program
