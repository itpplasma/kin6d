module microscopic_cloud_demo_callbacks
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

    use fortplot_animation, only: animation_callback_t
    implicit none
    private
    public :: run_microscopic_cloud_demo
    integer, parameter :: fps = 24
    real(dp), parameter :: two_pi = 6.2831853071795864769_dp
    real(dp), parameter :: blue(3) = [0.0_dp, 0.4471_dp, 0.6980_dp]
    real(dp), parameter :: box_aspect = 1.594_dp
    type :: run_microscopic_cloud_demo_context_t
        type(figure_t), pointer :: canvas => null()
        real(dp), pointer :: cdf_x(:) => null()
        real(dp), pointer :: cdf_y(:) => null()
        real(dp), pointer :: circle_x(:) => null()
        real(dp), pointer :: circle_y(:) => null()
        real(dp), pointer :: display_defect(:,:) => null()
        real(dp), pointer :: energy_view => null()
        real(dp), pointer :: fraction(:) => null()
        integer, pointer :: frame => null()
        integer, pointer :: frames => null()
        real(dp), pointer :: impurity_charge => null()
        integer, pointer :: j => null()
        integer, pointer :: k => null()
        real(dp), pointer :: max_view => null()
        integer, pointer :: p => null()
        integer, pointer :: particles => null()
        real(dp), pointer :: radii(:) => null()
        type(result_file_t), pointer :: result => null()
        integer, pointer :: run => null()
        real(dp), pointer :: scratch => null()
        real(dp), pointer :: time(:) => null()
        real(dp), pointer :: trajectory(:,:,:) => null()
        real(dp), pointer :: xp(:) => null()
        real(dp), pointer :: yp(:) => null()
    end type

    type, extends(animation_callback_t) :: run_microscopic_cloud_demo_update_frame_t
        type(run_microscopic_cloud_demo_context_t), pointer :: data => null()
    contains
        procedure :: update => run_microscopic_cloud_demo_update_frame_binding
    end type

contains
    subroutine run_microscopic_cloud_demo()
    implicit none
    ! Namelist &microscopic_cloud: electrons, output frames, substeps of dt per
    ! frame, encounter control eta, fixed central impurity charge of run 2 and
    ! the neutralizing sphere (defaults from kin6d_microscopic_cloud).
    ! 241 frames of 50 substeps: the same 12000 steps as before, shown in 10 s.
integer , target :: particles, frames, substeps
real(dp) , target :: dt, encounter_eta, impurity_charge
    real(dp) :: sphere_charge, sphere_radius
    namelist /microscopic_cloud/ particles, frames, substeps, dt, encounter_eta, &
        impurity_charge, sphere_charge, sphere_radius
    type(run_t) :: app
    type(result_file_t) , target :: result
    integer :: ios
    character(len=256) :: msg
    ! Width/height of the cloud axes boxes in pixels (measured), for equal aspect.
    type(cloud_model_t) :: model(2)
    type(nbody_system_t) :: system
    type(fortnum_status_t) :: status
    type(figure_t), pointer :: canvas
    type(animation_t) :: animation
    real(dp), allocatable , target :: state(:), trajectory(:, :, :)
    real(dp), allocatable , target :: time(:), defect(:, :), display_defect(:, :)
    real(dp) :: initial_energy, energy
    real(dp), allocatable :: candidate_cpu(:, :), diagnostic_cpu(:, :)
    real(dp) :: start, finish
    real(dp), allocatable , target :: xp(:), yp(:), radii(:), fraction(:), cdf_x(:), &
        cdf_y(:)
    real(dp) , target :: circle_x(101), circle_y(101), phase, scratch, max_view, &
        energy_view
    real(dp) :: limit
    integer , target :: run, frame, step, p, j, k, unit, snapshot_frame, encode_status
    integer :: substep, pieces, refined(2)
    logical :: ok, snapshot_only
    character(len=64) :: argument
        type(run_microscopic_cloud_demo_context_t), target :: ctx
        type(run_microscopic_cloud_demo_update_frame_t), target :: update_frame_operator
        particles = 32
        frames = 241
        substeps = 50
        dt = 0.001_dp
        encounter_eta = 0.05_dp
        impurity_charge = -2.0_dp
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
        call bind_run_microscopic_cloud_demo_context()
    write (unit, '(a)') '# q=-1,m=1,Coulomb prefactor=1; sphere Q='// &
        trim(run_microscopic_cloud_demo_number(ctx, &
            sphere_charge))//',R='//trim(run_microscopic_cloud_demo_number(ctx, &
            sphere_radius))
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
        call bind_run_microscopic_cloud_demo_context()
                call run_microscopic_cloud_demo_projection(trajectory(j:j + 2, &
                    frame, run), scratch, phase)
                max_view = max(max_view, abs(scratch)*1.05_dp, abs(phase)*1.05_dp)
            end do
        end do
    end do
    display_defect = 1.0e6_dp*defect
    energy_view = max(1.0e-6_dp, 1.15_dp*maxval(abs(display_defect)))
    call figure(figsize=[13.2_dp, 8.8_dp])
    canvas => get_global_figure()
    frame = snapshot_frame
        call bind_run_microscopic_cloud_demo_context()
    call run_microscopic_cloud_demo_draw_frame(ctx)
    call savefig('output/microscopic-cloud.png')
    if (.not. snapshot_only) then
        call bind_run_microscopic_cloud_demo_context()
        animation = animation_func(update_frame_operator, frames=frames + hold_frames, &
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
    subroutine bind_run_microscopic_cloud_demo_context()
        ctx%canvas => canvas
        nullify(ctx%cdf_x)
        if (allocated(cdf_x)) ctx%cdf_x => cdf_x
        nullify(ctx%cdf_y)
        if (allocated(cdf_y)) ctx%cdf_y => cdf_y
        ctx%circle_x => circle_x
        ctx%circle_y => circle_y
        nullify(ctx%display_defect)
        if (allocated(display_defect)) ctx%display_defect => display_defect
        ctx%energy_view => energy_view
        nullify(ctx%fraction)
        if (allocated(fraction)) ctx%fraction => fraction
        ctx%frame => frame
        ctx%frames => frames
        ctx%impurity_charge => impurity_charge
        ctx%j => j
        ctx%k => k
        ctx%max_view => max_view
        ctx%p => p
        ctx%particles => particles
        nullify(ctx%radii)
        if (allocated(radii)) ctx%radii => radii
        ctx%result => result
        ctx%run => run
        ctx%scratch => scratch
        nullify(ctx%time)
        if (allocated(time)) ctx%time => time
        nullify(ctx%trajectory)
        if (allocated(trajectory)) ctx%trajectory => trajectory
        nullify(ctx%xp)
        if (allocated(xp)) ctx%xp => xp
        nullify(ctx%yp)
        if (allocated(yp)) ctx%yp => yp
        update_frame_operator%data => ctx
    end subroutine bind_run_microscopic_cloud_demo_context

    !> Integral values as i0 (header "Q=32,R=3"), others as g0.

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


    end subroutine run_microscopic_cloud_demo

    function run_microscopic_cloud_demo_number(ctx, value) result(text)
        type(run_microscopic_cloud_demo_context_t), intent(in) :: ctx
        real(dp), intent(in) :: value
        character(len=32) :: text
        if (value == anint(value) .and. abs(value) < 1.0e9_dp) then
            write (text, '(i0)') nint(value)
        else
            write (text, '(g0)') value
        end if
    end function run_microscopic_cloud_demo_number

    pure subroutine run_microscopic_cloud_demo_projection(r, x, y)
        real(dp), intent(in) :: r(3)
        real(dp), intent(out) :: x, y
        x = r(1) + r(3)/4; y = r(2) + r(3)/6
    end subroutine run_microscopic_cloud_demo_projection

    subroutine run_microscopic_cloud_demo_update_frame(ctx, index)
        type(run_microscopic_cloud_demo_context_t), intent(in) :: ctx
        integer, intent(in) :: index
        if (index < 1 .or. index > ctx%frames + hold_frames) error stop &
            'invalid animation frame'
        ctx%frame = min(index, ctx%frames)
        call ctx%canvas%clear()
        call run_microscopic_cloud_demo_draw_frame(ctx)
        call ctx%canvas%set_rendered(.false.)
    end subroutine run_microscopic_cloud_demo_update_frame

    subroutine run_microscopic_cloud_demo_update_frame_binding(self, frame)
        class(run_microscopic_cloud_demo_update_frame_t), intent(inout) :: self
        integer, intent(in) :: frame
        call run_microscopic_cloud_demo_update_frame(self%data, frame)
    end subroutine

    subroutine run_microscopic_cloud_demo_draw_frame(ctx)
        type(run_microscopic_cloud_demo_context_t), intent(in) :: ctx
        integer, pointer :: local_p
        integer, pointer :: local_run
        local_p => ctx%p
        local_run => ctx%run
        call suptitle('A fixed negative impurity pushes the nearest electrons outward; '// &
            'the radial distribution shifts (exploratory, N = '// &
            trim(run_microscopic_cloud_demo_number(ctx, real(ctx%particles, dp)))//')')
        do local_run = 1, 2
            call subplot(2, 2, local_run)
            do local_p = 1, ctx%particles
                ctx%j = 3*local_p - 2
                call &
                    run_microscopic_cloud_demo_projection(ctx%trajectory(ctx%j:ctx%j &
                    + 2, ctx%frame, local_run), ctx%xp(local_p), ctx%yp(local_p))
            end do
            ! Dotted circle: sphere equator, a neutralizing background but no wall.
            call plot(ctx%circle_x, ctx%circle_y, linestyle=':', color=grey, &
                linewidth=1.5_dp)
            if (local_run == 1) then
                call &
                    title('Coulomb electrons, control  |  t = '//num(ctx%time(ctx%frame), 2))
                call scatter(ctx%xp, ctx%yp, s=36.0_dp, color=blue, marker='o')
            else
                call title('Same start, fixed impurity q = '// &
                    trim(run_microscopic_cloud_demo_number(ctx, &
                        ctx%impurity_charge))//' (x) at the centre')
                call scatter(ctx%xp, ctx%yp, s=36.0_dp, color=vermillion, marker='o')
                call scatter([0.0_dp], [0.0_dp], s=90.0_dp, color=black, marker='x')
            end if
            call xlabel('projection X = x + z/4')
            call ylabel('projection Y = y + z/6')
            call xlim(-box_aspect*ctx%max_view, box_aspect*ctx%max_view)
            call ylim(-ctx%max_view, ctx%max_view)
        end do
        call subplot(2, 2, 3)
        do local_run = 1, 2
            do local_p = 1, ctx%particles
                ctx%j = 3*local_p - 2
                ctx%radii(local_p) = norm2(ctx%trajectory(ctx%j:ctx%j + 2, &
                    ctx%frame, local_run))
                ctx%fraction(local_p) = real(local_p, dp)/real(ctx%particles, dp)
            end do
            ! Sort radii: the empirical CDF is an observable, not a smooth density.
            do local_p = 2, ctx%particles
                ctx%scratch = ctx%radii(local_p); ctx%k = local_p - 1
                do while (ctx%k >= 1)
                    if (ctx%radii(ctx%k) <= ctx%scratch) exit
                    ctx%radii(ctx%k + 1) = ctx%radii(ctx%k); ctx%k = ctx%k - 1
                end do
                ctx%radii(ctx%k + 1) = ctx%scratch
            end do
            ! Post-step staircase drawn as a polyline (fortplot step() does not
            ! route to subplots at the pinned revision).
            ctx%cdf_x(1) = 0.0_dp; ctx%cdf_y(1) = 0.0_dp
            do local_p = 1, ctx%particles
                ctx%cdf_x(2*local_p) = ctx%radii(local_p); ctx%cdf_y(2*local_p) = &
                    ctx%cdf_y(2*local_p - 1)
                ctx%cdf_x(2*local_p + 1) = ctx%radii(local_p); ctx%cdf_y(2*local_p + &
                    1) = ctx%fraction(local_p)
            end do
            ctx%cdf_x(2*ctx%particles + 2) = 5.0_dp; ctx%cdf_y(2*ctx%particles + 2) &
                = 1.0_dp
            if (local_run == 1) then
                call plot(ctx%cdf_x, ctx%cdf_y, color=blue, linewidth=2.0_dp, &
                    label='control')
            else
                call plot(ctx%cdf_x, ctx%cdf_y, color=vermillion, linestyle='--', &
                    linewidth=2.0_dp, label='negative impurity')
            end if
        end do
        call title('Empirical radial CDF (no statistical band)')
        call xlabel('3D radius r')
        call ylabel('fraction of electrons inside r')
        call xlim(0.0_dp, 5.0_dp); call ylim(0.0_dp, 1.05_dp)
        call legend(loc='lower right')
        call subplot(2, 2, 4)
        call plot(ctx%time(:ctx%frame), ctx%display_defect(:ctx%frame, 1), color=blue, &
                  linewidth=2.0_dp, label='control')
        call plot(ctx%time(:ctx%frame), ctx%display_defect(:ctx%frame, 2), &
            color=vermillion, &
                  linestyle='--', linewidth=2.0_dp, label='negative impurity')
        call title('Energy drift diagnostic (no trajectory bound)')
        call xlabel('t'); call ylabel('$10^6 (H(t) - H(0))/|H(0)|$')
        call xlim(0.0_dp, ctx%time(ctx%frames)); call ylim(-ctx%energy_view, &
            ctx%energy_view)
        call legend(loc='upper left')
        call tight_layout()
    end subroutine run_microscopic_cloud_demo_draw_frame
end module microscopic_cloud_demo_callbacks

program microscopic_cloud_demo
    use microscopic_cloud_demo_callbacks, only: run_microscopic_cloud_demo
    implicit none
    call run_microscopic_cloud_demo()
end program microscopic_cloud_demo
