module streaming_moments_demo_callbacks
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_interval, only: interval_t, interval, hull
    use kin6d_streaming_moments, only: streaming_moments_t, &
        initialize_streaming_moments, advance_streaming_moments, &
        sample_streaming_moments, sample_streaming_moments_box, &
        sample_streaming_moments_distribution
    use kin6d_streaming_moments_leaves, only: streaming_moments_initial
    use kin6d_streaming_comparison, only: streaming_state_t, streaming_particles_t, &
        init_streaming, advance_streaming, init_streaming_particles, &
        get_streaming_particles
    use kin6d_streaming_leaves, only: streaming_plot_x, streaming_velocity_coordinate, &
        streaming_frame_time, streaming_characteristic, streaming_display_shift, &
        streaming_particle_midpoint, streaming_particle_glyph_size
    use fortplot, only: figure, figure_t, get_global_figure, subplot, scatter, &
        pcolormesh, plot, fill_between, &
        title, xlabel, ylabel, xlim, ylim, legend, tight_layout, savefig, colorbar, &
        suptitle
    use kin6d_demo_style, only: blue, movie_fps, hold_frames, num
    use fortplot_animation, only: animation_t, animation_func => FuncAnimation, &
        save_animation
    use fortplot_system_runtime, only: create_directory_runtime
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use kin6d_result_io, only: run_t, result_file_t

    use fortplot_animation, only: animation_callback_t
    implicit none
    private
    public :: run_streaming_moments_demo
    character(len=*), parameter :: bound = 'certified-continuous whole-cell '// &
        'moment enclosure (density strictly positive)'
    type :: run_streaming_moments_demo_context_t
        integer, pointer :: bands => null()
        type(streaming_state_t), pointer :: baseline => null()
        real(dp), pointer :: bulk(:) => null()
        type(figure_t), pointer :: canvas => null()
        real(dp), pointer :: complete_s => null()
        real(dp), pointer :: current(:) => null()
        real(dp), pointer :: density(:) => null()
        integer, pointer :: distribution_data => null()
        real(dp), pointer :: duration => null()
        integer, pointer :: evidence => null()
        real(dp), pointer :: f_record(:,:,:,:) => null()
        integer, pointer :: frames => null()
        real(dp), pointer :: glyph(:) => null()
        real(dp), pointer :: image(:,:) => null()
        real(dp), pointer :: image_hi(:,:) => null()
        real(dp), pointer :: image_lo(:,:) => null()
        type(interval_t), pointer :: jcells(:) => null()
        real(dp), pointer :: jhi(:) => null()
        real(dp), pointer :: jlo(:) => null()
        real(dp), pointer :: moment_record(:,:,:) => null()
        type(streaming_moments_t), pointer :: moments => null()
        type(interval_t), pointer :: ncells(:) => null()
        real(dp), pointer :: nhi(:) => null()
        real(dp), pointer :: nlo(:) => null()
        integer, pointer :: np => null()
        integer, pointer :: nv => null()
        integer, pointer :: nx => null()
        logical, pointer :: ok => null()
        real(dp), pointer :: period => null()
        real(dp), pointer :: post_s => null()
        real(dp), pointer :: pv(:) => null()
        real(dp), pointer :: px(:) => null()
        logical, pointer :: recorded(:) => null()
        integer, pointer :: repeats => null()
        real(dp), pointer :: solve_s => null()
        real(dp), pointer :: t => null()
        real(dp), pointer :: times(:) => null()
        type(interval_t), pointer :: ucells(:) => null()
        real(dp), pointer :: uhi(:) => null()
        real(dp), pointer :: ulo(:) => null()
        real(dp), pointer :: v(:) => null()
        real(dp), pointer :: x(:) => null()
        real(dp), pointer :: x0(:) => null()
    end type

    type, extends(animation_callback_t) :: run_streaming_moments_demo_update_t
        type(run_streaming_moments_demo_context_t), pointer :: data => null()
    contains
        procedure :: update => run_streaming_moments_demo_update_binding
    end type

contains
    subroutine run_streaming_moments_demo()
    implicit none
integer , target :: nx, nv, particle_nx, particle_nv
    ! 193 frames span T = 0..8 in 8 s at 24 fps; a one-second hold follows.
integer , target :: np, frames, repeats
real(dp) , target :: duration
    namelist /streaming_moments/ nx, nv, particle_nx, particle_nv, frames, &
        repeats, duration
    type(run_t) :: run
    type(result_file_t) :: result
    integer :: unit, ios
    character(len=256) :: msg
    type(streaming_moments_t) , target :: moments
    type(streaming_state_t) , target :: baseline, particle_source
    type(streaming_particles_t) :: particles
    type(animation_t) :: animation
    type(figure_t), pointer :: canvas
    real(dp), allocatable , target :: x(:), v(:), a0(:, :), b0(:, :)
real(dp) :: vp(2), ap(0:1, 2), bp(0:1, 2)
    real(dp), allocatable , target :: x0(:), pv(:), weight(:), px(:), glyph(:)
    real(dp) , target :: period, measure
    real(dp), allocatable , target :: image(:, :), density(:), current(:), bulk(:)
    real(dp), allocatable , target :: image_lo(:, :), image_hi(:, :)
    real(dp), allocatable , target :: nlo(:), nhi(:), jlo(:), jhi(:), ulo(:), uhi(:)
    real(dp) , target :: t, solve_s, complete_s, post_s
    real(dp) :: ignored_x, ignored_v
    type(interval_t), allocatable , target :: ncells(:), jcells(:), ucells(:)
    ! Result record of every computed frame (NaN where a frame was not drawn).
    real(dp), allocatable , target :: times(:), moment_record(:, :, :), f_record(:, &
        :, :, :)
integer , target :: i, frame, snapshot_frame, status, evidence, bands
    integer , target :: distribution_data
logical , target :: ok, snapshot
    logical, allocatable , target :: recorded(:)
    character(len=64) :: argument
        type(run_streaming_moments_demo_context_t), target :: ctx
        type(run_streaming_moments_demo_update_t), target :: update_operator

        nx = 96
        nv = 64
        particle_nx = 32
        particle_nv = 48
        frames = 193
        repeats = 20
        duration = 8.0_dp
        vp = [-1.0_dp, 1.0_dp]
        snapshot_frame = 41
        snapshot = .false.
    call run%start('streaming_moments_demo')
    ios = 0; msg = ''
    call run%open_config(unit)
    if (run%has_config()) read (unit, nml=streaming_moments, iostat=ios, iomsg=msg)
    call run%close_config(unit, 'streaming_moments', ios, msg)
    if (min(nx, nv, particle_nx, particle_nv, repeats) < 1 .or. frames < 2) &
        error stop 'streaming_moments_demo sizes must be positive (frames >= 2)'
    if (.not. (duration > 0.0_dp)) error stop 'duration must be positive'
    call run%open_resolved(unit)
    write (unit, nml=streaming_moments)
    call run%close_resolved(unit)
    np = particle_nx*particle_nv
    allocate (x(nx + 1), v(nv + 1), a0(0:1, nv + 1), b0(0:1, nv + 1))
    allocate (x0(np), pv(np), weight(np), px(np), glyph(np))
    allocate (image(nv + 1, nx + 1), density(nx + 1), current(nx + 1), bulk(nx + 1))
    allocate (image_lo(nv + 1, nx + 1), image_hi(nv + 1, nx + 1))
    allocate (nlo(nx + 1), nhi(nx + 1), jlo(nx + 1), jhi(nx + 1), ulo(nx + 1), &
        uhi(nx + 1))
    allocate (ncells(nx), jcells(nx), ucells(nx), recorded(frames))
    allocate (times(frames), moment_record(nx + 1, frames, 9), &
        f_record(nv + 1, nx + 1, frames, 3))
    recorded = .false.
    times = ieee_value(1.0_dp, ieee_quiet_nan)
    moment_record = times(1); f_record = times(1)

    if (run%argument_count() > 2) error stop 'usage: demo [--config FILE] [--snapshot [frame]]'
    if (run%argument_count() >= 1) then
        call run%argument(1, argument)
        if (trim(argument) /= '--snapshot') error stop 'expected --snapshot'
        snapshot = .true.
    end if
    if (run%argument_count() == 2) then
        call run%argument(2, argument)
        read (argument, *, iostat=status) snapshot_frame
        if (status /= 0) error stop 'snapshot frame must be integer'
    end if
    if (snapshot_frame < 1 .or. snapshot_frame > frames) &
        error stop 'frame outside 1..frames'
    call create_directory_runtime('output', ok)
    if (.not. ok) error stop 'output directory unavailable'
    do i = 1, nx + 1
        call streaming_plot_x(real(i - 1, dp), real(nx + 1, dp), x(i))
    end do
    period = x(nx + 1)
    do i = 1, nv + 1
        call streaming_velocity_coordinate(real(i - 1, dp), real(nv, dp), v(i))
        call streaming_moments_initial(v(i), a0(0, i), a0(1, i))
    end do
    b0 = 0.0_dp
    call initialize_streaming_moments(moments, v, a0, b0, ok)
    if (.not. ok) error stop 'moment packet initialization failed'
    call init_streaming(baseline, v, a0, b0, ok)
    if (.not. ok) error stop 'matching baseline initialization failed'
    do i = 1, 2
        call streaming_moments_initial(vp(i), ap(0, i), ap(1, i))
    end do
    bp = 0.0_dp
    call init_streaming(particle_source, vp, ap, bp, ok)
    if (.not. ok) error stop 'particle source initialization failed'
    call init_streaming_particles(particle_source, particle_nx, particle_nv, &
        particles, ok)
    if (.not. ok) error stop 'particle initialization failed'
    call get_streaming_particles(particles, x0, pv, weight, ok)
    if (.not. ok) error stop 'particle retrieval failed'
    call streaming_particle_midpoint(0.0_dp, real(particle_nx, dp), 0.0_dp, &
        real(particle_nv, dp), ignored_x, ignored_v, measure)
    do i = 1, np
        call streaming_particle_glyph_size(weight(i), measure, glyph(i))
    end do
    print '(a)', 'Chapter 1: weighted particles -> distribution -> moments'
    print '(a)', 'Model: f_t+v*f_x=0, E=0; x period 2*pi, v[-1,1]; normalized.'
    print '(a)', 'Exact input: f0=(1+0.5*cos(x))*(1+0.5*v); positive finite Fourier/P1.'
    print '(a)', 'Scales [1]: X=2*pi*x/L, V=v/v_ref, T=2*pi*v_ref*t/L, F=2*v_ref*f/n0.'
    print '(a)', 'n=0.5*integral(f dv); j=0.5*integral(v*f dv); u=j/n; all normalized.'
    print '(a)', 'Bands: rigorous whole-cell moments; heatmap is candidate f.'
    print '(a)', 'Heatmap pointwise f bounds: particles-to-moments-distribution.csv.'
    print '(a)', 'No Maxwellian/beam, collisions, feedback, or certified sensitivities.'
    print '(a)', 'Costs: 64-cell primal solve vs moment packet; 20 CPU repeats.'
    open (newunit=evidence, file='output/particles-to-moments-cost.csv', &
        status='replace')
    write (evidence, '(a)') &
        'frame,t,solve_seconds,complete_seconds,extra_seconds,factor,'// &
        'postprocess_seconds'
    open (newunit=bands, file='output/particles-to-moments-bands.csv', status='replace')
    write (bands, '(a)') &
        'frame,t,x,n,n_lower,n_upper,j,j_lower,j_upper,u,u_lower,u_upper'
    open (newunit=distribution_data, &
        file='output/particles-to-moments-distribution.csv', status='replace')
    write (distribution_data, '(a)') 'frame,T,x,v,f,lower,upper'
    call figure(figsize=[13.2_dp, 8.8_dp])
    canvas => get_global_figure()
    frame = 1
    if (snapshot) frame = snapshot_frame
        call bind_run_streaming_moments_demo_context()
    call run_streaming_moments_demo_compute_frame(ctx, frame)
        call bind_run_streaming_moments_demo_context()
    call run_streaming_moments_demo_draw_frame(ctx)
    call savefig('output/particles-to-moments.png')
    if (.not. snapshot) then
        call bind_run_streaming_moments_demo_context()
        animation = animation_func(update_operator, frames=frames + hold_frames, &
            interval=42, &
            fig=canvas)
        call save_animation(animation, 'output/particles-to-moments.mp4', fps=movie_fps, &
            status=status)
        if (status /= 0) error stop 'moments animation encoding failed'
    end if
    close (evidence)
    close (bands)
    close (distribution_data)
    call result%create(run)
    call result%attribute('frame_record', 'frames not computed in this run are NaN')
    call result%write('time', times, '1', 'normalized time T', dims='frame')
    call result%write('x', x, '1', '2*pi*x/L', dims='x')
    call result%write('v', v, '1', 'v/v_ref', dims='v')
    call write_moment(1, 'n', 'density n/n0', 'candidate')
    call write_moment(2, 'n_lower', 'density lower bound', bound)
    call write_moment(3, 'n_upper', 'density upper bound', bound)
    call write_moment(4, 'j', 'current j', 'candidate')
    call write_moment(5, 'j_lower', 'current lower bound', bound)
    call write_moment(6, 'j_upper', 'current upper bound', bound)
    call write_moment(7, 'u', 'bulk velocity u = j/n', 'candidate')
    call write_moment(8, 'u_lower', 'bulk velocity lower bound', bound)
    call write_moment(9, 'u_upper', 'bulk velocity upper bound', bound)
    call result%write('f', f_record(:, :, :, 1), '1', '2*v_ref*f/n0', &
        dims='v,x,frame', status='candidate')
    call result%write('f_lower', f_record(:, :, :, 2), '1', 'pointwise f lower bound', &
        dims='v,x,frame', status='certified-continuous pointwise bound')
    call result%write('f_upper', f_record(:, :, :, 3), '1', 'pointwise f upper bound', &
        dims='v,x,frame', status='certified-continuous pointwise bound')
    call result%close()
contains
    subroutine bind_run_streaming_moments_demo_context()
        ctx%bands => bands
        ctx%baseline => baseline
        nullify(ctx%bulk)
        if (allocated(bulk)) ctx%bulk => bulk
        ctx%canvas => canvas
        ctx%complete_s => complete_s
        nullify(ctx%current)
        if (allocated(current)) ctx%current => current
        nullify(ctx%density)
        if (allocated(density)) ctx%density => density
        ctx%distribution_data => distribution_data
        ctx%duration => duration
        ctx%evidence => evidence
        nullify(ctx%f_record)
        if (allocated(f_record)) ctx%f_record => f_record
        ctx%frames => frames
        nullify(ctx%glyph)
        if (allocated(glyph)) ctx%glyph => glyph
        nullify(ctx%image)
        if (allocated(image)) ctx%image => image
        nullify(ctx%image_hi)
        if (allocated(image_hi)) ctx%image_hi => image_hi
        nullify(ctx%image_lo)
        if (allocated(image_lo)) ctx%image_lo => image_lo
        nullify(ctx%jcells)
        if (allocated(jcells)) ctx%jcells => jcells
        nullify(ctx%jhi)
        if (allocated(jhi)) ctx%jhi => jhi
        nullify(ctx%jlo)
        if (allocated(jlo)) ctx%jlo => jlo
        nullify(ctx%moment_record)
        if (allocated(moment_record)) ctx%moment_record => moment_record
        ctx%moments => moments
        nullify(ctx%ncells)
        if (allocated(ncells)) ctx%ncells => ncells
        nullify(ctx%nhi)
        if (allocated(nhi)) ctx%nhi => nhi
        nullify(ctx%nlo)
        if (allocated(nlo)) ctx%nlo => nlo
        ctx%np => np
        ctx%nv => nv
        ctx%nx => nx
        ctx%ok => ok
        ctx%period => period
        ctx%post_s => post_s
        nullify(ctx%pv)
        if (allocated(pv)) ctx%pv => pv
        nullify(ctx%px)
        if (allocated(px)) ctx%px => px
        nullify(ctx%recorded)
        if (allocated(recorded)) ctx%recorded => recorded
        ctx%repeats => repeats
        ctx%solve_s => solve_s
        ctx%t => t
        nullify(ctx%times)
        if (allocated(times)) ctx%times => times
        nullify(ctx%ucells)
        if (allocated(ucells)) ctx%ucells => ucells
        nullify(ctx%uhi)
        if (allocated(uhi)) ctx%uhi => uhi
        nullify(ctx%ulo)
        if (allocated(ulo)) ctx%ulo => ulo
        nullify(ctx%v)
        if (allocated(v)) ctx%v => v
        nullify(ctx%x)
        if (allocated(x)) ctx%x => x
        nullify(ctx%x0)
        if (allocated(x0)) ctx%x0 => x0
        update_operator%data => ctx
    end subroutine bind_run_streaming_moments_demo_context

    subroutine write_moment(k, name, long_name, status)
        integer, intent(in) :: k
        character(len=*), intent(in) :: name, long_name, status
        call result%write(name, moment_record(:, :, k), '1', long_name, &
            dims='x,frame', status=status)
    end subroutine

    end subroutine run_streaming_moments_demo

    subroutine run_streaming_moments_demo_compute_frame(ctx, number)
        type(run_streaming_moments_demo_context_t), intent(in) :: ctx
        integer, intent(in) :: number
        real(dp) :: start, finish, value, lo, hi, unwrapped, factor
        type(interval_t) :: nb, jb, ub
        integer :: ix, iv, p, repeat, turns, left, right
        integer, pointer :: local_bands
        integer, pointer :: local_distribution_data
        integer, pointer :: local_evidence
        local_bands => ctx%bands
        local_distribution_data => ctx%distribution_data
        local_evidence => ctx%evidence

        call streaming_frame_time(real(number - 1, dp), real(ctx%frames, dp), &
            ctx%duration, ctx%t)
        call cpu_time(start)
        do repeat = 1, ctx%repeats
            call advance_streaming(ctx%baseline, ctx%t, ctx%ok, certify=.false.)
        end do
        call cpu_time(finish)
        ctx%solve_s = (finish - start)/ctx%repeats
        if (.not. ctx%ok) error stop 'primal baseline failed'
        call cpu_time(start)
        do repeat = 1, ctx%repeats
            call advance_streaming_moments(ctx%moments, ctx%t, ctx%ok)
        end do
        call cpu_time(finish)
        ctx%complete_s = (finish - start)/ctx%repeats
        if (.not. ctx%ok) error stop 'moment certification failed'
        call cpu_time(start)
        do p = 1, ctx%np
            call streaming_characteristic(ctx%x0(p), ctx%pv(p), ctx%t, unwrapped)
            turns = floor(unwrapped/ctx%period)
            call streaming_display_shift(unwrapped, real(turns, dp), ctx%px(p))
        end do
        do ix = 1, ctx%nx + 1
            do iv = 1, ctx%nv + 1
                call sample_streaming_moments_distribution(ctx%moments, ctx%x(ix), &
                    ctx%v(iv), &
                    value, lo, hi, ctx%ok)
                if (.not. ctx%ok) error stop 'distribution sampling failed'
                ctx%image(iv, ix) = value
                ctx%image_lo(iv, ix) = lo
                ctx%image_hi(iv, ix) = hi
            end do
            call sample_streaming_moments(ctx%moments, ctx%x(ix), ctx%density(ix), &
                ctx%current(ix), &
                ctx%bulk(ix), nb, jb, ub, ctx%ok)
            if (.not. ctx%ok) error stop 'positive-density point moments unavailable'
            ctx%nlo(ix) = nb%lo; ctx%nhi(ix) = nb%hi
            ctx%jlo(ix) = jb%lo; ctx%jhi(ix) = jb%hi
            ctx%ulo(ix) = ub%lo; ctx%uhi(ix) = ub%hi
        end do
        do ix = 1, ctx%nx
            call sample_streaming_moments_box(ctx%moments, interval(ctx%x(ix), &
                ctx%x(ix + 1)), &
                ctx%ncells(ix), ctx%jcells(ix), ctx%ucells(ix), ctx%ok)
            if (.not. ctx%ok) error stop 'whole-cell moment bands unavailable'
        end do
        ! Each vertex covers both adjacent complete cells. Every convex
        ! combination of its lower/upper endpoints therefore covers that cell.
        do ix = 1, ctx%nx + 1
            nb = interval(ctx%nlo(ix), ctx%nhi(ix))
            jb = interval(ctx%jlo(ix), ctx%jhi(ix))
            ub = interval(ctx%ulo(ix), ctx%uhi(ix))
            left = max(1, ix - 1)
            right = min(ctx%nx, ix)
            nb = hull(nb, hull(ctx%ncells(left), ctx%ncells(right)))
            jb = hull(jb, hull(ctx%jcells(left), ctx%jcells(right)))
            ub = hull(ub, hull(ctx%ucells(left), ctx%ucells(right)))
            ctx%nlo(ix) = nb%lo; ctx%nhi(ix) = nb%hi
            ctx%jlo(ix) = jb%lo; ctx%jhi(ix) = jb%hi
            ctx%ulo(ix) = ub%lo; ctx%uhi(ix) = ub%hi
        end do
        call cpu_time(finish)
        ctx%post_s = finish - start
        if (.not. ctx%recorded(number)) then
            factor = -1.0_dp
            if (ctx%solve_s > 0.0_dp) factor = ctx%complete_s/ctx%solve_s
            write (local_evidence, '(i0,6(",",es24.16e3))') number, ctx%t, &
                ctx%solve_s, ctx%complete_s, &
                ctx%complete_s - ctx%solve_s, factor, ctx%post_s
            do ix = 1, ctx%nx + 1
                write (local_bands, '(i0,11(",",es24.16e3))') number, ctx%t, &
                    ctx%x(ix), &
                    ctx%density(ix), ctx%nlo(ix), ctx%nhi(ix), ctx%current(ix), &
                        ctx%jlo(ix), ctx%jhi(ix), &
                    ctx%bulk(ix), ctx%ulo(ix), ctx%uhi(ix)
                do iv = 1, ctx%nv + 1
                    write (local_distribution_data, '(i0,6(",",es24.16e3))') &
                        number, ctx%t, ctx%x(ix), ctx%v(iv), ctx%image(iv, ix), &
                        ctx%image_lo(iv, ix), ctx%image_hi(iv, ix)
                end do
            end do
            ctx%recorded(number) = .true.
            ctx%times(number) = ctx%t
            ctx%moment_record(:, number, :) = reshape([ctx%density, ctx%nlo, &
                ctx%nhi, ctx%current, ctx%jlo, &
                ctx%jhi, ctx%bulk, ctx%ulo, ctx%uhi], [ctx%nx + 1, 9])
            ctx%f_record(:, :, number, 1) = ctx%image
            ctx%f_record(:, :, number, 2) = ctx%image_lo
            ctx%f_record(:, :, number, 3) = ctx%image_hi
        end if
    end subroutine run_streaming_moments_demo_compute_frame

    subroutine run_streaming_moments_demo_draw_frame(ctx)
        type(run_streaming_moments_demo_context_t), intent(in) :: ctx
        call suptitle('Velocity moments compress f(x,v) into density and bulk '// &
            'velocity: they follow the flow but discard the velocity shape')
        call subplot(2, 2, 1)
        call title('Weighted particles (marker area ~ weight)  |  T = '//num(ctx%t, 2))
        call scatter(ctx%px, ctx%pv, s=0.6_dp*ctx%glyph, color=blue, alpha=0.5_dp)
        call xlabel('$2\pi x/L$'); call ylabel('$v/v_{ref}$')
        call xlim(0.0_dp, ctx%period); call ylim(-1.0_dp, 1.0_dp)
        call subplot(2, 2, 2)
        call title('Eulerian Fourier/P1 distribution (candidate)')
        call pcolormesh(ctx%x, ctx%v, ctx%image, cmap='viridis', vmin=0.25_dp, &
            vmax=2.25_dp)
        call colorbar(label='$2 v_{ref} f/n_0$')
        call xlabel('$2\pi x/L$'); call ylabel('$v/v_{ref}$')
        call xlim(0.0_dp, ctx%period); call ylim(-1.0_dp, 1.0_dp)
        call subplot(2, 2, 3)
        call title('Density: zeroth velocity moment')
        call fill_between(ctx%x, ctx%nlo, ctx%nhi, color=blue, alpha=0.25_dp)
        call plot(ctx%x, ctx%density, color=blue, linewidth=2.0_dp, &
            label='n ± whole-cell bound')
        call xlabel('$2\pi x/L$'); call ylabel('$n/n_0$')
        call xlim(0.0_dp, ctx%period); call ylim(0.35_dp, 1.65_dp)
        call legend(loc='upper center')
        call subplot(2, 2, 4)
        call title('Bulk velocity u = j/n')
        call fill_between(ctx%x, ctx%ulo, ctx%uhi, color=blue, alpha=0.25_dp)
        call plot(ctx%x, ctx%bulk, color=blue, linewidth=2.0_dp, &
            label='u ± whole-cell bound')
        call xlabel('$2\pi x/L$')
        call ylabel('$u/v_{ref}$')
        call xlim(0.0_dp, ctx%period); call ylim(-0.15_dp, 0.40_dp)
        call legend(loc='upper right')
        call tight_layout()
    end subroutine run_streaming_moments_demo_draw_frame

    subroutine run_streaming_moments_demo_update(ctx, index)
        type(run_streaming_moments_demo_context_t), intent(in) :: ctx
        integer, intent(in) :: index
        call run_streaming_moments_demo_compute_frame(ctx, min(index, ctx%frames))
        call ctx%canvas%clear()
        call run_streaming_moments_demo_draw_frame(ctx)
        call ctx%canvas%set_rendered(.false.)
    end subroutine run_streaming_moments_demo_update

    subroutine run_streaming_moments_demo_update_binding(self, frame)
        class(run_streaming_moments_demo_update_t), intent(inout) :: self
        integer, intent(in) :: frame
        call run_streaming_moments_demo_update(self%data, frame)
    end subroutine
end module streaming_moments_demo_callbacks

program streaming_moments_demo
    use streaming_moments_demo_callbacks, only: run_streaming_moments_demo
    implicit none
    call run_streaming_moments_demo()
end program streaming_moments_demo
