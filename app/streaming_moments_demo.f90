program streaming_moments_demo
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
        title, xlabel, ylabel, xlim, ylim, legend, tight_layout, savefig
    use fortplot_animation, only: animation_t, animation_func => FuncAnimation, &
        save_animation
    use fortplot_system_runtime, only: create_directory_runtime
    implicit none
    integer, parameter :: nx = 96, nv = 64, particle_nx = 32, particle_nv = 48
    integer, parameter :: np = particle_nx*particle_nv, frames = 81, repeats = 20
    real(dp), parameter :: duration = 8.0_dp
    type(streaming_moments_t) :: moments
    type(streaming_state_t) :: baseline, particle_source
    type(streaming_particles_t) :: particles
    type(animation_t) :: animation
    type(figure_t), pointer :: canvas
    real(dp) :: x(nx + 1), v(nv + 1), a0(0:1, nv + 1), b0(0:1, nv + 1)
    real(dp) :: vp(2) = [-1.0_dp, 1.0_dp], ap(0:1, 2), bp(0:1, 2)
    real(dp) :: x0(np), pv(np), weight(np), px(np), glyph(np), period, measure
    real(dp) :: image(nv + 1, nx + 1), density(nx + 1), current(nx + 1), bulk(nx + 1)
    real(dp) :: image_lo(nv + 1, nx + 1), image_hi(nv + 1, nx + 1)
    real(dp) :: nlo(nx + 1), nhi(nx + 1), jlo(nx + 1), jhi(nx + 1)
    real(dp) :: ulo(nx + 1), uhi(nx + 1), t, solve_s, complete_s, post_s
    real(dp) :: ignored_x, ignored_v
    type(interval_t) :: ncells(nx), jcells(nx), ucells(nx)
    integer :: i, frame, snapshot_frame = 41, status, evidence, bands
    integer :: distribution_data
    logical :: ok, snapshot = .false., recorded(frames) = .false.
    character(len=64) :: argument

    if (command_argument_count() > 2) error stop 'usage: demo [--snapshot [frame]]'
    if (command_argument_count() >= 1) then
        call get_command_argument(1, argument)
        if (trim(argument) /= '--snapshot') error stop 'expected --snapshot'
        snapshot = .true.
    end if
    if (command_argument_count() == 2) then
        call get_command_argument(2, argument)
        read (argument, *, iostat=status) snapshot_frame
        if (status /= 0) error stop 'snapshot frame must be integer'
    end if
    if (snapshot_frame < 1 .or. snapshot_frame > frames) &
        error stop 'frame outside 1..81'
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
    call figure(figsize=[12.0_dp, 8.0_dp])
    canvas => get_global_figure()
    frame = 1
    if (snapshot) frame = snapshot_frame
    call compute_frame(frame)
    call draw_frame()
    call savefig('output/particles-to-moments.png')
    if (.not. snapshot) then
        animation = animation_func(update, frames=frames, interval=62, fig=canvas)
        call save_animation(animation, 'output/particles-to-moments.mp4', fps=16, &
            status=status)
        if (status /= 0) error stop 'moments animation encoding failed'
    end if
    close (evidence)
    close (bands)
    close (distribution_data)
contains
    subroutine compute_frame(number)
        integer, intent(in) :: number
        real(dp) :: start, finish, value, lo, hi, unwrapped, factor
        type(interval_t) :: nb, jb, ub
        integer :: ix, iv, p, repeat, turns, left, right

        call streaming_frame_time(real(number - 1, dp), real(frames, dp), duration, t)
        call cpu_time(start)
        do repeat = 1, repeats
            call advance_streaming(baseline, t, ok, certify=.false.)
        end do
        call cpu_time(finish)
        solve_s = (finish - start)/repeats
        if (.not. ok) error stop 'primal baseline failed'
        call cpu_time(start)
        do repeat = 1, repeats
            call advance_streaming_moments(moments, t, ok)
        end do
        call cpu_time(finish)
        complete_s = (finish - start)/repeats
        if (.not. ok) error stop 'moment certification failed'
        call cpu_time(start)
        do p = 1, np
            call streaming_characteristic(x0(p), pv(p), t, unwrapped)
            turns = floor(unwrapped/period)
            call streaming_display_shift(unwrapped, real(turns, dp), px(p))
        end do
        do ix = 1, nx + 1
            do iv = 1, nv + 1
                call sample_streaming_moments_distribution(moments, x(ix), v(iv), &
                    value, lo, hi, ok)
                if (.not. ok) error stop 'distribution sampling failed'
                image(iv, ix) = value
                image_lo(iv, ix) = lo
                image_hi(iv, ix) = hi
            end do
            call sample_streaming_moments(moments, x(ix), density(ix), current(ix), &
                bulk(ix), nb, jb, ub, ok)
            if (.not. ok) error stop 'positive-density point moments unavailable'
            nlo(ix) = nb%lo; nhi(ix) = nb%hi
            jlo(ix) = jb%lo; jhi(ix) = jb%hi
            ulo(ix) = ub%lo; uhi(ix) = ub%hi
        end do
        do ix = 1, nx
            call sample_streaming_moments_box(moments, interval(x(ix), x(ix + 1)), &
                ncells(ix), jcells(ix), ucells(ix), ok)
            if (.not. ok) error stop 'whole-cell moment bands unavailable'
        end do
        ! Each vertex covers both adjacent complete cells. Every convex
        ! combination of its lower/upper endpoints therefore covers that cell.
        do ix = 1, nx + 1
            nb = interval(nlo(ix), nhi(ix))
            jb = interval(jlo(ix), jhi(ix))
            ub = interval(ulo(ix), uhi(ix))
            left = max(1, ix - 1)
            right = min(nx, ix)
            nb = hull(nb, hull(ncells(left), ncells(right)))
            jb = hull(jb, hull(jcells(left), jcells(right)))
            ub = hull(ub, hull(ucells(left), ucells(right)))
            nlo(ix) = nb%lo; nhi(ix) = nb%hi
            jlo(ix) = jb%lo; jhi(ix) = jb%hi
            ulo(ix) = ub%lo; uhi(ix) = ub%hi
        end do
        call cpu_time(finish)
        post_s = finish - start
        if (.not. recorded(number)) then
            factor = -1.0_dp
            if (solve_s > 0.0_dp) factor = complete_s/solve_s
            write (evidence, '(i0,6(",",es24.16e3))') number, t, solve_s, complete_s, &
                complete_s - solve_s, factor, post_s
            do ix = 1, nx + 1
                write (bands, '(i0,11(",",es24.16e3))') number, t, x(ix), &
                    density(ix), nlo(ix), nhi(ix), current(ix), jlo(ix), jhi(ix), &
                    bulk(ix), ulo(ix), uhi(ix)
                do iv = 1, nv + 1
                    write (distribution_data, '(i0,6(",",es24.16e3))') &
                        number, t, x(ix), v(iv), image(iv, ix), &
                        image_lo(iv, ix), image_hi(iv, ix)
                end do
            end do
            recorded(number) = .true.
        end if
    end subroutine

    subroutine draw_frame()
        character(len=96) :: heading
        call subplot(2, 2, 1)
        write (heading, '(a,f4.1,a)') 'Weighted particles | T [1] = ', t, &
            ' | area ~ weight'
        call title(trim(heading))
        call scatter(px, pv, s=glyph, color=[0.0_dp, 0.4471_dp, 0.6980_dp], &
            alpha=0.65_dp)
        call xlabel('2*pi*x/L [1]'); call ylabel('v/v_ref [1]')
        call xlim(0.0_dp, period); call ylim(-1.0_dp, 1.0_dp)
        call subplot(2, 2, 2)
        call title('2*v_ref*f/n0 [1] | candidate; color [0.25,2.25]')
        call pcolormesh(x, v, image, cmap='viridis', vmin=0.25_dp, vmax=2.25_dp)
        call xlabel('2*pi*x/L [1]'); call ylabel('v/v_ref [1]')
        call xlim(0.0_dp, period); call ylim(-1.0_dp, 1.0_dp)
        call subplot(2, 2, 3)
        call title('Density n = (1/2) integral f dv')
        call fill_between(x, nlo, nhi, color='lightblue', alpha=0.5_dp)
        call plot(x, density, color=[0.0_dp, 0.4471_dp, 0.6980_dp], &
            label='n + whole-cell bound')
        call xlabel('2*pi*x/L [1]'); call ylabel('n/n0 [1]')
        call xlim(0.0_dp, period); call ylim(0.35_dp, 1.65_dp)
        call legend(loc='lower right')
        call subplot(2, 2, 4)
        call title('Bulk velocity u = j/n | moments omit shape')
        call fill_between(x, ulo, uhi, color='moccasin', alpha=0.5_dp)
        call plot(x, bulk, color=[0.8353_dp, 0.3686_dp, 0.0_dp], linestyle='--', &
            label='u + whole-cell bound')
        call xlabel('2*pi*x/L [1]')
        call ylabel('u/v_ref [1]')
        call xlim(0.0_dp, period); call ylim(-0.15_dp, 0.40_dp)
        call legend(loc='lower right')
        call tight_layout()
    end subroutine

    subroutine update(index)
        integer, intent(in) :: index
        call compute_frame(index)
        call canvas%clear()
        call draw_frame()
        call canvas%set_rendered(.false.)
    end subroutine
end program streaming_moments_demo
