module kin6d_demo_application
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_is_nan
    use kin6d_laboratory, only: laboratory_t, initialize, advance, sample, &
                                solve_steady, certificate_radius, certificate_accepts
    use kin6d_leaves, only: coordinate, reconstruct, interpolate, &
                            characteristic_reference, interval_reconstruct, &
                            interval_band
    use fortnum_interval, only: interval_t, interval, interval_pi
    use fortplot, only: figure_t, figure, get_global_figure, subplot, title, &
                   xlabel, ylabel, xlim, ylim, pcolormesh, plot, fill_between, legend, &
                        savefig, tight_layout, colorbar, suptitle, set_yscale
    use kin6d_demo_style, only: blue, vermillion, black, grey, movie_fps, hold_frames, num
    use fortplot_animation, only: animation_t, FuncAnimation, save_animation
    use fortplot_system_runtime, only: create_directory_runtime
    use kin6d_result_io, only: run_t, result_file_t
    implicit none
    integer, parameter :: nx = 192
    ! 193 frames of dt = 5/128 span t = 0..7.5 in 8 s at 24 fps.
    integer :: nv = 256, frames = 193
    real(dp) :: time_step = 0.0390625_dp, steady_nu = 0.75_dp
    namelist /kin6d_demo/ nv, frames, time_step, steady_nu
    type(run_t) :: run
    type(result_file_t) :: result
    type(laboratory_t) :: state, coarse, steady
    type(figure_t), pointer :: fig
    type(animation_t) :: animation
    type(interval_t) :: pi, center, lower, upper
    real(dp) :: xedges(nx + 1)
    real(dp), allocatable :: image(:, :), profile(:), exact(:)
    real(dp), allocatable :: lower_line(:), upper_line(:), coarse_line(:)
    real(dp), allocatable :: error_line(:), error_lower(:), &
        error_upper(:)
    real(dp), allocatable :: frame_time(:), frame_bound(:)
    real(dp) :: x, normalized_x, velocity, m, a, b, radius, residual, &
                solve_start, solve_end
    real(dp) :: bound_max, final_time
    integer :: i, status, ledger, step, unit, ios
    logical :: directory_ok, movie
    character(96) :: header
    character(32) :: argument
    character(256) :: msg
contains
    subroutine run_demo()
        pi = interval_pi()
        call run%start('kin6d_demo')
        ios = 0; msg = ''
        call run%open_config(unit)
        if (run%has_config()) read (unit, nml=kin6d_demo, iostat=ios, iomsg=msg)
        call run%close_config(unit, 'kin6d_demo', ios, msg)
        if (nv < 4) error stop 'nv must be at least 4 (coarse grid uses nv/2)'
        if (frames < 1) error stop 'frames must be positive'
        if (.not. (time_step > 0.0_dp .and. steady_nu > 0.0_dp)) &
            error stop 'time_step and steady_nu must be positive'
        call run%open_resolved(unit)
        write (unit, nml=kin6d_demo)
        call run%close_resolved(unit)
        allocate (image(nv, nx), profile(nv + 1), exact(nv + 1), lower_line(nv + 1), &
            upper_line(nv + 1), coarse_line(nv + 1), frame_time(frames), &
            frame_bound(frames), error_line(nv + 1), &
            error_lower(nv + 1), error_upper(nv + 1))
        ! Frames not evaluated (--snapshot renders only the last) stay NaN.
        frame_time = ieee_value(0.0_dp, ieee_quiet_nan); frame_bound = frame_time
        movie = .true.
        if (run%argument_count() > 0) then
            call run%argument(1, argument)
            if (trim(argument) == '--snapshot') movie = .false.
        end if
        call create_directory_runtime('output', directory_ok)
        if (.not. directory_ok) error stop 'cannot create output directory'
        call initialize(state, nv)
        call initialize(coarse, nv/2)
        call initialize(steady, nv)
        call cpu_time(solve_start)
        call solve_steady(steady, steady_nu, residual)
        call cpu_time(solve_end)
        print '(a,es12.4,a,f10.6)', 'Direct steady residual=', residual, &
            ' CPU seconds=', solve_end - solve_start
        do i = 1, nx + 1
            call coordinate(0.0_dp, 1.0_dp, real(i - 1, dp), real(nx, dp), xedges(i))
        end do
        open (newunit=ledger, file='output/streaming.csv', status='replace')
        write (ledger, '(a)') 'frame,time,continuous_sup_error_bound'
        call figure(figsize=[13.2_dp, 8.8_dp])
        fig => get_global_figure()
        bound_max = 0.0_dp
        final_time = 0.0_dp
        if (movie) then
            animation = FuncAnimation(update_frame, frames=frames + hold_frames, &
                                      interval=50, fig=fig)
            call save_animation(animation, 'output/free-streaming.mp4', movie_fps, status)
            if (status /= 0) error stop 'fortplot animation encoding failed'
        else
            do step = 1, frames - 1
                call advance(state, time_step)
                call advance(coarse, time_step)
            end do
            call render_frame(frames)
        end if
        call savefig('output/final-snapshot.png')
        close (ledger)
        call write_result()
        print '(a,f7.3,a,es12.4)', 'Final time=', final_time, &
            ' uniform snapshot bound=', radius
       print '(a)', 'Bounds use reviewed interpolation analysis and FortNum arithmetic.'
        print '(a)', 'This streaming bound has no Lean or velocity-tail claim.'
    end subroutine run_demo

    subroutine update_frame(frame)
        integer, intent(in) :: frame
        ! Fortplot requests frames sequentially, starting at one. Frames past
        ! the last physical frame hold the final state for about one second.
        if (frame > frames) then
            call draw_frame(frames)
            return
        end if
        if (frame > 1) then
            call advance(state, time_step)
            call advance(coarse, time_step)
        end if
        call render_frame(frame)
    end subroutine

    subroutine render_frame(frame)
        integer, intent(in) :: frame
        integer :: ix, iv
        radius = certificate_radius(state)
        if (.not. certificate_accepts(state, &
                                      radius)) error stop 'snapshot bound rejected'
        bound_max = max(bound_max, radius)
        final_time = state%time
        write (ledger, '(i0,a,f12.6,a,es24.16)') frame, ',', state%time, ',', radius
        frame_time(frame) = state%time; frame_bound(frame) = radius
        ! One rendered tile: no stored trajectory and no duplicate phase-space state.
        do ix = 1, nx
            call coordinate(xedges(ix), xedges(ix + 1), 1.0_dp, 2.0_dp, normalized_x)
            call coordinate(0.0_dp, 2.0_dp*pi%lo, normalized_x, 1.0_dp, x)
            do iv = 1, nv
                call interpolate(state%m(iv), state%m(iv + 1), 0.5_dp, m)
                call interpolate(state%a(iv), state%a(iv + 1), 0.5_dp, a)
                call interpolate(state%b(iv), state%b(iv + 1), 0.5_dp, b)
                call reconstruct(m, a, b, x, image(iv, ix))
            end do
        end do
        do iv = 1, nv + 1
            velocity = state%v(iv)
            call reconstruct(state%m(iv), state%a(iv), state%b(iv), 0.0_dp, profile(iv))
            call characteristic_reference(0.0_dp, velocity, state%time, &
                                          state%amplitude, exact(iv))
            call sample(coarse, 0.0_dp, velocity, coarse_line(iv))
            call interval_reconstruct(interval(state%m(iv)), interval(state%a(iv)), &
                                      interval(state%b(iv)), interval(0.0_dp), center)
            call interval_band(center, interval(radius), lower, upper)
            lower_line(iv) = lower%lo
            upper_line(iv) = upper%hi
        end do
        error_line = profile - exact
        error_lower = lower_line - exact
        error_upper = upper_line - exact
        call draw_frame(frame)
    end subroutine

    subroutine draw_frame(frame)
        integer, intent(in) :: frame
        call figure(figsize=[13.2_dp, 8.8_dp])
        call suptitle('Free streaming shears f(x,v) into ever finer velocity '// &
            'stripes; the certified error bound grows with them (E = 0)')
        call subplot(2, 2, 1)
        call title('Phase-space density  |  t = '//num(state%time, 2))
        call xlabel('x / L  (periodic)')
        call ylabel('$v/v_{th}$')
        call xlim(0.0_dp, 1.0_dp)
        call ylim(-4.0_dp, 4.0_dp)
        call pcolormesh(xedges, state%v, image, cmap='viridis', vmin=0.0_dp, &
                        vmax=1.5_dp)
        call colorbar(label='$f/F_M(0)$')
        call subplot(2, 2, 2)
        call title('Velocity slice at x = 0')
        call xlabel('$v/v_{th}$')
        call ylabel('$f/F_M(0)$')
        call xlim(-4.0_dp, 4.0_dp)
        ! Retain the raw enclosure below zero in the near-zero velocity tails.
        call ylim(-0.05_dp, 1.65_dp)
        call fill_between(state%v, lower_line, upper_line, color=blue, alpha=0.25_dp)
        call plot(state%v, profile, label=itoa(nv)//' cells', color=blue, &
                  linewidth=2.0_dp)
        call plot(state%v, coarse_line, label=itoa(nv/2)//' cells', &
                  linestyle=':', color=vermillion, linewidth=2.0_dp)
        call plot(state%v, exact, label='exact characteristics', linestyle='--', &
                  color=black, linewidth=1.2_dp)
        call legend(loc='upper right')
        call subplot(2, 2, 3)
        call title('Error of the 256-cell slice at x = 0 (band invisible above)')
        call xlabel('$v/v_{th}$')
        call ylabel('$(f - f_{exact})/F_M(0)$')
        call xlim(-4.0_dp, 4.0_dp)
        call ylim(-5.0e-3_dp, 5.0e-3_dp)
        call fill_between(state%v, error_lower, error_upper, color=blue, alpha=0.25_dp)
        call plot(state%v, error_line, label='256 cells ± certified bound', color=blue, &
                  linewidth=2.0_dp)
        call legend(loc='lower right')
        call subplot(2, 2, 4)
        call title('Certified sup-error bound over the whole (x,v) domain')
        call xlabel('t $v_{th}$ / L')
        call ylabel('bound / $F_M(0)$')
        call xlim(0.0_dp, real(frames - 1, dp)*time_step)
        call ylim(1.0e-4_dp, 1.0e-2_dp)
        call set_yscale('log')
        ! --snapshot evaluates only the last frame; unevaluated frames are NaN.
        call plot(pack(frame_time(1:frame), .not. ieee_is_nan(frame_bound(1:frame))), &
                  pack(frame_bound(1:frame), .not. ieee_is_nan(frame_bound(1:frame))), &
                  color=blue, linewidth=2.0_dp)
        call plot([frame_time(frame)], [frame_bound(frame)], color=blue, &
                  linestyle='none', marker='o', markersize=6.0_dp)
        call tight_layout()
    end subroutine

    pure function itoa(i) result(text)
        integer, intent(in) :: i
        character(:), allocatable :: text
        character(16) :: buffer
        write (buffer, '(i0)') i
        text = trim(buffer)
    end function

    subroutine write_result()
        character(len=*), parameter :: bound = 'certified-continuous pointwise '// &
            'sup-error enclosure (reviewed interpolation analysis); no velocity-tail claim'
        call result%create(run)
        call result%write('steady_residual', residual, '1', 'direct steady residual', &
            status='candidate')
        call result%write('time', frame_time, 'L/v_th', 'frame time', dims='frame')
        call result%write('sup_error_bound', frame_bound, 'f/F_M(0)', &
            'continuous sup error bound of the streaming snapshot (NaN: frame '// &
            'not evaluated)', dims='frame', &
            status=bound)
        call result%write('velocity', state%v, 'v_th', 'velocity nodes', dims='v')
        call result%write('final_profile', profile, 'f/F_M(0)', &
            'final slice x=0, fine grid', dims='v', status='candidate center')
        call result%write('final_lower', lower_line, 'f/F_M(0)', &
            'final slice lower band', dims='v', status=bound)
        call result%write('final_upper', upper_line, 'f/F_M(0)', &
            'final slice upper band', dims='v', status=bound)
        call result%write('final_exact', exact, 'f/F_M(0)', &
            'exact characteristic reference at x=0', dims='v', status='exact oracle')
        call result%write('final_coarse', coarse_line, 'f/F_M(0)', &
            'coarse-grid slice x=0', dims='v', status='candidate')
        call result%close()
    end subroutine
end module

program kin6d_demo
    use kin6d_demo_application, only: run_demo
    implicit none
    call run_demo()
end program
