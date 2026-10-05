program debye_sphere_demo
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_interval, only: interval_t, interval, rad
    use kin6d_debye_sphere_reference, only: debye_sphere_values, debye_sphere_bands
    use kin6d_debye_sphere_leaves, only: sphere_inside, sphere_outside, &
        sphere_bare_inside, sphere_bare_outside, sphere_coordinate, sphere_parameter, &
        sphere_response, sphere_response_interval
    use fortplot, only: figure, get_global_figure, figure_t, subplot, scatter, &
        plot, fill_between, title, xlabel, ylabel, xlim, ylim, legend, &
        tight_layout, savefig
    use fortplot_animation, only: animation_t, animation_func => FuncAnimation, &
        save_animation
    use fortplot_system_runtime, only: create_directory_runtime
    implicit none
    real(dp), parameter :: pi = 3.14159265358979323846_dp
    integer, parameter :: cells = 192, points = cells + 1
    integer, parameter :: frames = 121, fps = 24, repeats = 3
    integer, parameter :: circle_points = 129, max_markers = 121
    real(dp), parameter :: chi = 0.02_dp
    real(dp) :: r(points), phi(points), bare(points), tangent(points)
    real(dp) :: lower(3, points), upper(3, points), work(3, cells)
    real(dp) :: ne(points), ne_lo(points), ne_hi(points)
    real(dp) :: circle_x(circle_points), circle_y(circle_points)
    real(dp) :: source_x(max_markers), source_y(max_markers)
    real(dp) :: electron_x(max_markers), electron_y(max_markers)
    real(dp) :: k, angle, marker_x, marker_y, marker_r, start, finish
    real(dp) :: primal_cpu, derivative_cpu, complete_cpu, max_phi_radius, max_g_radius
    type(interval_t) :: ne_box
    type(figure_t), pointer :: canvas
    type(animation_t) :: animation
    integer :: i, j, frame, repeat, snapshot_frame, status, data, evidence
    integer :: source_count, electron_count
    logical :: ok, snapshot_only, response_view
    character(len=64) :: argument

    snapshot_only = .false.; snapshot_frame = 61; response_view = .false.
    if (command_argument_count() > 2) error stop 'usage: demo [--snapshot [frame]]'
    if (command_argument_count() >= 1) then
        call get_command_argument(1, argument)
        if (trim(argument) /= '--snapshot') error stop 'expected --snapshot'
        snapshot_only = .true.
    end if
    if (command_argument_count() == 2) then
        call get_command_argument(2, argument)
        read (argument, *, iostat=status) snapshot_frame
        if (status /= 0) error stop 'snapshot frame must be integer'
    end if
    if (snapshot_frame < 1 .or. snapshot_frame > frames) &
        error stop 'frame outside 1..121'
    call create_directory_runtime('output', ok)
    if (.not. ok) error stop 'cannot create output directory'
    do i = 1, points
        call sphere_coordinate(real(i - 1, dp), real(cells, dp), r(i))
    end do
    do i = 1, circle_points
        angle = 2*pi*real(i - 1, dp)/real(circle_points - 1, dp)
        circle_x(i) = cos(angle); circle_y(i) = sin(angle)
    end do
    call prepare_illustration()
    open (newunit=data, file='output/debye-potential-reduction.csv', status='replace')
    write (data, '(a)') 'frame,kappa,r_over_R,screened,'// &
        'screened_lower,screened_upper,'// &
        'bare,bare_lower,bare_upper,dphi_dk,dphi_dk_lower,dphi_dk_upper,'// &
        'electron_density,electron_lower,electron_upper'
    open (newunit=evidence, file='output/debye-sphere-evidence.csv', status='replace')
    write (evidence, '(a)') 'frame,kappa,primal_s,primal_tangent_s,'// &
        'complete_bands_s,'// &
        'extra_over_primal_s,factor_over_primal,extra_over_tangent_s,'// &
        'factor_over_tangent,max_potential_band_radius,max_tangent_band_radius'
    if (snapshot_only) then
        call compute_frame(snapshot_frame, timed=.true.)
        call write_evidence(snapshot_frame)
    else
        do frame = 1, frames
            call compute_frame(frame, timed=.true.)
            call write_evidence(frame)
        end do
    end if
    close (data); close (evidence)
    call compute_frame(snapshot_frame, timed=.false.)
    call figure(figsize=[11.0_dp, 5.0_dp])
    canvas => get_global_figure()
    call draw_frame()
    call savefig('output/debye-potential-reduction.png')
    response_view = .true.
    call figure(figsize=[11.0_dp, 5.0_dp])
    canvas => get_global_figure()
    call draw_frame()
    call savefig('output/debye-sphere-response.png')
    if (.not. snapshot_only) then
        response_view = .false.
        call figure(figsize=[11.0_dp, 5.0_dp])
        canvas => get_global_figure()
        animation = animation_func(update_frame, frames=frames, interval=42, fig=canvas)
        call save_animation(animation, 'output/debye-potential-reduction.mp4', fps=fps)
        response_view = .true.
        call figure(figsize=[11.0_dp, 5.0_dp])
        canvas => get_global_figure()
        animation = animation_func(update_frame, frames=frames, interval=42, fig=canvas)
        call save_animation(animation, 'output/debye-sphere-response.mp4', fps=fps)
    end if
    write (*, '(a)') 'ANALYTICAL REFERENCE: matched permeable uniform 3D sphere.'
    write (*, '(a)') 'Equilibrium kappa sweep; no transient or N-body screening.'
    write (*, '(a)') 'R=1; Phi unit Q/(4*pi*epsilon0*R); r/R [1]; kappa=R/lambdaD [1].'
    write (*, '(a)') 'Electron n/n0=1+chi*Phi, chi=.02; max bare ePhi/kBT=.03.'
    write (*, '(a)') 'Strict analytical evaluation and full radial chord bands.'
    write (*, '(a)') 'Model error unavailable; dPhi/dk interval certified.'
    write (*, '(a)') 'R/Q/input/model sensitivities not reported.'
    write (*, '(a)') 'CPU timings: cold frame bands and repeated references.'
    write (*, '(a)') 'Rendering and encoding excluded from scientific timers.'
contains
    !> Illustration marker grid on [-2, 2]: -2 + 4 index/(count - 1).
    pure real(dp) function marker(index, count)
        integer, intent(in) :: index, count
        marker = -2 + 4*real(index, dp)/real(count - 1, dp)
    end function

    subroutine prepare_illustration()
        source_count = 0; electron_count = 0
        do j = 0, 10
            marker_y = marker(j, 11)
            do i = 0, 10
                marker_x = marker(i, 11)
                marker_r = sqrt(marker_x**2 + marker_y**2)
                if (marker_r < 0.95_dp) then
                    source_count = source_count + 1
                    source_x(source_count) = marker_x
                    source_y(source_count) = marker_y
                end if
            end do
        end do
        do j = 0, 8
            marker_y = marker(j, 9)
            do i = 0, 8
                marker_x = marker(i, 9)
                electron_count = electron_count + 1
                electron_x(electron_count) = marker_x
                electron_y(electron_count) = marker_y
            end do
        end do
    end subroutine prepare_illustration

    subroutine compute_frame(number, timed)
        integer, intent(in) :: number
        logical, intent(in) :: timed
        call sphere_parameter(real(number - 1, dp), real(frames, dp), k)
        if (timed) then
            call cpu_time(start)
            do repeat = 1, repeats
                do i = 1, points
                    if (r(i) <= 1.0_dp) then
                        call sphere_inside(r(i), k, phi(i))
                        call sphere_bare_inside(r(i), bare(i))
                    else
                        call sphere_outside(r(i), k, phi(i))
                        call sphere_bare_outside(r(i), bare(i))
                    end if
                end do
            end do
            call cpu_time(finish)
            primal_cpu = (finish - start)/repeats
            call cpu_time(start)
            do repeat = 1, repeats
                do i = 1, points
                    call debye_sphere_values(r(i), k, phi(i), bare(i), tangent(i), ok)
                    if (.not. ok) error stop 'reference derivative evaluation failed'
                end do
            end do
            call cpu_time(finish)
            derivative_cpu = (finish - start)/repeats
        end if
        call cpu_time(start)
        call debye_sphere_bands(r, k, phi, bare, tangent, lower, upper, work, ok)
        call cpu_time(finish)
        if (.not. ok) error stop 'reference whole-cell bands failed'
        complete_cpu = finish - start
        max_phi_radius = 0.0_dp; max_g_radius = 0.0_dp
        do i = 1, points
            max_phi_radius = max(max_phi_radius, &
                rad(interval(lower(1, i), upper(1, i))))
            max_g_radius = max(max_g_radius, rad(interval(lower(3, i), upper(3, i))))
            call sphere_response(phi(i), chi, ne(i))
            call sphere_response_interval(interval(lower(1, i), upper(1, i)), &
                interval(chi), ne_box)
            ne_lo(i) = ne_box%lo; ne_hi(i) = ne_box%hi
        end do
    end subroutine compute_frame

    subroutine write_evidence(number)
        integer, intent(in) :: number
        do i = 1, points
            write (data, '(i0,14(",",es24.16e3))') number, k, r(i), phi(i), &
                lower(1, i), upper(1, i), bare(i), lower(2, i), upper(2, i), &
                tangent(i), lower(3, i), upper(3, i), ne(i), ne_lo(i), ne_hi(i)
        end do
        write (evidence, '(i0,10(",",es24.16e3))') number, k, primal_cpu, &
            derivative_cpu, complete_cpu, complete_cpu - primal_cpu, &
            timer_factor(complete_cpu, primal_cpu), complete_cpu - derivative_cpu, &
            timer_factor(complete_cpu, derivative_cpu), max_phi_radius, max_g_radius
    end subroutine write_evidence

    real(dp) function timer_factor(total, baseline) result(value)
        real(dp), intent(in) :: total, baseline
        value = -1.0_dp
        if (baseline > 0.0_dp) value = total/baseline
    end function timer_factor

    subroutine update_frame(number)
        integer, intent(in) :: number
        call compute_frame(number, timed=.false.)
        call canvas%clear()
        call draw_frame()
        call canvas%set_rendered(.false.)
    end subroutine update_frame

    subroutine draw_frame()
        character(len=120) :: heading
        call subplot(1, 2, 1)
        call title('Permeable 3D sphere | cross-section illustration')
        call xlabel('x/R [1]'); call ylabel('y/R [1]')
        call xlim(-2.25_dp, 2.25_dp); call ylim(-2.25_dp, 2.25_dp)
        call scatter(electron_x(:electron_count), electron_y(:electron_count), &
            s=8.0_dp, color=[0.4_dp, 0.4_dp, 0.4_dp], marker='o', &
            label='mobile electrons (illustration)')
        call plot(circle_x, circle_y, color=[0.8353_dp, 0.3686_dp, 0.0_dp])
        call scatter(source_x(:source_count), source_y(:source_count), &
            s=45.0_dp, color=[0.8353_dp, 0.3686_dp, 0.0_dp], marker='+', &
            label='fixed uniform source; R=1')
        call legend(loc='upper left')
        call subplot(1, 2, 2)
        write (heading, '(a,f5.2)') 'Equilibrium reference | R/lambdaD = ', k
        call title(trim(heading))
        call xlabel('r/R [1]'); call xlim(0.0_dp, 6.0_dp)
        if (response_view) then
            call ylabel('electron n/n0 [1]'); call ylim(0.999_dp, 1.035_dp)
            call envelope(ne_lo, ne_hi, 'lightblue', [0.0_dp, 0.4471_dp, 0.6980_dp])
            call plot(r, ne, color=[0.0_dp, 0.4471_dp, 0.6980_dp], &
                label='linear response + full bound')
            call plot(r, [(1.0_dp, i=1, points)], linestyle='--', &
                color=[0.4_dp, 0.4_dp, 0.4_dp], label='background n/n0=1')
        else
            call ylabel('Phi / [Q/(4 pi epsilon0 R)] [1]'); call ylim(-0.02_dp, 1.80_dp)
            call envelope(lower(2, :), upper(2, :), 'moccasin', &
                [0.8353_dp, 0.3686_dp, 0.0_dp])
            call envelope(lower(1, :), upper(1, :), 'lightblue', &
                [0.0_dp, 0.4471_dp, 0.6980_dp])
            call plot(r, bare, linestyle='--', color=[0.8353_dp, 0.3686_dp, 0.0_dp], &
                label='bare finite source + full bound')
            call plot(r, phi, color=[0.0_dp, 0.4471_dp, 0.6980_dp], &
                label='screened equilibrium + full bound')
        end if
        call legend(loc='upper right')
        call tight_layout()
    end subroutine draw_frame

    subroutine envelope(lo, hi, fill, color)
        real(dp), intent(in) :: lo(:), hi(:), color(3)
        character(*), intent(in) :: fill
        call fill_between(r, lo, hi, color=fill, alpha=0.5_dp)
        call plot(r, lo, linestyle=':', color=color)
        call plot(r, hi, linestyle=':', color=color)
    end subroutine envelope
end program debye_sphere_demo
