module debye_sphere_demo_callbacks
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_interval, only: interval_t, interval, rad
    use kin6d_debye_sphere_reference, only: debye_sphere_values, debye_sphere_bands
    use kin6d_debye_sphere_leaves, only: sphere_inside, sphere_outside, &
        sphere_bare_inside, sphere_bare_outside, sphere_coordinate, sphere_parameter, &
        sphere_response, sphere_response_interval
    use fortplot, only: figure, get_global_figure, figure_t, subplot, &
        plot, fill_between, title, xlabel, ylabel, xlim, ylim, legend, &
        tight_layout, savefig, pcolormesh, colorbar, suptitle, set_yscale
    use kin6d_demo_style, only: blue, vermillion, green, black, grey, hold_frames, num
    use fortplot_animation, only: animation_t, animation_func => FuncAnimation, &
        save_animation
    use fortplot_system_runtime, only: create_directory_runtime
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use kin6d_result_io, only: run_t, result_file_t

    use fortplot_animation, only: animation_callback_t
    implicit none
    private
    public :: run_debye_sphere_demo
    real(dp), parameter :: pi = 3.14159265358979323846_dp
    integer, parameter :: fps = 24
    integer, parameter :: circle_points = 129, map_cells = 120
    real(dp), parameter :: map_x = 3.23_dp  ! equal aspect in the measured axes box
    character(len=*), parameter :: bound = 'certified-continuous strict analytical '// &
        'reference band; model error unavailable'
    type :: run_debye_sphere_demo_context_t
        real(dp), pointer :: bare(:) => null()
        type(figure_t), pointer :: canvas => null()
        real(dp), pointer :: chi => null()
        real(dp), pointer :: circle_x(:) => null()
        real(dp), pointer :: circle_y(:) => null()
        real(dp), pointer :: complete_cpu => null()
        real(dp), pointer :: derivative_cpu => null()
        real(dp), pointer :: finish => null()
        integer, pointer :: frames => null()
        integer, pointer :: i => null()
        real(dp), pointer :: k => null()
        real(dp), pointer :: lower(:,:) => null()
        real(dp), pointer :: map(:,:) => null()
        real(dp), pointer :: map_xe(:) => null()
        real(dp), pointer :: map_ye(:) => null()
        real(dp), pointer :: max_g_radius => null()
        real(dp), pointer :: max_phi_radius => null()
        real(dp), pointer :: ne(:) => null()
        type(interval_t), pointer :: ne_box => null()
        real(dp), pointer :: ne_hi(:) => null()
        real(dp), pointer :: ne_lo(:) => null()
        logical, pointer :: ok => null()
        real(dp), pointer :: phi(:) => null()
        integer, pointer :: points => null()
        real(dp), pointer :: primal_cpu => null()
        real(dp), pointer :: r(:) => null()
        integer, pointer :: repeat => null()
        integer, pointer :: repeats => null()
        real(dp), pointer :: start => null()
        real(dp), pointer :: tangent(:) => null()
        real(dp), pointer :: upper(:,:) => null()
        real(dp), pointer :: work(:,:) => null()
    end type

    type, extends(animation_callback_t) :: run_debye_sphere_demo_update_frame_t
        type(run_debye_sphere_demo_context_t), pointer :: data => null()
    contains
        procedure :: update => run_debye_sphere_demo_update_frame_binding
    end type

contains
    subroutine run_debye_sphere_demo()
    implicit none
    ! Cross-section map: map_cells^2 display cells on [-map_x, map_x] x [-2.5, 2.5].
integer , target :: cells, points
    ! 193 equilibria (kappa = 0.5..2) play in 8 s at 24 fps, then a 1 s hold.
integer , target :: frames, repeats
real(dp) , target :: chi
    namelist /debye_sphere/ cells, frames, repeats, chi
    type(run_t) :: run
    type(result_file_t) :: result
    integer :: unit, ios
    character(len=256) :: msg
    real(dp), allocatable , target :: r(:), phi(:), bare(:), tangent(:)
    real(dp), allocatable , target :: lower(:, :), upper(:, :), work(:, :)
    real(dp), allocatable , target :: ne(:), ne_lo(:), ne_hi(:)
    ! Result record (point, frame, quantity); NaN where a frame was not computed.
    real(dp), allocatable :: kappas(:), record(:, :, :)
    real(dp) , target :: circle_x(circle_points), circle_y(circle_points)
    real(dp) , target :: map_xe(map_cells + 1), map_ye(map_cells + 1), &
        map(map_cells, map_cells)
    real(dp) , target :: k, angle, start, finish
    real(dp) , target :: primal_cpu, derivative_cpu, complete_cpu, max_phi_radius, &
        max_g_radius
    type(interval_t) , target :: ne_box
    type(figure_t), pointer :: canvas
    type(animation_t) :: animation
    integer , target :: i, j, frame, repeat, snapshot_frame, status, data, evidence
    logical , target :: ok, snapshot_only
    character(len=64) :: argument
        type(run_debye_sphere_demo_context_t), target :: ctx
        type(run_debye_sphere_demo_update_frame_t), target :: update_frame_operator

        cells = 192
        frames = 193
        repeats = 3
        chi = 0.02_dp
    call run%start('debye_sphere_demo')
    ios = 0; msg = ''
    call run%open_config(unit)
    if (run%has_config()) read (unit, nml=debye_sphere, iostat=ios, iomsg=msg)
    call run%close_config(unit, 'debye_sphere', ios, msg)
    if (cells < 1 .or. frames < 2 .or. repeats < 1) &
        error stop 'debye_sphere needs cells >= 1, frames >= 2, repeats >= 1'
    if (.not. (chi >= 0.0_dp)) error stop 'debye_sphere chi must be nonnegative'
    call run%open_resolved(unit)
    write (unit, nml=debye_sphere)
    call run%close_resolved(unit)
    points = cells + 1
    allocate (r(points), phi(points), bare(points), tangent(points))
    allocate (lower(3, points), upper(3, points), work(3, cells))
    allocate (ne(points), ne_lo(points), ne_hi(points))
    allocate (kappas(frames), record(points, frames, 12))
    kappas = ieee_value(1.0_dp, ieee_quiet_nan); record = kappas(1)
    snapshot_only = .false.; snapshot_frame = (frames + 1)/2
    if (run%argument_count() > 2) error stop 'usage: demo [--config FILE] [--snapshot [frame]]'
    if (run%argument_count() >= 1) then
        call run%argument(1, argument)
        if (trim(argument) /= '--snapshot') error stop 'expected --snapshot'
        snapshot_only = .true.
    end if
    if (run%argument_count() == 2) then
        call run%argument(2, argument)
        read (argument, *, iostat=status) snapshot_frame
        if (status /= 0) error stop 'snapshot frame must be integer'
    end if
    if (snapshot_frame < 1 .or. snapshot_frame > frames) &
        error stop 'frame outside 1..frames'
    call create_directory_runtime('output', ok)
    if (.not. ok) error stop 'cannot create output directory'
    do i = 1, points
        call sphere_coordinate(real(i - 1, dp), real(cells, dp), r(i))
    end do
    do i = 1, circle_points
        angle = 2*pi*real(i - 1, dp)/real(circle_points - 1, dp)
        circle_x(i) = cos(angle); circle_y(i) = sin(angle)
    end do
    do i = 1, map_cells + 1
        map_xe(i) = -map_x + 2*map_x*real(i - 1, dp)/real(map_cells, dp)
        map_ye(i) = -2.5_dp + 5.0_dp*real(i - 1, dp)/real(map_cells, dp)
    end do
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
        call bind_run_debye_sphere_demo_context()
        call run_debye_sphere_demo_compute_frame(ctx, snapshot_frame, timed=.true.)
        call write_evidence(snapshot_frame)
    else
        do frame = 1, frames
        call bind_run_debye_sphere_demo_context()
            call run_debye_sphere_demo_compute_frame(ctx, frame, timed=.true.)
            call write_evidence(frame)
        end do
    end if
    close (data); close (evidence)
    call write_result()
        call bind_run_debye_sphere_demo_context()
    call run_debye_sphere_demo_compute_frame(ctx, snapshot_frame, timed=.false.)
    call figure(figsize=[13.2_dp, 8.8_dp])
    canvas => get_global_figure()
        call bind_run_debye_sphere_demo_context()
    call run_debye_sphere_demo_draw_frame(ctx)
    call savefig('output/debye-potential-reduction.png')
    if (.not. snapshot_only) then
        call bind_run_debye_sphere_demo_context()
        animation = animation_func(update_frame_operator, frames=frames + hold_frames, &
            interval=42, fig=canvas)
        call save_animation(animation, 'output/debye-potential-reduction.mp4', fps=fps)
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
    subroutine bind_run_debye_sphere_demo_context()
        nullify(ctx%bare)
        if (allocated(bare)) ctx%bare => bare
        ctx%canvas => canvas
        ctx%chi => chi
        ctx%circle_x => circle_x
        ctx%circle_y => circle_y
        ctx%complete_cpu => complete_cpu
        ctx%derivative_cpu => derivative_cpu
        ctx%finish => finish
        ctx%frames => frames
        ctx%i => i
        ctx%k => k
        nullify(ctx%lower)
        if (allocated(lower)) ctx%lower => lower
        ctx%map => map
        ctx%map_xe => map_xe
        ctx%map_ye => map_ye
        ctx%max_g_radius => max_g_radius
        ctx%max_phi_radius => max_phi_radius
        nullify(ctx%ne)
        if (allocated(ne)) ctx%ne => ne
        ctx%ne_box => ne_box
        nullify(ctx%ne_hi)
        if (allocated(ne_hi)) ctx%ne_hi => ne_hi
        nullify(ctx%ne_lo)
        if (allocated(ne_lo)) ctx%ne_lo => ne_lo
        ctx%ok => ok
        nullify(ctx%phi)
        if (allocated(phi)) ctx%phi => phi
        ctx%points => points
        ctx%primal_cpu => primal_cpu
        nullify(ctx%r)
        if (allocated(r)) ctx%r => r
        ctx%repeat => repeat
        ctx%repeats => repeats
        ctx%start => start
        nullify(ctx%tangent)
        if (allocated(tangent)) ctx%tangent => tangent
        nullify(ctx%upper)
        if (allocated(upper)) ctx%upper => upper
        nullify(ctx%work)
        if (allocated(work)) ctx%work => work
        update_frame_operator%data => ctx
    end subroutine bind_run_debye_sphere_demo_context

    subroutine write_result()
        character(len=16), parameter :: names(12) = [character(len=16) :: &
            'screened', 'screened_lower', 'screened_upper', 'bare', 'bare_lower', &
            'bare_upper', 'dphi_dk', 'dphi_dk_lower', 'dphi_dk_upper', &
            'electron_density', 'electron_lower', 'electron_upper']
        character(len=24), parameter :: units(4) = [character(len=24) :: &
            'Q/(4 pi eps0 R)', 'Q/(4 pi eps0 R)', 'Q/(4 pi eps0 R)', 'n0']
        integer :: q
        call result%create(run)
        call result%attribute('frame_record', 'frames not computed in this run are NaN')
        call result%write('kappa', kappas, '1', 'kappa = R/lambda_D', dims='frame')
        call result%write('r', r, 'R', 'radius r/R', dims='r')
        do q = 1, 12
            if (mod(q, 3) == 1) then
                call result%write(trim(names(q)), record(:, :, q), &
                    trim(units((q + 2)/3)), trim(names(q))//' center', dims='r,frame', &
                    status='candidate center')
            else
                call result%write(trim(names(q)), record(:, :, q), &
                    trim(units((q + 2)/3)), trim(names(q))//' bound', dims='r,frame', &
                    status=bound)
            end if
        end do
        call result%close()
    end subroutine write_result


    subroutine write_evidence(number)
        integer, intent(in) :: number
        kappas(number) = k
        record(:, number, :) = reshape([phi, lower(1, :), upper(1, :), bare, &
            lower(2, :), upper(2, :), tangent, lower(3, :), upper(3, :), ne, ne_lo, &
            ne_hi], [points, 12])
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

    end subroutine run_debye_sphere_demo

    subroutine run_debye_sphere_demo_compute_frame(ctx, number, timed)
        type(run_debye_sphere_demo_context_t), intent(in) :: ctx
        integer, intent(in) :: number
        logical, intent(in) :: timed
        integer, pointer :: local_i
        integer, pointer :: local_repeat
        local_i => ctx%i
        local_repeat => ctx%repeat
        call sphere_parameter(real(number - 1, dp), real(ctx%frames, dp), ctx%k)
        if (timed) then
            call cpu_time(ctx%start)
            do local_repeat = 1, ctx%repeats
                do local_i = 1, ctx%points
                    if (ctx%r(local_i) <= 1.0_dp) then
                        call sphere_inside(ctx%r(local_i), ctx%k, ctx%phi(local_i))
                        call sphere_bare_inside(ctx%r(local_i), ctx%bare(local_i))
                    else
                        call sphere_outside(ctx%r(local_i), ctx%k, ctx%phi(local_i))
                        call sphere_bare_outside(ctx%r(local_i), ctx%bare(local_i))
                    end if
                end do
            end do
            call cpu_time(ctx%finish)
            ctx%primal_cpu = (ctx%finish - ctx%start)/ctx%repeats
            call cpu_time(ctx%start)
            do local_repeat = 1, ctx%repeats
                do local_i = 1, ctx%points
                    call debye_sphere_values(ctx%r(local_i), ctx%k, &
                        ctx%phi(local_i), ctx%bare(local_i), ctx%tangent(local_i), &
                        ctx%ok)
                    if (.not. ctx%ok) error stop &
                        'reference derivative evaluation failed'
                end do
            end do
            call cpu_time(ctx%finish)
            ctx%derivative_cpu = (ctx%finish - ctx%start)/ctx%repeats
        end if
        call cpu_time(ctx%start)
        call debye_sphere_bands(ctx%r, ctx%k, ctx%phi, ctx%bare, ctx%tangent, &
            ctx%lower, ctx%upper, ctx%work, ctx%ok)
        call cpu_time(ctx%finish)
        if (.not. ctx%ok) error stop 'reference whole-cell bands failed'
        ctx%complete_cpu = ctx%finish - ctx%start
        ctx%max_phi_radius = 0.0_dp; ctx%max_g_radius = 0.0_dp
        do local_i = 1, ctx%points
            ctx%max_phi_radius = max(ctx%max_phi_radius, &
                rad(interval(ctx%lower(1, local_i), ctx%upper(1, local_i))))
            ctx%max_g_radius = max(ctx%max_g_radius, rad(interval(ctx%lower(3, &
                local_i), ctx%upper(3, local_i))))
            call sphere_response(ctx%phi(local_i), ctx%chi, ctx%ne(local_i))
            call sphere_response_interval(interval(ctx%lower(1, local_i), &
                ctx%upper(1, local_i)), &
                interval(ctx%chi), ctx%ne_box)
            ctx%ne_lo(local_i) = ctx%ne_box%lo; ctx%ne_hi(local_i) = ctx%ne_box%hi
        end do
    end subroutine run_debye_sphere_demo_compute_frame

    subroutine run_debye_sphere_demo_update_frame(ctx, number)
        type(run_debye_sphere_demo_context_t), intent(in) :: ctx
        integer, intent(in) :: number
        call run_debye_sphere_demo_compute_frame(ctx, min(number, ctx%frames), &
            timed=.false.)
        call ctx%canvas%clear()
        call run_debye_sphere_demo_draw_frame(ctx)
        call ctx%canvas%set_rendered(.false.)
    end subroutine run_debye_sphere_demo_update_frame

    subroutine run_debye_sphere_demo_update_frame_binding(self, frame)
        class(run_debye_sphere_demo_update_frame_t), intent(inout) :: self
        integer, intent(in) :: frame
        call run_debye_sphere_demo_update_frame(self%data, frame)
    end subroutine

    subroutine run_debye_sphere_demo_draw_frame(ctx)
        type(run_debye_sphere_demo_context_t), intent(in) :: ctx
        integer :: ix, iy, cell
        real(dp) :: rho, w, dr
        integer, pointer :: local_i
        local_i => ctx%i
        ! Display only: linear interpolation of the radial reference on the z = 0 plane.
        dr = ctx%r(2) - ctx%r(1)
        do iy = 1, map_cells
            do ix = 1, map_cells
                rho = hypot(0.5_dp*(ctx%map_xe(ix) + ctx%map_xe(ix + 1)), &
                    0.5_dp*(ctx%map_ye(iy) + ctx%map_ye(iy + 1)))
                cell = min(ctx%points - 1, int(rho/dr) + 1)
                w = rho/dr - real(cell - 1, dp)
                ctx%map(iy, ix) = (1 - w)*ctx%phi(cell) + w*ctx%phi(cell + 1)
            end do
        end do
        call suptitle('Linear electron screening shrinks the sphere potential within '// &
            'a Debye length; the bare Coulomb tail disappears')
        call subplot(2, 2, 1)
        call &
            title('Potential at z = 0 (solid: r = R, dashed: r = $\lambda_D$)  |  $R/\lambda_D$ = '//num(ctx%k, 2))
        call xlabel('x / R'); call ylabel('y / R')
        call xlim(-map_x, map_x); call ylim(-2.5_dp, 2.5_dp)
        call pcolormesh(ctx%map_xe, ctx%map_ye, ctx%map, cmap='viridis', &
            vmin=0.0_dp, vmax=1.5_dp)
        call colorbar(label='4πε0 R Φ / Q')
        ! Solid circle: source sphere r = R; dashed circle: r = lambda_D.
        call plot(ctx%circle_x, ctx%circle_y, color=[1.0_dp, 1.0_dp, 1.0_dp], &
            linewidth=2.0_dp)
        call plot(ctx%circle_x/ctx%k, ctx%circle_y/ctx%k, color=[1.0_dp, 1.0_dp, &
            1.0_dp], &
            linestyle='--', linewidth=1.5_dp)
        call subplot(2, 2, 2)
        call title('Radial potential: bare source and screened equilibrium')
        call xlabel('r / R'); call xlim(0.0_dp, 6.0_dp)
        call ylabel('$4\pi\epsilon_0 R \Phi/Q$'); call ylim(-0.02_dp, 1.60_dp)
        call fill_between(ctx%r, ctx%lower(1, :), ctx%upper(1, :), color=blue, &
            alpha=0.25_dp)
        call plot(ctx%r, ctx%phi, color=blue, linewidth=2.0_dp, &
            label='screened ± certified bound')
        call plot(ctx%r, ctx%bare, linestyle='--', color=black, linewidth=1.2_dp, &
            label='bare uniform sphere')
        call plot([1.0_dp, 1.0_dp]/ctx%k, [-0.02_dp, 1.6_dp], color=grey, &
            linestyle=':', &
            linewidth=1.5_dp, label='$r = \lambda_D$')
        call legend(loc='upper right')
        call subplot(2, 2, 3)
        call title('Electron density response $n_e/n_0 = 1 + \chi\Phi$')
        call xlabel('r / R'); call xlim(0.0_dp, 6.0_dp)
        call ylabel('$n_e/n_0$'); call ylim(0.999_dp, 1.035_dp)
        call fill_between(ctx%r, ctx%ne_lo, ctx%ne_hi, color=blue, alpha=0.25_dp)
        call plot(ctx%r, ctx%ne, color=blue, linewidth=2.0_dp, &
            label='linear response ± bound')
        call plot(ctx%r, [(1.0_dp, local_i=1, ctx%points)], linestyle='--', &
            color=black, &
            linewidth=1.2_dp, label='uniform background')
        call legend(loc='upper right')
        call subplot(2, 2, 4)
        call title('Certified band half-widths (too thin to see in the curves)')
        call xlabel('r / R'); call xlim(0.0_dp, 6.0_dp)
        call ylabel('half-width'); call ylim(1.0e-12_dp, 1.0e-2_dp)
        call set_yscale('log')
        call plot(ctx%r, max(0.5_dp*(ctx%upper(1, :) - ctx%lower(1, :)), &
            1.0e-12_dp), color=blue, &
            linewidth=2.0_dp, label='screened potential')
        call plot(ctx%r, max(0.5_dp*(ctx%upper(2, :) - ctx%lower(2, :)), &
            1.0e-12_dp), color=black, &
            linestyle='--', linewidth=1.2_dp, label='bare potential')
        call plot(ctx%r, max(0.5_dp*(ctx%ne_hi - ctx%ne_lo), 1.0e-12_dp), &
            color=vermillion, &
            linestyle='-.', linewidth=1.5_dp, label='electron density')
        call legend(loc='upper right')
        call tight_layout()
    end subroutine run_debye_sphere_demo_draw_frame
end module debye_sphere_demo_callbacks

program debye_sphere_demo
    use debye_sphere_demo_callbacks, only: run_debye_sphere_demo
    implicit none
    call run_debye_sphere_demo()
end program debye_sphere_demo
