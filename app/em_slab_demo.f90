!> Plane electromagnetic wave on a cold electron plasma edge, below and above
!> the electron plasma frequency (docs/em-slab-reflection.md).
!>
!> Runs the two configured frequencies on the same plasma, compares the
!> lock-in reflection/transmission, skin depth and transmitted wavelength
!> with the linear cold-plasma theory (independent oracle leaves), writes
!> output/em_slab_demo.nc and, with --export DIR, the lecture-scene data
!> DIR/em-wave-plasma-reflection/ (kin6d animation export).
program em_slab_demo
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_em_slab, only: em_slab_config_t, em_slab_result_t, em_slab_frames_t, &
        em_slab_run, em_slab_reflection, em_slab_decay_length, &
        em_slab_mean_amplitude, em_slab_wavenumber
    use kin6d_em_slab_leaves, only: em_slab_fresnel, em_slab_evanescent, &
        em_slab_vacuum_envelope, em_slab_evanescent_envelope
    use kin6d_result_io, only: run_t, result_file_t
    use kin6d_animation_export, only: animation_export_t
    implicit none
    real(dp), parameter :: pi = 3.14159265358979323846_dp
    ! Configuration (&em_slab); run 1 below, run 2 above omega_pe.
    real(dp) :: omegas(2) = [0.5_dp, 2.0_dp]
    real(dp) :: amplitude = 5.0e-3_dp, ramp = 25.0_dp, dx = 0.04_dp
    real(dp) :: courant = 0.95_dp
    real(dp) :: x_left = -20.0_dp, x_source = -18.0_dp, x_edge = 0.0_dp
    real(dp) :: x_right = 135.0_dp, probe = -16.0_dp, probe_far = -12.0_dp
    integer :: per_cell = 8
    real(dp) :: window_starts(2) = [75.0_dp, 62.0_dp]
    integer :: window_periods(2) = [4, 8]
    real(dp) :: frame_intervals(2) = [0.25_dp, 0.1_dp]
    real(dp) :: display_min = -15.0_dp, display_max = 10.0_dp
    integer :: display_stride = 2
    real(dp) :: particle_min = 0.0_dp, particle_max = 4.0_dp
    integer :: particle_stride = 16
    namelist /em_slab/ omegas, amplitude, ramp, dx, courant, x_left, x_source, &
        x_edge, x_right, probe, probe_far, per_cell, window_starts, window_periods, &
        frame_intervals, display_min, display_max, display_stride, particle_min, &
        particle_max, particle_stride
    type(run_t) :: run
    type(result_file_t) :: result
    type(em_slab_config_t) :: config
    type(em_slab_result_t) :: res(2)
    type(em_slab_frames_t) :: frames(2)
    character(len=256) :: argument, export_dir, msg
    real(dp) :: theory(6, 2), measured(6, 2), balance(2), gauss(2), a0(2)
    real(dp), allocatable :: envelope(:, :), lockin(:, :)
    integer :: unit, ios, i, r

    call run%start('em_slab_demo')
    call run%open_config(unit)
    if (run%has_config()) read (unit, nml=em_slab, iostat=ios, iomsg=msg)
    call run%close_config(unit, 'em_slab', ios, msg)
    call run%open_resolved(unit)
    write (unit, nml=em_slab)
    call run%close_resolved(unit)
    export_dir = ''
    i = 1
    do while (i <= run%argument_count())
        call run%argument(i, argument)
        if (trim(argument) /= '--export' .or. i == run%argument_count()) &
            error stop 'usage: em_slab_demo [--config FILE] [--export DIR]'
        call run%argument(i + 1, export_dir); i = i + 2
    end do
    if (.not. (omegas(1) < 1 .and. omegas(2) > 1)) &
        error stop 'em_slab_demo: need omegas(1) < 1 < omegas(2) [omega_pe]'

    config%amplitude = amplitude; config%ramp = ramp; config%dx = dx
    config%courant = courant; config%x_left = x_left; config%x_source = x_source
    config%x_edge = x_edge; config%x_right = x_right; config%per_cell = per_cell
    do r = 1, 2
        config%omega = omegas(r)
        frames(r)%interval = frame_intervals(r)
        frames(r)%x_min = display_min; frames(r)%x_max = display_max
        frames(r)%stride = display_stride
        frames(r)%particle_min = particle_min; frames(r)%particle_max = particle_max
        frames(r)%particle_stride = particle_stride
        call em_slab_run(config, window_starts(r), window_periods(r), probe, &
                         probe_far, res(r), frames(r))
        balance(r) = res(r)%balance; gauss(r) = res(r)%gauss_residual
        a0(r) = amplitude/omegas(r)
    end do
    call compare()
    call write_result()
    if (len_trim(export_dir) > 0) call export_animation()

contains

    !> theory/measured rows: |R|, arg R, |T| (edge), skin depth or wavelength,
    !> |R|^2 + N|T|^2 (transmission only), k.
    subroutine compare()
        real(dp) :: kappa, delta, rre, rim, phase, tabs, n, rr, tt, k, lambda
        complex(dp) :: refl
        theory = 0; measured = 0
        call em_slab_evanescent(omegas(1), kappa, delta, rre, rim, phase, tabs)
        refl = em_slab_reflection(res(1), omegas(1), x_edge)
        theory(:, 1) = [hypot(rre, rim), phase, tabs, delta, 1.0_dp, 0.0_dp]
        measured(1:4, 1) = [abs(refl), atan2(aimag(refl), real(refl)), &
                            abs(res(1)%amplitude(node(x_edge)))/abs(res(1)%forward), &
                            em_slab_decay_length(res(1), x_edge + 0.25_dp, &
                                                 x_edge + 3.0_dp)]
        measured(5, 1) = abs(refl)**2
        call em_slab_fresnel(omegas(2), n, rr, tt, k, lambda)
        refl = em_slab_reflection(res(2), omegas(2), x_edge)
        theory(:, 2) = [rr, 0.0_dp, tt, lambda, 1.0_dp, k]
        measured(:, 2) = [abs(refl), atan2(aimag(refl), real(refl)), &
                          em_slab_mean_amplitude(res(2), x_edge + 1, x_edge + 10)/ &
                          abs(res(2)%forward), 0.0_dp, 0.0_dp, &
                          em_slab_wavenumber(res(2), x_edge + 1, x_edge + 10)]
        measured(4, 2) = 2*pi/measured(6, 2)
        measured(5, 2) = measured(1, 2)**2 + n*measured(3, 2)**2
        print '(a,f5.2,a)', 'omega = ', omegas(1), ' omega_pe (reflection)'
        print '(a,2f11.6)', '  |R|            measured, theory: ', measured(1, 1), &
            theory(1, 1)
        print '(a,2f11.6)', '  arg R          measured, theory: ', measured(2, 1), &
            theory(2, 1)
        print '(a,2f11.6)', '  |E(edge)|/E0   measured, theory: ', measured(3, 1), &
            theory(3, 1)
        print '(a,2f11.6)', '  skin depth     measured, theory: ', measured(4, 1), &
            theory(4, 1)
        print '(a,f5.2,a)', 'omega = ', omegas(2), ' omega_pe (transmission)'
        print '(a,2f11.6)', '  R              measured, theory: ', measured(1, 2), &
            theory(1, 2)
        print '(a,2f11.6)', '  T              measured, theory: ', measured(3, 2), &
            theory(3, 2)
        print '(a,2f11.6)', '  wavelength     measured, theory: ', measured(4, 2), &
            theory(4, 2)
        print '(a,2es11.3)', 'energy balance residual / incident: ', balance
        print '(a,2es11.3)', 'discrete Gauss residual:            ', gauss
    end subroutine

    integer function node(x) result(i)
        real(dp), intent(in) :: x
        i = nint((x - x_left)/dx)
    end function

    subroutine write_result()
        integer :: status
        call execute_command_line('mkdir -p output', exitstat=status)
        if (status /= 0) error stop 'em_slab_demo: cannot create output/'
        call result%create(run)
        call result%write('omega', omegas, 'omega_pe', 'wave frequencies', dims='run')
        call result%write('theory', theory, '1', 'linear cold-plasma oracle: '// &
            '|R|, arg R, |T|, skin depth or wavelength, flux sum, k', &
            dims='quantity,run')
        call result%write('measured', measured, '1', 'lock-in estimates '// &
            '(candidate), same rows as theory', dims='quantity,run', &
            status='candidate')
        call result%write('energy_balance', balance, '1', &
            '(incident - reflected - outflow - field - kinetic)/incident', dims='run')
        call result%write('gauss_residual', gauss, 'e n0', &
            'max discrete Gauss-law residual at the end', dims='run')
        call result%close()
    end subroutine

    !> Lecture-scene data (slug em-wave-plasma-reflection): E_y/E0 on the
    !> display grid per frame, sample electron sheets (x - x0, y), the steady
    !> linear-theory envelope and the measured lock-in envelope.
    subroutine export_animation()
        type(animation_export_t) :: export
        character(len=4), parameter :: tag(2) = ['low ', 'high']
        real(dp) :: kappa, delta, rre, rim, phase, tabs, n, rr, tt, k, lambda
        integer :: nx, j, nt, m
        call export%create(run, export_dir, 'em-wave-plasma-reflection')
        call export%note('model', 'relativistic Vlasov-Maxwell, 1D in x, normal '// &
            'incidence: E_x, E_y, B_z on a staggered Yee grid, first-order Mur '// &
            'ends, '// &
            'soft sheet-current source; cold electron macroparticle sheets (CIC, '// &
            'relativistic Boris, charge-conserving J_x), fixed neutralizing ions, '// &
            'uniform density n0 for x > x_edge')
        call export%note('units', 'x in c/omega_pe, t in 1/omega_pe, E_y in E0 '// &
            '(incident amplitude), electron displacement in c/omega_pe')
        call export%note('status', 'candidate PIC dynamics; R, T, skin depth and '// &
            'wavelength are lock-in estimates checked against linear cold-plasma '// &
            'theory (CTest em_slab_*), not certificates')
        call export%note('runs', 'low: omega = omegas(1) < omega_pe (reflection); '// &
            'high: omega = omegas(2) > omega_pe (transmission); same plasma')
        call export%parameter('omega_low', omegas(1), 'omega_pe', 'run 1 frequency')
        call export%parameter('omega_high', omegas(2), 'omega_pe', 'run 2 frequency')
        call export%parameter('amplitude', amplitude, 'm c omega_pe/e', &
            'incident E_y amplitude E0 (both runs)')
        call export%parameter('a0_low', a0(1), '1', 'e E0/(m omega c), run 1')
        call export%parameter('a0_high', a0(2), '1', 'e E0/(m omega c), run 2')
        call export%parameter('x_edge', x_edge, 'c/omega_pe', 'plasma edge')
        call export%parameter('x_source', x_source, 'c/omega_pe', 'source sheet')
        call export%parameter('x_left', x_left, 'c/omega_pe', 'left Mur boundary')
        call export%parameter('x_right', x_right, 'c/omega_pe', 'right Mur boundary')
        call export%parameter('dx', dx, 'c/omega_pe', 'cell width')
        call export%parameter('dt', courant*dx, '1/omega_pe', 'time step')
        call export%parameter('per_cell', real(per_cell, dp), '1', &
            'electron sheets per cell')
        call export%parameter('ramp', ramp, '1/omega_pe', 'sin^2 envelope rise time')
        call export%parameter('R_low_measured', measured(1, 1), '1', '|R|, run 1')
        call export%parameter('R_low_theory', theory(1, 1), '1', '|R| = 1')
        call export%parameter('R_phase_low_measured', measured(2, 1), 'rad', 'arg R')
        call export%parameter('R_phase_low_theory', theory(2, 1), 'rad', &
            '-2 atan(kappa)')
        call export%parameter('skin_depth_measured', measured(4, 1), 'c/omega_pe', &
            'fit of the lock-in |E_y| on [0.25, 3] c/omega_pe inside')
        call export%parameter('skin_depth_theory', theory(4, 1), 'c/omega_pe', &
            'c/sqrt(omega_pe^2 - omega^2)')
        call export%parameter('R_high_measured', measured(1, 2), '1', '|R|, run 2')
        call export%parameter('R_high_theory', theory(1, 2), '1', '(1 - N)/(1 + N)')
        call export%parameter('T_high_measured', measured(3, 2), '1', &
            'mean lock-in |E_y|/E0 on [1, 10] c/omega_pe inside')
        call export%parameter('T_high_theory', theory(3, 2), '1', '2/(1 + N)')
        call export%parameter('wavelength_measured', measured(4, 2), 'c/omega_pe', &
            'lock-in phase slope on [1, 10] c/omega_pe inside')
        call export%parameter('wavelength_theory', theory(4, 2), 'c/omega_pe', &
            '2 pi c/sqrt(omega^2 - omega_pe^2)')
        call export%parameter('energy_balance_low', balance(1), '1', &
            '(in - reflected - outflow - field - kinetic)/in')
        call export%parameter('energy_balance_high', balance(2), '1', &
            '(in - reflected - outflow - field - kinetic)/in')
        nx = size(frames(1)%x)
        call export%real32('x', frames(1)%x - x_edge, [nx], 'c/omega_pe', &
            'display grid, x - x_edge')
        call export%real32('x0', frames(1)%x0 - x_edge, [size(frames(1)%x0)], &
            'c/omega_pe', 'initial positions of the sample electron sheets')
        allocate (envelope(nx, 2), lockin(nx, 2))
        call em_slab_evanescent(omegas(1), kappa, delta, rre, rim, phase, tabs)
        call em_slab_fresnel(omegas(2), n, rr, tt, k, lambda)
        do j = 1, nx
            if (frames(1)%x(j) < x_edge) then
                call em_slab_vacuum_envelope(frames(1)%x(j) - x_edge, omegas(1), &
                                             rre, rim, envelope(j, 1))
                call em_slab_vacuum_envelope(frames(1)%x(j) - x_edge, omegas(2), &
                                             rr, 0.0_dp, envelope(j, 2))
            else
                call em_slab_evanescent_envelope(frames(1)%x(j) - x_edge, omegas(1), &
                                                 envelope(j, 1))
                envelope(j, 2) = tt
            end if
            do m = 1, 2
                lockin(j, m) = abs(res(m)%amplitude(node(frames(1)%x(j))))/ &
                    abs(res(m)%forward)
            end do
        end do
        do m = 1, 2
            nt = size(frames(m)%t)
            call export%parameter('window_start_'//trim(tag(m)), &
                res(m)%window_start, '1/omega_pe', 'lock-in window start')
            call export%parameter('window_end_'//trim(tag(m)), res(m)%window_end, &
                '1/omega_pe', 'lock-in window end')
            call export%real32('t_'//trim(tag(m)), frames(m)%t, [nt], '1/omega_pe', &
                'frame times')
            call export%real32('ey_'//trim(tag(m)), &
                reshape(frames(m)%ey/amplitude, [nx*nt]), [nx, nt], '1', &
                'E_y/E0 on the display grid per frame')
            call export%real32('xdisp_'//trim(tag(m)), &
                reshape(frames(m)%px - spread(frames(m)%x0, 2, nt), &
                        [size(frames(m)%px)]), shape(frames(m)%px), 'c/omega_pe', &
                'longitudinal displacement x - x0 of the sample sheets')
            call export%real32('y_'//trim(tag(m)), &
                reshape(frames(m)%py, [size(frames(m)%py)]), shape(frames(m)%py), &
                'c/omega_pe', 'transverse displacement y of the sample sheets')
            call export%real32('envelope_'//trim(tag(m)), envelope(:, m), [nx], '1', &
                'steady linear cold-plasma |E_y|/E0 (oracle)')
            call export%real32('lockin_'//trim(tag(m)), lockin(:, m), [nx], '1', &
                'measured steady |E_y|/|E_forward| from the lock-in window')
        end do
        call export%close()
        print '(2a)', 'animation data written to ', &
            trim(export_dir)//'/em-wave-plasma-reflection'
    end subroutine
end program em_slab_demo
