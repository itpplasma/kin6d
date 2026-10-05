program quantum_reference_demo
    use, intrinsic :: iso_fortran_env, only: dp => real64, i64 => int64, &
        compiler_version, &
        compiler_options
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortnum_interval, only: interval_t, interval
    use kin6d_quantum_fock, only: assemble_jellium_pair_sector
    use kin6d_species, only: units_t, periodic_pair_units
    use kin6d_quantum_evolution, only: quantum_evolution_t, &
        initialize_quantum_evolution, solve_quantum_evolution, &
        sample_quantum_evolution
    use kin6d_quantum_reference_gallery, only: quantum_observable_packet_t, &
        prepare_quantum_pair_sector, &
        sample_quantum_observables, &
        quantum_occupation_bands, quantum_initial_norm, quantum_matrix_defect, &
        quantum_physical_radius
    use kin6d_quantum_reference_gallery_leaves, only: qgallery_gap_interval, &
        qgallery_radius_interval, qgallery_time_interval
    use fortplot, only: figure, get_global_figure, figure_t, subplot, plot, &
        errorbar, fill_between, title, xlabel, ylabel, xlim, ylim, legend, &
        tight_layout, savefig
    use fortplot_animation, only: animation_t, animation_func => FuncAnimation, &
        save_animation
    use kin6d_result_io, only: run_t, result_file_t
    use kin6d_periodic_pair_config, only: periodic_pair_t, read_periodic_pair, &
        write_periodic_pair, periodic_pair_bins
    implicit none
    integer, parameter :: points = 193, bin_count = periodic_pair_bins
    integer :: momentum_cutoff = 2, particles = 2, frames, steps_per_frame = 1
    ! Box, mass, coupling and time grid come from the &periodic_pair group
    ! shared with classical_pair_demo; hbar stays quantum-only.
    real(dp) :: length, mass, hbar, coupling, time_max
    type(run_t) :: run
    type(result_file_t) :: result
    type(periodic_pair_t) :: pair
    type(units_t) :: units
    integer(i64) :: max_states = 10000_i64, max_bytes = 1073741824_i64
    integer :: initial_terms = 0
    integer, allocatable :: initial_momenta(:, :)
    real(dp), allocatable :: initial_real(:), initial_imag(:)
    namelist /quantum_reference/ momentum_cutoff, particles, steps_per_frame, hbar, &
        max_states, max_bytes
    namelist /quantum_initial/ initial_terms, initial_momenta, initial_real, &
        initial_imag
    type(quantum_evolution_t) :: state(2), candidate(2)
    type(quantum_observable_packet_t) :: packet
    type(interval_t) :: norm0, gap, radius_box, time_box
    type(figure_t), pointer :: canvas
    type(animation_t) :: animation
    integer, allocatable :: mode_momentum(:, :)
    real(dp), allocatable :: initial_re(:), initial_im(:), hr(:, :, :)
    real(dp), allocatable :: hi(:, :, :), re(:), im(:), physical_re(:), physical_im(:)
    type(interval_t), allocatable :: hri(:, :), hii(:, :)
    real(dp), allocatable :: values(:, :, :, :), lower(:, :, :, :), upper(:, :, :, :)
    real(dp), allocatable :: bin_value(:, :, :, :), bin_lower(:, :, :, :)
    real(dp), allocatable :: bin_upper(:, :, :, :)
    real(dp), allocatable :: occupation(:, :, :), occupation_lo(:, :, :)
    real(dp), allocatable :: occupation_hi(:, :, :), orbital(:), err(:)
    real(dp), allocatable :: times(:)
    real(dp) :: x(points), bin_centers(bin_count), widths(3), matrix_defect(2)
    real(dp) :: before, after, setup_s, strict_setup_s, candidate_s, certified_s
    real(dp) :: common_setup_s
    real(dp) :: post_s, radius, physical_radius, max_gap, max_radius, h, total_s
    real(dp) :: ylim_pair(2), ylim_connected(2)
    integer :: modes, bound, shell, kx, ky, kz, i, j, k, q
    integer :: f, model_index, cell, unit, status, argn, snapshot
    integer :: refinement
    logical :: ok, bounded, snapshot_only, accepted, solved
    character(512) :: argument
    character(256) :: message
    call run%start('quantum_reference_demo')
    snapshot_only = .false.; snapshot = -1
    argn = 1
    do while (argn <= run%argument_count())
        call run%argument(argn, argument)
        select case (trim(argument))
        case ('--snapshot')
            snapshot_only = .true.
            if (argn < run%argument_count()) then
                call run%argument(argn + 1, argument)
                if (argument(1:2) /= '--') then
                    read (argument, *, iostat=status) snapshot
                    if (status /= 0) error stop 'invalid snapshot frame'
                    argn = argn + 1
                end if
            end if
        case default
            error stop 'use --snapshot [frame] and/or --config filename'
        end select
        argn = argn + 1
    end do
    units = periodic_pair_units()
    hbar = units%hbar
    call read_periodic_pair(run, pair)
    length = pair%length; mass = pair%mass; coupling = pair%coupling
    time_max = pair%time_max; frames = pair%frames
    status = 0; message = ''
    call run%open_config(unit)
    if (run%has_config()) read (unit, nml=quantum_reference, iostat=status, &
        iomsg=message)
    call run%close_config(unit, 'quantum_reference', status, message)
    if (particles /= 2) error stop 'this matched gallery is the N=2 spin-pair sector'
    if (momentum_cutoff < 0 .or. momentum_cutoff > 16) &
        error stop 'invalid spatial momentum cutoff'
    if (frames < 2 .or. frames > 1001) error stop 'invalid frame budget'
    if (steps_per_frame < 1 .or. steps_per_frame > 1024) &
        error stop 'invalid time-step budget'
    if (.not. all(ieee_is_finite([length, mass, hbar, coupling, time_max]))) &
        error stop 'nonfinite parameters'
    if (min(length, mass, hbar, time_max) <= 0.0_dp .or. coupling < 0.0_dp) &
        error stop 'invalid physical parameters'
    if (snapshot < 0) snapshot = frames
    if (snapshot < 1 .or. snapshot > frames) error stop 'snapshot frame range'
    call execute_command_line('mkdir -p output', exitstat=status)
    if (status /= 0) error stop 'cannot create output'
    call cpu_time(before)
    bound = momentum_cutoff; modes = (2*bound + 1)**3
    if (int(modes, i64) > max_states) error stop 'sector state budget exceeded'
    if (64_i64*int(modes, i64)**2 > max_bytes) error stop 'matrix byte budget exceeded'
    allocate (mode_momentum(3, modes))
    k = 0
    do shell = 0, 3*bound**2
        do kz = -bound, bound
            do ky = -bound, bound
                do kx = -bound, bound
                    if (kx*kx + ky*ky + kz*kz /= shell) cycle
                    k = k + 1; mode_momentum(:, k) = [kx, ky, kz]
                end do
            end do
        end do
    end do
    allocate (initial_momenta(3, modes), initial_real(modes), initial_imag(modes))
    initial_momenta = 0; initial_real = 0.0_dp; initial_imag = 0.0_dp
    status = 0; message = ''
    call run%open_config(unit)
    if (run%has_config()) read (unit, nml=quantum_initial, iostat=status, &
        iomsg=message)
    call run%close_config(unit, 'quantum_initial', status, message)
    if (initial_terms < 0 .or. initial_terms > modes) error stop 'initial term count'
    ! Early capture rejects unknown groups before the expensive solve; the
    ! final capture in write_outputs records the adapted steps_per_frame.
    call capture_resolved()
    open (newunit=unit, file='output/quantum-input.nml', status='replace', &
        access='stream', form='unformatted')
    if (run%has_config()) then
        write (unit) run%literal
    else
        write (unit) run%resolved
    end if
    close (unit)
    allocate (initial_re(modes), initial_im(modes), re(modes), im(modes), &
        physical_re(modes), physical_im(modes))
    initial_re = 0.0_dp; initial_im = 0.0_dp
    if (initial_terms == 0) then
        initial_terms = 1; initial_momenta(:, 1) = 0; initial_real(1) = 1.0_dp
    end if
    do i = 1, initial_terms
        j = 0
        do k = 1, modes
            if (all(mode_momentum(:, k) == initial_momenta(:, i))) j = k
        end do
        if (j < 1) error stop 'initial momentum outside retained basis'
        if (.not. all(ieee_is_finite([initial_real(i), initial_imag(i)]))) &
            error stop 'nonfinite initial amplitude'
        do k = 1, i - 1
            if (all(initial_momenta(:, k) == initial_momenta(:, i))) error stop &
                'duplicate initial momentum'
        end do
        initial_re(j) = initial_real(i); initial_im(j) = initial_imag(i)
    end do
    call quantum_initial_norm(initial_re, initial_im, norm0, ok)
    if (.not. ok) error stop 'zero or invalid initial state'
    allocate (hr(modes, modes, 2), hi(modes, modes, 2), &
        hri(modes, modes), hii(modes, modes))
    call cpu_time(after)
    common_setup_s = after - before
    call cpu_time(before)
    do model_index = 1, 2
        radius = coupling
        if (model_index == 2) radius = 0.0_dp
        call assemble_jellium_pair_sector(mode_momentum, length, hbar, mass, &
            radius, max_bytes, hr(:, :, model_index), hi(:, :, model_index), &
            hri, hii, ok, message, certify=.false.)
        if (.not. ok) error stop 'candidate Hamiltonian construction failed'
    end do
    call cpu_time(after)
    setup_s = after - before
    call cpu_time(before)
    do model_index = 1, 2
        radius = coupling
        if (model_index == 2) radius = 0.0_dp
        call assemble_jellium_pair_sector(mode_momentum, length, hbar, mass, &
            radius, max_bytes, hr(:, :, model_index), hi(:, :, model_index), &
            hri, hii, ok, message)
        if (.not. ok) error stop trim(message)
        call quantum_matrix_defect(hr(:, :, model_index), hi(:, :, model_index), &
            hri, hii, matrix_defect(model_index), ok)
        if (.not. ok) error stop 'physical Hamiltonian rounding budget failed'
    end do
    call cpu_time(after)
    strict_setup_s = after - before
    ! A row-norm scale avoids obviously unstable initial time steps.
    radius = maxval(sum(abs(hr(:, :, 1)) + abs(hi(:, :, 1)), dim=2))/hbar
    radius = time_max/real(frames - 1, dp)*radius/0.5_dp
    if (.not. ieee_is_finite(radius)) error stop 'time-step scale overflow'
    if (radius > 1024.0_dp) error stop 'initial time-step work budget exceeded'
    steps_per_frame = max(steps_per_frame, ceiling(radius))
    allocate (times(frames), values(6, points, frames, 2), lower(6, points, &
        frames, 2), &
        upper(6, points, frames, 2), occupation(modes, frames, 2), &
        occupation_lo(modes, frames, 2), occupation_hi(modes, frames, 2), &
        orbital(modes), err(modes), bin_value(6, bin_count, frames, 2), &
        bin_lower(6, bin_count, frames, 2), bin_upper(6, bin_count, frames, 2))
    do i = 1, points
        x(i) = -0.5_dp + real(i - 1, dp)/real(points - 1, dp)
    end do
    do i = 1, bin_count
        bin_centers(i) = -0.5_dp + (real(i, dp) - 0.5_dp)/real(bin_count, dp)
    end do
    widths = [1.0_dp/real(bin_count, dp), 1.0_dp, 1.0_dp]
    max_gap = 0.0_dp
    do i = 1, points - 1
        call qgallery_gap_interval(interval(x(i)), interval(x(i + 1)), gap)
        max_gap = max(max_gap, gap%hi)
    end do
    do i = 1, modes
        orbital(i) = real(i, dp)
    end do
    accepted = .false.
    total_s = 0.0_dp
    do refinement = 0, 7
        if (steps_per_frame > 1024) exit
        if (64_i64*int(modes, i64)**2 + 1024_i64*int(modes, i64)* &
            int((frames - 1)*steps_per_frame, i64) > max_bytes) &
            error stop 'time-trajectory working-byte budget exceeded'
        h = time_max/real((frames - 1)*steps_per_frame, dp)
        do f = 1, frames
            times(f) = real((f - 1)*steps_per_frame, dp)*h
        end do
        call cpu_time(before)
        solved = .true.
        do model_index = 1, 2
            call initialize_quantum_evolution(state(model_index), &
                hr(:, :, model_index), hi(:, :, model_index), initial_re, &
                initial_im, &
                (frames - 1)*steps_per_frame, h, &
                hbar, ok)
            if (.not. ok) error stop 'time evolution initialization failed'
            call solve_quantum_evolution(state(model_index), ok)
            if (.not. ok) then
                solved = .false.; exit
            end if
        end do
        call cpu_time(after)
        certified_s = after - before; total_s = total_s + certified_s
        if (.not. solved) then
            steps_per_frame = 2*steps_per_frame
            cycle
        end if
        call cpu_time(before)
        max_radius = 0.0_dp
        do model_index = 1, 2
            do f = 1, frames
                call frame_packet(model_index, f, max_gap)
                call quantum_occupation_bands(packet, occupation(:, f, model_index), &
                    occupation_lo(:, f, model_index), occupation_hi(:, f, &
                    model_index), ok)
                if (.not. ok) error stop 'occupation enclosure failed'
                do i = 1, modes
                    call qgallery_radius_interval(interval(occupation_lo(i, f, &
                        model_index)), &
                        interval(occupation_hi(i, f, model_index)), &
                        radius_box)
                    max_radius = max(max_radius, radius_box%hi)
                end do
                do i = 1, points
                    call sample_quantum_observables(packet, x(i), &
                        values(:, i, f, model_index), lower(:, i, f, model_index), &
                        upper(:, i, f, model_index), ok)
                    if (.not. ok) error stop 'point/bin observable certificate failed'
                    do q = 1, 6
                        call qgallery_radius_interval(interval(lower(q, i, f, &
                            model_index)), &
                            interval(upper(q, i, f, model_index)), &
                            radius_box)
                        max_radius = max(max_radius, radius_box%hi)
                    end do
                end do
                do i = 1, bin_count
                    call sample_quantum_observables(packet, bin_centers(i), &
                        bin_value(:, i, f, model_index), bin_lower(:, i, f, &
                        model_index), &
                        bin_upper(:, i, f, model_index), ok)
                    if (.not. ok) error stop 'matched-bin certificate unavailable'
                end do
            end do
        end do
        call cpu_time(after)
        post_s = after - before; total_s = total_s + post_s
        if (max_radius <= 0.010000000001_dp) then
            accepted = .true.; exit
        end if
        steps_per_frame = 2*steps_per_frame
    end do
    if (.not. accepted) error stop 'observable 1% radius refinement budget exhausted'
    call cpu_time(before)
    do model_index = 1, 2
        call initialize_quantum_evolution(candidate(model_index), hr(:, :, &
            model_index), &
            hi(:, :, model_index), initial_re, initial_im, &
            (frames - 1)*steps_per_frame, h, hbar, ok)
        if (.not. ok) error stop 'candidate-only initialization failed'
        call solve_quantum_evolution(candidate(model_index), ok, certify=.false.)
        if (.not. ok) error stop 'candidate-only solve failed'
    end do
    call cpu_time(after)
    candidate_s = after - before
    call write_outputs()
    ylim_pair = [minval(lower(1, :, :, :)), maxval(upper(1, :, :, :))]
    ylim_pair = [min(0.0_dp, ylim_pair(1) - 0.05_dp), ylim_pair(2) + 0.4_dp]
    ylim_connected = [minval(lower([2, 4], :, :, :)), maxval(upper([2, 4], :, :, :))]
    ylim_connected = [ylim_connected(1) - 0.08_dp, ylim_connected(2) + 0.3_dp]
    call figure(figsize=[12.8_dp, 8.2_dp])
    canvas => get_global_figure()
    call draw_frame(snapshot)
    call savefig('output/quantum-correlations.png')
    if (.not. snapshot_only) then
        animation = animation_func(update, frames=frames, interval=83, fig=canvas)
call save_animation(animation, 'output/quantum-correlations.mp4', fps=12, status=status)
        if (status /= 0) error stop 'quantum movie encoding failed'
    end if
    print '(a,i0,a,i0,a,es12.4)', 'quantum3D modes=', modes, ', sector=', modes, &
        ', maximum observable radius=', max_radius
contains
    subroutine write_result()
        character(len=*), parameter :: bound = 'certified-discrete finite-basis '// &
            'enclosure; continuum cutoff/model/input uncertainty unavailable'
        character(len=*), parameter :: dims = 'x,frame,model'
        character(24), parameter :: names(6) = [character(24) :: 'pair_up_down', &
            'cumulant_up_down', 'pair_up_up', 'cumulant_up_up', 'density_up', &
            'density_down']
        integer :: qq
        call result%create(run)
        call result%attribute('model_index', '1=Coulomb, 2=noninteracting control')
        call result%write('time', times, 'periodic-pair time unit', 'frame time', &
            dims='frame')
        call result%write('x', x, '1', 'relative coordinate x/L', dims='x')
        call result%write('bin_center', bin_centers, '1', &
            'relative-x slab center x/L (width 1/16)', dims='bin')
        do qq = 1, 6
            call result%write(trim(names(qq)), values(qq, :, :, :), '1', &
                trim(names(qq))//' center', dims=dims, status='candidate center')
            call result%write(trim(names(qq))//'_lower', lower(qq, :, :, :), '1', &
                trim(names(qq))//' lower bound', dims=dims, status=bound)
            call result%write(trim(names(qq))//'_upper', upper(qq, :, :, :), '1', &
                trim(names(qq))//' upper bound', dims=dims, status=bound)
        end do
        call result%write('bin_g', bin_value(1, :, :, :), '1', 'slab pair correlation g', &
            dims='bin,frame,model', status='candidate center')
        call result%write('bin_g_lower', bin_lower(1, :, :, :), '1', 'g lower bound', &
            dims='bin,frame,model', status=bound)
        call result%write('bin_g_upper', bin_upper(1, :, :, :), '1', 'g upper bound', &
            dims='bin,frame,model', status=bound)
        call result%write('occupation', occupation, '1', 'paired-orbital occupation', &
            dims='orbital,frame,model', status='candidate center')
        call result%write('occupation_lower', occupation_lo, '1', 'occupation lower', &
            dims='orbital,frame,model', status=bound)
        call result%write('occupation_upper', occupation_hi, '1', 'occupation upper', &
            dims='orbital,frame,model', status=bound)
        call result%write('max_radius', max_radius, '1', &
            'maximum certified observable radius', status=bound)
        call result%close()
    end subroutine

    subroutine capture_resolved()
        integer :: u
        call run%open_resolved(u)
        call write_periodic_pair(u, pair)
        write (u, nml=quantum_reference)
        write (u, nml=quantum_initial)
        call run%close_resolved(u)
    end subroutine

    subroutine frame_packet(model_number, frame_number, display_gap)
        integer, intent(in) :: model_number, frame_number
        real(dp), intent(in) :: display_gap
        real(dp) :: fraction
        cell = 1; fraction = 0.0_dp
        if (frame_number > 1) then
            cell = (frame_number - 1)*steps_per_frame; fraction = 1.0_dp
        end if
        call sample_quantum_evolution(state(model_number), cell, fraction, &
            re, im, radius, ok, bounded)
        if (.not. ok .or. .not. bounded) error stop &
            'unavailable quantum state enclosure'
        call qgallery_time_interval(interval(real((frame_number - 1)*steps_per_frame, &
            dp)), &
            interval(h), time_box)
        call quantum_physical_radius(norm0, radius, matrix_defect(model_number), &
            time_box%hi, hbar, physical_radius, ok)
        if (.not. ok) error stop 'physical coefficient rounding propagation failed'
        physical_re = re; physical_im = im
        call prepare_quantum_pair_sector(packet, mode_momentum, physical_re, &
            physical_im, norm0, physical_radius, display_gap, &
            widths, ok)
        if (.not. ok) error stop 'quantum observable packet failed'
    end subroutine
    subroutine write_outputs()
        integer :: u, mi, ff, ii, qq, b
        character(24), parameter :: names(6) = [character(24) :: &
            'pair_up_down', 'cumulant_up_down', &
            'pair_up_up', &
            'cumulant_up_up', 'density_up', &
            'density_down']
        call capture_resolved()
        open (newunit=u, file='output/quantum-resolved.nml', status='replace', &
            access='stream', form='unformatted')
        write (u) run%resolved
        close (u)
        call write_result()
        open (newunit=u, file='output/quantum-observables.csv', status='replace')
        write (u, &
            '(a)') '# normalized exact3D finite-basis jellium spin-pair K=0 sector'
        write (u, '(a)') '# x=x/L; full-transverse relative-bin average, width_x/L=1/16'
    write (u, '(a)') '# certified finite-basis time/observable/physical-H rounding;'// &
        ' continuum cutoff/model/input uncertainty unavailable'
        write (u, '(a)') '# smooth generated leaf JVP/VJP;'// &
            ' observable sensitivity bands unavailable'
        write (u, '(a)') 'model,frame,t,x,observable,value,lower,upper'
        do mi = 1, 2
            do ff = 1, frames
                do ii = 1, points
                    do qq = 1, 6
                        write (u, &
                            '(i0,",",i0,2(",",es24.16e3),",",a,3(",",es24.16e3))') &
                            mi, ff, times(ff), x(ii), trim(names(qq)), &
                            values(qq, ii, ff, mi), lower(qq, ii, ff, mi), upper(qq, &
                                ii, ff, mi)
                    end do
                end do
            end do
        end do
        close (u)
        open (newunit=u, file='output/quantum-pair-bins.csv', status='replace')
        write (u, '(a)') '# identical classical/quantum3D bins: widths/L=(1/16,1,1)'
        write (u, '(a)') '# model1=Coulomb,model2=noninteracting;'// &
            ' deterministic finite-basis bands'
write (u, '(a)') 'model,frame,t,bin,x_center_over_L,g,g_lower,g_upper,C,C_lower,C_upper'
        do mi = 1, 2
            do ff = 1, frames
                do b = 1, bin_count
                    write (u, '(i0,",",i0,",",es24.16e3,",",i0,7(",",es24.16e3))') &
                        mi, ff, times(ff), b, bin_centers(b), bin_value(1, b, ff, mi), &
                        bin_lower(1, b, ff, mi), bin_upper(1, b, ff, mi), &
                        bin_value(2, b, ff, mi), bin_lower(2, b, ff, mi), bin_upper(2, &
                            b, &
                        ff, mi)
                end do
            end do
        end do
        close (u)
        open (newunit=u, file='output/quantum-occupations.csv', status='replace')
        write (u, &
            '(a)') '# paired orbitals k_up,-k_down have equal occupation; Pauli scale1'
        write (u, '(a)') 'model,frame,t,orbital_pair,kx,ky,kz,occupation,lower,upper'
        do mi = 1, 2
            do ff = 1, frames
                do ii = 1, modes
                    write (u, '(i0,",",i0,",",es24.16e3,4(",",i0),3(",",es24.16e3))') &
                        mi, ff, times(ff), ii, mode_momentum(:, ii), occupation(ii, &
                        ff, mi), &
                        occupation_lo(ii, ff, mi), occupation_hi(ii, ff, mi)
                end do
            end do
        end do
        close (u)
        open (newunit=u, file='output/quantum-performance.csv', status='replace')
        write (u, '(a)') '# single local CPU trial; both Coulomb/control;'// &
            ' matched final steps;'// &
            ' initializer+solve baseline; rendering/CSV excluded'
    write (u, '(a)') 'candidate_init_solve_s,certified_init_solve_s,extra_s,factor,'// &
        'observable_post_s,all_refinement_s,max_radius,steps_per_frame'
        write (u, '(7(es24.16e3,","),i0)') candidate_s, certified_s, &
            certified_s - candidate_s, certified_s/max(candidate_s, tiny(1.0_dp)), &
            post_s, total_s, max_radius, steps_per_frame
        write (u, &
            '(a,es24.16e3)') '# candidate_matrices_setup_s=', setup_s
        write (u, '(a,es24.16e3)') '# strict_matrices_setup_s=', strict_setup_s
        write (u, '(a,2es24.16e3)') '# physical_H_matrix_norm_defects=', matrix_defect
        write (u, '(a,es24.16e3)') '# final_certified_solve_and_all_bands_s=', &
            certified_s + post_s
        write (u, '(a,es24.16e3)') '# all_bands_extra_over_solve_s=', &
            certified_s + post_s - candidate_s
        if (candidate_s > 0.0_dp) write (u, '(a,es24.16e3)') &
            '# all_bands_over_solve_factor=', &
            (certified_s + post_s)/candidate_s
        write (u, '(a,es24.16e3)') '# cold_candidate_first_result_s=', &
            common_setup_s + setup_s + candidate_s
        write (u, '(a,es24.16e3)') '# cold_certified_first_result_with_bands_s=', &
            common_setup_s + strict_setup_s + certified_s + post_s
        write (u, '(a,es24.16e3)') '# adaptive_first_result_with_bands_s=', &
            common_setup_s + strict_setup_s + total_s
        write (u, &
            '(a)') '# requested derivatives not timed; gradient bounds unavailable'
        write (u, '(a)') '# compiler='//compiler_version()
        write (u, '(a)') '# flags='//compiler_options()
        close (u)
    end subroutine

    subroutine update(number)
        integer, intent(in) :: number
        call canvas%clear()
        call draw_frame(number)
        call canvas%set_rendered(.false.)
    end subroutine

    subroutine draw_frame(number)
        integer, intent(in) :: number
        character(120) :: heading
        real(dp), parameter :: blue(3) = [0.0_dp, 0.4471_dp, 0.6980_dp]
        real(dp), parameter :: orange(3) = [0.8353_dp, 0.3686_dp, 0.0_dp]
        real(dp), parameter :: gray(3) = [0.45_dp, 0.45_dp, 0.45_dp]
        call subplot(2, 2, 1)
        write (heading, '(a,f6.2,a,i0)') 'Coulomb quench | t/t0 = ', times(number), &
            ' | 3D modes = ', modes
        call title(trim(heading))
        call xlabel('Orbital pair: k up, -k down (discrete index)')
        call ylabel('Occupation per spin orbital')
        call xlim(0.0_dp, real(modes + 1, dp)); call ylim(-0.05_dp, 1.3_dp)
        err = max(occupation(:, number, 1) - occupation_lo(:, number, 1), &
            occupation_hi(:, number, 1) - occupation(:, number, 1))
        call errorbar(orbital, occupation(:, number, 1), yerr=err, fmt='o', &
            color=blue, ecolor=blue, label='Coulomb + deterministic bounds', &
                capsize=2.0_dp)
        err = max(occupation(:, number, 2) - occupation_lo(:, number, 2), &
            occupation_hi(:, number, 2) - occupation(:, number, 2))
        call errorbar(orbital, occupation(:, number, 2), yerr=err, fmt='x', &
            linestyle=':', color=gray, ecolor=gray, &
            label='noninteracting, same engine', capsize=2.0_dp)
        call legend(loc='upper right')
        call subplot(2, 2, 2)
        call title('Unlike-spin relative pair density')
        call xlabel('Relative bin center x / L; full transverse average')
        call ylabel('g(up,down) = V times relative-bin density')
        call xlim(-0.5_dp, 0.5_dp); call ylim(ylim_pair(1), ylim_pair(2))
        call band(1, number, blue, 'Coulomb + finite-basis time band', '-')
        call band(1, number, gray, 'noninteracting, same preparation', ':', 2)
        call legend(loc='upper right')
        call subplot(2, 2, 3)
        call title('Beyond-Wick pair cumulants')
        call xlabel('Relative bin center x / L; full transverse average')
        call ylabel('V^2 times connected pair cumulant')
        call xlim(-0.5_dp, 0.5_dp); call ylim(ylim_connected(1), ylim_connected(2))
        call band(2, number, blue, 'up,down: g - 1', '-')
        call band(4, number, orange, 'up,up: fixed-number correlation', '--')
        call band(2, number, gray, 'noninteracting Slater', ':', 2)
        call legend(loc='upper right')
        call subplot(2, 2, 4)
        call title('Homogeneous one-body density remains unchanged')
        call xlabel('x / L; periodic3D box')
        call ylabel('V times density per spin')
        call xlim(-0.5_dp, 0.5_dp); call ylim(0.85_dp, 1.25_dp)
        call band(5, number, blue, 'up: exact symmetry invariant', '-')
        call band(6, number, orange, 'down: exact symmetry invariant', '--')
        call legend(loc='upper right')
        call tight_layout()
    end subroutine

    subroutine band(channel, number, color, label, style, model_number)
        integer, intent(in) :: channel, number
        real(dp), intent(in) :: color(3)
        character(*), intent(in) :: label, style
        integer, intent(in), optional :: model_number
        integer :: mi
        mi = 1
        if (present(model_number)) mi = model_number
        call fill_between(x, lower(channel, :, number, mi), upper(channel, :, &
            number, mi), &
            color=color, alpha=0.25_dp)
        call plot(x, lower(channel, :, number, mi), color=color, linewidth=0.5_dp)
        call plot(x, upper(channel, :, number, mi), color=color, linewidth=0.5_dp)
        call plot(x, values(channel, :, number, mi), linestyle=style, color=color, &
            label=label)
    end subroutine
end program quantum_reference_demo
