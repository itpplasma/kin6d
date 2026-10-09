program gs_phase1_main
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64, input_unit, output_unit
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use kin6d_grad_shafranov, only: gs_state_t, initialize_gs, solve_gs, &
        write_gs_mesh, read_gs_mesh, sample_gs
    use kin6d_gs_field_readback, only: sample_gs_many
    use kin6d_curved_triangle, only: periodic_curve_t, initialize_periodic_curve, &
        curve_value
    use kin6d_gs_readback, only: gs_readback_t, read_gs_equilibrium, gs_profile_f, &
        gs_readback_resolution, gs_readback_q
    use kin6d_result_io, only: run_t, result_file_t, runtime_git
    use fortplot_system_runtime, only: create_directory_runtime
    implicit none
    type(gs_state_t) :: state
    type(periodic_curve_t) :: curve
    type(gs_readback_t) :: readback
    type(run_t) :: run
    type(result_file_t) :: result
    real(dp), allocatable :: samples(:, :), contour(:, :)
    real(dp) :: max_area = 0, pprime = 0, ffprime = 0, f_edge = 0
    real(dp) :: mu0 = 1.2566370614359173e-6_dp, wall_s, perimeter
    integer :: degree = 2
    integer :: boundary_nodes = 96, radial_points = 0, angular_points = 0
    integer :: unit, status, n, i, next
    integer(int64) :: start, finish, rate
    logical :: curved = .true., directory_ok
    character(512) :: boundary_file = '', mesh_file = 'mesh.dat'
    character(256) :: message
    character(512) :: argument
    character(:), allocatable :: commit, dirty
    namelist /gs_phase1/ boundary_file, mesh_file, boundary_nodes, max_area, &
        pprime, ffprime, f_edge, mu0, curved, radial_points, angular_points, degree
    call run%start('gs_phase1')
    if (run%argument_count() > 0) then
        call run%argument(1, argument)
        select case (trim(argument))
        case ('--evaluate', '--evaluate-stream')
            call evaluate_file(trim(argument) == '--evaluate-stream')
        case ('--q')
            call evaluate_q_file()
        case default
            error stop 'expected --config, --evaluate or --q'
        end select
        stop
    end if
    status = 0; message = ''
    call run%open_config(unit)
    if (run%has_config()) read (unit, nml=gs_phase1, iostat=status, iomsg=message)
    call run%close_config(unit, 'gs_phase1', status, message)
    if (len_trim(boundary_file) == 0) error stop 'boundary_file required'
    if (degree /= 2 .and. degree /= 3) error stop 'degree must be 2 or 3'
    if (boundary_nodes < 8) error stop 'at least eight boundary nodes required'
    call system_clock(start, rate)
    open (newunit=unit, file=trim(boundary_file), status='old', action='read')
    read (unit, *) n
    if (n < 4) error stop 'boundary requires four samples without closing duplicate'
    allocate (samples(2, n), contour(2, boundary_nodes))
    do i = 1, n
        read (unit, *) samples(:, i)
    end do
    close (unit)
    call initialize_periodic_curve(curve, samples, status)
    if (status /= 0) error stop 'invalid closed curve'
    do i = 1, boundary_nodes
        call boundary(real(i - 1, dp)/boundary_nodes, contour(:, i))
    end do
    if (max_area <= 0) then
        perimeter = 0
        do i = 1, boundary_nodes
            next = mod(i, boundary_nodes) + 1
            perimeter = perimeter + norm2(contour(:, next) - contour(:, i))
        end do
        max_area = (perimeter/boundary_nodes)**2/2
    end if
    if (curved) then
        call initialize_gs(state, contour, max_area, status, degree, curve=boundary)
    else
        call initialize_gs(state, contour, max_area, status, degree)
    end if
    if (status /= 0) error stop 'mesh/assembly failed'
    call solve_gs(state, mu0, pprime, ffprime, status)
    if (status /= 0) error stop 'linear equilibrium solve failed'
    call system_clock(finish)
    wall_s = real(finish - start, dp)/rate
    if (wall_s > 300) error stop 'solve exceeded 300 seconds'
    call gs_readback_resolution(state, radial_points, angular_points)
    call read_gs_equilibrium(state, f_edge, ffprime, radial_points, angular_points, &
        readback, status)
    if (status /= 0) then
        print *, 'readback status:', status
        error stop 'equilibrium contour readback failed'
    end if
    call write_gs_mesh(state, trim(mesh_file), status)
    if (status /= 0) error stop 'mesh export failed'
    call run%open_resolved(unit)
    write (unit, nml=gs_phase1)
    call run%close_resolved(unit)
    call create_directory_runtime('output', directory_ok)
    if (.not. directory_ok) error stop 'cannot create output directory'
    call result%create(run)
    call result%attribute('convention', 'COCOS 3; SI; psi in Wb/rad; psi_edge=0')
    call result%write('boundary_samples', samples, 'm', 'supplied closed R Z curve')
    call result%write('points', state%points, 'm', &
        'physical isoparametric geometry nodes')
    call result%write('psi', state%psi, 'Wb/rad', 'native FE coefficients')
    call result%write('triangles', real(state%triangles, dp), '1', &
        'isoparametric connectivity')
    call result%write('axis', readback%axis, 'm', 'physical magnetic axis')
    call result%write('volume', readback%volume, 'm3', 'revolved plasma volume')
    call result%write('phi_edge', readback%phi_edge, 'Wb', 'signed full toroidal flux')
    call result%write('s_tor', readback%s_tor, '1', 'normalized toroidal flux')
    call result%write('q', readback%q, '1', 'signed safety factor on FE contours')
    call result%write('profile_psi', readback%psi, 'Wb/rad', 'native contour levels')
    call result%write('contour_s', readback%contour_s, 's', 'contour readback time')
    call result%write('integral_s', readback%integral_s, 's', 'area integrals time')
    call result%write('wall_s', wall_s, 's', 'mesh assembly and solve wall time')
    call result%write('dof', real(count(.not. state%boundary), dp), '1', &
        'free unknowns')
    call result%write('iterations', real(state%iterations, dp), '1', 'CG iterations')
    call result%write('residual', state%residual, '1', 'relative weak residual')
    call result%write('flux_consistency', readback%flux_consistency, '1', &
        'relative contour versus area flux quadrature difference')
    call result%close()
    open (newunit=unit, file='profiles.dat', status='replace')
    write (unit, '(a)') '# psi[Wb/rad] s_tor q; COCOS 3'
    do i = 1, radial_points
        write (unit, '(3es25.16e3)') readback%psi(i), readback%s_tor(i), readback%q(i)
    end do
    close (unit)
    call runtime_git(commit, dirty, 'gs_phase1')
    open (newunit=unit, file='summary.dat', status='replace')
    write (unit, '(a,1x,i0)') 'degree', state%element_order
    write (unit, '(a,1x,i0)') 'iterations', state%iterations
    write (unit, '(a,1x,i0)') 'dof', count(.not. state%boundary)
    write (unit, '(a,1x,es25.16e3)') 'wall_s', wall_s
    write (unit, '(a,1x,es25.16e3)') 'contour_s', readback%contour_s
    write (unit, '(a,1x,es25.16e3)') 'integral_s', readback%integral_s
    write (unit, '(a,1x,i0)') 'radial_points', radial_points
    write (unit, '(a,1x,i0)') 'angular_points', angular_points
    write (unit, '(a,1x,es25.16e3)') 'axis_R', readback%axis(1)
    write (unit, '(a,1x,es25.16e3)') 'axis_Z', readback%axis(2)
    write (unit, '(a,1x,es25.16e3)') 'volume', readback%volume
    write (unit, '(a,1x,es25.16e3)') 'phi_edge', readback%phi_edge
    write (unit, '(a,1x,es25.16e3)') 'residual', state%residual
    write (unit, '(a,1x,es25.16e3)') 'flux_consistency', readback%flux_consistency
    write (unit, '(a,1x,a)') 'code_commit', commit
    write (unit, '(a,1x,a)') 'code_tree', dirty
    close (unit)
    print '(a,i0,a,es13.5)', 'DOF=', count(.not. state%boundary), ' wall_s=', wall_s
contains
    subroutine boundary(t, point)
        real(dp), intent(in) :: t
        real(dp), intent(out) :: point(2)
        call curve_value(curve, t, point)
    end subroutine

    subroutine evaluate_q_file()
        character(512) :: path
        real(dp) :: row(3), s
        integer :: n, i, unit, stat
        if (run%argument_count() /= 2) &
            error stop '--q profiles.dat; stdin: count, s_tor'
        call run%argument(2, path)
        open (newunit=unit, file=trim(path), status='old', action='read')
        read (unit, *)
        n = 0
        do
            read (unit, *, iostat=stat) row
            if (stat /= 0) exit
            n = n + 1
        end do
        if (n < 4) error stop 'incomplete q profile'
        rewind (unit)
        read (unit, *)
        allocate (readback%psi(n), readback%s_tor(n), readback%q(n))
        do i = 1, n
            read (unit, *) readback%psi(i), readback%s_tor(i), readback%q(i)
        end do
        close (unit)
        read (input_unit, *) n
        do i = 1, n
            read (input_unit, *) s
            write (*, '(es25.16e3)') gs_readback_q(readback, s)
        end do
    end subroutine

    subroutine evaluate_file(streaming)
        logical, intent(in) :: streaming
        character(512) :: path, text
        real(dp) :: f, nan
        real(dp), allocatable :: points(:, :), values(:, :)
        integer, allocatable :: statuses(:)
        integer :: count, i, stat
        if (run%argument_count() /= 4) &
            error stop '--evaluate mesh_file F_edge FFprime; stdin: count then R Z'
        call run%argument(2, path)
        call run%argument(3, text)
        read (text, *) f_edge
        call run%argument(4, text)
        read (text, *) ffprime
        call read_gs_mesh(state, path, stat)
        if (stat /= 0) error stop 'invalid native field export'
        nan = ieee_value(0.0_dp, ieee_quiet_nan)
        do
            read (input_unit, *, iostat=stat) count
            if (stat < 0) exit
            if (stat /= 0 .or. count < 0) error stop 'invalid point count'
            if (allocated(points)) then
                if (size(points, 2) /= count) deallocate (points, values, statuses)
            end if
            if (.not. allocated(points)) &
                allocate (points(2, count), values(4, count), statuses(count))
            do i = 1, count
                read (input_unit, *) points(:, i)
            end do
            call sample_gs_many(state, points, f_edge, values, statuses)
            do i = 1, count
                if (statuses(i) == 0) then
                    f = gs_profile_f(values(1, i), f_edge, ffprime)
                    if (abs(f) <= tiny(1.0_dp)) statuses(i) = 1
                    values(4, i) = f/points(1, i)
                end if
                if (statuses(i) /= 0) values(:, i) = nan
                write (*, '(4es25.16e3)') values(:, i)
            end do
            flush (output_unit)
            if (.not. streaming) exit
        end do
    end subroutine
end program
