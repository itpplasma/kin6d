program coil_box_vacuum
    !! Vacuum field of MAKEGRID coils in a domain bounded by a star-shaped
    !! Fourier wall with P field periods, from the tangential coil vector
    !! potential on the wall (mapped NN1), checked against direct
    !! Biot-Savart in the L2 norm over the domain.
    !! Usage: coil_box_vacuum input.nml   (namelist /box/)
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use coil_box_geometry
    use milca_makegrid_circuits, only: makegrid_circuits_read
    use milca_coils_3d, only: coil_mu0
    use milca_coil_quadrature_field, only: coil_quadrature_field_t, coil_quadrature_init
    use milca_island_diagnostics, only: island_trace_options_t, island_trace
    use milca_mhd_mesh, only: flux_mesh_t, flux_mesh_build
    use milca_equilibrium_nonnested, only: free_b_field_t, nn_build
    use milca_nonnested_mapped, only: nn_map_t, nn_map_build, nn_mapped_report_t, &
        nn_mapped_cartesian_data, nn_mapped_solve, nn_mapped_cartesian_error, &
        nn_map_point
    use fortfem_meridian_feec_n, only: meridian_qgrad
    implicit none
    character(len=512) :: coils_file = '', wall_file = '', arg
    real(dp) :: extcur(64) = 0, tolerance = 1.0e-10_dp
    integer :: period = 1, nring = 8, ncore = 8, nthmax = 256, order = 2, modes = 4
    integer :: map_modes = 8, nphi = 0, nphi_error = 0
    integer :: trace_lines = 0, trace_periods = 100
    real(dp) :: trace_u0 = 0.5_dp, trace_u1 = 0.95_dp
    character(len=512) :: trace_file = 'poincare.txt'
    ! A test gate: stop with an error if the L2 field error exceeds it.
    real(dp) :: max_error = huge(1.0_dp)
    namelist /box/ coils_file, wall_file, extcur, period, nring, ncore, nthmax, &
        order, modes, map_modes, nphi, nphi_error, tolerance, trace_lines, &
        trace_periods, trace_u0, trace_u1, trace_file, max_error
    type(flux_mesh_t) :: fm
    type(free_b_field_t), target :: field
    type(nn_map_t), target :: map
    type(nn_mapped_report_t) :: report
    integer, allocatable :: circuit(:)
    real(dp), allocatable :: turns(:), data(:), rings(:)
    real(dp) :: err, eh, ee, fvac, b(3), p(2), scale, quality
    integer :: unit, st, ncircuit, j, c0, c1, rate

    call get_command_argument(1, arg)
    open (newunit=unit, file=trim(arg), status='old', action='read')
    read (unit, nml=box)
    close (unit)
    if (nphi == 0) nphi = 4*modes + 4
    if (nphi_error == 0) nphi_error = 2*nphi + 1
    call system_clock(c0, rate)

    call makegrid_circuits_read(trim(coils_file), coils, circuit, turns, ncircuit, st)
    if (st /= 0) error stop 'coils'
    allocate (weights(size(circuit)))
    weights = 0
    where (circuit >= 1) weights = turns*extcur(max(circuit, 1))
    call wall_read(trim(wall_file), period, st)
    if (st /= 0) error stop 'wall'
    quality = star_quality(4*map_modes + 4)
    write (*, '(a,es10.3)') 'wall star quality (min sin angle)= ', quality
    if (quality <= 0) error stop 'wall not star-shaped about its centroids'

    allocate (rings(nring))
    rings = [(real(j, dp)/nring, j=1, nring)]
    call flux_mesh_build(rings, ncore, nthmax, 2, ref_polar, fm, st)
    if (st /= 0) error stop 'mesh'
    ! F = R B_phi on the mean wall at phi = 0 (vacuum: constant on any loop
    ! linking the axis once, up to the 3D ripple), the reference harmonic.
    fvac = 0
    do j = 0, 63
        call ref_polar(1.0_dp, 2*acos(-1.0_dp)*j/64, p)
        call box_field([p(1), 0.0_dp, p(2)], b)
        fvac = fvac + p(1)*b(2)/64
    end do
    scale = 0
    do j = 0, 63
        call ref_polar(1.0_dp, 2*acos(-1.0_dp)*j/64, p)
        scale = scale + norm2(p - ref_centre())/64
    end do
    call nn_build(fm, order, modes, ref_centre(), scale, fvac, coil_mu0, field, st)
    if (st /= 0) error stop 'field'
    field%period = period
    field%coordinates => box_coordinates
    call set_map_planes(2*map_modes + 2)
    call nn_map_build(field%fs, map_modes, 2*map_modes + 2, box_map, map, st, &
        period=period)
    if (st /= 0) error stop 'map'
    call map_jacobian_range(field, map, nphi)
    call nn_mapped_cartesian_data(field, map, nphi, box_potential, data, st)
    if (st /= 0) error stop 'data (map Jacobian not positive?)'
    call nn_mapped_solve(field, map, nphi, data, tolerance, report, st)
    if (st /= 0) error stop 'solve'
    call nn_mapped_cartesian_error(field, map, nphi_error, box_field, err, eh, ee, st)
    if (st /= 0) error stop 'error norm'
    call system_clock(c1)
    write (*, '(a,i0,a,i0,a,i0,a,i0,a,i0)') 'period=', period, ' rings=', nring, &
        ' order=', order, ' modes=', modes, ' ndof=', report%ndof
    write (*, '(a,es10.3,a,es10.3,a,i0,a,es9.2,a,es9.2,a,f9.2)') &
        'field error / |B| (L2)= ', err/sqrt(2*coil_mu0*ee), &
        ' energy rel= ', abs(eh - ee)/ee, ' fgmres= ', report%iterations, &
        ' residual= ', report%residual, ' div= ', report%divergence, &
        ' s= ', real(c1 - c0, dp)/rate
    if (err/sqrt(2*coil_mu0*ee) > max_error) error stop 'field error above max_error'
    if (trace_lines > 0) call poincare(field)
contains
    subroutine poincare(field)
        ! Lines from the outboard ray theta = 0 of the reference section,
        ! traced period by period in the box field (reference coordinates,
        ! written mapped to physical) and from the same physical start in the
        ! direct Biot-Savart field. Section phi = 0 mod 2 pi / P.
        type(free_b_field_t), intent(in) :: field
        type(coil_quadrature_field_t) :: cq
        type(island_trace_options_t) :: opt
        real(dp) :: x(2), y(2), xb(2), u, fin(2)
        integer :: unit, i, k, st, sb
        call coil_quadrature_init(cq, coils, weights, 4, 0.0_dp, st)
        if (st /= 0) error stop 'coil quadrature field'
        open (newunit=unit, file=trim(trace_file), status='replace', action='write')
        write (unit, '(a)') '# line period R_box Z_box R_biot_savart Z_biot_savart (m), phi = 0'
        !$omp parallel do default(shared) private(i, k, u, x, y, xb, fin, st, sb) &
        !$omp schedule(dynamic, 1)
        do i = 1, trace_lines
            u = trace_u0 + (trace_u1 - trace_u0)*(i - 1)/max(1, trace_lines - 1)
            call ref_polar(u, 0.0_dp, x)
            call box_map(x, 0.0_dp, xb)
            st = 0; sb = 0
            do k = 0, trace_periods
                y = -1
                if (st == 0) call box_map(x, 0.0_dp, y)
                !$omp critical
                write (unit, '(2i7,4es20.11)') i, k, merge(y, [-1.0_dp, -1.0_dp], st == 0), &
                    merge(xb, [-1.0_dp, -1.0_dp], sb == 0)
                !$omp end critical
                if (st == 0) then
                    call island_trace(field, x, 0.0_dp, 1.0_dp/period, opt, fin, st)
                    x = fin
                end if
                if (sb == 0) then
                    call island_trace(cq, xb, 0.0_dp, 1.0_dp/period, opt, fin, sb)
                    xb = fin
                end if
                if (st /= 0 .and. sb /= 0) exit
            end do
        end do
        !$omp end parallel do
        close (unit)
    end subroutine

    subroutine map_jacobian_range(field, map, nphi)
        ! Smallest and largest map Jacobian over the quadrature points of the
        ! nphi planes of one period; the data and solve need it positive.
        type(free_b_field_t), intent(in) :: field
        type(nn_map_t), intent(in) :: map
        integer, intent(in) :: nphi
        real(dp) :: x(2), w, val(field%fs%nl0), dval(2, field%fs%nl0), pos(3)
        real(dp) :: p1(3, 3), p2(3, 3), jac, jmin, jmax, xmin(3)
        integer :: t, q, j, st, bad
        jmin = huge(1.0_dp); jmax = -huge(1.0_dp); bad = 0
        do t = 1, field%fs%mm%topology%n_triangles
            do q = 1, field%fs%nq
                call meridian_qgrad(field%fs, t, q, x, w, val, dval, st)
                do j = 1, nphi
                    call nn_map_point(map, field%fs%sp(0)%dof(:, t), val, dval, x, &
                        2*acos(-1.0_dp)*(j - 1)/(nphi*field%period), pos, p1, p2, jac, st)
                    if (st /= 0) bad = bad + 1
                    if (jac < jmin) then
                        jmin = jac; xmin = [x, 2*acos(-1.0_dp)*(j - 1)/(nphi*field%period)]
                    end if
                    jmax = max(jmax, jac)
                end do
            end do
        end do
        write (*, '(a,2es11.3,a,i0,a,3f9.4)') 'map Jacobian min/max= ', jmin, jmax, &
            ' failed points= ', bad, ' min at (R,Z,phi)= ', xmin
    end subroutine
end program coil_box_vacuum
