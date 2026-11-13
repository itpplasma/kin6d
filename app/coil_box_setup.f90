module coil_box_setup
    !! The coil box shared by the vacuum and finite-pressure applications:
    !! coils and wall from files, the reference mesh on the phi-mean wall
    !! section, the NN1 field with P field periods and the box map, the map
    !! Jacobian report and Poincare tracing against direct Biot-Savart.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use coil_box_geometry
    use milca_makegrid_circuits, only: makegrid_circuits_read
    use milca_coils_3d, only: coil_mu0
    use milca_coil_quadrature_field, only: coil_quadrature_field_t, coil_quadrature_init
    use milca_island_diagnostics, only: island_trace_options_t, island_trace
    use milca_mhd_mesh, only: flux_mesh_t, flux_mesh_build
    use milca_equilibrium_nonnested, only: free_b_field_t, nn_build
    use milca_nonnested_mapped, only: nn_map_t, nn_map_build, nn_map_point
    use fortfem_meridian_feec_n, only: meridian_qgrad
    implicit none
    private
    public :: box_build, poincare, map_jacobian_range
contains
    subroutine box_build(coils_file, wall_file, extcur, period, nring, ncore, &
                         nthmax, order, modes, map_modes, nphi, field, map)
        ! F = R B_phi of the coils on the mean wall at phi = 0 is the
        ! reference harmonic F/R e_phi; the coil potential carries the rest.
        character(len=*), intent(in) :: coils_file, wall_file
        real(dp), intent(in) :: extcur(:)
        integer, intent(in) :: period, nring, ncore, nthmax, order, modes, map_modes, nphi
        type(free_b_field_t), intent(out) :: field
        type(nn_map_t), intent(out) :: map
        type(flux_mesh_t) :: fm
        integer, allocatable :: circuit(:)
        real(dp), allocatable :: turns(:), rings(:)
        real(dp) :: fvac, b(3), p(2), scale, quality
        integer :: st, ncircuit, j
        call makegrid_circuits_read(coils_file, coils, circuit, turns, ncircuit, st)
        if (st /= 0) error stop 'coils'
        allocate (weights(size(circuit)))
        weights = 0
        where (circuit >= 1) weights = turns*extcur(max(circuit, 1))
        call wall_read(wall_file, period, st)
        if (st /= 0) error stop 'wall'
        quality = star_quality(4*map_modes + 4)
        write (*, '(a,es10.3)') 'wall star quality (min sin angle)= ', quality
        if (quality <= 0) error stop 'wall not star-shaped about its centroids'
        allocate (rings(nring))
        rings = [(real(j, dp)/nring, j=1, nring)]
        call flux_mesh_build(rings, ncore, nthmax, 2, ref_polar, fm, st)
        if (st /= 0) error stop 'mesh'
        fvac = 0; scale = 0
        do j = 0, 63
            call ref_polar(1.0_dp, 2*acos(-1.0_dp)*j/64, p)
            call box_field([p(1), 0.0_dp, p(2)], b)
            fvac = fvac + p(1)*b(2)/64
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
    end subroutine

    subroutine poincare(field, period, trace_lines, trace_periods, trace_u0, &
                        trace_u1, trace_file)
        ! Lines from the outboard ray theta = 0 of the reference section,
        ! traced period by period in the box field (reference coordinates,
        ! written mapped to physical) and from the same physical start in the
        ! direct Biot-Savart field. Section phi = 0 mod 2 pi / P.
        type(free_b_field_t), intent(in) :: field
        integer, intent(in) :: period, trace_lines, trace_periods
        real(dp), intent(in) :: trace_u0, trace_u1
        character(len=*), intent(in) :: trace_file
        type(coil_quadrature_field_t) :: cq
        type(island_trace_options_t) :: opt
        ! Section points of each line in the box (kind 1) and in the direct
        ! field (kind 2), -1 after a failed trace, and the period of failure.
        real(dp), allocatable :: path(:, :, :, :)
        integer, allocatable :: last(:, :)
        real(dp) :: x(2), y(2), u, fin(2)
        integer :: unit, i, k, st, task, kind
        call coil_quadrature_init(cq, coils, weights, 4, 0.0_dp, st)
        if (st /= 0) error stop 'coil quadrature field'
        allocate (path(2, 0:trace_periods, 2, trace_lines), last(2, trace_lines))
        path = -1
        last = trace_periods
        ! The box line and the direct line are separate tasks, so that twice
        ! as many threads find work as there are lines.
        !$omp parallel do default(shared) private(task, i, kind, k, u, x, y, fin, st) &
        !$omp schedule(dynamic, 1)
        do task = 1, 2*trace_lines
            i = (task + 1)/2; kind = 2 - mod(task, 2)
            u = trace_u0 + (trace_u1 - trace_u0)*(i - 1)/max(1, trace_lines - 1)
            call ref_polar(u, 0.0_dp, x)
            if (kind == 2) then
                call box_map(x, 0.0_dp, y)
                x = y
            end if
            do k = 0, trace_periods
                if (kind == 1) then
                    call box_map(x, 0.0_dp, y)
                    path(:, k, kind, i) = y
                    call island_trace(field, x, 0.0_dp, 1.0_dp/period, opt, fin, st)
                else
                    path(:, k, kind, i) = x
                    call island_trace(cq, x, 0.0_dp, 1.0_dp/period, opt, fin, st)
                end if
                if (st /= 0) then
                    last(kind, i) = k
                    exit
                end if
                x = fin
            end do
        end do
        !$omp end parallel do
        open (newunit=unit, file=trim(trace_file), status='replace', action='write')
        write (unit, '(a)') '# line period R_box Z_box R_biot_savart Z_biot_savart (m), phi = 0'
        do i = 1, trace_lines
            ! A period is written while either line is still alive.
            do k = 0, maxval(last(:, i))
                write (unit, '(2i7,4es20.11)') i, k, path(:, k, 1, i), path(:, k, 2, i)
            end do
        end do
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
end module coil_box_setup
