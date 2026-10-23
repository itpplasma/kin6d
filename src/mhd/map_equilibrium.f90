module kin6d_map_equilibrium
    ! Read-only N=0 cylindrical map -> S adapter. The reference form is retained
    ! with the map, rather than reconstructed from truncated Fourier B modes.
    ! Inverting the represented map is a geometric root solve; the field leaf
    ! and its JVP/VJP are generated from the same Piola formula. This entry
    ! deliberately rejects nonaxisymmetric/curve charts until D3/D7 admit them. D7 label contours use the same Pk
    ! projection used by the GS geometry services. This is a numerical readback;
    ! the projection and quadrature errors are not enclosed.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use kin6d_equilibrium_state, only: flux_labelled_field_t, equilibrium_state_t, &
                   EQ_BR, EQ_BZ, EQ_BPHI, EQ_P, EQ_PSI, EQ_PSI_R, EQ_PSI_Z, EQ_VALUES, &
                                       EQ_UNAVAILABLE
    use kin6d_equilibrium_chart, only: cylindrical_chart_t
    use kin6d_ideal_map, only: ideal_map_t, ideal_frozen_t
    use kin6d_ideal_map_leaves, only: ideal_axisymmetric_field, ideal_reference_form, &
        ideal_region_weight, ideal_contour_density
    use kin6d_grad_shafranov, only: gs_state_t, gs_axis_candidate, gs_element_jet
    use kin6d_gs_coarea, only: gs_level_coarea
    use kin6d_gs_axis_recovery, only: gs_magnetic_axis
    use fortfem_triangle_duffy_quadrature, only: triangle_duffy_quadrature
    use kin6d_gs_equilibrium, only: gs_geometry_view, gs_readback_resolution
    use kin6d_triangle_basis, only: triangle_basis_jet, triangle_bernstein_controls, &
        TRIANGLE_MAX_NODES
    use kin6d_curved_triangle, only: triangle_geometry, triangle_locate
    use delaunay_types, only: mesh_t
    use point_location, only: locate_point
    implicit none
    private
    public :: map_field_t, map_equilibrium

    type, extends(flux_labelled_field_t) :: map_field_t
        private
        type(gs_state_t) :: geometry
        real(dp), allocatable :: reference(:, :), omega(:), boxes(:, :, :)
        class(ideal_frozen_t), allocatable :: frozen
        type(mesh_t) :: locator
        real(dp) :: psi_bound = 0, orientation = 1, margin = 0
    contains
        procedure :: sample => map_field_sample
        procedure :: label_axis => map_label_axis
        procedure :: magnetic_axis => map_magnetic_axis
        procedure :: label_coarea => map_label_coarea
        procedure :: region_integrals => map_region_integrals
        procedure :: resolution => map_resolution
        procedure :: toroidal_density => map_toroidal_density
    end type
contains
    subroutine map_equilibrium(map, frozen, mesh, eq, status, f_vac, psi_bound)
        type(ideal_map_t), intent(in) :: map
        class(ideal_frozen_t), intent(in) :: frozen
        type(gs_state_t), intent(in) :: mesh
        type(equilibrium_state_t), intent(out) :: eq
        integer, intent(out) :: status
        real(dp), intent(in) :: f_vac
        real(dp), intent(in), optional :: psi_bound
        type(map_field_t), allocatable :: field
        real(dp) :: controls(TRIANGLE_MAX_NODES), point(2), jac(2, 2)
        real(dp) :: curvature(2, 2, 2)
        real(dp) :: det, inverse(2, 2), chi, grad(2), b3, pressure
        real(dp) :: form(3), dj(2, 2)
        integer :: e, a, j, g, local
        status = EQ_UNAVAILABLE
        if (map%ntor /= 0 .or. abs(map%sigma + 1) > epsilon(1.0_dp)) return
        if (.not. allocated(map%chart)) return
        select type (chart => map%chart)
        type is (cylindrical_chart_t)
        class default
            return
        end select
        status = 1
        if (.not. mesh%ready .or. map%mu0 <= 0) return
        if (.not. allocated(map%p)) return
        if (.not. allocated(map%chi_node)) return
        if (.not. all(ieee_is_finite([f_vac, map%mu0]))) return
        if (map%nodes /= size(mesh%points, 2)) return
        if (map%cells /= size(mesh%triangles, 2)) return
        if (map%nloc /= size(mesh%triangles, 1)) return
        if (any(map%conn /= mesh%triangles)) return
        if (.not. all(ieee_is_finite(map%x))) return
        allocate (field)
        field%psi_bound = sum(map%chi_node, mask=map%boundary)/count(map%boundary)
        if (present(psi_bound)) then
            if (.not. ieee_is_finite(psi_bound)) return
            field%psi_bound = psi_bound
        end if
        if (maxval(abs(pack(map%chi_node, map%boundary) - field%psi_bound)) &
            > 1.0e-12_dp*max(1.0_dp, abs(field%psi_bound))) return
        call gs_geometry_view(mesh, 0.0_dp, field%geometry)
        field%geometry%points = map%x(:2, 1, :)
        field%margin = 64*epsilon(1.0_dp)*maxval(abs(field%geometry%points))
        field%geometry%mapped = .true.
        field%geometry%curved_cell = .true.
        field%reference = map%reference
        field%omega = map%x(3, 1, :)
        allocate (field%frozen, source=frozen)
        do a = 1, map%nodes
            call frozen%evaluate(map%reference(:, a), chi, grad, b3, pressure)
            if (abs(chi - map%chi_node(a)) > &
                1.0e-12_dp*max(1.0_dp, abs(chi))) return
            field%geometry%psi(a) = chi
        end do
        ! Check that this evaluator is the one used to freeze the solver.
        do e = 1, map%cells
            do g = 1, map%nq
                point = matmul(map%reference(:, map%conn(:, e)), map%shape(:, g))
                call frozen%evaluate(point, chi, grad, b3, pressure)
                if (abs(pressure - map%p(g, e)) > &
                    1.0e-12_dp*max(1.0_dp, abs(pressure))) return
                dj = matmul(map%reference(:, map%conn(:, e)), &
                            transpose(map%dshape(:, :, g)))
                call ideal_reference_form(dj(1, 1), dj(2, 1), dj(1, 2), dj(2, 2), &
                                        grad(1), grad(2), b3, form(1), form(2), form(3))
                if (maxval(abs(form - map%b(:, g, e))) > &
                    1.0e-12_dp*max(1.0_dp, maxval(abs(form)))) return
                jac = matmul(field%geometry%points(:, map%conn(:, e)), &
                             transpose(map%dshape(:, :, g)))
                det = jac(1, 1)*jac(2, 2) - jac(1, 2)*jac(2, 1)
                if (det <= 0 .or. .not. ieee_is_finite(det)) return
            end do
        end do
        allocate (field%boxes(2, 2, map%cells))
        field%locator%npoints = map%nodes
        field%locator%ntriangles = map%cells
        allocate (field%locator%points(map%nodes), &
                  field%locator%triangles(map%cells))
        field%locator%points%x = field%geometry%points(1, :)
        field%locator%points%y = field%geometry%points(2, :)
        do e = 1, map%cells
            field%locator%triangles(e)%vertices = map%conn(:3, e)
            ! Triangle compatibility uses -1 on the exterior; the walking
            ! locator uses zero (triangle_compat.f90:180, point_location.f90:302).
            field%locator%triangles(e)%neighbors = max(0, mesh%neighbors(:, e))
            do j = 1, 2
                call triangle_bernstein_controls(field%geometry%points(j, &
                                                   map%conn(:, e)), controls(:map%nloc))
                field%boxes(j, 1, e) = minval(controls(:map%nloc))
                field%boxes(j, 2, e) = maxval(controls(:map%nloc))
            end do
            call triangle_geometry(field%geometry%points(:, map%conn(:, e)), &
                       [1.0_dp/3, 1.0_dp/3], point, jac, curvature, det, inverse, local)
            if (local /= 0 .or. det <= 0) return
            call triangle_geometry(field%geometry%points(:, map%conn(:3, e)), &
                [1.0_dp/3, 1.0_dp/3], point, jac, curvature, det, inverse, local)
            if (local /= 0) return
            field%geometry%gradient(:, 2, e) = inverse(1, :)
            field%geometry%gradient(:, 3, e) = inverse(2, :)
            field%geometry%gradient(:, 1, e) = -sum(transpose(inverse), dim=2)
            field%geometry%area(e) = det/2
        end do
        field%geometry%psi = field%geometry%psi - field%psi_bound
        if (maxval(field%geometry%psi) > -minval(field%geometry%psi)) then
            field%orientation = -1
            field%geometry%psi = -field%geometry%psi
        end if
        eq%mu0 = map%mu0
        eq%psi_bound = field%psi_bound
        eq%f_vac = f_vac
        eq%lane = 'map'
        eq%claim = 'candidate'
        eq%seed = 'frozen reference form and represented N=0 map'
        allocate (eq%chart, source=map%chart)
        call move_alloc(field, eq%field)
        status = 0
    end subroutine

    subroutine map_field_sample(field, points, values, statuses)
        class(map_field_t), intent(in) :: field
        real(dp), intent(in) :: points(:, :)
        real(dp), intent(out) :: values(:, :)
        integer, intent(out) :: statuses(:)
        real(dp) :: xi(2)
        integer :: i, cell, seed, local
        values = 0; statuses = 1
        if (size(points, 1) < 2 .or. size(points, 1) > 3) return
        if (size(values, 1) /= EQ_VALUES) return
        if (size(values, 2) /= size(points, 2)) return
        if (size(statuses) /= size(points, 2)) return
        seed = 1
        !$omp parallel do default(shared) schedule(static) firstprivate(seed) &
        !$omp& private(i, cell, xi, local)
        do i = 1, size(points, 2)
            if (.not. all(ieee_is_finite(points(:, i)))) cycle
            if (points(1, i) <= 0) cycle
            call map_locate(field, points(:2, i), seed, cell, xi, local)
            if (local /= 0) cycle
            seed = cell
            call map_cell_sample(field, cell, xi, values(:, i), statuses(i))
        end do
        !$omp end parallel do
    end subroutine

    subroutine map_locate(field, point, seed, cell, xi, status)
        class(map_field_t), intent(in) :: field
        real(dp), intent(in) :: point(2)
        integer, intent(in) :: seed
        integer, intent(out) :: cell, status
        real(dp), intent(out) :: xi(2)
        real(dp) :: margin
        integer :: location, local, e
        margin = field%margin
        call locate_point(field%locator, point(1), point(2), &
                          cell, location, seed)
        local = 1
        if (cell > 0) then
            call triangle_locate(field%geometry%points(:, &
                       field%geometry%triangles(:, cell)), point, xi, local)
        end if
        ! The chord locator can miss the strip next to a curved face.
        ! Bernstein hulls restrict this fallback to possible cells.
        if (local /= 0) then
            do e = 1, field%locator%ntriangles
                if (any(point < field%boxes(:, 1, e) - margin)) cycle
                if (any(point > field%boxes(:, 2, e) + margin)) cycle
                call triangle_locate(field%geometry%points(:, &
                          field%geometry%triangles(:, e)), point, xi, local)
                if (local == 0) then
                    cell = e
                    exit
                end if
            end do
        end if
        status = local
    end subroutine

    subroutine map_cell_sample(field, e, xi, value, status)
        class(map_field_t), intent(in) :: field
        integer, intent(in) :: e
        real(dp), intent(in) :: xi(2)
        real(dp), intent(out) :: value(EQ_VALUES)
        integer, intent(out) :: status
        real(dp) :: shape(TRIANGLE_MAX_NODES), grad(2, TRIANGLE_MAX_NODES)
        real(dp) :: reference(2), dg(2, 2), a(2, 2)
        real(dp) :: omega(2), chi, cg(2), b3, pressure, r, jac
        integer :: n
        n = size(field%geometry%triangles, 1)
        call triangle_basis_jet(xi, shape(:n), grad(:, :n))
        reference = matmul(field%reference(:, field%geometry%triangles(:, e)), &
                           shape(:n))
        dg = matmul(field%reference(:, field%geometry%triangles(:, e)), &
                    transpose(grad(:, :n)))
        a = matmul(field%geometry%points(:, field%geometry%triangles(:, e)), &
                   transpose(grad(:, :n)))
        omega = matmul(grad(:, :n), field%omega(field%geometry%triangles(:, e)))
        r = dot_product(field%geometry%points(1, &
                                             field%geometry%triangles(:, e)), shape(:n))
        call field%frozen%evaluate(reference, chi, cg, b3, pressure)
        call ideal_axisymmetric_field(a(1, 1), a(2, 1), a(1, 2), a(2, 2), &
                           dg(1, 1), dg(2, 1), dg(1, 2), dg(2, 2), omega(1), omega(2), &
                      cg(1), cg(2), b3, r, value(EQ_BR), value(EQ_BZ), value(EQ_BPHI), &
                                      value(EQ_PSI_R), value(EQ_PSI_Z), jac)
        value(EQ_P) = pressure
        value(EQ_PSI) = chi
        status = 1
        if (jac <= 0) return
        if (.not. all(ieee_is_finite(value))) return
        status = 0
    end subroutine
    subroutine map_resolution(field, radial_points, angular_points)
        class(map_field_t), intent(in) :: field
        integer, intent(inout) :: radial_points, angular_points
        call gs_readback_resolution(field%geometry, radial_points, angular_points)
    end subroutine

    subroutine map_label_axis(field, axis, psi_axis, status)
        class(map_field_t), intent(in) :: field
        real(dp), intent(out) :: axis(2), psi_axis
        integer, intent(out) :: status
        real(dp) :: hessian(2, 2)
        integer :: element
        call gs_axis_candidate(field%geometry, axis, psi_axis, hessian, element, status)
        psi_axis = field%orientation*psi_axis + field%psi_bound
    end subroutine

    subroutine map_magnetic_axis(field, axis, status)
        class(map_field_t), intent(in) :: field
        real(dp), intent(out) :: axis(2)
        integer, intent(out) :: status
        call gs_magnetic_axis(field%geometry, axis, status)
    end subroutine

    subroutine map_label_coarea(field, levels, angular_points, coarea, status, work)
        class(map_field_t), intent(in) :: field
        real(dp), intent(in) :: levels(:)
        integer, intent(in) :: angular_points
        real(dp), intent(out) :: coarea(:)
        integer, intent(out) :: status
        integer, intent(out), optional :: work(:)
        real(dp) :: axis(2)
        call gs_level_coarea(field%geometry, &
            field%orientation*(levels - field%psi_bound), angular_points, &
            coarea, axis, status, work=work)
    end subroutine

    subroutine map_toroidal_density(field, levels, angular_points, density, &
                                    status, work)
        class(map_field_t), intent(in) :: field
        real(dp), intent(in) :: levels(:)
        integer, intent(in) :: angular_points
        real(dp), intent(out) :: density(:)
        integer, intent(out) :: status
        integer, intent(out), optional :: work(:)
        real(dp), allocatable :: coarea(:), pieces(:, :)
        integer, allocatable :: counts(:, :)
        real(dp) :: axis(2), psi_axis, value(EQ_VALUES), theta, upper, ray(2)
        integer :: i, j, local, failed, cell, used
        real(dp) :: xi(2)
        status = 1
        if (angular_points < 16) return
        if (size(levels) /= size(density)) return
        if (present(work)) then
            if (size(work) /= size(levels)) return
        end if
        allocate (coarea(size(levels)), pieces(angular_points, size(levels)), &
            counts(angular_points, size(levels)))
        call field%label_coarea(levels, angular_points, coarea, status, work)
        if (status /= 0) return
        call field%label_axis(axis, psi_axis, status)
        if (status /= 0) return
        call map_locate(field, axis, 1, cell, xi, status)
        if (status /= 0) return
        call map_cell_sample(field, cell, xi, value, status)
        if (status /= 0) return
        density(1) = axis(1)*value(EQ_BPHI)*coarea(1)
        upper = 2*maxval(abs(field%geometry%points - spread(axis, 2, &
            size(field%geometry%points, 2))))
        failed = 0
        !$omp parallel do collapse(2) default(shared) schedule(static) &
        !$omp& private(i, j, theta, ray, local, used) reduction(max:failed)
        do i = 2, size(levels)
            do j = 1, angular_points
                theta = 2*acos(-1.0_dp)*(j - 0.5_dp)/angular_points
                ray = [cos(theta), sin(theta)]
                call map_ray_density(field, axis, ray, upper, &
                    field%orientation*(levels(i) - field%psi_bound), &
                    pieces(j, i), used, local)
                counts(j, i) = used
                if (local /= 0) failed = 1
            end do
        end do
        !$omp end parallel do
        status = failed
        if (status /= 0) return
        do i = 2, size(levels)
            density(i) = 2*acos(-1.0_dp)*sum(pieces(:, i))/angular_points
            if (present(work)) work(i) = work(i) + sum(counts(:, i))
        end do
    end subroutine

    subroutine map_ray_density(field, axis, ray, upper, target, density, used, status)
        class(map_field_t), intent(in) :: field
        real(dp), intent(in) :: axis(2), ray(2), upper, target
        real(dp), intent(out) :: density
        integer, intent(out) :: used, status
        real(dp) :: lo, hi, rho, point(2), flux, grad(2), hessian(2, 2)
        real(dp) :: xi(2), value(EQ_VALUES), last(2), denom
        integer :: seed, cell, local, iteration
        logical :: found
        lo = 0; hi = upper; seed = 1; found = .false.; used = 0
        density = 0; status = 1
        last = axis
        do iteration = 1, 60
            rho = (lo + hi)/2
            point = axis + rho*ray
            call map_locate(field, point, seed, cell, xi, local)
            used = used + 1
            if (local /= 0) then
                hi = rho
            else
                seed = cell
                call gs_element_jet(field%geometry, cell, point(1), point(2), &
                    flux, grad, hessian, local)
                if (local /= 0) return
                if (flux > target) then
                    hi = rho
                else
                    lo = rho
                    last = point
                    found = .true.
                end if
            end if
            if (hi - lo < 2.0e-11_dp*max(1.0_dp, upper)) exit
        end do
        if (.not. found) return
        call map_locate(field, last, seed, cell, xi, local)
        if (local /= 0) return
        call gs_element_jet(field%geometry, cell, last(1), last(2), &
            flux, grad, hessian, local)
        if (local /= 0) return
        if (abs(flux - target) > 1.0e-7_dp*max(1.0_dp, abs(target))) return
        denom = dot_product(grad, ray)
        if (denom <= 0) return
        call map_cell_sample(field, cell, xi, value, local)
        if (local /= 0) return
        call ideal_contour_density(lo, value(EQ_BPHI), grad(1), grad(2), &
            ray(1), ray(2), density)
        status = 0
    end subroutine

    subroutine map_region_integrals(field, volume, phi_tor, status)
        class(map_field_t), intent(in) :: field
        real(dp), intent(out) :: volume, phi_tor
        integer, intent(out) :: status
        real(dp), allocatable :: x(:), y(:), weight(:)
        real(dp) :: point(2), jac(2, 2), curvature(2, 2, 2), det, inverse(2, 2)
        real(dp) :: values(EQ_VALUES), v, f
        integer :: e, g, local, failed
        volume = 0; phi_tor = 0
        call triangle_duffy_quadrature(field%geometry%element_order + 2, &
            x, y, weight, status)
        if (status /= 0) return
        failed = 0
        !$omp parallel do default(shared) schedule(static) &
        !$omp& private(g, point, jac, curvature, det, inverse, values, local, v, f) &
        !$omp& reduction(+:volume, phi_tor) reduction(max:failed)
        do e = 1, size(field%geometry%triangles, 2)
            do g = 1, size(weight)
                call triangle_geometry(field%geometry%points(:, &
                    field%geometry%triangles(:, e)), [x(g), y(g)], point, &
                    jac, curvature, det, inverse, local)
                if (local /= 0) then
                    failed = 1
                    cycle
                end if
                call map_cell_sample(field, e, [x(g), y(g)], values, local)
                if (local /= 0) then
                    failed = 1
                    cycle
                end if
                call ideal_region_weight(point(1), det, values(EQ_BPHI), v, f)
                volume = volume + weight(g)*v
                phi_tor = phi_tor + weight(g)*f
            end do
        end do
        !$omp end parallel do
        status = failed
    end subroutine
end module kin6d_map_equilibrium
