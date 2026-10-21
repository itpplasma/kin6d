module kin6d_map_equilibrium
    ! Read-only N=0 cylindrical map -> S adapter. The reference form is retained
    ! with the map, rather than reconstructed from truncated Fourier B modes.
    ! Inverting the represented map is a geometric root solve; the field leaf
    ! and its JVP/VJP are generated from the same Piola formula. This entry
    ! deliberately rejects nonaxisymmetric/curve charts until D3/D7 admit them.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use kin6d_equilibrium_state, only: equilibrium_field_t, equilibrium_state_t, &
                   EQ_BR, EQ_BZ, EQ_BPHI, EQ_P, EQ_PSI, EQ_PSI_R, EQ_PSI_Z, EQ_VALUES, &
                                       EQ_UNAVAILABLE
    use kin6d_equilibrium_chart, only: cylindrical_chart_t
    use kin6d_ideal_map, only: ideal_map_t, ideal_frozen_t
    use kin6d_ideal_map_leaves, only: ideal_axisymmetric_field, ideal_reference_form
    use kin6d_grad_shafranov, only: gs_state_t
    use kin6d_gs_equilibrium, only: gs_geometry_view
    use kin6d_triangle_basis, only: triangle_basis_jet, triangle_bernstein_controls
    use kin6d_curved_triangle, only: triangle_geometry, triangle_locate
    use delaunay_types, only: mesh_t
    use point_location, only: locate_point
    implicit none
    private
    public :: map_field_t, map_equilibrium

    type, extends(equilibrium_field_t) :: map_field_t
        private
        type(gs_state_t) :: geometry
        real(dp), allocatable :: reference(:, :), omega(:), boxes(:, :, :)
        class(ideal_frozen_t), allocatable :: frozen
        type(mesh_t) :: locator
        real(dp) :: psi_bound = 0
    contains
        procedure :: sample => map_field_sample
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
        real(dp) :: controls(10), point(2), jac(2, 2), curvature(2, 2, 2)
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
            field%geometry%area(e) = det/2
        end do
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
        real(dp) :: xi(2), margin
        integer :: i, e, cell, location, seed, local
        values = 0; statuses = 1
        if (size(points, 1) < 2 .or. size(points, 1) > 3) return
        if (size(values, 1) /= EQ_VALUES) return
        if (size(values, 2) /= size(points, 2)) return
        if (size(statuses) /= size(points, 2)) return
        margin = 64*epsilon(1.0_dp)*maxval(abs(field%geometry%points))
        seed = 1
        !$omp parallel do default(shared) schedule(static) firstprivate(seed) &
        !$omp& private(i, e, cell, location, xi, local)
        do i = 1, size(points, 2)
            if (.not. all(ieee_is_finite(points(:, i)))) cycle
            if (points(1, i) <= 0) cycle
            call locate_point(field%locator, points(1, i), points(2, i), &
                              cell, location, seed)
            local = 1
            if (cell > 0) then
                call triangle_locate(field%geometry%points(:, &
                           field%geometry%triangles(:, cell)), points(:2, i), xi, local)
            end if
            ! The chord locator can miss the strip next to a curved face.
            ! Bernstein hulls restrict this fallback to possible cells.
            if (local /= 0) then
                do e = 1, field%locator%ntriangles
                    if (any(points(:2, i) < field%boxes(:, 1, e) - margin)) cycle
                    if (any(points(:2, i) > field%boxes(:, 2, e) + margin)) cycle
                    call triangle_locate(field%geometry%points(:, &
                              field%geometry%triangles(:, e)), points(:2, i), xi, local)
                    if (local == 0) then
                        cell = e
                        exit
                    end if
                end do
            end if
            if (local /= 0) cycle
            seed = cell
            call map_cell_sample(field, cell, xi, values(:, i), statuses(i))
        end do
        !$omp end parallel do
    end subroutine

    subroutine map_cell_sample(field, e, xi, value, status)
        class(map_field_t), intent(in) :: field
        integer, intent(in) :: e
        real(dp), intent(in) :: xi(2)
        real(dp), intent(out) :: value(EQ_VALUES)
        integer, intent(out) :: status
        real(dp) :: shape(10), grad(2, 10), reference(2), dg(2, 2), a(2, 2)
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
end module kin6d_map_equilibrium
