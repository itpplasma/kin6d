module kin6d_gs_equilibrium
    ! psi lane -> core state adapter (docs/equilibrium-core.md, M1). The field
    ! of a solved gs_state_t is B = grad psi x grad phi + F(psi) grad phi with
    ! the profile law (or the linear p', FF' law) of the solve. The adapter
    ! keeps a read-only geometry view of the mesh (no solve operator) whose
    ! flux is u = psi - psi_bound, oriented so that the magnetic axis is a
    ! minimum; the sign is undone on every output. psi_bound is the flux of
    ! the plasma boundary in the lane gauge: 0 for fixed boundary (psi = 0 on
    ! the mesh boundary), the LCFS flux for free boundary (psi = 0 on the
    ! symmetry axis). plasma_cell restricts the view to the plasma cells of a
    ! free-boundary mesh (sub-mesh, compactly renumbered).
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use kin6d_grad_shafranov, only: gs_state_t, gs_profile_integral, &
        gs_axis_candidate
    use kin6d_gs_field_readback, only: sample_gs_many
    use kin6d_gs_coarea, only: gs_level_coarea
    use kin6d_gs_axis_recovery, only: gs_magnetic_axis
    use kin6d_curved_triangle, only: triangle_geometry
    use kin6d_triangle_basis, only: triangle_basis_jet
    use fortfem_triangle_duffy_quadrature, only: triangle_duffy_quadrature
    use kin6d_gs_profile_law, only: gs_profile_law_t, gs_law_primitives
    use kin6d_equilibrium_state, only: labelled_field_t, equilibrium_state_t, &
        EQ_BR, EQ_BZ, EQ_BPHI, EQ_P, EQ_PSI, EQ_PSI_R, EQ_PSI_Z, EQ_VALUES, &
        EQ_GAUGE_AXIS, EQ_GAUGE_BOUNDARY, EQ_UNAVAILABLE
    use kin6d_equilibrium_chart, only: cylindrical_chart_t
    implicit none
    private
    public :: gs_field_t, gs_equilibrium, gs_profile_f, gs_readback_resolution

    type, extends(labelled_field_t) :: gs_field_t
        type(gs_state_t) :: view
        real(dp) :: orientation = 1, psi_bound = 0
        real(dp) :: f_edge = 0, ffprime = 0, pprime = 0
        ! The view is the whole mesh: its boundary is the plasma boundary.
        logical :: whole = .true.
        type(gs_profile_law_t), allocatable :: law
    contains
        procedure :: sample => gs_field_sample
        procedure :: profiles => gs_field_profiles
        procedure :: label_axis => gs_field_label_axis
        procedure :: magnetic_axis => gs_field_magnetic_axis
        procedure :: label_coarea => gs_field_label_coarea
        procedure :: region_integrals => gs_field_region_integrals
        procedure :: resolution => gs_field_resolution
    end type
contains
    subroutine gs_readback_resolution(state, radial_points, angular_points)
        type(gs_state_t), intent(in) :: state
        integer, intent(inout) :: radial_points, angular_points
        ! Resolve radial changes on the element scale. Cubic profile
        ! interpolation/integration is at least the P2/P3 field order;
        ! each smooth cell arc also shrinks with the mesh.
        if (radial_points == 0) radial_points = max(17, &
            ceiling(sqrt(real(size(state%triangles, 2), dp))) + 1)
        if (angular_points == 0) angular_points = 16*(state%element_order + 2)
    end subroutine

    subroutine gs_field_resolution(field, radial_points, angular_points)
        class(gs_field_t), intent(in) :: field
        integer, intent(inout) :: radial_points, angular_points
        call gs_readback_resolution(field%view, radial_points, angular_points)
    end subroutine

    subroutine gs_field_label_axis(field, axis, psi_axis, status)
        ! Represented flux extremum (the contour and axis-q reference).
        class(gs_field_t), intent(in) :: field
        real(dp), intent(out) :: axis(2), psi_axis
        integer, intent(out) :: status
        real(dp) :: hessian(2, 2)
        integer :: element
        call gs_axis_candidate(field%view, axis, psi_axis, hessian, element, status)
        psi_axis = field%orientation*psi_axis + field%psi_bound
    end subroutine

    subroutine gs_field_magnetic_axis(field, axis, status)
        class(gs_field_t), intent(in) :: field
        real(dp), intent(out) :: axis(2)
        integer, intent(out) :: status
        call gs_magnetic_axis(field%view, axis, status)
    end subroutine

    subroutine gs_field_label_coarea(field, levels, angular_points, coarea, status, &
                                     work)
        class(gs_field_t), intent(in) :: field
        real(dp), intent(in) :: levels(:)
        integer, intent(in) :: angular_points
        real(dp), intent(out) :: coarea(:)
        integer, intent(out) :: status
        integer, intent(out), optional :: work(:)
        real(dp) :: axis(2)
        call gs_level_coarea(field%view, field%orientation*(levels - field%psi_bound), &
            angular_points, coarea, axis, status, work=work)
    end subroutine

    subroutine gs_field_region_integrals(field, volume, phi_tor, status)
        ! Whole-mesh views only: a plasma sub-mesh holds cut cells whose
        ! plasma part needs the solver's level-cut quadrature.
        class(gs_field_t), intent(in) :: field
        real(dp), intent(out) :: volume, phi_tor
        integer, intent(out) :: status
        real(dp), allocatable :: x(:), y(:), w(:)
        real(dp) :: point(2), jac(2, 2), curvature(2, 2, 2), det, inverse(2, 2)
        real(dp) :: shape(10), gradients(2, 10), hessians(2, 2, 10)
        real(dp) :: flux, f, pressure, two_pi, volume_sum, flux_sum
        real(dp), allocatable :: shapes(:, :)
        integer :: e, node, ndof, local_status, failed_element, failed_status
        integer :: snapshot
        volume = 0; phi_tor = 0
        status = EQ_UNAVAILABLE
        if (.not. field%whole) return
        ndof = size(field%view%triangles, 1)
        two_pi = 2*acos(-1.0_dp)
        call triangle_duffy_quadrature(field%view%element_order + 2, x, y, w, status)
        if (status /= 0) return
        allocate (shapes(ndof, size(w)))
        do node = 1, size(w)
            call triangle_basis_jet([x(node), y(node)], shape(:ndof), &
                gradients(:, :ndof), hessians(:, :, :ndof))
            shapes(:, node) = shape(:ndof)
        end do
        ! Element-parallel volume and toroidal-flux sums: static partition per
        ! thread, reduction over the thread partial sums.
        volume_sum = 0; flux_sum = 0
        failed_element = huge(1); failed_status = 1
        !$omp parallel do default(shared) schedule(static) &
        !$omp& private(node, point, jac, curvature, det, inverse, flux, f, pressure) &
        !$omp& private(local_status, snapshot) &
        !$omp& reduction(+: volume_sum, flux_sum)
        do e = 1, size(field%view%triangles, 2)
            !$omp atomic read
            snapshot = failed_element
            if (e > snapshot) cycle
            do node = 1, size(w)
                if (field%view%curved_cell(e)) then
                    call triangle_geometry(field%view%points(:, &
                        field%view%triangles(:, e)), [x(node), y(node)], point, &
                        jac, curvature, det, inverse, local_status)
                    if (local_status /= 0) then
                        call record_failure(e, local_status)
                        exit
                    end if
                else
                    point = matmul(field%view%points(:, field%view%triangles(1:3, e)), &
                        [1 - x(node) - y(node), x(node), y(node)])
                    det = 2*field%view%area(e)
                end if
                flux = dot_product(shapes(:, node), &
                    field%view%psi(field%view%triangles(:, e)))
                call field%profiles(field%orientation*flux, f, pressure)
                if (abs(f) <= tiny(1.0_dp)) then
                    call record_failure(e, 1)
                    exit
                end if
                volume_sum = volume_sum + two_pi*w(node)*det*point(1)
                flux_sum = flux_sum + w(node)*det*f/point(1)
            end do
        end do
        !$omp end parallel do
        if (failed_element /= huge(1)) then
            status = failed_status
            return
        end if
        volume = volume_sum
        phi_tor = flux_sum
        status = 0
    contains
        subroutine record_failure(element, code)
            ! The lowest failing element decides the status, as in a serial pass.
            integer, intent(in) :: element, code
            !$omp critical(gs_equilibrium_failure)
            if (element < failed_element) then
                failed_element = element
                failed_status = code
            end if
            !$omp end critical(gs_equilibrium_failure)
        end subroutine
    end subroutine

    subroutine gs_equilibrium(state, eq, status, f_edge, ffprime, pprime, law, mu0, &
                              psi_bound, plasma_cell, gauge)
        ! Without psi_bound the state is a fixed-boundary solution and its
        ! boundary flux must be 0 (to 1e-12); with it, the boundary flux of
        ! the whole-mesh view must equal psi_bound. gauge defaults to the
        ! boundary gauge without psi_bound and to the axis gauge with it.
        type(gs_state_t), intent(in) :: state
        type(equilibrium_state_t), intent(out) :: eq
        integer, intent(out) :: status
        real(dp), intent(in) :: f_edge, ffprime
        real(dp), intent(in), optional :: pprime, mu0, psi_bound
        type(gs_profile_law_t), intent(in), optional :: law
        logical, intent(in), optional :: plasma_cell(:)
        integer, intent(in), optional :: gauge
        type(gs_field_t), allocatable :: field
        status = 1
        if (.not. state%ready .or. state%element_order < 2) return
        if (.not. all(ieee_is_finite([f_edge, ffprime]))) return
        allocate (field)
        field%f_edge = f_edge; field%ffprime = ffprime
        if (present(pprime)) field%pprime = pprime
        if (present(law)) field%law = law
        if (present(psi_bound)) then
            if (.not. ieee_is_finite(psi_bound)) return
            field%psi_bound = psi_bound
            eq%psi_gauge = EQ_GAUGE_AXIS
        else
            eq%psi_gauge = EQ_GAUGE_BOUNDARY
        end if
        if (present(gauge)) eq%psi_gauge = gauge
        if (present(mu0)) eq%mu0 = mu0
        if (present(plasma_cell)) then
            if (size(plasma_cell) /= size(state%triangles, 2)) return
            field%whole = all(plasma_cell)
        end if
        if (field%whole) then
            call geometry_view(state, field%psi_bound, field%view)
            if (any(abs(pack(field%view%psi, field%view%boundary)) > 1.0e-12_dp)) &
                return
        else
            call plasma_view(state, field%psi_bound, plasma_cell, field%view, status)
            if (status /= 0) return
            status = 1
        end if
        if (maxval(field%view%psi) > -minval(field%view%psi)) then
            field%orientation = -1
            field%view%psi = -field%view%psi
        end if
        eq%psi_bound = field%psi_bound
        eq%f_vac = f_edge
        eq%lane = 'psi'
        eq%claim = state%claim
        allocate (cylindrical_chart_t :: eq%chart)
        call move_alloc(field, eq%field)
        status = 0
    end subroutine

    subroutine gs_field_sample(field, points, values, statuses)
        class(gs_field_t), intent(in) :: field
        real(dp), intent(in) :: points(:, :)
        real(dp), intent(out) :: values(:, :)
        integer, intent(out) :: statuses(:)
        real(dp), allocatable :: raw(:, :)
        real(dp) :: u, f, pressure
        integer :: i
        values = 0
        if (size(values, 1) /= EQ_VALUES .or. size(values, 2) /= size(points, 2)) then
            statuses = 1
            return
        end if
        allocate (raw(4, size(points, 2)))
        call sample_gs_many(field%view, points, 1.0_dp, raw, statuses)
        do i = 1, size(points, 2)
            if (statuses(i) /= 0) cycle
            u = field%orientation*raw(1, i)
            ! Beyond the plasma boundary the profiles keep their edge values
            ! (vacuum F; the oriented flux there is positive).
            call field%profiles(field%orientation*min(raw(1, i), 0.0_dp), f, &
                pressure)
            values(EQ_BR, i) = field%orientation*raw(2, i)
            values(EQ_BZ, i) = field%orientation*raw(3, i)
            values(EQ_BPHI, i) = f/points(1, i)
            values(EQ_P, i) = pressure
            values(EQ_PSI, i) = u + field%psi_bound
            values(EQ_PSI_R, i) = points(1, i)*values(EQ_BZ, i)
            values(EQ_PSI_Z, i) = -points(1, i)*values(EQ_BR, i)
        end do
    end subroutine

    subroutine gs_field_profiles(field, u, f, pressure)
        ! F and p of the flux u = psi - psi_bound measured from the boundary.
        class(gs_field_t), intent(in) :: field
        real(dp), intent(in) :: u
        real(dp), intent(out) :: f, pressure
        integer :: status
        if (allocated(field%law)) then
            call gs_law_primitives(field%law, u, pressure, f, status)
            if (status /= 0) f = 0
            return
        end if
        f = gs_profile_f(u, field%f_edge, field%ffprime)
        pressure = field%pprime*u
    end subroutine

    real(dp) function gs_profile_f(psi, f_edge, ffprime, law) result(f)
        real(dp), intent(in) :: psi, f_edge, ffprime
        type(gs_profile_law_t), intent(in), optional :: law
        real(dp) :: squared, pressure
        integer :: status
        if (present(law)) then
            call gs_law_primitives(law, psi, pressure, f, status)
            if (status /= 0) f = 0
            return
        end if
        squared = f_edge**2 + 2*gs_profile_integral( &
            [-1.0_dp, 0.0_dp], [ffprime, ffprime], psi)
        f = 0
        if (squared <= 0) return
        f = sign(sqrt(squared), f_edge)
    end function

    subroutine geometry_view(state, psi_bound, view)
        ! Geometry and flux only; no stiffness, operator or factor.
        type(gs_state_t), intent(in) :: state
        real(dp), intent(in) :: psi_bound
        type(gs_state_t), intent(out) :: view
        view%points = state%points
        view%psi = state%psi - psi_bound
        view%triangles = state%triangles
        view%neighbors = state%neighbors
        view%boundary = state%boundary
        view%area = state%area
        view%curved_cell = state%curved_cell
        if (allocated(state%gradient)) view%gradient = state%gradient
        view%mapped = state%mapped
        view%element_order = state%element_order
        view%claim = state%claim
        view%ready = state%ready
    end subroutine

    subroutine plasma_view(state, psi_bound, plasma_cell, view, status)
        ! The sub-mesh of the flagged cells, compactly renumbered, with the
        ! flux measured from psi_bound. Nodes on faces without a flagged
        ! neighbour are its boundary.
        type(gs_state_t), intent(in) :: state
        real(dp), intent(in) :: psi_bound
        logical, intent(in) :: plasma_cell(:)
        type(gs_state_t), intent(out) :: view
        integer, intent(out) :: status
        integer, allocatable :: cells(:), cell_map(:), node_map(:), nodes(:)
        integer :: e, k, j, ne, nn, order, neighbor
        status = 1
        ne = size(state%triangles, 2)
        nn = size(state%points, 2)
        if (.not. any(plasma_cell)) return
        cells = pack([(e, e=1, ne)], plasma_cell)
        allocate (cell_map(ne), node_map(nn))
        cell_map = 0
        cell_map(cells) = [(e, e=1, size(cells))]
        node_map = 0
        do e = 1, size(cells)
            node_map(state%triangles(:, cells(e))) = 1
        end do
        nodes = pack([(j, j=1, nn)], node_map > 0)
        node_map(nodes) = [(j, j=1, size(nodes))]
        view%points = state%points(:, nodes)
        view%psi = state%psi(nodes) - psi_bound
        view%triangles = state%triangles(:, cells)
        do e = 1, size(cells)
            view%triangles(:, e) = node_map(view%triangles(:, e))
        end do
        view%neighbors = state%neighbors(:, cells)
        view%area = state%area(cells)
        view%curved_cell = state%curved_cell(cells)
        if (allocated(state%gradient)) view%gradient = state%gradient(:, :, cells)
        view%mapped = state%mapped
        view%element_order = state%element_order
        view%claim = state%claim
        order = state%element_order
        allocate (view%boundary(size(nodes)))
        view%boundary = .false.
        do e = 1, size(cells)
            do k = 1, 3
                neighbor = view%neighbors(k, e)
                if (neighbor > 0) neighbor = cell_map(neighbor)
                view%neighbors(k, e) = neighbor
                if (neighbor > 0) cycle
                view%boundary(view%triangles([k, mod(k, 3) + 1], e)) = .true.
                do j = 1, order - 1
                    view%boundary(view%triangles(3 + (k - 1)*(order - 1) + j, e)) &
                        = .true.
                end do
            end do
        end do
        view%ready = .true.
        status = 0
    end subroutine
end module kin6d_gs_equilibrium
