module kin6d_gs_nonlinear
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use kin6d_grad_shafranov, only: gs_state_t, solve_gs, apply_gs, &
                                    gs_cell_data, gs_axis_candidate, &
                                    gs_linear_response, gs_factor_ready, &
                                    gs_apply_inverse
    use kin6d_gs_csr, only: gs_csr_dirichlet_apply, gs_csr_fill_masked
    use kin6d_gs_options, only: gs_use_csr
    use kin6d_gs_profile_law, only: gs_profile_law_t, gs_law_sources, gs_law_kink, &
                                  gs_law_c1
    use kin6d_grad_shafranov_leaves, only: gs_source, gs_source_jvp, &
                                      gs_cut_discriminant, gs_cut_quadratic_discriminant
    use kin6d_triangle_basis, only: triangle_basis_jet, cubic_segment_roots, &
                                    triangle_bernstein_controls, &
                                    cubic_segment_stationary
    use kin6d_curved_triangle, only: triangle_locate
    use fortfem_triangle_level_quadrature, only: triangle_level_workspace_t, &
                       initialize_triangle_level_workspace, integrate_triangle_levels, &
                                                 triangle_level_coefficients
    use fortfem_triangle_duffy_quadrature, only: triangle_duffy_quadrature
    use fortnum_krylov, only: real_gmres_operator, real_conjugate_gradient_operator, &
                               KRYLOV_OK, KRYLOV_MAX_ITERATIONS
    use fortnum_linalg, only: det2
    use kin6d_gs_timers, only: gs_wall_s, gs_timer_add, gs_timer_iterations, &
                               gs_progress, gs_progress_enabled, GS_T_ASSEMBLY, &
                               GS_T_JAC_CG, GS_T_JAC_GMRES, GS_T_PRECOND_CG
    use fortnum_quadrature, only: gauss_legendre_ab
    implicit none
    private
    public :: gs_forcing_t, gs_normalize_current
    public :: solve_gs_cubic_profiles, gs_cubic_load, gs_profile_axis
    public :: gs_profile_jacobian_t, apply_gs_profile_jacobian, gs_profile_response
    type :: gs_profile_jacobian_t
        real(dp), allocatable :: mass(:, :, :), axis_load(:), axis_gradient(:)
        ! `mass` summed onto the CSR pattern of the stiffness (interior rows and
        ! columns), refilled with every Jacobian. Absent: element-wise products.
        real(dp), allocatable :: mval(:)
    end type
    type :: source_rule_t
        real(dp), allocatable :: weights(:), shapes(:, :), gradients(:, :, :)
    end type
    ! Shared read-only quadrature tables of gs_cubic_load.
    type :: load_tables_t
        real(dp), allocatable :: w(:), wh(:), shapes(:, :), grads(:, :, :)
        real(dp), allocatable :: shapes_h(:, :), grads_h(:, :, :)
    end type
    ! Per-thread scratch of gs_cubic_load: the FortFEM workspace and the cut
    ! event and value buffers.
    type :: load_scratch_t
        type(triangle_level_workspace_t) :: workspace
        type(source_rule_t) :: rules(4)
        real(dp), allocatable :: events(:), inner(:), epsabs(:), rtol(:), values(:)
        real(dp), allocatable :: errors(:), cell_values(:), high_values(:)
        logical, allocatable :: inner_active(:)
        logical :: ready = .false.
    end type
    abstract interface
        function gs_forcing_t(r, z) result(source)
            import :: dp
            real(dp), intent(in) :: r, z
            real(dp) :: source
        end function
    end interface
    interface
        subroutine dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
            import :: dp
            integer, intent(in) :: n, nrhs, lda, ldb
            integer, intent(out) :: ipiv(*), info
            real(dp), intent(inout) :: a(lda, *), b(ldb, *)
        end subroutine
        subroutine dgeev(jobvl, jobvr, n, a, lda, wr, wi, vl, ldvl, vr, ldvr, &
                         work, lwork, info)
            import :: dp
            character, intent(in) :: jobvl, jobvr
            integer, intent(in) :: n, lda, ldvl, ldvr, lwork
            real(dp), intent(inout) :: a(lda, *)
            real(dp), intent(out) :: wr(*), wi(*), vl(ldvl, *), vr(ldvr, *), work(*)
            integer, intent(out) :: info
        end subroutine
    end interface
contains
    subroutine gs_profile_axis(state, axis, gradient, status)
        type(gs_state_t), intent(inout) :: state
        real(dp), intent(out) :: axis, gradient(:)
        integer, intent(out) :: status
        real(dp) :: point(2), hessian(2, 2), local(2), shape(10)
        real(dp) :: g(2, 10), h(2, 2, 10), orientation
        integer :: cell, ndof
        gradient = 0
        orientation = 1
        if (maxval(state%psi) > -minval(state%psi)) orientation = -1
        state%psi = orientation*state%psi
        call gs_axis_candidate(state, point, axis, hessian, cell, status)
        state%psi = orientation*state%psi
        axis = orientation*axis
        if (status /= 0) return
        if (abs(axis) <= tiny(1.0_dp)) then
            status = 1
            return
        end if
        ndof = size(state%triangles, 1)
        call triangle_locate(state%points(:, state%triangles(:, cell)), &
                             point, local, status)
        if (status /= 0) return
        call triangle_basis_jet(local, shape(:ndof), g(:, :ndof), h(:, :, :ndof))
        ! At a stationary point the spatial displacement term vanishes.
        gradient(state%triangles(:, cell)) = shape(:ndof)
    end subroutine

    subroutine gs_cubic_load(state, mu0, law, rhs, jacobian, status, tolerance, &
                             forcing, jacobian_tolerance, load_only, &
                             quadrature_evaluations, source_integral, plasma_cell, &
                             cell_psi)
        !! Free-boundary mode (`plasma_cell` present): the flux handed in is
        !! measured from the plasma boundary flux, law%axis is the axis flux
        !! in the same measure, and the source is the law inside the cells
        !! flagged in `plasma_cell` where the normalised flux is below one
        !! and zero elsewhere. The cut at the boundary level (knot s = 1)
        !! uses the same exact level-cut integration as every other knot.
        !! No Dirichlet rows are applied and the axis rank-one gradient is
        !! left to the caller (axis_gradient is returned zero); axis_load is
        !! the derivative of the load with respect to law%axis.
        type(gs_state_t), intent(inout) :: state
        real(dp), intent(in) :: mu0
        type(gs_profile_law_t), intent(inout) :: law
        real(dp), intent(out) :: rhs(:)
        type(gs_profile_jacobian_t), intent(out) :: jacobian
        integer, intent(out) :: status
        real(dp), intent(in), optional :: tolerance
        procedure(gs_forcing_t), optional :: forcing
        real(dp), intent(out), optional :: source_integral
        logical, intent(in), optional :: plasma_cell(:)
        !> Free-boundary mode only: per-cell nodal flux (ndof, cells) that replaces
        !! state%psi in the level cuts (face-coincident cut handling).
        real(dp), intent(in), optional :: cell_psi(:, :)
        real(dp), intent(in), optional :: jacobian_tolerance
        logical, intent(in), optional :: load_only
        integer, intent(out), optional :: quadrature_evaluations
        real(dp), allocatable :: levels(:), x(:), y(:), w(:)
        real(dp), allocatable :: xh(:), yh(:), wh(:)
        real(dp), allocatable :: rhs_local(:, :), axis_local(:, :), source_local(:)
        integer, allocatable :: element_status(:), element_evaluations(:)
        type(load_tables_t) :: tables
        logical :: smooth_law, need_jacobian, free_mode, scratch_failed
        real(dp) :: jtol, tol, h(2, 2, 10)
        integer :: e, ndof, nc, i, j, k, q, nel, capacity, first_failure
        integer :: evaluations
        status = 1; rhs = 0
        if (present(source_integral)) source_integral = 0
        need_jacobian = .true.
        if (present(load_only)) need_jacobian = .not. load_only
        if (present(quadrature_evaluations)) quadrature_evaluations = 0
        evaluations = 0
        if (.not. state%ready .or. mu0 <= 0) return
        if (.not. ieee_is_finite(mu0)) return
        if (.not. law%valid()) return
        smooth_law = gs_law_c1(law)
        free_mode = present(plasma_cell)
        if (free_mode) then
            if (.not. law%normalized .or. size(plasma_cell) /= size(state%triangles, 2)) &
                return
            ! The source jumps at the boundary level: never the smooth rule.
            smooth_law = .false.
        end if
        if (.not. all(ieee_is_finite(state%psi))) return
        tol = 1.0e-9_dp
        if (present(tolerance)) tol = tolerance
        if (.not. ieee_is_finite(tol)) return
        if (tol <= 0 .or. tol >= 1) return
        jtol = tol
        if (present(jacobian_tolerance)) jtol = jacobian_tolerance
        if (.not. ieee_is_finite(jtol)) return
        if (jtol <= 0 .or. jtol >= 1) return
        ndof = size(state%triangles, 1); nc = ndof*(ndof + 2)
        if (.not. need_jacobian) nc = ndof
        allocate (jacobian%axis_gradient(size(rhs)))
        jacobian%axis_gradient = 0
        if (need_jacobian) then
            allocate (jacobian%mass(ndof, ndof, size(state%triangles, 2)), &
                      jacobian%axis_load(size(rhs)))
            jacobian%mass = 0; jacobian%axis_load = 0
        end if
        if (law%normalized .and. .not. free_mode) then
            call gs_profile_axis(state, law%axis, jacobian%axis_gradient, status)
            if (status /= 0) return
        end if
        levels = law%knots
        if (law%normalized) levels = law%axis*(1 - levels)
        if (levels(1) > levels(size(levels))) levels = levels(size(levels):1:-1)
        capacity = 12*size(levels) + 3
        if (state%element_order == 3 .and. .not. smooth_law) capacity = 27
        nel = size(state%triangles, 2)
        call triangle_duffy_quadrature(2*state%element_order + 2, x, y, w, status)
        if (status /= 0) return
        call triangle_duffy_quadrature(2*state%element_order + 4, xh, yh, wh, status)
        if (status /= 0) return
        allocate (tables%w(size(w)), tables%wh(size(wh)))
        tables%w = w; tables%wh = wh
        allocate (tables%shapes(ndof, size(w)), tables%grads(2, ndof, size(w)), &
                  tables%shapes_h(ndof, size(wh)), tables%grads_h(2, ndof, size(wh)))
        do q = 1, size(w)
            call triangle_basis_jet([x(q), y(q)], tables%shapes(:, q), &
                                    tables%grads(:, :, q), h(:, :, :ndof))
        end do
        do q = 1, size(wh)
            call triangle_basis_jet([xh(q), yh(q)], tables%shapes_h(:, q), &
                                    tables%grads_h(:, :, q), h(:, :, :ndof))
        end do
        ! Element-parallel assembly. Every element writes only its own slice;
        ! the node gather below runs in element order, so the result does not
        ! depend on the thread count. Each thread owns one quadrature workspace.
        allocate (element_status(nel), element_evaluations(nel), rhs_local(ndof, nel), &
                  source_local(nel))
        allocate (axis_local(ndof, merge(nel, 0, need_jacobian)))
        element_status = 0; element_evaluations = 0; rhs_local = 0; source_local = 0
        if (need_jacobian) axis_local = 0
        first_failure = huge(1)
        scratch_failed = .false.
        !$omp parallel default(shared)
        block
            type(load_scratch_t) :: scratch
            integer :: status_local, local_evaluations, ie, snapshot
            call init_load_scratch(scratch, nc, capacity, ndof, tol, jtol, status_local)
            if (status_local /= 0) then
                !$omp atomic write
                scratch_failed = .true.
            end if
            !$omp do schedule(dynamic, 4)
            do ie = 1, nel
                if (free_mode) then
                    if (.not. plasma_cell(ie)) cycle
                end if
                if (.not. scratch%ready) then
                    element_status(ie) = 1
                    cycle
                end if
                !$omp atomic read
                snapshot = first_failure
                if (ie > snapshot) cycle
                call gs_load_element(state, mu0, law, ie, free_mode, smooth_law, &
                                     need_jacobian, nc, tol, jtol, levels, tables, &
                                     forcing, scratch, local_evaluations, &
                                     status_local, cell_psi)
                element_evaluations(ie) = local_evaluations
                if (status_local /= 0) then
                    element_status(ie) = status_local
                    !$omp critical(gs_load_failure)
                    first_failure = min(first_failure, ie)
                    !$omp end critical(gs_load_failure)
                    cycle
                end if
                rhs_local(:, ie) = scratch%values(:ndof)
                source_local(ie) = sum(scratch%values(:ndof))
                if (need_jacobian) then
                    axis_local(:, ie) = scratch%values(nc - ndof + 1:nc)
                    jacobian%mass(:, :, ie) = &
                        reshape(scratch%values(ndof + 1:nc - ndof), [ndof, ndof])
                end if
            end do
            !$omp end do
        end block
        !$omp end parallel
        if (scratch_failed) then
            status = 1
            if (first_failure == huge(1)) return
        end if
        if (first_failure /= huge(1)) then
            status = element_status(first_failure)
            rhs = 0
            return
        end if
        evaluations = sum(element_evaluations)
        do e = 1, nel
            if (free_mode) then
                if (.not. plasma_cell(e)) cycle
            end if
            if (present(source_integral)) &
                source_integral = source_integral + source_local(e)
            do i = 1, ndof
                k = state%triangles(i, e)
                rhs(k) = rhs(k) + rhs_local(i, e)
                if (need_jacobian) jacobian%axis_load(k) = jacobian%axis_load(k) &
                    + axis_local(i, e)
            end do
            do j = 1, ndof
                k = state%triangles(j, e)
                if (free_mode) exit
                if (.not. state%boundary(k)) cycle
                do i = 1, ndof
                    rhs(state%triangles(i, e)) = rhs(state%triangles(i, e)) &
                                                 - state%stiffness(i, j, e)*state%psi(k)
                end do
            end do
        end do
        if (.not. free_mode) then
            where (state%boundary)
                rhs = state%psi
                jacobian%axis_gradient = 0
            end where
            if (need_jacobian) then
                where (state%boundary) jacobian%axis_load = 0
            end if
        end if
        if (need_jacobian .and. allocated(state%op%rowptr) .and. gs_use_csr()) then
            allocate (jacobian%mval(size(state%op%val)))
            call gs_csr_fill_masked(state%op, state%triangles, state%boundary, &
                                    jacobian%mass, jacobian%mval)
        end if
        status = 1
        if (all(ieee_is_finite(rhs))) status = 0
        if (present(quadrature_evaluations)) quadrature_evaluations = evaluations
    end subroutine

    subroutine init_load_scratch(scr, nc, capacity, ndof, tol, jtol, status)
        type(load_scratch_t), intent(out) :: scr
        integer, intent(in) :: nc, capacity, ndof
        real(dp), intent(in) :: tol, jtol
        integer, intent(out) :: status
        allocate (scr%events(capacity), scr%inner(capacity), &
                  scr%inner_active(capacity), &
                  scr%epsabs(nc), scr%rtol(nc), scr%values(nc), scr%errors(nc), &
                  scr%cell_values(nc), scr%high_values(nc))
        scr%rtol(:ndof) = tol; scr%rtol(ndof + 1:) = jtol
        call initialize_triangle_level_workspace(scr%workspace, nc, capacity, &
                   status, max_panels=max(256, capacity), max_evaluations=2000000, &
                                                 max_points=1000000)
        scr%ready = status == 0
    end subroutine

    subroutine gs_load_element(state, mu0, law, e, free_mode, smooth_law, &
                               need_jacobian, nc, tol, jtol, levels, tables, forcing, &
                               scr, evaluations, status, cell_psi)
        !! One element of gs_cubic_load. All mutable state is local or in the
        !! per-thread scratch `scr`, so concurrent calls on distinct elements are
        !! independent. The integrals are returned in scr%values.
        type(gs_state_t), intent(in) :: state
        real(dp), intent(in) :: mu0, tol, jtol, levels(:)
        type(gs_profile_law_t), intent(in) :: law
        integer, intent(in) :: e, nc
        logical, intent(in) :: free_mode, smooth_law, need_jacobian
        type(load_tables_t), intent(in) :: tables
        procedure(gs_forcing_t), optional :: forcing
        type(load_scratch_t), intent(inout) :: scr
        integer, intent(out) :: evaluations, status
        real(dp), intent(in), optional :: cell_psi(:, :)
        real(dp) :: coefficients(6), local_psi(10), controls(10), low, high, scale
        real(dp) :: edge_values(4), roots(3), shape(10), g(2, 10)
        real(dp) :: last_x, actual(2), width, samples(4), band_lower, band_upper
        real(dp) :: cell_points(2, 10)
        integer :: ndof, i, j, nr, ne, ni, nlevels, q
        integer :: ix, left, right, band, cell_left, cell_right, rule
        logical :: integrated, closed_cut
        character(160) :: message
        ndof = size(state%triangles, 1)
        evaluations = 0
        status = 0
        message = ''
        scr%values = 0
            cell_points(:, :ndof) = state%points(:, state%triangles(:, e))
            local_psi(:ndof) = state%psi(state%triangles(:, e))
            if (present(cell_psi)) local_psi(:ndof) = cell_psi(:ndof, e)
            if (maxval(abs(local_psi(:ndof) - local_psi(1))) <= 0) then
                if (gs_law_kink(law, local_psi(1))) then
                    status = 4
                    return
                end if
            end if
            call triangle_bernstein_controls(local_psi(:ndof), controls(:ndof))
            low = minval(controls(:ndof)); high = maxval(controls(:ndof))
            left = count(levels <= low) + 1; right = count(levels < high)
            nlevels = max(0, right - left + 1)
            closed_cut = .false.
            if (nlevels > 0 .and. smooth_law) closed_cut = possible_closed_cut()
            scr%values = 0
            scale = max(1.0_dp, mu0*maxval(abs(law%pp(1, :))) &
                        + maxval(abs(law%ff(1, :))))*2*state%area(e)
            scr%epsabs = scr%rtol*scale
            integrated = .false.
            if (nlevels > 0 .and. smooth_law .and. .not. closed_cut) then
                scr%values = 0
                do rule = 1, size(scr%rules)
                    call smooth_rule(rule, scr%high_values)
                    if (status /= 0) return
                    if (rule > 1) then
                        if (all(abs(scr%high_values - scr%values) <= scr%epsabs + &
                            scr%rtol*abs(scr%high_values))) then
                            integrated = .true.
                            exit
                        end if
                    end if
                    scr%values = scr%high_values
                end do
                scr%values = scr%high_values
            end if
            if (integrated) then
                continue
            else if (nlevels == 0) then
                do q = 1, size(tables%w)
                    call source_values(tables%shapes(:, q), tables%grads(:, :, q), &
                                       scr%errors, status)
                    if (status /= 0) return
                    scr%values = scr%values + tables%w(q)*scr%errors
                end do
                scr%high_values = 0
                do q = 1, size(tables%wh)
                    call source_values(tables%shapes_h(:, q), &
                                       tables%grads_h(:, :, q), scr%errors, status)
                    if (status /= 0) return
                    scr%high_values = scr%high_values + tables%wh(q)*scr%errors
                end do
                if (any(abs(scr%high_values - scr%values) > scr%epsabs + tol*0.1_dp &
                        *abs(scr%high_values))) then
                    coefficients = 0
                    call integrate_triangle_levels(coefficients, [real(dp) ::], &
                        source_products, scr%epsabs, tol*0.1_dp, scr%workspace, &
                        scr%values, scr%errors, status, message)
                    if (status /= 0) then
                        print *, 'smooth profile quadrature cell:', e, trim(message)
                        return
                    end if
                else
                    scr%values = scr%high_values
                end if
            else
                if (state%element_order == 3 .and. smooth_law .and. &
                    .not. closed_cut) then
                    band_lower = -huge(1.0_dp); band_upper = huge(1.0_dp)
                    ne = 2; scr%events(:2) = [0.0_dp, 1.0_dp]
                    coefficients = [0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
                    last_x = -1
                    call integrate_triangle_levels(coefficients, scr%events(:ne), &
                        cut_products, scr%epsabs, tol, scr%workspace, scr%values, &
                        scr%errors, status, message)
                else if (state%element_order <= 2) then
                    call triangle_level_coefficients(state%element_order, &
                                                 local_psi(:ndof), coefficients, status)
                    if (status /= 0) return
                    call integrate_triangle_levels(coefficients, levels(left:right), &
                        source_products, scr%epsabs, tol*0.1_dp, scr%workspace, &
                        scr%values, scr%errors, status, message)
                else
                    ! One flux band at a time: work grows with the number of
                    ! knots rather than repeating all bands at every event.
                    ! FortFEM owns the positive paired adaptive integration.
                    cell_left = left; cell_right = right; scr%cell_values = 0
                    do band = 0, nlevels
                        band_lower = -huge(1.0_dp); band_upper = huge(1.0_dp)
                        if (band > 0) band_lower = levels(cell_left + band - 1)
                        if (band < nlevels) band_upper = levels(cell_left + band)
                        left = max(cell_left, cell_left + band - 1)
                        right = min(cell_right, cell_left + band)
                        scr%events(:2) = [0.0_dp, 1.0_dp]; ne = 2
                        do ix = 1, 2
                            do i = 1, 4
                                actual = [real(i - 1, dp)/3, 0.0_dp]
                                if (ix == 2) actual(2) = 1 - actual(1)
                                edge_values(i) = flux_value(actual)
                            end do
                            do j = left, right
                             call cubic_segment_roots(edge_values, levels(j), roots, nr)
                                scr%events(ne + 1:ne + nr) = roots(:nr); ne = ne + nr
                            end do
                        end do
                        call tangency_events()
                        if (status /= 0) return
                        call sort_cuts(scr%events, ne)
                        call cut_rule(4, scr%values)
                        if (status /= 0) return
                        call cut_rule(8, scr%high_values)
                        if (status /= 0) return
                        scr%errors = abs(scr%high_values - scr%values)
                        scr%epsabs(:ndof) = max(tol*scale*0.01_dp/(nlevels + 1), &
                            tol*0.1_dp*maxval(abs(scr%high_values(:ndof))))
                        if (need_jacobian) then
                            scr%epsabs(ndof + 1:) = &
                                max(jtol*scale*0.01_dp/(nlevels + 1), &
                                    jtol*0.1_dp*maxval(abs(scr%high_values(ndof + 1:))))
                        end if
                        coefficients = [0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
                        last_x = -1
                        if (all(scr%errors <= scr%epsabs &
                                + tol*0.1_dp*abs(scr%high_values))) then
                            scr%values = scr%high_values
                        else
                            call integrate_triangle_levels(coefficients, &
                                  scr%events(:ne), cut_products, scr%epsabs, &
                                  tol*0.1_dp, scr%workspace, scr%values, &
                                  scr%errors, status, message)
                        end if
                        if (status /= 0) exit
                        scr%cell_values = scr%cell_values + scr%values
                    end do
                    scr%values = scr%cell_values
                end if
                if (status /= 0) then
                   print *, 'profile quadrature cell/status: ', e, status, trim(message)
                    print *, 'cuts, panels, evaluations, depth:', nlevels, &
                        scr%workspace%npanels, scr%workspace%nevaluations, &
                        scr%workspace%inner_depth
                    return
                end if
            end if
    contains
        logical function possible_closed_cut() result(closed)
            real(dp) :: edge_min, edge_max, t, value, stationary(2), v(4)
            integer :: edge, node, count
            edge_min = huge(1.0_dp); edge_max = -huge(1.0_dp)
            do edge = 1, 3
                do node = 1, 4
                    t = real(node - 1, dp)/3
                    v(node) = edge_flux(edge, t)
                end do
                edge_min = min(edge_min, v(1), v(4))
                edge_max = max(edge_max, v(1), v(4))
                call cubic_segment_stationary(v, stationary, count)
                do node = 1, count
                    value = edge_flux(edge, stationary(node))
                    edge_min = min(edge_min, value)
                    edge_max = max(edge_max, value)
                end do
            end do
            closed = any(levels(left:right) < edge_min) .or. &
                     any(levels(left:right) > edge_max)
        end function

        real(dp) function edge_flux(edge, t) result(value)
            integer, intent(in) :: edge
            real(dp), intent(in) :: t
            select case (edge)
            case (1)
                value = flux_value([t, 0.0_dp])
            case (2)
                value = flux_value([1 - t, t])
            case default
                value = flux_value([0.0_dp, 1 - t])
            end select
        end function

        subroutine smooth_rule(index, integral)
            integer, intent(in) :: index
            real(dp), intent(out) :: integral(:)
            integer, parameter :: orders(4) = [8, 12, 24, 48]
            real(dp), allocatable :: rx(:), ry(:)
            real(dp) :: piece(nc)
            integer :: node
            if (.not. allocated(scr%rules(index)%weights)) then
                call triangle_duffy_quadrature(2*orders(index) - 2, rx, ry, &
                    scr%rules(index)%weights, status)
                if (status /= 0) return
                allocate (scr%rules(index)%shapes(ndof, size(rx)), &
                    scr%rules(index)%gradients(2, ndof, size(rx)))
                do node = 1, size(rx)
                    call triangle_basis_jet([rx(node), ry(node)], &
                        scr%rules(index)%shapes(:, node), &
                        scr%rules(index)%gradients(:, :, node))
                end do
            end if
            integral = 0
            do node = 1, size(scr%rules(index)%weights)
                call source_values(scr%rules(index)%shapes(:, node), &
                    scr%rules(index)%gradients(:, :, node), piece, status)
                if (status /= 0) return
                integral = integral + scr%rules(index)%weights(node)*piece
            end do
        end subroutine

        subroutine cut_rule(order, integral)
            integer, intent(in) :: order
            real(dp), intent(out) :: integral(:)
            real(dp) :: nodes(order), weights(order), piece(nc), xi, eta, scale
            integer :: panel, a, b
            call gauss_legendre_ab(order, 0.0_dp, 1.0_dp, nodes, weights)
            integral = 0; last_x = -1
            do panel = 1, ne - 1
                scale = scr%events(panel + 1) - scr%events(panel)
                do a = 1, order
                    xi = scr%events(panel) + scale*nodes(a)
                    do b = 1, order
                        eta = (1 - xi)*nodes(b)
                        call cut_products([xi, eta], piece, status)
                        if (status /= 0) return
                        integral = integral + scale*(1 - xi) &
                                   *weights(a)*weights(b)*piece
                    end do
                end do
            end do
        end subroutine
        subroutine tangency_events()
            real(dp) :: slices(4, 7), v(4), d1, d2, d3, t, span, c(4)
            real(dp) :: vandermonde(7, 7), matrix(7, 7), polynomial(7, 1)
            real(dp) :: companion(6, 6), wr(6), wi(6), vl(1, 1), vr(1, 1), work(256)
            integer :: a, b, level, pivots(7), degree, info
            logical :: quadratic
            span = max(high - low, tiny(1.0_dp))
            do a = 1, 7
                t = real(a - 1, dp)/3 - 1
                vandermonde(a, 1) = 1
                do b = 2, 7
                    vandermonde(a, b) = vandermonde(a, b - 1)*t
                end do
                do b = 1, 4
                    v(b) = (flux_value([(t + 1)/2, real(b - 1, dp)/3]) &
                            - local_psi(1))/span
                end do
                d1 = v(2) - v(1); d2 = v(3) - 2*v(2) + v(1)
                d3 = v(4) - 3*v(3) + 3*v(2) - v(1)
                slices(:, a) = [v(1), 3*d1 - 1.5_dp*d2 + d3, &
                                4.5_dp*(d2 - d3), 4.5_dp*d3]
            end do
            quadratic = maxval(abs(slices(4, :))) < 128*epsilon(1.0_dp)
            do level = left, right
                do a = 1, 7
                    c = slices(:, a)
                    c(1) = c(1) + (local_psi(1) - levels(level))/span
                    if (quadratic) then
                        call gs_cut_quadratic_discriminant(c(1), c(2), c(3), &
                                                           polynomial(a, 1))
                    else
                        call gs_cut_discriminant(c(1), c(2), c(3), c(4), &
                                                 polynomial(a, 1))
                    end if
                end do
                matrix = vandermonde
                call dgesv(7, 1, matrix, 7, pivots, polynomial, 7, info)
                if (info /= 0) then
                    status = 3
                    return
                end if
                degree = 6
                do while (degree > 0)
                    if (abs(polynomial(degree + 1, 1)) > &
                        1.0e-11_dp*maxval(abs(polynomial))) exit
                    degree = degree - 1
                end do
                if (degree == 0) cycle
                companion = 0
                do a = 2, degree
                    companion(a, a - 1) = 1
                end do
                companion(:degree, degree) = -polynomial(:degree, 1) &
                                             /polynomial(degree + 1, 1)
                call dgeev('N', 'N', degree, companion, 6, wr, wi, vl, 1, &
                           vr, 1, work, size(work), info)
                if (info /= 0) then
                    status = 3
                    return
                end if
                do a = 1, degree
                    if (abs(wi(a)) > 1.0e-7_dp) cycle
                    if (wr(a) <= -1 .or. wr(a) >= 1) cycle
                    ne = ne + 1; scr%events(ne) = (wr(a) + 1)/2
                end do
            end do
        end subroutine
        real(dp) function flux_value(point) result(flux)
            real(dp), intent(in) :: point(2)
            call triangle_basis_jet(point, shape(:ndof), g(:, :ndof))
            flux = local_psi(1) + dot_product(shape(:ndof), &
                                              local_psi(:ndof) - local_psi(1))
        end function
        subroutine source_products(point, output, stat)
            real(dp), intent(in) :: point(2)
            real(dp), intent(out) :: output(:)
            integer, intent(out) :: stat
            real(dp) :: s(10), grad(2, 10)
            call triangle_basis_jet(point, s(:ndof), grad(:, :ndof))
            call source_values(s(:ndof), grad(:, :ndof), output, stat)
        end subroutine
        subroutine source_values(s, grad, output, stat)
            real(dp), intent(in) :: s(:), grad(:, :)
            real(dp), intent(out) :: output(:)
            integer, intent(out) :: stat
            real(dp) :: physical(2), det, psi, pp, ff, dpd, dfd, app, aff, f, df, af
            real(dp) :: map_jac(2, 2)
            integer :: a, b
            evaluations = evaluations + 1
            stat = 1
            physical = matmul(cell_points(:, :ndof), s)
            det = 2*state%area(e)
            if (state%curved_cell(e)) then
                map_jac = matmul(cell_points(:, :ndof), &
                                 transpose(grad))
                det = det2(map_jac)
            end if
            if (det <= 0 .or. physical(1) <= 0) return
            psi = local_psi(1) + dot_product(s(:ndof), &
                                             local_psi(:ndof) - local_psi(1))
            if (free_mode) then
                ! Outside the plasma boundary (normalised flux above one).
                if (psi*law%axis < 0) then
                    output = 0
                    stat = 0
                    return
                end if
            end if
            if (need_jacobian) then
                call gs_law_sources(law, psi, pp, ff, dpd, dfd, app, aff)
            else
                call gs_law_sources(law, psi, pp, ff)
            end if
            call gs_source(physical(1), mu0, pp, ff, f)
            if (present(forcing)) f = f + forcing(physical(1), physical(2))
            output(:ndof) = det*f*s(:ndof)
            if (.not. need_jacobian) then
                if (all(ieee_is_finite(output))) stat = 0
                return
            end if
            call gs_source_jvp(physical(1), mu0, pp, ff, 0.0_dp, 0.0_dp, dpd, dfd, df)
            call gs_source_jvp(physical(1), mu0, pp, ff, 0.0_dp, 0.0_dp, app, aff, af)
            do b = 1, ndof
                do a = 1, ndof
                    output(ndof + (b - 1)*ndof + a) = det*df*s(a)*s(b)
                end do
            end do
            output(ndof*(ndof + 1) + 1:) = det*af*s(:ndof)
            if (all(ieee_is_finite(output))) stat = 0
        end subroutine
        subroutine cut_products(point, output, stat)
            real(dp), intent(in) :: point(2)
            real(dp), intent(out) :: output(:)
            integer, intent(out) :: stat
            real(dp) :: piece(nc), eta, upper, midpoint_flux
            integer :: a, b, root_count
            if (abs(point(1) - last_x) > 0) then
                last_x = point(1)
                upper = 1 - last_x
                do a = 1, 4
                    samples(a) = flux_value([last_x, upper*real(a - 1, dp)/3])
                end do
                scr%inner(:2) = [0.0_dp, 1.0_dp]; ni = 2
                do b = left, right
                    call cubic_segment_roots(samples, levels(b), roots, root_count)
                    scr%inner(ni + 1:ni + root_count) = roots(:root_count)
                    ni = ni + root_count
                end do
                call sort_cuts(scr%inner, ni)
                do a = 1, ni - 1
                    midpoint_flux = flux_value([last_x, &
                                            upper*(scr%inner(a) + scr%inner(a + 1))/2])
                    scr%inner_active(a) = midpoint_flux >= band_lower .and. &
                                      midpoint_flux <= band_upper
                end do
            end if
            output = 0; stat = 0
            do a = 1, ni - 1
                if (.not. scr%inner_active(a)) cycle
                width = scr%inner(a + 1) - scr%inner(a)
                eta = (1 - point(1))*scr%inner(a) + width*point(2)
                call source_products([point(1), eta], piece, stat)
                if (stat /= 0) return
                output = output + width*piece
            end do
        end subroutine
    end subroutine


    subroutine apply_gs_profile_jacobian(state, jacobian, input, output, transpose)
        type(gs_state_t), intent(in) :: state
        type(gs_profile_jacobian_t), intent(in) :: jacobian
        real(dp), intent(in) :: input(:)
        real(dp), intent(out) :: output(:)
        logical, intent(in), optional :: transpose
        real(dp) :: local(size(state%triangles, 1)), product(size(local))
        integer :: e, i, j, node
        logical :: adjoint
        adjoint = .false.
        if (present(transpose)) adjoint = transpose
        if (allocated(jacobian%mval) .and. allocated(state%op%rowptr)) then
            call gs_csr_dirichlet_apply(state%op, state%boundary, input, output, &
                                        jacobian%mval)
            if (adjoint) then
                output = output - jacobian%axis_gradient*dot_product( &
                         jacobian%axis_load, input)
            else
                output = output - jacobian%axis_load*dot_product( &
                         jacobian%axis_gradient, input)
            end if
            return
        end if
        call apply_gs(state, input, output)
        do e = 1, size(state%triangles, 2)
            local = input(state%triangles(:, e))
            where (state%boundary(state%triangles(:, e))) local = 0
            ! The local source mass matrix is symmetric. Only the solved-axis
            ! normalization contributes a nonsymmetric, rank-one term.
            ! Explicit column-ordered product: avoids a libgfortran matmul call
            ! on a 10x10 or 20x20 block with runtime extents.
            product = 0
            do j = 1, size(local)
                do i = 1, size(local)
                    product(i) = product(i) + jacobian%mass(i, j, e)*local(j)
                end do
            end do
            do i = 1, size(local)
                node = state%triangles(i, e)
                if (state%boundary(node)) cycle
                output(node) = output(node) - product(i)
            end do
        end do
        if (adjoint) then
            output = output - jacobian%axis_gradient*dot_product( &
                     jacobian%axis_load, input)
        else
            output = output - jacobian%axis_load*dot_product( &
                     jacobian%axis_gradient, input)
        end if
    end subroutine

    subroutine gs_profile_response(state, jacobian, rhs, solution, status, transpose)
        type(gs_state_t), intent(in) :: state
        type(gs_profile_jacobian_t), intent(in) :: jacobian
        real(dp), intent(in) :: rhs(:)
        real(dp), intent(out) :: solution(:)
        integer, intent(out) :: status
        logical, intent(in), optional :: transpose
        real(dp) :: residual, direct_tolerance, denominator, coupling
        real(dp) :: preconditioned_rhs(size(rhs)), action_work(size(rhs))
        real(dp) :: diagonal(size(rhs)), axis_response(size(rhs))
        real(dp) :: rank_rhs(size(rhs)), rank_test(size(rhs)), clock
        integer :: iterations, precondition_status, e, j, node
        logical :: adjoint, use_inverse
        adjoint = .false.
        if (present(transpose)) adjoint = transpose
        precondition_status = 0
        solution = 0
        if (norm2(rhs) <= tiny(1.0_dp)) then
            status = 0
            return
        end if
        ! The axis-normalization term is rank one; the remaining block is
        ! symmetric. Try the existing CG solve on that block, rejecting
        ! negative curvature/breakdown and checking the original J residual.
        ! Positive diagonal alone is not a definiteness claim.
        rank_rhs = jacobian%axis_load; rank_test = jacobian%axis_gradient
        if (adjoint) then
            rank_rhs = jacobian%axis_gradient; rank_test = jacobian%axis_load
        end if
        ! With the cached Cholesky factor of A, A^-1 preconditions both Krylov
        ! solves: A^-1 (A - M) = I - A^-1 M is a compact perturbation of the
        ! identity, so the iteration counts do not grow with the mesh.
        use_inverse = gs_factor_ready(state)
        direct_tolerance = 1.0e-11_dp
        if (maxval(abs(rank_rhs))*maxval(abs(rank_test)) > tiny(1.0_dp)) &
            direct_tolerance = 1.0e-9_dp
        diagonal = state%diagonal
        do e = 1, size(state%triangles, 2)
            do j = 1, size(state%triangles, 1)
                node = state%triangles(j, e)
                diagonal(node) = diagonal(node) - jacobian%mass(j, j, e)
            end do
        end do
        where (state%boundary) diagonal = 1
        if (all(diagonal > 0)) then
            clock = gs_wall_s()
            call real_conjugate_gradient_operator(symmetric_action, rhs, solution, &
                direct_tolerance*min(1.0_dp, norm2(rhs)), &
                max(500, 10*ceiling(sqrt(real(size(rhs), dp)))), &
                status, iterations, residual, symmetric_precondition)
            call gs_timer_add(GS_T_JAC_CG, gs_wall_s() - clock, iterations)
            if (status == KRYLOV_OK) then
                call symmetric_action(solution, action_work)
                if (norm2(action_work - rhs) > 1.0e-9_dp*norm2(rhs)) status = 1
            end if
            if (status == KRYLOV_OK .and. &
                maxval(abs(rank_rhs))*maxval(abs(rank_test)) > tiny(1.0_dp)) then
                axis_response = 0
                clock = gs_wall_s()
                call real_conjugate_gradient_operator(symmetric_action, rank_rhs, &
                    axis_response, direct_tolerance*min(1.0_dp, norm2(rank_rhs)), &
                    max(500, 10*ceiling(sqrt(real(size(rhs), dp)))), &
                    status, iterations, residual, symmetric_precondition)
                call gs_timer_add(GS_T_JAC_CG, gs_wall_s() - clock, iterations)
                if (status == KRYLOV_OK) then
                    call symmetric_action(axis_response, action_work)
                    if (norm2(action_work - rank_rhs) > &
                        1.0e-9_dp*norm2(rank_rhs)) status = 1
                end if
                if (status == KRYLOV_OK) then
                    coupling = dot_product(rank_test, axis_response)
                    denominator = 1 - coupling
                    if (abs(denominator) <= sqrt(epsilon(1.0_dp))* &
                        max(1.0_dp, abs(coupling))) then
                        status = 1
                    else
                        solution = solution + axis_response* &
                            (dot_product(rank_test, solution)/denominator)
                    end if
                end if
            end if
            if (status == KRYLOV_OK) then
                call jacobian_action(solution, action_work)
                if (norm2(action_work - rhs) <= 1.0e-9_dp*norm2(rhs)) return
            end if
            solution = 0
        end if
        clock = gs_wall_s()
        call precondition(rhs, preconditioned_rhs)
        call real_gmres_operator(action, preconditioned_rhs, solution, &
                           1.0e-11_dp*min(1.0_dp, norm2(preconditioned_rhs)), 300, 30, &
                                 status, iterations, residual)
        call gs_timer_add(GS_T_JAC_GMRES, gs_wall_s() - clock, iterations)
        if (precondition_status /= 0) status = precondition_status
        ! GMRES certifies only the preconditioned system. Like the direct-CG
        ! branch above, confirm the physical Jacobian residual independently
        ! before accepting a finite stalled iterate; any other solver status
        ! remains a failure.
        if (status == KRYLOV_MAX_ITERATIONS .and. all(ieee_is_finite(solution))) then
            call jacobian_action(solution, action_work)
            if (norm2(action_work - rhs) <= 1.0e-9_dp*norm2(rhs)) status = KRYLOV_OK
        end if
    contains
        subroutine jacobian_action(input, output)
            real(dp), intent(in) :: input(:)
            real(dp), intent(out) :: output(:)
            call apply_gs_profile_jacobian(state, jacobian, input, output, adjoint)
        end subroutine
        subroutine symmetric_action(input, output)
            real(dp), intent(in) :: input(:)
            real(dp), intent(out) :: output(:)
            call jacobian_action(input, output)
            output = output + rank_rhs*dot_product(rank_test, input)
        end subroutine
        subroutine symmetric_precondition(input, output)
            ! A^-1 when the factor is cached, the Jacobi diagonal otherwise.
            real(dp), intent(in) :: input(:)
            real(dp), intent(out) :: output(:)
            logical :: ok
            if (use_inverse) then
                call gs_apply_inverse(state, input, output, ok)
                if (ok) return
            end if
            output = input/diagonal
        end subroutine
        subroutine action(input, output)
            real(dp), intent(in) :: input(:)
            real(dp), intent(out) :: output(:)
            call apply_gs_profile_jacobian(state, jacobian, input, action_work, adjoint)
            call precondition(action_work, output)
        end subroutine
        subroutine precondition(input, output)
            real(dp), intent(in) :: input(:)
            real(dp), intent(out) :: output(:)
            real(dp), allocatable :: candidate(:)
            integer, parameter :: precondition_nonfinite = 4
            integer :: stat
            logical :: ok
            if (use_inverse) then
                call gs_apply_inverse(state, input, output, ok)
                if (ok) return
            end if
            ! A preconditioner is not the solve: it may only fall back, never
            ! silently accept an unusable candidate. The diagonal-CG inner
            ! loop keeps the baseline 1e-11 accuracy so the inverse decks
            ! that reach it converge GMRES honestly; only the n>=320 TC24
            ! roundoff floor (~1.1e-10 relative at 95k DOF) that no
            ! iteration budget reaches down to 1e-11 retains its finite
            ! stalled iterate (stagnation_ok). The outer GMRES residual plus
            ! the Newton line search stay the authority for the step.
            call gs_linear_response(state, input, candidate, stat, &
                                    relative_tolerance=1.0e-11_dp, stagnation_ok=.true., &
                                    timer_kind=GS_T_PRECOND_CG)
            ! .and. does not short-circuit: allocation is checked on its own
            ! before any elementwise test may touch the candidate.
            if (.not. allocated(candidate)) then
                precondition_status = merge(stat, 1, stat /= 0)
                output = input
                return
            end if
            if (.not. all(ieee_is_finite(candidate))) then
                ! A nonfinite response is never usable, even when the inner
                ! solver reported success: only finite stagnation is tolerated.
                precondition_status = merge(stat, precondition_nonfinite, stat /= 0)
                output = input
                return
            end if
            if (stat /= KRYLOV_OK .and. stat /= KRYLOV_MAX_ITERATIONS) then
                ! Geometry, breakdown or error responses are never usable,
                ! finite or not: fall back to identity and flag the status.
                precondition_status = stat
                output = input
                return
            end if
            output = candidate
        end subroutine
    end subroutine

    subroutine solve_gs_cubic_profiles(state, mu0, law, status, relative_tolerance, &
                                       forcing, quadrature_tolerance)
        type(gs_state_t), intent(inout) :: state
        real(dp), intent(in) :: mu0
        type(gs_profile_law_t), intent(inout) :: law
        integer, intent(out) :: status
        real(dp), intent(in), optional :: relative_tolerance
        real(dp), intent(in), optional :: quadrature_tolerance
        procedure(gs_forcing_t), optional :: forcing
        type(gs_profile_jacobian_t) :: jacobian, trial_jacobian
        real(dp), allocatable :: rhs(:), action(:), step(:), previous(:), trial_rhs(:)
        real(dp) :: pp, ff, tol, damping, old_norm, new_norm, load_tolerance
        real(dp) :: work_tolerance, next_tolerance, previous_axis, clock, start_s
        character(160) :: line
        integer :: iteration, backtrack
        logical :: reuse_jacobian
        status = 1
        if (.not. state%ready .or. mu0 <= 0) return
        if (.not. ieee_is_finite(mu0)) return
        if (.not. law%valid()) return
        tol = 1.0e-9_dp
        if (present(relative_tolerance)) tol = relative_tolerance
        if (.not. ieee_is_finite(tol)) return
        if (tol <= 0 .or. tol >= 1) return
        if (maxval(abs(state%psi)) <= tiny(1.0_dp)) then
            call gs_law_sources(law, 0.0_dp, pp, ff)
            if (abs(pp) + abs(ff) <= tiny(1.0_dp)) &
                call gs_law_sources(law, 0.5_dp*law%axis, pp, ff)
            call solve_gs(state, mu0, pp, ff, status)
            if (status /= 0) return

        end if
        allocate (rhs(size(state%psi)), action(size(state%psi)), &
           step(size(state%psi)), previous(size(state%psi)), trial_rhs(size(state%psi)))
        state%iterations = 0
        ! Algebraic convergence and source integration are separate controls.
        load_tolerance = 1.0e-8_dp
        if (present(quadrature_tolerance)) load_tolerance = quadrature_tolerance
        if (.not. ieee_is_finite(load_tolerance)) return
        if (load_tolerance <= 0 .or. load_tolerance >= 1) return
        work_tolerance = max(1.0e-4_dp, load_tolerance)
        start_s = gs_wall_s()
        clock = start_s
        call gs_cubic_load(state, mu0, law, rhs, jacobian, status, work_tolerance, &
                           forcing, max(work_tolerance, 1.0e-3_dp))
        call gs_timer_add(GS_T_ASSEMBLY, gs_wall_s() - clock)
        if (status /= 0) return
        do iteration = 0, 25
            call apply_gs(state, state%psi, action)
            action = rhs - action
            old_norm = norm2(action)
            state%residual = old_norm/max(norm2(rhs), tiny(1.0_dp))
            next_tolerance = max(load_tolerance, &
                                 min(work_tolerance, 0.1_dp*state%residual**2))
            if (next_tolerance < 10*load_tolerance) next_tolerance = load_tolerance
            if (next_tolerance < 0.2_dp*work_tolerance .or. &
                (next_tolerance <= load_tolerance .and. work_tolerance > &
                 load_tolerance)) then
                work_tolerance = next_tolerance
                clock = gs_wall_s()
                call gs_cubic_load(state, mu0, law, rhs, jacobian, status, &
                                   work_tolerance, forcing, &
                                   max(work_tolerance, 1.0e-3_dp))
                call gs_timer_add(GS_T_ASSEMBLY, gs_wall_s() - clock)
                if (status /= 0) return
                call apply_gs(state, state%psi, action)
                action = rhs - action
                old_norm = norm2(action)
                state%residual = old_norm/max(norm2(rhs), tiny(1.0_dp))
            end if
            state%iterations = iteration
            if (gs_progress_enabled()) then
                write (line, '(a,i0,a,es10.3,a,es9.2,a,i0,a,i0,a,i0,a,f9.1)') &
                    'newton it=', iteration, ' residual=', state%residual, &
                    ' work_tol=', work_tolerance, ' jac_cg_its=', &
                    gs_timer_iterations(GS_T_JAC_CG), ' gmres_its=', &
                    gs_timer_iterations(GS_T_JAC_GMRES), ' precond_cg_its=', &
                    gs_timer_iterations(GS_T_PRECOND_CG), ' elapsed_s=', &
                    gs_wall_s() - start_s
                call gs_progress(line)
            end if
            if (state%residual <= tol .and. work_tolerance <= load_tolerance) then
                status = 0
                return
            end if
            if (iteration == 25) exit
            call gs_profile_response(state, jacobian, action, step, status)
            if (status /= 0) return
            previous = state%psi; previous_axis = law%axis; damping = 1
            reuse_jacobian = state%residual < 1.0e-3_dp
            do backtrack = 0, 14
                state%psi = previous + damping*step
                clock = gs_wall_s()
                call gs_cubic_load(state, mu0, law, trial_rhs, trial_jacobian, &
                                   status, work_tolerance, forcing, &
                                   max(work_tolerance, 1.0e-3_dp), &
                                   load_only=reuse_jacobian)
                call gs_timer_add(GS_T_ASSEMBLY, gs_wall_s() - clock)
                if (status == 0) then
                    call apply_gs(state, state%psi, action)
                    new_norm = norm2(trial_rhs - action)
                    if (new_norm <= (1 - 1.0e-4_dp*damping)*old_norm) exit
                end if
                damping = damping/2
            end do
            if (backtrack > 14) then
                state%psi = previous; law%axis = previous_axis
                status = 2
                return
            end if
            rhs = trial_rhs
            if (.not. reuse_jacobian) then
                call move_alloc(trial_jacobian%mass, jacobian%mass)
                call move_alloc(trial_jacobian%axis_load, jacobian%axis_load)
                call move_alloc(trial_jacobian%axis_gradient, jacobian%axis_gradient)
                call move_alloc(trial_jacobian%mval, jacobian%mval)
            end if
        end do
        status = 2
    end subroutine

    subroutine gs_normalize_current(state, mu0, law, target, current, status, tolerance)
        ! For normalized derivative shapes, psi -> a psi and (pprime,FFprime)
        ! -> a (pprime,FFprime) preserve s and the GS equation exactly. The
        ! primitives relative to the edge scale as a**2. F_edge stays fixed.
        ! This is an input current constraint, never a fit to another solution.
        type(gs_state_t), intent(inout) :: state
        real(dp), intent(in) :: mu0, target
        type(gs_profile_law_t), intent(inout) :: law
        real(dp), intent(out) :: current
        integer, intent(out) :: status
        real(dp), intent(in), optional :: tolerance
        type(gs_profile_jacobian_t) :: jac
        real(dp), allocatable :: rhs(:)
        real(dp) :: integral, scale
        current = 0; status = 1
        if (.not. law%normalized .or. .not. state%ready) return
        if (.not. ieee_is_finite(target) .or. abs(target) <= tiny(1.0_dp)) return
        if (any(abs(pack(state%psi, state%boundary)) > tiny(1.0_dp))) return
        allocate (rhs(size(state%psi)))
        call gs_cubic_load(state, mu0, law, rhs, jac, status, tolerance, &
                           source_integral=integral)
        if (status /= 0) return
        current = integral/mu0
        status = 1
        if (abs(current) <= tiny(1.0_dp)) return
        scale = target/current
        if (.not. ieee_is_finite(scale) .or. scale <= 0) return
        state%psi = scale*state%psi
        law%axis = scale*law%axis
        law%pp = scale*law%pp
        law%ff = scale*law%ff
        current = target
        status = 0
    end subroutine

    subroutine sort_cuts(cuts, n)
        real(dp), intent(inout) :: cuts(:)
        integer, intent(inout) :: n
        real(dp) :: value
        integer :: i, j, kept
        do i = 2, n
            value = cuts(i); j = i - 1
            do while (j >= 1)
                if (cuts(j) <= value) exit
                cuts(j + 1) = cuts(j); j = j - 1
            end do
            cuts(j + 1) = value
        end do
        kept = 0
        do i = 1, n
            if (kept > 0) then
                if (abs(cuts(i) - cuts(kept)) < 32*epsilon(1.0_dp)) cycle
            end if
            kept = kept + 1; cuts(kept) = cuts(i)
        end do
        n = kept
    end subroutine
end module
