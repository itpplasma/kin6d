module kin6d_gs_nonlinear
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use kin6d_grad_shafranov, only: gs_state_t, solve_gs, apply_gs, &
                                    gs_cell_data, gs_axis_candidate, gs_linear_response
    use kin6d_gs_profile_law, only: gs_profile_law_t, gs_law_sources, gs_law_kink
    use kin6d_grad_shafranov_leaves, only: gs_source, gs_source_jvp, &
                                      gs_cut_discriminant, gs_cut_quadratic_discriminant
    use kin6d_triangle_basis, only: triangle_basis_jet, cubic_segment_roots, &
                                    triangle_bernstein_controls
    use kin6d_curved_triangle, only: triangle_locate
    use fortfem_triangle_level_quadrature, only: triangle_level_workspace_t, &
                       initialize_triangle_level_workspace, integrate_triangle_levels, &
                                                 triangle_level_coefficients
    use fortfem_triangle_duffy_quadrature, only: triangle_duffy_quadrature
    use fortnum_krylov, only: real_gmres_operator
    use fortnum_linalg, only: det2
    use fortnum_quadrature, only: gauss_legendre_ab
    implicit none
    private
    public :: gs_forcing_t, gs_normalize_current
    public :: solve_gs_cubic_profiles, gs_cubic_load, gs_profile_axis
    public :: gs_profile_jacobian_t, apply_gs_profile_jacobian, gs_profile_response
    type :: gs_profile_jacobian_t
        real(dp), allocatable :: mass(:, :, :), axis_load(:), axis_gradient(:)
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
                             forcing, source_integral)
        type(gs_state_t), intent(inout) :: state
        real(dp), intent(in) :: mu0
        type(gs_profile_law_t), intent(inout) :: law
        real(dp), intent(out) :: rhs(:)
        type(gs_profile_jacobian_t), intent(out) :: jacobian
        integer, intent(out) :: status
        real(dp), intent(in), optional :: tolerance
        procedure(gs_forcing_t), optional :: forcing
        real(dp), intent(out), optional :: source_integral
        type(triangle_level_workspace_t) :: workspace
        real(dp), allocatable :: levels(:), events(:), inner(:), epsabs(:)
        real(dp), allocatable :: values(:), errors(:), x(:), y(:), w(:)
        real(dp), allocatable :: xh(:), yh(:), wh(:), shapes(:, :), grads(:, :, :)
        real(dp), allocatable :: shapes_h(:, :), grads_h(:, :, :)
        real(dp), allocatable :: high_values(:), cell_values(:)
        logical, allocatable :: inner_active(:)
        real(dp) :: coefficients(6), local_psi(10), controls(10), low, high, scale, tol
        real(dp) :: edge_values(4), roots(3), shape(10), g(2, 10), h(2, 2, 10)
        real(dp) :: last_x, actual(2), width, samples(4), band_lower, band_upper
        integer :: e, ndof, nc, i, j, k, nr, ne, ni, nlevels, q
        integer :: ix, left, right, band, cell_left, cell_right, capacity
        character(160) :: message
        status = 1; rhs = 0
        if (present(source_integral)) source_integral = 0
        if (.not. state%ready .or. mu0 <= 0) return
        if (.not. ieee_is_finite(mu0)) return
        if (.not. law%valid()) return
        if (.not. all(ieee_is_finite(state%psi))) return
        tol = 1.0e-9_dp
        if (present(tolerance)) tol = tolerance
        if (.not. ieee_is_finite(tol)) return
        if (tol <= 0 .or. tol >= 1) return
        ndof = size(state%triangles, 1); nc = ndof*(ndof + 2)
        allocate (jacobian%mass(ndof, ndof, size(state%triangles, 2)), &
                  jacobian%axis_load(size(rhs)), jacobian%axis_gradient(size(rhs)))
        jacobian%mass = 0; jacobian%axis_load = 0; jacobian%axis_gradient = 0
        if (law%normalized) then
            call gs_profile_axis(state, law%axis, jacobian%axis_gradient, status)
            if (status /= 0) return
        end if
        levels = law%knots
        if (law%normalized) levels = law%axis*(1 - levels)
        if (levels(1) > levels(size(levels))) levels = levels(size(levels):1:-1)
        capacity = 12*size(levels) + 3
        if (state%element_order == 3) capacity = 27
        allocate (events(capacity), inner(capacity), inner_active(capacity), &
                  epsabs(nc), values(nc), errors(nc), cell_values(nc))
        call initialize_triangle_level_workspace(workspace, nc, size(events), &
                   status, max_panels=max(256, size(events)), max_evaluations=2000000, &
                                                 max_points=1000000)
        if (status /= 0) return
        call triangle_duffy_quadrature(2*state%element_order + 2, x, y, w, status)
        if (status /= 0) return
        call triangle_duffy_quadrature(2*state%element_order + 4, xh, yh, wh, status)
        if (status /= 0) return
        allocate (shapes(ndof, size(w)), grads(2, ndof, size(w)), &
                  shapes_h(ndof, size(wh)), grads_h(2, ndof, size(wh)), high_values(nc))
        do q = 1, size(w)
            call triangle_basis_jet([x(q), y(q)], shapes(:, q), grads(:, :, q), &
                                    h(:, :, :ndof))
        end do
        do q = 1, size(wh)
            call triangle_basis_jet([xh(q), yh(q)], shapes_h(:, q), grads_h(:, :, q), &
                                    h(:, :, :ndof))
        end do
        do e = 1, size(state%triangles, 2)
            local_psi(:ndof) = state%psi(state%triangles(:, e))
            if (maxval(abs(local_psi(:ndof) - local_psi(1))) <= 0) then
                if (gs_law_kink(law, local_psi(1))) then
                    status = 4
                    rhs = 0
                    return
                end if
            end if
            call triangle_bernstein_controls(local_psi(:ndof), controls(:ndof))
            low = minval(controls(:ndof)); high = maxval(controls(:ndof))
            left = count(levels <= low) + 1; right = count(levels < high)
            nlevels = max(0, right - left + 1)
            values = 0
            scale = max(1.0_dp, mu0*maxval(abs(law%pp(1, :))) &
                        + maxval(abs(law%ff(1, :))))*2*state%area(e)
            epsabs = tol*scale*0.01_dp
            if (nlevels == 0) then
                do q = 1, size(w)
                    call source_values(shapes(:, q), grads(:, :, q), errors, status)
                    if (status /= 0) return
                    values = values + w(q)*errors
                end do
                high_values = 0
                do q = 1, size(wh)
                    call source_values(shapes_h(:, q), grads_h(:, :, q), errors, status)
                    if (status /= 0) return
                    high_values = high_values + wh(q)*errors
                end do
                if (any(abs(high_values - values) > epsabs + tol*0.1_dp &
                        *abs(high_values))) then
                    coefficients = 0
                    call integrate_triangle_levels(coefficients, [real(dp) ::], &
                               source_products, epsabs, tol*0.1_dp, workspace, values, &
                                                   errors, status, message)
                    if (status /= 0) then
                        print *, 'smooth profile quadrature cell:', e, trim(message)
                        return
                    end if
                else
                    values = high_values
                end if
            else
                if (state%element_order <= 2) then
                    call triangle_level_coefficients(state%element_order, &
                                                 local_psi(:ndof), coefficients, status)
                    if (status /= 0) return
                    call integrate_triangle_levels(coefficients, levels(left:right), &
                               source_products, epsabs, tol*0.1_dp, workspace, values, &
                                                   errors, status, message)
                else
                    ! One flux band at a time: work grows with the number of
                    ! knots rather than repeating all bands at every event.
                    ! FortFEM owns the positive paired adaptive integration.
                    cell_left = left; cell_right = right; cell_values = 0
                    do band = 0, nlevels
                        band_lower = -huge(1.0_dp); band_upper = huge(1.0_dp)
                        if (band > 0) band_lower = levels(cell_left + band - 1)
                        if (band < nlevels) band_upper = levels(cell_left + band)
                        left = max(cell_left, cell_left + band - 1)
                        right = min(cell_right, cell_left + band)
                        events(:2) = [0.0_dp, 1.0_dp]; ne = 2
                        do ix = 1, 2
                            do i = 1, 4
                                actual = [real(i - 1, dp)/3, 0.0_dp]
                                if (ix == 2) actual(2) = 1 - actual(1)
                                edge_values(i) = flux_value(actual)
                            end do
                            do j = left, right
                             call cubic_segment_roots(edge_values, levels(j), roots, nr)
                                events(ne + 1:ne + nr) = roots(:nr); ne = ne + nr
                            end do
                        end do
                        call tangency_events()
                        if (status /= 0) return
                        call sort_cuts(events, ne)
                        call cut_rule(4, values)
                        if (status /= 0) return
                        call cut_rule(8, high_values)
                        if (status /= 0) return
                        errors = abs(high_values - values)
                        epsabs(:ndof) = max(tol*scale*0.01_dp/(nlevels + 1), &
                                            tol*0.1_dp*maxval(abs(high_values(:ndof))))
                        epsabs(ndof + 1:) = &
                            max(tol*scale*0.01_dp/(nlevels + 1), &
                                tol*0.1_dp*maxval(abs(high_values(ndof + 1:))))
                        coefficients = [0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
                        last_x = -1
                        if (all(errors <= epsabs + tol*0.1_dp*abs(high_values))) then
                            values = high_values
                        else
                            call integrate_triangle_levels(coefficients, events(:ne), &
                                  cut_products, epsabs, tol*0.1_dp, workspace, values, &
                                                           errors, status, message)
                        end if
                        if (status /= 0) exit
                        cell_values = cell_values + values
                    end do
                    values = cell_values
                end if
                if (status /= 0) then
                   print *, 'profile quadrature cell/status: ', e, status, trim(message)
                    print *, 'cuts, panels, evaluations, depth:', nlevels, &
                        workspace%npanels, workspace%nevaluations, workspace%inner_depth
                    return
                end if
            end if
            if (present(source_integral)) &
                source_integral = source_integral + sum(values(:ndof))
            do i = 1, ndof
                k = state%triangles(i, e)
                rhs(k) = rhs(k) + values(i)
                jacobian%axis_load(k) = jacobian%axis_load(k) &
                                        + values(nc - ndof + i)
            end do
            jacobian%mass(:, :, e) = reshape(values(ndof + 1:nc - ndof), &
                                             [ndof, ndof])
            do j = 1, ndof
                k = state%triangles(j, e)
                if (.not. state%boundary(k)) cycle
                do i = 1, ndof
                    rhs(state%triangles(i, e)) = rhs(state%triangles(i, e)) &
                                                 - state%stiffness(i, j, e)*state%psi(k)
                end do
            end do
        end do
        where (state%boundary)
            rhs = state%psi
            jacobian%axis_load = 0
            jacobian%axis_gradient = 0
        end where
        status = 1
        if (all(ieee_is_finite(rhs))) status = 0
    contains
        subroutine cut_rule(order, integral)
            integer, intent(in) :: order
            real(dp), intent(out) :: integral(:)
            real(dp) :: nodes(order), weights(order), piece(nc), xi, eta, scale
            integer :: panel, a, b
            call gauss_legendre_ab(order, 0.0_dp, 1.0_dp, nodes, weights)
            integral = 0; last_x = -1
            do panel = 1, ne - 1
                scale = events(panel + 1) - events(panel)
                do a = 1, order
                    xi = events(panel) + scale*nodes(a)
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
                    ne = ne + 1; events(ne) = (wr(a) + 1)/2
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
            stat = 1
            physical = matmul(state%points(:, state%triangles(:, e)), s)
            det = 2*state%area(e)
            if (state%curved_cell(e)) then
                map_jac = matmul(state%points(:, state%triangles(:, e)), &
                                 transpose(grad))
                det = det2(map_jac)
            end if
            if (det <= 0 .or. physical(1) <= 0) return
            psi = local_psi(1) + dot_product(s(:ndof), &
                                             local_psi(:ndof) - local_psi(1))
            call gs_law_sources(law, psi, pp, ff, dpd, dfd, app, aff)
            call gs_source(physical(1), mu0, pp, ff, f)
            if (present(forcing)) f = f + forcing(physical(1), physical(2))
            output(:ndof) = det*f*s(:ndof)
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
                inner(:2) = [0.0_dp, 1.0_dp]; ni = 2
                do b = left, right
                    call cubic_segment_roots(samples, levels(b), roots, root_count)
                    inner(ni + 1:ni + root_count) = roots(:root_count)
                    ni = ni + root_count
                end do
                call sort_cuts(inner, ni)
                do a = 1, ni - 1
                    midpoint_flux = flux_value([last_x, &
                                                upper*(inner(a) + inner(a + 1))/2])
                    inner_active(a) = midpoint_flux >= band_lower .and. &
                                      midpoint_flux <= band_upper
                end do
            end if
            output = 0; stat = 0
            do a = 1, ni - 1
                if (.not. inner_active(a)) cycle
                width = inner(a + 1) - inner(a)
                eta = (1 - point(1))*inner(a) + width*point(2)
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
        integer :: e, i, node
        logical :: adjoint
        adjoint = .false.
        if (present(transpose)) adjoint = transpose
        call apply_gs(state, input, output)
        do e = 1, size(state%triangles, 2)
            local = input(state%triangles(:, e))
            where (state%boundary(state%triangles(:, e))) local = 0
            ! The local source mass matrix is symmetric. Only the solved-axis
            ! normalization contributes a nonsymmetric, rank-one term.
            product = matmul(jacobian%mass(:, :, e), local)
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
        real(dp) :: residual
        real(dp) :: preconditioned_rhs(size(rhs)), action_work(size(rhs))
        integer :: iterations, precondition_status
        logical :: adjoint
        adjoint = .false.
        if (present(transpose)) adjoint = transpose
        precondition_status = 0
        solution = 0
        if (norm2(rhs) <= tiny(1.0_dp)) then
            status = 0
            return
        end if
        call precondition(rhs, preconditioned_rhs)
        call real_gmres_operator(action, preconditioned_rhs, solution, &
                           1.0e-11_dp*min(1.0_dp, norm2(preconditioned_rhs)), 300, 30, &
                                 status, iterations, residual)
        if (precondition_status /= 0) status = precondition_status
    contains
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
            integer :: stat
            call gs_linear_response(state, input, candidate, stat, &
                                    relative_tolerance=1.0e-11_dp)
            if (stat /= 0) precondition_status = stat
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
        real(dp) :: work_tolerance, next_tolerance, previous_axis
        integer :: iteration, backtrack
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
        call gs_cubic_load(state, mu0, law, rhs, jacobian, status, work_tolerance, &
                           forcing)
        if (status /= 0) return
        do iteration = 0, 25
            call apply_gs(state, state%psi, action)
            action = rhs - action
            old_norm = norm2(action)
            state%residual = old_norm/max(norm2(rhs), tiny(1.0_dp))
            next_tolerance = max(load_tolerance, &
                                 min(work_tolerance, 0.1_dp*state%residual**2))
            if (next_tolerance < 0.2_dp*work_tolerance .or. &
                (next_tolerance <= load_tolerance .and. work_tolerance > &
                 load_tolerance)) then
                work_tolerance = next_tolerance
                call gs_cubic_load(state, mu0, law, rhs, jacobian, status, &
                                   work_tolerance, forcing)
                if (status /= 0) return
                call apply_gs(state, state%psi, action)
                action = rhs - action
                old_norm = norm2(action)
                state%residual = old_norm/max(norm2(rhs), tiny(1.0_dp))
            end if
            state%iterations = iteration
            if (state%residual <= tol .and. work_tolerance <= load_tolerance) then
                status = 0
                return
            end if
            if (iteration == 25) exit
            call gs_profile_response(state, jacobian, action, step, status)
            if (status /= 0) return
            previous = state%psi; previous_axis = law%axis; damping = 1
            do backtrack = 0, 14
                state%psi = previous + damping*step
                call gs_cubic_load(state, mu0, law, trial_rhs, trial_jacobian, &
                                   status, work_tolerance, forcing)
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
            call move_alloc(trial_jacobian%mass, jacobian%mass)
            call move_alloc(trial_jacobian%axis_load, jacobian%axis_load)
            call move_alloc(trial_jacobian%axis_gradient, jacobian%axis_gradient)
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
