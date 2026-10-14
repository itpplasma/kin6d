module kin6d_gs_ampere_loop
    !! Total signed toroidal current on the represented outer FEM boundary.
    !! Positive reference-edge orientation is CCW in (R,Z); I=-int B.dl/mu0.
    !! Integrate each smooth curved edge separately, using its owning cell's
    !! one-sided physical gradient. No field interpolation or inward offset.
    !! Fixed-geometry nodal weights provide the linear JVP and transpose VJP.
    !! Successive quadrature differences are indicators, not rigorous bounds.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use kin6d_grad_shafranov, only: gs_state_t, gs_cell_data
    use kin6d_triangle_basis, only: triangle_basis_jet
    use kin6d_grad_shafranov_leaves, only: gs_current_sheet
    use kin6d_gs_profile_law, only: gs_profile_law_t
    use kin6d_gs_nonlinear, only: gs_cubic_load, gs_profile_jacobian_t, gs_forcing_t
    use kin6d_real_identity, only: same_real_value
    use fortfem_triangle_duffy_quadrature, only: triangle_duffy_quadrature
    use fortnum_quadrature, only: gauss_legendre_ab
    implicit none
    private
    public :: gs_boundary_current, gs_weak_boundary_current
contains
    subroutine gs_weak_boundary_current(state, mu0, law, order, current, status, &
                                        nodal_weights, source_current, lift, &
                                        tolerance, forcing)
        !! Green's identity: I=(load(v)-a(psi,v))/mu0, trace(v)=1.
        !! The default FE lift is one on boundary DOFs and zero elsewhere.
        !! Interior residuals control dependence on the choice of lift. For a
        !! Galerkin solution this equals the volume source, unlike the raw trace.
        !! Fixed-geometry/law coefficient weights give JVP and dual VJP; for a
        !! normalized law they include the equation-level stationary-axis term.
        !! Smooth-law/unique-axis premises are inherited from gs_cubic_load.
        !! Geometry derivatives and rigorous observable bounds are unavailable.
        type(gs_state_t), intent(inout) :: state
        real(dp), intent(in) :: mu0
        type(gs_profile_law_t), intent(in) :: law
        integer, intent(in) :: order
        real(dp), intent(out) :: current
        integer, intent(out) :: status
        real(dp), intent(out), optional :: nodal_weights(:), source_current
        real(dp), intent(in), optional :: lift(:), tolerance
        procedure(gs_forcing_t), optional :: forcing
        type(gs_profile_law_t) :: executed
        type(gs_profile_jacobian_t) :: jacobian
        real(dp), allocatable :: rhs(:), test(:), x(:), y(:), w(:)
        real(dp) :: point(2), det, shape(10), g(2, 10), h(2, 2, 10)
        real(dp) :: grad(2), test_grad(2), local(10), density, integral
        integer :: e, q, j, node, ndof, st

        status = 1; current = 0
        if (present(nodal_weights)) nodal_weights = 0
        if (present(source_current)) source_current = 0
        if (state%element_order < 2 .or. state%element_order > 3) return
        if (order < 1) return
        if (mu0 <= 0) return
        if (.not. ieee_is_finite(mu0)) return
        if (.not. any(state%boundary)) return
        if (present(nodal_weights)) then
            if (size(nodal_weights) /= size(state%psi)) return
        end if
        allocate (test(size(state%psi)), rhs(size(state%psi)))
        test = merge(1.0_dp, 0.0_dp, state%boundary)
        if (present(lift)) then
            if (size(lift) /= size(test)) return
            if (.not. all(ieee_is_finite(lift))) return
            if (any(abs(pack(lift, state%boundary) - 1.0_dp) > 0)) return
            test = lift
        end if
        executed = law
        call gs_cubic_load(state, mu0, executed, rhs, jacobian, st, &
                           tolerance=tolerance, forcing=forcing, &
                           load_only=.not. present(nodal_weights), &
                           source_integral=integral, uneliminated=.true.)
        if (st /= 0) then
            status = st
            return
        end if
        current = dot_product(test, rhs)/mu0
        if (present(source_current)) source_current = integral/mu0
        ndof = size(state%triangles, 1)
        if (present(nodal_weights)) then
            nodal_weights = jacobian%axis_gradient &
                            *dot_product(test, jacobian%axis_load)/mu0
            do e = 1, size(state%triangles, 2)
                local(:ndof) = test(state%triangles(:, e))
                do j = 1, ndof
                    node = state%triangles(j, e)
                    density = dot_product(local(:ndof), jacobian%mass(:, j, e))/mu0
                    nodal_weights(node) = nodal_weights(node) + density
                end do
            end do
        end if
        call triangle_duffy_quadrature(order, x, y, w, st)
        if (st /= 0) return
        do e = 1, size(state%triangles, 2)
            local(:ndof) = test(state%triangles(:, e))
            if (all(same_real_value(local(:ndof), local(1)))) cycle
            do q = 1, size(w)
                call gs_cell_data(state, e, [x(q), y(q)], point, det, &
                                  shape(:ndof), g(:, :ndof), h(:, :, :ndof), st)
                if (st /= 0 .or. point(1) <= 0) return
                grad = matmul(g(:, :ndof), state%psi(state%triangles(:, e)) &
                              - state%psi(state%triangles(1, e)))
                test_grad = matmul(g(:, :ndof), local(:ndof) - local(1))
                call gs_current_sheet(point(1), mu0, grad(1), grad(2), &
                                      test_grad(1), test_grad(2), density)
                current = current + w(q)*det*density
                if (.not. present(nodal_weights)) cycle
                do j = 1, ndof
                    node = state%triangles(j, e)
                    call gs_current_sheet(point(1), mu0, g(1, j), g(2, j), &
                                          test_grad(1), test_grad(2), density)
                    nodal_weights(node) = nodal_weights(node) + w(q)*det*density
                end do
            end do
        end do
        if (.not. ieee_is_finite(current)) return
        status = 0
    end subroutine

    subroutine gs_boundary_current(state, mu0, order, current, status, nodal_weights)
        type(gs_state_t), intent(in) :: state
        real(dp), intent(in) :: mu0
        integer, intent(in) :: order
        real(dp), intent(out) :: current
        integer, intent(out) :: status
        real(dp), intent(out), optional :: nodal_weights(:)
        real(dp), parameter :: vertex(2, 3) = reshape( &
                               [0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2, 3])
        real(dp), allocatable :: x(:), w(:)
        real(dp) :: point(2), det, shape(10), g(2, 10), h(2, 2, 10)
        real(dp) :: ref_g(2, 10), xi(2), direction(2), tangent(2), grad(2)
        real(dp) :: density, contribution, basis_current
        integer :: e, k, next, first, last, q, j, ndof, st, edges

        status = 1; current = 0
        if (present(nodal_weights)) nodal_weights = 0
        if (state%element_order < 2 .or. state%element_order > 3) return
        if (order < 1) return
        if (mu0 <= 0) return
        if (.not. ieee_is_finite(mu0)) return
        if (present(nodal_weights)) then
            if (size(nodal_weights) /= size(state%psi)) return
        end if
        ndof = size(state%triangles, 1)
        allocate (x(order), w(order))
        call gauss_legendre_ab(order, 0.0_dp, 1.0_dp, x, w)
        edges = 0
        do e = 1, size(state%triangles, 2)
            do k = 1, 3
                next = mod(k, 3) + 1
                if (.not. state%boundary(state%triangles(k, e))) cycle
                if (.not. state%boundary(state%triangles(next, e))) cycle
                first = 4 + (k - 1)*(state%element_order - 1)
                last = 3 + k*(state%element_order - 1)
                if (.not. all(state%boundary(state%triangles(first:last, e)))) cycle
                edges = edges + 1
                direction = vertex(:, next) - vertex(:, k)
                do q = 1, order
                    xi = vertex(:, k) + x(q)*direction
                    call gs_cell_data(state, e, xi, point, det, shape(:ndof), &
                                      g(:, :ndof), h(:, :, :ndof), st)
                    if (st /= 0 .or. point(1) <= 0) return
                    call triangle_basis_jet(xi, shape(:ndof), ref_g(:, :ndof))
                    if (state%curved_cell(e)) then
                        tangent = matmul(state%points(:, state%triangles(:, e)), &
                                         matmul(transpose(ref_g(:, :ndof)), direction))
                    else
                        tangent = state%points(:, state%triangles(next, e)) &
                                  - state%points(:, state%triangles(k, e))
                    end if
                    grad = matmul(g(:, :ndof), state%psi(state%triangles(:, e)) &
                                  - state%psi(state%triangles(1, e)))
                    call gs_current_sheet(point(1), mu0, grad(1), grad(2), &
                                          tangent(2), -tangent(1), density)
                    contribution = w(q)*density
                    if (.not. ieee_is_finite(contribution)) return
                    current = current + contribution
                    if (.not. present(nodal_weights)) cycle
                    do j = 1, ndof
                        call gs_current_sheet(point(1), mu0, g(1, j), g(2, j), &
                                              tangent(2), -tangent(1), basis_current)
                        nodal_weights(state%triangles(j, e)) = &
                            nodal_weights(state%triangles(j, e)) + w(q)*basis_current
                    end do
                end do
            end do
        end do
        if (edges == 0) return
        status = 0
    end subroutine
end module
