module kin6d_gs_majorant
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use kin6d_grad_shafranov, only: gs_state_t, gs_cell_data, apply_gs, &
                                  gs_linear_response
    use kin6d_grad_shafranov_leaves, only: gs_source, gs_source_jvp
    use kin6d_triangle_basis, only: triangle_bernstein_controls
    use kin6d_gs_profile_law, only: gs_profile_law_t, gs_law_sources
    use kin6d_gs_nonlinear, only: gs_forcing_t
    use fortfem_triangle_duffy_quadrature, only: triangle_duffy_quadrature
    use fortnum_krylov, only: real_conjugate_gradient_operator, KRYLOV_OK, &
        KRYLOV_MAX_ITERATIONS
    use kin6d_gs_csr, only: csr_matrix_t, gs_csr_build, gs_csr_apply_full
    use kin6d_gs_options, only: gs_use_direct
    use kin6d_gs_direct, only: gs_direct_t, gs_direct_available, gs_direct_factor, &
        gs_direct_solve, gs_direct_free
    implicit none
    private
    public :: gs_functional_majorant
contains
    subroutine gs_functional_majorant(state, mu0, pprime, ffprime, majorant, &
                                     flux_defect, balance_defect, status, &
                                     recovery_converged, recovery_iterations, law, &
                                     stability_lower_bound, forcing, &
                                     estimate_stability, stability_estimate)
        ! Repin functional majorant for -div(grad(psi)/R)=source, homogeneous
        ! Dirichlet trace, on the represented domain. Recover y in continuous
        ! vector Pk, k=2 or 3 (a subspace of H(div)) by minimizing D^2+C_F^2*B^2.
        ! Returns D+C_F*B, with ordinary quadrature; no rigorous enclosure.
        type(gs_state_t), intent(in) :: state
        real(dp), intent(in) :: mu0, pprime, ffprime
        real(dp), intent(out) :: majorant, flux_defect, balance_defect
        integer, intent(out) :: status
        logical, intent(out), optional :: recovery_converged
        integer, intent(out), optional :: recovery_iterations
        type(gs_profile_law_t), intent(in), optional :: law
        real(dp), intent(in), optional :: stability_lower_bound
        procedure(gs_forcing_t), optional :: forcing
        logical, intent(in), optional :: estimate_stability
        real(dp), intent(out), optional :: stability_estimate
        real(dp), allocatable :: x(:), y(:), w(:), matrix(:, :, :)
        real(dp), allocatable :: rhs(:), recovered(:), diagonal(:)
        integer, allocatable :: indices(:, :)
        real(dp), allocatable :: shape(:), gradients(:, :), hessians(:, :, :)
        real(dp), allocatable :: local(:, :), d(:), bezier(:, :)
        real(dp), allocatable :: source_mass(:, :, :)
        real(dp) :: point(2), det
        real(dp) :: grad(2), source, weight, cf2, lo(2), hi(2)
        real(dp) :: flux(2), divergence, residual, pp, ff, psi
        real(dp) :: dpp, dff, derivative, alpha
        logical :: practical, solved
        integer :: e, q, j, k, a, b, row, col, n, iterations, ndof
        majorant = huge(1.0_dp); flux_defect = 0; balance_defect = 0
        status = 1
        practical = .false.
        if (present(estimate_stability)) practical = estimate_stability
        if (.not. present(law)) practical = .false.
        if (present(law)) then
            if (law%normalized) practical = .false.
        end if
        if (present(stability_lower_bound)) practical = .false.
        if (present(stability_estimate)) stability_estimate = 0
        if (present(recovery_converged)) recovery_converged = .false.
        if (present(recovery_iterations)) recovery_iterations = 0
        if (.not. state%ready) return
        if (state%element_order /= 2 .and. state%element_order /= 3) return
        ndof = size(state%triangles, 1)
        if (practical) then
            allocate (source_mass(ndof, ndof, size(state%triangles, 2)))
            source_mass = 0
        end if
        allocate (shape(ndof), gradients(2, ndof), hessians(2, 2, ndof), &
            local(2, ndof), d(ndof), bezier(2, ndof))
        if (any(abs(pack(state%psi, state%boundary)) > tiny(1.0_dp))) return
        if (mu0 <= 0) return
        if (.not. all(ieee_is_finite([mu0, pprime, ffprime]))) return
        lo = huge(1.0_dp); hi = -huge(1.0_dp)
        do e = 1, size(state%triangles, 2)
            bezier = state%points(:, state%triangles(:, e))
            do j = 1, 2
                call triangle_bernstein_controls( &
                    state%points(j, state%triangles(:, e)), bezier(j, :))
            end do
            lo = min(lo, minval(bezier, dim=2))
            hi = max(hi, maxval(bezier, dim=2))
        end do
        if (lo(1) <= 0 .or. any(hi <= lo)) return
        ! Zero extension to the enclosing rectangle gives this Poincare bound.
        cf2 = hi(1)/(acos(-1.0_dp)**2*sum(1/(hi - lo)**2))
        n = 2*size(state%psi)
        allocate (rhs(n), recovered(n), diagonal(n), &
            matrix(2*ndof, 2*ndof, size(state%triangles, 2)), &
            indices(2*ndof, size(state%triangles, 2)))
        rhs = 0; diagonal = 0; matrix = 0; recovered = 0
        call triangle_duffy_quadrature(8, x, y, w, status)
        if (status /= 0) return
        do e = 1, size(state%triangles, 2)
            do j = 1, ndof
                indices(2*j - 1:2*j, e) = &
                    [2*state%triangles(j, e) - 1, 2*state%triangles(j, e)]
            end do
            d = state%psi(state%triangles(:, e)) - state%psi(state%triangles(1, e))
            do q = 1, size(w)
                call gs_cell_data(state, e, [x(q), y(q)], point, det, shape, &
                    gradients, hessians, status)
                if (status /= 0) return
                weight = w(q)*det
                grad = matmul(gradients, d)
                pp = pprime; ff = ffprime
                if (present(law)) then
                    psi = state%psi(state%triangles(1, e)) + dot_product(shape, d)
                    if (practical) then
                        call gs_law_sources(law, psi, pp, ff, dpp, dff)
                        call gs_source_jvp(point(1), mu0, pp, ff, &
                            0.0_dp, 0.0_dp, dpp, dff, derivative)
                        do k = 1, ndof
                            do j = 1, ndof
                                source_mass(j, k, e) = source_mass(j, k, e) &
                                    + weight*max(0.0_dp, derivative)*shape(j)*shape(k)
                            end do
                        end do
                    else
                        call gs_law_sources(law, psi, pp, ff)
                    end if
                end if
                call gs_source(point(1), mu0, pp, ff, source)
                if (present(forcing)) source = source + forcing(point(1), point(2))
                do j = 1, ndof
                    do a = 1, 2
                        row = 2*(j - 1) + a
                        rhs(indices(row, e)) = rhs(indices(row, e)) + weight &
                            *(shape(j)*grad(a) - cf2*source*gradients(a, j))
                        do k = 1, ndof
                            do b = 1, 2
                                col = 2*(k - 1) + b
                                matrix(row, col, e) = matrix(row, col, e) &
                                    + weight*cf2*gradients(a, j)*gradients(b, k)
                                if (a == b) matrix(row, col, e) = matrix(row, col, e) &
                                    + weight*point(1)*shape(j)*shape(k)
                            end do
                        end do
                    end do
                end do
            end do
            do j = 1, 2*ndof
                diagonal(indices(j, e)) = diagonal(indices(j, e)) + matrix(j, j, e)
            end do
        end do
        if (present(recovery_converged)) recovery_converged = .true.
        if (norm2(rhs) > tiny(1.0_dp)) then
            solved = .false.
            if (gs_use_direct() .and. gs_direct_available()) &
                call direct_recovery(solved)
            if (solved) then
                status = KRYLOV_OK
            else
                recovered = 0
                call real_conjugate_gradient_operator(action, rhs, recovered, &
                    1.0e-11_dp*min(1.0_dp, norm2(rhs)), 6000, status, iterations, &
                    residual, precondition)
            end if
            if (present(recovery_converged)) recovery_converged = status == KRYLOV_OK
            if (present(recovery_iterations)) recovery_iterations = iterations
            ! The functional estimate is valid for every continuous recovered
            ! vector field, not only the minimizer. An iteration cap can make
            ! it less sharp, but must not discard an otherwise valid estimate.
            ! Algebraic breakdown and invalid arguments remain failures.
            if (status /= KRYLOV_OK .and. status /= KRYLOV_MAX_ITERATIONS) return
            if (.not. all(ieee_is_finite(recovered))) then
                status = 1
                return
            end if
        end if
        do e = 1, size(state%triangles, 2)
            local = reshape(recovered(indices(:, e)), [2, ndof])
            d = state%psi(state%triangles(:, e)) - state%psi(state%triangles(1, e))
            do q = 1, size(w)
                call gs_cell_data(state, e, [x(q), y(q)], point, det, shape, &
                    gradients, hessians, status)
                if (status /= 0) return
                weight = w(q)*det
                grad = matmul(gradients, d)
                flux = matmul(local, shape)
                divergence = sum(local*gradients)
                pp = pprime; ff = ffprime
                if (present(law)) then
                    psi = state%psi(state%triangles(1, e)) + dot_product(shape, d)
                    call gs_law_sources(law, psi, pp, ff)
                end if
                call gs_source(point(1), mu0, pp, ff, source)
                if (present(forcing)) source = source + forcing(point(1), point(2))
                flux_defect = flux_defect + weight*point(1) &
                    *sum((flux - grad/point(1))**2)
                balance_defect = balance_defect + weight*(source + divergence)**2
            end do
        end do
        flux_defect = sqrt(flux_defect)
        balance_defect = sqrt(cf2*balance_defect)
        majorant = flux_defect + balance_defect
        if (.not. ieee_is_finite(majorant)) then
            status = 1
            return
        end if
        if (present(law)) then
            ! Without branch stability this is the frozen-source defect, not
            ! a nonlinear error bound. A caller-supplied alpha is conditional
            ! on the documented strong-monotonicity / normalized-axis assumptions.
            status = 4
            if (practical) then
                call estimate_coercivity(alpha)
                if (present(stability_estimate)) stability_estimate = alpha
                if (alpha > 0) then
                    majorant = majorant/alpha
                    status = 5
                end if
            end if
            if (present(stability_lower_bound)) then
                if (.not. ieee_is_finite(stability_lower_bound)) return
                if (stability_lower_bound <= 0) return
                majorant = majorant/stability_lower_bound
                status = 0
            end if
        end if
    contains
        subroutine direct_recovery(converged)
            ! The recovery system is symmetric positive definite: one Cholesky
            ! factor, one back-substitution and iterative refinement. When the
            ! factor or the refinement fails, the caller runs the CG path.
            logical, intent(out) :: converged
            type(csr_matrix_t) :: system
            type(gs_direct_t) :: factor
            real(dp), allocatable :: residual_vector(:), correction(:)
            real(dp) :: rnorm, previous, target
            logical :: ok
            integer :: sweep
            converged = .false.
            rnorm = huge(1.0_dp)
            call gs_csr_build(system, indices, matrix, n)
            call gs_direct_factor(factor, n, system%rowptr, system%colidx, &
                                  system%val, ok)
            if (.not. ok) then
                call gs_direct_free(factor)
                return
            end if
            allocate (residual_vector(n), correction(n))
            call gs_direct_solve(factor, rhs, recovered, ok)
            iterations = 1
            previous = huge(1.0_dp)
            target = 1.0e-11_dp*min(1.0_dp, norm2(rhs))
            do sweep = 0, 4
                if (.not. ok) exit
                call gs_csr_apply_full(system, recovered, residual_vector)
                residual_vector = rhs - residual_vector
                rnorm = norm2(residual_vector)
                if (rnorm <= 1.0e-2_dp*target .or. sweep == 4) exit
                if (rnorm > 0.5_dp*previous) exit
                previous = rnorm
                call gs_direct_solve(factor, residual_vector, correction, ok)
                recovered = recovered + correction
                iterations = iterations + 1
            end do
            converged = ok .and. rnorm <= target
            call gs_direct_free(factor)
        end subroutine

        subroutine estimate_coercivity(alpha)
            ! Numerical positive-part comparison: J >= A - M_+.
            ! Power iteration is self-adjoint in the A pairing. Its result is
            ! a practical discrete stability estimate, not a continuous bound.
            real(dp), intent(out) :: alpha
            real(dp), allocatable :: v(:), av(:), mv(:), z(:), defect(:), ad(:)
            real(dp) :: theta, radius, energy
            real(dp) :: element(ndof)
            integer :: iteration, cell, node, i, j, stat, count
            count = size(state%psi)
            allocate (v(count), av(count), mv(count), defect(count), ad(count))
            v = state%psi
            if (norm2(v) <= tiny(1.0_dp)) v = 1
            where (state%boundary) v = 0
            alpha = 0
            do iteration = 1, 60
                call apply_gs(state, v, av)
                energy = dot_product(v, av)
                if (energy <= 0) return
                v = v/sqrt(energy); av = av/sqrt(energy)
                mv = 0
                do cell = 1, size(state%triangles, 2)
                    element = 0
                    do j = 1, ndof
                        do i = 1, ndof
                            element(i) = element(i) + source_mass(i, j, cell)* &
                                v(state%triangles(j, cell))
                        end do
                    end do
                    do j = 1, ndof
                        node = state%triangles(j, cell)
                        if (state%boundary(node)) cycle
                        mv(node) = mv(node) + element(j)
                    end do
                end do
                theta = dot_product(v, mv)
                if (norm2(mv) <= tiny(1.0_dp)) then
                    alpha = 1
                    return
                end if
                call gs_linear_response(state, mv, z, stat, &
                    relative_tolerance=1.0e-9_dp)
                if (stat /= 0) return
                defect = z - theta*v
                call apply_gs(state, defect, ad)
                radius = sqrt(max(0.0_dp, dot_product(defect, ad)))
                if (radius <= 1.0e-5_dp*max(1.0_dp, abs(theta))) then
                    alpha = 1 - theta - radius
                    return
                end if
                v = z
            end do
        end subroutine

        subroutine action(input, output)
            real(dp), intent(in) :: input(:)
            real(dp), intent(out) :: output(:)
            real(dp) :: value(2*ndof)
            integer :: cell, i, j
            output = 0
            do cell = 1, size(indices, 2)
                value = 0
                do j = 1, 2*ndof
                    do i = 1, 2*ndof
                        value(i) = value(i) + matrix(i, j, cell)*input(indices(j, cell))
                    end do
                end do
                do i = 1, 2*ndof
                    output(indices(i, cell)) = output(indices(i, cell)) + value(i)
                end do
            end do
        end subroutine
        subroutine precondition(input, output)
            real(dp), intent(in) :: input(:)
            real(dp), intent(out) :: output(:)
            output = input/diagonal
        end subroutine
    end subroutine
end module
