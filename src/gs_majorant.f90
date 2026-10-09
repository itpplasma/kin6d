module kin6d_gs_majorant
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use kin6d_grad_shafranov, only: gs_state_t, gs_cell_data
    use kin6d_grad_shafranov_leaves, only: gs_source
    use fortfem_triangle_duffy_quadrature, only: triangle_duffy_quadrature
    use fortnum_krylov, only: real_conjugate_gradient_operator, KRYLOV_OK
    implicit none
    private
    public :: gs_functional_majorant
contains
    subroutine gs_functional_majorant(state, mu0, pprime, ffprime, majorant, &
                                     flux_defect, balance_defect, status)
        ! Repin functional majorant for -div(grad(psi)/R)=source, homogeneous
        ! Dirichlet trace, on the represented domain. Recover y in continuous
        ! vector P2 (a subspace of H(div)) by minimizing D^2+C_F^2*B^2.
        ! Returns D+C_F*B, with ordinary quadrature; no rigorous enclosure.
        type(gs_state_t), intent(in) :: state
        real(dp), intent(in) :: mu0, pprime, ffprime
        real(dp), intent(out) :: majorant, flux_defect, balance_defect
        integer, intent(out) :: status
        real(dp), allocatable :: x(:), y(:), w(:), matrix(:, :, :)
        real(dp), allocatable :: rhs(:), recovered(:), diagonal(:)
        integer, allocatable :: indices(:, :)
        real(dp) :: shape(6), gradients(2, 6), hessians(2, 2, 6), point(2), det
        real(dp) :: grad(2), source, weight, cf2, lo(2), hi(2), bezier(2, 6)
        real(dp) :: flux(2), divergence, residual, local(2, 6), d(6)
        integer :: e, q, j, k, a, b, row, col, n, iterations, next
        majorant = huge(1.0_dp); flux_defect = 0; balance_defect = 0
        status = 1
        if (.not. state%ready .or. state%element_order /= 2) return
        if (any(abs(pack(state%psi, state%boundary)) > tiny(1.0_dp))) return
        if (mu0 <= 0) return
        if (.not. all(ieee_is_finite([mu0, pprime, ffprime]))) return
        lo = huge(1.0_dp); hi = -huge(1.0_dp)
        do e = 1, size(state%triangles, 2)
            bezier = state%points(:, state%triangles(:, e))
            do j = 1, 3
                next = mod(j, 3) + 1
                bezier(:, j + 3) = 2*bezier(:, j + 3) &
                    - (bezier(:, j) + bezier(:, next))/2
            end do
            lo = min(lo, minval(bezier, dim=2))
            hi = max(hi, maxval(bezier, dim=2))
        end do
        if (lo(1) <= 0 .or. any(hi <= lo)) return
        ! Zero extension to the enclosing rectangle gives this Poincare bound.
        cf2 = hi(1)/(acos(-1.0_dp)**2*sum(1/(hi - lo)**2))
        n = 2*size(state%psi)
        allocate (rhs(n), recovered(n), diagonal(n), &
            matrix(12, 12, size(state%triangles, 2)), &
            indices(12, size(state%triangles, 2)))
        rhs = 0; diagonal = 0; matrix = 0; recovered = 0
        call triangle_duffy_quadrature(8, x, y, w, status)
        if (status /= 0) return
        do e = 1, size(state%triangles, 2)
            do j = 1, 6
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
                call gs_source(point(1), mu0, pprime, ffprime, source)
                do j = 1, 6
                    do a = 1, 2
                        row = 2*(j - 1) + a
                        rhs(indices(row, e)) = rhs(indices(row, e)) + weight &
                            *(shape(j)*grad(a) - cf2*source*gradients(a, j))
                        do k = 1, 6
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
            do j = 1, 12
                diagonal(indices(j, e)) = diagonal(indices(j, e)) + matrix(j, j, e)
            end do
        end do
        if (norm2(rhs) > tiny(1.0_dp)) then
            call real_conjugate_gradient_operator(action, rhs, recovered, &
                1.0e-11_dp*min(1.0_dp, norm2(rhs)), 6000, status, iterations, &
                residual, precondition)
            if (status /= KRYLOV_OK) return
        end if
        do e = 1, size(state%triangles, 2)
            local = reshape(recovered(indices(:, e)), [2, 6])
            d = state%psi(state%triangles(:, e)) - state%psi(state%triangles(1, e))
            do q = 1, size(w)
                call gs_cell_data(state, e, [x(q), y(q)], point, det, shape, &
                    gradients, hessians, status)
                if (status /= 0) return
                weight = w(q)*det
                grad = matmul(gradients, d)
                flux = matmul(local, shape)
                divergence = sum(local*gradients)
                call gs_source(point(1), mu0, pprime, ffprime, source)
                flux_defect = flux_defect + weight*point(1) &
                    *sum((flux - grad/point(1))**2)
                balance_defect = balance_defect + weight*(source + divergence)**2
            end do
        end do
        flux_defect = sqrt(flux_defect)
        balance_defect = sqrt(cf2*balance_defect)
        majorant = flux_defect + balance_defect
        if (.not. ieee_is_finite(majorant)) status = 1
    contains
        subroutine action(input, output)
            real(dp), intent(in) :: input(:)
            real(dp), intent(out) :: output(:)
            real(dp) :: value(12)
            integer :: cell, i
            output = 0
            do cell = 1, size(indices, 2)
                value = matmul(matrix(:, :, cell), input(indices(:, cell)))
                do i = 1, 12
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
