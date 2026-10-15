module kin6d_gs_axis_recovery
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_grad_shafranov, only: gs_state_t, gs_axis_candidate
    implicit none
    private
    public :: gs_magnetic_axis
    interface
        subroutine dgels(trans, m, n, nrhs, a, lda, b, ldb, work, lwork, info)
            import :: dp
            character, intent(in) :: trans
            integer, intent(in) :: m, n, nrhs, lda, ldb, lwork
            integer, intent(out) :: info
            real(dp), intent(inout) :: a(lda, *), b(ldb, *)
            real(dp), intent(out) :: work(*)
        end subroutine
    end interface
contains
    subroutine gs_magnetic_axis(state, axis, status)
        type(gs_state_t), intent(in) :: state
        real(dp), intent(out) :: axis(2)
        integer, intent(out) :: status
        integer, parameter :: terms = 15
        real(dp), allocatable :: matrix(:, :), rhs(:, :), work(:)
        real(dp) :: center(2), scale, local(2), hessian(2, 2), gradient(2)
        real(dp) :: flux, coefficients(terms), step(2), determinant, weight
        integer :: seed, element, m, row, i, j, k, degree, px(terms), py(terms)
        integer :: iteration, info

        ! The represented minimum remains the contour/axis-q reference. Recover
        ! a smooth local flux for the physical axis-location diagnostic: its
        ! stationary point is less sensitive to which P2 cell contains the axis.
        call gs_axis_candidate(state, axis, flux, hessian, element, status)
        if (status /= 0) return
        ! P3 already resolves an interior stationary point with positive Hessian.
        ! A patch fit can move it in response to flux outside its containing cell.
        if (state%element_order == 3) return
        status = 1
        seed = minloc(state%psi, dim=1)
        center = state%points(:, seed)
        scale = 6*sqrt(2*state%area(element))
        m = 0
        do i = 1, size(state%psi)
            if (norm2(state%points(:, i) - center) < scale) m = m + 1
        end do
        if (m < 2*terms) return
        allocate (matrix(m, terms), rhs(m, 1), work(64*m))
        k = 0
        do degree = 0, 4
            do j = 0, degree
                k = k + 1
                px(k) = j; py(k) = degree - j
            end do
        end do
        row = 0
        do i = 1, size(state%psi)
            local = (state%points(:, i) - center)/scale
            if (norm2(local) >= 1) cycle
            row = row + 1
            weight = (1 - sum(local**2))**2
            do k = 1, terms
                matrix(row, k) = weight*local(1)**px(k)*local(2)**py(k)
            end do
            rhs(row, 1) = weight*(state%psi(i) - state%psi(seed))
        end do
        ! QR least squares, in coordinates scaled to the shrinking local patch.
        call dgels('N', m, terms, 1, matrix, m, rhs, m, work, size(work), info)
        if (info /= 0) return
        coefficients = rhs(:terms, 1)
        local = (axis - center)/scale
        do iteration = 1, 12
            gradient = 0; hessian = 0
            do k = 1, terms
                i = px(k); j = py(k)
                if (i > 0) gradient(1) = gradient(1) &
                    + coefficients(k)*i*local(1)**(i - 1)*local(2)**j
                if (j > 0) gradient(2) = gradient(2) &
                    + coefficients(k)*j*local(1)**i*local(2)**(j - 1)
                if (i > 1) hessian(1, 1) = hessian(1, 1) &
                    + coefficients(k)*i*(i - 1)*local(1)**(i - 2)*local(2)**j
                if (j > 1) hessian(2, 2) = hessian(2, 2) &
                    + coefficients(k)*j*(j - 1)*local(1)**i*local(2)**(j - 2)
                if (i > 0 .and. j > 0) hessian(1, 2) = hessian(1, 2) &
                    + coefficients(k)*i*j*local(1)**(i - 1)*local(2)**(j - 1)
            end do
            hessian(2, 1) = hessian(1, 2)
            determinant = hessian(1, 1)*hessian(2, 2) - hessian(1, 2)**2
            if (hessian(1, 1) <= 0 .or. determinant <= 0) return
            step = [hessian(2, 2)*gradient(1) - hessian(1, 2)*gradient(2), &
                hessian(1, 1)*gradient(2) - hessian(1, 2)*gradient(1)]/determinant
            local = local - step
            if (norm2(local) >= 1) return
            if (norm2(step) < 32*epsilon(1.0_dp)) exit
        end do
        if (iteration > 12) return
        axis = center + scale*local
        status = 0
    end subroutine
end module
