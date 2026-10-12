module kin6d_gs_csr
    ! Row-parallel CSR kernels for the Grad--Shafranov operators, built on the
    ! fortfem element-to-CSR assembly. The pattern is the full symmetric one
    ! (diagonal stored). The stiffness is applied in exchange form,
    !   (A x)_i = sum_{j/=i} a_ij (x_j - x_i),
    ! so the zero row sum of every element block survives rounding; the stored
    ! diagonal then multiplies exactly zero. Dirichlet nodes follow
    ! apply_gs: their input is treated as zero inside interior rows and their
    ! output row is the identity.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortfem_element_to_csr, only: csr_matrix_t, csr_build_pattern, &
        csr_assemble_blocks, csr_find_entry
    implicit none
    private
    public :: csr_matrix_t, gs_csr_build, gs_csr_dirichlet_apply, gs_csr_fill_masked
    public :: gs_csr_spd_matrix, gs_csr_apply_full
    integer, parameter :: parallel_threshold = 2048
contains
    subroutine gs_csr_build(op, conn, blocks, n_dofs)
        ! Pattern and values of the element blocks blocks(a, b, e).
        type(csr_matrix_t), intent(out) :: op
        integer, intent(in) :: conn(:, :), n_dofs
        real(dp), intent(in) :: blocks(:, :, :)
        call csr_build_pattern(n_dofs, conn, .false., op)
        call csr_assemble_blocks(op, conn, blocks)
    end subroutine

    subroutine gs_csr_dirichlet_apply(op, boundary, x, y, mval)
        ! y = A x (exchange form, Dirichlet rows identity) minus M x when the
        ! mass values mval (same pattern, boundary rows and columns zero) are given.
        type(csr_matrix_t), intent(in) :: op
        logical, intent(in) :: boundary(:)
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: y(:)
        real(dp), intent(in), optional :: mval(:)
        real(dp), allocatable :: xm(:)
        real(dp) :: acc, mass
        integer :: i, k
        allocate (xm(op%n))
        !$omp parallel do default(none) shared(op, boundary, x, xm) private(i) &
        !$omp schedule(static) if(op%n >= parallel_threshold)
        do i = 1, op%n
            if (boundary(i)) then
                xm(i) = 0
            else
                xm(i) = x(i)
            end if
        end do
        !$omp end parallel do
        if (present(mval)) then
            !$omp parallel do default(none) shared(op, boundary, x, xm, y, mval) &
            !$omp private(i, k, acc, mass) schedule(static) &
            !$omp if(op%n >= parallel_threshold)
            do i = 1, op%n
                if (boundary(i)) then
                    y(i) = x(i)
                    cycle
                end if
                acc = 0; mass = 0
                do k = op%rowptr(i), op%rowptr(i + 1) - 1
                    acc = acc + op%val(k)*(xm(op%colidx(k)) - xm(i))
                    mass = mass + mval(k)*xm(op%colidx(k))
                end do
                y(i) = acc - mass
            end do
            !$omp end parallel do
        else
            !$omp parallel do default(none) shared(op, boundary, x, xm, y) &
            !$omp private(i, k, acc) schedule(static) &
            !$omp if(op%n >= parallel_threshold)
            do i = 1, op%n
                if (boundary(i)) then
                    y(i) = x(i)
                    cycle
                end if
                acc = 0
                do k = op%rowptr(i), op%rowptr(i + 1) - 1
                    acc = acc + op%val(k)*(xm(op%colidx(k)) - xm(i))
                end do
                y(i) = acc
            end do
            !$omp end parallel do
        end if
    end subroutine

    subroutine gs_csr_apply_full(op, x, y)
        ! y = M x for a full-form matrix.
        type(csr_matrix_t), intent(in) :: op
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: y(:)
        real(dp) :: acc
        integer :: i, k
        !$omp parallel do default(none) shared(op, x, y) private(i, k, acc) &
        !$omp schedule(static) if(op%n >= parallel_threshold)
        do i = 1, op%n
            acc = 0
            do k = op%rowptr(i), op%rowptr(i + 1) - 1
                acc = acc + op%val(k)*x(op%colidx(k))
            end do
            y(i) = acc
        end do
        !$omp end parallel do
    end subroutine

    subroutine gs_csr_fill_masked(op, conn, boundary, blocks, mval)
        ! Values of sum_e blocks(:, :, e) on the pattern of op, restricted to
        ! interior rows and columns. Each row is summed by one thread in
        ! ascending element order.
        type(csr_matrix_t), intent(in) :: op
        integer, intent(in) :: conn(:, :)
        logical, intent(in) :: boundary(:)
        real(dp), intent(in) :: blocks(:, :, :)
        real(dp), intent(inout) :: mval(:)
        integer :: i, k, e, a, b, j, p
        !$omp parallel do default(none) shared(op, conn, boundary, blocks, mval) &
        !$omp private(i, k, e, a, b, j, p) schedule(static) &
        !$omp if(op%n >= parallel_threshold)
        do i = 1, op%n
            mval(op%rowptr(i):op%rowptr(i + 1) - 1) = 0
            if (boundary(i)) cycle
            do k = op%elem_ptr(i), op%elem_ptr(i + 1) - 1
                e = op%elem_idx(k)
                do a = 1, size(conn, 1)
                    if (conn(a, e) /= i) cycle
                    do b = 1, size(conn, 1)
                        j = conn(b, e)
                        if (boundary(j)) cycle
                        p = csr_find_entry(op, i, j)
                        mval(p) = mval(p) + blocks(a, b, e)
                    end do
                end do
            end do
        end do
        !$omp end parallel do
    end subroutine

    subroutine gs_csr_spd_matrix(op, boundary, rowptr, colidx, val)
        ! The matrix that apply realises: interior rows keep interior columns
        ! only, with the diagonal -sum_{j/=i} a_ij over all neighbours (the
        ! exchange-form diagonal); Dirichlet rows are the identity.
        type(csr_matrix_t), intent(in) :: op
        logical, intent(in) :: boundary(:)
        integer, allocatable, intent(out) :: rowptr(:), colidx(:)
        real(dp), allocatable, intent(out) :: val(:)
        real(dp) :: diagonal
        integer(8) :: total
        integer :: i, k, q, count_row
        allocate (rowptr(op%n + 1))
        rowptr(1) = 1
        !$omp parallel do default(none) shared(op, boundary, rowptr) &
        !$omp private(i, k, count_row) schedule(static) &
        !$omp if(op%n >= parallel_threshold)
        do i = 1, op%n
            if (boundary(i)) then
                count_row = 1
            else
                count_row = 0
                do k = op%rowptr(i), op%rowptr(i + 1) - 1
                    if (.not. boundary(op%colidx(k)) .or. op%colidx(k) == i) &
                        count_row = count_row + 1
                end do
            end if
            rowptr(i + 1) = count_row
        end do
        !$omp end parallel do
        total = 1
        do i = 1, op%n
            total = total + rowptr(i + 1)
            if (total > int(huge(1), 8)) error stop 'gs_csr_spd_matrix: nnz overflow'
            rowptr(i + 1) = int(total)
        end do
        allocate (colidx(rowptr(op%n + 1) - 1), val(rowptr(op%n + 1) - 1))
        !$omp parallel do default(none) shared(op, boundary, rowptr, colidx, val) &
        !$omp private(i, k, q, diagonal) schedule(static) &
        !$omp if(op%n >= parallel_threshold)
        do i = 1, op%n
            q = rowptr(i)
            if (boundary(i)) then
                colidx(q) = i; val(q) = 1
                cycle
            end if
            diagonal = 0
            do k = op%rowptr(i), op%rowptr(i + 1) - 1
                if (op%colidx(k) /= i) diagonal = diagonal - op%val(k)
            end do
            do k = op%rowptr(i), op%rowptr(i + 1) - 1
                if (op%colidx(k) == i) then
                    colidx(q) = i; val(q) = diagonal
                    q = q + 1
                else if (.not. boundary(op%colidx(k))) then
                    colidx(q) = op%colidx(k); val(q) = op%val(k)
                    q = q + 1
                end if
            end do
        end do
        !$omp end parallel do
    end subroutine
end module
