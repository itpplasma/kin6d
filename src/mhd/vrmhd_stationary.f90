module kin6d_vrmhd_stationary
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_vrmhd_parameters, only: vrmhd_parameters_t
    use kin6d_vrmhd_section, only: vrmhd_section_t, vrmhd_section_basis, &
        vrmhd_section_residual, vrmhd_section_jacobian
    use kin6d_vrmhd_section_leaves, only: vrmhd_pairing
    use fortsparse, only: csc_t, csc_from_triplet, csc_matvec, csc_transpose, &
        fortsparse_status_t, status_ok
    implicit none
    private
    type, public :: vrmhd_stationary_t
        integer :: ndof = 0, ngauge = 0
        type(csc_t) :: gradient, transpose
    end type
    public :: vrmhd_stationary_build, vrmhd_stationary_residual
    public :: vrmhd_stationary_jacobian
contains
    subroutine vrmhd_stationary_build(op, steady, status)
        type(vrmhd_section_t), intent(in) :: op
        type(vrmhd_stationary_t), intent(out) :: steady
        integer, intent(out) :: status
        type(fortsparse_status_t) :: st
        integer, allocatable :: index(:), rows(:), cols(:)
        real(dp), allocatable :: values(:)
        real(dp) :: basis(36, op%nlocal), block(op%nlocal, op%scalar_fs%nl0)
        real(dp) :: x(2), weight, value
        integer :: i, j, t, q, pos, count, row, col, nt
        nt = op%fs%mm%topology%n_triangles
        allocate (index(op%scalar_fs%n0))
        index = 0
        do i = 1, op%scalar_fs%n0
            if (op%scalar_fs%sp(0)%bnd(i)) cycle
            steady%ngauge = steady%ngauge + 1
            index(i) = steady%ngauge
        end do
        steady%ndof = op%ndof + steady%ngauge
        count = nt*op%nlocal*op%scalar_fs%nl0
        allocate (rows(count), cols(count), values(count))
        pos = 0
        do t = 1, nt
            block = 0
                do q = 1, op%fs%nq
                    call vrmhd_section_basis(op, t, q, x, weight, basis, status)
                    if (status /= 0) return
                    do j = 1, op%scalar_fs%nl0
                        do i = 1, op%nlocal
                            if (op%row(i) /= 3) cycle
                            call vrmhd_pairing(basis(5, i), basis(6, i), basis(7, i), &
                                basis(13, j), basis(14, j), basis(15, j), value)
                            block(i, j) = block(i, j) + weight*value
                        end do
                    end do
                end do
                do j = 1, op%scalar_fs%nl0
                    col = index(op%scalar_fs%sp(0)%dof(j, t))
                    if (col == 0) cycle
                    do i = 1, op%nlocal
                        if (op%row(i) /= 3) cycle
                        row = op%map(i, t)
                        pos = pos + 1
                        rows(pos) = row; cols(pos) = col; values(pos) = block(i, j)
                    end do
                end do
            end do
            call csc_from_triplet(op%ndof, steady%ngauge, rows(:pos), cols(:pos), &
                values(:pos), steady%gradient, st)
            status = 1
            if (.not. status_ok(st)) return
            call csc_transpose(steady%gradient, steady%transpose, st)
            if (status_ok(st)) status = 0
        end subroutine

        subroutine vrmhd_stationary_residual(op, steady, p, u, du, boundary, source, &
                load, r, jv, status)
            type(vrmhd_section_t), intent(in) :: op
            type(vrmhd_stationary_t), intent(in) :: steady
            type(vrmhd_parameters_t), intent(in) :: p
            real(dp), intent(in) :: u(steady%ndof), du(steady%ndof), boundary(op%ndof)
            real(dp), intent(in) :: source(6, op%fs%nq, op%fs%mm%topology%n_triangles)
            real(dp), intent(in) :: load(op%ndof)
            real(dp), intent(out) :: r(steady%ndof), jv(steady%ndof)
            integer, intent(out) :: status
            real(dp) :: trial(op%ndof), tangent(op%ndof), zero(op%ndof)
            integer :: n
            n = op%ndof
            trial = u(:n); tangent = du(:n); zero = 0
            where (op%fixed)
                trial = boundary
                tangent = 0
            end where
            call vrmhd_section_residual(op, p, trial, zero, tangent, zero, source, &
                r(:n), jv(:n), status)
            if (status /= 0) return
            r(:n) = r(:n) - load + csc_matvec(steady%gradient, u(n + 1:))
            jv(:n) = jv(:n) + csc_matvec(steady%gradient, du(n + 1:))
            r(n + 1:) = csc_matvec(steady%transpose, trial)
            jv(n + 1:) = csc_matvec(steady%transpose, tangent)
            where (op%fixed)
                r(:n) = u(:n) - boundary
                jv(:n) = du(:n)
            end where
        end subroutine

        subroutine vrmhd_stationary_jacobian(op, steady, p, u, boundary, source, &
                jac, status)
            type(vrmhd_section_t), intent(in) :: op
            type(vrmhd_stationary_t), intent(in) :: steady
            type(vrmhd_parameters_t), intent(in) :: p
            real(dp), intent(in) :: u(steady%ndof), boundary(op%ndof)
            real(dp), intent(in) :: source(6, op%fs%nq, op%fs%mm%topology%n_triangles)
            type(csc_t), intent(out) :: jac
            integer, intent(out) :: status
            type(csc_t) :: volume
            type(fortsparse_status_t) :: st
            integer, allocatable :: rows(:), cols(:)
            real(dp), allocatable :: values(:)
            real(dp) :: trial(op%ndof), zero(op%ndof)
            integer :: i, j, k, pos, capacity, n
            n = op%ndof; trial = u(:n); zero = 0
            where (op%fixed) trial = boundary
            call vrmhd_section_jacobian(op, p, trial, zero, source, 0.0_dp, volume, &
                status)
            if (status /= 0) return
            capacity = size(volume%val) + 2*size(steady%gradient%val) + count(op%fixed)
            allocate (rows(capacity), cols(capacity), values(capacity))
            pos = 0
            do j = 1, n
                if (op%fixed(j)) then
                    pos = pos + 1; rows(pos) = j; cols(pos) = j; values(pos) = 1
                    cycle
                end if
                do k = volume%col_ptr(j), volume%col_ptr(j + 1) - 1
                    i = volume%row_idx(k)
                    if (op%fixed(i)) cycle
                    pos = pos + 1
                    rows(pos) = i; cols(pos) = j; values(pos) = volume%val(k)
                end do
            end do
            do j = 1, steady%ngauge
                do k = steady%gradient%col_ptr(j), steady%gradient%col_ptr(j + 1) - 1
                    i = steady%gradient%row_idx(k)
                    if (op%fixed(i)) cycle
                    pos = pos + 1
                    rows(pos) = i; cols(pos) = n + j
                    values(pos) = steady%gradient%val(k)
                    pos = pos + 1
                    rows(pos) = n + j; cols(pos) = i
                    values(pos) = steady%gradient%val(k)
                end do
            end do
            call csc_from_triplet(steady%ndof, steady%ndof, rows(:pos), cols(:pos), &
                values(:pos), jac, st)
            status = 1
            if (status_ok(st)) status = 0
        end subroutine
    end module kin6d_vrmhd_stationary
