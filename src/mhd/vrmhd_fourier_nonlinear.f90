module kin6d_vrmhd_fourier_nonlinear
    use, intrinsic :: iso_fortran_env, only: dp => real64, wp => real128
    use kin6d_vrmhd_parameters, only: vrmhd_parameters_t
    use kin6d_vrmhd_fourier, only: vrmhd_fourier_t, vrmhd_fourier_data
    use kin6d_vrmhd_section, only: vrmhd_section_test
    use kin6d_vrmhd_point, only: vrmhd_point, vrmhd_point_matrix
    use fortsparse, only: csc_t, csc_from_triplet, fortsparse_status_t, status_ok
    implicit none
    private
    public :: vrmhd_fourier_nonlinear_residual, vrmhd_fourier_nonlinear_jacobian
    public :: vrmhd_fourier_nonlinear_basis
contains
    ! theta is the dimensionless orbit angle; the physical derivative is n/r0.
    subroutine vrmhd_fourier_nonlinear_basis(op, t, q, theta, x, weight, basis, &
            map, status)
        type(vrmhd_fourier_t), intent(in) :: op
        integer, intent(in) :: t, q
        real(dp), intent(in) :: theta
        real(dp), intent(out) :: x(2), weight
        real(dp), intent(out) :: basis(36, op%section%nlocal*sum(op%copies))
        integer, intent(out) :: map(op%section%nlocal*sum(op%copies)), status
        real(dp) :: base(36, op%section%nlocal), orbit(36, op%section%nlocal)
        real(dp) :: phi, derivative, angle, kappa
        integer :: b, c, lo, nl, ns
        nl = op%section%nlocal; ns = op%section%ndof; lo = 0
        do b = 1, size(op%harmonic)
            do c = 1, op%copies(b)
                map(lo+1:lo+nl) = op%offset(b) + (c-1)*ns + op%section%map(:, t)
                lo = lo + nl
            end do
        end do
        call vrmhd_fourier_data(op, t, q, x, weight, base, orbit, status)
        if (status /= 0) return
        lo = 0
        do b = 1, size(op%harmonic)
            angle = op%harmonic(b)*theta
            kappa = op%harmonic(b)/op%orbit_scale
            do c = 1, op%copies(b)
                if (op%harmonic(b) == 0) then
                    phi = 1; derivative = 0
                else if (c == 1) then
                    phi = cos(angle); derivative = -kappa*sin(angle)
                else
                    phi = sin(angle); derivative = kappa*cos(angle)
                end if
                basis(:, lo+1:lo+nl) = phi*base + derivative*orbit
                lo = lo + nl
            end do
        end do
    end subroutine

    ! Orbit averaging uses the same angles for every test and nonlinear term.
    ! Nonpolynomial transport requires an angular-refinement check by the caller.
    subroutine vrmhd_fourier_nonlinear_residual(op, p, u, ud, du, dud, source, &
            residual, tangent, status)
        type(vrmhd_fourier_t), intent(in) :: op
        type(vrmhd_parameters_t), intent(in) :: p
        real(dp), intent(in) :: u(op%ndof), ud(op%ndof), du(op%ndof), dud(op%ndof)
        real(dp), intent(in) :: source(:, :, :, :)
        real(dp), intent(out) :: residual(op%ndof), tangent(op%ndof)
        integer, intent(out) :: status
        real(dp) :: basis(36, op%section%nlocal*sum(op%copies))
        real(dp) :: local_r(op%section%nlocal*sum(op%copies))
        real(dp) :: local_j(op%section%nlocal*sum(op%copies))
        real(dp) :: x(2), weight, v(36), vd(12), dv(36), dvd(12), test(19)
        real(dp) :: value, derivative, theta
        integer :: map(op%section%nlocal*sum(op%copies))
        integer :: t, q, a, i, row, st, fail, na, nl, section_column
        residual = 0; tangent = 0
        status = 2
        if (.not. valid_source(op, source)) return
        na = size(source, 2); nl = size(map); fail = 0
        !$omp parallel do default(shared) schedule(static) &
        !$omp private(t,q,a,i,row,st,basis,map,local_r,local_j,x,weight,v,vd, &
        !$omp dv,dvd,test,value,derivative,theta,section_column) reduction(max:fail)
        do t = 1, op%section%fs%mm%topology%n_triangles
            local_r = 0; local_j = 0
            do q = 1, op%section%fs%nq
                do a = 1, na
                    theta = 2*acos(-1.0_dp)*real(a-1, dp)/na
                    call vrmhd_fourier_nonlinear_basis(op, t, q, theta, x, &
                        weight, basis, map, st)
                    if (st /= 0) then
                        fail = 1; cycle
                    end if
                    v = matmul(basis, u(map)); vd = matmul(basis(:12, :), ud(map))
                    dv = matmul(basis, du(map)); dvd = matmul(basis(:12, :), dud(map))
                    do i = 1, nl
                        section_column = modulo(i-1, op%section%nlocal) + 1
                        row = op%section%row(section_column)
                        call vrmhd_section_test(row, basis(:, i), test)
                        call vrmhd_point(p, row, v(:33), dv(:33), vd, dvd, test, &
                            source(:, a, q, t), value, derivative)
                        local_r(i) = local_r(i) + weight/na*value
                        local_j(i) = local_j(i) + weight/na*derivative
                    end do
                end do
            end do
            do i = 1, nl
                !$omp atomic update
                residual(map(i)) = residual(map(i)) + local_r(i)
                !$omp atomic update
                tangent(map(i)) = tangent(map(i)) + local_j(i)
            end do
        end do
        !$omp end parallel do
        status = fail
    end subroutine

    subroutine vrmhd_fourier_nonlinear_jacobian(op, p, u, ud, source, shift, jac, &
            status)
        type(vrmhd_fourier_t), intent(in) :: op
        type(vrmhd_parameters_t), intent(in) :: p
        real(dp), intent(in) :: u(op%ndof), ud(op%ndof), shift
        real(dp), intent(in) :: source(:, :, :, :)
        type(csc_t), intent(out) :: jac
        integer, intent(out) :: status
        real(dp) :: basis(36, op%section%nlocal*sum(op%copies))
        real(dp) :: trial(45, op%section%nlocal*sum(op%copies))
        real(dp) :: action(19, op%section%nlocal*sum(op%copies), 6)
        real(wp), allocatable :: local(:, :)
        real(dp) :: x(2), weight, v(36), vd(12), test(19), kmat(19, 45, 6), theta
        integer :: map(op%section%nlocal*sum(op%copies))
        integer, allocatable :: rows(:), columns(:)
        real(dp), allocatable :: values(:)
        type(fortsparse_status_t) :: sparse_status
        integer :: t, q, a, i, j, row, st, fail, na, nl, nt, pos, section_column
        status = 2
        if (.not. valid_source(op, source)) return
        na = size(source, 2); nl = size(map)
        nt = op%section%fs%mm%topology%n_triangles
        allocate (rows(nl*nl*nt), columns(nl*nl*nt), values(nl*nl*nt))
        fail = 0
        !$omp parallel default(shared) &
        !$omp private(t,q,a,i,j,row,st,basis,trial,action,local,x,weight,v,vd,test, &
        !$omp kmat,theta,map,pos,section_column) reduction(max:fail)
        ! Harmonic pairs make this block too large for a worker thread stack.
        allocate (local(nl, nl))
        !$omp do schedule(static)
        do t = 1, nt
            local = 0
            do q = 1, op%section%fs%nq
                do a = 1, na
                    theta = 2*acos(-1.0_dp)*real(a-1, dp)/na
                    call vrmhd_fourier_nonlinear_basis(op, t, q, theta, x, &
                        weight, basis, map, st)
                    if (st /= 0) then
                        fail = 1; cycle
                    end if
                    v = matmul(basis, u(map)); vd = matmul(basis(:12, :), ud(map))
                    call vrmhd_point_matrix(p, v(:33), vd, source(:, a, q, t), &
                        shift, kmat)
                    trial(:33, :) = basis(:33, :); trial(34:, :) = basis(:12, :)
                    do row = 1, 6
                        action(:, :, row) = matmul(kmat(:, :, row), trial)
                    end do
                    do i = 1, nl
                        section_column = modulo(i-1, op%section%nlocal) + 1
                        row = op%section%row(section_column)
                        call vrmhd_section_test(row, basis(:, i), test)
                        local(i, :) = local(i, :) + real(weight/na, wp) &
                            *real(matmul(test, action(:, :, row)), wp)
                    end do
                end do
            end do
            pos = (t-1)*nl*nl
            do j = 1, nl
                do i = 1, nl
                    pos = pos + 1
                    rows(pos) = map(i); columns(pos) = map(j)
                    values(pos) = real(local(i, j), dp)
                end do
            end do
        end do
        !$omp end do
        deallocate (local)
        !$omp end parallel
        status = fail
        if (status /= 0) return
        call csc_from_triplet(op%ndof, op%ndof, rows, columns, values, jac, &
            sparse_status)
        if (.not. status_ok(sparse_status)) status = 1
    end subroutine

    logical function valid_source(op, source) result(valid)
        type(vrmhd_fourier_t), intent(in) :: op
        real(dp), intent(in) :: source(:, :, :, :)
        valid = size(source, 1) == 6 .and. size(source, 2) >= 4*maxval(op%harmonic)+1 &
            .and. size(source, 3) == op%section%fs%nq &
            .and. size(source, 4) == op%section%fs%mm%topology%n_triangles &
            .and. any(op%harmonic == 0)
    end function
end module kin6d_vrmhd_fourier_nonlinear
