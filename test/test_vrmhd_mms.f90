program test_vrmhd_mms
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use milca_gs_symmetry, only: gs_symmetry_t, GS_ROTATION, GS_SCREW, &
        GS_CHART_RZ, GS_CHART_CARTESIAN
    use kin6d_vrmhd_parameters, only: vrmhd_parameters_t
    use kin6d_vrmhd_section
    use vrmhd_mms_leaves, only: mms_rotation, mms_screw
    use fortfem_mesh_2d, only: mesh_2d_t
    use fortfem_meridian_curved_mesh, only: meridian_mesh_t, meridian_mesh_straight
    use fortsparse, only: csc_t, csc_from_triplet, sparse_solve_once, &
        fortsparse_status_t, status_ok
    implicit none
    real(dp) :: errors(8, 3), defect(3), rates(8)
    integer :: kind, level, degree
    real(dp) :: target(8)
    character(16) :: mode
    call get_command_argument(1, mode)
    degree = 0
    if (mode == 'mixed1') degree = 1
    if (mode == 'mixed2') degree = 2
    target = 1
    if (degree > 0) target = [degree, degree + 1, degree + 1, degree + 1, &
        degree + 1, degree + 1, degree, degree]
    do kind = 1, 3
        do level = 1, 3
            call manufactured(kind, 2**(level + 1), errors(:, level), defect(level))
            print '(a,i0,i0,8es12.4)', 'natural-norm errors ', kind, level, &
                errors(:, level)
        end do
        do level = 2, 3
            rates = log(errors(:, level - 1)/errors(:, level))/log(2.0_dp)
            print '(a,i0,a,i0,8f7.3)', 'natural-norm rates symmetry ', kind, &
                ' level ', level, rates
            call require(all(rates > target - 0.2_dp), 'natural-norm design orders')
            print '(a,es12.4)', 'scaled strong-source weak defect: ', defect(level)
            call require(defect(level) < 0.8_dp*defect(level - 1), &
                'strong-source residual consistency')
        end do
    end do
contains
    subroutine require(ok, label)
        logical, intent(in) :: ok
        character(*), intent(in) :: label
        if (.not. ok) then
            print '(a)', label
            error stop 'manufactured section consistency'
        end if
    end subroutine

    pure function norms(u) result(e)
        real(dp), intent(in) :: u(36)
        real(dp) :: e(8)
        integer :: i
        e(1) = u(1)**2 + sum(u(13:15)**2)
        do i = 2, 4
            e(i) = u(i)**2 + u(i + 14)**2 + u(i + 17)**2 + u(i + 20)**2
        end do
        e(5) = sum(u(5:7)**2) + sum(u(25:27)**2)
        e(6) = sum(u(8:10)**2) + sum(u(34:36)**2)
        e(7) = u(11)**2 + sum(u(28:30)**2)
        e(8) = u(12)**2 + sum(u(31:33)**2)
    end function

    subroutine manufactured(kind, n, error, defect)
        integer, intent(in) :: kind, n
        real(dp), intent(out) :: error(8), defect
        type(vrmhd_section_t) :: op
        type(vrmhd_parameters_t) :: p
        type(mesh_2d_t) :: mesh
        type(meridian_mesh_t) :: mm
        type(gs_symmetry_t) :: symmetry
        type(csc_t) :: mass
        type(fortsparse_status_t) :: sparse_status
        real(dp), allocatable :: basis(:, :), rhs(:), u(:), zero(:), source(:, :, :)
        real(dp), allocatable :: load(:), r(:), jv(:), diagonal(:), vals(:), &
            block(:, :)
        integer, allocatable :: rows(:), cols(:)
        logical, allocatable :: free(:)
        real(dp) :: x(2), weight, exact(48), v(36), shift
        integer :: t, q, i, j, pos, status, nlocal, nt
        shift = 0
        symmetry = gs_symmetry_t(kind, GS_CHART_CARTESIAN, 0.0_dp)
        if (kind == GS_ROTATION) then
            shift = 2; symmetry%chart = GS_CHART_RZ
        end if
        if (kind == GS_SCREW) symmetry%k = 0.7_dp
        call mesh%create_rectangular(n + 1, n + 1, shift + 0.1_dp, shift + 1.1_dp, &
            0.2_dp, 1.2_dp)
        call meridian_mesh_straight(mesh, 1, mm, status)
        call require(status == 0, 'manufactured mesh')
        if (degree > 0) then
            call vrmhd_section_build(mm, symmetry, degree + 1, 1, op, status, &
                scalar_degree=degree, velocity_degree=degree + 1)
        else
            call vrmhd_section_build(mm, symmetry, 1, 1, op, status)
        end if
        call require(status == 0, 'manufactured section')
        p%mi = 1; p%mu0 = 1; p%eta0 = 0.2_dp; p%nu0 = 0.2_dp
        p%chi_parallel0 = 0.2_dp; p%chi_perpendicular = 0.05_dp
        p%d_parallel = 0.1_dp; p%d_perpendicular = 0.1_dp; p%epsilon_b = 0.1_dp
        nt = mm%topology%n_triangles; nlocal = op%nlocal
        allocate (basis(36, nlocal), rhs(op%ndof), u(op%ndof), zero(op%ndof), &
            source(6, op%fs%nq, nt), load(op%ndof), r(op%ndof), jv(op%ndof), &
            diagonal(op%ndof), free(op%ndof), block(nlocal, nlocal), &
            rows(nlocal*nlocal*nt), cols(nlocal*nlocal*nt), &
            vals(nlocal*nlocal*nt))
        rhs = 0; load = 0; diagonal = 0; zero = 0; source = 0
        do t = 1, nt
            block = 0
                do q = 1, op%fs%nq
                    call vrmhd_section_basis(op, t, q, x, weight, basis, status)
                    call require(status == 0, 'manufactured basis')
                    if (kind == GS_ROTATION) then
                        call mms_rotation(x(1), x(2), symmetry%k, exact)
                    else
                        call mms_screw(x(1), x(2), symmetry%k, exact)
                    end if
                    block = block + weight*matmul(transpose(basis), basis)
                        do i = 1, nlocal
                            j = op%map(i, t)
                            rhs(j) = rhs(j) + weight*dot_product(basis(:, i), exact(:36))
                            load(j) = load(j) + weight*dot_product(basis(:12, i), &
                                exact(37:48))
                            diagonal(j) = diagonal(j) + weight*sum(norms(basis(:, i)))
                        end do
                    end do
                    pos = (t - 1)*nlocal*nlocal
                    do j = 1, nlocal
                        do i = 1, nlocal
                            pos = pos + 1
                            rows(pos) = op%map(i, t)
                            cols(pos) = op%map(j, t)
                            vals(pos) = block(i, j)
                        end do
                    end do
                end do
                call csc_from_triplet(op%ndof, op%ndof, rows, cols, vals, mass, sparse_status)
                call require(status_ok(sparse_status), 'projection mass assembly')
                call sparse_solve_once(mass, rhs, u, sparse_status)
                call require(status_ok(sparse_status), 'independent natural-norm projection')
                error = 0
                do t = 1, nt
                    do q = 1, op%fs%nq
                        call vrmhd_section_basis(op, t, q, x, weight, basis, status)
                        if (kind == GS_ROTATION) then
                            call mms_rotation(x(1), x(2), symmetry%k, exact)
                        else
                            call mms_screw(x(1), x(2), symmetry%k, exact)
                        end if
                        v = matmul(basis, u(op%map(:, t)))
                        error = error + weight*norms(v - exact(:36))
                    end do
                end do
                error = sqrt(error)
                call vrmhd_section_residual(op, p, u, zero, zero, zero, source, r, jv, status)
                call require(status == 0, 'manufactured residual')
                free = .not. op%fixed
                free(op%offset(8) + 1:) = .not. op%fs%sp(1)%bnd
                defect = sqrt(sum((r - load)**2/diagonal, mask=free))
            end subroutine
        end program test_vrmhd_mms
