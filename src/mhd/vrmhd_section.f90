module kin6d_vrmhd_section
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_gs_symmetry, only: gs_symmetry_t, GS_ROTATION
    use kin6d_vrmhd_parameters, only: vrmhd_parameters_t
    use kin6d_vrmhd_point, only: vrmhd_point
    use kin6d_vrmhd_section_leaves, only: vrmhd_grad_rotation, vrmhd_grad_screw, &
        vrmhd_velocity_rotation, vrmhd_velocity_screw, &
        vrmhd_curl_rotation, vrmhd_curl_screw
    use fortfem_meridian_curved_mesh, only: meridian_mesh_t
    use fortfem_meridian_feec_n, only: meridian_feec_t, meridian_feec_build, &
        meridian_qpoint, meridian_qgrad, meridian_point_basis_t, &
        meridian_point_basis
    use fortsparse, only: csc_t, csc_from_triplet, csc_transpose, csc_matvec, &
        fortsparse_status_t, status_ok
    implicit none
    private
    type, public :: vrmhd_section_t
        type(meridian_feec_t) :: fs, scalar_fs
        type(gs_symmetry_t) :: symmetry
        integer :: harmonic_count = 1, ndof = 0, nlocal = 0
        integer, allocatable :: map(:, :), row(:)
        logical, allocatable :: fixed(:)
    end type
    public :: vrmhd_section_build, vrmhd_section_basis, vrmhd_section_residual
    public :: vrmhd_section_jacobian, vrmhd_section_vjp, vrmhd_section_mass
    public :: vrmhd_section_trace
contains
    subroutine vrmhd_section_build(mm, symmetry, degree, harmonic_count, op, status, &
            scalar_degree, study_pairing)
        type(meridian_mesh_t), intent(in) :: mm
        type(gs_symmetry_t), intent(in) :: symmetry
        integer, intent(in) :: degree, harmonic_count
        type(vrmhd_section_t), intent(out) :: op
        integer, intent(out) :: status
        logical, intent(in), optional :: study_pairing
        logical :: study
        integer, intent(in), optional :: scalar_degree
        integer :: i, t, j, b, n0, n1, l0, l1, offset, lo, ks, qdegree
        status = 1
        if (.not. symmetry%valid() .or. harmonic_count /= 1) return
        ks = degree + 1
        if (present(scalar_degree)) ks = scalar_degree
        study = .false.
        if (present(study_pairing)) study = study_pairing
        if (degree < 1 .or. ks < 1) return
        if (study) then
            if (ks < degree - 1) return
        else
            if (ks < degree + 1) return
        end if
        op%symmetry = symmetry; op%harmonic_count = harmonic_count
        qdegree = 4*max(ks, degree) + 2*mm%geo_order + 2
        if (symmetry%kind == GS_ROTATION) then
            call meridian_feec_build(mm, degree, op%fs, status, &
                quad_degree=qdegree)
            if (status /= 0) return
            call meridian_feec_build(mm, ks, op%scalar_fs, status, &
                quad_degree=qdegree)
        else
            call meridian_feec_build(mm, degree, op%fs, status, &
                quad_degree=qdegree, cylinder_r0=1.0_dp)
            if (status /= 0) return
            call meridian_feec_build(mm, ks, op%scalar_fs, status, &
                quad_degree=qdegree, cylinder_r0=1.0_dp)
        end if
        if (status /= 0) return
        call op%fs%mm%topology%build_edge_connectivity()
        n0 = op%scalar_fs%n0; n1 = op%fs%sp(1)%ndof
        l0 = op%scalar_fs%nl0; l1 = op%fs%sp(1)%nloc
        op%ndof = harmonic_count*(6*n0 + 2*n1)
        op%nlocal = harmonic_count*(6*l0 + 2*l1)
        allocate (op%map(op%nlocal, mm%topology%n_triangles), &
            op%row(op%nlocal), op%fixed(op%ndof))
        op%fixed = .false.
        do b = 1, harmonic_count
            offset = (b - 1)*(6*n0 + 2*n1); lo = (b - 1)*(6*l0 + 2*l1)
            do i = 1, 6
                do t = 1, mm%topology%n_triangles
                    op%map(lo + (i - 1)*l0 + 1:lo + i*l0, t) = &
                        offset + (i - 1)*n0 + op%scalar_fs%sp(0)%dof(:, t)
                end do
                op%row(lo + (i - 1)*l0 + 1:lo + i*l0) = 2
                if (i == 1) op%row(lo + 1:lo + l0) = 1
                if (i >= 5) op%row(lo + (i - 1)*l0 + 1:lo + i*l0) = i
                op%fixed(offset + (i - 1)*n0 + 1:offset + i*n0) = op%scalar_fs%sp(0)%bnd
            end do
            do j = 1, 2
                i = lo + 6*l0 + (j - 1)*l1
                do t = 1, mm%topology%n_triangles
                    op%map(i + 1:i + l1, &
                        t) = offset + 6*n0 + (j - 1)*n1 + op%fs%sp(1)%dof(:, t)
                end do
                op%row(i + 1:i + l1) = j + 2
                ! A and the free current use the same zero-tangential tests.
                op%fixed(offset + 6*n0 + (j - 1)*n1 + 1: &
                    offset + 6*n0 + j*n1) = op%fs%sp(1)%bnd
            end do
        end do
    end subroutine

    subroutine vrmhd_section_basis(op, t, q, x, weight, basis, status)
        type(vrmhd_section_t), intent(in) :: op
        integer, intent(in) :: t, q
        real(dp), intent(out) :: x(2), weight, basis(36, op%nlocal)
        integer, intent(out) :: status
        real(dp) :: val(op%scalar_fs%nl0), grad(2, op%scalar_fs%nl0)
        real(dp) :: aval(op%fs%nl0), agrad(2, op%fs%nl0), g(3), v(3), vx(3), vy(3)
        real(dp) :: gv(3, 3), phi(3, op%fs%sp(1)%nloc), psi(3, op%fs%sp(2)%nloc)
        real(dp) :: curl(3), av(3), bv(3), a, ax, ay, c, scale
        integer :: i, j, l, col, p, r, local, s, b, offset
        call meridian_qgrad(op%scalar_fs, t, q, x, weight, val, grad, status)
        if (status /= 0) return
        call meridian_qgrad(op%fs, t, q, x, weight, aval, agrad, status)
        if (status /= 0) return
        call meridian_qpoint(op%fs, 1, t, q, x, weight, phi, status)
        if (status /= 0) return
        call meridian_qpoint(op%fs, 2, t, q, x, weight, psi, status)
        if (status /= 0) return
        basis = 0
        do b = 1, op%harmonic_count
            offset = (b - 1)*(6*op%scalar_fs%nl0 + 2*op%fs%sp(1)%nloc)
            do i = 1, op%scalar_fs%nl0
                call scalar_gradient(op, x, grad(:, i), g)
                basis(1, offset + i) = val(i); basis(13:15, offset + i) = g
                do j = 1, 3
                    local = offset + j*op%scalar_fs%nl0 + i
                    v = 0; vx = 0; vy = 0
                    v(j) = val(i); vx(j) = grad(1, i); vy(j) = grad(2, i)
                    call vector_gradient(op, x, v, vx, vy, gv)
                    basis(2:4, local) = v; basis(16:24, local) = reshape(gv, [9])
                end do
                do j = 1, 2
                    local = offset + (3 + j)*op%scalar_fs%nl0 + i
                    basis(10 + j, local) = val(i)
                    basis(25 + 3*j:27 + 3*j, local) = g
                end do
            end do
            do i = 1, op%fs%sp(1)%nloc
                col = op%fs%sp(1)%dof(i, t); curl = 0
                do p = op%fs%curl%a%col_ptr(col), op%fs%curl%a%col_ptr(col + 1) - 1
                    r = op%fs%curl%a%row_idx(p)
                    do l = 1, op%fs%sp(2)%nloc
                        if (op%fs%sp(2)%dof(l, t) /= r) cycle
                        curl = curl + op%fs%curl%a%val(p)*psi(:, l) &
                            *real(op%fs%sp(2)%sgn(l, t), dp)
                    end do
                end do
                scale = real(op%fs%sp(1)%sgn(i, t), dp)
                a = 0; ax = 0; ay = 0
                if (i > op%fs%nln) then
                    j = i - op%fs%nln
                    a = aval(j); ax = agrad(1, j); ay = agrad(2, j)
                end if
                c = -curl(2)
                call vector_curl(op, x, scale*phi(1, i), scale*phi(3, i), c, a, &
                    ax, ay, av, bv)
                do j = 1, 2
                    local = offset + 6*op%scalar_fs%nl0 + (j - 1)*op%fs%sp(1)%nloc + i
                    s = 2 + 3*j
                    basis(s:s + 2, local) = av
                    if (j == 1) basis(25:27, local) = bv
                    if (j == 2) basis(34:36, local) = bv
                end do
            end do
        end do
    end subroutine

    subroutine scalar_gradient(op, x, f, g)
        type(vrmhd_section_t), intent(in) :: op
        real(dp), intent(in) :: x(2), f(2)
        real(dp), intent(out) :: g(3)
        if (op%symmetry%kind == GS_ROTATION) then
            call vrmhd_grad_rotation(x(1), x(2), op%symmetry%k, f(1), f(2), g(1), &
                g(2), g(3))
        else
            call vrmhd_grad_screw(x(1), x(2), op%symmetry%k, f(1), f(2), g(1), &
                g(2), g(3))
        end if
    end subroutine

    subroutine vrmhd_section_trace(op, pb, t, ref, x, trace, status)
        type(vrmhd_section_t), intent(in) :: op
        type(meridian_point_basis_t), intent(in) :: pb
        integer, intent(in) :: t
        real(dp), intent(in) :: ref(2)
        real(dp), intent(out) :: x(2), trace(3, op%fs%sp(1)%nloc)
        integer, intent(out) :: status
        real(dp) :: phi(3, op%fs%sp(1)%nloc), a, unused(3), radius, scale
        integer :: i
        call meridian_point_basis(op%fs, pb, 1, t, ref(1), ref(2), x, phi, status)
        if (status /= 0) return
        radius = 1
        if (op%symmetry%kind == GS_ROTATION) radius = x(1)
        do i = 1, op%fs%sp(1)%nloc
            scale = real(op%fs%sp(1)%sgn(i, t), dp)
            a = radius*phi(2, i)
            call vector_curl(op, x, scale*phi(1, i), scale*phi(3, i), 0.0_dp, &
                scale*a, 0.0_dp, 0.0_dp, trace(:, i), unused)
        end do
    end subroutine

    subroutine vector_gradient(op, x, v, vx, vy, g)
        type(vrmhd_section_t), intent(in) :: op
        real(dp), intent(in) :: x(2), v(3), vx(3), vy(3)
        real(dp), intent(out) :: g(3, 3)
        if (op%symmetry%kind == GS_ROTATION) then
            call vrmhd_velocity_rotation(x(1), x(2), op%symmetry%k, &
                v(1), v(2), v(3), vx(1), vx(2), vx(3), &
                vy(1), vy(2), vy(3), &
                g(1, 1), g(2, 1), g(3, 1), &
                g(1, 2), g(2, 2), g(3, 2), &
                g(1, 3), g(2, 3), g(3, 3))
        else
            call vrmhd_velocity_screw(x(1), x(2), op%symmetry%k, &
                v(1), v(2), v(3), vx(1), vx(2), vx(3), &
                vy(1), vy(2), vy(3), &
                g(1, 1), g(2, 1), g(3, 1), g(1, 2), g(2, 2), &
                g(3, 2), g(1, 3), g(2, 3), g(3, 3))
        end if
    end subroutine

    subroutine vector_curl(op, x, u1, u2, c, a, ax, ay, v, b)
        type(vrmhd_section_t), intent(in) :: op
        real(dp), intent(in) :: x(2), u1, u2, c, a, ax, ay
        real(dp), intent(out) :: v(3), b(3)
        if (op%symmetry%kind == GS_ROTATION) then
            call vrmhd_curl_rotation(x(1), x(2), op%symmetry%k, u1, u2, c, a, ax, &
                ay, &
                v(1), v(2), v(3), b(1), b(2), b(3))
        else
            call vrmhd_curl_screw(x(1), x(2), op%symmetry%k, u1, u2, c, a, ax, ay, &
                v(1), v(2), v(3), b(1), b(2), b(3))
        end if
    end subroutine

    pure subroutine test_data(row, basis, test)
        integer, intent(in) :: row
        real(dp), intent(in) :: basis(36)
        real(dp), intent(out) :: test(19)
        test = 0
        select case (row)
        case (1)
            test(1) = basis(1); test(2:4) = basis(13:15)
        case (2)
            test(5:7) = basis(2:4); test(8:16) = basis(16:24)
        case (3)
            test(5:7) = basis(5:7); test(17:19) = basis(25:27)
        case (4)
            test(5:7) = basis(8:10); test(17:19) = basis(34:36)
        case (5, 6)
            test(1) = basis(6 + row)
            test(2:4) = basis(13 + 3*row:15 + 3*row)
        end select
    end subroutine

    subroutine vrmhd_section_residual(op, p, u, ud, du, dud, source, r, jv, status)
        type(vrmhd_section_t), intent(in) :: op
        type(vrmhd_parameters_t), intent(in) :: p
        real(dp), intent(in) :: u(op%ndof), ud(op%ndof), du(op%ndof), dud(op%ndof)
        real(dp), intent(in) :: source(6, op%fs%nq, op%fs%mm%topology%n_triangles)
        real(dp), intent(out) :: r(op%ndof), jv(op%ndof)
        integer, intent(out) :: status
        real(dp) :: local(op%nlocal), dlocal(op%nlocal), dotlocal(op%nlocal)
        real(dp) :: ddotlocal(op%nlocal), lr(op%nlocal), lj(op%nlocal)
        real(dp) :: basis(36, op%nlocal), x(2), weight, v(36), dv(36), vd(12), dvd(12)
        real(dp) :: test(19), value, tangent
        integer :: t, q, i, st, fail
        r = 0; jv = 0; fail = 0
        !$omp parallel do default(shared) schedule(static) &
        !$omp private(t,q,i,st,local,dlocal,dotlocal,ddotlocal,lr,lj, &
        !$omp basis,x,weight,v,dv,vd,dvd,test,value,tangent) reduction(max:fail)
        do t = 1, op%fs%mm%topology%n_triangles
            local = u(op%map(:, t)); dlocal = du(op%map(:, t))
            dotlocal = ud(op%map(:, t)); ddotlocal = dud(op%map(:, t))
            lr = 0; lj = 0
            do q = 1, op%fs%nq
                call vrmhd_section_basis(op, t, q, x, weight, basis, st)
                if (st /= 0) then
                    fail = 1; cycle
                end if
                v = matmul(basis, local); dv = matmul(basis, dlocal)
                vd = matmul(basis(:12, :), dotlocal)
                dvd = matmul(basis(:12, :), ddotlocal)
                do i = 1, op%nlocal
                    call test_data(op%row(i), basis(:, i), test)
                    call vrmhd_point(p, op%row(i), v(:33), dv(:33), vd, dvd, &
                        test, source(:, q, t), value, tangent)
                    lr(i) = lr(i) + weight*value; lj(i) = lj(i) + weight*tangent
                end do
            end do
            do i = 1, op%nlocal
                !$omp atomic update
                r(op%map(i, t)) = r(op%map(i, t)) + lr(i)
                !$omp atomic update
                jv(op%map(i, t)) = jv(op%map(i, t)) + lj(i)
            end do
        end do
        !$omp end parallel do
        status = fail
    end subroutine

    subroutine vrmhd_section_mass(op, p, u, ud, mass, status)
        type(vrmhd_section_t), intent(in) :: op
        type(vrmhd_parameters_t), intent(in) :: p
        real(dp), intent(in) :: u(op%ndof), ud(op%ndof)
        real(dp), intent(out) :: mass(op%ndof)
        integer, intent(out) :: status
        real(dp) :: zero(op%ndof), unused(op%ndof)
        real(dp) :: source(6, op%fs%nq, op%fs%mm%topology%n_triangles)
        zero = 0; source = 0
        call vrmhd_section_residual(op, p, u, zero, zero, ud, source, unused, &
            mass, status)
    end subroutine

    subroutine vrmhd_section_jacobian(op, p, u, ud, source, shift, jac, status)
        type(vrmhd_section_t), intent(in) :: op
        type(vrmhd_parameters_t), intent(in) :: p
        real(dp), intent(in) :: u(op%ndof), ud(op%ndof), shift
        real(dp), intent(in) :: source(6, op%fs%nq, op%fs%mm%topology%n_triangles)
        type(csc_t), intent(out) :: jac
        integer, intent(out) :: status
        real(dp) :: basis(36, op%nlocal), x(2), weight, v(36), vd(12), test(19), &
            value, dvalue
        real(dp) :: block(op%nlocal, op%nlocal)
        integer, allocatable :: rows(:), cols(:)
        real(dp), allocatable :: vals(:)
        type(fortsparse_status_t) :: sparse_status
        integer :: t, q, i, j, pos, st, fail
        allocate (rows(op%nlocal**2*op%fs%mm%topology%n_triangles), &
            cols(op%nlocal**2*op%fs%mm%topology%n_triangles), &
            vals(op%nlocal**2*op%fs%mm%topology%n_triangles))
        fail = 0
        !$omp parallel do default(shared) schedule(static) &
        !$omp private(t,q,i,j,pos,st,basis,x,weight,v,vd,test,value,dvalue,block) &
        !$omp reduction(max:fail)
        do t = 1, op%fs%mm%topology%n_triangles
            block = 0
                do q = 1, op%fs%nq
                    call vrmhd_section_basis(op, t, q, x, weight, basis, st)
                    if (st /= 0) then
                        fail = 1; cycle
                    end if
                    v = matmul(basis, u(op%map(:, t)))
                    vd = matmul(basis(:12, :), ud(op%map(:, t)))
                    do i = 1, op%nlocal
                        call test_data(op%row(i), basis(:, i), test)
                        do j = 1, op%nlocal
                            call vrmhd_point(p, op%row(i), v(:33), basis(:33, j), vd, &
                                shift*basis(:12, j), test, source(:, q, t), &
                                value, dvalue)
                            block(i, j) = block(i, j) + weight*dvalue
                        end do
                    end do
                end do
                pos = (t - 1)*op%nlocal**2
                do j = 1, op%nlocal
                    do i = 1, op%nlocal
                        pos = pos + 1
                        rows(pos) = op%map(i, t)
                        cols(pos) = op%map(j, t)
                        vals(pos) = block(i, j)
                    end do
                end do
            end do
            !$omp end parallel do
            status = fail
            if (status /= 0) return
            call csc_from_triplet(op%ndof, op%ndof, rows, cols, vals, jac, sparse_status)
            if (.not. status_ok(sparse_status)) status = 1
        end subroutine

        subroutine vrmhd_section_vjp(jac, dual, result, status)
            type(csc_t), intent(in) :: jac
            real(dp), intent(in) :: dual(jac%nrow)
            real(dp), intent(out) :: result(jac%ncol)
            integer, intent(out) :: status
            type(csc_t) :: transpose
            type(fortsparse_status_t) :: st
            call csc_transpose(jac, transpose, st)
            status = 1
            result = 0
            if (.not. status_ok(st)) return
            result = csc_matvec(transpose, dual)
            status = 0
        end subroutine
    end module kin6d_vrmhd_section
