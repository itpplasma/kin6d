submodule(kin6d_nonnested_coupled) kin6d_nc_matrix_free
    ! Opt-in Fourier pseudo-spectral Jacobian action and mean-n-block solve.
    ! All constitutive derivatives come from the retained generated leaves.
    use fortnum_fft, only: fortnum_fft_plan_t, fortnum_fft_workspace_t, &
        fft_c2c_plan_init, fft_workspace_init, fft_c2c
    implicit none
contains
    module subroutine nc_jvp(field, problem, x, direction, output, status)
        ! Exact directional derivative, with the same phi quadrature as the
        ! dense reference. Fixed rows are zero and floating rows are summed.
        type(free_b_field_t), intent(in) :: field
        type(nc_problem_t), intent(in) :: problem
        real(dp), intent(in) :: x(:), direction(:)
        real(dp), allocatable, intent(out) :: output(:)
        integer, intent(out) :: status
        type(nc_ops_t) :: ops
        call operators(field, problem, ops, status)
        if (status /= 0) return
        call jvp_ops(field, problem, ops, x, direction, output, status)
    end subroutine

    module subroutine nc_vjp(field, problem, x, weight, output, status)
        ! Transpose of nc_jvp in the Euclidean coefficient pairing. Reverse
        ! the row reduction before testing and retain boundary input columns.
        type(free_b_field_t), intent(in) :: field
        type(nc_problem_t), intent(in) :: problem
        real(dp), intent(in) :: x(:), weight(:)
        real(dp), allocatable, intent(out) :: output(:)
        integer, intent(out) :: status
        type(nc_ops_t) :: ops
        type(nc_layout_t) :: lay
        type(csc_t) :: ct, st
        type(fortsparse_status_t) :: sparse_status
        type(fortnum_fft_plan_t) :: plan
        type(fortnum_fft_workspace_t) :: work
        real(dp), allocatable :: b(:), w(:), rb(:, :), rs(:, :), dual_b(:)
        integer, allocatable :: rep(:)
        integer :: t, nt, bk, j, failure, local_status, n0, n2
        status = 1; lay = nc_layout(field)
        if (size(x) /= lay%n .or. size(weight) /= lay%n) return
        if (.not. all(ieee_is_finite(x))) return
        if (.not. all(ieee_is_finite(weight))) return
        call operators(field, problem, ops, status)
        if (status /= 0) return
        status = 1
        call csc_transpose(ops%c, ct, sparse_status)
        if (.not. status_ok(sparse_status)) return
        call csc_transpose(ops%s, st, sparse_status)
        if (.not. status_ok(sparse_status)) return
        rep = representative(field, problem)
        w = weight(rep)
        where (nc_fixed(field, problem)) w = 0
        b = csc_matvec(ops%c, x(:lay%na))
        call fft_c2c_plan_init(plan, ops%nphi)
        call fft_workspace_init(work, plan)
        nt = field%fs%mm%topology%n_triangles
        n0 = field%fs%nl0; n2 = field%fs%sp(2)%nloc
        allocate (rb(n2*lay%nb, nt), rs(2*n0*lay%nb + 1, nt))
        failure = 0
        !$omp parallel do default(shared) private(t, local_status) &
        !$omp& firstprivate(work) schedule(static) reduction(max:failure)
        do t = 1, nt
            call element_vjp(field, problem, lay, ops%nphi, b, x, w, t, &
                rb(:, t), rs(:, t), plan, work, ops%nfp, local_status)
            failure = max(failure, local_status)
        end do
        !$omp end parallel do
        if (failure /= 0) return
        allocate (output(lay%n), dual_b(lay%nb2*lay%nb))
        output = 0; dual_b = 0
        output(:lay%na) = csc_matvec(st, w(:lay%na))/field%mu0 &
            + csc_matvec(ops%mg, w(lay%il + 1:lay%ip))
        output(lay%il + 1:lay%ip) = csc_matvec(ops%mgt, w(:lay%na))
        do t = 1, nt
            do bk = 0, lay%nb - 1
                do j = 1, n2
                    local_status = bk*lay%nb2 + field%fs%sp(2)%dof(j, t)
                    dual_b(local_status) = dual_b(local_status) &
                        + field%fs%sp(2)%sgn(j, t)*rb(bk*n2 + j, t)
                end do
                do j = 1, n0
                    local_status = bk*lay%n0 + field%fs%sp(0)%dof(j, t)
                    output(lay%ip + local_status) = output(lay%ip + local_status) &
                        + rs(bk*n0 + j, t)
                    output(lay%iphi + local_status) = &
                        output(lay%iphi + local_status) + rs((lay%nb + bk)*n0 + j, t)
                end do
            end do
            output(lay%iv) = output(lay%iv) + rs(size(rs, 1), t)
        end do
        output(:lay%na) = output(:lay%na) + csc_matvec(ct, dual_b)
        if (.not. all(ieee_is_finite(output))) return
        status = 0
    end subroutine

    subroutine element_vjp(field, problem, lay, nphi, b, x, weight, t, rb, rs, &
                           plan, work, nfp, status)
        type(free_b_field_t), intent(in) :: field
        type(nc_problem_t), intent(in) :: problem
        type(nc_layout_t), intent(in) :: lay
        integer, intent(in) :: nphi, t, nfp
        real(dp), intent(in) :: b(:), x(:), weight(:)
        real(dp), intent(out) :: rb(:), rs(:)
        type(fortnum_fft_plan_t), intent(in) :: plan
        type(fortnum_fft_workspace_t), intent(inout) :: work
        integer, intent(out) :: status
        real(dp) :: v1(3, field%fs%sp(1)%nloc), v2(3, field%fs%sp(2)%nloc)
        real(dp) :: val(field%fs%nl0), grad(2, field%fs%nl0), xq(2), wq
        real(dp) :: u(10, nphi), w(10, nphi), gu(10, nphi), gq(8)
        real(dp) :: am(3, lay%nb), wa(3, nphi), gm(9, lay%nb)
        real(dp) :: al(field%fs%sp(1)%nloc, lay%nb), dummy(2), phi, gpar
        real(dp) :: dpmean(lay%nb), dfmean(lay%nb)
        integer :: q, j, k, bk, n, n0, n2
        n0 = field%fs%nl0
        n2 = field%fs%sp(2)%nloc
        rb = 0; rs = 0
        do bk = 1, lay%nb
            al(:, bk) = weight((bk - 1)*lay%na1 + field%fs%sp(1)%dof(:, t)) &
                *field%fs%sp(1)%sgn(:, t)
        end do
        do q = 1, field%fs%nq
            call point_series(field, lay, t, q, b, x, u, xq, wq, &
                v1, v2, val, grad, plan, work, nfp, status)
            if (status /= 0) return
            call point_series(field, lay, t, q, b, weight, w, xq, wq, &
                v1, v2, val, grad, plan, work, nfp, status)
            if (status /= 0) return
            am = matmul(v1, al)
            do k = 1, 3
                call synthesize(am(k, :), wa(k, :), plan, work, nfp)
            end do
            u(2, :) = u(2, :) + field%f_vac/xq(1)
            w(7:9, :) = w(7:9, :) - wa
            w(8, :) = w(8, :) + weight(lay%iv)/xq(1)
            do j = 1, nphi
                status = 1
                if (dot_product(u(:3, j), u(:3, j)) <= tiny(1.0_dp)) return
                phi = 2*acos(-1.0_dp)*(j - 1)/(nphi*nfp); gpar = problem%g_par
                if (associated(problem%g_par_shape)) then
                    gpar = gpar*problem%g_par_shape(xq)
                else if (associated(problem%g_par_shape_3d)) then
                    gpar = gpar*problem%g_par_shape_3d(xq, phi)
                end if
                call coupled_current_vjp(u(1, j), u(2, j), u(3, j), u(4, j), &
                    u(5, j), u(6, j), u(7, j), u(8, j), u(9, j), u(10, j), &
                    gpar, problem%g_perp, w(7, j), w(8, j), w(9, j), &
                    gu(1, j), gu(2, j), gu(3, j), gu(4, j), gu(5, j), gu(6, j), &
                    gu(7, j), gu(8, j), gu(9, j), gu(10, j), dummy(1), dummy(2))
                call anisotropic_flux_vjp(u(1, j), u(2, j), u(3, j), u(4, j), &
                    u(5, j), u(6, j), problem%chi_par, problem%chi_perp, &
                    w(4, j), w(5, j), w(6, j), gq(1), gq(2), gq(3), &
                    gq(4), gq(5), gq(6), gq(7), gq(8))
                gu(:6, j) = gu(:6, j) + gq(:6)
            end do
            do k = 1, 9
                call project(gu(k, :), gm(k, :), plan, work, nfp)
            end do
            dpmean = 0; dfmean = 0
            do n = 1, lay%modes
                dpmean(2*n) = -n*gm(5, 2*n + 1)
                dpmean(2*n + 1) = n*gm(5, 2*n)
                dfmean(2*n) = -n*gm(8, 2*n + 1)
                dfmean(2*n + 1) = n*gm(8, 2*n)
            end do
            do bk = 1, lay%nb
                rb((bk - 1)*n2 + 1:bk*n2) = rb((bk - 1)*n2 + 1:bk*n2) &
                    + wq*matmul(gm(:3, bk), v2)
                rs((bk - 1)*n0 + 1:bk*n0) = rs((bk - 1)*n0 + 1:bk*n0) &
                    + wq*(grad(1, :)*gm(4, bk) + grad(2, :)*gm(6, bk) &
                          + val*dpmean(bk)/xq(1))
                rs((lay%nb + bk - 1)*n0 + 1:(lay%nb + bk)*n0) = &
                    rs((lay%nb + bk - 1)*n0 + 1:(lay%nb + bk)*n0) &
                    + wq*(grad(1, :)*gm(7, bk) + grad(2, :)*gm(9, bk) &
                          + val*dfmean(bk)/xq(1))
            end do
            rs(size(rs)) = rs(size(rs)) &
                + wq*sum(gu(10, :))/(nphi*2*acos(-1.0_dp)*xq(1))
        end do
        status = 0
    end subroutine

    subroutine synthesize(coeff, values, plan, work, nfp, derivative)
        real(dp), intent(in) :: coeff(:)
        real(dp), intent(out) :: values(:)
        type(fortnum_fft_plan_t), intent(in) :: plan
        type(fortnum_fft_workspace_t), intent(inout) :: work
        integer, intent(in) :: nfp
        logical, intent(in), optional :: derivative
        complex(dp) :: z(size(values)), a
        integer :: n
        z = 0
        z(1) = coeff(1)
        if (present(derivative)) then
            if (derivative) z(1) = 0
        end if
        do n = 1, (size(coeff) - 1)/(2*nfp)
            a = cmplx(coeff(2*n*nfp), -coeff(2*n*nfp + 1), dp)/2
            if (present(derivative)) then
                if (derivative) a = cmplx(0.0_dp, real(n*nfp, dp), dp)*a
            end if
            z(n + 1) = a
            z(size(z) - n + 1) = conjg(a)
        end do
        call fft_c2c(z, 1, plan, work)
        values = real(z, dp)
    end subroutine

    subroutine project(values, coeff, plan, work, nfp)
        ! Test coefficients are phi means, not Fourier amplitudes: no factor 2.
        real(dp), intent(in) :: values(:)
        real(dp), intent(out) :: coeff(:)
        type(fortnum_fft_plan_t), intent(in) :: plan
        type(fortnum_fft_workspace_t), intent(inout) :: work
        complex(dp) :: z(size(values))
        integer, intent(in) :: nfp
        integer :: n
        z = cmplx(values, 0.0_dp, dp)
        call fft_c2c(z, -1, plan, work)
        z = z/size(values)
        coeff = 0
        coeff(1) = real(z(1), dp)
        do n = 1, (size(coeff) - 1)/(2*nfp)
            coeff(2*n*nfp) = real(z(n + 1), dp)
            coeff(2*n*nfp + 1) = -aimag(z(n + 1))
        end do
    end subroutine

    subroutine point_series(field, lay, t, q, b, x, values, xq, wq, &
                            v1, v2, val, grad, plan, work, nfp, status)
        ! Geometry/basis evaluation belongs to FortFEM. Here we only gather
        ! coefficients and schedule Fourier transforms of B, grad p, grad Phi.
        type(free_b_field_t), intent(in) :: field
        type(nc_layout_t), intent(in) :: lay
        integer, intent(in) :: t, q, nfp
        real(dp), intent(in) :: b(:), x(:)
        real(dp), intent(out) :: values(:, :), xq(2), wq
        real(dp), intent(out) :: v1(:, :), v2(:, :), val(:), grad(:, :)
        type(fortnum_fft_plan_t), intent(in) :: plan
        type(fortnum_fft_workspace_t), intent(inout) :: work
        integer, intent(out) :: status
        real(dp) :: bl(field%fs%sp(2)%nloc, lay%nb)
        real(dp) :: pl(field%fs%nl0, lay%nb), fl(field%fs%nl0, lay%nb)
        real(dp) :: bm(3, lay%nb), pg(2, lay%nb), fg(2, lay%nb)
        integer :: bk, k
        call meridian_qpoint(field%fs, 2, t, q, xq, wq, v2, status)
        if (status /= 0) return
        call meridian_qpoint(field%fs, 1, t, q, xq, wq, v1, status)
        if (status /= 0) return
        call meridian_qgrad(field%fs, t, q, xq, wq, val, grad, status)
        if (status /= 0) return
        do bk = 1, lay%nb
            bl(:, bk) = b((bk - 1)*lay%nb2 + field%fs%sp(2)%dof(:, t)) &
                        *field%fs%sp(2)%sgn(:, t)
            pl(:, bk) = x(lay%ip + (bk - 1)*lay%n0 + field%fs%sp(0)%dof(:, t))
            fl(:, bk) = x(lay%iphi + (bk - 1)*lay%n0 + field%fs%sp(0)%dof(:, t))
        end do
        bm = matmul(v2, bl); pg = matmul(grad, pl); fg = matmul(grad, fl)
        do k = 1, 3
            call synthesize(bm(k, :), values(k, :), plan, work, nfp)
        end do
        call synthesize(pg(1, :), values(4, :), plan, work, nfp)
        call synthesize(matmul(val, pl), values(5, :), plan, work, nfp, .true.)
        values(5, :) = values(5, :)/xq(1)
        call synthesize(pg(2, :), values(6, :), plan, work, nfp)
        call synthesize(fg(1, :), values(7, :), plan, work, nfp)
        call synthesize(matmul(val, fl), values(8, :), plan, work, nfp, .true.)
        values(8, :) = values(8, :)/xq(1)
        call synthesize(fg(2, :), values(9, :), plan, work, nfp)
        values(10, :) = x(lay%iv)/(2*acos(-1.0_dp)*xq(1))
    end subroutine

    subroutine jvp_ops(field, problem, ops, x, direction, output, status)
        type(free_b_field_t), intent(in) :: field
        type(nc_problem_t), intent(in) :: problem
        type(nc_ops_t), intent(in) :: ops
        real(dp), intent(in) :: x(:), direction(:)
        real(dp), allocatable, intent(out) :: output(:)
        integer, intent(out) :: status
        type(nc_layout_t) :: lay
        type(fortnum_fft_plan_t) :: plan
        type(fortnum_fft_workspace_t) :: work
        real(dp), allocatable :: b(:), db(:), rl(:, :)
        integer :: t, nt, st, failure
        status = 1; lay = nc_layout(field)
        if (size(x) /= lay%n .or. size(direction) /= lay%n) return
        if (.not. all(ieee_is_finite(x))) return
        if (.not. all(ieee_is_finite(direction))) return
        call fft_c2c_plan_init(plan, ops%nphi)
        call fft_workspace_init(work, plan)
        b = csc_matvec(ops%c, x(:lay%na))
        db = csc_matvec(ops%c, direction(:lay%na))
        allocate (output(lay%n)); output = 0
        output(:lay%na) = csc_matvec(ops%s, direction(:lay%na))/field%mu0 &
                          + csc_matvec(ops%mg, direction(lay%il + 1:lay%ip))
        output(lay%il + 1:lay%ip) = csc_matvec(ops%mgt, direction(:lay%na))
        nt = field%fs%mm%topology%n_triangles
        allocate (rl(local_rows(field, lay), nt)); failure = 0
        !$omp parallel do default(shared) private(t, st) firstprivate(work) &
        !$omp& schedule(static) reduction(max:failure)
        do t = 1, nt
            call element_jvp(field, problem, lay, ops%nphi, b, db, x, direction, &
                             t, rl(:, t), plan, work, ops%nfp, st)
            failure = max(failure, st)
        end do
        !$omp end parallel do
        if (failure /= 0) return
        do t = 1, nt
            call scatter_residual(field, lay, t, rl(:, t), output)
        end do
        where (nc_fixed(field, problem)) output = 0
        call nc_reduce(field, problem, output)
        if (.not. all(ieee_is_finite(output))) return
        status = 0
    end subroutine

    subroutine element_jvp(field, problem, lay, nphi, b, db, x, dx, t, rl, &
                           plan, work, nfp, status)
        type(free_b_field_t), intent(in) :: field
        type(nc_problem_t), intent(in) :: problem
        type(nc_layout_t), intent(in) :: lay
        integer, intent(in) :: nphi, t, nfp
        real(dp), intent(in) :: b(:), db(:), x(:), dx(:)
        real(dp), intent(out) :: rl(:)
        type(fortnum_fft_plan_t), intent(in) :: plan
        type(fortnum_fft_workspace_t), intent(inout) :: work
        integer, intent(out) :: status
        real(dp) :: v1(3, field%fs%sp(1)%nloc), v2(3, field%fs%sp(2)%nloc)
        real(dp) :: val(field%fs%nl0), grad(2, field%fs%nl0), xq(2), wq
        real(dp) :: u(10, nphi), du(10, nphi), cj(3, nphi), fq(3, nphi)
        real(dp) :: cm(3, lay%nb), fm(3, lay%nb), dc(lay%nb), df(lay%nb)
        real(dp) :: rr, phi, gpar, twopi
        integer :: q, j, k, bk, n1, n0, rp, rf, n
        n1 = field%fs%sp(1)%nloc; n0 = field%fs%nl0
        rp = n1*lay%nb; rf = rp + n0*lay%nb
        twopi = 2*acos(-1.0_dp); rl = 0
        do q = 1, field%fs%nq
            call point_series(field, lay, t, q, b, x, u, xq, wq, &
                              v1, v2, val, grad, plan, work, nfp, status)
            if (status /= 0) return
            call point_series(field, lay, t, q, db, dx, du, xq, wq, &
                              v1, v2, val, grad, plan, work, nfp, status)
            if (status /= 0) return
            rr = xq(1); u(2, :) = u(2, :) + field%f_vac/rr
            do j = 1, nphi
                status = 1
                if (dot_product(u(:3, j), u(:3, j)) <= tiny(1.0_dp)) return
                phi = twopi*(j - 1)/(nphi*nfp); gpar = problem%g_par
                if (associated(problem%g_par_shape)) then
                    gpar = gpar*problem%g_par_shape(xq)
                else if (associated(problem%g_par_shape_3d)) then
                    gpar = gpar*problem%g_par_shape_3d(xq, phi)
                end if
                call coupled_current_jvp(u(1, j), u(2, j), u(3, j), u(4, j), u(5, j), &
                    u(6, j), u(7, j), u(8, j), u(9, j), u(10, j), &
                    gpar, problem%g_perp, du(1, j), du(2, j), du(3, j), &
                    du(4, j), du(5, j), du(6, j), du(7, j), &
                    du(8, j), du(9, j), du(10, j), 0.0_dp, 0.0_dp, &
                    cj(1, j), cj(2, j), cj(3, j))
                call anisotropic_flux_jvp(u(1, j), u(2, j), u(3, j), u(4, j), u(5, j), &
                    u(6, j), problem%chi_par, problem%chi_perp, du(1, j), du(2, j), &
                    du(3, j), du(4, j), du(5, j), du(6, j), 0.0_dp, 0.0_dp, &
                    fq(1, j), fq(2, j), fq(3, j))
            end do
            do k = 1, 3
                call project(cj(k, :), cm(k, :), plan, work, nfp)
                call project(fq(k, :), fm(k, :), plan, work, nfp)
            end do
            dc = 0; df = 0
            do n = 1, lay%modes
                dc(2*n) = -n*cm(2, 2*n + 1); dc(2*n + 1) = n*cm(2, 2*n)
                df(2*n) = -n*fm(2, 2*n + 1); df(2*n + 1) = n*fm(2, 2*n)
            end do
            do bk = 1, lay%nb
                rl((bk - 1)*n1 + 1:bk*n1) = rl((bk - 1)*n1 + 1:bk*n1) &
                                            - wq*matmul(cm(:, bk), v1)
                rl(rp + (bk - 1)*n0 + 1:rp + bk*n0) = &
                    rl(rp + (bk - 1)*n0 + 1:rp + bk*n0) &
                    + wq*(grad(1, :)*fm(1, bk) + grad(2, :)*fm(3, bk) &
                          + val*df(bk)/rr)
                rl(rf + (bk - 1)*n0 + 1:rf + bk*n0) = &
                    rl(rf + (bk - 1)*n0 + 1:rf + bk*n0) &
                    + wq*(grad(1, :)*cm(1, bk) + grad(2, :)*cm(3, bk) &
                          + val*dc(bk)/rr)
            end do
            rl(size(rl)) = rl(size(rl)) + wq*cm(2, 1)/rr
        end do
        status = 0
    end subroutine

    subroutine mean_leaves(field, problem, ops, x, mean, status)
        type(free_b_field_t), intent(in) :: field
        type(nc_problem_t), intent(in) :: problem
        type(nc_ops_t), intent(in) :: ops
        real(dp), intent(in) :: x(:)
        type(nc_mean_t), intent(out) :: mean
        integer, intent(out) :: status
        type(nc_layout_t) :: lay
        type(fortnum_fft_plan_t) :: plan
        type(fortnum_fft_workspace_t) :: work
        real(dp), allocatable :: b(:)
        integer :: t, st, failure, nt
        lay = nc_layout(field); nt = field%fs%mm%topology%n_triangles
        b = csc_matvec(ops%c, x(:lay%na))
        allocate (mean%dj(3, 10, field%fs%nq, nt), &
                  mean%dq(3, 6, field%fs%nq, nt))
        call fft_c2c_plan_init(plan, ops%nphi)
        call fft_workspace_init(work, plan)
        failure = 0
        !$omp parallel do default(shared) private(t, st) firstprivate(work) &
        !$omp& schedule(static) reduction(max:failure)
        do t = 1, nt
            call element_mean(field, problem, lay, ops%nphi, b, x, t, &
                      mean%dj(:, :, :, t), mean%dq(:, :, :, t), plan, work, ops%nfp, st)
            failure = max(failure, st)
        end do
        !$omp end parallel do
        status = failure
    end subroutine

    subroutine element_mean(field, problem, lay, nphi, b, x, t, dj, dq, &
                            plan, work, nfp, status)
        type(free_b_field_t), intent(in) :: field
        type(nc_problem_t), intent(in) :: problem
        type(nc_layout_t), intent(in) :: lay
        integer, intent(in) :: nphi, t, nfp
        real(dp), intent(in) :: b(:), x(:)
        real(dp), intent(out) :: dj(:, :, :), dq(:, :, :)
        type(fortnum_fft_plan_t), intent(in) :: plan
        type(fortnum_fft_workspace_t), intent(inout) :: work
        integer, intent(out) :: status
        real(dp) :: v1(3, field%fs%sp(1)%nloc), v2(3, field%fs%sp(2)%nloc)
        real(dp) :: val(field%fs%nl0), grad(2, field%fs%nl0), xq(2), wq
        real(dp) :: u(10, nphi), e(12), cj(3), fq(3), gpar, phi
        integer :: q, j, k
        dj = 0; dq = 0
        do q = 1, field%fs%nq
            call point_series(field, lay, t, q, b, x, u, xq, wq, &
                              v1, v2, val, grad, plan, work, nfp, status)
            if (status /= 0) return
            u(2, :) = u(2, :) + field%f_vac/xq(1)
            do j = 1, nphi
                status = 1
                if (dot_product(u(:3, j), u(:3, j)) <= tiny(1.0_dp)) return
                phi = 2*acos(-1.0_dp)*(j - 1)/(nphi*nfp); gpar = problem%g_par
                if (associated(problem%g_par_shape)) then
                    gpar = gpar*problem%g_par_shape(xq)
                else if (associated(problem%g_par_shape_3d)) then
                    gpar = gpar*problem%g_par_shape_3d(xq, phi)
                end if
                do k = 1, 10
                    e = 0; e(k) = 1
                    call coupled_current_jvp(u(1, j), u(2, j), u(3, j), &
                        u(4, j), u(5, j), u(6, j), u(7, j), u(8, j), u(9, j), &
                        u(10, j), gpar, problem%g_perp, e(1), e(2), e(3), &
                        e(4), e(5), e(6), e(7), e(8), e(9), e(10), &
                        e(11), e(12), cj(1), cj(2), cj(3))
                    dj(:, k, q) = dj(:, k, q) + cj/nphi
                end do
                do k = 1, 6
                    e = 0; e(k) = 1
                    call anisotropic_flux_jvp(u(1, j), u(2, j), u(3, j), u(4, j), &
                        u(5, j), u(6, j), problem%chi_par, problem%chi_perp, &
                        e(1), e(2), e(3), e(4), e(5), e(6), 0.0_dp, 0.0_dp, &
                        fq(1), fq(2), fq(3))
                    dq(:, k, q) = dq(:, k, q) + fq/nphi
                end do
            end do
        end do
        status = 0
    end subroutine

    function mode_layout(lay, n) result(one)
        type(nc_layout_t), intent(in) :: lay
        integer, intent(in) :: n
        type(nc_layout_t) :: one
        one = lay
        one%nb = merge(1, 2, n == 0)
        one%na = one%na1*one%nb; one%nv = one%n0*one%nb
        one%il = one%na; one%ip = one%na + one%nv
        one%iphi = one%na + 2*one%nv; one%iv = one%na + 3*one%nv + 1
        one%n = one%iv
    end function

    function mode_indices(lay, n) result(idx)
        type(nc_layout_t), intent(in) :: lay
        integer, intent(in) :: n
        integer, allocatable :: idx(:)
        type(nc_layout_t) :: one
        integer :: j, b0
        one = mode_layout(lay, n); b0 = max(0, 2*n - 1)
        idx = [(b0*lay%na1 + j, j=1, one%na), &
               (lay%il + b0*lay%n0 + j, j=1, one%nv), &
               (lay%ip + b0*lay%n0 + j, j=1, one%nv), &
               (lay%iphi + b0*lay%n0 + j, j=1, one%nv), lay%iv]
    end function

    subroutine mode_slice(a, row0, nrow, col0, ncol, result, status)
        ! A contiguous Fourier sub-block, used only for the constant FEEC forms.
        type(csc_t), intent(in) :: a
        integer, intent(in) :: row0, nrow, col0, ncol
        type(csc_t), intent(out) :: result
        integer, intent(out) :: status
        type(fortsparse_status_t) :: st
        integer, allocatable :: rows(:), cols(:)
        real(dp), allocatable :: vals(:)
        integer :: j, p, k, count
        count = 0
        do j = col0 + 1, col0 + ncol
            do p = a%col_ptr(j), a%col_ptr(j + 1) - 1
                k = a%row_idx(p) - row0
                if (k >= 1 .and. k <= nrow) count = count + 1
            end do
        end do
        allocate (rows(count), cols(count), vals(count)); count = 0
        do j = col0 + 1, col0 + ncol
            do p = a%col_ptr(j), a%col_ptr(j + 1) - 1
                k = a%row_idx(p) - row0
                if (k < 1 .or. k > nrow) cycle
                count = count + 1
                rows(count) = k; cols(count) = j - col0; vals(count) = a%val(p)
            end do
        end do
        call csc_from_triplet(nrow, ncol, rows, cols, vals, result, st)
        status = 1
        if (status_ok(st)) status = 0
    end subroutine

    subroutine averaged_mode(field, ops, mean, lay, n, jac, status)
        type(free_b_field_t), intent(in) :: field
        type(nc_ops_t), intent(in) :: ops
        type(nc_mean_t), intent(in) :: mean
        type(nc_layout_t), intent(in) :: lay
        integer, intent(in) :: n
        type(csc_t), intent(out) :: jac
        integer, intent(out) :: status
        type(nc_layout_t) :: one
        type(csc_t) :: c, xb, xbc, a
        type(fortsparse_status_t) :: st
        integer, allocatable :: rows(:), cols(:), brows(:), bcols(:)
        real(dp), allocatable :: vals(:), bvals(:), rl(:), jb(:, :), jd(:, :)
        integer :: t, b0, nr, count, bcount
        one = mode_layout(lay, n); b0 = max(0, 2*n - 1)
        nr = local_rows(field, one)
        allocate (rl(nr), jb(nr, field%fs%sp(2)%nloc*one%nb), &
                  jd(nr, 2*field%fs%nl0*one%nb + 1))
        allocate (rows(16), cols(16), vals(16), brows(16), bcols(16), bvals(16))
        count = 0; bcount = 0
        do t = 1, field%fs%mm%topology%n_triangles
            call averaged_element(field, mean, one, n, t, jb, jd, status)
            if (status /= 0) return
            call scatter_jacobian(field, one, t, jb, jd, rows, cols, vals, &
                                  count, brows, bcols, bvals, bcount)
        end do
        call csc_from_triplet(one%n, one%nb2*one%nb, brows(:bcount), &
                              bcols(:bcount), bvals(:bcount), xb, st)
        status = 1
        if (.not. status_ok(st)) return
        deallocate (brows, bcols, bvals)
        call mode_slice(ops%c, b0*lay%nb2, one%nb2*one%nb, &
                        b0*lay%na1, one%na, c, status)
        if (status /= 0) return
        call csc_matmul(xb, c, xbc, st)
        status = 1
        if (.not. status_ok(st)) return
        call append(xbc, 0, 0, 1.0_dp, rows, cols, vals, count)
        call mode_slice(ops%s, b0*lay%na1, one%na, b0*lay%na1, one%na, a, status)
        if (status /= 0) return
        call append(a, 0, 0, 1/field%mu0, rows, cols, vals, count)
        call mode_slice(ops%mg, b0*lay%na1, one%na, b0*lay%n0, one%nv, a, status)
        if (status /= 0) return
        call append(a, 0, one%il, 1.0_dp, rows, cols, vals, count)
        call mode_slice(ops%mgt, b0*lay%n0, one%nv, b0*lay%na1, one%na, a, status)
        if (status /= 0) return
        call append(a, one%il, 0, 1.0_dp, rows, cols, vals, count)
        if (n > 0) then
            call reserve(rows, cols, vals, count + 1)
            count = count + 1
            rows(count) = one%iv; cols(count) = one%iv; vals(count) = 1
        end if
        call csc_from_triplet(one%n, one%n, rows(:count), cols(:count), &
                              vals(:count), jac, st)
        status = 1
        if (.not. status_ok(st)) return
        if (.not. all(ieee_is_finite(jac%val))) return
        status = 0
    end subroutine

    subroutine averaged_element(field, mean, one, n, t, jb, jd, status)
        ! Three angular points integrate the quadratic sine/cosine products
        ! exactly after averaging the leaf. Only this n's two real blocks
        ! are ever formed. This is a mean-coefficient preconditioner, not
        ! the diagonal of the full phi-dependent Jacobian.
        type(free_b_field_t), intent(in) :: field
        type(nc_mean_t), intent(in) :: mean
        type(nc_layout_t), intent(in) :: one
        integer, intent(in) :: n, t
        real(dp), intent(out) :: jb(:, :), jd(:, :)
        integer, intent(out) :: status
        real(dp) :: v1(3, field%fs%sp(1)%nloc), v2(3, field%fs%sp(2)%nloc)
        real(dp) :: val(field%fs%nl0), grad(2, field%fs%nl0), xq(2), wq
        real(dp) :: bop(3, field%fs%sp(2)%nloc*one%nb)
        real(dp) :: dop(3, field%fs%nl0*one%nb), w1(3, field%fs%sp(1)%nloc*one%nb)
        real(dp) :: jc(3, size(jb, 2) + size(jd, 2)), qc(3, size(jc, 2))
        real(dp) :: tr(one%nb), dtr(one%nb), angle, w
        integer :: q, j, bk, n1, n2, n0, cb, rp, rf, np
        n1 = field%fs%sp(1)%nloc; n2 = field%fs%sp(2)%nloc
        n0 = field%fs%nl0; cb = size(jb, 2)
        rp = n1*one%nb; rf = rp + n0*one%nb
        np = merge(1, 3, n == 0); jb = 0; jd = 0
        do q = 1, field%fs%nq
            call meridian_qpoint(field%fs, 2, t, q, xq, wq, v2, status)
            if (status /= 0) return
            call meridian_qpoint(field%fs, 1, t, q, xq, wq, v1, status)
            if (status /= 0) return
            call meridian_qgrad(field%fs, t, q, xq, wq, val, grad, status)
            if (status /= 0) return
            do j = 1, np
                angle = 2*acos(-1.0_dp)*(j - 1)/np
                tr = 1; dtr = 0
                if (n > 0) then
                    tr = [cos(angle), sin(angle)]
                    dtr = n*[-sin(angle), cos(angle)]
                end if
                do bk = 1, one%nb
                    bop(:, (bk - 1)*n2 + 1:bk*n2) = tr(bk)*v2
                    w1(:, (bk - 1)*n1 + 1:bk*n1) = tr(bk)*v1
                    dop(1, (bk - 1)*n0 + 1:bk*n0) = tr(bk)*grad(1, :)
                    dop(2, (bk - 1)*n0 + 1:bk*n0) = dtr(bk)*val/xq(1)
                    dop(3, (bk - 1)*n0 + 1:bk*n0) = tr(bk)*grad(2, :)
                end do
                jc(:, :cb) = matmul(mean%dj(:, 1:3, q, t), bop)
                jc(:, cb + 1:cb + n0*one%nb) = &
                    matmul(mean%dj(:, 4:6, q, t), dop)
                jc(:, cb + n0*one%nb + 1:cb + 2*n0*one%nb) = &
                    matmul(mean%dj(:, 7:9, q, t), dop)
                jc(:, size(jc, 2)) = 0
                if (n == 0) jc(:, size(jc, 2)) = &
                    mean%dj(:, 10, q, t)/(2*acos(-1.0_dp)*xq(1))
                qc = 0
                qc(:, :cb) = matmul(mean%dq(:, 1:3, q, t), bop)
                qc(:, cb + 1:cb + n0*one%nb) = &
                    matmul(mean%dq(:, 4:6, q, t), dop)
                w = wq/np
                jb(:rp, :) = jb(:rp, :) - w*matmul(transpose(w1), jc(:, :cb))
                jd(:rp, :) = jd(:rp, :) - w*matmul(transpose(w1), jc(:, cb + 1:))
                jb(rp + 1:rf, :) = jb(rp + 1:rf, :) &
                                   + w*matmul(transpose(dop), qc(:, :cb))
                jd(rp + 1:rf, :) = jd(rp + 1:rf, :) &
                                   + w*matmul(transpose(dop), qc(:, cb + 1:))
                jb(rf + 1:size(jb, 1) - 1, :) = jb(rf + 1:size(jb, 1) - 1, :) &
                                                + w*matmul(transpose(dop), jc(:, :cb))
                jd(rf + 1:size(jd, 1) - 1, :) = jd(rf + 1:size(jd, 1) - 1, :) &
                                              + w*matmul(transpose(dop), jc(:, cb + 1:))
                if (n == 0) then
                    jb(size(jb, 1), :) = jb(size(jb, 1), :) + w*jc(2, :cb)/xq(1)
                    jd(size(jd, 1), :) = jd(size(jd, 1), :) + w*jc(2, cb + 1:)/xq(1)
                end if
            end do
        end do
        status = 0
    end subroutine

    module subroutine prepare_modes(field, problem, ops, x, factors, status)
        type(free_b_field_t), intent(in) :: field
        type(nc_problem_t), intent(in) :: problem
        type(nc_ops_t), intent(in) :: ops
        real(dp), intent(in) :: x(:)
        type(nc_mode_factor_t), intent(inout) :: factors(0:)
        integer, intent(out) :: status
        type(nc_mean_t) :: mean
        type(nc_layout_t) :: lay, one
        type(csc_t), allocatable :: reduced(:)
        logical, allocatable :: fixed(:)
        integer, allocatable :: rep(:)
        integer :: n, b, j, first, st, failure
        real(dp) :: tick
        lay = nc_layout(field); tick = -1
        call nc_phase('', tick)
        call mean_leaves(field, problem, ops, x, mean, status)
        call nc_phase('nn2-mean-leaves', tick)
        if (status /= 0) return
        allocate (reduced(0:lay%modes))
        failure = 0
        !$omp parallel do default(shared) private(n, one, fixed, rep, &
        !$omp& b, j, first, st) schedule(static) reduction(max:failure)
        do n = 0, lay%modes
            if (mod(n, problem%nfp) /= 0) cycle
            if (factors(n)%ready) call sparse_destroy(factors(n)%solver)
            factors(n)%ready = .false.
            factors(n)%indices = mode_indices(lay, n)
            one = mode_layout(lay, n)
            if (allocated(factors(n)%rhs)) &
                deallocate (factors(n)%rhs, factors(n)%solution)
            allocate (factors(n)%rhs(one%n), factors(n)%solution(one%n))
            call averaged_mode(field, ops, mean, lay, n, factors(n)%jac, st)
            if (st /= 0) then
                failure = 1
                cycle
            end if
            allocate (fixed(one%n)); fixed = .false.
            rep = [(j, j=1, one%n)]
            do b = 0, one%nb - 1
                fixed(b*one%na1 + 1:(b + 1)*one%na1) = field%fs%sp(1)%bnd
                fixed(one%il + b*one%n0 + 1:one%il + (b + 1)*one%n0) = &
                    field%fs%sp(0)%bnd
                fixed(one%ip + b*one%n0 + 1:one%ip + (b + 1)*one%n0) = &
                    field%fs%sp(0)%bnd
                fixed(one%iphi + b*one%n0 + 1:one%iphi + (b + 1)*one%n0) = &
                    field%fs%sp(0)%bnd
            end do
            if (n == 0 .and. allocated(problem%floating)) then
                first = findloc(problem%floating, .true., 1)
                if (first > 0) then
                    do j = 1, one%n0
                        if (.not. problem%floating(j)) cycle
                        fixed(one%il + j) = .false.
                        fixed(one%iphi + j) = .false.
                        rep(one%il + j) = one%il + first
                        rep(one%iphi + j) = one%iphi + first
                    end do
                end if
            end if
            call reduced_system(factors(n)%jac, fixed, rep, reduced(n), st)
            deallocate (fixed, rep)
            if (st /= 0) then
                failure = 1
                cycle
            end if
        end do
        !$omp end parallel do
        call nc_phase('nn2-mode-assembly', tick)
        status = failure
        if (failure /= 0) return
        !$omp parallel do default(shared) private(n, st) schedule(static) &
        !$omp& reduction(max:failure)
        do n = 0, lay%modes
            if (mod(n, problem%nfp) /= 0) cycle
            call factor_reduced(reduced(n), factors(n)%solver, st)
            if (st /= 0) then
                failure = 1
                cycle
            end if
            factors(n)%ready = .true.
        end do
        !$omp end parallel do
        call nc_phase('nn2-mode-factor', tick)
        status = failure
    end subroutine

    module function mode_scale(factors, x) result(scale)
        type(nc_mode_factor_t), intent(in) :: factors(0:)
        real(dp), intent(in) :: x(:)
        real(dp) :: scale(size(x))
        real(dp), allocatable :: local(:)
        integer :: n, m
        scale = 0
        do n = 0, size(factors) - 1
            if (.not. factors(n)%ready) cycle
            local = abs_matvec(factors(n)%jac, x(factors(n)%indices))
            m = size(local)
            if (n > 0) m = m - 1 ! dummy voltage belongs only to n=0
            scale(factors(n)%indices(:m)) = local(:m)
        end do
    end function

    module subroutine matrix_free_step(field, problem, ops, x, fixed, rep, rhs, &
                                      factors, dx, status, iterations, residual, &
                                      krylov_status)
        type(free_b_field_t), intent(in) :: field
        type(nc_problem_t), intent(in) :: problem
        type(nc_ops_t), intent(in) :: ops
        real(dp), intent(in) :: x(:), rhs(:)
        logical, intent(in) :: fixed(:)
        integer, intent(in) :: rep(:)
        type(nc_mode_factor_t), intent(inout) :: factors(0:)
        real(dp), allocatable, intent(out) :: dx(:)
        integer, intent(out) :: status, iterations, krylov_status
        real(dp), intent(out) :: residual
        real(dp), allocatable :: b(:), y(:), row_scale(:), test(:)
        real(dp) :: bnorm
        status = 1; iterations = 0; krylov_status = -1
        residual = huge(1.0_dp)
        if (problem%krylov_restart < 1) return
        if (problem%krylov_max < problem%krylov_restart) return
        row_scale = mode_scale(factors, spread(1.0_dp, 1, size(x)))
        call nc_reduce(field, problem, row_scale)
        where (row_scale <= tiny(1.0_dp)) row_scale = 1
        row_scale = 1/row_scale
        b = rhs*row_scale
        allocate (y(size(b))); y = 0
        bnorm = norm2(b)
        if (bnorm <= tiny(1.0_dp)) then
            dx = y; residual = 0; status = 0; krylov_status = KRYLOV_OK
            return
        end if
        call real_fgmres_operator(apply, precondition, b/bnorm, y, &
                 problem%linear_tolerance, problem%krylov_max, problem%krylov_restart, &
                                  krylov_status, iterations, residual)
        if (krylov_status /= KRYLOV_OK) return
        if (.not. all(ieee_is_finite(y))) return
        allocate (test(size(y)))
        call apply(y, test)
        residual = norm2(test - b/bnorm)
        if (residual > problem%linear_tolerance) return
        dx = bnorm*y(rep)
        where (fixed) dx = 0
        status = 0
    contains
        subroutine apply(input, output)
            real(dp), intent(in) :: input(:)
            real(dp), intent(out) :: output(:)
            real(dp), allocatable :: tangent(:), action(:)
            integer :: st, i
            tangent = input(rep)
            where (fixed) tangent = 0
            call jvp_ops(field, problem, ops, x, tangent, action, st)
            if (st /= 0) then
                output = ieee_value(1.0_dp, ieee_quiet_nan)
                return
            end if
            output = row_scale*action
            do i = 1, size(input)
                if (fixed(i) .or. rep(i) /= i) output(i) = input(i)
            end do
        end subroutine

        subroutine precondition(input, output)
            real(dp), intent(in) :: input(:)
            real(dp), intent(out) :: output(:)
            type(fortsparse_status_t) :: st
            integer :: n, m
            output = 0
            !$omp parallel do default(shared) private(n, m, st) schedule(static)
            do n = 0, size(factors) - 1
                if (.not. factors(n)%ready) cycle
                factors(n)%rhs = &
                    input(factors(n)%indices)/row_scale(factors(n)%indices)
                m = size(factors(n)%rhs)
                if (n > 0) factors(n)%rhs(m) = 0
                call sparse_solve(factors(n)%solver, factors(n)%rhs, &
                                  factors(n)%solution, st)
                if (.not. status_ok(st)) &
                    factors(n)%solution = ieee_value(1.0_dp, ieee_quiet_nan)
                if (n > 0) m = m - 1
                output(factors(n)%indices(:m)) = factors(n)%solution(:m)
            end do
            !$omp end parallel do
        end subroutine
    end subroutine

end submodule kin6d_nc_matrix_free
