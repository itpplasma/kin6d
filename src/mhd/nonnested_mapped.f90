module kin6d_nonnested_mapped
    ! NN1 on a phi-dependent domain (docs/nn1-mapped-3d.md). The reference
    ! meridional FEEC complex times real Fourier blocks is mapped by
    ! X(x, phi) = x + sum_b d_b(x) tr_b(phi), which keeps the toroidal angle.
    ! The 2-form transforms by the contravariant Piola map B = DF b / J, so
    ! the complex and div b = 0 are unchanged; the map enters only through
    ! the Hodge tensor K = DF^T DF / J of the upstream Fourier V2 mass.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_meridian_feec_n, only: meridian_feec_t, meridian_qpoint, &
        meridian_qgrad, meridian_lagrange_nodes
    use fortfem_meridian_fourier_scalar, only: meridian_fourier_tensor_t
    use fortfem_meridian_fourier_curl, only: meridian_fourier_mass, &
        meridian_fourier_curl_curl, meridian_fourier_curl, meridian_fourier_gradient
    use fortsparse, only: csc_t, csc_matvec, csc_transpose, csc_from_triplet, &
        fortsparse_status_t, status_ok
    use fortnum_krylov, only: real_fgmres_operator, KRYLOV_OK
    use kin6d_equilibrium_nonnested, only: free_b_field_t, nn_factor_t, nn_factor, &
        nn_factor_solve_real, nn_factor_free, nn_divergence
    implicit none
    private
    public :: nn_map_t, nn_map_build, nn_map_jet, nn_mapped_report_t
    public :: nn_mapped_data, nn_mapped_solve, nn_mapped_error
    public :: nn_map_function_i, nn_vector_field_i

    integer, parameter :: METRIC_MAPPED = 0, METRIC_MEAN = 1, METRIC_IDENTITY = 2
    real(dp), parameter :: pi = acos(-1.0_dp)

    type :: nn_map_t
        ! Displacement of (R, Z): Lagrange V0 times real Fourier blocks
        ! b = 0, 2n - 1 (cos n phi), 2n (sin n phi), as the N1 map unknowns.
        integer :: modes = 0
        real(dp), allocatable :: d(:, :, :) ! (n0, 2, 0:2*modes)
    end type

    type, extends(meridian_fourier_tensor_t) :: map_metric_t
        type(meridian_feec_t), pointer :: fs => null()
        type(nn_map_t), pointer :: map => null()
        real(dp) :: f_vac = 0
        integer :: kind = METRIC_MAPPED, nphi_mean = 0
    contains
        procedure :: tensor => metric_tensor
    end type

    type :: nn_mapped_report_t
        real(dp) :: energy = 0, residual = 0, divergence = 0
        integer :: iterations = 0, info = 0, ndof = 0
    end type

    abstract interface
        subroutine nn_map_function_i(x, phi, y)
            ! Physical (R, Z) of the reference point x = (R, Z) at angle phi.
            import :: dp
            real(dp), intent(in) :: x(2), phi
            real(dp), intent(out) :: y(2)
        end subroutine
        subroutine nn_vector_field_i(x, v)
            ! Physical components (R, phi, Z) at the point x = (R, Z, phi).
            import :: dp
            real(dp), intent(in) :: x(3)
            real(dp), intent(out) :: v(3)
        end subroutine
    end interface
contains
    subroutine nn_map_build(fs, modes, nphi, fun, map, status)
        ! Nodal interpolation in the plane, trigonometric interpolation on
        ! nphi uniform planes (exact for degree < nphi/2).
        type(meridian_feec_t), intent(in) :: fs
        integer, intent(in) :: modes, nphi
        procedure(nn_map_function_i) :: fun
        type(nn_map_t), intent(out) :: map
        integer, intent(out) :: status
        real(dp), allocatable :: nodes(:, :)
        real(dp) :: y(2), phi, tr(0:2*modes)
        integer :: i, j
        status = 1
        if (modes < 0 .or. nphi < 2*modes + 1) return
        call meridian_lagrange_nodes(fs, nodes)
        map%modes = modes
        allocate (map%d(fs%n0, 2, 0:2*modes))
        map%d = 0
        !$omp parallel do default(shared) private(i, j, y, phi, tr) schedule(static)
        do i = 1, fs%n0
            do j = 1, nphi
                phi = 2*pi*(j - 1)/nphi
                call fun(nodes(:, i), phi, y)
                call trig(modes, phi, tr)
                tr(1:) = 2*tr(1:)
                map%d(i, 1, :) = map%d(i, 1, :) + (y(1) - nodes(1, i))*tr/nphi
                map%d(i, 2, :) = map%d(i, 2, :) + (y(2) - nodes(2, i))*tr/nphi
            end do
        end do
        !$omp end parallel do
        if (all(ieee_is_finite(map%d))) status = 0
    end subroutine

    pure subroutine nn_map_jet(map, dofs, val, grad, x, phi, y, df, jac)
        ! Point y = X(x, phi), the Jacobian df in orthonormal cylindrical
        ! frames (rows physical (R, phi, Z), columns reference (R, phi, Z))
        ! and jac = det df, from the V0 values/gradients of one element.
        type(nn_map_t), intent(in) :: map
        integer, intent(in) :: dofs(:)
        real(dp), intent(in) :: val(:), grad(:, :), x(2), phi
        real(dp), intent(out) :: y(2), df(3, 3), jac
        real(dp) :: tr(0:2*map%modes), dtr(0:2*map%modes), c(2), g(2, 2), u(2)
        real(dp) :: up(2)
        integer :: i, n
        call trig(map%modes, phi, tr)
        dtr = 0
        do n = 1, map%modes
            dtr(2*n - 1) = -n*tr(2*n)
            dtr(2*n) = n*tr(2*n - 1)
        end do
        u = 0; up = 0
        g = reshape([1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2, 2])
        do i = 1, size(dofs)
            c = matmul(map%d(dofs(i), :, :), tr)
            u = u + val(i)*c
            g(:, 1) = g(:, 1) + grad(1, i)*c
            g(:, 2) = g(:, 2) + grad(2, i)*c
            up = up + val(i)*matmul(map%d(dofs(i), :, :), dtr)
        end do
        y = x + u
        df = 0
        df(1, :) = [g(1, 1), up(1)/x(1), g(1, 2)]
        df(2, 2) = y(1)/x(1)
        df(3, :) = [g(2, 1), up(2)/x(1), g(2, 2)]
        jac = df(2, 2)*(g(1, 1)*g(2, 2) - g(1, 2)*g(2, 1))
    end subroutine

    subroutine nn_mapped_data(field, map, nphi, potential, data, status)
        ! Real-block Fourier-V1 boundary data from the physical potential A:
        ! L2 projection of the pulled-back 1-form DF^T A + F ln R e_Z (the
        ! reference harmonic F/R e_phi is carried by field%f_vac), then the
        ! exact discrete gauge R A_phi = 0 for n > 0 (A - grad_F lambda).
        type(free_b_field_t), intent(in), target :: field
        type(nn_map_t), intent(in), target :: map
        integer, intent(in) :: nphi
        procedure(nn_vector_field_i) :: potential
        real(dp), allocatable, intent(out) :: data(:)
        integer, intent(out) :: status
        type(map_metric_t) :: unit
        type(csc_t) :: mass, grad
        type(nn_factor_t) :: factor
        real(dp), allocatable :: loads(:, :), load(:), lambda(:)
        integer :: modes, nb, n1, n0, nn, nl, t, b, n, failure, st
        status = 1
        modes = ubound(field%b, 2)
        if (nphi < 2*modes + 1) return
        nb = 2*modes + 1; n1 = field%fs%sp(1)%ndof; n0 = field%fs%n0
        nn = field%fs%nn; nl = field%fs%sp(1)%nloc
        allocate (loads(nl*nb, field%fs%mm%topology%n_triangles))
        failure = 0
        !$omp parallel do default(shared) private(t, st) schedule(dynamic, 4) &
        !$omp& reduction(max:failure)
        do t = 1, size(loads, 2)
            call element_load(t, loads(:, t), st)
            failure = max(failure, st)
        end do
        !$omp end parallel do
        if (failure /= 0) return
        allocate (load(nb*n1))
        load = 0
        do t = 1, size(loads, 2)
            do b = 0, nb - 1
                load(b*n1 + field%fs%sp(1)%dof(:, t)) = &
                    load(b*n1 + field%fs%sp(1)%dof(:, t)) + loads(b*nl + 1:(b + 1)*nl, t)
            end do
        end do
        unit%fs => field%fs; unit%kind = METRIC_IDENTITY
        ! The trigonometric basis is orthogonal on the nphi-point grid (nphi >=
        ! 2 modes + 1) and the identity metric is phi-independent, so the
        ! Fourier-V1 mass is blockdiag(c_b M0) with c_0 = 1, c_b > 0 = 1/2.
        ! Factor the single meridian mass M0 once: assembling all nb blocks
        ! stores O(nb^2) roundoff cross-mode entries that make the direct
        ! factorisation fail at fine meshes.
        call meridian_fourier_mass(field%fs, 1, 0, nphi, unit, mass, st)
        if (st /= 0) return
        call nn_factor(mass, spread(.false., 1, n1), factor, st)
        if (st /= 0) return
        allocate (data(nb*n1))
        do b = 0, nb - 1
            call nn_factor_solve_real(factor, load(b*n1 + 1:(b + 1)*n1), &
                data(b*n1 + 1:(b + 1)*n1), st)
            if (st /= 0) exit
            data(b*n1 + 1:(b + 1)*n1) = data(b*n1 + 1:(b + 1)*n1)* &
                merge(1.0_dp, 2.0_dp, b == 0)
        end do
        call nn_factor_free(factor)
        if (st /= 0) return
        call meridian_fourier_gradient(field%fs, modes, grad, st)
        if (st /= 0) return
        allocate (lambda(nb*n0))
        lambda = 0
        do n = 1, modes
            lambda((2*n - 1)*n0 + 1:2*n*n0) = -data(2*n*n1 + nn + 1:(2*n + 1)*n1)/n
            lambda(2*n*n0 + 1:(2*n + 1)*n0) = data((2*n - 1)*n1 + nn + 1:2*n*n1)/n
        end do
        data = data - csc_matvec(grad, lambda)
        do b = 1, nb - 1
            data(b*n1 + nn + 1:(b + 1)*n1) = 0
        end do
        status = 0
    contains
        subroutine element_load(t, local, status)
            integer, intent(in) :: t
            real(dp), intent(out) :: local(:)
            integer, intent(out) :: status
            real(dp) :: x(2), w, w0, phi(3, nl), val(field%fs%nl0)
            real(dp) :: dval(2, field%fs%nl0), tr(0:nb - 1), y(2), df(3, 3), jac
            real(dp) :: a(3), ahat(3), angle
            integer :: q, j, i, bb
            local = 0
            do q = 1, field%fs%nq
                call meridian_qpoint(field%fs, 1, t, q, x, w, phi, status)
                if (status /= 0) return
                call meridian_qgrad(field%fs, t, q, x, w0, val, dval, status)
                if (status /= 0) return
                do i = 1, nl
                    phi(:, i) = phi(:, i)*field%fs%sp(1)%sgn(i, t)
                end do
                do j = 1, nphi
                    angle = 2*pi*(j - 1)/nphi
                    call nn_map_jet(map, field%fs%sp(0)%dof(:, t), val, dval, x, &
                        angle, y, df, jac)
                    status = 1
                    if (jac <= 0) return
                    call potential([y, angle], a)
                    ahat = matmul(transpose(df), a)
                    ahat(3) = ahat(3) + field%f_vac*log(x(1))
                    call trig(modes, angle, tr)
                    do bb = 0, nb - 1
                        local(bb*nl + 1:(bb + 1)*nl) = local(bb*nl + 1:(bb + 1)*nl) &
                            + (w/nphi)*tr(bb)*matmul(ahat, phi)
                    end do
                end do
            end do
            status = 0
        end subroutine
    end subroutine

    subroutine nn_mapped_solve(field, map, nphi, data, tolerance, report, status)
        ! Minimise int |h + curl A|^2 / (2 mu0) over the mapped domain with
        ! the fixed traces of data: S A = -C^T M_2(K) h, S = C^T M_2(K) C.
        ! FGMRES right-preconditioned by the mode-diagonal operator of the
        ! phi-averaged metric (the per-n axisymmetric solves). B is stored in
        ! field%b in the NN1 complex convention, so nn_point and
        ! nn_divergence apply unchanged (reference components).
        type(free_b_field_t), intent(inout), target :: field
        type(nn_map_t), intent(in), target :: map
        integer, intent(in) :: nphi
        real(dp), intent(in) :: data(:), tolerance
        type(nn_mapped_report_t), intent(out) :: report
        integer, intent(out) :: status
        type(map_metric_t) :: metric, mean
        type(csc_t) :: s, sm, c, ct, p
        type(nn_factor_t) :: factor
        type(fortsparse_status_t) :: fst
        logical, allocatable :: fixed(:)
        real(dp), allocatable :: g(:), gl(:, :), hkh(:), rhs(:), x(:), r(:), y(:)
        real(dp), allocatable :: bh(:)
        real(dp) :: scale, rnorm
        integer :: modes, nb, n1, n2, nn, nl, t, b, n, failure, st
        status = 1
        modes = ubound(field%b, 2)
        if (nphi < 2*modes + 1 .or. tolerance <= 0) return
        nb = 2*modes + 1; n1 = field%fs%sp(1)%ndof; n2 = field%fs%sp(2)%ndof
        nn = field%fs%nn; nl = field%fs%sp(2)%nloc
        if (size(data) /= nb*n1) return
        metric%fs => field%fs; metric%map => map
        mean = metric; mean%kind = METRIC_MEAN; mean%nphi_mean = nphi
        call meridian_fourier_curl_curl(field%fs, modes, nphi, metric, s, st)
        if (st /= 0) return
        call meridian_fourier_curl(field%fs, modes, c, st)
        if (st /= 0) return
        allocate (gl(nl*nb, field%fs%mm%topology%n_triangles), &
            hkh(field%fs%mm%topology%n_triangles))
        failure = 0
        !$omp parallel do default(shared) private(t, st) schedule(dynamic, 4) &
        !$omp& reduction(max:failure)
        do t = 1, size(gl, 2)
            call harmonic_load(t, gl(:, t), hkh(t), st)
            failure = max(failure, st)
        end do
        !$omp end parallel do
        if (failure /= 0) return
        allocate (g(nb*n2))
        g = 0
        do t = 1, size(gl, 2)
            do b = 0, nb - 1
                g(b*n2 + field%fs%sp(2)%dof(:, t)) = &
                    g(b*n2 + field%fs%sp(2)%dof(:, t)) + gl(b*nl + 1:(b + 1)*nl, t)
            end do
        end do
        call csc_transpose(c, ct, fst)
        if (.not. status_ok(fst)) return
        rhs = -csc_matvec(ct, g)
        allocate (fixed(nb*n1))
        do b = 0, nb - 1
            fixed(b*n1 + 1:(b + 1)*n1) = field%fs%sp(1)%bnd
            if (b > 0) fixed(b*n1 + nn + 1:(b + 1)*n1) = .true.
        end do
        report%ndof = count(.not. fixed)
        x = merge(data, 0.0_dp, fixed)
        r = rhs - csc_matvec(s, x)
        where (fixed) r = 0
        call meridian_fourier_curl_curl(field%fs, modes, nphi, mean, sm, st)
        if (st /= 0) return
        call gauge_regularised(sm, p, st)
        if (st /= 0) return
        call nn_factor(p, fixed, factor, st)
        if (st /= 0) return
        scale = norm2(r)
        allocate (y(size(r)))
        y = 0
        if (scale > 0) then
            call real_fgmres_operator(apply, precondition, r/scale, y, tolerance, &
                400, 100, report%info, report%iterations, rnorm)
        end if
        call nn_factor_free(factor)
        if (scale > 0 .and. report%info /= KRYLOV_OK) return
        x = x + scale*y
        r = rhs - csc_matvec(s, x)
        where (fixed) r = 0
        report%residual = norm2(r)/max(tiny(1.0_dp), norm2(merge(rhs, 0.0_dp, &
            .not. fixed)), norm2(merge(csc_matvec(s, x), 0.0_dp, .not. fixed)))
        report%energy = (dot_product(x, csc_matvec(s, x)) - 2*dot_product(x, rhs) &
            + sum(hkh))/(2*field%mu0)
        bh = csc_matvec(c, x)
        field%b(:, 0) = cmplx(bh(:n2), 0.0_dp, dp)
        do n = 1, modes
            field%b(:, n) = cmplx(bh((2*n - 1)*n2 + 1:2*n*n2), &
                -bh(2*n*n2 + 1:(2*n + 1)*n2), dp)/2
        end do
        call nn_divergence(field, report%divergence, st)
        if (st /= 0) return
        status = 0
    contains
        subroutine apply(input, output)
            real(dp), intent(in) :: input(:)
            real(dp), intent(out) :: output(:)
            output = csc_matvec(s, merge(0.0_dp, input, fixed))
            where (fixed) output = input
        end subroutine

        subroutine precondition(input, output)
            real(dp), intent(in) :: input(:)
            real(dp), intent(out) :: output(:)
            integer :: info
            call nn_factor_solve_real(factor, input, output, info)
            if (info /= 0) output = input
        end subroutine

        subroutine gauge_regularised(a, out, status)
            ! Mean-metric operator plus eps M_1 on the n = 0 Nedelec block,
            ! whose curl kernel is grad V0. On residuals orthogonal to that
            ! kernel the preconditioner adds no kernel component.
            type(csc_t), intent(in) :: a
            type(csc_t), intent(out) :: out
            integer, intent(out) :: status
            type(fortsparse_status_t) :: fst
            integer, allocatable :: row(:), col(:)
            real(dp), allocatable :: val(:)
            real(dp) :: trace_a, trace_m, eps
            integer :: j, q, k, cnt
            associate (m => field%fs%mass(1))
                trace_a = 0; trace_m = 0
                do j = 1, nn
                    do q = a%col_ptr(j), a%col_ptr(j + 1) - 1
                        if (a%row_idx(q) == j) trace_a = trace_a + a%val(q)
                    end do
                    do q = m%col_ptr(j), m%col_ptr(j + 1) - 1
                        if (m%row_idx(q) == j) trace_m = trace_m + m%val(q)
                    end do
                end do
                eps = 1.0e-8_dp*trace_a/max(trace_m, tiny(1.0_dp))
                allocate (row(a%nnz + m%nnz), col(a%nnz + m%nnz), val(a%nnz + m%nnz))
                cnt = 0
                do j = 1, a%ncol
                    do q = a%col_ptr(j), a%col_ptr(j + 1) - 1
                        cnt = cnt + 1
                        row(cnt) = a%row_idx(q); col(cnt) = j; val(cnt) = a%val(q)
                    end do
                end do
                do j = 1, nn
                    do q = m%col_ptr(j), m%col_ptr(j + 1) - 1
                        k = m%row_idx(q)
                        if (k > nn) cycle
                        cnt = cnt + 1
                        row(cnt) = k; col(cnt) = j; val(cnt) = eps*m%val(q)
                    end do
                end do
            end associate
            call csc_from_triplet(a%nrow, a%ncol, row(:cnt), col(:cnt), val(:cnt), &
                out, fst)
            status = merge(0, 1, status_ok(fst))
        end subroutine

        subroutine harmonic_load(t, local, hk, status)
            ! g_i = < psi_i . K h > and < h . K h > for h = F/R e_phi.
            integer, intent(in) :: t
            real(dp), intent(out) :: local(:), hk
            integer, intent(out) :: status
            real(dp) :: xq(2), w, phi(3, nl), pts(3, nphi), kt(3, 3, nphi), h(3)
            real(dp) :: tr(0:nb - 1), kh(3)
            integer :: q, j, i, bb
            local = 0; hk = 0
            do q = 1, field%fs%nq
                call meridian_qpoint(field%fs, 2, t, q, xq, w, phi, status)
                if (status /= 0) return
                do i = 1, nl
                    phi(:, i) = phi(:, i)*field%fs%sp(2)%sgn(i, t)
                end do
                do j = 1, nphi
                    pts(:, j) = [xq, 2*pi*(j - 1)/nphi]
                end do
                call metric%tensor(t, pts, kt, status)
                if (status /= 0) return
                h = [0.0_dp, field%f_vac/xq(1), 0.0_dp]
                do j = 1, nphi
                    kh = matmul(kt(:, :, j), h)
                    hk = hk + (w/nphi)*dot_product(h, kh)
                    call trig(modes, pts(3, j), tr)
                    do bb = 0, nb - 1
                        local(bb*nl + 1:(bb + 1)*nl) = local(bb*nl + 1:(bb + 1)*nl) &
                            + (w/nphi)*tr(bb)*matmul(kh, phi)
                    end do
                end do
            end do
            status = 0
        end subroutine
    end subroutine

    subroutine nn_mapped_error(field, map, nphi, exact, error, deviation, &
                               energy_h, energy_exact, status)
        ! Physical L2 norms over the mapped domain on nphi planes (choose
        ! nphi different from the assembly rule): ||B_h - B||, ||B - F/R e_phi||
        ! and both energies int |B|^2 / (2 mu0).
        type(free_b_field_t), intent(in) :: field
        type(nn_map_t), intent(in) :: map
        integer, intent(in) :: nphi
        procedure(nn_vector_field_i) :: exact
        real(dp), intent(out) :: error, deviation, energy_h, energy_exact
        integer, intent(out) :: status
        real(dp) :: x(2), w, w0, phi(3, field%fs%sp(2)%nloc), val(field%fs%nl0)
        real(dp) :: dval(2, field%fs%nl0), y(2), df(3, 3), jac, angle, bhat(3)
        real(dp) :: bh(3), be(3), wt
        complex(dp) :: bn(3)
        integer :: t, q, j, n, i, st, failure
        error = 0; deviation = 0; energy_h = 0; energy_exact = 0
        failure = 0
        !$omp parallel do default(shared) schedule(dynamic, 4) &
        !$omp& private(t, q, j, n, i, st, x, w, w0, phi, val, dval, y, df, jac) &
        !$omp& private(angle, bhat, bh, be, wt, bn) reduction(max:failure) &
        !$omp& reduction(+:error, deviation, energy_h, energy_exact)
        do t = 1, field%fs%mm%topology%n_triangles
            do q = 1, field%fs%nq
                call meridian_qpoint(field%fs, 2, t, q, x, w, phi, st)
                if (st == 0) call meridian_qgrad(field%fs, t, q, x, w0, val, dval, st)
                if (st /= 0) then
                    failure = 1
                    exit
                end if
                do i = 1, size(phi, 2)
                    phi(:, i) = phi(:, i)*field%fs%sp(2)%sgn(i, t)
                end do
                do j = 1, nphi
                    angle = 2*pi*(j - 0.5_dp)/nphi
                    call nn_map_jet(map, field%fs%sp(0)%dof(:, t), val, dval, x, &
                        angle, y, df, jac)
                    if (jac <= 0) failure = 1
                    bhat = [0.0_dp, field%f_vac/x(1), 0.0_dp]
                    do n = 0, ubound(field%b, 2)
                        bn = matmul(phi, field%b(field%fs%sp(2)%dof(:, t), n))
                        bhat = bhat + merge(1, 2, n == 0) &
                            *real(bn*exp(cmplx(0.0_dp, n*angle, dp)), dp)
                    end do
                    bh = matmul(df, bhat)/jac
                    call exact([y, angle], be)
                    wt = w*jac/nphi
                    error = error + wt*sum((bh - be)**2)
                    deviation = deviation + wt*(be(1)**2 + be(3)**2 &
                        + (be(2) - field%f_vac/y(1))**2)
                    energy_h = energy_h + wt*sum(bh**2)/(2*field%mu0)
                    energy_exact = energy_exact + wt*sum(be**2)/(2*field%mu0)
                end do
            end do
        end do
        !$omp end parallel do
        error = sqrt(error); deviation = sqrt(deviation)
        status = failure
    end subroutine

    subroutine metric_tensor(self, t, x, k, status)
        class(map_metric_t), intent(in) :: self
        integer, intent(in) :: t
        real(dp), intent(in) :: x(:, :)
        real(dp), intent(out) :: k(:, :, :)
        integer, intent(out) :: status
        real(dp) :: xq(2, self%fs%nq), w, val(self%fs%nl0, self%fs%nq)
        real(dp) :: dval(2, self%fs%nl0, self%fs%nq), kj(3, 3)
        real(dp) :: kmean(3, 3, self%fs%nq)
        logical :: have_mean(self%fs%nq)
        integer :: p, q, j, nm
        status = 1
        k = 0
        have_mean = .false.
        if (self%kind == METRIC_IDENTITY) then
            do p = 1, size(x, 2)
                k(1, 1, p) = 1; k(2, 2, p) = 1; k(3, 3, p) = 1
            end do
            status = 0
            return
        end if
        do q = 1, self%fs%nq
            call meridian_qgrad(self%fs, t, q, xq(:, q), w, val(:, q), dval(:, :, q), &
                status)
            if (status /= 0) return
        end do
        status = 1
        do p = 1, size(x, 2)
            q = minloc(sum((xq - spread(x(:2, p), 2, self%fs%nq))**2, 1), 1)
            if (sum((xq(:, q) - x(:2, p))**2) > 1.0e-24_dp*(1 + sum(x(:2, p)**2))) return
            if (self%kind == METRIC_MEAN) then
                ! The phi mean is the same at every angle of a meridian point.
                if (.not. have_mean(q)) then
                    nm = self%nphi_mean
                    kmean(:, :, q) = 0
                    do j = 1, nm
                        if (.not. pulled(q, 2*pi*(j - 1)/nm, kj)) return
                        kmean(:, :, q) = kmean(:, :, q) + kj/nm
                    end do
                    have_mean(q) = .true.
                end if
                k(:, :, p) = kmean(:, :, q)
            else
                if (.not. pulled(q, x(3, p), k(:, :, p))) return
            end if
        end do
        status = 0
    contains
        logical function pulled(q, angle, kk) result(ok)
            integer, intent(in) :: q
            real(dp), intent(in) :: angle
            real(dp), intent(out) :: kk(3, 3)
            real(dp) :: y(2), df(3, 3), jac
            call nn_map_jet(self%map, self%fs%sp(0)%dof(:, t), val(:, q), &
                dval(:, :, q), xq(:, q), angle, y, df, jac)
            ok = jac > 0
            kk = 0
            if (ok) kk = matmul(transpose(df), df)/jac
        end function
    end subroutine

    pure subroutine trig(modes, phi, tr)
        ! tr = 1, cos phi, sin phi, cos 2 phi, ... (upstream block order).
        integer, intent(in) :: modes
        real(dp), intent(in) :: phi
        real(dp), intent(out) :: tr(0:)
        integer :: n
        tr(0) = 1
        do n = 1, modes
            tr(2*n - 1) = cos(n*phi)
            tr(2*n) = sin(n*phi)
        end do
    end subroutine
end module kin6d_nonnested_mapped
