module kin6d_nonnested_mapped
    ! Mapped NN1 operators (docs/nn1-mapped-3d.md). The same reference FEEC
    ! complex supports the retained phi-preserving R/Z displacement and an
    ! owned full Cartesian chart H(R,Z,phi). Incidence/curl are unchanged;
    ! existing generated Piola maps give the physical fields and K1/K2.
    ! Geometry is fixed in this operator slice, not an NN2 moving interface.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_meridian_feec_n, only: meridian_feec_t, meridian_qpoint, &
        meridian_qgrad, meridian_lagrange_nodes
    use fortfem_meridian_fourier_scalar, only: meridian_fourier_tensor_t
    use fortfem_meridian_fourier_curl, only: meridian_fourier_mass, &
        meridian_fourier_curl_curl, meridian_fourier_curl, meridian_fourier_gradient
    use fortfem_tetra_piola_maps, only: map_tetra_nedelec_covariant
    use fortnum_linalg, only: det3
    use kin6d_equilibrium_chart, only: ideal_chart_t, cylindrical_chart_t
    use fortsparse, only: csc_t, csc_matvec, csc_transpose, csc_from_triplet, &
        fortsparse_status_t, status_ok
    use fortnum_krylov, only: real_fgmres_operator, KRYLOV_OK
    use kin6d_equilibrium_nonnested, only: free_b_field_t, nn_factor_t, nn_factor, &
        nn_factor_solve_real, nn_factor_free, nn_divergence
    implicit none
    private
    public :: map_metric_t, METRIC_IDENTITY, nn_harmonic_forms
    public :: nn_map_t, nn_map_build, nn_map_jet, nn_mapped_report_t
    public :: nn_mapped_data, nn_mapped_solve, nn_mapped_error
    public :: nn_map_function_i, nn_vector_field_i
    public :: nn_map_set_chart, nn_map_point
    public :: nn_mapped_cartesian_data, nn_mapped_cartesian_error

    integer, parameter :: METRIC_MAPPED = 0, METRIC_MEAN = 1, METRIC_IDENTITY = 2
    real(dp), parameter :: pi = acos(-1.0_dp)

    type :: nn_map_t
        ! Displacement of (R, Z): Lagrange V0 times real Fourier blocks
        ! b = 0, 2n - 1 (cos n phi), 2n (sin n phi), as the N1 map unknowns.
        integer :: modes = 0
        real(dp), allocatable :: d(:, :, :) ! (n0, 2, 0:2*modes)
        ! Alternative full Cartesian H(R,Z,phi); owned, fixed geometry.
        ! The two representations are exclusive through their constructors.
        class(ideal_chart_t), allocatable :: chart
    end type

    type, extends(meridian_fourier_tensor_t) :: map_metric_t
        type(meridian_feec_t), pointer :: fs => null()
        type(nn_map_t), pointer :: map => null()
        real(dp) :: f_vac = 0
        integer :: kind = METRIC_MAPPED, nphi_mean = 0, space = 2
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

    subroutine nn_map_set_chart(chart, map)
        ! Copy the chart, not a global regularity or embedding certificate.
        ! Actual quadrature points are admitted by nn_map_point.
        class(ideal_chart_t), intent(in) :: chart
        type(nn_map_t), intent(out) :: map
        allocate (map%chart, source=chart)
    end subroutine

    pure subroutine nn_map_point(map, dofs, val, grad, x, phi, position, &
                                 p1, p2, jac, status)
        ! Physical Cartesian point and Piola maps from orthonormal reference
        ! cylindrical (R,phi,Z) components. q=(R,Z,phi) is left handed:
        ! D=[H_R,H_phi/R,H_Z], det(D)=-det(dH/dq)/R must be positive.
        ! Existing generated upstream maps own inverse/determinant algebra.
        type(nn_map_t), intent(in) :: map
        integer, intent(in) :: dofs(:)
        real(dp), intent(in) :: val(:), grad(:, :), x(2), phi
        real(dp), intent(out) :: position(3), p1(3, 3), p2(3, 3), jac
        integer, intent(out) :: status
        type(cylindrical_chart_t) :: cylinder
        real(dp) :: j(3, 3), d2(3, 3, 3), d3(3, 3, 3, 3)
        real(dp) :: d(3, 3), identity(3, 3), frame(3, 3), y(2), df(3, 3)
        integer :: i
        position = 0; p1 = 0; p2 = 0; jac = 0; status = 1
        if (.not. all(ieee_is_finite(x)) .or. .not. ieee_is_finite(phi)) return
        if (x(1) <= 0) return
        if (allocated(map%chart)) then
            if (allocated(map%d)) return
            call map%chart%jet([x, phi], position, j, d2, d3)
            d(:, 1) = j(:, 1); d(:, 2) = j(:, 3)/x(1); d(:, 3) = j(:, 2)
        else
            if (.not. allocated(map%d)) return
            if (size(val) /= size(dofs)) return
            if (size(grad, 1) /= 2 .or. size(grad, 2) /= size(dofs)) return
            if (size(dofs) < 1) return
            if (minval(dofs) < 1 .or. maxval(dofs) > size(map%d, 1)) return
            if (size(map%d, 2) /= 2 .or. size(map%d, 3) /= 2*map%modes + 1) return
            call nn_map_jet(map, dofs, val, grad, x, phi, y, df, jac)
            if (.not. all(ieee_is_finite(y)) .or. y(1) <= 0) return
            call cylinder%jet([y, phi], position, j, d2, d3)
            frame(:, 1) = j(:, 1); frame(:, 2) = j(:, 3)/y(1)
            frame(:, 3) = j(:, 2)
            d = matmul(frame, df)
        end if
        if (.not. all(ieee_is_finite(position))) return
        if (.not. all(ieee_is_finite(d))) return
        jac = det3(d)
        if (.not. ieee_is_finite(jac)) return
        identity = 0
        do i = 1, 3
            identity(i, i) = 1
        end do
        call map_tetra_nedelec_covariant(d, identity, identity, p1, p2, status)
        if (status /= 0) return
        if (.not. all(ieee_is_finite(p1)) .or. .not. all(ieee_is_finite(p2))) status = 1
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
        y = 0; df = 0; jac = 0
        ! This legacy entry has cylindrical physical outputs. Full charts use
        ! nn_map_point, whose Cartesian contract does not reinterpret them.
        if (allocated(map%chart) .or. .not. allocated(map%d)) return
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

    subroutine nn_harmonic_forms(field, metric, nphi, g, energy, volume, status)
        ! Unit harmonic h=e_phi/R: <V2,K h>, <h,K h>, and int dR dZ/R.
        ! The same forms serve the vacuum and Taylor closures.
        type(free_b_field_t), intent(in) :: field
        type(map_metric_t), intent(in) :: metric
        integer, intent(in) :: nphi
        real(dp), allocatable, intent(out) :: g(:)
        real(dp), intent(out) :: energy, volume
        integer, intent(out) :: status
        real(dp), allocatable :: local(:, :), hk(:), vol(:)
        integer :: modes, nb, nl, n2, nt, t, b, failure, st
        modes = ubound(field%b, 2); nb = 2*modes + 1
        nl = field%fs%sp(2)%nloc; n2 = field%fs%sp(2)%ndof
        nt = field%fs%mm%topology%n_triangles
        allocate (local(nl*nb, nt), hk(nt), vol(nt))
        failure = 0
        !$omp parallel do default(shared) private(t, st) schedule(static) &
        !$omp& reduction(max:failure)
        do t = 1, nt
            call element(t, local(:, t), hk(t), vol(t), st)
            failure = max(failure, st)
        end do
        !$omp end parallel do
        status = failure
        if (status /= 0) return
        allocate (g(nb*n2)); g = 0
        do t = 1, nt
            do b = 0, nb - 1
                g(b*n2 + field%fs%sp(2)%dof(:, t)) = &
                    g(b*n2 + field%fs%sp(2)%dof(:, t)) + local(b*nl + 1:(b + 1)*nl, t)
            end do
        end do
        energy = sum(hk); volume = sum(vol)
    contains
        subroutine element(t, local, hk, vol, status)
            integer, intent(in) :: t
            real(dp), intent(out) :: local(:), hk, vol
            integer, intent(out) :: status
            real(dp) :: x(2), w, basis(3, nl), pts(3, nphi), kt(3, 3, nphi)
            real(dp) :: tr(0:nb - 1), h(3), kh(3), radius
            integer :: q, j, b
            local = 0; hk = 0; vol = 0
            do q = 1, field%fs%nq
                call meridian_qpoint(field%fs, 2, t, q, x, w, basis, status)
                if (status /= 0) return
                radius = x(1)
                if (field%fs%cyl_r0 > 0) radius = field%fs%cyl_r0
                vol = vol + w/radius**2
                do j = 1, nphi
                    pts(:, j) = [x, 2*pi*(j - 1)/nphi]
                end do
                call metric%tensor(t, pts, kt, status)
                if (status /= 0) return
                h = [0.0_dp, 1/radius, 0.0_dp]
                do j = 1, nphi
                    call trig(modes, pts(3, j), tr)
                    kh = matmul(kt(:, :, j), h)
                    hk = hk + (w/nphi)*dot_product(h, kh)
                    do b = 0, nb - 1
                        local(b*nl + 1:(b + 1)*nl) = local(b*nl + 1:(b + 1)*nl) &
                            + (w/nphi)*tr(b)*matmul(kh, basis)*field%fs%sp(2)%sgn(:, t)
                    end do
                end do
            end do
            status = 0
        end subroutine
    end subroutine

    subroutine nn_mapped_data(field, map, nphi, potential, data, status)
        ! Legacy physical (R,Z,phi) samples and cylindrical vector components.
        type(free_b_field_t), intent(in), target :: field
        type(nn_map_t), intent(in), target :: map
        integer, intent(in) :: nphi
        procedure(nn_vector_field_i) :: potential
        real(dp), allocatable, intent(out) :: data(:)
        integer, intent(out) :: status
        call mapped_data(field, map, nphi, potential, data, status, .false.)
    end subroutine

    subroutine nn_mapped_cartesian_data(field, map, nphi, potential, data, status)
        ! Explicit physical Cartesian XYZ samples and Cartesian A components.
        type(free_b_field_t), intent(in), target :: field
        type(nn_map_t), intent(in), target :: map
        integer, intent(in) :: nphi
        procedure(nn_vector_field_i) :: potential
        real(dp), allocatable, intent(out) :: data(:)
        integer, intent(out) :: status
        call mapped_data(field, map, nphi, potential, data, status, .true.)
    end subroutine

    subroutine mapped_data(field, map, nphi, potential, data, status, cartesian)
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
        logical, intent(in) :: cartesian
        type(map_metric_t) :: unit
        type(csc_t) :: mass, grad
        type(nn_factor_t) :: factor
        real(dp), allocatable :: loads(:, :), load(:), lambda(:)
        integer :: modes, nb, n1, n0, nn, nl, t, b, n, failure, st
        status = 1
        if (.not. cartesian .and. allocated(map%chart)) return
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
            real(dp) :: a(3), ahat(3), angle, position(3), p1(3, 3), p2(3, 3)
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
                    if (cartesian) then
                        call nn_map_point(map, field%fs%sp(0)%dof(:, t), val, dval, &
                                          x, angle, position, p1, p2, jac, status)
                        if (status /= 0) return
                        call potential(position, a)
                        ahat = jac*matmul(transpose(p2), a)
                    else
                        call nn_map_jet(map, field%fs%sp(0)%dof(:, t), val, dval, x, &
                                        angle, y, df, jac)
                        status = 1
                        if (.not. ieee_is_finite(jac) .or. jac <= 0) return
                        call potential([y, angle], a)
                        ahat = matmul(transpose(df), a)
                    end if
                    status = 1
                    if (.not. all(ieee_is_finite(ahat))) return
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
        real(dp), allocatable :: g(:), rhs(:), x(:), r(:), y(:)
        real(dp), allocatable :: bh(:)
        real(dp) :: scale, rnorm, energy_h, volume
        integer :: modes, nb, n1, n2, nn, b, n, st
        status = 1
        modes = ubound(field%b, 2)
        if (nphi < 2*modes + 1 .or. tolerance <= 0) return
        nb = 2*modes + 1; n1 = field%fs%sp(1)%ndof; n2 = field%fs%sp(2)%ndof
        nn = field%fs%nn
        if (size(data) /= nb*n1) return
        metric%fs => field%fs; metric%map => map
        mean = metric; mean%kind = METRIC_MEAN; mean%nphi_mean = nphi
        call meridian_fourier_curl_curl(field%fs, modes, nphi, metric, s, st)
        if (st /= 0) return
        call meridian_fourier_curl(field%fs, modes, c, st)
        if (st /= 0) return
        call nn_harmonic_forms(field, metric, nphi, g, energy_h, volume, st)
        if (st /= 0) return
        g = field%f_vac*g
        energy_h = field%f_vac**2*energy_h
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
        ! Relative to the lifted right-hand side rhs - S x_D (norm scale), the
        ! quantity FGMRES reduced. With zero forcing (f_vac = 0) the free part of
        ! S x equals -r, so a denominator without scale is identically |r|.
        report%residual = norm2(r)/max(tiny(1.0_dp), scale, norm2(merge(rhs, &
            0.0_dp, .not. fixed)), norm2(merge(csc_matvec(s, x), 0.0_dp, .not. fixed)))
        report%energy = (dot_product(x, csc_matvec(s, x)) - 2*dot_product(x, rhs) &
            + energy_h)/(2*field%mu0)
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

    end subroutine

    subroutine nn_mapped_error(field, map, nphi, exact, error, deviation, &
                               energy_h, energy_exact, status)
        type(free_b_field_t), intent(in) :: field
        type(nn_map_t), intent(in) :: map
        integer, intent(in) :: nphi
        procedure(nn_vector_field_i) :: exact
        real(dp), intent(out) :: error, deviation, energy_h, energy_exact
        integer, intent(out) :: status
        call mapped_error(field, map, nphi, exact, error, deviation, &
                          energy_h, energy_exact, status, .false.)
    end subroutine

    subroutine nn_mapped_cartesian_error(field, map, nphi, exact, error, &
                                         energy_h, energy_exact, status)
        ! Same toroidal-MEAN L2 norm/energy as NN1 (no extra 2*pi factor).
        type(free_b_field_t), intent(in) :: field
        type(nn_map_t), intent(in) :: map
        integer, intent(in) :: nphi
        procedure(nn_vector_field_i) :: exact
        real(dp), intent(out) :: error, energy_h, energy_exact
        integer, intent(out) :: status
        real(dp) :: unused_deviation
        call mapped_error(field, map, nphi, exact, error, unused_deviation, &
                          energy_h, energy_exact, status, .true.)
    end subroutine

    subroutine mapped_error(field, map, nphi, exact, error, deviation, &
                            energy_h, energy_exact, status, cartesian)
        ! Physical L2 norms over the mapped domain on nphi planes (choose
        ! nphi different from the assembly rule): ||B_h - B||, ||B - F/R e_phi||
        ! and both energies int |B|^2 / (2 mu0).
        type(free_b_field_t), intent(in) :: field
        type(nn_map_t), intent(in) :: map
        integer, intent(in) :: nphi
        procedure(nn_vector_field_i) :: exact
        real(dp), intent(out) :: error, deviation, energy_h, energy_exact
        integer, intent(out) :: status
        logical, intent(in) :: cartesian
        real(dp) :: x(2), w, w0, phi(3, field%fs%sp(2)%nloc), val(field%fs%nl0)
        real(dp) :: dval(2, field%fs%nl0), y(2), df(3, 3), jac, angle, bhat(3)
        real(dp) :: bh(3), be(3), wt, position(3), p1(3, 3), p2(3, 3)
        complex(dp) :: bn(3)
        integer :: t, q, j, n, i, st, failure
        error = 0; deviation = 0; energy_h = 0; energy_exact = 0
        status = 1
        if (nphi < 1) return
        if (.not. cartesian .and. allocated(map%chart)) return
        failure = 0
        !$omp parallel do default(shared) schedule(dynamic, 4) &
        !$omp& private(t, q, j, n, i, st, x, w, w0, phi, val, dval, y, df, jac) &
        !$omp& private(angle, bhat, bh, be, wt, bn, position, p1, p2) reduction(max:failure) &
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
                    if (cartesian) then
                        call nn_map_point(map, field%fs%sp(0)%dof(:, t), val, dval, &
                                          x, angle, position, p1, p2, jac, st)
                    else
                        call nn_map_jet(map, field%fs%sp(0)%dof(:, t), val, dval, x, &
                                        angle, y, df, jac)
                        st = 0
                        if (.not. ieee_is_finite(jac) .or. jac <= 0) st = 1
                    end if
                    if (st /= 0) then
                        failure = 1
                        exit
                    end if
                    bhat = [0.0_dp, field%f_vac/x(1), 0.0_dp]
                    do n = 0, ubound(field%b, 2)
                        bn = matmul(phi, field%b(field%fs%sp(2)%dof(:, t), n))
                        bhat = bhat + merge(1, 2, n == 0) &
                            *real(bn*exp(cmplx(0.0_dp, n*angle, dp)), dp)
                    end do
                    if (cartesian) then
                        bh = matmul(p2, bhat)
                        call exact(position, be)
                    else
                        bh = matmul(df, bhat)/jac
                        call exact([y, angle], be)
                    end if
                    if (.not. all(ieee_is_finite(be)) .or. &
                        .not. all(ieee_is_finite(bh))) then
                        failure = 1
                        exit
                    end if
                    wt = w*jac/nphi
                    error = error + wt*sum((bh - be)**2)
                    if (.not. cartesian) deviation = deviation + &
                        wt*(be(1)**2 + be(3)**2 + (be(2) - field%f_vac/y(1))**2)
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
        if (self%space /= 1 .and. self%space /= 2) return
        if (self%kind == METRIC_MEAN .and. self%nphi_mean < 1) return
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
            real(dp) :: position(3), p1(3, 3), p2(3, 3), jac
            integer :: st
            call nn_map_point(self%map, self%fs%sp(0)%dof(:, t), val(:, q), &
                              dval(:, :, q), xq(:, q), angle, &
                              position, p1, p2, jac, st)
            ok = st == 0
            kk = 0
            if (.not. ok) return
            if (self%space == 1) then
                kk = jac*matmul(transpose(p1), p1)
            else
                kk = jac*matmul(transpose(p2), p2)
            end if
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
