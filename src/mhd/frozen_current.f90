module kin6d_frozen_current
    ! Level 2 of the device island field hierarchy: the magnetic field of a
    ! frozen current density J on the NN1 FEEC complex times a real Fourier
    ! series in phi (modes 0..N, FortFEM's real block layout as in
    ! kin6d_nonnested_coupled). Unknowns x = [A (V1 blocks) | lambda (V0)]:
    !     <curl dA, curl A>/mu0 + <dA, grad lambda> = <dA, J>,
    !     <A, grad chi> = 0,
    ! for test functions with zero trace: the Coulomb-gauge saddle point of
    ! the coupled solver with its current frozen. Testing Ampere with grad
    ! chi shows that lambda vanishes iff the discrete J is divergence free;
    ! ||grad lambda||/||J_h|| (J_h the L2 projection of J into V1, R dR dZ
    ! times the phi mean) is the removed non-solenoidal fraction.
    ! Boundary data are the exterior field of J itself, exact because J is
    ! frozen and the boundary lies outside its support (fc_exterior): R A_phi
    ! of mode 0 and, for n > 0, the A_pol trace in the gauge A_phi = 0 from
    ! the B_R, B_Z modes (vacuum_gauge), projected by FortFEM's boundary
    ! trace. The mode-0 A_pol trace carries the toroidal flux of the
    ! poloidal current. Its circulation c is one extra unknown, fixed by the
    ! natural condition that the field of J has no toroidal component on the
    ! boundary: the unreduced Ampere residual r vanishes along the
    ! unit-circulation trace d. With x_a (data, c = 0, J) and x_b (c = 1,
    ! no data, no J), c = -d.r(x_a)/d.(r(x_b) + load) and x = x_a + c x_b.
    ! Only B = curl A of the frozen current is computed; a coil field is
    ! added by the consumer. The phi rule (nphi points) must resolve the
    ! harmonics of J up to N without aliasing those above N. With
    ! field%period = P, block n is the toroidal mode n P.
    !
    ! Exterior data from field samples (fc_casing): with a field F that
    ! satisfies curl F = mu0 J and div F = 0 in a volume V containing J
    ! (an equilibrium field minus the coil field), the field of J outside
    ! V follows from F on dV, and its potential from F on dV and in V, by
    ! integration by parts. The kernels are smooth because the domain
    ! boundary lies outside V, unlike the volume Biot-Savart sum of
    ! fc_exterior over a current that jumps at dV.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_meridian_feec_n, only: meridian_qpoint, &
        meridian_boundary_trace_points, meridian_boundary_trace_values
    use fortfem_meridian_fourier_curl, only: meridian_fourier_identity_forms, &
        meridian_fourier_diagonal_mass, meridian_fourier_gradient
    use fortfem_axisymmetric_gs_kernel, only: gs_green
    use fortsparse, only: csc_t, csc_matvec, csc_from_triplet, sparse_solver_t, &
        sparse_solve, sparse_destroy, fortsparse_status_t, status_ok
    use kin6d_equilibrium_nonnested, only: free_b_field_t, nn_factor_t, nn_factor, &
        nn_factor_solve_real, nn_factor_free
    use kin6d_nonnested_coupled, only: nc_layout_t, nc_layout, nc_store, &
        nc_factor_reduced
    use milca_nonnested_leaves, only: vacuum_gauge
    use milca_island_leaves, only: biot_savart_integrand
    implicit none
    private
    public :: fc_current_i, fc_exterior_t, fc_report_t, fc_exterior, fc_solve
    public :: fc_source_t, fc_casing

    abstract interface
        subroutine fc_current_i(x, phi, j)
            ! Physical (J_R, J_phi, J_Z) at (R, Z) = x and toroidal angle phi.
            import :: dp
            real(dp), intent(in) :: x(2), phi
            real(dp), intent(out) :: j(3)
        end subroutine
    end interface

    ! Exterior field of J at the boundary trace points xy (FortFEM order):
    ! psi = R A_phi of mode 0; bc, bs = cosine and sine amplitudes of
    ! (B_R, B_Z) for modes 1..N, B = B_0 + sum_n (bc cos n phi + bs sin n phi).
    type :: fc_exterior_t
        real(dp), allocatable :: xy(:, :), psi(:), bc(:, :, :), bs(:, :, :)
    end type

    ! Samples of F (see the module header). Points are (R, Z, phi), vectors
    ! physical (R, phi, Z). Weights (dV, and n dS with n the outward
    ! normal) count the whole torus: samples of one period of a P-periodic
    ! source carry the factor P, and fc_casing rotates them P times.
    type :: fc_source_t
        integer :: period = 1
        real(dp), allocatable :: xv(:, :), wv(:), fv(:, :)
        real(dp), allocatable :: xs(:, :), ns(:, :), fs(:, :)
    end type

    type :: fc_report_t
        real(dp) :: circulation = 0 ! mode-0 A_pol circulation c (= -toroidal flux)
        real(dp) :: projection = 0 ! ||grad lambda|| / ||J_h||
        real(dp) :: current_norm = 0, residual = 0
    end type
contains
    subroutine fc_exterior(field, current, nphi, ext, status)
        ! Biot-Savart field of J from its volume quadrature on field's mesh
        ! (FEEC rule times nphi uniform angles). Rotational symmetry turns
        ! the phi sum into a convolution: B_n = 2 pi K_{-n} J_n per source
        ! point, with the kernel modes K_m from an nphi-point DFT in the
        ! relative angle. Mode 0 R A_phi uses the exact ring Green function.
        type(free_b_field_t), intent(in) :: field
        procedure(fc_current_i) :: current
        integer, intent(in) :: nphi
        type(fc_exterior_t), intent(out) :: ext
        integer, intent(out) :: status
        real(dp), parameter :: pi = acos(-1.0_dp)
        real(dp), allocatable :: src(:, :), jm(:, :, :)
        complex(dp), allocatable :: jn(:, :, :), bn(:, :), kn(:, :, :)
        real(dp) :: v1(3, field%fs%sp(1)%nloc), xq(2), wq, phi, j(3), d(3)
        real(dp) :: kern(3, 3), c, s, coef
        integer :: modes, nt, ns, t, q, k, i, n, m, st, failure, np, per
        status = 1
        modes = ubound(field%b, 2); per = field%period
        if (nphi < 2*modes*per + 1) return
        coef = field%mu0/(4*pi)
        nt = field%fs%mm%topology%n_triangles
        ns = nt*field%fs%nq
        allocate (src(3, ns), jm(3, nphi, ns))
        failure = 0
        !$omp parallel do default(shared) private(t, q, k, i, xq, wq, v1, st, j, phi) &
        !$omp& reduction(max:failure) schedule(static)
        do t = 1, nt
            do q = 1, field%fs%nq
                i = (t - 1)*field%fs%nq + q
                call meridian_qpoint(field%fs, 1, t, q, xq, wq, v1, st)
                failure = max(failure, st)
                src(:, i) = [xq, wq]
                do k = 1, nphi
                    phi = 2*pi*(k - 1)/nphi
                    call current(xq, phi, j)
                    jm(:, k, i) = j
                end do
            end do
        end do
        !$omp end parallel do
        if (failure /= 0) return
        if (.not. all(ieee_is_finite(jm))) return
        ! Source modes J_n = mean_k J(phi_k) e^{-i n phi_k}, n = 0..N.
        allocate (jn(3, 0:modes, ns))
        do i = 1, ns
            do n = 0, modes
                jn(:, n, i) = 0
                do k = 1, nphi
                    phi = 2*pi*(k - 1)/nphi
                    jn(:, n, i) = jn(:, n, i) + jm(:, k, i) &
                        *cmplx(cos(n*per*phi), -sin(n*per*phi), dp)/nphi
                end do
            end do
        end do
        call meridian_boundary_trace_points(field%fs, ext%xy)
        np = size(ext%xy, 2)
        allocate (ext%psi(np), ext%bc(2, np, modes), ext%bs(2, np, modes))
        !$omp parallel do default(shared) private(i, k, n, m, d, kern, c, s, &
        !$omp& kn, bn, phi) schedule(dynamic)
        do i = 1, np
            allocate (bn(3, modes), kn(3, 3, modes))
            bn = 0
            ext%psi(i) = 0
            do k = 1, ns
                if (maxval(abs(jn(:, :, k))) <= tiny(1.0_dp)) cycle
                ! Ring Green function per ampere: R A_phi = mu0 G I.
                ext%psi(i) = ext%psi(i) + field%mu0*gs_green(ext%xy(1, i), &
                    ext%xy(2, i), src(1, k), src(2, k))*real(jn(2, 0, k), dp) &
                    *src(3, k)/src(1, k)
                if (modes == 0) cycle
                kn = 0
                do n = 1, nphi
                    phi = 2*pi*(n - 1)/nphi
                    c = cos(phi); s = sin(phi)
                    d = [ext%xy(1, i) - src(1, k)*c, -src(1, k)*s, &
                        ext%xy(2, i) - src(2, k)]
                    ! Columns: unit J_R, J_phi, J_Z at the source angle.
                    call biot_savart_integrand(d(1), d(2), d(3), c, s, 0.0_dp, &
                        kern(1, 1), kern(2, 1), kern(3, 1))
                    call biot_savart_integrand(d(1), d(2), d(3), -s, c, 0.0_dp, &
                        kern(1, 2), kern(2, 2), kern(3, 2))
                    call biot_savart_integrand(d(1), d(2), d(3), 0.0_dp, &
                        0.0_dp, 1.0_dp, kern(1, 3), kern(2, 3), kern(3, 3))
                    ! K_{-m} = mean K(phi) e^{+i m phi}
                    do m = 1, modes
                        kn(:, :, m) = kn(:, :, m) &
                            + kern*cmplx(cos(m*per*phi), sin(m*per*phi), dp)/nphi
                    end do
                end do
                do m = 1, modes
                    bn(:, m) = bn(:, m) + 2*pi*coef*src(3, k) &
                        *matmul(kn(:, :, m), jn(:, m, k))
                end do
            end do
            do m = 1, modes
                ! Cartesian (x, y, z) at phi = 0 is (R, phi, Z).
                ext%bc(:, i, m) = 2*real(bn([1, 3], m), dp)
                ext%bs(:, i, m) = -2*aimag(bn([1, 3], m))
            end do
            deallocate (bn, kn)
        end do
        !$omp end parallel do
        if (.not. all(ieee_is_finite(ext%psi))) return
        if (.not. all(ieee_is_finite(ext%bc))) return
        if (.not. all(ieee_is_finite(ext%bs))) return
        status = 0
    end subroutine

    subroutine fc_solve(field, current, ext, x, report, status, nphi)
        ! x = [A | lambda]; B = curl A is stored into field (nc_store).
        type(free_b_field_t), intent(inout) :: field
        procedure(fc_current_i) :: current
        type(fc_exterior_t), intent(in) :: ext
        real(dp), allocatable, intent(out) :: x(:)
        type(fc_report_t), intent(out) :: report
        integer, intent(out) :: status
        integer, intent(in), optional :: nphi
        type(nc_layout_t) :: lay
        type(csc_t) :: c, s, mg, k, reduced, mass, grad
        type(sparse_solver_t) :: solver
        type(nn_factor_t) :: mfactor
        type(fortsparse_status_t) :: st
        real(dp), allocatable :: load(:), xa(:), xb(:), da(:), db(:), r(:), ra(:)
        real(dp), allocatable :: rb(:), dvec(:), full(:), jh(:), g(:), mg_l(:)
        logical, allocatable :: fixed(:)
        real(dp) :: denominator
        integer :: n, na, nv, np, modes, ip, istat
        status = 1
        if (.not. allocated(field%b)) return
        if (.not. allocated(ext%xy)) return
        lay = nc_layout(field)
        modes = lay%modes; na = lay%na; nv = lay%nv; n = na + nv
        np = 4*modes*field%period + 1
        if (present(nphi)) np = nphi
        if (np < 2*modes*field%period + 1) return
        if (size(ext%bc, 3) /= modes) return
        call meridian_fourier_identity_forms(field%fs, modes, c, s, mg, istat, &
            field%period)
        if (istat /= 0) return
        call current_load(field, current, lay, np, load, istat)
        if (istat /= 0) return
        allocate (fixed(n))
        call fixed_mask(field, lay, fixed)
        call saddle_matrix(s, mg, field%mu0, na, nv, k, istat)
        if (istat /= 0) return
        call lifts(field, ext, lay, xa, xb, istat)
        if (istat /= 0) return
        dvec = xb ! the unit-circulation trace in the mode-0 A rows
        call reduce(k, fixed, reduced, istat)
        if (istat /= 0) return
        call nc_factor_reduced(reduced, solver, istat)
        if (istat /= 0) return
        allocate (ra(n), rb(n), da(n), db(n))
        ra = 0; ra(:na) = load
        ra = ra - csc_matvec(k, xa)
        rb = -csc_matvec(k, xb)
        where (fixed) ra = 0
        where (fixed) rb = 0
        call sparse_solve(solver, ra, da, st)
        if (.not. status_ok(st)) return
        call sparse_solve(solver, rb, db, st)
        call sparse_destroy(solver)
        if (.not. status_ok(st)) return
        xa = xa + da; xb = xb + db
        ! Unreduced residuals r = K x - [load; 0].
        r = csc_matvec(k, xa); r(:na) = r(:na) - load
        rb = csc_matvec(k, xb)
        denominator = dot_product(dvec, rb)
        if (.not. abs(denominator) > 0) return
        report%circulation = -dot_product(dvec, r)/denominator
        x = xa + report%circulation*xb
        r = csc_matvec(k, x); r(:na) = r(:na) - load
        where (fixed) r = 0
        report%residual = maxval(abs(r))/max(maxval(abs(load)), tiny(1.0_dp))
        ! Removed non-solenoidal fraction.
        call meridian_fourier_diagonal_mass(field%fs, 1, modes, mass, istat)
        if (istat /= 0) return
        call meridian_fourier_gradient(field%fs, modes, grad, istat, field%period)
        if (istat /= 0) return
        allocate (jh(na))
        call nn_factor(mass, [(.false., ip=1, na)], mfactor, istat)
        if (istat /= 0) return
        call nn_factor_solve_real(mfactor, load, jh, istat)
        call nn_factor_free(mfactor)
        if (istat /= 0) return
        g = csc_matvec(grad, x(na + 1:))
        mg_l = csc_matvec(mass, g)
        report%current_norm = sqrt(max(dot_product(load, jh), 0.0_dp))
        report%projection = sqrt(max(dot_product(g, mg_l), 0.0_dp)) &
            /max(report%current_norm, tiny(1.0_dp))
        if (.not. all(ieee_is_finite(x))) return
        allocate (full(lay%n))
        full = 0; full(:na) = x(:na)
        call nc_store(field, full, istat)
        if (istat /= 0) return
        status = 0
    end subroutine

    subroutine current_load(field, current, lay, nphi, load, status)
        ! load = <w, J> on every V1 block with the phi rule of nphi points.
        type(free_b_field_t), intent(in) :: field
        procedure(fc_current_i) :: current
        type(nc_layout_t), intent(in) :: lay
        integer, intent(in) :: nphi
        real(dp), allocatable, intent(out) :: load(:)
        integer, intent(out) :: status
        real(dp), allocatable :: local(:, :)
        real(dp) :: v1(3, field%fs%sp(1)%nloc), xq(2), wq, phi, j(3), tr
        integer :: t, q, k, bk, n1, nt, st, failure, mode, i
        n1 = field%fs%sp(1)%nloc
        nt = field%fs%mm%topology%n_triangles
        allocate (local(n1*lay%nb, nt))
        failure = 0
        !$omp parallel do default(shared) private(t, q, k, bk, xq, wq, v1, st, &
        !$omp& phi, j, tr, mode) reduction(max:failure) schedule(static)
        do t = 1, nt
            local(:, t) = 0
            do q = 1, field%fs%nq
                call meridian_qpoint(field%fs, 1, t, q, xq, wq, v1, st)
                failure = max(failure, st)
                do k = 1, nphi
                    phi = 2*acos(-1.0_dp)*(k - 1)/nphi
                    call current(xq, phi, j)
                    do bk = 1, lay%nb
                        mode = bk/2*field%period
                        if (bk == 1) then
                            tr = 1
                        else if (mod(bk, 2) == 0) then
                            tr = cos(mode*phi)
                        else
                            tr = sin(mode*phi)
                        end if
                        local((bk - 1)*n1 + 1:bk*n1, t) = &
                            local((bk - 1)*n1 + 1:bk*n1, t) &
                            + wq/nphi*tr*matmul(j, v1)
                    end do
                end do
            end do
        end do
        !$omp end parallel do
        status = failure
        if (status /= 0) return
        allocate (load(lay%na))
        load = 0
        do t = 1, nt
            do bk = 1, lay%nb
                do i = 1, n1
                    associate (g => lay%ia + (bk - 1)*lay%na1 &
                            + field%fs%sp(1)%dof(i, t))
                        load(g) = load(g) + field%fs%sp(1)%sgn(i, t) &
                            *local((bk - 1)*n1 + i, t)
                    end associate
                end do
            end do
        end do
        if (.not. all(ieee_is_finite(load))) status = 1
    end subroutine

    subroutine fixed_mask(field, lay, fixed)
        ! A trace and lambda on the boundary, in every block.
        type(free_b_field_t), intent(in) :: field
        type(nc_layout_t), intent(in) :: lay
        logical, intent(out) :: fixed(:)
        integer :: b
        do b = 0, lay%nb - 1
            fixed(b*lay%na1 + 1:(b + 1)*lay%na1) = field%fs%sp(1)%bnd
            fixed(lay%na + b*lay%n0 + 1:lay%na + (b + 1)*lay%n0) = &
                field%fs%sp(0)%bnd
        end do
    end subroutine

    subroutine lifts(field, ext, lay, xa, xb, status)
        ! xa: the exterior data on the fixed A dofs; xb: a mode-0 A_pol
        ! trace of unit circulation, (A_R, A_Z) = (-(Z - Z_c), R - R_c)/(2 S)
        ! with S the cross-section area (Stokes).
        type(free_b_field_t), intent(in) :: field
        type(fc_exterior_t), intent(in) :: ext
        type(nc_layout_t), intent(in) :: lay
        real(dp), allocatable, intent(out) :: xa(:), xb(:)
        integer, intent(out) :: status
        complex(dp), allocatable :: v(:, :), values(:)
        real(dp) :: ar, az, area, xq(2), wq, v0(3, field%fs%nl0)
        integer :: np, i, m, t, q, o
        np = size(ext%xy, 2)
        allocate (xa(lay%na + lay%nv), xb(lay%na + lay%nv), v(3, np))
        xa = 0; xb = 0
        v = 0
        do i = 1, np
            v(2, i) = ext%psi(i)/ext%xy(1, i)
        end do
        call meridian_boundary_trace_values(field%fs, v, values, status)
        if (status /= 0) return
        xa(field%fs%nn + 1:lay%na1) = real(values(field%fs%nn + 1:), dp)
        do m = 1, lay%modes
            ! cos block <- gauge(B^s), sin block <- -gauge(B^c) (rotation of
            ! A = i gauge(B) for B = Re[(c - i s) e^{i n phi}]).
            do i = 1, np
                call vacuum_gauge(ext%xy(1, i), ext%bs(1, i, m), ext%bs(2, i, m), &
                    real(m*field%period, dp), ar, az)
                v(:, i) = [cmplx(ar, 0.0_dp, dp), (0.0_dp, 0.0_dp), &
                    cmplx(az, 0.0_dp, dp)]
                call vacuum_gauge(ext%xy(1, i), ext%bc(1, i, m), ext%bc(2, i, m), &
                    real(m*field%period, dp), ar, az)
                v(:, i) = v(:, i) + cmplx(0.0_dp, -1.0_dp, dp)*[ar, 0.0_dp, az]
            end do
            call meridian_boundary_trace_values(field%fs, v, values, status)
            if (status /= 0) return
            o = (2*m - 1)*lay%na1
            xa(o + 1:o + field%fs%nn) = real(values(:field%fs%nn), dp)
            o = 2*m*lay%na1
            xa(o + 1:o + field%fs%nn) = aimag(values(:field%fs%nn))
        end do
        area = 0
        do t = 1, field%fs%mm%topology%n_triangles
            do q = 1, field%fs%nq
                call meridian_qpoint(field%fs, 0, t, q, xq, wq, v0, status)
                if (status /= 0) return
                area = area + wq/xq(1)
            end do
        end do
        do i = 1, np
            v(:, i) = [cmplx(-(ext%xy(2, i) - field%centre(2)), 0.0_dp, dp), &
                (0.0_dp, 0.0_dp), cmplx(ext%xy(1, i) - field%centre(1), 0.0_dp, &
                dp)]/(2*area)
        end do
        call meridian_boundary_trace_values(field%fs, v, values, status)
        if (status /= 0) return
        xb(:field%fs%nn) = real(values(:field%fs%nn), dp)
    end subroutine

    subroutine fc_casing(field, src, nphi, ext, status)
        ! Exterior data of the frozen current from its samples at the
        ! boundary trace points, which lie outside V. With d = x - x' and
        ! a = d/|d|^3, by parts from mu0 J = curl F and div F = 0:
        !   B_J(x) = (1/4 pi) <(n x F) x a + (n . F) a>_dV,
        !   A_J(x) = (1/4 pi) [<(n x F)/|d|>_dV + (F x a)_V].
        ! B is sampled at nphi angles of one period and transformed to the
        ! modes k P, k = 1..N; psi is R times the phi mean of A_phi.
        type(free_b_field_t), intent(in) :: field
        type(fc_source_t), intent(in) :: src
        integer, intent(in) :: nphi
        type(fc_exterior_t), intent(out) :: ext
        integer, intent(out) :: status
        real(dp), parameter :: pi = acos(-1.0_dp)
        real(dp), allocatable :: ps(:, :), qs(:, :), fs(:), pv(:, :), qv(:, :)
        complex(dp), allocatable :: bn(:, :)
        real(dp) :: x(3), d(3), b(3), a(3), r2, ri, ang, c, s, aphi
        integer :: modes, per, np, ns, nv, i, l, k, m
        status = 1
        modes = ubound(field%b, 2); per = field%period
        if (src%period /= per) return
        if (nphi < 2*modes + 1) return
        call rotated(src%xs, src%ns, src%fs, ps, qs, fs)
        call rotated(src%xv, spread(src%wv, 1, 3), src%fv, pv, qv)
        call meridian_boundary_trace_points(field%fs, ext%xy)
        np = size(ext%xy, 2); ns = size(ps, 2); nv = size(pv, 2)
        allocate (ext%psi(np), ext%bc(2, np, modes), ext%bs(2, np, modes))
        !$omp parallel do default(shared) private(i, l, k, m, x, d, b, a, r2, ri, &
        !$omp& ang, c, s, aphi, bn) schedule(dynamic)
        do i = 1, np
            allocate (bn(2, modes))
            bn = 0
            ext%psi(i) = 0
            do l = 0, nphi - 1
                ang = 2*pi*l/(nphi*per)
                c = cos(ang); s = sin(ang)
                x = [ext%xy(1, i)*c, ext%xy(1, i)*s, ext%xy(2, i)]
                b = 0; a = 0
                do k = 1, ns
                    ! ps: position, qs(1:3) n dS, qs(4:6) n x F dS, fs n.F dS
                    d = x - ps(:, k)
                    r2 = dot_product(d, d); ri = 1/sqrt(r2)
                    b = b + (cross(qs(4:6, k), d) + fs(k)*d)*ri**3
                    a = a + qs(4:6, k)*ri
                end do
                do k = 1, nv
                    d = x - pv(:, k)
                    ri = 1/sqrt(dot_product(d, d))
                    a = a + cross(qv(4:6, k), d)*ri**3
                end do
                b = b/(4*pi); a = a/(4*pi)
                aphi = -a(1)*s + a(2)*c
                ext%psi(i) = ext%psi(i) + ext%xy(1, i)*aphi/nphi
                do m = 1, modes
                    bn(:, m) = bn(:, m) + [b(1)*c + b(2)*s, b(3)] &
                        *cmplx(cos(m*l*2*pi/nphi), -sin(m*l*2*pi/nphi), dp)/nphi
                end do
            end do
            do m = 1, modes
                ext%bc(:, i, m) = 2*real(bn(:, m), dp)
                ext%bs(:, i, m) = -2*aimag(bn(:, m))
            end do
            deallocate (bn)
        end do
        !$omp end parallel do
        if (.not. all(ieee_is_finite(ext%psi))) return
        if (.not. all(ieee_is_finite(ext%bc))) return
        if (.not. all(ieee_is_finite(ext%bs))) return
        status = 0
    contains
        subroutine rotated(xc, wvec, f, p, q, nf)
            ! The per rotated copies in Cartesian form, weights / per:
            ! p position, q(1:3) the weight vector, q(4:6) its cross product
            ! with F (surface) or W F (volume), nf = weight . F.
            real(dp), intent(in) :: xc(:, :), wvec(:, :), f(:, :)
            real(dp), allocatable, intent(out) :: p(:, :), q(:, :)
            real(dp), allocatable, intent(out), optional :: nf(:)
            real(dp) :: cs, sn, wc(3), fc(3)
            integer :: n, j, cp, o
            n = size(xc, 2)
            allocate (p(3, n*per), q(6, n*per))
            if (present(nf)) allocate (nf(n*per))
            do cp = 0, per - 1
                do j = 1, n
                    o = cp*n + j
                    cs = cos(xc(3, j) + 2*pi*cp/per); sn = sin(xc(3, j) + 2*pi*cp/per)
                    p(:, o) = [xc(1, j)*cs, xc(1, j)*sn, xc(2, j)]
                    fc = [f(1, j)*cs - f(2, j)*sn, f(1, j)*sn + f(2, j)*cs, f(3, j)]
                    if (present(nf)) then
                        wc = [wvec(1, j)*cs - wvec(2, j)*sn, &
                            wvec(1, j)*sn + wvec(2, j)*cs, wvec(3, j)]/per
                        q(1:3, o) = wc
                        q(4:6, o) = cross(wc, fc)
                        nf(o) = dot_product(wc, fc)
                    else
                        q(1:3, o) = 0
                        q(4:6, o) = wvec(1, j)/per*fc ! volume: the scalar dV
                    end if
                end do
            end do
        end subroutine
    end subroutine

    pure function cross(a, b) result(c)
        real(dp), intent(in) :: a(3), b(3)
        real(dp) :: c(3)
        c = [a(2)*b(3) - a(3)*b(2), a(3)*b(1) - a(1)*b(3), a(1)*b(2) - a(2)*b(1)]
    end function

    subroutine saddle_matrix(s, mg, mu0, na, nv, k, status)
        ! K = [[S/mu0, MG], [MG^T, 0]].
        type(csc_t), intent(in) :: s, mg
        real(dp), intent(in) :: mu0
        integer, intent(in) :: na, nv
        type(csc_t), intent(out) :: k
        integer, intent(out) :: status
        type(fortsparse_status_t) :: st
        integer, allocatable :: row(:), col(:)
        real(dp), allocatable :: val(:)
        integer :: j, p, count
        status = 1
        allocate (row(s%nnz + 2*mg%nnz), col(s%nnz + 2*mg%nnz), &
            val(s%nnz + 2*mg%nnz))
        count = 0
        do j = 1, na
            do p = s%col_ptr(j), s%col_ptr(j + 1) - 1
                count = count + 1
                row(count) = s%row_idx(p); col(count) = j; val(count) = s%val(p)/mu0
            end do
        end do
        do j = 1, nv
            do p = mg%col_ptr(j), mg%col_ptr(j + 1) - 1
                count = count + 1
                row(count) = mg%row_idx(p); col(count) = na + j; val(count) = mg%val(p)
                count = count + 1
                row(count) = na + j; col(count) = mg%row_idx(p); val(count) = mg%val(p)
            end do
        end do
        call csc_from_triplet(na + nv, na + nv, row(:count), col(:count), &
            val(:count), k, st)
        if (.not. status_ok(st)) return
        status = 0
    end subroutine

    subroutine reduce(a, fixed, reduced, status)
        ! Unit rows and columns for the fixed dofs.
        type(csc_t), intent(in) :: a
        logical, intent(in) :: fixed(:)
        type(csc_t), intent(out) :: reduced
        integer, intent(out) :: status
        type(fortsparse_status_t) :: st
        integer, allocatable :: row(:), col(:)
        real(dp), allocatable :: val(:)
        integer :: n, j, p, count
        n = size(fixed); status = 1
        allocate (row(a%nnz + n), col(a%nnz + n), val(a%nnz + n))
        count = 0
        do j = 1, n
            if (fixed(j)) then
                count = count + 1
                row(count) = j; col(count) = j; val(count) = 1
                cycle
            end if
            do p = a%col_ptr(j), a%col_ptr(j + 1) - 1
                if (fixed(a%row_idx(p))) cycle
                count = count + 1
                row(count) = a%row_idx(p); col(count) = j; val(count) = a%val(p)
            end do
        end do
        call csc_from_triplet(n, n, row(:count), col(:count), val(:count), &
            reduced, st)
        if (.not. status_ok(st)) return
        status = 0
    end subroutine
end module
