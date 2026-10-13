module nn_vacuum_fixture
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_coils_3d, only: coil_set_t, coil_set_from_polylines, coil_set_b_cart, &
        coil_set_modes, coil_mu0
    use kin6d_equilibrium_nonnested, only: free_b_field_t, nn_velocity
    use kin6d_equilibrium_state, only: EQ_VALUES, EQ_BR, EQ_BZ, EQ_BPHI
    use kin6d_nonnested_leaves, only: segment_potential, vacuum_gauge
    use fortnum_ode_cash_karp, only: cash_karp_step
    use fortnum_multiroot, only: multiroot_hybrids
    use fortnum_status, only: fortnum_status_t
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp), r0 = 6, ftor = 100
    real(dp), parameter :: current = 5.0e6_dp, excursion = 0.02_dp
    integer, parameter :: segments = 192, nphi = 32, nmodes = 3
    type(coil_set_t) :: coils
    type(free_b_field_t) :: field
    integer :: active_mode = 1
contains
    subroutine make_coil()
        real(dp) :: points(3, segments), phi, r
        integer :: i, st
        do i = 1, segments
            phi = 2*pi*(i - 1)/segments
            r = r0 + excursion*cos(phi)
            points(:, i) = [r*cos(phi), r*sin(phi), -excursion*sin(phi)]
        end do
        call coil_set_from_polylines(coils, [segments], points, status=st)
        if (st /= 0) error stop 'coil'
    end subroutine
    subroutine circle(u, theta, x)
        real(dp), intent(in) :: u, theta
        real(dp), intent(out) :: x(2)
        x = [r0 + u*cos(theta), u*sin(theta)]
    end subroutine
    subroutine gauge_potential(x, a)
        real(dp), intent(in) :: x(2)
        complex(dp), intent(out) :: a(3)
        real(dp) :: rz(2, 1), ar, az, br, bz
        complex(dp) :: bn(3, 1, 1)
        integer :: st, part
        rz(:, 1) = x
        call coil_set_modes(coils, [current], rz, nphi, [active_mode], bn, st)
        if (st /= 0) error stop 'coil Fourier'
        a = 0
        do part = 1, 2
            if (part == 1) then
                br = real(bn(1, 1, 1), dp); bz = real(bn(3, 1, 1), dp)
            else
                br = aimag(bn(1, 1, 1)); bz = aimag(bn(3, 1, 1))
            end if
            call vacuum_gauge(x(1), br, bz, real(active_mode, dp), ar, az)
            if (part == 1) then
                a = cmplx(0.0_dp, [ar, 0.0_dp, az], dp)
            else
                a = a - [ar, 0.0_dp, az]
            end if
        end do
    end subroutine
    real(dp) function mean_psi(x) result(psi)
        real(dp), intent(in) :: x(2)
        real(dp) :: phi, p(3), a(3), coefficient, len, r1, r2
        integer :: j, s
        psi = 0
        do j = 1, nphi
            phi = 2*pi*(j - 1)/nphi
            p = [x(1)*cos(phi), x(1)*sin(phi), x(2)]
            a = 0
            do s = 1, coils%nseg
                len = norm2(coils%l(:, s))
                r1 = norm2(p - coils%a(:, s))
                r2 = norm2(p - coils%a(:, s) - coils%l(:, s))
                call segment_potential(r1, r2, len, coefficient)
                a = a + coefficient*coils%l(:, s)/len
            end do
            psi = psi + x(1)*(-a(1)*sin(phi) + a(2)*cos(phi))/nphi
        end do
        psi = psi*coil_mu0*current/(4*pi)
    end function
    subroutine direct_field(x, phi, b)
        real(dp), intent(in) :: x(2), phi
        real(dp), intent(out) :: b(3) ! (R,phi,Z)
        real(dp) :: cart(3), p(3)
        p = [x(1)*cos(phi), x(1)*sin(phi), x(2)]
        call coil_set_b_cart(coils, [current], p, cart)
        b = [cart(1)*cos(phi) + cart(2)*sin(phi), &
            -cart(1)*sin(phi) + cart(2)*cos(phi) + ftor/x(1), cart(3)]
    end subroutine
    subroutine rhs(phi, y, dy, ctx)
        real(dp), intent(in) :: phi, y(:)
        real(dp), intent(out) :: dy(:)
        class(*), intent(in), optional :: ctx
        real(dp) :: b(3)
        logical :: direct
        integer :: st
        direct = .false.
        if (present(ctx)) then
            select type (ctx)
            type is (logical)
                direct = ctx
            end select
        end if
        if (direct) then
            call direct_field(y, phi, b)
        else
            call nn_velocity(field, [y, phi], dy, st)
            if (st /= 0) then
                print '(a,3es14.6)', 'trace outside at R,Z,phi ', y, phi
                error stop 'trace left vacuum mesh'
            end if
            return
        end if
        dy = y(1)*[b(1), b(3)]/b(2)
    end subroutine
    subroutine orbit(start, turns, steps, direct, finish, span)
        real(dp), intent(in) :: start(2), turns
        integer, intent(in) :: steps
        logical, intent(in) :: direct
        real(dp), intent(out) :: finish(2), span
        real(dp) :: k1(2), k2(2), k3(2), k4(2), k5(2), k6(2), tmp(2), next(2)
        real(dp) :: err(2), phi, h, theta, oldtheta, angle, lo, hi
        integer :: i, nfev, count
        finish = start
        oldtheta = atan2(start(2), start(1) - r0)
        angle = oldtheta; lo = angle; hi = angle
        count = nint(turns*steps); h = 2*pi*turns/count; nfev = 0
        do i = 1, count
            phi = h*(i - 1)
            call cash_karp_step(rhs, phi, finish, h, .false., k1, k2, k3, k4, &
                k5, k6, tmp, next, err, nfev, ctx=direct)
            finish = next
            theta = atan2(finish(2), finish(1) - r0)
            angle = angle + atan2(sin(theta - oldtheta), cos(theta - oldtheta))
            oldtheta = theta
            lo = min(lo, angle + i*h); hi = max(hi, angle + i*h)
            ! The range cannot shrink: a completed relative rotation decides
            ! the finite-time outside predicate immediately.
            if (turns > 1 .and. hi - lo > 2*pi) exit
        end do
        span = hi - lo
    end subroutine
    real(dp) function periodic_radius(left, right, direct, steps) result(root)
        real(dp), intent(in) :: left, right
        logical, intent(in) :: direct
        integer, intent(in) :: steps
        real(dp) :: a, b, fa, fm, finish(2), span
        integer :: i
        a = left; b = right
        call orbit([a, 0.0_dp], 0.5_dp, steps, direct, finish, span)
        fa = finish(2)
        call orbit([b, 0.0_dp], 0.5_dp, steps, direct, finish, span)
        if (fa*finish(2) >= 0) error stop 'periodic orbit not bracketed'
        do i = 1, 22
            root = (a + b)/2
            call orbit([root, 0.0_dp], 0.5_dp, steps, direct, finish, span)
            fm = finish(2)
            if (fa*fm <= 0) then
                b = root
            else
                a = root; fa = fm
            end if
        end do
        root = (a + b)/2
    end function
    subroutine full_periodic(seed, root, direct)
        real(dp), intent(in) :: seed
        real(dp), intent(out) :: root(2)
        logical, intent(in) :: direct
        type(fortnum_status_t) :: st
        real(dp) :: residual(2)
        call multiroot_hybrids(map_residual, 2, [seed, 0.0_dp], root, st, &
            xtol=1.0e-8_dp, ftol=2.0e-6_dp, max_iter=30, ctx=direct)
        call map_residual(root, residual, direct)
        if (maxval(abs(residual)) > 2.0e-5_dp) error stop 'periodic map root'
    end subroutine
    subroutine map_residual(x, residual, ctx)
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: residual(:)
        class(*), intent(in), optional :: ctx
        real(dp) :: finish(2), span
        logical :: direct
        direct = .false.
        if (present(ctx)) then
            select type (ctx)
            type is (logical)
                direct = ctx
            end select
        end if
        call orbit(x, 1.0_dp, 512, direct, finish, span)
        residual = finish - x
    end subroutine
    subroutine stability(root, direct, trace, defect)
        real(dp), intent(in) :: root(2)
        logical, intent(in) :: direct
        real(dp), intent(out) :: trace, defect
        real(dp) :: start(2), plus(2), minus(2), end(2), span, jac(2, 2)
        real(dp), parameter :: delta = 1.0e-4_dp
        integer :: j
        start = root
        call orbit(start, 1.0_dp, 512, direct, end, span)
        defect = norm2(end - start)
        do j = 1, 2
            start = root; start(j) = start(j) + delta
            call orbit(start, 1.0_dp, 512, direct, plus, span)
            start(j) = start(j) - 2*delta
            call orbit(start, 1.0_dp, 512, direct, minus, span)
            jac(:, j) = (plus - minus)/(2*delta)
        end do
        trace = jac(1, 1) + jac(2, 2)
    end subroutine
    real(dp) function island_edge(outside, inside, direct, turns, steps) result(edge)
        real(dp), intent(in) :: outside, inside, turns
        logical, intent(in) :: direct
        integer, intent(in) :: steps
        real(dp) :: a, b, middle, finish(2), span
        integer :: i
        a = outside; b = inside
        call orbit([r0 - a, 0.0_dp], turns, steps, direct, finish, span)
        if (span <= 2*pi) error stop 'island outside bracket'
        call orbit([r0 - b, 0.0_dp], turns, steps, direct, finish, span)
        if (span > 2*pi) error stop 'island inside bracket'
        do i = 1, 9
            middle = (a + b)/2
            call orbit([r0 - middle, 0.0_dp], turns, steps, direct, finish, span)
            if (span > 2*pi) then
                a = middle
            else
                b = middle
            end if
        end do
        edge = (a + b)/2
    end function
end module

program test_nonnested_vacuum
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use nn_vacuum_fixture, only: pi, r0, ftor, nmodes, field, active_mode, &
        coil_mu0, &
        make_coil, circle, mean_psi, gauge_potential, direct_field, orbit, &
        periodic_radius, full_periodic, stability, island_edge
    use kin6d_mhd_mesh, only: flux_mesh_t, flux_mesh_build
    use fortfem_mesh_2d, only: mesh_2d_t
    use fortfem_meridian_feec_n, only: meridian_lagrange_nodes
    use kin6d_equilibrium_nonnested, only: nn_report_t, nn_build, nn_solve, &
        nn_vacuum_mode, nn_point
    use kin6d_equilibrium_state, only: EQ_VALUES, EQ_BR, EQ_BZ, EQ_BPHI
    use kin6d_nonnested_residual, only: nn_residual
    use kin6d_equilibrium_residual, only: equilibrium_residual_t
    implicit none
    type(flux_mesh_t) :: fm
    type(mesh_2d_t) :: topo
    type(nn_report_t) :: report
    type(equilibrium_residual_t) :: stress
    real(dp), allocatable :: rings(:), nodes(:, :)
    complex(dp), allocatable :: bc(:), seed(:)
    real(dp) :: eb(3), em(3), d3(3), points(3, 64), exact(3, 64), values(3, 64)
    real(dp) :: val(EQ_VALUES), finish(2), reference(2), span, width, width_ref
    real(dp) :: oref, xref, oh, xh, trace_o, trace_x, defect, lower, upper
    real(dp) :: opoint(2), xpoint(2), oexact(2), xexact(2)
    real(dp) :: lower_ref, upper_ref, width_long, integration_error, rate
    real(dp) :: feec_step_error, width_ref_long, scan_resolution, map_error
    integer :: lev, nr, nt, fan, st, i, n
    call make_coil()
    oref = periodic_radius(5.49_dp, 5.54_dp, .true., 512)
    xref = periodic_radius(6.69_dp, 6.72_dp, .true., 512)
    call full_periodic(oref, oexact, .true.)
    call stability(oexact, .true., trace_o, defect)
    print '(a,3es14.6)', 'direct O, trace, closure ', oref, trace_o, defect
    if (abs(trace_o) >= 2 .or. defect > 1.0e-6_dp) error stop 'direct O point'
    call full_periodic(xref, xexact, .true.)
    call stability(xexact, .true., trace_x, defect)
    print '(a,3es14.6)', 'direct X, trace, closure ', xref, trace_x, defect
    if (trace_x <= 2 .or. defect > 1.0e-6_dp) error stop 'direct X point'
    !$omp parallel sections default(shared)
    !$omp section
    lower_ref = island_edge(0.36_dp, 0.48_dp, .true., 32.0_dp, 256)
    !$omp section
    upper_ref = island_edge(0.63_dp, 0.48_dp, .true., 32.0_dp, 256)
    !$omp end parallel sections
    width_ref = upper_ref - lower_ref
    !$omp parallel sections default(shared)
    !$omp section
    lower = island_edge(0.36_dp, 0.48_dp, .true., 64.0_dp, 512)
    !$omp section
    upper = island_edge(0.63_dp, 0.48_dp, .true., 64.0_dp, 512)
    !$omp end parallel sections
    width_long = upper - lower
    width_ref_long = width_long
    print '(a,4es14.6)', 'direct island edges, width, long-run width ', &
        lower_ref, upper_ref, width_ref, width_long
    if (abs(width_long - width_ref) > 0.003_dp) error stop 'unresolved island width'
    call orbit([5.52_dp, 0.0_dp], 1.0_dp, 256, .true., reference, span)
    call orbit([5.52_dp, 0.0_dp], 1.0_dp, 512, .true., finish, span)
    integration_error = norm2(reference - finish)
    print '(a,es14.6)', 'direct integration step-doubling ', integration_error
    if (integration_error > 1.0e-6_dp) error stop 'direct integration'
    do i = 1, 64
        call circle(0.4_dp + 0.07_dp*mod(i, 8), 2*pi*(i - 0.31_dp)/64, points(:2, i))
        points(3, i) = 2*pi*(i - 0.17_dp)/19
        call direct_field(points(:2, i), points(3, i), exact(:, i))
        exact(2, i) = exact(2, i) - ftor/points(1, i)
    end do
    do lev = 1, 3
        nr = 4*2**(lev - 1)
        if (allocated(rings)) deallocate (rings, bc, seed)
        allocate (rings(nr + 1))
        do i = 1, nr + 1
            rings(i) = 0.3_dp + 0.9_dp*(i - 1)/nr
        end do
        call flux_mesh_build(rings, 8*nr, 8*nr, 2, circle, fm, st)
        if (st /= 0) error stop 'vacuum mesh'
        ! Extract the annular strips, excluding the current-bearing cavity.
        fan = fm%tlast(1)
        nt = fm%mm%topology%n_triangles - fan
        topo = mesh_2d_t()
        topo%n_triangles = nt; topo%n_vertices = fm%mm%topology%n_vertices - 1
        topo%has_triangles = .true.
        topo%vertices = fm%mm%topology%vertices(:, 2:)
        topo%triangles = fm%mm%topology%triangles(:, fan + 1:) - 1
        fm%mm%topology = topo
        fm%mm%nodes = fm%mm%nodes(:, :, fan + 1:)
        fm%tlast = fm%tlast - fan; fm%vfirst = fm%vfirst - 1
        fm%ring_tri = fm%ring_tri - fan
        call nn_build(fm, 2, nmodes, [r0, 0.0_dp], 1.0_dp, ftor, coil_mu0, &
            field, st)
        if (st /= 0) error stop 'vacuum FEEC'
        call meridian_lagrange_nodes(field%fs, nodes)
        allocate (bc(field%fs%n0), seed(field%fs%n0))
        bc = 0; seed = 0
        !$omp parallel do default(shared) private(i) schedule(static)
        do i = 1, size(bc)
            if (field%fs%sp(0)%bnd(i)) bc(i) = mean_psi(nodes(:, i))
        end do
        !$omp end parallel do
        call nn_solve(field, 0, bc, seed, report, st)
        if (st /= 0) error stop 'vacuum n0 solve'
        do n = 1, nmodes
            active_mode = n
            call nn_vacuum_mode(field, n, gauge_potential, report, st)
            if (st /= 0) error stop 'vacuum mode solve'
            if (report%divergence > 1.0e-12_dp) error stop 'vacuum closed form'
        end do
        call nn_solve(field, 0, bc, seed, report, st, pprime=1.0_dp)
        if (st == 0) error stop '3D finite pressure must be rejected'
        do i = 1, 64
            call nn_point(field, points(:, i), val, st)
            if (st /= 0) error stop 'vacuum sampling'
            values(:, i) = [val(EQ_BR), val(EQ_BPHI) - ftor/points(1, i), val(EQ_BZ)]
        end do
        eb(lev) = sqrt(sum((values - exact)**2)/sum(exact**2))
        map_error = 0
        !$omp parallel do default(shared) private(i, reference, finish, span) &
        !$omp& reduction(max:map_error) schedule(static)
        do i = 1, 6
            call orbit([5.52_dp + 0.008_dp*i, 0.0_dp], 1.0_dp, 256, &
                .true., reference, span)
            call orbit([5.52_dp + 0.008_dp*i, 0.0_dp], 1.0_dp, 256, &
                .false., finish, span)
            map_error = max(map_error, norm2(finish - reference))
        end do
        !$omp end parallel do
        em(lev) = map_error
        call nn_residual(field, 1, stress, st)
        if (st /= 0) error stop 'vacuum D3'
        d3(lev) = stress%normalized
        print '(a,i0,a,3es14.6)', 'nr=', nr, ' B / map metres / D3 ', &
            eb(lev), em(lev), d3(lev)
    end do
    rate = log(eb(1)/eb(3))/log(4.0_dp)
    print '(a,3f9.3)', 'orders field/map/D3 ', rate, &
        log(em(1)/em(3))/log(4.0_dp), log(d3(1)/d3(3))/log(4.0_dp)
    if (rate < 1.5_dp .or. eb(3) > 0.01_dp) error stop 'vacuum field convergence'
    if (em(3) > 0.01_dp .or. em(3) >= em(1)) error stop 'field-line map convergence'
    if (d3(3) >= d3(1)) error stop 'vacuum D3 convergence'
    oh = periodic_radius(5.49_dp, 5.54_dp, .false., 512)
    xh = periodic_radius(6.69_dp, 6.72_dp, .false., 512)
    call full_periodic(oh, opoint, .false.)
    call stability(opoint, .false., trace_o, defect)
    print '(a,4es14.6)', 'FEEC O(R,Z), trace, closure ', opoint, trace_o, defect
    if (abs(trace_o) >= 2) error stop 'FEEC O point'
    call full_periodic(xh, xpoint, .false.)
    call stability(xpoint, .false., trace_x, defect)
    print '(a,4es14.6)', 'FEEC X(R,Z), trace, closure ', xpoint, trace_x, defect
    if (trace_x <= 2) error stop 'FEEC X point'
    if (max(norm2(opoint - oexact), norm2(xpoint - xexact)) > 0.002_dp) &
        error stop 'O/X location'
    !$omp parallel sections default(shared)
    !$omp section
    lower = island_edge(0.36_dp, 0.48_dp, .false., 32.0_dp, 256)
    !$omp section
    upper = island_edge(0.63_dp, 0.48_dp, .false., 32.0_dp, 256)
    !$omp end parallel sections
    width = upper - lower
    print '(a,3es14.6)', 'FEEC edges and width ', lower, upper, width
    scan_resolution = (abs(0.36_dp - 0.48_dp) + abs(0.63_dp - 0.48_dp))/2**9
    print '(a,2es14.6)', 'width difference and scan resolution ', &
        abs(width - width_ref), scan_resolution
    if (abs(width - width_ref) > 0.004_dp) error stop 'island width'
    !$omp parallel sections default(shared)
    !$omp section
    lower = island_edge(0.36_dp, 0.48_dp, .false., 64.0_dp, 512)
    !$omp section
    upper = island_edge(0.63_dp, 0.48_dp, .false., 64.0_dp, 512)
    !$omp end parallel sections
    width_long = upper - lower
    print '(a,2es14.6)', '64-turn FEEC width and oracle difference ', &
        width_long, abs(width_long - width_ref_long)
    if (abs(width_long - width) > 0.003_dp) error stop 'FEEC width resolution'
    if (abs(width_long - width_ref_long) > 0.004_dp) error stop 'long-run width'
    call orbit([5.52_dp, 0.0_dp], 1.0_dp, 512, .false., reference, span)
    call orbit([5.52_dp, 0.0_dp], 1.0_dp, 1024, .false., finish, span)
    feec_step_error = norm2(reference - finish)
    print '(a,es14.6)', 'FEEC integration step-doubling ', feec_step_error
    if (feec_step_error > em(3)) error stop 'integration dominates map comparison'
    print '(a)', 'PASS: NN1 vacuum island against direct Biot-Savart tracing'
end program
