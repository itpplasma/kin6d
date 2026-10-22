module nn_vacuum_fixture
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_coils_3d, only: coil_set_t, coil_set_from_polylines, coil_set_b_cart, &
        coil_set_modes, coil_mu0
    use kin6d_equilibrium_nonnested, only: free_b_field_t, nn_report_t, &
                                           nn_build, nn_solve, nn_vacuum_mode
    use fortfem_meridian_feec_n, only: meridian_lagrange_nodes
    use kin6d_equilibrium_state, only: EQ_VALUES, EQ_BR, EQ_BZ, EQ_BPHI
    use kin6d_nonnested_leaves, only: segment_potential, vacuum_gauge
    use kin6d_equilibrium_state, only: equilibrium_field_t
    use kin6d_island_diagnostics, only: island_trace_options_t, island_point_t, &
                           island_trace, island_symmetry_point, island_periodic_point, &
                                        island_stability, island_libration_edge
    use kin6d_mhd_mesh, only: flux_mesh_t, flux_mesh_build
    use fortfem_mesh_2d, only: mesh_2d_t
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp), r0 = 6, ftor = 100
    real(dp), parameter :: current = 5.0e6_dp, excursion = 0.02_dp
    integer, parameter :: segments = 192, nphi = 32, nmodes = 3
    type(coil_set_t) :: coils
    type(free_b_field_t) :: field
    integer :: active_mode = 1
    type, extends(equilibrium_field_t) :: direct_coil_field_t
    contains
        procedure :: sample => direct_sample
    end type
    type(direct_coil_field_t) :: direct_coil
contains
    subroutine solve_vacuum(nr, order, modes)
        integer, intent(in) :: nr, order, modes
        type(flux_mesh_t) :: fm
        type(nn_report_t) :: report
        real(dp), allocatable :: nodes(:, :)
        complex(dp), allocatable :: bc(:), seed(:)
        integer :: j, n, st
        call annulus_mesh(nr, order, fm)
        call nn_build(fm, order, modes, [r0, 0.0_dp], 1.0_dp, ftor, &
                      coil_mu0, field, st)
        if (st /= 0) error stop 'vacuum build'
        call meridian_lagrange_nodes(field%fs, nodes)
        allocate (bc(field%fs%n0), seed(field%fs%n0))
        bc = 0; seed = 0
        !$omp parallel do default(shared) private(j) schedule(static)
        do j = 1, size(bc)
            if (field%fs%sp(0)%bnd(j)) bc(j) = mean_psi(nodes(:, j))
        end do
        !$omp end parallel do
        call nn_solve(field, 0, bc, seed, report, st)
        if (st /= 0) error stop 'vacuum n0'
        do n = 1, modes
            active_mode = n
            call nn_vacuum_mode(field, n, gauge_potential, report, st)
            if (st /= 0) then
                print '(a,4i10,a,3es16.7)', 'vacuum failure P/nr/mode/status ', &
                    order, nr, n, st, ' residual/energy-before/after ', &
                    report%residual, report%energy_before, report%energy_after
                error stop 'vacuum Fourier mode'
            end if
        end do
    end subroutine
    subroutine annulus_mesh(nr, order, fm)
        ! Rings u in [0.3, 1.2] about the coil, without the current-bearing
        ! cavity (the axis fan is removed).
        integer, intent(in) :: nr, order
        type(flux_mesh_t), intent(out) :: fm
        type(mesh_2d_t) :: topo
        real(dp) :: rings(nr + 1)
        integer :: i, st, fan, nt
        do i = 1, nr + 1
            rings(i) = 0.3_dp + 0.9_dp*(i - 1)/nr
        end do
        call flux_mesh_build(rings, 8*nr, 8*nr, order, circle, fm, st)
        if (st /= 0) error stop 'vacuum mesh'
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
    end subroutine
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
    subroutine direct_sample(field, points, values, statuses)
        class(direct_coil_field_t), intent(in) :: field
        real(dp), intent(in) :: points(:, :)
        real(dp), intent(out) :: values(:, :)
        integer, intent(out) :: statuses(:)
        real(dp) :: b(3), phi
        integer :: i
        values = 0; statuses = 0
        do i = 1, size(points, 2)
            phi = 0
            if (size(points, 1) >= 3) phi = points(3, i)
            call direct_field(points(:2, i), phi, b)
            values(EQ_BR, i) = b(1); values(EQ_BZ, i) = b(3)
            values(EQ_BPHI, i) = b(2)
        end do
    end subroutine
    subroutine orbit(start, turns, steps, direct, finish, span)
        real(dp), intent(in) :: start(2), turns
        integer, intent(in) :: steps
        logical, intent(in) :: direct
        real(dp), intent(out) :: finish(2), span
        type(island_trace_options_t) :: opt
        integer :: st
        opt%fixed_steps = steps
        if (direct) then
            call island_trace(direct_coil, start, 0.0_dp, turns, opt, finish, &
                              st, [r0, 0.0_dp], 1, -1, span)
        else
            call island_trace(field, start, 0.0_dp, turns, opt, finish, &
                              st, [r0, 0.0_dp], 1, -1, span)
        end if
        if (st /= 0) error stop 'island trace domain or integration failure'
    end subroutine
    real(dp) function periodic_radius(left, right, direct, steps) result(root)
        real(dp), intent(in) :: left, right
        logical, intent(in) :: direct
        integer, intent(in) :: steps
        type(island_trace_options_t) :: opt
        real(dp) :: point(2)
        integer :: st
        opt%fixed_steps = steps; opt%root_tol = 1.0e-8_dp
        if (direct) then
            call island_symmetry_point(direct_coil, [left, right], 0.0_dp, &
                                       1.0_dp, opt, point, st)
        else
            call island_symmetry_point(field, [left, right], 0.0_dp, &
                                       1.0_dp, opt, point, st)
        end if
        if (st /= 0) error stop 'periodic orbit not bracketed'
        root = point(1)
    end function
    subroutine full_periodic(seed, root, direct)
        real(dp), intent(in) :: seed
        real(dp), intent(out) :: root(2)
        logical, intent(in) :: direct
        type(island_trace_options_t) :: opt
        type(island_point_t) :: report
        opt%fixed_steps = 512; opt%root_tol = 2.0e-5_dp
        if (direct) then
            call island_periodic_point(direct_coil, [seed, 0.0_dp], 0.0_dp, &
                                       1, [r0, 0.0_dp], 1, -1, opt, report)
        else
            call island_periodic_point(field, [seed, 0.0_dp], 0.0_dp, &
                                       1, [r0, 0.0_dp], 1, -1, opt, report)
        end if
        print '(a,i0,a,es12.4)', 'periodic status ', report%status, &
            ' closure ', report%closure
        if (report%status /= 0) error stop 'periodic map root'
        root = report%x
    end subroutine
    subroutine stability(root, direct, trace, defect)
        real(dp), intent(in) :: root(2)
        logical, intent(in) :: direct
        real(dp), intent(out) :: trace, defect
        type(island_trace_options_t) :: opt
        real(dp) :: jac(2, 2)
        integer :: st
        opt%fixed_steps = 512; opt%map_delta = 1.0e-4_dp
        if (direct) then
            call island_stability(direct_coil, root, 0.0_dp, 1.0_dp, opt, &
                                  jac, defect, st)
        else
            call island_stability(field, root, 0.0_dp, 1.0_dp, opt, &
                                  jac, defect, st)
        end if
        if (st /= 0) error stop 'stability trace'
        trace = jac(1, 1) + jac(2, 2)
    end subroutine
    real(dp) function island_edge(outside, inside, direct, turns, steps) result(edge)
        real(dp), intent(in) :: outside, inside, turns
        logical, intent(in) :: direct
        integer, intent(in) :: steps
        type(island_trace_options_t) :: opt
        integer :: st
        opt%fixed_steps = steps
        if (direct) then
            call island_libration_edge(direct_coil, outside, inside, &
                                [r0, 0.0_dp], [-1.0_dp, 0.0_dp], 0.0_dp, turns, 1, -1, &
                                       opt, edge, st, scans=9)
        else
            call island_libration_edge(field, outside, inside, &
                                [r0, 0.0_dp], [-1.0_dp, 0.0_dp], 0.0_dp, turns, 1, -1, &
                                       opt, edge, st, scans=9)
        end if
        if (st /= 0) error stop 'island edge bracket or trace'
    end function
end module
