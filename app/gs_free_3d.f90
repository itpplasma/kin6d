module gs_free_3d_support
    ! Reference data of gs_free_3d: the 2D free-boundary result (gs_free
    ! mesh.dat, flux measured from psi_b, and the executed law), a bucket
    ! locator for it, and the frozen map-lane data pushed forward by the
    ! offset map T (offset start):
    !     y = T(x):  y = (R + dR, Z + dZ)  or  y_2 = Z_a + (Z - Z_a)(1 + e),
    !     chi(y) = psi(T^-1 y),  b3(y) = det(DT^-1) F(psi)/R(T^-1 y),
    !     p(y) = p(psi) - p_edge,
    ! so the 2D equilibrium itself is a stationary point of the map lane
    ! (X = T^-1) and the seed X = id is displaced from it by T.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_grad_shafranov, only: gs_state_t, gs_element_jet
    use kin6d_gs_profile_law, only: gs_profile_law_t, gs_law_primitives
    use kin6d_ideal_map, only: ideal_frozen_t
    implicit none
    private
    public :: reference_t, reference_locate_jet, offset_field_t
    public :: offset_forward, offset_backward

    type :: reference_t
        type(gs_state_t) :: mesh
        type(gs_profile_law_t) :: law
        real(dp) :: psi_axis_rel = 0, f_edge = 0, p_edge = 0, axis(2) = 0
        real(dp) :: lo(2) = 0, cell(2) = 1
        integer :: nb(2) = 1
        integer, allocatable :: first(:), list(:)
    end type

    type, extends(ideal_frozen_t) :: offset_field_t
        type(reference_t), pointer :: ref => null()
        real(dp) :: shift(2) = 0, elongation = 0
    contains
        procedure :: evaluate => offset_evaluate
    end type
contains
    pure function offset_forward(field, x) result(y)
        class(offset_field_t), intent(in) :: field
        real(dp), intent(in) :: x(2)
        real(dp) :: y(2)
        y = x + field%shift
        y(2) = field%ref%axis(2) + field%shift(2) &
            + (x(2) - field%ref%axis(2))*(1 + field%elongation)
    end function

    pure function offset_backward(field, y) result(x)
        class(offset_field_t), intent(in) :: field
        real(dp), intent(in) :: y(2)
        real(dp) :: x(2)
        x(1) = y(1) - field%shift(1)
        x(2) = field%ref%axis(2) + (y(2) - field%shift(2) - field%ref%axis(2)) &
            /(1 + field%elongation)
    end function

    subroutine reference_locate_jet(ref, x, psi, grad, status)
        ! Flux (from psi_b) and gradient of the 2D result at x.
        type(reference_t), intent(in) :: ref
        real(dp), intent(in) :: x(2)
        real(dp), intent(out) :: psi, grad(2)
        integer, intent(out) :: status
        real(dp) :: hess(2, 2)
        integer :: i, j, b, p
        status = 1; psi = 0; grad = 0
        i = int((x(1) - ref%lo(1))/ref%cell(1)) + 1
        j = int((x(2) - ref%lo(2))/ref%cell(2)) + 1
        if (i < 1 .or. j < 1 .or. i > ref%nb(1) .or. j > ref%nb(2)) return
        b = i + ref%nb(1)*(j - 1)
        do p = ref%first(b), ref%first(b + 1) - 1
            call gs_element_jet(ref%mesh, ref%list(p), x(1), x(2), psi, grad, hess, &
                status)
            if (status == 0) return
        end do
    end subroutine

    subroutine offset_evaluate(field, point, chi, chi_gradient, b3, pressure)
        class(offset_field_t), intent(in) :: field
        real(dp), intent(in) :: point(2)
        real(dp), intent(out) :: chi, chi_gradient(2), b3, pressure
        real(dp) :: x(2), grad(2), s, p, f
        integer :: st
        x = offset_backward(field, point)
        call reference_locate_jet(field%ref, x, chi, grad, st)
        if (st /= 0) error stop 'frozen point outside the 2D reference mesh'
        chi_gradient = [grad(1), grad(2)/(1 + field%elongation)]
        ! Profiles: clamp to the plasma side of the 2D boundary.
        s = max(0.0_dp, chi/field%ref%psi_axis_rel)
        call gs_law_primitives(field%ref%law, s*field%ref%psi_axis_rel, p, f, st)
        if (st /= 0) error stop 'profile law evaluation'
        b3 = f/x(1)/(1 + field%elongation)
        pressure = p - field%ref%p_edge
    end subroutine
end module gs_free_3d_support

program gs_free_3d_main
    ! N4 driver: 3D free boundary of the map lane against a 2D gs_free result
    ! (stage NL-3a). Reads &gs_free_3d from --config; writes summary-3d.dat and
    ! boundary-3d.dat in the working directory.
    use, intrinsic :: iso_fortran_env, only: dp => real64, output_unit
    use kin6d_grad_shafranov, only: gs_state_t, initialize_gs, read_gs_mesh
    use kin6d_gs_profile_law, only: read_gs_profile_law, gs_law_primitives
    use kin6d_gs_coils, only: gs_circuit_set_t, gs_read_coil_table, &
        gs_circuit_flux_vectors
    use kin6d_gs_free_boundary, only: gs_fb_trace_lcfs
    use kin6d_curved_triangle, only: initialize_periodic_curve, periodic_curve_t, &
        curve_value
    use kin6d_coils_3d, only: coil_set_t, coil_set_read, coil_read_currents, &
        coil_set_modes
    use kin6d_ideal_map, only: cylindrical_chart_t
    use kin6d_ideal_free_boundary, only: ideal_fb_t, ideal_fb_report_t, &
        ideal_fb_init, ideal_fb_set_drive, ideal_fb_solve_vacuum, ideal_fb_solve, &
        ideal_fb_solve_multipliers, ideal_fb_shell_field, ideal_fb_gradient, &
        ideal_fb_solve_shell_map, ideal_fb_rigid_stiffness, ideal_fb_mode_stability, &
        IDEAL_FB_OK
    use fortfem_axisymmetric_gs_kernel, only: gs_green_gradient
    use gs_free_3d_support, only: reference_t, reference_locate_jet, &
        offset_field_t, offset_forward, offset_backward
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp)
    integer, parameter :: max_circuits = 64
    type(reference_t), target :: ref
    type(offset_field_t) :: field
    type(gs_circuit_set_t) :: coils
    type(coil_set_t) :: rmp
    type(periodic_curve_t) :: gamma_curve
    type(gs_state_t) :: mesh
    type(cylindrical_chart_t) :: chart
    type(ideal_fb_t) :: fb
    type(ideal_fb_report_t) :: report
    real(dp), allocatable :: samples(:, :), contour(:, :), lcfs(:, :), chain(:, :)
    real(dp), allocatable :: rmp_currents(:)
    real(dp) :: mu0 = 4.0e-7_dp*acos(-1.0_dp), max_area = 0, offset_r = 0
    real(dp) :: offset_z = 0, elongation = 0, rmp_scale = 1, tolerance = 1.0e-9_dp
    real(dp) :: step_tolerance = 1.0e-5_dp, currents(max_circuits) = 0
    real(dp) :: gamma_c_scale = 0, residual_tolerance = 1.0e-9_dp, rigid(2, 2) = 0, seed_amplitude = 0
    real(dp) :: circulation, closure, perimeter, psi_b_abs, ip_2d, t0, t1
    integer :: boundary_nodes = 96, degree = 3, chain_nodes = 64, nfp = 1, ntor = 0
    integer :: nphi = 0, max_iterations = 40, rmp_nphi = 128, axis_sign = 1
    integer :: panels_per_edge = 0, st, i, nc, unit
    logical :: verbose = .true., start_at_target = .false., stability = .false.
    real(dp) :: stab(64) = 0, shifts(8) = [-1.0e-8_dp, -1.0e-7_dp, -1.0e-6_dp, &
        -1.0e-5_dp, -1.0e-4_dp, -1.0e-3_dp, -1.0e-2_dp, -1.0e-1_dp]
    character(512) :: mesh_file = 'mesh.dat', law_file = 'profile-law.dat'
    character(512) :: summary_file = 'summary.dat', coil_file = ''
    character(512) :: boundary_file = '', rmp_coil_file = '', rmp_current_file = ''
    character(512) :: config
    character(256) :: message
    namelist /gs_free_3d/ mesh_file, law_file, summary_file, coil_file, &
        boundary_file, boundary_nodes, max_area, degree, chain_nodes, mu0, nfp, ntor, &
        nphi, offset_r, offset_z, elongation, rmp_coil_file, rmp_current_file, &
        rmp_scale, rmp_nphi, tolerance, step_tolerance, max_iterations, axis_sign, &
        panels_per_edge, verbose, start_at_target, residual_tolerance, seed_amplitude, &
        stability, shifts, gamma_c_scale

    call get_command_argument(1, config)
    if (trim(config) /= '--config') error stop 'usage: gs_free_3d --config input.nml'
    call get_command_argument(2, config)
    open (newunit=unit, file=trim(config), status='old', action='read')
    read (unit, nml=gs_free_3d)
    close (unit)
    call cpu_time(t0)

    ! 2D reference.
    call read_gs_mesh(ref%mesh, trim(mesh_file), st)
    if (st /= 0) error stop 'cannot read the 2D mesh'
    call read_gs_profile_law(trim(law_file), ref%law, st)
    if (st /= 0) error stop 'cannot read the executed law'
    call read_summary()
    ref%psi_axis_rel = ref%law%axis
    call gs_law_primitives(ref%law, 0.0_dp, ref%p_edge, ref%f_edge, st)
    if (st /= 0) error stop 'edge law'
    call build_locator(ref)
    field%ref => ref
    field%shift = [offset_r, offset_z]; field%elongation = elongation

    ! Coils.
    call gs_read_coil_table(trim(coil_file), coils, st, message)
    if (st /= 0) error stop 'invalid coil table'
    nc = size(coils%incidence, 2)
    if (len_trim(rmp_coil_file) > 0) then
        call coil_set_read(rmp, trim(rmp_coil_file), st)
        if (st /= 0) error stop 'invalid RMP coil table'
        allocate (rmp_currents(rmp%ncoil))
        call coil_read_currents(trim(rmp_current_file), rmp%ncoil, rmp_currents, st)
        if (st /= 0) error stop 'invalid RMP currents'
        rmp_currents = rmp_scale*rmp_currents
    end if

    ! Gamma_c: the boundary_file curve, or (gamma_c_scale > 0) the 2D
    ! boundary scaled about the 2D axis, e.g. to keep in-vessel coils outside.
    if (gamma_c_scale > 0) then
        call gs_fb_trace_lcfs(ref%mesh, ref%mesh%psi, ref%axis, 0.0_dp, &
            real(axis_sign, dp), 1024, samples, st)
        if (st /= 0) error stop 'cannot trace the 2D boundary for Gamma_c'
        do i = 1, size(samples, 2)
            samples(:, i) = ref%axis + gamma_c_scale*(samples(:, i) - ref%axis)
        end do
    else
        call read_curve(trim(boundary_file), samples)
    end if
    call initialize_periodic_curve(gamma_curve, samples, st)
    if (st /= 0) error stop 'invalid Gamma_c'
    if (allocated(rmp_currents)) call check_coils_outside()
    allocate (contour(2, boundary_nodes))
    do i = 1, boundary_nodes
        call gamma(real(i - 1, dp)/boundary_nodes, contour(:, i))
    end do
    perimeter = sum(norm2(contour - cshift(contour, 1, dim=2), dim=1))
    if (max_area <= 0) max_area = 0.4_dp*(perimeter/boundary_nodes)**2
    call gs_fb_trace_lcfs(ref%mesh, ref%mesh%psi, ref%axis, 0.0_dp, &
        real(axis_sign, dp), chain_nodes, lcfs, st)
    if (st /= 0) error stop 'cannot trace the 2D boundary'
    allocate (chain(2, chain_nodes))
    do i = 1, chain_nodes
        chain(:, i) = offset_forward(field, lcfs(:, i))
    end do
    circulation = poloidal_circulation()
    write (output_unit, '(a,es14.6,a,es14.6)') 'circulation (T m) ', circulation, &
        '  -mu0 Ip ', -mu0*ip_2d
    call initialize_gs(mesh, contour, max_area, st, degree, curve=gamma, &
        interior_chain=chain)
    if (st /= 0) error stop 'mesh'
    if (panels_per_edge <= 0) panels_per_edge = merge(1, 2, degree == 2)
    call ideal_fb_init(fb, mesh, chain, chart, field, nfp, ntor, ref%f_edge, &
        circulation, st, mu0=mu0, panels_per_edge=panels_per_edge)
    if (st /= IDEAL_FB_OK) error stop 'free-boundary init'
    if (nphi > 0 .and. nphi /= fb%plasma%nphi) &
        write (output_unit, '(a)') 'note: nphi is set by ideal_map_init (3 ntor + 1)'
    fb%energy_scale = energy_scale()
    call ideal_fb_set_drive(fb, sampler, st, closure)
    write (output_unit, '(a,i0,a,i0,a,i0,a,es10.3)') 'plasma nodes ', &
        fb%plasma%nodes, ' shell nodes ', fb%ns, ' unknowns ', fb%n, &
        ' drive closure ', closure
    if (start_at_target) then
        ! Diagnostic: the exact target X = T^-1 in full coefficients.
        do i = 1, fb%plasma%nodes
            fb%plasma%x(1:2, 1, i) = offset_backward(field, fb%plasma%reference(:, i))
        end do
        call ideal_fb_solve_shell_map(fb, st)
        if (st /= IDEAL_FB_OK) error stop 'shell map projection'
    end if
    if (seed_amplitude /= 0 .and. ntor > 0) then
        ! Leakage seed: an interior n = nfp displacement (cos, R) weighted by
        ! the axis-expansion weight T (0 on Gamma_p, 1 at the axis node); the
        ! boundary is left axisymmetric. Without RMP drive it must decay.
        do i = 1, fb%plasma%nodes
            fb%plasma%x(1, 2, i) = seed_amplitude*fb%plasma%axis_t(i)
        end do
    end if
    call ideal_fb_solve_vacuum(fb, st)
    if (st /= IDEAL_FB_OK) error stop 'initial vacuum solve'
    write (output_unit, '(a,es10.3)') 'vacuum shell field vs 2D field, rel L2 ', &
        shell_mismatch()
    call xrow_norm('before multipliers')
    call ideal_fb_solve_multipliers(fb, st)
    if (st /= IDEAL_FB_OK) error stop 'initial multiplier solve'
    call xrow_norm('after multipliers')
    call ideal_fb_solve(fb, tolerance, max_iterations, report, st, verbose, &
        step_tolerance, residual_tolerance)
    call cpu_time(t1)
    block
        real(dp) :: f, k
        integer :: c, s2
        do c = 1, 2
            call ideal_fb_rigid_stiffness(fb, c, 1.0e-3_dp, f, k, s2)
            rigid(:, c) = [f, k]
            write (output_unit, '(a,a,a,es12.4,a,es12.4,a,i0)') 'rigid ', &
                trim(merge('R', 'Z', c == 1)), ' shift: force ', f, &
                ' stiffness (J/m^2/period) ', k, ' status ', s2
        end do
    end block
    ! Second variation per real mode (meaningful at an axisymmetric state).
    if (stability) then
        block
            real(dp) :: lam
            integer :: m, s2, j
            do m = 1, fb%nmode
                if (m > 1 .and. mod(m, 2) == 1) cycle
                stab(m) = huge(1.0_dp)
                do j = 0, size(shifts)
                    if (j > 0) then
                        if (.not. shifts(j) < 0) cycle
                        call ideal_fb_mode_stability(fb, m, lam, s2, shift=shifts(j))
                    else
                        call ideal_fb_mode_stability(fb, m, lam, s2)
                    end if
                    stab(m) = min(stab(m), lam)
                    write (output_unit, '(a,i0,a,i0,a,es10.2,a,es12.4,a,i0)') &
                        'stability mode ', m, ' (n=', (m/2)*nfp, ') shift ', &
                        merge(0.0_dp, shifts(max(j, 1)), j == 0), ' eigenvalue ', lam, &
                        ' status ', s2
                end do
            end do
        end block
    end if
    write (output_unit, '(a,i0,a,i0,a,es10.3)') 'newton status ', st, &
        ' iterations ', report%iterations, ' residual ', report%residual
    call write_results(st)
contains
    subroutine gamma(t, point)
        real(dp), intent(in) :: t
        real(dp), intent(out) :: point(2)
        call curve_value(gamma_curve, t, point)
    end subroutine

    subroutine check_coils_outside()
        ! The formulation needs every coil outside Gamma_c (coil fields enter
        ! only through their traces): reject RMP polyline points inside it,
        ! and report the smallest clearance.
        real(dp) :: r, z, clear
        integer :: k, i, j, inside, n
        n = size(samples, 2)
        clear = huge(1.0_dp)
        do k = 1, rmp%nseg
            r = hypot(rmp%a(1, k), rmp%a(2, k)); z = rmp%a(3, k)
            inside = 0
            j = n
            do i = 1, n
                if ((samples(2, i) > z) .neqv. (samples(2, j) > z)) then
                    if (r < (samples(1, j) - samples(1, i))*(z - samples(2, i)) &
                        /(samples(2, j) - samples(2, i)) + samples(1, i)) &
                        inside = 1 - inside
                end if
                j = i
            end do
            if (inside == 1) error stop 'RMP coil inside Gamma_c'
            clear = min(clear, minval(hypot(samples(1, :) - r, samples(2, :) - z)))
        end do
        write (output_unit, '(a,f8.4)') 'RMP coil clearance from Gamma_c (m) ', clear
    end subroutine

    subroutine read_curve(path, points)
        character(*), intent(in) :: path
        real(dp), allocatable, intent(out) :: points(:, :)
        integer :: u, n, k
        open (newunit=u, file=path, status='old', action='read')
        read (u, *) n
        allocate (points(2, n))
        do k = 1, n
            read (u, *) points(:, k)
        end do
        close (u)
    end subroutine

    subroutine read_summary()
        character(256) :: line, key
        real(dp) :: value
        integer :: u, ios, k
        open (newunit=u, file=trim(summary_file), status='old', action='read')
        do
            read (u, '(a)', iostat=ios) line
            if (ios /= 0) exit
            read (line, *, iostat=ios) key, value
            if (ios /= 0) cycle
            select case (trim(key))
            case ('axis_R'); ref%axis(1) = value
            case ('axis_Z'); ref%axis(2) = value
            case ('psi_bound'); psi_b_abs = value
            case ('plasma_current'); ip_2d = value
            end select
            ! The 2D result's circuit currents (A) drive the 3D solve.
            if (key(1:16) == 'circuit_current_') then
                read (key(17:), *) k
                currents(k) = value
            end if
        end do
        close (u)
    end subroutine

    subroutine build_locator(r)
        ! Uniform buckets over the 2D mesh: every cell is listed in the
        ! buckets its node bounding box (padded) touches.
        type(reference_t), intent(inout) :: r
        integer, allocatable :: count(:)
        real(dp) :: lo(2), hi(2), cmin(2), cmax(2)
        integer :: e, i, j, nbk, pass, i0, i1, j0, j1, total
        lo = minval(r%mesh%points, dim=2); hi = maxval(r%mesh%points, dim=2)
        nbk = max(1, int(sqrt(real(size(r%mesh%triangles, 2), dp)/2)))
        r%nb = nbk
        r%lo = lo - 1.0e-9_dp*(hi - lo)
        r%cell = (hi - lo)*(1 + 2.0e-9_dp)/nbk
        allocate (count(nbk*nbk + 1))
        do pass = 1, 2
            count = 0
            if (pass == 2) total = 0
            do e = 1, size(r%mesh%triangles, 2)
                cmin = minval(r%mesh%points(:, r%mesh%triangles(:, e)), dim=2)
                cmax = maxval(r%mesh%points(:, r%mesh%triangles(:, e)), dim=2)
                cmin = cmin - 0.05_dp*(cmax - cmin); cmax = cmax + 0.05_dp*(cmax - cmin)
                i0 = max(1, int((cmin(1) - r%lo(1))/r%cell(1)) + 1)
                i1 = min(nbk, int((cmax(1) - r%lo(1))/r%cell(1)) + 1)
                j0 = max(1, int((cmin(2) - r%lo(2))/r%cell(2)) + 1)
                j1 = min(nbk, int((cmax(2) - r%lo(2))/r%cell(2)) + 1)
                do j = j0, j1
                    do i = i0, i1
                        count(i + nbk*(j - 1)) = count(i + nbk*(j - 1)) + 1
                        if (pass == 2) then
                            r%list(r%first(i + nbk*(j - 1)) + count(i + nbk*(j - 1)) &
                                - 1) = e
                        end if
                    end do
                end do
            end do
            if (pass == 1) then
                allocate (r%first(nbk*nbk + 1), r%list(sum(count)))
                r%first(1) = 1
                do i = 1, nbk*nbk
                    r%first(i + 1) = r%first(i) + count(i)
                end do
            end if
        end do
        if (total < 0) continue
    end subroutine

    real(dp) function poloidal_circulation() result(c)
        ! ccw circulation of the 2D field on Gamma_c pulled 10 % toward the
        ! axis (in vacuum, enclosing the plasma).
        integer, parameter :: n = 2048
        real(dp) :: x(2), xp(2), psi, grad(2), xm(2)
        integer :: k, s
        real(dp) :: loop
        loop = 0.9_dp
        if (gamma_c_scale > 0) loop = (1 + gamma_c_scale)/(2*gamma_c_scale)
        c = 0
        do k = 1, n
            call gamma(real(k - 1, dp)/n, x)
            call gamma(real(k, dp)/n, xp)
            x = ref%axis + loop*(x - ref%axis); xp = ref%axis + loop*(xp - ref%axis)
            xm = 0.5_dp*(x + xp)
            call reference_locate_jet(ref, xm, psi, grad, s)
            if (s /= 0) error stop 'circulation loop outside the 2D mesh'
            c = c + dot_product([-grad(2), grad(1)]/xm(1), xp - x)
        end do
    end function

    subroutine sampler(r, z, k, bcos, bsin)
        ! PF/CS circuits and the unit-circulation filament at the 2D axis
        ! (n = 0), RMP coils (all n in nfp Z, including any n = 0 part).
        real(dp), intent(in) :: r, z
        integer, intent(in) :: k
        real(dp), intent(out) :: bcos(3), bsin(3)
        real(dp) :: g(1, max_circuits), gr(1, max_circuits), gz(1, max_circuits)
        real(dp) :: gf, gg(4)
        real(dp), save :: last(2) = -1
        complex(dp), allocatable, save :: modes(:, :, :)
        integer, allocatable :: list(:)
        integer :: j
        bcos = 0; bsin = 0
        if (k == 0) then
            call gs_circuit_flux_vectors(coils, [r], [z], g(:, :nc), gr(:, :nc), &
                gz(:, :nc))
            bcos(1) = -dot_product(gz(1, :nc), currents(:nc))/r
            bcos(3) = dot_product(gr(1, :nc), currents(:nc))/r
            call gs_green_gradient(r, z, ref%axis(1), ref%axis(2), gf, gg)
            bcos([1, 3]) = bcos([1, 3]) + circulation*[gg(2), -gg(1)]/r
        end if
        if (.not. allocated(rmp_currents)) return
        if (any([r, z] /= last) .or. .not. allocated(modes)) then
            list = [(j*nfp, j=0, ntor)]
            if (allocated(modes)) deallocate (modes)
            allocate (modes(3, ntor + 1, 1))
            call coil_set_modes(rmp, rmp_currents, reshape([r, z], [2, 1]), rmp_nphi, &
                list, modes, j)
            if (j /= 0) error stop 'RMP coil modes'
            last = [r, z]
        end if
        if (k == 0) then
            bcos = bcos + real(modes(:, 1, 1))
        else
            bcos = bcos + 2*real(modes(:, k + 1, 1))
            bsin = bsin - 2*aimag(modes(:, k + 1, 1))
        end if
    end subroutine

    subroutine xrow_norm(label)
        character(*), intent(in) :: label
        real(dp), allocatable :: g(:)
        integer :: s2
        allocate (g(fb%n))
        call ideal_fb_gradient(fb, g, s2)
        write (output_unit, '(a,a,4es11.3)') label, ' |g| p x l f ', &
            norm2(g(1:fb%np)), norm2(g(fb%np + 1:fb%np + fb%nx)), &
            norm2(g(fb%np + fb%nx + 1:fb%np + 2*fb%nx)), norm2(g(fb%np + 2*fb%nx + 1:))
    end subroutine

    real(dp) function shell_mismatch() result(rel)
        ! Shell field of the vacuum solve against the 2D field (B_pol from
        ! psi, B_phi = F_vac/R); meaningful for the unshifted start, where
        ! Gamma_p is the 2D boundary.
        real(dp) :: x(3), b(3), psi, grad(2), w, e2, r2, xi(2)
        real(dp), parameter :: pts(2, 3) = reshape([1.0_dp/6, 1.0_dp/6, &
            2.0_dp/3, 1.0_dp/6, 1.0_dp/6, 2.0_dp/3], [2, 3])
        integer :: c, g, s
        e2 = 0; r2 = 0
        do c = 1, fb%ncell
            do g = 1, 3
                xi = pts(:, g)
                call ideal_fb_shell_field(fb, c, xi, 0.0_dp, x, b, s)
                call reference_locate_jet(ref, x(1:2), psi, grad, s)
                if (s /= 0) cycle
                w = x(1)
                e2 = e2 + w*((b(1) + grad(2)/x(1))**2 + (b(3) - grad(1)/x(1))**2)
                r2 = r2 + w*sum(grad**2)/x(1)**2
            end do
        end do
        rel = sqrt(e2/r2)
    end function

    real(dp) function energy_scale() result(e)
        ! sigma V = int (|B_pol|^2/mu0 + p) dV per field period, from the
        ! frozen data on the plasma reference (toroidal field removed).
        real(dp) :: y(2), chi, grad(2), b3, p, dg(2, 2)
        integer :: c, g
        e = 0
        do c = 1, fb%plasma%cells
            do g = 1, fb%plasma%nq
                y = matmul(fb%plasma%reference(:, fb%plasma%conn(:, c)), &
                    fb%plasma%shape(:, g))
                dg = matmul(fb%plasma%reference(:, fb%plasma%conn(:, c)), &
                    transpose(fb%plasma%dshape(:, :, g)))
                call field%evaluate(y, chi, grad, b3, p)
                e = e + fb%plasma%weight(g)*abs(dg(1, 1)*dg(2, 2) - dg(1, 2)*dg(2, 1)) &
                    *y(1)*(sum(grad**2)/(y(1)**2*mu0) + p)
            end do
        end do
        e = e*2*pi/nfp
    end function

    subroutine write_results(status)
        ! Return distances of Gamma_p from the 2D boundary (signed normal
        ! distance psi/|grad psi| of the 2D result at the mapped nodes, every
        ! phi plane), interior surface defect, axis position and the n /= 0
        ! normal displacement of Gamma_p per harmonic.
        integer, intent(in) :: status
        real(dp) :: x(2), psi, grad(2), d, dmax, drms, defect, chi, cg(2), b3, p
        real(dp) :: axis_x(2), amp(0:16), nrm(2), phi, interior(16)
        integer :: a, j, m, k, s, u, cnt
        dmax = 0; drms = 0; cnt = 0; amp = 0
        open (newunit=u, file='boundary-3d.dat', status='replace', action='write')
        write (u, '(a)') '# phi_hat R Z distance_to_2D_boundary (m), Gamma_p nodes'
        do j = 1, fb%plasma%nphi
            phi = fb%plasma%phi(j)
            do a = 1, fb%plasma%nodes
                if (.not. fb%plasma%boundary(a)) cycle
                x = 0
                do m = 1, fb%plasma%nmode
                    x = x + fb%plasma%x(1:2, m, a)*fb%plasma%trig(m, j)
                end do
                call reference_locate_jet(ref, x, psi, grad, s)
                if (s /= 0) cycle
                d = psi/norm2(grad)
                dmax = max(dmax, abs(d)); drms = drms + d**2; cnt = cnt + 1
                write (u, '(4es18.9)') phi, x, d
            end do
        end do
        close (u)
        drms = sqrt(drms/max(cnt, 1))
        defect = 0
        do a = 1, fb%plasma%nodes
            x = fb%plasma%x(1:2, 1, a)
            call reference_locate_jet(ref, x, psi, grad, s)
            if (s /= 0) cycle
            call field%evaluate(fb%plasma%reference(:, a), chi, cg, b3, p)
            defect = max(defect, abs(psi - chi)/norm2(grad))
        end do
        axis_x = fb%plasma%x(1:2, 1, fb%plasma%axis_node)
        block
            real(dp) :: mean(2)
            integer :: nb
            mean = 0; nb = 0
            do a = 1, fb%plasma%nodes
                if (.not. fb%plasma%boundary(a)) cycle
                mean = mean + fb%plasma%x(1:2, 1, a) - fb%plasma%reference(:, a)
                nb = nb + 1
            end do
            write (output_unit, '(a,2es11.3,a,2es11.3)') 'mean Gamma_p shift ', &
                mean/nb, '  axis node shift ', axis_x &
                - fb%plasma%reference(:, fb%plasma%axis_node)
        end block
        do a = 1, fb%plasma%nodes
            if (.not. fb%plasma%boundary(a)) cycle
            nrm = fb%plasma%direction(:, a)
            do k = 1, ntor
                amp(k) = max(amp(k), sqrt(dot_product(nrm, fb%plasma%x(1:2, 2*k, a))**2 &
                    + dot_product(nrm, fb%plasma%x(1:2, 2*k + 1, a))**2))
            end do
        end do
        do k = 1, ntor
            interior(k) = maxval(abs(fb%plasma%x(1:2, 2*k:2*k + 1, :)))
        end do
        open (newunit=u, file='summary-3d.dat', status='replace', action='write')
        write (u, '(a,1x,i0)') 'status', status
        write (u, '(a,1x,i0)') 'iterations', report%iterations
        write (u, '(a,1x,es16.8)') 'residual_sigma', report%residual
        write (u, '(a,1x,es16.8)') 'merit_ratio', sqrt(report%merit/report%merit0)
        write (u, '(a,1x,es16.8)') 'last_step_m', report%step
        write (u, '(a,1x,es16.8)') 'constraint', report%constraint
        write (u, '(a,1x,l1)') 'roundoff_floor', report%floor
        write (u, '(a,2(1x,es16.8))') 'rigid_R_force_stiffness', rigid(:, 1)
        write (u, '(a,2(1x,es16.8))') 'rigid_Z_force_stiffness', rigid(:, 2)
        if (stability) then
            do k = 1, fb%nmode
                if (k > 1 .and. mod(k, 2) == 1) cycle
                write (u, '(a,i0,1x,es16.8)') 'stability_min_lambda_mode', k, stab(k)
            end do
        end if
        write (u, '(a,1x,es16.8)') 'offset_r', offset_r
        write (u, '(a,1x,es16.8)') 'offset_z', offset_z
        write (u, '(a,1x,es16.8)') 'elongation', elongation
        write (u, '(a,1x,es16.8)') 'boundary_distance_max_m', dmax
        write (u, '(a,1x,es16.8)') 'boundary_distance_rms_m', drms
        write (u, '(a,1x,es16.8)') 'surface_defect_max_m', defect
        write (u, '(a,2(1x,es16.8))') 'axis_node_position', axis_x
        write (u, '(a,2(1x,es16.8))') 'axis_2d', ref%axis
        write (u, '(a,1x,es16.8)') 'circulation', circulation
        write (u, '(a,1x,es16.8)') 'ip_2d', ip_2d
        write (u, '(a,1x,es16.8)') 'f_vac', ref%f_edge
        write (u, '(a,1x,es16.8)') 'drive_closure', closure
        write (u, '(a,1x,es16.8)') 'rmp_scale', rmp_scale
        do k = 1, ntor
            write (u, '(a,i0,1x,es16.8)') 'gamma_p_normal_amplitude_n', k*nfp, amp(k)
            write (u, '(a,i0,1x,es16.8)') 'plasma_max_coefficient_n', k*nfp, interior(k)
        end do
        write (u, '(a,1x,i0)') 'unknowns', fb%n
        write (u, '(a,1x,es16.8)') 'assemble_s', report%assemble_s
        write (u, '(a,1x,es16.8)') 'factor_s', report%factor_s
        write (u, '(a,1x,es16.8)') 'cpu_s', t1 - t0
        close (u)
        write (output_unit, '(a,es10.3,a,es10.3,a,es10.3)') 'boundary distance max ', &
            dmax, ' rms ', drms, ' surface defect ', defect
        write (output_unit, '(a,2f12.6,a,2f12.6)') 'axis node ', axis_x, ' 2D axis ', &
            ref%axis
        do k = 1, ntor
            write (output_unit, '(a,i0,a,es12.4,a,es12.4)') 'Gamma_p normal amplitude n=', &
                k*nfp, ' (m) ', amp(k), '  plasma max coefficient ', interior(k)
        end do
    end subroutine
end program gs_free_3d_main
