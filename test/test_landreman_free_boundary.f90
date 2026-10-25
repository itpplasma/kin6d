module landreman_circuit_basis
    ! Landreman drive as circuits: ring source k (unit current, the
    ! fortfem ring Green function, as drive%external), then the toroidal
    ! field B_phi = 1/R (unit F), then optional coils_3d circuits.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortfem_axisymmetric_gs_kernel, only: gs_green_gradient
    use kin6d_ideal_free_boundary, only: ideal_fb_field_basis_t
    use kin6d_ideal_fb_circuits, only: coil_circuit_basis_t
    implicit none
    private
    type, extends(ideal_fb_field_basis_t), public :: ring_basis_t
        real(dp), allocatable :: source(:, :)
        type(coil_circuit_basis_t), allocatable :: coils
    contains
        procedure :: nfield => ring_nfield
        procedure :: sample => ring_sample
    end type
contains
    integer function ring_nfield(self)
        class(ring_basis_t), intent(in) :: self
        ring_nfield = size(self%source, 2) + 1
        if (allocated(self%coils)) ring_nfield = ring_nfield + self%coils%ncircuit
    end function

    subroutine ring_sample(self, rz, nfp, ntor, bcos, bsin, status)
        class(ring_basis_t), intent(in) :: self
        real(dp), intent(in) :: rz(:, :)
        integer, intent(in) :: nfp, ntor
        real(dp), intent(out) :: bcos(:, 0:, :, :), bsin(:, 0:, :, :)
        integer, intent(out) :: status
        real(dp) :: g, gg(4), r, z
        integer :: q, k, ns
        ns = size(self%source, 2)
        bcos = 0; bsin = 0
        do q = 1, size(rz, 2)
            r = rz(1, q); z = rz(2, q)
            do k = 1, ns
                call gs_green_gradient(r, z, self%source(1, k), self%source(2, k), &
                                       g, gg)
                bcos(:, 0, k, q) = [-gg(2)/r, 0.0_dp, gg(1)/r]
            end do
            bcos(2, 0, ns + 1, q) = 1/r
        end do
        status = 0
        if (allocated(self%coils)) call self%coils%sample(rz, nfp, ntor, &
            bcos(:, :, ns + 2:, :), bsin(:, :, ns + 2:, :), status)
    end subroutine
end module landreman_circuit_basis

program test_landreman_free_boundary
    use, intrinsic :: iso_fortran_env, only: dp => real64, compiler_version, &
                                                                        compiler_options
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use kin6d_landreman_leaves, only: landreman_axisymmetric, landreman_surface, &
                                landreman_project_surface, landreman_axisymmetric_jvp, &
                                      landreman_axisymmetric_vjp
    use kin6d_landreman_benchmark, only: landreman_frozen_t, landreman_drive_t, &
                                         landreman_volume_boundary, LANDREMAN_DELTA
    use kin6d_grad_shafranov, only: gs_state_t, initialize_gs, gs_element_jet
    use kin6d_gs_coils, only: gs_circuit_set_t, gs_coil_filaments
    use kin6d_gs_profile_law, only: gs_profile_law_t
    use kin6d_gs_free_boundary, only: gs_fb_system_t, gs_fb_config_t, &
             gs_fb_eval_t, gs_fb_history_t, gs_fb_setup, gs_fb_newton, gs_fb_evaluate, &
                                      gs_fb_locate, GS_FB_TARGET_MEAN
    use fortfem_axisymmetric_gs_kernel, only: gs_mu0, gs_green_gradient
    use kin6d_equilibrium_chart, only: cylindrical_chart_t
    use kin6d_ideal_free_boundary, only: ideal_fb_t, ideal_fb_init, &
             ideal_fb_set_drive, ideal_fb_solve, ideal_fb_report_t, ideal_fb_shell_field
    use kin6d_ideal_fb_circuits, only: ideal_fb_circuits_t, ideal_fb_circuits_build, &
             ideal_fb_set_currents, coil_circuit_basis_t, coil_circuit_basis_init
    use kin6d_coils_3d, only: coil_set_t, coil_set_from_polylines, coil_mu0
    use landreman_circuit_basis, only: ring_basis_t
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp), delta = LANDREMAN_DELTA
    type(landreman_drive_t) :: drive
    real(dp) :: distance(2), field_error(2), pressure_error(2)
    integer :: st, l
    ! Modes: (none) set_drive with the fitted sampler; circuits: the same
    ! drive through circuit currents (kin6d_ideal_fb_circuits); helical:
    ! nfp = 3, N = 1 with an added 3-fold window-coil circuit.
    character(16) :: mode
    logical :: use_circuits
    integer :: case_nfp = 1, case_ntor = 0
    real(dp) :: window_current = 0, helical_amplitude = 0
    mode = 'sampler'
    if (command_argument_count() > 0) call get_command_argument(1, mode)
    if (all(trim(mode) /= [character(16) :: 'sampler', 'circuits', 'helical'])) &
        error stop 'unknown mode'
    use_circuits = trim(mode) /= 'sampler'
    print '(a)', compiler_version()
    print '(a)', compiler_options()
    print '(a)', 'mu0=1; epsilon=0; delta=0.01; sources=24; fit=96; quad=48/96'
    call local_checks()
    call derivative_checks()
    call drive%fit(24, 96, 48, st)
    if (st /= 0) error stop 'external source fit'
    call drive_checks()
    if (trim(mode) == 'helical') then
        call helical_case()
        stop
    end if
    do l = 1, 2
        call solve_case(16*2**(l - 1), distance(l), field_error(l), pressure_error(l))
    end do
    ! Analytic interpolation scaling: P3 geometry O(h^4), B O(h^3).
    ! Tolerances are engineering budgets, not a posteriori certificates:
    ! 1% of minor radius, poloidal field scale 2sqrt(delta), pressure 2delta.
    if (distance(2) > 0.01_dp*sqrt(delta)) error stop 'surface return budget'
    if (field_error(2) > 0.01_dp*2*sqrt(delta)) error stop 'field budget'
    if (pressure_error(2) > 0.01_dp*2*delta) error stop 'pressure budget'
    if (distance(2) > 0.8_dp*distance(1)) error stop 'surface refinement'
    if (field_error(2) > 0.8_dp*field_error(1)) error stop 'field refinement'
    if (pressure_error(2) > 0.8_dp*pressure_error(1)) error stop 'pressure refinement'
    print '(a)', 'candidate; numerical evidence only; continuous bounds unavailable'
contains
    subroutine exact(r, z, b, j, psi, p, gp)
        real(dp), intent(in) :: r, z
        real(dp), intent(out) :: b(3), j(3), psi, p, gp(2)
        call landreman_axisymmetric(r, z, delta, psi, p, b(1), b(2), b(3), &
                                    j(1), j(2), j(3), gp(1), gp(2))
    end subroutine

    subroutine local_checks()
        ! Independent explicit curl and equilibrium identities, and FD oracle
        ! for derivatives of the emitted primal. mu0=1.
        real(dp) :: r, z, b(3), j(3), psi, p, gp(2), h, brp(3), brm(3), bzp(3), bzm(3)
        real(dp) :: junk(3), g(2), s, t, err, div, curl(3), force(3), tr(10), tz(10)
        integer :: i
        h = 1.0e-5_dp; err = 0
        do i = 1, 32
            r = 1 + 0.07_dp*cos(real(i, dp)); z = 0.04_dp*sin(2.0_dp*i)
            call exact(r, z, b, j, psi, p, gp)
            force = [j(2)*b(3) - j(3)*b(2), j(3)*b(1) - j(1)*b(3), &
                     j(1)*b(2) - j(2)*b(1)]
            if (maxval(abs(force + [2*gp(1), 0.0_dp, 2*gp(2)])) > 5e-14_dp) &
                error stop 'force balance'
            if (abs(b(1)*gp(1) + b(3)*gp(2)) > 5e-14_dp) error stop 'flux invariant'
            call exact(r + h, z, brp, junk, s, t, g)
            call exact(r - h, z, brm, junk, s, t, g)
            call exact(r, z + h, bzp, junk, s, t, g)
            call exact(r, z - h, bzm, junk, s, t, g)
            div = (brp(1) - brm(1) + bzp(3) - bzm(3))/(2*h) + b(1)/r
            curl = [-(bzp(2) - bzm(2))/(2*h), &
                    (bzp(1) - bzm(1) - brp(3) + brm(3))/(2*h), &
                    (brp(2) - brm(2))/(2*h) + b(2)/r]
            err = max(err, abs(div), maxval(abs(curl - j)))
            call tangent_values([r, z, delta], [1.0_dp, 0.0_dp, 0.0_dp], tr)
            call tangent_values([r, z, delta], [0.0_dp, 1.0_dp, 0.0_dp], tz)
            if (abs(tr(3) + b(1)/r + tz(5)) > 5e-14_dp) &
                error stop 'generated divergence to roundoff'
            ! Closed-form independent cylindrical current.
            if (maxval(abs(j - [4*z/(r*sqrt(1 - 4*psi)), 2/r + 2*r, &
                        -2*(r*r - 1)/sqrt(1 - 4*psi)])) > 5e-14_dp) error stop 'current'
        end do
        print '(a,es12.4)', 'field FD curl/div residual ', err
        if (err > 2e-8_dp) error stop 'independent curl/div'
    end subroutine

    subroutine values(x, v)
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: v(10)
        call landreman_axisymmetric(x(1), x(2), x(3), v(1), v(2), v(3), v(4), &
                                    v(5), v(6), v(7), v(8), v(9), v(10))
    end subroutine

    subroutine tangent_values(x, v, tv)
        real(dp), intent(in) :: x(3), v(3)
        real(dp), intent(out) :: tv(10)
        call landreman_axisymmetric_jvp(x(1), x(2), x(3), v(1), v(2), v(3), &
                  tv(1), tv(2), tv(3), tv(4), tv(5), tv(6), tv(7), tv(8), tv(9), tv(10))
    end subroutine

    subroutine derivative_checks()
        real(dp) :: x(3), v(3), bar(10), tangent(10), adj(3), vp(10), vm(10), h, err
        integer :: i
        x = [1.03_dp, 0.02_dp, delta]; v = [0.13_dp, -0.23_dp, 0.01_dp]
        bar = [(sin(real(i, dp)), i=1, 10)]; h = 1e-5_dp
        call tangent_values(x, v, tangent)
        call landreman_axisymmetric_vjp(x(1), x(2), x(3), bar(1), bar(2), bar(3), &
                              bar(4), bar(5), bar(6), bar(7), bar(8), bar(9), bar(10), &
                                        adj(1), adj(2), adj(3))
        call values(x + h*v, vp); call values(x - h*v, vm)
        err = maxval(abs((vp - vm)/(2*h) - tangent))
        if (err > 1e-8_dp) error stop 'generated JVP FD'
        if (abs(dot_product(bar, tangent) - dot_product(adj, v)) > 1e-13_dp) &
            error stop 'Euclidean leaf JVP/VJP pairing'
        print '(a,es12.4)', 'generated JVP FD error ', err
    end subroutine

    subroutine drive_checks()
        real(dp) :: theta, r, z, dr, dz, b(3), j(3), psi, p, gp(2), bj(3), bj2(3), bx(3)
        real(dp) :: ip, ip2, quad, flux, traction, jump, circ, extcirc, h
       real(dp) :: bp(3), bm(3), zp(3), zm(3), vacerr, scale, der(3, 2), normal(3), d(3)
        integer :: i, l
        quad = 0; flux = 0; traction = 0; jump = 0; circ = 0; extcirc = 0
        do i = 1, 96
            theta = 2*pi*(i - 0.371_dp)/96
            call landreman_surface(theta, delta, r, z, dr, dz)
            call exact(r, z, b, j, psi, p, gp)
            call landreman_volume_boundary(theta, delta, 48, bj, ip)
            call landreman_volume_boundary(theta, delta, 96, bj2, ip2)
            quad = max(quad, norm2(bj - bj2))
            call drive%external(r, z, bx)
            if (abs(r*b(2) - sqrt(1 - 4*delta)) > 1e-14_dp) &
                error stop 'interior toroidal period'
            if (abs(r*bx(2) - sqrt(1 - 4*delta)) > 1e-14_dp) &
                error stop 'external toroidal period (lk=0)'
            d = bj2 + bx - b; normal = [gp(1), 0.0_dp, gp(2)]/norm2(gp)
            flux = max(flux, abs(dot_product(d, normal)))
            jump = max(jump, norm2(d - dot_product(d, normal)*normal))
            traction = max(traction, abs(p + (sum(b*b) - sum((bj2 + bx)**2))/2))
            circ = circ + 2*pi*(b(1)*dr + b(3)*dz)/96
            extcirc = extcirc + 2*pi*(bx(1)*dr + bx(3)*dz)/96
            if (abs(p) > 1e-15_dp) error stop 'edge pressure'
        end do
        print '(a,5es12.4)', 'fit, quadrature, flux, tangent jump, traction ', &
            drive%fit_error, quad, flux, jump, traction
        print '(a,f12.8)', 'F from both toroidal periods (lk=0) ', sqrt(1 - 4*delta)
        print '(a,3es12.4)', 'c=-Ip, boundary circulation, external circulation ', &
            drive%circulation, circ, extcirc
        ! For coaxial longitude gamma and axis, the Gauss linking integrand
        ! vanishes after a disjoint isotopy to coplanar coaxial circles.
        ! A resolved Gauss sum below also guards the actual cycle choice.
        call linking_check()
        if (max(quad, flux, jump) > 2e-5_dp) error stop 'virtual casing defect'
        if (traction > 2e-5_dp) error stop 'pressure balance defect'
        if (abs(circ - drive%circulation) > 1e-10_dp) error stop 'Ampere sign/period'
        if (abs(extcirc) > 1e-10_dp) error stop 'external poloidal period'
        h = 1e-5_dp
        do l = 1, 3
            scale = 1 + 0.2_dp*l; vacerr = 0
            do i = 1, 32
                theta = 2*pi*(i - 0.2_dp)/32
                call landreman_surface(theta, scale*scale*delta, r, z, dr, dz)
                call drive%external(r, z, b)
                call drive%external(r + h, z, bp); call drive%external(r - h, z, bm)
                call drive%external(r, z + h, zp); call drive%external(r, z - h, zm)
                der(:, 1) = (bp - bm)/(2*h); der(:, 2) = (zp - zm)/(2*h)
                vacerr = max(vacerr, abs(der(1, 1) + b(1)/r + der(3, 2)), &
                             abs(der(1, 2) - der(3, 1)), abs(der(2, 2)), &
                             abs(der(2, 1) + b(2)/r))
            end do
            print '(a,f4.1,a,es12.4)', 'G4 numerical vacuum to shell rho/a=', &
                scale, ' curl/div FD ', vacerr
            if (vacerr > 2e-7_dp) error stop 'vacuum continuation residual'
        end do
    end subroutine

    subroutine linking_check()
        real(dp) :: a, b, x(3), y(3), dx(3), dy(3), cross(3), s, total, r, z, dr, dz
        integer :: i, j
        call landreman_surface(0.3_dp, delta, r, z, dr, dz)
        total = 0
        do i = 1, 512
            a = 2*pi*(i - 0.5_dp)/512
            x = [r*cos(a), r*sin(a), z]; dx = [-r*sin(a), r*cos(a), 0.0_dp]
            do j = 1, 512
                b = 2*pi*(j - 0.2_dp)/512
                y = [cos(b), sin(b), 0.0_dp]; dy = [-sin(b), cos(b), 0.0_dp]
                cross = [dx(2)*dy(3) - dx(3)*dy(2), dx(3)*dy(1) - dx(1)*dy(3), &
                         dx(1)*dy(2) - dx(2)*dy(1)]
                s = norm2(x - y)
                total = total + dot_product(x - y, cross)/s**3*pi/512**2
            end do
        end do
        print '(a,es12.4)', 'longitude-axis Gauss linking number ', total
        if (abs(total) > 1e-8_dp) error stop 'linking condition'
    end subroutine

    subroutine gamma_c(t, point)
        real(dp), intent(in) :: t
        real(dp), intent(out) :: point(2)
        real(dp) :: dr, dz
        call landreman_surface(2*pi*t, 2.56_dp*delta, point(1), point(2), dr, dz)
    end subroutine

    subroutine sampler(r, z, k, bc, bs)
        real(dp), intent(in) :: r, z
        integer, intent(in) :: k
        real(dp), intent(out) :: bc(3), bs(3)
        call drive%sample(r, z, k, bc, bs)
    end subroutine

    subroutine filament(r, z, k, bc, bs)
        ! Fixed part of the drive: the unit-circulation filament times c.
        real(dp), intent(in) :: r, z
        integer, intent(in) :: k
        real(dp), intent(out) :: bc(3), bs(3)
        real(dp) :: g, gg(4)
        bc = 0; bs = 0
        if (k /= 0) return
        call gs_green_gradient(r, z, 1.0_dp, 0.0_dp, g, gg)
        bc(1) = drive%circulation*gg(2)/r
        bc(3) = -drive%circulation*gg(1)/r
    end subroutine

    subroutine windows(basis)
        ! Six rectangular window coils at R = 1.3, phi = 0, 60, ..., 300 deg,
        ! half-width 0.2 rad, Z in [-0.15, 0.15], one circuit with alternating
        ! winding signs: the field has n = 3(2j + 1) only. (Same-sign windows
        ! add an n = 0 radial field that drives the vertical mode.)
        ! mu0 = 1 here: field_scale removes the SI mu0 of coils_3d.
        type(coil_circuit_basis_t), intent(out) :: basis
        type(coil_set_t) :: coils
        real(dp) :: pts(3, 24), corner(2, 4), phi0
        integer :: c, j, st
        corner = reshape([-0.2_dp, -0.15_dp, 0.2_dp, -0.15_dp, 0.2_dp, 0.15_dp, &
            -0.2_dp, 0.15_dp], [2, 4])
        do c = 1, 6
            phi0 = 2*pi*(c - 1)/6
            do j = 1, 4
                pts(:, 4*(c - 1) + j) = [1.3_dp*cos(phi0 + corner(1, j)), &
                    1.3_dp*sin(phi0 + corner(1, j)), corner(2, j)]
            end do
        end do
        call coil_set_from_polylines(coils, [4, 4, 4, 4, 4, 4], pts, status=st)
        if (st /= 0) error stop 'window coils'
        call coil_circuit_basis_init(basis, coils, [1, 1, 1, 1, 1, 1], &
            [1.0_dp, -1.0_dp, 1.0_dp, -1.0_dp, 1.0_dp, -1.0_dp], 1, 128, st, &
            field_scale=1/coil_mu0)
        if (st /= 0) error stop 'window circuit'
    end subroutine

    subroutine circuit_drive(fb)
        ! Replace the sampler traces by the circuit-current path:
        ! currents = (fitted ring currents, F(delta)[, window current]).
        ! The axisymmetric part must reproduce the sampler traces to roundoff.
        type(ideal_fb_t), intent(inout) :: fb
        type(ring_basis_t) :: basis
        type(ideal_fb_circuits_t) :: circuits
        real(dp), allocatable :: d(:, :), l(:, :), currents(:)
        real(dp) :: err
        integer :: st
        d = fb%drive; l = fb%load
        basis%source = drive%source
        if (case_ntor > 0) then
            allocate (basis%coils)
            call windows(basis%coils)
        end if
        call ideal_fb_circuits_build(fb, basis, circuits, st, background=filament)
        if (st /= 0) error stop 'circuit traces'
        ! Same 1e-8 closure as the sampler for the axisymmetric circuits; a
        ! window entry holds the n = 3 P2 curl-free indicator (printed).
        if (maxval(circuits%closure(:size(drive%current) + 1)) > 1e-8_dp) &
            error stop 'circuit trace closure'
        currents = [drive%current, sqrt(1 - 4*delta)]
        if (case_ntor > 0) currents = [currents, 0.0_dp]
        call ideal_fb_set_currents(fb, circuits, currents, st)
        if (st /= 0) error stop 'circuit currents'
        err = max(maxval(abs(fb%drive - d))/maxval(abs(d)), &
            maxval(abs(fb%load - l))/maxval(abs(l)), &
            abs(fb%f_vac - sqrt(1 - 4*delta))/sqrt(1 - 4*delta))
        print '(a,es11.3)', 'circuit traces vs sampler traces (relative) ', err
        if (err > 1e-12_dp) error stop 'circuit traces differ from sampler traces'
        if (case_ntor > 0) then
            currents(size(currents)) = window_current
            call ideal_fb_set_currents(fb, circuits, currents, st)
            if (st /= 0) error stop 'window current'
            print '(a,es11.3,a,es11.3)', 'window current ', window_current, &
                '  closure ', circuits%closure(size(currents))
        end if
    end subroutine

    subroutine helical_case()
        ! Nonaxisymmetric synthetic case: Landreman drive plus a 3-fold
        ! window-coil circuit, nfp = 3, N = 1, chain 16. The window current
        ! gives |B_n=3| = 1e-4 on the outer midplane of Gamma_c. The solve
        ! must meet the same absolute residual and budgets; the n = 3
        ! boundary response must be linear in the current (I vs I/2).
        type(coil_circuit_basis_t) :: basis
        real(dp) :: bcos(3, 0:1, 1, 1), bsin(3, 0:1, 1, 1), unit, amp(2), d, f, pe
        integer :: st
        call windows(basis)
        call basis%sample(reshape([1.149_dp, 0.0_dp], [2, 1]), 3, 1, bcos, bsin, st)
        if (st /= 0) error stop 'window sample'
        unit = norm2([bcos(:, 1, 1, 1), bsin(:, 1, 1, 1)])
        case_nfp = 3; case_ntor = 1
        window_current = 1e-4_dp/unit
        call solve_case(16, d, f, pe)
        amp(1) = helical_amplitude
        window_current = 0.5_dp*window_current
        call solve_case(16, d, f, pe)
        amp(2) = helical_amplitude
        print '(a,2es12.4,a,f9.5)', 'n=3 boundary amplitude (I, I/2) ', amp, &
            '  ratio ', amp(1)/amp(2)
        if (.not. amp(2) > 1e-9_dp) error stop 'no n = 3 response'
        if (abs(amp(1)/amp(2) - 2) > 0.05_dp) error stop 'n = 3 response not linear'
        print '(a)', 'helical circuit case converged; numerical evidence only'
    end subroutine

    subroutine solve_case(nchain, distance, field_error, pressure_error)
        integer, intent(in) :: nchain
        real(dp), intent(out) :: distance, field_error, pressure_error
        type(gs_state_t) :: mesh
        type(cylindrical_chart_t) :: chart
        type(landreman_frozen_t) :: frozen
        type(ideal_fb_t) :: fb
        type(ideal_fb_report_t) :: report
        real(dp), allocatable :: chain(:, :), contour(:, :)
        real(dp) :: dr, dz, h, closure, b(3), j(3), psi, p, gp(2), q(3), qd(3, 3), det
        real(dp) :: mapped(3), r, z, projection(2), xi(2), point(3), vac(3)
        real(dp) :: vertex(2, 3), tau, vac_error, traction, edge_pressure
        integer :: a, c, k
        integer :: i, e, g, st, nc
        vertex = reshape([0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2, 3])
        nc = 3*nchain/2
        allocate (chain(2, nchain), contour(2, nc))
        do i = 1, nchain
            call landreman_surface(2*pi*(i - 1)/nchain, delta, &
                                   chain(1, i), chain(2, i), dr, dz)
        end do
        do i = 1, nc
            call gamma_c(real(i - 1, dp)/nc, contour(:, i))
        end do
        h = 2*pi*1.6_dp*sqrt(delta)/nc
        call initialize_gs(mesh, contour, 0.5_dp*h*h, st, 3, curve=gamma_c, &
                           interior_chain=chain)
        if (st /= 0) error stop 'Landreman mesh'
        call ideal_fb_init(fb, mesh, chain, chart, frozen, case_nfp, case_ntor, &
                           sqrt(1 - 4*delta), drive%circulation, st, mu0=1.0_dp)
        if (st /= 0) error stop 'Landreman init'
        call ideal_fb_set_drive(fb, sampler, st, closure)
        if (st /= 0 .or. closure > 1e-8_dp) error stop 'drive trace closure'
        if (use_circuits) call circuit_drive(fb)
        ! Family 1 has closed lines: target-start stationarity/refinement
        ! test only. No artificial residual shift or manufactured load.
        call ideal_fb_solve(fb, 1e-9_dp, 50, report, st, verbose=.false., &
                            residual_tolerance=1e-8_dp)
        print '(a,i0)', 'KKT unknowns ', fb%n
        print '(a,i0,a,i0,a,i0,a,es12.4)', 'solve chain=', nchain, ' status=', st, &
            ' iterations=', report%iterations, ' residual=', report%residual
        if (st /= 0) error stop 'Landreman solve did not meet absolute residual'
        helical_amplitude = 0
        do i = 1, fb%plasma%nodes
            if (.not. fb%plasma%boundary(i)) cycle
            if (fb%nmode > 1) helical_amplitude = max(helical_amplitude, &
                norm2(fb%plasma%x(:2, 2:3, i)))
        end do
        distance = 0; field_error = 0; pressure_error = 0
        do i = 1, fb%plasma%nodes
            if (.not. fb%plasma%boundary(i)) cycle
            call exact(fb%plasma%x(1, 1, i), fb%plasma%x(2, 1, i), b, j, psi, p, gp)
            call landreman_project_surface(fb%plasma%x(1, 1, i), &
                              fb%plasma%x(2, 1, i), delta, projection(1), projection(2))
            distance = max(distance, norm2(fb%plasma%x(:2, 1, i) - projection))
        end do
        do e = 1, fb%plasma%cells
            do g = 1, fb%plasma%nq
                q = matmul(fb%plasma%x(:, 1, fb%plasma%conn(:, e)), &
                           fb%plasma%shape(:, g))
                qd = 0
                qd(:, :2) = matmul(fb%plasma%x(:, 1, fb%plasma%conn(:, e)), &
                                   transpose(fb%plasma%dshape(:, :, g)))
                r = q(1); z = q(2)
                det = qd(1, 1)*qd(2, 2) - qd(1, 2)*qd(2, 1)
                mapped(1) = dot_product(qd(1, :2), fb%plasma%b(:2, g, e))/(r*det)
                mapped(3) = dot_product(qd(2, :2), fb%plasma%b(:2, g, e))/(r*det)
                mapped(2) = fb%plasma%b(3, g, e)/det
                call exact(r, z, b, j, psi, p, gp)
                field_error = max(field_error, norm2(mapped - b))
                pressure_error = max(pressure_error, abs(fb%plasma%p(g, e) - p))
            end do
        end do
        ! Trace probes include edge interiors; nodal agreement alone misses
        ! curved-interface interpolation error. Compare actual solved vacuum.
        vac_error = 0; traction = 0; edge_pressure = 0
        do e = 1, fb%ncell
            do a = 1, 3
                c = mod(a, 3) + 1
                if (fb%role(fb%conn(a, e)) /= 1) cycle
                if (fb%role(fb%conn(c, e)) /= 1) cycle
                do k = 0, 4
                    tau = real(k, dp)/4
                    xi = (1 - tau)*vertex(:, a) + tau*vertex(:, c)
                    call ideal_fb_shell_field(fb, e, xi, 0.0_dp, point, vac, st)
                    if (st /= 0) error stop 'solved vacuum trace'
                    call exact(point(1), point(2), b, j, psi, p, gp)
                    call landreman_project_surface(point(1), point(2), delta, &
                                                   projection(1), projection(2))
                    distance = max(distance, norm2(point(:2) - projection))
                    vac_error = max(vac_error, norm2(vac - b))
                    traction = max(traction, abs(p + (sum(b*b) - sum(vac*vac))/2))
                    edge_pressure = max(edge_pressure, abs(p))
                end do
            end do
        end do
        print '(a,3es12.4)', 'solved vacuum trace B, traction, edge p defects ', &
            vac_error, traction, edge_pressure
        ! Stress tolerance follows |B.dB|+|dB|^2/2+|dp| with B scale 1.2,
        ! field allowance 1% of 2a and pressure allowance 1% of 2delta.
        if (vac_error > 0.01_dp*2*sqrt(delta)) error stop 'solved vacuum field'
        if (traction > 1.2_dp*0.002_dp + 0.002_dp**2/2 + 0.0002_dp) &
            error stop 'solved traction budget'
        print '(a,3es12.4)', 'S0 distance, max B error, max pressure error ', &
            distance, field_error, pressure_error
        if (.not. all(ieee_is_finite([distance, field_error, pressure_error]))) &
            error stop 'nonfinite result'
        call psi_crosscheck(mesh, chain)

    end subroutine
    subroutine psi_crosscheck(state, chain)
        ! Same external ring sources through gs_free's ordinary coil input.
        ! GS uses chi=-psi; its SI mu0 is absorbed in p and circuit currents.
        type(gs_state_t), intent(inout) :: state
        real(dp), intent(in) :: chain(:, :)
        type(gs_circuit_set_t) :: coils
        type(gs_profile_law_t) :: law
        type(gs_fb_config_t) :: cfg
        type(gs_fb_system_t) :: sys
        type(gs_fb_eval_t) :: ev
        type(gs_fb_history_t) :: history
        real(dp) :: b(3), j(3), psi, p, gp(2), chi, grad(2), hess(2, 2), xi(2)
        real(dp) :: dist, berr, r, z, point(2), p_error
        integer :: k, st, e, i, ns
        logical :: found
        ns = size(drive%current)
        allocate (coils%coils(ns), coils%incidence(ns, 1))
        coils%incidence(:, 1) = drive%current
        do k = 1, ns
            coils%coils(k)%kind = gs_coil_filaments
            coils%coils(k)%fil_r = [drive%source(1, k)]
            coils%coils(k)%fil_z = [drive%source(2, k)]
            coils%coils(k)%fil_turns = [1.0_dp]
        end do
        allocate (law%knots(2), law%pp(4, 1), law%ff(4, 1))
        law%knots = [0.0_dp, 1.0_dp]; law%normalized = .true.; law%axis = delta
        law%p_edge = 0; law%f_edge = sqrt(1 - 4*delta)
        law%pp = 0; law%ff = 0; law%pp(1, 1) = 2/gs_mu0; law%ff(1, 1) = 2
        cfg%mode = GS_FB_TARGET_MEAN; cfg%points = chain
        cfg%axis_sign = 1; cfg%axis_guess = [1.0_dp, 0.0_dp]
        cfg%search_radius = 0.05_dp; cfg%tolerance = 1e-9_dp
        cfg%quadrature_tolerance = 1e-10_dp; cfg%jacobian_tolerance = 1e-10_dp
        cfg%max_iterations = 30
        call gs_fb_setup(state, sys, law, cfg, coils, st)
        if (st /= 0) error stop 'GS Landreman setup'
        do i = 1, size(state%points, 2)
            call exact(state%points(1, i), state%points(2, i), b, j, psi, p, gp)
            state%psi(i) = -psi
        end do
        call gs_fb_newton(state, sys, [1/gs_mu0], st, history)
        print '(a,i0,a,i0)', 'GS cross-check status=', st, &
            ' iterations=', history%iterations
        if (st /= 0) error stop 'GS Landreman solve'
        call gs_fb_evaluate(state, sys, state%psi, [1/gs_mu0], ev, st, jacobian=.false.)
        if (st /= 0) error stop 'GS Landreman evaluate'
        dist = 0; berr = 0; p_error = 0
        do i = 1, size(chain, 2)
            r = chain(1, i); z = chain(2, i)
            call gs_fb_locate(state, r, z, e, xi, found)
            if (.not. found) error stop 'GS surface locate'
            call gs_element_jet(state, e, r, z, chi, grad, hess, st)
            if (st /= 0) error stop 'GS surface jet'
            call exact(r, z, b, j, psi, p, gp)
            dist = max(dist, abs(chi - ev%psi_bound)/norm2(gp))
        end do
        do e = 1, size(state%triangles, 2)
            point = sum(state%points(:, state%triangles(:3, e)), dim=2)/3
            r = point(1); z = point(2)
            call exact(r, z, b, j, psi, p, gp)
            if (psi > 0.9_dp*delta) cycle
            call gs_element_jet(state, e, r, z, chi, grad, hess, st)
            if (st /= 0) error stop 'GS interior jet'
            berr = max(berr, norm2([-grad(2)/r, grad(1)/r] - b([1, 3])))
            p_error = max(p_error, abs(2*(chi - ev%psi_bound) - p))
        end do
        print '(a,3es12.4)', 'GS S0 distance indicator, poloidal B, p errors ', &
            dist, berr, p_error
        if (dist > 0.01_dp*sqrt(delta)) error stop 'GS S0 return'
        if (berr > 0.01_dp*2*sqrt(delta)) error stop 'GS field oracle'
        if (p_error > 0.01_dp*2*delta) error stop 'GS pressure oracle'
    end subroutine
end program test_landreman_free_boundary
