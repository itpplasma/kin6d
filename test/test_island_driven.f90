module island_driven_case
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_island_leaves, only: driven_cylinder_data
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp), r0 = 20
    real(dp) :: q0 = 1.2_dp, qc = 1.8_dp
    real(dp) :: amplitude = 0
contains
    function zero_source(x) result(s)
        real(dp), intent(in) :: x(2)
        real(dp) :: s
        s = 0
    end function
    subroutine mapping(u, theta, x)
        real(dp), intent(in) :: u, theta
        real(dp), intent(out) :: x(2)
        x = [r0 + u*cos(theta), u*sin(theta)]
    end subroutine
    function conductivity(x) result(w)
        real(dp), intent(in) :: x(2)
        real(dp) :: w, q, psi
        call driven_cylinder_data(norm2(x - [r0, 0.0_dp]), q0, qc, q, w, psi)
    end function
    subroutine boundary_potential(x, a)
        real(dp), intent(in) :: x(2)
        complex(dp), intent(out) :: a(3)
        complex(dp) :: z
        ! a_1=eps (x-iZ)^2/2 in R A_phi; gauge Lambda=i a_1
        ! removes A_phi and preserves the imposed normal trace of curl A.
        z = cmplx(x(1) - r0, -x(2), dp)
        a = amplitude*[cmplx(0.0_dp, 1.0_dp, dp)*z, &
                       cmplx(0.0_dp, 0.0_dp, dp), z]
    end subroutine
end module

program test_island_driven
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use island_driven_case, only: pi, r0, q0, qc, amplitude, mapping, &
                                  conductivity, boundary_potential, zero_source
    use island_outer_reference, only: outer_response
    use kin6d_island_leaves, only: driven_cylinder_data
    use kin6d_mhd_mesh, only: flux_mesh_t, flux_mesh_build
    use kin6d_equilibrium_nonnested, only: free_b_field_t, nn_build, nn_report_t, &
                                           nn_vacuum_mode
    use fortfem_meridian_feec_n, only: meridian_lagrange_nodes
    use kin6d_nonnested_coupled, only: nc_problem_t, nc_layout_t, nc_report_t, &
                                       nc_layout, nc_continue
    use kin6d_island_diagnostics, only: island_trace_options_t, island_point_t, &
                       island_width_t, island_trace, island_periodic_point, island_width
    use kin6d_equilibrium_state, only: EQ_VALUES, EQ_BR, EQ_BZ, EQ_BPHI
    use kin6d_nonnested_residual, only: nn_residual
    use kin6d_equilibrium_residual, only: equilibrium_residual_t
    implicit none
    type(free_b_field_t) :: field
    type(flux_mesh_t) :: mesh
    type(nc_problem_t) :: problem
    type(nc_layout_t) :: lay
    type(nc_report_t) :: report
    type(nn_report_t) :: vacuum_report
    type(island_trace_options_t) :: opt
    type(island_point_t) :: op, xp
    type(island_width_t) :: width
    type(equilibrium_residual_t) :: stress
    real(dp), allocatable :: rings(:), nodes(:, :), x(:), seed(:)
    complex(dp), allocatable :: avec(:)
    real(dp) :: response(8), tearing, drive, delta, theory, rs, shear, wq, psi, q
    real(dp) :: eps(3), widths(3), exponent, phase, centre(2), direction(2), inside
    real(dp) :: points(3, 2), values(EQ_VALUES, 2), qmeasured(2), r
    integer :: k, nr, modes, j, st, statuses(2), lev
    real(dp) :: gperp = 1.0e-4_dp
    logical :: stable_case, prefactor_failed
    character(32) :: arg, mode
    call get_command_argument(1, mode)
    stable_case = index(trim(mode), 'stable') == 1
    if (stable_case) then
        q0 = 0.5_dp; qc = 2.5_dp; gperp = 1.0e-8_dp
    end if
    call get_command_argument(5, arg)
    if (len_trim(arg) > 0) read (arg, *) q0
    call get_command_argument(6, arg)
    if (len_trim(arg) > 0) read (arg, *) qc
    call get_command_argument(7, arg)
    if (len_trim(arg) > 0) read (arg, *) gperp
    do lev = 1, 8
        delta = 1.0e-4_dp/2**(lev - 1)
        call outer_response(q0, q0 + qc, delta, response(lev), tearing, drive, st)
        if (st /= 0) error stop 'outer ODE'
        print '(a,4es14.6)', 'cut / response / Delta-prime / drive ', &
            delta, response(lev), tearing, drive
    end do
    if (abs(response(8) - response(7)) > 0.01_dp*abs(response(8))) &
        error stop 'unresolved outer Frobenius limit'
    if (stable_case .and. tearing >= 0) error stop 'profile is not tearing stable'
    call outer_response(0.0_dp, 3.0_dp, delta, q, tearing, drive, st)
    print '(a,i0,a,3es14.6)', 'vacuum reduction status ', st, &
        ' response / Delta-prime / drive ', q, tearing, drive
    if (st /= 0 .or. abs(q - 2.0_dp/3) > 1.0e-6_dp) &
        error stop 'outer vacuum r^m oracle'
    if (trim(mode) == 'outer' .or. trim(mode) == 'stable_outer') then
        call check_radial_reference()
        call check_signed_projection()
        print '(a)', 'PASS: cylindrical driven outer response and vacuum reduction'
        stop
    end if
    call get_command_argument(2, arg); read (arg, *) k
    call get_command_argument(3, arg); read (arg, *) nr
    call get_command_argument(4, arg); read (arg, *) modes
    if (k < 2 .or. k > 3 .or. nr < 4 .or. modes < 1) error stop 'arguments'
    centre = [r0, 0.0_dp]; rs = sqrt((2 - q0)/qc); shear = qc*rs*rs
    allocate (rings(nr))
    rings = [(real(j, dp)/nr, j=1, nr)]
    call flux_mesh_build(rings, 8, 8*nr, k, mapping, mesh, st)
    if (st /= 0) error stop 'mesh'
    call nn_build(mesh, k, modes, centre, 1.0_dp, r0, 1.0_dp, field, st)
    if (st /= 0) error stop 'build'
    call meridian_lagrange_nodes(field%fs, nodes)
    lay = nc_layout(field)
    allocate (x(lay%n)); x = 0
    do j = 1, lay%n0
        call driven_cylinder_data(norm2(nodes(:, j) - centre), q0, qc, q, wq, psi)
        x(lay%ia + field%fs%nn + j) = psi
    end do
    x(lay%iv) = -4*pi/q0
    problem%current = -2*pi/(r0*(q0 + qc))
    problem%g_par_shape => conductivity
    problem%source => zero_source
    problem%g_par = 1; problem%g_perp = gperp
    problem%chi_par = 1; problem%chi_perp = 1
    call nc_continue(field, problem, x, report, st)
    print '(a,i0,a,i0,a,i0,a,i0,a,es12.4)', 'P', k, ' nr=', nr, ' N=', modes, &
        ' base status ', st, ' voltage ', report%voltage
    print '(a,2i4,a,41es11.3)', 'base stages/steps ', report%stages, &
        report%total_steps, ' history ', report%history(:report%newton_steps)
    print '(a,5es12.4)', 'base block backward errors ', report%residual
    if (st /= 0) error stop 'axisymmetric Ohmic seed'
    call nn_residual(field, modes, stress, st)
    if (st /= 0) error stop 'base force residual'
    print '(a,es13.5)', 'base D3 ', stress%normalized
    points(:, 1) = [r0 + 0.10_dp, 0.0_dp, 0.0_dp]
    points(:, 2) = [r0 + 0.95_dp, 0.0_dp, 0.0_dp]
    call field%sample(points, values, statuses)
    if (any(statuses /= 0)) error stop 'q sampling'
    do j = 1, 2
        r = points(1, j) - r0
        qmeasured(j) = r*values(EQ_BPHI, j)/(r0*values(EQ_BZ, j))
    end do
    print '(a,2es14.6)', 'q at r=0.10,0.95 ', qmeasured
    if (abs(qmeasured(1) - q0) > 0.2_dp .or. &
        abs(qmeasured(2) - (q0 + qc*0.95_dp**2)) > 0.2_dp) &
        error stop 'Ohmic q profile'
    if (maxval(abs(x(lay%ip + 1:lay%ip + lay%nv))) > 0) &
        error stop 'zero-pressure transport invariant lost'
    if (trim(mode) == 'base') then
        print '(a)', 'PASS: zero-beta Ohmic solve preserves homogeneous pressure'
        stop
    end if
    seed = x
    eps = [0.002_dp, 0.0005_dp, 0.000125_dp]
    if (stable_case) eps = eps/4
    prefactor_failed = .false.
    opt%root_tol = 2.0e-6_dp; opt%rtol = 1.0e-8_dp
    do lev = 1, 3
        amplitude = eps(lev); x = seed
        call nn_vacuum_mode(field, 1, boundary_potential, vacuum_report, st, &
                            solution=avec)
        if (st /= 0) error stop 'boundary drive'
        x(lay%ia + lay%na1 + 1:lay%ia + lay%na1 + field%fs%nn) = &
            2*real(avec, dp)
        x(lay%ia + 2*lay%na1 + 1:lay%ia + 2*lay%na1 + field%fs%nn) = &
            -2*aimag(avec)
        call nc_continue(field, problem, x, report, st)
        print '(a,es13.5,a,i0)', 'eps ', amplitude, ' solve status ', st
        if (st /= 0) error stop 'driven Ohmic solve'
        call nn_residual(field, modes, stress, st)
        if (st /= 0) error stop 'driven force residual'
        print '(a,es13.5)', 'driven D3 ', stress%normalized
        call response_samples()
        if (trim(mode) == 'response' .or. trim(mode) == 'stable_response') stop
        theory = 4*sqrt(abs(amplitude*response(8)/r0)/(rs/(2*r0)*shear/rs))
        phase = 0
        if (response(8) < 0) phase = pi/2
        call island_periodic_point(field, centre + rs*[cos(phase), sin(phase)], &
                                   0.0_dp, 2, centre, 2, 1, opt, op)
        if (stable_case) then
            ! The stable response puts X off Z=0. No Z=0 symmetry bracket applies.
            call island_periodic_point(field, &
                centre + rs*[cos(phase + pi/2), sin(phase + pi/2)], &
                0.0_dp, 2, centre, 2, 1, opt, xp)
        else
            call island_periodic_point(field, &
                centre + rs*[cos(phase + pi/2), sin(phase + pi/2)], &
                0.0_dp, 2, centre, 2, 1, opt, xp, [r0 + 0.4_dp, r0 + 0.9_dp])
        end if
        print '(a,2i3,a,4es13.5)', 'O/X status ', op%status, xp%status, &
            ' closure/residue ', op%closure, xp%closure, op%residue, xp%residue
        if (op%status /= 0 .or. xp%status /= 0) error stop 'driven O/X unresolved'
        if (op%residue <= 0 .or. xp%residue >= 0) error stop 'driven O/X type'
        direction = (op%x - centre)/norm2(op%x - centre)
        inside = norm2(op%x - centre)
        ! Stronger fixed-step ladders are explicit in the error ledger.
        opt%fixed_steps = 512
        call island_width(field, [max(0.1_dp, inside - 1.5_dp*theory), &
                     min(0.98_dp, inside + 1.5_dp*theory)], inside, centre, direction, &
                          0.0_dp, 64.0_dp, 2, 1, opt, minor_radius, 'r [m]', width)
        opt%fixed_steps = 0
        if (width%status /= 0) error stop 'driven width unresolved'
        widths(lev) = width%metres
        print '(a,6es14.6)', 'eps / W / theory / trace / duration / phase ', &
            amplitude, widths(lev), theory, width%trace_error, &
            width%duration_error, op%phase
        if (abs(atan2(sin(op%phase - 2*phase), cos(op%phase - 2*phase))) &
            > 0.1_dp + op%phase_error) error stop 'phase not locked to drive'
        ! Mesh/mode uncertainty is NOT supplied by a single run. This
        ! provisional necessary gate includes only tracing and O(a/R0).
        if (abs(widths(lev) - theory) > width%trace_error + width%duration_error &
            + width%scan_error + theory/r0) prefactor_failed = .true.
    end do
    exponent = log(widths(1)/widths(3))/log(eps(1)/eps(3))
    print '(a,f10.5)', 'NN2 driven width exponent ', exponent
    if (abs(exponent - 0.5_dp) > 0.05_dp) error stop 'NN2 driven scaling'
    if (prefactor_failed) error stop 'linear prefactor mismatch'
    print '(a)', 'PASS: necessary driven gates; mesh/mode ladders still required'
contains
    subroutine response_grid(rad, xyz)
        real(dp), intent(in) :: rad
        real(dp), intent(out) :: xyz(3, 128)
        real(dp) :: theta, phi
        integer :: it, ip, idx
        idx = 0
        do ip = 0, 7
            phi = 2*pi*ip/8
            do it = 0, 15
                theta = 2*pi*it/16
                idx = idx + 1
                xyz(:, idx) = [r0 + rad*cos(theta), rad*sin(theta), phi]
            end do
        end do
    end subroutine

    real(dp) function flux_projection(rad, xyz, bv, eps) result(f)
        real(dp), intent(in) :: rad, xyz(3, 128), bv(EQ_VALUES, 128), eps
        real(dp) :: theta, phi, br
        integer :: it, ip, idx
        f = 0; idx = 0
        do ip = 0, 7
            phi = 2*pi*ip/8
            do it = 0, 15
                theta = 2*pi*it/16
                idx = idx + 1
                br = bv(EQ_BR, idx)*cos(theta) + bv(EQ_BZ, idx)*sin(theta)
                f = f + xyz(1, idx)*rad*br*sin(2*theta - phi)
            end do
        end do
        f = f/(128*eps)
    end function

    subroutine check_signed_projection()
        ! Independent analytic curl of a_1=eps(i z,0,z), z=R-R0-i Z.
        ! This is a signed extraction oracle, not a native vacuum-solve claim.
        real(dp), parameter :: eps0 = 0.001_dp, rr(3) = [0.25_dp, 2.0_dp/3, 1.0_dp]
        real(dp) :: xyz(3, 128), bv(EQ_VALUES, 128), xr, zz, phi, measured
        integer :: ir, j, polarity
        do ir = 1, size(rr)
            call response_grid(rr(ir), xyz)
            do polarity = -1, 1, 2
                bv = 0
                do j = 1, 128
                    xr = xyz(1, j) - r0; zz = xyz(2, j); phi = xyz(3, j)
                    bv(EQ_BR, j) = polarity*2*eps0 &
                        *(zz*cos(phi) - xr*sin(phi))/xyz(1, j)
                    bv(EQ_BZ, j) = polarity*2*eps0 &
                        *(xr*cos(phi) + zz*sin(phi))/xyz(1, j)
                end do
                measured = flux_projection(rr(ir), xyz, bv, eps0)
                if (abs(measured - polarity*rr(ir)**2) > 2.0e-13_dp) &
                    error stop 'signed analytic radial-flux projection'
            end do
        end do
        print '(a)', 'PASS: signed analytic curl projection +/-r^2; edge +/-1'
    end subroutine

    subroutine check_radial_reference()
        real(dp) :: rt(5), radial(5), ref, dt, drv, resonant
        integer :: stat
        resonant = sqrt((2 - q0)/qc)
        rt = [0.0_dp, resonant/2, resonant, (1 + resonant)/2, 1.0_dp]
        call outer_response(q0, q0 + qc, delta, ref, dt, drv, stat, rt, radial)
        if (stat /= 0) error stop 'radial outer reference'
        if (abs(radial(5) - 1) > 1.0e-12_dp) error stop 'outer edge normalization'
        if (abs(radial(3) - response(8)) > 1.0e-12_dp) &
            error stop 'outer resonant normalization'
        ! With J'=0, the independent regular solution is f(r)/f_b=r^2.
        rt = [0.0_dp, 0.25_dp, sqrt(2.0_dp/3), 0.99_dp, 1.0_dp]
        call outer_response(0.0_dp, 3.0_dp, delta, ref, dt, drv, stat, rt, radial)
        if (stat /= 0) error stop 'radial vacuum reference'
        if (maxval(abs(radial - rt**2)) > 1.0e-6_dp) &
            error stop 'outer vacuum radial r^2 oracle'
        rt(1) = resonant + delta/2
        call outer_response(q0, q0 + qc, delta, ref, dt, drv, stat, rt, radial)
        if (stat /= 4) error stop 'unresolved Frobenius cut must be rejected'
        print '(a)', 'PASS: signed outer radial r^2/edge/resonance; ' // &
            'unresolved cut rejected'
    end subroutine

    subroutine response_samples()
        real(dp) :: xyz(3, 128), bv(EQ_VALUES, 128), radii(5), radial(5), rad, f
        real(dp) :: ref, dt, drv
        integer :: stat(128), ir, status
        radii = [(rs + (0.99_dp - rs)*ir/4, ir=0, 4)]
        call outer_response(q0, q0 + qc, delta, ref, dt, drv, status, radii, radial)
        if (status /= 0) error stop 'response radial outer reference'
        print '(a,2es16.7)', 'resonant radius / outer flux response ', rs, ref
        do ir = 0, 4
            rad = radii(ir + 1)
            call response_grid(rad, xyz)
            call field%sample(xyz, bv, stat)
            if (any(stat /= 0)) error stop 'response sampling'
            f = flux_projection(rad, xyz, bv, amplitude)
            print '(a,3es16.7)', 'radius / measured flux response / radial outer ', &
                rad, f, radial(ir + 1)
        end do
    end subroutine
    function minor_radius(x) result(u)
        real(dp), intent(in) :: x(2)
        real(dp) :: u
        u = norm2(x - centre)
    end function
end program
