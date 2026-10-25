module nn_exterior_operator_data
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
contains
    function source(x) result(value)
        real(dp), intent(in) :: x(2)
        real(dp) :: value
        value = 0*x(1)
    end function
end module

program test_nonnested_exterior_operators
    ! Independent unit-sphere l=1 harmonic oracle for the frozen composition.
    ! The supplied normal-map columns are synthetic FEEC-to-trace inputs.
    ! This does not qualify their geometric construction or full BEM assembly.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use nn_solovev_fixture, only: r0, minor, ftor, mapping
    use nn_exterior_operator_data, only: source
    use kin6d_mhd_mesh, only: flux_mesh_t, flux_mesh_build
    use kin6d_equilibrium_nonnested, only: free_b_field_t, nn_build
    use kin6d_nonnested_coupled, only: nc_problem_t, nc_layout_t, nc_layout
    use kin6d_nonnested_exterior, only: nn_exterior_t, nn_attach_exterior_operators
    use fortsparse, only: csc_matvec
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp), mass = 4*pi/3
    type(flux_mesh_t) :: mesh
    type(free_b_field_t) :: field
    type(nc_problem_t) :: problem, invalid
    type(nc_layout_t) :: lay
    type(nn_exterior_t) :: exterior, rejected
    real(dp), allocatable :: normal(:, :), volume(:), a(:), direction(:), gradient(:)
    real(dp) :: s(3, 3), coil(3), f(3, 2), av(2), v(2), g(3), c(3), dc(3)
    real(dp) :: u(3), actual_energy, expected_energy, expected_gradient(2)
    real(dp) :: axis(3), response(3), bn, boundary_energy, angular_energy
    real(dp) :: full_energy, net_flux, period, energy_plus, energy_minus, h
    real(dp) :: tangent(3), theta, dtheta, point(3), sign_error, pairing_error
    integer :: status, i, j, first, last, interior, ids(2)

    call flux_mesh_build([1.0_dp], 4, 128, 2, mapping, mesh, status)
    call require(status == 0, 'mesh')
    call nn_build(mesh, 2, 1, [r0, 0.0_dp], minor, ftor, 1.0_dp, field, status)
    call require(status == 0, 'FEEC field')
    lay = nc_layout(field)
    problem%source => source
    first = findloc(field%fs%sp(1)%bnd, .true., 1)
    last = findloc(field%fs%sp(1)%bnd, .true., 1, back=.true.)
    call require(first > 0 .and. last > first, 'distinct boundary columns')
    ids = [first, last]
    allocate (normal(3, lay%na), volume(lay%na), a(lay%na), direction(lay%na))
    ! Trace basis n_x,n_y,n_z on the unit sphere, not orthonormalized.
    ! The exterior u=(c.x)/r**3 has B_response=-grad(u) and weak DtN S=2 M.
    s = 0
    do j = 1, 3
        s(j, j) = 2*mass
    end do
    f(:, 1) = [1.0_dp, 0.3_dp, -0.2_dp]
    f(:, 2) = [-0.4_dp, 0.2_dp, 0.7_dp]
    normal = 0; normal(:, ids) = mass*f
    g = [0.2_dp, -0.3_dp, 0.4_dp]; coil = mass*g
    volume = 0
    call nn_attach_exterior_operators(field, problem, s, normal, coil, volume, &
        exterior, status)
    call require(status == 0, 'frozen exterior attachment')
    call require(.not. allocated(exterior%flux_pairing), 'no fabricated BEM trace')

    av = [0.7_dp, -0.2_dp]; v = [-0.3_dp, 0.6_dp]
    a = 0; a(ids) = av
    direction = 0; direction(ids) = v
    c = (matmul(f, av) - g)/2
    call require(norm2(c) > 0.1_dp, 'nonzero analytic response')
    u = matmul(exterior%potential_map(:, :size(exterior%magnetic_dofs)), &
        a(exterior%magnetic_dofs)) &
        - exterior%potential_map(:, size(exterior%magnetic_dofs) + 1)
    sign_error = norm2(u - c)/norm2(c)
    call require(sign_error <= 1.0e-12_dp, 'response sign and trace coefficients')
    call require(norm2(matmul(s, u) - (matmul(normal, a) - coil)) &
        <= 1.0e-12_dp*norm2(coil), 'weak normal matching')

    ! Six axes integrate every degree-two sphere polynomial exactly. Integrate
    ! the analytic exterior field independently of the supplied Schur matrix.
    boundary_energy = 0; angular_energy = 0; net_flux = 0
    do j = 1, 3
        do i = -1, 1, 2
            axis = 0; axis(j) = real(i, dp)
            response = 3*dot_product(c, axis)*axis - c
            bn = dot_product(response, axis)
            call require(abs(bn - dot_product(matmul(f, av) - g, axis)) &
                < 1.0e-12_dp, 'physical normal response')
            boundary_energy = boundary_energy + (4*pi/6)*dot_product(c, axis)*bn/2
            angular_energy = angular_energy + (4*pi/6)*dot_product(response, response)
            net_flux = net_flux + (4*pi/6)*bn
        end do
    end do
    ! B scales as r**(-3), so integral_1^infinity r**2 r**(-6) dr = 1/3.
    full_energy = angular_energy/(6*field%mu0)
    expected_energy = full_energy/(2*pi)
    actual_energy = (dot_product(a, csc_matvec(problem%exterior, a))/2 &
        - dot_product(a, problem%coil_drive) + dot_product(coil, g)/(8*pi)) &
        /field%mu0
    call require(expected_energy > 0, 'positive exterior magnetic energy')
    call require(abs(actual_energy - expected_energy) &
        <= 1.0e-12_dp*expected_energy, 'once-only toroidal-mean normalization')
    call require(abs(boundary_energy/field%mu0 - full_energy) &
        <= 1.0e-12_dp*full_energy, 'independent volume and boundary energies')
    call require(abs(net_flux) <= 1.0e-12_dp*norm2(c), 'zero net response flux')

    gradient = (csc_matvec(problem%exterior, a) - problem%coil_drive)/field%mu0
    expected_gradient = matmul(transpose(f), matmul(f, av) - g)/(3*field%mu0)
    call require(norm2(gradient(ids) - expected_gradient) &
        <= 1.0e-12_dp*norm2(expected_gradient), 'independent Schur gradient')
    h = 1.0e-5_dp
    dc = matmul(f, v)/2
    energy_plus = mass*sum((c + h*dc)**2)/(2*pi*field%mu0)
    energy_minus = mass*sum((c - h*dc)**2)/(2*pi*field%mu0)
    pairing_error = abs((energy_plus - energy_minus)/(2*h) &
        - dot_product(gradient, direction))/abs(dot_product(gradient, direction))
    call require(pairing_error <= 1.0e-9_dp, 'independent energy derivative')

    ! The response is a single-valued scalar gradient: integrate one closed
    ! equatorial period. The tangential field has the same -grad(u) sign.
    period = 0; dtheta = 2*pi/64
    do j = 1, 64
        theta = (j - 0.5_dp)*dtheta
        point = [cos(theta), sin(theta), 0.0_dp]
        tangent = [-sin(theta), cos(theta), 0.0_dp]
        response = 3*dot_product(c, point)*point - c
        period = period + dot_product(response, tangent)*dtheta
        call require(abs(dot_product(response, tangent) + dot_product(c, tangent)) &
            <= 1.0e-12_dp, 'tangential response sign')
    end do
    call require(abs(period) <= 1.0e-12_dp*norm2(c), 'zero gradient period')

    ! Counterexamples are genuinely inadmissible inputs, not empty responses.
    invalid%source => source
    invalid%current = 1
    call nn_attach_exterior_operators(field, invalid, s, normal, coil, volume, &
        rejected, status)
    call require(status /= 0, 'nonzero current needs harmonic response')
    invalid%current = 0
    call nn_attach_exterior_operators(field, invalid, s(:2, :2), normal, coil, volume, &
        rejected, status)
    call require(status /= 0, 'inconsistent trace basis size')
    interior = findloc(field%fs%sp(1)%bnd, .false., 1)
    call require(interior > 0, 'interior trace counterexample exists')
    normal(1, interior) = 0.1_dp
    call nn_attach_exterior_operators(field, invalid, s, normal, coil, volume, &
        rejected, status)
    call require(status /= 0, 'interior normal column rejected')
    normal(1, interior) = 0
    s(1, 2) = 0.1_dp
    call nn_attach_exterior_operators(field, invalid, s, normal, coil, volume, &
        rejected, status)
    call require(status /= 0, 'nonsymmetric energy map rejected')
    s(1, 2) = 0
    coil(1) = ieee_value(0.0_dp, ieee_quiet_nan)
    call nn_attach_exterior_operators(field, invalid, s, normal, coil, volume, &
        rejected, status)
    call require(status /= 0, 'nonfinite coil trace rejected')
    print '(a,4es16.7)', 'frozen sphere energy/sign/energy-JVP/netflux ', &
        actual_energy, sign_error, pairing_error, net_flux
    print '(a)', 'PASS: frozen exterior composition; geometric BEM and G8 unqualified'
contains
    subroutine require(ok, label)
        logical, intent(in) :: ok
        character(*), intent(in) :: label
        if (ok) return
        print '(2a,i0)', trim(label), ': status ', status
        error stop 'frozen exterior operator gate'
    end subroutine
end program
