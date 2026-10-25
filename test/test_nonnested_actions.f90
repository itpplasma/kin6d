module nn_action_data
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    implicit none
contains
    function zero_source(x) result(value)
        real(dp), intent(in) :: x(2)
        real(dp) :: value
        value = 0*x(1)
    end function

    function bad_source(x) result(value)
        real(dp), intent(in) :: x(2)
        real(dp) :: value
        value = ieee_value(x(1), ieee_quiet_nan)
    end function

    function negative_shape(x) result(value)
        real(dp), intent(in) :: x(2)
        real(dp) :: value
        value = -1 + 0*x(1)
    end function

    function conductivity(x, phi) result(value)
        real(dp), intent(in) :: x(2), phi
        real(dp) :: value
        value = 1 + 0.2_dp*cos(phi) + 0.01_dp*x(2)**2
    end function
end module

program test_nonnested_actions
    ! Operator checks at nonstationary states, including boundary columns.
    ! These forced/incompatible data do not claim to be an equilibrium.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use nn_solovev_fixture, only: r0, minor, ftor, pp, mapping, exact
    use nn_action_data, only: zero_source, bad_source, negative_shape, conductivity
    use kin6d_mhd_mesh, only: flux_mesh_t, flux_mesh_build
    use fortfem_meridian_feec_n, only: meridian_lagrange_nodes
    use fortsparse, only: csc_t, csc_matvec, csc_transpose, fortsparse_status_t, &
        status_ok
    use kin6d_equilibrium_nonnested, only: free_b_field_t, nn_build
    use kin6d_nonnested_coupled, only: nc_problem_t, nc_layout_t, nc_layout, &
        nc_fixed, nc_reduce, nc_residual, nc_jacobian, nc_jvp, nc_vjp
    use kin6d_nonnested_exterior, only: nn_exterior_t, nn_attach_exterior_operators
    use kin6d_island_closure, only: ISLAND_INVALID, ISLAND_UNSUPPORTED, &
        ISLAND_NULL_FIELD, ISLAND_BAD_DATA
    implicit none
    type(flux_mesh_t) :: mesh
    type(free_b_field_t) :: field
    type(nc_problem_t) :: problem
    type(nc_layout_t) :: lay
    type(csc_t) :: jac, transposed
    type(fortsparse_status_t) :: sparse_status
    real(dp), allocatable :: nodes(:, :), x(:), direction(:), weight(:)
    real(dp), allocatable :: action(:), adjoint(:), reference(:), raw_weight(:)
    real(dp), allocatable :: plus(:), minus(:)
    logical, allocatable :: fixed(:)
    real(dp) :: psi, b(2), voltage, error, pairing(2), h
    real(dp) :: fd_errors(3), fd_noise(3)
    character(24) :: mode, arg
    integer :: degree, i, first, status, level

    call get_command_argument(1, mode)
    call get_command_argument(2, arg)
    read (arg, *) degree
    if (degree < 2 .or. degree > 3) error stop 'degree'
    call flux_mesh_build([0.5_dp, 1.0_dp], 8, 1024, degree, mapping, mesh, status)
    call require(status == 0, 'mesh')
    call nn_build(mesh, degree, 2, [r0, 0.0_dp], minor, ftor, 1.0_dp, field, status)
    call require(status == 0, 'FEEC field')
    call meridian_lagrange_nodes(field%fs, nodes)
    lay = nc_layout(field)
    allocate (x(lay%n), direction(lay%n), weight(lay%n))
    x = 0
    do i = 1, lay%n0
        call exact(nodes(:, i), psi, b)
        x(field%fs%nn + i) = psi
        x(lay%ip + i) = pp*psi
    end do
    voltage = 2*acos(-1.0_dp)*r0*r0*pp
    x(lay%iv) = voltage
    direction = [(sin(1.3_dp*i), i=1, lay%n)]
    direction(:lay%na) = direction(:lay%na)*0.01_dp*maxval(abs(x(:lay%na)))
    direction(lay%il + 1:lay%ip) = direction(lay%il + 1:lay%ip)*1.0e-3_dp
    direction(lay%ip + 1:lay%iphi) = direction(lay%ip + 1:lay%iphi) &
        *0.03_dp*maxval(abs(x(lay%ip + 1:lay%iphi)))
    direction(lay%iphi + 1:) = direction(lay%iphi + 1:)*0.01_dp*abs(voltage)
    x = x + 0.2_dp*direction
    weight = [(cos(0.73_dp*i), i=1, lay%n)]
    problem%source => zero_source
    problem%g_par_shape_3d => conductivity
    problem%chi_par = 100; problem%g_perp = 0.01_dp
    ! An even FFT grid also checks that no odd-grid convention is assumed.
    problem%nphi = 8
    select case (trim(mode))
    case ('fixed')
    case ('exterior')
        call attach_exterior()
    case ('floating')
        allocate (problem%floating(lay%n0)); problem%floating = .false.
        first = 0
        do i = 1, lay%n0
            if (.not. field%fs%sp(0)%bnd(i)) cycle
            problem%floating(i) = .true.
            first = first + 1
            if (first == 3) exit
        end do
        call require(first == 3, 'three floating nodes')
    case ('voltage')
        problem%closure%drive = 'V'
        problem%voltage = x(lay%iv)
    case ('invalid')
        call invalid_inputs()
        print '(a)', 'PASS: I3 invalid input and residual domain'
        stop
    case default
        error stop 'mode'
    end select
    fixed = nc_fixed(field, problem)
    call nc_jacobian(field, problem, x, jac, status)
    call require(status == 0, 'assembled Jacobian')
    reference = csc_matvec(jac, direction)
    where (fixed) reference = 0
    call nc_reduce(field, problem, reference)
    call require(norm2(reference) > tiny(1.0_dp), 'nonzero derivative')
    call nc_jvp(field, problem, x, direction, action, status)
    call require(status == 0, 'matrix-free JVP')
    error = norm2(action - reference)/norm2(reference)
    print '(a,es12.4)', 'JVP / assembled relative error ', error
    call require(error <= 1.0e-10_dp, 'JVP / assembled')

    ! A coefficient direction has larger pointwise curl at higher order.
    ! Resolve central-difference truncation at the SAME state/direction; do
    ! not confuse O(h**2) oracle error with a constitutive derivative defect.
    do level = 1, 3
        h = 1.0e-4_dp/2**(level - 1)
        call nc_residual(field, problem, x + h*direction, plus, status)
        call require(status == 0, 'positive residual perturbation')
        call nc_residual(field, problem, x - h*direction, minus, status)
        call require(status == 0, 'negative residual perturbation')
        fd_errors(level) = norm2((plus - minus)/(2*h) - action)/norm2(action)
        fd_noise(level) = 32*epsilon(h)*(norm2(plus) + norm2(minus)) &
            /(2*h*norm2(action))
    end do
    print '(a,3es12.4)', 'JVP / central difference h,h/2,h/4 ', fd_errors

    call nc_vjp(field, problem, x, weight, adjoint, status)
    call require(status == 0, 'matrix-free VJP')
    pairing = [dot_product(weight, action), dot_product(adjoint, direction)]
    error = abs(pairing(1) - pairing(2))/max(tiny(1.0_dp), maxval(abs(pairing)))
    print '(a,es12.4)', 'Euclidean JVP / VJP pairing error ', error
    call require(error <= 1.0e-10_dp, 'adjoint pairing')
    ! Independently reverse the prescribed row map for the assembled oracle.
    raw_weight = weight
    if (allocated(problem%floating)) then
        first = findloc(problem%floating, .true., 1)
        do i = 1, lay%n0
            if (.not. problem%floating(i)) cycle
            raw_weight(lay%il + i) = weight(lay%il + first)
            raw_weight(lay%iphi + i) = weight(lay%iphi + first)
        end do
    end if
    where (fixed) raw_weight = 0
    call csc_transpose(jac, transposed, sparse_status)
    call require(status_ok(sparse_status), 'assembled transpose')
    reference = csc_matvec(transposed, raw_weight)
    error = norm2(adjoint - reference)/norm2(reference)
    print '(a,es12.4)', 'VJP / assembled transpose relative error ', error
    call require(error <= 1.0e-10_dp, 'VJP / assembled transpose')
    call require(fd_errors(3) <= 1.0e-7_dp, 'finite-difference JVP')
    call require(fd_errors(3) <= 0.12_dp*fd_errors(1) + maxval(fd_noise), &
        'central-difference truncation convergence')
    print '(a,i0,2a)', 'PASS: NN2 action P', degree, ' ', trim(mode)
contains
    subroutine attach_exterior()
        ! Frozen positive energy/trace blocks: no full BEM qualification.
        type(nn_exterior_t) :: exterior
        real(dp), allocatable :: normal(:, :), drive(:)
        real(dp) :: steklov(3, 3), coil(3), mass
        integer :: j, first, last
        mass = 4*acos(-1.0_dp)/3
        steklov = 0
        do j = 1, 3
            steklov(j, j) = 2*mass
        end do
        allocate (normal(3, lay%na), drive(lay%na))
        normal = 0; drive = 0
        first = findloc(field%fs%sp(1)%bnd, .true., 1)
        last = findloc(field%fs%sp(1)%bnd, .true., 1, back=.true.)
        call require(first > 0 .and. last > first, 'boundary trace columns')
        normal(:, first) = mass*[1.0_dp, 0.3_dp, -0.2_dp]
        normal(:, last) = mass*[-0.4_dp, 0.2_dp, 0.7_dp]
        coil = mass*[0.2_dp, -0.3_dp, 0.4_dp]
        call nn_attach_exterior_operators(field, problem, steklov, normal, &
            coil, drive, exterior, status)
        call require(status == 0, 'frozen exterior attachment')
        call require(norm2(csc_matvec(problem%exterior, direction(:lay%na))) > 0, &
            'nonzero exterior action')
    end subroutine

    subroutine require(ok, label)
        logical, intent(in) :: ok
        character(*), intent(in) :: label
        if (ok) return
        print '(2a,i0)', trim(label), ': status ', status
        error stop 'I3 action gate'
    end subroutine

    subroutine reject(expected, label)
        integer, intent(in) :: expected
        character(*), intent(in) :: label
        call nc_jvp(field, problem, x, direction, action, status)
        call require(status == expected, trim(label)//' JVP')
        call nc_vjp(field, problem, x, weight, adjoint, status)
        call require(status == expected, trim(label)//' VJP')
    end subroutine

    subroutine invalid_inputs()
        problem%chi_perp = 0
        call reject(ISLAND_INVALID, 'zero conductivity')
        problem%chi_perp = 1
        problem%closure%target_class = 'ideal'
        call reject(ISLAND_UNSUPPORTED, 'unsupported ideal claim')
        problem%closure%target_class = 'regularised'
        nullify (problem%g_par_shape_3d)
        problem%g_par_shape => negative_shape
        call reject(ISLAND_BAD_DATA, 'negative pointwise conductivity')
        nullify (problem%g_par_shape)
        problem%source => bad_source
        call reject(ISLAND_BAD_DATA, 'nonfinite declared source')
        problem%source => zero_source
        x = 0; field%f_vac = 0
        call reject(ISLAND_NULL_FIELD, 'magnetic null')
    end subroutine
end program
