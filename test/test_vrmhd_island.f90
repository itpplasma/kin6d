program test_vrmhd_island
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use vrmhd_periodic_slab, only: periodic_slab_t, periodic_slab_build, &
        periodic_scalar
    use vrmhd_island_leaves, only: island_shear_equilibrium, island_seed
    use kin6d_vrmhd_section, only: vrmhd_section_basis, vrmhd_section_trace, &
        vrmhd_section_residual
    use kin6d_vrmhd_parameters, only: vrmhd_parameters_t
    use kin6d_vrmhd_evolve, only: vrmhd_evolve_t, vrmhd_evolve_build, &
        vrmhd_evolve_consistent, vrmhd_evolve_expand
    use fortfem_meridian_feec_n, only: meridian_lagrange_nodes, meridian_load, &
        meridian_point_basis_t, meridian_point_basis_init
    use fortfem_triangle_vector_interpolation, only: &
        interpolate_triangle_nedelec_second_kind
    use fortsparse, only: csc_t, csc_from_triplet, csc_transpose, csc_matvec, &
        fortsparse_status_t, status_ok, &
        FORTSPARSE_LU_SCALAR, FORTSPARSE_LU_SUPERNODAL, sparse_factor_stats, &
        sparse_factor_stats_t
    use fortnum_ode_bdf, only: bdf_options_t, bdf_state_t, bdf_report_t, bdf_init, &
        bdf_integrate_to
    use fortnum_quadrature, only: gauss_legendre_ab
    use fortnum_status, only: fortnum_status_t, FORTNUM_OK
    implicit none
    type(periodic_slab_t), target :: slab
    type(meridian_point_basis_t) :: pb
    type(vrmhd_parameters_t) :: p
    type(csc_t) :: reduction
    real(dp), allocatable :: lift(:), y(:), u(:), source(:, :, :), current_wall(:)
    real(dp), allocatable :: radial_grid(:)
    logical, allocatable :: magnetic_wall(:)
    integer, parameter :: MASS = 1, KINETIC = 2, MAGNETIC = 3, INTERNAL = 4, &
        OHMIC = 5, VISCOUS = 6, COMPRESSION = 7, WALL_POWER = 8, &
        MASS_RATE = 9, KINETIC_RATE = 10, MAGNETIC_RATE = 11, &
        INTERNAL_RATE = 12, NLEDGER = 12
    integer, allocatable :: kept(:)
    real(dp), parameter :: pi = acos(-1.0_dp), half_width = 1.0_dp
    character(32) :: gate, family, arg, preconditioner = 'supernodal'
    character(1024) :: mesh_path
    integer :: ny = 4, nx = 4, steps = 20
    real(dp) :: eta = 0.01_dp, dt = 0.1_dp, seed_amplitude = 0.01_dp
    integer :: status, mesh_unit, io_status
    namelist /radial_mesh/ radial_grid
    call get_command_argument(1, gate)
    call get_command_argument(2, family)
    if (family == '') family = 'mixed'
    call require(family == 'equal' .or. family == 'mixed', 'equal or mixed family')
    call require(any(command_argument_count() == [1, 2, 4, 7, 8, 9, 10]), &
        'arguments: MODE FAMILY [NY NX [ETA DT STEPS [EPSILON [scalar|supernodal [MESH]]]]]')
    if (command_argument_count() >= 3) then
        call get_command_argument(3, arg); read(arg, *) ny
        call get_command_argument(4, arg); read(arg, *) nx
    end if
    if (command_argument_count() >= 5) then
        call get_command_argument(5, arg); read(arg, *) eta
        call get_command_argument(6, arg); read(arg, *) dt
        call get_command_argument(7, arg); read(arg, *) steps
    end if
    if (command_argument_count() >= 8) then
        call get_command_argument(8, arg); read(arg, *) seed_amplitude
    end if
    if (command_argument_count() >= 9) &
        call get_command_argument(9, preconditioner)
    call require(ny > 0 .and. nx > 0, 'positive mesh cell counts')
    if (command_argument_count() == 10) then
        call get_command_argument(10, mesh_path)
        allocate (radial_grid(nx + 1)); radial_grid = -1
        open (newunit=mesh_unit, file=trim(mesh_path), status='old', iostat=io_status)
        call require(io_status == 0, 'open radial mesh namelist')
        read (mesh_unit, nml=radial_mesh, iostat=io_status)
        close (mesh_unit)
        call require(io_status == 0, 'read radial mesh namelist')
    end if
    if (gate == 'graded') then
        call require(command_argument_count() <= 2, 'graded fixture takes only the family')
        nx = 6
        radial_grid = [0.0_dp, 0.125_dp, 0.875_dp, 1.0_dp, 1.0625_dp, 1.5_dp, 2.0_dp]
        call invalid_grid_gate()
    end if
    call require(preconditioner == 'scalar' .or. preconditioner == 'supernodal', &
        'known direct right preconditioner')
    call require(eta > 0 .and. dt > 0 .and. steps > 0, 'positive run parameters')
    call require(seed_amplitude > 0 .and. seed_amplitude < 0.5_dp, &
        'positive seed with separatrix inside the walls')
    call periodic_slab_build(ny, 2, slab, status, radial_cells=nx, &
        radial_length=2*half_width, radial_walls=.true., mixed=family /= 'equal', &
        radial_grid=radial_grid)
    call require(status == 0, 'wall slab construction')
    if (allocated(radial_grid)) &
        print '(a,*(1x,es24.16))', 'radial breakpoints [1]:', slab%radial_nodes
    call wall_reduction()
    call equilibrium()
    call meridian_point_basis_init(slab%op%fs, pb, status)
    call require(status == 0, 'point evaluation basis')
    select case (trim(gate))
    case ('boundary')
        call boundary_gate()
    case ('equilibrium')
        call equilibrium_gate()
    case ('width')
        call width_gate()
    case ('ledger')
        call ledger_oracle()
        call evolve_seed()
    case ('preconditioner')
        call preconditioner_gate()
    case ('graded')
        call boundary_gate()
        call equilibrium_gate()
        call width_gate()
        call ledger_oracle()
    case ('probe', 'induction')
        call evolve_seed()
    case default
        error stop 'unknown island gate'
    end select
contains
    subroutine invalid_grid_gate()
        type(periodic_slab_t) :: invalid
        real(dp) :: nodes(3, 6)
        integer :: failure, i
        nodes(:, 1) = [0.0_dp, 1.0_dp, 2.0_dp]
        call periodic_slab_build(2, 2, invalid, failure, radial_cells=3, &
            radial_length=2.0_dp, radial_grid=nodes(:, 1))
        call require(failure /= 0, 'reject radial grid with wrong size')
        nodes(:, 2) = [0.25_dp, 1.0_dp, 2.0_dp]
        nodes(:, 3) = [0.0_dp, 1.0_dp, 1.75_dp]
        nodes(:, 4) = [0.0_dp, 0.0_dp, 2.0_dp]
        nodes(:, 5) = [0.0_dp, 3.0_dp, 2.0_dp]
        nodes(:, 6) = [0.0_dp, ieee_value(0.0_dp, ieee_quiet_nan), 2.0_dp]
        do i = 2, size(nodes, 2)
            call periodic_slab_build(2, 2, invalid, failure, radial_cells=2, &
                radial_length=2.0_dp, radial_grid=nodes(:, i))
            call require(failure /= 0, 'reject invalid radial breakpoints')
        end do
    end subroutine

    subroutine require(ok, message)
        logical, intent(in) :: ok
        character(*), intent(in) :: message
        if (ok) return
        print '(a)', message
        error stop 'island benchmark failed'
    end subroutine

    subroutine wall_reduction()
        !! No-slip walls and fixed tangential A/J; the y faces are identified.
        type(fortsparse_status_t) :: st
        real(dp), allocatable :: xy(:, :), vals(:)
        integer, allocatable :: rows(:), cols(:)
        logical, allocatable :: fixed(:), magnetic(:), active(:)
        real(dp) :: start(2), finish(2)
        integer :: n1, i, j, t, edge, q, pos, col, first, last, field
        n1 = slab%op%fs%sp(1)%ndof
        allocate (fixed(slab%op%ndof), magnetic(n1), active(slab%nred))
        fixed = .false.; magnetic = .false.
        if (allocated(slab%op%velocity_fs)) then
            call meridian_lagrange_nodes(slab%op%velocity_fs, xy)
        else
            call meridian_lagrange_nodes(slab%op%scalar_fs, xy)
        end if
        do i = 1, size(xy, 2)
            if (.not. wall(xy(1, i))) cycle
            do field = 2, 4
                fixed(slab%op%offset(field) + i) = .true.
            end do
        end do
        deallocate (xy)
        call meridian_lagrange_nodes(slab%op%fs, xy)
        do i = 1, slab%op%fs%n0
            magnetic(slab%op%fs%nn + i) = wall(xy(1, i))
        end do
        do t = 1, slab%op%fs%mm%topology%n_triangles
            do edge = 1, 3
                j = modulo(edge, 3) + 1
                start = slab%op%fs%mm%topology%vertices(:, &
                    slab%op%fs%mm%topology%triangles(edge, t))
                finish = slab%op%fs%mm%topology%vertices(:, &
                    slab%op%fs%mm%topology%triangles(j, t))
                if (.not. wall(start(1)) .or. .not. wall(finish(1))) cycle
                if (abs(start(1) - finish(1)) > 1.0e-12_dp) cycle
                do q = 1, slab%op%fs%ke
                    i = (edge - 1)*slab%op%fs%ke + q
                    magnetic(slab%op%fs%sp(1)%dof(i, t)) = .true.
                end do
            end do
        end do
        magnetic_wall = magnetic
        fixed(slab%op%offset(7) + 1:slab%op%offset(8)) = magnetic
        fixed(slab%op%offset(8) + 1:) = magnetic
        do j = 1, slab%nred
            first = slab%reduction%col_ptr(j)
            last = slab%reduction%col_ptr(j + 1) - 1
            active(j) = .not. any(fixed(slab%reduction%row_idx(first:last)))
            call require(all(fixed(slab%reduction%row_idx(first:last))) &
                .or. active(j), 'periodic class cannot mix wall and free rows')
        end do
        kept = pack([(i, i=1, slab%nred)], active)
        allocate (rows(size(slab%reduction%val)), &
            cols(size(slab%reduction%val)), vals(size(slab%reduction%val)))
        pos = 0
        do col = 1, size(kept)
            j = kept(col)
            do q = slab%reduction%col_ptr(j), slab%reduction%col_ptr(j + 1) - 1
                pos = pos + 1
                rows(pos) = slab%reduction%row_idx(q); cols(pos) = col
                vals(pos) = slab%reduction%val(q)
            end do
        end do
        call csc_from_triplet(slab%op%ndof, size(kept), rows(:pos), cols(:pos), &
            vals(:pos), reduction, st)
        call require(status_ok(st), 'wall reduction')
    end subroutine

    logical function wall(x)
        real(dp), intent(in) :: x
        wall = min(abs(x), abs(x - 2*half_width)) < 1.0e-12_dp
    end function

    subroutine equilibrium()
        complex(dp), allocatable :: load(:)
        real(dp), allocatable :: xy(:, :), basis(:, :), moments(:)
        real(dp) :: rho, ay, az, by, jz, bz, point(2), weight
        integer :: n0, i, t, q, a0
        n0 = slab%op%scalar_fs%n0
        allocate (lift(slab%op%ndof), u(slab%op%ndof), y(size(kept)), &
            source(6, slab%op%fs%nq, slab%op%fs%mm%topology%n_triangles), &
            basis(36, slab%op%nlocal))
        lift = 0; y = 0; source = 0
        call meridian_lagrange_nodes(slab%op%scalar_fs, xy)
        do i = 1, n0
            call island_shear_equilibrium(xy(1, i) - half_width, &
                rho, ay, az, by, jz, bz)
            lift(i) = rho
        end do
        lift(slab%op%offset(5) + 1:slab%op%offset(7)) = 1
        call meridian_load(slab%op%fs, 1, equilibrium_potential, load, status)
        call require(status == 0, 'equilibrium potential load')
        a0 = slab%op%offset(7)
        do t = 1, slab%op%fs%mm%topology%n_triangles
            call interpolate_triangle_nedelec_second_kind( &
                slab%op%fs%mm%topology%vertices(:, &
                slab%op%fs%mm%topology%triangles(:, t)), &
                slab%op%fs%k, 2*slab%op%fs%k + 2, &
                equilibrium_section_potential, moments, status)
            call require(status == 0, 'equilibrium defining moments')
            lift(a0 + slab%op%fs%sp(1)%dof(:slab%op%fs%nln, t)) = &
                moments*slab%op%fs%sp(1)%sgn(:slab%op%fs%nln, t)
        end do
        deallocate (xy)
        call meridian_lagrange_nodes(slab%op%fs, xy)
        do i = 1, size(xy, 2)
            call island_shear_equilibrium(xy(1, i) - half_width, &
                rho, ay, az, by, jz, bz)
            lift(a0 + slab%op%fs%nn + i) = az
        end do
        call require(maxval(abs(csc_matvec(slab%op%fs%mass(1), &
            lift(slab%op%offset(7) + 1:slab%op%offset(8))) - real(load, dp))) &
            < 1.0e-11_dp*maxval(abs(load)), 'full projection mass residual')
        ! Constant axial current is an orbit-node coefficient, not an edge moment.
        lift(slab%op%offset(8) + slab%op%fs%nn + 1:) = 1
        current_wall = lift(slab%op%offset(8) + 1:)
        where (.not. magnetic_wall) current_wall = 0
        source(3, :, :) = 1
        p%mi = 1; p%mu0 = 1; p%eta0 = eta; p%nu0 = eta
        u = lift
        do t = 1, slab%op%fs%mm%topology%n_triangles
            do q = 1, slab%op%fs%nq
                call vrmhd_section_basis(slab%op, t, q, point, weight, basis, status)
                call require(status == 0, 'equilibrium basis')
                call island_shear_equilibrium(point(1) - half_width, &
                    rho, ay, az, by, jz, bz)
                call require(maxval(abs(matmul(basis([1, 6, 7, 26, 10, 27], :), &
                    lift(slab%op%map(:, t))) - [rho, ay, az, by, jz, bz])) &
                    < 1.0e-11_dp, 'physical equilibrium represented exactly')
            end do
        end do
    end subroutine

    subroutine equilibrium_section_potential(x, y, a)
        real(dp), intent(in) :: x, y
        real(dp), intent(out) :: a(2)
        real(dp) :: rho, ay, az, by, jz, bz
        call island_shear_equilibrium(x - half_width, rho, ay, az, by, jz, bz)
        a = [0.0_dp, ay]
    end subroutine

    subroutine equilibrium_potential(x, a)
        real(dp), intent(in) :: x(2)
        complex(dp), intent(out) :: a(3)
        real(dp) :: rho, ay, az, by, jz, bz
        call island_shear_equilibrium(x(1) - half_width, rho, ay, az, by, jz, bz)
        a = cmplx([0.0_dp, az, ay], 0.0_dp, dp)
    end subroutine

    subroutine boundary_gate()
        real(dp), allocatable :: full(:), reduced(:), scalar(:), xy(:, :)
        real(dp) :: a(3), b(3), x, nodes(4), weights(4)
        integer :: i, cell, n0
        n0 = slab%op%scalar_fs%n0
        allocate (reduced(slab%nred), scalar(slab%ns))
        call periodic_scalar(slab, radial_coordinate, scalar)
        reduced = 0; reduced(:slab%ns) = scalar
        full = csc_matvec(slab%reduction, reduced)
        call meridian_lagrange_nodes(slab%op%scalar_fs, xy)
        call require(maxval(abs(full(:n0) - xy(1, :))) < 1.0e-13_dp, &
            'opposite radial walls must retain independent density traces')
        deallocate (reduced)
        allocate (reduced(reduction%ncol))
        reduced = [(sin(real(37*i, dp)), i=1, size(reduced))]
        full = csc_matvec(reduction, reduced)
        call gauss_legendre_ab(4, 0.0_dp, 1.0_dp, nodes, weights)
        ! H(curl) tangential traces are defined on facet interiors; the vector
        ! value at a shared vertex need not agree across adjacent facets.
        do cell = 1, slab%nx
            do i = 1, 4
                x = slab%radial_nodes(cell) + nodes(i) &
                    *(slab%radial_nodes(cell + 1) - slab%radial_nodes(cell))
                a = potential([x, 0.0_dp], full)
                b = potential([x, 2*pi], full)
                call require(maxval(abs(a([1, 3]) - b([1, 3]))) < 1.0e-10_dp, &
                    'periodic tangential potential traces and edge orientation')
            end do
        end do
        do cell = 1, slab%n
            do i = 1, 4
                x = 2*pi*(cell - 1 + nodes(i))/slab%n
                a = potential([0.0_dp, x], full)
                b = potential([2*half_width, x], full)
                call require(maxval(abs([a([2, 3]), b([2, 3])])) < 1.0e-10_dp, &
                    'conducting-wall tangential perturbation traces')
            end do
        end do
    end subroutine

    pure real(dp) function radial_coordinate(x) result(value)
        real(dp), intent(in) :: x(2)
        value = x(1)
    end function

    subroutine equilibrium_gate()
        type(vrmhd_evolve_t) :: problem
        type(bdf_options_t) :: options
        type(bdf_state_t) :: state
        type(bdf_report_t) :: report
        type(fortnum_status_t) :: st
        real(dp), allocatable :: derivative(:), residual(:), zero(:)
        integer :: step
        allocate (derivative(size(y)), residual(size(y)), zero(size(y)))
        zero = 0
        call vrmhd_evolve_build(problem, slab%op, p, source, lift, reduction, status)
        call require(status == 0, 'wall evolution build')
        call problem%evaluate(0.0_dp, y, zero, residual, status)
        call require(status == 0, 'equilibrium residual')
        print '(a,es14.5)', 'equilibrium max residual [1]: ', maxval(abs(residual))
        call require(maxval(abs(residual)) < 1.0e-11_dp, &
            'undriven momentum and source-balanced Ohm equilibrium')
        call vrmhd_evolve_consistent(problem, 0.0_dp, y, derivative, status)
        call require(status == 0, 'consistent equilibrium derivative')
        print '(a,es14.5)', 'equilibrium max derivative [1]: ', &
            maxval(abs(derivative))
        call require(maxval(abs(derivative)) < 1.0e-8_dp, 'stationary initial state')
        options%adaptive = .false.; options%initial_step = 0.1_dp
        options%max_step = 0.1_dp; options%rtol = 1.0e-9_dp
        options%atol = 1.0e-11_dp
        call bdf_init(state, 0.0_dp, y, options, st, derivative=derivative)
        call require(st%code == FORTNUM_OK, 'equilibrium BDF start')
        do step = 1, 2
            call bdf_integrate_to(problem, state, step*0.1_dp, options, report, st)
            call require(st%code == FORTNUM_OK, 'equilibrium BDF step')
        end do
        call vrmhd_evolve_expand(problem, state%y, u)
        print '(a,es14.5)', 'equilibrium BDF max drift [1]: ', maxval(abs(u - lift))
        call require(maxval(abs(u - lift)) < 1.0e-9_dp, 'equilibrium remains steady')
    end subroutine

    subroutine evolve_seed()
        type(vrmhd_evolve_t) :: problem
        type(bdf_options_t) :: options
        type(bdf_state_t) :: state
        type(bdf_report_t) :: report
        type(fortnum_status_t) :: st
        real(dp), allocatable :: xy(:, :), derivative(:), full_dot(:), rhs(:)
        real(dp), allocatable :: normal(:)
        type(csc_t) :: restriction
        type(fortsparse_status_t) :: sparse_status
        type(sparse_factor_stats_t) :: factor_stats
        integer :: newton, krylov, accepted, rejected
        real(dp) :: az, bx, by, jz, rate, expected, width, delta, relative
        real(dp) :: totals(NLEDGER), initial(NLEDGER), previous(NLEDGER)
        real(dp) :: wall_integral, mechanical_integral, mechanical_power, old_power
        integer :: i, j, k, step
        u = lift
        if (gate == 'ledger') then
            call meridian_lagrange_nodes(slab%op%scalar_fs, xy)
            u(slab%op%offset(5) + 1:slab%op%offset(6)) = &
                1 + 0.1_dp*cos(xy(2, :))
            deallocate (xy)
        end if
        call meridian_lagrange_nodes(slab%op%fs, xy)
        do i = 1, size(xy, 2)
            call island_seed(xy(1, i) - half_width, xy(2, i), half_width, &
                1.0_dp, seed_amplitude, az, bx, by, jz)
            u(slab%op%offset(7) + slab%op%fs%nn + i) = az
        end do
        call csc_transpose(reduction, restriction, sparse_status)
        call require(status_ok(sparse_status), 'seed restriction')
        rhs = csc_matvec(restriction, u - lift)
        allocate (normal(size(y)), derivative(size(y)), full_dot(slab%op%ndof))
        do j = 1, reduction%ncol
            normal(j) = 0
            do k = reduction%col_ptr(j), reduction%col_ptr(j + 1) - 1
                normal(j) = normal(j) + reduction%val(k)**2
            end do
        end do
        y = rhs/normal
        call vrmhd_evolve_build(problem, slab%op, p, source, lift, reduction, status)
        call require(status == 0, 'island evolution build')
        if (preconditioner == 'scalar') &
            problem%solver%options%lu_method = FORTSPARSE_LU_SCALAR
        problem%lag = 8
        problem%shift_tolerance = 0.5_dp
        call vrmhd_evolve_consistent(problem, 0.0_dp, y, derivative, status)
        call require(status == 0, 'consistent island current and derivative')
        full_dot = csc_matvec(reduction, derivative)
        rate = flux([half_width, 0.0_dp], full_dot) &
            - flux([half_width, pi], full_dot)
        ! At rest, induction gives A_z,t = eta Laplacian(delta A_z).
        expected = -2*eta*seed_amplitude*(6/half_width**2 + 1)
        relative = abs(rate/expected - 1)
        if (gate /= 'ledger') then
            print '(a,2i5,3es15.6)', 'initial induction: ny,nx,rate,exact,relative [1] ', &
                ny, nx, rate, expected, relative
        end if
        call vrmhd_evolve_expand(problem, y, u)
        call measure(u, full_dot, initial)
        if (gate == 'ledger') then
            call check_ledger(u, full_dot, initial)
            return
        end if
        if (gate == 'induction') then
            call require(relative < 0.02_dp, 'initial resistive flux derivative')
            return
        end if
        options%adaptive = .false.; options%initial_step = dt
        options%max_step = dt; options%rtol = 1.0e-8_dp
        options%atol = 1.0e-10_dp
        call bdf_init(state, 0.0_dp, y, options, st, derivative=derivative)
        call require(st%code == FORTNUM_OK, 'island BDF start')
        delta = -2/tanh(half_width)
        print '(a,es14.5)', 'outer tearing index [1]: ', delta
        print '(a,es14.5)', 'Rutherford slope, not an admitted gate [1]: ', &
            1.22_dp*eta*delta
        newton = 0; krylov = 0; accepted = 0; rejected = 0
        wall_integral = 0; mechanical_integral = 0
        previous = initial
        print '(a)', 'step,time,full_cut_width,flux_ox,max_velocity,mass,'// &
            'total_energy,wall_power,energy_defect,mechanical_power,mechanical_defect'
        flush(6)
        do step = 0, steps
            if (step > 0) then
                call bdf_integrate_to(problem, state, step*dt, options, report, st)
                if (st%code /= FORTNUM_OK) print *, report
                call require(st%code == FORTNUM_OK, 'island BDF step')
                newton = newton + report%newton_iterations
                krylov = krylov + report%krylov_iterations
                accepted = accepted + report%accepted
                rejected = rejected + report%rejected_newton
                if (report%accepted /= 1 .or. report%rejected_newton /= 0) &
                    print '(a,2i9)', 'accepted steps, nonlinear retries: ', &
                    report%accepted, report%rejected_newton
                call require(report%accepted == 1 .and. report%rejected_newton == 0, &
                    'fixed time grid without nonlinear retries')
            end if
            call vrmhd_evolve_expand(problem, state%y, u)
            width = symmetric_width()
            full_dot = csc_matvec(reduction, state%derivative)
            call measure(u, full_dot, totals)
            mechanical_power = totals(COMPRESSION) - totals(OHMIC) &
                - totals(VISCOUS) + totals(WALL_POWER)
            if (step > 0) then
                wall_integral = wall_integral + dt/2 &
                    *(previous(WALL_POWER) + totals(WALL_POWER))
                old_power = previous(COMPRESSION) - previous(OHMIC) &
                    - previous(VISCOUS) + previous(WALL_POWER)
                mechanical_integral = mechanical_integral + dt/2 &
                    *(old_power + mechanical_power)
            end if
            previous = totals
            print '(i0,10(a,es24.16))', step, ',', step*dt, ',', width, ',', &
                flux([half_width, 0.0_dp]) - flux([half_width, pi]), ',', &
                maxval(abs(u(slab%op%offset(2) + 1:slab%op%offset(5)))), &
                ',', totals(MASS), ',', sum(totals(KINETIC:INTERNAL)), ',', &
                totals(WALL_POWER), ',', &
                sum(totals(KINETIC:INTERNAL) - initial(KINETIC:INTERNAL)) &
                - wall_integral, ',', mechanical_power, ',', &
                sum(totals(KINETIC:MAGNETIC) - initial(KINETIC:MAGNETIC)) &
                - mechanical_integral
            flush(6)
        end do
        print '(a,3i9)', 'Newton, Krylov, factor setups: ', &
            newton, krylov, problem%factorizations
        print '(a,2i9)', 'accepted steps, nonlinear retries: ', accepted, rejected
        call sparse_factor_stats(problem%solver, factor_stats, sparse_status)
        call require(status_ok(sparse_status), 'direct factor statistics')
        print '(a,a,a,l1,1x,a)', 'last factor: ', trim(factor_stats%native_factor_family), &
            ', fallback ', factor_stats%lu_fallback, trim(factor_stats%lu_fallback_reason)
        print '(a,2i14)', 'last factor entries, retained owned bytes: ', &
            factor_stats%factor_nnz, factor_stats%host_current_bytes
        if (factor_stats%backward_error >= 0) &
            print '(a,es17.9)', 'last solve backward error: ', &
            factor_stats%backward_error
        if (factor_stats%perturbed_pivots >= 0) &
            print '(a,i9)', 'last factor perturbed pivots: ', factor_stats%perturbed_pivots
        if (gate == 'preconditioner' .and. preconditioner == 'supernodal') &
            call require(factor_stats%lu_method_used == FORTSPARSE_LU_SUPERNODAL, &
            'exercise supernodal factor rather than scalar fallback')
    end subroutine

    subroutine preconditioner_gate()
        real(dp), allocatable :: reference(:), basis(:, :)
        real(dp) :: point(2), weight, difference(36), perturbation(36)
        real(dp) :: error, scale, reference_width, reference_flux
        integer :: t, q
        dt = 0.125_dp; steps = 4
        preconditioner = 'scalar'
        call evolve_seed()
        reference = u
        reference_width = symmetric_width()
        reference_flux = flux([half_width, 0.0_dp]) - flux([half_width, pi])
        preconditioner = 'supernodal'
        call evolve_seed()
        allocate (basis(36, slab%op%nlocal))
        error = 0; scale = 0
        do t = 1, slab%op%fs%mm%topology%n_triangles
            do q = 1, slab%op%fs%nq
                call vrmhd_section_basis(slab%op, t, q, point, weight, basis, status)
                call require(status == 0, 'preconditioner physical norm')
                difference = matmul(basis, u(slab%op%map(:, t)) &
                    - reference(slab%op%map(:, t)))
                perturbation = matmul(basis, reference(slab%op%map(:, t)) &
                    - lift(slab%op%map(:, t)))
                error = error + weight*(sum(difference(:12)**2) &
                    + sum(difference(25:27)**2))
                scale = scale + weight*(sum(perturbation(:12)**2) &
                    + sum(perturbation(25:27)**2))
            end do
        end do
        call require(scale > 0, 'nontrivial reference perturbation')
        print '(a,es17.9)', 'relative physical state difference: ', sqrt(error/scale)
        call require(sqrt(error/scale) < 1.0e-6_dp, 'full nonlinear physical state')
        call require(abs(symmetric_width()/reference_width - 1) < 1.0e-7_dp, &
            'preconditioner-independent full cut width')
        call require(abs((flux([half_width, 0.0_dp]) - flux([half_width, pi])) &
            /reference_flux - 1) < 1.0e-6_dp, 'preconditioner-independent island flux')
    end subroutine

    subroutine measure(field, derivative, totals)
        !! Independent physical integrals, including the prescribed current port.
        real(dp), intent(in) :: field(:), derivative(:)
        real(dp), intent(out) :: totals(NLEDGER)
        real(dp), allocatable :: basis(:, :)
        real(dp) :: point(2), weight, v(36), vd(36), g(3, 3), strain(3, 3)
        real(dp) :: eta, nu, pressure, divergence, wall_field(3), vxB(3)
        integer :: t, q, i, first
        allocate (basis(36, slab%op%nlocal))
        totals = 0
        first = slab%op%local_offset(8) + 1
        do t = 1, slab%op%fs%mm%topology%n_triangles
            do q = 1, slab%op%fs%nq
                call vrmhd_section_basis(slab%op, t, q, point, weight, basis, status)
                call require(status == 0, 'island ledger basis')
                v = matmul(basis, field(slab%op%map(:, t)))
                vd = matmul(basis, derivative(slab%op%map(:, t)))
                g = reshape(v(16:24), [3, 3])
                divergence = g(1, 1) + g(2, 2) + g(3, 3)
                strain = (g + transpose(g))/2
                do i = 1, 3
                    strain(i, i) = strain(i, i) - divergence/3
                end do
                eta = p%eta0*(v(11)/p%t0)**(-1.5_dp)
                nu = p%nu0*(v(12)/p%t0)**(-1.5_dp)
                pressure = v(1)*(v(11) + v(12))/p%mi
                wall_field = matmul(basis(8:10, first:), &
                    current_wall(slab%op%map(first:, t) - slab%op%offset(8)))
                vxB = [v(3)*v(27) - v(4)*v(26), &
                    v(4)*v(25) - v(2)*v(27), v(2)*v(26) - v(3)*v(25)]
                totals(MASS) = totals(MASS) + weight*v(1)
                totals(KINETIC) = totals(KINETIC) + weight*v(1)*sum(v(2:4)**2)/2
                totals(MAGNETIC) = totals(MAGNETIC) + weight*sum(v(25:27)**2)/(2*p%mu0)
                totals(INTERNAL) = totals(INTERNAL) + weight*1.5_dp*pressure
                totals(OHMIC) = totals(OHMIC) &
                    + weight*eta*dot_product(v(8:10), v(8:10) - source(:3, q, t))
                totals(VISCOUS) = totals(VISCOUS) + weight*2*nu*sum(strain*strain)
                totals(COMPRESSION) = totals(COMPRESSION) + weight*pressure*divergence
                totals(WALL_POWER) = totals(WALL_POWER) + weight*dot_product(wall_field, &
                    vd(5:7) - vxB + eta*(v(8:10) - source(:3, q, t)))
                totals(MASS_RATE) = totals(MASS_RATE) + weight*vd(1)
                totals(KINETIC_RATE) = totals(KINETIC_RATE) + weight &
                    *(vd(1)*sum(v(2:4)**2)/2 + v(1)*dot_product(v(2:4), vd(2:4)))
                totals(MAGNETIC_RATE) = totals(MAGNETIC_RATE) &
                    + weight*dot_product(v(25:27), vd(25:27))/p%mu0
                totals(INTERNAL_RATE) = totals(INTERNAL_RATE) + weight*1.5_dp/p%mi &
                    *(vd(1)*(v(11) + v(12)) + v(1)*(vd(11) + vd(12)))
            end do
        end do
    end subroutine

    subroutine ledger_oracle()
        real(dp) :: values(NLEDGER), exact(4), rates(4)
        real(dp), allocatable :: field(:), derivative(:)
        integer :: first, last, j
        field = lift
        do j = 2, 4
            first = slab%op%offset(j) + 1; last = slab%op%offset(j + 1)
            field(first:last) = real(j, dp)/10
        end do
        derivative = field
        derivative(:slab%op%offset(2)) = 0.1_dp*field(:slab%op%offset(2))
        derivative(slab%op%offset(2) + 1:slab%op%offset(5)) = &
            0.3_dp*field(slab%op%offset(2) + 1:slab%op%offset(5))
        derivative(slab%op%offset(5) + 1:slab%op%offset(7)) = &
            0.2_dp*field(slab%op%offset(5) + 1:slab%op%offset(7))
        derivative(slab%op%offset(7) + 1:) = &
            0.4_dp*field(slab%op%offset(7) + 1:)
        call measure(field, derivative, values)
        exact = [11*pi/3, 0.145_dp*11*pi/3, 26*pi/3, 11*pi]
        rates = exact*[0.1_dp, 0.7_dp, 0.8_dp, 0.3_dp]
        call require(maxval(abs(values(MASS:INTERNAL)/exact - 1)) < 1.0e-10_dp, &
            'analytic mass and energy integrals')
        call require(maxval(abs(values(MASS_RATE:INTERNAL_RATE)/rates - 1)) &
            < 1.0e-10_dp, 'independent physical chain-rule rates')
    end subroutine

    subroutine check_ledger(field, derivative, values)
        real(dp), intent(in) :: field(:), derivative(:), values(NLEDGER)
        real(dp), allocatable :: residual(:), unused(:), zero(:)
        real(dp) :: scale, defect, reaction, mechanical, thermal
        integer :: a0, j0
        allocate (residual(size(field)), unused(size(field)), zero(size(field)))
        zero = 0
        call vrmhd_section_residual(slab%op, p, field, derivative, zero, zero, &
            source, residual, unused, status)
        call require(status == 0, 'independent wall reaction residual')
        a0 = slab%op%offset(7); j0 = slab%op%offset(8)
        reaction = dot_product(current_wall, residual(a0 + 1:j0))
        scale = max(sum(abs(values(OHMIC:WALL_POWER))), 1.0e-10_dp)
        defect = sum(values(KINETIC_RATE:INTERNAL_RATE)) - values(WALL_POWER)
        print '(a,4es17.9)', 'mass rate, energy rate, wall power, defect [1]: ', &
            values(MASS_RATE), sum(values(KINETIC_RATE:INTERNAL_RATE)), &
            values(WALL_POWER), defect
        call require(abs(values(MASS_RATE)) < 1.0e-10_dp, 'closed mass balance')
        call require(abs(values(WALL_POWER) - reaction) < 1.0e-10_dp*scale, &
            'wall power from physical pairings matches constrained reaction')
        mechanical = sum(values(KINETIC_RATE:MAGNETIC_RATE)) &
            - values(COMPRESSION) + values(OHMIC) + values(VISCOUS) &
            - values(WALL_POWER)
        thermal = values(INTERNAL_RATE) + values(COMPRESSION) &
            - values(OHMIC) - values(VISCOUS)
        call require(max(abs(defect), abs(mechanical), abs(thermal)) &
            < 1.0e-8_dp*scale, 'source-accounted energy exchanges')
        call require(abs(values(WALL_POWER)) > 1.0e-10_dp, 'nontrivial wall port')
    end subroutine

    subroutine width_gate()
        real(dp), allocatable :: xy(:, :)
        real(dp) :: rho, ay, az, by, jz, bz, width, amplitude, expected
        integer :: i
        call meridian_lagrange_nodes(slab%op%fs, xy)
        amplitude = 0.01_dp
        u = lift
        do i = 1, slab%op%fs%n0
            ! A quadratic x profile with a nodal y perturbation has an exact
            ! separatrix cut, independently of the interpolant between y nodes.
            call island_shear_equilibrium(xy(1, i) - half_width, &
                rho, ay, az, by, jz, bz)
            u(slab%op%offset(7) + slab%op%fs%nn + i) = az + amplitude*cos(xy(2, i))
        end do
        width = symmetric_width()
        expected = 4*sqrt(amplitude)
        print '(a,2es14.5)', 'full cut width, analytic [1]: ', width, expected
        call require(abs(width - expected) < 1.0e-11_dp, 'full-width convention')
        u(slab%op%offset(7) + slab%op%fs%nn + 1:slab%op%offset(8)) = &
            u(slab%op%offset(7) + slab%op%fs%nn + 1:slab%op%offset(8)) + 7.0_dp
        call require(abs(symmetric_width() - width) < 1.0e-11_dp, &
            'flux gauge offset cannot change width')
    end subroutine

    real(dp) function flux(point, field) result(value)
        real(dp), intent(in) :: point(2)
        real(dp), intent(in), optional :: field(:)
        real(dp) :: a(3)
        a = potential(point, field)
        value = a(3)
    end function

    function potential(point, field) result(value)
        real(dp), intent(in) :: point(2)
        real(dp), intent(in), optional :: field(:)
        real(dp) :: value(3)
        real(dp) :: start(2), e1(2), e2(2), ref(2), det, x(2)
        real(dp) :: trace(3, slab%op%fs%sp(1)%nloc)
        integer :: t, i(3)
        do t = 1, slab%op%fs%mm%topology%n_triangles
            i = slab%op%fs%mm%topology%triangles(:, t)
            start = slab%op%fs%mm%topology%vertices(:, i(1))
            e1 = slab%op%fs%mm%topology%vertices(:, i(2)) - start
            e2 = slab%op%fs%mm%topology%vertices(:, i(3)) - start
            det = e1(1)*e2(2) - e1(2)*e2(1)
            ref(1) = ((point(1) - start(1))*e2(2) &
                - (point(2) - start(2))*e2(1))/det
            ref(2) = (e1(1)*(point(2) - start(2)) &
                - e1(2)*(point(1) - start(1)))/det
            if (minval(ref) < -1.0e-12_dp .or. sum(ref) > 1 + 1.0e-12_dp) cycle
            call vrmhd_section_trace(slab%op, pb, t, ref, x, trace, status)
            call require(status == 0, 'flux point evaluation')
            if (present(field)) then
                value = matmul(trace, &
                    field(slab%op%offset(7) + slab%op%fs%sp(1)%dof(:, t)))
            else
                value = matmul(trace, &
                    u(slab%op%offset(7) + slab%op%fs%sp(1)%dof(:, t)))
            end if
            return
        end do
        error stop 'point outside slab'
    end function

    real(dp) function symmetric_width() result(width)
        !! Full separatrix width on y=0, for the stated reflection-symmetric seed.
        real(dp) :: separatrix, lo, hi, middle, crossing(2)
        integer :: side, iteration
        separatrix = flux([half_width, pi])
        call require(flux([half_width, 0.0_dp]) > separatrix, 'positive island flux')
        do side = 1, 2
            lo = 0; hi = half_width
            call require(flux([half_width + (2*side - 3)*hi, 0.0_dp]) &
                < separatrix, 'separatrix crossing bracket')
            do iteration = 1, 48
                middle = (lo + hi)/2
                if (flux([half_width + (2*side - 3)*middle, 0.0_dp]) &
                    > separatrix) then
                    lo = middle
                else
                    hi = middle
                end if
            end do
            crossing(side) = (lo + hi)/2
        end do
        width = sum(crossing)
    end function
end program test_vrmhd_island
