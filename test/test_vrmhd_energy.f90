program test_vrmhd_energy
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use milca_gs_symmetry, only: gs_symmetry_t, GS_ROTATION, GS_SCREW, &
        GS_CHART_RZ, GS_CHART_CARTESIAN
    use milca_vrmhd_parameters, only: vrmhd_parameters_t
    use milca_vrmhd_section, only: vrmhd_section_t, vrmhd_section_build, &
        vrmhd_section_basis, vrmhd_section_residual, vrmhd_section_jacobian
    use vrmhd_wall_oracle, only: current_wall_dofs
    use milca_vrmhd_boundary, only: vrmhd_current_lift
    use milca_vrmhd_drive, only: vrmhd_drive_build, VRMHD_CURRENT_SOURCE, &
        VRMHD_LOOP_VOLTAGE, VRMHD_LOOP_VOLTAGE_SOURCE
    use fortfem_mesh_2d, only: mesh_2d_t
    use fortfem_meridian_curved_mesh, only: meridian_mesh_t, meridian_mesh_straight
    use fortsparse, only: csc_t
    implicit none
    integer :: kind, seed_size, i, convention, family_degree
    character(16) :: mode
    integer, allocatable :: seed(:)
    logical :: passed
    real(dp) :: defect

    call get_command_argument(1, mode)
    family_degree = 0
    if (mode == 'mixed1') family_degree = 1
    if (mode == 'mixed2') family_degree = 2
    call random_seed(size=seed_size)
    allocate (seed(seed_size))
    seed = [(104729 + 37*i, i=1, seed_size)]
    call random_seed(put=seed)
    passed = .true.
    do kind = GS_ROTATION, GS_SCREW
        call hydrodynamic_ledger(kind, defect)
        passed = passed .and. defect <= 1.0e-12_dp
    end do
    call require(passed, 'closed hydrodynamic energy ledger exceeds 1e-12')
    do kind = GS_ROTATION, GS_SCREW
        do convention = VRMHD_CURRENT_SOURCE, VRMHD_LOOP_VOLTAGE_SOURCE
            call full_ledger(kind, convention, defect)
            passed = passed .and. defect <= 1.0e-12_dp
        end do
    end do
    call require(passed, 'full semidiscrete energy ledger exceeds 1e-12')
contains
    subroutine require(ok, message)
        logical, intent(in) :: ok
        character(*), intent(in) :: message
        if (.not. ok) then
            print '(a)', message
            error stop 'RM2 energy gate failed'
        end if
    end subroutine

    subroutine hydrodynamic_ledger(kind, relative_defect)
        integer, intent(in) :: kind
        real(dp), intent(out) :: relative_defect
        type(vrmhd_section_t) :: op
        type(vrmhd_parameters_t) :: p
        type(mesh_2d_t) :: mesh
        type(meridian_mesh_t) :: mm
        type(gs_symmetry_t) :: symmetry
        real(dp), allocatable :: u(:), ud(:), zero(:), r(:), unused(:)
        real(dp), allocatable :: mass(:, :), original(:, :), rhs(:, :)
        real(dp), allocatable :: source(:, :, :), basis(:, :)
        integer, allocatable :: free_index(:), pivot(:)
        logical, allocatable :: active(:)
        real(dp) :: x(2), weight, v(36), vd(36), gv(3, 3), pressure
        real(dp) :: kinetic, internal, advection, continuity, pressure_work
        real(dp) :: scale, solve_defect, ledger_defect, shift
        integer :: status, n0, nf, j, q, t, info
        interface
            subroutine dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
                import :: dp
                integer, intent(in) :: n, nrhs, lda, ldb
                integer, intent(out) :: ipiv(*), info
                real(dp), intent(inout) :: a(lda, *), b(ldb, *)
            end subroutine
        end interface

        symmetry = gs_symmetry_t(kind, GS_CHART_CARTESIAN, 0.0_dp)
        shift = 0
        if (kind == GS_ROTATION) then
            shift = 2
            symmetry%chart = GS_CHART_RZ
        end if
        if (kind == GS_SCREW) symmetry%k = 0.7_dp
        call mesh%create_rectangular(3, 3, shift - 0.4_dp, shift + 0.6_dp, &
            -0.3_dp, 0.7_dp)
        call meridian_mesh_straight(mesh, 1, mm, status)
        call require(status == 0, 'energy mesh')
        if (family_degree > 0) then
            call vrmhd_section_build(mm, symmetry, family_degree + 1, 1, op, status, &
                scalar_degree=family_degree, velocity_degree=family_degree + 1)
        else
            call vrmhd_section_build(mm, symmetry, 1, 1, op, status)
        end if
        call require(status == 0, 'energy section')
        n0 = op%scalar_fs%n0
        p%mi = 1.2_dp; p%mu0 = 0.8_dp
        allocate (u(op%ndof), ud(op%ndof), zero(op%ndof), r(op%ndof), &
            unused(op%ndof), active(op%ndof), &
            source(6, op%fs%nq, op%fs%mm%topology%n_triangles), &
            basis(36, op%nlocal))
        u = 0; ud = 0; zero = 0; source = 0
        call random_number(u(:op%offset(7)))
        u(:n0) = 1 + 0.4_dp*u(:n0)
        u(n0 + 1:op%offset(5)) = 0.4_dp*(u(n0 + 1:op%offset(5)) - 0.5_dp)
        u(op%offset(5) + 1:op%offset(7)) = 1 + 0.2_dp*u(op%offset(5) + 1:op%offset(7))
        where (op%fixed(n0 + 1:op%offset(5))) u(n0 + 1:op%offset(5)) = 0

        ! No-flux scalars retain constants; only velocity has wall constraints.
        active = .false.
        active(:op%offset(7)) = .true.
        active(n0 + 1:op%offset(5)) = .not. op%fixed(n0 + 1:op%offset(5))
        nf = count(active)
        allocate (free_index(nf), pivot(nf), mass(nf, nf), original(nf, nf), &
            rhs(nf, 1))
        free_index = pack([(j, j=1, op%ndof)], active)
        call vrmhd_section_residual(op, p, u, zero, zero, zero, source, &
            r, unused, status)
        call require(status == 0, 'energy spatial residual')
        rhs(:, 1) = -r(free_index)
        call free_mass(op, p, u, free_index, mass, status)
        call require(status == 0, 'energy mass assembly')
        original = mass
        call dgesv(nf, 1, mass, nf, pivot, rhs, nf, info)
        call require(info == 0, 'energy mass solve')
        ud(free_index) = rhs(:, 1)
        solve_defect = maxval(abs(matmul(original, rhs(:, 1)) + r(free_index))) &
            /maxval(abs(r(free_index)))
        call require(solve_defect < 1.0e-13_dp, 'energy mass solve accuracy')

        kinetic = 0; internal = 0; advection = 0; continuity = 0
        pressure_work = 0
        do t = 1, op%fs%mm%topology%n_triangles
            do q = 1, op%fs%nq
                call vrmhd_section_basis(op, t, q, x, weight, basis, status)
                call require(status == 0, 'energy quadrature basis')
                v = matmul(basis, u(op%map(:, t)))
                vd = matmul(basis, ud(op%map(:, t)))
                gv = reshape(v(16:24), [3, 3])
                pressure = v(1)*(v(11) + v(12))/p%mi
                ! Differentiate the physical energy, independently of the leaves.
                kinetic = kinetic + weight*( &
                    0.5_dp*vd(1)*sum(v(2:4)**2) + v(1)*dot_product(v(2:4), vd(2:4)))
                internal = internal + weight*1.5_dp/p%mi*( &
                    vd(1)*(v(11) + v(12)) + v(1)*(vd(11) + vd(12)))
                advection = advection - weight*0.5_dp*v(1)* &
                    dot_product(v(2:4), matmul(gv, v(2:4)))
                continuity = continuity + weight*0.5_dp*v(1)* &
                    dot_product(v(2:4), matmul(gv, v(2:4)))
                pressure_work = pressure_work + weight*pressure* &
                    (gv(1, 1) + gv(2, 2) + gv(3, 3))
            end do
        end do
        scale = abs(advection) + abs(continuity) + 2*abs(pressure_work)
        call require(scale > tiny(scale), 'nontrivial energy exchanges')
        ledger_defect = max(abs(kinetic - advection - continuity - pressure_work), &
            abs(internal + pressure_work), abs(advection + continuity))/scale
        call require(ledger_defect < 1.0e-12_dp, 'independent energy accounting')
        relative_defect = abs(kinetic + internal)/scale
        print '(a,i2)', 'symmetry (1 rotation, 2 translation, 3 screw): ', kind
        if (kind == GS_ROTATION) then
            print '(a)', 'power in W/rad (per unit rotation angle)'
        else
            print '(a)', 'power in W/m (per unit axial orbit coordinate)'
        end if
        print '(a,es24.16)', 'forward skew advection work:       ', advection
        print '(a,es24.16)', 'reverse skew advection work:       ', continuity
        print '(a,es24.16)', 'advection pair defect:             ', &
            advection + continuity
        print '(a,es24.16)', 'pressure work in kinetic energy:   ', pressure_work
        print '(a,es24.16)', 'pressure work in internal energy:  ', internal
        print '(a,es24.16)', 'pressure pair defect:              ', &
            pressure_work + internal
        print '(a,es24.16)', 'd(kinetic + internal)/dt:          ', kinetic + internal
        print '(a,es24.16)', 'explicit source and boundary power:', 0.0_dp
        print '(a,es24.16)', 'relative total energy defect:      ', relative_defect
        print '(a,es24.16)', 'relative mass solve residual:      ', solve_defect
        print '(a,es24.16)', 'relative independent ledger defect:', ledger_defect
        print '(a)', 'B=J=0; eta=nu=chi=D=Srho=Qe=Qi=V=0; '// &
            'no magnetic or dissipative power'
    end subroutine
    subroutine full_ledger(kind, convention, relative_defect)
        integer, intent(in) :: kind, convention
        real(dp), intent(out) :: relative_defect
        type(vrmhd_section_t) :: op
        type(vrmhd_parameters_t) :: p
        type(mesh_2d_t) :: mesh
        type(meridian_mesh_t) :: mm
        type(gs_symmetry_t) :: symmetry
        real(dp), allocatable :: u(:), ud(:), zero(:), r(:), unused(:)
        real(dp), allocatable :: mass(:, :), original(:, :), rhs(:, :), mj(:, :)
        real(dp), allocatable :: source(:, :, :), js(:, :, :), basis(:, :)
        real(dp), allocatable :: harmonic(:), load(:), wall(:), recovery(:, :)
        integer, allocatable :: free_index(:), pivot(:), current_index(:)
        logical, allocatable :: active(:)
        real(dp) :: wall_field(3), wall_terms(4), wall_power
        real(dp) :: x(2), weight, v(36), vd(36), gv(3, 3), strain(3, 3), tau(3, 3)
        real(dp) :: b(3), projector(3, 3), flux(3), identity(3, 3), pressure, divv
        real(dp) :: kinetic, internal, magnetic, ampere, ideal, lorentz, viscous
        real(dp) :: pressure_work, ohmic, current_drive, voltage_drive, heat_source
        real(dp) :: density_port, source_port, diffusion_port, conduction, chi, g(3)
        real(dp) :: scale, solve_defect, ledger_defect, shift, voltage, h(3), eta, nu
        real(dp) :: skew(2), reaction_power, magnetic_balance, kinetic_balance
        real(dp) :: internal_balance, energy_rate, external_power, removed
        integer :: status, n0, n1, nf, j, i, q, t, info, a0, j0, l1, s
        interface
            subroutine dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
                import :: dp
                integer, intent(in) :: n, nrhs, lda, ldb
                integer, intent(out) :: ipiv(*), info
                real(dp), intent(inout) :: a(lda, *), b(ldb, *)
            end subroutine
        end interface

        symmetry = gs_symmetry_t(kind, GS_CHART_CARTESIAN, 0.0_dp)
        shift = 0
        if (kind == GS_ROTATION) then
            shift = 2
            symmetry%chart = GS_CHART_RZ
        end if
        if (kind == GS_SCREW) symmetry%k = 0.7_dp
        call mesh%create_rectangular(4, 4, shift - 0.4_dp, shift + 0.6_dp, &
            -0.3_dp, 0.7_dp)
        call meridian_mesh_straight(mesh, 1, mm, status)
        call require(status == 0, 'magnetic energy mesh')
        if (family_degree > 0) then
            call vrmhd_section_build(mm, symmetry, family_degree + 1, 1, op, status, &
                scalar_degree=family_degree, velocity_degree=family_degree + 1)
        else
            call vrmhd_section_build(mm, symmetry, 1, 1, op, status)
        end if
        call require(status == 0, 'magnetic energy section')
        n0 = op%scalar_fs%n0; n1 = op%fs%sp(1)%ndof
        l1 = op%fs%sp(1)%nloc
        a0 = op%offset(7); j0 = a0 + n1
        p%mi = 1.2_dp; p%mu0 = 0.8_dp; p%epsilon_b = 0.02_dp
        p%eta0 = 0.015_dp; p%nu0 = 0.02_dp; p%nu_parallel = 0.01_dp
        p%chi_parallel0 = [0.04_dp, 0.03_dp]
        p%chi_perpendicular = [0.02_dp, 0.01_dp]
        p%d_parallel = 0.01_dp; p%d_perpendicular = 0.005_dp
        p%equilibration_rate = 0.02_dp
        if (convention == VRMHD_LOOP_VOLTAGE_SOURCE) then
            ! Non-classical closures and a partial Ohmic deposit.
            p%eta_exponent = -1; p%nu_exponent = 0.5_dp
            p%chi_parallel_exponent = 1.5_dp; p%ohmic_fraction = 0.25_dp
        end if
        allocate (u(op%ndof), ud(op%ndof), zero(op%ndof), r(op%ndof), &
            unused(op%ndof), active(op%ndof), &
            source(6, op%fs%nq, op%fs%mm%topology%n_triangles), &
            js(3, op%fs%nq, op%fs%mm%topology%n_triangles), &
            basis(36, op%nlocal), harmonic(n1), load(op%ndof), mj(n1, n1), wall(n1))
        u = 0; ud = 0; zero = 0; source = 0
        call random_number(u(:j0))
        u(:n0) = 1 + 0.4_dp*u(:n0)
        u(n0 + 1:op%offset(5)) = 0.4_dp*(u(n0 + 1:op%offset(5)) - 0.5_dp)
        u(op%offset(5) + 1:op%offset(7)) = 1 + 0.2_dp*u(op%offset(5) + 1:op%offset(7))
        u(a0 + 1:j0) = 0.03_dp*(u(a0 + 1:j0) - 0.5_dp)
        where (op%fixed(n0 + 1:op%offset(5))) u(n0 + 1:op%offset(5)) = 0
        ! Compact support isolates the prescribed current trace from magnetic trace data.
        do t = 1, op%fs%mm%topology%n_triangles
            if (.not. any(op%scalar_fs%sp(0)%bnd(op%scalar_fs%sp(0)%dof(:, t)))) cycle
            u(op%map(op%local_offset(7) + 1:op%local_offset(8), t)) = 0
        end do
        call require(maxval(abs(u(a0 + 1:j0))) > 0, 'nonzero random potential')
        call random_number(js); js = 0.02_dp*(js - 0.5_dp)
        do t = 1, op%fs%mm%topology%n_triangles
            if (any(op%scalar_fs%sp(0)%bnd(op%scalar_fs%sp(0)%dof(:, t)))) js(:, :, t) = 0
        end do
        where (op%scalar_fs%sp(0)%bnd) u(op%offset(5) + 1:op%offset(6)) = 1
        source(4, :, :) = 0.003_dp
        source(5, :, :) = 0.001_dp; source(6, :, :) = -0.0003_dp
        voltage = 0
        if (convention /= VRMHD_CURRENT_SOURCE) voltage = 0.004_dp
        call vrmhd_drive_build(op, convention, 1.7_dp, voltage, js, harmonic, &
            source, load, status)
        call require(status == 0, 'energy drive')

        call vrmhd_section_residual(op, p, u, zero, zero, zero, source, &
            r, unused, status)
        call require(status == 0, 'Ampere recovery residual')
        mj = 0
        do t = 1, op%fs%mm%topology%n_triangles
            do q = 1, op%fs%nq
                call vrmhd_section_basis(op, t, q, x, weight, basis, status)
                call require(status == 0, 'Ampere recovery quadrature')
                do j = op%local_offset(8) + 1, op%nlocal
                    do i = op%local_offset(8) + 1, op%nlocal
                        mj(op%map(i, t) - j0, op%map(j, t) - j0) = &
                            mj(op%map(i, t) - j0, op%map(j, t) - j0) &
                            + weight*p%mu0*dot_product(basis(8:10, i), basis(8:10, j))
                    end do
                end do
            end do
        end do
        allocate (rhs(n1, 1))
        if (convention /= VRMHD_CURRENT_SOURCE) then
            call current_wall_dofs(op, voltage_current, &
                [real(kind, dp), voltage/p%eta0], rhs(:, 1))
        else
            call random_number(rhs)
            rhs = 0.02_dp*(rhs - 0.5_dp)
        end if
        call vrmhd_current_lift(op, rhs(:, 1), wall)
        call require(maxval(abs(wall)) > 1.0e-5_dp, 'prescribed wall current active')
        u(j0 + 1:) = wall
        if (convention == VRMHD_CURRENT_SOURCE) then
            do t = 1, op%fs%mm%topology%n_triangles
                do q = 1, op%fs%nq
                    call vrmhd_section_basis(op, t, q, x, weight, basis, status)
                    call require(status == 0, 'wall source basis')
                    source(:3, q, t) = source(:3, q, t) &
                        + matmul(basis(8:10, op%local_offset(8) + 1:), &
                        wall(op%map(op%local_offset(8) + 1:, t) - j0))
                end do
            end do
        end if
        deallocate (rhs)
        nf = count(.not. op%fs%sp(1)%bnd)
        allocate (current_index(nf), rhs(nf, 1), pivot(nf), recovery(nf, nf))
        current_index = pack([(j, j=1, n1)], .not. op%fs%sp(1)%bnd)
        call vrmhd_section_residual(op, p, u, zero, zero, zero, source, &
            r, unused, status)
        call require(status == 0, 'lifted Ampere residual')
        recovery = mj(current_index, current_index)
        rhs(:, 1) = -r(j0 + current_index)
        call dgesv(nf, 1, recovery, nf, pivot, rhs, nf, info)
        call require(info == 0, 'H0(curl) weak Ampere current solve')
        u(j0 + current_index) = rhs(:, 1)
        ! Refine against the residual action, whose summation differs from the matrix.
        do i = 1, 2
            call vrmhd_section_residual(op, p, u, zero, zero, zero, source, &
                r, unused, status)
            call require(status == 0, 'Ampere refinement residual')
            recovery = mj(current_index, current_index)
            rhs(:, 1) = -r(j0 + current_index)
            call dgesv(nf, 1, recovery, nf, pivot, rhs, nf, info)
            call require(info == 0, 'Ampere refinement solve')
            u(j0 + current_index) = u(j0 + current_index) + rhs(:, 1)
        end do
        deallocate (rhs, pivot, recovery)

        active = .true.; active(j0 + 1:) = .false.
        active(n0 + 1:op%offset(5)) = .not. op%fixed(n0 + 1:op%offset(5))
        active(a0 + 1:j0) = .not. op%fixed(a0 + 1:j0)
        nf = count(active)
        allocate (free_index(nf), pivot(nf), mass(nf, nf), original(nf, nf), rhs(nf, 1))
        free_index = pack([(j, j=1, op%ndof)], active)
        call vrmhd_section_residual(op, p, u, zero, zero, zero, source, &
            r, unused, status)
        call require(status == 0, 'full energy spatial residual')
        r = r - load
        call require(maxval(abs(r(j0 + current_index))) < 1.0e-13_dp, 'weak Ampere residual')
        rhs(:, 1) = -r(free_index)
        call free_mass(op, p, u, free_index, mass, status)
        call require(status == 0, 'full energy mass assembly')
        original = mass
        call dgesv(nf, 1, mass, nf, pivot, rhs, nf, info)
        call require(info == 0, 'full energy mass solve')
        ud(free_index) = rhs(:, 1)
        solve_defect = maxval(abs(matmul(original, rhs(:, 1)) + r(free_index))) &
            /maxval(abs(r(free_index)))
        call require(solve_defect < 1.0e-13_dp, 'full energy mass solve accuracy')
        call vrmhd_section_residual(op, p, u, ud, zero, zero, source, &
            r, unused, status)
        call require(status == 0, 'full energy rate residual')
        r = r - load
        reaction_power = dot_product(wall, r(a0 + 1:j0))

        kinetic = 0; internal = 0; magnetic = 0; ampere = 0
        pressure_work = 0; lorentz = 0; ideal = 0; viscous = 0
        ohmic = 0; current_drive = 0; voltage_drive = 0; heat_source = 0
        source_port = 0; diffusion_port = 0; conduction = 0; skew = 0; wall_terms = 0
        identity = 0
        do i = 1, 3
            identity(i, i) = 1
        end do
        do t = 1, op%fs%mm%topology%n_triangles
            do q = 1, op%fs%nq
                call vrmhd_section_basis(op, t, q, x, weight, basis, status)
                call require(status == 0, 'full energy quadrature')
                v = matmul(basis, u(op%map(:, t)))
                vd = matmul(basis, ud(op%map(:, t)))
                gv = reshape(v(16:24), [3, 3])
                divv = gv(1, 1) + gv(2, 2) + gv(3, 3)
                pressure = v(1)*(v(11) + v(12))/p%mi
                b = v(25:27)/sqrt(sum(v(25:27)**2) + p%epsilon_b**2)
                strain = (gv + transpose(gv))/2 - identity*divv/3
                projector = spread(b, 2, 3)*spread(b, 1, 3) - identity*sum(b*b)/3
                eta = p%eta0*(v(11)/p%t0)**p%eta_exponent
                nu = p%nu0*(v(12)/p%t0)**p%nu_exponent
                tau = 2*nu*strain + 3*p%nu_parallel*sum(projector*strain)*projector
                flux = p%d_perpendicular*v(13:15) &
                    + (p%d_parallel - p%d_perpendicular)*b*dot_product(b, v(13:15))
                h = [0.0_dp, 0.0_dp, 1/1.7_dp]
                if (kind == GS_ROTATION) h = [0.0_dp, 1/(2*acos(-1.0_dp)*x(1)), 0.0_dp]
                wall_field = matmul(basis(8:10, op%local_offset(8) + 1:), &
                    wall(op%map(op%local_offset(8) + 1:, t) - j0))
                wall_terms = wall_terms + weight*[ &
                    dot_product(wall_field, vd(5:7)), &
                    -dot_product(wall_field, cross(v(2:4), v(25:27))), &
                    eta*dot_product(wall_field, v(8:10) - source(:3, q, t)), &
                    -voltage*dot_product(wall_field, h)]
                kinetic = kinetic + weight*(vd(1)*sum(v(2:4)**2)/2 &
                    + v(1)*dot_product(v(2:4), vd(2:4)))
                internal = internal + weight*1.5_dp/p%mi*( &
                    vd(1)*(v(11) + v(12)) + v(1)*(vd(11) + vd(12)))
                magnetic = magnetic + weight*dot_product(v(25:27), vd(25:27))/p%mu0
                ampere = ampere + weight*dot_product(v(8:10), vd(5:7))
                pressure_work = pressure_work + weight*pressure*divv
                lorentz = lorentz + weight*dot_product(v(2:4), cross(v(8:10), v(25:27)))
                ideal = ideal + weight*dot_product(v(8:10), cross(v(2:4), v(25:27)))
                viscous = viscous + weight*sum(tau*gv)
                ohmic = ohmic + weight*eta*sum(v(8:10)**2)
                current_drive = current_drive &
                    + weight*eta*dot_product(v(8:10), source(:3, q, t))
                voltage_drive = voltage_drive + weight*voltage*dot_product(h, v(8:10))
                heat_source = heat_source + weight*sum(source(5:6, q, t))
                source_port = source_port - weight*source(4, q, t)*sum(v(2:4)**2)/2
                diffusion_port = diffusion_port &
                    + weight*dot_product(flux, matmul(transpose(gv), v(2:4)))
                skew(1) = skew(1) &
                    - weight*v(1)*dot_product(v(2:4), matmul(gv, v(2:4)))/2
                skew(2) = skew(2) &
                    + weight*v(1)*dot_product(v(2:4), matmul(gv, v(2:4)))/2
                do s = 1, 2
                    g = v(25 + 3*s:27 + 3*s)
                    chi = p%chi_parallel0(s)*(v(10 + s)/p%t0)**p%chi_parallel_exponent
                    conduction = conduction + weight*v(1)/p%mi &
                        *(p%chi_perpendicular(s)*sum(g*g) &
                        + (chi - p%chi_perpendicular(s))*dot_product(b, g)**2)
                end do
            end do
        end do
        wall_power = sum(wall_terms)
        density_port = source_port + diffusion_port
        kinetic_balance = kinetic - pressure_work - lorentz + viscous - density_port
        removed = (1 - p%ohmic_fraction)*(ohmic - current_drive)
        internal_balance = internal + pressure_work - viscous &
            - p%ohmic_fraction*(ohmic - current_drive) - heat_source
        magnetic_balance = magnetic - ideal + ohmic - current_drive - voltage_drive
        scale = 2*(abs(pressure_work) + abs(lorentz) + abs(viscous) + abs(ohmic)) &
            + sum(abs(skew)) + abs(source_port) + abs(diffusion_port) &
            + 2*abs(current_drive) + abs(voltage_drive) + abs(heat_source) &
            + sum(abs(wall_terms))
        call require(min(viscous, ohmic, conduction) > 1.0e-8_dp, &
            'all dissipative pairs active')
        call require(abs(lorentz) > 1.0e-10_dp, 'magnetic exchange active')
        ledger_defect = max(abs(kinetic_balance), abs(internal_balance), &
            abs(magnetic - ampere), abs(lorentz + ideal), abs(sum(skew)), &
            abs(magnetic_balance - wall_power), abs(wall_power - reaction_power))/scale
        call require(ledger_defect < 1.0e-12_dp, 'independent full energy accounting')
        energy_rate = kinetic + internal + magnetic
        call require(abs(wall_power) > 1.0e-10_dp, 'nonzero independent wall power')
        if (convention /= VRMHD_CURRENT_SOURCE) &
            call require(abs(voltage_drive) > 1.0e-10_dp, 'nonzero loop power')
        if (convention == VRMHD_LOOP_VOLTAGE_SOURCE) then
            call require(abs(current_drive) > 1.0e-10_dp, 'nonzero source power')
            call require(abs(removed) > 1.0e-10_dp, 'nonzero removed Ohmic power')
        end if
        external_power = density_port + heat_source + voltage_drive + wall_power &
            - removed
        relative_defect = abs(energy_rate - external_power)/scale
        print '(a,2i3)', 'full ledger symmetry, drive (1 current, 2 voltage, 3 both): ', &
            kind, convention
        if (kind == GS_ROTATION) then
            print '(a)', 'power in W/rad; prescribed current wall power included'
        else
            print '(a)', 'power in W/m; prescribed current wall power included'
        end if
        print '(a,3es24.16)', 'skew forward / reverse / defect: ', skew, sum(skew)
        print '(a,3es24.16)', 'pressure kinetic / internal / defect: ', &
            pressure_work, internal - viscous - ohmic + current_drive - heat_source, &
            internal_balance
        print '(a,3es24.16)', 'Lorentz / ideal induction / defect: ', &
            lorentz, ideal, lorentz + ideal
        print '(a,3es24.16)', 'magnetic derivative / Ampere pairing / defect: ', &
            magnetic, ampere, magnetic - ampere
        print '(a,2es24.16)', 'Ohmic magnetic loss / electron heating: ', -ohmic, ohmic
        print '(a,2es24.16)', 'current-source magnetic / thermal withdrawal: ', &
            current_drive, -current_drive
        print '(a,2es24.16)', 'viscous kinetic loss / ion heating: ', -viscous, viscous
        print '(a,2es24.16)', 'heat conduction activity / net boundary flux: ', &
            conduction, 0.0_dp
        print '(a,3es24.16)', 'mass source / diffusion / total kinetic port: ', &
            source_port, diffusion_port, density_port
        print '(a,es24.16)', 'Ohmic power removed from the electrons: ', removed
        print '(a,2es24.16)', 'loop-voltage / external heat power: ', &
            voltage_drive, heat_source
        print '(a,4es24.16)', 'wall rate / ideal / resistive / voltage power: ', wall_terms
        print '(a,es24.16)', 'explicit current wall power: ', wall_power
        print '(a,3es24.16)', 'kinetic / internal / magnetic rate: ', &
            kinetic, internal, magnetic
        print '(a,3es24.16)', 'energy rate / explicit power / defect: ', &
            energy_rate, external_power, energy_rate - external_power
        print '(a,2es24.16)', 'magnetic wall work / independent row pairing: ', &
            magnetic_balance, reaction_power
        print '(a,es24.16)', 'maximum current coefficient on fixed A trace: ', &
            maxval(abs(u(j0 + 1:)), mask=op%fixed(a0 + 1:j0))
        print '(a,es24.16)', 'applied loop voltage in V: ', voltage
        print '(a,3es24.16)', 'relative energy / solve / accounting defects: ', &
            relative_defect, solve_defect, ledger_defect
    end subroutine

    subroutine free_mass(op, p, u, free_index, mass, status)
        !! Dense time mass on the free coefficients from one sparse assembly.
        type(vrmhd_section_t), intent(in) :: op
        type(vrmhd_parameters_t), intent(in) :: p
        real(dp), intent(in) :: u(:)
        integer, intent(in) :: free_index(:)
        real(dp), intent(out) :: mass(:, :)
        integer, intent(out) :: status
        type(csc_t) :: time_mass
        real(dp), allocatable :: zero(:), source(:, :, :)
        integer, allocatable :: position(:)
        integer :: j, k
        allocate (zero(op%ndof), position(op%ndof), &
            source(6, op%fs%nq, op%fs%mm%topology%n_triangles))
        zero = 0; source = 0; mass = 0
        call vrmhd_section_jacobian(op, p, u, zero, source, 1.0_dp, time_mass, &
            status, mass_only=.true.)
        if (status /= 0) return
        position = 0
        position(free_index) = [(j, j=1, size(free_index))]
        do j = 1, op%ndof
            if (position(j) == 0) cycle
            do k = time_mass%col_ptr(j), time_mass%col_ptr(j + 1) - 1
                if (position(time_mass%row_idx(k)) == 0) cycle
                mass(position(time_mass%row_idx(k)), position(j)) = time_mass%val(k)
            end do
        end do
    end subroutine

    subroutine voltage_current(x, data, current)
        real(dp), intent(in) :: x(2), data(:)
        real(dp), intent(out) :: current(3)
        current = [0.0_dp, 0.0_dp, data(2)/1.7_dp]
        if (nint(data(1)) == GS_ROTATION) &
            current = [0.0_dp, data(2)/(2*acos(-1.0_dp)*x(1)), 0.0_dp]
    end subroutine

    function cross(a, b) result(c)
        real(dp), intent(in) :: a(3), b(3)
        real(dp) :: c(3)
        c = [a(2)*b(3) - a(3)*b(2), a(3)*b(1) - a(1)*b(3), &
            a(1)*b(2) - a(2)*b(1)]
    end function
end program test_vrmhd_energy
