program test_vrmhd_energy
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_gs_symmetry, only: gs_symmetry_t, GS_ROTATION, GS_SCREW, &
        GS_CHART_RZ, GS_CHART_CARTESIAN
    use kin6d_vrmhd_parameters, only: vrmhd_parameters_t
    use kin6d_vrmhd_section, only: vrmhd_section_t, vrmhd_section_build, &
        vrmhd_section_basis, vrmhd_section_residual, vrmhd_section_mass
    use fortfem_mesh_2d, only: mesh_2d_t
    use fortfem_meridian_curved_mesh, only: meridian_mesh_t, meridian_mesh_straight
    implicit none
    integer :: kind, seed_size, i
    integer, allocatable :: seed(:)
    logical :: passed
    real(dp) :: defect

    call random_seed(size=seed_size)
    allocate (seed(seed_size))
    seed = [(104729 + 37*i, i=1, seed_size)]
    call random_seed(put=seed)
    passed = .true.
    do kind = GS_ROTATION, GS_SCREW
        call hydrodynamic_ledger(kind, defect)
        passed = passed .and. defect <= 1.0e-12_dp
    end do
    call require(passed, 'closed semidiscrete energy ledger exceeds 1e-12')
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
        real(dp), allocatable :: u(:), ud(:), zero(:), r(:), unused(:), probe(:)
        real(dp), allocatable :: mass(:, :), original(:, :), rhs(:, :)
        real(dp), allocatable :: source(:, :, :), basis(:, :)
        integer, allocatable :: free_index(:), pivot(:)
        logical, allocatable :: active(:)
        real(dp) :: x(2), weight, v(36), vd(36), gv(3, 3), pressure
        real(dp) :: kinetic, internal, advection, continuity, pressure_work
        real(dp) :: strong_product, scale, solve_defect, ledger_defect, shift
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
        call vrmhd_section_build(mm, symmetry, 1, 1, op, status)
        call require(status == 0, 'energy section')
        n0 = op%fs%n0
        p%mi = 1.2_dp; p%mu0 = 0.8_dp
        allocate (u(op%ndof), ud(op%ndof), zero(op%ndof), r(op%ndof), &
            unused(op%ndof), probe(op%ndof), active(op%ndof), &
            source(6, op%fs%nq, op%fs%mm%topology%n_triangles), &
            basis(36, op%nlocal))
        u = 0; ud = 0; zero = 0; source = 0
        call random_number(u(:6*n0))
        u(:n0) = 1 + 0.4_dp*u(:n0)
        u(n0 + 1:4*n0) = 0.4_dp*(u(n0 + 1:4*n0) - 0.5_dp)
        u(4*n0 + 1:6*n0) = 1 + 0.2_dp*u(4*n0 + 1:6*n0)
        where (op%fixed(n0 + 1:4*n0)) u(n0 + 1:4*n0) = 0

        ! No-flux scalars retain constants; only velocity has wall constraints.
        active = .false.
        active(:6*n0) = .true.
        active(n0 + 1:4*n0) = .not. op%fixed(n0 + 1:4*n0)
        nf = count(active)
        allocate (free_index(nf), pivot(nf), mass(nf, nf), original(nf, nf), &
            rhs(nf, 1))
        free_index = pack([(j, j=1, op%ndof)], active)
        call vrmhd_section_residual(op, p, u, zero, zero, zero, source, &
            r, unused, status)
        call require(status == 0, 'energy spatial residual')
        rhs(:, 1) = -r(free_index)
        do j = 1, nf
            probe = 0; probe(free_index(j)) = 1
            call vrmhd_section_mass(op, p, u, probe, unused, status)
            call require(status == 0, 'energy mass assembly')
            mass(:, j) = unused(free_index)
        end do
        original = mass
        call dgesv(nf, 1, mass, nf, pivot, rhs, nf, info)
        call require(info == 0, 'energy mass solve')
        ud(free_index) = rhs(:, 1)
        solve_defect = maxval(abs(matmul(original, rhs(:, 1)) + r(free_index))) &
            /maxval(abs(r(free_index)))
        call require(solve_defect < 1.0e-13_dp, 'energy mass solve accuracy')

        kinetic = 0; internal = 0; advection = 0; continuity = 0
        pressure_work = 0; strong_product = 0
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
                advection = advection + weight*v(1)* &
                    dot_product(v(2:4), matmul(gv, v(2:4)))
                continuity = continuity - weight*0.5_dp*vd(1)*sum(v(2:4)**2)
                pressure_work = pressure_work + weight*pressure* &
                    (gv(1, 1) + gv(2, 2) + gv(3, 3))
                strong_product = strong_product + weight*0.5_dp*sum(v(2:4)**2)* &
                    (vd(1) + dot_product(v(2:4), v(13:15)) &
                    + v(1)*(gv(1, 1) + gv(2, 2) + gv(3, 3)))
            end do
        end do
        scale = abs(advection) + abs(continuity) + 2*abs(pressure_work)
        call require(scale > tiny(scale), 'nontrivial energy exchanges')
        ledger_defect = max(abs(kinetic - advection - continuity - pressure_work), &
            abs(internal + pressure_work), &
            abs(strong_product + advection + continuity))/scale
        call require(ledger_defect < 1.0e-12_dp, 'independent energy accounting')
        relative_defect = abs(kinetic + internal)/scale
        print '(a,i2)', 'symmetry (1 rotation, 2 translation, 3 screw): ', kind
        if (kind == GS_ROTATION) then
            print '(a)', 'power in W/rad (per unit rotation angle)'
        else
            print '(a)', 'power in W/m (per unit axial orbit coordinate)'
        end if
        print '(a,es24.16)', 'momentum advection work:           ', advection
        print '(a,es24.16)', 'kinetic continuity partner:        ', continuity
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
end program test_vrmhd_energy
