program test_vrmhd_stationary
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_gs_symmetry, only: gs_symmetry_t, GS_ROTATION, GS_SCREW, &
        GS_CHART_RZ, GS_CHART_CARTESIAN
    use kin6d_vrmhd_parameters, only: vrmhd_parameters_t
    use kin6d_vrmhd_section
    use kin6d_vrmhd_stationary
    use kin6d_vrmhd_boundary, only: vrmhd_ampere_boundary
    use vrmhd_mms_leaves, only: mms_rotation, mms_screw
    use fortfem_mesh_2d, only: mesh_2d_t
    use fortfem_meridian_curved_mesh, only: meridian_mesh_t, meridian_mesh_straight
    use fortsparse, only: csc_t, csc_from_triplet, csc_matvec, sparse_solve_once, &
        fortsparse_status_t, status_ok
    implicit none
    character(32) :: gate
    type(gs_symmetry_t) :: trace_symmetry
    integer :: kind, level
    real(dp) :: errors(8, 3), rates(8)
    call get_command_argument(1, gate)
    do kind = 1, 3
        select case (trim(gate))
        case ('boundary')
            call check_boundary(kind)
        case ('derivatives')
            call check_derivatives(kind)
        case ('solve')
            do level = 1, 3
                call solve_manufactured(kind, 2**(level + 1), errors(:, level))
                print '(a,2i3,8es13.5)', 'solved natural-norm errors ', &
                    kind, level, errors(:, level)
                if (level == 1) cycle
                rates = log(errors(:, level - 1)/errors(:, level))/log(2.0_dp)
                print '(a,2i3,8f8.3)', 'solved natural-norm rates ', kind, level, rates
                call require(minval(rates) > 0.8_dp, 'P1 solved design order')
            end do
        case default
            error stop 'unknown stationary gate'
        end select
    end do
contains
    subroutine require(ok, message)
        logical, intent(in) :: ok
        character(*), intent(in) :: message
        if (ok) return
        print '(a)', message
        error stop 'stationary MHD gate'
    end subroutine

    subroutine exact_at(op, x, exact)
        type(vrmhd_section_t), intent(in) :: op
        real(dp), intent(in) :: x(2)
        real(dp), intent(out) :: exact(48)
        if (op%symmetry%kind == GS_ROTATION) then
            call mms_rotation(x(1), x(2), op%symmetry%k, exact)
        else
            call mms_screw(x(1), x(2), op%symmetry%k, exact)
        end if
    end subroutine

    subroutine magnetic_trace(x, b)
        real(dp), intent(in) :: x(2)
        real(dp), intent(out) :: b(3)
        real(dp) :: exact(48)
        if (trace_symmetry%kind == GS_ROTATION) then
            call mms_rotation(x(1), x(2), trace_symmetry%k, exact)
        else
            call mms_screw(x(1), x(2), trace_symmetry%k, exact)
        end if
        b = exact(25:27)
    end subroutine

    pure function norms(u) result(e)
        real(dp), intent(in) :: u(36)
        real(dp) :: e(8)
        integer :: i
        e(1) = u(1)**2 + sum(u(13:15)**2)
        do i = 2, 4
            e(i) = u(i)**2 + u(i + 14)**2 + u(i + 17)**2 + u(i + 20)**2
        end do
        e(5) = sum(u(5:7)**2) + sum(u(25:27)**2)
        e(6) = sum(u(8:10)**2) + sum(u(34:36)**2)
        e(7) = u(11)**2 + sum(u(28:30)**2)
        e(8) = u(12)**2 + sum(u(31:33)**2)
    end function

    subroutine prepare(kind, n, op, steady, p, u, load)
        integer, intent(in) :: kind, n
        type(vrmhd_section_t), intent(out) :: op
        type(vrmhd_stationary_t), intent(out) :: steady
        type(vrmhd_parameters_t), intent(out) :: p
        real(dp), allocatable, intent(out) :: u(:), load(:)
        type(mesh_2d_t) :: mesh
        type(meridian_mesh_t) :: mm
        type(gs_symmetry_t) :: symmetry
        real(dp) :: shift
        integer :: status
        symmetry = gs_symmetry_t(kind, GS_CHART_CARTESIAN, 0.0_dp)
        shift = 0
        if (kind == GS_ROTATION) then
            shift = 2; symmetry%chart = GS_CHART_RZ
        end if
        if (kind == GS_SCREW) symmetry%k = 0.7_dp
        call mesh%create_rectangular(n + 1, n + 1, shift + 0.1_dp, shift + 1.1_dp, &
            0.2_dp, 1.2_dp)
        call meridian_mesh_straight(mesh, 1, mm, status)
        call require(status == 0, 'mesh')
        call vrmhd_section_build(mm, symmetry, 1, 1, op, status)
        call require(status == 0, 'section')
        call vrmhd_stationary_build(op, steady, status)
        call require(status == 0, 'Coulomb gauge')
        p%mi = 1; p%mu0 = 1; p%eta0 = 0.2_dp; p%nu0 = 0.2_dp
        p%chi_parallel0 = 0.2_dp; p%chi_perpendicular = 0.05_dp
        p%d_parallel = 0.1_dp; p%d_perpendicular = 0.1_dp; p%epsilon_b = 0.1_dp
        allocate (u(steady%ndof), load(op%ndof))
        u = 0
        call project(op, u(:op%ndof), load)
    end subroutine

    subroutine project(op, u, load)
        type(vrmhd_section_t), intent(in) :: op
        real(dp), intent(out) :: u(op%ndof), load(op%ndof)
        type(csc_t) :: mass
        type(fortsparse_status_t) :: st
        integer, allocatable :: rows(:), cols(:)
        real(dp), allocatable :: values(:), rhs(:)
        real(dp) :: basis(36, op%nlocal), block(op%nlocal, op%nlocal)
        real(dp) :: exact(48), x(2), weight
        integer :: i, j, t, q, pos, status, nt
        nt = op%fs%mm%topology%n_triangles
        allocate (rows(nt*op%nlocal**2), cols(nt*op%nlocal**2), &
            values(nt*op%nlocal**2), rhs(op%ndof))
        rhs = 0; load = 0
        do t = 1, nt
            block = 0
                do q = 1, op%fs%nq
                    call vrmhd_section_basis(op, t, q, x, weight, basis, status)
                    call require(status == 0, 'projection basis')
                    call exact_at(op, x, exact)
                    block = block + weight*matmul(transpose(basis), basis)
                        do i = 1, op%nlocal
                            j = op%map(i, t)
                            rhs(j) = rhs(j) + weight*dot_product(basis(:, i), exact(:36))
                            load(j) = load(j) + weight*dot_product(basis(:12, i), exact(37:48))
                        end do
                    end do
                    pos = (t - 1)*op%nlocal**2
                    do j = 1, op%nlocal
                        do i = 1, op%nlocal
                            pos = pos + 1
                            rows(pos) = op%map(i, t); cols(pos) = op%map(j, t)
                            values(pos) = block(i, j)
                        end do
                    end do
                end do
                call csc_from_triplet(op%ndof, op%ndof, rows, cols, values, mass, st)
                call require(status_ok(st), 'projection assembly')
                call sparse_solve_once(mass, rhs, u, st)
                call require(status_ok(st), 'projection solve')
            end subroutine

            subroutine check_boundary(kind)
                integer, intent(in) :: kind
                type(vrmhd_section_t) :: op
                type(vrmhd_stationary_t) :: steady
                type(vrmhd_parameters_t) :: p
                real(dp), allocatable :: u(:), load(:), surface(:), volume(:)
                real(dp), allocatable :: basis(:, :)
                real(dp) :: exact(48), x(2), weight, value
                real(dp) :: defect
                integer :: t, q, i, row, status
                call prepare(kind, 2, op, steady, p, u, load)
                allocate (surface(op%ndof), volume(op%ndof), basis(36, op%nlocal))
                trace_symmetry = op%symmetry
                call vrmhd_ampere_boundary(op, magnetic_trace, surface, status)
                call require(status == 0, 'Ampere trace assembly')
                volume = 0
                do t = 1, op%fs%mm%topology%n_triangles
                    do q = 1, op%fs%nq
                        call vrmhd_section_basis(op, t, q, x, weight, basis, status)
                        call exact_at(op, x, exact)
                        do i = 1, op%nlocal
                            if (op%row(i) /= 4) cycle
                            value = dot_product(exact(8:10), basis(8:10, i)) &
                                - dot_product(exact(25:27), basis(34:36, i))
                            row = op%map(i, t)
                            volume(row) = volume(row) + weight*value
                        end do
                    end do
                end do
                defect = maxval(abs(volume + surface))
                print '(a,i2,2es14.6)', 'Ampere volume+surface, surface norm ', &
                    kind, defect, maxval(abs(surface))
                call require(maxval(abs(surface)) > 1.0e-3_dp, 'nonzero magnetic trace')
                call require(defect < 2.0e-8_dp, 'independent Ampere integration by parts')
            end subroutine

            subroutine check_derivatives(kind)
                integer, intent(in) :: kind
                type(vrmhd_section_t) :: op
                type(vrmhd_stationary_t) :: steady
                type(vrmhd_parameters_t) :: p
                type(csc_t) :: jac
                real(dp), allocatable :: u(:), load(:), source(:, :, :), du(:), zero(:)
                real(dp), allocatable :: r(:), jv(:), plus(:), minus(:), unused(:), boundary(:)
                real(dp) :: eps, error, previous
                integer :: i, status
                call prepare(kind, 2, op, steady, p, u, load)
                allocate (source(6, op%fs%nq, op%fs%mm%topology%n_triangles))
                allocate (du(steady%ndof), zero(steady%ndof), r(steady%ndof), &
                    jv(steady%ndof), plus(steady%ndof), minus(steady%ndof), &
                    unused(steady%ndof), boundary(op%ndof))
                source = 0; zero = 0; boundary = u(:op%ndof)
                do i = 1, steady%ndof
                    du(i) = sin(real(i, dp))/10
                end do
                u(op%ndof + 1:) = 0.03_dp
                call vrmhd_stationary_residual(op, steady, p, u, du, boundary, source, &
                    load, r, jv, status)
                call require(status == 0, 'gauged residual')
                call vrmhd_stationary_jacobian(op, steady, p, u, boundary, source, jac, status)
                call require(status == 0, 'gauged sparse Jacobian')
                error = maxval(abs(csc_matvec(jac, du) - jv))
                call require(error < 1.0e-11_dp, 'sparse and matrix-free gauged Jacobian')
                previous = huge(1.0_dp)
                do i = 1, 3
                    eps = 0.02_dp/2**i
                    call vrmhd_stationary_residual(op, steady, p, u + eps*du, zero, &
                        boundary, source, load, plus, unused, status)
                    call vrmhd_stationary_residual(op, steady, p, u - eps*du, zero, &
                        boundary, source, load, minus, unused, status)
                    error = maxval(abs((plus - minus)/(2*eps) - jv))
                    print '(a,2i3,es14.6)', 'gauged centered difference ', kind, i, error
                    call require(error < 0.3_dp*previous, 'quadratic centered differences')
                    previous = error
                end do
            end subroutine

            subroutine solve_manufactured(kind, n, error)
                integer, intent(in) :: kind, n
                real(dp), intent(out) :: error(8)
                type(vrmhd_section_t) :: op
                type(vrmhd_stationary_t) :: steady
                type(vrmhd_parameters_t) :: p
                type(csc_t) :: jac
                type(fortsparse_status_t) :: st
                real(dp), allocatable :: u(:), load(:), source(:, :, :), boundary(:)
                real(dp), allocatable :: r(:), unused(:), zero(:), step(:), surface(:)
                real(dp), allocatable :: basis(:, :), projected(:)
                real(dp) :: history(10), x(2), weight, exact(48), field(36), order
                real(dp) :: current_error(2), projection_error(2), projection(36)
                integer :: i, it, status, t, q
                logical :: quadratic
                call prepare(kind, n, op, steady, p, u, load)
                allocate (source(6, op%fs%nq, op%fs%mm%topology%n_triangles), &
                    boundary(op%ndof), surface(op%ndof), basis(36, op%nlocal), &
                    projected(op%ndof))
                allocate (r(steady%ndof), unused(steady%ndof), zero(steady%ndof), &
                    step(steady%ndof))
                source = 0; zero = 0; boundary = u(:op%ndof)
                projected = u(:op%ndof)
                trace_symmetry = op%symmetry
                call vrmhd_ampere_boundary(op, magnetic_trace, surface, status)
                call require(status == 0, 'manufactured Ampere trace')
                load = load - surface
                do i = 1, op%ndof
                    if (op%fixed(i)) cycle
                    u(i) = u(i) + 0.02_dp*sin(real(i, dp))
                end do
                quadratic = .false.
                do it = 1, size(history)
                    call vrmhd_stationary_residual(op, steady, p, u, zero, boundary, &
                        source, load, r, unused, status)
                    call require(status == 0, 'Newton residual')
                    history(it) = sqrt(sum(r*r))
                    print '(a,3i4,es16.8)', 'Newton residual ', kind, n, it - 1, history(it)
                    if (it > 2) then
                        if (history(it) > 1.0e-12_dp .and. &
                            history(it - 1) < 0.2_dp*history(it - 2)) then
                            order = log(history(it)/history(it - 1)) &
                                /log(history(it - 1)/history(it - 2))
                            print '(a,f9.4)', 'Newton observed order ', order
                            quadratic = quadratic .or. order > 1.8_dp
                        end if
                    end if
                    if (history(it) < 1.0e-10_dp) exit
                    call vrmhd_stationary_jacobian(op, steady, p, u, boundary, source, &
                        jac, status)
                    call require(status == 0, 'Newton assembled sparse Jacobian')
                    call sparse_solve_once(jac, -r, step, st)
                    if (.not. status_ok(st)) print '(a,i0)', 'direct solve status ', st%code
                    call require(status_ok(st), 'Newton direct solve')
                    u = u + step
                end do
                call require(it <= size(history), 'manufactured Newton convergence')
                call require(quadratic, 'quadratic Newton convergence')
                print '(a,es14.6)', 'Coulomb multiplier maximum ', maxval(abs(u(op%ndof + 1:)))
                error = 0; current_error = 0; projection_error = 0
                do t = 1, op%fs%mm%topology%n_triangles
                    do q = 1, op%fs%nq
                        call vrmhd_section_basis(op, t, q, x, weight, basis, status)
                        call exact_at(op, x, exact)
                        field = matmul(basis, u(op%map(:, t)))
                        projection = matmul(basis, projected(op%map(:, t)))
                        error = error + weight*norms(field - exact(:36))
                        current_error = current_error + weight*[ &
                            sum((field(8:10) - exact(8:10))**2), &
                            sum((field(34:36) - exact(34:36))**2)]
                        projection_error = projection_error + weight*[ &
                            sum((projection(8:10) - exact(8:10))**2), &
                            sum((projection(34:36) - exact(34:36))**2)]
                    end do
                end do
                error = sqrt(error)
                print '(a,2es14.6)', 'solved J L2, curl L2 errors ', sqrt(current_error)
                print '(a,2es14.6)', 'projected J L2, curl L2 errors ', sqrt(projection_error)
            end subroutine
        end program test_vrmhd_stationary
