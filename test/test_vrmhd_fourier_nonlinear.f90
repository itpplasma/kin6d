program test_vrmhd_fourier_nonlinear
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_gs_symmetry, only: gs_symmetry_t, GS_ROTATION, GS_TRANSLATION, GS_SCREW, &
        GS_CHART_RZ, GS_CHART_CARTESIAN
    use kin6d_vrmhd_parameters, only: vrmhd_parameters_t
    use kin6d_vrmhd_fourier
    use kin6d_vrmhd_fourier_nonlinear
    use fortfem_mesh_2d, only: mesh_2d_t
    use fortfem_meridian_curved_mesh, only: meridian_mesh_t, meridian_mesh_straight
    use fortfem_meridian_feec_n, only: meridian_qgrad, meridian_lagrange_nodes
    use fortsparse, only: csc_t, csc_matvec
    implicit none
    character(24) :: gate
    integer :: kind, family
    call get_command_argument(1, gate)
    select case (gate)
    case ('coupling')
        call coupling()
    case ('energy')
        do family = 0, 1
            do kind = GS_TRANSLATION, GS_SCREW
                call energy(kind, family == 1, .false.)
                call energy(kind, family == 1, .true.)
            end do
        end do
    case ('derivatives')
        do family = 0, 1
            do kind = GS_ROTATION, GS_SCREW
                call derivatives(kind, family == 1)
            end do
        end do
    case default
        error stop 'unknown nonlinear Fourier gate'
    end select
contains
    subroutine require(ok, message)
        logical, intent(in) :: ok
        character(*), intent(in) :: message
        if (.not. ok) then
            print '(a)', message
            error stop 'nonlinear Fourier gate failed'
        end if
    end subroutine

    subroutine build(kind, harmonics, op, cells, mixed)
        integer, intent(in) :: kind, harmonics(:)
        integer, intent(in), optional :: cells
        logical, intent(in), optional :: mixed
        type(vrmhd_fourier_t), intent(out) :: op
        type(mesh_2d_t) :: mesh
        type(meridian_mesh_t) :: mm
        type(gs_symmetry_t) :: symmetry
        real(dp) :: shift, scale
        integer :: status, count
        count = 1
        if (present(cells)) count = cells
        shift = 0; scale = 0.7_dp
        symmetry = gs_symmetry_t(kind, GS_CHART_CARTESIAN, 0.0_dp)
        if (kind == GS_ROTATION) then
            shift = 2; scale = 1; symmetry%chart = GS_CHART_RZ
        end if
        if (kind == GS_SCREW) symmetry%k = 1/scale
        call mesh%create_rectangular(count+1, count+1, shift, shift+1, 0.0_dp, 1.0_dp)
        call meridian_mesh_straight(mesh, 1, mm, status)
        call require(status == 0, 'nonlinear Fourier mesh')
        if (present(mixed)) then
            if (mixed) then
                call vrmhd_fourier_build(mm, symmetry, 2, harmonics, scale, &
                    op, status, scalar_degree=1, velocity_degree=2, &
                    coordinate_coframe=kind == GS_SCREW)
            else
                call vrmhd_fourier_build(mm, symmetry, 1, harmonics, scale, op, &
                    status, &
                    coordinate_coframe=kind == GS_SCREW)
            end if
        else
            call vrmhd_fourier_build(mm, symmetry, 1, harmonics, scale, op, status, &
                coordinate_coframe=kind == GS_SCREW)
        end if
        call require(status == 0, 'nonlinear Fourier layout')
    end subroutine

    subroutine parameters(p)
        type(vrmhd_parameters_t), intent(out) :: p
        p%mi = 1.2_dp; p%mu0 = 0.8_dp
        p%eta0 = 0.02_dp; p%nu0 = 0.03_dp; p%epsilon_b = 0.1_dp
        p%chi_parallel0 = [0.07_dp, 0.04_dp]
        p%chi_perpendicular = [0.005_dp, 0.003_dp]
        p%d_parallel = 0.02_dp; p%d_perpendicular = 0.002_dp
        p%nu_parallel = 0.01_dp; p%equilibration_rate = 0.04_dp
        p%electron_mass = 0.04_dp; p%gradient_floor = 0.03_dp
        p%flux_limited = .true.
    end subroutine

    subroutine derivatives(kind, mixed)
        integer, intent(in) :: kind
        logical, intent(in) :: mixed
        type(vrmhd_fourier_t) :: op
        type(vrmhd_parameters_t) :: p
        type(csc_t) :: jac
        real(dp), allocatable :: u(:), ud(:), du(:), dud(:), r(:), jv(:), unused(:)
        real(dp), allocatable :: plus(:), minus(:), source(:, :, :, :), s0(:, :, :)
        real(dp), allocatable :: u0(:), ud0(:), expected(:)
        real(dp) :: shift, step, defect, errors(2), norm, dual
        integer :: ns, n, i, j, a, status, angles
        call build(kind, [2, 0, 1], op, mixed=mixed)
        call parameters(p)
        ns = op%section%ndof; n = op%ndof; shift = 0.37_dp
        allocate (u(n), ud(n), du(n), dud(n), r(n), jv(n), unused(n), plus(n), &
            minus(n), expected(n), u0(ns), ud0(ns), &
            source(6, 65, op%section%fs%nq, &
            op%section%fs%mm%topology%n_triangles), &
            s0(6, op%section%fs%nq, op%section%fs%mm%topology%n_triangles))
        do i = 1, ns
            u0(i) = 0.01_dp*sin(real(i, dp)); ud0(i) = 0.005_dp*cos(real(i, dp))
        end do
        u0(:op%section%scalar_fs%n0) = u0(:op%section%scalar_fs%n0) + 2
        u0(op%section%offset(5)+1:op%section%offset(7)) = &
            u0(op%section%offset(5)+1:op%section%offset(7)) + 1.5_dp
        do i = 1, n
            du(i) = 0.01_dp*cos(real(3*i, dp)); dud(i) = shift*du(i)
        end do
        s0 = 0.003_dp
        do a = 1, size(source, 2)
            source(:, a, :, :) = s0
        end do
        u = 0; ud = 0
        u(op%offset(2)+1:op%offset(2)+ns) = u0
        ud(op%offset(2)+1:op%offset(2)+ns) = ud0
        call vrmhd_fourier_nonlinear_residual(op, p, u, ud, du, dud, source, &
            r, jv, status)
        call require(status == 0, 'nonlinear invariant tangent')
        call vrmhd_fourier_residual(op, p, u0, ud0, du, dud, s0, plus, expected, status)
        call require(status == 0, 'independent admitted linear operator')
        defect = norm2(r-plus)/max(norm2(plus), tiny(1.0_dp))
        call require(defect < 1.0e-12_dp, 'invariant nonlinear residual equals section')
        defect = norm2(jv-expected)/norm2(expected)
        print '(a,i2,es16.8)', 'symmetry / invariant tangent defect: ', kind, defect
        call require(defect < 1.0e-12_dp, &
            'nonlinear tangent equals linear harmonic operator')
        ! Excite every retained harmonic before differentiating the nonlinear terms.
        u = u + du; ud = ud + dud
        call vrmhd_fourier_nonlinear_residual(op, p, u, ud, du, dud, source, &
            r, jv, status)
        call require(status == 0, 'nonlinear coupled residual')
        do j = 1, 2
            step = 0.02_dp/2**(j-1)
            call vrmhd_fourier_nonlinear_residual(op, p, u+step*du, ud+step*dud, &
                0*du, 0*dud, source, plus, unused, status)
            call require(status == 0, 'positive difference state')
            call vrmhd_fourier_nonlinear_residual(op, p, u-step*du, ud-step*dud, &
                0*du, 0*dud, source, minus, unused, status)
            call require(status == 0, 'negative difference state')
            errors(j) = norm2((plus-minus)/(2*step)-jv)/norm2(jv)
        end do
        print '(a,2es16.8)', 'centered directional derivative errors: ', errors
        call require(errors(2) < 1.0e-6_dp .and. errors(1)/errors(2) > 3.5_dp, &
            'nonlinear directional derivative has second-order finite differences')
        call vrmhd_fourier_nonlinear_jacobian(op, p, u, ud, source, shift, jac, status)
        call require(status == 0, 'nonlinear assembled Jacobian')
        expected = csc_matvec(jac, du)
        defect = norm2(expected-jv)/norm2(jv)
        print '(a,es16.8)', 'assembled / matrix-free nonlinear tangent defect: ', defect
        call require(defect < 1.0e-12_dp, 'assembled nonlinear tangent equals JVP')
        ! Transpose the actual sparse matrix, with an independent dual vector.
        plus = [(sin(real(i, dp)), i=1, n)]; minus = 0
        do j = 1, n
            do i = jac%col_ptr(j), jac%col_ptr(j+1)-1
                minus(j) = minus(j) + jac%val(i)*plus(jac%row_idx(i))
            end do
        end do
        dual = abs(dot_product(plus, jv)-dot_product(minus, du))
        norm = norm2(plus)*norm2(jv)
        call require(dual < 1.0e-12_dp*norm, 'nonlinear sparse transpose duality')
        ! Nonpolynomial closures must converge as the angular rule is refined.
        do angles = 1, 4
            a = 2*size(source, 2)-1
            deallocate (source)
            allocate (source(6, a, op%section%fs%nq, &
                op%section%fs%mm%topology%n_triangles))
            do a = 1, size(source, 2)
                source(:, a, :, :) = s0
            end do
            call vrmhd_fourier_nonlinear_residual(op, p, u, ud, du, dud, source, &
                plus, expected, status)
            call require(status == 0, 'refined angular rule')
            defect = max(norm2(plus-r)/norm2(plus), norm2(expected-jv)/norm2(expected))
            print '(a,i6,es16.8)', 'angular samples / residual and tangent change: ', &
                size(source, 2), defect
            if (defect < 1.0e-10_dp) exit
            r = plus; jv = expected
        end do
        call require(defect < 1.0e-10_dp, 'nonpolynomial angular refinement')
    end subroutine

    subroutine energy(kind, mixed, dissipative)
        integer, intent(in) :: kind
        logical, intent(in) :: mixed, dissipative
        type(vrmhd_fourier_t) :: op
        type(vrmhd_parameters_t) :: p
        type(csc_t) :: jac, shifted
        real(dp), allocatable :: u(:), ud(:), zero(:), r(:), unused(:), s0(:, :, :)
        real(dp), allocatable :: source(:, :, :, :), dense(:, :), mass(:, :), rhs(:, :)
        real(dp), allocatable :: basis(:, :), u0(:)
        integer, allocatable :: algebraic(:), dynamic(:), pivot(:), map(:)
        logical, allocatable :: active(:), current(:)
        real(dp) :: x(2), weight, theta, v(36), vd(36), kinetic, internal, magnetic
        real(dp) :: scale, defect, original, shift
        integer :: n, ns, n0, a0, j0, b, c, lo, i, j, q, t, a, status, info, nf
        interface
            subroutine dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
                import :: dp
                integer, intent(in) :: n, nrhs, lda, ldb
                integer, intent(out) :: ipiv(*)
                integer, intent(out) :: info
                real(dp), intent(inout) :: a(lda, *), b(ldb, *)
            end subroutine
        end interface
        call build(kind, [0, 1], op, cells=2, mixed=mixed)
        p%mi = 1.2_dp; p%mu0 = 0.8_dp
        if (dissipative) then
            p%eta0 = 0.02_dp; p%nu0 = 0.03_dp
            p%chi_parallel0 = [0.07_dp, 0.04_dp]
            p%chi_perpendicular = [0.005_dp, 0.003_dp]
            p%epsilon_b = 0.1_dp
        end if
        n = op%ndof; ns = op%section%ndof; n0 = op%section%scalar_fs%n0
        a0 = op%section%offset(7); j0 = op%section%offset(8)
        allocate (u(n), ud(n), zero(n), r(n), unused(n), u0(ns), &
            active(n), current(n), dense(n, n), &
            source(6, 9, op%section%fs%nq, op%section%fs%mm%topology%n_triangles), &
            s0(6, op%section%fs%nq, op%section%fs%mm%topology%n_triangles), &
            basis(36, op%section%nlocal*sum(op%copies)), &
            map(op%section%nlocal*sum(op%copies)))
        zero = 0; ud = 0; source = 0; s0 = 0; active = .false.; current = .false.
        u = [(0.001_dp*sin(real(17*i, dp)), i=1, n)]
        do b = 1, size(op%harmonic)
            do c = 1, op%copies(b)
                lo = op%offset(b)+(c-1)*ns
                u(lo+j0+1:lo+ns) = 0
                where (op%fixed(lo+n0+1:lo+j0)) u(lo+n0+1:lo+j0) = 0
                active(lo+1:lo+j0) = .true.
                active(lo+n0+1:lo+j0) = .not. op%fixed(lo+n0+1:lo+j0)
                active(lo+op%section%offset(5)+1:lo+a0) = .true.
                current(lo+j0+1:lo+ns) = .not. op%fixed(lo+j0+1:lo+ns)
            end do
        end do
        u(:n0) = u(:n0) + 2
        u(op%section%offset(5)+1:a0) = u(op%section%offset(5)+1:a0) + 1.5_dp
        u0 = u(:ns)
        call vrmhd_fourier_jacobian(op, p, u0, 0*u0, s0, 0.0_dp, jac, status)
        call require(status == 0, 'independent linear Ampere mass')
        dense = 0
        do j = 1, n
            do i = jac%col_ptr(j), jac%col_ptr(j+1)-1
                dense(jac%row_idx(i), j) = jac%val(i)
            end do
        end do
        algebraic = pack([(i, i=1, n)], current)
        dynamic = pack([(i, i=1, n)], active)
        call vrmhd_fourier_nonlinear_residual(op, p, u, zero, zero, zero, source, &
            r, unused, status)
        call require(status == 0, 'Ampere recovery residual')
        nf = size(algebraic)
        allocate (mass(nf, nf), rhs(nf, 1), pivot(nf))
        mass = dense(algebraic, algebraic); rhs(:, 1) = -r(algebraic)
        call dgesv(nf, 1, mass, nf, pivot, rhs, nf, info)
        call require(info == 0, 'interior current mass solve')
        u(algebraic) = rhs(:, 1)
        deallocate (mass, rhs, pivot)
        shift = 2.0_dp**30
        call vrmhd_fourier_nonlinear_jacobian(op, p, u, zero, source, 0.0_dp, &
            jac, status)
        call require(status == 0, 'nonlinear energy spatial Jacobian')
        call vrmhd_fourier_nonlinear_jacobian(op, p, u, zero, source, shift, &
            shifted, status)
        call require(status == 0, 'nonlinear energy shifted Jacobian')
        dense = 0
        do j = 1, n
            do i = shifted%col_ptr(j), shifted%col_ptr(j+1)-1
                dense(shifted%row_idx(i), j) = shifted%val(i)/shift
            end do
            do i = jac%col_ptr(j), jac%col_ptr(j+1)-1
                dense(jac%row_idx(i), j) = dense(jac%row_idx(i), j)-jac%val(i)/shift
            end do
        end do
        nf = size(dynamic)
        allocate (mass(nf, nf), rhs(nf, 1), pivot(nf))
        do a = 1, 3
            call vrmhd_fourier_nonlinear_residual(op, p, u, ud, zero, zero, &
                source, r, unused, status)
            call require(status == 0, 'nonlinear energy mass residual')
            if (a == 1) original = norm2(r(dynamic))
            mass = dense(dynamic, dynamic); rhs(:, 1) = -r(dynamic)
            call dgesv(nf, 1, mass, nf, pivot, rhs, nf, info)
            call require(info == 0, 'nonlinear energy mass correction')
            ud(dynamic) = ud(dynamic)+rhs(:, 1)
        end do
        call vrmhd_fourier_nonlinear_residual(op, p, u, ud, zero, zero, &
            source, r, unused, status)
        call require(status == 0, 'corrected nonlinear energy residual')
        call require(norm2(r(dynamic))/original < 1.0e-12_dp, &
            'dynamic mass solve accuracy')
        call require(norm2(r(algebraic))/original < 1.0e-12_dp, 'weak Ampere accuracy')
        kinetic = 0; internal = 0; magnetic = 0; scale = 0
        do t = 1, op%section%fs%mm%topology%n_triangles
            do q = 1, op%section%fs%nq
                do a = 1, size(source, 2)
                    theta = 2*acos(-1.0_dp)*real(a-1, dp)/size(source, 2)
                    call vrmhd_fourier_nonlinear_basis(op, t, q, theta, x, &
                        weight, basis, map, status)
                    call require(status == 0, 'physical energy readback')
                    v = matmul(basis, u(map)); vd = matmul(basis, ud(map))
                    weight = weight/size(source, 2)
                    kinetic = kinetic + weight*(0.5_dp*vd(1)*sum(v(2:4)**2) &
                        + v(1)*dot_product(v(2:4), vd(2:4)))
                    internal = internal + weight*1.5_dp/p%mi &
                        *(vd(1)*(v(11)+v(12))+v(1)*(vd(11)+vd(12)))
                    magnetic = magnetic + weight*dot_product(v(25:27), vd(25:27))/p%mu0
                    scale = scale + weight*(abs(0.5_dp*vd(1)*sum(v(2:4)**2)) &
                        + abs(v(1)*dot_product(v(2:4), vd(2:4))) &
                        + abs(1.5_dp/p%mi*(vd(1)*(v(11)+v(12)) &
                        + v(1)*(vd(11)+vd(12)))) &
                        + abs(dot_product(v(25:27), vd(25:27))/p%mu0))
                end do
            end do
        end do
        call require(scale > tiny(scale), 'nontrivial nonlinear energy exchange')
        defect = abs(kinetic+internal+magnetic)/scale
        print '(a,3es16.8)', 'kinetic / internal / magnetic power [W/m]: ', &
            kinetic, internal, magnetic
        print '(a,es16.8)', 'nonlinear Fourier total energy defect: ', defect
        call require(defect < 1.0e-12_dp, 'closed nonlinear Fourier energy ledger')
    end subroutine

    subroutine coupling()
        type(vrmhd_fourier_t) :: op
        type(vrmhd_parameters_t) :: p
        real(dp), allocatable :: u(:), zero(:), r(:), unused(:), expected(:)
        real(dp), allocatable :: source(:, :, :, :), nodes(:, :), values(:), grad(:, :)
        real(dp) :: x(2), weight, rho0, density, velocity, defect
        integer :: n, ns, n0, status, i, t, q, row, b, lo
        call build(GS_TRANSLATION, [2, 0, 1], op)
        p%mi = 1; p%mu0 = 1
        n = op%ndof; ns = op%section%ndof; n0 = op%section%scalar_fs%n0
        allocate (u(n), zero(n), r(n), unused(n), expected(n), &
            source(6, 9, op%section%fs%nq, op%section%fs%mm%topology%n_triangles), &
            values(op%section%scalar_fs%nl0), grad(2, op%section%scalar_fs%nl0))
        call meridian_lagrange_nodes(op%section%scalar_fs, nodes)
        rho0 = 2; density = 0.2_dp; velocity = 0.3_dp
        u = 0; zero = 0; source = 0; expected = 0
        lo = op%offset(2)
        u(lo+1:lo+n0) = rho0
        u(lo+op%section%offset(5)+1:lo+op%section%offset(7)) = 1
        lo = op%offset(3)
        u(lo+1:lo+n0) = density*nodes(1, :)
        u(lo+op%section%offset(2)+1:lo+op%section%offset(3)) = velocity
        call vrmhd_fourier_nonlinear_residual(op, p, u, zero, zero, zero, source, &
            r, unused, status)
        call require(status == 0, 'quadratic harmonic coupling residual')
        ! rho = rho0 + b x cos(theta), vx = a cos(theta).
        ! Angular averages give rho0*a/2 at n=1, a*b*x/2 at n=0,
        ! and a*b*x/4 at n=2; sine density rows all vanish.
        do t = 1, op%section%fs%mm%topology%n_triangles
            do q = 1, op%section%fs%nq
                call meridian_qgrad(op%section%scalar_fs, t, q, x, weight, &
                    values, grad, status)
                call require(status == 0, 'analytic density weak integration')
                do i = 1, size(values)
                    row = op%section%scalar_fs%sp(0)%dof(i, t)
                    expected(op%offset(1)+row) = expected(op%offset(1)+row) &
                        - weight*velocity*density*x(1)*grad(1, i)/4
                    expected(op%offset(2)+row) = expected(op%offset(2)+row) &
                        - weight*velocity*density*x(1)*grad(1, i)/2
                    expected(op%offset(3)+row) = expected(op%offset(3)+row) &
                        - weight*velocity*rho0*grad(1, i)/2
                end do
            end do
        end do
        defect = 0
        do b = 1, size(op%harmonic)
            do i = 1, op%copies(b)
                lo = op%offset(b)+(i-1)*ns
                defect = max(defect, norm2(r(lo+1:lo+n0)-expected(lo+1:lo+n0)))
            end do
        end do
        defect = defect/norm2(expected)
        print '(a,es16.8)', 'analytic density harmonic product defect: ', defect
        call require(defect < 1.0e-13_dp, &
            'sum and difference harmonics equal analytic products')
        source = 0
        call vrmhd_fourier_nonlinear_residual(op, p, u, zero, zero, zero, &
            source(:, :8, :, :), r, unused, status)
        call require(status /= 0, 'underresolved angular product rule is rejected')
    end subroutine
end program
