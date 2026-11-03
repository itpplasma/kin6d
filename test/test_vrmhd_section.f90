program test_vrmhd_section
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_gs_symmetry, only: gs_symmetry_t, GS_ROTATION, GS_TRANSLATION, &
        GS_SCREW, GS_CHART_RZ, GS_CHART_CARTESIAN
    use kin6d_vrmhd_parameters, only: vrmhd_parameters_t
    use kin6d_vrmhd_section
    use kin6d_vrmhd_evolve, only: vrmhd_evolve_t, vrmhd_evolve_build
    use fortfem_meridian_feec_n, only: meridian_lagrange_nodes
    use kin6d_vrmhd_section_leaves, only: vrmhd_grad_rotation, vrmhd_grad_screw, &
        vrmhd_velocity_rotation, vrmhd_velocity_screw, &
        vrmhd_curl_rotation, vrmhd_curl_screw
    use kin6d_vrmhd_leaves, only: vrmhd_chi_flux_limited, &
        vrmhd_chi_flux_limited_jvp, &
        vrmhd_chi_flux_limited_vjp
    use fortfem_mesh_2d, only: mesh_2d_t
    use fortfem_meridian_curved_mesh, only: meridian_mesh_t, meridian_mesh_straight
    use fortsparse, only: csc_t, csc_matvec, csc_from_triplet, &
        fortsparse_status_t, status_ok
    use omp_lib, only: omp_set_num_threads
    implicit none
    character(24) :: gate
    integer :: kind, degree
    call get_command_argument(1, gate)
    if (gate == 'mixed') then
        do degree = 1, 2
            call interpolation(degree, .false.)
            call interpolation(degree, .true.)
        end do
    end if
    if (gate == 'closure') then
        call closure()
    else
        do kind = GS_ROTATION, GS_SCREW
            select case (gate)
            case ('mixed')
                do degree = 1, 2
                    call assembled(kind, .false., mixed_degree=degree)
                    call assembled(kind, .true., oracle=.true., mixed_degree=degree)
                end do
            case ('geometry')
                do degree = 1, 2
                    call geometry(kind, degree)
                end do
            case ('derivatives')
                call assembled(kind, .false.)
                call assembled(kind, .true.)
            case ('volume')
                call assembled(kind, .true., oracle=.true.)
            case default
                error stop 'unknown RM2 gate'
            end select
        end do
    end if
contains
    subroutine require(ok, message)
        logical, intent(in) :: ok
        character(*), intent(in) :: message
        if (.not. ok) then
            print '(a)', message
            error stop 'RM2 oracle failed'
        end if
    end subroutine

    subroutine interpolation(degree, mixed)
        integer, intent(in) :: degree
        logical, intent(in) :: mixed
        type(vrmhd_section_t), target :: op
        type(vrmhd_parameters_t) :: p
        type(vrmhd_evolve_t) :: problem
        type(mesh_2d_t) :: mesh
        type(meridian_mesh_t) :: mm
        type(gs_symmetry_t) :: symmetry
        type(csc_t) :: identity
        type(fortsparse_status_t) :: sparse_status
        real(dp), allocatable :: xy(:, :), u(:), basis(:, :), source(:, :, :)
        real(dp), allocatable :: zero(:), probe(:), mass(:)
        integer, allocatable :: indices(:)
        real(dp) :: x(2), weight, value(36), expected(36), defect
        integer :: status, t, q, j, i, n, field
        call mesh%create_rectangular(3, 3, 0.1_dp, 1.1_dp, 0.2_dp, 1.2_dp)
        call meridian_mesh_straight(mesh, 1, mm, status)
        call require(status == 0, 'mixed interpolation mesh')
        symmetry = gs_symmetry_t(GS_TRANSLATION, GS_CHART_CARTESIAN, 0.0_dp)
        if (mixed) then
            call vrmhd_section_build(mm, symmetry, degree + 1, 1, op, status, &
                scalar_degree=degree, velocity_degree=degree + 1)
        else
            call vrmhd_section_build(mm, symmetry, degree + 1, 1, op, status, &
                scalar_degree=degree, study_pairing=.true.)
        end if
        call require(status == 0, 'mixed interpolation section')
        n = op%ndof
        p%mi = 1; p%mu0 = 1
        allocate (u(n), basis(36, op%nlocal), source(6, op%fs%nq, &
            mm%topology%n_triangles), zero(n), probe(n), mass(n), indices(n))
        u = 0; zero = 0; source = 0
        call meridian_lagrange_nodes(op%scalar_fs, xy)
        u(:op%offset(2)) = 2 + xy(1, :)**degree
        u(op%offset(5) + 1:op%offset(6)) = 3 + xy(2, :)**degree
        u(op%offset(6) + 1:op%offset(7)) = 4 + xy(1, :)**degree
        deallocate (xy)
        if (allocated(op%velocity_fs)) then
            call meridian_lagrange_nodes(op%velocity_fs, xy)
        else
            call meridian_lagrange_nodes(op%scalar_fs, xy)
        end if
        u(op%offset(2) + 1:op%offset(3)) = xy(1, :)**(degree + 1)
        u(op%offset(3) + 1:op%offset(4)) = xy(2, :)**(degree + 1)
        u(op%offset(4) + 1:op%offset(5)) = xy(1, :)**degree*xy(2, :)
        defect = 0
        do t = 1, mm%topology%n_triangles
            do q = 1, op%fs%nq
                call vrmhd_section_basis(op, t, q, x, weight, basis, status)
                call require(status == 0, 'mixed interpolation basis')
                value = matmul(basis, u(op%map(:, t)))
                expected = 0
                expected(1) = 2 + x(1)**degree
                expected(2:4) = [x(1)**(degree + 1), x(2)**(degree + 1), &
                    x(1)**degree*x(2)]
                expected(11:12) = [3 + x(2)**degree, 4 + x(1)**degree]
                expected(13) = degree*x(1)**(degree - 1)
                expected(16) = (degree + 1)*x(1)**degree
                expected(18) = degree*x(1)**(degree - 1)*x(2)
                expected(20) = (degree + 1)*x(2)**degree
                expected(21) = x(1)**degree
                expected(29) = degree*x(2)**(degree - 1)
                expected(31) = degree*x(1)**(degree - 1)
                defect = max(defect, maxval(abs(value - expected)))
            end do
        end do
        print '(a,i2,l2,es16.8)', 'polynomial interpolation degree/mixed/error ', &
            degree, mixed, defect
        if (.not. mixed) then
            call require(defect > 1.0e-3_dp, &
                'old velocity space rejects degree k+1 data')
            return
        end if
        call require(defect < 5.0e-12_dp, 'mixed polynomial values and gradients')
        indices = [(j, j=1, n)]
        call csc_from_triplet(n, n, indices, indices, spread(1.0_dp, 1, n), &
            identity, sparse_status)
        call require(status_ok(sparse_status), 'mixed identity reduction')
        call vrmhd_evolve_build(problem, op, p, source, zero, identity, status)
        call require(status == 0, 'mixed BDF consumer')
        ! The last A coordinate has mass; the first current coordinate has none.
        do field = 7, 8
            i = op%offset(8) + field - 7
            probe = 0; probe(i) = 1
            call vrmhd_section_mass(op, p, u, probe, mass, status)
            call require(status == 0, 'mixed independent mass action')
            if (problem%algebraic(i)) then
                call require(norm2(mass) < 1.0e-13_dp, &
                    'BDF algebraic coordinate has zero mass')
            else
                call require(norm2(mass) > 1.0e-8_dp, &
                    'BDF dynamic coordinate has nonzero mass')
            end if
        end do
    end subroutine

    subroutine setup(kind, op, degree, mixed_degree)
        integer, intent(in) :: kind
        type(vrmhd_section_t), intent(out) :: op
        type(mesh_2d_t) :: mesh
        type(meridian_mesh_t) :: mm
        type(gs_symmetry_t) :: symmetry
        integer, intent(in), optional :: degree, mixed_degree
        integer :: status, order
        real(dp) :: shift
        shift = 0
        symmetry = gs_symmetry_t(kind, GS_CHART_CARTESIAN, 0.0_dp)
        if (kind == GS_ROTATION) then
            shift = 2
            symmetry%chart = GS_CHART_RZ
        end if
        if (kind == GS_SCREW) symmetry%k = 0.7_dp
        call mesh%create_rectangular(3, 3, -0.4_dp + shift, 0.6_dp + shift, &
            -0.3_dp, 0.7_dp)
        call meridian_mesh_straight(mesh, 1, mm, status)
        call require(status == 0, 'straight section mesh')
        order = 1
        if (present(degree)) order = degree
        if (present(mixed_degree)) then
            call vrmhd_section_build(mm, symmetry, mixed_degree + 1, 1, op, status, &
                scalar_degree=mixed_degree, velocity_degree=mixed_degree + 1)
        else
            call vrmhd_section_build(mm, symmetry, order, 1, op, status)
        end if
        call require(status == 0, 'section build')
    end subroutine

    subroutine parameters(p, limited)
        type(vrmhd_parameters_t), intent(out) :: p
        logical, intent(in) :: limited
        p%mi = 1.2_dp; p%mu0 = 0.8_dp; p%eta0 = 0.2_dp; p%nu0 = 0.3_dp
        p%chi_parallel0 = [0.7_dp, 0.4_dp]; p%chi_perpendicular = [0.05_dp, 0.03_dp]
        p%d_parallel = 0.2_dp; p%d_perpendicular = 0.02_dp
        p%nu_parallel = 0.1_dp; p%equilibration_rate = 0.4_dp; p%epsilon_b = 0.1_dp
        p%flux_limited = limited; p%electron_mass = 0.04_dp
        p%gradient_floor = 0.03_dp
    end subroutine

    subroutine assembled(kind, limited, oracle, mixed_degree)
        integer, intent(in) :: kind
        logical, intent(in) :: limited
        logical, intent(in), optional :: oracle
        integer, intent(in), optional :: mixed_degree
        type(vrmhd_section_t) :: op
        type(vrmhd_parameters_t) :: p
        type(csc_t) :: jac
        real(dp), allocatable :: u(:), ud(:), du(:), dud(:), r(:), jv(:), plus(:), &
            minus(:)
        real(dp), allocatable :: unused(:), zero(:), dual(:), vjp(:), source(:, :, &
            :), ref(:)
        real(dp) :: h, errors(2), scale, shift, defect
        integer :: i, j, status, n, n0
        call setup(kind, op, mixed_degree=mixed_degree); call parameters(p, limited)
        n = op%ndof; n0 = op%scalar_fs%n0
        allocate (u(n), ud(n), du(n), dud(n), r(n), jv(n), plus(n), minus(n), &
            unused(n), zero(n), dual(n), vjp(n), ref(n), source(6, op%fs%nq, &
            op%fs%mm%topology%n_triangles))
        do i = 1, n
            u(i) = 0.05_dp*sin(real(i, dp)); du(i) = 0.1_dp*cos(real(2*i, dp))
            ud(i) = 0.02_dp*cos(real(i, dp)); dual(i) = sin(real(3*i, dp))
        end do
        u(:n0) = u(:n0) + 2
        u(op%offset(5) + 1:op%offset(7)) = u(op%offset(5) + 1:op%offset(7)) + 1.5_dp
        shift = 0.31_dp; dud = shift*du; zero = 0
        source(:3, :, :) = 0.13_dp; source(4:6, :, :) = 0.07_dp
        call vrmhd_section_residual(op, p, u, ud, du, dud, source, r, jv, status)
        call require(status == 0, 'assembled residual')
        if (present(oracle)) then
            call reference(op, p, u, ud, source, ref)
            defect = maxval(abs(r - ref))/max(1.0_dp, maxval(abs(ref)))
            print '(a,i0,a,es12.4)', 'volume symmetry ', kind, ': ', defect
            call require(defect < 2.0e-13_dp, 'independent weak volume assembly')
            call omp_set_num_threads(4)
            call vrmhd_section_residual(op, p, u, ud, du, dud, source, plus, &
                unused, status)
            call omp_set_num_threads(1)
            call require(status == 0, 'parallel residual')
            call require(norm2(plus - r) < 1.0e-13_dp*max(1.0_dp, norm2(r)), &
                'four-thread residual agrees')
            call vrmhd_section_residual(op, p, u, zero, zero, zero, source, minus, &
                unused, status)
            call require(status == 0, 'zero-rate residual')
            call vrmhd_section_mass(op, p, u, ud, plus, status)
            call require(status == 0, 'separate mass action')
            call require(norm2(r - minus - plus) < 1.0e-13_dp*max(1.0_dp, norm2(r)), &
                'mass action equals independent conservative mass difference')

            return
        end if
        do j = 1, 2
            h = 1.0e-4_dp/2**(j - 1)
            call vrmhd_section_residual(op, p, u + h*du, ud + h*dud, zero, zero, &
                source, plus, unused, status)
            call require(status == 0, 'positive finite difference')
            call vrmhd_section_residual(op, p, u - h*du, ud - h*dud, zero, zero, &
                source, minus, unused, status)
            call require(status == 0, 'negative finite difference')
            errors(j) = norm2((plus - minus)/(2*h) - jv)
        end do
        print '(a,i0,l2,2es12.4)', 'assembled JVP ', kind, limited, errors
        call require(errors(2) < 0.3_dp*errors(1), &
            'quadratic assembled finite differences')
        call require(errors(2) < 2.0e-7_dp*max(1.0_dp, norm2(jv)), &
            'assembled JVP accuracy')
        call vrmhd_section_jacobian(op, p, u, ud, source, shift, jac, status)
        call require(status == 0, 'sparse Jacobian')
        plus = csc_matvec(jac, du)
        call require(norm2(plus - jv) < 2.0e-13_dp*max(1.0_dp, norm2(jv)), &
            'assembled sparse versus matrix-free JVP')
        call vrmhd_section_vjp(jac, dual, vjp, status)
        call require(status == 0, 'assembled transpose')
        scale = max(1.0_dp, norm2(dual)*norm2(jv))
        call require(abs(dot_product(dual, jv) - dot_product(vjp, &
            du)) < 2.0e-13_dp*scale, &
            'assembled VJP duality')
    end subroutine

    subroutine reference(op, p, u, ud, source, result)
        type(vrmhd_section_t), intent(in) :: op
        type(vrmhd_parameters_t), intent(in) :: p
        real(dp), intent(in) :: u(:), ud(:), source(:, :, :)
        real(dp), intent(out) :: result(:)
        real(dp) :: basis(36, op%nlocal), x(2), weight, v(36), vd(12), test(36)
        real(dp) :: b(3), gv(3, 3), gw(3, 3), strain(3, 3), stress(3, 3), tensor(3, 3)
        real(dp) :: eta, nu, pressure, n, divv, tr, chi, temperature, bg, mass, &
            heat(2), value, density_flux(3), product_gradient(3)
        integer :: t, q, i, j, s, st
        result = 0
        do t = 1, op%fs%mm%topology%n_triangles
            do q = 1, op%fs%nq
                call vrmhd_section_basis(op, t, q, x, weight, basis, st)
                call require(st == 0, 'reference geometry')
                v = matmul(basis, u(op%map(:, t)))
                vd = matmul(basis(:12, :), ud(op%map(:, t)))
                n = v(1)/p%mi; pressure = n*(v(11) + v(12))
                eta = p%eta0*(v(11)/p%t0)**(-1.5_dp)
                nu = p%nu0*(v(12)/p%t0)**(-1.5_dp)
                b = v(25:27)/sqrt(sum(v(25:27)**2) + p%epsilon_b**2)
                density_flux = p%d_perpendicular*v(13:15) &
                    + (p%d_parallel - p%d_perpendicular)*b*dot_product(b, v(13:15))
                gv = reshape(v(16:24), [3, 3]); divv = gv(1, 1) + gv(2, 2) + gv(3, 3)
                strain = (gv + transpose(gv))/2
                do j = 1, 3
                    strain(j, j) = strain(j, j) - divv/3
                end do
                tensor = spread(b, 2, 3)*spread(b, 1, 3)
                tr = sum(tensor*strain)
                stress = 2*nu*strain + 3*p%nu_parallel*tr*tensor
                do j = 1, 3
                    stress(j, j) = stress(j, j) - p%nu_parallel*tr*sum(b*b)
                end do
                heat(1) = eta*dot_product(v(8:10), v(8:10) - source(:3, q, t)) &
                    + n*p%equilibration_rate*(v(12) - v(11)) + source(5, q, t)
                heat(2) = sum(stress*gv) &
                    - n*p%equilibration_rate*(v(12) - v(11)) + source(6, q, t)
                do i = 1, op%nlocal
                    test = basis(:, i); gw = reshape(test(16:24), [3, 3])
                    select case (op%row(i))
                    case (1)
                        value = (vd(1) - source(4, q, t))*test(1) &
                            - v(1)*dot_product(v(2:4), test(13:15)) &
                            + p%d_perpendicular*dot_product(v(13:15), &
                            test(13:15)) &
                            + (p%d_parallel - p%d_perpendicular) &
                            *dot_product(b, v(13:15))*dot_product(b, test(13:15))
                    case (2)
                        product_gradient = matmul(transpose(gv), test(2:4)) &
                            + matmul(transpose(gw), v(2:4))
                        value = dot_product(v(1)*vd(2:4) + vd(1)*v(2:4)/2, test(2:4)) &
                            + v(1)*(dot_product(matmul(gv, v(2:4)), test(2:4)) &
                            - dot_product(matmul(gw, v(2:4)), v(2:4)))/2 &
                            + source(4, q, t)*dot_product(v(2:4), test(2:4))/2 &
                            - dot_product(density_flux, product_gradient)/2 &
                            - pressure*(gw(1, 1) + gw(2, 2) + gw(3, 3)) &
                            + sum(stress*gw) &
                            - dot_product(cross(v(8:10), v(25:27)), test(2:4))
                    case (3)
                        value = dot_product(vd(5:7) - cross(v(2:4), v(25:27)) &
                            + eta*(v(8:10) - source(:3, q, t)), &
                            test(5:7))
                    case (4)
                        value = p%mu0*dot_product(v(8:10), test(8:10)) &
                            - dot_product(v(25:27), test(34:36))
                    case (5, 6)
                        s = op%row(i) - 4; temperature = v(10 + s)
                        chi = p%chi_parallel0(s)*(temperature/p%t0)**2.5_dp
                        if (p%flux_limited(s)) then
                            bg = dot_product(b, v(25 + 3*s:27 + 3*s))
                            mass = p%electron_mass
                            if (s == 2) mass = p%mi
                            chi = chi/(1 + chi*sqrt(bg*bg + p%gradient_floor(s)**2) &
                                /(p%flux_limit(s)*temperature*sqrt(2*temperature/mass)))
                        end if
                        value = (1.5_dp*(v(1)*vd(10 + s) + vd(1)*temperature)/p%mi &
                            + n*temperature*divv - heat(s))*test(10 + s) &
                            - 1.5_dp*n*temperature*dot_product(v(2:4), &
                            test(25 + 3*s:27 + 3*s)) &
                            + n*p%chi_perpendicular(s) &
                            *dot_product(v(25 + 3*s:27 + 3*s), &
                            test(25 + 3*s:27 + 3*s)) &
                            + n*(chi - p%chi_perpendicular(s))*dot_product(b, &
                            v(25 + 3*s:27 + 3*s)) &
                            *dot_product(b, test(25 + 3*s:27 + 3*s))
                    end select
                    result(op%map(i, t)) = result(op%map(i, t)) + weight*value
                end do
            end do
        end do
    end subroutine

    function cross(a, b) result(c)
        real(dp), intent(in) :: a(3), b(3)
        real(dp) :: c(3)
        c = [a(2)*b(3) - a(3)*b(2), a(3)*b(1) - a(1)*b(3), a(1)*b(2) - a(2)*b(1)]
    end function

    subroutine geometry(kind, degree)
        integer, intent(in) :: kind, degree
        type(vrmhd_section_t) :: op
        real(dp), allocatable :: basis(:, :)
        real(dp) :: x(2), w, g(3), v(3), b(3), values(12), plus(12), minus(12)
        real(dp) :: d(3, 12), point(3), h, k, u1, u2, c, a, ax, ay, expected(3), &
            grad(3, 3), volume, vx(3), vy(3)
        integer :: i, j, t, q, st, n0, n1
        real(dp), allocatable :: coefficients(:), scalar(:)
        call setup(kind, op, degree); k = op%symmetry%k
        allocate (basis(36, op%nlocal))
        point = [0.2_dp, 0.1_dp, -0.3_dp]
        if (kind == GS_ROTATION) point = [2.2_dp, 0.0_dp, 0.1_dp]
        if (kind /= GS_ROTATION) point(3) = 0
        h = 2.0e-5_dp
        call analytic(kind, k, point, values)
        do j = 1, 3
            point(j) = point(j) + h; call analytic(kind, k, point, plus)
            point(j) = point(j) - 2*h; call analytic(kind, k, point, minus)
            point(j) = point(j) + h; d(j, :) = (plus - minus)/(2*h)
        end do
        x = [point(1), point(2)]
        if (kind == GS_ROTATION) x(2) = point(3)
        u1 = sin(x(1) + x(2)); u2 = cos(x(1) - x(2))
        c = -sin(x(1) - x(2)) - cos(x(1) + x(2))
        a = 1 + sin(x(1))*cos(x(2)); ax = cos(x(1))*cos(x(2))
        ay = -sin(x(1))*sin(x(2))
        if (kind == GS_ROTATION) then
            call vrmhd_curl_rotation(x(1), x(2), k, u1, u2, c, a, ax, ay, &
                v(1), v(2), v(3), b(1), b(2), b(3))
            call vrmhd_grad_rotation(x(1), x(2), k, ax, ay, g(1), g(2), g(3))
        else
            call vrmhd_curl_screw(x(1), x(2), k, u1, u2, c, a, ax, ay, &
                v(1), v(2), v(3), b(1), b(2), b(3))
            call vrmhd_grad_screw(x(1), x(2), k, ax, ay, g(1), g(2), g(3))
        end if
        expected = [d(2, 7) - d(3, 6), d(3, 5) - d(1, 7), d(1, 6) - d(2, 5)]
        call require(norm2(v - values(5:7)) < 1.0e-13_dp, 'invariant physical vector')
        call require(norm2(b - expected) < 2.0e-9_dp, &
            'curl from independent Cartesian flow')
        call require(norm2(g - d(:, 1)) < 2.0e-9_dp, &
            'scalar gradient from Cartesian flow')
        v = [sin(x(1)), cos(x(2)), sin(x(1))*cos(x(2))]
        vx = [cos(x(1)), 0.0_dp, cos(x(1))*cos(x(2))]
        vy = [0.0_dp, -sin(x(2)), -sin(x(1))*sin(x(2))]
        if (kind == GS_ROTATION) then
            call vrmhd_velocity_rotation(x(1), x(2), k, v(1), v(2), v(3), &
                vx(1), vx(2), vx(3), vy(1), vy(2), vy(3), &
                grad(1, 1), grad(2, 1), grad(3, 1), &
                grad(1, 2), grad(2, 2), grad(3, 2), &
                grad(1, 3), grad(2, 3), grad(3, 3))
        else
            call vrmhd_velocity_screw(x(1), x(2), k, v(1), v(2), v(3), &
                vx(1), vx(2), vx(3), vy(1), vy(2), vy(3), &
                grad(1, 1), grad(2, 1), grad(3, 1), &
                grad(1, 2), grad(2, 2), grad(3, 2), &
                grad(1, 3), grad(2, 3), grad(3, 3))
        end if
        call require(norm2(reshape(grad - transpose(d(:, 2:4)), [9])) < 2.0e-9_dp, &
            'vector gradient from Cartesian flow')
        n0 = op%scalar_fs%n0; n1 = op%fs%sp(1)%ndof
        allocate (coefficients(op%ndof), scalar(op%fs%n0)); coefficients = 0
        do i = 1, op%fs%n0
            scalar(i) = sin(real(i, dp))
        end do
        coefficients(6*n0 + 1:6*n0 + n1) = csc_matvec(op%fs%grad%a, scalar)
        volume = 0
        do t = 1, op%fs%mm%topology%n_triangles
            do q = 1, op%fs%nq
                call vrmhd_section_basis(op, t, q, x, w, basis, st)
                call require(st == 0, 'basis quadrature')
                b = matmul(basis(25:27, :), coefficients(op%map(:, t)))
                call require(norm2(b) < 1.0e-12_dp, 'assembled curl grad is zero')
                volume = volume + w
            end do
        end do
        expected(1) = 1
        if (kind == GS_ROTATION) expected(1) = 2.1_dp
        call require(abs(volume - expected(1)) < 1.0e-13_dp, 'orbit-volume measure')
        print '(a,i0,a,es12.4)', 'section symmetry ', kind, ' volume: ', volume
    end subroutine

    subroutine analytic(kind, k, point, u)
        integer, intent(in) :: kind
        real(dp), intent(in) :: k, point(3)
        real(dp), intent(out) :: u(12)
        real(dp) :: x, y, theta, c, s, a, xi(3), v(3)
        if (kind == GS_ROTATION) then
            theta = atan2(point(2), point(1)); x = norm2(point(:2)); y = point(3)
            xi = [0.0_dp, x, 0.0_dp]
            a = 1 + sin(x)*cos(y); v = [sin(x + y), a/x, cos(x - y)]
        else
            theta = k*point(3)
            c = cos(theta); s = sin(theta)
            x = c*point(1) + s*point(2); y = -s*point(1) + c*point(2)
            xi = [-k*y, k*x, 1.0_dp]; a = 1 + sin(x)*cos(y)
            v = [sin(x + y), cos(x - y), &
                k*(y*sin(x + y) - x*cos(x - y))] + a*xi/sum(xi*xi)
        end if
        c = cos(theta); s = sin(theta)
        u = 0; u(1) = a
        u(5:7) = [c*v(1) - s*v(2), s*v(1) + c*v(2), v(3)]
        v = [sin(x), cos(y), sin(x)*cos(y)]
        u(2:4) = [c*v(1) - s*v(2), s*v(1) + c*v(2), v(3)]
    end subroutine

    subroutine closure()
        real(dp) :: reference, t, t0, mass, alpha, g, floor, chi, expected, dchi, &
            bars(7), h, yp, ym
        integer :: i
        reference = 4; t = 2; t0 = 1; mass = 0.5_dp; alpha = 0.3_dp; floor = 1.0e-8_dp
        do i = -8, 8
            g = 10.0_dp**i
            call vrmhd_chi_flux_limited(reference, t, t0, mass, alpha, g, floor, chi)
            expected = 1/(1/(reference*(t/t0)**2.5_dp) &
                + sqrt(g*g + floor*floor)/(alpha*t*sqrt(2*t/mass)))
            call require(abs(chi - expected) < 1.0e-13_dp*expected, &
                'harmonic-mean closure')
            call require(chi*g <= alpha*t*sqrt(2*t/mass), 'free-streaming flux cap')
        end do
        g = 0.7_dp; h = 1.0e-5_dp
        call vrmhd_chi_flux_limited_jvp(reference, t, t0, mass, alpha, g, floor, &
            0.0_dp, 0.2_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
            0.1_dp, 0.0_dp, dchi)
        call vrmhd_chi_flux_limited(reference, t + 0.2_dp*h, t0, mass, alpha, &
            g + 0.1_dp*h, floor, yp)
        call vrmhd_chi_flux_limited(reference, t - 0.2_dp*h, t0, mass, alpha, &
            g - 0.1_dp*h, floor, ym)
        call require(abs((yp - ym)/(2*h) - dchi) < 1.0e-9_dp, 'limited closure JVP')
        call vrmhd_chi_flux_limited_vjp(reference, t, t0, mass, alpha, g, floor, &
            1.0_dp, &
            bars(1), bars(2), bars(3), bars(4), &
            bars(5), bars(6), bars(7))
        call require(abs(0.2_dp*bars(2) + 0.1_dp*bars(6) - dchi) < 1.0e-13_dp, &
            'closure VJP')
    end subroutine
end program test_vrmhd_section
