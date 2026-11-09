program test_vrmhd_leaves
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use milca_vrmhd_leaves
    use milca_vrmhd_parameters, only: vrmhd_parameters_t
    use vrmhd_leaf_dispatch, only: nleaf, maxargs, nargs, names, evaluate, derivative
    implicit none
    character(32) :: gate
    integer, allocatable :: seed(:)
    integer :: n, seed_index
    real(dp) :: start, finish

    call random_seed(size=n)
    seed = [(104729 + 37*seed_index, seed_index=1, n)]
    call random_seed(put=seed)
    call get_command_argument(1, gate)
    call cpu_time(start)
    select case (gate)
    case ('derivatives')
        call derivatives()
    case ('physics')
        call physics()
    case ('energy')
        call energy()
    case ('covariance')
        call covariance()
    case default
        error stop 'choose derivatives, physics, energy or covariance'
    end select
    call cpu_time(finish)
    print '(a,a,a,f10.6)', 'PASS: ', trim(gate), '; CPU seconds = ', finish - start
contains
    subroutine close (actual, expected, label, tolerance)
        real(dp), intent(in) :: actual, expected
        character(*), intent(in) :: label
        real(dp), intent(in), optional :: tolerance
        real(dp) :: tol
        tol = 3.0e-13_dp
        if (present(tolerance)) tol = tolerance
        if (.not. abs(actual - expected) <= tol*max(1.0_dp, abs(expected))) then
            print *, trim(label), ': actual, expected = ', actual, expected
            error stop 'VRMHD independent oracle'
        end if
    end subroutine

    subroutine sample(x)
        real(dp), intent(out) :: x(:)
        call random_number(x)
        x = 0.5_dp + x
    end subroutine

    function cross(a, b) result(c)
        real(dp), intent(in) :: a(3), b(3)
        real(dp) :: c(3)
        c = [a(2)*b(3) - a(3)*b(2), a(3)*b(1) - a(1)*b(3), a(1)*b(2) - a(2)*b(1)]
    end function

    function stress(nu, nup, eps, b, g) result(tau)
        real(dp), intent(in) :: nu, nup, eps, b(3), g(3, 3)
        real(dp) :: tau(3, 3), s(3, 3), tensor(3, 3), identity(3, 3), bn(3)
        integer :: i
        identity = 0
        do i = 1, 3
            identity(i, i) = 1
        end do
        s = (g + transpose(g))/2 - identity*(g(1, 1) + g(2, 2) + g(3, 3))/3
        bn = b/sqrt(sum(b*b) + eps**2)
        tensor = spread(bn, 2, 3)*spread(bn, 1, 3) - identity*sum(bn*bn)/3
        tau = 2*nu*s + 3*nup*sum(tensor*s)*tensor
    end function

    function oracle(leaf, x) result(y)
        integer, intent(in) :: leaf
        real(dp), intent(in) :: x(maxargs)
        real(dp) :: y, tensor(3, 3), identity(3, 3), b(3), g(3, 3), h(3, 3)
        integer :: i, shift
        select case (leaf)
        case (1, 2, 3)
            y = x(1)*exp(x(4)*log(x(2)/x(3)))
        case (4)
            y = x(1)
        case (5)
            y = x(1)*(x(2) + x(3))/x(4)
        case (6)
            y = x(1)*x(2)
        case (7)
            y = dot_product(x(1)*x(6:8) + x(2)*x(3:5)/2, x(9:11))
        case (8)
            y = 1.5_dp*(x(1)*x(4) + x(2)*x(3))*x(6)/x(5)
        case (9)
            y = dot_product(x(1:3), x(4:6))
        case (10)
            y = x(1)*dot_product(x(2:4), x(5:7)) - dot_product(x(8:10), x(11:13))
        case (11)
            y = -x(1)*dot_product(x(2:4), x(5:7))
        case (12, 13)
            shift = 0
            if (leaf == 13) shift = 2
            identity = 0
            do i = 1, 3
                identity(i, i) = 1
            end do
            b = x(shift + 4:shift + 6) &
                /sqrt(sum(x(shift + 4:shift + 6)**2) + x(shift + 3)**2)
            tensor = x(shift + 2)*identity &
                + (x(shift + 1) - x(shift + 2))*spread(b, 2, 3)*spread(b, 1, 3)
            y = dot_product(x(shift + 10:shift + 12), &
                matmul(tensor, x(shift + 7:shift + 9)))
            if (leaf == 13) y = y*x(1)/x(2)
        case (14)
            g = reshape(x(8:16), [3, 3]); h = reshape(x(17:25), [3, 3])
            y = x(1)*(dot_product(matmul(g, x(2:4)), x(5:7)) &
                - dot_product(matmul(h, x(2:4)), x(2:4)))/2
        case (15)
            y = -x(1)*x(2)
        case (16, 17)
            y = dot_product(cross(x(1:3), x(4:6)), x(7:9))
        case (18)
            y = -x(1)*dot_product(x(2:4) - x(5:7), x(8:10))
        case (19)
            y = x(1)*dot_product(x(2:4), x(2:4) - x(5:7))*x(8)
        case (20, 21)
            g = reshape(x(7:15), [3, 3])
            tensor = stress(x(1), x(2), x(3), x(4:6), g)
            if (leaf == 20) then
                h = reshape(x(16:24), [3, 3])
                y = -sum(tensor*h)
            else
                y = sum(tensor*g)*x(16)
            end if
        case (22)
            y = -1.5_dp*x(1)*x(3)*dot_product(x(4:6), x(7:9))/x(2)
        case (23)
            y = x(1)*x(3)*x(4)*x(5)/x(2)
        case (24)
            y = x(1)*x(3)*(x(5) - x(4))*x(6)/x(2)
        case (25)
            y = x(1)*x(2)
        case (26)
            y = x(1)*dot_product(x(2:4), x(5:7))/2
        case (27)
            b = x(4:6)/sqrt(sum(x(4:6)**2) + x(3)**2)
            g = reshape(x(16:24), [3, 3]); h = reshape(x(25:33), [3, 3])
            y = -dot_product(x(2)*x(7:9) &
                + (x(1) - x(2))*b*dot_product(b, x(7:9)), &
                matmul(transpose(g), x(13:15)) &
                + matmul(transpose(h), x(10:12)))/2
        case default
            error stop 'missing independent primal oracle'
        end select
    end function

    subroutine derivatives()
        real(dp) :: x(maxargs), d(maxargs), grad(maxargs), weight, tangent
        real(dp) :: y, plus, minus, err(2), fine, scale, adjoint_error, max_fd
        real(dp), parameter :: step(2) = [0.004_dp, 0.002_dp]
        integer :: leaf, trial, k, coordinate
        max_fd = 0
        do leaf = 1, nleaf
            do trial = 1, 32
                call sample(x)
                call random_number(d)
                d = d - 0.5_dp
                d(nargs(leaf) + 1:) = 0
                call random_number(weight)
                weight = 2*weight - 1
                call evaluate(leaf, x, y)
                call close (y, oracle(leaf, x), trim(names(leaf))//' primal')
                call derivative(leaf, x, d, weight, tangent, grad)
                scale = max(1.0_dp, abs(tangent))
                do k = 1, 2
                    call evaluate(leaf, x + step(k)*d, plus)
                    call evaluate(leaf, x - step(k)*d, minus)
                    err(k) = abs(tangent - (plus - minus)/(2*step(k)))/scale
                end do
                if (err(2) > 0.3_dp*err(1) + 2.0e-11_dp) then
                    print *, trim(names(leaf)), err
                    error stop 'central difference must converge quadratically'
                end if
                call evaluate(leaf, x + 1.0e-5_dp*d, plus)
                call evaluate(leaf, x - 1.0e-5_dp*d, minus)
                fine = abs(tangent - (plus - minus)/2.0e-5_dp)/scale
                max_fd = max(max_fd, fine)
                call close (fine, 0.0_dp, trim(names(leaf))//' fine JVP', 2.0e-8_dp)
                adjoint_error = weight*tangent - dot_product(grad, d)
                scale = max(1.0_dp, abs(weight*tangent), sum(abs(grad*d)))
                call close (adjoint_error, 0.0_dp, &
                    trim(names(leaf))//' dual pairing', 2.0e-13_dp*scale)
            end do
            ! Coordinate directions prevent cancellation from hiding a missing partial.
            do coordinate = 1, nargs(leaf)
                d = 0; d(coordinate) = 1
                call derivative(leaf, x, d, 1.0_dp, tangent, grad)
                call evaluate(leaf, x + 1.0e-5_dp*d, plus)
                call evaluate(leaf, x - 1.0e-5_dp*d, minus)
                call close (tangent, (plus - minus)/2.0e-5_dp, &
                    trim(names(leaf))//' coordinate JVP', 2.0e-8_dp)
            end do
        end do
        print '(a,es12.4)', 'Maximum scaled centered-difference error: ', max_fd
    end subroutine

    subroutine physics()
        type(vrmhd_parameters_t) :: theta
        real(dp) :: x(maxargs), y, z, point(3), a(3), b(3), j(3), k(3), ck(3)
        real(dp) :: ampere_sum, left, right, volume, pi, temperature, value
        integer :: ix, iy, iz, trial
        integer, parameter :: points = 8
        theta%mi = 2.0_dp; theta%mu0 = 0.7_dp
        theta%eta0 = 0.13_dp; theta%t0 = 0.8_dp
        theta%nu0 = 0.27_dp; theta%chi_parallel0 = [1.7_dp, 2.3_dp]
        temperature = 4*theta%t0
        call vrmhd_eta(theta%eta0, temperature, theta%t0, theta%eta_exponent, value)
        call close (value, theta%eta0/8, 'Spitzer eta at four reference temperatures')
        call vrmhd_nu(theta%nu0, temperature, theta%t0, theta%nu_exponent, value)
        call close (value, theta%nu0/8, 'viscosity closure')
        call vrmhd_chi_parallel(theta%chi_parallel0(1), temperature, theta%t0, &
            theta%chi_parallel_exponent, value)
        call close (value, 32*theta%chi_parallel0(1), 'parallel chi closure')
        ! Prescribed constant coefficients, as in the JOREK comparisons.
        call vrmhd_chi_parallel(theta%chi_parallel0(1), temperature, theta%t0, &
            0.0_dp, value)
        call close (value, theta%chi_parallel0(1), 'temperature-independent chi')
        call vrmhd_nu(theta%nu0, temperature, theta%t0, 0.0_dp, value)
        call close (value, theta%nu0, 'temperature-independent viscosity')
        pi = acos(-1.0_dp); volume = (2*pi)**3
        ampere_sum = 0; left = 0; right = 0
        ! A=(sin z,sin x,sin y); curl curl A=A, independently differentiated.
        do iz = 0, points - 1
            do iy = 0, points - 1
                do ix = 0, points - 1
                    point = 2*pi*real([ix, iy, iz], dp)/points
                    a = sin(point([3, 1, 2]))
                    b = cos(point([2, 3, 1]))
                    j = a/theta%mu0
                    k = [a(1) + 0.2_dp*b(1), 0.7_dp*a(2) + 0.3_dp*b(2), &
                        0.5_dp*a(3) + 0.4_dp*b(3)]
                    ck = [0.5_dp*b(1) + 0.3_dp*a(1), b(2) + 0.4_dp*a(2), &
                        0.7_dp*b(3) + 0.2_dp*a(3)]
                    x = 0
                    x(:13) = [theta%mu0, j, k, b, ck]
                    call evaluate(10, x, y)
                    ampere_sum = ampere_sum + y
                    x(8:10) = 0
                    call evaluate(10, x, y)
                    left = left + y
                    x(2:4) = 0; x(8:10) = b
                    call evaluate(10, x, y)
                    right = right - y
                end do
            end do
        end do
        call close (ampere_sum/points**3, 0.0_dp, 'periodic weak Ampere')
        call close (left*volume/points**3, 1.1_dp*volume, &
            'analytic Ampere mass integral')
        call close (right*volume/points**3, 1.1_dp*volume, &
            'analytic Ampere curl integral')
        do trial = 1, 64
            call sample(x)
            b = x(4:6)
            x(2) = 0
            x(7:9) = cross(b, [0.3_dp, -0.7_dp, 1.1_dp])
            call evaluate(12, x, y)
            call close (y, 0.0_dp, 'parallel conduction, B dot grad T = 0')
            x(1) = 0; x(2) = 0.4_dp
            call evaluate(12, x, y)
            call close (y, 0.4_dp*dot_product(x(7:9), x(10:12)), &
                'perpendicular conduction for transverse gradient')
            x(4:6) = 0
            call evaluate(12, x, y)
            call close (y, 0.4_dp*dot_product(x(7:9), x(10:12)), &
                'regularized null field')
            call sample(x)
            x(2) = 0
            x(7:15) = [0.2_dp, 0.4_dp, -0.7_dp, -0.4_dp, 0.2_dp, 0.5_dp, &
                0.7_dp, -0.5_dp, 0.2_dp]
            call evaluate(21, x, y)
            call close (y, 0.0_dp, 'isotropic dilation and rigid rotation do not heat')
            call sample(x)
            call evaluate(21, x, y)
            if (y < 0) error stop 'nonnegative viscous dissipation'
            ! Thermal pressure work equals minus kinetic pressure work for w=v.
            call vrmhd_pressure(x(1), x(2), x(3), theta%mi, y)
            call vrmhd_pressure_work(y, x(4), z)
            call vrmhd_compression(x(1), theta%mi, x(2), x(4), 1.0_dp, value)
            call vrmhd_compression(x(1), theta%mi, x(3), x(4), 1.0_dp, y)
            call close (z + value + y, 0.0_dp, 'two-temperature pressure exchange')
        end do
    end subroutine

    subroutine energy()
        real(dp) :: x(maxargs), v(3), b(3), j(3), js(3), eta, lorentz, ideal, loss, heat
        real(dp) :: g(3, 3), strain(3, 3), identity(3, 3), projector(3, 3)
        real(dp) :: bn(3), expected, rho, rd, vd(3), mi, kinetic, density, momentum
        real(dp) :: advection, thermal
        integer :: trial, driven
        identity = 0
        do trial = 1, 3
            identity(trial, trial) = 1
        end do
        do trial = 1, 256
            call sample(x)
            v = 2*x(1:3) - 2; b = 2*x(4:6) - 2; j = 2*x(7:9) - 2
            eta = x(10)
            x(:9) = [j, b, v]
            call evaluate(16, x, lorentz)
            x(:9) = [v, b, j]
            call evaluate(17, x, ideal)
            call close (lorentz + ideal, 0.0_dp, 'Lorentz and magnetic work cancel')
            do driven = 0, 1
                js = driven*[0.2_dp, -0.4_dp, 0.7_dp]
                x(:10) = [eta, j, js, j]
                call evaluate(18, x, loss)
                x(:8) = [eta, j, js, 1.0_dp]
                call evaluate(19, x, heat)
                call close (heat, eta*dot_product(j, j - js), 'analytic Ohmic power')
                call close (loss + heat, 0.0_dp, &
                    'resistive loss equals electron heating')
            end do
            call sample(x)
            g = reshape(2*x(7:15) - 2, [3, 3])
            strain = (g + transpose(g))/2 &
                - identity*(g(1, 1) + g(2, 2) + g(3, 3))/3
            bn = b/sqrt(sum(b*b) + x(3)**2)
            projector = spread(bn, 2, 3)*spread(bn, 1, 3) &
                - identity*sum(bn*bn)/3
            expected = 2*x(1)*sum(strain*strain) &
                + 3*x(2)*sum(projector*strain)**2
            x(4:6) = b; x(7:15) = reshape(g, [9])
            x(16:24) = x(7:15)
            call evaluate(20, x, loss)
            x(16) = 1
            call evaluate(21, x, heat)
            call close (heat, expected, 'positive strain-norm viscous heating')
            call close (loss + heat, 0.0_dp, 'viscous force and heating cancel')

            rho = x(1); rd = x(2); vd = x(4:6); mi = x(3)
            kinetic = sum(v*v)/2
            x(:11) = [rho, rd, v, vd, v]
            call evaluate(7, x, momentum)
            x(:2) = [rd, kinetic]
            call evaluate(6, x, density)
            call close (density, rd*kinetic, 'density mass pairing')
            call close (momentum, &
                rho*dot_product(v, vd) + rd*kinetic, 'kinetic mass chain rule')
            x(:25) = [rho, v, v, reshape(g, [9]), reshape(g, [9])]
            call evaluate(14, x, advection)
            call close (advection, 0.0_dp, 'skew advection does no work')
            x(:9) = [rho, mi, 1.2_dp, v, 0.0_dp, 0.0_dp, 0.0_dp]
            call evaluate(22, x, thermal)
            call close (thermal, 0.0_dp, 'constant thermal test has no advection')
            x(:14) = [rho, mi, 1.7_dp, 0.3_dp, 0.1_dp, b, vd, &
                0.0_dp, 0.0_dp, 0.0_dp]
            call evaluate(13, x, thermal)
            call close (thermal, 0.0_dp, 'constant thermal test has no conduction')
        end do
    end subroutine

    subroutine covariance()
        real(dp) :: x(maxargs), rotated(maxargs), rotation(3, 3), a(3), b(3)
        real(dp) :: g(3, 3), y, yr
        integer :: trial, leaf, first, last, shift
        do trial = 1, 64
            call sample(a); call sample(b)
            a = a/sqrt(sum(a*a))
            b = b - dot_product(a, b)*a
            b = b/sqrt(sum(b*b))
            rotation(:, 1) = a; rotation(:, 2) = b; rotation(:, 3) = cross(a, b)
            do leaf = 1, nleaf
                call sample(x)
                rotated = x
                select case (leaf)
                case (7)
                    first = 3; last = 9
                case (9, 16, 17)
                    first = 1; last = 7
                case (10)
                    first = 2; last = 11
                case (11)
                    first = 2; last = 5
                case (12, 13)
                    shift = 0
                    if (leaf == 13) shift = 2
                    first = 4 + shift; last = 10 + shift
                case (14)
                    first = 2; last = 5
                    g = matmul(rotation, &
                        matmul(reshape(x(8:16), [3, 3]), transpose(rotation)))
                    rotated(8:16) = reshape(g, [9])
                    g = matmul(rotation, &
                        matmul(reshape(x(17:25), [3, 3]), transpose(rotation)))
                    rotated(17:25) = reshape(g, [9])
                case (26)
                    first = 2; last = 5
                case (27)
                    first = 4; last = 13
                    g = matmul(rotation, &
                        matmul(reshape(x(16:24), [3, 3]), transpose(rotation)))
                    rotated(16:24) = reshape(g, [9])
                    g = matmul(rotation, &
                        matmul(reshape(x(25:33), [3, 3]), transpose(rotation)))
                    rotated(25:33) = reshape(g, [9])
                case (18)
                    first = 2; last = 8
                case (19)
                    first = 2; last = 5
                case (20, 21)
                    first = 4; last = 4
                    g = matmul(rotation, &
                        matmul(reshape(x(7:15), [3, 3]), transpose(rotation)))
                    rotated(7:15) = reshape(g, [9])
                    if (leaf == 20) then
                        g = reshape(x(16:24), [3, 3])
                        g = matmul(rotation, matmul(g, transpose(rotation)))
                        rotated(16:24) = reshape(g, [9])
                    end if
                case (22)
                    first = 4; last = 7
                case default
                    cycle
                end select
                do shift = first, last, 3
                    rotated(shift:shift + 2) = matmul(rotation, x(shift:shift + 2))
                end do
                call evaluate(leaf, x, y)
                call evaluate(leaf, rotated, yr)
                call close (yr, y, trim(names(leaf))//' proper rotation')
            end do
        end do
    end subroutine
end program test_vrmhd_leaves
