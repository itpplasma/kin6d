module nn_mapped_fixture
    ! Closed-form vacuum oracle: B = F grad phi + B_v e_Z + sum_m grad u_m,
    ! u_m = c_m R0 (R/R0)^m cos(m phi) = Re of a holomorphic function of
    ! x + i y, so curl B = 0 and div B = 0 everywhere. Potential
    ! A = (0, B_v R/2, -F ln R + sum_m c_m R0 (R/R0)^m sin(m phi)).
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_equilibrium_chart, only: ideal_chart_t, curve_chart_t
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp), r0 = 3, minor = 1, ftor = 3
    real(dp), parameter :: bv = 0.1_dp, coef(3) = [0.05_dp, 0.03_dp, 0.02_dp]
    ! Rotating ellipse: semi-axes e1, e2 (equal area), turning with angle phi
    ! (the cross-section has period pi, two field periods).
    real(dp), parameter :: e1 = 1.25_dp, e2 = 0.8_dp
    ! Independent analytic test charts. Production geometry remains generated.
    real(dp), parameter :: shear(3, 3) = reshape([1.0_dp, 0.0_dp, 0.0_dp, &
        0.2_dp, 1.1_dp, 0.0_dp, 0.1_dp, 0.15_dp, 0.9_dp], [3, 3])
    type, extends(ideal_chart_t) :: analytic_chart_t
        real(dp) :: amplitude = 0.08_dp, orientation = 1.0_dp
        logical :: invalid_position = .false.
    contains
        procedure :: jet => analytic_jet
    end type
    type, extends(ideal_chart_t) :: wrapped_curve_t
        type(curve_chart_t) :: curve
    contains
        procedure :: jet => wrapped_jet
    end type
contains
    pure subroutine analytic_jet(chart, q, x, d1, d2, d3)
        class(analytic_chart_t), intent(in) :: chart
        real(dp), intent(in) :: q(3)
        real(dp), intent(out) :: x(3), d1(3, 3), d2(3, 3, 3), d3(3, 3, 3, 3)
        real(dp) :: a, c, sn, v(3), first(3, 3), second(3, 3, 3)
        real(dp) :: third(3, 3, 3, 3)
        integer :: i, j, k
        a = chart%amplitude; c = cos(q(3)); sn = sin(q(3))
        v = [q(1)*c, q(1)*sn, q(2) + a*sin(3*q(3))]
        first(:, 1) = [c, sn, 0.0_dp]
        first(:, 2) = [0.0_dp, 0.0_dp, 1.0_dp]
        first(:, 3) = [-q(1)*sn, q(1)*c, 3*a*cos(3*q(3))]
        second = 0; third = 0
        second(:, 1, 3) = [-sn, c, 0.0_dp]
        second(:, 3, 1) = second(:, 1, 3)
        second(:, 3, 3) = [-q(1)*c, -q(1)*sn, -9*a*sin(3*q(3))]
        third(:, 1, 3, 3) = [-c, -sn, 0.0_dp]
        third(:, 3, 1, 3) = third(:, 1, 3, 3)
        third(:, 3, 3, 1) = third(:, 1, 3, 3)
        third(:, 3, 3, 3) = [q(1)*sn, -q(1)*c, -27*a*cos(3*q(3))]
        x = matmul(shear, v); d1 = matmul(shear, first)
        do i = 1, 3
            do j = 1, 3
                d2(:, i, j) = matmul(shear, second(:, i, j))
                do k = 1, 3
                    d3(:, i, j, k) = matmul(shear, third(:, i, j, k))
                end do
            end do
        end do
        x(1) = chart%orientation*x(1)
        d1(1, :) = chart%orientation*d1(1, :)
        d2(1, :, :) = chart%orientation*d2(1, :, :)
        d3(1, :, :, :) = chart%orientation*d3(1, :, :, :)
        if (chart%invalid_position) x(1) = ieee_value(0.0_dp, ieee_quiet_nan)
    end subroutine

    pure subroutine wrapped_jet(chart, q, x, d1, d2, d3)
        class(wrapped_curve_t), intent(in) :: chart
        real(dp), intent(in) :: q(3)
        real(dp), intent(out) :: x(3), d1(3, 3), d2(3, 3, 3), d3(3, 3, 3, 3)
        integer, parameter :: p(3) = [2, 1, 3]
        real(dp) :: a(3, 3), b(3, 3, 3), c(3, 3, 3, 3)
        integer :: i, j, k
        call chart%curve%jet([q(2), q(1) - r0, q(3)], x, a, b, c)
        do i = 1, 3
            d1(:, i) = a(:, p(i))
            do j = 1, 3
                d2(:, i, j) = b(:, p(i), p(j))
                do k = 1, 3
                    d3(:, i, j, k) = c(:, p(i), p(j), p(k))
                end do
            end do
        end do
    end subroutine

    subroutine expected_point(kind, q, x, d, jac)
        ! Closed Cartesian geometry, independent of chart%jet and Piola kernels.
        integer, intent(in) :: kind
        real(dp), intent(in) :: q(3)
        real(dp), intent(out) :: x(3), d(3, 3), jac
        real(dp) :: er(3), ep(3), ez(3), n(3), dn(3), amp, g, gp, speed
        er = [cos(q(3)), sin(q(3)), 0.0_dp]
        ep = [-sin(q(3)), cos(q(3)), 0.0_dp]; ez = [0.0_dp, 0.0_dp, 1.0_dp]
        if (kind <= 2) then
            x = q(1)*er + q(2)*ez
            d(:, 1) = er; d(:, 2) = ep; d(:, 3) = ez
            jac = 1
            if (kind == 2) then
                x = matmul(shear, x + 0.08_dp*sin(3*q(3))*ez)
                d(:, 2) = ep + (0.24_dp*cos(3*q(3))/q(1))*ez
                d = matmul(shear, d)
                jac = 0.99_dp
            end if
        else
            amp = 0
            if (kind == 4) amp = 0.08_dp
            g = 3*amp*cos(3*q(3)); gp = -9*amp*sin(3*q(3))
            speed = sqrt(r0**2 + g**2)
            n = -(g/speed)*ep + (r0/speed)*ez
            dn = (g/speed)*er - (r0**2*gp/speed**3)*ep &
                - (r0*g*gp/speed**3)*ez
            x = q(1)*er + q(2)*n + amp*sin(3*q(3))*ez
            d(:, 1) = er
            d(:, 2) = ep + (g*ez + q(2)*dn)/q(1)
            d(:, 3) = n
            jac = (r0*q(1) + g**2)/(q(1)*speed) &
                - q(2)*r0*gp/(q(1)*speed**2)
        end if
    end subroutine

    pure function cross(a, b) result(c)
        real(dp), intent(in) :: a(3), b(3)
        real(dp) :: c(3)
        c = [a(2)*b(3) - a(3)*b(2), a(3)*b(1) - a(1)*b(3), &
             a(1)*b(2) - a(2)*b(1)]
    end function

    subroutine cartesian_potential(x, a)
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: a(3)
        complex(dp) :: z
        integer :: m
        z = cmplx(x(1), x(2), dp)
        a = [-bv*x(2)/2, bv*x(1)/2, -ftor*log(sum(x(:2)**2))/2]
        do m = 1, size(coef)
            a(3) = a(3) + coef(m)*r0**(1 - m)*aimag(z**m)
        end do
    end subroutine

    subroutine cartesian_field(x, b)
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: b(3)
        complex(dp) :: z
        real(dp) :: c
        integer :: m
        z = cmplx(x(1), x(2), dp)
        b = [-ftor*x(2)/sum(x(:2)**2), ftor*x(1)/sum(x(:2)**2), bv]
        do m = 1, size(coef)
            c = coef(m)*m*r0**(1 - m)
            b(1) = b(1) + c*real(z**(m - 1), dp)
            b(2) = b(2) - c*aimag(z**(m - 1))
        end do
    end subroutine

    subroutine constant_potential(x, a)
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: a(3)
        a = [-bv*x(2)/2, bv*x(1)/2, 0.0_dp]
    end subroutine

    subroutine constant_field(x, b)
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: b(3)
        b = [0.0_dp, 0.0_dp, bv] + 0*x
    end subroutine

    subroutine circle(u, theta, x)
        real(dp), intent(in) :: u, theta
        real(dp), intent(out) :: x(2)
        x = [r0 + minor*u*cos(theta), minor*u*sin(theta)]
    end subroutine

    subroutine identity_map(x, phi, y)
        real(dp), intent(in) :: x(2), phi
        real(dp), intent(out) :: y(2)
        y = x + 0*phi
    end subroutine

    subroutine ellipse_map(x, phi, y)
        real(dp), intent(in) :: x(2), phi
        real(dp), intent(out) :: y(2)
        real(dp) :: q(2, 2), d(2)
        q = reshape([cos(phi), sin(phi), -sin(phi), cos(phi)], [2, 2])
        d = matmul(transpose(q), x - [r0, 0.0_dp])
        y = [r0, 0.0_dp] + matmul(q, [e1*d(1), e2*d(2)])
    end subroutine

    subroutine potential(x, a)
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: a(3)
        integer :: m
        a = [0.0_dp, bv*x(1)/2, -ftor*log(x(1))]
        do m = 1, size(coef)
            a(3) = a(3) + coef(m)*r0*(x(1)/r0)**m*sin(m*x(3))
        end do
    end subroutine

    subroutine exact_field(x, b)
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: b(3)
        integer :: m
        b = [0.0_dp, ftor/x(1), bv]
        do m = 1, size(coef)
            b(1) = b(1) + coef(m)*m*(x(1)/r0)**(m - 1)*cos(m*x(3))
            b(2) = b(2) - coef(m)*m*(x(1)/r0)**(m - 1)*sin(m*x(3))
        end do
    end subroutine
end module nn_mapped_fixture

program test_nonnested_mapped
    ! Usage: identity | hladder N nr... | modes nr N... | zero_forcing
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use nn_mapped_fixture
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use kin6d_equilibrium_chart, only: cylindrical_chart_t
    use fortfem_meridian_feec_n, only: meridian_qpoint, meridian_qgrad
    use kin6d_mhd_mesh, only: flux_mesh_t, flux_mesh_build
    use kin6d_equilibrium_nonnested, only: free_b_field_t, nn_report_t, nn_build, &
        nn_solve
    use kin6d_nonnested_mapped, only: nn_map_t, nn_map_build, nn_mapped_data, &
        nn_mapped_solve, nn_mapped_error, nn_mapped_report_t, nn_map_set_chart, &
        nn_map_point, map_metric_t, nn_mapped_cartesian_data, &
        nn_mapped_cartesian_error, nn_vector_field_i
    implicit none
    character(32) :: mode, arg, env
    integer, allocatable :: list(:)
    integer :: i, fixed_arg, st
    real(dp) :: min_order = 1.7_dp
    real(dp), allocatable :: errs(:)
    call get_command_argument(1, mode)
    if (command_argument_count() > 1) then
        call get_command_argument(2, arg)
        read (arg, *) fixed_arg
        allocate (list(command_argument_count() - 2))
        do i = 1, size(list)
            call get_command_argument(i + 2, arg)
            read (arg, *) list(i)
        end do
    end if
    select case (trim(mode))
    case ('point_oracles')
        call point_oracles()
    case ('cartesian_identity')
        call cartesian_identity()
    case ('cartesian_quadrature')
        call cartesian_quadrature()
    case ('cartesian_vacuum')
        call cartesian_vacuum()
    case ('identity')
        call identity_regression()
    case ('hladder')
        allocate (errs(size(list)))
        do i = 1, size(list)
            call ellipse_case(list(i), fixed_arg, errs(i))
            if (i > 1) then
                write (*, '(a,i0,a,f6.3)') 'order nr=', list(i), ': ', &
                    log(errs(i - 1)/errs(i))/log(real(list(i), dp)/list(i - 1))
            end if
        end do
        if (size(errs) > 1) then
            call get_environment_variable('KIN6D_MAPPED_MIN_ORDER', env, status=st)
            if (st == 0) read (env, *) min_order
            if (log(errs(size(errs) - 1)/errs(size(errs))) &
                /log(real(list(size(list)), dp)/list(size(list) - 1)) < min_order) &
                error stop 'mapped h order below the gate'
        end if
    case ('modes')
        allocate (errs(size(list)))
        do i = 1, size(list)
            call ellipse_case(fixed_arg, list(i), errs(i))
        end do
        if (size(errs) > 1) then
            if (errs(size(errs)) > 0.5_dp*errs(1)) error stop 'mode ladder stalls'
        end if
    case ('zero_forcing')
        call zero_forcing_residual()
    case default
        error stop 'usage: identity | hladder | modes | point_oracles | cartesian_* | &
            &zero_forcing'
    end select
contains
    subroutine small(label, value, gate)
        character(*), intent(in) :: label
        real(dp), intent(in) :: value, gate
        if (.not. ieee_is_finite(value)) then
            write (*, '(a,a)') 'nonfinite oracle: ', label
            error stop 'mapped Cartesian oracle is nonfinite'
        end if
        if (value > gate) then
            write (*, '(a,2es18.9)') label, value, gate
            error stop 'mapped Cartesian oracle exceeds gate'
        end if
    end subroutine

    subroutine point_oracles()
        type(free_b_field_t), target :: field
        type(nn_map_t), target :: map
        type(cylindrical_chart_t) :: cylinder
        type(analytic_chart_t) :: analytic
        type(wrapped_curve_t) :: wrapper
        real(dp) :: ac(3, 0:3), as(3, 0:3), amp
        integer :: kind, st
        call make_field(1, 3, field)
        do kind = 0, 4
            select case (kind)
            case (0)
                call nn_map_build(field%fs, 0, 1, identity_map, map, st)
                if (st /= 0) error stop 'legacy identity map'
            case (1)
                call nn_map_set_chart(cylinder, map)
            case (2)
                call nn_map_set_chart(analytic, map)
            case (3, 4)
                amp = 0
                if (kind == 4) amp = 0.08_dp
                ac = 0; as = 0
                ac(1, 1) = r0; as(2, 1) = r0; as(3, 3) = amp
                call wrapper%curve%init(ac, as, 3, st)
                if (st /= 0) error stop 'circle/closed-helix chart init'
                call nn_map_set_chart(wrapper, map)
            end select
            call point_case(field, map, kind)
        end do
        call vacuum_cartesian_oracle()
        call physical_piola_fd()
    end subroutine

    subroutine point_case(field, map, kind)
        type(free_b_field_t), intent(in), target :: field
        type(nn_map_t), intent(in), target :: map
        integer, intent(in) :: kind
        real(dp), parameter :: angles(4) = [0.37_dp, 1.1_dp, 2.7_dp, 5.8_dp]
        real(dp), parameter :: steps(3) = [1.0e-4_dp, 5.0e-5_dp, 2.5e-5_dp]
        type(map_metric_t) :: metric
        type(analytic_chart_t) :: bad
        type(nn_map_t) :: badmap
        real(dp) :: x(2), val(field%fs%nl0), grad(2, field%fs%nl0), w
        real(dp) :: pos(3), want(3), d(3, 3), p1(3, 3), p2(3, 3), jac, jref
        real(dp) :: cof(3, 3), expected1(3, 3), expected2(3, 3), points(3, 4)
        real(dp) :: k(3, 3, 4), km(3, 3), a(3), b(3), ap(3), bp(3)
        real(dp) :: plus(3), minus(3), dummy1(3, 3), dummy2(3, 3), jd, xp(2), angle
        real(dp) :: fd(3, 3), errors(3), h, rotated(3), turn
        integer :: st, j, space, l, c
        call meridian_qgrad(field%fs, 1, 1, x, w, val, grad, st)
        if (st /= 0) error stop 'point oracle quadrature'
        do j = 1, size(angles)
            points(:, j) = [x, angles(j)]
        end do
        metric%fs => field%fs; metric%map => map
        do space = 1, 2
            metric%space = space
            call metric%tensor(1, points, k, st)
            if (st /= 0) error stop 'Cartesian metric oracle'
            do j = 1, size(angles)
                call expected_point(kind, points(:, j), want, d, jref)
                call nn_map_point(map, field%fs%sp(0)%dof(:, 1), val, grad, x, &
                    angles(j), pos, p1, p2, jac, st)
                if (st /= 0) error stop 'Cartesian point oracle'
                cof(:, 1) = cross(d(:, 2), d(:, 3))
                cof(:, 2) = cross(d(:, 3), d(:, 1))
                cof(:, 3) = cross(d(:, 1), d(:, 2))
                expected1 = cof/jref; expected2 = d/jref
                call small('position', maxval(abs(pos - want)), 2.0e-12_dp)
                call small('positive analytic Jacobian', abs(jac - jref), 2.0e-12_dp)
                call small('covariant Piola', maxval(abs(p1 - expected1)), 2.0e-12_dp)
                call small('contravariant Piola', maxval(abs(p2 - expected2)), &
                    2.0e-12_dp)
                if (space == 1) then
                    km = jref*matmul(transpose(expected1), expected1)
                else
                    km = jref*matmul(transpose(expected2), expected2)
                end if
                call small('pulled-back Hodge', maxval(abs(k(:, :, j) - km)), &
                    5.0e-12_dp)
                call cartesian_potential(want, a); call cartesian_field(want, b)
                ap = matmul(transpose(d), a); bp = matmul(transpose(cof), b)
                call small('1-form reconstruction', maxval(abs(matmul(p1, ap) - a)), &
                    5.0e-12_dp)
                call small('2-form flux reconstruction', &
                    maxval(abs(matmul(p2, bp) - b)), 5.0e-12_dp)
                call nn_map_point(map, field%fs%sp(0)%dof(:, 1), val, grad, x, &
                    angles(j) + 2*pi, plus, dummy1, dummy2, jd, st)
                if (st /= 0) error stop 'full-turn point oracle'
                call small('full-turn position', maxval(abs(plus - pos)), 2.0e-12_dp)
                call small('full-turn P1', maxval(abs(dummy1 - p1)), 2.0e-12_dp)
                call small('full-turn P2', maxval(abs(dummy2 - p2)), 2.0e-12_dp)
                if (kind >= 3) then
                    turn = 2*pi/3
                    rotated = [cos(turn)*pos(1) - sin(turn)*pos(2), &
                        sin(turn)*pos(1) + cos(turn)*pos(2), pos(3)]
                    call nn_map_point(map, field%fs%sp(0)%dof(:, 1), val, grad, &
                        x, angles(j) + turn, plus, dummy1, dummy2, jd, st)
                    if (st /= 0) error stop 'three-period curve point'
                    call small('three-period Cartesian equivariance', &
                        maxval(abs(plus - rotated)), 2.0e-12_dp)
                end if
                if (space == 2 .and. j == 1) then
                    do l = 1, size(steps)
                        h = steps(l)
                        do c = 1, 3
                            xp = x; angle = angles(j)
                            if (c == 1) xp(1) = x(1) + h
                            if (c == 2) angle = angles(j) + h
                            if (c == 3) xp(2) = x(2) + h
                            call nn_map_point(map, field%fs%sp(0)%dof(:, 1), &
                                val, grad, xp, angle, plus, dummy1, dummy2, jd, st)
                            if (st /= 0) error stop 'positive point FD'
                            xp = x; angle = angles(j)
                            if (c == 1) xp(1) = x(1) - h
                            if (c == 2) angle = angles(j) - h
                            if (c == 3) xp(2) = x(2) - h
                            call nn_map_point(map, field%fs%sp(0)%dof(:, 1), &
                                val, grad, xp, angle, minus, dummy1, dummy2, jd, st)
                            if (st /= 0) error stop 'negative point FD'
                            fd(:, c) = (plus - minus)/(2*h)
                            if (c == 2) fd(:, c) = fd(:, c)/x(1)
                        end do
                        errors(l) = maxval(abs(fd - d))/max(1.0_dp, maxval(abs(d)))
                    end do
                    call small('position FD finest', errors(3), 2.0e-8_dp)
                end if
            end do
        end do
        write (*, '(a,i0,a,3es13.5)') 'point kind=', kind, ' FD ladder=', errors
        if (kind >= 3) call integrated_period_trace(field, map, kind)
        if (kind == 2) then
            bad%orientation = -1
            call nn_map_set_chart(bad, badmap)
            call nn_map_point(badmap, field%fs%sp(0)%dof(:, 1), val, grad, x, &
                angles(1), pos, p1, p2, jac, st)
            if (st == 0) error stop 'negative orientation accepted as abs(detD)'
            bad%orientation = 1; bad%invalid_position = .true.
            call nn_map_set_chart(bad, badmap)
            call nn_map_point(badmap, field%fs%sp(0)%dof(:, 1), val, grad, x, &
                angles(1), pos, p1, p2, jac, st)
            if (st == 0) error stop 'nonfinite physical chart point accepted'
        end if
    end subroutine

    subroutine integrated_period_trace(field, map, kind)
        type(free_b_field_t), intent(in) :: field
        type(nn_map_t), intent(in) :: map
        integer, intent(in) :: kind
        integer, parameter :: nloop = 37, ntheta = 17, nphi = 19
        real(dp) :: x(2), refx(2), w, val(field%fs%nl0), grad(2, field%fs%nl0)
        real(dp) :: pos(3), want(3), d(3, 3), p1(3, 3), p2(3, 3), jac, jr
        real(dp) :: angle, theta, tangent(3), area(3), bref(3), pairing
        real(dp) :: period, flux, expected_flux
        integer :: i, j, st
        call meridian_qgrad(field%fs, 1, 1, refx, w, val, grad, st)
        if (st /= 0) error stop 'period/trace reference basis'
        x = [r0 + 0.25_dp, 0.1_dp]; period = 0
        do j = 1, nloop
            angle = 2*pi*(j - 0.27_dp)/nloop
            call expected_point(kind, [x, angle], want, d, jr)
            call nn_map_point(map, field%fs%sp(0)%dof(:, 1), val, grad, x, &
                angle, pos, p1, p2, jac, st)
            if (st /= 0) error stop 'physical loop-period point'
            period = period + dot_product(p1(:, 2)/x(1), x(1)*d(:, 2))
        end do
        period = (2*pi/nloop)*period
        call small('integrated physical cohomology period', &
            abs(period - 2*pi), 2.0e-11_dp)
        flux = 0
        do i = 1, ntheta
            theta = 2*pi*(i - 0.31_dp)/ntheta
            x = [r0 + minor*cos(theta), minor*sin(theta)]
            tangent = minor*[-sin(theta), 0.0_dp, cos(theta)]
            ! Manufactured reference flux vector, not a vacuum solution.
            bref = [cos(theta), 0.2_dp, sin(theta)]
            do j = 1, nphi
                angle = 2*pi*(j - 0.27_dp)/nphi
                call expected_point(kind, [x, angle], want, d, jr)
                area = cross(x(1)*d(:, 2), matmul(d, tangent))
                call nn_map_point(map, field%fs%sp(0)%dof(:, 1), val, grad, x, &
                    angle, pos, p1, p2, jac, st)
                if (st /= 0) error stop 'physical boundary-trace point'
                pairing = dot_product(matmul(p2, bref), area)
                call small('independent outward flux pairing', &
                    abs(pairing - x(1)*minor), 2.0e-11_dp)
                flux = flux + pairing
            end do
        end do
        flux = flux*(2*pi)**2/(ntheta*nphi)
        expected_flux = (2*pi)**2*r0*minor
        call small('integrated physical boundary flux', &
            abs(flux - expected_flux)/expected_flux, 2.0e-12_dp)
        write (*, '(a,i0,a,2es18.9)') 'chart kind=', kind, &
            ' integrated period/normal flux=', period, flux
    end subroutine

    subroutine vacuum_cartesian_oracle()
        real(dp), parameter :: h(3) = [1.0e-4_dp, 5.0e-5_dp, 2.5e-5_dp]
        real(dp) :: x(3), xp(3), ap(3), am(3), bp(3), bm(3), exact(3)
        real(dp) :: da(3, 3), db(3, 3), curl(3), errors(3), div
        integer :: l, c
        x = [2.7_dp, -0.4_dp, 0.3_dp]
        call cartesian_field(x, exact)
        do l = 1, 3
            do c = 1, 3
                xp = x; xp(c) = x(c) + h(l)
                call cartesian_potential(xp, ap); call cartesian_field(xp, bp)
                xp(c) = x(c) - h(l)
                call cartesian_potential(xp, am); call cartesian_field(xp, bm)
                da(:, c) = (ap - am)/(2*h(l)); db(:, c) = (bp - bm)/(2*h(l))
            end do
            curl = [da(3, 2) - da(2, 3), da(1, 3) - da(3, 1), &
                da(2, 1) - da(1, 2)]
            errors(l) = maxval(abs(curl - exact))
            div = abs(db(1, 1) + db(2, 2) + db(3, 3))
            call small('Cartesian vacuum divergence', div, 1.0e-7_dp)
            curl = [db(3, 2) - db(2, 3), db(1, 3) - db(3, 1), &
                db(2, 1) - db(1, 2)]
            call small('Cartesian vacuum curl', maxval(abs(curl)), 1.0e-7_dp)
        end do
        call small('Cartesian potential curl finest', errors(3), 1.0e-7_dp)
        write (*, '(a,3es13.5)') 'Cartesian curl(A)-B FD ladder=', errors
    end subroutine

    subroutine physical_piola_fd()
        ! Curl/divergence stencils perturb physical Cartesian coordinates.
        real(dp), parameter :: steps(3) = [2.0e-4_dp, 1.0e-4_dp, 5.0e-5_dp]
        type(free_b_field_t) :: field
        type(nn_map_t) :: map
        type(analytic_chart_t) :: chart
        real(dp) :: x(3), d(3, 3), jac, xp(3), ap(3), am(3), bp(3), bm(3)
        real(dp) :: da(3, 3), db(3, 3), curl(3), target(3), errors(2, 3)
        integer :: l, c
        call make_field(1, 1, field)
        call nn_map_set_chart(chart, map)
        call expected_point(2, [3.2_dp, 0.17_dp, 0.37_dp], x, d, jac)
        target = (d(:, 2) + d(:, 3))/jac
        do l = 1, 3
            do c = 1, 3
                xp = x; xp(c) = x(c) + steps(l)
                call physical_piola_value(field, map, xp, ap, bp)
                xp(c) = x(c) - steps(l)
                call physical_piola_value(field, map, xp, am, bm)
                da(:, c) = (ap - am)/(2*steps(l))
                db(:, c) = (bp - bm)/(2*steps(l))
            end do
            curl = [da(3, 2) - da(2, 3), da(1, 3) - da(3, 1), &
                da(2, 1) - da(1, 2)]
            errors(1, l) = maxval(abs(curl - target))
            errors(2, l) = abs(db(1, 1) + db(2, 2) + db(3, 3) - 2/jac)
        end do
        call small('physical Piola curl finest', errors(1, 3), 2.0e-7_dp)
        call small('physical Piola divergence finest', errors(2, 3), 2.0e-7_dp)
        write (*, '(a,3es13.5)') 'physical Piola curl FD ladder=', errors(1, :)
        write (*, '(a,3es13.5)') 'physical Piola div FD ladder=', errors(2, :)
    end subroutine

    subroutine physical_piola_value(field, map, x, a, b)
        type(free_b_field_t), intent(in) :: field
        type(nn_map_t), intent(in) :: map
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: a(3), b(3)
        real(dp) :: u(3), q(3), refx(2), w, val(field%fs%nl0)
        real(dp) :: grad(2, field%fs%nl0), pos(3), p1(3, 3), p2(3, 3), jac
        integer :: st
        ! Explicit inverse of the triangular shear and vertical ripple.
        u(3) = x(3)/0.9_dp
        u(2) = (x(2) - 0.15_dp*u(3))/1.1_dp
        u(1) = x(1) - 0.2_dp*u(2) - 0.1_dp*u(3)
        q(1) = sqrt(sum(u(:2)**2)); q(3) = atan2(u(2), u(1))
        q(2) = u(3) - 0.08_dp*sin(3*q(3))
        call meridian_qgrad(field%fs, 1, 1, refx, w, val, grad, st)
        if (st /= 0) error stop 'physical Piola reference basis'
        call nn_map_point(map, field%fs%sp(0)%dof(:, 1), val, grad, q(:2), &
            q(3), pos, p1, p2, jac, st)
        if (st /= 0) error stop 'physical Piola map'
        call small('physical Cartesian inverse chart', maxval(abs(pos - x)), &
            5.0e-12_dp)
        a = matmul(p1, [q(2)/2, q(1)/2, -q(1)/2])
        b = q(1)*p2(:, 1)
    end subroutine

    subroutine cartesian_identity()
        type(free_b_field_t), target :: legacy, cart
        type(nn_map_t), target :: oldmap, newmap
        type(cylindrical_chart_t) :: chart
        type(nn_mapped_report_t) :: report
        real(dp), allocatable :: olddata(:), newdata(:)
        real(dp) :: eo, en, deviation, ho, hn, exacto, exactn, scale
        integer :: st
        call make_field(2, 3, legacy); call make_field(2, 3, cart)
        call nn_map_build(legacy%fs, 0, 1, identity_map, oldmap, st)
        if (st /= 0) error stop 'Cartesian identity reference map'
        call nn_map_set_chart(chart, newmap)
        call nn_mapped_data(legacy, oldmap, 16, potential, olddata, st)
        if (st /= 0) error stop 'legacy identity data'
        call nn_mapped_cartesian_data(cart, newmap, 16, cartesian_potential, &
            newdata, st)
        if (st /= 0) error stop 'Cartesian identity data'
        scale = max(1.0_dp, maxval(abs(olddata)))
        call small('Cartesian/cylindrical data', &
            maxval(abs(newdata - olddata))/scale, 1.0e-10_dp)
        call nn_mapped_solve(legacy, oldmap, 16, olddata, 1.0e-12_dp, report, st)
        if (st /= 0) error stop 'legacy identity solve'
        call nn_mapped_solve(cart, newmap, 16, newdata, 1.0e-12_dp, report, st)
        if (st /= 0) error stop 'Cartesian identity solve'
        call small('Cartesian/cylindrical field coefficients', &
            maxval(abs(cart%b - legacy%b))/maxval(abs(legacy%b)), 1.0e-10_dp)
        call nn_mapped_error(legacy, oldmap, 33, exact_field, eo, deviation, ho, &
            exacto, st)
        if (st /= 0) error stop 'legacy identity norm'
        call nn_mapped_cartesian_error(cart, newmap, 33, cartesian_field, en, hn, &
            exactn, st)
        if (st /= 0) error stop 'Cartesian identity norm'
        call small('Cartesian/cylindrical L2 error', abs(en - eo), 1.0e-10_dp)
        call small('Cartesian/cylindrical energy', abs(hn - ho)/ho, 1.0e-10_dp)
        call small('Cartesian/cylindrical exact energy', &
            abs(exactn - exacto)/exacto, 1.0e-10_dp)
        call nn_mapped_data(cart, newmap, 16, potential, olddata, st)
        if (st == 0) error stop 'legacy data accepts Cartesian chart'
        call nn_mapped_error(cart, newmap, 33, exact_field, eo, deviation, ho, &
            exacto, st)
        if (st == 0) error stop 'legacy norm accepts Cartesian chart'
        write (*, '(a,3es13.5)') 'Cartesian identity L2/energy/exact=', en, hn, exactn
    end subroutine

    subroutine cartesian_quadrature()
        ! Shared spatial rule, independent Cartesian geometry/field and phi rule.
        integer, parameter :: np(3) = [12, 24, 48], reference_planes = 193
        type(free_b_field_t) :: field
        type(nn_map_t) :: map
        type(analytic_chart_t) :: chart
        real(dp) :: ref(3), ref_coarse(3), got(3), errs(3), err, eh, ee
        integer :: l, st
        call make_field(1, 1, field)
        field%b = 0
        call nn_map_set_chart(chart, map)
        call independent_norm(field, 97, ref_coarse, cartesian_field)
        call independent_norm(field, reference_planes, ref, cartesian_field)
        call small('independent phi reference refinement', &
            maxval(abs(ref - ref_coarse)/max(abs(ref), 1.0e-12_dp)), 1.0e-11_dp)
        do l = 1, 3
            call nn_mapped_cartesian_error(field, map, np(l), cartesian_field, &
                err, eh, ee, st)
            if (st /= 0) error stop 'Cartesian quadrature error status'
            got = [err, eh, ee]
            errs(l) = maxval(abs(got - ref)/max(abs(ref), 1.0e-12_dp))
        end do
        call small('independent phi refinement finest', errs(3), 1.0e-8_dp)
        call small('independent phi refinement reduction', errs(3), &
            max(2.0e-12_dp, 0.25_dp*errs(1)))
        write (*, '(a,3es13.5)') 'Cartesian independent phi12/24/48 vs193=', errs
    end subroutine

    subroutine cartesian_vacuum()
        ! Same physical chart/field on quadratic geometry, with p enrichment.
        ! R*A_phi contains R**2: after the curved P2 cell map it needs P4.
        ! The poloidal covariant pullback needs vector P3, contained in N4.
        integer, parameter :: orders(3) = [2, 3, 4], planes = 24
        type(free_b_field_t), target :: field
        type(nn_map_t), target :: map
        type(analytic_chart_t) :: chart
        type(nn_mapped_report_t) :: report
        real(dp), allocatable :: data(:)
        real(dp) :: error, eh, ee, ref(3), norm, discrepancy, errors(3)
        integer :: st, l, points
        do l = 1, size(orders)
            call make_field(2, 4, field, orders(l))
            field%f_vac = 0
            call nn_map_set_chart(chart, map)
            call nn_mapped_cartesian_data(field, map, planes, constant_potential, data, st)
            if (st /= 0) error stop 'full-chart constant vacuum data'
            call nn_mapped_solve(field, map, planes, data, 1.0e-12_dp, report, st)
            if (st /= 0) error stop 'full-chart constant vacuum solve'
            call nn_mapped_cartesian_error(field, map, 61, constant_field, error, eh, &
                ee, st)
            if (st /= 0) error stop 'full-chart constant vacuum norm'
            call independent_norm(field, 97, ref, constant_field)
            norm = sqrt(2*field%mu0*ref(3))
            errors(l) = ref(1)/norm
            discrepancy = max(abs(error - ref(1))/norm, abs(eh - ref(2))/ref(3), &
                abs(ee - ref(3))/ref(3))
            points = field%fs%mm%topology%n_triangles*field%fs%nq*planes
            ! Retain all metrics before assertions, including approximation error.
            write (*, '(a,i0,a,4es18.9,a,i0)') 'full-chart vacuum k=', orders(l), &
                ' L2/res/div/norm delta=', errors(l), report%residual, &
                report%divergence, discrepancy, ' iterations=', report%iterations
            write (*, '(a,6i10)') 'full-chart k/geometry/ndof/nq/nphi/points=', &
                orders(l), field%fs%g, report%ndof, field%fs%nq, planes, points
            call small('full-chart vacuum independent norm', discrepancy, 1.0e-10_dp)
            call small('full-chart vacuum residual', report%residual, 1.0e-9_dp)
            call small('full-chart vacuum divergence', report%divergence, 1.0e-12_dp)
        end do
        ! P2/P3 are explicit under-resolution controls, not exactness claims.
        ! Same original physical error gate, at the sufficient polynomial order.
        call small('full-chart vacuum relative L2', errors(3), 1.0e-7_dp)
    end subroutine

    subroutine independent_norm(field, nphi, result, exact)
        ! Toroidal MEAN: w already carries R dR dZ. No 2pi or nfp multiplier.
        type(free_b_field_t), intent(in) :: field
        integer, intent(in) :: nphi
        real(dp), intent(out) :: result(3)
        procedure(nn_vector_field_i) :: exact
        real(dp) :: x(2), w, basis(3, field%fs%sp(2)%nloc), point(3), d(3, 3)
        real(dp) :: jac, angle, b(3), bh(3), bhat(3), wt, volume, want_volume
        real(dp) :: carry_v, carry_w
        complex(dp) :: bn(3)
        integer :: t, q, j, st, n, i
        result = 0; volume = 0; want_volume = 0; carry_v = 0; carry_w = 0
        do t = 1, field%fs%mm%topology%n_triangles
            do q = 1, field%fs%nq
                call meridian_qpoint(field%fs, 2, t, q, x, w, basis, st)
                if (st /= 0) error stop 'independent norm spatial quadrature'
                do i = 1, size(basis, 2)
                    basis(:, i) = basis(:, i)*field%fs%sp(2)%sgn(i, t)
                end do
                ! Compensated sums: the measure check compares two sums of
                ! ~1e6 terms at 2e-13, below plain accumulation roundoff.
                call kahan_add(want_volume, carry_w, 0.99_dp*w)
                do j = 1, nphi
                    ! Distinct phase from the production midpoint rule.
                    angle = 2*pi*(j - 0.271_dp)/nphi
                    call expected_point(2, [x, angle], point, d, jac)
                    call exact(point, b)
                    bhat = [0.0_dp, field%f_vac/x(1), 0.0_dp]
                    do n = 0, ubound(field%b, 2)
                        bn = matmul(basis, field%b(field%fs%sp(2)%dof(:, t), n))
                        bhat = bhat + merge(1, 2, n == 0) &
                            *real(bn*cmplx(cos(n*angle), sin(n*angle), dp), dp)
                    end do
                    bh = matmul(d, bhat)/jac
                    wt = w*jac/nphi
                    result(1) = result(1) + wt*sum((bh - b)**2)
                    result(2) = result(2) + wt*sum(bh**2)/(2*field%mu0)
                    result(3) = result(3) + wt*sum(b**2)/(2*field%mu0)
                    call kahan_add(volume, carry_v, wt)
                end do
            end do
        end do
        result(1) = sqrt(result(1))
        call small('independent physical mean measure', &
            abs(volume - want_volume)/want_volume, 2.0e-13_dp)
    end subroutine

    subroutine kahan_add(sum, carry, term)
        real(dp), intent(inout) :: sum, carry
        real(dp), intent(in) :: term
        real(dp) :: y, t
        y = term - carry
        t = sum + y
        carry = (t - sum) - y
        sum = t
    end subroutine

    subroutine make_field(nr, modes, field, approximation_order)
        integer, intent(in) :: nr, modes
        integer, intent(in), optional :: approximation_order
        type(free_b_field_t), intent(out) :: field
        type(flux_mesh_t) :: fm
        real(dp) :: rings(nr)
        integer :: j, st, k
        k = 2
        if (present(approximation_order)) k = approximation_order
        do j = 1, nr
            rings(j) = real(j, dp)/nr
        end do
        call flux_mesh_build(rings, 8, 4096, 2, circle, fm, st)
        if (st /= 0) error stop 'mesh'
        call nn_build(fm, k, modes, [r0, 0.0_dp], minor, ftor, 1.0_dp, field, st)
        if (st /= 0) error stop 'build'
    end subroutine

    subroutine identity_regression()
        ! Identity phi-dependence: the coupled FGMRES solve must reproduce
        ! the per-n NN1 solves for the same traces to rounding.
        integer, parameter :: nr = 2, modes = 3, nphi = 4*modes + 4
        type(free_b_field_t), target :: mapped, ref
        type(nn_map_t), target :: map
        type(nn_mapped_report_t) :: report
        type(nn_report_t) :: r1
        real(dp), allocatable :: data(:)
        complex(dp), allocatable :: bc(:), seed(:)
        real(dp) :: diff, err, dev, eh, ee, err_ref
        integer :: st, n, n1, nn, n0
        call make_field(nr, modes, mapped)
        call make_field(nr, modes, ref)
        call nn_map_build(mapped%fs, 0, 1, identity_map, map, st)
        if (st /= 0) error stop 'identity map'
        call nn_mapped_data(mapped, map, nphi, potential, data, st)
        if (st /= 0) error stop 'identity data'
        call nn_mapped_solve(mapped, map, nphi, data, 1.0e-12_dp, report, st)
        if (st /= 0) error stop 'identity mapped solve'
        n1 = mapped%fs%sp(1)%ndof; nn = mapped%fs%nn; n0 = mapped%fs%n0
        bc = cmplx(data(nn + 1:n1), 0.0_dp, dp)
        seed = 0*bc
        call nn_solve(ref, 0, bc, seed, r1, st)
        if (st /= 0) error stop 'NN1 n=0'
        do n = 1, modes
            bc = cmplx(data((2*n - 1)*n1 + 1:(2*n - 1)*n1 + nn), &
                -data(2*n*n1 + 1:2*n*n1 + nn), dp)/2
            seed = 0*bc
            call nn_solve(ref, n, bc, seed, r1, st)
            if (st /= 0) error stop 'NN1 mode'
        end do
        diff = maxval(abs(mapped%b - ref%b))/maxval(abs(ref%b))
        call nn_mapped_error(mapped, map, 2*nphi + 1, exact_field, err, dev, eh, ee, st)
        if (st /= 0) error stop 'identity error'
        call nn_mapped_error(ref, map, 2*nphi + 1, exact_field, err_ref, dev, eh, ee, st)
        if (st /= 0) error stop 'identity error ref'
        write (*, '(a,es10.3,a,i0,a,es10.3,a,es10.3,a,es10.3)') &
            'identity: max|b-b_NN1|/max|b|=', diff, ' iterations=', &
            report%iterations, ' divergence=', report%divergence, &
            ' field error=', err/dev, ' NN1=', err_ref/dev
        if (diff > 1.0e-10_dp) error stop 'identity map does not reproduce NN1'
        if (report%iterations > 3) error stop 'identity preconditioner not exact'
        if (report%divergence > 1.0e-12_dp) error stop 'identity divergence'
        if (abs(err - err_ref) > 1.0e-9_dp*dev) error stop 'identity error mismatch'
    end subroutine

    subroutine ellipse_case(nr, modes, rel)
        integer, intent(in) :: nr, modes
        real(dp), intent(out) :: rel
        type(free_b_field_t), target :: field
        type(nn_map_t), target :: map
        type(nn_mapped_report_t) :: report
        real(dp), allocatable :: data(:)
        real(dp) :: err, dev, eh, ee
        integer :: st, nphi, c0, c1, rate
        nphi = 4*modes + 4
        call system_clock(c0, rate)
        call make_field(nr, modes, field)
        call nn_map_build(field%fs, 2, 8, ellipse_map, map, st)
        if (st /= 0) error stop 'ellipse map'
        call nn_mapped_data(field, map, nphi, potential, data, st)
        if (st /= 0) error stop 'ellipse data'
        call nn_mapped_solve(field, map, nphi, data, 1.0e-14_dp, report, st)
        if (st /= 0) error stop 'ellipse solve'
        call nn_mapped_error(field, map, 2*nphi + 1, exact_field, err, dev, eh, ee, st)
        if (st /= 0) error stop 'ellipse error'
        call system_clock(c1)
        rel = err/dev
        write (*, '(a,i0,a,i0,a,i0,a,es10.3,a,es10.3,a,i0,a,es9.2,a,es9.2,a,f8.2)') &
            'ellipse nr=', nr, ' N=', modes, ' ndof=', report%ndof, &
            ' field error=', rel, ' energy rel=', abs(report%energy - ee)/ee, &
            ' fgmres=', report%iterations, ' residual=', report%residual, &
            ' div=', report%divergence, ' s=', real(c1 - c0, dp)/rate
        if (report%divergence > 1.0e-12_dp) error stop 'mapped divergence identity'
        if (report%residual > 1.0e-9_dp) error stop 'mapped residual'
    end subroutine

    subroutine zero_forcing_residual()
        ! f_vac = 0 leaves only the trace data as forcing (rhs = 0). The
        ! reported residual must follow the FGMRES tolerance there too; a
        ! normalisation by rhs and S x alone is identically 1 in this case.
        integer, parameter :: nr = 2, modes = 1, nphi = 4*modes + 4
        real(dp), parameter :: loose = 1.0e-4_dp, tight = 1.0e-12_dp
        type(free_b_field_t), target :: field
        type(nn_map_t), target :: map
        type(nn_mapped_report_t) :: report
        real(dp), allocatable :: data(:)
        real(dp) :: res_loose, res_tight
        integer :: st
        call make_field(nr, modes, field)
        field%f_vac = 0
        call nn_map_build(field%fs, 2, 8, ellipse_map, map, st)
        if (st /= 0) error stop 'zero-forcing map'
        call nn_mapped_data(field, map, nphi, potential, data, st)
        if (st /= 0) error stop 'zero-forcing data'
        if (maxval(abs(data)) <= 0) error stop 'zero-forcing data vanish'
        call nn_mapped_solve(field, map, nphi, data, loose, report, st)
        if (st /= 0) error stop 'zero-forcing loose solve'
        res_loose = report%residual
        call nn_mapped_solve(field, map, nphi, data, tight, report, st)
        if (st /= 0) error stop 'zero-forcing tight solve'
        res_tight = report%residual
        write (*, '(a,2es12.4,a,i0)') 'zero forcing residual loose/tight=', &
            res_loose, res_tight, ' iterations=', report%iterations
        if (res_loose > 2*loose) error stop 'zero-forcing loose residual'
        if (res_tight > 1.0e-9_dp) error stop 'zero-forcing residual'
        if (res_tight >= res_loose) error stop 'zero-forcing residual not reduced'
    end subroutine
end program test_nonnested_mapped
