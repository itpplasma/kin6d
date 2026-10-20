program test_curve_chart
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan, ieee_value, ieee_quiet_nan
    use kin6d_equilibrium_chart, only: curve_chart_t, cylindrical_chart_t
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp)
    type(curve_chart_t) :: chart
    real(dp) :: c(3, 0:4), s(3, 0:4)
    character(40) :: mode
    integer :: status
    call get_command_argument(1, mode)
    select case (trim(mode))
    case ('circle')
        call circle()
    case ('jets')
        call jets()
    case ('tube')
        call tube()
    case ('periodic')
        call periodic()
    case ('invalid')
        call invalid()
    case default
        error stop 'test/test_curve_chart.f90: unknown test mode'
    end select
contains
    subroutine require(ok, label, line)
        logical, intent(in) :: ok
        character(*), intent(in) :: label
        integer, intent(in) :: line
        if (ok) return
        write (error_unit, '(a,i0,2a)') 'test/test_curve_chart.f90:', line, ': ', label
        error stop 1
    end subroutine

    subroutine fixture(family, nfp)
        integer, intent(in) :: family
        integer, intent(out) :: nfp
        real(dp) :: gc(3, 0:1), gs(3, 0:1)
        c = 0; s = 0; nfp = 1
        select case (family)
        case (1) ! Circle.
            c(1, 1) = 3; s(2, 1) = 3
        case (2, 5) ! Closed toroidal helix: (3+.2 cos 3t)e_R+.2 sin 3t e_z.
            nfp = 3
            c(1, 1) = 3; s(2, 1) = 3
            c(1, 2) = .1_dp; s(2, 2) = -.1_dp
            c(1, 4) = .1_dp; s(2, 4) = .1_dp; s(3, 3) = .2_dp
        case (3) ! Analytic lifted Gerono figure-8, not a published equilibrium.
            s(1, 1) = 1; s(2, 2) = 1; c(3, 1) = .3_dp
        case (4) ! Planar Gerono: exact inflections but crossing axis (zero reach).
            s(1, 1) = 1; s(2, 2) = 1
        end select
        if (family == 5) then
            gc = 0; gs = 0; gc(3, 0) = 1; gc(1, 1) = .1_dp; gs(2, 1) = .1_dp
            call chart%init(c, s, nfp, status, gc, gs)
        else
            call chart%init(c, s, nfp, status)
        end if
        call require(status == 0, 'initialize regular projected-guide chart', __LINE__)
    end subroutine

    subroutine circle()
        type(cylindrical_chart_t) :: cylinder
        real(dp) :: q(3), x(3), d1(3, 3), d2(3, 3, 3), d3(3, 3, 3, 3)
        real(dp) :: y(3), e1(3, 3), e2(3, 3, 3), e3(3, 3, 3, 3), err
        integer, parameter :: perm(3) = [2, 1, 3]
        integer :: a, b, d, i, nfp
        call fixture(1, nfp)
        err = 0
        do i = 0, 31
            q = [.13_dp, -.27_dp, 2*pi*i/32]
            call chart%jet(q, x, d1, d2, d3)
            call cylinder%jet([3 + q(2), q(1), q(3)], y, e1, e2, e3)
            err = max(err, maxval(abs(x - y)))
            do a = 1, 3
                err = max(err, maxval(abs(d1(:, a) - e1(:, perm(a)))))
                do b = 1, 3
                    err = max(err, maxval(abs(d2(:, a, b) - e2(:, perm(a), perm(b)))))
                    do d = 1, 3
                        err = max(err, maxval(abs(d3(:, a, b, d) &
                                                  - e3(:, perm(a), perm(b), perm(d)))))
                    end do
                end do
            end do
            call require(abs(det(d1) - (3 + q(2))) < 2e-13_dp, &
                         'positive circular determinant R0+q2', __LINE__)
        end do
        call chart%jet([0.0_dp, -3.0_dp, 0.0_dp], x, d1, d2, d3)
        call require(abs(det(d1)) < 1e-13_dp, 'curvature limit is singular', __LINE__)
        call chart%jet([0.0_dp, -3.1_dp, 0.0_dp], x, d1, d2, d3)
        call require(det(d1) < 0, &
            'beyond the circular tube limit reverses orientation', __LINE__)
        write (*, '(a,es12.4)') 'VERIFIED circle H through d3 max error: ', err
        call require(err < 2e-13_dp, &
            'fixed relabelling of complete circular jet', __LINE__)
    end subroutine

    function position(q, family) result(x)
        ! Independent direct Cartesian formula, no production Fourier/frame leaves.
        real(dp), intent(in) :: q(3)
        integer, intent(in) :: family
        real(dp) :: x(3), axis(3), v(3), t(3), n(3), b(3), a, guide(3)
        a = q(3)
        select case (family)
        case (2, 5)
            axis = [(3 + .2_dp*cos(3*a))*cos(a), &
                    (3 + .2_dp*cos(3*a))*sin(a), .2_dp*sin(3*a)]
            v = [-.6_dp*sin(3*a)*cos(a) - (3 + .2_dp*cos(3*a))*sin(a), &
                 -.6_dp*sin(3*a)*sin(a) + (3 + .2_dp*cos(3*a))*cos(a), &
                 .6_dp*cos(3*a)]
        case (3, 4)
            axis = [sin(a), sin(2*a), .3_dp*cos(a)]
            v = [cos(a), 2*cos(2*a), -.3_dp*sin(a)]
            if (family == 4) then
                axis(3) = 0; v(3) = 0
            end if
        end select
        t = v/sqrt(sum(v*v))
        guide = [0.0_dp, 0.0_dp, 1.0_dp]
        if (family == 5) guide = [.1_dp*cos(a), .1_dp*sin(a), 1.0_dp]
        n = guide - dot_product(guide, t)*t
        n = n/sqrt(sum(n*n))
        b = cross(t, n)
        x = axis + q(1)*n + q(2)*b
    end function

    subroutine weights(order, w)
        ! Independent moment equations for a centered 9-point FD stencil.
        integer, intent(in) :: order
        real(dp), intent(out) :: w(-4:4)
        real(dp) :: mat(0:8, 0:9), row(0:9), pivot
        integer :: i, j, k, p
        do i = 0, 8
            do j = 0, 8
                mat(i, j) = real(j - 4, dp)**i
            end do
        end do
        mat(:, 9) = 0
        mat(order, 9) = 1
        do i = 2, order
            mat(order, 9) = mat(order, 9)*i
        end do
        do i = 0, 8
            p = i + maxloc(abs(mat(i:8, i)), dim=1) - 1
            row = mat(i, :); mat(i, :) = mat(p, :); mat(p, :) = row
            pivot = mat(i, i); mat(i, :) = mat(i, :)/pivot
            do k = 0, 8
                if (k /= i) mat(k, :) = mat(k, :) - mat(k, i)*mat(i, :)
            end do
        end do
        w = mat(:, 9)
    end subroutine

    function difference(q, family, orders, h) result(value)
        real(dp), intent(in) :: q(3), h
        integer, intent(in) :: family, orders(3)
        real(dp) :: value(3), w(-4:4, 3), p(3), term, steps(3)
        integer :: i, j, k, lo(3), hi(3), d
        w = 0; w(0, :) = 1; lo = 0; hi = 0
        steps = [.1_dp, .1_dp, h]
        do d = 1, 3
            if (orders(d) == 0) cycle
            call weights(orders(d), w(:, d))
            lo(d) = -4; hi(d) = 4
        end do
        value = 0
        do k = lo(3), hi(3)
            do j = lo(2), hi(2)
                do i = lo(1), hi(1)
                    p = q + steps*[i, j, k]
                    term = w(i, 1)*w(j, 2)*w(k, 3)
                    value = value + term*position(p, family)
                end do
            end do
        end do
        value = value/product(steps**orders)
    end function

    subroutine jets()
        real(dp) :: q(3), x(3), d1(3, 3), d2(3, 3, 3), d3(3, 3, 3, 3), fd(3)
        real(dp) :: err(0:3), h, scale
        integer :: family, nfp, i, j, k, l, orders(3), pass
        do family = 2, 5
            call fixture(family, nfp)
            do pass = 1, 2
                err = 0; h = .008_dp/pass
                do i = 0, 10
                    q = [.017_dp, -.023_dp, 2*pi*i/11]
                    call chart%jet(q, x, d1, d2, d3)
                    err(0) = max(err(0), maxval(abs(x - position(q, family))))
                    do j = 1, 3
                        orders = 0; orders(j) = 1
                        fd = difference(q, family, orders, h)
                        scale = max(1.0_dp, maxval(abs(d1(:, j))))
                        err(1) = max(err(1), maxval(abs(fd - d1(:, j)))/scale)
                        do k = 1, 3
                            orders = 0; orders(j) = orders(j) + 1
                            orders(k) = orders(k) + 1
                            fd = difference(q, family, orders, h)
                            scale = max(1.0_dp, maxval(abs(d2(:, j, k))))
                            err(2) = max(err(2), maxval(abs(fd - d2(:, j, k)))/scale)
                            do l = 1, 3
                                orders = 0; orders(j) = orders(j) + 1
                                orders(k) = orders(k) + 1; orders(l) = orders(l) + 1
                                fd = difference(q, family, orders, h)
                                scale = max(1.0_dp, maxval(abs(d3(:, j, k, l))))
                                err(3) = max(err(3), &
                                    maxval(abs(fd - d3(:, j, k, l)))/scale)
                            end do
                        end do
                    end do
                end do
                write (*, '(a,2i3,4es12.3)') 'VERIFIED FD family/step, errors: ', &
                    family, pass, err
                call require(err(0) < 2e-14_dp, 'independent chart position', __LINE__)
                call require(err(1) < 5e-9_dp .and. err(2) < 2e-7_dp &
                             .and. err(3) < 2e-5_dp, &
                             '9-point differences of H for every jet entry', __LINE__)
            end do
        end do
    end subroutine

    subroutine tube()
        real(dp) :: rc, rt, r, q(3), x(3), d1(3, 3), d2(3, 3, 3), d3(3, 3, 3, 3)
        real(dp) :: minimum, identity, tangent(3), k(3), g(3, 3)
        integer :: family, nfp, i, j
        do family = 1, 4
            call fixture(family, nfp)
            call chart%tube_radii(2048, rc, rt, status)
            call require(status == 0, 'tube diagnostic resolved', __LINE__)
            write (*, '(a,i2,2es14.5)') &
                'VERIFIED numerical curvature ceiling / nonoverlap radius family ', &
                family, rc, rt
            if (family == 4) then
                call require(abs(rt) < tiny(rt), &
                    'self-crossing axis has no admitted tube', __LINE__)
                cycle
            end if
            call require(rt > 0, &
                'embedded axis has positive nonoverlap radius', __LINE__)
            if (family == 1) call require(abs(rc - 3) < 1e-12_dp, &
                'circle curvature radius equals major radius', __LINE__)
            r = .8_dp*rt; minimum = huge(r)
            do i = 0, 127
                do j = 0, 15
                    q = [r*cos(2*pi*j/16), r*sin(2*pi*j/16), 2*pi*i/128]
                    call chart%jet([0.0_dp, 0.0_dp, q(3)], x, d1, d2, d3)
                    tangent = d1(:, 3)/sqrt(sum(d1(:, 3)**2))
                    k = (d2(:, 3, 3) - tangent*dot_product(tangent, d2(:, 3, 3))) &
                        /sum(d1(:, 3)**2)
                    identity = sqrt(sum(d1(:, 3)**2)) &
                               *(1 - q(1)*dot_product(k, d1(:, 1)) &
                                 - q(2)*dot_product(k, d1(:, 2)))
                    call chart%jet(q, x, d1, d2, d3)
                    minimum = min(minimum, det(d1))
                    call require(abs(det(d1) - identity) < 1e-12_dp, &
                                 'tube determinant curvature identity', __LINE__)
                    g = matmul(transpose(d1), d1)
                    call require(abs(det(g) - det(d1)**2) < 2e-11_dp, &
                                 'metric determinant equals squared Jacobian', __LINE__)
                end do
            end do
            call require(minimum > 0, 'positive tube orientation and metric', __LINE__)
            write (*, '(a,es14.5)') 'VERIFIED sampled minimum det DH: ', minimum
        end do
    end subroutine

    subroutine periodic()
        real(dp) :: q(3), x(3), d1(3, 3), d2(3, 3, 3), d3(3, 3, 3, 3)
        real(dp) :: y(3), e1(3, 3), e2(3, 3, 3), e3(3, 3, 3, 3), rot(3, 3), angle
        real(dp) :: gc(3, 0:1), gs(3, 0:1), err
        integer :: family, nfp, i, j, k, l
        err = 0
        do family = 1, 3
            call fixture(family, nfp)
            ! A nonconstant guide exercises its derivatives and equivariance.
            gc = 0; gs = 0; gc(3, 0) = 1
            gc(1, 1) = .1_dp; gs(2, 1) = .1_dp
            call chart%init(c, s, nfp, status, gc, gs)
            call require(status == 0, 'admit periodic nonconstant guide', __LINE__)
            do l = 1, 2
                angle = 2*pi
                if (l == 2) angle = 2*pi/nfp
                rot = 0
                rot(1, :) = [cos(angle), -sin(angle), 0.0_dp]
                rot(2, :) = [sin(angle), cos(angle), 0.0_dp]; rot(3, 3) = 1
                do i = 0, 16
                    q = [.03_dp, -.04_dp, 2*pi*i/17]
                    call chart%jet(q, x, d1, d2, d3)
                    q(3) = q(3) + angle
                    call chart%jet(q, y, e1, e2, e3)
                    err = max(err, maxval(abs(y - matmul(rot, x))))
                    do j = 1, 3
                        err = max(err, maxval(abs(e1(:, j) - matmul(rot, d1(:, j)))))
                        do k = 1, 3
                            err = max(err, maxval(abs(e2(:, j, k) &
                                                      - matmul(rot, d2(:, j, k)))))
                        end do
                    end do
                    do j = 1, 3
                        do k = 1, 3
                            err = max(err, maxval(abs(e3(:, j, k, 1) &
                                                      - matmul(rot, d3(:, j, k, 1)))))
                            err = max(err, maxval(abs(e3(:, j, k, 2) &
                                                      - matmul(rot, d3(:, j, k, 2)))))
                            err = max(err, maxval(abs(e3(:, j, k, 3) &
                                                      - matmul(rot, d3(:, j, k, 3)))))
                        end do
                    end do
                end do
            end do
        end do
        write (*, '(a,es12.4)') 'VERIFIED periodic/full nfp jet equivariance: ', err
        call require(err < 2e-10_dp, &
            'full jet closure and rotational nfp symmetry', __LINE__)
    end subroutine

    subroutine invalid()
        real(dp) :: gc(3, 0:0), gs(3, 0:0), x(3), a(3, 3)
        real(dp) :: b(3, 3, 3), d(3, 3, 3, 3), r, t
        integer :: nfp
        call fixture(1, nfp)
        call chart%init(c, s, 3, status)
        call require(status == 0, 'circle has rotational symmetry for nfp=3', __LINE__)
        c(1, 2) = .1_dp
        call chart%init(c, s, 3, status)
        call require(status == 1, 'reject incompatible nfp coefficients', __LINE__)
        call chart%jet([0.0_dp, 0.0_dp, 0.0_dp], x, a, b, d)
        call require(all(ieee_is_nan(x)), &
            'failed init invalidates previous chart', __LINE__)
        c = 0; s = 0
        call chart%init(c, s, 1, status)
        call require(status == 2, 'reject stationary axis', __LINE__)
        call chart%tube_radii(256, r, t, status)
        call require(status /= 0, 'no radius on invalid chart', __LINE__)
        call fixture(1, nfp)
        gc = 0; gs = 0; gc(1, 0) = 1
        call chart%init(c, s, 1, status, gc, gs)
        call require(status == 2, 'reject guide tangent to axis', __LINE__)
        call chart%init(c(:, 0:-1), s(:, 0:-1), 1, status)
        call require(status == 1, 'reject empty Fourier series', __LINE__)
        call fixture(1, nfp)
        call chart%init(c, s, 0, status)
        call require(status == 1, 'reject zero field periods', __LINE__)
        call chart%init(c, s, 1, status, guide_c=gc)
        call require(status == 1, 'reject an unpaired guide array', __LINE__)
        call fixture(1, nfp)
        call chart%tube_radii(16, r, t, status)
        call require(status == 1, 'reject unresolved tube sampling', __LINE__)
        call fixture(1, nfp)
        c(1, 0) = ieee_value(0.0_dp, ieee_quiet_nan)
        call chart%init(c, s, 1, status)
        call require(status == 1, 'reject nonfinite coefficients', __LINE__)
    end subroutine

    pure function cross(a, b) result(c)
        real(dp), intent(in) :: a(3), b(3)
        real(dp) :: c(3)
        c = [a(2)*b(3) - a(3)*b(2), a(3)*b(1) - a(1)*b(3), a(1)*b(2) - a(2)*b(1)]
    end function
    pure function det(a) result(d)
        real(dp), intent(in) :: a(3, 3)
        real(dp) :: d
        d = dot_product(a(:, 1), cross(a(:, 2), a(:, 3)))
    end function
end program
