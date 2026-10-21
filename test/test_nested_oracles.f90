program test_nested_oracles
    ! Independent complex formulas, spatial differences and volume integration.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_nested_oracles, only: nested_oracle_t, NESTED_INTEGER, NESTED_SHEARED
    use kin6d_nested_oracle_leaves, only: nested_integer_averages, &
                                   nested_sheared_field_jvp, nested_sheared_field_vjp, &
                         nested_sheared_field_interval, nested_integer_field_interval, &
                                      nested_integer_field_jvp, nested_integer_field_vjp
    use fortnum_interval, only: interval_t, interval
    use fortnum_quadrature, only: gauss_legendre_ab
    implicit none
    type(nested_oracle_t) :: o
    real(dp), parameter :: pi = acos(-1.0_dp), step = 1e-4_dp
    real(dp) :: x(3), xx(3), v(20), ref(5), vr(5), vp(5), vm(5), vpp(5), vmm(5)
real(dp) :: db(3, 3), gp(3), curl(3), force(3), jac, rho, theta, t, error, divmax, fbmax
   real(dp) :: moments(4), fine(4), vol, pavg, b2, beta, iota, flux, lo, hi, limit, axis
    real(dp) :: pa
    real(dp) :: a(8), direction(8), weights(20), tangent(20), adjoint(8), fd(20)
    real(dp) :: plus(20), minus(20), jmat(3, 3), xp(3), xm(3), unused
    type(interval_t) :: ai(8), vi(20)
    integer :: family, asym, k, n, status
    error = 0; divmax = 0; fbmax = 0
    do family = 1, 2
        do asym = 0, 1
            o%family = family
            o%eps = 0.5_dp
            o%tau = asym*0.4_dp
            if (family == 2) o%eps = 1.2_dp
            o%delta = 1.0_dp/64
            if (family == 2) o%delta = 0.14_dp
            if (.not. o%valid()) error stop 'test_nested_oracles: valid case rejected'
            do n = 1, 24
                theta = 2*pi*(n - 0.3_dp)/24; t = 2*pi*(7*n - 0.2_dp)/31
                rho = sqrt(o%delta)*(0.2_dp + 0.8_dp*n/24)
                if (family == 2) rho = sqrt(2.0_dp)*rho
                call o%surface(rho, theta, t, x, jac)
                call o%point(x, v, status)
                if (status /= 0) error stop 'test_nested_oracles: point rejected'
                call reference(o, x, ref)
                error = max(error, maxval(abs(v(1:5) - ref)))
                limit = rho*rho
                if (family == 2) limit = limit/2
                if (abs(v(1) - limit) > 2e-12_dp) error stop 'surface label identity'
                do k = 1, 3
                    xx = x; xx(k) = x(k) + step; call reference(o, xx, vp)
                    xx = x; xx(k) = x(k) - step; call reference(o, xx, vm)
                    xx = x; xx(k) = x(k) + 2*step; call reference(o, xx, vpp)
                    xx = x; xx(k) = x(k) - 2*step; call reference(o, xx, vmm)
                    db(:, k) = (-vpp(3:5) + 8*vp(3:5) - 8*vm(3:5) + vmm(3:5))/(12*step)
                    gp(k) = (-vpp(2) + 8*vp(2) - 8*vm(2) + vmm(2))/(12*step)
                end do
                curl = [db(3, 2) - db(2, 3), db(1, 3) - db(3, 1), db(2, 1) - db(1, 2)]
                if (maxval(abs(curl - v(6:8))) > 2e-7_dp) error stop 'independent curl'
                if (maxval(abs(gp - v(9:11))) > 2e-7_dp) error stop 'independent grad p'
                if (maxval(abs(reshape(db, [9]) - v(12:20))) > 2e-7_dp) &
                    error stop 'independent grad B'
                force = cross(curl, ref(3:5)) - gp
                divmax = max(divmax, abs(db(1, 1) + db(2, 2) + db(3, 3)))
                fbmax = max(fbmax, maxval(abs(force)))
                ! Independent Jacobian of the retained physical surface chart.
                do k = 1, 3
                    xx = [rho, theta, t]; xx(k) = xx(k) + step
                    call o%surface(xx(1), xx(2), xx(3), xp, unused)
                    xx = [rho, theta, t]; xx(k) = xx(k) - step
                    call o%surface(xx(1), xx(2), xx(3), xm, unused)
                    jmat(:, k) = (xp - xm)/(2*step)
                end do
                if (abs(abs(dot_product(jmat(:, 1), cross(jmat(:, 2), jmat(:, 3)))) &
                        - jac) > 2e-7_dp) error stop 'surface volume Jacobian'
            end do
            call o%boundary(0.8_dp, 1.2_dp, x, status)
            if (status /= 0) error stop 'geometric phi boundary root'
            if (abs(atan2(x(2), x(1)) - 1.2_dp) > 2e-12_dp) error stop 'geometric phi'
            call o%point(x, v, status)
            if (abs(v(1) - o%delta) > 2e-12_dp) error stop 'export boundary label'
            xx = [x(1), -x(2), -x(3)]
            call reference(o, xx, vr)
     if (asym == 0 .and. abs(vr(2) - v(2)) > 2e-12_dp) error stop 'stellarator symmetry'
         if (asym == 1 .and. abs(vr(2) - v(2)) < 1e-5_dp) error stop 'missing asymmetry'
            call o%averages(12, moments)
            call o%averages(20, fine)
            if (maxval(abs(moments - fine)) > 2e-5_dp) error stop 'volume refinement'
            if (family == 1) then
                pa = 2*o%delta/(1 - o%tau**2)
                call nested_integer_averages(o%eps, o%tau, pa, o%delta, vol, &
                                             pavg, b2, beta, iota, flux)
                if (maxval(abs(fine - [vol, pavg, b2, beta])) > 2e-11_dp) &
                    error stop 'independent integer averages'
            else
                vol = independent_volume(o)
               if (abs(vol - fine(1)) > 2e-6_dp) error stop 'independent sheared volume'
            end if
            print '(a,2i2,4es15.6)', '[V] family,asym,V,<p>,<B^2>,beta ', &
                family, asym, fine
        end do
    end do
    if (error > 2e-12_dp .or. divmax > 2e-7_dp .or. fbmax > 2e-7_dp) &
        error stop 'independent Cartesian identities'
    print '(a,3es14.5)', '[V] max field,div,force differences ', error, divmax, fbmax
    ! Axis transform from the independent published numerical value.
    o%family = NESTED_SHEARED; o%eps = 2; o%saxis = 1; o%delta = 0.005_dp; o%tau = 0
    call o%rotation(0.0_dp, 8, 1e-7_dp, axis, status)
    if (status /= 0 .or. abs(axis - 2.28690_dp) > 1e-5_dp) error stop 'axis iota oracle'
    call o%rotation(o%delta, 64, 1e-7_dp, lo, status)
    call o%rotation(o%delta, 128, 3e-8_dp, hi, status)
    if (status /= 0 .or. abs(lo - hi) > 0.01_dp) error stop 'finite turn rotation'
    print '(a,3f12.7)', '[V] sheared iota: axis,64,128 turns ', axis, lo, hi
    ! Generated field parameter JVP/VJP and strict point enclosure.
    a = [1.1_dp, 0.3_dp, 0.05_dp, 1.2_dp, 0.4_dp, 1.4_dp, 1.6_dp, 0.1_dp]
    direction = [0.1_dp, -0.2_dp, 0.3_dp, 0.15_dp, 0.1_dp, -0.12_dp, 0.2_dp, 0.4_dp]
    weights = [(sin(real(k, dp)), k=1, 20)]
    call jvp_leaf(a, direction, tangent)
    call vjp_leaf(a, weights, adjoint)
  if (abs(dot_product(weights, tangent) - dot_product(adjoint, direction)) > 1e-11_dp) &
        error stop 'JVP VJP pairing'
    o%eps = a(4); o%tau = a(5); o%saxis = a(6); o%lam = a(7)
    ! FD checks the direct generated leaf, with pa independent of all parameters.
    call field_leaf(a + 1e-6_dp*direction, plus)
    call field_leaf(a - 1e-6_dp*direction, minus)
    fd = (plus - minus)/2e-6_dp
    if (maxval(abs(fd - tangent)) > 2e-7_dp) error stop 'JVP independent differences'
    ai = interval(a)
    call interval_leaf(ai, vi)
    call field_leaf(a, plus)
    do k = 1, 20
        if (plus(k) < vi(k)%lo .or. plus(k) > vi(k)%hi) error stop 'strict enclosure'
    end do
    a(1:6) = [1.0_dp, 0.3_dp, 0.05_dp, 0.5_dp, 0.4_dp, 0.1_dp]
    call jvp_leaf(a(1:6), direction(1:6), tangent)
    call vjp_leaf(a(1:6), weights, adjoint(1:6))
    if (abs(dot_product(weights, tangent) - &
  dot_product(adjoint(1:6), direction(1:6))) > 1e-11_dp) error stop 'integer AD pairing'
    call field_leaf(a(1:6) + 1e-6_dp*direction(1:6), plus)
    call field_leaf(a(1:6) - 1e-6_dp*direction(1:6), minus)
if (maxval(abs((plus - minus)/2e-6_dp - tangent)) > 2e-7_dp) error stop 'integer JVP FD'
    ai(1:6) = interval(a(1:6))
    call interval_leaf(ai(1:6), vi)
    call field_leaf(a(1:6), plus)
    do k = 1, 20
        if (plus(k) < vi(k)%lo .or. plus(k) > vi(k)%hi) error stop 'integer interval'
    end do
    print '(a)', '[V] JVP/VJP, finite differences and strict point containment PASS'
contains
    subroutine reference(o, x, v)
        type(nested_oracle_t), intent(in)::o
        real(dp), intent(in)::x(3)
        real(dp), intent(out)::v(5)
        complex(dp)::omega, kk, xi, ww, bxy
        real(dp)::a, b, s, f, om2, q, pa, mag(3)
        if (o%family == NESTED_INTEGER) then
            a = sqrt(1 + o%eps); b = sqrt(1 - o%eps)
            s = (x(1)/a)**2 + (x(2)/b)**2
            f = sqrt(1 - (1 - s)**2 - 4*x(3)**2 - 4*o%tau*s*x(3))
            mag = [((2*x(3) + o%tau*s)*x(1) - a*f*x(2)/b)/s, &
                   ((2*x(3) + o%tau*s)*x(2) + b*f*x(1)/a)/s, 1 - s - 2*o%tau*x(3)]
            om2 = 1 - o%tau**2
            q = om2*(x(1)**2 + x(2)**2 + 4*x(3)**2)/2 + 2*o%tau*x(3)
            v(1) = om2/2*(q + sum(mag**2)/2 - 1 + o%tau**2/(2*om2) + om2*o%eps**2/2)
            pa = o%edge_pressure + 2*o%delta/om2
            v(2) = pa - 2*v(1)/om2
        else
            omega = cmplx(x(1), x(2), dp)
            kk = conjg(omega)*sqrt(1 + o%eps/conjg(omega)**2)
            xi = omega*kk + cmplx(pi/2 - o%saxis, o%tau, dp)
            ww = cmplx(0.0_dp, 1.0_dp, dp)*sin(xi)/(2*kk)
            bxy = exp(cmplx(0.0_dp, -o%lam*x(3), dp))*ww
            mag = [real(bxy, dp), aimag(bxy), &
                   real(exp(cmplx(0.0_dp, -o%lam*x(3), dp))*cos(xi), dp)/o%lam]
            v(1) = (sin(o%lam*x(3))**2 + (o%lam*mag(3))**2)/2
            v(2) = o%edge_pressure + (o%delta - v(1))/o%lam**2
        end if
        v(3:5) = mag
    end subroutine
    function cross(a, b) result(c)
        real(dp), intent(in)::a(3), b(3)
        real(dp)::c(3)
        c = [a(2)*b(3) - a(3)*b(2), a(3)*b(1) - a(1)*b(3), a(1)*b(2) - a(2)*b(1)]
    end function
    function independent_volume(o) result(vol)
        type(nested_oracle_t), intent(in)::o
        real(dp)::vol, d, nu, xi, u, t
        real(dp), allocatable::q(:), w(:)
        integer::i, j
        d = asin(sqrt(2*o%delta)); vol = 0
        allocate (q(32), w(32))
        call gauss_legendre_ab(32, -pi/2, pi/2, q, w)
        do j = 0, 63
            t = 2*pi*j/64; nu = o%eps*sin(2*t)/2 + o%tau
            do i = 1, 32
                xi = d*sin(q(i))
                u = asin(sqrt(max(0.0_dp, 2*o%delta - sin(xi)**2) &
                              /(cosh(nu)**2 - sin(xi)**2)))
                vol = vol + 2*pi/64*w(i)*u*d*cos(q(i))/o%lam
            end do
        end do
    end function
    subroutine field_leaf(a, v)
        use kin6d_nested_oracle_leaves, only: nested_sheared_field, nested_integer_field
        real(dp), intent(in)::a(:)
        real(dp), intent(out)::v(20)
        if (size(a) == 6) then
            call nested_integer_field(a(1), a(2), a(3), a(4), a(5), a(6), v(1), &
                                      v(2), v(3), v(4), v(5), v(6), v(7), v(8), &
                                      v(9), v(10), v(11), v(12), v(13), v(14), v(15), &
                                      v(16), v(17), v(18), v(19), v(20))
            return
        end if
        call nested_sheared_field(a(1), a(2), a(3), a(4), a(5), a(6), a(7), a(8), &
                   v(1), v(2), v(3), v(4), v(5), v(6), v(7), v(8), v(9), v(10), v(11), &
                          v(12), v(13), v(14), v(15), v(16), v(17), v(18), v(19), v(20))
    end subroutine
    subroutine jvp_leaf(a, d, v)
        real(dp), intent(in)::a(:), d(:)
        real(dp), intent(out)::v(20)
        if (size(a) == 6) then
            call nested_integer_field_jvp(a(1), a(2), a(3), a(4), a(5), a(6), d(1), &
                                          d(2), d(3), d(4), d(5), d(6), v(1), v(2), &
                                          v(3), v(4), v(5), v(6), v(7), v(8), v(9), &
                                      v(10), v(11), v(12), v(13), v(14), v(15), v(16), &
                                          v(17), v(18), v(19), v(20))
            return
        end if
        call nested_sheared_field_jvp(a(1), a(2), a(3), a(4), a(5), a(6), a(7), a(8), &
               d(1), d(2), d(3), d(4), d(5), d(6), d(7), d(8), v(1), v(2), v(3), v(4), &
               v(5), v(6), v(7), v(8), v(9), v(10), v(11), v(12), v(13), v(14), v(15), &
                                      v(16), v(17), v(18), v(19), v(20))
    end subroutine
    subroutine vjp_leaf(a, b, v)
        real(dp), intent(in)::a(:), b(20)
        real(dp), intent(out)::v(:)
        if (size(a) == 6) then
            call nested_integer_field_vjp(a(1), a(2), a(3), a(4), a(5), a(6), b(1), &
                                          b(2), b(3), b(4), b(5), b(6), b(7), b(8), &
                                       b(9), b(10), b(11), b(12), b(13), b(14), b(15), &
                                        b(16), b(17), b(18), b(19), b(20), v(1), v(2), &
                                          v(3), v(4), v(5), v(6))
            return
        end if
        call nested_sheared_field_vjp(a(1), a(2), a(3), a(4), a(5), a(6), a(7), a(8), &
                   b(1), b(2), b(3), b(4), b(5), b(6), b(7), b(8), b(9), b(10), b(11), &
                        b(12), b(13), b(14), b(15), b(16), b(17), b(18), b(19), b(20), &
                                      v(1), v(2), v(3), v(4), v(5), v(6), v(7), v(8))
    end subroutine
    subroutine interval_leaf(a, v)
        type(interval_t), intent(in)::a(:)
        type(interval_t), intent(out)::v(20)
        if (size(a) == 6) then
          call nested_integer_field_interval(a(1), a(2), a(3), a(4), a(5), a(6), v(1), &
                                             v(2), v(3), v(4), v(5), v(6), v(7), v(8), &
                                       v(9), v(10), v(11), v(12), v(13), v(14), v(15), &
                                               v(16), v(17), v(18), v(19), v(20))
            return
        end if
    call nested_sheared_field_interval(a(1), a(2), a(3), a(4), a(5), a(6), a(7), a(8), &
                   v(1), v(2), v(3), v(4), v(5), v(6), v(7), v(8), v(9), v(10), v(11), &
                          v(12), v(13), v(14), v(15), v(16), v(17), v(18), v(19), v(20))
    end subroutine
end program
