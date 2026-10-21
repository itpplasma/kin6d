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
    use fortnum_ode, only: ode_solution_t
    implicit none
    type(nested_oracle_t) :: o, ode_model
    type(ode_solution_t) :: fieldline_solution
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
    o%eps = 0
    ! Independent symmetry limit: fixed conventional theta gives the same
    ! (R,Z) at every phi. A material-angle export instead winds twice.
    do asym = 0, 1
        o%tau = 0.4_dp*asym
        do n = 1, 17
            theta = 2*pi*n/17
            call o%boundary(theta, 0.0_dp, x, status)
            if (status /= 0) error stop 'axisymmetric boundary at phi=0'
            do k = 1, 13
                t = 2*pi*(k - 0.3_dp)/13
                call o%boundary(theta, t, xx, status)
                if (status /= 0) error stop 'axisymmetric boundary at phi'
                if (abs(hypot(x(1), x(2)) - hypot(xx(1), xx(2))) > 2e-13_dp &
                    .or. abs(x(3) - xx(3)) > 2e-13_dp) &
                    error stop 'poloidal angle must not wind in axisymmetric limit'
            end do
        end do
    end do
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
    call check_matched_profiles()
    call check_flux_derivatives()
    ! Axis transform from the independent published numerical value.
    o%family = NESTED_SHEARED; o%eps = 2; o%saxis = 1; o%delta = 0.005_dp; o%tau = 0
    call o%rotation(0.0_dp, 8, 1e-7_dp, axis, status)
    if (status /= 0 .or. abs(axis - 2.28690_dp) > 1e-5_dp) error stop 'axis iota oracle'
    call o%rotation(o%delta, 64, 1e-7_dp, lo, status)
    call o%rotation(o%delta, 128, 3e-8_dp, hi, status)
    if (status /= 0 .or. abs(lo - hi) > 1e-7_dp) error stop 'flux-ratio tolerance'
    print '(a,3f12.7)', '[V] sheared iota: axis,coarse,fine quadrature ', axis, lo, hi
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
    subroutine check_matched_profiles()
        use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
        use fortnum_ode, only: ode_problem_t, ode_workspace_t, ode_solution_t, &
                               ode_integrate
        use fortnum_status, only: fortnum_status_t
!$      use omp_lib, only: omp_set_num_threads
        type(nested_oracle_t) :: model
        type(ode_problem_t) :: problem
        type(ode_workspace_t) :: work
        type(fortnum_status_t) :: ode_status
        real(dp), parameter :: snodes(5) = [0.0_dp, 0.125_dp, 0.5_dp, 0.875_dp, 1.0_dp]
        ! gvecshear's independently integrated flux periods, rounded to 8 places.
        real(dp), parameter :: axis_values(2) = [2.79964725_dp, 2.64613001_dp]
        real(dp), parameter :: edge_values(2) = [2.81043883_dp, 2.63668554_dp]
        real(dp), parameter :: ode_values(2) = [2.8100648097671277_dp, &
                                                2.6369692267616878_dp]
        real(dp) :: labels(5), phis(5), transforms(5), pressures(5), serial(4, 5)
        real(dp) :: phiedge, coarse, tight, estimate, finite_turn, xyz(3), j, refv(5)
        real(dp) :: invalid_s(3), tmp(3), tmp2(3), tmp3(3), tmp4(3), nan
        integer :: codes(5), bad_codes(3), code, branch, k
        model%edge_pressure = 0.03_dp
        do branch = 1, 2
            model%tau = 0.4_dp*(branch - 1)
            model%family = NESTED_INTEGER
            do k = 1, 5
                call model%rotation(snodes(k)*model%delta, 1, 1e-10_dp, coarse, code)
                if (code /= 0 .or. coarse /= 2.0_dp) error stop 'integer exact iota'
            end do
        call model%profiles(snodes, 1e-8_dp, labels, phis, transforms, pressures, codes)
            if (any(codes /= 0)) error stop 'integer matched profiles'
            if (any(transforms /= 2.0_dp)) error stop 'integer profile exact iota'
            if (maxval(abs(labels - snodes*model%delta)) > 1e-15_dp) &
                error stop 'integer linear flux label'
            if (maxval(abs(pressures - model%edge_pressure &
                           - 2*(model%delta - labels)/(1 - model%tau**2))) > 1e-14_dp) &
                error stop 'integer profile pressure'
            model%family = NESTED_SHEARED
            model%eps = 1.2_dp; model%delta = 0.14_dp
            call model%flux(model%delta, 1e-9_dp, phiedge, code, estimate)
            if (code /= 0) error stop 'sheared flux convergence'
            if (abs(phiedge - 0.09622812651621984_dp) > 1e-10_dp) &
                error stop 'gvecshear matched toroidal flux'
            if (estimate > 1e-9_dp*phiedge) error stop 'flux tolerance budget'
!$          call omp_set_num_threads(1)
        call model%profiles(snodes, 1e-8_dp, labels, phis, transforms, pressures, codes)
            if (any(codes /= 0)) error stop 'sheared matched profiles'
            serial(1, :) = labels; serial(2, :) = phis
            serial(3, :) = transforms; serial(4, :) = pressures
!$          call omp_set_num_threads(4)
        call model%profiles(snodes, 1e-8_dp, labels, phis, transforms, pressures, codes)
            if (any(codes /= 0)) error stop 'parallel sheared profiles'
            if (maxval(abs(labels - serial(1, :))) > 1e-12_dp .or. &
                maxval(abs(phis - serial(2, :))) > 1e-12_dp .or. &
                maxval(abs(transforms - serial(3, :))) > 1e-12_dp .or. &
                maxval(abs(pressures - serial(4, :))) > 1e-12_dp) &
                error stop 'profile thread parity'
            if (abs(transforms(1) - axis_values(branch)) > 2e-8_dp .or. &
                abs(transforms(5) - edge_values(branch)) > 2e-8_dp) &
                error stop 'gvecshear axis/edge transforms'
            if (maxval(abs(labels - snodes*model%delta)) < 1e-5_dp) &
                error stop 'sheared toroidal flux is not psi/delta'
            do k = 1, 5
                call model%flux(labels(k), 1e-10_dp, tight, code)
                if (code /= 0 .or. abs(tight/phiedge - snodes(k)) > 1e-8_dp) &
                    error stop 'matched flux inverse residual'
                call model%surface(sqrt(2*labels(k)), 0.73_dp, 0.29_dp, xyz, j)
                call reference(model, xyz, refv)
                if (abs(refv(2) - pressures(k)) > 2e-10_dp) &
                    error stop 'matched pressure versus Cartesian field'
            end do
            tight = independent_toroidal_flux(model, model%delta)
            if (abs(tight - phiedge) > 2e-8_dp) error stop 'physical B dot area flux'
            call model%rotation(model%delta, 1, 1e-7_dp, coarse, code)
            if (code /= 0) error stop 'coarse flux-ratio convergence'
            call model%rotation(model%delta, 128, 1e-10_dp, tight, code, estimate)
            if (code /= 0 .or. abs(coarse - tight) > 1e-7_dp*abs(tight)) &
                error stop 'flux ratio within quadrature stopping tolerance'
            if (estimate > 1e-10_dp*abs(tight)) error stop 'iota tolerance budget'
            ode_model = model
            problem%rhs => independent_rhs
            problem%t1 = 256*pi
            problem%y0 = [0.0_dp]
            problem%rtol = 1e-10_dp; problem%atol = 1e-10_dp
            problem%hmax = 0.1_dp
            call ode_integrate(problem, work, fieldline_solution, ode_status)
            if (ode_status%code /= 0) error stop 'independent 128-turn integration'
            finite_turn = fieldline_solution%y(1, size(fieldline_solution%t))/(256*pi)
            if (abs(finite_turn - ode_values(branch)) > 3e-7_dp) &
                error stop 'independent 128-turn endpoint'
            ! Positive toroidal density maps each poloidal turn monotonically
            ! to a flux angle: its bounded endpoint correction is <1 turn.
            if (abs(tight - finite_turn) > 1.0_dp/128) &
                error stop '128-turn winding versus flux-period ratio'
            if (abs(tight - finite_turn) < 1e-4_dp) &
                error stop 'finite-turn endpoint must not define exact iota'
            print '(a,f4.1,4es18.10)', '[V] tau,Phi_edge,iota_axis,edge,128turn ', &
                model%tau, phiedge, transforms(1), transforms(5), finite_turn
            model%family = NESTED_INTEGER
            model%eps = 0.5_dp; model%delta = 1.0_dp/64
        end do
        model%family = NESTED_SHEARED
        model%eps = 1.2_dp; model%delta = 0.14_dp
        nan = ieee_value(0.0_dp, ieee_quiet_nan)
        call model%rotation(nan, 1, 1e-8_dp, tight, code)
        if (code == 0) error stop 'NaN flux label accepted'
        call model%rotation(0.1_dp, 1, nan, tight, code)
        if (code == 0) error stop 'NaN tolerance accepted'
        call model%rotation(0.1_dp, 1, 1e-20_dp, tight, code)
        if (code == 0) error stop 'unattainable quadrature tolerance accepted'
        invalid_s = [-0.1_dp, 1.1_dp, nan]
        call model%profiles(invalid_s, 1e-8_dp, tmp, tmp2, tmp3, tmp4, bad_codes)
        if (any(bad_codes == 0)) error stop 'invalid normalized flux accepted'
    end subroutine

    subroutine independent_rhs(time, state, rate, ctx)
        real(dp), intent(in) :: time, state(:)
        real(dp), intent(out) :: rate(:)
        class(*), intent(in), optional :: ctx
        real(dp) :: rr, y, nu, ch, sigma, h
        rr = sqrt(2*ode_model%delta)
        y = rr*sin(state(1)); nu = ode_model%eps*sin(2*time)/2 + ode_model%tau
        ch = cosh(nu)
        sigma = ode_model%saxis + atan(tanh(nu)*y/sqrt(1 - y*y)) &
                + asin(rr*cos(state(1))/sqrt(ch*ch - y*y))
        h = sqrt(4*sigma*sigma + ode_model%eps**2)
        rate(1) = (h + ode_model%eps*cos(2*time))*sqrt((1 - y*y)/(ch*ch - rr*rr))
        if (present(ctx)) continue
    end subroutine

    function independent_toroidal_flux(model, psi) result(value)
        type(nested_oracle_t), intent(in) :: model
        real(dp), intent(in) :: psi
        real(dp) :: value, radii(24), rw(24), r, angle, chart(3), xp(3), xm(3)
        real(dp) :: dr(3), dt(3), values(5), jac
        real(dp), parameter :: dh = 1e-5_dp
        integer :: ir, it
        call gauss_legendre_ab(24, 0.0_dp, sqrt(2*psi), radii, rw)
        value = 0
        do it = 1, 128
            angle = 2*pi*(it - 0.37_dp)/128
            do ir = 1, 24
                r = radii(ir)
                call model%surface(r, angle, 0.39_dp, chart, jac)
                call model%surface(r + dh, angle, 0.39_dp, xp, jac)
                call model%surface(r - dh, angle, 0.39_dp, xm, jac)
                dr = (xp - xm)/(2*dh)
                call model%surface(r, angle + dh, 0.39_dp, xp, jac)
                call model%surface(r, angle - dh, 0.39_dp, xm, jac)
                dt = (xp - xm)/(2*dh)
                call reference(model, chart, values)
               value = value + rw(ir)*(2*pi/128)*dot_product(values(3:5), cross(dr, dt))
            end do
        end do
    end function

    subroutine check_flux_derivatives()
        use kin6d_nested_oracle_leaves, only: nested_sheared_flux_density, &
                        nested_sheared_flux_density_jvp, nested_sheared_flux_density_vjp
        real(dp) :: a(7), d(7), b(2), jv(2), jt(7), plus(2), minus(2)
        a = [0.4_dp, 0.7_dp, 0.3_dp, 1.2_dp, 0.4_dp, 1.4_dp, 1.6_dp]
        d = [0.1_dp, -0.2_dp, 0.3_dp, 0.15_dp, 0.1_dp, -0.12_dp, 0.2_dp]
        b = [0.7_dp, -0.3_dp]
        call nested_sheared_flux_density_jvp(a(1), a(2), a(3), a(4), a(5), a(6), a(7), &
                                 d(1), d(2), d(3), d(4), d(5), d(6), d(7), jv(1), jv(2))
        call nested_sheared_flux_density_vjp(a(1), a(2), a(3), a(4), a(5), a(6), a(7), &
                            b(1), b(2), jt(1), jt(2), jt(3), jt(4), jt(5), jt(6), jt(7))
        if (abs(dot_product(b, jv) - dot_product(jt, d)) > 1e-12_dp) &
            error stop 'flux density JVP/VJP pairing'
        a = a + 1e-6_dp*d
        call nested_sheared_flux_density(a(1), a(2), a(3), a(4), a(5), a(6), a(7), &
                                         plus(1), plus(2))
        a = a - 2e-6_dp*d
        call nested_sheared_flux_density(a(1), a(2), a(3), a(4), a(5), a(6), a(7), &
                                         minus(1), minus(2))
        if (maxval(abs((plus - minus)/2e-6_dp - jv)) > 1e-8_dp) &
            error stop 'flux density independent directional difference'
        print '(a)', '[V] flux-density JVP/VJP and directional differences PASS'
    end subroutine

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
