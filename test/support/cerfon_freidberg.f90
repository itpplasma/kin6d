module cerfon_freidberg_fixture
    ! Published synthetic shapes only. No machine equilibrium data.
    ! Cerfon & Freidberg, Phys. Plasmas 17, 032502 (2010),
    ! https://doi.org/10.1063/1.3328818, eqs. (5)-(11), Secs. V-VI.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_cerfon_freidberg_leaves, only: cerfon_freidberg_jet, &
        cerfon_freidberg_jet_jvp, cerfon_freidberg_jet_vjp
    use kin6d_equilibrium_bulk_leaves, only: equilibrium_bulk_density, &
        equilibrium_bulk_density_jvp, equilibrium_bulk_density_vjp, &
        equilibrium_bulk_normalize, equilibrium_bulk_normalize_jvp, &
        equilibrium_bulk_normalize_vjp
    use fortnum_linalg, only: dense_solve
    use fortnum_quadrature, only: gauss_legendre_ab
    use fortnum_roots, only: root_brent
    use fortnum_status, only: fortnum_status_t, FORTNUM_OK

    use fortnum_roots, only: scalar_root_function_t

    use fortnum_roots, only: scalar_root_function_t
    implicit none
    private
    public :: cf_case_t, cf_setup, cf_jet, cf_surface, cf_oracle, cf_q
    real(dp), parameter :: pi = acos(-1.0_dp)
    type :: cf_case_t
        real(dp) :: eps, kappa, delta, a, c(7), axis, psi_axis
    end type
    type :: cf_setup_context_t
        type(cf_case_t), pointer :: model => null()
    end type

    type, extends(scalar_root_function_t) :: cf_setup_axis_derivative_t
        type(cf_setup_context_t), pointer :: data => null()
    contains
        procedure :: evaluate => cf_setup_axis_derivative_binding
    end type

    type :: cf_surface_context_t
        real(dp), pointer :: level => null()
        type(cf_case_t), pointer :: model => null()
        real(dp), pointer :: ray(:) => null()
    end type

    type, extends(scalar_root_function_t) :: cf_surface_residual_t
        type(cf_surface_context_t), pointer :: data => null()
    contains
        procedure :: evaluate => cf_surface_residual_binding
    end type

contains
    pure subroutine cf_jet(model, point, jet)
        type(cf_case_t), intent(in) :: model
        real(dp), intent(in) :: point(2)
        real(dp), intent(out) :: jet(6)
        call cerfon_freidberg_jet(point(1), point(2), model%a, model%c(1), &
            model%c(2), model%c(3), model%c(4), model%c(5), model%c(6), model%c(7), &
            jet(1), jet(2), jet(3), jet(4), jet(5), jet(6))
    end subroutine

    subroutine cf_setup(name, model)
        character(*), intent(in) :: name
        type(cf_case_t), intent(out) , target :: model
        type(cf_case_t) :: basis
        real(dp) :: matrix(7, 7), rhs(7), residual(7), jet(6), point(2)
        real(dp) :: alpha, curvature(3), step
        type(fortnum_status_t) :: root_status
        integer :: i, j, status
        type(cf_setup_context_t), target :: ctx
        type(cf_setup_axis_derivative_t), target :: axis_derivative_operator
        model%eps = 0.32_dp; model%kappa = 1.7_dp; model%delta = 0.33_dp
        model%a = 0
        select case (name)
        case ('iter-0')
        case ('iter-155')
            model%a = -0.155_dp ! Sec. V, Fig. 2.
        case ('nstx-0', 'nstx-1')
            model%eps = 0.78_dp; model%kappa = 2; model%delta = 0.35_dp
            if (name == 'nstx-1') model%a = 1 ! Sec. VI, Fig. 3(a,b).
        case default
            error stop 'unknown Cerfon-Freidberg case'
        end select
        alpha = asin(model%delta)
        ! Eq. (11). Curvatures of the Miller model constrain psi, but that
        ! model is NOT used as the boundary of the exact GS solution.
        curvature = [-(1 + alpha)**2/(model%eps*model%kappa**2), &
            (1 - alpha)**2/(model%eps*model%kappa**2), &
            -model%kappa/(model%eps*cos(alpha)**2)]
        model%c = 0
        call constraints(model, rhs)
        do j = 1, 7
            basis = model
            basis%c(j) = 1
            call constraints(basis, residual)
            matrix(:, j) = residual - rhs
        end do
        call dense_solve(matrix, -rhs, model%c, status)
        if (status /= 0) error stop 'CF seven boundary conditions dense solve'
        call constraints(model, residual)
        if (maxval(abs(residual)) > 2e-11_dp) error stop 'CF boundary residual'
        ! Up-down symmetry makes the axis lie on y=0. Bracket the derivative;
        ! no KIN6D axis search or readback enters this oracle.
        call bind_cf_setup_context()
        call root_brent(axis_derivative_operator, 1 - model%eps, 1 + model%eps, &
            model%axis, root_status, 2e-14_dp, 1e-15_dp, 128)
        if (root_status%code /= FORTNUM_OK) error stop 'CF exact axis root'
        call cf_jet(model, [model%axis, 0.0_dp], jet)
        model%psi_axis = jet(1)
        if (model%psi_axis >= 0 .or. min(jet(4), jet(6)) <= 0) &
            error stop 'CF regular minimum'
        ! Structural roundoff check of the GENERATED jet, away from the axis
        ! and the interpolation points (the generator also proves eq. (5)).
        do i = 1, 31
            point = [0.25_dp + 1.5_dp*i/31, -0.7_dp + 1.4_dp*i/31]
            call cf_jet(model, point, jet)
            if (abs(jet(4) - jet(2)/point(1) + jet(6) &
                - ((1 - model%a)*point(1)**2 + model%a)) > 1e-11_dp) &
                error stop 'CF generated GS residual'
            ! Independent centred directional derivative, not an AD pairing.
            step = 2e-5_dp
            call check_derivative(point, jet, step)
        end do
        call check_generated_products(model)
    contains
    subroutine bind_cf_setup_context()
        ctx%model => model
        axis_derivative_operator%data => ctx
    end subroutine bind_cf_setup_context

        subroutine constraints(m, result)
            type(cf_case_t), intent(in) :: m
            real(dp), intent(out) :: result(7)
            real(dp) :: outer(6), inner(6), top(6)
            call cf_jet(m, [1 + m%eps, 0.0_dp], outer)
            call cf_jet(m, [1 - m%eps, 0.0_dp], inner)
            call cf_jet(m, [1 - m%delta*m%eps, m%kappa*m%eps], top)
            ! Eq. (10): three zeros, horizontal top, three curvatures.
            result = [outer(1), inner(1), top(1), top(2), &
                outer(6) + curvature(1)*outer(2), &
                inner(6) + curvature(2)*inner(2), &
                top(4) + curvature(3)*top(3)]
        end subroutine
        subroutine check_derivative(point, jet, step)
            real(dp), intent(in) :: point(2), jet(6), step
            real(dp) :: plus(6), minus(6), direction(2)
            direction = [0.6_dp, -0.8_dp]
            call cf_jet(model, point + step*direction, plus)
            call cf_jet(model, point - step*direction, minus)
            if (abs((plus(1) - minus(1))/(2*step) &
                - dot_product(jet(2:3), direction)) > 2e-8_dp) &
                error stop 'CF directional derivative'
            if (norm2((plus(2:3) - minus(2:3))/(2*step) &
                - matmul(reshape([jet(4), jet(5), jet(5), jet(6)], [2, 2]), &
                direction)) > 2e-7_dp) error stop 'CF Hessian derivative'
        end subroutine
    end subroutine

    subroutine check_generated_products(model)
        type(cf_case_t), intent(in) :: model
        real(dp), parameter :: h = 1e-5_dp
        real(dp) :: x(10), d(10), b(6), tangent(6), adjoint(10), fd(6)
        real(dp) :: u(6), du(6), weight(2), out(2), bar(6), difference(2)
        real(dp) :: v(4), dv(4), gv(4)
        integer :: i
        x = [1.1_dp, 0.2_dp, model%a, model%c]
        d = [(real(i, dp)/100, i=1, 10)]
        b = [(real(4 - i, dp)/7, i=1, 6)]
        call cerfon_freidberg_jet_jvp(x(1), x(2), x(3), x(4), x(5), x(6), x(7), &
            x(8), x(9), x(10), d(1), d(2), d(3), d(4), d(5), d(6), d(7), d(8), &
            d(9), d(10), tangent(1), tangent(2), tangent(3), tangent(4), &
            tangent(5), tangent(6))
        call cerfon_freidberg_jet_vjp(x(1), x(2), x(3), x(4), x(5), x(6), x(7), &
            x(8), x(9), x(10), b(1), b(2), b(3), b(4), b(5), b(6), adjoint(1), &
            adjoint(2), adjoint(3), adjoint(4), adjoint(5), adjoint(6), adjoint(7), &
            adjoint(8), adjoint(9), adjoint(10))
        fd = (jet_array(x + h*d) - jet_array(x - h*d))/(2*h)
        call check_products(fd, tangent, d, adjoint, b)
        u = [3.0_dp, 0.2_dp, 0.1_dp, -0.5_dp, 3.0_dp, 1.0_dp]
        du = [0.1_dp, -0.2_dp, 0.3_dp, 0.2_dp, -0.1_dp, 0.4_dp]
        weight = [0.7_dp, -0.2_dp]
        call equilibrium_bulk_normalize_jvp(u(1), u(2), u(3), u(4), u(5), u(6), &
            du(1), du(2), du(3), du(4), du(5), du(6), out(1), out(2))
        call equilibrium_bulk_normalize_vjp(u(1), u(2), u(3), u(4), u(5), u(6), &
            weight(1), weight(2), bar(1), bar(2), bar(3), bar(4), bar(5), bar(6))
        difference = (bulk_array(u + h*du) - bulk_array(u - h*du))/(2*h)
        call check_products(difference, out, du, bar, weight)
        v = [1.3_dp, 0.2_dp, -0.1_dp, 0.4_dp]
        dv = [0.2_dp, -0.3_dp, 0.1_dp, 0.2_dp]
        call equilibrium_bulk_density_jvp(v(1), v(2), v(3), v(4), dv(1), dv(2), &
            dv(3), dv(4), out(1), out(2))
        call equilibrium_bulk_density_vjp(v(1), v(2), v(3), v(4), weight(1), &
            weight(2), gv(1), gv(2), gv(3), gv(4))
        difference = (density_array(v + h*dv) - density_array(v - h*dv))/(2*h)
        call check_products(difference, out, dv, gv, weight)
    end subroutine

    subroutine check_products(fd, tangent, direction, adjoint, weight)
        real(dp), intent(in) :: fd(:), tangent(:), direction(:), adjoint(:), weight(:)
        if (norm2(fd - tangent) > 1e-7_dp*max(1.0_dp, norm2(tangent))) &
            error stop 'CF/bulk generated JVP finite difference'
        if (abs(dot_product(tangent, weight) - dot_product(direction, adjoint)) &
            > 2e-12_dp*max(1.0_dp, norm2(tangent)*norm2(weight))) &
            error stop 'CF/bulk generated transpose pairing'
    end subroutine

    pure function jet_array(x) result(y)
        real(dp), intent(in) :: x(10)
        real(dp) :: y(6)
        call cerfon_freidberg_jet(x(1), x(2), x(3), x(4), x(5), x(6), x(7), x(8), &
            x(9), x(10), y(1), y(2), y(3), y(4), y(5), y(6))
    end function

    pure function bulk_array(x) result(y)
        real(dp), intent(in) :: x(6)
        real(dp) :: y(2)
        call equilibrium_bulk_normalize(x(1), x(2), x(3), x(4), x(5), x(6), y(1), &
            y(2))
    end function

    pure function density_array(x) result(y)
        real(dp), intent(in) :: x(4)
        real(dp) :: y(2)
        call equilibrium_bulk_density(x(1), x(2), x(3), x(4), y(1), y(2))
    end function

    subroutine cf_surface(model, angle, level, point, radius, slope)
        type(cf_case_t), intent(in) , target :: model
        real(dp), intent(in) , target :: angle, level
        real(dp), intent(out) :: point(2), radius, slope
        real(dp) , target :: ray(2), limit, lo, hi, value, jet(6)
        type(fortnum_status_t) :: status
        integer :: k
        type(cf_surface_context_t), target :: ctx
        type(cf_surface_residual_t), target :: residual_operator
        ray = model%eps*[cos(angle), model%kappa*sin(angle)]
        limit = 4
        if (ray(1) < 0) limit = min(limit, -0.999_dp*model%axis/ray(1))
        lo = 0
        hi = 0
        ! First crossing is the closed surface around this magnetic axis.
        do k = 1, 256
            hi = limit*k/256
        call bind_cf_surface_context()
            value = cf_surface_residual(ctx, hi)
            if (value >= 0) exit
            lo = hi
        end do
        if (value < 0) error stop 'CF surface bracket'
        call bind_cf_surface_context()
        call root_brent(residual_operator, lo, hi, radius, status, 2e-14_dp, &
            1e-15_dp, 128)
        if (status%code /= FORTNUM_OK) error stop 'CF surface root'
        point = [model%axis, 0.0_dp] + radius*ray
        call cf_jet(model, point, jet)
        slope = dot_product(jet(2:3), ray)
        if (slope <= 0) error stop 'CF non-star-shaped contour'
    contains
    subroutine bind_cf_surface_context()
        ctx%level => level
        ctx%model => model
        ctx%ray => ray
        residual_operator%data => ctx
    end subroutine bind_cf_surface_context

    end subroutine

    function cf_q(model, psi_n, order) result(q)
        type(cf_case_t), intent(in) :: model
        real(dp), intent(in) :: psi_n
        integer, intent(in) :: order
        real(dp) :: q, theta(order), weight(order), point(2), radius, slope, level
        integer :: i
        level = model%psi_axis*(1 - psi_n)
        call gauss_legendre_ab(order, 0.0_dp, 2*pi, theta, weight)
        q = 0
        do i = 1, order
            call cf_surface(model, theta(i), level, point, radius, slope)
            ! Direct delta-function area quadrature of the exact level set:
            ! q = F/(2pi) int s eps^2 kappa/(R dpsi/ds) dtheta.
            q = q + weight(i)*radius/(point(1)*slope)
        end do
        q = sqrt(1 - 2*model%a*level)*model%eps**2*model%kappa*q/(2*pi)
    end function

    subroutine cf_oracle(model, order, moments)
        type(cf_case_t), intent(in) :: model
        integer, intent(in) :: order
        ! volume, beta_p, li, current, perimeter, toroidal flux.
        real(dp), intent(out) :: moments(6)
        real(dp) :: theta(order), wt(order), s(order/4), ws(order/4), jet(6), point(2)
        real(dp) :: ray(2), ray_theta(2), sb, slope, sb_theta, jac, pressure
        real(dp) :: p_int, b_int, area_weight, vol_weight, length, current, phi
        real(dp) :: r, bp2
        integer :: i, j
        call gauss_legendre_ab(order, 0.0_dp, 2*pi, theta, wt)
        call gauss_legendre_ab(order/4, 0.0_dp, 1.0_dp, s, ws)
        moments = 0; p_int = 0; b_int = 0; length = 0; current = 0; phi = 0
        do i = 1, order
            ray = model%eps*[cos(theta(i)), model%kappa*sin(theta(i))]
            ray_theta = model%eps*[-sin(theta(i)), model%kappa*cos(theta(i))]
            call cf_surface(model, theta(i), 0.0_dp, point, sb, slope)
            call cf_jet(model, point, jet)
            sb_theta = -sb*dot_product(jet(2:3), ray_theta)/slope
            length = length + wt(i)*norm2(sb_theta*ray + sb*ray_theta)
            do j = 1, size(s)
                point = [model%axis, 0.0_dp] + sb*s(j)*ray
                call cf_jet(model, point, jet)
                r = point(1)
                jac = sb**2*s(j)*model%eps**2*model%kappa
                area_weight = wt(i)*ws(j)*jac
                vol_weight = 2*pi*r*area_weight
                pressure = (model%a - 1)*jet(1) ! mu0=R0=Psi0=F_edge=1.
                bp2 = sum(jet(2:3)**2)/r**2
                moments(1) = moments(1) + vol_weight
                p_int = p_int + vol_weight*pressure
                b_int = b_int + vol_weight*bp2
                current = current - area_weight*((1 - model%a)*r + model%a/r)
                phi = phi + area_weight*sqrt(1 - 2*model%a*jet(1))/r
            end do
        end do
        moments(2) = 2*p_int/moments(1)/(current/length)**2
        moments(3) = b_int/moments(1)/(current/length)**2
        moments(4:6) = [current, length, phi]
    end subroutine

    pure function cf_setup_axis_derivative(ctx, r) result(value)
        type(cf_setup_context_t), intent(in) :: ctx
            real(dp), intent(in) :: r
            real(dp) :: value, local(6)
            call cf_jet(ctx%model, [r, 0.0_dp], local)
            value = local(2)
    end function cf_setup_axis_derivative

    pure function cf_setup_axis_derivative_binding(self, x) result(result_value)
        class(cf_setup_axis_derivative_t), intent(in) :: self
            real(dp), intent(in) :: x
            real(dp) :: result_value, local(6)
        result_value = cf_setup_axis_derivative(self%data, x)
    end function

    pure function cf_surface_residual(ctx, s) result(value)
        type(cf_surface_context_t), intent(in) :: ctx
            real(dp), intent(in) :: s
            real(dp) :: value, local(6)
            call cf_jet(ctx%model, [ctx%model%axis, 0.0_dp] + s*ctx%ray, local)
            value = local(1) - ctx%level
    end function cf_surface_residual

    pure function cf_surface_residual_binding(self, x) result(result_value)
        class(cf_surface_residual_t), intent(in) :: self
            real(dp), intent(in) :: x
            real(dp) :: result_value, local(6)
        result_value = cf_surface_residual(self%data, x)
    end function
end module
