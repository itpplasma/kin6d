module kin6d_nested_oracles
    ! Exact nested ideal-MHD references, mu0=1. No solver uses these as forcing.
    ! Local mathematics (including derivatives) is in gen/gen_nested_oracles.f90.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
    use kin6d_nested_oracle_leaves, only: nested_integer_field, &
                 nested_sheared_field, nested_integer_surface, nested_sheared_surface, &
                                          nested_setup, nested_integer_domain, &
     nested_ray, nested_reduce, nested_rotation_number, nested_integer_field_interval, &
                            nested_sheared_field_interval, nested_pressure_axis_interval
    use fortnum_interval, only: interval_t, interval
    use fortnum_quadrature, only: gauss_legendre_ab
    use fortnum_ode, only: ode_problem_t, ode_workspace_t, ode_solution_t, ode_integrate
    use fortnum_status, only: fortnum_status_t
    implicit none
    private
    public :: nested_oracle_t, NESTED_INTEGER, NESTED_SHEARED
    integer, parameter :: NESTED_INTEGER = 1, NESTED_SHEARED = 2
    real(dp), parameter :: pi = acos(-1.0_dp)
    type :: nested_oracle_t
        integer :: family = NESTED_INTEGER
      real(dp) :: eps = 0.5_dp, tau = 0, delta = 1.0_dp/64, saxis = 1.4_dp, lam = 1.6_dp
        real(dp) :: edge_pressure = 0
    contains
        procedure, non_overridable :: valid
        procedure, non_overridable :: point
        procedure, non_overridable :: point_enclosure
        procedure, non_overridable :: sample
        procedure, non_overridable :: surface
        procedure, non_overridable :: boundary
        procedure, non_overridable :: averages
        procedure, non_overridable :: rotation
    end type
contains
    logical function valid(o)
        class(nested_oracle_t), intent(in) :: o
        real(dp) :: limit, extent, pa1, pa2
        valid = .false.
        if (.not. all(ieee_is_finite([o%eps, o%tau, o%delta, o%saxis, o%lam, &
                                      o%edge_pressure]))) return
        if (o%delta <= 0 .or. o%lam <= 0 .or. o%eps < 0) return
        select case (o%family)
        case (NESTED_INTEGER)
            if (o%eps >= 1 .or. abs(o%tau) >= 1) return
            call nested_setup(o%eps, o%tau, o%delta, o%lam, o%edge_pressure, &
                              limit, extent, pa1, pa2)
            if (o%delta >= limit) return
        case (NESTED_SHEARED)
            if (o%delta >= 0.5_dp) return
            call nested_setup(o%eps, 0.0_dp, o%delta, o%lam, o%edge_pressure, &
                              limit, extent, pa1, pa2)
            if (o%saxis <= extent) return
        case default
            return
        end select
        valid = .true.
    end function

    subroutine point(o, x, v, status)
        class(nested_oracle_t), intent(in) :: o
        real(dp), intent(in) :: x(3)
        ! psi,p,B(3),curl B(3),grad p(3),column-major grad B(3,3).
        real(dp), intent(out) :: v(20)
        integer, intent(out) :: status
        real(dp) :: s, f2, limit, extent, pa1, pa2
        status = 1
        v = ieee_value(0.0_dp, ieee_quiet_nan)
        if (.not. all(ieee_is_finite(x))) return
        if (o%family == NESTED_INTEGER) then
            call nested_integer_domain(x(1), x(2), x(3), o%eps, o%tau, s, f2)
            if (s <= 0 .or. f2 <= 0) return
        else
            if (x(1) == 0 .and. x(2)**2 <= o%eps) return
        end if
        call nested_setup(o%eps, merge(o%tau, 0.0_dp, o%family == NESTED_INTEGER), &
                          o%delta, o%lam, o%edge_pressure, limit, extent, pa1, pa2)
        if (o%family == NESTED_INTEGER) then
            call nested_integer_field(x(1), x(2), x(3), o%eps, o%tau, pa1, &
                          v(1), v(2), v(3), v(4), v(5), v(6), v(7), v(8), v(9), v(10), &
                   v(11), v(12), v(13), v(14), v(15), v(16), v(17), v(18), v(19), v(20))
        else
            call nested_sheared_field(x(1), x(2), x(3), o%eps, o%tau, o%saxis, &
                     o%lam, pa2, v(1), v(2), v(3), v(4), v(5), v(6), v(7), v(8), v(9), &
            v(10), v(11), v(12), v(13), v(14), v(15), v(16), v(17), v(18), v(19), v(20))
        end if
        if (all(ieee_is_finite(v))) status = 0
    end subroutine

    subroutine point_enclosure(o, x, v)
        class(nested_oracle_t), intent(in) :: o
        ! Encloses evaluation at the exact stored Cartesian point and parameters.
        ! Caller must first admit the point; this is not an enclosure of the map.
        real(dp), intent(in) :: x(3)
        type(interval_t), intent(out) :: v(20)
        type(interval_t) :: pa, pai1, pai2
        call nested_pressure_axis_interval( &
        interval(merge(o%tau, 0.0_dp, o%family == NESTED_INTEGER)), interval(o%delta), &
            interval(o%lam), interval(o%edge_pressure), pai1, pai2)
        if (o%family == NESTED_INTEGER) then
            pa = pai1
        else
            pa = pai2
        end if
        if (o%family == NESTED_INTEGER) then
            call nested_integer_field_interval(interval(x(1)), interval(x(2)), &
                                 interval(x(3)), interval(o%eps), interval(o%tau), pa, &
                          v(1), v(2), v(3), v(4), v(5), v(6), v(7), v(8), v(9), v(10), &
                   v(11), v(12), v(13), v(14), v(15), v(16), v(17), v(18), v(19), v(20))
        else
            call nested_sheared_field_interval(interval(x(1)), interval(x(2)), &
                  interval(x(3)), interval(o%eps), interval(o%tau), interval(o%saxis), &
            interval(o%lam), pa, v(1), v(2), v(3), v(4), v(5), v(6), v(7), v(8), v(9), &
            v(10), v(11), v(12), v(13), v(14), v(15), v(16), v(17), v(18), v(19), v(20))
        end if
    end subroutine

    subroutine sample(o, x, v, status)
        class(nested_oracle_t), intent(in) :: o
        real(dp), intent(in) :: x(:, :)
        real(dp), intent(out) :: v(:, :)
        integer, intent(out) :: status(:)
        integer :: k
        if (.not. o%valid()) error stop 'nested_oracles: invalid parameters'
        if (size(x, 1) /= 3 .or. size(v, 1) /= 20) error stop 'nested_oracles: shape'
        if (size(v, 2) /= size(x, 2) .or. size(status) /= size(x, 2)) &
            error stop 'nested_oracles: batch length'
        !$omp parallel do schedule(static)
        do k = 1, size(x, 2)
            call o%point(x(:, k), v(:, k), status(k))
        end do
        !$omp end parallel do
    end subroutine

    pure subroutine surface(o, rho, theta, t, x, jac, rhs)
        class(nested_oracle_t), intent(in) :: o
        ! rho^2=psi (integer); rho^2=2 psi (sheared). t is NOT geometric phi.
        real(dp), intent(in) :: rho, theta, t
        real(dp), intent(out) :: x(3), jac
        real(dp), intent(out), optional :: rhs
        real(dp) :: rate
        if (o%family == NESTED_INTEGER) then
         call nested_integer_surface(rho, theta, t, o%eps, o%tau, x(1), x(2), x(3), jac)
            rate = 2
        else
            call nested_sheared_surface(rho, theta, t, o%eps, o%tau, o%saxis, o%lam, &
                                        x(1), x(2), x(3), jac, rate)
        end if
        if (present(rhs)) rhs = rate
    end subroutine

    subroutine boundary(o, theta, phi, x, status)
        class(nested_oracle_t), intent(in) :: o
        real(dp), intent(in) :: theta, phi
        real(dp), intent(out) :: x(3)
        integer, intent(out) :: status
        real(dp) :: lo, hi, t, rho, angle, jac, ray, radius
        integer :: k
        status = 1
        if (.not. o%valid()) return
        rho = sqrt(o%delta)
        angle = -theta ! Clockwise poloidal angle, positive geometric phi.
        if (o%family == NESTED_SHEARED) then
            rho = sqrt(2*o%delta)
            angle = theta
        end if
        lo = phi - pi/2; hi = phi + pi/2
        do k = 1, 55
            t = (lo + hi)/2
            call o%surface(rho, angle, t, x, jac)
            call nested_ray(x(1), x(2), phi, ray, radius)
            if (ray < 0) then
                lo = t
            else
                hi = t
            end if
        end do
        if (radius <= 0 .or. abs(ray) > 1e-12_dp) return
        status = 0
    end subroutine

    subroutine averages(o, n, values)
        class(nested_oracle_t), intent(in) :: o
        integer, intent(in) :: n
        ! V,<p>,<B^2>,beta. Quadrature values: no continuous enclosure.
        real(dp), intent(out) :: values(4)
        real(dp), allocatable :: r(:), w(:)
        real(dp) :: theta, t, x(3), jac, v(20), rho, acc(3), tmp(3), factor
        integer :: ir, it, ip, status, bad
        if (.not. o%valid() .or. n < 2) error stop 'nested_oracles: quadrature input'
        rho = sqrt(o%delta)
        if (o%family == NESTED_SHEARED) rho = sqrt(2*o%delta)
        allocate (r(n), w(n))
        call gauss_legendre_ab(n, 0.0_dp, rho, r, w)
        acc = 0; bad = 0
        !$omp parallel do collapse(2) schedule(static) &
        !$omp private(theta,t,x,jac,v,tmp,factor,status,ir) reduction(+:acc,bad)
        do ip = 0, 2*n - 1
            do it = 0, 2*n - 1
                theta = pi*it/n; t = pi*ip/n
                do ir = 1, n
                    call o%surface(r(ir), theta, t, x, jac)
                    call o%point(x, v, status)
                    bad = bad + status
                    factor = w(ir)*(pi/n)**2*jac
              call nested_reduce(v(2), v(3), v(4), v(5), factor, tmp(1), tmp(2), tmp(3))
                    acc = acc + tmp
                end do
            end do
        end do
        !$omp end parallel do
        if (bad /= 0) error stop 'nested_oracles: inadmissible quadrature point'
        values(1) = acc(1)
        values(2:3) = acc(2:3)/acc(1)
        values(4) = 2*acc(2)/acc(3)
    end subroutine

    subroutine rotation(o, psi, turns, tolerance, iota, status)
        class(nested_oracle_t), intent(in) :: o
        real(dp), intent(in) :: psi, tolerance
        integer, intent(in) :: turns
        real(dp), intent(out) :: iota
        integer, intent(out) :: status
        type(ode_problem_t) :: problem
        type(ode_workspace_t) :: work
        type(ode_solution_t) :: solution
        type(fortnum_status_t) :: ode_status
        real(dp) :: rho
        status = 1
        iota = ieee_value(0.0_dp, ieee_quiet_nan)
        if (.not. o%valid()) return
        if (psi < 0 .or. psi > o%delta .or. turns < 1 .or. tolerance <= 0) return
        if (o%family == NESTED_INTEGER) then
            iota = 2; status = 0
            return
        end if
        rho = sqrt(2*psi)
        problem%rhs => rhs
        problem%t1 = 2*pi*turns
        problem%y0 = [0.0_dp]
        problem%rtol = tolerance; problem%atol = tolerance
        problem%hmax = 0.1_dp
        call ode_integrate(problem, work, solution, ode_status)
        if (ode_status%code /= 0) return
        call nested_rotation_number(solution%y(1, size(solution%t)), &
                                    real(turns, dp), iota)
        status = 0
    contains
        subroutine rhs(t, y, dy, ctx)
            real(dp), intent(in) :: t, y(:)
            real(dp), intent(out) :: dy(:)
            class(*), intent(in), optional :: ctx
            real(dp) :: x(3), jac
            call o%surface(rho, y(1), t, x, jac, dy(1))
            ! Context is intentionally host-associated, with no global state.
            if (present(ctx)) continue
        end subroutine
    end subroutine
end module
