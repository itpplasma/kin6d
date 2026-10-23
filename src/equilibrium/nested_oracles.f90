module kin6d_nested_oracles
    ! Exact nested ideal-MHD references, mu0=1. No solver uses these as forcing.
    ! Local mathematics (including derivatives) is in gen/gen_nested_oracles.f90.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
    use kin6d_nested_oracle_leaves, only: nested_integer_field, &
                 nested_sheared_field, nested_integer_surface, nested_sheared_surface, &
                                          nested_setup, nested_integer_domain, &
     nested_ray, nested_reduce, nested_integer_field_interval, &
                            nested_sheared_field_interval, nested_pressure_axis_interval, &
                            nested_sheared_flux_density, nested_integer_flux, &
                            nested_flux_ratio, nested_integer_pressure, nested_sheared_pressure, &
                            nested_profile_residual
    use fortnum_interval, only: interval_t, interval
    use fortnum_quadrature, only: gauss_legendre_ab
    use fortnum_integrate_gk, only: integrate_gk
    use fortnum_multiroot, only: multiroot_hybrid
    use fortnum_status, only: fortnum_status_t
    use fortnum_integrate_gk, only: gk_function_t
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
        procedure, non_overridable :: flux
        procedure, non_overridable :: profiles
    end type
    type, extends(gk_function_t) :: flux_periods_root_t
        integer, pointer :: k => null()
        class(nested_oracle_t), pointer :: o => null()
        real(dp), pointer :: rho => null()
    contains
        procedure :: evaluate => flux_periods_root_evaluate
    end type
    type :: profile_at_s_residual_context_t
        real(dp), pointer :: edge => null()
        integer, pointer :: failed => null()
        class(nested_oracle_t), pointer :: o => null()
        real(dp), pointer :: root(:) => null()
        real(dp), pointer :: s => null()
        real(dp), pointer :: tolerance => null()
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
        ! The integer family's material angle winds twice relative to the
        ! conventional poloidal angle. Untwist before Fourier export so
        ! iota=2 and an axisymmetric limit refer to the same coordinate.
        angle = 2*phi - theta
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
    subroutine rotation(o, psi, turns, tolerance, iota, status, error_estimate)
        class(nested_oracle_t), intent(in) :: o
        real(dp), intent(in) :: psi, tolerance
        ! Retained for existing callers; the flux ratio has no finite-turn bias.
        integer, intent(in) :: turns
        real(dp), intent(out) :: iota
        integer, intent(out) :: status
        real(dp), intent(out), optional :: error_estimate
        real(dp) :: periods(2), errors(2), estimate
        status = 1
        iota = ieee_value(0.0_dp, ieee_quiet_nan)
        if (present(error_estimate)) error_estimate = huge(1.0_dp)
        if (.not. admitted_profile(o, psi, tolerance) .or. turns < 1) return
        if (o%family == NESTED_INTEGER) then
            iota = 2; status = 0
            if (present(error_estimate)) error_estimate = 0
            return
        end if
        call flux_periods(o, psi, tolerance/8, periods, errors, status)
        if (status /= 0) return
        call nested_flux_ratio(periods(2), periods(1), iota)
        estimate = (errors(2) + abs(iota)*errors(1))/periods(1)
        if (present(error_estimate)) error_estimate = estimate
        status = 2
        if (.not. ieee_is_finite(iota)) return
        if (estimate <= tolerance*max(1.0_dp, abs(iota))) status = 0
    end subroutine
    logical function admitted_profile(o, psi, tolerance) result(ok)
        class(nested_oracle_t), intent(in) :: o
        real(dp), intent(in) :: psi, tolerance
        ok = .false.
        if (.not. o%valid()) return
        if (.not. all(ieee_is_finite([psi, tolerance]))) return
        if (psi < 0 .or. psi > o%delta .or. tolerance <= 0) return
        ok = .true.
    end function
    subroutine flux_periods(o, psi, tolerance, periods, errors, status)
        class(nested_oracle_t), intent(in), target :: o
        real(dp), intent(in) :: psi, tolerance
        ! dPhi_tor/dpsi, dpsi_pol/dpsi. Smooth axis limits need no rho division.
        real(dp), intent(out) :: periods(2), errors(2)
        integer, intent(out) :: status
        real(dp), target :: rho
        integer , target :: k, ierr
        type(flux_periods_root_t) :: root_function
        root_function%k => k
        root_function%o => o
        root_function%rho => rho
        rho = sqrt(2*psi)
        status = 2
        do k = 1, 2
            call integrate_gk(root_function, 0.0_dp, 2*pi, 0.0_dp, tolerance, &
                              periods(k), errors(k), ierr)
            if (ierr /= 0) return
        end do
        if (.not. all(ieee_is_finite(periods))) return
        if (any(periods <= 0)) return
        status = 0
    end subroutine
    subroutine flux(o, psi, tolerance, phi, status, error_estimate)
        class(nested_oracle_t), intent(in) :: o
        real(dp), intent(in) :: psi, tolerance
        real(dp), intent(out) :: phi
        integer, intent(out) :: status
        ! Relative quadrature tolerance; estimates exclude floating rounding.
        ! No continuous enclosure is available. Status 0 means E1 convergence,
        ! 1 invalid input, 2 quadrature failure (including a roundoff floor).
        real(dp), intent(out), optional :: error_estimate
        integer, parameter :: max_radial = 128
        real(dp) :: nodes(max_radial), weights(max_radial)
        real(dp) :: periods(2), errors(2), total, previous, estimate, angular
        integer :: n, j, ierr
        phi = ieee_value(0.0_dp, ieee_quiet_nan)
        status = 1
        if (present(error_estimate)) error_estimate = huge(1.0_dp)
        if (.not. admitted_profile(o, psi, tolerance)) return
        if (o%family == NESTED_INTEGER) then
            call nested_integer_flux(psi, o%eps, o%tau, phi)
            status = 0
            if (present(error_estimate)) error_estimate = 0
            return
        end if
        if (psi == 0) then
            phi = 0; status = 0
            if (present(error_estimate)) error_estimate = 0
            return
        end if
        previous = 0
        n = 8
        status = 2
        do
            call gauss_legendre_ab(n, 0.0_dp, psi, nodes(:n), weights(:n))
            total = 0; angular = 0
            do j = 1, n
                call flux_periods(o, nodes(j), tolerance/8, periods, errors, ierr)
                if (ierr /= 0) return
                total = total + weights(j)*periods(1)
                angular = angular + weights(j)*errors(1)
            end do
            estimate = abs(total - previous) + angular
            if (n > 8 .and. estimate <= tolerance*abs(total)) then
                phi = total; status = 0
                if (present(error_estimate)) error_estimate = estimate
                return
            end if
            if (n == max_radial) return
            previous = total
            n = 2*n
        end do
    end subroutine
    subroutine profiles(o, s, tolerance, psi, phi, iota, pressure, status)
        class(nested_oracle_t), intent(in) :: o
        ! Matched profiles at normalized TOROIDAL flux s=Phi(psi)/Phi(delta).
        ! mu0=1; Phi is the full positive toroidal flux (not flux/2pi).
        ! Shape-preserving batches; no interpolation and no guessed psi=delta*s.
        ! E1 quadrature/root convergence, continuous bounds unavailable.
        real(dp), intent(in) :: s(:), tolerance
        real(dp), intent(out) :: psi(:), phi(:), iota(:), pressure(:)
        integer, intent(out) :: status(:)
        real(dp) :: edge
        integer :: k, ierr
        if (size(psi) /= size(s) .or. size(phi) /= size(s) .or. &
            size(iota) /= size(s) .or. size(pressure) /= size(s) .or. &
            size(status) /= size(s)) error stop 'nested_oracles: profile shape'
        psi = ieee_value(0.0_dp, ieee_quiet_nan)
        phi = psi; iota = psi; pressure = psi
        status = 1
        call o%flux(o%delta, tolerance/16, edge, ierr)
        if (ierr /= 0) then
            status = ierr
            return
        end if
        !$omp parallel do schedule(static)
        do k = 1, size(s)
            call profile_at_s(o, s(k), tolerance, edge, psi(k), phi(k), &
                              iota(k), pressure(k), status(k))
        end do
        !$omp end parallel do
    end subroutine
    subroutine profile_at_s(o, s, tolerance, edge, psi, phi, iota, pressure, status)
        class(nested_oracle_t), intent(in), target :: o
        real(dp), intent(in), target :: s, tolerance, edge
        real(dp), intent(out) :: psi, phi, iota, pressure
        integer, intent(out) :: status
        type(fortnum_status_t) :: root_status
        real(dp), target :: check_flux, root(1)
        integer , target :: failed, ierr
        type(profile_at_s_residual_context_t) :: residual_context
        status = 1
        psi = ieee_value(0.0_dp, ieee_quiet_nan)
        phi = psi; iota = psi; pressure = psi
        if (.not. ieee_is_finite(s)) return
        if (s < 0 .or. s > 1) return
        failed = 0
        if (s == 0 .or. s == 1 .or. o%family == NESTED_INTEGER) then
            psi = s*o%delta
        else
            ! Equation-level derivative: dpsi=(d(s*edge)-dPhi|psi)/Phi'(psi).
            ! Differentiate the equation, not Newton or adaptive schedules.
        residual_context%edge => edge
        residual_context%failed => failed
        residual_context%o => o
        residual_context%root => root
        residual_context%s => s
        residual_context%tolerance => tolerance
            call multiroot_hybrid(profile_at_s_residual_callback, 1, [s*o%delta], &
                root, root_status, &
                                  xtol=tolerance/16, ftol=tolerance/16, &
                                      ctx=residual_context)
            psi = root(1)
            status = 3
            if (root_status%code /= 0 .or. failed /= 0) return
        end if
        call o%flux(psi, tolerance/16, check_flux, ierr)
        status = 3
        if (ierr /= 0) return
        if (abs(check_flux/edge - s) > tolerance) return
        phi = s*edge
        call o%rotation(psi, 1, tolerance, iota, status)
        if (status /= 0) return
        if (o%family == NESTED_INTEGER) then
            call nested_integer_pressure(psi, o%tau, o%delta, o%edge_pressure, pressure)
        else
            call nested_sheared_pressure(psi, o%delta, o%lam, o%edge_pressure, pressure)
        end if
    end subroutine
    function flux_periods_root_evaluate(self, x) result(y)
        class(flux_periods_root_t), intent(in) :: self
            real(dp), intent(in) :: x
            real(dp) :: y, unused
            if (self%k == 1) then
                call nested_sheared_flux_density(self%rho, x, 0.0_dp, self%o%eps, &
                                                 self%o%tau, self%o%saxis, &
                                                     self%o%lam, y, unused)
            else
                call nested_sheared_flux_density(self%rho, 0.0_dp, x, self%o%eps, &
                                                 self%o%tau, self%o%saxis, &
                                                     self%o%lam, unused, y)
            end if
    end function flux_periods_root_evaluate
        subroutine profile_at_s_residual_callback(x, f, jac, ctx)
            real(dp), intent(in) :: x(:)
            real(dp), intent(out) :: f(:), jac(:, :)
            class(*), intent(in), optional  :: ctx
            real(dp) :: label, flux_value, periods(2), errors(2)
            integer :: flux_status
        type(profile_at_s_residual_context_t) :: callback_context
        if (.not. present(ctx)) &
            error stop "missing private callback context"
        select type(selected_context => ctx)
        type is (profile_at_s_residual_context_t)
            callback_context = selected_context
        class default
            error stop "invalid private callback context"
        end select
            ! Monotone linear continuation outside the physical domain is
            ! solely for rejected Newton trials; the admitted root is inside.
            label = max(0.0_dp, min(callback_context%o%delta, x(1)))
            call callback_context%o%flux(label, callback_context%tolerance/16, &
                flux_value, flux_status)
            if (flux_status /= 0) callback_context%failed = flux_status
            call flux_periods(callback_context%o, label, &
                callback_context%tolerance/16, periods, errors, flux_status)
            if (flux_status /= 0) callback_context%failed = flux_status
            call nested_profile_residual(flux_value, callback_context%edge, &
                callback_context%s, x(1), label, &
                                         periods(1), f(1), jac(1, 1))
        end subroutine
end module
