module rm_mhd_oracles
    ! Test references only; no dependence on the evolving MHD residual.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use kin6d_rm_oracle_leaves, only: rm_wave_roots, rm_fkr, rm_internal_kink, &
        rm_rutherford, rm_symmetric_saturation, rm_cylinder_saturation, &
        rm_wesson_radial, rm_wesson_series, rm_wesson_curvature, &
        rm_outer_rhs, rm_frobenius_fit, &
        rm_newcomb_coefficients, rm_reduced_coefficients
    use fortnum_ode, only: ode_problem_t, ode_workspace_t, ode_solution_t, &
        ode_integrate
    use fortnum_roots, only: root_brent
    use fortnum_status, only: fortnum_status_t
    implicit none
    private
    real(dp), parameter, public :: oracle_mu0 = 4e-7_dp*acos(-1.0_dp)
    integer, parameter, public :: oracle_ok = 0, oracle_domain = -1, &
        oracle_no_resonance = -2, oracle_singular = -3
    public :: uniform_waves, tearing_rates, rutherford_speed, symmetric_width
    public :: shoot_outer, wesson_t, wesson_reference, newcomb_coefficients
    public :: cylinder_saturated_width, outer_result_t, outer_coefficients
    type :: outer_result_t
        real(dp) :: rs = 0, delta = 0, slope_minus = 0, slope_plus = 0
        real(dp) :: flux_minus = 0, flux_plus = 0
        integer :: status = oracle_domain
    end type
    type :: wesson_t
        real(dp) :: a = 1, q0 = 0.9_dp, qa = 3
        integer :: m = 2, n = 1
    contains
        procedure :: evaluate => wesson_evaluate
    end type
    abstract interface
        subroutine outer_coefficients(r, p, u)
            import :: dp
            real(dp), intent(in) :: r
            real(dp), intent(out) :: p, u
        end subroutine
    end interface
contains
    subroutine uniform_waves(k, b, rho, pressure, adiabatic, omega2, status)
        real(dp), intent(in) :: k(3), b(3), rho, pressure, adiabatic
        ! Squared angular frequencies: Alfven, fast, slow, in s^-2.
        real(dp), intent(out) :: omega2(3)
        integer, intent(out) :: status
        real(dp) :: k2, b2, kb, cs2, cross(3)
        omega2 = 0; status = oracle_domain
        if (.not. all(ieee_is_finite([k, b, rho, pressure, adiabatic]))) return
        if (rho <= 0 .or. pressure < 0 .or. adiabatic <= 0) return
        k2 = dot_product(k, k); b2 = dot_product(b, b); kb = dot_product(k, b)
        cs2 = adiabatic*pressure/rho
        status = oracle_ok
        if (k2 <= 0) return
        if (b2 <= 0 .and. cs2 <= 0) return
        cross = [k(2)*b(3) - k(3)*b(2), k(3)*b(1) - k(1)*b(3), &
            k(1)*b(2) - k(2)*b(1)]
        call rm_wave_roots(k2, kb, b2, dot_product(cross, cross), &
            rho, oracle_mu0, cs2, &
            omega2(1), omega2(2), omega2(3))
    end subroutine

    subroutine tearing_rates(eta, rho, fprime, delta, fkr, kink, status)
        ! eta: Ohm m; rho: kg/m^3; F'=d(k.B)/dr: T/m^2; Delta': 1/m.
        real(dp), intent(in) :: eta, rho, fprime, delta
        real(dp), intent(out) :: fkr, kink
        integer, intent(out) :: status
        real(dp) :: constant
        status = oracle_domain; fkr = 0; kink = 0
        if (.not. all(ieee_is_finite([eta, rho, fprime, delta]))) return
        if (eta < 0 .or. rho <= 0 .or. abs(fprime) <= 0 .or. delta < 0) return
        constant = 2*acos(-1.0_dp)*gamma(0.75_dp)/gamma(0.25_dp)
        call rm_fkr(delta, constant, eta/oracle_mu0, fprime, oracle_mu0, rho, fkr)
        call rm_internal_kink(eta/oracle_mu0, fprime, oracle_mu0, rho, kink)
        status = oracle_ok
    end subroutine

    pure real(dp) function rutherford_speed(eta, delta) result(speed)
        real(dp), intent(in) :: eta, delta
        call rm_rutherford(eta/oracle_mu0, delta, speed)
    end function

    subroutine symmetric_width(delta, curvature, width, status)
        real(dp), intent(in) :: delta, curvature
        real(dp), intent(out) :: width
        integer, intent(out) :: status
        status = oracle_domain; width = 0
        if (.not. all(ieee_is_finite([delta, curvature]))) return
        if (delta <= 0 .or. curvature >= 0) return
        call rm_symmetric_saturation(delta, curvature, width)
        status = oracle_ok
    end subroutine

    subroutine cylinder_saturated_width(delta, a, alpha, curvature, asymmetry, &
            rs, bracket, width, status)
        real(dp), intent(in) :: delta, a, alpha, curvature, asymmetry, rs
        real(dp), intent(in) :: bracket(2)
        real(dp), intent(out) :: width
        integer, intent(out) :: status
        type(fortnum_status_t) :: stat
        status = oracle_domain; width = 0
        if (.not. all(ieee_is_finite([delta, a, alpha, curvature, asymmetry, &
            rs, bracket]))) return
        if (a <= 0 .or. rs <= 0 .or. delta <= 0) return
        if (bracket(1) <= 0 .or. bracket(2) <= bracket(1)) return
        call root_brent(residual, bracket(1), bracket(2), width, stat, &
            xtol=1e-12_dp*a)
        status = stat%code
    contains
        pure real(dp) function residual(w) result(value)
            real(dp), intent(in) :: w
            call rm_cylinder_saturation(w, delta, a, alpha, curvature, &
                asymmetry, rs, value)
        end function
    end subroutine

    pure subroutine wesson_evaluate(model, r, f, fp, fpp, jp, jpp, jstar)
        class(wesson_t), intent(in) :: model
        real(dp), intent(in) :: r
        real(dp), intent(out) :: f, fp, fpp, jp, jstar
        real(dp), intent(out), optional :: jpp
        real(dp) :: nu
        nu = model%qa/model%q0 - 1
        if (r/model%a < 0.02_dp) then
            call rm_wesson_series(r, model%a, model%q0, nu, real(model%m, dp), &
                real(model%n, dp), f, fp, fpp, jp, jstar)
        else
            call rm_wesson_radial(r, model%a, model%q0, nu, real(model%m, dp), &
                real(model%n, dp), f, fp, fpp, jp, jstar)
        end if
        ! Noninteger nu can have divergent curvature at the wall; ODEs need J' only.
        if (present(jpp)) call rm_wesson_curvature(r, model%a, model%q0, nu, jpp)
    end subroutine

    subroutine wesson_reference(model, cut, tolerance, result)
        type(wesson_t), intent(in) :: model
        real(dp), intent(in) :: cut, tolerance
        type(outer_result_t), intent(out) :: result
        real(dp) :: rs, f, fp, fpp, jp, jpp, jstar
        type(fortnum_status_t) :: stat
        result%status = oracle_domain
        if (.not. all(ieee_is_finite([model%a, model%q0, model%qa]))) return
        if (model%a <= 0 .or. model%q0 <= 0 .or. model%qa <= 2*model%q0) return
        if (model%m < 1 .or. model%n < 1) return
        result%status = oracle_no_resonance
        if (real(model%m, dp)/model%n <= model%q0) return
        if (real(model%m, dp)/model%n >= model%qa) return
        call root_brent(resonance, 0.0_dp, model%a, rs, stat, &
            xtol=1e-13_dp*model%a)
        if (stat%code /= 0) then
            result%status = stat%code
            return
        end if
        call model%evaluate(rs, f, fp, fpp, jp, jpp, jstar)
        call shoot_outer(coefficients, model%m, model%a, rs, &
            model%m*jp/(rs*fp), cut, tolerance, result)
    contains
        pure real(dp) function resonance(r) result(value)
            real(dp), intent(in) :: r
            real(dp) :: v1, v2, v3, v5
            call model%evaluate(r, value, v1, v2, v3, jstar=v5)
        end function
        subroutine coefficients(r, p, u)
            real(dp), intent(in) :: r
            real(dp), intent(out) :: p, u
            real(dp) :: ff, v1, v2, jj, v4
            call model%evaluate(r, ff, v1, v2, jj, jstar=v4)
            call rm_reduced_coefficients(r, real(model%m, dp), ff, jj, p, u)
        end subroutine
    end subroutine

    pure subroutine newcomb_coefficients(r, m, kz, bt, bz, f, fp, fpp, p, u, pole)
        real(dp), intent(in) :: r, m, kz, bt, bz, f, fp, fpp
        real(dp), intent(out) :: p, u, pole
        real(dp) :: numerator
        ! U=regular+numerator/F. At F=0 only its pole residue is requested.
        call rm_newcomb_coefficients(r, m, kz, bt, bz, fp, fpp, p, u, numerator)
        pole = huge(1.0_dp)
        if (abs(fp) > 0) pole = numerator/fp
        if (abs(f) > 0) u = u + numerator/f
    end subroutine

    subroutine shoot_outer(coefficients, m, wall, rs, pole, cut, tolerance, result, &
            axis_start, wall_logderivative, axis_logderivative)
        procedure(outer_coefficients) :: coefficients
        integer, intent(in) :: m
        real(dp), intent(in) :: wall, rs, pole, cut, tolerance
        type(outer_result_t), intent(out) :: result
        real(dp), intent(in), optional :: axis_start, wall_logderivative
        real(dp), intent(in), optional :: axis_logderivative
        real(dp) :: rmin, initial(2), left(2), right(2), h, pm, pp, um, up, p0, u0
        result%status = oracle_domain; result%rs = rs
        if (.not. all(ieee_is_finite([wall, rs, pole, cut, tolerance]))) return
        if (m < 1 .or. wall <= 0 .or. tolerance <= 0 .or. cut <= 0) return
        rmin = 1e-5_dp*wall
        if (present(axis_start)) rmin = axis_start
        if (rmin <= 0 .or. rmin >= rs - cut .or. rs + cut >= wall) return
        h = min(1e-4_dp*wall, (rs - rmin)/4, (wall - rs)/4)
        call coefficients(rs - h, pm, um)
        call coefficients(rs + h, pp, up)
        p0 = (pm + pp)/2; u0 = (um + up)/2
        initial = [1.0_dp, real(m, dp)/rmin]
        if (present(axis_logderivative)) initial(2) = axis_logderivative
        call integrate(rmin, rs - cut, initial, left)
        if (result%status /= oracle_ok) return
        initial = [0.0_dp, 1.0_dp]
        if (present(wall_logderivative)) initial = [1.0_dp, wall_logderivative]
        call integrate(wall, rs + cut, initial, right)
        if (result%status /= oracle_ok) return
        result%status = oracle_singular
        ! A vanishing large-solution amplitude is a kink, not a huge Delta'.
        if (abs(left(1) + cut*left(2)) < &
            1e-7_dp*max(abs(left(1)), wall*abs(left(2)))) return
        if (abs(right(1) - cut*right(2)) < &
            1e-7_dp*max(abs(right(1)), wall*abs(right(2)))) return
        call rm_frobenius_fit(left(1), left(2), -cut, wall, pole, p0, u0, &
            result%flux_minus, result%slope_minus)
        call rm_frobenius_fit(right(1), right(2), cut, wall, pole, p0, u0, &
            result%flux_plus, result%slope_plus)
        result%delta = result%slope_plus - result%slope_minus
        if (.not. ieee_is_finite(result%delta)) return
        result%status = oracle_ok
    contains
        subroutine integrate(a, b, initial, final)
            real(dp), intent(in) :: a, b, initial(2)
            real(dp), intent(out) :: final(2)
            type(ode_problem_t) :: problem
            type(ode_workspace_t) :: work
            type(ode_solution_t) :: sol
            type(fortnum_status_t) :: stat
            problem%rhs => rhs
            problem%t0 = a; problem%t1 = b; problem%y0 = initial
            problem%rtol = tolerance; problem%atol = tolerance*1e-2_dp
            problem%hmax = 0.02_dp*wall
            call ode_integrate(problem, work, sol, stat)
            result%status = stat%code
            final = initial
            if (stat%code == 0) final = sol%y(:, size(sol%t))
        end subroutine
        subroutine rhs(r, y, dy, ctx)
            real(dp), intent(in) :: r, y(:)
            real(dp), intent(out) :: dy(:)
            class(*), intent(in), optional :: ctx
            real(dp) :: p, u
            call coefficients(r, p, u)
            dy(1) = y(2)
            call rm_outer_rhs(p, u, y(1), y(2), dy(2))
        end subroutine
    end subroutine
end module rm_mhd_oracles
