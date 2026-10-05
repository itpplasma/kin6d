program test_fermi_gas
    use, intrinsic :: iso_fortran_env, only: dp => real64, qp => real128
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_is_finite
    use fortnum_interval, only: interval_t, interval
    use kin6d_fermi_gas, only: fermi_workspace_t, fermi_state_t, fermi_si_state_t, &
        initialize_fermi_workspace, solve_fermi_gas, solve_fermi_gas_si, &
        occupation_fermi_gas, fermi_gas_jvp, fermi_gas_vjp
    use kin6d_fermi_gas_leaves, only: fermi_positive, fermi_negative, &
        fermi_positive_jvp, fermi_positive_vjp, fermi_si_parameters
    implicit none
    type(fermi_workspace_t) :: work
    type(fermi_state_t) :: state, corrupted
    type(fermi_si_state_t) :: si
    real(dp), parameter :: temperatures(3) = [0.05_dp, 1.0_dp, 10.0_dp]
    real(dp) :: value, lower, upper, tangent(6), reverse, cotangent(6)
    real(dp) :: pf, ef, theta, nan, derivative, adjoint
    real(qp) :: reference(6), gradient(6), refined(6), refined_gradient(6)
    real(qp) :: exact_pf, exact_ef, exact_theta
    type(interval_t) :: bands(6), tangent_band(6), reverse_band
    integer :: case_index, j
    logical :: ok

    call initialize_fermi_workspace(work, 16384, ok)
    call require(ok, 'workspace')
    call fermi_positive(0.0_dp, value)
    call require(value <= 0.5_dp .and. value >= 0.5_dp, 'exact half occupation')
    call fermi_positive_jvp(0.0_dp, 2.0_dp, derivative)
    call fermi_positive_vjp(0.0_dp, 2.0_dp, adjoint)
    call require(derivative <= -0.5_dp .and. derivative >= -0.5_dp .and. &
                 adjoint <= -0.5_dp .and. adjoint >= -0.5_dp, 'chart derivative')
    call fermi_positive(800.0_dp, value)
    call require(value <= 0.0_dp .and. value >= 0.0_dp, 'positive stable chart')
    call fermi_negative(-800.0_dp, value)
    call require(value <= 1.0_dp .and. value >= 1.0_dp, 'negative stable chart')
    cotangent = [0.25_dp, -0.5_dp, 0.125_dp, 0.75_dp, -0.25_dp, 0.5_dp]
    do case_index = 1, size(temperatures)
        call initialize_fermi_workspace(work, 16384, ok)
        call solve_fermi_gas(work, temperatures(case_index), state, ok)
        call require(ok .and. state%certified, 'general EOS certificate')
        call oracle(real(state%theta, qp), 2048, reference, gradient)
        call oracle(real(state%theta, qp), 4096, refined, refined_gradient)
        call require(maxval(abs(reference - refined)) < 1.0e-10_qp, &
                     'independent oracle refinement')
        call require(maxval(abs(gradient - refined_gradient)) < 1.0e-9_qp, &
                     'independent derivative oracle refinement')
        bands = [state%mu_bound, state%pressure_bound, state%energy_bound, &
            state%susceptibility_bound, state%screening_ratio_bound, &
            state%bulk_modulus_bound]
        do j = 1, 6
            call require(encloses(bands(j), refined(j)), 'EOS independent containment')
            call require(encloses(state%theta_gradient_bound(j), refined_gradient(j)), &
                         'implicit derivative containment')
        end do
        call require(encloses(state%density_bound, 1.0_qp), 'fixed density')
        call require(abs(state%energy - 1.5_dp*state%pressure) < 2.0e-14_dp, &
                     'ideal gas energy pressure identity')
        if (case_index == 1) then
            call require(abs(state%mu - 1.0_dp) < 0.003_dp .and. &
                abs(state%pressure - 0.4_dp) < 0.005_dp .and. &
                abs(state%susceptibility - 1.5_dp) < 0.004_dp, 'cold limit approach')
        end if
        if (case_index == 3) then
            call require(abs(state%pressure/state%theta - 1.0_dp) < 0.005_dp .and. &
                abs(state%susceptibility*state%theta - 1.0_dp) < 0.009_dp, &
                'classical limit approach')
        end if
        call fermi_gas_jvp(state, 0.5_dp, tangent, tangent_band, ok)
        call require(ok, 'equation JVP')
        call fermi_gas_vjp(state, cotangent, reverse, reverse_band, ok)
        call require(ok, 'equation VJP')
        call require(abs(dot_product(cotangent, tangent) - 0.5_dp*reverse) &
                     < 2.0e-12_dp, 'Euclidean output pairing')
        call require(encloses(reverse_band, &
            dot_product(real(cotangent, qp), refined_gradient)), 'independent VJP')
        call occupation_fermi_gas(state, 0.875_dp, value, lower, upper, ok)
        call require(ok, 'occupation query')
        call require(encloses(interval(lower, upper), &
            1.0_qp/(1.0_qp + exp((0.875_qp**2 - refined(1))/ &
                                real(state%theta, qp)))), 'occupation containment')
        corrupted = state
        corrupted%mu = corrupted%mu + 0.01_dp
        call occupation_fermi_gas(corrupted, 0.5_dp, value, lower, upper, ok)
        call require(.not. ok, 'mutated chemical potential rejected')
        corrupted = state
        corrupted%theta_bound = interval(0.25_dp)
        call fermi_gas_jvp(corrupted, 1.0_dp, tangent, tangent_band, ok)
        call require(.not. ok, 'mutated input packet rejected')
    end do

    call solve_fermi_gas(work, 1.0_dp, state, ok, certify=.false.)
    call require(ok .and. .not. state%certified, 'candidate alone')
    call occupation_fermi_gas(state, 0.5_dp, value, lower, upper, ok)
    call require(.not. ok, 'unprepared candidate rejected')
    corrupted = fermi_state_t()
    corrupted%certified = .true.; corrupted%derivative_certified = .true.
    call occupation_fermi_gas(corrupted, 0.0_dp, value, lower, upper, ok)
    call require(.not. ok, 'forged certificate rejected')
    call solve_fermi_gas(work, 0.0_dp, state, ok)
    call require(.not. ok, 'zero temperature explicit unsupported')
    nan = ieee_value(0.0_dp, ieee_quiet_nan)
    call solve_fermi_gas(work, nan, state, ok)
    call require(.not. ok, 'nonfinite input')
    call solve_fermi_gas(work, 1.0_dp, state, ok, &
                         theta_interval=interval(0.0_dp, 1.0_dp))
    call require(.not. ok, 'input family crossing zero')
    call initialize_fermi_workspace(work, 32, ok)
    call require(.not. ok, 'invalid workspace')
    call initialize_fermi_workspace(work, 16384, ok)
    ! Supplied exact binary64 constants define this dimensionless SI-map fixture.
    call fermi_si_parameters(1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, pf, ef, theta)
    call solve_fermi_gas_si(work, 1.0_dp, ef, 1.0_dp, 1.0_dp, 1.0_dp, &
                            2.0_dp, 1.0_dp, si, ok)
    call require(ok .and. si%certified, 'SI normalized closure')
    call require(encloses(si%fermi_momentum_bound, &
        (3.0_qp*acos(-1.0_qp)**2)**(1.0_qp/3.0_qp)), 'SI momentum oracle')
    exact_pf = (3.0_qp*acos(-1.0_qp)**2)**(1.0_qp/3.0_qp)
    exact_ef = exact_pf**2/2.0_qp
    exact_theta = real(ef, qp)/exact_ef
    call require(encloses(si%normalized%theta_bound, exact_theta), &
                 'SI temperature rounding ledger')
    call oracle(exact_theta, 4096, refined, refined_gradient)
    call require(encloses(si%chemical_potential_bound, exact_ef*refined(1)), &
                 'SI root normalization containment')
    call require(encloses(si%pressure_bound, exact_ef*refined(2)), &
                 'SI physical pressure oracle')
    call require(encloses(si%energy_density_bound, exact_ef*refined(3)), &
                 'SI physical energy oracle')
    call require(encloses(si%compressibility_bound, refined(4)/exact_ef), &
                 'SI physical compressibility oracle')
    call require(encloses(si%screening_wavenumber_bound, &
        sqrt(4.0_qp*refined(4)/exact_ef)), 'SI physical screening oracle')
    call require(encloses(si%bulk_modulus_bound, exact_ef/refined(4)), &
                 'SI physical bulk modulus oracle')
    call require(encloses(si%pressure_bound, &
        real(si%fermi_energy, qp)*real(si%normalized%pressure, qp)), 'SI pressure map')
    call require(abs(si%energy_density - 1.5_dp*si%pressure) < 1.0e-12_dp, &
                 'physical ideal energy identity')
    print '(a)', 'PASS: Fermi EOS, strict integrals/root, implicit AD, SI input ledger'

contains

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(*), intent(in) :: message
        if (.not. condition) then
            print '(a)', 'FAIL: '//message
            error stop 1
        end if
    end subroutine require

    logical function encloses(bound, expected) result(ok)
        type(interval_t), intent(in) :: bound
        real(qp), intent(in) :: expected
        ok = .false.
        if (.not. all(ieee_is_finite([bound%lo, bound%hi]))) return
        if (bound%lo > bound%hi) return
        ok = real(bound%lo, qp) <= expected .and. expected <= real(bound%hi, qp)
    end function encloses

    subroutine oracle(theta, panels, values, gradients)
        real(qp), intent(in) :: theta
        integer, intent(in) :: panels
        real(qp), intent(out) :: values(6), gradients(6)
        real(qp) :: lower, upper, mu, moments(8), dmu, dchi, radius
        integer :: iteration
        lower = -theta*(4.0_qp + 2.0_qp*log(1.0_qp + theta))
        upper = 4.0_qp + 4.0_qp*theta
        radius = sqrt(4.0_qp + 36.0_qp*theta)
        do iteration = 1, 70
            mu = (lower + upper)/2.0_qp
            call simpson_moments(theta, mu, radius, panels, moments)
            if (3.0_qp*moments(1) < 1.0_qp) then
                lower = mu
            else
                upper = mu
            end if
        end do
        mu = (lower + upper)/2.0_qp
        call simpson_moments(theta, mu, radius, panels, moments)
        values = [mu, 2.0_qp*moments(2), 3.0_qp*moments(2), &
            3.0_qp*moments(3), sqrt(2.0_qp*moments(3)), &
            1.0_qp/(3.0_qp*moments(3))]
        dmu = -moments(4)/moments(3)
        dchi = 3.0_qp*(moments(8) + moments(7)*dmu)
        gradients = [dmu, 2.0_qp*(moments(6) + moments(5)*dmu), &
            3.0_qp*(moments(6) + moments(5)*dmu), dchi, &
            dchi/(3.0_qp*values(5)), -dchi/values(4)**2]
    end subroutine oracle

    subroutine simpson_moments(theta, mu, radius, panels, moments)
        real(qp), intent(in) :: theta, mu, radius
        integer, intent(in) :: panels
        real(qp), intent(out) :: moments(8)
        real(qp) :: r, h, z, f, w, logistic, excess, integrand(8)
        integer :: node
        h = radius/real(panels, qp)
        moments = 0.0_qp
        do node = 0, panels
            r = real(node, qp)*h
            z = (r*r - mu)/theta
            ! Test-only direct high-precision formula; chart split is unnecessary
            ! for these bounded oracle fixtures, independent of production DAGs.
            f = 1.0_qp/(1.0_qp + exp(z))
            logistic = f*(1.0_qp - f)
            excess = r*r - mu
            integrand = [r*r*f, r**4*f, r*r*logistic/theta, &
                r*r*excess*logistic/theta**2, r**4*logistic/theta, &
                r**4*excess*logistic/theta**2, &
                r*r*logistic*(1.0_qp - 2.0_qp*f)/theta**2, &
                r*r*logistic*(excess*(1.0_qp - 2.0_qp*f)/theta**3 &
                             - 1.0_qp/theta**2)]
            w = 2.0_qp
            if (mod(node, 2) == 1) w = 4.0_qp
            if (node == 0 .or. node == panels) w = 1.0_qp
            moments = moments + w*integrand
        end do
        moments = moments*h/3.0_qp
    end subroutine simpson_moments
end program test_fermi_gas
