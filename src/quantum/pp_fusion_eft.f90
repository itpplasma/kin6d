! Pionless EFT threshold pp fusion. Natural lengths: fm; energies: MeV.
module kin6d_pp_fusion_eft
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use kin6d_pp_fusion, only: pp_hc, pp_hbar, pp_alpha, pp_mp, pp_mn, &
        pp_phase_space, axial, ft0
    use kin6d_pp_fusion_leaves, only: pp_sfactor, pp_sfactor_jvp
    use kin6d_pp_fusion_eft_leaves, only: eft_scales, eft_scales_jvp, &
        eft_orders, eft_orders_jvp, eft_coulomb_node, eft_coulomb_slope, &
        eft_rates, eft_rates_jvp, eft_square_error, eft_coulomb_moment, &
        eft_energy_coefficients, eft_energy_slope
    use fortnum_special_coulomb, only: coulomb_c0_squared
    use fortnum_quadrature, only: gauss_legendre
    implicit none
    private
    real(dp), parameter :: euler = 0.57721566490153286061_dp
    real(dp), parameter, public :: pp_eft_input_errors(4) = &
        [9e-6_dp, .0026_dp, .014_dp, .002_dp]
    real(dp), parameter, public :: pp_eft_breakdown_mev = 140._dp
    ! NPLQCD arXiv:1610.04545v2 (14), not the superseded v1 uncertainties.
    real(dp), parameter, public :: pp_eft_lattice_errors(4) = [.2_dp, 1._dp, &
        .4_dp, .9_dp]
    type, public :: pp_eft_config_t
        ! AV18 Tables VIII,X: measured binding and fitted Coulomb ERE inputs.
        real(dp) :: binding = 2.224575_dp, a_pp = -7.8063_dp, r_pp = 2.794_dp
        ! Deuteron pole effective range, BC (7); not the zero-energy triplet range.
        real(dp) :: rho_d = 1.764_dp, l1a = 3.9_dp
        ! PDS subtraction convention, not a pion degree of freedom.
        real(dp) :: mu_mev = pp_eft_breakdown_mev
    end type
    type, public :: pp_eft_result_t
        ! Entries: LO, NLO, NNLO, range-resummed NNLO diagnostic.
        real(dp) :: lambda(4) = 0, lambda2(4) = 0, sfactor(4) = 0, slope(4) = 0
        real(dp) :: truncation_lambda(4) = 0, truncation_lambda2(4) = 0
        real(dp) :: truncation_s(4) = 0
        real(dp) :: coefficients(3) = 0, lbar = 0, gamma = 0, chi = 0
        real(dp) :: phase_space = 0, expansion_q = 0, coulomb_integral = 0
        real(dp) :: quadrature_indicator = 0
        character(16) :: status = 'candidate'
        character(48) :: numerical_bound = 'unavailable; quadrature indicator only'
        character(48) :: model_bound = 'unavailable; EFT truncation estimate only'
    end type
    public :: pp_eft_evaluate, pp_eft_jvp, pp_eft_vjp, pp_eft_coulomb_integral
contains
    ! J(chi)=chi*I(chi)=1-chi*exp(chi)*E1(chi), KR (22),(28).
    ! Evaluate the Coulomb integral, using the owning special-function library.
    subroutine pp_eft_coulomb_integral(chi, value, indicator, ok, moment)
        real(dp), intent(in) :: chi
        real(dp), intent(out) :: value, indicator
        logical, intent(out) :: ok
        real(dp), intent(out), optional :: moment
        real(dp) :: estimates(2), moments(2), x(1024), w(1024), eta, factor, fac2, c2
        integer :: level, n, i
        value = 0; indicator = 0; ok = .false.
        if (present(moment)) moment = 0
        if (.not. ieee_is_finite(chi)) return
        if (chi < 0 .or. chi > 1) return
        if (chi <= 0) then
            if (present(moment)) return ! the neutral p^2 moment diverges
            value = 1; ok = .true.
            return
        end if
        ! The endpoint layer for chi<0.01 needs a different numerical schedule.
        if (chi < .01_dp) return
        estimates = 0; moments = 0
        do level = 1, 4
            n = 64*2**level
            call gauss_legendre(n, x(:n), w(:n))
            estimates(1) = estimates(2)
            estimates(2) = 0
            moments(1) = moments(2); moments(2) = 0
            !$omp parallel do reduction(+:estimates,moments) private(eta,factor,fac2,c2)
            do i = 1, n
                call eft_coulomb_node(x(i), chi, eta, factor)
                c2 = coulomb_c0_squared(eta)
                estimates(2) = estimates(2) + w(i)*factor*c2
                if (.not. present(moment)) cycle
                call eft_coulomb_moment(x(i), chi, fac2)
                moments(2) = moments(2) + w(i)*fac2*c2
            end do
            !$omp end parallel do
            if (level < 2) cycle
            if (abs(estimates(2) - estimates(1)) >= 1e-13_dp) cycle
            if (present(moment)) then
                if (abs(moments(2) - moments(1)) >= 1e-11_dp*max(1._dp, &
                    abs(moments(2)))) cycle
            end if
            exit
        end do
        value = estimates(2); indicator = abs(estimates(2) - estimates(1))
        ok = ieee_is_finite(value) .and. indicator < 1e-13_dp
        if (present(moment)) then
            moment = moments(2)
            ok = ok .and. abs(moments(2) - moments(1)) < 1e-11_dp*max(1._dp, &
                abs(moment))
        end if
    end subroutine

    ! Jacobian rows: lambda2(LO,NLO,NNLO,resummed), then S in the same order.
    ! Columns: B_d, a_pp, r_pp, rho_d, L1A. Weak constants and PDS scale frozen.
    subroutine pp_eft_evaluate(cfg, answer, ok, jacobian)
        type(pp_eft_config_t), intent(in) :: cfg
        type(pp_eft_result_t), intent(out) :: answer
        logical, intent(out) :: ok
        real(dp), intent(out), optional :: jacobian(8, 5)
        real(dp) :: gamma, chi, screen, mu, jc, djc, df, ceiling, zeros(5), v(5)
        real(dp) :: dg, dx, dscr, dc(3), dres, dlb, drates(4), unused, moment
        real(dp) :: ep(3), epr, energy_derivative(4)
        integer :: i, j, k
        answer = pp_eft_result_t(); ok = .false.
        if (present(jacobian)) jacobian = 0
        if (.not. ieee_is_finite(cfg%binding)) return
        if (.not. ieee_is_finite(cfg%a_pp)) return
        if (.not. ieee_is_finite(cfg%r_pp)) return
        if (.not. ieee_is_finite(cfg%rho_d)) return
        if (.not. ieee_is_finite(cfg%l1a)) return
        if (.not. ieee_is_finite(cfg%mu_mev)) return
        if (cfg%binding <= 0 .or. cfg%a_pp >= 0 .or. cfg%r_pp < 0 .or. &
            cfg%rho_d < 0 .or. cfg%mu_mev <= 0) return
        mu = cfg%mu_mev/pp_hc
        call eft_scales(cfg%binding, pp_mp, pp_mn, pp_hc, pp_alpha, mu, euler, &
            gamma, chi, screen)
        if (gamma*cfg%rho_d >= 1 .or. gamma >= mu) return
        if (abs(1/cfg%a_pp - mu + screen) < 1e-8_dp) return
        call pp_eft_coulomb_integral(chi, jc, answer%quadrature_indicator, ok, moment)
        if (.not. ok) return
        ok = .false.
        call eft_orders(gamma, cfg%a_pp, cfg%r_pp, cfg%rho_d, cfg%l1a, axial, &
            mu, chi, jc, screen, answer%coefficients(1), answer%coefficients(2), &
            answer%coefficients(3), answer%lambda(4), answer%lbar)
        do i = 1, 3
            answer%lambda(i) = sum(answer%coefficients(:i))
        end do
        if (.not. all(ieee_is_finite(answer%lambda))) return
        if (any(answer%lambda <= 0)) return
        call eft_rates(answer%coefficients(1), answer%coefficients(2), &
            answer%coefficients(3), answer%lambda(4), answer%lambda2(1), &
            answer%lambda2(2), answer%lambda2(3), answer%lambda2(4))
        if (any(answer%lambda2 <= 0)) return
        answer%gamma = gamma; answer%chi = chi; answer%coulomb_integral = jc
        answer%phase_space = pp_phase_space(0._dp, cfg%binding, df)
        if (answer%phase_space <= 0) return
        answer%expansion_q = gamma/(pp_eft_breakdown_mev/pp_hc)
        ceiling = abs(answer%coefficients(1))
        do i = 1, 4
            k = min(i, 3)
            ceiling = max(ceiling, &
                abs(answer%coefficients(k))/answer%expansion_q**(k - 1))
            answer%truncation_lambda(i) = ceiling*answer%expansion_q**k
            call eft_square_error(answer%lambda(i), answer%truncation_lambda(i), &
                answer%truncation_lambda2(i))
            call pp_sfactor(pp_mp, pp_alpha, pp_hbar, pp_hc, axial, gamma, &
                answer%lambda2(i), answer%phase_space, ft0, answer%sfactor(i))
            call pp_sfactor(pp_mp, pp_alpha, pp_hbar, pp_hc, axial, gamma, &
                answer%truncation_lambda2(i), answer%phase_space, ft0, &
                answer%truncation_s(i))
        end do
        call eft_energy_coefficients(gamma, cfg%a_pp, cfg%r_pp, cfg%rho_d, &
            cfg%l1a, axial, mu, chi, jc, screen, moment/gamma**2, pp_alpha, &
            pp_mp/pp_hc, ep(1), ep(2), ep(3), epr)
        call eft_rates_jvp(answer%coefficients(1), answer%coefficients(2), &
            answer%coefficients(3), answer%lambda(4), ep(1), ep(2), ep(3), epr, &
            energy_derivative(1), energy_derivative(2), energy_derivative(3), &
            energy_derivative(4))
        do i = 1, 4
            call eft_energy_slope(energy_derivative(i), answer%lambda2(i), &
                answer%phase_space, df, pp_mp/pp_hc, pp_hc, answer%slope(i))
        end do
        if (.not. all(ieee_is_finite(answer%slope))) return
        ok = .true.
        if (.not. present(jacobian)) return
        call eft_coulomb_slope(chi, jc, djc)
        zeros = 0
        do j = 1, 5
            v = zeros; v(j) = 1
            call eft_scales_jvp(cfg%binding, pp_mp, pp_mn, pp_hc, pp_alpha, &
                mu, euler, v(1), 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, &
                dg, dx, dscr)
            call eft_orders_jvp(gamma, cfg%a_pp, cfg%r_pp, cfg%rho_d, cfg%l1a, &
                axial, mu, chi, jc, screen, dg, v(2), v(3), v(4), v(5), &
                0._dp, 0._dp, dx, djc*dx, dscr, dc(1), dc(2), dc(3), dres, dlb)
            call eft_rates_jvp(answer%coefficients(1), answer%coefficients(2), &
                answer%coefficients(3), answer%lambda(4), dc(1), dc(2), dc(3), &
                dres, drates(1), drates(2), drates(3), drates(4))
            do i = 1, 4
                jacobian(i, j) = drates(i)
                call pp_sfactor_jvp(pp_mp, pp_alpha, pp_hbar, pp_hc, axial, &
                    gamma, answer%lambda2(i), answer%phase_space, ft0, &
                    0._dp, 0._dp, 0._dp, 0._dp, 0._dp, dg, drates(i), df*v(1), &
                    0._dp, unused)
                jacobian(i + 4, j) = unused
            end do
        end do
    end subroutine

    subroutine pp_eft_jvp(cfg, direction, answer, product, ok)
        type(pp_eft_config_t), intent(in) :: cfg
        real(dp), intent(in) :: direction(5)
        type(pp_eft_result_t), intent(out) :: answer
        real(dp), intent(out) :: product(8)
        logical, intent(out) :: ok
        real(dp) :: jacobian(8, 5)
        call pp_eft_evaluate(cfg, answer, ok, jacobian)
        product = matmul(jacobian, direction)
    end subroutine

    subroutine pp_eft_vjp(cfg, weights, answer, product, ok)
        type(pp_eft_config_t), intent(in) :: cfg
        real(dp), intent(in) :: weights(8)
        type(pp_eft_result_t), intent(out) :: answer
        real(dp), intent(out) :: product(5)
        logical, intent(out) :: ok
        real(dp) :: jacobian(8, 5)
        call pp_eft_evaluate(cfg, answer, ok, jacobian)
        product = matmul(transpose(jacobian), weights)
    end subroutine
end module
