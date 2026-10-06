module kin6d_fermi_gas
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortnum_interval, only: interval_t, interval, hull, intersect, mag, is_empty
    use fortnum_quadrature, only: gauss_legendre_ab
    use fortnum_roots, only: root_newton
    use fortnum_status, only: fortnum_status_t, status_ok
    use fortnum_verified_quadrature, only: composite_midpoint_batch
    use kin6d_fermi_gas_leaves, only: fermi_add_interval, fermi_apply_direction, &
        fermi_apply_direction_interval, fermi_apply_direction_vjp, &
        fermi_apply_direction_vjp_interval, fermi_argument, fermi_argument_interval, &
        fermi_curvature_interval, fermi_density_curvature_interval, &
        fermi_density_integrand, fermi_density_integrand_interval, &
        fermi_difference_interval, fermi_eos, fermi_eos_interval, fermi_grid, &
        fermi_half, fermi_integrands, fermi_integrands_interval, fermi_negative, &
        fermi_negative_interval, fermi_positive, fermi_positive_interval, &
        fermi_residual, fermi_residual_interval, fermi_root_band, &
        fermi_root_band_interval, &
        fermi_root_bracket, fermi_scaled_tolerance, fermi_si_eos, &
        fermi_si_eos_interval, fermi_si_parameters, fermi_si_parameters_interval, &
        fermi_tail_interval, fermi_theta_derivative, &
        fermi_theta_derivative_interval, fermi_weighted_add, &
        fermi_relativity_indicator, fermi_relativity_indicator_interval
    use kin6d_relativistic_fermi_leaves, only: fermi_rel_root, fermi_rel_root_interval, &
        fermi_rel_argument, fermi_rel_argument_interval, fermi_rel_integrands, &
        fermi_rel_integrands_interval, fermi_rel_density_integrand, &
        fermi_rel_density_integrand_interval, fermi_rel_density_curvature_interval, &
        fermi_rel_curvature_interval, fermi_rel_eos, fermi_rel_eos_interval, &
        fermi_rel_theta_derivative, fermi_rel_theta_derivative_interval, &
        fermi_rel_root_bracket, fermi_rel_tail_interval, fermi_rel_window, &
        fermi_rel_grid
    implicit none
    private
    public :: fermi_workspace_t, fermi_state_t, initialize_fermi_workspace
    public :: solve_fermi_gas, occupation_fermi_gas, fermi_gas_jvp, fermi_gas_vjp
    public :: fermi_si_state_t, solve_fermi_gas_si

    type :: fermi_workspace_t
        private
        real(dp), allocatable :: edges(:), nodes(:), weights(:)
        integer :: candidate_nodes = 128, maximum_nodes = 4096, cells = 64
        real(dp) :: cached_theta = -1.0_dp, cached_radius = -1.0_dp
        real(dp) :: cached_relativity = 0.0_dp, nodes_relativity = 0.0_dp
        integer :: cached_nodes = 0
        logical :: initialized = .false.
    end type fermi_workspace_t

    type :: fermi_state_t
        real(dp) :: theta = 0.0_dp, mu = 0.0_dp, pressure = 0.0_dp
        ! Relativity parameter x=pF/(m*c); zero selects the nonrelativistic model.
        real(dp) :: relativity = 0.0_dp
        type(interval_t) :: relativity_bound
        real(dp) :: energy = 0.0_dp, susceptibility = 0.0_dp
        real(dp) :: screening_ratio = 0.0_dp, bulk_modulus = 0.0_dp
        type(interval_t) :: theta_bound, mu_bound, pressure_bound, energy_bound
        type(interval_t) :: susceptibility_bound, screening_ratio_bound
        type(interval_t) :: bulk_modulus_bound, density_bound
        real(dp) :: theta_gradient(6) = 0.0_dp
        type(interval_t) :: theta_gradient_bound(6)
        integer :: cells = 0, candidate_nodes = 0
        logical :: valid = .false., certified = .false.
        logical :: derivative_certified = .false.
        character(48) :: status = 'uninitialized'
        logical, private :: prepared = .false.
        real(dp), private :: prepared_theta = 0.0_dp, prepared_center(6) = 0.0_dp
        real(dp), private :: prepared_gradient(6) = 0.0_dp
        type(interval_t), private :: prepared_value(6), prepared_derivative(6)
        type(interval_t), private :: prepared_theta_bound, prepared_relativity
    end type fermi_state_t

    type :: fermi_si_state_t
        type(fermi_state_t) :: normalized
        real(dp) :: fermi_momentum = 0.0_dp, fermi_energy = 0.0_dp
        real(dp) :: chemical_potential = 0.0_dp, pressure = 0.0_dp
        real(dp) :: energy_density = 0.0_dp, compressibility = 0.0_dp
        real(dp) :: screening_wavenumber = 0.0_dp, bulk_modulus = 0.0_dp
        type(interval_t) :: fermi_momentum_bound, fermi_energy_bound
        type(interval_t) :: chemical_potential_bound, pressure_bound
        type(interval_t) :: energy_density_bound, compressibility_bound
        type(interval_t) :: screening_wavenumber_bound, bulk_modulus_bound
        logical :: certified = .false.
    end type fermi_si_state_t

contains

    subroutine solve_fermi_gas_si(work, density, temperature, mass, hbar, kb, &
                                 charge, epsilon, state, ok, speed)
        type(fermi_workspace_t), intent(inout) :: work
        real(dp), intent(in) :: density, temperature, mass, hbar, kb, charge, epsilon
        type(fermi_si_state_t), intent(out) :: state
        logical, intent(out) :: ok
        !> Speed of light; when present the relativistic kinetic energy is used.
        real(dp), intent(in), optional :: speed
        real(dp) :: theta, x
        type(interval_t) :: theta_bound, result(6), x_bound
        integer :: j
        state = fermi_si_state_t()
        ok = .false.
        if (.not. all(ieee_is_finite([density, temperature, mass, hbar, kb, &
                                     charge, epsilon]))) return
        if (min(density, temperature, mass, hbar, kb, charge, epsilon) <= 0.0_dp) &
            return
        call fermi_si_parameters(density, temperature, mass, hbar, kb, &
                                 state%fermi_momentum, state%fermi_energy, theta)
        call fermi_si_parameters_interval(interval(density), interval(temperature), &
            interval(mass), interval(hbar), interval(kb), &
            state%fermi_momentum_bound, state%fermi_energy_bound, theta_bound)
        if (.not. finite_interval(state%fermi_momentum_bound)) return
        if (.not. finite_interval(state%fermi_energy_bound)) return
        if (.not. finite_interval(theta_bound)) return
        theta_bound = hull(theta_bound, interval(theta))
        if (present(speed)) then
            if (.not. ieee_is_finite(speed) .or. speed <= 0.0_dp) return
            call fermi_relativity_indicator(state%fermi_momentum, mass, speed, x)
            call fermi_relativity_indicator_interval(state%fermi_momentum_bound, &
                interval(mass), interval(speed), x_bound)
            if (.not. ieee_is_finite(x) .or. .not. finite_interval(x_bound)) return
            call solve_fermi_gas(work, theta, state%normalized, ok, &
                theta_interval=theta_bound, relativity=x, &
                relativity_interval=hull(x_bound, interval(x)))
        else
            call solve_fermi_gas(work, theta, state%normalized, ok, &
                                 theta_interval=theta_bound)
        end if
        if (.not. ok) return
        associate (n => state%normalized)
            call fermi_si_eos(density, state%fermi_energy, n%mu, n%pressure, &
                n%energy, n%susceptibility, charge, epsilon, &
                state%chemical_potential, state%pressure, state%energy_density, &
                state%compressibility, state%screening_wavenumber, state%bulk_modulus)
            call fermi_si_eos_interval(interval(density), state%fermi_energy_bound, &
                n%mu_bound, n%pressure_bound, n%energy_bound, n%susceptibility_bound, &
                interval(charge), interval(epsilon), result(1), result(2), result(3), &
                result(4), result(5), result(6))
        end associate
        ok = .false.
        if (.not. all(ieee_is_finite([state%chemical_potential, state%pressure, &
            state%energy_density, state%compressibility, state%screening_wavenumber, &
            state%bulk_modulus]))) return
        do j = 1, 6
            if (.not. finite_interval(result(j))) return
        end do
        state%chemical_potential_bound = result(1); state%pressure_bound = result(2)
        state%energy_density_bound = result(3); state%compressibility_bound = result(4)
        state%screening_wavenumber_bound = result(5)
        state%bulk_modulus_bound = result(6)
        state%certified = .true.
        ok = .true.
    end subroutine solve_fermi_gas_si

    subroutine initialize_fermi_workspace(work, max_cells, ok)
        type(fermi_workspace_t), intent(out) :: work
        integer, intent(in) :: max_cells
        logical, intent(out) :: ok

        ok = .false.
        if (max_cells < 64 .or. max_cells >= huge(max_cells)) return
        allocate (work%edges(max_cells + 1), work%nodes(work%maximum_nodes), &
                  work%weights(work%maximum_nodes))
        work%initialized = .true.
        ok = .true.
    end subroutine initialize_fermi_workspace

    subroutine solve_fermi_gas(work, theta, state, ok, certify, theta_interval, &
                               relativity, relativity_interval)
        type(fermi_workspace_t), intent(inout) :: work
        real(dp), intent(in) :: theta
        type(fermi_state_t), intent(out) :: state
        logical, intent(out) :: ok
        logical, intent(in), optional :: certify
        type(interval_t), intent(in), optional :: theta_interval
        !> x=pF/(m*c) >= 0 and an optional enclosing interval. Absent or zero
        !> keeps the nonrelativistic leaves bit for bit.
        real(dp), intent(in), optional :: relativity
        type(interval_t), intent(in), optional :: relativity_interval
        real(dp) :: lower, upper, radius, xtol, step, mu, x0, f, stats(11), sums(11)
        real(dp) :: unused_upper, x
        type(interval_t) :: lo, hi, mu_band, left_residual, right_residual
        type(interval_t) :: enclosed(11), values(6), gradients(6), residual
        type(fortnum_status_t) :: root_status
        integer :: j, k, maximum_cells, moments, expansion
        logical :: check, need_candidate, rel

        state = fermi_state_t()
        ok = .false.
        state%status = 'invalid-input'
        if (.not. work%initialized) return
        if (.not. ieee_is_finite(theta)) return
        if (theta <= 0.0_dp) return
        state%theta = theta
        state%theta_bound = interval(theta)
        if (present(theta_interval)) then
            if (.not. finite_interval(theta_interval)) return
            if (theta_interval%lo <= 0.0_dp) return
            if (theta_interval%lo > theta .or. theta_interval%hi < theta) return
            state%theta_bound = theta_interval
        end if
        x = 0.0_dp
        if (present(relativity)) x = relativity
        if (.not. ieee_is_finite(x) .or. x < 0.0_dp) return
        state%relativity = x
        state%relativity_bound = interval(x)
        if (present(relativity_interval)) then
            if (.not. finite_interval(relativity_interval)) return
            if (relativity_interval%lo < 0.0_dp) return
            if (relativity_interval%lo > x .or. relativity_interval%hi < x) return
            state%relativity_bound = relativity_interval
        end if
        rel = state%relativity_bound%hi > 0.0_dp
        moments = 8
        if (rel) moments = 11
        call fermi_root_bracket(theta, lower, upper, radius)
        if (rel) call fermi_rel_root_bracket(theta, x, upper, radius)
        if (.not. all(ieee_is_finite([lower, upper, radius]))) return
        call fermi_scaled_tolerance(theta, 1.0e-7_dp, xtol)
        call fermi_scaled_tolerance(theta, 0.005_dp, step)
        if (.not. ieee_is_finite(step) .or. step <= 0.0_dp) return
        maximum_cells = size(work%edges) - 1
        if (work%cached_theta < theta .or. work%cached_theta > theta) work%cells = 64
        if (work%cached_relativity < x .or. work%cached_relativity > x) then
            work%cells = 64
            work%cached_relativity = x
        end if
        do
            call prepare_nodes(work, theta, radius, x)
            if (rel) then
                ! Relativistic density exceeds the nonrelativistic one at equal
                ! mu; extend the candidate lower bracket until the sign changes.
                state%status = 'candidate-bracket-failed'
                do expansion = 1, 64
                    call density_equation(lower, f, unused_upper)
                    if (f < 0.0_dp) exit
                    if (expansion == 64 .or. .not. ieee_is_finite(f)) return
                    lower = lower - (abs(lower) + 4.0_dp*theta)
                end do
            end if
            call fermi_half(lower, upper, x0)
            call root_newton(density_equation, lower, upper, x0, mu, root_status, &
                             xtol=xtol, ftol=1.0e-12_dp, max_iter=200, &
                             deriv_floor=1.0e-280_dp)
            state%status = 'candidate-root-failed'
            if (.not. status_ok(root_status)) return
            if (.not. ieee_is_finite(mu)) return
            sums = 0.0_dp
            stats = 0.0_dp
            do j = 1, work%candidate_nodes
                f = occupation_value(work%nodes(j), theta, mu, x, rel)
                if (rel) then
                    call fermi_rel_integrands(work%nodes(j), &
                        rel_root(work%nodes(j), x), theta, mu, x, f, stats(1), &
                        stats(2), stats(3), stats(4), stats(5), stats(6), stats(7), &
                        stats(8), stats(9), stats(10), stats(11))
                else
                    call fermi_integrands(work%nodes(j), theta, mu, f, stats(1), &
                        stats(2), stats(3), stats(4), stats(5), stats(6), stats(7), &
                        stats(8))
                end if
                do k = 1, moments
                    call fermi_weighted_add(sums(k), stats(k), work%weights(j), f)
                    sums(k) = f
                end do
            end do
            if (.not. all(ieee_is_finite(sums))) return
            if (rel) then
                call fermi_rel_eos(mu, sums(2), sums(3), sums(4), state%mu, &
                    state%pressure, state%energy, state%susceptibility, &
                    state%screening_ratio, state%bulk_modulus)
                call fermi_rel_theta_derivative(sums(1), sums(2), sums(3), sums(4), &
                    sums(5), sums(6), sums(7), sums(8), sums(9), sums(10), sums(11), &
                    state%theta_gradient(1), state%theta_gradient(2), &
                    state%theta_gradient(3), state%theta_gradient(4), &
                    state%theta_gradient(5), state%theta_gradient(6))
            else
                call fermi_eos(mu, sums(2), sums(3), state%mu, state%pressure, &
                    state%energy, state%susceptibility, state%screening_ratio, &
                    state%bulk_modulus)
                call fermi_theta_derivative(sums(1), sums(2), sums(3), sums(4), &
                    sums(5), sums(6), sums(7), sums(8), state%theta_gradient(1), &
                    state%theta_gradient(2), state%theta_gradient(3), &
                    state%theta_gradient(4), state%theta_gradient(5), &
                    state%theta_gradient(6))
            end if
            if (.not. all(ieee_is_finite(center_values(state)))) return
            if (.not. all(ieee_is_finite(state%theta_gradient))) return
            state%valid = .true.
            state%candidate_nodes = work%candidate_nodes
            state%status = 'candidate'
            ok = .true.
            if (present(certify)) then
                if (.not. certify) return
            end if
            ok = .false.
            call fermi_root_band_interval(interval(mu), interval(step), lo, hi)
            if (.not. finite_interval(lo) .or. .not. finite_interval(hi)) return
            mu_band = interval(lo%lo, hi%hi)
            need_candidate = .false.
            do
                call prepare_edges(work, radius, theta, x)
                call verified_moments(state%theta_bound, interval(mu_band%lo), &
                    state%relativity_bound, rel, radius, &
                    work%edges(:work%cells + 1), enclosed(:1), check)
                if (.not. check) exit
                call fermi_residual_interval(enclosed(1), left_residual)
                call verified_moments(state%theta_bound, interval(mu_band%hi), &
                    state%relativity_bound, rel, radius, &
                    work%edges(:work%cells + 1), enclosed(:1), check)
                if (.not. check) exit
                call fermi_residual_interval(enclosed(1), right_residual)
                if (left_residual%lo > 0.0_dp .or. right_residual%hi < 0.0_dp) then
                    need_candidate = .true.
                    exit
                end if
                if (left_residual%hi < 0.0_dp .and. right_residual%lo > 0.0_dp) then
                    call verified_moments(state%theta_bound, mu_band, &
                        state%relativity_bound, rel, radius, &
                        work%edges(:work%cells + 1), enclosed(:moments), check)
                    if (check) then
                        call enclosed_eos(mu_band, enclosed, rel, values, gradients, &
                                          check)
                        if (check) then
                            if (accuracy_satisfied(state, values)) then
                                call set_bounds(state, values, gradients)
                                call fermi_residual_interval(enclosed(1), residual)
                                call fermi_add_interval(residual, interval(1.0_dp), &
                                                        state%density_bound)
                                state%cells = work%cells
                                state%status = 'certified-continuous'
                                ok = .true.
                                return
                            end if
                            if (candidate_outside_budget(state, values) .or. &
                                work%cells > maximum_cells/2) then
                                need_candidate = .true.
                                exit
                            end if
                        end if
                    end if
                end if
                if (work%cells > maximum_cells/2) exit
                work%cells = work%cells*2
            end do
            state%status = 'certificate-budget-failed'
            if (.not. need_candidate) return
            if (work%candidate_nodes > work%maximum_nodes/2) return
            work%candidate_nodes = work%candidate_nodes*2
            work%cells = 64
        end do

    contains

        pure subroutine density_equation(candidate, value, derivative)
            real(dp), intent(in) :: candidate
            real(dp), intent(out) :: value, derivative
            real(dp) :: occupation, i2, i2mu, sum2, sum_mu, next
            integer :: node
            sum2 = 0.0_dp; sum_mu = 0.0_dp
            do node = 1, work%candidate_nodes
                occupation = occupation_value(work%nodes(node), theta, candidate, x, rel)
                if (rel) then
                    call fermi_rel_density_integrand(work%nodes(node), &
                        rel_root(work%nodes(node), x), theta, candidate, x, &
                        occupation, i2, i2mu)
                else
                    call fermi_density_integrand(work%nodes(node), theta, candidate, &
                                                 occupation, i2, i2mu)
                end if
                call fermi_weighted_add(sum2, i2, work%weights(node), next)
                sum2 = next
                call fermi_weighted_add(sum_mu, i2mu, work%weights(node), next)
                sum_mu = next
            end do
            call fermi_residual(sum2, value)
            call fermi_scaled_tolerance(sum_mu, 3.0_dp, derivative)
        end subroutine density_equation

    end subroutine solve_fermi_gas

    subroutine prepare_nodes(work, theta, radius, x)
        type(fermi_workspace_t), intent(inout) :: work
        real(dp), intent(in) :: theta, radius, x
        real(dp) :: window(5)
        integer :: count, first, j
        logical :: graded
        if (work%cached_nodes == work%candidate_nodes) then
            if (work%cached_theta <= theta .and. work%cached_theta >= theta) then
                if (work%cached_radius <= radius .and. work%cached_radius >= radius &
                    .and. work%nodes_relativity <= x .and. work%nodes_relativity >= x) &
                    return
            end if
        end if
        call graded_window(theta, radius, x, window, graded)
        if (graded) then
            count = work%candidate_nodes/4
            first = 1
            do j = 1, 4
                call gauss_legendre_ab(count, window(j), window(j + 1), &
                    work%nodes(first:first + count - 1), &
                    work%weights(first:first + count - 1))
                first = first + count
            end do
        else
            call gauss_legendre_ab(work%candidate_nodes, 0.0_dp, radius, &
                work%nodes(:work%candidate_nodes), work%weights(:work%candidate_nodes))
        end if
        work%cached_theta = theta
        work%cached_radius = radius
        work%cached_nodes = work%candidate_nodes
        work%nodes_relativity = x
    end subroutine prepare_nodes

    subroutine prepare_edges(work, radius, theta, x)
        type(fermi_workspace_t), intent(inout) :: work
        real(dp), intent(in) :: radius, theta, x
        real(dp) :: window(5)
        integer :: j, segment, count, first
        logical :: graded
        call graded_window(theta, radius, x, window, graded)
        if (graded) then
            count = work%cells/4
            first = 1
            do segment = 1, 4
                do j = 0, count - 1
                    call fermi_rel_grid(real(j, dp), real(count, dp), &
                        window(segment), window(segment + 1), work%edges(first + j))
                end do
                first = first + count
            end do
        else
            do j = 1, work%cells + 1
                call fermi_grid(real(j - 1, dp), real(work%cells, dp), radius, &
                                work%edges(j))
            end do
        end if
        work%edges(1) = 0.0_dp
        work%edges(work%cells + 1) = radius
    end subroutine prepare_edges

    !> Degenerate states (relativistic or not) use a four-segment mesh with
    !> equal counts, split at the cold Fermi surface r=1 inside a window of
    !> fermi_rel_window halfwidth. One Gauss rule across the whole window
    !> leaves the occupation-derivative peak at its sparse centre (relative
    !> susceptibility error 1.6e-3 for 64 nodes); splitting at the peak puts
    !> it where the nodes cluster. Only cell placement changes, not bounds.
    pure subroutine graded_window(theta, radius, x, window, use)
        real(dp), intent(in) :: theta, radius, x
        real(dp), intent(out) :: window(5)
        logical, intent(out) :: use
        real(dp) :: halfwidth
        use = .false.
        window = 0.0_dp
        call fermi_rel_window(theta, x, halfwidth)
        if (.not. ieee_is_finite(halfwidth) .or. halfwidth >= 1.0_dp) return
        call fermi_root_band(1.0_dp, halfwidth, window(2), window(4))
        window(3) = 1.0_dp
        window(5) = radius
        if (.not. (window(2) > 0.0_dp .and. window(4) > window(2) .and. &
                   radius > window(4))) return
        use = .true.
    end subroutine graded_window

    pure real(dp) function rel_root(r, x) result(s)
        real(dp), intent(in) :: r, x
        call fermi_rel_root(r, x, s)
    end function rel_root

    pure real(dp) function occupation_value(r, theta, mu, x, rel) result(value)
        real(dp), intent(in) :: r, theta, mu, x
        logical, intent(in) :: rel
        real(dp) :: argument
        if (rel) then
            call fermi_rel_argument(r, rel_root(r, x), theta, mu, argument)
        else
            call fermi_argument(r, theta, mu, argument)
        end if
        if (argument >= 0.0_dp) then
            call fermi_positive(argument, value)
        else
            call fermi_negative(argument, value)
        end if
    end function occupation_value

    pure subroutine occupation_packet(r, theta, mu, x, rel, value, ok)
        type(interval_t), intent(in) :: r, theta, mu, x
        logical, intent(in) :: rel
        type(interval_t), intent(out) :: value
        logical, intent(out) :: ok
        type(interval_t) :: argument, lower, upper, s
        ok = .false.
        if (rel) then
            call fermi_rel_root_interval(r, x, s)
            if (.not. finite_interval(s)) return
            call fermi_rel_argument_interval(r, s, theta, mu, argument)
        else
            call fermi_argument_interval(r, theta, mu, argument)
        end if
        if (.not. finite_interval(argument)) return
        if (argument%hi >= 0.0_dp) then
            call fermi_positive_interval(interval(argument%hi), lower)
        else
            call fermi_negative_interval(interval(argument%hi), lower)
        end if
        if (argument%lo >= 0.0_dp) then
            call fermi_positive_interval(interval(argument%lo), upper)
        else
            call fermi_negative_interval(interval(argument%lo), upper)
        end if
        if (.not. finite_interval(lower) .or. .not. finite_interval(upper)) return
        value = intersect(interval(lower%lo, upper%hi), interval(0.0_dp, 1.0_dp))
        ok = finite_interval(value)
    end subroutine occupation_packet

    subroutine verified_moments(theta, mu, x, rel, radius, edges, integral, ok)
        real(dp), intent(in) :: radius, edges(:)
        type(interval_t), intent(in) :: theta, mu, x
        logical, intent(in) :: rel
        type(interval_t), intent(out) :: integral(:)
        logical, intent(out) :: ok
        type(interval_t) :: tails(11), addition, next, s
        integer :: j
        logical :: positive
        call composite_midpoint_batch(evaluate_cell, edges, integral, ok)
        if (.not. ok) return
        if (rel) then
            call fermi_rel_root_interval(interval(radius), x, s)
            call fermi_rel_tail_interval(theta, mu, interval(radius), s, x, &
                tails(1), tails(2), tails(3), tails(4), tails(5), tails(6), &
                tails(7), tails(8), tails(9), tails(10), tails(11))
        else
            call fermi_tail_interval(theta, mu, interval(radius), &
                tails(1), tails(2), tails(3), tails(4), tails(5), tails(6), &
                tails(7), tails(8))
        end if
        do j = 1, size(integral)
            if (.not. finite_interval(tails(j))) then
                ok = .false.
                return
            end if
            if (rel) then
                positive = any(j == [1, 2, 3, 4, 6, 8])
            else
                positive = any(j == [1, 2, 3, 5])
            end if
            if (positive) then
                addition = interval(0.0_dp, tails(j)%hi)
            else
                call fermi_root_band_interval(interval(0.0_dp), &
                    interval(tails(j)%hi), addition, next)
                addition = hull(addition, next)
            end if
            call fermi_add_interval(integral(j), addition, next)
            integral(j) = next
            if (.not. finite_interval(integral(j))) then
                ok = .false.
                return
            end if
        end do

    contains

        subroutine evaluate_cell(cell, midpoint, value, curvature, callback_ok)
            type(interval_t), intent(in) :: cell, midpoint
            type(interval_t), intent(out) :: value(:), curvature(:)
            logical, intent(out) :: callback_ok
            type(interval_t) :: midpoint_f, cell_f, raw(11), unused, s_mid, s_cell
            integer :: component
            call occupation_packet(midpoint, theta, mu, x, rel, midpoint_f, &
                                   callback_ok)
            if (.not. callback_ok) return
            call occupation_packet(cell, theta, mu, x, rel, cell_f, callback_ok)
            if (.not. callback_ok) return
            if (rel) then
                call fermi_rel_root_interval(midpoint, x, s_mid)
                call fermi_rel_root_interval(cell, x, s_cell)
            end if
            if (size(value) == 1) then
                if (rel) then
                    call fermi_rel_density_integrand_interval(midpoint, s_mid, theta, &
                        mu, x, midpoint_f, value(1), unused)
                    call fermi_rel_density_curvature_interval(cell, s_cell, theta, &
                        mu, x, cell_f, raw(1))
                else
                    call fermi_density_integrand_interval(midpoint, theta, &
                        mu, midpoint_f, value(1), unused)
                    call fermi_density_curvature_interval(cell, theta, &
                                                         mu, cell_f, raw(1))
                end if
            else if (rel) then
                call fermi_rel_integrands_interval(midpoint, s_mid, theta, mu, x, &
                    midpoint_f, value(1), value(2), value(3), value(4), value(5), &
                    value(6), value(7), value(8), value(9), value(10), value(11))
                call fermi_rel_curvature_interval(cell, s_cell, theta, mu, x, cell_f, &
                    raw(1), raw(2), raw(3), raw(4), raw(5), raw(6), raw(7), raw(8), &
                    raw(9), raw(10), raw(11))
            else
                call fermi_integrands_interval(midpoint, theta, mu, &
                    midpoint_f, value(1), value(2), value(3), value(4), value(5), &
                    value(6), value(7), value(8))
                call fermi_curvature_interval(cell, theta, mu, cell_f, &
                    raw(1), raw(2), raw(3), raw(4), raw(5), raw(6), raw(7), raw(8))
            end if
            do component = 1, size(value)
                if (.not. finite_interval(value(component))) then
                    callback_ok = .false.
                    return
                end if
                if (.not. finite_interval(raw(component))) then
                    callback_ok = .false.
                    return
                end if
                curvature(component) = interval(0.0_dp, mag(raw(component)))
            end do
        end subroutine evaluate_cell

    end subroutine verified_moments

    pure subroutine enclosed_eos(mu, integral, rel, value, gradient, ok)
        type(interval_t), intent(in) :: mu, integral(11)
        logical, intent(in) :: rel
        type(interval_t), intent(out) :: value(6), gradient(6)
        logical, intent(out) :: ok
        integer :: j
        ok = .false.
        if (rel) then
            if (integral(4)%lo <= 0.0_dp) return
            call fermi_rel_eos_interval(mu, integral(2), integral(3), integral(4), &
                value(1), value(2), value(3), value(4), value(5), value(6))
            call fermi_rel_theta_derivative_interval(integral(1), integral(2), &
                integral(3), integral(4), integral(5), integral(6), integral(7), &
                integral(8), integral(9), integral(10), integral(11), gradient(1), &
                gradient(2), gradient(3), gradient(4), gradient(5), gradient(6))
        else
            if (integral(3)%lo <= 0.0_dp) return
            call fermi_eos_interval(mu, integral(2), integral(3), value(1), value(2), &
                value(3), value(4), value(5), value(6))
            call fermi_theta_derivative_interval(integral(1), integral(2), &
                integral(3), integral(4), integral(5), integral(6), integral(7), &
                integral(8), gradient(1), gradient(2), gradient(3), gradient(4), &
                gradient(5), gradient(6))
        end if
        do j = 1, 6
            if (.not. finite_interval(value(j))) return
            if (.not. finite_interval(gradient(j))) return
        end do
        ok = .true.
    end subroutine enclosed_eos

    pure function center_values(state) result(value)
        type(fermi_state_t), intent(in) :: state
        real(dp) :: value(6)
        value = [state%mu, state%pressure, state%energy, state%susceptibility, &
                 state%screening_ratio, state%bulk_modulus]
    end function center_values

    pure logical function accuracy_satisfied(state, value) result(ok)
        type(fermi_state_t), intent(in) :: state
        type(interval_t), intent(in) :: value(6)
        type(interval_t) :: difference
        real(dp) :: center(6), scale, budget
        integer :: j
        ok = .false.
        center = center_values(state)
        do j = 1, 6
            call fermi_difference_interval(value(j), interval(center(j)), difference)
            scale = mag(interval(center(j)))
            if (scale < 1.0_dp) scale = 1.0_dp
            call fermi_scaled_tolerance(scale, 0.01_dp, budget)
            if (mag(difference) > budget) return
        end do
        ok = .true.
    end function accuracy_satisfied

    pure logical function candidate_outside_budget(state, value) result(outside)
        type(fermi_state_t), intent(in) :: state
        type(interval_t), intent(in) :: value(6)
        type(interval_t) :: difference
        real(dp) :: center(6), scale, budget
        integer :: j
        outside = .false.
        center = center_values(state)
        do j = 1, 6
            call fermi_difference_interval(value(j), interval(center(j)), difference)
            scale = mag(interval(center(j)))
            if (scale < 1.0_dp) scale = 1.0_dp
            call fermi_scaled_tolerance(scale, 0.01_dp, budget)
            if (difference%lo > budget .or. difference%hi < -budget) then
                outside = .true.
                return
            end if
        end do
    end function candidate_outside_budget

    subroutine set_bounds(state, value, gradient)
        type(fermi_state_t), intent(inout) :: state
        type(interval_t), intent(in) :: value(6), gradient(6)
        type(interval_t) :: displayed(6), displayed_gradient(6)
        real(dp) :: center(6)
        integer :: j
        center = center_values(state)
        do j = 1, 6
            displayed(j) = hull(value(j), interval(center(j)))
            displayed_gradient(j) = hull(gradient(j), interval(state%theta_gradient(j)))
        end do
        state%mu_bound = displayed(1); state%pressure_bound = displayed(2)
        state%energy_bound = displayed(3); state%susceptibility_bound = displayed(4)
        state%screening_ratio_bound = displayed(5)
        state%bulk_modulus_bound = displayed(6)
        state%theta_gradient_bound = displayed_gradient
        state%certified = .true.
        state%derivative_certified = .true.
        state%prepared_theta = state%theta
        state%prepared_theta_bound = state%theta_bound
        state%prepared_relativity = state%relativity_bound
        state%prepared_center = center_values(state)
        state%prepared_gradient = state%theta_gradient
        state%prepared_value = displayed
        state%prepared_derivative = displayed_gradient
        state%prepared = .true.
    end subroutine set_bounds

    subroutine occupation_fermi_gas(state, q, value, lower, upper, ok)
        type(fermi_state_t), intent(in) :: state
        real(dp), intent(in) :: q
        real(dp), intent(out) :: value, lower, upper
        logical, intent(out) :: ok
        type(interval_t) :: packet
        ok = .false.
        value = 0.0_dp; lower = 0.0_dp; upper = 0.0_dp
        if (.not. packet_valid(state)) return
        if (.not. ieee_is_finite(q)) return
        if (q < 0.0_dp) return
        call occupation_packet(interval(q), state%theta_bound, state%mu_bound, &
            state%relativity_bound, state%relativity_bound%hi > 0.0_dp, packet, ok)
        if (.not. ok) return
        value = occupation_value(q, state%theta, state%mu, state%relativity, &
                                 state%relativity_bound%hi > 0.0_dp)
        if (.not. ieee_is_finite(value)) then
            ok = .false.
            return
        end if
        packet = hull(packet, interval(value))
        lower = packet%lo; upper = packet%hi
    end subroutine occupation_fermi_gas

    subroutine fermi_gas_jvp(state, dtheta, tangent, tangent_bound, ok)
        type(fermi_state_t), intent(in) :: state
        real(dp), intent(in) :: dtheta
        real(dp), intent(out) :: tangent(6)
        type(interval_t), intent(out) :: tangent_bound(6)
        logical, intent(out) :: ok
        real(dp) :: g(6)
        type(interval_t) :: b(6)
        integer :: j
        ok = .false.
        if (.not. packet_valid(state)) return
        if (.not. ieee_is_finite(dtheta)) return
        g = state%theta_gradient; b = state%theta_gradient_bound
        call fermi_apply_direction(dtheta, g(1), g(2), g(3), g(4), g(5), g(6), &
            tangent(1), tangent(2), tangent(3), tangent(4), tangent(5), tangent(6))
        call fermi_apply_direction_interval(interval(dtheta), b(1), b(2), b(3), &
            b(4), b(5), b(6), tangent_bound(1), tangent_bound(2), tangent_bound(3), &
            tangent_bound(4), tangent_bound(5), tangent_bound(6))
        if (.not. all(ieee_is_finite(tangent))) return
        do j = 1, 6
            if (.not. finite_interval(tangent_bound(j))) return
            tangent_bound(j) = hull(tangent_bound(j), interval(tangent(j)))
        end do
        ok = .true.
    end subroutine fermi_gas_jvp

    subroutine fermi_gas_vjp(state, cotangent, gtheta, gtheta_bound, ok)
        type(fermi_state_t), intent(in) :: state
        real(dp), intent(in) :: cotangent(6)
        real(dp), intent(out) :: gtheta
        type(interval_t), intent(out) :: gtheta_bound
        logical, intent(out) :: ok
        real(dp) :: g(6), ignored(6)
        type(interval_t) :: b(6), weights(6), unused(6)
        ok = .false.
        if (.not. packet_valid(state)) return
        if (.not. all(ieee_is_finite(cotangent))) return
        g = state%theta_gradient; b = state%theta_gradient_bound
        weights = interval(cotangent)
        call fermi_apply_direction_vjp(0.0_dp, g(1), g(2), g(3), g(4), g(5), g(6), &
            cotangent(1), cotangent(2), cotangent(3), cotangent(4), cotangent(5), &
            cotangent(6), gtheta, ignored(1), ignored(2), ignored(3), ignored(4), &
            ignored(5), ignored(6))
        call fermi_apply_direction_vjp_interval(interval(0.0_dp), b(1), b(2), b(3), &
            b(4), b(5), b(6), weights(1), weights(2), weights(3), weights(4), &
            weights(5), weights(6), gtheta_bound, unused(1), unused(2), unused(3), &
            unused(4), unused(5), unused(6))
        if (.not. ieee_is_finite(gtheta)) return
        if (.not. finite_interval(gtheta_bound)) return
        gtheta_bound = hull(gtheta_bound, interval(gtheta))
        ok = .true.
    end subroutine fermi_gas_vjp

    pure logical function packet_valid(state) result(ok)
        type(fermi_state_t), intent(in) :: state
        type(interval_t) :: bounds(6)
        real(dp) :: center(6)
        integer :: j
        ok = .false.
        if (.not. state%prepared) return
        if (.not. state%certified .or. .not. state%derivative_certified) return
        if (.not. ieee_is_finite(state%theta)) return
        if (state%theta < state%prepared_theta .or. &
            state%theta > state%prepared_theta) return
        if (.not. same_interval(state%theta_bound, state%prepared_theta_bound)) return
        if (.not. same_interval(state%relativity_bound, state%prepared_relativity)) &
            return
        if (state%relativity < state%relativity_bound%lo .or. &
            state%relativity > state%relativity_bound%hi) return
        center = center_values(state)
        if (.not. all(ieee_is_finite(center))) return
        if (.not. all(ieee_is_finite(state%theta_gradient))) return
        bounds = [state%mu_bound, state%pressure_bound, state%energy_bound, &
                  state%susceptibility_bound, state%screening_ratio_bound, &
                  state%bulk_modulus_bound]
        do j = 1, 6
            if (center(j) < state%prepared_center(j) .or. &
                center(j) > state%prepared_center(j)) return
            if (state%theta_gradient(j) < state%prepared_gradient(j) .or. &
                state%theta_gradient(j) > state%prepared_gradient(j)) return
            if (.not. same_interval(bounds(j), state%prepared_value(j))) return
            if (.not. same_interval(state%theta_gradient_bound(j), &
                                    state%prepared_derivative(j))) return
        end do
        ok = .true.
    end function packet_valid

    pure logical function same_interval(a, b) result(ok)
        type(interval_t), intent(in) :: a, b
        ok = .false.
        if (.not. finite_interval(a) .or. .not. finite_interval(b)) return
        if (a%lo < b%lo .or. a%lo > b%lo) return
        if (a%hi < b%hi .or. a%hi > b%hi) return
        ok = .true.
    end function same_interval

    pure logical function finite_interval(value) result(ok)
        type(interval_t), intent(in) :: value
        ok = .false.
        if (is_empty(value)) return
        if (.not. ieee_is_finite(value%lo) .or. .not. ieee_is_finite(value%hi)) return
        ok = .true.
    end function finite_interval

end module kin6d_fermi_gas
