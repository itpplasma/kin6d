!> Q7 part 1: static screening hierarchy of an ideal electron gas with a
!> neutralizing background (Debye, finite-T Fermi/Thomas-Fermi, Lindhard/RPA).
!> Units kF = 1, EF = 1, x = q/(2 kF); G(x) = -chi0(q)/N(0). Static ideal-gas
!> response only: no local-field correction, no correlations.
module kin6d_screening_hierarchy
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortnum_interval, only: interval_t, interval
    use fortnum_quadrature, only: gauss_legendre
    use kin6d_species, only: units_t, species_t, plasma_scales_t, plasma_scales, &
        pair_coupling
    use kin6d_fermi_gas, only: fermi_workspace_t, fermi_state_t, fermi_si_state_t, &
        initialize_fermi_workspace, solve_fermi_gas, solve_fermi_gas_si
    use kin6d_fermi_gas_leaves, only: fermi_argument, fermi_positive, fermi_negative
    use kin6d_screening_leaves, only: screening_lindhard_inner, &
        screening_lindhard_outer, screening_lindhard_small, screening_lindhard_large, &
        screening_momentum_below, screening_momentum_above, screening_maldague_weight, &
        screening_maldague_term, screening_potential_term, screening_potential_term_jvp, &
        screening_momentum_tail_interval, screening_maldague_tail_interval, &
        screening_inversion_tail_interval
    implicit none
    private
    public :: screening_state_t, screening_value_t, screening_setup, &
        screening_setup_normalized, lindhard_zero, lindhard_momentum, &
        lindhard_maldague, lindhard_response, screening_q_table_t, &
        screening_q_table, screened_ratio, screened_ratio_from_table

    integer, parameter :: low_order = 12, high_order = 20, max_levels = 48
    real(dp), parameter :: grading = 0.25_dp, pi = acos(-1.0_dp)

    !> Static ideal-gas screening state at (theta, kappa2). theta = 0 is the
    !> exact T=0 limit (closed-form Lindhard); theta > 0 uses mu from kin6d_fermi_gas.
    type :: screening_state_t
        real(dp) :: theta = 0.0_dp, mu = 1.0_dp, kappa2 = 0.0_dp
        !> ks/kTF0 = sqrt(G(0)) from kin6d_fermi_gas; kD/kTF0 = sqrt(2/(3 theta)).
        real(dp) :: ks_ratio = 1.0_dp, kd_ratio = huge(1.0_dp)
        !> Certified kin6d_fermi_gas enclosure of ks/kTF0 (1% output budget).
        type(interval_t) :: ks_ratio_bound
        !> Physical scales in the caller's reduced units (0 if normalized setup).
        real(dp) :: kf = 0.0_dp, ktf0 = 0.0_dp, ks = 0.0_dp, kd = 0.0_dp
        type(plasma_scales_t) :: scales
        real(dp) :: xl(low_order) = 0.0_dp, wl(low_order) = 0.0_dp
        real(dp) :: xh(high_order) = 0.0_dp, wh(high_order) = 0.0_dp
        logical :: ready = .false., zero_temperature = .true.
    end type

    !> Candidate value (high-order composite Gauss), rigorous truncation tail
    !> bound and quadrature indicator |Q20 - Q12| (an indicator, not a bound).
    type :: screening_value_t
        real(dp) :: value = 0.0_dp, tail = 0.0_dp, indicator = 0.0_dp
    end type

    !> Composite q rule for Fourier inversion with G tabulated on its nodes.
    type :: screening_q_table_t
        real(dp), allocatable :: ql(:), wl(:), gl(:), qh(:), wh(:), gh(:)
        real(dp) :: qmax = 0.0_dp, tail = 0.0_dp, g_indicator = 0.0_dp
    end type

contains

    subroutine prepare_rules(state)
        type(screening_state_t), intent(inout) :: state
        call gauss_legendre(low_order, state%xl, state%wl)
        call gauss_legendre(high_order, state%xh, state%wh)
        state%ready = .true.
    end subroutine

    !> Physical setup from density and temperature (energy units, k_B = 1)
    !> through the shared species record and the finite-T ideal Fermi EOS.
    subroutine screening_setup(units, s, density, temperature, state, ok)
        type(units_t), intent(in) :: units
        type(species_t), intent(in) :: s
        real(dp), intent(in) :: density, temperature
        type(screening_state_t), intent(out) :: state
        logical, intent(out) :: ok
        type(fermi_workspace_t) :: work
        type(fermi_si_state_t) :: si
        real(dp) :: g
        ok = .false.
        state%scales = plasma_scales(units, s, density, temperature)
        g = pair_coupling(units, s, s)
        call initialize_fermi_workspace(work, 4096, ok)
        if (.not. ok) return
        ! k_s^2 = (q^2/eps0) dn/dmu with q^2/eps0 = 4 pi g: charge sqrt(g), eps0 1/(4 pi).
        call solve_fermi_gas_si(work, density, temperature, s%mass, units%hbar, &
                                1.0_dp, sqrt(abs(g)), 1.0_dp/(4*pi), si, ok)
        if (.not. ok) return
        state%theta = si%normalized%theta
        state%zero_temperature = .false.
        state%mu = si%normalized%mu
        state%ks_ratio = si%normalized%screening_ratio
        state%ks_ratio_bound = si%normalized%screening_ratio_bound
        state%kf = si%fermi_momentum/units%hbar
        state%ks = si%screening_wavenumber
        state%ktf0 = state%ks/state%ks_ratio
        state%kd = 1.0_dp/state%scales%debye_length
        state%kd_ratio = state%kd/state%ktf0
        state%kappa2 = (state%ktf0/state%kf)**2
        call prepare_rules(state)
        ok = all(ieee_is_finite([state%kappa2, state%mu, state%ks_ratio]))
    end subroutine

    !> Normalized setup: theta = T/EF >= 0 and kappa2 = (kTF0/kF)^2.
    subroutine screening_setup_normalized(theta, kappa2, state, ok)
        real(dp), intent(in) :: theta, kappa2
        type(screening_state_t), intent(out) :: state
        logical, intent(out) :: ok
        type(fermi_workspace_t) :: work
        type(fermi_state_t) :: fs
        ok = .false.
        if (theta < 0.0_dp .or. kappa2 <= 0.0_dp) return
        state%theta = theta
        state%kappa2 = kappa2
        state%zero_temperature = .not. (theta > 0.0_dp)
        if (theta > 0.0_dp) then
            call initialize_fermi_workspace(work, 4096, ok)
            if (.not. ok) return
            call solve_fermi_gas(work, theta, fs, ok)
            if (.not. ok) return
            state%mu = fs%mu
            state%ks_ratio = fs%screening_ratio
            state%ks_ratio_bound = fs%screening_ratio_bound
            state%kd_ratio = sqrt(2.0_dp/(3.0_dp*theta))
        end if
        call prepare_rules(state)
        ok = .true.
    end subroutine

    !> Exact T=0 Lindhard function F(x), x = q/(2 kF), from four generated charts.
    pure real(dp) function lindhard_zero(x) result(f)
        real(dp), intent(in) :: x
        if (x <= 0.1_dp) then
            call screening_lindhard_small(x, f)
        else if (x < 1.0_dp) then
            call screening_lindhard_inner(x, f)
        else if (x > 1.0_dp .and. x < 10.0_dp) then
            call screening_lindhard_outer(x, f)
        else if (x >= 10.0_dp) then
            call screening_lindhard_large(x, f)
        else
            f = 0.5_dp
        end if
    end function

    !> Fermi occupation at r = p/pF from the authoritative kin6d_fermi_gas leaves.
    pure real(dp) function occupation(r, theta, mu) result(f)
        real(dp), intent(in) :: r, theta, mu
        real(dp) :: z
        call fermi_argument(r, theta, mu, z)
        if (z >= 0.0_dp) then
            call fermi_positive(z, f)
        else
            call fermi_negative(z, f)
        end if
    end function

    !> G(x) by the default route: closed form at theta = 0, Maldague otherwise.
    function lindhard_response(state, x) result(res)
        type(screening_state_t), intent(in) :: state
        real(dp), intent(in) :: x
        type(screening_value_t) :: res
        if (state%zero_temperature) then
            res%value = lindhard_zero(x)
        else
            res = lindhard_maldague(state, x)
        end if
    end function

    !> Route A: G = integral r f(r) log|(r+x)/(r-x)| dr/(2x), log singularity
    !> at r = x resolved by geometric panel grading. theta = 0: step occupation.
    function lindhard_momentum(state, x) result(res)
        type(screening_state_t), intent(in) :: state
        real(dp), intent(in) :: x
        type(screening_value_t) :: res
        real(dp) :: radius, h, q(2)
        type(interval_t) :: bound
        if (state%zero_temperature) then
            radius = 1.0_dp
            h = 0.05_dp
        else
            radius = max(2*x, sqrt(max(state%mu, 0.0_dp) + 40*state%theta))
            h = min(0.05_dp, state%theta)
            call screening_momentum_tail_interval(interval(state%theta), &
                interval(state%mu), interval(radius), bound)
            res%tail = bound%hi
        end if
        q = 0.0_dp
        if (x < radius) then
            call graded_sum(state, 0.0_dp, x, 2, h, integrand, q)
            call graded_sum(state, x, radius, 1, h, integrand, q)
        else
            call graded_sum(state, 0.0_dp, radius, 2, h, integrand, q)
        end if
        res%value = q(2)
        res%indicator = abs(q(2) - q(1))
    contains
        real(dp) function integrand(r) result(value)
            real(dp), intent(in) :: r
            real(dp) :: f
            if (state%zero_temperature) then
                f = 1.0_dp
            else
                f = occupation(r, state%theta, state%mu)
            end if
            if (r < x) then
                call screening_momentum_below(r, x, f, value)
            else
                call screening_momentum_above(r, x, f, value)
            end if
        end function
    end function

    !> Route B (Maldague): G = integral w(e) sqrt(e) F(x/sqrt(e)) de over the
    !> exact T=0 response, e = s^2; the kink of F at s = x is graded.
    function lindhard_maldague(state, x) result(res)
        type(screening_state_t), intent(in) :: state
        real(dp), intent(in) :: x
        type(screening_value_t) :: res
        real(dp) :: smax, h, q(2)
        type(interval_t) :: bound
        if (state%theta <= 0.0_dp) error stop 'lindhard_maldague: theta must be > 0'
        smax = sqrt(max(state%mu, 0.0_dp) + 40*state%theta)
        h = min(0.05_dp, state%theta)
        call screening_maldague_tail_interval(interval(state%theta), &
            interval(state%mu), interval(smax), bound)
        res%tail = bound%hi
        q = 0.0_dp
        if (x < smax) then
            call graded_sum(state, 0.0_dp, x, 2, h, integrand, q)
            call graded_sum(state, x, smax, 1, h, integrand, q)
        else
            call graded_sum(state, 0.0_dp, smax, 0, h, integrand, q)
        end if
        res%value = q(2)
        res%indicator = abs(q(2) - q(1))
    contains
        real(dp) function integrand(s) result(value)
            real(dp), intent(in) :: s
            real(dp) :: w
            call screening_maldague_weight(s*s, state%mu, state%theta, w)
            call screening_maldague_term(s, w, lindhard_zero(x/s), value)
        end function
    end function

    !> Composite Gauss sums (12 and 20 nodes per panel) on [a,b]; toward = 1/2
    !> grades panels geometrically toward a/b, 0 none; panels are split to
    !> width <= step.
    subroutine graded_sum(state, a, b, toward, step, f, total)
        type(screening_state_t), intent(in) :: state
        real(dp), intent(in) :: a, b, step
        integer, intent(in) :: toward
        real(dp), intent(inout) :: total(2)
        interface
            real(dp) function f(r)
                import :: dp
                real(dp), intent(in) :: r
            end function
        end interface
        real(dp) :: lo, hi, width
        integer :: level, levels
        if (b <= a) return
        width = b - a
        levels = grading_levels(a, b, toward)
        do level = 1, levels
            if (toward == 0) then
                lo = a; hi = b
            else if (toward == 1) then
                lo = a + width*grading**level; hi = a + width*grading**(level - 1)
                if (level == levels) lo = a
            else
                lo = b - width*grading**(level - 1); hi = b - width*grading**level
                if (level == levels) hi = b
            end if
            call uniform_sum(lo, hi)
        end do
    contains
        subroutine uniform_sum(p, r)
            real(dp), intent(in) :: p, r
            integer :: cells, c, j
            real(dp) :: left, half, mid
            cells = max(1, ceiling((r - p)/step))
            do c = 1, cells
                left = p + (r - p)*(c - 1)/cells
                half = 0.5_dp*(r - p)/cells
                mid = left + half
                do j = 1, low_order
                    total(1) = total(1) + half*state%wl(j)*f(mid + half*state%xl(j))
                end do
                do j = 1, high_order
                    total(2) = total(2) + half*state%wh(j)*f(mid + half*state%xh(j))
                end do
            end do
        end subroutine
    end subroutine

    !> Number of geometric levels: the innermost panel keeps a width of at
    !> least 4096 ulp of its singular endpoint so Gauss nodes stay distinct.
    pure integer function grading_levels(a, b, toward) result(levels)
        real(dp), intent(in) :: a, b
        integer, intent(in) :: toward
        real(dp) :: floor_width
        if (toward == 0) then
            levels = 1
            return
        end if
        floor_width = 4096*epsilon(1.0_dp)*max(abs(a), abs(b), 1.0_dp)
        levels = 1
        do while (levels < max_levels .and. (b - a)*grading**levels > floor_width)
            levels = levels + 1
        end do
    end function

    !> Nodes of the composite q rule on [0, qmax] (kF = 1): panels graded toward
    !> the 2kF kink, width <= near/far, and G tabulated on them by the default
    !> route. The inversion tail bound uses G(q/2) <= 8/(3 q^2).
    subroutine screening_q_table(state, qmax, near, far, table, constant_g)
        type(screening_state_t), intent(in) :: state
        real(dp), intent(in) :: qmax, near, far
        type(screening_q_table_t), intent(out) :: table
        !> Optional synthetic G (oracle hook: constant G gives exact Yukawa).
        real(dp), intent(in), optional :: constant_g
        integer :: nl, nh, j
        type(screening_value_t) :: res
        type(interval_t) :: bound
        call count_nodes(nl, nh)
        allocate (table%ql(nl), table%wl(nl), table%gl(nl))
        allocate (table%qh(nh), table%wh(nh), table%gh(nh))
        nl = 0; nh = 0
        call build(0.0_dp, 2.0_dp, 2, near, .true.)
        call build(2.0_dp, 8.0_dp, 1, near, .true.)
        call build(8.0_dp, qmax, 0, far, .true.)
        do j = 1, size(table%ql)
            table%gl(j) = g_of(table%ql(j))
        end do
        do j = 1, size(table%qh)
            table%gh(j) = g_of(table%qh(j))
        end do
        table%qmax = qmax
        call screening_inversion_tail_interval(interval(state%kappa2), &
            interval(qmax), bound)
        table%tail = bound%hi
    contains
        real(dp) function g_of(q) result(g)
            real(dp), intent(in) :: q
            if (present(constant_g)) then
                g = constant_g
            else
                res = lindhard_response(state, 0.5_dp*q)
                g = res%value
                table%g_indicator = max(table%g_indicator, res%indicator + res%tail)
            end if
        end function
        subroutine count_nodes(cl, ch)
            integer, intent(out) :: cl, ch
            nl = 0; nh = 0
            call build(0.0_dp, 2.0_dp, 2, near, .false.)
            call build(2.0_dp, 8.0_dp, 1, near, .false.)
            call build(8.0_dp, qmax, 0, far, .false.)
            cl = nl; ch = nh
        end subroutine
        subroutine build(a, b, toward, step, store)
            real(dp), intent(in) :: a, b, step
            integer, intent(in) :: toward
            logical, intent(in) :: store
            real(dp) :: lo, hi, width, left, half
            integer :: level, levels, cells, c
            width = b - a
            levels = grading_levels(a, b, toward)
            do level = 1, levels
                if (toward == 0) then
                    lo = a; hi = b
                else if (toward == 1) then
                    lo = a + width*grading**level; hi = a + width*grading**(level - 1)
                    if (level == levels) lo = a
                else
                    lo = b - width*grading**(level - 1); hi = b - width*grading**level
                    if (level == levels) hi = b
                end if
                cells = max(1, ceiling((hi - lo)/step))
                do c = 1, cells
                    half = 0.5_dp*(hi - lo)/cells
                    left = lo + (hi - lo)*(c - 1)/cells
                    if (store) then
                        table%ql(nl + 1:nl + low_order) = left + half*(1 + state%xl)
                        table%wl(nl + 1:nl + low_order) = half*state%wl
                        table%qh(nh + 1:nh + high_order) = left + half*(1 + state%xh)
                        table%wh(nh + 1:nh + high_order) = half*state%wh
                    end if
                    nl = nl + low_order; nh = nh + high_order
                end do
            end do
        end subroutine
    end subroutine

    !> phi(r)/phi_bare(r) = 1 - integral_0^inf 2 sin(qr) kappa2 G/(pi q (q^2+kappa2 G)) dq
    !> for a point test charge (kF = 1). tail: rigorous q-truncation bound;
    !> indicator: |rule20 - rule12| plus the G indicator propagated to first order.
    subroutine screened_ratio_from_table(state, table, r, ratio, tail, indicator)
        type(screening_state_t), intent(in) :: state
        type(screening_q_table_t), intent(in) :: table
        real(dp), intent(in) :: r(:)
        real(dp), intent(out) :: ratio(:), tail, indicator
        real(dp) :: low, high, term, sensitivity, dterm
        integer :: i, j
        indicator = 0.0_dp
        do i = 1, size(r)
            low = 0.0_dp; high = 0.0_dp; sensitivity = 0.0_dp
            do j = 1, size(table%ql)
                call screening_potential_term(table%ql(j), r(i), state%kappa2, &
                                              table%gl(j), term)
                low = low + table%wl(j)*term
            end do
            do j = 1, size(table%qh)
                call screening_potential_term(table%qh(j), r(i), state%kappa2, &
                                              table%gh(j), term)
                high = high + table%wh(j)*term
                call screening_potential_term_jvp(table%qh(j), r(i), state%kappa2, &
                    table%gh(j), 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, dterm)
                sensitivity = sensitivity + table%wh(j)*abs(dterm)
            end do
            ratio(i) = 1.0_dp - high
            ! First-order propagation of the G indicator (indicator, not bound).
            indicator = max(indicator, abs(high - low) + table%g_indicator*sensitivity)
        end do
        tail = table%tail
    end subroutine

    !> Convenience: build the table and invert in one call.
    subroutine screened_ratio(state, r, ratio, tail, indicator, qmax)
        type(screening_state_t), intent(in) :: state
        real(dp), intent(in) :: r(:)
        real(dp), intent(out) :: ratio(:), tail, indicator
        real(dp), intent(in), optional :: qmax
        type(screening_q_table_t) :: table
        real(dp) :: qm
        qm = 200.0_dp
        if (present(qmax)) qm = qmax
        call screening_q_table(state, qm, min(0.05_dp, 1.0_dp/maxval(r)), &
                               min(0.25_dp, 4.0_dp/maxval(r)), table)
        call screened_ratio_from_table(state, table, r, ratio, tail, indicator)
    end subroutine
end module kin6d_screening_hierarchy
