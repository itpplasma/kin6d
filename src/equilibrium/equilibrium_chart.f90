module kin6d_equilibrium_chart
    ! Chart of the common equilibrium state (docs/equilibrium-core.md, D1):
    ! X = H(q) from chart coordinates q into Cartesian space, with the jet H,
    ! dH, d2H, d3H. first_jet returns only H and dH, with a full-jet
    ! fallback for charts that do not specialize it. A curve-following frame
    ! extends ideal_chart_t.
    ! This is the chart of the N1 map lane (ideal_map), lifted unchanged into
    ! the core; ideal_map uses it.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, &
        ieee_quiet_nan
    use kin6d_equilibrium_core_leaves, only: cylindrical_chart_position, &
        cylindrical_chart_first, cylindrical_chart_second, cylindrical_chart_third
    use kin6d_curve_chart_leaves, only: curve_fourier_mode, curve_normalize, &
        curve_cross, curve_compose, curve_geometry, curve_norm, curve_rotate, &
        curve_mode_bounds, curve_tube_bounds, curve_admission_bounds, &
        curve_fourier_shift
    implicit none
    private
    public :: ideal_chart_t, cylindrical_chart_t, curve_chart_t

    ! jet returns H, dH(a,k) = dH_a/dq_k and the second and third
    ! q-derivatives d2(a,k,l), d3(a,k,l,m). A consumer fixes the orientation
    ! sigma = sign(det dH) at its reference.
    type, abstract :: ideal_chart_t
    contains
        procedure(chart_jet), deferred :: jet
        procedure :: first_jet => chart_first_jet
    end type
    abstract interface
        pure subroutine chart_jet(chart, q, x, d1, d2, d3)
            import :: ideal_chart_t, dp
            class(ideal_chart_t), intent(in) :: chart
            real(dp), intent(in) :: q(3)
            real(dp), intent(out) :: x(3), d1(3, 3), d2(3, 3, 3), d3(3, 3, 3, 3)
        end subroutine
    end interface

    ! q = (R, Z, phi); X = (R cos phi, R sin phi, Z). (R, Z, phi) is left-handed.
    type, extends(ideal_chart_t) :: cylindrical_chart_t
    contains
        procedure :: jet => cylindrical_jet
        procedure :: first_jet => cylindrical_first_jet
    end type
    ! Orthonormal projected-guide G-frame (not Frenet or parallel transport).
    ! All coefficients use full-turn Cartesian Fourier modes, not nfp*m.
    type, extends(ideal_chart_t) :: curve_chart_t
        private
        real(dp), allocatable :: axis_c(:, :), axis_s(:, :)
        real(dp), allocatable :: guide_c(:, :), guide_s(:, :)
        integer :: nfp = 1
        real(dp) :: speed_lower = 0, acceleration_upper = 0, speed_upper = 0
    contains
        procedure :: init => curve_init
        procedure :: jet => curve_jet
        procedure :: tube_radii => curve_tube_radii
        procedure :: field_periods => curve_field_periods
    end type
contains
    pure subroutine chart_first_jet(chart, q, x, d1)
        class(ideal_chart_t), intent(in) :: chart
        real(dp), intent(in) :: q(3)
        real(dp), intent(out) :: x(3), d1(3, 3)
        real(dp) :: d2(3, 3, 3), d3(3, 3, 3, 3)
        call chart%jet(q, x, d1, d2, d3)
    end subroutine

    pure subroutine cylindrical_first_jet(chart, q, x, d1)
        class(cylindrical_chart_t), intent(in) :: chart
        real(dp), intent(in) :: q(3)
        real(dp), intent(out) :: x(3), d1(3, 3)
        select type (chart)
        type is (cylindrical_chart_t)
            call cylindrical_chart_position(q(1), q(2), q(3), x(1), x(2), x(3))
            call cylindrical_chart_first(q(1), q(2), q(3), d1(1, 1), d1(2, 1), &
                d1(3, 1), d1(1, 2), d1(2, 2), d1(3, 2), d1(1, 3), d1(2, 3), d1(3, 3))
        class default
            call chart_first_jet(chart, q, x, d1)
        end select
    end subroutine

    pure integer function curve_field_periods(chart) result(nfp)
        class(curve_chart_t), intent(in) :: chart
        nfp = 0
        if (allocated(chart%axis_c)) nfp = chart%nfp
    end function

    subroutine curve_init(chart, axis_c, axis_s, nfp, status, guide_c, guide_s, &
                          samples)
        ! status: 0 admitted, 1 invalid input/symmetry, 2 unresolved frame.
        ! A failed initialization leaves the chart uninitialized (jet -> NaNs).
        class(curve_chart_t), intent(inout) :: chart
        real(dp), intent(in) :: axis_c(:, 0:), axis_s(:, 0:)
        integer, intent(in) :: nfp
        integer, intent(out) :: status
        real(dp), intent(in), optional :: guide_c(:, 0:), guide_s(:, 0:)
        integer, intent(in), optional :: samples
        real(dp) :: a(3, 0:4), g(3, 0:4), v(3, 0:3), u(3, 0:3), z(3, 0:3)
        real(dp) :: ab(0:2), gb(0:2), speed, gap, min_speed, min_gap, step, gap_lower
        integer :: i, count
        if (allocated(chart%axis_c)) then
            deallocate (chart%axis_c, chart%axis_s, chart%guide_c, chart%guide_s)
        end if
        status = 1
        if (nfp < 1) return
        if (.not. valid_series(axis_c, axis_s, nfp)) return
        if (present(guide_c) .neqv. present(guide_s)) return
        if (present(guide_c)) then
            if (.not. valid_series(guide_c, guide_s, nfp)) return
            chart%guide_c = guide_c; chart%guide_s = guide_s
        else
            allocate (chart%guide_c(3, 0:0), chart%guide_s(3, 0:0))
            chart%guide_c = 0; chart%guide_s = 0
            chart%guide_c(3, 0) = 1
        end if
        chart%axis_c = axis_c; chart%axis_s = axis_s; chart%nfp = nfp
        count = max(1024, 16*size(axis_c, 2), 16*size(chart%guide_c, 2))
        if (present(samples)) count = max(count, samples)
        step = 2*acos(-1.0_dp)/count
        call series_bounds(chart%axis_c, chart%axis_s, ab)
        call series_bounds(chart%guide_c, chart%guide_s, gb)
        min_speed = huge(1.0_dp); min_gap = huge(1.0_dp)
        do i = 0, count - 1
            call series(chart%axis_c, chart%axis_s, i*step, a)
            call series(chart%guide_c, chart%guide_s, i*step, g)
            v = 0; u = 0
            v(:, 0) = a(:, 1); u(:, 0) = g(:, 0)
            call curve_cross(v, u, z)
            call curve_norm(a(1, 1), a(2, 1), a(3, 1), speed)
            call curve_norm(z(1, 0), z(2, 0), z(3, 0), gap)
            min_speed = min(min_speed, speed); min_gap = min(min_gap, gap)
        end do
        ! Lipschitz padding includes points between samples. Floating-point
        ! bounds are admission diagnostics, not outward-rounded certificates.
        call curve_admission_bounds(min_speed, min_gap, ab(2), ab(1), gb(0), gb(1), &
                                    step, chart%speed_lower, gap_lower)
        chart%acceleration_upper = ab(2); chart%speed_upper = ab(1)
        status = 2
        if (chart%speed_lower > 64*epsilon(step)*max(1.0_dp, ab(1))) then
            if (gap_lower > 64*epsilon(step)*max(1.0_dp, ab(1)*gb(0))) status = 0
        end if
        if (status /= 0) &
            deallocate (chart%axis_c, chart%axis_s, chart%guide_c, chart%guide_s)
    end subroutine

    logical function valid_series(c, s, nfp) result(ok)
        real(dp), intent(in) :: c(:, 0:), s(:, 0:)
        integer, intent(in) :: nfp
        real(dp) :: angle, phase, rc(3), rs(3), cs(3), ss(3), tol
        integer :: m
        ok = .false.
        if (size(c, 1) /= 3 .or. size(s, 1) /= 3) return
        if (any(shape(c) /= shape(s))) return
        if (size(c, 2) < 1) return
        if (.not. all(ieee_is_finite(c))) return
        if (.not. all(ieee_is_finite(s))) return
        if (any(abs(s(:, 0)) > 0)) return
        angle = 2*acos(-1.0_dp)/nfp
        tol = 128*epsilon(angle)*max(1.0_dp, maxval(abs(c)), maxval(abs(s)))
        do m = 0, ubound(c, 2)
            phase = mod(m, nfp)*angle
            call curve_rotate(c(1, m), c(2, m), c(3, m), angle, rc(1), rc(2), rc(3))
            call curve_rotate(s(1, m), s(2, m), s(3, m), angle, rs(1), rs(2), rs(3))
            call curve_fourier_shift(c(:, m), s(:, m), phase, cs, ss)
            if (maxval(abs(rc - cs)) > tol) return
            if (maxval(abs(rs - ss)) > tol) return
        end do
        ok = .true.
    end function

    pure subroutine series(c, s, t, a)
        real(dp), intent(in) :: c(:, 0:), s(:, 0:), t
        real(dp), intent(out) :: a(3, 0:4)
        real(dp) :: f(0:4)
        integer :: i, m
        a = 0
        do m = 0, ubound(c, 2)
            do i = 1, 3
                call curve_fourier_mode(c(i, m), s(i, m), real(m, dp), t, &
                                        f(0), f(1), f(2), f(3), f(4))
                a(i, :) = a(i, :) + f
            end do
        end do
    end subroutine

    subroutine series_bounds(c, s, bounds)
        real(dp), intent(in) :: c(:, 0:), s(:, 0:)
        real(dp), intent(out) :: bounds(0:2)
        real(dp) :: cn, sn, b(0:2)
        integer :: m
        bounds = 0
        do m = 0, ubound(c, 2)
            call curve_norm(c(1, m), c(2, m), c(3, m), cn)
            call curve_norm(s(1, m), s(2, m), s(3, m), sn)
            call curve_mode_bounds(cn, sn, real(m, dp), b(0), b(1), b(2))
            bounds = bounds + b
        end do
    end subroutine

    pure subroutine curve_jet(chart, q, x, d1, d2, d3)
        class(curve_chart_t), intent(in) :: chart
        real(dp), intent(in) :: q(3)
        real(dp), intent(out) :: x(3), d1(3, 3), d2(3, 3, 3), d3(3, 3, 3, 3)
        real(dp) :: a(3, 0:4), g(3, 0:4), t(3, 0:3), c(3, 0:3), n(3, 0:3), b(3, 0:3)
        real(dp) :: h13, h23, h33, h133, h233, h333, bad
        integer :: i, k
        if (.not. allocated(chart%axis_c)) then
            bad = ieee_value(0.0_dp, ieee_quiet_nan)
            x = bad; d1 = bad; d2 = bad; d3 = bad
            return
        end if
        call series(chart%axis_c, chart%axis_s, q(3), a)
        call series(chart%guide_c, chart%guide_s, q(3), g)
        call curve_normalize(a(:, 1:4), t)
        call curve_cross(t, g(:, 0:3), c)
        call curve_normalize(c, b)
        call curve_cross(b, t, n)
        d2 = 0; d3 = 0
        do i = 1, 3
            call curve_compose(a(i, 0), a(i, 1), a(i, 2), a(i, 3), &
               n(i, 0), n(i, 1), n(i, 2), n(i, 3), b(i, 0), b(i, 1), b(i, 2), b(i, 3), &
                               q(1), q(2), x(i), d1(i, 1), d1(i, 2), d1(i, 3), &
                               h13, h23, h33, h133, h233, h333)
            d2(i, 1, 3) = h13; d2(i, 3, 1) = h13
            d2(i, 2, 3) = h23; d2(i, 3, 2) = h23; d2(i, 3, 3) = h33
            do k = 1, 2
                if (k == 1) then
                    d3(i, k, 3, 3) = h133; d3(i, 3, k, 3) = h133; d3(i, 3, 3, k) = h133
                else
                    d3(i, k, 3, 3) = h233; d3(i, 3, k, 3) = h233; d3(i, 3, 3, k) = h233
                end if
            end do
            d3(i, 3, 3, 3) = h333
        end do
    end subroutine

    subroutine curve_tube_radii(chart, samples, curvature_radius, tube_radius, status)
        ! curvature_radius: sampled 1/max kappa, local ceiling only.
        ! tube_radius: conservative nonoverlap lower bound in real arithmetic;
        ! rounding is NOT enclosed. status 0: numerical diagnostic available.
        class(curve_chart_t), intent(in) :: chart
        integer, intent(in) :: samples
        real(dp), intent(out) :: curvature_radius, tube_radius
        integer, intent(out) :: status
        real(dp), allocatable :: points(:, :)
        real(dp) :: a(3, 0:4), step, delta, kmax, speed, kappa, far, dist, local, apart
        integer :: i, j, separation
        status = 1; curvature_radius = 0; tube_radius = 0
        if (.not. allocated(chart%axis_c)) return
        if (samples < 16) return
        step = 2*acos(-1.0_dp)/samples
        delta = min(acos(-1.0_dp), chart%speed_lower/chart%acceleration_upper)
        if (step >= delta/2) return
        allocate (points(3, 0:samples - 1))
        kmax = 0
        do i = 0, samples - 1
            call series(chart%axis_c, chart%axis_s, i*step, a)
            points(:, i) = a(:, 0)
            call curve_geometry(a(1, 1), a(2, 1), a(3, 1), a(1, 2), a(2, 2), &
                                a(3, 2), speed, kappa)
            kmax = max(kmax, kappa)
        end do
        far = huge(1.0_dp)
        do j = 0, samples - 1
            do i = 0, j - 1
                separation = min(j - i, samples - (j - i))
                if (separation*step < delta - step) cycle
                call curve_norm(points(1, i) - points(1, j), &
                         points(2, i) - points(2, j), points(3, i) - points(3, j), dist)
                far = min(far, dist)
            end do
        end do
        curvature_radius = 1/kmax
        call curve_tube_bounds(chart%speed_lower, chart%acceleration_upper, &
                               chart%speed_upper, step, far, local, apart)
        tube_radius = max(0.0_dp, min(local, apart))
        status = 0
    end subroutine

    pure subroutine cylindrical_jet(chart, q, x, d1, d2, d3)
        class(cylindrical_chart_t), intent(in) :: chart
        real(dp), intent(in) :: q(3)
        real(dp), intent(out) :: x(3), d1(3, 3), d2(3, 3, 3), d3(3, 3, 3, 3)
        call cylindrical_chart_position(q(1), q(2), q(3), x(1), x(2), x(3))
        call cylindrical_chart_first(q(1), q(2), q(3), d1(1, 1), d1(2, 1), &
            d1(3, 1), d1(1, 2), d1(2, 2), d1(3, 2), d1(1, 3), d1(2, 3), d1(3, 3))
        call jet2(d2)
        call jet3(d3)
    contains
        pure subroutine jet2(t)
            real(dp), intent(out) :: t(27)
            call cylindrical_chart_second(q(1), q(2), q(3), t(1), t(2), t(3), t(4), &
                t(5), t(6), t(7), t(8), t(9), t(10), t(11), t(12), t(13), t(14), &
                t(15), t(16), t(17), t(18), t(19), t(20), t(21), t(22), t(23), &
                t(24), t(25), t(26), t(27))
        end subroutine
        pure subroutine jet3(t)
            real(dp), intent(out) :: t(81)
            call cylindrical_chart_third(q(1), q(2), q(3), t(1), t(2), t(3), t(4), &
                t(5), t(6), t(7), t(8), t(9), t(10), t(11), t(12), t(13), t(14), &
                t(15), t(16), t(17), t(18), t(19), t(20), t(21), t(22), t(23), &
                t(24), t(25), t(26), t(27), t(28), t(29), t(30), t(31), t(32), &
                t(33), t(34), t(35), t(36), t(37), t(38), t(39), t(40), t(41), &
                t(42), t(43), t(44), t(45), t(46), t(47), t(48), t(49), t(50), &
                t(51), t(52), t(53), t(54), t(55), t(56), t(57), t(58), t(59), &
                t(60), t(61), t(62), t(63), t(64), t(65), t(66), t(67), t(68), &
                t(69), t(70), t(71), t(72), t(73), t(74), t(75), t(76), t(77), &
                t(78), t(79), t(80), t(81))
        end subroutine
    end subroutine
end module kin6d_equilibrium_chart
