module kin6d_equilibrium_chart
    ! Chart of the common equilibrium state (docs/equilibrium-core.md, D1):
    ! X = H(q) from chart coordinates q into Cartesian space, with the jet H,
    ! dH, d2H, d3H. Lanes see the chart only through `jet`; a curve-following
    ! frame is a further extension of ideal_chart_t and changes nothing else.
    ! This is the chart of the N1 map lane (ideal_map), lifted unchanged into
    ! the core; ideal_map uses it.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_equilibrium_core_leaves, only: cylindrical_chart_position, &
        cylindrical_chart_first, cylindrical_chart_second, cylindrical_chart_third
    implicit none
    private
    public :: ideal_chart_t, cylindrical_chart_t

    ! jet returns H, dH(a,k) = dH_a/dq_k and the second and third
    ! q-derivatives d2(a,k,l), d3(a,k,l,m). A consumer fixes the orientation
    ! sigma = sign(det dH) at its reference.
    type, abstract :: ideal_chart_t
    contains
        procedure(chart_jet), deferred :: jet
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
    end type
contains
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
