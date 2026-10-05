program test_streaming_moments
    use, intrinsic :: iso_fortran_env, only: dp => real64, qp => real128
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use fortnum_interval, only: interval_t, interval
    use kin6d_streaming_moments, only: streaming_moments_t, &
        initialize_streaming_moments, advance_streaming_moments, &
        sample_streaming_moments, sample_streaming_moments_box
    use kin6d_streaming_moments_leaves, only: streaming_moments_current_cell, &
        streaming_moments_current_cell_interval, streaming_moments_current_cell_jvp, &
        streaming_moments_current_cell_vjp
    implicit none
    type(streaming_moments_t) :: packet
    real(dp) :: v(65), a(0:1, 65), b(0:1, 65), t, x, n, j, u, nan
    type(interval_t) :: nb, jb, ub
    real(qp) :: nq, jq, uq
    integer :: i, frame, point
    logical :: ok

    call cell_oracle()
    do i = 1, 65
        v(i) = -1.0_dp + real(i - 1, dp)/32
        a(0, i) = 1.0_dp + v(i)/2
        a(1, i) = a(0, i)/2
    end do
    b = 0.0_dp
    call initialize_streaming_moments(packet, v, a, b, ok)
    call require(ok, 'moment initialization')
    do frame = 0, 3
        t = real(frame, dp)*2
        call advance_streaming_moments(packet, t, ok)
        call require(ok, 'moment evolution')
        do point = 0, 20
            x = real(point, dp)/4
            call analytic_moments(real(x, qp), real(t, qp), nq, jq)
            uq = jq/nq
            call sample_streaming_moments(packet, x, n, j, u, nb, jb, ub, ok)
            call require(ok, 'positive density point')
            call contains(nb, nq, 'independent density integral')
            call contains(jb, jq, 'independent current integral')
            call contains(ub, uq, 'independent bulk quotient')
            call sample_streaming_moments_box(packet, interval(x, x + 0.25_dp), &
                nb, jb, ub, ok)
            call require(ok, 'whole-cell moment box')
            do i = 0, 8
                call analytic_moments(real(x, qp) + real(i, qp)/32, &
                    real(t, qp), nq, jq)
                call contains(nb, nq, 'whole-cell density')
                call contains(jb, jq, 'whole-cell current')
                call contains(ub, jq/nq, 'whole-cell bulk')
            end do
        end do
    end do
    nan = ieee_value(0.0_dp, ieee_quiet_nan)
    call advance_streaming_moments(packet, nan, ok)
    call require(.not. ok, 'nonfinite time rejected')
    call sample_streaming_moments_box(packet, interval(0.0_dp), nb, jb, ub, ok)
    call require(.not. ok, 'failed preparation invalidates old moment bands')
    a(0, :) = -1.0_dp; a(1, :) = 0.0_dp
    call initialize_streaming_moments(packet, v, a, b, ok)
    call require(ok, 'signed source still valid for distribution')
    call sample_streaming_moments_box(packet, interval(0.0_dp), nb, jb, ub, ok)
    call require(.not. ok, 'bulk requires strictly positive density')
    print '(a)', 'Streaming mass/current/bulk and generated derivative oracles pass'
contains
    subroutine analytic_moments(x, t, n, j)
        real(qp), intent(in) :: x, t
        real(qp), intent(out) :: n, j
        real(qp) :: sinc, odd, even

        ! Independent exact integrals of (1+.5 cos(x-vt))*(1+.5v).
        if (abs(t) < 1.0e-30_qp) then
            sinc = 1.0_qp; odd = 0.0_qp; even = 1.0_qp/3
        else
            sinc = sin(t)/t
            odd = (sin(t) - t*cos(t))/t**2
            even = ((t**2 - 2)*sin(t) + 2*t*cos(t))/t**3
        end if
        n = 1 + sinc*cos(x)/2 + odd*sin(x)/4
        j = 1.0_qp/6 + odd*sin(x)/2 + even*cos(x)/4
    end subroutine

    subroutine cell_oracle()
        real(dp) :: args(5), direction(5), value, tangent, gradient(5)
        complex(qp) :: z(5), vm, am, exact
        real(qp), parameter :: epsilon = 1.0e-25_qp
        type(interval_t) :: box
        integer :: i

        args = [-0.75_dp, 0.25_dp, 1.25_dp, -0.5_dp, 0.375_dp]
        direction = [0.0625_dp, -0.125_dp, 0.25_dp, 0.5_dp, -0.375_dp]
        call streaming_moments_current_cell(args(1), args(2), args(3), &
            args(4), args(5), value)
        call streaming_moments_current_cell_interval(interval(args(1)), &
            interval(args(2)), interval(args(3)), interval(args(4)), &
            interval(args(5)), box)
        do i = 1, 5
            z(i) = cmplx(real(args(i), qp), epsilon*real(direction(i), qp), qp)
        end do
        vm = (z(1) + z(2))/2; am = (z(3) + z(4))/2
        ! Simpson integrates v times its affine reconstruction exactly.
        exact = z(5) + (z(2) - z(1))*(z(1)*z(3) + 4*vm*am + z(2)*z(4))/12
        call contains(box, real(exact, qp), 'independent quadratic Simpson integral')
        call require(abs(real(value, qp) - real(exact, qp)) < 1.0e-14_qp, 'cell center')
        call streaming_moments_current_cell_jvp(args(1), args(2), args(3), args(4), &
            args(5), direction(1), direction(2), direction(3), direction(4), &
            direction(5), tangent)
        call require(abs(real(tangent, qp) - aimag(exact)/epsilon) < 1.0e-14_qp, &
            'independent complex-step JVP')
        call streaming_moments_current_cell_vjp(args(1), args(2), args(3), args(4), &
            args(5), 0.75_dp, gradient(1), gradient(2), gradient(3), &
            gradient(4), gradient(5))
        call require(abs(dot_product(gradient, direction) - 0.75_dp*tangent) < &
            1.0e-14_dp, 'Euclidean transpose pairing')
    end subroutine

    subroutine contains(box, value, label)
        type(interval_t), intent(in) :: box
        real(qp), intent(in) :: value
        character(*), intent(in) :: label
        call require(real(box%lo, qp) <= value .and. value <= real(box%hi, qp), label)
    end subroutine

    subroutine require(condition, label)
        logical, intent(in) :: condition
        character(*), intent(in) :: label
        if (condition) return
        print '(a)', label
        error stop 1
    end subroutine
end program test_streaming_moments
