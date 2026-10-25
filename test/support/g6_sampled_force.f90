module g6_sampled_force
    ! External physical rows only: no KIN state import, normalization or gate.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, &
        ieee_quiet_nan
    implicit none
    private
    public :: g6_force_accumulator_t, g6_force_vector, g6_force_add, &
        g6_force_result
    type :: g6_force_accumulator_t
        private
        integer :: samples = 0
        real(dp) :: scale = 0, ssq = 1
        logical :: failed = .false.
    end type
contains
    pure subroutine g6_force_vector(b, mu0j, mu0gradp, force, status)
        real(dp), intent(in) :: b(3), mu0j(3), mu0gradp(3)
        real(dp), intent(out) :: force(3)
        integer, intent(out) :: status
        real(dp) :: product1, product2, cross, candidate(3)
        integer :: i, j, k
        logical :: ok
        force = ieee_value(0.0_dp, ieee_quiet_nan)
        status = 1
        if (.not. all(ieee_is_finite(b))) return
        if (.not. all(ieee_is_finite(mu0j))) return
        if (.not. all(ieee_is_finite(mu0gradp))) return
        do i = 1, 3
            j = modulo(i, 3) + 1
            k = modulo(i + 1, 3) + 1
            call finite_product(mu0j(j), b(k), product1, ok)
            if (.not. ok) return
            call finite_product(mu0j(k), b(j), product2, ok)
            if (.not. ok) return
            call finite_difference(product1, product2, cross, ok)
            if (.not. ok) return
            call finite_difference(cross, mu0gradp(i), candidate(i), ok)
            if (.not. ok) return
        end do
        force = candidate
        status = 0
    end subroutine

    pure subroutine g6_force_add(acc, b, mu0j, mu0gradp, status)
        type(g6_force_accumulator_t), intent(inout) :: acc
        real(dp), intent(in) :: b(3), mu0j(3), mu0gradp(3)
        integer, intent(out) :: status
        real(dp) :: force(3), a
        integer :: i
        status = 1
        if (acc%failed) return
        call g6_force_vector(b, mu0j, mu0gradp, force, status)
        if (status /= 0 .or. acc%samples == huge(acc%samples)) then
            acc%failed = .true.
            status = 1
            return
        end if
        do i = 1, 3
            a = abs(force(i))
            if (a == 0) cycle
            if (acc%scale < a) then
                acc%ssq = 1 + acc%ssq*(acc%scale/a)**2
                acc%scale = a
            else
                acc%ssq = acc%ssq + (a/acc%scale)**2
            end if
        end do
        acc%samples = acc%samples + 1
        status = 0
    end subroutine

    pure subroutine g6_force_result(acc, rms, status)
        type(g6_force_accumulator_t), intent(in) :: acc
        real(dp), intent(out) :: rms
        integer, intent(out) :: status
        real(dp) :: factor
        rms = ieee_value(0.0_dp, ieee_quiet_nan)
        status = 1
        if (acc%failed .or. acc%samples == 0) return
        if (acc%scale == 0) then
            rms = 0
            status = 0
            return
        end if
        ! Divide by the number of physical rows, not the component count.
        factor = sqrt(acc%ssq/real(acc%samples, dp))
        if (factor > 1) then
            if (acc%scale > huge(rms)/factor) return
        end if
        rms = acc%scale*factor
        if (.not. ieee_is_finite(rms)) then
            rms = ieee_value(0.0_dp, ieee_quiet_nan)
            return
        end if
        status = 0
    end subroutine

    pure subroutine finite_product(a, b, product, ok)
        real(dp), intent(in) :: a, b
        real(dp), intent(out) :: product
        logical, intent(out) :: ok
        product = 0
        ok = .false.
        if (abs(a) > 1) then
            if (abs(b) > huge(a)/abs(a)) return
        end if
        product = a*b
        ok = ieee_is_finite(product)
    end subroutine

    pure subroutine finite_difference(a, b, delta, ok)
        real(dp), intent(in) :: a, b
        real(dp), intent(out) :: delta
        logical, intent(out) :: ok
        delta = 0
        ok = .false.
        if (b < 0 .and. a > 0) then
            if (a > huge(a) + b) return
        else if (b > 0 .and. a < 0) then
            if (a < -huge(a) + b) return
        end if
        delta = a - b
        ok = ieee_is_finite(delta)
    end subroutine
end module
