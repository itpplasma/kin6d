module kin6d_response_calibration
    !! Direct spot readback. No singular powers or asymptotic transport.
    !! Fixed sample coordinates have a linear storage JVP/VJP. Changes in
    !! equilibrium, spot and sample coordinates need a separate derivative.
    !! All outputs are candidates; continuous enclosures are unavailable.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortnum_polynomial, only: lagrange_deriv_weights
    use kin6d_response_calibration_leaves, only: sheet_field, &
        distance_leaf => spot_distance
    implicit none
    private
    public :: spot_weights, spot_jump, spot_jump_vjp, equivalent_sheet_field
    public :: spot_distance
contains
    subroutine spot_distance(n, dq, spot, xs, status)
        integer, intent(in) :: n
        real(dp), intent(in) :: dq, spot
        real(dp), intent(out) :: xs
        integer, intent(out) :: status
        xs = 0.0_dp
        status = 1
        if (n == 0) return
        if (.not. ieee_is_finite(dq) .or. .not. ieee_is_finite(spot)) return
        if (abs(dq) <= tiny(1.0_dp) .or. spot <= 0.0_dp) return
        call distance_leaf(abs(real(n, dp)), abs(dq), spot, xs)
        if (.not. ieee_is_finite(xs)) return
        if (xs <= 0.0_dp) return
        status = 0
    end subroutine spot_distance

    subroutine spot_weights(psin, center, n, dq, spot, weights, status)
        real(dp), intent(in) :: psin(:, :), center, dq, spot
        integer, intent(in) :: n
        real(dp), intent(out) :: weights(size(psin, 1), 2)
        integer, intent(out) :: status
        real(dp) :: xs, at
        integer :: side, i, nn
        weights = 0.0_dp
        status = 1
        if (size(psin, 2) /= 2) return
        nn = size(psin, 1)
        if (nn < 3 .or. n == 0) return
        if (.not. all(ieee_is_finite(psin))) return
        if (.not. all(ieee_is_finite([center, dq, spot]))) return
        call spot_distance(n, dq, spot, xs, status)
        if (status /= 0) return
        status = 1
        do side = 1, 2
            at = center + merge(-xs, xs, side == 1)
            if (at < psin(1, side) .or. at > psin(nn, side)) return
            do i = 2, nn
                if (psin(i, side) <= psin(i - 1, side)) return
            end do
            ! Each stencil stays on its side of the sheet, including at a cut.
            if (side == 1) then
                if (psin(nn, side) >= center) return
            else
                if (psin(1, side) <= center) return
            end if
            call lagrange_deriv_weights(nn, at, psin(:, side), weights(:, side))
        end do
        status = 0
    end subroutine spot_weights

    pure function spot_jump(weights, flux) result(jump)
        real(dp), intent(in) :: weights(:, :)
        complex(dp), intent(in) :: flux(:, :)
        complex(dp) :: jump
        jump = sum(weights(:, 2)*flux(:, 2)) - sum(weights(:, 1)*flux(:, 1))
    end function spot_jump

    pure subroutine spot_jump_vjp(weights, seed, adjoint)
        real(dp), intent(in) :: weights(:, :)
        complex(dp), intent(in) :: seed
        complex(dp), intent(out) :: adjoint(size(weights, 1), 2)
        adjoint(:, 1) = -weights(:, 1)*seed
        adjoint(:, 2) = weights(:, 2)*seed
    end subroutine spot_jump_vjp

    subroutine equivalent_sheet_field(jump, unit_jump, unit_flux, area, b, status)
        complex(dp), intent(in) :: jump, unit_jump, unit_flux
        real(dp), intent(in) :: area
        complex(dp), intent(out) :: b
        integer, intent(out) :: status
        real(dp) :: br, bi
        status = 1
        b = (0.0_dp, 0.0_dp)
        if (.not. all(ieee_is_finite([real(jump), aimag(jump), real(unit_jump), &
            aimag(unit_jump), real(unit_flux), aimag(unit_flux), area]))) return
        if (area <= 0.0_dp .or. abs(unit_jump) <= tiny(1.0_dp)) return
        if (abs(unit_flux) <= tiny(1.0_dp)) return
        call sheet_field(real(jump), aimag(jump), real(unit_jump), aimag(unit_jump), &
            real(unit_flux), aimag(unit_flux), area, br, bi)
        b = cmplx(br, bi, dp)
        if (.not. all(ieee_is_finite([br, bi]))) return
        status = 0
    end subroutine equivalent_sheet_field
end module kin6d_response_calibration
