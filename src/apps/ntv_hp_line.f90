!> Thin-orbit NTV geometry of a general field line for the spectral-element
!! kernel (milca_dk_hp): |B| from the trigonometric interpolant of the line
!! (hp_field_t, which also fits the element rows to its orbit classes), the
!! precession from the line provider's drift profiles. The local equation
!! is that of milca_dk_ntv with the full drift weight,
!!   i (n/v) (omega_E + pre (drift (1 + xi^2) + drift2 xi^2)),
!! in the measure 1/(b k) of the weak form, with Boozer streaming
!! b.grad theta = kappa b (operator k = kappa, k_power = 1). The drift
!! profiles enter through their trigonometric interpolants on the same
!! samples as b.
module milca_ntv_hp_line
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use milca_dk_hp, only: hp_field_t, hp_field, hp_operator_t, HP_OK
    use milca_dk_hp_leaves, only: dk_hp_point, dk_hp_ntv
    use milca_dk_ntv, only: ntv_line_t
    implicit none
    private
    public :: hp_line_field_t, hp_line_field

    real(dp), parameter :: pi = 4*atan(1.0_dp)

    type, extends(hp_field_t) :: hp_line_field_t
        ! Trigonometric interpolants of the drift profiles on the samples of b:
        ! a quadrature point then costs no Boozer sum.
        real(dp), allocatable :: da(:, :), dc(:, :)
    contains
        procedure :: point => line_point
    end type

contains

    !> Geometry of line on n_sample (odd) points of b = 2/(3 source); kappa
    !> = b.grad theta/b of the line.
    subroutine hp_line_field(line, n_sample, g, kappa, status)
        class(ntv_line_t), intent(in) :: line
        integer, intent(in) :: n_sample
        type(hp_line_field_t), intent(out) :: g
        real(dp), intent(out) :: kappa
        integer, intent(out) :: status
        real(dp), dimension(n_sample) :: theta, stream, mirror, drift, drift2, source, &
            measure
        integer :: j
        theta = [(2*pi*(j - 1)/n_sample, j=1, n_sample)]
        call line%sample(theta, stream, mirror, drift, drift2, source, measure)
        kappa = stream(1)*3*source(1)/2
        call hp_field(2/(3*source), kappa, g%hp_field_t, status)
        if (status /= HP_OK) return
        allocate (g%da(0:g%m, 2), g%dc(g%m, 2))
        g%da(0, 1) = sum(drift)/n_sample
        g%da(0, 2) = sum(drift2)/n_sample
        do j = 1, g%m
            g%da(j, 1) = 2*sum(drift*cos(j*theta))/n_sample
            g%dc(j, 1) = 2*sum(drift*sin(j*theta))/n_sample
            g%da(j, 2) = 2*sum(drift2*cos(j*theta))/n_sample
            g%dc(j, 2) = 2*sum(drift2*sin(j*theta))/n_sample
        end do
    end subroutine hp_line_field

    subroutine line_point(g, th, xi, op, collide, drive_r, drive_e, density, phase, &
            precess, drive, torque)
        class(hp_line_field_t), intent(in) :: g
        real(dp), intent(in) :: th, xi
        type(hp_operator_t), intent(in) :: op
        real(dp), intent(out) :: collide, drive_r, drive_e, density, phase, precess, drive, &
            torque
        real(dp) :: bb, db, k, drift(2), cm, sm
        complex(dp) :: e, z
        integer :: m
        bb = g%a(0)
        db = 0
        drift = g%da(0, :)
        ! cos(m th), sin(m th) by angle addition from one sincos.
        e = cmplx(cos(th), sin(th), dp)
        z = 1
        do m = 1, g%m
            z = z*e
            cm = real(z, dp)
            sm = aimag(z)
            bb = bb + g%a(m)*cm + g%c(m)*sm
            db = db + m*(g%c(m)*cm - g%a(m)*sm)
            drift = drift + g%da(m, :)*cm + g%dc(m, :)*sm
        end do
        k = op%k*bb**op%k_power
        call dk_hp_point(bb, db, xi, k, op%nu, collide, drive_r, drive_e, density)
        call dk_hp_ntv(bb, xi, k, 0.0_dp, 0.0_dp, op%wsign, phase, precess, drive, torque)
        ! The drift profiles per unit pre replace the circular cos(theta)/b law.
        precess = density*(drift(1)*(1 + op%wsign*xi*xi) + drift(2)*xi*xi)
        torque = torque/bb**op%k_power
    end subroutine line_point

end module milca_ntv_hp_line
