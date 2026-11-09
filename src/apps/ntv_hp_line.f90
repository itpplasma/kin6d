!> Thin-orbit NTV geometry of a general field line for the spectral-element
!! kernel (milca_dk_hp): |B| from the trigonometric interpolant of the line
!! (hp_field_t, which also fits the element rows to its orbit classes), the
!! precession from the line provider's drift profiles. The local equation
!! is that of milca_dk_ntv with the full drift weight,
!!   i (n/v) (omega_E + pre (drift (1 + xi^2) + drift2 xi^2)),
!! in the measure 1/(b k) of the weak form, with Boozer streaming
!! b.grad theta = kappa b (operator k = kappa, k_power = 1).
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
        class(ntv_line_t), allocatable :: line
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
        allocate (g%line, source=line)
    end subroutine hp_line_field

    subroutine line_point(g, th, xi, op, collide, drive_r, drive_e, density, phase, &
            precess, drive, torque)
        class(hp_line_field_t), intent(in) :: g
        real(dp), intent(in) :: th, xi
        type(hp_operator_t), intent(in) :: op
        real(dp), intent(out) :: collide, drive_r, drive_e, density, phase, precess, drive, &
            torque
        real(dp) :: bb, db, k, s(1), stream(1), mirror(1), drift(1), drift2(1), source(1), &
            measure(1)
        integer :: m
        bb = g%a(0)
        db = 0
        do m = 1, g%m
            bb = bb + g%a(m)*cos(m*th) + g%c(m)*sin(m*th)
            db = db + m*(g%c(m)*cos(m*th) - g%a(m)*sin(m*th))
        end do
        k = op%k*bb**op%k_power
        call dk_hp_point(bb, db, xi, k, op%nu, collide, drive_r, drive_e, density)
        call dk_hp_ntv(bb, xi, k, 0.0_dp, 0.0_dp, op%wsign, phase, precess, drive, torque)
        s(1) = th
        call g%line%sample(s, stream, mirror, drift, drift2, source, measure)
        ! The drift profiles per unit pre replace the circular cos(theta)/b law.
        precess = density*(drift(1)*(1 + op%wsign*xi*xi) + drift2(1)*xi*xi)
        torque = torque/bb**op%k_power
    end subroutine line_point

end module milca_ntv_hp_line
