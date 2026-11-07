!> Single-ion force balance in the Gaussian units of KAMEL's profile files.
module kin6d_dk_force_balance
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use kin6d_dk_energy, only: ENERGY_OK, ENERGY_BAD_INPUT, ENERGY_SINGULAR
    use kin6d_dk_energy_leaves, only: dk_energy_force_balance, dk_energy_force_balance_jvp
    implicit none
    private
    public :: energy_balance_input_t, energy_balance_t, energy_force_balance_cgs

    !> T in erg, n in cm^-3, B in G, charge in statC, r and r0 in cm,
    !> vz in cm/s; nr and tr differentiate n and T with respect to r.
    !> k is the rotation coefficient; q is signed and bz is the signed toroidal field.
    type :: energy_balance_input_t
        real(dp) :: k = 0, n = 0, t = 0, nr = 0, tr = 0, b = 0, bz = 0
        real(dp) :: r = 0, vz = 0, q = 0, r0 = 0, charge = 0
    end type

    type :: energy_balance_t
        real(dp) :: er = 0, vth = 0 ! statV/cm, cm/s
        real(dp) :: der = 0, dvth = 0
        character(len=11) :: certification = 'unavailable'
    end type
contains
    subroutine energy_force_balance_cgs(p, value, status, direction)
        type(energy_balance_input_t), intent(in) :: p
        type(energy_balance_t), intent(out) :: value
        integer, intent(out) :: status
        type(energy_balance_input_t), intent(in), optional :: direction
        status = ENERGY_BAD_INPUT
        if (.not. all(ieee_is_finite(parameters(p)))) return
        if (min(p%n, p%t, p%b, p%r0) <= 0 .or. p%r < 0) return
        if (min(abs(p%q), abs(p%charge)) <= tiny(1.0_dp)) return
        if (present(direction)) then
            if (.not. all(ieee_is_finite(parameters(direction)))) return
        end if
        status = ENERGY_SINGULAR
        call dk_energy_force_balance(p%k, p%n, p%t, p%nr, p%tr, p%b, p%bz, &
            p%r, p%vz, p%q, p%r0, p%charge, value%er, value%vth)
        if (present(direction)) then
            call dk_energy_force_balance_jvp(p%k, p%n, p%t, p%nr, p%tr, p%b, p%bz, &
                p%r, p%vz, p%q, p%r0, p%charge, direction%k, direction%n, direction%t, &
                direction%nr, direction%tr, direction%b, direction%bz, direction%r, &
                direction%vz, direction%q, direction%r0, direction%charge, value%der, value%dvth)
        end if
        if (.not. all(ieee_is_finite([value%er, value%vth, value%der, value%dvth]))) return
        status = ENERGY_OK
    end subroutine energy_force_balance_cgs

    pure function parameters(p) result(values)
        type(energy_balance_input_t), intent(in) :: p
        real(dp) :: values(12)
        values = [p%k, p%n, p%t, p%nr, p%tr, p%b, p%bz, p%r, p%vz, p%q, p%r0, p%charge]
    end function parameters
end module kin6d_dk_force_balance
