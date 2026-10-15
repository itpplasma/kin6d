!> C1: two electrons in a periodic cube with neutralizing background, the
!> classical twin of the Q0 quantum pair quench. Units, species and cell edge
!> come from kin6d_species (periodic_pair_units, periodic_pair_length), the
!> same declarations the quantum reference uses; the Ewald parameters are
!> selected from a requested truncation tolerance and the dynamics are the
!> shared kin6d_nbody KDK path on the full two-particle state.
module kin6d_classical_pair
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use kin6d_species, only: units_t, periodic_pair_units, periodic_pair_length, &
                             electron
    use kin6d_ewald, only: ewald_parameters_t, ewald_select
    use kin6d_nbody, only: nbody_system_t, nbody_periodic
    implicit none
    private
    public :: classical_pair_model_t, initialize_classical_pair
    public :: classical_pair_sample, classical_pair_wrap, classical_pair_state
    public :: classical_pair_bin_center

    type :: classical_pair_model_t
        type(units_t) :: units
        real(dp) :: length = periodic_pair_length
        !> Ewald truncation tolerance relative to |g|/L (potential) and
        !> |g|/L^2 (force components).
        real(dp) :: tolerance = 1.0e-13_dp
        type(ewald_parameters_t) :: ewald
        type(nbody_system_t) :: system
        logical :: ready = .false.
    end type
contains
    subroutine initialize_classical_pair(model, ok)
        type(classical_pair_model_t), intent(inout) :: model
        logical, intent(out) :: ok
        model%ready = .false.
        model%units = periodic_pair_units()
        call ewald_select(model%length, model%tolerance, model%ewald, ok)
        if (.not. ok) return
        call nbody_periodic(model%system, model%units, [electron(), electron()], &
                            model%ewald, ok)
        model%ready = ok
    end subroutine

    !> Wrap into the chart [-L/2, L/2) by the integer cell offset.
    pure subroutine classical_pair_wrap(position, length, wrapped)
        real(dp), intent(in) :: position(3), length
        real(dp), intent(out) :: wrapped(3)
        wrapped = position - length*real(floor(position/length + 0.5_dp), dp)
    end subroutine

    !> Two-particle state [x1, x2, v1, v2] at rest from centre and relative
    !> displacement: x1 = c + r/2, x2 = c - r/2.
    pure subroutine classical_pair_state(center, relative, state)
        real(dp), intent(in) :: center(3), relative(3)
        real(dp), intent(out) :: state(12)
        state = 0
        state(:3) = center + relative/2
        state(4:6) = center - relative/2
    end subroutine

    !> Uniform centre and independent uniform relative displacement on the
    !> torus (Park--Miller; exact integer arithmetic, one stream).
    subroutine classical_pair_sample(seed, length, center, relative)
        integer(int64), intent(inout) :: seed
        real(dp), intent(in) :: length
        real(dp), intent(out) :: center(3), relative(3)
        real(dp) :: draws(6)
        integer :: j
        do j = 1, 6
            seed = modulo(16807_int64*seed, 2147483647_int64)
            draws(j) = length*(real(seed, dp)/2147483647.0_dp)
        end do
        center = draws(:3); relative = draws(4:)
    end subroutine

    !> Centre of histogram bin b of 16 equal bins on x/L in [-1/2, 1/2].
    pure real(dp) function classical_pair_bin_center(b) result(x)
        integer, intent(in) :: b
        x = (2*b - 17)/32.0_dp
    end function
end module
