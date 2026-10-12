module kin6d_equilibrium_state
    ! One physical equilibrium state S = (B, p, Lambda, C, pi) shared by the
    ! psi, map and free-B lanes (docs/equilibrium-core.md, D1). The core holds
    ! the field through an evaluator, never a vector potential:
    !   B, p      equilibrium_field_t (lane adapters extend it);
    !   Lambda    flux invariants and the poloidal-flux gauge;
    !   C         the chart X = H(q);
    !   pi        provenance (lane, seed, claim label).
    ! The poloidal flux is psi = (1/2 pi) int_disk B.dS = R A_phi with psi = 0
    ! on the symmetry axis. A lane whose domain does not reach the axis (fixed
    ! boundary) knows psi only up to a constant; it declares that by the gauge
    ! EQ_GAUGE_BOUNDARY, and psi_bound is then its value on the boundary.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_equilibrium_chart, only: ideal_chart_t
    implicit none
    private
    public :: equilibrium_field_t, equilibrium_state_t

    ! Rows of the batched point evaluation values(EQ_VALUES, n): cylindrical
    ! B components (R, Z, phi), pressure, absolute poloidal flux and its
    ! (R, Z) gradient.
    integer, parameter, public :: EQ_BR = 1, EQ_BZ = 2, EQ_BPHI = 3, EQ_P = 4, &
        EQ_PSI = 5, EQ_PSI_R = 6, EQ_PSI_Z = 7, EQ_VALUES = 7
    integer, parameter, public :: EQ_GAUGE_AXIS = 1, EQ_GAUGE_BOUNDARY = 2

    type, abstract :: equilibrium_field_t
    contains
        procedure(field_sample), deferred :: sample
    end type
    abstract interface
        ! Point values at points(2, n) = (R, Z); statuses(i) /= 0 marks a point
        ! outside the represented domain (its values are unspecified).
        subroutine field_sample(field, points, values, statuses)
            import :: equilibrium_field_t, dp
            class(equilibrium_field_t), intent(in) :: field
            real(dp), intent(in) :: points(:, :)
            real(dp), intent(out) :: values(:, :)
            integer, intent(out) :: statuses(:)
        end subroutine
    end interface

    type :: equilibrium_state_t
        class(equilibrium_field_t), allocatable :: field
        class(ideal_chart_t), allocatable :: chart
        real(dp) :: mu0 = 1.2566370614359173e-6_dp
        ! Lambda: flux on the plasma boundary in the declared gauge, vacuum
        ! F = R B_phi, toroidal current and toroidal flux (0 when unknown).
        integer :: psi_gauge = EQ_GAUGE_BOUNDARY
        real(dp) :: psi_bound = 0, f_vac = 0, current = 0, phi_tor = 0
        ! pi: producing lane, seed and claim label.
        character(32) :: lane = '', claim = 'candidate'
        character(256) :: seed = ''
    end type
end module kin6d_equilibrium_state
