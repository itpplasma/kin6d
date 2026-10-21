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
    public :: equilibrium_field_t, labelled_field_t, equilibrium_state_t
    public :: flux_labelled_field_t

    ! Rows of the batched point evaluation values(EQ_VALUES, n): cylindrical
    ! B components (R, Z, phi), pressure, absolute poloidal flux and its
    ! (R, Z) gradient.
    integer, parameter, public :: EQ_BR = 1, EQ_BZ = 2, EQ_BPHI = 3, EQ_P = 4, &
        EQ_PSI = 5, EQ_PSI_R = 6, EQ_PSI_Z = 7, EQ_VALUES = 7
    integer, parameter, public :: EQ_GAUGE_AXIS = 1, EQ_GAUGE_BOUNDARY = 2
    ! Status of a label primitive the lane does not provide.
    integer, parameter, public :: EQ_UNAVAILABLE = -1

    ! The evaluator of B and p.
    type, abstract :: equilibrium_field_t
    contains
        procedure(field_sample), deferred :: sample
    end type

    ! A field with the label primitives of the D7 readback of axisymmetric
    ! nested states, the poloidal flux psi being the label.
    type, abstract, extends(equilibrium_field_t) :: labelled_field_t
    contains
        ! O-point of the label: point and absolute flux there.
        procedure(field_label_axis), deferred :: label_axis
        ! Magnetic axis: the fixed point of the poloidal field-line map,
        ! B_pol = 0, resolved beyond the label O-point where the lane can.
        procedure(field_magnetic_axis), deferred :: magnetic_axis
        ! coarea(j) = oint dl/(R |grad psi|) on the level psi = levels(j),
        ! levels(1) the axis flux, levels monotone towards the boundary;
        ! work(j) is the lane's work count for level j.
        procedure(field_label_coarea), deferred :: label_coarea
        ! Plasma volume and toroidal flux int B_phi dA over the plasma;
        ! EQ_UNAVAILABLE when the lane cannot integrate over its plasma.
        procedure(field_region_integrals), deferred :: region_integrals
        ! Default radial levels and angular quadrature budget for 0 inputs.
        procedure(field_resolution), deferred :: resolution
    end type
    ! A labelled field whose toroidal flux density must be integrated around
    ! the contour (map candidates need not have R B_phi constant on a label).
    type, abstract, extends(labelled_field_t) :: flux_labelled_field_t
    contains
        procedure(field_toroidal_density), deferred :: toroidal_density
    end type
    abstract interface
        ! Point values at points(2 or 3, n) = (R, Z[, phi]); omitted phi is
        ! zero. Components remain cylindrical; Cartesian conversion uses C.
        ! statuses(i) /= 0 marks a point
        ! outside the represented domain (its values are unspecified).
        subroutine field_sample(field, points, values, statuses)
            import :: equilibrium_field_t, dp
            class(equilibrium_field_t), intent(in) :: field
            real(dp), intent(in) :: points(:, :)
            real(dp), intent(out) :: values(:, :)
            integer, intent(out) :: statuses(:)
        end subroutine

        subroutine field_label_axis(field, axis, psi_axis, status)
            import :: labelled_field_t, dp
            class(labelled_field_t), intent(in) :: field
            real(dp), intent(out) :: axis(2), psi_axis
            integer, intent(out) :: status
        end subroutine

        subroutine field_magnetic_axis(field, axis, status)
            import :: labelled_field_t, dp
            class(labelled_field_t), intent(in) :: field
            real(dp), intent(out) :: axis(2)
            integer, intent(out) :: status
        end subroutine

        subroutine field_label_coarea(field, levels, angular_points, coarea, &
                                      status, work)
            import :: labelled_field_t, dp
            class(labelled_field_t), intent(in) :: field
            real(dp), intent(in) :: levels(:)
            integer, intent(in) :: angular_points
            real(dp), intent(out) :: coarea(:)
            integer, intent(out) :: status
            integer, intent(out), optional :: work(:)
        end subroutine

        subroutine field_region_integrals(field, volume, phi_tor, status)
            import :: labelled_field_t, dp
            class(labelled_field_t), intent(in) :: field
            real(dp), intent(out) :: volume, phi_tor
            integer, intent(out) :: status
        end subroutine

        subroutine field_toroidal_density(field, levels, angular_points, density, &
                                          status, work)
            import :: flux_labelled_field_t, dp
            class(flux_labelled_field_t), intent(in) :: field
            real(dp), intent(in) :: levels(:)
            integer, intent(in) :: angular_points
            real(dp), intent(out) :: density(:)
            integer, intent(out) :: status
            integer, intent(out), optional :: work(:)
        end subroutine

        subroutine field_resolution(field, radial_points, angular_points)
            import :: labelled_field_t
            class(labelled_field_t), intent(in) :: field
            integer, intent(inout) :: radial_points, angular_points
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
