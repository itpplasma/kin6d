module kin6d_island_closure
    ! I1 declaration only. B fitting/projection belongs to I5, Taylor regions to I2/I8.
    implicit none
    private
    public :: island_closure_t, island_closure_check
    public :: ISLAND_INVALID, ISLAND_UNSUPPORTED, ISLAND_DELEGATE_B, ISLAND_DELEGATE_C
    public :: ISLAND_NULL_FIELD, ISLAND_BAD_DATA
    integer, parameter :: ISLAND_NULL_FIELD = 14, ISLAND_BAD_DATA = 15
    integer, parameter :: ISLAND_INVALID = 10, ISLAND_UNSUPPORTED = 11
    integer, parameter :: ISLAND_DELEGATE_B = 12, ISLAND_DELEGATE_C = 13
    type :: island_closure_t
        character(24) :: pressure = 'A', current = 'D', drive = 'I_phi'
        character(24) :: heat_boundary = 'dirichlet', potential_boundary = 'dirichlet'
        character(24) :: target_class = 'regularised'
    end type
contains
    subroutine island_closure_check(closure, status, message)
        type(island_closure_t), intent(in) :: closure
        integer, intent(out) :: status
        character(*), intent(out) :: message
        status = ISLAND_INVALID
        message = 'Select pressure A/B_frozen_fit/C and compatible current D/C'
        select case (trim(closure%pressure))
        case ('A', 'B_frozen_fit', 'B_projection')
            if (closure%current /= 'D') return
        case ('C')
            if (closure%current /= 'C') return
        case default
            return
        end select
        message = 'Prescribe exactly one drive: I_phi or V, never both'
        if (closure%drive /= 'I_phi' .and. closure%drive /= 'V') return
        status = ISLAND_UNSUPPORTED
        message = 'I1 supports Dirichlet heat data in the seed only'
        if (closure%heat_boundary /= 'dirichlet') return
        message = 'I1 supports Dirichlet potential gauge and floating conductors only'
        if (closure%potential_boundary /= 'dirichlet') return
        message = 'Finite-ratio A+D yields a regularised candidate, not an ideal claim'
        if (closure%target_class /= 'regularised') return
        select case (trim(closure%pressure))
        case ('B_frozen_fit', 'B_projection')
            status = ISLAND_DELEGATE_B
            message = 'Closure B delegated to I5: fit/projection required'
        case ('C')
            status = ISLAND_DELEGATE_C
            message = 'Closure C delegated to I2/I8: regional invariants required'
        case default
            status = 0
            message = 'A+D: transport pressure and resistive Ohm/current continuity'
        end select
    end subroutine
end module
