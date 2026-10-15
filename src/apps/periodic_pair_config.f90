!> The one `&periodic_pair` namelist group shared by quantum_reference_demo
!> and classical_pair_demo: box, interaction and time grid that both
!> experiments must agree on. Defaults come from kin6d_species. A demo that
!> cannot realize a value must reject it explicitly.
module kin6d_periodic_pair_config
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_result_io, only: run_t
    use kin6d_species, only: units_t, species_t, periodic_pair_units, &
        periodic_pair_length, electron, pair_coupling
    implicit none
    private
    public :: periodic_pair_t, read_periodic_pair, write_periodic_pair, &
        declared_periodic_pair, periodic_pair_bins

    !> Relative-x slabs of width L/16 on [-L/2, L/2) shared by both demos.
    integer, parameter :: periodic_pair_bins = 16

    type :: periodic_pair_t
        !> Cube edge, particle mass and pair coupling q^2/(4 pi eps0) in the
        !> declared periodic-pair units; equally spaced frames on [0, time_max].
        real(dp) :: length, mass, coupling, time_max
        integer :: frames
    end type

contains

    !> The declared Q0/C1 electron pair (single source of defaults).
    function declared_periodic_pair() result(pair)
        type(periodic_pair_t) :: pair
        type(units_t) :: units
        type(species_t) :: species
        units = periodic_pair_units()
        species = electron()
        pair = periodic_pair_t(length=periodic_pair_length, mass=species%mass, &
            coupling=pair_coupling(units, species, species), &
            time_max=16.0_dp, frames=61)
    end function

    subroutine read_periodic_pair(run, pair)
        type(run_t), intent(in) :: run
        type(periodic_pair_t), intent(out) :: pair
        real(dp) :: length, mass, coupling, time_max
        integer :: frames, unit, ios
        character(len=256) :: msg
        namelist /periodic_pair/ length, mass, coupling, time_max, frames
        pair = declared_periodic_pair()
        length = pair%length; mass = pair%mass; coupling = pair%coupling
        time_max = pair%time_max
        frames = pair%frames
        ios = 0; msg = ''
        call run%open_config(unit)
        if (run%has_config()) read (unit, nml=periodic_pair, iostat=ios, iomsg=msg)
        call run%close_config(unit, 'periodic_pair', ios, msg)
        pair = periodic_pair_t(length, mass, coupling, time_max, frames)
    end subroutine

    subroutine write_periodic_pair(unit, pair)
        integer, intent(in) :: unit
        type(periodic_pair_t), intent(in) :: pair
        real(dp) :: length, mass, coupling, time_max
        integer :: frames
        namelist /periodic_pair/ length, mass, coupling, time_max, frames
        length = pair%length; mass = pair%mass; coupling = pair%coupling
        time_max = pair%time_max
        frames = pair%frames
        write (unit, nml=periodic_pair)
    end subroutine
end module kin6d_periodic_pair_config
