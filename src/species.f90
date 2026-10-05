!> Species and declared reduced units shared by every microscopic case.
!>
!> Reduced units declare the values of hbar and of the Coulomb constant
!> 1/(4 pi eps0); masses and charges are in declared SI mass/charge units.
!> The SI length and time of one reduced unit then follow (species_unit_map).
!> Pair coupling g = coulomb q_i q_j; every physical scale is a generated leaf
!> from gen/kin6d_symbolic_physics.f90, never retyped here.
module kin6d_species
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_species_leaves, only: species_pair_scales, species_sommerfeld, &
                                    species_plasma_scales, species_unit_map
    implicit none
    private
    public :: units_t, species_t, pair_scales_t, plasma_scales_t
    public :: reduced_units, atomic_units, periodic_pair_units, electron
    public :: pair_coupling, pair_scales, sommerfeld_parameter, plasma_scales
    public :: periodic_pair_length
    public :: hbar_si, coulomb_si, electron_mass_si, elementary_charge_si

    ! CODATA 2018 (exact h and e; m_e and 1/(4 pi eps0) measured).
    real(dp), parameter :: hbar_si = 6.62607015e-34_dp/(2*3.14159265358979323846_dp)
    real(dp), parameter :: elementary_charge_si = 1.602176634e-19_dp
    real(dp), parameter :: electron_mass_si = 9.1093837015e-31_dp
    real(dp), parameter :: coulomb_si = 8.9875517923e9_dp
    !> Cell edge of the matched classical/quantum periodic pair (C1, Q0).
    real(dp), parameter :: periodic_pair_length = 6.2831853071795864769_dp

    type :: units_t
        real(dp) :: hbar = 1.0_dp
        !> Coulomb constant 1/(4 pi eps0) in reduced units.
        real(dp) :: coulomb = 1.0_dp
        !> SI value of one reduced mass, charge, length and time unit.
        real(dp) :: mass_si = 0.0_dp, charge_si = 0.0_dp
        real(dp) :: length_si = 0.0_dp, time_si = 0.0_dp
    end type

    type :: species_t
        character(16) :: name = ''
        real(dp) :: charge = 0.0_dp, mass = 0.0_dp
    end type

    type :: pair_scales_t
        real(dp) :: coupling = 0.0_dp, reduced_mass = 0.0_dp
        real(dp) :: bohr_radius = 0.0_dp, rydberg = 0.0_dp
    end type

    !> One responsive species at number density n and temperature T (k_B=1).
    type :: plasma_scales_t
        real(dp) :: density = 0.0_dp, temperature = 0.0_dp
        real(dp) :: wigner_seitz = 0.0_dp, gamma = 0.0_dp, fermi = 0.0_dp
        real(dp) :: theta = 0.0_dp, rs = 0.0_dp, debye_length = 0.0_dp
        real(dp) :: plasma_frequency = 0.0_dp
    end type
contains
    !> Units with declared hbar and Coulomb constant; mass_si/charge_si are
    !> the SI values of the reduced mass and charge units.
    pure function reduced_units(hbar, coulomb, mass_si, charge_si) result(units)
        real(dp), intent(in) :: hbar, coulomb, mass_si, charge_si
        type(units_t) :: units
        units%hbar = hbar; units%coulomb = coulomb
        units%mass_si = mass_si; units%charge_si = charge_si
        call species_unit_map(hbar, coulomb, hbar_si, coulomb_si, charge_si, &
                              mass_si, units%length_si, units%time_si)
    end function

    !> Hartree atomic units: hbar = m_e = e = 1/(4 pi eps0) = 1.
    pure function atomic_units() result(units)
        type(units_t) :: units
        units = reduced_units(1.0_dp, 1.0_dp, electron_mass_si, &
                              elementary_charge_si)
    end function

    !> Units of the matched classical/quantum periodic pair (C1, Q0):
    !> hbar = m_e = e = 1 and e^2/(4 pi eps0) = 2, i.e. length unit 2 a_0.
    pure function periodic_pair_units() result(units)
        type(units_t) :: units
        units = reduced_units(1.0_dp, 2.0_dp, electron_mass_si, &
                              elementary_charge_si)
    end function

    !> Electron in units whose mass/charge units are m_e and e.
    pure function electron() result(s)
        type(species_t) :: s
        s = species_t('electron', -1.0_dp, 1.0_dp)
    end function

    pure real(dp) function pair_coupling(units, a, b) result(g)
        type(units_t), intent(in) :: units
        type(species_t), intent(in) :: a, b
        g = units%coulomb*a%charge*b%charge
    end function

    pure function pair_scales(units, a, b) result(s)
        type(units_t), intent(in) :: units
        type(species_t), intent(in) :: a, b
        type(pair_scales_t) :: s
        s%coupling = pair_coupling(units, a, b)
        call species_pair_scales(units%hbar, a%mass, b%mass, s%coupling, &
                                 s%reduced_mass, s%bohr_radius, s%rydberg)
    end function

    !> Sommerfeld parameter g/(hbar v) at relative speed v.
    pure real(dp) function sommerfeld_parameter(units, a, b, v) result(eta)
        type(units_t), intent(in) :: units
        type(species_t), intent(in) :: a, b
        real(dp), intent(in) :: v
        call species_sommerfeld(pair_coupling(units, a, b), units%hbar, v, eta)
    end function

    pure function plasma_scales(units, s, density, temperature) result(p)
        type(units_t), intent(in) :: units
        type(species_t), intent(in) :: s
        real(dp), intent(in) :: density, temperature
        type(plasma_scales_t) :: p
        p%density = density; p%temperature = temperature
        call species_plasma_scales(units%hbar, s%mass, pair_coupling(units, s, s), &
                                   density, temperature, p%wigner_seitz, p%gamma, &
                                   p%fermi, p%theta, p%rs, p%debye_length, &
                                   p%plasma_frequency)
    end function
end module kin6d_species
