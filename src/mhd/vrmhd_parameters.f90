module milca_vrmhd_parameters
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    private
    public :: vrmhd_parameters_t

    ! Temperatures are energies in joules; see docs/vrmhd-weak-form.md.
    type :: vrmhd_parameters_t
        real(dp) :: mi, mu0
        real(dp) :: t0 = 1.0_dp
        real(dp) :: eta0 = 0.0_dp, nu0 = 0.0_dp
        real(dp) :: chi_parallel0(2) = 0.0_dp, chi_perpendicular(2) = 0.0_dp
        real(dp) :: d_parallel = 0.0_dp, d_perpendicular = 0.0_dp
        real(dp) :: nu_parallel = 0.0_dp, equilibration_rate = 0.0_dp
        real(dp) :: epsilon_b = 1.0e-12_dp
        logical :: flux_limited(2) = .false.
        real(dp) :: flux_limit(2) = 0.3_dp
        real(dp) :: electron_mass = 9.1093837139e-31_dp
        real(dp) :: gradient_floor(2) = 1.0e-30_dp
    end type
end module milca_vrmhd_parameters
