module vrmhd_leaf_dispatch
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use milca_vrmhd_leaves
    implicit none
    integer, parameter :: nleaf = 27, maxargs = 33
    integer, parameter :: nargs(nleaf) = [ &
        3, 3, 3, 1, 4, 2, 11, 6, 6, 13, 7, 12, 14, &
        25, 2, 9, 9, 10, 8, 24, 16, 9, 5, 6, 2, 7, 33]
    character(32), parameter :: names(nleaf) = [character(32) :: &
        'eta', &
        'nu', &
        'chi_parallel', &
        'chi_perpendicular', &
        'pressure', &
        'mass_density', &
        'mass_momentum', &
        'mass_temperature', &
        'mass_induction', &
        'ampere', &
        'density_advection', &
        'transport', &
        'conduction', &
        'momentum_advection', &
        'pressure_work', &
        'lorentz', &
        'induction_ideal', &
        'induction_resistive', &
        'ohmic', &
        'viscosity', &
        'viscous_heat', &
        'temperature_advection', &
        'compression', &
        'equipartition', &
        'source', &
        'momentum_mass_source', &
        'momentum_density_diffusion']
contains
    subroutine evaluate(leaf, x, y)
        integer, intent(in) :: leaf
        real(dp), intent(in) :: x(maxargs)
        real(dp), intent(out) :: y
        select case (leaf)
        case (1)
            call vrmhd_eta( &
                x(1), x(2), x(3), y)
        case (2)
            call vrmhd_nu( &
                x(1), x(2), x(3), y)
        case (3)
            call vrmhd_chi_parallel( &
                x(1), x(2), x(3), y)
        case (4)
            call vrmhd_chi_perpendicular( &
                x(1), y)
        case (5)
            call vrmhd_pressure( &
                x(1), x(2), x(3), x(4), y)
        case (6)
            call vrmhd_mass_density( &
                x(1), x(2), y)
        case (7)
            call vrmhd_mass_momentum( &
                x(1), x(2), x(3), x(4), x(5), &
                x(6), x(7), x(8), x(9), x(10), &
                x(11), y)
        case (8)
            call vrmhd_mass_temperature( &
                x(1), x(2), x(3), x(4), x(5), &
                x(6), y)
        case (9)
            call vrmhd_mass_induction( &
                x(1), x(2), x(3), x(4), x(5), &
                x(6), y)
        case (10)
            call vrmhd_ampere( &
                x(1), x(2), x(3), x(4), x(5), &
                x(6), x(7), x(8), x(9), x(10), &
                x(11), x(12), x(13), y)
        case (11)
            call vrmhd_density_advection( &
                x(1), x(2), x(3), x(4), x(5), &
                x(6), x(7), y)
        case (12)
            call vrmhd_transport( &
                x(1), x(2), x(3), x(4), x(5), &
                x(6), x(7), x(8), x(9), x(10), &
                x(11), x(12), y)
        case (13)
            call vrmhd_conduction( &
                x(1), x(2), x(3), x(4), x(5), &
                x(6), x(7), x(8), x(9), x(10), &
                x(11), x(12), x(13), x(14), y)
        case (14)
            call vrmhd_momentum_advection( &
                x(1), x(2), x(3), x(4), x(5), &
                x(6), x(7), x(8), x(9), x(10), &
                x(11), x(12), x(13), x(14), x(15), &
                x(16), x(17), x(18), x(19), x(20), &
                x(21), x(22), x(23), x(24), x(25), &
                y)
        case (15)
            call vrmhd_pressure_work( &
                x(1), x(2), y)
        case (16)
            call vrmhd_lorentz( &
                x(1), x(2), x(3), x(4), x(5), &
                x(6), x(7), x(8), x(9), y)
        case (17)
            call vrmhd_induction_ideal( &
                x(1), x(2), x(3), x(4), x(5), &
                x(6), x(7), x(8), x(9), y)
        case (18)
            call vrmhd_induction_resistive( &
                x(1), x(2), x(3), x(4), x(5), &
                x(6), x(7), x(8), x(9), x(10), &
                y)
        case (19)
            call vrmhd_ohmic( &
                x(1), x(2), x(3), x(4), x(5), &
                x(6), x(7), x(8), y)
        case (20)
            call vrmhd_viscosity( &
                x(1), x(2), x(3), x(4), x(5), &
                x(6), x(7), x(8), x(9), x(10), &
                x(11), x(12), x(13), x(14), x(15), &
                x(16), x(17), x(18), x(19), x(20), &
                x(21), x(22), x(23), x(24), y)
        case (21)
            call vrmhd_viscous_heat( &
                x(1), x(2), x(3), x(4), x(5), &
                x(6), x(7), x(8), x(9), x(10), &
                x(11), x(12), x(13), x(14), x(15), &
                x(16), y)
        case (22)
            call vrmhd_temperature_advection( &
                x(1), x(2), x(3), x(4), x(5), &
                x(6), x(7), x(8), x(9), y)
        case (23)
            call vrmhd_compression( &
                x(1), x(2), x(3), x(4), x(5), &
                y)
        case (24)
            call vrmhd_equipartition( &
                x(1), x(2), x(3), x(4), x(5), &
                x(6), y)
        case (25)
            call vrmhd_source( &
                x(1), x(2), y)
        case (26)
            call vrmhd_momentum_mass_source( &
                x(1), x(2), x(3), x(4), x(5), &
                x(6), x(7), y)
        case (27)
            call vrmhd_momentum_density_diffusion( &
                x(1), x(2), x(3), x(4), x(5), &
                x(6), x(7), x(8), x(9), x(10), &
                x(11), x(12), x(13), x(14), x(15), &
                x(16), x(17), x(18), x(19), x(20), &
                x(21), x(22), x(23), x(24), x(25), &
                x(26), x(27), x(28), x(29), x(30), &
                x(31), x(32), x(33), y)
        case default
            error stop 'unknown VRMHD leaf'
        end select
    end subroutine

    subroutine derivative(leaf, x, d, weight, tangent, gradient)
        integer, intent(in) :: leaf
        real(dp), intent(in) :: x(maxargs), d(maxargs), weight
        real(dp), intent(out) :: tangent, gradient(maxargs)
        gradient = 0
        select case (leaf)
        case (1)
            call vrmhd_eta_jvp( &
                x(1), x(2), x(3), d(1), d(2), &
                d(3), tangent)
            call vrmhd_eta_vjp( &
                x(1), x(2), x(3), weight, gradient(1), &
                gradient(2), gradient(3))
        case (2)
            call vrmhd_nu_jvp( &
                x(1), x(2), x(3), d(1), d(2), &
                d(3), tangent)
            call vrmhd_nu_vjp( &
                x(1), x(2), x(3), weight, gradient(1), &
                gradient(2), gradient(3))
        case (3)
            call vrmhd_chi_parallel_jvp( &
                x(1), x(2), x(3), d(1), d(2), &
                d(3), tangent)
            call vrmhd_chi_parallel_vjp( &
                x(1), x(2), x(3), weight, gradient(1), &
                gradient(2), gradient(3))
        case (4)
            call vrmhd_chi_perpendicular_jvp( &
                x(1), d(1), tangent)
            call vrmhd_chi_perpendicular_vjp( &
                x(1), weight, gradient(1))
        case (5)
            call vrmhd_pressure_jvp( &
                x(1), x(2), x(3), x(4), d(1), &
                d(2), d(3), d(4), tangent)
            call vrmhd_pressure_vjp( &
                x(1), x(2), x(3), x(4), weight, &
                gradient(1), gradient(2), gradient(3), gradient(4))
        case (6)
            call vrmhd_mass_density_jvp( &
                x(1), x(2), d(1), d(2), tangent)
            call vrmhd_mass_density_vjp( &
                x(1), x(2), weight, gradient(1), gradient(2))
        case (7)
            call vrmhd_mass_momentum_jvp( &
                x(1), x(2), x(3), x(4), x(5), &
                x(6), x(7), x(8), x(9), x(10), &
                x(11), d(1), d(2), d(3), d(4), &
                d(5), d(6), d(7), d(8), d(9), &
                d(10), d(11), tangent)
            call vrmhd_mass_momentum_vjp( &
                x(1), x(2), x(3), x(4), x(5), &
                x(6), x(7), x(8), x(9), x(10), &
                x(11), weight, gradient(1), gradient(2), gradient(3), &
                gradient(4), gradient(5), gradient(6), gradient(7), gradient(8), &
                gradient(9), gradient(10), gradient(11))
        case (8)
            call vrmhd_mass_temperature_jvp( &
                x(1), x(2), x(3), x(4), x(5), &
                x(6), d(1), d(2), d(3), d(4), &
                d(5), d(6), tangent)
            call vrmhd_mass_temperature_vjp( &
                x(1), x(2), x(3), x(4), x(5), &
                x(6), weight, gradient(1), gradient(2), gradient(3), &
                gradient(4), gradient(5), gradient(6))
        case (9)
            call vrmhd_mass_induction_jvp( &
                x(1), x(2), x(3), x(4), x(5), &
                x(6), d(1), d(2), d(3), d(4), &
                d(5), d(6), tangent)
            call vrmhd_mass_induction_vjp( &
                x(1), x(2), x(3), x(4), x(5), &
                x(6), weight, gradient(1), gradient(2), gradient(3), &
                gradient(4), gradient(5), gradient(6))
        case (10)
            call vrmhd_ampere_jvp( &
                x(1), x(2), x(3), x(4), x(5), &
                x(6), x(7), x(8), x(9), x(10), &
                x(11), x(12), x(13), d(1), d(2), &
                d(3), d(4), d(5), d(6), d(7), &
                d(8), d(9), d(10), d(11), d(12), &
                d(13), tangent)
            call vrmhd_ampere_vjp( &
                x(1), x(2), x(3), x(4), x(5), &
                x(6), x(7), x(8), x(9), x(10), &
                x(11), x(12), x(13), weight, gradient(1), &
                gradient(2), gradient(3), gradient(4), gradient(5), gradient(6), &
                gradient(7), gradient(8), gradient(9), gradient(10), gradient(11), &
                gradient(12), gradient(13))
        case (11)
            call vrmhd_density_advection_jvp( &
                x(1), x(2), x(3), x(4), x(5), &
                x(6), x(7), d(1), d(2), d(3), &
                d(4), d(5), d(6), d(7), tangent)
            call vrmhd_density_advection_vjp( &
                x(1), x(2), x(3), x(4), x(5), &
                x(6), x(7), weight, gradient(1), gradient(2), &
                gradient(3), gradient(4), gradient(5), gradient(6), gradient(7))
        case (12)
            call vrmhd_transport_jvp( &
                x(1), x(2), x(3), x(4), x(5), &
                x(6), x(7), x(8), x(9), x(10), &
                x(11), x(12), d(1), d(2), d(3), &
                d(4), d(5), d(6), d(7), d(8), &
                d(9), d(10), d(11), d(12), tangent)
            call vrmhd_transport_vjp( &
                x(1), x(2), x(3), x(4), x(5), &
                x(6), x(7), x(8), x(9), x(10), &
                x(11), x(12), weight, gradient(1), gradient(2), &
                gradient(3), gradient(4), gradient(5), gradient(6), gradient(7), &
                gradient(8), gradient(9), gradient(10), gradient(11), gradient(12))
        case (13)
            call vrmhd_conduction_jvp( &
                x(1), x(2), x(3), x(4), x(5), &
                x(6), x(7), x(8), x(9), x(10), &
                x(11), x(12), x(13), x(14), d(1), &
                d(2), d(3), d(4), d(5), d(6), &
                d(7), d(8), d(9), d(10), d(11), &
                d(12), d(13), d(14), tangent)
            call vrmhd_conduction_vjp( &
                x(1), x(2), x(3), x(4), x(5), &
                x(6), x(7), x(8), x(9), x(10), &
                x(11), x(12), x(13), x(14), weight, &
                gradient(1), gradient(2), gradient(3), gradient(4), gradient(5), &
                gradient(6), gradient(7), gradient(8), gradient(9), gradient(10), &
                gradient(11), gradient(12), gradient(13), gradient(14))
        case (14)
            call vrmhd_momentum_advection_jvp( &
                x(1), x(2), x(3), x(4), x(5), &
                x(6), x(7), x(8), x(9), x(10), &
                x(11), x(12), x(13), x(14), x(15), &
                x(16), x(17), x(18), x(19), x(20), &
                x(21), x(22), x(23), x(24), x(25), &
                d(1), d(2), d(3), d(4), d(5), &
                d(6), d(7), d(8), d(9), d(10), &
                d(11), d(12), d(13), d(14), d(15), &
                d(16), d(17), d(18), d(19), d(20), &
                d(21), d(22), d(23), d(24), d(25), &
                tangent)
            call vrmhd_momentum_advection_vjp( &
                x(1), x(2), x(3), x(4), x(5), &
                x(6), x(7), x(8), x(9), x(10), &
                x(11), x(12), x(13), x(14), x(15), &
                x(16), x(17), x(18), x(19), x(20), &
                x(21), x(22), x(23), x(24), x(25), &
                weight, gradient(1), gradient(2), gradient(3), gradient(4), &
                gradient(5), gradient(6), gradient(7), gradient(8), gradient(9), &
                gradient(10), gradient(11), gradient(12), gradient(13), gradient(14), &
                gradient(15), gradient(16), gradient(17), gradient(18), gradient(19), &
                gradient(20), gradient(21), gradient(22), gradient(23), gradient(24), &
                gradient(25))
        case (15)
            call vrmhd_pressure_work_jvp( &
                x(1), x(2), d(1), d(2), tangent)
            call vrmhd_pressure_work_vjp( &
                x(1), x(2), weight, gradient(1), gradient(2))
        case (16)
            call vrmhd_lorentz_jvp( &
                x(1), x(2), x(3), x(4), x(5), &
                x(6), x(7), x(8), x(9), d(1), &
                d(2), d(3), d(4), d(5), d(6), &
                d(7), d(8), d(9), tangent)
            call vrmhd_lorentz_vjp( &
                x(1), x(2), x(3), x(4), x(5), &
                x(6), x(7), x(8), x(9), weight, &
                gradient(1), gradient(2), gradient(3), gradient(4), gradient(5), &
                gradient(6), gradient(7), gradient(8), gradient(9))
        case (17)
            call vrmhd_induction_ideal_jvp( &
                x(1), x(2), x(3), x(4), x(5), &
                x(6), x(7), x(8), x(9), d(1), &
                d(2), d(3), d(4), d(5), d(6), &
                d(7), d(8), d(9), tangent)
            call vrmhd_induction_ideal_vjp( &
                x(1), x(2), x(3), x(4), x(5), &
                x(6), x(7), x(8), x(9), weight, &
                gradient(1), gradient(2), gradient(3), gradient(4), gradient(5), &
                gradient(6), gradient(7), gradient(8), gradient(9))
        case (18)
            call vrmhd_induction_resistive_jvp( &
                x(1), x(2), x(3), x(4), x(5), &
                x(6), x(7), x(8), x(9), x(10), &
                d(1), d(2), d(3), d(4), d(5), &
                d(6), d(7), d(8), d(9), d(10), &
                tangent)
            call vrmhd_induction_resistive_vjp( &
                x(1), x(2), x(3), x(4), x(5), &
                x(6), x(7), x(8), x(9), x(10), &
                weight, gradient(1), gradient(2), gradient(3), gradient(4), &
                gradient(5), gradient(6), gradient(7), gradient(8), gradient(9), &
                gradient(10))
        case (19)
            call vrmhd_ohmic_jvp( &
                x(1), x(2), x(3), x(4), x(5), &
                x(6), x(7), x(8), d(1), d(2), &
                d(3), d(4), d(5), d(6), d(7), &
                d(8), tangent)
            call vrmhd_ohmic_vjp( &
                x(1), x(2), x(3), x(4), x(5), &
                x(6), x(7), x(8), weight, gradient(1), &
                gradient(2), gradient(3), gradient(4), gradient(5), gradient(6), &
                gradient(7), gradient(8))
        case (20)
            call vrmhd_viscosity_jvp( &
                x(1), x(2), x(3), x(4), x(5), &
                x(6), x(7), x(8), x(9), x(10), &
                x(11), x(12), x(13), x(14), x(15), &
                x(16), x(17), x(18), x(19), x(20), &
                x(21), x(22), x(23), x(24), d(1), &
                d(2), d(3), d(4), d(5), d(6), &
                d(7), d(8), d(9), d(10), d(11), &
                d(12), d(13), d(14), d(15), d(16), &
                d(17), d(18), d(19), d(20), d(21), &
                d(22), d(23), d(24), tangent)
            call vrmhd_viscosity_vjp( &
                x(1), x(2), x(3), x(4), x(5), &
                x(6), x(7), x(8), x(9), x(10), &
                x(11), x(12), x(13), x(14), x(15), &
                x(16), x(17), x(18), x(19), x(20), &
                x(21), x(22), x(23), x(24), weight, &
                gradient(1), gradient(2), gradient(3), gradient(4), gradient(5), &
                gradient(6), gradient(7), gradient(8), gradient(9), gradient(10), &
                gradient(11), gradient(12), gradient(13), gradient(14), gradient(15), &
                gradient(16), gradient(17), gradient(18), gradient(19), gradient(20), &
                gradient(21), gradient(22), gradient(23), gradient(24))
        case (21)
            call vrmhd_viscous_heat_jvp( &
                x(1), x(2), x(3), x(4), x(5), &
                x(6), x(7), x(8), x(9), x(10), &
                x(11), x(12), x(13), x(14), x(15), &
                x(16), d(1), d(2), d(3), d(4), &
                d(5), d(6), d(7), d(8), d(9), &
                d(10), d(11), d(12), d(13), d(14), &
                d(15), d(16), tangent)
            call vrmhd_viscous_heat_vjp( &
                x(1), x(2), x(3), x(4), x(5), &
                x(6), x(7), x(8), x(9), x(10), &
                x(11), x(12), x(13), x(14), x(15), &
                x(16), weight, gradient(1), gradient(2), gradient(3), &
                gradient(4), gradient(5), gradient(6), gradient(7), gradient(8), &
                gradient(9), gradient(10), gradient(11), gradient(12), gradient(13), &
                gradient(14), gradient(15), gradient(16))
        case (22)
            call vrmhd_temperature_advection_jvp( &
                x(1), x(2), x(3), x(4), x(5), &
                x(6), x(7), x(8), x(9), d(1), &
                d(2), d(3), d(4), d(5), d(6), &
                d(7), d(8), d(9), tangent)
            call vrmhd_temperature_advection_vjp( &
                x(1), x(2), x(3), x(4), x(5), &
                x(6), x(7), x(8), x(9), weight, &
                gradient(1), gradient(2), gradient(3), gradient(4), gradient(5), &
                gradient(6), gradient(7), gradient(8), gradient(9))
        case (23)
            call vrmhd_compression_jvp( &
                x(1), x(2), x(3), x(4), x(5), &
                d(1), d(2), d(3), d(4), d(5), &
                tangent)
            call vrmhd_compression_vjp( &
                x(1), x(2), x(3), x(4), x(5), &
                weight, gradient(1), gradient(2), gradient(3), gradient(4), &
                gradient(5))
        case (24)
            call vrmhd_equipartition_jvp( &
                x(1), x(2), x(3), x(4), x(5), &
                x(6), d(1), d(2), d(3), d(4), &
                d(5), d(6), tangent)
            call vrmhd_equipartition_vjp( &
                x(1), x(2), x(3), x(4), x(5), &
                x(6), weight, gradient(1), gradient(2), gradient(3), &
                gradient(4), gradient(5), gradient(6))
        case (25)
            call vrmhd_source_jvp( &
                x(1), x(2), d(1), d(2), tangent)
            call vrmhd_source_vjp( &
                x(1), x(2), weight, gradient(1), gradient(2))
        case (26)
            call vrmhd_momentum_mass_source_jvp( &
                x(1), x(2), x(3), x(4), x(5), &
                x(6), x(7), d(1), d(2), d(3), &
                d(4), d(5), d(6), d(7), tangent)
            call vrmhd_momentum_mass_source_vjp( &
                x(1), x(2), x(3), x(4), x(5), &
                x(6), x(7), weight, gradient(1), gradient(2), &
                gradient(3), gradient(4), gradient(5), gradient(6), gradient(7))
        case (27)
            call vrmhd_momentum_density_diffusion_jvp( &
                x(1), x(2), x(3), x(4), x(5), &
                x(6), x(7), x(8), x(9), x(10), &
                x(11), x(12), x(13), x(14), x(15), &
                x(16), x(17), x(18), x(19), x(20), &
                x(21), x(22), x(23), x(24), x(25), &
                x(26), x(27), x(28), x(29), x(30), &
                x(31), x(32), x(33), d(1), d(2), &
                d(3), d(4), d(5), d(6), d(7), &
                d(8), d(9), d(10), d(11), d(12), &
                d(13), d(14), d(15), d(16), d(17), &
                d(18), d(19), d(20), d(21), d(22), &
                d(23), d(24), d(25), d(26), d(27), &
                d(28), d(29), d(30), d(31), d(32), &
                d(33), tangent)
            call vrmhd_momentum_density_diffusion_vjp( &
                x(1), x(2), x(3), x(4), x(5), &
                x(6), x(7), x(8), x(9), x(10), &
                x(11), x(12), x(13), x(14), x(15), &
                x(16), x(17), x(18), x(19), x(20), &
                x(21), x(22), x(23), x(24), x(25), &
                x(26), x(27), x(28), x(29), x(30), &
                x(31), x(32), x(33), weight, gradient(1), &
                gradient(2), gradient(3), gradient(4), gradient(5), gradient(6), &
                gradient(7), gradient(8), gradient(9), gradient(10), gradient(11), &
                gradient(12), gradient(13), gradient(14), gradient(15), gradient(16), &
                gradient(17), gradient(18), gradient(19), gradient(20), gradient(21), &
                gradient(22), gradient(23), gradient(24), gradient(25), gradient(26), &
                gradient(27), gradient(28), gradient(29), gradient(30), gradient(31), &
                gradient(32), gradient(33))
        case default
            error stop 'unknown VRMHD derivative'
        end select
    end subroutine
end module vrmhd_leaf_dispatch
