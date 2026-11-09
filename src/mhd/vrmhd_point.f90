module milca_vrmhd_point
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use milca_vrmhd_parameters, only: vrmhd_parameters_t
    use milca_vrmhd_dispatch, only: vrmhd_leaf
    use milca_vrmhd_leaf_ids, only: &
        leaf_ampere, leaf_chi_flux_limited, leaf_chi_parallel, &
        leaf_compression, leaf_conduction, leaf_density_advection, &
        leaf_equipartition, leaf_eta, leaf_induction_ideal, &
        leaf_induction_resistive, leaf_lorentz, leaf_mass_density, &
        leaf_mass_induction, leaf_mass_momentum, leaf_mass_temperature, &
        leaf_momentum_advection, leaf_momentum_density_diffusion, &
        leaf_momentum_mass_source, leaf_nu, leaf_ohmic, &
        leaf_parallel_gradient, leaf_pressure, leaf_pressure_work, &
        leaf_source, leaf_temperature_advection, leaf_trace, &
        leaf_transport, leaf_viscosity, leaf_viscous_heat
    implicit none
    private
    public :: vrmhd_point, vrmhd_point_matrix
    integer, parameter :: ntest(6) = [4, 12, 6, 6, 4, 4]
    !> vrmhd_point evaluations of one vrmhd_point_matrix call.
    integer, parameter, public :: VRMHD_POINT_MATRIX_CALLS = 45*sum(ntest)
    integer, parameter :: test_slots(12, 6) = reshape([ &
        1, 2, 3, 4, 0, 0, 0, 0, 0, 0, 0, 0, &
        5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, &
        5, 6, 7, 17, 18, 19, 0, 0, 0, 0, 0, 0, &
        5, 6, 7, 17, 18, 19, 0, 0, 0, 0, 0, 0, &
        1, 2, 3, 4, 0, 0, 0, 0, 0, 0, 0, 0, &
        1, 2, 3, 4, 0, 0, 0, 0, 0, 0, 0, 0], [12, 6])

contains
    ! Pointwise bilinear form of the JVP: kmat(:, j, row) is the tangent
    ! against unit tests for unit trial data j (33 values, 12 rates).
    subroutine vrmhd_point_matrix(p, v, vd, source, shift, kmat)
        type(vrmhd_parameters_t), intent(in) :: p
        real(dp), intent(in) :: v(33), vd(12), source(6), shift
        real(dp), intent(out) :: kmat(19, 45, 6)
        real(dp) :: du(33), dud(12), test(19), value, tangent
        integer :: row, it, j
        kmat = 0
        do row = 1, 6
            do it = 1, ntest(row)
                test = 0; test(test_slots(it, row)) = 1
                do j = 1, 33
                    du = 0; du(j) = 1; dud = 0
                    call vrmhd_point(p, row, v, du, vd, dud, test, source, value, &
                        tangent)
                    kmat(test_slots(it, row), j, row) = tangent
                end do
                do j = 1, 12
                    du = 0; dud = 0; dud(j) = shift
                    call vrmhd_point(p, row, v, du, vd, dud, test, source, value, &
                        tangent)
                    kmat(test_slots(it, row), 33 + j, row) = tangent
                end do
            end do
        end do
    end subroutine

    subroutine vrmhd_point(p, row, u, du, ud, dud, test, source, value, tangent)
        type(vrmhd_parameters_t), intent(in) :: p
        integer, intent(in) :: row
        real(dp), intent(in) :: u(33), du(33), ud(12), dud(12)
        real(dp), intent(in) :: test(19), source(6)
        real(dp), intent(out) :: value, tangent
        real(dp) :: eta, deta, nu, dnu, chi, dchi, pressure, dpressure
        real(dp) :: divv, ddivv, divw, h, dh, mass, bg, dbg
        real(dp) :: zero(24)
        integer :: species, it
        zero = 0; value = 0; tangent = 0
        ! Transport coefficients and traces only for the rows that use them.
        if (row == 3 .or. row == 5) call vrmhd_leaf(leaf_eta, [p%eta0, u(11), p%t0, &
            p%eta_exponent], [0.0_dp, du(11), 0.0_dp, 0.0_dp], eta, deta)
        if (row == 2 .or. row == 6) call vrmhd_leaf(leaf_nu, [p%nu0, u(12), p%t0, &
            p%nu_exponent], [0.0_dp, du(12), 0.0_dp, 0.0_dp], nu, dnu)
        if (row >= 5) call vrmhd_leaf(leaf_trace, [u(16), u(20), u(24)], [du(16), &
            du(20), du(24)], divv, ddivv)
        if (row == 2) call vrmhd_leaf(leaf_trace, [test(8), test(12), test(16)], &
            zero(:3), divw, h)
        select case (row)
        case (1)
            call add(leaf_mass_density, [ud(1), test(1)], [dud(1), 0.0_dp], 1.0_dp)
            call add(leaf_density_advection, [u(1:4), test(2:4)], [du(1:4), &
                zero(:3)], 1.0_dp)
            call add(leaf_transport, [p%d_parallel, p%d_perpendicular, &
                p%epsilon_b, u(25:27), &
                u(13:15), test(2:4)], [zero(:3), du(25:27), &
                du(13:15), zero(:3)], 1.0_dp)
            call add(leaf_source, [source(4), test(1)], zero(:2), -1.0_dp)
        case (2)
            call add(leaf_mass_momentum, [u(1), ud(1), u(2:4), ud(2:4), test(5:7)], &
                [du(1), dud(1), du(2:4), dud(2:4), zero(:3)], 1.0_dp)
            call add(leaf_momentum_advection, [u(1:4), test(5:7), u(16:24), &
                test(8:16)], [du(1:4), zero(:3), du(16:24), zero(:9)], 1.0_dp)
            call add(leaf_momentum_mass_source, [source(4), u(2:4), test(5:7)], &
                [0.0_dp, du(2:4), zero(:3)], 1.0_dp)
            call add(leaf_momentum_density_diffusion, &
                [p%d_parallel, p%d_perpendicular, p%epsilon_b, u(25:27), &
                u(13:15), u(2:4), test(5:7), u(16:24), test(8:16)], &
                [zero(:3), du(25:27), du(13:15), du(2:4), zero(:3), &
                du(16:24), zero(:9)], 1.0_dp)
            call vrmhd_leaf(leaf_pressure, [u(1), u(11:12), p%mi], &
                [du(1), du(11:12), 0.0_dp], pressure, dpressure)
            call add(leaf_pressure_work, [pressure, divw], [dpressure, 0.0_dp], &
                1.0_dp)
            call add(leaf_lorentz, [u(8:10), u(25:27), test(5:7)], &
                [du(8:10), du(25:27), zero(:3)], -1.0_dp)
            call add(leaf_viscosity, [nu, p%nu_parallel, p%epsilon_b, u(25:27), &
                u(16:24), test(8:16)], &
                [dnu, zero(:2), du(25:27), du(16:24), zero(:9)], -1.0_dp)
        case (3)
            call add(leaf_mass_induction, [ud(5:7), test(5:7)], [dud(5:7), &
                zero(:3)], 1.0_dp)
            call add(leaf_induction_ideal, [u(2:4), u(25:27), test(5:7)], &
                [du(2:4), du(25:27), zero(:3)], -1.0_dp)
            call add(leaf_induction_resistive, [eta, u(8:10), source(:3), &
                test(5:7)], &
                [deta, du(8:10), zero(:6)], -1.0_dp)
        case (4)
            call add(leaf_ampere, [p%mu0, u(8:10), test(5:7), u(25:27), &
                test(17:19)], &
                [0.0_dp, du(8:10), zero(:3), du(25:27), zero(:3)], 1.0_dp)
        case (5, 6)
            species = row - 4; it = 10 + species
            if (p%flux_limited(species)) then
                call vrmhd_leaf(leaf_parallel_gradient, [u(25:27), &
                    u(25 + 3*species:27 + 3*species), p%epsilon_b], &
                    [du(25:27), du(25 + 3*species:27 + 3*species), &
                    0.0_dp], bg, dbg)
                mass = p%electron_mass
                if (species == 2) mass = p%mi
                call vrmhd_leaf(leaf_chi_flux_limited, [p%chi_parallel0(species), &
                    u(it), p%t0, &
                    mass, &
                    p%flux_limit(species), &
                    bg, &
                    p%gradient_floor(species)], &
                    [0.0_dp, du(it), zero(:3), dbg, 0.0_dp], chi, dchi)
            else
                call vrmhd_leaf(leaf_chi_parallel, [p%chi_parallel0(species), &
                    u(it), p%t0, p%chi_parallel_exponent], &
                    [0.0_dp, du(it), 0.0_dp, 0.0_dp], chi, dchi)
            end if
            call add(leaf_mass_temperature, [u(1), ud(1), u(it), ud(it), p%mi, &
                test(1)], &
                [du(1), dud(1), du(it), dud(it), zero(:2)], 1.0_dp)
            call add(leaf_temperature_advection, [u(1), p%mi, u(it), u(2:4), &
                test(2:4)], &
                [du(1), 0.0_dp, du(it), du(2:4), zero(:3)], 1.0_dp)
            call add(leaf_compression, [u(1), p%mi, u(it), divv, test(1)], &
                [du(1), 0.0_dp, du(it), ddivv, 0.0_dp], 1.0_dp)
            call add(leaf_conduction, [u(1), p%mi, chi, &
                p%chi_perpendicular(species), p%epsilon_b, &
                u(25:27), u(25 + 3*species:27 + 3*species), &
                test(2:4)], &
                [du(1), 0.0_dp, dchi, zero(:2), du(25:27), &
                du(25 + 3*species:27 + 3*species), zero(:3)], 1.0_dp)
            if (species == 1) then
                call add(leaf_ohmic, [eta, u(8:10), source(:3), test(1)], &
                    [deta, du(8:10), zero(:4)], -p%ohmic_fraction)
            else
                call add(leaf_viscous_heat, [nu, p%nu_parallel, p%epsilon_b, &
                    u(25:27), u(16:24), test(1)], &
                    [dnu, zero(:2), du(25:27), du(16:24), 0.0_dp], -1.0_dp)
            end if
            call add(leaf_equipartition, [u(1), p%mi, p%equilibration_rate, &
                u(11:12), test(1)], &
                [du(1), zero(:2), du(11:12), 0.0_dp], real(2*species - 3, dp))
            call add(leaf_source, [source(4 + species), test(1)], zero(:2), -1.0_dp)
        end select
    contains
        subroutine add(id, x, dx, sign)
            integer, intent(in) :: id
            real(dp), intent(in) :: x(:), dx(:), sign
            call vrmhd_leaf(id, x, dx, h, dh)
            value = value + sign*h; tangent = tangent + sign*dh
        end subroutine
    end subroutine
end module milca_vrmhd_point
