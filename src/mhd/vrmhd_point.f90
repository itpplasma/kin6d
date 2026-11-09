module milca_vrmhd_point
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use milca_vrmhd_parameters, only: vrmhd_parameters_t
    use milca_vrmhd_dispatch, only: vrmhd_leaf, vrmhd_leaf_vjp
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
    !> Cost of one vrmhd_point_matrix call in vrmhd_point evaluations: one
    !> reverse sweep, about two evaluations, per test slot.
    integer, parameter, public :: VRMHD_POINT_MATRIX_CALLS = 2*sum(ntest)
    integer, parameter :: test_slots(12, 6) = reshape([ &
        1, 2, 3, 4, 0, 0, 0, 0, 0, 0, 0, 0, &
        5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, &
        5, 6, 7, 17, 18, 19, 0, 0, 0, 0, 0, 0, &
        5, 6, 7, 17, 18, 19, 0, 0, 0, 0, 0, 0, &
        1, 2, 3, 4, 0, 0, 0, 0, 0, 0, 0, 0, &
        1, 2, 3, 4, 0, 0, 0, 0, 0, 0, 0, 0], [12, 6])

contains
    ! Pointwise bilinear form of the JVP: kmat(:, j, row) is the tangent
    ! against unit tests for unit trial data j (33 values, 12 rates). The
    ! residual is linear in the test, so one reverse sweep per test slot gives
    ! every trial column at once.
    subroutine vrmhd_point_matrix(p, v, vd, source, shift, kmat)
        type(vrmhd_parameters_t), intent(in) :: p
        real(dp), intent(in) :: v(33), vd(12), source(6), shift
        real(dp), intent(out) :: kmat(19, 45, 6)
        real(dp) :: test(19), gu(33), gud(12)
        integer :: row, it, slot
        kmat = 0
        do row = 1, 6
            do it = 1, ntest(row)
                slot = test_slots(it, row)
                test = 0; test(slot) = 1
                call point_gradient(p, row, v, vd, test, source, gu, gud)
                kmat(slot, :33, row) = gu
                kmat(slot, 34:, row) = shift*gud
            end do
        end do
    end subroutine

    subroutine point_gradient(p, row, u, ud, test, source, gu, gud)
        !! Gradient of the vrmhd_point value with respect to u and ud: the
        !! reverse sweep of the same leaf composition, with generated VJPs.
        type(vrmhd_parameters_t), intent(in) :: p
        integer, intent(in) :: row
        real(dp), intent(in) :: u(33), ud(12), test(19), source(6)
        real(dp), intent(out) :: gu(33), gud(12)
        integer, parameter :: KETA = 101, KNU = 102, KDIVV = 103, KPRESS = 104, &
            KCHI = 105
        real(dp) :: eta, nu, chi, pressure, divv, divw, bg, mass, unused
        real(dp) :: geta, gnu, gchi, gpress, gdivv, gbg, gx(33)
        real(dp) :: zero(24)
        integer :: species, it, g1, i
        integer :: m(33)
        zero = 0; gu = 0; gud = 0
        geta = 0; gnu = 0; gchi = 0; gpress = 0; gdivv = 0; gbg = 0
        if (row == 3 .or. row == 5) call vrmhd_leaf(leaf_eta, [p%eta0, u(11), p%t0, &
            p%eta_exponent], zero(:4), eta, unused)
        if (row == 2 .or. row == 6) call vrmhd_leaf(leaf_nu, [p%nu0, u(12), p%t0, &
            p%nu_exponent], zero(:4), nu, unused)
        if (row >= 5) call vrmhd_leaf(leaf_trace, [u(16), u(20), u(24)], zero(:3), &
            divv, unused)
        if (row == 2) call vrmhd_leaf(leaf_trace, [test(8), test(12), test(16)], &
            zero(:3), divw, unused)
        select case (row)
        case (1)
            call add(leaf_mass_density, [ud(1), test(1)], [-1, 0], 1.0_dp)
            call add(leaf_density_advection, [u(1:4), test(2:4)], &
                [1, 2, 3, 4, 0, 0, 0], 1.0_dp)
            call add(leaf_transport, [p%d_parallel, p%d_perpendicular, &
                p%epsilon_b, u(25:27), u(13:15), test(2:4)], &
                [0, 0, 0, 25, 26, 27, 13, 14, 15, 0, 0, 0], 1.0_dp)
        case (2)
            call add(leaf_mass_momentum, [u(1), ud(1), u(2:4), ud(2:4), test(5:7)], &
                [1, -1, 2, 3, 4, -2, -3, -4, 0, 0, 0], 1.0_dp)
            m = 0; m(:4) = [1, 2, 3, 4]; m(8:16) = [(i, i=16, 24)]
            call add(leaf_momentum_advection, [u(1:4), test(5:7), u(16:24), &
                test(8:16)], m(:25), 1.0_dp)
            call add(leaf_momentum_mass_source, [source(4), u(2:4), test(5:7)], &
                [0, 2, 3, 4, 0, 0, 0], 1.0_dp)
            m = 0; m(4:12) = [25, 26, 27, 13, 14, 15, 2, 3, 4]
            m(16:24) = [(i, i=16, 24)]
            call add(leaf_momentum_density_diffusion, &
                [p%d_parallel, p%d_perpendicular, p%epsilon_b, u(25:27), &
                u(13:15), u(2:4), test(5:7), u(16:24), test(8:16)], m(:33), 1.0_dp)
            call vrmhd_leaf(leaf_pressure, [u(1), u(11:12), p%mi], zero(:4), &
                pressure, unused)
            call add(leaf_pressure_work, [pressure, divw], [KPRESS, 0], 1.0_dp)
            call add(leaf_lorentz, [u(8:10), u(25:27), test(5:7)], &
                [8, 9, 10, 25, 26, 27, 0, 0, 0], -1.0_dp)
            m = 0; m(1) = KNU; m(4:15) = [25, 26, 27, (i, i=16, 24)]
            call add(leaf_viscosity, [nu, p%nu_parallel, p%epsilon_b, u(25:27), &
                u(16:24), test(8:16)], m(:24), -1.0_dp)
            call vrmhd_leaf_vjp(leaf_pressure, [u(1), u(11:12), p%mi], gpress, gx)
            gu(1) = gu(1) + gx(1); gu(11:12) = gu(11:12) + gx(2:3)
        case (3)
            call add(leaf_mass_induction, [ud(5:7), test(5:7)], &
                [-5, -6, -7, 0, 0, 0], 1.0_dp)
            call add(leaf_induction_ideal, [u(2:4), u(25:27), test(5:7)], &
                [2, 3, 4, 25, 26, 27, 0, 0, 0], -1.0_dp)
            call add(leaf_induction_resistive, [eta, u(8:10), source(:3), &
                test(5:7)], [KETA, 8, 9, 10, 0, 0, 0, 0, 0, 0], -1.0_dp)
        case (4)
            call add(leaf_ampere, [p%mu0, u(8:10), test(5:7), u(25:27), &
                test(17:19)], [0, 8, 9, 10, 0, 0, 0, 25, 26, 27, 0, 0, 0], 1.0_dp)
        case (5, 6)
            species = row - 4; it = 10 + species; g1 = 25 + 3*species
            if (p%flux_limited(species)) then
                call vrmhd_leaf(leaf_parallel_gradient, [u(25:27), &
                    u(g1:g1 + 2), p%epsilon_b], zero(:7), bg, unused)
                mass = p%electron_mass
                if (species == 2) mass = p%mi
                call vrmhd_leaf(leaf_chi_flux_limited, [p%chi_parallel0(species), &
                    u(it), p%t0, mass, p%flux_limit(species), bg, &
                    p%gradient_floor(species)], zero(:7), chi, unused)
            else
                call vrmhd_leaf(leaf_chi_parallel, [p%chi_parallel0(species), &
                    u(it), p%t0, p%chi_parallel_exponent], zero(:4), chi, unused)
            end if
            call add(leaf_mass_temperature, [u(1), ud(1), u(it), ud(it), p%mi, &
                test(1)], [1, -1, it, -it, 0, 0], 1.0_dp)
            call add(leaf_temperature_advection, [u(1), p%mi, u(it), u(2:4), &
                test(2:4)], [1, 0, it, 2, 3, 4, 0, 0, 0], 1.0_dp)
            call add(leaf_compression, [u(1), p%mi, u(it), divv, test(1)], &
                [1, 0, it, KDIVV, 0], 1.0_dp)
            call add(leaf_conduction, [u(1), p%mi, chi, &
                p%chi_perpendicular(species), p%epsilon_b, u(25:27), &
                u(g1:g1 + 2), test(2:4)], &
                [1, 0, KCHI, 0, 0, 25, 26, 27, g1, g1 + 1, g1 + 2, 0, 0, 0], 1.0_dp)
            if (species == 1) then
                call add(leaf_ohmic, [eta, u(8:10), source(:3), test(1)], &
                    [KETA, 8, 9, 10, 0, 0, 0, 0], -p%ohmic_fraction)
            else
                m = 0; m(1) = KNU; m(4:15) = [25, 26, 27, (i, i=16, 24)]
                call add(leaf_viscous_heat, [nu, p%nu_parallel, p%epsilon_b, &
                    u(25:27), u(16:24), test(1)], m(:16), -1.0_dp)
            end if
            call add(leaf_equipartition, [u(1), p%mi, p%equilibration_rate, &
                u(11:12), test(1)], [1, 0, 0, 11, 12, 0], real(2*species - 3, dp))
            gu(16) = gu(16) + gdivv; gu(20) = gu(20) + gdivv
            gu(24) = gu(24) + gdivv
            if (p%flux_limited(species)) then
                call vrmhd_leaf_vjp(leaf_chi_flux_limited, &
                    [p%chi_parallel0(species), u(it), p%t0, mass, &
                    p%flux_limit(species), bg, p%gradient_floor(species)], &
                    gchi, gx)
                gu(it) = gu(it) + gx(2); gbg = gx(6)
                call vrmhd_leaf_vjp(leaf_parallel_gradient, [u(25:27), &
                    u(g1:g1 + 2), p%epsilon_b], gbg, gx)
                gu(25:27) = gu(25:27) + gx(1:3)
                gu(g1:g1 + 2) = gu(g1:g1 + 2) + gx(4:6)
            else
                call vrmhd_leaf_vjp(leaf_chi_parallel, &
                    [p%chi_parallel0(species), u(it), p%t0, &
                    p%chi_parallel_exponent], gchi, gx)
                gu(it) = gu(it) + gx(2)
            end if
        end select
        if (row == 3 .or. row == 5) then
            call vrmhd_leaf_vjp(leaf_eta, [p%eta0, u(11), p%t0, p%eta_exponent], &
                geta, gx)
            gu(11) = gu(11) + gx(2)
        end if
        if (row == 2 .or. row == 6) then
            call vrmhd_leaf_vjp(leaf_nu, [p%nu0, u(12), p%t0, p%nu_exponent], gnu, gx)
            gu(12) = gu(12) + gx(2)
        end if
    contains
        subroutine add(id, x, map, sign)
            !! map(k) > 0: u(map(k)); < 0: ud(-map(k)); intermediate codes; 0: data.
            integer, intent(in) :: id, map(:)
            real(dp), intent(in) :: x(:), sign
            real(dp) :: g(size(x))
            integer :: k
            call vrmhd_leaf_vjp(id, x, sign, g)
            do k = 1, size(map)
                select case (map(k))
                case (1:33)
                    gu(map(k)) = gu(map(k)) + g(k)
                case (-12:-1)
                    gud(-map(k)) = gud(-map(k)) + g(k)
                case (KETA)
                    geta = geta + g(k)
                case (KNU)
                    gnu = gnu + g(k)
                case (KDIVV)
                    gdivv = gdivv + g(k)
                case (KPRESS)
                    gpress = gpress + g(k)
                case (KCHI)
                    gchi = gchi + g(k)
                end select
            end do
        end subroutine
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
