!> Axisymmetric drift-kinetic solve on the energy grid (#216, #184).
!!   rigid_shift   with the momentum-conserving operator the D31 drive has
!!                 the exact solution h = x xi/b, so D31 = sum w x^2/3; the
!!                 test-particle operator misses it;
program test_dk_energy
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_dk_energy
    use kin6d_dk_collop, only: collop_species_t, collop_grid_t, collop_make_grid, &
        COLLOP_FULL, COLLOP_TEST_PARTICLE, COLLOP_LORENTZ, collop_qe
    implicit none
    real(dp), parameter :: pi = 4*atan(1.0_dp)
    character(64) :: check
    logical :: ok
    call get_command_argument(1, check)
    select case (trim(check))
    case ('rigid_shift'); ok = rigid_shift()
    case default
        print '(a)', 'unknown check '//trim(check)
        ok = .false.
    end select
    if (.not. ok) error stop 1
contains
    pure function field(n, eps) result(b)
        integer, intent(in) :: n
        real(dp), intent(in) :: eps
        real(dp) :: b(n)
        integer :: j
        b = [(1/(1 + eps*cos(2*pi*(j - 1)/n)), j=1, n)]
    end function field

    !> Deuterons at 1 keV with density n [m^-3].
    function deuterons(n) result(s)
        real(dp), intent(in) :: n
        type(collop_species_t) :: s
        s = collop_species_t(3.3436e-27_dp, 1.0_dp, n, 1e3_dp*collop_qe)
    end function deuterons

    logical function rigid_shift() result(ok)
        type(collop_grid_t) :: grid
        type(energy_axi_t) :: full, tp
        integer :: status
        call collop_make_grid(6.0_dp, 4, 5, grid, status)
        call energy_axi_coefficients(field(33, 0.1_dp), 1.0_dp, deuterons(1e20_dp), 17.0_dp, &
            grid, 24, COLLOP_FULL, full, status)
        if (status /= ENERGY_OK) error stop 'full solve failed'
        call energy_axi_coefficients(field(33, 0.1_dp), 1.0_dp, deuterons(1e20_dp), 17.0_dp, &
            grid, 24, COLLOP_LORENTZ, tp, status)
        if (status /= ENERGY_OK) error stop 'test-particle solve failed'
        print '(a,3es16.8)', 'D31 full, test particle, rigid shift:', full%d31, tp%d31, &
            full%d31_rigid
        print '(a,2f12.6,a,2es9.1)', 'k full, test particle:', full%k, tp%k, '  residuals', &
            full%residual
        ok = abs(full%d31/full%d31_rigid - 1) < 1e-8_dp .and. all(full%residual < 1e-8_dp)
    end function rigid_shift
end program test_dk_energy
