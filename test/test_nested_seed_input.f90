program test_nested_seed_input
    ! Seed input from the exact nested references: boundary reconstruction
    ! off the projection grid, and profile fits against independent oracle
    ! profile evaluations at points that are neither fit nor check nodes.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_nested_oracles, only: nested_oracle_t, NESTED_INTEGER, NESTED_SHEARED
    use kin6d_nested_oracle_leaves, only: nested_setup
    use kin6d_ideal_seed_3d, only: toroidal_profile_t
    use kin6d_nested_seed_input, only: nested_seed_boundary, nested_seed_profile, &
        NESTED_SEED_OK, NESTED_SEED_FIT, NESTED_SEED_BAD_INPUT
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp)
    type(nested_oracle_t) :: o
    type(toroidal_profile_t) :: profile
    real(dp) :: flux, exact_flux, fit_error, error(2), limit, extent, pa1, pa2
    real(dp) :: s(7), psi(7), phi(7), iota(7), pressure(7), fi, fp, unused, worst
    integer :: status, ierr, k, point_status(7)
    ! Integer family: exact constant transform and linear pressure.
    o%family = NESTED_INTEGER; o%eps = 0.5_dp; o%delta = 1.0_dp/64
    call nested_seed_profile(o, 4, 1e-9_dp, profile, flux, fit_error, status)
    if (status /= NESTED_SEED_OK) error stop 'integer profile status'
    call nested_setup(o%eps, o%tau, o%delta, o%lam, o%edge_pressure, limit, extent, &
        pa1, pa2)
    if (size(profile%iota) /= 1 .or. abs(profile%iota(1) - 2) > 1e-14_dp) &
        error stop 'integer iota must be the constant 2'
    if (size(profile%pressure) /= 2 .or. abs(profile%pressure(1) - pa1) > 1e-14_dp &
        .or. abs(profile%pressure(1) + profile%pressure(2) - o%edge_pressure) &
        > 1e-14_dp) error stop 'integer pressure must be linear in s'
    if (fit_error > 1e-13_dp) error stop 'integer profile check'
    call boundary_error(o, 8, 4, error(1))
    call boundary_error(o, 16, 8, error(2))
    print '(a,2es11.3)', 'integer boundary error m8n4 m16n8:', error
    if (error(2) > 1e-6_dp .or. error(2) > error(1)) error stop 'integer boundary'
    ! Sheared family (examples/nested-oracles/landreman-sheared.nml).
    o%family = NESTED_SHEARED; o%eps = 1.2_dp; o%delta = 0.14_dp
    call o%flux(o%delta, 1e-12_dp, exact_flux, ierr)
    if (ierr /= 0) error stop 'sheared reference flux'
    call nested_seed_profile(o, 2, 1e-9_dp, profile, flux, fit_error, status)
    if (status /= NESTED_SEED_FIT) error stop 'a quadratic sheared fit must be rejected'
    call nested_seed_profile(o, 12, 1e-9_dp, profile, flux, fit_error, status)
    print '(a,es11.3)', 'sheared profile fit error degree 12:', fit_error
    if (status /= NESTED_SEED_OK) error stop 'sheared profile status'
    if (abs(flux/exact_flux - 1) > 1e-11_dp) error stop 'sheared toroidal flux'
    s = [0.013_dp, 0.11_dp, 0.29_dp, 0.47_dp, 0.61_dp, 0.83_dp, 0.977_dp]
    call o%profiles(s, 1e-11_dp, psi, phi, iota, pressure, point_status)
    if (any(point_status /= 0)) error stop 'independent profile evaluation'
    worst = 0
    do k = 1, size(s)
        call profile%evaluate(s(k), fi, fp, unused)
        worst = max(worst, abs(fi/iota(k) - 1), abs(fp/pressure(k) - 1))
    end do
    print '(a,es11.3,a,2f10.6)', 'sheared off-node profile error:', worst, &
        ' iota axis/edge:', profile%iota(1), sum(profile%iota)
    if (worst > 1e-9_dp) error stop 'sheared off-node profile'
    call boundary_error(o, 8, 4, error(1))
    call boundary_error(o, 16, 8, error(2))
    print '(a,2es11.3)', 'sheared boundary error m8n4 m16n8:', error
    if (error(2) > 1e-4_dp .or. error(2) > error(1)/8) error stop 'sheared boundary'
    call nested_seed_profile(o, 0, 1e-9_dp, profile, flux, fit_error, status)
    if (status /= NESTED_SEED_BAD_INPUT) error stop 'degree 0 must be rejected'
    print '(a)', 'PASS nested seed input'
contains
    subroutine boundary_error(o, mpol, ntor, error)
        ! Independent reconstruction in the VMEC convention, off-grid angles,
        ! relative to the exact boundary's minor extent.
        type(nested_oracle_t), intent(in) :: o
        integer, intent(in) :: mpol, ntor
        real(dp), intent(out) :: error
        real(dp), allocatable :: c(:, :, :)
        real(dp) :: theta, phi, x(3), r, z, q, rmin, rmax
        integer :: k, m, n, st
        call nested_seed_boundary(o, mpol, ntor, c, st)
        if (st /= NESTED_SEED_OK) error stop 'boundary projection status'
        error = 0; rmin = huge(1.0_dp); rmax = 0
        do k = 1, 211
            theta = 2*pi*(k - 0.31_dp)/211
            phi = 2*pi*mod(37*k, 211)/211 + 0.017_dp
            call o%boundary(theta, phi, x, st)
            if (st /= 0) error stop 'boundary root'
            r = 0; z = 0
            do m = 0, mpol
                do n = -ntor, ntor
                    q = m*theta - 2*n*phi
                    r = r + c(1, n, m)*cos(q) + c(2, n, m)*sin(q)
                    z = z + c(3, n, m)*cos(q) + c(4, n, m)*sin(q)
                end do
            end do
            error = max(error, hypot(r - hypot(x(1), x(2)), z - x(3)))
            rmin = min(rmin, hypot(x(1), x(2))); rmax = max(rmax, hypot(x(1), x(2)))
        end do
        error = error/(rmax - rmin)
    end subroutine
end program
