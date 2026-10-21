program test_nested_boundary_input
    ! Independent reconstruction in the external VMEC Fourier convention.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_nested_oracles, only: nested_oracle_t
    implicit none
    integer, parameter::maxmode = 32
    type(nested_oracle_t)::o
    real(dp)::rbc(-maxmode:maxmode, 0:maxmode), rbs(-maxmode:maxmode, 0:maxmode)
    real(dp)::zbc(-maxmode:maxmode, 0:maxmode), zbs(-maxmode:maxmode, 0:maxmode)
    real(dp)::theta, phi, x(3), r, z, phase, error, pi
    real(dp)::phiedge, pres_scale, am(2), ai(1), ftol_array(3)
    real(dp)::xp(3), xm(3), dr(3), dt(3), v(20), jac, rho, step, integral
    integer::nfp, mpol, ntor, k, m, n, u, status, family, ncurr
    integer::ns_array(3), niter_array(3), ir, it
    logical::lfreeb, lasym
    character(512)::path, arg, mode
    character(32)::pmass_type, piota_type
    namelist /indata/ nfp, mpol, ntor, lfreeb, lasym, rbc, rbs, zbc, zbs, &
        phiedge, pres_scale, am, ai, ncurr, pmass_type, piota_type, &
        ns_array, ftol_array, niter_array
    rbc = 0; rbs = 0; zbc = 0; zbs = 0
    call get_command_argument(1, path)
    call get_command_argument(2, arg); read (arg, *) family
    call get_command_argument(3, mode)
    o%family = family; o%eps = 0.5_dp; o%tau = 0.4_dp
    if (family == 2) then
        o%eps = 1.2_dp; o%delta = 0.14_dp
    end if
    open (newunit=u, file=trim(path), status='old')
    read (u, nml=indata)
    close (u)
    if (nfp /= 2 .or. .not. lasym .or. lfreeb) error stop 'boundary symmetry metadata'
    pi = acos(-1.0_dp); error = 0
    do k = 1, 157
        theta = 2*pi*k/157; phi = 2*pi*mod(11*k, 157)/157
        call o%boundary(theta, phi, x, status)
        if (status /= 0) error stop 'boundary root'
        r = 0; z = 0
        do m = 0, mpol - 1
            do n = -ntor, ntor
                phase = m*theta - n*nfp*phi
                r = r + rbc(n, m)*cos(phase) + rbs(n, m)*sin(phase)
                z = z + zbc(n, m)*cos(phase) + zbs(n, m)*sin(phase)
            end do
        end do
        error = max(error, norm2([r - sqrt(x(1)**2 + x(2)**2), z - x(3)]))
    end do
    if (error > 0.001_dp) error stop 'external boundary convention error'
if (maxval(abs(rbs)) + maxval(abs(zbc)) < 1e-5_dp) error stop 'asymmetry lost in export'
    print '(a,es12.3)', '[V] independent VMEC boundary reconstruction max error ', error
    if (trim(mode) == 'matched') then
        if (family /= 1 .or. ncurr /= 0) error stop 'matched input model'
        if (trim(pmass_type) /= 'power_series' .or. &
            trim(piota_type) /= 'power_series') error stop 'matched profile type'
        if (abs(ai(1) - 2) > 1e-14_dp) error stop 'integer transform'
        if (abs(pres_scale*4*pi*1e-7_dp - 1) > 1e-14_dp) error stop 'SI pressure units'
        ! Independent physical flux integral on the material t=0 cross-section.
        ! Finite-difference tangents deliberately do not use the generated jacobian.
        integral = 0; step = 1e-5_dp
        do ir = 1, 16
            rho = sqrt(o%delta)*(ir - 0.5_dp)/16
            do it = 0, 63
                theta = 2*pi*(it + 0.5_dp)/64
                call o%surface(rho + step, theta, 0.0_dp, xp, jac)
                call o%surface(rho - step, theta, 0.0_dp, xm, jac)
                dr = (xp - xm)/(2*step)
                call o%surface(rho, theta + step, 0.0_dp, xp, jac)
                call o%surface(rho, theta - step, 0.0_dp, xm, jac)
                dt = (xp - xm)/(2*step)
                call o%surface(rho, theta, 0.0_dp, x, jac)
                call o%point(x, v, status)
                if (status /= 0) error stop 'flux quadrature domain'
                integral = integral + dot_product(v(3:5), &
                    [dr(2)*dt(3) - dr(3)*dt(2), dr(3)*dt(1) - dr(1)*dt(3), &
                     dr(1)*dt(2) - dr(2)*dt(1)])
                if (abs(am(1) + am(2)*rho**2/o%delta - v(2)) > 1e-13_dp) &
                    error stop 'matched pressure differs from Cartesian oracle'
            end do
        end do
        integral = abs(integral)*sqrt(o%delta)/16*(2*pi/64)
        if (abs(integral/phiedge - 1) > 1e-7_dp) error stop 'matched toroidal flux'
        print '(a,es12.3)', '[V] independent flux integral relative error ', &
            abs(integral/phiedge - 1)
    end if
end program
