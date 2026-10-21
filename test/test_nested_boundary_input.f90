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
    integer::nfp, mpol, ntor, k, m, n, u, status, family
    logical::lfreeb, lasym
    character(512)::path, arg
    namelist /indata/ nfp, mpol, ntor, lfreeb, lasym, rbc, rbs, zbc, zbs
    rbc = 0; rbs = 0; zbc = 0; zbs = 0
    call get_command_argument(1, path)
    call get_command_argument(2, arg); read (arg, *) family
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
end program
