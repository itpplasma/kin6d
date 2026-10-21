! Standalone recipe for the public regression inputs; no machine data.
! Compile with gfortran -o generate_gs_synthetic generate_gs_synthetic.f90;
! run in test/fixtures/gs_synthetic. Reference outputs come from gs_phase1.
program generate_gs_synthetic
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp), mu0 = 4*pi*1.0e-7_dp
    real(dp), parameter :: r0 = 3, a = 0.6_dp, c = 0.25_dp, f = 6
    integer, parameter :: n = 512, nq = 16384, np = 12
    real(dp) :: t, r, z, dz, integral, current
    integer :: u, j

    open (newunit=u, file='boundary.dat', status='replace')
    write (u, *) n
    do j = 0, n - 1
        t = 2*pi*j/n
        r = sqrt(r0*r0 + 2*r0*a*cos(t))
        z = r0*a*sin(t)/r
        write (u, '(2es25.16e3)') r, z
    end do
    close (u)
    ! Green's theorem: integral R dR dZ = (1/2) integral R**2 dZ.
    ! Smooth periodic midpoint quadrature of the exact boundary, independent
    ! of the production mesh, supplies the prescribed current in amperes.
    integral = 0
    do j = 0, nq - 1
        t = 2*pi*(j + 0.5_dp)/nq
        r = sqrt(r0*r0 + 2*r0*a*cos(t))
        dz = r0*a*cos(t)/r + r0*r0*a*a*sin(t)**2/r**3
        integral = integral + r*r*dz/2
    end do
    integral = integral*2*pi/nq
    current = -16*c*integral/mu0
    open (newunit=u, file='input.nml', status='replace')
    write (u, '(a)') '&gs_phase1', &
        " boundary_file='boundary.dat', profile_file='source.dat',", &
        ' degree=3, boundary_nodes=64, nonlinear_tolerance=1e-12,', &
        ' quadrature_tolerance=1e-5, radial_points=65, angular_points=64,', &
        ' solve_wall_limit_s=240,'
    write (u, '(a,es25.16e3,a)') ' total_current=', current, ','
    write (u, '(a)') '/'
    close (u)
    open (newunit=u, file='source.dat', status='replace')
    write (u, *) n + 1, 1
    write (u, '(3es25.16e3)') 0.0_dp, f, -4*c*r0*r0*a*a
    write (u, '(*(es25.16e3,1x))') (real(j, dp)/n, j=0, n)
    do j = 1, n
        write (u, '(8es25.16e3)') -16*c/mu0, 0.0_dp, 0.0_dp, 0.0_dp, &
            0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp
    end do
    close (u)
    open (newunit=u, file='points.dat', status='replace')
    write (u, *) np
    do j = 0, np - 1
        t = 2*pi*j/np
        r = sqrt(r0*r0 + r0*a*cos(t))
        z = r0*a*sin(t)/(2*r)
        write (u, '(2es25.16e3)') r, z
    end do
    close (u)
end program
