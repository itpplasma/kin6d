!> NTV flux kernel of the Legendre reference (#183).
!!   jvp          tangents in nu and omega_E against centred differences;
!!   kernel FILE  kernels of a surface at the speeds listed in FILE (the
!!                stage-1 gate driver feeds this mode).
program test_dk_ntv
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_dk_ntv, only: ntv_surface_t, ntv_kernel_t, ntv_kernel, NTV_OK
    implicit none
    character(256) :: check, path
    logical :: ok

    call get_command_argument(1, check)
    select case (trim(check))
    case ('jvp'); ok = jvp()
    case ('kernel')
        call get_command_argument(2, path)
        call kernel_table(trim(path))
        ok = .true.
    case default
        print '(a)', 'unknown check '//trim(check)
        ok = .false.
    end select
    if (.not. ok) error stop 1

contains

    !> One (m, n) = (6, 3) harmonic on eps = 0.05, q = 2.2 with collisions
    !> and precession of the stage-1 order.
    function sample() result(sf)
        type(ntv_surface_t) :: sf
        sf%eps = 0.05_dp
        sf%q = 2.2_dp
        sf%k = 1/(sf%q*6.2_dp)
        sf%offset = 0.01_dp
        sf%n = 3
        sf%m = [6]
        sf%c = [(1e-3_dp, 2e-4_dp)]
    end function sample

    logical function jvp() result(ok)
        real(dp), parameter :: v = 1.0_dp, nu = 2e-3_dp, omega = 5e-3_dp, pre = 2e-2_dp
        real(dp) :: h, fd
        type(ntv_kernel_t) :: k0, up, down
        integer :: status
        call ntv_kernel(sample(), v, nu, omega, pre, 33, 80, k0, status)
        if (status /= NTV_OK) error stop 'kernel failed'
        h = 1e-6_dp*nu
        call ntv_kernel(sample(), v, nu + h, omega, pre, 33, 80, up, status)
        call ntv_kernel(sample(), v, nu - h, omega, pre, 33, 80, down, status)
        fd = (up%r - down%r)/(2*h)
        print '(a,2es14.6,es9.1)', 'dR/dnu      ', k0%dr_dnu, fd, abs(k0%dr_dnu/fd - 1)
        ok = abs(k0%dr_dnu/fd - 1) < 1e-6_dp
        h = 1e-6_dp*omega
        call ntv_kernel(sample(), v, nu, omega + h, pre, 33, 80, up, status)
        call ntv_kernel(sample(), v, nu, omega - h, pre, 33, 80, down, status)
        fd = (up%r - down%r)/(2*h)
        print '(a,2es14.6,es9.1)', 'dR/domega_E ', k0%dr_domega, fd, abs(k0%dr_domega/fd - 1)
        ok = ok .and. abs(k0%dr_domega/fd - 1) < 1e-6_dp
    end function jvp

    !> FILE: eps k q offset n nm / nm lines m Re(c) Im(c) / n_theta n_leg /
    !> nx / nx lines v nu_d omega_E pre. Prints per node R, the theta and
    !> Legendre indicators, dR/dnu_D and dR/domega_E.
    subroutine kernel_table(file)
        character(*), intent(in) :: file
        type(ntv_surface_t) :: sf
        type(ntv_kernel_t) :: k
        integer :: u, nm, j, n_theta, n_leg, nx, status
        real(dp) :: re, im, v, nu_d, omega, pre
        open (newunit=u, file=file, status='old', action='read')
        read (u, *) sf%eps, sf%k, sf%q, sf%offset, sf%n, nm
        allocate (sf%m(nm), sf%c(nm))
        do j = 1, nm
            read (u, *) sf%m(j), re, im
            sf%c(j) = cmplx(re, im, dp)
        end do
        read (u, *) n_theta, n_leg
        read (u, *) nx
        do j = 1, nx
            read (u, *) v, nu_d, omega, pre
            call ntv_kernel(sf, v, nu_d, omega, pre, n_theta, n_leg, k, status)
            if (status /= NTV_OK) error stop 'kernel failed'
            print '(5es24.15)', k%r, k%indicator_theta, k%indicator_leg, k%dr_dnu, k%dr_domega
        end do
        close (u)
    end subroutine kernel_table

end program test_dk_ntv
