program coil_modes
    !! K-C1 driver: toroidal Fourier components of the 3D coil field at
    !! meridional points. Usage: coil_modes coils.tab currents.dat points.dat out.dat
    !! [nphi] [nmax]. points.dat: "np" then np lines "R Z". Output (one line per
    !! point and mode n = 0..nmax): n R Z, then Re/Im of B_R, B_phi, B_Z of the
    !! two-sided component B_n (physical n>0 harmonic = 2 Re[B_n e^{i n phi}]), in T.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_coils_3d
    implicit none
    type(coil_set_t) :: set
    character(len=512) :: arg
    real(dp), allocatable :: w(:), rz(:, :)
    complex(dp), allocatable :: bn(:, :, :)
    integer, allocatable :: modes(:)
    integer :: nphi = 512, nmax = 12, np, unit, i, im, st

    if (command_argument_count() < 4) &
        error stop 'usage: coil_modes coils.tab currents points out [nphi nmax]'
    call get_command_argument(1, arg)
    call coil_set_read(set, trim(arg), st)
    if (st /= 0) error stop 'cannot read coil table'
    allocate (w(set%ncoil))
    call get_command_argument(2, arg)
    call coil_read_currents(trim(arg), set%ncoil, w, st)
    if (st /= 0) error stop 'cannot read currents'
    if (command_argument_count() >= 5) then
        call get_command_argument(5, arg); read (arg, *) nphi
    end if
    if (command_argument_count() >= 6) then
        call get_command_argument(6, arg); read (arg, *) nmax
    end if
    call get_command_argument(3, arg)
    open (newunit=unit, file=trim(arg), status='old', action='read')
    read (unit, *) np
    allocate (rz(2, np))
    do i = 1, np
        read (unit, *) rz(:, i)
    end do
    close (unit)
    modes = [(i, i=0, nmax)]
    allocate (bn(3, nmax + 1, np))
    call coil_set_modes(set, w, rz, nphi, modes, bn, st)
    if (st /= 0) error stop 'projection failed'
    call get_command_argument(4, arg)
    open (newunit=unit, file=trim(arg), status='replace', action='write')
    write (unit, '(a)') '# n R Z ReBR ImBR ReBphi ImBphi ReBZ ImBZ  (two-sided B_n, T)'
    do i = 1, np
        do im = 1, nmax + 1
            write (unit, '(i4,2f14.8,6es24.15e3)') modes(im), rz(:, i), &
                bn(1, im, i), bn(2, im, i), bn(3, im, i)
        end do
    end do
    close (unit)
end program
