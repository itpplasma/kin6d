!> D31 of mono_coefficients_3d on a surface grid of surface_grid.py over a list
!! of nu* = pi R nu, and its ratio to a reference value (README.md).
program qh_scan
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use milca_dk_mono, only: mono_coefficients_t, mono_coefficients_3d
    implicit none
    real(dp), parameter :: pi = 4*atan(1.0_dp)
    type(mono_coefficients_t) :: c
    real(dp), allocatable :: b(:, :)
    real(dp) :: s, iota, g, ci, b0, rmaj, d31lc, nustar
    integer :: nt, nz, nfp, nl, u, i, j, st, k, nn
    character(256) :: fname, arg
    call get_command_argument(1, fname)
    call get_command_argument(2, arg); read (arg, *) nl
    call get_command_argument(3, arg); read (arg, *) rmaj
    call get_command_argument(4, arg); read (arg, *) d31lc
    open (newunit=u, file=fname, status='old')
    read (u, *) nt, nz, nfp
    read (u, *) s, iota, g, ci, b0
    allocate (b(nt, nz))
    read (u, *) ((b(i, j), i=1, nt), j=1, nz)
    close (u)
    nn = command_argument_count() - 4
    do k = 1, nn
        call get_command_argument(4 + k, arg); read (arg, *) nustar
        call mono_coefficients_3d(b, nfp, 1/iota, ci/g, iota*b0/(g + iota*ci), &
            nustar/(pi*rmaj), nl, c, st)
        print '(a,es9.2,a,i3,a,i4,a,es20.12,a,f12.8,a,es9.2)', 'nu*', nustar, ' st', st, &
            ' L', nl, ' D31', c%d31, ' ratio', c%d31/d31lc, ' D13+D31', c%d13 + c%d31
    end do
end program
