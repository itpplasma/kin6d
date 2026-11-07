program test_nested_equilibrium_app
    ! Checks output/nested_equilibrium.nc of examples/nested-equilibrium/
    ! circular-tokamak.nml: exact volume 2 pi^2 R a^2, and the axis and magnetic
    ! energy of the same input solved by PARVMEC (ns 101, mpol 8: raxis
    ! 3.03634 m, wb 0.0750949) and GVEC (raxis 3.03683 m).
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use netcdf, only: nf90_open, nf90_close, nf90_inq_varid, nf90_get_var, &
        nf90_get_att, NF90_NOWRITE, NF90_NOERR, NF90_GLOBAL
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp)
    character(256) :: path
    character(8) :: converged
    real(dp) :: volume, raxis, wb
    integer :: ncid
    call get_command_argument(1, path)
    if (nf90_open(trim(path), NF90_NOWRITE, ncid) /= NF90_NOERR) error stop 'open result'
    if (nf90_get_att(ncid, NF90_GLOBAL, 'converged', converged) /= NF90_NOERR) &
        error stop 'converged attribute'
    volume = scalar('volume_p'); raxis = scalar('raxis_cc'); wb = scalar('wb')
    if (nf90_close(ncid) /= NF90_NOERR) error stop 'close result'
    print '(a,3es16.8)', 'volume, raxis, wb: ', volume, raxis, wb
    if (trim(converged) /= 'true') error stop 'not converged'
    if (abs(volume/(2*pi**2*3) - 1) > 1e-4_dp) error stop 'volume'
    if (abs(raxis - 3.0366_dp) > 1.5e-3_dp) error stop 'axis'
    if (abs(wb/0.0750949_dp - 1) > 2e-5_dp) error stop 'magnetic energy'
    print '(a)', 'PASS nested equilibrium circular tokamak'
contains
    real(dp) function scalar(name)
        character(*), intent(in) :: name
        integer :: varid
        if (nf90_inq_varid(ncid, name, varid) /= NF90_NOERR) error stop name
        if (nf90_get_var(ncid, varid, scalar) /= NF90_NOERR) error stop name
    end function
end program
