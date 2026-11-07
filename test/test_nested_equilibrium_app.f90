program test_nested_equilibrium_app
    ! Checks output/nested_equilibrium.nc against the same input solved by
    ! PARVMEC: usage RESULT RAXIS WB [VOLUME]. The axis must agree to 1.5 mm,
    ! the magnetic energy to 2e-5 and the volume (when given) to 1e-4.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use netcdf, only: nf90_open, nf90_close, nf90_inq_varid, nf90_get_var, &
        nf90_get_att, NF90_NOWRITE, NF90_NOERR, NF90_GLOBAL
    implicit none
    character(256) :: path, arg
    character(8) :: converged
    real(dp) :: volume, raxis, wb, expected(3)
    integer :: ncid, k
    call get_command_argument(1, path)
    expected = -1
    do k = 2, min(command_argument_count(), 4)
        call get_command_argument(k, arg); read (arg, *) expected(k - 1)
    end do
    if (nf90_open(trim(path), NF90_NOWRITE, ncid) /= NF90_NOERR) error stop 'open result'
    if (nf90_get_att(ncid, NF90_GLOBAL, 'converged', converged) /= NF90_NOERR) &
        error stop 'converged attribute'
    volume = scalar('volume_p'); raxis = scalar('raxis_cc'); wb = scalar('wb')
    if (nf90_close(ncid) /= NF90_NOERR) error stop 'close result'
    print '(a,3es16.8)', 'volume, raxis, wb: ', volume, raxis, wb
    if (trim(converged) /= 'true') error stop 'not converged'
    if (abs(raxis - expected(1)) > 1.5e-3_dp) error stop 'axis'
    if (abs(wb/expected(2) - 1) > 2e-5_dp) error stop 'magnetic energy'
    if (expected(3) > 0 .and. abs(volume/expected(3) - 1) > 1e-4_dp) error stop 'volume'
    print '(a)', 'PASS nested equilibrium'
contains
    real(dp) function scalar(name)
        character(*), intent(in) :: name
        integer :: varid
        if (nf90_inq_varid(ncid, name, varid) /= NF90_NOERR) error stop name
        if (nf90_get_var(ncid, varid, scalar) /= NF90_NOERR) error stop name
    end function
end program
