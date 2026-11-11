!> Checks output/hasegawa_wakatani.nc of the saturated benchmark case
!> alpha = kappa = 1 on the box 2 pi/0.15 at 128^2 against the band recorded
!> in benchmarks/drift-turbulence/README.md.
program test_hasegawa_wakatani_benchmark
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortio, only: fortio_file_t, fortio_status_t
    implicit none
    !> Flux of the 256^2 and 512^2 runs and the half-width that covers the
    !> spread between resolutions, seeds and dissipation settings.
    real(dp), parameter :: flux_reference = 0.584_dp, flux_band = 0.10_dp
    type(fortio_file_t) :: file
    type(fortio_status_t) :: status
    real(dp) :: flux, flux_error, sink, hyper, zonal, residual
    integer :: failures

    failures = 0
    call file%open('output/hasegawa_wakatani.nc', status); call require('open')
    call file%read('Gamma_n', flux, status); call require('Gamma_n')
    call file%read('Gamma_n_error', flux_error, status); call require('Gamma_n_error')
    call file%read('D_alpha', sink, status); call require('D_alpha')
    call file%read('D_E', hyper, status); call require('D_E')
    call file%read('zonal_fraction', zonal, status); call require('zonal_fraction')
    call file%read('energy_balance_residual', residual, status); call require('residual')
    call file%close(status)
    write (*, '(a,2es12.4)') 'Gamma_n and its standard error: ', flux, flux_error
    call check('flux inside the benchmark band', abs(flux - flux_reference), flux_band)
    call check('standard error below 5 % of the flux', flux_error/flux, 0.05_dp)
    call check('mean energy balance', abs(residual)/flux, 0.01_dp)
    ! Numata, Ball and Dewar (2007), Sec. III: Gamma_n ~ D_alpha >> D_E.
    call check('flux balanced by the resistive sink', abs(sink/flux - 1), 0.1_dp)
    call check('hyperviscous sink small', hyper/flux, 0.1_dp)
    ! The original model keeps no dominant zonal flow.
    call check('zonal share of the kinetic energy', zonal, 0.2_dp)
    if (failures > 0) error stop 'test_hasegawa_wakatani_benchmark failed'
    write (*, '(a)') 'test_hasegawa_wakatani_benchmark passed'
contains
    subroutine require(what)
        character(len=*), intent(in) :: what
        if (.not. status%ok()) then
            print '(a)', what//': '//status%message
            error stop 'FortIO read failed'
        end if
    end subroutine

    subroutine check(name, value, tolerance)
        character(*), intent(in) :: name
        real(dp), intent(in) :: value, tolerance
        if (value <= tolerance) then
            write (*, '(a,a,es10.3)') 'ok   ', name//': ', value
        else
            write (*, '(a,a,es10.3,a,es10.3)') 'FAIL ', name//': ', value, ' > ', tolerance
            failures = failures + 1
        end if
    end subroutine
end program test_hasegawa_wakatani_benchmark
