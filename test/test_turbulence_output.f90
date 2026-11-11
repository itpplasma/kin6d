!> Reads output/hasegawa_wakatani.nc of a small saturated run through
!> FortIO's HDF5 reader and checks the output contract: axes, units,
!> Parseval identity of the spectrum, the particle flux recovered from the
!> cross spectrum, and the standard errors. A cross spectrum of the opposite
!> sign convention is rejected.
program test_turbulence_output
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortio, only: fortio_file_t, fortio_status_t, hdf5_attribute_t
    implicit none
    type(fortio_file_t) :: file
    type(fortio_status_t) :: status
    type(hdf5_attribute_t), allocatable :: attributes(:)
    real(dp), allocatable :: kx(:), ky(:), omega(:), s(:, :, :), cross(:, :, :), c_t(:)
    real(dp) :: variance, flux, flux_error, sink, hyper, value, cell, spectral_flux, tau, ell
    integer :: failures, j

    failures = 0
    call file%open('output/hasegawa_wakatani.nc', status); call require('open')
    call file%read('kx', kx, status); call require('kx')
    call file%read('ky', ky, status); call require('ky')
    call file%read('omega', omega, status); call require('omega')
    call file%read('S_phi_phi', s, status); call require('S_phi_phi')
    call file%read('S_n_phi_im', cross, status); call require('S_n_phi_im')
    call file%read('C_phi_t', c_t, status); call require('C_phi_t')
    call file%read('phi_variance', variance, status); call require('phi_variance')
    call file%read('Gamma_n', flux, status); call require('Gamma_n')
    call file%read('Gamma_n_error', flux_error, status); call require('Gamma_n_error')
    call file%read('D_alpha', sink, status); call require('D_alpha')

    if (any(shape(s) /= [size(kx), size(ky), size(omega)])) failures = failures + 1
    if (size(kx) /= 21 .or. size(ky) /= 21 .or. size(omega) /= 256) failures = failures + 1
    call check('kx axis symmetric', abs(kx(1) + kx(21)) + abs(kx(11)) + abs(kx(12) - 0.3_dp), &
               1e-14_dp)
    call check('ky axis symmetric', abs(ky(1) + ky(21)) + abs(ky(11)), 1e-14_dp)
    call check('omega spacing', abs(omega(2) - omega(1) - 2*acos(-1.0_dp)/(256*0.4_dp)), 1e-12_dp)
    cell = 0.3_dp**2*(omega(2) - omega(1))
    call check('Parseval: spectrum sum equals variance', abs(sum(s)*cell/variance - 1), 1e-12_dp)
    call check('autocorrelation at zero lag equals variance', abs(c_t(1)/variance - 1), 1e-12_dp)
    if (.not. (minval(s) >= 0 .and. variance > 0.1_dp)) failures = failures + 1

    ! Gamma_n = kappa <n v_x>, v_x = -d phi/dy: -kappa sum ky Im <n conjg(phi)>.
    spectral_flux = 0
    do j = 1, size(ky)
        spectral_flux = spectral_flux - ky(j)*sum(cross(:, j, :))*cell
    end do
    write (*, '(a,3es12.4)') 'Gamma_n: series, error, from cross spectrum ', flux, flux_error, &
        spectral_flux
    if (.not. (flux > 0.05_dp .and. flux_error > 0 .and. flux_error < 0.2_dp*flux)) &
        failures = failures + 1
    call check('flux from the cross spectrum', abs(spectral_flux/flux - 1), 0.1_dp)
    ! Mutant: the opposite sign convention of the cross spectrum.
    call check_fails('cross spectrum of opposite sign is rejected', &
                     abs(-spectral_flux/flux - 1), 0.1_dp)
    call file%read('D_E', hyper, status); call require('D_E')
    if (.not. (sink > 0.5_dp*flux .and. hyper > 0)) failures = failures + 1
    call check('saturated state: flux balances the sinks', &
               abs(flux - sink - hyper)/flux, 0.02_dp)

    call file%read('tau_c_euler', tau, status); call require('tau_c_euler')
    call file%read('ell_x', ell, status); call require('ell_x')
    if (.not. (tau > 0 .and. ell > 0)) failures = failures + 1
    call file%read('tau_c_euler_error', value, status); call require('tau_c_euler_error')
    if (.not. (value > 0 .and. value < tau)) failures = failures + 1
    call file%read('unit_length', value, status); call require('unit_length')
    call check('unit_length', abs(value - 0.2_dp), 0.0_dp)
    call file%read('unit_potential', value, status); call require('unit_potential')
    call check('unit_potential', abs(value - 2.5_dp), 0.0_dp)
    call file%get_attributes('S_phi_phi', attributes, status); call require('attributes')
    call expect('units', '(T_e/e)^2 (rho_s)^2 1/Omega_i')
    call file%get_attributes('unit_time', attributes, status); call require('attributes')
    call expect('units', 's')
    call file%get_attributes('Gamma_n_error', attributes, status); call require('attributes')
    call expect('milca_status', 'candidate; statistical standard error')
    call file%get_attributes('/', attributes, status); call require('globals')
    call expect('turbulence_model', 'Hasegawa-Wakatani')
    call expect('turbulence_provider', 'nonlinear')
    call expect('derivative_status', 'unavailable: time averages of a chaotic state')
    call file%close(status)
    if (failures > 0) error stop 'test_turbulence_output failed'
    write (*, '(a)') 'test_turbulence_output passed'
contains
    subroutine require(what)
        character(len=*), intent(in) :: what
        if (.not. status%ok()) then
            print '(a)', what//': '//status%message
            error stop 'FortIO read failed'
        end if
    end subroutine

    subroutine expect(name, value)
        character(len=*), intent(in) :: name, value
        integer :: n
        do n = 1, size(attributes)
            if (attributes(n)%name /= name) cycle
            if (.not. allocated(attributes(n)%value_text)) error stop 'not text: '//name
            if (trim(attributes(n)%value_text) /= value) then
                print '(a)', name//' = ['//attributes(n)%value_text//']'
                failures = failures + 1
            end if
            return
        end do
        print '(a)', 'missing attribute '//name
        failures = failures + 1
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

    subroutine check_fails(name, value, tolerance)
        character(*), intent(in) :: name
        real(dp), intent(in) :: value, tolerance
        if (value > tolerance) then
            write (*, '(a,a,es10.3)') 'ok   ', name//': ', value
        else
            write (*, '(a,a,es10.3)') 'FAIL ', name//': ', value
            failures = failures + 1
        end if
    end subroutine
end program test_turbulence_output
