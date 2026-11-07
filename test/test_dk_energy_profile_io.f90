program test_dk_energy_profile_io
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
    use fortio, only: fortio_file_t, fortio_status_t, hdf5_attribute_t
    use kin6d_dk_energy, only: energy_axi_t, ENERGY_OK, ENERGY_BAD_INPUT
    use kin6d_dk_force_balance, only: energy_balance_input_t, energy_balance_t, &
        energy_force_balance_cgs
    use kin6d_dk_energy_profile_io
    use kin6d_result_io, only: run_t
    implicit none
    type(energy_balance_input_t) :: p(2), input
    type(energy_axi_t) :: response(2)
    type(energy_balance_t) :: expected(2)
    type(run_t) :: run
    type(fortio_file_t) :: file
    type(fortio_status_t) :: read_status
    type(hdf5_attribute_t), allocatable :: attributes(:)
    real(dp), allocatable :: back(:), err(:, :)
    real(dp) :: radius, value
    integer :: unit, j, status, solver_status(2), rows = 2
    logical :: exists
    namelist /profile_probe/ rows

    call run%start('test_dk_energy_profile_io')
    call run%open_resolved(unit)
    write (unit, nml=profile_probe)
    call run%close_resolved(unit)
    call execute_command_line('mkdir -p output')
    do j = 1, 2
        p(j) = energy_balance_input_t(k=-999.0_dp, n=1e13_dp, t=1e-9_dp, &
            nr=-1e11_dp, tr=-j*1e-12_dp, b=2e4_dp, bz=-1.9e4_dp, r=real(j, dp), &
            vz=1e5_dp, q=2.0_dp, r0=165.0_dp, charge=4.8032047e-10_dp)
        response(j)%k = 0.5_dp*j
        response(j)%d31 = 0.5_dp
        response(j)%d32 = 0.5_dp*(2.5_dp - response(j)%k)
        response(j)%dk_ds = 0.01_dp*j
        response(j)%residual = [1e-14_dp, 2e-14_dp]
        input = p(j); input%k = response(j)%k
        call energy_force_balance_cgs(input, expected(j), status)
        if (status /= ENERGY_OK) error stop 'force-balance fixture'
    end do
    solver_status = ENERGY_OK
    call energy_profile_write_cgs(run, p, response, solver_status, status)
    if (status /= PROFILE_OK) error stop 'profile export failed'
    call check_legacy('Er.dat', expected%er)
    call check_legacy('Vth.dat', expected%vth)
    call file%open('output/test_dk_energy_profile_io.nc', read_status)
    call require('open')
    call file%read('r_eff', back, read_status); call require('r_eff')
    if (any(back /= p%r)) error stop 'radius round trip'
    call file%read('Er', back, read_status); call require('Er')
    if (any(back /= expected%er)) error stop 'electric field round trip'
    call file%read('Vth', back, read_status); call require('Vth')
    if (any(back /= expected%vth)) error stop 'velocity round trip'
    call file%read('k', back, read_status); call require('k')
    if (any(back /= response%k)) error stop 'stale force-balance k used'
    call file%read('solve_backward_error', err, read_status); call require('backward error')
    if (any(shape(err) /= [2, 2])) error stop 'backward-error shape'
    if (any(err(:, 1) /= response(1)%residual)) error stop 'backward-error values'
    call file%get_attributes('Er', attributes, read_status); call require('Er attributes')
    call expect('units', 'statV/cm'); call expect('kin6d_status', 'unavailable')
    call file%get_attributes('Vth', attributes, read_status); call require('Vth attributes')
    call expect('units', 'cm/s')
    call file%get_attributes('/', attributes, read_status); call require('provenance')
    call expect('kin6d_config_literal', run%literal)
    call expect('kin6d_config_resolved', run%resolved)
    call expect('kamel_export', 'complete')
    call file%close(read_status); call require('close')

    solver_status(1) = ENERGY_BAD_INPUT
    call energy_profile_write_cgs(run, p, response, solver_status, status)
    if (status /= PROFILE_INCOMPLETE) error stop 'failed solver row accepted'
    inquire (file='output/Er.dat', exist=exists)
    if (exists) error stop 'stale Er profile survived failed row'
    inquire (file='output/Vth.dat', exist=exists)
    if (exists) error stop 'stale Vth profile survived failed row'
    call file%open('output/test_dk_energy_profile_io.nc', read_status); call require('failed profile')
    call file%read('Er', back, read_status); call require('failed row')
    if (.not. ieee_is_nan(back(1)) .or. back(2) /= expected(2)%er) &
        error stop 'failed row discarded valid row'
    call file%read('solver_status', back, read_status); call require('solver statuses')
    if (any(back /= real(solver_status, dp))) error stop 'solver statuses lost'
    call file%get_attributes('/', attributes, read_status); call require('failed provenance')
    call expect('kamel_export', 'incomplete')
    call file%close(read_status); call require('failed close')

    solver_status = ENERGY_OK
    p(1)%n = -1
    call energy_profile_write_cgs(run, p, response, solver_status, status)
    if (status /= PROFILE_INCOMPLETE) error stop 'failed force balance accepted'
    call file%open('output/test_dk_energy_profile_io.nc', read_status); call require('balance profile')
    call file%read('force_balance_status', back, read_status); call require('balance statuses')
    if (any(back /= real([ENERGY_BAD_INPUT, ENERGY_OK], dp))) error stop 'balance statuses lost'
    call file%close(read_status); call require('balance close')
    call energy_profile_write_cgs(run, p, response(:1), solver_status, status)
    if (status /= PROFILE_BAD_INPUT) error stop 'shape mismatch accepted'
    p(2)%r = p(1)%r
    call energy_profile_write_cgs(run, p, response, solver_status, status)
    if (status /= PROFILE_BAD_INPUT) error stop 'duplicate radius accepted'
    print '(a)', 'energy profile export PASS'
contains
    subroutine check_legacy(name, expected_values)
        character(*), intent(in) :: name
        real(dp), intent(in) :: expected_values(:)
        integer :: ios
        open (newunit=unit, file='output/'//name, status='old', action='read')
        do j = 1, size(expected_values)
            read (unit, *, iostat=ios) radius, value
            if (ios /= 0) error stop 'legacy row missing'
            if (radius /= p(j)%r .or. value /= expected_values(j)) error stop 'legacy value or units'
        end do
        read (unit, *, iostat=ios) radius, value
        if (ios >= 0) error stop 'extra legacy row'
        close (unit)
    end subroutine check_legacy

    subroutine require(what)
        character(*), intent(in) :: what
        if (.not. read_status%ok()) then
            print '(a)', what//': '//read_status%message
            error stop 'FortIO read failed'
        end if
    end subroutine require

    subroutine expect(name, expected_value)
        character(*), intent(in) :: name, expected_value
        integer :: i
        do i = 1, size(attributes)
            if (attributes(i)%name /= name) cycle
            if (.not. allocated(attributes(i)%value_text)) error stop 'not text: '//name
            if (attributes(i)%value_text /= expected_value) error stop 'attribute mismatch: '//name
            return
        end do
        error stop 'missing attribute: '//name
    end subroutine expect
end program test_dk_energy_profile_io
