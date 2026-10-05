!> Round trip for the F5 configuration/result contract. The probe owns a
!> namelist group, writes output/test_result_io.nc and output/
!> test_result_io.resolved.nml, then reads the file back through FortIO's
!> independent HDF5 reader (not the netcdf writer adapter) and checks
!> values, dimensions and provenance attributes. test/result_io_check.cmake
!> adds ncdump/h5dump and the rejection cases.
program test_result_io
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortio, only: fortio_file_t, fortio_status_t, hdf5_attribute_t
    use kin6d_result_io, only: run_t, result_file_t
    implicit none
    real(dp) :: scale = 1.5_dp
    integer :: count = 3
    real(dp) :: weights(3) = [1.0_dp, 2.0_dp, 4.0_dp]
    character(len=32) :: label = 'probe'
    namelist /probe/ scale, count, weights, label
    type(run_t) :: run
    type(result_file_t) :: result
    type(fortio_file_t) :: file
    type(fortio_status_t) :: status
    type(hdf5_attribute_t), allocatable :: attributes(:)
    real(dp), allocatable :: back1(:), back2(:, :), back3(:, :, :)
    real(dp) :: back0, cube(2, 3, 4), grid(3, 2)
    integer :: unit, ios, i, j, k
    character(len=256) :: msg
    character(len=64) :: argument

    call run%start('test_result_io')
    if (run%argument_count() > 1) error stop 'usage: test_result_io [--config F] [tag]'
    argument = 'none'
    if (run%argument_count() == 1) call run%argument(1, argument)
    ios = 0; msg = ''
    call run%open_config(unit)
    if (run%has_config()) read (unit, nml=probe, iostat=ios, iomsg=msg)
    call run%close_config(unit, 'probe', ios, msg)
    if (count < 1) error stop 'count must be positive'
    call run%open_resolved(unit)
    write (unit, nml=probe)
    call run%close_resolved(unit)
    call execute_command_line('mkdir -p output')
    open (newunit=unit, file='output/test_result_io.resolved.nml', status='replace', &
        access='stream', form='unformatted')
    write (unit) run%resolved
    close (unit)

    do k = 1, 4
        do j = 1, 3
            do i = 1, 2
                cube(i, j, k) = scale*(i + 10*j + 100*k)
            end do
        end do
    end do
    grid = reshape([(real(i, dp), i=1, 6)], [3, 2])
    call result%create(run)
    call result%attribute('probe_argument', trim(argument))
    call result%write('scale', scale, 'm', 'probe scale', status='certified-discrete')
    call result%write('weights', weights*count, 'm', 'weighted probe', dims='w')
    call result%write('grid', grid, '1', 'column-major grid', dims='w,pair')
    call result%write('cube', cube, 's', 'rank-3 probe')
    call result%close()

    call file%open('output/test_result_io.nc', status)
    call require(status, 'open')
    call file%read('scale', back0, status); call require(status, 'scale')
    call file%read('weights', back1, status); call require(status, 'weights')
    call file%read('grid', back2, status); call require(status, 'grid')
    call file%read('cube', back3, status); call require(status, 'cube')
    if (back0 /= scale) error stop 'scalar round trip'
    if (any(back1 /= weights*count)) error stop 'rank-1 round trip'
    if (any(shape(back2) /= [3, 2])) error stop 'rank-2 shape'
    if (any(back2 /= grid)) error stop 'rank-2 round trip'
    if (any(shape(back3) /= [2, 3, 4])) error stop 'rank-3 shape'
    if (any(back3 /= cube)) error stop 'rank-3 round trip'
    call file%get_attributes('scale', attributes, status); call require(status, 'attrs')
    call expect(attributes, 'units', 'm')
    call expect(attributes, 'kin6d_status', 'certified-discrete')
    call file%get_attributes('/', attributes, status); call require(status, 'globals')
    call expect(attributes, 'kin6d_application', 'test_result_io')
    call expect(attributes, 'kin6d_config_resolved', run%resolved)
    call expect(attributes, 'kin6d_config_literal', run%literal)
    call expect(attributes, 'probe_argument', trim(argument))
    call file%close(status)
    call require(status, 'close')
    print '(a)', 'result_io round trip PASS'

contains

    subroutine require(status, what)
        type(fortio_status_t), intent(in) :: status
        character(len=*), intent(in) :: what
        if (.not. status%ok()) then
            print '(a)', what//': '//status%message
            error stop 'FortIO read failed'
        end if
    end subroutine

    subroutine expect(attributes, name, value)
        type(hdf5_attribute_t), intent(in) :: attributes(:)
        character(len=*), intent(in) :: name, value
        integer :: n
        do n = 1, size(attributes)
            if (attributes(n)%name /= name) cycle
            if (.not. allocated(attributes(n)%value_text)) error stop 'not text: '//name
            if (attributes(n)%value_text /= value) then
                print '(a)', name//' = ['//attributes(n)%value_text//']'
                error stop 'attribute mismatch'
            end if
            return
        end do
        error stop 'missing attribute '//name
    end subroutine
end program test_result_io
