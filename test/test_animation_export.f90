!> Animation-export contract: binary32 round trip, 8/16-bit quantization
!> decoded by the documented rule (independent decode here) within half a
!> quantization step, below-range code 0 and saturation, and JSON string
!> escaping of quotes, backslashes and newlines in provenance.json.
program test_animation_export
    use, intrinsic :: iso_fortran_env, only: dp => real64, sp => real32, &
        int8, int16
    use kin6d_result_io, only: run_t
    use kin6d_animation_export, only: animation_export_t
    implicit none
    character(len=*), parameter :: root = 'output/test_animation_export'
    type(run_t) :: run
    type(animation_export_t) :: export
    real(dp), parameter :: values(6) = [-3.0_dp, -1.0_dp, -0.123456_dp, 0.5_dp, &
                                        2.0_dp, 7.0_dp]
    real(dp), parameter :: lo = -1.0_dp, hi = 2.0_dp
    real(sp) :: back(6)
    integer(int8) :: bytes(6)
    integer(int16) :: words(6)
    character(len=4096) :: json
    integer :: unit, i, unit_size

    call run%start('test_animation_export')
    call run%open_resolved(unit)
    write (unit, '(a)') '&probe text=''say "hi"\path'' /'
    call run%close_resolved(unit)
    call execute_command_line('rm -rf '//root//' && mkdir -p output')
    call export%create(run, root, 'probe')
    call export%note('quote', 'a "quoted" \ word'//new_line('a')//'next line')
    call export%parameter('scale', 1.25_dp, 'lambda_D', 'probe scale')
    call export%real32('values', values, [2, 3], '1', 'probe values')
    call export%quantized('coarse', values, [6], lo, hi, 8, '1', 'probe 8 bit')
    call export%quantized('fine', values, [6], lo, hi, 16, '1', 'probe 16 bit')
    call export%close()

    open (newunit=unit, file=root//'/probe/values.f32', access='stream', &
        form='unformatted', status='old', action='read')
    read (unit) back
    close (unit)
    if (any(back /= real(values, sp))) error stop 'binary32 round trip failed'

    open (newunit=unit, file=root//'/probe/coarse.u8', access='stream', &
        form='unformatted', status='old', action='read')
    inquire (unit=unit, size=unit_size)
    if (unit_size /= 6) error stop '8-bit file size'
    read (unit) bytes
    close (unit)
    call check_codes([(iand(int(bytes(i)), 255), i = 1, 6)], 255)

    open (newunit=unit, file=root//'/probe/fine.u16', access='stream', &
        form='unformatted', status='old', action='read')
    inquire (unit=unit, size=unit_size)
    if (unit_size /= 12) error stop '16-bit file size'
    read (unit) words
    close (unit)
    call check_codes([(iand(int(words(i)), 65535), i = 1, 6)], 65535)

    json = ''
    open (newunit=unit, file=root//'/probe/provenance.json', access='stream', &
        form='unformatted', status='old', action='read')
    inquire (unit=unit, size=unit_size)
    if (unit_size > len(json)) error stop 'provenance.json unexpectedly large'
    read (unit) json(:unit_size)
    close (unit)
    if (index(json, '"a \"quoted\" \\ word\nnext line"') == 0) &
        error stop 'note not escaped'
    if (index(json, 'say \"hi\"\\path') == 0) error stop 'resolved text not escaped'
    if (index(json, '"shape": [2, 3]') == 0) error stop 'shape missing'
    if (index(json, '"schema": "kin6d-animation-export-v1"') == 0) &
        error stop 'schema missing'
    print '(a)', 'animation export contract: ok'
contains
    !> Independent decode: v = lo + (code - 1) (hi - lo)/(cmax - 1).
    subroutine check_codes(codes, cmax)
        integer, intent(in) :: codes(:), cmax
        real(dp) :: step, decoded
        integer :: k
        step = (hi - lo)/(cmax - 1)
        do k = 1, size(codes)
            if (values(k) < lo) then
                if (codes(k) /= 0) error stop 'below-range value not coded 0'
            else if (values(k) > hi) then
                if (codes(k) /= cmax) error stop 'above-range value not saturated'
            else
                decoded = lo + (codes(k) - 1)*step
                if (abs(decoded - values(k)) > step/2*(1 + 1.0e-12_dp)) &
                    error stop 'quantized value outside half a step'
            end if
        end do
    end subroutine
end program test_animation_export
