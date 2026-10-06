!> Compact animation-data export for external visualizations (lecture Manim
!> scenes). A visualization consumes the exported arrays unchanged; the
!> scientific source remains the kin6d application that wrote them.
!>
!> Layout of an export directory DIR/<slug>/:
!> - one raw little-endian file per array, `<name>.f32` (IEEE binary32) or
!>   `<name>.u8`/`<name>.u16` (unsigned codes of a linear quantization), stored in
!>   Fortran (column-major) order with the shape recorded in the header;
!> - `provenance.json`: application, command line, literal and resolved
!>   namelist, kin6d commit/tree, compiler, dependency revisions, declared
!>   physical parameters with units, and one entry per array (file, dtype,
!>   shape, order, units, description, decode rule for quantized arrays).
!> No timestamp is written, so equal inputs give byte-identical exports.
module kin6d_animation_export
    use, intrinsic :: iso_fortran_env, only: dp => real64, sp => real32, &
        int8, int16, int32, compiler_version
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use kin6d_result_io, only: run_t, runtime_git
    implicit none
    private
    include 'kin6d_build_info.inc'

    character(len=*), parameter, public :: animation_export_schema = &
        'kin6d-animation-export-v1'

    !> One export directory; entries accumulate until close().
    type, public :: animation_export_t
        character(len=:), allocatable :: path, parameters, arrays, notes
        type(run_t) :: run
    contains
        procedure :: create => export_create
        procedure :: parameter => export_parameter
        procedure :: note => export_note
        procedure :: real32 => export_real32
        procedure :: quantized => export_quantized
        procedure :: close => export_close
    end type

contains

    !> Create DIR/slug (parents included) for the export of `run`.
    subroutine export_create(this, run, directory, slug)
        class(animation_export_t), intent(inout) :: this
        type(run_t), intent(in) :: run
        character(len=*), intent(in) :: directory, slug
        integer :: status
        if (transfer(1_int32, 0_int8) /= 1_int8) &
            error stop 'animation export requires a little-endian host'
        this%path = trim(directory)//'/'//trim(slug)
        call execute_command_line("mkdir -p '"//this%path//"'", exitstat=status)
        if (status /= 0) error stop 'cannot create animation export directory'
        this%parameters = ''; this%arrays = ''; this%notes = ''
        this%run = run
    end subroutine

    !> Declared physical/numerical parameter with units.
    subroutine export_parameter(this, name, value, units, description)
        class(animation_export_t), intent(inout) :: this
        character(len=*), intent(in) :: name, units, description
        real(dp), intent(in) :: value
        call append(this%parameters, quoted(name)//': {"value": '//number(value)// &
            ', "units": '//quoted(units)//', "description": '//quoted(description)//'}')
    end subroutine

    !> Free-text statement (model, status, claim limits).
    subroutine export_note(this, name, text)
        class(animation_export_t), intent(inout) :: this
        character(len=*), intent(in) :: name, text
        call append(this%notes, quoted(name)//': '//quoted(text))
    end subroutine

    !> Array of any rank, passed flattened in Fortran order with its shape,
    !> written as binary32.
    subroutine export_real32(this, name, values, shape, units, description)
        class(animation_export_t), intent(inout) :: this
        character(len=*), intent(in) :: name, units, description
        real(dp), intent(in) :: values(:)
        integer, intent(in) :: shape(:)
        integer :: unit
        if (product(shape) /= size(values)) error stop 'export shape mismatch: '//name
        open (newunit=unit, file=this%path//'/'//name//'.f32', access='stream', &
            form='unformatted', status='replace', action='write')
        write (unit) real(values, sp)
        close (unit)
        call append(this%arrays, quoted(name)//': {"file": '//quoted(name//'.f32')// &
            ', "dtype": "float32", "endianness": "little", "order": "F", "shape": '// &
            int_list(shape)//', "units": '//quoted(units)//', "description": '// &
            quoted(description)//'}')
    end subroutine

    !> Linear quantization of values in [lo, hi] to `bits` = 8 or 16 bit
    !> unsigned codes: code 0 marks values below lo (or not finite), codes
    !> 1..cmax (cmax = 2**bits - 1) decode to lo + (code - 1) (hi - lo)/(cmax - 1);
    !> values above hi saturate at cmax.
    subroutine export_quantized(this, name, values, shape, lo, hi, bits, units, &
                                description)
        class(animation_export_t), intent(inout) :: this
        character(len=*), intent(in) :: name, units, description
        real(dp), intent(in) :: values(:), lo, hi
        integer, intent(in) :: shape(:), bits
        integer(int16), allocatable :: words(:)
        integer(int8), allocatable :: bytes(:)
        integer :: unit, i, code, cmax
        character(len=:), allocatable :: suffix, dtype
        character(len=16) :: buffer
        if (product(shape) /= size(values)) error stop 'export shape mismatch: '//name
        if (.not. hi > lo) error stop 'export quantization needs hi > lo: '//name
        if (bits /= 8 .and. bits /= 16) error stop 'export quantization: 8 or 16 bits'
        cmax = 2**bits - 1
        suffix = merge('.u8 ', '.u16', bits == 8); suffix = trim(suffix)
        dtype = merge('uint8 ', 'uint16', bits == 8); dtype = trim(dtype)
        allocate (words(size(values)))
        do i = 1, size(values)
            if (.not. values(i) >= lo) then
                code = 0
            else
                code = 1 + nint(min(1.0_dp, (values(i) - lo)/(hi - lo))*(cmax - 1))
            end if
            ! Two's-complement storage of the unsigned code.
            if (code > cmax/2) code = code - cmax - 1
            words(i) = int(code, int16)
        end do
        open (newunit=unit, file=this%path//'/'//name//suffix, access='stream', &
            form='unformatted', status='replace', action='write')
        if (bits == 8) then
            bytes = int(words, int8)
            write (unit) bytes
        else
            write (unit) words
        end if
        close (unit)
        write (buffer, '(i0)') cmax - 1
        call append(this%arrays, quoted(name)//': {"file": '//quoted(name//suffix)// &
            ', "dtype": '//quoted(dtype)//', "endianness": "little", "order": "F", '// &
            '"shape": '//int_list(shape)//', "units": '//quoted(units)// &
            ', "description": '//quoted(description)// &
            ', "decode": {"rule": "code 0: below lo or not finite; code >= 1: '// &
            'lo + (code - 1) * (hi - lo) / '//trim(buffer)//'", "lo": '//number(lo)// &
            ', "hi": '//number(hi)//'}}')
    end subroutine

    !> Write provenance.json.
    subroutine export_close(this)
        class(animation_export_t), intent(inout) :: this
        character(len=:), allocatable :: commit, tree
        integer :: unit
        call runtime_git(commit, tree, this%run%name)
        open (newunit=unit, file=this%path//'/provenance.json', status='replace', &
            action='write')
        write (unit, '(a)') '{'
        call field('schema', quoted(animation_export_schema))
        call field('application', quoted(this%run%name))
        call field('command_line', quoted(this%run%command_line))
        call field('config_path', quoted(this%run%config_path))
        call field('config_literal', quoted(this%run%literal))
        call field('config_resolved', quoted(this%run%resolved))
        call field('kin6d_git_commit_build', quoted(kin6d_git_commit))
        call field('kin6d_git_tree_build', quoted(kin6d_git_dirty))
        call field('kin6d_git_commit_runtime', quoted(commit))
        call field('kin6d_git_tree_runtime', quoted(tree))
        call field('kin6d_build_type', quoted(kin6d_build_type))
        call field('compiler_version', quoted(compiler_version()))
        call field('fortnum_revision', quoted(kin6d_fortnum_revision))
        call field('fortsym_revision', quoted(kin6d_fortsym_revision))
        call field('notes', '{'//this%notes//'}')
        call field('parameters', '{'//this%parameters//'}')
        write (unit, '(a)') '  "arrays": {'//this%arrays//'}'
        write (unit, '(a)') '}'
        close (unit)
    contains
        subroutine field(key, value)
            character(len=*), intent(in) :: key, value
            write (unit, '(a)') '  '//quoted(key)//': '//value//','
        end subroutine
    end subroutine

    subroutine append(list, entry)
        character(len=:), allocatable, intent(inout) :: list
        character(len=*), intent(in) :: entry
        if (len(list) > 0) list = list//','
        list = list//new_line('a')//'    '//entry
    end subroutine

    !> JSON string literal with escapes for quote, backslash and controls.
    function quoted(text) result(json)
        character(len=*), intent(in) :: text
        character(len=:), allocatable :: json
        character(len=6) :: hex
        integer :: i, c
        json = '"'
        do i = 1, len(text)
            c = iachar(text(i:i))
            select case (c)
            case (34); json = json//'\"'
            case (92); json = json//'\\'
            case (10); json = json//'\n'
            case (9); json = json//'\t'
            case (13); json = json//'\r'
            case (0:8, 11:12, 14:31)
                write (hex, '(a,z4.4)') '\u', c
                json = json//hex
            case default; json = json//text(i:i)
            end select
        end do
        json = json//'"'
    end function

    !> Shortest round-trip-safe decimal (17 significant digits).
    function number(value) result(text)
        real(dp), intent(in) :: value
        character(len=:), allocatable :: text
        character(len=32) :: buffer
        if (.not. ieee_is_finite(value)) &
            error stop 'animation export parameter is not finite'
        write (buffer, '(es24.16e3)') value
        text = trim(adjustl(buffer))
    end function

    function int_list(values) result(text)
        integer, intent(in) :: values(:)
        character(len=:), allocatable :: text
        character(len=16) :: buffer
        integer :: i
        text = '['
        do i = 1, size(values)
            write (buffer, '(i0)') values(i)
            text = text//trim(buffer)
            if (i < size(values)) text = text//', '
        end do
        text = text//']'
    end function
end module kin6d_animation_export
