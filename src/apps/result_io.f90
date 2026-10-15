!> Configuration and result contract for KIN6D applications (campaign F5).
!>
!> Input: an optional Fortran namelist file named by `--config FILE`. Each
!> application owns its namelist groups and their defaults; this module only
!> opens the file, checks the read status, captures the literal text and the
!> resolved (defaults materialized) text, and rejects groups no reader owns.
!>
!> Output: `output/<application>.nc`, a NetCDF-4/HDF5 file written by FortIO
!> with named real variables (scalar, 1-3D), units/long_name and an optional
!> per-variable error/certificate status, plus global provenance attributes.
module kin6d_result_io
    use, intrinsic :: iso_fortran_env, only: dp => real64, compiler_version, &
        compiler_options
    use netcdf, only: nf90_create, nf90_close, nf90_def_dim, nf90_def_var, &
        nf90_put_att, nf90_put_var, nf90_enddef, nf90_redef, nf90_strerror, &
        NF90_NOERR, NF90_NETCDF4, NF90_CLOBBER, NF90_DOUBLE, NF90_GLOBAL
    implicit none
    private
    public :: runtime_git
    include 'kin6d_build_info.inc'

    !> Directory that receives every application result.
    character(len=*), parameter, public :: result_directory = 'output'
    character(len=*), parameter, public :: result_schema = 'kin6d-result-v1'
    character(len=*), parameter :: no_config = &
        '(no configuration file; all defaults)'

    type :: text_t
        character(len=:), allocatable :: value
    end type

    !> One application run: command line, literal and resolved configuration.
    type, public :: run_t
        character(len=:), allocatable :: name, config_path, literal, resolved
        character(len=:), allocatable :: command_line
        type(text_t), allocatable :: arguments(:)
    contains
        procedure :: start => run_start
        procedure :: argument_count => run_argument_count
        procedure :: argument => run_argument
        procedure :: has_config => run_has_config
        procedure :: open_config => run_open_config
        procedure :: close_config => run_close_config
        procedure :: open_resolved => run_open_resolved
        procedure :: close_resolved => run_close_resolved
    end type

    !> NetCDF-4 result file with named dimensions shared by name.
    type, public :: result_file_t
        integer :: ncid = -1
        character(len=:), allocatable :: path
        integer :: dimension_count = 0
        character(len=64) :: dimension_names(64) = ''
        integer :: dimension_lengths(64) = 0, dimension_ids(64) = -1
    contains
        procedure :: create => result_create
        procedure :: close => result_close
        procedure :: attribute => result_attribute
        procedure, private :: write_scalar, write_rank1, write_rank2, write_rank3
        generic :: write => write_scalar, write_rank1, write_rank2, write_rank3
    end type

contains

    !> Parse `--config FILE` (or `--config=FILE`) and keep the remaining
    !> arguments for the application. Reads the literal configuration text.
    subroutine run_start(this, name)
        class(run_t), intent(inout) :: this
        character(len=*), intent(in) :: name
        character(len=:), allocatable :: arg
        integer :: i, n, kept
        this%name = name
        this%config_path = ''
        this%literal = no_config
        this%resolved = ''
        this%command_line = command_text()
        n = command_argument_count()
        allocate (this%arguments(n))
        kept = 0
        i = 1
        do while (i <= n)
            arg = argument_text(i)
            if (arg == '--config' .or. index(arg, '--config=') == 1) then
                if (len(this%config_path) > 0) &
                    call fail(name, '--config given more than once')
                if (arg == '--config') then
                    if (i == n) call fail(name, '--config needs a file name')
                    i = i + 1
                    this%config_path = argument_text(i)
                else
                    this%config_path = arg(len('--config=') + 1:)
                end if
                if (len(this%config_path) == 0) call fail(name, 'empty --config path')
            else
                kept = kept + 1
                this%arguments(kept)%value = arg
            end if
            i = i + 1
        end do
        this%arguments = this%arguments(:kept)
        if (len(this%config_path) > 0) &
            this%literal = file_text(this%config_path, name)
    end subroutine

    integer function run_argument_count(this) result(n)
        class(run_t), intent(in) :: this
        n = size(this%arguments)
    end function

    !> Application argument i, excluding the `--config` pair.
    subroutine run_argument(this, i, value)
        class(run_t), intent(in) :: this
        integer, intent(in) :: i
        character(len=*), intent(out) :: value
        if (i < 1 .or. i > size(this%arguments)) call fail(this%name, 'argument index')
        value = this%arguments(i)%value
    end subroutine

    !> True when `--config FILE` was given (otherwise all defaults).
    logical function run_has_config(this) result(has)
        class(run_t), intent(in) :: this
        has = len(this%config_path) > 0
    end function

    !> Open the configuration file for one namelist read when has_config();
    !> each call starts at the file head. Read only if has_config().
    subroutine run_open_config(this, unit)
        class(run_t), intent(in) :: this
        integer, intent(out) :: unit
        integer :: ios
        unit = 0
        if (.not. this%has_config()) return
        open (newunit=unit, file=this%config_path, status='old', action='read', &
            iostat=ios)
        if (ios /= 0) call fail(this%name, 'cannot open '//this%config_path)
    end subroutine

    !> Close after `read (unit, nml=group, iostat=ios, iomsg=msg)`. A missing
    !> group keeps its defaults; any other read error (unknown entry, bad
    !> value) is fatal and names the group.
    subroutine run_close_config(this, unit, group, ios, msg)
        class(run_t), intent(in) :: this
        integer, intent(in) :: unit
        character(len=*), intent(in) :: group, msg
        integer, intent(in) :: ios
        if (.not. this%has_config()) return
        close (unit)
        if (ios == 0 .or. is_iostat_end(ios)) return
        call fail(this%name, 'invalid namelist &'//group//' in '// &
            this%config_path//': '//trim(msg))
    end subroutine

    !> Scratch unit that receives `write (unit, nml=group)` for every group
    !> the application owns, after defaults and overrides are resolved.
    subroutine run_open_resolved(this, unit)
        class(run_t), intent(in) :: this
        integer, intent(out) :: unit
        integer :: ios
        open (newunit=unit, status='scratch', action='readwrite', form='formatted', &
            iostat=ios)
        if (ios /= 0) call fail(this%name, 'cannot open resolved scratch file')
    end subroutine

    !> Capture the resolved text and reject literal groups no reader owns.
    subroutine run_close_resolved(this, unit)
        class(run_t), intent(inout) :: this
        integer, intent(in) :: unit
        character(len=4096) :: line
        character(len=:), allocatable :: owned
        integer :: ios, size_read
        rewind (unit)
        this%resolved = ''
        do
            read (unit, '(a)', advance='no', iostat=ios, size=size_read) line
            if (is_iostat_end(ios)) exit
            if (ios /= 0 .and. .not. is_iostat_eor(ios)) &
                call fail(this%name, 'cannot read resolved configuration')
            this%resolved = this%resolved//line(:size_read)
            if (is_iostat_eor(ios)) this%resolved = this%resolved//new_line('a')
        end do
        close (unit)
        owned = group_names(this%resolved)
        call check_groups(this%name, group_names(this%literal), owned)
    end subroutine

    !> Create `output/<run name>.nc` and write the global provenance.
    subroutine result_create(this, run)
        class(result_file_t), intent(inout) :: this
        type(run_t), intent(in) :: run
        character(len=:), allocatable :: commit, dirty
        if (len(run%resolved) == 0) call fail(run%name, &
            'resolved configuration must be captured before the result file')
        this%path = result_directory//'/'//run%name//'.nc'
        this%dimension_count = 0
        call check(nf90_create(this%path, ior(NF90_NETCDF4, NF90_CLOBBER), &
            this%ncid), this%path)
        call runtime_git(commit, dirty, run%name)
        call this%attribute('kin6d_schema', result_schema)
        call this%attribute('kin6d_application', run%name)
        call this%attribute('kin6d_command_line', run%command_line)
        call this%attribute('kin6d_config_path', nonempty(run%config_path))
        call this%attribute('kin6d_config_literal', nonempty(run%literal))
        call this%attribute('kin6d_config_resolved', run%resolved)
        call this%attribute('kin6d_git_commit', kin6d_git_commit)
        call this%attribute('kin6d_git_tree', kin6d_git_dirty)
        call this%attribute('kin6d_git_runtime_commit', commit)
        call this%attribute('kin6d_git_runtime_tree', dirty)
        call this%attribute('kin6d_source_directory', kin6d_source_path)
        call this%attribute('kin6d_build_type', nonempty(kin6d_build_type))
        call this%attribute('kin6d_build_flags', nonempty(kin6d_build_flags))
        call this%attribute('kin6d_build_compiler', nonempty(kin6d_build_compiler))
        call this%attribute('compiler_version', compiler_version())
        call this%attribute('compiler_options', compiler_options())
        call this%attribute('fortsym_revision', kin6d_fortsym_revision)
        call this%attribute('fortnum_revision', kin6d_fortnum_revision)
        call this%attribute('fortplot_revision', kin6d_fortplot_revision)
        call this%attribute('fortio_revision', kin6d_fortio_revision)
        call this%attribute('fortfem_revision', kin6d_fortfem_revision)
        call this%attribute('created_utc', utc_timestamp())
    end subroutine

    subroutine result_close(this)
        class(result_file_t), intent(inout) :: this
        if (this%ncid < 0) return
        call check(nf90_close(this%ncid), this%path)
        this%ncid = -1
    end subroutine

    !> Global text attribute (e.g. an application-level certificate summary).
    subroutine result_attribute(this, name, value)
        class(result_file_t), intent(inout) :: this
        character(len=*), intent(in) :: name, value
        call check(nf90_put_att(this%ncid, NF90_GLOBAL, name, value), this%path)
    end subroutine

    subroutine write_scalar(this, name, value, units, long_name, status)
        class(result_file_t), intent(inout) :: this
        character(len=*), intent(in) :: name, units, long_name
        real(dp), intent(in) :: value
        character(len=*), intent(in), optional :: status
        integer :: varid
        call define(this, name, [integer ::], '', units, long_name, status, varid)
        call check(nf90_put_var(this%ncid, varid, value), this%path)
        call check(nf90_redef(this%ncid), this%path)
    end subroutine

    !> `dims` is a comma-separated list of dimension names (default
    !> `<name>_dim<k>`); equal names must have equal lengths.
    subroutine write_rank1(this, name, value, units, long_name, dims, status)
        class(result_file_t), intent(inout) :: this
        character(len=*), intent(in) :: name, units, long_name
        real(dp), intent(in) :: value(:)
        character(len=*), intent(in), optional :: dims, status
        integer :: varid
        call define(this, name, shape(value), optional_text(dims), units, &
            long_name, status, varid)
        call check(nf90_put_var(this%ncid, varid, value), this%path)
        call check(nf90_redef(this%ncid), this%path)
    end subroutine

    subroutine write_rank2(this, name, value, units, long_name, dims, status)
        class(result_file_t), intent(inout) :: this
        character(len=*), intent(in) :: name, units, long_name
        real(dp), intent(in) :: value(:, :)
        character(len=*), intent(in), optional :: dims, status
        integer :: varid
        call define(this, name, shape(value), optional_text(dims), units, &
            long_name, status, varid)
        call check(nf90_put_var(this%ncid, varid, value), this%path)
        call check(nf90_redef(this%ncid), this%path)
    end subroutine

    subroutine write_rank3(this, name, value, units, long_name, dims, status)
        class(result_file_t), intent(inout) :: this
        character(len=*), intent(in) :: name, units, long_name
        real(dp), intent(in) :: value(:, :, :)
        character(len=*), intent(in), optional :: dims, status
        integer :: varid
        call define(this, name, shape(value), optional_text(dims), units, &
            long_name, status, varid)
        call check(nf90_put_var(this%ncid, varid, value), this%path)
        call check(nf90_redef(this%ncid), this%path)
    end subroutine

    !> Define a double variable and its attributes, then leave define mode.
    subroutine define(this, name, extents, dims, units, long_name, status, varid)
        type(result_file_t), intent(inout) :: this
        character(len=*), intent(in) :: name, dims, units, long_name
        integer, intent(in) :: extents(:)
        character(len=*), intent(in), optional :: status
        integer, intent(out) :: varid
        integer :: ids(size(extents)), k, first, last
        character(len=64) :: dim_name
        first = 1
        do k = 1, size(extents)
            if (len(dims) == 0) then
                write (dim_name, '(a,"_dim",i0)') name, k
            else
                last = index(dims(first:), ',')
                if (last == 0) then
                    last = len(dims) + 1
                else
                    last = first + last - 1
                end if
                if (last <= first) call fail(this%path, 'empty dimension name for '//name)
                dim_name = adjustl(dims(first:last - 1))
                first = last + 1
            end if
            ids(k) = dimension_id(this, trim(dim_name), extents(k))
        end do
        if (len(dims) > 0 .and. first <= len(dims)) &
            call fail(this%path, 'too many dimension names for '//name)
        call check(nf90_def_var(this%ncid, name, NF90_DOUBLE, ids, varid), this%path)
        call check(nf90_put_att(this%ncid, varid, 'units', nonempty(units)), this%path)
        call check(nf90_put_att(this%ncid, varid, 'long_name', long_name), this%path)
        if (present(status)) call check(nf90_put_att(this%ncid, varid, &
            'kin6d_status', status), this%path)
        call check(nf90_enddef(this%ncid), this%path)
    end subroutine

    integer function dimension_id(this, name, extent) result(id)
        type(result_file_t), intent(inout) :: this
        character(len=*), intent(in) :: name
        integer, intent(in) :: extent
        integer :: k
        do k = 1, this%dimension_count
            if (this%dimension_names(k) /= name) cycle
            if (this%dimension_lengths(k) /= extent) call fail(this%path, &
                'dimension '//name//' reused with a different length')
            id = this%dimension_ids(k)
            return
        end do
        if (this%dimension_count == size(this%dimension_names)) &
            call fail(this%path, 'too many dimensions')
        if (extent < 1) call fail(this%path, 'empty dimension '//name)
        call check(nf90_def_dim(this%ncid, name, extent, id), this%path)
        this%dimension_count = this%dimension_count + 1
        this%dimension_names(this%dimension_count) = name
        this%dimension_lengths(this%dimension_count) = extent
        this%dimension_ids(this%dimension_count) = id
    end function

    subroutine check(code, context)
        integer, intent(in) :: code
        character(len=*), intent(in) :: context
        if (code /= NF90_NOERR) call fail(context, trim(nf90_strerror(code)))
    end subroutine

    subroutine fail(context, message)
        character(len=*), intent(in) :: context, message
        error stop 'kin6d '//context//': '//message
    end subroutine

    function nonempty(text) result(value)
        character(len=*), intent(in) :: text
        character(len=:), allocatable :: value
        value = text
        if (len_trim(text) == 0) value = '(none)'
    end function

    function optional_text(text) result(value)
        character(len=*), intent(in), optional :: text
        character(len=:), allocatable :: value
        value = ''
        if (present(text)) value = trim(text)
    end function

    function command_text() result(text)
        character(len=:), allocatable :: text
        integer :: length
        call get_command(length=length)
        allocate (character(len=max(length, 0)) :: text)
        if (length > 0) call get_command(text)
    end function

    function argument_text(i) result(text)
        integer, intent(in) :: i
        character(len=:), allocatable :: text
        integer :: length
        call get_command_argument(i, length=length)
        allocate (character(len=length) :: text)
        if (length > 0) call get_command_argument(i, text)
    end function

    !> Exact bytes of a text file.
    function file_text(path, context) result(text)
        character(len=*), intent(in) :: path, context
        character(len=:), allocatable :: text
        integer :: unit, ios, bytes
        open (newunit=unit, file=path, status='old', action='read', &
            access='stream', form='unformatted', iostat=ios)
        if (ios /= 0) call fail(context, 'cannot open configuration '//path)
        inquire (unit=unit, size=bytes)
        if (bytes < 0) call fail(context, 'cannot size configuration '//path)
        allocate (character(len=bytes) :: text)
        if (bytes > 0) read (unit, iostat=ios) text
        close (unit)
        if (ios /= 0) call fail(context, 'cannot read configuration '//path)
    end function

    !> Lower-case namelist group names `&name` that start a line, as
    !> ",a,b,"; comments after '!' and the `&end` terminator are ignored.
    function group_names(text) result(names)
        character(len=*), intent(in) :: text
        character(len=:), allocatable :: names, line, name
        integer :: first, last, k
        names = ','
        first = 1
        do while (first <= len(text))
            last = index(text(first:), new_line('a'))
            if (last == 0) then
                last = len(text) + 1
            else
                last = first + last - 1
            end if
            line = adjustl(text(first:last - 1))
            first = last + 1
            if (len_trim(line) < 2) cycle
            if (line(1:1) /= '&') cycle
            k = scan(line(2:), ' /!,'//achar(9)//achar(13))
            if (k == 0) then
                name = lower(trim(line(2:)))
            else
                name = lower(line(2:k))
            end if
            if (name == 'end' .or. len(name) == 0) cycle
            if (index(names, ','//name//',') == 0) names = names//name//','
        end do
    end function

    subroutine check_groups(context, literal, owned)
        character(len=*), intent(in) :: context, literal, owned
        integer :: first, last
        first = 2
        do while (first < len(literal))
            last = first + index(literal(first:), ',') - 1
            if (index(owned, ','//literal(first:last - 1)//',') == 0) &
                call fail(context, 'unknown namelist group &'// &
                literal(first:last - 1)//'; this application reads '//owned)
            first = last + 1
        end do
    end subroutine

    pure function lower(text) result(value)
        character(len=*), intent(in) :: text
        character(len=len(text)) :: value
        integer :: k
        value = text
        do k = 1, len(text)
            if (text(k:k) >= 'A' .and. text(k:k) <= 'Z') &
                value(k:k) = achar(iachar(text(k:k)) + 32)
        end do
    end function

    !> Runtime recheck of the source revision (one short git call); records
    !> "unavailable" when git or the build-time source tree is absent.
    subroutine runtime_git(commit, tree, name)
        character(len=:), allocatable, intent(out) :: commit, tree
        character(len=*), intent(in) :: name
        character(len=:), allocatable :: scratch, source
        character(len=128) :: line
        integer :: unit, ios, exit_status, changed
        commit = 'unavailable'; tree = 'unavailable'
        scratch = result_directory//'/.'//name//'.git-check'
        source = "'"//kin6d_source_path//"'"
        call execute_command_line('git -C '//source//' rev-parse HEAD > '//scratch// &
            ' 2>/dev/null && git -C '//source//' status --porcelain '// &
            '--untracked-files=no 2>/dev/null | wc -l >> '//scratch, &
            exitstat=exit_status, cmdstat=ios)
        if (ios /= 0) return
        open (newunit=unit, file=scratch, status='old', action='read', iostat=ios)
        if (ios /= 0) return
        if (exit_status == 0) then
            read (unit, '(a)', iostat=ios) line
            if (ios == 0) commit = trim(line)
            read (unit, *, iostat=ios) changed
            if (ios == 0) then
                tree = 'clean'
                if (changed > 0) tree = 'dirty'
            end if
        end if
        close (unit, status='delete')
    end subroutine

    !> ISO 8601 UTC time from the local clock and its zone offset.
    function utc_timestamp() result(text)
        character(len=:), allocatable :: text
        character(len=20) :: buffer
        integer :: v(8), y, m, d, era, yoe, doy, doe, mp
        integer(8) :: days, minutes
        call date_and_time(values=v)
        y = v(1); m = v(2); d = v(3)
        if (m <= 2) y = y - 1
        era = floor(real(y)/400.0)
        yoe = y - era*400
        mp = mod(m + 9, 12)
        doy = (153*mp + 2)/5 + d - 1
        doe = yoe*365 + yoe/4 - yoe/100 + doy
        days = int(era, 8)*146097 + doe - 719468
        minutes = days*1440 + v(5)*60 + v(6) - v(4)
        days = floor(real(minutes, dp)/1440.0_dp, 8)
        minutes = minutes - days*1440
        ! Inverse civil-from-days (Hinnant).
        days = days + 719468
        era = int(floor(real(days, dp)/146097.0_dp))
        doe = int(days - int(era, 8)*146097)
        yoe = (doe - doe/1460 + doe/36524 - doe/146096)/365
        y = yoe + era*400
        doy = doe - (365*yoe + yoe/4 - yoe/100)
        mp = (5*doy + 2)/153
        d = doy - (153*mp + 2)/5 + 1
        m = mp + 3
        if (m > 12) m = m - 12
        if (m <= 2) y = y + 1
        write (buffer, '(i4.4,"-",i2.2,"-",i2.2,"T",i2.2,":",i2.2,":",i2.2,"Z")') &
            y, m, d, int(minutes/60), int(mod(minutes, 60_8)), v(7)
        text = buffer
    end function
end module kin6d_result_io
