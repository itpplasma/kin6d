module kin6d_gs_timers
    ! Per-phase wall clocks for the Grad--Shafranov path. Phases are accumulated
    ! by the master thread only, outside parallel regions. A module-level
    ! accumulator keeps the intent(in) numerical interfaces unchanged.
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64, output_unit
    implicit none
    private
    public :: gs_wall_s, gs_timer_reset, gs_timer_add, gs_timer_total
    public :: gs_timer_calls, gs_timer_iterations, gs_timer_event_count
    public :: gs_timer_event_get, gs_progress_enabled, gs_progress_set
    public :: gs_progress, gs_phase_tick
    integer, parameter, public :: GS_T_ASSEMBLY = 1   ! gs_cubic_load
    integer, parameter, public :: GS_T_FIXED_CG = 2   ! solve_gs linear CG
    integer, parameter, public :: GS_T_JAC_CG = 3     ! Jacobian symmetric CG
    integer, parameter, public :: GS_T_JAC_GMRES = 4  ! Jacobian GMRES (incl. nested)
    integer, parameter, public :: GS_T_PRECOND_CG = 5 ! nested CG inside GMRES
    integer, parameter, public :: GS_T_FACTOR = 6    ! Cholesky factorisation of A
    integer, parameter, public :: GS_T_DIRECT = 7    ! back-substitutions with A
    integer, parameter, public :: GS_T_LOAD_CELLS = 8 ! element loop of gs_cubic_load
    integer, parameter, public :: GS_T_COUNT = 8
    integer, parameter :: max_events = 512
    real(dp) :: seconds(GS_T_COUNT) = 0
    integer(int64) :: calls(GS_T_COUNT) = 0, iterations(GS_T_COUNT) = 0
    integer :: events = 0
    character(16) :: event_label(max_events) = ''
    real(dp) :: event_seconds(max_events) = 0
    integer :: event_iterations(max_events) = 0
    logical :: progress_on = .false.
contains
    subroutine gs_phase_tick(label, t0)
        !! Wall time since the previous tick, written to stderr as
        !! `PROF label seconds` when KIN6D_PROFILE is set. t0 < 0 starts a
        !! sequence. Serial sections only.
        character(*), intent(in) :: label
        real(dp), intent(inout) :: t0
        logical, save :: known = .false., enabled = .false.
        character(len=8) :: value
        integer :: length, stat
        real(dp) :: t

        if (.not. known) then
            call get_environment_variable('KIN6D_PROFILE', value, length, stat)
            enabled = stat == 0 .and. length > 0
            known = .true.
        end if
        if (.not. enabled) return
        t = gs_wall_s()
        if (t0 > 0) write (0, '(a,a,f10.3)') 'PROF ', label, t - t0
        t0 = t
    end subroutine

    function gs_wall_s() result(t)
        real(dp) :: t
        integer(int64) :: count, rate
        call system_clock(count, rate)
        t = real(count, dp)/real(rate, dp)
    end function

    subroutine gs_timer_reset()
        seconds = 0; calls = 0; iterations = 0; events = 0
    end subroutine

    subroutine gs_timer_add(kind, elapsed, count)
        ! Accumulate one timed call. Top-level solves and assemblies are also
        ! logged one by one; nested preconditioner solves are aggregated only.
        integer, intent(in) :: kind
        real(dp), intent(in) :: elapsed
        integer, intent(in), optional :: count
        integer :: n
        n = 0
        if (present(count)) n = count
        seconds(kind) = seconds(kind) + elapsed
        calls(kind) = calls(kind) + 1
        iterations(kind) = iterations(kind) + n
        if (kind == GS_T_PRECOND_CG .or. kind == GS_T_DIRECT) return
        if (events >= max_events) return
        events = events + 1
        event_label(events) = kind_name(kind)
        event_seconds(events) = elapsed
        event_iterations(events) = n
    end subroutine

    function kind_name(kind) result(name)
        integer, intent(in) :: kind
        character(16) :: name
        select case (kind)
        case (GS_T_ASSEMBLY); name = 'assembly'
        case (GS_T_FIXED_CG); name = 'fixed_cg'
        case (GS_T_JAC_CG); name = 'jacobian_cg'
        case (GS_T_JAC_GMRES); name = 'jacobian_gmres'
        case (GS_T_FACTOR); name = 'factor'
        case default; name = 'other'
        end select
    end function

    real(dp) function gs_timer_total(kind) result(t)
        integer, intent(in) :: kind
        t = seconds(kind)
    end function

    integer(int64) function gs_timer_calls(kind) result(n)
        integer, intent(in) :: kind
        n = calls(kind)
    end function

    integer(int64) function gs_timer_iterations(kind) result(n)
        integer, intent(in) :: kind
        n = iterations(kind)
    end function

    integer function gs_timer_event_count() result(n)
        n = events
    end function

    subroutine gs_timer_event_get(i, label, elapsed, count)
        integer, intent(in) :: i
        character(*), intent(out) :: label
        real(dp), intent(out) :: elapsed
        integer, intent(out) :: count
        label = event_label(i); elapsed = event_seconds(i); count = event_iterations(i)
    end subroutine

    subroutine gs_progress_set(on)
        logical, intent(in) :: on
        progress_on = on
    end subroutine

    logical function gs_progress_enabled() result(on)
        on = progress_on
    end function

    subroutine gs_progress(text)
        character(*), intent(in) :: text
        if (.not. progress_on) return
        write (output_unit, '(a)') trim(text)
        flush (output_unit)
    end subroutine
end module
