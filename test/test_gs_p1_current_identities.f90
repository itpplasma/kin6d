program test_gs_p1_current_identities
    use fortsym_arena, only: arena_t
    use fortsym_expr, only: expr_t, sym, operator(-)
    use fortsym_parse, only: parse_expr
    use fortsym_diff, only: diff
    use fortsym_engine_native, only: native_engine_t, make_native_engine
    use fortsym_check, only: suite_t, suite_begin, check_identity, suite_end
    implicit none
    type(arena_t), target :: arena
    type(native_engine_t) :: eng
    type(suite_t) :: suite
    type(expr_t) :: r, t
    integer :: i
    character(16) :: label
    character(100), parameter :: affine_i(5) = [character(100) :: &
        '1/r**2', '1/r', '1/(2*r**2)', '1/(2*r)', &
        '(-1+3/r-2/r**2)/6']
    character(100), parameter :: affine_a(5) = [character(100) :: &
        '-1/r', 'log(r)', '-1/(2*r)', 'log(r)/2', &
        '(-r+3*log(r)+2/r)/6']
    character(100), parameter :: split_i(5) = [character(100) :: &
        '(2-r)/r**2', '(2-r)/r', '1/r-1/2', '1-r/2', &
        '(-2*r**3+15*r**2-43*r+58-36/r+8/r**2)/6']
    character(100), parameter :: split_a(5) = [character(100) :: &
        '-2/r-log(r)', '2*log(r)-r', 'log(r)-r/2', 'r-r**2/4', &
        '(-r**4/2+5*r**3-43*r**2/2+58*r-36*log(r)-8/r)/6']
    character(100), parameter :: edge_i(5) = [character(100) :: &
        '2/(1+t)', '2', '2*t/(1+t)', '2*t', '2*t**2*(1-t)**2/(1+t)']
    character(100), parameter :: edge_a(5) = [character(100) :: &
        '2*log(1+t)', '2*t', '2*(t-log(1+t))', 't**2', &
        '2*(t**4/4-t**3+2*t**2-4*t+4*log(1+t))']
    call arena%init()
    eng = make_native_engine(arena)
    r = sym(arena, 'r'); t = sym(arena, 't')
    call suite_begin(suite, 'P1 weak current exact antiderivatives; r>0, 0<=t<=1')
    do i = 1, 5
        write (label, '(a,i1)') 'affine cell ', i
        call check_identity(suite, eng, trim(label), &
            diff(parsed(affine_a(i)), r) - parsed(affine_i(i)))
        write (label, '(a,i1)') 'split cell ', i
        call check_identity(suite, eng, trim(label), &
            diff(parsed(split_a(i)), r) - parsed(split_i(i)))
        write (label, '(a,i1)') 'split edge ', i
        call check_identity(suite, eng, trim(label), &
            diff(parsed(edge_a(i)), t) - parsed(edge_i(i)))
    end do
    call suite_end(suite, 'gs-p1-current-antiderivatives.json')
contains
    function parsed(text) result(e)
        character(*), intent(in) :: text
        type(expr_t) :: e
        logical :: ok
        character(:), allocatable :: message
        e = parse_expr(arena, trim(text), ok, message)
        if (.not. ok) then
            print *, message
            error stop 'P1 oracle expression parse failed'
        end if
    end function parsed
end program test_gs_p1_current_identities
