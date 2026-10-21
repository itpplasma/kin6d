program check_response_strip_current
    !! Independent cylinder oracle for R-B2 physical-current readback.
    !! Positive radii and mu0; integer harmonic m>=1; no numeric zero probes.
    use fortsym_arena, only: arena_t
    use fortsym_expr, only: expr_t, sym, operator(+), operator(-), &
        operator(*), operator(/), operator(**)
    use fortsym_diff, only: diff
    use fortsym_subs, only: subs
    use fortsym_engine, only: engine_result_t, VERDICT_TRUE
    use fortsym_engine_native, only: native_engine_t, make_native_engine
    implicit none
    type(arena_t), target :: arena
    type(native_engine_t) :: engine
    type(expr_t) :: r, rs, major, mu, amplitude, inner, outer, jump
    type(expr_t) :: a, b, h, determinant, li, lo, symmetric
    integer :: m, count
    call arena%init()
    engine = make_native_engine(arena)
    r = sym(arena, 'r'); rs = sym(arena, 'rs')
    major = sym(arena, 'major'); mu = sym(arena, 'mu')
    amplitude = sym(arena, 'amplitude')
    a = sym(arena, 'a'); b = sym(arena, 'b')
    h = sym(arena, 'h')
    count = 0
    do m = 1, 6
        ! Magnetic scalar potential has jump Lambda; B=mu0 grad U and
        ! Phi=2 R0 r B_r (the production peak-flux convention).
        inner = -mu*m*major*amplitude*(r/rs)**m
        outer = -mu*m*major*amplitude*(rs/r)**m
        call zero(diff(r*diff(inner, r), r)/r - m*m*inner/(r*r))
        call zero(diff(r*diff(outer, r), r)/r - m*m*outer/(r*r))
        call zero(subs(outer - inner, r, rs))
        jump = subs(diff(outer - inner, r), r, rs)
        call zero(rs*jump - 2*mu*m*m*major*amplitude)
        ! Artificial edge sheets reproduce the inner/outer r^(+-m) moments
        ! of a sheet on rs. Their unweighted sum need not equal its amplitude.
        determinant = a**m*b**(-m) - a**(-m)*b**m
        li = (b**(-m) - b**m)/determinant
        lo = (a**m - a**(-m))/determinant
        call zero(li*a**m + lo*b**m - 1)
        call zero(li*a**(-m) + lo*b**(-m) - 1)
        if (m == 2) then
            call zero(li + lo - (a*a*b*b + 1)/(a*a + b*b))
            symmetric = subs(subs(li + lo - 1, a, 1 - h), b, 1 + h)
            call zero(symmetric - h*h*(h*h - 4)/(2*(1 + h*h)))
        end if
    end do
    print '(a,i0,a)', 'PASS: ', count, ' exact FortSym cylinder-current identities (no probes)'
contains
    subroutine zero(expression)
        type(expr_t), intent(in) :: expression
        type(engine_result_t) :: result
        result = engine%zero_test(expression)
        if (result%verdict /= VERDICT_TRUE) then
            print *, 'Unproved cylinder-current identity', count + 1, 'm=', m
            error stop 1
        end if
        count = count + 1
    end subroutine zero
end program check_response_strip_current
