program test_gs_coupled_majorant_identities
    !! Exact FortSym check of the coupled residual identity behind the P10
    !! majorant (DESIGN section 4): for psi solving -div(grad psi/R) = f in
    !! Omega with the natural condition d_nu psi + S psi = g on its boundary,
    !! every psi_h, test function v and flux y satisfy
    !!
    !!   a(e, v) + <S e, v> = (y - grad psi_h/R, grad v) + (f + div y, v)
    !!                        + <g - S psi_h - y.n, v>,     e = psi - psi_h,
    !!
    !! d_nu = (1/R) d_n, n the outward normal. Checked with exact definite
    !! integrals on the box [1,2] x [0,1] for a polynomial psi and S the
    !! multiplication by sigma (a symmetric S enters only through g). The
    !! residual is affine in (psi_h, y) and linear in v, so checking it for
    !! every basis monomial of v with psi_h = y = 0, with each basis monomial
    !! of psi_h, and with each basis monomial of either component of y covers
    !! the general polynomial case of those degrees. Also Ampere,
    !! -<d_nu psi, 1> = (f, 1), which is the Ip goal.
    !! Negative controls: dropping the Neumann defect or flipping the normal
    !! must leave a nonzero remainder.
    use fortsym_arena, only: arena_t
    use fortsym_expr, only: expr_t, sym, num, operator(+), operator(-), &
        operator(*), operator(/)
    use fortsym_parse, only: parse_expr
    use fortsym_diff, only: diff
    use fortsym_subs, only: subs
    use fortsym_defint, only: definite_integral
    use fortsym_engine_native, only: native_engine_t, make_native_engine
    use fortsym_engine, only: engine_result_t, VERDICT_FALSE
    use fortsym_check, only: suite_t, suite_begin, check_identity, suite_end
    implicit none
    type(arena_t), target :: arena
    type(native_engine_t) :: engine
    type(suite_t) :: suite
    type(expr_t) :: r, z, psi, f, sigma, zero, ampere
    character(len=*), parameter :: v_basis(5) = [character(8) :: '1', 'R', 'Z', &
        'R*Z', 'Z**2']
    character(len=*), parameter :: h_basis(4) = [character(8) :: 'R**2*Z', 'R*Z', &
        'Z**2', '1']
    character(len=*), parameter :: y_basis(4) = [character(8) :: '1', 'Z', 'R', &
        'R*Z']
    integer :: i, j
    call arena%init()
    engine = make_native_engine(arena)
    r = sym(arena, 'R'); z = sym(arena, 'Z')
    zero = num(arena, 0)
    psi = parsed('R**2*Z + Z**3 - R*Z**2 + 2*R**3')
    sigma = parsed('3')
    f = zero - diff(diff(psi, r)/r, r) - diff(diff(psi, z)/r, z)
    call suite_begin(suite, 'Coupled FEM-BEM majorant residual identity')
    do i = 1, size(v_basis)
        call check_identity(suite, engine, 'identity, v = '//trim(v_basis(i)), &
            residual(zero, zero, zero, parsed(v_basis(i)), .false.))
        do j = 1, size(h_basis)
            call check_identity(suite, engine, 'identity, v = '//trim(v_basis(i)) &
                //', psi_h = '//trim(h_basis(j)), &
                residual(parsed(h_basis(j)), zero, zero, parsed(v_basis(i)), .false.))
        end do
        do j = 1, size(y_basis)
            call check_identity(suite, engine, 'identity, v = '//trim(v_basis(i)) &
                //', y_R = '//trim(y_basis(j)), &
                residual(zero, parsed(y_basis(j)), zero, parsed(v_basis(i)), .false.))
            call check_identity(suite, engine, 'identity, v = '//trim(v_basis(i)) &
                //', y_Z = '//trim(y_basis(j)), &
                residual(zero, zero, parsed(y_basis(j)), parsed(v_basis(i)), .false.))
        end do
    end do
    ampere = zero - side_sum(diff(psi, r)/r, zero - diff(psi, r)/r, &
                             diff(psi, z)/r, zero - diff(psi, z)/r)
    call check_identity(suite, engine, 'Ampere: -<d_nu psi, 1> = (f, 1)', &
                        ampere - box(f))
    call require_nonzero('omitting the Neumann defect', &
        residual(zero, parsed('Z'), zero, parsed('R'), .false.) &
        + side_sum(neumann(zero, parsed('Z'), zero, parsed('R'), 2, .false.), &
                   neumann(zero, parsed('Z'), zero, parsed('R'), 1, .false.), &
                   neumann(zero, parsed('Z'), zero, parsed('R'), 3, .false.), &
                   neumann(zero, parsed('Z'), zero, parsed('R'), 4, .false.)))
    call require_nonzero('flipping the boundary normal', &
        residual(zero, parsed('Z'), parsed('R'), parsed('1+R'), .true.))
    call suite_end(suite, 'gs-coupled-majorant-identities.json')
contains
    function residual(psih, yr, yz, v, flip) result(value)
        !! a(e,v) + <S e,v> - (y - grad psi_h/R, grad v) - (f + div y, v)
        !! - <g - S psi_h - y.n, v>; flip uses -n in the last term.
        type(expr_t), intent(in) :: psih, yr, yz, v
        logical, intent(in) :: flip
        type(expr_t) :: value, e
        e = psi - psih
        value = box((diff(e, r)*diff(v, r) + diff(e, z)*diff(v, z))/r) &
            + side_sum(sigma*e*v, sigma*e*v, sigma*e*v, sigma*e*v) &
            - box((yr - diff(psih, r)/r)*diff(v, r) &
                  + (yz - diff(psih, z)/r)*diff(v, z) &
                  + (f + diff(yr, r) + diff(yz, z))*v) &
            - side_sum(neumann(psih, yr, yz, v, 2, flip), &
                       neumann(psih, yr, yz, v, 1, flip), &
                       neumann(psih, yr, yz, v, 3, flip), &
                       neumann(psih, yr, yz, v, 4, flip))
    end function residual

    function neumann(psih, yr, yz, v, side, flip) result(value)
        !! (g - S psi_h - y.n) v on side 2: R=2, 1: R=1, 3: Z=1, 4: Z=0, with
        !! g = d_nu psi + S psi.
        type(expr_t), intent(in) :: psih, yr, yz, v
        integer, intent(in) :: side
        logical, intent(in) :: flip
        type(expr_t) :: value, dn, yn
        select case (side)
        case (2)
            dn = diff(psi, r)/r; yn = yr
        case (1)
            dn = zero - diff(psi, r)/r; yn = zero - yr
        case (3)
            dn = diff(psi, z)/r; yn = yz
        case default
            dn = zero - diff(psi, z)/r; yn = zero - yz
        end select
        if (flip) yn = zero - yn
        value = (dn + sigma*psi - sigma*psih - yn)*v
    end function neumann

    function box(integrand) result(value)
        !! int_1^2 int_0^1 integrand dZ dR, exact.
        type(expr_t), intent(in) :: integrand
        type(expr_t) :: value, inner
        inner = integral(integrand, z, num(arena, 0), num(arena, 1))
        value = integral(inner, r, num(arena, 1), num(arena, 2))
    end function box

    function side_sum(on_r2, on_r1, on_z1, on_z0) result(value)
        !! Boundary integral of a side-dependent integrand, ds = dZ or dR.
        type(expr_t), intent(in) :: on_r2, on_r1, on_z1, on_z0
        type(expr_t) :: value
        value = integral(subs(on_r2, r, num(arena, 2)), z, num(arena, 0), &
                         num(arena, 1)) &
            + integral(subs(on_r1, r, num(arena, 1)), z, num(arena, 0), num(arena, 1)) &
            + integral(subs(on_z1, z, num(arena, 1)), r, num(arena, 1), num(arena, 2)) &
            + integral(subs(on_z0, z, num(arena, 0)), r, num(arena, 1), num(arena, 2))
    end function side_sum

    function integral(integrand, var, lo, hi) result(value)
        type(expr_t), intent(in) :: integrand, var, lo, hi
        type(expr_t) :: value
        logical :: ok
        character(:), allocatable :: why
        call definite_integral(arena, integrand, var, lo, hi, value, ok, why)
        if (.not. ok) then
            print *, why
            error stop 'FortSym definite integral refused'
        end if
    end function integral

    subroutine require_nonzero(label, expression)
        character(*), intent(in) :: label
        type(expr_t), intent(in) :: expression
        type(engine_result_t) :: result
        result = engine%zero_test(expression)
        if (.not. result%ok) error stop 'Native nonzero control undecided'
        if (result%verdict /= VERDICT_FALSE) error stop 'Native omission control failed'
        print *, 'PASS exact negative control: ', label
    end subroutine require_nonzero

    function parsed(text) result(expression)
        character(*), intent(in) :: text
        type(expr_t) :: expression
        logical :: ok
        character(:), allocatable :: message
        expression = parse_expr(arena, text, ok, message)
        if (.not. ok) then
            print *, message
            error stop 'Native FortSym parse failed'
        end if
    end function parsed
end program test_gs_coupled_majorant_identities
