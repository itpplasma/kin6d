program test_gs_forward_majorant_identities
    use fortsym_arena, only: arena_t
    use fortsym_expr, only: expr_t, sym, operator(+), operator(-), &
        operator(*), operator(/), operator(**)
    use fortsym_parse, only: parse_expr
    use fortsym_diff, only: diff
    use fortsym_engine_native, only: native_engine_t, make_native_engine
    use fortsym_engine, only: engine_result_t, VERDICT_FALSE
    use fortsym_check, only: suite_t, suite_begin, check_identity, suite_end
    implicit none
    type(arena_t), target :: arena
    type(native_engine_t) :: engine
    type(suite_t) :: suite
    type(expr_t) :: r, z, alpha, lambda, phi, eps, bubble, exact, candidate
    type(expr_t) :: qr, qz, source, residual, density, perturbed
    character(64) :: rejected_premise
    call arena%init()
    engine = make_native_engine(arena)
    r = sym(arena, 'R'); z = sym(arena, 'Z')
    alpha = sym(arena, 'alpha'); lambda = sym(arena, 'lambda')
    phi = parsed('c0+c1*R+c2*Z+c3*R*Z+c4*R**2*Z**2')
    eps = sym(arena, 'eps')
    bubble = parsed('(R-1)*(2-R)*Z*(1-Z)')
    exact = parsed('-R**4+5*R**2-4')
    candidate = exact + alpha*bubble
    qr = diff(exact, r)/r; qz = diff(exact, z)/r
    source = parsed('8*R')
    call suite_begin(suite, 'Exact weighted GS energy majorant controls')
    call get_command_argument(1, rejected_premise)
    if (len_trim(rejected_premise) > 0) then
        select case (trim(rejected_premise))
        case ('jump')
            call check_identity(suite, engine, 'invalid omitted normal jump', &
                parsed('4/3'))
        case ('boundary')
            call check_identity(suite, engine, 'invalid matching boundary trace', &
                parsed('1/2'))
        case ('quadrature')
            call check_identity(suite, engine, &
                'invalid zero integral from Gauss pair', parsed('11/15120'))
        case default
            error stop 'Unknown rejected premise'
        end select
        call suite_end(suite, 'gs-forward-majorant-rejected-premise.json')
        stop
    end if
    call check_identity(suite, engine, 'actual fixed source divergence', &
        diff(qr, r) + diff(qz, z) + source)
    call check_identity(suite, engine, 'same trace bubble R1', &
        parsed('(1-1)*(2-1)*Z*(1-Z)'))
    call check_identity(suite, engine, 'same trace bubble R2', &
        parsed('(2-1)*(2-2)*Z*(1-Z)'))
    call check_identity(suite, engine, 'same trace bubble Z0', &
        parsed('(R-1)*(2-R)*0*(1-0)'))
    call check_identity(suite, engine, 'same trace bubble Z1', &
        parsed('(R-1)*(2-R)*1*(1-1)'))
    call check_identity(suite, engine, 'majorant optimal flux R mismatch', &
        qr - diff(candidate, r)/r + alpha*diff(bubble, r)/r)
    call check_identity(suite, engine, 'majorant optimal flux Z mismatch', &
        qz - diff(candidate, z)/r + alpha*diff(bubble, z)/r)
    call check_identity(suite, engine, 'R gradient square primitive', &
        diff(parsed('2*R**2-12*R+9*log(R)'), r) - parsed('(3-2*R)**2/R'))
    call check_identity(suite, engine, 'R bubble square primitive', &
        diff(parsed('R**4/4-2*R**3+13*R**2/2-12*R+4*log(R)'), r) &
        - parsed('((R-1)*(2-R))**2/R'))
    call check_identity(suite, engine, 'Z bubble square primitive', &
        diff(parsed('Z**3/3-Z**4/2+Z**5/5'), z) - parsed('(Z*(1-Z))**2'))
    call check_identity(suite, engine, 'Z gradient square primitive', &
        diff(parsed('Z-2*Z**2+4*Z**3/3'), z) - parsed('(1-2*Z)**2'))
    call check_identity(suite, engine, 'Green product identity with boundary flux', &
        diff(phi*diff(bubble, r)/r, r) + diff(phi*diff(bubble, z)/r, z) &
        - (diff(phi, r)*diff(bubble, r) + diff(phi, z)*diff(bubble, z))/r &
        - phi*(diff(diff(bubble, r)/r, r)+diff(diff(bubble, z)/r, z)))
    ! Nonlinear manufactured E0 source: spatial datum plus lambda*psi.
    qr = diff(bubble, r)/r; qz = diff(bubble, z)/r
    source = -diff(qr, r)-diff(qz, z)-lambda*bubble/r+lambda*alpha*bubble/r
    residual = source + diff(qr, r) + diff(qz, z)
    call check_identity(suite, engine, 'nonlinear source residual exact', &
        residual - lambda*(alpha-1)*bubble/r)
    perturbed = bubble + eps*phi
    density = (diff(perturbed, r)**2+diff(perturbed, z)**2)/(2*r) &
        -lambda*perturbed**2/(2*r)
    call check_identity(suite, engine, 'nonlinear second variation sign', &
        diff(diff(density, eps), eps) &
        -(diff(phi, r)**2+diff(phi, z)**2-lambda*phi**2)/r)
    ! Continuous tent across R=3/2 has opposite R gradients on its sides.
    call check_identity(suite, engine, 'interior outward jump is not omitted', &
        parsed('1/R-(-1/R)') - parsed('2/R'))
    ! u=Z is zero-source but has nonzero trace and nonzero energy on the box.
    call check_identity(suite, engine, 'boundary omission harmonic witness', &
        diff(diff(z, r)/r, r)+diff(diff(z, z)/r, z))
    call check_identity(suite, engine, 'boundary witness energy primitive', &
        diff(parsed('log(R)'), r) - parsed('1/R'))
    call check_identity(suite, engine, 'Gauss pair missed integral primitive', &
        diff(parsed('Z**7/7-Z**6/2+43*Z**5/60-13*Z**4/24'// &
        '+25*Z**3/108-Z**2/18+Z/144'), z) &
        - parsed('(Z-1/2)**2*(Z**2-Z+1/6)**2'))
    call check_identity(suite, engine, 'equal Gauss rules actual integral', &
        parsed('1/7-1/2+43/60-13/24+25/108-1/18+1/144-11/15120'))
    call require_nonzero('reject broken flux as globally Hdiv', parsed('4/3'))
    call require_nonzero('reject boundary trace omission', parsed('1/2'))
    call require_nonzero('reject quadrature zero-integral claim', parsed('11/15120'))
    call suite_end(suite, 'gs-forward-majorant-identities.json')
contains
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
end program test_gs_forward_majorant_identities
