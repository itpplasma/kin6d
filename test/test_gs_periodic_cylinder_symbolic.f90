program periodic_cylinder_symbolic
    use fortsym_arena, only: arena_t
    use fortsym_expr, only: expr_t, sym, operator(+), operator(-), &
        operator(*), operator(/), operator(**)
    use fortsym_parse, only: parse_expr
    use fortsym_diff, only: diff
    use fortsym_engine_native, only: native_engine_t, make_native_engine
    use fortsym_check, only: suite_t, suite_begin, check_identity, suite_end
    implicit none
    type(arena_t), target :: arena
    type(native_engine_t) :: eng
    type(suite_t) :: suite
    type(expr_t) :: r, tau, length, twopi, q, major, mu0
    type(expr_t) :: bz, bt, lambda, phi, pol, curltheta, curlz, az
    type(expr_t) :: qtorus, iota, radial_force, qgh, qmap
    call arena%init()
    eng = make_native_engine(arena)
    r = sym(arena, 'r')
    tau = sym(arena, 'tau')
    length = sym(arena, 'L')
    twopi = sym(arena, 'P')
    q = sym(arena, 'q')
    major = sym(arena, 'Rmajor')
    mu0 = sym(arena, 'mu0')
    ! P denotes the exact symbolic constant 2*pi. L,q,tau are nonzero.
    ! r>0 is the punctured-cylinder algebra domain; axis limits are numeric/series.
    bz = parsed('Baxis/(1+tau**2*r**2)')
    bt = tau*r*bz
    lambda = parsed('2*tau/(1+tau**2*r**2)')
    phi = parsed('P*Baxis*log(1+tau**2*r**2)/(2*tau**2)')
    pol = parsed('L*Baxis*log(1+tau**2*r**2)/(2*tau)')
    az = parsed('-Baxis*log(1+tau**2*r**2)/(2*tau)')
    call suite_begin(suite, 'exact periodic-cylinder and leading-order map identities')
    curltheta = -diff(bz, r)
    curlz = diff(r*bt, r)/r
    call check_identity(suite, eng, 'cylinder curl theta', curltheta-lambda*bt)
    call check_identity(suite, eng, 'cylinder curl z', curlz-lambda*bz)
    radial_force = (curltheta*bz-curlz*bt)/mu0
    call check_identity(suite, eng, 'cylinder radial force', radial_force)
    call check_identity(suite, eng, 'radial magnetic stress', &
        diff((bz**2+bt**2)/(2*mu0), r)+bt**2/(mu0*r))
    call check_identity(suite, eng, 'axial flux primitive', diff(phi, r)-twopi*r*bz)
    call check_identity(suite, eng, 'meridional flux primitive', diff(pol, r)-length*bt)
    call check_identity(suite, eng, 'q is signed flux derivative ratio', &
        diff(phi, r)/diff(pol, r)-twopi/(length*tau))
    call check_identity(suite, eng, 'q is axial winding ratio', &
        twopi*r*bz/(length*bt)-twopi/(length*tau))
    call check_identity(suite, eng, 'explicit L prescribed twist', &
        (twopi/(length*q))*length*q-twopi)
    call check_identity(suite, eng, 'axial vector potential', diff(az, r)+bt)
    call check_identity(suite, eng, 'meridional flux potential', pol+length*az)
    qtorus = length*q/(twopi*major)
    iota = twopi*major/(length*q)
    call check_identity(suite, eng, 'finite-A RH winding scaling', &
        qtorus*twopi*major-length*q)
    call check_identity(suite, eng, 'finite-A inverse winding', qtorus*iota-1)
    call check_identity(suite, eng, 'finite-A pitch is held at fixed L', &
        qtorus*major-length*q/twopi)
    qgh = twopi/(length*tau)
    qmap = length*qgh/(twopi*major)
    call check_identity(suite, eng, 'finite-A physical field-line pitch', &
        r*bz/(major*bt)-qmap)
    call check_identity(suite, eng, 'canonical signed flux-derivative winding', &
        diff(phi, r)/(twopi*diff(major*az, r))+qmap)
    call check_identity(suite, eng, 'canonical minus RH transverse orientation', &
        diff(major*az, r)+major*bt)
    call suite_end(suite, 'periodic-cylinder-identities.json')
contains
    function parsed(text) result(e)
        character(*), intent(in) :: text
        type(expr_t) :: e
        character(:), allocatable :: message
        logical :: ok
        e = parse_expr(arena, text, ok, message)
        if (.not. ok) then
            print *, message
            error stop 1
        end if
    end function
end program
