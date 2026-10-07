program test_gs_variational_identities
    use fortsym_arena, only: arena_t
    use fortsym_expr, only: expr_t, sym, operator(+), operator(-), &
        operator(*), operator(/), operator(**)
    use fortsym_parse, only: parse_expr, parse_expr_in
    use fortsym_fparse, only: logical_lines, find_assignment
    use fortsym_string, only: str_t, chars
    use fortsym_dialect, only: dialect, DIA_FORTRAN
    use fortsym_diff, only: diff
    use fortsym_engine_native, only: native_engine_t, make_native_engine
    use fortsym_check, only: suite_t, suite_begin, check_identity, suite_end
    implicit none
    type(arena_t), target :: arena
    type(native_engine_t) :: eng
    type(suite_t) :: suite
    type(expr_t) :: r, z, eps, u, v, w, p, ff2, density, green
    type(expr_t) :: pr, pz, prr, pzz, f, fp, mu0, dstar, coef
    type(expr_t) :: b(3), j(3), force(3), coord(3), xi(3), a(3), cb(3)
    type(expr_t) :: db(3), surface(3), pp, divxi, gamma, internal
    type(expr_t) :: rho, prof, prof_s, prof_ss, native
    character(1024) :: source_path
    call arena%init()
    eng = make_native_engine(arena)
    r = sym(arena, 'R'); z = sym(arena, 'Z'); eps = sym(arena, 'eps')
    mu0 = sym(arena, 'mu0')
    call suite_begin(suite, 'GS strong weak and ideal-variation algebra')
    ! Independent symbolic jets retain arbitrary first and second derivatives.
    pr = sym(arena, 'psi_R'); pz = sym(arena, 'psi_Z')
    prr = sym(arena, 'psi_RR')
    pzz = sym(arena, 'psi_ZZ'); f = sym(arena, 'F'); fp = sym(arena, 'Fp')
    dstar = prr - pr/r + pzz
    b = [-pz/r, f/r, pr/r]
    ! Cylindrical curl in RH ordering (R,phi,Z), with F=F(psi).
    j = [-fp*pz/r, -dstar/r, fp*pr/r]
    force = cross(j, b)
    coef = -(dstar + f*fp)/r**2
    call check_identity(suite, eng, 'force R equals GS coefficient', force(1)-coef*pr)
    call check_identity(suite, eng, 'force Z equals GS coefficient', force(3)-coef*pz)
    call check_identity(suite, eng, 'toroidal force vanishes', force(2))
    ! Arbitrary-degree samples with exact symbolic coefficients: no float probes.
    u = parsed('c0+c1*R+c2*Z+c3*R**2+c4*R*Z+c5*Z**2+c6*R**3*Z')
    v = parsed('d0+d1*R+d2*Z+d3*R*Z+d4*R**2*Z**2')
    call check_identity(suite, eng, 'divergence form of Delta_star', &
        diff(diff(u,r)/r,r)+diff(diff(u,z)/r,z)-delta(u)/r)
    green = diff(v*diff(u,r)/r,r)+diff(v*diff(u,z)/r,z)
    call check_identity(suite, eng, 'pointwise Green product rule', green &
        -(diff(u,r)*diff(v,r)+diff(u,z)*diff(v,z))/r-v*delta(u)/r)
    ! Differentiated polynomial source families, including nonlinear sources.
    w = u + eps*v
    p = parsed('p0')+parsed('pa')*w+parsed('pb')*w**2
    ff2 = parsed('F0sq')+parsed('fa')*w+parsed('fb')*w**2
    density = (diff(w,r)**2+diff(w,z)**2)/(2*r)-mu0*r*p-ff2/(2*r)
    pp = parsed('pa')+2*parsed('pb')*w
    fp = (parsed('fa')+2*parsed('fb')*w)/2
    green = diff(v*diff(w,r)/r,r)+diff(v*diff(w,z)/r,z)
    call check_identity(suite, eng, 'GS functional first variation', &
        diff(density,eps)-green+v*(delta(w)/r+mu0*r*pp+fp/r))
    ! Three-dimensional ideal displacement: separate from scalar GS functional.
    coord = [sym(arena,'x'),sym(arena,'y'),sym(arena,'z')]
    a = [parsed('x*y*z'),parsed('x*y*z'),parsed('x*y*z')]
    b = curl(a)
    xi = [parsed('x**2+y'),parsed('y**2+z'),parsed('z**2+x')]
    a = cross(xi,b); db = curl(a); cb = curl(b)
    force = cross(cb,b); surface = cross(a,b)
    call check_identity(suite, eng, 'ideal magnetic variation local balance', &
        dot(b,db)-div(surface)+dot(xi,force))
    pp = parsed('p0+x**2+3*y**2+5*z**2')
    gamma = sym(arena,'gamma'); divxi = div(xi)
    internal = (-gamma*pp*divxi+pp*divxi)/(gamma-1)
    surface = [pp*xi(1), pp*xi(2), pp*xi(3)]
    call check_identity(suite, eng, 'material pressure variation local balance', &
        internal+div(surface)- &
        (xi(1)*diff(pp,coord(1))+xi(2)*diff(pp,coord(2))+xi(3)*diff(pp,coord(3))))
    rho = sym(arena,'rho')
    prof = parsed('a0+a1*rho**2+a2*rho**4+a3*rho**6')
    prof_s = parsed('a1+2*a2*rho**2+3*a3*rho**4')
    prof_ss = parsed('2*a2+6*a3*rho**2')
    call check_identity(suite, eng, 'toroidal flux to rho first derivative', &
        diff(prof,rho)-2*rho*prof_s)
    call check_identity(suite, eng, 'toroidal flux to rho second derivative', &
        diff(diff(prof,rho),rho)-2*prof_s-4*rho**2*prof_ss)
    call get_command_argument(1,source_path)
    if (len_trim(source_path)==0) error stop 'Require actual generated KIN source'
    native = actual('gs_source','source')
    call check_identity(suite, eng, 'executed KIN weak source weight', &
        native-parsed('mu0*r*pprime+ffprime/r'))
    native = actual('gs_field','b_r')
    call check_identity(suite, eng, 'executed KIN field R',native+parsed('psi_z/r'))
    native = actual('gs_field','b_z')
    call check_identity(suite, eng, 'executed KIN field Z',native-parsed('psi_r/r'))
    native = actual('gs_field','b_phi')
    call check_identity(suite, eng, 'executed KIN field phi',native-parsed('f/r'))
    call suite_end(suite, 'gs-variational-identities.json')
contains
    function actual(routine,symbol) result(e)
        character(*), intent(in) :: routine,symbol
        type(expr_t) :: e
        type(str_t), allocatable :: lines(:)
        logical :: good
        character(:), allocatable :: message,rhs,text
        integer :: n,k,first,last
        call logical_lines(trim(source_path),lines,n,good,message)
        if (.not. good) error stop 'Cannot read actual source'
        first=0; last=0
        do k=1,n
            text=trim(adjustl(chars(lines(k))))
            if (index(text,'subroutine '//routine//'(')>0) first=k
            if (first>0) then
                if (index(text,'end subroutine '//routine)==1) then
                    last=k
                    exit
                end if
            end if
        end do
        if (first==0) error stop 'Actual source routine absent'
        if (last<=first) error stop 'Actual source routine boundary absent'
        call find_assignment(lines(first:last),last-first+1,symbol,rhs,good,message)
        if (.not. good) error stop 'Actual source assignment absent'
        e=parse_expr_in(arena,rhs,dialect(DIA_FORTRAN),good,message)
        if (.not. good) then
            print *, message
            error stop 1
        end if
    end function
    function parsed(text) result(e)
        character(*), intent(in) :: text
        type(expr_t) :: e
        character(:), allocatable :: message
        logical :: good
        e = parse_expr(arena,text,good,message)
        if (.not. good) then
            print *, message
            error stop 1
        end if
    end function
    function delta(e) result(out)
        type(expr_t), intent(in) :: e
        type(expr_t) :: out
        out = diff(diff(e,r),r)-diff(e,r)/r+diff(diff(e,z),z)
    end function
    function cross(left,right) result(out)
        type(expr_t), intent(in) :: left(3),right(3)
        type(expr_t) :: out(3)
        out(1)=left(2)*right(3)-left(3)*right(2)
        out(2)=left(3)*right(1)-left(1)*right(3)
        out(3)=left(1)*right(2)-left(2)*right(1)
    end function
    function dot(left,right) result(out)
        type(expr_t), intent(in) :: left(3),right(3)
        type(expr_t) :: out
        out=left(1)*right(1)+left(2)*right(2)+left(3)*right(3)
    end function
    function curl(e) result(out)
        type(expr_t), intent(in) :: e(3)
        type(expr_t) :: out(3)
        out(1)=diff(e(3),coord(2))-diff(e(2),coord(3))
        out(2)=diff(e(1),coord(3))-diff(e(3),coord(1))
        out(3)=diff(e(2),coord(1))-diff(e(1),coord(2))
    end function
    function div(e) result(out)
        type(expr_t), intent(in) :: e(3)
        type(expr_t) :: out
        out=diff(e(1),coord(1))+diff(e(2),coord(2))+diff(e(3),coord(3))
    end function
end program
