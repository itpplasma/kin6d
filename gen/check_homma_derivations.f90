program check_homma_derivations
    !! Symbolic checks of the thermal-force derivations in Homma et al.,
    !! Contrib. Plasma Phys. 54 (2014), Eqs. (2), (3), (5), (6), (8), with the
    !! normalisation of thermal-force-benchmark/docs/conventions.md.
    !! Every claim is a FortSym zero test; integrals are exact FortSym definite
    !! integrals. No point evaluation.
    !! Dimensionless speed x = v/a, a^2 = 2 T_b/m_b, r = |x|, s = r^2.
    use, intrinsic :: iso_fortran_env, only: int64
    use fortsym_arena, only: arena_t
    use fortsym_expr, only: expr_t, sym, rat, func, exp, pi_expr, sqrt, oo_expr, operator(+), operator(-), &
        operator(*), operator(/), operator(**)
    use fortsym_diff, only: diff
    use fortsym_subs, only: subs
    use fortsym_engine, only: engine_result_t, VERDICT_TRUE
    use fortsym_engine_native, only: native_engine_t, make_native_engine
    use fortsym_defint, only: definite_integral
    use fortsym_print, only: print_expr
    use fortsym_string, only: chars
    implicit none
    type(arena_t), target :: arena
    type(native_engine_t) :: engine
    type(expr_t) :: x(3), q(3), s, r, sq, Ph, Xs, Ee, Xe, dPh, g, gam, coef, Qc, qx, br, hs
    type(expr_t) :: rr, Phr, Xr, dPr, coefr, Qcr, Th, closed
    type(expr_t) :: m, T, sm, sT, mb, ma, n, qa, qb, e0, lnl, B, kap, Fd, Fp, Kdrift
    type(expr_t) :: z, na, s2, cC, gab, ab, pref, Kt, tauT, tauB
    type(engine_result_t) :: rs
    integer :: i, j, k, nchk

    call arena%init()
    engine = make_native_engine(arena)
    nchk = 0
    do i = 1, 3
        x(i) = sym(arena, 'x' // achar(48 + i))
        q(i) = sym(arena, 'q' // achar(48 + i))
    end do
    Ee = sym(arena, 'Eerf'); Xe = sym(arena, 'Xgau')
    s = x(1)*x(1) + x(2)*x(2) + x(3)*x(3)
    r = sqrt(s); sq = sqrt(pi_expr(arena))
    Ph = func('erf', [r]); Xs = exp(-s); dPh = (2/sq)*Xs
    qx = q(1)*x(1) + q(2)*x(2) + q(3)*x(3)

    ! Identity 1: Poisson equations of the Maxwellian potentials. With
    ! H^M = n erf(r)/(a r) and G^M = n a [(r + 1/(2r)) erf(r) + exp(-r^2)/sqrt(pi)]
    ! these are lap H = -4 pi f and lap G = 2 H for f = n pi^(-3/2) a^(-3) exp(-r^2).
    call split(lap(Ph/r) + (4/sq)*Xs, 'Identity 1: lap(erf(r)/r) = -4 exp(-s)/sqrt(pi)')
    g = (r + 1/(2*r))*Ph + Xs/sq
    call split(lap(g) - 2*Ph/r, 'Identity 1: lap g = 2 erf(r)/r')

    ! Identity 2: Eq. (5). Its prefactor n/(16 pi) sqrt(m/2T) L equals
    ! (Gamma/2)(n/a) times the bracket/2, because L = 4 pi Gamma; the
    ! Hessian of g is that bracket/2.
    do i = 1, 3
        do j = 1, 3
            br = ((2/r - 1/(r*r*r))*Ph + dPh/(r*r))*kron(i, j) &
                + ((-2/(r*r*r) + 3/(r*r*r*r*r))*Ph - 3/(r*r*r*r)*dPh)*x(i)*x(j)
            call split(diff(diff(g, x(i)), x(j)) - br/2, 'Identity 2: Eq. (5) bracket = 2 Hess(g)')
        end do
    end do
    call zero(rat(arena, 1_int64, 16_int64)*(4*pi_expr(arena))/pi_expr(arena) &
        - rat(arena, 1_int64, 4_int64), 'Identity 2: (1/(16 pi)) L = Gamma/4 with L = 4 pi Gamma')

    ! Identity 3: Eq. (6) bracket is the Hessian of (q.x) gam(r) with
    ! gam = E/(2r^3) - Phi'/(2r^2), and gam' = r coef.
    gam = Ph/(2*r*r*r) - dPh/(2*r*r)
    coef = (-3/(2*r*r*r*r*r))*Ph + (1/(r*r) + 3/(2*r*r*r*r))*dPh
    Qc = (15/(2*r*r*r*r*r*r*r))*Ph + (-2/(r*r) - 5/(r*r*r*r) - 15/(2*r*r*r*r*r*r))*dPh
    do i = 1, 3
        do j = 1, 3
            hs = diff(diff(qx*gam, x(i)), x(j))
            br = coef*(q(i)*x(j) + q(j)*x(i) + qx*kron(i, j)) + Qc*qx*x(i)*x(j)
            call split(hs - br, 'Identity 3: Eq. (6) bracket = Hess[(q.x) gam]')
        end do
    end do

    ! Identities 4b, 4c: trace of the Eq. (6) Hessian is qx*Th(r) with
    ! Th = 5 coef + r^2 Qc. Th is proved as a one-variable identity in rr,
    ! which holds at rr = r(x) for all x, so the trace equals -(4/sqrt pi)(q.x)exp(-s).
    rr = sym(arena, 'rr')
    Phr = func('erf', [rr]); Xr = exp(-(rr*rr)); dPr = (2/sq)*Xr
    coefr = (-3/(2*rr*rr*rr*rr*rr))*Phr + (1/(rr*rr) + 3/(2*rr*rr*rr*rr))*dPr
    Qcr = (15/(2*rr*rr*rr*rr*rr*rr*rr))*Phr + (-2/(rr*rr) - 5/(rr*rr*rr*rr) - 15/(2*rr*rr*rr*rr*rr*rr))*dPr
    Th = 5*coefr + (rr*rr)*Qcr
    call split_r(Th + (4/sq)*Xr, 'Theta(r) = -(4/sqrt pi) exp(-r^2) for all r')
    closed = -(4/sq)*qx*Xs
    do k = 1, 3
        call zero(diff(closed, x(k)) + (4/sq)*Xs*(q(k) - 2*qx*x(k)), &
            'Identity 4b: d_k Tr D = -(4/sqrt pi) exp(-s)(q_k - 2 (q.x) x_k)')
    end do
    call zero(lap(closed) + (16/sq)*qx*(s - rat(arena, 5_int64, 2_int64))*Xs, &
        'Identity 4c: lap Tr D = -(16/sqrt pi)(q.x)(s - 5/2) exp(-s)')

    ! Identity 4d: prefactor of Eq. (6) is m L/(40 pi T^2) = m Gamma/(10 T^2).
    ! The Poisson chain lap^2 G^T = -8 pi f^T with f^T in the q-linear
    ! part of Eq. (1) fixes c_G = 2/(5T) in G^T = c_G (q.x) gam, and
    ! D^T = (Gamma/2) c_G/a^2 Hess(...).
    m = sym(arena, 'mb'); T = sym(arena, 'Tb')
    call zero((2/(5*T))*(m/(2*T))*rat(arena, 1_int64, 2_int64) - m/(10*T*T), &
        'Identity 4d: (Gamma/2) c_G/a^2 = m Gamma/(10 T^2)')
    call zero(m*4*pi_expr(arena)/(40*pi_expr(arena)*T*T) - m/(10*T*T), &
        'Identity 4d: m L/(40 pi T^2) = m Gamma/(10 T^2)')
    call zero((2/(5*T))*(m*m/(4*T*T))*(-16/sq) - (-8*pi_expr(arena))*(m*m/(5*T*T*T))/(pi_expr(arena)*sq), &
        'Identity 4d: biharmonic chain, c_G = 2/(5T)')
    ! The Poisson potential H^T = -(2m/(5 sqrt(pi) T^2)) (q.x) exp(-s) satisfies
    ! lap H^T = -4 pi f^T (same constants).
    call zero(lap(qx*Xs) - 4*qx*(s - rat(arena, 5_int64, 2_int64))*Xs, &
        'Identity 4d: lap[(q.x) exp(-s)] = 4 (q.x)(s - 5/2) exp(-s)')

    ! Identity 5: distorted-Maxwellian bracket 1 - (m/(n T^2))(1 - w^2/(5 v_th^2)) with
    ! w^2/(5 v_th^2) = 2 x^2/5 gives the conventions bracket (2 m/(5 n T^2))(x^2 - 5/2).
    call zero((1 - rat(arena, 2_int64, 5_int64)*s) + rat(arena, 2_int64, 5_int64)*(s - rat(arena, 5_int64, 2_int64)), &
        'Identity 5: 1 - 2x^2/5 = -(2/5)(x^2 - 5/2)')

    ! Identity 6: the Eq. (6) brackets are regular at r = 0. Their Laurent
    ! series from r**(-5) up shows every pole coefficient cancels, and the
    ! regular part fixes the r = 0 values; the Hessian of (q.x) exp(-s) at 0
    ! gives A^T(0) parallel to q.
    z = rat(arena, 0_int64, 1_int64)
    call laurent_zero(coefr, -(4/(5*sq)) + (4/(7*sq))*rr*rr, &
        'Identity 6: Eq. (6) scalar bracket = -4/(5 sqrt pi) + 4 r^2/(7 sqrt pi) + O(r^3)')
    call laurent_zero(Qcr, 8/(7*sq) - (8/(9*sq))*rr*rr, &
        'Identity 6: Eq. (6) x_i x_j bracket = 8/(7 sqrt pi) - 8 r^2/(9 sqrt pi) + O(r^3)')
    do k = 1, 3
        call zero(at0(diff(qx*Xs, x(k))) - q(k), 'Identity 6: grad_k[(q.x) exp(-s)] at x=0 equals q_k')
    end do

    ! Identity 7a: both terms of Eq. (2) follow from Eq. (3). sm = sqrt(m_b),
    ! sT = sqrt(T_b), a = sqrt(2 T_b/m_b), Gamma = L/(4 pi), grad_v = grad_x/a,
    ! H^M = n erf(r)/(a r) and H^T = Kt (q.x) exp(-s) (Identities 1, 4d).
    sm = sym(arena, 'sm'); sT = sym(arena, 'sT')
    mb = sm*sm; T = sT*sT; ma = sym(arena, 'ma'); n = sym(arena, 'n')
    qa = sym(arena, 'qa'); qb = sym(arena, 'qb'); e0 = sym(arena, 'e0')
    lnl = sym(arena, 'lnl'); B = sym(arena, 'B'); na = sym(arena, 'na')
    s2 = sqrt(rat(arena, 2_int64, 1_int64))
    kap = 1 + mb/ma
    cC = lnl*(qa*qb/e0)*(qa*qb/e0)*kap
    gab = qa*qa*qb*qb*lnl/(4*pi_expr(arena)*e0*e0*ma*ma)
    ab = s2*sT/sm
    pref = ma*(1 + ma/mb)*gab/ab
    Kt = -2*mb/(5*sq*T*T)
    do k = 1, 3
        Fd = pref*diff((n/ab)*Ph/r, x(k))
        Fp = -(cC/(4*pi_expr(arena)))*(n/(2*T))*(Ph - r*dPh)/(r*r*r)*x(k)
        call split(Fd - Fp, 'Identity 7a: printed F^0 of Eq. (2) = m_a A^M of Eq. (3)')
        Fd = pref*diff(Kt*qx*Xs, x(k))
        Fp = -(cC/(10*pi_expr(arena)*sq))*(sm/(s2*sT))/(T*T)*Xs*(q(k) - 2*qx*x(k))
        call zero(Fd - Fp, 'Identity 7a: printed F^gradT of Eq. (2) = m_a A^T of Eq. (3)')
        ! Control: C/(10 sqrt pi) without the factor pi is not Eq. (3).
        rs = engine%zero_test(Fd - pi_expr(arena)*Fp)
        if (rs%verdict == VERDICT_TRUE) error stop 'control: C/(10 sqrt pi) proved'
    end do
    ! Slow-limit drift: F^T(0) = pref Kt q, q_perp = kappa_wedge b x grad T with
    ! kappa_wedge = 5 n T/(2 q_b B), v = F x b/(q_a B) = Kdrift grad_perp T.
    Kdrift = pref*Kt*5*n*T/(2*qb*B*qa*B)
    tauT = 8*s2*pi_expr(arena)*e0*e0*sm*T*sT/(na*qa*qa*qb*qb*lnl)
    tauB = (3*sq/4)*tauT
    call zero(Kdrift - kap*v8(tauB), &
        'Identity 7a: slow-limit gc drift = (1 + m_b/m_a) x Eq. (8) with tau_B = (3 sqrt pi/4) tau_T')
    call zero(Kdrift - kap*(4/(3*sq))*v8(tauT), &
        'Identity 8: slow-limit gc drift = (1 + m_b/m_a) (4/(3 sqrt pi)) x Eq. (8) with tau_T')
    rs = engine%zero_test(Kdrift - kap*v8(tauT))
    if (rs%verdict == VERDICT_TRUE) error stop 'control: Eq. (8) with tau_T proved'

    call flr_ratio()
    print '(a,i0,a)', 'PASS: ', nchk, ' symbolic zero tests'

contains

    function v8(tau) result(v)
        !! Eq. (8) coefficient of grad_perp T_b for a collision time tau.
        type(expr_t), intent(in) :: tau
        type(expr_t) :: v
        v = -rat(arena, 3_int64, 2_int64)*mb*n/(qa*qb*B*B*na*tau)
    end function v8

    subroutine flr_ratio()
        !! Identity 7b. Test ion at x = X + rho, rho = (b x v)/Omega_a,
        !! T_b(x) = T_b + g.rho with g = (g1, g2, 0), b = e_z. The friction A^M
        !! taken at T_b(x) adds d_T A^M (g.rho); its mean over gyrophase and the
        !! test Maxwellian f_a (thermal speed aa) gives the FLR drift, compared
        !! with the gc drift of <F^T> from Eq. (2). Speeds use bb = sqrt(2 T_b/m_b),
        !! T_b = m_b bb^2/2, T_a = m_a aa^2/2. Exact for every mass and temperature.
        type(expr_t) :: aa, bb, vn, yv, Xv, nu, Tb, Ta, dnu, closed, g1, g2, ka
        type(expr_t) :: tau, sig, th, cc, v(3), vsq, gauss, expo, fM, jac, dT
        type(expr_t) :: kw, qv(3), qdv, w, fflr(2), fth(2), vflr(2), vth(2), slow
        integer :: l

        aa = sym(arena, 'aa'); bb = sym(arena, 'bb'); vn = sym(arena, 'vn')
        g1 = sym(arena, 'g1'); g2 = sym(arena, 'g2')
        tau = sym(arena, 'tau'); sig = sym(arena, 'sig'); th = sym(arena, 'th')
        Tb = mb*bb*bb/2; Ta = ma*aa*aa/2
        yv = vn/bb; Xv = exp(-yv**2)
        ! A^M = -nu(|v|) v, nu = (1 + m_a/m_b) Gamma n [erf(y) - y Phi'(y)]/|v|^3.
        nu = (1 + ma/mb)*gab*n*(func('erf', [yv]) - yv*(2/sq)*Xv)/(vn*vn*vn)
        ! d_T A^M = -(d nu/d bb)/(dT_b/d bb) v; the erf terms cancel.
        dnu = -diff(nu, bb)/diff(Tb, bb)
        closed = (1 + ma/mb)*gab*n*(4/sq)*Xv/(mb*bb*bb*bb*bb*bb)
        call split_atoms(dnu - closed, func('erf', [yv]), Xv, &
            'Identity 7b: d_T A^M = (1 + m_a/m_b) Gamma n Phiprime(y) v/(bb^3 T_b)')

        ! v = cc (tau cos th, tau sin th, sig), cc^2 = aa^2 bb^2/(aa^2 + bb^2),
        ! so exp(-|v|^2/aa^2) exp(-|v|^2/bb^2) = exp(-tau^2 - sig^2).
        cc = aa*bb/sqrt(aa*aa + bb*bb)
        v(1) = cc*tau*func('cos', [th]); v(2) = cc*tau*func('sin', [th]); v(3) = cc*sig
        vsq = v(1)*v(1) + v(2)*v(2) + v(3)*v(3)
        expo = -vsq/(aa*aa) - vsq/(bb*bb)
        call zero(expo + tau*tau + sig*sig, 'Identity 7b: cc scales both Gaussians to exp(-tau^2 - sig^2)')
        gauss = exp(-(tau*tau) - sig*sig)
        ! f_a without its Gaussian, times the Jacobian cc^3 tau of d^3v.
        fM = cc*cc*cc*tau/(pi_expr(arena)*sq*aa*aa*aa)
        ! With aa in place of cc the same measure covers f_a alone.
        call zero(average(tau*gauss/(pi_expr(arena)*sq)) - 1, &
            'Identity 7b: the test Maxwellian is normalised')

        ! FLR force m_a <closed v_k (g.rho)>, rho = (-v2, v1, 0) m_a/(q_a B).
        dT = (-g1*v(2) + g2*v(1))*ma/(qa*B)
        do l = 1, 2
            fflr(l) = average(ma*fM*gauss*subs(closed, Xv, z + 1)*v(l)*dT)
        end do
        ! <F^T>, F^T_k = pref Kt exp(-|v|^2/bb^2)(q_k - 2 (q.v) v_k/bb^2) from 7a,
        ! with q = kappa_wedge b x g.
        kw = 5*n*Tb/(2*qb*B)
        qv(1) = -kw*g2; qv(2) = kw*g1; qv(3) = z
        qdv = qv(1)*v(1) + qv(2)*v(2)
        w = (ma*(1 + ma/mb)*gab/bb)*(-2*mb/(5*sq*Tb*Tb))
        do l = 1, 2
            fth(l) = average(fM*gauss*w*(qv(l) - 2*qdv*v(l)/(bb*bb)))
        end do
        ! Drift F x b/(q_a B) = (F_2, -F_1)/(q_a B).
        vflr(1) = fflr(2)/(qa*B); vflr(2) = -fflr(1)/(qa*B)
        vth(1) = fth(2)/(qa*B); vth(2) = -fth(1)/(qa*B)
        ka = aa*aa + bb*bb
        do l = 1, 2
            call zero(vflr(l)*qa*Tb - qb*Ta*vth(l), &
                'Identity 7b: v_flr = (q_b/q_a)(T_a/T_b) v_gc, exact')
            ! Finite-mass Maxwellian average of the slow-limit drift:
            ! v_gc = v_gc(slow) (1 + T_a m_b/(T_b m_a))**(-5/2).
            if (l == 1) then
                slow = w*kw*g1/(qa*B)
            else
                slow = w*kw*g2/(qa*B)
            end if
            ! Multiplied by ka**(5/2) > 0 so that the half-integer powers of ka cancel.
            call zero(vth(l)*ka**rat(arena, 5_int64, 2_int64) - slow*bb**5, &
                'Identity 7b: <v_gc> = v_gc(0) (bb^2/(aa^2 + bb^2))^(5/2)')
        end do
    end subroutine flr_ratio

    function average(f) result(v)
        !! int_0^oo dtau int_-oo^oo dsig int_0^2pi dth f.
        type(expr_t), intent(in) :: f
        type(expr_t) :: v, th, sig, tau, u
        type(engine_result_t) :: ex
        logical :: ok
        character(:), allocatable :: why
        th = sym(arena, 'th'); sig = sym(arena, 'sig'); tau = sym(arena, 'tau')
        ex = engine%expand(f)
        if (.not. ex%ok) error stop 'average: expand failed'
        call definite_integral(arena, ex%value, th, z, 2*pi_expr(arena), u, ok, why)
        if (.not. ok) then
            print '(a,a)', 'gyrophase integral: ', why
            error stop 1
        end if
        call definite_integral(arena, u, sig, -oo_expr(arena), oo_expr(arena), v, ok, why)
        if (.not. ok) then
            print '(a,a)', 'parallel integral: ', why
            error stop 1
        end if
        u = v
        call definite_integral(arena, u, tau, z, oo_expr(arena), v, ok, why)
        if (.not. ok) then
            print '(a,a)', 'perpendicular integral: ', why
            error stop 1
        end if
    end function average

    subroutine laurent_zero(f, expected, label)
        !! Laurent coefficients of f at rr = 0 from order -5 to 2 equal expected.
        type(expr_t), intent(in) :: f, expected
        character(*), intent(in) :: label
        type(engine_result_t) :: sf, ls
        sf = engine%simplify(f)
        if (.not. sf%ok) error stop 'laurent_zero: simplify failed'
        ls = engine%laurent_series(sf%value, rr, z, -5, 2)
        if (.not. ls%ok) error stop 'laurent_zero: Laurent series failed'
        call zero(ls%value - expected, label)
    end subroutine laurent_zero

    subroutine split_atoms(f, atom_e, atom_x, label)
        !! f is linear in the atoms erf(.) and exp(.); each coefficient and the rest vanish.
        type(expr_t), intent(in) :: f, atom_e, atom_x
        character(*), intent(in) :: label
        type(expr_t) :: f1, cE, cX
        f1 = subs(subs(f, atom_e, Ee), atom_x, Xe)
        cE = diff(f1, Ee); cX = diff(f1, Xe)
        call zero(cE, label)
        call zero(cX, label)
        call zero(f1 - cE*Ee - cX*Xe, label)
    end subroutine split_atoms

    function kron(i, j) result(v)
        integer, intent(in) :: i, j
        type(expr_t) :: v
        if (i == j) then
            v = rat(arena, 1_int64, 1_int64)
        else
            v = rat(arena, 0_int64, 1_int64)
        end if
    end function kron

    function at0(f) result(v)
        type(expr_t), intent(in) :: f
        type(expr_t) :: v
        integer :: l
        v = f
        do l = 1, 3
            v = subs(v, x(l), rat(arena, 0_int64, 1_int64))
        end do
    end function at0

    function lap(f) result(v)
        type(expr_t), intent(in) :: f
        type(expr_t) :: v
        v = diff(diff(f, x(1)), x(1)) + diff(diff(f, x(2)), x(2)) + diff(diff(f, x(3)), x(3))
    end function lap

    subroutine split(f, label)
        !! f is linear in erf(r) and exp(-s); each coefficient and the remainder must be zero.
        type(expr_t), intent(in) :: f
        character(*), intent(in) :: label
        type(expr_t) :: f1, cE, cX, rest
        f1 = subs(subs(f, Ph, Ee), Xs, Xe)
        cE = diff(f1, Ee); cX = diff(f1, Xe)
        rest = f1 - cE*Ee - cX*Xe
        call zero(cE, label)
        call zero(cX, label)
        call zero(rest, label)
    end subroutine split

    subroutine split_r(f, label)
        !! Same split for a function of rr alone.
        type(expr_t), intent(in) :: f
        character(*), intent(in) :: label
        type(expr_t) :: f1, cE, cX, rest
        f1 = subs(subs(f, Phr, Ee), Xr, Xe)
        cE = diff(f1, Ee); cX = diff(f1, Xe)
        rest = f1 - cE*Ee - cX*Xe
        call zero(cE, label)
        call zero(cX, label)
        call zero(rest, label)
    end subroutine split_r

    subroutine zero(f, label)
        type(expr_t), intent(in) :: f
        character(*), intent(in) :: label
        type(engine_result_t) :: rs
        rs = engine%zero_test(f)
        if (rs%verdict /= VERDICT_TRUE) then
            print '(a,a)', 'Unproved: ', label
            print '(a,a)', '  residual: ', chars(print_expr(f))
            error stop 1
        end if
        nchk = nchk + 1
    end subroutine zero

end program check_homma_derivations
