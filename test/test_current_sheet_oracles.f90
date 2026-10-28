program test_current_sheet_oracles
    ! Checks of the linear current-sheet references
    ! (docs/current-sheet-equilibrium.md section 5).
    !   slab_identities: ODE, wall, interface conditions; Hahm-Kulsrud sheet
    !                    current; symmetric-jump independence; Delta -> 0 rate.
    !   slab_fe:         independent P2 finite elements of the reduced energy
    !                    (F^2 xi'^2/k^2 + F^2 xi^2)/2, interface node at x = 0:
    !                    order >= 2.7 in xi, agreement of xi(0).
    !   screw_pinch:     Bessel Wronskian, Laplace residual, interface and wall
    !                    conditions, continuous-iota limit at rate ~1.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_current_sheet_oracles, only: slab_sheet_t, pinch_sheet_t, &
        modified_bessel_ik
    implicit none
    character(32) :: arg
    call get_command_argument(1, arg)
    select case (trim(arg))
    case ('slab_identities')
        call slab_identities()
    case ('slab_fe')
        call slab_fe()
    case ('screw_pinch')
        call screw_pinch()
    case default
        error stop 'usage: test_current_sheet_oracles slab_identities|slab_fe|screw_pinch'
    end select
contains
    subroutine check(ok, what)
        logical, intent(in) :: ok
        character(*), intent(in) :: what
        if (.not. ok) then
            write (*, '(a)') 'FAIL '//what
            error stop 1
        end if
    end subroutine

    subroutine slab_identities()
        type(slab_sheet_t) :: s, s2
        real(dp) :: x, h, r, psi, dpsi, xi, dxi, pp, pm, d1, d2, xp, xm
        real(dp) :: k1, k1hk, e(3), pi
        integer :: side, i
        pi = acos(-1.0_dp)
        s%a = 1; s%k = pi; s%b0 = 1; s%dplus = 0.2_dp; s%dminus = 0.05_dp
        s%wp = 0.01_dp; s%wm = 0.02_dp
        call s%init()
        call check(s%status == 0, 'slab init')
        ! ODE psi'' = k^2 psi by centred differences on both sides.
        h = 1e-3_dp
        do side = -1, 1, 2
            do i = 1, 4
                x = side*0.2_dp*i
                call s%eval(x + h, side, pp, d1, xi, dxi)
                call s%eval(x - h, side, pm, d2, xi, dxi)
                call s%eval(x, side, psi, dpsi, xi, dxi)
                r = (pp - 2*psi + pm)/h**2 - s%k**2*psi
                call check(abs(r) < 1e-5_dp*s%k**2*abs(psi), 'slab ODE')
                call check(abs((pp - pm)/(2*h) - dpsi) < 1e-5_dp*abs(dpsi) + 1e-12_dp, &
                           'slab derivative')
            end do
        end do
        ! Walls.
        call s%eval(s%a, 1, psi, dpsi, xi, dxi)
        call check(abs(xi - s%wp) < 1e-14_dp, 'slab wall +')
        call s%eval(-s%a, -1, psi, dpsi, xi, dxi)
        call check(abs(xi - s%wm) < 1e-14_dp, 'slab wall -')
        ! Interface: xi continuous, [[F^2 xi']] = 0.
        call s%eval(0.0_dp, 1, psi, dpsi, xp, d1)
        call s%eval(0.0_dp, -1, psi, dpsi, xm, d2)
        call check(abs(xp - xm) < 1e-14_dp .and. abs(xp - s%xi0) < 1e-14_dp, &
                   'slab xi continuity')
        call check(abs((s%k*s%by(0.0_dp, 1))**2*d1 - (s%k*s%by(0.0_dp, -1))**2*d2) &
                   < 1e-13_dp*abs(s%k**2*d1), 'slab total pressure')
        write (*, '(a,es23.15,a,es23.15)') 'jump case xi0 =', s%xi0, ' K1 =', s%sheet()
        ! Hahm-Kulsrud: K1 = 2 k b0 delta/sinh(ka), symmetric walls.
        s2%a = 1; s2%k = pi; s2%b0 = 1; s2%wp = -1e-2_dp; s2%wm = 1e-2_dp
        call s2%init()
        k1hk = 2*s2%k*s2%b0*1e-2_dp/sinh(s2%k*s2%a)
        k1 = s2%sheet()
        write (*, '(a,2es23.15)') 'Hahm-Kulsrud K1, closed form:', k1, k1hk
        call check(abs(k1 - k1hk) < 1e-14_dp*abs(k1hk), 'Hahm-Kulsrud sheet current')
        call s2%eval(0.0_dp, 1, psi, dpsi, xp, d1)
        call s2%eval(0.0_dp, -1, psi, dpsi, xm, d2)
        call check(xp*xm < 0 .and. abs(xp + xm) < 1e-14_dp, 'Hahm-Kulsrud step')
        ! Symmetric jump: K1 = 2 k b0 (1 + Delta) delta/sinh(ka) for every Delta.
        do i = 1, 3
            s2%dplus = 0.1_dp**i; s2%dminus = s2%dplus
            call s2%init()
            call check(abs(s2%xi0) < 1e-16_dp, 'symmetric xi0')
            call check(abs(s2%sheet() - k1hk*(1 + s2%dplus)) < 1e-13_dp*k1hk, &
                       'symmetric sheet independent of Delta')
        end do
        ! General walls: psi(0+-) = F(0+-) xi0 = O(Delta), shielding recovered.
        do i = 1, 3
            s2%dplus = 0.1_dp**(i + 1); s2%dminus = 0.5_dp*s2%dplus
            s2%wp = 0.01_dp; s2%wm = 0.02_dp
            call s2%init()
            call s2%eval(0.0_dp, 1, psi, dpsi, xi, dxi)
            e(i) = abs(psi)
        end do
        write (*, '(a,3es11.3)') 'psi(0+) for Delta = 1e-2,1e-3,1e-4:', e
        call check(abs(log10(e(1)/e(2)) - 1) < 0.05_dp .and. &
                   abs(log10(e(2)/e(3)) - 1) < 0.05_dp, 'shielding rate O(Delta)')
        ! Limit k a (Delta+ w+ + Delta- w-)/((Delta+ + Delta-) sinh ka).
        x = s2%k*s2%a*(s2%dplus*s2%wp + s2%dminus*s2%wm)/ &
            ((s2%dplus + s2%dminus)*sinh(s2%k*s2%a))
        write (*, '(a,2es23.15)') 'xi0 and its Delta -> 0 limit:', s2%xi0, x
        call check(abs(s2%xi0 - x) < 1e-3_dp*abs(x), 'xi0 limit')
        write (*, '(a)') 'PASS slab_identities'
    end subroutine

    subroutine slab_fe()
        type(slab_sheet_t) :: s
        integer, parameter :: nlev = 6
        integer :: lev, ne
        real(dp) :: err(nlev), x0h(nlev), order, pi
        pi = acos(-1.0_dp)
        s%a = 1; s%k = pi; s%b0 = 1; s%dplus = 0.2_dp; s%dminus = 0.05_dp
        s%wp = 0.01_dp; s%wm = 0.02_dp
        call s%init()
        do lev = 1, nlev
            ne = 8*2**(lev - 1)
            call fe_solve(s, ne, err(lev), x0h(lev))
            write (*, '(a,i4,a,es11.3,a,es23.15)') 'ne', ne, ' L2(xi) rel', err(lev), &
                ' xi(0)', x0h(lev)
        end do
        order = log(err(nlev - 1)/err(nlev))/log(2.0_dp)
        write (*, '(a,f6.3,a,es23.15)') 'order', order, ' oracle xi0', s%xi0
        call check(order > 2.7_dp, 'P2 order with interface node')
        call check(abs(x0h(nlev) - s%xi0) < 1e-4_dp*abs(s%xi0), 'interface displacement')
        write (*, '(a)') 'PASS slab_fe'
    end subroutine

    !> Independent P2 Galerkin solve; node at x = 0 when ne is even.
    subroutine fe_solve(s, ne, err, xi0h)
        type(slab_sheet_t), intent(in) :: s
        integer, intent(in) :: ne
        real(dp), intent(out) :: err, xi0h
        real(dp), allocatable :: kmat(:, :), u(:), rhs(:)
        real(dp) :: g(6), w(6), x, h, t, nb(3), db(3), f, xl, e2, n2
        real(dp) :: psi, dpsi, xie, dxie
        integer :: e, q, nd, idx(3), side, st
        call gauss6(g, w)
        nd = 2*ne + 1
        allocate (kmat(nd, nd), u(nd), rhs(nd))
        kmat = 0
        h = 2*s%a/ne
        do e = 1, ne
            xl = -s%a + (e - 1)*h
            side = merge(1, -1, xl + h/2 > 0)
            idx = [2*e - 1, 2*e, 2*e + 1]
            do q = 1, 6
                t = (g(q) + 1)/2; x = xl + t*h
                call shape(t, h, nb, db)
                f = s%k*s%by(x, side)
                kmat(idx, idx) = kmat(idx, idx) + w(q)*h/2*f**2* &
                    (spread(db, 2, 3)*spread(db, 1, 3)/s%k**2 + &
                     spread(nb, 2, 3)*spread(nb, 1, 3))
            end do
        end do
        u = 0; u(1) = s%wm; u(nd) = s%wp
        rhs = -kmat(:, 1)*u(1) - kmat(:, nd)*u(nd)
        call dense_solve(kmat(2:nd - 1, 2:nd - 1), rhs(2:nd - 1), st)
        call check(st == 0, 'fe solve')
        u(2:nd - 1) = rhs(2:nd - 1)
        xi0h = u(ne + 1)
        e2 = 0; n2 = 0
        do e = 1, ne
            xl = -s%a + (e - 1)*h
            side = merge(1, -1, xl + h/2 > 0)
            idx = [2*e - 1, 2*e, 2*e + 1]
            do q = 1, 6
                t = (g(q) + 1)/2; x = xl + t*h
                call shape(t, h, nb, db)
                call s%eval(x, side, psi, dpsi, xie, dxie)
                e2 = e2 + w(q)*h/2*(sum(nb*u(idx)) - xie)**2
                n2 = n2 + w(q)*h/2*xie**2
            end do
        end do
        err = sqrt(e2/n2)
    end subroutine

    pure subroutine shape(t, h, nb, db)
        real(dp), intent(in) :: t, h
        real(dp), intent(out) :: nb(3), db(3)
        nb = [2*(t - 0.5_dp)*(t - 1), -4*t*(t - 1), 2*t*(t - 0.5_dp)]
        db = [4*t - 3, -8*t + 4, 4*t - 1]/h
    end subroutine

    pure subroutine gauss6(g, w)
        real(dp), intent(out) :: g(6), w(6)
        g = [-0.9324695142031521_dp, -0.6612093864662645_dp, -0.2386191860831969_dp, &
             0.2386191860831969_dp, 0.6612093864662645_dp, 0.9324695142031521_dp]
        w = [0.1713244923791704_dp, 0.3607615730481386_dp, 0.4679139345726910_dp, &
             0.4679139345726910_dp, 0.3607615730481386_dp, 0.1713244923791704_dp]
    end subroutine

    subroutine dense_solve(a, b, status)
        real(dp), intent(inout) :: a(:, :), b(:)
        integer, intent(out) :: status
        integer :: n, i, j, p
        real(dp) :: f, row(size(b)), t
        n = size(b); status = 1
        do i = 1, n
            p = maxloc(abs(a(i:n, i)), 1) + i - 1
            if (a(p, i) == 0) return
            if (p /= i) then
                row = a(i, :); a(i, :) = a(p, :); a(p, :) = row
                t = b(i); b(i) = b(p); b(p) = t
            end if
            do j = i + 1, n
                f = a(j, i)/a(i, i)
                a(j, i:n) = a(j, i:n) - f*a(i, i:n)
                b(j) = b(j) - f*b(i)
            end do
        end do
        do i = n, 1, -1
            b(i) = (b(i) - sum(a(i, i + 1:n)*b(i + 1:n)))/a(i, i)
        end do
        status = 0
    end subroutine

    subroutine screw_pinch()
        type(pinch_sheet_t) :: c, s
        real(dp) :: x, bi, dbi, bk, dbk, chi, dchi, xi, cp, cm, xp, xm, h, r
        real(dp) :: lap, chip, chim, jt0, jz0, jt, jz, err(3), fp, fm, dpp, dpm
        integer :: m, i
        ! Wronskian I_m K_m' - I_m' K_m = -1/x.
        do m = 0, 4
            do i = 1, 8
                x = 0.25_dp*i
                call modified_bessel_ik(m, x, bi, dbi, bk, dbk)
                call check(abs(x*(bi*dbk - dbi*bk) + 1) < 1e-13_dp, 'Wronskian')
            end do
        end do
        ! I_1(1), K_1(1) reference values (A&S table 9.8).
        call modified_bessel_ik(1, 1.0_dp, bi, dbi, bk, dbk)
        call check(abs(bi - 0.565159103992485_dp) < 1e-14_dp, 'I_1(1)')
        call check(abs(bk - 0.601907230197235_dp) < 1e-14_dp, 'K_1(1)')
        ! Continuous iota, q = m/n at sigma: R0 = 3, m = 2, n = 1.
        c%m = 2; c%k = -1.0_dp/3; c%cm = 0.1_dp; c%cp = c%cm; c%bzm = 1
        c%r0 = 0.2_dp; c%b = 1; c%delta = 1e-3_dp
        c%sigma = c%resonant_radius()
        call c%init()
        call check(c%status == 0, 'continuous init')
        call c%eval(0.5_dp*(c%r0 + c%sigma), -1, chi, dchi, xi)
        call check(abs(chi) == 0, 'inner region shielded')
        call c%eval(c%sigma, 1, chip, dpp, xp)
        call c%sheet(jt0, jz0)
        write (*, '(a,f8.5,a,es12.4,a,2es14.6)') 'continuous: sigma', c%sigma, &
            ' xi step', xp, ' [[B1]] (theta,z)', jt0, jz0
        call c%eval(c%b, 1, chi, dchi, xi)
        call check(abs(xi - c%delta) < 1e-13_dp, 'continuous outer wall')
        ! Laplace residual of chi+ by differences.
        h = 1e-3_dp
        r = 0.9_dp
        call c%eval(r + h, 1, cp, dchi, xi)
        call c%eval(r - h, 1, cm, dchi, xi)
        call c%eval(r, 1, chi, dchi, xi)
        lap = (cp - 2*chi + cm)/h**2 + dchi/r - (c%m**2/r**2 + c%k**2)*chi
        call check(abs(lap) < 1e-5_dp*abs(chi)/r**2, 'Laplace residual')
        ! Jump family c+- = c (1 +- eps/2) about the continuous resonance,
        ! B_z+ from pressure balance; F(sigma+-) brackets zero.
        s = c
        do i = 1, 3
            s%cm = c%cm*(1 - 0.5_dp*0.1_dp**(i + 1))
            s%cp = c%cm*(1 + 0.5_dp*0.1_dp**(i + 1))
            call s%init()
            call check(s%status == 0 .and. .not. s%continuous, 'jump init')
            call s%eval(s%r0, -1, chi, dchi, xi)
            call check(abs(xi) < 1e-14_dp*c%delta, 'inner wall')
            call s%eval(s%b, 1, chi, dchi, xi)
            call check(abs(xi - s%delta) < 1e-13_dp*s%delta, 'outer wall')
            call s%eval(s%sigma, -1, chim, dpm, xm)
            call s%eval(s%sigma, 1, chip, dpp, xp)
            call check(abs(xp - xm) < 1e-10_dp*abs(xp), 'xi continuity')
            fp = s%f(s%sigma, 1); fm = s%f(s%sigma, -1)
            call check(fp*fm < 0, 'iota jump brackets the rational')
            call check(abs(fp*chip - fm*chim - xp*(s%cp**2 - s%cm**2)/s%sigma**3) &
                       < 1e-9_dp*(abs(fp*chip) + abs(fm*chim)), 'Lagrangian total pressure')
            call s%sheet(jt, jz)
            err(i) = abs(jz - jz0)/abs(jz0) + abs(jt - jt0)/abs(jt0)
            write (*, '(a,es9.2,a,es12.4,a,2es14.6,a,es10.3)') 'eps', 0.1_dp**(i + 1), &
                ' xi_s', xp, ' [[B_L]]', jt, jz, ' rel.diff', err(i)
        end do
        call check(abs(log10(err(1)/err(2)) - 1) < 0.1_dp .and. &
                   abs(log10(err(2)/err(3)) - 1) < 0.1_dp, &
                   'sheet current converges to shielding at rate O(Delta)')
        write (*, '(a)') 'PASS screw_pinch'
    end subroutine
end program test_current_sheet_oracles
