program test_sheet_slab_map
    ! Nonlinear frozen-in slab map with a current sheet (stage 2 of
    ! docs/current-sheet-equilibrium.md). Cases:
    !   derivatives  gradient/Hessian (incl. interface-shape modes) against
    !                differences; the reference state is stationary (also with
    !                a pressure jump)
    !   linear_p2/p3 sliding interface, asymmetric jump with pressure jump:
    !                linear limit (central difference in the wall amplitude)
    !                against the slab oracle, design-order gates
    !   symmetric    C0 map (no sliding) with a symmetric jump: design order
    !   floor        comparators: C0 map with the asymmetric jump, interface
    !                inside an element, continuous iota (Hahm-Kulsrud)
    !   nonlinear    finite amplitude: self-convergence of the sheet current
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_sheet_slab_map, only: sheet_slab_map_t
    use kin6d_current_sheet_oracles, only: slab_sheet_t
    implicit none
    character(32) :: arg
    real(dp) :: e1(3), e2(3), e3(3)
    call get_command_argument(1, arg)
    select case (trim(arg))
    case ('derivatives')
        call derivatives(.false.)
        call derivatives(.true.)
        write (*, '(a)') 'PASS derivatives'
    case ('linear_p2')
        call linear_ladder(2, [8, 16, 32, 64], 0.0_dp, 0.6_dp, 0.4_dp, .true., e1, e2, e3)
        call gates(2, e1, e2, e3)
    case ('linear_p3')
        call linear_ladder(3, [8, 16, 32, 64, 128], 0.0_dp, 0.6_dp, 0.4_dp, .true., e1, e2, e3)
        call gates(3, e1, e2, e3)
    case ('symmetric')
        call linear_ladder(2, [8, 16, 32, 64], 0.0_dp, 0.5_dp, 0.5_dp, .false., e1, e2, e3)
        call gates(2, e1, e2, e3)
    case ('floor')
        call floor_cases()
    case ('nonlinear')
        call nonlinear()
    case ('margin')
        call margin()

    case default
        error stop 'usage: test_sheet_slab_map derivatives|linear_p2|linear_p3|symmetric|floor|nonlinear'
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

    subroutine setup(m, order, ne, modes, shift, dplus, dminus, amp, slide)
        type(sheet_slab_map_t), intent(inout) :: m
        integer, intent(in) :: order, ne, modes
        real(dp), intent(in) :: shift, dplus, dminus, amp
        logical, intent(in) :: slide
        m%order = order; m%ne = ne; m%modes = modes; m%shift = shift
        m%dplus = dplus; m%dminus = dminus; m%slide = slide
        m%a = 1; m%k = acos(-1.0_dp); m%b0 = 1; m%bz = 2; m%p0 = 1
        m%wp = -amp; m%wm = 0.5_dp*amp
        call m%init()
        call check(m%status == 0, 'init')
    end subroutine

    real(dp) function hentry(m, ab, c, d, i, j) result(h)
        type(sheet_slab_map_t), intent(in) :: m
        real(dp), intent(in) :: ab(:, :), c(:, :), d(:, :)
        integer, intent(in) :: i, j
        h = 0
        if (i <= m%nb .and. j <= m%nb) then
            if (abs(i - j) <= m%band) h = ab(2*m%band + 1 + i - j, j)
        else if (i <= m%nb) then
            h = c(i, j - m%nb)
        else if (j <= m%nb) then
            h = c(j, i - m%nb)
        else
            h = d(i - m%nb, j - m%nb)
        end if
    end function

    subroutine derivatives(slide)
        logical, intent(in) :: slide
        type(sheet_slab_map_t) :: m
        real(dp), allocatable :: g(:), ab(:, :), c(:, :), d(:, :), u0(:), gp(:), gm(:)
        real(dp), allocatable :: du(:), abx(:, :), cx(:, :), dx(:, :)
        real(dp) :: wp, wm, h, err, scale, hv, sym
        integer :: i, j
        logical :: ok
        call setup(m, 2, 4, 2, 0.0_dp, 0.6_dp, 0.4_dp, 0.0_dp, slide)
        allocate (g(m%ndof), ab(3*m%band + 1, m%nb), c(m%nb, m%ng), d(m%ng, m%ng))
        allocate (abx, mold=ab); allocate (cx, mold=c); allocate (dx, mold=d)
        call m%assemble(m%u, g, ab, c, d)
        write (*, '(a,l2,a,es10.3)') 'slide', slide, ' reference-state gradient:', norm2(g)
        call check(norm2(g) < 1e-13_dp, 'reference state stationary')
        call setup(m, 2, 4, 2, 0.0_dp, 0.6_dp, 0.4_dp, 0.05_dp, slide)
        allocate (u0(m%ndof), gp(m%ndof), gm(m%ndof), du(m%ndof))
        do i = 1, m%ndof
            u0(i) = 0.01_dp*sin(1.7_dp*i)
        end do
        call m%assemble(u0, g, ab, c, d)
        do i = 1, m%nb      ! keep Dirichlet values at zero
            if (hentry(m, ab, c, d, i, i) == 1 .and. g(i) == 0) u0(i) = 0
        end do
        call m%assemble(u0, g, ab, c, d)
        h = 1e-5_dp; err = 0; scale = maxval(abs(g))
        do i = 1, m%ndof
            if (g(i) == 0) cycle
            du = 0; du(i) = h
            call m%total_energy(u0 + du, wp, ok)
            call m%total_energy(u0 - du, wm, ok)
            err = max(err, abs((wp - wm)/(2*h) - g(i)))
        end do
        write (*, '(a,es10.3,a,es10.3)') ' gradient vs differences:', err, ' of', scale
        call check(err < 1e-8_dp*scale, 'gradient')
        err = 0; scale = 0; sym = 0
        do j = 1, m%ndof
            if (g(j) == 0) cycle
            du = 0; du(j) = h
            call m%assemble(u0 + du, gp, abx, cx, dx)
            call m%assemble(u0 - du, gm, abx, cx, dx)
            do i = 1, m%ndof
                if (g(i) == 0) cycle
                hv = hentry(m, ab, c, d, i, j)
                scale = max(scale, abs(hv))
                err = max(err, abs((gp(i) - gm(i))/(2*h) - hv))
                sym = max(sym, abs(hv - hentry(m, ab, c, d, j, i)))
            end do
        end do
        write (*, '(a,es10.3,a,es10.3,a,es10.3)') ' Hessian vs differences:', err, &
            ' of', scale, ' asymmetry', sym
        call check(err < 1e-7_dp*scale .and. sym < 1e-12_dp*scale, 'Hessian')
    end subroutine

    !> Linear-limit errors of one mesh against the oracle.
    subroutine linear_errors(order, ne, shift, dplus, dminus, slide, exi, eb, ek, steps)
        integer, intent(in) :: order, ne
        real(dp), intent(in) :: shift, dplus, dminus
        logical, intent(in) :: slide
        real(dp), intent(out) :: exi, eb, ek
        integer, intent(out) :: steps
        type(sheet_slab_map_t) :: mp, mm
        type(slab_sheet_t) :: o
        real(dp), parameter :: eps = 2e-5_dp
        real(dp), allocatable :: s1(:), s2(:), wt(:), y1p(:), bp(:, :), y1m(:), bm(:, :)
        real(dp) :: psi, dpsi, xi, dxi, by, xie, bex(3), r1, n1, r2, n2, jyp, jym, jx
        integer :: i, side, st
        call setup(mp, order, ne, 2, shift, dplus, dminus, eps, slide)
        call setup(mm, order, ne, 2, shift, dplus, dminus, -eps, slide)
        call mp%solve(); call mm%solve()
        call check(mp%status == 0 .and. mm%status == 0, 'linear solves converge')
        steps = max(mp%iterations, mm%iterations)
        call mp%sample(s1, s2, wt, y1p, bp)
        call mm%sample(s1, s2, wt, y1m, bm)
        o%a = 1; o%k = mp%k; o%b0 = 1; o%dplus = dplus; o%dminus = dminus
        o%wp = -1; o%wm = 0.5_dp
        call o%init()
        call check(o%status == 0, 'oracle')
        r1 = 0; n1 = 0; r2 = 0; n2 = 0
        do i = 1, size(s1)
            side = merge(1, -1, s1(i) > 0)
            call o%eval(s1(i), side, psi, dpsi, xi, dxi)
            by = o%by(s1(i), side)
            xie = xi*cos(o%k*s2(i))
            bex = [-o%k*by*xi*sin(o%k*s2(i)), -by*dxi*cos(o%k*s2(i)), 0.0_dp]
            r1 = r1 + wt(i)*((y1p(i) - y1m(i))/(2*eps) - xie)**2
            n1 = n1 + wt(i)*xie**2
            r2 = r2 + wt(i)*sum(((bp(:, i) - bm(:, i))/(2*eps) - bex)**2)
            n2 = n2 + wt(i)*sum(bex**2)
        end do
        exi = sqrt(r1/n1); eb = sqrt(r2/n2)
        call mp%sheet_jump(jyp, jx, st)
        if (st == 0) then
            call mm%sheet_jump(jym, jx, st)
            ek = abs((jyp - jym)/(2*eps) - o%sheet())/abs(o%sheet())
        else
            ek = -1
        end if
    end subroutine

    subroutine linear_ladder(order, nes, shift, dplus, dminus, slide, exi, eb, ek)
        integer, intent(in) :: order, nes(:)
        real(dp), intent(in) :: shift, dplus, dminus
        logical, intent(in) :: slide
        real(dp), intent(out) :: exi(3), eb(3), ek(3)
        real(dp) :: x(size(nes)), b(size(nes)), s(size(nes)), oxi, ob, os
        integer :: i, steps, n
        write (*, '(a,i2,a,f5.2,a,2f6.3,a,l2)') 'P', order, ' shift', shift, ' jump', &
            dplus, dminus, ' slide', slide
        write (*, '(a)') '   ne  L2(xi)     L2(B)      |dK1|/K1   orders xi,B,K1  newton'
        n = size(nes)
        do i = 1, n
            call linear_errors(order, nes(i), shift, dplus, dminus, slide, x(i), b(i), &
                               s(i), steps)
            if (i == 1) then
                write (*, '(i5,3es11.3,20x,i3)') nes(i), x(i), b(i), s(i), steps
            else
                oxi = log(x(i - 1)/x(i))/log(2.0_dp)
                ob = log(b(i - 1)/b(i))/log(2.0_dp)
                os = 0
                if (s(i) > 0) os = log(s(i - 1)/s(i))/log(2.0_dp)
                write (*, '(i5,3es11.3,3f6.2,2x,i3)') nes(i), x(i), b(i), s(i), oxi, ob, &
                    os, steps
            end if
        end do
        exi = x(n - 2:n); eb = b(n - 2:n); ek = s(n - 2:n)
    end subroutine

    subroutine gates(order, exi, eb, ek)
        integer, intent(in) :: order
        real(dp), intent(in) :: exi(3), eb(3), ek(3)
        call check(log(exi(2)/exi(3))/log(2.0_dp) > order + 1 - 0.4_dp, &
                   'displacement order k+1')
        call check(log(eb(2)/eb(3))/log(2.0_dp) > order - 0.3_dp, 'field order k')
        call check(log(ek(2)/ek(3))/log(2.0_dp) > order - 0.3_dp, 'sheet current order')
        write (*, '(a)') 'PASS linear ladder'
    end subroutine

    subroutine floor_cases()
        real(dp) :: xa(3), ba(3), ka(3), x(3), b(3), k(3)
        call linear_ladder(2, [8, 16, 32, 64], 0.0_dp, 0.6_dp, 0.4_dp, .true., xa, ba, ka)
        call linear_ladder(2, [8, 16, 32, 64], 0.0_dp, 0.6_dp, 0.4_dp, .false., x, b, k)
        call check(log(b(2)/b(3))/log(2.0_dp) < 1 .and. b(3) > 10*ba(3), &
                   'C0 map without sliding: floor')
        call linear_ladder(2, [8, 16, 32, 64], 0.37_dp, 0.6_dp, 0.4_dp, .false., x, b, k)
        call check(log(b(2)/b(3))/log(2.0_dp) < 1 .and. b(3) > 10*ba(3), &
                   'interface inside an element: floor')
        call linear_ladder(2, [8, 16, 32, 64], 0.0_dp, 0.0_dp, 0.0_dp, .true., x, b, k)
        call check(log(b(2)/b(3))/log(2.0_dp) < 1 .and. b(3) > 10*ba(3), &
                   'continuous iota: floor')
        write (*, '(a)') 'PASS floor'
    end subroutine

    subroutine nonlinear()
        type(sheet_slab_map_t) :: m
        type(slab_sheet_t) :: o
        integer, parameter :: nl = 5
        integer :: nes(nl), i, st
        real(dp) :: jy(nl), jx, amp, w(nl), o1, o2
        nes = [8, 16, 32, 64, 128]
        amp = 0.1_dp
        o%a = 1; o%k = acos(-1.0_dp); o%b0 = 1; o%dplus = 0.6_dp; o%dminus = 0.4_dp
        o%wp = -amp; o%wm = 0.5_dp*amp
        call o%init()
        do i = 1, nl
            call setup(m, 2, nes(i), 6, 0.0_dp, 0.6_dp, 0.4_dp, amp, .true.)
            call m%solve()
            call check(m%status == 0, 'nonlinear solve')
            call m%sheet_jump(jy(i), jx, st)
            w(i) = m%energy
            write (*, '(a,i4,a,es23.15,a,es23.15,a,i3)') 'ne', nes(i), ' K1', jy(i), &
                ' W', w(i), ' newton', m%iterations
        end do
        o1 = log(abs(jy(nl - 2) - jy(nl - 1))/abs(jy(nl - 1) - jy(nl)))/log(2.0_dp)
        o2 = log(abs(w(nl - 2) - w(nl - 1))/abs(w(nl - 1) - w(nl)))/log(2.0_dp)
        write (*, '(a,2f7.3,a,es12.4)') 'self-convergence orders (K1, W):', o1, o2, &
            '  linear K1:', o%sheet()
        write (*, '(a,es12.4)') 'nonlinear/linear - 1:', jy(nl)/o%sheet() - 1
        call check(o1 > 1.8_dp .and. o2 > 3.5_dp, 'nonlinear self-convergence orders')
        write (*, '(a)') 'PASS nonlinear'
    end subroutine
    !> Diagnostic (not registered): largest wall amplitude reached by
    !> continuation before the map folds, symmetric jumps Delta.
    subroutine margin()
        type(sheet_slab_map_t) :: m
        real(dp) :: dl(4), amp, last, jy, jx
        integer :: i, ne, st, j
        dl = [0.4_dp, 0.2_dp, 0.1_dp, 0.05_dp]
        do j = 1, 2
            ne = 32*j
            do i = 1, 4
                call setup(m, 2, ne, 6, 0.0_dp, dl(i), dl(i), 0.0_dp, .true.)
                m%maxit = 40
                last = 0; amp = 0
                do
                    amp = amp + 0.02_dp
                    m%wp = -amp; m%wm = amp
                    call m%solve()
                    if (m%status /= 0) write (*, '(a,i3,a,i4,a,es10.3,a,f6.3)') '  status', &
                        m%status, ' it', m%iterations, ' |g|', m%gnorm, ' amp', amp
                    if (m%status /= 0 .or. amp > 2) exit
                    last = amp
                    call m%sheet_jump(jy, jx, st)
                end do
                write (*, '(a,i4,a,f6.3,a,f6.3,a,es11.3)') 'ne', ne, ' Delta', dl(i), &
                    ' last converged wall amplitude', last, ' K1 there', jy
            end do
        end do
    end subroutine
end program test_sheet_slab_map
