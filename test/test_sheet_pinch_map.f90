program test_sheet_pinch_map
    ! Helical screw pinch with a sliding current-sheet interface at the
    ! (m, n) = (2, 1) rational of the continuous profile
    ! (docs/current-sheet-equilibrium.md sections 4 and 5.2). Cases:
    !   derivatives  gradient/Hessian against differences; reference state
    !                stationary (with the B_theta and B_z jumps)
    !   linear_p2/p3 linear limit (central difference in the wall amplitude)
    !                against the screw-pinch oracle: r displacement, field and
    !                sheet current, design-order gates
    !   nonlinear    finite amplitude: self-convergence of the sheet current
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_sheet_slab_map, only: sheet_slab_map_t, SHEET_PINCH
    use kin6d_current_sheet_oracles, only: pinch_sheet_t
    implicit none
    character(32) :: arg
    real(dp), parameter :: jump = 0.5_dp
    call get_command_argument(1, arg)
    select case (trim(arg))
    case ('derivatives')
        call derivatives()
    case ('linear_p2')
        call ladder(2, [16, 32, 64, 128])
    case ('linear_p3')
        call ladder(3, [16, 32, 64, 128, 256])
    case ('nonlinear')
        call nonlinear()
    case ('profile')
        call profile()
    case ('spectrum')
        call spectrum()

    case default
        error stop 'usage: test_sheet_pinch_map derivatives|linear_p2|linear_p3|nonlinear'
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

    !> Oracle with c+- = c (1 -+ jump/2) about the continuous (2,1) resonance
    !> (q(sigma-) < 2 < q(sigma+), no resonance inside either region).
    subroutine oracle(o, delta, jmp)
        type(pinch_sheet_t), intent(out) :: o
        real(dp), intent(in) :: delta
        real(dp), intent(in), optional :: jmp
        type(pinch_sheet_t) :: c
        real(dp) :: jj
        jj = jump
        if (present(jmp)) jj = jmp
        c%m = 2; c%k = -1.0_dp/3; c%cm = 0.1_dp; c%cp = 0.1_dp; c%bzm = 1
        c%r0 = 0.2_dp; c%b = 1
        o = c
        o%sigma = c%resonant_radius()
        o%cm = c%cm*(1 + jj/2); o%cp = c%cm*(1 - jj/2)
        o%delta = delta
        call o%init()
        call check(o%status == 0, 'oracle')
    end subroutine

    subroutine setup(m, o, order, ne, modes, amp)
        type(sheet_slab_map_t), intent(inout) :: m
        type(pinch_sheet_t), intent(in) :: o
        integer, intent(in) :: order, ne, modes
        real(dp), intent(in) :: amp
        m%geometry = SHEET_PINCH; m%slide = .true.
        m%order = order; m%ne = ne; m%modes = modes
        m%r0 = o%r0; m%b = o%b; m%sigma = o%sigma; m%mh = o%m; m%kz = o%k
        m%cm = o%cm; m%cp = o%cp; m%bzm = o%bzm; m%bzp = o%bzp; m%pm = 0; m%pp = 0
        m%wp = amp; m%wm = 0
        call m%init()
        call check(m%status == 0, 'init')
    end subroutine

    subroutine derivatives()
        type(sheet_slab_map_t) :: m
        type(pinch_sheet_t) :: o
        real(dp), allocatable :: g(:), ab(:, :), c(:, :), d(:, :), u0(:), gp(:), gm(:)
        real(dp), allocatable :: du(:), abx(:, :), cx(:, :), dx(:, :)
        real(dp) :: wp, wm, h, err, scale, hv
        integer :: i, j
        logical :: ok
        call oracle(o, 1.0_dp)
        call setup(m, o, 2, 4, 2, 0.0_dp)
        allocate (g(m%ndof), ab(3*m%band + 1, m%nb), c(m%nb, m%ng), d(m%ng, m%ng))
        allocate (abx, mold=ab); allocate (cx, mold=c); allocate (dx, mold=d)
        ! The reference state is an exact equilibrium; its discrete gradient is
        ! the quadrature error of the non-polynomial r-dependence and must
        ! vanish quickly with the quadrature order.
        call m%assemble(m%u, g, ab, c, d)
        hv = norm2(g)
        m%nq = 10
        call m%init()
        call m%assemble(m%u, g, ab, c, d)
        write (*, '(a,2es10.3)') 'reference-state gradient (nq 5, 10):', hv, norm2(g)
        call check(norm2(g) < 1e-6_dp*hv .or. norm2(g) < 1e-14_dp, &
                   'reference state stationary')
        call setup(m, o, 2, 4, 2, 0.02_dp)
        allocate (u0(m%ndof), gp(m%ndof), gm(m%ndof), du(m%ndof))
        do i = 1, m%ndof
            u0(i) = 0.003_dp*sin(1.7_dp*i)
        end do
        call m%assemble(u0, g, ab, c, d)
        where (g == 0) u0 = 0
        call m%assemble(u0, g, ab, c, d)
        h = 1e-6_dp; err = 0; scale = maxval(abs(g))
        do i = 1, m%ndof
            if (g(i) == 0) cycle
            du = 0; du(i) = h
            call m%total_energy(u0 + du, wp, ok)
            call m%total_energy(u0 - du, wm, ok)
            err = max(err, abs((wp - wm)/(2*h) - g(i)))
        end do
        write (*, '(a,es10.3,a,es10.3)') 'gradient vs differences:', err, ' of', scale
        call check(err < 1e-7_dp*scale, 'gradient')
        err = 0; scale = 0
        do j = 1, m%ndof
            if (g(j) == 0) cycle
            du = 0; du(j) = h
            call m%assemble(u0 + du, gp, abx, cx, dx)
            call m%assemble(u0 - du, gm, abx, cx, dx)
            do i = 1, m%ndof
                if (g(i) == 0) cycle
                if (i <= m%nb .and. j <= m%nb) then
                    if (abs(i - j) > m%band) cycle
                    hv = ab(2*m%band + 1 + i - j, j)
                else if (i <= m%nb) then
                    hv = c(i, j - m%nb)
                else if (j <= m%nb) then
                    hv = c(j, i - m%nb)
                else
                    hv = d(i - m%nb, j - m%nb)
                end if
                scale = max(scale, abs(hv))
                err = max(err, abs((gp(i) - gm(i))/(2*h) - hv))
            end do
        end do
        write (*, '(a,es10.3,a,es10.3)') 'Hessian vs differences:', err, ' of', scale
        call check(err < 1e-6_dp*scale, 'Hessian')
        write (*, '(a)') 'PASS derivatives'
    end subroutine

    subroutine errors(order, ne, exi, eb, ek, steps)
        integer, intent(in) :: order, ne
        real(dp), intent(out) :: exi, eb, ek
        integer, intent(out) :: steps
        type(sheet_slab_map_t) :: mp, mm
        type(pinch_sheet_t) :: o
        real(dp), parameter :: eps = 2e-5_dp
        real(dp), allocatable :: s1(:), s2(:), wt(:), y1p(:), bp(:, :), y1m(:), bm(:, :)
        real(dp) :: chi, dchi, xi, c, xie, bex(3), r1, n1, r2, n2, jp(3), jm(3), jt, jz
        integer :: i, side, st
        call oracle(o, 1.0_dp)
        call setup(mp, o, order, ne, 2, eps)
        call setup(mm, o, order, ne, 2, -eps)
        call mp%solve(); call mm%solve()
        call check(mp%status == 0 .and. mm%status == 0, 'linear solves converge')
        steps = max(mp%iterations, mm%iterations)
        call mp%sample(s1, s2, wt, y1p, bp)
        call mm%sample(s1, s2, wt, y1m, bm)
        r1 = 0; n1 = 0; r2 = 0; n2 = 0
        do i = 1, size(s1)
            side = merge(1, -1, s1(i) > o%sigma)
            call o%eval(s1(i), side, chi, dchi, xi)
            c = merge(o%cp, o%cm, side > 0)
            xie = xi*cos(s2(i))
            bex = [dchi*sin(s2(i)), (o%m*chi/s1(i) - xi*c/s1(i)**2)*cos(s2(i)), &
                   o%k*chi*cos(s2(i))]
            r1 = r1 + wt(i)*s1(i)*((y1p(i) - y1m(i))/(2*eps) - xie)**2
            n1 = n1 + wt(i)*s1(i)*xie**2
            r2 = r2 + wt(i)*s1(i)*sum(((bp(:, i) - bm(:, i))/(2*eps) - bex)**2)
            n2 = n2 + wt(i)*s1(i)*sum(bex**2)
        end do
        exi = sqrt(r1/n1); eb = sqrt(r2/n2)
        call mp%sheet_jump(jp, st)
        call check(st == 0, 'interface vertex')
        call mm%sheet_jump(jm, st)
        call o%sheet(jt, jz)
        ek = sqrt(((jp(2) - jm(2))/(2*eps) - jt)**2 + ((jp(3) - jm(3))/(2*eps) - jz)**2) &
             /sqrt(jt**2 + jz**2)
    end subroutine

    subroutine ladder(order, nes)
        integer, intent(in) :: order, nes(:)
        real(dp) :: x(size(nes)), b(size(nes)), k(size(nes))
        integer :: i, n, steps
        n = size(nes)
        write (*, '(a,i2,a,f5.2)') 'pinch P', order, ' jump', jump
        write (*, '(a)') '   ne  L2(xi_r)   L2(B)      |dK|/K     orders xi,B,K   newton'
        do i = 1, n
            call errors(order, nes(i), x(i), b(i), k(i), steps)
            if (i == 1) then
                write (*, '(i5,3es11.3,20x,i3)') nes(i), x(i), b(i), k(i), steps
            else
                write (*, '(i5,3es11.3,3f6.2,2x,i3)') nes(i), x(i), b(i), k(i), &
                    log(x(i - 1)/x(i))/log(2.0_dp), log(b(i - 1)/b(i))/log(2.0_dp), &
                    log(k(i - 1)/k(i))/log(2.0_dp), steps
            end if
        end do
        call check(log(x(n - 1)/x(n))/log(2.0_dp) > order + 1 - 0.4_dp, 'xi order')
        call check(log(b(n - 1)/b(n))/log(2.0_dp) > order - 0.3_dp, 'field order')
        call check(log(k(n - 1)/k(n))/log(2.0_dp) > order - 0.3_dp, 'sheet order')
        write (*, '(a)') 'PASS linear ladder'
    end subroutine

    subroutine nonlinear()
        type(sheet_slab_map_t) :: m
        type(pinch_sheet_t) :: o, ol
        integer, parameter :: nl = 4
        integer :: nes(nl), i, st, j, steps
        real(dp) :: jt(nl), jv(3), w(nl), amp, o1, o2, jtl, jzl
        nes = [16, 32, 64, 128]
        amp = 0.05_dp
        call oracle(o, 1.0_dp, 1.0_dp)
        call oracle(ol, amp, 1.0_dp)
        call ol%sheet(jtl, jzl)
        do i = 1, nl
            ! Amplitude continuation from the reference state.
            call setup(m, o, 2, nes(i), 4, 0.0_dp)
            steps = 0
            do j = 1, 5
                m%wp = amp*j/5
                call m%solve()
                steps = steps + m%iterations
                if (m%status /= 0) write (*, '(a,2i4,es10.3)') 'status', m%status, &
                    m%iterations, m%gnorm
                call check(m%status == 0, 'nonlinear solve')
            end do
            call m%sheet_jump(jv, st)
            jt(i) = jv(2); w(i) = m%energy
            write (*, '(a,i4,a,es23.15,a,es23.15,a,i4)') 'ne', nes(i), ' [[B_theta]]', &
                jt(i), ' W', w(i), ' newton', steps
        end do
        o1 = log(abs(jt(nl - 2) - jt(nl - 1))/abs(jt(nl - 1) - jt(nl)))/log(2.0_dp)
        o2 = log(abs(w(nl - 2) - w(nl - 1))/abs(w(nl - 1) - w(nl)))/log(2.0_dp)
        write (*, '(a,2f7.3,a,es12.4,a,es12.4)') 'self-convergence orders:', o1, o2, &
            ' linear:', jtl, ' nonlinear/linear - 1:', jt(nl)/jtl - 1
        call check(o1 > 1.7_dp .and. o2 > 3.5_dp, 'nonlinear self-convergence')
        write (*, '(a)') 'PASS nonlinear'
    end subroutine

    !> Diagnostic (not registered): cos(u) profiles against the oracle.
    subroutine profile()
        type(sheet_slab_map_t) :: m
        type(pinch_sheet_t) :: o
        real(dp), allocatable :: s1(:), s2(:), wt(:), y1(:), b(:, :)
        real(dp) :: chi, dchi, xi, cs, cb, amp, c
        integer :: i, j, side
        amp = 1e-5_dp
        call oracle(o, 1.0_dp)
        call setup(m, o, 2, 16, 2, amp)
        call m%solve()
        call m%sample(s1, s2, wt, y1, b)
        do i = 1, size(s1), m%n2*m%nq
            cs = 0; cb = 0
            do j = 0, m%n2 - 1
                cs = cs + 2*y1(i + j)*cos(s2(i + j))/m%n2
                cb = cb + 2*b(3, i + j)*cos(s2(i + j))/m%n2
            end do
            side = merge(1, -1, s1(i) > o%sigma)
            call o%eval(s1(i), side, chi, dchi, xi)
            c = merge(o%cp, o%cm, side > 0)
            write (*, '(f8.4,4es12.4)') s1(i), cs/amp, xi, cb/amp, o%k*chi
        end do
    end subroutine
    !> Diagnostic (not registered): lowest Hessian eigenvalues and their
    !> vectors at the reference state (near-marginal interface slip modes).
    subroutine spectrum()
        type(sheet_slab_map_t) :: m
        type(pinch_sheet_t) :: o
        real(dp), allocatable :: g(:), ab(:, :), c(:, :), d(:, :), h(:, :), ev(:), work(:)
        integer :: i, j, n, info, ne
        interface
            subroutine dsyev(jobz, uplo, n, a, lda, w, work, lwork, info)
                import :: dp
                character, intent(in) :: jobz, uplo
                integer, intent(in) :: n, lda, lwork
                real(dp), intent(inout) :: a(lda, *)
                real(dp), intent(out) :: w(*), work(*)
                integer, intent(out) :: info
            end subroutine
        end interface
        call oracle(o, 1.0_dp)
        do ne = 8, 16, 8
        call setup(m, o, 2, ne, 3, 0.0_dp)
        n = m%ndof
        allocate (g(n), ab(3*m%band + 1, m%nb), c(m%nb, m%ng), d(m%ng, m%ng), h(n, n), &
                  ev(n), work(10*n))
        call m%assemble(m%u, g, ab, c, d)
        h = 0
        do j = 1, m%nb
            do i = max(1, j - m%band), min(m%nb, j + m%band)
                h(i, j) = ab(2*m%band + 1 + i - j, j)
            end do
        end do
        h(1:m%nb, m%nb + 1:) = c; h(m%nb + 1:, 1:m%nb) = transpose(c)
        h(m%nb + 1:, m%nb + 1:) = d
        call dsyev('V', 'U', n, h, n, ev, work, 10*n, info)
        write (*, '(a,i4,a,6es11.3)') 'ne', ne, ' lowest eigenvalues', ev(1:6)
        do j = 1, 2
            write (*, '(a,i2,a,es10.3)') ' mode', j, ' largest entries (node, slot):', ev(j)
            do i = 1, n
                if (abs(h(i, j)) > 0.3_dp*maxval(abs(h(:, j)))) then
                    if (i <= m%nb) then
                        write (*, '(4x,i5,i4,es11.3)') (i - 1)/m%per + 1, mod(i - 1, m%per) + 1, h(i, j)
                    else
                        write (*, '(4x,a,i3,es11.3)') 'G', i - m%nb, h(i, j)
                    end if
                end if
            end do
        end do
        deallocate (g, ab, c, d, h, ev, work)
        end do
    end subroutine
end program test_sheet_pinch_map
