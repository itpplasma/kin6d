module kin6d_sheet_slab_map
    ! Nonlinear frozen-in ideal map of the sheared slab with a current sheet
    ! (docs/current-sheet-equilibrium.md section 4, stage 2), mu0 = 1.
    !
    ! Labels s = (s1, s2, s3), -a < s1 < a, s2 periodic with length 2 pi/k,
    ! z invariant. Frozen b = (0, b2(s1), bz), b2 = b0 (s1/a + dplus) for
    ! s1 > 0 and b0 (s1/a - dminus) for s1 < 0 (frozen iota jump at the
    ! interface label s1 = 0), p = p0 - b2^2/2, so [[p + B^2/2]] = 0 in the
    ! reference state (a pressure jump when dplus /= dminus).
    !
    ! Physical position Y = Phi(X) with chart
    !   Phi = (x + (l0 + l1 x) G(y) + (m0 + m1 x) cos(k y), y)
    ! (gen/gen_current_sheet.f90). Walls are at x = +-a + w(+-) cos(k y).
    ! slide = .true. (sliding interface): per side sigma the chart maps x = 0
    ! to the interface curve x = G(y), G = sum_m G_m cos(m k y) unknown, and
    ! the wall; the computational interface node is duplicated with xi = 0 on
    ! both copies and independent eta, so the two sides slide along the shared
    ! interface (tangential discontinuity of the displacement). slide =
    ! .false.: one smooth chart (G = 0), X continuous everywhere (C0 map).
    !
    ! Discretization: continuous P_k Lagrange elements in s1 (`shift` moves
    ! interior vertices, non-sliding only), parity Fourier series in s2:
    ! xi = sum_{m=0..M} xi_m cos(m k s2), eta = sum_{m=1..M} eta_m sin(m k s2)
    ! (up-down symmetric walls; the parity removes the s2-relabelling null).
    ! Dirichlet xi = 0 at walls (and interface copies), eta free. Damped
    ! Newton with the exact Hessian: banded node block bordered by the G
    ! modes (LAPACK dgbsv, dgesv). [[p + B^2/2]] = 0 is the natural interface
    ! condition of the energy.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_current_sheet_leaves, only: sheet_map_energy, sheet_map_gradient, &
        sheet_map_hessian, sheet_map_position, sheet_map_field
    implicit none
    private
    public :: sheet_slab_map_t

    interface
        subroutine dgbsv(n, kl, ku, nrhs, ab, ldab, ipiv, b, ldb, info)
            import :: dp
            integer, intent(in) :: n, kl, ku, nrhs, ldab, ldb
            real(dp), intent(inout) :: ab(ldab, *), b(ldb, *)
            integer, intent(out) :: ipiv(*), info
        end subroutine
        subroutine dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
            import :: dp
            integer, intent(in) :: n, nrhs, lda, ldb
            real(dp), intent(inout) :: a(lda, *), b(ldb, *)
            integer, intent(out) :: ipiv(*), info
        end subroutine
    end interface

    type :: sheet_slab_map_t
        real(dp) :: a = 1, k = acos(-1.0_dp), b0 = 1, dplus = 0, dminus = 0
        real(dp) :: bz = 1, p0 = 1, wp = 0, wm = 0
        integer :: order = 2, ne = 16, modes = 2, nq = 0, n2 = 0
        real(dp) :: shift = 0
        logical :: slide = .false.
        integer :: maxit = 60
        real(dp) :: rtol = 1e-12_dp
        ! layout, state and results
        real(dp), allocatable :: vert(:), u(:)
        integer :: nnode = 0, per = 0, nb = 0, ng = 0, ndof = 0, band = 0
        integer :: status = -1, iterations = 0
        real(dp) :: energy = 0, gnorm = 0
    contains
        procedure, non_overridable :: init
        procedure, non_overridable :: solve
        procedure, non_overridable :: total_energy
        procedure, non_overridable :: assemble
        procedure, non_overridable :: sample
        procedure, non_overridable :: sheet_jump
        procedure, non_overridable :: interface_element
    end type
contains
    subroutine init(m)
        class(sheet_slab_map_t), intent(inout) :: m
        integer :: i
        real(dp) :: h
        m%status = 1
        if (m%order < 1 .or. m%order > 4 .or. m%ne < 2 .or. m%modes < 1) return
        if (abs(m%shift) >= 0.5_dp) return
        if (m%slide .and. (mod(m%ne, 2) /= 0 .or. m%shift /= 0)) return
        if (m%nq <= 0) m%nq = m%order + 3
        if (m%n2 <= 0) m%n2 = 4*m%modes + 4
        h = 2*m%a/m%ne
        if (allocated(m%vert)) deallocate (m%vert)
        allocate (m%vert(0:m%ne))
        do i = 0, m%ne
            m%vert(i) = -m%a + i*h
            if (i > 0 .and. i < m%ne) m%vert(i) = m%vert(i) + m%shift*h
        end do
        if (mod(m%ne, 2) == 0 .and. m%shift == 0) m%vert(m%ne/2) = 0
        m%per = 2*m%modes + 1
        m%nnode = m%order*m%ne + 1 + merge(1, 0, m%slide)
        m%nb = m%nnode*m%per
        m%ng = merge(m%modes + 1, 0, m%slide)
        m%ndof = m%nb + m%ng
        m%band = (m%order + 1)*m%per - 1
        if (allocated(m%u)) deallocate (m%u)
        allocate (m%u(m%ndof))
        m%u = 0
        m%status = 0
    end subroutine

    !> Element to the left of the interface vertex s1 = 0, or -1.
    integer function interface_element(m) result(e)
        class(sheet_slab_map_t), intent(in) :: m
        integer :: i
        e = -1
        do i = 1, m%ne - 1
            if (m%vert(i) == 0) e = i
        end do
    end function

    pure integer function node_of(m, e, l) result(j)
        class(sheet_slab_map_t), intent(in) :: m
        integer, intent(in) :: e, l
        j = (e - 1)*m%order + l + 1
        if (m%slide .and. e > m%ne/2) j = j + 1
    end function

    pure subroutine chart(m, e, l0, l1, m0, m1, side)
        class(sheet_slab_map_t), intent(in) :: m
        integer, intent(in) :: e
        real(dp), intent(out) :: l0, l1, m0, m1
        integer, intent(out) :: side
        side = merge(1, -1, m%vert(e - 1) + m%vert(e) > 0)
        if (m%slide) then
            l0 = 1; m0 = 0
            l1 = -side/m%a
            m1 = merge(m%wp, -m%wm, side > 0)/m%a
        else
            l0 = 0; l1 = 0
            m0 = (m%wp + m%wm)/2; m1 = (m%wp - m%wm)/(2*m%a)
        end if
    end subroutine

    pure subroutine lagrange(order, t, nb, db)
        integer, intent(in) :: order
        real(dp), intent(in) :: t
        real(dp), intent(out) :: nb(0:order), db(0:order)
        integer :: i, j, l
        real(dp) :: ti, tj, prod, s
        do i = 0, order
            ti = real(i, dp)/order
            prod = 1
            do j = 0, order
                if (j /= i) prod = prod*(t - real(j, dp)/order)/(ti - real(j, dp)/order)
            end do
            nb(i) = prod
            s = 0
            do l = 0, order
                if (l == i) cycle
                prod = 1/(ti - real(l, dp)/order)
                do j = 0, order
                    if (j == i .or. j == l) cycle
                    tj = real(j, dp)/order
                    prod = prod*(t - tj)/(ti - tj)
                end do
                s = s + prod
            end do
            db(i) = s
        end do
    end subroutine

    pure subroutine gauss(n, x, w)
        integer, intent(in) :: n
        real(dp), intent(out) :: x(n), w(n)
        integer :: i, it, j
        real(dp) :: z, p1, p2, p3, pp, pi
        pi = acos(-1.0_dp)
        do i = 1, n
            z = cos(pi*(i - 0.25_dp)/(n + 0.5_dp))
            do it = 1, 100
                p1 = 1; p2 = 0
                do j = 1, n
                    p3 = p2; p2 = p1
                    p1 = ((2*j - 1)*z*p2 - (j - 1)*p3)/j
                end do
                pp = n*(z*p1 - p2)/(z*z - 1)
                z = z - p1/pp
                if (abs(p1/pp) < 1e-16_dp) exit
            end do
            x(i) = -z
            w(i) = 2/((1 - z*z)*pp*pp)
        end do
    end subroutine

    !> Frozen b2 and p on the given side at label s1.
    pure subroutine frozen(m, s1, side, b2, p)
        class(sheet_slab_map_t), intent(in) :: m
        real(dp), intent(in) :: s1
        integer, intent(in) :: side
        real(dp), intent(out) :: b2, p
        if (side > 0) then
            b2 = m%b0*(s1/m%a + m%dplus)
        else
            b2 = m%b0*(s1/m%a - m%dminus)
        end if
        p = m%p0 - b2**2/2
    end subroutine

    !> Point state z (8) at element e, local t, label s2, and the dof
    !> Jacobian jz(local node dofs + G modes, 8); g2, g3 = G'', G''' at X2.
    subroutine point_state(m, u, e, t, s2, s1, z, jz, g2, g3, hlen)
        class(sheet_slab_map_t), intent(in) :: m
        real(dp), intent(in) :: u(:), t, s2
        integer, intent(in) :: e
        real(dp), intent(out) :: s1, z(8), jz(:, :), g2, g3, hlen
        real(dp) :: nb(0:m%order), db(0:m%order), cm, sm, mk
        integer :: l, mm, base, ld, nl
        nl = (m%order + 1)*m%per
        hlen = m%vert(e) - m%vert(e - 1)
        s1 = m%vert(e - 1) + t*hlen
        call lagrange(m%order, t, nb, db)
        db = db/hlen
        z = [s1, s2, 1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp]
        jz = 0
        do l = 0, m%order
            base = (node_of(m, e, l) - 1)*m%per
            do mm = 0, m%modes
                cm = cos(mm*m%k*s2); sm = sin(mm*m%k*s2); mk = mm*m%k
                ld = l*m%per + mm + 1
                jz(ld, 1) = nb(l)*cm; jz(ld, 3) = db(l)*cm; jz(ld, 5) = -nb(l)*mk*sm
                if (mm > 0) then
                    ld = l*m%per + m%modes + mm + 1
                    jz(ld, 2) = nb(l)*sm; jz(ld, 4) = db(l)*sm; jz(ld, 6) = nb(l)*mk*cm
                end if
            end do
            z(1:6) = z(1:6) + matmul(u(base + 1:base + m%per), &
                                     jz(l*m%per + 1:(l + 1)*m%per, 1:6))
        end do
        g2 = 0; g3 = 0
        if (m%ng > 0) then
            do mm = 0, m%modes
                mk = mm*m%k
                cm = cos(mk*z(2)); sm = sin(mk*z(2))
                z(7) = z(7) + u(m%nb + mm + 1)*cm
                z(8) = z(8) - u(m%nb + mm + 1)*mk*sm
                g2 = g2 - u(m%nb + mm + 1)*mk**2*cm
                g3 = g3 + u(m%nb + mm + 1)*mk**3*sm
                jz(nl + mm + 1, 7) = cm
                jz(nl + mm + 1, 8) = -mk*sm
            end do
            jz(1:nl, 7) = z(8)*jz(1:nl, 2)
            jz(1:nl, 8) = g2*jz(1:nl, 2)
        end if
    end subroutine

    !> Energy of state u; ok = .false. if det DY <= 0 at a quadrature point.
    subroutine total_energy(m, u, w, ok)
        class(sheet_slab_map_t), intent(in) :: m
        real(dp), intent(in) :: u(:)
        real(dp), intent(out) :: w
        logical, intent(out) :: ok
        real(dp) :: xq(m%nq), wq(m%nq), z(8), jz((m%order + 1)*m%per + m%ng, 8)
        real(dp) :: s1, s2, hl, b2, p, dens, jac, wt, l0, l1, m0, m1, g2, g3, per2
        integer :: e, q, j, side
        call gauss(m%nq, xq, wq)
        per2 = 2*acos(-1.0_dp)/m%k
        w = 0; ok = .true.
        do e = 1, m%ne
            call chart(m, e, l0, l1, m0, m1, side)
            do q = 1, m%nq
                do j = 1, m%n2
                    s2 = (j - 1)*per2/m%n2
                    call point_state(m, u, e, (xq(q) + 1)/2, s2, s1, z, jz, g2, g3, hl)
                    call frozen(m, s1, side, b2, p)
                    call sheet_map_energy(z(1), z(2), z(3), z(4), z(5), z(6), z(7), &
                        z(8), b2, m%bz, p, m%k, l0, l1, m0, m1, dens, jac)
                    if (jac <= 0) ok = .false.
                    wt = wq(q)*hl/2*per2/m%n2
                    w = w + wt*dens
                end do
            end do
        end do
    end subroutine

    !> Gradient, banded node Hessian ab (dgbsv layout), border c (nb x ng)
    !> and corner d (ng x ng); Dirichlet rows and columns are identity.
    subroutine assemble(m, u, grad, ab, c, d)
        class(sheet_slab_map_t), intent(in) :: m
        real(dp), intent(in) :: u(:)
        real(dp), intent(out) :: grad(:), ab(:, :), c(:, :), d(:, :)
        real(dp) :: xq(m%nq), wq(m%nq), z(8), jz((m%order + 1)*m%per + m%ng, 8)
        real(dp) :: s1, s2, hl, b2, p, g8(8), h8(36), hm(8, 8), wt, per2
        real(dp) :: jh((m%order + 1)*m%per + m%ng, 8), l0, l1, m0, m1, g2, g3, hv
        real(dp) :: cnode, cg, mk
        integer :: e, q, j, i, l, nl, kl, ii, jj, gi, gj, side, nloc, mm
        integer :: dof((m%order + 1)*m%per + m%ng)
        nl = (m%order + 1)*m%per
        nloc = nl + m%ng
        kl = m%band
        per2 = 2*acos(-1.0_dp)/m%k
        call gauss(m%nq, xq, wq)
        grad = 0; ab = 0; c = 0; d = 0
        do e = 1, m%ne
            call chart(m, e, l0, l1, m0, m1, side)
            do l = 0, m%order
                do ii = 1, m%per
                    dof(l*m%per + ii) = (node_of(m, e, l) - 1)*m%per + ii
                end do
            end do
            do ii = 1, m%ng
                dof(nl + ii) = m%nb + ii
            end do
            do q = 1, m%nq
                do j = 1, m%n2
                    s2 = (j - 1)*per2/m%n2
                    call point_state(m, u, e, (xq(q) + 1)/2, s2, s1, z, jz, g2, g3, hl)
                    call frozen(m, s1, side, b2, p)
                    call sheet_map_gradient(z(1), z(2), z(3), z(4), z(5), z(6), z(7), &
                        z(8), b2, m%bz, p, m%k, l0, l1, m0, m1, g8(1), g8(2), g8(3), &
                        g8(4), g8(5), g8(6), g8(7), g8(8))
                    call sheet_map_hessian(z(1), z(2), z(3), z(4), z(5), z(6), z(7), &
                        z(8), b2, m%bz, p, m%k, l0, l1, m0, m1, h8(1), h8(2), h8(3), &
                        h8(4), h8(5), h8(6), h8(7), h8(8), h8(9), h8(10), h8(11), &
                        h8(12), h8(13), h8(14), h8(15), h8(16), h8(17), h8(18), h8(19), &
                        h8(20), h8(21), h8(22), h8(23), h8(24), h8(25), h8(26), h8(27), &
                        h8(28), h8(29), h8(30), h8(31), h8(32), h8(33), h8(34), h8(35), &
                        h8(36))
                    l = 0
                    do jj = 1, 8
                        do ii = 1, jj
                            l = l + 1
                            hm(ii, jj) = h8(l); hm(jj, ii) = h8(l)
                        end do
                    end do
                    wt = wq(q)*hl/2*per2/m%n2
                    jh(1:nloc, :) = wt*matmul(jz(1:nloc, :), hm)
                    ! Second derivatives of gv, gd through X2 and the G modes.
                    cnode = wt*(g8(7)*g2 + g8(8)*g3)
                    do jj = 1, nloc
                        gj = dof(jj)
                        grad(gj) = grad(gj) + wt*dot_product(jz(jj, :), g8)
                        do ii = 1, nloc
                            gi = dof(ii)
                            hv = dot_product(jh(ii, :), jz(jj, :))
                            if (m%ng > 0) then
                                if (ii <= nl .and. jj <= nl) then
                                    hv = hv + cnode*jz(ii, 2)*jz(jj, 2)
                                else if (ii <= nl .and. jj > nl) then
                                    mm = jj - nl - 1; mk = mm*m%k
                                    cg = wt*(-g8(7)*mk*sin(mk*z(2)) &
                                             - g8(8)*mk**2*cos(mk*z(2)))
                                    hv = hv + cg*jz(ii, 2)
                                else if (ii > nl .and. jj <= nl) then
                                    mm = ii - nl - 1; mk = mm*m%k
                                    cg = wt*(-g8(7)*mk*sin(mk*z(2)) &
                                             - g8(8)*mk**2*cos(mk*z(2)))
                                    hv = hv + cg*jz(jj, 2)
                                end if
                            end if
                            if (ii <= nl .and. jj <= nl) then
                                ab(2*kl + 1 + gi - gj, gj) = ab(2*kl + 1 + gi - gj, gj) + hv
                            else if (ii <= nl) then
                                c(gi, gj - m%nb) = c(gi, gj - m%nb) + hv
                            else if (jj > nl) then
                                d(gi - m%nb, gj - m%nb) = d(gi - m%nb, gj - m%nb) + hv
                            end if
                        end do
                    end do
                end do
            end do
        end do
        do i = 1, m%nnode
            if (.not. dirichlet(i)) cycle
            do l = 1, m%modes + 1
                gi = (i - 1)*m%per + l
                grad(gi) = 0
                do jj = max(1, gi - kl), min(m%nb, gi + kl)
                    ab(2*kl + 1 + gi - jj, jj) = 0
                    ab(2*kl + 1 + jj - gi, gi) = 0
                end do
                ab(2*kl + 1, gi) = 1
                if (m%ng > 0) c(gi, :) = 0
            end do
        end do
    contains
        logical function dirichlet(node)
            integer, intent(in) :: node
            integer :: left
            dirichlet = node == 1 .or. node == m%nnode
            if (m%slide) then
                left = m%order*m%ne/2 + 1
                dirichlet = dirichlet .or. node == left .or. node == left + 1
            end if
        end function
    end subroutine

    !> Newton direction for the bordered system; info /= 0 on failure.
    subroutine newton_direction(m, grad, ab, c, d, dir, info)
        class(sheet_slab_map_t), intent(in) :: m
        real(dp), intent(in) :: grad(:), ab(:, :), c(:, :), d(:, :)
        real(dp), intent(out) :: dir(:)
        integer, intent(out) :: info
        real(dp), allocatable :: abt(:, :), rhs(:, :), s(:, :), r2(:, :)
        integer, allocatable :: piv(:), piv2(:)
        allocate (abt, source=ab)
        allocate (rhs(m%nb, 1 + m%ng), piv(m%nb))
        rhs(:, 1) = -grad(1:m%nb)
        if (m%ng > 0) rhs(:, 2:) = c
        call dgbsv(m%nb, m%band, m%band, 1 + m%ng, abt, size(ab, 1), piv, rhs, &
                   m%nb, info)
        if (info /= 0) return
        dir(1:m%nb) = rhs(:, 1)
        if (m%ng == 0) return
        allocate (s(m%ng, m%ng), r2(m%ng, 1), piv2(m%ng))
        s = d - matmul(transpose(c), rhs(:, 2:))
        r2(:, 1) = -grad(m%nb + 1:) - matmul(transpose(c), rhs(:, 1))
        call dgesv(m%ng, 1, s, m%ng, piv2, r2, m%ng, info)
        if (info /= 0) return
        dir(m%nb + 1:) = r2(:, 1)
        dir(1:m%nb) = rhs(:, 1) - matmul(rhs(:, 2:), r2(:, 1))
    end subroutine

    !> Damped Newton from the current state (identity after init).
    subroutine solve(m)
        class(sheet_slab_map_t), intent(inout) :: m
        real(dp), allocatable :: grad(:), ab(:, :), c(:, :), d(:, :), dir(:), trial(:)
        real(dp) :: w, wt, g0, t, slope, gprev
        integer :: it, info, ls, stalls
        logical :: ok, roundoff
        allocate (grad(m%ndof), ab(3*m%band + 1, m%nb), c(m%nb, m%ng), &
                  d(m%ng, m%ng), dir(m%ndof), trial(m%ndof))
        m%status = 2
        call m%total_energy(m%u, w, ok)
        if (.not. ok) return
        g0 = -1; gprev = huge(1.0_dp); stalls = 0
        do it = 1, m%maxit
            call m%assemble(m%u, grad, ab, c, d)
            m%gnorm = norm2(grad)
            if (g0 < 0) g0 = max(m%gnorm, tiny(1.0_dp))
            m%iterations = it - 1
            m%energy = w
            if (m%gnorm <= m%rtol*g0 .or. m%gnorm < 1e-15_dp) then
                m%status = 0
                return
            end if
            call newton_direction(m, grad, ab, c, d, dir, info)
            if (info /= 0) then
                m%status = 3
                return
            end if
            slope = dot_product(grad, dir)
            if (slope >= 0) then
                dir = -grad
                slope = dot_product(grad, dir)
            end if
            ! Roundoff band of W: the predicted change is below the rounding
            ! of W itself; take full Newton steps until the gradient stalls.
            roundoff = abs(slope) < 1e3_dp*epsilon(1.0_dp)*max(abs(w), 1.0_dp)
            if (roundoff) then
                if (m%gnorm > 0.5_dp*gprev) then
                    stalls = stalls + 1
                    if (stalls > 1) then
                        m%status = 0
                        return
                    end if
                end if
            end if
            gprev = m%gnorm
            t = 1
            do ls = 1, 40
                trial = m%u + t*dir
                call m%total_energy(trial, wt, ok)
                if (ok) then
                    if (roundoff .or. wt <= w + 1e-4_dp*t*slope) exit
                end if
                t = t/2
            end do
            if (ls > 40) then
                m%status = 4
                return
            end if
            m%u = trial
            w = wt
        end do
        m%iterations = m%maxit
    end subroutine

    !> Lagrangian samples at all quadrature points: s1, s2, weights,
    !> Y1 - s1 and B (Cartesian, at the labelled point).
    subroutine sample(m, s1, s2, wt, dy1, b)
        class(sheet_slab_map_t), intent(in) :: m
        real(dp), allocatable, intent(out) :: s1(:), s2(:), wt(:), dy1(:), b(:, :)
        real(dp) :: xq(m%nq), wq(m%nq), y(2), u22, bl(3), hl, per2
        integer :: e, q, j, n
        call gauss(m%nq, xq, wq)
        per2 = 2*acos(-1.0_dp)/m%k
        n = m%ne*m%nq*m%n2
        allocate (s1(n), s2(n), wt(n), dy1(n), b(3, n))
        n = 0
        do e = 1, m%ne
            hl = m%vert(e) - m%vert(e - 1)
            do q = 1, m%nq
                do j = 1, m%n2
                    n = n + 1
                    s2(n) = (j - 1)*per2/m%n2
                    call trace(m, e, (xq(q) + 1)/2, s2(n), s1(n), y, u22, bl)
                    wt(n) = wq(q)*hl/2*per2/m%n2
                    dy1(n) = y(1) - s1(n); b(:, n) = bl
                end do
            end do
        end do
    end subroutine

    !> Physical position y, dX2/ds2 and B at element e, local t, label s2.
    subroutine trace(m, e, t, s2, s1, y, u22, b)
        class(sheet_slab_map_t), intent(in) :: m
        integer, intent(in) :: e
        real(dp), intent(in) :: t, s2
        real(dp), intent(out) :: s1, y(2), u22, b(3)
        real(dp) :: z(8), jz((m%order + 1)*m%per + m%ng, 8), hl, b2, p, g2, g3
        real(dp) :: l0, l1, m0, m1
        integer :: side
        call chart(m, e, l0, l1, m0, m1, side)
        call point_state(m, m%u, e, t, s2, s1, z, jz, g2, g3, hl)
        call frozen(m, s1, side, b2, p)
        call sheet_map_position(z(1), z(2), z(7), m%k, l0, l1, m0, m1, y(1))
        y(2) = z(2); u22 = z(6)
        call sheet_map_field(z(1), z(2), z(3), z(4), z(5), z(6), z(7), z(8), b2, &
            m%bz, p, m%k, l0, l1, m0, m1, b(1), b(2), b(3))
    end subroutine

    !> Sheet current at the interface: cos(k y) amplitude of [[B_y]] and
    !> sin(k y) amplitude of [[B_x]] in the physical coordinate y, from the
    !> one-sided traces (each side projected on its own trace points).
    !> status 1 if the mesh has no interface vertex.
    subroutine sheet_jump(m, jy, jx, status)
        class(sheet_slab_map_t), intent(in) :: m
        real(dp), intent(out) :: jy, jx
        integer, intent(out) :: status
        real(dp) :: s1, s2, y(2), u22, b(3), per2, f
        integer :: j, e, side
        status = 1; jy = 0; jx = 0
        e = m%interface_element()
        if (e < 0) return
        per2 = 2*acos(-1.0_dp)/m%k
        do side = -1, 1, 2
            do j = 1, m%n2
                s2 = (j - 1)*per2/m%n2
                if (side < 0) then
                    call trace(m, e, 1.0_dp, s2, s1, y, u22, b)
                else
                    call trace(m, e + 1, 0.0_dp, s2, s1, y, u22, b)
                end if
                f = side*2*u22/m%n2
                jy = jy + f*b(2)*cos(m%k*y(2))
                jx = jx + f*b(1)*sin(m%k*y(2))
            end do
        end do
        status = 0
    end subroutine
end module kin6d_sheet_slab_map
