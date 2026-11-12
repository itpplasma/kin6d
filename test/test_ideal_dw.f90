program test_ideal_dw
    ! Oracles for the delta-W leaves that do not reuse their symbolic proofs:
    ! integral equality of the Eulerian and Lagrangian forms on a force-balanced
    ! screw pinch (and its failure off equilibrium), the gamma = 0 limit against
    ! the independent xi form, Hermitian symmetry and JVP/VJP pairing.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_quadrature, only: gauss_legendre
    use milca_ideal_dw, only: dw_bernstein_form, dw_lagrange_form
    use milca_ideal_dw_leaves, only: dw_bernstein, dw_bernstein_jvp, &
        dw_bernstein_vjp
    use milca_ideal_xi, only: xi_bilinear
    implicit none
    real(dp), parameter :: mu = 1.25_dp, cur = 0.7_dp, bz0 = 1.3_dp, &
        p0 = 2.0_dp, kz = 1.5_dp, gam = 5.0_dp/3
    real(dp) :: balanced, unbalanced, scale
    call integrals(.true., balanced, scale)
    print '(a,es10.2)', 'force-balanced pinch: (L-E)/|L| = ', balanced/scale
    if (abs(balanced) > 1e-12_dp*scale) error stop 'Eulerian/Lagrangian forms differ'
    call integrals(.false., unbalanced, scale)
    print '(a,es10.2)', 'uniform pressure (no balance): (L-E)/|L| = ', unbalanced/scale
    if (abs(unbalanced) < 1e-3_dp*scale) error stop 'off-equilibrium test is blind'
    call gamma_zero_limit()
    call hermitian_symmetry()
    call derivative_pairing()
    print '(a)', 'PASS ideal delta-W leaves'
contains
    subroutine integrals(equilibrium, difference, magnitude)
        ! B = (-c y, c x, B0), J = 2c/mu e_z, p = p0 - c^2 r^2/mu balances J x B.
        logical, intent(in) :: equilibrium
        real(dp), intent(out) :: difference, magnitude
        integer, parameter :: nq = 24, nz = 8
        real(dp) :: t(nq), wt(nq), b(3), g(3, 3), p, gp(3), xi(3), a(3, 3), &
            eta(3), c(3, 3), x, y, z, w, lag, eul
        real(dp), parameter :: twopi = 2*acos(-1.0_dp)
        integer :: i, j, l
        call gauss_legendre(nq, t, wt)
        difference = 0; magnitude = 0
        g = 0; g(1, 2) = -cur; g(2, 1) = cur
        do l = 1, nz
            z = twopi/kz*(l - 1)/nz
            do j = 1, nq
                do i = 1, nq
                    x = t(i); y = t(j)
                    w = wt(i)*wt(j)*twopi/(kz*nz)
                    b = [-cur*y, cur*x, bz0]
                    if (equilibrium) then
                        p = p0 - cur**2*(x*x + y*y)/mu
                        gp = [-2*cur**2*x/mu, -2*cur**2*y/mu, 0.0_dp]
                    else
                        p = p0; gp = 0
                    end if
                    call field(x, y, z, [0.3_dp, -0.5_dp, 0.8_dp], &
                        [0.7_dp, 0.2_dp, -0.4_dp], [-0.6_dp, 0.9_dp, 0.1_dp], &
                        [0.0_dp, 0.4_dp, 1.1_dp], xi, a)
                    call field(x, y, z, [-0.2_dp, 0.6_dp, 0.5_dp], &
                        [0.1_dp, -0.8_dp, 0.3_dp], [0.5_dp, 0.3_dp, -0.7_dp], &
                        [0.9_dp, -0.3_dp, 0.2_dp], eta, c)
                    lag = dw_lagrange_form(b, p, gam, mu, a, c)
                    eul = dw_bernstein_form(b, g, p, gp, gam, mu, xi, a, eta, c)
                    difference = difference + w*(lag - eul)
                    magnitude = magnitude + w*abs(lag)
                end do
            end do
        end do
    end subroutine integrals

    pure subroutine field(x, y, z, alpha, beta, gamma, phase, v, dv)
        ! v_i = w(x,y) (alpha_i + beta_i x + gamma_i y) cos(kz z + phase_i),
        ! w = (1-x^2)^3 (1-y^2)^3 vanishes with its gradient on the box edge.
        real(dp), intent(in) :: x, y, z, alpha(3), beta(3), gamma(3), phase(3)
        real(dp), intent(out) :: v(3), dv(3, 3)
        real(dp) :: w, wx, wy, poly, cs, sn
        integer :: k
        w = (1 - x*x)**3*(1 - y*y)**3
        wx = -6*x*(1 - x*x)**2*(1 - y*y)**3
        wy = -6*y*(1 - x*x)**3*(1 - y*y)**2
        do k = 1, 3
            poly = alpha(k) + beta(k)*x + gamma(k)*y
            cs = cos(kz*z + phase(k)); sn = sin(kz*z + phase(k))
            v(k) = w*poly*cs
            dv(k, 1) = (wx*poly + w*beta(k))*cs
            dv(k, 2) = (wy*poly + w*gamma(k))*cs
            dv(k, 3) = -kz*w*poly*sn
        end do
    end subroutine field

    subroutine gamma_zero_limit()
        real(dp) :: b(3), a(3, 3), c(3, 3), lag
        complex(dp) :: ref
        integer :: i
        do i = 1, 5
            call sample(i, b, a, c)
            lag = dw_lagrange_form(b, 0.4_dp*i, 0.0_dp, mu, a, c)
            ref = xi_bilinear(b, 0.4_dp*i, mu, cmplx(a, 0, dp), cmplx(c, 0, dp))
            if (abs(lag - real(ref)) > 1e-13_dp*(1 + abs(lag))) &
                error stop 'gamma = 0 differs from the xi form'
        end do
    end subroutine gamma_zero_limit

    subroutine hermitian_symmetry()
        real(dp) :: b(3), g(3, 3), a(3, 3), c(3, 3), gp(3), ra(3, 3), rc(3, 3), &
            spare(3, 3)
        complex(dp) :: xi(3), eta(3), za(3, 3), zc(3, 3), f1, f2
        call sample(1, b, g, spare)
        call sample(2, gp, a, c)
        call sample(3, spare(:, 1), ra, rc)
        xi = cmplx(gp, b, dp); eta = cmplx(b, -gp, dp)
        za = cmplx(a, ra, dp); zc = cmplx(c, rc, dp)
        f1 = dw_bernstein_form(b, g, 1.7_dp, gp, gam, mu, xi, za, eta, zc)
        f2 = dw_bernstein_form(b, g, 1.7_dp, gp, gam, mu, eta, zc, xi, za)
        if (abs(f1 - conjg(f2)) > 1e-13_dp*abs(f1)) error stop 'not Hermitian'
    end subroutine hermitian_symmetry

    subroutine derivative_pairing()
        ! JVP against central differences and VJP . direction against JVP.
        real(dp) :: u(42), du(42), gr(42), up(42), um(42), jv, fp, fm, h
        integer :: i
        do i = 1, 42
            u(i) = 0.3_dp + 0.1_dp*sin(1.7_dp*i)
            du(i) = cos(0.9_dp*i)
        end do
        u(18) = 1.1_dp
        call dw_bernstein_jvp(u(1), u(2), u(3), u(4), u(5), u(6), u(7), u(8), &
            u(9), u(10), u(11), u(12), u(13), u(14), u(15), u(16), u(17), u(18), &
            u(19), u(20), u(21), u(22), u(23), u(24), u(25), u(26), u(27), u(28), &
            u(29), u(30), u(31), u(32), u(33), u(34), u(35), u(36), u(37), u(38), &
            u(39), u(40), u(41), u(42), du(1), du(2), du(3), du(4), du(5), du(6), &
            du(7), du(8), du(9), du(10), du(11), du(12), du(13), du(14), du(15), &
            du(16), du(17), du(18), du(19), du(20), du(21), du(22), du(23), &
            du(24), du(25), du(26), du(27), du(28), du(29), du(30), du(31), &
            du(32), du(33), du(34), du(35), du(36), du(37), du(38), du(39), &
            du(40), du(41), du(42), jv)
        call dw_bernstein_vjp(u(1), u(2), u(3), u(4), u(5), u(6), u(7), u(8), &
            u(9), u(10), u(11), u(12), u(13), u(14), u(15), u(16), u(17), u(18), &
            u(19), u(20), u(21), u(22), u(23), u(24), u(25), u(26), u(27), u(28), &
            u(29), u(30), u(31), u(32), u(33), u(34), u(35), u(36), u(37), u(38), &
            u(39), u(40), u(41), u(42), 1.0_dp, gr(1), gr(2), gr(3), gr(4), &
            gr(5), gr(6), gr(7), gr(8), gr(9), gr(10), gr(11), gr(12), gr(13), &
            gr(14), gr(15), gr(16), gr(17), gr(18), gr(19), gr(20), gr(21), &
            gr(22), gr(23), gr(24), gr(25), gr(26), gr(27), gr(28), gr(29), &
            gr(30), gr(31), gr(32), gr(33), gr(34), gr(35), gr(36), gr(37), &
            gr(38), gr(39), gr(40), gr(41), gr(42))
        if (abs(dot_product(gr, du) - jv) > 1e-12_dp*(1 + abs(jv))) &
            error stop 'VJP . direction /= JVP'
        h = 1e-5_dp
        up = u + h*du; um = u - h*du
        fp = leaf(up); fm = leaf(um)
        if (abs((fp - fm)/(2*h) - jv) > 1e-8_dp*(1 + abs(jv))) &
            error stop 'JVP differs from central difference'
    end subroutine derivative_pairing

    real(dp) function leaf(u)
        real(dp), intent(in) :: u(42)
        call dw_bernstein(u(1), u(2), u(3), u(4), u(5), u(6), u(7), u(8), u(9), &
            u(10), u(11), u(12), u(13), u(14), u(15), u(16), u(17), u(18), u(19), &
            u(20), u(21), u(22), u(23), u(24), u(25), u(26), u(27), u(28), u(29), &
            u(30), u(31), u(32), u(33), u(34), u(35), u(36), u(37), u(38), u(39), &
            u(40), u(41), u(42), leaf)
    end function leaf

    pure subroutine sample(seed, b, a, c)
        integer, intent(in) :: seed
        real(dp), intent(out) :: b(3), a(3, 3), c(3, 3)
        integer :: i, j
        do i = 1, 3
            b(i) = sin(1.3_dp*seed + i)
            do j = 1, 3
                a(i, j) = cos(0.7_dp*seed*i + 0.3_dp*j)
                c(i, j) = sin(0.5_dp*seed + 1.1_dp*i*j)
            end do
        end do
    end subroutine sample
end program test_ideal_dw
