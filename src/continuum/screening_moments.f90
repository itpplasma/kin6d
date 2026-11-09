! First-order (in grad T) spatial moments of a test-ion population in B = B e_z and
! a background with T_b(x) = T_0 + G.x and heat flux q (benchmarks B05-B07).
! With N = int f d^3x and M = int x f d^3x, to first order in G:
!   d_t N = L N + C^q[N] + G.d_T C[M - rho N]*,  d_t M = L M + v N,
!   L = C_T0 - Omega (v x b).grad_v,  N(0) = f_M(T_0),  M(0) = 0,
! where *rho N is subtracted only in the guiding-centre variant (collisions at
! T_b(X), X = x - rho, rho = b x v/Omega). d<x>/dt = int v N. All sources are l = 1,
! so the state is f = c(v).n on the B04 radial cells; zeta = c_x - i c_y carries the
! Lorentz rotation as +i Omega, c_z does not rotate (gen/gen_thermal_screening.f90).
! Pointwise mathematics comes from the FortSym leaves; this module holds the cell
! bookkeeping and the linear algebra.
module milca_screening_moments
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_quadrature, only: gauss_legendre_ab
    use milca_rosenbluth_reference, only: homma_pair_t, homma_radial_quad, homma_eps0
    use milca_fp_relaxation, only: fp_relax_t, fp_relax_setup
    use kin6d_fp_relaxation_leaves, only: fp_rel_w3, fp_rel_face_lo, fp_rel_face_hi, &
        fp_rel_sink, fp_rel_mom_u, fp_rel_trap_c, fp_rel_bdf2_a0, fp_rel_bdf2_a1, &
        fp_rel_bdf2_a2
    use kin6d_thermal_screening_leaves, only: scr_friction_dt, scr_diff_par_dt, &
        scr_diff_perp_dt, scr_cross_dt, scr_vfm_shell, scr_cq_shell
    implicit none
    private
    public :: scr_t, scr_setup, scr_stencil, scr_steady, scr_rate0, scr_transient

    integer, parameter :: nquad = 4
    real(dp), parameter :: oracle_xmax = 8.0_dp
    ! Slots of homma_radial_quad: h1', g1, g1', g1''.
    integer, parameter :: i_dh1 = 7, i_g1 = 8, i_dg1 = 9, i_d2g1 = 10
    complex(dp), parameter :: iu = (0.0_dp, 1.0_dp)

    ! Cell i spans [vface(i-1), vface(i)]; (K x)_i = kl_i x_(i-1) + kd_i x_i + ku_i x_(i+1)
    ! is the l = 1 collision operator at T_0 integrated over the cell, (Kt x) its T_b
    ! derivative. src_m and src_p are the cell integrals of v^2 (v f_M) and v^2 P,
    ! P per unit q.
    type :: scr_t
        integer :: n = 0
        real(dp) :: h = 0.0_dp
        real(dp), allocatable :: vface(:), vc(:), vol(:), w3(:)
        real(dp), allocatable :: kl(:), kd(:), ku(:), klt(:), kdt(:), kut(:)
        real(dp), allocatable :: src_m(:), src_p(:)
    end type scr_t

contains

    ! Test ions (ma, qa) Maxwellian at the background temperature t0 (J), so that
    ! L f_M = 0; ppu oracle panels per v_t,b for the heat-flux source.
    subroutine scr_setup(s, ncell, ma, mb, qa, qb, nb, t0, lnl, tail, ppu)
        type(scr_t), intent(out) :: s
        integer, intent(in) :: ncell
        real(dp), intent(in) :: ma, mb, qa, qb, nb, t0, lnl, tail, ppu
        type(fp_relax_t) :: f
        type(homma_pair_t) :: pair
        real(dp), allocatable :: rf(:), at(:), et(:), dpt(:), dperpt(:), pft(:), qft(:)
        real(dp) :: xv(nquad), wv(nquad), val, rv(10)
        integer :: i, j, n

        call fp_relax_setup(f, ncell, ma, mb, qa, qb, nb, t0, lnl, homma_eps0, t0, t0, &
            0.0_dp, 1.0_dp, tail)
        n = f%ncell
        s%n = n; s%h = f%h
        allocate (s%vface(0:n), s%vc(n), s%vol(n), s%w3(n), s%src_m(n), s%src_p(n))
        s%vface = f%vface; s%vc = f%vc; s%vol = f%vol
        call fp_rel_w3(f%vface(0:n - 1), f%vface(1:n), s%w3)
        call scr_stencil(f%pf, f%qf, f%h, f%dperp, 2.0_dp, s%kl, s%kd, s%ku)

        ! The face flux is linear in (A_r, E, D_par) and the sink in D_perp, so the
        ! B04 stencil of the T_b-derivative coefficients is d_T K.
        allocate (rf(n - 1), at(n - 1), et(n - 1), dpt(n), dperpt(n), pft(0:n), qft(0:n))
        rf = f%vface(1:n - 1)
        call scr_friction_dt(rf, t0, mb, nb, f%gam, f%mr, at)
        call scr_cross_dt(rf, t0, mb, nb, f%gam, et)
        call scr_diff_par_dt(f%vc, t0, mb, nb, f%gam, dpt)
        call scr_diff_perp_dt(f%vc, t0, mb, nb, f%gam, dperpt)
        pft = 0.0_dp; qft = 0.0_dp
        call fp_rel_face_lo(rf, at, et, dpt(1:n - 1), f%h, pft(1:n - 1))
        call fp_rel_face_hi(rf, at, et, dpt(2:n), f%h, qft(1:n - 1))
        call scr_stencil(pft, qft, f%h, dperpt, 2.0_dp, s%klt, s%kdt, s%kut)

        pair%n_b = nb; pair%t_b = t0; pair%m_b = mb; pair%q_b = qb
        pair%m_a = ma; pair%q_a = qa; pair%lnlambda = lnl
        do i = 1, n
            call gauss_legendre_ab(nquad, f%vface(i - 1), f%vface(i), xv, wv)
            s%src_m(i) = 0.0_dp
            s%src_p(i) = 0.0_dp
            do j = 1, nquad
                call scr_vfm_shell(xv(j), ma, t0, val)
                s%src_m(i) = s%src_m(i) + wv(j)*val
                call homma_radial_quad(pair, xv(j), ppu, oracle_xmax, rv)
                call scr_cq_shell(xv(j), ma, t0, t0, mb, f%gam, f%mr, rv(i_dh1), rv(i_g1), &
                    rv(i_dg1), rv(i_d2g1), val)
                s%src_p(i) = s%src_p(i) + wv(j)*val
            end do
        end do
    end subroutine scr_setup

    ! B04 per-l stencil from the face weights pf (own cell) and qf (upper cell) and
    ! the angular sink l(l+1) D_perp; faces 0 and n carry zero flux.
    subroutine scr_stencil(pf, qf, h, dperp, ll, kl, kd, ku)
        real(dp), intent(in) :: pf(0:), qf(0:), h, dperp(:), ll
        real(dp), allocatable, intent(out) :: kl(:), kd(:), ku(:)
        real(dp) :: sink(size(dperp))
        integer :: i, n

        n = size(dperp)
        allocate (kl(n), kd(n), ku(n))
        call fp_rel_sink(h, dperp, ll, sink)
        do i = 1, n
            kl(i) = -pf(i - 1)
            kd(i) = pf(i) - qf(i - 1) - sink(i)
            ku(i) = qf(i)
        end do
    end subroutine scr_stencil

    ! Steady first-order mean velocity (m/s) for grad T_b = grad_t (J/m) and heat
    ! flux q (W/m^2): L m = -v f_M, L n1 = -(C^q[f_M] + G.d_T C[m - rho f_M]*).
    subroutine scr_steady(s, omega, grad_t, q, gc, drift)
        type(scr_t), intent(in) :: s
        real(dp), intent(in) :: omega, grad_t(3), q(3)
        logical, intent(in) :: gc
        real(dp), intent(out) :: drift(3)
        complex(dp) :: zm(s%n), cm(s%n), zs(s%n), cs(s%n), zn(s%n), cn(s%n)

        call check_variant(omega, gc)
        call solve(s, 0.0_dp, -1.0_dp, omega, cmplx(-s%src_m, 0.0_dp, dp), zm)
        call solve(s, 0.0_dp, -1.0_dp, 0.0_dp, cmplx(-s%src_m, 0.0_dp, dp), cm)
        call n_source(s, omega, grad_t, q, gc, zm, cm, zs, cs)
        call solve(s, 0.0_dp, -1.0_dp, omega, -zs, zn)
        call solve(s, 0.0_dp, -1.0_dp, 0.0_dp, -cs, cn)
        call mean_velocity(s, zn, cn, drift)
    end subroutine scr_steady

    ! d<v>/dt at t = 0 (m/s^2) of the semi-discrete system: V dn1/dt = S(0).
    subroutine scr_rate0(s, omega, grad_t, q, gc, rate)
        type(scr_t), intent(in) :: s
        real(dp), intent(in) :: omega, grad_t(3), q(3)
        logical, intent(in) :: gc
        real(dp), intent(out) :: rate(3)
        complex(dp) :: zero_m(s%n), zs(s%n), cs(s%n)

        call check_variant(omega, gc)
        zero_m = (0.0_dp, 0.0_dp)
        call n_source(s, omega, grad_t, q, gc, zero_m, zero_m, zs, cs)
        call mean_velocity(s, zs/s%vol, cs/s%vol, rate)
    end subroutine scr_rate0

    ! TR-BDF2 (gamma = 2 - sqrt(2), B04 weights) from N = f_M, M = 0 over nsteps
    ! equal steps to tend (s); vel(:, k) = <v>(k tend/nsteps) in m/s. The M block is
    ! advanced first because the N source depends on M at the stage times.
    subroutine scr_transient(s, omega, grad_t, q, gc, tend, nsteps, vel)
        type(scr_t), intent(in) :: s
        real(dp), intent(in) :: omega, grad_t(3), q(3), tend
        logical, intent(in) :: gc
        integer, intent(in) :: nsteps
        real(dp), intent(out) :: vel(3, 0:nsteps)
        complex(dp), dimension(s%n) :: zm, cm, zn, cn, zm1, cm1, zn1, cn1, zm2, cm2
        complex(dp), dimension(s%n) :: zs0, cs0, zs1, cs1, zs2, cs2, msrc
        real(dp) :: dt, c, a0, a1, a2
        integer :: k

        call check_variant(omega, gc)
        if (nsteps < 1 .or. .not. tend > 0.0_dp) error stop 'scr_transient: invalid time grid'
        dt = tend/real(nsteps, dp)
        call fp_rel_trap_c(dt, c)
        call fp_rel_bdf2_a0(dt, a0)
        call fp_rel_bdf2_a1(dt, a1)
        call fp_rel_bdf2_a2(dt, a2)
        msrc = cmplx(s%src_m, 0.0_dp, dp)
        zm = 0.0_dp; cm = 0.0_dp; zn = 0.0_dp; cn = 0.0_dp
        vel(:, 0) = 0.0_dp
        call n_source(s, omega, grad_t, q, gc, zm, cm, zs0, cs0)
        do k = 1, nsteps
            call stage_tr(s, omega, c, zm, msrc, msrc, zm1)
            call stage_tr(s, 0.0_dp, c, cm, msrc, msrc, cm1)
            call stage_bdf(s, omega, a0, a1, a2, zm, zm1, msrc, zm2)
            call stage_bdf(s, 0.0_dp, a0, a1, a2, cm, cm1, msrc, cm2)
            call n_source(s, omega, grad_t, q, gc, zm1, cm1, zs1, cs1)
            call n_source(s, omega, grad_t, q, gc, zm2, cm2, zs2, cs2)
            call stage_tr(s, omega, c, zn, zs0, zs1, zn1)
            call stage_tr(s, 0.0_dp, c, cn, cs0, cs1, cn1)
            call stage_bdf(s, omega, a0, a1, a2, zn, zn1, zs2, zn)
            call stage_bdf(s, 0.0_dp, a0, a1, a2, cn, cn1, cs2, cn)
            zm = zm2; cm = cm2; zs0 = zs2; cs0 = cs2
            call mean_velocity(s, zn, cn, vel(:, k))
        end do
    end subroutine scr_transient

    subroutine check_variant(omega, gc)
        real(dp), intent(in) :: omega
        logical, intent(in) :: gc
        if (gc .and. .not. abs(omega) > 0.0_dp) error stop 'guiding-centre variant needs B /= 0'
    end subroutine check_variant

    ! Cell-integrated N source for the M state (zm: zeta, cm: c_z):
    ! zeta: (q_x - i q_y) P + Kt (G_x - i G_y)(zm - [gc] i v f_M/Omega), c_z: q_z P + Kt G_z cm.
    subroutine n_source(s, omega, grad_t, q, gc, zm, cm, zs, cs)
        type(scr_t), intent(in) :: s
        real(dp), intent(in) :: omega, grad_t(3), q(3)
        logical, intent(in) :: gc
        complex(dp), intent(in) :: zm(:), cm(:)
        complex(dp), intent(out) :: zs(:), cs(:)
        complex(dp) :: zg(s%n)

        zg = cmplx(grad_t(1), -grad_t(2), dp)*zm
        if (gc) zg = zg - cmplx(grad_t(1), -grad_t(2), dp)*iu*s%src_m/(omega*s%vol)
        call apply(s%klt, s%kdt, s%kut, zg, zs)
        zs = zs + cmplx(q(1), -q(2), dp)*s%src_p
        call apply(s%klt, s%kdt, s%kut, grad_t(3)*cm, cs)
        cs = cs + q(3)*s%src_p
    end subroutine n_source

    ! <v> from cell values: (4 pi/3) int v^3 c dv per component, c_y = -Im zeta.
    subroutine mean_velocity(s, zeta, cz, vel)
        type(scr_t), intent(in) :: s
        complex(dp), intent(in) :: zeta(:), cz(:)
        real(dp), intent(out) :: vel(3)
        real(dp) :: m(s%n)

        call fp_rel_mom_u(s%w3, real(zeta, dp), m)
        vel(1) = sum(m)
        call fp_rel_mom_u(s%w3, -aimag(zeta), m)
        vel(2) = sum(m)
        call fp_rel_mom_u(s%w3, real(cz, dp), m)
        vel(3) = sum(m)
    end subroutine mean_velocity

    ! Trapezoid over gamma dt: (V - c A) y1 = V y0 + c (A y0 + s0 + s1), A = K + i w V.
    subroutine stage_tr(s, w, c, y0, s0, s1, y1)
        type(scr_t), intent(in) :: s
        real(dp), intent(in) :: w, c
        complex(dp), intent(in) :: y0(:), s0(:), s1(:)
        complex(dp), intent(out) :: y1(:)
        complex(dp) :: ay(s%n)

        call apply(s%kl, s%kd, s%ku, y0, ay)
        ay = ay + iu*w*s%vol*y0
        call solve(s, 1.0_dp, c, w, s%vol*y0 + c*(ay + s0 + s1), y1)
    end subroutine stage_tr

    ! BDF2: (a0 V - A) y2 = -V (a1 y1 + a2 y0) + s2.
    subroutine stage_bdf(s, w, a0, a1, a2, y0, y1, s2, y2)
        type(scr_t), intent(in) :: s
        real(dp), intent(in) :: w, a0, a1, a2
        complex(dp), intent(in) :: y0(:), y1(:), s2(:)
        complex(dp), intent(out) :: y2(:)

        call solve(s, a0, 1.0_dp, w, -s%vol*(a1*y1 + a2*y0) + s2, y2)
    end subroutine stage_bdf

    subroutine apply(kl, kd, ku, x, y)
        real(dp), intent(in) :: kl(:), kd(:), ku(:)
        complex(dp), intent(in) :: x(:)
        complex(dp), intent(out) :: y(:)
        integer :: n

        n = size(x)
        y = kd*x
        y(2:n) = y(2:n) + kl(2:n)*x(1:n - 1)
        y(1:n - 1) = y(1:n - 1) + ku(1:n - 1)*x(2:n)
    end subroutine apply

    ! Solves (alpha V - beta (K + i w V)) x = rhs by the Thomas algorithm.
    subroutine solve(s, alpha, beta, w, rhs, x)
        type(scr_t), intent(in) :: s
        real(dp), intent(in) :: alpha, beta, w
        complex(dp), intent(in) :: rhs(:)
        complex(dp), intent(out) :: x(:)
        complex(dp) :: di(s%n), cp(s%n), dq(s%n), den
        integer :: i, n

        n = s%n
        di = alpha*s%vol - beta*(s%kd + iu*w*s%vol)
        cp(1) = -beta*s%ku(1)/di(1)
        dq(1) = rhs(1)/di(1)
        do i = 2, n
            den = di(i) + beta*s%kl(i)*cp(i - 1)
            cp(i) = -beta*s%ku(i)/den
            dq(i) = (rhs(i) + beta*s%kl(i)*dq(i - 1))/den
        end do
        x(n) = dq(n)
        do i = n - 1, 1, -1
            x(i) = dq(i) - cp(i)*x(i + 1)
        end do
    end subroutine solve

end module milca_screening_moments
