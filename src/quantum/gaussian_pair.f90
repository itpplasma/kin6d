! Closed-form matrix elements between complex isotropic 3D Gaussians
!   g_k(r) = exp(-alpha_k |r - q_k|^2 + i kappa_k.(r - q_k)),
! with real centre q_k and wave vector kappa_k = p_k/hbar in the x-z plane
! (y components zero) and Re alpha_k > 0. With the holomorphic coefficients
! a = alpha, b = 2 alpha q + i kappa, c = -alpha q.q - i kappa.q the product
! conj(g_k) g_l = exp(-A r^2 + B.r + C), A = conj(a_k) + a_l, etc., so
!   S = <g_k|g_l> = (pi/A)^(3/2) exp(T + C),  T = B.B/(4A)  (bilinear dot),
!   <g_k|1/r|g_l> = (2 pi/A) exp(C) E_0(T),   E_n(T) = int_0^1 u^2n e^{T(1-u^2)} du,
! the complex Boys function scaled by e^T. Polynomial moments are taken
! about the complex product centre c = B/(2A). Every element is exact; the
! only approximation is the evaluation of E_n (Gauss-Legendre for |T| <= 36,
! endpoint asymptotic series beyond, both verified against binary128
! quadrature in test/test_packet_models.f90).
!
! Elementary functions per Gaussian: e1 = g, e2 = u_x g, e3 = u_z g,
! e4 = |u|^2 g with u = r - q. These complex closed forms are written by hand
! because FortSym does not emit complex arithmetic (docs/packet-models.md).
module kin6d_gaussian_pair
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_quadrature, only: gauss_legendre
    implicit none
    private
    public :: gaussian_t, pair_rule_t, pair_rule_init, boys_scaled
    public :: pair_elements, gaussian_value, gaussian_norm2

    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp), parameter :: t_switch = 36

    type :: gaussian_t
        real(dp) :: q(2) = 0, kappa(2) = 0   ! (x, z) components
        complex(dp) :: alpha = (1, 0)
        complex(dp) :: d = (1, 0)            ! linear amplitude (not used by pair)
    end type gaussian_t

    ! Gauss-Legendre rule on [0,1] for the Boys integrals (explicit state).
    type :: pair_rule_t
        real(dp), allocatable :: u(:), w(:)
    end type pair_rule_t
contains
    subroutine pair_rule_init(rule, n)
        type(pair_rule_t), intent(out) :: rule
        integer, intent(in) :: n
        real(dp), allocatable :: x(:), w(:)
        allocate (x(n), w(n))
        call gauss_legendre(n, x, w)
        rule%u = (x + 1)/2; rule%w = w/2
    end subroutine pair_rule_init

    ! E_n(T) = exp(tau) f(n), n = 0..2, with tau = T for Re T > 0, else 0,
    ! so that |f| stays bounded (|f_n| <= 1/(2n+1) for the returned scaling).
    pure subroutine boys_scaled(rule, t, tau, f)
        type(pair_rule_t), intent(in) :: rule
        complex(dp), intent(in) :: t
        complex(dp), intent(out) :: tau, f(0:2)
        complex(dp) :: z, ser, term, e
        integer :: j, n
        if (real(t) > 0) then
            tau = t
        else
            tau = 0
        end if
        if (abs(t) <= t_switch) then
            f = 0
            do j = 1, size(rule%u)
                e = rule%w(j)*exp(t*(1 - rule%u(j)**2) - tau)
                do n = 0, 2
                    f(n) = f(n) + e*rule%u(j)**(2*n)
                end do
            end do
            return
        end if
        ! Large |T|: F_0 = sqrt(pi)/(2 z) - e^{-T}/(2T) sum_j (-1)^j (2j-1)!!/(2T)^j
        ! (erfc asymptotics, principal z = sqrt(T), valid for |arg z| <= pi/2),
        ! then E_{n+1} = ((2n+1) E_n - 1)/(2T), stable for |T| > 36.
        z = sqrt(t)
        ser = 1; term = 1
        do j = 1, 30
            term = -term*real(2*j - 1, dp)/(2*t)
            ser = ser + term
        end do
        ser = ser/(2*t)
        ! Scaled by exp(-tau): E_0 e^{-tau} = e^{T - tau} sqrt(pi)/(2z) - e^{-tau} ser.
        f(0) = exp(t - tau)*sqrt(pi)/(2*z) - exp(-tau)*ser
        do n = 0, 1
            f(n + 1) = (real(2*n + 1, dp)*f(n) - exp(-tau))/(2*t)
        end do
    end subroutine boys_scaled

    ! Pair elements between bra k and ket l:
    !   s      = <g_k|g_l>
    !   sig    = <e^k_i|e^l_j>/s           (4x4 normalized moments)
    !   hel    = <e^k_i|H|g_l>             (4), H = -hbar^2/(2 mass) lap + g/r + kap r^2/2
    ! A non-positive mass drops the kinetic term (used for pure overlaps).
    pure subroutine pair_elements(rule, gk, gl, hbar, mass, g, kap, s, sig, hel)
        type(pair_rule_t), intent(in) :: rule
        type(gaussian_t), intent(in) :: gk, gl
        real(dp), intent(in) :: hbar, mass, g, kap
        complex(dp), intent(out) :: s, sig(4, 4), hel(4)
        complex(dp) :: a, bb(2), cc, t, cen(2), dk(2), dl(2), h, dk2, dl2, tau, f(0:2)
        complex(dp) :: fac, w0, w1(2), wr2, tcoef(4), v(4), dkl
        real(dp) :: qk2
        a = conjg(gk%alpha) + gl%alpha
        bb = conjg(2*gk%alpha*gk%q + cmplx(0, 1, dp)*gk%kappa) + &
             2*gl%alpha*gl%q + cmplx(0, 1, dp)*gl%kappa
        cc = conjg(-gk%alpha*sum(gk%q**2) - cmplx(0, 1, dp)*sum(gk%kappa*gk%q)) &
             - gl%alpha*sum(gl%q**2) - cmplx(0, 1, dp)*sum(gl%kappa*gl%q)
        t = sum(bb*bb)/(4*a)
        s = exp(1.5_dp*log(pi/a) + t + cc)
        cen = bb/(2*a)
        dk = cen - gk%q; dl = cen - gl%q
        h = 1/(2*a)
        dk2 = sum(dk*dk); dl2 = sum(dl*dl); dkl = sum(dk*dl)
        sig(1, 1) = 1
        sig(1, 2:3) = dl; sig(2:3, 1) = dk
        sig(2, 2) = h + dk(1)*dl(1); sig(2, 3) = dk(1)*dl(2)
        sig(3, 2) = dk(2)*dl(1); sig(3, 3) = h + dk(2)*dl(2)
        sig(1, 4) = 3*h + dl2; sig(4, 1) = 3*h + dk2
        sig(2:3, 4) = 2*h*dl + dk*(3*h + dl2)
        sig(4, 2:3) = 2*h*dk + dl*(3*h + dk2)
        sig(4, 4) = 15*h**2 + 3*h*(dk2 + dl2) + 4*h*dkl + dk2*dl2
        ! Kinetic and harmonic actions on g_l in the elementary basis of l.
        tcoef = 0
        if (mass > 0) then
            tcoef = -(hbar**2/(2*mass))*[-6*gl%alpha - sum(gl%kappa**2), &
                    -4*cmplx(0, 1, dp)*gl%alpha*gl%kappa(1), &
                    -4*cmplx(0, 1, dp)*gl%alpha*gl%kappa(2), 4*gl%alpha**2]
        end if
        tcoef = tcoef + (kap/2)*[cmplx(sum(gl%q**2), 0, dp), cmplx(2*gl%q(1), 0, dp), &
                                 cmplx(2*gl%q(2), 0, dp), (1.0_dp, 0.0_dp)]
        hel = s*matmul(sig, tcoef)
        if (abs(g) <= 0.0_dp) return
        ! Coulomb moments W0 = int e^{..}/r, W_a = int r_a e^{..}/r, W_r2.
        call boys_scaled(rule, t, tau, f)
        fac = (2*pi/a)*exp(tau + cc)
        w0 = fac*f(0)
        w1 = fac*(f(0) - f(1))*bb/(2*a)
        wr2 = fac*(f(0) + t*(f(0) - f(1)))/a
        qk2 = sum(gk%q**2)
        v(1) = w0
        v(2:3) = w1 - gk%q*w0
        v(4) = wr2 - 2*sum(gk%q*w1) + qk2*w0
        hel = hel + g*v
    end subroutine pair_elements

    ! g(r) at Cartesian points (y enters through |u|^2 only).
    elemental complex(dp) function gaussian_value(gk, x, y, z) result(v)
        type(gaussian_t), intent(in) :: gk
        real(dp), intent(in) :: x, y, z
        real(dp) :: ux, uz
        ux = x - gk%q(1); uz = z - gk%q(2)
        v = exp(-gk%alpha*(ux**2 + y**2 + uz**2) + &
                cmplx(0, gk%kappa(1)*ux + gk%kappa(2)*uz, dp))
    end function gaussian_value

    ! ||g||^2 = (pi/(2 Re alpha))^(3/2).
    elemental real(dp) function gaussian_norm2(gk) result(n2)
        type(gaussian_t), intent(in) :: gk
        n2 = (pi/(2*real(gk%alpha)))**1.5_dp
    end function gaussian_norm2
end module kin6d_gaussian_pair
