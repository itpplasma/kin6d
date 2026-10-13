module kin6d_mhd_layer_n
    !! L5: resistive closure of the per-n response at rational surfaces by
    !! matching to the resistive-inertial(-viscous) inner layer in a
    !! rotating plasma (the MARS-F single-fluid layer physics, Liu et al.,
    !! Phys. Plasmas 7 (2000) 3681, with the layer reduced to its matched
    !! asymptotic form).  The outer region is the ideal static response of
    !! kin6d_mhd_ideal_n (D5: the shared per-n Ampere operator and the
    !! MEPHIT-form volume current); this module only supplies the layer
    !! closure, a separate current model on the layer edge sheets of
    !! kin6d_mhd_closure_n.
    !!
    !! Layer equations.  Near q = m/n, x the normal distance, F(x) = k.B =
    !! F' x, a perturbation exp(i(n phi - omega t)) in a plasma rotating
    !! toroidally with Omega has the Doppler rate gamma = -i (omega - n
    !! Omega) (gamma = i n Omega for the static coil field).  With b the
    !! resonant flux and v the radial velocity stream function, the slab
    !! layer (Furth, Killeen and Rosenbluth 1963; Fitzpatrick, Phys. Plasmas
    !! 5 (1998) 3325) reads
    !!     gamma b = i F' x v + D b'',                  D = eta/mu0,
    !!     gamma v'' - nu v'''' = i (F' x/(mu0 rho)) b''.
    !! Scaling x = L s, L^4 = |gamma| rho D mu0/F'^2, v = (|gamma|/(F' L))
    !! vt gives, with g = gamma/|gamma|, eps = D/(L^2 |gamma|) (= Q^(-3/2)
    !! in Fitzpatrick's notation) and P = nu/(|gamma| L^2),
    !!     eps b'' = g b - i s vt,
    !!     g vt'' - P vt'''' - s^2 vt = i g s b,
    !! b even, vt odd.  Far from the layer b = b0 (1 + Dhat |s|/2) and the
    !! layer response is Delta = [b']/b0 = Dhat/L (per metre); b(0)/b0 is
    !! the flux at the rational relative to the outer (matching) value, 1
    !! in the constant-psi regime and small in the inertial one.  The
    !! constant-psi resistive-inertial limit (eps -> infinity) is
    !!     Dhat = 2 pi Gamma(3/4)/Gamma(1/4) g^(5/4)/eps = 2.1236 g^(5/4)/eps.
    !! The general case (finite eps: non-constant psi; P > 0: viscous) is
    !! solved numerically (finite differences on a sinh-graded grid,
    !! Richardson extrapolated).
    !!
    !! Toroidal coefficients.  For a rational at psin_s with the flux label
    !! X = psin - psin_s: F' x = n |dq/dpsin| <B.grad theta> X with the
    !! chart-independent volume average <B.grad theta> = 4 pi^2 psi_o/V'
    !! (psi_o = |psi_edge - psi_axis| per radian, V' = dV/dpsin), and X =
    !! G x with G = <|grad psin|^2>^(1/2); so F' = n |q'| <B.grad theta> G
    !! and Delta_psin = Delta/G.  This is the cylinder-equivalent layer
    !! (exact in the large-aspect-ratio limit); curvature (Glasser) terms
    !! and the toroidal inertia enhancement are not included (rho may be
    !! scaled by the caller, e.g. 1 + 2 q^2).
    !!
    !! Closure.  On the two layer-edge sheets of an L6 layer
    !! (kin6d_mhd_closure_n) the resistive constraints are
    !!     a1(in) - a1(out) = 0        (continuous reconnected flux),
    !!     J - Delta_psin (a1(in) + a1(out))/2 = 0,
    !! with a1 the large-solution coefficients of the outer fits and J =
    !! junit (lambda_in + lambda_out) the jump of d Phi/d psin across the
    !! layer (net layer current).  As Delta -> infinity (eta -> 0 at fixed
    !! rotation) both tend to the ideal L6 constraints a1 = 0.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_vacuum_surfaces, only: psi_grid_t, surface_t, surface_build, &
        psi_eval
    use kin6d_mhd_ampere_n, only: mu0
    use kin6d_mhd_response_n, only: sheet_closure_t
    use kin6d_mhd_closure_n, only: rational_layer_t
    implicit none
    private

    real(dp), parameter :: pi = 3.14159265358979323846_dp
    !> 2 pi Gamma(3/4)/Gamma(1/4)
    real(dp), parameter, public :: ri_constant = 2.1236482381522474_dp

    !> Physical layer parameters (SI).
    type, public :: layer_params_t
        complex(dp) :: gamma = (0.0_dp, 0.0_dp)  !! -i (omega - n Omega), 1/s
        real(dp) :: diff = 0.0_dp                !! eta/mu0, m^2/s
        real(dp) :: rho = 0.0_dp                 !! mass density, kg/m^3
        real(dp) :: nu = 0.0_dp                  !! viscosity, m^2/s
        real(dp) :: fprime = 0.0_dp              !! d(k.B)/dx, T/m^2
    end type layer_params_t

    public :: layer_scaled_delta
    public :: layer_ri_constant_psi
    public :: layer_scales
    public :: layer_delta
    public :: layer_metrics
    public :: resistive_closure

contains

    pure function layer_ri_constant_psi(g, eps) result(d)
        !! Constant-psi resistive-inertial Dhat = 2.1236 g^(5/4)/eps
        !! (principal branch, |arg g| <= pi/2).
        complex(dp), intent(in) :: g
        real(dp), intent(in) :: eps
        complex(dp) :: d
        d = ri_constant*exp(1.25_dp*log(g))/eps
    end function layer_ri_constant_psi

    subroutine layer_scales(par, l, eps, p, g)
        !! Scale L (m), eps, P and the unit Doppler factor g of the layer.
        type(layer_params_t), intent(in) :: par
        real(dp), intent(out) :: l, eps, p
        complex(dp), intent(out) :: g
        real(dp) :: ag
        ag = abs(par%gamma)
        l = (ag*par%rho*par%diff*mu0/par%fprime**2)**0.25_dp
        eps = par%diff/(l*l*ag)
        p = par%nu/(ag*l*l)
        g = par%gamma/ag
    end subroutine layer_scales

    function layer_delta(par, status, center) result(d)
        !! Layer response Delta = [b']/b0 per metre; center = b(0)/b0.
        type(layer_params_t), intent(in) :: par
        integer, intent(out) :: status
        complex(dp), intent(out), optional :: center
        complex(dp) :: d
        real(dp) :: l, eps, p
        complex(dp) :: g, c
        call layer_scales(par, l, eps, p, g)
        d = layer_scaled_delta(g, eps, p, status, center=c)/l
        if (present(center)) center = c
    end function layer_delta

    function layer_scaled_delta(g, eps, p, status, npts, center) result(dhat)
        !! Dhat of the scaled layer equations (finite differences, two grids,
        !! Richardson extrapolation of the second-order error); center =
        !! b(0)/b0, the flux at the rational over the extrapolated outer flux.
        complex(dp), intent(in) :: g
        real(dp), intent(in) :: eps, p
        integer, intent(out) :: status
        integer, intent(in), optional :: npts
        complex(dp), intent(out), optional :: center
        complex(dp) :: dhat
        complex(dp) :: d1, d2, c1, c2
        integer :: nn

        nn = 2000
        if (present(npts)) nn = npts
        d1 = solve_grid(nn, status, c1)
        if (status /= 0) return
        d2 = solve_grid(2*nn, status, c2)
        if (status /= 0) return
        dhat = (4*d2 - d1)/3
        if (present(center)) center = (4*c2 - c1)/3
    contains
        function solve_grid(n, st, cen) result(dh)
            integer, intent(in) :: n
            integer, intent(out) :: st
            complex(dp), intent(out) :: cen
            complex(dp) :: dh
            real(dp), allocatable :: s(:)
            complex(dp), allocatable :: ab(:, :), rhs(:, :), fi(:)
            integer, allocatable :: ipiv(:)
            integer :: nv, nu, kl, ku, ldab, j, row
            real(dp) :: smax, h0, c, hm, hp, a, cm, c0, cp
            complex(dp) :: integral, bprime, b0

            ! grid: s = smax sinh(c t)/sinh(c), t uniform, spacing at 0 ~ h0
            smax = 40.0_dp*max(1.0_dp, 1.0_dp/sqrt(eps), p**0.25_dp)
            h0 = 0.02_dp*min(1.0_dp, sqrt(eps), max(p**0.25_dp, 1.0e-3_dp))
            if (p <= 0.0_dp) h0 = 0.02_dp*min(1.0_dp, sqrt(eps))
            h0 = h0*2000.0_dp/n
            c = grading(smax, h0, n)
            allocate (s(0:n))
            do j = 0, n
                if (c > 1.0e-8_dp) then
                    s(j) = smax*sinh(c*j/real(n, dp))/sinh(c)
                else
                    s(j) = smax*j/real(n, dp)
                end if
            end do
            nv = 2
            if (p > 0.0_dp) nv = 3
            nu = nv*(n + 1)
            kl = 2*nv - 1
            ku = 2*nv - 1
            ldab = 2*kl + ku + 1
            allocate (ab(ldab, nu), rhs(nu, 1), ipiv(nu))
            ab = (0.0_dp, 0.0_dp)
            rhs = (0.0_dp, 0.0_dp)
            do j = 0, n
                if (j == 0) then
                    ! b'(0) = 0 (even): symmetric stencil; v = 0, w = 0 (odd)
                    hp = s(1)
                    a = 2.0_dp/(hp*hp)
                    row = idx(0, 1)
                    call put_band(ab, kl, ku, row, idx(0, 1), -eps*a - g)
                    call put_band(ab, kl, ku, row, idx(1, 1), cmplx(eps*a, 0.0_dp, dp))
                    call put_band(ab, kl, ku, idx(0, 2), idx(0, 2), (1.0_dp, 0.0_dp))
                    if (nv == 3) call put_band(ab, kl, ku, idx(0, 3), idx(0, 3), (1.0_dp, 0.0_dp))
                    cycle
                end if
                if (j == n) then
                    ! b = 1, ideal outer relation vt = -i g b/s
                    call put_band(ab, kl, ku, idx(n, 1), idx(n, 1), (1.0_dp, 0.0_dp))
                    rhs(idx(n, 1), 1) = (1.0_dp, 0.0_dp)
                    call put_band(ab, kl, ku, idx(n, 2), idx(n, 2), (1.0_dp, 0.0_dp))
                    rhs(idx(n, 2), 1) = cmplx(0.0_dp, -1.0_dp, dp)*g/s(n)
                    if (nv == 3) then
                        call put_band(ab, kl, ku, idx(n, 3), idx(n, 3), (1.0_dp, 0.0_dp))
                        rhs(idx(n, 3), 1) = cmplx(0.0_dp, -2.0_dp, dp)*g/s(n)**3
                    end if
                    cycle
                end if
                hm = s(j) - s(j - 1)
                hp = s(j + 1) - s(j)
                cm = 2.0_dp/(hm*(hm + hp))
                cp = 2.0_dp/(hp*(hm + hp))
                c0 = -cm - cp
                ! eps b'' - g b + i s vt = 0
                row = idx(j, 1)
                call put_band(ab, kl, ku, row, idx(j - 1, 1), cmplx(eps*cm, 0.0_dp, dp))
                call put_band(ab, kl, ku, row, idx(j, 1), eps*c0 - g)
                call put_band(ab, kl, ku, row, idx(j + 1, 1), cmplx(eps*cp, 0.0_dp, dp))
                call put_band(ab, kl, ku, row, idx(j, 2), cmplx(0.0_dp, s(j), dp))
                if (nv == 2) then
                    ! g vt'' - s^2 vt - i g s b = 0
                    row = idx(j, 2)
                    call put_band(ab, kl, ku, row, idx(j - 1, 2), g*cm)
                    call put_band(ab, kl, ku, row, idx(j, 2), g*c0 - s(j)**2)
                    call put_band(ab, kl, ku, row, idx(j + 1, 2), g*cp)
                    call put_band(ab, kl, ku, row, idx(j, 1), cmplx(0.0_dp, -s(j), dp)*g)
                else
                    ! w = vt'':  vt'' - w = 0;  g w - P w'' - s^2 vt - i g s b = 0
                    row = idx(j, 2)
                    call put_band(ab, kl, ku, row, idx(j - 1, 2), cmplx(cm, 0.0_dp, dp))
                    call put_band(ab, kl, ku, row, idx(j, 2), cmplx(c0, 0.0_dp, dp))
                    call put_band(ab, kl, ku, row, idx(j + 1, 2), cmplx(cp, 0.0_dp, dp))
                    call put_band(ab, kl, ku, row, idx(j, 3), (-1.0_dp, 0.0_dp))
                    row = idx(j, 3)
                    call put_band(ab, kl, ku, row, idx(j - 1, 3), cmplx(-p*cm, 0.0_dp, dp))
                    call put_band(ab, kl, ku, row, idx(j, 3), g - p*c0)
                    call put_band(ab, kl, ku, row, idx(j + 1, 3), cmplx(-p*cp, 0.0_dp, dp))
                    call put_band(ab, kl, ku, row, idx(j, 2), cmplx(-s(j)**2, 0.0_dp, dp))
                    call put_band(ab, kl, ku, row, idx(j, 1), cmplx(0.0_dp, -s(j), dp)*g)
                end if
            end do
            call zgbsv(nu, kl, ku, 1, ab, ldab, ipiv, rhs, nu, st)
            if (st /= 0) return
            ! b'(smax) = (1/eps) int_0^smax (g b - i s vt) ds (trapezoid)
            allocate (fi(0:n))
            do j = 0, n
                fi(j) = g*rhs(idx(j, 1), 1) - cmplx(0.0_dp, s(j), dp)*rhs(idx(j, 2), 1)
            end do
            integral = (0.0_dp, 0.0_dp)
            do j = 1, n
                integral = integral + 0.5_dp*(s(j) - s(j - 1))*(fi(j) + fi(j - 1))
            end do
            bprime = integral/eps
            b0 = rhs(idx(n, 1), 1) - s(n)*bprime
            dh = 2.0_dp*bprime/b0
            cen = rhs(idx(0, 1), 1)/b0
            st = 0
        end function solve_grid

        pure integer function idx(j, iv)
            integer, intent(in) :: j, iv
            integer :: nvv
            nvv = 2
            if (p > 0.0_dp) nvv = 3
            idx = nvv*j + iv
        end function idx
    end function layer_scaled_delta

    subroutine put_band(ab, kl, ku, i, j, v)
        complex(dp), intent(inout) :: ab(:, :)
        integer, intent(in) :: kl, ku, i, j
        complex(dp), intent(in) :: v
        ab(kl + ku + 1 + i - j, j) = ab(kl + ku + 1 + i - j, j) + v
    end subroutine put_band

    real(dp) function grading(smax, h0, n) result(c)
        !! c with smax c/(n sinh c) = h0 (first spacing), c >= 0.
        real(dp), intent(in) :: smax, h0
        integer, intent(in) :: n
        real(dp) :: lo, hi, mid, target
        integer :: it
        target = n*h0/smax
        if (target >= 1.0_dp) then
            c = 0.0_dp
            return
        end if
        lo = 1.0e-8_dp
        hi = 200.0_dp
        do it = 1, 200
            mid = 0.5_dp*(lo + hi)
            if (mid/sinh(mid) > target) then
                lo = mid
            else
                hi = mid
            end if
        end do
        c = 0.5_dp*(lo + hi)
    end function grading

    subroutine layer_metrics(g, psin, nth, vprime, bth, gpsi, status)
        !! V' = dV/dpsin, <B.grad theta> = 4 pi^2 psi_o/V' and G = <|grad
        !! psin|^2>^(1/2) (volume average) of the surface psin.
        type(psi_grid_t), intent(in) :: g
        real(dp), intent(in) :: psin
        integer, intent(in) :: nth
        real(dp), intent(out) :: vprime, bth, gpsi
        integer, intent(out) :: status
        type(surface_t) :: s
        real(dp) :: ps, pr, pz, gs, w, wsum, gsum, psio
        integer :: j

        psio = abs(g%psi_edge - g%psi_axis)
        call surface_build(g, psin, nth, s, status)
        if (status /= 0) return
        wsum = 0.0_dp
        gsum = 0.0_dp
        do j = 1, nth
            call psi_eval(g, s%r(j), s%z(j), ps, pr, pz)
            gs = sqrt(pr**2 + pz**2)/psio
            w = s%r(j)*s%dl(j)/gs
            wsum = wsum + w
            gsum = gsum + w*gs*gs
        end do
        vprime = 2*pi*wsum*2*pi/nth
        bth = 4*pi*pi*psio/vprime
        gpsi = sqrt(gsum/wsum)
    end subroutine layer_metrics

    subroutine resistive_closure(lays, dlayer, junit, cl)
        !! Replace the constraints of the L6 layer closure cl (built by
        !! layer_closure for lays) by the resistive matching conditions with
        !! the layer responses dlayer (per psin) and the jumps of d Phi/d
        !! psin per unit sheet junit.
        type(rational_layer_t), intent(in) :: lays(:)
        complex(dp), intent(in) :: dlayer(:)
        complex(dp), intent(in) :: junit(:)
        type(sheet_closure_t), intent(inout) :: cl
        integer :: l, side, k, kf, r1, r2, j
        complex(dp) :: sc

        if (allocated(cl%wl)) deallocate (cl%wl)
        allocate (cl%wl(cl%ns, cl%ns))
        cl%wl = (0.0_dp, 0.0_dp)
        cl%w = (0.0_dp, 0.0_dp)
        cl%target = (0.0_dp, 0.0_dp)
        kf = 0
        do l = 1, size(lays)
            r1 = 2*(l - 1) + 1
            r2 = 2*(l - 1) + 2
            ! scale the current row so that it tends to the ideal row
            sc = (1.0_dp, 0.0_dp)
            if (abs(dlayer(l)) > 1.0_dp) sc = 1.0_dp/dlayer(l)
            cl%wl(r2, r1) = sc*junit(l)
            cl%wl(r2, r2) = sc*junit(l)
            do side = 1, 2
                do k = 1, lays(l)%nfit
                    kf = kf + 1
                    j = (side - 1)*lays(l)%nfit + k
                    cl%w(r1, kf) = lays(l)%pa1(j, 1) - lays(l)%pa1(j, 2)
                    cl%w(r2, kf) = -0.5_dp*sc*dlayer(l)*(lays(l)%pa1(j, 1) &
                                                         + lays(l)%pa1(j, 2))
                end do
            end do
        end do
    end subroutine resistive_closure

end module kin6d_mhd_layer_n
