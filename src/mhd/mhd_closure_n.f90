module kin6d_mhd_closure_n
    !! L6: ideal closure of the per-n response at rational surfaces by
    !! inner-outer matching with the Frobenius solutions of the outer
    !! (ideal) region (Glasser, Greene and Johnson, Phys. Fluids 18 (1975)
    !! 875; the ideal crossing of DCON).
    !!
    !! Outer solutions.  Near a rational surface q = m/n at psin_s, with
    !! x = psin - psin_s, the resonant flux harmonic Phi_m of an ideal outer
    !! solution (kin6d_mhd_response_n; Phi_m = chi' i (n q - m) xi_m in any
    !! straight-field-line chart) is, on each side,
    !!     Phi_m = A |x|^(1/2 - mu) (1 + O(x)) + B |x|^(1/2 + mu) (1 + O(x))
    !!             + C x + O(x^2),
    !! mu = sqrt(-D_I) with the Mercier index D_I of the surface; A is the
    !! large, B the small solution and C x the resonant part of the regular
    !! solutions driven by the sidebands.  With eps = mu - 1/2 and t = |x|/
    !! delta the fit basis on each side is
    !!     b1 = t^(-eps),   b2 = t,   b3 = (t^(1 + eps) - t)/eps,
    !! which spans {large, regular, small} and tends to {1, t, t ln t} as
    !! eps -> 0, the cylindrical zero-beta case where the large solution
    !! carries the x ln|x| term (D_I = -1/4 exactly when p' = 0).
    !!
    !! Closure.  The layer |x| < delta around each rational carries no volume
    !! current (kin6d_mhd_ideal_n sample rings inactive); it is replaced by
    !! two resonant sheets on its edge rings, x = -delta and x = +delta
    !! (kin6d_mhd_response_n sheet_closure_t), whose amplitudes are fixed by
    !! the large-solution coefficients a1 of the fits on nfit rings outside
    !! the layer on each side: ideal shielding sets a1 = 0 on both sides;
    !! prescribed a1 gives the Delta' (matching) data.  The two sheets carry
    !! the net current and the dipole moment of the layer, so the outer
    !! field is that of the ideal solution up to O(delta).
    !!
    !! Shielding current.  GPEC (gpout_singfld) reports the ideal sheet
    !! current as the jump of d Phi_m/d psin between psin_s -+ spot/(n
    !! |dq/dpsin|).  Here the net current of the layer (the sum of the edge
    !! sheets) gives the jump between the layer edges, J(delta), and the
    !! small-solution power carries it to spot: J(spot) = J(delta) (spot/
    !! delta)^(mu - 1/2), up to O(delta) (the regular part C x has no jump).
    !! For mu /= 1/2 GPEC's value therefore depends on spot.  The equivalent
    !! sheet amplitude lambda_eq is that of the sheet on the rational ring
    !! with the same jump (calibrated with the vacuum field of the unit
    !! sheet), and its resonant self flux per area is GPEC's singflx.
    !! With prescribed large coefficients (Delta' data) the jump between the
    !! layer edges also contains A (L'(delta) - L'(-delta)) = A (2 K^2 delta
    !! ln delta + O(delta)), removed by extrapolation in delta.
    !!
    !! Mercier index (DCON mercier_scan, flux-surface averages over the
    !! volume element):
    !!     D_I = -1/4 + T (1 - T) + p1 (v1/(q1 chi1^2))^2 <B^2/|grad s|^2>
    !!           (p1 (<1/B^2> + (2 pi F/chi1)^2 <1/(B^2 |grad s|^2)>) - v2/v1),
    !!     T = 2 pi F p1 v1 <1/|grad s|^2>/(q1 chi1^3),
    !! s = psin, chi1 = 2 pi |psi_edge - psi_axis|, p1 = mu0 dp/ds, v1 = dV/ds,
    !! v2 = d^2V/ds^2, q1 = d|q|/ds, F = |R B_phi|.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_vacuum_surfaces, only: psi_grid_t, surface_t, surface_build, &
        psi_eval, profiles_of_psin
    use kin6d_mhd_response_n, only: ring_angles_t, ring_angles_set, &
        sheet_closure_t
    use kin6d_mhd_ampere_n, only: mu0
    implicit none
    private

    real(dp), parameter :: pi = 3.14159265358979323846_dp

    !> One rational surface and its layer.
    !> Largest admissible |w|_1 of the a1 fit weights (ratio 1.5 gives about
    !! 22, ratio 2 about 8, ratio 1.25 about 70 for nfit = 3).
    real(dp), parameter, public :: layer_wmax = 30.0_dp

    type, public :: rational_layer_t
        integer :: m = 0, nfit = 3
        real(dp) :: psin = 0.0_dp, q = 0.0_dp, dq = 0.0_dp, di = -0.25_dp
        real(dp) :: eps = 0.0_dp, delta = 0.0_dp
        integer :: ring_rat = 0                  !! ring of the rational surface
        integer, allocatable :: ring(:, :)       !! (nfit, 2) fit rings; side 1 x < 0
        real(dp), allocatable :: t(:, :)         !! (nfit, 2) |x|/delta of the rings
        real(dp), allocatable :: pinv(:, :, :)   !! (3, nfit, 2) a = pinv Phi
        real(dp), allocatable :: preg(:, :, :)   !! (3, nfit, 2) smooth fit 1, t, t^2
        real(dp), allocatable :: pjump(:)        !! (2 nfit) two-sided kink fit
        !> (2 nfit) kink of the small solution |x|^(1+eps) at |x| = delta
        !! (layer_kink); pkinks: the same after removing the large-solution
        !! part a1 t^-eps of each side with its Frobenius coefficient a1.
        real(dp), allocatable :: pkink(:), pkinks(:)
        !> (2 nfit, 2) closure weights: a1 of side s = sum_j pa1(j, s) Phi_j
        !! over both sides' fit rings (j = (side - 1) nfit + k).  eps = 0: the
        !! per-side Frobenius fit (rows 1 of pinv, block diagonal).  eps > 0:
        !! joint fit with columns t^-eps (side 1), t^-eps (side 2), sigma t,
        !! t^(1+eps) (side 1), t^(1+eps) (side 2); the even regular kink |x|
        !! of the per-side basis is no local solution for eps > 0 and is
        !! excluded, the odd regular x (smooth) is kept.
        real(dp), allocatable :: pa1(:, :)
    end type rational_layer_t

    public :: mercier_di
    public :: frobenius_basis
    public :: layer_labels
    public :: layer_xmax
    public :: layer_init
    public :: layer_closure
    public :: layer_coefficients
    public :: layer_slope
    public :: layer_regular_slope
    public :: layer_jump
    public :: layer_kink
    public :: layer_large
    public :: layer_inside

contains

    subroutine mercier_di(g, psin, nth, di, q, dq, status)
        !! Mercier index D_I, |q| and d|q|/dpsin of the surface psin.
        type(psi_grid_t), intent(in) :: g
        real(dp), intent(in) :: psin
        integer, intent(in) :: nth
        real(dp), intent(out) :: di, q, dq
        integer, intent(out) :: status
        real(dp) :: h, v1m, v1p, v1, v2, qm, qp, avg(5), f, dfdx, p, dpdx, p1
        real(dp) :: psio, chi1, twopif, term

        status = 1
        if (psin <= 0.0_dp .or. psin >= 1.0_dp) return
        ! centred differences in psin, with the step kept inside the plasma
        ! (a rational closer than 1e-3 to the boundary has no surface at
        ! psin + 1e-3)
        h = min(1.0e-3_dp, 0.5_dp*psin, 0.5_dp*(1 - psin))
        psio = abs(g%psi_edge - g%psi_axis)
        call surface_average(psin - h, v1m, qm, avg, status)
        if (status /= 0) return
        call surface_average(psin + h, v1p, qp, avg, status)
        if (status /= 0) return
        call surface_average(psin, v1, q, avg, status)
        if (status /= 0) return
        v2 = (v1p - v1m)/(2*h)
        dq = (qp - qm)/(2*h)
        call profiles_of_psin(g, psin, f, dfdx, p, dpdx)
        chi1 = 2*pi*psio
        twopif = 2*pi*abs(f)
        p1 = mu0*dpdx
        term = twopif*p1*v1/(dq*chi1**3)*avg(2)
        di = -0.25_dp + term*(1 - term) + p1*(v1/(dq*chi1**2))**2*avg(1) &
            *(p1*(avg(3) + (twopif/chi1)**2*avg(4)) - v2/v1)
    contains
        subroutine surface_average(x, vp, qa, av, st)
            !! dV/dpsin, |q| and the volume averages of B^2/|grad s|^2,
            !! 1/|grad s|^2, 1/B^2, 1/(B^2 |grad s|^2), B^2 on surface x.
            real(dp), intent(in) :: x
            real(dp), intent(out) :: vp, qa, av(5)
            integer, intent(out) :: st
            type(surface_t) :: s
            real(dp) :: ps, pr, pz, gs2, b2, w, wsum, ff, df, pp, dp1
            integer :: j
            call surface_build(g, x, nth, s, st)
            if (st /= 0) return
            call profiles_of_psin(g, x, ff, df, pp, dp1)
            qa = abs(s%q)
            av = 0.0_dp
            wsum = 0.0_dp
            do j = 1, nth
                call psi_eval(g, s%r(j), s%z(j), ps, pr, pz)
                gs2 = (pr**2 + pz**2)/psio**2
                b2 = ((pr**2 + pz**2) + ff**2)/s%r(j)**2
                ! volume element R dl/|grad s| per d theta_g
                w = s%r(j)*s%dl(j)/sqrt(gs2)
                wsum = wsum + w
                av = av + w*[b2/gs2, 1/gs2, 1/b2, 1/(b2*gs2), b2]
            end do
            av = av/wsum
            vp = 2*pi*wsum*2*pi/nth
        end subroutine surface_average
    end subroutine mercier_di

    pure subroutine frobenius_basis(eps, t, f, df)
        !! b1 = t^(-eps), b2 = t, b3 = (t^(1 + eps) - t)/eps and their
        !! t-derivatives (eps -> 0: b3 = t ln t).
        real(dp), intent(in) :: eps, t
        real(dp), intent(out) :: f(3), df(3)
        real(dp) :: lt, e1

        lt = log(t)
        if (abs(eps*lt) > 1.0e-8_dp) then
            e1 = expm1(eps*lt)/eps
        else
            e1 = lt*(1 + 0.5_dp*eps*lt)
        end if
        f(1) = exp(-eps*lt)
        df(1) = -eps*f(1)/t
        f(2) = t
        df(2) = 1.0_dp
        f(3) = t*e1
        df(3) = e1 + exp(eps*lt)
    contains
        pure real(dp) function expm1(z)
            real(dp), intent(in) :: z
            if (abs(z) < 1.0e-5_dp) then
                expm1 = z*(1 + z*(0.5_dp + z/6))
            else
                expm1 = exp(z) - 1
            end if
        end function expm1
    end subroutine frobenius_basis

    pure function layer_xmax(psin_s, others, delta, ratio, nfit, nrad) result(xmax)
        !! Largest layer-ring offset: grade up to the uniform spacing
        !! 2 sqrt(psin_s)/nrad, but never below the outermost fit ring
        !! delta ratio^(nfit-1), and at most 0.45 of the gap to the other
        !! rationals (others).
        real(dp), intent(in) :: psin_s, others(:), delta, ratio
        integer, intent(in) :: nfit, nrad
        real(dp) :: xmax
        integer :: k

        xmax = 2*sqrt(psin_s)/(nrad*(ratio - 1))
        xmax = max(xmax, delta*ratio**(nfit - 1)*(1 + 1.0e-9_dp))
        do k = 1, size(others)
            xmax = min(xmax, 0.45_dp*abs(others(k) - psin_s))
        end do
    end function layer_xmax

    subroutine layer_labels(psin_s, delta, ratio, xmax, lo, hi, psin)
        !! psin labels of the layer rings: psin_s -+ delta ratio^k, k = 0, 1,
        !! ..., while delta ratio^k <= xmax and the label lies in (lo, hi).
        real(dp), intent(in) :: psin_s, delta, ratio, xmax, lo, hi
        real(dp), allocatable, intent(out) :: psin(:)
        real(dp) :: x, buf(400)
        integer :: c

        c = 0
        x = delta
        do while (x <= xmax .and. c < size(buf) - 1)
            if (psin_s - x > lo) then
                c = c + 1
                buf(c) = psin_s - x
            end if
            if (psin_s + x < hi) then
                c = c + 1
                buf(c) = psin_s + x
            end if
            x = x*ratio
        end do
        psin = buf(1:c)
    end subroutine layer_labels

    subroutine layer_init(ul, psin_s, q, dq, di, m, delta, ratio, nfit, lay, status)
        !! Layer of the rational at psin_s on the mesh with ring labels ul
        !! (u = sqrt(psin)): fit rings at psin_s -+ delta ratio^(k-1),
        !! k = 1 .. nfit (each must be a ring), and the fit operators.
        real(dp), intent(in) :: ul(:), psin_s, q, dq, di, delta, ratio
        integer, intent(in) :: m, nfit
        type(rational_layer_t), intent(out) :: lay
        integer, intent(out) :: status
        real(dp) :: bm(nfit, 3), br(nfit, 3), f(3), df(3), x, mu
        real(dp) :: bj(2*nfit, 5), pj(5, 2*nfit), bk(2*nfit, 5), pk(5, 2*nfit)
        real(dp) :: wk(nfit), ba(2*nfit, 5), pa(5, 2*nfit), sg, tk
        integer :: k, side, j

        status = 1
        if (nfit < 3 .or. di >= 0.0_dp) return
        lay%m = m
        lay%nfit = nfit
        lay%psin = psin_s
        lay%q = q
        lay%dq = dq
        lay%di = di
        mu = sqrt(-di)
        lay%eps = mu - 0.5_dp
        lay%delta = delta
        lay%ring_rat = minloc(abs(ul - sqrt(psin_s)), 1)
        allocate (lay%ring(nfit, 2), lay%t(nfit, 2), lay%pinv(3, nfit, 2), &
            lay%preg(3, nfit, 2), lay%pjump(2*nfit), lay%pkink(2*nfit), lay%pkinks(2*nfit), &
            lay%pa1(2*nfit, 2))
        do side = 1, 2
            do k = 1, nfit
                x = delta*ratio**(k - 1)
                j = minloc(abs(ul**2 - (psin_s + (2*side - 3)*x)), 1)
                if (abs(ul(j)**2 - (psin_s + (2*side - 3)*x)) > 1.0e-6_dp*x) return
                lay%ring(k, side) = j
                lay%t(k, side) = abs(ul(j)**2 - psin_s)/delta
                call frobenius_basis(lay%eps, lay%t(k, side), f, df)
                bm(k, :) = f
                br(k, :) = [1.0_dp, lay%t(k, side), lay%t(k, side)**2]
            end do
            call pseudo_inverse(bm, lay%pinv(:, :, side), status)
            if (status /= 0) return
            ! a1 is a difference of fluxes on neighbouring rings; the flux error
            ! is amplified by |w|_1 / delta in the sheet dipole (ratio 1.25: 70)
            if (sum(abs(lay%pinv(1, :, side))) > layer_wmax) then
                status = 2
                return
            end if
            call pseudo_inverse(br, lay%preg(:, :, side), status)
            if (status /= 0) return
        end do
        lay%pa1 = 0.0_dp
        if (lay%eps > 0.0_dp) then
            do side = 1, 2
                sg = real(2*side - 3, dp)
                do k = 1, nfit
                    tk = lay%t(k, side)
                    ba((side - 1)*nfit + k, :) = 0.0_dp
                    ba((side - 1)*nfit + k, side) = tk**(-lay%eps)
                    ba((side - 1)*nfit + k, 3) = sg*tk
                    ba((side - 1)*nfit + k, 3 + side) = tk**(1 + lay%eps)
                end do
            end do
            call pseudo_inverse_n(ba, pa, status)
            if (status /= 0) return
            lay%pa1 = transpose(pa(1:2, :))
        else
            do side = 1, 2
                lay%pa1((side - 1)*nfit + 1:side*nfit, side) = lay%pinv(1, :, side)
            end do
        end if
        if (maxval(sum(abs(lay%pa1), 1)) > layer_wmax) then
            status = 2
            return
        end if
        ! symmetric two-sided fit of the regular part with a slope kink at x = 0
        ! (columns 1, sigma t, t, t^2, sigma t^2; the kink is column 3)
        do side = 1, 2
            do k = 1, nfit
                bj((side - 1)*nfit + k, :) = [1.0_dp, (2*side - 3)*lay%t(k, side), &
                    lay%t(k, side), lay%t(k, side)**2, &
                    (2*side - 3)*lay%t(k, side)**2]
            end do
        end do
        call pseudo_inverse_n(bj, pj, status)
        if (status /= 0) return
        lay%pjump = pj(3, :)
        ! the same fit with the small solution |x|^(1+eps) of the Mercier
        ! exponent as the kink column and the even large solution t^-eps in
        ! place of the constant (columns t^-eps, sigma t, t^(1+eps), t^2,
        ! sigma t^(2+eps); eps = 0: the polynomial kink fit).  With the
        ! constant column a residual large solution A|x|^-eps (eps > 0) leaks
        ! into the kink as A delta^(-1-2 eps); the odd large part sigma t^-eps
        ! drops out on the symmetric rings.  For pkinks the large solution
        ! t^-eps of each side is first removed with its closure coefficient
        ! a1 (pa1);
        ! for a1 = 0 (ideal) the removal would only inject the a1 error of the
        ! small-solution corrections, so pkink uses the raw fluxes
        do side = 1, 2
            do k = 1, nfit
                x = lay%t(k, side)
                bk((side - 1)*nfit + k, :) = [x**(-lay%eps), (2*side - 3)*x, &
                    x**(1 + lay%eps), x**2, (2*side - 3)*x**(2 + lay%eps)]
            end do
        end do
        call pseudo_inverse_n(bk, pk, status)
        if (status /= 0) return
        lay%pkink = pk(3, :)
        ! w (I - sum_side f1_side pa1_side^T): weights on the raw fluxes
        lay%pkinks = lay%pkink
        do side = 1, 2
            wk = pk(3, (side - 1)*nfit + 1:side*nfit)
            lay%pkinks = lay%pkinks - sum(wk*lay%t(:, side)**(-lay%eps))*lay%pa1(:, side)
        end do
        status = 0
    end subroutine layer_init

    subroutine pseudo_inverse_n(a, p, status)
        !! p = (a^T a)^{-1} a^T for a (nrow, ncol) of full column rank.
        real(dp), intent(in) :: a(:, :)
        real(dp), intent(out) :: p(:, :)
        integer, intent(out) :: status
        real(dp), allocatable :: nn(:, :), rhs(:, :)
        integer, allocatable :: ipiv(:)
        integer :: nc

        nc = size(a, 2)
        allocate (nn(nc, nc), rhs(nc, size(a, 1)), ipiv(nc))
        nn = matmul(transpose(a), a)
        rhs = transpose(a)
        call dgesv(nc, size(a, 1), nn, nc, ipiv, rhs, nc, status)
        p = rhs
    end subroutine pseudo_inverse_n

    subroutine pseudo_inverse(a, p, status)
        !! p = (a^T a)^{-1} a^T for a (nfit, 3) of full column rank.
        real(dp), intent(in) :: a(:, :)
        real(dp), intent(out) :: p(:, :)
        integer, intent(out) :: status
        real(dp) :: n3(3, 3), rhs(3, size(a, 1))
        integer :: ipiv(3)

        n3 = matmul(transpose(a), a)
        rhs = transpose(a)
        call dgesv(3, size(a, 1), n3, 3, ipiv, rhs, 3, status)
        p = rhs
    end subroutine pseudo_inverse

    subroutine layer_closure(g, ul, lays, nfine, angle, cl, status)
        !! Sheet closure of the layers: sheets on the inner and outer edge
        !! rings, flux rings = fit rings, constraints a1 = target (0 for
        !! ideal shielding; set cl%target for Delta' data).
        type(psi_grid_t), intent(in) :: g
        real(dp), intent(in) :: ul(:)
        type(rational_layer_t), intent(in) :: lays(:)
        integer, intent(in) :: nfine, angle
        type(sheet_closure_t), intent(out) :: cl
        integer, intent(out) :: status
        type(surface_t) :: s
        integer :: l, side, k, nl, nfit, kf, st

        status = 1
        nl = size(lays)
        cl%ns = 2*nl
        cl%nf = 0
        do l = 1, nl
            cl%nf = cl%nf + 2*lays(l)%nfit
        end do
        allocate (cl%sheet_ring(cl%ns), cl%sheet_m(cl%ns), cl%sheet_ra(cl%ns), &
            cl%flux_ring(cl%nf), cl%flux_m(cl%nf), cl%flux_ra(cl%nf), &
            cl%w(cl%ns, cl%nf), cl%target(cl%ns))
        cl%w = 0.0_dp
        cl%target = (0.0_dp, 0.0_dp)
        kf = 0
        do l = 1, nl
            nfit = lays(l)%nfit
            do side = 1, 2
                do k = 1, nfit
                    kf = kf + 1
                    cl%flux_ring(kf) = lays(l)%ring(k, side)
                    cl%flux_m(kf) = lays(l)%m
                    call surface_build(g, ul(lays(l)%ring(k, side))**2, nfine, s, st)
                    if (st /= 0) return
                    call ring_angles_set(s%q, s%theta_s(:, angle), s%nu, cl%flux_ra(kf))
                    cl%w(2*(l - 1) + 1, kf) = lays(l)%pa1((side - 1)*nfit + k, 1)
                    cl%w(2*(l - 1) + 2, kf) = lays(l)%pa1((side - 1)*nfit + k, 2)
                    if (k == 1) then
                        cl%sheet_ring(2*(l - 1) + side) = lays(l)%ring(1, side)
                        cl%sheet_m(2*(l - 1) + side) = lays(l)%m
                        cl%sheet_ra(2*(l - 1) + side) = cl%flux_ra(kf)
                    end if
                end do
            end do
        end do
        status = 0
    end subroutine layer_closure

    pure function layer_coefficients(lay, phi) result(a)
        !! Fit coefficients a(:, side) of the basis b1, b2, b3 from the
        !! resonant fluxes phi(k, side) on the fit rings.
        type(rational_layer_t), intent(in) :: lay
        complex(dp), intent(in) :: phi(:, :)
        complex(dp) :: a(3, 2)
        integer :: side

        do side = 1, 2
            a(:, side) = matmul(lay%pinv(:, :, side), phi(1:lay%nfit, side))
        end do
    end function layer_coefficients

    pure function layer_slope(lay, a, ax) result(d)
        !! d Phi/d|x| (per unit psin) of the fitted outer solution with
        !! coefficients a (one side) at |x| = ax.
        type(rational_layer_t), intent(in) :: lay
        complex(dp), intent(in) :: a(3)
        real(dp), intent(in) :: ax
        complex(dp) :: d
        real(dp) :: f(3), df(3)

        call frobenius_basis(lay%eps, ax/lay%delta, f, df)
        d = sum(a*df)/lay%delta
    end function layer_slope

    pure function layer_regular_slope(lay, phi) result(d)
        !! One-sided slopes d Phi/d|x| at x = 0 of a field smooth on each
        !! side (quadratic fit), summed over both sides: the jump of
        !! d Phi/d psin across the rational surface.
        type(rational_layer_t), intent(in) :: lay
        complex(dp), intent(in) :: phi(:, :)
        complex(dp) :: d
        integer :: side

        d = (0.0_dp, 0.0_dp)
        do side = 1, 2
            d = d + sum(lay%preg(2, :, side)*phi(1:lay%nfit, side))/lay%delta
        end do
    end function layer_regular_slope

    pure function layer_jump(lay, phi) result(d)
        !! Jump of d Phi/d psin across the rational surface from a joint
        !! two-sided least-squares fit of all fit rings (inner and outer) by
        !! c0 + c1 x + (J/2) |x| + c2 x^2 + c3 x |x|: the kink J is one
        !! coefficient of the joint fit, not a difference of one-sided fits.
        type(rational_layer_t), intent(in) :: lay
        complex(dp), intent(in) :: phi(:, :)
        complex(dp) :: d
        integer :: n

        n = lay%nfit
        d = 2*(sum(lay%pjump(1:n)*phi(1:n, 1)) &
               + sum(lay%pjump(n + 1:2*n)*phi(1:n, 2)))/lay%delta
    end function layer_jump

    pure function layer_kink(lay, phi, subtract_large) result(d)
        !! Jump of d Phi/d psin across the rational surface carried by the
        !! small solution |x|^(1/2 + mu) = |x|^(1+eps), evaluated at |x| = delta:
        !! 2 (1 + eps) c/delta with c the coefficient of t^(1+eps) in the joint
        !! two-sided fit (columns t^-eps, sigma t, t^(1+eps), t^2, sigma t^(2+eps)).
        !! Multiply by (x/delta)^eps for the jump at |x| = x; at eps = 0 this is
        !! layer_jump.  With a1 /= 0 (resistive layer) subtract_large removes
        !! the large solution a1 t^-eps of each side first; for a1 = 0 (ideal
        !! shielding) leave it off, the removal would only add the a1 error.
        type(rational_layer_t), intent(in) :: lay
        complex(dp), intent(in) :: phi(:, :)
        logical, intent(in), optional :: subtract_large
        complex(dp) :: d
        real(dp), allocatable :: w(:)
        integer :: n

        n = lay%nfit
        w = lay%pkink
        if (present(subtract_large)) then
            if (subtract_large) w = lay%pkinks
        end if
        d = 2*(1 + lay%eps)*(sum(w(1:n)*phi(1:n, 1)) &
                             + sum(w(n + 1:2*n)*phi(1:n, 2)))/lay%delta
    end function layer_kink

    pure function layer_large(lay, phi) result(a1)
        !! Large-solution coefficients a1(side) of the closure (pa1) from the
        !! resonant fluxes phi(k, side) on the fit rings of both sides.
        type(rational_layer_t), intent(in) :: lay
        complex(dp), intent(in) :: phi(:, :)
        complex(dp) :: a1(2)
        integer :: side, n

        n = lay%nfit
        do side = 1, 2
            a1(side) = sum(lay%pa1(1:n, side)*phi(1:n, 1)) &
                + sum(lay%pa1(n + 1:2*n, side)*phi(1:n, 2))
        end do
    end function layer_large

    pure logical function layer_inside(lays, ul, u)
        !! u lies strictly between the edge rings of one of the layers.
        type(rational_layer_t), intent(in) :: lays(:)
        real(dp), intent(in) :: ul(:), u
        integer :: l

        layer_inside = .false.
        do l = 1, size(lays)
            if (u > ul(lays(l)%ring(1, 1)) .and. u < ul(lays(l)%ring(1, 2))) &
                layer_inside = .true.
        end do
    end function layer_inside

end module kin6d_mhd_closure_n
