!> QG6 binary Coulomb collision kernel in the relative coordinate.
!>
!> Model: two nonrelativistic point charges, reduced mass mu, coupling g > 0
!> (repulsive; electron-electron), relative speed v, k = mu v/hbar,
!> eta = g/(hbar v), b90 = g/(mu v^2) = eta/k, lambdabar = 1/k. Screening is
!> either a sharp impact-parameter (classical) / partial-wave (quantum) cut at
!> lambda_D, kappa = k lambda_D, or the static Yukawa potential
!> g exp(-r/lambda_D)/r (gen/kin6d_symbolic_physics.f90).
!>
!> Every Coulomb logarithm is the momentum-transfer cross section normalized
!> as lnL = sigma_tr/(4 pi b90^2); the Landau/Spitzer rates use it linearly.
!>
!> Pure Coulomb: the quantum differential cross section of distinguishable
!> particles equals Rutherford (ce_rutherford) at every angle; identical
!> particles add the exchange and Mott interference terms (ce_mott). In
!> partial waves sigma_{l+1} - sigma_l = atan(eta/(l+1)) exactly, so the
!> quantum deflection of ring l is the Rutherford deflection at the quantized
!> impact parameter (l+1)/k; the momentum-transfer sum is local in l:
!>   k^2 sigma_tr/(4 pi) = sum_l (l+1) sin^2(delta_l - delta_{l+1}).
!> The quantum diffraction correction to lnL (classical b90 -> lambdabar) is
!> therefore carried by the near rings l < O(max(eta,1)).
!>
!> Yukawa: phase shifts delta_l from the shared radial FEM-DVR
!> (src/radial_dvr.f90) at fixed energy with the screened radial operator,
!> Dirichlet data deep inside the centrifugal barrier and at a bridge node
!> beyond r_turn + tail*lambda_D, matched to Riccati-Bessel functions by least
!> squares. d delta_l/d eta from the same factorization (implicit equation-level
!> tangent (H - E) da = -(dV/dg) a), not by differentiating iterations.
module kin6d_quantum_collision_kernel
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_special_gamma, only: digamma_complex, trigamma_complex
    use fortnum_special_riccati_bessel, only: riccati_bessel
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortnum_quadrature, only: gauss_legendre, gauss_legendre_ab
    use kin6d_radial_dvr, only: radial_mesh_t, build_radial_mesh, mesh_spacing
    use kin6d_coulomb_encounter_leaves, only: ce_mott
    use kin6d_quantum_collision_leaves
    implicit none
    private
    public :: SPIN_DISTINGUISHABLE, SPIN_SINGLET, SPIN_TRIPLET, SPIN_UNPOLARIZED
    public :: SPIN_IDENTICAL_CLASSICAL
    public :: complex_digamma, quantum_cut_log, classical_cut_log, born_yukawa_log
    public :: mott_transfer_log, coulomb_dsdo
    public :: classical_yukawa_deflection, classical_yukawa_log
    public :: yukawa_config_t, yukawa_plan_t, yukawa_ledger_t
    public :: yukawa_plan, yukawa_phases, yukawa_transfer_log
    public :: transfer_from_phases, viscosity_from_phases, riccati_bessel
    public :: coulomb_ring_viscosity, coulomb_viscosity_exact

    integer, parameter :: SPIN_DISTINGUISHABLE = 0, SPIN_SINGLET = 1, SPIN_TRIPLET = 2
    integer, parameter :: SPIN_UNPOLARIZED = 3, SPIN_IDENTICAL_CLASSICAL = 4
    real(dp), parameter :: pi = 3.14159265358979323846_dp
    real(dp), parameter :: euler_gamma = 0.57721566490153286061_dp

    !> Yukawa partial-wave solver design (dimensionless units k = hbar = mu = 1,
    !> so the coupling equals eta and lambda_D equals kappa).
    type :: yukawa_config_t
        integer :: order = 12                  ! GLL points per element
        real(dp) :: nodes_per_wavelength = 10  ! resolution of k = 1
        real(dp) :: tail = 28                  ! matching beyond r_turn, units lambda_D
        real(dp) :: barrier = 36               ! WKB depth of the inner Dirichlet node
        real(dp) :: tolerance = 1.0e-10_dp     ! Born-asymptotic l-tail of lnL
        integer :: lcap = 200000
    end type

    !> Per-solve error ledger (all indicators; no rigorous continuum bound).
    type :: yukawa_ledger_t
        real(dp) :: tail_estimate = 0      ! Born-asymptotic omitted l >= lmax, lnL units
        real(dp) :: fit_residual = 0       ! max relative least-squares matching residual
        real(dp) :: potential_at_match = 0 ! max |V/E| over the matching elements
        real(dp) :: barrier_weight = 0     ! max exp(-2 WKB depth) at the inner node
        integer :: lmax = 0, nodes = 0
    end type

    !> Frozen design: mesh and, per partial wave, the inner and matching nodes.
    !> Finite differences in eta must reuse the plan.
    type :: yukawa_plan_t
        type(yukawa_config_t) :: cfg
        type(radial_mesh_t) :: mesh
        real(dp) :: kappa = 0, eta_design = 0
        integer :: lmax = 0
        integer, allocatable :: inner(:), outer(:)
        type(yukawa_ledger_t) :: ledger
    end type

    interface
        subroutine dgbtrf(m, n, kl, ku, ab, ldab, ipiv, info)
            import :: dp
            integer, intent(in) :: m, n, kl, ku, ldab
            real(dp), intent(inout) :: ab(ldab, *)
            integer, intent(out) :: ipiv(*), info
        end subroutine dgbtrf
        subroutine dgbtrs(trans, n, kl, ku, nrhs, ab, ldab, ipiv, b, ldb, info)
            import :: dp
            character, intent(in) :: trans
            integer, intent(in) :: n, kl, ku, nrhs, ldab, ldb
            real(dp), intent(in) :: ab(ldab, *)
            integer, intent(in) :: ipiv(*)
            real(dp), intent(inout) :: b(ldb, *)
            integer, intent(out) :: info
        end subroutine dgbtrs
    end interface
contains
    !> Complex digamma psi(z) and trigamma psi'(z) (FortNum).
    elemental subroutine complex_digamma(z, psi, dpsi)
        complex(dp), intent(in) :: z
        complex(dp), intent(out) :: psi, dpsi
        psi = digamma_complex(z); dpsi = trigamma_complex(z)
    end subroutine complex_digamma

    !> Classical Rutherford transfer with a sharp impact-parameter cut at
    !> lambda_D: lnL = ln(1 + kappa^2/eta^2)/2 and d lnL/d eta.
    elemental subroutine classical_cut_log(eta, kappa, lnl, dlnl)
        real(dp), intent(in) :: eta, kappa
        real(dp), intent(out) :: lnl, dlnl
        call qc_cut_log_classical(eta, kappa, lnl)
        call qc_cut_log_classical_jvp(eta, kappa, 1.0_dp, 0.0_dp, dlnl)
    end subroutine classical_cut_log

    !> Quantum pure-Coulomb transfer with the matching partial-wave cut: Langer
    !> rings sqrt(l(l+1)) <= k b < sqrt((l+1)(l+2)) inside k b <= kappa, the last
    !> ring weighted by its area fraction. lnL = sum_n n/(n^2+eta^2) (n = l+1);
    !> long sums use Re[psi(N+1+i eta) - psi(1+i eta)]. dlnl = d lnL/d eta.
    elemental subroutine quantum_cut_log(eta, kappa, lnl, dlnl)
        real(dp), intent(in) :: eta, kappa
        real(dp), intent(out) :: lnl, dlnl
        complex(dp) :: p1, d1, p2, d2
        real(dp) :: w, t, dt
        integer :: nfull, n
        nfull = floor((sqrt(1 + 4*kappa**2) - 1)/2)
        do while (real(nfull + 1, dp)*(nfull + 2) <= kappa**2)
            nfull = nfull + 1
        end do
        do while (nfull > 0 .and. real(nfull, dp)*(nfull + 1) > kappa**2)
            nfull = nfull - 1
        end do
        lnl = 0; dlnl = 0
        if (nfull <= 4000) then
            do n = nfull, 1, -1
                call qc_ring_transfer(eta, real(n, dp), t)
                call qc_ring_transfer_jvp(eta, real(n, dp), 1.0_dp, 0.0_dp, dt)
                lnl = lnl + t; dlnl = dlnl + dt
            end do
        else
            call complex_digamma(cmplx(nfull + 1, eta, dp), p2, d2)
            call complex_digamma(cmplx(1, eta, dp), p1, d1)
            lnl = real(p2 - p1)
            dlnl = -aimag(d2) + aimag(d1)
        end if
        call qc_partial_ring(kappa, real(nfull + 1, dp), w)
        call qc_ring_transfer(eta, real(nfull + 1, dp), t)
        call qc_ring_transfer_jvp(eta, real(nfull + 1, dp), 1.0_dp, 0.0_dp, dt)
        lnl = lnl + w*t; dlnl = dlnl + w*dt
    end subroutine quantum_cut_log

    !> First Born Yukawa transfer, closed form (eta-independent normalized).
    elemental real(dp) function born_yukawa_log(kappa) result(lnl)
        real(dp), intent(in) :: kappa
        call qc_born_yukawa_log(kappa, lnl)
    end function born_yukawa_log

    !> Pure-Coulomb differential cross section in units of 1/k^2 (k = 1):
    !> distinguishable Rutherford, or identical singlet/triplet/unpolarized Mott
    !> (laboratory convention: both outgoing particles counted, 0 <= theta <= pi).
    elemental real(dp) function coulomb_dsdo(theta, eta, spin) result(dsdo)
        real(dp), intent(in) :: theta, eta
        integer, intent(in) :: spin
        real(dp) :: d, x, c
        call ce_mott(theta, eta, 1.0_dp, d, x, c)
        select case (spin)
        case (SPIN_SINGLET)
            dsdo = d + x + c
        case (SPIN_TRIPLET)
            dsdo = d + x - c
        case (SPIN_UNPOLARIZED)
            dsdo = d + x - c/2
        case (SPIN_IDENTICAL_CLASSICAL)
            dsdo = d + x
        case default
            dsdo = d
        end select
    end function coulomb_dsdo

    !> Momentum-transfer logarithm of the pure-Coulomb kernel with the small-angle
    !> cut theta_c = 2 atan(eta/kappa) (the classical impact cut at lambda_D).
    !> Identical particles follow the forward-going (tagged) particle:
    !> theta in [theta_c, pi/2] of the symmetrized cross section; distinguishable
    !> particles use [theta_c, pi]. Quadrature in t = ln tan(theta/2) on panels
    !> resolving the Mott phase cos(2 eta t); nodes per panel declared.
    function mott_transfer_log(eta, kappa, spin, panel_nodes) result(lnl)
        real(dp), intent(in) :: eta, kappa
        integer, intent(in) :: spin
        integer, intent(in), optional :: panel_nodes
        real(dp) :: lnl, t0, t1, width, a, b, theta, jac
        real(dp), allocatable :: x(:), w(:)
        integer :: np, npanel, ip, i
        np = 24
        if (present(panel_nodes)) np = panel_nodes
        allocate (x(np), w(np))
        t0 = log(eta/kappa)
        t1 = 0
        if (spin == SPIN_DISTINGUISHABLE) t1 = 6   ! tan(theta/2) up to e^6 (theta ~ pi)
        width = min(0.5_dp, pi/(4*max(eta, 1.0e-12_dp)))
        npanel = max(1, ceiling((t1 - t0)/width))
        lnl = 0
        do ip = 1, npanel
            a = t0 + (t1 - t0)*(ip - 1)/npanel; b = t0 + (t1 - t0)*ip/npanel
            call gauss_legendre_ab(np, a, b, x, w)
            do i = 1, np
                theta = 2*atan(exp(x(i)))
                jac = sin(theta)     ! d theta/d t = sin(theta)
                lnl = lnl + w(i)*jac*2*sin(theta/2)**2*coulomb_dsdo(theta, eta, spin)* &
                      2*pi*sin(theta)
            end do
        end do
        if (spin == SPIN_DISTINGUISHABLE) lnl = lnl + tail_back(t1)
        lnl = lnl/(4*pi*eta**2)
    contains
        ! Rutherford remainder theta in [2 atan(e^t1), pi], closed form:
        ! 4 pi b90^2 [ln sin(theta/2)] with b90 = eta (k = 1).
        real(dp) function tail_back(t) result(s)
            real(dp), intent(in) :: t
            s = -4*pi*eta**2*log(sin(atan(exp(t))))
        end function tail_back
    end function mott_transfer_log

    !> Classical deflection in the repulsive Yukawa potential, lengths in units
    !> of b90 (so V/E = 2 exp(-r/lam)/r, lam = lambda_D/b90):
    !> Theta = pi - 2 b int_0^{u0} du/sqrt(F(u)), F = 1 - b^2 u^2 - 2u exp(-1/(u lam)),
    !> with u = u0 (1 - s^2) removing the turning-point singularity; nq
    !> Gauss-Legendre nodes in s.
    function classical_yukawa_deflection(b, lam, nq) result(theta)
        real(dp), intent(in) :: b, lam
        integer, intent(in) :: nq
        real(dp) :: theta, lo, hi, mid, f, u0, s(nq), w(nq), integral, u
        integer :: it, i
        lo = 0; hi = 1
        call qc_yukawa_radial(hi, b, 1.0_dp, 0.5_dp, lam, f)
        do while (f > 0)
            hi = 2*hi
            call qc_yukawa_radial(hi, b, 1.0_dp, 0.5_dp, lam, f)
        end do
        do it = 1, 200
            mid = (lo + hi)/2
            if (.not. (mid > lo .and. mid < hi)) exit
            call qc_yukawa_radial(mid, b, 1.0_dp, 0.5_dp, lam, f)
            if (f > 0) then
                lo = mid
            else
                hi = mid
            end if
        end do
        u0 = lo
        call gauss_legendre_ab(nq, 0.0_dp, 1.0_dp, s, w)
        integral = 0
        do i = 1, nq
            u = u0*(1 - s(i)**2)
            call qc_yukawa_radial(u, b, 1.0_dp, 0.5_dp, lam, f)
            integral = integral + w(i)*2*u0*s(i)/sqrt(max(f, tiny(1.0_dp)))
        end do
        theta = pi - 2*b*integral
    end function classical_yukawa_deflection

    !> Classical Yukawa transfer log (1/2) int (1 - cos Theta) b db in b90 units,
    !> lam = lambda_D/b90 = kappa/eta; ln b panels on [1e-6, 30 lam].
    function classical_yukawa_log(lam, nq, per_unit) result(lnl)
        real(dp), intent(in) :: lam
        integer, intent(in) :: nq, per_unit
        real(dp) :: lnl, t0, t1, x(8), w(8), b, th
        integer :: np, ip, i
        t0 = log(1.0e-6_dp); t1 = log(30*lam)
        np = max(1, ceiling((t1 - t0)*per_unit/8))
        lnl = (1.0e-6_dp)**2/2   ! b < 1e-6: Theta ~ pi, (1/2) int 2 b db
        do ip = 1, np
            call gauss_legendre_ab(8, t0 + (t1 - t0)*(ip - 1)/np, t0 + (t1 - t0)*ip/np, x, w)
            do i = 1, 8
                b = exp(x(i))
                th = classical_yukawa_deflection(b, lam, nq)
                lnl = lnl + w(i)*b*b*sin(th/2)**2   ! (1/2)(1-cos) b db = sin^2 b^2 dt
            end do
        end do
    end function classical_yukawa_log

    !> k^2 sigma_tr/(4 pi) = sum_l (l+1) sin^2(delta_l - delta_{l+1}).
    pure real(dp) function transfer_from_phases(delta) result(s)
        real(dp), intent(in) :: delta(0:)
        integer :: l
        s = 0
        do l = ubound(delta, 1) - 1, 0, -1
            s = s + (l + 1)*sin(delta(l) - delta(l + 1))**2
        end do
    end function transfer_from_phases

    !> k^2 sigma^(2)/(4 pi), sigma^(2) = int sin^2(theta) dsigma:
    !> sum_l (l+1)(l+2)/(2l+3) sin^2(delta_l - delta_{l+2}).
    pure real(dp) function viscosity_from_phases(delta) result(s)
        real(dp), intent(in) :: delta(0:)
        integer :: l
        s = 0
        do l = ubound(delta, 1) - 2, 0, -1
            s = s + real((l + 1)*(l + 2), dp)/(2*l + 3)*sin(delta(l) - delta(l + 2))**2
        end do
    end function viscosity_from_phases

    !> Pure Coulomb, rings n = 1..nmax, sigma^(2) in units 4 pi b90^2:
    !> exact partial waves (delta_{l+2}-delta_l = atan(eta/(l+1))+atan(eta/(l+2)))
    !> and the l-local ring kernel (deflection chi_n = 2 atan(eta/n), ring area
    !> 2 pi n/k^2).
    pure real(dp) function coulomb_viscosity_exact(eta, nmax) result(s)
        real(dp), intent(in) :: eta
        integer, intent(in) :: nmax
        integer :: l
        s = 0
        do l = nmax - 1, 0, -1
            s = s + real((l + 1)*(l + 2), dp)/(2*l + 3)* &
                sin(atan(eta/(l + 1)) + atan(eta/(l + 2)))**2
        end do
        s = s/eta**2
    end function coulomb_viscosity_exact

    pure real(dp) function coulomb_ring_viscosity(eta, nmax) result(s)
        real(dp), intent(in) :: eta
        integer, intent(in) :: nmax
        real(dp) :: chi
        integer :: n
        s = 0
        do n = nmax, 1, -1
            call qc_ring_deflection(eta, real(n, dp), chi)
            s = s + n*sin(chi)**2
        end do
        s = s/(2*eta**2)
    end function coulomb_ring_viscosity

    !> Effective radial potential (k = hbar = mu = 1, E = 1/2) of the Yukawa pair.
    elemental real(dp) function veff(r, l, eta, kappa) result(v)
        real(dp), intent(in) :: r, eta, kappa
        integer, intent(in) :: l
        call qc_screened_radial_potential(r, real(l, dp), eta, 1.0_dp, 1.0_dp, kappa, v)
    end function veff

    !> Classical turning radius veff(r) = 1/2 (largest root; veff decreasing).
    real(dp) function turning_radius(l, eta, kappa) result(rt)
        integer, intent(in) :: l
        real(dp), intent(in) :: eta, kappa
        real(dp) :: lo, hi, mid
        integer :: it
        lo = 0; hi = eta + sqrt(eta**2 + real(l, dp)*(l + 1)) + 1
        do it = 1, 200
            mid = (lo + hi)/2
            if (.not. (mid > lo .and. mid < hi)) exit
            if (veff(mid, l, eta, kappa) > 0.5_dp) then
                lo = mid
            else
                hi = mid
            end if
        end do
        rt = hi
    end function turning_radius

    !> Frozen design for a Yukawa solve at (eta_design, kappa).
    subroutine yukawa_plan(plan, eta, kappa, cfg, ok)
        type(yukawa_plan_t), intent(out) :: plan
        real(dp), intent(in) :: eta, kappa
        type(yukawa_config_t), intent(in) :: cfg
        logical, intent(out) :: ok
        real(dp) :: alpha, q, term, tail, h, rmax, rt, depth, r1, r2, k1, k2
        integer :: l, deg, nelem, j, jt
        ok = .false.
        if (.not. (eta > 0 .and. kappa > 0 .and. cfg%order >= 4)) return
        plan%cfg = cfg; plan%kappa = kappa; plan%eta_design = eta
        ! Born-asymptotic l tail: Q_l(cosh alpha) ~ sqrt(pi/(2(l+1/2) sinh alpha))
        ! exp(-(l+1/2) alpha), term (l+1) Q_l^2 (1-e^-alpha)^2, geometric remainder.
        alpha = 2*asinh(1/(2*kappa))
        l = 1
        do
            q = sqrt(pi/(2*(l + 0.5_dp)*sinh(alpha)))*exp(-(l + 0.5_dp)*alpha)
            term = (l + 1)*(q*(1 - exp(-alpha)))**2
            tail = term/(1 - exp(-2*alpha))
            if ((tail < cfg%tolerance .and. eta*q < 1.0e-3_dp) .or. l >= cfg%lcap) exit
            l = l + 1
        end do
        plan%lmax = l
        plan%ledger%tail_estimate = tail; plan%ledger%lmax = l
        deg = cfg%order - 1
        h = mesh_spacing(1.0_dp, cfg%order, cfg%nodes_per_wavelength)
        rmax = turning_radius(plan%lmax, eta, kappa) + cfg%tail*kappa + 3*h
        nelem = ceiling(rmax/h)
        call build_radial_mesh(plan%mesh, nelem*h, cfg%order, h, h, 1.0_dp, ok)
        if (.not. ok) return
        plan%ledger%nodes = plan%mesh%n
        allocate (plan%inner(0:plan%lmax), plan%outer(0:plan%lmax))
        plan%ledger%barrier_weight = 0; plan%ledger%potential_at_match = 0
        do l = 0, plan%lmax
            rt = turning_radius(l, eta, kappa)
            j = deg*ceiling((rt + cfg%tail*kappa)/h)
            if (j + deg > plan%mesh%n) then
                ok = .false.; return
            end if
            plan%outer(l) = j
            plan%ledger%potential_at_match = max(plan%ledger%potential_at_match, &
                abs(veff(plan%mesh%r(j - deg), 0, eta, kappa))/0.5_dp)
            ! March inward from the turning point accumulating the WKB depth.
            jt = min(plan%mesh%n, max(1, count(plan%mesh%r < rt)))
            depth = 0; plan%inner(l) = 0
            do j = jt, 2, -1
                r1 = plan%mesh%r(j - 1); r2 = plan%mesh%r(j)
                k1 = sqrt(max(0.0_dp, 2*veff(r1, l, eta, kappa) - 1))
                k2 = sqrt(max(0.0_dp, 2*veff(r2, l, eta, kappa) - 1))
                depth = depth + (r2 - r1)*(k1 + k2)/2
                if (depth >= cfg%barrier .and. mod(j - 1, deg) == 0) then
                    plan%inner(l) = j - 1
                    plan%ledger%barrier_weight = max(plan%ledger%barrier_weight, &
                                                     exp(-2*depth))
                    exit
                end if
            end do
        end do
        ok = .true.
    end subroutine yukawa_plan

    !> Phase shifts delta_l (l = 0..lmax) of the repulsive Yukawa pair at
    !> Sommerfeld eta on the frozen plan, and d delta_l/d eta (equation-level
    !> tangent of the same banded system). Fixed energy E = 1/2 (k = 1):
    !> (H_l - E) a = 0 on rows inner+1..outer-1, a(inner) = 0, a(outer) = 1.
    subroutine yukawa_phases(plan, eta, delta, ddelta, ok)
        type(yukawa_plan_t), intent(inout) :: plan
        real(dp), intent(in) :: eta
        real(dp), allocatable, intent(out) :: delta(:), ddelta(:)
        logical, intent(out) :: ok
        real(dp), allocatable :: ab(:, :), a(:, :), jh(:), nh(:)
        integer, allocatable :: ipiv(:)
        real(dp) :: v, dv, t, cj, cn, dcj, dcn, res
        integer :: l, kd, i0, i1, m, i, j, info, ldab, deg, lmax
        ok = .false.
        lmax = plan%lmax; deg = plan%cfg%order - 1; kd = deg
        allocate (delta(0:lmax), ddelta(0:lmax), jh(0:lmax), nh(0:lmax))
        ldab = 3*kd + 1
        plan%ledger%fit_residual = 0
        do l = 0, lmax
            i0 = plan%inner(l); i1 = plan%outer(l)
            m = i1 - i0 - 1
            allocate (ab(ldab, m), a(m, 2), ipiv(m))
            ab = 0; a = 0
            do j = 1, m
                do i = max(1, j - kd), min(m, j + kd)
                    t = band(i0 + i, i0 + j)/2
                    if (i == j) then
                        call qc_screened_radial_potential(plan%mesh%r(i0 + i), real(l, dp), &
                            eta, 1.0_dp, 1.0_dp, plan%kappa, v)
                        t = t + v - 0.5_dp
                    end if
                    ab(2*kd + 1 + i - j, j) = t
                end do
            end do
            do i = max(1, m - kd + 1), m
                a(i, 1) = -band(i0 + i, i1)/2
            end do
            call dgbtrf(m, m, kd, kd, ab, ldab, ipiv, info)
            if (info /= 0) return
            call dgbtrs('N', m, kd, kd, 1, ab, ldab, ipiv, a(:, 1), m, info)
            if (info /= 0) return
            ! Implicit tangent: (H - E) da = -(dV/dg) a, boundary data fixed.
            do i = 1, m
                call qc_screened_radial_potential_jvp(plan%mesh%r(i0 + i), real(l, dp), &
                    eta, 1.0_dp, 1.0_dp, plan%kappa, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, &
                    0.0_dp, 0.0_dp, dv)
                a(i, 2) = -dv*a(i, 1)
            end do
            call dgbtrs('N', m, kd, kd, 1, ab, ldab, ipiv, a(:, 2), m, info)
            if (info /= 0) return
            call match(l, i0, i1, a, cj, cn, dcj, dcn, res)
            delta(l) = atan2(-cn, cj)
            delta(l) = delta(l) - pi*nint(delta(l)/pi)   ! modulo pi, |delta| <= pi/2
            ddelta(l) = (-cj*dcn + cn*dcj)/(cj**2 + cn**2)
            plan%ledger%fit_residual = max(plan%ledger%fit_residual, res)
            deallocate (ab, a, ipiv)
        end do
        ok = all(ieee_is_finite(delta)) .and. all(ieee_is_finite(ddelta))
    contains
        ! Entry (i,j) of the normalized DVR -d^2/dr^2 (symmetric lower band).
        real(dp) function band(ii, jj) result(x)
            integer, intent(in) :: ii, jj
            x = 0
            if (abs(ii - jj) > kd) return
            if (ii >= jj) then
                x = plan%mesh%tband(1 + ii - jj, jj)
            else
                x = plan%mesh%tband(1 + jj - ii, ii)
            end if
        end function band
        ! Least squares u = cj jhat + cn nhat over the last element's nodes
        ! (physical u = a/sqrt(w)); also for the tangent.
        subroutine match(l, i0, i1, a, cj, cn, dcj, dcn, res)
            integer, intent(in) :: l, i0, i1
            real(dp), intent(in) :: a(:, :)
            real(dp), intent(out) :: cj, cn, dcj, dcn, res
            real(dp) :: g11, g12, g22, r1, r2, s1, s2, det, u, du, x, fit, nrm
            integer :: idx
            g11 = 0; g12 = 0; g22 = 0; r1 = 0; r2 = 0; s1 = 0; s2 = 0; nrm = 0
            do idx = i1 - deg, i1
                x = plan%mesh%r(idx)
                call riccati_bessel(l, x, jh(0:l), nh(0:l))
                if (idx == i1) then
                    u = 1; du = 0
                else
                    u = a(idx - i0, 1); du = a(idx - i0, 2)
                end if
                u = u/sqrt(plan%mesh%w(idx)); du = du/sqrt(plan%mesh%w(idx))
                g11 = g11 + jh(l)**2; g12 = g12 + jh(l)*nh(l); g22 = g22 + nh(l)**2
                r1 = r1 + jh(l)*u; r2 = r2 + nh(l)*u
                s1 = s1 + jh(l)*du; s2 = s2 + nh(l)*du
                nrm = nrm + u**2
            end do
            det = g11*g22 - g12**2
            cj = (g22*r1 - g12*r2)/det; cn = (g11*r2 - g12*r1)/det
            dcj = (g22*s1 - g12*s2)/det; dcn = (g11*s2 - g12*s1)/det
            res = 0
            do idx = i1 - deg, i1
                x = plan%mesh%r(idx)
                call riccati_bessel(l, x, jh(0:l), nh(0:l))
                u = 1
                if (idx < i1) u = a(idx - i0, 1)
                u = u/sqrt(plan%mesh%w(idx))
                fit = cj*jh(l) + cn*nh(l)
                res = res + (u - fit)**2
            end do
            res = sqrt(res/nrm)
        end subroutine match
    end subroutine yukawa_phases

    !> Yukawa momentum-transfer log lnL = sum (l+1) sin^2(delta_l - delta_{l+1})/eta^2
    !> and d lnL/d eta (implicit phase tangents; plan frozen).
    subroutine yukawa_transfer_log(plan, eta, lnl, dlnl, ok, viscosity)
        type(yukawa_plan_t), intent(inout) :: plan
        real(dp), intent(in) :: eta
        real(dp), intent(out) :: lnl, dlnl
        logical, intent(out) :: ok
        real(dp), intent(out), optional :: viscosity
        real(dp), allocatable :: delta(:), ddelta(:)
        real(dp) :: d
        integer :: l
        call yukawa_phases(plan, eta, delta, ddelta, ok)
        lnl = 0; dlnl = 0
        if (.not. ok) return
        do l = plan%lmax - 1, 0, -1
            d = delta(l) - delta(l + 1)
            lnl = lnl + (l + 1)*sin(d)**2
            dlnl = dlnl + (l + 1)*sin(2*d)*(ddelta(l) - ddelta(l + 1))
        end do
        lnl = lnl/eta**2
        dlnl = dlnl/eta**2 - 2*lnl/eta
        if (present(viscosity)) viscosity = viscosity_from_phases(delta)/eta**2
    end subroutine yukawa_transfer_log
end module kin6d_quantum_collision_kernel
