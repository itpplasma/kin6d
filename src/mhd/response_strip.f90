module kin6d_response_strip
    !! R-B1: resolved resistive-inertial layer strip in the periodic cylinder
    !! (docs/response-unified.md, treatment (B); theory/derivations/
    !! response-r2-r4.md, (R2)-(R3) single fluid, no rotation, zero beta,
    !! tokamak ordering).
    !!
    !! Units: minor radius a = 1, B_z/R0 = 1, mass density rho = 1, time in
    !! tau_A = R0 sqrt(mu0 rho)/B_z, eta = 1/S.  One helicity (m, n),
    !! perturbation psi(r), phi(r) exp(i(m theta - n z/R0) + gamma t); gamma is
    !! complex (gamma = -i omega for a rotating drive in the lab frame).
    !!
    !! Strip |r - r_s| <= X_s, unknowns psi (flux) and phi (stream function):
    !!     gamma psi - i F phi - eta L psi = 0                         (Ohm, R3)
    !!     gamma L phi - i F L psi + i (m/r) J0' psi = 0          (vorticity, R2)
    !! L f = (1/r)(r f')' - m^2 f/r^2, F = m/q - n, J0 = mu0 J_z.  The ideal
    !! static limit gamma, eta -> 0 is the outer tearing equation L psi =
    !! Lambda psi, Lambda = m J0'/(r F) (cyl_response.py, Fitzpatrick 1993).
    !!
    !! Coupling.  Outside the strip the plasma is ideal and static.  Each
    !! interface carries the outer trace relation alpha psi' + beta psi = s
    !! (the Dirichlet-to-Neumann map of the outer region, i.e. the Schur
    !! complement of the outer system on that interface ring) and the ideal
    !! Ohm law i F phi = gamma psi.  No Delta', mu or layer coefficient is
    !! fitted or extrapolated; the layer is resolved on a radial mesh with a
    !! uniform core of a few layer widths and geometric grading to X_s.
    !!
    !! Equilibrium: the sharp current step (uniform J_z inside r1, none
    !! outside, cyl_response.py profile 'step' with w -> 0), whose outer
    !! traces are closed form; r_s Delta' = -2m/(1 + rho),
    !! rho = (r1/r_s)^(2m) g/(2m - g), g = 2/(1 - n q0/m) (m = 1: Delta' =
    !! infinity, the resistive internal kink).
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    private

    complex(dp), parameter :: iu = (0.0_dp, 1.0_dp)
    integer, parameter :: kl = 4, ku = 4

    type, public :: cyl_step_t
        !! Current step: q = q0 inside r1, q0 (r/r1)^2 outside, plasma r <= 1,
        !! vacuum beyond (no wall), external drive psi_ext r^m.
        integer :: m = 2, n = 1
        real(dp) :: q0 = 1.2_dp, r1 = 0.5_dp
        complex(dp) :: psi_ext = (1.0_dp, 0.0_dp)
    end type cyl_step_t

    type, public :: strip_options_t
        real(dp) :: xs = 0.08_dp       !! strip half width X_s
        real(dp) :: core = 4.0_dp      !! uniform core half width / layer width
        real(dp) :: hcore = 0.05_dp    !! core spacing / layer width
        real(dp) :: ratio = 1.05_dp    !! geometric grading ratio beyond the core
        real(dp) :: hmax = 1.0e-3_dp   !! spacing cap (outer ring spacing)
    end type strip_options_t

    type, public :: trace_t
        !! Outer trace relation alpha psi' + beta psi = s on one interface.
        complex(dp) :: alpha = (1.0_dp, 0.0_dp), beta = (0.0_dp, 0.0_dp)
        complex(dp) :: s = (0.0_dp, 0.0_dp)
    end type trace_t

    type, public :: strip_solution_t
        real(dp), allocatable :: r(:)
        complex(dp), allocatable :: psi(:), phi(:), j(:)
        complex(dp) :: psi_s = (0.0_dp, 0.0_dp)   !! psi at r_s
        complex(dp) :: current = (0.0_dp, 0.0_dp) !! int r j dr over the strip
        real(dp) :: width = 0.0_dp                 !! HWHM of |j| about r_s
        integer :: info = 0
    end type strip_solution_t

    public :: step_rational, step_fprime, step_delta_prime, step_q, step_dj
    public :: step_traces, layer_width, strip_mesh, strip_solve
    public :: strip_growth_rate, strip_dirichlet_map
    public :: ri_layer_dhat, ri_dispersion_gamma

contains

    pure real(dp) function step_rational(eq) result(rs)
        type(cyl_step_t), intent(in) :: eq
        rs = eq%r1*sqrt(real(eq%m, dp)/(eq%n*eq%q0))
    end function step_rational

    pure real(dp) function step_q(eq, r) result(q)
        type(cyl_step_t), intent(in) :: eq
        real(dp), intent(in) :: r
        q = eq%q0*max(1.0_dp, (r/eq%r1)**2)
    end function step_q

    pure real(dp) function step_dj(eq, r) result(dj)
        !! J0' off the step (the step itself lies outside every strip).
        type(cyl_step_t), intent(in) :: eq
        real(dp), intent(in) :: r
        dj = 0.0_dp*(eq%q0 + r)
    end function step_dj

    pure real(dp) function step_fprime(eq) result(fp)
        !! |dF/dr| at r_s; F = m/q - n.
        type(cyl_step_t), intent(in) :: eq
        fp = 2.0_dp*eq%n/step_rational(eq)
    end function step_fprime

    pure real(dp) function step_delta_prime(eq) result(dp_rs)
        !! r_s Delta' of the sharp step (huge for the marginal m = 1).
        type(cyl_step_t), intent(in) :: eq
        real(dp) :: g, rho, rs
        rs = step_rational(eq)
        g = 2.0_dp/(1.0_dp - eq%n*eq%q0/eq%m)
        rho = (eq%r1/rs)**(2*eq%m)*g/(2*eq%m - g)
        if (abs(1.0_dp + rho) < 1.0e-12_dp) then
            dp_rs = huge(1.0_dp)
        else
            dp_rs = -2.0_dp*eq%m/(1.0_dp + rho)
        end if
    end function step_delta_prime

    subroutine step_traces(eq, rm, rp, tin, tout)
        !! Exact outer traces of the step at the strip edges rm < r_s < rp.
        !! Inside: psi = r^m (r < r1), [psi'] = -g psi/r1 at r1, A r^m + B r^-m
        !! beyond; trace psi_in psi' - psi_in' psi = 0.  Outside: vacuum-like
        !! A r^m + B r^-m with psi'(1) = -m psi(1) + 2m psi_ext, i.e.
        !! psi' + (m/r) psi = 2m psi_ext r^(m-1).
        type(cyl_step_t), intent(in) :: eq
        real(dp), intent(in) :: rm, rp
        type(trace_t), intent(out) :: tin, tout
        real(dp) :: g, a, b, pin, dpin
        integer :: m
        m = eq%m
        g = 2.0_dp/(1.0_dp - eq%n*eq%q0/m)
        a = 1.0_dp - g/(2*m)
        b = g*eq%r1**(2*m)/(2*m)
        pin = a*rm**m + b*rm**(-m)
        dpin = m*(a*rm**(m - 1) - b*rm**(-m - 1))
        tin%alpha = pin
        tin%beta = -dpin
        tin%s = 0.0_dp
        tout%alpha = 1.0_dp
        tout%beta = real(m, dp)/rp
        tout%s = 2*m*eq%psi_ext*rp**(m - 1)
    end subroutine step_traces

    pure real(dp) function layer_width(eta, gamma, fp) result(w)
        !! Mesh scale: the ideal-inertial width |gamma|/F' or, with
        !! resistivity, the larger of it and min(constant-psi RI width
        !! (|gamma| eta/F'^2)^(1/4), non-constant-psi width (eta^2/F'^2)^(1/6)).
        real(dp), intent(in) :: eta, fp
        complex(dp), intent(in) :: gamma
        w = abs(gamma)/fp
        if (eta > 0.0_dp) w = max(w, min((abs(gamma)*eta/fp**2)**0.25_dp, &
            (eta**2/fp**2)**(1.0_dp/6.0_dp)))
    end function layer_width

    subroutine strip_mesh(rs, w, opt, r)
        !! Symmetric mesh on [rs - X_s, rs + X_s] with a node at rs: uniform
        !! spacing hcore w for |x| <= core w, then geometric growth by ratio up
        !! to the spacing cap hmax.
        real(dp), intent(in) :: rs, w
        type(strip_options_t), intent(in) :: opt
        real(dp), allocatable, intent(out) :: r(:)
        real(dp), allocatable :: x(:)
        real(dp) :: h, xc
        integer :: k, nside
        h = opt%hcore*w
        allocate (x(0))
        xc = 0.0_dp
        do
            if (xc + h >= opt%xs) exit
            if (xc + h > opt%core*w) h = min(h*opt%ratio, max(opt%hmax, opt%hcore*w))
            xc = xc + h
            x = [x, xc]
        end do
        ! The last interval would be shorter than h/2: merge it into the edge.
        if (size(x) > 0) then
            if (opt%xs - x(size(x)) < 0.5_dp*h) x = x(1:size(x) - 1)
        end if
        x = [x, opt%xs]
        nside = size(x)
        allocate (r(2*nside + 1))
        r(nside + 1) = rs
        do k = 1, nside
            r(nside + 1 + k) = rs + x(k)
            r(nside + 1 - k) = rs - x(k)
        end do
    end subroutine strip_mesh

    subroutine strip_solve(eq, eta, gamma, r, sol, traces)
        !! Driven strip response at complex gamma on the mesh r (from
        !! strip_mesh, node of r_s in the middle), coupled to the exact outer
        !! traces of eq, or to the given traces(1:2) (inner, outer) of another
        !! outer region (R-B2: the FEEC torus).  Banded LU (LAPACK zgbsv).
        type(cyl_step_t), intent(in) :: eq
        real(dp), intent(in) :: eta
        complex(dp), intent(in) :: gamma
        real(dp), intent(in) :: r(:)
        type(strip_solution_t), intent(out) :: sol
        type(trace_t), intent(in), optional :: traces(2)
        complex(dp), allocatable :: ab(:, :), rhs(:, :)
        integer, allocatable :: ipiv(:)
        type(trace_t) :: tin, tout
        complex(dp) :: lc(3), d1(3), fnode
        real(dp) :: hm, hp
        integer :: nn, ntot, i, row, ldab, info, mid
        external :: zgbsv

        nn = size(r)
        ntot = 2*nn
        ldab = 2*kl + ku + 1
        allocate (ab(ldab, ntot), rhs(ntot, 1), ipiv(ntot))
        ab = (0.0_dp, 0.0_dp)
        rhs = (0.0_dp, 0.0_dp)
        if (present(traces)) then
            tin = traces(1)
            tout = traces(2)
        else
            call step_traces(eq, r(1), r(nn), tin, tout)
        end if

        do i = 2, nn - 1
            call lap_stencil(eq%m, r(i - 1:i + 1), lc)
            fnode = eq%m/step_q(eq, r(i)) - eq%n
            row = 2*i - 1
            ! Ohm: gamma psi_i - i F phi_i - eta (L psi)_i = 0
            call put(row, 2*i - 3, -eta*lc(1))
            call put(row, 2*i - 1, gamma - eta*lc(2))
            call put(row, 2*i + 1, -eta*lc(3))
            call put(row, 2*i, -iu*fnode)
            row = 2*i
            ! vorticity: gamma (L phi)_i - i F (L psi)_i + i (m/r) J0' psi_i = 0
            call put(row, 2*i - 2, gamma*lc(1))
            call put(row, 2*i, gamma*lc(2))
            call put(row, 2*i + 2, gamma*lc(3))
            call put(row, 2*i - 3, -iu*fnode*lc(1))
            call put(row, 2*i - 1, -iu*fnode*lc(2) + iu*eq%m/r(i)*step_dj(eq, r(i)))
            call put(row, 2*i + 1, -iu*fnode*lc(3))
        end do

        ! inner interface: trace relation with a one-sided second-order psi'
        hm = r(2) - r(1)
        hp = r(3) - r(2)
        d1 = [-(2*hm + hp)/(hm*(hm + hp)), (hm + hp)/(hm*hp), -hm/(hp*(hm + hp))]
        call put(1, 1, tin%alpha*d1(1) + tin%beta)
        call put(1, 3, tin%alpha*d1(2))
        call put(1, 5, tin%alpha*d1(3))
        rhs(1, 1) = tin%s
        fnode = eq%m/step_q(eq, r(1)) - eq%n
        call put(2, 1, -gamma)
        call put(2, 2, iu*fnode)
        ! outer interface
        hm = r(nn - 1) - r(nn - 2)
        hp = r(nn) - r(nn - 1)
        d1 = [hp/(hm*(hm + hp)), -(hm + hp)/(hm*hp), (2*hp + hm)/(hp*(hm + hp))]
        row = 2*nn - 1
        call put(row, 2*nn - 5, tout%alpha*d1(1))
        call put(row, 2*nn - 3, tout%alpha*d1(2))
        call put(row, 2*nn - 1, tout%alpha*d1(3) + tout%beta)
        rhs(row, 1) = tout%s
        fnode = eq%m/step_q(eq, r(nn)) - eq%n
        call put(2*nn, 2*nn - 1, -gamma)
        call put(2*nn, 2*nn, iu*fnode)

        call zgbsv(ntot, kl, ku, 1, ab, ldab, ipiv, rhs, ntot, info)
        sol%info = info
        sol%r = r
        sol%psi = rhs(1:ntot:2, 1)
        sol%phi = rhs(2:ntot:2, 1)
        if (info /= 0) return
        mid = (nn + 1)/2
        sol%psi_s = sol%psi(mid)
        allocate (sol%j(nn))
        do i = 2, nn - 1
            call lap_stencil(eq%m, r(i - 1:i + 1), lc)
            sol%j(i) = sum(lc*sol%psi(i - 1:i + 1))
        end do
        sol%j(1) = sol%j(2)
        sol%j(nn) = sol%j(nn - 1)
        ! layer current: int r j dr = [r psi'] - int m^2 psi/r dr, in the
        ! conservative form of the stencil (exact telescoping sum)
        sol%current = (0.0_dp, 0.0_dp)
        do i = 2, nn - 1
            hm = r(i) - r(i - 1)
            hp = r(i + 1) - r(i)
            sol%current = sol%current + r(i)*sol%j(i)*0.5_dp*(hm + hp)
        end do
        sol%width = hwhm(r, sol%j, mid)

    contains

        subroutine put(i_row, j_col, v)
            integer, intent(in) :: i_row, j_col
            complex(dp), intent(in) :: v
            ab(kl + ku + 1 + i_row - j_col, j_col) = &
                ab(kl + ku + 1 + i_row - j_col, j_col) + v
        end subroutine put

    end subroutine strip_solve

    subroutine strip_dirichlet_map(eq, eta, gamma, r, dmap, info)
        !! Dirichlet-to-Neumann map of the strip on the mesh r: the edge
        !! slopes psi'(r(1)), psi'(r(n)) = dmap psi(r(1)), psi(r(n)) (one-sided
        !! second-order differences, as in the trace rows), with the ideal Ohm
        !! law on both edges.  This is the strip's Schur complement on its two
        !! interface rings (R-B2 coupling to the FEEC outer region).
        type(cyl_step_t), intent(in) :: eq
        real(dp), intent(in) :: eta
        complex(dp), intent(in) :: gamma
        real(dp), intent(in) :: r(:)
        complex(dp), intent(out) :: dmap(2, 2)
        integer, intent(out) :: info
        type(trace_t) :: tr(2)
        type(strip_solution_t) :: sol
        real(dp) :: hm, hp, d1(3), d2(3)
        integer :: k, nn

        nn = size(r)
        hm = r(2) - r(1)
        hp = r(3) - r(2)
        d1 = [-(2*hm + hp)/(hm*(hm + hp)), (hm + hp)/(hm*hp), -hm/(hp*(hm + hp))]
        hm = r(nn - 1) - r(nn - 2)
        hp = r(nn) - r(nn - 1)
        d2 = [hp/(hm*(hm + hp)), -(hm + hp)/(hm*hp), (2*hp + hm)/(hp*(hm + hp))]
        info = 0
        do k = 1, 2
            tr%alpha = (0.0_dp, 0.0_dp)
            tr%beta = (1.0_dp, 0.0_dp)
            tr%s = (0.0_dp, 0.0_dp)
            tr(k)%s = (1.0_dp, 0.0_dp)
            call strip_solve(eq, eta, gamma, r, sol, tr)
            if (sol%info /= 0) then
                info = sol%info
                return
            end if
            dmap(1, k) = sum(d1*sol%psi(1:3))
            dmap(2, k) = sum(d2*sol%psi(nn - 2:nn))
        end do
    end subroutine strip_dirichlet_map

    pure subroutine lap_stencil(m, rr, lc)
        !! Conservative three-point (1/r)(r f')' - m^2 f/r^2 on a nonuniform mesh.
        integer, intent(in) :: m
        real(dp), intent(in) :: rr(3)
        complex(dp), intent(out) :: lc(3)
        real(dp) :: hm, hp, rmh, rph, c
        hm = rr(2) - rr(1)
        hp = rr(3) - rr(2)
        rmh = 0.5_dp*(rr(1) + rr(2))
        rph = 0.5_dp*(rr(2) + rr(3))
        c = 2.0_dp/((hm + hp)*rr(2))
        lc(1) = c*rmh/hm
        lc(3) = c*rph/hp
        lc(2) = -lc(1) - lc(3) - real(m*m, dp)/rr(2)**2
    end subroutine lap_stencil

    pure real(dp) function hwhm(r, j, mid) result(w)
        !! Half width at half maximum of |j| about the node mid (mean of the
        !! two sides, linear interpolation).
        real(dp), intent(in) :: r(:)
        complex(dp), intent(in) :: j(:)
        integer, intent(in) :: mid
        real(dp) :: half, side(2), a0, a1
        integer :: k, dir, s
        half = 0.5_dp*abs(j(mid))
        side = r(size(r)) - r(1)
        do s = 1, 2
            dir = merge(1, -1, s == 1)
            k = mid
            do
                if (k + dir < 1 .or. k + dir > size(r)) exit
                a0 = abs(j(k))
                a1 = abs(j(k + dir))
                if (a1 <= half) then
                    side(s) = abs(r(k) - r(mid)) + abs(r(k + dir) - r(k)) &
                        *(a0 - half)/max(a0 - a1, tiny(1.0_dp))
                    exit
                end if
                k = k + dir
            end do
        end do
        w = 0.5_dp*sum(side)
    end function hwhm

    subroutine strip_growth_rate(eq, eta, guess, opt, gamma, sol, iterations, status)
        !! Tearing eigenvalue of the coupled problem: zero of psi_ext/psi_s(gamma)
        !! (pole of the driven response), complex secant from guess and
        !! 1.02 guess on the mesh of the guess.  status 0 on convergence
        !! (relative step < 1e-11).
        type(cyl_step_t), intent(in) :: eq
        real(dp), intent(in) :: eta
        complex(dp), intent(in) :: guess
        type(strip_options_t), intent(in) :: opt
        complex(dp), intent(out) :: gamma
        type(strip_solution_t), intent(out) :: sol
        integer, intent(out) :: iterations, status
        real(dp), allocatable :: r(:)
        complex(dp) :: g0, g1, f0, f1, g2
        real(dp) :: rs, fp
        integer :: it
        rs = step_rational(eq)
        fp = step_fprime(eq)
        call strip_mesh(rs, layer_width(eta, guess, fp), opt, r)
        g0 = guess
        g1 = 1.02_dp*guess
        f0 = response_inverse(g0)
        f1 = response_inverse(g1)
        status = 1
        do it = 1, 60
            iterations = it
            g2 = g1 - f1*(g1 - g0)/(f1 - f0)
            g0 = g1
            f0 = f1
            g1 = g2
            f1 = response_inverse(g1)
            if (abs(g1 - g0) <= 1.0e-11_dp*abs(g1)) then
                status = 0
                exit
            end if
        end do
        gamma = g1

    contains

        complex(dp) function response_inverse(g) result(f)
            complex(dp), intent(in) :: g
            call strip_solve(eq, eta, g, r, sol)
            f = eq%psi_ext/sol%psi_s
        end function response_inverse

    end subroutine strip_growth_rate

    pure real(dp) function ri_layer_dhat(q) result(d)
        !! Oracle, not used by the solver: resistive-inertial slab layer
        !! Delta' d = -(pi/8) Q^(5/4) Gamma((Q^(3/2) - 1)/4)/Gamma((Q^(3/2) + 5)/4)
        !! with d = (rho eta^2/F'^2)^(1/6), Q = gamma d^2/eta (Ara et al. 1978;
        !! constants fixed by the constant-psi limit 2 pi Gamma(3/4)/Gamma(1/4)
        !! Q^(5/4) of FKR and by an independent BVP solve of the scaled layer
        !! to 1e-4 at Q = 0.01-2).  Q = 1 is the Delta' = infinity root.
        real(dp), intent(in) :: q
        real(dp) :: z
        z = q**1.5_dp
        d = -0.125_dp*acos(-1.0_dp)*q**1.25_dp*gamma(0.25_dp*(z - 1.0_dp)) &
            /gamma(0.25_dp*(z + 5.0_dp))
    end function ri_layer_dhat

    subroutine ri_dispersion_gamma(eta, fp, dprime, gamma_ri, status)
        !! Oracle growth rate from Delta' d = ri_layer_dhat(Q) (rho = 1);
        !! dprime >= huge/2 means Delta' = infinity (Q = 1).
        use fortnum_roots, only: root_brent
        use fortnum_status, only: fortnum_status_t
        real(dp), intent(in) :: eta, fp, dprime
        real(dp), intent(out) :: gamma_ri
        integer, intent(out) :: status
        type(fortnum_status_t) :: st
        real(dp) :: d, target, q
        d = (eta**2/fp**2)**(1.0_dp/6.0_dp)
        status = 0
        if (dprime >= 0.5_dp*huge(1.0_dp)) then
            q = 1.0_dp
        else
            target = dprime*d
            call root_brent(residual, 1.0e-14_dp, 1.0_dp - 1.0e-12_dp, q, st, &
                xtol=1.0e-15_dp)
            if (st%code /= 0) status = 1
        end if
        gamma_ri = q*eta/d**2

    contains

        pure real(dp) function residual(x) result(y)
            real(dp), intent(in) :: x
            y = ri_layer_dhat(x) - target
        end function residual

    end subroutine ri_dispersion_gamma

end module kin6d_response_strip
