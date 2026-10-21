program test_mhd_closure
    !! V5: ideal closure at rational surfaces (kin6d_mhd_closure_n) against
    !! cylindrical tearing Delta'.
    !!
    !! A large-aspect-ratio torus (circular surfaces about (R0, 0), minor
    !! radius 1, zero beta) carries the current profile mu0 J_z = (2 B0/(R0
    !! q0)) j(r).  Cylindrical outer region (k -> 0):
    !!     psi'' + psi'/r - m^2 psi/r^2 - G psi = 0,
    !!     G = mu0 J_z'/(B_theta (1 - n q/m)) = j' r/(I(r) (1 - n q/m)),
    !! I(r) = int_0^r j r dr, q = q0 r^2/(2 I), no wall (psi ~ r^-m outside
    !! r = 1, where j = 0).  Near the rational, psi = A (1 + K x ln|x|) + B_+-
    !! x and Delta' = (B_+ - B_-)/A.
    !!
    !! Oracles:
    !!  1. Step current j = 1 (r < r1), 0 (r > r1): closed form
    !!         r_s Delta' = -2m (m - 1 - n q0)/(m - 1 - n q0 + (n q0/m)^m)
    !!     checks the independent ODE solver (fortnum Cash-Karp, adaptive,
    !!     Richardson in the matching distance) on tanh-smoothed steps.
    !!  2. The FEEC response with the Frobenius closure in Delta' mode
    !!     (unit large-solution coefficients on both sides of q = m/n, the
    !!     other rationals shielded, no external field) on the torus against
    !!     the ODE: Wesson j = (1 - r^2)^2 (resonant current gradient, x ln|x|
    !!     large solution) at R0/a = 10 and 40 with the layer remainder
    !!     extrapolated (large-aspect-ratio limit), and a tanh step of width
    !!     0.08 at R0/a = 40.
    !!  3. Edge current jump (uniform current to the plasma edge, vacuum
    !!     outside, q0 = 1.2, m = 3, n = 1 drive): plasma response 1 +
    !!     1/(m - n q0 - 1) in the cylinder limit (field-aligned boundary
    !!     sheet = moved current xi_n J0 at zero beta).
    !!  4. Mercier index: D_I = -1/4 exactly at zero pressure gradient.
    !!  4b. Driven ideal shielding of q = 3/1 in the GPEC lar default cylinder
    !!     (q = 2.2 (1 + r^2), edge current j(1) /= 0) against the cylinder
    !!     ODE with the edge-current jump (cyl_response.py S_edge_current =
    !!     2.2967, edge total/applied 0.6507) at ring ratios 2 and 3.
    !! Argument "resistive" runs the L5 suite instead (kin6d_mhd_layer_n):
    !!  5. Scaled layer solver: constant-psi resistive-inertial limit 2 pi
    !!     Gamma(3/4)/Gamma(1/4) g^(5/4)/eps, resolution, g -> conj(g),
    !!     viscous-resistive scaling g P^(1/6)/eps.
    !!  6. Cylinder (Wesson profile, q = 2/1, static drive r^m in a plasma
    !!     rotating with gamma = i n Omega): the matched response (ODE Delta'
    !!     and drive jump E, numerical layer Delta and b(0)/b0) against the
    !!     full linear resistive reduced-MHD BVP (no matching), constant-psi
    !!     and inertial (Q up to 12) layers; Omega -> -Omega conjugates.
    !!  7. FEEC resistive closure on the torus R0/a = 40 against the BVP.
    !!  8. Ideal limit: eta -> 0 at fixed rotation gives the ideal-inertial
    !!     layer, whose current tends to the L6 closure as the inertial width
    !!     |gamma| (mu0 rho)^(1/2)/F' -> 0 (order 1/2 in rho).
    !! KIN6D_CLOSURE_STUDY="prof R0 nrad delta ratio maxit nfit" runs one FEEC
    !! Delta' case (convergence studies), "edge" the edge oracle alone.
    !$ use omp_lib, only: omp_get_max_threads, omp_set_num_threads
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_ode_cash_karp, only: cash_karp_step
    use fortfem_meridian_feec_n, only: meridian_feec_t, meridian_feec_build, &
        meridian_point_basis_t, meridian_point_basis_init
    use kin6d_vacuum_surfaces, only: psi_grid_t, psi_grid_init, flux_ray_point, &
        surface_t, surface_build, find_surface_q
    use kin6d_mhd_mesh, only: flux_mesh_t, flux_mesh_build, ring_labels_merge, &
        ring_labels_smooth
    use kin6d_gs_direct, only: gs_direct_available
    use kin6d_mhd_ampere_n, only: ampere_n_t, ampere_n_build, ampere_n_solve, &
        ampere_n_curl, mu0
    use kin6d_mhd_response_n, only: ring_sheet_load, layer_closure_t, closure_fluxes, &
        ring_angles_t, ring_angles_set, ring_flux_harmonic, ring_flux_from_values, &
        ring_points
    use kin6d_mhd_ideal_n, only: sample_ring_t, sample_ring_build, &
        ideal_options_t, ideal_response
    use kin6d_mhd_closure_n, only: rational_layer_t, mercier_di, layer_labels, &
        layer_xmax, layer_fit_widths, layer_jump, layer_kink, layer_init, &
        layer_closure, layer_wmax, &
        layer_coefficients, layer_regular_slope, layer_large
    use kin6d_mhd_layer_n, only: layer_params_t, layer_scales, layer_delta, &
        layer_scaled_delta, layer_ri_constant_psi, layer_metrics, resistive_closure
    use kin6d_response_strip, only: cyl_step_t, strip_options_t, strip_mesh, &
        layer_width, step_fprime, step_delta_prime, step_rational, &
        strip_dirichlet_map, strip_growth_rate, strip_solution_t, &
        ri_dispersion_gamma
    implicit none

    real(dp), parameter :: pi = 3.14159265358979323846_dp
    real(dp), parameter :: b0 = 2.0_dp
    integer, parameter :: ntab = 14001
    real(dp), parameter :: rtab_max = 1.4_dp
    integer :: prof = 1                  ! 1 Wesson, 2 tanh step, 3 uniform,
                                         ! 4 GPEC lar default j = 1/(1 + r^2)^2
    real(dp) :: q0 = 1.2_dp, r1 = 0.5_dp, wstep = 0.08_dp
    integer :: nmode = 1, mres = 2
    real(dp) :: r0, msgn = 1.0_dp
    real(dp) :: tab_r(ntab), tab_i(ntab), tab_p(ntab), tab_b(ntab)
    type(psi_grid_t) :: g
    logical :: all_passed
    integer :: nrad_fe = 32, gmres_max = 400, nfit_fe = 3
    real(dp) :: shield_tol = 1.0e-6_dp
    character(len=32) :: suite
    ! L5 resistive layer (SI): mass density, magnetic diffusivity eta/mu0,
    ! Doppler rate gamma = i n Omega of the static drive
    real(dp) :: rho_l = 1.0e-7_dp, diff_l = 8.8e-5_dp
    complex(dp) :: gam_l = (0.0_dp, 300.0_dp)
    ! R-B2 strip coupling state (fe_strip and its helpers)
    integer, parameter :: np = 256
    real(dp), parameter :: ratio = 1.5_dp
    type(flux_mesh_t) :: fm
    type(meridian_feec_t) :: fs
    type(meridian_point_basis_t) :: pb
    type(ampere_n_t) :: amp
    type(ideal_options_t) :: opt
    type(sample_ring_t), allocatable :: sr(:)
    type(layer_closure_t) :: cl
    type(surface_t) :: s
    type(ring_angles_t), allocatable :: ras(:)
    type(cyl_step_t) :: eq
    type(strip_options_t) :: sopt
    type(strip_solution_t) :: ssol
    real(dp), allocatable :: extra(:), ul(:), psin_s(:), rr(:), xq(:, :), nds(:, :)
    real(dp), allocatable :: th(:), tau(:), wq(:)
    integer, allocatable :: tt(:)
    integer, allocatable :: mr(:), je_all(:)
    complex(dp), allocatable :: b2(:), f(:), a(:), b2cl(:, :), lam(:), phit(:)
    complex(dp), allocatable :: phe(:), bc(:)
    complex(dp) :: mo(2, 2), mfe(2, 2), p0(2), dmap(2, 2), vmap(2, 2), smap(2, 2), jr(2)
    complex(dp) :: current_map(2)
    complex(dp) :: lam2(2), lam_l4, ga, gb, fa, fb, gn, bb(3), amat(2, 2)
    real(dp) :: merr, qsgn, qt, x, pa, rs, re(2), fp, dpr, gr, rb(2, 2), wl, tol
    integer :: k, st, nrat, itest, j, je(2), side, nr, ns, info, its, it
    logical :: found

    all_passed = .true.
    call study()
    block
        character(len=16) :: mode
        call get_command_argument(1, mode)
        if (trim(mode) == 'strip') then
            !$ call omp_set_num_threads(min(4, omp_get_max_threads()))
            call strip_suite()
        end if
        if (trim(mode) == 'setup') then
            call layer_all_rational_setup()
            if (.not. all_passed) error stop 1
            stop
        end if
    end block
    !$ call omp_set_num_threads(min(4, omp_get_max_threads()))
    if (command_argument_count() > 0) then
        call get_command_argument(1, suite)
        if (trim(suite) == 'resistive') then
            call resistive_suite()
            if (.not. all_passed) error stop 1
            print '(a)', "test_mhd_closure resistive: all checks passed"
            stop
        end if
    end if
    call step_closed_form()
    call mercier_zero_beta()
    call layer_rings_fine_mesh()
    call edge_current_jump()
    prof = 1
    call layer_conditioning()
    call layer_ring_cap()
    call layer_jump_exact()
    call layer_kink_finite_mu()
    call layer_closure_finite_mu()
    call shield_cylinder()
    prof = 1
    call fe_ladder("Wesson j = (1 - r^2)^2", [10.0_dp, 40.0_dp], 2.0e-2_dp)
    prof = 2
    wstep = 0.08_dp
    call fe_single("tanh step, width 0.08", 40.0_dp, 1.0e-3_dp, 3.0e-2_dp)
    if (.not. all_passed) error stop 1
    print '(a)', "test_mhd_closure: all checks passed"

contains

    subroutine layer_all_rational_setup()
        !! Unsorted, crowded and near-boundary rationals must all retain
        !! two complete fit stencils without crossing a neighbor or boundary.
        real(dp), parameter :: centers(4) = &
            [0.9996_dp, 0.49_dp, 0.4902_dp, 0.0001_dp]
        real(dp), parameter :: lo = 1.0e-6_dp, hi = 0.9998_dp
        real(dp) :: widths(4), reach, separation
        real(dp), allocatable :: extra(:), labels(:), grid(:)
        type(rational_layer_t) :: lay
        integer :: i, j, st

        call layer_fit_widths(centers, 1.0e-3_dp, 2.0_dp, 3, &
            lo, hi, widths, st)
        call record(st == 0, 'all-rational widths accept crowded and edge layers')
        if (st /= 0) return
        extra = sqrt(centers)
        do i = 1, size(centers)
            reach = 4*widths(i)
            separation = min(centers(i) - lo, hi - centers(i))
            do j = 1, size(centers)
                if (j == i) cycle
                separation = min(separation, abs(centers(i) - centers(j)))
            end do
            call record(reach < 0.5_dp*separation, &
                'fit stencil remains inside its rational neighborhood')
            call layer_labels(centers(i), widths(i), 2.0_dp, &
                reach*(1 + 1.0e-9_dp), lo, hi, labels)
            call record(size(labels) == 6, 'both complete fit stencils retained')
            extra = [extra, sqrt(labels)]
        end do
        call ring_labels_merge(128, extra, 0.3_dp, grid)
        do i = 1, size(centers)
            call layer_init(grid, centers(i), 2.0_dp, 1.0_dp, -0.28_dp, &
                i, widths(i), 2.0_dp, 3, lay, st)
            call record(st == 0, 'capped fit rings initialize on the merged mesh')
        end do
        call layer_fit_widths([0.5_dp, 0.5_dp], 1.0e-3_dp, 2.0_dp, 3, &
            lo, hi, widths(:2), st)
        call record(st /= 0, 'duplicate rational centers refused')
        call layer_fit_widths([hi], 1.0e-3_dp, 2.0_dp, 3, &
            lo, hi, widths(:1), st)
        call record(st /= 0, 'center outside the fit domain refused')
        call layer_fit_widths([0.5_dp], 1.0e-3_dp, 1.0_dp, 3, &
            lo, hi, widths(:1), st)
        call record(st /= 0, 'non-growing fit grading refused')
    end subroutine layer_all_rational_setup

    subroutine record(c, d)
        logical, intent(in) :: c
        character(len=*), intent(in) :: d
        all_passed = all_passed .and. c
        if (c) then
            print '(a,a)', "  PASS ", d
        else
            print '(a,a)', "  FAIL ", d
        end if
    end subroutine record

    ! ------------------------------------------------------------ profiles

    pure real(dp) function jprof(r)
        real(dp), intent(in) :: r
        if (prof == 3) then
            jprof = 1.0_dp
        else if (prof == 4) then
            jprof = 1.0_dp/(1.0_dp + r*r)**2
        else if (prof == 1) then
            jprof = max(0.0_dp, 1.0_dp - r*r)**2
        else
            jprof = (1.0_dp - tanh((r - r1)/wstep))/(1.0_dp + tanh(r1/wstep))
        end if
    end function jprof

    pure real(dp) function djprof(r)
        real(dp), intent(in) :: r
        if (prof == 3) then
            djprof = 0.0_dp
        else if (prof == 4) then
            djprof = -4.0_dp*r/(1.0_dp + r*r)**3
        else if (prof == 1) then
            djprof = -4.0_dp*r*max(0.0_dp, 1.0_dp - r*r)
        else
            djprof = -1.0_dp/(wstep*cosh((r - r1)/wstep)**2*(1.0_dp + tanh(r1/wstep)))
        end if
    end function djprof

    subroutine profile_rhs(t, y, dydt, ctx)
        !! y = (I, int I/r, int I j/r)
        real(dp), intent(in) :: t
        real(dp), intent(in) :: y(:)
        real(dp), intent(out) :: dydt(:)
        class(*), intent(in), optional :: ctx
        dydt(1) = jprof(t)*t
        if (t > 0.0_dp) then
            dydt(2) = y(1)/t
            dydt(3) = y(1)*jprof(t)/t
        else
            dydt(2:3) = 0.0_dp
        end if
    end subroutine profile_rhs

    subroutine make_tables()
        !! I(r), psi-hat(r) = int I/r and int I j/r on a uniform grid (RK5).
        real(dp) :: y(3), y5(3), yerr(3), k(3, 6), yt(3), h, t
        integer :: i, nfev
        h = rtab_max/(ntab - 1)
        y = 0.0_dp
        t = 0.0_dp
        nfev = 0
        tab_r(1) = 0.0_dp
        tab_i(1) = 0.0_dp
        tab_p(1) = 0.0_dp
        tab_b(1) = 0.0_dp
        do i = 2, ntab
            call cash_karp_step(profile_rhs, t, y, h, .false., k(:, 1), k(:, 2), &
                k(:, 3), k(:, 4), k(:, 5), k(:, 6), yt, y5, yerr, nfev)
            y = y5
            t = (i - 1)*h
            tab_r(i) = t
            tab_i(i) = y(1)
            tab_p(i) = y(2)
            tab_b(i) = y(3)
        end do
    end subroutine make_tables

    real(dp) function hermite(tab, dfun, r)
        !! Cubic Hermite interpolation of a table with known derivative.
        real(dp), intent(in) :: tab(:), r
        interface
            real(dp) function dfun(x)
                import :: dp
                real(dp), intent(in) :: x
            end function dfun
        end interface
        real(dp) :: h, s, ra, rb, h00, h10, h01, h11
        integer :: i
        h = tab_r(2)
        i = min(ntab - 1, max(1, int(r/h) + 1))
        ra = tab_r(i)
        rb = tab_r(i + 1)
        s = (r - ra)/h
        h00 = (1 + 2*s)*(1 - s)**2
        h10 = s*(1 - s)**2
        h01 = s*s*(3 - 2*s)
        h11 = s*s*(s - 1)
        hermite = h00*tab(i) + h10*h*dfun(ra) + h01*tab(i + 1) + h11*h*dfun(rb)
    end function hermite

    real(dp) function d_i(r)
        real(dp), intent(in) :: r
        d_i = jprof(r)*r
    end function d_i

    real(dp) function d_p(r)
        real(dp), intent(in) :: r
        if (r > 0.0_dp) then
            d_p = hermite(tab_i, d_i, r)/r
        else
            d_p = 0.0_dp
        end if
    end function d_p

    real(dp) function d_b(r)
        real(dp), intent(in) :: r
        if (r > 0.0_dp) then
            d_b = hermite(tab_i, d_i, r)*jprof(r)/r
        else
            d_b = 0.0_dp
        end if
    end function d_b

    real(dp) function ival(r)
        real(dp), intent(in) :: r
        ival = hermite(tab_i, d_i, r)
    end function ival

    real(dp) function qcyl(r)
        real(dp), intent(in) :: r
        qcyl = q0*r*r/(2.0_dp*ival(r))
    end function qcyl

    real(dp) function psi_of_r(r)
        !! psi = R0 int B_theta = (2 B0/q0) int I/r
        real(dp), intent(in) :: r
        psi_of_r = 2.0_dp*b0/q0*hermite(tab_p, d_p, r)
    end function psi_of_r

    real(dp) function bz_of_r(r)
        !! B_z^2/2 = B0^2/2 + c^2 int_r^rmax I j/r, c = 2 B0/(R0 q0)
        real(dp), intent(in) :: r
        real(dp) :: c
        c = 2.0_dp*b0/(r0*q0)
        bz_of_r = sqrt(b0**2 + 2.0_dp*c*c*(tab_b(ntab) - hermite(tab_b, d_b, r)))
    end function bz_of_r

    real(dp) function r_of_psi(psi)
        real(dp), intent(in) :: psi
        real(dp) :: lo, hi, mid
        integer :: it
        lo = 0.0_dp
        hi = rtab_max
        do it = 1, 80
            mid = 0.5_dp*(lo + hi)
            if (psi_of_r(mid) < psi) then
                lo = mid
            else
                hi = mid
            end if
        end do
        r_of_psi = 0.5_dp*(lo + hi)
    end function r_of_psi

    real(dp) function r_of_q(qt)
        !! q increases outward in both profiles
        real(dp), intent(in) :: qt
        real(dp) :: lo, hi, mid
        integer :: it
        lo = 1.0e-3_dp
        hi = 1.0_dp
        do it = 80, 1, -1
            mid = 0.5_dp*(lo + hi)
            if (qcyl(mid) < qt) then
                lo = mid
            else
                hi = mid
            end if
        end do
        r_of_q = 0.5_dp*(lo + hi)
    end function r_of_q

    ! ------------------------------------------------------- ODE Delta'

    subroutine tear_rhs(t, y, dydt, ctx)
        real(dp), intent(in) :: t
        real(dp), intent(in) :: y(:)
        real(dp), intent(out) :: dydt(:)
        class(*), intent(in), optional :: ctx
        real(dp) :: gg
        gg = djprof(t)*t/(ival(t)*(1.0_dp - nmode*qcyl(t)/mres))
        dydt(1) = y(2)
        dydt(2) = -y(2)/t + (mres*mres/(t*t) + gg)*y(1)
    end subroutine tear_rhs

    subroutine integrate(ra, rb, y)
        !! Adaptive Cash-Karp from ra to rb (either direction).
        real(dp), intent(in) :: ra, rb
        real(dp), intent(inout) :: y(2)
        real(dp) :: t, h, y5(2), yerr(2), k(2, 6), yt(2), err, sc
        integer :: nfev, it
        t = ra
        h = sign(1.0e-4_dp, rb - ra)
        nfev = 0
        do it = 1, 10000000
            if ((t + h - rb)*sign(1.0_dp, rb - ra) > 0) h = rb - t
            call cash_karp_step(tear_rhs, t, y, h, .false., k(:, 1), k(:, 2), &
                k(:, 3), k(:, 4), k(:, 5), k(:, 6), yt, y5, yerr, nfev)
            sc = 1.0e-13_dp*(abs(y(1)) + abs(y(2))*abs(h)) + 1.0e-300_dp
            err = maxval(abs(yerr))/sc
            if (err <= 1.0_dp) then
                t = t + h
                y = y5
                if (abs(t - rb) <= 1.0e-15_dp*abs(rb)) exit
            end if
            h = h*min(4.0_dp, max(0.1_dp, 0.9_dp*err**(-0.2_dp)))
        end do
    end subroutine integrate

    real(dp) function ode_delta(rs)
        !! r_s Delta' with Richardson extrapolation in the matching distance
        !! eta (error terms eta ln^2 eta, eta ln eta, eta).
        real(dp), intent(in) :: rs
        integer, parameter :: ne = 6
        real(dp) :: eta(ne), d(ne), yl(2), yr(2), am(ne, 4), rhs(ne), ra
        real(dp) :: work(64), le
        integer :: i, info
        ra = 1.0e-3_dp
        do i = 1, ne
            eta(i) = 1.0e-3_dp*0.25_dp**(i - 1)
            yl = [ra**mres, mres*ra**(mres - 1)]
            call integrate(ra, rs - eta(i), yl)
            yr = [1.0_dp, -real(mres, dp)]
            call integrate(1.0_dp, rs + eta(i), yr)
            d(i) = rs*(yr(2)/yr(1) - yl(2)/yl(1))
            le = log(eta(i))
            am(i, :) = [1.0_dp, eta(i)*le*le, eta(i)*le, eta(i)]
            rhs(i) = d(i)
        end do
        call dgels('N', ne, 4, 1, am, ne, rhs, ne, work, size(work), info)
        ode_delta = rhs(1)
    end function ode_delta

    subroutine step_closed_form()
        real(dp) :: rs, exact, dl(3), ws(3)
        integer :: i
        prof = 2
        q0 = 1.2_dp
        r1 = 0.5_dp
        nmode = 1
        mres = 2
        exact = -2.0_dp*mres*(mres - 1 - nmode*q0) &
            /(mres - 1 - nmode*q0 + (nmode*q0/mres)**mres)
        ws = [0.01_dp, 0.005_dp, 0.0025_dp]
        print '(a,f12.6)', "  step closed form r_s Delta' =", exact
        do i = 1, 3
            wstep = ws(i)
            call make_tables()
            rs = r_of_q(real(mres, dp)/nmode)
            dl(i) = ode_delta(rs)
            print '(a,f7.4,a,f12.6,a,es10.2)', "    ODE, width", wstep, ":", dl(i), &
                "  error", abs(dl(i) - exact)
        end do
        call record(abs(dl(3) - exact) < 2.0e-2_dp*abs(exact) .and. &
            abs(dl(3) - exact) < abs(dl(1) - exact), &
            "ODE Delta' of the smoothed step approaches the closed form")
    end subroutine step_closed_form

    ! ------------------------------------------------------- equilibrium

    subroutine make_equilibrium()
        integer, parameter :: ng = 401, nf = 1025
        real(dp) :: psirz(ng, ng), fpol(nf), r, z, x, pa
        integer :: i, j, st
        do j = 1, ng
            do i = 1, ng
                r = r0 - rtab_max + 2*rtab_max*(i - 1)/(ng - 1)
                z = -rtab_max + 2*rtab_max*(j - 1)/(ng - 1)
                psirz(i, j) = psi_of_r(min(rtab_max, sqrt((r - r0)**2 + z**2)))
            end do
        end do
        pa = psi_of_r(1.0_dp)
        do i = 1, nf
            x = (i - 1.0_dp)/(nf - 1)
            fpol(i) = r0*bz_of_r(r_of_psi(x*pa))
        end do
        call psi_grid_init(g, psirz, r0 - rtab_max, 2*rtab_max, 0.0_dp, 2*rtab_max, &
            r0 + 0.01_dp, 0.01_dp, pa, 1.0e6_dp, fpol, st)
        if (st /= 0) error stop "equilibrium"
    end subroutine make_equilibrium

    subroutine polar_map(u, theta, x)
        real(dp), intent(in) :: u, theta
        real(dp), intent(out) :: x(2)
        real(dp) :: rho
        logical :: ok
        if (u <= 0.0_dp) then
            x = [g%r_axis, g%z_axis]
            return
        end if
        call flux_ray_point(g, u*u, theta, rho, ok)
        if (.not. ok) error stop "ray"
        x = [g%r_axis + rho*cos(theta), g%z_axis + rho*sin(theta)]
    end subroutine polar_map

    subroutine mercier_zero_beta()
        real(dp) :: di, qa, dq
        integer :: st
        prof = 1
        q0 = 1.2_dp
        r0 = 10.0_dp
        call make_tables()
        call make_equilibrium()
        call mercier_di(g, 0.5_dp, 512, di, qa, dq, st)
        print '(a,es12.4,a,f8.4,a,f8.4)', "  D_I + 1/4 =", di + 0.25_dp, &
            ",  |q| =", qa, ", cylinder q =", qcyl(r_of_psi(0.5_dp*psi_of_r(1.0_dp)))
        call record(st == 0 .and. abs(di + 0.25_dp) < 1.0e-12_dp, &
            "Mercier D_I = -1/4 at zero pressure gradient")
        ! a rational closer to the boundary than the difference step
        call mercier_di(g, 0.9996_dp, 512, di, qa, dq, st)
        call record(st == 0 .and. abs(di + 0.25_dp) < 1.0e-12_dp, &
            "Mercier D_I = -1/4 at psin = 0.9996 (step inside the plasma)")
    end subroutine mercier_zero_beta

    subroutine layer_rings_fine_mesh()
        !! Layer rings graded with layer_xmax on a fine mesh (uniform u
        !! spacing below the outermost fit offset) contain every fit ring
        !! that layer_init requires.
        real(dp), parameter :: delta = 5.0e-4_dp, ratio = 2.0_dp
        integer, parameter :: nfit = 5
        real(dp) :: ps(2), xmax
        real(dp), allocatable :: extra(:), lab(:), ul(:)
        type(rational_layer_t) :: lay
        integer :: nrad, l, st
        logical :: ok
        ps = [0.49_dp, 0.76_dp]
        ok = .true.
        do nrad = 96, 512, 32
            extra = sqrt(ps)
            do l = 1, 2
                xmax = layer_xmax(ps(l), ps(3 - l:3 - l), &
                    delta, ratio, nfit, nrad)
                call layer_labels(ps(l), delta, ratio, xmax, (0.5_dp/nrad)**2, &
                    1 - 2.0e-4_dp, lab)
                extra = [extra, sqrt(lab)]
            end do
            call ring_labels_merge(nrad, extra, 0.3_dp, ul)
            call ring_labels_smooth(ul, 2*ratio)
            do l = 1, 2
                call layer_init(ul, ps(l), 1.5_dp, 2.0_dp, -0.3_dp, 3 + l, delta, &
                    ratio, nfit, lay, st)
                ok = ok .and. st == 0
            end do
        end do
        call record(ok, "layer fit rings present for nrad = 96 .. 512")
    end subroutine layer_rings_fine_mesh

    ! ----------------------------------------------------------- FE ladder

    ! ------------------------------------------------ R-B2 strip coupling

    subroutine strip_suite()
        !! R-B2 first increment (docs/response-unified.md (B)): the resolved
        !! resistive strip of kin6d_response_strip (single fluid, no
        !! rotation, resonant harmonic resolved in the strip, cylinder-
        !! equivalent local equations) coupled to the FEEC torus outer region
        !! by trace continuity on the two strip edge rings: the strip's
        !! Dirichlet-to-Neumann map D(gamma), the current-free interior that
        !! the FEEC sees between the edges (map V) and the jump of
        !! d Phi_m/dr per unit edge sheet give the edge sheets lam = S Phi,
        !! imposed as the closure rows S Phi - lam = 0 (layer_closure_t w, wl).
        !! Tanh step (q0 1.2, r1 0.5, width 0.02), 2/1 strip, 3/1 and 4/1 L4
        !! sheets.  Pre-declared rules (RB2 handoff, before the first run):
        !!  (a) eta = 0, gamma = 1e-6: strip net current / L4 sheet current on
        !!      the rational ring, same mesh and drive: |ratio - 1| <= 1e-3;
        !!  (b) as first declared: R0/a = 40: torus growth rate / R-B1 sharp-
        !!      step cylinder rate (same eta, r1 = sqrt(2 I_tot)), S = 1e6, 1e7,
        !!      1e8: |ratio - 1| <= 2 %;
        !!  consistency: the coupled closure solve equals the Schur form
        !!      (I - S M)^-1 S P0 to the GMRES tolerance (5e-4 relative).
        !! (b) re-declared 2026-10-14 (controller): the sharp-step cylinder
        !! (r_s Delta' = 5.000) is a different problem from the tanh profile
        !! solved here (width 0.02, ODE r_s Delta' = 5.247, about +3.9 % in
        !! gamma ~ Delta'^(4/5)), so it does not test the torus coupling.  (b)
        !! is now: torus growth rate / cylinder growth rate of the SAME tanh
        !! profile (strip coupling with the cylinder tearing-ODE outer map,
        !! step_outer_map(.true.)), same S values, |ratio - 1| <= 2 % (limit
        !! unchanged), with the gap shrinking under nrad refinement (ladder in
        !! docs/response-unified.md, R-B2 increment 1).  The gated run uses
        !! the smallest nrad of that ladder at which (b) passes.  The old
        !! ratio to the sharp-step rate stays printed ("rel").
        !! (a) retains the 2026-10-14 interface-sheet diagnostics: their sum
        !! agrees with L4 to 1e-3 at xd=0.01 and its error has order >=1.8
        !! between xd=0.02 and 0.01. That sum is the magnetic replacement's
        !! current, not the resolved strip current. Exact cylinder moments
        !! explain its finite-width error (response-strip-current.md).
        !! The physical-current row now supplies the original gate as well:
        !! actual strip current / L4 within 1e-3 at BOTH xd=0.02 and 0.01.
        !! No fitted correction is applied to the interface amplitudes.
        !! KIN6D_STRIP_LADDER="r0 xd nrad" runs one ladder case (no gates).
        character(len=128) :: env
        real(dp) :: r0l, xdl, ierr, cerr, ierr2, cerr2, order, physical, physical2
        integer :: nr, i
        integer, parameter :: nrad_gate = 64
        real(dp), parameter :: etas(3) = [1.0e-6_dp, 1.0e-7_dp, 1.0e-8_dp]
        real(dp), parameter :: xd_gate = 1.0e-2_dp
        real(dp) :: none(0)
        complex(dp) :: gt(3), gc(3), gs(3), g0(0)
        call get_environment_variable('KIN6D_STRIP_LADDER', env)
        if (len_trim(env) > 0) then
            read (env, *) r0l, xdl, nr
            call fe_strip(r0l, xdl, nr, etas, gt, gc, gs, ierr, cerr)
            stop
        end if
        ! (a) at twice the gated width, no growth rates: the xd-order reference
        call fe_strip(40.0_dp, 2*xd_gate, nrad_gate, none, g0, g0, g0, ierr2, &
            cerr2, physical2)
        call fe_strip(40.0_dp, xd_gate, nrad_gate, etas, gt, gc, gs, ierr, cerr, &
            physical)
        order = log(ierr2/ierr)/log(2.0_dp)
        print '(a,es10.3,a,es10.3,a,f6.3)', "  (a) xd 0.02 rel", ierr2, &
            ", xd 0.01 rel", ierr, ", xd order", order
        call record(ierr <= 1.0e-3_dp, "(a) ideal limit: strip current = L4 sheet")
        call record(order >= 1.8_dp, "(a) strip -> sheet error converges as xd^2")
        call record(physical2 <= 1.0e-3_dp .and. physical <= 1.0e-3_dp, &
            '(a) physical strip current = L4 sheet at original and half width')
        do i = 1, 3
            call record(abs(gt(i)/gs(i) - 1) <= 2.0e-2_dp, &
                "(b) torus growth rate = same-profile cylinder rate")
        end do
        call record(cerr <= 5.0e-4_dp, "coupled closure = Schur complement form")
        call record(merr <= 1.0e-3_dp, &
            "coupling with the sharp-step outer map = R-B1 cylinder rate")
        if (.not. all_passed) error stop 1
        print '(a)', "test_mhd_closure strip: all checks passed"
        stop
    end subroutine strip_suite

    subroutine fe_strip(r0a, xd, nrad, etas, gtor, gcyl, gsame, ierr, cerr, current_error)
        !! gcyl: R-B1 sharp-step cylinder rate; gsame: cylinder rate of the
        !! same (tanh) profile, the reference of gate (b).
        real(dp), intent(in) :: r0a, xd, etas(:)
        integer, intent(in) :: nrad
        complex(dp), intent(out) :: gtor(:), gcyl(:), gsame(:)
        real(dp), intent(out) :: ierr, cerr
        real(dp), intent(out), optional :: current_error
        complex(dp) :: physical_current
        nmode = 1
        mres = 2
        q0 = 1.2_dp
        prof = 2
        wstep = 0.02_dp
        call make_tables()
        r0 = r0a
        call make_equilibrium()
        pa = psi_of_r(1.0_dp)
        call surface_build(g, 0.25_dp, 256, s, st)
        qsgn = sign(1.0_dp, s%q)
        ! fe_strip is called twice by strip_suite: release the previous call
        if (allocated(psin_s)) deallocate (psin_s)
        if (allocated(mr)) deallocate (mr)
        if (allocated(ras)) deallocate (ras)
        if (allocated(f)) deallocate (f)
        if (allocated(a)) deallocate (a)
        if (allocated(sr)) deallocate (sr)
        if (allocated(lam)) deallocate (lam)
        if (allocated(phit)) deallocate (phit)
        if (allocated(phe)) deallocate (phe)
        if (allocated(b2cl)) deallocate (b2cl)
        allocate (psin_s(0), mr(0))
        do k = 1, 20
            qt = real(k, dp)/nmode
            found = find_surface_q(g, qsgn*qt, 0.02_dp, 0.98_dp, 256, 1.0e-13_dp, &
                x, st)
            if (.not. found) cycle
            psin_s = [psin_s, x]
            mr = [mr, nint(qsgn*qt*nmode)]
            if (k == mres) itest = size(mr)
        end do
        nrat = size(mr)
        msgn = sign(1.0_dp, real(mr(itest), dp))
        ! rings: rationals, strip edges and two graded rings on each side
        extra = sqrt(psin_s)
        do k = 0, 2
            extra = [extra, sqrt(psin_s(itest) - xd*(1 + 0.25_dp*k*(k + 1))), &
                sqrt(psin_s(itest) + xd*(1 + 0.25_dp*k*(k + 1)))]
        end do
        do k = 1, 2
            extra = [extra, sqrt(psin_s(itest) - xd*(1 - 0.25_dp*k*(k + 1)/2)), &
                sqrt(psin_s(itest) + xd*(1 - 0.25_dp*k*(k + 1)/2))]
        end do
        block
            character(len=16) :: ev2
            integer :: nsr
            nsr = 0
            call get_environment_variable('KIN6D_STRIP_STEPRINGS', ev2)
            if (len_trim(ev2) > 0) read (ev2, *) nsr
            ! rings across the tanh current step, spacing wstep/2
            if (nsr > 0) then
                do k = -nsr, nsr
                    extra = [extra, sqrt(psi_of_r(r1 + 0.5_dp*k*wstep)/pa)]
                end do
            end if
        end block
        call ring_labels_merge(nrad, extra, 0.3_dp, ul)
        nr = size(ul)
        do side = 1, 2
            je(side) = minloc(abs(ul**2 - (psin_s(itest) + (2*side - 3)*xd)), 1)
            if (abs(ul(je(side))**2 - psin_s(itest) - (2*side - 3)*xd) > 1.0e-8_dp) &
                error stop "strip edge ring"
        end do
        if (je(1) + 2 >= je(2) - 2) error stop "strip interior rings"
        call flux_mesh_build(ul, 8, 256, 2, polar_map, fm, st)
        if (st /= 0) error stop "mesh"
        call meridian_feec_build(fm%mm, 2, fs, st)
        call ampere_n_build(fs, nmode, amp, st)
        if (st /= 0) error stop "ampere"
        call meridian_point_basis_init(fs, pb, st)
        allocate (ras(nr))
        do j = 1, nr
            call surface_build(g, min(ul(j), 0.999_dp)**2, 1024, s, st)
            call ring_angles_set(s%q, s%theta_s(:, 2), s%nu, ras(j))
        end do
        ! jump of d Phi/dr per unit sheet on each edge ring (vacuum sheet field,
        ! one-sided quadratic slopes in psin on two rings each side)
        allocate (f(amp%nn), a(amp%nn))
        do side = 1, 2
            j = je(side)
            f = (0.0_dp, 0.0_dp)
            call ring_sheet_load(fs, pb, fm, j, ras(j), nmode, mr(itest), &
                (1.0_dp, 0.0_dp), f)
            call ampere_n_solve(amp, f, a, st)
            call ampere_n_curl(amp, a, b2)
            jr(side) = slope3(j, j + 1, j + 2) - slope3(j, j - 1, j - 2)
            re(side) = r_of_psi(ul(j)**2*pa)
            jr(side) = jr(side)*(psi_of_r(re(side) + 1.0e-5_dp) &
                - psi_of_r(re(side) - 1.0e-5_dp))/(2.0e-5_dp*pa)
        end do
        ! closures: 3/1, 4/1 (and 1/1 if present) L4 sheets; 2/1 strip rows
        ns = nrat + 1
        allocate (sr(nr), lam(ns), phit(ns), phe(ns), b2cl(fs%sp(2)%ndof, ns))
        call build_closure(.true.)
        do j = 1, ns
            f = (0.0_dp, 0.0_dp)
            call ring_sheet_load(fs, pb, fm, cl%sheet_ring(j), cl%sheet_ra(j), nmode, &
                cl%sheet_m(j), (1.0_dp, 0.0_dp), f)
            call ampere_n_solve(amp, f, a, st)
            call ampere_n_curl(amp, a, b2)
            b2cl(:, j) = b2
        end do
        call build_rings(.true.)
        ! One GMRES tolerance on every backend: with the SuperLU K solve the old
        ! 3e-5 fallback left a solver floor (direct/Schur 3.7e-4) above the
        ! xd 0.01 discretisation error 1.55e-4; 1e-7 converges within 400
        ! iterations on both SuperLU and MUMPS and gives the same (a) values.
        tol = 1.0e-7_dp
        block
            character(len=32) :: ev
            ! KIN6D_STRIP_TOL="tol maxit" overrides the GMRES settings (studies)
            call get_environment_variable('KIN6D_STRIP_TOL', ev)
            if (len_trim(ev) > 0) read (ev, *) tol, gmres_max
        end block
        ! outer response to prescribed edge sheets (no drive): M
        cl%wl(ns - 1:ns, ns - 1:ns) = reshape([1, 0, 0, 1], [2, 2])
        do k = 1, 2
            cl%target = (0.0_dp, 0.0_dp)
            cl%target(ns - 2 + k) = (1.0_dp, 0.0_dp)
            call solve(.false.)
            mo(:, k) = phit(ns - 1:ns)
        end do
        ! driven, lam = 0: P0
        cl%target = (0.0_dp, 0.0_dp)
        call solve(.true.)
        p0 = phit(ns - 1:ns)
        ! linearity of the prescribed-sheet path: P0 + M lx at arbitrary lx
        block
            complex(dp) :: lx(2)
            lx = [(0.3_dp, 0.7_dp), (-1.1_dp, 0.2_dp)]
            cl%target(ns - 1:ns) = lx
            call solve(.true.)
            call record(maxval(abs(phit(ns - 1:ns) - p0 - matmul(mo, lx))) <= &
                1.0e-3_dp*maxval(abs(matmul(mo, lx))), &
                "prescribed edge sheets: superposition P0 + M lam")
            cl%target = (0.0_dp, 0.0_dp)
        end block
        ! cylinder-equivalent strip
        eq%m = mres
        eq%n = nmode
        eq%q0 = q0
        eq%r1 = sqrt(2*ival(1.0_dp))
        rs = step_rational(eq)
        fp = step_fprime(eq)
        dpr = step_delta_prime(eq)
        rb(1, :) = [re(1)**mres, re(1)**(-mres)]
        rb(2, :) = [re(2)**mres, re(2)**(-mres)]
        vmap(1, :) = mres*[re(1)**(mres - 1), -re(1)**(-mres - 1)]
        vmap(2, :) = mres*[re(2)**(mres - 1), -re(2)**(-mres - 1)]
        vmap = matmul(vmap, inv2(cmplx(rb, 0.0_dp, dp)))
        print '(a,f6.1,a,es9.2,a,i0,a,2f9.5,a,f9.5,a,f9.4,a,2es11.3)', 'strip R0/a', &
            r0, ' xd', xd, ' rings ', nr, ' r_e', re, ' r_s', rs, ' rs Delta''', &
            dpr, ' |junit|', abs(jr)
        print '(a,f9.5)', '  same-profile (tanh) cylinder ODE r_s Delta''', &
            ode_delta(rs)
        ! (a) ideal limit
        wl = layer_width(0.0_dp, (1.0e-6_dp, 0.0_dp), fp)
        call smat(0.0_dp, (1.0e-6_dp, 0.0_dp), wl)
        lam2 = matmul(inv2(eye2() - matmul(smap, mo)), matmul(smap, p0))
        ! The Schur sheets replace the interior magnetic field for the outer
        ! solve; their sum is not the current in the resolved strip. Apply the
        ! strip's physical-current row to its actual two-ring flux traces.
        ! In the peak-flux convention Phi=2 R0 r B_r, a cylindrical unit
        ! stream-function sheet has int r L(Phi) dr = 2 mu0 m^2 R0 Lambda.
        ! gen/check_response_strip_current.f90 checks this normalization exactly.
        physical_current = sum(current_map*(p0 + matmul(mo, lam2))) &
            /(2*mu0*mres**2*r0)
        ! direct coupled closure solve
        cl%w(ns - 1:ns, ns - 1:ns) = smap
        cl%wl(ns - 1:ns, ns - 1:ns) = -eye2()
        cl%target = (0.0_dp, 0.0_dp)
        call solve(.true.)
        cerr = maxval(abs(lam(ns - 1:ns) - lam2))/maxval(abs(lam2))
        print '(a,4es12.4)', '  direct: |S Phi - lam|, |Phi - P0 - M lam|', &
            abs(matmul(smap, phit(ns - 1:ns)) - lam(ns - 1:ns)), &
            abs(phit(ns - 1:ns) - p0 - matmul(mo, lam(ns - 1:ns)))
        ! L4 reference: one sheet on the rational ring, all half rings active
        call build_closure(.false.)
        do j = ns - 1, ns
            f = (0.0_dp, 0.0_dp)
            call ring_sheet_load(fs, pb, fm, cl%sheet_ring(j), cl%sheet_ra(j), nmode, &
                cl%sheet_m(j), (1.0_dp, 0.0_dp), f)
            call ampere_n_solve(amp, f, a, st)
            call ampere_n_curl(amp, a, b2)
            b2cl(:, j) = b2
        end do
        call build_rings(.false.)
        call solve(.true.)
        lam_l4 = lam(ns - 1) + lam(ns)
        ierr = abs((lam2(1) + lam2(2))/lam_l4 - 1)
        if (present(current_error)) current_error = abs(physical_current/lam_l4 - 1)
        print '(a,2es13.5,a,es10.2)', '  (a) physical strip current', &
            physical_current, ' rel to L4', abs(physical_current/lam_l4 - 1)
        print '(a,2es13.5,a,2es13.5,a,es10.2,a,es10.2)', '  (a) strip current', &
            lam2(1) + lam2(2), ' L4', lam_l4, ' rel', ierr, ' direct/Schur', cerr
        ! (b) growth rates: det(I - S(gamma) M) = 0, mesh from the cylinder rate
        merr = 0.0_dp
        do k = 1, size(etas)
            call ri_dispersion_gamma(etas(k), fp, dpr, gr, st)
            sopt%xs = max(rs - re(1), re(2) - rs)
            call strip_growth_rate(eq, etas(k), cmplx(1.3_dp*gr, 0.0_dp, dp), sopt, &
                gcyl(k), ssol, it, st)
            wl = layer_width(etas(k), gcyl(k), fp)
            call rate_secant(etas(k), gcyl(k), wl, gtor(k), it)
            ! same coupling with the exact sharp-step outer map: must give gcyl
            mfe = mo
            mo = step_outer_map(.false.)
            call rate_secant(etas(k), 1.05_dp*gcyl(k), wl, gn, it)
            ! diagnostic: cylinder map of the smoothed (tanh) profile itself
            mo = step_outer_map(.true.)
            call rate_secant(etas(k), 1.05_dp*gcyl(k), wl, gb, it)
            gsame(k) = gb
            mo = mfe
            merr = max(merr, abs(gn/gcyl(k) - 1))
            print '(a,es9.2,a,2es13.5,a,2es13.5,a,es10.2,a,es10.2)', '  (b) S', &
                1/etas(k), ' gamma torus', gtor(k), ' cylinder', gcyl(k), ' rel', &
                abs(gtor(k)/gcyl(k) - 1), ' step map', abs(gn/gcyl(k) - 1)
            print '(a,es10.2,a,es10.2)', '      same-profile cylinder / step', &
                abs(gb/gcyl(k) - 1), ', torus / same-profile cylinder (b)', &
                abs(gtor(k)/gb - 1)
        end do

    end subroutine fe_strip

    subroutine rate_secant(eta, g0, w, gam, it2)
        !! Root of det(I - S(gamma) mo) by the secant method from g0.
        real(dp), intent(in) :: eta, w
        complex(dp), intent(in) :: g0
        complex(dp), intent(out) :: gam
        integer, intent(out) :: it2
        complex(dp) :: g1, g2, f1, f2, g3
        g1 = g0
        g2 = 1.05_dp*g0
        f1 = disp(eta, g1, w)
        f2 = disp(eta, g2, w)
        do it2 = 1, 40
            g3 = g2 - f2*(g2 - g1)/(f2 - f1)
            g1 = g2
            f1 = f2
            g2 = g3
            f2 = disp(eta, g2, w)
            if (abs(g2 - g1) <= 1.0e-10_dp*abs(g2)) exit
        end do
        gam = g2
    end subroutine rate_secant

    function step_outer_map(smooth) result(mc)
        !! Edge fluxes per unit edge sheet of the cylinder: regular solution
        !! f inside the sheet (sharp step of eq, a r^m + b r^-m; smooth: the
        !! tearing ODE of the tabulated tanh profile up to the inner edge and
        !! its current-free continuation across the strip), r^-m outside, jump jr;
        !! the cylinder counterpart of the FEEC map mo.
        logical, intent(in) :: smooth
        complex(dp) :: mc(2, 2)
        real(dp) :: gs, as, bs, fin(2), dfin(2), p, y(2)
        integer :: i2, k2
        gs = 2.0_dp/(1.0_dp - eq%n*eq%q0/eq%m)
        as = 1.0_dp - gs/(2*eq%m)
        bs = gs*eq%r1**(2*eq%m)/(2*eq%m)
        do i2 = 1, 2
            fin(i2) = as*re(i2)**mres + bs*re(i2)**(-mres)
            dfin(i2) = mres*(as*re(i2)**(mres - 1) - bs*re(i2)**(-mres - 1))
            if (smooth .and. i2 == 1) then
                y = [1.0e-3_dp**mres, mres*1.0e-3_dp**(mres - 1)]
                call integrate(1.0e-3_dp, re(1), y)
                fin(1) = y(1)
                dfin(1) = y(2)
            else if (smooth) then
                ! current-free continuation across the strip (as the FEEC
                ! interior map V); the ODE must not be integrated through r_s
                as = (fin(1) + re(1)*dfin(1)/mres)/(2*re(1)**mres)
                bs = (fin(1) - re(1)*dfin(1)/mres)*re(1)**mres/2
                fin(2) = as*re(2)**mres + bs*re(2)**(-mres)
                dfin(2) = mres*(as*re(2)**(mres - 1) - bs*re(2)**(-mres - 1))
            end if
        end do
        do k2 = 1, 2
            p = real(jr(k2), dp)/(-mres/re(k2) - dfin(k2)/fin(k2))
            do i2 = 1, 2
                if (i2 <= k2) then
                    mc(i2, k2) = p*fin(i2)/fin(k2)
                else
                    mc(i2, k2) = p*(re(k2)/re(i2))**mres
                end if
            end do
        end do
    end function step_outer_map

    complex(dp) function slope3(j0, j1, j2) result(d)
        integer, intent(in) :: j0, j1, j2
        real(dp) :: x0, x1, x2
        complex(dp) :: f0, f1, f2
        x0 = ul(j0)**2
        x1 = ul(j1)**2
        x2 = ul(j2)**2
        f0 = ring_flux_harmonic(fs, pb, fm, j0, ras(j0), nmode, b2, mr(itest))
        f1 = ring_flux_harmonic(fs, pb, fm, j1, ras(j1), nmode, b2, mr(itest))
        f2 = ring_flux_harmonic(fs, pb, fm, j2, ras(j2), nmode, b2, mr(itest))
        d = f0*(2*x0 - x1 - x2)/((x0 - x1)*(x0 - x2)) &
            + f1*(x0 - x2)/((x1 - x0)*(x1 - x2)) &
            + f2*(x0 - x1)/((x2 - x0)*(x2 - x1))
    end function slope3

    subroutine build_closure(strip)
        logical, intent(in) :: strip
        integer :: l2, q2
        if (allocated(cl%sheet_ring)) deallocate (cl%sheet_ring, cl%sheet_m, &
            cl%sheet_ra, cl%flux_ring, cl%flux_m, cl%flux_ra, cl%w, cl%target, &
            cl%wl)
        cl%ns = ns
        cl%nf = ns
        allocate (cl%sheet_ring(ns), cl%sheet_m(ns), cl%sheet_ra(ns), &
            cl%flux_ring(ns), cl%flux_m(ns), cl%flux_ra(ns), cl%w(ns, ns), &
            cl%target(ns), cl%wl(ns, ns))
        cl%w = (0.0_dp, 0.0_dp)
        cl%wl = (0.0_dp, 0.0_dp)
        cl%target = (0.0_dp, 0.0_dp)
        q2 = 0
        do l2 = 1, nrat
            if (l2 == itest) cycle
            q2 = q2 + 1
            cl%sheet_ring(q2) = minloc(abs(ul**2 - psin_s(l2)), 1)
            cl%sheet_m(q2) = mr(l2)
        end do
        if (strip) then
            cl%sheet_ring(ns - 1:ns) = je
            cl%sheet_m(ns - 1:ns) = mr(itest)
        else
            ! L4: the rational ring carries the sheet; ns - 1 is a dummy
            ! sheet on the inner edge fixed to zero amplitude
            cl%sheet_ring(ns - 1) = je(1)
            cl%sheet_ring(ns) = minloc(abs(ul**2 - psin_s(itest)), 1)
            cl%sheet_m(ns - 1:ns) = mr(itest)
        end if
        cl%flux_ring = cl%sheet_ring
        cl%flux_m = cl%sheet_m
        do l2 = 1, ns
            cl%sheet_ra(l2) = ras(cl%sheet_ring(l2))
            cl%flux_ra(l2) = ras(cl%flux_ring(l2))
            if (l2 < ns - 1) cl%w(l2, l2) = (1.0_dp, 0.0_dp)
        end do
        if (.not. strip) then
            cl%wl(ns - 1, ns - 1) = (1.0_dp, 0.0_dp)
            cl%w(ns, ns) = (1.0_dp, 0.0_dp)
        end if
    end subroutine build_closure

    subroutine build_rings(strip)
        logical, intent(in) :: strip
        integer :: j2, k2, st2
        do j2 = 2, nr
            call sample_ring_build(g, opt, fm, fs, pb, 0.5_dp*(ul(j2 - 1) + ul(j2)), &
                j2, np, sr(j2 - 1), st2)
            if (st2 /= 0) error stop "half ring"
            if (strip .and. sr(j2 - 1)%u > ul(je(1)) .and. sr(j2 - 1)%u < ul(je(2))) &
                sr(j2 - 1)%active = .false.
        end do
        call sample_ring_build(g, opt, fm, fs, pb, 1.0_dp, nr, np, sr(nr), st2)
        do j2 = 1, nr
            do k2 = 1, np
                call drive(sr(j2)%x(:, k2), sr(j2)%bc(:, k2))
            end do
        end do
        ! resonant drive fluxes on the closure flux rings
        do k2 = 1, ns
            call ring_points(fm, cl%flux_ring(k2), xq, nds, th, tt, tau, wq)
            if (allocated(bc)) deallocate (bc)
            allocate (bc(size(th)))
            do j2 = 1, size(th)
                call drive(xq(:, j2), bb)
                bc(j2) = (bb(1)*nds(1, j2) + bb(3)*nds(2, j2))*xq(1, j2)
            end do
            phe(k2) = ring_flux_from_values(cl%flux_ra(k2), nmode, cl%flux_m(k2), &
                th, bc, wq)
        end do
    end subroutine build_rings

    subroutine solve(driven)
        logical, intent(in) :: driven
        complex(dp) :: ph(ns)
        real(dp) :: tl
        integer :: j2
        type(sample_ring_t), allocatable :: s0(:)
        s0 = sr
        if (.not. driven) then
            do j2 = 1, nr
                s0(j2)%bc = (0.0_dp, 0.0_dp)
            end do
            ph = (0.0_dp, 0.0_dp)
        else
            ph = phe
        end if
        a = (0.0_dp, 0.0_dp)
        ! the edge-sheet responses M enter (I - S M)^-1 amplified: solve them
        ! three orders tighter than the driven problems
        tl = tol
        if (.not. driven) tl = 1.0e-3_dp*tol
        call ideal_response(fs, pb, fm, g, amp, opt, nmode, s0, cl, b2cl, ph, a, &
            lam, phit, tl, gmres_max, info, its, x)
        call record(info == 0, "GMRES converged")
    end subroutine solve

    subroutine smat(eta, gam, w)
        real(dp), intent(in) :: eta, w
        complex(dp), intent(in) :: gam
        integer :: n2, h2, st2
        sopt%xs = max(rs - re(1), re(2) - rs)
        call strip_mesh(rs, w, sopt, rr)
        n2 = size(rr)
        h2 = (n2 + 1)/2
        rr(1:h2) = rs + (rr(1:h2) - rs)*(rs - re(1))/sopt%xs
        rr(h2:n2) = rs + (rr(h2:n2) - rs)*(re(2) - rs)/sopt%xs
        call strip_dirichlet_map(eq, eta, gam, rr, dmap, st2, current_map)
        if (st2 /= 0) error stop "strip LU"
        smap(1, :) = (vmap(1, :) - dmap(1, :))/jr(1)
        smap(2, :) = (dmap(2, :) - vmap(2, :))/jr(2)
    end subroutine smat

    complex(dp) function disp(eta, gam, w) result(d)
        real(dp), intent(in) :: eta, w
        complex(dp), intent(in) :: gam
        call smat(eta, gam, w)
        amat = eye2() - matmul(smap, mo)
        d = amat(1, 1)*amat(2, 2) - amat(1, 2)*amat(2, 1)
    end function disp

    pure function eye2() result(e)
        complex(dp) :: e(2, 2)
        e = reshape([(1.0_dp, 0.0_dp), (0.0_dp, 0.0_dp), (0.0_dp, 0.0_dp), &
            (1.0_dp, 0.0_dp)], [2, 2])
    end function eye2

    pure function inv2(m2) result(v)
        complex(dp), intent(in) :: m2(2, 2)
        complex(dp) :: v(2, 2), det
        det = m2(1, 1)*m2(2, 2) - m2(1, 2)*m2(2, 1)
        v = reshape([m2(2, 2), -m2(2, 1), -m2(1, 2), m2(1, 1)], [2, 2])/det
    end function inv2

    subroutine study()
        !! KIN6D_CLOSURE_STUDY="prof r0 nrad delta ratio": one FE Delta'.
        character(len=128) :: env
        real(dp) :: delta, ratio, rs, dfe, dode
        call get_environment_variable('KIN6D_CLOSURE_STUDY', env)
        if (len_trim(env) == 0) return
        if (trim(env) == 'edge') then
            call edge_current_jump()
            stop
        end if
        if (trim(env) == 'shieldcyl') then
            call shield_cylinder()
            if (.not. all_passed) error stop 1
            stop
        end if
        if (trim(env) == 'reconnect') then
            call reconnect_cylinder()
            if (.not. all_passed) error stop 1
            stop
        end if
        if (env(1:6) == 'shield') then
            call shield_study(env(7:))
            stop
        end if
        read (env, *) prof, r0, nrad_fe, delta, ratio, gmres_max, nfit_fe
        nmode = 1
        mres = 2
        q0 = 1.2_dp
        wstep = 0.08_dp
        call make_tables()
        rs = r_of_q(real(mres, dp)/nmode)
        dode = ode_delta(rs)
        call make_equilibrium()
        call fe_delta(delta, ratio, rs, dfe)
        print '(a,i2,f6.1,i5,es10.2,f6.2,a,2f12.6,es10.2)', 'study', prof, r0, &
            nrad_fe, &
            delta, ratio, ' ODE FE rel', dode, dfe, (dfe - dode)/dode
        stop
    end subroutine study

    subroutine edge_current_jump()
        !! Uniform current j = 1 up to the plasma edge r = 1 (q = q0, vacuum
        !! outside), driven by the vacuum harmonic r^m e^{-i m theta}: inside
        !! psi = alpha r^m, outside psi = d (r^m + beta r^-m), and the edge
        !! current jump (moved current xi_r [[J_z]]) gives [psi'] = -2 psi(1)/
        !! (1 - n q0/m), so beta = 1/(m - n q0 - 1) and the resonant flux of
        !! the total field is (1 + beta) times that of the drive on every
        !! surface.  In the ideal model this is the field-aligned boundary
        !! sheet (charge conservation of the volume current sigma B leaving
        !! the plasma; at zero beta it equals the moved current xi_n J0).
        real(dp), parameter :: r0s(3) = [10.0_dp, 20.0_dp, 40.0_dp]
        real(dp) :: ratio(3), exact, err(3), rate, rno
        integer :: i
        prof = 3
        q0 = 1.2_dp
        nmode = 1
        mres = 3
        exact = 1.0_dp + 1.0_dp/(mres - nmode*q0 - 1.0_dp)
        call make_tables()
        print '(a,f10.6)', &
            "  edge current jump, m = 3, n = 1, q0 = 1.2: exact (1 + beta) =", &
            exact
        do i = 1, 3
            r0 = r0s(i)
            call make_equilibrium()
            call fe_forced(.true., ratio(i))
            err(i) = abs(ratio(i) - exact)/exact
            print '(a,f5.1,a,f10.6,a,es10.2)', "    R0/a =", r0, &
                ": FEEC ideal response", &
                ratio(i), "  rel. error", err(i)
        end do
        call fe_forced(.false., rno)
        print '(a,f10.6)', "    R0/a = 40 without the boundary sheet:", rno
        rate = log(err(1)/err(3))/log(4.0_dp)
        print '(a,f6.2)', "    error order in a/R0", rate
        call record(err(3) < 1.0e-2_dp .and. err(3) < err(1), &
            "edge current jump: response approaches 1 + 1/(m - n q0 - 1)")
        call record(abs(rno - exact) > 0.5_dp, &
            "edge current jump: the response without the boundary sheet is wrong")
    end subroutine edge_current_jump

    subroutine drive(x, b)
        !! grad_n of u = (x - i s z)^m, s = sign of the harmonic (e^{-i m_s theta})
        real(dp), intent(in) :: x(2)
        complex(dp), intent(out) :: b(3)
        complex(dp) :: w, u
        w = cmplx(x(1) - r0, -msgn*x(2), dp)
        u = w**mres
        b(1) = mres*w**(mres - 1)
        b(3) = cmplx(0.0_dp, -msgn, dp)*mres*w**(mres - 1)
        b(2) = cmplx(0.0_dp, nmode, dp)*u/x(1)
        b = 1.0e-3_dp*b
    end subroutine drive

    subroutine fe_forced(sheet, ratio)
        !! Resonant-harmonic flux of the total field over that of the drive on
        !! the boundary ring (no rational surface in the plasma).
        use kin6d_mhd_response_n, only: ring_angles_t, ring_angles_set, &
            ring_flux_harmonic, ring_flux_from_values, ring_points
        logical, intent(in) :: sheet
        real(dp), intent(out) :: ratio
        integer, parameter :: nrad = 32, np = 256
        type(flux_mesh_t) :: fm
        type(meridian_feec_t) :: fs
        type(meridian_point_basis_t) :: pb
        type(ampere_n_t) :: amp
        type(ideal_options_t) :: opt
        type(sample_ring_t), allocatable :: sr(:)
        type(layer_closure_t) :: cl
        type(ring_angles_t) :: rac
        type(surface_t) :: s
        complex(dp), allocatable :: a(:), b2(:), bc(:)
        complex(dp) :: lam(0), b2s(0, 0), phe(0), pht(0), phv, pht1, bb(3)
        real(dp), allocatable :: ul(:), x(:, :), nds(:, :), th(:), tau(:), wq(:)
        integer, allocatable :: tt(:)
        real(dp) :: res
        integer :: j, k, st, info, its, mk
        allocate (ul(nrad))
        do j = 1, nrad
            ul(j) = real(j, dp)/nrad
        end do
        call flux_mesh_build(ul, 8, 256, 2, polar_map, fm, st)
        call meridian_feec_build(fm%mm, 2, fs, st)
        call ampere_n_build(fs, nmode, amp, st)
        call meridian_point_basis_init(fs, pb, st)
        call surface_build(g, 1.0_dp, 1024, s, st)
        mk = nint(sign(real(mres, dp), s%q))
        msgn = sign(1.0_dp, real(mk, dp))
        opt%boundary_sheet = sheet
        allocate (sr(nrad))
        do j = 2, nrad
            call sample_ring_build(g, opt, fm, fs, pb, 0.5_dp*(ul(j - 1) + ul(j)), j, &
                np, sr(j - 1), st)
        end do
        call sample_ring_build(g, opt, fm, fs, pb, 1.0_dp, nrad, np, sr(nrad), st)
        do j = 1, nrad
            do k = 1, np
                call drive(sr(j)%x(:, k), sr(j)%bc(:, k))
            end do
        end do
        allocate (a(amp%nn))
        a = (0.0_dp, 0.0_dp)
        call ideal_response(fs, pb, fm, g, amp, opt, nmode, sr, cl, b2s, phe, a, lam, &
            pht, 1.0e-7_dp, 300, info, its, res)
        call record(info == 0, "GMRES converged")
        call ampere_n_curl(amp, a, b2)
        call ring_angles_set(s%q, s%theta_s(:, 2), s%nu, rac)
        pht1 = ring_flux_harmonic(fs, pb, fm, nrad, rac, nmode, b2, mk)
        call ring_points(fm, nrad, x, nds, th, tt, tau, wq)
        allocate (bc(size(th)))
        do k = 1, size(th)
            call drive(x(:, k), bb)
            bc(k) = (bb(1)*nds(1, k) + bb(3)*nds(2, k))*x(1, k)
        end do
        phv = ring_flux_from_values(rac, nmode, mk, th, bc, wq)
        ratio = real((pht1 + phv)/phv, dp)
    end subroutine fe_forced

    subroutine fe_single(name, r0a, delta, tol)
        !! FEEC closure Delta' at one layer half-width (profiles with a small
        !! resonant current gradient, where the layer remainder is small).
        character(len=*), intent(in) :: name
        real(dp), intent(in) :: r0a, delta, tol
        real(dp) :: dfe, dode, rs, err
        nmode = 1
        mres = 2
        q0 = 1.2_dp
        call make_tables()
        rs = r_of_q(real(mres, dp)/nmode)
        dode = ode_delta(rs)
        r0 = r0a
        call make_equilibrium()
        call fe_delta(delta, 1.5_dp, rs, dfe)
        err = abs(dfe - dode)/abs(dode)
        print '(a,a,a,f12.6,a,f5.1,a,es9.2,a,f12.6,a,es10.2)', "  ", name, &
            ": ODE r_s Delta' =", dode, ", R0/a =", r0, ", delta", delta, ": FEEC", &
            dfe, "  rel. error", err
        call record(err < tol, "closure Delta' against the cylinder ODE ("//name//")")
    end subroutine fe_single

    subroutine layer_conditioning()
        !! The a1 fit of three rings at ratio 1.25 amplifies flux errors by
        !! |w|_1 = 70 (cond 200) and makes the sheet dipole, a1_in - a1_out of
        !! O(delta), noise dominated: layer_init must reject it, and accept
        !! ratios 1.5 and 2.
        real(dp), parameter :: ratios(3) = [1.25_dp, 1.5_dp, 2.0_dp]
        real(dp), parameter :: ps = 0.5_dp, delta = 1.0e-3_dp
        real(dp) :: ul(7), x
        type(rational_layer_t) :: lay
        integer :: i, k, st(3)

        do i = 1, 3
            do k = 1, 3
                x = delta*ratios(i)**(k - 1)
                ul(4 - k) = sqrt(ps - x)
                ul(4 + k) = sqrt(ps + x)
            end do
            ul(4) = sqrt(ps)
            call layer_init(ul, ps, 1.5_dp, 1.0_dp, -0.2_dp, 3, delta, ratios(i), 3, &
                lay, st(i))
        end do
        call record(st(1) /= 0 .and. st(2) == 0 .and. st(3) == 0, &
            "layer fit rejects ill-conditioned ring ratio 1.25 (|w|_1 > layer_wmax)")
        if (layer_wmax <= 0.0_dp) error stop
    end subroutine layer_conditioning

    subroutine layer_ring_cap()
        !! Ratio 3 with a rational at small psin: the uniform-spacing cap
        !! 2u/(nrad (ratio-1)) = 0.0073 lies below the third fit ring
        !! delta ratio^2 = 9e-3, which used to be dropped (layer_init status 1).
        real(dp), parameter :: ps = 0.49_dp, delta = 1.0e-3_dp, ratio = 3.0_dp
        real(dp), allocatable :: lab(:), ul(:)
        real(dp) :: xm, t
        type(rational_layer_t) :: lay
        integer :: st, i, j

        xm = layer_xmax(ps, [0.9_dp], delta, ratio, 3, 96)
        call layer_labels(ps, delta, ratio, xm, 1.0e-6_dp, 0.999_dp, lab)
        ul = sqrt([lab, ps])
        do i = 2, size(ul)
            t = ul(i)
            j = i - 1
            do while (j >= 1)
                if (ul(j) <= t) exit
                ul(j + 1) = ul(j)
                j = j - 1
            end do
            ul(j + 1) = t
        end do
        call layer_init(ul, ps, 1.5_dp, 1.0_dp, -0.2_dp, 3, delta, ratio, 3, lay, st)
        call record(xm >= delta*ratio**2 .and. st == 0, &
            "layer rings reach the outermost fit ring at ratio 3 (layer_xmax)")
    end subroutine layer_ring_cap

    subroutine layer_jump_exact()
        !! Manufactured layer solution c0 + c1 x + (J/2)|x| + c2 x^2 + c3 x|x|
        !! on the fit rings: the joint two-sided fit returns J exactly, also
        !! with an asymmetric (c1, c3) part that would corrupt a difference of
        !! one-sided fits with a common value constraint.
        real(dp), parameter :: ps = 0.49_dp, delta = 1.0e-3_dp
        real(dp), parameter :: c0 = 0.7_dp, c1 = -3.1_dp, c2 = 40.0_dp, c3 = -25.0_dp
        real(dp), parameter :: jtrue = 0.47_dp
        integer, parameter :: nf(2) = [3, 4]
        real(dp), allocatable :: lab(:), ul(:)
        real(dp) :: xm, x, err
        complex(dp), allocatable :: phi(:, :)
        type(rational_layer_t) :: lay
        integer :: st, k, side, ir

        err = 0.0_dp
        do ir = 1, 2
            xm = layer_xmax(ps, [0.9_dp], delta, 2.0_dp, nf(ir), 96)
            call layer_labels(ps, delta, 2.0_dp, xm, 1.0e-6_dp, 0.999_dp, lab)
            ul = sqrt(sorted([lab, ps]))
            call layer_init(ul, ps, 1.5_dp, 1.0_dp, -0.2_dp, 3, delta, 2.0_dp, nf(ir), &
                lay, st)
            if (st /= 0) error stop "layer_jump_exact init"
            allocate (phi(nf(ir), 2))
            do side = 1, 2
                do k = 1, nf(ir)
                    x = ul(lay%ring(k, side))**2 - ps
                    phi(k, side) = cmplx(c0 + c1*x + 0.5_dp*jtrue*abs(x) + c2*x*x &
                                         + c3*x*abs(x), 0.0_dp, dp)
                end do
            end do
            err = max(err, abs(layer_jump(lay, phi) - jtrue))
            deallocate (phi)
        end do
        call record(err < 1.0e-9_dp, &
            "layer_jump recovers the kink of a manufactured layer solution")
    end subroutine layer_jump_exact

    subroutine layer_kink_finite_mu()
        !! Finite-beta layer exponent mu = sqrt(-D_I) /= 1/2 (Glasser, Greene
        !! and Johnson, Phys. Fluids 18, 875 (1975)): near the rational surface
        !! the outer flux is Phi = A|x|^(1/2 - mu) + B|x|^(1/2 + mu) + C x
        !! (the Frobenius solutions of (x^2 xi')' = (mu^2 - 1/4) xi, xi = Phi/x).
        !! layer_kink recovers the jump 2 (1 + eps) B delta^eps of the small
        !! solution at |x| = delta, eps = mu - 1/2, from the fit rings for every
        !! ring ratio and fit size, without and with the large solution
        !! (a1 /= 0 on both sides), against the polynomial-kink estimate that
        !! is biased by (xbar/delta)^eps.
        real(dp), parameter :: ps = 0.49_dp, delta = 1.0e-3_dp
        real(dp), parameter :: mus(3) = [0.5_dp, 0.6_dp, 0.75_dp]
        real(dp), parameter :: ratios(3) = [1.5_dp, 2.0_dp, 1.5_dp]
        integer, parameter :: nfs(3) = [3, 3, 5]
        real(dp), parameter :: bb = 0.8_dp, cc = -2.3_dp, b1 = 3.0_dp, b2 = 200.0_dp
        real(dp), parameter :: c2 = 40.0_dp, c4 = 50.0_dp
        real(dp), allocatable :: lab(:), ul(:)
        real(dp) :: xm, x, eps, jtrue, err, errl, errb, a1(2), ax
        complex(dp), allocatable :: phi(:, :)
        type(rational_layer_t) :: lay
        integer :: st, k, side, im, ir, ia

        err = 0.0_dp
        errl = 0.0_dp
        errb = 0.0_dp
        do im = 1, 3
            eps = mus(im) - 0.5_dp
            jtrue = 2.0_dp*(1.0_dp + eps)*bb*delta**eps
            do ir = 1, 3
                xm = layer_xmax(ps, [0.9_dp], delta, ratios(ir), nfs(ir), 96)
                call layer_labels(ps, delta, ratios(ir), xm, 1.0e-6_dp, 0.999_dp, lab)
                ul = sqrt(sorted([lab, ps]))
                call layer_init(ul, ps, 1.5_dp, 1.0_dp, -mus(im)**2, 3, delta, &
                    ratios(ir), nfs(ir), lay, st)
                if (st /= 0) error stop "layer_kink_finite_mu init"
                allocate (phi(nfs(ir), 2))
                do ia = 1, 2
                    ! ia = 1: small solution with analytic corrections; ia = 2:
                    ! exact Frobenius span with different a1 on the two sides
                    a1 = 0.0_dp
                    if (ia == 2) a1 = [0.3_dp, -0.2_dp]
                    do side = 1, 2
                        do k = 1, nfs(ir)
                            x = ul(lay%ring(k, side))**2 - ps
                            ax = abs(x)
                            if (ia == 1) then
                                phi(k, side) = cmplx(cc*x &
                                    + bb*ax**(1.0_dp + eps)*(1.0_dp + b1*x + b2*x*x) &
                                    + c2*x*x + c4*x**3, 0.0_dp, dp)
                            else
                                phi(k, side) = cmplx(a1(side)*(ax/delta)**(-eps) &
                                    + cc*x + bb*ax**(1.0_dp + eps), 0.0_dp, dp)
                            end if
                        end do
                    end do
                    if (ia == 1) then
                        err = max(err, abs(layer_kink(lay, phi) - jtrue)/jtrue)
                        if (im > 1) errb = max(errb, abs(layer_jump(lay, phi) &
                            *delta**eps - jtrue)/jtrue)
                    else
                        errl = max(errl, abs(layer_kink(lay, phi, .true.) - jtrue)/jtrue)
                        ! the joint fit carries the large solution itself
                        errl = max(errl, abs(layer_kink(lay, phi) - jtrue)/jtrue)
                    end if
                end do
                deallocate (phi)
            end do
        end do
        print '(a,3es10.2)', "  layer_kink finite mu: error small / with a1 / "// &
            "polynomial kink ", err, errl, errb
        call record(err < 2.0e-2_dp .and. errl < 2.0e-2_dp, &
            "layer_kink recovers the small-solution jump for mu = 0.5, 0.6, 0.75 "// &
            "(ratio 1.5/2, nfit 3/5, also with a large solution)")
    end subroutine layer_kink_finite_mu

    subroutine layer_closure_finite_mu()
        !! Ideal edge-sheet closure for mu > 1/2 on an exact outer model.  Each
        !! side carries Phi_s = A_s |x|^-eps + C_s |x|^(1+eps) + sigma g |x|
        !! + f x^2 (Frobenius span plus the smooth odd x and even x^2 parts),
        !! with Phi_s(|x| = xb) prescribed; the closure a1_s = 0 (layer_large)
        !! fixes A_s.  Oracle (ideal: A = 0): the net layer current (outer
        !! slope jump) is (1 + eps)(C_1 + C_2) delta^eps + 4 f delta, so its
        !! fitted delta exponent must be eps, and the extracted jump at the
        !! spot must converge to (1 + eps)(C_1 + C_2) xs^eps as delta -> 0.
        real(dp), parameter :: ps = 0.49_dp, xb = 0.3_dp, xs = 5.0e-3_dp
        real(dp), parameter :: mus(3) = [0.553_dp, 0.6_dp, 0.75_dp]
        real(dp), parameter :: ratios(2) = [1.5_dp, 2.0_dp]
        real(dp), parameter :: deltas(4) = [1.0e-3_dp, 5.0e-4_dp, 2.5e-4_dp, 1.25e-4_dp]
        real(dp), parameter :: pb(2) = [0.7_dp, 1.3_dp], gp = -2.0_dp, fp = 1.0_dp
        integer, parameter :: nfit = 5
        real(dp), allocatable :: lab(:), ul(:)
        real(dp) :: eps, am(4, 4), rhs(4), a(nfit, 2), cex(2), sinf, net(4), sk(4)
        real(dp) :: lx(4), ly(4), slope, errp, errs, sg, w
        complex(dp) :: phi(nfit, 2), a1(2)
        type(rational_layer_t) :: lay
        integer :: im, ir, id, st, side, s2, k, ipiv(4)

        errp = 0.0_dp
        errs = 0.0_dp
        do im = 1, 3
            eps = mus(im) - 0.5_dp
            do side = 1, 2
                sg = real(2*side - 3, dp)
                cex(side) = (pb(side) - sg*gp*xb - fp*xb**2)/xb**(1 + eps)
            end do
            sinf = (1 + eps)*sum(cex)*xs**eps
            do ir = 1, 2
                do id = 1, 4
                    call layer_labels(ps, deltas(id), ratios(ir), &
                        layer_xmax(ps, [0.9_dp], deltas(id), ratios(ir), nfit, 96), &
                        1.0e-6_dp, 0.999_dp, lab)
                    ul = sqrt(sorted([lab, ps]))
                    call layer_init(ul, ps, 1.5_dp, 1.0_dp, -mus(im)**2, 3, deltas(id), &
                        ratios(ir), nfit, lay, st)
                    if (st /= 0) error stop "layer_closure_finite_mu init"
                    a = lay%t*deltas(id)
                    ! unknowns A_1, A_2, C_1, C_2: boundary rows, closure rows
                    am = 0.0_dp
                    do side = 1, 2
                        sg = real(2*side - 3, dp)
                        am(side, side) = xb**(-eps)
                        am(side, 2 + side) = xb**(1 + eps)
                        rhs(side) = pb(side) - sg*gp*xb - fp*xb**2
                        rhs(2 + side) = 0.0_dp
                        do s2 = 1, 2
                            sg = real(2*s2 - 3, dp)
                            do k = 1, nfit
                                w = lay%pa1((s2 - 1)*nfit + k, side)
                                am(2 + side, s2) = am(2 + side, s2) + w*a(k, s2)**(-eps)
                                am(2 + side, 2 + s2) = am(2 + side, 2 + s2) &
                                    + w*a(k, s2)**(1 + eps)
                                rhs(2 + side) = rhs(2 + side) &
                                    - w*(sg*gp*a(k, s2) + fp*a(k, s2)**2)
                            end do
                        end do
                    end do
                    call dgesv(4, 1, am, 4, ipiv, rhs, 4, st)
                    if (st /= 0) error stop "layer_closure_finite_mu solve"
                    do side = 1, 2
                        sg = real(2*side - 3, dp)
                        phi(:, side) = cmplx(rhs(side)*a(:, side)**(-eps) &
                            + rhs(2 + side)*a(:, side)**(1 + eps) + sg*gp*a(:, side) &
                            + fp*a(:, side)**2, 0.0_dp, dp)
                    end do
                    a1 = layer_large(lay, phi)
                    if (maxval(abs(a1)) > 1.0e-10_dp) error stop &
                        "layer_closure_finite_mu: closure a1 not zero"
                    ! net current: outer slope jump d Phi/d|x| summed over both sides
                    net(id) = sum(-eps*rhs(1:2)*deltas(id)**(-1 - eps) &
                        + (1 + eps)*rhs(3:4)*deltas(id)**eps) + 4*fp*deltas(id)
                    sk(id) = real(layer_kink(lay, phi), dp)*(xs/deltas(id))**eps
                end do
                lx = log(deltas)
                ly = log(abs(net))
                slope = sum((lx - sum(lx)/4)*(ly - sum(ly)/4))/sum((lx - sum(lx)/4)**2)
                errp = max(errp, abs(slope - eps))
                errs = max(errs, abs(sk(4) - sinf)/abs(sinf))
                if (abs(sk(4) - sinf) > abs(sk(1) - sinf) + 1.0e-12_dp*abs(sinf)) &
                    errs = huge(1.0_dp)
            end do
        end do
        print '(a,2es10.2)', "  closure finite mu: |net-current exponent - eps|, "// &
            "spot-jump error ", errp, errs
        call record(errp < 2.0e-2_dp .and. errs < 1.0e-2_dp, &
            "finite-mu edge-sheet closure: net layer current ~ delta^eps and the "// &
            "spot jump converges (mu 0.553/0.6/0.75, ratio 1.5/2, nfit 5)")
    end subroutine layer_closure_finite_mu

    pure function sorted(v) result(w)
        real(dp), intent(in) :: v(:)
        real(dp) :: w(size(v)), t
        integer :: a, b
        w = v
        do a = 2, size(w)
            t = w(a)
            b = a - 1
            do while (b >= 1)
                if (w(b) <= t) exit
                w(b + 1) = w(b)
                b = b - 1
            end do
            w(b + 1) = t
        end do
    end function sorted

    subroutine fe_ladder(name, r0s, tol)
        !! FEEC closure Delta' at the layer half-widths delta = 2, 1, 0.5e-3
        !! (psin), extrapolated with the Frobenius remainder of the large
        !! solution inside the layer, c1 delta ln delta + c2 delta (the jump
        !! of L' = K (ln|x| + 1) + K^2 x ln|x| + ... across the layer), on the
        !! torus R0/a = r0s, against the cylindrical ODE.
        character(len=*), intent(in) :: name
        real(dp), intent(in) :: r0s(:), tol
        real(dp), parameter :: deltas(3) = [2.0e-3_dp, 1.0e-3_dp, 5.0e-4_dp]
        real(dp) :: dfe(3), dode, rs, err(size(r0s)), am(3, 3), rhs(3)
        integer :: i, k, ipiv(3), info
        nmode = 1
        mres = 2
        q0 = 1.2_dp
        call make_tables()
        rs = r_of_q(real(mres, dp)/nmode)
        dode = ode_delta(rs)
        print '(a,a,a,f12.6)', "  ", name, ": ODE r_s Delta' =", dode
        do i = 1, size(r0s)
            r0 = r0s(i)
            call make_equilibrium()
            do k = 1, 3
                call fe_delta(deltas(k), 1.5_dp, rs, dfe(k))
                am(k, :) = [1.0_dp, deltas(k)*log(deltas(k)), deltas(k)]
            end do
            rhs = dfe
            call dgesv(3, 1, am, 3, ipiv, rhs, 3, info)
            err(i) = abs(rhs(1) - dode)/abs(dode)
            print '(a,f5.1,a,3f10.5,a,f10.5,a,es10.2)', "    R0/a =", r0, &
                ": FEEC closure at delta 2, 1, 0.5e-3:", dfe, ", extrapolated", &
                rhs(1), &
                "  rel. error", err(i)
        end do
        call record(err(size(r0s)) < tol, &
            "closure Delta' against the cylinder ODE ("// &
            name//")")
        if (size(r0s) > 1) call record(err(size(r0s)) < err(1), &
            "closure Delta' approaches the cylinder as a/R0 -> 0 ("//name//")")
    end subroutine fe_ladder

    subroutine fe_delta(delta, ratio, rs_cyl, dprime)
        !! r_s Delta' of q = mres/nmode from the FEEC response with the
        !! Frobenius closure in Delta' mode.
        real(dp), intent(in) :: delta, ratio, rs_cyl
        real(dp), intent(out) :: dprime
        real(dp) :: gmres_tol
        integer, parameter :: np = 256
        integer :: nfit
        integer :: nrad
        type(flux_mesh_t) :: fm
        type(meridian_feec_t) :: fs
        type(meridian_point_basis_t) :: pb
        type(ampere_n_t) :: amp
        type(ideal_options_t) :: opt
        type(sample_ring_t), allocatable :: sr(:)
        type(rational_layer_t), allocatable :: lays(:)
        type(layer_closure_t) :: cl
        type(surface_t) :: s
        real(dp), allocatable :: ul(:), extra(:), lab(:), psin_s(:)
        complex(dp), allocatable :: a(:), f(:), b2(:), b2cl(:, :), lam(:), phit(:)
        complex(dp), allocatable :: phe(:)
        complex(dp) :: cf(3, 2)
        integer, allocatable :: mr(:)
        real(dp) :: qsgn, qt, di, qa, dq, res, dpsindr, x, dfit
        complex(dp) :: junit
        integer :: nr, j, l, st, info, its, nrat, k, base, itest
        logical :: found

        nrad = nrad_fe
        nfit = nfit_fe
        block
            character(len=4) :: rsh
            call get_environment_variable('KIN6D_RES_SHAPE', rsh)
            if (rsh(1:1) == '0') opt%strip_shape = .false.
        end block
        ! rationals q = m/n inside the plasma (FE equilibrium)
        call surface_build(g, 0.25_dp, 256, s, st)
        qsgn = sign(1.0_dp, s%q)
        allocate (psin_s(0), mr(0))
        do k = 1, 20
            qt = real(k, dp)/nmode
            found = find_surface_q(g, qsgn*qt, 0.02_dp, 0.98_dp, 256, 1.0e-13_dp, &
                x, st)
            if (.not. found) cycle
            psin_s = [psin_s, x]
            mr = [mr, nint(qsgn*qt*nmode)]
            if (k == mres) itest = size(mr)
        end do
        nrat = size(mr)
        extra = sqrt(psin_s)
        do l = 1, nrat
            qa = 2*sqrt(psin_s(l))/(nrad*(ratio - 1))
            do k = 1, nrat
                if (k /= l) qa = min(qa, 0.45_dp*abs(psin_s(k) - psin_s(l)))
            end do
            call layer_labels(psin_s(l), delta, ratio, qa, (0.5_dp/nrad)**2, &
                (1 - 0.5_dp/nrad)**2, lab)
            extra = [extra, sqrt(lab)]
        end do
        call ring_labels_merge(nrad, extra, 0.3_dp, ul)
        call ring_labels_smooth(ul, 2*ratio)
        nr = size(ul)
        call flux_mesh_build(ul, 8, 256, 2, polar_map, fm, st)
        if (st /= 0) error stop "mesh"
        call meridian_feec_build(fm%mm, 2, fs, st)
        call ampere_n_build(fs, nmode, amp, st)
        if (st /= 0) error stop "ampere"
        call meridian_point_basis_init(fs, pb, st)
        allocate (lays(nrat))
        do l = 1, nrat
            call mercier_di(g, psin_s(l), 512, di, qa, dq, st)
            if (st /= 0) error stop "Mercier"
            qa = 2*sqrt(psin_s(l))/(nrad*(ratio - 1))
            do k = 1, nrat
                if (k /= l) qa = min(qa, 0.45_dp*abs(psin_s(k) - psin_s(l)))
            end do
            call layer_init(ul, psin_s(l), real(abs(mr(l)), dp)/nmode, dq, di, mr(l), &
                delta, ratio, nfit, lays(l), st)
            if (st /= 0) error stop "layer"
        end do
        opt%res_m = mr(1:nrat)
        opt%res_eps = lays(1:nrat)%eps
        call layer_closure(g, ul, lays, 1024, 2, cl, st)
        if (st /= 0) error stop "closure"
        cl%target(2*itest - 1:2*itest) = (1.0_dp, 0.0_dp)
        allocate (a(amp%nn), f(amp%nn), b2cl(fs%sp(2)%ndof, cl%ns))
        do j = 1, cl%ns
            f = (0.0_dp, 0.0_dp)
            call ring_sheet_load(fs, pb, fm, cl%sheet_ring(j), cl%sheet_ra(j), nmode, &
                cl%sheet_m(j), (1.0_dp, 0.0_dp), f)
            call ampere_n_solve(amp, f, a, st)
            call ampere_n_curl(amp, a, b2)
            b2cl(:, j) = b2
        end do
        allocate (sr(nr))
        do j = 2, nr
            call sample_ring_build(g, opt, fm, fs, pb, 0.5_dp*(ul(j - 1) + ul(j)), j, &
                np, sr(j - 1), st)
            if (st /= 0) error stop "half ring"
            do l = 1, nrat
                if (sr(j - 1)%u > ul(lays(l)%ring(1, 1)) .and. &
                    sr(j - 1)%u < ul(lays(l)%ring(1, 2))) sr(j - 1)%active = .false.
            end do
        end do
        call sample_ring_build(g, opt, fm, fs, pb, 1.0_dp, nr, np, sr(nr), st)
        allocate (lam(cl%ns), phit(cl%nf), phe(cl%nf))
        phe = (0.0_dp, 0.0_dp)
        a = (0.0_dp, 0.0_dp)
        ! GMRES wraps K^-1.  The MUMPS Cholesky inverts K to roundoff; the complex
        ! SuperLU fallback leaves a relative residual of ~1e-7 at R0/a = 40, which
        ! stalls GMRES at 2e-6 .. 1.5e-5 (400 iterations, measured; plain or
        ! GMRES-polished refinement did not pay for itself), so allow 3e-5 there.
        gmres_tol = 1.0e-6_dp
        if (.not. gs_direct_available()) gmres_tol = 3.0e-5_dp
        call ideal_response(fs, pb, fm, g, amp, opt, nmode, sr, cl, b2cl, phe, a, &
            lam, phit, gmres_tol, gmres_max, info, its, res)
        call record(info == 0, "GMRES converged")
        base = 0
        do l = 1, itest - 1
            base = base + 2*lays(l)%nfit
        end do
        cf = layer_coefficients(lays(itest), reshape(phit(base + 1:base + 2*nfit), &
            [nfit, 2]))
        ! jump of d Phi/d psin per unit sheet: vacuum field of the unit sheet
        ! on the rational ring, one-sided slopes on the fit rings
        f = (0.0_dp, 0.0_dp)
        call ring_sheet_load(fs, pb, fm, lays(itest)%ring_rat, &
            cl%sheet_ra(2*itest - 1), &
            nmode, mr(itest), (1.0_dp, 0.0_dp), f)
        call ampere_n_solve(amp, f, a, st)
        call ampere_n_curl(amp, a, b2)
        phe = closure_fluxes(fs, pb, fm, cl, nmode, b2)
        junit = layer_regular_slope(lays(itest), &
            reshape(phe(base + 1:base + 2*nfit), [nfit, 2]))
        ! Delta' in psin: (a2+ + a2-)/(delta a1); to r with dpsin/dr at r_s
        dpsindr = 2.0_dp*b0/q0*ival(rs_cyl)/rs_cyl/psi_of_r(1.0_dp)
        dfit = real((cf(2, 1) + cf(2, 2))/(lays(itest)%delta*cf(1, 1)), dp) &
            *dpsindr*rs_cyl
        dprime = real((lam(2*itest - 1) + lam(2*itest))*junit/cf(1, 1), dp) &
            *dpsindr*rs_cyl
        print '(a,f12.6,a,4es12.4,a,2es12.4)', "      Delta' from the fitted slopes", &
            dfit, &
            " lam+-", lam(2*itest - 1), lam(2*itest), " junit", junit
        print '(a,i0,a,i0,a,es9.2,a,2es10.2,a,es10.2)', "      rings ", nr, &
            ", GMRES ", &
            its, " res", res, &
            ", a1 in/out ", abs(cf(1, 1)), abs(cf(1, 2)), ", log-term asymmetry ", &
            abs(cf(3, 1) + cf(3, 2))/max(abs(cf(3, 1)), tiny(1.0_dp))
    end subroutine fe_delta

    ! ============================================== driven ideal shielding

    subroutine ode_shield(rs, s_edge, edge, dp_edge)
        !! Ideal shielding of the vacuum drive (r/a)^m in the cylinder with
        !! the edge current: the outer solution u(r_s) = 0, u'(r_s) = 1, the
        !! edge jump psi'(1-) = psi'(1+) + kappa psi(1), kappa = j(1)/(I(1)
        !! (1 - n q(1)/m)), matched to the drive plus r^-m outside (no wall).
        !! s_edge = r_s [psi']/psi_vac(r_s), edge = psi(1)/psi_vac(1)
        !! (cyl_response.py S_edge_current, edge_total_over_applied), and
        !! dp_edge = r_s Delta' with the same edge jump (psi'(1-) = (kappa -
        !! m) psi(1) for the free tearing solution).
        real(dp), intent(in) :: rs
        real(dp), intent(out) :: s_edge, edge, dp_edge
        integer, parameter :: ne = 6
        real(dp) :: eta(ne), sv(ne), ev(ne), y(2), am(ne, 4), am2(ne, 4), kap
        real(dp) :: work(64), le, rhs(ne, 3), yl(2)
        integer :: i, info
        kap = jprof(1.0_dp)/(ival(1.0_dp)*(1.0_dp - nmode*qcyl(1.0_dp)/mres))
        do i = 1, ne
            eta(i) = 1.0e-3_dp*0.25_dp**(i - 1)
            y = [eta(i), 1.0_dp]
            call integrate(rs + eta(i), 1.0_dp, y)
            ev(i) = 2*mres/(mres - kap + y(2)/y(1))
            sv(i) = rs*(ev(i)/y(1))/rs**mres
            yl = [1.0e-3_dp**mres, mres*1.0e-3_dp**(mres - 1)]
            call integrate(1.0e-3_dp, rs - eta(i), yl)
            y = [1.0_dp, kap - mres]
            call integrate(1.0_dp, rs + eta(i), y)
            rhs(i, 3) = rs*(y(2)/y(1) - yl(2)/yl(1))
            le = log(eta(i))
            am(i, :) = [1.0_dp, eta(i)*le*le, eta(i)*le, eta(i)]
        end do
        rhs(:, 1) = sv
        rhs(:, 2) = ev
        am2 = am
        call dgels('N', ne, 4, 3, am2, ne, rhs, ne, work, size(work), info)
        s_edge = rhs(1, 1)
        edge = rhs(1, 2)
        dp_edge = rhs(1, 3)
    end subroutine ode_shield

    subroutine shield_cylinder()
        !! Driven ideal shielding of q = 3/1 in the GPEC lar default cylinder
        !! (q = q0 (1 + r^2), q0 = 2.2, j(1) /= 0: edge current) at R0/a = 40:
        !! the jump of d Phi/d psin across the layer over the vacuum resonant
        !! flux, r_s [psi']/psi_vac, against the cylinder ODE with the edge
        !! current (cyl_response.py S_edge_current = 2.2967), at the ring
        !! ratios 2 and 3 of the layer fit.
        real(dp) :: rs, sref, eref, dpe, sfe(2), sl(2), efe(2), err
        integer :: i
        prof = 4
        q0 = 2.2_dp
        nmode = 1
        mres = 3
        r0 = 40.0_dp
        call make_tables()
        rs = r_of_q(real(mres, dp)/nmode)
        call ode_shield(rs, sref, eref, dpe)
        print '(a,f10.6,a,f10.6,a,f10.6)', "  shielding q = 3/1, q0 = 2.2: r_s =", rs, &
            ", ODE r_s [psi']/psi_vac =", sref, ", edge total/applied =", eref
        call make_equilibrium()
        err = 0.0_dp
        do i = 1, 2
            call fe_shield(1.0e-3_dp, real(i + 1, dp), 32, rs, sfe(i), sl(i), efe(i))
            err = max(err, abs(sfe(i) - sref)/sref, abs(sl(i) - sref)/sref)
            print '(a,f4.1,a,3f10.6)', "    ratio", real(i + 1, dp), &
                ": FEEC S (fit, sheets), edge ratio", sfe(i), sl(i), efe(i)
        end do
        call record(err < 2.0e-2_dp .and. abs(efe(1) - eref) < 2.0e-2_dp*eref, &
            "driven ideal shielding with edge current against the cylinder ODE")
    end subroutine shield_cylinder

    subroutine shield_study(args)
        !! KIN6D_CLOSURE_STUDY="shield q0 m r0 nrad delta ratio bs mc [ReD ImD
        !! [nfit]]" (GPEC lar default profile, R0/a = r0, boundary sheet bs and
        !! moved edge current mc as 0/1).  KIN6D_SHIELD_DPRIME=1 replaces the
        !! drive by unit large-solution targets (Delta' mode, ReD >= 900) and
        !! prints r_s Delta' from the fitted slopes and from the sheets;
        !! KIN6D_SHIELD_TOL sets the GMRES tolerance (default 1e-6)."
        !! (with r_s Delta_layer = ReD + i ImD the layer is resistive)
        character(len=*), intent(in) :: args
        real(dp) :: rs, sref, eref, dpe, delta, ratio, sfe, sl, efe, dd(2)
        complex(dp) :: pr, pex
        integer :: nrad, ibs, imc, ios
        read (args, *) q0, mres, r0, nrad, delta, ratio, ibs, imc
        block
            character(len=16) :: tv
            call get_environment_variable('KIN6D_SHIELD_TOL', tv)
            if (len_trim(tv) > 0) read (tv, *) shield_tol
        end block
        read (args, *, iostat=ios) q0, mres, r0, nrad, delta, ratio, ibs, imc, dd, &
            nfit_fe
        if (ios /= 0) nfit_fe = 3
        read (args, *, iostat=ios) q0, mres, r0, nrad, delta, ratio, ibs, imc, dd
        prof = 4
        nmode = 1
        call make_tables()
        rs = r_of_q(real(mres, dp)/nmode)
        call ode_shield(rs, sref, eref, dpe)
        call make_equilibrium()
        if (ios == 0 .and. dd(1) < 900.0_dp) then
            call fe_shield(delta, ratio, nrad, rs, sfe, sl, efe, ibs == 1, imc == 1, &
                cmplx(dd(1), dd(2), dp), pr)
            pex = sref/(cmplx(dd(1), dd(2), dp) - dpe)
            print '(a,2f9.4,a,2f11.6,a,2f11.6,a,f10.6)', 'resist D', dd, &
                ' FE A/psi_vac', pr, ' exact S/(D - Delta_p)', pex, ' rs Delta_p', dpe
        else
            call fe_shield(delta, ratio, nrad, rs, sfe, sl, efe, ibs == 1, imc == 1)
        end if
        print '(a,f5.2,i3,f6.1,i5,es10.2,f5.2,2i2,a,3f11.6,a,2f11.6)', 'shield', &
            q0, mres, r0, nrad, delta, ratio, ibs, imc, ' FE S fit/sheets edge', &
            sfe, sl, efe, ' ODE', sref, eref
    end subroutine shield_study

    subroutine fe_shield(delta, ratio, nrad, rs_cyl, s_fit, s_sheet, edge, bsheet, &
            moved, dlay, prec, dprime, dpr)
        !! Driven ideal response of the torus with the Frobenius closure
        !! (a1 = 0 on both sides of every rational, the gs_response_3d
        !! 'frobenius' path): s_fit = r_s [psi']/psi_vac from layer_jump of
        !! the total fluxes on the fit rings, s_sheet the same from the edge
        !! sheet amplitudes (lam_in + lam_out) jump_unit, edge the total over
        !! applied resonant flux on the boundary ring.  With dlay the layer of
        !! q = mres/nmode is resistive (resistive_closure) with the prescribed
        !! layer response r_s Delta_layer = dlay (the other layers stay ideal)
        !! and prec returns the reconnected flux A/psi_vac (a1 average).
        !! dprime (or KIN6D_SHIELD_DPRIME) replaces the drive by unit large-
        !! solution targets; dpr returns r_s Delta' from the edge sheets.
        use kin6d_mhd_response_n, only: ring_angles_t, ring_angles_set, &
            ring_flux_harmonic
        real(dp), intent(in) :: delta, ratio, rs_cyl
        integer, intent(in) :: nrad
        real(dp), intent(out) :: s_fit, s_sheet, edge
        logical, intent(in), optional :: bsheet, moved
        complex(dp), intent(in), optional :: dlay
        complex(dp), intent(out), optional :: prec, dpr
        logical, intent(in), optional :: dprime
        integer, parameter :: np = 256
        type(flux_mesh_t) :: fm
        type(meridian_feec_t) :: fs
        type(meridian_point_basis_t) :: pb
        type(ampere_n_t) :: amp
        type(ideal_options_t) :: opt
        type(sample_ring_t), allocatable :: sr(:)
        type(rational_layer_t), allocatable :: lays(:)
        type(layer_closure_t) :: cl
        type(surface_t) :: s
        type(ring_angles_t) :: rac
        real(dp), allocatable :: ul(:), extra(:), lab(:), psin_s(:)
        complex(dp), allocatable :: a(:), f(:), b2(:), b2cl(:, :), lam(:), phit(:)
        complex(dp), allocatable :: phe(:), phu(:), jun(:), dll(:)
        complex(dp) :: ph(nfit_fe, 2), cf(3, 2), junit, jfit, phv, pht1, phb
        integer, allocatable :: mr(:)
        real(dp) :: qsgn, qt, di, qa, dq, res, dpsindr, xx
        integer :: nr, j, l, st, info, its, nrat, k, base, itest, mk, nfit
        logical :: found

        nfit = nfit_fe
        block
            character(len=4) :: rsh
            call get_environment_variable('KIN6D_RES_SHAPE', rsh)
            if (rsh(1:1) == '0') opt%strip_shape = .false.
        end block
        if (present(bsheet)) opt%boundary_sheet = bsheet
        if (present(moved)) opt%moved_current = moved
        call surface_build(g, 0.25_dp, 256, s, st)
        qsgn = sign(1.0_dp, s%q)
        msgn = qsgn
        allocate (psin_s(0), mr(0))
        do k = 1, 20
            qt = real(k, dp)/nmode
            found = find_surface_q(g, qsgn*qt, 0.02_dp, 0.98_dp, 256, 1.0e-13_dp, &
                xx, st)
            if (.not. found) cycle
            psin_s = [psin_s, xx]
            mr = [mr, nint(qsgn*qt*nmode)]
            if (k == mres) itest = size(mr)
        end do
        nrat = size(mr)
        extra = sqrt(psin_s)
        do l = 1, nrat
            qa = layer_xmax(psin_s(l), [psin_s(:l - 1), psin_s(l + 1:)], delta, ratio, &
                nfit, nrad)
            call layer_labels(psin_s(l), delta, ratio, qa, (0.5_dp/nrad)**2, &
                1 - 2.0e-4_dp, lab)
            extra = [extra, sqrt(lab)]
        end do
        call ring_labels_merge(nrad, extra, 0.3_dp, ul)
        call ring_labels_smooth(ul, 2*ratio)
        nr = size(ul)
        call flux_mesh_build(ul, 8, 256, 2, polar_map, fm, st)
        if (st /= 0) error stop "mesh"
        call meridian_feec_build(fm%mm, 2, fs, st)
        call ampere_n_build(fs, nmode, amp, st)
        if (st /= 0) error stop "ampere"
        call meridian_point_basis_init(fs, pb, st)
        allocate (lays(nrat))
        do l = 1, nrat
            call mercier_di(g, psin_s(l), 512, di, qa, dq, st)
            if (st /= 0) error stop "Mercier"
            call layer_init(ul, psin_s(l), real(abs(mr(l)), dp)/nmode, dq, di, mr(l), &
                delta, ratio, nfit, lays(l), st)
            if (st /= 0) error stop "layer"
        end do
        opt%res_m = mr(1:nrat)
        opt%res_eps = lays(1:nrat)%eps
        call layer_closure(g, ul, lays, 1024, 2, cl, st)
        if (st /= 0) error stop "closure"
        allocate (a(amp%nn), f(amp%nn), b2cl(fs%sp(2)%ndof, cl%ns))
        do j = 1, cl%ns
            f = (0.0_dp, 0.0_dp)
            call ring_sheet_load(fs, pb, fm, cl%sheet_ring(j), cl%sheet_ra(j), nmode, &
                cl%sheet_m(j), (1.0_dp, 0.0_dp), f)
            call ampere_n_solve(amp, f, a, st)
            call ampere_n_curl(amp, a, b2)
            b2cl(:, j) = b2
        end do
        dpsindr = 2.0_dp*b0/q0*ival(rs_cyl)/rs_cyl/psi_of_r(1.0_dp)
        ! jump of d Phi/d psin per unit rational-ring sheet of every layer
        ! (gs_response_3d jump_unit)
        allocate (jun(nrat), dll(nrat))
        base = 0
        do l = 1, nrat
            f = (0.0_dp, 0.0_dp)
            call surface_build(g, ul(lays(l)%ring_rat)**2, 1024, s, st)
            call ring_angles_set(s%q, s%theta_s(:, 2), s%nu, rac)
            call ring_sheet_load(fs, pb, fm, lays(l)%ring_rat, rac, nmode, mr(l), &
                (1.0_dp, 0.0_dp), f)
            call ampere_n_solve(amp, f, a, st)
            call ampere_n_curl(amp, a, b2)
            phu = closure_fluxes(fs, pb, fm, cl, nmode, b2)
            jun(l) = layer_jump(lays(l), reshape(phu(base + 1:base + 2*nfit), &
                [nfit, 2]))
            base = base + 2*nfit
        end do
        if (present(dlay)) then
            dll = (1.0e12_dp, 0.0_dp)
            dll(itest) = dlay/(rs_cyl*dpsindr)
            call resistive_closure(lays, dll, jun, cl)
        end if
        allocate (sr(nr))
        do j = 2, nr
            call sample_ring_build(g, opt, fm, fs, pb, 0.5_dp*(ul(j - 1) + ul(j)), j, &
                np, sr(j - 1), st)
            if (st /= 0) error stop "half ring"
            do l = 1, nrat
                if (sr(j - 1)%u > ul(lays(l)%ring(1, 1)) .and. &
                    sr(j - 1)%u < ul(lays(l)%ring(1, 2))) sr(j - 1)%active = .false.
            end do
        end do
        call sample_ring_build(g, opt, fm, fs, pb, 1.0_dp, nr, np, sr(nr), st)
        do j = 1, nr
            do k = 1, np
                call drive(sr(j)%x(:, k), sr(j)%bc(:, k))
            end do
        end do
        ! drive fluxes on the closure flux rings
        allocate (lam(cl%ns), phit(cl%nf), phe(cl%nf))
        do k = 1, cl%nf
            phe(k) = drive_flux(fm, cl%flux_ring(k), cl%flux_ra(k), cl%flux_m(k))
        end do
        block
            character(len=8) :: dbg
            logical :: dpm
            dpm = .false.
            if (present(dprime)) dpm = dprime
            call get_environment_variable('KIN6D_SHIELD_DPRIME', dbg)
            if (len_trim(dbg) > 0 .or. dpm) then
                phe = (0.0_dp, 0.0_dp)
                do j = 1, nr
                    sr(j)%bc = (0.0_dp, 0.0_dp)
                end do
                cl%target(2*itest - 1:2*itest) = (1.0_dp, 0.0_dp)
            end if
        end block
        a = (0.0_dp, 0.0_dp)
        call ideal_response(fs, pb, fm, g, amp, opt, nmode, sr, cl, b2cl, phe, a, &
            lam, phit, shield_tol, gmres_max, info, its, res)
        call record(info == 0, "GMRES converged")
        base = 0
        do l = 1, itest - 1
            base = base + 2*lays(l)%nfit
        end do
        ph = reshape(phit(base + 1:base + 2*nfit), [nfit, 2])
        jfit = layer_kink(lays(itest), ph)
        cf = layer_coefficients(lays(itest), ph)
        ! total over applied resonant flux on the boundary ring
        call ampere_n_curl(amp, a, b2)
        do j = 1, cl%ns
            b2 = b2 + lam(j)*b2cl(:, j)
        end do
        call surface_build(g, 1.0_dp, 1024, s, st)
        mk = nint(sign(real(mres, dp), s%q))
        call ring_angles_set(s%q, s%theta_s(:, 2), s%nu, rac)
        phb = drive_flux(fm, nr, rac, mk)
        pht1 = ring_flux_harmonic(fs, pb, fm, nr, rac, nmode, b2, mk)
        edge = real((pht1 + phb)/phb, dp)
        junit = jun(itest)
        if (abs(cl%target(2*itest)) > 0.0_dp) &
            print '(a,2f11.5,a,2f11.5)', '      DPRIME fit (a2+ + a2-)/(delta a1)', &
            (cf(2, 1) + cf(2, 2))/(delta*0.5_dp*(cf(1, 1) + cf(1, 2)))*dpsindr*rs_cyl, &
            ' sheets', (lam(2*itest - 1) + lam(2*itest))*junit &
            /(0.5_dp*(cf(1, 1) + cf(1, 2)))*dpsindr*rs_cyl
        ! vacuum resonant flux on the rational ring
        call surface_build(g, ul(lays(itest)%ring_rat)**2, 1024, s, st)
        call ring_angles_set(s%q, s%theta_s(:, 2), s%nu, rac)
        phv = drive_flux(fm, lays(itest)%ring_rat, rac, mr(itest))
        if (present(prec)) prec = 0.5_dp*(cf(1, 1) + cf(1, 2))/phv
        if (present(dpr)) dpr = (lam(2*itest - 1) + lam(2*itest))*junit &
            /(0.5_dp*(cf(1, 1) + cf(1, 2)))*dpsindr*rs_cyl
        s_fit = real(jfit/phv, dp)*dpsindr*rs_cyl
        s_sheet = real((lam(2*itest - 1) + lam(2*itest))*junit/phv, dp)*dpsindr*rs_cyl
        print '(a,i0,a,i0,a,es9.2,a,2es10.2,a,2es11.3,a,2es11.3)', "      rings ", nr, &
            ", GMRES ", its, " res", res, ", |a1/Phi| in/out ", &
            abs(cf(1, :))/abs(ph(1, :)), ", lam+- ", abs(lam(2*itest - 1)), &
            abs(lam(2*itest)), ", Im S fit/sheet ", aimag(jfit/phv)*dpsindr*rs_cyl, &
            aimag((lam(2*itest - 1) + lam(2*itest))*junit/phv)*dpsindr*rs_cyl
    end subroutine fe_shield

    complex(dp) function drive_flux(fm, ring, rang, m)
        !! Resonant flux harmonic m of the vacuum drive on mesh ring ring.
        use kin6d_mhd_response_n, only: ring_angles_t, ring_flux_from_values, &
            ring_points
        type(flux_mesh_t), intent(in) :: fm
        integer, intent(in) :: ring, m
        type(ring_angles_t), intent(in) :: rang
        real(dp), allocatable :: x(:, :), nds(:, :), th(:), tau(:), wq(:)
        integer, allocatable :: tt(:)
        complex(dp), allocatable :: bc(:)
        complex(dp) :: bb(3)
        integer :: kk
        call ring_points(fm, ring, x, nds, th, tt, tau, wq)
        allocate (bc(size(th)))
        do kk = 1, size(th)
            call drive(x(:, kk), bb)
            bc(kk) = (bb(1)*nds(1, kk) + bb(3)*nds(2, kk))*x(1, kk)
        end do
        drive_flux = ring_flux_from_values(rang, nmode, m, th, bc, wq)
    end function drive_flux

    ! ===================================================== L5 resistive layer

    subroutine resistive_suite()
        !! L5: resistive closure with rotation (kin6d_mhd_layer_n).
        real(dp) :: rs, dprime, ejump, diffs(3), errc(4), fp, dfe(5), rfe_rho(5)
        complex(dp) :: rm(4), rb(4), dl(4), rfe(5), jfe(5), cfe(5), cen(4)
        complex(dp) :: rid, jid
        type(layer_params_t) :: par
        real(dp) :: lw, eps, pv, errf(4), errl(3)
        complex(dp) :: gu
        integer :: i, st
        prof = 1
        q0 = 1.2_dp
        nmode = 1
        mres = 2
        r0 = 40.0_dp
        call make_tables()
        call layer_checks()
        ! cylinder: matched asymptotics against the full resistive BVP
        rs = r_of_q(real(mres, dp)/nmode)
        dprime = ode_delta(rs)/rs
        ejump = ode_drive_jump(rs)
        fp = cyl_fprime(rs)
        diffs = 8.8e-5_dp*[1.0_dp, 1.0_dp/16, 1.0_dp/256]
        print '(a,f10.6,a,f12.6,a,f12.6,a,es11.3)', "  cylinder q = 2 at r_s =", rs, &
            ": r_s Delta' =", dprime*rs, ", E =", ejump, ", F' =", fp
        print '(a)', "    flux at the rational over the drive, psi(r_s)/r_s^m: "// &
            "matched (outer Delta', E; layer Delta, b(0)/b0) vs full BVP"
        do i = 1, 4
            if (i == 4) gam_l = gam_l/10
            par = layer_params_t(gam_l, diffs(min(i, 3)), rho_l, 0.0_dp, fp)
            if (i == 4) par%diff = diffs(1)
            call layer_scales(par, lw, eps, pv, gu)
            dl(i) = layer_delta(par, st, cen(i))
            rm(i) = ejump/(dl(i) - dprime)/rs**mres*cen(i)
            diff_l = par%diff
            rb(i) = bvp_ratio(rs, lw)
            errc(i) = abs(rm(i) - rb(i))/abs(rb(i))
            print '(a,es9.2,a,es9.2,a,es9.2,a,f7.3,a,2es11.3,a,2es10.2,a,2es11.3,a,'// &
                '2es11.3,a,es9.2)', "    gamma/i =", aimag(gam_l), " eta/mu0 =", &
                par%diff, " L =", lw, " Q =", eps**(-2.0_dp/3), " r_s Delta_layer", &
                dl(i)*rs, " b(0)/b0", cen(i), " matched", rm(i), " BVP", rb(i), &
                " rel", errc(i)
            if (i == 4) gam_l = gam_l*10
        end do
        call record(maxval(errc) < 1.0e-2_dp, &
            "cylinder: matched resistive layer (constant- and non-constant-psi) "// &
            "against the full resistive BVP")
        ! conjugate symmetry Omega -> -Omega
        diff_l = diffs(1)
        gam_l = conjg(gam_l)
        rid = bvp_ratio(rs, (abs(gam_l)*rho_l*diff_l*1.25663706212e-6_dp/fp**2)**0.25_dp)
        gam_l = conjg(gam_l)
        call record(abs(rid - conjg(rb(1))) < 1.0e-8_dp*abs(rb(1)), &
            "cylinder BVP: Omega -> -Omega conjugates the response")
        ! FEEC torus R0/a = 40 against the cylinder BVP; ideal limit
        call make_equilibrium()
        dfe = [diffs(1), diffs(3), 8.8e-11_dp, 8.8e-11_dp, 8.8e-11_dp]
        rfe_rho = rho_l*[1.0_dp, 1.0_dp, 1.0_dp, 1.0e-2_dp, 1.0e-4_dp]
        call fe_resistive(dfe, rfe_rho, rs, rfe, jfe, cfe, rid, jid)
        do i = 1, 2
            errf(i) = abs(rfe(i)*cfe(i) - rb(2*i - 1))/abs(rb(2*i - 1))
            errf(i + 2) = abs(rfe(i)*rs**mres - ejump/(dl(2*i - 1) - dprime)) &
                /abs(ejump/(dl(2*i - 1) - dprime))
            print '(a,es9.2,a,2es11.3,a,2es11.3,a,es9.2,a,es9.2)', "    FEEC eta/mu0 =", &
                dfe(i), ": b(0)/drive", rfe(i)*cfe(i), " BVP", rb(2*i - 1), " rel", &
                errf(i), "; outer A vs matched cylinder rel", errf(i + 2)
        end do
        call record(maxval(errf(1:2)) < 3.0e-2_dp, &
            "FEEC resistive closure against the cylinder BVP (R0/a = 40)")
        do i = 1, 3
            errl(i) = abs(jfe(2 + i) - jid)/abs(jid)
            print '(a,es9.2,a,es9.2,a,2es11.3,a,2es11.3,a,es9.2,a,es9.2)', &
                "    eta/mu0 =", dfe(2 + i), " rho =", rfe_rho(2 + i), &
                ": layer jump", jfe(2 + i), " L6", jid, " rel", errl(i), &
                ", |b(0)/drive|", abs(rfe(2 + i)*cfe(2 + i))
        end do
        print '(a,f6.2)', "    order in rho:", log(errl(1)/errl(3))/log(1.0e4_dp)
        call record(errl(3) < 1.0e-3_dp .and. errl(3) < errl(2) .and. &
            errl(2) < errl(1), &
            "ideal limit: eta -> 0 gives the ideal-inertial layer; its current "// &
            "tends to the L6 closure as the inertial width -> 0")
        call reconnect_cylinder()
    end subroutine resistive_suite

    subroutine reconnect_cylinder()
        !! Resistive cylinder q = 3/1 (GPEC lar default profile, q0 = 2.2, R0/a
        !! = 40): r_s Delta' of the edge sheets with unit large-solution
        !! targets, and the fully reconnected flux A/psi_vac of a resistive
        !! layer with r_s Delta_layer = 0, against the cylinder ODE
        !! (cyl_response.py: r_s Delta' = -0.5862, P = S/(-r_s Delta') =
        !! 3.918).  The omitted layer current biases the edge-sheet value by
        !! C delta (ln delta - 1); both are extrapolated from delta = 1e-3,
        !! 5e-4 (1/P is linear in Delta').  Layer ratio 1.5, nfit 5, nrad 32.
        !! A strip-constant 1/x current (strip_shape off) gives r_s Delta'
        !! -0.667 here.  Im(A)/Re(A) is a discretisation error (exact: 0).
        real(dp), parameter :: dl(2) = [1.0e-3_dp, 5.0e-4_dp]
        real(dp) :: rs, sref, eref, dpe, sfe, sl, efe, fd(2), pref, edp, ep, tol0
        complex(dp) :: dpr(2), pr(2), dinf, pinv
        integer :: i, nfit0
        prof = 4
        q0 = 2.2_dp
        nmode = 1
        mres = 3
        r0 = 40.0_dp
        call make_tables()
        rs = r_of_q(real(mres, dp)/nmode)
        call ode_shield(rs, sref, eref, dpe)
        pref = sref/(-dpe)
        call make_equilibrium()
        nfit0 = nfit_fe
        nfit_fe = 5
        ! GMRES stagnates near 1.5e-6 at ratio 1.5, delta 5e-4 (with or
        ! without the strip shape)
        tol0 = shield_tol
        shield_tol = 1.0e-5_dp
        do i = 1, 2
            call fe_shield(dl(i), 1.5_dp, 32, rs, sfe, sl, efe, dprime=.true., &
                dpr=dpr(i))
            call fe_shield(dl(i), 1.5_dp, 32, rs, sfe, sl, efe, &
                dlay=(0.0_dp, 0.0_dp), prec=pr(i))
            fd(i) = dl(i)*(log(dl(i)) - 1)
            print '(a,es9.2,a,2f10.5,a,2f10.5)', "    delta", dl(i), ": r_s Delta'", &
                dpr(i), ", A/psi_vac", pr(i)
        end do
        nfit_fe = nfit0
        shield_tol = tol0
        dinf = dpr(2) - (dpr(1) - dpr(2))/(fd(1) - fd(2))*fd(2)
        pinv = 1/pr(2) - (1/pr(1) - 1/pr(2))/(fd(1) - fd(2))*fd(2)
        edp = abs(real(dinf, dp)/dpe - 1)
        ep = abs(real(1/pinv, dp)/pref - 1)
        print '(a,2f10.5,a,f10.5,a,2f10.5,a,f9.5)', "    delta -> 0: r_s Delta'", &
            dinf, " (ODE", dpe, "), A/psi_vac", 1/pinv, " (ODE", pref
        call record(edp < 1.0e-2_dp .and. ep < 1.0e-2_dp .and. &
            abs(aimag(1/pinv)) < 2.0e-2_dp*abs(1/pinv), &
            "resistive cylinder: delta-extrapolated Delta' and reconnected flux "// &
            "within 1 % of the ODE, |Im A| < 2 %")
    end subroutine reconnect_cylinder

    subroutine layer_checks()
        !! Scaled layer solver against the constant-psi resistive-inertial
        !! closed form, Richardson resolution, symmetry, viscous scaling.
        complex(dp) :: gs(3), d, dc, d1, d2, dv(2)
        real(dp) :: eps(3), err(3, 3)
        integer :: i, k, st
        gs = [(1.0_dp, 0.0_dp), (0.0_dp, 1.0_dp), &
            cmplx(1.0_dp, -1.0_dp, dp)/sqrt(2.0_dp)]
        eps = [1.0e2_dp, 1.0e3_dp, 1.0e4_dp]
        do k = 1, 3
            do i = 1, 3
                d = layer_scaled_delta(gs(k), eps(i), 0.0_dp, st)
                dc = layer_ri_constant_psi(gs(k), eps(i))
                err(i, k) = abs(d/dc - 1)
            end do
            print '(a,2f7.3,a,3es10.2)', "  layer g =", gs(k), &
                ": |Dhat/constant-psi - 1| at eps 1e2, 1e3, 1e4:", err(:, k)
        end do
        print '(a,f6.2)', "    order in 1/eps:", log(err(1, 1)/err(3, 1))/log(100.0_dp)
        call record(maxval(err(3, :)) < 2.0e-3_dp .and. all(err(3, :) < err(1, :)), &
            "layer: constant-psi limit 2 pi Gamma(3/4)/Gamma(1/4) g^(5/4)/eps")
        d1 = layer_scaled_delta(gs(2), 1.0_dp, 0.0_dp, st, 1000)
        d2 = layer_scaled_delta(gs(2), 1.0_dp, 0.0_dp, st, 4000)
        print '(a,2es14.6,a,es9.2)', "  layer eps = 1, g = i:", d2, &
            "  Richardson 1000 vs 4000:", abs(d1 - d2)/abs(d2)
        call record(abs(d1 - d2) < 1.0e-6_dp*abs(d2), "layer: resolution (Q = 1)")
        d1 = layer_scaled_delta(conjg(gs(2)), 1.0_dp, 0.0_dp, st, 4000)
        call record(abs(d1 - conjg(d2)) < 1.0e-10_dp*abs(d2), &
            "layer: g -> conj(g) conjugates Dhat")
        ! viscous, constant psi, inertia-free: Dhat = c_vr g P^(1/6)/eps
        do i = 1, 2
            dv(i) = layer_scaled_delta(gs(1), 1.0e4_dp, 1.0e4_dp*100.0_dp**(i - 1), st) &
                *1.0e4_dp/(1.0e4_dp*100.0_dp**(i - 1))**(1.0_dp/6)
        end do
        d = layer_scaled_delta(gs(2), 1.0e4_dp, 1.0e6_dp, st)*1.0e4_dp/1.0e6_dp**(1.0_dp/6)
        print '(a,2f10.5,a,2f10.5)', "  viscous-resistive constant: P = 1e4, 1e6:", &
            real(dv), "; g = i:", d
        call record(abs(dv(2) - dv(1)) < 2.0e-2_dp*abs(dv(2)) .and. &
            abs(d/dv(2) - (0.0_dp, 1.0_dp)) < 2.0e-2_dp, &
            "layer: viscous-resistive scaling g P^(1/6)/eps")
    end subroutine layer_checks

    real(dp) function cyl_fprime(rs)
        !! |dF/dr| at r_s, F = (m/r) B_theta (1 - n q/m), B_theta = (2 B0/(R0
        !! q0)) I/r: |F'| = (B_theta/r) n q'.
        real(dp), intent(in) :: rs
        real(dp) :: bt, dq
        bt = 2.0_dp*b0/(r0*q0)*ival(rs)/rs
        dq = qcyl(rs)*(2.0_dp/rs - jprof(rs)*rs/ival(rs))
        cyl_fprime = bt/rs*nmode*dq
    end function cyl_fprime

    real(dp) function ode_drive_jump(rs)
        !! E: slope at r_s+ of the outer solution with zero resonant flux
        !! (small solution only) and the vacuum field r^m + c r^-m outside
        !! (no wall), extrapolated in the matching distance.
        real(dp), intent(in) :: rs
        integer, parameter :: ne = 5
        real(dp) :: eta(ne), y1(2), y2(2), am(ne, 3), rhs(ne), work(64), c
        integer :: i, info
        do i = 1, ne
            eta(i) = 1.0e-3_dp*0.25_dp**(i - 1)
            y1 = [1.0_dp, real(mres, dp)]
            y2 = [1.0_dp, -real(mres, dp)]
            call integrate(1.0_dp, rs + eta(i), y1)
            call integrate(1.0_dp, rs + eta(i), y2)
            c = -y1(1)/y2(1)
            rhs(i) = y1(2) + c*y2(2)
            am(i, :) = [1.0_dp, eta(i)*log(eta(i)), eta(i)]
        end do
        call dgels('N', ne, 3, 1, am, ne, rhs, ne, work, size(work), info)
        ode_drive_jump = rhs(1)
    end function ode_drive_jump

    complex(dp) function bvp_ratio(rs, lw)
        !! psi(r_s)/r_s^m of the full cylindrical resistive response to the
        !! external field r^m (Richardson in the grid).
        real(dp), intent(in) :: rs, lw
        complex(dp) :: p1, p2
        p1 = bvp_solve(rs, lw, 1500)
        p2 = bvp_solve(rs, lw, 3000)
        bvp_ratio = (4*p2 - p1)/3/rs**mres
    end function bvp_ratio

    complex(dp) function bvp_solve(rs, lw, nside)
        !! Linear resistive reduced MHD in the periodic cylinder (uniform
        !! rho, eta/mu0 = diff_l, Doppler rate gam_l), harmonic (m, n):
        !!   D L psi - gamma psi + i F phi = 0,
        !!   gamma rho mu0 L phi - i (F L psi - (m/r) mu0 J' psi) = 0,
        !! L = d^2/dr^2 + (1/r) d/dr - m^2/r^2, regular at r = 0, vacuum
        !! outside r = 1 with the field r^m + c r^-m (psi' + m psi = 2 m),
        !! phi = -i gamma psi/F at r = 1.  Second-order finite differences
        !! on a grid graded at r_s to the layer scale lw.
        real(dp), intent(in) :: rs, lw
        integer, intent(in) :: nside
        real(dp), allocatable :: r(:)
        complex(dp), allocatable :: ab(:, :), rhs(:, :)
        integer, allocatable :: ipiv(:)
        integer, parameter :: kl = 4, ku = 4
        integer :: nn, nu, j, k, info, ncol, cols(6)
        real(dp) :: c1, c2, hm, hp, d1(3), d2(3), lap(3), f, rr, bt, mj, h1, h2
        complex(dp) :: vals(6)
        real(dp), parameter :: mu0l = 1.25663706212e-6_dp

        nn = 2*nside
        allocate (r(0:nn))
        c1 = grading_t(rs, lw/20, nside)
        c2 = grading_t(1.0_dp - rs, lw/20, nside)
        do k = 0, nside
            r(nside - k) = rs - rs*sinh(c1*k/real(nside, dp))/sinh(c1)
            r(nside + k) = rs + (1.0_dp - rs)*sinh(c2*k/real(nside, dp))/sinh(c2)
        end do
        r(0) = 0.0_dp
        nu = 2*(nn + 1)
        allocate (ab(2*kl + ku + 1, nu), rhs(nu, 1), ipiv(nu))
        ab = (0.0_dp, 0.0_dp)
        rhs = (0.0_dp, 0.0_dp)
        call setrow(ab, rhs, kl, ku, 1, [1], [(1.0_dp, 0.0_dp)], (0.0_dp, 0.0_dp))
        call setrow(ab, rhs, kl, ku, 2, [2], [(1.0_dp, 0.0_dp)], (0.0_dp, 0.0_dp))
        do j = 1, nn - 1
            hm = r(j) - r(j - 1)
            hp = r(j + 1) - r(j)
            rr = r(j)
            d1 = [-hp/(hm*(hm + hp)), (hp - hm)/(hm*hp), hm/(hp*(hm + hp))]
            d2 = [2/(hm*(hm + hp)), -2/(hm*hp), 2/(hp*(hm + hp))]
            lap = d2 + d1/rr
            lap(2) = lap(2) - mres**2/rr**2
            bt = 2.0_dp*b0/(r0*q0)*ival(rr)/rr
            f = mres*bt/rr - nmode*b0/r0
            mj = mres/rr*2.0_dp*b0/(r0*q0)*djprof(rr)
            cols = [ix(j - 1, 1), ix(j, 1), ix(j + 1, 1), ix(j, 2), 0, 0]
            vals(1:4) = [cmplx(diff_l*lap(1), 0.0_dp, dp), diff_l*lap(2) - gam_l, &
                cmplx(diff_l*lap(3), 0.0_dp, dp), cmplx(0.0_dp, f, dp)]
            call setrow(ab, rhs, kl, ku, ix(j, 1), cols(1:4), vals(1:4), (0.0_dp, 0.0_dp))
            cols = [ix(j - 1, 2), ix(j, 2), ix(j + 1, 2), ix(j - 1, 1), ix(j, 1), &
                ix(j + 1, 1)]
            vals(1:3) = gam_l*rho_l*mu0l*lap
            vals(4:6) = cmplx(0.0_dp, -f, dp)*lap
            vals(5) = vals(5) + cmplx(0.0_dp, mj, dp)
            call setrow(ab, rhs, kl, ku, ix(j, 2), cols, vals, (0.0_dp, 0.0_dp))
        end do
        ! r = 1: psi' + m psi = 2 m, phi = -i gamma psi/F
        h1 = r(nn) - r(nn - 1)
        h2 = r(nn - 1) - r(nn - 2)
        d1 = [h1/(h2*(h1 + h2)), -(h1 + h2)/(h1*h2), (2*h1 + h2)/(h1*(h1 + h2))]
        cols(1:3) = [ix(nn - 2, 1), ix(nn - 1, 1), ix(nn, 1)]
        vals(1:3) = cmplx(d1, 0.0_dp, dp)
        vals(3) = vals(3) + mres
        call setrow(ab, rhs, kl, ku, ix(nn, 1), cols(1:3), vals(1:3), cmplx(2.0_dp*mres, 0.0_dp, dp))
        bt = 2.0_dp*b0/(r0*q0)*ival(1.0_dp)
        f = mres*bt - nmode*b0/r0
        call setrow(ab, rhs, kl, ku, ix(nn, 2), [ix(nn, 2), ix(nn, 1)], &
            [(1.0_dp, 0.0_dp), cmplx(0.0_dp, 1.0_dp, dp)*gam_l/f], (0.0_dp, 0.0_dp))
        call zgbsv(nu, kl, ku, 1, ab, size(ab, 1), ipiv, rhs, nu, info)
        if (info /= 0) error stop "bvp"
        bvp_solve = rhs(ix(nside, 1), 1)
    end function bvp_solve

    pure integer function ix(jj, iv)
        integer, intent(in) :: jj, iv
        ix = 2*jj + iv
    end function ix

    subroutine setrow(ab, rhs, kl, ku, i, cc, vv, b)
        !! Row i of a banded system, equilibrated by its largest entry.
        complex(dp), intent(inout) :: ab(:, :), rhs(:, :)
        integer, intent(in) :: kl, ku, i, cc(:)
        complex(dp), intent(in) :: vv(:), b
        real(dp) :: mx
        integer :: q
        mx = maxval(abs(vv))
        do q = 1, size(cc)
            ab(kl + ku + 1 + i - cc(q), cc(q)) = ab(kl + ku + 1 + i - cc(q), cc(q)) &
                + vv(q)/mx
        end do
        rhs(i, 1) = b/mx
    end subroutine setrow

    real(dp) function grading_t(len, h0, n) result(c)
        !! c with len c/(n sinh c) = h0.
        real(dp), intent(in) :: len, h0
        integer, intent(in) :: n
        real(dp) :: lo, hi, mid
        integer :: it
        lo = 1.0e-8_dp
        hi = 300.0_dp
        do it = 1, 200
            mid = 0.5_dp*(lo + hi)
            if (len*mid/(n*sinh(mid)) > h0) then
                lo = mid
            else
                hi = mid
            end if
        end do
        c = 0.5_dp*(lo + hi)
    end function grading_t

    subroutine fe_resistive(diffs, rhos, rs_cyl, rfe, jfe, cfe, rid, jid)
        !! FEEC response to the drive r^m e^{-i m theta} with the resistive
        !! closure at q = mres/nmode (other rationals ideal) for the
        !! diffusivities diffs, and the ideal L6 closure: Phi_s/Phi_vac
        !! (reconnected over drive resonant flux) and the layer jump of d Phi/
        !! d psin.
        use kin6d_mhd_response_n, only: ring_angles_t, ring_angles_set, &
            ring_flux_from_values, ring_points
        real(dp), intent(in) :: diffs(:), rhos(:), rs_cyl
        complex(dp), intent(out) :: rfe(:), jfe(:), cfe(:), rid, jid
        integer, parameter :: np = 256
        real(dp), parameter :: delta = 1.0e-3_dp, ratio = 1.5_dp
        type(flux_mesh_t) :: fm
        type(meridian_feec_t) :: fs
        type(meridian_point_basis_t) :: pb
        type(ampere_n_t) :: amp
        type(ideal_options_t) :: opt
        type(sample_ring_t), allocatable :: sr(:)
        type(rational_layer_t), allocatable :: lays(:)
        type(layer_closure_t) :: cl0, cl
        type(surface_t) :: s
        type(ring_angles_t) :: rar
        type(layer_params_t) :: par
        real(dp), allocatable :: ul(:), extra(:), lab(:), psin_s(:)
        real(dp), allocatable :: xg(:, :), nds(:, :), th(:), tau(:), wq(:)
        integer, allocatable :: tg(:)
        complex(dp), allocatable :: a(:), f(:), b2(:), b2cl(:, :), lam(:), phit(:)
        complex(dp), allocatable :: phe(:), junit(:), dlay(:), bcv(:)
        complex(dp) :: cf(3, 2), phivac, bb(3)
        integer, allocatable :: mr(:)
        real(dp) :: qsgn, qt, di, qa, dq, res, x, vp, bth, gp
        integer :: nr, j, l, st, info, its, nrat, k, base, itest, ic, kk
        logical :: found

        call surface_build(g, 0.25_dp, 256, s, st)
        qsgn = sign(1.0_dp, s%q)
        msgn = qsgn
        allocate (psin_s(0), mr(0))
        do k = 1, 20
            qt = real(k, dp)/nmode
            found = find_surface_q(g, qsgn*qt, 0.02_dp, 0.98_dp, 256, 1.0e-13_dp, &
                x, st)
            if (.not. found) cycle
            psin_s = [psin_s, x]
            mr = [mr, nint(qsgn*qt*nmode)]
            if (k == mres) itest = size(mr)
        end do
        nrat = size(mr)
        extra = sqrt(psin_s)
        do l = 1, nrat
            qa = 2*sqrt(psin_s(l))/(nrad_fe*(ratio - 1))
            do k = 1, nrat
                if (k /= l) qa = min(qa, 0.45_dp*abs(psin_s(k) - psin_s(l)))
            end do
            call layer_labels(psin_s(l), delta, ratio, qa, (0.5_dp/nrad_fe)**2, &
                (1 - 0.5_dp/nrad_fe)**2, lab)
            extra = [extra, sqrt(lab)]
        end do
        call ring_labels_merge(nrad_fe, extra, 0.3_dp, ul)
        call ring_labels_smooth(ul, 2*ratio)
        nr = size(ul)
        call flux_mesh_build(ul, 8, 256, 2, polar_map, fm, st)
        if (st /= 0) error stop "mesh"
        call meridian_feec_build(fm%mm, 2, fs, st)
        call ampere_n_build(fs, nmode, amp, st)
        if (st /= 0) error stop "ampere"
        call meridian_point_basis_init(fs, pb, st)
        allocate (lays(nrat), junit(nrat), dlay(nrat))
        do l = 1, nrat
            call mercier_di(g, psin_s(l), 512, di, qa, dq, st)
            if (st /= 0) error stop "Mercier"
            if (l == itest) then
                call layer_metrics(g, psin_s(l), 1024, vp, bth, gp, st)
                par = layer_params_t(gam_l, diffs(1), rho_l, 0.0_dp, &
                    nmode*abs(dq)*bth*gp)
                print '(a,es12.4,a,es12.4)', "    torus F' =", par%fprime, &
                    ", cylinder F' =", cyl_fprime(rs_cyl)
            end if
            call layer_init(ul, psin_s(l), real(abs(mr(l)), dp)/nmode, dq, di, mr(l), &
                delta, ratio, nfit_fe, lays(l), st)
            if (st /= 0) error stop "layer"
        end do
        call layer_closure(g, ul, lays, 1024, 2, cl0, st)
        if (st /= 0) error stop "closure"
        allocate (a(amp%nn), f(amp%nn), b2cl(fs%sp(2)%ndof, cl0%ns))
        do j = 1, cl0%ns
            f = (0.0_dp, 0.0_dp)
            call ring_sheet_load(fs, pb, fm, cl0%sheet_ring(j), cl0%sheet_ra(j), &
                nmode, cl0%sheet_m(j), (1.0_dp, 0.0_dp), f)
            call ampere_n_solve(amp, f, a, st)
            call ampere_n_curl(amp, a, b2)
            b2cl(:, j) = b2
        end do
        ! jump of d Phi/d psin per unit sheet on each rational ring
        base = 0
        do l = 1, nrat
            f = (0.0_dp, 0.0_dp)
            call ring_sheet_load(fs, pb, fm, lays(l)%ring_rat, cl0%sheet_ra(2*l - 1), &
                nmode, mr(l), (1.0_dp, 0.0_dp), f)
            call ampere_n_solve(amp, f, a, st)
            call ampere_n_curl(amp, a, b2)
            phe = closure_fluxes(fs, pb, fm, cl0, nmode, b2)
            junit(l) = layer_regular_slope(lays(l), &
                reshape(phe(base + 1:base + 2*nfit_fe), [nfit_fe, 2]))
            base = base + 2*nfit_fe
        end do
        allocate (sr(nr))
        do j = 2, nr
            call sample_ring_build(g, opt, fm, fs, pb, 0.5_dp*(ul(j - 1) + ul(j)), j, &
                np, sr(j - 1), st)
            if (st /= 0) error stop "half ring"
            do l = 1, nrat
                if (sr(j - 1)%u > ul(lays(l)%ring(1, 1)) .and. &
                    sr(j - 1)%u < ul(lays(l)%ring(1, 2))) sr(j - 1)%active = .false.
            end do
        end do
        call sample_ring_build(g, opt, fm, fs, pb, 1.0_dp, nr, np, sr(nr), st)
        do j = 1, nr
            do k = 1, np
                call drive(sr(j)%x(:, k), sr(j)%bc(:, k))
            end do
        end do
        ! drive fluxes on the closure rings and on the rational ring
        if (allocated(phe)) deallocate (phe)
        allocate (phe(cl0%nf))
        do kk = 1, cl0%nf
            call ring_points(fm, cl0%flux_ring(kk), xg, nds, th, tg, tau, wq)
            allocate (bcv(size(th)))
            do k = 1, size(th)
                call drive(xg(:, k), bb)
                bcv(k) = (bb(1)*nds(1, k) + bb(3)*nds(2, k))*xg(1, k)
            end do
            phe(kk) = ring_flux_from_values(cl0%flux_ra(kk), nmode, cl0%flux_m(kk), &
                th, bcv, wq)
            deallocate (bcv)
        end do
        call surface_build(g, ul(lays(itest)%ring_rat)**2, 1024, s, st)
        call ring_angles_set(s%q, s%theta_s(:, 2), s%nu, rar)
        call ring_points(fm, lays(itest)%ring_rat, xg, nds, th, tg, tau, wq)
        allocate (bcv(size(th)))
        do k = 1, size(th)
            call drive(xg(:, k), bb)
            bcv(k) = (bb(1)*nds(1, k) + bb(3)*nds(2, k))*xg(1, k)
        end do
        phivac = ring_flux_from_values(rar, nmode, mr(itest), th, bcv, wq)
        base = 2*nfit_fe*(itest - 1)
        allocate (lam(cl0%ns), phit(cl0%nf))
        do ic = 0, size(diffs)
            cl = cl0
            if (ic > 0) then
                dlay = (1.0e30_dp, 0.0_dp)
                par%diff = diffs(ic)
                par%rho = rhos(ic)
                dlay(itest) = layer_delta(par, st, cfe(ic))/gp
                if (st /= 0) error stop "layer solve"
                call resistive_closure(lays, dlay, junit, cl)
            end if
            a = (0.0_dp, 0.0_dp)
            call ideal_response(fs, pb, fm, g, amp, opt, nmode, sr, cl, b2cl, phe, a, &
                lam, phit, 1.0e-7_dp, gmres_max, info, its, res)
            print '(a,i3,a,i3,a,i5,a,es9.2)', "    FEEC case", ic, ": info", info, &
                ", GMRES", its, ", residual", res
            call record(info == 0, "GMRES converged")
            cf = layer_coefficients(lays(itest), &
                reshape(phit(base + 1:base + 2*nfit_fe), [nfit_fe, 2]))
            if (ic == 0) then
                rid = 0.5_dp*(cf(1, 1) + cf(1, 2))/phivac
                jid = (lam(2*itest - 1) + lam(2*itest))*junit(itest)
            else
                rfe(ic) = 0.5_dp*(cf(1, 1) + cf(1, 2))/phivac
                jfe(ic) = (lam(2*itest - 1) + lam(2*itest))*junit(itest)
                print '(a,es9.2,a,2es11.3,a,es9.2,a,i4)', "    FEEC eta/mu0", &
                    diffs(ic), ": Delta_psin", dlay(itest), " |a1 in - out|/|a1|", &
                    abs(cf(1, 1) - cf(1, 2))/abs(cf(1, 1)), " GMRES", its
            end if
        end do
    end subroutine fe_resistive

end program test_mhd_closure
