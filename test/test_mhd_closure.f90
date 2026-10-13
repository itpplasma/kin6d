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
    !! KIN6D_CLOSURE_STUDY="prof R0 nrad delta ratio maxit nfit" runs one FEEC
    !! Delta' case (convergence studies), "edge" the edge oracle alone.
    !$ use omp_lib, only: omp_get_max_threads, omp_set_num_threads
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_ode_cash_karp, only: cash_karp_step
    use fortfem_meridian_feec_n, only: meridian_feec_t, meridian_feec_build, &
        meridian_point_basis_t, meridian_point_basis_init
    use kin6d_vacuum_surfaces, only: psi_grid_t, psi_grid_init, flux_ray_point, &
        surface_t, surface_build, find_surface_q
    use kin6d_mhd_mesh, only: flux_mesh_t, flux_mesh_build, ring_labels_merge
    use kin6d_mhd_ampere_n, only: ampere_n_t, ampere_n_build, ampere_n_solve, &
        ampere_n_curl
    use kin6d_mhd_response_n, only: ring_sheet_load, sheet_closure_t, closure_fluxes
    use kin6d_mhd_ideal_n, only: sample_ring_t, sample_ring_build, &
        ideal_options_t, ideal_response
    use kin6d_mhd_closure_n, only: rational_layer_t, mercier_di, layer_labels, &
        layer_init, layer_closure, layer_coefficients, layer_regular_slope
    implicit none

    real(dp), parameter :: pi = 3.14159265358979323846_dp
    real(dp), parameter :: b0 = 2.0_dp
    integer, parameter :: ntab = 14001
    real(dp), parameter :: rtab_max = 1.4_dp
    integer :: prof = 1                  ! 1 Wesson, 2 tanh step, 3 uniform
    real(dp) :: q0 = 1.2_dp, r1 = 0.5_dp, wstep = 0.08_dp
    integer :: nmode = 1, mres = 2
    real(dp) :: r0, msgn = 1.0_dp
    real(dp) :: tab_r(ntab), tab_i(ntab), tab_p(ntab), tab_b(ntab)
    type(psi_grid_t) :: g
    logical :: all_passed
    integer :: nrad_fe = 32, gmres_max = 400, nfit_fe = 3

    all_passed = .true.
    call study()
    !$ call omp_set_num_threads(min(4, omp_get_max_threads()))
    call step_closed_form()
    call mercier_zero_beta()
    call edge_current_jump()
    prof = 1
    call fe_ladder("Wesson j = (1 - r^2)^2", [10.0_dp, 40.0_dp], 2.0e-2_dp)
    prof = 2
    wstep = 0.08_dp
    call fe_single("tanh step, width 0.08", 40.0_dp, 1.0e-3_dp, 3.0e-2_dp)
    if (.not. all_passed) error stop 1
    print '(a)', "test_mhd_closure: all checks passed"

contains

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
    end subroutine mercier_zero_beta

    ! ----------------------------------------------------------- FE ladder

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
        type(sheet_closure_t) :: cl
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
        type(sheet_closure_t) :: cl
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
                if (k /= l) qa = min(qa, 0.4_dp*abs(psin_s(k) - psin_s(l)))
            end do
            call layer_labels(psin_s(l), delta, ratio, qa, (0.5_dp/nrad)**2, &
                (1 - 0.5_dp/nrad)**2, lab)
            extra = [extra, sqrt(lab)]
        end do
        call ring_labels_merge(nrad, extra, 0.3_dp, ul)
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
                if (k /= l) qa = min(qa, 0.4_dp*abs(psin_s(k) - psin_s(l)))
            end do
            call layer_init(ul, psin_s(l), real(abs(mr(l)), dp)/nmode, dq, di, mr(l), &
                delta, ratio, nfit, lays(l), st)
            if (st /= 0) error stop "layer"
        end do
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
        call ideal_response(fs, pb, fm, g, amp, opt, nmode, sr, cl, b2cl, phe, a, &
            lam, phit, 1.0e-6_dp, gmres_max, info, its, res)
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

end program test_mhd_closure
