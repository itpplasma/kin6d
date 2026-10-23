module test_mhd_newcomb_callbacks
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_ode_cash_karp, only: cash_karp_step
    use fortfem_meridian_feec_n, only: meridian_feec_t, meridian_feec_build, &
        meridian_point_basis_t, meridian_point_basis_init
    use kin6d_vacuum_surfaces, only: psi_grid_t, psi_grid_init, flux_ray_point, &
        surface_t, surface_build
    use kin6d_mhd_mesh, only: flux_mesh_t, flux_mesh_build
    use kin6d_mhd_ampere_n, only: ampere_n_t, ampere_n_build, ampere_n_curl, mu0
    use kin6d_mhd_response_n, only: ring_angles_t, ring_angles_set, &
        ring_flux_harmonic, ring_flux_from_values, ring_points, layer_closure_t
    use kin6d_mhd_ideal_n, only: sample_ring_t, sample_ring_build, &
        ideal_options_t, ideal_response

    use kin6d_mhd_mesh, only: polar_map_t

    implicit none
    private
    public :: run_test_mhd_newcomb
    real(dp), parameter :: pi = 3.14159265358979323846_dp
    real(dp), parameter :: b0 = 2.0_dp
    integer, parameter :: nmode = 2, mabs = 3, nrad = 32, np = 256
    integer, parameter :: nprof = 2001
    type, abstract :: callback_function_t
    contains
        procedure(callback_evaluate), deferred :: evaluate
    end type
    abstract interface
        subroutine callback_evaluate(self, t, y, dydt, ctx)
            import :: callback_function_t, dp
            class(callback_function_t), intent(in) :: self
            real(dp), intent(in) :: t, y(:)
            real(dp), intent(out) :: dydt(:)
            class(*), optional, intent(in) :: ctx
        end subroutine
    end interface

    type :: run_test_mhd_newcomb_newcomb_rhs_context_t
        real(dp), pointer :: beta => null()
        real(dp), pointer :: bzprof(:) => null()
        real(dp), pointer :: r0 => null()
        real(dp), pointer :: rprof(:) => null()
    end type

    type, extends(callback_function_t) :: run_test_mhd_newcomb_newcomb_rhs_t
        type(run_test_mhd_newcomb_newcomb_rhs_context_t), pointer :: data => null()
    contains
        procedure :: evaluate => run_test_mhd_newcomb_newcomb_rhs_binding
    end type

    type :: run_test_mhd_newcomb_polar_map_context_t
        type(psi_grid_t), pointer :: g => null()
    end type

    type, extends(polar_map_t) :: run_test_mhd_newcomb_polar_map_t
        type(run_test_mhd_newcomb_polar_map_context_t), pointer :: data => null()
    contains
        procedure :: evaluate => run_test_mhd_newcomb_polar_map_binding
    end type

contains
    subroutine run_test_mhd_newcomb()
    !! V4(ii): ideal response with pressure and current against Newcomb's
    !! equation in the screw-pinch limit (independent ODE, fortnum Cash-Karp).
    !!
    !! A large-aspect-ratio torus (circular surfaces about (R0, 0), minor
    !! radius 1) carries the screw pinch q(r) = 1.1 + 0.3 r^2 (B_theta =
    !! r B0/(R0 q)), mu0 p = beta (B0/R0)^2 (1 - r^2)^2 and B_z from radial force
    !! balance, at the tokamak ordering (B_theta, p^(1/2) ~ 1/R0).  It is driven
    !! by the field of u = (x -+ i z)^3 (one poloidal harmonic, |m| = 3, n = 2,
    !! no rational surface).  Inside the plasma the perturbed equilibrium
    !! satisfies Newcomb's equation (f xi')' - g xi = 0 with
    !!     f = r F^2/k0^2,  g = (2k^2/k0^2) mu0 p' + (k0^2 r^2 - 1)/(k0^2 r^2) r F^2
    !!         + 2k^2/(r k0^4) (k^2 B_z^2 - m^2 B_theta^2/r^2),
    !! F = k B_z + m B_theta/r, k0^2 = k^2 + m^2/r^2, k = n/R0 (Freidberg, Ideal
    !! MHD, the screw-pinch delta-W; checked separately against the vacuum
    !! I_m and the Lundquist force-free perturbation).  Regularity at the axis
    !! leaves one solution, so the radial shape of the resonant-harmonic flux
    !! Phi_m(r) ~ r F xi is fixed independently of the edge.  The torus shape
    !! Phi_m(r)/Phi_m(r_ref) must approach it as a/R0 -> 0.
    !! At the tokamak ordering the cylindrical pressure term is O(eps^2) and
    !! toroidal curvature adds O(eps beta_p): beta = 0 checks the current
    !! terms (observed second order, 6e-3 at R0/a = 40); beta = 0.5
    !! (beta_p about 2) checks that the pressure path converges to the same
    !! limit (first order).  A pressure-sensitive Newcomb test needs the
    !! periodic-cylinder chart, which the per-n FEEC complex does not have.
    !$ use omp_lib, only: omp_get_max_threads, omp_set_num_threads
    implicit none

real(dp) , target :: beta
    real(dp) , target :: r0, rprof(nprof), bzprof(nprof)
    integer :: mk
    type(psi_grid_t) , target :: g
    logical :: all_passed
        type(run_test_mhd_newcomb_newcomb_rhs_context_t), target :: newcomb_rhs_data
        type(run_test_mhd_newcomb_newcomb_rhs_t), target :: newcomb_rhs_operator
        type(run_test_mhd_newcomb_polar_map_context_t), target :: polar_map_data
        type(run_test_mhd_newcomb_polar_map_t), target :: polar_map_operator

        beta = 2.0_dp
    all_passed = .true.
    !$ call omp_set_num_threads(min(4, omp_get_max_threads()))
    beta = 0.0_dp
    call run(1.8_dp, 1.0e-2_dp)
    beta = 0.5_dp
    call run(0.9_dp, 3.0e-2_dp)
    if (.not. all_passed) error stop 1
    print '(a)', "test_mhd_newcomb: all checks passed"

contains
    subroutine bind_run_test_mhd_newcomb_polar_map_context()
        polar_map_data%g => g
        polar_map_operator%data => polar_map_data
    end subroutine bind_run_test_mhd_newcomb_polar_map_context

    subroutine bind_run_test_mhd_newcomb_newcomb_rhs_context()
        newcomb_rhs_data%beta => beta
        newcomb_rhs_data%bzprof => bzprof
        newcomb_rhs_data%r0 => r0
        newcomb_rhs_data%rprof => rprof
        newcomb_rhs_operator%data => newcomb_rhs_data
    end subroutine bind_run_test_mhd_newcomb_newcomb_rhs_context


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

    real(dp) function mu0p(r)
        real(dp), intent(in) :: r
        mu0p = beta*(b0/r0)**2*max(0.0_dp, 1.0_dp - r*r)**2
    end function mu0p


    pure real(dp) function psi_of_r(r)
        !! B_theta = psi'/R0
        real(dp), intent(in) :: r
        psi_of_r = b0/0.6_dp*log((1.1_dp + 0.3_dp*r*r)/1.1_dp)
    end function psi_of_r

    subroutine make_bz()
        !! d/dr (mu0 p + B_z^2/2) = -(B_theta/r) d(r B_theta)/dr, B_z(1.3) = B0;
        !! trapezoidal integration of B_z^2 on a fine grid.
        real(dp) :: h, s1, s2, x
        integer :: i
        h = 1.3_dp/(nprof - 1)
        do i = 1, nprof
            rprof(i) = (i - 1)*h
        end do
        bzprof(nprof) = b0**2/2 + mu0p(1.3_dp)
        do i = nprof - 1, 1, -1
            s1 = src(rprof(i))
            s2 = src(rprof(i + 1))
            bzprof(i) = bzprof(i + 1) + 0.5_dp*h*(s1 + s2)
        end do
        do i = 1, nprof
            x = rprof(i)
            bzprof(i) = sqrt(2.0_dp*(bzprof(i) - mu0p(x)))
        end do
    end subroutine make_bz

    real(dp) function src(r)
        real(dp), intent(in) :: r
        real(dp) :: dd
        dd = 1.0e-6_dp
        if (r < dd) then
            src = 0.0_dp
        else
        call bind_run_test_mhd_newcomb_newcomb_rhs_context()
            src = run_test_mhd_newcomb_btheta(newcomb_rhs_data, r)/r*((r + &
                dd)*run_test_mhd_newcomb_btheta(newcomb_rhs_data, r + dd) - (r - &
                dd)*run_test_mhd_newcomb_btheta(newcomb_rhs_data, r - dd)) &
                /(2*dd)
        end if
    end function src


    subroutine make_equilibrium()
        integer, parameter :: ng = 201, nf = 129
        real(dp) :: psirz(ng, ng), fpol(nf), r, z, rr, x, pa
        integer :: i, j, st
        do j = 1, ng
            do i = 1, ng
                r = r0 - 1.4_dp + 2.8_dp*(i - 1)/(ng - 1)
                z = -1.4_dp + 2.8_dp*(j - 1)/(ng - 1)
                psirz(i, j) = psi_of_r(sqrt((r - r0)**2 + z**2))
            end do
        end do
        pa = psi_of_r(1.0_dp)
        do i = 1, nf
            x = (i - 1.0_dp)/(nf - 1)
            rr = r_of_psin(x)
        call bind_run_test_mhd_newcomb_newcomb_rhs_context()
            fpol(i) = r0*run_test_mhd_newcomb_bz(newcomb_rhs_data, rr)
        end do
        call psi_grid_init(g, psirz, r0 - 1.4_dp, 2.8_dp, 0.0_dp, 2.8_dp, &
            r0 + 0.01_dp, 0.01_dp, pa, 1.0e6_dp, fpol, st)
        if (st /= 0) error stop "equilibrium"
        allocate (g%pcoef(0:nf + 1))
        ! cubic B-spline coefficients reproducing the samples to O(h^2): the
        ! profile used by the code is this spline; the ODE uses the same
        ! values up to that error
        do i = 0, nf + 1
            x = min(1.0_dp, max(0.0_dp, (i - 1.0_dp)/(nf - 1)))
            g%pcoef(i) = mu0p(r_of_psin(x))/mu0
        end do
    end subroutine make_equilibrium

    real(dp) function r_of_psin(x)
        real(dp), intent(in) :: x
        real(dp) :: t
        ! psi = (B0/0.6) ln(1 + 0.3 r^2/1.1) inverted
        t = x*psi_of_r(1.0_dp)*0.6_dp/b0
        r_of_psin = sqrt(max(0.0_dp, (exp(t) - 1.0_dp)*1.1_dp/0.3_dp))
    end function r_of_psin


    subroutine drive(x, b)
        !! grad_n of u = (x - i s z)^3, s = sign(mk) (harmonic e^{-i mk theta}).
        real(dp), intent(in) :: x(2)
        complex(dp), intent(out) :: b(3)
        complex(dp) :: w, u
        real(dp) :: s
        s = sign(1.0_dp, real(mk, dp))
        w = cmplx(x(1) - r0, -s*x(2), dp)
        u = w**mabs
        b(1) = mabs*w**(mabs - 1)
        b(3) = cmplx(0.0_dp, -s, dp)*mabs*w**(mabs - 1)
        b(2) = cmplx(0.0_dp, nmode, dp)*u/x(1)
        b = 1.0e-3_dp*b
    end subroutine drive

    ! -------------------------------------------------------------- Newcomb

    subroutine newcomb_shape(rs, shape)
        !! r F xi at the radii rs (increasing), regular solution.
        real(dp), intent(in) :: rs(:)
        real(dp), intent(out) :: shape(:)
        real(dp) :: y(2), y5(2), yerr(2), k1(2), k2(2), k3(2), k4(2), k5(2)
        real(dp) :: k6(2), yt(2), t, h, f, gg, ff, bt
        integer :: i, nfev, j
        t = 1.0e-3_dp
        y(1) = t**(mabs - 1)
        call bind_run_test_mhd_newcomb_newcomb_rhs_context()
        call run_test_mhd_newcomb_newcomb_coeff(newcomb_rhs_data, t, f, gg)
        y(2) = f*(mabs - 1)*t**(mabs - 2)
        nfev = 0
        j = 1
        h = 1.0e-4_dp
        do i = 1, 200000
            if (j > size(rs)) exit
            if (t + h >= rs(j)) then
        call bind_run_test_mhd_newcomb_newcomb_rhs_context()
                call cash_karp_step(newcomb_rhs_bridge, t, y, rs(j) - t, .false., k1, &
                    k2, k3, k4, k5, k6, yt, y5, yerr, nfev, ctx=newcomb_rhs_operator)
                t = rs(j)
                y = y5
        call bind_run_test_mhd_newcomb_newcomb_rhs_context()
                bt = run_test_mhd_newcomb_btheta(newcomb_rhs_data, t)
        call bind_run_test_mhd_newcomb_newcomb_rhs_context()
                ff = bt/t*(nmode*run_test_mhd_newcomb_qloc(newcomb_rhs_data, t, bt, &
                    run_test_mhd_newcomb_bz(newcomb_rhs_data, t)) - mabs)
                shape(j) = t*ff*y(1)
                j = j + 1
                cycle
            end if
        call bind_run_test_mhd_newcomb_newcomb_rhs_context()
            call cash_karp_step(newcomb_rhs_bridge, t, y, h, .false., k1, k2, k3, k4, &
                k5, k6, yt, y5, yerr, nfev, ctx=newcomb_rhs_operator)
            t = t + h
            y = y5
        end do
    end subroutine newcomb_shape

    ! ----------------------------------------------------------------- run

    subroutine run(order, tol)
        real(dp), intent(in) :: order, tol
        real(dp), parameter :: r0s(3) = [10.0_dp, 20.0_dp, 40.0_dp]
        real(dp) :: err(3), rate
        character(len=24) :: tag
        integer :: i
        write (tag, '(a,f4.2)') ", beta ", beta
        do i = 1, 3
            r0 = r0s(i)
            call make_bz()
            call make_equilibrium()
            call level(err(i))
        end do
        rate = log(err(1)/err(3))/log(4.0_dp)
        print '(a,f6.2)', "  shape error order in a/R0", rate
        call record(rate > order, "torus response approaches the Newcomb screw pinch" &
            //trim(tag))
        call record(err(3) < tol, "Newcomb shape at R0/a = 40"//trim(tag))
    end subroutine run

    subroutine level(err)
        real(dp), intent(out) :: err
        integer, parameter :: nchk = 6
        integer, parameter :: rings(nchk) = [8, 11, 14, 17, 20, 23]
        integer, parameter :: jref = 4
        type(flux_mesh_t) :: fm
        type(meridian_feec_t) :: fs
        type(meridian_point_basis_t) :: pb
        type(ampere_n_t) :: amp
        type(ideal_options_t) :: opt
        type(sample_ring_t), allocatable :: sr(:)
        type(ring_angles_t) :: rac
        type(layer_closure_t) :: cl
        type(surface_t) :: s
        complex(dp), allocatable :: a(:), b2(:), bc(:)
        complex(dp) :: lam(0), b2s(0, 0), phe(0), pht(0), phi(nchk)
        real(dp), allocatable :: ul(:), x(:, :), nds(:, :), th(:), tau(:), wq(:)
        integer, allocatable :: tt(:)
        real(dp) :: rs(nchk), shp(nchk), res
        complex(dp) :: sh(nchk), bb(3)
        integer :: j, k, st, info, its

        allocate (ul(nrad))
        do j = 1, nrad
            ul(j) = real(j, dp)/nrad
        end do
        call bind_run_test_mhd_newcomb_polar_map_context()
        call flux_mesh_build(ul, 8, 256, 2, polar_map_operator, fm, st)
        call meridian_feec_build(fm%mm, 2, fs, st)
        call ampere_n_build(fs, nmode, amp, st)
        call meridian_point_basis_init(fs, pb, st)
        call surface_build(g, 0.25_dp, 256, s, st)
        mk = nint(sign(real(mabs, dp), s%q))
        allocate (sr(nrad))
        do j = 2, nrad
            call sample_ring_build(g, opt, fm, fs, pb, 0.5_dp*(ul(j - 1) + ul(j)), &
                j, np, sr(j - 1), st)
            if (st /= 0) error stop "half ring"
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
            pht, 1.0e-10_dp, 300, info, its, res)
        call record(info == 0, "GMRES converged")
        call ampere_n_curl(amp, a, b2)
        do j = 1, nchk
            call surface_build(g, ul(rings(j))**2, 512, s, st)
            call ring_angles_set(s%q, s%theta_s(:, 1), s%nu, rac)
            phi(j) = ring_flux_harmonic(fs, pb, fm, rings(j), rac, nmode, b2, mk)
            call ring_points(fm, rings(j), x, nds, th, tt, tau, wq)
            allocate (bc(size(th)))
            do k = 1, size(th)
                call drive(x(:, k), bb)
                bc(k) = (bb(1)*nds(1, k) + bb(3)*nds(2, k))*x(1, k)
            end do
            phi(j) = phi(j) + ring_flux_from_values(rac, nmode, mk, th, bc, wq)
            deallocate (bc)
            rs(j) = r_of_psin(ul(rings(j))**2)
        end do
        call newcomb_shape(rs, shp)
        sh = phi/phi(jref)
        err = maxval(abs(sh - shp/shp(jref)))
        print '(a,f5.1,a,i3)', "  R0/a =", r0, ", GMRES iterations", its
        do j = 1, nchk
            print '(a,f6.3,a,2f10.5,a,f10.5)', "    r", rs(j), "  torus", sh(j), &
                "  Newcomb", shp(j)/shp(jref)
        end do
        print '(a,es10.2)', "    max shape error", err
    end subroutine level


    end subroutine run_test_mhd_newcomb

    pure real(dp) function run_test_mhd_newcomb_qcyl(r)
        real(dp), intent(in) :: r
        run_test_mhd_newcomb_qcyl = 1.1_dp + 0.3_dp*r*r
    end function run_test_mhd_newcomb_qcyl

    pure real(dp) function run_test_mhd_newcomb_btheta(newcomb_rhs_data, r)
        type(run_test_mhd_newcomb_newcomb_rhs_context_t), intent(in), target :: &
            newcomb_rhs_data
        real(dp), intent(in) :: r
        run_test_mhd_newcomb_btheta = &
            r*b0/(newcomb_rhs_data%r0*run_test_mhd_newcomb_qcyl(r))
    end function run_test_mhd_newcomb_btheta

    real(dp) function run_test_mhd_newcomb_dmu0p(newcomb_rhs_data, r)
        type(run_test_mhd_newcomb_newcomb_rhs_context_t), intent(in), target :: &
            newcomb_rhs_data
        real(dp), intent(in) :: r
        run_test_mhd_newcomb_dmu0p = &
            -4.0_dp*r*newcomb_rhs_data%beta*(b0/newcomb_rhs_data%r0)**2*max(0.0_dp, &
            1.0_dp - r*r)
    end function run_test_mhd_newcomb_dmu0p

    real(dp) function run_test_mhd_newcomb_bz(newcomb_rhs_data, r)
        type(run_test_mhd_newcomb_newcomb_rhs_context_t), intent(in), target :: &
            newcomb_rhs_data
        real(dp), intent(in) :: r
        real(dp) :: s
        integer :: i
        s = r/newcomb_rhs_data%rprof(2)
        i = min(nprof - 1, int(s) + 1)
        s = s - (i - 1)
        run_test_mhd_newcomb_bz = (1.0_dp - s)*newcomb_rhs_data%bzprof(i) + &
            s*newcomb_rhs_data%bzprof(i + 1)
    end function run_test_mhd_newcomb_bz

    subroutine run_test_mhd_newcomb_newcomb_coeff(newcomb_rhs_data, r, f, gg)
        type(run_test_mhd_newcomb_newcomb_rhs_context_t), intent(in), target :: &
            newcomb_rhs_data
        real(dp), intent(in) :: r
        real(dp), intent(out) :: f, gg
        real(dp) :: k, k02, ff, bt, bzz, m
        k = nmode/newcomb_rhs_data%r0
        m = real(mabs, dp)
        bt = run_test_mhd_newcomb_btheta(newcomb_rhs_data, r)
        bzz = run_test_mhd_newcomb_bz(newcomb_rhs_data, r)
        ! F = (B_theta/r)(n |q| - |m|): the helicity of the kin6d harmonic
        ff = bt/r*(nmode*run_test_mhd_newcomb_qloc(newcomb_rhs_data, r, bt, bzz) - m)
        k02 = k*k + m*m/(r*r)
        f = r*ff**2/k02
        gg = 2*k*k/k02*run_test_mhd_newcomb_dmu0p(newcomb_rhs_data, r) + (k02*r*r - &
            1)/(k02*r*r)*r*ff**2 &
            + 2*k*k/(r*k02**2)*(k*k*bzz**2 - m*m*bt**2/(r*r))
    end subroutine run_test_mhd_newcomb_newcomb_coeff

    pure real(dp) function run_test_mhd_newcomb_qloc(newcomb_rhs_data, r, bt, bzz)
        type(run_test_mhd_newcomb_newcomb_rhs_context_t), intent(in), target :: &
            newcomb_rhs_data
        real(dp), intent(in) :: r, bt, bzz
        run_test_mhd_newcomb_qloc = r*bzz/(newcomb_rhs_data%r0*bt)
    end function run_test_mhd_newcomb_qloc

    subroutine run_test_mhd_newcomb_newcomb_rhs(newcomb_rhs_data, t, y, dydt, ctx)
        type(run_test_mhd_newcomb_newcomb_rhs_context_t), intent(in), target :: &
            newcomb_rhs_data
        real(dp), intent(in) :: t
        real(dp), intent(in) :: y(:)
        real(dp), intent(out) :: dydt(:)
        class(*), intent(in), optional :: ctx
        real(dp) :: f, gg
        call run_test_mhd_newcomb_newcomb_coeff(newcomb_rhs_data, t, f, gg)
        dydt(1) = y(2)/f
        dydt(2) = gg*y(1)
    end subroutine run_test_mhd_newcomb_newcomb_rhs

    subroutine run_test_mhd_newcomb_newcomb_rhs_binding(self, t, y, dydt, ctx)
        class(run_test_mhd_newcomb_newcomb_rhs_t), intent(in) :: self
        real(dp), intent(in) :: t
        real(dp), intent(in) :: y(:)
        real(dp), intent(out) :: dydt(:)
        class(*), intent(in), optional :: ctx
        call run_test_mhd_newcomb_newcomb_rhs(self%data, t, y, dydt, ctx)
    end subroutine
    subroutine newcomb_rhs_bridge(t, y, dydt, ctx)
        real(dp), intent(in) :: t, y(:)
        real(dp), intent(out) :: dydt(:)
        class(*), optional, intent(in) :: ctx
        if (.not. present(ctx)) error stop 'missing callback object'
        select type (ctx)
        class is (callback_function_t)
            call ctx%evaluate(t, y, dydt)
        class default
            error stop 'invalid callback object'
        end select
    end subroutine

    subroutine run_test_mhd_newcomb_polar_map(polar_map_data, u, theta, x)
        type(run_test_mhd_newcomb_polar_map_context_t), intent(in), target :: &
            polar_map_data
        real(dp), intent(in) :: u, theta
        real(dp), intent(out) :: x(2)
        real(dp) :: rho
        logical :: ok
        if (u <= 0.0_dp) then
            x = [polar_map_data%g%r_axis, polar_map_data%g%z_axis]
            return
        end if
        call flux_ray_point(polar_map_data%g, u*u, theta, rho, ok)
        if (.not. ok) error stop "ray"
        x = [polar_map_data%g%r_axis + rho*cos(theta), polar_map_data%g%z_axis + &
            rho*sin(theta)]
    end subroutine run_test_mhd_newcomb_polar_map

    subroutine run_test_mhd_newcomb_polar_map_binding(self, u, theta, x)
        class(run_test_mhd_newcomb_polar_map_t), intent(in) :: self
        real(dp), intent(in) :: u, theta
        real(dp), intent(out) :: x(2)
        call run_test_mhd_newcomb_polar_map(self%data, u, theta, x)
    end subroutine
end module test_mhd_newcomb_callbacks

program test_mhd_newcomb
    use test_mhd_newcomb_callbacks, only: run_test_mhd_newcomb
    implicit none
    call run_test_mhd_newcomb()
end program test_mhd_newcomb
