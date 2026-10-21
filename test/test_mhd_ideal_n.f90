program test_mhd_ideal_n
    !! L4 stage 3: volume current of the ideal perturbed equilibrium.
    !!
    !! 1. Magnetic differential equation on a sample surface: for f given,
    !!    r = B0 . grad f is formed independently (finite differences of f
    !!    along the field line in (theta_g, phi)), and mde_solve(r) returns f
    !!    (spectrally accurate).
    !! 2. Force-free identity (independent oracle).  With p0' = 0 and the
    !!    equilibrium current replaced by J0 = sigma B0/mu0 (sigma constant),
    !!    the model gives mu0 J = sigma B for every field B (J_perp =
    !!    (sigma/mu0)(B - g B0), B0 . grad h = (sigma/mu0) B0 . grad g).  The
    !!    response then solves the linear problem curl curl A = sigma (B_ext +
    !!    curl A) in the plasma (outside the axis fan), which is assembled and
    !!    solved directly with the FEEC operators, without magnetic
    !!    differential equations, half rings or angle interpolation.  The
    !!    MEPHIT-type response must converge to it in the ring spacing (the
    !!    strip current is that of its mid surface; observed second order).
    !$ use omp_lib, only: omp_get_max_threads, omp_set_num_threads
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortfem_meridian_feec_n, only: meridian_feec_t, meridian_feec_build, &
        meridian_point_basis_t, meridian_point_basis_init, meridian_qpoint
    use fortnum_krylov, only: complex_gmres_operator
    use kin6d_vacuum_surfaces, only: psi_grid_t, psi_grid_init, flux_ray_point
    use kin6d_mhd_mesh, only: flux_mesh_t, flux_mesh_build
    use kin6d_mhd_ampere_n, only: ampere_n_t, ampere_n_build, ampere_n_solve, &
        ampere_n_curl, ampere_n_error, ampere_n_error_t, mu0
    use kin6d_mhd_response_n, only: layer_closure_t
    use kin6d_mhd_ideal_n, only: sample_ring_t, sample_ring_build, &
        sample_field, mde_solve, ideal_options_t, ideal_response, &
        ideal_response_error_t, ideal_current, ideal_load
    use fortsparse, only: csc_matvec
    implicit none

    real(dp), parameter :: pi = 3.14159265358979323846_dp
    real(dp), parameter :: r0 = 3.0_dp, cpsi = 0.6_dp, fval = 4.5_dp
    real(dp), parameter :: sigma = 0.5_dp
    integer, parameter :: nmode = 2
    type(psi_grid_t) :: g
    logical :: all_passed
    character(len=16) :: mode
    type(flux_mesh_t), pointer :: pfm => null()
    type(meridian_feec_t), pointer :: pfs => null()
    type(ampere_n_t), pointer :: pamp => null()

    all_passed = .true.
    !$ call omp_set_num_threads(min(4, omp_get_max_threads()))
    call make_equilibrium()
    call get_command_argument(1, mode)
    if (trim(mode) == 'error') then
        call error_check()
        if (.not. all_passed) error stop 1
        stop
    end if
    call mde_check()
    call force_free_check()
    call divergence_check()
    call load_thread_check()
    if (.not. all_passed) error stop 1
    print '(a)', "test_mhd_ideal_n: all checks passed"

contains

    subroutine error_check()
        !! Final state versus a tightly solved discrete reference. This tests
        !! an algebraic indicator, not containment of the continuous solution.
        type(flux_mesh_t) :: fm
        type(meridian_feec_t) :: fs
        type(meridian_point_basis_t) :: pb
        type(ampere_n_t) :: amp
        type(ideal_options_t) :: opt
        type(sample_ring_t), allocatable :: rings(:)
        type(layer_closure_t) :: cl
        type(ideal_response_error_t) :: estimate
        type(ampere_n_error_t) :: independent
        complex(dp), allocatable :: a(:), reference(:), delta(:), b2(:), f(:)
        complex(dp), allocatable :: bf(:, :, :), jc(:, :, :), kap(:)
        complex(dp) :: lam(0), b2s(0, 0), phe(0), pht(0)
        real(dp), allocatable :: ul(:)
        real(dp) :: res, truth, efficiency, elapsed, tsolve, tol
        integer :: lev, nrad, j, k, st, info, its, c0, c1, cr, trial
        integer, parameter :: np = 128

        call error_leaf_check()
        do lev = 1, 3
            nrad = 4*2**(lev - 1)
            allocate (ul(nrad), rings(nrad))
            do j = 1, nrad
                ul(j) = real(j, dp)/nrad
            end do
            call flux_mesh_build(ul, 8, 64, 2, polar_map, fm, st)
            if (st /= 0) error stop 'error mesh'
            call meridian_feec_build(fm%mm, 2, fs, st)
            if (st /= 0) error stop 'error FEEC'
            call ampere_n_build(fs, nmode, amp, st)
            if (st /= 0) error stop 'error Ampere build'
            call meridian_point_basis_init(fs, pb, st)
            if (st /= 0) error stop 'error point basis'
            opt%force_free = .true.
            opt%sigma = sigma
            opt%boundary_sheet = .false.
            opt%strip_shape = .false.
            do j = 2, nrad
                call sample_ring_build(g, opt, fm, fs, pb, &
                    0.5_dp*(ul(j - 1) + ul(j)), j, np, rings(j - 1), st)
                if (st /= 0) error stop 'error half ring'
            end do
            call sample_ring_build(g, opt, fm, fs, pb, 1.0_dp, nrad, np, &
                rings(nrad), st)
            if (st /= 0) error stop 'error boundary ring'
            do j = 1, nrad
                do k = 1, np
                    call bext_fun(rings(j)%x(:, k), rings(j)%bc(:, k))
                end do
            end do
            allocate (a(amp%nn), reference(amp%nn), f(amp%nn), &
                bf(3, np, nrad), jc(3, np, nrad), kap(np))
            reference = (0.0_dp, 0.0_dp)
            call ideal_response(fs, pb, fm, g, amp, opt, nmode, rings, cl, &
                b2s, phe, reference, lam, pht, 1.0e-12_dp, 60, info, its, res)
            call record(info == 0, 'tight discrete reference converged')
            tol = 10.0_dp**(-lev - 1)
            do trial = 1, 3
                a = (0.0_dp, 0.0_dp)
                call system_clock(c0, cr)
                call ideal_response(fs, pb, fm, g, amp, opt, nmode, rings, cl, &
                    b2s, phe, a, lam, pht, tol, 60, info, its, res, error=estimate)
                call system_clock(c1)
                elapsed = real(c1 - c0, dp)/cr
                tsolve = elapsed - estimate%seconds
                delta = a - reference
                truth = sqrt(real(dot_product(delta, &
                    csc_matvec(amp%matrix, delta)), dp))
                efficiency = estimate%algebraic_estimate/truth
                print '(a,4i6,7es13.5)', 'ideal error: ', nrad, amp%nn, its, trial, &
                    truth, estimate%defect_energy, estimate%contraction_probe, &
                    efficiency, tsolve, estimate%seconds, &
                    estimate%seconds/max(tsolve, 1.0e-12_dp)
                call record(info == 0 .and. estimate%status == 0, &
                    'final-state diagnostic succeeds')
                call record(efficiency >= 1.0_dp .and. efficiency < 2.0_dp, &
                    'practical estimate contains this discrete reference, index < 2')
            end do
            call record(estimate%continuous_claim == 'unavailable-ideal-stability', &
                'no continuous certificate from directional contraction')
            call record(estimate%sampling_claim == 'unavailable-current-sampling', &
                'current-sampling remainder is explicit')
            call ampere_n_curl(amp, a, b2)
            do j = 1, nrad
                call sample_field(fs, rings(j), b2, bf(:, :, j))
                bf(:, :, j) = bf(:, :, j) + rings(j)%bc
            end do
            call ideal_current(rings, nmode, g%sgn_dpsi, opt%qlow, bf, jc, kap, 0)
            call ideal_load(fs, pb, fm, g, jc(:, :, 1:nrad - 1), kap, 0, f)
            call ampere_n_error(amp, f, a, independent, st)
            call record(st == 0, 'independent final frozen-current audit succeeds')
            call record(abs(independent%discrete_energy_error &
                - estimate%defect_energy) < 1.0e-9_dp*truth, &
                'reported defect uses the returned state and its actual current')
            call record(abs(amp%error%discrete_energy_error &
                - estimate%defect_energy) < 1.0e-9_dp*truth, &
                'amp error restored after internal probe solves')
            ! Matched production-tolerance cost, separate from the deliberately
            ! loose candidates used to resolve algebraic truth above rounding.
            do trial = 1, 3
                a = (0.0_dp, 0.0_dp)
                call system_clock(c0, cr)
                call ideal_response(fs, pb, fm, g, amp, opt, nmode, rings, cl, &
                    b2s, phe, a, lam, pht, 1.0e-10_dp, 60, info, its, res, &
                    error=estimate)
                call system_clock(c1)
                elapsed = real(c1 - c0, dp)/cr
                tsolve = elapsed - estimate%seconds
                print '(a,4i6,3es13.5)', 'ideal cost: ', nrad, amp%nn, its, trial, &
                    tsolve, estimate%seconds, &
                    estimate%seconds/max(tsolve, 1.0e-12_dp)
                call record(info == 0 .and. estimate%status == 0, &
                    'production-tolerance final-state diagnostic succeeds')
            end do
            if (lev == 1) then
                do j = 1, nrad
                    rings(j)%j0 = 100.0_dp*rings(j)%j0
                end do
                a = (0.0_dp, 0.0_dp)
                call ideal_response(fs, pb, fm, g, amp, opt, nmode, rings, cl, &
                    b2s, phe, a, lam, pht, 1.0e-14_dp, 1, info, its, res, &
                    error=estimate)
                call record(info == 1 .and. its == 1 .and. estimate%status == 0, &
                    'nonconverged candidate still has a final-state defect')
                call record(estimate%contraction_probe >= 1.0_dp, &
                    'strong-current probe detects a noncontraction')
                call record(estimate%algebraic_claim == &
                    'unavailable-noncontractive-probe', &
                    'noncontractive probe never supplies a finite tail estimate')
            end if
            ! Exact-recovery case: no equilibrium current or pressure gradient.
            do j = 1, nrad
                rings(j)%j0 = 0.0_dp
                rings(j)%gsig = 0.0_dp
            end do
            a = (0.0_dp, 0.0_dp)
            call ideal_response(fs, pb, fm, g, amp, opt, nmode, rings, cl, b2s, &
                phe, a, lam, pht, 1.0e-12_dp, 60, info, its, res, error=estimate)
            call record(info == 0 .and. estimate%status == 0, &
                'zero-current diagnostic succeeds')
            call record(estimate%defect_energy < 1.0e-14_dp .and. &
                estimate%algebraic_estimate < 1.0e-14_dp, &
                'exact recovery of zero response')
            deallocate (ul, rings, a, reference, f, bf, jc, kap)
        end do
    end subroutine error_check

    subroutine error_leaf_check()
        use kin6d_ideal_response_error_leaves, only: ideal_ratio, ideal_tail, &
            ideal_ratio_jvp, ideal_ratio_vjp, ideal_tail_jvp, ideal_tail_vjp
        real(dp) :: q, tail, dq, dtail, gd, ga, bd, bq

        call ideal_ratio(2.0_dp, 0.5_dp, q)
        call ideal_tail(2.0_dp, q, tail)
        call record(abs(tail - 8.0_dp/3) < 1.0e-14_dp, &
            'geometric-series formula against rational arithmetic')
        call ideal_ratio_jvp(2.0_dp, 0.5_dp, 0.3_dp, -0.2_dp, dq)
        call ideal_ratio_vjp(2.0_dp, 0.5_dp, 1.0_dp, gd, ga)
        call record(abs(dq + 0.1375_dp) < 1.0e-14_dp .and. &
            abs(dq - (0.3_dp*gd - 0.2_dp*ga)) < 1.0e-14_dp, &
            'ratio JVP and VJP against independent rational derivative')
        call ideal_tail_jvp(2.0_dp, q, 0.3_dp, 0.2_dp, dtail)
        call ideal_tail_vjp(2.0_dp, q, 1.0_dp, bd, bq)
        call record(abs(dtail - 10.0_dp/9) < 1.0e-14_dp .and. &
            abs(dtail - (0.3_dp*bd + 0.2_dp*bq)) < 1.0e-14_dp, &
            'tail JVP and VJP against independent rational derivative')
    end subroutine error_leaf_check

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

    subroutine make_equilibrium()
        !! Circular surfaces psi = cpsi ((R - R0)^2 + Z^2), constant F, no
        !! pressure: q = F/(2 cpsi R0) (1 + O(eps^2)) = 1.25, n q = 2.5 (no
        !! rational surface for n = 2).
        integer, parameter :: ng = 161
        real(dp) :: psirz(ng, ng), fpol(65), r, z
        integer :: i, j, st
        do j = 1, ng
            do i = 1, ng
                r = r0 - 1.6_dp + 3.2_dp*(i - 1)/(ng - 1)
                z = -1.6_dp + 3.2_dp*(j - 1)/(ng - 1)
                psirz(i, j) = cpsi*((r - r0)**2 + z**2)
            end do
        end do
        fpol = fval
        call psi_grid_init(g, psirz, r0 - 1.6_dp, 3.2_dp, 0.0_dp, 3.2_dp, &
            r0 + 0.01_dp, 0.01_dp, cpsi, 1.0e6_dp, fpol, st)
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

    subroutine mde_check()
        !! f = exp(i k theta_g) cos(theta_g); B0 . grad f by differences.
        type(flux_mesh_t) :: fm
        type(meridian_feec_t) :: fs
        type(meridian_point_basis_t) :: pb
        type(sample_ring_t) :: sr
        type(ideal_options_t) :: opt
        integer, parameter :: np = 256
        complex(dp) :: f(np), r(np), fr(np)
        real(dp) :: th, h, err, ul(8), bt, dfdth
        integer :: k, st
        do k = 1, 8
            ul(k) = k/8.0_dp
        end do
        call flux_mesh_build(ul, 8, 64, 2, polar_map, fm, st)
        call meridian_feec_build(fm%mm, 2, fs, st)
        call meridian_point_basis_init(fs, pb, st)
        call sample_ring_build(g, opt, fm, fs, pb, 0.55_dp, 5, np, sr, st)
        if (st /= 0) error stop "sample ring"
        h = 1.0e-5_dp
        do k = 1, np
            th = 2*pi*(k - 1)/np
            f(k) = ftest(th)
            dfdth = real((ftest(th + h) - ftest(th - h))/(2*h), dp)
            ! B0 . grad f = (B_t/|x'|) df/dtheta + i n (F/R^2) f, with
            ! |x'|/B_t = fac and F/R^2 = B0_phi/R
            bt = 1.0_dp/sr%fac(k)
            r(k) = bt*(ftest(th + h) - ftest(th - h))/(2*h) &
                + cmplx(0.0_dp, nmode, dp)*sr%b0(2, k)/sr%x(1, k)*f(k)
        end do
        call mde_solve(sr, nmode, r, fr)
        err = maxval(abs(fr - f))/maxval(abs(f))
        print '(a,es10.2)', "  MDE inverse error", err
        call record(err < 1.0e-6_dp, "magnetic differential equation inverse")
    end subroutine mde_check

    complex(dp) function ftest(th)
        real(dp), intent(in) :: th
        ftest = exp(cmplx(0.0_dp, 3*th, dp))*cos(th) + 0.3_dp
    end function ftest

    subroutine divergence_check()
        use kin6d_mhd_ideal_n, only: ideal_current
        type(flux_mesh_t) :: fm
        type(meridian_feec_t) :: fs
        type(meridian_point_basis_t) :: pb
        type(ideal_options_t) :: opt
        type(sample_ring_t) :: rings(3)
        type(psi_grid_t) :: g0
        integer, parameter :: np = 256
        complex(dp) :: b(3, np, 3), jc(3, np, 3), kap(np), dj(np), djp(np)
        complex(dp) :: jperp(3, np, 3)
        real(dp) :: ul(16), du, uc, jac(np, 3), xx, gth(2)
        integer :: i, k, st, nf
        g0 = g
        ! profiles with p' /= 0 and F' /= 0 (any smooth profiles: the identity
        ! does not need an equilibrium)
        nf = size(g%fcoef) - 2
        do i = 0, nf + 1
            xx = real(i - 1, dp)/(nf - 1)
            g%fcoef(i) = fval*(1.0_dp - 0.08_dp*xx**2)
        end do
        allocate (g%pcoef(0:nf + 1))
        do i = 0, nf + 1
            xx = real(i - 1, dp)/(nf - 1)
            g%pcoef(i) = 2.0e5_dp*(1.0_dp - xx**2)**2
        end do
        do k = 1, 16
            ul(k) = k/16.0_dp
        end do
        call flux_mesh_build(ul, 8, 64, 2, polar_map, fm, st)
        call meridian_feec_build(fm%mm, 2, fs, st)
        call meridian_point_basis_init(fs, pb, st)
        uc = 0.6_dp
        du = 2.0e-3_dp
        do i = 1, 3
            call sample_ring_build(g, opt, fm, fs, pb, uc + (i - 2)*du, 10, np, &
                rings(i), st)
            if (st /= 0) error stop "divergence rings"
            do k = 1, np
                call btest_fun(rings(i)%x(:, k), b(:, k, i))
                gth = rings(i)%gth(:, k)
                jac(k, i) = rings(i)%x(1, k)/(rings(i)%gpsi(2, k)*gth(1) &
                    - rings(i)%gpsi(1, k)*gth(2))
            end do
        end do
        call ideal_current(rings, nmode, g%sgn_dpsi, 0.0_dp, b, jc, kap, 0)
        ! J_perp alone: remove the parallel part h B0 (h = J . B0/B0^2 minus
        ! the parallel part of J_perp, which vanishes)
        do i = 1, 3
            do k = 1, np
                jperp(:, k, i) = jc(:, k, i) - dot_product(rings(i)%b0(:, k), &
                    jc(:, k, i))/rings(i)%b0sq(k)*rings(i)%b0(:, k)
            end do
        end do
        call fd_div(rings, jac, jc, dj)
        call fd_div(rings, jac, jperp, djp)
        print '(a,es10.2)', "  max|div J|/max|div J_perp|", maxval(abs(dj))/maxval(abs(djp))
        call record(maxval(abs(dj)) < 1.0e-3_dp*maxval(abs(djp)), &
            "charge conservation with pressure and poloidal current")
        g = g0
    end subroutine divergence_check

    subroutine load_thread_check()
        !! ideal_load with the near-rational strip shape (psi and q at every
        !! quadrature point in the OpenMP region) is independent of the
        !! thread count: 1 against 4 threads to roundoff.
        use kin6d_mhd_ideal_n, only: ideal_load
        integer, parameter :: np = 128, nu = 16
        type(flux_mesh_t) :: fm
        type(meridian_feec_t) :: fs
        type(meridian_point_basis_t) :: pb
        type(ideal_options_t) :: opt
        type(sample_ring_t) :: rings(nu - 1)
        complex(dp) :: jc(3, np, nu - 1), kap(np)
        complex(dp), allocatable :: f1(:), f4(:)
        real(dp) :: ul(nu), th
        integer :: i, k, st, nth0
        do k = 1, nu
            ul(k) = k/real(nu, dp)
        end do
        call flux_mesh_build(ul, 8, 64, 2, polar_map, fm, st)
        call meridian_feec_build(fm%mm, 2, fs, st)
        call meridian_point_basis_init(fs, pb, st)
        do i = 2, nu
            call sample_ring_build(g, opt, fm, fs, pb, 0.5_dp*(ul(i - 1) + ul(i)), i, &
                np, rings(i - 1), st)
            if (st /= 0) error stop "load rings"
            do k = 1, np
                th = 2*pi*(k - 1)/np
                jc(:, k, i - 1) = [cmplx(cos(3*th), sin(2*th), dp), &
                    cmplx(i*0.1_dp, cos(th), dp), cmplx(sin(th), 0.3_dp, dp)]/ &
                    (1.0_dp + 0.2_dp*i)
            end do
        end do
        kap = (0.0_dp, 0.0_dp)
        allocate (f1(fs%sp(1)%ndof), f4(fs%sp(1)%ndof))
        nth0 = 1
        !$ nth0 = omp_get_max_threads()
        !$ call omp_set_num_threads(1)
        call ideal_load(fs, pb, fm, g, jc, kap, 0, f1, rings, nmode)
        !$ call omp_set_num_threads(4)
        call ideal_load(fs, pb, fm, g, jc, kap, 0, f4, rings, nmode)
        !$ call omp_set_num_threads(nth0)
        print '(a,es10.2)', "  ideal_load 1 vs 4 threads, max rel diff", &
            maxval(abs(f4 - f1))/maxval(abs(f1))
        call record(maxval(abs(f4 - f1)) <= 1.0e-13_dp*maxval(abs(f1)), &
            "ideal_load (strip shape) independent of the thread count")
    end subroutine load_thread_check

    subroutine fd_div(rings, jac, v, d)
        !! div v at the middle of three surfaces: (1/J)[d_psi(J v^psi) +
        !! d_theta(J v^theta)] + i n v_phi/R, J the (psi, theta_g, phi) Jacobian.
        type(sample_ring_t), intent(in) :: rings(3)
        real(dp), intent(in) :: jac(:, :)
        complex(dp), intent(in) :: v(:, :, :)
        complex(dp), intent(out) :: d(:)
        complex(dp) :: fpsi(size(d), 3), jth(size(d)), djth(size(d))
        integer :: np
        real(dp) :: dpsi
        integer :: ii, kk
        np = size(d)
        do ii = 1, 3
            do kk = 1, np
                fpsi(kk, ii) = jac(kk, ii)*(v(1, kk, ii)*rings(ii)%gpsi(1, kk) &
                    + v(3, kk, ii)*rings(ii)%gpsi(2, kk))
            end do
        end do
        dpsi = rings(3)%psi - rings(1)%psi
        do kk = 1, np
            jth(kk) = jac(kk, 2)*(v(1, kk, 2)*rings(2)%gth(1, kk) &
                + v(3, kk, 2)*rings(2)%gth(2, kk))
        end do
        call spectral_dtheta(jth, djth)
        do kk = 1, np
            d(kk) = ((fpsi(kk, 3) - fpsi(kk, 1))/dpsi + djth(kk))/jac(kk, 2) &
                + cmplx(0.0_dp, nmode, dp)*v(2, kk, 2)/rings(2)%x(1, kk)
        end do
    end subroutine fd_div
    subroutine spectral_dtheta(f, df)
        complex(dp), intent(in) :: f(:)
        complex(dp), intent(out) :: df(:)
        complex(dp) :: c
        integer :: kk, jj, m, np
        np = size(f)
        df = (0.0_dp, 0.0_dp)
        do m = -np/2 + 1, np/2 - 1
            c = (0.0_dp, 0.0_dp)
            do jj = 1, np
                c = c + f(jj)*exp(cmplx(0.0_dp, -m*2*pi*(jj - 1)/np, dp))
            end do
            c = c/np
            do kk = 1, np
                df(kk) = df(kk) + cmplx(0.0_dp, m, dp)*c &
                    *exp(cmplx(0.0_dp, m*2*pi*(kk - 1)/np, dp))
            end do
        end do
    end subroutine spectral_dtheta

    subroutine btest_fun(x, b)
        !! Vacuum per-n field grad_n (R^n (1 + Z/2)) (curl and divergence free).
        real(dp), intent(in) :: x(2)
        complex(dp), intent(out) :: b(3)
        real(dp) :: u
        u = x(1)**nmode*(1.0_dp + 0.5_dp*x(2))
        b(1) = nmode*x(1)**(nmode - 1)*(1.0_dp + 0.5_dp*x(2))
        b(2) = cmplx(0.0_dp, nmode, dp)*u/x(1)
        b(3) = 0.5_dp*x(1)**nmode
        b = b*1.0e-3_dp
    end subroutine btest_fun

    subroutine bext_fun(x, b)
        !! Exact per-n vacuum field grad_n (R^n).
        real(dp), intent(in) :: x(2)
        complex(dp), intent(out) :: b(3)
        b = [cmplx(nmode*x(1)**(nmode - 1), 0.0_dp, dp), &
            cmplx(0.0_dp, nmode*x(1)**(nmode - 1), dp), (0.0_dp, 0.0_dp)]
        b = b*1.0e-3_dp
    end subroutine bext_fun

    subroutine force_free_check()
        integer, parameter :: nlev = 3
        real(dp) :: diff(nlev), rate
        integer :: lev
        do lev = 1, nlev
            call force_free_level(8*2**lev, diff(lev))
        end do
        rate = log(diff(1)/diff(nlev))/(log(2.0_dp)*(nlev - 1))
        print '(a,f6.2)', "  force-free identity: order in ring spacing", rate
        call record(rate > 1.5_dp, "force-free response converges to the direct solve")
        call record(diff(nlev) < 2.0e-4_dp, "force-free response within 2e-4 on the finest level")
    end subroutine force_free_check

    subroutine force_free_level(nrad, diff)
        integer, intent(in) :: nrad
        real(dp), intent(out) :: diff
        type(flux_mesh_t), target :: fm
        type(meridian_feec_t), target :: fs
        type(meridian_point_basis_t) :: pb
        type(ampere_n_t), target :: amp
        type(ideal_options_t) :: opt
        type(sample_ring_t), allocatable :: rings(:)
        type(layer_closure_t) :: cl
        complex(dp), allocatable :: a(:), ad(:), b2(:), b2d(:)
        complex(dp), allocatable :: bs(:, :), bsd(:, :)
        complex(dp) :: lam(0), b2s(0, 0), phe(0), pht(0)
        real(dp), allocatable :: ul(:)
        integer :: j, st, info, its, np, k
        real(dp) :: res

        allocate (ul(nrad))
        do j = 1, nrad
            ul(j) = real(j, dp)/nrad
        end do
        call flux_mesh_build(ul, 8, 256, 2, polar_map, fm, st)
        call meridian_feec_build(fm%mm, 2, fs, st)
        call ampere_n_build(fs, nmode, amp, st)
        call meridian_point_basis_init(fs, pb, st)
        opt%force_free = .true.
        opt%sigma = sigma
        opt%boundary_sheet = .false.
        np = 256
        allocate (rings(nrad))
        do j = 2, nrad
            call sample_ring_build(g, opt, fm, fs, pb, 0.5_dp*(ul(j - 1) + ul(j)), &
                j, np, rings(j - 1), st)
            if (st /= 0) error stop "half ring"
        end do
        call sample_ring_build(g, opt, fm, fs, pb, 1.0_dp, nrad, np, rings(nrad), st)
        if (st /= 0) error stop "boundary ring"
        do j = 1, nrad
            do k = 1, np
                call bext_fun(rings(j)%x(:, k), rings(j)%bc(:, k))
            end do
        end do
        allocate (a(amp%nn))
        a = (0.0_dp, 0.0_dp)
        call ideal_response(fs, pb, fm, g, amp, opt, nmode, rings, cl, b2s, phe, a, &
            lam, pht, 1.0e-10_dp, 400, info, its, res)
        print '(a,i4,a,i4,a,es10.2)', "  nrad", nrad, ": GMRES iterations", its, &
            " residual", res
        call record(info == 0, "GMRES converged")
        ! direct solve: (I - sigma K^{-1} G) a = sigma K^{-1} b_ext
        pfm => fm
        pfs => fs
        pamp => amp
        call direct(fm, fs, amp, ad)
        call ampere_n_curl(amp, a, b2)
        call ampere_n_curl(amp, ad, b2d)
        allocate (bs(3, np), bsd(3, np))
        diff = 0.0_dp
        res = 0.0_dp
        do j = nrad/4, nrad - 1, max(1, nrad/8)
            call sample_field(fs, rings(j), b2, bs)
            call sample_field(fs, rings(j), b2d, bsd)
            diff = diff + sum(abs(bs - bsd)**2)
            res = res + sum(abs(bsd)**2)
        end do
        diff = sqrt(diff/res)
        print '(a,i4,a,es10.3,a,es10.3)', "  nrad", nrad, &
            ": |B_mde - B_direct|/|B_direct|", diff, "  |B_resp|", sqrt(res)
    end subroutine force_free_level

    subroutine direct(fm, fs, amp, a)
        !! Assemble G a = (A', curl A) and b = (A', B_ext) on the strips outside
        !! the fan and solve a - sigma K^{-1} (G a + b) = 0 by GMRES.
        type(flux_mesh_t), intent(in) :: fm
        type(meridian_feec_t), intent(in) :: fs
        type(ampere_n_t), intent(inout) :: amp
        complex(dp), allocatable, intent(out) :: a(:)
        complex(dp), allocatable :: rhs(:), bext(:), y(:)
        real(dp), allocatable :: phi1(:, :)
        real(dp) :: x(2), wq
        complex(dp) :: be(3)
        integer :: t, q, i, st, info, its
        real(dp) :: res

        allocate (bext(amp%nn), rhs(amp%nn), a(amp%nn), y(amp%nn))
        allocate (phi1(3, fs%sp(1)%nloc))
        bext = (0.0_dp, 0.0_dp)
        do t = fm%tlast(1) + 1, fs%mm%topology%n_triangles
            do q = 1, fs%nq
                call meridian_qpoint(fs, 1, t, q, x, wq, phi1, st)
                call bext_fun(x, be)
                do i = 1, fs%nln
                    bext(fs%sp(1)%dof(i, t)) = bext(fs%sp(1)%dof(i, t)) + wq* &
                        real(fs%sp(1)%sgn(i, t), dp)*sum(phi1(:, i)*be)
                end do
            end do
        end do
        call ampere_n_solve(amp, mu0*sigma/mu0*bext, rhs, st)
        a = (0.0_dp, 0.0_dp)
        call complex_gmres_operator(op, rhs, a, 1.0e-11_dp, 400, 60, info, its, res)
        if (info /= 0) error stop "direct GMRES"
    end subroutine direct

    subroutine op(xin, yout)
        complex(dp), intent(in) :: xin(:)
        complex(dp), intent(out) :: yout(:)
        complex(dp), allocatable :: b2(:), g(:), kg(:)
        complex(dp) :: bq(3)
        integer :: tt, qq, ii, ll, s2
        real(dp) :: xx(2), w2
        real(dp), allocatable :: phi1(:, :), phi2(:, :)
        call ampere_n_curl(pamp, xin, b2)
        allocate (g(size(xin)), kg(size(xin)), phi2(3, pfs%sp(2)%nloc), &
            phi1(3, pfs%sp(1)%nloc))
        g = (0.0_dp, 0.0_dp)
        do tt = pfm%tlast(1) + 1, pfs%mm%topology%n_triangles
            do qq = 1, pfs%nq
                call meridian_qpoint(pfs, 2, tt, qq, xx, w2, phi2, s2)
                bq = (0.0_dp, 0.0_dp)
                do ll = 1, pfs%sp(2)%nloc
                    bq = bq + real(pfs%sp(2)%sgn(ll, tt), dp)* &
                        b2(pfs%sp(2)%dof(ll, tt))*phi2(:, ll)
                end do
                call meridian_qpoint(pfs, 1, tt, qq, xx, w2, phi1, s2)
                do ii = 1, pfs%nln
                    g(pfs%sp(1)%dof(ii, tt)) = g(pfs%sp(1)%dof(ii, tt)) + w2* &
                        real(pfs%sp(1)%sgn(ii, tt), dp)*sum(phi1(:, ii)*bq)
                end do
            end do
        end do
        call ampere_n_solve(pamp, sigma*g, kg, s2)
        yout = xin - kg
    end subroutine op

end program test_mhd_ideal_n
