!> Spectral elements on the orbit-class mesh against the Legendre-Fourier
!! reference milca_dk_mono and against limits derived independently by
!! averaging the same equation (as in test_dk_mono): collisionless, banana
!! and Pfirsch-Schlueter, at collisionalities below the reach of the
!! reference. Field b = 1/(1 + eps cos theta), kappa = 1.
!! The modes scan, adapt and mono are timing drivers, not registered tests:
!!   test_dk_hp scan nu p n_col n_bulk delta   uniform mesh
!!   test_dk_hp adapt nu p tol                 adaptive D31 to relative tol
!!   test_dk_hp mono nu n_theta n_leg          Legendre-Fourier reference
program test_dk_hp
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use milca_dk_hp, only: hp_field_t, hp_field, hp_mesh_t, hp_default_mesh, &
        hp_coefficients, HP_OK
    use milca_dk_mono, only: mono_coefficients_t, mono_coefficients, MONO_OK
    implicit none
    real(dp), parameter :: pi = 4*atan(1.0_dp)
    ! nu_D/v (1/m) at x = 1 of the energy-coupled cases.
    real(dp), parameter :: energy_kappa = 3e-2_dp
    character(32) :: check
    logical :: ok

    call get_command_argument(1, check)
    select case (trim(check))
    case ('reference'); ok = reference()
    case ('onsager'); ok = onsager()
    case ('collisionless'); ok = collisionless()
    case ('banana'); ok = banana()
    case ('pfirsch_schlueter'); ok = pfirsch_schlueter()
    case ('estimate'); ok = estimate()
    case ('ntv_reference'); ok = ntv_reference()
    case ('ntv_superbanana_plateau'); ok = ntv_superbanana_plateau()
    case ('ntv_one_over_nu'); ok = ntv_one_over_nu()
    case ('ntv_rotation'); ok = ntv_rotation()
    case ('terms_recombine'); ok = terms_recombine()
    case ('line_moments'); ok = line_moments()
    case ('transfer'); ok = transfer()
    case ('energy_lorentz'); ok = energy_lorentz()
    case ('energy_test_particle'); ok = energy_test_particle()
    case ('energy_two_mesh'); ok = energy_two_mesh()
    case ('energy_density'); ok = energy_density()
    case ('energy_iterations'); ok = energy_iterations()
    case ('energy_full'); ok = energy_full()
    case ('collision_identities'); ok = collision_identities()
    case ('probe'); ok = probe()
    case ('ntv_scan'); ok = ntv_scan()
    case ('ntv_legendre'); ok = ntv_legendre()
    case ('scan'); ok = scan()
    case ('mono'); ok = mono_scan()
    case ('adapt'); ok = adapt_scan()
    case default
        print '(a)', 'unknown check '//trim(check)
        ok = .false.
    end select
    if (.not. ok) error stop 1

contains

    pure function field(n, eps) result(b)
        integer, intent(in) :: n
        real(dp), intent(in) :: eps
        real(dp) :: b(n)
        integer :: j
        b = [(1/(1 + eps*cos(2*pi*(j - 1)/n)), j=1, n)]
    end function field

    function hp_solve(eps, nu, p, n_col, n_bulk, delta) result(c)
        real(dp), intent(in) :: eps, nu, delta
        integer, intent(in) :: p, n_col, n_bulk
        type(mono_coefficients_t) :: c
        type(hp_field_t) :: f
        type(hp_mesh_t) :: mesh
        integer :: status
        call hp_field(field(33, eps), 1.0_dp, f, status)
        if (status /= HP_OK) error stop 'hp_field failed'
        call hp_default_mesh(f, p, n_col, n_bulk, delta, 0.3_dp, mesh)
        call hp_coefficients(f, mesh, nu, c, status)
        if (status /= HP_OK) error stop 'hp_coefficients failed'
    end function hp_solve

    logical function close_to(name, got, want, tol) result(ok)
        character(*), intent(in) :: name
        real(dp), intent(in) :: got, want, tol
        ok = abs(got - want) <= tol*abs(want) .and. abs(want) < huge(want) &
            .and. abs(got) < huge(got)
        print '(a,es14.6,a,es14.6,a,es9.2,a,l2)', name//' ', got, ' oracle ', want, &
            ' rel ', abs(got - want)/abs(want), ' pass', ok
    end function close_to

    !> Same equation, independent discretisation, at moderate collisionality
    !> where the reference is converged to rounding (65 theta, 200 Legendre).
    logical function reference() result(ok)
        real(dp), parameter :: eps = 0.1_dp
        real(dp) :: nu
        type(mono_coefficients_t) :: c, m
        integer :: status, i
        ok = .true.
        do i = 1, 3
            nu = 10.0_dp**(-i)
            call mono_coefficients(field(65, eps), 1.0_dp, nu, 200, m, status)
            if (status /= MONO_OK) error stop 'mono failed'
            c = hp_solve(eps, nu, 6, 8, 4, 1e-2_dp)
            ok = close_to('D11', c%d11, m%d11, 1e-8_dp) .and. ok
            ok = close_to('D31', c%d31, m%d31, 1e-8_dp) .and. ok
            ok = close_to('D33', c%d33, m%d33, 1e-8_dp) .and. ok
        end do
    end function reference

    !> D31 = -D13 by the parity of drives and operator, and the density
    !> sources vanish; the streaming bracket is exactly skew, so both hold to
    !> rounding amplified by the conditioning, which grows like 1/nu.
    logical function onsager() result(ok)
        type(mono_coefficients_t) :: c
        real(dp) :: scale, nu
        integer :: i
        ok = .true.
        do i = -6, 2, 2
            nu = 10.0_dp**i
            c = hp_solve(0.1_dp, nu, 4, 8, 4, min(0.1_dp, sqrt(nu)))
            scale = sqrt(c%d11*c%d33)
            ok = ok .and. abs(c%d31 + c%d13) <= 1e-9_dp*scale .and. &
                maxval(abs(c%border)) <= 1e-10_dp*scale
            print '(a,es9.1,a,es9.2,a,es9.2)', 'nu ', nu, &
                ' |D31 + D13|/sqrt(D11 D33) ', abs(c%d31 + c%d13)/scale, &
                ' |mu|/sqrt(D11 D33) ', maxval(abs(c%border))/scale
        end do
    end function onsager

    pure real(dp) function average(x, b)
        real(dp), intent(in) :: x(:), b(:)
        average = sum(x/b**2)/sum(1/b**2)
    end function average

    ! int_0^{1/b_max} lam dlam/<sqrt(1 - lam b)>, with lam = lam_c (1 - t^2).
    real(dp) function passing_integral(eps) result(q)
        use fortnum_quadrature, only: gauss_legendre_ab
        real(dp), intent(in) :: eps
        integer, parameter :: m = 40001, nq = 64
        real(dp) :: b(m), t(nq), w(nq), lam_c, lam
        integer :: i
        b = field(m, eps)
        lam_c = 1/maxval(b)
        call gauss_legendre_ab(nq, 0.0_dp, 1.0_dp, t, w)
        q = 0
        do i = 1, nq
            lam = lam_c*(1 - t(i)**2)
            q = q + w(i)*2*lam_c*t(i)*lam/average(sqrt(max(1 - lam*b, 0.0_dp)), b)
        end do
    end function passing_integral

    ! Value at nu -> 0 from two collisionalities with X = X0 + c sqrt(nu).
    pure real(dp) function richardson(x1, nu1, x2, nu2)
        real(dp), intent(in) :: x1, nu1, x2, nu2
        richardson = (x2*sqrt(nu1) - x1*sqrt(nu2))/(sqrt(nu1) - sqrt(nu2))
    end function richardson

    !> 3 nu D33/<b^2> -> f_c and 3 D31 -> 1 - f_c. At nu = 1e-6 and 2.5e-7
    !> the next term of the expansion leaves about 1e-6 after Richardson; the
    !> Legendre-Fourier reference needs 1e-4 and 2.5e-5 for its 1e-4 check.
    logical function collisionless() result(ok)
        real(dp), parameter :: eps = 0.1_dp, nu1 = 1e-6_dp, nu2 = 2.5e-7_dp
        type(mono_coefficients_t) :: c1, c2
        real(dp) :: b(40001), b2, f_c
        b = field(size(b), eps)
        b2 = average(b**2, b)
        f_c = 0.75_dp*b2*passing_integral(eps)
        c1 = hp_solve(eps, nu1, 6, 16, 4, 1e-5_dp)
        c2 = hp_solve(eps, nu2, 6, 16, 4, 5e-6_dp)
        ok = close_to('3 nu D33/<b^2>', richardson(3*nu1*c1%d33/b2, nu1, &
            3*nu2*c2%d33/b2, nu2), f_c, 1e-5_dp)
        ok = close_to('3 D31', richardson(3*c1%d31, nu1, 3*c2%d31, nu2), &
            1 - f_c, 1e-5_dp) .and. ok
    end function collisionless

    logical function banana() result(ok)
        real(dp), parameter :: eps = 0.1_dp, nu1 = 1e-6_dp, nu2 = 2.5e-7_dp
        type(mono_coefficients_t) :: c1, c2
        real(dp) :: b(40001), want
        b = field(size(b), eps)
        want = average(1/b**2, b)/3 - passing_integral(eps)/4
        c1 = hp_solve(eps, nu1, 6, 16, 4, 1e-5_dp)
        c2 = hp_solve(eps, nu2, 6, 16, 4, 5e-6_dp)
        ok = close_to('D11/nu', richardson(c1%d11/nu1, nu1, c2%d11/nu2, nu2), &
            want, 1e-5_dp)
    end function banana

    ! D11/nu -> <s_0 A> with d_theta A = -(b^-2 - mean(b^-2))/kappa and the
    ! P_0 drive s_0 = -(2/3) kappa b_theta/b (derivation in test_dk_mono).
    ! The correction falls like 1/nu^2 (5e-6 at nu = 1e2); far above nu = 1e3
    ! rounding grows with the conditioning of the collision-dominated system.
    logical function pfirsch_schlueter() result(ok)
        real(dp), parameter :: eps = 0.1_dp, nu = 1e3_dp
        integer, parameter :: m = 200001
        type(mono_coefficients_t) :: c
        real(dp), allocatable :: b(:), theta(:), s0(:), integrand(:), big_a(:)
        real(dp) :: h
        integer :: j
        allocate (b(m), theta(m), s0(m), integrand(m), big_a(m))
        b = field(m, eps)
        h = 2*pi/m
        theta = [(h*(j - 1), j=1, m)]
        s0 = -(2.0_dp/3)*eps*sin(theta)*b
        integrand = -(1/b**2 - sum(1/b**2)/m)
        big_a(1) = 0
        do j = 2, m
            big_a(j) = big_a(j - 1) + h*(integrand(j - 1) + integrand(j))/2
        end do
        c = hp_solve(eps, nu, 6, 8, 2, 0.25_dp)
        ok = close_to('D11/nu', c%d11/nu, average(s0*big_a, b), 1e-6_dp)
    end function pfirsch_schlueter

    !> The dual-weighted residual estimate of D31 on a coarse mesh predicts
    !> its error against the reference (converged to 2e-10 at this nu) within
    !> a factor 0.5 to 2, and one adaptive refinement step reduces the error
    !> tenfold.
    logical function estimate() result(ok)
        use milca_dk_hp, only: hp_estimate, hp_refine
        real(dp), parameter :: eps = 0.1_dp, nu = 1e-4_dp
        type(hp_field_t) :: f
        type(hp_mesh_t) :: mesh, next
        type(mono_coefficients_t) :: fine
        real(dp) :: value, eta, value2, eta2, err, err2
        real(dp), allocatable :: eta_el(:)
        integer :: status
        call hp_field(field(33, eps), 1.0_dp, f, status)
        call hp_default_mesh(f, 4, 4, 2, 1e-2_dp, 0.3_dp, mesh)
        call hp_estimate(f, mesh, nu, 1, 2, value, eta, eta_el, status)
        if (status /= HP_OK) error stop 'hp_estimate failed'
        call mono_coefficients(field(129, eps), 1.0_dp, nu, 400, fine, status)
        err = fine%d31 - value
        call hp_refine(mesh, eta_el, 0.5_dp, next)
        call hp_estimate(f, next, nu, 1, 2, value2, eta2, eta_el, status)
        err2 = fine%d31 - value2
        print '(a,3es12.4)', 'error, estimate, effectivity ', err, eta, eta/err
        print '(a,2es12.4)', 'error after refinement, estimate ', err2, eta2
        ok = eta/err >= 0.5_dp .and. eta/err <= 2 .and. abs(err2) < abs(err)/10
    end function estimate

    !> One toroidal harmonic with streaming shift, precession and drive,
    !> against the Legendre-Fourier kernel milca_dk_ntv on the circular
    !> stage-1 surface b = 1 - eps cos theta.
    logical function ntv_reference() result(ok)
        use milca_dk_ntv, only: ntv_surface_t, ntv_kernel_t, ntv_kernel, NTV_OK, &
            NTV_WEIGHT_FULL
        use milca_dk_hp, only: hp_operator_t, hp_ntv_kernel
        type(ntv_surface_t) :: sf
        type(ntv_kernel_t) :: k1
        type(hp_operator_t) :: op
        type(hp_field_t) :: f
        type(hp_mesh_t) :: mesh
        real(dp), parameter :: nus(3) = [1e-3_dp, 1e-4_dp, 1e-4_dp], &
            pres(3) = [0.0_dp, 1e-2_dp, 1e-2_dp], omegas(3) = [0.0_dp, 0.0_dp, 3e-3_dp]
        real(dp) :: r, nu, pre, omega
        integer :: status, i, j
        sf%eps = 0.05_dp; sf%q = 2.2_dp; sf%k = 1/(sf%q*6.2_dp); sf%offset = 0
        sf%n = 3; sf%m = [6]; sf%c = [(1e-3_dp, 0.0_dp)]
        op%k = sf%k; op%k_power = 0; op%q = sf%q; op%n = sf%n; op%m = sf%m; op%c = sf%c
        call hp_field([(1 - sf%eps*cos(2*pi*(j - 1)/33), j=1, 33)], 1.0_dp, f, status)
        ok = .true.
        do i = 1, 3
            nu = nus(i); pre = pres(i); omega = omegas(i)
            op%nu = nu; op%pre = pre; op%omega_e = omega
            ! Unchanged to all printed digits at 193 x 900.
            call ntv_kernel(sf, 1.0_dp, nu, omega, pre, 129, 600, k1, status, &
                indicators=.false.)
            if (status /= NTV_OK) error stop 'ntv_kernel failed'
            call hp_default_mesh(f, 6, 16, 4, sqrt(nu), 0.3_dp, mesh)
            call hp_ntv_kernel(f, mesh, op, r, status)
            if (status /= HP_OK) error stop 'hp_ntv_kernel failed'
            ok = close_to('R', r, k1%r, 1e-6_dp) .and. ok
        end do
        ! Boozer streaming k b with the full drift weight: the torque carries
        ! the surface measure 1/b of that streaming.
        sf%boozer_stream = .true.; sf%weight = NTV_WEIGHT_FULL
        op%k_power = 1; op%wsign = 1; op%nu = 1e-3_dp; op%pre = 1e-4_dp; op%omega_e = 0
        call ntv_kernel(sf, 1.0_dp, op%nu, 0.0_dp, op%pre, 65, 80, k1, status, &
            indicators=.false.)
        if (status /= NTV_OK) error stop 'ntv_kernel failed'
        call hp_default_mesh(f, 6, 16, 8, 1e-4_dp, 0.3_dp, mesh)
        call hp_ntv_kernel(f, mesh, op, r, status)
        if (status /= HP_OK) error stop 'hp_ntv_kernel failed'
        ok = close_to('R boozer', r, k1%r, 1e-6_dp) .and. ok
    end function ntv_reference

    !> The stage-1 harmonic (eps 0.05, q 2.2, n 3, m 6) with the operator of
    !> collisionality nu, precession pre and rotation omega_E.
    function stage1(nu, pre, omega) result(op)
        use milca_dk_hp, only: hp_operator_t
        real(dp), intent(in) :: nu, pre, omega
        type(hp_operator_t) :: op
        op%k = 1/(2.2_dp*6.2_dp); op%k_power = 0; op%q = 2.2_dp; op%n = 3
        op%m = [6]; op%c = [(1e-3_dp, 0.0_dp)]
        op%nu = nu; op%pre = pre; op%omega_e = omega
    end function stage1

    real(dp) function stage1_kernel(op, p, n_col, n_bulk, delta) result(r)
        use milca_dk_hp, only: hp_operator_t, hp_ntv_kernel
        type(hp_operator_t), intent(in) :: op
        integer, intent(in) :: p, n_col, n_bulk
        real(dp), intent(in) :: delta
        type(hp_field_t) :: f
        type(hp_mesh_t) :: mesh
        integer :: status, j
        call hp_field([(1 - 0.05_dp*cos(2*pi*(j - 1)/33), j=1, 33)], 1.0_dp, f, status)
        call hp_default_mesh(f, p, n_col, n_bulk, delta, 0.3_dp, mesh)
        call hp_ntv_kernel(f, mesh, op, r, status)
        if (status /= HP_OK) error stop 'hp_ntv_kernel failed'
    end function stage1_kernel

    !> As nu -> 0 with a precession reversal inside the trapped region the
    !> kernel tends to the resonant limit of the bounce-averaged equation
    !> (sbp_limit). It approaches it slowly from below, then overshoots by
    !> about 2 % at nu = 1e-8 (the bounce-averaged equation, solved
    !> numerically, gives 1.021 there); the Legendre reference cannot reach
    !> this collisionality. The check bounds the ratio at nu = 1e-8.
    logical function ntv_superbanana_plateau() result(ok)
        real(dp) :: r, limit
        r = stage1_kernel(stage1(1e-8_dp, 1e-4_dp, 0.0_dp), 6, 16, 8, 1e-4_dp)
        limit = sbp_limit(0.05_dp, stage1(1e-8_dp, 1e-4_dp, 0.0_dp))
        ok = close_to('R/R_sbp at nu = 1e-8', r/limit, 1.0_dp, 0.03_dp)
    end function ntv_superbanana_plateau

    !> R nu at omega_E = pre = 0 from three collisionalities extrapolated as
    !> A + C sqrt(nu) + B nu against the 1/nu limit
    !>   R nu -> (1/2) int dlambda |Q|^2/(lambda I),
    !>   Q(lambda) = int_lambda^{1/b_min} lambda' S_b dlambda',
    !> derived in test_dk_ntv (one_over_nu_limit), where the Legendre
    !> reference reaches 1e-2 from nu = 1e-5 to 1e-6.
    logical function ntv_one_over_nu() result(ok)
        use milca_dk_hp, only: hp_operator_t
        integer, parameter :: ns = 96
        real(dp) :: nu(3), y(3), m(3, 3), lc, lm, xs(ns), ws(ns), xt(ns), wt(ns), &
            s, lam, lam2, outer, limit, qq, fit
        integer :: i, j, piv(3), info
        type(hp_operator_t) :: op
        interface
            subroutine dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
                import :: dp
                integer, intent(in) :: n, nrhs, lda, ldb
                real(dp), intent(inout) :: a(lda, *), b(ldb, *)
                integer, intent(out) :: ipiv(*), info
            end subroutine dgesv
        end interface
        nu = [1e-6_dp, 3e-7_dp, 1e-7_dp]
        do i = 1, 3
            y(i) = nu(i)*stage1_kernel(stage1(nu(i), 0.0_dp, 0.0_dp), 6, 16, 8, &
                0.3_dp*sqrt(nu(i)))
            m(i, :) = [1.0_dp, sqrt(nu(i)), nu(i)]
        end do
        call dgesv(3, 1, m, 3, piv, y, 3, info)
        fit = y(1)
        op = stage1(1.0_dp, 0.0_dp, 0.0_dp)
        lc = 1/(1 + 0.05_dp)
        lm = 1/(1 - 0.05_dp)
        call gauss_legendre_local(ns, xs, ws)
        xs = (xs + 1)/2
        ws = ws/2
        outer = 0
        do i = 1, ns
            lam = lc + (lm - lc)*xs(i)**2
            call gauss_legendre_local(ns, xt, wt)
            qq = 0
            do j = 1, ns
                s = xs(i) + (1 - xs(i))*(xt(j) + 1)/2
                lam2 = lc + (lm - lc)*s**2
                qq = qq + wt(j)*(1 - xs(i))/2*2*(lm - lc)*s*lam2 &
                    *op%n*abs(op%c(1))*well(0.05_dp, op, lam2, 1)
            end do
            outer = outer + ws(i)*2*(lm - lc)*xs(i)*qq**2/(lam*speed(0.05_dp, lam))
        end do
        limit = outer/2
        ok = close_to('R nu at nu -> 0', fit, limit, 1e-3_dp)
    end function ntv_one_over_nu

    !> sum_sigma oint dtheta |xi|/b over the well of b = 1 - eps cos theta.
    real(dp) function speed(eps, lam) result(r)
        real(dp), intent(in) :: eps, lam
        integer, parameter :: nq = 200
        real(dp) :: u(nq), w(nq), t0, th, a
        integer :: k
        call gauss_legendre_local(nq, u, w)
        t0 = acos(max(-1.0_dp, min(1.0_dp, (1 - 1/lam)/eps)))
        r = 0
        do k = 1, nq
            a = (pi/4)*(u(k) + 1)
            th = t0*sin(a)
            r = r + w(k)*(pi/4)*xi_well(eps, lam, t0, th)/(1 - eps*cos(th))*t0*cos(a)
        end do
        r = 4*r
    end function speed

    !> |xi| = sqrt(lam eps (cos theta - cos theta_t)) without cancellation.
    pure real(dp) function xi_well(eps, lam, t0, th)
        real(dp), intent(in) :: eps, lam, t0, th
        xi_well = sqrt(lam*eps*2*sin((t0 + th)/2)*sin((t0 - th)/2))
    end function xi_well


    !> Rapid rotation, omega_E large against streaming and pre = 0:
    !>   R -> (16/15) nu |c|^2/omega_E^2 int_0^2pi dtheta/b^3,
    !>   int_0^2pi dtheta/(1 - eps cos theta)^3 = pi (2 + eps^2)/(1 - eps^2)^(5/2),
    !> approached as omega_E^-2 (derivation in test_dk_ntv, rotation_limit).
    logical function ntv_rotation() result(ok)
        real(dp), parameter :: eps = 0.05_dp, nu = 1e-3_dp, omega = 1.0_dp
        real(dp) :: r, limit
        r = stage1_kernel(stage1(nu, 0.0_dp, omega), 4, 8, 2, 0.1_dp)
        limit = 16/15.0_dp*nu*1e-6_dp/omega**2*pi*(2 + eps**2)/(1 - eps**2)**2.5_dp
        ok = close_to('R/limit at omega_E = 1', r/limit, 1.0_dp, 2e-3_dp)
    end function ntv_rotation

    !> The system recombined from its split terms and solved through the
    !> exported factor against the monoenergetic kernel on the same mesh, at
    !> three (nu, v, omega_E, pre) and with both quasilinear sources, where
    !> f_10 enters through the operators mir and ls.
    logical function terms_recombine() result(ok)
        use milca_dk_hp, only: hp_operator_t, hp_ntv_kernel, hp_terms_t, hp_ntv_terms, &
            hp_system_t, hp_assemble, hp_factor, hp_solve, hp_free
        real(dp), parameter :: par(4, 3) = reshape([1e-2_dp, 1.0_dp, 0.0_dp, 1e-3_dp, &
            1e-3_dp, 2.5_dp, 0.3_dp, -2e-3_dp, 1e-4_dp, 0.5_dp, -0.1_dp, 5e-4_dp], [4, 3])
        type(hp_operator_t) :: op
        type(hp_field_t) :: f
        type(hp_mesh_t) :: mesh
        type(hp_terms_t) :: terms
        type(hp_system_t) :: ax
        real(dp), allocatable :: f10(:)
        complex(dp), allocatable :: rhs(:)
        real(dp) :: r, mu
        integer :: status, i, j, e, n_loc

        ok = .true.
        call hp_field([(1 - 0.05_dp*cos(2*pi*(j - 1)/33), j=1, 33)], 1.0_dp, f, status)
        if (status /= HP_OK) error stop 'hp_field failed'
        call hp_default_mesh(f, 4, 8, 2, 0.1_dp, 0.3_dp, mesh)
        ! nu = 0 and v = 1 of the default operator must not enter the terms.
        op = stage1(0.0_dp, 0.0_dp, 0.0_dp)
        call hp_ntv_terms(f, mesh, op, terms, status)
        if (status /= HP_OK) error stop 'hp_ntv_terms failed'
        do i = 1, 3
            op%nu = par(1, i); op%v = par(2, i); op%omega_e = par(3, i); op%pre = par(4, i)
            call hp_ntv_kernel(f, mesh, op, r, status)
            if (status /= HP_OK) error stop 'hp_ntv_kernel failed'
            ok = close_to('R from terms', solved(terms, op, terms%load3/op%v), r, 1e-10_dp) .and. ok
            ! Without refinement the equilibrated factor alone is accurate here.
            if (i == 1) ok = close_to('R, no refinement', &
                solved(terms, op, terms%load3/op%v, 0), r, 1e-6_dp) .and. ok
        end do

        ! Both sources on a field with b.grad theta = kappa b; f_10 as the
        ! kernel builds it, -2/(v q kappa) times the radial-drive solution.
        op = stage1(1e-3_dp, 1e-3_dp, 0.1_dp)
        op%k = f%kappa; op%k_power = 1; op%v = 1.5_dp
        op%source_dpar = .true.; op%source_dlf10 = .true.
        call hp_ntv_kernel(f, mesh, op, r, status)
        if (status /= HP_OK) error stop 'hp_ntv_kernel with sources failed'
        call hp_ntv_terms(f, mesh, op, terms, status)
        if (status /= HP_OK) error stop 'hp_ntv_terms with sources failed'
        call hp_assemble(f, mesh, op%nu, ax, status)
        if (status == HP_OK) call hp_factor(ax, status)
        allocate (f10(ax%n))
        if (status == HP_OK) call hp_solve(ax, terms%load1, f10, mu, status)
        call hp_free(ax)
        if (status /= HP_OK) error stop 'axisymmetric solve failed'
        f10 = -2/(op%v*op%q*f%kappa)*f10
        rhs = terms%load3/op%v
        n_loc = size(terms%loc, 1)
        do e = 1, size(terms%loc, 2)
            rhs(terms%loc(:, e)) = rhs(terms%loc(:, e)) + matmul(reshape(terms%mir(:, e) &
                + op%nu*terms%ls(:, e), [n_loc, n_loc]), f10(terms%loc(:, e)))
        end do
        ok = close_to('R with sources', solved(terms, op, rhs), r, 1e-10_dp) .and. ok
        ! The f_10 load is not a rounding-level part of this R.
        ok = below('R without the f_10 load, -relative change', &
            -abs(solved(terms, op, terms%load3/op%v)/r - 1), -1e-4_dp) .and. ok
    end function terms_recombine

    !> Re(load4 . x) of the node system of op built from terms, A x = b.
    real(dp) function solved(terms, op, b, refine) result(rt)
        use milca_dk_hp, only: hp_operator_t, hp_terms_t, hp_node_matrix, hp_ntv_factor_t, &
            hp_ntv_factor, hp_ntv_solve, hp_ntv_free
        use fortsparse, only: csc_z_t
        type(hp_terms_t), intent(in) :: terms
        type(hp_operator_t), intent(in) :: op
        complex(dp), intent(in) :: b(:)
        integer, intent(in), optional :: refine
        type(csc_z_t) :: a
        type(hp_ntv_factor_t) :: fac
        complex(dp), allocatable :: x(:)
        integer :: st
        call hp_node_matrix(terms, op%nu, op%v, op%omega_e, op%pre, op%n, op%q, a)
        call hp_ntv_factor(a, fac, st)
        if (st /= HP_OK) error stop 'hp_ntv_factor failed'
        allocate (x(terms%n))
        call hp_ntv_solve(fac, b, x, .false., st, refine)
        call hp_ntv_free(fac)
        if (st /= HP_OK) error stop 'hp_ntv_solve failed'
        rt = real(sum(terms%load4*x), dp)
    end function solved

    pure function csc_apply(a, x) result(y)
        use fortsparse, only: csc_t
        type(csc_t), intent(in) :: a
        real(dp), intent(in) :: x(:)
        real(dp) :: y(a%nrow)
        integer :: j, k
        y = 0
        do j = 1, a%ncol
            do k = a%col_ptr(j), a%col_ptr(j + 1) - 1
                y(a%row_idx(k)) = y(a%row_idx(k)) + a%val(k)*x(j)
            end do
        end do
    end function csc_apply

    pure function csc_column_sums(a) result(y)
        use fortsparse, only: csc_t
        type(csc_t), intent(in) :: a
        real(dp) :: y(a%ncol)
        integer :: j
        y = [(sum(a%val(a%col_ptr(j):a%col_ptr(j + 1) - 1)), j=1, a%ncol)]
    end function csc_column_sums

    logical function below(name, got, bound) result(ok)
        character(*), intent(in) :: name
        real(dp), intent(in) :: got, bound
        ok = got <= bound
        print '(a,es10.2,a,es10.2,a,l2)', name//' ', got, ' bound ', bound, ' pass', ok
    end function below

    !> Mass, pitch moments and their loads on b = 1 + eps cos theta with
    !> b.grad theta = kappa b, measure 1/(kappa b^2), against closed forms:
    !>   sum M = 2 int dtheta/(kappa b^2) = 4 pi (1 - eps^2)^(-3/2)/kappa;
    !>   G_l of P_3(xi)(1 + cos theta) is delta_l3 (1 + cos theta_t) up to
    !>     the interpolation error in theta;
    !>   column t of D_l sums to w_t int P_l dxi/(kappa b^2), and with the
    !>     mesh's own quadrature angles 1^T D_0 G_0 g = d^T g for every g;
    !>   on a node angle the interpolant of xi^3 is xi^3 along the line, whose
    !>     moments are 3/5 (l = 1) and 2/5 (l = 3).
    logical function line_moments() result(ok)
        use milca_dk_hp, only: hp_operator_t, hp_terms_t, hp_ntv_terms, hp_mass, hp_lines_t, &
            hp_lines, hp_line_ops
        use fortsparse, only: csc_t
        use fortnum_quadrature, only: gauss_lobatto_legendre
        real(dp), parameter :: eps = 0.1_dp, kappa = 1.3_dp
        type(hp_operator_t) :: op
        type(hp_field_t) :: f
        type(hp_mesh_t) :: mesh
        type(hp_terms_t) :: terms
        type(hp_lines_t) :: lines
        type(csc_t) :: mass
        type(csc_t), allocatable :: g(:), d(:)
        real(dp), allocatable :: data(:), want(:), b2(:)
        real(dp) :: err(2), e_all, node(4), wn(4)
        integer :: status, j, l, ip, ec

        ok = .true.
        call hp_field([(1 + eps*cos(2*pi*(j - 1)/33), j=1, 33)], kappa, f, status)
        if (status /= HP_OK) error stop 'hp_field failed'
        op%k = kappa; op%k_power = 1; op%n = 1; op%m = [1]; op%c = [(1e-3_dp, 0.0_dp)]
        do ip = 1, 2
            call hp_default_mesh(f, 3*ip, 8, 2, 0.1_dp, 0.3_dp, mesh)
            call hp_ntv_terms(f, mesh, op, terms, status)
            if (status /= HP_OK) error stop 'hp_ntv_terms failed'
            call hp_lines([mesh], lines)
            call hp_line_ops(f, mesh, op, lines, 4, g, d, status)
            if (status /= HP_OK) error stop 'hp_line_ops failed'
            data = (5*terms%xi**3 - 3*terms%xi)/2*(1 + cos(terms%theta))
            err(ip) = 0
            do l = 0, 4
                want = merge(1.0_dp, 0.0_dp, l == 3)*(1 + cos(lines%theta))
                err(ip) = max(err(ip), maxval(abs(csc_apply(g(l), data) - want)))
            end do
        end do
        ok = below('P_3 moments, p = 3', err(1), 1e-3_dp) .and. ok
        ok = below('P_3 moments, p = 6', err(2), 1e-3_dp*err(1)) .and. ok

        ! On the mesh of degree 6.
        call hp_mass(terms, mass)
        ok = close_to('sum of mass', sum(mass%val), 4*pi*(1 - eps**2)**(-1.5_dp)/kappa, &
            1e-12_dp) .and. ok
        ok = close_to('sum of line weights', sum(lines%w), 2*pi, 1e-14_dp) .and. ok
        b2 = (1 + eps*cos(lines%theta))**2
        ok = below('load of P_0', maxval(abs(csc_column_sums(d(0))*kappa*b2/(2*lines%w) &
            - 1)), 1e-12_dp) .and. ok
        ok = below('load of P_2', maxval(abs(csc_column_sums(d(2))*kappa*b2/lines%w)), &
            1e-12_dp) .and. ok
        call random_number(data)
        ok = close_to('1^T D_0 G_0 g', sum(csc_column_sums(d(0))*csc_apply(g(0), data)), &
            sum(terms%density*data), 1e-12_dp) .and. ok

        ! Lines on the interior node angles of every column, degree 3.
        call hp_default_mesh(f, 3, 8, 2, 0.1_dp, 0.3_dp, mesh)
        call hp_ntv_terms(f, mesh, op, terms, status)
        if (status /= HP_OK) error stop 'hp_ntv_terms failed'
        call gauss_lobatto_legendre(4, node, wn)
        lines%theta = [((mesh%theta(ec) + (mesh%theta(ec + 1) - mesh%theta(ec)) &
            *(1 + node(j))/2, j=2, 3), ec=1, size(mesh%theta) - 1)]
        lines%w = 0*lines%theta + 1
        call hp_line_ops(f, mesh, op, lines, 20, g, d, status)
        if (status /= HP_OK) error stop 'hp_line_ops failed'
        e_all = 0
        do l = 0, 20
            want = 0*lines%theta
            if (l == 1) want = 0.6_dp
            if (l == 3) want = 0.4_dp
            e_all = max(e_all, maxval(abs(csc_apply(g(l), terms%xi**3) - want)))
        end do
        ok = below('moments of xi^3 to l = 20', e_all, 1e-10_dp) .and. ok
    end function line_moments

    !> Transfer of exp(xi) cos(theta), given on the nodes of one mesh, to the
    !> nodes of another mesh of the same field with other columns, bulk rows
    !> and separatrix layers, none nested: the error is that of the
    !> interpolant and falls spectrally with p. Onto the mesh of another
    !> surface (theta_x, lambda_c) the nodes are taken at equal angle from the
    !> X point and equal row labels. The constant must be reproduced, equal
    !> meshes must give the identity, and the mass matrix applied to 1 is the
    !> density weight because the basis sums to 1.
    logical function transfer() result(ok)
        use milca_dk_hp, only: hp_operator_t, hp_terms_t, hp_ntv_terms, hp_mass, hp_transfer, &
            hp_mesh_join
        use fortsparse, only: csc_t
        type(hp_operator_t) :: op
        type(hp_field_t) :: f1, f2
        type(hp_mesh_t) :: mesh1, mesh2, mesh3, joined
        type(hp_terms_t) :: t1, t2, t3
        type(csc_t) :: t, mass
        real(dp), allocatable :: one(:), th(:), xi(:)
        real(dp) :: err(3), err_xi(3)
        integer :: status, j, k, p
        logical :: same

        ok = .true.
        call hp_field([(1 + 0.1_dp*cos(2*pi*(j - 1)/33), j=1, 33)], 1.0_dp, f1, status)
        if (status == HP_OK) call hp_field([(1 + 0.2_dp*cos(2*pi*(j - 1)/33 - 1), &
            j=1, 33)], 1.0_dp, f2, status)
        if (status /= HP_OK) error stop 'hp_field failed'
        op%k = 1; op%k_power = 1; op%n = 1; op%m = [1]; op%c = [(1e-3_dp, 0.0_dp)]

        do k = 1, 3
            p = 2*k + 1
            call hp_default_mesh(f1, p, 6, 3, 0.1_dp, 0.3_dp, mesh1)
            call hp_default_mesh(f1, p, 9, 4, 0.01_dp, 0.3_dp, mesh3)
            call hp_ntv_terms(f1, mesh1, op, t1, status)
            if (status == HP_OK) call hp_ntv_terms(f1, mesh3, op, t3, status)
            if (status /= HP_OK) error stop 'hp_ntv_terms failed'
            call hp_transfer(f1, mesh3, mesh1, t3%loc, t1%loc, t, status)
            if (status /= HP_OK) error stop 'hp_transfer failed'
            err(k) = maxval(abs(csc_apply(t, smooth(t1%theta, t1%xi)) &
                - smooth(t3%theta, t3%xi)))
            err_xi(k) = maxval(abs(csc_apply(t, t1%xi) - t3%xi))
        end do
        ok = below('exp(xi) cos(theta), p = 3', err(1), 5e-3_dp) .and. ok
        ok = below('exp(xi) cos(theta), p = 5', err(2), 3e-2_dp*err(1)) .and. ok
        ok = below('exp(xi) cos(theta), p = 7', err(3), 3e-2_dp*err(2)) .and. ok
        ok = below('xi, p = 7', err_xi(3), 3e-2_dp*err_xi(2)) .and. ok

        call hp_default_mesh(f1, 4, 6, 3, 0.05_dp, 0.3_dp, mesh1)
        call hp_default_mesh(f2, 5, 9, 2, 0.01_dp, 0.3_dp, mesh2)
        call hp_ntv_terms(f1, mesh1, op, t1, status)
        if (status == HP_OK) call hp_ntv_terms(f2, mesh2, op, t2, status)
        if (status /= HP_OK) error stop 'hp_ntv_terms failed'
        call hp_transfer(f1, mesh2, mesh1, t2%loc, t1%loc, t, status)
        if (status /= HP_OK) error stop 'hp_transfer failed'
        ok = t%nrow == t2%n .and. t%ncol == t1%n .and. ok
        call labelled_nodes(f1, mesh2, t2, th, xi)
        ok = below('other surface, equal labels', maxval(abs(csc_apply(t, &
            smooth(t1%theta, t1%xi)) - smooth(th, xi))), 5e-4_dp) .and. ok
        allocate (one(t1%n))
        one = 1
        ok = below('constant', maxval(abs(csc_apply(t, one) - 1)), 1e-13_dp) .and. ok

        ! Into the union with a finer mesh of the same surface.
        call hp_default_mesh(f1, 4, 9, 2, 0.01_dp, 0.3_dp, mesh3)
        joined = mesh1
        call hp_mesh_join(joined, mesh3)
        ok = holds(joined%theta, mesh1%theta) .and. holds(joined%theta, mesh3%theta) &
            .and. holds(joined%lam, mesh1%lam) .and. holds(joined%lam, mesh3%lam) &
            .and. holds(joined%zeta, mesh1%zeta) .and. holds(joined%zeta, mesh3%zeta) &
            .and. size(joined%theta) < size(mesh1%theta) + size(mesh3%theta) .and. ok
        print '(a,3i4,a,l2)', 'joined breakpoints ', size(joined%theta), size(joined%lam), &
            size(joined%zeta), ' pass', ok
        call hp_ntv_terms(f1, joined, op, t3, status)
        if (status /= HP_OK) error stop 'hp_ntv_terms failed'
        call hp_transfer(f1, joined, mesh1, t3%loc, t1%loc, t, status)
        if (status /= HP_OK) error stop 'hp_transfer failed'
        ok = below('joined mesh', maxval(abs(csc_apply(t, smooth(t1%theta, t1%xi)) &
            - smooth(t3%theta, t3%xi))), 5e-4_dp) .and. ok

        call hp_transfer(f1, mesh1, mesh1, t1%loc, t1%loc, t, status)
        if (status /= HP_OK) error stop 'hp_transfer failed'
        same = t%nnz == t1%n .and. t%nrow == t1%n .and. t%ncol == t1%n
        if (same) same = all(t%val(:t%nnz) == 1) .and. &
            all(t%row_idx(:t%nnz) == [(k, k=1, t1%n)])
        print '(a,l2)', 'equal meshes give the identity pass', same
        ok = same .and. ok

        call hp_mass(t1, mass)
        ok = below('mass times 1 against density', maxval(abs(csc_apply(mass, one) &
            - t1%density))/maxval(t1%density), 1e-13_dp) .and. ok
    end function transfer

    elemental real(dp) function smooth(th, xi)
        real(dp), intent(in) :: th, xi
        smooth = exp(xi)*cos(th)
    end function smooth

    !> whole increases and holds every point of part, to the tolerance of the union.
    pure logical function holds(whole, part)
        real(dp), intent(in) :: whole(:), part(:)
        integer :: k
        holds = all(whole(2:) > whole(:size(whole) - 1))
        do k = 1, size(part)
            holds = holds .and. minval(abs(whole - part(k))) <= 1e-12_dp
        end do
    end function holds

    !> The nodes of mesh, in the numbering of terms, placed in the geometry f
    !> of another surface: at their angle from the X point and at the xi of
    !> their row labels (passing levels as fractions of lambda_c, trapped
    !> rows zeta), linear between the two curves of a row.
    subroutine labelled_nodes(f, mesh, terms, th, xi)
        use milca_dk_hp, only: hp_terms_t
        use fortnum_quadrature, only: gauss_lobatto_legendre
        type(hp_field_t), intent(in) :: f
        type(hp_mesh_t), intent(in) :: mesh
        type(hp_terms_t), intent(in) :: terms
        real(dp), allocatable, intent(out) :: th(:), xi(:)
        real(dp), allocatable :: val(:)
        integer, allocatable :: kind(:)
        real(dp) :: node(mesh%p + 1), wn(mesh%p + 1), x(2), d
        integer :: n_col, n_lev, ec, er, a1, a2, k, c
        n_col = size(mesh%theta) - 1
        n_lev = size(mesh%lam)
        val = [mesh%lam(:n_lev - 1)/mesh%lam(n_lev)*f%lambda_c, mesh%zeta, &
            mesh%lam(n_lev - 1:1:-1)/mesh%lam(n_lev)*f%lambda_c]
        kind = [(-1, c=1, n_lev - 1), (0, c=1, size(mesh%zeta)), (1, c=1, n_lev - 1)]
        call gauss_lobatto_legendre(mesh%p + 1, node, wn)
        allocate (th(terms%n), xi(terms%n))
        do er = 1, size(val) - 1
            do ec = 1, n_col
                do a2 = 0, mesh%p
                    do a1 = 0, mesh%p
                        k = terms%loc(1 + a1 + (mesh%p + 1)*a2, ec + (er - 1)*n_col)
                        th(k) = f%theta_x + mesh%theta(ec) + (mesh%theta(ec + 1) &
                            - mesh%theta(ec))*(1 + node(a1 + 1))/2 - mesh%theta(n_col + 1)
                        do c = 1, 2
                            call f%level(kind(er + c - 1), val(er + c - 1), th(k), x(c), d)
                        end do
                        xi(k) = x(1) + (x(2) - x(1))*(1 + node(a2 + 1))/2
                    end do
                end do
            end do
        end do
    end subroutine labelled_nodes

    !> R of the stage-1 harmonic (eps 0.05, q 2.2, n 3, m 6) on a uniform
    !> mesh, with its superbanana-plateau limit.
    logical function ntv_scan() result(ok)
        use milca_dk_hp, only: hp_operator_t, hp_ntv_kernel
        character(16) :: arg
        type(hp_operator_t) :: op
        type(hp_field_t) :: f
        type(hp_mesh_t) :: mesh
        real(dp) :: r, eps, delta
        integer :: status, j, p, ncol, nb, n
        integer(8) :: t0, t1, rate
        eps = 0.05_dp
        op%k = 1/(2.2_dp*6.2_dp); op%k_power = 0; op%q = 2.2_dp; op%n = 3
        op%m = [6]; op%c = [(1e-3_dp, 0.0_dp)]
        call get_command_argument(2, arg); read (arg, *) op%nu
        call get_command_argument(3, arg); read (arg, *) op%pre
        call get_command_argument(4, arg); read (arg, *) op%omega_e
        call get_command_argument(5, arg); read (arg, *) p
        call get_command_argument(6, arg); read (arg, *) ncol
        call get_command_argument(7, arg); read (arg, *) nb
        call get_command_argument(8, arg); read (arg, *) delta
        call hp_field([(1 - eps*cos(2*pi*(j - 1)/33), j=1, 33)], 1.0_dp, f, status)
        call hp_default_mesh(f, p, ncol, nb, delta, 0.3_dp, mesh)
        call system_clock(t0, rate)
        call hp_ntv_kernel(f, mesh, op, r, status, n)
        call system_clock(t1)
        print '(a,3es10.2,3i4,es10.2,i9,es22.13,es14.6,f9.2)', 'ntv', op%nu, op%pre, &
            op%omega_e, p, ncol, nb, delta, n, r, r/sbp_limit(eps, op), &
            real(t1 - t0)/rate
        ok = status == HP_OK
    end function ntv_scan

    !> The same harmonic with the Legendre-Fourier kernel milca_dk_ntv.
    logical function ntv_legendre() result(ok)
        use milca_dk_ntv, only: ntv_surface_t, ntv_kernel_t, ntv_kernel, NTV_OK
        character(16) :: arg
        type(ntv_surface_t) :: sf
        type(ntv_kernel_t) :: k
        real(dp) :: nu, pre, omega
        integer :: status, nt, nl
        integer(8) :: t0, t1, rate
        sf%eps = 0.05_dp; sf%q = 2.2_dp; sf%k = 1/(sf%q*6.2_dp); sf%offset = 0
        sf%n = 3; sf%m = [6]; sf%c = [(1e-3_dp, 0.0_dp)]
        call get_command_argument(2, arg); read (arg, *) nu
        call get_command_argument(3, arg); read (arg, *) pre
        call get_command_argument(4, arg); read (arg, *) omega
        call get_command_argument(5, arg); read (arg, *) nt
        call get_command_argument(6, arg); read (arg, *) nl
        call system_clock(t0, rate)
        call ntv_kernel(sf, 1.0_dp, nu, omega, pre, nt, nl, k, status, indicators=.false.)
        call system_clock(t1)
        print '(a,3es10.2,2i6,es22.13,f9.2)', 'legendre', nu, pre, omega, nt, nl, k%r, &
            real(t1 - t0)/rate
        ok = status == NTV_OK
    end function ntv_legendre

    !> Superbanana-plateau limit of the bounce-averaged stage-1 equation
    !>   i n Omega_b g - 2 nu d_lambda(lambda I d_lambda g) = lambda S_b,
    !>   Omega_b = sum_sigma oint dtheta Omega/|xi|, Omega = pre lambda cos theta,
    !>   S_b = sum_sigma oint dtheta e^{-i n q theta} S/|xi|,
    !> as nu -> 0: R = int dlambda lambda Re(S_b^* g) = pi lambda_r^2
    !> |S_b|^2/|n Omega_b'| at the zero lambda_r of Omega_b (offset = 0,
    !> omega_E = 0). Bounce integrals by Gauss quadrature in
    !> theta = theta_t sin(u), which removes the turning-point singularity.
    real(dp) function sbp_limit(eps, op) result(r)
        use milca_dk_hp, only: hp_operator_t
        real(dp), intent(in) :: eps
        type(hp_operator_t), intent(in) :: op
        real(dp) :: lo, hi, mid, h, d
        integer :: it
        lo = 1/(1 + eps)*(1 + 1e-9_dp)
        hi = 1/(1 - eps)*(1 - 1e-9_dp)
        do it = 1, 200
            mid = (lo + hi)/2
            if (lo*well(eps, op, lo, 0)*mid*well(eps, op, mid, 0) <= 0) then
                hi = mid
            else
                lo = mid
            end if
        end do
        h = 1e-6_dp*mid
        d = ((mid + h)*well(eps, op, mid + h, 0) - (mid - h)*well(eps, op, mid - h, 0)) &
            /(2*h)
        r = pi*mid**2*(op%n*abs(op%c(1))*well(eps, op, mid, 1))**2/abs(op%n*op%pre*d)
    end function sbp_limit

    !> 4 int_0^theta_t w(theta) dtheta/|xi|, w = cos theta (kind 0) or
    !> cos((m - n q) theta) (kind 1); |xi|^2 = lam eps (cos - cos theta_t).
    real(dp) function well(eps, op, lam, kind) result(s)
        use milca_dk_hp, only: hp_operator_t
        real(dp), intent(in) :: eps, lam
        type(hp_operator_t), intent(in) :: op
        integer, intent(in) :: kind
        integer, parameter :: nq = 400
        real(dp) :: u(nq), w(nq), t0, th, a
        integer :: j
        call gauss_legendre_local(nq, u, w)
        t0 = acos(max(-1.0_dp, min(1.0_dp, (1 - 1/lam)/eps)))
        s = 0
        do j = 1, nq
            a = (pi/4)*(u(j) + 1)
            th = t0*sin(a)
            if (kind == 0) then
                s = s + w(j)*(pi/4)*cos(th)/xi_well(eps, lam, t0, th)*t0*cos(a)
            else
                s = s + w(j)*(pi/4)*cos((op%m(1) - op%n*op%q)*th) &
                    /xi_well(eps, lam, t0, th)*t0*cos(a)
            end if
        end do
        s = 4*s
    end function well

    subroutine gauss_legendre_local(n, x, w)
        use fortnum_quadrature, only: gauss_legendre
        integer, intent(in) :: n
        real(dp), intent(out) :: x(n), w(n)
        call gauss_legendre(n, x, w)
    end subroutine gauss_legendre_local

    !> One deuterium species at 2 keV, nu_D/v about kappa_coll (1/m) at x
    !> = 1, on a speed grid of nelem elements of order 4 to x = 5.
    subroutine energy_plasma(kappa_coll, nelem, sp, lnl, grid)
        use milca_dk_collop, only: collop_species_t, collop_grid_t, collop_make_grid, &
            collop_nuhat, collop_vth, COLLOP_OK
        real(dp), intent(in) :: kappa_coll
        integer, intent(in) :: nelem
        type(collop_species_t), intent(out) :: sp(1)
        real(dp), intent(out) :: lnl(1, 1)
        type(collop_grid_t), intent(out) :: grid
        real(dp), parameter :: kev = 1.602176634e-16_dp
        integer :: st
        sp(1) = collop_species_t(3.3435837768e-27_dp, 1.0_dp, 3e19_dp, 2*kev)
        lnl = 3*sqrt(pi)/8*collop_vth(sp(1))*kappa_coll/collop_nuhat(sp(1), sp(1), 1.0_dp)
        call collop_make_grid(5.0_dp, nelem, 4, grid, st)
        if (st /= COLLOP_OK) error stop 'collop_make_grid failed'
    end subroutine energy_plasma

    !> The stage-1 harmonic at the speed nodes 2..n: pre = (T/(Z e)) x_k^2,
    !> as the Legendre kernel scales its precession with the node energy.
    function energy_ops(sp, x, omega) result(ops)
        use milca_dk_hp, only: hp_operator_t
        use milca_dk_collop, only: collop_species_t
        type(collop_species_t), intent(in) :: sp
        real(dp), intent(in) :: x(:), omega
        type(hp_operator_t), allocatable :: ops(:, :)
        real(dp), parameter :: charge = 1.602176634e-19_dp
        integer :: k
        allocate (ops(size(x) - 1, 1))
        do k = 1, size(x) - 1
            ops(k, 1) = stage1(0.0_dp, sp%t/(sp%z*charge)*x(k + 1)**2, omega)
        end do
    end function energy_ops

    subroutine stage1_field(f)
        type(hp_field_t), intent(out) :: f
        integer :: status, j
        call hp_field([(1 - 0.05_dp*cos(2*pi*(j - 1)/33), j=1, 33)], 1.0_dp, f, status)
        if (status /= HP_OK) error stop 'hp_field failed'
    end subroutine stage1_field

    !> K and q of the energy-coupled hp kernel on the stage-1 surface, node k
    !> on meshes(mid(k)); J_0 = int dtheta/b^2 = 2 pi (1 - eps^2)^(-3/2).
    subroutine energy_run(kappa_coll, model, meshes, mid, kc, qc, node_r, its, status, &
            max_its, budget, omega, factors, setup_only)
        use milca_dk_hp_energy, only: hp_energy_system_t, hp_energy_setup, &
            hp_energy_coefficients, hp_energy_free
        use milca_dk_collop, only: collop_species_t, collop_grid_t
        real(dp), intent(in) :: kappa_coll
        integer, intent(in) :: model, mid(:)
        type(hp_mesh_t), intent(in) :: meshes(:)
        real(dp), intent(out) :: kc(1, 1, 2), qc(1, 1, 2)
        real(dp), allocatable, intent(out) :: node_r(:, :, :, :)
        integer, intent(out) :: its, status
        integer, intent(in), optional :: max_its
        integer(8), intent(in), optional :: budget
        real(dp), intent(in), optional :: omega
        ! factors: n_factors, factor_bytes and factor_peak after setup.
        integer(8), intent(out), optional :: factors(3)
        logical, intent(in), optional :: setup_only
        type(hp_field_t) :: f
        type(hp_field_t), allocatable :: geom(:, :)
        type(collop_species_t) :: sp(1)
        type(collop_grid_t) :: grid
        type(hp_energy_system_t) :: sys
        real(dp) :: lnl(1, 1), om
        integer(8) :: b

        call stage1_field(f)
        call energy_plasma(kappa_coll, 3, sp, lnl, grid)
        allocate (geom(grid%n - 1, 1))
        geom = f
        b = huge(b)
        if (present(budget)) b = budget
        om = 0
        if (present(omega)) om = omega
        its = 0
        call hp_energy_setup(geom, energy_ops(sp(1), grid%x, om), meshes, &
            reshape(mid, [grid%n - 1, 1]), grid, sp, lnl, model, b, sys, status)
        if (present(factors)) factors = [int(sys%n_factors, 8), sys%factor_bytes, &
            sys%factor_peak]
        if (status /= HP_OK) return
        if (present(setup_only)) then
            if (setup_only) then
                call hp_energy_free(sys)
                return
            end if
        end if
        call hp_energy_coefficients(sys, 2*pi*(1 - 0.05_dp**2)**(-1.5_dp), kc, qc, node_r, &
            its, status, max_iterations=max_its)
        call hp_energy_free(sys)
    end subroutine energy_run

    !> Pitch-angle scattering alone: every node is the monoenergetic kernel
    !> on its own mesh at nu = sum_b nuhat_ab (erf(y) - G(y))/x^3/(v_a x) and
    !> speed v_a x, with the drive tz x^2 (A_1) and tz x^4 (A_2), so that
    !> K(j) = tz/(2 sqrt(pi) J_0) sum_k w_k x_k^4 exp(-x_k^2) tz x_k^(2j) R_k.
    !> Nodes on two band meshes; omega_E and precession on.
    logical function energy_lorentz() result(ok)
        use milca_dk_hp, only: hp_operator_t, hp_ntv_kernel
        use milca_dk_hp_energy, only: HP_NO_MEMORY
        use milca_dk_collop, only: collop_species_t, collop_grid_t, collop_vth, &
            collop_nuhat, COLLOP_LORENTZ
        real(dp), parameter :: charge = 1.602176634e-19_dp, kappa_coll = 1e-2_dp, &
            omega = 300.0_dp
        type(hp_field_t) :: f
        type(hp_mesh_t) :: meshes(2)
        type(hp_operator_t), allocatable :: ops(:, :)
        type(hp_operator_t) :: op
        type(collop_species_t) :: sp(1)
        type(collop_grid_t) :: grid
        real(dp), allocatable :: node_r(:, :, :, :), r(:), x(:), d4(:), d6(:)
        real(dp) :: lnl(1, 1), kc(1, 1, 2), qc(1, 1, 2), vt, tz, y, g, j0, want(3), err(2)
        integer :: its, status, k, mid(9), nx
        integer(8) :: whole(3), half(3)
        logical :: fits

        call stage1_field(f)
        call hp_default_mesh(f, 5, 12, 3, 1e-2_dp, 0.3_dp, meshes(1))
        call hp_default_mesh(f, 5, 16, 4, 3e-3_dp, 0.3_dp, meshes(2))
        mid = [1, 1, 1, 1, 2, 2, 2, 2, 2]
        call energy_run(kappa_coll, COLLOP_LORENTZ, meshes, mid, kc, qc, node_r, its, &
            status, omega=omega)
        if (status /= HP_OK) error stop 'hp energy kernel failed'
        print '(a,i0)', 'GMRES iterations ', its

        call energy_plasma(kappa_coll, 3, sp, lnl, grid)
        nx = grid%n - 1
        x = grid%x(2:)
        vt = collop_vth(sp(1))
        tz = sp(1)%t/(sp(1)%z*charge)
        ops = energy_ops(sp(1), grid%x, omega)
        allocate (r(nx))
        do k = 1, nx
            y = x(k)
            g = (erf(y) - 2*y*exp(-y*y)/sqrt(pi))/(2*y*y)
            op = ops(k, 1)
            op%nu = collop_nuhat(sp(1), sp(1), lnl(1, 1))*(erf(y) - g)/x(k)**3/(vt*x(k))
            op%v = vt*x(k)
            call hp_ntv_kernel(f, meshes(mid(k)), op, r(k), status)
            if (status /= HP_OK) error stop 'hp_ntv_kernel failed'
        end do
        err(1) = maxval(abs(node_r(:, 1, 1, 1) - tz*x**2*r))/maxval(abs(tz*x**2*r))
        err(2) = maxval(abs(node_r(:, 1, 1, 2) - tz*x**4*r))/maxval(abs(tz*x**4*r))
        ok = below('node responses A_1, relative', err(1), 1e-8_dp)
        ok = below('node responses A_2, relative', err(2), 1e-8_dp) .and. ok
        j0 = 2*pi*(1 - 0.05_dp**2)**(-1.5_dp)
        d4 = grid%w(2:)*x**4*exp(-x**2)
        d6 = grid%w(2:)*x**6*exp(-x**2)
        want = tz/(2*sqrt(pi)*j0)*[sum(d4*tz*x**2*r), sum(d4*tz*x**4*r), sum(d6*tz*x**2*r)]
        ok = close_to('K11', kc(1, 1, 1), want(1), 1e-9_dp) .and. ok
        ok = close_to('K12', kc(1, 1, 2), want(2), 1e-9_dp) .and. ok
        ok = close_to('q11', qc(1, 1, 1), want(3), 1e-9_dp) .and. ok

        ! Half the bytes of all factors: refused from one factor per mesh,
        ! before the others are built.
        call energy_run(kappa_coll, COLLOP_LORENTZ, meshes, mid, kc, qc, node_r, its, &
            status, omega=omega, factors=whole, setup_only=.true.)
        call energy_run(kappa_coll, COLLOP_LORENTZ, meshes, mid, kc, qc, node_r, its, &
            status, omega=omega, budget=whole(2)/2, factors=half, setup_only=.true.)
        print '(a,i0,a,i0,a,i0,a,i0,a,i0)', 'budget ', whole(2)/2, ' of ', whole(2), &
            ': status ', status, ', factors ', half(1), ', peak bytes ', half(3)
        fits = status == HP_NO_MEMORY .and. whole(1) == nx .and. half(1) == nx .and. &
            half(2) == 0 .and. half(3) < whole(2)/2
        print '(a,l2)', 'refused within the budget pass', fits
        ok = ok .and. fits

        call energy_run(kappa_coll, COLLOP_LORENTZ, meshes, mid, kc, qc, node_r, its, &
            status, budget=1_8)
        print '(a,i0,a,l2)', 'factors beyond a budget of one byte: status ', status, &
            ' pass', status == HP_NO_MEMORY
        ok = ok .and. status == HP_NO_MEMORY
    end function energy_lorentz

    !> Test-particle operator against the Legendre-Fourier energy-coupled
    !> kernel ntv_energy_coefficients on the same speed grid and surface
    !> (sources off): the same equation in another angular discretisation,
    !> so the difference must fall when both resolutions rise.
    logical function energy_test_particle() result(ok)
        use milca_dk_collop, only: COLLOP_TEST_PARTICLE
        real(dp) :: diff(2), kc(1, 1, 2), qc(1, 1, 2), kl(1, 1, 2), ql(1, 1, 2)
        integer :: ir
        ok = .true.
        do ir = 1, 2
            call energy_pair(ir, kc, qc, kl, ql)
            diff(ir) = max(maxval(abs(kc - kl)/abs(kl)), maxval(abs(qc - ql)/abs(ql)))
            print '(a,i0,a,4es14.6)', 'resolution ', ir, ' hp K11 K12 q11 q12', kc, qc
            print '(a,i0,a,4es14.6)', 'resolution ', ir, ' Legendre         ', kl, ql
        end do
        ok = below('largest relative difference', diff(2), 1e-5_dp)
        ok = below('falls with resolution', diff(2) - diff(1), 0.0_dp) .and. ok
    end function energy_test_particle

    !> The hp and Legendre coefficients of the stage-1 surface at resolution
    !> ir (1 coarse, 2 fine), test-particle operator.
    subroutine energy_pair(ir, kc, qc, kl, ql)
        use milca_dk_ntv, only: ntv_surface_t, ntv_energy_coefficients, NTV_OK
        use milca_dk_collop, only: collop_species_t, collop_grid_t, COLLOP_TEST_PARTICLE
        integer, intent(in) :: ir
        real(dp), intent(out) :: kc(1, 1, 2), qc(1, 1, 2), kl(1, 1, 2), ql(1, 1, 2)
        type(hp_field_t) :: f
        type(hp_mesh_t) :: mesh(1)
        type(ntv_surface_t) :: sf
        type(collop_species_t) :: sp(1)
        type(collop_grid_t) :: grid
        real(dp), allocatable :: node_r(:, :, :, :)
        real(dp) :: lnl(1, 1)
        integer :: its, status
        call stage1_field(f)
        mesh(1) = energy_mesh(f, ir)
        call energy_run(energy_kappa, COLLOP_TEST_PARTICLE, mesh, [1, 1, 1, 1, 1, 1, 1, 1, 1], kc, qc, &
            node_r, its, status)
        print '(a,i0,a,i0)', 'hp status ', status, ' GMRES iterations ', its
        if (status /= HP_OK) error stop 'hp energy kernel failed'
        sf%eps = 0.05_dp; sf%q = 2.2_dp; sf%k = 1/(sf%q*6.2_dp); sf%offset = 0
        sf%n = 3; sf%m = [6]; sf%c = [(1e-3_dp, 0.0_dp)]
        call energy_plasma(energy_kappa, 3, sp, lnl, grid)
        call ntv_energy_coefficients(sf, sp, lnl, grid, COLLOP_TEST_PARTICLE, 0.0_dp, &
            merge(17, 33, ir == 1), merge(24, 48, ir == 1), kl, its, status, ql)
        print '(a,i0,a,i0)', 'Legendre status ', status, ' GMRES iterations ', its
        if (status /= NTV_OK) error stop 'ntv_energy_coefficients failed'
    end subroutine energy_pair

    !> Mesh of resolution ir: 1 coarse, 2 fine, 3 the fine mesh with more
    !> rows, 4 with more columns.
    function energy_mesh(f, ir) result(mesh)
        type(hp_field_t), intent(in) :: f
        integer, intent(in) :: ir
        type(hp_mesh_t) :: mesh
        select case (ir)
        case (1); call hp_default_mesh(f, 4, 8, 2, 3e-2_dp, 0.3_dp, mesh)
        case (2); call hp_default_mesh(f, 6, 12, 3, 1e-2_dp, 0.3_dp, mesh)
        case (3); call hp_default_mesh(f, 6, 12, 5, 5e-3_dp, 0.3_dp, mesh)
        case default; call hp_default_mesh(f, 6, 20, 3, 1e-2_dp, 0.3_dp, mesh)
        end select
    end function energy_mesh

    !> Nodes alternating between a mesh refined in rows and one refined in
    !> columns, against all nodes on the fine common mesh: energy scattering
    !> then couples every node to the other mesh through the transfers.
    logical function energy_two_mesh() result(ok)
        use milca_dk_collop, only: COLLOP_TEST_PARTICLE
        type(hp_field_t) :: f
        type(hp_mesh_t) :: common(1), two(2)
        real(dp), allocatable :: node_r(:, :, :, :)
        real(dp) :: kc(1, 1, 2), qc(1, 1, 2), k2(1, 1, 2), q2(1, 1, 2), diff
        integer :: its, status
        call stage1_field(f)
        common(1) = energy_mesh(f, 2)
        two(1) = energy_mesh(f, 3)
        two(2) = energy_mesh(f, 4)
        call energy_run(energy_kappa, COLLOP_TEST_PARTICLE, common, [1, 1, 1, 1, 1, 1, 1, 1, 1], kc, &
            qc, node_r, its, status)
        if (status /= HP_OK) error stop 'hp energy kernel failed'
        call energy_run(energy_kappa, COLLOP_TEST_PARTICLE, two, [1, 2, 1, 2, 1, 2, 1, 2, 1], k2, q2, &
            node_r, its, status)
        if (status /= HP_OK) error stop 'hp energy kernel on two meshes failed'
        print '(a,i0)', 'GMRES iterations on two meshes ', its
        print '(a,4es14.6)', 'common mesh K11 K12 q11 q12', kc, qc
        print '(a,4es14.6)', 'two meshes  K11 K12 q11 q12', k2, q2
        diff = max(maxval(abs(k2 - kc)/abs(kc)), maxval(abs(q2 - qc)/abs(qc)))
        ok = below('largest relative difference', diff, 5e-4_dp)
    end function energy_two_mesh

    !> The density moment of the collision operator, sum_k v_a d3_k 1^T y_k,
    !> on a random vector over three meshes: pitch-angle scattering and
    !> energy scattering conserve it except for the dropped x = 0 node,
    !> whose test row a_0(1, k') (energy0_w of the collop matrices) gives
    !>   -sum_p sum_k' energy0_w(p, k') <P_{p<-k'} h_k'>_p,
    !> <g>_p the density moment on the mesh of grid node p (built here from
    !> the terms and transfers of each mesh).
    logical function energy_density() result(ok)
        use milca_dk_hp, only: hp_operator_t, hp_terms_t, hp_ntv_terms, hp_transfer
        use milca_dk_hp_energy, only: hp_energy_system_t, hp_energy_setup, &
            hp_energy_collide, hp_energy_free
        use milca_dk_ntv, only: ntv_collision_rows, NTV_OK
        use milca_dk_collop, only: collop_species_t, collop_grid_t, collop_vth, &
            COLLOP_TEST_PARTICLE
        use fortsparse, only: csc_t
        type(hp_field_t) :: f
        type(hp_field_t), allocatable :: geom(:, :)
        type(hp_mesh_t) :: meshes(3)
        type(hp_terms_t) :: terms(3)
        type(hp_operator_t), allocatable :: ops(:, :)
        type(hp_energy_system_t) :: sys
        type(collop_species_t) :: sp(1)
        type(collop_grid_t) :: grid
        type(csc_t) :: t
        real(dp), allocatable :: rows(:, :, :, :, :), ratio(:, :), d3(:), d4(:), d6(:), &
            vt(:), tz(:), e0w(:, :, :), re(:), im(:), moved(:)
        complex(dp), allocatable :: h(:), y(:)
        real(dp) :: lnl(1, 1), got, want, size_
        integer :: status, k, j, p, nx, mid(9), m, mp
        integer, allocatable :: off(:)

        call stage1_field(f)
        do m = 1, 3
            meshes(m) = energy_mesh(f, m + 1)
        end do
        mid = [1, 2, 3, 1, 2, 3, 1, 2, 3]
        call energy_plasma(energy_kappa, 3, sp, lnl, grid)
        nx = grid%n - 1
        ops = energy_ops(sp(1), grid%x, 0.0_dp)
        do m = 1, 3
            call hp_ntv_terms(f, meshes(m), ops(1, 1), terms(m), status)
            if (status /= HP_OK) error stop 'hp_ntv_terms failed'
        end do
        allocate (geom(nx, 1), off(nx + 1))
        geom = f
        call hp_energy_setup(geom, ops, meshes, reshape(mid, [nx, 1]), grid, sp, lnl, &
            COLLOP_TEST_PARTICLE, huge(1_8), sys, status)
        if (status /= HP_OK) error stop 'hp_energy_setup failed'
        off(1) = 0
        do k = 1, nx
            off(k + 1) = off(k) + terms(mid(k))%n
        end do
        ok = sys%n_dof == off(nx + 1)
        allocate (re(off(nx + 1)), im(off(nx + 1)), y(off(nx + 1)))
        call random_number(re)
        call random_number(im)
        h = cmplx(re - 0.5_dp, im - 0.5_dp, dp)
        call hp_energy_collide(sys, h, y)
        call hp_energy_free(sys)

        call ntv_collision_rows(sp, lnl, grid, COLLOP_TEST_PARTICLE, 2, rows, ratio, d3, d4, &
            d6, vt, tz, status, energy0_w=e0w)
        if (status /= NTV_OK) error stop 'ntv_collision_rows failed'
        d3 = grid%w(2:)*grid%x(2:)**3*exp(-grid%x(2:)**2)
        got = 0
        size_ = 0
        do k = 1, nx
            got = got + collop_vth(sp(1))*d3(k)*real(sum(y(off(k) + 1:off(k + 1))), dp)
            size_ = size_ + collop_vth(sp(1))*d3(k)*sum(abs(y(off(k) + 1:off(k + 1))))
        end do
        want = 0
        do p = 1, grid%n
            mp = mid(max(p - 1, 1))
            do j = 1, nx
                if (e0w(p, j, 1) == 0) cycle
                if (mid(j) == mp) then
                    moved = real(h(off(j) + 1:off(j + 1)), dp)
                else
                    call hp_transfer(f, meshes(mp), meshes(mid(j)), terms(mp)%loc, &
                        terms(mid(j))%loc, t, status)
                    if (status /= HP_OK) error stop 'hp_transfer failed'
                    moved = csc_apply(t, real(h(off(j) + 1:off(j + 1)), dp))
                end if
                want = want - e0w(p, j, 1)*sum(terms(mp)%density*moved)
            end do
        end do
        print '(a,es14.6,a,es14.6,a,es9.2)', 'density moment ', got, ' defect ', want, &
            ' scale ', size_
        ok = below('density moment minus defect, relative', abs(got - want)/size_, &
            1e-10_dp) .and. ok
        ! The defect is far above the error: the bound has teeth.
        ok = below('error over defect', abs(got - want)/abs(want), 1e-6_dp) .and. ok
    end function energy_density

    !> GMRES iterations of the test-particle and full kernels at the
    !> collisionality of the comparison case and 100 times above and below
    !> it, and the iteration limit.
    logical function energy_iterations() result(ok)
        use milca_dk_hp, only: HP_NO_CONVERGENCE
        use milca_dk_collop, only: COLLOP_TEST_PARTICLE, COLLOP_FULL
        real(dp), parameter :: factor(3) = [1.0_dp, 1e2_dp, 1e-2_dp]
        type(hp_field_t) :: f
        type(hp_mesh_t) :: mesh(1)
        real(dp), allocatable :: node_r(:, :, :, :)
        real(dp) :: kc(1, 1, 2), qc(1, 1, 2)
        integer :: its, status, i, model
        call stage1_field(f)
        mesh(1) = energy_mesh(f, 2)
        ok = .true.
        do model = COLLOP_TEST_PARTICLE, COLLOP_FULL
            do i = 1, 3
                call energy_run(energy_kappa*factor(i), model, mesh, &
                    [1, 1, 1, 1, 1, 1, 1, 1, 1], kc, qc, node_r, its, status)
                print '(a,i0,a,es9.2,a,i0,a,i0,a,es14.6)', 'model ', model, &
                    ' nu_D/v at x = 1 ', energy_kappa*factor(i), ' status ', status, &
                    ' iterations ', its, ' K11 ', kc(1, 1, 1)
                ok = ok .and. status == HP_OK .and. its <= 60
            end do
        end do
        call energy_run(energy_kappa, COLLOP_TEST_PARTICLE, mesh, [1, 1, 1, 1, 1, 1, 1, 1, 1], kc, qc, &
            node_r, its, status, max_its=2)
        print '(a,i0,a,l2)', 'max_iterations = 2: status ', status, ' pass', &
            status == HP_NO_CONVERGENCE
        ok = ok .and. status == HP_NO_CONVERGENCE
    end function energy_iterations

    !> Deuterium ions and electrons at 2 keV and 3e19 m^-3, one Coulomb
    !> logarithm for all pairs: at equal temperatures nu_D/v at x = 1 is
    !> kappa_coll for both species. Speed grid of nelem elements of order
    !> order to x = 5.
    subroutine two_species(kappa_coll, nelem, order, sp, lnl, grid)
        use milca_dk_collop, only: collop_species_t, collop_grid_t, collop_make_grid, &
            collop_nuhat, collop_vth, COLLOP_OK
        real(dp), intent(in) :: kappa_coll
        integer, intent(in) :: nelem, order
        type(collop_species_t), intent(out) :: sp(2)
        real(dp), intent(out) :: lnl(2, 2)
        type(collop_grid_t), intent(out) :: grid
        real(dp), parameter :: kev = 1.602176634e-16_dp
        integer :: st
        sp(1) = collop_species_t(3.3435837768e-27_dp, 1.0_dp, 3e19_dp, 2*kev)
        sp(2) = collop_species_t(9.1093837015e-31_dp, -1.0_dp, 3e19_dp, 2*kev)
        lnl = 3*sqrt(pi)/8*collop_vth(sp(1))*kappa_coll/collop_nuhat(sp(1), sp(1), 1.0_dp)
        call collop_make_grid(5.0_dp, nelem, order, grid, st)
        if (st /= COLLOP_OK) error stop 'collop_make_grid failed'
    end subroutine two_species

    !> The stage-1 operators of both species at the speed nodes 2..n.
    function two_ops(sp, x) result(ops)
        use milca_dk_hp, only: hp_operator_t
        use milca_dk_collop, only: collop_species_t
        type(collop_species_t), intent(in) :: sp(2)
        real(dp), intent(in) :: x(:)
        type(hp_operator_t), allocatable :: ops(:, :), one(:, :)
        integer :: a
        allocate (ops(size(x) - 1, 2))
        do a = 1, 2
            one = energy_ops(sp(a), x, 0.0_dp)
            ops(:, a) = one(:, 1)
        end do
    end function two_ops

    !> K and q of ions and electrons on the stage-1 surface, all nodes on
    !> mesh.
    subroutine two_run(model, l_field, mesh, kc, qc, its, status)
        use milca_dk_hp_energy, only: hp_energy_system_t, hp_energy_setup, &
            hp_energy_coefficients, hp_energy_free
        use milca_dk_collop, only: collop_species_t, collop_grid_t
        integer, intent(in) :: model, l_field
        type(hp_mesh_t), intent(in) :: mesh
        real(dp), intent(out) :: kc(2, 2, 2), qc(2, 2, 2)
        integer, intent(out) :: its, status
        type(hp_field_t) :: f
        type(hp_field_t), allocatable :: geom(:, :)
        type(collop_species_t) :: sp(2)
        type(collop_grid_t) :: grid
        type(hp_energy_system_t) :: sys
        real(dp), allocatable :: node_r(:, :, :, :)
        real(dp) :: lnl(2, 2)
        integer, allocatable :: mid(:, :)

        call stage1_field(f)
        call two_species(energy_kappa, 3, 4, sp, lnl, grid)
        allocate (geom(grid%n - 1, 2), mid(grid%n - 1, 2))
        geom = f
        mid = 1
        its = 0
        call hp_energy_setup(geom, two_ops(sp, grid%x), [mesh], mid, grid, sp, lnl, model, &
            huge(1_8), sys, status, l_field=l_field)
        if (status /= HP_OK) return
        call hp_energy_coefficients(sys, 2*pi*(1 - 0.05_dp**2)**(-1.5_dp), kc, qc, node_r, &
            its, status)
        call hp_energy_free(sys)
    end subroutine two_run

    !> Full linearised operator, ions and electrons, against the Legendre
    !> kernel ntv_energy_coefficients (sources off, field part to the
    !> Legendre degree of its angular basis): every K and q, including the
    !> cross-species entries, which vanish without the field part. The field
    !> part truncated at l = 3 and 7 agrees within the tolerance, at l = 1
    !> the cross-species entries move by more.
    logical function energy_full() result(ok)
        use milca_dk_ntv, only: ntv_surface_t, ntv_energy_coefficients, NTV_OK
        use milca_dk_collop, only: collop_species_t, collop_grid_t, COLLOP_FULL
        real(dp), parameter :: tol = 2e-2_dp
        type(hp_field_t) :: f
        type(hp_mesh_t) :: mesh
        type(ntv_surface_t) :: sf
        type(collop_species_t) :: sp(2)
        type(collop_grid_t) :: grid
        real(dp) :: lnl(2, 2), k3(2, 2, 2), q3(2, 2, 2), k7(2, 2, 2), q7(2, 2, 2), &
            k1(2, 2, 2), q1(2, 2, 2), kl(2, 2, 2), ql(2, 2, 2), err, cross1
        integer :: its, status, a, b, j

        call stage1_field(f)
        mesh = energy_mesh(f, 1)
        call two_run(COLLOP_FULL, 3, mesh, k3, q3, its, status)
        print '(a,i0,a,i0)', 'l_field 3: status ', status, ' GMRES iterations ', its
        if (status /= HP_OK) error stop 'hp energy kernel failed'
        call two_run(COLLOP_FULL, 7, mesh, k7, q7, its, status)
        print '(a,i0,a,i0)', 'l_field 7: status ', status, ' GMRES iterations ', its
        if (status /= HP_OK) error stop 'hp energy kernel failed'
        call two_run(COLLOP_FULL, 1, mesh, k1, q1, its, status)
        print '(a,i0,a,i0)', 'l_field 1: status ', status, ' GMRES iterations ', its
        if (status /= HP_OK) error stop 'hp energy kernel failed'

        sf%eps = 0.05_dp; sf%q = 2.2_dp; sf%k = 1/(sf%q*6.2_dp); sf%offset = 0
        sf%n = 3; sf%m = [6]; sf%c = [(1e-3_dp, 0.0_dp)]
        call two_species(energy_kappa, 3, 4, sp, lnl, grid)
        call ntv_energy_coefficients(sf, sp, lnl, grid, COLLOP_FULL, 0.0_dp, 33, 48, kl, its, &
            status, ql)
        print '(a,i0,a,i0)', 'Legendre status ', status, ' GMRES iterations ', its
        if (status /= NTV_OK) error stop 'ntv_energy_coefficients failed'

        ok = .true.
        err = 0
        cross1 = huge(1.0_dp)
        do j = 1, 2
            do b = 1, 2
                do a = 1, 2
                    print '(a,3i2,a,4es14.6)', 'K(a,b,j)', a, b, j, ' hp l=3, l=7, l=1, Legendre', &
                        k3(a, b, j), k7(a, b, j), k1(a, b, j), kl(a, b, j)
                    print '(a,3i2,a,4es14.6)', 'q(a,b,j)', a, b, j, ' hp l=3, l=7, l=1, Legendre', &
                        q3(a, b, j), q7(a, b, j), q1(a, b, j), ql(a, b, j)
                    err = max(err, abs(k3(a, b, j) - kl(a, b, j))/abs(kl(a, b, j)), &
                        abs(q3(a, b, j) - ql(a, b, j))/abs(ql(a, b, j)))
                    if (a /= b) cross1 = min(cross1, abs(k1(a, b, j) - k3(a, b, j)) &
                        /abs(k3(a, b, j)), abs(q1(a, b, j) - q3(a, b, j))/abs(q3(a, b, j)))
                end do
            end do
        end do
        ok = below('largest relative difference to Legendre', err, tol)
        ok = below('l_field 7 against 3, relative', max(maxval(abs(k7 - k3)/abs(k3)), &
            maxval(abs(q7 - q3)/abs(q3))), tol) .and. ok
        print '(a,es10.3,a,es10.3)', 'smallest cross-species change at l_field 1 ', cross1, &
            ' above ', tol
        ok = ok .and. cross1 > tol
    end function energy_full

    !> Moments of hp_energy_collide, ions and electrons, on smooth random
    !> vectors, thin geometry: total energy sum_a n_a T_a sum_k x_k^2 v_a d3_k
    !> 1^T y_ak (the dropped node x = 0 has weight zero), and the density of
    !> each species, sum_k v_a d3_k 1^T y_ak, against the defect of the
    !> dropped x = 0 row, on one mesh and on three meshes. Total momentum
    !> sum_a n_a m_a v_a sum_k x_k v_a d3_k w^T y_ak, w the nodal values of xi,
    !> relative to the friction of the test-particle part of the same
    !> vector: what the speed discretisation leaves (lumped pitch-angle
    !> scattering; independent of the mesh, from the Legendre rows) is
    !> removed, the rest falls with the mesh degree.
    logical function collision_identities() result(ok)
        real(dp) :: en, dens(2), mom(3), m3, floor
        integer :: p
        ok = .true.
        do p = 2, 4
            call identities(p, [1, 1, 1], en, dens, mom(p - 1), floor)
            print '(a,i0,a,es10.2,a,2es10.2,a,2es10.2)', 'p = ', p, ': energy ', en, &
                ' density ', dens, ' momentum and its speed part over friction ', &
                mom(p - 1), floor
            ok = below('energy, relative', en, 1e-10_dp) .and. ok
            ok = below('density minus defect, relative', maxval(dens), 1e-10_dp) .and. ok
        end do
        call identities(3, [1, 2, 3], en, dens, m3, floor)
        print '(a,2es10.2)', 'three meshes: density ', dens
        ok = below('density minus defect on three meshes, relative', maxval(dens), 1e-10_dp) &
            .and. ok
        ok = below('momentum p = 3 over p = 2', mom(2)/mom(1), 0.05_dp) .and. ok
        ok = below('momentum p = 4 over p = 3', mom(3)/mom(2), 0.05_dp) .and. ok
    end function collision_identities

    !> The moments of collision_identities at mesh degree p, node (k, a) on
    !> mesh pattern(mod(k + a - 2, 3) + 1): en the relative energy moment,
    !> dens(a) the relative density error after the x = 0 defect; over the
    !> test-particle friction, mom the momentum moment minus that of the
    !> Legendre rows, floor the latter (one mesh only).
    subroutine identities(p, pattern, en, dens, mom, floor)
        use milca_dk_hp, only: hp_operator_t, hp_terms_t, hp_ntv_terms, hp_transfer, &
            hp_lines_t, hp_lines, hp_line_ops
        use milca_dk_hp_energy, only: hp_energy_system_t, hp_energy_setup, &
            hp_energy_collide, hp_energy_free
        use milca_dk_ntv, only: ntv_collision_rows, NTV_OK
        use milca_dk_collop, only: collop_species_t, collop_grid_t, collop_vth, &
            COLLOP_TEST_PARTICLE, COLLOP_FULL
        use fortsparse, only: csc_t
        integer, intent(in) :: p, pattern(3)
        real(dp), intent(out) :: en, dens(2), mom, floor
        integer, parameter :: l_field = 3
        type :: line_ops_t
            type(csc_t), allocatable :: g(:), d(:)
        end type
        type(line_ops_t) :: lo(3)
        type(hp_field_t) :: f
        type(hp_field_t), allocatable :: geom(:, :)
        type(hp_mesh_t) :: meshes(3)
        type(hp_terms_t) :: terms(3)
        type(hp_lines_t) :: lines
        type(hp_operator_t), allocatable :: ops(:, :)
        type(hp_energy_system_t) :: sys
        type(collop_species_t) :: sp(2)
        type(collop_grid_t) :: grid
        type(csc_t) :: t
        real(dp), allocatable :: rows(:, :, :, :, :), ratio(:, :), d3(:), d4(:), d6(:), &
            vt(:), tz(:), e0w(:, :, :), f0w(:, :, :, :, :), c(:, :), hr(:), moved(:), x(:)
        complex(dp), allocatable :: h(:), y3(:), y2(:)
        real(dp) :: lnl(2, 2), want, got, size_, scale, en_sum
        complex(dp) :: mom3, mom2, mom_l, h1
        integer, allocatable :: mid(:, :), off(:, :)
        integer :: status, k, j, a, b, m, mp, q, l, nx, n, model

        call stage1_field(f)
        call hp_default_mesh(f, p, 12, 3, 1e-2_dp, 0.3_dp, meshes(1))
        call hp_default_mesh(f, p, 16, 4, 5e-3_dp, 0.3_dp, meshes(2))
        call hp_default_mesh(f, p, 8, 2, 2e-2_dp, 0.3_dp, meshes(3))
        call two_species(energy_kappa, 3, 4, sp, lnl, grid)
        nx = grid%n - 1
        x = grid%x(2:)
        ops = two_ops(sp, grid%x)
        do m = 1, 3
            call hp_ntv_terms(f, meshes(m), ops(1, 1), terms(m), status)
            if (status /= HP_OK) error stop 'hp_ntv_terms failed'
        end do
        allocate (geom(nx, 2), mid(nx, 2), off(nx + 1, 2))
        geom = f
        n = 0
        do a = 1, 2
            do k = 1, nx
                mid(k, a) = pattern(mod(k + a - 2, 3) + 1)
                off(k, a) = n
                n = n + terms(mid(k, a))%n
            end do
            off(nx + 1, a) = n
        end do
        ! Smooth: random multiples of (1 + c x^2) xi^l exp(i m theta), l, m <=
        ! 2, the same at every resolution.
        call random_seed(size=q)
        call random_seed(put=[(17*k + 3, k=1, q)])
        allocate (h(n), c(9, 6))
        do a = 1, 2
            call random_number(c)
            c = c - 0.5_dp
            do k = 1, nx
                associate (tm => terms(mid(k, a)))
                    h(off(k, a) + 1:off(k + 1, a)) = 0
                    do q = 0, 8
                        h(off(k, a) + 1:off(k + 1, a)) = h(off(k, a) + 1:off(k + 1, a)) &
                            + cmplx(c(q + 1, 1) + c(q + 1, 2)*x(k)**2, c(q + 1, 3) &
                            + c(q + 1, 4)*x(k)**2, dp) &
                            *tm%xi**(q/3)*exp(cmplx(0.0_dp, mod(q, 3)*tm%theta, dp))
                    end do
                end associate
            end do
        end do
        do model = COLLOP_TEST_PARTICLE, COLLOP_FULL
            call hp_energy_setup(geom, ops, meshes, mid, grid, sp, lnl, model, huge(1_8), &
                sys, status, l_field=l_field)
            if (status /= HP_OK) error stop 'hp_energy_setup failed'
            if (sys%n_dof /= n) error stop 'unexpected unknowns'
            if (model == COLLOP_FULL) then
                allocate (y3(n))
                call hp_energy_collide(sys, h, y3)
            else
                allocate (y2(n))
                call hp_energy_collide(sys, h, y2)
            end if
            call hp_energy_free(sys)
        end do

        call ntv_collision_rows(sp, lnl, grid, COLLOP_FULL, l_field + 1, rows, ratio, d3, d4, &
            d6, vt, tz, status, energy0_w=e0w, field0_w=f0w)
        if (status /= NTV_OK) error stop 'ntv_collision_rows failed'
        d3 = grid%w(2:)*x**3*exp(-x**2)
        en = 0
        en_sum = 0
        mom3 = 0
        mom2 = 0
        call hp_lines(meshes, lines)
        do m = 1, 3
            call hp_line_ops(f, meshes(m), ops(1, 1), lines, l_field, lo(m)%g, lo(m)%d, status)
            if (status /= HP_OK) error stop 'hp_line_ops failed'
        end do
        do a = 1, 2
            associate (va => collop_vth(sp(a)))
                got = 0
                size_ = 0
                do k = 1, nx
                    associate (y => y3(off(k, a) + 1:off(k + 1, a)), tm => terms(mid(k, a)))
                        got = got + va*d3(k)*real(sum(y), dp)
                        size_ = size_ + va*d3(k)*sum(abs(y))
                        scale = sp(a)%n*sp(a)%t*x(k)**2*va*d3(k)
                        en = en + scale*real(sum(y), dp)
                        en_sum = en_sum + scale*sum(abs(y))
                        scale = sp(a)%n*sp(a)%m*va*x(k)*va*d3(k)
                        mom3 = mom3 + scale*sum(tm%xi*y)
                        mom2 = mom2 + scale*sum(tm%xi*y2(off(k, a) + 1:off(k + 1, a)))
                    end associate
                end do
                ! The x = 0 row: energy scattering and field part.
                want = 0
                do q = 1, grid%n
                    mp = mid(max(q - 1, 1), a)
                    do j = 1, nx
                        hr = real(h(off(j, a) + 1:off(j + 1, a)), dp)
                        if (mid(j, a) == mp) then
                            moved = hr
                        else
                            call hp_transfer(f, meshes(mp), meshes(mid(j, a)), terms(mp)%loc, &
                                terms(mid(j, a))%loc, t, status)
                            if (status /= HP_OK) error stop 'hp_transfer failed'
                            moved = csc_apply(t, hr)
                        end if
                        want = want - e0w(q, j, a)*sum(terms(mp)%density*moved)
                    end do
                    do b = 1, 2
                        do j = 1, nx
                            hr = real(h(off(j, b) + 1:off(j + 1, b)), dp)
                            do l = 0, l_field
                                ! 1^T D^l of the measure mesh times G^l h.
                                want = want - f0w(q, j, a, b, l)*sum(csc_column_sums( &
                                    lo(mp)%d(l))*csc_apply(lo(mid(j, b))%g(l), hr))
                            end do
                        end do
                    end do
                end do
                dens(a) = abs(got - want)/size_
                print '(a,i0,a,es14.6,a,es14.6)', 'species ', a, ': density moment ', got, &
                    ' x = 0 defect ', want
            end associate
        end do
        en = abs(en)/en_sum
        ! The momentum of the Legendre rows (lumped pitch-angle scattering,
        ! energy and field part at l = 1) on the moments int measure xi h =
        ! (1/3) (1^T D^0)(G^1 h) of thin geometry, measure independent of xi.
        mom_l = 0
        do a = 1, 2
            do k = 1, nx
                scale = sp(a)%n*sp(a)%m*vt(a)*x(k)*vt(a)*d3(k)/ratio(k, a)
                do b = 1, 2
                    do j = 1, nx
                        mp = mid(j, b)
                        h1 = cmplx(sum(csc_column_sums(lo(mp)%d(0))*csc_apply(lo(mp)%g(1), &
                            real(h(off(j, b) + 1:off(j + 1, b)), dp))), &
                            sum(csc_column_sums(lo(mp)%d(0))*csc_apply(lo(mp)%g(1), &
                            aimag(h(off(j, b) + 1:off(j + 1, b))))), dp)/3
                        mom_l = mom_l + scale*rows(k, j, a, b, 1)*h1
                    end do
                end do
            end do
        end do
        mom = abs(mom3 - mom_l)/abs(mom2)
        floor = abs(mom_l)/abs(mom2)
    end subroutine identities

    logical function probe() result(ok)
        use milca_dk_hp_energy, only: hp_energy_system_t, hp_energy_setup, &
            hp_energy_coefficients, hp_energy_free
        use milca_dk_collop, only: collop_species_t, collop_grid_t, COLLOP_TEST_PARTICLE, COLLOP_FULL
        type(collop_species_t) :: s1(1)
        type(collop_grid_t) :: grid
        type(hp_field_t) :: f
        type(hp_field_t), allocatable :: geom(:, :)
        type(hp_energy_system_t) :: sys
        type(hp_mesh_t) :: mesh(1)
        real(dp), allocatable :: node_r(:, :, :, :)
        real(dp) :: l1(1,1), kc(1,1,2), qc(1,1,2), kap
        integer :: its, status, i, model
        character(16) :: arg
        call get_command_argument(2, arg); read (arg, *) kap
        call get_command_argument(3, arg); read (arg, *) i
        call get_command_argument(4, arg); read (arg, *) model
        call stage1_field(f)
        call energy_plasma(kap, 3, s1, l1, grid)
        allocate (geom(grid%n - 1, 1))
        geom = f
        mesh(1) = energy_mesh(f, i)
        call hp_energy_setup(geom, energy_ops(s1(1), grid%x, 0.0_dp), mesh, &
            reshape([(1, its=1, grid%n - 1)], [grid%n - 1, 1]), grid, s1, l1, model, &
            huge(1_8), sys, status, l_field=3)
        print *, 'unknowns', sys%n_dof
        call hp_energy_coefficients(sys, 2*pi*(1 - 0.05_dp**2)**(-1.5_dp), kc, qc, node_r, &
            its, status)
        call hp_energy_free(sys)
        print *, 'hp kappa', kap, 'mesh', i, 'model', model, status, its, kc(1,1,1)
        ok = .true.
    end function probe

    logical function scan() result(ok)
        real(dp), parameter :: eps = 0.1_dp
        character(16) :: arg
        real(dp) :: nu, delta
        integer :: p, ncol, nb
        type(mono_coefficients_t) :: c
        call get_command_argument(2, arg); read (arg, *) nu
        call get_command_argument(3, arg); read (arg, *) p
        call get_command_argument(4, arg); read (arg, *) ncol
        call get_command_argument(5, arg); read (arg, *) nb
        call get_command_argument(6, arg); read (arg, *) delta
        c = hp_solve(eps, nu, p, ncol, nb, delta)
        print '(a,es9.1,3i4,es9.1,4es20.12)', 'scan', nu, p, ncol, nb, delta, c%d11/nu, c%d31, c%d13, c%d33*nu
        ok = .true.
    end function scan

    logical function adapt_scan() result(ok)
        use milca_dk_hp, only: hp_estimate, hp_refine
        real(dp), parameter :: eps = 0.1_dp
        character(16) :: arg
        real(dp) :: nu, tol, value, eta
        real(dp), allocatable :: eta_el(:)
        integer :: p, it, status, n
        integer(8) :: t0, t1, rate
        type(hp_field_t) :: f
        type(hp_mesh_t) :: mesh, next
        call get_command_argument(2, arg); read (arg, *) nu
        call get_command_argument(3, arg); read (arg, *) p
        call get_command_argument(4, arg); read (arg, *) tol
        call hp_field(field(33, eps), 1.0_dp, f, status)
        call hp_default_mesh(f, p, 4, 2, sqrt(nu), 0.3_dp, mesh)
        call system_clock(t0, rate)
        do it = 0, 30
            call hp_estimate(f, mesh, nu, 1, 2, value, eta, eta_el, status, n)
            call system_clock(t1)
            print '(a,i3,i8,3es22.13,f10.2)', 'adapt', it, n, value, eta, value + eta, &
                real(t1 - t0)/rate
            if (abs(eta) <= tol*abs(value)) exit
            call hp_refine(mesh, eta_el, 0.5_dp, next)
            mesh = next
        end do
        ok = status == HP_OK
    end function adapt_scan

    logical function mono_scan() result(ok)
        real(dp), parameter :: eps = 0.1_dp
        character(16) :: arg
        real(dp) :: nu
        integer :: nt, nl, status
        type(mono_coefficients_t) :: c
        call get_command_argument(2, arg); read (arg, *) nu
        call get_command_argument(3, arg); read (arg, *) nt
        call get_command_argument(4, arg); read (arg, *) nl
        call mono_coefficients(field(nt, eps), 1.0_dp, nu, nl, c, status)
        print '(a,es9.1,2i6,4es20.12)', 'mono', nu, nt, nl, c%d11/nu, c%d31, c%d13, c%d33*nu
        ok = status == MONO_OK
    end function mono_scan

end program test_dk_hp
