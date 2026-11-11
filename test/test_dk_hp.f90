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
