!> Spectral elements on the orbit-class mesh against the Legendre-Fourier
!! reference milca_dk_mono and against limits derived independently by
!! averaging the same equation (as in test_dk_mono): collisionless, banana
!! and Pfirsch-Schlueter, at collisionalities below the reach of the
!! reference. Field b = 1/(1 + eps cos theta), kappa = 1.
!! The modes scan, adapt and mono are timing drivers, not registered tests:
!!   test_dk_hp scan nu p n_col n_bulk delta   uniform mesh
!!   test_dk_hp adapt nu p tol                 adaptive D31 to relative tol
!!   test_dk_hp mono nu n_theta n_leg          Legendre-Fourier reference
!> The thin field with the fallbacks of the radial deposit forced: points
!! with cos theta > 1/2 have no orbit, orbits of sign +1 have no dwell time
!! and orbits of sign -1 spend all their time in box away.
module hp_fallback_field
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use milca_dk_hp, only: hp_field_t, hp_operator_t, HP_OK
    implicit none
    type, extends(hp_field_t) :: fallback_t
        integer :: away = 1
    contains
        procedure :: point_source => fallback_source
        procedure :: dwell => fallback_dwell
    end type
contains
    subroutine fallback_source(g, th, xi, op, src, collide, drive_r, drive_e, density, &
            phase, precess, drive, torque, level, valid)
        class(fallback_t), intent(in) :: g
        real(dp), intent(in) :: th, xi
        type(hp_operator_t), intent(in) :: op
        complex(dp), intent(out) :: src
        real(dp), intent(out) :: collide, drive_r, drive_e, density, phase, precess, drive, &
            torque
        real(dp), intent(out), optional :: level
        logical, intent(out), optional :: valid
        call g%hp_field_t%point_source(th, xi, op, src, collide, drive_r, drive_e, &
            density, phase, precess, drive, torque, level, valid)
        if (present(valid)) valid = cos(th) <= 0.5_dp
    end subroutine fallback_source

    subroutine fallback_dwell(g, kind, level, panels, edges, time, status)
        class(fallback_t), intent(in) :: g
        integer, intent(in) :: kind
        real(dp), intent(in) :: level, panels(:), edges(:)
        real(dp), intent(out) :: time(:)
        integer, intent(out) :: status
        select case (kind)
        case (1)
            time = 0
            status = HP_OK + 1
        case (-1)
            time = 0
            time(g%away) = 1
            status = HP_OK
        case default
            call g%hp_field_t%dwell(kind, level, panels, edges, time, status)
        end select
    end subroutine fallback_dwell
end module hp_fallback_field

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
    case ('ntv_deposit_conserves'); ok = ntv_deposit_conserves()
    case ('ntv_deposit_fallbacks'); ok = ntv_deposit_fallbacks()
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

    !> The radial deposit of a thin kernel: the band torques and the box
    !> shares both sum to the kernel functional, which is unchanged by asking
    !> for them, and thin orbits leave all of it in the box of offset zero.
    !> Sums are compared on the scale of the band torques, which can cancel.
    !> hp_ntv_estimate deposits the same solution.
    logical function ntv_deposit_conserves() result(ok)
        use milca_dk_hp, only: hp_operator_t, hp_ntv_kernel, hp_ntv_estimate, hp_bands_t
        type(hp_operator_t) :: op
        type(hp_field_t) :: f
        type(hp_mesh_t) :: mesh
        type(hp_bands_t) :: bands
        real(dp), parameter :: edges(4) = [-1.0_dp, -0.1_dp, 0.2_dp, 1.0_dp]
        real(dp) :: r, r0, dep(3), dep_e(3), scale, eta
        real(dp), allocatable :: eta_el(:)
        integer :: status, j
        op = stage1(1e-3_dp, 0.3_dp, 0.0_dp)
        call hp_field([(1 - 0.05_dp*cos(2*pi*(j - 1)/33), j=1, 33)], 1.0_dp, f, status)
        call hp_default_mesh(f, 4, 16, 6, 1e-3_dp, 0.3_dp, mesh)
        call hp_ntv_kernel(f, mesh, op, r0, status)
        if (status /= HP_OK) error stop 'hp_ntv_kernel failed'
        call hp_ntv_estimate(f, mesh, op, r, eta, eta_el, status, edges=edges, dep=dep_e)
        if (status /= HP_OK) error stop 'hp_ntv_estimate with deposit failed'
        call hp_ntv_kernel(f, mesh, op, r, status, edges=edges, dep=dep, bands=bands)
        if (status /= HP_OK) error stop 'hp_ntv_kernel with deposit failed'
        scale = sum(abs(bands%torque))
        print '(a,es20.12,a,3es11.3,a,es10.3)', 'R ', r, ' sums-R, off-centre ', &
            sum(dep) - r, sum(bands%torque) - r, maxval(abs(dep([1, 3]))), ' scale ', scale
        print '(a,i0,a,i0,a,es10.3)', 'bands ', size(bands%torque), ' unplaced ', &
            bands%unplaced, ' stray ', bands%stray
        ok = r == r0 .and. all(dep_e == dep) .and. abs(sum(dep) - r) < 1e-10_dp*scale &
            .and. abs(sum(bands%torque) - r) < 1e-10_dp*scale .and. all(dep([1, 3]) == 0) &
            .and. bands%unplaced == 0 .and. bands%stray == 0 .and. bands%home == 2 &
            .and. count(abs(bands%torque) > 1e-6_dp*scale) > size(bands%torque)/2
    end function ntv_deposit_conserves

    !> Torque of points without an orbit (stray) and of bands without time
    !> (unplaced, here all passing rows of sign +1) goes to the home box,
    !> bands of sign -1 go to box 3, and the deposit still sums to R.
    logical function ntv_deposit_fallbacks() result(ok)
        use milca_dk_hp, only: hp_operator_t, hp_ntv_kernel, hp_bands_t
        use hp_fallback_field, only: fallback_t
        type(hp_operator_t) :: op
        type(fallback_t) :: f
        type(hp_mesh_t) :: mesh
        type(hp_bands_t) :: bands
        real(dp), parameter :: edges(4) = [-1.0_dp, -0.1_dp, 0.2_dp, 1.0_dp]
        real(dp) :: r, dep(3), scale, home, away
        integer :: status, j
        op = stage1(1e-3_dp, 0.3_dp, 0.0_dp)
        call hp_field([(1 - 0.05_dp*cos(2*pi*(j - 1)/33), j=1, 33)], 1.0_dp, f%hp_field_t, &
            status)
        f%away = 3
        call hp_default_mesh(f, 4, 16, 6, 1e-3_dp, 0.3_dp, mesh)
        call hp_ntv_kernel(f, mesh, op, r, status, edges=edges, dep=dep, bands=bands)
        if (status /= HP_OK) error stop 'hp_ntv_kernel with deposit failed'
        scale = sum(abs(bands%torque)) + abs(bands%stray)
        home = bands%stray + sum(bands%torque, mask=bands%kind >= 0)
        away = sum(bands%torque, mask=bands%kind < 0)
        print '(a,es10.3,a,i0,a,i0,a,i0)', 'stray/scale ', abs(bands%stray)/scale, &
            ' unplaced ', bands%unplaced, ' of ', count(bands%kind == 1), ' home ', bands%home
        print '(a,4es11.3)', 'sum-R, dep(1), home and away errors / scale ', &
            (sum(dep) - r)/scale, dep(1)/scale, (dep(2) - home)/scale, (dep(3) - away)/scale
        ok = abs(bands%stray) > 1e-3_dp*scale .and. bands%home == 2 &
            .and. bands%unplaced == count(bands%kind == 1) .and. bands%unplaced > 0 &
            .and. abs(sum(dep) - r) < 1e-10_dp*scale .and. dep(1) == 0 &
            .and. abs(dep(2) - home) < 1e-12_dp*scale .and. abs(dep(3) - away) < 1e-12_dp*scale &
            .and. abs(away) > 1e-3_dp*scale
    end function ntv_deposit_fallbacks

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
