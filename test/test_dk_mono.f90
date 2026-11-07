!> Monoenergetic coefficients against limits derived independently by
!! bounce or flux-surface averaging of the same equation (u = I/Omega_0 = 1,
!! v = 1, b.grad theta = kappa b, flux-surface weight 1/b^2):
!!   collisionless  3 nu D33/<b^2> -> f_c,  3 D31 -> f_t = 1 - f_c,
!!   banana         D11/nu -> <b^-2>/3 - (1/4) int_0^{1/b_max} lam dlam/<sqrt(1 - lam b)>,
!!   Pfirsch-Schlueter  D11/nu -> c_PS from the theta problem of the P_0, P_1 rows,
!!   plateau        D11 -> pi kappa eps^2/16 for b = 1/(1 + eps cos theta), eps -> 0,
!!   Onsager        D31 = -D13, and the density-source multipliers vanish,
!! with f_c = (3/4) <b^2> int_0^{1/b_max} lam dlam/<sqrt(1 - lam b)>.
!! The finite-collisionality corrections of the first two limits scale like
!! sqrt(nu); two collisionalities eliminate them (Richardson).
program test_dk_mono
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_dk_mono, only: mono_coefficients_t, mono_coefficients, MONO_OK, &
        mono_coefficients_jvp, mono_estimate_t, mono_estimate, mono_factor_t, mono_factor, &
        mono_solve, mono_apply, mono_drives
    use fortnum_quadrature, only: gauss_legendre_ab
    use kin6d_boozer_file, only: boozer_file_t, boozer_surface_t, load_boozer_file, &
        boozer_surface, BOOZER_OK
    implicit none
    real(dp), parameter :: pi = 4*atan(1.0_dp)
    character(32) :: check
    logical :: ok

    call get_command_argument(1, check)
    select case (trim(check))
    case ('collisionless'); ok = collisionless()
    case ('banana'); ok = banana()
    case ('pfirsch_schlueter'); ok = pfirsch_schlueter()
    case ('plateau'); ok = plateau()
    case ('onsager'); ok = onsager()
    case ('neo2_tok02'); ok = neo2_tok02()
    case ('jvp'); ok = jvp()
    case ('estimate'); ok = estimate()
    case ('factor_api'); ok = factor_api()
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

    function solve(n, eps, nu, n_leg) result(c)
        integer, intent(in) :: n, n_leg
        real(dp), intent(in) :: eps, nu
        type(mono_coefficients_t) :: c
        integer :: status
        call mono_coefficients(field(n, eps), 1.0_dp, nu, n_leg, c, status)
        if (status /= MONO_OK) error stop 'mono_coefficients failed'
    end function solve

    ! Flux-surface average with weight 1/b^2 on a fine uniform grid.
    pure real(dp) function average(x, b)
        real(dp), intent(in) :: x(:), b(:)
        average = sum(x/b**2)/sum(1/b**2)
    end function average

    ! int_0^{1/b_max} lam dlam/<sqrt(1 - lam b)>, with lam = lam_c (1 - t^2)
    ! so that the square-root edge at lam_c becomes smooth in t.
    real(dp) function passing_integral(eps) result(q)
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

    real(dp) function passing_fraction(eps) result(f_c)
        real(dp), intent(in) :: eps
        real(dp) :: b(40001)
        b = field(size(b), eps)
        f_c = 0.75_dp*average(b**2, b)*passing_integral(eps)
    end function passing_fraction

    logical function close_to(name, got, want, tol) result(ok)
        character(*), intent(in) :: name
        real(dp), intent(in) :: got, want, tol
        ok = abs(got - want) <= tol*abs(want)
        print '(a,es14.6,a,es14.6,a,es9.2,a,l2)', name//' ', got, ' oracle ', want, &
            ' rel ', abs(got - want)/abs(want), ' pass', ok
    end function close_to


    ! Value at nu -> 0 from two collisionalities with X = X0 + c sqrt(nu).
    pure real(dp) function richardson(x1, nu1, x2, nu2)
        real(dp), intent(in) :: x1, nu1, x2, nu2
        richardson = (x2*sqrt(nu1) - x1*sqrt(nu2))/(sqrt(nu1) - sqrt(nu2))
    end function richardson

    logical function collisionless() result(ok)
        ! Resolutions converged to 1e-5 relative at these collisionalities.
        real(dp), parameter :: eps = 0.1_dp, nu1 = 1e-4_dp, nu2 = 2.5e-5_dp
        type(mono_coefficients_t) :: c1, c2
        real(dp) :: b2, f_c
        c1 = solve(129, eps, nu1, 400)
        c2 = solve(257, eps, nu2, 800)
        b2 = average(field(40001, eps)**2, field(40001, eps))
        f_c = passing_fraction(eps)
        ok = close_to('3 nu D33/<b^2>', richardson(3*nu1*c1%d33/b2, nu1, &
            3*nu2*c2%d33/b2, nu2), f_c, 1e-4_dp)
        ok = close_to('3 D31', richardson(3*c1%d31, nu1, 3*c2%d31, nu2), &
            1 - f_c, 1e-4_dp) .and. ok
    end function collisionless

    logical function banana() result(ok)
        ! Resolutions converged to 1e-5 relative at these collisionalities.
        real(dp), parameter :: eps = 0.1_dp, nu1 = 1e-4_dp, nu2 = 2.5e-5_dp
        type(mono_coefficients_t) :: c1, c2
        real(dp) :: b(40001), want
        b = field(size(b), eps)
        want = average(1/b**2, b)/3 - passing_integral(eps)/4
        c1 = solve(129, eps, nu1, 400)
        c2 = solve(257, eps, nu2, 800)
        ok = close_to('D11/nu', richardson(c1%d11/nu1, nu1, c2%d11/nu2, nu2), &
            want, 1e-4_dp)
    end function banana

    ! The P_0 row gives a_1 = 1/b + c b, the O(nu) P_1 row k d_theta a_0 =
    ! -nu a_1, whose periodicity fixes c = -mean(b^-2), so that
    ! D11/nu -> <s_0 A> with d_theta A = -(b^-2 - mean(b^-2))/kappa and the
    ! P_0 drive s_0 = -(2/3) kappa b_theta/b.
    logical function pfirsch_schlueter() result(ok)
        real(dp), parameter :: eps = 0.1_dp, nu = 1e4_dp
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
        c = solve(33, eps, nu, 50)
        ok = close_to('D11/nu', c%d11/nu, average(s0*big_a, b), 1e-6_dp)
    end function pfirsch_schlueter

    ! The plateau value is approached linearly in nu*; eps = 1e-5 keeps the
    ! banana correction eps^(3/2)/nu* below 1e-5.
    logical function plateau() result(ok)
        real(dp), parameter :: eps = 1e-5_dp, nu1 = 1e-2_dp, nu2 = 1e-3_dp
        type(mono_coefficients_t) :: c1, c2
        real(dp) :: d_p
        d_p = pi*eps**2/16
        c1 = solve(129, eps, nu1, 400)
        c2 = solve(129, eps, nu2, 400)
        ok = close_to('D11/D_p', (nu1*c2%d11 - nu2*c1%d11)/((nu1 - nu2)*d_p), &
            1.0_dp, 1e-4_dp)
    end function plateau

    ! D31 = -D13 from the parity of the drives and of the operator, and the
    ! density-source multipliers vanish with the density moments of both drives;
    ! both measured against sqrt(D11 D33), the scale of the off-diagonal.
    logical function onsager() result(ok)
        type(mono_coefficients_t) :: c
        real(dp) :: scale
        integer :: i
        ok = .true.
        do i = -4, 3
            c = solve(33, 0.1_dp, 10.0_dp**i, 400)
            scale = sqrt(c%d11*c%d33)
            ok = ok .and. abs(c%d31 + c%d13) <= 1e-10_dp*scale .and. &
                maxval(abs(c%border)) <= 1e-10_dp*scale
            print '(a,es9.1,a,es9.2,a,es9.2)', 'nu ', 10.0_dp**i, &
                ' |D31 + D13|/sqrt(D11 D33) ', abs(c%d31 + c%d13)/scale, &
                ' |mu|/sqrt(D11 D33) ', maxval(abs(c%border))/scale
        end do
    end function onsager

    ! b(theta) and kappa = iota B0/(G + iota I) of an axisymmetric surface.
    subroutine surface_field(path, s, n, b, kappa)
        character(*), intent(in) :: path
        real(dp), intent(in) :: s
        integer, intent(in) :: n
        real(dp), intent(out) :: b(n), kappa
        type(boozer_file_t) :: bf
        type(boozer_surface_t) :: surf
        real(dp) :: b0, theta
        integer :: status, j, k
        call load_boozer_file(path, bf, status)
        if (status /= BOOZER_OK) error stop 'cannot read the Boozer file'
        call boozer_surface(bf, s, surf, status)
        if (status /= BOOZER_OK) error stop 'surface outside the file'
        if (any(surf%n /= 0)) error stop 'surface is not axisymmetric'
        b0 = sum(surf%bmn, mask=surf%m == 0)
        do j = 1, n
            theta = 2*pi*(j - 1)/n
            b(j) = 0
            do k = 1, size(surf%m)
                b(j) = b(j) + surf%bmn(k)*cos(surf%m(k)*theta)
            end do
        end do
        b = b/b0
        kappa = abs(surf%iota*b0/(surf%g_cov + surf%iota*surf%i_cov))
    end subroutine surface_field

    ! NEO-2-PAR in Lorentz mode (built from main on 2026-11-05) on tok02 at
    ! s = 0.50781, at four times its default resolution (nstep, eta_part,
    ! ripple_solver_accurfac x 4, bsfunc errors / 4). Its collpar per cm is
    ! 2 nu_D/v, and its gamma_E is (pi/4) 3 nu D33/<b^2>. R13 =
    ! D13 D31/(D11 D33) is free of the drive normalisations. The tolerance is
    ! NEO-2's resolution error: its last refinement still moves R13 by 3.5e-4.
    ! At collpar 3.5e-2 R13 ~ 6e-7 is below NEO-2's resolution and is not used.
    logical function neo2_tok02() result(ok)
        real(dp), parameter :: collpar(3) = [3.5154e-6_dp, 3.5154e-4_dp, 3.5154e-2_dp]
        real(dp), parameter :: gamma_e(3) = [0.357024_dp, 0.577930_dp, 0.785227_dp]
        real(dp), parameter :: r13(2) = [-1.068652_dp, -0.2858438_dp]
        character(512) :: path
        real(dp) :: b(65), kappa, nu
        type(mono_coefficients_t) :: c
        integer :: i, status
        call get_command_argument(2, path)
        call surface_field(trim(path), 0.50781_dp, size(b), b, kappa)
        ok = .true.
        do i = 1, 3
            nu = collpar(i)*100/2
            call mono_coefficients(b, kappa, nu, 200, c, status)
            if (status /= MONO_OK) error stop 'mono_coefficients failed'
            ok = close_to('gamma_E', pi/4*3*nu*c%d33/average(b**2, b), gamma_e(i), &
                2e-3_dp) .and. ok
            if (i <= 2) ok = close_to('R13', c%d13*c%d31/(c%d11*c%d33), r13(i), &
                2e-3_dp) .and. ok
        end do
    end function neo2_tok02


    !> Tangents through the factorisation against centred differences in nu,
    !> kappa and the field shape b -> b + h cos(2 theta).
    logical function jvp() result(ok)
        integer, parameter :: n = 33, n_leg = 40
        real(dp), parameter :: eps = 0.2_dp, nu = 1e-2_dp, kappa = 0.7_dp, h = 1e-6_dp
        real(dp) :: b(n), shape(n), dirs(3, 3)
        type(mono_coefficients_t) :: c, dc, up, down
        integer :: j, k, status
        b = field(n, eps)
        shape = [(cos(4*pi*(j - 1)/n), j=1, n)]
        ! Rows: d(nu), d(kappa), d(shape amplitude).
        dirs = reshape([1, 0, 0, 0, 1, 0, 0, 0, 1], [3, 3])
        ok = .true.
        do k = 1, 3
            call mono_coefficients_jvp(b, dirs(3, k)*shape, kappa, dirs(2, k), nu, &
                dirs(1, k), n_leg, c, dc, status)
            if (status /= MONO_OK) error stop 'jvp failed'
            call mono_coefficients(b + h*dirs(3, k)*shape, kappa + h*dirs(2, k), &
                nu + h*dirs(1, k), n_leg, up, status)
            call mono_coefficients(b - h*dirs(3, k)*shape, kappa - h*dirs(2, k), &
                nu - h*dirs(1, k), n_leg, down, status)
            ok = close_to('dD11', dc%d11, (up%d11 - down%d11)/(2*h), 1e-6_dp) .and. ok
            ok = close_to('dD31', dc%d31, (up%d31 - down%d31)/(2*h), 1e-6_dp) .and. ok
            ok = close_to('dD13', dc%d13, (up%d13 - down%d13)/(2*h), 1e-6_dp) .and. ok
            ok = close_to('dD33', dc%d33, (up%d33 - down%d33)/(2*h), 1e-6_dp) .and. ok
        end do
    end function jvp

    !> The discretisation indicator against the error measured with a reference
    !> at four times the resolution (effectivity index within [0.5, 2] where the
    !> error is resolved), and the algebraic indicator at the level of rounding.
    logical function estimate() result(ok)
        real(dp), parameter :: eps = 0.2_dp, nu = 3e-3_dp
        type(mono_coefficients_t) :: c, ref
        type(mono_estimate_t) :: e
        integer :: status
        real(dp) :: eff
        call mono_estimate(field(9, eps), 1.0_dp, nu, 24, c, e, status)
        if (status /= MONO_OK) error stop 'estimate failed'
        ref = solve(37, eps, nu, 96)
        eff = abs(e%discretisation%d11)/abs(ref%d11 - c%d11)
        print '(a,es10.2,a,es10.2,a,f6.2)', 'D11 indicator ', e%discretisation%d11, &
            ' error ', ref%d11 - c%d11, ' effectivity ', eff
        ok = eff >= 0.5_dp .and. eff <= 2
        eff = abs(e%discretisation%d33)/abs(ref%d33 - c%d33)
        print '(a,es10.2,a,es10.2,a,f6.2)', 'D33 indicator ', e%discretisation%d33, &
            ' error ', ref%d33 - c%d33, ' effectivity ', eff
        ok = ok .and. eff >= 0.5_dp .and. eff <= 2
        print '(a,2es10.2)', 'algebraic indicators D11, D33 ', e%algebraic%d11, e%algebraic%d33
        ok = ok .and. abs(e%algebraic%d11) <= 1e-12_dp*abs(c%d11) &
            .and. abs(e%algebraic%d33) <= 1e-12_dp*abs(c%d33)
    end function estimate

    !> The public factorisation, bordered solve and operator action agree:
    !> A x + mu e_0 = s and <w x_0> = 0 for both drives.
    logical function factor_api() result(ok)
        type(mono_factor_t) :: f
        real(dp), allocatable :: s(:, :, :), x(:, :, :), r(:, :)
        real(dp) :: mu(2), err
        integer :: status, k
        call mono_factor(field(33, 0.2_dp), 0.7_dp, 1e-2_dp, 40, f, status)
        if (status /= MONO_OK) error stop 'mono_factor failed'
        s = mono_drives(f)
        call mono_solve(f, s, [0.0_dp, 0.0_dp], x, mu, status)
        if (status /= MONO_OK) error stop 'mono_solve failed'
        err = 0
        do k = 1, 2
            r = mono_apply(f, 1e-2_dp, x(:, :, k)) - s(:, :, k)
            r(:, 0) = r(:, 0) + mu(k)
            err = max(err, maxval(abs(r))/maxval(abs(s(:, :, k))), &
                abs(sum(f%p%weight*x(:, 0, k)))/maxval(abs(x(:, 0, k))))
        end do
        print '(a,es10.2)', 'residual of the bordered solve through mono_apply:', err
        ok = err < 1e-12_dp
    end function factor_api

end program test_dk_mono
