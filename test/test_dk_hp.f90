!> Spectral elements on the orbit-class mesh against the Legendre-Fourier
!! reference kin6d_dk_mono and against limits derived independently by
!! averaging the same equation (as in test_dk_mono): collisionless, banana
!! and Pfirsch-Schlueter, at collisionalities below the reach of the
!! reference. Field b = 1/(1 + eps cos theta), kappa = 1.
!! The modes scan, adapt and mono are timing drivers, not registered tests:
!!   test_dk_hp scan nu p n_col n_bulk delta   uniform mesh
!!   test_dk_hp adapt nu p tol                 adaptive D31 to relative tol
!!   test_dk_hp mono nu n_theta n_leg          Legendre-Fourier reference
program test_dk_hp
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_dk_hp, only: hp_field_t, hp_field, hp_mesh_t, hp_default_mesh, &
        hp_coefficients, HP_OK
    use kin6d_dk_mono, only: mono_coefficients_t, mono_coefficients, MONO_OK
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
        ok = abs(got - want) <= tol*abs(want)
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
        use kin6d_dk_hp, only: hp_estimate, hp_refine
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
        use kin6d_dk_hp, only: hp_estimate, hp_refine
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
