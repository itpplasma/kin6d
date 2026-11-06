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
    use kin6d_dk_mono, only: mono_coefficients_t, mono_coefficients, MONO_OK
    use fortnum_quadrature, only: gauss_legendre_ab
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

end program test_dk_mono
