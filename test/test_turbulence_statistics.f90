!> Spectral estimator against synthetic Ornstein-Uhlenbeck modes with known
!> Lorentzian spectrum, correlation time and correlation length: exact
!> Parseval identity, recovery within the reported standard errors, the
!> symmetry of the full-plane spectrum and the series standard error. A
!> spectrum normalised without the window power is rejected.
program test_turbulence_statistics
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use fortnum_rng, only: rng_t, rng_seed, rng_normal
    use fortnum_status, only: fortnum_status_t
    use milca_turbulence_statistics
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp)
    integer :: failures

    failures = 0
    call ornstein_uhlenbeck()
    call series()
    if (failures > 0) error stop 'test_turbulence_statistics failed'
    write (*, '(a)') 'test_turbulence_statistics passed'
contains
    subroutine ornstein_uhlenbeck()
        integer, parameter :: nkx = 5, nky = 3, nseg = 256, samples = 256*130
        real(dp), parameter :: dt = 0.25_dp, gamma = 0.2_dp, omega0 = 0.8_dp
        real(dp), parameter :: sigma2 = 0.7_dp, kx(nkx) = 0.3_dp*[-2, -1, 0, 1, 2]
        real(dp), parameter :: ky(nky) = 0.5_dp*[0, 1, 2]
        type(turbulence_accumulator_t) :: acc
        type(turbulence_spectrum_t) :: s
        type(rng_t) :: rng
        type(fortnum_status_t) :: status
        complex(dp) :: phi(nkx, nky), n(nkx, nky), x, decay
        complex(dp), allocatable :: series(:)
        real(dp) :: a, b, weighted, cell, lorentz, tau_exact, peak
        real(dp) :: worst, cross_re, cross_im, expected
        integer :: t, w, pos, segments, i

        call rng_seed(rng, 20261111_int64, status)
        call turbulence_accumulator_init(acc, kx, ky, dt, nseg)
        decay = exp(cmplx(-gamma*dt, -omega0*dt, dp))
        call rng_normal(rng, a); call rng_normal(rng, b)
        x = sqrt(sigma2/2)*cmplx(a, b, dp)
        allocate (series(samples))
        do t = 1, samples
            call rng_normal(rng, a); call rng_normal(rng, b)
            x = x*decay + sqrt(sigma2/2*(1 - exp(-2*gamma*dt)))*cmplx(a, b, dp)
            ! One mode (kx, ky) = (0.3, 0.5); n = (2 - i) phi.
            phi = 0
            phi(4, 2) = x
            n = (2.0_dp, -1.0_dp)*phi
            call turbulence_accumulate(acc, phi, n)
            series(t) = x
        end do
        call turbulence_spectrum(acc, s)
        segments = samples/(nseg/2) - 1
        if (s%segments /= segments) failures = failures + 1
        if (s%batches < 8) failures = failures + 1
        cell = (kx(2) - kx(1))*(ky(2) - ky(1))*(s%omega(2) - s%omega(1))
        ! Independent window-weighted variance over the same segments; the
        ! mirror mode doubles it.
        weighted = 0
        do pos = 1, samples - nseg + 1, nseg/2
            do t = 1, nseg
                weighted = weighted + 2*sin(pi*(t - 1)/nseg)**4*abs(series(pos + t - 1))**2
            end do
        end do
        weighted = weighted/(segments*sum(sin(pi*[(i, i=0, nseg - 1)]/nseg)**4))
        call check('Parseval: spectrum sum equals window-weighted variance', &
                   abs(sum(s%s_phi_phi)*cell - weighted)/weighted, 1e-12_dp)
        ! Mutant: a periodogram divided by nseg^2 instead of nseg times the window power.
        call check_fails('spectrum without window-power normalisation is rejected', &
                         abs(sum(s%s_phi_phi)*cell*sum(sin(pi*[(i, i=0, nseg - 1)]/nseg)**4)/nseg &
                             - weighted)/weighted, 0.5_dp)
        call check('variance equals spectrum sum', abs(s%variance%value - sum(s%s_phi_phi)*cell) &
                   /weighted, 1e-12_dp)
        call within('variance of both mirror modes', s%variance, 2*sigma2, 4.0_dp)

        ! Lorentzian at (kx, ky) = (0.3, 0.5) and its mirror at (-0.3, -0.5, -omega).
        worst = 0; peak = 0
        do w = 1, nseg
            lorentz = sigma2/pi*gamma/((s%omega(w) - omega0)**2 + gamma**2)/((kx(2) - kx(1))*0.5_dp)
            peak = max(peak, lorentz)
            worst = max(worst, abs(s%s_phi_phi(4, 4, w) - lorentz))
            if (abs(s%s_phi_phi(2, 2, nseg + 2 - max(w, 2)) - s%s_phi_phi(4, 4, max(w, 2))) > 0) &
                failures = failures + 1
        end do
        ! Peak estimate has relative standard error of about sqrt(11/(9 segments)).
        call check('Lorentzian spectrum', worst/peak, 5*sqrt(11.0_dp/(9*segments)) + 0.02_dp)
        if (abs(sum(s%s_phi_phi) - 2*sum(s%s_phi_phi(4, 4, :))) > 1e-12_dp*sum(s%s_phi_phi)) &
            failures = failures + 1
        ! Cross spectrum <n conjg(phi)> = (2 - i) S_phi_phi.
        cross_re = sum(s%s_n_phi_re(4, 4, :)); cross_im = sum(s%s_n_phi_im(4, 4, :))
        expected = sum(s%s_phi_phi(4, 4, :))
        call check('cross spectrum', abs(cross_re - 2*expected) + abs(cross_im + expected), &
                   1e-12_dp*expected)
        call check('cross spectrum of the mirror mode', &
                   abs(sum(s%s_n_phi_im(2, 2, :)) - expected), 1e-12_dp*expected)

        ! Correlation scales: |c(s)| = exp(-gamma s); c(s) = exp(-gamma s) cos(omega0 s);
        ! c(r) = cos(k r).
        call within('mode correlation time 1/gamma', s%tau_mode, 1/gamma, 4.0_dp)
        tau_exact = crossing_time(gamma, omega0)
        call within('single-point correlation time', s%tau_euler, tau_exact, 4.0_dp)
        call check('correlation length in x', abs(s%ell_x%value - acos(exp(-1.0_dp))/0.3_dp), &
                   0.02_dp)
        call check('correlation length in y', abs(s%ell_y%value - acos(exp(-1.0_dp))/0.5_dp), &
                   0.02_dp)
        call within('Kubo number', s%kubo, tau_exact*0.5_dp*sqrt(2*sigma2) &
                    /(acos(exp(-1.0_dp))/0.3_dp), 4.0_dp)
        if (.not. (s%tau_mode%error > 0 .and. s%tau_mode%error < 0.1_dp/gamma)) &
            failures = failures + 1
    end subroutine

    !> First s with exp(-gamma s) cos(omega0 s) = 1/e, by bisection.
    real(dp) function crossing_time(gamma, omega0)
        real(dp), intent(in) :: gamma, omega0
        real(dp) :: lo, hi, mid
        integer :: it
        lo = 0; hi = min(1/gamma, 0.5_dp*pi/omega0)
        do it = 1, 60
            mid = (lo + hi)/2
            if (exp(-gamma*mid)*cos(omega0*mid) > exp(-1.0_dp)) then
                lo = mid
            else
                hi = mid
            end if
        end do
        crossing_time = (lo + hi)/2
    end function

    !> AR(1) series: mean, standard error and integrated autocorrelation time.
    subroutine series()
        integer, parameter :: n = 200000
        real(dp), parameter :: rho = 0.9_dp, mean = 0.6_dp, dt = 0.5_dp
        type(turbulence_scalar_t) :: scalar
        type(rng_t) :: rng
        type(fortnum_status_t) :: status
        real(dp), allocatable :: x(:)
        real(dp) :: a, exact_error
        integer :: t

        allocate (x(n))
        call rng_seed(rng, 7_int64, status)
        x(1) = mean
        do t = 2, n
            call rng_normal(rng, a)
            x(t) = mean + rho*(x(t - 1) - mean) + sqrt(1 - rho**2)*a
        end do
        scalar = turbulence_series_scalar(x, dt)
        ! Variance of the mean: (1 + rho)/(1 - rho)/n for unit variance.
        exact_error = sqrt((1 + rho)/(1 - rho)/n)
        call check('series standard error', abs(scalar%error/exact_error - 1), 0.1_dp)
        call check('series mean within four standard errors', &
                   abs(scalar%value - mean)/scalar%error, 4.0_dp)
        call check('integrated autocorrelation time', &
                   abs(scalar%tau/(dt*(1 + rho)/(2*(1 - rho))) - 1), 0.1_dp)
    end subroutine

    subroutine within(name, scalar, exact, sigmas)
        character(*), intent(in) :: name
        type(turbulence_scalar_t), intent(in) :: scalar
        real(dp), intent(in) :: exact, sigmas
        write (*, '(a,3es12.4)') name//': value, error, exact ', scalar%value, scalar%error, exact
        if (.not. scalar%error > 0) failures = failures + 1
        call check(name//' within standard errors', abs(scalar%value - exact)/scalar%error, sigmas)
        call check(name//' within 5 %', abs(scalar%value/exact - 1), 0.05_dp)
    end subroutine

    subroutine check(name, value, tolerance)
        character(*), intent(in) :: name
        real(dp), intent(in) :: value, tolerance
        if (value <= tolerance) then
            write (*, '(a,a,es10.3)') 'ok   ', name//': ', value
        else
            write (*, '(a,a,es10.3,a,es10.3)') 'FAIL ', name//': ', value, ' > ', tolerance
            failures = failures + 1
        end if
    end subroutine

    subroutine check_fails(name, value, tolerance)
        character(*), intent(in) :: name
        real(dp), intent(in) :: value, tolerance
        if (value > tolerance) then
            write (*, '(a,a,es10.3)') 'ok   ', name//': ', value
        else
            write (*, '(a,a,es10.3)') 'FAIL ', name//': ', value
            failures = failures + 1
        end if
    end subroutine
end program test_turbulence_statistics
