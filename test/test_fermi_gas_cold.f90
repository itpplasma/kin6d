! Degenerate ideal Fermi gas: candidate accuracy and certificate coverage at
! small theta = kB T/EF against an independent binary128 Sommerfeld oracle.
! Oracle: integral H(eps)(-df/deps) deps = sum_n c_2n theta^2n H^(2n)(mu)
! for power laws H = a eps^p (c0 = 1, c2 = pi^2/6, c4 = 7 pi^4/360,
! c6 = 31 pi^6/15120); exponentially small lower-limit terms ~exp(-mu/theta)
! are below 1e-11 and the theta^8 remainder below 3e-10 for theta <= 0.025.
! Density 3 int r^2 f dr = 1 uses G = eps^(3/2); pressure G = (2/5) eps^(5/2);
! susceptibility chi = 3 int r^2 f(1-f)/theta dr uses H = (3/2) eps^(1/2);
! screening ratio k_s/k_TF = sqrt(2 chi/3).
program test_fermi_gas_cold
    use, intrinsic :: iso_fortran_env, only: dp => real64, qp => real128
    use fortnum_interval, only: interval_t
    use kin6d_fermi_gas, only: fermi_workspace_t, fermi_state_t, &
        initialize_fermi_workspace, solve_fermi_gas
    implicit none
    real(dp), parameter :: thetas(7) = [0.025_dp, 0.02_dp, 0.01_dp, 3.0e-3_dp, &
                                        1.0e-3_dp, 3.0e-4_dp, 1.0e-4_dp]
    type(fermi_workspace_t) :: work
    type(fermi_state_t) :: state
    real(qp) :: mu, pressure, chi, ratio, errors(4)
    integer :: j
    logical :: ok

    do j = 1, size(thetas)
        call initialize_fermi_workspace(work, 16384, ok)
        call require(ok, 'workspace')
        call sommerfeld_state(real(thetas(j), qp), mu, pressure, chi, ratio)
        call solve_fermi_gas(work, thetas(j), state, ok, certify=.false.)
        call require(ok, 'cold candidate')
        errors = abs([real(state%mu, qp) - mu, real(state%pressure, qp) - pressure, &
                      real(state%susceptibility, qp) - chi, &
                      real(state%screening_ratio, qp) - ratio]/ &
                     [mu, pressure, chi, ratio])
        print '(a,es9.2,a,4es10.2)', 'theta=', thetas(j), &
            ' candidate rel. error mu/P/chi/k_s:', errors
        call require(maxval(errors) < 1.0e-9_qp, 'candidate matches Sommerfeld')
        call solve_fermi_gas(work, thetas(j), state, ok)
        call require(ok .and. state%status == 'certified-continuous', &
                     'cold certificate does not fail closed')
        call require(encloses(state%mu_bound, mu) .and. &
                     encloses(state%pressure_bound, pressure) .and. &
                     encloses(state%susceptibility_bound, chi) .and. &
                     encloses(state%screening_ratio_bound, ratio), &
                     'cold certificate contains Sommerfeld values')
    end do
    print '(a)', 'PASS: cold Fermi candidate and certificate vs Sommerfeld oracle'

contains

    subroutine sommerfeld_state(theta, mu, pressure, chi, ratio)
        real(qp), intent(in) :: theta
        real(qp), intent(out) :: mu, pressure, chi, ratio
        real(qp) :: lower, upper
        integer :: iteration
        lower = 0.5_qp
        upper = 1.5_qp
        do iteration = 1, 200
            mu = (lower + upper)/2
            if (sommerfeld(1.0_qp, 1.5_qp, mu, theta) < 1.0_qp) then
                lower = mu
            else
                upper = mu
            end if
        end do
        mu = (lower + upper)/2
        pressure = sommerfeld(0.4_qp, 2.5_qp, mu, theta)
        chi = sommerfeld(1.5_qp, 0.5_qp, mu, theta)
        ratio = sqrt(2*chi/3)
    end subroutine sommerfeld_state

    !> sum_n c_2n theta^2n d^2n/deps^2n (a eps^p) at eps = mu, n = 0..3.
    pure real(qp) function sommerfeld(a, p, mu, theta) result(value)
        real(qp), intent(in) :: a, p, mu, theta
        real(qp) :: pi, coefficient(0:3), falling
        integer :: n, m
        pi = acos(-1.0_qp)
        coefficient = [1.0_qp, pi**2/6, 7*pi**4/360, 31*pi**6/15120]
        value = 0
        do n = 0, 3
            falling = 1
            do m = 0, 2*n - 1
                falling = falling*(p - m)
            end do
            value = value + coefficient(n)*theta**(2*n)*a*falling*mu**(p - 2*n)
        end do
    end function sommerfeld

    logical function encloses(bound, expected) result(inside)
        type(interval_t), intent(in) :: bound
        real(qp), intent(in) :: expected
        inside = real(bound%lo, qp) <= expected .and. expected <= real(bound%hi, qp)
    end function encloses

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(*), intent(in) :: message
        if (condition) return
        print '(a)', 'FAIL: '//message
        error stop 1
    end subroutine require
end program test_fermi_gas_cold
