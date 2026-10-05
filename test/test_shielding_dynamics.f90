program test_shielding_dynamics
    use, intrinsic :: iso_fortran_env, only: dp => real64, qp => real128
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use kin6d_shielding_dynamics, only: shielding_state_t, init_shielding, &
     advance_shielding, sample_shielding, shielding_coefficients, shielding_plot_padding
    use kin6d_shielding_dynamics_leaves, only: shielding_mode, shielding_mode_jvp, &
                                               shielding_mode_vjp, shielding_equilibrium
    implicit none
    call continuous_oracle()
    call derivative_oracle()
    call failure_oracle()
   write (*, '(a)') 'Shielding continuous, Poisson, derivative and failure oracles pass'
contains
    subroutine continuous_oracle()
        type(shielding_state_t) :: state, coarse, fine
      real(dp), parameter :: source(4) = [0.125_dp, -0.0625_dp, 0.03125_dp, 0.015625_dp]
        real(dp) :: value(3), lo(3), hi(3), vl(3), ll(3), ul(3), vr(3), lr(3), ur(3)
       real(dp) :: x, times(3), padding(3), coarse_a(4), fine_a(4), r(4), src(4), steady
        real(qp) :: exact(3), a, kk, xx, mu, phi, error_coarse, error_fine
        logical :: ok
        integer :: n, j, q, p
        times = [0.0_dp, 0.125_dp, 0.875_dp]
        call init_shielding(state, source, ok)
        if (.not. ok) error stop 'init'
        do n = 1, size(times)
            call advance_shielding(state, times(n), 0.0005_dp, ok)
            if (.not. ok) error stop 'advance'
            do p = 0, 64
                x = -3.125_dp + real(p, dp)/10.0_dp
                call sample_shielding(state, x, value, lo, hi, ok)
                if (.not. ok) error stop 'sample'
                exact = 0.0_qp
                xx = real(x, qp)
                do j = 1, 4
                    kk = real(j, qp)
                    mu = kk**2 + 1.0_qp
                    a = real(source(j), qp)/mu*(1.0_qp - exp(-mu*real(times(n), qp)))
                    phi = (real(source(j), qp) - a)/kk**2
                    exact = exact + [a*cos(kk*xx), phi*cos(kk*xx), kk*phi*sin(kk*xx)]
                end do
                do q = 1, 3
                    if (exact(q) < real(lo(q), qp)) error stop 'lower containment'
                    if (exact(q) > real(hi(q), qp)) error stop 'upper containment'
                end do
                if (p < 64) then
                    call sample_shielding(state, x - 0.05_dp, vl, ll, ul, ok)
                    call sample_shielding(state, x + 0.05_dp, vr, lr, ur, ok)
               call shielding_plot_padding(state, x - 0.05_dp, x + 0.05_dp, padding, ok)
                    if (.not. ok) error stop 'chord padding'
                    do q = 1, 3
                        if (exact(q) < real((ll(q) + lr(q))/2 - padding(q), qp)) &
                            error stop 'chord lower'
                        if (exact(q) > real((ul(q) + ur(q))/2 + padding(q), qp)) &
                            error stop 'chord upper'
                    end do
                end if
            end do
        end do
        call init_shielding(coarse, source, ok)
        call init_shielding(fine, source, ok)
        call advance_shielding(coarse, 0.25_dp, 0.002_dp, ok)
        call advance_shielding(fine, 0.25_dp, 0.001_dp, ok)
        call shielding_coefficients(coarse, src, coarse_a, r, ok)
        call shielding_coefficients(fine, src, fine_a, r, ok)
        error_coarse = 0.0_qp
        error_fine = 0.0_qp
        do j = 1, 4
            mu = real(j, qp)**2 + 1.0_qp
            a = real(source(j), qp)/mu*(1.0_qp - exp(-mu/4.0_qp))
            error_coarse = error_coarse + abs(real(coarse_a(j), qp) - a)
            error_fine = error_fine + abs(real(fine_a(j), qp) - a)
            call shielding_equilibrium(source(j), real(j, dp), steady)
            if (abs((real(j, dp)**2 + 1)*steady - source(j)) > 1e-15_dp) &
                error stop 'Helmholtz equilibrium'
        end do
        if (error_fine >= 0.6_qp*error_coarse) error stop 'Euler refinement'
    end subroutine

    subroutine derivative_oracle()
        real(dp) :: primal(3), tangent(3), plus(3), minus(3), gradient(4)
        real(dp), parameter :: a = 0.125_dp, k = 2.0_dp, rho = 0.5_dp, x = 0.75_dp
        real(dp), parameter :: direction(4) = [0.5_dp, -0.25_dp, 0.125_dp, 0.375_dp]
        real(dp), parameter :: weight(3) = [-0.25_dp, 0.5_dp, 0.75_dp], eps = 1e-6_dp
        call shielding_mode(a, k, rho, x, primal(1), primal(2), primal(3))
       call shielding_mode_jvp(a, k, rho, x, direction(1), direction(2), direction(3), &
                                direction(4), tangent(1), tangent(2), tangent(3))
        call shielding_mode_vjp(a, k, rho, x, weight(1), weight(2), weight(3), &
                                gradient(1), gradient(2), gradient(3), gradient(4))
        call shielding_mode(a + eps*direction(1), k + eps*direction(2), &
                rho + eps*direction(3), x + eps*direction(4), plus(1), plus(2), plus(3))
        call shielding_mode(a - eps*direction(1), k - eps*direction(2), &
             rho - eps*direction(3), x - eps*direction(4), minus(1), minus(2), minus(3))
        if (maxval(abs(tangent - (plus - minus)/(2*eps))) > 2e-10_dp) &
            error stop 'independent mode finite difference'
        if (abs(dot_product(weight, tangent) - dot_product(direction, gradient)) > &
            1e-14_dp) error stop 'JVP/VJP pairing'
    end subroutine

    subroutine failure_oracle()
        type(shielding_state_t) :: state
        real(dp) :: val(3), lo(3), hi(3), nan
        logical :: ok
        nan = ieee_value(0.0_dp, ieee_quiet_nan)
        call init_shielding(state, [0.01_dp, huge(0.0_dp)], ok)
        if (.not. ok) error stop 'finite initial huge source'
        call advance_shielding(state, 1.0_dp, 1.0_dp, ok)
        if (ok) error stop 'overflow accepted'
        call sample_shielding(state, 0.0_dp, val, lo, hi, ok)
        if (ok) error stop 'partially mutated failed packet usable'
        call init_shielding(state, [nan], ok)
        if (ok) error stop 'NaN source accepted'
        call init_shielding(state, [0.1_dp], ok)
        call advance_shielding(state, 0.1_dp, 0.001_dp, ok, certify=.false.)
        if (.not. ok) error stop 'candidate-only advance'
        call sample_shielding(state, 0.0_dp, val, lo, hi, ok)
        if (ok) error stop 'candidate-only packet certified'
    end subroutine
end program test_shielding_dynamics
