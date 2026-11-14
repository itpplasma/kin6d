!> Cylindrical kink and tearing against closed forms: uniform current (exact
!> xi = x^(m-1) and the classic m-1 < nq_a < m window) and a smoothed uniform
!> current channel whose Delta' converges to the step-current formula.
program test_cylinder_stability
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use milca_cylinder_stability
    implicit none
    integer, parameter :: nx = 2001
    real(dp) :: q(nx), cur(nx), dcur(nx), x, dw, g, fa, lam, ref, xs, rd, rd2
    real(dp) :: q0, x1, w, i_enc
    integer :: i, rc, m

    ! Uniform current, no wall: dW = (2F/q)(n q - m + 1) with F = n/m - 1/q.
    do m = 2, 3
        q = 0.97_dp*real(m, dp) - 0.4_dp
        call cyl_kink_energy(nx, q, m, 1, 1e6_dp, dw, g, rc)
        if (rc /= CYL_OK) error stop 'uniform current kink status'
        fa = 1.0_dp/m - 1/q(1)
        lam = 1
        ref = fa*fa*(m-1) + 2*fa/q(1) + (1 + m*lam)*fa*fa
        if (abs(dw - ref) > 1e-6_dp*abs(ref)) error stop 'uniform current dW'
        if (dw >= 0 .or. g <= 0) error stop 'window m-1<q<m must be unstable'
    end do
    q = 2.2_dp
    call cyl_kink_energy(nx, q, 2, 1, 1e6_dp, dw, g, rc)
    if (rc /= CYL_OK .or. dw <= 0) error stop 'q_a above m/n must be stable'
    ! A close wall stabilises the external kink.
    q = 1.7_dp
    call cyl_kink_energy(nx, q, 2, 1, 1.05_dp, dw, g, rc)
    if (rc /= CYL_OK .or. dw <= 0) error stop 'close wall must stabilise'
    ! Peaked current with q = 2 inside: ideal problem is resonant.
    call peaked(1.0_dp, 3.5_dp)
    call cyl_kink_energy(nx, q, 2, 1, 1.2_dp, dw, g, rc)
    if (rc /= CYL_RESONANT_INSIDE) error stop 'resonant m=2'

    ! Smoothed current channel: J = 2/q0 for x < x1, converging to the step.
    q0 = 1.5_dp; x1 = 0.5_dp
    do i = 1, nx
        x = real(i-1, dp)/real(nx-1, dp)
        w = 0.004_dp
        i_enc = channel(x, x1, w)
        if (i == 1) then
            q(i) = q0
        else
            q(i) = q0*x*x/i_enc
        end if
    end do
    call cyl_current_from_q(nx, q, cur, dcur)
    call cyl_tearing_delta_prime(nx, q, dcur, 2, 1, 1.5_dp, 1e-4_dp, xs, rd, rc)
    if (rc /= CYL_OK) error stop 'step tearing status'
    ref = cyl_tearing_delta_prime_step(2, 1, q0, x1, xs, 1.5_dp)
    print '(a,3f12.5)', 'step channel x_s, r_s Delta'' shooting, closed form', xs, rd, ref
    if (abs(rd - ref) > 0.03_dp*abs(ref)) error stop 'step-current Delta'''
    call cyl_tearing_delta_prime(nx, q, dcur, 2, 1, 1.5_dp, 5e-5_dp, xs, rd2, rc)
    if (abs(rd2 - rd) > 0.01_dp*abs(rd)) error stop 'Delta'' eps convergence'

    ! Peaked current: Delta'(2/1) converges with the q grid (half the points).
    call peaked(1.0_dp, 3.0_dp)
    call cyl_current_from_q(nx, q, cur, dcur)
    call cyl_tearing_delta_prime(nx, q, dcur, 2, 1, 1.2_dp, 1e-4_dp, xs, rd, rc)
    if (rc /= CYL_OK) error stop 'peaked tearing status'
    call coarse_delta_prime(rd2)
    print '(a,2f12.5)', 'peaked nu=2 r_s Delta''(2/1), fine and coarse grid', rd, rd2
    if (abs(rd2 - rd) > 0.01_dp*abs(rd)) error stop 'Delta'' grid convergence'
    print '(a)', 'PASS cylindrical kink and tearing'
contains
    subroutine coarse_delta_prime(r)
        real(dp), intent(out) :: r
        integer, parameter :: nc = (nx + 1)/2
        real(dp) :: qc(nc), cc(nc), dc(nc), xc
        integer :: rcc
        qc = q(1:nx:2)
        call cyl_current_from_q(nc, qc, cc, dc)
        call cyl_tearing_delta_prime(nc, qc, dc, 2, 1, 1.2_dp, 1e-4_dp, xc, r, rcc)
        if (rcc /= CYL_OK) error stop 'coarse status'
    end subroutine
    subroutine peaked(qaxis, qedge)
        real(dp), intent(in) :: qaxis, qedge
        real(dp) :: nu, xx
        integer :: k
        nu = qedge/qaxis - 1
        do k = 1, nx
            xx = real(k-1, dp)/real(nx-1, dp)
            if (k == 1) then
                q(k) = qaxis
            else
                q(k) = qedge*xx*xx/(1 - (1 - xx*xx)**(nu + 1))
            end if
        end do
    end subroutine
    real(dp) function channel(xx, xc, width) result(ie)
        real(dp), intent(in) :: xx, xc, width
        real(dp) :: s, t
        integer :: k
        ! Enclosed current of J = (2/q0) H(xc - x) smoothed, normalised so
        ! that q = q0 inside: ie = int_0^x 2 s H ds / 1 (trapezoid, fine grid).
        ie = 0
        do k = 1, 400
            s = xx*(real(k, dp) - 0.5_dp)/400
            t = 0.5_dp*(1 - tanh((s - xc)/width))
            ie = ie + 2*s*t*xx/400
        end do
    end function
end program
