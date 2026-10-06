!> Linear Vlasov (collisionless continuum) response of a Maxwellian electron
!> OCP with a fixed neutralizing background: the oracles of the microscopic
!> plasma oscillation (C3). Plasma units omega_p = v_th = lambda_D = 1.
!>
!> - Plasma dispersion function Z(zeta) = i sqrt(pi) exp(-zeta^2) - 2 D(zeta)
!>   with the entire Dawson series D(z) = sum_n (-2)^n n!/(2n+1)! 2^n ... in
!>   quadruple precision (valid on the whole plane, analytic continuation
!>   included; |zeta| <= 6 keeps the series cancellation below 1e-18).
!> - Longitudinal dielectric epsilon(k, omega) = 1 + (1 + zeta Z(zeta))/k^2,
!>   zeta = omega/(sqrt(2) k), and its least-damped Landau root by Newton.
!> - Initial-value response of a density perturbation with Maxwellian
!>   velocities, the Volterra equation
!>   N(t) = exp(-k^2 t^2/2) - int_0^t (t-s) exp(-k^2 (t-s)^2/2) N(s) ds,
!>   solved by the trapezoidal product rule with Richardson extrapolation
!>   (h, h/2). This is the full finite-k response including the ballistic
!>   transient; its late time is the Landau root.
module kin6d_vlasov_response
    use, intrinsic :: iso_fortran_env, only: dp => real64, qp => real128
    implicit none
    private
    public :: plasma_dispersion, landau_root, vlasov_density_response
    public :: fit_damped_cosine

    real(qp), parameter :: pi_q = 3.14159265358979323846264338327950288_qp
contains
    !> Z(zeta) and Z'(zeta) = -2 (1 + zeta Z).
    subroutine plasma_dispersion(zeta, z, dz)
        complex(dp), intent(in) :: zeta
        complex(dp), intent(out) :: z, dz
        complex(qp) :: x, term, total, x2
        integer :: n
        x = cmplx(real(zeta, qp), real(aimag(zeta), qp), qp)
        x2 = x*x
        ! D(x) = sum_n (-1)^n 2^n x^(2n+1)/(2n+1)!!
        term = x; total = x
        do n = 1, 400
            term = -term*2*x2/(2*n + 1)
            total = total + term
            if (abs(term) <= 1.0e-30_qp*abs(total)) exit
        end do
        total = cmplx(0.0_qp, sqrt(pi_q), qp)*exp(-x2) - 2*total
        z = cmplx(real(total, dp), real(aimag(total), dp), dp)
        dz = -2*(1 + zeta*z)
    end subroutine

    !> Least-damped root omega of epsilon(k, omega) = 0 (Newton from the
    !> Bohm-Gross start); ok if converged to 1e-13.
    subroutine landau_root(k, omega, ok)
        real(dp), intent(in) :: k
        complex(dp), intent(out) :: omega
        logical, intent(out) :: ok
        complex(dp) :: zeta, z, dz, eps, deps
        integer :: it
        omega = cmplx(sqrt(1 + 3*k*k), -0.01_dp, dp)
        ok = .false.
        do it = 1, 100
            zeta = omega/(sqrt(2.0_dp)*k)
            call plasma_dispersion(zeta, z, dz)
            eps = 1 + (1 + zeta*z)/k**2
            deps = (z + zeta*dz)/k**2/(sqrt(2.0_dp)*k)
            omega = omega - eps/deps
            if (abs(eps/deps) < 1.0e-13_dp*abs(omega)) then
                ok = .true.; return
            end if
        end do
    end subroutine

    !> N(t_j) on t_j = (j-1) h, j = 1..size(n): normalized density response.
    subroutine vlasov_density_response(k, h, n)
        real(dp), intent(in) :: k, h
        real(dp), intent(out) :: n(:)
        real(dp), allocatable :: coarse(:), fine(:)
        integer :: m, j
        m = size(n)
        allocate (coarse(m), fine(2*m - 1))
        call trapezoid(k, h, coarse)
        call trapezoid(k, h/2, fine)
        do j = 1, m
            n(j) = (4*fine(2*j - 1) - coarse(j))/3
        end do
    end subroutine

    subroutine trapezoid(k, h, n)
        real(dp), intent(in) :: k, h
        real(dp), intent(out) :: n(:)
        real(dp), allocatable :: kernel(:)
        real(dp) :: s
        integer :: i, j, m
        m = size(n)
        allocate (kernel(0:m - 1))
        do i = 0, m - 1
            kernel(i) = i*h*exp(-(k*i*h)**2/2)
        end do
        n(1) = 1
        do i = 2, m
            s = kernel(i - 1)*n(1)/2
            do j = 2, i - 1
                s = s + kernel(i - j)*n(j)
            end do
            ! kernel(0) = 0: the scheme is explicit.
            n(i) = exp(-(k*(i - 1)*h)**2/2) - h*s
        end do
    end subroutine

    !> Least-squares fit y(t) ~ exp(-gamma t) (a cos(omega t) + b sin(omega t))
    !> over the samples: omega on [omin, omax] by a grid then golden section,
    !> gamma by golden section on [-0.2, 1] for each omega.
    subroutine fit_damped_cosine(t, y, omin, omax, omega, gamma, amplitude)
        real(dp), intent(in) :: t(:), y(:), omin, omax
        real(dp), intent(out) :: omega, gamma, amplitude
        real(dp) :: best, r, w, lo, hi, a, b, g
        integer :: i
        best = huge(1.0_dp); omega = omin
        do i = 0, 200
            w = omin + (omax - omin)*i/200
            call inner(w, g, r)
            if (r < best) then
                best = r; omega = w
            end if
        end do
        lo = omega - (omax - omin)/200; hi = omega + (omax - omin)/200
        do i = 1, 80
            a = hi - 0.618033988749895_dp*(hi - lo); b = lo + 0.618033988749895_dp*(hi - lo)
            call inner(a, g, r)
            call inner(b, g, best)
            if (r < best) then
                hi = b
            else
                lo = a
            end if
        end do
        omega = (lo + hi)/2
        call inner(omega, gamma, r, amplitude)
    contains
        subroutine inner(w, gbest, rbest, amp)
            real(dp), intent(in) :: w
            real(dp), intent(out) :: gbest, rbest
            real(dp), intent(out), optional :: amp
            real(dp) :: glo, ghi, ga, gb, ra, rb
            integer :: it
            glo = -0.2_dp; ghi = 1.0_dp
            do it = 1, 80
                ga = ghi - 0.618033988749895_dp*(ghi - glo)
                gb = glo + 0.618033988749895_dp*(ghi - glo)
                ra = residual(w, ga); rb = residual(w, gb)
                if (ra < rb) then
                    ghi = gb
                else
                    glo = ga
                end if
            end do
            gbest = (glo + ghi)/2
            rbest = residual(w, gbest, amp)
        end subroutine

        real(dp) function residual(w, g, amp) result(res)
            real(dp), intent(in) :: w, g
            real(dp), intent(out), optional :: amp
            real(dp) :: c(size(t)), s(size(t)), m11, m12, m22, r1, r2, det, ca, cb
            c = exp(-g*t)*cos(w*t); s = exp(-g*t)*sin(w*t)
            m11 = sum(c*c); m12 = sum(c*s); m22 = sum(s*s)
            r1 = sum(c*y); r2 = sum(s*y)
            det = m11*m22 - m12*m12
            ca = (m22*r1 - m12*r2)/det; cb = (m11*r2 - m12*r1)/det
            res = sum((y - ca*c - cb*s)**2)
            if (present(amp)) amp = sqrt(ca*ca + cb*cb)
        end function
    end subroutine
end module kin6d_vlasov_response
