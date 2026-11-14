!> Low-beta, large-aspect-ratio cylindrical stability of a tokamak q profile.
!>
!> Two classical calculations that need only q(x), x = r/a:
!> - ideal external kink: Newcomb shooting of the Euler-Lagrange equation of
!>   the cylindrical energy (Wesson, Tokamaks, 4th ed., sec. 6.5) with a
!>   conducting wall at x_b = b/a;
!> - tearing stability index Delta' of a resonant surface q(x_s) = m/n from
!>   the cylindrical outer-region equation (Furth, Rutherford and Selberg 1973),
!>   with the current gradient written through q as
!>   mu0 R0 j/B0 = (1/x) d(x^2/q)/dx.
!> Shaping and toroidicity enter only through q, so results are teaching-level
!> reductions, not toroidal stability limits.
module milca_cylinder_stability
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    private

    integer, parameter, public :: CYL_OK = 0, CYL_INVALID = 1, &
        CYL_NO_RESONANCE = 2, CYL_RESONANT_INSIDE = 3

    public :: cyl_current_from_q, cyl_resonant_surface, cyl_kink_energy, &
        cyl_tearing_delta_prime, cyl_tearing_delta_prime_step

contains

    !> Normalised current J = mu0 R0 j/B0 = (1/x) d(x^2/q)/dx and dJ/dx on the
    !> uniform grid x(i) = (i-1)/(nx-1), by second-order differences.
    pure subroutine cyl_current_from_q(nx, q, cur, dcur)
        integer, intent(in) :: nx
        real(dp), intent(in) :: q(nx)
        real(dp), intent(out) :: cur(nx), dcur(nx)
        real(dp) :: h, g(nx), x
        integer :: i
        h = 1.0_dp/real(nx-1, dp)
        do i = 1, nx
            x = real(i-1, dp)*h
            g(i) = x*x/q(i)
        end do
        cur(1) = 2/q(1)
        do i = 2, nx-1
            x = real(i-1, dp)*h
            cur(i) = (g(i+1) - g(i-1))/(2*h*x)
        end do
        cur(nx) = (3*g(nx) - 4*g(nx-1) + g(nx-2))/(2*h)
        dcur(1) = 0
        do i = 2, nx-1
            dcur(i) = (cur(i+1) - cur(i-1))/(2*h)
        end do
        dcur(nx) = (3*cur(nx) - 4*cur(nx-1) + cur(nx-2))/(2*h)
    end subroutine

    !> Outermost x in (0,1] with q(x) = m/n, by linear interpolation.
    pure subroutine cyl_resonant_surface(nx, q, m, n, xs, found)
        integer, intent(in) :: nx, m, n
        real(dp), intent(in) :: q(nx)
        real(dp), intent(out) :: xs
        logical, intent(out) :: found
        real(dp) :: qs, h, f0, f1
        integer :: i
        qs = real(m, dp)/real(n, dp)
        h = 1.0_dp/real(nx-1, dp)
        found = .false.
        xs = 0
        do i = nx-1, 1, -1
            f0 = q(i) - qs
            f1 = q(i+1) - qs
            if (f0*f1 < 0 .or. abs(f0) <= tiny(f0)) then
                xs = (real(i-1, dp) + f0/(f0 - f1))*h
                found = .true.
                return
            end if
        end do
    end subroutine

    pure real(dp) function lin(nx, f, x) result(v)
        integer, intent(in) :: nx
        real(dp), intent(in) :: f(nx), x
        real(dp) :: s, w
        integer :: i
        s = min(max(x, 0.0_dp), 1.0_dp)*real(nx-1, dp)
        i = min(int(s) + 1, nx-1)
        w = s - real(i-1, dp)
        v = (1-w)*f(i) + w*f(i+1)
    end function

    !> Ideal external kink energy for poloidal mode m and toroidal mode n, with
    !> xi(a) = 1 and the common factor pi^2 B0^2 a^2/(mu0 R0) removed:
    !>   dW = F_a^2 a xi'/xi + 2 F_a/q_a + (1 + m lambda) F_a^2,
    !>   F = n/m - 1/q, lambda = (1 + x_b^-2m)/(1 - x_b^-2m),
    !> where xi solves (x^3 F^2 xi')' = (m^2-1) x F^2 xi from the axis.
    !> gamma_hat = sqrt(-dW/(2K)) estimates gamma*tau_A, tau_A = R0 sqrt(mu0 rho)/B0,
    !> with K = int x (xi^2 + ((x xi)')^2/m^2) dx the incompressible inertia.
    !> A resonance q = m/n inside the plasma makes this mode ideally stable in
    !> the cylinder (status CYL_RESONANT_INSIDE): it is then a tearing problem.
    pure subroutine cyl_kink_energy(nx, q, m, n, wall_ratio, dw, gamma_hat, status)
        integer, intent(in) :: nx, m, n
        real(dp), intent(in) :: q(nx), wall_ratio
        real(dp), intent(out) :: dw, gamma_hat
        integer, intent(out) :: status
        real(dp) :: qs, x, h, y(2), k1(2), k2(2), k3(2), k4(2), fa, lam, c
        real(dp) :: kin, xi, dxi, prev_k, cur_k
        integer :: i, nstep
        logical :: found
        dw = 0
        gamma_hat = 0
        status = CYL_INVALID
        if (nx < 8 .or. m < 1 .or. n < 1 .or. wall_ratio <= 1) return
        if (any(q <= 0)) return
        qs = real(m, dp)/real(n, dp)
        call cyl_resonant_surface(nx, q, m, n, x, found)
        if (found .and. x < 1) then
            status = CYL_RESONANT_INSIDE
            dw = huge(1.0_dp)
            return
        end if
        nstep = 4*(nx-1)
        h = (1.0_dp - 1e-3_dp)/real(nstep, dp)
        x = 1e-3_dp
        y(1) = x**(m-1)
        y(2) = x**3*fsq(x)*real(m-1, dp)*x**max(m-2, 0)
        if (m == 1) y(2) = 0
        kin = 0
        prev_k = inertia(x, y)
        do i = 1, nstep
            call el_rhs(x, y, k1)
            call el_rhs(x + h/2, y + h*k1/2, k2)
            call el_rhs(x + h/2, y + h*k2/2, k3)
            call el_rhs(x + h, y + h*k3, k4)
            y = y + (h/6)*(k1 + 2*k2 + 2*k3 + k4)
            x = x + h
            cur_k = inertia(x, y)
            kin = kin + h*(prev_k + cur_k)/2
            prev_k = cur_k
        end do
        xi = y(1)
        dxi = y(2)/fsq(1.0_dp)
        fa = real(n, dp)/real(m, dp) - 1/q(nx)
        c = (1/wall_ratio)**(2*m)
        lam = (1 + c)/(1 - c)
        dw = fa*fa*dxi/xi + 2*fa/q(nx) + (1 + m*lam)*fa*fa
        kin = kin/(xi*xi)
        if (dw < 0 .and. kin > 0) gamma_hat = sqrt(-dw/(2*kin))
        status = CYL_OK
    contains
        pure real(dp) function fsq(xx)
            real(dp), intent(in) :: xx
            real(dp) :: f
            f = real(n, dp)/real(m, dp) - 1/lin(nx, q, xx)
            fsq = f*f
        end function
        pure subroutine el_rhs(xx, yy, dy)
            real(dp), intent(in) :: xx, yy(2)
            real(dp), intent(out) :: dy(2)
            real(dp) :: f2
            f2 = fsq(xx)
            dy(1) = yy(2)/(xx**3*f2)
            dy(2) = real(m*m - 1, dp)*xx*f2*yy(1)
        end subroutine
        pure real(dp) function inertia(xx, yy)
            real(dp), intent(in) :: xx, yy(2)
            real(dp) :: d, dxx
            d = yy(2)/(xx**3*fsq(xx))
            dxx = yy(1) + xx*d
            inertia = xx*(yy(1)**2 + dxx*dxx/real(m*m, dp))
        end function
    end subroutine

    !> Tearing index r_s Delta' of the m/n surface for the q profile on the
    !> uniform grid, with a vacuum gap and an ideal wall at x_b = b/a. The
    !> outer equation psi'' + psi'/x - m^2 psi/x^2 - G psi = 0,
    !> G = q J'/(x (1 - n q/m)), is integrated to x_s -+ eps from both sides.
    pure subroutine cyl_tearing_delta_prime(nx, q, dcur, m, n, wall_ratio, &
            eps, xs, rs_delta_prime, status)
        integer, intent(in) :: nx, m, n
        real(dp), intent(in) :: q(nx), dcur(nx), wall_ratio, eps
        real(dp), intent(out) :: xs, rs_delta_prime
        integer, intent(out) :: status
        real(dp) :: y(2), c, din, dout
        logical :: found
        xs = 0
        rs_delta_prime = 0
        status = CYL_INVALID
        if (nx < 8 .or. m < 1 .or. n < 1 .or. wall_ratio <= 1 .or. eps <= 0) return
        if (any(q <= 0)) return
        call cyl_resonant_surface(nx, q, m, n, xs, found)
        if (.not. found .or. xs >= 1 - 2*eps .or. xs <= 2*eps) then
            status = CYL_NO_RESONANCE
            return
        end if
        y(1) = 1e-3_dp**m
        y(2) = real(m, dp)*1e-3_dp**(m-1)
        call shoot(1e-3_dp, xs - eps, y)
        din = y(2)/y(1)
        c = wall_ratio**(2*m)
        y(1) = 1
        y(2) = real(m, dp)*(1 + c)/(1 - c)
        call shoot(1.0_dp, xs + eps, y)
        dout = y(2)/y(1)
        rs_delta_prime = xs*(dout - din)
        status = CYL_OK
    contains
        pure subroutine shoot(x0, x1, yy)
            real(dp), intent(in) :: x0, x1
            real(dp), intent(inout) :: yy(2)
            real(dp) :: x, h, hmax, dir, k1(2), k2(2), k3(2), k4(2), s
            integer :: it
            x = x0
            dir = sign(1.0_dp, x1 - x0)
            hmax = 0.25_dp/real(nx-1, dp)
            do it = 1, 200000
                if (dir*(x1 - x) <= 0) exit
                h = min(hmax, 0.1_dp*abs(x - xs), abs(x1 - x))
                h = dir*h
                call rhs(x, yy, k1)
                call rhs(x + h/2, yy + h*k1/2, k2)
                call rhs(x + h/2, yy + h*k2/2, k3)
                call rhs(x + h, yy + h*k3, k4)
                yy = yy + (h/6)*(k1 + 2*k2 + 2*k3 + k4)
                x = x + h
                s = abs(yy(1)) + abs(yy(2))
                if (s > 1e100_dp) yy = yy/s
            end do
        end subroutine
        pure subroutine rhs(x, yy, dy)
            real(dp), intent(in) :: x, yy(2)
            real(dp), intent(out) :: dy(2)
            real(dp) :: qx, g
            qx = lin(nx, q, x)
            g = qx*lin(nx, dcur, x)/(x*(1 - real(n, dp)*qx/real(m, dp)))
            dy(1) = yy(2)
            dy(2) = -yy(2)/x + real(m*m, dp)*yy(1)/(x*x) + g*yy(1)
        end subroutine
    end subroutine

    !> Closed form of r_s Delta' for a uniform current channel x < x1 inside a
    !> current-free region with the resonance at x_s in (x1, 1): the
    !> independent reference of the shooting solver (Wesson sec. 6.8).
    pure real(dp) function cyl_tearing_delta_prime_step(m, n, q0, x1, xs, &
            wall_ratio) result(rd)
        integer, intent(in) :: m, n
        real(dp), intent(in) :: q0, x1, xs, wall_ratio
        real(dp) :: k, a, b, c, pin, dpin, pout, dpout
        k = -2/(1 - real(n, dp)*q0/real(m, dp))
        a = 1 + k/(2*m)
        b = -k*x1**(2*m)/(2*m)
        pin = a*xs**m + b*xs**(-m)
        dpin = m*(a*xs**(m-1) - b*xs**(-m-1))
        c = wall_ratio**(2*m)
        pout = xs**m - c*xs**(-m)
        dpout = m*(xs**(m-1) + c*xs**(-m-1))
        rd = xs*(dpout/pout - dpin/pin)
    end function

end module milca_cylinder_stability
