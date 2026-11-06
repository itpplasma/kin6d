!> Monoenergetic NTV flux kernel of one toroidal harmonic n with Lorentz
!! pitch-angle scattering (#183), the simpler reference discretisation of
!! drift-kinetics.md: Legendre polynomials P_l(xi) and Fourier collocation in
!! theta, as for the axisymmetric problem (dk_mono).
!!
!! With f = F(theta, xi) exp(-i n zeta), the local equation per unit speed is
!!   xi k (d_theta - i n q) F - (1 - xi^2)/2 k d_theta(ln b) d_xi F
!!     + i (n/v) Omega F - (nu/2) d_xi((1 - xi^2) d_xi F) = lambda S/v,
!! b = |B|/B0, lambda = (1 - xi^2)/b, nu = nu_D/v,
!! Omega = omega_E + pre lambda (cos(theta) - offset) (grad-B precession,
!! bounce average -(mu B0/(e chi')) eps' (2E/K - 1 - offset)) and
!! S = i n sum_m c_m exp(i m theta) for B1 = Re sum_m c_m exp(i(m theta - n zeta)).
!! Its bounce average over a trapped orbit is the equation of the stage-1
!! bounce reference, i n Omega I_th f - 2 nu_D d_lambda(lambda I_xi d_lambda f)
!! = lambda S_hat, and the flux kernel is that reference's
!!   R = int dlambda oint dtheta/|xi| lambda Re(S^* f)
!!     = int dtheta int dxi (2/b) lambda Re(S^* F),
!! here over trapped and passing particles. (1 - xi^2) couples l to l +- 2,
!! so the rows l = 2j, 2j + 1 form pair blocks of size 2 n_theta, block
!! tridiagonal in j, factorised once by block elimination from the highest j.
!! Tangents with respect to nu and omega_E reuse the factorisation.
module kin6d_dk_ntv
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_dk_ntv_leaves, only: dk_ntv_point
    implicit none
    private
    public :: ntv_surface_t, ntv_kernel_t, ntv_kernel
    integer, parameter, public :: NTV_OK = 0, NTV_BAD_INPUT = 1, NTV_SINGULAR = 2

    !> Surface of the frozen-field benchmark and the perturbation spectrum.
    type :: ntv_surface_t
        real(dp) :: eps = 0, k = 0, q = 0, offset = 0
        integer :: n = 0
        integer, allocatable :: m(:)
        complex(dp), allocatable :: c(:)
    end type

    !> Flux kernel, its tangents and its two-resolution indicator.
    type :: ntv_kernel_t
        real(dp) :: r = 0, dr_dnu = 0, dr_domega = 0
        ! Directional two-resolution indicators: theta points 2 n + 1, and
        ! Legendre polynomials 2 L; their sum indicates the error of r.
        real(dp) :: indicator_theta = 0, indicator_leg = 0
    end type

    type :: ntv_system_t
        integer :: n = 0, nblock = 0
        real(dp) :: nu = 0, nq = 0
        complex(dp) :: shift = 0, omega = 0, precess = 0
        real(dp), allocatable :: d(:, :), stream(:), mirror(:), drift(:), source(:), &
            measure(:)
        complex(dp), allocatable :: s(:), s_lu(:, :, :), y(:, :, :)
        integer, allocatable :: piv(:, :)
    end type

    real(dp), parameter :: pi = 4*atan(1.0_dp)
    complex(dp), parameter :: iu = (0.0_dp, 1.0_dp)

    interface
        subroutine zgetrf(m, n, a, lda, ipiv, info)
            import :: dp
            integer, intent(in) :: m, n, lda
            complex(dp), intent(inout) :: a(lda, *)
            integer, intent(out) :: ipiv(*), info
        end subroutine zgetrf
        subroutine zgetrs(trans, n, nrhs, a, lda, ipiv, b, ldb, info)
            import :: dp
            character, intent(in) :: trans
            integer, intent(in) :: n, nrhs, lda, ldb, ipiv(*)
            complex(dp), intent(in) :: a(lda, *)
            complex(dp), intent(inout) :: b(ldb, *)
            integer, intent(out) :: info
        end subroutine zgetrs
    end interface

contains

    !> Kernel R at speed v [m/s], deflection frequency nu_d [1/s], omega_E
    !> [rad/s] and precession prefactor pre [rad/s], on n_theta (odd) points
    !> and n_leg (even) Legendre polynomials, with the directional
    !> indicators at 2 n_theta + 1 points and at 2 n_leg polynomials.
    subroutine ntv_kernel(sf, v, nu_d, omega_e, pre, n_theta, n_leg, kern, status)
        type(ntv_surface_t), intent(in) :: sf
        real(dp), intent(in) :: v, nu_d, omega_e, pre
        integer, intent(in) :: n_theta, n_leg
        type(ntv_kernel_t), intent(out) :: kern
        integer, intent(out) :: status
        type(ntv_system_t) :: a
        complex(dp), allocatable :: x(:, :), dx(:, :), rhs(:, :)
        real(dp) :: fine

        call factor(sf, v, nu_d, omega_e, pre, n_theta, n_leg, a, status)
        if (status /= NTV_OK) return
        allocate (rhs(2*a%n, 0:a%nblock - 1))
        rhs = 0
        rhs(1:a%n, 0) = a%source*a%s/v
        rhs(a%n + 1:, 0) = 0
        if (a%nblock > 1) rhs(1:a%n, 1) = -a%source*a%s/v
        call solve(a, rhs, x, status)
        if (status /= NTV_OK) return
        kern%r = flux(a, x)

        ! A dx = -dA x with dA = diag(l(l + 1)/2) for nu and i n/v for omega_E.
        rhs = -lorentz(a, x)
        call solve(a, rhs, dx, status)
        if (status /= NTV_OK) return
        kern%dr_dnu = flux(a, dx)/v
        rhs = -a%shift/v*x
        call solve(a, rhs, dx, status)
        if (status /= NTV_OK) return
        kern%dr_domega = flux(a, dx)

        call refined(sf, v, nu_d, omega_e, pre, 2*n_theta + 1, n_leg, fine, status)
        if (status /= NTV_OK) return
        kern%indicator_theta = fine - kern%r
        call refined(sf, v, nu_d, omega_e, pre, n_theta, 2*n_leg, fine, status)
        if (status /= NTV_OK) return
        kern%indicator_leg = fine - kern%r
    end subroutine ntv_kernel

    subroutine refined(sf, v, nu_d, omega_e, pre, n_theta, n_leg, r, status)
        type(ntv_surface_t), intent(in) :: sf
        real(dp), intent(in) :: v, nu_d, omega_e, pre
        integer, intent(in) :: n_theta, n_leg
        real(dp), intent(out) :: r
        integer, intent(out) :: status
        type(ntv_system_t) :: a
        complex(dp), allocatable :: x(:, :), rhs(:, :)
        call factor(sf, v, nu_d, omega_e, pre, n_theta, n_leg, a, status)
        if (status /= NTV_OK) return
        allocate (rhs(2*a%n, 0:a%nblock - 1))
        rhs = 0
        rhs(1:a%n, 0) = a%source*a%s/v
        rhs(1:a%n, 1) = -a%source*a%s/v
        call solve(a, rhs, x, status)
        if (status /= NTV_OK) return
        r = flux(a, x)
    end subroutine refined

    subroutine factor(sf, v, nu_d, omega_e, pre, n_theta, n_leg, a, status)
        type(ntv_surface_t), intent(in) :: sf
        real(dp), intent(in) :: v, nu_d, omega_e, pre
        integer, intent(in) :: n_theta, n_leg
        type(ntv_system_t), intent(out) :: a
        integer, intent(out) :: status
        real(dp), allocatable :: theta(:), b(:)
        complex(dp), allocatable :: m(:, :)
        integer :: n, j, i, info

        status = NTV_BAD_INPUT
        if (mod(n_theta, 2) == 0 .or. n_theta < 3 .or. mod(n_leg, 2) /= 0 .or. &
            n_leg < 4 .or. v <= 0 .or. nu_d <= 0 .or. sf%n == 0 .or. &
            .not. allocated(sf%m)) return
        n = n_theta
        a%n = n
        a%nblock = n_leg/2
        a%nu = nu_d/v
        a%shift = iu*sf%n/v
        a%nq = sf%n*sf%q
        a%omega = omega_e
        a%precess = pre
        allocate (theta(n), b(n), a%d(n, n), a%stream(n), a%mirror(n), a%drift(n), &
            a%source(n), a%measure(n), a%s(n))
        theta = [(2*pi*(i - 1)/n, i=1, n)]
        call fourier_derivative(n, a%d)
        ! Streaming along the field line in the frame exp(-i n zeta): d_theta - i n q.
        a%stream = sf%k
        call dk_ntv_point(cos(theta), sin(theta), sf%eps, sf%k, sf%offset, b, &
            a%mirror, a%drift, a%source, a%measure)
        a%s = 0
        do j = 1, size(sf%m)
            a%s = a%s + iu*sf%n*sf%c(j)*exp(iu*sf%m(j)*theta)
        end do

        ! S_J = Di_J; S_j = Di_j - Up_j Y_{j+1}, Y_j = S_j^-1 Lo_j.
        allocate (a%s_lu(2*n, 2*n, 0:a%nblock - 1), a%y(2*n, 2*n, a%nblock - 1), &
            a%piv(2*n, 0:a%nblock - 1), m(2*n, 2*n))
        status = NTV_SINGULAR
        do j = a%nblock - 1, 0, -1
            m = block(a, j, j)
            if (j < a%nblock - 1) m = m - matmul(block(a, j, j + 1), a%y(:, :, j + 1))
            call zgetrf(2*n, 2*n, m, 2*n, a%piv(:, j), info)
            if (info /= 0) return
            a%s_lu(:, :, j) = m
            if (j > 0) then
                a%y(:, :, j) = block(a, j, j - 1)
                call zgetrs('N', 2*n, 2*n, m, 2*n, a%piv(:, j), a%y(:, :, j), 2*n, info)
            end if
        end do
        status = NTV_OK
    end subroutine factor

    !> Solve A x = rhs, pair block j holding l = 2j (rows 1:n) and 2j + 1.
    subroutine solve(a, rhs, x, status)
        type(ntv_system_t), intent(in) :: a
        complex(dp), intent(in) :: rhs(:, 0:)
        complex(dp), allocatable, intent(out) :: x(:, :)
        integer, intent(out) :: status
        complex(dp), allocatable :: z(:, :)
        integer :: j, info, nn
        nn = 2*a%n
        allocate (z(nn, 0:a%nblock - 1), x(nn, 0:a%nblock - 1))
        status = NTV_SINGULAR
        do j = a%nblock - 1, 0, -1
            z(:, j) = rhs(:, j)
            if (j < a%nblock - 1) z(:, j) = z(:, j) - matmul(block(a, j, j + 1), z(:, j + 1))
            call zgetrs('N', nn, 1, a%s_lu(:, :, j), nn, a%piv(:, j), z(:, j), nn, info)
            if (info /= 0) return
        end do
        x(:, 0) = z(:, 0)
        do j = 1, a%nblock - 1
            x(:, j) = z(:, j) - matmul(a%y(:, :, j), x(:, j - 1))
        end do
        status = NTV_OK
    end subroutine solve

    !> Rows of pair block j, columns of pair block jc: the couplings of
    !> rows l to columns l' with |l - l'| <= 2.
    function block(a, j, jc) result(m)
        type(ntv_system_t), intent(in) :: a
        integer, intent(in) :: j, jc
        complex(dp) :: m(2*a%n, 2*a%n)
        integer :: r, c, l, lc
        m = 0
        do r = 0, 1
            do c = 0, 1
                l = 2*j + r
                lc = 2*jc + c
                m(r*a%n + 1:(r + 1)*a%n, c*a%n + 1:(c + 1)*a%n) = entry(a, l, lc)
            end do
        end do
    end function block

    !> Coefficient matrix of x_lc in row l.
    function entry(a, l, lc) result(m)
        type(ntv_system_t), intent(in) :: a
        integer, intent(in) :: l, lc
        complex(dp) :: m(a%n, a%n)
        integer :: i
        real(dp) :: f
        complex(dp) :: p(a%n)
        m = 0
        if (abs(l - lc) > 2) return
        p = a%shift*a%precess*a%drift
        select case (lc - l)
        case (-1)
            f = l/real(2*l - 1, dp)
            do i = 1, a%n
                m(i, :) = a%stream(i)*f*a%d(i, :)
                m(i, i) = m(i, i) - a%stream(i)*f*iu*a%nq &
                    + a%mirror(i)*(l - 1)*l/real(2*l - 1, dp)
            end do
        case (1)
            f = (l + 1)/real(2*l + 3, dp)
            do i = 1, a%n
                m(i, :) = a%stream(i)*f*a%d(i, :)
                m(i, i) = m(i, i) - a%stream(i)*f*iu*a%nq &
                    - a%mirror(i)*(l + 1)*(l + 2)/real(2*l + 3, dp)
            end do
        case (0)
            do i = 1, a%n
                m(i, i) = a%nu*l*(l + 1)/2.0_dp + a%shift*a%omega + p(i)*(1 - c0(l))
            end do
        case (-2)
            do i = 1, a%n
                m(i, i) = -p(i)*cp(lc)
            end do
        case (2)
            do i = 1, a%n
                m(i, i) = -p(i)*cm(lc)
            end do
        end select
    end function entry

    ! xi^2 P_l = cp P_{l+2} + c0 P_l + cm P_{l-2}, so (1 - xi^2) P_l has
    ! -cp, 1 - c0 and -cm.
    pure real(dp) function cp(l)
        integer, intent(in) :: l
        cp = (l + 1)*(l + 2)/real((2*l + 1)*(2*l + 3), dp)
    end function cp

    pure real(dp) function c0(l)
        integer, intent(in) :: l
        c0 = (l + 1)**2/real((2*l + 1)*(2*l + 3), dp)
        if (l > 0) c0 = c0 + l**2/real((2*l + 1)*(2*l - 1), dp)
    end function c0

    pure real(dp) function cm(l)
        integer, intent(in) :: l
        cm = 0
        if (l > 1) cm = l*(l - 1)/real((2*l + 1)*(2*l - 1), dp)
    end function cm

    function lorentz(a, x) result(y)
        type(ntv_system_t), intent(in) :: a
        complex(dp), intent(in) :: x(:, 0:)
        complex(dp) :: y(size(x, 1), 0:ubound(x, 2))
        integer :: j, l
        do j = 0, ubound(x, 2)
            l = 2*j
            y(1:a%n, j) = l*(l + 1)/2.0_dp*x(1:a%n, j)
            y(a%n + 1:, j) = (l + 1)*(l + 2)/2.0_dp*x(a%n + 1:, j)
        end do
    end function lorentz

    !> R = (2 pi/n) sum (2/b^2) Re(S^* ((4/3) a_0 - (4/15) a_2)).
    real(dp) function flux(a, x)
        type(ntv_system_t), intent(in) :: a
        complex(dp), intent(in) :: x(:, 0:)
        complex(dp) :: moment(a%n)
        moment = 4*x(1:a%n, 0)/3
        if (ubound(x, 2) >= 1) moment = moment - 4*x(1:a%n, 1)/15
        flux = 2*pi/a%n*sum(a%measure*real(conjg(a%s)*moment, dp))
    end function flux

    subroutine fourier_derivative(n, d)
        integer, intent(in) :: n
        real(dp), intent(out) :: d(n, n)
        integer :: i, j
        do j = 1, n
            do i = 1, n
                if (i == j) then
                    d(i, j) = 0
                else
                    d(i, j) = 0.5_dp*(-1)**(i - j)/sin(pi*(i - j)/n)
                end if
            end do
        end do
    end subroutine fourier_derivative

end module kin6d_dk_ntv
