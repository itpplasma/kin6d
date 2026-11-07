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
!! Tangents with respect to nu and omega_E, the perturbation amplitudes c_m
!! and a direction in (eps, k, offset, c_J) reuse the factorisation: with
!! A F = H, R = <H, F>_W and the FortSym JVP of the pointwise leaves,
!! dR = <dH, F>_W + <H, F>_dW + <H, A^-1 (dH - dA F)>_W.
!!
!! Three model switches, each with a measurable defect, select the equation of
!! NEO-2-QL and NEO-RT instead of the stage-1 bounce reference:
!!   weight  NTV_WEIGHT_MU (default): drift and drive carry v_perp^2/(2 v^2 b)
!!           = lambda, the grad-B part of the reference; NTV_WEIGHT_FULL: they
!!           carry (v_par^2 + v_perp^2/2)/(v^2 b) = (1 + xi^2)/b, grad-B plus
!!           curvature at low beta, equal to lambda at xi = 0;
!!   current c_J  adds the parallel-current part of the curvature drift,
!!           pre c_J xi^2, c_J = 2 B0 (I' + q G')/((G + iota I) q eps') with the
!!           covariant Boozer components G, I and ' = d/ds (zero by default);
!!   boozer_stream  b.grad theta = k b (Boozer angles, Jacobian ~ 1/B^2)
!!           instead of the constant k, with the surface measure
!!           dtheta/(b.grad theta B) ~ 1/b^2 in the flux.
!!
!! A general field line enters through an ntv_line_t provider: its sample
!! returns, at the theta points of the solver, the coefficients the circular
!! leaf returns otherwise (b.grad theta, the mirror factor
!! b.grad theta d_theta(ln b)/2, the drift multiplying (1 + wsign xi^2), the
!! pure xi^2 drift, the source factor 2/(3 b) and the flux measure), with b
!! = |B|/B_ref. The perturbation spectrum (any m, n) and q stay in the
!! surface. Tangents along (eps, k, offset, c_J) exist only for the circular
!! surface.
module kin6d_dk_ntv
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_dk_ntv_leaves, only: dk_ntv_point, dk_ntv_point_jvp
    implicit none
    private
    public :: ntv_surface_t, ntv_kernel_t, ntv_kernel, ntv_line_t, ntv_circular_line_t
    integer, parameter, public :: NTV_OK = 0, NTV_BAD_INPUT = 1, NTV_SINGULAR = 2
    integer, parameter, public :: NTV_WEIGHT_MU = 0, NTV_WEIGHT_FULL = 1

    !> Field-line coefficients of a general surface (see the module header).
    type, abstract :: ntv_line_t
    contains
        procedure(line_sample), deferred :: sample
    end type

    abstract interface
        subroutine line_sample(self, theta, stream, mirror, drift, drift2, source, &
                measure)
            import :: ntv_line_t, dp
            class(ntv_line_t), intent(in) :: self
            real(dp), intent(in) :: theta(:)
            real(dp), intent(out) :: stream(:), mirror(:), drift(:), drift2(:), &
                source(:), measure(:)
        end subroutine line_sample
    end interface

    !> The circular stage-1 field line as a provider: the same leaf and
    !> arguments as a surface without line.
    type, extends(ntv_line_t) :: ntv_circular_line_t
        real(dp) :: eps = 0, k = 0, offset = 0, current = 0
        logical :: boozer_stream = .false.
    contains
        procedure :: sample => circular_sample
    end type

    !> Surface of the frozen-field benchmark and the perturbation spectrum;
    !> with line allocated, eps, k, offset, current and boozer_stream are
    !> replaced by the provider.
    type :: ntv_surface_t
        real(dp) :: eps = 0, k = 0, q = 0, offset = 0
        integer :: n = 0
        integer :: weight = NTV_WEIGHT_MU
        real(dp) :: current = 0
        logical :: boozer_stream = .false.
        class(ntv_line_t), allocatable :: line
        integer, allocatable :: m(:)
        complex(dp), allocatable :: c(:)
    end type

    !> Flux kernel, its tangents and its two-resolution indicator. dr_dc(j)
    !> = dR/dRe(c_j) + i dR/dIm(c_j); dr_dsurface is the derivative along
    !> the (eps, k, offset, current) components of the optional direction.
    type :: ntv_kernel_t
        real(dp) :: r = 0, dr_dnu = 0, dr_domega = 0, dr_dsurface = 0
        complex(dp), allocatable :: dr_dc(:)
        ! Directional two-resolution indicators: theta points 2 n + 1, and
        ! Legendre polynomials 2 L; their sum indicates the error of r.
        real(dp) :: indicator_theta = 0, indicator_leg = 0
    end type

    type :: ntv_system_t
        integer :: n = 0, nblock = 0
        ! Legendre factors of the drift and drive weight in xi^2: (1 - xi^2)
        ! has sign -1, (1 + xi^2) sign +1.
        real(dp) :: wsign = -1
        real(dp) :: nu = 0, nq = 0
        complex(dp) :: shift = 0, omega = 0, precess = 0
        real(dp), allocatable :: d(:, :), stream(:), mirror(:), drift(:), drift2(:), &
            source(:), measure(:)
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
    subroutine ntv_kernel(sf, v, nu_d, omega_e, pre, n_theta, n_leg, kern, status, &
            indicators, direction)
        type(ntv_surface_t), intent(in) :: sf
        ! Tangent direction: its eps, k, offset and current components.
        type(ntv_surface_t), intent(in), optional :: direction
        real(dp), intent(in) :: v, nu_d, omega_e, pre
        integer, intent(in) :: n_theta, n_leg
        type(ntv_kernel_t), intent(out) :: kern
        integer, intent(out) :: status
        ! A resolution ladder supplies its own two-resolution estimate and
        ! can skip the refined solves, which dominate memory.
        logical, intent(in), optional :: indicators
        type(ntv_system_t) :: a, da
        complex(dp), allocatable :: x(:, :), dx(:, :), rhs(:, :), ds(:)
        real(dp) :: fine
        integer :: j, i

        if (present(direction)) then
            status = NTV_BAD_INPUT
            if (allocated(sf%line)) return
        end if
        call factor(sf, v, nu_d, omega_e, pre, n_theta, n_leg, a, status)
        if (status /= NTV_OK) return
        call drive(a, v, rhs)
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

        ! Amplitudes: s is linear in c_j, ds = i n exp(i m_j theta).
        allocate (kern%dr_dc(size(sf%m)), ds(a%n))
        da = coefficients(a)
        do j = 1, size(sf%m)
            ds = [(iu*sf%n*exp(iu*sf%m(j)*2*pi*(i - 1)/a%n), i=1, a%n)]
            da%s = ds
            call drive(da, v, rhs)
            call solve(a, rhs, dx, status)
            if (status /= NTV_OK) return
            kern%dr_dc(j) = cmplx(pair(a, ds, x) + flux(a, dx), &
                pair(a, iu*ds, x) + flux(a, iu*dx), dp)
        end do

        if (present(direction)) then
            call tangent_system(sf, direction, a, da)
            call drive(da, v, rhs)
            rhs = rhs - apply(da, x)
            call solve(a, rhs, dx, status)
            if (status /= NTV_OK) return
            kern%dr_dsurface = flux(a, dx) &
                + 2*pi/a%n*sum(da%measure*real(conjg(a%s)*moment(a, x), dp))
        end if
        if (present(indicators)) then
            if (.not. indicators) return
        end if

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
        call drive(a, v, rhs)
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
        real(dp), allocatable :: theta(:), b(:), stream_b(:), mirror_b(:), measure_b(:)
        complex(dp), allocatable :: m(:, :)
        integer :: n, j, i, info

        status = NTV_BAD_INPUT
        if (mod(n_theta, 2) == 0 .or. n_theta < 3 .or. mod(n_leg, 2) /= 0 .or. &
            n_leg < 4 .or. v <= 0 .or. nu_d <= 0 .or. sf%n == 0 .or. &
            .not. allocated(sf%m) .or. &
            (sf%weight /= NTV_WEIGHT_MU .and. sf%weight /= NTV_WEIGHT_FULL)) return
        n = n_theta
        a%n = n
        a%nblock = n_leg/2
        a%nu = nu_d/v
        a%shift = iu*sf%n/v
        a%nq = sf%n*sf%q
        a%omega = omega_e
        a%precess = pre
        allocate (theta(n), b(n), a%d(n, n), a%stream(n), a%mirror(n), a%drift(n), &
            a%drift2(n), a%source(n), a%measure(n), a%s(n))
        theta = [(2*pi*(i - 1)/n, i=1, n)]
        call fourier_derivative(n, a%d)
        if (allocated(sf%line)) then
            call sf%line%sample(theta, a%stream, a%mirror, a%drift, a%drift2, &
                a%source, a%measure)
        else
            ! Streaming along the field line in the frame exp(-i n zeta):
            ! d_theta - i n q.
            a%stream = sf%k
            allocate (stream_b(n), mirror_b(n), measure_b(n))
            call dk_ntv_point(cos(theta), sin(theta), sf%eps, sf%k, sf%offset, b, &
                a%mirror, a%drift, a%source, a%measure, stream_b, mirror_b, measure_b)
            if (sf%boozer_stream) then
                a%stream = stream_b
                a%mirror = mirror_b
                a%measure = measure_b
            end if
            a%drift2 = sf%current
        end if
        if (sf%weight == NTV_WEIGHT_FULL) a%wsign = 1
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
        complex(dp) :: p(a%n), pj(a%n)
        m = 0
        if (abs(l - lc) > 2) return
        p = a%shift*a%precess*a%drift
        pj = a%shift*a%precess*a%drift2
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
                m(i, i) = a%nu*l*(l + 1)/2.0_dp + a%shift*a%omega &
                    + p(i)*(1 + a%wsign*c0(l)) + pj(i)*c0(l)
            end do
        case (-2)
            do i = 1, a%n
                m(i, i) = (a%wsign*p(i) + pj(i))*cp(lc)
            end do
        case (2)
            do i = 1, a%n
                m(i, i) = (a%wsign*p(i) + pj(i))*cm(lc)
            end do
        end select
    end function entry

    ! xi^2 P_l = cp P_{l+2} + c0 P_l + cm P_{l-2}, so (1 + wsign xi^2) P_l
    ! has wsign cp, 1 + wsign c0 and wsign cm.
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

    subroutine circular_sample(self, theta, stream, mirror, drift, drift2, source, &
            measure)
        class(ntv_circular_line_t), intent(in) :: self
        real(dp), intent(in) :: theta(:)
        real(dp), intent(out) :: stream(:), mirror(:), drift(:), drift2(:), &
            source(:), measure(:)
        real(dp) :: b(size(theta)), stream_b(size(theta)), mirror_b(size(theta)), &
            measure_b(size(theta))
        stream = self%k
        call dk_ntv_point(cos(theta), sin(theta), self%eps, self%k, self%offset, b, &
            mirror, drift, source, measure, stream_b, mirror_b, measure_b)
        if (self%boozer_stream) then
            stream = stream_b
            mirror = mirror_b
            measure = measure_b
        end if
        drift2 = self%current
    end subroutine circular_sample

    !> <ds, x>_W: the flux functional with the conjugated drive replaced.
    real(dp) function pair(a, ds, x)
        type(ntv_system_t), intent(in) :: a
        complex(dp), intent(in) :: ds(:), x(:, 0:)
        pair = 2*pi/a%n*sum(a%measure*real(conjg(ds)*moment(a, x), dp))
    end function pair

    !> int (1 + wsign xi^2) F dxi per theta point.
    function moment(a, x) result(mom)
        type(ntv_system_t), intent(in) :: a
        complex(dp), intent(in) :: x(:, 0:)
        complex(dp) :: mom(a%n)
        mom = (2 + 2*a%wsign/3)*x(1:a%n, 0)
        if (ubound(x, 2) >= 1) mom = mom + 4*a%wsign*x(1:a%n, 1)/15
    end function moment

    !> Tangent of the operator and the pointwise leaves along a direction in
    !> (eps, k, offset, current): da holds dA (no collision or E x B part,
    !> which do not depend on these), dsource and dmeasure.
    subroutine tangent_system(sf, dir, a, da)
        type(ntv_surface_t), intent(in) :: sf, dir
        type(ntv_system_t), intent(in) :: a
        type(ntv_system_t), intent(out) :: da
        real(dp), allocatable :: theta(:), db(:), dstream_b(:), dmirror_b(:), &
            dmeasure_b(:)
        integer :: i, n
        n = a%n
        da = coefficients(a)
        da%nu = 0
        da%omega = 0
        da%drift2 = dir%current
        allocate (theta(n), db(n), dstream_b(n), dmirror_b(n), dmeasure_b(n))
        theta = [(2*pi*(i - 1)/n, i=1, n)]
        call dk_ntv_point_jvp(cos(theta), sin(theta), sf%eps, sf%k, sf%offset, &
            0*theta, 0*theta, dir%eps + 0*theta, dir%k + 0*theta, dir%offset + 0*theta, &
            db, da%mirror, da%drift, da%source, da%measure, dstream_b, dmirror_b, &
            dmeasure_b)
        da%stream = dir%k
        if (sf%boozer_stream) then
            da%stream = dstream_b
            da%mirror = dmirror_b
            da%measure = dmeasure_b
        end if
    end subroutine tangent_system

    !> The operator coefficients of a without its factors.
    function coefficients(a) result(c)
        type(ntv_system_t), intent(in) :: a
        type(ntv_system_t) :: c
        c%n = a%n; c%nblock = a%nblock; c%wsign = a%wsign; c%nu = a%nu; c%nq = a%nq
        c%drift2 = a%drift2; c%shift = a%shift; c%omega = a%omega
        c%precess = a%precess
        c%d = a%d; c%stream = a%stream; c%mirror = a%mirror; c%drift = a%drift
        c%source = a%source; c%measure = a%measure; c%s = a%s
    end function coefficients

    !> y = A x for the pair-block tridiagonal operator.
    function apply(a, x) result(y)
        type(ntv_system_t), intent(in) :: a
        complex(dp), intent(in) :: x(:, 0:)
        complex(dp) :: y(size(x, 1), 0:ubound(x, 2))
        integer :: j
        do j = 0, a%nblock - 1
            y(:, j) = matmul(block(a, j, j), x(:, j))
            if (j > 0) y(:, j) = y(:, j) + matmul(block(a, j, j - 1), x(:, j - 1))
            if (j < a%nblock - 1) y(:, j) = y(:, j) + matmul(block(a, j, j + 1), x(:, j + 1))
        end do
    end function apply

    !> R = (2 pi/n) sum (2/b^2) Re(S^* int (1 + wsign xi^2) F dxi), with
    !> int (1 + wsign xi^2) F dxi = (2 + 2 wsign/3) a_0 + (4 wsign/15) a_2.
    real(dp) function flux(a, x)
        type(ntv_system_t), intent(in) :: a
        complex(dp), intent(in) :: x(:, 0:)
        flux = pair(a, a%s, x)
    end function flux

    !> Drive (1 + wsign xi^2) S/(b v): source = 2/(3 b) per unit of the
    !> Legendre components 1 + wsign/3 (P_0) and 2 wsign/3 (P_2).
    subroutine drive(a, v, rhs)
        type(ntv_system_t), intent(in) :: a
        real(dp), intent(in) :: v
        complex(dp), allocatable, intent(out) :: rhs(:, :)
        allocate (rhs(2*a%n, 0:a%nblock - 1))
        rhs = 0
        rhs(1:a%n, 0) = (1.5_dp + a%wsign/2)*a%source*a%s/v
        rhs(1:a%n, 1) = a%wsign*a%source*a%s/v
    end subroutine drive

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
