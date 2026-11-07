!> Monoenergetic drift-kinetic coefficients D11, D31, D13, D33 on an
!! axisymmetric surface with Lorentz pitch-angle scattering (#182).
!!
!! Equation for g(theta, xi) per drive s, v = 1, B = B0 b(theta):
!!   xi k d_theta g - (1 - xi^2)/2 k d_theta(ln b) d_xi g
!!     - (nu/2) d_xi((1 - xi^2) d_xi g) = s,       k = b.grad theta = kappa b,
!! with s1 = v_d.grad psi (radial drift, u = I/Omega_0 = 1) and s3 = xi b.
!! This is the simpler reference discretisation that drift-kinetics.md asks
!! DK1 to run: Legendre polynomials P_l(xi), in which the Lorentz operator is
!! diagonal and its null space is P_0, and Fourier collocation in theta,
!! exact for a field given by finitely many harmonics. The system is block
!! tridiagonal in l. It is factorised once by block elimination from the
!! highest l down (mono_factor); the remaining l = 0 block is singular by the
!! density null space and is bordered: a density source mu enters the l = 0
!! row and <g_0> = c closes the system (c = 0 for the coefficients). mu is
!! zero for a drive with vanishing density moment, so |mu| measures the
!! discrete solvability defect. Every further right-hand side (drives,
!! tangents, refinement corrections) reuses the factorisation (mono_solve).
!!
!! Coefficients with the flux-surface average <.> (weight 1/b^2) and the
!! pitch average int dxi/2:
!!   D11 = <s1 g1>, D31 = <xi b g1>, D13 = <s1 g3>, D33 = <xi b g3>.
!!
!! Derivatives (mono_coefficients_jvp) follow the implicit rule
!! A dg = ds - dA g through the same factorisation; A depends linearly on the
!! leaf outputs and on nu, and the leaf tangents come from the FortSym JVP.
!! Error estimates (mono_estimate) are indicators in the sense of
!! certification.md: a two-resolution difference for the discretisation and
!! one refinement step on the true residual for the algebraic error.
module kin6d_dk_mono
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_dk_mono_leaves, only: dk_mono_point, dk_mono_point_jvp
    implicit none
    private
    public :: mono_coefficients_t, mono_coefficients, mono_coefficients_jvp
    public :: mono_estimate_t, mono_estimate, trig_resample
    ! Per-speed building blocks for solvers that couple speeds through
    ! collisions: factorisation, bordered solve, operator action, drives.
    public :: mono_point_t, mono_factor_t, mono_factor, mono_solve, mono_apply, mono_drives
    public :: mono_coupling
    integer, parameter, public :: MONO_OK = 0, MONO_BAD_INPUT = 1, &
        MONO_SINGULAR = 2

    type :: mono_coefficients_t
        real(dp) :: d11 = 0, d31 = 0, d13 = 0, d33 = 0
        ! Density-source multipliers of the radial and parallel drives.
        real(dp) :: border(2) = 0
    end type

    !> Error indicators of the coefficients at the requested resolution:
    !> discretisation (value at doubled resolution minus value) and algebraic
    !> (functional of one refinement step on the true residual).
    type :: mono_estimate_t
        type(mono_coefficients_t) :: discretisation, algebraic
    end type

    !> Leaf outputs on the theta points, the flux-surface weight normalised
    !> to sum 1 and its unnormalised sum.
    type :: mono_point_t
        real(dp), allocatable :: stream(:), mirror(:), drift0(:), drift2(:), &
            drive_e(:), weight(:)
        real(dp) :: weight_sum = 0
    end type

    !> Block elimination of A: LU of S_l and Y_l = S_l^-1 Lo_l for l >= 1,
    !> and the LU of the bordered l = 0 system.
    type :: mono_factor_t
        integer :: n = 0, lmax = 0
        real(dp) :: nu = 0
        real(dp), allocatable :: d(:, :), s_lu(:, :, :), y(:, :, :), border_lu(:, :)
        integer, allocatable :: s_piv(:, :), border_piv(:)
        type(mono_point_t) :: p
    end type

    real(dp), parameter :: pi = 4*atan(1.0_dp)

    interface
        subroutine dgetrf(m, n, a, lda, ipiv, info)
            import :: dp
            integer, intent(in) :: m, n, lda
            real(dp), intent(inout) :: a(lda, *)
            integer, intent(out) :: ipiv(*), info
        end subroutine dgetrf
        subroutine dgetrs(trans, n, nrhs, a, lda, ipiv, b, ldb, info)
            import :: dp
            character, intent(in) :: trans
            integer, intent(in) :: n, nrhs, lda, ldb, ipiv(*)
            real(dp), intent(in) :: a(lda, *)
            real(dp), intent(inout) :: b(ldb, *)
            integer, intent(out) :: info
        end subroutine dgetrs
    end interface

contains

    !> b on n_theta (odd) uniform points theta_j = 2 pi (j - 1)/n_theta,
    !! kappa = iota B0/(G + iota I) and the deflection frequency nu = nu_D/v,
    !! both per unit length; n_leg Legendre polynomials P_0 ... P_{n_leg-1}.
    subroutine mono_coefficients(b, kappa, nu, n_leg, c, status)
        real(dp), intent(in) :: b(:), kappa, nu
        integer, intent(in) :: n_leg
        type(mono_coefficients_t), intent(out) :: c
        integer, intent(out) :: status
        type(mono_factor_t) :: f
        real(dp), allocatable :: g(:, :, :)

        call mono_factor(b, kappa, nu, n_leg, f, status)
        if (status /= MONO_OK) return
        call drive_solutions(f, g, c%border, status)
        if (status /= MONO_OK) return
        c = functionals(f%p, g, c%border)
    end subroutine mono_coefficients

    !> Coefficients c and their directional derivatives dc for the tangents
    !> db (on the same points), dkappa and dnu.
    subroutine mono_coefficients_jvp(b, db, kappa, dkappa, nu, dnu, n_leg, c, dc, &
            status)
        real(dp), intent(in) :: b(:), db(:), kappa, dkappa, nu, dnu
        integer, intent(in) :: n_leg
        type(mono_coefficients_t), intent(out) :: c, dc
        integer, intent(out) :: status
        type(mono_factor_t) :: f
        type(mono_point_t) :: dp_
        real(dp), allocatable :: g(:, :, :), dg(:, :, :), rhs(:, :, :), dbt(:)
        real(dp) :: dmu(2), constraint(2)
        integer :: k

        status = MONO_BAD_INPUT
        if (size(db) /= size(b)) return
        call mono_factor(b, kappa, nu, n_leg, f, status)
        if (status /= MONO_OK) return
        call drive_solutions(f, g, c%border, status)
        if (status /= MONO_OK) return
        c = functionals(f%p, g, c%border)

        ! Leaf tangents, then the tangent of the normalised weight.
        allocate (dbt(f%n))
        dbt = matmul(f%d, db)
        call allocate_point(dp_, f%n)
        call dk_mono_point_jvp(b, matmul(f%d, b), kappa, db, dbt, dkappa, &
            dp_%stream, dp_%mirror, dp_%drift0, dp_%drift2, dp_%drive_e, dp_%weight)
        dp_%weight = (dp_%weight - f%p%weight*sum(dp_%weight))/f%p%weight_sum

        ! A dg = ds - dA g and <w dg_0> = -<dw g_0>, with the same factors.
        allocate (rhs(f%n, 0:f%lmax, 2))
        rhs = drives(dp_, f%lmax)
        do k = 1, 2
            rhs(:, :, k) = rhs(:, :, k) - apply_operator(f%d, dp_, dnu, g(:, :, k))
            constraint(k) = -sum(dp_%weight*g(:, 0, k))
        end do
        call mono_solve(f, rhs, constraint, dg, dmu, status)
        if (status /= MONO_OK) return

        block
            type(mono_coefficients_t) :: part
            type(mono_point_t) :: fixed
            fixed = f%p
            fixed%weight = dp_%weight
            ! Product rule: weight tangent, then drive tangents, then dg.
            dc = functionals_weighted(fixed, f%p, g)
            part = functionals_weighted(f%p, dp_, g)
            call add(dc, part)
            part = functionals(f%p, dg, dmu)
            call add(dc, part)
        end block
        dc%border = dmu
    end subroutine mono_coefficients_jvp

    !> Coefficients c at (n_theta, n_leg) with error indicators: the
    !> discretisation indicator compares with (2 n_theta + 1, 2 n_leg), b
    !> resampled by its trigonometric interpolant; the algebraic indicator is
    !> the functional of one solve with the true residual.
    subroutine mono_estimate(b, kappa, nu, n_leg, c, e, status)
        real(dp), intent(in) :: b(:), kappa, nu
        integer, intent(in) :: n_leg
        type(mono_coefficients_t), intent(out) :: c
        type(mono_estimate_t), intent(out) :: e
        integer, intent(out) :: status
        type(mono_factor_t) :: f
        type(mono_coefficients_t) :: fine
        real(dp), allocatable :: g(:, :, :), r(:, :, :), dg(:, :, :)
        real(dp) :: dmu(2), constraint(2)
        integer :: k

        call mono_factor(b, kappa, nu, n_leg, f, status)
        if (status /= MONO_OK) return
        call drive_solutions(f, g, c%border, status)
        if (status /= MONO_OK) return
        c = functionals(f%p, g, c%border)

        allocate (r(f%n, 0:f%lmax, 2))
        r = drives(f%p, f%lmax)
        do k = 1, 2
            r(:, :, k) = r(:, :, k) - apply_operator(f%d, f%p, nu, g(:, :, k))
            r(:, 0, k) = r(:, 0, k) - c%border(k)
            constraint(k) = -sum(f%p%weight*g(:, 0, k))
        end do
        call mono_solve(f, r, constraint, dg, dmu, status)
        if (status /= MONO_OK) return
        e%algebraic = functionals(f%p, dg, dmu)

        call mono_coefficients(trig_resample(b, 2*size(b) + 1), kappa, nu, 2*n_leg, &
            fine, status)
        if (status /= MONO_OK) return
        e%discretisation = fine
        call add(e%discretisation, c, -1.0_dp)
    end subroutine mono_estimate

    subroutine mono_factor(b, kappa, nu, n_leg, f, status)
        real(dp), intent(in) :: b(:), kappa, nu
        integer, intent(in) :: n_leg
        type(mono_factor_t), intent(out) :: f
        integer, intent(out) :: status
        integer :: n, l, j, info
        real(dp), allocatable :: s(:, :), up(:, :), weight(:)

        n = size(b)
        status = MONO_BAD_INPUT
        if (mod(n, 2) == 0 .or. n < 3 .or. n_leg < 3 .or. nu <= 0 .or. &
            abs(kappa) <= tiny(kappa) .or. any(b <= 0)) return
        f%n = n
        f%lmax = n_leg - 1
        f%nu = nu
        allocate (f%d(n, n), weight(n))
        call fourier_derivative(n, f%d)
        call allocate_point(f%p, n)
        call dk_mono_point(b, matmul(f%d, b), kappa, f%p%stream, f%p%mirror, &
            f%p%drift0, f%p%drift2, f%p%drive_e, weight)
        f%p%weight_sum = sum(weight)
        f%p%weight = weight/f%p%weight_sum

        ! S_lmax = nu lmax (lmax + 1)/2, S_l = Di_l - Up_l Y_{l+1}.
        allocate (f%s_lu(n, n, f%lmax), f%s_piv(n, f%lmax), f%y(n, n, f%lmax), &
            s(n, n), up(n, n))
        status = MONO_SINGULAR
        do l = f%lmax, 1, -1
            s = 0
            do j = 1, n
                s(j, j) = nu*l*(l + 1)/2.0_dp
            end do
            if (l < f%lmax) then
                up = upper(f%d, f%p, l)
                s = s - matmul(up, f%y(:, :, l + 1))
            end if
            call dgetrf(n, n, s, n, f%s_piv(:, l), info)
            if (info /= 0) return
            f%s_lu(:, :, l) = s
            f%y(:, :, l) = lower(f%d, f%p, l)
            call dgetrs('N', n, n, s, n, f%s_piv(:, l), f%y(:, :, l), n, info)
        end do

        ! l = 0: S_0 = -Up_0 Y_1, bordered by the density source and <w g_0>.
        allocate (f%border_lu(n + 1, n + 1), f%border_piv(n + 1))
        f%border_lu = 0
        f%border_lu(1:n, 1:n) = -matmul(upper(f%d, f%p, 0), f%y(:, :, 1))
        f%border_lu(1:n, n + 1) = 1
        f%border_lu(n + 1, 1:n) = f%p%weight
        call dgetrf(n + 1, n + 1, f%border_lu, n + 1, f%border_piv, info)
        if (info /= 0) return
        status = MONO_OK
    end subroutine mono_factor

    !> Solve A x + mu e_0 = rhs, <w x_0> = constraint for each right-hand side.
    subroutine mono_solve(f, rhs, constraint, x, mu, status)
        type(mono_factor_t), intent(in) :: f
        real(dp), intent(in) :: rhs(:, 0:, :), constraint(:)
        real(dp), allocatable, intent(out) :: x(:, :, :)
        real(dp), intent(out) :: mu(:)
        integer, intent(out) :: status
        real(dp), allocatable :: z(:, :, :), r(:, :), b0(:, :)
        integer :: l, n, m, info

        n = f%n
        m = size(rhs, 3)
        allocate (z(n, m, 0:f%lmax), r(n, m), x(n, 0:f%lmax, m), b0(n + 1, m))
        status = MONO_SINGULAR
        do l = f%lmax, 1, -1
            r = rhs(:, l, :)
            if (l < f%lmax) r = r - matmul(upper(f%d, f%p, l), z(:, :, l + 1))
            call dgetrs('N', n, m, f%s_lu(:, :, l), n, f%s_piv(:, l), r, n, info)
            if (info /= 0) return
            z(:, :, l) = r
        end do
        b0(1:n, :) = rhs(:, 0, :) - matmul(upper(f%d, f%p, 0), z(:, :, 1))
        b0(n + 1, :) = constraint
        call dgetrs('N', n + 1, m, f%border_lu, n + 1, f%border_piv, b0, n + 1, info)
        if (info /= 0) return
        z(:, :, 0) = b0(1:n, :)
        mu = b0(n + 1, :)
        do l = 1, f%lmax
            z(:, :, l) = z(:, :, l) - matmul(f%y(:, :, l), z(:, :, l - 1))
        end do
        do l = 0, f%lmax
            x(:, l, :) = z(:, :, l)
        end do
        status = MONO_OK
    end subroutine mono_solve

    !> A x without the border, with deflection frequency nu per unit length.
    pure function mono_apply(f, nu, x) result(y)
        type(mono_factor_t), intent(in) :: f
        real(dp), intent(in) :: nu, x(:, 0:)
        real(dp) :: y(size(x, 1), 0:ubound(x, 2))
        y = apply_operator(f%d, f%p, nu, x)
    end function mono_apply

    !> Coefficients of x_{l-1} (lower) and x_{l+1} (upper) in row l of the
    !> operator: the theta blocks of streaming and mirror force, for solvers
    !> that eliminate in l with several speeds per block.
    pure subroutine mono_coupling(f, l, lower_l, upper_l)
        type(mono_factor_t), intent(in) :: f
        integer, intent(in) :: l
        real(dp), intent(out) :: lower_l(f%n, f%n), upper_l(f%n, f%n)
        lower_l = 0
        upper_l = 0
        if (l > 0) lower_l = lower(f%d, f%p, l)
        upper_l = upper(f%d, f%p, l)
    end subroutine mono_coupling

    !> Radial (1) and parallel (2) drives of the factorised surface.
    pure function mono_drives(f) result(s)
        type(mono_factor_t), intent(in) :: f
        real(dp) :: s(f%n, 0:f%lmax, 2)
        s = drives(f%p, f%lmax)
    end function mono_drives

    subroutine drive_solutions(f, g, mu, status)
        type(mono_factor_t), intent(in) :: f
        real(dp), allocatable, intent(out) :: g(:, :, :)
        real(dp), intent(out) :: mu(2)
        integer, intent(out) :: status
        call mono_solve(f, drives(f%p, f%lmax), [0.0_dp, 0.0_dp], g, mu, status)
    end subroutine drive_solutions

    !> Legendre components of the radial and parallel drives, row by row.
    pure function drives(p, lmax) result(s)
        type(mono_point_t), intent(in) :: p
        integer, intent(in) :: lmax
        real(dp) :: s(size(p%stream), 0:lmax, 2)
        s = 0
        s(:, 0, 1) = p%drift0
        s(:, 2, 1) = p%drift2
        s(:, 1, 2) = p%drive_e
    end function drives

    !> A x without the border: Lo_l x_{l-1} + nu l(l + 1)/2 x_l + Up_l x_{l+1}.
    !> A is linear in (stream, mirror, nu), so the same routine applies dA.
    pure function apply_operator(d, p, nu, x) result(y)
        real(dp), intent(in) :: d(:, :), nu, x(:, 0:)
        type(mono_point_t), intent(in) :: p
        real(dp) :: y(size(x, 1), 0:ubound(x, 2))
        integer :: l, lmax
        lmax = ubound(x, 2)
        do l = 0, lmax
            y(:, l) = nu*l*(l + 1)/2.0_dp*x(:, l)
            if (l > 0) y(:, l) = y(:, l) + matmul(lower(d, p, l), x(:, l - 1))
            if (l < lmax) y(:, l) = y(:, l) + matmul(upper(d, p, l), x(:, l + 1))
        end do
    end function apply_operator

    !> Coefficient of x_{l+1} in row l: xi P_{l+1} and (1 - xi^2) P'_{l+1}
    !> have P_l components (l + 1)/(2l + 3) and (l + 1)(l + 2)/(2l + 3).
    pure function upper(d, p, l) result(m)
        real(dp), intent(in) :: d(:, :)
        type(mono_point_t), intent(in) :: p
        integer, intent(in) :: l
        real(dp) :: m(size(d, 1), size(d, 1))
        integer :: i
        do i = 1, size(d, 1)
            m(i, :) = p%stream(i)*(l + 1)/real(2*l + 3, dp)*d(i, :)
            m(i, i) = m(i, i) - p%mirror(i)*(l + 1)*(l + 2)/real(2*l + 3, dp)
        end do
    end function upper

    !> Coefficient of x_{l-1} in row l: components l/(2l - 1) and
    !> -(l - 1) l/(2l - 1).
    pure function lower(d, p, l) result(m)
        real(dp), intent(in) :: d(:, :)
        type(mono_point_t), intent(in) :: p
        integer, intent(in) :: l
        real(dp) :: m(size(d, 1), size(d, 1))
        integer :: i
        do i = 1, size(d, 1)
            m(i, :) = p%stream(i)*l/real(2*l - 1, dp)*d(i, :)
            m(i, i) = m(i, i) + p%mirror(i)*(l - 1)*l/real(2*l - 1, dp)
        end do
    end function lower

    pure function functionals(p, g, mu) result(c)
        type(mono_point_t), intent(in) :: p
        real(dp), intent(in) :: g(:, 0:, :), mu(2)
        type(mono_coefficients_t) :: c
        c = functionals_weighted(p, p, g)
        c%border = mu
    end function functionals

    !> D_ij with the weight of w and the drive factors of s.
    pure function functionals_weighted(w, s, g) result(c)
        type(mono_point_t), intent(in) :: w, s
        real(dp), intent(in) :: g(:, 0:, :)
        type(mono_coefficients_t) :: c
        c%d11 = sum(w%weight*(s%drift0*g(:, 0, 1) + s%drift2*g(:, 2, 1)/5))
        c%d31 = sum(w%weight*s%drive_e*g(:, 1, 1))/3
        c%d13 = sum(w%weight*(s%drift0*g(:, 0, 2) + s%drift2*g(:, 2, 2)/5))
        c%d33 = sum(w%weight*s%drive_e*g(:, 1, 2))/3
    end function functionals_weighted

    subroutine add(a, b, factor)
        type(mono_coefficients_t), intent(inout) :: a
        type(mono_coefficients_t), intent(in) :: b
        real(dp), intent(in), optional :: factor
        real(dp) :: s
        s = 1
        if (present(factor)) s = factor
        a%d11 = a%d11 + s*b%d11
        a%d31 = a%d31 + s*b%d31
        a%d13 = a%d13 + s*b%d13
        a%d33 = a%d33 + s*b%d33
        a%border = a%border + s*b%border
    end subroutine add

    subroutine allocate_point(p, n)
        type(mono_point_t), intent(out) :: p
        integer, intent(in) :: n
        allocate (p%stream(n), p%mirror(n), p%drift0(n), p%drift2(n), &
            p%drive_e(n), p%weight(n))
    end subroutine allocate_point

    !> Samples on m (odd) uniform points of the trigonometric interpolant of
    !> the n (odd) samples b, exact for b band-limited to |k| <= (n - 1)/2.
    pure function trig_resample(b, m) result(c)
        real(dp), intent(in) :: b(:)
        integer, intent(in) :: m
        real(dp) :: c(m)
        integer :: n, k, j, i
        real(dp) :: ak, bk, t
        n = size(b)
        c = sum(b)/n
        do k = 1, (n - 1)/2
            ak = 0
            bk = 0
            do j = 1, n
                t = 2*pi*(j - 1)/n
                ak = ak + b(j)*cos(k*t)
                bk = bk + b(j)*sin(k*t)
            end do
            ak = 2*ak/n
            bk = 2*bk/n
            do i = 1, m
                t = 2*pi*(i - 1)/m
                c(i) = c(i) + ak*cos(k*t) + bk*sin(k*t)
            end do
        end do
    end function trig_resample

    !> Fourier collocation derivative on n (odd) uniform points of [0, 2 pi).
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

end module kin6d_dk_mono
