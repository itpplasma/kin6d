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
!! tridiagonal in l and is solved directly by block elimination from the
!! highest l down; the remaining l = 0 block is singular by the density null
!! space and is bordered: a density source mu enters the l = 0 row and
!! <g_0> = 0 closes the system. mu is zero for a drive with vanishing
!! density moment, so |mu| measures the discrete solvability defect.
!!
!! Coefficients with the flux-surface average <.> (weight 1/b^2) and the
!! pitch average int dxi/2:
!!   D11 = <s1 g1>, D31 = <xi b g1>, D13 = <s1 g3>, D33 = <xi b g3>.
module kin6d_dk_mono
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_dk_mono_leaves, only: dk_mono_point
    implicit none
    private
    public :: mono_coefficients_t, mono_coefficients
    integer, parameter, public :: MONO_OK = 0, MONO_BAD_INPUT = 1, &
        MONO_SINGULAR = 2

    type :: mono_coefficients_t
        real(dp) :: d11 = 0, d31 = 0, d13 = 0, d33 = 0
        ! Density-source multipliers of the radial and parallel drives.
        real(dp) :: border(2) = 0
    end type

    real(dp), parameter :: pi = 4*atan(1.0_dp)

    interface
        subroutine dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
            import :: dp
            integer, intent(in) :: n, nrhs, lda, ldb
            real(dp), intent(inout) :: a(lda, *), b(ldb, *)
            integer, intent(out) :: ipiv(*), info
        end subroutine dgesv
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
        integer :: n, l, lmax, info, j
        integer, allocatable :: ipiv(:)
        real(dp), allocatable :: d(:, :), db(:), stream(:), mirror(:), drift0(:), &
            drift2(:), drive_e(:), weight(:), y(:, :, :), z(:, :, :), s(:, :), &
            lo(:, :), up(:, :), r(:, :), a(:, :, :), border(:, :), rhs(:, :)

        n = size(b)
        status = MONO_BAD_INPUT
        if (mod(n, 2) == 0 .or. n < 3 .or. n_leg < 3 .or. nu <= 0 .or. &
            abs(kappa) <= tiny(kappa) .or. any(b <= 0)) return
        lmax = n_leg - 1
        allocate (d(n, n), db(n), stream(n), mirror(n), drift0(n), drift2(n), &
            drive_e(n), weight(n), ipiv(n + 1))
        call fourier_derivative(n, d)
        db = matmul(d, b)
        call dk_mono_point(b, db, kappa, stream, mirror, drift0, drift2, drive_e, &
            weight)
        weight = weight/sum(weight)

        ! Y_l = S_l^-1 Lo_l and Z_l = S_l^-1 r_l from l = lmax down to l = 1,
        ! with S_lmax = nu lmax (lmax + 1)/2 and
        ! S_l = Di_l - Up_l Y_{l+1}, r_l = s_l - Up_l Z_{l+1}.
        allocate (y(n, n, lmax), z(n, 2, 0:lmax), s(n, n), lo(n, n), up(n, n), &
            r(n, 2))
        do l = lmax, 1, -1
            s = 0
            do j = 1, n
                s(j, j) = nu*l*(l + 1)/2.0_dp
            end do
            r = source(l)
            if (l < lmax) then
                up = upper(l)
                s = s - matmul(up, y(:, :, l + 1))
                r = r - matmul(up, z(:, :, l + 1))
            end if
            lo = lower(l)
            allocate (rhs(n, n + 2))
            rhs(:, 1:n) = lo
            rhs(:, n + 1:n + 2) = r
            call dgesv(n, n + 2, s, n, ipiv, rhs, n, info)
            if (info /= 0) then
                status = MONO_SINGULAR
                return
            end if
            y(:, :, l) = rhs(:, 1:n)
            z(:, :, l) = rhs(:, n + 1:n + 2)
            deallocate (rhs)
        end do

        ! l = 0: S_0 = -Up_0 Y_1, bordered by the density source and <g_0> = 0.
        up = upper(0)
        allocate (border(n + 1, n + 1), rhs(n + 1, 2), a(n, 2, 0:lmax))
        border = 0
        border(1:n, 1:n) = -matmul(up, y(:, :, 1))
        border(1:n, n + 1) = 1
        border(n + 1, 1:n) = weight
        rhs = 0
        rhs(1:n, :) = source(0) - matmul(up, z(:, :, 1))
        call dgesv(n + 1, 2, border, n + 1, ipiv, rhs, n + 1, info)
        if (info /= 0) then
            status = MONO_SINGULAR
            return
        end if
        a(:, :, 0) = rhs(1:n, :)
        c%border = rhs(n + 1, :)
        do l = 1, lmax
            a(:, :, l) = z(:, :, l) - matmul(y(:, :, l), a(:, :, l - 1))
        end do

        c%d11 = sum(weight*(drift0*a(:, 1, 0) + drift2*a(:, 1, 2)/5))
        c%d31 = sum(weight*drive_e*a(:, 1, 1))/3
        c%d13 = sum(weight*(drift0*a(:, 2, 0) + drift2*a(:, 2, 2)/5))
        c%d33 = sum(weight*drive_e*a(:, 2, 1))/3
        status = MONO_OK
    contains
        ! Legendre components of the two drives in row l.
        function source(l) result(sl)
            integer, intent(in) :: l
            real(dp) :: sl(n, 2)
            sl = 0
            if (l == 0) sl(:, 1) = drift0
            if (l == 2) sl(:, 1) = drift2
            if (l == 1) sl(:, 2) = drive_e
        end function source

        ! Coefficient of a_{l+1} in row l: xi P_{l+1} and (1 - xi^2) P'_{l+1}
        ! have P_l components (l + 1)/(2l + 3) and (l + 1)(l + 2)/(2l + 3).
        function upper(l) result(m)
            integer, intent(in) :: l
            real(dp) :: m(n, n)
            integer :: i
            do i = 1, n
                m(i, :) = stream(i)*(l + 1)/real(2*l + 3, dp)*d(i, :)
                m(i, i) = m(i, i) - mirror(i)*(l + 1)*(l + 2)/real(2*l + 3, dp)
            end do
        end function upper

        ! Coefficient of a_{l-1} in row l: components l/(2l - 1) and
        ! -(l - 1) l/(2l - 1).
        function lower(l) result(m)
            integer, intent(in) :: l
            real(dp) :: m(n, n)
            integer :: i
            do i = 1, n
                m(i, :) = stream(i)*l/real(2*l - 1, dp)*d(i, :)
                m(i, i) = m(i, i) + mirror(i)*(l - 1)*l/real(2*l - 1, dp)
            end do
        end function lower
    end subroutine mono_coefficients

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
