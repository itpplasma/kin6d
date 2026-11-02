module kin6d_nested_seed_input
    ! Boundary/profile input of the fixed-boundary 3D map seed
    ! (kin6d_ideal_seed_3d) from an exact nested reference. Only the boundary
    ! surface, the total toroidal flux and the profiles iota(s), p(s) in the
    ! normalized toroidal flux s enter. Interior oracle fields do not, so they
    ! remain independent error oracles for the solved state.
    ! Both oracle families have two field periods (nested_boundary_harmonic).
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use kin6d_nested_oracles, only: nested_oracle_t, NESTED_INTEGER, NESTED_SHEARED
    use kin6d_nested_oracle_leaves, only: nested_boundary_harmonic, &
        nested_boundary_radius
    use kin6d_ideal_seed_3d, only: toroidal_profile_t
    implicit none
    private
    public :: nested_seed_boundary, nested_seed_profile
    integer, parameter, public :: NESTED_SEED_NFP = 2
    integer, parameter, public :: NESTED_SEED_OK = 0, NESTED_SEED_BAD_INPUT = 1, &
        NESTED_SEED_ORACLE = 2, NESTED_SEED_FIT = 3
    real(dp), parameter :: pi = acos(-1.0_dp)
contains
    subroutine nested_seed_boundary(oracle, mpol, ntor, boundary, status)
        ! Fourier projection of the exact boundary on one field period:
        ! R,Z = sum C cos(m theta-n nfp phi)+S sin(m theta-n nfp phi),
        ! stored as boundary(1:4,n,m) = (rc, rs, zc, zs).
        type(nested_oracle_t), intent(in) :: oracle
        integer, intent(in) :: mpol, ntor
        real(dp), allocatable, intent(out) :: boundary(:, :, :)
        integer, intent(out) :: status
        real(dp), allocatable :: rz(:, :, :)
        real(dp) :: theta, phi, x(3), term(4)
        integer :: angles, j, k, m, n, local, failed
        status = NESTED_SEED_BAD_INPUT
        if (.not. oracle%valid()) return
        if (min(mpol, ntor) < 1) return
        angles = max(64, 4*max(mpol, ntor) + 1)
        allocate (rz(2, 0:angles - 1, 0:angles - 1))
        allocate (boundary(4, -ntor:ntor, 0:mpol))
        failed = 0
        !$omp parallel do collapse(2) schedule(static) &
        !$omp private(j,k,theta,phi,x,local) reduction(+:failed)
        do j = 0, angles - 1
            do k = 0, angles - 1
                theta = 2*pi*k/angles; phi = 2*pi*j/(NESTED_SEED_NFP*angles)
                call oracle%boundary(theta, phi, x, local)
                failed = failed + local
                call nested_boundary_radius(x(1), x(2), rz(1, k, j))
                rz(2, k, j) = x(3)
            end do
        end do
        !$omp end parallel do
        status = NESTED_SEED_ORACLE
        if (failed /= 0) return
        boundary = 0
        !$omp parallel do collapse(2) schedule(static) private(m,n,k,j,theta,phi,term)
        do m = 0, mpol
            do n = -ntor, ntor
                if (m == 0 .and. n < 0) cycle
                do j = 0, angles - 1
                    do k = 0, angles - 1
                        theta = 2*pi*k/angles; phi = 2*pi*j/(NESTED_SEED_NFP*angles)
                        call nested_boundary_harmonic(rz(1, k, j), rz(2, k, j), &
                            theta, phi, real(m, dp), real(n, dp), &
                            term(1), term(2), term(3), term(4))
                        boundary(:, n, m) = boundary(:, n, m) + term
                    end do
                end do
                boundary(:, n, m) = boundary(:, n, m)/angles**2
                if (m /= 0 .or. n /= 0) boundary(:, n, m) = 2*boundary(:, n, m)
            end do
        end do
        !$omp end parallel do
        status = NESTED_SEED_OK
    end subroutine

    subroutine nested_seed_profile(oracle, degree, tolerance, profile, &
            toroidal_flux, fit_error, status)
        ! Matched profiles of the seed: the full toroidal flux Phi(delta) and
        ! iota(s), p(s) with s=Phi(psi)/Phi(delta) (mu0=1 pressure units).
        ! The integer family is exactly polynomial (iota=2, p linear in s).
        ! The sheared family is a least-squares monomial fit of the given
        ! degree on Chebyshev nodes. fit_error is the larger relative maximum
        ! deviation of iota and p on a staggered check set; status
        ! NESTED_SEED_FIT is returned when it exceeds tolerance.
        type(nested_oracle_t), intent(in) :: oracle
        integer, intent(in) :: degree
        real(dp), intent(in) :: tolerance
        type(toroidal_profile_t), intent(out) :: profile
        real(dp), intent(out) :: toroidal_flux, fit_error
        integer, intent(out) :: status
        real(dp), allocatable :: s(:), psi(:), phi(:), iota(:), pressure(:)
        real(dp), allocatable :: a(:, :), rhs(:, :), work(:)
        real(dp) :: fit_iota, fit_pressure, unused
        integer, allocatable :: point_status(:)
        integer :: nodes, checks, k, j, info, ierr, lwork
        status = NESTED_SEED_BAD_INPUT
        fit_error = huge(1.0_dp); toroidal_flux = 0
        if (.not. oracle%valid()) return
        if (degree < 1 .or. degree > 16) return
        if (.not. (tolerance > 0 .and. tolerance < 1)) return
        status = NESTED_SEED_ORACLE
        call oracle%flux(oracle%delta, 1e-3_dp*tolerance, toroidal_flux, ierr)
        if (ierr /= 0 .or. .not. toroidal_flux > 0) return
        nodes = 4*(degree + 1); checks = 2*nodes + 1
        if (oracle%family == NESTED_INTEGER) nodes = 2
        allocate (s(nodes + checks))
        if (oracle%family == NESTED_INTEGER) then
            s(1:2) = [0.0_dp, 1.0_dp]
        else
            do k = 1, nodes
                s(k) = (1 - cos(pi*(k - 0.5_dp)/nodes))/2
            end do
        end if
        do k = 1, checks
            s(nodes + k) = (k - 1)/real(checks - 1, dp)
        end do
        allocate (psi(size(s)), phi(size(s)), iota(size(s)), pressure(size(s)))
        allocate (point_status(size(s)))
        call oracle%profiles(s, 1e-2_dp*tolerance, psi, phi, iota, pressure, &
            point_status)
        if (any(point_status /= 0)) return
        if (.not. all(ieee_is_finite([iota, pressure]))) return
        status = NESTED_SEED_FIT
        if (oracle%family == NESTED_INTEGER) then
            profile%iota = [iota(1)]
            profile%pressure = [pressure(1), pressure(2) - pressure(1)]
        else
            allocate (a(nodes, degree + 1), rhs(nodes, 2), work(1))
            do k = 1, nodes
                a(k, 1) = 1
                do j = 2, degree + 1
                    a(k, j) = a(k, j - 1)*s(k)
                end do
            end do
            rhs(:, 1) = iota(1:nodes); rhs(:, 2) = pressure(1:nodes)
            call dgels('N', nodes, degree + 1, 2, a, nodes, rhs, nodes, work, -1, info)
            lwork = max(1, int(work(1))); deallocate (work); allocate (work(lwork))
            call dgels('N', nodes, degree + 1, 2, a, nodes, rhs, nodes, work, &
                size(work), info)
            if (info /= 0) return
            profile%iota = rhs(1:degree + 1, 1)
            profile%pressure = rhs(1:degree + 1, 2)
        end if
        fit_error = 0
        do k = nodes + 1, nodes + checks
            call profile%evaluate(s(k), fit_iota, fit_pressure, unused)
            fit_error = max(fit_error, &
                abs(fit_iota - iota(k))/maxval(abs(iota(nodes + 1:))), &
                abs(fit_pressure - pressure(k))/maxval(abs(pressure(nodes + 1:))))
        end do
        if (.not. ieee_is_finite(fit_error)) return
        if (fit_error > tolerance) return
        status = NESTED_SEED_OK
    end subroutine
end module
