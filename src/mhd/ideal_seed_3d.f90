module kin6d_ideal_seed_3d
    ! Fixed-boundary seed only: no exact interior equilibrium is used.
    ! Boundary convention: R,Z = sum C_mn cos(m theta-n*nfp*phi)
    !                                  + S_mn sin(m theta-n*nfp*phi).
    ! The reference disk is (Rref-2,Zref)=(u,v), s=u^2+v^2. Profiles use
    ! normalized TOROIDAL flux, unlike ideal_seed's poloidal-flux profiles.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use kin6d_ideal_map, only: ideal_map_t, ideal_frozen_t, ideal_map_init, &
        ideal_map_freeze, ideal_map_energy, ideal_map_rdof, IDEAL_MAP_OK, &
        IDEAL_MAP_BAD_INPUT
    use kin6d_equilibrium_chart, only: cylindrical_chart_t
    use kin6d_grad_shafranov, only: gs_state_t, initialize_gs
    use kin6d_ideal_seed_3d_leaves, only: seed_complex_product, seed_boundary_mode, &
        seed_horner, seed_toroidal_form
    implicit none
    private
    public :: toroidal_profile_t, toroidal_seed_t, ideal_seed_3d_init
    public :: ideal_seed_disk_mesh, ideal_seed_boundary_coefficients
    real(dp), parameter :: pi = acos(-1.0_dp)
    type :: toroidal_profile_t
        ! Ascending monomial coefficients in s, pressure in solver units.
        real(dp), allocatable :: iota(:), pressure(:)
    contains
        procedure :: evaluate => profile_evaluate
    end type
    type, extends(ideal_frozen_t) :: toroidal_seed_t
        type(toroidal_profile_t) :: profile
        real(dp) :: flux = 0
    contains
        procedure :: evaluate => frozen_evaluate
    end type
contains
    pure subroutine profile_evaluate(profile, s, iota, pressure, integral)
        class(toroidal_profile_t), intent(in) :: profile
        real(dp), intent(in) :: s
        real(dp), intent(out) :: iota, pressure, integral
        real(dp) :: next
        integer :: k
        iota = 0; pressure = 0; integral = 0
        do k = size(profile%iota), 1, -1
            call seed_horner(s, profile%iota(k), iota, next)
            iota = next
            call seed_horner(s, profile%iota(k)/k, integral, next)
            integral = next
        end do
        integral = s*integral
        do k = size(profile%pressure), 1, -1
            call seed_horner(s, profile%pressure(k), pressure, next)
            pressure = next
        end do
    end subroutine

    subroutine frozen_evaluate(field, point, chi, chi_gradient, b3, pressure)
        class(toroidal_seed_t), intent(in) :: field
        real(dp), intent(in) :: point(2)
        real(dp), intent(out) :: chi, chi_gradient(2), b3, pressure
        real(dp) :: u, v, s, iota, p, integral, edge, qi, pi_unused
        u = point(1) - 2; v = point(2); s = u*u + v*v
        call field%profile%evaluate(s, iota, p, integral)
        call field%profile%evaluate(1.0_dp, qi, pi_unused, edge)
        call seed_toroidal_form(u, v, field%flux, pi, integral, edge, iota, p, &
            chi, chi_gradient(1), chi_gradient(2), b3, pressure)
    end subroutine

    subroutine ideal_seed_disk_mesh(mesh, boundary_nodes, degree, status)
        type(gs_state_t), intent(out) :: mesh
        integer, intent(in) :: boundary_nodes, degree
        integer, intent(out) :: status
        real(dp) :: points(2, boundary_nodes), h
        integer :: j
        status = IDEAL_MAP_BAD_INPUT
        if (boundary_nodes < 8) return
        if (degree < 2 .or. degree > 3) return
        do j = 1, boundary_nodes
            call disk_boundary(real(j - 1, dp)/boundary_nodes, points(:, j))
        end do
        h = 2*pi/boundary_nodes
        call initialize_gs(mesh, points, h*h/2, status, degree, curve=disk_boundary)
    end subroutine

    subroutine disk_boundary(t, point)
        real(dp), intent(in) :: t
        real(dp), intent(out) :: point(2)
        point = [2 + cos(2*pi*t), sin(2*pi*t)]
    end subroutine

    subroutine ideal_seed_boundary_coefficients(point, coefficients, values)
        ! Harmonic extension: rho^m exp(i m theta)=(u+i v)^m is regular
        ! at the axis. Mode sums and their storage are orchestration; the
        ! complex recurrence and real Fourier conversion are generated.
        real(dp), intent(in) :: point(2), coefficients(:, :, 0:)
        real(dp), intent(out) :: values(:, :)
        real(dp) :: u, v, re, im, nr, ni, term(4)
        integer :: m, n, ntor, mode
        ntor = (size(coefficients, 2) - 1)/2
        values = 0; re = 1; im = 0
        u = point(1) - 2; v = point(2)
        do m = 0, ubound(coefficients, 3)
            do n = -ntor, ntor
                call seed_boundary_mode(re, im, &
                    coefficients(1, n + ntor + 1, m), &
                    coefficients(2, n + ntor + 1, m), &
                    coefficients(3, n + ntor + 1, m), &
                    coefficients(4, n + ntor + 1, m), &
                    term(1), term(2), term(3), term(4))
                mode = 2*abs(n)
                if (n == 0) then
                    values(1:2, 1) = values(1:2, 1) + term([1, 3])
                else
                    values(1:2, mode) = values(1:2, mode) + term([1, 3])
                    values(1:2, mode + 1) = values(1:2, mode + 1) &
                        + sign(1.0_dp, real(n, dp))*term([2, 4])
                end if
            end do
            call seed_complex_product(re, im, u, v, nr, ni)
            re = nr; im = ni
        end do
    end subroutine

    subroutine ideal_seed_3d_init(map, seed, mesh, coefficients, nfp, flux, &
            profile, status, nphi)
        type(ideal_map_t), intent(out) :: map
        type(toroidal_seed_t), intent(out) :: seed
        type(gs_state_t), intent(in) :: mesh
        real(dp), intent(in) :: coefficients(:, :, 0:), flux
        integer, intent(in) :: nfp
        type(toroidal_profile_t), intent(in) :: profile
        integer, intent(out) :: status
        integer, intent(in), optional :: nphi
        type(cylindrical_chart_t) :: chart
        real(dp) :: iota, pressure, integral, previous, energy
        integer :: a, m, ntor, angles
        status = IDEAL_MAP_BAD_INPUT
        if (nfp < 1 .or. .not. flux > 0) return
        if (.not. ieee_is_finite(flux)) return
        if (size(coefficients, 1) /= 4) return
        if (mod(size(coefficients, 2), 2) /= 1) return
        if (.not. all(ieee_is_finite(coefficients))) return
        if (.not. allocated(profile%iota)) return
        if (.not. allocated(profile%pressure)) return
        if (min(size(profile%iota), size(profile%pressure)) < 1) return
        if (.not. all(ieee_is_finite(profile%iota))) return
        if (.not. all(ieee_is_finite(profile%pressure))) return
        previous = 0
        do a = 0, 128
            call profile%evaluate(a/128.0_dp, iota, pressure, integral)
            if (.not. all(ieee_is_finite([iota, pressure, integral]))) return
            if (abs(iota) <= tiny(iota)) return
            if (a > 0) then
                if ((iota > 0) .neqv. (previous > 0)) return
            end if
            previous = iota
        end do
        ntor = (size(coefficients, 2) - 1)/2
        angles = max(1, 4*ntor + 1)
        if (present(nphi)) angles = nphi
        call ideal_map_init(map, mesh, chart, nfp, ntor, status, nphi=angles)
        if (status /= IDEAL_MAP_OK) return
        seed%profile = profile; seed%flux = flux
        call ideal_map_freeze(map, seed, status)
        if (status /= IDEAL_MAP_OK) return
        !$omp parallel do schedule(static) private(a,m)
        do a = 1, map%nodes
            call ideal_seed_boundary_coefficients(map%reference(:, a), &
                coefficients, map%x(:, :, a))
            if (.not. map%boundary(a)) cycle
            ! A nonaxisymmetric fixed wall is parametrized by phi-hat. A
            ! boundary omega would move its points off that prescribed wall.
            do m = 1, map%nmode
                map%fixed(ideal_map_rdof(map, 2, m, a)) = .true.
            end do
        end do
        !$omp end parallel do
        ! Sampled positive-Jacobian admission, not a whole-cell certificate.
        call ideal_map_energy(map, map%x, energy, status)
    end subroutine
end module
