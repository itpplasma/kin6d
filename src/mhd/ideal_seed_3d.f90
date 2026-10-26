module kin6d_ideal_seed_3d
    ! Fixed-boundary seed only: no exact interior equilibrium is used.
    ! Boundary convention: R,Z = sum C_mn cos(m theta-n*nfp*phi)
    !                                  + S_mn sin(m theta-n*nfp*phi).
    ! The reference disk is (Rref-2,Zref)=(u,v), s=u^2+v^2. Profiles use
    ! normalized TOROIDAL flux, unlike ideal_seed's poloidal-flux profiles.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use kin6d_ideal_map, only: ideal_map_t, ideal_frozen_t, ideal_wall_t, &
        ideal_map_init, ideal_map_freeze, ideal_map_energy, ideal_map_wall_slide, &
        IDEAL_MAP_OK, IDEAL_MAP_BAD_INPUT, IDEAL_MAP_FOLDED
    use kin6d_equilibrium_chart, only: cylindrical_chart_t
    use kin6d_triangle_basis, only: TRIANGLE_MAX_DEGREE
    use kin6d_grad_shafranov, only: gs_state_t, initialize_gs
    use kin6d_ideal_seed_3d_leaves, only: seed_complex_product, seed_boundary_mode, &
        seed_horner, seed_toroidal_form, seed_disk_coordinates, seed_disk_boundary, &
        seed_primitive_coefficient, seed_primitive_finish, seed_wall_mode
    implicit none
    private
    public :: toroidal_profile_t, toroidal_seed_t, toroidal_wall_t, ideal_seed_3d_init
    public :: ideal_seed_disk_mesh, ideal_seed_boundary_coefficients
    real(dp), parameter :: pi = acos(-1.0_dp)
    type :: toroidal_profile_t
        ! Ascending monomial coefficients in s, pressure in solver units.
        real(dp), allocatable :: iota(:), pressure(:)
    contains
        procedure :: evaluate => profile_evaluate
    end type
    type, extends(ideal_wall_t) :: toroidal_wall_t
        ! The boundary series (module header convention) as a sliding wall.
        real(dp), allocatable :: coefficients(:, :, :)
        integer :: nfp = 1
    contains
        procedure :: evaluate => wall_evaluate
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
        real(dp) :: next, coefficient
        integer :: k
        iota = 0; pressure = 0; integral = 0
        do k = size(profile%iota), 1, -1
            call seed_horner(s, profile%iota(k), iota, next)
            iota = next
            call seed_primitive_coefficient(profile%iota(k), real(k, dp), coefficient)
            call seed_horner(s, coefficient, integral, next)
            integral = next
        end do
        call seed_primitive_finish(s, integral, next)
        integral = next
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
        call seed_disk_coordinates(point(1), point(2), u, v, s)
        call field%profile%evaluate(s, iota, p, integral)
        call field%profile%evaluate(1.0_dp, qi, pi_unused, edge)
        call seed_toroidal_form(u, v, field%flux, pi, integral, edge, iota, p, &
            chi, chi_gradient(1), chi_gradient(2), b3, pressure)
    end subroutine

    subroutine wall_evaluate(wall, theta, phi, x, xu, xuu)
        ! Mode sums are orchestration; each term and its theta-derivatives are
        ! the generated seed_wall_mode leaf.
        class(toroidal_wall_t), intent(in) :: wall
        real(dp), intent(in) :: theta, phi
        real(dp), intent(out) :: x(2), xu(2), xuu(2)
        real(dp) :: t(6)
        integer :: m, n, ntor
        ntor = (size(wall%coefficients, 2) - 1)/2
        x = 0; xu = 0; xuu = 0
        do m = 0, ubound(wall%coefficients, 3)
            do n = -ntor, ntor
                call seed_wall_mode(theta, phi, real(m, dp), real(n, dp), &
                    real(wall%nfp, dp), wall%coefficients(1, n + ntor + 1, m), &
                    wall%coefficients(2, n + ntor + 1, m), &
                    wall%coefficients(3, n + ntor + 1, m), &
                    wall%coefficients(4, n + ntor + 1, m), &
                    t(1), t(2), t(3), t(4), t(5), t(6))
                x = x + t([1, 4]); xu = xu + t([2, 5]); xuu = xuu + t([3, 6])
            end do
        end do
    end subroutine

    subroutine ideal_seed_disk_mesh(mesh, boundary_nodes, degree, status)
        type(gs_state_t), intent(out) :: mesh
        integer, intent(in) :: boundary_nodes, degree
        integer, intent(out) :: status
        real(dp) :: points(2, boundary_nodes), h
        integer :: j
        status = IDEAL_MAP_BAD_INPUT
        if (degree < 2 .or. degree > TRIANGLE_MAX_DEGREE) return
        ! Fewer, larger cells are the point of the high orders; the periodic
        ! boundary spline needs four nodes.
        if (boundary_nodes < merge(8, 4, degree <= 3)) return
        do j = 1, boundary_nodes
            call disk_boundary(real(j - 1, dp)/boundary_nodes, points(:, j))
        end do
        h = 2*pi/boundary_nodes
        call initialize_gs(mesh, points, h*h/2, status, degree, curve=disk_boundary, &
            geometry_only=.true.)
    end subroutine

    subroutine disk_boundary(t, point)
        real(dp), intent(in) :: t
        real(dp), intent(out) :: point(2)
        call seed_disk_boundary(t, pi, point(1), point(2))
    end subroutine

    subroutine ideal_seed_boundary_coefficients(point, coefficients, values)
        ! Harmonic extension: rho^m exp(i m theta)=(u+i v)^m is regular
        ! at the axis. Mode sums and their storage are orchestration; the
        ! complex recurrence and real Fourier conversion are generated.
        real(dp), intent(in) :: point(2), coefficients(:, :, 0:)
        real(dp), intent(out) :: values(:, :)
        real(dp) :: u, v, s, re, im, nr, ni, term(4)
        integer :: m, n, ntor, mode
        ntor = (size(coefficients, 2) - 1)/2
        values = 0; re = 1; im = 0
        call seed_disk_coordinates(point(1), point(2), u, v, s)
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
            profile, status, nphi, sliding_wall)
        type(ideal_map_t), intent(out) :: map
        type(toroidal_seed_t), intent(out) :: seed
        type(gs_state_t), intent(in) :: mesh
        real(dp), intent(in) :: coefficients(:, :, 0:), flux
        integer, intent(in) :: nfp
        type(toroidal_profile_t), intent(in) :: profile
        integer, intent(out) :: status
        integer, intent(in), optional :: nphi
        ! Default (pinned): every wall node is held at its input angle, which
        ! is exact only when that angle is the straight-field-line angle of
        ! the target. sliding_wall = .true.: the wall is a shape and its nodes
        ! slide along it (ideal_map_wall_slide). Experimental: the P_k
        ! boundary edges between sliding nodes are not the wall, and the
        ! minimiser enlarges the enclosed area through them (docs/ideal-seed-3d.md).
        logical, intent(in), optional :: sliding_wall
        type(cylindrical_chart_t) :: chart
        type(toroidal_wall_t) :: wall
        real(dp) :: iota, pressure, integral, previous, energy, minimum_radius
        real(dp) :: u, v, s
        real(dp), allocatable :: angle(:)
        integer :: a, j, ntor, angles
        logical :: sliding
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
        ! A fixed physical toroidal coordinate admits both R/Z increments.
        ! Reference-disk radial directions can become tangent to the seeded
        ! surfaces (already for a clockwise circular boundary). The frozen
        ! b3 is positive, so the phi-fixed representative is admissible.
        map%position_gauge = .true.
        call ideal_map_freeze(map, seed, status)
        if (status /= IDEAL_MAP_OK) return
        minimum_radius=huge(1.0_dp)
        !$omp parallel do schedule(static) private(a,j) reduction(min:minimum_radius)
        do a = 1, map%nodes
            call ideal_seed_boundary_coefficients(map%reference(:, a), &
                coefficients, map%x(:, :, a))
            do j=1,map%nphi
                minimum_radius=min(minimum_radius, &
                    dot_product(map%x(1,:,a),map%trig(:,j)))
            end do
        end do
        !$omp end parallel do
        if (.not. minimum_radius>0) then
            status=IDEAL_MAP_BAD_INPUT
            return
        end if
        ! The position gauge holds both (q1, q2) dofs of every wall node and
        ! there is no omega dof. The frozen field has straight field lines in
        ! the reference angle, so the pinned wall is exact only for a boundary
        ! given in the straight-field-line angle.
        sliding = .false.
        if (present(sliding_wall)) sliding = sliding_wall
        if (sliding) then
            wall%coefficients = coefficients; wall%nfp = nfp
            allocate (angle(map%nodes))
            do a = 1, map%nodes
                call seed_disk_coordinates(map%reference(1, a), map%reference(2, a), &
                    u, v, s)
                angle(a) = atan2(v, u)
            end do
            call ideal_map_wall_slide(map, wall, angle, status)
            if (status /= IDEAL_MAP_OK) return
        end if
        ! Sampled positive-Jacobian admission, not a whole-cell certificate.
        call ideal_map_energy(map, map%x, energy, status)
        if (status == IDEAL_MAP_FOLDED) then
            ! The supplied poloidal angle may be clockwise. Orientation is
            ! fixed by the initialized physical map, not the identity disk
            ! used temporarily by ideal_map_init. Mixed signs still fail.
            map%sigma = -map%sigma
            call ideal_map_energy(map, map%x, energy, status)
        end if
    end subroutine
end module
