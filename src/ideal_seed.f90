module kin6d_ideal_seed
    ! Nested seed of the frozen-in energy map from a boundary and profiles
    ! (DESIGN 7.1, VMEC's interpolated initial guess made regular at the axis).
    ! The boundary b(t) = R + iZ, t in [0, 1), counter-clockwise, is extended
    ! harmonically into the unit disk, zeta = xi1 + i xi2:
    !     Y(zeta) = sum_{m>=0} a_m zeta^m + sum_{m>=1} a_{-m} conj(zeta)^m
    !               + (c - a_0)(1 - |zeta|^2),
    ! with a_m the Fourier coefficients of b; the last term moves the seed axis
    ! Y(0) to a chosen centre c. Y is polynomial in (xi1, xi2), so the labels
    ! are smooth at the seed axis (a ray map y = c + rho (b - c) is only C1
    ! there unless the boundary is a circle about c, which caps the map's
    ! convergence order). Profiles are q(s) and p(s) in the normalised poloidal
    ! flux s = |zeta|^2 with the total toroidal flux Phi_edge, as in KIN6D's 2D
    ! prescribed-q solve:
    !     psi_axis = -Phi_edge/(2 pi int_0^1 q ds),   chi = psi_axis (1 - s),
    !     b3 = 2 |psi_axis| q(s) / det DY,
    ! i.e. toroidal flux 2 pi |psi_axis| q ds between seed surfaces. Y is a
    ! diffeomorphism for convex boundaries (Rado-Kneser-Choquet); init checks
    ! det DY > 0 on a polar grid.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_ideal_map, only: ideal_frozen_t
    use kin6d_curved_triangle, only: boundary_curve
    use fortnum_quadrature, only: gauss_legendre_ab
    implicit none
    private
    public :: ideal_profile_t, nested_seed_t, ideal_seed_init

    ! q(s) and p(s) in normalised poloidal flux.
    type, abstract :: ideal_profile_t
    contains
        procedure(profile_evaluate), deferred :: evaluate
    end type
    abstract interface
        subroutine profile_evaluate(profile, s, q, pressure)
            import :: ideal_profile_t, dp
            class(ideal_profile_t), intent(in) :: profile
            real(dp), intent(in) :: s
            real(dp), intent(out) :: q, pressure
        end subroutine
    end interface

    type, extends(ideal_frozen_t) :: nested_seed_t
        class(ideal_profile_t), allocatable :: profile
        ! a(-modes:modes) boundary Fourier coefficients; shift = c - a_0.
        complex(dp), allocatable :: a(:)
        complex(dp) :: shift = 0
        integer :: modes = 0
        real(dp) :: psi_axis = 0, scale = 1
    contains
        procedure :: evaluate => seed_evaluate
        procedure :: map => seed_map
        procedure :: inverse => seed_inverse
    end type
    real(dp), parameter :: pi = acos(-1.0_dp)
    integer, parameter :: samples = 1024
contains
    subroutine ideal_seed_init(seed, boundary, phi_edge, profile, status, centre)
        type(nested_seed_t), intent(out) :: seed
        procedure(boundary_curve) :: boundary
        real(dp), intent(in) :: phi_edge
        class(ideal_profile_t), intent(in) :: profile
        integer, intent(out) :: status
        real(dp), intent(in), optional :: centre(2)
        real(dp) :: nodes(64), weights(64), q, p, integral, point(2), y(2), dy(2, 2)
        complex(dp) :: b(0:samples - 1)
        integer :: j, m, k
        status = 1
        allocate (seed%profile, source=profile)
        do j = 0, samples - 1
            call boundary(real(j, dp)/samples, point)
            b(j) = cmplx(point(1), point(2), dp)
        end do
        seed%scale = maxval(abs(b - sum(b)/samples))
        ! Keep the modes above roundoff; the rest of the spectrum is dropped.
        allocate (seed%a(-samples/2 + 1:samples/2 - 1))
        do m = lbound(seed%a, 1), ubound(seed%a, 1)
            seed%a(m) = sum(b*exp(cmplx(0, -2*pi*m*[(k, k=0, samples - 1)]/samples, &
                dp)))/samples
        end do
        seed%modes = 0
        do m = 1, ubound(seed%a, 1)
            if (max(abs(seed%a(m)), abs(seed%a(-m))) > 1.0e-13_dp*seed%scale) &
                seed%modes = m
        end do
        status = 2
        if (seed%modes > samples/4) return
        block
            complex(dp), allocatable :: kept(:)
            allocate (kept(-seed%modes:seed%modes))
            kept = seed%a(-seed%modes:seed%modes)
            call move_alloc(kept, seed%a)
        end block
        seed%shift = 0
        if (present(centre)) seed%shift = cmplx(centre(1), centre(2), dp) - seed%a(0)
        call gauss_legendre_ab(64, 0.0_dp, 1.0_dp, nodes, weights)
        integral = 0
        do j = 1, 64
            call profile%evaluate(nodes(j), q, p)
            integral = integral + weights(j)*q
        end do
        status = 3
        if (.not. integral > 0 .or. .not. phi_edge > 0) return
        status = 4
        seed%psi_axis = -phi_edge/(2*pi*integral)
        do j = 0, 32
            do k = 0, 63
                call seed%map(j/32.0_dp*[cos(2*pi*k/64), sin(2*pi*k/64)], y, dy)
                if (.not. dy(1, 1)*dy(2, 2) - dy(1, 2)*dy(2, 1) > 0) return
            end do
        end do
        status = 0
    end subroutine

    subroutine seed_map(seed, xi, y, dy)
        ! y = Y(xi) and dy(:, i) = dY/dxi_i.
        class(nested_seed_t), intent(in) :: seed
        real(dp), intent(in) :: xi(2)
        real(dp), intent(out) :: y(2), dy(2, 2)
        complex(dp) :: zeta, value, dz, dzbar, power, powerbar
        integer :: m
        zeta = cmplx(xi(1), xi(2), dp)
        value = seed%a(0); dz = 0; dzbar = 0
        power = 1; powerbar = 1
        do m = 1, seed%modes
            dz = dz + m*seed%a(m)*power
            dzbar = dzbar + m*seed%a(-m)*powerbar
            power = power*zeta; powerbar = powerbar*conjg(zeta)
            value = value + seed%a(m)*power + seed%a(-m)*powerbar
        end do
        value = value + seed%shift*(1 - xi(1)**2 - xi(2)**2)
        ! d/dxi1 = d/dzeta + d/dzetabar, d/dxi2 = i (d/dzeta - d/dzetabar).
        y = [real(value), aimag(value)]
        dy(:, 1) = cplx(dz + dzbar - 2*xi(1)*seed%shift)
        dy(:, 2) = cplx(cmplx(0, 1, dp)*(dz - dzbar) - 2*xi(2)*seed%shift)
    contains
        pure function cplx(z) result(v)
            complex(dp), intent(in) :: z
            real(dp) :: v(2)
            v = [real(z), aimag(z)]
        end function
    end subroutine

    subroutine seed_inverse(seed, y, xi, ok)
        ! Damped Newton for Y(xi) = y from the seed axis.
        class(nested_seed_t), intent(in) :: seed
        real(dp), intent(in) :: y(2)
        real(dp), intent(out) :: xi(2)
        logical, intent(out) :: ok
        real(dp) :: point(2), dy(2, 2), step(2), r, rt, trial(2), t, det
        integer :: it, k
        xi = 0
        call seed%map(xi, point, dy)
        r = norm2(point - y)
        ok = .false.
        do it = 1, 100
            if (r <= 1.0e-15_dp*seed%scale) then
                ok = .true.; return
            end if
            det = dy(1, 1)*dy(2, 2) - dy(1, 2)*dy(2, 1)
            step = [dy(2, 2)*(y(1) - point(1)) - dy(1, 2)*(y(2) - point(2)), &
                -dy(2, 1)*(y(1) - point(1)) + dy(1, 1)*(y(2) - point(2))]/det
            t = 1
            do k = 1, 30
                trial = xi + t*step
                call seed%map(trial, point, dy)
                rt = norm2(point - y)
                if (rt < r) exit
                t = t/2
            end do
            if (.not. rt < r) then
                ok = r <= 1.0e-13_dp*seed%scale
                call seed%map(xi, point, dy)
                return
            end if
            xi = trial; r = rt
        end do
        ok = r <= 1.0e-13_dp*seed%scale
    end subroutine

    subroutine seed_evaluate(field, point, chi, chi_gradient, b3, pressure)
        class(nested_seed_t), intent(in) :: field
        real(dp), intent(in) :: point(2)
        real(dp), intent(out) :: chi, chi_gradient(2), b3, pressure
        real(dp) :: xi(2), y(2), dy(2, 2), det, s, q, grad_xi(2)
        logical :: ok
        call field%inverse(point, xi, ok)
        if (.not. ok) error stop 'nested seed: point outside the seed map'
        call field%map(xi, y, dy)
        s = min(1.0_dp, sum(xi**2))
        det = dy(1, 1)*dy(2, 2) - dy(1, 2)*dy(2, 1)
        chi = field%psi_axis*(1 - s)
        ! grad_y chi = DY^-T grad_xi chi, grad_xi chi = -2 psi_axis xi.
        grad_xi = -2*field%psi_axis*xi
        chi_gradient = [dy(2, 2)*grad_xi(1) - dy(2, 1)*grad_xi(2), &
            -dy(1, 2)*grad_xi(1) + dy(1, 1)*grad_xi(2)]/det
        call field%profile%evaluate(s, q, pressure)
        b3 = 2*abs(field%psi_axis)*q/det
    end subroutine
end module kin6d_ideal_seed
