module ideal_map_cases
    ! Exact Solovev states for the N1 tests (test support; module state holds
    ! the active case because the mesh boundary callback has no context).
    !   A3:     psi = c((R^2-R0^2)^2 + 4R^2Z^2 - level), FF' = 0
    !   shaped: psi/c = (R^2-R0^2)^2 + (4/kappa^2) R^2 Z^2
    !           + tau (R^4 - 4 R^2 Z^2 - R0^4) + dd Z^2 - level  (FF' /= 0)
    ! pullback_field_t: the exact state pulled back through an interior shift
    ! Y (Y = id on the LCFS). solovev_profile_t: the exact q(s), p(s) in
    ! normalised poloidal flux (q by spectral contour quadrature, tabulated
    ! on Chebyshev nodes) for the nested-seed builder.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_ideal_map, only: ideal_frozen_t
    use kin6d_ideal_seed, only: ideal_profile_t
    use fortnum_quadrature, only: gauss_legendre_ab
    implicit none
    private
    public :: setup_case, psi_exact, f_of_psi, boundary_point, label, phi_edge
    public :: f_edge, seed_map_point
    public :: pullback_field_t, solovev_profile_t
    public :: r0, pi, a, raxis, psi_axis, pprime, ffprime, shaped
    real(dp), parameter :: r0 = 6.2_dp, f0 = 32.86_dp, c = f0/(12*r0**3)
    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp), parameter :: kappa = 2.2_dp, tau = -0.1_dp, dd = 5.0_dp
    integer, parameter :: ncheb = 48
    logical :: shaped = .false.
    real(dp) :: a, level, raxis, psi_axis, pprime, ffprime, shift(2)
    real(dp) :: qcheb(ncheb)
    type, extends(ideal_frozen_t) :: pullback_field_t
    contains
        procedure :: evaluate => pullback_evaluate
    end type
    type, extends(ideal_profile_t) :: solovev_profile_t
    contains
        procedure :: evaluate => profile_evaluate
    end type
contains
    subroutine setup_case(shaped_case)
        logical, intent(in) :: shaped_case
        real(dp) :: g(2), s
        integer :: j
        shaped = shaped_case
        a = r0/3
        if (shaped) then
            ! axis: d/dR on Z=0 of (R^2-R0^2)^2 + tau R^4 vanishes
            raxis = r0/sqrt(1 + tau)
            level = 0
            call psi_raw((raxis + a), 0.0_dp, level)
            pprime = -(8 + 2*(4/kappa**2))*c
            ffprime = -2*dd*c
        else
            raxis = r0
            level = 4*r0*r0*a*a
            pprime = -16*c
            ffprime = 0
        end if
        call psi_exact(raxis, 0.0_dp, psi_axis, g)
        if (.not. psi_axis < 0) error stop 'axis flux sign'
        shift = [0.10_dp, 0.06_dp]*a
        do j = 1, ncheb
            s = (1 - cos(pi*(j - 0.5_dp)/ncheb))/2
            qcheb(j) = contour_q(s)
        end do
    end subroutine

    subroutine psi_raw(r, z, value)
        real(dp), intent(in) :: r, z
        real(dp), intent(out) :: value
        value = (r*r - r0*r0)**2 + 4/kappa**2*r*r*z*z &
            + tau*(r**4 - 4*r*r*z*z - r0**4) + dd*z*z
    end subroutine

    pure subroutine psi_exact(r, z, value, grad)
        real(dp), intent(in) :: r, z
        real(dp), intent(out) :: value, grad(2)
        if (shaped) then
            value = c*((r*r - r0*r0)**2 + 4/kappa**2*r*r*z*z &
                + tau*(r**4 - 4*r*r*z*z - r0**4) + dd*z*z - level)
            grad = c*[4*r*(r*r - r0*r0) + 8/kappa**2*r*z*z + tau*(4*r**3 - 8*r*z*z), &
                8/kappa**2*r*r*z - 8*tau*r*r*z + 2*dd*z]
        else
            value = c*((r*r - r0*r0)**2 + 4*r*r*z*z - level)
            grad = c*[4*r*(r*r - r0*r0 + 2*z*z), 8*r*r*z]
        end if
    end subroutine

    pure real(dp) function f_of_psi(psi)
        real(dp), intent(in) :: psi
        f_of_psi = sqrt(f0*f0 + 2*ffprime*psi)
    end function

    subroutine boundary_point(t, point)
        ! Exact LCFS: psi = 0 along the ray from the axis at angle 2*pi*t.
        real(dp), intent(in) :: t
        real(dp), intent(out) :: point(2)
        real(dp) :: angle, lo, hi, rho, v, g(2)
        integer :: it
        if (.not. shaped) then
            angle = 2*pi*t
            point(1) = sqrt(r0*r0 + 2*r0*a*cos(angle))
            point(2) = r0*a*sin(angle)/point(1)
            return
        end if
        angle = 2*pi*t
        lo = 1.0e-3_dp*a; hi = 3*a
        do it = 1, 200
            rho = (lo + hi)/2
            call psi_exact(raxis + rho*cos(angle), rho*sin(angle), v, g)
            if (v < 0) then
                lo = rho
            else
                hi = rho
            end if
            if (hi - lo <= 2*spacing(hi)) exit
        end do
        rho = (lo + hi)/2
        point = [raxis + rho*cos(angle), rho*sin(angle)]
    end subroutine

    subroutine seed_map(y, x, dy)
        ! Y(y) = y + shift*psi_exact(y)/psi_axis (identity on the LCFS).
        real(dp), intent(in) :: y(2)
        real(dp), intent(out) :: x(2), dy(2, 2)
        real(dp) :: v, g(2)
        call psi_exact(y(1), y(2), v, g)
        x = y + shift*v/psi_axis
        dy(1, :) = shift(1)*g/psi_axis
        dy(2, :) = shift(2)*g/psi_axis
        dy(1, 1) = dy(1, 1) + 1; dy(2, 2) = dy(2, 2) + 1
    end subroutine

    subroutine pullback_evaluate(field, point, chi, chi_gradient, b3, pressure)
        ! Pull-back of the exact state through Y: chi = psi o Y,
        ! b3 = det(DY) F(psi(Y))/R(Y), p = p(psi(Y)).
        class(pullback_field_t), intent(in) :: field
        real(dp), intent(in) :: point(2)
        real(dp), intent(out) :: chi, chi_gradient(2), b3, pressure
        real(dp) :: x(2), dy(2, 2), v, g(2)
        call seed_map(point, x, dy)
        call psi_exact(x(1), x(2), v, g)
        chi_gradient = matmul(transpose(dy), g)
        b3 = (dy(1, 1)*dy(2, 2) - dy(1, 2)*dy(2, 1))*f_of_psi(v)/x(1)
        pressure = pprime*v
        chi = v
    end subroutine

    real(dp) function label(point)
        real(dp), intent(in) :: point(2)
        real(dp) :: x(2), dy(2, 2), g(2)
        call seed_map(point, x, dy)
        call psi_exact(x(1), x(2), label, g)
    end function

    real(dp) function contour_q(s) result(q)
        ! q = F/(2 pi) oint dl/(R |grad psi|) on psi = psi_axis (1 - s), with the
        ! contour star-shaped about the axis: dl/|grad psi| = rho dtheta/psi_rho.
        real(dp), intent(in) :: s
        real(dp) :: theta, lo, hi, rho, v, g(2), target, dpsi
        integer :: j, it
        target = psi_axis*(1 - s)
        q = 0
        do j = 1, 256
            theta = 2*pi*(j - 1)/256
            lo = 0; hi = 3*a
            do it = 1, 200
                rho = (lo + hi)/2
                call psi_exact(raxis + rho*cos(theta), rho*sin(theta), v, g)
                if (v < target) then
                    lo = rho
                else
                    hi = rho
                end if
                if (hi - lo <= 2*spacing(hi)) exit
            end do
            rho = (lo + hi)/2
            call psi_exact(raxis + rho*cos(theta), rho*sin(theta), v, g)
            dpsi = g(1)*cos(theta) + g(2)*sin(theta)
            q = q + rho/((raxis + rho*cos(theta))*dpsi)
        end do
        q = f_of_psi(target)*q/256
    end function

    subroutine profile_evaluate(profile, s, q, pressure)
        ! Chebyshev (first-kind nodes) barycentric interpolation of q; exact p.
        class(solovev_profile_t), intent(in) :: profile
        real(dp), intent(in) :: s
        real(dp), intent(out) :: q, pressure
        real(dp) :: x, xj, w, num, den
        integer :: j
        x = 2*s - 1
        num = 0; den = 0
        do j = 1, ncheb
            xj = -cos(pi*(j - 0.5_dp)/ncheb)
            if (abs(x - xj) < 4*epsilon(1.0_dp)) then
                q = qcheb(j); exit
            end if
            w = (-1)**(j - 1)*sin(pi*(j - 0.5_dp)/ncheb)/(x - xj)
            num = num + w*qcheb(j); den = den + w
        end do
        if (j > ncheb) q = num/den
        pressure = pprime*psi_axis*(1 - s)
    end subroutine

    real(dp) function f_edge()
        ! Edge (vacuum) F: psi = 0 on the LCFS.
        f_edge = f_of_psi(0.0_dp)
    end function

    subroutine seed_map_point(y, x, dy)
        ! The pullback shift Y (exact minimiser of the pullback case).
        real(dp), intent(in) :: y(2)
        real(dp), intent(out) :: x(2), dy(2, 2)
        call seed_map(y, x, dy)
    end subroutine

    real(dp) function phi_edge()
        ! Total toroidal flux 2 pi |psi_axis| int_0^1 q ds of the exact state.
        type(solovev_profile_t) :: profile
        real(dp) :: nodes(64), weights(64), q, p
        integer :: j
        call gauss_legendre_ab(64, 0.0_dp, 1.0_dp, nodes, weights)
        phi_edge = 0
        do j = 1, 64
            call profile%evaluate(nodes(j), q, p)
            phi_edge = phi_edge + weights(j)*q
        end do
        phi_edge = 2*pi*abs(psi_axis)*phi_edge
    end function
end module ideal_map_cases
