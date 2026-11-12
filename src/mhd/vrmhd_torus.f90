module milca_vrmhd_torus
    !! Axisymmetric state of the circular-torus flux-pumping case: the
    !! cylinder case's p(psi_N) and FF'(psi_N) in a fixed-boundary
    !! Grad-Shafranov solve. See docs/vrmhd-torus.md.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use milca_vrmhd_cylinder, only: vrmhd_cylinder_case_t, vrmhd_cylinder_profile_t
    use milca_gs_profile_law, only: gs_profile_law_t, gs_law_primitives
    use milca_grad_shafranov, only: gs_state_t, initialize_gs, gs_axis_candidate, &
        gs_axis_safety_factor
    use milca_gs_nonlinear, only: solve_gs_cubic_profiles
    implicit none
    private
    real(dp), parameter :: pi = acos(-1.0_dp)

    !> Axis of the solved state; shift = R_axis - R0.
    type, public :: vrmhd_torus_axis_t
        real(dp) :: position(2) = 0, flux = 0, shift = 0, q = 0
    end type

    ! The boundary callback of initialize_gs takes no context.
    real(dp), save :: circle_major = 0

    public :: vrmhd_torus_law, vrmhd_torus_equilibrium
contains
    subroutine vrmhd_torus_law(c, prof, knots, law, status)
        !! Normalised law, piecewise linear in s = psi_N on `knots` uniform
        !! knots. psi = R0 int B_theta dr is the cylinder flux at R = R0, so
        !! dp/dpsi = T'/(R0 B_theta) and, from radial force balance and
        !! Ampere's law, FF' = -R0 (J_z + T'/B_theta).
        type(vrmhd_cylinder_case_t), intent(in) :: c
        type(vrmhd_cylinder_profile_t), intent(in) :: prof
        integer, intent(in) :: knots
        type(gs_profile_law_t), intent(out) :: law
        integer, intent(out) :: status
        real(dp) :: pp(knots), ff(knots), s, w, dtr, btr, jz
        integer :: i, j, n
        status = 1
        n = size(prof%x)
        if (knots < 2 .or. n < 2) return
        j = 1
        do i = 1, knots
            s = real(i - 1, dp)/(knots - 1)
            do while (j < n - 1 .and. prof%psin(j + 1) < s)
                j = j + 1
            end do
            w = (s - prof%psin(j))/(prof%psin(j + 1) - prof%psin(j))
            dtr = (1 - w)*prof%dt_over_r(j) + w*prof%dt_over_r(j + 1)
            btr = (1 - w)*prof%btheta_over_r(j) + w*prof%btheta_over_r(j + 1)
            jz = (1 - w)*prof%jz(j) + w*prof%jz(j + 1)
            pp(i) = dtr/(c%r0*btr)
            ff(i) = -c%r0*(jz + dtr/btr)
        end do
        allocate (law%knots(knots), law%pp(4, knots - 1), law%ff(4, knots - 1))
        law%knots = [(real(i - 1, dp)/(knots - 1), i=1, knots)]
        law%pp = 0; law%ff = 0
        law%pp(1, :) = pp(:knots - 1); law%ff(1, :) = ff(:knots - 1)
        law%pp(2, :) = (pp(2:) - pp(:knots - 1))*(knots - 1)
        law%ff(2, :) = (ff(2:) - ff(:knots - 1))*(knots - 1)
        law%normalized = .true.
        law%p_edge = prof%t(n)
        law%f_edge = c%r0*prof%bz(n)
        law%axis = -c%r0*prof%psi(n)
        if (law%valid()) status = 0
    end subroutine

    subroutine vrmhd_torus_equilibrium(c, prof, edge, law, state, axis, status)
        !! Fixed-boundary solve in the circle of radius a = 1 about R0 with
        !! cubic elements of edge length about `edge`; mu0 = 1. Newton starts
        !! from the cylinder flux; from the solver's default start (psi = 0)
        !! it fails for these profiles.
        type(vrmhd_cylinder_case_t), intent(in) :: c
        type(vrmhd_cylinder_profile_t), intent(in) :: prof
        real(dp), intent(in) :: edge
        type(gs_profile_law_t), intent(inout) :: law
        type(gs_state_t), intent(out) :: state
        type(vrmhd_torus_axis_t), intent(out) :: axis
        integer, intent(out) :: status
        real(dp), allocatable :: polygon(:, :)
        real(dp) :: hessian(2, 2), pressure, f, r, w
        integer :: nodes, i, j, cell, n
        status = 1
        if (edge <= 0 .or. edge > 1) return
        circle_major = c%r0
        nodes = max(8, nint(2*pi/edge))
        allocate (polygon(2, nodes))
        do i = 1, nodes
            call circle(real(i - 1, dp)/nodes, polygon(:, i))
        end do
        call initialize_gs(state, polygon, edge*edge/2, status, 3, curve=circle)
        if (status /= 0) return
        n = size(prof%x)
        do i = 1, size(state%psi)
            if (state%boundary(i)) then
                state%psi(i) = 0
                cycle
            end if
            r = min(norm2(state%points(:, i) - [c%r0, 0.0_dp]), 1.0_dp)*(n - 1)
            j = min(int(r) + 1, n - 1); w = r - (j - 1)
            state%psi(i) = c%r0*((1 - w)*prof%psi(j) + w*prof%psi(j + 1) - prof%psi(n))
        end do
        call solve_gs_cubic_profiles(state, 1.0_dp, law, status, 1.0e-11_dp, &
            quadrature_tolerance=1.0e-11_dp)
        if (status /= 0) return
        call gs_axis_candidate(state, axis%position, axis%flux, hessian, cell, status)
        if (status /= 0) return
        call gs_law_primitives(law, axis%flux, pressure, f, status)
        if (status /= 0) return
        call gs_axis_safety_factor(axis%position(1), f, hessian, axis%q, status)
        axis%shift = axis%position(1) - c%r0
    end subroutine

    pure subroutine circle(t, point)
        real(dp), intent(in) :: t
        real(dp), intent(out) :: point(2)
        point = [circle_major + cos(2*pi*t), sin(2*pi*t)]
    end subroutine
end module milca_vrmhd_torus
