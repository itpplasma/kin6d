module gs_fb_fixture
    !! Shared pieces of the free-boundary Newton tests (DESIGN 5.1, T3, T4,
    !! T7): normalised profile laws with and without an edge jump, compactly
    !! supported analytic flux bumps (a main plasma bump and an optional
    !! second one that creates an X point), and nodal sampling.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_grad_shafranov, only: gs_state_t
    use kin6d_gs_profile_law, only: gs_profile_law_t
    implicit none
    private

    public :: make_law, bump, bump_laplacian, sample_nodes, smooth_direction

contains

    subroutine make_law(law, edge_jump, scale)
        !! s in [0, 1], p'(s) = scale * 1e5 * q(s), FF'(s) = scale * 1.2 * r(s).
        !! Without a jump both vanish at s = 1 (quadratically); with a jump
        !! they stay finite there.
        type(gs_profile_law_t), intent(out) :: law
        logical, intent(in) :: edge_jump
        real(dp), intent(in) :: scale

        allocate (law%knots(2), law%pp(4, 1), law%ff(4, 1))
        law%knots = [0.0_dp, 1.0_dp]
        law%normalized = .true.
        law%axis = -1.0_dp
        law%p_edge = 0.0_dp
        law%f_edge = 5.0_dp
        if (edge_jump) then
            law%pp(:, 1) = scale*1.0e5_dp*[1.0_dp, -1.4_dp, 0.3_dp, 0.25_dp]
            law%ff(:, 1) = scale*1.2_dp*[1.0_dp, -0.9_dp, -0.4_dp, 0.5_dp]
        else
            law%pp(:, 1) = scale*1.0e5_dp*[1.0_dp, -2.0_dp, 1.0_dp, 0.0_dp]
            law%ff(:, 1) = scale*1.2_dp*[1.0_dp, -2.0_dp, 1.0_dp, 0.0_dp]
        end if
    end subroutine make_law

    pure function bump(r, z, c, semi, amplitude) result(w)
        !! amplitude * (1 - rho^2)^4 inside the ellipse of semi axes `semi`
        !! about c, zero outside (C^3).
        real(dp), intent(in) :: r, z, c(2), semi(2), amplitude
        real(dp) :: w, rho2

        rho2 = ((r - c(1))/semi(1))**2 + ((z - c(2))/semi(2))**2
        w = 0.0_dp
        if (rho2 < 1.0_dp) w = amplitude*(1.0_dp - rho2)**4
    end function bump

    pure function bump_laplacian(r, z, c, semi, amplitude) result(f)
        !! f = -d_R(psi_R / R) - d_Z(psi_Z / R) of the bump: the source
        !! mu0 j_phi that the bump alone would carry.
        real(dp), intent(in) :: r, z, c(2), semi(2), amplitude
        real(dp) :: f, u, ur, uz, urr, uzz, w1, w2, wr, wz, wrr, wzz

        f = 0.0_dp
        u = ((r - c(1))/semi(1))**2 + ((z - c(2))/semi(2))**2
        if (u >= 1.0_dp) return
        ur = 2.0_dp*(r - c(1))/semi(1)**2
        uz = 2.0_dp*(z - c(2))/semi(2)**2
        urr = 2.0_dp/semi(1)**2
        uzz = 2.0_dp/semi(2)**2
        w1 = -4.0_dp*amplitude*(1.0_dp - u)**3
        w2 = 12.0_dp*amplitude*(1.0_dp - u)**2
        wr = w1*ur
        wz = w1*uz
        wrr = w2*ur*ur + w1*urr
        wzz = w2*uz*uz + w1*uzz
        f = -(wrr/r - wr/r**2) - wzz/r
    end function bump_laplacian

    subroutine sample_nodes(state, field, values)
        type(gs_state_t), intent(in) :: state
        interface
            function field(r, z) result(v)
                import :: dp
                real(dp), intent(in) :: r, z
                real(dp) :: v
            end function field
        end interface
        real(dp), allocatable, intent(out) :: values(:)
        integer :: i

        allocate (values(size(state%points, 2)))
        do i = 1, size(values)
            values(i) = field(state%points(1, i), state%points(2, i))
        end do
    end subroutine sample_nodes

    subroutine smooth_direction(state, direction)
        !! A smooth nodal perturbation of unit size.
        type(gs_state_t), intent(in) :: state
        real(dp), allocatable, intent(out) :: direction(:)
        integer :: i
        real(dp) :: r, z

        allocate (direction(size(state%points, 2)))
        do i = 1, size(direction)
            r = state%points(1, i)
            z = state%points(2, i)
            direction(i) = sin(2.3_dp*r + 0.4_dp)*cos(1.7_dp*z) &
                           + 0.3_dp*cos(1.1_dp*r*z) + 0.2_dp*r - 0.5_dp
        end do
    end subroutine smooth_direction

end module gs_fb_fixture
