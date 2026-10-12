module gs_ci_fixture
    !! Twin-experiment fixture for the coil-current inverse tests (T6): eight
    !! circuits around the artificial curve of the coupling fixture, and a
    !! manufactured free-boundary problem (analytic flux bump plus the vacuum
    !! flux of the true currents I*; the psi-independent load `extra` makes
    !! psi* an exact discrete solution at I*, as in the T3 test). The forward
    !! problem keeps `extra` fixed, so it is a genuine free-boundary problem
    !! in the circuit currents.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortfem_axisymmetric_gs_kernel, only: gs_mu0, gs_rect_coil_t
    use kin6d_grad_shafranov, only: gs_state_t
    use kin6d_gs_coils, only: gs_circuit_set_t, gs_circuit_flux_vectors
    use kin6d_gs_coupling, only: gs_coupling_source_load, GS_COUPLING_OK
    use kin6d_gs_profile_law, only: gs_profile_law_t
    use kin6d_gs_free_boundary, only: gs_fb_config_t, gs_fb_system_t, gs_fb_eval_t, &
                      gs_fb_setup, gs_fb_evaluate, gs_fb_newton, GS_FB_OK, GS_FB_LIMITER
    use gs_coupling_fixture, only: build_state, gp_r0, gp_a, gp_b
    use gs_fb_fixture, only: make_law, bump, bump_laplacian
    implicit none
    private

    integer, parameter, public :: n_circ = 8
    real(dp), parameter, public :: istar(n_circ) = [1.2e4_dp, -0.8e4_dp, 0.9e4_dp, &
                                     -1.1e4_dp, 0.5e4_dp, -0.4e4_dp, 1.0e4_dp, 1.0e4_dp]
    real(dp), parameter, public :: current_limit = 4.0e4_dp
    real(dp), parameter :: centre(2) = [gp_r0 + 0.02_dp, 0.03_dp]
    real(dp), parameter :: semi(2) = [gp_a, gp_b]
    real(dp), parameter :: amplitude = 1.5_dp
    type(gs_circuit_set_t), save :: module_set
    !> Plasma source amplitude (Ip about 0.6 MA).
    real(dp), public :: law_scale = 5.0_dp

    public :: build_circuits, build_twin, target_points

contains

    subroutine build_circuits(set)
        !! Eight single-coil circuits, 10 turns, +-4e4 A; circuits 7 and 8
        !! (top and bottom) are the pair that a series tie couples.
        type(gs_circuit_set_t), intent(out) :: set
        real(dp), parameter :: pos(2, n_circ) = reshape([4.7_dp, 1.9_dp, &
                     4.8_dp, -1.8_dp, 1.2_dp, 2.2_dp, 1.3_dp, -2.1_dp, 4.9_dp, 0.0_dp, &
                          1.1_dp, 0.0_dp, 3.0_dp, 2.6_dp, 3.0_dp, -2.6_dp], [2, n_circ])
        integer :: c

        allocate (set%coils(n_circ), set%incidence(n_circ, n_circ), &
                  set%i_min(n_circ), set%i_max(n_circ))
        set%incidence = 0.0_dp
        do c = 1, n_circ
            set%incidence(c, c) = 1.0_dp
            set%coils(c)%rect = gs_rect_coil_t(pos(1, c), pos(2, c), 0.3_dp, 0.3_dp, &
                                               10.0_dp)
        end do
        set%i_min = -current_limit
        set%i_max = current_limit
    end subroutine build_circuits

    subroutine target_points(points)
        !! Eight points on an ellipse inside the plasma.
        real(dp), intent(out) :: points(2, 8)
        real(dp) :: angle
        integer :: k

        do k = 1, 8
            angle = 6.283185307179586_dp*(k - 0.5_dp)/8.0_dp
          points(:, k) = [centre(1) + 0.3_dp*cos(angle), centre(2) + 0.48_dp*sin(angle)]
        end do
    end subroutine target_points

    function vacuum(r, z) result(psi)
        real(dp), intent(in) :: r, z
        real(dp) :: psi, g(1, n_circ)

        call gs_circuit_flux_vectors(module_set, [r], [z], g)
        psi = dot_product(g(1, :), istar)
    end function vacuum

    function density(r, z) result(j)
        real(dp), intent(in) :: r, z
        real(dp) :: j

        j = bump_laplacian(r, z, centre, semi, amplitude)/gs_mu0
    end function density

    subroutine build_twin(n_outer, order, state, sys, set, psi_star, status)
        !! The mesh, the free-boundary system (limiter mode, edge jump) with the
        !! manufactured load, and psi* = the Newton solution at I*.
        integer, intent(in) :: n_outer, order
        type(gs_state_t), intent(out) :: state
        type(gs_fb_system_t), intent(out) :: sys
        type(gs_circuit_set_t), intent(out) :: set
        real(dp), allocatable, intent(out) :: psi_star(:)
        integer, intent(out) :: status
        type(gs_profile_law_t) :: law
        type(gs_fb_config_t) :: cfg
        type(gs_fb_eval_t) :: ev
        real(dp), allocatable :: exact(:), seed(:), f_load(:)
        integer :: i, n, stat

        status = 1
        call build_circuits(set)
        module_set = set
        call build_state(state, n_outer, order, stat)
        if (stat .ne. 0) return
        n = size(state%points, 2)
        allocate (exact(n), seed(n))
        do i = 1, n
            exact(i) = vacuum(state%points(1, i), state%points(2, i)) &
                       + bump(state%points(1, i), state%points(2, i), centre, semi, &
                              amplitude)
            seed(i) = vacuum(state%points(1, i), state%points(2, i)) &
                      + bump(state%points(1, i), state%points(2, i), &
                             centre + [0.06_dp, -0.05_dp], semi, 0.5_dp*amplitude)
        end do
        call make_law(law, .true., law_scale)
        cfg%mode = GS_FB_LIMITER
        cfg%axis_sign = 1
        cfg%axis_guess = centre
        cfg%search_radius = 0.4_dp
        cfg%quadrature_tolerance = 1.0e-11_dp
        cfg%jacobian_tolerance = 1.0e-11_dp
        cfg%tolerance = 1.0e-10_dp
        allocate (cfg%points(2, 1))
        cfg%points(:, 1) = [gp_r0 + 0.33_dp, 0.0_dp]
        call gs_fb_setup(state, sys, law, cfg, set, stat)
        if (stat .ne. GS_FB_OK) return
        call gs_coupling_source_load(state, density, f_load, stat)
        if (stat .ne. GS_COUPLING_OK) return
        call gs_fb_evaluate(state, sys, exact, istar, ev, stat, jacobian=.false.)
        if (stat .ne. GS_FB_OK) return
        sys%extra_load = f_load - ev%source
        state%psi = seed
        call gs_fb_newton(state, sys, istar, stat)
        if (stat .ne. GS_FB_OK) return
        psi_star = state%psi
        status = 0
    end subroutine build_twin

end module gs_ci_fixture
