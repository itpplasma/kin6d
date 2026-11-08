program test_gs_fb_law_tangent
    !! Implicit tangent of the converged free-boundary state with respect to
    !! the profile-law amplitudes (pp and ff scaled separately) against
    !! central differences of re-converged Newton solves: the nodal flux, the
    !! axis and boundary flux levels and the plasma current. A negative
    !! control drops the plasma response (dpsi = 0) and must miss the
    !! differences by far more. The manufactured problem of T3 supplies a
    !! genuine free-boundary equilibrium with an edge jump (interface term)
    !! and the separatrix with a vanishing edge source (the FreeGS setting).
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use milca_grad_shafranov, only: gs_state_t
    use milca_gs_coils, only: gs_circuit_set_t, gs_circuit_flux_vectors
    use milca_gs_coupling, only: gs_coupling_source_load, GS_COUPLING_OK
    use milca_gs_profile_law, only: gs_profile_law_t
    use milca_gs_free_boundary, only: gs_fb_config_t, gs_fb_system_t, gs_fb_eval_t, &
        gs_fb_setup, gs_fb_evaluate, gs_fb_newton, gs_fb_law_tangent, GS_FB_OK, &
        GS_FB_LIMITER, GS_FB_XPOINT
    use gs_coupling_fixture, only: build_state, build_coils, gp_r0, gp_a, gp_b
    use gs_fb_fixture, only: make_law, bump, bump_laplacian
    implicit none

    real(dp), parameter :: currents(2) = [1.5e4_dp, -1.0e4_dp]
    real(dp), parameter :: centre(2) = [gp_r0 + 0.02_dp, 0.03_dp]
    real(dp), parameter :: semi(2) = [gp_a, gp_b]
    real(dp), parameter :: low_c(2) = [gp_r0, -1.15_dp], low_s(2) = [0.5_dp, 0.45_dp]
    type(gs_circuit_set_t) :: coils
    logical :: with_second
    integer :: failures

    failures = 0
    call build_coils(coils)
    with_second = .false.
    call run_case('limiter, edge jump, P3', GS_FB_LIMITER, .true.)
    with_second = .true.
    call run_case('separatrix, vanishing edge, P3', GS_FB_XPOINT, .false.)
    if (failures /= 0) error stop 'gs_fb_law_tangent failures'
    print '(a)', 'gs_fb_law_tangent oracles passed'

contains

    subroutine require(ok, what)
        logical, intent(in) :: ok
        character(len=*), intent(in) :: what
        if (.not. ok) then
            failures = failures + 1
            print '(2a)', 'FAIL: ', what
        end if
    end subroutine require

    function exact_flux(r, z) result(psi)
        real(dp), intent(in) :: r, z
        real(dp) :: psi, g(1, 2)

        call gs_circuit_flux_vectors(coils, [r], [z], g)
        psi = g(1, 1)*currents(1) + g(1, 2)*currents(2) &
            + bump(r, z, centre, semi, 1.5_dp)
        if (with_second) psi = psi + bump(r, z, low_c, low_s, 0.9_dp)
    end function exact_flux

    function density(r, z) result(j)
        use fortfem_axisymmetric_gs_kernel, only: gs_mu0
        real(dp), intent(in) :: r, z
        real(dp) :: j

        j = bump_laplacian(r, z, centre, semi, 1.5_dp)
        if (with_second) j = j + bump_laplacian(r, z, low_c, low_s, 0.9_dp)
        j = j/gs_mu0
    end function density

    subroutine run_case(label, mode, jump)
        character(len=*), intent(in) :: label
        integer, intent(in) :: mode
        logical, intent(in) :: jump
        real(dp), parameter :: eps = 1.0e-3_dp, tol = 1.0e-7_dp
        type(gs_state_t) :: state
        type(gs_profile_law_t) :: law, directions(2)
        type(gs_fb_config_t) :: cfg
        type(gs_fb_system_t) :: sys
        type(gs_fb_eval_t) :: ev
        real(dp), allocatable :: exact(:), f_load(:), root(:), dpsi(:, :), fd(:, :)
        real(dp) :: daxis(2), dbound(2), dcurrent(2), levels(3, 2), plus(3), minus(3)
        real(dp) :: error_flux, error_levels, error_current, control, frozen(2)
        real(dp), allocatable :: probe_flux(:)
        integer :: status, i, k

        call build_state(state, 32, 3, status)
        call require(status == 0, 'mesh')
        allocate (exact(size(state%points, 2)))
        do i = 1, size(exact)
            exact(i) = exact_flux(state%points(1, i), state%points(2, i))
        end do
        call make_law(law, jump, 0.7_dp)
        cfg%mode = mode
        cfg%axis_sign = 1
        cfg%axis_guess = centre
        cfg%search_radius = 0.4_dp
        cfg%quadrature_tolerance = 1.0e-12_dp
        cfg%jacobian_tolerance = 1.0e-12_dp
        cfg%tolerance = 1.0e-13_dp
        if (mode == GS_FB_XPOINT) then
            cfg%xpoint_guess = [gp_r0, -0.62_dp]
            cfg%search_radius = 0.45_dp
        else
            allocate (cfg%points(2, 1))
            cfg%points(:, 1) = [gp_r0 + 0.33_dp, 0.0_dp]
        end if
        call gs_fb_setup(state, sys, law, cfg, coils, status)
        call require(status == GS_FB_OK, 'setup')
        call gs_coupling_source_load(state, density, f_load, status)
        call require(status == GS_COUPLING_OK, 'forcing load')
        call gs_fb_evaluate(state, sys, exact, currents, ev, status, jacobian=.false.)
        call require(status == GS_FB_OK, 'evaluate at the exact flux')
        if (failures /= 0) return
        sys%extra_load = f_load - ev%source
        state%psi = exact
        call gs_fb_newton(state, sys, currents, status)
        call require(status == GS_FB_OK, 'reference solve')
        root = state%psi

        directions = law
        directions(1)%ff = 0.0_dp
        directions(2)%pp = 0.0_dp
        allocate (dpsi(size(root), 2), fd(size(root), 2), probe_flux(size(root)))
        call gs_fb_law_tangent(state, sys, root, currents, directions, dpsi, &
            daxis, dbound, dcurrent, status)
        call require(status == GS_FB_OK, 'law tangent')
        if (failures /= 0) return

        do k = 1, 2
            call solve_scaled(state, sys, law, root, label, k, 1.0_dp + eps, .true., plus, fd(:, k))
            call solve_scaled(state, sys, law, root, label, k, 1.0_dp - eps, .true., minus, probe_flux)
            levels(:, k) = (plus - minus)/(2.0_dp*eps)
            fd(:, k) = (fd(:, k) - probe_flux)/(2.0_dp*eps)
            ! Negative control: the same law change at frozen flux.
            call solve_scaled(state, sys, law, root, label, k, 1.0_dp + eps, .false., plus, probe_flux)
            call solve_scaled(state, sys, law, root, label, k, 1.0_dp - eps, .false., minus, probe_flux)
            frozen(k) = (plus(3) - minus(3))/(2.0_dp*eps)
        end do
        error_flux = maxval(abs(dpsi - fd))/maxval(abs(fd))
        error_levels = maxval(abs([daxis, dbound] - [levels(1, :), levels(2, :)])) &
            /maxval(abs(levels(1:2, :)))
        error_current = maxval(abs(dcurrent - levels(3, :))/abs(levels(3, :)))
        control = maxval(abs(frozen - levels(3, :))/abs(levels(3, :)))
        print '(2a,3es10.2)', label, ': flux, levels, Ip tangent errors ', &
            error_flux, error_levels, error_current
        print '(a,es10.2)', '   frozen-flux Ip control error ', control
        call require(error_flux < tol, label//': flux tangent')
        call require(error_levels < tol, label//': flux-level tangents')
        call require(error_current < tol, label//': plasma-current tangent')
        call require(control > 100.0_dp*tol, label//': plasma response matters')

    end subroutine run_case

    subroutine solve_scaled(state, sys, law, root, label, component, factor, &
            resolve, observed, flux)
        !! Amplitude `component` scaled by `factor`; with `resolve` the
        !! state is re-converged, otherwise the flux stays at the root.
        type(gs_state_t), intent(inout) :: state
        type(gs_fb_system_t), intent(inout) :: sys
        type(gs_profile_law_t), intent(in) :: law
        real(dp), intent(in) :: root(:)
        character(len=*), intent(in) :: label
        integer, intent(in) :: component
        real(dp), intent(in) :: factor
        logical, intent(in) :: resolve
        real(dp), intent(out) :: observed(3), flux(:)
        type(gs_fb_eval_t) :: probe
        integer :: stat

        sys%law = law
        if (component == 1) sys%law%pp = factor*law%pp
        if (component == 2) sys%law%ff = factor*law%ff
        state%psi = root
        if (resolve) then
            call gs_fb_newton(state, sys, currents, stat)
            call require(stat == GS_FB_OK, label//': perturbed solve')
        end if
        flux = state%psi
        call gs_fb_evaluate(state, sys, flux, currents, probe, stat, &
            jacobian=.false.)
        call require(stat == GS_FB_OK, label//': perturbed evaluation')
        observed = [probe%psi_axis, probe%psi_bound, probe%plasma_current]
        sys%law = law
        state%psi = root
    end subroutine solve_scaled

end program test_gs_fb_law_tangent
