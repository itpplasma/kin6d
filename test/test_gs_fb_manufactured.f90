program test_gs_fb_manufactured
    !! T3 (DESIGN 5.1): a manufactured free-boundary solution and the
    !! quadratic convergence of the Newton iteration.
    !!
    !! psi* = sum_k I_k g_k + w with w a compactly supported analytic bump
    !! inside Gamma_c. The coil part is Delta*-harmonic inside Gamma_c, so the
    !! source of psi* is f*(R, Z) = -d_R(w_R/R) - d_Z(w_Z/R), known in closed
    !! form. The solver is given the nodal load
    !!
    !!     extra = int f* v dA - b(Pi psi*),
    !!
    !! with b the plasma source with the cut at psi_b of the interpolant of
    !! psi*; the exact LCFS {psi* = psi_b*} is fixed by the limiter point (or
    !! the mean of the target points). Newton then solves
    !! A psi - b(psi) - sum I load - extra = 0 from a seed that has a different
    !! amplitude and centre, and must reproduce psi* up to the discretisation
    !! error of the coupled operator. Run with an edge jump of the source and
    !! with a law that vanishes at the edge, for the limiter and the
    !! mean-at-targets modes, P2 and P3, on a mesh sequence.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_grad_shafranov, only: gs_state_t
    use kin6d_gs_coils, only: gs_circuit_set_t, gs_circuit_flux_vectors
    use kin6d_gs_coupling, only: gs_coupling_source_load, GS_COUPLING_OK
    use kin6d_gs_profile_law, only: gs_profile_law_t
    use kin6d_gs_free_boundary, only: gs_fb_config_t, gs_fb_system_t, gs_fb_eval_t, &
                 gs_fb_history_t, gs_fb_setup, gs_fb_evaluate, gs_fb_newton, GS_FB_OK, &
                                      GS_FB_LIMITER, GS_FB_TARGET_MEAN, GS_FB_XPOINT
    use gs_coupling_fixture, only: build_state, build_coils, gp_r0, gp_a, gp_b, fe_value
    use gs_fb_fixture, only: make_law, bump, bump_laplacian
    implicit none

    real(dp), parameter :: currents(2) = [1.5e4_dp, -1.0e4_dp]
    real(dp), parameter :: centre(2) = [gp_r0 + 0.02_dp, 0.03_dp]
    real(dp), parameter :: semi(2) = [gp_a, gp_b]
    real(dp), parameter :: amplitude = 1.5_dp
    real(dp), parameter :: low_c(2) = [gp_r0, -1.15_dp], low_s(2) = [0.5_dp, 0.45_dp]
    type(gs_circuit_set_t) :: coils
    logical :: with_second = .false.
    real(dp) :: case_fraction = 0.8_dp
    integer :: failures

    failures = 0
    call build_coils(coils)
    call sweep('limiter, edge jump, P2', 2, [24, 48, 96], GS_FB_LIMITER, .true., 2.0_dp)
    call sweep('limiter, vanishing edge, P2', 2, [24, 48, 96], GS_FB_LIMITER, &
               .false., 2.0_dp)
    call sweep('targets, edge jump, P3', 3, [24, 48], GS_FB_TARGET_MEAN, .true., 3.0_dp)
    call sweep('targets, vanishing edge, P3', 3, [24, 48], GS_FB_TARGET_MEAN, &
               .false., 3.0_dp)
    ! A second flux bump below the plasma creates a saddle: FB-cut mode with
    ! psi_b = psi_a + c (psi_X - psi_a), c < 1, and the X-point rank-one term.
    with_second = .true.
    call sweep('FB-cut c=0.8 X point, edge jump, P3', 3, [32, 64], GS_FB_XPOINT, &
               .true., 3.0_dp)
    ! The separatrix itself (c = 1): the source must vanish at the edge there.
    case_fraction = 1.0_dp
    call sweep('separatrix X point, vanishing edge, P2', 2, [32, 64, 128], &
               GS_FB_XPOINT, .false., 2.0_dp)
    case_fraction = 0.8_dp
    with_second = .false.
    if (failures /= 0) error stop 'gs_fb_manufactured failures'
    print '(a)', 'gs_fb_manufactured (T3) oracles passed'

contains

    subroutine require(ok, what)
        logical, intent(in) :: ok
        character(len=*), intent(in) :: what
        if (.not. ok) then
            failures = failures + 1
            print '(2a)', 'FAIL: ', what
        end if
    end subroutine require

    function vacuum(r, z) result(psi)
        real(dp), intent(in) :: r, z
        real(dp) :: psi, g(1, 2)

        call gs_circuit_flux_vectors(coils, [r], [z], g)
        psi = g(1, 1)*currents(1) + g(1, 2)*currents(2)
    end function vacuum

    function exact_flux(r, z) result(psi)
        real(dp), intent(in) :: r, z
        real(dp) :: psi

        psi = vacuum(r, z) + bump(r, z, centre, semi, amplitude)
        if (with_second) psi = psi + bump(r, z, low_c, low_s, 0.9_dp)
    end function exact_flux

    function density(r, z) result(j)
        !! mu0 j_phi of psi*: the source of the bump.
        use fortfem_axisymmetric_gs_kernel, only: gs_mu0
        real(dp), intent(in) :: r, z
        real(dp) :: j

        j = bump_laplacian(r, z, centre, semi, amplitude)
        if (with_second) j = j + bump_laplacian(r, z, low_c, low_s, 0.9_dp)
        j = j/gs_mu0
    end function density

    subroutine sweep(label, order, levels, mode, jump, min_rate)
        character(len=*), intent(in) :: label
        integer, intent(in) :: order, levels(:), mode
        logical, intent(in) :: jump
        real(dp), intent(in) :: min_rate
        real(dp) :: errors(3, size(levels)), rate
        integer :: i, iterations

        do i = 1, size(levels)
          call solve_case(label, order, levels(i), mode, jump, errors(:, i), iterations)
        end do
        if (size(levels) > 1) then
            rate = log(errors(1, size(levels) - 1)/errors(1, size(levels)))/log(2.0_dp)
      print '(a,a,f6.2)', trim(label), ': observed nodal rate on the finest pair ', rate
            call require(rate > min_rate, trim(label)//': nodal convergence rate')
            call require(all(errors(1, 2:) < errors(1, :size(levels) - 1)), &
                         trim(label)//': error decreases with refinement')
        end if
    end subroutine sweep

    subroutine solve_case(label, order, n_outer, mode, jump, errors, iterations)
        !! errors = [relative max nodal error, relative psi_b error, relative
        !! residual of the seed-independent solution].
        character(len=*), intent(in) :: label
        integer, intent(in) :: order, n_outer, mode
        logical, intent(in) :: jump
        real(dp), intent(out) :: errors(3)
        integer, intent(out) :: iterations
        type(gs_state_t) :: state
        type(gs_profile_law_t) :: law
        type(gs_fb_config_t) :: cfg
        type(gs_fb_system_t) :: sys
        type(gs_fb_eval_t) :: ev
        type(gs_fb_history_t) :: history
        real(dp), allocatable :: exact(:), seed(:), f_load(:), vac(:)
        real(dp) :: scale, psib_exact, angle, g(1, 2)
        integer :: status, i, n, m, k

        call build_state(state, n_outer, order, status)
        call require(status == 0, 'mesh')
        n = size(state%points, 2)
        allocate (exact(n), seed(n), vac(n))
        do i = 1, n
            exact(i) = exact_flux(state%points(1, i), state%points(2, i))
            vac(i) = vacuum(state%points(1, i), state%points(2, i))
            seed(i) = vac(i) + bump(state%points(1, i), state%points(2, i), &
                                   centre + [0.10_dp, -0.08_dp], semi, 0.3_dp*amplitude)
            if (with_second) seed(i) = seed(i) + bump(state%points(1, i), &
                          state%points(2, i), low_c + [0.03_dp, 0.02_dp], low_s, 0.7_dp)
        end do
        call make_law(law, jump, 0.7_dp)
        cfg%mode = mode
        cfg%axis_sign = 1
        cfg%axis_guess = centre
        cfg%search_radius = 0.4_dp
        cfg%quadrature_tolerance = 1.0e-11_dp
        cfg%jacobian_tolerance = 1.0e-11_dp
        cfg%tolerance = 1.0e-10_dp
        if (mode == GS_FB_XPOINT) then
            cfg%xpoint_guess = [gp_r0, -0.62_dp]
            cfg%cut_fraction = case_fraction
            cfg%search_radius = 0.45_dp
            psib_exact = 0.0_dp
        else if (mode == GS_FB_LIMITER) then
            allocate (cfg%points(2, 1))
            cfg%points(:, 1) = [gp_r0 + 0.33_dp, 0.0_dp]
            psib_exact = exact_flux(cfg%points(1, 1), cfg%points(2, 1))
        else
            m = 8
            allocate (cfg%points(2, m))
            psib_exact = 0.0_dp
            do k = 1, m
                angle = 6.283185307179586_dp*(k - 0.5_dp)/m
                cfg%points(:, k) = [gp_r0 + 0.3_dp*cos(angle), 0.48_dp*sin(angle)]
              psib_exact = psib_exact + exact_flux(cfg%points(1, k), cfg%points(2, k))/m
            end do
        end if
        call gs_fb_setup(state, sys, law, cfg, coils, status)
        call require(status == GS_FB_OK, 'setup')
        call gs_coupling_source_load(state, density, f_load, status)
        call require(status == GS_COUPLING_OK, 'forcing load')
        ! extra = int f* v - b(Pi psi*).
        call gs_fb_evaluate(state, sys, exact, currents, ev, status, jacobian=.false.)
        call require(status == GS_FB_OK, 'evaluate at the exact flux')
        sys%extra_load = f_load - ev%source
        ! The X-point level is a property of the sampled flux itself.
        if (mode == GS_FB_XPOINT) psib_exact = ev%psi_bound
        errors(3) = norm2(ev%residual + ev%source - f_load + ev%source - ev%source) &
                    /max(norm2(f_load), tiny(1.0_dp))
        ! Residual of the exact flux in the manufactured problem = discretisation
        ! error of the coupled operator and the quadrature of f*.
        call gs_fb_evaluate(state, sys, exact, currents, ev, status, jacobian=.false.)
        errors(3) = norm2(ev%residual)/norm2(f_load)
        state%psi = seed
        call gs_fb_newton(state, sys, currents, status, history)
        iterations = history%iterations
        if (status /= GS_FB_OK) print '(a,i4)', '   Newton status ', status
        call require(status == GS_FB_OK, trim(label)//': Newton converged')
        scale = maxval(abs(exact))
        errors(1) = maxval(abs(state%psi - exact))/scale
      call gs_fb_evaluate(state, sys, state%psi, currents, ev, status, jacobian=.false.)
        errors(2) = abs(ev%psi_bound - psib_exact)/scale
    print '(a,a,i4,a,i2,a,3es10.2)', trim(label), ' edges ', n_outer, '  Newton its ', &
            iterations, '  nodal, psi_b, exact-residual ', errors
        print '(a,12es10.2)', '   residual history ', history%residual
        print '(a,12es10.2)', '   damping history  ', history%damping
        print '(a,12es10.2)', '   step history     ', history%step
        print '(a,12es10.2)', '   Ip history       ', history%ip
        call check_quadratic(history)
    end subroutine solve_case

    subroutine check_quadratic(history)
        !! Quadratic Newton convergence: once the relative residual is below
        !! 1e-2, each step squares it up to a modest constant, until the
        !! quadrature floor.
        type(gs_fb_history_t), intent(in) :: history
        integer :: k
        logical :: seen

        seen = .false.
        do k = 1, history%iterations
 if (history%residual(k - 1) < 1.0e-1_dp .and. history%residual(k - 1) > 1.0e-6_dp) then
                seen = .true.
               call require(history%residual(k) <= 30.0_dp*history%residual(k - 1)**2, &
                             'quadratic residual decay')
                print '(a,i2,a,es10.2)', '   iteration ', k, ' r_k / r_{k-1}^2 ', &
                    history%residual(k)/history%residual(k - 1)**2
            end if
        end do
        call require(seen, 'quadratic phase observed')
    end subroutine check_quadratic

end program test_gs_fb_manufactured
