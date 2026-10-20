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
    use, intrinsic :: iso_fortran_env, only: dp => real64, output_unit
    use kin6d_grad_shafranov, only: gs_state_t
    use kin6d_gs_coils, only: gs_circuit_set_t, gs_circuit_flux_vectors
    use kin6d_gs_coupling, only: gs_coupling_source_load, GS_COUPLING_OK
    use kin6d_gs_profile_law, only: gs_profile_law_t
    use kin6d_gs_free_boundary, only: gs_fb_config_t, gs_fb_system_t, gs_fb_eval_t, &
                 gs_fb_history_t, gs_fb_setup, gs_fb_evaluate, gs_fb_newton, GS_FB_OK, &
                                      GS_FB_LIMITER, GS_FB_TARGET_MEAN, GS_FB_XPOINT
    use gs_coupling_fixture, only: build_state, build_coils, gp_r0, gp_a, gp_b
    use gs_fb_fixture, only: make_law, bump, bump_laplacian
    use kin6d_equilibrium_harness, only: equilibrium_case_t, &
        equilibrium_resolution_t, equilibrium_start_t, equilibrium_run_report_t, &
        equilibrium_start, run_equilibrium_case, write_equilibrium_report
    implicit none
    type :: psi_case_state_t
        type(gs_state_t) :: gs
        type(gs_fb_system_t) :: system
        type(gs_fb_history_t) :: history
        type(gs_fb_eval_t) :: evaluation
        real(dp), allocatable :: exact(:), forcing(:)
        real(dp) :: scale, seed_error
    end type
    type(equilibrium_start_t) :: selected_start
    character(16) :: command, start_name
    integer :: input_status

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
    start_name = 'legacy'
    if (command_argument_count() > 0) then
        call get_command_argument(1, command)
        if (command /= 'offset') error stop 'expected offset start-name'
        call get_command_argument(2, start_name)
    end if
    call equilibrium_start(trim(start_name), selected_start, input_status)
    if (input_status /= 0) error stop 'unknown D8 start'
    if (start_name /= 'legacy') then
        call sweep('D8 limiter offset P2', 2, [24, 48, 96], &
            GS_FB_LIMITER, .true., 2.0_dp)
    else
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
    end if
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
        real(dp) :: scale, psib_exact, angle, seed_semi(2)
        type(psi_case_state_t) :: result
        type(equilibrium_run_report_t) :: verification
        character(32) :: case_name
        integer :: status, i, n, m, k

        call build_state(state, n_outer, order, status)
        call require(status == 0, 'mesh')
        n = size(state%points, 2)
        allocate (exact(n), seed(n), vac(n))
        seed_semi = semi
        seed_semi(2) = semi(2)*(1 + selected_start%elongation)
        do i = 1, n
            exact(i) = exact_flux(state%points(1, i), state%points(2, i))
            vac(i) = vacuum(state%points(1, i), state%points(2, i))
            seed(i) = vac(i) + bump(state%points(1, i), state%points(2, i), &
                                   centre + selected_start%displacement, seed_semi, &
                                   selected_start%amplitude*amplitude)
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
        scale = maxval(abs(exact))
        result%gs = state
        result%system = sys
        result%exact = exact
        result%forcing = f_load
        result%scale = scale
        result%seed_error = maxval(abs(seed - exact))/scale
        select case (mode)
        case (GS_FB_LIMITER)
            case_name = 'limiter-'
        case (GS_FB_TARGET_MEAN)
            case_name = 'targets-'
        case (GS_FB_XPOINT)
            case_name = 'xpoint-'
            if (case_fraction >= 1.0_dp) case_name = 'separatrix-'
        end select
        if (jump) then
            case_name = trim(case_name)//'jump'
        else
            case_name = trim(case_name)//'smooth'
        end if
        call run_equilibrium_case(trim(case_name), 'psi', &
            equilibrium_resolution_t(order, n_outer), trim(start_name), &
            run_psi_case, result, verification)
        call write_equilibrium_report(output_unit, verification)
        call require(verification%status == 0, trim(label)//': D8 solve/report')
        if (verification%status /= 0) error stop 'D8 psi case'
        state = result%gs
        history = result%history
        ev = result%evaluation
        iterations = verification%iterations
        errors(1) = verification%oracle_error
        if (start_name /= 'legacy' .and. n_outer == 96) then
            call require(errors(1) < 0.03_dp, 'offset return within 3 percent')
            if (start_name /= 'reference') &
                call require(errors(1) < 0.1_dp*verification%seed_error, &
                             'offset corrected beyond seed stationarity')
        end if
        errors(2) = abs(ev%psi_bound - psib_exact)/scale
    print '(a,a,i4,a,i2,a,3es10.2)', trim(label), ' edges ', n_outer, '  Newton its ', &
            iterations, '  nodal, psi_b, exact-residual ', errors
        print '(a,12es10.2)', '   residual history ', history%residual
        print '(a,12es10.2)', '   damping history  ', history%damping
        print '(a,12es10.2)', '   step history     ', history%step
        print '(a,12es10.2)', '   Ip history       ', history%ip
        if (start_name == 'legacy') call check_quadratic(history)
    end subroutine solve_case

    subroutine run_psi_case(spec, resolution, start, state, verification)
        type(equilibrium_case_t), intent(in) :: spec
        type(equilibrium_resolution_t), intent(in) :: resolution
        type(equilibrium_start_t), intent(in) :: start
        class(*), intent(inout) :: state
        type(equilibrium_run_report_t), intent(inout) :: verification
        if (spec%lane /= 'psi') return
        if (resolution%degree < 2 .or. start%name /= start_name) return
        select type (state)
        type is (psi_case_state_t)
            call gs_fb_newton(state%gs, state%system, currents, &
                verification%status, state%history)
            if (verification%status /= GS_FB_OK) return
            verification%iterations = state%history%iterations
            verification%oracle_error = &
                maxval(abs(state%gs%psi - state%exact))/state%scale
            verification%seed_error = state%seed_error
            call gs_fb_evaluate(state%gs, state%system, state%gs%psi, &
                currents, state%evaluation, verification%status, jacobian=.false.)
            if (verification%status /= GS_FB_OK) return
            verification%residual = norm2(state%evaluation%residual) &
                /max(norm2(state%forcing), tiny(1.0_dp))
            verification%residual_norm = 'forced discrete l2 / l2(manufactured load)'
        end select
    end subroutine

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
