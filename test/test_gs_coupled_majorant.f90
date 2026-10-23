module test_gs_coupled_majorant_callbacks
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortfem_axisymmetric_gs_kernel, only: gs_mu0
    use kin6d_grad_shafranov, only: gs_state_t, gs_cell_data
    use kin6d_grad_shafranov_leaves, only: gs_source
    use kin6d_gs_coils, only: gs_circuit_set_t, gs_circuit_flux_vectors
    use kin6d_gs_coupling, only: gs_coupling_t, gs_coupling_build, &
        gs_coupling_coil_load, gs_coupling_source_load, gs_coupled_solve, &
        gs_coupled_matrix, GS_COUPLING_OK
    use kin6d_gs_profile_law, only: gs_profile_law_t, gs_law_sources
    use kin6d_gs_free_boundary, only: gs_fb_config_t, gs_fb_system_t, gs_fb_eval_t, &
        gs_fb_history_t, gs_fb_setup, gs_fb_evaluate, gs_fb_newton, &
        gs_fb_jacobian_apply, GS_FB_OK, GS_FB_LIMITER
    use kin6d_gs_majorant, only: gs_coupled_estimator_t, gs_coupled_majorant_t, &
        gs_coupled_goal_t, gs_fb_majorant_t, gs_coupled_estimator_build, &
        gs_coupled_majorant, gs_coupled_goal, gs_coupled_ip_moments, &
        gs_fb_majorant, gs_fb_inf_sup, gs_fb_goal, gs_coil_current_bounds, GS_CM_OK
    use fortfem_triangle_duffy_quadrature, only: triangle_duffy_quadrature
    use gs_coupling_fixture, only: build_state, build_coils, gp_r0, gp_a, gp_b
    use gs_fb_fixture, only: make_law
    use kin6d_gs_certificate, only: gs_fb_certifier_t
    use kin6d_equilibrium_certificate, only: equilibrium_certificate_t, certify, &
        EQ_CLAIM_PRACTICAL

    use kin6d_gs_nonlinear, only: gs_forcing_function_t
    use kin6d_gs_coupling, only: gs_density_function_t
    use kin6d_gs_majorant, only: gs_fb_majorant_object
    implicit none
    private
    public :: run_test_gs_coupled_majorant
    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp), parameter :: currents(2) = [1.5e4_dp, -1.0e4_dp]
    real(dp), parameter :: centre(2) = [gp_r0 + 0.02_dp, 0.03_dp]
    real(dp), parameter :: semi(2) = [gp_a, gp_b]
    real(dp), parameter :: amplitude = 1.5_dp
    real(dp), parameter :: goal_point(2) = [gp_r0 + 0.25_dp, 0.35_dp]
    real(dp), parameter :: goal_radius = 0.12_dp
    real(dp), parameter :: limiter(2) = [gp_r0 + 0.33_dp, 0.0_dp]
    type :: run_test_gs_coupled_majorant_context_t
        type(gs_circuit_set_t), pointer :: coils => null()
        type(gs_profile_law_t), pointer :: law_exact => null()
        real(dp), pointer :: psi_bound_exact => null()
    end type

    type, extends(gs_forcing_function_t) :: run_test_gs_coupled_majorant_forcing_t
        type(run_test_gs_coupled_majorant_context_t), pointer :: data => null()
    contains
        procedure :: evaluate => run_test_gs_coupled_majorant_forcing_binding
    end type

    type, extends(gs_density_function_t) :: &
        run_test_gs_coupled_majorant_forcing_density_t
        type(run_test_gs_coupled_majorant_context_t), pointer :: data => null()
    contains
        procedure :: evaluate => run_test_gs_coupled_majorant_forcing_density_binding
    end type

contains
    subroutine run_test_gs_coupled_majorant()
    !! P10 (DESIGN section 4): coupled FEM-BEM majorant, inf-sup estimate and
    !! goal functionals on manufactured free-boundary problems with exact
    !! solutions psi* = sum_k I_k g_k + w (w a compact analytic bump).
    !!
    !! 1. Frozen (linear) coupled problem a(psi,v) + <S psi,v> = (f*,v) + <g,v>:
    !!    the majorant D + C_F B + N must exceed the exact energy error
    !!    |||psi* - psi_h||| (interior energy by quadrature of the analytic
    !!    gradient, exterior energy <S e, e> on the reference BEM) with an
    !!    efficiency index in [1, 3] that does not grow under refinement, P2 and
    !!    P3; the shape and Ip goal enclosures must contain the exact values.
    !! 2. Nonlinear free boundary (limiter mode, law vanishing at the edge)
    !!    with the pointwise forcing f* - f_law(psi*) 1{plasma*}: the practical
    !!    bound M0 / beta against the exact error, beta against a dense SVD of
    !!    J in the energy pairing on the coarsest mesh, and the dual-weighted
    !!    shape goal against the exact goal error.
    !! 3. Coil-current bounds: attained by the sign vector, never exceeded.
    implicit none

    type(gs_circuit_set_t) , target :: coils
    type(gs_profile_law_t) , target :: law_exact
    real(dp) , target :: psi_axis_exact, psi_bound_exact
    integer :: failures
    character(len=16) :: mode
    real(dp) :: quick_errors(2), quick_bounds(2)
        type(run_test_gs_coupled_majorant_context_t), target :: ctx
        type(run_test_gs_coupled_majorant_forcing_t), target :: forcing_operator
        type(run_test_gs_coupled_majorant_forcing_density_t), target :: &
            forcing_density_operator

    failures = 0
    call build_coils(coils)
    call get_command_argument(1, mode)
    if (trim(mode) == 'quick') then
        call linear_case('P2 quick', 2, 12, quick_errors(1), quick_bounds(1))
        call linear_case('P2 quick', 2, 24, quick_errors(2), quick_bounds(2))
        call require(all(quick_bounds >= quick_errors), &
            'RT bound contains energy error')
        call require(quick_bounds(2) < quick_bounds(1), &
            'RT bound decreases on refinement')
        call linear_case('P3 quick', 3, 12, quick_errors(1), quick_bounds(1))
        call require(quick_bounds(1) >= quick_errors(1), &
            'P3 bound contains energy error')
        call nonlinear_sweep('P2 limiter quick', 2, [24])
    else
        call linear_sweep('P2', 2, [24, 48, 96])
        call linear_sweep('P3', 3, [24, 48])
        call nonlinear_sweep('P2 limiter', 2, [24, 48, 96])
    end if
    call coil_bound_check()
    if (failures /= 0) error stop 'gs_coupled_majorant failures'
    print '(a)', 'gs_coupled_majorant (P10) checks passed'

contains
    subroutine bind_run_test_gs_coupled_majorant_context()
        ctx%coils => coils
        ctx%law_exact => law_exact
        ctx%psi_bound_exact => psi_bound_exact
        forcing_operator%data => ctx
        forcing_density_operator%data => ctx
    end subroutine bind_run_test_gs_coupled_majorant_context


    subroutine require(ok, what)
        logical, intent(in) :: ok
        character(len=*), intent(in) :: what
        if (.not. ok) then
            failures = failures + 1
            print '(2a)', 'FAIL: ', what
        end if
    end subroutine require

    pure subroutine bump_gradient(r, z, gr, gz)
        real(dp), intent(in) :: r, z
        real(dp), intent(out) :: gr, gz
        real(dp) :: u, w1
        gr = 0.0_dp
        gz = 0.0_dp
        u = ((r - centre(1))/semi(1))**2 + ((z - centre(2))/semi(2))**2
        if (u >= 1.0_dp) return
        w1 = -6.0_dp*amplitude*(1.0_dp - u)**5
        gr = w1*2.0_dp*(r - centre(1))/semi(1)**2
        gz = w1*2.0_dp*(z - centre(2))/semi(2)**2
    end subroutine bump_gradient

    function goal_weight(r, z) result(q)
        !! Mollifier 3/(pi rho^2) (1 - |x - p|^2/rho^2)^2, unit integral.
        real(dp), intent(in) :: r, z
        real(dp) :: q, u
        u = ((r - goal_point(1))**2 + (z - goal_point(2))**2)/goal_radius**2
        q = 0.0_dp
        if (u < 1.0_dp) q = 3.0_dp/(pi*goal_radius**2)*(1.0_dp - u)**2
    end function goal_weight

    function zero_weight(r, z) result(q)
        real(dp), intent(in) :: r, z
        real(dp) :: q
        q = 0.0_dp*(r + z)
    end function zero_weight

    function exact_goal() result(q)
        !! int q psi* over the goal disk, polar Gauss-Legendre 24 x 64.
        real(dp) :: q, rho, theta, xr(24), xw(24)
        integer :: i, j
        call gauss_nodes(24, xr, xw)
        q = 0.0_dp
        do i = 1, 24
            rho = goal_radius*xr(i)
            do j = 1, 64
                theta = 2.0_dp*pi*(j - 0.5_dp)/64
        call bind_run_test_gs_coupled_majorant_context()
                q = q + xw(i)*goal_radius*rho*(2.0_dp*pi/64) &
                    *goal_weight(goal_point(1) + rho*cos(theta), &
                                 goal_point(2) + rho*sin(theta)) &
                    *run_test_gs_coupled_majorant_exact_flux(ctx, goal_point(1) + &
                        rho*cos(theta), &
                                goal_point(2) + rho*sin(theta))
            end do
        end do
    end function exact_goal

    function absolute_current() result(total)
        !! int |f*| / mu0 dA, midpoint rule on 400 x 400 cells of the bump box.
        real(dp) :: total, dr, dz
        integer :: i, j
        dr = 2*semi(1)/400
        dz = 2*semi(2)/400
        total = 0.0_dp
        do j = 1, 400
            do i = 1, 400
        call bind_run_test_gs_coupled_majorant_context()
                total = total + abs(run_test_gs_coupled_majorant_bump_laplacian( &
                    centre(1) - semi(1) + (i - 0.5_dp)*dr, &
                    centre(2) - semi(2) + (j - 0.5_dp)*dz))
            end do
        end do
        total = total*dr*dz/gs_mu0
    end function absolute_current

    subroutine gauss_nodes(n, x, w)
        use fortnum_quadrature, only: gauss_legendre
        integer, intent(in) :: n
        real(dp), intent(out) :: x(n), w(n)
        call gauss_legendre(n, x, w)
        x = 0.5_dp*(x + 1.0_dp)
        w = 0.5_dp*w
    end subroutine gauss_nodes

    function energy_error(state, est, psi) result(err)
        !! |||psi* - psi_h|||: interior energy of the analytic gradient error
        !! (Duffy rule of degree 12) and the exterior energy on the reference BEM.
        type(gs_state_t), intent(in) :: state
        type(gs_coupled_estimator_t), intent(in) :: est
        real(dp), intent(in) :: psi(:)
        real(dp) :: err, interior, exterior, shape(10), grads(2, 10), hess(2, 2, 10)
        real(dp) :: pt(2), det, grad(2), br, bz
        real(dp), allocatable :: x(:), y(:), w(:), g(:, :), gr(:, :), gz(:, :), e(:)
        real(dp), allocatable :: pts(:, :, :)
        integer :: c, q, nd, stat, i
        call triangle_duffy_quadrature(12, x, y, w, stat)
        nd = size(state%triangles, 1)
        interior = 0.0_dp
        allocate (g(size(w), 2), gr(size(w), 2), gz(size(w), 2), pts(2, size(w), 1))
        do c = 1, size(state%triangles, 2)
            do q = 1, size(w)
                call gs_cell_data(state, c, [x(q), y(q)], pts(:, q, 1), det, &
                                  shape(:nd), grads(:, :nd), hess(:, :, :nd), stat)
            end do
            call gs_circuit_flux_vectors(coils, pts(1, :, 1), pts(2, :, 1), g, gr, gz)
            do q = 1, size(w)
                call gs_cell_data(state, c, [x(q), y(q)], pt, det, shape(:nd), &
                                  grads(:, :nd), hess(:, :, :nd), stat)
                grad = matmul(grads(:, :nd), psi(state%triangles(:, c)))
                call bump_gradient(pt(1), pt(2), br, bz)
                interior = interior + w(q)*det*((dot_product(gr(q, :), currents) &
                    + br - grad(1))**2 + (dot_product(gz(q, :), currents) + bz &
                    - grad(2))**2)/pt(1)
            end do
        end do
        allocate (e(est%nx))
        e = matmul(est%restriction, psi(est%coupling%trace_node))
        do i = 1, est%nx
        call bind_run_test_gs_coupled_majorant_context()
            e(i) = run_test_gs_coupled_majorant_exact_flux(ctx, est%curve%nodes(1, &
                i), est%curve%nodes(2, i)) - e(i)
        end do
        exterior = dot_product(e, matmul(est%s_ref, e))
        err = sqrt(interior + exterior)
    end function energy_error

    subroutine linear_sweep(label, order, levels)
        character(len=*), intent(in) :: label
        integer, intent(in) :: order, levels(:)
        real(dp) :: errors(size(levels)), bounds(size(levels)), eff(size(levels))
        real(dp) :: rate_e, rate_m
        integer :: i, nl
        nl = size(levels)
        do i = 1, nl
            call linear_case(label, order, levels(i), errors(i), bounds(i))
            eff(i) = bounds(i)/errors(i)
            call require(eff(i) >= 1.0_dp .and. eff(i) <= 3.0_dp, &
                         label//': linear efficiency index in [1, 3]')
        end do
        rate_e = log(errors(nl - 1)/errors(nl))/log(2.0_dp)
        rate_m = log(bounds(nl - 1)/bounds(nl))/log(2.0_dp)
        print '(a,a,2f7.3)', label, ' linear: energy-error rate, majorant rate ', &
            rate_e, rate_m
        ! Sharp: the majorant converges at least at the rate of the error and
        ! its efficiency index decreases towards 1 under refinement.
        call require(rate_m >= rate_e - 0.1_dp, label//': majorant rate >= error rate')
        call require(all(eff(2:) < eff(:nl - 1)), label//': efficiency decreasing')
    end subroutine linear_sweep

    subroutine linear_case(label, order, n_outer, error, majorant)
        character(len=*), intent(in) :: label
        integer, intent(in) :: order, n_outer
        real(dp), intent(out) :: error, majorant
        type(gs_state_t) :: state
        type(gs_coupling_t) :: coupling
        type(gs_coupled_estimator_t) :: est
        type(gs_coupled_majorant_t) :: res, shifted
        type(gs_coupled_goal_t) :: shape_goal, ip_goal
        real(dp), allocatable :: coil_load(:, :), f_load(:), rhs(:), psi(:), bl(:)
        real(dp) :: q_exact, ip_h
        integer :: status
        integer(8) :: t0, t1, t2, rate
        call build_state(state, n_outer, order, status)
        call require(status == 0, 'mesh')
        call gs_coupling_build(state, coupling, status)
        call require(status == GS_COUPLING_OK, 'coupling')
        call gs_coupling_coil_load(state, coupling, coils, coil_load, status)
        call gs_coupling_source_load(state, linear_density, f_load, status, n_gauss=10)
        rhs = f_load + matmul(coil_load, currents)
        call gs_coupled_solve(coupling, rhs, psi, status, relative_tolerance=1.0e-13_dp)
        call require(status == GS_COUPLING_OK, 'coupled solve')
        call system_clock(t0, rate)
        call gs_coupled_estimator_build(state, coupling, est, status, set=coils)
        if (status /= GS_CM_OK) print '(a,i3)', '   estimator status ', status
        call require(status == GS_CM_OK, 'estimator build')
        bl = matmul(est%coil_moments, currents)
        call system_clock(t1)
        call gs_coupled_majorant(state, est, psi, linear_source, bl, res, status)
        call system_clock(t2)
        print '(a,i7,a,f8.3,a,f8.3)', '   nodes ', size(psi), '  estimator build s ', &
            real(t1 - t0, dp)/rate, '  majorant s ', real(t2 - t1, dp)/rate
        call require(status == GS_CM_OK, 'coupled majorant')
        call require(res%exterior_status == 'unavailable-bem-remainder', &
            'measured interface norm does not bound the continuous BEM remainder')
        if (trim(mode) == 'quick') then
            call gs_coupled_majorant(state, est, psi, linear_source, &
                bl + 1.0_dp, shifted, status)
            call require(status == GS_CM_OK, 'perturbed boundary estimate')
            call require(shifted%neumann_defect > res%neumann_defect, &
                'arbitrary boundary-load mismatch increases interface defect')
            call require(abs(shifted%flux_defect - res%flux_defect) &
                <= 1.0e-12_dp*max(1.0_dp, res%flux_defect), &
                'boundary load does not change local interior RT recovery')
        end if
        error = energy_error(state, est, psi)
        majorant = res%majorant
        print '(a,a,i4,a,es10.3,a,es10.3,a,f6.3,a,3es10.2,a,es10.3,a,i4)', label, &
            ' edges', n_outer, ' error ', error, ' majorant ', majorant, ' I_eff ', &
            majorant/error, ' D,CfB,N ', res%flux_defect, res%balance_defect, &
            res%neumann_defect, ' C_F^2 ', est%cf2, ' PCG its', res%iterations
        ! Shape goal: mollified point value; Ip goal: Ampere on Gamma_c.
        call gs_coupled_goal(state, est, psi, linear_source, bl, goal_weight, &
                             0.0_dp*bl, majorant, shape_goal, status)
        call require(status == GS_CM_OK, 'shape goal')
        q_exact = exact_goal()
        print '(a,es12.4,a,es10.2,a,es10.2,a,es10.2)', '   shape goal Q* ', q_exact, &
            '  Q*-Q_h ', q_exact - shape_goal%value, '  rho ', shape_goal%correction, &
            '  bound ', shape_goal%bound
        call require(q_exact >= shape_goal%lower .and. q_exact <= shape_goal%upper, &
                     label//': shape goal enclosure')
        call gs_coupled_goal(state, est, psi, linear_source, bl, zero_weight, &
                             gs_coupled_ip_moments(est, gs_mu0), majorant, ip_goal, &
                             status)
        call require(status == GS_CM_OK, 'Ip goal')
        ! The bump carries no net current: Ip* = 0. The adjoint z = 1/mu0 lies
        ! in V_h, so the bound vanishes and Ip_h + rho is the volume current
        ! int f/mu0 up to the quadrature of f: tolerance 1e-7 of int |f|/mu0.
        ip_h = ip_goal%value - sum(bl)/gs_mu0
        print '(a,es10.2,a,es10.2,a,es10.2,a,es10.2)', '   Ip goal Ip_h ', ip_h, &
            '  rho ', ip_goal%correction, '  bound ', ip_goal%bound, &
            '  int|j| ', absolute_current()
        call require(abs(ip_h + ip_goal%correction) <= ip_goal%bound + 1.0e-7_dp &
                     *absolute_current(), label//': Ip goal enclosure')
    end subroutine linear_case

    ! ---------------------------------------------------------------
    ! Nonlinear manufactured free boundary.
    ! ---------------------------------------------------------------

    subroutine exact_levels()
        !! psi_a*: Newton on the analytic flux (central-difference gradient
        !! and Hessian), psi_b*: the limiter point.
        real(dp) :: x(2), h, g(2), hm(2, 2), det, dx(2)
        integer :: it
        x = centre
        h = 1.0e-4_dp
        do it = 1, 30
        call bind_run_test_gs_coupled_majorant_context()
            g(1) = (run_test_gs_coupled_majorant_exact_flux(ctx, x(1) + h, x(2)) - &
                run_test_gs_coupled_majorant_exact_flux(ctx, x(1) - h, x(2)))/(2*h)
        call bind_run_test_gs_coupled_majorant_context()
            g(2) = (run_test_gs_coupled_majorant_exact_flux(ctx, x(1), x(2) + h) - &
                run_test_gs_coupled_majorant_exact_flux(ctx, x(1), x(2) - h))/(2*h)
        call bind_run_test_gs_coupled_majorant_context()
            hm(1, 1) = (run_test_gs_coupled_majorant_exact_flux(ctx, x(1) + h, x(2)) &
                - 2*run_test_gs_coupled_majorant_exact_flux(ctx, x(1), x(2)) &
                        + run_test_gs_coupled_majorant_exact_flux(ctx, x(1) - h, &
                            x(2)))/h**2
        call bind_run_test_gs_coupled_majorant_context()
            hm(2, 2) = (run_test_gs_coupled_majorant_exact_flux(ctx, x(1), x(2) + h) &
                - 2*run_test_gs_coupled_majorant_exact_flux(ctx, x(1), x(2)) &
                        + run_test_gs_coupled_majorant_exact_flux(ctx, x(1), x(2) - &
                            h))/h**2
        call bind_run_test_gs_coupled_majorant_context()
            hm(1, 2) = (run_test_gs_coupled_majorant_exact_flux(ctx, x(1) + h, x(2) &
                + h) &
                        - run_test_gs_coupled_majorant_exact_flux(ctx, x(1) + h, &
                            x(2) - h) &
                        - run_test_gs_coupled_majorant_exact_flux(ctx, x(1) - h, &
                            x(2) + h) &
                        + run_test_gs_coupled_majorant_exact_flux(ctx, x(1) - h, &
                            x(2) - h))/(4*h*h)
            hm(2, 1) = hm(1, 2)
            det = hm(1, 1)*hm(2, 2) - hm(1, 2)**2
            dx = -[hm(2, 2)*g(1) - hm(1, 2)*g(2), hm(1, 1)*g(2) - hm(1, 2)*g(1)]/det
            x = x + dx
            if (norm2(dx) < 1.0e-11_dp) exit
        end do
        call bind_run_test_gs_coupled_majorant_context()
        psi_axis_exact = run_test_gs_coupled_majorant_exact_flux(ctx, x(1), x(2))
        call bind_run_test_gs_coupled_majorant_context()
        psi_bound_exact = run_test_gs_coupled_majorant_exact_flux(ctx, limiter(1), &
            limiter(2))
        call make_law(law_exact, .false., 0.7_dp)
        law_exact%axis = psi_axis_exact - psi_bound_exact
    end subroutine exact_levels

    subroutine nonlinear_sweep(label, order, levels)
        character(len=*), intent(in) :: label
        integer, intent(in) :: order, levels(:)
        real(dp) :: errors(size(levels)), bounds(size(levels))
        integer :: i
        call exact_levels()
        print '(a,2es14.6)', 'exact psi_a*, psi_b* ', psi_axis_exact, psi_bound_exact
        do i = 1, size(levels)
            call nonlinear_case(label, order, levels(i), i == 1, errors(i), bounds(i))
            call require(bounds(i) >= errors(i), label//': M0/beta above the error')
            call require(bounds(i) <= 4.0_dp*errors(i), label//': M0/beta efficiency')
        end do
    end subroutine nonlinear_sweep

    subroutine nonlinear_case(label, order, n_outer, dense, error, bound)
        character(len=*), intent(in) :: label
        integer, intent(in) :: order, n_outer
        logical, intent(in) :: dense
        real(dp), intent(out) :: error, bound
        type(gs_state_t), target :: state
        type(gs_profile_law_t) :: law
        type(gs_fb_config_t) :: cfg
        type(gs_fb_system_t), target :: sys
        type(gs_fb_eval_t) :: ev
        type(gs_fb_history_t) :: history
        type(gs_coupled_estimator_t), target :: est
        type(gs_fb_certifier_t) :: certifier
        type(equilibrium_certificate_t) :: certificate
        type(gs_fb_majorant_t) :: res
        type(gs_coupled_goal_t) :: goal
        real(dp), allocatable :: f_load(:), exact(:)
        real(dp) :: beta_dense, q_exact
        integer :: status, i, its
        call build_state(state, n_outer, order, status)
        call make_law(law, .false., 0.7_dp)
        cfg%mode = GS_FB_LIMITER
        cfg%axis_guess = centre
        cfg%search_radius = 0.4_dp
        cfg%quadrature_tolerance = 1.0e-11_dp
        cfg%jacobian_tolerance = 1.0e-11_dp
        allocate (cfg%points(2, 1))
        cfg%points(:, 1) = limiter
        call gs_fb_setup(state, sys, law, cfg, coils, status)
        call require(status == GS_FB_OK, 'fb setup')
        call bind_run_test_gs_coupled_majorant_context()
        call gs_coupling_source_load(state, forcing_density_operator, f_load, &
            status, n_gauss=12)
        sys%extra_load = f_load
        allocate (exact(size(state%points, 2)))
        do i = 1, size(exact)
        call bind_run_test_gs_coupled_majorant_context()
            exact(i) = run_test_gs_coupled_majorant_exact_flux(ctx, state%points(1, &
                i), state%points(2, i))
        end do
        state%psi = exact
        call gs_fb_newton(state, sys, currents, status, history)
        call require(status == GS_FB_OK, label//': Newton converged')
        call gs_coupled_estimator_build(state, sys%coupling, est, status, set=coils)
        if (status /= GS_CM_OK) print '(a,i3)', '   estimator status ', status
        call require(status == GS_CM_OK, 'fb estimator build')
        call bind_run_test_gs_coupled_majorant_context()
        call gs_fb_majorant_object(state, sys, est, currents, res, status, &
            forcing_operator)
        call require(status == GS_CM_OK, 'fb majorant')
        error = energy_error(state, est, state%psi)
        bound = res%bound
        print '(a,a,i4,a,es10.3,a,es10.3,a,es10.3,a,f7.4,a,i4,a,f6.3)', label, &
            ' edges', n_outer, ' error ', error, ' M0 ', res%frozen%majorant, &
            ' M0/beta ', bound, ' beta ', res%beta, ' its', res%beta_iterations, &
            ' I_eff ', bound/error
        if (dense) then
            call gs_fb_evaluate(state, sys, state%psi, currents, ev, status, &
                                jacobian=.true.)
            beta_dense = dense_inf_sup(state, sys, ev)
            print '(a,f9.5,a,f9.5)', '   beta (power iteration) ', res%beta, &
                '  dense SVD ', beta_dense
            call require(abs(res%beta - beta_dense) <= 1.0e-3_dp*beta_dense, &
                         label//': beta matches the dense SVD')
            ! D4: the psi free-boundary certifier reports the same M0 and
            ! beta, and the core labels M0/beta a practical estimate.
            certifier%state => state
            certifier%sys => sys
            certifier%est => est
            certifier%currents = currents
        call bind_run_test_gs_coupled_majorant_context()
            certifier%forcing_function => forcing_operator
            call certify(certifier, certificate, status)
            call require(status == 0 .and. certificate%lane_status == GS_CM_OK, &
                         label//': certifier status')
            call require(certificate%claim == EQ_CLAIM_PRACTICAL, &
                         label//': free-boundary claim is a practical estimate')
            call require(abs(certificate%bound - res%bound) <= 1.0e-10_dp*res%bound, &
                         label//': certificate bound equals M0/beta')
            call gs_fb_inf_sup(state, sys, ev, beta_dense, its, status)
            call require(status == GS_CM_OK, 'inf-sup direct call')
        end if
        call gs_fb_goal(state, sys, est, currents, goal_weight, &
                        0.0_dp*est%coil_moments(:, 1), res%frozen%majorant, goal, &
                        status)
        call require(status == GS_CM_OK, 'fb goal')
        q_exact = exact_goal()
        print '(a,es10.2,a,es10.2,a,es10.2)', '   shape goal Q*-Q_h ', &
            q_exact - goal%value, '  rho ', goal%correction, '  |||z_h||| M0 ', &
            goal%bound
        call require(abs(q_exact - goal%value) <= goal%bound, &
                     label//': nonlinear shape goal within |||z_h||| M0')
    end subroutine nonlinear_case

    function dense_inf_sup(state, sys, ev) result(beta)
        !! Smallest singular value of L^-1 J L^-T, A = L L^T (dense oracle).
        type(gs_state_t), intent(in) :: state
        type(gs_fb_system_t), intent(in) :: sys
        type(gs_fb_eval_t), intent(in) :: ev
        real(dp) :: beta
        real(dp), allocatable :: a(:, :), j(:, :), e(:), s(:), work(:)
        real(dp), allocatable :: u(:, :), vt(:, :)
        integer :: n, i, info, lwork
        interface
            subroutine dpotrf(uplo, n, a, lda, info)
                import :: dp
                character(len=1), intent(in) :: uplo
                integer, intent(in) :: n, lda
                real(dp), intent(inout) :: a(lda, *)
                integer, intent(out) :: info
            end subroutine dpotrf
            subroutine dtrsm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
                import :: dp
                character(len=1), intent(in) :: side, uplo, transa, diag
                integer, intent(in) :: m, n, lda, ldb
                real(dp), intent(in) :: alpha, a(lda, *)
                real(dp), intent(inout) :: b(ldb, *)
            end subroutine dtrsm
            subroutine dgesvd(jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, &
                              lwork, info)
                import :: dp
                character(len=1), intent(in) :: jobu, jobvt
                integer, intent(in) :: m, n, lda, ldu, ldvt, lwork
                real(dp), intent(inout) :: a(lda, *)
                real(dp), intent(out) :: s(*), u(ldu, *), vt(ldvt, *), work(*)
                integer, intent(out) :: info
            end subroutine dgesvd
        end interface
        n = sys%n
        call gs_coupled_matrix(sys%coupling, a)
        allocate (j(n, n), e(n))
        do i = 1, n
            e = 0.0_dp
            e(i) = 1.0_dp
            call gs_fb_jacobian_apply(state, sys, ev, e, j(:, i))
        end do
        call dpotrf('L', n, a, n, info)
        call require(info == 0, 'dense Cholesky of A')
        call dtrsm('L', 'L', 'N', 'N', n, n, 1.0_dp, a, n, j, n)
        call dtrsm('R', 'L', 'T', 'N', n, n, 1.0_dp, a, n, j, n)
        allocate (s(n), u(1, 1), vt(1, 1), work(1))
        call dgesvd('N', 'N', n, n, j, n, s, u, 1, vt, 1, work, -1, info)
        lwork = int(work(1))
        deallocate (work)
        allocate (work(lwork))
        call dgesvd('N', 'N', n, n, j, n, s, u, 1, vt, 1, work, lwork, info)
        call require(info == 0, 'dense SVD')
        beta = s(n)
    end function dense_inf_sup

    subroutine coil_bound_check()
        real(dp) :: d(7, 3), sens(3, 7), bounds(3), rows(7), di(3), dr(7)
        integer :: status, k, trial
        logical :: free(3)
        do k = 1, 7
            d(k, :) = [sin(1.3_dp*k), cos(0.7_dp*k) + 0.2_dp, 0.5_dp*sin(2.1_dp*k + 1)]
            rows(k) = 0.1_dp + 0.05_dp*k
        end do
        free = [.true., .true., .false.]
        call gs_coil_current_bounds(d, 0.3_dp, [1.0_dp, 2.0_dp, 1.0_dp], free, rows, &
                                    bounds, status, sens)
        call require(status == GS_CM_OK, 'coil bounds status')
        call require(bounds(3) == 0.0_dp, 'fixed circuit has no bound')
        do k = 1, 2
            dr = rows*sign(1.0_dp, sens(k, :))
            di = matmul(sens, dr)
            call require(abs(abs(di(k)) - bounds(k)) <= 1.0e-12_dp*bounds(k), &
                         'coil bound attained')
        end do
        do trial = 1, 50
            dr = rows*sin(3.7_dp*trial + [(1.1_dp*k, k=1, 7)])
            di = matmul(sens, dr)
            call require(all(abs(di) <= bounds*(1 + 1.0e-12_dp)), 'coil bound holds')
        end do
        print '(a,3es11.3)', 'coil current bounds ', bounds
    end subroutine coil_bound_check


    end subroutine run_test_gs_coupled_majorant

    function run_test_gs_coupled_majorant_exact_flux(ctx, r, z) result(psi)
        type(run_test_gs_coupled_majorant_context_t), intent(in) :: ctx
        real(dp), intent(in) :: r, z
        real(dp) :: psi, g(1, 2)
        call gs_circuit_flux_vectors(ctx%coils, [r], [z], g)
        psi = dot_product(g(1, :), currents) + run_test_gs_coupled_majorant_bump(r, z)
    end function run_test_gs_coupled_majorant_exact_flux

    pure function run_test_gs_coupled_majorant_bump(r, z) result(w)
        !! amplitude (1 - u)^6, u the elliptic radius squared (C^5).
        real(dp), intent(in) :: r, z
        real(dp) :: w, u
        u = ((r - centre(1))/semi(1))**2 + ((z - centre(2))/semi(2))**2
        w = 0.0_dp
        if (u < 1.0_dp) w = amplitude*(1.0_dp - u)**6
    end function run_test_gs_coupled_majorant_bump

    pure function run_test_gs_coupled_majorant_bump_laplacian(r, z) result(f)
        !! f = -d_R(w_R/R) - d_Z(w_Z/R): the source the bump alone carries.
        real(dp), intent(in) :: r, z
        real(dp) :: f, u, ur, uz, w1, w2
        f = 0.0_dp
        u = ((r - centre(1))/semi(1))**2 + ((z - centre(2))/semi(2))**2
        if (u >= 1.0_dp) return
        ur = 2.0_dp*(r - centre(1))/semi(1)**2
        uz = 2.0_dp*(z - centre(2))/semi(2)**2
        w1 = -6.0_dp*amplitude*(1.0_dp - u)**5
        w2 = 30.0_dp*amplitude*(1.0_dp - u)**4
        f = -((w2*ur*ur + w1*2.0_dp/semi(1)**2)/r - w1*ur/r**2) &
            - (w2*uz*uz + w1*2.0_dp/semi(2)**2)/r
    end function run_test_gs_coupled_majorant_bump_laplacian

    function run_test_gs_coupled_majorant_forcing(ctx, r, z) result(f)
        type(run_test_gs_coupled_majorant_context_t), intent(in) :: ctx
        !! f* - f_law(psi*) 1{plasma*}: psi* solves the nonlinear problem.
        real(dp), intent(in) :: r, z
        real(dp) :: f, p, pp, ff, s
        f = run_test_gs_coupled_majorant_bump_laplacian(r, z)
        if (((r - centre(1))/semi(1))**2 + ((z - centre(2))/semi(2))**2 >= 1) return
        p = run_test_gs_coupled_majorant_exact_flux(ctx, r, z)
        if (p <= ctx%psi_bound_exact) return
        call gs_law_sources(ctx%law_exact, p - ctx%psi_bound_exact, pp, ff)
        call gs_source(r, gs_mu0, pp, ff, s)
        f = f - s
    end function run_test_gs_coupled_majorant_forcing

    function run_test_gs_coupled_majorant_forcing_binding(self, r, z) &
        result(result_value)
        class(run_test_gs_coupled_majorant_forcing_t), intent(in) :: self
        real(dp), intent(in) :: r, z
        real(dp) :: result_value, p, pp, ff, s
        result_value = run_test_gs_coupled_majorant_forcing(self%data, r, z)
    end function

    function run_test_gs_coupled_majorant_forcing_density(ctx, r, z) result(j)
        type(run_test_gs_coupled_majorant_context_t), intent(in) :: ctx
        real(dp), intent(in) :: r, z
        real(dp) :: j
        j = run_test_gs_coupled_majorant_forcing(ctx, r, z)/gs_mu0
    end function run_test_gs_coupled_majorant_forcing_density

    function run_test_gs_coupled_majorant_forcing_density_binding(self, r, z) &
        result(result_value)
        class(run_test_gs_coupled_majorant_forcing_density_t), intent(in) :: self
        real(dp), intent(in) :: r, z
        real(dp) :: result_value
        result_value = run_test_gs_coupled_majorant_forcing_density(self%data, r, z)
    end function
    function linear_source(cell, r, z, psi) result(f)
        integer, intent(in) :: cell
        real(dp), intent(in) :: r, z, psi
        real(dp) :: f
        f = run_test_gs_coupled_majorant_bump_laplacian(r, z) + 0.0_dp*(cell + psi)
    end function linear_source

    function linear_density(r, z) result(j)
        real(dp), intent(in) :: r, z
        real(dp) :: j
        j = run_test_gs_coupled_majorant_bump_laplacian(r, z)/gs_mu0
    end function linear_density

end module test_gs_coupled_majorant_callbacks

program test_gs_coupled_majorant
    use test_gs_coupled_majorant_callbacks, only: run_test_gs_coupled_majorant
    implicit none
    call run_test_gs_coupled_majorant()
end program test_gs_coupled_majorant
