program test_gs_fb_refit
    !! P8 (DESIGN 1.2, step 2): the outer refit loop of the free-boundary solver.
    !!
    !! Part 1, boundary faces. The flux equals the boundary level exactly on
    !! the nodes of the Gamma_p chain (limiter point on the chain), so every
    !! chain edge is a face of the plasma boundary. The evaluation must not
    !! fail, must assign the cells by the sign at their centroid, must be
    !! continuous against a sliver (chain nodes moved by rounding-level and by
    !! 1e-9 relative amounts) and must carry the face interface term in the
    !! Jacobian (one-sided finite differences; the term is visible).
    !!
    !! Part 2, T3 manufactured case (limiter, edge jump) solved on the fixed
    !! mesh and with the refit loop; the nodal error orders are compared. Same
    !! manufactured solution as test_gs_fb_manufactured: the nodal load `extra`
    !! is rebuilt after every mesh change.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_grad_shafranov, only: gs_state_t, initialize_gs, gs_element_jet
    use kin6d_gs_coils, only: gs_circuit_set_t, gs_circuit_flux_vectors
  use kin6d_gs_coupling, only: gs_coupling_source_load, gs_coupled_apply, GS_COUPLING_OK
    use kin6d_gs_profile_law, only: gs_profile_law_t
    use kin6d_gs_free_boundary, only: gs_fb_config_t, gs_fb_system_t, gs_fb_eval_t, &
             gs_fb_refit_config_t, gs_fb_refit_report_t, gs_fb_history_t, gs_fb_setup, &
  gs_fb_evaluate, gs_fb_jacobian_apply, gs_fb_newton, gs_fb_refit_solve, gs_fb_locate, &
                                      gs_fb_trace_lcfs, &
                                      GS_FB_OK, GS_FB_LIMITER, GS_FB_REFIT_NOT_CONVERGED
    use gs_coupling_fixture, only: build_state, build_coils, outer_curve, pi, gp_r0, &
                                   gp_a, gp_b, gc_a, gc_b
    use gs_fb_fixture, only: make_law, bump, bump_laplacian, smooth_direction
    implicit none

    real(dp), parameter :: currents(2) = [1.5e4_dp, -1.0e4_dp]
    real(dp), parameter :: centre(2) = [gp_r0 + 0.02_dp, 0.03_dp]
    real(dp), parameter :: semi(2) = [gp_a, gp_b]
    real(dp), parameter :: amplitude = 1.5_dp
    type(gs_circuit_set_t) :: coils
    real(dp) :: law_scale = 5.0_dp, seed_amplitude = 1.5_dp
    integer :: failures
    character(len=16) :: part, arg

    failures = 0
    call build_coils(coils)
    part = 'all'
    if (command_argument_count() > 0) call get_command_argument(1, part)
    if (command_argument_count() > 2) then
        call get_command_argument(2, arg)
        read (arg, *) law_scale
        call get_command_argument(3, arg)
        read (arg, *) seed_amplitude
    end if
    if (part == 'all' .or. part == 'faces') then
        call face_case('P2', 2)
        call face_case('P3', 3)
    end if
    if (part == 'all' .or. part == 'sweep') then
        call refit_sweep('P2', 2, [24, 48, 96], 2.5_dp)
        call refit_sweep('P3', 3, [24, 48], 3.0_dp)
    end if
    if (part == 'all' .or. part == 'study') call reference_study()
    if (failures /= 0) error stop 'gs_fb_refit failures'
    print '(a)', 'gs_fb_refit (P8) oracles passed'

contains

    subroutine require(ok, what)
        logical, intent(in) :: ok
        character(len=*), intent(in) :: what
        if (.not. ok) then
            failures = failures + 1
            print '(2a)', 'FAIL: ', what
        end if
    end subroutine require

    ! ------------------------------------------------------------------
    ! Part 1: faces on the boundary.
    ! ------------------------------------------------------------------

    subroutine face_case(label, order)
        character(len=*), intent(in) :: label
        integer, intent(in) :: order
        integer, parameter :: n_outer = 36, nc = n_outer/2
        type(gs_state_t) :: state
        type(gs_profile_law_t) :: law
        type(gs_fb_config_t) :: cfg
        type(gs_fb_system_t) :: sys
        type(gs_fb_eval_t) :: ev, ev_a, ev_b
        real(dp), allocatable :: psi(:), delta(:), jd(:), none(:), ad(:), db(:)
        real(dp), allocatable :: forward(:), backward(:)
        real(dp) :: chain(2, nc), level, rho2, ip0, shift, plus, minus, err_one
        real(dp) :: err_no_face, scale
        real(dp), parameter :: shifts(3) = [1.0e-13_dp, 1.0e-9_dp, 1.0e-7_dp]
        logical, allocatable :: on_chain(:)
        integer :: status, i, n, j, step

        call build_state(state, n_outer, order, status)
        call require(status == 0, label//': mesh')
        do j = 1, nc
            chain(:, j) = [gp_r0 + gp_a*cos(2.0_dp*pi*(j - 1)/nc), &
                           gp_b*sin(2.0_dp*pi*(j - 1)/nc)]
        end do
        call chain_nodes(state, chain, on_chain)
        call make_law(law, .true., 1.0_dp)
        cfg%mode = GS_FB_LIMITER
        cfg%axis_sign = 1
        cfg%axis_guess = [gp_r0, 0.0_dp]
        cfg%search_radius = 0.35_dp
        cfg%quadrature_tolerance = 1.0e-12_dp
        cfg%jacobian_tolerance = 1.0e-12_dp
        allocate (cfg%points(2, 2))
        cfg%points(:, 1) = chain(:, 1)
        cfg%points(:, 2) = chain(:, nc/2 + 1)
        call gs_fb_setup(state, sys, law, cfg, status=status)
        call require(status == GS_FB_OK, label//': setup')
        n = sys%n
        allocate (psi(n), none(0))
        level = 0.3_dp
        scale = 1.0_dp
        do i = 1, n
            rho2 = ((state%points(1, i) - gp_r0)/gp_a)**2 + (state%points(2, i)/gp_b)**2
            psi(i) = level + scale*(1.0_dp - rho2)
            if (on_chain(i)) psi(i) = level
        end do
        call gs_fb_evaluate(state, sys, psi, none, ev, status)
        call require(status == GS_FB_OK, label//': evaluate on faces')
        if (status /= GS_FB_OK) return
        print '(a,a,i4,a,es12.4)', label, ': boundary faces ', ev%face_cells, &
            '  Ip ', ev%plasma_current
        call require(ev%face_cells == nc, label//': one plasma face per chain edge')
        ip0 = ev%plasma_current

        ! Continuity against rounding noise and slivers.
        do step = 1, 3
            shift = shifts(step)
            call gs_fb_evaluate(state, sys, merge(psi + shift, psi, on_chain), none, &
                                ev_a, status)
            call require(status == GS_FB_OK, label//': evaluate, chain nodes raised')
            call gs_fb_evaluate(state, sys, merge(psi - shift, psi, on_chain), none, &
                                ev_b, status)
            call require(status == GS_FB_OK, label//': evaluate, chain nodes lowered')
            plus = abs(ev_a%plasma_current - ip0)/ip0
            minus = abs(ev_b%plasma_current - ip0)/ip0
            print '(a,es9.1,a,2es10.2,a,2i5)', '   chain shift ', shift, &
                '  relative change of Ip ', plus, minus, '  face cells ', &
                ev_a%face_cells, ev_b%face_cells
            call require(max(plus, minus) < 50.0_dp*shift + 1.0e-12_dp, &
                         label//': source continuous across the face')
        end do

        ! Jacobian of the source: J d = A d - db. One-sided differences of the
        ! source into the plasma and into the vacuum side of the faces.
        call smooth_direction(state, delta)
        allocate (jd(n), ad(n))
        call gs_fb_jacobian_apply(state, sys, ev, delta, jd)
        call gs_coupled_apply(sys%coupling, delta, ad)
        shift = 1.0e-6_dp
        call gs_fb_evaluate(state, sys, psi + shift*delta, none, ev_a, status, &
                            jacobian=.false.)
        call require(status == GS_FB_OK, label//': evaluate forward')
        call gs_fb_evaluate(state, sys, psi - shift*delta, none, ev_b, status, &
                            jacobian=.false.)
        call require(status == GS_FB_OK, label//': evaluate backward')
        db = ad - jd
        forward = (ev_a%source - ev%source)/shift
        backward = (ev%source - ev_b%source)/shift
        err_one = min(norm2(db - forward), norm2(db - backward))/norm2(db)
        call gs_fb_jacobian_apply(state, sys, strip_faces(state, ev), delta, jd)
        db = ad - jd
        err_no_face = min(norm2(db - forward), norm2(db - backward))/norm2(db)
        print '(a,es10.2,a,es10.2)', '   source derivative vs one-sided differences ', &
            err_one, '  without the face term ', err_no_face
        call require(err_one < 0.05_dp, label//': face Jacobian one-sided difference')
        call require(err_no_face > 3.0_dp*err_one, label//': the face term matters')
    end subroutine face_case

    function strip_faces(state, ev) result(control)
        !! A copy of the evaluation whose interface mass (only the faces carry
        !! it in these states: every cut is a face) is removed from J.
        type(gs_state_t), intent(in) :: state
        type(gs_fb_eval_t), intent(in) :: ev
        type(gs_fb_eval_t) :: control
        integer :: e, k

        control = ev
        do e = 1, size(state%triangles, 2)
            do k = 1, size(state%triangles, 1)
                control%ub(state%triangles(k, e)) = control%ub(state%triangles(k, e)) &
                                                    + sum(ev%interface_mass(k, :, e))
            end do
        end do
        control%mass = ev%mass - ev%interface_mass
    end function strip_faces

    subroutine chain_nodes(state, chain, on_chain)
        !! Nodes of the chain edges: chain vertices and the edge nodes of the
        !! interior edges that join two consecutive chain vertices.
        type(gs_state_t), intent(in) :: state
        real(dp), intent(in) :: chain(:, :)
        logical, allocatable, intent(out) :: on_chain(:)
        integer, allocatable :: index(:)
        integer :: e, k, a, b, i, j, order, nc, node

        nc = size(chain, 2)
        order = merge(2, 3, size(state%triangles, 1) == 6)
        allocate (on_chain(size(state%points, 2)), index(size(state%points, 2)))
        on_chain = .false.
        index = 0
        do i = 1, size(state%points, 2)
            do j = 1, nc
                if (maxval(abs(state%points(:, i) - chain(:, j))) < 1.0e-12_dp) then
                    index(i) = j
                    on_chain(i) = .true.
                end if
            end do
        end do
        do e = 1, size(state%triangles, 2)
            do k = 1, 3
                if (state%neighbors(k, e) <= 0) cycle
                a = index(state%triangles(k, e))
                b = index(state%triangles(mod(k, 3) + 1, e))
                if (a == 0 .or. b == 0) cycle
                if (mod(a, nc) + 1 /= b .and. mod(b, nc) + 1 /= a) cycle
                do j = 1, order - 1
                    node = state%triangles(3 + (k - 1)*(order - 1) + j, e)
                    on_chain(node) = .true.
                end do
            end do
        end do
    end subroutine chain_nodes

    ! ------------------------------------------------------------------
    ! Part 2: T3 with and without the refit.
    ! ------------------------------------------------------------------

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
    end function exact_flux

    function density(r, z) result(j)
        use fortfem_axisymmetric_gs_kernel, only: gs_mu0
        real(dp), intent(in) :: r, z
        real(dp) :: j

        j = bump_laplacian(r, z, centre, semi, amplitude)/gs_mu0
    end function density

    subroutine manufactured_extra(state, sys)
        !! extra = int f* v - b(Pi psi*) on the current mesh.
        type(gs_state_t), intent(inout) :: state
        type(gs_fb_system_t), intent(inout) :: sys
        type(gs_fb_eval_t) :: ev
        real(dp), allocatable :: exact(:), f_load(:)
        integer :: i, status

        allocate (exact(size(state%points, 2)))
        do i = 1, size(exact)
            exact(i) = exact_flux(state%points(1, i), state%points(2, i))
        end do
        call gs_coupling_source_load(state, density, f_load, status)
        call require(status == GS_COUPLING_OK, 'forcing load')
        sys%extra_load = 0.0_dp
        call gs_fb_evaluate(state, sys, exact, currents, ev, status, jacobian=.false.)
        call require(status == GS_FB_OK, 'evaluate at the exact flux')
        sys%extra_load = f_load - ev%source
    end subroutine manufactured_extra

    subroutine smooth_extra(state, sys)
        !! extra = int f* v: the smooth load of the bump flux, no cancellation of
        !! the source jump.
        type(gs_state_t), intent(inout) :: state
        type(gs_fb_system_t), intent(inout) :: sys
        real(dp), allocatable :: f_load(:)
        integer :: status

        call gs_coupling_source_load(state, density, f_load, status)
        call require(status == GS_COUPLING_OK, 'forcing load')
        sys%extra_load = f_load
    end subroutine smooth_extra

    subroutine refit_sweep(label, order, levels, min_rate)
        character(len=*), intent(in) :: label
        integer, intent(in) :: order, levels(:)
        real(dp), intent(in) :: min_rate
        real(dp) :: fixed(size(levels)), refit(size(levels)), bound(size(levels))
        real(dp) :: bound_fixed(size(levels)), rate_fixed, rate_refit, mv
        integer :: i, passes(size(levels)), n

        n = size(levels)
        do i = 1, n
            call solve_case(order, levels(i), .false., fixed(i), bound_fixed(i), &
                            passes(i), mv)
            call solve_case(order, levels(i), .true., refit(i), bound(i), passes(i), mv)
            print '(a,a,i4,a,2es11.3,a,2es11.3,a,i2,a,es9.2)', label, ' edges ', &
                levels(i), '  nodal error fixed/refit ', fixed(i), refit(i), &
                '  psi_b error fixed/refit ', bound_fixed(i), bound(i), '  refits ', &
                passes(i), '  final chain distance/h ', mv
        end do
        rate_fixed = log(fixed(n - 1)/fixed(n))/log(2.0_dp)
        rate_refit = log(refit(n - 1)/refit(n))/log(2.0_dp)
        print '(a,a,2f7.2)', label, ': finest-pair nodal order fixed, refit ', &
            rate_fixed, rate_refit
        call require(rate_refit > min_rate, label//': refit nodal convergence order')
        call require(all(refit(2:) < refit(:n - 1)), label//': refit error decreases')
        call require(refit(n) <= 1.5_dp*fixed(n), label//': refit not worse than fixed')
    end subroutine refit_sweep

    subroutine solve_case(order, n_outer, refit, error, bound_error, refits, moved)
        integer, intent(in) :: order, n_outer
        logical, intent(in) :: refit
        real(dp), intent(out) :: error, bound_error, moved
        integer, intent(out) :: refits
        type(gs_state_t) :: state
        type(gs_profile_law_t) :: law
        type(gs_fb_config_t) :: cfg
        type(gs_fb_system_t) :: sys
        type(gs_fb_eval_t) :: ev
        type(gs_fb_refit_config_t) :: rcfg
        type(gs_fb_refit_report_t) :: report
        real(dp), allocatable :: contour(:, :), chain(:, :), seed(:), exact(:)
        real(dp) :: perimeter, h, psib_exact, scale
        integer :: j, nc, status, i

        nc = n_outer/2
        allocate (contour(2, n_outer), chain(2, nc))
        do j = 1, n_outer
            call outer_curve(real(j - 1, dp)/n_outer, contour(:, j))
        end do
        do j = 1, nc
            chain(:, j) = [gp_r0 + gp_a*cos(2.0_dp*pi*(j - 1)/nc), &
                           gp_b*sin(2.0_dp*pi*(j - 1)/nc)]
        end do
        perimeter = pi*(3.0_dp*(gc_a + gc_b) &
                        - sqrt((3.0_dp*gc_a + gc_b)*(gc_a + 3.0_dp*gc_b)))
        h = perimeter/n_outer
        call initialize_gs(state, contour, 0.4_dp*h*h, status, element_order=order, &
                           curved=.true., curve=outer_curve, interior_chain=chain)
        call require(status == 0, 'mesh')
        allocate (seed(size(state%points, 2)))
        do i = 1, size(seed)
            seed(i) = vacuum(state%points(1, i), state%points(2, i)) &
                      + bump(state%points(1, i), state%points(2, i), &
                             centre + [0.10_dp, -0.08_dp], semi, 0.3_dp*amplitude)
        end do
        call make_law(law, .true., 0.7_dp)
        cfg%mode = GS_FB_LIMITER
        cfg%axis_sign = 1
        cfg%axis_guess = centre
        cfg%search_radius = 0.4_dp
        cfg%quadrature_tolerance = 1.0e-11_dp
        cfg%jacobian_tolerance = 1.0e-11_dp
        cfg%tolerance = 1.0e-10_dp
        allocate (cfg%points(2, 1))
        cfg%points(:, 1) = [gp_r0 + 0.33_dp, 0.0_dp]
        psib_exact = exact_flux(cfg%points(1, 1), cfg%points(2, 1))
        state%psi = seed
        rcfg%max_refits = merge(4, 0, refit)
        rcfg%tolerance = 1.0e-3_dp
        call gs_fb_refit_solve(state, sys, law, cfg, currents, contour, 0.4_dp*h*h, &
                               chain, rcfg, report, status, set=coils, &
                               curve=outer_curve, extra_cb=manufactured_extra)
        call require(status == GS_FB_OK .or. .not. refit, 'refit solve')
        if (status /= GS_FB_OK .and. refit) print '(a,i4)', '   refit status ', status
        allocate (exact(size(state%points, 2)))
        do i = 1, size(exact)
            exact(i) = exact_flux(state%points(1, i), state%points(2, i))
        end do
        scale = maxval(abs(exact))
        error = maxval(abs(state%psi - exact))/scale
      call gs_fb_evaluate(state, sys, state%psi, currents, ev, status, jacobian=.false.)
        bound_error = abs(ev%psi_bound - psib_exact)/scale
        refits = report%refits
        moved = report%movement(report%refits)
    end subroutine solve_case

    ! ------------------------------------------------------------------
    ! Part 3: a free-boundary problem with a genuine source jump at the
    ! boundary (no manufactured cancellation), errors against a fine
    ! reference solution on a refitted P3 mesh.
    ! ------------------------------------------------------------------

    subroutine jump_law(law)
        !! Edge values (s = 1) are 70 % of the axis values: a strong jump.
        type(gs_profile_law_t), intent(out) :: law

        allocate (law%knots(2), law%pp(4, 1), law%ff(4, 1))
        law%knots = [0.0_dp, 1.0_dp]
        law%normalized = .true.
        law%axis = -1.0_dp
        law%p_edge = 0.0_dp
        law%f_edge = 5.0_dp
        law%pp(:, 1) = law_scale*1.0e5_dp*[1.0_dp, -0.2_dp, -0.1_dp, 0.0_dp]
        law%ff(:, 1) = law_scale*1.2_dp*[1.0_dp, -0.2_dp, -0.1_dp, 0.0_dp]
    end subroutine jump_law

    subroutine free_solve(order, n_outer, refit, summary, shape_error_point, status)
        !! summary = [psi_b, psi_a, Ip, psi at 6 probes, chain(1) radius].
        integer, intent(in) :: order, n_outer
        logical, intent(in) :: refit
        real(dp), intent(out) :: summary(10), shape_error_point(2)
        integer, intent(out) :: status
        type(gs_state_t) :: state
        type(gs_profile_law_t) :: law
        type(gs_fb_config_t) :: cfg
        type(gs_fb_system_t) :: sys
        type(gs_fb_eval_t) :: ev
        type(gs_fb_refit_config_t) :: rcfg
        type(gs_fb_refit_report_t) :: report
        type(gs_fb_history_t) :: hist
        real(dp), allocatable :: contour(:, :), chain(:, :), seed(:), trace(:, :)
        real(dp) :: perimeter, h, g(2), hess(2, 2), value
        integer :: j, nc, i, element, stat
        logical :: found
        real(dp), parameter :: probes(2, 6) = reshape([gp_r0, 0.0_dp, gp_r0 + 0.2_dp, &
            0.4_dp, gp_r0 - 0.25_dp, -0.5_dp, gp_r0 + 0.5_dp, -0.3_dp, gp_r0 + 0.8_dp, &
                                                0.0_dp, gp_r0 - 0.4_dp, 1.0_dp], [2, 6])

        nc = n_outer/2
        allocate (contour(2, n_outer), chain(2, nc))
        do j = 1, n_outer
            call outer_curve(real(j - 1, dp)/n_outer, contour(:, j))
        end do
        do j = 1, nc
            chain(:, j) = [gp_r0 + 0.97_dp*gp_a*cos(2.0_dp*pi*(j - 1)/nc), &
                           0.97_dp*gp_b*sin(2.0_dp*pi*(j - 1)/nc)]
        end do
        perimeter = pi*(3.0_dp*(gc_a + gc_b) &
                        - sqrt((3.0_dp*gc_a + gc_b)*(gc_a + 3.0_dp*gc_b)))
        h = perimeter/n_outer
        call initialize_gs(state, contour, 0.4_dp*h*h, status, element_order=order, &
                           curved=.true., curve=outer_curve, interior_chain=chain)
        if (status /= 0) return
        allocate (seed(size(state%points, 2)))
        do i = 1, size(seed)
            seed(i) = vacuum(state%points(1, i), state%points(2, i)) &
                      + bump(state%points(1, i), state%points(2, i), centre, semi, &
                             seed_amplitude)
        end do
        call jump_law(law)
        cfg%mode = GS_FB_LIMITER
        cfg%axis_sign = 1
        cfg%axis_guess = centre
        cfg%search_radius = 0.4_dp
        cfg%quadrature_tolerance = 1.0e-9_dp
        cfg%jacobian_tolerance = 1.0e-9_dp
        cfg%tolerance = 1.0e-9_dp
        allocate (cfg%points(2, 1))
        cfg%points(:, 1) = [gp_r0 + 0.33_dp, 0.0_dp]
        state%psi = seed
        rcfg%max_refits = merge(5, 0, refit)
        call gs_fb_refit_solve(state, sys, law, cfg, currents, contour, 0.4_dp*h*h, &
                            chain, rcfg, report, status, set=coils, curve=outer_curve, &
                               extra_cb=smooth_extra, history=hist)
        if (status == GS_FB_REFIT_NOT_CONVERGED) status = GS_FB_OK
        if (status /= GS_FB_OK) &
            print '(a,i4,a,i4,a,i4)', '   free solve status ', status, ' order ', &
            order, ' edges ', n_outer
        if (status /= GS_FB_OK) return
      call gs_fb_evaluate(state, sys, state%psi, currents, ev, status, jacobian=.false.)
        if (status /= GS_FB_OK) return
        summary(1) = ev%psi_bound
        summary(2) = ev%psi_axis
        summary(3) = ev%plasma_current
        do i = 1, 6
            call gs_fb_locate(state, probes(1, i), probes(2, i), element, g, found)
            call gs_element_jet(state, element, probes(1, i), probes(2, i), value, g, &
                                hess, stat)
            summary(3 + i) = value
        end do
     call gs_fb_trace_lcfs(state, state%psi, ev%axis, ev%psi_bound, 1.0_dp, 16, trace, &
                              stat)
        if (stat /= GS_FB_OK) then
            status = stat
            return
        end if
        summary(10) = trace(1, 1)
        shape_error_point = [real(report%refits, dp), report%movement(report%refits)]
    end subroutine free_solve

    subroutine reference_study()
        real(dp) :: ref(10), run(10), info(2), errors(4, 3), dummy(2)
        real(dp) :: fixed_finest(4), rates(4), fixed_rates(3)
        integer :: status, k, order, level, n, j
        character(len=2) :: tag
        integer, parameter :: sizes(3) = [18, 36, 72]

        call free_solve(3, 144, .true., ref, dummy, status)
        call require(status == GS_FB_OK, 'reference solve')
        if (status /= GS_FB_OK) return
        print '(a,10es12.4)', 'reference (P3, 144 edges, refit): ', ref
        do order = 2, 3
            do k = 1, 2
                do level = 1, 3
                    n = sizes(level)
                    call free_solve(order, n, k == 2, run, info, status)
                    call require(status == GS_FB_OK, 'jump-study solve')
                    if (status /= GS_FB_OK) cycle
                errors(1, level) = maxval(abs(run(4:9) - ref(4:9)))/abs(ref(2) - ref(1))
                    errors(2, level) = abs(run(1) - ref(1))/abs(ref(2) - ref(1))
                    errors(3, level) = abs(run(3) - ref(3))/abs(ref(3))
                    errors(4, level) = abs(run(10) - ref(10))
                    write (tag, '(a,i1)') 'P', order
             print '(a,a,a,i4,a,4es11.3,a,i2,a,es9.2)', tag, merge(' refit', ' fixed', &
                     k == 2), ' edges ', n, ' errors psi(probes), psi_b, Ip, R_lcfs ', &
                       errors(:, level), '  refits ', int(info(1)), '  dist/h ', info(2)
                end do
                rates = [(log(errors(j, 2)/errors(j, 3))/log(2.0_dp), j=1, 4)]
                print '(a,a,a,4f7.2)', tag, merge(' refit', ' fixed', k == 2), &
                    ' orders (finest pair) ', rates
                if (k == 1) fixed_finest = errors(:, 3)
                if (k == 2) then
                    ! The refit recovers the order of the flux level quantities
                    ! and the boundary position; the fixed mesh does not.
                    call require(all(errors(1:4, 3) < fixed_finest), &
                                 tag//': refit more accurate than the fixed mesh')
                    call require(rates(3) > 2.5_dp .and. rates(4) > 2.5_dp, &
                                 tag//': refit order of Ip and boundary position')
                    call require(rates(2) > 1.9_dp, tag//': refit order of psi_b')
                    call require(maxval(fixed_rates) < 2.5_dp, &
                                 tag//': fixed mesh loses the order')
                end if
                if (k == 1) fixed_rates = rates(2:4)
            end do
        end do
    end subroutine reference_study

end program test_gs_fb_refit
