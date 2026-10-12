program test_gs_fb_greens
    !! T4 (DESIGN 5.1): cross-check of the FEM-BEM free-boundary solution
    !! against a von Hagenow-style Green's-function evaluation that shares
    !! only the kernel and the coil model with the solver.
    !!
    !! Part A, the cut source. For the nodal flux psi* of the manufactured
    !! problem and a normalised law with a finite edge value, the plasma load
    !! b_i = int f(psi*) 1{plasma} v_i dA of the solver (exact level-cut
    !! quadrature, flood-fill plasma, jump at s = 1) is compared through
    !! sum_i b_i g_i = int f I[g] dA for smooth g with an independent polar
    !! cubature about the axis whose LCFS radius is found by bisection on each
    !! ray, so the jump is resolved by construction.
    !!
    !! Part B, the converged solution. The equilibrium of a law that vanishes
    !! at the edge, with a continuous forcing so that the problem has a known
    !! solution, is computed by Newton. The total source of the solution,
    !! f(psi_h) 1{plasma} + forcing, is integrated against the axisymmetric
    !! Green's function on the polar grid (graded about the axis; pieces
    !! split at the LCFS ray radius and at the edge of the forcing support):
    !!   * psi_h = sum_k I_k g_k + int G f_total dA on every Gamma_c trace
    !!     node and at points between the LCFS and Gamma_c, and on the
    !!     magnetic axis (logarithmically singular integrand);
    !!   * Ampere's law: the exterior flux of S_h through Gamma_c is mu0 Ip.
    !! The agreement is the discretisation error of the FEM-BEM solution and
    !! must fall under mesh refinement.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortfem_axisymmetric_gs_kernel, only: gs_green, gs_mu0
    use fortnum_quadrature, only: gauss_legendre
    use kin6d_grad_shafranov, only: gs_state_t
    use kin6d_gs_coils, only: gs_circuit_set_t, gs_circuit_flux_vectors
    use kin6d_gs_coupling, only: gs_coupling_source_load, GS_COUPLING_OK
    use kin6d_gs_profile_law, only: gs_profile_law_t, gs_law_sources
    use kin6d_gs_free_boundary, only: gs_fb_config_t, gs_fb_system_t, gs_fb_eval_t, &
     gs_fb_history_t, gs_fb_setup, gs_fb_evaluate, gs_fb_newton, GS_FB_OK, GS_FB_LIMITER
    use gs_coupling_fixture, only: build_state, build_coils, gp_r0, fe_value
    use gs_fb_fixture, only: make_law, bump, bump_laplacian
    implicit none

    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp), parameter :: currents(2) = [1.5e4_dp, -1.0e4_dp]
    real(dp), parameter :: centre(2) = [gp_r0 + 0.02_dp, 0.03_dp]
    real(dp), parameter :: semi(2) = [0.65_dp, 1.0_dp]
    real(dp), parameter :: amplitude = 1.5_dp
    real(dp), parameter :: law_scale = 0.25_dp
    real(dp), parameter :: limiter(2) = [gp_r0 + 0.33_dp, 0.0_dp]
    type(gs_circuit_set_t) :: coils
    real(dp), allocatable :: yr(:), yz(:), yw(:), yf(:)
    ! Reference levels of the forcing (any fixed numbers; the problem is
    ! defined by them), set once from the analytic flux.
    type(gs_profile_law_t) :: law_edge, law_zero
    real(dp) :: ref_axis, ref_bound
    integer :: failures

    failures = 0
    call build_coils(coils)
    call reference_levels()
    call part_a(3, 40)
    call part_a(2, 40)
    call part_b(3, [24, 48])
    call part_b(2, [24, 48])
    if (failures /= 0) error stop 'gs_fb_greens failures'
    print '(a)', 'gs_fb_greens (T4) oracles passed'

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
    end function exact_flux

    subroutine reference_levels()
        !! psi_a of psi* by a few Newton steps on finite-difference gradients,
        !! psi_b at the limiter point.
        real(dp) :: p(2), g(2), h(2, 2), d, det, step(2)
        integer :: it

        d = 1.0e-4_dp
        p = centre
        do it = 1, 20
            g = [exact_flux(p(1) + d, p(2)) - exact_flux(p(1) - d, p(2)), &
                 exact_flux(p(1), p(2) + d) - exact_flux(p(1), p(2) - d)]/(2.0_dp*d)
            h(1, 1) = (exact_flux(p(1) + d, p(2)) - 2.0_dp*exact_flux(p(1), p(2)) &
                       + exact_flux(p(1) - d, p(2)))/d**2
            h(2, 2) = (exact_flux(p(1), p(2) + d) - 2.0_dp*exact_flux(p(1), p(2)) &
                       + exact_flux(p(1), p(2) - d))/d**2
            h(1, 2) = (exact_flux(p(1) + d, p(2) + d) - exact_flux(p(1) + d, p(2) - d) &
                    - exact_flux(p(1) - d, p(2) + d) + exact_flux(p(1) - d, p(2) - d)) &
                      /(4.0_dp*d*d)
            h(2, 1) = h(1, 2)
            det = h(1, 1)*h(2, 2) - h(1, 2)**2
            step = -[h(2, 2)*g(1) - h(1, 2)*g(2), -h(2, 1)*g(1) + h(1, 1)*g(2)]/det
            p = p + step
            if (maxval(abs(step)) < 1.0e-9_dp) exit
        end do
        ref_axis = exact_flux(p(1), p(2))
        ref_bound = exact_flux(limiter(1), limiter(2))
        call make_law(law_edge, .true., law_scale)
        call make_law(law_zero, .false., law_scale)
    end subroutine reference_levels

    function law_source(law, psi, r) result(f)
        !! f = mu0 R p' + FF'/R for the normalised law at flux psi, with the
        !! reference levels, zero for s >= 1.
        type(gs_profile_law_t), intent(in) :: law
        real(dp), intent(in) :: psi, r
        real(dp) :: f, pp, ff
        type(gs_profile_law_t) :: shifted

        f = 0.0_dp
        if ((psi - ref_bound)*(ref_axis - ref_bound) <= 0.0_dp) return
        shifted = law
        shifted%axis = ref_axis - ref_bound
        call gs_law_sources(shifted, psi - ref_bound, pp, ff)
        f = ff/r + gs_mu0*pp*r
    end function law_source

    function forcing(r, z) result(f)
        !! Continuous forcing of part B: f* - f_law(psi*).
        real(dp), intent(in) :: r, z
        real(dp) :: f

        f = bump_laplacian(r, z, centre, semi, amplitude) &
            - law_source(law_zero, exact_flux(r, z), r)
    end function forcing

    function forcing_density(r, z) result(j)
        real(dp), intent(in) :: r, z
        real(dp) :: j

        j = forcing(r, z)/gs_mu0
    end function forcing_density

    subroutine ray_radius(state, psi, psi_bound, origin, theta, radius)
        !! First radius along the ray where psi falls to psi_bound.
        type(gs_state_t), intent(inout) :: state
        real(dp), intent(in) :: psi(:), psi_bound, origin(2), theta
        real(dp), intent(out) :: radius
        real(dp) :: lo, hi, mid, value
        integer :: bis
        logical :: found

        lo = 0.0_dp
        hi = 1.6_dp
        do bis = 1, 60
            mid = 0.5_dp*(lo + hi)
            call fe_value(state, psi, origin(1) + mid*cos(theta), &
                          origin(2) + mid*sin(theta), value, found)
            if (found .and. value > psi_bound) then
                lo = mid
            else
                hi = mid
            end if
        end do
        radius = 0.5_dp*(lo + hi)
    end subroutine ray_radius

    function support_radius(origin, theta) result(radius)
        !! Radius along the ray at which the forcing bump ends (its ellipse).
        real(dp), intent(in) :: origin(2), theta
        real(dp) :: radius, a, b, c, ex, ey

        ex = cos(theta)
        ey = sin(theta)
        a = (ex/semi(1))**2 + (ey/semi(2))**2
        b = 2.0_dp * ((origin(1) - centre(1)) * ex / semi(1)**2 &
                      + (origin(2) - centre(2)) * ey / semi(2)**2)
c = ((origin(1) - centre(1))/semi(1))**2 + ((origin(2) - centre(2))/semi(2))**2 - 1.0_dp
        radius = (-b + sqrt(b*b - 4.0_dp*a*c))/(2.0_dp*a)
    end function support_radius

    subroutine part_a(order, n_outer)
        !! Cut source against an independent polar cubature, edge-jump law.
        integer, intent(in) :: order, n_outer
        type(gs_state_t) :: state
        type(gs_fb_config_t) :: cfg
        type(gs_fb_system_t) :: sys
        type(gs_fb_eval_t) :: ev
        real(dp), allocatable :: exact(:), nodal(:), gx(:), gw(:), tests(:, :)
        real(dp) :: radius, theta, value, fval, gval, total(4), ref(4), solver(4), r, z
        real(dp) :: weight, worst
        integer, parameter :: ntheta = 256, nrad = 96
        integer :: status, i, k, q, n, j
        logical :: found

        call build_state(state, n_outer, order, status)
        n = size(state%points, 2)
        allocate (exact(n), nodal(n), tests(n, 4))
        do i = 1, n
            exact(i) = exact_flux(state%points(1, i), state%points(2, i))
            do j = 1, 4
                tests(i, j) = test_function(j, state%points(1, i), state%points(2, i))
            end do
        end do
        cfg%mode = GS_FB_LIMITER
        cfg%axis_guess = centre
        cfg%search_radius = 0.4_dp
        cfg%quadrature_tolerance = 1.0e-11_dp
        allocate (cfg%points(2, 1))
        cfg%points(:, 1) = limiter
        call gs_fb_setup(state, sys, law_edge, cfg, coils, status)
        call require(status == GS_FB_OK, 'part A setup')
        call gs_fb_evaluate(state, sys, exact, currents, ev, status, jacobian=.false.)
        call require(status == GS_FB_OK, 'part A evaluate')
        ! Test functions 1, R, Z, R Z, interpolated at the nodes.
        solver = 0.0_dp
        do j = 1, 4
            solver(j) = dot_product(ev%source, tests(:, j))
        end do
        allocate (gx(nrad), gw(nrad))
        call gauss_legendre(nrad, gx, gw)
        gx = 0.5_dp*(gx + 1.0_dp)
        gw = 0.5_dp*gw
        ref = 0.0_dp
        do k = 1, ntheta
            theta = 2.0_dp*pi*(k - 0.5_dp)/ntheta
            call ray_radius(state, exact, ev%psi_bound, ev%axis, theta, radius)
            do q = 1, nrad
                r = ev%axis(1) + radius*gx(q)*cos(theta)
                z = ev%axis(2) + radius*gx(q)*sin(theta)
                call fe_value(state, exact, r, z, value, found)
                fval = law_value(value, r, ev%psi_axis, ev%psi_bound)
                weight = 2.0_dp*pi/ntheta*gw(q)*radius**2*gx(q)
                do j = 1, 4
                    call fe_value(state, tests(:, j), r, z, gval, found)
                    ref(j) = ref(j) + weight*fval*gval
                end do
            end do
        end do
        worst = maxval(abs(solver - ref)/abs(ref))
        total = ref
        print '(a,i1,a,i3,a,4es12.4)', 'part A P', order, ' edges ', n_outer, &
            '  polar cubature ', total
        print '(a,4es12.4)', '          cut-cell load   ', solver
        print '(a,es10.2)', '          max relative difference ', worst
        call require(worst < merge(1.0e-6_dp, 5.0e-6_dp, order == 3), &
                     'cut source load against polar cubature')
    end subroutine part_a

    function test_function(j, r, z) result(g)
        integer, intent(in) :: j
        real(dp), intent(in) :: r, z
        real(dp) :: g

        select case (j)
        case (1)
            g = 1.0_dp
        case (2)
            g = r
        case (3)
            g = z + 1.0_dp
        case default
            g = r*(z + 1.0_dp)
        end select
    end function test_function

    function law_value(psi, r, psi_axis, psi_bound) result(f)
        real(dp), intent(in) :: psi, r, psi_axis, psi_bound
        real(dp) :: f, pp, ff
        type(gs_profile_law_t) :: shifted

        f = 0.0_dp
        if ((psi - psi_bound)*(psi_axis - psi_bound) <= 0.0_dp) return
        shifted = law_edge
        shifted%axis = psi_axis - psi_bound
        call gs_law_sources(shifted, psi - psi_bound, pp, ff)
        f = ff/r + gs_mu0*pp*r
    end function law_value

    subroutine part_b(order, levels)
        integer, intent(in) :: order, levels(:)
        real(dp) :: errors(3, size(levels))
        integer :: i

        do i = 1, size(levels)
            call compare(order, levels(i), errors(:, i))
        end do
        print '(a,i1,a,3f7.2)', 'part B P', order, &
            ' observed rates (trace, exterior+axis, Ampere) ', &
            log(errors(:, 1)/errors(:, 2))/log(2.0_dp)
     call require(all(errors(:, 2) < errors(:, 1)), 'Green cross-check error decreases')
        call require(all(log(errors(1:2, 1)/errors(1:2, 2))/log(2.0_dp) > 2.0_dp), &
                     'Green cross-check convergence rate above 2')
        call require(errors(1, 2) < merge(1.0e-5_dp, 5.0e-5_dp, order == 3), &
                     'trace nodes against Green integral')
        call require(errors(2, 2) < merge(2.0e-5_dp, 3.0e-3_dp, order == 3), &
                     'exterior samples and axis against Green integral')
        call require(errors(3, 2) < 1.0e-3_dp, 'Ampere law with the Green-grid current')
    end subroutine part_b

    subroutine compare(order, n_outer, errors)
        !! errors = [max trace-node difference / flux range, max exterior and
        !! axis difference / flux range, relative Ampere error].
        integer, intent(in) :: order, n_outer
        real(dp), intent(out) :: errors(3)
        type(gs_state_t) :: state
        type(gs_fb_config_t) :: cfg
        type(gs_fb_system_t) :: sys
        type(gs_fb_eval_t) :: ev
        type(gs_fb_history_t) :: history
        integer, parameter :: ntheta = 96, nrad = 36
        real(dp), allocatable :: seed(:), trace_psi(:), ref(:), gx(:), gw(:), f_load(:)
        real(dp), allocatable :: ring(:)
        real(dp) :: radius(ntheta), flux_range, value, ampere, plasma_ip, sample(2, 5)
        real(dp) :: theta, lo, hi, piece, x0, w0, xr, xz
        integer :: status, i, k, n, m, np, piece_index
        logical :: found

        call build_state(state, n_outer, order, status)
        call require(status == 0, 'mesh')
        n = size(state%points, 2)
        cfg%mode = GS_FB_LIMITER
        cfg%axis_guess = centre
        cfg%search_radius = 0.4_dp
        cfg%quadrature_tolerance = 1.0e-11_dp
        cfg%jacobian_tolerance = 1.0e-11_dp
        allocate (cfg%points(2, 1))
        cfg%points(:, 1) = limiter
        call gs_fb_setup(state, sys, law_zero, cfg, coils, status)
        call require(status == GS_FB_OK, 'setup')
        call gs_coupling_source_load(state, forcing_density, f_load, status)
        call require(status == GS_COUPLING_OK, 'forcing load')
        sys%extra_load = f_load
        allocate (seed(n))
        do i = 1, n
            seed(i) = vacuum(state%points(1, i), state%points(2, i)) &
                      + bump(state%points(1, i), state%points(2, i), &
                             centre + [0.05_dp, -0.04_dp], semi, 0.6_dp*amplitude)
        end do
        state%psi = seed
        call gs_fb_newton(state, sys, currents, status, history)
        if (status /= GS_FB_OK) print '(a,i4)', '   Newton status ', status
        call require(status == GS_FB_OK, 'Newton converged')
      call gs_fb_evaluate(state, sys, state%psi, currents, ev, status, jacobian=.false.)
        call require(status == GS_FB_OK, 'final evaluation')

        ! Polar grid about the axis: pieces [0, rho_lcfs] and [rho_lcfs, rho_support]
        ! (or reversed), Gauss points per piece, first piece graded.
        allocate (gx(nrad), gw(nrad))
        call gauss_legendre(nrad, gx, gw)
        gx = 0.5_dp*(gx + 1.0_dp)
        gw = 0.5_dp*gw
        np = 2*ntheta*nrad
        if (allocated(yr)) deallocate (yr, yz, yw, yf)
        allocate (yr(np), yz(np), yw(np), yf(np))
        plasma_ip = 0.0_dp
        m = 0
        do k = 1, ntheta
            theta = 2.0_dp*pi*(k - 0.5_dp)/ntheta
            call ray_radius(state, state%psi, ev%psi_bound, ev%axis, theta, radius(k))
            lo = min(radius(k), support_radius(ev%axis, theta))
            hi = max(radius(k), support_radius(ev%axis, theta))
            do piece_index = 1, 2
                do i = 1, nrad
                    m = m + 1
                    if (piece_index == 1) then
                        x0 = lo*gx(i)**2
                        w0 = gw(i)*2.0_dp*gx(i)*lo
                    else
                        x0 = lo + (hi - lo)*gx(i)
                        w0 = gw(i)*(hi - lo)
                    end if
                    xr = ev%axis(1) + x0*cos(theta)
                    xz = ev%axis(2) + x0*sin(theta)
                    yr(m) = xr
                    yz(m) = xz
                    yw(m) = 2.0_dp*pi/ntheta*w0*x0
                    call fe_value(state, state%psi, xr, xz, value, found)
                    call require(found, 'grid point located')
                    yf(m) = law_value_zero(value, xr, ev%psi_axis, ev%psi_bound) &
                            + forcing(xr, xz)
                end do
            end do
        end do
        piece = sum(yw*yf)/gs_mu0
        ! Ip of the law part alone: the law current is what the exterior sees
        ! together with the forcing; Ampere uses the total.
        plasma_ip = piece

        np = sys%coupling%n_trace
        allocate (trace_psi(np), ref(np))
        flux_range = ev%psi_axis - minval(state%psi(sys%coupling%trace_node))
        !$omp parallel do default(shared) private(i) schedule(static)
        do i = 1, np
            ref(i) = vacuum(state%points(1, sys%coupling%trace_node(i)), &
                            state%points(2, sys%coupling%trace_node(i))) &
                     + green_sum(state%points(1, sys%coupling%trace_node(i)), &
                                 state%points(2, sys%coupling%trace_node(i)))
        end do
        !$omp end parallel do
        trace_psi = state%psi(sys%coupling%trace_node)
        errors(1) = maxval(abs(trace_psi - ref))/flux_range

        sample(:, 1) = [gp_r0 + 0.92_dp, 0.3_dp]
        sample(:, 2) = [gp_r0 - 0.85_dp, -0.4_dp]
        sample(:, 3) = [gp_r0 + 0.1_dp, 1.25_dp]
        sample(:, 4) = [gp_r0 + 0.2_dp, -1.3_dp]
        sample(:, 5) = ev%axis
        errors(2) = 0.0_dp
        do i = 1, 5
            call fe_value(state, state%psi, sample(1, i), sample(2, i), value, found)
            call require(found, 'sample point located')
            errors(2) = max(errors(2), abs(value - vacuum(sample(1, i), sample(2, i)) &
                                    - green_sum(sample(1, i), sample(2, i)))/flux_range)
        end do
        allocate (ring(np))
        ! The Steklov-Poincare map acts on the exterior-decaying plasma part
        ! of the flux; the coil flux is regular outside Gamma_c and has no
        ! share in it.
        do i = 1, np
          ring(i) = trace_psi(i) - vacuum(state%points(1, sys%coupling%trace_node(i)), &
                                            state%points(2, sys%coupling%trace_node(i)))
        end do
        ring = matmul(sys%coupling%s_trace, ring)
        ampere = sum(ring)/gs_mu0
        ! The manufactured flux bump carries no net current, so the error is
        ! measured against the absolute current content int |f| dA / mu0.
        errors(3) = abs(ampere - plasma_ip)/(sum(abs(yw*yf))/gs_mu0)
    print '(a,i1,a,i3,a,i2,a,es11.3,a,es11.3,a,es11.3)', 'part B P', order, ' edges ', &
            n_outer, '  Newton its ', history%iterations, '  psi_axis ', ev%psi_axis, &
            '  psi_b ', ev%psi_bound, '  total current ', plasma_ip
        print '(a,3es10.2)', '   trace, exterior+axis, Ampere errors ', errors
    end subroutine compare

    function law_value_zero(psi, r, psi_axis, psi_bound) result(f)
        real(dp), intent(in) :: psi, r, psi_axis, psi_bound
        real(dp) :: f, pp, ff
        type(gs_profile_law_t) :: shifted

        f = 0.0_dp
        if ((psi - psi_bound)*(psi_axis - psi_bound) <= 0.0_dp) return
        shifted = law_zero
        shifted%axis = psi_axis - psi_bound
        call gs_law_sources(shifted, psi - psi_bound, pp, ff)
        f = ff/r + gs_mu0*pp*r
    end function law_value_zero

    function green_sum(r, z) result(v)
        real(dp), intent(in) :: r, z
        real(dp) :: v
        integer :: q

        v = 0.0_dp
        do q = 1, size(yw)
            v = v + yw(q)*yf(q)*gs_green(r, z, yr(q), yz(q))
        end do
    end function green_sum

end program test_gs_fb_greens
