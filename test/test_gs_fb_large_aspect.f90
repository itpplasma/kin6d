module test_gs_fb_large_aspect_callbacks
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortfem_axisymmetric_gs_kernel, only: gs_mu0, gs_rect_coil_t
    use kin6d_grad_shafranov, only: gs_state_t, initialize_gs, gs_cell_data
    use kin6d_gs_coils, only: gs_circuit_set_t, gs_circuit_flux_vectors
    use kin6d_gs_profile_law, only: gs_profile_law_t, gs_law_primitives
    use kin6d_gs_free_boundary, only: gs_fb_config_t, gs_fb_system_t, gs_fb_eval_t, &
        gs_fb_setup, gs_fb_evaluate, gs_fb_newton, gs_fb_trace_lcfs, GS_FB_OK, &
        GS_FB_TARGET_MEAN
    use fortfem_triangle_duffy_quadrature, only: triangle_duffy_quadrature

    use kin6d_curved_triangle, only: boundary_curve_t
    use kin6d_grad_shafranov, only: initialize_gs_object
    implicit none
    private
    public :: run_test_gs_fb_large_aspect
    real(dp), parameter :: pi = acos(-1.0_dp), a = 0.5_dp
    real(dp), parameter :: eps_list(4) = [0.25_dp, 0.125_dp, 0.0625_dp, 0.03125_dp]
    type :: run_test_gs_fb_large_aspect_context_t
        real(dp), pointer :: r0 => null()
    end type

    type, extends(boundary_curve_t) :: run_test_gs_fb_large_aspect_outer_curve_t
        type(run_test_gs_fb_large_aspect_context_t), pointer :: data => null()
    contains
        procedure :: evaluate => run_test_gs_fb_large_aspect_outer_curve_binding
    end type

contains
    subroutine run_test_gs_fb_large_aspect()
    !! T5 (DESIGN 5.1, step P13): large-aspect-ratio free-boundary asymptotics.
    !!
    !! A circular plasma of minor radius a = 0.5 m at major radius R0 = a/eps
    !! inside the circle Gamma_c of radius 2a; the boundary flux is the mean
    !! over 16 points of the circle of radius a about (R0, 0) (FreeGS-style
    !! targets: the plasma may shift, R_g and a_g are measured). The
    !! vertical field comes from a ring pair at (3 R0, +-1.5 R0) (Helmholtz
    !! spacing, uniform near the plasma). For each eps the free-boundary
    !! solution gives Ip, the boundary centre R_g and half width a_g, the
    !! poloidal beta beta_p = 2 mu0 <p> / B_pa^2 and the internal inductance
    !! l_i = <B_p^2> / B_pa^2 (cross-section averages, B_pa = mu0 Ip / L_p).
    !! Oracle (Shafranov):
    !!
    !!     B_v = mu0 Ip / (4 pi R_g) (ln(8 R_g / a_g) + beta_p + l_i/2 - 3/2),
    !!
    !! against the coil field at (R_g, 0); the vertical-field current is
    !! rescaled twice to the Shafranov value at R0 so that the plasma is
    !! centred; the relative mismatch of the third pass is an asymptotic
    !! error that must fall as eps -> 0 (observed order printed). A
    !! quadrupole (coils at (R0, +-3a) carrying +I_q, at (R0 +- 3a, 0)
    !! carrying -I_q, ramped in two continuation steps) with I_q of the sign
    !! of Ip pulls the plasma vertically and raises the elongation kappa;
    !! the opposite sign lowers it.
    implicit none

    real(dp) , target :: r0, mismatch(size(eps_list)), kappa(3), &
        kappa_eps(size(eps_list))
    real(dp) :: rate
    integer :: failures, i
        type(run_test_gs_fb_large_aspect_context_t), target :: ctx
        type(run_test_gs_fb_large_aspect_outer_curve_t), target :: outer_curve_operator

    failures = 0
    do i = 1, size(eps_list)
        r0 = a/eps_list(i)
        call shafranov_case(0.0_dp, mismatch(i), kappa_eps(i))
    end do
    kappa(1) = kappa_eps(2)
    do i = 2, size(eps_list)
        rate = log(abs(mismatch(i - 1)/mismatch(i)))/log(2.0_dp)
        print '(a,f8.5,a,f6.2)', 'eps ', eps_list(i), '  observed order ', rate
        call require(abs(mismatch(i)) < abs(mismatch(i - 1)), &
                     'Shafranov mismatch falls with eps')
    end do
    call require(rate > 1.5_dp, 'Shafranov mismatch order above 1.5 (finest pair)')
    call require(abs(mismatch(size(eps_list))) < 2.0e-3_dp, &
                 'Shafranov mismatch below 0.2 % at eps = 1/32')
    r0 = a/eps_list(2)
    call shafranov_case(0.1_dp, mismatch(1), kappa(2))
    call shafranov_case(-0.1_dp, mismatch(1), kappa(3))
    print '(a,3f9.5)', 'kappa without, with +, with - quadrupole ', kappa
    call require(kappa(2) > kappa(1), 'like-signed quadrupole elongates')
    call require(kappa(3) < kappa(1), 'opposite quadrupole squeezes')
    if (failures /= 0) error stop 'gs_fb_large_aspect failures'
    print '(a)', 'gs_fb_large_aspect (T5) checks passed'

contains
    subroutine bind_run_test_gs_fb_large_aspect_context()
        ctx%r0 => r0
        outer_curve_operator%data => ctx
    end subroutine bind_run_test_gs_fb_large_aspect_context


    subroutine require(ok, what)
        logical, intent(in) :: ok
        character(len=*), intent(in) :: what
        if (.not. ok) then
            failures = failures + 1
            print '(2a)', 'FAIL: ', what
        end if
    end subroutine require


    subroutine shafranov_case(quadrupole, mismatch, kappa)
        !! quadrupole: current of each quadrupole coil per ampere of Ip.
        real(dp), intent(in) :: quadrupole
        real(dp), intent(out) :: mismatch, kappa
        integer, parameter :: n_outer = 96
        type(gs_state_t) :: state
        type(gs_circuit_set_t) :: set
        type(gs_profile_law_t) :: law
        type(gs_fb_config_t) :: cfg
        type(gs_fb_system_t) :: sys
        type(gs_fb_eval_t) :: ev
        real(dp) :: contour(2, n_outer), h, currents(2), g(1, 2)
        real(dp) :: gr(1, 2), gz(1, 2), ip, b_shafranov, b_coil, beta_p, l_i
        real(dp) :: r_geo, a_geo, rho2
        integer :: status, j, pass
        do j = 1, n_outer
        call bind_run_test_gs_fb_large_aspect_context()
            call run_test_gs_fb_large_aspect_outer_curve(ctx, real(j - 1, &
                dp)/n_outer, contour(:, j))
        end do
        h = 2*pi*2*a/n_outer
        call bind_run_test_gs_fb_large_aspect_context()
        call initialize_gs_object(state, contour, 0.4_dp*h*h, status, element_order=2, &
                           curved=.true., curve=outer_curve_operator)
        call require(status == 0, 'mesh')
        allocate (set%coils(6), set%incidence(6, 2))
        set%incidence = 0.0_dp
        set%incidence(1:2, 1) = 1.0_dp
        set%incidence(3:4, 2) = 1.0_dp
        set%incidence(5:6, 2) = -1.0_dp
        set%coils(1)%rect = gs_rect_coil_t(3*r0, 1.5_dp*r0, 0.1_dp, 0.1_dp, 1.0_dp)
        set%coils(2)%rect = gs_rect_coil_t(3*r0, -1.5_dp*r0, 0.1_dp, 0.1_dp, 1.0_dp)
        set%coils(3)%rect = gs_rect_coil_t(r0, 3*a, 0.1_dp, 0.1_dp, 1.0_dp)
        set%coils(4)%rect = gs_rect_coil_t(r0, -3*a, 0.1_dp, 0.1_dp, 1.0_dp)
        set%coils(5)%rect = gs_rect_coil_t(r0 + 3*a, 0.0_dp, 0.1_dp, 0.1_dp, 1.0_dp)
        set%coils(6)%rect = gs_rect_coil_t(r0 - 3*a, 0.0_dp, 0.1_dp, 0.1_dp, 1.0_dp)
        ! Normalised law, source ~ 1.2 (mu0 j) at every R0: p' and FF' vanish at
        ! the edge.
        allocate (law%knots(2), law%pp(4, 1), law%ff(4, 1))
        law%knots = [0.0_dp, 1.0_dp]
        law%normalized = .true.
        law%axis = -1.0_dp
        law%p_edge = 0.0_dp
        law%f_edge = 5.0_dp*r0
        law%pp(:, 1) = 0.6_dp/(gs_mu0*r0)*[1.0_dp, -1.0_dp, 0.0_dp, 0.0_dp]
        law%ff(:, 1) = 0.6_dp*r0*[1.0_dp, -1.0_dp, 0.0_dp, 0.0_dp]
        cfg%mode = GS_FB_TARGET_MEAN
        cfg%axis_guess = [r0, 0.0_dp]
        cfg%search_radius = 0.4_dp
        allocate (cfg%points(2, 16))
        do j = 1, 16
            cfg%points(:, j) = [r0 + a*cos(2*pi*(j - 0.5_dp)/16), &
                                a*sin(2*pi*(j - 0.5_dp)/16)]
        end do
        call gs_fb_setup(state, sys, law, cfg, set, status)
        call require(status == GS_FB_OK, 'setup')
        ! First guess: Ip ~ 1.2 pi a^2 / (2 mu0), beta_p + l_i/2 ~ 1.
        ip = 1.2_dp*pi*a*a/(2*gs_mu0)
        call gs_circuit_flux_vectors(set, [r0], [0.0_dp], g, gr, gz)
        b_shafranov = gs_mu0*ip/(4*pi*r0)*(log(8*r0/a) - 0.5_dp)
        currents = [-b_shafranov/(gr(1, 1)/r0), 0.0_dp]
        do j = 1, size(state%psi)
            rho2 = ((state%points(1, j) - r0)**2 + state%points(2, j)**2)/a**2
            call gs_circuit_flux_vectors(set, state%points(1, j:j), &
                                         state%points(2, j:j), g)
            state%psi(j) = dot_product(g(1, :), currents) &
                           + 0.075_dp*r0*max(0.0_dp, 1 - rho2)**2
        end do
        ! Three centring passes without the quadrupole, then the quadrupole
        ! current ramped in two continuation steps (last pass measured).
        do pass = 1, merge(3, 5, quadrupole == 0.0_dp)
            if (pass > 3) currents(2) = quadrupole*ip*(pass - 3)/2.0_dp
            call gs_fb_newton(state, sys, currents, status)
            if (status /= GS_FB_OK) print '(a,i3)', '   Newton status ', status
            call require(status == GS_FB_OK, 'Newton converged')
            call measure(state, sys, currents, ev, ip, r_geo, a_geo, kappa, beta_p, l_i)
            call gs_circuit_flux_vectors(set, [r_geo], [0.0_dp], g, gr, gz)
            b_coil = -dot_product(gr(1, :), currents)/r_geo
            b_shafranov = gs_mu0*ip/(4*pi*r_geo)*(log(8*r_geo/a_geo) + beta_p &
                                                  + l_i/2 - 1.5_dp)
            mismatch = (b_coil - b_shafranov)/b_shafranov
            print '(a,f8.5,a,f5.2,a,i1,a,es11.4,a,f8.4,a,f7.4,a,2f7.4,a,f7.4,' &
                //'a,es10.2)', &
                'eps ', a/r0, ' q ', quadrupole, ' pass ', pass, ' Ip ', ip, &
                ' R_g-R0 ', r_geo - r0, ' a_g ', a_geo, ' beta_p,l_i ', beta_p, &
                l_i, ' kappa ', kappa, ' B_v mismatch ', mismatch
            ! Rescale the vertical field to the Shafranov value at R0 so that
            ! the next pass is centred.
            if (pass >= 3) cycle
            call gs_circuit_flux_vectors(set, [r0], [0.0_dp], g, gr, gz)
            currents(1) = currents(1)*gs_mu0*ip/(4*pi*r0)*(log(8*r0/a_geo) + beta_p &
                + l_i/2 - 1.5_dp)/(-dot_product(gr(1, :), currents)/r0)
        end do
    end subroutine shafranov_case

    subroutine measure(state, sys, currents, ev, ip, r_geo, a_geo, kappa, beta_p, l_i)
        !! Boundary geometry from the traced LCFS; cross-section averages of p
        !! and B_p^2 over {s < 1} in the plasma cells (Duffy rule of degree 12).
        type(gs_state_t), intent(inout) :: state
        type(gs_fb_system_t), intent(in) :: sys
        real(dp), intent(in) :: currents(:)
        type(gs_fb_eval_t), intent(inout) :: ev
        real(dp), intent(out) :: ip, r_geo, a_geo, kappa, beta_p, l_i
        type(gs_profile_law_t) :: law
        real(dp), allocatable :: chain(:, :), x(:), y(:), w(:), psi(:)
        real(dp) :: shape(6), grads(2, 6), hess(2, 2, 6), pt(2), det, area, p_int
        real(dp) :: b2_int, perimeter, b_pa, ph, grad(2), pressure, f
        integer :: status, e, q, k
        psi = state%psi
        call gs_fb_evaluate(state, sys, psi, currents, ev, status, jacobian=.false.)
        call require(status == GS_FB_OK, 'evaluate')
        ip = ev%plasma_current
        call gs_fb_trace_lcfs(state, psi, ev%axis, ev%psi_bound, 1.0_dp, 256, chain, &
                              status)
        call require(status == GS_FB_OK, 'LCFS trace')
        if (status /= GS_FB_OK) then
            r_geo = huge(1.0_dp); a_geo = 0; kappa = 0; beta_p = 0; l_i = 0
            return
        end if
        r_geo = 0.5_dp*(maxval(chain(1, :)) + minval(chain(1, :)))
        a_geo = 0.5_dp*(maxval(chain(1, :)) - minval(chain(1, :)))
        kappa = (maxval(chain(2, :)) - minval(chain(2, :)))/(2*a_geo)
        perimeter = 0.0_dp
        do k = 1, size(chain, 2)
            perimeter = perimeter + norm2(chain(:, mod(k, size(chain, 2)) + 1) &
                                          - chain(:, k))
        end do
        law = sys%law
        law%axis = ev%psi_axis - ev%psi_bound
        call triangle_duffy_quadrature(12, x, y, w, status)
        area = 0.0_dp
        p_int = 0.0_dp
        b2_int = 0.0_dp
        do e = 1, size(state%triangles, 2)
            if (.not. ev%plasma_cell(e)) cycle
            do q = 1, size(w)
                call gs_cell_data(state, e, [x(q), y(q)], pt, det, shape, grads, &
                                  hess, status)
                ph = dot_product(shape, psi(state%triangles(:, e))) - ev%psi_bound
                if (ph <= 0.0_dp) cycle
                grad = matmul(grads, psi(state%triangles(:, e)))
                call gs_law_primitives(law, ph, pressure, f, status)
                area = area + w(q)*det
                p_int = p_int + w(q)*det*pressure
                b2_int = b2_int + w(q)*det*sum(grad**2)/pt(1)**2
            end do
        end do
        b_pa = gs_mu0*ip/perimeter
        beta_p = 2*gs_mu0*(p_int/area)/b_pa**2
        l_i = (b2_int/area)/b_pa**2
    end subroutine measure


    end subroutine run_test_gs_fb_large_aspect

    subroutine run_test_gs_fb_large_aspect_outer_curve(ctx, t, point)
        type(run_test_gs_fb_large_aspect_context_t), intent(in) :: ctx
        real(dp), intent(in) :: t
        real(dp), intent(out) :: point(2)
        point = [ctx%r0 + 2*a*cos(2*pi*t), 2*a*sin(2*pi*t)]
    end subroutine run_test_gs_fb_large_aspect_outer_curve

    subroutine run_test_gs_fb_large_aspect_outer_curve_binding(self, t, point)
        class(run_test_gs_fb_large_aspect_outer_curve_t), intent(in) :: self
        real(dp), intent(in) :: t
        real(dp), intent(out) :: point(2)
        call run_test_gs_fb_large_aspect_outer_curve(self%data, t, point)
    end subroutine
end module test_gs_fb_large_aspect_callbacks

program test_gs_fb_large_aspect
    use test_gs_fb_large_aspect_callbacks, only: run_test_gs_fb_large_aspect
    implicit none
    call run_test_gs_fb_large_aspect()
end program test_gs_fb_large_aspect
