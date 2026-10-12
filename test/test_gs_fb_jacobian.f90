program test_gs_fb_jacobian
    !! T7 (DESIGN 5.1): the free-boundary Newton Jacobian against central
    !! differences of the residual, including the interface term of the moving
    !! source cut and the rank-two axis and boundary terms; transpose pairing
    !! <J^T l, d> = <l, J d>; factor/solve consistency of J and J^T; and the
    !! negative controls (dropping the boundary rank-one term or the interface
    !! term must be visible). The flux is analytic bumps sampled at the nodes,
    !! so no equilibrium is needed.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_grad_shafranov, only: gs_state_t
    use kin6d_gs_profile_law, only: gs_profile_law_t
    use kin6d_gs_free_boundary, only: gs_fb_config_t, gs_fb_system_t, gs_fb_eval_t, &
      gs_fb_factor_t, gs_fb_setup, gs_fb_evaluate, gs_fb_jacobian_apply, gs_fb_factor, &
                    gs_fb_solve, gs_fb_release, GS_FB_OK, GS_FB_XPOINT, GS_FB_LIMITER, &
                                      GS_FB_TARGET_MEAN
    use gs_coupling_fixture, only: build_state, gp_r0, gp_a, gp_b
    use gs_fb_fixture, only: make_law, bump, smooth_direction
    implicit none

    real(dp), parameter :: main_c(2) = [gp_r0, 0.0_dp], main_s(2) = [gp_a, gp_b]
    real(dp), parameter :: low_c(2) = [gp_r0, -1.15_dp], low_s(2) = [0.5_dp, 0.45_dp]
    integer :: failures

    failures = 0
    call run_case('limiter, edge jump, P2', 2, GS_FB_LIMITER, .true., .false., 1.0_dp)
    call run_case('target mean, edge jump, P3', 3, GS_FB_TARGET_MEAN, .true., &
                  .false., 1.0_dp)
 call run_case('FB-cut c=0.85, edge jump, P3', 3, GS_FB_XPOINT, .true., .true., 0.85_dp)
    call run_case('separatrix, vanishing edge, P2', 2, GS_FB_XPOINT, .false., &
                  .true., 1.0_dp)
    if (failures /= 0) error stop 'gs_fb_jacobian failures'
    print '(a)', 'gs_fb_jacobian (T7) oracles passed'

contains

    subroutine require(ok, what)
        logical, intent(in) :: ok
        character(len=*), intent(in) :: what
        if (.not. ok) then
            failures = failures + 1
            print '(2a)', 'FAIL: ', what
        end if
    end subroutine require

    function flux(r, z, second) result(psi)
        real(dp), intent(in) :: r, z
        logical, intent(in) :: second
        real(dp) :: psi

        psi = bump(r, z, main_c, main_s, 1.5_dp)
        if (second) psi = psi + bump(r, z, low_c, low_s, 0.9_dp)
    end function flux

    subroutine run_case(label, order, mode, jump, second, fraction)
        character(len=*), intent(in) :: label
        integer, intent(in) :: order, mode
        logical, intent(in) :: jump, second
        real(dp), intent(in) :: fraction
        type(gs_state_t) :: state
        type(gs_profile_law_t) :: law
        type(gs_fb_config_t) :: cfg
        type(gs_fb_system_t) :: sys
        type(gs_fb_eval_t) :: ev, ev_plus, ev_minus, ev_control
        type(gs_fb_factor_t) :: fac
        real(dp), allocatable :: psi(:), delta(:), jd(:), fd(:), lambda(:), jtl(:)
        real(dp), allocatable :: x(:), r(:), none(:), probe(:)
        real(dp) :: eps, err, norm_jd, pair_a, pair_b, err_rank, err_noint
        integer :: status, i, n, m
        real(dp) :: angle

        call build_state(state, 36, order, status)
        call require(status == 0, 'mesh')
        call make_law(law, jump, 1.0_dp)
        cfg%mode = mode
        cfg%axis_sign = 1
        cfg%axis_guess = [gp_r0 + 0.02_dp, 0.02_dp]
        cfg%xpoint_guess = [gp_r0, -0.62_dp]
        cfg%search_radius = 0.35_dp
        cfg%cut_fraction = fraction
        cfg%quadrature_tolerance = 1.0e-12_dp
        cfg%jacobian_tolerance = 1.0e-12_dp
        if (mode == GS_FB_LIMITER) then
            allocate (cfg%points(2, 2))
            cfg%points(:, 1) = [gp_r0 + 0.5_dp, 0.0_dp]
            cfg%points(:, 2) = [gp_r0 - 0.45_dp, 0.1_dp]
        else if (mode == GS_FB_TARGET_MEAN) then
            m = 8
            allocate (cfg%points(2, m))
            do i = 1, m
                angle = 6.283185307179586_dp*(i - 0.5_dp)/m
                cfg%points(:, i) = [gp_r0 + 0.38_dp*cos(angle), 0.62_dp*sin(angle)]
            end do
        end if
        call gs_fb_setup(state, sys, law, cfg, status=status)
        call require(status == GS_FB_OK, 'setup')
        n = sys%n
        allocate (psi(n), none(0))
        do i = 1, n
            psi(i) = flux(state%points(1, i), state%points(2, i), second)
        end do
        call smooth_direction(state, delta)
        call gs_fb_evaluate(state, sys, psi, none, ev, status)
        call require(status == GS_FB_OK, 'evaluate')
        if (status /= GS_FB_OK) return

        allocate (jd(n), fd(n), jtl(n), lambda(n), probe(n))
        call gs_fb_jacobian_apply(state, sys, ev, delta, jd)
        norm_jd = norm2(jd)
        eps = 1.0e-4_dp
        call gs_fb_evaluate(state, sys, psi + eps*delta, none, ev_plus, status)
        call require(status == GS_FB_OK, 'evaluate plus')
        call gs_fb_evaluate(state, sys, psi - eps*delta, none, ev_minus, status)
        call require(status == GS_FB_OK, 'evaluate minus')
        fd = (ev_plus%residual - ev_minus%residual)/(2.0_dp*eps)
        err = norm2(jd - fd)/norm_jd
        print '(a,a,a,es10.2,a,es10.2)', trim(label), ': Jd vs central FD ', &
            'relative ', err, '  psi_b ', ev%psi_bound
        print '(a,3es12.4)', '   psi_axis, psi_bound, Ip ', ev%psi_axis, ev%psi_bound, &
            ev%plasma_current
     call require(err < 2.0e-6_dp, trim(label)//': Jacobian against finite differences')

        ! Negative controls: the terms under test change J visibly.
        ev_control = ev
        ev_control%ub = 0.0_dp
        call gs_fb_jacobian_apply(state, sys, ev_control, delta, probe)
        err_rank = norm2(probe - fd)/norm_jd
        print '(a,es10.2)', '   without the boundary rank-one term ', err_rank
  call require(err_rank > 100.0_dp*err, trim(label)//': boundary rank-one term matters')
        if (jump) then
            ev_control = ev
            call strip_interface(state, ev_control)
            call gs_fb_jacobian_apply(state, sys, ev_control, delta, probe)
            err_noint = norm2(probe - fd)/norm_jd
            print '(a,es10.2)', '   without the interface term ', err_noint
         call require(err_noint > 100.0_dp*err, trim(label)//': interface term matters')
        end if

        ! Transpose pairing and solves.
        allocate (x(n), r(n))
        do i = 1, n
            lambda(i) = cos(0.9_dp*state%points(1, i) + 1.3_dp*state%points(2, i))
        end do
        call gs_fb_jacobian_apply(state, sys, ev, lambda, jtl, transpose=.true.)
        pair_a = dot_product(jtl, delta)
        pair_b = dot_product(lambda, jd)
        print '(a,es10.2)', '   transpose pairing relative gap ', &
            abs(pair_a - pair_b)/max(abs(pair_a), abs(pair_b))
        call require(abs(pair_a - pair_b) <= 1.0e-12_dp*max(abs(pair_a), abs(pair_b)), &
                     trim(label)//': transpose pairing')
        call gs_fb_factor(state, sys, ev, fac, status)
        call require(status == GS_FB_OK, 'factor')
        call gs_fb_solve(sys, ev, fac, lambda, x, status)
        call require(status == GS_FB_OK, 'solve')
        call gs_fb_jacobian_apply(state, sys, ev, x, r)
        print '(a,es10.2)', '   J x = rhs residual ', norm2(r - lambda)/norm2(lambda)
    call require(norm2(r - lambda) <= 1.0e-9_dp*norm2(lambda), trim(label)//': J solve')
        call gs_fb_solve(sys, ev, fac, lambda, x, status, transpose=.true.)
        call require(status == GS_FB_OK, 'transpose solve')
        call gs_fb_jacobian_apply(state, sys, ev, x, r, transpose=.true.)
        print '(a,es10.2)', '   J^T x = rhs residual ', norm2(r - lambda)/norm2(lambda)
        call require(norm2(r - lambda) <= 1.0e-9_dp*norm2(lambda), &
                     trim(label)//': J^T solve')
        call gs_fb_release(fac)
    end subroutine run_case

    subroutine strip_interface(state, control)
        !! Remove the interface part of M from a copy of the evaluation.
        type(gs_state_t), intent(in) :: state
        type(gs_fb_eval_t), intent(inout) :: control
        integer :: e, k

        do e = 1, size(state%triangles, 2)
            do k = 1, size(state%triangles, 1)
                control%ub(state%triangles(k, e)) = control%ub(state%triangles(k, e)) &
                                                  + sum(control%interface_mass(k, :, e))
            end do
        end do
        control%mass = control%mass - control%interface_mass
    end subroutine strip_interface

end program test_gs_fb_jacobian
