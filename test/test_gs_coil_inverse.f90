program test_gs_coil_inverse
    !! P9 / T6 (DESIGN 3.4 and 5.1): the inverse coil-current problem on a twin
    !! free-boundary problem with known currents I*.
    !!
    !!  1. forward sensitivities against central differences of the converged
    !!     flux, and thread-count independence;
    !!  2. T6: constraints extracted from the solution at I*, Gauss-Newton from
    !!     perturbed currents recovers I* (psi_b free: up to the GSVD null
    !!     mode, which the report identifies as the flux gauge; with noisy
    !!     targets the error on each mode obeys |c_j| gsv_j <= |noise|; psi_b
    !!     fixed: the null mode is removed);
    !!  3. the gradient D^T r + lambda^2 L^2 (I - I_ref) against central
    !!     differences of the converged objective;
    !!  4. the three fit variants, bounds and the series tie honoured exactly;
    !!  5. FreeGS-emulation mode reaches a higher objective than Gauss-Newton.
    use, intrinsic :: iso_fortran_env, only: dp => real64
!$  use omp_lib
    use kin6d_grad_shafranov, only: gs_state_t
    use kin6d_gs_coils, only: gs_circuit_set_t, gs_circuit_flux_vectors
    use kin6d_gs_free_boundary, only: gs_fb_system_t, gs_fb_eval_t, gs_fb_factor_t, &
                     gs_fb_evaluate, gs_fb_factor, gs_fb_release, gs_fb_newton, GS_FB_OK
    use kin6d_gs_coil_inverse
    use gs_ci_fixture
    implicit none

    type(gs_state_t) :: state
    type(gs_fb_system_t) :: sys
    type(gs_circuit_set_t) :: set
    type(gs_ci_model_t) :: model
    real(dp), allocatable :: psi_star(:)
    real(dp) :: psib_star, pts(2, 8)
    integer :: failures, status, n_outer
    character(len=32) :: arg

    failures = 0
    n_outer = 24
    if (command_argument_count() .gt. 0) then
        call get_command_argument(1, arg)
        read (arg, *) n_outer       ! mesh size for the thread-scaling measurement
    end if
    call build_twin(n_outer, 2, state, sys, set, psi_star, status)
    call require(status .eq. 0, 'twin problem')
    if (status .ne. 0) error stop 'gs_coil_inverse: no twin'
    call target_points(pts)
    call make_model(.false., model)
    if (command_argument_count() .gt. 1) then
        call test_scaling()
        stop
    end if
    call test_sensitivities()
    call test_recovery()
    call test_gradient()
    call test_variants()
    call test_emulation()
    if (failures .ne. 0) error stop 'gs_coil_inverse failures'
    print '(a)', 'gs_coil_inverse (P9, T6) oracles passed'

contains

    subroutine require(ok, what)
        logical, intent(in) :: ok
        character(len=*), intent(in) :: what
        if (.not. ok) then
            failures = failures + 1
            print '(2a)', 'FAIL: ', what
        end if
    end subroutine require

    subroutine make_model(with_psib, mdl)
        !! Rows: psi(p) - psi_b at eight points, B_R and B_Z at two points, and
        !! with `with_psib` the flux state psi_b. Targets from the twin.
        logical, intent(in) :: with_psib
        type(gs_ci_model_t), intent(out) :: mdl
        type(gs_ci_row_t), allocatable :: rows(:)
        type(gs_fb_eval_t) :: ev
        integer :: k, nr, stat
        real(dp), parameter :: tol = 1.0e-3_dp

        nr = 12
        if (with_psib) nr = 13
        allocate (rows(nr))
        do k = 1, 8
            rows(k) = gs_ci_row_t(GS_CI_ISOFLUX_LCFS, pts(:, k), 0.0_dp, 0.0_dp, tol)
        end do
        rows(9) = gs_ci_row_t(GS_CI_XPOINT_BR, [3.05_dp, -0.7_dp], 0.0_dp, 0.0_dp, tol)
        rows(10) = gs_ci_row_t(GS_CI_XPOINT_BZ, [3.05_dp, -0.7_dp], 0.0_dp, 0.0_dp, tol)
        rows(11) = gs_ci_row_t(GS_CI_XPOINT_BR, [2.6_dp, 0.6_dp], 0.0_dp, 0.0_dp, tol)
        rows(12) = gs_ci_row_t(GS_CI_XPOINT_BZ, [3.5_dp, 0.6_dp], 0.0_dp, 0.0_dp, tol)
     if (with_psib) rows(13) = gs_ci_row_t(GS_CI_PSI_BOUND, 0.0_dp, 0.0_dp, 0.0_dp, tol)
        call gs_ci_build_model(state, sys, set, rows, mdl, stat)
        call require(stat .eq. GS_CI_OK, 'model')
        call gs_fb_evaluate(state, sys, psi_star, istar, ev, stat, jacobian=.false.)
        psib_star = ev%psi_bound
        call gs_ci_set_targets(mdl, psi_star, psib_star)
    end subroutine make_model

    function start_currents() result(i0)
        real(dp) :: i0(n_circ)
        i0 = istar*[1.25_dp, 0.8_dp, 1.2_dp, 0.75_dp, 1.3_dp, 0.7_dp, 1.1_dp, 0.9_dp]
    end function start_currents

    subroutine add_noise(mdl, amplitude)
        !! Deterministic noise of size amplitude * scale on every target.
        type(gs_ci_model_t), intent(inout) :: mdl
        real(dp), intent(in) :: amplitude
        integer :: i

        do i = 1, mdl%n_rows
            mdl%rows(i)%target = mdl%rows(i)%target &
                                 + amplitude*mdl%rows(i)%scale*sin(1.7_dp*i + 0.3_dp)
        end do
    end subroutine add_noise

    subroutine test_scaling()
        !! Strong scaling of the sensitivity solves (manual: test <n_outer> x).
        type(gs_fb_eval_t) :: ev
        type(gs_fb_factor_t) :: fac
        real(dp), allocatable :: ds(:, :), ref(:, :)
        real(dp) :: t, best
        integer :: stat, nt, rep

        call gs_fb_evaluate(state, sys, psi_star, istar, ev, stat)
        call gs_fb_factor(state, sys, ev, fac, stat)
        do nt = 1, 8
!$          call omp_set_num_threads(nt)
            call gs_ci_flux_sensitivities(sys, ev, fac, ds, stat)
            best = huge(1.0_dp)
            do rep = 1, 5
                t = wtime()
                call gs_ci_flux_sensitivities(sys, ev, fac, ds, stat)
                best = min(best, wtime() - t)
            end do
            if (nt .eq. 1) ref = ds
            print '(a,i2,a,es10.3,a,es10.2,a,i7)', 'threads ', nt, '  time ', best, &
                '  max diff to serial ', maxval(abs(ds - ref)), '  nodes ', sys%n
        end do
        call gs_fb_release(fac)
    end subroutine test_scaling

    subroutine test_sensitivities()
        type(gs_fb_eval_t) :: ev
        type(gs_fb_factor_t) :: fac
        real(dp), allocatable :: ds(:, :), ds1(:, :), plus(:), minus(:)
        real(dp) :: h, err, t0, t1, tn, t_one
        integer :: stat, k, kk(2)

        call gs_fb_evaluate(state, sys, psi_star, istar, ev, stat)
        call gs_fb_factor(state, sys, ev, fac, stat)
        call require(stat .eq. GS_FB_OK, 'factor')
        call gs_ci_flux_sensitivities(sys, ev, fac, ds, stat)    ! warm-up
        t0 = wtime()
        call gs_ci_flux_sensitivities(sys, ev, fac, ds, stat)
        t1 = wtime()
        tn = t1 - t0
        call require(stat .eq. GS_CI_OK, 'sensitivities')
!$      call omp_set_num_threads(1)
        t0 = wtime()
        call gs_ci_flux_sensitivities(sys, ev, fac, ds1, stat)
        t_one = wtime() - t0
!$      call omp_set_num_threads(omp_get_num_procs())
        err = maxval(abs(ds - ds1))/maxval(abs(ds))
        print '(a,es10.2,a,f7.4,a,f7.4,a)', 'sensitivities: thread-count difference ', &
            err, ';  1 thread ', t_one, ' s, all threads ', tn, ' s'
      call require(err .lt. 1.0e-12_dp, 'sensitivities independent of the thread count')
        kk = [1, 8]
        do k = 1, 2
            h = 20.0_dp
            plus = psi_star
            call solve_at(istar + h*unit(kk(k)), plus)
            minus = psi_star
            call solve_at(istar - h*unit(kk(k)), minus)
            err = maxval(abs((plus - minus)/(2.0_dp*h) - ds(:, kk(k)))) &
                  /maxval(abs(ds(:, kk(k))))
            print '(a,i2,a,es10.2)', 'd psi / d I_', kk(k), ' against central FD ', err
      call require(err .lt. 2.0e-6_dp, 'forward sensitivity against finite differences')
        end do
        call gs_fb_release(fac)
    end subroutine test_sensitivities

    function unit(k) result(e)
        integer, intent(in) :: k
        real(dp) :: e(n_circ)
        e = 0.0_dp
        e(k) = 1.0_dp
    end function unit

    function wtime() result(t)
        real(dp) :: t
        integer :: c, rate
        call system_clock(c, rate)
        t = real(c, dp)/real(rate, dp)
!$      t = omp_get_wtime()
    end function wtime

    subroutine solve_at(cur, psi)
        real(dp), intent(in) :: cur(:)
        real(dp), intent(inout) :: psi(:)
        integer :: stat

        state%psi = psi
        call gs_fb_newton(state, sys, cur, stat)
        call require(stat .eq. GS_FB_OK, 'forward solve')
        psi = state%psi
    end subroutine solve_at

    subroutine test_recovery()
        type(gs_ci_model_t) :: noisy, fixed
        type(gs_ci_config_t) :: cfg
        type(gs_ci_result_t) :: res, free
        real(dp), allocatable :: coef(:)
        real(dp) :: flux(8, n_circ), mean, spread(2), noise_norm, worst, c_free
        integer :: status, j

        cfg%use_bounds = .false.
        ! Noise-free data, psi_b free.
       call gs_ci_solve(state, sys, model, cfg, start_currents(), psi_star, res, status)
        print '(a,i2,a,i3,a,es10.2,a,es10.2)', 'T6 noise-free: status ', status, &
            '  iterations ', res%iterations, '  Phi ', res%phi, '  |I - I*|/|I*| ', &
            norm2(res%currents - istar)/norm2(istar)
        print '(a,8es10.2)', '   GSVD spectrum ', res%gsvd%gsv
        call require(status .eq. GS_CI_OK, 'T6 converged')
        call require(res%phi .lt. 1.0e-9_dp, 'T6 constraints reproduced')
        call require(norm2(res%currents - istar) .lt. 1.0e-4_dp*norm2(istar &
                                                             ), 'T6 currents recovered')
        ! The flux-gauge null mode: smallest generalised singular value, a
        ! coil combination whose vacuum flux is nearly uniform over the targets.
        call gs_circuit_flux_vectors(set, pts(1, :), pts(2, :), flux)
        do j = 1, 2
            mean = sum(matmul(flux, res%gsvd%modes(:, j)))/8.0_dp
            spread(j) = maxval(abs(matmul(flux, res%gsvd%modes(:, j)) - mean)) &
                        /max(abs(mean), tiny(1.0_dp))
        end do
        print '(a,2es10.2,a,es10.2)', '   flux non-uniformity of modes 1, 2 ', spread, &
            '  gsv(2)/gsv(1) ', res%gsvd%gsv(2)/res%gsvd%gsv(1)
       call require(res%gsvd%gsv(2) .gt. 50.0_dp*res%gsvd%gsv(1), 'null mode separated')
        call require(spread(1) .lt. 0.1_dp*spread(2), 'null mode is the flux gauge')
        call require(all(abs(res%gsvd%gsv - [(norm2(matmul(dmat_of(res), &
            res%gsvd%modes(:, j))), j=1, n_circ)]) .lt. 1.0e-8_dp*res%gsvd%gsv(n_circ) &
                         ), 'modes reproduce the generalised singular values')
       call require(abs(norm2(res%gsvd%modes(:, 1)/set%i_max) - 1.0_dp) .lt. 1.0e-12_d &
                     p, 'modes normalised in the L norm')

        ! Noisy targets: |c_j| gsv_j <= |noise| mode by mode (exact for the
        ! linearised problem), with the error concentrated in the null mode.
        noisy = model
        call add_noise(noisy, 1.0_dp)
        noise_norm = sqrt(real(noisy%n_rows, dp))
        call gs_ci_solve(state, sys, noisy, cfg, start_currents(), psi_star, free, &
                         status)
        call require(status .eq. GS_CI_OK, 'noisy fit converged')
        coef = gs_ci_mode_coefficients(free%gsvd, free%currents - istar)
        worst = maxval(abs(coef)*free%gsvd%gsv)
        c_free = abs(coef(1))
        print '(a,es10.2,a,es10.2,a,8es9.1)', '   noisy: max |c_j| gsv_j ', worst, &
            '  bound ', noise_norm, '  |c_j| ', abs(coef)
        call require(worst .le. 1.2_dp*noise_norm, 'mode errors bounded by the noise')
        call require(c_free .gt. 5.0_dp*maxval(abs(coef(3:)) &
                                         ), 'noise error concentrates in the null mode')
        ! Fixed psi_b: the extra row removes the null mode.
        call make_model(.true., fixed)
        call add_noise(fixed, 1.0_dp)
        call gs_ci_solve(state, sys, fixed, cfg, start_currents(), psi_star, res, &
                         status)
        coef = gs_ci_mode_coefficients(res%gsvd, res%currents - istar)
        print '(a,i2,a,es10.2,a,es10.2)', '   psi_b fixed: status ', status, &
            '  gsv(1) ', res%gsvd%gsv(1), '  max |c_j| ', maxval(abs(coef))
        call require(status .eq. GS_CI_OK, 'fixed-psi_b fit converged')
        call require(res%gsvd%gsv(1) .gt. 20.0_dp*free%gsvd%gsv(1 &
                                                     ), 'psi_b row lifts the null mode')
        call require(maxval(abs(coef)) .lt. 0.2_dp*c_fre &
                     e, 'psi_b row removes the null-mode error')
    end subroutine test_recovery

    function dmat_of(res) result(d)
        !! Jacobian D at the solution of `res` (full response).
        type(gs_ci_result_t), intent(in) :: res
        real(dp), allocatable :: d(:, :), grad(:), psi(:)
        type(gs_ci_config_t) :: cfg
        real(dp) :: phi
        integer :: status

        psi = res%psi
        call gs_ci_objective(state, sys, model, cfg, res%currents, psi, phi, grad, &
                             status, dmat=d)
    end function dmat_of

    subroutine test_gradient()
        type(gs_ci_config_t) :: cfg
        real(dp), allocatable :: grad(:), psi(:), gfd(:), dummy(:)
        real(dp) :: phi, pp, pm, h, err
        real(dp) :: cur(n_circ)
        integer :: status, k

        cfg%lambda = 0.05_dp
        allocate (cfg%i_ref(n_circ))
        cfg%i_ref = 0.8_dp*istar
        cur = start_currents()
        psi = psi_star
        call gs_ci_objective(state, sys, model, cfg, cur, psi, phi, grad, status)
        call require(status .eq. GS_CI_OK, 'objective')
        allocate (gfd(n_circ))
        do k = 1, n_circ
            h = 1.0e-3_dp*set%i_max(k)
            psi = psi_star
            call gs_ci_objective(state, sys, model, cfg, cur + h*unit(k), psi, pp, &
                                 dummy, status)
            psi = psi_star
            call gs_ci_objective(state, sys, model, cfg, cur - h*unit(k), psi, pm, &
                                 dummy, status)
            gfd(k) = (pp - pm)/(2.0_dp*h)
        end do
        err = maxval(abs(grad - gfd))/maxval(abs(gfd))
        print '(a,es10.2,a,es10.2)', 'gradient vs central FD of Phi: ', err, &
            '  Phi ', phi
        print '(a,8es10.2)', '   analytic ', grad
        print '(a,8es10.2)', '   FD       ', gfd
    call require(err .lt. 1.0e-5_dp, 'Gauss-Newton gradient against finite differences')
    end subroutine test_gradient

    subroutine test_variants()
        type(gs_ci_model_t) :: lim
        type(gs_ci_config_t) :: cfg
        type(gs_ci_result_t) :: unc, bnd, fix
        integer :: status

        lim = model
        ! Limits tighter than I* on circuits 1 and 4; the series tie I7 = I8.
        lim%lower = -current_limit
        lim%upper = current_limit
        lim%upper(1) = 0.8e4_dp
        lim%lower(4) = -0.9e4_dp
        allocate (cfg%tie_c(1, n_circ), cfg%tie_d(1))
        cfg%tie_c = 0.0_dp
        cfg%tie_c(1, 7) = 1.0_dp
        cfg%tie_c(1, 8) = -1.0_dp
        cfg%tie_d = 0.0_dp
        cfg%lambda = 1.0e-3_dp
        call gs_ci_fit(state, sys, lim, cfg, GS_CI_UNCONSTRAINED, start_currents(), &
                       psi_star, unc, status)
        call require(status .eq. GS_CI_OK, 'unconstrained fit')
        call gs_ci_fit(state, sys, lim, cfg, GS_CI_BOUNDED, start_currents(), &
                       psi_star, bnd, status)
        call require(status .eq. GS_CI_OK, 'bounded fit')
        call gs_ci_fit(state, sys, lim, cfg, GS_CI_FIXED_AT_LIMIT, start_currents(), &
                       psi_star, fix, status)
        call require(status .eq. GS_CI_OK, 'fixed-at-limit fit')
        print '(a,3es12.4)', 'variants Phi (unconstrained, bounded, fixed) ', unc%phi, &
            bnd%phi, fix%phi
        print '(a,8f9.1)', '   unconstrained ', unc%currents
        print '(a,8f9.1)', '   bounded       ', bnd%currents
        print '(a,8f9.1)', '   fixed-at-limit', fix%currents
        print '(a,8i3,a,i2,a,8l2)', '   bounded active set ', bnd%active, &
            '   passes ', fix%passes, '  pinned ', fix%pinned
     call require(unc%currents(1) .gt. lim%upper(1) .or. unc%currents(4) .lt. lim%lowe &
                     r(4), 'the unconstrained fit violates a limit')
    call require(all(bnd%currents .le. lim%upper) .and. all(bnd%currents .ge. lim%lowe &
                                           r), 'bounded fit honours the bounds exactly')
  call require(bnd%currents(1) .eq. lim%upper(1) .and. bnd%currents(4) .eq. lim%lower( &
                     4), 'bounded fit sits on the active limits')
       call require(abs(bnd%currents(7) - bnd%currents(8)) .lt. 1.0e-9_dp*current_limi &
                     t, 'bounded fit honours the series tie')
        call require(fix%pinned(1) .and. fix%pinned(4) .and. count(fix%pinned) .eq. 2 &
                     , 'fixed-at-limit pins the over-limit circuits')
  call require(fix%currents(1) .eq. lim%upper(1) .and. fix%currents(4) .eq. lim%lower( &
                     4), 'pinned circuits sit exactly at the limit')
    call require(all(fix%currents .le. lim%upper) .and. all(fix%currents .ge. lim%lowe &
                                                 r), 'fixed-at-limit fit within limits')
       call require(abs(fix%currents(7) - fix%currents(8)) .lt. 1.0e-9_dp*current_limi &
                     t, 'fixed-at-limit fit honours the series tie')
        call require(unc%phi .le. bnd%phi*(1.0_dp + 1.0e-9_dp) + 1.0e-12_dp &
                     , 'bounds cannot lower Phi')
        call require(bnd%phi .le. fix%phi*(1.0_dp + 1.0e-6_dp) &
                     , 'bounded optimum no worse than fixed-at-limit')
    end subroutine test_variants

    subroutine test_emulation()
        !! The measurable "at least as good as FreeGS" claim (DESIGN 3.4): both
        !! solutions are scored by the same full-response objective.
        type(gs_ci_model_t) :: noisy
        type(gs_ci_config_t) :: cfg
        type(gs_ci_result_t) :: gn, fg, zero_ref, scen_ref
        real(dp), allocatable :: d_full(:, :), d_vac(:, :), grad(:), psi(:)
        real(dp) :: phi, response, gn_scaled, fg_scaled
        integer :: status

        noisy = model
        call add_noise(noisy, 3.0_dp)
        cfg%use_bounds = .false.
        cfg%lambda = 0.05_dp
        allocate (cfg%i_ref(n_circ))
        cfg%i_ref = 0.0_dp
        psi = psi_star
        call gs_ci_objective(state, sys, noisy, cfg, istar, psi, phi, grad, status, &
                             dmat=d_full)
        cfg%emulate_freegs = .true.
        psi = psi_star
        call gs_ci_objective(state, sys, noisy, cfg, istar, psi, phi, grad, status, &
                             dmat=d_vac)
        response = norm2(d_full - d_vac)/norm2(d_full)
        print '(a,es10.2)', 'plasma response |D - D_vac| / |D| ', response
        call require(response .gt. 0.05_dp, 'the plasma response is significant')
        cfg%emulate_freegs = .false.
        call gs_ci_solve(state, sys, noisy, cfg, start_currents(), psi_star, gn, status)
        call require(status .eq. GS_CI_OK, 'Gauss-Newton converged')
        cfg%emulate_freegs = .true.
        cfg%max_iterations = 80
        call gs_ci_solve(state, sys, noisy, cfg, start_currents(), psi_star, fg, status)
        print '(a,i2,a,i3)', 'emulation status ', status, '  iterations ', fg%iterations
        call require(status .eq. GS_CI_OK, 'FreeGS emulation converged')
        cfg%emulate_freegs = .false.
        psi = psi_star
       call gs_ci_objective(state, sys, noisy, cfg, fg%currents, psi, phi, grad, status)
        print '(a,2es14.6,a,es10.2)', 'Phi Gauss-Newton, FreeGS emulation ', gn%phi, &
            phi, '  |grad Phi| at the emulation point ', norm2(grad)
        call require(phi .gt. gn%phi*(1.0_dp + 1.0e-4_dp &
                         ), 'FreeGS emulation has a higher objective than Gauss-Newton')
        gn_scaled = norm2(gn%gradient*set%i_max)
        fg_scaled = norm2(grad*set%i_max)
        print '(a,2es10.2)', '   |L^-1 grad Phi| Gauss-Newton, emulation ', gn_scaled, &
            fg_scaled
        call require(fg_scaled .gt. 100.0_dp*gn_scale &
                     d, 'Gauss-Newton is stationary, the emulation is not')

        ! Reference currents: a scenario and zero. A dominant regulariser pulls
        ! the fit to the reference.
        cfg%emulate_freegs = .false.
        cfg%lambda = 1.0e4_dp
        cfg%i_ref = 0.5_dp*istar
        call gs_ci_solve(state, sys, noisy, cfg, start_currents(), psi_star, scen_ref, &
                         status)
        call require(status .eq. GS_CI_OK, 'scenario-reference fit')
        deallocate (cfg%i_ref)
        call gs_ci_solve(state, sys, noisy, cfg, start_currents(), psi_star, zero_ref, &
                         status)
        call require(status .eq. GS_CI_OK, 'zero-reference fit')
     print '(a,2es10.2)', 'dominant regulariser: |I - I_ref|/|I_ref| scenario, zero ', &
            norm2(scen_ref%currents - 0.5_dp*istar)/norm2(0.5_dp*istar), &
            norm2(zero_ref%currents)/norm2(istar)
        call require(norm2(scen_ref%currents - 0.5_dp*istar) .lt. 0.05_dp*norm2(istar &
                                                       ), 'scenario reference honoured')
        call require(norm2(zero_ref%currents) .lt. 0.05_dp*norm2(istar &
                                                           ), 'zero reference honoured')
    end subroutine test_emulation

end program test_gs_coil_inverse
