!> Q7 part 2 demo: 1X1V Wigner-Poisson vs Vlasov-Poisson on one grid.
!> Quantum Langmuir (Bohm-Gross) dispersion, quantum two-stream growth,
!> hbar -> 0 convergence, unclipped Wigner negativity, error/cost ledger,
!> figures and the f_W | f_V movie. Units: omega_p, v_ref, m = 1.
program wigner_demo
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_wigner_poisson
    use kin6d_wigner_dispersion, only: beam_t, dispersion_root
    use kin6d_wigner_leaves, only: wigner_bohm_gross, wigner_cold_two_stream
    use kin6d_wigner_plots, only: plot_summary, plot_movie
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp)
    integer, parameter :: nlk = 6, nts = 15, nh = 3, nhb = 5, frames = 192
    integer, parameter :: mx = 64, mv = 256, hold = 24
    real(dp) :: hl(nh) = [0.0_dp, 1.0_dp, 2.0_dp], ht(nh) = [0.0_dp, 1.0_dp, 2.0_dp]
    real(dp) :: kl(nlk), wl_fit(nlk, nh), wl_root(nlk, nh), wl_bg(nlk, nh)
    real(dp) :: kt(nts), g_fit(nts, nh), g_root(nts, nh), g_cold(nts, nh)
    real(dp) :: hb(nhb), err_hb(nhb, 2), rate(2), errors(3, 3), cost(6), t0, t1
    real(dp) :: neg(frames, 3), amp(frames, 3), times(frames)
    real(dp), allocatable :: fw(:, :, :), fv(:, :, :)
    real(dp) :: ledger(3, 2), dgamma(3), u, bg
    complex(dp) :: root, droot, fit, guess
    type(beam_t) :: cool(1), streams(2)
    integer :: i, j, unit, status
    logical :: ok

    call execute_command_line('mkdir -p output', exitstat=status)
    cool(1) = beam_t(1.0_dp, 0.0_dp, 0.1_dp)
    streams = [beam_t(0.5_dp, 1.0_dp, 0.2_dp), beam_t(0.5_dp, -1.0_dp, 0.2_dp)]

    ! 1. Quantum Langmuir dispersion of a cool plasma (sigma = 0.1).
    call cpu_time(t0)
    open (newunit=unit, file='output/wigner-langmuir.csv', status='replace')
    write (unit, '(a)') '# k,hbar,omega_fit,gamma_fit,omega_root,gamma_root,bohm_gross'
    do j = 1, nh
        do i = 1, nlk
            kl(i) = 0.4_dp + 0.2_dp*(i - 1)
            call wigner_bohm_gross(kl(i), 0.1_dp, hl(j), bg)
            call dispersion_root(kl(i), hl(j), cool, cmplx(bg, 0, dp), root, droot, ok)
            call run_fit([1.0_dp], [0.0_dp], [0.1_dp], kl(i), hl(j), &
                         1.5_dp + hl(j)*kl(i), 256, 50.0_dp, 110.0_dp, 1e-9_dp, fit)
            wl_fit(i, j) = real(fit); wl_root(i, j) = real(root); wl_bg(i, j) = bg
            write (unit, '(f5.2,",",f4.1,5(",",es24.16e3))') kl(i), hl(j), fit, root, bg
        end do
    end do
    close (unit)
    call cpu_time(t1); cost(1) = t1 - t0

    ! 2. Two-stream growth rate scan; quantum suppression at large k.
    call cpu_time(t0)
    open (newunit=unit, file='output/wigner-two-stream.csv', status='replace')
    write (unit, '(a)') '# k,hbar,gamma_fit,gamma_root,gamma_cold_quartic'
    do j = 1, nh
        guess = (0.0_dp, 0.3_dp)
        do i = 1, nts
            kt(i) = 0.2_dp + 0.1_dp*(i - 1)
            call dispersion_root(kt(i), ht(j), streams, guess, root, droot, ok)
            if (.not. ok .or. aimag(root) < 0.01_dp) &
                call dispersion_root(kt(i), ht(j), streams, (0.0_dp, 0.3_dp), root, &
                                     droot, ok)
            g_root(i, j) = merge(aimag(root), 0.0_dp, ok .and. aimag(root) > 0)
            if (ok .and. aimag(root) > 0.01_dp) guess = root
            call wigner_cold_two_stream(kt(i), 1.0_dp, ht(j), u)
            g_cold(i, j) = sqrt(max(-u, 0.0_dp))
            g_fit(i, j) = -1
            if (g_root(i, j) > 0.03_dp) then
                call run_fit([0.5_dp, 0.5_dp], [1.0_dp, -1.0_dp], [0.2_dp, 0.2_dp], &
                    kt(i), ht(j), 4.0_dp, 256, 24.0_dp, &
                    min(60.0_dp, 24.0_dp + 16.0_dp*0.33_dp/g_root(i, j)), 1e-12_dp, fit)
                g_fit(i, j) = aimag(fit)
            end if
            write (unit, '(f5.2,",",f4.1,3(",",es24.16e3))') kt(i), ht(j), &
                g_fit(i, j), g_root(i, j), g_cold(i, j)
        end do
    end do
    close (unit)
    call cpu_time(t1); cost(2) = t1 - t0

    ! 3. Growth-rate derivative wrt hbar at k = 0.6, hbar = 1: discrete
    ! tangent (JVP) of log|n_1| over the linear window, equation-level
    ! implicit derivative of the dispersion root, and finite differences.
    call cpu_time(t0)
    call growth_tangent(0.6_dp, 1.0_dp, dgamma(1))
    call dispersion_root(0.6_dp, 1.0_dp, streams, (0.0_dp, 0.3_dp), root, droot, ok)
    dgamma(2) = aimag(droot)
    call dispersion_root(0.6_dp, 1.001_dp, streams, root, fit, droot, ok)
    dgamma(3) = aimag(fit)
    call dispersion_root(0.6_dp, 0.999_dp, streams, root, fit, droot, ok)
    dgamma(3) = (dgamma(3) - aimag(fit))/0.002_dp
    call cpu_time(t1); cost(3) = t1 - t0

    ! 4. hbar -> 0: Wigner -> Vlasov on the same grid (two-stream, t = 20,
    ! and Landau alpha = 0.05, t = 10).
    call cpu_time(t0)
    do i = 1, nhb
        hb(i) = 0.8_dp/2**(i - 1)
    end do
    call hbar_limit(err_hb, rate)
    call cpu_time(t1); cost(4) = t1 - t0

    ! 5. Error indicators: Strang dt, nv, nx (two-stream hbar = 1, t = 20).
    call cpu_time(t0)
    call error_ledger(errors)
    call cpu_time(t1); cost(5) = t1 - t0

    ! 6. Movie run and negativity histories (hbar = 1, 0.5, 0).
    call cpu_time(t0)
    allocate (fw(mv + 1, mx + 1, frames), fv(mv + 1, mx + 1, frames))
    call movie_runs()
    call cpu_time(t1); cost(6) = t1 - t0
    call write_summary()
    call plot_summary(kl, hl, wl_fit, wl_root, wl_bg, kt, g_fit, g_root, g_cold, &
                      hb, err_hb, times, neg)
    call plot_movie(fw, fv, times, neg, amp, 2*pi/0.6_dp, 4.0_dp, frames, hold)
contains
    include 'wigner_demo_runs.inc'
end program wigner_demo
