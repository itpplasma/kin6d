!> Figures and movie of the Q7 Wigner demo (presentation only; every plotted
!> number is also written to output/wigner_demo.nc and the CSV files).
module kin6d_wigner_plots
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortplot, only: figure, subplot, plot, title, suptitle, xlabel, ylabel, &
        xlim, ylim, legend, tight_layout, savefig, pcolormesh, colorbar, &
        set_xscale, set_yscale, get_global_figure, figure_t
    use fortplot_animation, only: animation_t, animation_func => FuncAnimation, &
        save_animation
    implicit none
    private
    public :: plot_summary, plot_movie
    real(dp), parameter :: blue(3) = [0.000_dp, 0.447_dp, 0.698_dp]
    real(dp), parameter :: verm(3) = [0.835_dp, 0.369_dp, 0.000_dp]
    real(dp), parameter :: green(3) = [0.000_dp, 0.620_dp, 0.451_dp]
    real(dp), parameter :: black(3) = [0.0_dp, 0.0_dp, 0.0_dp]
    real(dp), parameter :: grey(3) = [0.45_dp, 0.45_dp, 0.45_dp]
    ! Movie data (module-held for the frame callback of FuncAnimation).
    real(dp), allocatable :: mfw(:, :, :), mfv(:, :, :), mt(:), mneg(:, :), mamp(:, :)
    real(dp), allocatable :: mx(:), mv(:)
    real(dp) :: mrange, mnegmax, mresolved
    logical :: still_frame = .false.
    integer :: mframes
    class(figure_t), pointer :: canvas => null()
contains

    subroutine plot_summary(kl, hl, wfit, wroot, wbg, kt, gfit, groot, gcold, hb, &
                            errhb, times, neg)
        real(dp), intent(in) :: kl(:), hl(:), wfit(:, :), wroot(:, :), wbg(:, :)
        real(dp), intent(in) :: kt(:), gfit(:, :), groot(:, :), gcold(:, :)
        real(dp), intent(in) :: hb(:), errhb(:, :), times(:), neg(:, :)
        real(dp) :: colors(3, 3)
        character(16) :: lab
        integer :: j
        colors(:, 1) = blue; colors(:, 2) = verm; colors(:, 3) = green
        call figure(figsize=[13.2_dp, 8.8_dp])
        call subplot(2, 2, 1)
        call title('Langmuir waves, $\sigma = 0.1 v_{ref}$')
        do j = 1, size(hl)
            write (lab, '(a,f3.1)') 'ħ = ', hl(j)
            call plot(kl, wroot(:, j), color=colors(:, j), label='root, '//trim(lab))
            call plot(kl, wfit(:, j), color=colors(:, j), linestyle='none', &
                marker='o', markersize=5.0_dp, label='Wigner run, '//trim(lab))
            call plot(kl, wbg(:, j), color=black, linestyle='--', linewidth=1.2_dp)
        end do
        call xlabel('$k v_{ref}/\omega_p$'); call ylabel('$\omega/\omega_p$')
        call legend(loc='upper left')
        call subplot(2, 2, 2)
        call title('Two-stream growth, beams $\pm v_{ref}$, $\sigma = 0.2 v_{ref}$')
        do j = 1, size(hl)
            write (lab, '(a,f3.1)') 'ħ = ', hl(j)
            call plot(pack(kt, groot(:, j) > 0), pack(groot(:, j), groot(:, j) > 0), &
                color=colors(:, j), label='root, '//trim(lab))
            call plot(pack(kt, gfit(:, j) >= 0), pack(gfit(:, j), gfit(:, j) >= 0), &
                color=colors(:, j), linestyle='none', marker='o', markersize=5.0_dp)
            call plot(pack(kt, gcold(:, j) > 0), pack(gcold(:, j), gcold(:, j) > 0), &
                color=black, linestyle='--', linewidth=1.2_dp)
        end do
        call xlabel('$k v_{ref}/\omega_p$'); call ylabel('$\gamma/\omega_p$')
        call ylim(0.0_dp, 0.42_dp)
        call legend(loc='upper right')
        call subplot(2, 2, 3)
        call title('Wigner to Vlasov as ħ to 0')
        call plot(hb, errhb(:, 1), color=blue, marker='o', label='two-stream, t = 20')
        call plot(hb, errhb(:, 2), color=verm, marker='s', label='Landau, t = 10')
        call plot(hb, errhb(1, 1)*(hb/hb(1))**2, color=black, linestyle='--', &
            linewidth=1.2_dp, label='slope 2')
        call set_xscale('log'); call set_yscale('log')
        call xlabel('ħ $\omega_p/(m v_{ref}^2)$'); call ylabel('$||f_W - f_V||_2$')
        call legend(loc='upper left')
        call subplot(2, 2, 4)
        call title('Negative mass, two-stream (not clipped)')
        call plot(times, max(neg(:, 1), 1e-16_dp), color=blue, label='Wigner ħ = 1')
        call plot(times, max(neg(:, 2), 1e-16_dp), color=verm, label='Wigner ħ = 0.5')
        call plot(times, max(neg(:, 3), 1e-16_dp), color=grey, linestyle=':', &
            label='Vlasov (numerical indicator)')
        call set_yscale('log'); call ylim(1e-16_dp, 10.0_dp)
        call xlim(times(1), times(size(times)))
        call xlabel('$\omega_p t$'); call ylabel('$\int\int \max(-f, 0) dx dv$')
        call legend(loc='lower right')
        call suptitle('1X1V Wigner-Poisson vs Vlasov-Poisson; markers: split-step '// &
            'runs, lines: Z-function roots, dashed: Bohm-Gross / cold quartic')
        call tight_layout()
        call savefig('output/wigner-dispersion.png')
    end subroutine

    !> Movie f_W | f_V with a diverging map centred at zero (symmetric limits),
    !> plus first-harmonic amplitude and negative-mass histories.
    subroutine plot_movie(fw, fv, times, neg, amp, length, vmax, frames, hold, snapshot)
        real(dp), intent(in) :: fw(:, :, :), fv(:, :, :), times(:), neg(:, :)
        real(dp), intent(in) :: amp(:, :), length, vmax
        integer, intent(in) :: frames, hold, snapshot
        type(animation_t) :: animation
        integer :: j, status
        mfw = fw; mfv = fv; mt = times; mneg = max(neg, 1e-16_dp)
        mamp = max(amp, 1e-16_dp); mframes = frames
        mrange = max(maxval(abs(fw)), maxval(abs(fv)))
        mnegmax = max(maxval(-fw), 1e-12_dp)
        ! First frame at which the Vlasov negative mass (pure numerics) exceeds
        ! 1e-6: beyond it the filamentation is no longer resolved by the grid.
        mresolved = huge(1.0_dp)
        do j = 1, size(neg, 1)
            if (neg(j, 3) > 1e-6_dp) then
                mresolved = times(j); exit
            end if
        end do
        allocate (mx(size(fw, 2)), mv(size(fw, 1)))
        do j = 1, size(mx)
            mx(j) = length*(j - 1)/(size(mx) - 1)
        end do
        do j = 1, size(mv)
            mv(j) = -vmax + 2*vmax*(j - 1)/(size(mv) - 1)
        end do
        call figure(figsize=[13.2_dp, 8.8_dp])
        canvas => get_global_figure()
        still_frame = .true.
        call draw(snapshot)
        still_frame = .false.
        call savefig('output/wigner-two-stream.png')
        animation = animation_func(update, frames=frames + hold, interval=42, fig=canvas)
        call save_animation(animation, 'output/wigner-two-stream.mp4', fps=24, &
                            status=status)
        if (status /= 0) write (*, '(a)') 'warning: movie encoding failed'
    end subroutine

    subroutine update(number)
        integer, intent(in) :: number
        call canvas%clear()
        call draw(min(number, mframes))
        call canvas%set_rendered(.false.)
    end subroutine

    subroutine draw(n)
        integer, intent(in) :: n
        character(64) :: heading
        call subplot(2, 2, 1)
        write (heading, '(a,f5.1)') 'Wigner $f_W$, ħ = 1, $\omega_p t$ = ', mt(n)
        call title(trim(heading))
        call pcolormesh(mx, mv, mfw(:, :, n), cmap='coolwarm', vmin=-mrange, vmax=mrange)
        call colorbar(label='$f_W v_{ref}/n_0$')
        call xlim(mx(1), mx(size(mx))); call ylim(mv(1), mv(size(mv)))
        call xlabel('$x \omega_p/v_{ref}$'); call ylabel('$v/v_{ref}$')
        call subplot(2, 2, 2)
        write (heading, '(a,f5.1)') 'Vlasov $f$, ħ = 0, $\omega_p t$ = ', mt(n)
        call title(trim(heading))
        call pcolormesh(mx, mv, mfv(:, :, n), cmap='coolwarm', vmin=-mrange, vmax=mrange)
        call colorbar(label='$f v_{ref}/n_0$')
        call xlim(mx(1), mx(size(mx))); call ylim(mv(1), mv(size(mv)))
        call xlabel('$x \omega_p/v_{ref}$'); call ylabel('$v/v_{ref}$')
        call subplot(2, 2, 3)
        call title('Negative part of $f_W$, ħ = 1')
        call pcolormesh(mx, mv, max(-mfw(:, :, n), 0.0_dp), cmap='viridis', &
                        vmin=0.0_dp, vmax=mnegmax)
        call colorbar(label='$-f_W v_{ref}/n_0$ where $f_W < 0$')
        call xlim(mx(1), mx(size(mx))); call ylim(mv(1), mv(size(mv)))
        call xlabel('$x \omega_p/v_{ref}$'); call ylabel('$v/v_{ref}$')
        call subplot(2, 2, 4)
        call title('Mode growth, saturation and negative mass')
        call plot(mt(:n), mamp(:n, 1), color=blue, label='$|n_1|/n_0$, ħ = 1')
        call plot(mt(:n), mamp(:n, 3), color=verm, linestyle='--', &
                  label='$|n_1|/n_0$, Vlasov')
        call plot(mt(:n), mneg(:n, 1), color=green, linestyle='-.', &
                  label='negative mass, ħ = 1')
        call plot(mt(:n), mneg(:n, 3), color=grey, linestyle=':', &
                  label='negative mass, Vlasov')
        if (mresolved < merge(mt(n), mt(mframes), still_frame)) call plot([mresolved, mresolved], [1e-16_dp, 10.0_dp], &
            color=black, linestyle='--', linewidth=1.0_dp, &
            label='Vlasov under-resolved')
        call set_yscale('log')
        call xlim(0.0_dp, merge(mt(n), mt(mframes), still_frame)); call ylim(1e-16_dp, 10.0_dp)
        call xlabel('$\omega_p t$')
        call ylabel('$|n_1|/n_0$ and negative mass')
        call legend(loc='lower right')
        call suptitle('Two-stream instability: linear growth, saturation and a '// &
            'phase-space vortex; quantum $f_W$ turns negative where the beams overlap')
        call tight_layout()
    end subroutine
end module kin6d_wigner_plots

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
    use kin6d_result_io, only: run_t, result_file_t
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp)
    integer, parameter :: nlk = 6, nts = 15, nh = 3, nhb = 5, frames = 288
    integer, parameter :: mx = 128, mv = 512, hold = 24
    real(dp) :: hl(nh) = [0.0_dp, 1.0_dp, 2.0_dp], ht(nh) = [0.0_dp, 1.0_dp, 2.0_dp]
    real(dp) :: kl(nlk), wl_fit(nlk, nh), wl_root(nlk, nh), wl_bg(nlk, nh)
    real(dp) :: kt(nts), g_fit(nts, nh), g_root(nts, nh), g_cold(nts, nh)
    real(dp) :: hb(nhb), err_hb(nhb, 2), rate(2), errors(3, 3), cost(6), t0, t1
    real(dp) :: neg(frames, 3), amp(frames, 3), times(frames)
    real(dp), allocatable :: fw(:, :, :), fv(:, :, :)
    real(dp) :: ledger(3, 2), dgamma(3), u, bg
    complex(dp) :: root, droot, fit, guess
    type(beam_t) :: cool(1), streams(2)
    integer :: i, j, unit, status, config_unit, ios
    logical :: ok
    type(run_t) :: run
    type(result_file_t) :: result
    character(256) :: msg
    !> Scan time step, movie end time and movie quantum parameter.
    real(dp) :: dt = 0.05_dp, t_movie = 45.0_dp, hbar_movie = 1.0_dp
    namelist /wigner/ dt, t_movie, hbar_movie

    call run%start('wigner_demo')
    if (run%argument_count() > 0) error stop 'usage: wigner_demo [--config FILE]'
    ios = 0; msg = ''
    call run%open_config(config_unit)
    if (run%has_config()) read (config_unit, nml=wigner, iostat=ios, iomsg=msg)
    call run%close_config(config_unit, 'wigner', ios, msg)
    if (.not. (dt > 0 .and. dt <= 0.1_dp .and. t_movie > 0 .and. hbar_movie > 0)) &
        error stop 'wigner: invalid configuration'
    call run%open_resolved(config_unit)
    write (config_unit, nml=wigner)
    call run%close_resolved(config_unit)
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
            ! Window after ballistic phase mixing exp(-(k sigma t)^2/2) has decayed.
            call run_fit([1.0_dp], [0.0_dp], [0.1_dp], kl(i), hl(j), &
                         1.5_dp + hl(j)*kl(i), 512, 60/kl(i), 60/kl(i) + 60, &
                         1e-9_dp, fit)
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
                    kt(i), ht(j), 4.0_dp + ht(j)*kt(i)/2, 256, 24.0_dp, 80.0_dp, &
                    1e-12_dp, fit)
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
    call write_result()
    call plot_movie(fw, fv, times, neg, amp, 2*pi/0.6_dp, 4.0_dp, frames, hold, &
                    nint(min(32.0_dp, t_movie)/t_movie*(frames - 1)) + 1)
contains
    !> Run a linear seeded mode and fit its dominant complex frequency on the
    !> window [max(ta, tend - (tb - ta)), tend]; tend = tb or the first time the
    !> first-harmonic amplitude exceeds 1e-5 (linear regime).
    subroutine run_fit(w, drift, sigma, k, hbar, vmax, nv, ta, tb, alpha, omega)
        real(dp), intent(in) :: w(:), drift(:), sigma(:), k, hbar, vmax, ta, tb, alpha
        integer, intent(in) :: nv
        complex(dp), intent(out) :: omega
        type(wigner_grid_t) :: g
        type(wigner_work_t) :: wk
        type(wigner_state_t) :: s
        type(wigner_diagnostics_t) :: d
        real(dp), allocatable :: signal(:), at(:)
        integer :: n, m, stride, first
        logical :: fit_ok
        call wigner_grid_init(g, 16, nv, 2*pi/k, vmax)
        call wigner_work_init(g, hbar, wk)
        call wigner_initial_beams(g, w, drift, sigma, alpha, 1, hbar, s, .false.)
        stride = max(1, nint(0.25_dp/dt)); m = 0
        allocate (signal(nint(tb/dt)/stride + 1), at(nint(tb/dt)/stride + 1))
        do n = 1, nint(tb/dt)
            call wigner_step(g, wk, s, dt)
            if (mod(n, stride) == 0) then
                call wigner_diagnostics(g, wk, s, d)
                m = m + 1; signal(m) = real(d%mode, dp); at(m) = s%time
                if (abs(d%mode) > 1e-5_dp .and. s%time > ta) exit
            end if
        end do
        first = 1
        do while (at(first) < max(ta, at(m) - (tb - ta)) - 1e-9_dp)
            first = first + 1
        end do
        call wigner_fit_mode(signal(first:m), stride*dt, omega, fit_ok)
        if (.not. fit_ok) omega = cmplx(-1, -1, dp)
    end subroutine

    !> d gamma/d hbar from the discrete forward tangent: slope of
    !> d log|n_1|/d hbar over the late linear window [24, 40].
    subroutine growth_tangent(k, hbar, slope)
        real(dp), intent(in) :: k, hbar
        real(dp), intent(out) :: slope
        type(wigner_grid_t) :: g
        type(wigner_work_t) :: wk
        type(wigner_state_t) :: s
        type(wigner_diagnostics_t) :: d
        real(dp) :: la
        integer :: n
        call wigner_grid_init(g, 16, 256, 2*pi/k, 4.0_dp + hbar*k/2)
        call wigner_work_init(g, hbar, wk)
        call wigner_initial_beams(g, [0.5_dp, 0.5_dp], [1.0_dp, -1.0_dp], &
            [0.2_dp, 0.2_dp], 1e-12_dp, 1, hbar, s, .true.)
        la = 0
        do n = 1, nint(40/dt)
            call wigner_step(g, wk, s, dt)
            if (n == nint(24/dt)) then
                call wigner_diagnostics(g, wk, s, d)
                la = real(conjg(d%mode)*d%dmode, dp)/abs(d%mode)**2
            end if
        end do
        call wigner_diagnostics(g, wk, s, d)
        slope = (real(conjg(d%mode)*d%dmode, dp)/abs(d%mode)**2 - la)/16
    end subroutine

    !> Nonlinear final state f(x, v, tend) on a fixed grid.
    subroutine final_state(case, hbar, step, nx, nv, f, diag, d0)
        integer, intent(in) :: case, nx, nv
        real(dp), intent(in) :: hbar, step
        real(dp), intent(out) :: f(:, :)
        type(wigner_diagnostics_t), intent(out) :: diag, d0
        type(wigner_grid_t) :: g
        type(wigner_work_t) :: wk
        type(wigner_state_t) :: s
        real(dp) :: tend
        integer :: n
        if (case == 1) then
            call wigner_grid_init(g, nx, nv, 2*pi/0.6_dp, 4.0_dp)
            call wigner_work_init(g, hbar, wk)
            call wigner_initial_beams(g, [0.5_dp, 0.5_dp], [1.0_dp, -1.0_dp], &
                [0.2_dp, 0.2_dp], 1e-3_dp, 1, hbar, s, .false.)
            tend = 20
        else
            call wigner_grid_init(g, nx, nv, 4*pi, 8.0_dp)
            call wigner_work_init(g, hbar, wk)
            call wigner_initial_beams(g, [1.0_dp], [0.0_dp], [1.0_dp], 0.05_dp, 1, &
                hbar, s, .false.)
            tend = 10
        end if
        call wigner_diagnostics(g, wk, s, d0)
        do n = 1, nint(tend/step)
            call wigner_step(g, wk, s, step)
        end do
        call wigner_diagnostics(g, wk, s, diag)
        f = s%f*sqrt(g%dx*g%dv)
    end subroutine

    subroutine hbar_limit(errors, rates)
        real(dp), intent(out) :: errors(:, :), rates(2)
        real(dp) :: fref(64, 256), f(64, 256)
        type(wigner_diagnostics_t) :: d, d0
        integer :: c, i
        do c = 1, 2
            call final_state(c, 0.0_dp, 0.025_dp, 64, 256, fref, d, d0)
            do i = 1, size(hb)
                call final_state(c, hb(i), 0.025_dp, 64, 256, f, d, d0)
                errors(i, c) = sqrt(sum((f - fref)**2))
            end do
            rates(c) = log(errors(size(hb) - 1, c)/errors(size(hb), c))/log(2.0_dp)
        end do
    end subroutine

    !> Indicators (not bounds) for |n_1(t=20)| of the hbar = 1 two-stream:
    !> rows dt (0.1, 0.05, 0.025 vs 0.0125), nv (128, 256, 512 vs 1024),
    !> nx (16, 32, 64 vs 128); differences to the finest member.
    subroutine error_ledger(errors)
        real(dp), intent(out) :: errors(3, 3)
        real(dp) :: ref
        integer :: i
        ref = amplitude(0.0125_dp, 64, 256)
        do i = 1, 3
            errors(1, i) = abs(amplitude(0.2_dp/2**i, 64, 256) - ref)
        end do
        ref = amplitude(0.025_dp, 64, 1024)
        do i = 1, 3
            errors(2, i) = abs(amplitude(0.025_dp, 64, 64*2**i) - ref)
        end do
        ref = amplitude(0.025_dp, 128, 256)
        do i = 1, 3
            errors(3, i) = abs(amplitude(0.025_dp, 8*2**i, 256) - ref)
        end do
    end subroutine

    real(dp) function amplitude(step, nx, nv) result(a)
        real(dp), intent(in) :: step
        integer, intent(in) :: nx, nv
        real(dp), allocatable :: f(:, :)
        type(wigner_diagnostics_t) :: d, d0
        allocate (f(nx, nv))
        call final_state(1, 1.0_dp, step, nx, nv, f, d, d0)
        a = abs(d%mode)
    end function

    subroutine movie_runs()
        type(wigner_grid_t) :: g
        type(wigner_work_t) :: wk
        type(wigner_state_t) :: s
        type(wigner_diagnostics_t) :: d, d0
        real(dp) :: hbars(3), step
        integer :: c, n, frame, sub
        hbars = [hbar_movie, hbar_movie/2, 0.0_dp]
        ! Substeps per frame so that the split step stays at or below 0.05.
        sub = ceiling(t_movie/(frames - 1)/0.05_dp - 1e-12_dp)
        step = t_movie/(sub*(frames - 1))
        do c = 1, 3
            call wigner_grid_init(g, mx, mv, 2*pi/0.6_dp, 4.0_dp)
            call wigner_work_init(g, hbars(c), wk)
            call wigner_initial_beams(g, [0.5_dp, 0.5_dp], [1.0_dp, -1.0_dp], &
                [0.2_dp, 0.2_dp], 1e-3_dp, 1, hbars(c), s, .false.)
            call wigner_diagnostics(g, wk, s, d0)
            frame = 0
            do n = 0, sub*(frames - 1)
                if (n > 0) call wigner_step(g, wk, s, step)
                if (mod(n, sub) /= 0) cycle
                frame = frame + 1
                call wigner_diagnostics(g, wk, s, d)
                times(frame) = s%time; neg(frame, c) = d%negative_mass
                amp(frame, c) = abs(d%mode)
                if (c == 1) call store(fw(:, :, frame), s%f)
                if (c == 3) call store(fv(:, :, frame), s%f)
            end do
            if (c /= 2) then
                j = merge(1, 2, c == 1)
                ledger(1, j) = d%mass - d0%mass
                ledger(2, j) = d%momentum - d0%momentum
                ledger(3, j) = (d%energy - d0%energy)/d0%energy
            end if
        end do
    end subroutine

    !> Image (v, x) with periodic duplicates of the first row/column.
    subroutine store(image, f)
        real(dp), intent(out) :: image(:, :)
        real(dp), intent(in) :: f(:, :)
        image(1:mv, 1:mx) = transpose(f)
        image(mv + 1, 1:mx) = f(:, 1)
        image(:, mx + 1) = image(:, 1)
    end subroutine

    subroutine write_summary()
        open (newunit=unit, file='output/wigner-summary.txt', status='replace')
        write (unit, '(a)') 'Q7 Wigner-Poisson vs Vlasov-Poisson (units omega_p, v_ref, m)'
        write (unit, '(a,es10.3)') 'Langmuir max |omega_fit - omega_root|: ', &
            maxval(abs(wl_fit - wl_root))
        write (unit, '(a,es10.3)') 'Langmuir max |Bohm-Gross/root - 1|:   ', &
            maxval(abs(wl_bg/wl_root - 1))
        write (unit, '(a,es10.3)') 'two-stream max |gamma_fit - gamma_root|: ', &
            maxval(abs(g_fit - g_root), mask=g_fit >= 0)
        write (unit, '(a,3es12.4)') 'dgamma/dhbar k=0.6 hbar=1 [tangent, implicit, FD]:', &
            dgamma
        write (unit, '(a,5es10.3)') 'hbar->0 two-stream L2 errors: ', err_hb(:, 1)
        write (unit, '(a,5es10.3)') 'hbar->0 Landau L2 errors:     ', err_hb(:, 2)
        write (unit, '(a,2f7.3)') 'observed hbar rates (two-stream, Landau): ', rate
        write (unit, '(a,3es10.3)') '|n1| indicator dt 0.1/0.05/0.025: ', errors(1, :)
        write (unit, '(a,3es10.3)') '|n1| indicator nv 128/256/512:    ', errors(2, :)
        write (unit, '(a,3es10.3)') '|n1| indicator nx 16/32/64:       ', errors(3, :)
        write (unit, '(a,3es10.3)') 'Wigner ledger dmass dmomentum dE/E: ', ledger(:, 1)
        write (unit, '(a,3es10.3)') 'Vlasov ledger dmass dmomentum dE/E: ', ledger(:, 2)
        write (unit, '(a,3es10.3)') 'negative mass at t=20 (hbar 1, 0.5, 0): ', &
            neg(nint(20/t_movie*(frames - 1)) + 1, :)
        write (unit, '(a,6f8.2)') 'cpu seconds langmuir/two-stream/jvp/hbar/errors/movie:', &
            cost
        close (unit)
    end subroutine

    subroutine write_result()
        call result%create(run)
        call result%write('k_langmuir', kl, '1', 'k v_ref/omega_p', dims='k_langmuir')
        call result%write('hbar_scan', hl, '1', 'hbar omega_p/(m v_ref^2)', dims='hbar_scan')
        call result%write('omega_fit', wl_fit, '1', 'Langmuir frequency, split-step run', &
            dims='k_langmuir,hbar_scan', status='candidate; linear-prediction fit')
        call result%write('omega_root', wl_root, '1', 'Z-function dispersion root', &
            dims='k_langmuir,hbar_scan', status='independent oracle')
        call result%write('omega_bohm_gross', wl_bg, '1', 'small-k Bohm-Gross law', &
            dims='k_langmuir,hbar_scan', status='asymptotic reference')
        call result%write('k_two_stream', kt, '1', 'k v_ref/omega_p', dims='k_two_stream')
        call result%write('gamma_fit', g_fit, '1', 'growth rate, split-step run (-1: '// &
            'not run)', dims='k_two_stream,hbar_scan', status='candidate')
        call result%write('gamma_root', g_root, '1', 'growth rate, Z-function root', &
            dims='k_two_stream,hbar_scan', status='independent oracle')
        call result%write('gamma_cold', g_cold, '1', 'cold-beam quartic growth rate', &
            dims='k_two_stream,hbar_scan', status='cold-limit reference')
        call result%write('dgamma_dhbar', dgamma, '1', 'tangent, implicit, FD', &
            dims='derivative_route')
        call result%write('hbar_limit', hb, '1', 'hbar values', dims='hbar_limit')
        call result%write('wigner_vlasov_l2', err_hb, '1', '||f_W-f_V||_2 two-stream, '// &
            'Landau', dims='hbar_limit,limit_case', status='candidate')
        call result%write('error_indicators', errors, '1', '|n1| indicators rows dt,nv,nx', &
            dims='error_axis,error_level', status='indicator, not a bound')
        call result%write('time', times, '1', 'omega_p t of movie frames', dims='time')
        call result%write('negative_mass', neg, '1', 'negative mass hbar=1,0.5,0', &
            dims='time,movie_case', status='diagnostic; f not clipped')
        call result%write('mode_amplitude', amp, '1', '|n_1|/n_0', dims='time,movie_case')
        call result%write('conservation', ledger, '1', 'dmass,dmomentum,dE/E Wigner/Vlasov', &
            dims='ledger,ledger_case')
        call result%write('cpu_seconds', cost, 's', 'section costs', dims='cost')
        call result%close()
    end subroutine
end program wigner_demo
