!> Q7 part 2 oracles: Wigner-Poisson / Vlasov-Poisson split-step solver
!> against the independent Z-function dispersion relation (Landau, quantum
!> Langmuir, two-stream), exact cold-beam limits, conservation, splitting
!> order, hbar -> 0 convergence, tangent-vs-FD, and unclipped negativity.
program test_wigner
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_wigner_poisson
    use kin6d_wigner_dispersion, only: beam_t, dispersion_root
    use kin6d_wigner_leaves, only: wigner_bohm_gross, wigner_cold_two_stream
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp)
    type(beam_t) :: maxwell(1), cold(1), streams(2)
    complex(dp) :: fit, root, droot
    real(dp) :: bg, u, e(4), ref(32, 128), fv(32, 128), ratio, plus, minus, tangent
    real(dp) :: neg_w, neg_v, mass_change
    type(wigner_diagnostics_t) :: d
    integer :: failures, j
    logical :: ok
    failures = 0
    maxwell(1) = beam_t(1.0_dp, 0.0_dp, 1.0_dp)
    streams = [beam_t(0.5_dp, 1.0_dp, 0.2_dp), beam_t(0.5_dp, -1.0_dp, 0.2_dp)]

    ! 1. Classical Landau root vs literature (Canosa 1973 table: 1.4157-0.1534i).
    call dispersion_root(0.5_dp, 0.0_dp, maxwell, (1.4_dp, -0.15_dp), root, droot, ok)
    call check('Landau root vs literature', abs(root - (1.4157_dp, -0.1534_dp)), 1e-4_dp)
    ! 2. Simulated Landau mode, classical and quantum (hbar = 2).
    call landau(0.0_dp, fit)
    call check('Vlasov Landau vs root', abs(fit - root), 3e-4_dp)
    write (*, '(a,2f11.6,a,2f11.6)') 'Vlasov  k=0.5 fit ', fit, ' root ', root
    call dispersion_root(0.5_dp, 2.0_dp, maxwell, root, root, droot, ok)
    call landau(2.0_dp, fit)
    write (*, '(a,2f11.6,a,2f11.6)') 'Wigner  k=0.5 hbar=2 fit ', fit, ' root ', root
    call check('Wigner Landau vs quantum root', abs(fit - root), 3e-4_dp)

    ! 3. Quantum Langmuir wave in a cool plasma (sigma = 0.1, k = 1, hbar = 1).
    cold(1) = beam_t(1.0_dp, 0.0_dp, 0.1_dp)
    call wigner_bohm_gross(1.0_dp, 0.1_dp, 1.0_dp, bg)
    call dispersion_root(1.0_dp, 1.0_dp, cold, cmplx(bg, 0, dp), root, droot, ok)
    call run_fit([1.0_dp], [0.0_dp], [0.1_dp], 1.0_dp, 1.0_dp, 1.5_dp, 256, &
                 0.05_dp, 50.0_dp, 110.0_dp, fit)
    write (*, '(a,f11.6,a,f11.6,a,f11.6)') 'cool quantum Langmuir fit ', real(fit), &
        ' root ', real(root), ' Bohm-Gross ', bg
    call check('cool quantum Langmuir vs root', abs(fit - root), 3e-4_dp)
    call check('Bohm-Gross asymptotic (1%)', abs(bg/real(root) - 1), 1e-2_dp)
    ! Exact cold limit omega^2 = 1 + hbar^2 k^4/4 of the oracle itself.
    cold(1)%sigma = 0.02_dp
    call dispersion_root(1.0_dp, 1.0_dp, cold, (1.118_dp, 0.0_dp), root, droot, ok)
    call wigner_bohm_gross(1.0_dp, 0.02_dp, 1.0_dp, bg)
    call check('cold oracle vs Bohm-Gross sigma=0.02', abs(real(root) - bg), 2e-4_dp)

    ! 4. Two-stream growth (k = 0.6), late linear window [24, 40].
    do j = 0, 1
        call dispersion_root(0.6_dp, 1.0_dp*j, streams, (0.0_dp, 0.3_dp), root, &
                             droot, ok)
        call run_fit([0.5_dp, 0.5_dp], [1.0_dp, -1.0_dp], [0.2_dp, 0.2_dp], 0.6_dp, &
                     1.0_dp*j, 4.0_dp, 256, 0.05_dp, 24.0_dp, 40.0_dp, fit)
        write (*, '(a,f4.1,a,f11.6,a,f11.6)') 'two-stream hbar=', 1.0_dp*j, &
            ' growth fit ', aimag(fit), ' root ', aimag(root)
        call check('two-stream growth vs root', abs(aimag(fit) - aimag(root)), 3e-4_dp)
    end do
    ! Cold-beam quartic vs oracle at small beam temperature.
    streams%sigma = 0.02_dp
    call dispersion_root(0.6_dp, 1.0_dp, streams, (0.0_dp, 0.34_dp), root, droot, ok)
    call wigner_cold_two_stream(0.6_dp, 1.0_dp, 1.0_dp, u)
    call check('cold quantum two-stream quartic', abs(aimag(root) - sqrt(-u)), 2e-3_dp)
    streams%sigma = 0.2_dp

    ! 5./8. Conservation and Strang order (Landau, alpha = 0.05, t = 10).
    call nonlinear(0.0_dp, 0.0125_dp, 10.0_dp, ref, d, mass_change)
    do j = 1, 3
        call nonlinear(0.0_dp, 0.2_dp/2**j, 10.0_dp, fv, d, mass_change)
        e(j) = sqrt(sum((fv - ref)**2))
        call check('mass conservation', abs(mass_change), 1e-12_dp)
    end do
    write (*, '(a,3es11.3)') 'Strang dt-errors (0.1, 0.05, 0.025) ', e(1:3)
    call check('Strang order 2', abs(log(e(1)/e(2))/log(2.0_dp) - 2), 0.15_dp)

    ! 6. hbar -> 0: ||f_W - f_V|| = O(hbar^2).
    call nonlinear(0.0_dp, 0.05_dp, 10.0_dp, ref, d, mass_change)
    do j = 1, 3
        call nonlinear(0.4_dp/2**(j - 1), 0.05_dp, 10.0_dp, fv, d, mass_change)
        e(j) = sqrt(sum((fv - ref)**2))
    end do
    ratio = log(e(2)/e(3))/log(2.0_dp)
    write (*, '(a,3es11.3,a,f6.3)') 'hbar 0.4/0.2/0.1 Wigner-Vlasov ', e(1:3), &
        ' rate ', ratio
    call check('hbar^2 convergence', abs(ratio - 2), 0.05_dp)

    ! 7. Discrete tangent d(mode)/d hbar vs central finite difference.
    call tangent_check(tangent, plus, minus)
    write (*, '(a,es14.6,a,es14.6)') 'JVP d mode/d hbar ', tangent, &
        ' FD ', (plus - minus)/2e-4_dp
    call check('tangent vs FD', abs(tangent - (plus - minus)/2e-4_dp)/abs(tangent), &
        1e-6_dp)

    ! 9. Negativity is retained, not clipped: nonlinear two-stream.
    call two_stream_nonlinear(1.0_dp, neg_w, mass_change)
    call two_stream_nonlinear(0.0_dp, neg_v, mass_change)
    write (*, '(a,es11.3,a,es11.3)') 'negative mass Wigner ', neg_w, ' Vlasov ', neg_v
    call check('Wigner negative mass present', -neg_w, -1e-3_dp)
    call check('mass conserved with negativity', abs(mass_change), 1e-12_dp)

    if (failures > 0) error stop 'wigner oracles failed'
    write (*, '(a)') 'wigner oracles passed'
contains
    subroutine check(label, value, tolerance)
        character(*), intent(in) :: label
        real(dp), intent(in) :: value, tolerance
        if (.not. (value <= tolerance)) then
            write (*, '(a,a,es12.4,a,es10.3)') 'FAIL ', label, value, ' > ', tolerance
            failures = failures + 1
        end if
    end subroutine

    subroutine landau(hbar, omega)
        real(dp), intent(in) :: hbar
        complex(dp), intent(out) :: omega
        call run_fit([1.0_dp], [0.0_dp], [1.0_dp], 0.5_dp, hbar, 8.0_dp, 128, &
                     0.05_dp, 8.0_dp, 35.0_dp, omega)
    end subroutine

    subroutine run_fit(w, drift, sigma, k, hbar, vmax, nv, dt, ta, tb, omega)
        real(dp), intent(in) :: w(:), drift(:), sigma(:), k, hbar, vmax, dt, ta, tb
        integer, intent(in) :: nv
        complex(dp), intent(out) :: omega
        type(wigner_grid_t) :: g
        type(wigner_work_t) :: wk
        type(wigner_state_t) :: s
        real(dp) :: signal(1000)
        integer :: n, m, stride
        logical :: fit_ok
        call wigner_grid_init(g, 16, nv, 2*pi/k, vmax)
        call wigner_work_init(g, hbar, wk)
        call wigner_initial_beams(g, w, drift, sigma, 1e-9_dp, 1, hbar, s, .false.)
        stride = nint(0.25_dp/dt); m = 0
        do n = 1, nint(tb/dt)
            call wigner_step(g, wk, s, dt)
            if (mod(n, stride) == 0 .and. s%time >= ta - 1e-9_dp) then
                call wigner_diagnostics(g, wk, s, d)
                m = m + 1; signal(m) = real(d%mode, dp)
            end if
        end do
        call wigner_fit_mode(signal(:m), stride*dt, omega, fit_ok)
        if (.not. fit_ok) failures = failures + 1
    end subroutine

    subroutine nonlinear(hbar, dt, tend, f, diag, dmass)
        real(dp), intent(in) :: hbar, dt, tend
        real(dp), intent(out) :: f(:, :), dmass
        type(wigner_diagnostics_t), intent(out) :: diag
        type(wigner_grid_t) :: g
        type(wigner_work_t) :: wk
        type(wigner_state_t) :: s
        type(wigner_diagnostics_t) :: d0
        integer :: n
        call wigner_grid_init(g, 32, 128, 4*pi, 8.0_dp)
        call wigner_work_init(g, hbar, wk)
        call wigner_initial_beams(g, [1.0_dp], [0.0_dp], [1.0_dp], 0.05_dp, 1, hbar, &
                                  s, .false.)
        call wigner_diagnostics(g, wk, s, d0)
        do n = 1, nint(tend/dt)
            call wigner_step(g, wk, s, dt)
        end do
        call wigner_diagnostics(g, wk, s, diag)
        f = s%f
        dmass = diag%mass - d0%mass
    end subroutine

    subroutine tangent_check(dmode, plus, minus)
        real(dp), intent(out) :: dmode, plus, minus
        call short(1.0_dp, .true., dmode)
        call short(1.0001_dp, .false., plus)
        call short(0.9999_dp, .false., minus)
    end subroutine

    subroutine short(hbar, with_tangent, value)
        real(dp), intent(in) :: hbar
        logical, intent(in) :: with_tangent
        real(dp), intent(out) :: value
        type(wigner_grid_t) :: g
        type(wigner_work_t) :: wk
        type(wigner_state_t) :: s
        integer :: n
        call wigner_grid_init(g, 16, 128, 2*pi/0.6_dp, 4.0_dp)
        call wigner_work_init(g, hbar, wk)
        call wigner_initial_beams(g, [0.5_dp, 0.5_dp], [1.0_dp, -1.0_dp], &
            [0.2_dp, 0.2_dp], 0.01_dp, 1, hbar, s, with_tangent)
        do n = 1, 200
            call wigner_step(g, wk, s, 0.05_dp)
        end do
        call wigner_diagnostics(g, wk, s, d)
        value = real(d%mode, dp)
        if (with_tangent) value = real(d%dmode, dp)
    end subroutine

    subroutine two_stream_nonlinear(hbar, negative, dmass)
        real(dp), intent(in) :: hbar
        real(dp), intent(out) :: negative, dmass
        type(wigner_grid_t) :: g
        type(wigner_work_t) :: wk
        type(wigner_state_t) :: s
        type(wigner_diagnostics_t) :: d0
        integer :: n
        call wigner_grid_init(g, 32, 256, 2*pi/0.6_dp, 4.0_dp)
        call wigner_work_init(g, hbar, wk)
        call wigner_initial_beams(g, [0.5_dp, 0.5_dp], [1.0_dp, -1.0_dp], &
            [0.2_dp, 0.2_dp], 0.01_dp, 1, hbar, s, .false.)
        call wigner_diagnostics(g, wk, s, d0)
        do n = 1, 600
            call wigner_step(g, wk, s, 0.05_dp)
        end do
        call wigner_diagnostics(g, wk, s, d)
        negative = d%negative_mass
        dmass = d%mass - d0%mass
    end subroutine
end program test_wigner
