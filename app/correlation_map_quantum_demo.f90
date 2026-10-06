!> H1b demo: quantum effects (diffraction, exchange, degeneracy) on pair and
!> three-body correlations against the classical parent at identical finite N,
!> periodic cube, neutralizing background, density and temperature, over a
!> (r_s, Theta) grid (docs/quantum-classical-correlations.md).
!>
!> N = 2: exact relative-coordinate quantum cell (kin6d_quantum_pair_cell)
!> vs classical quadrature of the Boltzmann factor, Q3 ideal fermions,
!> Wigner--Kirkwood + leading exchange prediction, Kelbg and Kelbg+Pauli
!> pair models, and the Theta derivative. N = 3, 4: Q4 cell, every sector
!> fully diagonalized (two cutoffs), vs classical Metropolis Monte Carlo
!> (statistical intervals), pair and connected three-body Fourier moments.
program correlation_map_quantum_demo
    use, intrinsic :: iso_fortran_env, only: dp => real64, i64 => int64
    use kin6d_result_io, only: run_t, result_file_t
    use kin6d_classical_cell
    use kin6d_quantum_pair_cell
    use kin6d_quantum_classical_map
    use kin6d_jellium_cell
    use kin6d_quantum_fock, only: fock_hamiltonian_t
    use fortplot, only: figure, subplot, plot, title, suptitle, xlabel, ylabel, &
        xlim, ylim, legend, tight_layout, savefig, pcolormesh, colorbar, xscale, &
        yscale
    use kin6d_demo_style, only: blue, vermillion, green, orange, purple, sky, &
        black, grey, num
    implicit none
    integer, parameter :: max_rs = 8, max_theta = 8, max_bins = 32
    integer, parameter :: np2 = 16, nm2 = 22, npn = 12, nmn = 28
    ! Configuration (&correlation_map_quantum).
    integer :: nrs = 5, ntheta = 5, bins = 12, ngrid = 40, nshell = 8
    integer :: sweeps = 40000, model_sweeps = 16000, equilibration = 2000
    integer :: batches = 40, cutoff_n3 = 3, cutoff_n4 = 2
    real(dp) :: rs_values(max_rs) = [0.5_dp, 1.0_dp, 2.0_dp, 4.0_dp, 8.0_dp, &
                                     0.0_dp, 0.0_dp, 0.0_dp]
    real(dp) :: theta_values(max_theta) = [0.5_dp, 1.0_dp, 2.0_dp, 4.0_dp, 10.0_dp, &
                                           0.0_dp, 0.0_dp, 0.0_dp]
    real(dp) :: tau = 36, tau_check = 25, derivative_rs = 2, derivative_theta = 2
    real(dp) :: fd_step = 1.0e-3_dp
    ! Q4 thermal rows only where the cutoff holds all but q4_missing_max of
    ! the ideal canonical occupation; per-N state weight thresholds.
    real(dp) :: q4_missing_max = 0.05_dp, threshold_n3 = 1.0e-10_dp
    real(dp) :: threshold_n4 = 1.0e-6_dp
    integer(i64) :: seed = 20261006_i64
    namelist /correlation_map_quantum/ nrs, ntheta, bins, ngrid, nshell, sweeps, &
        model_sweeps, equilibration, batches, cutoff_n3, cutoff_n4, rs_values, &
        theta_values, tau, tau_check, derivative_rs, derivative_theta, fd_step, seed, &
        q4_missing_max, threshold_n3, threshold_n4
    type(run_t) :: run
    type(result_file_t) :: result
    integer :: unit, ios, r, t, n
    character(256) :: msg
    ! pair2(bin, item, rs, theta), N = 2 items:
    !  1 g quantum, 2 g_same q, 3 g_opp q, 4 g quantum (check cutoff),
    !  5 g ideal Q3, 6 g_same ideal, 7 g_opp ideal, 8 g classical,
    !  9 g WK, 10 g exchange, 11 g WK+exchange, 12 g Kelbg, 13 g Kelbg+Pauli,
    ! 14 dg_q/dTheta, 15 dg_cl/dTheta, 16 bin centre r/a.
    real(dp), allocatable :: pair2(:, :, :, :), met2(:, :, :)
    ! pairn(bin, item, rs, theta 0:ntheta, N): 1 g q (cut c), 2 g q (c-1),
    !  3 g q ideal (Q4 basis), 4 g ideal Q3, 5 g classical MC, 6 se,
    !  7 g WK+exchange MC, 8 g Kelbg MC, 9 g Kelbg+Pauli MC, 10 g_same q,
    ! 11 g_opp q, 12 bin centre r/a.
    real(dp), allocatable :: pairn(:, :, :, :, :), metn(:, :, :, :)
    real(dp), allocatable :: deriv(:)
    real(dp) :: t0, t1
    ! Figure state shared by the plotting helpers.
    real(dp), allocatable :: xe(:), ye(:)
    real(dp) :: palette(3, 6)
    character(80) :: label

    call run%start('correlation_map_quantum_demo')
    call run%open_config(unit)
    if (run%has_config()) read (unit, nml=correlation_map_quantum, iostat=ios, iomsg=msg)
    call run%close_config(unit, 'correlation_map_quantum', ios, msg)
    if (nrs < 1 .or. nrs > max_rs .or. any(rs_values(1:nrs) <= 0)) error stop 'rs list'
    if (ntheta < 1 .or. ntheta > max_theta .or. any(theta_values(1:ntheta) <= 0)) &
        error stop 'theta list'
    if (bins < 2 .or. bins > max_bins .or. batches < 4) error stop 'bins/batches'
    if (cutoff_n3 < 2 .or. cutoff_n3 > 3 .or. cutoff_n4 < 2 .or. cutoff_n4 > 3) &
        error stop 'Q4 cutoffs (2..3)'
    call run%open_resolved(unit)
    write (unit, nml=correlation_map_quantum)
    call run%close_resolved(unit)
    call execute_command_line('mkdir -p output', exitstat=ios)
    if (ios /= 0) error stop 'cannot create output'
    allocate (pair2(bins, np2, nrs, ntheta), met2(nm2, nrs, ntheta))
    allocate (pairn(bins, npn, nrs, 0:ntheta, 3:4), metn(nmn, nrs, 0:ntheta, 3:4))
    allocate (deriv(8))
    pair2 = 0; met2 = 0; pairn = 0; metn = 0; deriv = 0
    open (newunit=unit, file='output/qcmap-cost.csv', status='replace')
    write (unit, '(a)') 'stage,N,rs,seconds'
    do r = 1, nrs
        call cpu_time(t0)
        call pair_case(r)
        call cpu_time(t1)
        write (unit, '(a,",",i0,",",f6.2,",",es11.4)') 'N2', 2, rs_values(r), t1 - t0
        flush (unit)
    end do
    call derivative_case()
    do n = 3, 4
        do r = 1, nrs
            call cpu_time(t0)
            call many_case(n, r)
            call cpu_time(t1)
            write (unit, '(a,",",i0,",",f6.2,",",es11.4)') 'Q4+MC', n, rs_values(r), t1 - t0
            flush (unit)
        end do
    end do
    close (unit)
    call write_tables()
    call write_netcdf()
    call make_figures()
    print '(a)', 'correlation_map_quantum_demo: wrote output/qcmap-*.csv, '// &
        'output/correlation_map_quantum_demo.nc and output/qcmap-*.png'
contains
    function shell_edges(length) result(edges)
        real(dp), intent(in) :: length
        real(dp) :: edges(bins + 1)
        integer :: i
        edges = [(i*length/(2*bins), i=0, bins)]
    end function

    !> Relative RMS ||a - b|| / ||c|| guarded against a vanishing scale.
    real(dp) function ratio(edges, a, b, c) result(x)
        real(dp), intent(in) :: edges(:), a(:), b(:), c(:)
        x = shell_rms(edges, a - b)/max(shell_rms(edges, c), 1.0e-300_dp)
    end function

    !> N = 2 deterministic row at every Theta of one r_s.
    subroutine pair_case(r)
        integer, intent(in) :: r
        type(quantum_pair_cell_t) :: pc, pcc
        type(classical_cell_t) :: cc
        integer, allocatable :: q(:, :)
        real(dp), dimension(bins) :: g, gs, go, gchk, gi, gis, gio, gcl, gwk, gx, gp, &
            dgcl, gk, gkp, dgq, d1, d2, d3, d4
        real(dp) :: f2q(26), f2c(26), f2d(26), theta, rs, edges(bins + 1), lam_t
        real(dp) :: num_, den
        integer :: t, b
        logical :: ok
        rs = rs_values(r)
        call cell_q_set(q)
        call init_quantum_pair_cell(pc, rs, quantum_pair_cutoff(rs, &
            maxval(theta_values(1:ntheta)), tau, 16), ok)
        if (.not. ok) error stop 'quantum pair cell'
        call init_quantum_pair_cell(pcc, rs, quantum_pair_cutoff(rs, &
            maxval(theta_values(1:ntheta)), tau_check, 12), ok)
        if (.not. ok) error stop 'quantum pair cell (check cutoff)'
        edges = shell_edges(pc%length)
        do t = 1, ntheta
            theta = theta_values(t)
            call quantum_pair_thermal(pc, theta, edges, q, g, gs, go, f2q, ok, dgq)
            if (.not. ok) error stop 'quantum pair thermal'
            call quantum_pair_thermal(pcc, theta, edges, q, gchk, d1, d2, f2d, ok)
            call ideal_cell_pair(2, theta, pc%length, edges, gi, gis, gio, ok)
            if (.not. ok) error stop 'ideal pair'
            call init_classical_cell(cc, 2, rs, theta, model_coulomb, ok)
            call classical_pair_quadrature(cc, edges, ngrid, nshell, gcl, f2c, gwk, gx, &
                                           gp, dgcl, ok)
            if (.not. ok) error stop 'classical quadrature'
            dgcl = dgcl*(-cc%beta/theta)
            lam_t = sqrt(2*acos(-1.0_dp))*cc%lambda
            met2(1, r, t) = cc%gamma; met2(2, r, t) = cc%theta
            met2(3, r, t) = lam_t/rs
            call init_classical_cell(cc, 2, rs, theta, model_kelbg, ok)
            call classical_pair_quadrature(cc, edges, ngrid, nshell, gk, f2d, d1, d2, &
                                           d3, d4, ok)
            call init_classical_cell(cc, 2, rs, theta, model_kelbg_pauli, ok)
            call classical_pair_quadrature(cc, edges, ngrid, nshell, gkp, f2d, d1, d2, &
                                           d3, d4, ok)
            if (.not. ok) error stop 'Kelbg quadrature'
            pair2(:, 1, r, t) = g; pair2(:, 2, r, t) = gs; pair2(:, 3, r, t) = go
            pair2(:, 4, r, t) = gchk; pair2(:, 5, r, t) = gi; pair2(:, 6, r, t) = gis
            pair2(:, 7, r, t) = gio; pair2(:, 8, r, t) = gcl; pair2(:, 9, r, t) = gwk
            pair2(:, 10, r, t) = gx; pair2(:, 11, r, t) = gp; pair2(:, 12, r, t) = gk
            pair2(:, 13, r, t) = gkp; pair2(:, 14, r, t) = dgq; pair2(:, 15, r, t) = dgcl
            pair2(:, 16, r, t) = (edges(1:bins) + edges(2:))/(2*rs)
            met2(4, r, t) = shell_rms(edges, g - gcl)
            met2(5, r, t) = shell_rms(edges, gi - 1)
            met2(6, r, t) = shell_rms(edges, (g - gi) - (gcl - 1))
            met2(7, r, t) = shell_rms(edges, gs - gcl/2)
            met2(8, r, t) = shell_rms(edges, go - gcl/2)
            met2(9, r, t) = shell_rms(edges, g - gchk)
            met2(10, r, t) = ratio(edges, gp, gcl, g - gcl)
            met2(11, r, t) = ratio(edges, gp, g, g - gcl)
            met2(12, r, t) = ratio(edges, gwk, g, g - gcl)
            met2(13, r, t) = ratio(edges, gx, g, g - gcl)
            met2(14, r, t) = ratio(edges, gk, g, g - gcl)
            met2(15, r, t) = ratio(edges, gkp, g, g - gcl)
            met2(16, r, t) = sqrt(sum(f2q**2)/26); met2(17, r, t) = sqrt(sum(f2c**2)/26)
            met2(18, r, t) = sqrt(sum((f2q - f2c)**2)/26)
            met2(19, r, t) = pc%kmax2
            met2(20, r, t) = g(1); met2(21, r, t) = gcl(1)
            ! Quantum-classical difference outside the thermal wavelength.
            num_ = 0; den = 0
            do b = 1, bins
                if (edges(b) < lam_t) cycle
                num_ = num_ + (edges(b + 1)**3 - edges(b)**3)*(g(b) - gcl(b))**2
                den = den + edges(b + 1)**3 - edges(b)**3
            end do
            met2(22, r, t) = 0
            if (den > 0) met2(22, r, t) = sqrt(num_/den)
        end do
    end subroutine

    !> d||g_q - g_cl||/dTheta at one point: exact covariance derivatives of
    !> both sides vs central differences of the deterministic pipeline.
    subroutine derivative_case()
        type(quantum_pair_cell_t) :: pc
        type(classical_cell_t) :: cc
        integer, allocatable :: q(:, :)
        real(dp), dimension(bins) :: g, gs, go, gcl, dgq, dgcl, d1, d2, d3
        real(dp) :: edges(bins + 1), f2(26), dist(-1:1), dd, h, theta, vol(bins)
        integer :: k
        logical :: ok
        call cell_q_set(q)
        call init_quantum_pair_cell(pc, derivative_rs, quantum_pair_cutoff( &
            derivative_rs, 2*derivative_theta, tau, 16), ok)
        edges = shell_edges(pc%length)
        vol = edges(2:)**3 - edges(1:bins)**3
        h = fd_step*derivative_theta
        do k = -1, 1
            theta = derivative_theta + k*h
            call quantum_pair_thermal(pc, theta, edges, q, g, gs, go, f2, ok, dgq)
            call init_classical_cell(cc, 2, derivative_rs, theta, model_coulomb, ok)
            call classical_pair_quadrature(cc, edges, ngrid, nshell, gcl, f2, d1, d2, &
                                           d3, dgcl, ok)
            if (.not. ok) error stop 'derivative case'
            dist(k) = shell_rms(edges, g - gcl)
            if (k == 0) then
                dgcl = dgcl*(-cc%beta/theta)
                dd = sum(vol*(g - gcl)*(dgq - dgcl))/(sum(vol)*dist(0))
                deriv(1) = dist(0); deriv(2) = dd
                deriv(5) = maxval(abs(dgq)); deriv(6) = maxval(abs(dgcl))
            end if
        end do
        deriv(3) = (dist(1) - dist(-1))/(2*h)
        deriv(4) = abs(deriv(2) - deriv(3))/abs(deriv(3))
        deriv(7) = derivative_rs; deriv(8) = derivative_theta
        print '(a,f5.2,a,f5.2,a,es12.5,a,es12.5,a,es9.2)', 'dD/dTheta at rs=', &
            derivative_rs, ' Theta=', derivative_theta, ' analytic ', dd, ' FD ', &
            deriv(3), ' rel ', deriv(4)
    end subroutine

    !> N = 3, 4 at one r_s: Q4 cells at cutoffs c-1 and c (all sectors dense),
    !> the Q4 ideal cell (g = 0) at c, Q3 ideal pair, classical Monte Carlo.
    subroutine many_case(n, r)
        integer, intent(in) :: n, r
        type(jellium_cell_t) :: cell(3)
        type(jellium_sector_t) :: sector
        type(jellium_spectrum_t) :: spec(3)
        type(fock_hamiltonian_t) :: model
        type(classical_cell_t) :: cc
        type(classical_mc_t) :: mc, mk, mkp
        integer, allocatable :: modes(:, :), q(:, :), triple(:, :)
        real(dp), allocatable :: w(:), f2(:, :), f3(:, :)
        real(dp) :: edges(bins + 1), gsum(bins, 3), gsame(bins, 3), gopp(bins, 3)
        real(dp) :: gi(bins), gis(bins), gio(bins), trunc(3), ground, theta, rs, dtot
        real(dp) :: missing, omitted(3)
        logical :: quantum
        integer :: k, t, c, cmax
        logical :: ok
        rs = rs_values(r)
        cmax = merge(cutoff_n3, cutoff_n4, n == 3)
        call cell_q_set(q); call cell_triple_set(q, triple)
        allocate (f2(size(q, 2), 3), f3(size(triple, 2), 3))
        do k = 1, 3
            c = merge(cmax - 1, cmax, k == 1)
            call jellium_modes('spherical', c, modes)
            if (k < 3) then
                call init_jellium_cell(cell(k), n, modes, rs, ok)
            else
                call init_jellium_cell(cell(k), n, modes, rs, ok, 0.0_dp)
            end if
            if (.not. ok) error stop 'jellium cell'
            call jellium_sector(cell(k), n/2, reshape([0, 0, 0], [3, 1]), sector, ok)
            call jellium_hamiltonian(cell(k), sector%basis, model, ok)
            if (.not. ok) error stop 'jellium model'
            call jellium_full_spectrum(cell(k), model, huge(1.0_dp), spec(k), ok)
            if (.not. ok) error stop 'full spectrum'
        end do
        edges = shell_edges(cell(2)%length)
        do t = 0, ntheta
            theta = 0
            if (t > 0) theta = theta_values(t)
            missing = ideal_missing_fraction(n, theta, cmax)
            metn(27, r, t, n) = missing
            quantum = missing <= q4_missing_max
            if (quantum) then
                do k = 1, 3
                    if (allocated(w)) deallocate (w)
                    allocate (w(spec(k)%count))
                    call jellium_thermal_weights(cell(k), spec(k), theta, w, trunc(k), &
                                                 ground)
                    call jellium_diagonal_moments(cell(k), spec(k), w, edges, q, triple, &
                        gsame(:, k), gopp(:, k), f2(:, k), f3(:, k), ok, &
                        merge(threshold_n3, threshold_n4, n == 3), omitted(k))
                    if (.not. ok) error stop 'diagonal moments'
                    gsum(:, k) = (gsame(:, k) + gopp(:, k))*n/(2*(n - 1.0_dp))
                end do
                pairn(:, 1, r, t, n) = gsum(:, 2); pairn(:, 2, r, t, n) = gsum(:, 1)
                pairn(:, 3, r, t, n) = gsum(:, 3); pairn(:, 10, r, t, n) = gsame(:, 2)
                pairn(:, 11, r, t, n) = gopp(:, 2)
                metn(10, r, t, n) = trunc(2); metn(28, r, t, n) = maxval(omitted)
                metn(11, r, t, n) = sqrt(sum(f2(:, 2)**2)/size(f2, 1))
                metn(14, r, t, n) = sqrt(sum(f3(:, 2)**2)/size(f3, 1))
                metn(15, r, t, n) = sqrt(sum(f3(:, 1)**2)/size(f3, 1))
                metn(16, r, t, n) = sqrt(sum(f3(:, 3)**2)/size(f3, 1))
                metn(8, r, t, n) = shell_rms(edges, gsum(:, 2) - gsum(:, 1))
            else
                ! Basis saturated: no quantum row (sentinel -1).
                metn([3, 5, 6, 7, 8, 9, 11, 14, 15, 16, 19, 20, 21, 22], r, t, n) = -1
            end if
            pairn(:, 12, r, t, n) = (edges(1:bins) + edges(2:))/(2*rs)
            if (t == 0) cycle
            call ideal_cell_pair(n, theta, cell(2)%length, edges, gi, gis, gio, ok)
            pairn(:, 4, r, t, n) = gi
            call init_classical_cell(cc, n, rs, theta, model_coulomb, ok)
            call classical_cell_mc(cc, edges, sweeps, equilibration, batches, &
                                   seed + 1000*r + 10*t + n, mc, ok, .true.)
            if (.not. ok) error stop 'classical MC'
            call init_classical_cell(cc, n, rs, theta, model_kelbg, ok)
            call classical_cell_mc(cc, edges, model_sweeps, equilibration, batches, &
                                   seed + 2000*r + 10*t + n, mk, ok)
            call init_classical_cell(cc, n, rs, theta, model_kelbg_pauli, ok)
            call classical_cell_mc(cc, edges, model_sweeps, equilibration, batches, &
                                   seed + 3000*r + 10*t + n, mkp, ok)
            if (.not. ok) error stop 'Kelbg MC'
            pairn(:, 5, r, t, n) = mc%g; pairn(:, 6, r, t, n) = mc%g_err
            pairn(:, 7, r, t, n) = mc%g_pred; pairn(:, 8, r, t, n) = mk%g
            pairn(:, 9, r, t, n) = mkp%g
            metn(1, r, t, n) = cc%gamma
            metn(2, r, t, n) = sqrt(2*acos(-1.0_dp))*cc%lambda/rs
            metn(4, r, t, n) = shell_rms(edges, gi - 1)
            metn(12, r, t, n) = mc%c2; metn(13, r, t, n) = mc%c2_err
            metn(17, r, t, n) = mc%c3; metn(18, r, t, n) = mc%c3_err
            metn(24, r, t, n) = mc%acceptance; metn(25, r, t, n) = maxval(mc%g_err)
            metn(26, r, t, n) = spec(2)%count
            if (.not. quantum) cycle
            dtot = shell_rms(edges, gsum(:, 2) - mc%g)
            metn(3, r, t, n) = dtot
            metn(5, r, t, n) = shell_rms(edges, (gsum(:, 2) - gsum(:, 3)) - (mc%g - 1))
            metn(6, r, t, n) = shell_rms(edges, gsame(:, 2) - mc%g*(n - 1.0_dp)/n)
            metn(7, r, t, n) = shell_rms(edges, gopp(:, 2) - mc%g*(n - 1.0_dp)/n)
            metn(9, r, t, n) = shell_rms(edges, gsum(:, 3) - gi)
            metn(19, r, t, n) = sqrt(sum((f3(:, 2) - mc%f3)**2)/size(f3, 1))
            metn(20, r, t, n) = ratio(edges, mc%g_pred, gsum(:, 2), gsum(:, 2) - mc%g)
            metn(21, r, t, n) = ratio(edges, mk%g, gsum(:, 2), gsum(:, 2) - mc%g)
            metn(22, r, t, n) = ratio(edges, mkp%g, gsum(:, 2), gsum(:, 2) - mc%g)
            ! Reference row: both basis indicators and the statistical error
            ! well below the quantum-classical difference.
            metn(23, r, t, n) = merge(1.0_dp, 0.0_dp, metn(8, r, t, n) <= 0.2_dp*dtot &
                .and. metn(9, r, t, n) <= 0.2_dp*dtot .and. trunc(2) < 1.0e-3_dp &
                .and. maxval(mc%g_err) <= 0.2_dp*dtot .and. metn(28, r, t, n) < 1.0e-3_dp)
        end do
    end subroutine

    subroutine write_tables()
        integer :: lu, r, t, n, b
        open (newunit=lu, file='output/qcmap-n2.csv', status='replace')
        write (lu, '(a)') 'rs,theta,Gamma,theta_resolved,lambdaT_over_a,D_tot,D_ideal,D_interaction,'// &
            'D_same,D_opp,D_cutoff,WK_efficiency,WK_X_residual,WK_only_residual,'// &
            'X_only_residual,Kelbg_error,KelbgPauli_error,c2_quantum,c2_classical,'// &
            'f2_difference,kmax2,g0_quantum,g0_classical,D_outside_lambdaT'
        do r = 1, nrs
            do t = 1, ntheta
                write (lu, '(f6.2,",",f6.2,22(",",es12.5))') rs_values(r), &
                    theta_values(t), met2(:, r, t)
            end do
        end do
        close (lu)
        open (newunit=lu, file='output/qcmap-n2-g.csv', status='replace')
        write (lu, '(a)') 'rs,theta,r_over_a,g_q,g_same_q,g_opp_q,g_q_check,g_ideal,'// &
            'g_same_ideal,g_opp_ideal,g_cl,g_WK,g_X,g_WK_X,g_Kelbg,g_KelbgPauli,'// &
            'dgq_dtheta,dgcl_dtheta'
        do r = 1, nrs
            do t = 1, ntheta
                do b = 1, bins
                    write (lu, '(f6.2,",",f6.2,",",f8.4,15(",",es13.6))') rs_values(r), &
                        theta_values(t), pair2(b, 16, r, t), pair2(b, 1:15, r, t)
                end do
            end do
        end do
        close (lu)
        open (newunit=lu, file='output/qcmap-many.csv', status='replace')
        write (lu, '(a)') 'N,rs,theta,Gamma,lambdaT_over_a,D_tot,D_ideal,'// &
            'D_interaction,D_same,D_opp,D_cutoff,ideal_truncation,thermal_truncation,'// &
            'c2_q,c2_cl,c2_cl_se,c3_q,c3_q_lower_cutoff,c3_q_ideal,c3_cl,c3_cl_se,'// &
            'f3_difference,WK_X_residual,Kelbg_error,KelbgPauli_error,reference,'// &
            'acceptance,g_cl_se_max,stored_states,ideal_missing_fraction,omitted_weight'
        do n = 3, 4
            do r = 1, nrs
                do t = 0, ntheta
                    write (lu, '(i0,",",f6.2,",",f6.2,28(",",es12.5))') n, rs_values(r), &
                        merge(0.0_dp, theta_values(max(1, t)), t == 0), metn(:, r, t, n)
                end do
            end do
        end do
        close (lu)
        open (newunit=lu, file='output/qcmap-derivative.csv', status='replace')
        write (lu, '(a)') 'rs,theta,D,dD_dtheta_covariance,dD_dtheta_FD,relative_error,'// &
            'max_dgq_dtheta,max_dgcl_dtheta'
        write (lu, '(f6.2,",",f6.2,6(",",es14.7))') deriv(7), deriv(8), deriv(1:6)
        close (lu)
    end subroutine

    subroutine write_netcdf()
        character(8) :: tag
        integer :: n
        call result%create(run)
        call result%attribute('model', 'N electrons in a periodic cube with uniform '// &
            'neutralizing background, bare Coulomb (Ewald), canonical at fixed (n, T); '// &
            'quantum: exact N=2 relative-coordinate cell, Q4 plane-wave cell (N=3,4); '// &
            'classical: same cell, quadrature (N=2) and Metropolis Monte Carlo (N=3,4)')
        call result%attribute('status', 'N=2 deterministic finite-basis candidate with '// &
            'cutoff indicator; N=3,4 quantum rows are finite-basis indicators unless '// &
            'reference=1; classical Monte Carlo intervals are statistical (batch means)')
        call result%write('rs', rs_values(1:nrs), 'a0', 'Wigner-Seitz radius', dims='rs')
        call result%write('theta', theta_values(1:ntheta), '1', 'kT/E_F', dims='theta')
        call result%write('pair_N2', reshape(pair2, [bins, np2, nrs*ntheta]), '1', &
                          'g items (see app header); last index rs + nrs*(theta-1)', &
                          dims='bin,pair_item,rs_theta', status='candidate')
        call result%write('metrics_N2', met2, '1', 'see output/qcmap-n2.csv header', &
                          dims='metric2,rs,theta', status='candidate')
        do n = 3, 4
            write (tag, '(a,i0)') '_N', n
            call result%write('pair'//trim(tag), reshape(pairn(:, :, :, :, n), &
                [bins, npn, nrs*(ntheta + 1)]), '1', 'g items (see app header); last '// &
                'index rs + nrs*theta, theta index 0 = ground multiplet', &
                dims='bin,pairn_item,rs_theta0', status='finite-basis indicator')
            call result%write('metrics'//trim(tag), metn(:, :, :, n), '1', &
                'see output/qcmap-many.csv header', dims='metricn,rs,theta0', &
                status='finite-basis indicator; Monte Carlo statistical')
        end do
        call result%write('theta_derivative', deriv, '1', &
            'D, dD/dTheta covariance, FD, rel error, max dgq, max dgcl, rs, theta', &
            dims='derivative_item', status='candidate')
        call result%close()
    end subroutine

    subroutine make_figures()
        integer, parameter :: pick(2, 4) = reshape([2, 1, 2, 5, 4, 2, 5, 1], [2, 4])
        integer, parameter :: wk_rs(3) = [4, 3, 2], wk_theta(3) = [5, 4, 2]
        real(dp) :: z(ntheta, nrs)
        integer :: k, r, t, n
        palette(:, 1) = blue; palette(:, 2) = vermillion; palette(:, 3) = green
        palette(:, 4) = orange; palette(:, 5) = purple; palette(:, 6) = sky
        ! 1. Pair functions at selected points (N = 2).
        call figure(figsize=[13.2_dp, 8.8_dp])
        do k = 1, 4
            r = min(pick(1, k), nrs); t = min(pick(2, k), ntheta)
            call subplot(2, 2, k)
            call title('$r_s$ = '//num(rs_values(r), 1)//', Θ = '// &
                       num(theta_values(t), 1)//', Γ = '//num(met2(1, r, t), 1))
            call xlabel('r / a'); call ylabel('g(r), N = 2')
            call plot(pair2(:, 16, r, t), pair2(:, 1, r, t), color=blue, marker='o', &
                      label='quantum (exact cell)')
            call plot(pair2(:, 16, r, t), pair2(:, 5, r, t), color=green, linestyle=':', &
                      label='ideal fermions (Q3)')
            call plot(pair2(:, 16, r, t), pair2(:, 13, r, t), color=vermillion, &
                      linestyle='-.', label='classical, Kelbg + Pauli')
            call plot(pair2(:, 16, r, t), pair2(:, 11, r, t), color=orange, &
                      linestyle='--', label='Wigner-Kirkwood + exchange')
            call plot(pair2(:, 16, r, t), pair2(:, 8, r, t), color=black, linestyle='--', &
                      label='classical, Coulomb')
            call ylim(0.0_dp, 1.6_dp)
            if (k == 1) call legend(loc='lower right')
        end do
        call suptitle('Spin-summed pair correlation of two electrons in the periodic '// &
                      'cell: quantum vs classical parent')
        call tight_layout()
        call savefig('output/qcmap-pair-g.png')
        ! Log-spaced cell edges for the (r_s, Theta) heat maps.
        allocate (xe(nrs + 1), ye(ntheta + 1))
        call log_edges(rs_values(1:nrs), xe); call log_edges(theta_values(1:ntheta), ye)
        call figure(figsize=[13.2_dp, 8.8_dp])
        do k = 1, 4
            call subplot(2, 2, k)
            do r = 1, nrs
                do t = 1, ntheta
                    z(t, r) = log10(max(1.0e-4_dp, met2(3 + k + merge(1, 0, k == 4), r, t)))
                end do
            end do
            select case (k)
            case (1); call title('total $\|g_q - g_{cl}\|$')
            case (2); call title('ideal part $\|g_{q,id} - 1\|$ (exchange, degeneracy)')
            case (3); call title('interaction part (diffraction)')
            case (4); call title('opposite-spin channel $\|g_{q,\uparrow\downarrow} - g_{cl}/2\|$')
            end select
            call heat(z, '$\log_{10}$ RMS difference', -3.0_dp, 0.0_dp)
        end do
        call suptitle('N = 2 quantum-classical pair difference; dashed: Γ = 0.1, 1, 10; '// &
                      'dotted: $\lambda_T = a$')
        call tight_layout()
        call savefig('output/qcmap-regime-pair.png')
        ! 3. Three-body connected moments (N = 3, 4).
        call figure(figsize=[13.2_dp, 8.8_dp])
        do k = 1, 4
            call subplot(2, 2, k)
            n = merge(3, 4, k == 1 .or. k == 3)
            do r = 1, nrs
                do t = 1, ntheta
                    if (k <= 2) then
                        ! Skipped (basis-saturated) quantum rows sit below vmin.
                        z(t, r) = -9
                        if (metn(14, r, t, n) >= 0) z(t, r) = log10(max(1.0e-4_dp, &
                            abs(metn(14, r, t, n) - metn(17, r, t, n))))
                    else if (k == 3) then
                        z(t, r) = log10(max(1.0e-4_dp, metn(17, r, t, n)))
                    else
                        z(t, r) = log10(max(1.0e-4_dp, metn(17, r, t, 4)))
                    end if
                end do
            end do
            if (k <= 2) call title('$|c_3^q - c_3^{cl}|$, N = '//num(real(n, dp), 0))
            if (k == 3) call title('classical $c_3$, N = 3')
            if (k == 4) call title('classical $c_3$, N = 4')
            call heat(z, '$\log_{10}$ RMS connected moment', -3.0_dp, 0.0_dp)
        end do
        call suptitle('Connected three-body Fourier moments, Q4 cell vs classical '// &
                      'Monte Carlo (darkest: Q4 basis saturated, not computed)')
        call tight_layout()
        call savefig('output/qcmap-regime-three-body.png')
        ! 4. Wigner--Kirkwood prediction vs actual difference.
        call figure(figsize=[12.8_dp, 6.4_dp])
        call subplot(1, 2, 1)
        call title('N = 2: $g_q - g_{cl}$ (solid) and leading prediction (dashed)')
        call xlabel('r / a'); call ylabel('$\Delta g$')
        do k = 1, 3
            r = min(wk_rs(k), nrs); t = min(wk_theta(k), ntheta)
            write (label, '(a,a,a,a)') '$r_s$ = ', num(rs_values(r), 1), ', Θ = ', &
                num(theta_values(t), 1)
            call plot(pair2(:, 16, r, t), pair2(:, 1, r, t) - pair2(:, 8, r, t), &
                      color=palette(:, k), marker='o', label=trim(label))
            call plot(pair2(:, 16, r, t), pair2(:, 11, r, t) - pair2(:, 8, r, t), &
                      color=palette(:, k), linestyle='--')
        end do
        call legend(loc='upper right')
        call subplot(1, 2, 2)
        call title('relative residual of the WK + exchange prediction')
        call xlabel('$\lambda_T / a$'); call ylabel('$\|\Delta_{pred} - \Delta\| / \|\Delta\|$')
        call xscale('log'); call yscale('log')
        call rs_lines(11, 'WK + exchange')
        call plot([0.5_dp, 3.0_dp], [1.0_dp, 1.0_dp], color=black, linestyle='--', &
                  label='no skill (residual = 1)')
        call legend(loc='upper left')
        call tight_layout()
        call savefig('output/qcmap-wigner-kirkwood.png')
        ! 5. Kelbg model error.
        call figure(figsize=[12.8_dp, 6.4_dp])
        do k = 1, 2
            call subplot(1, 2, k)
            if (k == 1) call title('Kelbg diffraction potential')
            if (k == 2) call title('Kelbg + spin-averaged Pauli potential')
            call xlabel('$\lambda_T / a$')
            call ylabel('$\|g_{model} - g_q\| / \|g_{cl} - g_q\|$')
            call xscale('log'); call yscale('log')
            call rs_lines(13 + k, 'model')
            call plot([0.5_dp, 3.0_dp], [1.0_dp, 1.0_dp], color=black, linestyle='--', &
                      label='classical Coulomb (ratio 1)')
            call legend(loc='upper left')
        end do
        call suptitle('Semiclassical pair-potential model error vs the exact N = 2 cell')
        call tight_layout()
        call savefig('output/qcmap-kelbg.png')
    end subroutine

    subroutine log_edges(v, e)
        real(dp), intent(in) :: v(:)
        real(dp), intent(out) :: e(:)
        integer :: i, m
        m = size(v)
        if (m == 1) then
            e = [log10(v(1)) - 0.5_dp, log10(v(1)) + 0.5_dp]; return
        end if
        do i = 2, m
            e(i) = (log10(v(i - 1)) + log10(v(i)))/2
        end do
        e(1) = 2*log10(v(1)) - e(2); e(m + 1) = 2*log10(v(m)) - e(m)
    end subroutine
    subroutine heat(field, clabel, lo, hi)
        real(dp), intent(in) :: field(:, :), lo, hi
        character(*), intent(in) :: clabel
        real(dp) :: gx(2), g
        integer :: j
        call pcolormesh(xe, ye, field, cmap='viridis', vmin=lo, vmax=hi)
        call colorbar(label=clabel)
        gx = [xe(1), xe(nrs + 1)]
        do j = -1, 1
            g = 10.0_dp**j
            call plot(gx, gx - log10(1.8416_dp*g), color=black, linestyle='--')
        end do
        call plot(gx, [1.0_dp, 1.0_dp]*log10(2*acos(-1.0_dp)/1.8416_dp), color=grey, &
                  linestyle=':')
        call xlim(xe(1), xe(nrs + 1)); call ylim(ye(1), ye(ntheta + 1))
        call xlabel('$\log_{10} r_s$'); call ylabel('$\log_{10}$ Θ')
    end subroutine
    subroutine rs_lines(item, what)
        integer, intent(in) :: item
        character(*), intent(in) :: what
        integer :: i
        do i = 1, nrs
            write (label, '(a,a,a)') what//', $r_s$ = ', num(rs_values(i), 1)
            call plot(met2(3, i, 1:ntheta), max(1.0e-4_dp, met2(item, i, 1:ntheta)), &
                      color=palette(:, 1 + mod(i - 1, 6)), marker='o', label=trim(label))
        end do
    end subroutine
end program
