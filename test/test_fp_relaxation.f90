program test_fp_relaxation
    ! One executable, one check per CTest name: maxwellian, mass, slow_rate,
    ! convergence_v, convergence_dt.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_fp_relaxation, only: fp_relax_t, fp_relax_setup, fp_relax_initial, &
        fp_relax_step, fp_relax_run, fp_relax_moments, fp_relax_cells, fp_relax_nmom, &
        fp_relax_lmax
    use kin6d_fp_relaxation_leaves, only: fp_rel_maxwell, fp_rel_gauss, &
        fp_rel_rate_at, fp_rel_nu_s
    implicit none
    ! B04 parameters (thermal-force-benchmark/cases/B04_relaxation.nml, conventions.md).
    real(dp), parameter :: amu = 1.66053906660e-27_dp, e_charge = 1.602176634e-19_dp
    real(dp), parameter :: eps0 = 8.8541878128e-12_dp, mp = 1.67262192369e-27_dp
    real(dp), parameter :: ma = 183.84_dp*amu, mb = mp, qa = 3.0_dp*e_charge
    real(dp), parameter :: qb = e_charge, nb = 1.0e20_dp, tb = 50.0_dp*e_charge
    real(dp), parameter :: lnl = 15.0_dp, tpar = 200.0_dp*e_charge
    real(dp), parameter :: tperp = 50.0_dp*e_charge, drift = 7.2e3_dp, tail = 1.0e-14_dp
    character(len=32) :: name

    call get_command_argument(1, name)
    select case (trim(name))
    case ('maxwellian')
        call check_maxwellian()
    case ('mass')
        call check_mass()
    case ('slow_rate')
        call check_slow_rate()
    case ('convergence_v')
        call check_convergence_v()
    case ('convergence_dt')
        call check_convergence_dt()
    case default
        error stop 'usage: test_fp_relaxation maxwellian|mass|slow_rate|convergence_v|convergence_dt'
    end select
    print '(a)', 'PASS '//trim(name)

contains

    subroutine setup_b04(s, cells_per_unit)
        type(fp_relax_t), intent(out) :: s
        integer, intent(in) :: cells_per_unit
        integer :: n

        n = fp_relax_cells(cells_per_unit, ma, tb, drift, tpar, tail)
        call fp_relax_setup(s, n, ma, mb, qa, qb, nb, tb, lnl, eps0, tpar, tperp, &
            drift, 1.0_dp, tail)
    end subroutine setup_b04

    subroutine check_maxwellian()
        ! A Maxwellian at T_b with mass m_a is stationary for the test-particle
        ! operator in a Maxwellian background at T_b. The cell-discrete state is
        ! stationary only to O(h^2), so the drift of p/(n T_b) - 1 must shrink
        ! at least twice as fast as h is refined.
        type(fp_relax_t) :: s
        real(dp) :: mom0(fp_relax_nmom), mom1(fp_relax_nmom), times(1), drift(3)
        real(dp) :: mom_out(fp_relax_nmom, 1), vc_val, exact_v4
        integer :: i, k, cpu_list(3)

        cpu_list = [32, 64, 128]
        times = [5.0_dp]
        exact_v4 = 15.0_dp*(tb/ma)**2
        do k = 1, 3
            call setup_b04(s, cpu_list(k))
            do i = 1, s%ncell
                call fp_rel_maxwell(s%vc(i), ma, tb, vc_val)
                s%f(i, 0) = vc_val
            end do
            s%f(:, 1:) = 0.0_dp
            call fp_relax_moments(s, mom0)
            call fp_relax_run(s, times, 0.01_dp, mom_out)
            mom1 = mom_out(:, 1)
            drift(k) = max(abs(mom1(3)/(mom1(1)*tb) - 1.0_dp), &
                abs(mom1(4)/(mom1(1)*tb) - 1.0_dp), &
                abs(mom1(5)/exact_v4 - 1.0_dp))
            print '(a,i4,a,es12.4,a,es12.4)', 'maxwellian cells_per_unit=', cpu_list(k), &
                ' drift=', drift(k), ' density change=', mom1(1)/mom0(1) - 1.0_dp
            print '(a,3es12.4)', '  change t5/t0-1 (p_par, p_perp, v4):', &
                mom1(3)/mom0(3) - 1.0_dp, mom1(4)/mom0(4) - 1.0_dp, mom1(5)/mom0(5) - 1.0_dp
            print '(a,3es12.4)', '  t0 offset (p_par, p_perp, v4):', &
                mom0(3)/(mom0(1)*tb) - 1.0_dp, mom0(4)/(mom0(1)*tb) - 1.0_dp, mom0(5)/exact_v4 - 1.0_dp
            if (abs(mom1(2)) > 1.0e-3_dp*sqrt(exact_v4)) error stop 'Maxwellian drift appeared'
            if (abs(mom1(1)/mom0(1) - 1.0_dp) > 1.0e-10_dp) error stop 'Maxwellian density drifted'
        end do
        if (drift(2) > 1.0e-2_dp) error stop 'Maxwellian drift above 1e-2 at 64 cells'
        if (drift(3) > drift(2)/2.0_dp) error stop 'Maxwellian drift not refining'
    end subroutine check_maxwellian

    subroutine check_mass()
        ! Telescoping face fluxes and the stage weights a0+a1+a2 = 0 conserve
        ! the l=0 cell mass exactly, whatever the initial drifting state.
        type(fp_relax_t) :: s
        real(dp) :: mom0(fp_relax_nmom), mom1(fp_relax_nmom), times(1)
        real(dp) :: mom_out(fp_relax_nmom, 1)

        call setup_b04(s, 64)
        call fp_relax_initial(s)
        call fp_relax_moments(s, mom0)
        times = [5.0_dp]
        call fp_relax_run(s, times, 0.01_dp, mom_out)
        mom1 = mom_out(:, 1)
        print '(a,es24.16,a,es24.16)', 'mass before ', mom0(1), ' after ', mom1(1)
        ! Round-off accumulates over the 500 steps; an O(dt) defect in the stage
        ! weights would show at 1e-4 or more.
        if (abs(mom1(1)/mom0(1) - 1.0_dp) > 1.0e-10_dp) error stop 'FP l=0 mass not conserved'
    end subroutine check_mass

    subroutine check_slow_rate()
        ! A narrow l=1 pulse decays at the v^3-weighted slowing rate nu(v), and
        ! nu(v) = nu_s (1 - 3/5 x^2 + ...) with x = v/v_th,b for small x. The
        ! discrete identity holds to O(h^2), so the mismatch must fall under refinement.
        type(fp_relax_t) :: s
        real(dp), parameter :: v0_frac = 0.2_dp, width_frac = 0.12_dp, dt_tau = 1.0e-4_dp
        real(dp) :: v0, width, m1_0, m1_1, rate_num, nu_bar, nu_s, x2_bar
        real(dp) :: nu_at_v0, x2_v0, nu_point, wsum, mismatch(3)
        real(dp), allocatable :: g(:), nu_cell(:)
        integer :: i, n, k, cpu_list(3)

        cpu_list = [64, 128, 256]
        do k = 1, 3
            call setup_b04(s, cpu_list(k))
            n = s%ncell
            allocate (g(n), nu_cell(n))
            v0 = v0_frac*sqrt(2.0_dp*tb/ma)
            width = width_frac*sqrt(2.0_dp*tb/ma)
            do i = 1, n
                call fp_rel_gauss(s%vc(i), v0, width, g(i))
            end do
            s%f = 0.0_dp
            s%f(:, 1) = g
            m1_0 = sum(s%w3*s%f(:, 1))
            call fp_relax_step(s, dt_tau*s%tau)
            m1_1 = sum(s%w3*s%f(:, 1))
            rate_num = -(m1_1 - m1_0)/(dt_tau*s%tau*m1_0)
            call fp_rel_rate_at(s%vc, s%vt, nb, s%gam, s%mr, nu_cell)
            wsum = sum(s%w3*g)
            nu_bar = sum(s%w3*g*nu_cell)/wsum
            call fp_rel_nu_s(s%gam, s%mr, nb, s%vt, nu_s)
            x2_bar = sum(s%w3*g*(s%vc/s%vt)**2)/wsum
            mismatch(k) = rate_num/nu_bar - 1.0_dp
            print '(a,i4,a,3es16.8,a,es10.2)', 'cells_per_unit=', cpu_list(k), &
                ' nu_s, nu_bar, rate_num =', nu_s, nu_bar, rate_num, ' mismatch=', mismatch(k)
            if (abs(nu_bar/nu_s - (1.0_dp - 0.6_dp*x2_bar)) > 1.0e-4_dp) &
                error stop 'nu_bar misses small-x correction'
            deallocate (g, nu_cell)
        end do
        ! Expansion at the pulse centre, independent of the grid.
        call setup_b04(s, 64)
        x2_v0 = (v0_frac*sqrt(2.0_dp*tb/ma)/s%vt)**2
        call fp_rel_rate_at(v0_frac*sqrt(2.0_dp*tb/ma), s%vt, nb, s%gam, s%mr, nu_point)
        call fp_rel_nu_s(s%gam, s%mr, nb, s%vt, nu_s)
        nu_at_v0 = nu_s*(1.0_dp - 0.6_dp*x2_v0)
        print '(a,es16.8,a,es16.8)', 'nu(v0) =', nu_point, ' expansion =', nu_at_v0
        if (abs(nu_point/nu_at_v0 - 1.0_dp) > 1.0e-6_dp) error stop 'nu(v) misses small-x expansion'
        ! Second order: the mismatch falls by about 4 per doubling, and its
        ! Richardson extrapolation to h -> 0 must vanish.
        if (mismatch(2)/mismatch(3) < 3.0_dp) error stop 'l=1 decay mismatch not second order'
        if (abs(mismatch(3) - (mismatch(2) - mismatch(3))/3.0_dp) > 1.0e-4_dp) &
            error stop 'l=1 decay rate not nu(v) in the limit'
    end subroutine check_slow_rate

    subroutine check_convergence_v()
        ! Spatial order at fixed time step, against the finest grid (128 cells per u_a).
        type(fp_relax_t) :: s
        integer, parameter :: nlev = 4
        integer :: cells(nlev), k
        real(dp) :: ref(fp_relax_nmom), mom(fp_relax_nmom), err(nlev), times(1), order

        cells = [8, 16, 32, 64]
        times = [1.0_dp]
        call run_levels(128, 1.0e-3_dp, times, ref)
        do k = 1, nlev
            call run_levels(cells(k), 1.0e-3_dp, times, mom)
            err(k) = abs(mom(3)/ref(3) - 1.0_dp) + abs(mom(4)/ref(4) - 1.0_dp)
            print '(a,i4,a,es14.6)', 'cells_per_unit=', cells(k), ' error=', err(k)
        end do
        order = log(err(nlev - 1)/err(nlev))/log(2.0_dp)
        print '(a,f8.4)', 'observed spatial order (last pair) =', order
        if (order < 1.5_dp) error stop 'spatial order below 1.5'
    end subroutine check_convergence_v

    subroutine check_convergence_dt()
        ! Time order at fixed grid. TR-BDF2 is second order.
        real(dp) :: ref(fp_relax_nmom), mom(fp_relax_nmom), err(3), times(1), order
        real(dp) :: dts(3)
        integer :: k

        times = [1.0_dp]
        dts = [0.1_dp, 0.05_dp, 0.025_dp]
        call run_levels(32, 0.00625_dp, times, ref)
        do k = 1, 3
            call run_levels(32, dts(k), times, mom)
            err(k) = abs(mom(3)/ref(3) - 1.0_dp) + abs(mom(4)/ref(4) - 1.0_dp)
            print '(a,f8.4,a,es14.6)', 'dt_tau=', dts(k), ' error=', err(k)
        end do
        order = log(err(2)/err(3))/log(2.0_dp)
        print '(a,f8.4)', 'observed time order (last pair) =', order
        if (order < 1.8_dp) error stop 'time order below 1.8'
    end subroutine check_convergence_dt

    subroutine run_levels(cells_per_unit, dt_tau, times, mom)
        ! Returns the moments of the B04 drifting state at times (in tau_s).
        integer, intent(in) :: cells_per_unit
        real(dp), intent(in) :: dt_tau, times(:)
        real(dp), intent(out) :: mom(fp_relax_nmom)
        type(fp_relax_t) :: s
        real(dp) :: out(fp_relax_nmom, size(times))

        call setup_b04(s, cells_per_unit)
        call fp_relax_initial(s)
        call fp_relax_run(s, times, dt_tau, out)
        mom = out(:, size(times))
    end subroutine run_levels

end program test_fp_relaxation
