!> C4: emergent Debye shielding of a weak fixed repulsive point impurity in the
!> classical electron OCP (periodic 3D3V cube, uniform neutralizing background
!> that also compensates the impurity, bare Coulomb via PME, encounter-
!> resolved Strang step; no softening, screened potential or collision
!> operator). Responsive species: electrons only, lambda_De = 1. The impurity
!> charge Q = -kappa 4 pi Lambda e gives the nonlinear radius
!> r_nl = |qQ|/(4 pi eps0 k_B T) = kappa lambda_D (excluded from the fit).
!> Estimator: time- and ensemble-averaged electron density in spherical
!> shells around the impurity, delta n/n = counts/(n V_shell samples) - 1.
!> Reference: linear response -kappa G(r) with the periodic screened-Poisson
!> Green function G (image sum minus zero mode, shell- and bin-averaged);
!> the isolated Yukawa -kappa exp(-r)/r is an overlay with a finite-box model
!> difference. Statistical 95% Student-t intervals over independent
!> realizations are reported separately from the numerical settings.
module ocp_shielding_app
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use fortnum_status, only: fortnum_status_t, FORTNUM_OK
    use kin6d_ocp, only: ocp_model_t, ocp_system, ocp_prepare, ocp_advance, &
        ocp_thermostat, ocp_kinetic_temperature, ocp_mode, ocp_radial_histogram, &
        ocp_debye_shell
    use kin6d_nbody, only: nbody_system_t, nbody_energy
    use kin6d_ocp_leaves, only: debye_structure
    use kin6d_result_io, only: run_t, result_file_t
    use kin6d_animation_export, only: animation_export_t
    use kin6d_series_statistics, only: series_autocorrelation
    use fortplot, only: figure, figure_t, get_global_figure, subplot, scatter, plot, &
        fill_between, title, suptitle, xlabel, ylabel, xlim, ylim, legend, &
        tight_layout, savefig, colorbar
    use fortplot_animation, only: animation_t, animation_func => FuncAnimation, &
                                  save_animation
    use fortplot_system_runtime, only: create_directory_runtime
    implicit none
    private
    public :: main
    real(dp), parameter :: pi = 3.14159265358979323846_dp, two_pi = 2*pi
    real(dp), parameter :: blue(3) = [0.000_dp, 0.447_dp, 0.698_dp]
    real(dp), parameter :: vermillion(3) = [0.835_dp, 0.369_dp, 0.000_dp]
    real(dp), parameter :: green(3) = [0.000_dp, 0.620_dp, 0.451_dp]
    real(dp), parameter :: orange(3) = [0.902_dp, 0.624_dp, 0.000_dp]
    real(dp), parameter :: black(3) = [0.0_dp, 0.0_dp, 0.0_dp]
    integer, parameter :: max_cases = 6, table_points = 41
    ! Configuration (&ocp_shielding).
    character(64) :: label = 'ocp-shielding'
    real(dp) :: debye_volumes(max_cases) = [100.0_dp, 300.0_dp, 1000.0_dp, 100.0_dp, &
                                            0.0_dp, 0.0_dp]
    real(dp) :: impurities(max_cases) = [0.05_dp, 0.05_dp, 0.05_dp, 0.1_dp, &
                                         0.0_dp, 0.0_dp]
    real(dp) :: average_times(max_cases) = [40.0_dp, 30.0_dp, 20.0_dp, 40.0_dp, &
                                            0.0_dp, 0.0_dp]
    integer :: members(max_cases) = [8, 6, 4, 8, 0, 0]
    real(dp) :: box = 4.0_dp, tolerance = 1.0e-4_dp, equilibration_time = 5.0_dp
    real(dp) :: sample_interval = 0.2_dp, snapshot_interval = 10.0_dp
    real(dp) :: fit_min = 0.3_dp, fit_max = 1.8_dp, slice = 0.25_dp
    integer :: bins = 20, directions = 64
    integer(int64) :: seed = 4242_int64
    namelist /ocp_shielding/ label, debye_volumes, impurities, average_times, &
        members, box, tolerance, equilibration_time, sample_interval, &
        snapshot_interval, fit_min, fit_max, slice, bins, directions, seed
    type(run_t) :: run
    type(result_file_t) :: result
    type(figure_t), pointer :: canvas
    type(animation_t) :: animation
    integer :: cases, case_index, samples, per, frame, selected
    real(dp) :: rmax
    real(dp), allocatable :: rb(:), periodic(:), yukawa(:), table(:, :), kappas(:)
    ! Per case: profiles (bin, member), statistics and gates.
    real(dp), allocatable :: profile(:, :, :), pstat(:, :, :), fitted(:, :)
    real(dp), allocatable :: gates(:, :), cost(:, :, :), autocorr(:, :)
    real(dp), allocatable :: sel_x(:, :), sel_y(:, :), running(:, :), dts(:)
    integer, allocatable :: picked(:), electrons(:)
    ! Animation export (--export DIR): case 1 only; all electron positions of
    ! member 1 and the cumulative shell counts of every member per sample.
    character(256) :: export_dir = ''
    real(dp), allocatable :: trajectory(:, :, :), member_counts(:, :, :)
    character(160) :: heading
contains
    subroutine main()
        integer :: unit, ios, c
        logical :: ok
        character(256) :: msg
        call run%start('ocp_shielding_demo')
        call run%open_config(unit)
        if (run%has_config()) read (unit, nml=ocp_shielding, iostat=ios, iomsg=msg)
        call run%close_config(unit, 'ocp_shielding', ios, msg)
        call run%open_resolved(unit)
        write (unit, nml=ocp_shielding)
        call run%close_resolved(unit)
        c = 1
        do while (c <= run%argument_count())
            call run%argument(c, msg)
            if (trim(msg) /= '--export' .or. c == run%argument_count()) &
                error stop 'usage: ocp_shielding_demo [--config FILE] [--export DIR]'
            call run%argument(c + 1, export_dir); c = c + 2
        end do
        cases = count(debye_volumes > 0)
        if (cases < 1 .or. any(members(:cases) < 2)) error stop 'need cases with >= 2 members'
        if (fit_max > box/2 .or. fit_min <= 0) error stop 'fit window outside (0, L/2]'
        call create_directory_runtime('output', ok)
        rmax = box/2
        allocate (rb(bins), periodic(bins), yukawa(bins), kappas(table_points))
        allocate (table(bins, table_points), electrons(cases), dts(cases))
        allocate (profile(bins, maxval(members), cases), pstat(bins, 3, cases))
        allocate (fitted(maxval(members) + 1, cases), gates(12, cases))
        allocate (cost(maxval(members), 5, cases), autocorr(maxval(members), cases))
        call references()
        do c = 1, cases
            case_index = c
            call run_case(c)
            print '(a,f7.1,a,f5.3,a,f7.4,a,f7.4,a,f7.4)', 'Lambda = ', debye_volumes(c), &
                ' kappa = ', impurities(c), ': lambda_fit = ', gates(1, c), ' +- ', &
                gates(2, c), ' (mean-profile fit ', gates(3, c), ')'
        end do
        if (len_trim(export_dir) > 0) call export_animation()
        call write_outputs()
        call draw_figures()
    end subroutine

    !> Bin-volume average sum_q w_q r_q^2 f(r_q)/sum_q w_q r_q^2 (4-point
    !> Gauss-Legendre) of the shell-averaged periodic Green function and of
    !> the isolated Yukawa, per unit impurity coupling kappa; table over the
    !> screening wave number for the fits.
    subroutine references()
        real(dp), parameter :: xg(4) = [-0.861136311594053_dp, -0.339981043584856_dp, &
            0.339981043584856_dp, 0.861136311594053_dp]
        real(dp), parameter :: wg(4) = [0.347854845137454_dp, 0.652145154862546_dp, &
            0.652145154862546_dp, 0.347854845137454_dp]
        real(dp) :: dr
        integer :: b, i
        dr = rmax/bins
        do b = 1, bins
            rb(b) = (b - 0.5_dp)*dr
        end do
        do i = 1, table_points
            kappas(i) = 0.5_dp + 1.5_dp*(i - 1)/(table_points - 1)
        end do
        !$omp parallel do schedule(dynamic, 1) private(i, b)
        do i = 1, table_points
            do b = 1, bins
                table(b, i) = average(b, kappas(i), .true.)
            end do
        end do
        !$omp end parallel do
        do b = 1, bins
            periodic(b) = -average(b, 1.0_dp, .true.)
            yukawa(b) = -average(b, 1.0_dp, .false.)
        end do
    contains
        real(dp) function average(b, kappa, lattice) result(v)
            integer, intent(in) :: b
            real(dp), intent(in) :: kappa
            logical, intent(in) :: lattice
            real(dp) :: r, f, wsum
            integer :: q
            v = 0; wsum = 0
            do q = 1, 4
                r = rb(b) + xg(q)*dr/2
                if (lattice) then
                    f = ocp_debye_shell(r, kappa, box, directions)
                else
                    f = exp(-kappa*r)/r
                end if
                v = v + wg(q)*r*r*f; wsum = wsum + wg(q)*r*r
            end do
            v = v/wsum
        end function
    end subroutine

    !> Periodic reference at screening wave number kappa (cubic Lagrange
    !> interpolation of the table), per unit impurity coupling, negative.
    real(dp) function reference(b, kappa) result(v)
        integer, intent(in) :: b
        real(dp), intent(in) :: kappa
        real(dp) :: p, x, w(4)
        integer :: base, q
        p = (kappa - kappas(1))/(kappas(2) - kappas(1))
        base = min(table_points - 4, max(0, int(p) - 1))
        x = p - base
        w(1) = -(x - 1)*(x - 2)*(x - 3)/6; w(2) = x*(x - 2)*(x - 3)/2
        w(3) = -x*(x - 1)*(x - 3)/2; w(4) = x*(x - 1)*(x - 2)/6
        v = 0
        do q = 1, 4
            v = v + w(q)*table(b, base + q)
        end do
        v = -v
    end function

    !> Screening length 1/kappa_fit minimizing sum over the fit window of
    !> (dn/n - kappa_Q reference(kappa))^2, golden section on [0.5, 2].
    real(dp) function fit_length(dn, coupling) result(length)
        real(dp), intent(in) :: dn(:), coupling
        real(dp) :: lo, hi, a, b, fa, fb
        integer :: it
        lo = kappas(1); hi = kappas(table_points)
        do it = 1, 100
            a = hi - 0.618033988749895_dp*(hi - lo); b = lo + 0.618033988749895_dp*(hi - lo)
            fa = misfit(a); fb = misfit(b)
            if (fa < fb) then
                hi = b
            else
                lo = a
            end if
        end do
        length = 2/(lo + hi)
    contains
        real(dp) function misfit(kappa) result(s)
            real(dp), intent(in) :: kappa
            integer :: q
            s = 0
            do q = 1, bins
                if (rb(q) < fit_min .or. rb(q) > fit_max) cycle
                s = s + (dn(q) - coupling*reference(q, kappa))**2
            end do
        end function
    end function

    !> Two-sided 95% Student-t quantile and the mean/half-width of x.
    pure subroutine interval(x, mean, half)
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: mean, half
        real(dp), parameter :: table95(10) = [12.706_dp, 4.303_dp, 3.182_dp, 2.776_dp, &
            2.571_dp, 2.447_dp, 2.365_dp, 2.306_dp, 2.262_dp, 2.228_dp]
        real(dp) :: t
        integer :: n
        n = size(x); mean = sum(x)/n
        t = 1.96_dp
        if (n - 1 <= 10) t = table95(max(1, n - 1))
        if (n - 1 > 10 .and. n - 1 <= 30) t = 2.228_dp - (n - 11)*(2.228_dp - 2.042_dp)/20
        half = t*sqrt(sum((x - mean)**2)/(n - 1)/n)
    end subroutine

    subroutine run_case(c)
        integer, intent(in) :: c
        type(ocp_model_t) :: model
        type(nbody_system_t) :: probe
        real(dp) :: dt, mean, half, sk(maxval(members)), z, chi
        integer :: r, nm, b, nwin
        logical :: ok
        nm = members(c)
        model%debye_volume = debye_volumes(c); model%box = box
        model%tolerance = tolerance; model%impurity = impurities(c)
        dt = 0.01_dp*(100/model%debye_volume)**(1.0_dp/3)
        per = max(1, nint(sample_interval/dt)); dt = sample_interval/per
        model%dt = dt; dts(c) = dt
        call ocp_system(model, probe, ok)
        if (.not. ok) error stop 'OCP impurity system setup failed'
        electrons(c) = model%electrons
        samples = nint(average_times(c)/sample_interval)
        if (c == 1) allocate (sel_x(1, 1), sel_y(1, 1), running(bins, samples), picked(1))
        if (c == 1 .and. len_trim(export_dir) > 0) allocate ( &
            trajectory(3, model%electrons, samples), member_counts(bins, samples, nm))
        print '(a,i0,a,f6.4,a,i0,a,es9.2)', 'case: N = ', model%electrons, ', dt = ', dt, &
            ', PME K = ', model%pme%mesh, ', estimated rms force error/F0 = ', &
            sqrt(model%pme%real_estimate**2 + model%pme%reciprocal_estimate**2)/ &
            model%pme%force_scale
        !$omp parallel do schedule(dynamic, 1) private(r)
        do r = 1, nm
            call member(model, c, r, sk(r))
        end do
        !$omp end parallel do
        do b = 1, bins
            call interval(profile(b, :nm, c), mean, half)
            pstat(b, :, c) = [mean, mean - half, mean + half]
        end do
        do r = 1, nm
            fitted(r, c) = fit_length(profile(:, r, c), impurities(c))
        end do
        fitted(nm + 1, c) = fit_length(pstat(:, 1, c), impurities(c))
        call interval(fitted(:nm, c), gates(1, c), gates(2, c))
        gates(3, c) = fitted(nm + 1, c)
        call interval(sk(:nm), gates(4, c), gates(5, c))
        call debye_structure((two_pi/box)**2, 1.0_dp, gates(6, c))
        ! Geometry gate: periodic vs isolated reference over the fit window,
        ! relative to the reference magnitude at the window start.
        gates(7, c) = 0; chi = 0; nwin = 0
        do b = 1, bins
            if (rb(b) < fit_min .or. rb(b) > fit_max) cycle
            gates(7, c) = max(gates(7, c), abs(periodic(b) - yukawa(b)))
            ! Mean z-score of the measured mean against the periodic reference.
            half = (pstat(b, 3, c) - pstat(b, 1, c))/1.96_dp
            z = (pstat(b, 1, c) - impurities(c)*periodic(b))/max(half, tiny(1.0_dp))
            chi = chi + z*z; nwin = nwin + 1
        end do
        gates(7, c) = gates(7, c)/abs(periodic(first_bin()))
        gates(8, c) = sum(autocorr(:nm, c))/nm
        gates(9, c) = nm*samples*sample_interval/max(2*gates(8, c), sample_interval)
        gates(10, c) = maxval(abs(cost(:nm, 5, c)))
        gates(11, c) = sqrt(model%pme%real_estimate**2 + model%pme%reciprocal_estimate**2)/ &
                       model%pme%force_scale
        gates(12, c) = chi/max(1, nwin)
    end subroutine

    integer function first_bin() result(b)
        do b = 1, bins
            if (rb(b) >= fit_min) return
        end do
        b = bins
    end function

    !> One member: prepared ensemble with the bare-impurity Boltzmann factor,
    !> thermostat equilibration (switched off), then Hamiltonian sampling.
    subroutine member(template, c, r, skmin)
        type(ocp_model_t), intent(in) :: template
        integer, intent(in) :: c, r
        real(dp), intent(out) :: skmin
        type(ocp_model_t) :: m
        type(nbody_system_t) :: sys
        type(fortnum_status_t) :: st
        real(dp), allocatable :: s(:), counts(:), series(:), acc(:)
        real(dp) :: e0, e1, shell, density, kv(3), mean, err, tau, ess, rho, sampling, &
            wall0
        integer(int64) :: c0, c1, c2, rate
        integer :: rejected, step, j, b, q
        logical :: okm
        m = template
        call ocp_system(m, sys, okm)
        if (.not. okm) error stop 'member system failed'
        allocate (s(6*m%particles), counts(bins), series(samples), acc(bins))
        call system_clock(c0, rate)
        call ocp_prepare(m, seed + 1000_int64*c + r, s, rejected, st)
        do step = 1, nint(equilibration_time/m%dt)
            call ocp_advance(sys, s, m%dt, st)
            if (st%code /= FORTNUM_OK) error stop 'equilibration failed: '//trim(st%msg)
            call ocp_thermostat(m, s)
        end do
        call system_clock(c1); wall0 = real(c1 - c0, dp)/rate
        call nbody_energy(sys, s, e0, st)
        counts = 0; skmin = 0; sampling = 0
        if (c == 1 .and. r == 1) call select_slice(s)
        do j = 1, samples
            do step = 1, per
                call ocp_advance(sys, s, m%dt, st)
                if (st%code /= FORTNUM_OK) error stop 'OCP step failed: '//trim(st%msg)
            end do
            call system_clock(c1)
            acc = 0
            call ocp_radial_histogram(m, s, [0.0_dp, 0.0_dp, 0.0_dp], rmax, acc)
            counts = counts + acc
            series(j) = sum(acc(:max(1, nint(bins/rmax))))
            do q = 1, 3
                kv = 0; kv(q) = two_pi/box
                rho = abs(ocp_mode(m, s, kv))
                skmin = skmin + rho*rho/m%electrons/3
            end do
            if (c == 1 .and. r == 1) then
                call store_slice(s, j)
                running(:, j) = counts
            end if
            if (c == 1 .and. allocated(member_counts)) then
                member_counts(:, j, r) = counts
                if (r == 1) trajectory(:, :, j) = reshape(s(:3*m%electrons), &
                    [3, m%electrons]) - box*anint(reshape(s(:3*m%electrons), &
                    [3, m%electrons])/box)
            end if
            if (mod(j, max(1, nint(snapshot_interval/sample_interval))) == 0) call dump(s, j)
            call system_clock(c2); sampling = sampling + real(c2 - c1, dp)/rate
        end do
        call nbody_energy(sys, s, e1, st)
        skmin = skmin/samples
        density = m%electrons/box**3
        do b = 1, bins
            shell = 4*pi/3*((b*rmax/bins)**3 - ((b - 1)*rmax/bins)**3)
            profile(b, r, c) = counts(b)/(density*shell*samples) - 1
            if (c == 1 .and. r == 1) running(b, :) = running(b, :)/(density*shell* &
                [(real(j, dp), j = 1, samples)]) - 1
        end do
        ! Integrated autocorrelation time of the inner-shell count series
        ! (kin6d_series_statistics, Sokal window), in time units.
        call series_autocorrelation(series, mean, err, tau, ess)
        autocorr(r, c) = tau*sample_interval
        call system_clock(c2)
        cost(r, :, c) = [wall0, real(c2 - c0, dp)/rate - wall0, sampling, &
                         real(rejected, dp), (e1 - e0)/m%electrons]
    contains
        !> Explicit selection: electrons with |z| < slice (minimum image) at
        !> the start of sampling, followed by identity.
        subroutine select_slice(x)
            real(dp), intent(in) :: x(:)
            integer :: p, n
            n = 0
            do p = 1, m%electrons
                if (abs(x(3*p) - box*anint(x(3*p)/box)) < slice) n = n + 1
            end do
            deallocate (sel_x, sel_y, picked)
            allocate (sel_x(n, samples), sel_y(n, samples), picked(n))
            selected = n; n = 0
            do p = 1, m%electrons
                if (abs(x(3*p) - box*anint(x(3*p)/box)) < slice) then
                    n = n + 1; picked(n) = p
                end if
            end do
        end subroutine

        subroutine store_slice(x, j)
            real(dp), intent(in) :: x(:)
            integer, intent(in) :: j
            integer :: q, p
            do q = 1, selected
                p = picked(q)
                sel_x(q, j) = x(3*p - 2) - box*anint(x(3*p - 2)/box)
                sel_y(q, j) = x(3*p - 1) - box*anint(x(3*p - 1)/box)
            end do
        end subroutine

        !> Snapshot (stream): N, t, Lambda, L, impurity charge, x(3N), v(3N).
        subroutine dump(x, j)
            real(dp), intent(in) :: x(:)
            integer, intent(in) :: j
            integer :: u
            character(160) :: path
            write (path, '(a,a,a,i0,a,i0,a,i0,a)') 'output/', trim(label), '-snapshot-c', c, &
                '-r', r, '-s', j, '.bin'
            open (newunit=u, file=trim(path), access='stream', form='unformatted', &
                  status='replace', action='write')
            write (u) m%electrons, j*sample_interval, m%debye_volume, box, &
                m%impurity_charge, x(:3*m%electrons), &
                x(3*m%particles + 1:3*m%particles + 3*m%electrons)
            close (u)
        end subroutine
    end subroutine

    !> Animation data for the lecture scene (slug debye-shielding-particles),
    !> case 1: positions of all electrons of member 1 (minimum image around the
    !> fixed impurity at the origin, 16-bit over [-L/2, L/2]), the running
    !> ensemble mean of the time-averaged shell deficit delta n/n up to each
    !> sample, the final per-bin 95% Student-t interval over members, and the
    !> bin-averaged references per unit kappa (periodic Debye, isolated
    !> Yukawa, unscreened 1/r). Statistical status only; no deterministic
    !> certificate is claimed.
    subroutine export_animation()
        type(animation_export_t) :: export
        real(dp), allocatable :: mean_profile(:, :), bare(:)
        real(dp) :: density, shell
        integer :: b, j
        allocate (mean_profile(bins, samples), bare(bins))
        density = electrons(1)/box**3
        do b = 1, bins
            shell = 4*pi/3*((b*rmax/bins)**3 - ((b - 1)*rmax/bins)**3)
            do j = 1, samples
                mean_profile(b, j) = sum(member_counts(b, j, :members(1)))/ &
                    (members(1)*density*shell*j) - 1
            end do
            bare(b) = -bin_mean_inverse_r(b)
        end do
        call export%create(run, export_dir, 'debye-shielding-particles')
        call export%note('model', 'classical electron OCP, periodic cube, uniform '// &
            'neutralizing background, bare Coulomb (PME), fixed repulsive point '// &
            'impurity Q = -kappa 4 pi Lambda e compensated in the background (C4)')
        call export%note('units', 'k_B T = m = 1, lambda_D = 1, omega_p = 1')
        call export%note('status', 'candidate dynamics; profile intervals are '// &
            'statistical (95% Student-t over members), not certificates')
        call export%note('references', 'delta n/n reference = kappa * periodic '// &
            '(linear response with the periodic screened-Poisson Green function); '// &
            'yukawa: isolated kappa exp(-r)/r overlay; bare: unscreened linear '// &
            'Boltzmann response kappa/r; all bin-averaged, sign included')
        call export%parameter('Lambda', debye_volumes(1), '1', 'n lambda_D^3')
        call export%parameter('kappa', impurities(1), '1', &
            'impurity coupling, nonlinear radius r_nl = kappa lambda_D')
        call export%parameter('box', box, 'lambda_D', 'periodic cube edge L')
        call export%parameter('electrons', real(electrons(1), dp), '1', 'N')
        call export%parameter('members', real(members(1), dp), '1', &
            'independent realizations in the profile')
        call export%parameter('dt', dts(1), '1/omega_p', 'time step')
        call export%parameter('sample_interval', sample_interval, '1/omega_p', &
            'time between frames')
        call export%parameter('equilibration_time', equilibration_time, '1/omega_p', &
            'thermostat equilibration before t = 0 (switched off afterwards)')
        call export%real32('t', [(j*sample_interval, j = 1, samples)], [samples], &
            '1/omega_p', 'time since the end of equilibration')
        call export%quantized('positions', reshape(trajectory, [size(trajectory)]), &
            [3, electrons(1), samples], -box/2, box/2, 16, 'lambda_D', &
            'electron positions (x, y, z) of member 1, impurity at the origin')
        call export%real32('r', rb, [bins], 'lambda_D', 'shell-bin centres')
        call export%real32('dn_running', reshape(mean_profile, [size(mean_profile)]), &
            [bins, samples], '1', 'delta n/n, ensemble mean of the time average '// &
            'up to each frame')
        call export%real32('dn_interval', reshape(pstat(:, :, 1), [3*bins]), [bins, 3], &
            '1', 'final delta n/n: mean, lower, upper 95% Student-t over members')
        call export%real32('periodic', periodic, [bins], '1', &
            'periodic Debye reference per unit kappa (delta n/(n kappa))')
        call export%real32('yukawa', yukawa, [bins], '1', &
            'isolated -exp(-r)/r per unit kappa')
        call export%real32('bare', bare, [bins], '1', 'unscreened -1/r per unit kappa')
        call export%close()
        print '(2a)', 'animation data written to ', &
            trim(export_dir)//'/debye-shielding-particles'
    contains
        !> Bin-volume average of 1/r (exact: (3/2)(b^2 - a^2)/(b^3 - a^3)).
        real(dp) function bin_mean_inverse_r(b) result(v)
            integer, intent(in) :: b
            real(dp) :: lo, hi
            lo = (b - 1)*rmax/bins; hi = b*rmax/bins
            v = 1.5_dp*(hi**2 - lo**2)/(hi**3 - lo**3)
        end function
    end subroutine

    subroutine write_outputs()
        integer :: u, c, b, r
        open (newunit=u, file='output/'//trim(label)//'-profile.csv', status='replace')
        write (u, '(a)') '# C4 OCP Debye shielding; electron-only lambda_De = 1; '// &
            '95% Student-t statistical intervals over realizations'
        write (u, '(a)') 'case,lambda,kappa_nl,r_over_lambda_d,dn_mean,dn_lo,dn_hi,'// &
            'periodic_debye,isolated_yukawa'
        do c = 1, cases
            do b = 1, bins
                write (u, '(i0,",",*(g0.8,:,","))') c, debye_volumes(c), impurities(c), &
                    rb(b), pstat(b, :, c), impurities(c)*periodic(b), &
                    impurities(c)*yukawa(b)
            end do
        end do
        close (u)
        open (newunit=u, file='output/'//trim(label)//'-gates.csv', status='replace')
        write (u, '(a)') 'case,lambda,kappa_nl,n,members,dt,lambda_fit_mean,'// &
            'lambda_fit_half95,lambda_fit_mean_profile,s_kmin,s_kmin_half95,'// &
            's_kmin_debye_hueckel,geometry_periodic_vs_yukawa,autocorrelation_time,'// &
            'effective_samples,max_energy_defect_per_particle,pme_rms_force_estimate,'// &
            'mean_z2_vs_periodic'
        do c = 1, cases
            write (u, '(i0,",",g0.6,",",g0.6,",",i0,",",i0,",",*(g0.8,:,","))') c, &
                debye_volumes(c), impurities(c), electrons(c), members(c), dts(c), &
                gates(:, c)
        end do
        close (u)
        open (newunit=u, file='output/'//trim(label)//'-members.csv', status='replace')
        write (u, '(a)') 'case,member,lambda_fit,equilibration_s,sampling_run_s,'// &
            'observable_s,rejected_draws,energy_defect_per_particle,autocorrelation_time'
        do c = 1, cases
            do r = 1, members(c)
                write (u, '(i0,",",i0,",",*(g0.8,:,","))') c, r, fitted(r, c), &
                    cost(r, :, c), autocorr(r, c)
            end do
        end do
        close (u)
        call result%create(run)
        call result%attribute('model', 'classical electron OCP with fixed repulsive '// &
            'point impurity; periodic 3D3V; bare Coulomb (PME); no softening')
        call result%attribute('responsive_species', 'electrons only; lambda_De = 1')
        call result%attribute('preparation', 'uniform x Boltzmann factor exp(-kappa/r) '// &
            'of the bare impurity, Maxwellian T=1, isokinetic thermostat during '// &
            'equilibration only')
        call result%write('r', rb, 'lambda_D', 'bin centre', dims='bin')
        call result%write('periodic_debye', periodic, '1', &
            'periodic screened-Poisson reference per unit kappa', dims='bin')
        call result%write('isolated_yukawa', yukawa, '1', &
            'isolated Yukawa overlay per unit kappa', dims='bin')
        call result%write('profile_statistics', pstat(:, :, :cases), '1', &
            'delta n/n mean/lo/hi', dims='bin,statistic3,case', &
            status='statistical 95% Student-t over realizations')
        call result%write('gates', gates(:, :cases), '1', &
            'see output csv header for the gate columns', dims='gate,case')
        call result%write('member_fits', fitted(:, :cases), 'lambda_D', &
            'per-member screening length; last row: mean-profile fit', &
            dims='member_plus_mean,case')
        call result%close()
    end subroutine

    subroutine draw_figures()
        real(dp) :: lam(cases), lo(cases), hi(cases), view
        integer :: c, frames
        call figure(figsize=[13.2_dp, 8.8_dp])
        do c = 1, min(cases, 4)
            call subplot(2, 2, c)
            call fill_between(rb, -pstat(:, 2, c)/impurities(c), -pstat(:, 3, c)/ &
                              impurities(c), color=blue, alpha=0.25_dp)
            call plot(rb, -pstat(:, 1, c)/impurities(c), color=blue, linewidth=1.8_dp, &
                      label='microscopic (95% stat.)')
            call plot(rb, -yukawa, color=green, linestyle=':', linewidth=1.2_dp, &
                      label='isolated Yukawa exp(-r)/r')
            call plot(rb, -periodic, color=black, linestyle='--', linewidth=1.2_dp, &
                      label='periodic Debye G(r)')
            call plot([fit_min, fit_min], [-1.0_dp, 6.0_dp], color=[0.6_dp, 0.6_dp, &
                      0.6_dp], linestyle=':', linewidth=1.0_dp, label='fit window start')
            write (heading, '(a,f6.0,a,f5.2,a,f6.3,a,f5.3)') '$\Lambda$ = ', debye_volumes(c), &
                ', $\kappa$ = ', impurities(c), ': $\lambda_{fit}$ = ', gates(1, c), ' ± ', gates(2, c)
            call title(trim(heading))
            call xlabel('$r/\lambda_D$'); call ylabel('$-\delta n/(n\kappa)$')
            call xlim(0.0_dp, rmax); call ylim(-0.5_dp, 5.0_dp)
            call legend(loc='upper right')
        end do
        call suptitle('C4 emergent Debye shielding (electron-only $\lambda_{De}$); bands '// &
                      'statistical, not certified')
        call tight_layout()
        call savefig('output/'//trim(label)//'-profiles.png')
        lam = debye_volumes(:cases)
        lo = gates(1, :cases) - gates(2, :cases); hi = gates(1, :cases) + gates(2, :cases)
        call figure(figsize=[12.8_dp, 6.4_dp])
        call subplot(1, 2, 1)
        call scatter(lam, gates(1, :cases), color=blue, s=36.0_dp, &
                     label='member mean $\lambda_{fit}$')
        do c = 1, cases
            call plot([lam(c), lam(c)], [lo(c), hi(c)], color=blue, linewidth=1.8_dp)
        end do
        call scatter(lam, gates(3, :cases), color=vermillion, s=36.0_dp, marker='s', &
                     label='fit of the mean profile')
        call plot([0.5_dp*minval(lam), 2*maxval(lam)], [1.0_dp, 1.0_dp], color=black, &
                  linestyle='--', linewidth=1.2_dp, label='$\lambda_{De}$ = 1')
        call title('Fitted screening length (95% stat. bars)')
        call xlabel('$\Lambda = n\lambda_D^3$'); call ylabel('$\lambda_{fit}$ / $\lambda_{De}$')
        call xlim(0.0_dp, 1.1_dp*maxval(lam)); call ylim(0.5_dp, 1.5_dp)
        call legend(loc='upper right')
        call subplot(1, 2, 2)
        call scatter(lam, gates(4, :cases), color=blue, s=36.0_dp, &
                     label='measured $S(k_{min})$')
        do c = 1, cases
            call plot([lam(c), lam(c)], [gates(4, c) - gates(5, c), gates(4, c) + &
                      gates(5, c)], color=blue, linewidth=1.8_dp)
        end do
        call plot([0.5_dp*minval(lam), 2*maxval(lam)], [gates(6, 1), gates(6, 1)], &
                  color=black, linestyle='--', linewidth=1.2_dp, label='Debye-Hückel')
        call title('Fluctuation gate: $S(2\pi/L)$ with impurity present')
        call xlabel('$\Lambda = n\lambda_D^3$'); call ylabel('$S(k_{min})$')
        call xlim(0.0_dp, 1.1_dp*maxval(lam)); call ylim(0.0_dp, 1.2_dp)
        call legend(loc='lower right')
        call tight_layout()
        call savefig('output/'//trim(label)//'-scan.png')
        view = box/2
        call figure(figsize=[12.8_dp, 6.4_dp])
        canvas => get_global_figure()
        samples = size(running, 2)
        frame = samples
        call draw_frame()
        call savefig('output/'//trim(label)//'-particles.png')
        frames = samples + 24
        animation = animation_func(update_frame, frames=frames, interval=42, fig=canvas)
        call save_animation(animation, 'output/'//trim(label)//'.mp4', fps=24, status=c)
        if (c /= 0) print '(a)', 'animation encoding failed'
    contains
        subroutine update_frame(index)
            integer, intent(in) :: index
            frame = min(samples, index)
            call canvas%clear()
            call draw_frame()
            call canvas%set_rendered(.false.)
        end subroutine

        subroutine draw_frame()
            call subplot(1, 2, 1)
            call scatter(sel_x(:, frame), sel_y(:, frame), color=blue, s=6.0_dp, &
                         alpha=0.5_dp, label='electrons, |z| < 0.25 $\lambda_D$')
            call scatter([0.0_dp], [0.0_dp], color=vermillion, s=80.0_dp, marker='x', &
                         label='fixed impurity')
            write (heading, '(a,f6.2)') 'Slice through the impurity, $\omega_p t$ = ', &
                frame*sample_interval
            call title(trim(heading))
            call xlabel('$x/\lambda_D$'); call ylabel('$y/\lambda_D$')
            call xlim(-view, view); call ylim(-view, view)
            call legend(loc='upper right')
            call subplot(1, 2, 2)
            call plot(rb, -running(:, frame)/impurities(1), color=blue, linewidth=1.8_dp, &
                      label='running average, one realization')
            call plot(rb, -periodic, color=black, linestyle='--', linewidth=1.2_dp, &
                      label='periodic Debye G(r)')
            write (heading, '(a,f6.1,a,f6.2)') 'Accumulated response, $\Lambda$ = ', &
                debye_volumes(1), ', $\omega_p t$ = ', frame*sample_interval
            call title(trim(heading))
            call xlabel('$r/\lambda_D$'); call ylabel('$-\delta n/(n\kappa)$')
            call xlim(0.0_dp, rmax); call ylim(-3.0_dp, 6.0_dp)
            call legend(loc='upper right')
            call tight_layout()
        end subroutine
    end subroutine
end module ocp_shielding_app

program ocp_shielding_demo
    use ocp_shielding_app, only: main
    implicit none
    call main()
end program ocp_shielding_demo
