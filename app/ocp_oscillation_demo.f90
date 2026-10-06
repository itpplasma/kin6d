module ocp_oscillation_app
    ! C3: microscopic plasma oscillation of the classical electron OCP in a
    ! periodic 3D3V cube (uniform neutralizing background, bare Coulomb via
    ! PME, encounter-resolved Strang step, no softening/collision operator).
    ! Each ensemble member: uniform positions + Maxwellian velocities, a
    ! declared isokinetic-thermostat equilibration (a different dynamics,
    ! switched off before measurement), then from the same equilibrated state
    ! (common random numbers) an unperturbed run, a run with the spatial-mode
    ! displacement x_a += (eps/k) sin(k x_a) on all three axes, and the step
    ! tangent with respect to eps at eps = 0. Observable: the normalized
    ! density-mode response R(t) = (2/N) Re rho_k(t)/(-eps), averaged over the
    ! three axes; ensemble mean with a 95% Student-t interval (statistical).
    ! Oracles: linear Vlasov initial-value response (Volterra), Landau root
    ! of the Maxwellian dielectric, Bohm-Gross and the cold rigid slab
    ! (kin6d_slab_oscillation). Snapshots and g(r)/S(k) feed C5/H1.
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use fortnum_status, only: fortnum_status_t, FORTNUM_OK
    use kin6d_ocp, only: ocp_model_t, ocp_system, ocp_prepare, ocp_advance, &
        ocp_advance_jvp, ocp_thermostat, ocp_kinetic_temperature, ocp_mode, &
        ocp_displace, ocp_displace_tangent
    use kin6d_nbody, only: nbody_system_t, nbody_energy
    use kin6d_vlasov_response, only: landau_root, vlasov_density_response, &
        fit_damped_cosine
    use kin6d_pair_correlation, only: pair_counts, pair_distribution, structure_factor
    use kin6d_slab_oscillation, only: slab_oscillation_state_t, &
        initialize_slab_oscillation, solve_slab_oscillation, sample_slab_oscillation
    use kin6d_ocp_leaves, only: bohm_gross, debye_structure
    use kin6d_result_io, only: run_t, result_file_t
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
    real(dp), parameter :: purple(3) = [0.800_dp, 0.475_dp, 0.655_dp]
    real(dp), parameter :: black(3) = [0.0_dp, 0.0_dp, 0.0_dp]
    real(dp), parameter :: grey(3) = [0.55_dp, 0.55_dp, 0.55_dp]
    integer, parameter :: profile_bins = 32, radial_bins = 40
    ! Configuration (&ocp_oscillation).
    character(64) :: label = 'ocp-oscillation-landau'
    real(dp) :: debye_volume = 20.0_dp, box = 12.566370614359172_dp
    real(dp) :: tolerance = 1.0e-4_dp, amplitude = 0.01_dp
    real(dp) :: equilibration_time = 5.0_dp, measure_time = 25.0_dp
    real(dp) :: sample_interval = 0.1_dp, snapshot_interval = 5.0_dp
    real(dp) :: fit_start = 5.0_dp, fit_end = 20.0_dp, slice = 0.25_dp
    integer :: realizations = 8, structure_shells = 0
    integer(int64) :: seed = 20261005_int64
    logical :: tangent = .true.
    namelist /ocp_oscillation/ label, debye_volume, box, tolerance, amplitude, &
        equilibration_time, measure_time, sample_interval, snapshot_interval, &
        fit_start, fit_end, slice, realizations, structure_shells, seed, tangent
    type(run_t) :: run
    type(result_file_t) :: result
    type(ocp_model_t) :: model
    type(nbody_system_t) :: probe
    type(slab_oscillation_state_t) :: slab
    type(figure_t), pointer :: canvas
    type(animation_t) :: animation
    real(dp), allocatable :: time(:), pair(:, :), linear(:, :), temperature(:, :)
    real(dp), allocatable :: energy(:, :), tnorm(:, :), volterra(:), landau(:)
    real(dp), allocatable :: bgross(:), cold(:), stat(:, :), fits(:, :), cost(:, :)
    real(dp), allocatable :: sk(:, :), gcount(:, :), sel_x(:, :), sel_y(:, :)
    real(dp), allocatable :: sel_v(:, :), profile(:, :), xprofile(:), rk(:), sdh(:)
    real(dp), allocatable :: rr(:), gr(:, :), gstat(:, :), skstat(:, :), dh_g(:)
    integer, allocatable :: vectors(:), snapshots(:), picked(:)
    complex(dp) :: root
    real(dp) :: k, dt, omega_bg, vfit(3), volterra_fit(3), t0, rmax, lyapunov(2)
    real(dp) :: view_profile, vmax_colour
    integer :: samples, per, r, j, unit, ios, nshell, selected, frame, encode_status
    integer :: nsnap
    logical :: ok, have_root
    character(256) :: msg
    character(160) :: heading

contains
    subroutine main()
        call run%start('ocp_oscillation_demo')
        call run%open_config(unit)
        if (run%has_config()) read (unit, nml=ocp_oscillation, iostat=ios, iomsg=msg)
        call run%close_config(unit, 'ocp_oscillation', ios, msg)
        call run%open_resolved(unit)
        write (unit, nml=ocp_oscillation)
        call run%close_resolved(unit)
        if (realizations < 2) error stop 'need at least two realizations for intervals'
        if (.not. (amplitude > 0 .and. measure_time > 0 .and. sample_interval > 0)) &
            error stop 'invalid amplitude/time configuration'
        call create_directory_runtime('output', ok)
        ! Outer step: the automatic value, rounded so samples fall on steps.
        model%debye_volume = debye_volume; model%box = box; model%tolerance = tolerance
        dt = 0.01_dp*(100/debye_volume)**(1.0_dp/3)
        per = max(1, nint(sample_interval/dt)); dt = sample_interval/per
        model%dt = dt
        call ocp_system(model, probe, ok)
        if (.not. ok) error stop 'OCP system setup failed'
        samples = nint(measure_time/sample_interval) + 1
        nsnap = nint(measure_time/snapshot_interval) + 1
        k = two_pi/box
        nshell = structure_shells
        if (nshell <= 0) nshell = max(1, (nint(3*box/two_pi))**2)
        rmax = min(2.0_dp, box/2)
        allocate (time(samples), pair(samples, realizations), linear(samples, realizations))
        allocate (temperature(samples, realizations), energy(samples, realizations))
        allocate (tnorm(samples, realizations), volterra(samples), landau(samples))
        allocate (bgross(samples), cold(samples), stat(samples, 6))
        allocate (fits(realizations, 6), cost(realizations, 6))
        allocate (sk(nshell, realizations), gcount(radial_bins, realizations))
        allocate (vectors(nshell), snapshots(realizations), rk(nshell), sdh(nshell))
        allocate (rr(radial_bins), gr(radial_bins, realizations), gstat(radial_bins, 3))
        allocate (skstat(nshell, 3), dh_g(radial_bins))
        allocate (profile(profile_bins, samples), xprofile(profile_bins))
        do j = 1, samples
            time(j) = (j - 1)*sample_interval
        end do
        call oracles()
        print '(a,i0,a,f6.3,a,f7.4,a,i0,a,i0)', 'C3 OCP: N = ', model%electrons, &
            ', k lambda_D = ', k, ', dt = ', dt, ', PME K = ', model%pme%mesh, &
            ', realizations = ', realizations
        call ensemble()
        call statistics()
        call write_csv()
        call write_result()
        call draw_figures()
        print '(a,2f9.5)', 'Landau root (omega, gamma):        ', real(root), -aimag(root)
        print '(a,2f9.5)', 'Volterra fit on window:            ', volterra_fit(1:2)
        print '(a,4f9.5)', 'ensemble pair fit omega +- , gamma +-:', fits_summary(1), &
            fits_summary(2), fits_summary(3), fits_summary(4)
        print '(a,2f9.5)', 'tangent Lyapunov rate (mean, spread):', lyapunov
    end subroutine

    subroutine oracles()
        real(dp), allocatable :: fine(:)
        real(dp) :: q, p, ql, qu, pl, pu, a, b, c(samples), s(samples), m11, m12, m22
        integer :: j, w1, w2
        call landau_root(k, root, have_root)
        if (.not. have_root) error stop 'Landau root did not converge'
        allocate (fine((samples - 1)*10 + 1))
        call vlasov_density_response(k, sample_interval/10, fine)
        volterra = fine(1::10)
        call bohm_gross(1.0_dp, k, 1.0_dp, omega_bg)
        bgross = cos(omega_bg*time)
        call initialize_slab_oscillation(slab, 4000, measure_time, 1.0_dp, 0.0_dp, ok)
        if (ok) call solve_slab_oscillation(slab, ok)
        if (.not. ok) error stop 'cold slab oracle failed'
        do j = 1, samples
            call sample_slab_oscillation(slab, min(time(j), measure_time), q, p, ql, qu, &
                                         pl, pu, ok)
            cold(j) = q
        end do
        w1 = window(fit_start); w2 = window(fit_end)
        call fit_damped_cosine(time(w1:w2), volterra(w1:w2), 0.5_dp*real(root), &
                               1.5_dp*real(root), volterra_fit(1), volterra_fit(2), &
                               volterra_fit(3))
        ! Landau mode: omega and damping from the root, amplitude/phase fitted
        ! to the Volterra response on the fit window.
        c = exp(aimag(root)*time)*cos(real(root)*time)
        s = exp(aimag(root)*time)*sin(real(root)*time)
        m11 = sum(c(w1:w2)**2); m12 = sum(c(w1:w2)*s(w1:w2)); m22 = sum(s(w1:w2)**2)
        a = (m22*sum(c(w1:w2)*volterra(w1:w2)) - m12*sum(s(w1:w2)*volterra(w1:w2)))/ &
            (m11*m22 - m12**2)
        b = (m11*sum(s(w1:w2)*volterra(w1:w2)) - m12*sum(c(w1:w2)*volterra(w1:w2)))/ &
            (m11*m22 - m12**2)
        landau = a*c + b*s
    end subroutine

    integer function window(t) result(j)
        real(dp), intent(in) :: t
        j = max(1, min(samples, 1 + nint(t/sample_interval)))
    end function

    subroutine ensemble()
        integer :: r
        !$omp parallel do schedule(dynamic, 1) default(shared) private(r)
        do r = 1, realizations
            call member(r)
        end do
        !$omp end parallel do
    end subroutine

    !> One ensemble member: preparation, equilibration, unperturbed,
    !> perturbed and tangent runs from the same equilibrated state.
    subroutine member(r)
        integer, intent(in) :: r
        type(ocp_model_t) :: m
        type(nbody_system_t) :: sys
        type(fortnum_status_t) :: st
        real(dp), allocatable :: base(:), s(:), d(:), au(:, :), hist(:, :), sk1(:)
        real(dp) :: a(3), e0, e, sampling, wall(4)
        integer(int64) :: c0, c1, c2, rate
        integer :: rejected, j, step, n3, nsel, i, vec(nshell), snap
        logical :: okm
        m = model
        call ocp_system(m, sys, okm)
        if (.not. okm) error stop 'member OCP system failed'
        n3 = 3*m%particles
        allocate (base(6*m%particles), s(6*m%particles), d(6*m%particles))
        allocate (au(3, samples), hist(profile_bins, samples), sk1(nshell))
        call system_clock(c0, rate)
        call ocp_prepare(m, seed + 7919_int64*r, base, rejected, st)
        do step = 1, nint(equilibration_time/m%dt)
            call ocp_advance(sys, base, m%dt, st)
            if (st%code /= FORTNUM_OK) error stop 'equilibration step failed'
            call ocp_thermostat(m, base)
        end do
        call system_clock(c1); wall(1) = real(c1 - c0, dp)/rate
        sampling = 0; sk(:, r) = 0; gcount(:, r) = 0; snap = 0
        ! Unperturbed run (diagnostics, snapshots, structure).
        s = base
        call nbody_energy(sys, s, e0, st)
        do j = 1, samples
            if (j > 1) call advance(per)
            call system_clock(c1)
            call modes(s, au(:, j))
            temperature(j, r) = ocp_kinetic_temperature(m, s)
            energy(j, r) = energy(max(1, j - 1), r)
            if (mod(j - 1, 10) == 0) then
                call nbody_energy(sys, s, e, st)
                energy(j, r) = (e - e0)/m%electrons
            end if
            if (r == 1) call histogram(s, hist(:, j))
            if (mod(j - 1, nint(snapshot_interval/sample_interval)) == 0) then
                snap = snap + 1
                call structure_factor(s, m%electrons, box, nint(sqrt(real(nshell))), &
                                      sk1, vec)
                sk(:, r) = sk(:, r) + sk1
                call pair_counts(s, m%electrons, box, rmax, gcount(:, r))
                call dump(s, j)
            end if
            call system_clock(c2); sampling = sampling + real(c2 - c1, dp)/rate
        end do
        snapshots(r) = snap; sk(:, r) = sk(:, r)/snap
        if (r == 1) vectors = vec
        call system_clock(c2); wall(2) = real(c2 - c0, dp)/rate - wall(1)
        ! Perturbed run with common random numbers.
        s = base
        call ocp_displace(m, s, amplitude, [.true., .true., .true.])
        if (r == 1) call select_slice(s, nsel)
        do j = 1, samples
            if (j > 1) call advance(per)
            call modes(s, a)
            pair(j, r) = sum(a - au(:, j))/3/(-amplitude)
            if (r == 1) then
                call histogram(s, xprofile)
                profile(:, j) = xprofile - hist(:, j)
                call store_slice(s, j)
            end if
        end do
        call system_clock(c1); wall(3) = real(c1 - c0, dp)/rate - sum(wall(1:2))
        ! Tangent with respect to eps at eps = 0 (cluster structure frozen).
        linear(:, r) = 0; tnorm(:, r) = 0
        if (tangent) then
            s = base
            call ocp_displace_tangent(m, s, [.true., .true., .true.], d)
            do j = 1, samples
                if (j > 1) then
                    do step = 1, per
                        call ocp_advance_jvp(sys, s, d, m%dt, st)
                        if (st%code /= FORTNUM_OK) error stop 'tangent step failed'
                    end do
                end if
                call mode_tangents(s, d, a)
                linear(j, r) = sum(a)/3/(-1.0_dp)
                tnorm(j, r) = norm2(d(:n3))
            end do
        end if
        call system_clock(c2); wall(4) = real(c2 - c0, dp)/rate - sum(wall(1:3))
        cost(r, 1:4) = wall; cost(r, 5) = sampling; cost(r, 6) = rejected
    contains
        subroutine advance(count)
            integer, intent(in) :: count
            integer :: q
            do q = 1, count
                call ocp_advance(sys, s, m%dt, st)
                if (st%code /= FORTNUM_OK) error stop 'OCP step failed: '//trim(st%msg)
            end do
        end subroutine

        subroutine modes(x, out)
            real(dp), intent(in) :: x(:)
            real(dp), intent(out) :: out(3)
            integer :: q
            real(dp) :: kv(3)
            do q = 1, 3
                kv = 0; kv(q) = k
                out(q) = 2*real(ocp_mode(m, x, kv))/m%electrons
            end do
        end subroutine

        subroutine mode_tangents(x, dx, out)
            real(dp), intent(in) :: x(:), dx(:)
            real(dp), intent(out) :: out(3)
            integer :: q, p
            out = 0
            do p = 1, m%electrons
                do q = 1, 3
                    out(q) = out(q) - sin(k*x(3*p - 3 + q))*k*dx(3*p - 3 + q)
                end do
            end do
            out = 2*out/m%electrons
        end subroutine

        !> Electron density along x in profile bins, relative to the mean.
        subroutine histogram(x, h)
            real(dp), intent(in) :: x(:)
            real(dp), intent(out) :: h(:)
            integer :: p, b
            h = 0
            do p = 1, m%electrons
                b = 1 + int((x(3*p - 2)/box - floor(x(3*p - 2)/box))*profile_bins)
                b = min(profile_bins, b); h(b) = h(b) + 1
            end do
            h = h/(real(m%electrons, dp)/profile_bins)
        end subroutine

        !> Explicit selection: electrons with |z - L/2| < slice at t = 0 of
        !> the perturbed run, followed by identity.
        subroutine select_slice(x, count)
            real(dp), intent(in) :: x(:)
            integer, intent(out) :: count
            real(dp) :: z
            integer :: p
            count = 0
            do p = 1, m%electrons
                z = x(3*p) - box*floor(x(3*p)/box)
                if (abs(z - box/2) < slice) count = count + 1
            end do
            selected = count
            allocate (sel_x(count, samples), sel_y(count, samples), sel_v(count, samples))
            allocate (picked(count))
            count = 0
            do p = 1, m%electrons
                z = x(3*p) - box*floor(x(3*p)/box)
                if (abs(z - box/2) < slice) then
                    count = count + 1; picked(count) = p
                end if
            end do
        end subroutine

        subroutine store_slice(x, j)
            real(dp), intent(in) :: x(:)
            integer, intent(in) :: j
            integer :: q, p
            do q = 1, selected
                p = picked(q)
                sel_x(q, j) = x(3*p - 2) - box*floor(x(3*p - 2)/box)
                sel_y(q, j) = x(3*p - 1) - box*floor(x(3*p - 1)/box)
                sel_v(q, j) = x(3*m%particles + 3*p - 2)
            end do
        end subroutine

        !> Snapshot for correlation analysis: stream file of N, t, x(3N), v(3N).
        subroutine dump(x, j)
            real(dp), intent(in) :: x(:)
            integer, intent(in) :: j
            integer :: u
            character(160) :: path
            write (path, '(a,a,a,i0,a,i0,a)') 'output/', trim(label), '-snapshot-r', r, &
                '-s', j, '.bin'
            open (newunit=u, file=trim(path), access='stream', form='unformatted', &
                  status='replace', action='write')
            write (u) m%electrons, time(j), m%debye_volume, box, x(:3*m%electrons), &
                x(3*m%particles + 1:3*m%particles + 3*m%electrons)
            close (u)
        end subroutine
    end subroutine

    !> Two-sided 95% Student-t quantile for nu degrees of freedom.
    pure real(dp) function student(nu) result(t)
        integer, intent(in) :: nu
        real(dp), parameter :: table(10) = [12.706_dp, 4.303_dp, 3.182_dp, 2.776_dp, &
            2.571_dp, 2.447_dp, 2.365_dp, 2.306_dp, 2.262_dp, 2.228_dp]
        if (nu <= 10) then
            t = table(max(1, nu))
        else if (nu <= 20) then
            t = 2.228_dp - (nu - 10)*(2.228_dp - 2.086_dp)/10
        else if (nu <= 30) then
            t = 2.086_dp - (nu - 20)*(2.086_dp - 2.042_dp)/10
        else
            t = 1.96_dp
        end if
    end function

    !> Mean and 95% half-width of the realization values x(:) (statistical).
    pure subroutine interval(x, mean, half)
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: mean, half
        integer :: n
        n = size(x); mean = sum(x)/n
        half = student(n - 1)*sqrt(sum((x - mean)**2)/(n - 1)/n)
    end subroutine

    subroutine statistics()
        real(dp) :: mean, half, slope, xs(samples), lyap(realizations)
        integer :: j, r, w1, w2, h
        do j = 1, samples
            call interval(pair(j, :), mean, half); stat(j, 1:3) = [mean, mean - half, mean + half]
            call interval(linear(j, :), mean, half)
            stat(j, 4:6) = [mean, mean - half, mean + half]
        end do
        w1 = window(fit_start); w2 = window(fit_end)
        do r = 1, realizations
            call fit_damped_cosine(time(w1:w2), pair(w1:w2, r), 0.5_dp*real(root), &
                1.5_dp*real(root), fits(r, 1), fits(r, 2), slope)
            fits(r, 3:4) = 0
            if (tangent) call fit_damped_cosine(time(w1:w2), linear(w1:w2, r), &
                0.5_dp*real(root), 1.5_dp*real(root), fits(r, 3), fits(r, 4), slope)
            fits(r, 5) = energy(samples, r); fits(r, 6) = sum(temperature(:, r))/samples
            ! Lyapunov-type growth of the position tangent over the second half.
            h = samples/2; xs = time
            slope = 0
            if (tangent) slope = (sum((xs(h:) - sum(xs(h:))/(samples - h + 1))* &
                log(tnorm(h:, r))))/sum((xs(h:) - sum(xs(h:))/(samples - h + 1))**2)
            lyap(r) = slope
        end do
        call interval(lyap, lyapunov(1), lyapunov(2))
        do j = 1, nshell
            rk(j) = two_pi*sqrt(real(j, dp))/box
            call debye_structure(rk(j)**2, 1.0_dp, sdh(j))
            call interval(sk(j, :), mean, half)
            skstat(j, :) = [mean, mean - half, mean + half]
        end do
        do r = 1, realizations
            call pair_distribution(gcount(:, r), model%electrons, box, rmax, &
                                   snapshots(r), rr, gr(:, r))
        end do
        do j = 1, radial_bins
            call interval(gr(j, :), mean, half)
            gstat(j, :) = [mean, mean - half, mean + half]
            dh_g(j) = 1 - exp(-rr(j))/(4*pi*debye_volume*rr(j))
        end do
    end subroutine

    !> Ensemble pair-fit [omega, 95% half-width, gamma, 95% half-width].
    real(dp) function fits_summary(i) result(v)
        integer, intent(in) :: i
        real(dp) :: mean, half
        call interval(fits(:, 1 + (i - 1)/2), mean, half)
        v = merge(mean, half, mod(i, 2) == 1)
    end function

    subroutine write_csv()
        integer :: u, j, r
        open (newunit=u, file='output/'//trim(label)//'-response.csv', status='replace')
        write (u, '(a)') '# C3 OCP plasma oscillation; statistical 95% Student-t '// &
            'intervals over realizations (not deterministic bounds)'
        write (u, '(a)') 't,paired_mean,paired_lo,paired_hi,tangent_mean,tangent_lo,'// &
            'tangent_hi,vlasov_volterra,landau_mode,bohm_gross,cold_slab,'// &
            'temperature_mean,energy_defect_per_particle_mean'
        do j = 1, samples
            write (u, '(*(g0.8,:,","))') time(j), stat(j, :), volterra(j), landau(j), &
                bgross(j), cold(j), sum(temperature(j, :))/realizations, &
                sum(energy(j, :))/realizations
        end do
        close (u)
        open (newunit=u, file='output/'//trim(label)//'-fits.csv', status='replace')
        write (u, '(a)') 'realization,omega_paired,gamma_paired,omega_tangent,'// &
            'gamma_tangent,energy_defect_per_particle,mean_temperature,'// &
            'equilibration_s,unperturbed_s,perturbed_s,tangent_s,sampling_s,rejected'
        do r = 1, realizations
            write (u, '(i0,",",*(g0.8,:,","))') r, fits(r, :), cost(r, :)
        end do
        close (u)
        open (newunit=u, file='output/'//trim(label)//'-structure.csv', status='replace')
        write (u, '(a)') 'k_lambda_d,s_mean,s_lo,s_hi,debye_hueckel,vectors'
        do j = 1, nshell
            if (vectors(j) > 0) write (u, '(*(g0.8,:,","))') rk(j), skstat(j, :), &
                sdh(j), real(vectors(j), dp)
        end do
        close (u)
        open (newunit=u, file='output/'//trim(label)//'-pair.csv', status='replace')
        write (u, '(a)') 'r_over_lambda_d,g_mean,g_lo,g_hi,debye_hueckel'
        do j = 1, radial_bins
            write (u, '(*(g0.8,:,","))') rr(j), gstat(j, :), dh_g(j)
        end do
        close (u)
    end subroutine

    subroutine write_result()
        character(*), parameter :: st = 'statistical 95% Student-t over realizations'
        call result%create(run)
        call result%attribute('model', 'classical electron OCP, periodic 3D3V, bare '// &
            'Coulomb (PME), encounter-resolved Strang step, no softening')
        call result%attribute('preparation', 'uniform positions, Maxwellian T=1, '// &
            'zero momentum, isokinetic thermostat during equilibration only')
        call result%write('time', time, '1/omega_p', 'sample time', dims='time')
        call result%write('response_paired', pair, '1', 'paired-difference response', &
            dims='time,realization', status='candidate per realization')
        call result%write('response_tangent', linear, '1', 'tangent response dR/deps', &
            dims='time,realization', status='candidate; chaotic growth limits validity')
        call result%write('response_statistics', stat, '1', &
            'paired mean/lo/hi, tangent mean/lo/hi', dims='time,statistic6', status=st)
        call result%write('vlasov_volterra', volterra, '1', &
            'linear Vlasov initial-value response', dims='time')
        call result%write('fits', fits, '1', &
            'omega/gamma paired, omega/gamma tangent, energy defect, temperature', &
            dims='realization,fit_column')
        call result%write('cost', cost, 's', &
            'wall s: equilibration, unperturbed, perturbed, tangent, sampling; rejected', &
            dims='realization,cost_column')
        call result%write('structure_factor', skstat, '1', 'S(k) mean/lo/hi', &
            dims='shell,statistic3', status=st)
        call result%write('pair_distribution', gstat, '1', 'g(r) mean/lo/hi', &
            dims='radial_bin,statistic3', status=st)
        call result%close()
    end subroutine

    subroutine draw_figures()
        real(dp) :: ylo, yhi, span, xc(profile_bins), pred(profile_bins), sc(nshell)
        real(dp) :: dlo, dhi
        integer :: j, m, frames
        logical :: has(nshell)
        ! Response against the continuum and analytic limits.
        call figure(figsize=[12.8_dp, 6.4_dp])
        call subplot(1, 2, 1)
        call fill_between(time, stat(:, 2), stat(:, 3), color=blue, alpha=0.25_dp)
        call plot(time, stat(:, 1), color=blue, linewidth=1.8_dp, &
                  label='microscopic, paired difference (95% stat.)')
        if (tangent) call plot(time, stat(:, 4), color=vermillion, linestyle='-.', &
                               linewidth=1.5_dp, label='microscopic step tangent $dR/d\epsilon$')
        call plot(time, cold, color=grey, linestyle=':', linewidth=1.2_dp, &
                  label='cold slab $\omega_p$')
        call plot(time, bgross, color=orange, linestyle=':', linewidth=1.2_dp, &
                  label='Bohm-Gross')
        call plot(time, landau, color=green, linestyle='--', linewidth=1.2_dp, &
                  label='Landau root mode')
        call plot(time, volterra, color=black, linestyle='--', linewidth=1.2_dp, &
                  label='linear Vlasov initial value')
        write (heading, '(a,f4.2)') 'Density mode response, $k\lambda_D$ = ', k
        call title(trim(heading))
        call xlabel('$\omega_p t$'); call ylabel('R(t) = (2/N) Re $\rho_k/(-\epsilon)$')
        call xlim(0.0_dp, measure_time); call ylim(-1.4_dp, 1.4_dp)
        call legend(loc='lower right')
        call subplot(1, 2, 2)
        call fill_between(time, stat(:, 2) - volterra, stat(:, 3) - volterra, color=blue, &
                          alpha=0.25_dp)
        call plot(time, stat(:, 1) - volterra, color=blue, linewidth=1.8_dp, &
                  label='paired - Vlasov (95% stat.)')
        if (tangent) call plot(time, stat(:, 4) - volterra, color=vermillion, &
                               linestyle='-.', linewidth=1.5_dp, label='tangent - Vlasov')
        call plot(time, landau - volterra, color=green, linestyle='--', linewidth=1.2_dp, &
                  label='Landau mode - Vlasov')
        span = max(0.05_dp, 1.1_dp*maxval(abs(stat(:, 2:3) - spread(volterra, 2, 2))))
        call title('Difference to the linear Vlasov response')
        call xlabel('$\omega_p t$'); call ylabel('$\Delta R$')
        call xlim(0.0_dp, measure_time); call ylim(-span, span)
        call legend(loc='lower left')
        write (heading, '(a,i0,a,f5.1,a,i0,a)') 'C3 periodic OCP, N = ', model%electrons, &
            ', $\Lambda$ = ', debye_volume, ', ', realizations, &
            ' realizations; bands statistical, not certified'
        call suptitle(trim(heading))
        call tight_layout()
        call savefig('output/'//trim(label)//'-response.png')
        ! Static structure against the Debye-Hueckel weak-coupling limit.
        has = vectors > 0
        m = count(has)
        call figure(figsize=[12.8_dp, 6.4_dp])
        call subplot(1, 2, 1)
        sc(:m) = pack(rk, has)
        call fill_between(sc(:m), pack(skstat(:, 2), has), pack(skstat(:, 3), has), &
                          color=blue, alpha=0.25_dp)
        call plot(sc(:m), pack(skstat(:, 1), has), color=blue, marker='o', &
                  linewidth=1.5_dp, label='microscopic S(k) (95% stat.)')
        call plot(sc(:m), pack(sdh, has), color=black, linestyle='--', linewidth=1.2_dp, &
                  label='Debye-Hückel $k^2/(k^2+k_D^2)$')
        call title('Static structure factor')
        call xlabel('$k\lambda_D$'); call ylabel('S(k)')
        call xlim(0.0_dp, 1.05_dp*maxval(sc(:m))); call ylim(0.0_dp, 1.4_dp)
        call legend(loc='lower right')
        call subplot(1, 2, 2)
        dlo = 4*pi*debye_volume
        call fill_between(rr, dlo*(gstat(:, 2) - 1), dlo*(gstat(:, 3) - 1), color=blue, &
                          alpha=0.25_dp)
        call plot(rr, dlo*(gstat(:, 1) - 1), color=blue, linewidth=1.5_dp, &
                  label='microscopic (95% stat.)')
        call plot(rr, dlo*(dh_g - 1), color=black, linestyle='--', linewidth=1.2_dp, &
                  label='Debye-Hückel -exp(-r)/r')
        call title('Pair correlation h(r) = g(r) - 1, scaled by $4\pi\Lambda$')
        call xlabel('$r/\lambda_D$'); call ylabel('$4\pi\Lambda h(r)$')
        dhi = max(1.0_dp, 1.1_dp*maxval(abs(dlo*(gstat(3:, 2:3) - 1))))
        call xlim(0.0_dp, rmax); call ylim(-dhi, 0.5_dp*dhi)
        call legend(loc='lower right')
        call suptitle('C3 equilibrium snapshots: candidate estimator, statistical bands')
        call tight_layout()
        call savefig('output/'//trim(label)//'-structure.png')
        ! Particles in an explicit thin slice and the paired density profile.
        do j = 1, profile_bins
            xc(j) = (j - 0.5_dp)*box/profile_bins
        end do
        view_profile = max(2*amplitude, 1.1_dp*maxval(abs(profile)))
        vmax_colour = 3
        call figure(figsize=[12.8_dp, 6.4_dp])
        canvas => get_global_figure()
        frame = (samples + 1)/2
        call draw_frame()
        call savefig('output/'//trim(label)//'-particles.png')
        frames = samples + 24
        animation = animation_func(update_frame, frames=frames, interval=42, fig=canvas)
        call save_animation(animation, 'output/'//trim(label)//'.mp4', fps=24, &
                            status=encode_status)
        if (encode_status /= 0) print '(a)', 'animation encoding failed'
        ylo = 0; yhi = 0; pred = 0
    contains
        subroutine update_frame(index)
            integer, intent(in) :: index
            frame = min(samples, index)
            call canvas%clear()
            call draw_frame()
            call canvas%set_rendered(.false.)
        end subroutine

        subroutine draw_frame()
            pred = -amplitude*volterra(frame)*cos(k*xc)
            call subplot(1, 2, 1)
            call scatter(sel_x(:, frame), sel_y(:, frame), c=sel_v(:, frame), &
                         color=blue, cmap='coolwarm', vmin=-vmax_colour, &
                         vmax=vmax_colour, s=9.0_dp, alpha=0.6_dp)
            call colorbar(label='$v_x/v_{th}$')
            write (heading, '(a,f6.2)') 'Electrons with |z - L/2| < 0.25 $\lambda_D$, $\omega_p t$ = ', &
                time(frame)
            call title(trim(heading))
            call xlabel('$x/\lambda_D$'); call ylabel('$y/\lambda_D$')
            call xlim(0.0_dp, box); call ylim(0.0_dp, box)
            call subplot(1, 2, 2)
            call plot(xc, profile(:, frame), color=blue, linewidth=1.8_dp, &
                      label='paired difference, one realization')
            call plot(xc, pred, color=black, linestyle='--', linewidth=1.2_dp, &
                      label='linear Vlasov $-\epsilon N(t)\cos(kx)$')
            write (heading, '(a,f6.2)') 'Density perturbation along x, $\omega_p t$ = ', time(frame)
            call title(trim(heading))
            call xlabel('$x/\lambda_D$'); call ylabel('$\delta n/n$')
            call xlim(0.0_dp, box); call ylim(-view_profile, view_profile)
            call legend(loc='upper right')
            call tight_layout()
        end subroutine
    end subroutine
end module ocp_oscillation_app

program ocp_oscillation_demo
    use ocp_oscillation_app, only: main
    implicit none
    call main()
end program ocp_oscillation_demo
