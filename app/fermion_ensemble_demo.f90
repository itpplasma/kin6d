program fermion_ensemble_demo
    ! Q3 evidence: finite-N ideal spin-1/2 fermions in a periodic box versus the
    ! grand-canonical thermodynamic limit and the classical ideal gas at the
    ! same density and temperature (a different preparation, not forced equal).
    ! Writes CSV and preliminary PNG figures to output/ (not committed).
    use, intrinsic :: iso_fortran_env, only: dp => real64, qp => real128
    use fortnum_interval, only: interval_t
    use kin6d_fermion_ensemble
    use kin6d_fermi_gas, only: fermi_workspace_t, fermi_state_t, &
        initialize_fermi_workspace, solve_fermi_gas, occupation_fermi_gas
    use kin6d_fermion_ensemble_leaves, only: fermion_cold_density_matrix, &
        fermion_cold_density_matrix_series, fermion_exchange_hole, &
        fermion_boltzmann_hole, fermion_classical_limits
    use fortplot, only: figure, subplot, plot, title, xlabel, ylabel, xlim, ylim, &
        legend, tight_layout, savefig, xscale, yscale, suptitle
    use fortplot_system_runtime, only: create_directory_runtime
    use kin6d_result_io, only: run_t, result_file_t
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp), parameter :: colors(3, 4) = reshape([0.0_dp, 0.447_dp, 0.698_dp, &
        0.835_dp, 0.369_dp, 0.0_dp, 0.0_dp, 0.62_dp, 0.451_dp, 0.8_dp, 0.475_dp, &
        0.655_dp], [3, 4])
    real(dp), parameter :: black(3) = [0.0_dp, 0.0_dp, 0.0_dp]
    real(dp), parameter :: gray(3) = [0.45_dp, 0.45_dp, 0.45_dp]
    ! Namelist &fermion_ensemble (theta = kT/E_F, N = particles per spin).
    integer :: workspace_cells = 16384
    ! Exchange hole: x = k_F r grid on [0, hole_x_max]; ground-state cutoff.
    integer :: hole_points = 181, hole_n = 123, ground_cutoff = 40
    real(dp) :: hole_thetas(3) = [0.25_dp, 1.0_dp, 2.0_dp], hole_x_max = 9.0_dp
    real(dp) :: pair_tolerance = 1.0e-9_dp
    integer :: occupation_sizes(2) = [33, 123]
    real(dp) :: occupation_thetas(2) = [0.25_dp, 1.0_dp]
    ! EOS: log theta grid eos_theta_min*eos_theta_ratio**[0,1]; the large
    ! size is solved on every eos_large_stride-th point.
    integer :: eos_points = 24, eos_fine = 60, eos_small = 33, eos_large = 123
    integer :: eos_large_stride = 3
    real(dp) :: eos_theta_min = 0.05_dp, eos_theta_ratio = 200.0_dp
    integer :: size_sizes(6) = [7, 19, 33, 57, 123, 251]
    real(dp) :: size_thetas(3) = [0.25_dp, 1.0_dp, 4.0_dp]
    namelist /fermion_ensemble/ workspace_cells, hole_points, hole_n, ground_cutoff, &
        hole_thetas, hole_x_max, pair_tolerance, occupation_sizes, occupation_thetas, &
        eos_points, eos_fine, eos_small, eos_large, eos_large_stride, eos_theta_min, &
        eos_theta_ratio, size_sizes, size_thetas
    type(run_t) :: run
    type(result_file_t) :: result
    type(fermi_workspace_t) :: work
    logical :: ok
    integer :: config_unit, ios
    character(256) :: msg
    ! Evidence gathered by the sections for the NetCDF result.
    real(dp), allocatable :: hole_x(:), hole_cold(:), hole_ground(:), hole_boltz(:)
    real(dp), allocatable :: hole_canonical(:, :), hole_gc(:, :)
    real(dp), allocatable :: occ_rows(:, :), size_rows(:, :)
    real(dp), allocatable :: eos_theta(:), eos_small_rows(:, :, :), eos_large_rows(:, :, :)
    real(dp), allocatable :: eos_small_cert(:), eos_large_cert(:), eos_done(:)
    real(dp), allocatable :: limit_theta(:), limit_values(:, :)

    call run%start('fermion_ensemble_demo')
    if (run%argument_count() > 0) error stop 'usage: fermion_ensemble_demo [--config FILE]'
    ios = 0; msg = ''
    call run%open_config(config_unit)
    if (run%has_config()) read (config_unit, nml=fermion_ensemble, iostat=ios, iomsg=msg)
    call run%close_config(config_unit, 'fermion_ensemble', ios, msg)
    if (workspace_cells < 1 .or. hole_points < 2 .or. hole_n < 1 .or. &
        ground_cutoff < 1 .or. eos_points < 2 .or. eos_fine < 2 .or. eos_small < 1 .or. &
        eos_large < 1 .or. eos_large_stride < 1 .or. any(occupation_sizes < 1) .or. &
        any(size_sizes < 1)) error stop 'fermion_ensemble sizes must be positive'
    if (.not. (all(hole_thetas > 0.0_dp) .and. all(occupation_thetas > 0.0_dp) .and. &
        all(size_thetas > 0.0_dp) .and. hole_x_max > 0.0_dp .and. &
        pair_tolerance >= 0.0_dp .and. eos_theta_min > 0.0_dp .and. &
        eos_theta_ratio > 1.0_dp)) error stop 'fermion_ensemble invalid theta/grid'
    call run%open_resolved(config_unit)
    write (config_unit, nml=fermion_ensemble)
    call run%close_resolved(config_unit)

    call create_directory_runtime('output', ok)
    call initialize_fermi_workspace(work, workspace_cells, ok)
    if (.not. ok) error stop 'workspace'
    allocate (occ_rows(8, 0), size_rows(8, 0))
    call hole_section()
    call occupation_section()
    call size_section()
    call eos_section()
    call write_result()

contains

    subroutine hole_section()
        type(fermion_spectrum_t) :: spec
        type(canonical_state_t) :: st
        type(fermi_state_t) :: fermi
        real(dp), allocatable :: x(:), r(:), cold(:), ground(:)
        real(dp), allocatable :: canonical(:, :), gc(:, :), boltz(:), rho(:)
        real(dp) :: ef, b, trunc, kfl, thetas(3), dev, span
        integer :: i, k, unit, points, n
        character(80) :: label
        points = hole_points; n = hole_n; thetas = hole_thetas
        allocate (x(points), r(points), cold(points), ground(points), &
            canonical(points, 3), gc(points, 3), boltz(points), rho(points))
        kfl = (6.0_dp*pi**2*real(n, dp))**(1.0_dp/3.0_dp)
        do i = 1, points
            x(i) = hole_x_max*real(i - 1, dp)/real(points - 1, dp)
            r(i) = x(i)/kfl
            if (x(i) < 0.05_dp) then
                call fermion_cold_density_matrix_series(x(i), rho(i))
            else
                call fermion_cold_density_matrix(x(i), rho(i))
            end if
            call fermion_exchange_hole(rho(i), cold(i))
            call fermion_boltzmann_hole(thetas(3), x(i), boltz(i))
        end do
        call box_spectrum(ground_cutoff, spec, ok)
        call solve_ground(spec, n, st, ok)
        call same_spin_pair(spec, st, r, ground, ok, 0.0_dp, trunc)
        if (.not. ok) error stop 'ground pair'
        print '(a,i4,a,es10.3)', 'N_sigma=', n, ' ground state: max|g-g_TL(T=0)|=', &
            maxval(abs(ground - cold))
        do k = 1, 3
            call box_parameters(2*n, thetas(k), ef, b)
            call box_spectrum(box_cutoff(b, n), spec, ok)
            call solve_canonical(spec, n, b, st, ok)
            if (.not. ok) error stop 'canonical state'
            call same_spin_pair(spec, st, r, canonical(:, k), ok, pair_tolerance, trunc)
            if (.not. ok) error stop 'canonical pair'
            call solve_fermi_gas(work, thetas(k), fermi, ok)
            call grand_canonical_density(fermi, x, rho, ok)
            if (.not. ok) error stop 'grand canonical density'
            do i = 1, points
                call fermion_exchange_hole(rho(i), gc(i, k))
            end do
            print '(a,f5.2,a,es10.3,a,es9.2,a,l2)', 'theta=', thetas(k), &
                ' max|g_N-g_GC|=', maxval(abs(canonical(:, k) - gc(:, k))), &
                ' truncation<=', trunc, ' certified state:', st%certified
        end do
        open (newunit=unit, file='output/fermion-exchange-hole.csv', status='replace')
        write (label, '(i0)') n
        write (unit, '(a)') '# same-spin g(r), N_sigma='//trim(label)//' per spin in a periodic box; '// &
            'x=kF*r; canonical=finite-N candidate; gc=grand-canonical thermodynamic '// &
            'limit; classical ideal gas g=1'
        write (unit, '(a)') 'x,analytic_T0,ground_N,'//trim(column('canonical_', &
            thetas(1)))//','//trim(column('gc_', thetas(1)))//','// &
            trim(column('canonical_', thetas(2)))//','//trim(column('gc_', thetas(2)))// &
            ','//trim(column('canonical_', thetas(3)))//','// &
            trim(column('gc_', thetas(3)))//','//trim(column('boltzmann_', thetas(3)))// &
            ',classical'
        do i = 1, points
            write (unit, '(*(es16.8,:,","))') x(i), cold(i), ground(i), &
                (canonical(i, k), gc(i, k), k=1, 3), boltz(i), 1.0_dp
        end do
        close (unit)
        hole_x = x; hole_cold = cold; hole_ground = ground; hole_boltz = boltz
        hole_canonical = canonical; hole_gc = gc
        call figure(figsize=[12.8_dp, 6.4_dp])
        call subplot(1, 2, 1)
        call plot([0.0_dp, hole_x_max], [1.0_dp, 1.0_dp], label='classical ideal gas', &
                  color=gray, linestyle=':')
        write (label, '(a,f4.2,a)') 'Boltzmann limit, $\Theta = ', thetas(3), '$'
        call plot(x, boltz, label=trim(label), color=colors(:, 4), linestyle=':', &
                  linewidth=2.0_dp)
        write (label, '(a,i0,a)') '$N = ', n, '$ ground state'
        call plot(x, ground, label=trim(label), color=gray, linewidth=2.0_dp)
        do k = 1, 3
            write (label, '(a,i0,a,f4.2,a)') '$N = ', n, '$, $\Theta = ', thetas(k), '$'
            call plot(x, canonical(:, k), label=trim(label), color=colors(:, k), &
                      linewidth=2.0_dp)
        end do
        ! Thermodynamic-limit references last, thin dashed black.
        call plot(x, cold, label='thermodynamic limit (grand canonical)', color=black, &
                  linestyle='--', linewidth=1.2_dp)
        do k = 1, 3
            call plot(x, gc(:, k), color=black, linestyle='--', linewidth=1.2_dp)
        end do
        call title('Same-spin pair correlation of free fermions')
        call xlabel('$k_F r$'); call ylabel('$g_{\sigma\sigma}(r)$')
        call xlim(0.0_dp, hole_x_max); call ylim(0.0_dp, 1.15_dp)
        call legend(loc='lower right')
        call subplot(1, 2, 2)
        write (label, '(a,i0,a)') '$N = ', n, '$ ground state $-$ $T = 0$ limit'
        call plot(x, ground - cold, label=trim(label), color=gray, linewidth=2.0_dp)
        do k = 1, 3
            write (label, '(a,f4.2,a)') '$\Theta = ', thetas(k), '$'
            call plot(x, canonical(:, k) - gc(:, k), label=trim(label), &
                      color=colors(:, k), linewidth=2.0_dp)
        end do
        call plot([0.0_dp, hole_x_max], [0.0_dp, 0.0_dp], color=black, linestyle='--', &
                  linewidth=1.2_dp)
        write (label, '(a,i0,a)') 'Finite-size deviation at $N = ', n, '$ per spin'
        call title(trim(label))
        call xlabel('$k_F r$'); call ylabel('$g_N - g_{\infty}$')
        dev = min(minval(ground - cold), minval(canonical - gc))
        span = max(maxval(ground - cold), maxval(canonical - gc)) - dev
        call xlim(0.0_dp, hole_x_max); call ylim(dev - 0.08_dp*span, dev + 1.5_dp*span)
        call legend(loc='upper right')
        call suptitle('Exchange hole without interaction (opposite spins: $g = 1$)')
        call tight_layout()
        call savefig('output/fermion-exchange-hole.png')
    end subroutine hole_section

    subroutine occupation_section()
        integer :: sizes(2)
        real(dp) :: thetas(2)
        type(fermion_spectrum_t) :: spec
        type(canonical_state_t) :: st
        type(fermi_state_t) :: fermi
        real(dp) :: ef, b, u(200), fd(200), mb(200), lower, upper, e, p, f, m, floor_occ
        real(dp), allocatable :: eps(:)
        integer :: i, j, k, unit, shells
        character(40) :: label
        sizes = occupation_sizes; thetas = occupation_thetas
        open (newunit=unit, file='output/fermion-occupations.csv', status='replace')
        write (unit, '(a)') '# per-orbital canonical occupation of box shells; '// &
            'fd=grand-canonical Fermi-Dirac at the certified thermodynamic-limit mu; '// &
            'mb=Maxwell-Boltzmann at the same density (can exceed one)'
        write (unit, '(a)') 'N_sigma,theta,eps_over_EF,occupation,lower,upper,fd,mb'
        call figure(figsize=[12.8_dp, 6.4_dp])
        do k = 1, 2
            call subplot(1, 2, k)
            call solve_fermi_gas(work, thetas(k), fermi, ok)
            do i = 1, 200
                u(i) = 2.0_dp*real(i - 1, dp)/199.0_dp
                call occupation_fermi_gas(fermi, u(i), fd(i), lower, upper, ok)
                call fermion_classical_limits(thetas(k), u(i), e, p, mb(i))
            end do
            ! Maxwell-Boltzmann exceeds one at low energy for Theta < 1; show it
            ! only inside the axes (fortplot does not clip lines).
            floor_occ = 0.5_dp*minval(fd)
            call plot(pack(u**2, mb <= 1.5_dp .and. mb >= floor_occ), &
                      pack(mb, mb <= 1.5_dp .and. mb >= floor_occ), &
                      label='Maxwell-Boltzmann, same $n$ and $T$', color=gray, &
                      linestyle=':', linewidth=2.0_dp)
            do j = 1, 2
                call box_parameters(2*sizes(j), thetas(k), ef, b)
                call box_spectrum(box_cutoff(b, sizes(j)), spec, ok)
                call solve_canonical(spec, sizes(j), b, st, ok)
                shells = count(spec%level/ef <= 4.0_dp)
                eps = spec%level(:shells)/ef
                do i = 1, shells
                    call occupation_fermi_gas(fermi, sqrt(eps(i)), f, lower, upper, ok)
                    call fermion_classical_limits(thetas(k), sqrt(eps(i)), e, p, m)
                    write (unit, '(i4,",",f6.3,6(",",es16.8))') sizes(j), thetas(k), &
                        eps(i), st%occupation(i), st%occupation_bound(i)%lo, &
                        st%occupation_bound(i)%hi, f, m
                    occ_rows = reshape([occ_rows, real(sizes(j), dp), thetas(k), eps(i), &
                        st%occupation(i), st%occupation_bound(i)%lo, &
                        st%occupation_bound(i)%hi, f, m], [8, size(occ_rows, 2) + 1])
                end do
                write (label, '(a,i0,a)') 'canonical, $N = ', sizes(j), '$ per spin'
                call plot(eps, st%occupation(:shells), label=trim(label), &
                          color=colors(:, j), linestyle='None', marker='o')
            end do
            call plot(u**2, fd, label='Fermi-Dirac (thermodynamic limit)', color=black, &
                      linestyle='--', linewidth=1.2_dp)
            write (label, '(a,f4.2,a)') 'Orbital occupation at $\Theta = ', thetas(k), '$'
            call title(trim(label)); call xlabel('$\epsilon/E_F$')
            call ylabel('occupation per orbital')
            call yscale('log')
            call xlim(0.0_dp, 4.0_dp); call ylim(floor_occ, 1.5_dp)
            call legend(loc='lower left')
        end do
        call suptitle('Finite-N canonical occupations against the Fermi-Dirac '// &
                      'and Maxwell-Boltzmann laws')
        close (unit)
        call tight_layout()
        call savefig('output/fermion-occupations.png')
    end subroutine occupation_section

    subroutine eos_section()
        type(fermi_state_t) :: fermi
        real(dp), allocatable :: theta(:), tl(:, :), tt(:)
        real(dp), allocatable :: small(:, :), large(:, :)
        logical, allocatable :: done(:)
        integer :: i, unit, points, fine
        character(24) :: label
        points = eos_points; fine = eos_fine
        allocate (theta(points), tl(fine, 4), tt(fine), small(points, 5), &
            large(points, 5), done(points), eos_small_rows(points, 4, 3), &
            eos_large_rows(points, 4, 3), eos_small_cert(points), &
            eos_large_cert(points))
        large = 0.0_dp; eos_large_rows = 0.0_dp; eos_large_cert = 0.0_dp
        open (newunit=unit, file='output/fermion-eos-scan.csv', status='replace')
        write (unit, '(a)') '# box EOS per spin-1/2 electron in units of the '// &
            'thermodynamic-limit EF; e=E/(N EF), p=P/(n EF)=2e/3, de/dtheta=C_V/(N kB), '// &
            'B=n dP/dn|T/(n EF); certified=outward enclosure within 1% budget'
        write (unit, '(a)') 'N_sigma,theta,e,e_lo,e_hi,p,p_lo,p_hi,dedtheta,'// &
            'dedtheta_lo,dedtheta_hi,bulk,bulk_lo,bulk_hi,certified'
        done = .false.
        do i = 1, points
            theta(i) = eos_theta_min*(eos_theta_ratio)**(real(i - 1, dp)/real(points - 1, dp))
            call row(unit, eos_small, theta(i), small(i, :), eos_small_rows(i, :, :), &
                eos_small_cert(i))
            if (mod(i - 1, eos_large_stride) == 0) then
                call row(unit, eos_large, theta(i), large(i, :), eos_large_rows(i, :, :), &
                    eos_large_cert(i))
                done(i) = .true.
            end if
        end do
        close (unit)
        open (newunit=unit, file='output/fermion-eos-limit.csv', status='replace')
        write (unit, '(a)') '# certified homogeneous thermodynamic-limit Fermi gas'
        write (unit, '(a)') 'theta,e,p,dedtheta,bulk'
        do i = 1, fine
            tt(i) = eos_theta_min*(eos_theta_ratio)**(real(i - 1, dp)/real(fine - 1, dp))
            call solve_fermi_gas(work, tt(i), fermi, ok)
            tl(i, :) = [fermi%energy, fermi%pressure, fermi%theta_gradient(3), &
                        fermi%bulk_modulus]
            write (unit, '(*(es16.8,:,","))') tt(i), tl(i, :)
        end do
        close (unit)
        eos_theta = theta; eos_done = merge(1.0_dp, 0.0_dp, done)
        limit_theta = tt; limit_values = tl
        call figure(figsize=[13.2_dp, 8.8_dp])
        call subplot(2, 2, 1)
        call plot([tt(1), tt(fine)], [1.0_dp, 1.0_dp], label='classical $P = n k_B T$', &
                  color=gray, linestyle='--')
        call plot(pack(tt, 0.4_dp/tt >= 0.5_dp), pack(0.4_dp/tt, 0.4_dp/tt >= 0.5_dp), &
                  label='$T = 0$ Fermi pressure $\frac{2}{5} n E_F$', &
                  color=gray, linestyle=':', linewidth=2.0_dp)
        write (label, '(a,i0,a)') '$N = ', eos_small, '$ per spin'
        call plot(theta, small(:, 2)/theta, label=trim(label), color=colors(:, 1), &
                  linestyle='None', marker='o')
        write (label, '(a,i0,a)') '$N = ', eos_large, '$ per spin'
        call plot(pack(theta, done), pack(large(:, 2), done)/pack(theta, done), &
                  label=trim(label), color=colors(:, 2), linestyle='None', marker='s')
        call plot(tt, tl(:, 2)/tt, label='thermodynamic limit', color=black, &
                  linestyle='--', linewidth=1.2_dp)
        call xscale('log'); call yscale('log')
        call title('Pressure of the ideal Fermi gas'); call xlabel('$\Theta = k_B T/E_F$')
        call ylabel('$P/(n k_B T)$')
        call xlim(0.8_dp*tt(1), 1.25_dp*tt(fine))
        call ylim(0.5_dp, 2.0_dp*maxval(tl(:, 2)/tt))
        call legend(loc='upper right')
        call subplot(2, 2, 2)
        call plot([tt(1), tt(fine)], [1.5_dp, 1.5_dp], label='classical $3/2$', &
                  color=gray, linestyle='--')
        write (label, '(a,i0,a)') '$N = ', eos_small, '$ per spin'
        call plot(theta, small(:, 3), label=trim(label), color=colors(:, 1), &
                  linestyle='None', marker='o')
        write (label, '(a,i0,a)') '$N = ', eos_large, '$ per spin'
        call plot(pack(theta, done), pack(large(:, 3), done), label=trim(label), &
                  color=colors(:, 2), linestyle='None', marker='s')
        call plot(tt, tl(:, 3), label='thermodynamic limit', color=black, &
                  linestyle='--', linewidth=1.2_dp)
        call xscale('log')
        call title('Heat capacity $C_V/(N k_B) = de/d\Theta$')
        call xlabel('$\Theta = k_B T/E_F$'); call ylabel('$C_V/(N k_B)$')
        call xlim(0.8_dp*tt(1), 1.25_dp*tt(fine)); call ylim(0.0_dp, 1.8_dp)
        call legend(loc='lower right')
        call size_panel(3, 4, 'energy', '$|e_N - e_\infty|$ (per particle, units of $E_F$)')
        call size_panel(4, 6, 'heat capacity', '$|C_N - C_\infty|/(N k_B)$')
        call suptitle('Ideal Fermi gas: finite-$N$ canonical box against the '// &
                      'thermodynamic limit')
        call tight_layout()
        call savefig('output/fermion-eos-scan.png')

    end subroutine eos_section

    subroutine row(unit, n, th, out, evidence, certified)
        integer, intent(in) :: unit, n
        real(dp), intent(in) :: th
        real(dp), intent(out) :: out(5), evidence(4, 3), certified
        type(fermion_spectrum_t) :: spec
        type(canonical_state_t) :: st
        type(interval_t) :: bounds(5)
        real(dp) :: ef, b
        logical :: cert
        call box_parameters(2*n, th, ef, b)
        call box_spectrum(box_cutoff(b, n), spec, ok)
        call solve_canonical(spec, n, b, st, ok)
        call box_eos(st, st, out, bounds, ok, cert)
        if (.not. ok) error stop 'box EOS'
        write (unit, '(i4,",",es14.6,12(",",es16.8),",",l1)') n, th, out(1), &
            bounds(1)%lo, bounds(1)%hi, out(2), bounds(2)%lo, bounds(2)%hi, &
            out(3), bounds(3)%lo, bounds(3)%hi, out(5), bounds(5)%lo, bounds(5)%hi, cert
        evidence(:, 1) = out([1, 2, 3, 5])
        evidence(:, 2) = bounds([1, 2, 3, 5])%lo
        evidence(:, 3) = bounds([1, 2, 3, 5])%hi
        certified = merge(1.0_dp, 0.0_dp, cert)
    end subroutine row

    subroutine size_section()
        integer :: sizes(6)
        real(dp) :: thetas(3)
        type(fermion_spectrum_t) :: spec
        type(canonical_state_t) :: st
        type(fermi_state_t) :: fermi
        real(dp) :: values(5), ef, b
        type(interval_t) :: bounds(5)
        real(qp) :: log_z, cancellation
        logical :: cert
        integer :: i, k, unit
        sizes = size_sizes; thetas = size_thetas
        open (newunit=unit, file='output/fermion-size-scan.csv', status='replace')
        write (unit, '(a)') '# finite-N minus thermodynamic limit (indicator, not a '// &
            'bound); cancellation=alternating power-sum recursion factor (binary128)'
        write (unit, '(a)') 'N_sigma,theta,e,e_minus_limit,p_minus_limit,'// &
            'dedtheta_minus_limit,certified,power_sum_cancellation'
        do k = 1, 3
            call solve_fermi_gas(work, thetas(k), fermi, ok)
            do i = 1, 6
                if (thetas(k) > 2.0_dp .and. sizes(i) > 123) cycle
                call box_parameters(2*sizes(i), thetas(k), ef, b)
                call box_spectrum(box_cutoff(b, sizes(i)), spec, ok)
                call solve_canonical(spec, sizes(i), b, st, ok)
                call box_eos(st, st, values, bounds, ok, cert)
                call power_sum_recursion(spec, sizes(i), b, log_z, cancellation, ok)
                write (unit, '(i4,",",f6.3,4(",",es14.6),",",l1,",",es10.3)') sizes(i), &
                    thetas(k), values(1), values(1) - fermi%energy, values(2) - &
                    fermi%pressure, values(3) - fermi%theta_gradient(3), cert, &
                    real(cancellation, dp)
                size_rows = reshape([size_rows, real(sizes(i), dp), thetas(k), values(1), &
                    values(1) - fermi%energy, values(2) - fermi%pressure, &
                    values(3) - fermi%theta_gradient(3), merge(1.0_dp, 0.0_dp, cert), &
                    real(cancellation, dp)], [8, size(size_rows, 2) + 1])
                print '(a,i4,a,f5.2,a,es10.3,a,es9.2)', 'N_sigma=', sizes(i), ' theta=', &
                    thetas(k), ' e_N-e_TL=', values(1) - fermi%energy, &
                    ' power-sum cancellation=', real(cancellation, dp)
            end do
        end do
        close (unit)
    end subroutine size_section

    !> Finite-size indicator |X_N - X_inf| against N per spin (log-log), one
    !> curve per size_thetas entry, from the size-scan rows (column col).
    subroutine size_panel(position, col, name, axis_label)
        integer, intent(in) :: position, col
        character(len=*), intent(in) :: name, axis_label
        character(len=80) :: label
        logical, allocatable :: sel(:)
        real(dp) :: lo, hi
        integer :: k
        character(len=1), parameter :: markers(3) = ['o', 's', '^']
        call subplot(2, 2, position)
        lo = huge(1.0_dp); hi = 0.0_dp
        do k = 1, 3
            sel = abs(size_rows(2, :) - size_thetas(k)) < 1.0e-12_dp .and. &
                  abs(size_rows(col, :)) > 0.0_dp
            if (.not. any(sel)) cycle
            write (label, '(a,f4.2,a)') '$\Theta = ', size_thetas(k), '$'
            call plot(pack(size_rows(1, :), sel), abs(pack(size_rows(col, :), sel)), &
                      label=trim(label), color=colors(:, k), marker=markers(k))
            lo = min(lo, minval(abs(pack(size_rows(col, :), sel))))
            hi = max(hi, maxval(abs(pack(size_rows(col, :), sel))))
        end do
        call xscale('log'); call yscale('log')
        call title('Finite-size '//name//' deviation (indicator)')
        call xlabel('particles per spin $N$'); call ylabel(axis_label)
        call xlim(0.8_dp*minval(size_rows(1, :)), 1.25_dp*maxval(size_rows(1, :)))
        call ylim(0.5_dp*lo, 3.0_dp*hi)
        call legend(loc='lower left')
    end subroutine size_panel

    !> CSV column name prefix//theta with trailing zeros removed (0.25, 1, 2).
    function column(prefix, theta) result(name)
        character(len=*), intent(in) :: prefix
        real(dp), intent(in) :: theta
        character(len=48) :: name
        character(len=32) :: text
        integer :: last
        write (text, '(f0.6)') theta
        if (text(1:1) == '.') text = '0'//text
        last = len_trim(text)
        do while (text(last:last) == '0')
            last = last - 1
        end do
        if (text(last:last) == '.') last = last - 1
        name = prefix//text(:last)
    end function column

    subroutine write_result()
        character(len=*), parameter :: bound = 'certified-discrete outward enclosure '// &
            'of the finite-N box state where certified=1'
        call result%create(run)
        call result%attribute('scope', 'ideal spin-1/2 fermions, no interaction; '// &
            'theta = kT/E_F; finite-N minus limit differences are indicators, not bounds')
        call result%write('hole_x', hole_x, '1', 'k_F r', dims='hole_x')
        call result%write('hole_thetas', hole_thetas, '1', 'theta of the canonical holes', &
            dims='hole_theta')
        call result%write('hole_analytic_cold', hole_cold, '1', &
            'T=0 thermodynamic-limit same-spin g', dims='hole_x', status='reference limit')
        call result%write('hole_ground', hole_ground, '1', 'finite-N ground-state g', &
            dims='hole_x', status='candidate')
        call result%write('hole_canonical', hole_canonical, '1', 'finite-N canonical g', &
            dims='hole_x,hole_theta', status='candidate')
        call result%write('hole_grand_canonical', hole_gc, '1', &
            'grand-canonical thermodynamic-limit g', dims='hole_x,hole_theta', &
            status='candidate')
        call result%write('hole_boltzmann', hole_boltz, '1', &
            'Boltzmann limit at hole_thetas(3)', dims='hole_x', status='reference limit')
        call result%write('occupation_rows', occ_rows, '1', 'rows: N_sigma, theta, '// &
            'eps/E_F, occupation, lower, upper, Fermi-Dirac, Maxwell-Boltzmann', &
            dims='occupation_field,occupation_row', status='occupation lower/upper: '// &
            'certified-discrete canonical enclosure; fd/mb candidate references')
        call result%write('eos_theta', eos_theta, '1', 'theta grid', dims='eos_theta')
        call result%write('eos_small', eos_small_rows, '1', 'N=eos_small box EOS [e,p,'// &
            'de/dtheta,bulk] x [center,lower,upper] in E_F units', &
            dims='eos_theta,eos_quantity,eos_bound', status=bound)
        call result%write('eos_small_certified', eos_small_cert, '1', &
            '1 if within the 1% budget', dims='eos_theta')
        call result%write('eos_large', eos_large_rows, '1', 'N=eos_large box EOS (zero '// &
            'where eos_large_done=0)', dims='eos_theta,eos_quantity,eos_bound', status=bound)
        call result%write('eos_large_certified', eos_large_cert, '1', &
            '1 if within the 1% budget', dims='eos_theta')
        call result%write('eos_large_done', eos_done, '1', '1 where N=eos_large was solved', &
            dims='eos_theta')
        call result%write('limit_theta', limit_theta, '1', 'theta grid', dims='limit_theta')
        call result%write('limit_eos', limit_values, '1', 'thermodynamic limit [e,p,'// &
            'de/dtheta,bulk]', dims='limit_theta,limit_quantity', &
            status='candidate centers of the certified limit solve; bounds not recorded')
        call result%write('size_rows', size_rows, '1', 'rows: N_sigma, theta, e, '// &
            'e-e_TL, p-p_TL, dedtheta-limit, certified, power-sum cancellation', &
            dims='size_field,size_row', status='indicator, not a bound')
        call result%close()
    end subroutine write_result
end program fermion_ensemble_demo
