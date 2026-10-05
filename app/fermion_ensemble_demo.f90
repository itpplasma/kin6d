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
        legend, tight_layout, savefig, xscale, yscale
    use fortplot_system_runtime, only: create_directory_runtime
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp), parameter :: colors(3, 4) = reshape([0.0_dp, 0.447_dp, 0.698_dp, &
        0.835_dp, 0.369_dp, 0.0_dp, 0.0_dp, 0.62_dp, 0.451_dp, 0.8_dp, 0.475_dp, &
        0.655_dp], [3, 4])
    type(fermi_workspace_t) :: work
    logical :: ok

    call create_directory_runtime('output', ok)
    call initialize_fermi_workspace(work, 16384, ok)
    if (.not. ok) error stop 'workspace'
    call hole_section()
    call occupation_section()
    call eos_section()
    call size_section()

contains

    subroutine hole_section()
        integer, parameter :: points = 181, n = 123
        real(dp), parameter :: thetas(3) = [0.25_dp, 1.0_dp, 2.0_dp]
        type(fermion_spectrum_t) :: spec
        type(canonical_state_t) :: st
        type(fermi_state_t) :: fermi
        real(dp) :: x(points), r(points), cold(points), ground(points)
        real(dp) :: canonical(points, 3), gc(points, 3), boltz(points), rho(points)
        real(dp) :: ef, b, trunc, kfl
        integer :: i, k, unit
        character(32) :: label
        kfl = (6.0_dp*pi**2*real(n, dp))**(1.0_dp/3.0_dp)
        do i = 1, points
            x(i) = 9.0_dp*real(i - 1, dp)/real(points - 1, dp)
            r(i) = x(i)/kfl
            if (x(i) < 0.05_dp) then
                call fermion_cold_density_matrix_series(x(i), rho(i))
            else
                call fermion_cold_density_matrix(x(i), rho(i))
            end if
            call fermion_exchange_hole(rho(i), cold(i))
            call fermion_boltzmann_hole(thetas(3), x(i), boltz(i))
        end do
        call box_spectrum(40, spec, ok)
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
            call same_spin_pair(spec, st, r, canonical(:, k), ok, 1.0e-9_dp, trunc)
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
        write (unit, '(a)') '# same-spin g(r), N_sigma=123 per spin in a periodic box; '// &
            'x=kF*r; canonical=finite-N candidate; gc=grand-canonical thermodynamic '// &
            'limit; classical ideal gas g=1'
        write (unit, '(a)') 'x,analytic_T0,ground_N,canonical_0.25,gc_0.25,'// &
            'canonical_1,gc_1,canonical_2,gc_2,boltzmann_2,classical'
        do i = 1, points
            write (unit, '(*(es16.8,:,","))') x(i), cold(i), ground(i), &
                (canonical(i, k), gc(i, k), k=1, 3), boltz(i), 1.0_dp
        end do
        close (unit)
        call figure(figsize=[8.0_dp, 5.5_dp])
        call plot(x, cold, label='T=0 limit 1-9(j1(x)/x)^2', color=[0.0_dp, 0.0_dp, &
                  0.0_dp], linewidth=2.0_dp)
        call plot(x, ground, label='N=123 ground state', color=[0.5_dp, 0.5_dp, &
                  0.5_dp], linestyle='--')
        do k = 1, 3
            write (label, '(a,f4.2,a)') 'N=123, Theta=', thetas(k), ' canonical'
            call plot(x, canonical(:, k), label=trim(label), color=colors(:, k))
            write (label, '(a,f4.2,a)') 'Theta=', thetas(k), ' grand canonical'
            call plot(x, gc(:, k), label=trim(label), color=colors(:, k), linestyle='--')
        end do
        call plot(x, boltz, label='Boltzmann limit 1-exp(-Theta x^2/2), Theta=2', &
                  color=colors(:, 4), linestyle=':')
        call plot([0.0_dp, 9.0_dp], [1.0_dp, 1.0_dp], label='classical ideal gas g=1', &
                  color=[0.3_dp, 0.3_dp, 0.3_dp], linestyle=':')
        call title('Same-spin exchange hole without interaction (opposite spin: g=1)')
        call xlabel('k_F r'); call ylabel('g_{sigma sigma}(r)')
        call xlim(0.0_dp, 9.0_dp); call ylim(0.0_dp, 1.2_dp)
        call legend(loc='lower right')
        call savefig('output/fermion-exchange-hole.png')
    end subroutine hole_section

    subroutine occupation_section()
        integer, parameter :: sizes(2) = [33, 123]
        real(dp), parameter :: thetas(2) = [0.25_dp, 1.0_dp]
        type(fermion_spectrum_t) :: spec
        type(canonical_state_t) :: st
        type(fermi_state_t) :: fermi
        real(dp) :: ef, b, u(200), fd(200), mb(200), lower, upper, e, p, f, m
        real(dp), allocatable :: eps(:)
        integer :: i, j, k, unit, shells
        character(40) :: label
        open (newunit=unit, file='output/fermion-occupations.csv', status='replace')
        write (unit, '(a)') '# per-orbital canonical occupation of box shells; '// &
            'fd=grand-canonical Fermi-Dirac at the certified thermodynamic-limit mu; '// &
            'mb=Maxwell-Boltzmann at the same density (can exceed one)'
        write (unit, '(a)') 'N_sigma,theta,eps_over_EF,occupation,lower,upper,fd,mb'
        call figure(figsize=[11.0_dp, 4.8_dp])
        do k = 1, 2
            call subplot(1, 2, k)
            call solve_fermi_gas(work, thetas(k), fermi, ok)
            do i = 1, 200
                u(i) = 2.0_dp*real(i - 1, dp)/199.0_dp
                call occupation_fermi_gas(fermi, u(i), fd(i), lower, upper, ok)
                call fermion_classical_limits(thetas(k), u(i), e, p, mb(i))
            end do
            call plot(u**2, fd, label='Fermi-Dirac (thermodynamic limit)', &
                      color=[0.0_dp, 0.0_dp, 0.0_dp])
            call plot(u**2, mb, label='Maxwell-Boltzmann, same n and T', &
                      color=[0.4_dp, 0.4_dp, 0.4_dp], linestyle='--')
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
                end do
                write (label, '(a,i3,a)') 'canonical N=', sizes(j), ' per spin'
                call plot(eps, st%occupation(:shells), label=trim(label), &
                          color=colors(:, j), linestyle='None', marker='o')
            end do
            write (label, '(a,f4.2)') 'Theta=', thetas(k)
            call title(trim(label)); call xlabel('epsilon/E_F'); call ylabel('occupation')
            call xlim(0.0_dp, 4.0_dp); call ylim(0.0_dp, 1.6_dp); call legend()
        end do
        close (unit)
        call tight_layout()
        call savefig('output/fermion-occupations.png')
    end subroutine occupation_section

    subroutine eos_section()
        integer, parameter :: points = 24, fine = 60
        type(fermi_state_t) :: fermi
        real(dp) :: theta(points), tl(fine, 4), tt(fine)
        real(dp) :: small(points, 5), large(points, 5)
        logical :: done(points)
        integer :: i, unit
        open (newunit=unit, file='output/fermion-eos-scan.csv', status='replace')
        write (unit, '(a)') '# box EOS per spin-1/2 electron in units of the '// &
            'thermodynamic-limit EF; e=E/(N EF), p=P/(n EF)=2e/3, de/dtheta=C_V/(N kB), '// &
            'B=n dP/dn|T/(n EF); certified=outward enclosure within 1% budget'
        write (unit, '(a)') 'N_sigma,theta,e,e_lo,e_hi,p,p_lo,p_hi,dedtheta,'// &
            'dedtheta_lo,dedtheta_hi,bulk,bulk_lo,bulk_hi,certified'
        done = .false.
        do i = 1, points
            theta(i) = 0.05_dp*(200.0_dp)**(real(i - 1, dp)/real(points - 1, dp))
            call row(unit, 33, theta(i), small(i, :))
            if (mod(i - 1, 3) == 0) then
                call row(unit, 123, theta(i), large(i, :))
                done(i) = .true.
            end if
        end do
        close (unit)
        open (newunit=unit, file='output/fermion-eos-limit.csv', status='replace')
        write (unit, '(a)') '# certified homogeneous thermodynamic-limit Fermi gas'
        write (unit, '(a)') 'theta,e,p,dedtheta,bulk'
        do i = 1, fine
            tt(i) = 0.05_dp*(200.0_dp)**(real(i - 1, dp)/real(fine - 1, dp))
            call solve_fermi_gas(work, tt(i), fermi, ok)
            tl(i, :) = [fermi%energy, fermi%pressure, fermi%theta_gradient(3), &
                        fermi%bulk_modulus]
            write (unit, '(*(es16.8,:,","))') tt(i), tl(i, :)
        end do
        close (unit)
        call figure(figsize=[11.0_dp, 4.8_dp])
        call subplot(1, 2, 1)
        call plot(tt, tl(:, 2)/tt, label='thermodynamic limit', color=[0.0_dp, &
                  0.0_dp, 0.0_dp])
        call plot(theta, small(:, 2)/theta, label='N=33 per spin', color=colors(:, 1), &
                  linestyle='None', marker='o')
        call plot(pack(theta, done), pack(large(:, 2), done)/pack(theta, done), &
                  label='N=123 per spin', color=colors(:, 2), linestyle='None', marker='s')
        call plot([0.05_dp, 10.0_dp], [1.0_dp, 1.0_dp], label='classical P=nkT', &
                  color=[0.4_dp, 0.4_dp, 0.4_dp], linestyle='--')
        call plot(tt, 0.4_dp/tt, label='T=0 Fermi pressure 0.4 n E_F', &
                  color=[0.4_dp, 0.4_dp, 0.4_dp], linestyle=':')
        call xscale('log')
        call title('Pressure without interaction'); call xlabel('Theta = kT/E_F')
        call ylabel('P/(n k_B T)'); call ylim(0.0_dp, 9.0_dp)
        call legend(loc='upper right')
        call subplot(1, 2, 2)
        call plot(tt, tl(:, 3), label='thermodynamic limit', color=[0.0_dp, 0.0_dp, &
                  0.0_dp])
        call plot(theta, small(:, 3), label='N=33 per spin', color=colors(:, 1), &
                  linestyle='None', marker='o')
        call plot(pack(theta, done), pack(large(:, 3), done), label='N=123 per spin', &
                  color=colors(:, 2), linestyle='None', marker='s')
        call plot([0.05_dp, 10.0_dp], [1.5_dp, 1.5_dp], label='classical 3/2', &
                  color=[0.4_dp, 0.4_dp, 0.4_dp], linestyle='--')
        call xscale('log')
        call title('Heat capacity (derivative de/dTheta)'); call xlabel('Theta = kT/E_F')
        call ylabel('C_V/(N k_B)'); call ylim(0.0_dp, 1.8_dp)
        call legend(loc='lower right')
        call tight_layout()
        call savefig('output/fermion-eos-scan.png')

    end subroutine eos_section

    subroutine row(unit, n, th, out)
        integer, intent(in) :: unit, n
        real(dp), intent(in) :: th
        real(dp), intent(out) :: out(5)
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
    end subroutine row

    subroutine size_section()
        integer, parameter :: sizes(6) = [7, 19, 33, 57, 123, 251]
        real(dp), parameter :: thetas(3) = [0.25_dp, 1.0_dp, 4.0_dp]
        type(fermion_spectrum_t) :: spec
        type(canonical_state_t) :: st
        type(fermi_state_t) :: fermi
        real(dp) :: values(5), ef, b
        type(interval_t) :: bounds(5)
        real(qp) :: log_z, cancellation
        logical :: cert
        integer :: i, k, unit
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
                print '(a,i4,a,f5.2,a,es10.3,a,es9.2)', 'N_sigma=', sizes(i), ' theta=', &
                    thetas(k), ' e_N-e_TL=', values(1) - fermi%energy, &
                    ' power-sum cancellation=', real(cancellation, dp)
            end do
        end do
        close (unit)
    end subroutine size_section
end program fermion_ensemble_demo
