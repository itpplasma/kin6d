!> Q7 part 1 demo: Debye -> finite-T Fermi/Thomas-Fermi -> Lindhard/RPA static
!> screening of an electron gas with neutralizing background at r_s = 1.
!> Writes CSV evidence and output/screening-hierarchy.png.
program screening_demo
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortplot, only: figure, subplot, plot, title, suptitle, xlabel, ylabel, &
        xlim, ylim, legend, tight_layout, savefig, set_xscale, set_yscale
    use kin6d_species, only: atomic_units, electron, plasma_scales, plasma_scales_t
    use kin6d_screening_hierarchy
    use kin6d_result_io, only: run_t, result_file_t
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp), parameter :: blue(3) = [0.000_dp, 0.447_dp, 0.698_dp]
    real(dp), parameter :: verm(3) = [0.835_dp, 0.369_dp, 0.000_dp]
    real(dp), parameter :: green(3) = [0.000_dp, 0.620_dp, 0.451_dp]
    real(dp), parameter :: orange(3) = [0.902_dp, 0.624_dp, 0.000_dp]
    real(dp), parameter :: black(3) = [0.0_dp, 0.0_dp, 0.0_dp]
    integer, parameter :: nt = 41, nq = 161, nr = 400, nf = 4
    real(dp) :: density, fermi, theta(nt), ks_kd(nt), lind_kd(nt), debye(nt), tf(nt)
    real(dp) :: q(nq), eps(nq, nf), eps_a(nq, nf), r(nr), ratio(nr, nf), friedel(nr, nf)
    real(dp) :: theta_eps(nf) = [0.0_dp, 0.1_dp, 1.0_dp, 10.0_dp]
    real(dp) :: theta_fr(nf) = [0.0_dp, 0.02_dp, 0.05_dp, 0.1_dp]
    real(dp) :: tail(nf), ind(nf), tail_f(nf), ind_f(nf), ktf, kd10, ks1
    real(dp) :: route_gap, period, period_near, t0, t1, cost(4), worst_debye(2)
    real(dp) :: rl(1601), ratio_long(1601)
    type(screening_state_t) :: st, sts(nf)
    type(screening_value_t) :: a, b
    type(plasma_scales_t) :: scales
    integer :: i, j, unit, status, config_unit, ios
    logical :: ok
    type(run_t) :: run
    type(result_file_t) :: result
    character(256) :: msg
    real(dp) :: rs = 1.0_dp, qmax = 200.0_dp
    namelist /screening/ rs, qmax, theta_eps, theta_fr

    call run%start('screening_demo')
    if (run%argument_count() > 0) error stop 'usage: screening_demo [--config FILE]'
    ios = 0; msg = ''
    call run%open_config(config_unit)
    if (run%has_config()) read (config_unit, nml=screening, iostat=ios, iomsg=msg)
    call run%close_config(config_unit, 'screening', ios, msg)
    if (.not. (rs > 0 .and. qmax > 20 .and. all(theta_eps >= 0) .and. &
        all(theta_fr >= 0))) error stop 'screening: invalid configuration'
    call run%open_resolved(config_unit)
    write (config_unit, nml=screening)
    call run%close_resolved(config_unit)
    call execute_command_line('mkdir -p output', exitstat=status)
    density = 3.0_dp/(4*pi*rs**3)
    scales = plasma_scales(atomic_units(), electron(), density, 1.0_dp)
    fermi = scales%fermi

    ! 1. ks/kD along the degeneracy axis: fermi_gas EOS vs Lindhard q -> 0.
    call cpu_time(t0)
    route_gap = 0
    open (newunit=unit, file='output/screening-ks.csv', status='replace')
    write (unit, '(a)') '# Q7 static screening, electron gas r_s=1 (atomic units).'
    write (unit, '(a)') '# theta,gamma,ks_over_kd_fermi_gas,ks_over_kd_lindhard_q0,'// &
        'ks_lo,ks_hi,ktf0_over_kd,lindhard_tail,lindhard_indicator'
    do i = 1, nt
        theta(i) = 10.0_dp**(-2 + 4*(i - 1)/real(nt - 1, dp))
        call screening_setup(atomic_units(), electron(), density, theta(i)*fermi, st, ok)
        if (.not. ok) error stop 'screening setup failed'
        b = lindhard_maldague(st, 1e-6_dp)
        ks_kd(i) = st%ks_ratio/st%kd_ratio
        lind_kd(i) = sqrt(b%value)/st%kd_ratio
        debye(i) = 1; tf(i) = 1/st%kd_ratio
        route_gap = max(route_gap, abs(b%value/st%ks_ratio**2 - 1))
        write (unit, '(f0.8,8(",",es24.16e3))') theta(i), st%scales%gamma, ks_kd(i), &
            lind_kd(i), st%ks_ratio_bound%lo/st%kd_ratio, &
            st%ks_ratio_bound%hi/st%kd_ratio, tf(i), b%tail, b%indicator
    end do
    close (unit)
    call cpu_time(t1); cost(1) = t1 - t0

    ! 2. Static dielectric function, both finite-T routes.
    call cpu_time(t0)
    open (newunit=unit, file='output/screening-epsilon.csv', status='replace')
    write (unit, '(a)') '# q/kF,theta,eps_maldague,eps_momentum,|A-B|,tail,indicator'
    do j = 1, nf
        call screening_setup_normalized(theta_eps(j), 4/(pi*kf_ab()), sts(j), ok)
        do i = 1, nq
            q(i) = 10.0_dp**(-1.5_dp + 2.5_dp*(i - 1)/real(nq - 1, dp))
            a = lindhard_momentum(sts(j), q(i)/2)
            b = lindhard_response(sts(j), q(i)/2)
            eps(i, j) = sts(j)%kappa2*b%value/q(i)**2
            eps_a(i, j) = sts(j)%kappa2*a%value/q(i)**2
            write (unit, '(es16.8,6(",",es24.16e3))') q(i), theta_eps(j), &
                1 + eps(i, j), 1 + eps_a(i, j), abs(a%value - b%value), &
                a%tail + b%tail, a%indicator + b%indicator
        end do
    end do
    close (unit)
    call cpu_time(t1); cost(2) = t1 - t0

    ! 3. Screened point-charge potential phi/phi_bare (default route).
    call cpu_time(t0)
    do i = 1, nr
        r(i) = 0.05_dp*i
    end do
    do j = 1, nf
        call screened_ratio(sts(j), r, ratio(:, j), tail(j), ind(j), qmax)
    end do
    ktf = sqrt(sts(1)%kappa2)
    kd10 = ktf*sts(4)%kd_ratio
    ks1 = ktf*sts(3)%ks_ratio
    worst_debye = 0
    do i = 1, nr
        if (r(i) >= 1/(ktf*sts(4)%ks_ratio) .and. r(i) <= 4/(ktf*sts(4)%ks_ratio)) &
            worst_debye(1) = max(worst_debye(1), &
            abs(ratio(i, 4)/exp(-ktf*sts(4)%ks_ratio*r(i)) - 1))
    end do
    call cpu_time(t1); cost(3) = t1 - t0
    call cpu_time(t0)
    do j = 1, nf
        call screening_setup_normalized(theta_fr(j), 4/(pi*kf_ab()), st, ok)
        call screened_ratio(st, r, friedel(:, j), tail_f(j), ind_f(j))
        friedel(:, j) = friedel(:, j)*r**3
    end do
    period_near = friedel_period(r, friedel(:, 1))
    ! Asymptotic period from a long window 20 <= r kF <= 60 at T = 0.
    do i = 1, size(rl)
        rl(i) = 20.0_dp + 40.0_dp*(i - 1)/real(size(rl) - 1, dp)
    end do
    call screened_ratio(sts(1), rl, ratio_long, tail_f(1), ind_f(1))
    period = friedel_period(rl, ratio_long)
    call cpu_time(t1); cost(4) = t1 - t0
    open (newunit=unit, file='output/screening-potential.csv', status='replace')
    write (unit, '(a)') '# r*kF, phi/phi_bare at theta=0,0.1,1,10; r^3 phi/phi_bare at '// &
        'theta=0,0.02,0.05,0.1'
    write (unit, '(a,4es11.3,a,4es11.3)') '# q-tail bounds ', tail, ' indicators ', ind
    do i = 1, nr
        write (unit, '(f0.4,8(",",es24.16e3))') r(i), ratio(i, :), friedel(i, :)
    end do
    close (unit)

    open (newunit=unit, file='output/screening-summary.txt', status='replace')
    write (unit, '(a)') 'Q7 screening hierarchy, r_s = 1, ideal-gas RPA (no LFC)'
    write (unit, '(a,f8.5)') 'kappa2 = (kTF0/kF)^2 = ', sts(1)%kappa2
    write (unit, '(a,es10.3)') 'max |G(0)/(ks/kTF0)^2_fermi_gas - 1| over theta: ', &
        route_gap
    write (unit, '(a,f9.6,a,f9.6,a,f9.6)') 'Friedel period kF*lambda, 20<=r kF<=60: ', &
        period, '; 8<=r kF<=20: ', period_near, '  (pi = ', pi
    write (unit, '(a,es10.3)') 'theta=10 max rel. deviation from exp(-ks r), '// &
        '1<=ks r<=4: ', worst_debye(1)
    write (unit, '(a,4es10.3)') 'inversion q-tail bounds: ', tail
    write (unit, '(a,4es10.3)') 'inversion indicators:   ', ind
    write (unit, '(a,4f8.3)') 'cpu seconds ks-scan/epsilon/potential/friedel: ', cost
    close (unit)
    write (*, '(a,f9.6,a,es10.3,a,es10.3)') 'Friedel period ', period, &
        ' | G(0)-fermi_gas gap ', route_gap, ' | Debye-limit dev ', worst_debye(1)

    call result%create(run)
    call result%write('theta', theta, '1', 'T/T_F of the k_s scan', dims='theta')
    call result%write('ks_over_kd_fermi_gas', ks_kd, '1', &
        'finite-T Fermi screening k_s/k_D (kin6d_fermi_gas candidate)', dims='theta', &
        status='candidate inside fermi_gas certified enclosure')
    call result%write('ks_over_kd_lindhard', lind_kd, '1', &
        'sqrt(G(q->0)) of the Maldague route over k_D', dims='theta', &
        status='candidate; rigorous tail, quadrature indicator')
    call result%write('q', q, '1', 'q/k_F', dims='q')
    call result%write('theta_eps', theta_eps, '1', 'T/T_F of the epsilon curves', &
        dims='theta_eps')
    call result%write('epsilon_minus_one', eps, '1', 'static RPA eps(q)-1, default route', &
        dims='q,theta_eps', status='candidate; routes agree to roundoff')
    call result%write('epsilon_minus_one_momentum', eps_a, '1', &
        'static RPA eps(q)-1, momentum route', dims='q,theta_eps', status='candidate')
    call result%write('r', r, '1', 'r k_F', dims='r')
    call result%write('potential_ratio', ratio, '1', 'phi/phi_bare at theta_eps', &
        dims='r,theta_eps', status='candidate; q-tail bound and indicator in attrs')
    call result%write('theta_friedel', theta_fr, '1', 'T/T_F of Friedel curves', &
        dims='theta_friedel')
    call result%write('friedel', friedel, '1', '(r k_F)^3 phi/phi_bare', &
        dims='r,theta_friedel', status='candidate')
    call result%write('inversion_tail_bound', tail, '1', 'rigorous q-truncation bound', &
        dims='theta_eps', status='bound')
    call result%write('inversion_indicator', ind, '1', 'quadrature indicator', &
        dims='theta_eps', status='indicator, not a bound')
    call result%write('friedel_period', period, '1/k_F', 'Friedel period, 20<=r k_F<=60')
    call result%write('cpu_seconds', cost, 's', 'ks scan, epsilon, potential, Friedel', &
        dims='cost')
    call result%close()

    call figure(figsize=[13.2_dp, 8.8_dp])
    call subplot(2, 2, 1)
    call title('Static screening wavenumber, $r_s = 1$')
    call plot(theta, ks_kd, color=blue, label='finite-T Fermi (fermi_gas)')
    call plot(theta(1:nt:4), lind_kd(1:nt:4), color=verm, &
        linestyle='none', marker='o', markersize=5.0_dp, label='Lindhard, small-$q$ limit')
    call plot(theta, debye, color=black, linestyle='--', label='Debye')
    call plot(theta, tf, color=black, linestyle=':', &
        label='T=0 Thomas-Fermi')
    call set_xscale('log'); call set_yscale('log')
    call xlim(1e-2_dp, 1e2_dp); call ylim(0.1_dp, 15.0_dp)
    call xlabel('$\Theta = T/T_F$')
    call ylabel('$k_s / k_D$')
    call legend(loc='upper left')
    call subplot(2, 2, 2)
    call title('Static RPA dielectric function')
    call plot(q, eps(:, 2), color=blue, label='$\Theta$ = 0.1')
    call plot(q, eps(:, 3), color=verm, label='$\Theta$ = 1')
    call plot(q, eps(:, 4), color=green, label='$\Theta$ = 10')
    call plot(q(1:nq:10), eps_a(1:nq:10, 3), color=verm, &
        linestyle='none', marker='o', markersize=4.0_dp, &
        label='momentum route, $\Theta$ = 1')
    call plot(q, eps(:, 1), color=black, linestyle='--', &
        label='T=0 closed form')
    call set_xscale('log'); call set_yscale('log')
    call xlim(10.0_dp**(-1.5_dp), 10.0_dp); call ylim(10.0_dp**(-4.5_dp), 1e3_dp)
    call xlabel('$q / k_F$'); call ylabel('$\epsilon(q) - 1$')
    call legend(loc='lower left')
    call subplot(2, 2, 3)
    call title('Screened point charge')
    call plot(r, floor_abs(ratio(:, 1)), color=blue, label='Lindhard T=0')
    call plot(r, floor_abs(ratio(:, 3)), color=verm, label='Lindhard $\Theta$ = 1')
    call plot(r, floor_abs(ratio(:, 4)), color=green, label='Lindhard $\Theta$ = 10')
    call plot(r, floor_abs(exp(-ktf*r)), color=black, linestyle='--', label='Thomas-Fermi')
    call plot(r, floor_abs(exp(-ks1*r)), color=black, linestyle='-.', &
        label='Yukawa $k_s(\Theta=1)$')
    call plot(r, floor_abs(exp(-kd10*r)), color=black, linestyle=':', &
        label='Debye $\Theta$ = 10')
    call set_yscale('log')
    call xlim(0.0_dp, 20.0_dp); call ylim(1e-8_dp, 1.5_dp)
    call xlabel('$r k_F$'); call ylabel('$|\phi / \phi_{bare}|$')
    call legend(loc='lower left')
    call subplot(2, 2, 4)
    call title('Friedel oscillations and thermal damping')
    call plot(r(160:), clip(friedel(160:, 1)), color=blue, label='T=0')
    call plot(r(160:), clip(friedel(160:, 2)), color=verm, label='$\Theta$ = 0.02')
    call plot(r(160:), clip(friedel(160:, 3)), color=green, label='$\Theta$ = 0.05')
    call plot(r(160:), clip(friedel(160:, 4)), color=orange, label='$\Theta$ = 0.1')
    call xlim(8.0_dp, 20.0_dp); call ylim(-0.8_dp, 1.5_dp)
    call xlabel('$r k_F$'); call ylabel('$(r k_F)^3 \phi / \phi_{bare}$')
    call legend(loc='upper left')
    call suptitle('Ideal electron gas, static RPA without local-field correction; '// &
        'period ' // fixed(period) // ' $/k_F$ ($\pi$ = 3.1416)')
    call tight_layout()
    call savefig('output/screening-hierarchy.png')
contains
    !> |y| floored at the lower y limit of panel 3 (display only).
    pure function floor_abs(y) result(z)
        real(dp), intent(in) :: y(:)
        real(dp) :: z(size(y))
        z = max(abs(y), 1e-8_dp)
    end function
    !> Clip into the plotted band of panel 4 (display only).
    pure function clip(y) result(z)
        real(dp), intent(in) :: y(:)
        real(dp) :: z(size(y))
        z = max(-0.8_dp, min(1.5_dp, y))
    end function
    !> kF a_B = (9 pi/4)^(1/3)/r_s (textbook map, cross-checked in the test).
    real(dp) function kf_ab()
        kf_ab = (9*pi/4)**(1.0_dp/3.0_dp)/rs
    end function
    real(dp) function friedel_period(x, y) result(p)
        real(dp), intent(in) :: x(:), y(:)
        real(dp) :: first, last, z
        integer :: k, count
        count = 0; first = 0; last = 0
        do k = 1, size(x) - 1
            if (x(k) < 8.0_dp) cycle
            if (y(k)*y(k + 1) < 0) then
                z = x(k) - y(k)*(x(k + 1) - x(k))/(y(k + 1) - y(k))
                if (count == 0) first = z
                last = z; count = count + 1
            end if
        end do
        p = 2*(last - first)/(count - 1)
    end function
    function fixed(value) result(text)
        real(dp), intent(in) :: value
        character(:), allocatable :: text
        character(16) :: buffer
        write (buffer, '(f7.4)') value
        text = trim(adjustl(buffer))
    end function
end program screening_demo
