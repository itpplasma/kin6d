program relativistic_fermi_demo
    ! Q8 evidence: certified relativistic ideal electron pressure versus density
    ! at T=1e7 K (white-dwarf-like, Ye=1/2), with the classical, cold
    ! nonrelativistic, cold ultrarelativistic and exact cold (Chandrasekhar)
    ! curves and the nonrelativistic finite-T EOS. Electrons only: no ions,
    ! Coulomb/correlation corrections or constants uncertainty. CSV and
    ! preliminary PNG figures go to output/ (not committed).
    use, intrinsic :: iso_fortran_env, only: dp => real64, qp => real128
    use kin6d_fermi_gas, only: fermi_workspace_t, fermi_state_t, fermi_si_state_t, &
        initialize_fermi_workspace, solve_fermi_gas_si, occupation_fermi_gas
    use kin6d_fermi_gas_leaves, only: fermi_si_parameters
    use kin6d_relativistic_fermi_leaves, only: fermi_rel_limits
    use fortplot, only: figure, plot, title, xlabel, ylabel, ylim, legend, savefig, &
        xscale, yscale
    use fortplot_system_runtime, only: create_directory_runtime
    implicit none
    ! CODATA 2018 values, supplied explicitly (exact SI where defined).
    real(dp), parameter :: mass = 9.1093837015e-31_dp, hbar = 1.054571817e-34_dp
    real(dp), parameter :: speed = 299792458.0_dp, kb = 1.380649e-23_dp
    real(dp), parameter :: charge = 1.602176634e-19_dp, eps0 = 8.8541878128e-12_dp
    real(dp), parameter :: atomic_mass = 1.66053906660e-27_dp, mu_e = 2.0_dp
    real(dp), parameter :: temperature = 1.0e7_dp, pi = acos(-1.0_dp)
    integer, parameter :: points = 25, shown(3) = [7, 13, 25]
    real(dp), parameter :: blue(3) = [0.0_dp, 0.447_dp, 0.698_dp]
    real(dp), parameter :: orange(3) = [0.835_dp, 0.369_dp, 0.0_dp]
    real(dp), parameter :: green(3) = [0.0_dp, 0.62_dp, 0.451_dp]
    real(dp), parameter :: gray(3) = [0.45_dp, 0.45_dp, 0.45_dp]
    real(dp), parameter :: palette(3, 3) = reshape([blue, orange, green], [3, 3])
    type(fermi_workspace_t) :: work
    type(fermi_si_state_t) :: rel, nonrel
    real(dp) :: x(points), n(points), rho(points), p(points), plo(points), phi(points)
    real(dp) :: pnr(points), pcl(points), pnr0(points), pur0(points), pch(points)
    real(dp) :: scale(points), theta(points), pf, ef, limits(3), t0, t1, u(161), f(161, 3), lo, hi
    logical :: ok, certified(points), nonrel_ok(points)
    integer :: i, j, unit
    character(24) :: label

    call create_directory_runtime('output', ok)
    call initialize_fermi_workspace(work, 65536, ok)
    if (.not. ok) error stop 'workspace'
    open (newunit=unit, file='output/relativistic-fermi-pressure.csv', status='replace')
    write (unit, '(a)') '# ideal electron gas, T=1e7 K, Ye=1/2, CODATA 2018 constants; '// &
        'P certified-continuous enclosure (electrons only); limits are references'
    write (unit, '(a)') 'xF,n_m3,rho_kg_m3,theta,P_Pa,P_lo,P_hi,certified,'// &
        'P_nonrel_T(0=certificate failed),P_classical,P_nonrel_T0,P_ultrarel_T0,P_chandrasekhar_T0,cells,seconds'
    do i = 1, points
        x(i) = 0.01_dp*1000.0_dp**(real(i - 1, dp)/real(points - 1, dp))
        n(i) = (mass*speed*x(i)/hbar)**3/(3.0_dp*pi**2)
        rho(i) = n(i)*mu_e*atomic_mass
        call cpu_time(t0)
        call solve_fermi_gas_si(work, n(i), temperature, mass, hbar, kb, charge, eps0, &
                                rel, ok, speed=speed)
        call cpu_time(t1)
        certified(i) = ok .and. rel%certified
        p(i) = rel%pressure; plo(i) = rel%pressure_bound%lo; phi(i) = rel%pressure_bound%hi
        call solve_fermi_gas_si(work, n(i), temperature, mass, hbar, kb, charge, eps0, &
                                nonrel, ok)
        ! The nonrelativistic certificate uses a uniform mesh and fails closed
        ! for very cold states; those points are omitted from its curve.
        nonrel_ok(i) = ok .and. nonrel%certified
        pnr(i) = nonrel%pressure
        call fermi_si_parameters(n(i), temperature, mass, hbar, kb, pf, ef, theta(i))
        call fermi_rel_limits(theta(i), x(i), limits(1), limits(2), limits(3))
        pnr0(i) = n(i)*ef*limits(1); pur0(i) = n(i)*ef*limits(2)
        pcl(i) = n(i)*ef*limits(3)
        pch(i) = chandrasekhar(x(i))
        write (unit, '(es12.5,6(",",es16.8),",",l1,5(",",es16.8),",",i6,",",f8.3)') &
            x(i), n(i), rho(i), theta(i), p(i), plo(i), phi(i), certified(i), pnr(i), &
            pcl(i), pnr0(i), pur0(i), pch(i), rel%normalized%cells, t1 - t0
        print '(a,es9.2,a,es9.2,a,l2,a,es12.5,a,f8.5,a,f8.5,a,f6.2)', 'xF=', x(i), &
            ' theta=', theta(i), ' certified', certified(i), ' P=', p(i), &
            ' P/Pch0=', p(i)/pch(i), ' P/Pnr=', merge(p(i)/pnr(i), 0.0_dp, &
            nonrel_ok(i)), ' s=', t1 - t0
    end do
    close (unit)
    ! Occupations versus p/pF at three densities (normalized certified state).
    do j = 1, 3
        call solve_fermi_gas_si(work, n(shown(j)), temperature, mass, hbar, kb, &
                                charge, eps0, rel, ok, speed=speed)
        do i = 1, 161
            u(i) = 2.0_dp*real(i - 1, dp)/160.0_dp
            call occupation_fermi_gas(rel%normalized, u(i), f(i, j), lo, hi, ok)
        end do
    end do
    ! Single-axes figures: fortplot subplot axis scales currently leak
    ! between panels of one figure.
    call figure(figsize=[7.5_dp, 5.5_dp])
    ! Pressure per rest-energy density P/(n m c^2) keeps the 12 decades of
    ! P(rho) readable: classical kT/mc^2, cold NR x^2/5, cold UR x/4.
    scale = n*mass*speed**2
    call plot(x, pch/scale, label='T=0 exact (Chandrasekhar)', &
              color=[0.0_dp, 0.0_dp, 0.0_dp])
    call plot(x, pnr0/scale, label='T=0 nonrelativistic', color=gray, linestyle='--')
    call plot(x, pur0/scale, label='T=0 ultrarelativistic', color=gray, linestyle=':')
    call plot(x, pcl/scale, label='classical nkT', color=green, linestyle='-.')
    call plot(pack(x, nonrel_ok), pack(pnr/scale, nonrel_ok), &
              label='nonrelativistic, T=1e7 K', color=orange)
    call plot(x, p/scale, label='relativistic, T=1e7 K (certified)', color=blue, &
              linestyle='None', marker='o')
    call xscale('log'); call yscale('log')
    call title('Ideal electrons, T=1e7 K, rho=1.95e9 x_F^3 kg/m^3')
    call xlabel('x_F = p_F/(m c)'); call ylabel('P/(n m c^2)')
    call ylim(1.0e-3_dp, 10.0_dp); call legend(loc='upper left')
    call savefig('output/relativistic-fermi-pressure.png')
    call figure(figsize=[7.5_dp, 5.5_dp])
    call plot(x, p/pch, label='relativistic / T=0 exact', color=blue, marker='o')
    call plot(pack(x, nonrel_ok), pack(p/merge(pnr, 1.0_dp, nonrel_ok), nonrel_ok), &
              label='relativistic / nonrelativistic (same T)', color=orange, marker='s')
    call plot(x, p/pcl, label='relativistic / classical nkT', color=green, marker='^')
    call xscale('log')
    call title('Pressure ratios at T=1e7 K'); call xlabel('x_F = p_F/(m c)')
    call ylabel('ratio'); call ylim(0.0_dp, 3.0_dp); call legend(loc='upper right')
    call savefig('output/relativistic-fermi-ratios.png')
    call figure(figsize=[7.5_dp, 5.5_dp])
    do j = 1, 3
        write (label, '(a,es8.1)') 'x_F=', x(shown(j))
        call plot(u, f(:, j), label=trim(label), color=palette(:, j))
    end do
    call title('Occupation at T=1e7 K'); call xlabel('p/p_F'); call ylabel('f')
    call ylim(0.0_dp, 1.05_dp); call legend(loc='upper right')
    call savefig('output/relativistic-fermi-occupation.png')

contains

    !> Exact cold pressure (m^4 c^5/(24 pi^2 hbar^3))[x(2x^2-3)s+3 asinh x],
    !> evaluated in binary128 to absorb the small-x cancellation.
    real(dp) function chandrasekhar(xf) result(pressure)
        real(dp), intent(in) :: xf
        real(qp) :: xq
        xq = real(xf, qp)
        pressure = real(real(mass, qp)**4*real(speed, qp)**5/(24.0_qp* &
            acos(-1.0_qp)**2*real(hbar, qp)**3)*(xq*(2.0_qp*xq*xq - 3.0_qp)* &
            sqrt(1.0_qp + xq*xq) + 3.0_qp*asinh(xq)), dp)
    end function chandrasekhar
end program relativistic_fermi_demo
