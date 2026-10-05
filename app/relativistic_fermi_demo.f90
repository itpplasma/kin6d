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
    use kin6d_species, only: electron_mass_si, hbar_si, elementary_charge_si, coulomb_si
    use fortplot, only: figure, plot, title, xlabel, ylabel, ylim, legend, savefig, &
        xscale, yscale
    use fortplot_system_runtime, only: create_directory_runtime
    use kin6d_result_io, only: run_t, result_file_t
    implicit none
    ! Electron constants from the shared species record (CODATA 2018); the
    ! remaining exact SI constants and the atomic mass unit are declared here.
    real(dp), parameter :: mass = electron_mass_si, hbar = hbar_si
    real(dp), parameter :: speed = 299792458.0_dp, kb = 1.380649e-23_dp
    real(dp), parameter :: charge = elementary_charge_si
    real(dp), parameter :: eps0 = 1.0_dp/(4.0_dp*acos(-1.0_dp)*coulomb_si)
    real(dp), parameter :: atomic_mass = 1.66053906660e-27_dp, pi = acos(-1.0_dp)
    ! Namelist &relativistic_fermi: nucleon masses per electron mu_e,
    ! temperature [K], x_F = p_F/(m c) log grid x_min*x_ratio**[0,1], three
    ! occupation densities (grid indices), occupation grid on [0, pmax] p_F.
    real(dp) :: mu_e = 2.0_dp, temperature = 1.0e7_dp
    real(dp) :: x_min = 0.01_dp, x_ratio = 1000.0_dp, occupation_pmax = 2.0_dp
    integer :: points = 25, shown(3) = [7, 13, 25], workspace_cells = 65536
    integer :: occupation_points = 161
    namelist /relativistic_fermi/ mu_e, temperature, x_min, x_ratio, points, &
        shown, workspace_cells, occupation_points, occupation_pmax
    type(run_t) :: run
    type(result_file_t) :: result
    integer :: ios
    character(256) :: msg
    real(dp), parameter :: blue(3) = [0.0_dp, 0.447_dp, 0.698_dp]
    real(dp), parameter :: orange(3) = [0.835_dp, 0.369_dp, 0.0_dp]
    real(dp), parameter :: green(3) = [0.0_dp, 0.62_dp, 0.451_dp]
    real(dp), parameter :: gray(3) = [0.45_dp, 0.45_dp, 0.45_dp]
    real(dp), parameter :: palette(3, 3) = reshape([blue, orange, green], [3, 3])
    type(fermi_workspace_t) :: work
    type(fermi_si_state_t) :: rel, nonrel
    real(dp), allocatable :: x(:), n(:), rho(:), p(:), plo(:), phi(:)
    real(dp), allocatable :: pnr(:), pcl(:), pnr0(:), pur0(:), pch(:)
    real(dp), allocatable :: scale(:), theta(:), u(:), f(:, :), cells(:), seconds(:)
    real(dp) :: pf, ef, limits(3), t0, t1, lo, hi
    logical, allocatable :: certified(:), nonrel_ok(:)
    logical :: ok
    integer :: i, j, unit
    character(24) :: label
    character(64) :: tlabel

    call run%start('relativistic_fermi_demo')
    if (run%argument_count() > 0) error stop 'usage: relativistic_fermi_demo [--config FILE]'
    ios = 0; msg = ''
    call run%open_config(unit)
    if (run%has_config()) read (unit, nml=relativistic_fermi, iostat=ios, iomsg=msg)
    call run%close_config(unit, 'relativistic_fermi', ios, msg)
    if (points < 2) error stop 'points must be at least 2'
    if (any(shown < 1) .or. any(shown > points)) error stop 'shown indices outside 1..points'
    if (.not. (mu_e > 0.0_dp .and. temperature > 0.0_dp .and. x_min > 0.0_dp .and. &
        x_ratio > 1.0_dp .and. occupation_pmax > 0.0_dp)) error stop 'invalid physical grid'
    if (workspace_cells < 1 .or. occupation_points < 2) error stop 'invalid resolution'
    call run%open_resolved(unit)
    write (unit, nml=relativistic_fermi)
    call run%close_resolved(unit)
    allocate (x(points), n(points), rho(points), p(points), plo(points), phi(points), &
        pnr(points), pcl(points), pnr0(points), pur0(points), pch(points), &
        scale(points), theta(points), cells(points), seconds(points), &
        certified(points), nonrel_ok(points), u(occupation_points), &
        f(occupation_points, 3))
    write (tlabel, '(a,es9.2,a)') 'T=', temperature, ' K'

    call create_directory_runtime('output', ok)
    call initialize_fermi_workspace(work, workspace_cells, ok)
    if (.not. ok) error stop 'workspace'
    open (newunit=unit, file='output/relativistic-fermi-pressure.csv', status='replace')
    write (unit, '(a,es9.2,a,es9.2,a)') '# ideal electron gas, T=', temperature, &
        ' K, mu_e=', mu_e, ', CODATA 2018 constants; '// &
        'P certified-continuous enclosure (electrons only); limits are references'
    write (unit, '(a)') 'xF,n_m3,rho_kg_m3,theta,P_Pa,P_lo,P_hi,certified,'// &
        'P_nonrel_T(0=certificate failed),P_classical,P_nonrel_T0,P_ultrarel_T0,P_chandrasekhar_T0,cells,seconds'
    do i = 1, points
        x(i) = x_min*x_ratio**(real(i - 1, dp)/real(points - 1, dp))
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
        cells(i) = real(rel%normalized%cells, dp); seconds(i) = t1 - t0
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
        do i = 1, occupation_points
            u(i) = occupation_pmax*real(i - 1, dp)/real(occupation_points - 1, dp)
            call occupation_fermi_gas(rel%normalized, u(i), f(i, j), lo, hi, ok)
        end do
    end do
    call write_result()
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
              label='nonrelativistic, '//trim(tlabel), color=orange)
    call plot(x, p/scale, label='relativistic, '//trim(tlabel)//' (certified)', color=blue, &
              linestyle='None', marker='o')
    call xscale('log'); call yscale('log')
    call title('Ideal electrons, '//trim(tlabel))
    call xlabel('x_F = p_F/(m c)'); call ylabel('P/(n m c^2)')
    call ylim(1.0e-3_dp, 10.0_dp); call legend(loc='upper left')
    call savefig('output/relativistic-fermi-pressure.png')
    call figure(figsize=[7.5_dp, 5.5_dp])
    call plot(x, p/pch, label='relativistic / T=0 exact', color=blue, marker='o')
    call plot(pack(x, nonrel_ok), pack(p/merge(pnr, 1.0_dp, nonrel_ok), nonrel_ok), &
              label='relativistic / nonrelativistic (same T)', color=orange, marker='s')
    call plot(x, p/pcl, label='relativistic / classical nkT', color=green, marker='^')
    call xscale('log')
    call title('Pressure ratios at '//trim(tlabel)); call xlabel('x_F = p_F/(m c)')
    call ylabel('ratio'); call ylim(0.0_dp, 3.0_dp); call legend(loc='upper right')
    call savefig('output/relativistic-fermi-ratios.png')
    call figure(figsize=[7.5_dp, 5.5_dp])
    do j = 1, 3
        write (label, '(a,es8.1)') 'x_F=', x(shown(j))
        call plot(u, f(:, j), label=trim(label), color=palette(:, j))
    end do
    call title('Occupation at '//trim(tlabel)); call xlabel('p/p_F'); call ylabel('f')
    call ylim(0.0_dp, 1.05_dp); call legend(loc='upper right')
    call savefig('output/relativistic-fermi-occupation.png')

contains

    subroutine write_result()
        character(len=*), parameter :: bound = 'certified-continuous enclosure '// &
            '(electrons only; constants uncertainty unavailable)'
        call result%create(run)
        call result%write('x_fermi', x, '1', 'relativity parameter p_F/(m c)', dims='point')
        call result%write('number_density', n, 'm-3', 'electron number density', dims='point')
        call result%write('mass_density', rho, 'kg m-3', 'mass density mu_e m_u n', dims='point')
        call result%write('theta', theta, '1', 'degeneracy kT/E_F (nonrelativistic E_F)', &
            dims='point')
        call result%write('pressure', p, 'Pa', 'relativistic ideal electron pressure', &
            dims='point', status='candidate center')
        call result%write('pressure_lower', plo, 'Pa', 'pressure lower bound', &
            dims='point', status=bound//' where certified=1')
        call result%write('pressure_upper', phi, 'Pa', 'pressure upper bound', &
            dims='point', status=bound//' where certified=1')
        call result%write('certified', merge(1.0_dp, 0.0_dp, certified), '1', &
            '1 if the relativistic pressure certificate succeeded', dims='point')
        call result%write('pressure_nonrelativistic', merge(pnr, 0.0_dp, nonrel_ok), 'Pa', &
            'nonrelativistic finite-T pressure (0 = certificate failed)', dims='point', &
            status='candidate center; certified enclosure exists where nonzero')
        call result%write('pressure_classical', pcl, 'Pa', 'classical n k T', dims='point', &
            status='reference limit')
        call result%write('pressure_nonrelativistic_cold', pnr0, 'Pa', &
            'cold nonrelativistic pressure', dims='point', status='reference limit')
        call result%write('pressure_ultrarelativistic_cold', pur0, 'Pa', &
            'cold ultrarelativistic pressure', dims='point', status='reference limit')
        call result%write('pressure_chandrasekhar', pch, 'Pa', &
            'exact cold (Chandrasekhar) pressure, binary128', dims='point', &
            status='reference limit')
        call result%write('cells', cells, '1', 'certificate mesh cells', dims='point')
        call result%write('solve_cpu', seconds, 's', 'relativistic solve CPU time', &
            dims='point')
        call result%write('occupation_momentum', u, '1', 'p/p_F', dims='momentum')
        call result%write('occupation', f, '1', 'occupation at the shown densities', &
            dims='momentum,shown', status='candidate center')
        call result%write('shown_x_fermi', x(shown), '1', 'x_F of the shown densities', &
            dims='shown')
        call result%close()
    end subroutine

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
