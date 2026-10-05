program steady_laboratory_demo
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use kin6d_fp_laboratory, only: fp_state_t, initialize_fp, solve_fp, &
        fp_certificate_radius, fp_certificate_accepts, fp_mass_interval, &
        fp_normalization_interval
    use kin6d_fp_leaves, only: fp_band, fp_reference, fp_mass_error_band
    use poisson_laboratory, only: coupled_steady_t, solve_coupled_steady, &
        coupled_sample, coupled_certificate_radius, &
        coupled_certificate_accepts
    use kin6d_poisson_leaves, only: poisson_reconstruct_interval, &
        coupled_distribution_interval
    use kin6d_leaves, only: coordinate, interval_band
    use fortnum_interval, only: interval_t, interval, interval_pi
    use fortplot, only: figure, subplot, plot, fill_between, title, xlabel, &
        ylabel, legend, savefig, tight_layout, xlim, ylim
    use fortplot_system_runtime, only: create_directory_runtime
    use kin6d_result_io, only: run_t, result_file_t
    implicit none
    type(fp_state_t) :: fp, coarse
    type(coupled_steady_t) :: coupled, krook, homogeneous
    type(interval_t) :: z, value, lower, upper, pi, potential, electric, energy
    integer :: cells = 64, nx = 96
    ! Coupled VPFP source drive, Krook limiting drive and collision rate.
    real(dp) :: coupled_drive = 1.875_dp, krook_drive = 0.375_dp, nu = 0.75_dp
    namelist /steady_laboratory/ cells, nx, coupled_drive, krook_drive, nu
    type(run_t) :: run
    type(result_file_t) :: result
    real(dp) :: radius, coarse_radius, residual, start, finish
    real(dp), allocatable :: lo(:), hi(:), clo(:), chi(:)
    real(dp), allocatable :: ref_x(:), ref_lo(:), ref_hi(:)
    real(dp), allocatable :: x(:), xf(:), flo(:), fhi(:)
    real(dp), allocatable :: elo(:), ehi(:), fline(:), eline(:)
    real(dp), allocatable :: limit_line(:), equilibrium_line(:)
    real(dp) :: fr, pr, er, e, p, source, defect, mass_lo, mass_hi
    real(dp) :: zlo, zhi, lower_mass, upper_mass
    real(dp) :: resolution(2), mass_lower(2), mass_upper(2), exact_mass(2)
    integer :: i, ledger, unit, ios
    logical :: ok
    character(256) :: msg
    character(64) :: label

    call run%start('steady_laboratory_demo')
    if (run%argument_count() > 0) error stop &
        'usage: steady_laboratory_demo [--config FILE]'
    ios = 0; msg = ''
    call run%open_config(unit)
    if (run%has_config()) read (unit, nml=steady_laboratory, iostat=ios, iomsg=msg)
    call run%close_config(unit, 'steady_laboratory', ios, msg)
    if (cells < 4 .or. mod(cells, 2) /= 0) error stop 'cells must be even and >= 4'
    if (nx < 1) error stop 'nx must be positive'
    if (.not. (nu > 0.0_dp)) error stop 'nu must be positive'
    call run%open_resolved(unit)
    write (unit, nml=steady_laboratory)
    call run%close_resolved(unit)
    allocate (lo(cells + 1), hi(cells + 1), clo(cells/2 + 1), chi(cells/2 + 1), &
        ref_x(2*cells), ref_lo(2*cells), ref_hi(2*cells), x(nx + 1), xf(2*nx), &
        flo(2*nx), fhi(2*nx), elo(2*nx), ehi(2*nx), fline(nx + 1), eline(nx + 1), &
        limit_line(nx + 1), equilibrium_line(nx + 1))
    call create_directory_runtime('output', ok)
    if (.not. ok) error stop 'cannot create output directory'
    pi = interval_pi()
    call initialize_fp(fp, cells)
    call initialize_fp(coarse, cells/2)
    call cpu_time(start)
    call solve_fp(fp, residual)
    call solve_fp(coarse, defect)
    call cpu_time(finish)
    print '(a,f10.6)', 'FP solve CPU seconds=', finish - start
    print '(a,es12.4,a,es12.4)', 'FP direct flux residual=', residual, &
        ' coarse residual=', defect
    call cpu_time(start)
    radius = fp_certificate_radius(fp)
    coarse_radius = fp_certificate_radius(coarse)
    if (.not. fp_certificate_accepts(fp, radius)) error stop 'FP certificate failed'
    if (.not. fp_certificate_accepts(coarse, coarse_radius)) &
        error stop 'coarse FP certificate failed'
    call cpu_time(finish)
    print '(a,f10.6)', 'FP certificate CPU seconds=', finish - start
    print '(a,i0,a,es12.4,a,i0,a,es12.4)', 'FP whole-cell sup bounds: ', cells, '=', &
        radius, ' ', cells/2, '=', coarse_radius
    call fp_normalization_interval(fp, zlo, zhi)
    z = interval(zlo, zhi)
    open (newunit=ledger, file='output/fp-candidate.csv', status='replace')
    write (ledger, '(a)') 'cells,v,f,lower,upper,v_binary64,f_binary64,' // &
        'lower_binary64,upper_binary64'
    do i = 1, cells + 1
        call fp_band(interval(fp%f(i)), interval(radius), lower, upper)
        lo(i) = lower%lo; hi(i) = upper%hi
        write (ledger, '(i0,4(",",es24.16),4(",",z16.16))') cells, &
            fp%v(i), fp%f(i), lo(i), hi(i), transfer(fp%v(i), 0_int64), &
            transfer(fp%f(i), 0_int64), transfer(lo(i), 0_int64), &
            transfer(hi(i), 0_int64)
    end do
    do i = 1, cells/2 + 1
        call fp_band(interval(coarse%f(i)), interval(coarse_radius), lower, upper)
        clo(i) = lower%lo; chi(i) = upper%hi
        write (ledger, '(i0,4(",",es24.16),4(",",z16.16))') cells/2, &
            coarse%v(i), coarse%f(i), clo(i), chi(i), &
            transfer(coarse%v(i), 0_int64), transfer(coarse%f(i), 0_int64), &
            transfer(clo(i), 0_int64), transfer(chi(i), 0_int64)
    end do
    close (ledger)
    do i = 1, cells
        ref_x(2*i - 1) = fp%v(i); ref_x(2*i) = fp%v(i + 1)
        call fp_reference(interval(fp%v(i), fp%v(i + 1)), interval(fp%theta), &
            interval(fp%mass), z, value)
        ref_lo(2*i - 1:2*i) = value%lo
        ref_hi(2*i - 1:2*i) = value%hi
    end do
    call fp_mass_interval(fp, mass_lo, mass_hi)
    call fp_mass_error_band(interval(mass_lo), interval(mass_hi), interval(radius), &
        lower, upper)
    lower_mass = lower%lo; upper_mass = upper%hi
    resolution = [real(cells/2, dp), real(cells, dp)]
    mass_lower(2) = lower_mass; mass_upper(2) = upper_mass
    call fp_mass_interval(coarse, mass_lo, mass_hi)
    call fp_mass_error_band(interval(mass_lo), interval(mass_hi), &
        interval(coarse_radius), lower, upper)
    mass_lower(1) = lower%lo; mass_upper(1) = upper%hi
    exact_mass = fp%mass
    print '(a,2es24.16)', 'FP continuous mass interval=', lower_mass, upper_mass
    call figure(figsize=[9.6_dp, 6.0_dp])
    call title('Bounded OU steady state | reviewed continuous bounds')
    call xlabel('v / thermal scale  (reflecting boundaries at +/-2)')
    call ylabel('normalized distribution f  (mass = 1)')
    call fill_between(coarse%v, clo, chi, color='moccasin', alpha=0.65_dp)
    call fill_between(fp%v, lo, hi, color='lightblue', alpha=0.75_dp)
    call plot(coarse%v, clo, linestyle=':', color=[0.8_dp, 0.4_dp, 0.0_dp])
    call plot(coarse%v, chi, linestyle=':', color=[0.8_dp, 0.4_dp, 0.0_dp])
    call plot(fp%v, lo, linestyle=':', color=[0.0_dp, 0.45_dp, 0.7_dp])
    call plot(fp%v, hi, linestyle=':', color=[0.0_dp, 0.45_dp, 0.7_dp])
    write (label, '(i0,a)') cells/2, ' cells + full band'
    call plot(coarse%v, coarse%f, label=trim(label), linestyle='--', &
        color=[0.8_dp, 0.4_dp, 0.0_dp])
    write (label, '(i0,a)') cells, ' cells + full band'
    call plot(fp%v, fp%f, label=trim(label), &
        color=[0.0_dp, 0.45_dp, 0.7_dp])
    call plot(ref_x, ref_lo, label='rigorous Gaussian cell edges', linestyle=':', &
        color=[0.0_dp, 0.0_dp, 0.0_dp])
    call plot(ref_x, ref_hi, linestyle=':', color=[0.0_dp, 0.0_dp, 0.0_dp])
    call legend(loc='upper right')
    call tight_layout()
    call savefig('output/fp-steady.png')

    ! The prescribed drive is a source density moment, not a supplied solution.
    call solve_coupled_steady(coupled_drive, nu, 1, coupled, ok, &
        fp_rate=1.0_dp, temperature=1.0_dp)
    if (.not. ok) error stop 'coupled FP steady solve failed'
    call solve_coupled_steady(krook_drive, nu, 1, krook, ok)
    if (.not. ok) error stop 'Krook limiting solve failed'
    call solve_coupled_steady(0.0_dp, nu, 1, homogeneous, ok)
    if (.not. ok) error stop 'homogeneous limiting solve failed'
    call coupled_certificate_radius(coupled, fr, pr, er, energy, ok)
    if (.not. ok) error stop 'coupled enclosure failed'
    if (.not. coupled_certificate_accepts(coupled, fr, pr, er)) &
        error stop 'coupled certificate rejected'
    print '(a,3es12.4)', 'Coupled pointwise bounds f,phi,E=', fr, pr, er
    print '(a,2es24.16)', 'Spatial mean electric energy interval=', energy%lo, energy%hi
    open (newunit=ledger, file='output/coupled-candidate.csv', status='replace')
    write (ledger, '(a)') 'x,f_at_v_zero,E,phi,source_at_v_zero,residual'
    do i = 1, nx + 1
        call coordinate(0.0_dp, 2.0_dp*pi%hi, real(i - 1, dp), real(nx, dp), x(i))
        call coupled_sample(coupled, x(i), 0.0_dp, fline(i), eline(i), p, &
            source, defect, ok)
        if (.not. ok) error stop 'coupled sample failed'
        write (ledger, '(es24.16,5(",",es24.16))') x(i), fline(i), eline(i), &
            p, source, defect
        call coupled_sample(krook, x(i), 0.0_dp, limit_line(i), e, p, source, &
            defect, ok)
        if (.not. ok) error stop 'Krook sample failed'
        call coupled_sample(homogeneous, x(i), 0.0_dp, equilibrium_line(i), e, p, &
            source, defect, ok)
        if (.not. ok) error stop 'homogeneous sample failed'
    end do
    close (ledger)
    do i = 1, nx
        xf(2*i - 1) = x(i); xf(2*i) = x(i + 1)
        call coupled_distribution_interval(interval(coupled%amplitude), &
            interval(real(coupled%k, dp)), interval(x(i), x(i + 1)), &
            interval(0.0_dp), value)
        call interval_band(value, interval(fr), lower, upper)
        flo(2*i - 1:2*i) = lower%lo; fhi(2*i - 1:2*i) = upper%hi
        call poisson_reconstruct_interval(interval(real(coupled%k, dp)), &
            interval(coupled%field%phi_cos), interval(coupled%field%phi_sin), &
            interval(coupled%field%e_cos), interval(coupled%field%e_sin), &
            interval(x(i), x(i + 1)), potential, electric)
        call interval_band(electric, interval(er), lower, upper)
        elo(2*i - 1:2*i) = lower%lo; ehi(2*i - 1:2*i) = upper%hi
    end do
    call figure(figsize=[9.6_dp, 7.2_dp])
    call subplot(2, 1, 1)
    call title('Source-balanced 1X1V Poisson-FP | whole-cell enclosures')
    call xlabel('x / slab scale  (period 2 pi)')
    call ylabel('f(x,v=0)  (velocity domain [-1,1])')
    call xlim(0.0_dp, 2.0_dp*pi%hi)
    call fill_between(xf, flo, fhi, color='lightblue', alpha=0.75_dp)
    call plot(xf, flo, linestyle=':', color=[0.0_dp, 0.45_dp, 0.7_dp])
    call plot(xf, fhi, linestyle=':', color=[0.0_dp, 0.45_dp, 0.7_dp])
    call plot(x, fline, label='FP candidate + full cell band', &
        color=[0.0_dp, 0.45_dp, 0.7_dp])
    call plot(x, limit_line, label='Krook limit (FP rate = 0)', linestyle='--', &
        color=[0.8_dp, 0.4_dp, 0.0_dp])
    call plot(x, equilibrium_line, label='zero drive / zero field', linestyle=':', &
        color=[0.0_dp, 0.0_dp, 0.0_dp])
    call legend(loc='upper right')
    call subplot(2, 1, 2)
    call title('Self-consistent electric field | reviewed continuous bounds')
    call xlabel('x / slab scale  (period 2 pi)')
    call ylabel('E / field scale')
    call xlim(0.0_dp, 2.0_dp*pi%hi)
    call fill_between(xf, elo, ehi, color='lightblue', alpha=0.75_dp)
    call plot(xf, elo, linestyle=':', color=[0.0_dp, 0.45_dp, 0.7_dp])
    call plot(xf, ehi, linestyle=':', color=[0.0_dp, 0.45_dp, 0.7_dp])
    call plot(x, eline, label='E from charge density + full cell band', &
        color=[0.0_dp, 0.45_dp, 0.7_dp])
    call legend(loc='upper right')
    call tight_layout(h_pad=3.0_dp)
    call savefig('output/poisson-fp-steady.png')
    call figure(figsize=[9.6_dp, 6.4_dp])
    call subplot(2, 1, 1)
    call title('Bounded OU mass | full continuous observable interval')
    call xlabel('velocity cells')
    call ylabel('integral f dv  (prescribed mass = 1)')
    call plot([resolution(1), resolution(1)], [mass_lower(1), mass_upper(1)], &
        label=trim(cell_label(cells/2)), color=[0.8_dp, 0.4_dp, 0.0_dp], marker='s')
    call plot([resolution(2), resolution(2)], [mass_lower(2), mass_upper(2)], &
        label=trim(cell_label(cells)), color=[0.0_dp, 0.45_dp, 0.7_dp], marker='o')
    call plot(resolution, exact_mass, label='exact prescribed mass', &
        linestyle='--', color=[0.0_dp, 0.0_dp, 0.0_dp])
    call legend(loc='upper right')
    call subplot(2, 1, 2)
    call title('Poisson electric energy | full continuous interval')
    call xlabel('observable 1: energy at normalized amplitude a = 1/2')
    call ylabel('spatial mean E squared / 2')
    call xlim(0.0_dp, 2.0_dp)
    call ylim(0.0_dp, 0.07_dp)
    call plot([1.0_dp, 1.0_dp], [energy%lo, energy%hi], &
        label='full interval (arithmetic width)', marker='o', &
        color=[0.0_dp, 0.45_dp, 0.7_dp])
    call legend(loc='upper right')
    call tight_layout(h_pad=3.0_dp)
    call savefig('output/steady-observables.png')
    open (newunit=ledger, file='output/steady-bounds.csv', status='replace')
    write (ledger, '(a)') 'case,quantity,lower,upper,proof_level'
    write (ledger, '(a,es24.16,a,es24.16,a)') 'bounded_ou,mass,', lower_mass, ',', &
        upper_mass, ',reviewed-continuous'
    write (ledger, '(a,es24.16,a,es24.16,a)') 'coupled_vpfp,mean_E_energy,', &
        energy%lo, ',', energy%hi, ',reviewed-continuous'
    close (ledger)
    open (newunit=ledger, file='output/steady-certificate.txt', status='replace')
    write (ledger, '(a)') 'schema_version=1'
    write (ledger, '(a)') 'proof_level=reviewed-continuous'
    write (ledger, '(a)') 'definitions=docs/steady-laboratory.md; gen/gen_fp.f90; ' // &
        'gen/gen_poisson.f90'
    write (ledger, '(a)') 'normalization=dimensionless; coupled q=m=epsilon0=1'
    write (ledger, '(a)') 'fp_domain=[-2,2]; no_flux; piecewise_linear; mass=theta=1'
    write (ledger, '(a)') 'fp_candidate_bytes=fp-candidate.csv IEEE binary64 columns'
    write (ledger, '(a,es24.16)') 'fp_64_sup_radius=', radius
    write (ledger, '(a,es24.16)') 'fp_32_sup_radius=', coarse_radius
    write (ledger, '(a,z16.16)') 'fp_64_radius_binary64=', transfer(radius, 0_int64)
    write (ledger, '(a,z16.16)') 'fp_32_radius_binary64=', transfer(coarse_radius, 0_int64)
    write (ledger, '(a)') 'coupled_domain=x periodic 2*pi,v[-1,1]; f=0 at v endpoints'
    write (ledger, '(a)') 'coupled_reconstruction=B(v)*(1+a*cos(x)); phi=p*cos(x)'
    write (ledger, '(a)') 'coupled_drive='//short(coupled_drive)//'; nu='//short(nu)// &
        '; beta=theta=k=1; fixed_charge=-1'
    write (ledger, '(a,z16.16)') 'coupled_amplitude_binary64=', &
        transfer(coupled%amplitude, 0_int64)
    write (ledger, '(a,z16.16)') 'coupled_phi_cos_binary64=', &
        transfer(coupled%field%phi_cos, 0_int64)
    write (ledger, '(a,z16.16)') 'coupled_e_sin_binary64=', &
        transfer(coupled%field%e_sin, 0_int64)
    write (ledger, '(a,3es24.16)') 'coupled_pointwise_radii_f_phi_E=', fr, pr, er
    write (ledger, '(a,3(z16.16,1x))') 'coupled_radii_binary64_f_phi_E=', &
        transfer(fr, 0_int64), transfer(pr, 0_int64), transfer(er, 0_int64)
    write (ledger, '(a,2(z16.16,1x))') 'mean_E_energy_binary64_lower_upper=', &
        transfer(energy%lo, 0_int64), transfer(energy%hi, 0_int64)
    write (ledger, '(a,2(z16.16,1x))') 'fp_mass_binary64_lower_upper=', &
        transfer(lower_mass, 0_int64), transfer(upper_mass, 0_int64)
    write (ledger, '(a)') 'formal_certificate=ComputedKrook.lean; ' // &
        'separate prescribed-field polynomial model'
    write (ledger, '(a)') 'exclusions=infinite_velocity_tails; general_nonlinear_uniqueness'
    write (ledger, '(a)') 'source_identity=source-base-commit.txt; source-manifest.sha256'
    close (ledger)
    call write_result()
    print '(a)', 'FP and coupled bounds use reviewed analysis and FortNum intervals.'
    print '(a)', 'The separate Lean Krook theorem has its own certificate check.'
contains
    function cell_label(n) result(text)
        integer, intent(in) :: n
        character(64) :: text
        write (text, '(i0,a)') n, ' cells: full interval'
    end function

    !> Shortest fixed decimal (up to 6 digits) that reads back exactly, else g0.
    function short(value) result(text)
        real(dp), intent(in) :: value
        character(len=:), allocatable :: text
        character(64) :: buffer
        real(dp) :: back
        integer :: k, ios
        write (buffer, '(f0.6)') value
        k = len_trim(buffer)
        do while (k > 1 .and. buffer(k:k) == '0')
            k = k - 1
        end do
        read (buffer(:k), *, iostat=ios) back
        if (ios == 0 .and. back == value) then
            text = buffer(:k)
        else
            write (buffer, '(g0)') value
            text = trim(buffer)
        end if
    end function

    subroutine write_result()
        character(len=*), parameter :: bound = 'certified-continuous '// &
            '(reviewed-continuous analysis, FortNum intervals); no tails/uniqueness claim'
        call result%create(run)
        call result%write('fp_v', fp%v, 'thermal scale', 'FP velocity nodes', dims='v')
        call result%write('fp_f', fp%f, '1', 'bounded OU/FP candidate', dims='v', &
            status='candidate center')
        call result%write('fp_lower', lo, '1', 'FP whole-cell lower band', dims='v', &
            status=bound)
        call result%write('fp_upper', hi, '1', 'FP whole-cell upper band', dims='v', &
            status=bound)
        call result%write('fp_coarse_v', coarse%v, 'thermal scale', &
            'coarse FP velocity nodes', dims='v_coarse')
        call result%write('fp_coarse_f', coarse%f, '1', 'coarse FP candidate', &
            dims='v_coarse', status='candidate center')
        call result%write('fp_coarse_lower', clo, '1', 'coarse FP lower band', &
            dims='v_coarse', status=bound)
        call result%write('fp_coarse_upper', chi, '1', 'coarse FP upper band', &
            dims='v_coarse', status=bound)
        call result%write('fp_sup_radius', radius, '1', 'FP whole-cell sup radius', &
            status=bound)
        call result%write('fp_coarse_sup_radius', coarse_radius, '1', &
            'coarse FP whole-cell sup radius', status=bound)
        call result%write('fp_mass_bounds', [lower_mass, upper_mass], '1', &
            'FP mass interval [lower, upper]', dims='bound', status=bound)
        call result%write('x', x, 'slab scale', 'coupled sample points', dims='x')
        call result%write('coupled_f_v0', fline, '1', 'coupled f(x, v=0) candidate', &
            dims='x', status='candidate center')
        call result%write('coupled_E', eline, 'field scale', 'coupled E(x) candidate', &
            dims='x', status='candidate center')
        call result%write('coupled_radii', [fr, pr, er], '1', &
            'pointwise radii of f, phi, E', dims='field', status=bound)
        call result%write('mean_E_energy_bounds', [energy%lo, energy%hi], '1', &
            'spatial mean electric energy interval', dims='bound', status=bound)
        call result%write('krook_f_v0', limit_line, '1', 'Krook limit f(x, v=0)', &
            dims='x', status='candidate')
        call result%close()
    end subroutine
end program
