program nested_oracle_export
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_nested_oracles, only: nested_oracle_t, NESTED_INTEGER
    use kin6d_nested_oracle_leaves, only: nested_boundary_radius, &
        nested_boundary_harmonic, nested_boundary_phase, nested_setup, &
        nested_integer_averages
    use kin6d_result_io, only: run_t, result_file_t
    use fortnum_interval, only: interval_t
    implicit none
    type(nested_oracle_t) :: oracle
    type(run_t) :: run
    type(result_file_t) :: result
    integer :: family, mmax, nmax, angles, quadrature, turns, unit, ios, k, l, m, n, status, bad
    integer :: index, jj
    real(dp) :: eps, tau, delta, saxis, lam, edge_pressure, requested_radius, pi
    real(dp) :: theta, phi, x(3), r, cs, sn, rec(2), error, term(4), iota(3), fine_iota
    real(dp) :: coarse(4), fine(4), radius
    real(dp) :: limit, extent, pa1, pa2, volume, pmean, bmean, beta, transform, flux
    real(dp), allocatable :: points(:, :), fields(:, :), lower(:, :), upper(:, :), rz(:, :)
    real(dp), allocatable :: coefficients(:, :, :)
    type(interval_t) :: enclosure(20)
    character(256) :: message
    namelist /nested_oracle/ family, eps, tau, delta, saxis, lam, edge_pressure, &
        mmax, nmax, angles, quadrature, turns, requested_radius
    family = oracle%family; eps = oracle%eps; tau = oracle%tau; delta = oracle%delta
    saxis = oracle%saxis; lam = oracle%lam; edge_pressure = oracle%edge_pressure
    mmax = 12; nmax = 12; angles = 64; quadrature = 16; turns = 128
    requested_radius = 0.01_dp
    pi = acos(-1.0_dp)
    call run%start('nested_oracle_export')
    call run%open_config(unit)
    ios = 0; message = ''
    if (run%has_config()) read (unit, nml=nested_oracle, iostat=ios, iomsg=message)
    call run%close_config(unit, 'nested_oracle', ios, message)
    if (run%argument_count() /= 0) error stop 'usage: nested_oracle_export [--config FILE]'
    oracle%family = family; oracle%eps = eps; oracle%tau = tau; oracle%delta = delta
    oracle%saxis = saxis; oracle%lam = lam; oracle%edge_pressure = edge_pressure
    if (.not. oracle%valid()) error stop 'nested_oracle_export: invalid domain'
    if (mmax < 1 .or. nmax < 1 .or. angles < 4*max(mmax, nmax) .or. quadrature < 2 &
        .or. turns < 1 .or. requested_radius <= 0) error stop 'invalid export resolution'
    call run%open_resolved(unit)
    write (unit, nml=nested_oracle)
    call run%close_resolved(unit)
    allocate (points(3, angles*angles), fields(20, angles*angles), &
              lower(20, angles*angles), upper(20, angles*angles), rz(2, angles*angles))
    allocate (coefficients(4, -nmax:nmax, 0:mmax))
    bad = 0
    !$omp parallel do collapse(2) schedule(static) &
    !$omp private(theta,phi,index,x,status,enclosure,jj) reduction(+:bad)
    do l = 0, angles - 1
        do k = 0, angles - 1
            theta = 2*pi*k/angles; phi = 2*pi*l/angles
            index = 1 + k + angles*l
            call oracle%boundary(theta, phi, x, status)
            bad = bad + status
            points(:, index) = x
            call oracle%point(x, fields(:, index), status)
            bad = bad + status
            call oracle%point_enclosure(x, enclosure)
            do jj = 1, 20
                lower(jj, index) = enclosure(jj)%lo
                upper(jj, index) = enclosure(jj)%hi
            end do
            call nested_boundary_radius(x(1), x(2), rz(1, index))
            rz(2, index) = x(3)
        end do
    end do
    !$omp end parallel do
    if (bad /= 0) error stop 'nested_oracle_export: rejected boundary point'
    coefficients = 0
    !$omp parallel do collapse(2) schedule(static) private(theta,phi,index,term,k,l)
    do m = 0, mmax
        do n = -nmax, nmax
            if (m == 0 .and. n < 0) cycle
            do l = 0, angles - 1
                do k = 0, angles - 1
                    theta = 2*pi*k/angles; phi = 2*pi*l/angles
                    index = 1 + k + angles*l
                    call nested_boundary_harmonic(rz(1, index), rz(2, index), &
                                                  theta, phi, real(m, dp), real(n, dp), term(1), term(2), term(3), term(4))
                    coefficients(:, n, m) = coefficients(:, n, m) + term
                end do
            end do
            coefficients(:, n, m) = coefficients(:, n, m)/real(angles*angles, dp)
            if (m /= 0 .or. n /= 0) coefficients(:, n, m) = 2*coefficients(:, n, m)
        end do
    end do
    !$omp end parallel do
    if (tau == 0) then
        coefficients(2, :, :) = 0
        coefficients(3, :, :) = 0
    end if
    ! Independent staggered grid checks the reconstructed physical boundary.
    ! This is a sampled truncation indicator, not a continuous surface bound.
    error = 0; bad = 0
    !$omp parallel do collapse(2) schedule(static) &
    !$omp private(theta,phi,x,status,r,rec,cs,sn,m,n) &
    !$omp reduction(max:error) reduction(+:bad)
    do l = 0, 2*angles - 1
        do k = 0, 2*angles - 1
            theta = pi*(k + 0.37_dp)/angles; phi = pi*(l + 0.23_dp)/angles
            call oracle%boundary(theta, phi, x, status)
            bad = bad + status
            call nested_boundary_radius(x(1), x(2), r)
            rec = 0
            do m = 0, mmax
                do n = -nmax, nmax
                    call nested_boundary_phase(theta, phi, real(m, dp), real(n, dp), cs, sn)
                    rec = rec + [coefficients(1, n, m)*cs + coefficients(2, n, m)*sn, &
                                 coefficients(3, n, m)*cs + coefficients(4, n, m)*sn]
                end do
            end do
            error = max(error, norm2(rec - [r, x(3)]))
        end do
    end do
    !$omp end parallel do
    if (bad /= 0) error stop 'nested_oracle_export: verification point rejected'
    radius = sqrt(delta)
    if (family /= NESTED_INTEGER) radius = sqrt(2*delta)/lam
    if (error > requested_radius*radius) &
        error stop 'boundary truncation exceeds requested radius'
    call oracle%averages(quadrature, coarse)
    call oracle%averages(2*quadrature, fine)
    call oracle%rotation(0.0_dp, turns, 1e-7_dp, iota(1), status)
    if (status /= 0) error stop 'axis rotation integration failed'
    call oracle%rotation(delta, turns, 1e-7_dp, iota(2), status)
    if (status /= 0) error stop 'edge rotation integration failed'
    call oracle%rotation(delta, 2*turns, 3e-8_dp, fine_iota, status)
    if (status /= 0) error stop 'rotation refinement failed'
    iota(3) = abs(fine_iota - iota(2))
    call execute_command_line('mkdir -p output', exitstat=status)
    if (status /= 0) error stop 'cannot create output directory'
    call result%create(run)
    call result%attribute('physical_model', 'smooth nested ideal MHD; mu0=1')
    call result%attribute('source', 'arXiv:2609.26742v1; arXiv:2610.07304v1')
    call result%attribute('symbolic_proof', 'eight exact FortSym divergence/tension lemmas')
    call result%attribute('coordinates', &
        'Cartesian vectors; theta clockwise; geometric phi positive')
    call result%attribute('field_rows', &
        'psi,p,Bx,By,Bz,Jx,Jy,Jz,px,py,pz,grad B column-major')
    call result%attribute('boundary_input_scope', &
                          'geometry only: profiles and flux constraints must be supplied separately')
    call result%attribute('model_input_stochastic_error', &
        'not applicable to exact named reference')
    call result%write('boundary_points', points, 'm', 'stored Cartesian evaluation points', &
                      'xyz,point', status='candidate geometry; map evaluation bound unavailable')
    call result%write('fields', fields, 'normalized', 'oracle values at stored points', &
                      'value,point', status='exact closed-form reference; floating evaluation')
    call result%write('fields_lower', lower, 'normalized', &
        'outward evaluation lower endpoints', &
                      'value,point', status='point evaluation enclosure at stored Cartesian points')
    call result%write('fields_upper', upper, 'normalized', &
        'outward evaluation upper endpoints', &
                      'value,point', status='point evaluation enclosure at stored Cartesian points')
    call result%write('moments', fine, 'normalized', 'V,mean p,mean B squared,beta', &
                      'moment', status='quadrature reference; continuous bound unavailable')
    call result%write('moment_refinement', abs(fine - coarse), 'normalized', &
                      'difference under doubled quadrature', 'moment', status='indicator')
    call result%write('iota', iota, '1', 'axis,edge,edge refinement difference', &
                      'rotation', status='rotation reference; ODE and finite turn bound unavailable')
    call result%write('boundary_sampled_error', error, 'm', &
        'staggered Fourier reconstruction error', &
                      status='sampled indicator; continuous boundary bound unavailable')
    call result%write('boundary_coefficients', coefficients, 'm', 'Rcos,Rsin,Zcos,Zsin', &
                      'coef,tor,pol', status='Fourier quadrature; truncation indicator reported')
    call result%close()
    open (newunit=unit, file='output/input.nested_oracle', status='replace')
    write (unit, '(a)') '! Geometry-only boundary: no pressure/current/flux matching claim.'
    write (unit, '(a)') '! R,Z = sum c cos(m theta - n*NFP*phi) + s sin(...).'
    write (unit, '(a)') '&INDATA'
    write (unit, '(a)') ' NFP=2, LFREEB=F,'
    write (unit, '(a,l1,a)') ' LASYM=', abs(tau) > 0, ','
    write (unit, '(a,i0,a,i0,a)') ' MPOL=', mmax + 1, ', NTOR=', nmax, ','
    do m = 0, mmax
        do n = -nmax, nmax
            call coefficient('RBC', 1)
            call coefficient('RBS', 2)
            call coefficient('ZBC', 3)
            call coefficient('ZBS', 4)
        end do
    end do
    write (unit, '(a)') '/'
    close (unit)
    if (family == NESTED_INTEGER) then
        ! Phi(psi) is linear, so s=Phi/Phi_edge=psi/delta. Reuse the
        ! retained generated pressure setup and exact flux, never a fit.
        call nested_setup(eps, tau, delta, lam, edge_pressure, &
                          limit, extent, pa1, pa2)
        call nested_integer_averages(eps, tau, pa1, delta, volume, pmean, &
                                     bmean, beta, transform, flux)
        open (newunit=unit, file='output/input.nested_equilibrium', status='replace')
        write (unit, '(a)') '! Matched integer-family equilibrium: s=Phi/Phi_edge=psi/delta.'
        write (unit, '(a)') '! B in tesla, length in metres; normalized pressure divided by SI mu0.'
        write (unit, '(a)') '&INDATA'
        write (unit, '(a)') ' NFP=2, LFREEB=F, NCURR=0,'
        write (unit, '(a,l1,a)') ' LASYM=', abs(tau) > 0, ','
        write (unit, '(a,i0,a,i0,a)') ' MPOL=', mmax + 1, ', NTOR=', nmax, ','
        write (unit, '(a,es25.16,a)') ' PHIEDGE=', flux, ','
        write (unit, '(a)') " PMASS_TYPE='power_series', PIOTA_TYPE='power_series',"
        write (unit, '(a,es25.16,a,es25.16,a)') ' AM=', pa1, ',', edge_pressure - pa1, ','
        write (unit, '(a,es25.16,a)') ' AI=', transform, ','
        write (unit, '(a,es25.16,a)') ' PRES_SCALE=', 1/(4*pi*1e-7_dp), ','
        write (unit, '(a)') ' NS_ARRAY=17,33,65, FTOL_ARRAY=1e-8,1e-9,1e-10,'
        write (unit, '(a)') ' NITER_ARRAY=3000,5000,10000,'
        do m = 0, mmax
            do n = -nmax, nmax
                call coefficient('RBC', 1)
                call coefficient('RBS', 2)
                call coefficient('ZBC', 3)
                call coefficient('ZBS', 4)
            end do
        end do
        write (unit, '(a)') '/'
        close (unit)
        print '(a)', '[V] matched integer-family input: output/input.nested_equilibrium'
    end if
    print '(a,4es16.7)', '[V] V,<p>,<B^2>,beta (quadrature reference) ', fine
    print '(a,4es12.3)', '[V] quadrature refinement indicators ', abs(fine - coarse)
    print '(a,3es16.7)', '[V] axis,edge iota; edge refinement indicator ', iota
    print '(a,es12.3)', '[V] sampled boundary error (m) ', error
    print '(a)', '[V] geometry-only input: output/input.nested_oracle'
    print '(a)', '[V] result: output/nested_oracle_export.nc'
contains
    subroutine coefficient(name, row)
        character(*), intent(in)::name
        integer, intent(in)::row
        write (unit, '(1x,a,a,i0,a,i0,a,es25.16,a)') name, '(', n, ',', m, ')=', &
            coefficients(row, n, m), ','
    end subroutine
end program
