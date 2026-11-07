program nested_equilibrium
    ! Fixed-boundary nested 3D equilibrium from a boundary, the toroidal flux and
    ! the profiles iota(s), p(s) in the normalized toroidal flux s, solved by the
    ! native map Newton (docs/nested-equilibrium-input.md). Writes the solved
    ! state's global quantities for cross-code comparison.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use kin6d_ideal_seed_3d, only: toroidal_profile_t, toroidal_seed_t, &
        ideal_seed_3d_init, ideal_seed_disk_mesh
    use kin6d_ideal_map, only: ideal_map_t, ideal_map_solve, ideal_map_solve_report_t, &
        ideal_map_energy, ideal_map_reduce, ideal_map_prolong_correction
    use kin6d_map_3d_equilibrium, only: map_3d_coordinates, map_3d_sample
    use kin6d_equilibrium_state, only: EQ_VALUES, EQ_BR, EQ_BZ, EQ_BPHI, EQ_P
    use kin6d_grad_shafranov, only: gs_state_t
    use kin6d_curved_triangle, only: triangle_locate
    use kin6d_result_io, only: run_t, result_file_t
    use fortnum_quadrature, only: gauss_legendre
    implicit none
    integer, parameter :: mmax = 64, nmax = 32, kmax = 32
    real(dp), parameter :: pi = acos(-1.0_dp), mu0 = 4e-7_dp*pi
    integer :: nfp, mpol, ntor, boundary_nodes, degree, max_iterations, angle_mpol
    integer :: unit, ios, status, solve_status, step, steps, nodes, order
    integer :: ladder_nodes(8), ladder_degree(8)
    real(dp) :: phiedge, tolerance
    real(dp) :: rbc(-nmax:nmax, 0:mmax), rbs(-nmax:nmax, 0:mmax)
    real(dp) :: zbc(-nmax:nmax, 0:mmax), zbs(-nmax:nmax, 0:mmax)
    real(dp) :: iota(0:kmax), pressure(0:kmax)
    real(dp), allocatable :: boundary(:, :, :), seed_x(:, :, :), coarse_seed(:, :, :)
    real(dp) :: residual, reference, volume, area, magnetic, kinetic, raxis
    real(dp) :: iota_axis, iota_edge, p_axis, unused
    type(run_t) :: run
    type(result_file_t) :: result
    type(gs_state_t) :: mesh
    type(ideal_map_t) :: map, coarse
    type(toroidal_seed_t) :: seed
    type(toroidal_profile_t) :: profile
    type(ideal_map_solve_report_t) :: report
    character(256) :: message
    namelist /nested_input/ nfp, mpol, ntor, rbc, rbs, zbc, zbs, phiedge, &
        iota, pressure, boundary_nodes, degree, tolerance, max_iterations, angle_mpol, &
        ladder_nodes, ladder_degree

    nfp = 1; mpol = 1; ntor = 0; phiedge = 1
    rbc = 0; rbs = 0; zbc = 0; zbs = 0; iota = 0; pressure = 0
    boundary_nodes = 16; degree = 3; tolerance = 1e-8_dp; max_iterations = 200
    angle_mpol = -1; ladder_nodes = 0; ladder_degree = 0
    call run%start('nested_equilibrium')
    call run%open_config(unit)
    ios = 0; message = ''
    if (run%has_config()) read (unit, nml=nested_input, iostat=ios, iomsg=message)
    call run%close_config(unit, 'nested_input', ios, message)
    if (run%argument_count() /= 0) error stop 'usage: nested_equilibrium [--config FILE]'
    if (nfp < 1 .or. mpol < 1 .or. mpol > mmax .or. ntor < 0 .or. ntor > nmax) &
        error stop 'nested_equilibrium: mpol, ntor or nfp out of range'
    if (.not. phiedge > 0 .or. .not. ieee_is_finite(phiedge)) &
        error stop 'nested_equilibrium: phiedge must be positive'
    call run%open_resolved(unit)
    write (unit, nml=nested_input)
    call run%close_resolved(unit)

    ! boundary(1:4, n, m) = (rc, rs, zc, zs) of cos/sin(m theta - n nfp phi)
    allocate (boundary(4, -ntor:ntor, 0:mpol))
    boundary(1, :, :) = rbc(-ntor:ntor, 0:mpol); boundary(2, :, :) = rbs(-ntor:ntor, 0:mpol)
    boundary(3, :, :) = zbc(-ntor:ntor, 0:mpol); boundary(4, :, :) = zbs(-ntor:ntor, 0:mpol)
    profile%iota = iota(0:last(iota))
    ! the solver uses mu0 = 1 pressure units
    profile%pressure = mu0*pressure(0:last(pressure))
    call profile%evaluate(0.0_dp, iota_axis, p_axis, unused)
    call profile%evaluate(1.0_dp, iota_edge, unused, unused)

    ! Resolution ladder (docs/equilibrium-harness.md): each level starts from its
    ! own seed plus the correction of the solved coarser level and stops at the
    ! same absolute gradient relative to its uncorrected seed.
    steps = count(ladder_nodes > 0)
    do step = 1, steps + 1
        nodes = boundary_nodes; order = degree
        if (step <= steps) then
            nodes = ladder_nodes(step); order = ladder_degree(step)
        end if
        call ideal_seed_disk_mesh(mesh, nodes, order, status)
        if (status /= 0) error stop 'nested_equilibrium: disk mesh'
        if (angle_mpol >= 0) then
            call ideal_seed_3d_init(map, seed, mesh, boundary, nfp, phiedge, profile, &
                status, angle_mpol=angle_mpol)
        else
            call ideal_seed_3d_init(map, seed, mesh, boundary, nfp, phiedge, profile, status)
        end if
        if (status /= 0) error stop 'nested_equilibrium: seed'
        reference = free_gradient()
        seed_x = map%x
        if (step > 1) then
            call ideal_map_prolong_correction(map, coarse, coarse_seed, status)
            if (status /= 0) error stop 'nested_equilibrium: ladder correction'
        end if
        call ideal_map_solve(map, tolerance, max_iterations, report, solve_status, &
            reference_gradient=reference)
        coarse = map; coarse_seed = seed_x
    end do
    residual = free_gradient()/max(tiny(1.0_dp), reference)
    call integrals(volume, area, magnetic, kinetic, status)
    if (status /= 0) error stop 'nested_equilibrium: integrals'
    call axis_radius(raxis, status)
    if (status /= 0) error stop 'nested_equilibrium: axis'

    call result%create(run)
    call result%attribute('physical_model', 'fixed-boundary nested ideal MHD, prescribed iota')
    call result%attribute('converged', merge('true ', 'false', &
        solve_status == 0 .and. residual <= tolerance))
    call result%write('volume_p', volume, 'm^3', 'plasma volume')
    call result%write('aminor_p', sqrt(area/pi), 'm', 'sqrt of mean cross-section area / pi')
    call result%write('rmajor_p', volume/(2*pi*area), 'm', 'volume / (2 pi mean area)')
    call result%write('aspect', volume/(2*pi*area)/sqrt(area/pi), '1', 'rmajor_p / aminor_p')
    call result%write('raxis_cc', raxis, 'm', 'toroidal mean of the axis major radius')
    call result%write('iotaf_axis', iota_axis, '1', 'rotational transform on axis')
    call result%write('iotaf_edge', iota_edge, '1', 'rotational transform at the edge')
    call result%write('pressure_axis', p_axis/mu0, 'Pa', 'pressure on axis')
    call result%write('wb', magnetic/(4*pi**2), 'T^2 m^3', &
        'VMEC normalisation int B^2/2 dV / (2 pi)^2 (mu0 = 1)')
    call result%write('betatotal', kinetic/magnetic, '1', 'int p dV / int B^2/2 dV')
    call result%write('force_residual', residual, '1', 'relative free reduced-gradient l2')
    call result%write('tolerance', tolerance, '1', 'requested force_residual')
    call result%write('iterations', real(report%iterations, dp), '1', 'Newton iterations')
    call result%close()

contains

    real(dp) function free_gradient()
        ! l2 norm of the free reduced energy gradient of the current map
        real(dp), allocatable :: gradient(:), reduced(:)
        real(dp) :: energy
        integer :: local
        allocate (gradient(map%ndof), reduced(map%nred))
        call ideal_map_energy(map, map%x, energy, local, gradient)
        if (local /= 0) error stop 'nested_equilibrium: energy'
        call ideal_map_reduce(map, gradient, reduced)
        where (map%fixed) reduced = 0
        free_gradient = norm2(reduced)
    end function

    integer function last(c)
        ! highest nonzero coefficient, at least the constant
        real(dp), intent(in) :: c(0:)
        do last = ubound(c, 1), 1, -1
            if (c(last) /= 0) return
        end do
    end function

    subroutine integrals(volume, area, magnetic, kinetic, status)
        ! Gauss-Legendre on the collapsed reference triangle and the trapezoidal
        ! rule in phi, exact for the periodic Fourier dependence of the map.
        real(dp), intent(out) :: volume, area, magnetic, kinetic
        integer, intent(out) :: status
        integer, parameter :: nq = 8
        real(dp) :: x(nq), w(nq), xi(2), weight, q(3), d(3, 3), ref(2), g(2, 2)
        real(dp) :: point(3), values(EQ_VALUES), phi, jac, sums(4)
        integer :: e, i, j, k, nphi, failed, local
        call gauss_legendre(nq, x, w)
        x = (x + 1)/2; w = w/2
        nphi = max(16, 4*map%ntor + 4)
        sums = 0; failed = 0
        !$omp parallel do collapse(2) schedule(static) &
        !$omp private(e,i,j,k,xi,weight,phi,q,d,ref,g,point,values,jac,local) &
        !$omp reduction(+:sums,failed)
        do e = 1, map%cells
            do k = 0, nphi - 1
                phi = 2*pi*k/(nfp*nphi)
                do j = 1, nq
                    do i = 1, nq
                        xi = [x(i)*(1 - x(j)), x(j)]
                        weight = w(i)*w(j)*(1 - x(j))
                        call map_3d_coordinates(map, e, xi, phi, q, d, ref, g)
                        ! dV = R |det dq/d(xi, phi)|; the map also displaces phi
                        jac = abs(d(1, 1)*(d(2, 2)*d(3, 3) - d(2, 3)*d(3, 2)) &
                            - d(1, 2)*(d(2, 1)*d(3, 3) - d(2, 3)*d(3, 1)) &
                            + d(1, 3)*(d(2, 1)*d(3, 2) - d(2, 2)*d(3, 1)))
                        call map_3d_sample(map, seed, e, xi, phi, point, values, local)
                        failed = failed + local
                        sums = sums + weight*[q(1)*jac, jac, &
                            q(1)*jac*sum(values([EQ_BR, EQ_BZ, EQ_BPHI])**2)/2, &
                            q(1)*jac*values(EQ_P)]
                    end do
                end do
            end do
        end do
        !$omp end parallel do
        status = failed
        ! mean cross-section area = int dV/R / (2 pi)
        volume = 2*pi*sums(1)/nphi; area = sums(2)/nphi
        magnetic = 2*pi*sums(3)/nphi; kinetic = 2*pi*sums(4)/nphi
    end subroutine

    subroutine axis_radius(raxis, status)
        ! image of the centre of the reference disc (the label's O-point),
        ! averaged over phi
        real(dp), intent(out) :: raxis
        integer, intent(out) :: status
        real(dp) :: centre(2), xi(2), q(3), d(3, 3), ref(2), g(2, 2)
        integer :: e, k, nphi, local
        centre = (maxval(map%reference, dim=2) + minval(map%reference, dim=2))/2
        status = 1
        do e = 1, map%cells
            call triangle_locate(map%reference(:, map%conn(:, e)), centre, xi, local)
            if (local /= 0) cycle
            if (minval([xi, 1 - sum(xi)]) < -1e-10_dp) cycle
            nphi = max(16, 4*map%ntor + 4)
            raxis = 0
            do k = 0, nphi - 1
                call map_3d_coordinates(map, e, xi, 2*pi*k/(nfp*nphi), q, d, ref, g)
                raxis = raxis + q(1)/nphi
            end do
            status = 0
            return
        end do
    end subroutine
end program nested_equilibrium
