program nested_equilibrium
    ! Fixed-boundary nested 3D equilibrium from a boundary, the toroidal flux and
    ! the profiles iota(s), p(s) in the normalized toroidal flux s, solved by the
    ! native map Newton (docs/nested-equilibrium-input.md). Writes the solved
    ! state's global quantities for cross-code comparison.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use milca_ideal_seed_3d, only: toroidal_profile_t, toroidal_seed_t, &
        ideal_seed_3d_init, ideal_seed_disk_mesh
    use milca_ideal_map, only: ideal_map_t, ideal_map_solve, ideal_map_solve_report_t, &
        ideal_map_energy, ideal_map_reduce, ideal_map_prolong_correction
    use milca_map_3d_equilibrium, only: map_3d_coordinates, map_3d_sample
    use milca_equilibrium_state, only: EQ_VALUES, EQ_BR, EQ_BZ, EQ_BPHI, EQ_P
    use milca_grad_shafranov, only: gs_state_t
    use milca_curved_triangle, only: triangle_locate
    use milca_result_io, only: run_t, result_file_t
    use fortnum_quadrature, only: gauss_legendre
    implicit none
    integer, parameter :: mmax = 64, nmax = 32, kmax = 32
    real(dp), parameter :: pi = acos(-1.0_dp), mu0 = 4e-7_dp*pi
    integer :: nfp, mpol, ntor, boundary_nodes, degree, max_iterations, angle_mpol
    integer :: unit, ios, status, solve_status, step, steps, nodes, order, n
    integer :: ladder_nodes(8), ladder_degree(8)
    real(dp) :: phiedge, tolerance
    real(dp) :: rbc(-nmax:nmax, 0:mmax), rbs(-nmax:nmax, 0:mmax)
    real(dp) :: zbc(-nmax:nmax, 0:mmax), zbs(-nmax:nmax, 0:mmax)
    real(dp) :: iota(0:kmax), pressure(0:kmax), current(0:kmax)
    real(dp) :: raxis_cc(0:nmax), raxis_cs(0:nmax), zaxis_cc(0:nmax), zaxis_cs(0:nmax)
    real(dp), allocatable :: axis(:, :)
    logical :: guessed = .false.
    logical :: prescribe_current
    integer :: max_current_iterations
    real(dp) :: current_tolerance, current_error, solve_tolerance
    integer :: current_iterations
    real(dp), allocatable :: boundary(:, :, :), seed_x(:, :, :), coarse_seed(:, :, :)
    real(dp) :: residual, reference, volume, area, magnetic, kinetic, raxis
    real(dp) :: iota_axis, iota_edge, p_axis, unused, ctor
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
        ladder_nodes, ladder_degree, prescribe_current, current, max_current_iterations, &
        current_tolerance, raxis_cc, raxis_cs, zaxis_cc, zaxis_cs

    nfp = 1; mpol = 1; ntor = 0; phiedge = 1
    rbc = 0; rbs = 0; zbc = 0; zbs = 0; iota = 0; pressure = 0
    boundary_nodes = 16; degree = 3; tolerance = 1e-8_dp; max_iterations = 200
    angle_mpol = -1; ladder_nodes = 0; ladder_degree = 0
    prescribe_current = .false.; current = 0; max_current_iterations = 30
    raxis_cc = 0; raxis_cs = 0; zaxis_cc = 0; zaxis_cs = 0
    current_tolerance = 1e-5_dp
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
    if (any([raxis_cc, raxis_cs, zaxis_cc, zaxis_cs] /= 0)) then
        ! VMEC's axis arrays follow its boundary's m = 0 convention, the
        ! layout cos/sin(-n nfp phi); negative n keep the boundary's values
        allocate (axis(4, -ntor:ntor), source=boundary(:, :, 0))
        do n = 0, ntor
            axis(:, n) = [raxis_cc(n), raxis_cs(n), zaxis_cc(n), zaxis_cs(n)]
        end do
        axis(2:4:2, 0) = 0
    end if
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
        call init_seed(status)
        if (status /= 0 .and. .not. guessed) then
            ! as VMEC does after a bad initial Jacobian: guess the axis
            call guess_axis()
            guessed = .true.
            call init_seed(status)
        end if
        if (status /= 0) then
            print '(a,i0)', 'nested_equilibrium: seed status ', status
            error stop 'nested_equilibrium: seed'
        end if
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
    current_error = 0; current_iterations = 0; solve_tolerance = tolerance
    if (prescribe_current) call match_current(current_error)
    call profile%evaluate(0.0_dp, iota_axis, unused, unused)
    call profile%evaluate(1.0_dp, iota_edge, unused, unused)
    residual = free_gradient()/max(tiny(1.0_dp), reference)
    call integrals(volume, area, magnetic, kinetic, status)
    if (status /= 0) error stop 'nested_equilibrium: integrals'
    call axis_radius(raxis, status)
    if (status /= 0) error stop 'nested_equilibrium: axis'
    ctor = enclosed_current(1.0_dp, status)
    if (status /= 0) error stop 'nested_equilibrium: current'

    call result%create(run)
    call result%attribute('physical_model', 'fixed-boundary nested ideal MHD, prescribed iota')
    call result%attribute('converged', merge('true ', 'false', &
        solve_status == 0 .and. residual <= tolerance .and. current_error <= current_tolerance))
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
    call result%write('ctor', ctor, 'A', 'enclosed toroidal current at the edge')
    if (prescribe_current) then
        call result%write('current_error', current_error, '1', &
            'largest enclosed-current error relative to the largest prescribed current')
        call result%write('current_iterations', real(current_iterations, dp), '1', &
            'iota updates to match the prescribed current')
    end if
    call result%write('force_residual', residual, '1', 'relative free reduced-gradient l2')
    call result%write('tolerance', tolerance, '1', 'requested force_residual')
    call result%write('iterations', real(report%iterations, dp), '1', 'Newton iterations')
    call result%write('energy_seconds', report%energy_s, 's', 'final solve: energy and gradient')
    call result%write('hessian_seconds', report%hessian_s, 's', 'final solve: Hessian assembly')
    call result%write('factor_seconds', report%factor_s, 's', 'final solve: factorisation and solves')
    call result%write('factorizations', real(report%factorizations, dp), '1', &
        'final solve: factorisations')
    call result%write('lu_fallbacks', real(report%lu_fallbacks, dp), '1', &
        'final solve: Cholesky failures refactored by pivoted LU')
    call result%write('analyses', real(report%analyses, dp), '1', 'final solve: symbolic analyses')
    call result%close()

contains

    subroutine match_current(error)
        ! Prescribed enclosed toroidal current I(s). For a fixed map I on each
        ! surface is affine in iota there; its slope is read on the solved map
        ! with the frozen form's iota shifted by delta, so every update is a
        ! Newton step on the current map, refitted and re-solved.
        real(dp), intent(out) :: error
        integer, parameter :: ns = 12
        real(dp), parameter :: delta = 1e-3_dp
        real(dp) :: s(ns), target(ns), measured(ns), shifted(ns), slope(ns)
        real(dp) :: iota_s(ns), scale, unused1, unused2, update(ns)
        real(dp), allocatable :: previous_iota(:)
        integer :: j, k, local
        do j = 1, ns
            s(j) = real(j, dp)/ns
            target(j) = sum(current*s(j)**[(k, k=0, kmax)])
        end do
        scale = maxval(abs(target))
        do k = 0, max_current_iterations
            call measure(s, measured)
            if (k == 0) then
                ! the current per unit iota of the starting map scales a
                ! current-free target, so that the error is an iota error there
                seed%profile%iota(1) = seed%profile%iota(1) + delta
                call measure(s, shifted)
                seed%profile%iota(1) = seed%profile%iota(1) - delta
                scale = max(scale, maxval(abs(shifted - measured))/delta, tiny(1.0_dp))
            end if
            error = maxval(abs(measured - target))/scale
            print '(a,i3,a,es10.3)', 'current iteration ', k, ': relative error ', error
            if (error <= current_tolerance .or. k == max_current_iterations) then
                ! the current of a loosely solved map is measured again on the
                ! fully solved one
                if (solve_tolerance <= tolerance .or. k == max_current_iterations) exit
                solve_tolerance = tolerance
                call ideal_map_solve(map, solve_tolerance, max_iterations, report, &
                    solve_status, reference_gradient=reference)
                call measure(s, measured)
                error = maxval(abs(measured - target))/scale
                print '(a,i3,a,es10.3)', 'current iteration ', k, ': relative error ', error
                if (error <= current_tolerance) exit
            end if
            current_iterations = k + 1
            seed%profile%iota(1) = seed%profile%iota(1) + delta
            call measure(s, shifted)
            seed%profile%iota(1) = seed%profile%iota(1) - delta
            slope = (shifted - measured)/delta
            do j = 1, ns
                call profile%evaluate(s(j), iota_s(j), unused1, unused2)
            end do
            ! KIN6D's frozen field needs iota without zeros: halve an update
            ! that would change its sign
            update = (target - measured)/slope
            do j = 1, 30
                if (all(iota_s*(iota_s + update) > 0)) exit
                update = update/2
            end do
            iota_s = iota_s + update
            previous_iota = profile%iota
            call refit(s, iota_s)
            if (.not. single_signed()) then
                profile%iota = previous_iota
                print '(a)', 'current iteration: the refitted iota changes sign'
                exit
            end if
            ! intermediate maps need only the accuracy of the current mismatch
            solve_tolerance = max(tolerance, min(1e-4_dp, 1e-3_dp*error))
            call resolve(local)
            if (local /= 0) error stop 'nested_equilibrium: current re-seed'
        end do
    end subroutine

    subroutine init_seed(status)
        integer, intent(out) :: status
        ! an unallocated axis is an absent optional argument
        if (angle_mpol >= 0) then
            call ideal_seed_3d_init(map, seed, mesh, boundary, nfp, phiedge, profile, &
                status, angle_mpol=angle_mpol, axis=axis)
        else
            call ideal_seed_3d_init(map, seed, mesh, boundary, nfp, phiedge, profile, &
                status, axis=axis)
        end if
    end subroutine

    subroutine guess_axis()
        ! In each toroidal plane the axis point that maximises the smallest
        ! Jacobian of the seed's harmonic extension (grid search, as VMEC's
        ! guess_axis), then fitted with the axis Fourier modes.
        integer, parameter :: nr = 16, nt = 96, ng = 41
        integer :: nplane, k, i, j, m, nn, a, b, info
        real(dp) :: phi, theta, rho, phase, xm(2, 0:mmax), dxm(2, 0:mmax), b0(2)
        real(dp) :: best, worst, cand(2), jac, rr, rt, zr, zt, lo(2), hi(2)
        real(dp), allocatable :: points(:, :), lsq(:, :), rhs(:, :), work(:)
        nplane = max(8, 2*ntor + 2)
        allocate (points(2, nplane), lsq(nplane, 2*ntor + 1), rhs(nplane, 2), work(64*nplane))
        do k = 1, nplane
            phi = 2*pi*(k - 1)/(nfp*nplane)
            lo = huge(1.0_dp); hi = -huge(1.0_dp)
            do j = 0, nt - 1
                theta = 2*pi*j/nt
                call plane_modes(theta, phi, xm, dxm)
                lo = min(lo, sum(xm, dim=2)); hi = max(hi, sum(xm, dim=2))
            end do
            call plane_modes(0.0_dp, phi, xm, dxm); b0 = xm(:, 0)
            best = -huge(1.0_dp)
            do a = 1, ng
                do b = 1, ng
                    cand = lo + (hi - lo)*[real(a, dp), real(b, dp)]/(ng + 1)
                    worst = huge(1.0_dp)
                    do j = 0, nt - 1
                        theta = 2*pi*j/nt
                        call plane_modes(theta, phi, xm, dxm)
                        do i = 1, nr
                            rho = real(i, dp)/nr
                            rr = -2*rho*(cand(1) - b0(1)); zr = -2*rho*(cand(2) - b0(2))
                            rt = 0; zt = 0
                            do m = 1, mpol
                                rr = rr + m*rho**(m - 1)*xm(1, m); zr = zr + m*rho**(m - 1)*xm(2, m)
                                rt = rt + rho**m*dxm(1, m); zt = zt + rho**m*dxm(2, m)
                            end do
                            jac = (rr*zt - rt*zr)/rho
                            worst = min(worst, jac)
                        end do
                    end do
                    if (worst > best) then
                        best = worst; points(:, k) = cand
                    end if
                end do
            end do
            ! cos and sin of (-n nfp phi), the m = 0 layout
            lsq(k, 1) = 1
            do nn = 1, ntor
                lsq(k, 2*nn) = cos(-nn*nfp*phi); lsq(k, 2*nn + 1) = sin(-nn*nfp*phi)
            end do
        end do
        rhs = transpose(points)
        call dgels('N', nplane, 2*ntor + 1, 2, lsq, nplane, rhs, nplane, work, size(work), info)
        if (info /= 0) error stop 'nested_equilibrium: axis fit'
        if (allocated(axis)) deallocate (axis)
        allocate (axis(4, -ntor:ntor), source=boundary(:, :, 0))
        axis(:, 0) = [rhs(1, 1), 0.0_dp, rhs(1, 2), 0.0_dp]
        do nn = 1, ntor
            axis(:, nn) = [rhs(2*nn, 1), rhs(2*nn + 1, 1), rhs(2*nn, 2), rhs(2*nn + 1, 2)]
            axis(:, -nn) = 0
        end do
    end subroutine

    subroutine plane_modes(theta, phi, xm, dxm)
        ! (R, Z) of each poloidal mode m of the boundary at (theta, phi) and its
        ! theta derivative
        real(dp), intent(in) :: theta, phi
        real(dp), intent(out) :: xm(2, 0:mmax), dxm(2, 0:mmax)
        real(dp) :: phase
        integer :: m, nn
        xm = 0; dxm = 0
        do m = 0, mpol
            do nn = -ntor, ntor
                phase = m*theta - nn*nfp*phi
                xm(:, m) = xm(:, m) + boundary([1, 3], nn, m)*cos(phase) &
                    + boundary([2, 4], nn, m)*sin(phase)
                dxm(:, m) = dxm(:, m) + m*(boundary([2, 4], nn, m)*cos(phase) &
                    - boundary([1, 3], nn, m)*sin(phase))
            end do
        end do
    end subroutine

    logical function single_signed()
        ! the seed's admission rule: iota without zeros on [0, 1]
        real(dp) :: value, first, unused1, unused2
        integer :: a
        call profile%evaluate(0.0_dp, first, unused1, unused2)
        single_signed = abs(first) > tiny(1.0_dp)
        do a = 1, 128
            call profile%evaluate(a/128.0_dp, value, unused1, unused2)
            single_signed = single_signed .and. value*first > 0
        end do
    end function

    subroutine measure(s, values)
        ! enclosed toroidal current on the surfaces s
        real(dp), intent(in) :: s(:)
        real(dp), intent(out) :: values(:)
        integer :: i, local
        do i = 1, size(s)
            values(i) = enclosed_current(s(i), local)
            if (local /= 0) error stop 'nested_equilibrium: current'
        end do
    end subroutine

    subroutine refit(s, values)
        ! least-squares iota monomials through the sampled surfaces and the value
        ! extrapolated to the axis
        real(dp), intent(in) :: s(:), values(:)
        integer, parameter :: order = 8
        real(dp) :: a(size(s) + 1, 0:order), rhs(size(s) + 1), work(64*(size(s) + 1))
        integer :: i, k, n, info
        n = size(s)
        do i = 1, n
            a(i, :) = s(i)**[(k, k=0, order)]
            rhs(i) = values(i)
        end do
        a(n + 1, :) = 0; a(n + 1, 0) = 1
        rhs(n + 1) = 2*values(1) - values(2)
        call dgels('N', n + 1, order + 1, 1, a, n + 1, rhs, n + 1, work, size(work), info)
        if (info /= 0) error stop 'nested_equilibrium: iota fit'
        profile%iota = rhs(1:order + 1)
    end subroutine

    subroutine resolve(local)
        ! re-seed the frozen form with the new iota and keep the solved map
        integer, intent(out) :: local
        real(dp), allocatable :: x_saved(:, :, :), angle_saved(:)
        x_saved = map%x
        if (allocated(map%angle)) angle_saved = map%angle
        call ideal_seed_disk_mesh(mesh, boundary_nodes, degree, local)
        call init_seed(local)
        if (local /= 0) return
        map%x = x_saved
        if (allocated(angle_saved)) map%angle = angle_saved
        ! the same absolute target as the first solve, relative to its seed
        call ideal_map_solve(map, solve_tolerance, max_iterations, report, solve_status, &
            reference_gradient=reference)
    end subroutine

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

    real(dp) function enclosed_current(s, status)
        ! (1/mu0) oint B.dl around a poloidal loop on the flux surface s at phi = 0;
        ! every such loop encloses the same current. Inside, the loop is the image
        ! of the reference circle; at s = 1 it runs along the discrete wall.
        real(dp), intent(in) :: s
        integer, intent(out) :: status
        integer, parameter :: nt = 512, ng = 8
        real(dp) :: centre(2), radius, y(2), xi(2), trial(2), point(3, 0:nt), b(3, 0:nt)
        real(dp) :: values(EQ_VALUES), dl(3), corner(2, 3), x(ng), w(ng), q(3), d(3, 3)
        real(dp) :: ref(2), g(2, 2), tangent(2), along(3), sense
        integer :: k, e, local, owner, i, a
        centre = (maxval(map%reference, dim=2) + minval(map%reference, dim=2))/2
        radius = sqrt(s)*(maxval(map%reference(1, :)) - minval(map%reference(1, :)))/2
        enclosed_current = 0; status = 1
        if (s >= 1) then
            corner = reshape([0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2, 3])
            call gauss_legendre(ng, x, w)
            x = (x + 1)/2; w = w/2
            do e = 1, map%cells
                do i = 1, 3
                    a = mod(i, 3) + 1
                    if (.not. (map%boundary(map%conn(i, e)) .and. map%boundary(map%conn(a, e)))) cycle
                    ! a wall edge, not a chord between two wall nodes
                    call map_3d_coordinates(map, e, (corner(:, i) + corner(:, a))/2, 0.0_dp, q, d, ref, g)
                    if (abs(norm2(ref - centre)/radius - 1) > 1e-3_dp) cycle
                    tangent = corner(:, a) - corner(:, i)
                    sense = sign(1.0_dp, (ref(1) - centre(1))*dot_product(g(2, :), tangent) &
                        - (ref(2) - centre(2))*dot_product(g(1, :), tangent))
                    do k = 1, ng
                        xi = corner(:, i) + x(k)*tangent
                        call map_3d_coordinates(map, e, xi, 0.0_dp, q, d, ref, g)
                        call map_3d_sample(map, seed, e, xi, 0.0_dp, point(:, 0), values, local)
                        if (local /= 0) return
                        along = matmul(d(:, 1:2), tangent); along(3) = along(3)*q(1)
                        enclosed_current = enclosed_current + sense*w(k)* &
                            dot_product(values([EQ_BR, EQ_BZ, EQ_BPHI]), along)
                    end do
                end do
            end do
        else
            do k = 0, nt
                y = centre + radius*[cos(2*pi*k/nt), sin(2*pi*k/nt)]
                owner = 0
                do e = 1, map%cells
                    call triangle_locate(map%reference(:, map%conn(:, e)), y, trial, local)
                    if (local /= 0) cycle
                    if (minval([trial, 1 - sum(trial)]) < -1e-10_dp) cycle
                    owner = e; xi = trial; exit
                end do
                if (owner == 0) return
                call map_3d_sample(map, seed, owner, xi, 0.0_dp, point(:, k), values, local)
                if (local /= 0) return
                b(:, k) = values([EQ_BR, EQ_BZ, EQ_BPHI])
            end do
            do k = 0, nt - 1
                ! (R, Z, phi) increments; the toroidal length is R dphi
                dl = point(:, k + 1) - point(:, k)
                dl(3) = dl(3)*(point(1, k + 1) + point(1, k))/2
                enclosed_current = enclosed_current + dot_product((b(:, k + 1) + b(:, k))/2, dl)
            end do
        end if
        enclosed_current = enclosed_current/mu0
        status = 0
    end function

    subroutine axis_radius(raxis, status)
        ! image of the centre of the reference disc (the label's O-point),
        ! averaged over the cylindrical angle
        real(dp), intent(out) :: raxis
        integer, intent(out) :: status
        real(dp) :: centre(2), xi(2), q(3), d(3, 3), ref(2), g(2, 2), weight
        integer :: e, k, nphi, local
        centre = (maxval(map%reference, dim=2) + minval(map%reference, dim=2))/2
        status = 1
        do e = 1, map%cells
            call triangle_locate(map%reference(:, map%conn(:, e)), centre, xi, local)
            if (local /= 0) cycle
            if (minval([xi, 1 - sum(xi)]) < -1e-10_dp) cycle
            nphi = max(16, 4*map%ntor + 4)
            ! mean over the cylindrical angle, which the map displaces: dphi_cyl = d(3, 3) dphi
            raxis = 0; weight = 0
            do k = 0, nphi - 1
                call map_3d_coordinates(map, e, xi, 2*pi*k/(nfp*nphi), q, d, ref, g)
                raxis = raxis + q(1)*d(3, 3); weight = weight + d(3, 3)
            end do
            raxis = raxis/weight
            status = 0
            return
        end do
    end subroutine
end program nested_equilibrium
