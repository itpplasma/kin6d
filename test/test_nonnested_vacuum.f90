program test_nonnested_vacuum
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use nn_vacuum_fixture, only: pi, r0, ftor, nmodes, field, active_mode, &
        coil_mu0, &
        make_coil, circle, mean_psi, gauge_potential, direct_field, orbit, &
        periodic_radius, full_periodic, stability, island_edge, annulus_mesh
    use kin6d_mhd_mesh, only: flux_mesh_t
    use fortfem_meridian_feec_n, only: meridian_lagrange_nodes
    use kin6d_equilibrium_nonnested, only: nn_report_t, nn_build, nn_solve, &
        nn_vacuum_mode, nn_point
    use kin6d_equilibrium_state, only: EQ_VALUES, EQ_BR, EQ_BZ, EQ_BPHI
    use kin6d_nonnested_residual, only: nn_residual
    use kin6d_equilibrium_residual, only: equilibrium_residual_t
    implicit none
    type(flux_mesh_t) :: fm
    type(nn_report_t) :: report
    type(equilibrium_residual_t) :: stress
    real(dp), allocatable :: nodes(:, :)
    complex(dp), allocatable :: bc(:), seed(:)
    real(dp) :: eb(3), em(3), d3(3), points(3, 64), exact(3, 64), values(3, 64)
    real(dp) :: val(EQ_VALUES), finish(2), reference(2), span, width, width_ref
    real(dp) :: oref, xref, oh, xh, trace_o, trace_x, defect, lower, upper
    real(dp) :: opoint(2), xpoint(2), oexact(2), xexact(2)
    real(dp) :: lower_ref, upper_ref, width_long, integration_error, rate
    real(dp) :: feec_step_error, width_ref_long, scan_resolution, map_error
    integer :: lev, nr, st, i, n
    call make_coil()
    oref = periodic_radius(5.49_dp, 5.54_dp, .true., 512)
    xref = periodic_radius(6.69_dp, 6.72_dp, .true., 512)
    call full_periodic(oref, oexact, .true.)
    call stability(oexact, .true., trace_o, defect)
    print '(a,3es14.6)', 'direct O, trace, closure ', oref, trace_o, defect
    if (abs(trace_o) >= 2 .or. defect > 1.0e-6_dp) error stop 'direct O point'
    call full_periodic(xref, xexact, .true.)
    call stability(xexact, .true., trace_x, defect)
    print '(a,3es14.6)', 'direct X, trace, closure ', xref, trace_x, defect
    if (trace_x <= 2 .or. defect > 1.0e-6_dp) error stop 'direct X point'
    !$omp parallel sections default(shared)
    !$omp section
    lower_ref = island_edge(0.36_dp, 0.48_dp, .true., 32.0_dp, 256)
    !$omp section
    upper_ref = island_edge(0.63_dp, 0.48_dp, .true., 32.0_dp, 256)
    !$omp end parallel sections
    width_ref = upper_ref - lower_ref
    !$omp parallel sections default(shared)
    !$omp section
    lower = island_edge(0.36_dp, 0.48_dp, .true., 64.0_dp, 512)
    !$omp section
    upper = island_edge(0.63_dp, 0.48_dp, .true., 64.0_dp, 512)
    !$omp end parallel sections
    width_long = upper - lower
    width_ref_long = width_long
    print '(a,4es14.6)', 'direct island edges, width, long-run width ', &
        lower_ref, upper_ref, width_ref, width_long
    if (abs(width_long - width_ref) > 0.003_dp) error stop 'unresolved island width'
    call orbit([5.52_dp, 0.0_dp], 1.0_dp, 256, .true., reference, span)
    call orbit([5.52_dp, 0.0_dp], 1.0_dp, 512, .true., finish, span)
    integration_error = norm2(reference - finish)
    print '(a,es14.6)', 'direct integration step-doubling ', integration_error
    if (integration_error > 1.0e-6_dp) error stop 'direct integration'
    do i = 1, 64
        call circle(0.4_dp + 0.07_dp*mod(i, 8), 2*pi*(i - 0.31_dp)/64, points(:2, i))
        points(3, i) = 2*pi*(i - 0.17_dp)/19
        call direct_field(points(:2, i), points(3, i), exact(:, i))
        exact(2, i) = exact(2, i) - ftor/points(1, i)
    end do
    do lev = 1, 3
        nr = 4*2**(lev - 1)
        if (allocated(bc)) deallocate (bc, seed)
        call annulus_mesh(nr, 2, fm)
        call nn_build(fm, 2, nmodes, [r0, 0.0_dp], 1.0_dp, ftor, coil_mu0, &
            field, st)
        if (st /= 0) error stop 'vacuum FEEC'
        call meridian_lagrange_nodes(field%fs, nodes)
        allocate (bc(field%fs%n0), seed(field%fs%n0))
        bc = 0; seed = 0
        !$omp parallel do default(shared) private(i) schedule(static)
        do i = 1, size(bc)
            if (field%fs%sp(0)%bnd(i)) bc(i) = mean_psi(nodes(:, i))
        end do
        !$omp end parallel do
        call nn_solve(field, 0, bc, seed, report, st)
        if (st /= 0) error stop 'vacuum n0 solve'
        do n = 1, nmodes
            active_mode = n
            call nn_vacuum_mode(field, n, gauge_potential, report, st)
            if (st /= 0) error stop 'vacuum mode solve'
            if (report%divergence > 1.0e-12_dp) error stop 'vacuum closed form'
        end do
        call nn_solve(field, 0, bc, seed, report, st, pprime=1.0_dp)
        if (st == 0) error stop '3D finite pressure must be rejected'
        do i = 1, 64
            call nn_point(field, points(:, i), val, st)
            if (st /= 0) error stop 'vacuum sampling'
            values(:, i) = [val(EQ_BR), val(EQ_BPHI) - ftor/points(1, i), val(EQ_BZ)]
        end do
        eb(lev) = sqrt(sum((values - exact)**2)/sum(exact**2))
        map_error = 0
        !$omp parallel do default(shared) private(i, reference, finish, span) &
        !$omp& reduction(max:map_error) schedule(static)
        do i = 1, 6
            call orbit([5.52_dp + 0.008_dp*i, 0.0_dp], 1.0_dp, 256, &
                .true., reference, span)
            call orbit([5.52_dp + 0.008_dp*i, 0.0_dp], 1.0_dp, 256, &
                .false., finish, span)
            map_error = max(map_error, norm2(finish - reference))
        end do
        !$omp end parallel do
        em(lev) = map_error
        call nn_residual(field, 1, stress, st)
        if (st /= 0) error stop 'vacuum D3'
        d3(lev) = stress%normalized
        print '(a,i0,a,3es14.6)', 'nr=', nr, ' B / map metres / D3 ', &
            eb(lev), em(lev), d3(lev)
    end do
    rate = log(eb(1)/eb(3))/log(4.0_dp)
    print '(a,3f9.3)', 'orders field/map/D3 ', rate, &
        log(em(1)/em(3))/log(4.0_dp), log(d3(1)/d3(3))/log(4.0_dp)
    if (rate < 1.5_dp .or. eb(3) > 0.01_dp) error stop 'vacuum field convergence'
    if (em(3) > 0.01_dp .or. em(3) >= em(1)) error stop 'field-line map convergence'
    if (d3(3) >= d3(1)) error stop 'vacuum D3 convergence'
    oh = periodic_radius(5.49_dp, 5.54_dp, .false., 512)
    xh = periodic_radius(6.69_dp, 6.72_dp, .false., 512)
    call full_periodic(oh, opoint, .false.)
    call stability(opoint, .false., trace_o, defect)
    print '(a,4es14.6)', 'FEEC O(R,Z), trace, closure ', opoint, trace_o, defect
    if (abs(trace_o) >= 2) error stop 'FEEC O point'
    call full_periodic(xh, xpoint, .false.)
    call stability(xpoint, .false., trace_x, defect)
    print '(a,4es14.6)', 'FEEC X(R,Z), trace, closure ', xpoint, trace_x, defect
    if (trace_x <= 2) error stop 'FEEC X point'
    if (max(norm2(opoint - oexact), norm2(xpoint - xexact)) > 0.002_dp) &
        error stop 'O/X location'
    !$omp parallel sections default(shared)
    !$omp section
    lower = island_edge(0.36_dp, 0.48_dp, .false., 32.0_dp, 256)
    !$omp section
    upper = island_edge(0.63_dp, 0.48_dp, .false., 32.0_dp, 256)
    !$omp end parallel sections
    width = upper - lower
    print '(a,3es14.6)', 'FEEC edges and width ', lower, upper, width
    scan_resolution = (abs(0.36_dp - 0.48_dp) + abs(0.63_dp - 0.48_dp))/2**9
    print '(a,2es14.6)', 'width difference and scan resolution ', &
        abs(width - width_ref), scan_resolution
    if (abs(width - width_ref) > 0.004_dp) error stop 'island width'
    !$omp parallel sections default(shared)
    !$omp section
    lower = island_edge(0.36_dp, 0.48_dp, .false., 64.0_dp, 512)
    !$omp section
    upper = island_edge(0.63_dp, 0.48_dp, .false., 64.0_dp, 512)
    !$omp end parallel sections
    width_long = upper - lower
    print '(a,2es14.6)', '64-turn FEEC width and oracle difference ', &
        width_long, abs(width_long - width_ref_long)
    ! The bisected edge is a finite-turn escape classification of a chaotic layer,
    ! so the width erodes with the turn count and its last 1e-3 depends on roundoff
    ! (the 64-turn FEEC width is 0.17271 with system gfortran/reference BLAS on a
    ! Ryzen 5950X and 0.17446 with conda gfortran/OpenBLAS on EPYC 9455, each
    ! independent of thread count).  Test that FEEC erodes like the direct oracle
    ! instead of bounding the absolute 32 -> 64 turn change (1.2e-3 for the oracle
    ! itself); the measured excess, 2.1e-3 and 0.3e-3, stays below 0.003.
    if (abs((width_long - width) - (width_ref_long - width_ref)) > 0.003_dp) &
        error stop 'FEEC width resolution'
    if (abs(width_long - width_ref_long) > 0.004_dp) error stop 'long-run width'
    call orbit([5.52_dp, 0.0_dp], 1.0_dp, 512, .false., reference, span)
    call orbit([5.52_dp, 0.0_dp], 1.0_dp, 1024, .false., finish, span)
    feec_step_error = norm2(reference - finish)
    print '(a,es14.6)', 'FEEC integration step-doubling ', feec_step_error
    if (feec_step_error > em(3)) error stop 'integration dominates map comparison'
    print '(a)', 'PASS: NN1 vacuum island against direct Biot-Savart tracing'
end program
