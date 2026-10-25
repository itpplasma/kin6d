program test_nonnested_relax
    ! Metriplectic relaxation oracles (docs/metriplectic-relaxation.md).
    !   taylor : periodic cylinder, perturbed Lundquist start at the Lundquist
    !            helicity; the global metric must return to the independently
    !            solved Beltrami state (LDLT + helicity bisection) and to the
    !            analytic Bessel field, then hand off to that Newton solve.
    !   ideal  : same start, pointwise metric; energy monotone, helicity and
    !            fluxes kept to roundoff, energy bounded below by the Taylor
    !            (Woltjer) minimum at that helicity.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_mhd_mesh, only: flux_mesh_t, flux_mesh_build
    use kin6d_equilibrium_nonnested, only: free_b_field_t, nn_build, nn_point, &
        nn_linear_solve
    use kin6d_nonnested_beltrami, only: beltrami_system_t, beltrami_report_t, &
        beltrami_prepare, beltrami_solve, beltrami_constraint
    use kin6d_nonnested_relax, only: relax_system_t, relax_options_t, &
        relax_report_t, relax_state_t, relax_prepare, relax_run, relax_state, &
        relax_store, relax_free, RELAX_IDEAL, RELAX_TAYLOR
    use fortfem_meridian_feec_n, only: meridian_load
    use kin6d_equilibrium_state, only: EQ_VALUES, EQ_BR, EQ_BPHI, EQ_BZ
    use kin6d_equilibrium_nonnested, only: nn_report_t, nn_solve, nn_vacuum_mode
    use kin6d_nonnested_relax, only: relax_pack_modes
    use nn_vacuum_fixture, only: vfield => field, ftor, nmodes, active_mode, &
        coil_mu0, make_coil, mean_psi, gauge_potential, annulus_mesh, &
        periodic_radius, full_periodic, stability, island_edge
    use fortfem_meridian_feec_n, only: meridian_lagrange_nodes
    use fortsparse, only: csc_matvec
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp), r0 = 6.0_dp, mu = 1.1_dp
    type(flux_mesh_t) :: mesh
    type(free_b_field_t) :: field, oracle
    type(beltrami_system_t) :: bsys
    type(beltrami_report_t) :: bref, bnewton
    type(relax_system_t) :: rs
    type(relax_options_t) :: options
    type(relax_report_t) :: report
    type(relax_state_t) :: s0, s1
    real(dp), allocatable :: rings(:), data(:), a(:), a0(:), da(:), v(:)
    real(dp) :: flux, kq(3), sigma, qa, qb, qc, e, scale, err, val(EQ_VALUES)
    real(dp) :: point(3), r, theta, bex(3), amp
    integer :: st, i, j, nr, modes, n1
    character(24) :: kind
    call get_command_argument(1, kind)
    if (trim(kind) == 'island') then
        call island_case()
        stop
    end if
    nr = 4; modes = 1
    allocate (rings(nr)); rings = [(real(i, dp)/nr, i=1, nr)]
    call flux_mesh_build(rings, 8, 256, 2, mapping, mesh, st)
    if (st /= 0) error stop 'mesh'
    call nn_build(mesh, 2, modes, [r0, 0.0_dp], 1.0_dp, 0.0_dp, 1.0_dp, field, st, &
        cylinder_r0=r0)
    if (st /= 0) error stop 'build'
    n1 = field%fs%sp(1)%ndof
    flux = 2*pi*bessel_j1(mu)/mu
    allocate (data((2*modes + 1)*n1)); data = 0
    call beltrami_prepare(field, flux, data, 2*modes + 3, bsys, st)
    if (st /= 0) error stop 'Beltrami prepare'
    call beltrami_solve(bsys, mu, field, bref, potential=a0)
    if (bref%status /= 0) error stop 'Beltrami oracle'
    oracle = field
    call relax_prepare(field, rs, st)
    if (st /= 0) error stop 'relax prepare'
    call relax_state(rs, RELAX_TAYLOR, a0, s0, st)
    if (st /= 0) error stop 'relax state'
    ! Same discrete forms as the Beltrami closure: energy and helicity agree.
    print '(a,4es24.15)', 'oracle W, K (Beltrami / relax) ', bref%energy, s0%energy, &
        bref%helicity, s0%helicity
    if (abs(s0%energy - bref%energy) > 1.0e-12_dp*abs(bref%energy)) &
        error stop 'energy form mismatch'
    if (abs(s0%helicity - bref%helicity) > 1.0e-12_dp*abs(bref%helicity)) &
        error stop 'helicity form mismatch'
    if (abs(s0%lambda - mu) > 1.0e-9_dp) error stop 'Beltrami multiplier'
    if (s0%force > 1.0e-8_dp) error stop 'oracle not force free'
    ! Smooth interior perturbation with n=0 and n=1 content.
    call perturbation(da)
    amp = 0.3_dp*maxval(abs(a0))/maxval(abs(da))
    a = a0 + amp*da
    ! Restore the Lundquist helicity exactly along v (quadratic in sigma).
    v = merge(0.0_dp, a0, rs%fixed)
    do i = 1, 3
        call relax_state(rs, RELAX_TAYLOR, a + (i - 2)*v, s1, st)
        kq(i) = s1%helicity - s0%helicity
    end do
    qa = (kq(3) + kq(1))/2 - kq(2); qb = (kq(3) - kq(1))/2; qc = kq(2)
    if (qb*qb - 4*qa*qc < 0) error stop 'no helicity-matching start'
    ! Smaller root of qa s^2 + qb s + qc = 0, cancellation-free.
    sigma = -2*qc/(qb + sign(sqrt(qb*qb - 4*qa*qc), qb))
    a = a + sigma*v
    call relax_state(rs, RELAX_TAYLOR, a, s1, st)
    print '(a,3es14.6)', 'start sigma, helicity defect, energy excess ', sigma, &
        (s1%helicity - s0%helicity)/s0%helicity, (s1%energy - s0%energy)/s0%energy
    if (abs(s1%helicity - s0%helicity) > 1.0e-12_dp*abs(s0%helicity)) &
        error stop 'helicity-matched start'
    if (s1%energy <= s0%energy) error stop 'start below the Taylor minimum'
    select case (trim(kind))
    case ('taylor')
        options%kind = RELAX_TAYLOR; options%dt = 1.0e-3_dp; options%dt_max = 1.0e2_dp
        options%handoff = 1.0e-5_dp; options%max_steps = 400
        call relax_run(rs, a, options, report)
        call summary()
        if (.not. report%handoff) error stop 'Taylor relaxation did not hand off'
        if (report%helicity_drift > 1.0e-12_dp) error stop 'helicity drift'
        if (report%flux_drift > 1.0e-13_dp .or. report%trace_drift > 1.0e-13_dp) &
            error stop 'flux drift'
        if (report%energy_increase > 1.0e-14_dp) error stop 'energy increased'
        if (report%identity_defect > 1.0e-8_dp .or. report%law_defect > 1.0e-13_dp) &
            error stop 'discrete energy law'
        if (abs(report%final%lambda - mu) > 1.0e-5_dp) error stop 'Taylor multiplier'
        call relax_store(rs, a, field, st)
        if (st /= 0) error stop 'store'
        e = maxval(abs(field%b - oracle%b))/maxval(abs(oracle%b))
        print '(a,es12.4)', 'max |B - B_Beltrami| / max |B| ', e
        if (e > 1.0e-4_dp) error stop 'relaxed field differs from Beltrami state'
        ! Newton hand-off: the helicity-constrained Beltrami solve started
        ! from the relaxed multiplier.
        call beltrami_constraint(bsys, 'helicity', report%final%helicity, &
            [report%final%lambda - 1.0e-3_dp, report%final%lambda + 1.0e-3_dp], &
            1.0e-12_dp, field, bnewton)
        print '(a,3es14.6)', 'Newton mu, |mu - lambda|, residual ', bnewton%mu, &
            abs(bnewton%mu - report%final%lambda), bnewton%residual
        if (bnewton%status /= 0 .or. bnewton%residual > 1.0e-9_dp) &
            error stop 'Newton hand-off'
        if (abs(bnewton%mu - mu) > 1.0e-9_dp) error stop 'Newton multiplier'
        ! Analytic Lundquist field B = c (0, J1(mu r), J0(mu r)) with the
        ! harmonic flux normalisation c = 1.
        err = 0; scale = 0
        do i = 1, 12
            do j = 1, 12
                r = 0.9_dp*(i - 0.5_dp)/12; theta = 2*pi*(j - 0.3_dp)/12
                point = [r0 + r*cos(theta), r*sin(theta), 0.7_dp]
                call nn_point(field, point, val, st)
                if (st /= 0) error stop 'sample'
                bex = [bessel_j1(mu*r)*sin(theta), bessel_j0(mu*r), &
                    -bessel_j1(mu*r)*cos(theta)]
                err = err + sum(([val(EQ_BR), val(EQ_BPHI), val(EQ_BZ)] - bex)**2)
                scale = scale + sum(bex**2)
            end do
        end do
        print '(a,es12.4)', 'relative RMS error to analytic Lundquist ', &
            sqrt(err/scale)
        if (sqrt(err/scale) > 1.0e-2_dp) error stop 'analytic Lundquist'
    case ('ideal')
        options%kind = RELAX_IDEAL; options%dt = 1.0e-3_dp; options%dt_max = 1.0_dp
        options%handoff = 1.0e-4_dp; options%max_steps = 40
        call relax_run(rs, a, options, report)
        call summary()
        if (report%status /= 0) error stop 'ideal relaxation failed'
        if (report%helicity_drift > 1.0e-12_dp) error stop 'helicity drift'
        if (report%flux_drift > 1.0e-13_dp .or. report%trace_drift > 1.0e-13_dp) &
            error stop 'flux drift'
        if (report%energy_increase > 1.0e-14_dp) error stop 'energy increased'
        ! The law holds to the B1 Picard tolerance (outer_tol), helicity to
        ! roundoff through the scalar projection.
        if (report%identity_defect > 1.0e-6_dp .or. report%law_defect > 1.0e-13_dp) &
            error stop 'discrete energy law'
        if (report%projection > 1.0e-6_dp) error stop 'helicity projection too large'
        if (report%final%force > 0.01_dp*report%initial%force) &
            error stop 'force not reduced a hundredfold'
        ! Woltjer bound: the Taylor state minimises W at fixed helicity/flux.
        if (report%final%energy < s0%energy*(1 - 1.0e-12_dp)) &
            error stop 'energy below the Taylor minimum'
    case default
        error stop 'usage: test_nonnested_relax taylor|ideal|island'
    end select
    call relax_free(rs)
    print '(a)', 'PASS'
contains
    subroutine island_case()
        ! NN1 vacuum helical-coil island (nonnested_vacuum_helical_coil_island
        ! setup at its finest mesh level). A small interior perturbation is
        ! relaxed with the ideal (topology-preserving) metric; the island
        ! chain must survive with its O/X points and width, the invariants
        ! must hold, and the vacuum Newton solve from the relaxed state must
        ! return to the vacuum within the NN1 residual gate.
        type(flux_mesh_t) :: fm
        type(nn_report_t) :: nrep
        type(free_b_field_t) :: vacuum
        complex(dp), allocatable :: bc(:), seed(:), sol0(:), soln(:), an(:, :)
        complex(dp), allocatable :: sol(:)
        real(dp), allocatable :: nodes(:, :)
        real(dp) :: o0(2), x0(2), o1(2), x1(2), os(2), xs(2), w0, w1, ws, beta
        integer :: i, n
        call make_coil()
        call annulus_mesh(16, 2, fm)
        call nn_build(fm, 2, nmodes, [r0, 0.0_dp], 1.0_dp, ftor, coil_mu0, vfield, st)
        if (st /= 0) error stop 'vacuum build'
        call meridian_lagrange_nodes(vfield%fs, nodes)
        allocate (bc(vfield%fs%n0), seed(vfield%fs%n0)); bc = 0; seed = 0
        do i = 1, size(bc)
            if (vfield%fs%sp(0)%bnd(i)) bc(i) = mean_psi(nodes(:, i))
        end do
        call nn_solve(vfield, 0, bc, seed, nrep, st, solution=sol0)
        if (st /= 0) error stop 'vacuum n0'
        allocate (an(vfield%fs%nn, nmodes))
        do n = 1, nmodes
            active_mode = n
            call nn_vacuum_mode(vfield, n, gauge_potential, nrep, st, solution=soln)
            if (st /= 0) error stop 'vacuum mode'
            an(:, n) = soln
        end do
        vacuum = vfield
        call island_geometry(o0, x0, w0)
        call relax_prepare(vfield, rs, st)
        if (st /= 0) error stop 'relax prepare'
        call relax_pack_modes(rs, sol0, an, a0)
        call relax_store(rs, a0, vfield, st)
        e = maxval(abs(vfield%b - vacuum%b))/maxval(abs(vacuum%b))
        print '(a,es12.4)', 'packed potential reproduces vacuum B to ', e
        if (st /= 0 .or. e > 1.0e-12_dp) error stop 'mode packing'
        call relax_state(rs, RELAX_IDEAL, a0, s0, st)
        print '(a,es12.4)', 'vacuum force ', s0%force
        ! Interior perturbation with max |curl dA| = 1e-3 max |curl A_vac|
        ! (the coil field; the toroidal harmonic is not part of curl A).
        n1 = rs%n1; modes = nmodes
        field = vfield
        call perturbation(da)
        a = a0 + 1.0e-3_dp*maxval(abs(csc_matvec(rs%curl, a0))) &
            /maxval(abs(csc_matvec(rs%curl, da)))*da
        ! Island of the perturbed start: the ideal flow must keep this one.
        call relax_store(rs, a, vfield, st)
        if (st /= 0) error stop 'store start'
        call island_geometry(os, xs, ws)
        beta = (ftor/r0)**2
        options%kind = RELAX_IDEAL; options%dt = 1.0e-4_dp/beta
        options%dt_max = 1.0e-2_dp/beta; options%handoff = 0.1_dp*s0%force
        options%max_steps = 12; options%verbose = .true.
        call relax_run(rs, a, options, report)
        call summary()
        if (report%status /= 0) error stop 'island relaxation'
        if (report%helicity_drift > 1.0e-12_dp) error stop 'helicity drift'
        if (report%flux_drift > 1.0e-13_dp .or. report%trace_drift > 1.0e-13_dp) &
            error stop 'flux drift'
        if (report%energy_increase > 1.0e-14_dp) error stop 'energy increased'
        if (report%final%force > 0.1_dp*report%initial%force) &
            error stop 'island force not reduced tenfold'
        call relax_store(rs, a, vfield, st)
        if (st /= 0) error stop 'store'
        call island_geometry(o1, x1, w1)
        print '(a,3es12.4)', 'start vs vacuum O/X shift, width change (m) ', &
            norm2(os - o0), norm2(xs - x0), ws - w0
        print '(a,3es12.4)', 'relaxed vs start O/X shift, width change (m) ', &
            norm2(o1 - os), norm2(x1 - xs), w1 - ws
        ! Gates of nonnested_vacuum_helical_coil_island: 2 mm O/X location,
        ! 3 mm between trace refinements of the width.
        if (max(norm2(o1 - os), norm2(x1 - xs)) > 2.0e-3_dp) error stop 'O/X moved'
        if (abs(w1 - ws) > 3.0e-3_dp) error stop 'island width changed'
        ! Newton hand-off: the NN1 vacuum solves seeded by the relaxed state.
        call nn_solve(vfield, 0, bc, cmplx(a(rs%nn + 1:rs%n1), 0.0_dp, dp), &
            nrep, st, solution=sol)
        if (st /= 0 .or. nrep%residual > 1.0e-9_dp) error stop 'Newton n0'
        print '(a,2es12.4)', 'Newton n=0 energy before/after ', nrep%energy_before, &
            nrep%energy_after
        do n = 1, nmodes
            active_mode = n
            call nn_vacuum_mode(vfield, n, gauge_potential, nrep, st)
            if (st /= 0 .or. nrep%residual > 1.0e-9_dp) error stop 'Newton mode'
        end do
        e = maxval(abs(vfield%b - vacuum%b))/maxval(abs(vacuum%b))
        print '(a,es12.4)', 'Newton state vs vacuum ', e
        if (e > 1.0e-10_dp) error stop 'Newton did not return to the vacuum'
        call relax_free(rs)
        print '(a)', 'PASS'
    end subroutine

    subroutine island_geometry(opoint, xpoint, width)
        real(dp), intent(out) :: opoint(2), xpoint(2), width
        real(dp) :: lower, upper, tr, defect
        call full_periodic(periodic_radius(5.49_dp, 5.54_dp, .false., 512), &
            opoint, .false.)
        call stability(opoint, .false., tr, defect)
        if (abs(tr) >= 2) error stop 'O point lost'
        call full_periodic(periodic_radius(6.69_dp, 6.72_dp, .false., 512), &
            xpoint, .false.)
        call stability(xpoint, .false., tr, defect)
        if (tr <= 2) error stop 'X point lost'
        !$omp parallel sections default(shared)
        !$omp section
        lower = island_edge(0.36_dp, 0.48_dp, .false., 32.0_dp, 256)
        !$omp section
        upper = island_edge(0.63_dp, 0.48_dp, .false., 32.0_dp, 256)
        !$omp end parallel sections
        width = upper - lower
        print '(a,5es14.6)', 'O(R,Z), X(R,Z), width ', opoint, xpoint, width
    end subroutine

    subroutine summary()
        print '(a,a)', 'status: ', trim(report%message)
        print '(a,i0,a,i0,a,i0,a,i0)', 'steps ', report%steps, ' rejected ', &
            report%rejected, ' outer ', report%outer, ' krylov ', report%krylov
        print '(a,2es14.6)', 'energy initial/final ', report%initial%energy, &
            report%final%energy
        print '(a,2es14.6)', 'force initial/final ', report%initial%force, &
            report%final%force
        print '(a,es14.6)', 'final lambda ', report%final%lambda
        print '(a,3es12.4)', 'helicity/flux/trace drift ', report%helicity_drift, &
            report%flux_drift, report%trace_drift
        print '(a,4es12.4)', 'energy increase, law defect rel/abs, projection ', &
            report%energy_increase, report%identity_defect, report%law_defect, &
            report%projection
    end subroutine

    subroutine perturbation(d)
        real(dp), allocatable, intent(out) :: d(:)
        complex(dp), allocatable :: load(:), x(:)
        integer :: b, st
        allocate (d((2*modes + 1)*n1)); d = 0
        do b = 0, 1
            j = b
            call meridian_load(field%fs, 1, bump, load, st)
            if (st /= 0) error stop 'perturbation load'
            allocate (x(size(load)))
            call nn_linear_solve(field%fs%mass(1), field%fs%sp(1)%bnd, load, x, st)
            if (st /= 0) error stop 'perturbation projection'
            d(b*n1 + 1:(b + 1)*n1) = real(x, dp)
            deallocate (x)
        end do
    end subroutine

    subroutine bump(x, w)
        real(dp), intent(in) :: x(2)
        complex(dp), intent(out) :: w(3)
        real(dp) :: y, z, q
        y = x(1) - r0; z = x(2); q = max(1 - y*y - z*z, 0.0_dp)
        if (j == 0) then
            w = cmplx([q*z*(1 + y), q*q*(1 + 0.5_dp*y), -q*y], 0.0_dp, dp)
        else
            w = cmplx([q*(z + 0.3_dp), q*y*z, q*(y - 0.2_dp)], 0.0_dp, dp)
        end if
    end subroutine

    subroutine mapping(u, theta, x)
        real(dp), intent(in) :: u, theta
        real(dp), intent(out) :: x(2)
        x = [r0 + u*cos(theta), u*sin(theta)]
    end subroutine
end program
