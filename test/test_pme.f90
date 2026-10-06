program test_pme
    ! Particle-mesh Ewald against independent oracles: B-spline partition of
    ! unity and derivative identities, the converged direct Ewald sum
    ! (kin6d_nbody, structure-factor reciprocal sum) for forces and energy,
    ! the cell-list and all-pairs real-space paths, the energy-gradient
    ! identity and finite-difference force tangents, fixed particles.
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use fortnum_status, only: fortnum_status_t, FORTNUM_OK
    use fortnum_rng, only: rng_t, rng_seed, rng_uniform
    use kin6d_species, only: units_t, species_t, reduced_units, electron
    use kin6d_ewald, only: ewald_parameters_t, ewald_select
    use kin6d_pme, only: pme_parameters_t, pme_select, pme_estimate
    use kin6d_ocp_leaves, only: pme_bspline4, pme_bspline6, pme_bspline8
    use kin6d_nbody, only: nbody_system_t, nbody_periodic, nbody_periodic_mesh, &
        nbody_force, nbody_force_jvp, nbody_step, nbody_fix, nbody_set_encounter, &
        nbody_energy, nbody_step_jvp
    implicit none
    real(dp), parameter :: pi = 3.14159265358979323846_dp
    type(fortnum_status_t) :: status
    type(units_t) :: units
    logical :: ok

    units = reduced_units(1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp)
    call bspline_identities()
    call against_direct(160, 5.0_dp, 1.0e-4_dp)
    call against_direct(160, 5.0_dp, 1.0e-6_dp)
    call against_direct(400, 3.0_dp, 1.0e-5_dp)
    call gradient_and_tangent()
    call fixed_particle()
    call encounter_split()
    print '(a)', 'PASS: B-splines, PME vs direct Ewald, gradient, tangent, fixed, encounters'
contains
    subroutine check(condition, message)
        logical, intent(in) :: condition
        character(*), intent(in) :: message
        if (.not. condition) then
            print '(a)', 'FAIL: '//message
            error stop 1
        end if
    end subroutine

    subroutine bspline_identities()
        real(dp) :: t, w(8), d(8), s(8), wp(8), dp8(8), sp(8), h
        integer :: k
        h = 1.0e-6_dp
        do k = 0, 9
            t = (k + 0.37_dp)/10
            call pme_bspline4(t, w(1), w(2), w(3), w(4), d(1), d(2), d(3), d(4), &
                              s(1), s(2), s(3), s(4))
            call check(abs(sum(w(1:4)) - 1) < 1e-14_dp, 'order 4 partition of unity')
            call check(abs(sum(d(1:4))) < 1e-13_dp .and. abs(sum(s(1:4))) < 1e-12_dp, &
                       'order 4 derivative sums vanish')
            call pme_bspline8(t, w(1), w(2), w(3), w(4), w(5), w(6), w(7), w(8), &
                d(1), d(2), d(3), d(4), d(5), d(6), d(7), d(8), &
                s(1), s(2), s(3), s(4), s(5), s(6), s(7), s(8))
            call pme_bspline8(t + h, wp(1), wp(2), wp(3), wp(4), wp(5), wp(6), wp(7), &
                wp(8), dp8(1), dp8(2), dp8(3), dp8(4), dp8(5), dp8(6), dp8(7), dp8(8), &
                sp(1), sp(2), sp(3), sp(4), sp(5), sp(6), sp(7), sp(8))
            call check(abs(sum(w) - 1) < 1e-13_dp, 'order 8 partition of unity')
            ! Centre of mass of the weights: sum_j (t + j) M_p(t + j) = p/2.
            call check(abs(sum([(real(k - 1, dp) + t, k = 1, 8)]*w) - 4) < 1e-12_dp, &
                       'order 8 first moment p/2')
            call check(maxval(abs((wp - w)/h - d)) < 1e-5_dp, 'order 8 dw/dt')
            call check(maxval(abs((dp8 - d)/h - s)) < 1e-5_dp, 'order 8 d2w/dt2')
        end do
        call pme_bspline6(0.5_dp, w(1), w(2), w(3), w(4), w(5), w(6), d(1), d(2), &
            d(3), d(4), d(5), d(6), s(1), s(2), s(3), s(4), s(5), s(6))
        ! M_6 at half-integers: 1/3840, 237/3840, 1682/3840 (symmetric).
        call check(abs(w(1) - 1/3840.0_dp) < 1e-15_dp .and. &
                   abs(w(3) - 1682/3840.0_dp) < 1e-14_dp, 'order 6 known values')
    end subroutine

    subroutine random_state(n, length, seed, state)
        integer, intent(in) :: n
        real(dp), intent(in) :: length
        integer(int64), intent(in) :: seed
        real(dp), intent(out) :: state(:)
        type(rng_t) :: g
        integer :: i
        call rng_seed(g, seed, status)
        do i = 1, size(state)
            call rng_uniform(g, state(i))
        end do
        state(:3*n) = length*state(:3*n)
        state(3*n + 1:) = state(3*n + 1:) - 0.5_dp
    end subroutine

    !> n electrons (plus, for n=400, one fixed charge -5e) in a cube: PME
    !> forces and energy against the direct Ewald sum at tolerance 1e-12.
    subroutine against_direct(n, length, tolerance)
        integer, intent(in) :: n
        real(dp), intent(in) :: length, tolerance
        type(species_t), allocatable :: particles(:)
        type(ewald_parameters_t) :: ep
        type(pme_parameters_t) :: pp
        type(nbody_system_t) :: direct, mesh
        real(dp), allocatable :: state(:), fd(:), fm(:)
        real(dp) :: ud, um, a, f0, rms, momentum(3)
        integer :: i
        allocate (particles(n), state(6*n), fd(3*n), fm(3*n))
        particles = electron()
        if (n == 400) particles(n)%charge = -5
        call random_state(n, length, 11_int64*n, state)
        call ewald_select(length, 1.0e-12_dp, ep, ok, particles=n)
        call check(ok, 'direct Ewald parameters')
        call nbody_periodic(direct, units, particles, ep, ok)
        call check(ok, 'direct Ewald system')
        call pme_select(length, n, units%coulomb, tolerance, pp, ok)
        call check(ok, 'PME parameters from tolerance')
        call nbody_periodic_mesh(mesh, units, particles, pp, ok)
        call check(ok, 'PME system')
        call nbody_force(direct, state, fd, ud, status)
        call check(status%code == FORTNUM_OK, 'direct force')
        call nbody_force(mesh, state, fm, um, status)
        call check(status%code == FORTNUM_OK, 'PME force')
        a = (3*length**3/(4*pi*n))**(1.0_dp/3); f0 = units%coulomb/a**2
        rms = sqrt(sum((fm - fd)**2)/n)
        do i = 1, 3
            momentum(i) = sum(fm(i::3))
        end do
        print '(a,i5,a,es9.2,a,f6.2,a,i4,a,i2,a,f5.2,a,es9.2,a,es9.2,a,es9.2,a,f9.5)', &
            'N=', n, ' tol=', tolerance, ' rc/a=', pp%cutoff/a, ' K=', pp%mesh, &
            ' p=', pp%order, ' alpha*rc=', pp%alpha*pp%cutoff, ' rms/F0=', rms/f0, &
            ' est/F0=', sqrt(pp%real_estimate**2 + pp%reciprocal_estimate**2)/f0, &
            ' dU/(N c/a)=', (um - ud)/(n*units%coulomb/a), ' U/(N c/a)=', &
            ud/(n*units%coulomb/a)
        call check(rms <= tolerance*f0, 'measured rms force error within tolerance')
        call check(rms <= 1.5_dp*sqrt(pp%real_estimate**2 + pp%reciprocal_estimate**2) &
                   .and. rms >= 0.1_dp*pp%reciprocal_estimate, &
                   "error estimate within a factor 1.5 of the measurement")
        call check(abs(um - ud) <= 10*tolerance*f0*a*n, 'PME energy')
        call check(maxval(abs(momentum)) <= 10*tolerance*f0*sqrt(real(n, dp)), &
                   'net mesh force (momentum defect) at the force tolerance')
    end subroutine

    !> F = -grad U and dF = (dF/dx) dx by central differences; the cell-list
    !> path (L/rc >= 3) against the direct sum.
    subroutine gradient_and_tangent()
        integer, parameter :: n = 1500
        type(species_t) :: particles(n)
        type(pme_parameters_t) :: pp
        type(ewald_parameters_t) :: ep
        type(nbody_system_t) :: mesh, direct
        real(dp) :: state(6*n), x(3*n), dx(3*n), f(3*n), fp(3*n), fmn(3*n), df(3*n)
        real(dp) :: up, um, u, h, length, fd(3*n), a
        integer :: i, j
        length = 6.0_dp
        particles = electron()
        call random_state(n, length, 7_int64, state)
        call pme_select(length, n, units%coulomb, 1.0e-5_dp, pp, ok, order=6, &
                        max_cutoff=length/3)
        call check(ok .and. length/pp%cutoff >= 3, 'cell-list parameters')
        call nbody_periodic_mesh(mesh, units, particles, pp, ok)
        call check(ok, 'PME system (cell list)')
        call ewald_select(length, 1.0e-12_dp, ep, ok, particles=n)
        call nbody_periodic(direct, units, particles, ep, ok)
        call nbody_force(mesh, state, f, u, status)
        call nbody_force(direct, state, fd, up, status)
        a = (3*length**3/(4*pi*n))**(1.0_dp/3)
        call check(sqrt(sum((f - fd)**2)/n) <= 1.0e-5_dp*units%coulomb/a**2, &
                   'cell-list PME vs direct Ewald')
        x = state(:3*n); h = 1.0e-5_dp
        do j = 1, 6
            i = 1 + 47*(j - 1)
            x(i) = x(i) + h; call nbody_force(mesh, x, fp, up, status)
            x(i) = x(i) - 2*h; call nbody_force(mesh, x, fmn, um, status)
            x(i) = x(i) + h
            call check(abs(-(up - um)/(2*h) - f(i)) <= 1.0e-6_dp*maxval(abs(f)), &
                       'PME force is the gradient of the PME energy')
        end do
        do i = 1, 3*n
            dx(i) = sin(1.3_dp*i)
        end do
        call nbody_force_jvp(mesh, x, dx, df, status)
        call check(status%code == FORTNUM_OK, 'PME tangent')
        h = 1.0e-6_dp
        call nbody_force(mesh, x + h*dx, fp, up, status)
        call nbody_force(mesh, x - h*dx, fmn, um, status)
        call check(maxval(abs((fp - fmn)/(2*h) - df)) <= 1.0e-5_dp*maxval(abs(df)), &
                   'PME force tangent vs central difference')
    end subroutine

    subroutine fixed_particle()
        integer, parameter :: n = 40
        type(species_t) :: particles(n)
        type(pme_parameters_t) :: pp
        type(nbody_system_t) :: mesh
        real(dp) :: state(6*n), start(3)
        integer :: k
        particles = electron(); particles(1)%charge = -3
        call random_state(n, 4.0_dp, 3_int64, state)
        state(3*n + 1:3*n + 3) = 0; start = state(1:3)
        call pme_select(4.0_dp, n, units%coulomb, 1.0e-4_dp, pp, ok)
        call nbody_periodic_mesh(mesh, units, particles, pp, ok)
        call nbody_fix(mesh, 1, ok)
        call check(ok, 'fix particle')
        do k = 1, 20
            call nbody_step(mesh, state, 1.0e-3_dp, status)
            call check(status%code == FORTNUM_OK, 'mesh step')
        end do
        call check(all(state(1:3) == start) .and. all(state(3*n + 1:3*n + 3) == 0), &
                   'fixed particle stays at rest')
    end subroutine

    !> Weakly coupled cell (64 electrons, coupling 1e-3, thermal speeds ~1)
    !> with one large-angle encounter (b = 1e-4 < b90 = 5e-4): the
    !> encounter-resolved split step (dt = 0.02) against plain KDK of the same
    !> Hamiltonian at dt = 1e-5, its energy, and the split-step tangent by
    !> central differences.
    subroutine encounter_split()
        integer, parameter :: n = 64
        type(species_t) :: particles(n)
        type(pme_parameters_t) :: pp
        type(nbody_system_t) :: split, plain
        type(units_t) :: weak
        real(dp) :: base(6*n), a(6*n), b(6*n), t(6*n), p(6*n), m(6*n), e0, e1, h
        integer :: k
        weak = reduced_units(1.0_dp, 1.0e-3_dp, 1.0_dp, 1.0_dp)
        particles = electron()
        call random_state(n, 4.0_dp, 21_int64, base)
        base(3*n + 1:) = 2*base(3*n + 1:)
        base(1:3) = [1.0_dp, 2.0_dp, 2.0_dp]; base(4:6) = [1.1_dp, 2.0001_dp, 2.0_dp]
        base(3*n + 1:3*n + 3) = [1.0_dp, 0.0_dp, 0.0_dp]
        base(3*n + 4:3*n + 6) = [-1.0_dp, 0.0_dp, 0.0_dp]
        call pme_select(4.0_dp, n, weak%coulomb, 1.0e-6_dp, pp, ok)
        call nbody_periodic_mesh(split, weak, particles, pp, ok)
        call nbody_periodic_mesh(plain, weak, particles, pp, ok)
        call nbody_set_encounter(split, 0.02_dp, 2.0_dp, ok)
        call check(ok, 'encounter split setup')
        a = base; b = base
        call nbody_energy(split, a, e0, status)
        do k = 1, 5
            call nbody_step(split, a, 0.02_dp, status)
            call check(status%code == FORTNUM_OK, 'split step')
            if (k == 3) call check(split%near%deepest > 100, 'the close pair is substepped')
        end do
        do k = 1, 10000
            call nbody_step(plain, b, 1.0e-5_dp, status)
        end do
        call nbody_energy(split, a, e1, status)
        print '(a,es9.2,a,es9.2,a,i5)', ' split vs fine KDK max|dx|=', &
            maxval(abs(a(:3*n) - b(:3*n))), ' dE=', e1 - e0, ' deepest=', split%near%deepest
        call check(maxval(abs(a(:3*n) - b(:3*n))) < 1e-4_dp, &
                   'encounter-resolved step reproduces the fine-step trajectory')
        call check(abs(e1 - e0) < 1e-5_dp*sum(base(3*n + 1:)**2)/2, 'split energy')
        ! Tangent of three split steps in a random direction.
        do k = 1, 6*n
            t(k) = cos(0.7_dp*k)
        end do
        a = base
        do k = 1, 3
            call nbody_step_jvp(split, a, t, 0.02_dp, status)
        end do
        h = 1.0e-7_dp
        p = base + h*[(cos(0.7_dp*k), k = 1, 6*n)]
        m = base - h*[(cos(0.7_dp*k), k = 1, 6*n)]
        do k = 1, 3
            call nbody_step(split, p, 0.02_dp, status)
            call nbody_step(split, m, 0.02_dp, status)
        end do
        call check(maxval(abs((p - m)/(2*h) - t)) < 1e-5_dp*maxval(abs(t)), &
                   'split-step tangent vs central difference')
    end subroutine
end program test_pme
