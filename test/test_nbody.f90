program test_nbody
    ! Independent oracles for the shared microscopic foundation: CODATA/hand
    ! scale values, lattice Madelung constants (literature and an Evjen direct
    ! sum), Ewald truncation bounds, Newton/momentum/energy identities and
    ! finite-difference derivative checks.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_status, only: fortnum_status_t, FORTNUM_OK
    use kin6d_species, only: units_t, species_t, atomic_units, periodic_pair_units, &
        electron, pair_scales, pair_scales_t, plasma_scales, plasma_scales_t, &
        sommerfeld_parameter, pair_coupling
    use kin6d_coulomb_leaves, only: coulomb_partial_wave, coulomb_pair, &
        coulomb_pair_vjp, coulomb_pair_jvp
    use kin6d_ewald, only: ewald_parameters_t, ewald_select, ewald_truncation_bound
    use kin6d_nbody, only: nbody_system_t, nbody_isolated, nbody_periodic, &
        nbody_add_source, nbody_set_sphere, nbody_force, nbody_force_jvp, &
        nbody_energy, nbody_step, nbody_step_jvp, nbody_encounter_step
    implicit none
    real(dp), parameter :: pi = 3.14159265358979323846_dp
    ! Literature: Madelung constants of rock salt and of the simple-cubic
    ! point lattice in a neutralizing background (Ewald/Wigner constant).
    real(dp), parameter :: madelung_nacl = 1.747564594633_dp
    real(dp), parameter :: wigner_sc = -2.837297479480619_dp
    type(units_t) :: au
    type(fortnum_status_t) :: status
    logical :: ok

    au = atomic_units()
    call scales()
    call madelung()
    call ewald_bound()
    call isolated_identities()
    call periodic_identities()
    call derivatives()
    call encounters()
    print '(a)', 'PASS: scales, Madelung, Ewald bound, Newton, momentum, energy, AD'
contains
    subroutine check(condition, message)
        logical, intent(in) :: condition
        character(*), intent(in) :: message
        if (.not. condition) then
            print '(a)', 'FAIL: '//message
            error stop 1
        end if
    end subroutine

    subroutine scales()
        type(pair_scales_t) :: s
        type(plasma_scales_t) :: p
        type(units_t) :: u
        real(dp) :: v, fermi
        ! CODATA 2018 Bohr radius and atomic time unit.
        call check(abs(au%length_si/5.29177210903e-11_dp - 1) < 1e-9_dp, &
                   'atomic length unit is the Bohr radius')
        call check(abs(au%time_si/2.4188843265857e-17_dp - 1) < 1e-9_dp, &
                   'atomic time unit')
        u = periodic_pair_units()
        call check(abs(u%length_si/(2*5.29177210903e-11_dp) - 1) < 1e-9_dp, &
                   'pair units: e^2/(4 pi eps0)=2 gives length unit 2 a0')
        call check(abs(u%time_si/(4*2.4188843265857e-17_dp) - 1) < 1e-9_dp, &
                   'pair units time')
        call check(abs(pair_coupling(u, electron(), electron()) - 2) < 1e-15_dp, &
                   'pair coupling g=2')
        s = pair_scales(au, electron(), electron())
        call check(abs(s%coupling - 1) + abs(s%reduced_mass - 0.5_dp) + &
                   abs(s%bohr_radius - 2) + abs(s%rydberg - 0.25_dp) < 1e-15_dp, &
                   'electron pair scales: mu=1/2, a=2, Ry=1/4')
        call check(abs(sommerfeld_parameter(au, electron(), electron(), 0.5_dp) - 2) &
                   < 1e-15_dp, 'Sommerfeld eta = g/(hbar v)')
        p = plasma_scales(au, electron(), 3/(4*pi), 0.5_dp)
        fermi = (9*pi/4)**(2.0_dp/3)/2
        call check(abs(p%wigner_seitz - 1) < 1e-14_dp, 'Wigner-Seitz radius')
        call check(abs(p%gamma - 2) < 1e-14_dp, 'Gamma = g/(a T)')
        call check(abs(p%fermi - fermi) < 1e-14_dp, 'Fermi energy')
        call check(abs(p%theta - 0.5_dp/fermi) < 1e-14_dp, 'Theta')
        call check(abs(p%rs - 1) < 1e-14_dp, 'r_s')
        call check(abs(p%debye_length - sqrt(0.5_dp/3)) < 1e-14_dp, 'Debye length')
        call check(abs(p%plasma_frequency - sqrt(3.0_dp)) < 1e-14_dp, 'omega_p')
        call coulomb_partial_wave(2.0_dp, 1.0_dp, 1.0_dp, 0.5_dp, 1.0_dp, v)
        call check(abs(v - 1) < 1e-15_dp, 'partial-wave potential g/r + l(l+1)/(2mu r^2)')
    end subroutine

    subroutine madelung()
        type(ewald_parameters_t) :: e
        type(nbody_system_t) :: system
        type(species_t) :: ions(8)
        real(dp) :: state(48), force(24), u, evjen
        integer :: i
        call ewald_select(1.0_dp, 1e-13_dp, e, ok)
        call check(ok, 'Ewald parameter selection')
        call nbody_periodic(system, au, [species_t('ion', 1.0_dp, 1.0_dp)], e, ok)
        call check(ok, 'single charge periodic system')
        state = 0
        call nbody_force(system, state(:6), force(:3), u, status)
        call check(status%code == FORTNUM_OK, 'single charge force')
        call check(abs(u - wigner_sc/2) < 1e-11_dp, 'simple-cubic Wigner constant')
        call check(maxval(abs(force(:3))) < 1e-12_dp, 'lattice symmetry: no self force')
        ! Rock salt conventional cell, a=1: four cations and four anions.
        ions(:4) = species_t('cation', 1.0_dp, 1.0_dp)
        ions(5:) = species_t('anion', -1.0_dp, 1.0_dp)
        state = 0
        state(:24) = [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, .5_dp, .5_dp, .5_dp, 0.0_dp, &
                      .5_dp, .5_dp, .5_dp, 0.0_dp, .5_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                      .5_dp, 0.0_dp, 0.0_dp, 0.0_dp, .5_dp, .5_dp, .5_dp, .5_dp]
        call nbody_periodic(system, au, ions, e, ok)
        call nbody_force(system, state, force, u, status)
        call check(status%code == FORTNUM_OK, 'rock salt force')
        evjen = evjen_madelung(24)
        call check(abs(evjen - madelung_nacl) < 1e-6_dp, 'Evjen direct sum oracle')
        call check(abs(u + 8*madelung_nacl) < 1e-10_dp, 'rock salt Madelung energy')
        do i = 1, 8
            call check(maxval(abs(force(3*i - 2:3*i))) < 1e-10_dp, &
                       'rock salt equilibrium forces vanish')
        end do
    end subroutine

    !> Evjen's charge-neutral cube sum of the rock-salt Madelung constant.
    pure real(dp) function evjen_madelung(r) result(m)
        integer, intent(in) :: r
        integer :: i, j, k, faces
        real(dp) :: w
        m = 0
        do i = -r, r
            do j = -r, r
                do k = -r, r
                    if (i == 0 .and. j == 0 .and. k == 0) cycle
                    faces = merge(1, 0, abs(i) == r) + merge(1, 0, abs(j) == r) + &
                            merge(1, 0, abs(k) == r)
                    w = 0.5_dp**faces
                    m = m - w*(-1)**(i + j + k)/sqrt(real(i*i + j*j + k*k, dp))
                end do
            end do
        end do
    end function

    !> Five charges of mixed sign; deterministic irregular positions.
    subroutine sample_state(state, length)
        real(dp), intent(out) :: state(30)
        real(dp), intent(in) :: length
        integer :: j
        do j = 1, 30
            state(j) = length*(modulo(0.6180339887_dp*j*j + 0.31_dp*j, 1.0_dp) - .5_dp)
        end do
        state(16:) = state(16:)/length
    end subroutine

    function mixed() result(s)
        type(species_t) :: s(5)
        s = [species_t('e', -1.0_dp, 1.0_dp), species_t('e', -1.0_dp, 1.0_dp), &
             species_t('i', 2.0_dp, 3.0_dp), species_t('e', -1.0_dp, 1.0_dp), &
             species_t('p', 1.0_dp, 2.0_dp)]
    end function

    subroutine ewald_bound()
        type(ewald_parameters_t) :: loose, tight
        type(nbody_system_t) :: a, b
        type(species_t) :: s(5)
        real(dp) :: state(30), fa(15), fb(15), ua, ub, gsum, gpair, l
        integer :: i, j
        l = 2.5_dp; s = mixed()
        call ewald_select(l, 1e-5_dp, loose, ok)
        call check(ok, 'loose Ewald selection')
        call ewald_select(l, 1e-14_dp, tight, ok)
        tight%alpha = 0.9_dp*tight%alpha; tight%real_cutoff = tight%real_cutoff + 1
        tight%reciprocal_cutoff = tight%reciprocal_cutoff + 4
        call ewald_truncation_bound(tight)
        call check(tight%potential_bound*l < 1e-15_dp, 'reference truncation')
        call sample_state(state, l)
        call nbody_periodic(a, au, s, loose, ok)
        call nbody_periodic(b, au, s, tight, ok)
        call nbody_force(a, state, fa, ua, status)
        call nbody_force(b, state, fb, ub, status)
        gsum = 0
        do i = 1, 5
            gsum = gsum + s(i)%charge**2/2
            gpair = 0
            do j = 1, 5
                if (j /= i) gpair = gpair + abs(s(i)%charge*s(j)%charge)
                if (j > i) gsum = gsum + abs(s(i)%charge*s(j)%charge)
            end do
            call check(maxval(abs(fa(3*i - 2:3*i) - fb(3*i - 2:3*i))) <= &
                       gpair*(loose%force_bound + tight%force_bound) + 1e-12_dp, &
                       'force truncation within Ewald bound')
        end do
        call check(abs(ua - ub) <= gsum*(loose%potential_bound + &
                   tight%potential_bound) + 1e-12_dp, 'energy within Ewald bound')
        call check(abs(ua - ub) > 1e-13_dp, 'comparison resolves the loose truncation')
        call check(loose%potential_bound*l <= 1e-5_dp, 'requested tolerance met')
    end subroutine

    subroutine isolated_identities()
        type(nbody_system_t) :: system
        real(dp) :: state(30), force(15), u, e0, e1, p0(3), p1(3), f2(6), u2
        integer :: k
        call nbody_isolated(system, au, mixed(), ok)
        call sample_state(state, 4.0_dp)
        call nbody_force(system, state, force, u, status)
        call check(status%code == FORTNUM_OK, 'isolated force')
        call check(maxval(abs(sum(reshape(force, [3, 5]), 2))) < 1e-13_dp, &
                   'Newton third law: total internal force vanishes')
        call nbody_isolated(system, au, [electron(), electron()], ok)
        call nbody_force(system, [0.3_dp, -.4_dp, 1.2_dp, -.7_dp, .2_dp, .1_dp], &
                         f2, u2, status)
        call check(maxval(abs(f2(:3) + f2(4:))) < 1e-15_dp, 'pair antisymmetry')
        call check(abs(u2 - 1/sqrt(1.0_dp + .36_dp + 1.21_dp)) < 1e-15_dp, &
                   'pair potential g/r')
        call nbody_isolated(system, au, mixed(), ok)
        call nbody_energy(system, state, e0, status)
        p0 = momentum(state)
        do k = 1, 2000
            call nbody_step(system, state, 1e-4_dp, status)
            call check(status%code == FORTNUM_OK, 'isolated step')
        end do
        call nbody_energy(system, state, e1, status)
        p1 = momentum(state)
        call check(maxval(abs(p1 - p0)) < 1e-12_dp, 'N=5 momentum conservation')
        call check(abs(e1 - e0) < 1e-6_dp*abs(e0), 'N=5 energy conservation')
        call gradient_identity(system, state)
    end subroutine

    function momentum(state) result(p)
        real(dp), intent(in) :: state(30)
        real(dp) :: p(3), m(5)
        integer :: i
        m = [1.0_dp, 1.0_dp, 3.0_dp, 1.0_dp, 2.0_dp]; p = 0
        do i = 1, 5
            p = p + m(i)*state(15 + 3*i - 2:15 + 3*i)
        end do
    end function

    !> F = -grad U by central differences.
    subroutine gradient_identity(system, state)
        type(nbody_system_t), intent(inout) :: system
        real(dp), intent(in) :: state(30)
        real(dp) :: x(30), force(15), f(15), u, up, um, h
        integer :: j
        h = 1e-5_dp
        call nbody_force(system, state, force, u, status)
        do j = 1, 15
            x = state; x(j) = x(j) + h
            call nbody_force(system, x, f, up, status)
            x(j) = x(j) - 2*h
            call nbody_force(system, x, f, um, status)
            call check(abs((up - um)/(2*h) + force(j)) < 1e-6_dp*(1 + abs(force(j))), &
                       'energy-gradient identity')
        end do
    end subroutine

    subroutine periodic_identities()
        type(ewald_parameters_t) :: e
        type(nbody_system_t) :: system
        real(dp) :: state(30), force(15), u, e0, e1, p0(3), p1(3)
        integer :: k
        call ewald_select(2.5_dp, 1e-12_dp, e, ok)
        call nbody_periodic(system, au, mixed(), e, ok)
        call check(ok, 'periodic system')
        call sample_state(state, 2.5_dp)
        call nbody_force(system, state, force, u, status)
        call check(maxval(abs(sum(reshape(force, [3, 5]), 2))) < 1e-11_dp, &
                   'periodic total force vanishes')
        call gradient_identity(system, state)
        call nbody_energy(system, state, e0, status)
        p0 = momentum(state)
        do k = 1, 500
            call nbody_step(system, state, 2e-4_dp, status)
            call check(status%code == FORTNUM_OK, 'periodic step')
        end do
        call nbody_energy(system, state, e1, status)
        p1 = momentum(state)
        call check(maxval(abs(p1 - p0)) < 1e-10_dp, 'periodic momentum conservation')
        call check(abs(e1 - e0) < 1e-6_dp*abs(e0), 'periodic energy conservation')
    end subroutine

    subroutine derivatives()
        type(ewald_parameters_t) :: e
        type(nbody_system_t) :: iso, per
        real(dp) :: state(30), dx(30), df(15), fp(15), fm(15), u, h
        real(dp) :: sp(30), sm(30), tangent(30), y(30), fd(30), d(4), b(4), gr(4)
        integer :: k
        call nbody_isolated(iso, au, mixed(), ok)
        call nbody_add_source(iso, [0.2_dp, -0.1_dp, 0.05_dp], -2.0_dp, ok)
        call check(ok, 'fixed source')
        call nbody_set_sphere(iso, 3.0_dp, 1.5_dp, ok)
        call check(ok, 'fixed sphere')
        call ewald_select(2.5_dp, 1e-12_dp, e, ok)
        call nbody_periodic(per, au, mixed(), e, ok)
        call sample_state(state, 2.5_dp)
        call sample_state(dx, 1.0_dp)
        h = 1e-6_dp
        call nbody_force_jvp(iso, state(:15), dx(:15), df, status)
        call check(status%code == FORTNUM_OK, 'isolated force JVP')
        call nbody_force(iso, state + h*dx*[(1, k=1, 15), (0, k=1, 15)], fp, u, status)
        call nbody_force(iso, state - h*dx*[(1, k=1, 15), (0, k=1, 15)], fm, u, status)
        call check(maxval(abs(df - (fp - fm)/(2*h))) < 1e-6_dp*maxval(abs(df)), &
                   'isolated force JVP vs finite difference')
        call nbody_force_jvp(per, state(:15), dx(:15), df, status)
        call nbody_force(per, state + h*dx*[(1, k=1, 15), (0, k=1, 15)], fp, u, status)
        call nbody_force(per, state - h*dx*[(1, k=1, 15), (0, k=1, 15)], fm, u, status)
        call check(maxval(abs(df - (fp - fm)/(2*h))) < 1e-6_dp*maxval(abs(df)), &
                   'periodic force JVP vs finite difference')
        ! Tangent of ten KDK steps vs finite difference of the discrete flow.
        y = state; tangent = dx
        do k = 1, 10
            call nbody_step_jvp(iso, y, tangent, 1e-3_dp, status)
            call check(status%code == FORTNUM_OK, 'step tangent')
        end do
        sp = state + h*dx; sm = state - h*dx
        do k = 1, 10
            call nbody_step(iso, sp, 1e-3_dp, status)
            call nbody_step(iso, sm, 1e-3_dp, status)
        end do
        fd = (sp - sm)/(2*h)
        call check(maxval(abs(fd - tangent)) < 1e-6_dp*maxval(abs(tangent)), &
                   'KDK step tangent vs finite difference')
        call check(maxval(abs(y - sp)) < 1e-4_dp, 'tangent run advances state')
        ! Leaf JVP/VJP Euclidean pairing of the single Coulomb definition.
        d = [0.2_dp, -0.3_dp, 0.1_dp, 0.4_dp]; b = [0.5_dp, 0.1_dp, -0.2_dp, 0.3_dp]
        call coulomb_pair_jvp(0.7_dp, -1.1_dp, 0.3_dp, 1.3_dp, d(1), d(2), d(3), d(4), &
                              fp(1), fp(2), fp(3), fp(4))
        call coulomb_pair_vjp(0.7_dp, -1.1_dp, 0.3_dp, 1.3_dp, b(1), b(2), b(3), b(4), &
                              gr(1), gr(2), gr(3), gr(4))
        call check(abs(dot_product(b, fp(:4)) - dot_product(d, gr)) < 1e-14_dp, &
                   'Coulomb JVP/VJP pairing')
    end subroutine

    subroutine encounters()
        type(nbody_system_t) :: system
        real(dp) :: dt, force(6), u
        call nbody_isolated(system, au, [electron(), electron()], ok)
        call nbody_encounter_step(system, [0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, &
            0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], 0.1_dp, dt, status)
        call check(status%code == FORTNUM_OK, 'encounter step')
        ! r=1, at rest, g=1, mu=1/2: Kepler time sqrt(r^3 mu/g).
        call check(abs(dt - 0.1_dp*sqrt(0.5_dp)) < 1e-15_dp, 'Kepler encounter time')
        call nbody_encounter_step(system, [0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, &
            0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, -100.0_dp, 0.0_dp, 0.0_dp], 0.1_dp, dt, &
            status)
        call check(abs(dt - 0.1_dp/100) < 1e-15_dp, 'crossing time controls fast pair')
        call nbody_encounter_step(system, [1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, &
            0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], 0.1_dp, dt, status)
        call check(status%code /= FORTNUM_OK, 'coincident encounter fails closed')
        call nbody_force(system, [1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp], &
                         force, u, status)
        call check(status%code /= FORTNUM_OK, 'coincident force fails closed')
        call coulomb_pair(3.0_dp, -4.0_dp, 12.0_dp, 2.0_dp, force(1), force(2), &
                          force(3), u)
        call check(abs(u - 2.0_dp/13) < 1e-15_dp .and. &
                   abs(force(3) - 24.0_dp/2197) < 1e-17_dp, 'Coulomb hand value')
    end subroutine
end program test_nbody
