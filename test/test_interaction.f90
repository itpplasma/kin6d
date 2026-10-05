program test_interaction
    use, intrinsic :: iso_fortran_env, only: dp => real64, qp => real128
    use fortnum_interval, only: interval_t, ipoint
    use kin6d_interaction_laboratory
    use kin6d_interaction_leaves
    use kin6d_coulomb_leaves, only: coulomb_pair, coulomb_pair_interval
    use kin6d_species, only: units_t, species_t, atomic_units, electron
    use kin6d_nbody, only: nbody_system_t, nbody_isolated, nbody_add_source, &
                           nbody_step, nbody_step_jvp
    use fortnum_status, only: fortnum_status_t, FORTNUM_OK
    implicit none
    type(units_t) :: au
    type(nbody_system_t) :: relative
    type(fortnum_status_t) :: status
    type(interaction_state_t) :: state
    type(interval_t) :: rb(3), ub(3), fi(3), pi
    real(dp) :: r0(3) = [-6._dp, .6_dp, .2_dp], u0(3) = [2._dp, .1_dp, -.15_dp]
    real(dp) :: r(3), u(3), f(3), fn(3), p, pn, err(3), initial_energy, energy
    real(dp) :: x(6), direction(6), j(6), yp(6), ym(6), fd(6), eps
    real(qp) :: exact_r(3), exact_u(3), force(3), norm, physical_time
    logical :: ok, bounded
    integer :: k, n, level
    au = atomic_units()
    call coulomb_pair(3._dp, -4._dp, 2._dp, 1.7_dp, f(1), f(2), f(3), p)
    call coulomb_pair(-3._dp, 4._dp, -2._dp, 1.7_dp, fn(1), fn(2), fn(3), pn)
    call check(maxval(abs(f + fn)) < 1.e-15_dp, 'pair force antisymmetry')
    call check(abs(p - pn) < 1.e-15_dp, 'pair potential exchange symmetry')
    norm = sqrt(29._qp); force = real(1.7_dp, qp)*[3._qp, -4._qp, 2._qp]/norm**3
    call check(maxval(abs(real(f, qp) - force)) < 2.e-16_qp, &
               'Cartesian Coulomb oracle')
    call coulomb_pair_interval(ipoint(3._dp), ipoint(-4._dp), ipoint(2._dp), &
                               ipoint(1.7_dp), fi(1), fi(2), fi(3), pi)
    call check(all(real(fi%lo, qp) <= force) .and. all(force <= real(fi%hi, qp)), &
               'independent force interval containment')
    call interaction_particle_coordinate_interval(ipoint(nearest(0._dp, &
                                                            1._dp)), ipoint(-1._dp), pi)
    call check(real(pi%lo, qp) <= -real(nearest(0._dp, 1._dp), qp)/2 .and. &
               -real(nearest(0._dp, 1._dp), qp)/2 <= real(pi%hi, qp), &
               'subnormal physical coordinate containment')
    ! Relative problem: reduced mass 1/2, charge -1, fixed partner charge -1.
    call nbody_isolated(relative, au, [species_t('relative', -1._dp, .5_dp)], ok)
    call nbody_add_source(relative, [0._dp, 0._dp, 0._dp], -1._dp, ok)
    call check(ok, 'relative system')
    x = [r0, u0]
    direction = [.2_dp, -.1_dp, .3_dp, -.4_dp, .2_dp, .1_dp]
    yp = x; j = direction
    call nbody_step_jvp(relative, yp, j, .013_dp, status)
    call check(status%code == FORTNUM_OK, 'KDK tangent step')
    eps = 2.e-5_dp
    yp = x + eps*direction; ym = x - eps*direction
    call nbody_step(relative, yp, .013_dp, status)
    call nbody_step(relative, ym, .013_dp, status)
    fd = (yp - ym)/(2*eps)
    call check(maxval(abs(fd - j)) < 3.e-10_dp, &
               'KDK step independent directional finite difference')
    call orbit(real(r0, qp), real(u0, qp), 2._qp, 6._qp, exact_r, exact_u)
    do level = 1, 3
        n = 1500*2**level
        call init_interaction(state, r0, u0, au, electron(), electron(), ok)
        call check(ok, 'initialize arbitrary true 3D scattering')
        do k = 1, n
            call step_interaction(state, 6._dp/real(n, dp), ok, certify=.false.)
            call check(ok, 'candidate step')
        end do
        call sample_interaction(state, r, u, rb, ub, ok, bounded)
        call check(.not. bounded, &
                   'uncertified advance cannot expose a certified packet')
        err(level) = real(maxval(abs(real(r, qp) - exact_r)), dp)
    end do
    call check(err(1)/err(2) > 3.8_dp .and. err(2)/err(3) > 3.8_dp, &
               'Rutherford orbit second order convergence')
    call check(err(3) < 2.e-6_dp, 'finite time Rutherford orbit accuracy')
    call init_interaction(state, r0, u0, au, electron(), electron(), ok)
    call interaction_energy(r0(1), r0(2), r0(3), u0(1), u0(2), u0(3), &
                            interaction_coupling(state), initial_energy)
    do k = 1, 12000
        call step_interaction(state, .0005_dp, ok)
        call check(ok, 'validated candidate step')
        call sample_interaction(state, r, u, rb, ub, ok, bounded)
        call check(bounded, 'well separated Picard continuation closes')
        if (mod(k, 100) /= 0) cycle
        physical_time = real(k, qp)*real(.0005_dp, qp)
        call orbit(real(r0, qp), real(u0, qp), 2._qp, physical_time, &
                   exact_r, exact_u)
        call check(all(real(rb%lo, &
                            qp) <= exact_r) .and. all(exact_r <= real(rb%hi, qp)), &
                   'Rutherford true position inside rigorous propagated box')
        call check(all(real(ub%lo, &
                            qp) <= exact_u) .and. all(exact_u <= real(ub%hi, qp)), &
                   'Rutherford true velocity inside rigorous propagated box')
    end do
    call interaction_energy(r(1), r(2), r(3), u(1), u(2), u(3), &
                            interaction_coupling(state), energy)
    call check(abs(energy - initial_energy) < 1.e-7_dp, 'energy diagnostic')
    call neutral_contact_trajectory([-6._dp, .6_dp, 0._dp], [2._dp, 0._dp, &
                                                         0._dp], 1._dp, 2._dp, r, u, ok)
    call check(ok .and. maxval(abs(r - [-2._dp, .6_dp, 0._dp])) < 1.e-14_dp, &
               'neutral free approach')
    call neutral_contact_trajectory([-6._dp, .6_dp, 0._dp], [2._dp, 0._dp, &
                                                         0._dp], 1._dp, 4._dp, r, u, ok)
    call check(ok .and. maxval(abs(u - [-.56_dp, 1.92_dp, &
                                0._dp])) < 2.e-14_dp, 'neutral specular contact oracle')
    call check(abs(sum(u*u) - 4._dp) < 2.e-14_dp, 'elastic neutral contact')
    call init_interaction(state, [0._dp, 0._dp, 0._dp], u0, au, electron(), &
                          electron(), ok)
    call check(.not. ok, 'singular initial state rejected')
    call init_interaction(state, r0, u0, au, electron(), &
                          species_t('proton', 1._dp, 1836.15267343_dp), ok)
    call check(.not. ok, 'nonrepulsive coupling rejected by this named model')
    call init_interaction(state, r0, u0, au, electron(), electron(), ok)
    call step_interaction(state, 0._dp, ok)
    call check(.not. ok, 'invalid step rejected')
    call step_interaction(state, 1.e5_dp, ok)
    call sample_interaction(state, r, u, rb, ub, ok, bounded)
    call check(.not. bounded, 'failed tube never reports a certified continuation')
    print '(a,3es12.3)', 'Rutherford position errors: ', err
    print '(a)', &
        'PASS: Coulomb oracle, convergent 3D Rutherford orbit, rigorous box'// &
        'es, AD and neutral contact'
contains
    subroutine check(condition, message)
        logical, intent(in)::condition
        character(*), intent(in)::message
        if (.not. condition) then
            print '(a)', 'FAIL: '//message
            error stop 1
        end if
    end subroutine
    pure function cross(a, b) result(c)
        real(qp), intent(in)::a(3), b(3)
        real(qp)::c(3)
        c = [a(2)*b(3) - a(3)*b(2), a(3)*b(1) - a(1)*b(3), a(1)*b(2) - a(2)*b(1)]
    end function
    subroutine orbit(r0, u0, c, t, r, u)
        ! Test-only exact repulsive Kepler hyperbola, independent of the generator
        ! and candidate. Exact binary64 initial parameters are promoted to real128.
        real(qp), intent(in)::r0(3), u0(3), c, t
        real(qp), intent(out)::r(3), u(3)
        real(qp)::energy, l(3), ln, a, e, p(3), q(3), h0, h, lo, hi, target, tau, &
                   s, den, angle
        integer::i
        energy = sum(u0*u0)/2 + c/sqrt(sum(r0*r0))
        l = cross(r0, u0); ln = sqrt(sum(l*l)); a = c/(2*energy)
        e = sqrt(1 + 2*energy*ln**2/c**2)
        p = (cross(u0, l) + c*r0/sqrt(sum(r0*r0)))/(c*e)
        q = cross(l/ln, p); s = sqrt(e**2 - 1)
        h0 = asinh(dot_product(r0, q)/(a*s))
        tau = sqrt(a**3/c); target = t + tau*(e*sinh(h0) + h0)
        lo = -30; hi = 30
        do i = 1, 200
            h = (lo + hi)/2
            if (tau*(e*sinh(h) + h) < target) then
                lo = h
            else
                hi = h
            end if
        end do
        h = (lo + hi)/2; den = tau*(e*cosh(h) + 1)
        r = a*(e + cosh(h))*p + a*s*sinh(h)*q
        u = (a*sinh(h)*p + a*s*cosh(h)*q)/den
        angle = 2*atan(c/(ln*sqrt(2*energy)))
        call check(abs(angle - 2*asin(1/e)) < 1.e-30_qp, &
                   'Rutherford impact/deflection identity')
    end subroutine
end program
