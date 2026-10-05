program test_transport
    use, intrinsic :: iso_fortran_env, only: dp => real64, qp => real128
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use fortnum_interval, only: interval_t, interval, operator(+), operator(/)
    use kin6d_transport_laboratory, only: transport_state_t, initialize_transport, &
        transport_valid, solve_transport, solve_transport_tangent, &
        solve_transport_vjp, sample_transport, sample_transport_arrays, &
        sample_transport_coefficients_interval
    use kin6d_transport_leaves, only: transport_node, transport_node_jvp, &
        transport_node_vjp, transport_source_affine, transport_residual_cell_interval, &
        transport_cosine_moment_cell_interval
    use kin6d_transport_certificate, only: transport_certificate_t, certify_transport, &
        transport_certificate_accepts, certify_transport_tangent, &
        certify_transport_sample, certify_transport_tangent_sample, &
        transport_observable_band
    implicit none

    type(transport_state_t) :: coarse, fine, invalid, empty
    type(transport_certificate_t) :: cc, cf, ct, audit, tangent_audit
    real(dp), allocatable :: dsc(:, :), dss(:, :), da(:, :), db(:, :)
    real(dp), allocatable :: ba(:, :), bb(:, :), bsc(:, :), bss(:, :)
    real(dp) :: a, b, ta, tb, gv, gnu, gk, gsc, gss, value, tangent, lhs, rhs
    real(dp) :: x, velocity, radius, lo, hi, tangent_radius
    real(qp) :: reference, derivative
    type(interval_t) :: range, bc0, bc1, bc2, bs0, bs1, bs2, goal, tangent_goal, term
    logical :: ok
    integer :: j, k

    ! Independently inverted exact rational matrix [[2,1],[-1,2]].
    call transport_node(0.5_dp, 2.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, a, b)
    if (abs(a - 0.4_dp) > 1.0e-14_dp .or. abs(b - 2.2_dp) > 1.0e-14_dp) &
        error stop 'transport rational node fixture'
    ! L*dU=[1,-2]-.5*[.4,2.2]=[.8,-3.1].
    call transport_node_jvp(0.5_dp, 2.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, &
                            0.0_dp, 0.5_dp, 0.0_dp, 1.0_dp, -2.0_dp, ta, tb)
    if (abs(ta - 0.94_dp) > 1.0e-14_dp .or. abs(tb + 1.08_dp) > 1.0e-14_dp) &
        error stop 'transport independent implicit tangent fixture'
    call transport_node_vjp(0.5_dp, 2.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, &
                            2.0_dp, -1.0_dp, gv, gnu, gk, gsc, gss)
    if (max(abs(gv + 3.28_dp), abs(gnu - 1.52_dp), abs(gk + 0.82_dp), &
            abs(gsc - 0.6_dp), abs(gss + 0.8_dp)) > 1.0e-14_dp) &
        error stop 'transport rational coefficient VJP fixture'

    ! Arbitrary unsolved data yield rc(t)=8-33t+24t^2 and
    ! rs(t)=-10+16t+12t^2 for v=-1+2t. Independent Bernstein fixtures:
    ! bc=[8,-8.5,-1], bs=[-10,-2,18], rc(.5)=-2.5, rs(.5)=1.
    call transport_residual_cell_interval(interval(-1.0_dp), interval(1.0_dp), &
               interval(2.0_dp), interval(3.0_dp), interval(1.0_dp), interval(3.0_dp), &
             interval(2.0_dp), interval(-2.0_dp), interval(4.0_dp), interval(-1.0_dp), &
                      interval(-3.0_dp), interval(5.0_dp), bc0, bc1, bc2, bs0, bs1, bs2)
    call contains_fixture(bc0, 8.0_dp)
    call contains_fixture(bc1, -8.5_dp)
    call contains_fixture(bc2, -1.0_dp)
    call contains_fixture(bs0, -10.0_dp)
    call contains_fixture(bs1, -2.0_dp)
    call contains_fixture(bs2, 18.0_dp)
    range = bc0/4.0_dp + bc1/2.0_dp + bc2/4.0_dp
    call contains_fixture(range, -2.5_dp)
    range = bs0/4.0_dp + bs1/2.0_dp + bs2/4.0_dp
    call contains_fixture(range, 1.0_dp)

    call supplied_case(coarse, 16)
    call supplied_case(fine, 64)
    call certify_transport(coarse, cc, ok)
    if (.not. ok) error stop 'coarse continuous transport certificate failed'
    call certify_transport(fine, cf, ok)
    if (.not. ok) error stop 'fine continuous transport certificate failed'
    if (cf%state_radius >= 0.1_dp*cc%state_radius) &
        error stop 'whole-cell transport bound failed refinement'
    if (cf%state_radius <= 0.0_dp) error stop 'finite-grid transport radius is zero'
    if (.not. transport_certificate_accepts(fine, cf%state_radius)) &
        error stop 'valid transport radius rejected'
    if (transport_certificate_accepts(fine, 0.0_dp)) error stop 'zero radius accepted'
    allocate (dsc(0:2, 65), dss(0:2, 65), da(0:2, 65), db(0:2, 65))
    allocate (ba(0:2, 65), bb(0:2, 65), bsc(0:2, 65), bss(0:2, 65))
    dsc = 0.0_dp
    dss = 0.0_dp
    call solve_transport_tangent(fine, 1.0_dp, dsc, dss, da, db, ok)
    if (.not. ok) error stop 'implicit nu tangent solve failed'
    call certify_transport_tangent(fine, 1.0_dp, dsc, dss, da, db, ct, ok)
    if (.not. ok) error stop 'continuous nu tangent certificate failed'
    call certify_transport(fine, audit, ok, strict=.true.)
    if (.not. ok) error stop 'strict transport audit failed'
    call certify_transport_tangent(fine, 1.0_dp, dsc, dss, da, db, &
                                   tangent_audit, ok, strict=.true.)
    if (.not. ok) error stop 'strict nu tangent audit failed'
    ! Oracle exists only in TEST: independent binary128 mode resolvents for
    ! the supplied source, evaluated between nodes and at nonspecial x values.
    do j = 1, 31
        x = real(j, dp)/13.0_dp
        velocity = -0.95_dp + real(j, dp)/17.0_dp
        call test_only_resolvent(real(x, qp), real(velocity, qp), &
                                 real(fine%nu, qp), reference, derivative)
        call sample_transport(fine, x, velocity, value, ok)
        if (.not. ok) error stop 'transport sample failed'
        if (abs(real(value, qp) - reference) > real(cf%state_radius, qp)) &
            error stop 'continuous state radius misses independent mode oracle'
        if (abs(real(value, qp) - reference) > real(audit%state_radius, qp)) &
            error stop 'strict state radius misses independent mode oracle'
        call sample_transport_arrays(fine%v, da, db, x, velocity, tangent, ok)
        if (.not. ok) error stop 'tangent sample failed'
        if (abs(real(tangent, qp) - derivative) > real(ct%tangent_radius, qp)) &
            error stop 'continuous tangent radius misses independent mode oracle'
        if (abs(real(tangent, qp) - derivative) > &
            real(tangent_audit%tangent_radius, qp)) &
            error stop 'strict tangent radius misses independent mode oracle'
    end do
    call sample_transport_coefficients_interval(fine, fine%cosine, fine%sine, &
                       interval(0.2_dp, 0.3_dp), interval(-0.02_dp, 0.02_dp), range, ok)
    if (.not. ok) error stop 'interval reconstruction across velocity cells failed'
    do j = 1, 9
        x = 0.2_dp + real(j, dp)/100.0_dp
        velocity = -0.02_dp + real(j, dp)/250.0_dp
        call sample_transport(fine, x, velocity, value, ok)
        if (value < range%lo .or. value > range%hi) &
            error stop 'interval reconstruction misses interior candidate'
        call test_only_resolvent(real(x, qp), real(velocity, qp), &
                                 real(fine%nu, qp), reference, derivative)
        call certify_transport_sample(fine, cf, interval(0.2_dp, 0.3_dp), &
                                      interval(-0.02_dp, 0.02_dp), lo, hi, ok)
        if (.not. ok) error stop 'certified state cell sample failed'
        if (reference < real(lo, qp) .or. reference > real(hi, qp)) &
            error stop 'state cell band misses independent continuous oracle'
        call certify_transport_tangent_sample(fine, ct, da, db, &
                      interval(0.2_dp, 0.3_dp), interval(-0.02_dp, 0.02_dp), lo, hi, ok)
        if (.not. ok) error stop 'certified tangent cell sample failed'
        if (derivative < real(lo, qp) .or. derivative > real(hi, qp)) &
            error stop 'tangent cell band misses independent continuous oracle'
    end do

    ! Q=(1/(2*pi))*integral cos(x)*f dx dv. Source cos1=(1-v)/4
    ! gives the independent identities Q=atan(1/nu)/4 and dQ=-1/(4*(nu^2+1)).
    goal = interval(0.0_dp)
    tangent_goal = interval(0.0_dp)
    do j = 1, size(fine%v) - 1
        call transport_cosine_moment_cell_interval(interval(fine%v(j)), &
                                 interval(fine%v(j + 1)), interval(fine%cosine(1, j)), &
                                                  interval(fine%cosine(1, j + 1)), term)
        goal = goal + term
        call transport_cosine_moment_cell_interval(interval(fine%v(j)), &
              interval(fine%v(j + 1)), interval(da(1, j)), interval(da(1, j + 1)), term)
        tangent_goal = tangent_goal + term
    end do
    call transport_observable_band(goal, interval(2.0_dp), cf%state_radius, lo, hi, ok)
    if (.not. ok) error stop 'observable enclosure failed'
    reference = atan(1.0_qp/real(fine%nu, qp))/4.0_qp
    if (reference < real(lo, qp) .or. reference > real(hi, qp)) &
        error stop 'observable band misses independent analytic integral'
    call transport_observable_band(tangent_goal, interval(2.0_dp), &
                                   ct%tangent_radius, lo, hi, ok)
    if (.not. ok) error stop 'observable tangent enclosure failed'
    derivative = -1.0_qp/(4.0_qp*(real(fine%nu, qp)**2 + 1.0_qp))
    if (derivative < real(lo, qp) .or. derivative > real(hi, qp)) &
        error stop 'observable tangent band misses independent analytic integral'
    tangent_radius = ct%tangent_radius
    da(1, 33) = da(1, 33) + 1.0_dp
    call certify_transport_tangent(fine, 1.0_dp, dsc, dss, da, db, tangent_audit, ok)
    if (.not. ok) error stop 'corrupt tangent should receive a wider bound'
    if (tangent_audit%tangent_radius <= tangent_radius) &
        error stop 'corrupt tangent did not increase certified uncertainty'

    ! Complete caller-level transpose identity on raw nodal coefficient duals.
    ! This is separate from the continuous stability/tangent enclosure.
    do j = 1, 65
        do k = 0, 2
            dsc(k, j) = real(mod(j + k, 7), dp)/16.0_dp
            dss(k, j) = real(mod(j + 2*k, 5), dp)/32.0_dp
            ba(k, j) = real(mod(j + k, 3), dp)/8.0_dp
            bb(k, j) = real(mod(j + 2*k, 4), dp)/16.0_dp
        end do
    end do
    dss(0, :) = 0.0_dp
    call solve_transport_tangent(fine, 0.25_dp, dsc, dss, da, db, ok)
    if (.not. ok) error stop 'mixed source/nu tangent failed'
    call solve_transport_vjp(fine, ba, bb, gnu, bsc, bss, ok)
    if (.not. ok) error stop 'implicit coefficient VJP failed'
    lhs = sum(ba*da + bb*db)
    rhs = 0.25_dp*gnu + sum(dsc*bsc + dss*bss)
    if (abs(lhs - rhs) > 1.0e-11_dp) error stop 'nodal dual pairing identity'

    radius = cf%state_radius
    invalid = fine
    invalid%cosine(1, 33) = invalid%cosine(1, 33) + 1.0_dp
    if (transport_certificate_accepts(invalid, radius)) &
        error stop 'corrupt reconstructed candidate accepted'
    if (transport_certificate_accepts(fine, ieee_value(0.0_dp, ieee_quiet_nan))) &
        error stop 'NaN radius accepted'
    if (transport_certificate_accepts(empty, radius)) &
        error stop 'unallocated state accepted'
    invalid = fine
    invalid%nu = 0.0_dp
    if (transport_valid(invalid)) error stop 'zero damping accepted'
    invalid = fine
    invalid%v(2) = invalid%v(1)
    if (transport_certificate_accepts(invalid, radius)) &
        error stop 'collapsed cell accepted'
    invalid = fine
    invalid%source_sin(0, 1) = 1.0_dp
    if (transport_valid(invalid)) error stop 'nonphysical mean sine accepted'
    invalid = fine
    invalid%sine(0, 1) = 1.0_dp
    if (transport_certificate_accepts(invalid, radius)) &
        error stop 'nonphysical candidate mean sine accepted'
    invalid = fine
    invalid%nu = ieee_value(0.0_dp, ieee_quiet_nan)
    if (transport_certificate_accepts(invalid, radius)) &
        error stop 'NaN damping accepted'
    invalid = fine
    invalid%v(1) = -2.0_dp
    if (transport_certificate_accepts(invalid, radius)) &
        error stop 'changed endpoint accepted'
    invalid = fine
    deallocate (invalid%v)
    allocate (invalid%v(0:64))
    invalid%v = fine%v
    if (transport_valid(invalid)) error stop 'unsupported grid lower bound accepted'
    invalid = fine
    invalid%source_cos(1, 2) = ieee_value(0.0_dp, ieee_quiet_nan)
    if (transport_certificate_accepts(invalid, radius)) error stop 'NaN source accepted'
    dss(0, 1) = 1.0_dp
    call solve_transport_tangent(fine, 1.0_dp, dsc, dss, da, db, ok)
    if (ok) error stop 'invalid mean-sine source direction accepted'
    write (*, '(a,2es15.6)') 'General transport continuous radii (16/64): ', &
        cc%state_radius, cf%state_radius

contains

    subroutine contains_fixture(enclosure, exact)
        type(interval_t), intent(in) :: enclosure
        real(dp), intent(in) :: exact

        if (exact < enclosure%lo .or. exact > enclosure%hi) &
            error stop 'independent quadratic Bernstein fixture missed'
    end subroutine contains_fixture

    subroutine supplied_case(state, cells)
        type(transport_state_t), intent(out) :: state
        integer, intent(in) :: cells
        integer :: n
        logical :: solved

        call initialize_transport(state, cells, 2, 0.875_dp, solved)
        if (.not. solved) error stop 'transport setup failed'
        do n = 1, size(state%v)
            call transport_source_affine(state%v(n), 1.0_dp, 0.125_dp, &
                                         state%source_cos(0, n))
            call transport_source_affine(state%v(n), 0.25_dp, -0.25_dp, &
                                         state%source_cos(1, n))
            call transport_source_affine(state%v(n), 0.0_dp, 0.125_dp, &
                                         state%source_sin(2, n))
        end do
        call solve_transport(state, solved)
        if (.not. solved) error stop 'transport direct solve failed'
    end subroutine supplied_case

    subroutine test_only_resolvent(xq, vq, nuq, valueq, tangentq)
        real(qp), intent(in) :: xq, vq, nuq
        real(qp), intent(out) :: valueq, tangentq
        complex(qp) :: forcing, denominator, mode
        integer :: m

        ! Complex convention f=Re(z*exp(i*k*x)), z=(sc-i*ss)/(nu+i*k*v).
        ! Independent implementation/precision, forbidden in production checker.
        valueq = (1.0_qp + vq/8.0_qp)/nuq
        tangentq = -(1.0_qp + vq/8.0_qp)/nuq**2
        do m = 1, 2
            if (m == 1) then
                forcing = cmplx((1.0_qp - vq)/4.0_qp, 0.0_qp, qp)
            else
                forcing = cmplx(0.0_qp, -vq/8.0_qp, qp)
            end if
            denominator = cmplx(nuq, real(m, qp)*vq, qp)
            mode = exp(cmplx(0.0_qp, real(m, qp)*xq, qp))
            valueq = valueq + real(forcing*mode/denominator, qp)
            tangentq = tangentq - real(forcing*mode/denominator**2, qp)
        end do
    end subroutine test_only_resolvent

end program test_transport
