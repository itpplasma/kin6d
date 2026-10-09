program test_solovev_geometry
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use kin6d_solovev_geometry, only: solovev_geometry_t, geometry_sample_t, &
        initialize_solovev_geometry, sample_solovev_geometry, &
        solovev_geometry_jvp, solovev_geometry_vjp, &
        geometry_ok, geometry_invalid, geometry_interface, geometry_arithmetic
    implicit none
    type(solovev_geometry_t) :: map, bad
    type(geometry_sample_t) :: sample, plus, minus
    real(dp) :: polygon(2, 8), r0, level, tau, boundary(2), point(2), physical(2)
    real(dp) :: want_j(2, 2), want_h(2, 2, 2), delta, q, a, radial(2)
    real(dp) :: m(2, 2), direction(2), step, v(2)
    real(dp) :: b, t, k, corner(2), gamma, tangent(2), adjoint(2), weight(2), radius
    real(dp) :: wrong(2, 8)
    integer :: status, i, j, passed
    passed = 0
    r0 = 3; boundary = [0.3_dp, 0.2_dp]
    level = 5.3145_dp; tau = 0.75_dp; q = 1-tau
    ! Rotated regular octagon: inradius .288, outradius < .313.
    ! Omega lies in the radius-.446 disk and H < level on radius-.313:
    ! q*Omega is inside P and delta <= 1-.288/.446 gives a > .17.
    ! These strict fixture premises are independent of callback guard checks.
    gamma = 1.44_dp*(sqrt(2.0_dp)-1)
    do i = 1, 8
        select case (i)
        case (1)
            v = [-gamma, -1.44_dp]
        case (2)
            v = [gamma, -1.44_dp]
        case (3)
            v = [1.44_dp, -gamma]
        case (4)
            v = [1.44_dp, gamma]
        case (5)
            v = [gamma, 1.44_dp]
        case (6)
            v = [-gamma, 1.44_dp]
        case (7)
            v = [-1.44_dp, gamma]
        case (8)
            v = [-1.44_dp, -gamma]
        end select
        polygon(:, i) = [r0, 0.0_dp] + &
            [4*v(1)-3*v(2), 3*v(1)+4*v(2)]/25
    end do
    call initialize_solovev_geometry(map, r0, level, tau, polygon, status)
    call require(status == geometry_ok, 'explicit noncircular geometry initializes')
    call sample_solovev_geometry(map, [r0, 0.0_dp], sample)
    call require(sample%status == geometry_ok, 'axis has an identity neighborhood')
    call require(sample%identity, 'axis branch tagged identity')
    call close(maxval(abs(sample%value-[r0, 0.0_dp])), 0.0_dp, 0.0_dp, 'axis value')
    call close(sample%jacobian(1, 1), 1.0_dp, 0.0_dp, 'axis identity Jacobian')
    call close(sample%jacobian(2, 2), 1.0_dp, 0.0_dp, 'axis identity Jacobian')
    call close(sum(abs(sample%hessian)), 0.0_dp, 0.0_dp, 'axis identity Hessian')
    call require(.not. sample%certified, 'native evaluation makes no root certificate')
    ! Known forward image on the quartic ray. Independent formal-power-series
    ! residual solving gives the full jets below, without the closed jet formula.
    point = [r0+128.0_dp/625, 256.0_dp/1875]
    physical = [r0, 0.0_dp] + 0.8_dp*boundary
    call sample_solovev_geometry(map, point, sample)
    call require(sample%status == geometry_ok, 'smooth annular sample')
    call require(.not. sample%identity, 'annular branch')
    call require(sample%facet == 3, 'unique active facet selected')
    call close(maxval(abs(sample%value-physical)), 0.0_dp, 2.0e-13_dp, &
        'smaller root value')
    want_j = reshape([6885.0_dp/4864, 395.0_dp/2432, &
        2145.0_dp/9728, 6415.0_dp/4864], [2, 2])
    want_h(1, :, :) = reshape([8399001875.0_dp/3034202112.0_dp, &
        12750974375.0_dp/4045602816.0_dp, 12750974375.0_dp/4045602816.0_dp, &
        -286488125.0_dp/674267136], [2, 2])
    want_h(2, :, :) = reshape([73758125.0_dp/284456448, &
        15614155625.0_dp/6068404224.0_dp, 15614155625.0_dp/6068404224.0_dp, &
        3782714375.0_dp/2022801408.0_dp], [2, 2])
    call close(maxval(abs(sample%jacobian-want_j)), 0.0_dp, 3.0e-12_dp, &
        'independent noncircular full Jacobian')
    call close(maxval(abs(sample%hessian-want_h)), 0.0_dp, 3.0e-10_dp, &
        'independent noncircular full Hessian')
    direction = [0.4_dp, -0.7_dp]; weight = [-0.3_dp, 0.8_dp]
    call solovev_geometry_jvp(map, point, direction, tangent, status)
    call require(status == geometry_ok, 'implicit spatial JVP status')
    call close(maxval(abs(tangent-matmul(want_j, direction))), &
        0.0_dp, 3.0e-12_dp, 'JVP independent Jacobian oracle')
    call solovev_geometry_vjp(map, point, weight, adjoint, status)
    call require(status == geometry_ok, 'Cartesian pairing VJP status')
    call close(maxval(abs(adjoint-matmul(transpose(want_j), weight))), &
        0.0_dp, 3.0e-12_dp, 'VJP independent Jacobian oracle')
    call close(dot_product(weight, tangent), dot_product(adjoint, direction), &
        2.0e-14_dp, 'spatial JVP VJP Euclidean pairing')
    call solovev_geometry_jvp(map, point, [huge(1.0_dp), huge(1.0_dp)], &
        tangent, status)
    call require(status == geometry_arithmetic, 'overflowed derivative product fails closed')
    call close(norm2(tangent), 0.0_dp, 0.0_dp, 'overflowed product cleared')
    delta = 0.2_dp; a = 16.0_dp/25
    do i = 1, 2
        radial(i) = dot_product(boundary, matmul(sample%hessian(i, :, :), boundary))
    end do
    call close(maxval(abs(radial-boundary*2*delta/(tau*a**3))), &
        0.0_dp, 2.0e-11_dp, 'inverse radial curvature is nonzero')
    ! Nonsymmetric affine directions exercise all Cartesian chain terms.
    m = reshape([2.0_dp, 0.25_dp, 1.0_dp/3, 1.5_dp], [2, 2])
    step = 2.0e-6_dp
    do j = 1, 2
        direction = m(:, j)
        call sample_solovev_geometry(map, point+step*direction, plus)
        call sample_solovev_geometry(map, point-step*direction, minus)
        call require(plus%status == geometry_ok, 'positive displaced smooth point')
        call require(minus%status == geometry_ok, 'negative displaced smooth point')
        call close(maxval(abs((plus%value-minus%value)/(2*step) - &
            matmul(want_j, direction))), 0.0_dp, 2.0e-8_dp, 'affine direction JVP')
        do i = 1, 2
            call close(maxval(abs((plus%jacobian(i, :)-minus%jacobian(i, :))/(2*step) &
                - matmul(want_h(i, :, :), direction))), &
                0.0_dp, 2.0e-7_dp, 'affine direction Hessian')
        end do
    end do
    ! Independent forward radial equation on several known quartic-ray radii.
    do i = 1, 5
        t = 0.3_dp+0.1_dp*i; k = 0.8_dp
        b = t*(1-(1-k)*(t-q)/tau)
        point = [r0, 0.0_dp]+b*boundary
        call sample_solovev_geometry(map, point, sample)
        call require(sample%status == geometry_ok, 'forward-image smooth status')
        call close(maxval(abs(sample%value-([r0, 0.0_dp]+t*boundary))), &
            0.0_dp, 3.0e-13_dp, 'forward inverse composition')
    end do
    ! Polygon boundary maps to the exact curved boundary along this ray.
    point = [r0, 0.0_dp]+k*boundary
    call sample_solovev_geometry(map, point, sample)
    call require(sample%status == geometry_ok, 'open polygon edge supported')
    call close(maxval(abs(sample%value-([r0, 0.0_dp]+boundary))), &
        0.0_dp, 3.0e-13_dp, 'quartic boundary ray')
    call close(oracle_h(sample%value), level, 2.0e-12_dp, 'physical level boundary')
    do i = -1, 1, 2
        t = q+i*1.0e-8_dp
        b = t
        if (i > 0) b = t*(1-0.2_dp*(t-q)/tau)
        call sample_solovev_geometry(map, [r0, 0.0_dp]+b*boundary, sample)
        call require(sample%status == geometry_ok, &
            'strict one-sided interface neighborhood')
        call require(sample%identity .eqv. (i < 0), 'correct strict interface side')
        call close(maxval(abs(sample%value-([r0, 0.0_dp]+t*boundary))), &
            0.0_dp, 3.0e-13_dp, 'one-sided value continuity')
    end do
    ! The horizontal ray has an independent analytic quartic exit root.
    radius = sqrt(r0*r0+sqrt(level))-r0
    k = 1.44_dp*sqrt(2.0_dp)/(7*radius)
    t = 0.85_dp
    b = t*(1-(1-k)*(t-q)/tau)
    point = [r0+b*radius, 0.0_dp]
    call sample_solovev_geometry(map, point, sample)
    call require(sample%status == geometry_ok, 'analytic horizontal-axis ray status')
    call close(maxval(abs(sample%value-[r0+t*radius, 0.0_dp])), &
        0.0_dp, 3.0e-13_dp, 'analytic horizontal-axis ray inverse')
    call sample_solovev_geometry(map, [r0, 0.0_dp]+q*boundary, sample)
    call require(sample%status == geometry_interface, 'inner derivative cut rejected')
    call require(sample%inner_interface, 'inner cut explicitly tagged')
    call solovev_geometry_jvp(map, [r0, 0.0_dp]+q*boundary, direction, tangent, status)
    call require(status == geometry_interface, 'JVP does not cross the inner cut')
    call close(norm2(tangent), 0.0_dp, 0.0_dp, &
        'unavailable cut derivative clears output')
    corner = polygon(:, 4)-[r0, 0.0_dp]
    call sample_solovev_geometry(map, [r0, 0.0_dp]+0.98_dp*corner, sample)
    call require(sample%status == geometry_interface, 'sector derivative cut rejected')
    call require(sample%sector_interface, 'sector cut explicitly tagged')
    call close(sum(abs(sample%jacobian))+sum(abs(sample%hessian)), &
        0.0_dp, 0.0_dp, 'failed interface sample clears every jet')
    call sample_solovev_geometry(map, [r0, 0.0_dp]+1.2_dp*corner, sample)
    call require(sample%status == geometry_invalid, 'outside polygon rejected')
    call sample_solovev_geometry(map, [ieee_value(0.0_dp, ieee_quiet_nan), &
        0.0_dp], sample)
    call require(sample%status == geometry_invalid, 'nonfinite point rejected')
    wrong = polygon(:, 8:1:-1)
    call initialize_solovev_geometry(bad, r0, level, tau, wrong, status)
    call require(status == geometry_invalid, 'clockwise polygon rejected')
    call initialize_solovev_geometry(bad, r0, level, 0.0_dp, polygon, status)
    call require(status == geometry_invalid, 'zero tau rejected')
    call sample_solovev_geometry(bad, point, sample)
    call require(sample%status == geometry_invalid, &
        'failed initialization stays closed')
    print *, 'Solovev geometry primitive checks:', passed
contains
    real(dp) function oracle_h(x) result(value)
        real(dp), intent(in) :: x(2)
        value = (x(1)**2-r0**2)**2+4*x(1)**2*x(2)**2
    end function
    subroutine require(ok, label)
        logical, intent(in) :: ok
        character(*), intent(in) :: label
        if (.not. ok) then
            print *, 'FAIL: ', label
            error stop 1
        end if
        passed = passed+1
    end subroutine
    subroutine close(value, expected, tolerance, label)
        real(dp), intent(in) :: value, expected, tolerance
        character(*), intent(in) :: label
        call require(abs(value-expected) <= tolerance, label)
    end subroutine
end program
