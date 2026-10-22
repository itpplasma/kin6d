program test_nonnested_projection
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_special_elliptic, only: ellip_k, ellip_e
    use kin6d_nonnested_projection, only: nn_projection_t, nn_projection_report_t, &
             nn_projection_build, NN_CLOSED_LINE, NN_INVARIANT_SURFACE, NN_MIXING_VOLUME
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp), eps = 0.02_dp
    type(nn_projection_t) :: projection
    type(nn_projection_report_t) :: report
    real(dp), allocatable :: target(:), weights(:), projected(:), again(:), expected(:)
    integer, allocatable :: component(:)
    real(dp) :: angle, m, k, x, h, error, previous_error, xmax, a
    real(dp) :: last_mean(3), means(3)
    real(dp), parameter :: abscissa(3) = [-sqrt(0.6_dp), 0.0_dp, sqrt(0.6_dp)]
    real(dp), parameter :: gw(3) = [5.0_dp/9, 8.0_dp/9, 5.0_dp/9]
    integer :: n, lev, i, j, c, st, index
    ! Integrable irrational tori: a tensor quadrature of invariant surface
    ! measure, with s_ref constant on each surface. Non-invariant modes average out.
    n = 64
    allocate (target(3*n*n), weights(3*n*n), component(3*n*n), &
              projected(3*n*n), again(3*n*n), expected(3*n*n))
    do c = 1, 3
        do j = 1, n
            do i = 1, n
                index = (c - 1)*n*n + (j - 1)*n + i
                component(index) = c; weights(index) = 1.0_dp/(n*n)
                target(index) = (0.2_dp*c)**2 + cos(2*pi*i/n) + sin(2*pi*j/n)
                expected(index) = (0.2_dp*c)**2
            end do
        end do
    end do
    call nn_projection_build(component, [NN_INVARIANT_SURFACE, NN_INVARIANT_SURFACE, &
                                         NN_INVARIANT_SURFACE], weights, projection, st)
    if (st /= 0) error stop 'integrable projection build'
    call check_projection('integrable surfaces')
    if (maxval(abs(projected - expected)) > 1.0e-6_dp) error stop 'integrable target'
    deallocate (target, weights, component, projected, again, expected)

    ! Closed rational lines with theta' = 1 + a cos theta, phi' = theta'/2.
    ! dl/|B| = dtheta/(1+a cos theta). A uniform-arclength average is wrong.
    previous_error = huge(1.0_dp)
    do lev = 1, 4
        n = 8*2**(lev - 1); a = 0.6_dp
     allocate (target(n), weights(n), component(n), projected(n), again(n), expected(n))
        do i = 1, n
            angle = 2*pi*(i - 0.5_dp)/n
            target(i) = 2 + cos(angle)
            weights(i) = 2*pi/n/(1 + a*cos(angle))
        end do
        component = 1
        expected = 2 + (sqrt(1 - a*a) - 1)/a
        call nn_projection_build(component, [NN_CLOSED_LINE], weights, projection, st)
        if (st /= 0) error stop 'rational projection build'
        call check_projection('closed rational dl/absB')
        error = maxval(abs(projected - expected))
        print '(a,i5,es13.4)', 'rational trace ladder ', n, error
        if (error > previous_error + 1.0e-12_dp) error stop 'rational trace refinement'
        previous_error = error
        deallocate (target, weights, component, projected, again, expected)
    end do
    if (error > 1.0e-6_dp) error stop 'rational target'

    ! Analytic helical island H=x^2/2+eps cos(alpha), alpha=theta-phi.
    ! alpha=pi+2 asin(k sin u), x=2 k sqrt(eps) cos u,
    ! dt/du=1/(sqrt(eps)*sqrt(1-k^2 sin^2 u)). Independent complete
    ! elliptic integrals give <x^2>=4 eps (E/K - 1 + k^2).
    ! Each k is a different internal invariant surface, not a mixing region.
    last_mean = 0
    do lev = 1, 5
        n = 16*2**(lev - 1)
        allocate (target(3*n), weights(3*n), component(3*n), projected(3*n), &
                  again(3*n), expected(3*n))
        do c = 1, 3
            select case (c)
            case (1)
                k = 0.3_dp
            case (2)
                k = 0.8_dp
            case (3)
                k = 0.97_dp
            end select
            m = k*k; h = eps*(2*m - 1)
            do i = 1, n
                index = (c - 1)*n + i
                angle = 2*pi*(i - 0.5_dp)/n
                x = 2*k*sqrt(eps)*cos(angle)
                component(index) = c
                weights(index) = 2*pi/n/sqrt(eps*(1 - m*sin(angle)**2))
                target(index) = h*h + x*x
                expected(index) = h*h + 4*eps*(ellip_e(m)/ellip_k(m) - 1 + m)
            end do
        end do
        call nn_projection_build(component, &
            [NN_CLOSED_LINE, NN_CLOSED_LINE, NN_CLOSED_LINE], &
                                 weights, projection, st)
        if (st /= 0) error stop 'island projection build'
        call check_projection('helical island internal surfaces')
        means = [projected(1), projected(n + 1), projected(2*n + 1)]
        error = maxval(abs(projected - expected))
        print '(a,i5,2es13.4)', 'island trace ladder / analytic error / change ', &
            n, error, maxval(abs(means - last_mean))
        if (lev == 5) then
            if (error > 1.0e-6_dp) error stop 'island invariant target'
        if (maxval(abs(means - last_mean)) > 1.0e-6_dp) error stop 'island trace defect'
      if (maxval(means) - minval(means) < 0.005_dp) error stop 'internal surfaces mixed'
        end if
        last_mean = means
        deallocate (target, weights, component, projected, again, expected)
    end do

    ! Explicit whole-island mixing uses volume measure dx d(alpha), not
    ! line-time weights. Three-point Gauss integrates x^2 on each chord.
    n = 1024
    allocate (target(3*n), weights(3*n), component(3*n), projected(3*n), &
              again(3*n), expected(3*n))
    do i = 1, n
        angle = 2*pi*(i - 0.5_dp)/n
        xmax = 2*sqrt(eps)*sin(angle/2)
        do j = 1, 3
            index = 3*(i - 1) + j
            x = xmax*abscissa(j)
            target(index) = x*x; weights(index) = 2*pi/n*xmax*gw(j)
        end do
    end do
    component = 1; expected = 8*eps/9
    call nn_projection_build(component, [NN_MIXING_VOLUME], weights, projection, st)
    if (st /= 0) error stop 'mixing projection build'
    call check_projection('declared whole-island volume mixing')
    if (maxval(abs(projected - expected)) > 1.0e-6_dp) error stop 'mixing volume mean'
    if (maxval(projected) - minval(projected) > 1.0e-12_dp) error stop 'mixing not flat'
    ! Open components and nonpositive measures are unsupported, not averaged.
    call nn_projection_build(component, [0], weights, projection, st)
    if (st == 0) error stop 'open projection accepted'
    weights(1) = 0
    call nn_projection_build(component, [NN_MIXING_VOLUME], weights, projection, st)
    if (st == 0) error stop 'nonpositive measure accepted'
    print *, 'PASS: G1b declared-component projection diagnostic only'
contains
    subroutine check_projection(name)
        character(*), intent(in) :: name
        real(dp) :: mean_before, pairing
        call projection%apply(target, projected, report, st)
        if (st /= 0) error stop 'projection'
        mean_before = sum(weights*target)/sum(weights)
        if (abs(sum(weights*projected)/sum(weights) - mean_before) > 1.0e-9_dp) &
            error stop 'projection mean'
        if (report%mean_defect > 1.0e-9_dp) error stop 'component mean'
        call projection%apply(projected, again, report, st, previous=projected)
        if (st /= 0) error stop 'second projection'
     if (maxval(abs(again - projected)) > 1.0e-9_dp) error stop 'projection idempotence'
        if (report%trace_change > 1.0e-9_dp) error stop 'projection repeat defect'
        ! Weighted orthogonality is independent of idempotence and mean tests.
        pairing = sum(weights*projected*(target - projected))
        if (abs(pairing) > 1.0e-9_dp*sum(weights)) error stop 'projection orthogonality'
        print '(a,a)', 'projection verified: ', name
    end subroutine
end program
