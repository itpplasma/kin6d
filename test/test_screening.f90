!> Q7 part 1 oracles: T=0 Lindhard closed form, two independent finite-T
!> routes, kin6d_fermi_gas long-wavelength limit, Sommerfeld and Boltzmann
!> (Dawson) limits, exact Yukawa inversion, tail bound and Friedel period.
program test_screening
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_quadrature, only: gauss_legendre_ab
    use kin6d_species, only: atomic_units, electron
    use kin6d_screening_hierarchy
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp)
    type(screening_state_t) :: st, zero
    type(screening_value_t) :: a, b
    type(screening_q_table_t) :: table
    real(dp) :: xs(7) = [0.05_dp, 0.3_dp, 0.9_dp, 0.999_dp, 1.001_dp, 1.5_dp, 4.0_dp]
    real(dp) :: thetas(4) = [0.01_dp, 0.1_dp, 1.0_dp, 10.0_dp]
    real(dp) :: worst, y, ya, yb, err, rs, n, kf, ef, tail, ind, period
    real(dp), allocatable :: r(:), ratio(:), ratio2(:)
    integer :: i, j, failures, zeros
    logical :: ok
    failures = 0

    ! 1. Chart continuity of the generated T=0 Lindhard function.
    call check('series/closed chart at x=0.1', abs(lindhard_zero(0.1_dp) - &
        lindhard_zero(0.1_dp + 1e-15_dp)), 1e-14_dp)
    call check('closed/asymptotic chart at y=10', abs(lindhard_zero(10.0_dp) - &
        lindhard_zero(10.0_dp - 1e-14_dp)), 1e-15_dp)
    call check('F(1) = 1/2 limit', abs(lindhard_zero(1.0_dp + 1e-12_dp) - 0.5_dp), &
        1e-10_dp)

    ! 2. Lemma F(y) <= min(1, 2/(3y^2)) used by the inversion tail bound.
    ! y >= sqrt(2): series F <= 1/(3(y^2-1)) <= 2/(3 y^2) (proof in doc).
    ! sqrt(2/3) <= y <= sqrt(2): F decreasing, so F(ya) <= 2/(3 yb^2) per cell.
    worst = -1.0_dp
    do i = 0, 399
        ya = sqrt(2.0_dp/3.0_dp) + (sqrt(2.0_dp) - sqrt(2.0_dp/3.0_dp))*i/400
        yb = sqrt(2.0_dp/3.0_dp) + (sqrt(2.0_dp) - sqrt(2.0_dp/3.0_dp))*(i + 1)/400
        worst = max(worst, lindhard_zero(ya) - 2.0_dp/(3.0_dp*yb**2))
    end do
    call check('tail lemma margin (must be < 0)', worst + 1.0_dp, 1.0_dp - 1e-3_dp)

    ! 3. Route A at T=0 against the closed form.
    call screening_setup_normalized(0.0_dp, 1.0_dp, zero, ok)
    worst = 0.0_dp
    do i = 1, size(xs)
        a = lindhard_momentum(zero, xs(i))
        worst = max(worst, abs(a%value - lindhard_zero(xs(i))))
    end do
    write (*, '(a,es10.3)') 'T=0 momentum route vs closed form, max error ', worst
    call check('T=0 momentum route vs closed form', worst, 1e-12_dp)

    ! 4. Momentum (A) vs Maldague (B) routes at finite theta.
    worst = 0.0_dp
    do j = 1, size(thetas)
        call screening_setup_normalized(thetas(j), 1.0_dp, st, ok)
        if (.not. ok) error stop 'setup failed'
        do i = 1, size(xs)
            a = lindhard_momentum(st, xs(i)); b = lindhard_maldague(st, xs(i))
            worst = max(worst, abs(a%value - b%value))
            if (max(a%indicator, b%indicator) > 1e-11_dp) then
                write (*, *) 'indicator', thetas(j), xs(i), a%indicator, b%indicator
                failures = failures + 1
            end if
            if (max(a%tail, b%tail) > 1e-15_dp) failures = failures + 1
        end do
        ! 5. Long-wavelength limit vs kin6d_fermi_gas screening ratio.
        ! G(0) must lie in the certified fermi_gas enclosure of (ks/kTF0)^2.
        b = lindhard_maldague(st, 1e-6_dp)
        write (*, '(a,f6.2,a,es10.3)') 'theta=', thetas(j), &
            ' G(0)/fermi_gas candidate - 1 = ', b%value/st%ks_ratio**2 - 1
        if (b%value < st%ks_ratio_bound%lo**2 .or. b%value > st%ks_ratio_bound%hi**2) &
            failures = failures + 1
    end do
    write (*, '(a,es10.3)') 'finite-T route A vs route B, max difference ', worst
    call check('finite-T route agreement', worst, 1e-11_dp)

    ! 6. Sommerfeld limit: G(0) = sqrt(mu)(1 - pi^2 theta^2/24) + O(theta^4) with
    ! mu = 1 - pi^2 theta^2/12; finite-x deviation from F scales as theta^2.
    call screening_setup_normalized(0.01_dp, 1.0_dp, st, ok)
    if (.not. ok) error stop 'setup failed'
    b = lindhard_maldague(st, 1e-6_dp)
    err = abs(b%value - sqrt(1 - pi**2*1e-4_dp/12)*(1 - pi**2*1e-4_dp/24))
    write (*, '(a,es10.3)') 'theta=0.01 G(0) vs Sommerfeld expansion ', err
    call check('Sommerfeld G(0), O(theta^4) remainder', err, 1e-7_dp)
    a = lindhard_response(st, 0.5_dp)
    y = abs(a%value - lindhard_zero(0.5_dp))
    call screening_setup_normalized(0.02_dp, 1.0_dp, st, ok)
    a = lindhard_response(st, 0.5_dp)
    err = abs(a%value - lindhard_zero(0.5_dp))/y
    write (*, '(a,f8.4)') 'degenerate-limit deviation ratio theta 0.02/0.01 ', err
    call check('theta^2 approach to T=0 (ratio-4)', abs(err - 4), 0.05_dp)

    ! 7. Boltzmann limit: (3 theta/2) G -> Dawson(y)/y, y = x/sqrt(theta).
    call screening_setup_normalized(400.0_dp, 1.0_dp, st, ok)
    worst = 0.0_dp
    do i = 1, 3
        y = 0.1_dp*i
        a = lindhard_maldague(st, y*sqrt(st%theta))
        worst = max(worst, abs(1.5_dp*st%theta*a%value - dawson(y)/y))
    end do
    write (*, '(a,es10.3)') 'theta=400 vs Boltzmann Dawson response ', worst
    call check('nondegenerate Dawson limit', worst, 2e-3_dp)

    ! 8. Physical setup through species + fermi_gas: kappa2 = 4/(pi kF aB).
    rs = 2.0_dp
    n = 3.0_dp/(4*pi*rs**3)
    kf = (3*pi**2*n)**(1.0_dp/3.0_dp)
    ef = 0.5_dp*kf**2
    call screening_setup(atomic_units(), electron(), n, 0.5_dp*ef, st, ok)
    if (.not. ok) error stop 'physical setup failed'
    call check('theta from species', abs(st%theta - 0.5_dp), 1e-13_dp)
    call check('rs from species', abs(st%scales%rs - rs), 1e-13_dp)
    call check('kappa2 = 4/(pi kF aB)', abs(st%kappa2 - 4/(pi*kf)), 1e-12_dp)
    call check('kD/kTF0 = sqrt(2/(3 theta))', &
        abs(st%kd_ratio - sqrt(2/(3*st%theta))), 1e-12_dp)

    ! 9. Inversion with constant G is exactly Yukawa exp(-sqrt(kappa2 G0) r).
    allocate (r(60), ratio(60), ratio2(60))
    do i = 1, size(r)
        r(i) = 0.1_dp*1.1_dp**(i - 1)
    end do
    call screening_q_table(zero, 200.0_dp, 0.02_dp, 0.1_dp, table, constant_g=0.7_dp)
    call screened_ratio_from_table(zero, table, r, ratio, tail, ind)
    worst = maxval(abs(ratio - exp(-sqrt(0.7_dp)*r)))
    ! Constant G violates the Lindhard tail lemma; its own q-tail is bounded by
    ! (2/pi) kappa2 G0 integral_qmax^inf dq/q^3 = kappa2 G0/(pi qmax^2).
    write (*, '(a,es10.3,a,es10.3)') 'constant-G inversion vs Yukawa ', worst, &
        ' own tail bound ', 0.7_dp/(pi*200.0_dp**2)
    call check('Yukawa inversion oracle', worst, 0.7_dp/(pi*200.0_dp**2) + 1e-10_dp)

    ! 10. q-truncation bound encloses the change from qmax = 50 to 200.
    call screening_q_table(zero, 50.0_dp, 0.02_dp, 0.1_dp, table)
    call screened_ratio_from_table(zero, table, r, ratio2, tail, ind)
    call screening_q_table(zero, 200.0_dp, 0.02_dp, 0.1_dp, table)
    call screened_ratio_from_table(zero, table, r, ratio, err, ind)
    write (*, '(a,es10.3,a,es10.3)') 'qmax 50->200 change ', &
        maxval(abs(ratio - ratio2)), ' bound(50) ', tail
    call check('truncation bound encloses change', &
        maxval(abs(ratio - ratio2)) - tail - 2*ind, 0.0_dp)

    ! 11. Friedel period pi/kF at T=0 from sign changes on 20 <= r <= 60.
    deallocate (r, ratio)
    allocate (r(1601), ratio(1601))
    do i = 1, size(r)
        r(i) = 20.0_dp + 40.0_dp*(i - 1)/1600
    end do
    call screened_ratio(zero, r, ratio, tail, ind)
    zeros = 0; ya = 0; yb = 0
    do i = 1, size(r) - 1
        if (ratio(i)*ratio(i + 1) < 0) then
            y = r(i) - ratio(i)*(r(i + 1) - r(i))/(ratio(i + 1) - ratio(i))
            if (zeros == 0) ya = y
            yb = y; zeros = zeros + 1
        end if
    end do
    period = 2*(yb - ya)/(zeros - 1)
    write (*, '(a,f10.6,a,i0,a)') 'Friedel period kF*lambda = ', period, &
        ' from ', zeros, ' zeros'
    call check('Friedel period pi/kF', abs(period/pi - 1), 2e-3_dp)

    if (failures > 0) error stop 'screening oracles failed'
    write (*, '(a)') 'screening oracles passed'
contains
    subroutine check(label, value, tolerance)
        character(*), intent(in) :: label
        real(dp), intent(in) :: value, tolerance
        if (.not. (value <= tolerance)) then
            write (*, '(a,a,es12.4,a,es10.3)') 'FAIL ', label, value, ' > ', tolerance
            failures = failures + 1
        end if
    end subroutine
    !> Dawson integral D(y) = exp(-y^2) integral_0^y exp(t^2) dt (oracle).
    real(dp) function dawson(y) result(d)
        real(dp), intent(in) :: y
        real(dp) :: t(40), w(40)
        call gauss_legendre_ab(40, 0.0_dp, y, t, w)
        d = exp(-y*y)*sum(w*exp(t*t))
    end function
end program test_screening
