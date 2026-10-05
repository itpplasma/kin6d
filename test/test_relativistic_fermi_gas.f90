program test_relativistic_fermi_gas
    ! Q8 relativistic ideal Fermi EOS. Independent oracle: binary128 composite
    ! Simpson moments and bisection written here from E=m*c**2*(sqrt(1+p^2/m^2c^2)-1);
    ! no production leaf is used for expected values.
    use, intrinsic :: iso_fortran_env, only: dp => real64, qp => real128
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortnum_interval, only: interval_t, interval
    use kin6d_fermi_gas, only: fermi_workspace_t, fermi_state_t, fermi_si_state_t, &
        initialize_fermi_workspace, solve_fermi_gas, solve_fermi_gas_si, &
        occupation_fermi_gas, fermi_gas_jvp
    implicit none
    type(fermi_workspace_t) :: work
    type(fermi_state_t) :: state, classic, corrupted, plus, minus
    type(fermi_si_state_t) :: si
    real(dp), parameter :: thetas(4) = [0.05_dp, 1.0_dp, 10.0_dp, 0.002_dp]
    real(dp), parameter :: xs(4) = [1.0_dp, 1.0_dp, 3.0_dp, 5.0_dp]
    real(qp) :: reference(6), gradient(6), refined(6), refined_gradient(6), chandra
    real(dp) :: value, lower, upper, tangent(6), h, fd
    type(interval_t) :: bands(6), tangent_band(6)
    integer :: case_index, j
    logical :: ok

    call initialize_fermi_workspace(work, 65536, ok)
    call require(ok, 'workspace')
    do case_index = 1, size(thetas)
        call solve_fermi_gas(work, thetas(case_index), state, ok, &
                             relativity=xs(case_index))
        call require(ok .and. state%certified, 'relativistic EOS certificate')
        call oracle(real(thetas(case_index), qp), real(xs(case_index), qp), 512, &
                    reference, gradient)
        call oracle(real(thetas(case_index), qp), real(xs(case_index), qp), 1024, &
                    refined, refined_gradient)
        call require(maxval(abs(reference - refined)) < 1.0e-10_qp, &
                     'relativistic oracle refinement')
        call require(maxval(abs(gradient - refined_gradient)) < 1.0e-8_qp, &
                     'relativistic derivative oracle refinement')
        bands = [state%mu_bound, state%pressure_bound, state%energy_bound, &
            state%susceptibility_bound, state%screening_ratio_bound, &
            state%bulk_modulus_bound]
        do j = 1, 6
            call require(encloses(bands(j), refined(j)), 'relativistic containment')
            call require(encloses(state%theta_gradient_bound(j), &
                refined_gradient(j)), 'relativistic implicit derivative containment')
        end do
        call require(encloses(state%density_bound, 1.0_qp), 'relativistic density')
        ! Per particle E_kin/(p*v)=s/(s+1) lies in [1/2,1): 1.5P<=U<3P.
        call require(state%energy_bound%lo >= 1.5_dp*state%pressure_bound%lo .and. &
            state%energy_bound%hi < 3.0_dp*state%pressure_bound%hi, &
            'relativistic energy-pressure window')
        call occupation_fermi_gas(state, 1.0_dp, value, lower, upper, ok)
        call require(ok .and. lower <= value .and. value <= upper, &
                     'relativistic occupation query')
        call require(encloses(interval(lower, upper), occupation(1.0_qp, &
            real(thetas(case_index), qp), refined(1), real(xs(case_index), qp))) &
            .or. abs(value - real(occupation(1.0_qp, real(thetas(case_index), qp), &
            refined(1), real(xs(case_index), qp)), dp)) < 1.0e-6_dp, &
            'relativistic occupation oracle')
    end do

    ! Nonrelativistic limit: x->0 reproduces the existing certified EOS.
    call solve_fermi_gas(work, 1.0_dp, classic, ok)
    call require(ok, 'nonrelativistic reference')
    call solve_fermi_gas(work, 1.0_dp, state, ok, relativity=1.0e-6_dp)
    call require(ok, 'small-x relativistic solve')
    call require(abs(state%pressure - classic%pressure) < 1.0e-9_dp .and. &
        abs(state%mu - classic%mu) < 1.0e-9_dp .and. &
        abs(state%energy - classic%energy) < 1.0e-9_dp, 'nonrelativistic limit')
    call require(overlap(state%pressure_bound, classic%pressure_bound), &
                 'nonrelativistic limit enclosure')
    call solve_fermi_gas(work, 1.0_dp, state, ok, relativity=0.0_dp)
    call require(ok .and. state%pressure <= classic%pressure .and. &
        state%pressure >= classic%pressure, 'x=0 is the nonrelativistic path')

    ! Cold ultrarelativistic limit P=(hbar*c/4)*(3*pi**2)**(1/3)*n**(4/3),
    ! normalized 1/(2x), and the exact Chandrasekhar T=0 curve.
    call solve_fermi_gas(work, 1.0e-4_dp, state, ok, relativity=50.0_dp)
    call require(ok, 'cold ultrarelativistic certificate')
    chandra = chandrasekhar(50.0_qp)
    call require(abs(state%pressure*100.0_dp - 1.0_dp) < 1.0e-3_dp, &
                 'ultrarelativistic cold limit')
    call require(abs(real(state%pressure, qp)/chandra - 1.0_qp) < 1.0e-4_qp, &
                 'Chandrasekhar cold curve')
    call solve_fermi_gas(work, 1.0e-4_dp, state, ok, relativity=0.1_dp)
    call require(ok .and. abs(real(state%pressure, qp)/chandrasekhar(0.1_qp) - &
                 1.0_qp) < 1.0e-6_qp, 'Chandrasekhar weakly relativistic')
    ! Classical limit P=n*kB*T for every x.
    call solve_fermi_gas(work, 300.0_dp, state, ok, relativity=3.0_dp)
    call require(ok .and. abs(state%pressure/300.0_dp - 1.0_dp) < 1.0e-3_dp, &
                 'classical relativistic limit')

    ! Equation-level derivative versus centered candidate differences.
    call solve_fermi_gas(work, 0.3_dp, state, ok, relativity=2.0_dp)
    call require(ok, 'derivative base state')
    h = 1.0e-4_dp
    call solve_fermi_gas(work, 0.3_dp + h, plus, ok, certify=.false., relativity=2.0_dp)
    call solve_fermi_gas(work, 0.3_dp - h, minus, ok, certify=.false., relativity=2.0_dp)
    fd = (plus%pressure - minus%pressure)/(2.0_dp*h)
    call require(abs(fd - state%theta_gradient(2)) < 1.0e-5_dp*max(1.0_dp, abs(fd)), &
                 'relativistic pressure derivative finite difference')
    call fermi_gas_jvp(state, 2.0_dp, tangent, tangent_band, ok)
    call require(ok .and. abs(tangent(2) - 2.0_dp*state%theta_gradient(2)) < &
                 1.0e-14_dp*abs(tangent(2)), 'relativistic jvp')
    corrupted = state
    corrupted%relativity = 2.5_dp
    call occupation_fermi_gas(corrupted, 1.0_dp, value, lower, upper, ok)
    call require(.not. ok, 'altered relativity rejected')

    ! SI electron gas at 1e9 kg/m**3, Ye=1/2, T=1e7 K: x=pF/(m*c) near one.
    call solve_fermi_gas_si(work, 1.0e9_dp*0.5_dp/1.66053906660e-27_dp, 1.0e7_dp, &
        9.1093837015e-31_dp, 1.054571817e-34_dp, 1.380649e-23_dp, &
        1.602176634e-19_dp, 8.8541878128e-12_dp, si, ok, speed=299792458.0_dp)
    call require(ok .and. si%certified, 'relativistic SI certificate')
    call require(abs(si%normalized%relativity - si%fermi_momentum/ &
        (9.1093837015e-31_dp*299792458.0_dp)) < 1.0e-14_dp, 'SI relativity parameter')
    call require(si%normalized%relativity > 0.8_dp .and. &
        si%normalized%relativity < 1.2_dp, 'SI relativity range')
    call require(si%pressure_bound%lo <= si%pressure .and. &
                 si%pressure <= si%pressure_bound%hi, 'SI pressure enclosure')
    print '(a)', 'relativistic Fermi EOS oracles passed'

contains

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(*), intent(in) :: message
        if (.not. condition) then
            print '(a)', 'FAIL: '//message
            error stop 1
        end if
    end subroutine require

    logical function encloses(bound, expected) result(ok)
        type(interval_t), intent(in) :: bound
        real(qp), intent(in) :: expected
        ok = .false.
        if (.not. all(ieee_is_finite([bound%lo, bound%hi]))) return
        if (bound%lo > bound%hi) return
        ok = real(bound%lo, qp) <= expected .and. expected <= real(bound%hi, qp)
    end function encloses

    logical function overlap(a, b)
        type(interval_t), intent(in) :: a, b
        overlap = a%lo <= b%hi .and. b%lo <= a%hi
    end function overlap

    real(qp) function energy(r, x)
        real(qp), intent(in) :: r, x
        energy = 2.0_qp*r*r/(sqrt(1.0_qp + x*x*r*r) + 1.0_qp)
    end function energy

    real(qp) function occupation(r, theta, mu, x)
        real(qp), intent(in) :: r, theta, mu, x
        occupation = 1.0_qp/(1.0_qp + exp((energy(r, x) - mu)/theta))
    end function occupation

    !> Exact T=0 P/(n*EF) = [x*(2x^2-3)*sqrt(1+x^2)+3*asinh(x)]/(4x^5).
    real(qp) function chandrasekhar(x)
        real(qp), intent(in) :: x
        chandrasekhar = (x*(2.0_qp*x*x - 3.0_qp)*sqrt(1.0_qp + x*x) + &
                         3.0_qp*asinh(x))/(4.0_qp*x**5)
    end function chandrasekhar

    subroutine oracle(theta, x, panels, values, gradients)
        real(qp), intent(in) :: theta, x
        integer, intent(in) :: panels
        real(qp), intent(out) :: values(6), gradients(6)
        real(qp) :: lower, upper, mu, m(11), dmu, dchi
        integer :: iteration
        lower = -1.0_qp
        do
            call simpson(theta, lower, x, panels, m, .true.)
            if (3.0_qp*m(1) < 1.0_qp) exit
            lower = 2.0_qp*lower - 4.0_qp*theta
        end do
        upper = 5.0_qp
        do iteration = 1, 90
            mu = (lower + upper)/2.0_qp
            call simpson(theta, mu, x, panels, m, .true.)
            if (3.0_qp*m(1) < 1.0_qp) then
                lower = mu
            else
                upper = mu
            end if
        end do
        mu = (lower + upper)/2.0_qp
        call simpson(theta, mu, x, panels, m, .false.)
        values = [mu, m(2), 3.0_qp*m(3), 3.0_qp*m(4), sqrt(2.0_qp*m(4)), &
                  1.0_qp/(3.0_qp*m(4))]
        dmu = -m(5)/m(4)
        dchi = 3.0_qp*(m(11) + m(10)*dmu)
        gradients = [dmu, m(7) + m(6)*dmu, 3.0_qp*(m(9) + m(8)*dmu), dchi, &
                     dchi/(3.0_qp*values(5)), -dchi/values(4)**2]
    end subroutine oracle

    !> Composite Simpson on [0,1-w],[1-w,1+w],[1+w,R] if degenerate. Otherwise
    !> one Simpson rule in the rapidity t, r=sinh(t)/x, whose even entire
    !> integrands remove the r=+-i/x branch-point limitation.
    subroutine simpson(theta, mu, x, panels, m, density_only)
        real(qp), intent(in) :: theta, mu, x
        integer, intent(in) :: panels
        real(qp), intent(out) :: m(11)
        logical, intent(in) :: density_only
        real(qp) :: edges(4), level, r, h, f, l, e, v, w, g(11), jacobian, tt
        integer :: segment, node, count
        logical :: rapidity
        level = 4.0_qp + 80.0_qp*theta
        edges = [0.0_qp, 0.0_qp, 0.0_qp, sqrt(level*(1.0_qp + x*x*level/4.0_qp))]
        w = 25.0_qp*theta*sqrt(1.0_qp + x*x)
        rapidity = w >= 0.5_qp
        if (rapidity) then
            edges(4) = asinh(x*edges(4))
            edges(2:3) = [edges(4)/3.0_qp, 2.0_qp*edges(4)/3.0_qp]
        else
            edges(2:3) = [1.0_qp - w, 1.0_qp + w]
        end if
        m = 0.0_qp
        do segment = 1, 3
            h = (edges(segment + 1) - edges(segment))/real(panels, qp)
            do node = 0, panels
                r = edges(segment) + real(node, qp)*h
                jacobian = 1.0_qp
                if (rapidity) then
                    tt = r
                    r = sinh(tt)/x
                    jacobian = cosh(tt)/x
                end if
                f = occupation(r, theta, mu, x)
                l = f*(1.0_qp - f)
                e = energy(r, x)
                v = 2.0_qp*r/sqrt(1.0_qp + x*x*r*r)
                g = 0.0_qp
                g(1) = r*r*f
                if (.not. density_only) then
                    g(2:11) = [r**3*v*f, r*r*e*f, r*r*l/theta, &
                        r*r*(e - mu)*l/theta**2, r**3*v*l/theta, &
                        r**3*v*(e - mu)*l/theta**2, r*r*e*l/theta, &
                        r*r*e*(e - mu)*l/theta**2, r*r*l*(1.0_qp - 2.0_qp*f)/theta**2, &
                        r*r*l*((e - mu)*(1.0_qp - 2.0_qp*f)/theta**3 - 1.0_qp/theta**2)]
                end if
                count = 2
                if (mod(node, 2) == 1) count = 4
                if (node == 0 .or. node == panels) count = 1
                m = m + real(count, qp)*h*jacobian/3.0_qp*g
            end do
        end do
    end subroutine simpson
end program test_relativistic_fermi_gas
