program test_ocp
    ! Oracles of the periodic OCP cases (C3, C4) and their estimators:
    ! plasma dispersion function values/derivative identity, literature Landau
    ! roots, the Volterra initial-value response (cold limit, late-time Landau
    ! root), the periodic screened-Poisson Green function (Helmholtz equation by
    ! finite differences, isolated Yukawa limit), g(r)/S(k) estimators on an
    ! ideal gas and a simple-cubic lattice, plasma units against the shared
    ! species scales, and the step tangent of a mode amplitude by central
    ! differences.
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use fortnum_status, only: fortnum_status_t, FORTNUM_OK
    use fortnum_rng, only: rng_t, rng_seed, rng_uniform
    use kin6d_vlasov_response, only: plasma_dispersion, landau_root, &
        vlasov_density_response, fit_damped_cosine
    use kin6d_pair_correlation, only: pair_counts, pair_distribution, structure_factor
    use kin6d_ocp, only: ocp_model_t, ocp_system, ocp_prepare, ocp_advance, &
        ocp_advance_jvp, ocp_mode, ocp_displace, ocp_displace_tangent, &
        ocp_debye_green, ocp_kinetic_temperature
    use kin6d_nbody, only: nbody_system_t, nbody_energy
    implicit none
    real(dp), parameter :: pi = 3.14159265358979323846_dp
    type(fortnum_status_t) :: status

    call dispersion()
    call volterra()
    call debye_green()
    call estimators()
    call ocp_dynamics()
    print '(a)', 'PASS: Z function, Landau roots, Volterra, Debye Green, g/S, OCP'
contains
    subroutine check(condition, message)
        logical, intent(in) :: condition
        character(*), intent(in) :: message
        if (.not. condition) then
            print '(a)', 'FAIL: '//message
            error stop 1
        end if
    end subroutine

    subroutine dispersion()
        complex(dp) :: z, dz, zp, zm, dummy, omega
        complex(dp), parameter :: zeta = (1.3_dp, -0.4_dp), h = (1.0e-5_dp, 0.0_dp)
        logical :: ok
        call plasma_dispersion((0.0_dp, 0.0_dp), z, dz)
        call check(abs(z - cmplx(0.0_dp, sqrt(pi), dp)) < 1e-15_dp, 'Z(0) = i sqrt(pi)')
        ! Z(1) = i sqrt(pi)/e - 2 D(1), Dawson D(1) = 0.5380795069127684.
        call plasma_dispersion((1.0_dp, 0.0_dp), z, dz)
        call check(abs(z - cmplx(-2*0.5380795069127684_dp, sqrt(pi)*exp(-1.0_dp), dp)) &
                   < 1e-14_dp, 'Z(1) from Dawson D(1)')
        ! Large argument: Z ~ -1/zeta (1 + 1/(2 zeta^2) + 3/(4 zeta^4)) + ...
        call plasma_dispersion((5.0_dp, 0.0_dp), z, dz)
        call check(abs(real(z) + 0.2_dp*(1 + 0.02_dp + 0.0012_dp + 15/8.0_dp/5**6)) &
                   < 5e-6_dp, 'Z asymptotic series')
        call plasma_dispersion(zeta + h, zp, dummy)
        call plasma_dispersion(zeta - h, zm, dummy)
        call plasma_dispersion(zeta, z, dz)
        call check(abs((zp - zm)/(2*h) - dz) < 1e-8_dp, "Z' = -2(1 + zeta Z)")
        ! Literature (e.g. Canosa 1973; Fried-Conte tables): least-damped roots.
        call landau_root(0.5_dp, omega, ok)
        call check(ok .and. abs(real(omega) - 1.4156_dp) < 1e-4_dp .and. &
                   abs(aimag(omega) + 0.1534_dp) < 2e-4_dp, 'Landau root k lambda_D = 0.5')
        call landau_root(0.3_dp, omega, ok)
        call check(ok .and. abs(real(omega) - 1.1598_dp) < 2e-4_dp .and. &
                   abs(aimag(omega) + 0.0126_dp) < 2e-4_dp, 'Landau root k lambda_D = 0.3')
    end subroutine

    subroutine volterra()
        integer, parameter :: m = 4001
        real(dp) :: n(m), t(m), omega, gamma, amplitude
        complex(dp) :: root
        logical :: ok
        integer :: j
        real(dp), parameter :: h = 0.01_dp
        do j = 1, m
            t(j) = (j - 1)*h
        end do
        ! Cold limit k -> 0: N'' = -N, N(0) = 1, N'(0) = 0.
        call vlasov_density_response(1.0e-4_dp, h, n(:2001))
        call check(maxval(abs(n(:2001) - cos(t(:2001)))) < 1e-5_dp, 'cold Volterra limit')
        ! Late time of the full response is the least-damped Landau root.
        call vlasov_density_response(0.5_dp, h, n)
        call landau_root(0.5_dp, root, ok)
        call fit_damped_cosine(t(1001:2001), n(1001:2001), 1.2_dp, 1.6_dp, omega, &
                               gamma, amplitude)
        call check(abs(omega - real(root)) < 1e-3_dp .and. &
                   abs(gamma + aimag(root)) < 1e-3_dp, 'Volterra late time = Landau root')
    end subroutine

    subroutine debye_green()
        real(dp) :: d(3), g0, lap, b, h, length, v
        integer :: a
        length = 4.0_dp; h = 1.0e-3_dp
        d = [0.7_dp, -0.4_dp, 1.1_dp]
        g0 = ocp_debye_green(d, 1.0_dp, length, b)
        call check(b < 1e-16_dp, 'Debye image truncation bound')
        lap = 0
        do a = 1, 3
            lap = lap + ocp_debye_green(d + h*unit(a), 1.0_dp, length, b) + &
                  ocp_debye_green(d - h*unit(a), 1.0_dp, length, b) - 2*g0
        end do
        lap = lap/h**2
        v = length**3
        call check(abs(-lap + g0 + 4*pi/v) < 1e-5_dp*4*pi/v, &
                   'periodic screened Poisson: (-lap + kappa^2) G = -4 pi/V')
        g0 = ocp_debye_green([1.0_dp, 0.0_dp, 0.0_dp], 1.0_dp, 20.0_dp, b)
        call check(abs(g0 - (exp(-1.0_dp) - 4*pi/20.0_dp**3)) < 1e-8_dp, &
                   'large box: isolated Yukawa minus zero mode')
    end subroutine

    function unit(a) result(e)
        integer, intent(in) :: a
        real(dp) :: e(3)
        e = 0; e(a) = 1
    end function

    subroutine estimators()
        integer, parameter :: n = 4000, bins = 20, lattice = 10
        type(rng_t) :: g
        real(dp) :: x(3*n), counts(bins), r(bins), gr(bins), s(9), y(3*lattice**3)
        real(dp) :: mean, length, c2(21)
        integer :: vectors(9), i, j, k, p
        length = 10.0_dp
        call rng_seed(g, 99_int64, status)
        do i = 1, 3*n
            call rng_uniform(g, x(i)); x(i) = length*x(i)
        end do
        counts = 0
        call pair_counts(x, n, length, 2.0_dp, counts)
        call pair_distribution(counts, n, length, 2.0_dp, 1, r, gr)
        mean = sum(counts(6:))/sum(counts(6:)/gr(6:))
        call check(abs(mean - 1) < 0.01_dp, 'ideal gas g(r) = 1 (mean)')
        call check(all(abs(gr(6:) - 1) < 5/sqrt(counts(6:))), 'ideal gas g(r) = 1 (bins)')
        call structure_factor(x, n, length, 3, s, vectors)
        mean = sum(s*vectors)/sum(vectors)
        call check(abs(mean - 1) < 3/sqrt(real(sum(vectors), dp)), 'ideal gas S(k) = 1')
        ! Simple-cubic lattice of spacing 1: S = 0 at every non-Bragg k, six
        ! nearest neighbours at distance 1.
        p = 0
        do k = 0, lattice - 1
            do j = 0, lattice - 1
                do i = 0, lattice - 1
                    y(p + 1:p + 3) = [i, j, k] + 0.5_dp; p = p + 3
                end do
            end do
        end do
        call structure_factor(y, lattice**3, length, 3, s, vectors)
        call check(maxval(s) < 1e-12_dp, 'lattice S(k) = 0 off Bragg')
        c2 = 0
        call pair_counts(y, lattice**3, length, 2.1_dp, c2)
        call check(nint(sum(c2(9:12))) == 3*lattice**3, 'lattice first shell: 6 neighbours')
    end subroutine

    !> Tiny OCP (Lambda = 2, N = 128): energy conservation of the
    !> encounter-resolved step and the tangent of Re rho_k(t) with respect
    !> to the displacement amplitude against central differences.
    subroutine ocp_dynamics()
        type(ocp_model_t) :: model
        type(nbody_system_t) :: system
        real(dp), allocatable :: base(:), state(:), tangent(:)
        real(dp) :: k(3), e0, e1, jvp, plus, minus, eps
        integer :: rejected, step, steps
        logical :: ok
        model%debye_volume = 2; model%box = 4; model%tolerance = 1e-6_dp
        model%dt = 0.01_dp
        call ocp_system(model, system, ok)
        call check(ok .and. model%electrons == 128, 'OCP system in plasma units')
        allocate (base(6*model%particles), state(6*model%particles))
        allocate (tangent(6*model%particles))
        call ocp_prepare(model, 5_int64, base, rejected, status)
        call check(status%code == FORTNUM_OK .and. &
                   abs(sum(base(3*128 + 1::3))) < 1e-12_dp, 'prepared, zero momentum')
        state = base
        call nbody_energy(system, state, e0, status)
        steps = nint(5/model%dt)
        do step = 1, steps
            call ocp_advance(system, state, model%dt, status)
            call check(status%code == FORTNUM_OK, 'OCP step')
        end do
        call nbody_energy(system, state, e1, status)
        print *, 'energy defect per particle', (e1 - e0)/model%electrons, ' T', &
            ocp_kinetic_temperature(model, state)
        call check(abs(e1 - e0) < 1e-3_dp*model%electrons, 'energy over 5/omega_p')
        k = [2*pi/model%box, 0.0_dp, 0.0_dp]
        state = base; tangent = 0
        call ocp_displace_tangent(model, state, [.true., .false., .false.], tangent)
        steps = nint(2/model%dt)
        do step = 1, steps
            call ocp_advance_jvp(system, state, tangent, model%dt, status)
            call check(status%code == FORTNUM_OK, 'OCP tangent step')
        end do
        jvp = mode_tangent(model, state, tangent, k)
        eps = 1.0e-6_dp
        plus = run(model, system, base, steps, k, eps)
        minus = run(model, system, base, steps, k, -eps)
        print *, 'jvp', jvp, ' fd', (plus - minus)/(2*eps)
        call check(abs((plus - minus)/(2*eps) - jvp) < 1e-5_dp*abs(jvp), &
                   'd Re rho_k(t)/d eps: tangent vs central difference')
    end subroutine

    real(dp) function run(model, system, base, steps, k, amplitude) result(mode)
        type(ocp_model_t), intent(in) :: model
        type(nbody_system_t), intent(inout) :: system
        real(dp), intent(in) :: base(:), k(3), amplitude
        integer, intent(in) :: steps
        real(dp) :: state(size(base))
        integer :: s
        state = base
        call ocp_displace(model, state, amplitude, [.true., .false., .false.])
        do s = 1, steps
            call ocp_advance(system, state, model%dt, status)
        end do
        mode = real(ocp_mode(model, state, k))
    end function

    real(dp) function mode_tangent(model, x, dx, kv) result(d)
        type(ocp_model_t), intent(in) :: model
        real(dp), intent(in) :: x(:), dx(:), kv(3)
        integer :: i
        d = 0
        do i = 1, model%electrons
            ! d/dx cos(k.x) = -sin(k.x) k
            d = d - sin(dot_product(kv, x(3*i - 2:3*i)))*dot_product(kv, dx(3*i - 2:3*i))
        end do
    end function
end program test_ocp
