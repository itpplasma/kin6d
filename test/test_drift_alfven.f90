!> Two-dimensional drift-Alfven model (stage DA, homogeneous step): single
!> modes against the three roots of the drift-Alfven dispersion relation,
!> conservation in the ideal limit, the energy and enstrophy balances of a
!> forced nonlinear run, the Hasegawa-Wakatani limit beta, mu -> 0 of a
!> nonlinear run, the E x B advection of the electron inertia, and the
!> tangent and adjoint of the right-hand side.
program test_drift_alfven
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_matrix_exponential, only: matrix_exponential
    use milca_drift_dispersion, only: hw_parameters_t, drift_alfven_parameters_t, &
        drift_alfven_roots, drift_root_ok
    use milca_drift_spectral, only: drift_to_spectral, drift_mean_product
    use milca_drift_fluid_leaves, only: drift_hyperviscosity
    use milca_hasegawa_wakatani, only: hw_state_t, hw_initialize, hw_step
    use milca_drift_alfven
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp)
    integer :: failures

    failures = 0
    call linear_modes()
    call ideal_invariants()
    call balances()
    call hasegawa_wakatani_limit()
    call inertia_advection()
    call derivatives()
    if (failures > 0) error stop 'test_drift_alfven failed'
    write (*, '(a)') 'test_drift_alfven passed'
contains
    !> One seeded eigenmode of each of the three roots follows
    !> exp(-i omega t - nu t), in the electromagnetic, inertial and nearly
    !> resistive (stiff) regimes.
    subroutine linear_modes()
        type(da_parameters_t) :: par
        real(dp) :: error(3), mutant(3)

        par = da_parameters_t(beta=1.0_dp, mu=0.05_dp, c=0.5_dp, kpar=0.8_dp, &
                              kappa=1.0_dp, d=1e-4_dp, order=2)
        call mode_error(par, 3, 5, .false., error)
        call check('electromagnetic regime, three roots', maxval(error), 1e-10_dp)
        ! Mutant: a propagator without the term -beta kappa d a/dy of Ohm's law.
        call mode_error(par, 3, 5, .true., mutant)
        call check_fails('propagator without beta d a/dy is rejected', maxval(mutant), &
                         1e-3_dp)
        par = da_parameters_t(beta=0.1_dp, mu=1.0_dp, c=0.01_dp, kpar=2.0_dp, &
                              kappa=1.0_dp, d=1e-6_dp, order=3)
        call mode_error(par, 7, 1, .false., error)
        call check('inertial regime, three roots', maxval(error), 1e-10_dp)
        par = da_parameters_t(beta=1e-4_dp, mu=1e-4_dp, c=1.0_dp, kpar=1.0_dp, &
                              kappa=1.0_dp, d=1e-4_dp, order=2)
        call mode_error(par, 0, 2, .false., error)
        call check('nearly resistive regime, kx = 0, three roots', maxval(error), 1e-10_dp)
    end subroutine

    !> Error of exp(-i omega T - nu T) for the three roots of
    !> drift_alfven_roots, absolute for a decayed and relative for a grown mode.
    subroutine mode_error(par, mx, my, mutate, error)
        type(da_parameters_t), intent(in) :: par
        integer, intent(in) :: mx, my
        logical, intent(in) :: mutate
        real(dp), intent(out) :: error(3)
        integer, parameter :: n = 32, steps = 40
        real(dp), parameter :: box = 2*pi/0.15_dp, amplitude = 1e-9_dp
        type(da_state_t) :: state
        complex(dp) :: omega(3), v0(3), l(3, 3), factor
        real(dp) :: kx, ky, k2, residual, nu, ind
        integer :: root, status, step, i, j, f, info

        kx = 0.15_dp*mx; ky = 0.15_dp*my; k2 = kx**2 + ky**2
        call drift_alfven_roots(drift_alfven_parameters_t(par%beta, par%mu, par%c), &
                                par%kpar, kx, ky, omega, residual, status)
        if (status /= drift_root_ok .or. residual > 1e-12_dp) failures = failures + 1
        call drift_hyperviscosity(par%d, real(par%order, dp), kx, ky, nu)
        ind = par%beta/k2 + par%mu
        i = mx + 1; j = my + 1
        do root = 1, 3
            call da_initialize(state, par, n, n, box, box, 0.05_dp)
            if (mutate) then
                l = da_linear_operator(par, kx, ky)
                l(3, 3) = real(l(3, 3), dp)
                do f = 1, 3
                    l(f, f) = l(f, f) - nu
                end do
                call matrix_exponential(l*(state%dt/2), state%half(:, :, i, j), info)
            end if
            ! phi = 1: zeta = -k^2, n = (kappa ky - omega k^2)/omega,
            ! j = i omega k^2/kpar, psi = inductance j.
            v0 = amplitude*[cmplx(-k2, 0, dp), (par%kappa*ky - omega(root)*k2)/omega(root), &
                            ind*(0, 1)*omega(root)*k2/par%kpar]
            call set_mode(state, i, j, v0)
            do step = 1, steps
                call da_step(state)
            end do
            factor = exp(-(0, 1)*omega(root)*state%time - nu*state%time)
            error(root) = (abs(state%zeta(i, j)/v0(1) - factor) &
                           + abs(state%n(i, j)/v0(2) - factor) &
                           + abs(state%psi(i, j)/v0(3) - factor))/max(1.0_dp, abs(factor))
        end do
    end subroutine

    subroutine set_mode(state, i, j, v)
        type(da_state_t), intent(inout) :: state
        integer, intent(in) :: i, j
        complex(dp), intent(in) :: v(3)
        integer :: ii, jj
        ii = state%grid%mx(i); jj = state%grid%my(j)
        state%zeta(i, j) = v(1); state%n(i, j) = v(2); state%psi(i, j) = v(3)
        state%zeta(ii, jj) = conjg(v(1)); state%n(ii, jj) = conjg(v(2))
        state%psi(ii, jj) = conjg(v(3))
    end subroutine

    !> kappa = c = d = 0: the linear propagator conserves energy and
    !> potential enstrophy exactly; they change only by the error of the
    !> third-order bracket substep.
    subroutine ideal_invariants()
        type(da_parameters_t) :: par
        type(da_state_t) :: state
        type(da_ledger_t) :: first, last
        real(dp) :: drift(2, 2), dt
        integer :: k, step

        par = da_parameters_t(beta=0.5_dp, mu=0.2_dp, c=0.0_dp, kpar=0.7_dp, &
                              kappa=0.0_dp, d=0.0_dp, order=2)
        do k = 1, 2
            dt = 0.1_dp/k
            call da_initialize(state, par, 32, 32, 2*pi, 2*pi, dt)
            call seed(state, 0.3_dp)
            call da_ledger(state, first)
            do step = 1, nint(4.0_dp/dt)
                call da_step(state)
            end do
            call da_ledger(state, last)
            drift(:, k) = abs([last%energy - first%energy, &
                               last%enstrophy - first%enstrophy]) &
                          /[first%energy, first%enstrophy]
        end do
        write (*, '(a,4es10.2)') 'ideal drift of E and W at dt, dt/2: ', drift
        call check('energy conserved', drift(1, 2), 1e-5_dp)
        call check('potential enstrophy conserved', drift(2, 2), 1e-5_dp)
        call check('third-order convergence of the energy drift', &
                   drift(1, 2)/drift(1, 1), 1/6.0_dp)
        call check('third-order convergence of the enstrophy drift', &
                   drift(2, 2)/drift(2, 1), 1/6.0_dp)
        if (.not. all(drift(:, 1) > 1e-9_dp)) failures = failures + 1
        if (.not. last%magnetic > 1e-3_dp*last%energy) failures = failures + 1
    end subroutine

    !> Centred time derivative of E and W against the ledger of a driven,
    !> dissipative, electromagnetic nonlinear state; the residual falls as dt^2.
    subroutine balances()
        type(da_parameters_t) :: par
        real(dp) :: residual(2, 2)
        integer :: k

        par = da_parameters_t(beta=0.3_dp, mu=0.1_dp, c=0.5_dp, kpar=1.0_dp, &
                              kappa=1.0_dp, d=2e-3_dp, order=2)
        do k = 1, 2
            call balance_residual(par, 0.04_dp/k, residual(:, k))
        end do
        write (*, '(a,4es10.2)') 'balance residuals: ', residual
        call check('energy balance closes', residual(1, 2), 2e-3_dp)
        call check('enstrophy balance closes', residual(2, 2), 2e-3_dp)
        call check('second-order convergence of the energy balance', &
                   residual(1, 2)/residual(1, 1), 1/3.0_dp)
        call check('second-order convergence of the enstrophy balance', &
                   residual(2, 2)/residual(2, 1), 1/3.0_dp)
    end subroutine

    subroutine balance_residual(par, dt, residual)
        type(da_parameters_t), intent(in) :: par
        real(dp), intent(in) :: dt
        real(dp), intent(out) :: residual(2)
        type(da_state_t) :: state
        type(da_ledger_t) :: before, at, after
        real(dp) :: scale
        integer :: step

        call da_initialize(state, par, 32, 32, 2*pi/0.3_dp, 2*pi/0.3_dp, dt)
        call seed(state, 1.5_dp)
        do step = 1, nint(2.0_dp/dt) - 1
            call da_step(state)
        end do
        call da_ledger(state, before)
        call da_step(state)
        call da_ledger(state, at)
        call da_step(state)
        call da_ledger(state, after)
        associate (g => state%grid)
            if (any(state%psi /= conjg(state%psi(g%mx, g%my)))) then
                write (*, '(a)') 'FAIL mirror symmetry of the coefficients'
                failures = failures + 1
            end if
        end associate
        scale = abs(at%flux) + at%d_c + at%d_e
        if (.not. (at%d_c > 0.05_dp*scale .and. abs(at%flux) > 0.05_dp*scale &
                   .and. at%d_e > 1e-3_dp*scale .and. at%magnetic > 1e-3_dp*at%energy)) &
            failures = failures + 1
        residual(1) = abs((after%energy - before%energy)/(2*dt) &
                          - (at%flux - at%d_c - at%d_e))/scale
        residual(2) = abs((after%enstrophy - before%enstrophy)/(2*dt) &
                          - (at%flux - at%d_w))/scale
    end subroutine

    !> beta = mu = eps -> 0 at fixed c: j -> -kpar (n - phi)/c and the model
    !> becomes Hasegawa-Wakatani with alpha = kpar^2/c. A nonlinear run that
    !> starts on that slow manifold approaches the Hasegawa-Wakatani run at the
    !> same time step with an error of order eps.
    subroutine hasegawa_wakatani_limit()
        type(da_parameters_t) :: par
        type(hw_parameters_t) :: hw
        type(da_state_t) :: state
        type(hw_state_t) :: reference
        complex(dp), allocatable :: phi(:, :)
        real(dp), parameter :: dt = 0.02_dp, box = 2*pi/0.3_dp
        real(dp) :: difference(3), eps
        integer :: k, step

        hw%alpha = 1; hw%kappa = 1; hw%d = 2e-3_dp; hw%order = 2
        call hw_initialize(reference, hw, 32, 32, box, box, dt)
        do k = 1, 3
            eps = 10.0_dp**(-k - 1)
            par = da_parameters_t(beta=eps, mu=eps, c=1.0_dp, kpar=1.0_dp, kappa=1.0_dp, &
                                  d=hw%d, order=hw%order)
            call da_initialize(state, par, 32, 32, box, box, dt)
            call seed(state, 1.5_dp)
            allocate (phi, mold=state%zeta)
            phi = merge(-state%zeta/state%grid%k2, (0.0_dp, 0.0_dp), state%grid%k2 > 0)
            state%psi = -par%kpar*(state%n - phi)/par%c
            state%psi = merge(state%psi/state%admittance, (0.0_dp, 0.0_dp), &
                              state%admittance > 0)
            deallocate (phi)
            if (k == 1) then
                reference%zeta = state%zeta; reference%n = state%n
                do step = 1, nint(2.0_dp/dt)
                    call hw_step(reference)
                end do
            end if
            do step = 1, nint(2.0_dp/dt)
                call da_step(state)
            end do
            difference(k) = (maxval(abs(state%zeta - reference%zeta)) &
                             + maxval(abs(state%n - reference%n))) &
                            /(maxval(abs(reference%zeta)) + maxval(abs(reference%n)))
        end do
        write (*, '(a,3es10.2)') 'distance to Hasegawa-Wakatani at eps = 1e-2, 1e-3, 1e-4: ', &
            difference
        call check('Hasegawa-Wakatani limit', difference(3), 1e-3_dp)
        call check('first-order approach to the limit', difference(3)/difference(2), 0.2_dp)
        call check('first-order approach to the limit (coarser)', &
                   difference(2)/difference(1), 0.2_dp)
    end subroutine

    !> Without linear terms (kpar = c = kappa = d = 0) the rate of psi is
    !> -mu {phi, j}: phi = a sin x and j = b cos y give mu a b cos x sin y.
    !> The electron inertia is energy-neutral, so the other tests do not see it.
    subroutine inertia_advection()
        type(da_parameters_t) :: par
        type(da_state_t) :: state
        complex(dp), allocatable :: rate(:, :, :), expected(:, :)
        real(dp), allocatable :: f(:, :)
        real(dp), parameter :: a = 0.7_dp, b = -1.3_dp
        real(dp) :: x, y
        integer :: i, k, n

        n = 16
        par = da_parameters_t(beta=0.4_dp, mu=0.3_dp, c=0.0_dp, kpar=0.0_dp, &
                              kappa=0.0_dp, d=0.0_dp, order=2)
        call da_initialize(state, par, n, n, 2*pi, 2*pi, 0.01_dp)
        allocate (rate(n, n, 3), expected(n, n), f(n, n))
        associate (g => state%grid)
            do k = 1, n
                do i = 1, n
                    x = 2*pi*(i - 1)/n; y = 2*pi*(k - 1)/n
                    f(i, k) = -a*sin(x)
                end do
            end do
            call drift_to_spectral(g, state%work, f, state%zeta)
            do k = 1, n
                do i = 1, n
                    y = 2*pi*(k - 1)/n
                    f(i, k) = (par%beta + par%mu)*b*cos(y)
                end do
            end do
            call drift_to_spectral(g, state%work, f, state%psi)
            do k = 1, n
                do i = 1, n
                    x = 2*pi*(i - 1)/n; y = 2*pi*(k - 1)/n
                    f(i, k) = par%mu*a*b*cos(x)*sin(y)
                end do
            end do
            call drift_to_spectral(g, state%work, f, expected)
        end associate
        call da_rhs(state, state%zeta, state%n, state%psi, rate)
        call check('inertia advected by E x B', maxval(abs(rate(:, :, 3) - expected)) &
                   /maxval(abs(expected)), 1e-13_dp)
        call check('vorticity and density at rest', maxval(abs(rate(:, :, 1:2))), 1e-14_dp)
    end subroutine

    subroutine derivatives()
        type(da_parameters_t) :: par
        type(da_state_t) :: state, tangent, cotangent
        complex(dp), allocatable :: d(:, :, :), p(:, :, :), m(:, :, :), g(:, :, :), b(:, :, :)
        real(dp), parameter :: h = 1e-5_dp
        real(dp) :: left, right
        integer :: n, f

        n = 24
        par = da_parameters_t(beta=0.4_dp, mu=0.3_dp, c=0.2_dp, kpar=0.9_dp, &
                              kappa=1.2_dp, d=1e-3_dp, order=2)
        call da_initialize(state, par, n, n, 9.0_dp, 7.0_dp, 0.01_dp)
        call da_initialize(tangent, par, n, n, 9.0_dp, 7.0_dp, 0.01_dp)
        call da_initialize(cotangent, par, n, n, 9.0_dp, 7.0_dp, 0.01_dp)
        call seed(state, 1.0_dp); call seed(tangent, 0.7_dp, 5)
        call seed(cotangent, 0.4_dp, 9)
        allocate (d(n, n, 3), p(n, n, 3), m(n, n, 3), g(n, n, 3), b(n, n, 3))
        associate (s => state, t => tangent)
            call da_rhs_jvp(s, s%zeta, s%n, s%psi, t%zeta, t%n, t%psi, d)
            call da_rhs(s, s%zeta + h*t%zeta, s%n + h*t%n, s%psi + h*t%psi, p)
            call da_rhs(s, s%zeta - h*t%zeta, s%n - h*t%n, s%psi - h*t%psi, m)
        end associate
        call check('tangent vs central difference', &
                   maxval(abs(d - (p - m)/(2*h)))/maxval(abs(d)), 1e-9_dp)
        b(:, :, 1) = cotangent%zeta; b(:, :, 2) = cotangent%n; b(:, :, 3) = cotangent%psi
        call da_rhs_vjp(state, state%zeta, state%n, state%psi, b, g)
        left = 0; right = 0
        do f = 1, 3
            left = left + drift_mean_product(b(:, :, f), d(:, :, f))
        end do
        right = drift_mean_product(g(:, :, 1), tangent%zeta) &
                + drift_mean_product(g(:, :, 2), tangent%n) &
                + drift_mean_product(g(:, :, 3), tangent%psi)
        if (.not. abs(left) > 1e-6_dp) failures = failures + 1
        call check('adjoint identity', abs(left - right)/abs(left), 1e-12_dp)
    end subroutine

    !> Smooth reproducible fields with many interacting modes; psi from a
    !> current of comparable size.
    subroutine seed(state, amplitude, shift)
        type(da_state_t), intent(inout) :: state
        real(dp), intent(in) :: amplitude
        integer, intent(in), optional :: shift
        real(dp), allocatable :: zeta(:, :), n(:, :), j(:, :)
        complex(dp), allocatable :: jh(:, :)
        real(dp) :: x, y, s
        integer :: i, k, m

        s = 0
        if (present(shift)) s = shift
        associate (g => state%grid)
            allocate (zeta(g%nx, g%ny), n(g%nx, g%ny), j(g%nx, g%ny), jh(g%nx, g%ny))
            zeta = 0; n = 0; j = 0
            do k = 1, g%ny
                do i = 1, g%nx
                    x = 2*pi*(i - 1)/g%nx; y = 2*pi*(k - 1)/g%ny
                    do m = 1, 4
                        zeta(i, k) = zeta(i, k) + amplitude*sin(m*x + 0.7_dp*m + s) &
                                     *cos((5 - m)*y + 0.3_dp*m*m - s)/m
                        n(i, k) = n(i, k) + amplitude*cos(m*x - 1.1_dp*m + s) &
                                  *cos(m*y + 0.2_dp*m + 2*s) + amplitude*sin((m + 1)*y - x + s)/m
                        j(i, k) = j(i, k) + amplitude*cos((m + 1)*x + 0.4_dp*m - s) &
                                  *sin(m*y - 0.9_dp*m + s)/m
                    end do
                end do
            end do
            call drift_to_spectral(g, state%work, zeta, state%zeta)
            call drift_to_spectral(g, state%work, n, state%n)
            call drift_to_spectral(g, state%work, j, jh)
            state%psi = merge(jh/state%admittance, (0.0_dp, 0.0_dp), state%admittance > 0)
        end associate
    end subroutine

    subroutine check(name, value, tolerance)
        character(*), intent(in) :: name
        real(dp), intent(in) :: value, tolerance
        if (value <= tolerance) then
            write (*, '(a,a,es10.3)') 'ok   ', name//': ', value
        else
            write (*, '(a,a,es10.3,a,es10.3)') 'FAIL ', name//': ', value, ' > ', tolerance
            failures = failures + 1
        end if
    end subroutine

    subroutine check_fails(name, value, tolerance)
        character(*), intent(in) :: name
        real(dp), intent(in) :: value, tolerance
        if (value > tolerance) then
            write (*, '(a,a,es10.3)') 'ok   ', name//': ', value
        else
            write (*, '(a,a,es10.3)') 'FAIL ', name//': ', value
            failures = failures + 1
        end if
    end subroutine
end program test_drift_alfven
