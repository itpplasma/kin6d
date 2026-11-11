!> Hasegawa-Wakatani on the periodic box: single modes against the dispersion
!> relation, conservation in the ideal limit, the energy and enstrophy
!> balances of a forced nonlinear run with their time-step convergence, the
!> modified model, and the tangent and adjoint of the right-hand side.
program test_hasegawa_wakatani
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use milca_drift_dispersion, only: hw_parameters_t, hw_linear_roots, drift_root_ok
    use milca_drift_spectral, only: drift_to_spectral, drift_mean_product
    use milca_hasegawa_wakatani
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp)
    integer :: failures

    failures = 0
    call linear_modes()
    call ideal_invariants()
    call balances()
    call modified_zonal()
    call derivatives()
    if (failures > 0) error stop 'test_hasegawa_wakatani failed'
    write (*, '(a)') 'test_hasegawa_wakatani passed'
contains
    !> One seeded eigenmode of each root follows exp(-i omega t). The box and
    !> the hyperviscosity are those of Numata, Ball and Dewar (2007).
    subroutine linear_modes()
        type(hw_parameters_t) :: par, wrong
        real(dp) :: error(2), mutant(2)
        integer :: model

        par%alpha = 1; par%kappa = 1; par%d = 1e-4_dp; par%order = 2
        do model = 1, 2
            par%modified = model == 2
            wrong = par
            wrong%kappa = -par%kappa
            call mode_error(par, par, 3, 5, error)
            call check('growing mode vs dispersion relation', error(1), 1e-11_dp)
            call check('damped mode vs dispersion relation', error(2), 1e-11_dp)
            ! Mutant: a propagator built with the opposite sign of kappa.
            call mode_error(wrong, par, 3, 5, mutant)
            call check_fails('propagator with wrong sign of kappa is rejected', &
                             mutant(1), 1e-3_dp)
        end do
        par%alpha = 0.01_dp
        call mode_error(par, par, 0, 2, error)
        call check('hydrodynamic regime, kx = 0', maxval(error), 1e-11_dp)
        par%alpha = 50; par%order = 3; par%d = 1e-6_dp
        call mode_error(par, par, 7, 1, error)
        call check('adiabatic regime, sixth-order hyperviscosity', maxval(error), 1e-11_dp)
    end subroutine

    !> Relative error of exp(-i omega T) for the two roots of `expected`,
    !> measured in a run with the parameters `run`.
    subroutine mode_error(run, expected, mx, my, error)
        type(hw_parameters_t), intent(in) :: run, expected
        integer, intent(in) :: mx, my
        real(dp), intent(out) :: error(2)
        integer, parameter :: n = 32, steps = 40
        real(dp), parameter :: box = 2*pi/0.15_dp, amplitude = 1e-9_dp
        type(hw_state_t) :: state
        complex(dp) :: omega(2), response(2), factor, zeta0
        real(dp) :: kx, ky
        integer :: root, status, step, i, j

        kx = 0.15_dp*mx; ky = 0.15_dp*my
        call hw_linear_roots(expected, kx, ky, omega, response, status)
        if (status /= drift_root_ok) failures = failures + 1
        i = mx + 1; j = my + 1
        do root = 1, 2
            call hw_initialize(state, run, n, n, box, box, 0.05_dp)
            zeta0 = -(kx**2 + ky**2)*amplitude
            call set_mode(state, i, j, zeta0, response(root)*amplitude)
            do step = 1, steps
                call hw_step(state)
            end do
            factor = exp(-(0, 1)*omega(root)*state%time)
            ! Absolute error for a decayed mode, relative for a grown one.
            error(root) = (abs(state%zeta(i, j)/zeta0 - factor) &
                           + abs(state%n(i, j)/(response(root)*amplitude) - factor)) &
                          /max(1.0_dp, abs(factor))
        end do
    end subroutine

    subroutine set_mode(state, i, j, zeta, n)
        type(hw_state_t), intent(inout) :: state
        integer, intent(in) :: i, j
        complex(dp), intent(in) :: zeta, n
        state%zeta(i, j) = zeta
        state%n(i, j) = n
        state%zeta(state%grid%mx(i), state%grid%my(j)) = conjg(zeta)
        state%n(state%grid%mx(i), state%grid%my(j)) = conjg(n)
    end subroutine

    !> alpha = kappa = D = 0: energy and potential enstrophy change only by
    !> the error of the third-order bracket substep.
    subroutine ideal_invariants()
        type(hw_parameters_t) :: par
        type(hw_ledger_t) :: first, last
        real(dp) :: drift(2, 2)
        integer :: k

        par%alpha = 0; par%kappa = 0; par%d = 0
        do k = 1, 2
            call run(par, 0.1_dp/k, 4.0_dp, 0.3_dp, first, last)
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
    end subroutine

    !> Centred time derivative of E and W against the ledger of a driven,
    !> dissipative, nonlinear state; the residual falls as dt^2.
    subroutine balances()
        type(hw_parameters_t) :: par
        real(dp) :: residual(2, 2)
        integer :: k, model

        par%alpha = 0.5_dp; par%kappa = 1; par%d = 2e-3_dp; par%order = 2
        do model = 1, 2
            par%modified = model == 2
            do k = 1, 2
                call balance_residual(par, 0.04_dp/k, residual(:, k))
            end do
            write (*, '(a,l2,4es10.2)') 'balance residuals, modified =', par%modified, residual
            call check('energy balance closes', residual(1, 2), 2e-3_dp)
            call check('enstrophy balance closes', residual(2, 2), 2e-3_dp)
            call check('second-order convergence of the energy balance', &
                       residual(1, 2)/residual(1, 1), 1/3.0_dp)
            call check('second-order convergence of the enstrophy balance', &
                       residual(2, 2)/residual(2, 1), 1/3.0_dp)
        end do
    end subroutine

    subroutine balance_residual(par, dt, residual)
        type(hw_parameters_t), intent(in) :: par
        real(dp), intent(in) :: dt
        real(dp), intent(out) :: residual(2)
        type(hw_state_t) :: state
        type(hw_ledger_t) :: before, at, after
        real(dp) :: scale
        integer :: step, steps

        call hw_initialize(state, par, 32, 32, 2*pi/0.3_dp, 2*pi/0.3_dp, dt)
        call seed(state, 1.5_dp)
        steps = nint(2.0_dp/dt)
        do step = 1, steps - 1
            call hw_step(state)
        end do
        call hw_ledger(state, before)
        call hw_step(state)
        call hw_ledger(state, at)
        call hw_step(state)
        call hw_ledger(state, after)
        scale = abs(at%flux) + at%d_alpha + at%d_e
        if (.not. (at%d_alpha > 0.05_dp*scale .and. abs(at%flux) > 0.05_dp*scale &
                   .and. at%d_e > 1e-3_dp*scale)) failures = failures + 1
        residual(1) = abs((after%energy - before%energy)/(2*dt) &
                          - (at%flux - at%d_alpha - at%d_e))/scale
        residual(2) = abs((after%enstrophy - before%enstrophy)/(2*dt) &
                          - (at%flux - at%d_w))/scale
    end subroutine

    !> Without drive the zonal modes of the modified model decay by
    !> hyperviscosity alone and feel no resistive sink; in the original model
    !> the coupling damps them.
    subroutine modified_zonal()
        type(hw_parameters_t) :: par
        type(hw_state_t) :: state
        type(hw_ledger_t) :: ledger
        complex(dp) :: zeta0
        real(dp) :: nu
        integer :: step

        par%alpha = 2; par%kappa = 1; par%d = 1e-2_dp; par%order = 2
        par%modified = .true.
        call hw_initialize(state, par, 16, 16, 2*pi, 2*pi, 0.05_dp)
        zeta0 = (0.3_dp, -0.2_dp)
        call set_mode(state, 3, 1, zeta0, (0.1_dp, 0.4_dp))
        call hw_ledger(state, ledger)
        call check('no resistive sink on zonal modes (modified)', ledger%d_alpha, 0.0_dp)
        call check('zonal kinetic fraction is one', &
                   abs(ledger%zonal_kinetic/ledger%kinetic - 1), 1e-15_dp)
        do step = 1, 20
            call hw_step(state)
        end do
        nu = par%d*2.0_dp**4
        call check('zonal decay by hyperviscosity alone', &
                   abs(state%zeta(3, 1)/zeta0 - exp(-nu*state%time)), 1e-13_dp)
        par%modified = .false.
        call hw_initialize(state, par, 16, 16, 2*pi, 2*pi, 0.05_dp)
        call set_mode(state, 3, 1, zeta0, (0.1_dp, 0.4_dp))
        call hw_ledger(state, ledger)
        if (.not. ledger%d_alpha > 0.1_dp) failures = failures + 1
    end subroutine

    subroutine derivatives()
        type(hw_parameters_t) :: par
        type(hw_state_t) :: state, tangent, cotangent
        complex(dp), allocatable, dimension(:, :) :: dz, dn, pz, pn, mz, mn, gz, gn
        real(dp), parameter :: h = 1e-5_dp
        real(dp) :: left, right
        integer :: n

        n = 24
        par%alpha = 0.7_dp; par%kappa = 1.2_dp; par%d = 1e-3_dp; par%modified = .true.
        call hw_initialize(state, par, n, n, 9.0_dp, 7.0_dp, 0.01_dp)
        call hw_initialize(tangent, par, n, n, 9.0_dp, 7.0_dp, 0.01_dp)
        call hw_initialize(cotangent, par, n, n, 9.0_dp, 7.0_dp, 0.01_dp)
        call seed(state, 1.0_dp); call seed(tangent, 0.7_dp, 5); call seed(cotangent, 0.4_dp, 9)
        allocate (dz(n, n), dn(n, n), pz(n, n), pn(n, n), mz(n, n), mn(n, n), gz(n, n), &
                  gn(n, n))
        call hw_rhs_jvp(state, state%zeta, state%n, tangent%zeta, tangent%n, dz, dn)
        call hw_rhs(state, state%zeta + h*tangent%zeta, state%n + h*tangent%n, pz, pn)
        call hw_rhs(state, state%zeta - h*tangent%zeta, state%n - h*tangent%n, mz, mn)
        call check('tangent vs central difference', &
                   (maxval(abs(dz - (pz - mz)/(2*h))) + maxval(abs(dn - (pn - mn)/(2*h)))) &
                   /(maxval(abs(dz)) + maxval(abs(dn))), 1e-9_dp)
        call hw_rhs_vjp(state, state%zeta, state%n, cotangent%zeta, cotangent%n, gz, gn)
        left = drift_mean_product(cotangent%zeta, dz) + drift_mean_product(cotangent%n, dn)
        right = drift_mean_product(gz, tangent%zeta) + drift_mean_product(gn, tangent%n)
        if (.not. abs(left) > 1e-6_dp) failures = failures + 1
        call check('adjoint identity', abs(left - right)/abs(left), 1e-12_dp)
    end subroutine

    subroutine run(par, dt, time, amplitude, first, last)
        type(hw_parameters_t), intent(in) :: par
        real(dp), intent(in) :: dt, time, amplitude
        type(hw_ledger_t), intent(out) :: first, last
        type(hw_state_t) :: state
        integer :: step

        call hw_initialize(state, par, 32, 32, 2*pi, 2*pi, dt)
        call seed(state, amplitude)
        call hw_ledger(state, first)
        do step = 1, nint(time/dt)
            call hw_step(state)
        end do
        call hw_ledger(state, last)
    end subroutine

    !> Smooth reproducible fields with many interacting modes.
    subroutine seed(state, amplitude, shift)
        type(hw_state_t), intent(inout) :: state
        real(dp), intent(in) :: amplitude
        integer, intent(in), optional :: shift
        real(dp), allocatable :: zeta(:, :), n(:, :)
        real(dp) :: x, y, s
        integer :: i, j, m

        s = 0
        if (present(shift)) s = shift
        associate (g => state%grid)
            allocate (zeta(g%nx, g%ny), n(g%nx, g%ny))
            zeta = 0; n = 0
            do j = 1, g%ny
                do i = 1, g%nx
                    x = 2*pi*(i - 1)/g%nx; y = 2*pi*(j - 1)/g%ny
                    do m = 1, 4
                        zeta(i, j) = zeta(i, j) + amplitude*sin(m*x + 0.7_dp*m + s) &
                                     *cos((5 - m)*y + 0.3_dp*m*m - s)/m
                        n(i, j) = n(i, j) + amplitude*cos(m*x - 1.1_dp*m + s) &
                                  *cos(m*y + 0.2_dp*m + 2*s) + amplitude*sin((m + 1)*y - x + s)/m
                    end do
                end do
            end do
            call drift_to_spectral(g, state%work, zeta, state%zeta)
            call drift_to_spectral(g, state%work, n, state%n)
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
end program test_hasegawa_wakatani
