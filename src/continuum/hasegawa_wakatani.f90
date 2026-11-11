!> Hasegawa-Wakatani drift-wave turbulence on a doubly periodic box, original
!> and modified model (docs/drift-turbulence.md).
!>
!> State: Fourier coefficients of the vorticity zeta = laplace(phi) and of the
!> density n. One step is the Strang composition of the exact linear
!> propagator over dt/2, the bracket terms over dt with the third-order
!> strong-stability-preserving Runge-Kutta scheme, and the propagator again.
!> The linear part (coupling, drive, hyperviscosity) is therefore free of
!> time-step error for a single mode. The state stays the transform of real
!> fields: its mirror symmetry is restored after every step.
module milca_hasegawa_wakatani
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use milca_drift_dispersion, only: hw_parameters_t
    use milca_drift_spectral, only: drift_grid_t, drift_work_t, drift_grid_init, &
        drift_work_init, drift_bracket, drift_poisson_inverse, drift_enforce_reality
    use milca_drift_fluid_leaves, only: drift_hyperviscosity, hw_linear_matrix, &
        hw_propagator, hw_propagator_decoupled, hw_mode_ledger
    implicit none
    private
    public :: hw_state_t, hw_ledger_t, hw_initialize, hw_step, hw_ledger
    public :: hw_rhs, hw_rhs_jvp, hw_rhs_vjp

    type :: hw_state_t
        type(hw_parameters_t) :: par
        type(drift_grid_t) :: grid
        type(drift_work_t) :: work
        real(dp) :: dt = 0, time = 0
        complex(dp), allocatable :: zeta(:, :), n(:, :)
        !> Propagator over dt/2, entries (1,1), (1,2), (2,1), (2,2) per mode.
        complex(dp), allocatable :: half(:, :, :)
        complex(dp), allocatable :: phi(:, :), stage(:, :, :), rate(:, :, :)
    end type

    !> Box averages; d energy/dt = flux - d_alpha - d_e and
    !> d enstrophy/dt = flux - d_w. kinetic = <|grad phi|^2>/2 and its part
    !> zonal_kinetic in the modes ky = 0.
    type :: hw_ledger_t
        real(dp) :: energy = 0, enstrophy = 0, flux = 0
        real(dp) :: d_alpha = 0, d_e = 0, d_w = 0
        real(dp) :: kinetic = 0, zonal_kinetic = 0
    end type
contains
    subroutine hw_initialize(state, par, nx, ny, lx, ly, dt)
        type(hw_state_t), intent(out) :: state
        type(hw_parameters_t), intent(in) :: par
        integer, intent(in) :: nx, ny
        real(dp), intent(in) :: lx, ly, dt
        real(dp) :: nu, e(8)
        integer :: i, j

        if (par%alpha < 0 .or. par%d < 0 .or. par%order < 0) &
            error stop 'hw_initialize: need alpha, d, order >= 0'
        if (.not. dt > 0) error stop 'hw_initialize: need dt > 0'
        state%par = par
        state%dt = dt
        call drift_grid_init(state%grid, nx, ny, lx, ly)
        call drift_work_init(state%work, state%grid)
        allocate (state%zeta(nx, ny), state%n(nx, ny), state%phi(nx, ny))
        allocate (state%half(4, nx, ny), state%stage(nx, ny, 2), state%rate(nx, ny, 2))
        state%zeta = 0; state%n = 0
        associate (g => state%grid)
            do j = 1, ny
                do i = 1, nx
                    state%half(:, i, j) = [1, 0, 0, 1]
                    if (.not. g%keep(i, j) .or. .not. g%k2(i, j) > 0) cycle
                    call drift_hyperviscosity(par%d, real(par%order, dp), g%kx(i), &
                                              g%ky(j), nu)
                    if (coupling(par, g%ky(j)) > 0) then
                        call hw_propagator(par%alpha, par%kappa, g%kx(i), g%ky(j), nu, &
                            dt/2, e(1), e(2), e(3), e(4), e(5), e(6), e(7), e(8))
                        state%half(:, i, j) = cmplx(e(1:7:2), e(2:8:2), dp)
                    else
                        call hw_propagator_decoupled(par%kappa, g%kx(i), g%ky(j), nu, &
                                                     dt/2, e(1), e(2))
                        state%half(:, i, j) = [cmplx(e(1), 0, dp), (0.0_dp, 0.0_dp), &
                                               cmplx(0, e(2), dp), cmplx(e(1), 0, dp)]
                    end if
                end do
            end do
        end associate
    end subroutine

    !> Adiabaticity acting on a mode: zero for ky = 0 in the modified model.
    pure real(dp) function coupling(par, ky)
        type(hw_parameters_t), intent(in) :: par
        real(dp), intent(in) :: ky
        coupling = par%alpha
        if (par%modified .and. .not. abs(ky) > 0) coupling = 0
    end function

    subroutine hw_step(state)
        type(hw_state_t), intent(inout) :: state
        real(dp) :: dt

        dt = state%dt
        call propagate_half(state)
        associate (u => state%stage, k => state%rate, z => state%zeta, n => state%n)
            call bracket_terms(state, z, n, k(:, :, 1), k(:, :, 2))
            u(:, :, 1) = z + dt*k(:, :, 1)
            u(:, :, 2) = n + dt*k(:, :, 2)
            call bracket_terms(state, u(:, :, 1), u(:, :, 2), k(:, :, 1), k(:, :, 2))
            u(:, :, 1) = 0.75_dp*z + 0.25_dp*(u(:, :, 1) + dt*k(:, :, 1))
            u(:, :, 2) = 0.75_dp*n + 0.25_dp*(u(:, :, 2) + dt*k(:, :, 2))
            call bracket_terms(state, u(:, :, 1), u(:, :, 2), k(:, :, 1), k(:, :, 2))
            z = (z + 2*(u(:, :, 1) + dt*k(:, :, 1)))/3
            n = (n + 2*(u(:, :, 2) + dt*k(:, :, 2)))/3
        end associate
        call propagate_half(state)
        call drift_enforce_reality(state%grid, state%zeta)
        call drift_enforce_reality(state%grid, state%n)
        state%time = state%time + dt
    end subroutine

    subroutine propagate_half(state)
        type(hw_state_t), intent(inout) :: state
        complex(dp) :: z
        integer :: i, j

        !$omp parallel do private(i, j, z) schedule(static)
        do j = 1, state%grid%ny
            do i = 1, state%grid%nx
                z = state%zeta(i, j)
                state%zeta(i, j) = state%half(1, i, j)*z + state%half(2, i, j)*state%n(i, j)
                state%n(i, j) = state%half(3, i, j)*z + state%half(4, i, j)*state%n(i, j)
            end do
        end do
        !$omp end parallel do
    end subroutine

    !> (-{phi, zeta}, -{phi, n}) with laplace(phi) = zeta.
    subroutine bracket_terms(state, zeta, n, dzeta, dn)
        type(hw_state_t), intent(inout) :: state
        complex(dp), intent(in) :: zeta(:, :), n(:, :)
        complex(dp), intent(out) :: dzeta(:, :), dn(:, :)

        call drift_poisson_inverse(state%grid, zeta, state%phi)
        call drift_bracket(state%grid, state%work, state%phi, zeta, dzeta)
        call drift_bracket(state%grid, state%work, state%phi, n, dn)
        dzeta = -dzeta
        dn = -dn
    end subroutine

    !> Full right-hand side d/dt (zeta, n) at the coefficients (zeta, n).
    subroutine hw_rhs(state, zeta, n, dzeta, dn)
        type(hw_state_t), intent(inout) :: state
        complex(dp), intent(in) :: zeta(:, :), n(:, :)
        complex(dp), intent(out) :: dzeta(:, :), dn(:, :)

        call bracket_terms(state, zeta, n, dzeta, dn)
        call add_linear(state, zeta, n, dzeta, dn, adjoint=.false.)
    end subroutine

    !> Tangent of hw_rhs at (zeta, n) along (tzeta, tn).
    subroutine hw_rhs_jvp(state, zeta, n, tzeta, tn, dzeta, dn)
        type(hw_state_t), intent(inout) :: state
        complex(dp), intent(in) :: zeta(:, :), n(:, :), tzeta(:, :), tn(:, :)
        complex(dp), intent(out) :: dzeta(:, :), dn(:, :)
        complex(dp), allocatable :: phi(:, :), tphi(:, :), part(:, :)

        allocate (phi, mold=zeta); allocate (tphi, mold=zeta); allocate (part, mold=zeta)
        call drift_poisson_inverse(state%grid, zeta, phi)
        call drift_poisson_inverse(state%grid, tzeta, tphi)
        call drift_bracket(state%grid, state%work, tphi, zeta, dzeta)
        call drift_bracket(state%grid, state%work, phi, tzeta, part)
        dzeta = -dzeta - part
        call drift_bracket(state%grid, state%work, tphi, n, dn)
        call drift_bracket(state%grid, state%work, phi, tn, part)
        dn = -dn - part
        call add_linear(state, tzeta, tn, dzeta, dn, adjoint=.false.)
    end subroutine

    !> Adjoint of hw_rhs_jvp under the box average: cotangent (bzeta, bn) of
    !> the rates to the cotangent (gzeta, gn) of the state.
    subroutine hw_rhs_vjp(state, zeta, n, bzeta, bn, gzeta, gn)
        type(hw_state_t), intent(inout) :: state
        complex(dp), intent(in) :: zeta(:, :), n(:, :), bzeta(:, :), bn(:, :)
        complex(dp), intent(out) :: gzeta(:, :), gn(:, :)
        complex(dp), allocatable :: phi(:, :), part(:, :), source(:, :)

        allocate (phi, mold=zeta); allocate (part, mold=zeta); allocate (source, mold=zeta)
        call drift_poisson_inverse(state%grid, zeta, phi)
        ! <b, -{phi_t, f}> = -<phi_t, {f, b}> and <b, -{phi, f_t}> = <f_t, {phi, b}>.
        call drift_bracket(state%grid, state%work, zeta, bzeta, source)
        call drift_bracket(state%grid, state%work, n, bn, part)
        source = -source - part
        call drift_poisson_inverse(state%grid, source, gzeta)
        call drift_bracket(state%grid, state%work, phi, bzeta, part)
        gzeta = gzeta + part
        call drift_bracket(state%grid, state%work, phi, bn, gn)
        call add_linear(state, bzeta, bn, gzeta, gn, adjoint=.true.)
    end subroutine

    !> (dzeta, dn) += (M - nu) (zeta, n), or its conjugate transpose.
    subroutine add_linear(state, zeta, n, dzeta, dn, adjoint)
        type(hw_state_t), intent(in) :: state
        complex(dp), intent(in) :: zeta(:, :), n(:, :)
        complex(dp), intent(inout) :: dzeta(:, :), dn(:, :)
        logical, intent(in) :: adjoint
        real(dp) :: m11, m12, m21r, m21i, m22, nu
        complex(dp) :: m21
        integer :: i, j

        associate (g => state%grid, par => state%par)
            !$omp parallel do private(i, j, m11, m12, m21r, m21i, m22, nu, m21) &
            !$omp schedule(static)
            do j = 1, g%ny
                do i = 1, g%nx
                    if (.not. g%keep(i, j) .or. .not. g%k2(i, j) > 0) cycle
                    call drift_hyperviscosity(par%d, real(par%order, dp), g%kx(i), &
                                              g%ky(j), nu)
                    call hw_linear_matrix(coupling(par, g%ky(j)), par%kappa, g%kx(i), &
                                          g%ky(j), m11, m12, m21r, m21i, m22)
                    m21 = cmplx(m21r, m21i, dp)
                    if (adjoint) then
                        dzeta(i, j) = dzeta(i, j) + (m11 - nu)*zeta(i, j) + conjg(m21)*n(i, j)
                        dn(i, j) = dn(i, j) + m12*zeta(i, j) + (m22 - nu)*n(i, j)
                    else
                        dzeta(i, j) = dzeta(i, j) + (m11 - nu)*zeta(i, j) + m12*n(i, j)
                        dn(i, j) = dn(i, j) + m21*zeta(i, j) + (m22 - nu)*n(i, j)
                    end if
                end do
            end do
            !$omp end parallel do
        end associate
    end subroutine

    subroutine hw_ledger(state, ledger)
        type(hw_state_t), intent(in) :: state
        type(hw_ledger_t), intent(out) :: ledger
        real(dp) :: nu, part(6), total(6), kinetic, zonal, ek
        integer :: i, j

        total = 0; kinetic = 0; zonal = 0
        associate (g => state%grid, par => state%par)
            !$omp parallel do private(i, j, nu, part, ek) schedule(static) &
            !$omp reduction(+:total, kinetic, zonal)
            do j = 1, g%ny
                do i = 1, g%nx
                    if (.not. g%keep(i, j) .or. .not. g%k2(i, j) > 0) cycle
                    call drift_hyperviscosity(par%d, real(par%order, dp), g%kx(i), &
                                              g%ky(j), nu)
                    call hw_mode_ledger(coupling(par, g%ky(j)), par%kappa, g%kx(i), &
                        g%ky(j), nu, real(state%zeta(i, j), dp), aimag(state%zeta(i, j)), &
                        real(state%n(i, j), dp), aimag(state%n(i, j)), part(1), part(2), &
                        part(3), part(4), part(5), part(6))
                    total = total + part
                    ek = abs(state%zeta(i, j))**2/(2*g%k2(i, j))
                    kinetic = kinetic + ek
                    if (.not. abs(g%ky(j)) > 0) zonal = zonal + ek
                end do
            end do
            !$omp end parallel do
        end associate
        ledger = hw_ledger_t(total(1), total(2), total(3), total(4), total(5), total(6), &
                             kinetic, zonal)
    end subroutine
end module milca_hasegawa_wakatani
