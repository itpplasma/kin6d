!> Two-dimensional electromagnetic drift-Alfven model with one parallel wave
!> number on a doubly periodic box: first step of stage DA in
!> docs/drift-turbulence.md.
!>
!> State: Fourier coefficients of the vorticity zeta = laplace(phi), the
!> density n and the generalised flux psi = beta a + mu j, j = -laplace(a).
!> One step is the Strang composition of the exact linear propagator of each
!> mode over dt/2 (coupling, drive, Alfven and inertial response, resistivity
!> and hyperviscosity), the bracket terms over dt with the third-order
!> strong-stability-preserving Runge-Kutta scheme, and the propagator again.
!> The linear block is stiff when beta and mu are small; its exponential has
!> no stability limit, so the Hasegawa-Wakatani limit beta, mu -> 0 can be
!> approached at the time step of that model.
module milca_drift_alfven
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_matrix_exponential, only: matrix_exponential
    use milca_drift_spectral, only: drift_grid_t, drift_work_t, drift_grid_init, &
        drift_work_init, drift_bracket, drift_poisson_inverse, drift_enforce_reality
    use milca_drift_fluid_leaves, only: drift_hyperviscosity, da_inductance, &
        da_linear_matrix, da_mode_ledger
    implicit none
    private
    public :: da_parameters_t, da_state_t, da_ledger_t
    public :: da_initialize, da_step, da_ledger, da_linear_operator, da_current
    public :: da_rhs, da_rhs_jvp, da_rhs_vjp

    !> Units: lengths rho_s, time L_n/c_s, fluctuations times L_n/rho_s.
    !> beta = (beta_e/2)(L_par/L_n)^2, mu = (m_e/m_i)(L_par/L_n)^2,
    !> c = 0.51 (nu_e L_n/c_s) mu, kpar in 1/L_par, kappa the drive (1 for the
    !> density gradient that sets L_n), damping d (-laplace)^order on every field.
    type :: da_parameters_t
        real(dp) :: beta = 0, mu = 0, c = 0, kpar = 1, kappa = 1, d = 0
        integer :: order = 2
    end type

    type :: da_state_t
        type(da_parameters_t) :: par
        type(drift_grid_t) :: grid
        type(drift_work_t) :: work
        real(dp) :: dt = 0, time = 0
        complex(dp), allocatable :: zeta(:, :), n(:, :), psi(:, :)
        !> exp[(L - nu) dt/2] of each mode, (row, column, kx, ky).
        complex(dp), allocatable :: half(:, :, :, :)
        !> 1/inductance of each mode, j = psi/inductance; zero outside the band.
        real(dp), allocatable :: admittance(:, :)
        complex(dp), allocatable :: phi(:, :), current(:, :)
        complex(dp), allocatable :: stage(:, :, :), rate(:, :, :)
    end type

    !> Box averages; d energy/dt = flux - d_c - d_e and
    !> d enstrophy/dt = flux - d_w. kinetic = <|grad phi|^2>/2 and its part
    !> zonal_kinetic in the modes ky = 0; magnetic = beta <|grad a|^2>/2.
    type :: da_ledger_t
        real(dp) :: energy = 0, enstrophy = 0, flux = 0
        real(dp) :: d_c = 0, d_e = 0, d_w = 0
        real(dp) :: kinetic = 0, zonal_kinetic = 0, magnetic = 0
    end type
contains
    subroutine da_initialize(state, par, nx, ny, lx, ly, dt)
        type(da_state_t), intent(out) :: state
        type(da_parameters_t), intent(in) :: par
        integer, intent(in) :: nx, ny
        real(dp), intent(in) :: lx, ly, dt
        complex(dp) :: l(3, 3)
        real(dp) :: nu, ind
        integer :: i, j, k, info

        if (par%beta < 0 .or. par%mu < 0 .or. .not. par%beta + par%mu > 0) &
            error stop 'da_initialize: need beta, mu >= 0 and beta + mu > 0'
        if (par%c < 0 .or. par%d < 0 .or. par%order < 0) &
            error stop 'da_initialize: need c, d, order >= 0'
        if (.not. dt > 0) error stop 'da_initialize: need dt > 0'
        state%par = par
        state%dt = dt
        call drift_grid_init(state%grid, nx, ny, lx, ly)
        call drift_work_init(state%work, state%grid)
        allocate (state%zeta(nx, ny), state%n(nx, ny), state%psi(nx, ny))
        allocate (state%phi(nx, ny), state%current(nx, ny))
        allocate (state%half(3, 3, nx, ny), state%admittance(nx, ny))
        allocate (state%stage(nx, ny, 3), state%rate(nx, ny, 3))
        state%zeta = 0; state%n = 0; state%psi = 0
        state%half = 0; state%admittance = 0
        associate (g => state%grid)
            do j = 1, ny
                do i = 1, nx
                    do k = 1, 3
                        state%half(k, k, i, j) = 1
                    end do
                    if (.not. g%keep(i, j) .or. .not. g%k2(i, j) > 0) cycle
                    call da_inductance(par%beta, par%mu, g%kx(i), g%ky(j), ind)
                    state%admittance(i, j) = 1/ind
                    call drift_hyperviscosity(par%d, real(par%order, dp), g%kx(i), &
                                              g%ky(j), nu)
                    l = da_linear_operator(par, g%kx(i), g%ky(j))
                    do k = 1, 3
                        l(k, k) = l(k, k) - nu
                    end do
                    call matrix_exponential(l*(dt/2), state%half(:, :, i, j), info)
                    if (info /= 0) error stop 'da_initialize: mode propagator failed'
                end do
            end do
        end associate
    end subroutine

    !> L of d/dt (zeta, n, psi) = L (zeta, n, psi) for one mode, without damping.
    pure function da_linear_operator(par, kx, ky) result(l)
        type(da_parameters_t), intent(in) :: par
        real(dp), intent(in) :: kx, ky
        complex(dp) :: l(3, 3)
        real(dp) :: l13, l21i, l31, l32, l33r, l33i

        call da_linear_matrix(par%beta, par%mu, par%c, par%kpar, par%kappa, kx, ky, &
                              l13, l21i, l31, l32, l33r, l33i)
        l = 0
        l(1, 3) = l13
        l(2, 1) = cmplx(0, l21i, dp)
        l(2, 3) = l13
        l(3, 1) = l31
        l(3, 2) = l32
        l(3, 3) = cmplx(l33r, l33i, dp)
    end function

    subroutine da_step(state)
        type(da_state_t), intent(inout) :: state
        real(dp) :: dt

        dt = state%dt
        call propagate_half(state)
        associate (u => state%stage, k => state%rate, z => state%zeta, n => state%n, &
                   p => state%psi)
            call bracket_terms(state, z, n, p, k)
            u(:, :, 1) = z + dt*k(:, :, 1)
            u(:, :, 2) = n + dt*k(:, :, 2)
            u(:, :, 3) = p + dt*k(:, :, 3)
            call bracket_terms(state, u(:, :, 1), u(:, :, 2), u(:, :, 3), k)
            u(:, :, 1) = 0.75_dp*z + 0.25_dp*(u(:, :, 1) + dt*k(:, :, 1))
            u(:, :, 2) = 0.75_dp*n + 0.25_dp*(u(:, :, 2) + dt*k(:, :, 2))
            u(:, :, 3) = 0.75_dp*p + 0.25_dp*(u(:, :, 3) + dt*k(:, :, 3))
            call bracket_terms(state, u(:, :, 1), u(:, :, 2), u(:, :, 3), k)
            z = (z + 2*(u(:, :, 1) + dt*k(:, :, 1)))/3
            n = (n + 2*(u(:, :, 2) + dt*k(:, :, 2)))/3
            p = (p + 2*(u(:, :, 3) + dt*k(:, :, 3)))/3
        end associate
        call propagate_half(state)
        call drift_enforce_reality(state%grid, state%zeta)
        call drift_enforce_reality(state%grid, state%n)
        call drift_enforce_reality(state%grid, state%psi)
        state%time = state%time + dt
    end subroutine

    subroutine propagate_half(state)
        type(da_state_t), intent(inout) :: state
        complex(dp) :: v(3)
        integer :: i, j

        !$omp parallel do private(i, j, v) schedule(static)
        do j = 1, state%grid%ny
            do i = 1, state%grid%nx
                v = [state%zeta(i, j), state%n(i, j), state%psi(i, j)]
                v = matmul(state%half(:, :, i, j), v)
                state%zeta(i, j) = v(1); state%n(i, j) = v(2); state%psi(i, j) = v(3)
            end do
        end do
        !$omp end parallel do
    end subroutine

    !> j = psi/inductance mode by mode.
    subroutine da_current(state, psi, j)
        type(da_state_t), intent(in) :: state
        complex(dp), intent(in) :: psi(:, :)
        complex(dp), intent(out) :: j(:, :)
        j = state%admittance*psi
    end subroutine

    !> (-{phi, zeta}, -{phi, n}, -mu {phi, j}) with laplace(phi) = zeta.
    subroutine bracket_terms(state, zeta, n, psi, rate)
        type(da_state_t), intent(inout) :: state
        complex(dp), intent(in) :: zeta(:, :), n(:, :), psi(:, :)
        complex(dp), intent(out) :: rate(:, :, :)

        call drift_poisson_inverse(state%grid, zeta, state%phi)
        call drift_bracket(state%grid, state%work, state%phi, zeta, rate(:, :, 1))
        call drift_bracket(state%grid, state%work, state%phi, n, rate(:, :, 2))
        rate(:, :, 1:2) = -rate(:, :, 1:2)
        if (state%par%mu > 0) then
            call da_current(state, psi, state%current)
            call drift_bracket(state%grid, state%work, state%phi, state%current, &
                               rate(:, :, 3))
            rate(:, :, 3) = -state%par%mu*rate(:, :, 3)
        else
            rate(:, :, 3) = 0
        end if
    end subroutine

    !> Full right-hand side d/dt (zeta, n, psi), rate(:, :, field).
    subroutine da_rhs(state, zeta, n, psi, rate)
        type(da_state_t), intent(inout) :: state
        complex(dp), intent(in) :: zeta(:, :), n(:, :), psi(:, :)
        complex(dp), intent(out) :: rate(:, :, :)

        call bracket_terms(state, zeta, n, psi, rate)
        call add_linear(state, zeta, n, psi, rate, adjoint=.false.)
    end subroutine

    !> Tangent of da_rhs at (zeta, n, psi) along (tzeta, tn, tpsi).
    subroutine da_rhs_jvp(state, zeta, n, psi, tzeta, tn, tpsi, rate)
        type(da_state_t), intent(inout) :: state
        complex(dp), intent(in) :: zeta(:, :), n(:, :), psi(:, :)
        complex(dp), intent(in) :: tzeta(:, :), tn(:, :), tpsi(:, :)
        complex(dp), intent(out) :: rate(:, :, :)
        complex(dp), allocatable :: phi(:, :), tphi(:, :), part(:, :), j(:, :), tj(:, :)

        allocate (phi, tphi, part, j, tj, mold=zeta)
        call drift_poisson_inverse(state%grid, zeta, phi)
        call drift_poisson_inverse(state%grid, tzeta, tphi)
        call drift_bracket(state%grid, state%work, tphi, zeta, rate(:, :, 1))
        call drift_bracket(state%grid, state%work, phi, tzeta, part)
        rate(:, :, 1) = -rate(:, :, 1) - part
        call drift_bracket(state%grid, state%work, tphi, n, rate(:, :, 2))
        call drift_bracket(state%grid, state%work, phi, tn, part)
        rate(:, :, 2) = -rate(:, :, 2) - part
        call da_current(state, psi, j)
        call da_current(state, tpsi, tj)
        call drift_bracket(state%grid, state%work, tphi, j, rate(:, :, 3))
        call drift_bracket(state%grid, state%work, phi, tj, part)
        rate(:, :, 3) = -state%par%mu*(rate(:, :, 3) + part)
        call add_linear(state, tzeta, tn, tpsi, rate, adjoint=.false.)
    end subroutine

    !> Adjoint of da_rhs_jvp under the box average: cotangent b(:, :, field)
    !> of the rates to the cotangent g(:, :, field) of the state.
    subroutine da_rhs_vjp(state, zeta, n, psi, b, g)
        type(da_state_t), intent(inout) :: state
        complex(dp), intent(in) :: zeta(:, :), n(:, :), psi(:, :), b(:, :, :)
        complex(dp), intent(out) :: g(:, :, :)
        complex(dp), allocatable :: phi(:, :), part(:, :), source(:, :), j(:, :)

        allocate (phi, part, source, j, mold=zeta)
        call drift_poisson_inverse(state%grid, zeta, phi)
        call da_current(state, psi, j)
        ! <b, -{phi_t, f}> = -<phi_t, {f, b}> and <b, -{phi, f_t}> = <f_t, {phi, b}>;
        ! phi = P zeta and j = psi/inductance are self-adjoint mode by mode.
        call drift_bracket(state%grid, state%work, zeta, b(:, :, 1), source)
        call drift_bracket(state%grid, state%work, n, b(:, :, 2), part)
        source = -source - part
        call drift_bracket(state%grid, state%work, j, b(:, :, 3), part)
        source = source - state%par%mu*part
        call drift_poisson_inverse(state%grid, source, g(:, :, 1))
        call drift_bracket(state%grid, state%work, phi, b(:, :, 1), part)
        g(:, :, 1) = g(:, :, 1) + part
        call drift_bracket(state%grid, state%work, phi, b(:, :, 2), g(:, :, 2))
        call drift_bracket(state%grid, state%work, phi, b(:, :, 3), part)
        call da_current(state, state%par%mu*part, g(:, :, 3))
        call add_linear(state, b(:, :, 1), b(:, :, 2), b(:, :, 3), g, adjoint=.true.)
    end subroutine

    !> rate += (L - nu) (zeta, n, psi), or its conjugate transpose.
    subroutine add_linear(state, zeta, n, psi, rate, adjoint)
        type(da_state_t), intent(in) :: state
        complex(dp), intent(in) :: zeta(:, :), n(:, :), psi(:, :)
        complex(dp), intent(inout) :: rate(:, :, :)
        logical, intent(in) :: adjoint
        complex(dp) :: l(3, 3), v(3)
        real(dp) :: nu
        integer :: i, j, k

        associate (g => state%grid, par => state%par)
            !$omp parallel do private(i, j, k, l, v, nu) schedule(static)
            do j = 1, g%ny
                do i = 1, g%nx
                    if (.not. g%keep(i, j) .or. .not. g%k2(i, j) > 0) cycle
                    call drift_hyperviscosity(par%d, real(par%order, dp), g%kx(i), &
                                              g%ky(j), nu)
                    l = da_linear_operator(par, g%kx(i), g%ky(j))
                    if (adjoint) l = conjg(transpose(l))
                    do k = 1, 3
                        l(k, k) = l(k, k) - nu
                    end do
                    v = matmul(l, [zeta(i, j), n(i, j), psi(i, j)])
                    rate(i, j, :) = rate(i, j, :) + v
                end do
            end do
            !$omp end parallel do
        end associate
    end subroutine

    subroutine da_ledger(state, ledger)
        type(da_state_t), intent(in) :: state
        type(da_ledger_t), intent(out) :: ledger
        real(dp) :: nu, part(6), total(6), kinetic, zonal, magnetic, ek
        integer :: i, j

        total = 0; kinetic = 0; zonal = 0; magnetic = 0
        associate (g => state%grid, par => state%par)
            !$omp parallel do private(i, j, nu, part, ek) schedule(static) &
            !$omp reduction(+:total, kinetic, zonal, magnetic)
            do j = 1, g%ny
                do i = 1, g%nx
                    if (.not. g%keep(i, j) .or. .not. g%k2(i, j) > 0) cycle
                    call drift_hyperviscosity(par%d, real(par%order, dp), g%kx(i), &
                                              g%ky(j), nu)
                    call da_mode_ledger(par%beta, par%mu, par%c, par%kappa, g%kx(i), &
                        g%ky(j), nu, real(state%zeta(i, j), dp), aimag(state%zeta(i, j)), &
                        real(state%n(i, j), dp), aimag(state%n(i, j)), &
                        real(state%psi(i, j), dp), aimag(state%psi(i, j)), part(1), &
                        part(2), part(3), part(4), part(5), part(6))
                    total = total + part
                    ek = abs(state%zeta(i, j))**2/(2*g%k2(i, j))
                    kinetic = kinetic + ek
                    if (.not. abs(g%ky(j)) > 0) zonal = zonal + ek
                    magnetic = magnetic + par%beta*abs(state%admittance(i, j) &
                               *state%psi(i, j))**2/(2*g%k2(i, j))
                end do
            end do
            !$omp end parallel do
        end associate
        ledger = da_ledger_t(total(1), total(2), total(3), total(4), total(5), total(6), &
                             kinetic, zonal, magnetic)
    end subroutine
end module milca_drift_alfven
