!> Q7 part 2: 1X1V electron Wigner-Poisson (hbar > 0) and Vlasov-Poisson
!> (hbar = 0) on one periodic Fourier(x) x Fourier(v) grid with neutralizing
!> background (units in kin6d_wigner_leaves). Strang splitting: exact free
!> streaming in (k, v) by the shared streaming_rotation leaf, exact nonlocal
!> Wigner kick in (x, eta) as the phase s_k(eta) Vhat_k, Poisson at the half
!> step. The forward tangent d/dhbar propagates through every substep.
!> The Wigner function is not clipped; negative mass is a diagnostic.
module kin6d_wigner_poisson
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_fft, only: fortnum_fft_plan_t, fft_c2c_plan_init, fft_c2c
    use kin6d_streaming_leaves, only: streaming_rotation
    use kin6d_wigner_leaves, only: wigner_shift_symbol, wigner_shift_symbol_jvp, &
        vlasov_shift_symbol, wigner_potential_mode, wigner_field_mode, &
        wigner_phase_coefficient, wigner_phase_coefficient_jvp, wigner_kick, &
        wigner_kick_tangent, wigner_field_energy_mode, wigner_maxwellian
    implicit none
    private
    public :: wigner_grid_t, wigner_state_t, wigner_work_t, wigner_diagnostics_t
    public :: wigner_grid_init, wigner_work_init, wigner_initial_beams, wigner_step
    public :: wigner_diagnostics, wigner_fit_mode

    real(dp), parameter :: pi = acos(-1.0_dp)

    type :: wigner_grid_t
        integer :: nx = 0, nv = 0
        real(dp) :: length = 0, vmax = 0, dx = 0, dv = 0
        real(dp), allocatable :: x(:), v(:), k(:), eta(:)
        type(fortnum_fft_plan_t) :: px, pv
    end type

    type :: wigner_state_t
        real(dp), allocatable :: f(:, :), tangent(:, :)
        real(dp) :: time = 0, hbar = 0
        logical :: with_tangent = .false.
    end type

    !> Preallocated step workspace; the step itself allocates nothing.
    type :: wigner_work_t
        real(dp) :: hbar = -1
        real(dp), allocatable :: s(:, :), ds(:, :), phase(:, :), dphase(:, :)
        real(dp), allocatable :: density(:), ddensity(:)
        complex(dp), allocatable :: vhat(:), dvhat(:), bx(:), tx(:), bv(:), tv(:)
    end type

    type :: wigner_diagnostics_t
        real(dp) :: mass = 0, momentum = 0, kinetic = 0, field = 0, energy = 0
        real(dp) :: negative_mass = 0, edge_mass = 0
        !> Complex amplitude of cos(k1 x) in the density and its hbar tangent.
        complex(dp) :: mode = 0, dmode = 0
    end type

contains

    subroutine wigner_grid_init(grid, nx, nv, length, vmax)
        type(wigner_grid_t), intent(out) :: grid
        integer, intent(in) :: nx, nv
        real(dp), intent(in) :: length, vmax
        integer :: j
        grid%nx = nx; grid%nv = nv; grid%length = length; grid%vmax = vmax
        grid%dx = length/nx; grid%dv = 2*vmax/nv
        allocate (grid%x(nx), grid%v(nv), grid%k(nx), grid%eta(nv))
        do j = 1, nx
            grid%x(j) = (j - 1)*grid%dx
            grid%k(j) = 2*pi*signed(j, nx)/length
        end do
        do j = 1, nv
            grid%v(j) = -vmax + (j - 1)*grid%dv
            grid%eta(j) = 2*pi*signed(j, nv)/(2*vmax)
        end do
        call fft_c2c_plan_init(grid%px, nx)
        call fft_c2c_plan_init(grid%pv, nv)
    end subroutine

    !> Signed DFT index of position j; the Nyquist index returns 0 (filtered).
    pure integer function signed(j, n) result(m)
        integer, intent(in) :: j, n
        m = j - 1
        if (m > n/2) m = m - n
        if (2*m == n) m = 0
    end function

    pure logical function nyquist(j, n)
        integer, intent(in) :: j, n
        nyquist = 2*(j - 1) == n
    end function

    !> Workspace and the shift-symbol table s_k(eta) (and d/dhbar) for hbar.
    subroutine wigner_work_init(grid, hbar, work)
        type(wigner_grid_t), intent(in) :: grid
        real(dp), intent(in) :: hbar
        type(wigner_work_t), intent(out) :: work
        integer :: m, l
        allocate (work%s(grid%nx, grid%nv), work%ds(grid%nx, grid%nv))
        allocate (work%phase(grid%nx, grid%nv), work%dphase(grid%nx, grid%nv))
        allocate (work%density(grid%nx), work%ddensity(grid%nx))
        allocate (work%vhat(grid%nx), work%dvhat(grid%nx), work%bx(grid%nx), &
                  work%tx(grid%nx), work%bv(grid%nv), work%tv(grid%nv))
        work%hbar = hbar
        do l = 1, grid%nv
            do m = 1, grid%nx
                if (hbar > 0) then
                    call wigner_shift_symbol(hbar, grid%k(m), grid%eta(l), work%s(m, l))
                    call wigner_shift_symbol_jvp(hbar, grid%k(m), grid%eta(l), &
                        1.0_dp, 0.0_dp, 0.0_dp, work%ds(m, l))
                else
                    call vlasov_shift_symbol(grid%k(m), grid%eta(l), work%s(m, l))
                    work%ds(m, l) = 0
                end if
            end do
        end do
    end subroutine

    !> f = sum_b Maxwellian_b(v) (1 + alpha cos(2 pi mode x/L)); zero tangent.
    subroutine wigner_initial_beams(grid, weight, drift, sigma, alpha, mode, hbar, &
                                    state, with_tangent)
        type(wigner_grid_t), intent(in) :: grid
        real(dp), intent(in) :: weight(:), drift(:), sigma(:), alpha, hbar
        integer, intent(in) :: mode
        type(wigner_state_t), intent(out) :: state
        logical, intent(in) :: with_tangent
        real(dp) :: value, f0
        integer :: j, l, b
        allocate (state%f(grid%nx, grid%nv))
        do l = 1, grid%nv
            f0 = 0
            do b = 1, size(weight)
                call wigner_maxwellian(grid%v(l), drift(b), sigma(b), weight(b), value)
                f0 = f0 + value
            end do
            do j = 1, grid%nx
                state%f(j, l) = f0*(1 + alpha*cos(2*pi*mode*grid%x(j)/grid%length))
            end do
        end do
        state%hbar = hbar
        state%with_tangent = with_tangent
        if (with_tangent) then
            allocate (state%tangent(grid%nx, grid%nv))
            state%tangent = 0
        end if
    end subroutine

    !> One Strang step: stream dt/2, Poisson, Wigner/Vlasov kick dt, stream dt/2.
    subroutine wigner_step(grid, work, state, dt)
        type(wigner_grid_t), intent(in) :: grid
        type(wigner_work_t), intent(inout) :: work
        type(wigner_state_t), intent(inout) :: state
        real(dp), intent(in) :: dt
        if (abs(work%hbar - state%hbar) > 0) error stop 'wigner_step: work hbar mismatch'
        call stream(grid, work, state, dt/2)
        call poisson(grid, work, state)
        call kick(grid, work, state, dt)
        call stream(grid, work, state, dt/2)
        state%time = state%time + dt
    end subroutine

    subroutine stream(grid, work, state, tau)
        type(wigner_grid_t), intent(in) :: grid
        type(wigner_work_t), intent(inout) :: work
        type(wigner_state_t), intent(inout) :: state
        real(dp), intent(in) :: tau
        integer :: l
        do l = 1, grid%nv
            call stream_column(state%f(:, l), work%bx, grid%v(l))
            if (state%with_tangent) call stream_column(state%tangent(:, l), work%tx, &
                                                       grid%v(l))
        end do
    contains
        subroutine stream_column(column, buffer, v)
            real(dp), intent(inout) :: column(:)
            complex(dp), intent(inout) :: buffer(:)
            real(dp), intent(in) :: v
            real(dp) :: a, b
            integer :: m
            buffer = cmplx(column, 0.0_dp, dp)
            call fft_c2c(buffer, -1, grid%px)
            do m = 1, grid%nx
                if (nyquist(m, grid%nx)) then
                    buffer(m) = 0
                    cycle
                end if
                ! Coefficients (a, b) of cos/sin rotate exactly as in the
                ! classical free-streaming case: f(x - v tau).
                call streaming_rotation(tau, v, grid%k(m), real(buffer(m), dp), &
                                        -aimag(buffer(m)), a, b)
                buffer(m) = cmplx(a, -b, dp)
            end do
            call fft_c2c(buffer, 1, grid%px)
            column = real(buffer, dp)/grid%nx
        end subroutine
    end subroutine

    !> Vhat_k = nhat_k/k^2 (electron potential energy V = -phi) and its tangent.
    subroutine poisson(grid, work, state)
        type(wigner_grid_t), intent(in) :: grid
        type(wigner_work_t), intent(inout) :: work
        type(wigner_state_t), intent(in) :: state
        integer :: j
        do j = 1, grid%nx
            work%density(j) = grid%dv*sum(state%f(j, :))
        end do
        call potential(work%density, work%vhat)
        if (state%with_tangent) then
            do j = 1, grid%nx
                work%ddensity(j) = grid%dv*sum(state%tangent(j, :))
            end do
            call potential(work%ddensity, work%dvhat)
        end if
    contains
        subroutine potential(n, vhat)
            real(dp), intent(in) :: n(:)
            complex(dp), intent(out) :: vhat(:)
            real(dp) :: vr, vi
            integer :: m
            work%bx = cmplx(n, 0.0_dp, dp)
            call fft_c2c(work%bx, -1, grid%px)
            do m = 1, grid%nx
                if (m == 1 .or. nyquist(m, grid%nx)) then
                    vhat(m) = 0
                else
                    call wigner_potential_mode(grid%k(m), real(work%bx(m), dp), &
                                               aimag(work%bx(m)), vr, vi)
                    vhat(m) = cmplx(vr, vi, dp)
                end if
            end do
        end subroutine
    end subroutine

    subroutine kick(grid, work, state, dt)
        type(wigner_grid_t), intent(in) :: grid
        type(wigner_work_t), intent(inout) :: work
        type(wigner_state_t), intent(inout) :: state
        real(dp), intent(in) :: dt
        real(dp) :: pr, pim, dpr, dpi, gr, gi, n
        integer :: j, l, m
        n = real(grid%nx, dp)
        ! Phase(x_j, eta_l) = sum_k (dt/N) i s_k(eta_l) Vhat_k exp(i k x_j).
        do l = 1, grid%nv
            do m = 1, grid%nx
                call wigner_phase_coefficient(dt, work%s(m, l), real(work%vhat(m), dp), &
                    aimag(work%vhat(m)), n, pr, pim)
                work%bx(m) = cmplx(pr, pim, dp)
                if (state%with_tangent) then
                    call wigner_phase_coefficient_jvp(dt, work%s(m, l), &
                        real(work%vhat(m), dp), aimag(work%vhat(m)), n, 0.0_dp, &
                        work%ds(m, l), real(work%dvhat(m), dp), aimag(work%dvhat(m)), &
                        0.0_dp, dpr, dpi)
                    work%tx(m) = cmplx(dpr, dpi, dp)
                end if
            end do
            call fft_c2c(work%bx, 1, grid%px)
            work%phase(:, l) = real(work%bx, dp)
            if (state%with_tangent) then
                call fft_c2c(work%tx, 1, grid%px)
                work%dphase(:, l) = real(work%tx, dp)
            end if
        end do
        do j = 1, grid%nx
            work%bv = cmplx(state%f(j, :), 0.0_dp, dp)
            call fft_c2c(work%bv, -1, grid%pv)
            if (state%with_tangent) then
                work%tv = cmplx(state%tangent(j, :), 0.0_dp, dp)
                call fft_c2c(work%tv, -1, grid%pv)
            end if
            do l = 1, grid%nv
                if (nyquist(l, grid%nv)) then
                    work%bv(l) = 0; work%tv(l) = 0
                    cycle
                end if
                if (state%with_tangent) then
                    call wigner_kick_tangent(real(work%bv(l), dp), aimag(work%bv(l)), &
                        real(work%tv(l), dp), aimag(work%tv(l)), work%phase(j, l), &
                        work%dphase(j, l), gr, gi)
                    work%tv(l) = cmplx(gr, gi, dp)
                end if
                call wigner_kick(real(work%bv(l), dp), aimag(work%bv(l)), &
                                 work%phase(j, l), gr, gi)
                work%bv(l) = cmplx(gr, gi, dp)
            end do
            call fft_c2c(work%bv, 1, grid%pv)
            state%f(j, :) = real(work%bv, dp)/grid%nv
            if (state%with_tangent) then
                call fft_c2c(work%tv, 1, grid%pv)
                state%tangent(j, :) = real(work%tv, dp)/grid%nv
            end if
        end do
    end subroutine

    !> Mass, momentum, kinetic and field energy, negative mass (not clipped),
    !> edge mass |v| > 0.9 vmax (velocity-truncation indicator), and the
    !> density amplitude of the first harmonic with its hbar tangent.
    subroutine wigner_diagnostics(grid, work, state, diag)
        type(wigner_grid_t), intent(in) :: grid
        type(wigner_work_t), intent(inout) :: work
        type(wigner_state_t), intent(in) :: state
        type(wigner_diagnostics_t), intent(out) :: diag
        real(dp) :: cell, er, ei, total
        integer :: j, l, m
        cell = grid%dx*grid%dv
        do l = 1, grid%nv
            diag%mass = diag%mass + cell*sum(state%f(:, l))
            diag%momentum = diag%momentum + cell*grid%v(l)*sum(state%f(:, l))
            diag%kinetic = diag%kinetic + cell*grid%v(l)**2/2*sum(state%f(:, l))
            diag%negative_mass = diag%negative_mass + cell*sum(max(-state%f(:, l), 0.0_dp))
            if (abs(grid%v(l)) > 0.9_dp*grid%vmax) &
                diag%edge_mass = diag%edge_mass + cell*sum(abs(state%f(:, l)))
        end do
        do j = 1, grid%nx
            work%density(j) = grid%dv*sum(state%f(j, :))
        end do
        work%bx = cmplx(work%density, 0.0_dp, dp)
        call fft_c2c(work%bx, -1, grid%px)
        do m = 2, grid%nx
            if (nyquist(m, grid%nx)) cycle
            call wigner_field_mode(grid%k(m), real(work%bx(m), dp), aimag(work%bx(m)), &
                                   er, ei)
            total = diag%field
            call wigner_field_energy_mode(er, ei, grid%length, real(grid%nx, dp), &
                                          total, diag%field)
        end do
        diag%energy = diag%kinetic + diag%field
        diag%mode = 2*work%bx(2)/grid%nx
        if (state%with_tangent) then
            do j = 1, grid%nx
                work%ddensity(j) = grid%dv*sum(state%tangent(j, :))
            end do
            work%tx = cmplx(work%ddensity, 0.0_dp, dp)
            call fft_c2c(work%tx, -1, grid%px)
            diag%dmode = 2*work%tx(2)/grid%nx
        end if
    end subroutine

    !> Dominant complex frequency of a real mode signal sampled at spacing tau
    !> from the least-squares two-term linear prediction s(n+1) = p s(n) - q s(n-1),
    !> exact for A exp(-i omega t) + c.c. or for two real exponentials. Returns
    !> omega with Im(omega) the growth rate (negative: damping).
    subroutine wigner_fit_mode(signal, tau, omega, ok)
        real(dp), intent(in) :: signal(:), tau
        complex(dp), intent(out) :: omega
        logical, intent(out) :: ok
        real(dp) :: a11, a12, a22, r1, r2, det, p, q, disc
        complex(dp) :: z1, z2
        integer :: n
        a11 = 0; a12 = 0; a22 = 0; r1 = 0; r2 = 0
        do n = 2, size(signal) - 1
            a11 = a11 + signal(n)**2
            a12 = a12 - signal(n)*signal(n - 1)
            a22 = a22 + signal(n - 1)**2
            r1 = r1 + signal(n + 1)*signal(n)
            r2 = r2 - signal(n + 1)*signal(n - 1)
        end do
        det = a11*a22 - a12**2
        ok = abs(det) > 0
        if (.not. ok) return
        p = (r1*a22 - a12*r2)/det
        q = (a11*r2 - a12*r1)/det
        disc = p**2 - 4*q
        z1 = (p + sqrt(cmplx(disc, 0.0_dp, dp)))/2
        z2 = (p - sqrt(cmplx(disc, 0.0_dp, dp)))/2
        if (disc < 0) then
            if (aimag(z2) < aimag(z1)) z1 = z2
        else
            if (abs(z2) > abs(z1)) z1 = z2
        end if
        ! z = exp(-i omega tau): omega = i log(z)/tau.
        omega = (0.0_dp, 1.0_dp)*log(z1)/tau
        if (real(omega, dp) < 0) omega = cmplx(-real(omega, dp), aimag(omega), dp)
    end subroutine
end module kin6d_wigner_poisson
