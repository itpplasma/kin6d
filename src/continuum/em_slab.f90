!> Plane electromagnetic wave at normal incidence on a cold electron plasma
!> with a sharp edge (docs/em-slab-reflection.md).
!>
!> Model: relativistic Vlasov--Maxwell in one spatial dimension x with
!> velocity components (v_x, v_y); fields E_x, E_y, B_z; fixed neutralizing
!> ions. Carrier: macroparticle sheets (PIC) with linear (CIC) shape on a
!> staggered Yee grid. Units: t omega_pe, x omega_pe/c, E e/(m c omega_pe),
!> c B e/(m c omega_pe), u = p/(m c), J/(e n0 c), energy n0 m c^2 c/omega_pe
!> per unit area. Electrons carry q/m = -1; the plasma density is n0 = 1.
!>
!> Time levels: x, E at integer steps n; u, B_z, J at n + 1/2. One step:
!> half Faraday to B^n, gather E^n and B^n at x^n, Boris push, move, deposit
!> J_y at the midpoint and J_x with the charge-conserving (continuity-exact)
!> CIC rule, half Faraday to B^(n+1/2), Ampere to E^(n+1), first-order Mur
!> outgoing condition for E_y at both ends. A soft sheet current at the
!> source node launches the wave. The discrete Gauss law is preserved by the
!> J_x deposit from the neutral initial state.
!>
!> Status: candidate. Reported reflection/transmission/skin-depth values are
!> lock-in estimates with indicator errors from independent refinement, not
!> certificates.
module kin6d_em_slab
    use, intrinsic :: iso_fortran_env, only: dp => real64
    !$ use omp_lib, only: omp_get_max_threads, omp_get_thread_num
    use kin6d_em_slab_leaves, only: em_slab_boris, em_slab_velocity, em_slab_cic, &
        em_slab_faraday, em_slab_ampere, em_slab_mur, em_slab_ramp, &
        em_slab_sheet_current, em_slab_field_energy, em_slab_kinetic_energy, &
        em_slab_characteristics
    implicit none
    private

    !> Geometry, wave and discretization of one run.
    type, public :: em_slab_config_t
        real(dp) :: omega = 2.0_dp         !< wave frequency [omega_pe]
        real(dp) :: amplitude = 5.0e-3_dp  !< incident E_y amplitude [m c omega_pe/e]
        real(dp) :: ramp = 25.0_dp         !< sin^2 envelope rise time [1/omega_pe]
        real(dp) :: dx = 0.04_dp           !< cell width [c/omega_pe]
        real(dp) :: courant = 0.95_dp      !< c dt/dx
        real(dp) :: x_left = -20.0_dp      !< left boundary [c/omega_pe]
        real(dp) :: x_source = -18.0_dp    !< source sheet
        real(dp) :: x_edge = 0.0_dp        !< plasma edge
        real(dp) :: x_right = 135.0_dp     !< right boundary
        integer :: per_cell = 8            !< electron sheets per cell
        logical :: plasma = .true.         !< .false. gives the vacuum run
    end type

    !> Field and particle state with preallocated per-thread scratch.
    type, public :: em_slab_t
        type(em_slab_config_t) :: config
        integer :: cells = 0, source_node = 0, steps = 0, threads = 1
        real(dp) :: dt = 0, time = 0, weight = 0
        !> Last completed step: kinetic energy at the old integer time level.
        real(dp) :: kinetic = 0
        real(dp), allocatable :: ey(:), jy(:), rho_ion(:)        ! nodes 0..cells
        real(dp), allocatable :: bz(:), ex(:), jx(:), bn(:)      ! i -> (i + 1/2) dx
        real(dp), allocatable :: x(:), y(:), ux(:), uy(:), x0(:)
        real(dp), allocatable :: jy_thread(:, :), jx_thread(:, :), k_thread(:)
    end type

    !> Steady-state lock-in and energy ledger of one run.
    type, public :: em_slab_result_t
        real(dp) :: window_start = 0, window_end = 0
        real(dp), allocatable :: x(:)                 !< node positions
        complex(dp), allocatable :: amplitude(:)      !< E_y lock-in per node
        complex(dp) :: forward = 0, backward = 0      !< characteristics at probe
        complex(dp) :: forward_far = 0                !< forward at second probe
        real(dp) :: probe = 0, probe_far = 0
        !> Time integrals [n0 m c^3/omega_pe]: incident (forward) and reflected
        !> (backward) Poynting flux at the probe, outflow at the right end;
        !> changes of field energy right of the probe and of particle energy.
        real(dp) :: incident = 0, reflected = 0, outflow = 0
        real(dp) :: field_change = 0, kinetic_change = 0, balance = 0
        real(dp) :: gauss_residual = 0, max_momentum = 0
    end type

    !> Optional frame recording for visualization.
    type, public :: em_slab_frames_t
        real(dp) :: interval = 0.25_dp, x_min = -15.0_dp, x_max = 10.0_dp
        integer :: stride = 2
        real(dp) :: particle_min = 0.0_dp, particle_max = 4.0_dp
        integer :: particle_stride = 16
        real(dp), allocatable :: t(:), x(:), ey(:, :), x0(:), px(:, :), py(:, :)
    end type

    public :: em_slab_init, em_slab_step, em_slab_magnetic_now, em_slab_run
    public :: em_slab_field_energy_range, em_slab_gauss_residual
    public :: em_slab_characteristics_at
    public :: em_slab_reflection, em_slab_decay_length, em_slab_mean_amplitude
    public :: em_slab_wavenumber

contains

    !> Allocate fields and load the quiet cold electron lattice. The ion
    !> background equals the initial electron charge on the grid, so the
    !> initial state is exactly neutral and E = B = 0.
    subroutine em_slab_init(s, config)
        type(em_slab_t), intent(out) :: s
        type(em_slab_config_t), intent(in) :: config
        integer :: plasma_cells, first, k, n
        real(dp) :: spacing
        if (.not. (config%dx > 0 .and. config%courant > 0 .and. config%courant < 1)) &
            error stop 'em_slab: need dx > 0 and 0 < courant < 1'
        if (.not. (config%x_left < config%x_source)) error stop 'em_slab: source'
        if (.not. (config%x_source < config%x_edge)) error stop 'em_slab: edge'
        if (.not. (config%x_edge < config%x_right)) error stop 'em_slab: right'
        if (config%per_cell < 1) error stop 'em_slab: per_cell >= 1'
        s%config = config
        s%cells = nint((config%x_right - config%x_left)/config%dx)
        s%dt = config%courant*config%dx
        s%source_node = nint((config%x_source - config%x_left)/config%dx)
        s%threads = 1
        !$ s%threads = omp_get_max_threads()
        allocate (s%ey(0:s%cells), s%jy(0:s%cells), s%rho_ion(0:s%cells))
        allocate (s%bz(0:s%cells - 1), s%ex(0:s%cells - 1), s%jx(0:s%cells - 1))
        allocate (s%bn(0:s%cells - 1))
        allocate (s%jy_thread(0:s%cells, s%threads))
        allocate (s%jx_thread(0:s%cells - 1, s%threads))
        allocate (s%k_thread(s%threads))
        s%ey = 0; s%jy = 0; s%bz = 0; s%ex = 0; s%jx = 0; s%bn = 0; s%rho_ion = 0
        ! Plasma from the edge node to two cells before the right boundary.
        first = nint((config%x_edge - config%x_left)/config%dx)
        plasma_cells = s%cells - 2 - first
        if (.not. config%plasma) plasma_cells = 0
        n = plasma_cells*config%per_cell
        allocate (s%x(n), s%y(n), s%ux(n), s%uy(n), s%x0(n))
        spacing = config%dx/config%per_cell
        s%weight = spacing                       ! n0 = 1 sheets per unit length
        do k = 1, n
            s%x(k) = config%x_left + first*config%dx + (k - 0.5_dp)*spacing
        end do
        s%x0 = s%x; s%y = 0; s%ux = 0; s%uy = 0
        call electron_density(s, s%rho_ion)
        s%rho_ion = -s%rho_ion
        s%time = 0; s%steps = 0; s%kinetic = 0
    end subroutine

    !> Advance one step (x^n, u^(n-1/2), E^n, B^(n-1/2)) -> level n + 1.
    subroutine em_slab_step(s)
        type(em_slab_t), intent(inout) :: s
        real(dp) :: r_half, r_full, envelope, sheet, edge_left, edge_right
        integer :: i
        r_half = s%dt/(2*s%config%dx); r_full = s%dt/s%config%dx
        call em_slab_magnetic_now(s, s%bn)
        call push_and_deposit(s)
        ! Second half of Faraday to B^(n+1/2).
        call em_slab_faraday(s%bn, s%ey(0:s%cells - 1), s%ey(1:s%cells), r_half, s%bz)
        ! Source sheet current at t^(n+1/2).
        envelope = 1
        if (s%time + s%dt/2 < s%config%ramp) &
            call em_slab_ramp(s%time + s%dt/2, s%config%ramp, envelope)
        call em_slab_sheet_current(s%time + s%dt/2, s%config%omega, &
                                   s%config%amplitude, envelope, sheet)
        s%jy(s%source_node) = s%jy(s%source_node) + sheet/s%config%dx
        edge_left = s%ey(1); edge_right = s%ey(s%cells - 1)
        do i = 1, s%cells - 1
            call em_slab_ampere(s%ey(i), s%bz(i - 1), s%bz(i), s%jy(i), r_full, &
                                s%dt, s%ey(i))
        end do
        call em_slab_ampere(s%ex, 0.0_dp, 0.0_dp, s%jx, r_full, s%dt, s%ex)
        call em_slab_mur(s%ey(0), edge_left, s%ey(1), s%config%courant, s%ey(0))
        call em_slab_mur(s%ey(s%cells), edge_right, s%ey(s%cells - 1), &
                         s%config%courant, s%ey(s%cells))
        s%steps = s%steps + 1
        s%time = s%steps*s%dt
    end subroutine

    !> B_z at the current integer level: B^(n-1/2) advanced by dt/2 with E^n.
    subroutine em_slab_magnetic_now(s, bn)
        type(em_slab_t), intent(in) :: s
        real(dp), intent(out) :: bn(0:)
        call em_slab_faraday(s%bz, s%ey(0:s%cells - 1), s%ey(1:s%cells), &
                             s%dt/(2*s%config%dx), bn(0:s%cells - 1))
    end subroutine

    !> Gather, Boris, move and deposit. Each thread owns a contiguous static
    !> particle block and a private current column; columns are summed in
    !> thread order, so a fixed thread count is bitwise reproducible.
    subroutine push_and_deposit(s)
        type(em_slab_t), intent(inout) :: s
        real(dp) :: xi, f, w0, w1, eyp, exp_, bzp, uxn, uyn, vx, vy, g_old, g_new
        real(dp) :: xn, shift, inv_dx, inv_dt, old(0:2), new(0:2), running
        integer :: k, j, j0, tid, m, low, high
        inv_dx = 1/s%config%dx; inv_dt = 1/s%dt
        low = 1; high = s%cells - 2
        s%jy_thread = 0; s%jx_thread = 0; s%k_thread = 0
        !$omp parallel do schedule(static) default(shared) private(k, j, j0, tid, &
        !$omp& m, xi, f, w0, w1, eyp, exp_, bzp, uxn, uyn, vx, vy, g_old, g_new, xn, &
        !$omp& shift, old, new, running)
        do k = 1, size(s%x)
            tid = 1
            !$ tid = omp_get_thread_num() + 1
            ! Gather E_y at nodes, E_x and B_z at half nodes.
            xi = (s%x(k) - s%config%x_left)*inv_dx
            j = floor(xi); f = xi - j
            if (j < low .or. j > high) error stop 'em_slab: electron left the domain'
            call em_slab_cic(f, w0, w1)
            eyp = w0*s%ey(j) + w1*s%ey(j + 1)
            shift = xi - 0.5_dp
            j = floor(shift); f = shift - j
            call em_slab_cic(f, w0, w1)
            exp_ = w0*s%ex(j) + w1*s%ex(j + 1)
            bzp = w0*s%bn(j) + w1*s%bn(j + 1)
            call em_slab_kinetic_energy(s%ux(k), s%uy(k), g_old)
            call em_slab_boris(s%ux(k), s%uy(k), exp_, eyp, bzp, s%dt/2, uxn, uyn)
            call em_slab_kinetic_energy(uxn, uyn, g_new)
            s%k_thread(tid) = s%k_thread(tid) + s%weight*(g_old + g_new)/2
            call em_slab_velocity(uxn, uyn, vx, vy, g_new)
            s%ux(k) = uxn; s%uy(k) = uyn
            xn = s%x(k) + s%dt*vx
            s%y(k) = s%y(k) + s%dt*vy
            ! J_y at the midpoint position (charge -1, weight w).
            xi = ((s%x(k) + xn)/2 - s%config%x_left)*inv_dx
            j = floor(xi); f = xi - j
            call em_slab_cic(f, w0, w1)
            s%jy_thread(j, tid) = s%jy_thread(j, tid) - s%weight*vy*w0*inv_dx
            s%jy_thread(j + 1, tid) = s%jy_thread(j + 1, tid) - s%weight*vy*w1*inv_dx
            ! Charge-conserving J_x: J_(i+1/2) - J_(i-1/2) = w (S_i^new - S_i^old)/dt.
            j0 = floor((min(s%x(k), xn) - s%config%x_left)*inv_dx)
            call shape3(s%x(k), j0, old)
            call shape3(xn, j0, new)
            running = 0
            do m = 0, 1
                running = running + s%weight*(new(m) - old(m))*inv_dt
                s%jx_thread(j0 + m, tid) = s%jx_thread(j0 + m, tid) + running
            end do
            s%x(k) = xn
        end do
        !$omp end parallel do
        !$omp parallel do schedule(static) private(m)
        do k = 0, s%cells
            s%jy(k) = 0
            do m = 1, s%threads
                s%jy(k) = s%jy(k) + s%jy_thread(k, m)
            end do
        end do
        !$omp end parallel do
        !$omp parallel do schedule(static) private(m)
        do k = 0, s%cells - 1
            s%jx(k) = 0
            do m = 1, s%threads
                s%jx(k) = s%jx(k) + s%jx_thread(k, m)
            end do
        end do
        !$omp end parallel do
        s%kinetic = 0
        do m = 1, s%threads
            s%kinetic = s%kinetic + s%k_thread(m)
        end do
    contains
        !> CIC node weights of a particle at x on nodes j0, j0 + 1, j0 + 2.
        pure subroutine shape3(xp, base, weights)
            real(dp), intent(in) :: xp
            integer, intent(in) :: base
            real(dp), intent(out) :: weights(0:2)
            real(dp) :: q, a, b
            integer :: jj
            q = (xp - s%config%x_left)*inv_dx
            jj = floor(q)
            call em_slab_cic(q - jj, a, b)
            weights = 0
            weights(jj - base) = a
            weights(jj - base + 1) = b
        end subroutine
    end subroutine

    !> Electron charge density -sum w S_i/dx on the nodes.
    subroutine electron_density(s, rho)
        type(em_slab_t), intent(in) :: s
        real(dp), intent(out) :: rho(0:)
        real(dp) :: xi, w0, w1
        integer :: k, j
        rho(0:s%cells) = 0
        do k = 1, size(s%x)
            xi = (s%x(k) - s%config%x_left)/s%config%dx
            j = floor(xi)
            call em_slab_cic(xi - j, w0, w1)
            rho(j) = rho(j) - s%weight*w0/s%config%dx
            rho(j + 1) = rho(j + 1) - s%weight*w1/s%config%dx
        end do
    end subroutine

    !> max_i |(E_x(i+1/2) - E_x(i-1/2))/dx - rho_i| over interior nodes.
    real(dp) function em_slab_gauss_residual(s) result(residual)
        type(em_slab_t), intent(in) :: s
        real(dp), allocatable :: rho(:)
        integer :: i
        allocate (rho(0:s%cells))
        call electron_density(s, rho)
        rho = rho + s%rho_ion
        residual = 0
        do i = 1, s%cells - 1
            residual = max(residual, abs((s%ex(i) - s%ex(i - 1))/s%config%dx - rho(i)))
        end do
    end function

    !> Field energy per area on nodes first..last (trapezoid for E_y, half
    !> nodes in between for E_x, B_z), with B_z at the integer level bn.
    real(dp) function em_slab_field_energy_range(s, bn, first, last) result(total)
        type(em_slab_t), intent(in) :: s
        real(dp), intent(in) :: bn(0:)
        integer, intent(in) :: first, last
        real(dp) :: w
        integer :: i
        total = 0
        do i = first, last
            call em_slab_field_energy(0.0_dp, s%ey(i), 0.0_dp, w)
            if (i == first .or. i == last) w = w/2
            total = total + w
        end do
        do i = first, last - 1
            call em_slab_field_energy(s%ex(i), 0.0_dp, bn(i), w)
            total = total + w
        end do
        total = total*s%config%dx
    end function

    !> Vacuum characteristics (E_y +- c B_z)/2 at node i, B_z averaged to the node.
    subroutine em_slab_characteristics_at(s, bn, i, forward, backward)
        type(em_slab_t), intent(in) :: s
        real(dp), intent(in) :: bn(0:)
        integer, intent(in) :: i
        real(dp), intent(out) :: forward, backward
        call em_slab_characteristics(s%ey(i), (bn(i - 1) + bn(i))/2, forward, backward)
    end subroutine

    !> Run to window_start + periods*2 pi/omega; lock in E_y on all nodes and
    !> the characteristics at two vacuum probes over the final window of
    !> `periods` whole periods; accumulate the energy ledger right of probe.
    subroutine em_slab_run(config, window_start, periods, probe, probe_far, result, &
                           frames)
        type(em_slab_config_t), intent(in) :: config
        real(dp), intent(in) :: window_start, probe, probe_far
        integer, intent(in) :: periods
        type(em_slab_result_t), intent(out) :: result
        type(em_slab_frames_t), intent(inout), optional :: frames
        type(em_slab_t) :: s
        real(dp), allocatable :: bn(:)
        real(dp) :: period, fwd, bwd, fwd2, bwd2, s_out, f_in, b_in, weight
        real(dp) :: field_start, energy_now, kinetic_start, time_now
        real(dp) :: f_prev, b_prev, out_prev
        complex(dp) :: rotor
        integer :: first_step, window_steps, total_steps, n, ip, ip2, i
        integer :: frame_count, next_frame
        call em_slab_init(s, config)
        allocate (bn(0:s%cells - 1))
        period = 8*atan(1.0_dp)/config%omega
        first_step = nint(window_start/s%dt)
        window_steps = nint(periods*period/s%dt)
        total_steps = first_step + window_steps
        if (window_steps < 8) error stop 'em_slab_run: window too short'
        if (total_steps*s%dt > config%x_right - config%x_source - 2*config%dx) &
            error stop 'em_slab_run: signal would reach the right boundary'
        ip = nint((probe - config%x_left)/config%dx)
        ip2 = nint((probe_far - config%x_left)/config%dx)
        if (ip <= s%source_node .or. ip2 <= ip .or. &
            ip2 >= nint((config%x_edge - config%x_left)/config%dx)) &
            error stop 'em_slab_run: need source < probe < probe_far < edge'
        allocate (result%x(0:s%cells), result%amplitude(0:s%cells))
        result%x = [(config%x_left + i*config%dx, i=0, s%cells)]
        result%amplitude = 0
        result%probe = result%x(ip); result%probe_far = result%x(ip2)
        result%window_start = first_step*s%dt
        result%window_end = total_steps*s%dt
        if (present(frames)) then
            frame_count = int(total_steps*s%dt/frames%interval) + 1
            call frame_setup(frames, s, frame_count)
            next_frame = 1
        end if
        call em_slab_magnetic_now(s, bn)
        field_start = em_slab_field_energy_range(s, bn, ip, s%cells - 1)
        kinetic_start = 0
        f_prev = 0; b_prev = 0; out_prev = 0
        do n = 0, total_steps
            time_now = n*s%dt
            call em_slab_magnetic_now(s, bn)
            call em_slab_characteristics_at(s, bn, ip, fwd, bwd)
            call em_slab_characteristics_at(s, bn, ip2, fwd2, bwd2)
            s_out = s%ey(s%cells - 1)*(bn(s%cells - 2) + bn(s%cells - 1))/2
            f_in = fwd**2; b_in = bwd**2
            if (n > 0) then
                result%incident = result%incident + s%dt*(f_prev + f_in)/2
                result%reflected = result%reflected + s%dt*(b_prev + b_in)/2
                result%outflow = result%outflow + s%dt*(out_prev + s_out)/2
            end if
            f_prev = f_in; b_prev = b_in; out_prev = s_out
            if (n >= first_step .and. n < total_steps) then
                ! Lock-in a = (2/T) sum f e^(i omega t) dt over whole periods.
                weight = 2*s%dt/(window_steps*s%dt)
                rotor = exp(cmplx(0, config%omega*time_now, dp))
                result%amplitude = result%amplitude + weight*rotor*s%ey
                result%forward = result%forward + weight*rotor*fwd
                result%backward = result%backward + weight*rotor*bwd
                result%forward_far = result%forward_far + weight*rotor*fwd2
            end if
            if (present(frames)) then
                if (next_frame <= size(frames%t)) then
                    if (time_now >= frames%t(next_frame) - s%dt/2) then
                        call frame_record(frames, s, next_frame)
                        next_frame = next_frame + 1
                    end if
                end if
            end if
            if (n == total_steps) exit
            call em_slab_step(s)
            if (size(s%x) > 0) result%max_momentum = max(result%max_momentum, &
                maxval(abs(s%ux)), maxval(abs(s%uy)))
        end do
        call em_slab_magnetic_now(s, bn)
        energy_now = em_slab_field_energy_range(s, bn, ip, s%cells - 1)
        result%field_change = energy_now - field_start
        ! Kinetic energy at the final level: one more half-level average.
        result%kinetic_change = final_kinetic(s) - kinetic_start
        result%balance = (result%incident - result%reflected - result%outflow - &
                          result%field_change - result%kinetic_change)/result%incident
        result%gauss_residual = em_slab_gauss_residual(s)
    contains
        !> Kinetic energy at the final integer level, centred like s%kinetic:
        !> half-level average using the push that would follow.
        real(dp) function final_kinetic(state) result(k)
            type(em_slab_t), intent(in) :: state
            type(em_slab_t) :: probe_state
            probe_state = state
            call em_slab_step(probe_state)
            k = probe_state%kinetic
        end function
    end subroutine

    !> Complex reflection coefficient referred to the plasma edge from the
    !> probe characteristics: backward/forward = R exp(-2 i omega (x_p - x_e)).
    complex(dp) function em_slab_reflection(res, omega, x_edge) result(r)
        type(em_slab_result_t), intent(in) :: res
        real(dp), intent(in) :: omega, x_edge
        r = res%backward/res%forward*exp(cmplx(0, 2*omega*(res%probe - x_edge), dp))
    end function

    !> Decay length -1/slope of a least-squares fit of log|a(x)| on [x0, x1].
    real(dp) function em_slab_decay_length(res, x0, x1) result(length)
        type(em_slab_result_t), intent(in) :: res
        real(dp), intent(in) :: x0, x1
        integer :: i0, i1
        i0 = result_node(res, x0); i1 = result_node(res, x1)
        length = -1/fit_slope(res%x(i0:i1), log(abs(res%amplitude(i0:i1))))
    end function

    !> Mean lock-in |a(x)| on [x0, x1].
    real(dp) function em_slab_mean_amplitude(res, x0, x1) result(value)
        type(em_slab_result_t), intent(in) :: res
        real(dp), intent(in) :: x0, x1
        integer :: i0, i1
        i0 = result_node(res, x0); i1 = result_node(res, x1)
        value = sum(abs(res%amplitude(i0:i1)))/(i1 - i0 + 1)
    end function

    !> Wave number from a fit of the unwrapped lock-in phase on [x0, x1].
    real(dp) function em_slab_wavenumber(res, x0, x1) result(k)
        type(em_slab_result_t), intent(in) :: res
        real(dp), intent(in) :: x0, x1
        real(dp), allocatable :: phase(:)
        complex(dp) :: ratio
        integer :: i0, i1, i
        i0 = result_node(res, x0); i1 = result_node(res, x1)
        allocate (phase(i1 - i0 + 1))
        phase(1) = atan2(aimag(res%amplitude(i0)), real(res%amplitude(i0)))
        do i = 2, size(phase)
            ratio = res%amplitude(i0 + i - 1)/res%amplitude(i0 + i - 2)
            phase(i) = phase(i - 1) + atan2(aimag(ratio), real(ratio))
        end do
        k = fit_slope(res%x(i0:i1), phase)
    end function

    integer function result_node(res, x) result(i)
        type(em_slab_result_t), intent(in) :: res
        real(dp), intent(in) :: x
        integer :: first
        first = lbound(res%x, 1)
        i = first + nint((x - res%x(first))/(res%x(first + 1) - res%x(first)))
        if (i < first .or. i > ubound(res%x, 1)) error stop 'em_slab: x outside grid'
    end function

    real(dp) function fit_slope(xs, ys) result(m)
        real(dp), intent(in) :: xs(:), ys(:)
        real(dp) :: xm, ym
        xm = sum(xs)/size(xs); ym = sum(ys)/size(ys)
        m = sum((xs - xm)*(ys - ym))/sum((xs - xm)**2)
    end function

    subroutine frame_setup(frames, s, count)
        type(em_slab_frames_t), intent(inout) :: frames
        type(em_slab_t), intent(in) :: s
        integer, intent(in) :: count
        integer :: lo, hi, i, k, m
        lo = nint((frames%x_min - s%config%x_left)/s%config%dx)
        hi = nint((frames%x_max - s%config%x_left)/s%config%dx)
        m = (hi - lo)/frames%stride + 1
        if (allocated(frames%t)) deallocate (frames%t, frames%x, frames%ey, &
                                             frames%x0, frames%px, frames%py)
        allocate (frames%t(count), frames%x(m), frames%ey(m, count))
        frames%t = [(i*frames%interval, i=0, count - 1)]
        frames%x = [(s%config%x_left + (lo + i*frames%stride)*s%config%dx, &
                     i=0, m - 1)]
        m = 0
        do k = 1, size(s%x0), frames%particle_stride
            if (s%x0(k) >= frames%particle_min .and. &
                s%x0(k) <= frames%particle_max) &
                m = m + 1
        end do
        allocate (frames%x0(m), frames%px(m, count), frames%py(m, count))
        m = 0
        do k = 1, size(s%x0), frames%particle_stride
            if (s%x0(k) >= frames%particle_min .and. &
                s%x0(k) <= frames%particle_max) then
                m = m + 1; frames%x0(m) = s%x0(k)
            end if
        end do
    end subroutine

    subroutine frame_record(frames, s, index)
        type(em_slab_frames_t), intent(inout) :: frames
        type(em_slab_t), intent(in) :: s
        integer, intent(in) :: index
        integer :: lo, i, k, m
        frames%t(index) = s%time
        lo = nint((frames%x_min - s%config%x_left)/s%config%dx)
        do i = 1, size(frames%x)
            frames%ey(i, index) = s%ey(lo + (i - 1)*frames%stride)
        end do
        m = 0
        do k = 1, size(s%x0), frames%particle_stride
            if (s%x0(k) >= frames%particle_min .and. &
                s%x0(k) <= frames%particle_max) then
                m = m + 1
                frames%px(m, index) = s%x(k); frames%py(m, index) = s%y(k)
            end if
        end do
    end subroutine
end module kin6d_em_slab
