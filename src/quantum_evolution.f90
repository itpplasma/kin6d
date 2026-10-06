module kin6d_quantum_evolution
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortnum_interval, only: interval_t, ipoint, interval, mag, intersect, entire
    use fortnum_ode_cash_karp, only: cash_karp_step
    use fortnum_rounding, only: add_up, mul_up, mul_down, sub_up, gamma_up
    use fortnum_verified_matvec, only: matvec_bound_t, prepare_matvec_bound, &
        enclose_matvec, matvec_row_bound
    use fortnum_rng, only: rng_t
    use fortnum_sampled_norm, only: sampled_norm_upper
    use fortnum_ode_residual_certificate, only: hermite_cubic_enclosure, &
        linear_residual_bound, residual_radius_step
    use kin6d_quantum_evolution_leaves, only: qe_rhs_term, qe_rhs_term_interval, &
        qe_rhs_term_vjp, qe_add, qe_add_interval, qe_scale, qe_hermite, &
        qe_evaluate, qe_evaluate_interval, qe_square_add_interval, &
        qe_norm_interval, qe_radius, qe_band_interval, &
        qe_fraction, qe_fraction_interval, qe_complement, &
        qe_action_allowance_interval, &
        qe_evaluate_rate_interval
    implicit none
    private
    public :: quantum_evolution_t, initialize_quantum_evolution
    public :: solve_quantum_evolution, sample_quantum_evolution
    public :: sample_quantum_evolution_box, quantum_evolution_info
    public :: quantum_evolution_radius, quantum_evolution_accepted
    public :: solve_quantum_evolution_tangent, sample_quantum_evolution_tangent
    public :: quantum_evolution_vjp
    public :: quantum_hamiltonian_action, quantum_rhs_action
    public :: certify_quantum_evolution, quantum_evolution_error, &
        quantum_evolution_indicator, quantum_evolution_sampled_residual
    integer, parameter, public :: QE_WORKHORSE = 1, QE_STRICT = 2, QE_APPROXIMATE = 3

    type :: quantum_evolution_t
        private
        integer :: n = 0, steps = 0, active_cell = 0, mode = 0, action_mode = 0
        real(dp) :: h = 0.0_dp, hbar = 0.0_dp, budget = 0.0_dp
        real(dp), allocatable :: hr(:, :), hi(:, :), y(:, :), f(:, :)
        real(dp), allocatable :: c(:, :, :), rk(:, :), dhr(:, :), dhi(:, :)
        real(dp), allocatable :: dy(:, :), df(:, :), dc(:, :, :)
        type(interval_t), allocatable :: ci(:, :, :), dci(:, :, :)
        type(interval_t), allocatable :: action(:, :), forcing(:, :)
        real(dp), allocatable :: radius(:), dradius(:)
        type(interval_t) :: initial_norm, delta_norm
        logical :: initialized = .false., ready = .false., bounded = .false.
        logical :: tangent_ready = .false., tangent_bounded = .false.
        integer :: error_mode = QE_WORKHORSE
        logical :: approximate_ready = .false.
        logical :: error_attempted = .false.
        real(dp) :: embedded_indicator = 0.0_dp
        real(dp), allocatable :: evaluation_allowance(:), cell_residual(:)
        real(dp), allocatable :: cell_allowance(:)
        real(dp) :: bound_relative_error = .05_dp
        type(matvec_bound_t) :: hr_bound, hi_bound, dhr_bound, dhi_bound
    end type
    type :: vjp_context_t
        type(quantum_evolution_t), pointer :: primal => null(), adjoint => null()
        integer :: reverse_cell = 0
    end type
contains
    subroutine initialize_quantum_evolution(state, hre, him, psi_re, psi_im, &
                                            steps, h, hbar, ok, budget)
        type(quantum_evolution_t), intent(out) :: state
        real(dp), intent(in) :: hre(:, :), him(:, :), psi_re(:), psi_im(:)
        integer, intent(in) :: steps
        real(dp), intent(in) :: h, hbar
        logical, intent(out) :: ok
        real(dp), intent(in), optional :: budget
        type(interval_t) :: acc, next
        integer :: n, j, status

        ok = .false.
        n = size(psi_re)
        if (n < 1 .or. steps < 1) return
        if (n > shiftr(huge(n), 1)) return
        if (n > shiftr(huge(n),3)-1) return
        if (steps >= huge(steps)) return
        if (size(psi_im) /= n) return
        if (.not. valid_hamiltonian(hre, him, n)) return
        if (.not. all(ieee_is_finite(psi_re))) return
        if (.not. all(ieee_is_finite(psi_im))) return
        if (.not. ieee_is_finite(h) .or. .not. ieee_is_finite(hbar)) return
        if (h <= 0.0_dp .or. hbar <= 0.0_dp) return
        if (present(budget)) then
            if (.not. ieee_is_finite(budget)) return
            if (budget < 0.0_dp) return
        end if
        allocate (state%hr(n, n), state%hi(n, n), state%y(2*n, 0:steps), &
                  state%f(2*n, 0:steps), state%c(2*n, 0:3, steps), &
                  state%ci(2*n, 0:3, steps), state%rk(2*n, 9), &
                  state%action(2*n, 0:3), &
                  state%radius(0:steps), stat=status)
        if (status /= 0) return
        allocate(state%evaluation_allowance(0:steps),state%cell_residual(steps), &
            state%cell_allowance(steps))
        state%evaluation_allowance = 0
        state%cell_allowance = 0
        state%n = n; state%steps = steps
        state%h = h; state%hbar = hbar
        state%hr = hre; state%hi = him
        call prepare_matvec_bound(hre, state%hr_bound, ok)
        if (.not. ok) return
        call prepare_matvec_bound(him, state%hi_bound, ok)
        if (.not. ok) return
        ok = .false.
        state%y(1:n, 0) = psi_re; state%y(n + 1:2*n, 0) = psi_im
        acc = ipoint(0.0_dp)
        do j = 1, 2*n
            call qe_square_add_interval(acc, ipoint(state%y(j, 0)), next)
            acc = next
        end do
        call squared_norm(acc, state%initial_norm)
        if (.not. finite_interval(state%initial_norm)) return
        call qe_radius(1.e-12_dp, 0.01_dp, state%initial_norm%hi, &
                       0.0_dp, 0.0_dp, state%budget)
        if (present(budget)) state%budget = budget
        if (.not. ieee_is_finite(state%budget)) return
        state%initialized = .true.
        ok = .true.
    end subroutine

    subroutine solve_quantum_evolution(state, ok, certify, error_mode)
        type(quantum_evolution_t), intent(inout) :: state
        logical, intent(out) :: ok
        logical, intent(in), optional :: certify
        integer, intent(in), optional :: error_mode
        integer :: cell, nfev
        logical :: check

        ok = .false.; state%ready = .false.; state%bounded = .false.
        state%tangent_ready = .false.; state%tangent_bounded = .false.
        state%approximate_ready = .false.; state%embedded_indicator = 0
        state%error_mode = QE_WORKHORSE
        if (present(error_mode)) state%error_mode = error_mode
        if (state%error_mode < QE_WORKHORSE .or. &
            state%error_mode > QE_APPROXIMATE) return
        if (.not. state%initialized) return
        check = .true.
        if (present(certify)) check = certify
        state%error_attempted = check
        state%mode = 0
        call evolution_rhs(0.0_dp, state%y(:, 0), state%f(:, 0), state)
        nfev = 0
        do cell = 1, state%steps
            call cash_karp_step(evolution_rhs, 0.0_dp, state%y(:, cell - 1), &
                                state%h, .false., state%rk(:, 1), state%rk(:, 2), &
                       state%rk(:, 3), state%rk(:, 4), state%rk(:, 5), state%rk(:, 6), &
                            state%rk(:, 7), state%rk(:, 8), state%rk(:, 9), nfev, state)
            if (.not. all(ieee_is_finite(state%rk(:, 8)))) return
            state%embedded_indicator = max(state%embedded_indicator, &
                norm2(state%rk(:, 9)))
            state%y(:, cell) = state%rk(:, 8)
            call evolution_rhs(0.0_dp, state%y(:, cell), state%f(:, cell), state)
            if (.not. all(ieee_is_finite(state%f(:, cell)))) return
            call make_cell(state, cell, .false., check)
        end do
        state%ready = .true.
        if (check) call certify_trajectory(state, .false.)
        ok = .true.
    end subroutine

    ! Deferred certification audits every residual cell of the immutable trace;
    ! reporting only at checkpoints does not omit the intervening integral.
    subroutine certify_quantum_evolution(state, ok, error_mode, epsilon_b)
        type(quantum_evolution_t), intent(inout) :: state
        logical, intent(out) :: ok
        integer, intent(in), optional :: error_mode
        real(dp), intent(in), optional :: epsilon_b
        integer :: cell
        ok = .false.
        state%bounded = .false.; state%approximate_ready = .false.
        state%error_attempted = .true.
        state%tangent_bounded = .false.
        if (.not. state%ready) return
        if (present(epsilon_b)) state%bound_relative_error = epsilon_b
        if (.not. ieee_is_finite(state%bound_relative_error)) return
        if (state%bound_relative_error <= 0 .or. &
            state%bound_relative_error > .1_dp) return
        if (present(error_mode)) state%error_mode = error_mode
        if (state%error_mode < QE_WORKHORSE .or. &
            state%error_mode > QE_APPROXIMATE) return
        do cell = 1, state%steps
            call make_cell(state, cell, .false., .true.)
        end do
        call certify_trajectory(state, .false.)
        ok = state%bounded .or. state%approximate_ready
    end subroutine

    subroutine quantum_evolution_indicator(state, value, status)
        type(quantum_evolution_t), intent(in) :: state
        real(dp), intent(out) :: value
        character(*), intent(out) :: status
        value = state%embedded_indicator; status = 'indicator'
        if (.not. state%ready) status = 'unavailable'
    end subroutine

    ! A one-time residual-norm statement, not the integral needed for a state
    ! certificate. The probability space is iid uniform coordinate draws.
    subroutine quantum_evolution_sampled_residual(state, cell, fraction, samples, &
        delta, draw, estimate, upper, ok)
        type(quantum_evolution_t), intent(inout) :: state
        integer, intent(in) :: cell, samples
        real(dp), intent(in) :: fraction, delta
        type(rng_t), intent(inout) :: draw
        real(dp), intent(out) :: estimate, upper
        logical, intent(out) :: ok
        type(interval_t), allocatable :: x(:), derivative(:)
        type(interval_t) :: envelope_box
        real(dp) :: envelope
        integer :: j,n
        ok = .false.; estimate = 0; upper = huge(0._dp)
        n = state%n
        if (.not. sample_guard(state,cell,fraction,n,n)) return
        call make_cell(state,cell,.false.,.true.)
        allocate(x(2*n),derivative(2*n))
        do j = 1,2*n
            call qe_evaluate_interval(state%ci(j,0,cell),state%ci(j,1,cell), &
                state%ci(j,2,cell),state%ci(j,3,cell),ipoint(fraction),x(j))
            call qe_evaluate_rate_interval(state%ci(j,1,cell),state%ci(j,2,cell), &
                state%ci(j,3,cell),ipoint(fraction),ipoint(state%h),derivative(j))
        end do
        ! A global envelope proved from prepared matrix row sums, not sampled maxima.
        call qe_action_allowance_interval(ipoint(matvec_row_bound(state%hr_bound)), &
            ipoint(matvec_row_bound(state%hi_bound)),ipoint(0._dp),ipoint(0._dp), &
            ipoint(maxval(mag(x))),ipoint(0._dp),ipoint(0._dp),ipoint(1._dp), &
            ipoint(state%hbar),envelope_box)
        envelope = add_up(envelope_box%hi,maxval(mag(derivative)))
        call sampled_norm_upper(2*n,samples,delta,envelope,draw,component, &
            estimate,upper,ok)
    contains
        subroutine component(index,value,valid)
            use fortnum_interval, only: operator(-)
            integer, intent(in) :: index
            type(interval_t), intent(out) :: value
            logical, intent(out) :: valid
            type(interval_t) :: ar,ai,total,next
            integer :: i,k
            i = index
            if (index > n) i = index-n
            total = ipoint(0._dp)
            do k = 1,n
                call qe_rhs_term_interval(ipoint(state%hr(i,k)),ipoint(state%hi(i,k)), &
                    x(k),x(n+k),ipoint(state%hbar),ar,ai)
                if (index <= n) then
                    call qe_add_interval(total,ar,next)
                else
                    call qe_add_interval(total,ai,next)
                end if
                total = next
            end do
            value = derivative(index)-total
            valid = finite_interval(value)
        end subroutine
    end subroutine

    ! The approximate action omits proved evaluation allowances. Each prefix
    ! must pass allowance <= epsilon_b*eta. The uninflated eta is not certified.
    subroutine quantum_evolution_error(state, cell, radius, status, epsilon_b, ok, &
                                      fraction, evaluation_error)
        type(quantum_evolution_t), intent(in) :: state
        integer, intent(in) :: cell
        real(dp), intent(out) :: radius, epsilon_b
        character(*), intent(out) :: status
        logical, intent(out) :: ok
        real(dp), intent(in), optional :: fraction
        real(dp), intent(out), optional :: evaluation_error
        real(dp) :: re(state%n), im(state%n), at, unused
        logical :: sampled, bounded
        radius = huge(0._dp); epsilon_b = 0; status = 'candidate'; ok = .false.
        if (present(evaluation_error)) evaluation_error = huge(0._dp)
        if (.not. state%ready) return
        if (cell < 0 .or. cell > state%steps) return
        if (.not. state%bounded .and. .not. state%approximate_ready) then
            if (state%error_attempted) status = 'certification-failed'
            return
        end if
        radius = state%radius(cell)
        if (cell > 0) then
            at = 1
            if (present(fraction)) at = fraction
            call sample_quantum_evolution(state,cell,at,re,im,unused,sampled,bounded)
            if (.not. sampled) return
            call sampled_radius(state,cell,at,re,im,.false.,radius,ok)
            if (.not. ok) return
        end if
        status = 'certified-discrete'
        if (state%approximate_ready) then
            epsilon_b = state%bound_relative_error
            write(status,'(a,es10.3,a)') 'approximate-bound; epsilon_b=', &
                epsilon_b,'; omission-controlled'
        end if
        if (present(evaluation_error)) then
            evaluation_error = state%evaluation_allowance(cell)
            if (cell > 0) evaluation_error = &
                add_up(state%evaluation_allowance(cell-1), &
                    mul_up(mul_up(state%h,at),state%cell_allowance(cell)))
        end if
        ok = .true.
    end subroutine

    subroutine evolution_rhs(t, y, dydt, ctx)
        real(dp), intent(in) :: t, y(:)
        real(dp), intent(out) :: dydt(:)
        class(*), intent(in), optional :: ctx
        real(dp) :: a, b, next, fraction, xr, xi, da, db
        integer :: i, j, n, cell

        dydt = 0.0_dp
        if (.not. ieee_is_finite(t)) return
        if (.not. present(ctx)) return
        select type (ctx)
        type is (quantum_evolution_t)
            n = ctx%n
            if (size(y) /= 2*n .or. size(dydt) /= 2*n) return
            cell = ctx%active_cell
            if (ctx%mode == 1) call qe_fraction(t, ctx%h, fraction)
            do i = 1, n
                do j = 1, n
                    call qe_rhs_term(ctx%hr(i, j), ctx%hi(i, j), &
                                     y(j), y(n + j), ctx%hbar, a, b)
                    if (ctx%mode == 1) then
                        call qe_evaluate(ctx%c(j, 0, cell), ctx%c(j, 1, cell), &
                                     ctx%c(j, 2, cell), ctx%c(j, 3, cell), fraction, xr)
                        call qe_evaluate(ctx%c(n + j, 0, cell), &
                                         ctx%c(n + j, 1, cell), ctx%c(n + j, 2, cell), &
                                         ctx%c(n + j, 3, cell), fraction, xi)
                        call qe_rhs_term(ctx%dhr(i, j), ctx%dhi(i, j), &
                                         xr, xi, ctx%hbar, da, db)
                        call qe_add(a, da, next); a = next
                        call qe_add(b, db, next); b = next
                    end if
                    call qe_add(dydt(i), a, next); dydt(i) = next
                    call qe_add(dydt(n + i), b, next); dydt(n + i) = next
                end do
            end do
        end select
    end subroutine

    subroutine make_cell(state, cell, tangent, strict)
        type(quantum_evolution_t), intent(inout) :: state
        integer, intent(in) :: cell
        logical, intent(in) :: tangent, strict
        integer :: j

        do j = 1, 2*state%n
            if (tangent) then
                call qe_hermite(state%dy(j, cell - 1), state%dy(j, cell), &
                                state%df(j, cell - 1), state%df(j, cell), state%h, &
                                state%dc(j, 0, cell), state%dc(j, 1, cell), &
                                state%dc(j, 2, cell), state%dc(j, 3, cell))
            else
                call qe_hermite(state%y(j, cell - 1), state%y(j, cell), &
                                state%f(j, cell - 1), state%f(j, cell), state%h, &
                                state%c(j, 0, cell), state%c(j, 1, cell), &
                                state%c(j, 2, cell), state%c(j, 3, cell))
            end if
        end do
        if (.not. strict) return
        if (tangent) then
            call hermite_cubic_enclosure(state%dy(:, cell - 1), state%dy(:, cell), &
                                         state%df(:, cell - 1), state%df(:, cell), &
                                         ipoint(state%h), state%dci(:, :, cell))
        else
            call hermite_cubic_enclosure(state%y(:, cell - 1), state%y(:, cell), &
                                         state%f(:, cell - 1), state%f(:, cell), &
                                         ipoint(state%h), state%ci(:, :, cell))
        end if
    end subroutine

    ! FortNum's Hermite/Bernstein residual and radius accumulation on the
    ! kin6d interval action -i H x/hbar. The generator is skew-adjoint in the
    ! packed real/imag Hilbert pairing, so the logarithmic norm is zero. The
    ! tangent residual adds the forcing -i dH y_h/hbar and the discrepancy
    ! ||dH||_F/hbar * radius(cell) of the primal reconstruction.
    subroutine certify_trajectory(state, tangent)
        type(quantum_evolution_t), intent(inout) :: state
        logical, intent(in) :: tangent
        real(dp) :: bound, next, omitted, xmax, width, row_r, row_i
        type(interval_t) :: allowance
        integer :: cell, k
        logical :: ok

        if (tangent) then
            state%dradius(0) = 0.0_dp
            state%tangent_bounded = .false.
        else
            state%radius(0) = 0.0_dp
            state%bounded = .false.
            state%approximate_ready = .false.
            state%evaluation_allowance = 0
            state%cell_allowance = 0
        end if
        do cell = 1, state%steps
            if (tangent) then
                state%action_mode = 1
                do k = 0, 3
                    call qe_action(state%ci(:, k, cell), state%forcing(:, k), ok, state)
                    if (.not. ok) return
                end do
                state%action_mode = 0
                call linear_residual_bound(qe_action, state%dci(:, :, cell), &
                                           ipoint(state%h), state%action, bound, ok, &
                                           state, state%forcing)
                if (.not. ok) return
                bound = add_up(bound, mul_up(state%delta_norm%hi, state%radius(cell)))
                call residual_radius_step(ipoint(0.0_dp), state%dradius(cell - 1), &
                                          ipoint(state%h), bound, next, ok)
                if (.not. ok) return
                state%dradius(cell) = next
            else
                state%action_mode = 0
                call linear_residual_bound(qe_action, state%ci(:, :, cell), &
                                           ipoint(state%h), state%action, bound, ok, state)
                if (.not. ok) return
                call residual_radius_step(ipoint(0.0_dp), state%radius(cell - 1), &
                                          ipoint(state%h), bound, next, ok)
                if (.not. ok) return
                state%radius(cell) = next
                state%cell_residual(cell) = bound
                if (state%error_mode == QE_APPROXIMATE) then
                    omitted = 0
                    row_r = matvec_row_bound(state%hr_bound)
                    row_i = matvec_row_bound(state%hi_bound)
                    do k = 0,3
                        xmax = maxval(abs(state%ci(:,k,cell)%lo))
                        width = maxval(sub_up(state%ci(:,k,cell)%hi, &
                            state%ci(:,k,cell)%lo))
                        call qe_action_allowance_interval( &
                            ipoint(row_r),ipoint(row_i), &
                            ipoint(gamma_up(8*state%n+4)),ipoint(xmax),ipoint(width), &
                            ipoint(real(8*state%n+4,dp)),ipoint(tiny(0._dp)), &
                            ipoint(real(2*state%n,dp)), &
                            ipoint(state%hbar),allowance)
                        if (.not. finite_interval(allowance)) return
                        omitted = add_up(omitted,allowance%hi)
                    end do
                    state%evaluation_allowance(cell) = &
                        add_up(state%evaluation_allowance(cell-1), &
                               mul_up(state%h,omitted))
                    state%cell_allowance(cell) = omitted
                    if (state%evaluation_allowance(cell) > &
                        mul_down(state%bound_relative_error,next)) return
                end if
            end if
        end do
        if (tangent) then
            state%tangent_bounded = .true.
        else
            state%bounded = state%error_mode /= QE_APPROXIMATE
            state%approximate_ready = state%error_mode == QE_APPROXIMATE
        end if
    end subroutine

    ! Interval action x -> -i H x/hbar (action_mode 0) or -i dH x/hbar
    ! (action_mode 1) over the packed real/imag split.
    subroutine qe_action(x, ax, ok, ctx)
        type(interval_t), intent(in) :: x(:)
        type(interval_t), intent(out) :: ax(:)
        logical, intent(out) :: ok
        class(*), intent(in), optional :: ctx
        type(interval_t) :: a, b, next
        real(dp) :: ar, ai, sr, si, nc
        integer :: i, j, n

        ok = .false.
        ax = ipoint(0.0_dp)
        if (.not. present(ctx)) return
        select type (ctx)
        type is (quantum_evolution_t)
            n = ctx%n
            if (size(x) /= 2*n .or. size(ax) /= 2*n) return
            if (ctx%error_mode == QE_WORKHORSE) then
                call bounded_action(ctx, x, ax, ok)
                return
            end if
            if (ctx%error_mode == QE_APPROXIMATE) then
                ! Unenclosed action at the lower coefficient endpoints. The
                ! omitted box widths and FP error distinguish A from D/S.
                do i = 1,n
                    sr = 0; si = 0
                    do j = 1,n
                        call qe_rhs_term(ctx%hr(i,j), ctx%hi(i,j), &
                            x(j)%lo, &
                            x(n+j)%lo, ctx%hbar, ar, ai)
                        call qe_add(sr,ar,nc); sr = nc
                        call qe_add(si,ai,nc); si = nc
                    end do
                    ax(i) = ipoint(sr); ax(n+i) = ipoint(si)
                end do
                ok = all_finite_intervals(ax)
                return
            end if
            do j = 1, n
                do i = 1, n
                    if (ctx%action_mode == 1) then
                        call qe_rhs_term_interval(ipoint(ctx%dhr(i, j)), &
                                                  ipoint(ctx%dhi(i, j)), x(j), x(n + j), &
                                                  ipoint(ctx%hbar), a, b)
                    else
                        call qe_rhs_term_interval(ipoint(ctx%hr(i, j)), &
                                                  ipoint(ctx%hi(i, j)), x(j), x(n + j), &
                                                  ipoint(ctx%hbar), a, b)
                    end if
                    call qe_add_interval(ax(i), a, next); ax(i) = next
                    call qe_add_interval(ax(n + i), b, next); ax(n + i) = next
                end do
            end do
            ok = all_finite_intervals(ax)
        end select
    end subroutine

    subroutine bounded_action(state, x, ax, ok)
        type(quantum_evolution_t), intent(in) :: state
        type(interval_t), intent(in) :: x(:)
        type(interval_t), intent(inout) :: ax(:)
        logical, intent(out) :: ok
        integer :: n,j
        type(interval_t) :: next
        n = state%n; ok = .false.
        if (state%action_mode == 1) then
            call enclose_matvec(state%dhr,state%dhr_bound,x(n+1:),ax(:n),ok)
            if (.not. ok) return
            call enclose_matvec(state%dhr,state%dhr_bound,x(:n),ax(n+1:),ok,sign=-1)
            if (.not. ok) return
            call enclose_matvec(state%dhi,state%dhi_bound,x(:n),ax(:n),ok, &
                accumulate=.true.)
            if (.not. ok) return
            call enclose_matvec(state%dhi,state%dhi_bound,x(n+1:),ax(n+1:),ok, &
                accumulate=.true.)
        else
            call enclose_matvec(state%hr,state%hr_bound,x(n+1:),ax(:n),ok)
            if (.not. ok) return
            call enclose_matvec(state%hr,state%hr_bound,x(:n),ax(n+1:),ok,sign=-1)
            if (.not. ok) return
            call enclose_matvec(state%hi,state%hi_bound,x(:n),ax(:n),ok, &
                accumulate=.true.)
            if (.not. ok) return
            call enclose_matvec(state%hi,state%hi_bound,x(n+1:),ax(n+1:),ok, &
                accumulate=.true.)
        end if
        if (.not. ok) return
        do j = 1,2*n
            call qe_fraction_interval(ax(j),ipoint(state%hbar),next)
            ax(j) = next
        end do
        ok = all_finite_intervals(ax)
    end subroutine

    subroutine sample_quantum_evolution(state, cell, fraction, re, im, radius, &
                                        ok, bounded)
        type(quantum_evolution_t), intent(in) :: state
        integer, intent(in) :: cell
        real(dp), intent(in) :: fraction
        real(dp), intent(out) :: re(:), im(:), radius
        logical, intent(out) :: ok, bounded
        integer :: j, n

        ok = .false.; bounded = .false.; radius = huge(0.0_dp)
        if (.not. sample_guard(state, cell, fraction, size(re), size(im))) return
        n = state%n
        do j = 1, n
            call qe_evaluate(state%c(j, 0, cell), state%c(j, 1, cell), &
                             state%c(j, 2, cell), state%c(j, 3, cell), fraction, re(j))
            call qe_evaluate(state%c(n + j, 0, cell), state%c(n + j, 1, cell), &
                      state%c(n + j, 2, cell), state%c(n + j, 3, cell), fraction, im(j))
        end do
        if (.not. all(ieee_is_finite(re))) return
        if (.not. all(ieee_is_finite(im))) return
        ok = .true.; bounded = state%bounded
        if (bounded) then
            call sampled_radius(state, cell, fraction, re, im, .false., radius, &
                                bounded)
        end if
    end subroutine

    subroutine sampled_radius(state, cell, fraction, re, im, tangent, radius, ok)
        type(quantum_evolution_t), intent(in) :: state
        integer, intent(in) :: cell
        real(dp), intent(in) :: fraction, re(:), im(:)
        logical, intent(in) :: tangent
        real(dp), intent(out) :: radius
        logical, intent(out) :: ok
        type(interval_t) :: yr, yi, error, sumsq, next, norm, total
        real(dp) :: partial
        integer :: j, n

        n = state%n; sumsq = ipoint(0.0_dp); ok = .false.
        radius = huge(0.0_dp)
        do j = 1, n
            call reconstructed_point(state, j, cell, fraction, tangent, yr)
            call reconstructed_point(state, n + j, cell, fraction, tangent, yi)
            call qe_add_interval(yr, ipoint(-re(j)), error)
            call qe_square_add_interval(sumsq, ipoint(mag(error)), next)
            sumsq = next
            call qe_add_interval(yi, ipoint(-im(j)), error)
            call qe_square_add_interval(sumsq, ipoint(mag(error)), next)
            sumsq = next
        end do
        call squared_norm(sumsq, norm)
        if (tangent) then
            call qe_add_interval(ipoint(state%dradius(cell)), norm, total)
        else
            partial = add_up(state%radius(cell-1), &
                mul_up(mul_up(state%h,fraction),state%cell_residual(cell)))
            call qe_add_interval(ipoint(partial), norm, total)
        end if
        if (.not. finite_interval(total)) return
        radius = total%hi; ok = .true.
    end subroutine

    subroutine reconstructed_point(state, index, cell, fraction, tangent, value)
        type(quantum_evolution_t), intent(in) :: state
        integer, intent(in) :: index, cell
        real(dp), intent(in) :: fraction
        logical, intent(in) :: tangent
        type(interval_t), intent(out) :: value

        if (tangent) then
            call qe_evaluate_interval(state%dci(index, 0, cell), &
                                 state%dci(index, 1, cell), state%dci(index, 2, cell), &
                                     state%dci(index, 3, cell), ipoint(fraction), value)
        else
            call qe_evaluate_interval(state%ci(index, 0, cell), &
                                   state%ci(index, 1, cell), state%ci(index, 2, cell), &
                                      state%ci(index, 3, cell), ipoint(fraction), value)
        end if
    end subroutine

    subroutine sample_quantum_evolution_box(state, cell, fraction, re, im, &
                                            radius, ok, bounded)
        type(quantum_evolution_t), intent(in) :: state
        integer, intent(in) :: cell
        real(dp), intent(in) :: fraction
        type(interval_t), intent(out) :: re(:), im(:)
        real(dp), intent(out) :: radius
        logical, intent(out) :: ok, bounded
        type(interval_t) :: value, lower, upper
        integer :: j, n

        ok = .false.; bounded = .false.; radius = huge(0.0_dp)
        if (.not. sample_guard(state, cell, fraction, size(re), size(im))) return
        if (.not. state%bounded) return
        n = state%n; radius = state%radius(cell)
        do j = 1, n
            call reconstructed_point(state, j, cell, fraction, .false., value)
            call qe_band_interval(value, ipoint(radius), lower, upper)
            re(j) = interval(lower%lo, upper%hi)
            call reconstructed_point(state, n + j, cell, fraction, .false., value)
            call qe_band_interval(value, ipoint(radius), lower, upper)
            im(j) = interval(lower%lo, upper%hi)
        end do
        ok = all_finite_intervals(re) .and. all_finite_intervals(im)
        bounded = ok
    end subroutine

    subroutine quantum_evolution_info(state, n, steps, h, hbar, ready, bounded)
        type(quantum_evolution_t), intent(in) :: state
        integer, intent(out) :: n, steps
        real(dp), intent(out) :: h, hbar
        logical, intent(out) :: ready, bounded
        n = state%n; steps = state%steps; h = state%h; hbar = state%hbar
        ready = state%ready; bounded = state%bounded
    end subroutine

    subroutine quantum_evolution_radius(state, radius, ok)
        type(quantum_evolution_t), intent(in) :: state
        real(dp), intent(out) :: radius
        logical, intent(out) :: ok
        radius = huge(0.0_dp); ok = state%bounded
        if (ok) radius = state%radius(state%steps)
    end subroutine

    logical function quantum_evolution_accepted(state)
        type(quantum_evolution_t), intent(in) :: state
        quantum_evolution_accepted = .false.
        if (.not. state%bounded) return
        quantum_evolution_accepted = state%radius(state%steps) <= state%budget
    end function

    subroutine solve_quantum_evolution_tangent(state, dhre, dhim, dpsi_re, &
                                               dpsi_im, ok, certify)
        type(quantum_evolution_t), intent(inout) :: state
        real(dp), intent(in) :: dhre(:, :), dhim(:, :), dpsi_re(:), dpsi_im(:)
        logical, intent(out) :: ok
        logical, intent(in), optional :: certify
        integer :: n, cell, i, j, nfev, status
        logical :: check
        type(interval_t) :: acc, next, norm

        ok = .false.; state%tangent_ready = .false.; state%tangent_bounded = .false.
        if (.not. state%ready) return
        n = state%n
        if (.not. valid_hamiltonian(dhre, dhim, n)) return
        if (size(dpsi_re) /= n .or. size(dpsi_im) /= n) return
        if (.not. all(ieee_is_finite(dpsi_re))) return
        if (.not. all(ieee_is_finite(dpsi_im))) return
        check = .true.
        if (present(certify)) check = certify
        if (check) then
            if (state%error_mode == QE_APPROXIMATE) return
            if (.not. state%bounded) then
                do cell = 1, state%steps
                    call make_cell(state, cell, .false., .true.)
                end do
                call certify_trajectory(state, .false.)
                if (.not. state%bounded) return
            end if
        end if
        if (.not. allocated(state%dy)) then
            allocate (state%dhr(n, n), state%dhi(n, n), &
                      state%dy(2*n, 0:state%steps), state%df(2*n, 0:state%steps), &
                state%dc(2*n, 0:3, state%steps), state%dci(2*n, 0:3, state%steps), &
                state%dradius(0:state%steps), state%forcing(2*n, 0:3), stat=status)
            if (status /= 0) then
                if (allocated(state%dhr)) deallocate (state%dhr)
                if (allocated(state%dhi)) deallocate (state%dhi)
                if (allocated(state%dy)) deallocate (state%dy)
                if (allocated(state%df)) deallocate (state%df)
                if (allocated(state%dc)) deallocate (state%dc)
                if (allocated(state%dci)) deallocate (state%dci)
                if (allocated(state%dradius)) deallocate (state%dradius)
                if (allocated(state%forcing)) deallocate (state%forcing)
                return
            end if
        end if
        state%dhr = dhre; state%dhi = dhim
        call prepare_matvec_bound(dhre, state%dhr_bound, ok)
        if (.not. ok) return
        call prepare_matvec_bound(dhim, state%dhi_bound, ok)
        if (.not. ok) return
        state%dy(1:n, 0) = dpsi_re; state%dy(n + 1:2*n, 0) = dpsi_im
        if (check) then
            acc = ipoint(0.0_dp)
            do j = 1, n
                do i = 1, n
                    call qe_square_add_interval(acc, ipoint(dhre(i, j)), next)
                    acc = next
                    call qe_square_add_interval(acc, ipoint(dhim(i, j)), next)
                    acc = next
                end do
            end do
            call squared_norm(acc, norm)
            call qe_fraction_interval(norm, ipoint(state%hbar), state%delta_norm)
            if (.not. finite_interval(state%delta_norm)) return
        end if
        state%mode = 1; state%active_cell = 1
        call evolution_rhs(0.0_dp, state%dy(:, 0), state%df(:, 0), state)
        nfev = 0
        do cell = 1, state%steps
            state%active_cell = cell
            call cash_karp_step(evolution_rhs, 0.0_dp, state%dy(:, cell - 1), &
                                state%h, .false., state%rk(:, 1), state%rk(:, 2), &
                       state%rk(:, 3), state%rk(:, 4), state%rk(:, 5), state%rk(:, 6), &
                            state%rk(:, 7), state%rk(:, 8), state%rk(:, 9), nfev, state)
            if (.not. all(ieee_is_finite(state%rk(:, 8)))) return
            state%dy(:, cell) = state%rk(:, 8)
            call evolution_rhs(state%h, state%dy(:, cell), state%df(:, cell), state)
            if (.not. all(ieee_is_finite(state%df(:, cell)))) return
            call make_cell(state, cell, .true., check)
        end do
        state%mode = 0; state%tangent_ready = .true.
        if (check) call certify_trajectory(state, .true.)
        ok = .true.
    end subroutine

    subroutine sample_quantum_evolution_tangent(state, cell, fraction, re, im, &
                                                radius, ok, bounded)
        type(quantum_evolution_t), intent(in) :: state
        integer, intent(in) :: cell
        real(dp), intent(in) :: fraction
        real(dp), intent(out) :: re(:), im(:), radius
        logical, intent(out) :: ok, bounded
        integer :: j, n

        ok = .false.; bounded = .false.; radius = huge(0.0_dp)
        if (.not. state%tangent_ready) return
        if (.not. sample_guard(state, cell, fraction, size(re), size(im))) return
        n = state%n
        do j = 1, n
            call qe_evaluate(state%dc(j, 0, cell), state%dc(j, 1, cell), &
                            state%dc(j, 2, cell), state%dc(j, 3, cell), fraction, re(j))
            call qe_evaluate(state%dc(n + j, 0, cell), state%dc(n + j, 1, cell), &
                    state%dc(n + j, 2, cell), state%dc(n + j, 3, cell), fraction, im(j))
        end do
        if (.not. all(ieee_is_finite(re))) return
        if (.not. all(ieee_is_finite(im))) return
        ok = .true.; bounded = state%tangent_bounded
        if (bounded) call sampled_radius(state, cell, fraction, re, im, .true., &
                                         radius, bounded)
    end subroutine

    subroutine quantum_rhs_action(hre, him, re, im, hbar, out_re, out_im, ok)
        real(dp), intent(in) :: hre(:, :), him(:, :), re(:), im(:), hbar
        real(dp), intent(out) :: out_re(:), out_im(:)
        logical, intent(out) :: ok
        real(dp) :: a, b, next
        integer :: n, i, j

        ok = .false.; n = size(re)
        if (size(im) /= n .or. size(out_re) /= n .or. size(out_im) /= n) return
        if (.not. valid_hamiltonian(hre, him, n)) return
        if (.not. all(ieee_is_finite(re)) .or. .not. all(ieee_is_finite(im))) return
        if (.not. ieee_is_finite(hbar)) return
        if (hbar <= 0.0_dp) return
        out_re = 0.0_dp; out_im = 0.0_dp
        do i = 1, n
            do j = 1, n
                call qe_rhs_term(hre(i, j), him(i, j), re(j), im(j), hbar, a, b)
                call qe_add(out_re(i), a, next); out_re(i) = next
                call qe_add(out_im(i), b, next); out_im(i) = next
            end do
        end do
        ok = all(ieee_is_finite(out_re)) .and. all(ieee_is_finite(out_im))
    end subroutine

    subroutine quantum_hamiltonian_action(hre, him, re, im, out_re, out_im, ok)
        real(dp), intent(in) :: hre(:, :), him(:, :), re(:), im(:)
        real(dp), intent(out) :: out_re(:), out_im(:)
        logical, intent(out) :: ok
        real(dp) :: a, b
        integer :: j

        call quantum_rhs_action(hre, him, re, im, 1.0_dp, out_re, out_im, ok)
        if (.not. ok) return
        do j = 1, size(re)
            a = out_re(j)
            call qe_scale(out_im(j), -1.0_dp, b)
            out_re(j) = b; out_im(j) = a
        end do
    end subroutine

    logical function sample_guard(state, cell, fraction, nr, ni)
        type(quantum_evolution_t), intent(in) :: state
        integer, intent(in) :: cell, nr, ni
        real(dp), intent(in) :: fraction
        sample_guard = .false.
        if (.not. state%ready) return
        if (cell < 1 .or. cell > state%steps) return
        if (nr /= state%n .or. ni /= state%n) return
        if (.not. ieee_is_finite(fraction)) return
        if (fraction < 0.0_dp .or. fraction > 1.0_dp) return
        sample_guard = .true.
    end function

    logical function valid_hamiltonian(hr, hi, n)
        real(dp), intent(in) :: hr(:, :), hi(:, :)
        integer, intent(in) :: n
        integer :: i, j
        real(dp) :: opposite
        valid_hamiltonian = .false.
        if (n < 1) return
        if (size(hr, 1) /= n .or. size(hr, 2) /= n) return
        if (size(hi, 1) /= n .or. size(hi, 2) /= n) return
        if (.not. all(ieee_is_finite(hr))) return
        if (.not. all(ieee_is_finite(hi))) return
        do j = 1, n
            do i = 1, n
                if (hr(i, j) < hr(j, i) .or. hr(i, j) > hr(j, i)) return
                call qe_scale(hi(j, i), -1.0_dp, opposite)
                if (hi(i, j) < opposite .or. hi(i, j) > opposite) return
            end do
        end do
        valid_hamiltonian = .true.
    end function

    logical function finite_interval(value)
        type(interval_t), intent(in) :: value
        finite_interval = ieee_is_finite(value%lo) .and. &
                          ieee_is_finite(value%hi) .and. value%lo <= value%hi
    end function

    logical function all_finite_intervals(values)
        type(interval_t), intent(in) :: values(:)
        all_finite_intervals = all(ieee_is_finite(values%lo)) .and. &
                        all(ieee_is_finite(values%hi)) .and. all(values%lo <= values%hi)
    end function

    subroutine squared_norm(sumsq, norm)
        type(interval_t), intent(in) :: sumsq
        type(interval_t), intent(out) :: norm
        ! Every caller constructed a sum of squares. Tighten only its domain;
        ! reject overflow before intersection can discard infinite endpoints.
        norm = entire()
        if (.not. finite_interval(sumsq)) return
        if (sumsq%hi < 0.0_dp) return
        call qe_norm_interval(intersect(sumsq, interval(0.0_dp, huge(0.0_dp))), norm)
    end subroutine

    subroutine quantum_evolution_vjp(state, terminal_re, terminal_im, &
                               initial_re, initial_im, ghre, ghim, initial_radius, ok, &
                                     initial_bounded, h_gradient_bounded)
        type(quantum_evolution_t), intent(in), target :: state
        real(dp), intent(in) :: terminal_re(:), terminal_im(:)
        real(dp), intent(out) :: initial_re(:), initial_im(:), ghre(:, :), ghim(:, :)
        real(dp), intent(out) :: initial_radius
        logical, intent(out) :: ok, initial_bounded, h_gradient_bounded
        type(quantum_evolution_t), target :: adjoint
        type(vjp_context_t) :: context
        real(dp), allocatable :: neg_hr(:, :), neg_hi(:, :), stages(:, :), gradient(:)
        integer :: n, length, nfev, reverse_cell, i, j, idx, status
        logical :: success

        ok = .false.; initial_bounded = .false.; h_gradient_bounded = .false.
        initial_radius = huge(0.0_dp)
        if (.not. state%ready) return
        n = state%n
        if (n > shiftr(huge(n)/n, 1)) return
        if (size(terminal_re) /= n .or. size(terminal_im) /= n) return
        if (size(initial_re) /= n .or. size(initial_im) /= n) return
        if (size(ghre, 1) /= n .or. size(ghre, 2) /= n) return
        if (size(ghim, 1) /= n .or. size(ghim, 2) /= n) return
        length = 2*n*n
        allocate (neg_hr(n, n), neg_hi(n, n), stages(length, 9), &
                  gradient(length), stat=status)
        if (status /= 0) return
        call qe_scale(state%hr, -1.0_dp, neg_hr)
        call qe_scale(state%hi, -1.0_dp, neg_hi)
        call initialize_quantum_evolution(adjoint, neg_hr, neg_hi, &
                    terminal_re, terminal_im, state%steps, state%h, state%hbar, success)
        if (.not. success) return
        call solve_quantum_evolution(adjoint, success)
        if (.not. success) return
        call sample_quantum_evolution(adjoint, state%steps, 1.0_dp, &
                       initial_re, initial_im, initial_radius, success, initial_bounded)
        if (.not. success) return
        gradient = 0.0_dp; nfev = 0
        context%primal => state; context%adjoint => adjoint
        do reverse_cell = 1, state%steps
            context%reverse_cell = reverse_cell
            call cash_karp_step(gradient_rhs, 0.0_dp, gradient, state%h, &
                      .false., stages(:, 1), stages(:, 2), stages(:, 3), stages(:, 4), &
                               stages(:, 5), stages(:, 6), stages(:, 7), stages(:, 8), &
                                stages(:, 9), nfev, context)
            if (.not. all(ieee_is_finite(stages(:, 8)))) return
            gradient = stages(:, 8)
        end do
        do j = 1, n
            do i = 1, n
                idx = i + (j - 1)*n
                ghre(i, j) = gradient(idx)
                ghim(i, j) = gradient(n*n + idx)
            end do
        end do
        ok = .true.
    end subroutine

    subroutine gradient_rhs(t, values, derivative, ctx)
        real(dp), intent(in) :: t, values(:)
        real(dp), intent(out) :: derivative(:)
        class(*), intent(in), optional :: ctx
        real(dp) :: s, forward_s, pr, pi, lr, li, gr, gi, gx, gy, gh
        integer :: i, j, idx, forward_cell, reverse_cell, n, length
        type(quantum_evolution_t), pointer :: state, adjoint

        derivative = 0.0_dp
        if (.not. present(ctx)) return
        select type (ctx)
        type is (vjp_context_t)
            if (.not. associated(ctx%primal)) return
            if (.not. associated(ctx%adjoint)) return
            state => ctx%primal; adjoint => ctx%adjoint
            n = state%n; length = 2*n*n; reverse_cell = ctx%reverse_cell
            if (size(values) /= length .or. size(derivative) /= length) return
            call qe_fraction(t, state%h, s)
            call qe_complement(s, forward_s)
            forward_cell = state%steps - reverse_cell + 1
            do i = 1, n
                call qe_evaluate(adjoint%c(i, 0, reverse_cell), &
                         adjoint%c(i, 1, reverse_cell), adjoint%c(i, 2, reverse_cell), &
                                 adjoint%c(i, 3, reverse_cell), s, lr)
                call qe_evaluate(adjoint%c(n + i, 0, reverse_cell), &
                                 adjoint%c(n + i, 1, reverse_cell), &
                                 adjoint%c(n + i, 2, reverse_cell), &
                                 adjoint%c(n + i, 3, reverse_cell), s, li)
                do j = 1, n
                    call qe_evaluate(state%c(j, 0, forward_cell), &
                             state%c(j, 1, forward_cell), state%c(j, 2, forward_cell), &
                                     state%c(j, 3, forward_cell), forward_s, pr)
                    call qe_evaluate(state%c(n + j, 0, forward_cell), &
                                     state%c(n + j, 1, forward_cell), &
                                     state%c(n + j, 2, forward_cell), &
                                     state%c(n + j, 3, forward_cell), forward_s, pi)
                    call qe_rhs_term_vjp(state%hr(i, j), state%hi(i, j), &
                                         pr, pi, state%hbar, lr, li, gr, gi, gx, gy, gh)
                    idx = i + (j - 1)*n
                    derivative(idx) = gr
                    derivative(n*n + idx) = gi
                end do
            end do
        end select
    end subroutine
end module kin6d_quantum_evolution
