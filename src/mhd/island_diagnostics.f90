module kin6d_island_diagnostics
    ! Poincare diagnostics for any equilibrium_field_t. Angles are radians,
    ! section coordinates (R,Z) are metres. Estimates are tolerance/refinement
    ! indicators, never continuous enclosures. No global field or workspace.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use kin6d_equilibrium_state, only: equilibrium_field_t, EQ_VALUES, &
                                       EQ_BR, EQ_BZ, EQ_BPHI
    use kin6d_nonnested_leaves, only: fieldline_velocity
    use kin6d_island_leaves, only: poincare_invariants, helical_phase
    use fortnum_ode, only: ode_problem_t, ode_workspace_t, ode_solution_t, &
                           ode_integrate
    use fortnum_ode_cash_karp, only: cash_karp_step
    use fortnum_multiroot, only: multiroot_hybrid
    use fortnum_status, only: fortnum_status_t
    implicit none
    private
    public :: island_trace_options_t, island_point_t, island_width_t
    public :: island_trace, island_periodic_point, island_stability
    public :: island_symmetry_point, island_libration_edge, island_width
    real(dp), parameter :: pi = acos(-1.0_dp)

    type :: island_trace_options_t
        real(dp) :: rtol = 1.0e-8_dp, atol = 1.0e-10_dp
        real(dp) :: max_step = 2*pi/128
        real(dp) :: map_delta = 1.0e-3_dp, root_tol = 1.0e-7_dp
        integer :: max_steps = 200000
        ! Nonzero selects the legacy fixed Cash-Karp schedule, per turn.
        ! This is useful for explicit step ladders, not tolerance control.
        integer :: fixed_steps = 0
    end type

    type :: island_point_t
        real(dp) :: x(2) = 0, closure = huge(1.0_dp)
        real(dp) :: residue = 0, determinant = 0, phase = 0
        real(dp) :: position_error = huge(1.0_dp), residue_error = huge(1.0_dp)
        real(dp) :: phase_error = huge(1.0_dp)
        logical :: used_symmetry = .false.
        ! 0: reference comparison; 1: input/domain; 2: trace; 3: unresolved root.
        integer :: status = 1
    end type

    type :: island_width_t
        real(dp) :: edge(2) = 0, metres = 0, label_width = 0
        real(dp) :: trace_error = 0, duration_error = 0, scan_error = 0
        character(80) :: radial_label = ''
        character(48) :: claim = 'finite-time-libration; enclosure unavailable'
        integer :: status = 1
    end type

    abstract interface
        function radial_label_t(x) result(label)
            import dp
            real(dp), intent(in) :: x(2)
            real(dp) :: label
        end function
    end interface
contains
    subroutine island_trace(field, start, phi0, turns, opt, finish, status, &
                            centre, m, n, span)
        class(equilibrium_field_t), intent(in) :: field
        real(dp), intent(in) :: start(2), phi0, turns
        type(island_trace_options_t), intent(in) :: opt
        real(dp), intent(out) :: finish(2)
        integer, intent(out) :: status
        real(dp), intent(in), optional :: centre(2)
        integer, intent(in), optional :: m, n
        real(dp), intent(out), optional :: span
        type(ode_problem_t) :: problem
        type(ode_workspace_t) :: work
        type(ode_solution_t) :: sol
        type(fortnum_status_t) :: stat
        real(dp) :: k1(2), k2(2), k3(2), k4(2), k5(2), k6(2), tmp(2)
        real(dp) :: next(2), err(2), h, phi, angle, oldtheta, theta, lo, hi
        integer :: j, count, nfev, fault
        finish = start; status = 1; fault = 0
        if (present(span)) span = 0
        if (abs(turns) <= tiny(1.0_dp) .or. opt%max_step <= 0 .or. opt%rtol <= 0 &
            .or. opt%atol <= 0 .or. opt%fixed_steps < 0) return
        if (present(span)) then
            if (.not. present(centre)) return
            if (.not. present(m)) return
            if (.not. present(n)) return
            if (norm2(start - centre) <= tiny(1.0_dp)) return
            oldtheta = atan2(start(2) - centre(2), start(1) - centre(1))
            angle = oldtheta
            lo = m*angle - n*phi0; hi = lo
        end if
        if (opt%fixed_steps > 0) then
            count = max(1, nint(abs(turns)*opt%fixed_steps))
            h = 2*pi*turns/count; nfev = 0
            do j = 1, count
                phi = phi0 + h*(j - 1)
                call cash_karp_step(rhs, phi, finish, h, .false., k1, k2, &
                                    k3, k4, k5, k6, tmp, next, err, nfev)
                if (fault /= 0) exit
                finish = next
                if (present(span)) then
                    call accumulate(phi + h, finish)
                    if (abs(turns) > 1 .and. hi - lo > 2*pi) exit
                end if
            end do
        else
            problem%rhs => rhs
            problem%t0 = phi0; problem%t1 = phi0 + 2*pi*turns
            problem%y0 = start; problem%rtol = opt%rtol; problem%atol = opt%atol
            problem%hmax = opt%max_step; problem%max_steps = opt%max_steps
            call ode_integrate(problem, work, sol, stat)
            if (stat%code /= 0) fault = 2
            if (fault == 0) then
                finish = sol%y(:, size(sol%t))
                if (present(span)) then
                    do j = 2, size(sol%t)
                        call accumulate(sol%t(j), sol%y(:, j))
                    end do
                end if
            end if
        end if
        status = 0
        if (fault /= 0 .or. .not. all(ieee_is_finite(finish))) status = 2
        if (present(span)) span = hi - lo
    contains
        subroutine rhs(t, y, dy, ctx)
            real(dp), intent(in) :: t, y(:)
            real(dp), intent(out) :: dy(:)
            class(*), intent(in), optional :: ctx
            real(dp) :: point(3, 1), value(EQ_VALUES, 1)
            integer :: st(1)
            dy = 0
            if (fault /= 0) return
            point(:, 1) = [y, t]
            call field%sample(point, value, st)
            if (st(1) /= 0 .or. y(1) <= 0) then
                fault = 1
                return
            end if
            if (.not. all(ieee_is_finite(value(:3, 1)))) then
                fault = 1
                return
            end if
            if (abs(value(EQ_BPHI, 1)) <= tiny(1.0_dp)) then
                fault = 1
                return
            end if
            call fieldline_velocity(y(1), value(EQ_BR, 1), value(EQ_BZ, 1), &
                                    value(EQ_BPHI, 1), 0.0_dp, dy(1), dy(2))
        end subroutine
        subroutine accumulate(t, y)
            real(dp), intent(in) :: t, y(2)
            theta = atan2(y(2) - centre(2), y(1) - centre(1))
            angle = angle + atan2(sin(theta - oldtheta), cos(theta - oldtheta))
            oldtheta = theta
            lo = min(lo, m*angle - n*t); hi = max(hi, m*angle - n*t)
        end subroutine
    end subroutine

    subroutine island_stability(field, root, phi0, turns, opt, jac, closure, status)
        class(equilibrium_field_t), intent(in) :: field
        real(dp), intent(in) :: root(2), phi0, turns
        type(island_trace_options_t), intent(in) :: opt
        real(dp), intent(out) :: jac(2, 2), closure
        integer, intent(out) :: status
        real(dp) :: finish(2), plus(2), minus(2), x(2)
        integer :: j, statuses(2)
        jac = 0; closure = huge(1.0_dp)
        status = 1
        if (opt%map_delta <= 0) return
        call island_trace(field, root, phi0, turns, opt, finish, status)
        if (status /= 0) return
        closure = norm2(finish - root)
        statuses = 0
        !$omp parallel do default(shared) private(j, x, plus, minus) schedule(static)
        do j = 1, 2
            x = root; x(j) = x(j) + opt%map_delta
            call island_trace(field, x, phi0, turns, opt, plus, statuses(j))
            if (statuses(j) /= 0) cycle
            x(j) = root(j) - opt%map_delta
            call island_trace(field, x, phi0, turns, opt, minus, statuses(j))
            if (statuses(j) /= 0) cycle
            jac(:, j) = (plus - minus)/(2*opt%map_delta)
        end do
        !$omp end parallel do
        status = maxval(statuses)
    end subroutine

    subroutine island_symmetry_point(field, bracket, phi0, turns, opt, root, status)
        ! ONLY for a reversible field on Z=0 at phi0: half-period map Z=0.
        ! The complete-map closure must still be checked by the caller.
        class(equilibrium_field_t), intent(in) :: field
        real(dp), intent(in) :: bracket(2), phi0, turns
        type(island_trace_options_t), intent(in) :: opt
        real(dp), intent(out) :: root(2)
        integer, intent(out) :: status
        real(dp) :: a, b, fa, fb, finish(2)
        integer :: i
        a = bracket(1); b = bracket(2); root = [a, 0.0_dp]
        call island_trace(field, root, phi0, turns/2, opt, finish, status)
        if (status /= 0) return
        fa = finish(2)
        call island_trace(field, [b, 0.0_dp], phi0, turns/2, opt, finish, status)
        if (status /= 0) return
        fb = finish(2); status = 3
        if (fa*fb > 0) return
        do i = 1, 50
            root = [(a + b)/2, 0.0_dp]
            call island_trace(field, root, phi0, turns/2, opt, finish, status)
            if (status /= 0) return
            if (abs(b - a) < opt%root_tol) return
            if (fa*finish(2) <= 0) then
                b = root(1)
            else
                a = root(1); fa = finish(2)
            end if
        end do
        status = 3
    end subroutine

    subroutine island_periodic_point(field, seed, phi0, turns, centre, m, n, &
                                     opt, report, symmetry_bracket)
        class(equilibrium_field_t), intent(in) :: field
        real(dp), intent(in) :: seed(2), phi0, centre(2)
        integer, intent(in) :: turns, m, n
        type(island_trace_options_t), intent(in) :: opt
        type(island_point_t), intent(out) :: report
        real(dp), intent(in), optional :: symmetry_bracket(2)
        type(island_trace_options_t) :: tight
        type(fortnum_status_t) :: stat
        real(dp) :: roots(2, 2), jac(2, 2), res(2), closure, residues(2), phases(2)
        integer :: lev, st, fault
        logical :: fallback
        report = island_point_t(); tight = opt
        if (turns < 1 .or. opt%map_delta <= 0 .or. opt%root_tol <= 0) return
        roots(:, 1) = seed; roots(:, 2) = seed
        do lev = 1, 2
            if (lev == 2) then
                tight%rtol = opt%rtol/8; tight%atol = opt%atol/8
                tight%max_step = opt%max_step/2
                tight%map_delta = opt%map_delta/2
                tight%fixed_steps = 2*opt%fixed_steps
                roots(:, 2) = roots(:, 1)
            end if
            fault = 0; fallback = .false.
            call multiroot_hybrid(map_jac, 2, roots(:, lev), report%x, stat, &
                               xtol=tight%root_tol/10, ftol=tight%root_tol, max_iter=30)
            call map_jac(report%x, res, jac)
            if (fault /= 0 .or. norm2(res) > tight%root_tol) then
                if (.not. present(symmetry_bracket)) then
                    report%status = 3
                    return
                end if
                call island_symmetry_point(field, symmetry_bracket, phi0, &
                                           real(turns, dp), tight, report%x, st)
                if (st /= 0) then
                    report%status = st
                    return
                end if
                fallback = .true.
                ! A symmetric seed need not close on a non-reflected FEEC mesh.
                ! Refine that seed with the complete map before accepting it.
                roots(:, lev) = report%x
                fault = 0
                call map_jac(report%x, res, jac)
                if (fault /= 0 .or. norm2(res) > tight%root_tol) then
                    fault = 0
                    call multiroot_hybrid(map_jac, 2, roots(:, lev), report%x, stat, &
                        xtol=tight%root_tol/10, ftol=tight%root_tol, max_iter=30)
                    call map_jac(report%x, res, jac)
                    if (fault /= 0 .or. norm2(res) > tight%root_tol) then
                        report%status = 3
                        return
                    end if
                end if
            end if
            report%used_symmetry = report%used_symmetry .or. fallback
            roots(:, lev) = report%x
            call island_stability(field, report%x, phi0, real(turns, dp), &
                                  tight, jac, closure, st)
            if (st /= 0) then
                report%status = st
                return
            end if
            report%closure = closure
            call poincare_invariants(jac(1, 1), jac(1, 2), jac(2, 1), &
                                     jac(2, 2), residues(lev), report%determinant)
            call helical_phase(report%x(1) - centre(1), &
                   report%x(2) - centre(2), real(m, dp), real(n, dp), phi0, phases(lev))
        end do
        report%residue = residues(2)
        report%phase = atan2(sin(phases(2)), cos(phases(2)))
        report%position_error = norm2(roots(:, 2) - roots(:, 1)) + opt%root_tol
        report%residue_error = abs(residues(2) - residues(1))
        report%phase_error = abs(atan2(sin(phases(2) - phases(1)), &
                                       cos(phases(2) - phases(1))))
        if (norm2(report%x - centre) > report%position_error) then
            report%phase_error = report%phase_error + abs(m)*report%position_error &
                /(norm2(report%x - centre) - report%position_error)
        else
            report%phase_error = huge(1.0_dp)
        end if
        report%status = 0
        if (report%closure > opt%root_tol) report%status = 3
    contains
        subroutine map_jac(x, f, jacobian, ctx)
            real(dp), intent(in) :: x(:)
            real(dp), intent(out) :: f(:), jacobian(:, :)
            class(*), intent(in), optional :: ctx
            real(dp) :: finish(2), defect
            integer :: j, ierr
            call island_trace(field, x, phi0, real(turns, dp), tight, finish, ierr)
            f = finish - x
            call island_stability(field, x, phi0, real(turns, dp), tight, &
                                  jacobian, defect, st)
            if (ierr /= 0 .or. st /= 0) fault = 1
            do j = 1, 2
                jacobian(j, j) = jacobian(j, j) - 1
            end do
        end subroutine
    end subroutine

    subroutine island_libration_edge(field, outside, inside, origin, direction, &
                                     phi0, turns, m, n, opt, edge, status, scans)
        ! Coordinate along a straight section chord, direction a unit vector.
        class(equilibrium_field_t), intent(in) :: field
        real(dp), intent(in) :: outside, inside, origin(2), direction(2)
        real(dp), intent(in) :: phi0, turns
        integer, intent(in) :: m, n
        type(island_trace_options_t), intent(in) :: opt
        real(dp), intent(out) :: edge
        integer, intent(out) :: status
        integer, intent(in), optional :: scans
        real(dp) :: a, b, finish(2), span
        integer :: j, count
        a = outside; b = inside; count = 12
        if (present(scans)) count = scans
        edge = a; status = 1
        if (abs(norm2(direction) - 1) > 1.0e-12_dp .or. count < 1) return
        call island_trace(field, origin + a*direction, phi0, turns, opt, &
                          finish, status, origin, m, n, span)
        if (status /= 0) return
        status = 3
        if (span <= 2*pi) return
        call island_trace(field, origin + b*direction, phi0, turns, opt, &
                          finish, status, origin, m, n, span)
        if (status /= 0) return
        status = 3
        if (span > 2*pi) return
        do j = 1, count
            edge = (a + b)/2
            call island_trace(field, origin + edge*direction, phi0, turns, opt, &
                              finish, status, origin, m, n, span)
            if (status /= 0) return
            if (span > 2*pi) then
                a = edge
            else
                b = edge
            end if
        end do
        edge = (a + b)/2
    end subroutine

    subroutine island_width(field, outside, inside, origin, direction, phi0, &
                            turns, m, n, opt, label, label_name, report, scans)
        class(equilibrium_field_t), intent(in) :: field
        real(dp), intent(in) :: outside(2), inside, origin(2), direction(2)
        real(dp), intent(in) :: phi0, turns
        integer, intent(in) :: m, n
        type(island_trace_options_t), intent(in) :: opt
        procedure(radial_label_t) :: label
        character(*), intent(in) :: label_name
        type(island_width_t), intent(out) :: report
        integer, intent(in), optional :: scans
        type(island_trace_options_t) :: tight
        real(dp) :: edges(2, 3), widths(3)
        integer :: lev, j, statuses(2), count
        report = island_width_t(); report%radial_label = label_name
        count = 12
        if (present(scans)) count = scans
        if (len_trim(label_name) == 0 .or. count < 1 .or. count > 24) return
        tight = opt
        do lev = 1, 3
            if (lev == 2) then
                tight%rtol = opt%rtol/8; tight%atol = opt%atol/8
                tight%max_step = opt%max_step/2
                tight%fixed_steps = 2*opt%fixed_steps
            end if
            !$omp parallel do default(shared) private(j) schedule(static)
            do j = 1, 2
                call island_libration_edge(field, outside(j), inside, origin, &
                            direction, phi0, turns*merge(2, 1, lev == 3), m, n, tight, &
                                           edges(j, lev), statuses(j), scans=count)
            end do
            !$omp end parallel do
            if (any(statuses /= 0)) then
                report%status = maxval(statuses)
                return
            end if
            widths(lev) = abs(edges(2, lev) - edges(1, lev))
        end do
        report%edge = edges(:, 2); report%metres = widths(2)
        report%label_width = abs(label(origin + edges(2, 2)*direction) &
                                 - label(origin + edges(1, 2)*direction))
        report%trace_error = abs(widths(2) - widths(1))
        report%duration_error = abs(widths(3) - widths(2))
        report%scan_error = sum(abs(outside - inside))/2**count
        report%status = 0
    end subroutine
end module
