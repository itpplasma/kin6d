module kin6d_equilibrium_readback
    ! One readback for every lane (docs/equilibrium-core.md, D7, migration
    ! M2). An axisymmetric nested state is read through its evaluator and its
    ! label primitives only:
    !   axis      O-point of the flux, then the magnetic axis (B_pol = 0);
    !   levels    psi_j = psi_b + (psi_a - psi_b)(1 - rho_j^2), rho uniform,
    !             with psi_b the plasma-boundary flux of the state (no psi = 0
    !             boundary assumption);
    !   q         from B: q = sigma oint B_phi dl/|grad psi| /(2 pi), with
    !             R B_phi taken from the evaluator at a point of the level
    !             (a surface quantity of an axisymmetric equilibrium) and the
    !             lane's coarea oint dl/(R |grad psi|); sigma = sign(psi_b -
    !             psi_a) gives the signed COCOS 3 q of the psi lane;
    !   Phi_tor   int (dPhi/drho) drho with dPhi/dpsi = 2 pi q, normalised to
    !             s_tor; volume and edge toroidal flux from region integrals
    !             of B_phi when the lane provides them, else the contour
    !             integral (volume and flux_consistency NaN).
    ! The flux values reported are absolute in the state's gauge.
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
    use kin6d_equilibrium_state, only: equilibrium_state_t, equilibrium_field_t, &
        labelled_field_t, EQ_BPHI, EQ_PSI, EQ_PSI_R, EQ_VALUES
    use fortnum_polynomial, only: lagrange_weights
    use fortnum_quadrature, only: gauss_legendre_ab
    implicit none
    private
    public :: equilibrium_readback_t, read_equilibrium, readback_q

    type :: equilibrium_readback_t
        real(dp) :: axis(2) = 0, psi_axis = 0, psi_bound = 0, volume = 0
        real(dp) :: phi_edge = 0, flux_consistency = 0, contour_s = 0, integral_s = 0
        integer, allocatable :: contour_work(:)
        real(dp), allocatable :: psi(:), s_tor(:), q(:)
    end type
contains
    real(dp) function readback_q(result, s) result(q)
        ! Cubic Lagrange interpolation of q in the normalised toroidal flux.
        type(equilibrium_readback_t), intent(in) :: result
        real(dp), intent(in) :: s
        real(dp) :: weights(4)
        integer :: i, first
        q = ieee_value(0.0_dp, ieee_quiet_nan)
        if (.not. ieee_is_finite(s)) return
        if (s < 0 .or. s > 1) return
        i = min(size(result%q) - 1, max(1, count(result%s_tor < s)))
        first = max(1, min(i - 1, size(result%q) - 3))
        call lagrange_weights(4, s, result%s_tor(first:first + 3), weights)
        q = dot_product(weights, result%q(first:first + 3))
    end function

    subroutine read_equilibrium(eq, radial_points, angular_points, result, status, &
                                rho_grid)
        !! rho_grid, when present, replaces the uniform rho levels: it must
        !! increase strictly from 0 to 1 (5..4097 points); radial_points is
        !! then ignored.
        type(equilibrium_state_t), intent(in) :: eq
        integer, intent(in) :: radial_points, angular_points
        type(equilibrium_readback_t), intent(out) :: result
        integer, intent(out) :: status
        real(dp), intent(in), optional :: rho_grid(:)
        status = 1
        if (present(rho_grid)) then
            if (size(rho_grid) < 5 .or. size(rho_grid) > 4097) return
            if (.not. all(ieee_is_finite(rho_grid))) return
            if (abs(rho_grid(1)) > 0) return
            if (abs(rho_grid(size(rho_grid)) - 1) > 0) return
            if (any(rho_grid(2:) <= rho_grid(:size(rho_grid) - 1))) return
        end if
        if (.not. allocated(eq%field)) return
        select type (field => eq%field)
        class is (labelled_field_t)
            call read_labelled(eq, field, radial_points, angular_points, result, &
                status, rho_grid)
        end select
    end subroutine

    subroutine read_labelled(eq, field, radial_points, angular_points, result, &
                             status, rho_grid)
        type(equilibrium_state_t), intent(in) :: eq
        class(labelled_field_t), intent(in) :: field
        integer, intent(in) :: radial_points, angular_points
        type(equilibrium_readback_t), intent(inout) :: result
        integer, intent(out) :: status
        real(dp), intent(in), optional :: rho_grid(:)
        real(dp), allocatable :: coarea(:), density(:), rho(:), f(:)
        real(dp) :: two_pi, integral, gx(2), gw(2), weights(4), r, sigma, depth
        real(dp) :: volume, phi_tor
        integer(int64) :: tick, tock, clock_rate
        integer :: i, nr, na, first, k, region_status
        status = 1
        nr = radial_points; na = angular_points
        call field%resolution(nr, na)
        if (present(rho_grid)) nr = size(rho_grid)
        if (nr < 5 .or. nr > 4097) return
        result%psi_bound = eq%psi_bound
        call field%label_axis(result%axis, result%psi_axis, status)
        if (status /= 0) return
        status = 1
        depth = result%psi_axis - eq%psi_bound
        if (abs(depth) <= tiny(1.0_dp) .or. .not. ieee_is_finite(depth)) return
        sigma = sign(1.0_dp, -depth)
        two_pi = 2*acos(-1.0_dp)
        allocate (result%psi(nr), result%s_tor(nr), result%q(nr), &
            coarea(nr), density(nr), rho(nr), f(nr), result%contour_work(nr))
        do i = 1, nr
            rho(i) = real(i - 1, dp)/(nr - 1)
            if (present(rho_grid)) rho(i) = rho_grid(i)
            result%psi(i) = eq%psi_bound + depth*(1 - rho(i)**2)
        end do
        call level_rbphi(field, result%axis, result%psi, f, status)
        if (status /= 0) return
        status = 1
        if (abs(f(1)) <= tiny(1.0_dp)) return
        call system_clock(tick, clock_rate)
        call field%label_coarea(result%psi, na, coarea, status, &
            work=result%contour_work)
        call system_clock(tock)
        result%contour_s = real(tock - tick, dp)/clock_rate
        if (status /= 0) return
        density = f*coarea
        result%q = sigma*density/two_pi
        result%s_tor(1) = 0
        ! Integrate the cubic interpolant of dPhi/drho. Working in sqrt(s_pol)
        ! resolves the axis cell without giving its Hessian a full O(h) weight.
        density = -2*(sigma*depth)*rho*density
        call gauss_legendre_ab(2, 0.0_dp, 1.0_dp, gx, gw)
        do i = 2, nr
            result%s_tor(i) = result%s_tor(i - 1)
            first = max(1, min(i - 2, nr - 3))
            do k = 1, 2
                r = rho(i - 1) + (rho(i) - rho(i - 1))*gx(k)
                call lagrange_weights(4, r, rho(first:first + 3), weights)
                result%s_tor(i) = result%s_tor(i) + (rho(i) - rho(i - 1)) &
                    *gw(k)*dot_product(weights, density(first:first + 3))
            end do
        end do
        integral = result%s_tor(nr)
        if (abs(integral) <= tiny(1.0_dp)) return
        result%s_tor = result%s_tor/integral
        if (.not. all(ieee_is_finite(result%s_tor))) return
        if (any(result%s_tor(2:) <= result%s_tor(:nr - 1))) return
        call system_clock(tick)
        call field%region_integrals(volume, phi_tor, region_status)
        call system_clock(tock)
        result%integral_s = real(tock - tick, dp)/clock_rate
        if (region_status == 0) then
            result%volume = volume
            result%phi_edge = phi_tor
            result%flux_consistency = abs(integral/result%phi_edge - 1)
        else if (region_status < 0) then
            result%phi_edge = integral
            result%volume = ieee_value(0.0_dp, ieee_quiet_nan)
            result%flux_consistency = ieee_value(0.0_dp, ieee_quiet_nan)
        else
            status = region_status
            return
        end if
        call field%magnetic_axis(result%axis, status)
    end subroutine

    subroutine level_rbphi(field, axis, levels, rbphi, status)
        ! R B_phi at one point of every level: the label O-point for the axis
        ! level, else the crossing of the outboard ray from the axis, located
        ! by bracketed Newton steps on the flux (batched evaluator calls). A
        ! level not reached inside the evaluator's domain (the boundary level
        ! of a curved mesh) takes the cubic extrapolation in psi of R B_phi
        ! from four ray points just inside the domain edge.
        class(equilibrium_field_t), intent(in) :: field
        real(dp), intent(in) :: axis(2), levels(:)
        real(dp), intent(out) :: rbphi(:)
        integer, intent(out) :: status
        integer, parameter :: table = 512, max_iterations = 200
        real(dp), allocatable :: points(:, :), values(:, :), t(:), lo(:), hi(:)
        real(dp), allocatable :: x(:)
        integer, allocatable :: statuses(:)
        logical, allocatable :: active(:), beyond(:)
        logical :: newton
        real(dp) :: inside, outside, step, target, tolerance, d, edge_psi(4)
        real(dp) :: edge_rbphi(4), weights(4)
        integer :: n, i, j, iteration
        status = 1
        n = size(levels)
        tolerance = 8*epsilon(1.0_dp)*max(maxval(abs(levels)), &
            abs(levels(n) - levels(1)))
        allocate (points(2, max(table + 1, n)), values(EQ_VALUES, max(table + 1, n)), &
            statuses(max(table + 1, n)))
        points(:, 1) = axis
        call field%sample(points(:, :1), values(:, :1), statuses(:1))
        if (statuses(1) /= 0) return
        rbphi(1) = axis(1)*values(EQ_BPHI, 1)
        ! Domain edge on the outboard ray: doubling, then bisection.
        inside = 0
        outside = 1.0e-3_dp*max(1.0_dp, abs(axis(1)))
        do i = 1, 64
            points(:, 1) = axis + [outside, 0.0_dp]
            call field%sample(points(:, :1), values(:, :1), statuses(:1))
            if (statuses(1) /= 0) exit
            inside = outside
            outside = 2*outside
        end do
        if (statuses(1) == 0) return
        do i = 1, 64
            points(:, 1) = axis + [(inside + outside)/2, 0.0_dp]
            call field%sample(points(:, :1), values(:, :1), statuses(:1))
            if (statuses(1) == 0) then
                inside = (inside + outside)/2
            else
                outside = (inside + outside)/2
            end if
        end do
        if (inside <= 0) return
        ! Normalised flux x = (psi - psi_a)/(psi_b - psi_a) on a ray table.
        allocate (t(0:table), x(0:table))
        t = [(inside*i/table, i=0, table)]
        do i = 0, table
            points(:, i + 1) = axis + [t(i), 0.0_dp]
        end do
        call field%sample(points(:, :table + 1), values(:, :table + 1), &
            statuses(:table + 1))
        if (any(statuses(:table + 1) /= 0)) return
        x = (values(EQ_PSI, :table + 1) - levels(1))/(levels(n) - levels(1))
        allocate (lo(n), hi(n), active(n), beyond(n))
        active = .false.; beyond = .false.
        do j = 2, n
            target = (levels(j) - levels(1))/(levels(n) - levels(1))
            i = 1
            do while (i <= table)
                if (x(i) >= target) exit
                i = i + 1
            end do
            if (i <= table) then
                lo(j) = t(i - 1); hi(j) = t(i)
                active(j) = .true.
                points(:, j) = axis + [(lo(j) + hi(j))/2, 0.0_dp]
            else
                beyond(j) = .true.
                points(:, j) = axis
            end if
        end do
        ! Bracketed Newton; a step leaving the bracket falls back to bisection.
        do iteration = 1, max_iterations
            if (.not. any(active)) exit
            call field%sample(points(:, :n), values(:, :n), statuses(:n))
            do j = 2, n
                if (.not. active(j)) cycle
                if (statuses(j) /= 0) return
                d = points(1, j) - axis(1)
                rbphi(j) = points(1, j)*values(EQ_BPHI, j)
                target = values(EQ_PSI, j) - levels(j)
                if (abs(target) <= tolerance .or. hi(j) - lo(j) <= &
                    4*epsilon(1.0_dp)*(abs(axis(1)) + inside)) then
                    active(j) = .false.
                    cycle
                end if
                if ((target < 0) .eqv. (levels(n) > levels(1))) then
                    lo(j) = d
                else
                    hi(j) = d
                end if
                newton = abs(values(EQ_PSI_R, j)) > tiny(1.0_dp)
                if (newton) then
                    step = -target/values(EQ_PSI_R, j)
                    newton = d + step > lo(j) .and. d + step < hi(j)
                end if
                if (newton) then
                    ! A Newton step at the coordinate roundoff has converged.
                    if (abs(step) <= 4*epsilon(1.0_dp)*(abs(axis(1)) + inside)) then
                        active(j) = .false.
                        cycle
                    end if
                    points(1, j) = axis(1) + d + step
                else
                    points(1, j) = axis(1) + (lo(j) + hi(j))/2
                end if
            end do
        end do
        if (any(active)) return
        if (any(beyond)) then
            do i = 1, 4
                points(:, i) = axis + [inside*(1 - (i - 1)*1.0e-4_dp), 0.0_dp]
            end do
            call field%sample(points(:, :4), values(:, :4), statuses(:4))
            if (any(statuses(:4) /= 0)) return
            edge_psi = values(EQ_PSI, :4)
            edge_rbphi = points(1, :4)*values(EQ_BPHI, :4)
            do j = 2, n
                if (.not. beyond(j)) cycle
                call lagrange_weights(4, levels(j), edge_psi, weights)
                rbphi(j) = dot_product(weights, edge_rbphi)
            end do
        end if
        if (.not. all(ieee_is_finite(rbphi))) return
        status = 0
    end subroutine
end module kin6d_equilibrium_readback
