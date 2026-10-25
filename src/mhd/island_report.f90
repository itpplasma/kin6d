module kin6d_island_report
    ! Explicit diagnostic pass, independent of Newton convergence and acceptance.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, &
        ieee_value, ieee_quiet_nan
    use fortfem_meridian_feec_n, only: meridian_qpoint, meridian_qgrad
    use kin6d_equilibrium_nonnested, only: free_b_field_t, nn_divergence
    use kin6d_nonnested_coupled, only: nc_problem_t, nc_layout_t, nc_layout, &
        nc_residual, nc_store, nc_source_data
    use kin6d_nonnested_residual, only: nn_residual, NN_RESIDUAL_ZERO_SCALE
    use kin6d_equilibrium_residual, only: equilibrium_residual_t
    use kin6d_nonnested_leaves, only: island_defect_density, island_defect_norm
    implicit none
    private
    public :: island_report_t, island_residual_report
    type :: island_report_t
        ! Euclidean norms of free weak rows, with NN2's phi-mean measure.
        real(dp) :: ampere = 0, gauge = 0, pressure = 0, ohm_continuity = 0
        real(dp) :: drive = 0, current = 0, voltage = 0
        ! Full-torus integrated balances. Heat flux is the Dirichlet reaction.
        real(dp) :: heating = 0, outward_heat = 0, heat_balance = 0
        ! Fitted-source reaction is for live minus seed flux, not physical heat.
        real(dp) :: weak_heat_reaction = 0, weak_heat_balance = 0
        character(64) :: heat_status = 'physical scalar-source Dirichlet balance'
        real(dp) :: outward_current = 0, charge_balance = 0
        real(dp) :: epsilon_parallel = 0, perpendicular_ohmic_force_rms = 0
        real(dp) :: algebraic_divergence = 0, length_scale = 0, pressure_scale = 0
        type(equilibrium_residual_t) :: force
        integer :: test_modes = 0, nphi = 0, quadrature_points = 0
        character(64) :: force_status = 'finite-test-space norm; no continuous bound'
        character(32) :: claim = 'candidate'
        character(32) :: interface_status = 'not applicable: single region'
        character(32) :: gamma_c_status = 'unavailable: prescribed trace'
        character(32) :: sampled_divergence_status = 'unavailable'
        character(32) :: certificate_status = 'unavailable: no stability bound'
        character(32) :: derivative_status = 'fixed branch, nonzero norms'
    end type
contains
    subroutine island_residual_report(field, problem, x, length_scale, &
        pressure_scale, report, status, test_modes)
        type(free_b_field_t), intent(inout) :: field
        type(nc_problem_t), intent(in) :: problem
        real(dp), intent(in) :: x(:), length_scale, pressure_scale
        type(island_report_t), intent(out) :: report
        integer, intent(out) :: status
        integer, intent(in), optional :: test_modes
        type(nc_layout_t) :: lay
        real(dp), allocatable :: r(:), raw(:), sums(:, :)
        real(dp) :: total(4), two_pi
        integer :: st, t, nt, failure
        status = 1
        if (.not. all(ieee_is_finite([length_scale, pressure_scale]))) return
        if (length_scale <= 0 .or. pressure_scale <= 0) return
        call nc_residual(field, problem, x, r, status)
        if (status /= 0) return
        call nc_residual(field, problem, x, raw, status, boundary_rows=.true.)
        if (status /= 0) return
        lay = nc_layout(field)
        report%ampere = norm2(r(1:lay%il))
        report%gauge = norm2(r(lay%il + 1:lay%ip))
        report%pressure = norm2(r(lay%ip + 1:lay%iphi))
        report%ohm_continuity = norm2(r(lay%iphi + 1:lay%iv - 1))
        report%current = raw(lay%iv)
        report%voltage = x(lay%iv)
        if (problem%closure%drive == 'V') then
            report%drive = x(lay%iv) - problem%voltage
        else
            report%drive = raw(lay%iv) - problem%current
        end if
        two_pi = 2*acos(-1.0_dp)
        report%outward_heat = -two_pi*sum(raw(lay%ip + 1:lay%ip + lay%n0), &
            mask=field%fs%sp(0)%bnd)
        report%outward_current = two_pi*sum(raw(lay%iphi + 1:lay%iphi + lay%n0), &
            mask=field%fs%sp(0)%bnd)
        report%charge_balance = abs(report%outward_current)
        report%length_scale = length_scale; report%pressure_scale = pressure_scale
        report%nphi = 4*lay%modes + 1
        if (problem%nphi > 0) report%nphi = problem%nphi
        report%quadrature_points = field%fs%nq
        report%test_modes = lay%modes
        if (present(test_modes)) report%test_modes = test_modes
        call nc_store(field, x, status)
        if (status /= 0) return
        nt = field%fs%mm%topology%n_triangles
        allocate (sums(4, nt))
        failure = 0
        !$omp parallel do default(shared) private(t, st) schedule(static) &
        !$omp& reduction(max:failure)
        do t = 1, nt
            call cell(t, sums(:, t), st)
            failure = max(failure, st)
        end do
        !$omp end parallel do
        status = failure
        if (status /= 0) return
        total = sum(sums, dim=2)
        status = 1
        if (total(3) <= 0) return
        if (total(1) <= 0 .or. total(2) <= 0) &
            report%derivative_status = 'unavailable: norm at zero'
        report%heating = two_pi*total(4)
        report%heat_balance = abs(report%outward_heat - report%heating)
        if (allocated(problem%fitted)) then
            report%weak_heat_reaction = report%outward_heat
            report%weak_heat_balance = abs(report%weak_heat_reaction)
            report%heat_status = 'unavailable: weak fit has no physical heat trace'
            report%heating = ieee_value(0.0_dp, ieee_quiet_nan)
            report%outward_heat = ieee_value(0.0_dp, ieee_quiet_nan)
            report%heat_balance = ieee_value(0.0_dp, ieee_quiet_nan)
        end if
        call island_defect_norm(max(0.0_dp, total(1)), max(0.0_dp, total(2)), &
            total(3), length_scale, pressure_scale, report%epsilon_parallel, &
            report%perpendicular_ohmic_force_rms)
        call nn_divergence(field, report%algebraic_divergence, status)
        if (status /= 0) return
        call nn_residual(field, report%test_modes, report%force, status)
        if (status == NN_RESIDUAL_ZERO_SCALE) then
            report%force_status = 'normalization unavailable: zero D3 scale'
            status = 0
        end if
        ! Other D3 failures propagate; its finite norm is not a majorant.
    contains
        subroutine cell(t, acc, status)
            integer, intent(in) :: t
            real(dp), intent(out) :: acc(4)
            integer, intent(out) :: status
            real(dp) :: v2(3, field%fs%sp(2)%nloc), val(field%fs%nl0)
            real(dp) :: grad(2, field%fs%nl0), xq(2), wq, phi, w
            real(dp) :: bm(3, lay%nb), pl(field%fs%nl0, lay%nb)
            real(dp) :: fl(field%fs%nl0, lay%nb), tr(lay%nb), dtr(lay%nb)
            real(dp) :: b(3), gp(3), gf(3), par2, force2, source
            real(dp) :: source_flux(3), drive
            integer :: q, j, n, bk
            acc = 0
            do bk = 1, lay%nb
                pl(:, bk) = x(lay%ip + (bk - 1)*lay%n0 + field%fs%sp(0)%dof(:, t))
                fl(:, bk) = x(lay%iphi + (bk - 1)*lay%n0 &
                    + field%fs%sp(0)%dof(:, t))
            end do
            do q = 1, field%fs%nq
                call meridian_qpoint(field%fs, 2, t, q, xq, wq, v2, status)
                if (status /= 0) return
                bm(:, 1) = real(matmul(v2, field%b(field%fs%sp(2)%dof(:, t), 0) &
                    *field%fs%sp(2)%sgn(:, t)), dp)
                do n = 1, lay%modes
                    bm(:, 2*n) = 2*real(matmul(v2, &
                        field%b(field%fs%sp(2)%dof(:, t), n) &
                        *field%fs%sp(2)%sgn(:, t)), dp)
                    bm(:, 2*n + 1) = -2*aimag(matmul(v2, &
                        field%b(field%fs%sp(2)%dof(:, t), n) &
                        *field%fs%sp(2)%sgn(:, t)))
                end do
                call meridian_qgrad(field%fs, t, q, xq, wq, val, grad, status)
                if (status /= 0) return
                w = wq/report%nphi
                do j = 1, report%nphi
                    phi = two_pi*(j - 1)/report%nphi
                    tr(1) = 1; dtr(1) = 0
                    do n = 1, lay%modes
                        tr(2*n) = cos(n*phi); tr(2*n + 1) = sin(n*phi)
                        dtr(2*n) = -n*sin(n*phi); dtr(2*n + 1) = n*cos(n*phi)
                    end do
                    b = matmul(bm, tr); b(2) = b(2) + field%f_vac/xq(1)
                    gp = [dot_product(matmul(grad(1, :), pl), tr), &
                        dot_product(matmul(val, pl), dtr)/xq(1), &
                        dot_product(matmul(grad(2, :), pl), tr)]
                    gf = [dot_product(matmul(grad(1, :), fl), tr), &
                        dot_product(matmul(val, fl), dtr)/xq(1), &
                        dot_product(matmul(grad(2, :), fl), tr)]
                    call island_defect_density(b(1), b(2), b(3), gp(1), gp(2), &
                        gp(3), gf(1), gf(2), gf(3), problem%g_perp, par2, force2)
                    call nc_source_data(problem, xq, phi, source, &
                        source_flux, drive, status)
                    if (status /= 0) return
                    acc = acc + w*[par2, force2, 1.0_dp, source]
                end do
            end do
            status = 1
            if (all(ieee_is_finite(acc))) status = 0
        end subroutine
    end subroutine
end module
