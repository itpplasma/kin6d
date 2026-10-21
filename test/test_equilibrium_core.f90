module test_equilibrium_core_fields
    ! Solov'ev: psi = c ((R^2 - R0^2)^2 + 4 R^2 Z^2) - 4 c R0^2 a^2 vanishes on
    ! the boundary R^2 = R0^2 + 2 R0 a cos t, R Z = R0 a sin t; mu0 = 1,
    ! mu0 p' = -16 c, F = f.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_equilibrium_state, only: equilibrium_field_t, EQ_BR, EQ_BZ, &
        EQ_BPHI, EQ_P, EQ_PSI, EQ_PSI_R, EQ_PSI_Z
    implicit none
    private
    public :: solovev_field_t
    real(dp), parameter, public :: r0 = 6.2_dp, minor = 0.62_dp, f_tor = 32.86_dp
    real(dp), parameter, public :: c_sol = f_tor/(12*r0**3)
    type, extends(equilibrium_field_t) :: solovev_field_t
    contains
        procedure :: sample => solovev_sample
    end type
contains
    subroutine solovev_sample(field, points, values, statuses)
        class(solovev_field_t), intent(in) :: field
        real(dp), intent(in) :: points(:, :)
        real(dp), intent(out) :: values(:, :)
        integer, intent(out) :: statuses(:)
        real(dp) :: r, z, psi_r, psi_z
        integer :: i
        values = 0; statuses = 0
        do i = 1, size(points, 2)
            r = points(1, i); z = points(2, i)
            psi_r = c_sol*(4*r*(r*r - r0*r0) + 8*r*z*z)
            psi_z = c_sol*8*r*r*z
            values(EQ_PSI, i) = c_sol*((r*r - r0*r0)**2 + 4*r*r*z*z &
                - 4*r0*r0*minor*minor)
            values(EQ_PSI_R, i) = psi_r; values(EQ_PSI_Z, i) = psi_z
            values(EQ_BR, i) = -psi_z/r; values(EQ_BZ, i) = psi_r/r
            values(EQ_BPHI, i) = f_tor/r
            values(EQ_P, i) = -16*c_sol*values(EQ_PSI, i)
        end do
    end subroutine
end module test_equilibrium_core_fields

program test_equilibrium_core
    ! Common equilibrium core (docs/equilibrium-core.md).
    !   chart              jet of the cylindrical chart against finite differences
    !   compare WORK FIX   K6-mini regression: gs_phase1 outputs in WORK against
    !                      the fixture FIX to 1e-12 relative
    !   readback K         D7 on P_K Solov'ev: flux gauge shift and sign flip
    !                      leave q and s_tor unchanged; a plasma sub-mesh
    !                      reproduces the whole-mesh q at its boundary level
    !   residual K         D3 of the exact Solov'ev state (quadrature level) and
    !                      of psi-lane P_K solutions on three meshes (order >= K;
    !                      observed K + 1)
    !   residual_bound     D3 equilibrated upper bound (P2): the exact Solov'ev
    !                      state has D3 = 0, so the bound is pure estimator
    !                      error and shrinks; for a coarse psi-lane state the
    !                      bound brackets D3: the lower bound of a finer
    !                      test space does not exceed it
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_equilibrium_chart, only: cylindrical_chart_t
    use kin6d_equilibrium_state, only: equilibrium_state_t
    use test_equilibrium_core_fields, only: solovev_field_t, r0, minor, f_tor, c_sol
    use kin6d_equilibrium_residual, only: equilibrium_residual_t, &
        equilibrium_residual
    use kin6d_grad_shafranov, only: gs_state_t, initialize_gs, solve_gs
    use kin6d_gs_equilibrium, only: gs_equilibrium
    use kin6d_equilibrium_readback, only: equilibrium_readback_t, read_equilibrium
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan, ieee_is_finite
    implicit none
    character(512) :: mode, work, fixtures
    call get_command_argument(1, mode)
    select case (trim(mode))
    case ('chart')
        call test_chart()
    case ('compare')
        call get_command_argument(2, work)
        call get_command_argument(3, fixtures)
        call compare_k6_mini(trim(work), trim(fixtures))
    case ('synthetic_analytic')
        call get_command_argument(2, work)
        call get_command_argument(3, fixtures)
        call check_synthetic_analytic(trim(work), trim(fixtures))
    case ('readback')
        call get_command_argument(2, work)
        call test_readback(merge(3, 2, trim(work) == '3'))
    case ('residual')
        call get_command_argument(2, work)
        call test_residual(merge(3, 2, trim(work) == '3'))
    case ('residual_bound')
        call test_residual_bound()
    case default
        error stop 'mode: chart | compare WORK FIXTURES | readback K | residual K'
    end select
contains
    subroutine test_residual(order)
        ! Gate of migration M1: D3 of the exact state is at quadrature level;
        ! for psi-lane solutions it decreases at least at order K.
        integer, intent(in) :: order
        type(equilibrium_state_t) :: exact, lane
        type(equilibrium_residual_t) :: report
        type(gs_state_t) :: mesh
        real(dp) :: normalized(3), rate
        integer :: level, status, n
        do level = 1, 3
            n = 12*2**level
            call solovev_mesh(n, order, mesh)
            if (level == 3) then
                allocate (solovev_field_t :: exact%field)
                exact%mu0 = 1; exact%f_vac = f_tor
                call equilibrium_residual(exact, mesh, report, status)
                if (status /= 0) error stop 'exact residual evaluation'
                print '(a,i0,a,es10.3,a,es10.3)', 'exact Solov''ev, n=', n, &
                    ': normalized D3 ', report%normalized, ' sigma ', report%sigma
                if (report%normalized > 1.0e-9_dp) &
                    error stop 'D3 of the exact state is not at quadrature level'
            end if
            call solve_gs(mesh, 1.0_dp, -16*c_sol, 0.0_dp, status)
            if (status /= 0) error stop 'Solov''ev solve'
            call gs_equilibrium(mesh, lane, status, f_tor, 0.0_dp, &
                pprime=-16*c_sol, mu0=1.0_dp)
            if (status /= 0) error stop 'psi adapter'
            call equilibrium_residual(lane, mesh, report, status)
            if (status /= 0) error stop 'psi-lane residual evaluation'
            normalized(level) = report%normalized
            print '(a,i0,a,i0,a,es10.3)', 'psi lane P', order, ', n=', n, &
                ': normalized D3 ', normalized(level)
        end do
        rate = log(normalized(2)/normalized(3))/log(2.0_dp)
        print '(a,2f7.3)', 'observed orders ', &
            log(normalized(1)/normalized(2))/log(2.0_dp), rate
        if (rate < order) error stop 'psi-lane D3 order below K'
        ! TODO(N1): the same gate for the map lane once ideal_map is on main.
        print '(a)', 'PASS: D3 residual'
    end subroutine

    subroutine test_residual_bound()
        ! D3 = 0 for the exact state: dual_upper = sqrt(0 + M^2) is estimator
        ! error only (practical claim, quadrature not enclosed) and must shrink
        ! under refinement. For the coarse psi-lane state, the lower bound on a
        ! finer test space (same state, finer mesh) is below the coarse upper
        ! bound, which is above the coarse lower bound.
        type(equilibrium_state_t) :: exact, lane
        type(equilibrium_residual_t) :: report, coarse, fine
        type(gs_state_t) :: mesh, fine_mesh
        real(dp) :: upper(2)
        integer :: level, status
        allocate (solovev_field_t :: exact%field)
        exact%mu0 = 1; exact%f_vac = f_tor
        do level = 1, 2
            call solovev_mesh(12*2**level, 2, mesh)
            call equilibrium_residual(exact, mesh, report, status, &
                error_estimate=.true.)
            if (status /= 0) error stop 'exact residual evaluation'
            if (report%bound_status /= 0) error stop 'exact upper bound unavailable'
            if (trim(report%claim) /= 'practical-upper-bound') &
                error stop 'upper-bound claim label'
            if (trim(report%quadrature_status) /= 'approximate') &
                error stop 'upper-bound quadrature status'
            if (report%dual_upper < report%dual) error stop 'upper below lower'
            upper(level) = report%normalized_upper
            print '(a,i0,a,es10.3,a,es10.3)', 'exact Solov''ev n=', 12*2**level, &
                ': normalized D3 lower ', report%normalized, ' upper ', upper(level)
        end do
        if (upper(2) > upper(1)/2.8_dp) error stop 'exact upper bound not shrinking'
        call solovev_mesh(24, 2, mesh)
        call solve_gs(mesh, 1.0_dp, -16*c_sol, 0.0_dp, status)
        if (status /= 0) error stop 'Solov''ev solve'
        call gs_equilibrium(mesh, lane, status, f_tor, 0.0_dp, pprime=-16*c_sol, &
            mu0=1.0_dp)
        if (status /= 0) error stop 'psi adapter'
        call equilibrium_residual(lane, mesh, coarse, status, error_estimate=.true.)
        if (status /= 0 .or. coarse%bound_status /= 0) error stop 'coarse bound'
        call solovev_mesh(72, 2, fine_mesh)
        call equilibrium_residual(lane, fine_mesh, fine, status)
        if (status /= 0) error stop 'fine lower bound'
        print '(a,3es11.3)', 'psi lane P2 n=24 dual: lower, upper; finer lower ', &
            coarse%dual, coarse%dual_upper, fine%dual
        if (coarse%dual > coarse%dual_upper) error stop 'coarse lower above upper'
        if (fine%dual > coarse%dual_upper) error stop 'finer lower bound above upper'
        print '(a)', 'PASS: D3 upper bound'
    end subroutine

    subroutine test_readback(order)
        ! Gate of migration M2 beyond the K6-mini fixture: the readback reads
        ! the state's own boundary flux and orientation, never psi = 0.
        integer, intent(in) :: order
        real(dp), parameter :: shift = 0.37_dp
        type(gs_state_t) :: mesh, moved
        type(equilibrium_state_t) :: eq
        type(equilibrium_readback_t) :: base, shifted, flipped, inner
        logical, allocatable :: plasma(:)
        real(dp) :: psi_b
        integer :: status, k, e
        call solovev_mesh(48, order, mesh)
        call solve_gs(mesh, 1.0_dp, -16*c_sol, 0.0_dp, status)
        if (status /= 0) error stop 'Solov''ev solve'
        call gs_equilibrium(mesh, eq, status, f_tor, 0.0_dp, pprime=-16*c_sol)
        if (status /= 0) error stop 'psi adapter'
        call read_equilibrium(eq, 65, 0, base, status)
        if (status /= 0) error stop 'fixed-boundary readback'
        ! Axis gauge: the same field with psi = 0 on the symmetry axis.
        moved = mesh
        moved%psi = mesh%psi + shift
        call gs_equilibrium(moved, eq, status, f_tor, 0.0_dp, pprime=-16*c_sol, &
            psi_bound=shift)
        if (status /= 0) error stop 'shifted adapter'
        call read_equilibrium(eq, 65, 0, shifted, status)
        if (status /= 0) error stop 'shifted readback'
        call check('gauge shift q', maxval(abs(shifted%q/base%q - 1)), 1.0e-12_dp)
        call check('gauge shift s_tor', maxval(abs(shifted%s_tor - base%s_tor)), &
            1.0e-12_dp)
        call check('gauge shift psi', maxval(abs(shifted%psi - shift - base%psi)), &
            1.0e-12_dp)
        call check('gauge shift volume', abs(shifted%volume/base%volume - 1), &
            1.0e-14_dp)
        call check('contour vs area toroidal flux', base%flux_consistency, 1.0e-5_dp)
        ! Reversed flux (and poloidal field): q changes sign, s_tor does not.
        moved%psi = -mesh%psi
        call gs_equilibrium(moved, eq, status, f_tor, 0.0_dp, pprime=16*c_sol)
        if (status /= 0) error stop 'flipped adapter'
        call read_equilibrium(eq, 65, 0, flipped, status)
        if (status /= 0) error stop 'flipped readback'
        call check('flip q', maxval(abs(flipped%q/base%q + 1)), 1.0e-12_dp)
        call check('flip s_tor', maxval(abs(flipped%s_tor - base%s_tor)), 1.0e-12_dp)
        call check('flip axis flux', abs(flipped%psi_axis + base%psi_axis), &
            1.0e-14_dp)
        ! An interior level as plasma boundary on the cells reaching it.
        k = 33
        psi_b = base%psi(k)
        allocate (plasma(size(mesh%triangles, 2)))
        do e = 1, size(plasma)
            plasma(e) = minval(mesh%psi(mesh%triangles(:, e))) &
                <= psi_b - 0.05_dp*base%psi_axis
        end do
        if (all(plasma)) error stop 'sub-mesh is the whole mesh'
        call gs_equilibrium(mesh, eq, status, f_tor, 0.0_dp, pprime=-16*c_sol, &
            psi_bound=psi_b, plasma_cell=plasma)
        if (status /= 0) error stop 'sub-mesh adapter'
        call read_equilibrium(eq, 65, 0, inner, status)
        if (status /= 0) error stop 'sub-mesh readback'
        call check('sub-mesh edge q', abs(inner%q(65)/base%q(k) - 1), 1.0e-10_dp)
        call check('sub-mesh edge level', abs(inner%psi(65) - psi_b), 1.0e-14_dp)
        if (.not. ieee_is_nan(inner%volume)) error stop 'sub-mesh volume must be NaN'
        print '(a)', 'PASS: D7 readback'
    end subroutine

    subroutine check(label, value, tolerance)
        character(*), intent(in) :: label
        real(dp), intent(in) :: value, tolerance
        print '(a30,es10.3)', label, value
        if (.not. value <= tolerance) error stop 'D7 readback check failed'
    end subroutine

    subroutine solovev_mesh(n, order, mesh)
        integer, intent(in) :: n, order
        type(gs_state_t), intent(out) :: mesh
        real(dp) :: contour(2, n), h
        integer :: j, status
        do j = 1, n
            call solovev_boundary(real(j - 1, dp)/n, contour(:, j))
        end do
        h = 2*acos(-1.0_dp)*minor/n
        call initialize_gs(mesh, contour, h*h/2, status, order, &
            curve=solovev_boundary)
        if (status /= 0) error stop 'Solov''ev mesh'
    end subroutine

    subroutine solovev_boundary(t, point)
        real(dp), intent(in) :: t
        real(dp), intent(out) :: point(2)
        real(dp) :: theta
        theta = 2*acos(-1.0_dp)*t
        point(1) = sqrt(r0*r0 + 2*r0*minor*cos(theta))
        point(2) = r0*minor*sin(theta)/point(1)
    end subroutine

    subroutine test_chart()
        ! Second and third q-derivatives are the derivatives of the lower ones.
        type(cylindrical_chart_t) :: chart
        real(dp) :: q(3), x(3), d1(3, 3), d2(3, 3, 3), d3(3, 3, 3, 3)
        real(dp) :: xp(3), d1p(3, 3), d2p(3, 3, 3), d3p(3, 3, 3, 3)
        real(dp) :: xm(3), d1m(3, 3), d2m(3, 3, 3), d3m(3, 3, 3, 3), dq(3)
        real(dp) :: error, h
        integer :: k, sample
        error = 0
        h = 1.0e-5_dp
        do sample = 1, 4
            q = [6.2_dp + 0.3_dp*sample, -0.5_dp + 0.4_dp*sample, 0.7_dp*sample]
            call chart%jet(q, x, d1, d2, d3)
            if (maxval(abs(x - [q(1)*cos(q(3)), q(1)*sin(q(3)), q(2)])) > 1.0e-14_dp) &
                error stop 'cylindrical chart position'
            do k = 1, 3
                dq = 0; dq(k) = h
                call chart%jet(q + dq, xp, d1p, d2p, d3p)
                call chart%jet(q - dq, xm, d1m, d2m, d3m)
                error = max(error, maxval(abs((xp - xm)/(2*h) - d1(:, k))), &
                    maxval(abs((d1p - d1m)/(2*h) - d2(:, :, k))), &
                    maxval(abs((d2p - d2m)/(2*h) - d3(:, :, :, k))))
            end do
        end do
        print '(a,es10.3)', 'cylindrical chart jet finite-difference defect: ', error
        if (error > 1.0e-8_dp) error stop 'cylindrical chart jet inconsistent'
        print '(a)', 'PASS: chart jet'
    end subroutine

    subroutine compare_k6_mini(work, fixtures)
        character(*), intent(in) :: work, fixtures
        real(dp), parameter :: tolerance = 1.0e-12_dp
        real(dp) :: worst
        worst = 0
        call compare_summary(work//'/summary.dat', &
            fixtures//'/reference_summary.dat', worst)
        call compare_table(work//'/profiles.dat', &
            fixtures//'/reference_profiles.dat', 7, 1, worst)
        call compare_table(work//'/field.dat', &
            fixtures//'/reference_field.dat', 4, 0, worst)
        call check_synthetic_analytic(work, fixtures)
        print '(a,es10.3)', 'K6-mini worst relative deviation: ', worst
        if (worst > tolerance) error stop 'K6-mini regression exceeds 1e-12'
        print '(a)', 'PASS: K6-mini regression'
    end subroutine

    subroutine check_synthetic_analytic(work, fixtures)
        ! Independent exact Solov'ev oracle, not the retained numerical output.
        ! psi=C*((R^2-R0^2)^2+4*R^2*Z^2)-4*C*R0^2*a^2,
        ! Delta*psi=16*C*R^2, p'=-16*C/mu0, FF'=0, F=6.
        character(*), intent(in) :: work, fixtures
        real(dp), parameter :: centre = 3, minor = 0.6_dp, amplitude = 0.25_dp
        real(dp), allocatable :: points(:, :), fields(:, :)
        real(dp) :: r, z, expected(4), scale(4), worst
        integer :: j
        call read_table(fixtures//'/points.dat', 2, 1, points)
        call read_table(work//'/field.dat', 4, 0, fields)
        if (size(points, 2) /= size(fields, 2)) error stop 'analytic query count'
        if (.not. all(ieee_is_finite(fields))) error stop 'nonfinite synthetic field'
        scale = [4*amplitude*centre**2*minor**2, &
            8*amplitude*centre*minor, 8*amplitude*centre*minor, 6/centre]
        worst = 0
        do j = 1, size(points, 2)
            r = points(1, j); z = points(2, j)
            expected = [amplitude*((r*r - centre**2)**2 + 4*r*r*z*z) &
                - 4*amplitude*centre**2*minor**2, -8*amplitude*r*z, &
                4*amplitude*(r*r - centre**2 + 2*z*z), 6/r]
            worst = max(worst, maxval(abs(fields(:, j) - expected)/scale))
        end do
        print '(a,es10.3)', 'synthetic Solov''ev analytic scaled error: ', worst
        ! Coarse P3 mesh/curved-boundary discretization gate (0.1% of the
        ! declared physical scales), separate from version regression drift.
        if (worst > 1.0e-3_dp) error stop 'synthetic analytic field mismatch'
    end subroutine

    subroutine compare_summary(path, reference, worst)
        ! Every key of the reference must be present; values relative.
        character(*), intent(in) :: path, reference
        real(dp), intent(inout) :: worst
        character(64) :: key, other
        character(256) :: line
        real(dp) :: expected, value, deviation
        integer :: unit, ref_unit, stat
        real(dp) :: psi_scale, radius_scale
        logical :: found
        open (newunit=ref_unit, file=reference, status='old', action='read')
        psi_scale = 0; radius_scale = 0
        do
            read (ref_unit, '(a)', iostat=stat) line
            if (stat /= 0) exit
            if (len_trim(line) == 0) cycle
            read (line, *) key, expected
            if (trim(key) == 'psi_axis') psi_scale = abs(expected)
            if (trim(key) == 'axis_R') radius_scale = abs(expected)
        end do
        if (psi_scale == 0) error stop 'K6-mini reference lacks psi_axis'
        rewind (ref_unit)
        do
            read (ref_unit, '(a)', iostat=stat) line
            if (stat /= 0) exit
            if (len_trim(line) == 0) cycle
            read (line, *) key, expected
            found = .false.
            open (newunit=unit, file=path, status='old', action='read')
            do
                read (unit, '(a)', iostat=stat) line
                if (stat /= 0) exit
                read (line, *) other
                if (trim(other) /= trim(key)) cycle
                read (line, *) other, value
                found = .true.
                exit
            end do
            close (unit)
            if (.not. found) then
                print '(2a)', 'missing summary key: ', trim(key)
                error stop 'K6-mini summary incomplete'
            end if
            ! flux_consistency is a difference |a/b - 1| ~ 1e-4: absolute.
            if (trim(key) == 'axis_Z') then
                ! Up/down symmetry makes Z_axis zero; use the major radius
                ! as its physical scale, not roundoff in the reference Z.
                deviation = abs(value - expected)/max(radius_scale, tiny(1.0_dp))
            else if (trim(key) == 'flux_consistency') then
                deviation = abs(value - expected)
            else if (any(trim(key) == [character(16) :: 'majorant', &
                'majorant_flux', 'majorant_balance'])) then
                ! A-posteriori estimators are differences of solution-sized
                ! quantities, so their
                ! absolute noise is the solve's, not relative to themselves.
                ! The solve stops at nonlinear_tolerance 1e-12, and the fast
                ! reductions/FMA move the iterate at that level. Measure
                ! against the solution scale |psi_axis|: the observed spread
                ! is ~1e-14 there, 100x below the gate, while a 1e-8 change
                ! of the solution still gives ~1e-8.
                deviation = abs(value - expected)/max(psi_scale, tiny(1.0_dp))
            else
                deviation = abs(value - expected)/max(abs(expected), tiny(1.0_dp))
            end if
            print '(a24,2es25.16e3,es10.3)', key, expected, value, deviation
            worst = max(worst, deviation)
        end do
        close (ref_unit)
    end subroutine

    subroutine compare_table(path, reference, columns, header, worst)
        ! Columnwise normwise deviation max|a - b|/max|b|: the fast (default)
        ! reductions leave 1e-14 run-to-run noise in the flux coefficients,
        ! which a pointwise ratio amplifies where a component crosses zero.
        character(*), intent(in) :: path, reference
        integer, intent(in) :: columns, header
        real(dp), intent(inout) :: worst
        real(dp), allocatable :: a(:, :), b(:, :)
        real(dp) :: deviation
        integer :: j
        call read_table(path, columns, header, a)
        call read_table(reference, columns, header, b)
        if (size(a, 2) /= size(b, 2)) error stop 'K6-mini table length differs'
        do j = 1, columns
            deviation = maxval(abs(a(j, :) - b(j, :)))/max(maxval(abs(b(j, :))), &
                tiny(1.0_dp))
            print '(a,a,i0,es10.3)', trim(reference(index(reference, '/', &
                back=.true.) + 1:)), ' column ', j, deviation
            worst = max(worst, deviation)
        end do
    end subroutine

    subroutine read_table(path, columns, header, table)
        character(*), intent(in) :: path
        integer, intent(in) :: columns, header
        real(dp), allocatable, intent(out) :: table(:, :)
        real(dp) :: row(columns)
        integer :: unit, stat, n, i
        open (newunit=unit, file=path, status='old', action='read')
        do i = 1, header
            read (unit, *)
        end do
        n = 0
        do
            read (unit, *, iostat=stat) row
            if (stat /= 0) exit
            n = n + 1
        end do
        rewind (unit)
        do i = 1, header
            read (unit, *)
        end do
        allocate (table(columns, n))
        do i = 1, n
            read (unit, *) table(:, i)
        end do
        close (unit)
        if (n == 0) error stop 'empty K6-mini table'
    end subroutine
end program test_equilibrium_core
