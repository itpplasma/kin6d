module kin6d_equilibrium_harness
    ! D8/M6: one registry and reporting seam for the existing lane drivers.
    ! A driver owns its native state (including frozen labels / coil data).
    ! The harness never changes physics, differentiates a solver iteration,
    ! or promotes an oracle comparison to a certificate.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use kin6d_equilibrium_state, only: equilibrium_state_t
    use kin6d_triangle_basis, only: TRIANGLE_MAX_DEGREE
    implicit none
    private
    public :: equilibrium_case_t, equilibrium_resolution_t, equilibrium_start_t
    public :: equilibrium_run_report_t, equilibrium_case, equilibrium_start
    public :: run_equilibrium_case, write_equilibrium_report
    public :: equilibrium_floor_report_t, equilibrium_residual_floor
    public :: write_residual_floor
    integer, parameter, public :: D8_OK = 0, D8_BAD_INPUT = -1
    integer, parameter, public :: D8_BAD_REPORT = -2

    type :: equilibrium_case_t
        character(32) :: name, lane
        character(80) :: oracle
    end type
    type(equilibrium_case_t), parameter :: cases(*) = [ &
        equilibrium_case_t('solovev', 'map', 'Solovev A3; independent 2D GS'), &
        equilibrium_case_t('cerfon-iter-0', 'psi', 'Cerfon-Freidberg symmetric ITER A=0'), &
        equilibrium_case_t('cerfon-iter-155', 'psi', 'Cerfon-Freidberg symmetric ITER A=-0.155'), &
        equilibrium_case_t('cerfon-nstx-0', 'psi', 'Cerfon-Freidberg symmetric NSTX A=0'), &
        equilibrium_case_t('cerfon-nstx-1', 'psi', 'Cerfon-Freidberg symmetric NSTX A=1'), &
        equilibrium_case_t('shaped', 'map', 'shaped Solovev; independent 2D GS'), &
        equilibrium_case_t('seed', 'map', 'nested A3; exact flux and 2D GS'), &
        equilibrium_case_t('seedshaped', 'map', 'nested shaped; exact flux and 2D GS'), &
        equilibrium_case_t('limiter-jump', 'psi', 'manufactured compact flux bump'), &
        equilibrium_case_t('limiter-smooth', 'psi', 'manufactured compact flux bump'), &
        equilibrium_case_t('targets-jump', 'psi', 'manufactured compact flux bump'), &
        equilibrium_case_t('targets-smooth', 'psi', 'manufactured compact flux bump'), &
        equilibrium_case_t('xpoint-jump', 'psi', 'manufactured two-bump FB-cut'), &
        equilibrium_case_t('separatrix-smooth', 'psi', 'manufactured two-bump separatrix')]

    type :: equilibrium_resolution_t
        integer :: degree = 2, boundary_nodes = 24
    end type
    type :: equilibrium_start_t
        character(16) :: name = 'legacy'
        real(dp) :: displacement(2) = [0.10_dp, -0.08_dp]
        real(dp) :: elongation = 0, amplitude = 0.3_dp
    end type
    type :: equilibrium_run_report_t
        integer :: status = D8_BAD_REPORT, iterations = 0
        character(32) :: case_name = '', lane = '', claim = 'candidate'
        character(16) :: start = ''
        character(80) :: oracle = ''
        integer :: degree = 0, boundary_nodes = 0
        real(dp) :: oracle_error = huge(1.0_dp), seed_error = huge(1.0_dp)
        real(dp) :: residual = huge(1.0_dp), sigma = 0
        character(64) :: residual_norm = 'unavailable'
        ! Bounds are unavailable unless a lane has actually evaluated them.
        ! Oracle errors and sampled residuals do not fill a bound field.
        character(32) :: state_bound_status = 'unavailable'
        character(32) :: derivative_bound_status = 'unavailable'
    end type
    ! Residual floor of a refinement ladder (docs/relaxation-comparators.md,
    ! item 5): current sheets at rational surfaces set a force-residual
    ! floor that h or p refinement does not lower. A ladder is 'decreasing'
    ! while every refinement reduces the residual by at least the factor
    ! ratio**(-min_order); otherwise the finest level is reported as the
    ! floor. The orders are observed rates log(r_i-1/r_i)/log(ratio), with
    ! ratio the nominal refinement factor per level (2 for h doubling; a p
    ! ladder states its own nominal factor). This is a diagnostic: a gate
    ! reports the floor instead of chasing it, and nothing is certified.
    type :: equilibrium_floor_report_t
        integer :: levels = 0
        real(dp) :: ratio = 2, min_order = 0.5_dp, floor = 0, last_order = 0
        real(dp), allocatable :: residuals(:), orders(:)
        character(16) :: status = 'unavailable' ! decreasing | floor
        character(64) :: residual_norm = 'unavailable'
    end type
    abstract interface
        subroutine lane_run(spec, resolution, start, state, report, equilibrium)
            import :: equilibrium_case_t, equilibrium_resolution_t
            import :: equilibrium_start_t, equilibrium_run_report_t
            import :: equilibrium_state_t
            type(equilibrium_case_t), intent(in) :: spec
            type(equilibrium_resolution_t), intent(in) :: resolution
            type(equilibrium_start_t), intent(in) :: start
            class(*), intent(inout) :: state
            type(equilibrium_run_report_t), intent(inout) :: report
            type(equilibrium_state_t), intent(out) :: equilibrium
        end subroutine
    end interface
contains
    subroutine equilibrium_case(name, lane, spec, status)
        character(*), intent(in) :: name, lane
        type(equilibrium_case_t), intent(out) :: spec
        integer, intent(out) :: status
        integer :: i
        status = D8_BAD_INPUT
        spec = equilibrium_case_t('', '', '')
        do i = 1, size(cases)
            if (name /= cases(i)%name .or. lane /= cases(i)%lane) cycle
            spec = cases(i)
            status = D8_OK
            return
        end do
    end subroutine

    subroutine equilibrium_start(name, start, status)
        character(*), intent(in) :: name
        type(equilibrium_start_t), intent(out) :: start
        integer, intent(out) :: status
        start = equilibrium_start_t()
        status = D8_OK
        if (name == 'legacy') return
        start%name = name
        start%displacement = 0
        start%amplitude = 1
        select case (name)
        case ('reference')
        case ('rp')
            start%displacement(1) = 0.015_dp
        case ('rm')
            start%displacement(1) = -0.015_dp
        case ('zp')
            start%displacement(2) = 0.015_dp
        case ('zm')
            start%displacement(2) = -0.015_dp
        case ('el')
            start%elongation = 0.02_dp
        case default
            status = D8_BAD_INPUT
        end select
    end subroutine

    subroutine run_equilibrium_case(name, lane, resolution, start_name, solve, &
            state, report, equilibrium)
        character(*), intent(in) :: name, lane, start_name
        type(equilibrium_resolution_t), intent(in) :: resolution
        procedure(lane_run) :: solve
        class(*), intent(inout) :: state
        type(equilibrium_run_report_t), intent(out) :: report
        type(equilibrium_state_t), intent(out) :: equilibrium
        type(equilibrium_case_t) :: spec
        type(equilibrium_start_t) :: start
        integer :: status
        report = equilibrium_run_report_t()
        report%case_name = name; report%lane = lane; report%start = start_name
        report%degree = resolution%degree
        report%boundary_nodes = resolution%boundary_nodes
        report%status = D8_BAD_INPUT
        call equilibrium_case(name, lane, spec, status)
        if (status /= D8_OK) return
        report%oracle = spec%oracle
        call equilibrium_start(start_name, start, status)
        if (status /= D8_OK) return
        ! The psi lane needs the GS solver (P2/P3); the isoparametric map lane
        ! admits every Lagrange degree of the triangle basis.
        if (resolution%degree < 2) return
        if (resolution%degree > merge(TRIANGLE_MAX_DEGREE, 3, lane == 'map')) return
        if (resolution%boundary_nodes < merge(8, 4, resolution%degree <= 3)) return
        ! Map seeds are the retained interior pullback or nested-axis seed;
        ! the free-boundary offset registry currently belongs to the psi lane.
        if (lane == 'map' .and. start_name /= 'legacy') return
        report%status = D8_BAD_REPORT
        call solve(spec, resolution, start, state, report, equilibrium)
        if (report%status /= D8_OK) return
        if (.not. allocated(equilibrium%field) &
            .or. .not. allocated(equilibrium%chart)) then
            report%status = D8_BAD_REPORT
            return
        end if
        if (equilibrium%lane /= lane) then
            report%status = D8_BAD_REPORT
            return
        end if
        equilibrium%claim = report%claim
        write (equilibrium%seed, '(4a,i0,a,i0)') name, '/', start_name, &
            ' P', resolution%degree, ' n', resolution%boundary_nodes
        if (.not. ieee_is_finite(report%oracle_error) &
            .or. .not. ieee_is_finite(report%seed_error) &
            .or. .not. ieee_is_finite(report%residual) &
            .or. .not. ieee_is_finite(report%sigma)) then
            report%status = D8_BAD_REPORT
            return
        end if
        if (min(report%oracle_error, report%seed_error, report%residual) < 0 &
            .or. report%iterations < 0 .or. report%residual_norm == 'unavailable') &
            report%status = D8_BAD_REPORT
        if (report%oracle_error >= huge(1.0_dp) &
            .or. report%seed_error >= huge(1.0_dp) &
            .or. report%residual >= huge(1.0_dp)) report%status = D8_BAD_REPORT
    end subroutine

    subroutine write_equilibrium_report(unit, report)
        integer, intent(in) :: unit
        type(equilibrium_run_report_t), intent(in) :: report
        write (unit, '(6a,3(a,i0))') '[V] D8 case=', trim(report%case_name), &
            ' lane=', trim(report%lane), ' start=', trim(report%start), &
            ' degree=', report%degree, ' nodes=', report%boundary_nodes, &
            ' status=', report%status
        write (unit, '(2a)') '[V] oracle=', trim(report%oracle)
        write (unit, '(a,i0,3(a,es16.8))') '[V] iterations=', report%iterations, &
            ' oracle_error=', report%oracle_error, ' seed_error=', report%seed_error, &
            ' residual=', report%residual
        write (unit, '(5a,es16.8)') '[V] claim=', trim(report%claim), &
            ' residual_norm=', trim(report%residual_norm), ' sigma=', report%sigma
        write (unit, '(4a)') '[V] state_bound=', trim(report%state_bound_status), &
            ' derivative_bound=', trim(report%derivative_bound_status)
    end subroutine

    subroutine equilibrium_residual_floor(residuals, ratio, residual_norm, report, &
            min_order)
        real(dp), intent(in) :: residuals(:), ratio
        character(*), intent(in) :: residual_norm
        type(equilibrium_floor_report_t), intent(out) :: report
        real(dp), intent(in), optional :: min_order
        integer :: i, m
        m = size(residuals)
        report%levels = m; report%ratio = ratio
        report%residual_norm = residual_norm
        if (present(min_order)) report%min_order = min_order
        report%residuals = residuals
        allocate (report%orders(max(m - 1, 0))); report%orders = 0
        if (m < 2 .or. .not. ratio > 1) return
        if (.not. all(ieee_is_finite(residuals)) .or. any(residuals <= 0)) return
        do i = 2, m
            report%orders(i - 1) = log(residuals(i - 1)/residuals(i))/log(ratio)
        end do
        report%last_order = report%orders(m - 1)
        report%status = 'decreasing'
        if (minval(report%orders) < report%min_order) then
            report%status = 'floor'
            report%floor = residuals(m)
        end if
    end subroutine

    subroutine write_residual_floor(unit, label, report)
        integer, intent(in) :: unit
        character(*), intent(in) :: label
        type(equilibrium_floor_report_t), intent(in) :: report
        character(32) :: fmt
        write (unit, '(5a,i0,a,f5.2,a,es12.4)') '[V] residual ladder ', label, &
            ' norm=', trim(report%residual_norm), ' levels=', report%levels, &
            ' min_order=', report%min_order, ' floor=', report%floor
        write (fmt, '(a,i0,a)') '(a,', max(report%levels, 1), 'es12.4)'
        write (unit, fmt) '[V]   residuals', report%residuals
        if (report%levels > 1) then
            write (fmt, '(a,i0,a)') '(a,', report%levels - 1, 'f8.3,2a)'
            write (unit, fmt) '[V]   orders', report%orders, '  status=', &
                trim(report%status)
        end if
    end subroutine
end module kin6d_equilibrium_harness
