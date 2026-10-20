program test_equilibrium_harness
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use kin6d_equilibrium_harness, only: equilibrium_case_t, &
        equilibrium_resolution_t, equilibrium_start_t, equilibrium_run_report_t, &
        equilibrium_case, equilibrium_start, run_equilibrium_case, &
        D8_BAD_INPUT, D8_BAD_REPORT
    implicit none
    type(equilibrium_case_t) :: spec
    type(equilibrium_start_t) :: start
    type(equilibrium_run_report_t) :: report
    integer :: state, status
    state = 0
    call equilibrium_case('solovev', 'psi', spec, status)
    if (status /= D8_BAD_INPUT) error stop 'unsupported lane admitted'
    call equilibrium_start('zp', start, status)
    if (status /= 0) error stop 'vertical offset missing'
    if (maxval(abs(start%displacement - [0.0_dp, 0.015_dp])) > 1e-15_dp) &
        error stop 'vertical displacement units'
    call equilibrium_start('el', start, status)
    if (abs(start%elongation - 0.02_dp) > 1e-15_dp) &
        error stop 'elongation fraction'
    call run_equilibrium_case('solovev', 'psi', equilibrium_resolution_t(), &
        'legacy', solve, state, report)
    if (state /= 0 .or. report%status /= D8_BAD_INPUT) &
        error stop 'invalid case ran the solver'
    call run_equilibrium_case('solovev', 'map', equilibrium_resolution_t(1, 24), &
        'legacy', solve, state, report)
    if (state /= 0 .or. report%status /= D8_BAD_INPUT) &
        error stop 'unsupported resolution ran the solver'
    call run_equilibrium_case('solovev', 'map', equilibrium_resolution_t(), &
        'legacy', solve, state, report)
    if (state /= 1 .or. report%status /= 0) error stop 'native state not returned'
    if (report%claim /= 'candidate' .or. report%state_bound_status /= 'unavailable') &
        error stop 'oracle comparison promoted to a certificate'
    call run_equilibrium_case('solovev', 'map', equilibrium_resolution_t(), &
        'legacy', incomplete, state, report)
    if (report%status /= D8_BAD_REPORT) error stop 'missing evidence accepted'
    state = -1
    call run_equilibrium_case('solovev', 'map', equilibrium_resolution_t(), &
        'legacy', solve, state, report)
    if (report%status /= D8_BAD_REPORT) error stop 'NaN accepted'
contains
    subroutine solve(spec, resolution, start, state, report)
        type(equilibrium_case_t), intent(in) :: spec
        type(equilibrium_resolution_t), intent(in) :: resolution
        type(equilibrium_start_t), intent(in) :: start
        class(*), intent(inout) :: state
        type(equilibrium_run_report_t), intent(inout) :: report
        if (spec%lane /= 'map' .or. resolution%degree /= 2) return
        if (start%name /= 'legacy') return
        select type (state)
            type is (integer)
            report%oracle_error = 0.01_dp
            if (state < 0) report%oracle_error = &
                ieee_value(0.0_dp, ieee_quiet_nan)
            state = state + 1
        end select
        report%seed_error = 1
        report%residual = 1e-8_dp
        report%residual_norm = 'test oracle'
        report%status = 0
    end subroutine

    subroutine incomplete(spec, resolution, start, state, report)
        type(equilibrium_case_t), intent(in) :: spec
        type(equilibrium_resolution_t), intent(in) :: resolution
        type(equilibrium_start_t), intent(in) :: start
        class(*), intent(inout) :: state
        type(equilibrium_run_report_t), intent(inout) :: report
        ! A solver returning success without a report must fail closed.
        call solve(spec, resolution, start, state, report)
        report%residual_norm = 'unavailable'
    end subroutine
end program
