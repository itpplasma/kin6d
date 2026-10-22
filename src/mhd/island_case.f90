module kin6d_island_case
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use kin6d_nonnested_coupled, only: nc_problem_t, nc_validate
    use kin6d_island_closure, only: ISLAND_INVALID
    implicit none
    private
    public :: island_case_read
contains
    subroutine island_case_read(unit, problem, status, message)
        ! Source hooks and boundary values are supplied by the case adapter.
        ! Read transactionally: a rejected declaration never changes problem.
        integer, intent(in) :: unit
        type(nc_problem_t), intent(inout) :: problem
        integer, intent(out) :: status
        character(*), intent(out) :: message
        type(nc_problem_t) :: trial
        character(24) :: pressure_closure(1), current_closure(1), drive
        character(24) :: heat_boundary, potential_boundary, target_class
        real(dp) :: chi_par, chi_perp, g_par, g_perp, I_phi, V
        real(dp), parameter :: absent = huge(1.0_dp)
        integer :: ios
        logical :: has_current, has_voltage
        namelist /island_case/ pressure_closure, current_closure, drive, &
            heat_boundary, potential_boundary, target_class, &
            chi_par, chi_perp, g_par, g_perp, I_phi, V
        trial = problem
        pressure_closure = trial%closure%pressure
        current_closure = trial%closure%current
        drive = 'auto'
        heat_boundary = trial%closure%heat_boundary
        potential_boundary = trial%closure%potential_boundary
        target_class = trial%closure%target_class
        chi_par = trial%chi_par; chi_perp = trial%chi_perp
        g_par = trial%g_par; g_perp = trial%g_perp
        I_phi = absent; V = absent
        read (unit, nml=island_case, iostat=ios, iomsg=message)
        status = ISLAND_INVALID
        if (ios /= 0) return
        ! Compare the exact reserved sentinel bit pattern, not physical values.
        has_current = transfer(I_phi, 0_int64) /= transfer(absent, 0_int64)
        has_voltage = transfer(V, 0_int64) /= transfer(absent, 0_int64)
        message = 'Do not prescribe both I_phi and V'
        if (has_current .and. has_voltage) return
        message = 'Drive selector conflicts with supplied I_phi or V'
        if (has_voltage .and. drive /= 'auto' .and. drive /= 'V') return
        if (has_current .and. drive /= 'auto' .and. drive /= 'I_phi') return
        if (has_voltage) then
            drive = 'V'; trial%voltage = V
        else if (has_current) then
            drive = 'I_phi'; trial%current = I_phi
        end if
        if (drive == 'auto') drive = trial%closure%drive
        trial%closure%pressure = pressure_closure(1)
        trial%closure%current = current_closure(1)
        trial%closure%drive = drive
        trial%closure%heat_boundary = heat_boundary
        trial%closure%potential_boundary = potential_boundary
        trial%closure%target_class = target_class
        trial%chi_par = chi_par; trial%chi_perp = chi_perp
        trial%g_par = g_par; trial%g_perp = g_perp
        call nc_validate(trial, status, message)
        if (status == 0) problem = trial
    end subroutine
end module
