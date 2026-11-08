!> KAMEL's CGS profile files and the canonical energy-flow result.
module kin6d_dk_energy_profile_io
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_is_finite
    use kin6d_dk_energy, only: energy_axi_t, energy_estimate_t, ENERGY_OK, ENERGY_BAD_INPUT
    use kin6d_dk_force_balance, only: energy_balance_input_t, energy_balance_t, &
        energy_force_balance_cgs
    use kin6d_result_io, only: run_t, result_file_t, result_directory
    implicit none
    private
    public :: energy_profile_write_cgs, PROFILE_OK, PROFILE_INCOMPLETE, &
        PROFILE_IO_ERROR, PROFILE_BAD_INPUT
    integer, parameter :: PROFILE_OK = 0, PROFILE_INCOMPLETE = 1, &
        PROFILE_IO_ERROR = 2, PROFILE_BAD_INPUT = 3
contains
    !> p%r is increasing effective radius [cm]; p%k is replaced by response%k.
    !> Failed rows remain in NetCDF with their status and NaN outputs. Legacy
    !> files are withheld after a failed row; I/O errors have a separate status.
    subroutine energy_profile_write_cgs(run, p, response, solver_status, status, estimate)
        type(run_t), intent(in) :: run
        type(energy_balance_input_t), intent(in) :: p(:)
        type(energy_axi_t), intent(in) :: response(:)
        integer, intent(in) :: solver_status(:)
        integer, intent(out) :: status
        type(energy_estimate_t), intent(in), optional :: estimate(:)
        type(energy_balance_input_t) :: input, direction
        type(energy_balance_t) :: balance, delta_balance
        type(result_file_t) :: file
        real(dp), allocatable :: er(:), vth(:), k(:), d31(:), d32(:), ds(:), &
            d31s(:), d32s(:), backward(:, :), delta(:, :)
        integer, allocatable :: balance_status(:), refinement_status(:)
        real(dp) :: nan
        integer :: n, j, ue, uv, ios, ios2, delta_status
        logical :: erased_er, erased_vth
        character(16) :: export_state
        status = PROFILE_BAD_INPUT
        n = size(p)
        if (n < 1 .or. size(response) /= n .or. size(solver_status) /= n) return
        if (present(estimate)) then
            if (size(estimate) /= n) return
        end if
        if (.not. all(ieee_is_finite(p%r)) .or. any(p%r < 0)) return
        if (any(p(2:)%r <= p(:n - 1)%r)) return
        nan = ieee_value(0.0_dp, ieee_quiet_nan)
        allocate (er(n), vth(n), k(n), d31(n), d32(n), ds(n), d31s(n), d32s(n), &
            backward(2, n), balance_status(n), delta(5, n), refinement_status(n))
        er = nan; vth = nan; k = nan; d31 = nan; d32 = nan; ds = nan
        d31s = nan; d32s = nan; backward = nan
        delta = nan
        balance_status = ENERGY_BAD_INPUT
        refinement_status = ENERGY_BAD_INPUT
        do j = 1, n
            if (solver_status(j) /= ENERGY_OK) cycle
            associate (v => response(j))
                if (.not. all(ieee_is_finite([v%k, v%d31, v%d32, v%dk_ds, &
                    v%dd31_ds, v%dd32_ds, v%residual]))) cycle
                k(j) = v%k; d31(j) = v%d31; d32(j) = v%d32; ds(j) = v%dk_ds
                d31s(j) = v%dd31_ds; d32s(j) = v%dd32_ds; backward(:, j) = v%residual
            end associate
            input = p(j)
            input%k = k(j)
            call energy_force_balance_cgs(input, balance, balance_status(j))
            if (balance_status(j) /= ENERGY_OK) cycle
            er(j) = balance%er; vth(j) = balance%vth
            if (.not. present(estimate)) cycle
            if (.not. all(ieee_is_finite([estimate(j)%k, estimate(j)%d31, estimate(j)%d32]))) cycle
            direction = energy_balance_input_t(k=estimate(j)%k)
            call energy_force_balance_cgs(input, delta_balance, delta_status, direction)
            if (delta_status /= ENERGY_OK) cycle
            delta(:, j) = [estimate(j)%k, estimate(j)%d31, estimate(j)%d32, &
                delta_balance%der, delta_balance%dvth]
            refinement_status(j) = ENERGY_OK
        end do
        status = PROFILE_INCOMPLETE
        export_state = 'incomplete'
        erased_er = remove_legacy('Er.dat')
        erased_vth = remove_legacy('Vth.dat')
        if (.not. erased_er .or. .not. erased_vth) then
            status = PROFILE_IO_ERROR
            export_state = 'io-error'
        else if (all(balance_status == ENERGY_OK)) then
            open (newunit=ue, file=result_directory//'/Er.dat', status='new', &
                action='write', iostat=ios)
            if (ios == 0) then
                open (newunit=uv, file=result_directory//'/Vth.dat', status='new', &
                    action='write', iostat=ios)
                if (ios == 0) then
                    do j = 1, n
                        write (ue, '(2es26.17)', iostat=ios) p(j)%r, er(j)
                        if (ios /= 0) exit
                        write (uv, '(2es26.17)', iostat=ios) p(j)%r, vth(j)
                        if (ios /= 0) exit
                    end do
                    close (uv, iostat=ios2)
                    if (ios2 /= 0) ios = ios2
                end if
                close (ue, iostat=ios2)
                if (ios2 /= 0) ios = ios2
            end if
            if (ios == 0) then
                status = PROFILE_OK
                export_state = 'complete'
            else
                erased_er = remove_legacy('Er.dat')
                erased_vth = remove_legacy('Vth.dat')
                status = PROFILE_IO_ERROR
                export_state = 'io-error'
            end if
        end if
        call file%create(run)
        call file%attribute('certification', 'unavailable')
        call file%attribute('collision_scope', 'single-ion energy-flow response')
        call file%attribute('kamel_export', trim(export_state))
        call file%attribute('radial_derivative', 'increasing effective radius in cm')
        if (present(estimate)) then
            call file%attribute('refinement_indicator', 'signed refined-minus-coarse; not a bound')
        else
            call file%attribute('refinement_indicator', 'unavailable')
        end if
        call file%write('r_eff', p%r, 'cm', 'effective minor radius', dims='surface')
        call file%write('n', p%n, 'cm^-3', 'ion number density', dims='surface')
        call file%write('T', p%t, 'erg', 'ion temperature', dims='surface')
        call file%write('dn_dr', p%nr, 'cm^-4', 'radial density derivative', dims='surface')
        call file%write('dT_dr', p%tr, 'erg/cm', 'radial temperature derivative', dims='surface')
        call file%write('B', p%b, 'G', 'magnetic-field strength', dims='surface')
        call file%write('Bz', p%bz, 'G', 'signed toroidal magnetic field', dims='surface')
        call file%write('Vz', p%vz, 'cm/s', 'signed toroidal velocity', dims='surface')
        call file%write('q', p%q, '1', 'signed safety factor', dims='surface')
        call file%write('R0', p%r0, 'cm', 'major radius', dims='surface')
        call file%write('charge', p%charge, 'statC', 'ion charge', dims='surface')
        call file%write('k', k, '1', 'poloidal-rotation coefficient', 'surface', 'unavailable')
        call file%write('D31', d31, '1', 'normalised density-driven flow', 'surface', 'unavailable')
        call file%write('D32', d32, '1', 'normalised temperature-driven flow', 'surface', 'unavailable')
        call file%write('dk_ds', ds, '1', 'tangent with s=ln(n/T^2)', 'surface', 'unavailable')
        call file%write('dD31_ds', d31s, '1', 'flow tangent with s=ln(n/T^2)', 'surface', 'unavailable')
        call file%write('dD32_ds', d32s, '1', 'flow tangent with s=ln(n/T^2)', 'surface', 'unavailable')
        call file%write('solve_backward_error', backward, '1', 'componentwise discrete backward error', &
            'drive,surface', 'not a continuous error bound')
        call file%write('Er', er, 'statV/cm', 'radial electric field', 'surface', 'unavailable')
        call file%write('Vth', vth, 'cm/s', 'poloidal velocity', 'surface', 'unavailable')
        call file%write('k_refinement_delta', delta(1, :), '1', 'refined-minus-coarse k', &
            'surface', 'empirical refinement difference')
        call file%write('D31_refinement_delta', delta(2, :), '1', 'refined-minus-coarse D31', &
            'surface', 'empirical refinement difference')
        call file%write('D32_refinement_delta', delta(3, :), '1', 'refined-minus-coarse D32', &
            'surface', 'empirical refinement difference')
        call file%write('Er_refinement_delta', delta(4, :), 'statV/cm', 'refined-minus-coarse Er', &
            'surface', 'empirical refinement difference')
        call file%write('Vth_refinement_delta', delta(5, :), 'cm/s', 'refined-minus-coarse poloidal velocity', &
            'surface', 'empirical refinement difference')
        call file%write('solver_status', real(solver_status, dp), '1', &
            'ENERGY_OK=0; other codes indicate failed energy solve', 'surface')
        call file%write('force_balance_status', real(balance_status, dp), '1', &
            'ENERGY_OK=0; other codes indicate unavailable force balance', 'surface')
        call file%write('refinement_status', real(refinement_status, dp), '1', &
            '0=empirical difference available; 1=unavailable', 'surface')
        call file%close()
    end subroutine energy_profile_write_cgs

    logical function remove_legacy(name) result(ok)
        character(*), intent(in) :: name
        logical :: exists
        integer :: unit, ios
        inquire (file=result_directory//'/'//name, exist=exists)
        ok = .not. exists
        if (.not. exists) return
        open (newunit=unit, file=result_directory//'/'//name, status='old', iostat=ios)
        if (ios /= 0) return
        close (unit, status='delete', iostat=ios)
        ok = ios == 0
    end function remove_legacy
end module kin6d_dk_energy_profile_io
