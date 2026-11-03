program test_vrmhd_spectra
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use vrmhd_spectrum_support, only: spectrum, require
    use rm_mhd_oracles, only: uniform_waves, oracle_mu0, oracle_ok
    implicit none
    integer, parameter :: max_levels = 4
    real(dp) :: frequencies(3, max_levels), errors(3, max_levels)
    real(dp) :: rates(3), exact(3), b0(3)
    real(dp) :: optical(max_levels)
    real(dp) :: mass_condition(max_levels), condition_growth, refinement
    integer, allocatable :: cells(:)
    integer :: level, previous, levels, status, unmatched, ka, ks, mesh, repeat
    integer :: rate_base
    character(32) :: selection
    logical :: fine
    call get_command_argument(1, selection)
    call require(command_argument_count() <= 1, 'one optional ladder selection')
    call require(selection == '' .or. selection == 'fine', 'known ladder selection')
    fine = selection == 'fine'
    b0 = sqrt(oracle_mu0)*[0.6_dp, 0.0_dp, 0.8_dp]
    call uniform_waves([1.0_dp, 0.0_dp, 0.0_dp], b0, 1.0_dp, 0.6_dp, &
        5.0_dp/3, exact, status)
    call require(status == oracle_ok, 'existing uniform-wave oracle')
    exact = sqrt(exact)
    print '(a,3es20.12)', 'oracle Alfven / fast / slow angular frequencies [1/s]: ', exact
    do mesh = 1, 2
        repeat = mesh
        print '(a,l1)', 'periodically perturbed mesh: ', mesh == 2
        do ka = merge(2, 1, fine), 2
            ks = ka
            if (fine) then
                cells = [32, 48, 64]
            else if (ka == 1) then
                cells = [4, 8, 16, 32]
            else
                cells = [8, 16, 32]
            end if
            levels = size(cells)
            print '(a,2i4)', 'edge order / scalar degree: ', ka, ks
            do level = 1, levels
                call spectrum(cells(level), ka, ks, &
                    repeat**2*(6*ks**2 + 2*(3*ka**2 + ka)), &
                    b0, exact, .true., &
                    frequencies(:, level), optical(level), unmatched, &
                    perturbed=mesh == 2, structural=.true., &
                    mass_condition=mass_condition(level))
                errors(:, level) = abs(frequencies(:, level)/exact - 1)
                print '(a,i4,3es20.12)', 'wave frequencies [1/s], cells per period: ', &
                    cells(level), frequencies(:, level)
                print '(a,i4,3es20.12)', 'relative frequency errors: ', &
                    cells(level), errors(:, level)
                if (level == 1) cycle
                previous = level - 1
                refinement = real(cells(level), dp)/cells(previous)
                condition_growth = (mass_condition(level)/mass_condition(previous)) &
                    **(log(2.0_dp)/log(refinement))
                print '(a,es20.12)', 'M_r condition growth per halving: ', &
                    condition_growth
                call require(condition_growth <= 2.0_dp, &
                    'bounded physical time-mass condition')
                rates = log(errors(:, previous)/errors(:, level))/log(refinement)
                print '(a,3f12.6)', 'Alfven / fast / slow frequency orders: ', rates
            end do
            if (ka == 2 .and. .not. fine) then
                print '(a)', &
                    'degree-two design order requires vrmhd_section_spectra_fine'
                cycle
            end if
            rate_base = merge(1, levels-1, fine)
            refinement = real(cells(levels), dp)/cells(rate_base)
            rates = log(errors(:, rate_base)/errors(:, levels))/log(refinement)
            print '(a,2i4,3f12.6)', 'design-order cells / Alfven / fast / slow: ', &
                cells(rate_base), cells(levels), rates
            call require(minval(rates) > 2*ka - 0.2_dp, 'frequency design order')
        end do
    end do
end program test_vrmhd_spectra
