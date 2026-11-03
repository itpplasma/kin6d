program test_vrmhd_spectra
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use vrmhd_spectrum_support, only: spectrum, require
    use rm_mhd_oracles, only: uniform_waves, oracle_mu0, oracle_ok
    implicit none
    integer, parameter :: levels = 4
    real(dp) :: frequencies(3, levels), errors(3, levels), rates(3), exact(3), b0(3)
    real(dp) :: optical(levels)
    integer :: level, status, unmatched, ka, ks, mesh, repeat
    b0 = sqrt(oracle_mu0)*[0.6_dp, 0.0_dp, 0.8_dp]
    call uniform_waves([1.0_dp, 0.0_dp, 0.0_dp], b0, 1.0_dp, 0.6_dp, &
        5.0_dp/3, exact, status)
    call require(status == oracle_ok, 'existing uniform-wave oracle')
    exact = sqrt(exact)
    print '(a,3es20.12)', 'oracle Alfven / fast / slow angular frequencies [1/s]: ', exact
    do mesh = 1, 2
        repeat = mesh
        print '(a,l1)', 'periodically perturbed mesh: ', mesh == 2
        do ka = 1, 2
            ks = ka
            print '(a,2i4)', 'edge order / scalar degree: ', ka, ks
            do level = 1, levels
                call spectrum(2**(level + 1), ka, ks, &
                    repeat**2*(6*ks**2 + 2*(3*ka**2 + ka)), &
                    b0, exact, .true., &
                    frequencies(:, level), optical(level), unmatched, &
                    perturbed=mesh == 2, structural=.true.)
                errors(:, level) = abs(frequencies(:, level)/exact - 1)
                print '(a,i4,3es20.12)', 'wave frequencies [1/s], cells per period: ', &
                    2**(level + 1), frequencies(:, level)
                print '(a,i4,3es20.12)', 'relative frequency errors: ', &
                    2**(level + 1), errors(:, level)
                if (level == 1) cycle
                rates = log(errors(:, level - 1)/errors(:, level))/log(2.0_dp)
                print '(a,3f12.6)', 'Alfven / fast / slow frequency orders: ', rates
                if (level == levels) &
                    call require(minval(rates) > 2*ka - 0.2_dp, 'frequency design order')
            end do
        end do
    end do
end program test_vrmhd_spectra
