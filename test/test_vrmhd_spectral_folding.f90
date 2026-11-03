program test_vrmhd_spectral_folding
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use vrmhd_structural_spectrum, only: match_folded
    use vrmhd_spectrum_support, only: require
    use rm_mhd_oracles, only: oracle_mu0
    implicit none
    real(dp) :: omega(25), exact(3), folded(3), frequencies(3), b0(3)
    integer :: cell, offset
    b0 = sqrt(oracle_mu0)*[0.6_dp,0.0_dp,0.8_dp]
    exact = [0.6_dp,sqrt(1.8_dp),sqrt(0.2_dp)]
    ! q=(+/-1,-2,0) has |q|^2=5 and (q.v_A)^2=0.36.
    folded = [0.6_dp,sqrt(5.0_dp+sqrt(23.2_dp)),sqrt(5.0_dp-sqrt(23.2_dp))]
    omega = 0
    do cell = 1, 4
        offset = 6*(cell-1)
        if (cell <= 2) then
            omega(offset+1:offset+3) = exact
        else
            omega(offset+1:offset+3) = folded
        end if
        omega(offset+4:offset+6) = -omega(offset+1:offset+3)
    end do
    call match_folded(omega,4,2,1,b0,exact,frequencies)
    call require(maxval(abs(frequencies-exact)) < 1.0e-14_dp, &
        'four-cell folded multiplicities retain all three base branches')
    ! At four cells the Alfven dispersion error can exceed the isolation window.
    omega = 0
    omega(1:3) = [0.7_dp,1.25_dp,0.42_dp]
    omega(4:6) = -omega(1:3)
    call match_folded(omega,4,1,1,b0,exact,frequencies,resolved=.false.)
    call require(maxval(abs(frequencies-[0.7_dp,1.25_dp,0.42_dp])) < 1.0e-14_dp, &
        'coarse ladder retains counts and frequency estimates without isolation')
end program
