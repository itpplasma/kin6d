program test_vrmhd_spectral_folding
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use vrmhd_structural_spectrum, only: match_folded
    use vrmhd_spectrum_support, only: require
    use rm_mhd_oracles, only: oracle_mu0
    implicit none
    real(dp) :: omega(25), exact(3), folded(3), frequencies(3), b0(3)
    real(dp) :: shares(4,25), perturbed_shares(16,25), captured(25)
    character(16) :: gate
    integer :: cell, offset
    b0 = sqrt(oracle_mu0)*[0.6_dp,0.0_dp,0.8_dp]
    exact = [0.6_dp,sqrt(1.8_dp),sqrt(0.2_dp)]
    call get_command_argument(1,gate)
    if (len_trim(gate) > 0) then
        omega = 0; shares = 0; captured = 1
        omega(1:3) = exact; omega(4:6) = -exact
        shares(4,:) = 1
        ! The alias q=(1,-16,0) has |q|^2=257 and (q.v_A)^2=0.36.
        omega(7:8) = [0.6_dp,-0.6_dp]
        omega(9) = sqrt(257.0_dp-sqrt(257.0_dp**2-0.36_dp*257.0_dp))
        omega(10) = -omega(9)
        captured(7:10) = 0.05_dp
        select case (gate)
        case ('spurious')
            captured(9:10) = 1
        case ('ambiguous')
            captured(1) = 0.5_dp
        case default
            error stop 'unknown resolved-power regression'
        end select
        call match_folded(omega,16,1,2,b0,exact,frequencies, &
            shares=shares,captured=captured)
        stop
    end if
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
    omega([2,8]) = exact(2)+1.0e-6_dp
    omega([5,11]) = -omega([2,8])
    call match_folded(omega,4,2,1,b0,exact,frequencies)
    call require(maxval(abs(frequencies-[exact(1),omega(2),exact(3)])) < 1.0e-14_dp, &
        'retained-band gap includes both displaced folded fast waves')
    ! A resolved fast wave can approach its oracle from above.
    omega = 0
    omega(1:3) = exact
    omega(2) = exact(2)+1.0e-6_dp
    omega(4:6) = -omega(1:3)
    shares = 1; captured = 1
    call match_folded(omega,16,1,1,b0,exact,frequencies, &
        shares=shares(:1,:),captured=captured)
    call require(maxval(abs(frequencies-omega(1:3))) < 1.0e-14_dp, &
        'gap cutoff includes a fast wave displaced above its oracle')
    omega = 0; shares = 0; captured = 1
    omega(1:3) = exact; omega(4:6) = -exact
    omega(7:9) = [0.6147_dp,1.1318_dp,2.6248_dp]
    omega(10:12) = -omega(7:9)
    shares(4,:) = 1
    captured(1:6) = 0.9_dp
    captured(7:12) = 0.1_dp
    call match_folded(omega,16,1,2,b0,exact,frequencies, &
        shares=shares,captured=captured)
    call require(maxval(abs(frequencies-exact)) < 1.0e-14_dp, &
        'total resolved power separates waves at both inclusive class thresholds')
    captured(1:6) = 1; captured(7:12) = 0.5_dp
    shares(1,7:12) = 0.99_dp; shares(4,7:12) = 0.01_dp
    call match_folded(omega,16,1,2,b0,exact,frequencies, &
        shares=shares,captured=captured)
    call require(maxval(abs(frequencies-exact)) < 1.0e-14_dp, &
        'resolved power includes the base-alias share and captured fraction')
    omega(9:) = 0
    perturbed_shares = 0; perturbed_shares(11,:) = 1
    perturbed_shares(11,7:8) = 0.01_dp; perturbed_shares(16,7:8) = 0.99_dp
    call match_folded(omega,16,2,2,b0,exact,frequencies, &
        shares=perturbed_shares,captured=captured)
    call require(maxval(abs(frequencies-exact)) < 1.0e-14_dp, &
        'two-by-two supercell uses its interior G=0 alias with unchanged thresholds')
    ! At four cells the Alfven dispersion error can exceed the isolation window.
    omega = 0
    omega(1:3) = [0.7_dp,1.25_dp,0.42_dp]
    omega(4:5) = -omega(1:2)
    call match_folded(omega,4,1,1,b0,exact,frequencies,resolved=.false.)
    call require(maxval(abs(frequencies-[0.7_dp,1.25_dp,0.42_dp])) < 1.0e-14_dp, &
        'coarse ladder reports unequal counts without asserting count or isolation')
end program
