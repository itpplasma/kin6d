program test_classical_pair
    use, intrinsic :: iso_fortran_env, only: dp=>real64
    use fortnum_status, only: fortnum_status_t, FORTNUM_OK
    use kin6d_classical_pair, only: classical_pair_model_t, initialize_classical_pair, &
        classical_pair_state, classical_pair_bin_center
    use kin6d_classical_pair_leaves, only: cp_statistic
    use kin6d_coulomb_leaves, only: ewald_real_pair_jvp, ewald_real_pair_vjp
    use kin6d_nbody, only: nbody_force, nbody_force_jvp, nbody_step, nbody_energy
    implicit none
    type(classical_pair_model_t) :: a, b
    type(fortnum_status_t) :: status
    real(dp) :: r(3), s(12), f(6), fm(6), df(6), u, um, up, h, weights(4), g(5)
    real(dp) :: x(12), e0, energy, dot1, dot2, tangent(4)
    integer :: j
    logical :: ok
    call initialize_classical_pair(a, ok)
    if (.not. ok) error stop 'model initialization'
    ! Independent splitting: a looser tolerance gives another alpha and cutoffs.
    b%tolerance = 1e-9_dp
    call initialize_classical_pair(b, ok)
    if (.not. ok) error stop 'second model initialization'
    if (abs(a%ewald%alpha - b%ewald%alpha) < 1e-3_dp) error stop 'distinct splitting'
    r = [0.7_dp, -1.1_dp, 0.3_dp]
    call classical_pair_state([0.4_dp, 2.0_dp, -1.0_dp], r, s)
    call nbody_force(a%system, s, f, u, status)
    if (status%code /= FORTNUM_OK) error stop 'force'
    if (maxval(abs(f(:3) + f(4:))) > 1e-12_dp) error stop 'pair antisymmetry'
    call classical_pair_state([0.4_dp, 2.0_dp, -1.0_dp], -r, s)
    call nbody_force(a%system, s, fm, um, status)
    if (maxval(abs(f(:3) + fm(:3))) > 1e-12_dp) error stop 'pair parity'
    if (abs(u - um) > 1e-12_dp) error stop 'pair energy parity'
    call classical_pair_state([3.0_dp, -0.2_dp, 5.0_dp], r, s)
    call nbody_force(b%system, s, fm, um, status)
    ! Translation invariance and splitting independence within both bounds
    ! (|g|=2 per pair; energy also carries two self terms of |g|/2 each).
    if (maxval(abs(f - fm)) > 2*(a%ewald%force_bound + b%ewald%force_bound) + &
        1e-12_dp) error stop 'independent Ewald splitting'
    if (abs(u - um) > 4*(a%ewald%potential_bound + b%ewald%potential_bound) + &
        1e-12_dp) error stop 'background energy splitting'
    h = 1e-5_dp
    call classical_pair_state([0.4_dp, 2.0_dp, -1.0_dp], r, s)
    do j = 1, 3
        x = s; x(j) = s(j) + h
        call nbody_force(a%system, x, fm, up, status)
        x(j) = s(j) - h
        call nbody_force(a%system, x, fm, um, status)
        if (abs((up - um)/(2*h) + f(j)) > 1e-8_dp) error stop 'energy gradient'
        x = 0; x(j) = 1
        call nbody_force_jvp(a%system, s(:6), x(:6), df, status)
        x = s; x(j) = s(j) + h
        call nbody_force(a%system, x, fm, up, status)
        x(j) = s(j) - h
        call nbody_force(a%system, x, f, up, status)
        if (maxval(abs(df - (fm - f)/(2*h))) > 1e-7_dp) error stop 'force JVP'
        call nbody_force(a%system, s, f, u, status)
    end do
    weights = [0.2_dp, -0.7_dp, 0.4_dp, 0.3_dp]
    call ewald_real_pair_jvp(r(1), r(2), r(3), a%ewald%alpha, 2.0_dp, &
        0.2_dp, -0.1_dp, 0.3_dp, 0.4_dp, -0.2_dp, tangent(1), tangent(2), &
        tangent(3), tangent(4))
    call ewald_real_pair_vjp(r(1), r(2), r(3), a%ewald%alpha, 2.0_dp, &
        weights(1), weights(2), weights(3), weights(4), g(1), g(2), g(3), g(4), g(5))
    dot1 = dot_product(weights, tangent)
    dot2 = dot_product(g, [0.2_dp, -0.1_dp, 0.3_dp, 0.4_dp, -0.2_dp])
    if (abs(dot1 - dot2) > 1e-12_dp) error stop 'leaf JVP/VJP pairing'
    call nbody_energy(a%system, s, e0, status)
    do j = 1, 100
        call nbody_step(a%system, s, 0.001_dp, status)
        if (status%code /= FORTNUM_OK) error stop 'step'
    end do
    call nbody_energy(a%system, s, energy, status)
    if (abs(energy - e0) > 1e-7_dp) error stop 'short-run Hamiltonian drift'
    if (abs(classical_pair_bin_center(1) + 15.0_dp/32) > 1e-14_dp) &
        error stop 'first bin center'
    if (abs(classical_pair_bin_center(16) - 15.0_dp/32) > 1e-14_dp) &
        error stop 'last bin center'
    call cp_statistic(256.0_dp, 4096.0_dp, 11.0_dp, u, um)
    if (abs(u - 1) > 1e-14_dp) error stop 'histogram normalization'
    if (abs(um - 16*sqrt(11.0_dp/(2*4096))) > 1e-14_dp) error stop 'sampling radius'
    if (2*16*61*exp(-11.0_dp) >= .05_dp) error stop 'simultaneous coverage allocation'
    print '(a)', 'PASS: Ewald splitting, parity, energy gradient, AD and Hamiltonian'
end program
