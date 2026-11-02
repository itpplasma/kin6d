program test_map_3d_equilibrium
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_ideal_seed_3d, only: toroidal_profile_t, toroidal_seed_t, &
                                   ideal_seed_3d_init, ideal_seed_disk_mesh
    use kin6d_ideal_map, only: ideal_map_t
    use kin6d_grad_shafranov, only: gs_state_t
    use kin6d_equilibrium_state, only: equilibrium_state_t, EQ_VALUES
    use kin6d_map_3d_equilibrium, only: map_3d_equilibrium, map_3d_sample, &
                                        map_3d_coordinates
    use kin6d_map_3d_readback_leaves, only: map_3d_field, &
        map_3d_field_jvp, map_3d_field_vjp
    use omp_lib, only: omp_set_num_threads
    implicit none
    type(ideal_map_t) :: map
    type(gs_state_t) :: mesh
    type(toroidal_profile_t) :: profile
    type(toroidal_seed_t) :: seed, wrong
    type(equilibrium_state_t) :: equilibrium
   real(dp) :: c(4,-1:1,0:1), points(3,31), values(EQ_VALUES,31), expected(EQ_VALUES,31)
    real(dp) :: parallel_values(EQ_VALUES, 31), xi(2), q(3), d(3, 3), g(2, 2), ref(2)
    real(dp) :: phi, cg(2), chi, b3, pressure, error
    real(dp) :: adinput(18), direction_ad(18), bar(7), tangent(7), adjoint(18)
    real(dp) :: yplus(7), yminus(7), h
    integer :: status, statuses(31), parallel_statuses(31), j, cell
    c = 0; c(1, 0, 0) = 2; c(1, 0, 1) = 1; c(4, 0, 1) = 1
    profile%iota = [0.4_dp, 0.2_dp]; profile%pressure = [0.3_dp, -0.2_dp]
    call ideal_seed_disk_mesh(mesh, 8, 2, status)
    if (status /= 0) error stop 'disk mesh'
    call ideal_seed_3d_init(map, seed, mesh, c, 3, 0.7_dp, profile, status)
    if (status /= 0) error stop '3D readback seed'
    ! Independent cylindrical identity-map formula.
    do j = 1, 31
        cell = 1 + mod(j - 1, map%cells); xi = [0.21_dp, 0.37_dp]; phi = 0.17_dp*j
      call map_3d_sample(map, seed, cell, xi, phi, points(:, j), expected(:, j), status)
        if (status /= 0) error stop 'identity map sample'
        call map_3d_coordinates(map, cell, xi, phi, q, d, ref, g)
        call seed%evaluate(ref, chi, cg, b3, pressure)
        if (maxval(abs(expected(:3, j) - [-cg(2)/q(1), cg(1)/q(1), b3])) > 1e-14_dp) &
            error stop 'independent cylindrical Piola'
    end do
    ! Interior nonzero Fourier modes, including a toroidal displacement.
    do j = 1, map%nodes
        map%x(1, 2, j) = 0.02_dp*map%axis_t(j)
        map%x(2, 3, j) = -0.01_dp*map%axis_t(j)
        map%x(3, 3, j) = 0.03_dp*map%axis_t(j)
    end do
    do j = 1, 31
        cell = 1 + mod(j - 1, map%cells); xi = [0.21_dp, 0.37_dp]; phi = 0.17_dp*j
      call map_3d_sample(map, seed, cell, xi, phi, points(:, j), expected(:, j), status)
        if (status /= 0) error stop 'nonzero Fourier map sample'
    end do
    wrong = seed; wrong%flux = 1.01_dp*seed%flux
    call map_3d_equilibrium(map, wrong, equilibrium, status)
    if (status == 0) error stop 'mismatched frozen form admitted'
    call map_3d_equilibrium(map, seed, equilibrium, status)
    if (status /= 0) error stop '3D core state'
    call omp_set_num_threads(1)
    call equilibrium%field%sample(points, values, statuses)
    call omp_set_num_threads(4)
    call equilibrium%field%sample(points, parallel_values, parallel_statuses)
if (any(statuses /= 0) .or. any(parallel_statuses /= 0)) error stop '3D point inversion'
    error = maxval(abs(values - expected))
    if (error > 2e-11_dp) error stop 'physical 3D map round trip'
 if (maxval(abs(values - parallel_values)) > 2e-11_dp) error stop 'thread sample parity'
    adinput=[1.2_dp,0.1_dp,0.02_dp,0.05_dp,0.9_dp,-0.01_dp, &
        0.03_dp,0.04_dp,1.01_dp,0.4_dp,0.02_dp,0.03_dp,0.3_dp, &
        0.2_dp,-0.3_dp,0.7_dp,2.8_dp,-1.0_dp]
    direction_ad=[(0.02_dp*j,j=1,18)];direction_ad(18)=0
    bar=[(0.03_dp*j,j=1,7)];h=1e-5_dp
    call map_3d_field_jvp(adinput(1), adinput(2), adinput(3), adinput(4), &
        adinput(5), adinput(6), adinput(7), adinput(8), adinput(9), adinput(10), &
        adinput(11), adinput(12), adinput(13), adinput(14), adinput(15), &
        adinput(16), adinput(17), adinput(18), direction_ad(1), direction_ad(2), &
        direction_ad(3), direction_ad(4), direction_ad(5), direction_ad(6), &
        direction_ad(7), direction_ad(8), direction_ad(9), direction_ad(10), &
        direction_ad(11), direction_ad(12), direction_ad(13), direction_ad(14), &
        direction_ad(15), direction_ad(16), direction_ad(17), direction_ad(18), &
        tangent(1), tangent(2), tangent(3), tangent(4), tangent(5), tangent(6), &
        tangent(7))
    call map_3d_field_vjp(adinput(1), adinput(2), adinput(3), adinput(4), &
        adinput(5), adinput(6), adinput(7), adinput(8), adinput(9), adinput(10), &
        adinput(11), adinput(12), adinput(13), adinput(14), adinput(15), &
        adinput(16), adinput(17), adinput(18), bar(1), bar(2), bar(3), bar(4), &
        bar(5), bar(6), bar(7), adjoint(1), adjoint(2), adjoint(3), adjoint(4), &
        adjoint(5), adjoint(6), adjoint(7), adjoint(8), adjoint(9), adjoint(10), &
        adjoint(11), adjoint(12), adjoint(13), adjoint(14), adjoint(15), &
        adjoint(16), adjoint(17), adjoint(18))
    if(abs(dot_product(tangent,bar)-dot_product(direction_ad,adjoint))>1e-12_dp) &
        error stop '3D field JVP/VJP pairing'
    call field_at(adinput+h*direction_ad,yplus)
    call field_at(adinput-h*direction_ad,yminus)
    if(maxval(abs((yplus-yminus)/(2*h)-tangent))>1e-8_dp) &
        error stop '3D field JVP finite difference'
    print '(a,es14.5)', '[V] 3D physical round-trip error: ', error
   print '(a)', '[V] candidate state and derivatives have unavailable continuous bounds'
contains
    subroutine field_at(a,y)
        real(dp),intent(in)::a(18)
        real(dp),intent(out)::y(7)
    call map_3d_field(a(1), a(2), a(3), a(4), a(5), a(6), a(7), a(8), a(9), a(10), &
        a(11), a(12), a(13), a(14), a(15), a(16), a(17), a(18), y(1), y(2), y(3), &
        y(4), y(5), y(6), y(7))
    end subroutine
end program
