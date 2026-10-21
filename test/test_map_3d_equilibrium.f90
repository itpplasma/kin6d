program test_map_3d_equilibrium
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_ideal_seed_3d, only: toroidal_profile_t, toroidal_seed_t, &
                                   ideal_seed_3d_init, ideal_seed_disk_mesh
    use kin6d_ideal_map, only: ideal_map_t
    use kin6d_grad_shafranov, only: gs_state_t
    use kin6d_equilibrium_state, only: equilibrium_state_t, EQ_VALUES
    use kin6d_map_3d_equilibrium, only: map_3d_equilibrium, map_3d_sample, &
                                        map_3d_coordinates
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
    print '(a,es14.5)', '[V] 3D physical round-trip error: ', error
   print '(a)', '[V] candidate state and derivatives have unavailable continuous bounds'
end program
