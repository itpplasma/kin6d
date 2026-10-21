program test_equilibrium_boundary_3d
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortfem_sphere_surface_mesh, only: generate_sphere_surface_mesh
    use fortfem_parametric_bem_high_order, only: bem_pair_options_t, bem_pair_report_t
    use kin6d_equilibrium_boundary_3d, only: gamma_c_3d_t, gamma_c_3d_build, &
                                             gamma_c_3d_energy
    use fortfem_torus_surface_mesh, only: generate_torus_surface_mesh
    use fortfem_parametric_bem_high_order, only: build_p2_trace_connectivity
    use kin6d_equilibrium_boundary_3d, only: gamma_c_3d_harmonics
    use parametric_bem_test_geometry, only: test_sphere_t, test_torus_t
    implicit none
    type(test_sphere_t) :: sphere
    type(gamma_c_3d_t) :: boundary
    type(bem_pair_options_t) :: options
    type(bem_pair_report_t) :: report
    real(dp), allocatable :: r(:, :), u(:), drive(:), h(:), load(:), d(:), lp(:), lm(:)
    real(dp) :: energy, ep, em, step
    integer :: i, n, status
    character(16) :: argument
    if (command_argument_count() > 0) then
        call get_command_argument(1, argument)
        if (trim(argument) == 'harmonic') then
            call harmonic_data()
            stop
        end if
    end if
    call generate_sphere_surface_mesh(1.0_dp, 1, sphere%vertices, sphere%faces)
    ! The P2 trace has vertices and shared edges: identity restriction is a
    ! valid D6 trace adapter and the exact constant spherical energy is 4 pi.
    n = size(sphere%vertices, 2) + 3*size(sphere%faces, 2)/2
    allocate (r(n, n), u(n), drive(n), h(n), load(n), d(n), lp(n), lm(n))
    r = 0
    do i = 1, n
        r(i, i) = 1
    end do
    call gamma_c_3d_build(sphere, sphere%faces, r, options, boundary, report, status)
    if (status /= 0) error stop '3D D6 build'
    u = 1; drive = 0; h = 0
    call gamma_c_3d_energy(boundary, u, drive, h, 1.0_dp, energy, load, status)
    if (status /= 0) error stop '3D D6 action'
    if (abs(energy/(2*acos(-1.0_dp)) - 1) > 0.01_dp) error stop 'sphere energy oracle'
    d = [(sin(real(i, dp)), i=1, n)]; step = 1.e-5_dp
    call gamma_c_3d_energy(boundary, u + step*d, drive, h, 1.0_dp, ep, lp, status)
    call gamma_c_3d_energy(boundary, u - step*d, drive, h, 1.0_dp, em, lm, status)
    if (abs((ep - em)/(2*step) - dot_product(load, d)) > 1.e-8_dp) &
        error stop 'D6 energy gradient'
    if (maxval(abs((lp - lm)/(2*step) - matmul(boundary%s, d))) > 1.e-8_dp) &
        error stop 'D6 Hessian action'
    print '(a,es16.7)', '3D D6 sphere energy ', energy
contains
    subroutine harmonic_data()
        type(test_torus_t) :: torus
        type(gamma_c_3d_t) :: boundary
        type(bem_pair_options_t) :: options
        type(bem_pair_report_t) :: report
        integer, allocatable :: td(:, :), fd(:, :)
        real(dp), allocatable :: vertices(:, :), r(:, :), u(:), a(:), h(:), load(:)
        real(dp) :: wire(3, 128), theta, periods(2, 1), energy
        integer :: j, n, status
        call generate_torus_surface_mesh(2.0_dp, 0.6_dp, 4, 8, &
                                         vertices, torus%faces, torus%parameters)
        call build_p2_trace_connectivity(torus%faces, td, fd, status)
        if (status /= 0) error stop 'harmonic data topology'
        n = maxval(td); allocate (r(n, n), u(n), a(n), h(n), load(n))
        r = 0; u = 0; a = 0; h = 0
        do j = 1, n
            r(j, j) = 1
        end do
        options%maximum_degree = 80
        call gamma_c_3d_build(torus, torus%faces, r, options, boundary, report, status)
        print '(a,3i10,es18.9)', 'Harmonic-data assembly status failed degree ratio ', &
            status, report%failed_pairs, report%maximum_degree_used, report%error_ratio
        if (status /= 0) error stop 'harmonic data BEM'
        do j = 1, 128
            theta = 2*acos(-1.0_dp)*real(j - 1, dp)/128
            wire(:, j) = [2*cos(theta), 2*sin(theta), 0.0_dp]
        end do
        ! Independent Ampere linking coordinates for the oriented core circle.
        periods(:, 1) = [-1.0_dp, 0.0_dp]
       call gamma_c_3d_harmonics(boundary, torus, torus%faces, td, 16, wire, [1, 129], &
                                  periods, [2.0_dp, 0.0_dp], [1], 1.e-8_dp, status)
        if (status /= 0) error stop 'compatible harmonic data'
    if(abs(boundary%coefficients(1)+2)>1.e-12_dp) error stop 'physical period amplitude'
        call gamma_c_3d_energy(boundary, u, a, h, 1.0_dp, energy, load, status)
        if (status /= 0) error stop 'harmonic data energy'
        if (abs(energy/(2*boundary%gram(1, 1)) - 1) > 1.e-12_dp) &
            error stop 'harmonic vacuum energy'
       if (maxval(abs(load)) < 1.e-6_dp) error stop 'harmonic scalar cross load missing'
       call gamma_c_3d_harmonics(boundary, torus, torus%faces, td, 16, wire, [1, 129], &
                                  periods, [2.0_dp, 1.0_dp], [1], 1.e-8_dp, status)
        if (status /= 3) error stop 'incompatible exterior harmonic data accepted'
        print '(a,es18.9)', '3D harmonic-data energy ', energy
    end subroutine
end program
