program test_equilibrium_boundary_3d
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortfem_sphere_surface_mesh, only: generate_sphere_surface_mesh
    use fortfem_parametric_bem_high_order, only: bem_pair_options_t, bem_pair_report_t
    use kin6d_equilibrium_boundary_3d, only: gamma_c_3d_t, gamma_c_3d_build, &
                                             gamma_c_3d_energy
    use parametric_bem_test_geometry, only: test_sphere_t
    implicit none
    type(test_sphere_t) :: sphere
    type(gamma_c_3d_t) :: boundary
    type(bem_pair_options_t) :: options
    type(bem_pair_report_t) :: report
    real(dp), allocatable :: r(:, :), u(:), drive(:), h(:), load(:), d(:), lp(:), lm(:)
    real(dp) :: energy, ep, em, step
    integer :: i, n, status
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
end program
