program test_ideal_seed_3d
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_ideal_seed_3d, only: toroidal_profile_t, toroidal_seed_t, &
        ideal_seed_3d_init, ideal_seed_disk_mesh, ideal_seed_boundary_coefficients
    use kin6d_ideal_map, only: ideal_map_t, ideal_map_rdof, ideal_map_energy, &
        ideal_map_init, ideal_map_freeze
    use kin6d_equilibrium_chart, only: cylindrical_chart_t
    use kin6d_grad_shafranov, only: gs_state_t
    use kin6d_ideal_seed_3d_leaves, only: seed_boundary_mode, &
        seed_boundary_mode_jvp, seed_boundary_mode_vjp
    use omp_lib, only: omp_set_num_threads
    implicit none
    type(toroidal_profile_t) :: profile
    type(toroidal_seed_t) :: seed
    type(ideal_map_t) :: map, threaded, gauge_one, gauge_four
    type(cylindrical_chart_t) :: chart
    type(gs_state_t) :: mesh
    real(dp) :: coefficients(4, -1:1, 0:3), modes(3, 3), q(3), exact(2)
    real(dp) :: point(2), grad(2), chi, b3, p, theta, phi, phase, error, e1, e4, normal(3), span, minimum_span
    real(dp) :: iota, integral, inputs(6), direction(6), bar(4), tangent(4), adj(6)
    real(dp) :: plus(4), minus(4), h
    integer :: status, j, m, n
    real(dp), parameter :: pi = acos(-1.0_dp), flux = 0.7_dp
    coefficients = 0
    coefficients(1, 0, 0) = 3
    coefficients(1, 0, 1) = 0.4_dp
    coefficients(4, 0, 1) = 0.5_dp
    ! All four coefficient channels, both signs of n, general field periods.
    coefficients(:, 1, 1) = [0.02_dp, 0.01_dp, -0.01_dp, 0.02_dp]
    coefficients(:, -1, 2) = [0.01_dp, -0.02_dp, 0.01_dp, 0.01_dp]
    coefficients(2, 0, 3) = 0.005_dp
    error = 0
    do j = 1, 71
        theta = 2*pi*j/71; phi = 2*pi*mod(13*j, 71)/71
        point = [2 + cos(theta), sin(theta)]
        call ideal_seed_boundary_coefficients(point, coefficients, modes)
        q = modes(:, 1) + modes(:, 2)*cos(3*phi) + modes(:, 3)*sin(3*phi)
        exact = 0
        do m = 0, 3
            do n = -1, 1
                phase = m*theta - 3*n*phi
                exact = exact + coefficients([1, 3], n, m)*cos(phase) &
                    + coefficients([2, 4], n, m)*sin(phase)
            end do
        end do
        error = max(error, maxval(abs(q(:2) - exact)))
    end do
    if (error > 3e-15_dp) error stop 'independent boundary reconstruction'
    ! Analytic integral of a sheared iota, including nonzero edge pressure.
    profile%iota = [0.4_dp, 0.2_dp, -0.05_dp]
    profile%pressure = [0.3_dp, -0.2_dp]
    call profile%evaluate(0.6_dp, iota, p, integral)
    if (abs(iota - (0.4_dp + 0.2_dp*0.6_dp - 0.05_dp*0.6_dp**2)) > 1e-15_dp) &
        error stop 'iota polynomial'
    if (abs(integral - (0.4_dp*0.6_dp + 0.1_dp*0.6_dp**2 &
        - 0.05_dp*0.6_dp**3/3)) > 1e-15_dp) error stop 'iota primitive'
    call ideal_seed_disk_mesh(mesh, 8, 2, status)
    if (status /= 0) error stop 'reference mesh'
    call omp_set_num_threads(1)
    call ideal_seed_3d_init(map, seed, mesh, coefficients, 3, flux, profile, status)
    if (status /= 0) error stop '3D boundary/profile seed'
    call ideal_map_energy(map, map%x, e1, status)
    if (status /= 0) error stop 'seed positive Jacobian'
    call omp_set_num_threads(4)
    call ideal_seed_3d_init(threaded, seed, mesh, coefficients, 3, flux, &
        profile, status)
    if (status /= 0) error stop 'threaded seed'
    call ideal_map_energy(threaded, threaded%x, e4, status)
    if (maxval(abs(map%x - threaded%x)) > 1e-15_dp) error stop 'thread seed parity'
    call omp_set_num_threads(1)
    call ideal_map_init(gauge_one,mesh,chart,3,1,status,nphi=7)
    if(status/=0) error stop 'serial gauge map'
    call ideal_map_freeze(gauge_one,seed,status)
    if(status/=0) error stop 'serial gauge freeze'
    call omp_set_num_threads(4)
    call ideal_map_init(gauge_four,mesh,chart,3,1,status,nphi=7)
    if(status/=0) error stop 'thread gauge map'
    call ideal_map_freeze(gauge_four,seed,status)
    if(status/=0) error stop 'thread gauge freeze'
    if (maxval(abs(gauge_one%omega_gauge-gauge_four%omega_gauge)) > &
        128*epsilon(1.0_dp)*maxval(abs(gauge_one%omega_gauge))) &
        error stop 'thread gauge assembly parity'
    if (abs(e1 - e4) > 1e-13_dp*abs(e1)) error stop 'thread energy parity'
    do j = 1, map%nodes
        if (.not. map%boundary(j)) cycle
        do m = 1, map%nmode
            if (.not. map%fixed(ideal_map_rdof(map, 2, m, j))) &
                error stop 'fixed 3D wall omega released'
        end do
    end do
    point = [2.3_dp, -0.2_dp]
    call seed%evaluate(point, chi, grad, b3, p)
    call profile%evaluate(0.13_dp, iota, p, integral)
    if (abs(pi*b3 - flux) > 1e-15_dp) error stop 'total toroidal flux'
    if (maxval(abs(grad - flux/pi*iota*[0.3_dp, -0.2_dp])) > 1e-15_dp) &
        error stop 'poloidal flux gradient'
    ! Independent finite differences and JVP/VJP physical pairings.
    inputs = [0.3_dp, -0.4_dp, 0.5_dp, 0.2_dp, -0.1_dp, 0.7_dp]
    direction = [0.1_dp, -0.2_dp, 0.3_dp, 0.4_dp, -0.5_dp, 0.6_dp]
    bar = [0.3_dp, -0.1_dp, 0.4_dp, 0.2_dp]; h = 1e-5_dp
    call seed_boundary_mode_jvp(inputs(1), inputs(2), inputs(3), inputs(4), &
        inputs(5), inputs(6), direction(1), direction(2), direction(3), &
        direction(4), direction(5), direction(6), &
        tangent(1), tangent(2), tangent(3), tangent(4))
    call seed_boundary_mode_vjp(inputs(1), inputs(2), inputs(3), inputs(4), &
        inputs(5), inputs(6), bar(1), bar(2), bar(3), bar(4), &
        adj(1), adj(2), adj(3), adj(4), adj(5), adj(6))
    if (abs(dot_product(tangent, bar) - dot_product(direction, adj)) > 1e-15_dp) &
        error stop 'boundary JVP/VJP pairing'
    call mode_at(inputs + h*direction, plus)
    call mode_at(inputs - h*direction, minus)
    if (maxval(abs((plus - minus)/(2*h) - tangent)) > 1e-10_dp) &
        error stop 'boundary JVP independent difference'
    ! Reversing Z reverses the supplied angle; it does not fold the map.
    coefficients(3:4, :, :) = -coefficients(3:4, :, :)
    call ideal_seed_3d_init(threaded, seed, mesh, coefficients, 3, flux, &
        profile, status)
    if (status /= 0) error stop 'clockwise boundary rejected'
    if (map%sigma*threaded%sigma >= 0) error stop 'boundary orientation not updated'
    minimum_span=huge(1.0_dp)
    do j=1,threaded%nodes
        if(threaded%boundary(j)) cycle
        normal=[threaded%reference(1,j)-2,-threaded%reference(2,j),0.0_dp]
        if(norm2(normal)<0.1_dp) cycle
        normal=normal/norm2(normal)
        span=dot_product(normal,threaded%basis(:,1,j))**2 &
            +dot_product(normal,threaded%basis(:,2,j))**2
        minimum_span=min(minimum_span,span)
    end do
    if(minimum_span<0.5_dp) error stop 'clockwise seed loses a physical normal direction'
    coefficients(1, 0, 2) = 0.8_dp
    call ideal_seed_3d_init(threaded, seed, mesh, coefficients, 3, flux, &
        profile, status)
    if (status == 0) error stop 'folded boundary admitted'
    coefficients=0
    coefficients(1,0,0)=-3
    coefficients(1,0,1)=0.4_dp;coefficients(4,0,1)=0.5_dp
    call ideal_seed_3d_init(threaded,seed,mesh,coefficients,3,flux,profile,status)
    if(status==0) error stop 'negative cylindrical radius admitted'
    print '(a,3es14.5)', '[V] seed boundary error, energies (1/4 threads): ', error, e1, e4
    print '(a)', '[V] flux, sheared profiles, fixed wall, derivatives and fold rejection pass'
    print '(a)', '[V] initializer candidate only; continuous bounds unavailable'
contains
    subroutine mode_at(x, y)
        real(dp), intent(in) :: x(6)
        real(dp), intent(out) :: y(4)
        call seed_boundary_mode(x(1), x(2), x(3), x(4), x(5), x(6), &
            y(1), y(2), y(3), y(4))
    end subroutine
end program
