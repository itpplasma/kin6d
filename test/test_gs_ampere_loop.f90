program test_gs_ampere_loop
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_grad_shafranov, only: gs_state_t, initialize_gs, solve_gs, &
                                    write_gs_mesh, read_gs_mesh
    use kin6d_gs_ampere_loop, only: gs_boundary_current
    implicit none
    real(dp), parameter :: r0 = 3.0_dp, a = 0.7_dp, c = 0.03_dp, mu0 = 1.0_dp
    real(dp), parameter :: exact = -16*c*r0*acos(-1.0_dp)*a*a/mu0
    type(gs_state_t) :: state, loaded
    real(dp), allocatable :: boundary(:, :), weights(:), direction(:), original(:)
    real(dp) :: current, refined, errors(3), rate, h, value, tangent
    integer :: degree, level, n, j, status
    character(32) :: argument
    degree = 2
    if (command_argument_count() > 0) then
        call get_command_argument(1, argument)
        read (argument, *) degree
    end if
    ! Exact Solov'ev psi=c*((R^2-R0^2)^2+4*R^2*Z^2), Delta*=16*c*R^2.
    ! On a disk centered at R0, integral j_phi dR dZ=-16*c*R0*pi*a^2/mu0.
    ! This volume oracle is independent of field readback and loop quadrature.
    do level = 1, 3
        n = merge(192, 48, degree == 2)*2**(level - 1)
        allocate (boundary(2, n))
        do j = 1, n
            call circle(real(j - 1, dp)/n, boundary(:, j))
        end do
        h = 2*acos(-1.0_dp)*a/n
        call initialize_gs(state, boundary, h*h/2, status, degree, curve=circle)
        if (status /= 0) error stop 'Solovev disk mesh'
        do j = 1, size(state%psi)
            if (.not. state%boundary(j)) cycle
            state%psi(j) = solovev(state%points(1, j), state%points(2, j))
        end do
        call solve_gs(state, mu0, -16*c/mu0, 0.0_dp, status)
        if (status /= 0) error stop 'Solovev solve'
        call gs_boundary_current(state, mu0, 8, current, status)
        if (status /= 0) error stop 'boundary integral'
        call gs_boundary_current(state, mu0, 16, refined, status)
        if (status /= 0) error stop 'refined boundary integral'
        if (abs(current - refined) > 1.0e-10_dp*abs(exact)) &
            error stop 'loop quadrature unresolved'
        errors(level) = abs(refined/exact - 1)
        print '(a,2i6,2es20.10)', 'degree n current relative_error ', &
            degree, n, refined, errors(level)
        if (level == 1) then
            call write_gs_mesh(state, 'ampere_mesh.dat', status)
            if (status /= 0) error stop 'mesh export'
            call read_gs_mesh(loaded, 'ampere_mesh.dat', status)
            if (status /= 0) error stop 'evaluation-only mesh readback'
            call gs_boundary_current(loaded, mu0, 16, value, status)
            if (status /= 0) error stop 'exported boundary current'
            if (abs(value - refined) > 1.0e-12_dp*abs(exact)) &
                error stop 'export readback current'
            allocate (weights(size(state%psi)), direction(size(state%psi)), &
                      original(size(state%psi)))
            original = state%psi
            call gs_boundary_current(state, mu0, 16, value, status, weights)
            if (status /= 0) error stop 'current nodal weights'
            if (abs(dot_product(weights, original) - value) &
                > 1.0e-11_dp*abs(exact)) error stop 'linear current weights'
            direction = sin(state%points(1, :)) + cos(state%points(2, :))
            tangent = dot_product(weights, direction)
            state%psi = original + 1.0e-4_dp*direction
            call gs_boundary_current(state, mu0, 16, value, status)
            if (status /= 0) error stop 'perturbed current'
            if (abs((value - refined)/1.0e-4_dp - tangent) &
                > 1.0e-7_dp*abs(exact)) error stop 'fixed-geometry current JVP'
            state%psi = original + 32.0_dp
            call gs_boundary_current(state, mu0, 16, value, status)
            if (status /= 0) error stop 'gauge shifted current'
            if (abs(value - refined) > 1.0e-10_dp*abs(exact)) &
                error stop 'constant flux gauge'
            state%psi = -original
            call gs_boundary_current(state, mu0, 16, value, status)
            if (status /= 0) error stop 'reversed current'
            if (abs(value + refined) > 1.0e-12_dp*abs(exact)) &
                error stop 'current polarity'
            deallocate (weights, direction, original)
        end if
        deallocate (boundary)
    end do
    do level = 2, 3
        rate = log(errors(level - 1)/errors(level))/log(2.0_dp)
        print '(a,f10.5)', 'current convergence rate ', rate
        if (rate < degree - 0.5_dp) error stop 'total current convergence rate'
    end do
    if (errors(3) > 1.0e-5_dp) error stop 'total current accuracy'
contains
    pure real(dp) function solovev(r, z)
        real(dp), intent(in) :: r, z
        solovev = c*((r*r - r0*r0)**2 + 4*r*r*z*z)
    end function
    subroutine circle(t, point)
        real(dp), intent(in) :: t
        real(dp), intent(out) :: point(2)
        point = [r0 + a*cos(2*acos(-1.0_dp)*t), a*sin(2*acos(-1.0_dp)*t)]
    end subroutine
end program
