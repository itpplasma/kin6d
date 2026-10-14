program test_gs_weak_ampere_loop
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_grad_shafranov, only: gs_state_t, initialize_gs, solve_gs, &
                                    write_gs_mesh, read_gs_mesh
    use kin6d_gs_profile_law, only: gs_profile_law_t
    use kin6d_gs_nonlinear, only: solve_gs_cubic_profiles
    use kin6d_gs_ampere_loop, only: gs_weak_boundary_current
    implicit none
    real(dp), parameter :: r0 = 3.0_dp, a = 0.7_dp, c = 0.03_dp, mu0 = 1.0_dp
    real(dp), parameter :: exact = -16*c*r0*acos(-1.0_dp)*a*a/mu0
    real(dp), parameter :: reaction = -0.5_dp
    real(dp), parameter :: rectangle_current = -8*c*r0*a*a/mu0
    type(gs_state_t) :: state, loaded
    type(gs_profile_law_t) :: law
    real(dp), allocatable :: boundary(:, :), weights(:), direction(:), original(:)
    real(dp), allocatable :: lift(:)
    real(dp), allocatable :: interior(:, :)
    real(dp) :: current, refined, source, errors(3), rate, h, value, tangent
    integer :: degree, level, n, j, status, side, k, node
    character(32) :: argument

    degree = 2
    if (command_argument_count() > 0) then
        call get_command_argument(1, argument)
        read (argument, *) degree
    end if
    law%knots = [0.0_dp, 1.0_dp]
    allocate (law%pp(4, 1), law%ff(4, 1))
    law%pp = 0; law%pp(1, 1) = -16*c/mu0; law%ff = 0
    ! Independent analytic volume integral of the exact Solov'ev source.
    ! No imposed-current constraint or exact solution is used by the diagnostic.
    do level = 1, 3
        n = 12*2**(level - 1)
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
        call gs_weak_boundary_current(state, mu0, law, 8, current, status)
        if (status /= 0) error stop 'weak current'
        call gs_weak_boundary_current(state, mu0, law, 16, refined, status, &
                                      source_current=source)
        if (status /= 0) error stop 'refined weak current'
        if (abs(current - refined) > 1.0e-10_dp*abs(exact)) &
            error stop 'weak quadrature unresolved'
        errors(level) = abs(refined/exact - 1)
        print '(a,2i6,3es20.10)', 'degree n weak_current relative_error balance ', &
            degree, n, refined, errors(level), refined - source
        if (level == 1) then
            call write_gs_mesh(state, 'weak_ampere_mesh.dat', status)
            if (status /= 0) error stop 'mesh export'
            call read_gs_mesh(loaded, 'weak_ampere_mesh.dat', status)
            if (status /= 0) error stop 'evaluation-only mesh readback'
            call gs_weak_boundary_current(loaded, mu0, law, 16, value, status)
            if (status /= 0) error stop 'exported weak current'
            if (loaded%ready .or. allocated(loaded%stiffness)) &
                error stop 'readback must remain evaluation-only'
            if (abs(value - refined) > 1.0e-12_dp*abs(exact)) &
                error stop 'weak export readback current'
            allocate (weights(size(state%psi)), direction(size(state%psi)), &
                      original(size(state%psi)), lift(size(state%psi)))
            original = state%psi
            lift = 1
            call gs_weak_boundary_current(state, mu0, law, 16, value, status, lift=lift)
            if (status /= 0) error stop 'constant-one lift'
            ! Galerkin balance: boundary residual agrees with the independent
            ! all-domain source integral, to the existing assembly accuracy.
            if (abs(value - refined) > 1.0e-6_dp*abs(exact)) &
                error stop 'weak lift dependence on solved state'
            call gs_weak_boundary_current(state, mu0, law, 16, value, status, weights)
            if (status /= 0) error stop 'weak current nodal weights'
            direction = sin(state%points(1, :)) + cos(state%points(2, :))
            tangent = dot_product(weights, direction)
            state%psi = original + 1.0e-4_dp*direction
            call gs_weak_boundary_current(state, mu0, law, 16, value, status)
            if (status /= 0) error stop 'perturbed weak current'
            if (abs((value - refined)/1.0e-4_dp - tangent) &
                > 1.0e-7_dp*abs(exact)) error stop 'fixed-geometry weak current JVP'
            if (abs(value - refined) < 1.0e-8_dp*abs(exact)) &
                error stop 'weak flux must detect off-equation gradient'
            state%psi = original + 32.0_dp
            call gs_weak_boundary_current(state, mu0, law, 16, value, status)
            if (status /= 0) error stop 'gauge shifted weak current'
            if (abs(value - refined) > 1.0e-10_dp*abs(exact)) &
                error stop 'constant flux gauge'
            state%psi = -original
            law%pp = -law%pp
            call gs_weak_boundary_current(state, mu0, law, 16, value, status)
            if (status /= 0) error stop 'reversed weak current'
            if (abs(value + refined) > 1.0e-12_dp*abs(exact)) &
                error stop 'weak current polarity'
            law%pp = -law%pp
            ! A normalized nonconstant source exercises both the mass term and
            ! the stationary-axis rank-one derivative on an evaluation-only state.
            loaded%psi = a*a - (loaded%points(1, :) - r0 - 0.023_dp)**2 &
                         - (loaded%points(2, :) - 0.017_dp)**2
            original = loaded%psi
            law%normalized = .true.
            law%axis = a*a
            law%pp(:, 1) = [-0.1_dp, 0.2_dp, 0.0_dp, 0.0_dp]
            call gs_weak_boundary_current(loaded, mu0, law, 16, current, status, &
                                          nodal_weights=weights)
            if (status /= 0) error stop 'normalized weak current weights'
            tangent = dot_product(weights, direction)
            loaded%psi = original + 1.0e-5_dp*direction
            call gs_weak_boundary_current(loaded, mu0, law, 16, value, status)
            if (status /= 0) error stop 'normalized weak positive perturbation'
            loaded%psi = original - 1.0e-5_dp*direction
            call gs_weak_boundary_current(loaded, mu0, law, 16, source, status)
            if (status /= 0) error stop 'normalized weak negative perturbation'
            if (abs((value - source)/2.0e-5_dp - tangent) &
                > 1.0e-7_dp*max(1.0_dp, abs(tangent))) &
                error stop 'normalized source and axis current derivative'
            law%normalized = .false.
            law%pp(:, 1) = [-16*c/mu0, 0.0_dp, 0.0_dp, 0.0_dp]
            deallocate (weights, direction, original, lift)
        end if
        deallocate (boundary)
    end do
    do level = 2, 3
        rate = log(errors(level - 1)/errors(level))/log(2.0_dp)
        print '(a,f10.5)', 'curved-domain weak current convergence rate ', rate
        ! Native equispaced cubic boundary nodes give area error O(h^4).
        ! Domain error is separate from the Galerkin functional error.
        if (rate < 3.5_dp) error stop 'curved-domain weak current rate'
    end do
    if (errors(3) > 1.0e-5_dp) error stop 'total weak current accuracy'
    ! psi=c*R^2*(Z^2-a^2) is Solov'ev: Delta*psi=2*c*R^2.
    ! The rectangle is exact and every Dirichlet edge trace is quadratic.
    ! Add a monotone linear source and an independent forcing cancelling it
    ! at the exact state. The total flux now depends on the numerical state;
    ! its smooth volume goal has nontrivial Galerkin O(h^(2p)) error.
    law%pp(:, 1) = [-2*c/mu0, reaction/mu0, 0.0_dp, 0.0_dp]
    do level = 1, 3
        n = 32*2**(level - 1)
        allocate (boundary(2, n))
        do j = 1, n
            call rectangle(real(j - 1, dp)/n, boundary(:, j))
        end do
        h = 8*a/n
        side = n/4
        allocate (interior(2, (side - 1)**2))
        node = 0
        do k = 1, side - 1
            do j = 1, side - 1
                node = node + 1
                interior(:, node) = [r0 - a + j*h, -a + k*h]
            end do
        end do
        call initialize_gs(state, boundary, 1.01_dp*h*h/2, status, degree, &
                           curved=.false., interior_points=interior)
        if (status /= 0) error stop 'exact rectangle mesh'
        do j = 1, size(state%psi)
            if (.not. state%boundary(j)) cycle
            state%psi(j) = rectangle_solovev(state%points(1, j), state%points(2, j))
        end do
        call solve_gs_cubic_profiles(state, mu0, law, status, &
                                     relative_tolerance=1.0e-14_dp, &
                                     quadrature_tolerance=1.0e-12_dp, forcing=forcing)
        if (status /= 0) error stop 'manufactured Solovev reaction solve'
        call gs_weak_boundary_current(state, mu0, law, 16, current, status, &
                                      tolerance=1.0e-12_dp, forcing=forcing, &
                                      source_current=source)
        if (status /= 0) error stop 'manufactured Solovev reaction flux'
        errors(level) = abs(current/rectangle_current - 1)
        print '(a,3i7,3es20.10)', 'degree n dof flux relative_error balance ', &
            degree, n, count(.not. state%boundary), current, errors(level), &
            current - source
        deallocate (boundary, interior)
    end do
    do level = 2, 3
        rate = log(errors(level - 1)/errors(level))/log(2.0_dp)
        print '(a,f10.5)', 'Galerkin weak current convergence rate ', rate
        if (rate < 2*degree - 0.5_dp) error stop 'weak current superconvergence rate'
    end do
    if (errors(3) > 1.0e-5_dp) error stop 'manufactured weak current accuracy'
contains
    pure real(dp) function rectangle_solovev(r, z) result(psi)
        real(dp), intent(in) :: r, z
        psi = c*r*r*(z*z - a*a)
    end function
    function forcing(r, z) result(source)
        real(dp), intent(in) :: r, z
        real(dp) :: source
        source = -reaction*r*rectangle_solovev(r, z)
    end function
    subroutine rectangle(t, point)
        real(dp), intent(in) :: t
        real(dp), intent(out) :: point(2)
        real(dp) :: u
        u = 4*t
        select case (int(u))
        case (0)
            point = [r0 - a + 2*a*u, -a]
        case (1)
            point = [r0 + a, -a + 2*a*(u - 1)]
        case (2)
            point = [r0 + a - 2*a*(u - 2), a]
        case default
            point = [r0 - a, a - 2*a*(u - 3)]
        end select
    end subroutine
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
