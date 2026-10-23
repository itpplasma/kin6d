module test_gs_axis_refinement_callbacks
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_grad_shafranov, only: gs_state_t, initialize_gs, solve_gs, &
        gs_axis_candidate, gs_axis_safety_factor, sample_gs
    use kin6d_gs_axis_recovery, only: gs_magnetic_axis
    use kin6d_gs_profile_law, only: gs_profile_law_t
    use kin6d_gs_nonlinear, only: solve_gs_cubic_profiles
    use kin6d_gs_readback, only: gs_readback_t, read_gs_equilibrium

    use kin6d_curved_triangle, only: boundary_curve_t
    use kin6d_grad_shafranov, only: initialize_gs_object
    implicit none
    private
    public :: run_test_gs_axis_refinement
    real(dp), parameter :: r0 = 6.2_dp, a = r0/3, c = 32.86_dp/(12*r0**3)
    real(dp), parameter :: rc = 62.0_dp, ac = 0.62_dp
    real(dp), parameter :: mu0 = 1.2566370614359173e-6_dp
    real(dp), parameter :: two_pi = 2*acos(-1.0_dp)
    integer, parameter :: nshift = 3
    real(dp), parameter :: scatter = 3.0_dp  ! measured realisation scatter 2.4
    type :: run_test_gs_axis_refinement_context_t
        real(dp), pointer :: phase => null()
    end type

    type, extends(boundary_curve_t) :: run_test_gs_axis_refinement_phased_boundary_t
        type(run_test_gs_axis_refinement_context_t), pointer :: data => null()
    contains
        procedure :: evaluate => run_test_gs_axis_refinement_phased_boundary_binding
    end type

contains
    subroutine run_test_gs_axis_refinement()
    implicit none
    real(dp) :: psi_coarse, psi_medium, psi_fine, shift(2)
    type(gs_state_t) :: state
    real(dp), allocatable :: contour(:, :)
    real(dp) :: axis(2), axis_coarse(2), axis_medium(2), axis_fine(2)
    real(dp) :: error(3), rate, h, q_coarse, q_medium, q_fine
    real(dp) , target :: shifted(nshift), phase, worst(3)
    integer :: mesh, n, i, k, status
        type(run_test_gs_axis_refinement_context_t), target :: ctx
        type(run_test_gs_axis_refinement_phased_boundary_t), target :: &
            phased_boundary_operator
    call represented_cubic_axis()
    ! The axis error of an unstructured mesh is a signed number that depends on the
    ! mesh realisation, not on the solver: the solution is unchanged when the CG
    ! tolerance is varied from 1e-8 to 1e-13, but the error at n = 384 is 0.8e-8 to
    ! 1.1e-8 (default flags) and 1.7e-8 to 4.1e-8 (-march=x86-64-v3), whose
    ! triangulation differs, and varies by a factor 2.4 between boundary-node phases
    ! of one flavour. Two consecutive errors therefore carry no rate. Each
    ! resolution is run at `nshift` node phases; the gate is that no realisation
    ! exceeds `scatter` times the h**2 extrapolation of the coarse rms error.
    do mesh = 1, 3
        n = 96*2**(mesh - 1)
        allocate (contour(2, n))
        h = two_pi*a/n
        error(mesh) = 0
        do k = 0, nshift - 1
            phase = real(k, dp)/(nshift*n)
            do i = 1, n
        call bind_run_test_gs_axis_refinement_context()
                call run_test_gs_axis_refinement_phased_boundary(ctx, real(i - 1, &
                    dp)/n, contour(:, i))
            end do
        call bind_run_test_gs_axis_refinement_context()
            call initialize_gs_object(state, contour, h*h/2, status, 2, &
                               curve=phased_boundary_operator)
            if (status /= 0) error stop 'axis mesh'
            call solve_gs(state, 1.0_dp, -16*c, 0.0_dp, status)
            if (status /= 0) error stop 'axis solve'
            call gs_magnetic_axis(state, axis, status)
            if (status /= 0) error stop 'axis candidate'
            shifted(k + 1) = norm2(axis - [r0, 0.0_dp])/r0
            error(mesh) = error(mesh) + shifted(k + 1)**2/nshift
        end do
        error(mesh) = sqrt(error(mesh))
        worst(mesh) = maxval(shifted)
        print *, n, error(mesh), shifted
        deallocate (contour)
    end do
    do mesh = 2, 3
        rate = log(error(mesh - 1)/error(mesh))/log(2.0_dp)
        print *, 'axis rate of the rms error', rate
        if (worst(mesh) > scatter*error(1)*(96.0_dp/(96*2**(mesh - 1)))**2) &
            error stop 'axis error above the second-order extrapolation'
    end do

    ! Cylinder at aspect 100, the PLAN Phase 5 Lundquist family. On the coarse
    ! 16-node curved P3 mesh the flux minimum sits deep inside a cell corner, so
    ! the first undamped axis step overshot that cell, no cell then reported a
    ! minimum, and the readback failed with gs_axis_candidate status 4. The axis
    ! must be recovered, must minimise the represented flux, and must converge.
    call cylinder(two_pi*ac/16, axis_coarse, psi_coarse, q_coarse)
    call cylinder(two_pi*ac/32, axis_medium, psi_medium, q_medium)
    call cylinder(two_pi*ac/64, axis_fine, psi_fine, q_fine)
    print '(a,2(es15.8,1x),es18.10,es12.4)', 'coarse P3 axis', axis_coarse, &
        psi_coarse, q_coarse
    ! Registered coarse run cyl20261009_kinp3_profiles_kin6d_lundquist_A100_00,
    ! its flux minimum confirmed by an independent monomial fit of the exported
    ! P3 field to -6.6303029203853043e+01.
    if (abs(psi_coarse + 6.6303029203853043e+01_dp) > 1.0e-6_dp) then
        error stop 'coarse axis flux moved'
    end if
    if (norm2(axis_fine - [rc, 0.0_dp]) > 0.01_dp*ac) then
        error stop 'fine axis off centre'
    end if
    shift = [abs(psi_coarse - psi_medium), norm2(axis_coarse - axis_medium)]
    if (shift(1) <= abs(psi_medium - psi_fine)) error stop 'flux not refining'
    if (shift(2) <= norm2(axis_medium - axis_fine)) error stop 'axis not refining'
    if (q_coarse*q_medium <= 0 .or. q_fine*q_medium <= 0) then
        error stop 'axis q sign changed across refinement'
    end if

    ! Coarse curved P2 cylinder, aspect 100. Here the represented flux minimum
    ! sits exactly on the axis node shared by the surrounding cells. Degree 2
    ! used to answer status 4: it folded every nodal value into one "edge is
    ! lower" rejection, so a minimum on a shared entity never became an axis.
    call coarse_axis_degree2
    ! One manufactured quadratic whose unique minimum sits strictly inside an
    ! edge that two cells own jointly: one axis, reported once, and a minimum
    ! that stays on the wall is still no axis at all.
    call shared_entity_axis
contains
    subroutine bind_run_test_gs_axis_refinement_context()
        ctx%phase => phase
        phased_boundary_operator%data => ctx
    end subroutine bind_run_test_gs_axis_refinement_context

    subroutine represented_cubic_axis()
        type(gs_state_t) :: field
        real(dp) :: points(2, 48), exact(2), found(2), flux, hessian(2, 2)
        real(dp) :: x, y, psi, br, bz, bt
        integer :: j, stat, element
        exact = [3.073_dp, 0.017_dp]
        do j = 1, 48
            points(:, j) = [3 + cos(two_pi*(j - 1)/48), sin(two_pi*(j - 1)/48)]
        end do
        call initialize_gs(field, points, 0.01_dp, stat, 3)
        if (stat /= 0) error stop 'represented cubic mesh'
        do j = 1, size(field%psi)
            x = field%points(1, j) - exact(1)
            y = field%points(2, j) - exact(2)
            field%psi(j) = -1 + x*x + 2*y*y + 0.1_dp*x*x*x
            ! Perturb only distant nodal coefficients. The cubic containing
            ! the axis still has the exact displaced minimum and Hessian diag(2,4).
            if (x > 0.25_dp .and. x < 0.6_dp .and. abs(y) < 0.25_dp) &
                field%psi(j) = field%psi(j) + 0.01_dp
        end do
        call gs_axis_candidate(field, found, flux, hessian, element, stat)
        if (stat /= 0) error stop 'represented cubic minimum'
        if (abs(flux + 1) > 1.0e-12_dp) error stop 'cubic minimum value'
        if (maxval(abs(hessian - reshape([2.0_dp, 0.0_dp, 0.0_dp, 4.0_dp], &
            [2, 2]))) > 1.0e-10_dp) error stop 'cubic minimum Hessian'
        call gs_magnetic_axis(field, found, stat)
        if (stat /= 0) error stop 'represented cubic axis'
        if (norm2(found - exact) > 1.0e-12_dp) error stop 'cubic axis moved'
        call sample_gs(field, found(1), found(2), 1.0_dp, psi, br, bz, bt, stat)
        if (stat /= 0) error stop 'cubic axis sample'
        if (norm2([br, bz]) > 1.0e-12_dp) error stop 'axis not stationary'
    end subroutine

    pure subroutine cylinder_boundary(t, point)
        real(dp), intent(in) :: t
        real(dp), intent(out) :: point(2)
        point = [rc + ac*cos(two_pi*t), ac*sin(two_pi*t)]
    end subroutine

    ! One resolution of the Lundquist cylinder: solve the registered cubic source,
    ! read the axis back, and test it against the flux it must minimise.
    subroutine cylinder(edge, axis_out, psi_axis, q_axis)
        real(dp), intent(in) :: edge
        real(dp), intent(out) :: axis_out(2), psi_axis, q_axis
        type(gs_profile_law_t) :: law
        type(gs_readback_t) :: readback
        real(dp), allocatable :: polygon(:, :)
        real(dp) :: hessian(2, 2), psi, br, bz, bt, radius, angle, psi_min
        integer :: nodes, cell, step, turn, sample
        nodes = int(nint(two_pi*ac/edge))
        allocate (polygon(2, nodes))
        do sample = 1, nodes
            call cylinder_boundary(real(sample - 1, dp)/nodes, polygon(:, sample))
        end do
        call initialize_gs(state, polygon, edge*edge/2, status, 3, &
            curve=cylinder_boundary)
        if (status /= 0) error stop 'cylinder mesh'
        call lundquist_law(law)
        if (.not. law%valid()) error stop 'cylinder source law'
        call solve_gs_cubic_profiles(state, mu0, law, status, 1.0e-11_dp, &
            quadrature_tolerance=1.0e-11_dp)
        if (status /= 0) error stop 'cylinder solve'
        call gs_axis_candidate(state, axis_out, psi_axis, hessian, cell, status)
        if (status /= 0) error stop 'cylinder lost its axis'
        if (psi_axis >= 0) error stop 'cylinder axis flux sign'
        if (norm2(axis_out - [rc, 0.0_dp]) > 0.02_dp*ac) then
            error stop 'axis off centre'
        end if
        call gs_axis_safety_factor(axis_out(1), law%f_edge, hessian, q_axis, &
            status)
        if (status /= 0) error stop 'cylinder axis q'
        call read_gs_equilibrium(state, law%f_edge, 0.0_dp, 0, 0, readback, &
            status, law)
        if (status /= 0) error stop 'cylinder readback failed'
        psi_min = huge(1.0_dp)
        do step = 0, 48
            do turn = 0, 35
                radius = ac*sqrt(real(step, dp)/48)
                angle = two_pi*real(turn, dp)/36
                call sample_gs(state, rc + radius*cos(angle), radius*sin(angle), &
                    law%f_edge, psi, br, bz, bt, sample)
                if (sample /= 0) cycle
                psi_min = min(psi_min, psi)
            end do
        end do
        if (psi_min < psi_axis - 1.0e-9_dp*abs(psi_axis)) then
            error stop 'axis is not the flux minimum'
        end if
        deallocate (polygon)
    end subroutine

    ! Registered Lundquist A=100 source: p'=0, FF' linear in psi, edge F.
    subroutine lundquist_law(law)
        type(gs_profile_law_t), intent(out) :: law
        law%normalized = .false.
        law%p_edge = 0
        law%f_edge = 1.68186572932422706e+02_dp
        law%axis = -6.63042165212652890e+01_dp
        allocate (law%knots(2), law%pp(4, 1), law%ff(4, 1))
        law%knots = [-9.94563247818979335e+01_dp, 0.0_dp]
        law%pp = 0
        law%ff(:, 1) = [-9.89048500484972578e+02_dp, 5.85327783558792980_dp, &
            0.0_dp, 0.0_dp]
        if (.not. law%valid()) error stop 'cylinder source law'
    end subroutine lundquist_law

    ! Degree 2 on the coarse 16-node cylinder. The axis is a mesh node, so the
    ! reported flux is that node's flux bit for bit, and the axis cannot sit
    ! above the flux it is supposed to minimise.
    subroutine coarse_axis_degree2
        type(gs_profile_law_t) :: law
        real(dp), allocatable :: polygon(:, :)
        real(dp) :: hessian(2, 2), axis_p2(2), psi_p2, q_p2
        integer :: nodes, sample, cell, node, status

        nodes = 16
        allocate (polygon(2, nodes))
        do sample = 1, nodes
            call cylinder_boundary(real(sample - 1, dp)/nodes, polygon(:, sample))
        end do
        call initialize_gs(state, polygon, (two_pi*ac/nodes)**2/2, status, 2, &
            curve=cylinder_boundary)
        if (status /= 0) error stop 'coarse P2 mesh'
        call lundquist_law(law)
        call solve_gs_cubic_profiles(state, mu0, law, status, 1.0e-11_dp, &
            quadrature_tolerance=1.0e-11_dp)
        if (status /= 0) error stop 'coarse P2 solve'
        call gs_axis_candidate(state, axis_p2, psi_p2, hessian, cell, status)
        print '(a,2(es15.8,1x),es18.10,1x,i4)', 'coarse P2 axis', axis_p2, &
            psi_p2, cell
        if (status /= 0) error stop 'coarse P2 rejected a shared-node axis'
        node = minloc(state%psi, 1)
        if (state%boundary(node)) error stop 'coarse P2 axis node is on the wall'
        if (psi_p2 /= minval(state%psi)) then
            error stop 'coarse P2 axis sits above the flux minimum'
        end if
        if (.not. all(axis_p2 == state%points(:, node))) then
            error stop 'coarse P2 axis moved off its node'
        end if
        if (norm2(axis_p2 - [rc, 0.0_dp]) > 1.0e-9_dp*ac) then
            error stop 'coarse P2 axis off centre'
        end if
        call gs_axis_safety_factor(axis_p2(1), law%f_edge, hessian, q_p2, status)
        if (status /= 0 .or. q_p2 <= 0) error stop 'coarse P2 axis q not an O-point'
        deallocate (polygon)
    end subroutine coarse_axis_degree2

    ! psi = ((R-r*)^2 + Z^2)/2 - 1 on the same mesh reproduces itself inside
    ! every P2 cell, so its minimum -1 at r* is known exactly. With r* strictly
    ! inside an interior edge the two adjacent cells see one entity, and the
    ! search must report it once rather than as a tie (status 5). With r* on a
    ! wall node the flux keeps falling to the boundary and there is no axis.
    subroutine shared_entity_axis
        real(dp), allocatable :: polygon(:, :)
        real(dp) :: hessian(2, 2), axis_q(2), psi_q, centre(2), q_q, h
        integer :: nodes, sample, cell, edge, j, k, node, status

        nodes = 16
        allocate (polygon(2, nodes))
        do sample = 1, nodes
            call cylinder_boundary(real(sample - 1, dp)/nodes, polygon(:, sample))
        end do
        h = (two_pi*ac/nodes)**2/2
        call initialize_gs(state, polygon, h, status, 2, curve=cylinder_boundary)
        if (status /= 0) error stop 'quadratic mesh'
        deallocate (polygon)
        centre = 0
        edge = 0
        do cell = 1, size(state%triangles, 2)
            do j = 1, 3
                k = mod(j, 3) + 1
                if (state%boundary(state%triangles(j + 3, cell))) cycle
                centre = state%points(:, state%triangles(j, cell)) + 0.25_dp* &
                    (state%points(:, state%triangles(k, cell)) &
                    - state%points(:, state%triangles(j, cell)))
                edge = 1
                exit
            end do
            if (edge == 1) exit
        end do
        if (edge == 0) error stop 'no interior edge in the coarse mesh'
        do node = 1, size(state%psi)
            state%psi(node) = ((state%points(1, node) - centre(1))**2 &
                + (state%points(2, node) - centre(2))**2)/2 - 1
        end do
        call gs_axis_candidate(state, axis_q, psi_q, hessian, cell, status)
        print '(a,2(es17.10,1x),es20.12)', 'shared edge axis', axis_q, psi_q
        if (status == 5) error stop 'shared edge axis counted as a tie'
        if (status /= 0) error stop 'interior edge minimum was not an axis'
        if (psi_q /= -1.0_dp) error stop 'edge axis flux off the exact minimum'
        if (psi_q > minval(state%psi)) then
            error stop 'edge axis sits above the flux minimum'
        end if
        if (norm2(axis_q - centre) > 1.0e-8_dp*sqrt(2*h)) then
            error stop 'edge axis off the exact minimum location'
        end if
        call gs_axis_safety_factor(axis_q(1), 1.0_dp, hessian, q_q, status)
        if (status /= 0 .or. q_q <= 0) error stop 'edge axis q not an O-point'
        node = 1
        do while (.not. state%boundary(node))
            node = node + 1
        end do
        centre = state%points(:, node)
        do j = 1, size(state%psi)
            state%psi(j) = ((state%points(1, j) - centre(1))**2 &
                + (state%points(2, j) - centre(2))**2)/2 - 1
        end do
        call gs_axis_candidate(state, axis_q, psi_q, hessian, cell, status)
        if (status /= 4) error stop 'wall-centred flux was called an enclosed axis'
    end subroutine shared_entity_axis

    end subroutine run_test_gs_axis_refinement

    subroutine run_test_gs_axis_refinement_phased_boundary(ctx, t, point)
        type(run_test_gs_axis_refinement_context_t), intent(in) :: ctx
        real(dp), intent(in) :: t
        real(dp), intent(out) :: point(2)
        call run_test_gs_axis_refinement_boundary(t + ctx%phase, point)
    end subroutine run_test_gs_axis_refinement_phased_boundary

    subroutine run_test_gs_axis_refinement_phased_boundary_binding(self, t, point)
        class(run_test_gs_axis_refinement_phased_boundary_t), intent(in) :: self
        real(dp), intent(in) :: t
        real(dp), intent(out) :: point(2)
        call run_test_gs_axis_refinement_phased_boundary(self%data, t, point)
    end subroutine

    subroutine run_test_gs_axis_refinement_boundary(t, point)
        real(dp), intent(in) :: t
        real(dp), intent(out) :: point(2)
        point(1) = sqrt(r0*r0 + 2*r0*a*cos(two_pi*t))
        point(2) = r0*a*sin(two_pi*t)/point(1)
    end subroutine run_test_gs_axis_refinement_boundary
end module test_gs_axis_refinement_callbacks

program test_gs_axis_refinement
    use test_gs_axis_refinement_callbacks, only: run_test_gs_axis_refinement
    implicit none
    call run_test_gs_axis_refinement()
end program test_gs_axis_refinement
