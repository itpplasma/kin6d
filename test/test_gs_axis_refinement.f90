program test_gs_axis_refinement
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_grad_shafranov, only: gs_state_t, initialize_gs, solve_gs, &
        gs_axis_candidate, gs_axis_safety_factor, sample_gs
    use kin6d_gs_axis_recovery, only: gs_magnetic_axis
    use kin6d_gs_profile_law, only: gs_profile_law_t
    use kin6d_gs_nonlinear, only: solve_gs_cubic_profiles
    use kin6d_gs_readback, only: gs_readback_t, read_gs_equilibrium
    implicit none
    real(dp), parameter :: r0 = 6.2_dp, a = r0/3, c = 32.86_dp/(12*r0**3)
    real(dp), parameter :: rc = 62.0_dp, ac = 0.62_dp
    real(dp), parameter :: mu0 = 1.2566370614359173e-6_dp
    real(dp), parameter :: two_pi = 2*acos(-1.0_dp)
    real(dp) :: psi_coarse, psi_medium, psi_fine, shift(2)
    type(gs_state_t) :: state
    real(dp), allocatable :: contour(:, :)
    real(dp) :: axis(2), axis_coarse(2), axis_medium(2), axis_fine(2)
    real(dp) :: error(3), rate, h, q_coarse, q_medium, q_fine
    integer :: mesh, n, i, status
    do mesh = 1, 3
        n = 96*2**(mesh - 1)
        allocate (contour(2, n))
        do i = 1, n
            call boundary(real(i - 1, dp)/n, contour(:, i))
        end do
        h = two_pi*a/n
        call initialize_gs(state, contour, h*h/2, status, 2, curve=boundary)
        if (status /= 0) error stop 'axis mesh'
        call solve_gs(state, 1.0_dp, -16*c, 0.0_dp, status)
        if (status /= 0) error stop 'axis solve'
        call gs_magnetic_axis(state, axis, status)
        if (status /= 0) error stop 'axis candidate'
        error(mesh) = norm2(axis - [r0, 0.0_dp])/r0
        print *, n, error(mesh)
        deallocate (contour)
    end do
    do mesh = 2, 3
        rate = log(error(mesh - 1)/error(mesh))/log(2.0_dp)
        print *, 'axis rate', rate
        if (rate < 1.7_dp) error stop 'axis below field order'
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
contains
    subroutine boundary(t, point)
        real(dp), intent(in) :: t
        real(dp), intent(out) :: point(2)
        point(1) = sqrt(r0*r0 + 2*r0*a*cos(two_pi*t))
        point(2) = r0*a*sin(two_pi*t)/point(1)
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
        ! Registered Lundquist A=100 source: p'=0, FF' linear in psi, edge F.
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
end program
