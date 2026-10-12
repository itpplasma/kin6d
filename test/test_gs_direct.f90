program test_gs_direct
    ! The CSR operator, the cached Cholesky factor and the A^-1-preconditioned
    ! Newton solves agree with the element-wise operator and the Jacobi-CG
    ! solves, and the Krylov work per Newton solve does not grow with the mesh.
    ! Manufactured P3 problem of test_gs_nonlinear on a circle-like boundary.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_grad_shafranov, only: gs_state_t, initialize_gs, apply_gs
    use kin6d_gs_profile_law, only: gs_profile_law_t
    use kin6d_gs_nonlinear, only: solve_gs_cubic_profiles
    use kin6d_gs_options, only: gs_set_options
    use kin6d_gs_direct, only: gs_direct_available
    use kin6d_gs_timers, only: gs_timer_reset, gs_timer_total, gs_timer_calls, &
        gs_timer_iterations, GS_T_JAC_CG, GS_T_JAC_GMRES, GS_T_PRECOND_CG, &
        GS_T_FACTOR, GS_T_DIRECT
    implicit none
    real(dp), parameter :: r0 = 6.2_dp, f = 32.86_dp
    real(dp), parameter :: c = f/(12*r0**3)
    real(dp) :: a, level
    type(gs_state_t) :: direct_state, iterative_state
    type(gs_profile_law_t) :: law, direct_law, iterative_law
    real(dp), allocatable :: boundary(:, :), input(:), element(:), csr(:)
    real(dp) :: h, difference, per_solve(3), solves
    integer :: mesh, n, j, status

    a = r0/3
    level = 4*r0*r0*a*a
    law%knots = [0.0_dp, 0.31_dp, 0.77_dp, 1.0_dp]
    allocate (law%pp(4, 3), law%ff(4, 3))
    law%normalized = .true.; law%axis = -c*level; law%p_edge = 0; law%f_edge = f
    do j = 1, 3
        law%pp(:, j) = -16*c*[1 + 0.2_dp*law%knots(j)**2, 0.4_dp*law%knots(j), &
                              0.2_dp, 0.0_dp]
    end do
    law%ff = 0
    print '(a)', 'n dof  direct: CG its/calls, GMRES its/calls, factor_s, its/solve'
    do mesh = 1, 3
        n = 32*2**(mesh - 1)
        allocate (boundary(2, n))
        do j = 1, n
            call exact_boundary(real(j - 1, dp)/n, boundary(:, j))
        end do
        h = 2*acos(-1.0_dp)*a/n
        call initialize_gs(direct_state, boundary, h*h/2, status, 3, &
                           curve=exact_boundary)
        if (status /= 0) error stop 'mesh initialization'
        iterative_state = direct_state
        deallocate (boundary)

        ! CSR and element-wise operators, nonzero boundary values included.
        j = size(direct_state%psi)
        allocate (input(j), element(j), csr(j))
        input = [(sin(0.37_dp*j) + 0.1_dp*cos(1.9_dp*j), j=1, size(input))]
        call gs_set_options(.true., .true.)
        call apply_gs(direct_state, input, csr)
        call gs_set_options(.false., .true.)
        call apply_gs(direct_state, input, element)
        difference = maxval(abs(csr - element))/maxval(abs(element))
        if (difference > 1.0e-12_dp) error stop 'CSR operator differs from element-wise'
        deallocate (input, element, csr)

        ! The direct oracles need the MUMPS backend (KIN6D_ENABLE_MUMPS); CTest
        ! reports exit code 77 as skipped (SKIP_RETURN_CODE).
        if (.not. gs_direct_available()) then
            print '(a)', 'SKIP: built without MUMPS; CSR operator check passed'
            stop 77
        end if

        direct_law = law
        call gs_set_options(.true., .true.)
        call gs_timer_reset()
        call solve_gs_cubic_profiles(direct_state, 1.0_dp, direct_law, status, &
            relative_tolerance=1.0e-10_dp, forcing=manufactured_source)
        if (status /= 0) error stop 'direct nonlinear solve'
        if (gs_timer_total(GS_T_FACTOR) <= 0 .or. gs_timer_calls(GS_T_DIRECT) < 1) &
            error stop 'direct path was not used'
        if (gs_timer_iterations(GS_T_PRECOND_CG) /= 0) error stop 'nested CG still used'
        solves = real(gs_timer_calls(GS_T_JAC_CG) + gs_timer_calls(GS_T_JAC_GMRES), dp)
        per_solve(mesh) = real(gs_timer_iterations(GS_T_JAC_CG) &
            + gs_timer_iterations(GS_T_JAC_GMRES), dp)/max(solves, 1.0_dp)
        print '(i5,i8,a,i6,a,i3,a,i6,a,i3,a,f7.3,a,f6.1)', n, size(direct_state%psi), &
            '  ', gs_timer_iterations(GS_T_JAC_CG), '/', gs_timer_calls(GS_T_JAC_CG), &
            '  ', gs_timer_iterations(GS_T_JAC_GMRES), '/', &
            gs_timer_calls(GS_T_JAC_GMRES), '  ', gs_timer_total(GS_T_FACTOR), &
            '  Krylov its per Newton solve:', per_solve(mesh)
        if (per_solve(mesh) > 50) error stop 'Krylov iterations per Newton solve'

        iterative_law = law
        call gs_set_options(.false., .false.)
        call gs_timer_reset()
        call solve_gs_cubic_profiles(iterative_state, 1.0_dp, iterative_law, status, &
            relative_tolerance=1.0e-10_dp, forcing=manufactured_source)
        if (status /= 0) error stop 'iterative nonlinear solve'
        print '(a,i8,a,i8)', '      iterative: CG its', &
            gs_timer_iterations(GS_T_JAC_CG), ' nested CG its', &
            gs_timer_iterations(GS_T_PRECOND_CG)
        difference = maxval(abs(direct_state%psi - iterative_state%psi)) &
            /maxval(abs(iterative_state%psi))
        print '(a,es10.2)', '      psi relative difference:', difference
        if (difference > 1.0e-8_dp) error stop 'direct and iterative flux differ'
        if (abs(direct_law%axis - iterative_law%axis) > 1.0e-8_dp*abs(law%axis)) &
            error stop 'direct and iterative axis flux differ'
    end do
    ! Mesh independence of the preconditioned iteration.
    if (per_solve(3) > 1.5_dp*per_solve(1) + 3) error stop 'iterations grow with mesh'
    print '(a)', 'PASS'
contains
    function manufactured_source(r, z) result(source)
        real(dp), intent(in) :: r, z
        real(dp) :: source, exact_s
        exact_s = ((r*r - r0*r0)**2 + 4*r*r*z*z)/level
        source = 16*c*r*(1 + 0.2_dp*exact_s*exact_s - 1)
    end function

    subroutine exact_boundary(t, point)
        real(dp), intent(in) :: t
        real(dp), intent(out) :: point(2)
        real(dp) :: angle
        angle = 2*acos(-1.0_dp)*t
        point(1) = sqrt(r0*r0 + 2*r0*a*cos(angle))
        point(2) = r0*a*sin(angle)/point(1)
    end subroutine
end program
