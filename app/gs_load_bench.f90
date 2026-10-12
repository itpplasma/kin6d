program gs_load_bench
    !! Timing and equivalence driver for gs_cubic_load on an elliptic cross
    !! section with a prescribed flux. Modes: smooth (no knot inside the flux
    !! range) and knot (one C1 knot inside: smooth rule on the crossed cells).
    !! Usage: gs_load_bench MODE MAX_AREA NREP [write|compare FILE [TOL]]
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_grad_shafranov, only: gs_state_t, initialize_gs
    use kin6d_gs_profile_law, only: gs_profile_law_t
    use kin6d_gs_nonlinear, only: gs_profile_jacobian_t, gs_cubic_load
    use kin6d_gs_timers, only: gs_wall_s, gs_timer_total, gs_timer_reset, GS_T_LOAD_CELLS
    implicit none
    type(gs_state_t) :: state
    type(gs_profile_law_t) :: law
    type(gs_profile_jacobian_t) :: jac
    real(dp), allocatable :: rhs(:), reference(:)
    real(dp) :: contour(2, 64), area, t0, best, rho2, u0, c(4), d(4), diff, norm
    real(dp) :: integral, tol, cells_best
    real(dp), allocatable :: mass_flat(:)
    integer :: j, status, i, n, nrep, rep, unit, evaluations, nread
    character(64) :: mode, text, action, path
    call get_command_argument(1, mode)
    call get_command_argument(2, text); read (text, *) area
    call get_command_argument(3, text); read (text, *) nrep
    action = ''; path = ''; tol = 1.0e-9_dp
    if (command_argument_count() >= 6) then
        call get_command_argument(6, text); read (text, *) tol
    end if
    if (command_argument_count() >= 5) then
        call get_command_argument(4, action)
        call get_command_argument(5, path)
    end if
    do i = 1, 64
        contour(:, i) = [3 + cos(2*acos(-1.0_dp)*real(i - 1, dp)/64), &
                         1.4_dp*sin(2*acos(-1.0_dp)*real(i - 1, dp)/64)]
    end do
    call initialize_gs(state, contour, area, status, 3)
    if (status /= 0) error stop 'mesh'
    n = size(state%psi)
    do i = 1, n
        rho2 = (state%points(1, i) - 3)**2 + (state%points(2, i)/1.4_dp)**2
        state%psi(i) = -0.5_dp*(1 - rho2)*(1 + 0.1_dp*(state%points(1, i) - 3))
    end do
    c = [-0.3_dp, -0.08_dp, -0.02_dp, -0.01_dp]
    d = [-0.2_dp, 0.05_dp, 0.03_dp, 0.01_dp]
    law%f_edge = 5; law%p_edge = 0; law%normalized = .false.
    select case (trim(mode))
    case ('smooth')
        law%knots = [-1.0_dp, 1.0_dp]
        allocate (law%pp(4, 1), law%ff(4, 1)); law%pp(:, 1) = c; law%ff(:, 1) = d
    case default
        law%knots = [-1.0_dp, -0.25_dp, 1.0_dp]
        allocate (law%pp(4, 2), law%ff(4, 2))
        law%pp(:, 1) = c; law%ff(:, 1) = d
        u0 = 0.75_dp
        law%pp(:, 2) = shifted(c, u0); law%ff(:, 2) = shifted(d, u0)
    end select
    allocate (rhs(n))
    best = huge(1.0_dp); cells_best = huge(1.0_dp)
    do rep = 1, nrep
        call gs_timer_reset()
        t0 = gs_wall_s()
        call gs_cubic_load(state, 1.0_dp, law, rhs, jac, status, tol, &
                           quadrature_evaluations=evaluations, &
                           source_integral=integral)
        best = min(best, gs_wall_s() - t0)
        cells_best = min(cells_best, gs_timer_total(GS_T_LOAD_CELLS))
        if (status /= 0) error stop 'load failed'
    end do
    print '(a,a,a,i0,a,i0,a,es12.5,a,es12.5,a,i0,a,es22.14)', 'mode ', trim(mode), &
        ' dof ', n, ' cells ', size(state%triangles, 2), ' best_s ', best, &
        ' cells_s ', cells_best, ' evals ', evaluations, ' source_integral ', integral
    if (trim(action) == 'write') then
        open (newunit=unit, file=trim(path), access='stream', status='replace')
        write (unit) n, rhs, jac%mass, jac%axis_load
        close (unit)
    else if (trim(action) == 'compare') then
        open (newunit=unit, file=trim(path), access='stream', status='old')
        read (unit) nread
        if (nread /= n) error stop 'reference size'
        allocate (reference(n)); read (unit) reference
        diff = maxval(abs(reference - rhs)); norm = maxval(abs(reference))
        print '(a,es10.3)', 'rhs max rel diff ', diff/norm
        deallocate (reference); allocate (reference(size(jac%mass)))
        read (unit) reference
        diff = maxval(abs(reference - reshape(jac%mass, [size(jac%mass)])))
        mass_flat = reshape(jac%mass, [size(jac%mass)])
        j = 0
        do i = 1, size(state%triangles, 2)
            if (maxval(abs(reference((i - 1)*100 + 1:i*100) &
                           - mass_flat((i - 1)*100 + 1:i*100))) > 1.0e-13_dp*norm) &
                j = j + 1
        end do
        print '(a,i0)', 'cells with mass diff > 1e-13 of max: ', j
        norm = maxval(abs(reference))
        print '(a,es10.3)', 'mass max rel diff ', diff/norm
        print '(a,es10.3)', 'mass 2-norm rel diff ', &
            norm2(reference - mass_flat)/norm2(reference)
        deallocate (reference); allocate (reference(n)); read (unit) reference
        diff = maxval(abs(reference - jac%axis_load)); norm = maxval(abs(reference))
        print '(a,es10.3)', 'axis_load max rel diff ', diff/norm
        close (unit)
    end if
contains
    pure function shifted(p, u) result(q)
        real(dp), intent(in) :: p(4), u
        real(dp) :: q(4)
        q(1) = p(1) + u*(p(2) + u*(p(3) + u*p(4)))
        q(2) = p(2) + u*(2*p(3) + 3*u*p(4))
        q(3) = p(3) + 3*u*p(4)
        q(4) = p(4)
    end function
end program
