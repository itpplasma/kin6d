program test_vrmhd_torus
    !! Circular-torus flux-pumping case (docs/vrmhd-torus.md): axisymmetric
    !! Grad-Shafranov state of the cylinder profiles and its stationary
    !! section state.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use milca_vrmhd_cylinder, only: vrmhd_cylinder_case_t, vrmhd_cylinder_profile_t, &
        vrmhd_cylinder_profiles, vrmhd_cylinder_parameters
    use milca_vrmhd_torus, only: vrmhd_torus_law, vrmhd_torus_equilibrium, &
        vrmhd_torus_axis_t, vrmhd_torus_mesh, vrmhd_torus_sources, vrmhd_torus_state
    use milca_vrmhd_parameters, only: vrmhd_parameters_t
    use milca_vrmhd_section, only: vrmhd_section_t, vrmhd_section_build, &
        vrmhd_section_basis
    use milca_vrmhd_stationary, only: vrmhd_stationary_t, vrmhd_stationary_build
    use milca_vrmhd_drive, only: vrmhd_drive_build, VRMHD_LOOP_VOLTAGE_SOURCE
    use milca_vrmhd_steady, only: vrmhd_steady_solve, vrmhd_steady_options_t, &
        vrmhd_steady_report_t
    use milca_vrmhd_evolve, only: vrmhd_evolve_t, vrmhd_evolve_build, &
        vrmhd_evolve_consistent, vrmhd_evolve_expand
    use fortsparse, only: csc_t, csc_from_triplet, fortsparse_status_t, status_ok
    use fortnum_ode_bdf, only: bdf_options_t, bdf_state_t, bdf_report_t, bdf_init, &
        bdf_integrate_to
    use fortnum_status, only: fortnum_status_t, FORTNUM_OK
    use milca_gs_profile_law, only: gs_profile_law_t, gs_law_primitives
    use milca_grad_shafranov, only: gs_state_t, gs_element_jet, gs_axis_safety_factor
    use fortfem_meridian_curved_mesh, only: meridian_mesh_t
    use rm_mhd_oracles, only: shafranov_formula
    implicit none
    character(64) :: gate, arg
    real(dp) :: edge = 0.1_dp, beta_poloidal = 1.0e-6_dp, dtau = 10
    !> 1 holds the initial poloidal current by a source, 0 lets it evolve.
    integer :: poloidal_source = 1
    real(dp) :: t_end = 1.0e4_dp, dt_out = 1.0e3_dp, rtol = 1.0e-4_dp, atol = 1.0e-12_dp
    integer :: max_steps = 100000
    integer :: knots = 129, degree = 1, k, e, failures = 0
    call get_command_argument(1, gate)
    do k = 2, command_argument_count()
        call get_command_argument(k, arg)
        e = index(arg, '=')
        if (e == 0) error stop 'arguments are key=value'
        select case (arg(:e - 1))
        case ('edge'); read (arg(e + 1:), *) edge
        case ('knots'); read (arg(e + 1:), *) knots
        case ('beta_p'); read (arg(e + 1:), *) beta_poloidal
        case ('degree'); read (arg(e + 1:), *) degree
        case ('dtau'); read (arg(e + 1:), *) dtau
        case ('poloidal_source'); read (arg(e + 1:), *) poloidal_source
        case ('t_end'); read (arg(e + 1:), *) t_end
        case ('dt_out'); read (arg(e + 1:), *) dt_out
        case ('rtol'); read (arg(e + 1:), *) rtol
        case ('atol'); read (arg(e + 1:), *) atol
        case ('max_steps'); read (arg(e + 1:), *) max_steps
        case default; error stop 'unknown argument'
        end select
    end do
    select case (trim(gate))
    case ('shafranov')
        call shafranov()
    case ('stationary')
        call stationary(.false.)
    case ('evolve')
        call stationary(.true.)
    case default
        error stop 'expected shafranov or stationary'
    end select
    if (failures /= 0) error stop 'torus checks failed'
contains
    subroutine check(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description
        if (condition) then
            print '(2a)', 'PASS ', description
        else
            print '(2a)', 'FAIL ', description
            failures = failures + 1
        end if
    end subroutine

    ! Axis shift against the large-aspect-ratio formula where the shifted
    ! circles stay well apart (max |dDelta/dr| <= 0.3), and elsewhere against
    ! the polar finite-difference Grad-Shafranov solve of docs/vrmhd-torus.md.
    ! The unstructured mesh is not up-down symmetric; its axis height is a
    ! discretisation error below 1e-4 at edge = 0.1.
    subroutine shafranov()
        type(vrmhd_cylinder_case_t) :: c
        type(vrmhd_cylinder_profile_t) :: prof
        type(gs_profile_law_t) :: law
        type(gs_state_t) :: state
        type(vrmhd_torus_axis_t) :: axis
        real(dp) :: formula, slope, reference, ratio
        real(dp) :: cases(3, 4)
        integer :: i, status
        ! R0/a, beta_p, Delta(0) of the finite-difference solve
        cases = reshape([3.0_dp, 0.02_dp, 0.0755_dp, 3.0_dp, 1.0_dp, 0.157_dp, &
            3.0_dp, 3.0_dp, 0.289_dp, 10.0_dp, 1.0_dp, 0.0476_dp], [3, 4])
        do i = 1, size(cases, 2)
            c%r0 = cases(1, i); c%beta_poloidal = cases(2, i)
            call vrmhd_cylinder_profiles(c, prof, status)
            call check(status == 0, 'cylinder profiles')
            call vrmhd_torus_law(c, prof, knots, law, status)
            call check(status == 0, 'torus profile law')
            call vrmhd_torus_equilibrium(c, prof, edge, law, state, axis, status)
            call check(status == 0, 'torus Grad-Shafranov solve')
            if (status /= 0) cycle
            call shafranov_formula(prof%x, prof%t, prof%x*prof%btheta_over_r, &
                c%r0, formula, slope)
            reference = cases(3, i)
            if (slope <= 0.3_dp) reference = formula
            ratio = axis%shift/reference
            print '(a,f5.1,a,f5.2,a,4f9.5,a,f8.5,a,es10.2)', 'R0/a', c%r0, ' beta_p', &
                c%beta_poloidal, ' shift, formula, FD, ratio', axis%shift, &
                formula, cases(3, i), ratio, ' q0', axis%q, ' Z', axis%position(2)
            call check(abs(axis%position(2)) < 1.0e-3_dp, 'axis on the midplane')
            if (slope <= 0.3_dp) then
                call check(abs(ratio - 1) < 0.03_dp, 'axis shift within 3 % of the formula')
            else
                call check(abs(ratio - 1) < 0.01_dp, &
                    'axis shift within 1 % of the finite-difference solve')
            end if
        end do
    end subroutine

    ! Steady solve from the interpolated Grad-Shafranov state at R0/a = 3 with
    ! isotropic conduction. At p' -> 0 the state is an exact static Ohmic
    ! equilibrium (docs/vrmhd-torus.md), so the solve may change it by
    ! discretisation error only. The section's curl of the interpolated A must
    ! reproduce the Grad-Shafranov field B = (-psi_Z/R, F/R, psi_R/R).
    subroutine stationary(evolve)
        logical, intent(in) :: evolve
        type(vrmhd_cylinder_case_t) :: c
        type(vrmhd_cylinder_profile_t) :: prof
        type(gs_profile_law_t) :: law
        type(gs_state_t) :: gs
        type(vrmhd_torus_axis_t) :: axis
        type(meridian_mesh_t) :: mm
        type(vrmhd_section_t) :: op
        type(vrmhd_stationary_t) :: steady
        type(vrmhd_parameters_t) :: p
        type(vrmhd_steady_options_t) :: options
        type(vrmhd_steady_report_t) :: report
        real(dp), allocatable :: source(:, :, :), js(:, :, :), load(:), harmonic(:)
        real(dp), allocatable :: u0(:), u(:)
        real(dp) :: voltage, field_error, change_te, speed
        integer :: status
        c%r0 = 3; c%beta_poloidal = beta_poloidal; c%chi_ratio = 1
        call vrmhd_cylinder_profiles(c, prof, status)
        call check(status == 0, 'cylinder profiles')
        call vrmhd_torus_law(c, prof, knots, law, status)
        call check(status == 0, 'torus profile law')
        call vrmhd_torus_equilibrium(c, prof, edge, law, gs, axis, status)
        call check(status == 0, 'torus Grad-Shafranov solve')
        call vrmhd_torus_mesh(gs, mm, status)
        call check(status == 0, 'section mesh')
        call vrmhd_section_build(mm, gs%symmetry, degree + 1, 1, op, status, &
            scalar_degree=degree, velocity_degree=degree + 1)
        call check(status == 0, 'section')
        call vrmhd_stationary_build(op, steady, status)
        call check(status == 0, 'gauge')
        if (failures /= 0) return
        call vrmhd_cylinder_parameters(c, prof, p)
        allocate (source(6, op%fs%nq, op%fs%mm%topology%n_triangles), &
            load(op%ndof), harmonic(op%fs%sp(1)%ndof), u0(op%ndof), u(steady%ndof))
        call vrmhd_torus_sources(op, c, prof, gs, law, source, voltage, status)
        call check(status == 0, 'sources')
        js = source(:3, :, :)
        if (poloidal_source == 0) js = 0
        call vrmhd_drive_build(op, VRMHD_LOOP_VOLTAGE_SOURCE, 1.0_dp, voltage, js, &
            harmonic, source, load, status)
        call check(status == 0, 'drive')
        call vrmhd_torus_state(op, gs, law, u0, status)
        call check(status == 0, 'interpolated state')
        if (failures /= 0) return
        if (evolve) then
            print '(a,f9.5)', 'Grad-Shafranov q0 ', axis%q
            call evolve_run(op, p, source, load, u0)
            return
        end if
        field_error = interpolated_field_error(op, gs, law, u0)
        print '(a,i0,a,i0,a,es10.2)', 'triangles ', mm%topology%n_triangles, &
            ' unknowns ', steady%ndof, ' max |B - B_GS|/B0 ', field_error
        call check(field_error < 1.0e-3_dp, 'section field of the interpolated state')
        u = 0; u(:op%ndof) = u0
        options%dtau = dtau; options%verbose = .true.; options%rtol = 1.0e-12_dp
        options%growth = 4
        call vrmhd_steady_solve(op, steady, p, u0, source, load, u, options, report)
        call check(report%status == 0, 'steady solve')
        change_te = maxval(abs(u(op%offset(5) + 1:op%offset(6)) &
            - u0(op%offset(5) + 1:op%offset(6))))/(prof%temperature_axis/2)
        speed = maxval(abs(u(op%offset(2) + 1:op%offset(5))))
        print '(a,es10.2,a,es10.2,a,i0)', 'max change of Te / Te(0) ', change_te, &
            ' max velocity ', speed, ' iterations ', report%iterations
        print '(a,f9.5)', 'Grad-Shafranov q0 ', axis%q
        call axis_report(op, u0, 0.0_dp)
        call axis_report(op, u(:op%ndof), huge(1.0_dp))
        call check(change_te < 1.0e-3_dp, 'stationary temperature')
        call check(speed < 1.0e-5_dp, 'stationary state at rest')
    end subroutine

    ! Axisymmetric BDF2 run from the interpolated state, no helical seed.
    subroutine evolve_run(op, p, source, load, u0)
        type(vrmhd_section_t), intent(in), target :: op
        type(vrmhd_parameters_t), intent(in) :: p
        real(dp), intent(in) :: source(:, :, :), load(:), u0(:)
        type(vrmhd_evolve_t) :: problem
        type(bdf_options_t) :: options
        type(bdf_state_t) :: state
        type(bdf_report_t) :: report
        type(fortnum_status_t) :: st
        type(csc_t) :: reduction
        type(fortsparse_status_t) :: sparse_status
        real(dp), allocatable :: lift(:), y(:), derivative(:), u(:)
        integer, allocatable :: kept(:)
        real(dp) :: t_next
        integer :: i, status, nfree
        kept = pack([(i, i=1, op%ndof)], .not. op%fixed)
        nfree = size(kept)
        call csc_from_triplet(op%ndof, nfree, kept, [(i, i=1, nfree)], &
            [(1.0_dp, i=1, nfree)], reduction, sparse_status)
        call check(status_ok(sparse_status), 'reduction')
        lift = u0
        y = lift(kept)
        lift(kept) = 0
        call vrmhd_evolve_build(problem, op, p, source, lift, reduction, status, &
            load=load, lag=8)
        call check(status == 0, 'evolution build')
        allocate (derivative(nfree), u(op%ndof))
        call vrmhd_evolve_consistent(problem, 0.0_dp, y, derivative, status)
        call check(status == 0, 'consistent start')
        options%adaptive = .true.; options%initial_step = 1
        options%rtol = rtol; options%atol = atol; options%max_steps = max_steps
        options%max_step = dt_out
        call bdf_init(state, 0.0_dp, y, options, st, derivative=derivative)
        call check(st%code == FORTNUM_OK, 'BDF start')
        if (failures /= 0) return
        t_next = 0
        do
            call vrmhd_evolve_expand(problem, state%y, u)
            call axis_report(op, u, state%t)
            if (state%t >= t_end) exit
            t_next = min(t_next + dt_out, t_end)
            call bdf_integrate_to(problem, state, t_next, options, report, st)
            print '(a,i0,a,i0,a,i0,a,i0,a,i0,a,es10.2,a,es10.2)', 'steps ', &
                report%accepted, ' rejected error ', report%rejected_error, &
                ' newton ', report%rejected_newton, ' newton its ', &
                report%newton_iterations, ' krylov ', report%krylov_iterations, &
                ' t ', state%t, ' next step ', state%next_step
            call check(st%code == FORTNUM_OK, 'BDF step')
            if (st%code /= FORTNUM_OK) return
        end do
    end subroutine

    ! Axis and q0 from a quadratic least-squares fit of psi = R A_phi over
    ! the quadrature points within 0.15 a of the flux extremum.
    subroutine axis_report(op, u, time)
        type(vrmhd_section_t), intent(in) :: op
        real(dp), intent(in) :: u(:), time
        interface
            subroutine dgels(trans, m, n, nrhs, a, lda, b, ldb, work, lwork, info)
                import :: dp
                character :: trans
                integer :: m, n, nrhs, lda, ldb, lwork, info
                real(dp) :: a(lda, *), b(ldb, *), work(*)
            end subroutine
        end interface
        real(dp) :: basis(36, op%nlocal), x(2), weight, v(36), best, centre(2)
        real(dp) :: f, hessian(2, 2), a(4000, 6), b(4000, 1), work(4000), d(2)
        real(dp) :: axis(2), q0, speed
        integer :: t, q, status, n, info
        best = 0; centre = 0; f = 0
        do t = 1, op%fs%mm%topology%n_triangles
            do q = 1, op%fs%nq
                call vrmhd_section_basis(op, t, q, x, weight, basis, status)
                v = matmul(basis, u(op%map(:, t)))
                if (abs(v(6)*x(1)) <= best) cycle
                best = abs(v(6)*x(1)); centre = x; f = v(26)*x(1)
            end do
        end do
        n = 0
        do t = 1, op%fs%mm%topology%n_triangles
            do q = 1, op%fs%nq
                call vrmhd_section_basis(op, t, q, x, weight, basis, status)
                d = x - centre
                if (norm2(d) > 0.15_dp .or. n == size(b, 1)) cycle
                v = matmul(basis, u(op%map(:, t)))
                n = n + 1
                a(n, :) = [1.0_dp, d(1), d(2), d(1)**2, d(1)*d(2), d(2)**2]
                b(n, 1) = v(6)*x(1)
            end do
        end do
        call dgels('N', n, 6, 1, a, size(a, 1), b, size(b, 1), work, size(work), info)
        hessian = reshape([2*b(4, 1), b(5, 1), b(5, 1), 2*b(6, 1)], [2, 2])
        d = -matmul(reshape([hessian(2, 2), -hessian(2, 1), -hessian(1, 2), &
            hessian(1, 1)], [2, 2]), b(2:3, 1))/(hessian(1, 1)*hessian(2, 2) &
            - hessian(1, 2)**2)
        axis = centre + d
        call gs_axis_safety_factor(axis(1), f, hessian, q0, status)
        speed = maxval(abs(u(op%offset(2) + 1:op%offset(5))))
        print '(a,f10.1,a,2f10.6,a,f9.5,a,es10.2,a,i0)', 't ', time, ' axis R, Z ', &
            axis, ' q0 ', q0, ' max |v| ', speed, ' fit ', info
    end subroutine

    real(dp) function interpolated_field_error(op, gs, law, u) result(error)
        type(vrmhd_section_t), intent(in) :: op
        type(gs_state_t), intent(in) :: gs
        type(gs_profile_law_t), intent(in) :: law
        real(dp), intent(in) :: u(:)
        real(dp) :: basis(36, op%nlocal), x(2), weight, v(36), psi, grad(2)
        real(dp) :: hessian(2, 2), pressure, f, exact(3)
        integer :: t, q, status
        error = huge(1.0_dp)
        do t = 1, op%fs%mm%topology%n_triangles
            do q = 1, op%fs%nq
                call vrmhd_section_basis(op, t, q, x, weight, basis, status)
                if (status /= 0) return
                call gs_element_jet(gs, t, x(1), x(2), psi, grad, hessian, status)
                if (status /= 0) return
                call gs_law_primitives(law, psi, pressure, f, status)
                if (status /= 0) return
                v = matmul(basis, u(op%map(:, t)))
                exact = [-grad(2), f, grad(1)]/x(1)
                if (t == 1 .and. q == 1) error = 0
                error = max(error, maxval(abs(v(25:27) - exact)))
            end do
        end do
    end function
end program test_vrmhd_torus
