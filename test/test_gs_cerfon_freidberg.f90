module test_gs_cerfon_freidberg_callbacks
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use cerfon_freidberg_fixture, only: cf_case_t, cf_setup, cf_jet, cf_surface, &
        cf_oracle, cf_q
    use kin6d_grad_shafranov, only: gs_state_t, initialize_gs, solve_gs, gs_cell_data
    use kin6d_gs_profile_law, only: gs_profile_law_t
    use kin6d_gs_equilibrium, only: gs_equilibrium
    use kin6d_equilibrium_state, only: equilibrium_state_t, EQ_VALUES, EQ_PSI, &
        EQ_BR, EQ_BZ, EQ_BPHI
    use kin6d_equilibrium_readback, only: equilibrium_readback_t, read_equilibrium
    use kin6d_equilibrium_harness, only: equilibrium_case_t, equilibrium_start_t, &
        equilibrium_resolution_t, equilibrium_run_report_t, run_equilibrium_case
    use fortfem_triangle_duffy_quadrature, only: triangle_duffy_quadrature

    use kin6d_curved_triangle, only: boundary_curve_t
    use kin6d_equilibrium_harness, only: equilibrium_solver_t
    use kin6d_grad_shafranov, only: initialize_gs_object
    implicit none
    private
    public :: run_test_gs_cerfon_freidberg
    real(dp), parameter :: pi = acos(-1.0_dp)
    type :: run_test_gs_cerfon_freidberg_context_t
        real(dp), pointer :: last_error(:) => null()
        type(gs_profile_law_t), pointer :: law => null()
        type(cf_case_t), pointer :: model => null()
        real(dp), pointer :: reference(:) => null()
    end type

    type, extends(boundary_curve_t) :: run_test_gs_cerfon_freidberg_boundary_t
        type(run_test_gs_cerfon_freidberg_context_t), pointer :: data => null()
    contains
        procedure :: evaluate => run_test_gs_cerfon_freidberg_boundary_binding
    end type

    type, extends(equilibrium_solver_t) :: run_test_gs_cerfon_freidberg_solve_case_t
        type(run_test_gs_cerfon_freidberg_context_t), pointer :: data => null()
    contains
        procedure :: evaluate => run_test_gs_cerfon_freidberg_solve_case_binding
    end type

contains
    subroutine run_test_gs_cerfon_freidberg()
    implicit none
    type(cf_case_t) , target :: model
    type(gs_state_t) :: state
    type(gs_profile_law_t) , target :: law
    type(equilibrium_state_t) :: eq
    type(equilibrium_readback_t) :: readback
    type(equilibrium_run_report_t) :: report
    real(dp) , target :: errors(3, 3), last_error(3), rates(2), reference(6), again(6)
    real(dp) :: rho(11), q_error, q_exact, q_coarse, actual(6), globals(6), axis_error
    real(dp) :: elapsed
    integer(int64) :: tick, tock, clock_rate
    integer :: degree, mesh, n, status, j
    character(32) :: name, argument
        type(run_test_gs_cerfon_freidberg_context_t), target :: ctx
        type(run_test_gs_cerfon_freidberg_boundary_t), target :: boundary_operator
        type(run_test_gs_cerfon_freidberg_solve_case_t), target :: solve_case_operator
    name = 'iter-155'; degree = 2
    if (command_argument_count() >= 1) call get_command_argument(1, name)
    if (command_argument_count() >= 2) then
        call get_command_argument(2, argument)
        read (argument, *) degree
    end if
    if (degree < 2 .or. degree > 3) error stop 'CF degree'
    call system_clock(tick, clock_rate)
    call cf_setup(trim(name), model)
    call cf_oracle(model, 256, again)
    call cf_oracle(model, 512, reference)
    print '(a,6es14.5)', '[V] oracle refinement differences=', &
        abs(reference - again)/max(abs(reference), 1.0_dp)
    if (maxval(abs(reference - again)/max(abs(reference), 1.0_dp)) > 1e-10_dp) &
        error stop 'independent analytic volume quadrature unresolved'
    ! Published secondary evidence (Table I) has two decimal places. The
    ! original table gives 1.07, versus 1.0649076 from resolved quadrature;
    ! its 1% comparison is separate from the 1e-3 numerical budgets below.
    if (name == 'nstx-0') then
        if (abs(reference(2)/1.07_dp - 1) > 0.01_dp) error stop 'Table I beta_p'
        if (abs((model%axis - 1)/model%eps - 0.34_dp) > 0.01_dp) &
            error stop 'Table I vacuum-field axis shift'
    else if (name == 'nstx-1') then
        if (abs(reference(2)) > 1e-14_dp) error stop 'Table I zero pressure'
        if (abs((model%axis - 1)/model%eps - 0.11_dp) > 0.01_dp) &
            error stop 'Table I force-free axis shift'
    end if
    print '(3a,i0,a,f7.3)', '[V] ', trim(name), ' P', degree, ' A=', model%a
    print '(a,7es16.7)', '[V] coefficients (eq. 8): ', model%c
    print '(a,6es16.7)', '[V] analytic volume,beta_p,li,Ip,Lp,Phi: ', reference
    ! Constant physical Solov'ev law, R0=Psi0=mu0=F_edge=1.
    allocate (law%knots(2), law%pp(4, 1), law%ff(4, 1))
    law%knots = [model%psi_axis, 0.0_dp]
    law%pp = 0; law%ff = 0
    law%pp(1, 1) = model%a - 1; law%ff(1, 1) = -model%a
    law%axis = model%psi_axis
    do mesh = 1, 3
        n = 64*2**(mesh - 1)
        call bind_run_test_gs_cerfon_freidberg_context()
        call run_equilibrium_case('cerfon-'//trim(name), 'psi', &
            equilibrium_resolution_t(degree, n), 'reference', solve_case_operator, &
            state, report, eq)
        if (report%status /= 0) then
            print *, 'D8 status', report%status
            error stop 'CF D8 solver'
        end if
        if (report%state_bound_status /= 'unavailable' .or. &
            report%claim /= 'candidate') error stop 'CF numerical evidence promoted'
        errors(:, mesh) = last_error
        call bind_run_test_gs_cerfon_freidberg_context()
        print '(a,i4,a,i7,a,3es14.5,a,es12.3)', '[V] n=', n, ' dof=', &
            count(.not. state%boundary), ' relative L2 psi/Bpol/Bphi=', last_error, &
            ' weak residual=', state%residual
    end do
    do mesh = 2, 3
        rates = log(errors(:2, mesh - 1)/errors(:2, mesh))/log(2.0_dp)
        print '(a,i0,a,2f10.5)', '[V] pair ', mesh - 1, ' orders psi/Bpol=', rates
        ! Budgets are frozen after the first measured ladder (see docs).
        if (minval(rates - [real(degree + 1, dp), real(degree, dp)]) < -0.4_dp) &
            error stop 'CF convergence rate'
    end do
    if (errors(1, 3) > 1e-5_dp .or. errors(2, 3) > 5e-4_dp &
        .or. errors(3, 3) > 1e-6_dp) error stop 'CF finest field budget'
    rho = [(sqrt(real(j - 1, dp)/10), j=1, 11)]
    call read_equilibrium(eq, 11, 64, readback, status, rho_grid=rho, &
        bulk_quantities=.true.)
    if (status /= 0) then
        print *, 'D7 status', status
        error stop 'CF D7 readback'
    end if
    if (readback%bulk_status /= 0) error stop 'CF bulk diagnostics unavailable'
    if (readback%bound_status /= 'unavailable' &
        .or. readback%bulk_bound_status /= 'unavailable') &
        error stop 'CF unavailable bounds lost'
    q_error = 0
    do j = 2, 10
        q_exact = cf_q(model, rho(j)**2, 512)
        q_coarse = cf_q(model, rho(j)**2, 256)
        if (abs(q_coarse/q_exact - 1) > 1e-10_dp) &
            error stop 'independent analytic q quadrature unresolved'
        print '(a,f5.2,a,2es16.7)', '[V] psi_n=', rho(j)**2, ' q D7/exact=', &
            readback%q(j), q_exact
        q_error = max(q_error, abs(readback%q(j)/q_exact - 1))
    end do
    actual = [readback%volume, readback%beta_p, readback%li, readback%current, &
        readback%perimeter, readback%phi_edge]
    ! Zero-pressure force-free beta uses a unit absolute floor.
    globals = abs(actual - reference)/max(abs(reference), 1e-2_dp)
    axis_error = norm2(readback%axis - [model%axis, 0.0_dp])/model%eps
    print '(a,6es14.5)', '[V] relative D7 volume,beta_p,li,Ip,Lp,Phi errors=', globals
    print '(a,2es14.5)', '[V] max relative q and axis/eps errors=', q_error, axis_error
    if (.not. all(ieee_is_finite([errors(:, 3), globals, q_error, axis_error]))) &
        error stop 'CF nonfinite diagnostics'
    if (q_error > 1e-3_dp .or. axis_error > 1e-3_dp &
        .or. maxval(globals) > 1e-3_dp) error stop 'CF finest D7 budget'
    call bulk_contract()
    call system_clock(tock)
    elapsed = real(tock - tick, dp)/clock_rate
    print '(a,f10.3)', '[V] total wall seconds=', elapsed
contains
    subroutine bind_run_test_gs_cerfon_freidberg_context()
        ctx%last_error => last_error
        ctx%law => law
        ctx%model => model
        ctx%reference => reference
        boundary_operator%data => ctx
        solve_case_operator%data => ctx
    end subroutine bind_run_test_gs_cerfon_freidberg_context

    subroutine bulk_contract()
        ! A second pressure/current unit scale and a reversed poloidal field
        ! must preserve dimensionless beta_p/li, and reverse signed Ip and q.
        type(gs_profile_law_t) :: reversed
        type(equilibrium_state_t) :: mirrored
        type(equilibrium_readback_t) :: check
        real(dp), parameter :: mu0 = 3.0_dp
        real(dp) :: values(6)
        integer :: stat
        reversed = law
        reversed%knots = [0.0_dp, -model%psi_axis]
        reversed%axis = -model%psi_axis
        reversed%pp = -law%pp/mu0
        reversed%ff = -law%ff
        state%psi = -state%psi
        call gs_equilibrium(state, mirrored, stat, 1.0_dp, model%a, &
            law=reversed, mu0=mu0)
        state%psi = -state%psi
        if (stat /= 0) error stop 'CF reversed core state'
        call read_equilibrium(mirrored, 11, 64, check, stat, rho_grid=rho, &
            bulk_quantities=.true.)
        if (stat /= 0 .or. check%bulk_status /= 0) error stop 'CF reversed D7'
        values = [check%volume, check%beta_p, check%li, -mu0*check%current, &
            check%perimeter, check%phi_edge]
        if (maxval(abs(values - actual)/max(abs(actual), 1.0_dp)) > 1e-11_dp) &
            error stop 'D7 polarity or mu0 normalization'
        if (maxval(abs(check%q + readback%q)/abs(readback%q)) > 1e-11_dp) &
            error stop 'D7 signed q polarity'
        call read_equilibrium(mirrored, 11, 64, check, stat, rho_grid=rho)
        if (stat /= 0 .or. check%bulk_status /= -1) error stop 'D7 bulk opt-in'
        if (ieee_is_finite(check%beta_p) .or. ieee_is_finite(check%li)) &
            error stop 'D7 unavailable bulk values'
    end subroutine

    end subroutine run_test_gs_cerfon_freidberg

    subroutine run_test_gs_cerfon_freidberg_boundary(ctx, t, point)
        type(run_test_gs_cerfon_freidberg_context_t), intent(in), target :: ctx
        real(dp), intent(in) :: t
        real(dp), intent(out) :: point(2)
        real(dp) :: radius, slope
        call cf_surface(ctx%model, 2*pi*t, 0.0_dp, point, radius, slope)
    end subroutine run_test_gs_cerfon_freidberg_boundary

    subroutine run_test_gs_cerfon_freidberg_boundary_binding(self, t, point)
        class(run_test_gs_cerfon_freidberg_boundary_t), intent(in) :: self
        real(dp), intent(in) :: t
        real(dp), intent(out) :: point(2)
        call run_test_gs_cerfon_freidberg_boundary(self%data, t, point)
    end subroutine

    subroutine run_test_gs_cerfon_freidberg_solve_case(ctx, spec, resolution, start, &
        native, report, equilibrium)
        type(run_test_gs_cerfon_freidberg_context_t), intent(in), target :: ctx
        type(equilibrium_case_t), intent(in) :: spec
        type(equilibrium_resolution_t), intent(in) :: resolution
        type(equilibrium_start_t), intent(in) :: start
        class(*), intent(inout) :: native
        type(equilibrium_run_report_t), intent(inout) :: report
        type(equilibrium_state_t), intent(out) :: equilibrium
        real(dp), allocatable :: contour(:, :)
        real(dp) :: h
        integer :: k, stat
        type(run_test_gs_cerfon_freidberg_boundary_t) :: boundary_operator
        boundary_operator%data => ctx
        if (spec%lane /= 'psi' .or. start%name /= 'reference') return
        select type (native)
        type is (gs_state_t)
            allocate (contour(2, resolution%boundary_nodes))
            do k = 1, size(contour, 2)
                call run_test_gs_cerfon_freidberg_boundary(ctx, real(k - 1, &
                    dp)/size(contour, 2), contour(:, k))
            end do
            h = ctx%reference(5)/size(contour, 2)
            call initialize_gs_object(native, contour, h*h/2, stat, resolution%degree, &
                curve=boundary_operator)
            if (stat /= 0) error stop 'CF initialize_gs'
            ! This forward solve receives ONLY the physical boundary and two
            ! constant source coefficients, never nodal values of exact psi.
            call solve_gs(native, 1.0_dp, ctx%model%a - 1, -ctx%model%a, stat)
            if (stat /= 0) error stop 'CF solve_gs'
            call gs_equilibrium(native, equilibrium, stat, 1.0_dp, -ctx%model%a, &
                law=ctx%law, mu0=1.0_dp)
            if (stat /= 0) error stop 'CF core adapter'
            call run_test_gs_cerfon_freidberg_field_errors(ctx, native, equilibrium, &
                ctx%last_error)
            report%status = 0
            report%oracle_error = maxval(ctx%last_error)
            ! Initial zero flux has exactly unit relative psi L2 error.
            report%seed_error = 1
            report%residual = native%residual
            report%iterations = native%iterations
            report%residual_norm = 'relative Euclidean free-row GS weak residual'
        end select
    end subroutine run_test_gs_cerfon_freidberg_solve_case

    subroutine run_test_gs_cerfon_freidberg_solve_case_binding(self, spec, &
        resolution, start, state, report, equilibrium)
        class(run_test_gs_cerfon_freidberg_solve_case_t), intent(in) :: self
        type(equilibrium_case_t), intent(in) :: spec
        type(equilibrium_resolution_t), intent(in) :: resolution
        type(equilibrium_start_t), intent(in) :: start
        class(*), intent(inout) :: state
        type(equilibrium_run_report_t), intent(inout) :: report
        type(equilibrium_state_t), intent(out) :: equilibrium
        call run_test_gs_cerfon_freidberg_solve_case(self%data, spec, resolution, &
            start, state, report, equilibrium)
    end subroutine

    subroutine run_test_gs_cerfon_freidberg_field_errors(ctx, native, equilibrium, &
        errors)
        type(run_test_gs_cerfon_freidberg_context_t), intent(in), target :: ctx
        type(gs_state_t), intent(in) :: native
        type(equilibrium_state_t), intent(in) :: equilibrium
        real(dp), intent(out) :: errors(3)
        real(dp), allocatable :: x(:), y(:), w(:), points(:, :), weights(:)
        real(dp), allocatable :: values(:, :)
        integer, allocatable :: statuses(:)
        real(dp) :: point(2), det, shape(10), gradient(2, 10), hessian(2, 2, 10)
        real(dp) :: jet(6), numerator(3), denominator(3), exact_b(3), error_b(3)
        integer :: e, q, ndof, stat, failures, k, count_points
        call triangle_duffy_quadrature(8, x, y, w, stat)
        if (stat /= 0) error stop 'CF error quadrature'
        ndof = size(native%triangles, 1)
        count_points = size(w)*size(native%triangles, 2)
        allocate (points(2, count_points), weights(count_points), &
            values(EQ_VALUES, count_points), statuses(count_points))
        failures = 0
        !$omp parallel do default(shared) schedule(static) &
        !$omp& private(q, k, point, det, shape, gradient, hessian, stat) &
        !$omp& reduction(+:failures)
        do e = 1, size(native%triangles, 2)
            do q = 1, size(w)
                call gs_cell_data(native, e, [x(q), y(q)], point, det, &
                    shape(:ndof), gradient(:, :ndof), hessian(:, :, :ndof), stat)
                if (stat /= 0) then
                    failures = failures + 1
                    cycle
                end if
                k = (e - 1)*size(w) + q
                points(:, k) = point
                weights(k) = w(q)*det
            end do
        end do
        !$omp end parallel do
        if (failures /= 0) error stop 'CF physical cell'
        ! Compare the actual D7 evaluator, including its field/profile/sign
        ! conventions, rather than reconstructing B from the solved nodal psi.
        call equilibrium%field%sample(points, values, statuses)
        if (any(statuses /= 0)) error stop 'CF D7 quadrature-point evaluation'
        numerator = 0; denominator = 0
        !$omp parallel do default(shared) schedule(static) &
        !$omp& private(point, jet, exact_b, error_b) reduction(+:numerator, denominator)
        do k = 1, count_points
            point = points(:, k)
            call cf_jet(ctx%model, point, jet)
            exact_b = [-jet(3), jet(2), sqrt(1 - 2*ctx%model%a*jet(1))]/point(1)
            error_b = values([EQ_BR, EQ_BZ, EQ_BPHI], k) - exact_b
            numerator = numerator + weights(k)*[(values(EQ_PSI, k) - jet(1))**2, &
                sum(error_b(:2)**2), error_b(3)**2]
            denominator = denominator + weights(k)*[jet(1)**2, &
                sum(exact_b(:2)**2), exact_b(3)**2]
        end do
        !$omp end parallel do
        errors = sqrt(numerator/denominator)
    end subroutine run_test_gs_cerfon_freidberg_field_errors
end module test_gs_cerfon_freidberg_callbacks

program test_gs_cerfon_freidberg
    use test_gs_cerfon_freidberg_callbacks, only: run_test_gs_cerfon_freidberg
    implicit none
    call run_test_gs_cerfon_freidberg()
end program test_gs_cerfon_freidberg
