program test_gs_verified_stability
    use, intrinsic :: iso_fortran_env, only: dp => real64, qp => real128
    use omp_lib, only: omp_get_wtime
    use kin6d_grad_shafranov, only: gs_state_t, initialize_gs, gs_cell_data
    use kin6d_gs_profile_law, only: gs_profile_law_t
    use kin6d_gs_nonlinear, only: solve_gs_cubic_profiles
    use kin6d_gs_equilibrated, only: gs_equilibrated_majorant, &
        equilibrated_flux_estimate_t
    use kin6d_gs_error_control, only: gs_stability_t, gs_fixed_stability, &
        gs_fixed_stability_verified, gs_state_error_bound, GS_ERROR_OK, &
        GS_ERROR_UNAVAILABLE, GS_ERROR_NONCOERCIVE
    use fortfem_triangle_duffy_quadrature, only: triangle_duffy_quadrature
    implicit none
    type(gs_state_t) :: state
    type(gs_profile_law_t) :: law
    type(gs_stability_t) :: stability
    type(equilibrated_flux_estimate_t) :: estimate
    real(dp), parameter :: contour(2, 4) = reshape([2.0_dp, -1.0_dp, &
        4.0_dp, -1.0_dp, 4.0_dp, 1.0_dp, 2.0_dp, 1.0_dp], [2, 4])
    real(dp), allocatable :: x(:), y(:), w(:)
    real(dp) :: errors(3), radii(3), start, solve_s, bound_s, radius, rate
    real(dp) :: point(2), det, shape(6), grad(2, 6), hess(2, 2, 6), dg(2)
    real(qp) :: exact_alpha
    integer :: status, mesh, e, q, repeat
    character(16) :: mode

    law%knots = [0.0_dp, 2.0_dp]
    allocate (law%pp(4, 1), law%ff(4, 1))
    law%pp = 0; law%ff = 0
    call initialize_gs(state, contour, 0.08_dp, status, 2)
    if (status /= 0) error stop 'mesh'
    call gs_fixed_stability_verified(state, law, 1.0_dp, stability, status)
    if (status /= GS_ERROR_OK) error stop 'fixed-source verified stability'
    if (.not. stability%arithmetic_rigorous) error stop 'missing arithmetic status'
    if (stability%alpha > 1 .or. stability%alpha < 1 - 1e-14_dp) &
        error stop 'exact-recovery coercivity'
    call gs_state_error_bound(0.0_dp, stability, radius, status)
    if (status /= 0 .or. radius > 1e-100_dp) error stop 'zero residual recovery'
    ! Independent high precision oracle for a derivative maximum away from
    ! the origin (completing the square, not evaluating the generated leaf).
    law%pp(:, 1) = [0.0_dp, -0.1_dp, 0.3_dp, -0.2_dp]
    call gs_fixed_stability_verified(state, law, 1.0_dp, stability, status)
    exact_alpha = 1 - 32*(real(law%pp(2, 1), qp) &
        - real(law%pp(3, 1), qp)**2/(3*real(law%pp(4, 1), qp))) &
        /acos(-1.0_qp)**2
    if (status /= GS_ERROR_OK) error stop 'concave cubic'
    if (real(stability%alpha, qp) > exact_alpha) error stop 'lower bound direction'
    if (exact_alpha - real(stability%alpha, qp) > 1e-13_qp) &
        error stop 'outward bound too loose'
    call gs_state_error_bound(0.2_dp, stability, radius, status)
    if (status /= GS_ERROR_OK) error stop 'verified composition'
    if (real(radius, qp) < real(0.2_dp, qp)/exact_alpha) &
        error stop 'radius rounded inward'
    law%pp = 0; law%pp(2, 1) = 1
    call gs_fixed_stability_verified(state, law, 1.0_dp, stability, status)
    if (status /= GS_ERROR_NONCOERCIVE) error stop 'noncoercive accepted'
    if (stability%available) error stop 'noncoercive availability'
    law%pp = 0; law%pp(4, 1) = 0.01_dp
    call gs_fixed_stability_verified(state, law, 1.0_dp, stability, status)
    if (status /= GS_ERROR_UNAVAILABLE) error stop 'positive cubic tail accepted'
    law%pp = 0; law%pp(3, 1) = -0.01_dp
    call gs_fixed_stability_verified(state, law, 1.0_dp, stability, status)
    if (status /= GS_ERROR_UNAVAILABLE) error stop 'unbounded quadratic accepted'
    law%pp = 0; law%normalized = .true.; law%knots = [0.0_dp, 1.0_dp]
    call gs_fixed_stability_verified(state, law, 1.0_dp, stability, status)
    if (status /= GS_ERROR_UNAVAILABLE) error stop 'normalized law accepted'
    law%normalized = .false.; law%knots = [0.0_dp, 2.0_dp]
    state%curved_cell(1) = .true.
    call gs_fixed_stability_verified(state, law, 1.0_dp, stability, status)
    if (status /= GS_ERROR_UNAVAILABLE) error stop 'curved geometry accepted'
    state%curved_cell(1) = .false.

    ! A genuine nonlinear problem with known polynomial exact flux, zero
    ! trace, and a globally bounded positive one-sided source slope.
    law%pp(:, 1) = [0.0_dp, 0.01_dp, 0.0_dp, -0.02_dp]
    call triangle_duffy_quadrature(14, x, y, w, status)
    if (status /= 0) error stop 'quadrature'
    print '(a)', 'dof energy_error radius efficiency solve_s stability_s ratio'
    do mesh = 1, 3
        call initialize_gs(state, contour, 0.08_dp/4**(mesh - 1), status, 2)
        if (status /= 0) error stop 'manufactured mesh'
        start = omp_get_wtime()
        call solve_gs_cubic_profiles(state, 1.0_dp, law, status, &
            relative_tolerance=1e-11_dp, forcing=forcing)
        solve_s = omp_get_wtime() - start
        if (status /= 0) error stop 'nonlinear solve'
        start = omp_get_wtime()
        do repeat = 1, 20
            call gs_fixed_stability(state, law, 1.0_dp, stability, status)
            if (status /= GS_ERROR_OK) error stop 'automatic verified path'
        end do
        bound_s = (omp_get_wtime() - start)/20
        if (.not. stability%arithmetic_rigorous) error stop 'automatic path status'
        call gs_equilibrated_majorant(state, state%psi, source, estimate, status)
        if (status /= 0) error stop 'equilibrated residual'
        call gs_state_error_bound(estimate%bound, stability, radii(mesh), status)
        if (status /= 0) error stop 'monotone nonlinear bound'
        errors(mesh) = 0
        do e = 1, size(state%triangles, 2)
            do q = 1, size(w)
                call gs_cell_data(state, e, [x(q), y(q)], point, det, &
                    shape, grad, hess, status)
                if (status /= 0) error stop 'cell'
                dg = matmul(grad, state%psi(state%triangles(:, e))) &
                    - [(6 - 2*point(1))*(1 - point(2)**2), &
                       -2*point(2)*(point(1) - 2)*(4 - point(1))]
                errors(mesh) = errors(mesh) + w(q)*det*sum(dg*dg)/point(1)
            end do
        end do
        errors(mesh) = sqrt(errors(mesh))
        if (radii(mesh) < errors(mesh)) error stop 'manufactured containment'
        if (radii(mesh)/errors(mesh) > 1.3_dp) error stop 'efficiency cap'
        print '(i7,6es13.4)', count(.not. state%boundary), errors(mesh), &
            radii(mesh), radii(mesh)/errors(mesh), solve_s, bound_s, bound_s/solve_s
    end do
    do mesh = 2, 3
        rate = log(radii(mesh - 1)/radii(mesh))/log(2.0_dp)
        if (rate < 1.5_dp .or. rate > 2.5_dp) error stop 'P2 majorant rate'
        print '(a,f8.4)', 'majorant rate ', rate
    end do
    ! Small upward jumps can pass the input validity tolerance but cannot
    ! have a finite one-sided Lipschitz constant. Never certify such a law.
    deallocate (law%knots, law%pp, law%ff)
    law%knots = [0.0_dp, 1.0_dp, 2.0_dp]
    allocate (law%pp(4, 2), law%ff(4, 2))
    law%pp = 0; law%ff = 0; law%pp(1, 2) = epsilon(1.0_dp)
    if (.not. law%valid()) error stop 'jump fixture'
    call gs_fixed_stability_verified(state, law, 1.0_dp, stability, status)
    if (status /= GS_ERROR_UNAVAILABLE) error stop 'profile jump certified'
    call get_command_argument(1, mode)
    if (trim(mode) == 'scan') call scan_benchmark()
contains
    subroutine scan_benchmark()
        ! Optional geometry-scan microbenchmark; no solve timings are inferred
        ! from it. Only the points, affine flags and ready marker are consumed.
        type(gs_state_t) :: domain
        type(gs_profile_law_t) :: profile
        type(gs_stability_t) :: checked
        integer :: n, level, trial, i, j, node, stat, rep
        real(dp) :: clock, samples(5)
        profile%knots = [0.0_dp, 2.0_dp]
        allocate (profile%pp(4, 1), profile%ff(4, 1))
        profile%pp(:, 1) = [0.0_dp, 0.01_dp, 0.0_dp, -0.02_dp]
        profile%ff = 0
        domain%ready = .true.
        do level = 1, 3
            n = 256*2**(level - 1)
            allocate (domain%points(2, (n + 1)**2), domain%curved_cell(2*n*n))
            domain%curved_cell = .false.
            do j = 0, n
                do i = 0, n
                    node = 1 + i + (n + 1)*j
                    domain%points(:, node) = &
                        [2 + 2*real(i, dp)/n, -1 + 2*real(j, dp)/n]
                end do
            end do
            do trial = 1, 5
                clock = omp_get_wtime()
                do rep = 1, 100
                    call gs_fixed_stability_verified(domain, profile, &
                        1.0_dp, checked, stat)
                    if (stat /= 0) error stop 'scan stability'
                end do
                samples(trial) = (omp_get_wtime() - clock)/100
            end do
            print '(a,i9,6es13.4)', 'scan nodes, alpha, seconds ', &
                size(domain%points, 2), checked%alpha, samples
            deallocate (domain%points, domain%curved_cell)
        end do
    end subroutine scan_benchmark
    function exact(r, z) result(p)
        real(dp), intent(in) :: r, z
        real(dp) :: p
        p = (r - 2)*(4 - r)*(1 - z*z)
    end function exact
    function forcing(r, z) result(f)
        real(dp), intent(in) :: r, z
        real(dp) :: f, p
        p = exact(r, z)
        f = (2*(1 - z*z) + 2*(r - 2)*(4 - r))/r &
            + (6 - 2*r)*(1 - z*z)/(r*r) - r*(0.01_dp*p - 0.02_dp*p**3)
    end function forcing
    function source(cell, r, z, p) result(f)
        integer, intent(in) :: cell
        real(dp), intent(in) :: r, z, p
        real(dp) :: f
        f = forcing(r, z) + r*(0.01_dp*p - 0.02_dp*p**3) + 0.0_dp*cell
    end function source
end program test_gs_verified_stability
