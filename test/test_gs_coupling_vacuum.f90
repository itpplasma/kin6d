program test_gs_coupling_vacuum
    !! T1 (DESIGN 5.1): coils outside Gamma_c, no plasma. The coupled solution
    !! of a(psi, v) + <S psi, v> = sum_k I_k <g_k, v> is psi = sum_k I_k g_k
    !! inside Gamma_c, with g_k the coil flux of kin6d_gs_coils evaluated
    !! directly. Errors are measured at all nodes and at interior sample
    !! points, and their convergence is checked under mesh refinement (P2 and
    !! P3) and under panel refinement of the BEM (fixed mesh).
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_grad_shafranov, only: gs_state_t
    use kin6d_gs_coils, only: gs_circuit_set_t, gs_circuit_flux_vectors
    use kin6d_gs_coupling, only: gs_coupling_t, gs_coupling_build, gs_coupled_apply, &
                                 gs_coupling_coil_load, gs_coupled_solve, GS_COUPLING_OK
    use gs_coupling_fixture, only: build_state, build_coils, fe_value
    implicit none

    real(dp), parameter :: currents(2) = [1.2e5_dp, -0.8e5_dp]
    real(dp), parameter :: samples(2, 6) = reshape([3.0_dp, 0.0_dp, 3.4_dp, 0.5_dp, &
              2.5_dp, -0.9_dp, 3.9_dp, -0.6_dp, 2.2_dp, 1.0_dp, 3.8_dp, 1.3_dp], [2, 6])
    integer :: failures

    failures = 0
    call mesh_convergence(2, [24, 48, 96, 192], 1, 2.9_dp, 2.0e-6_dp)
    call mesh_convergence(3, [24, 48, 96], 2, 3.8_dp, 5.0e-7_dp)
    call panel_refinement()
    if (failures /= 0) error stop 'gs_coupling vacuum oracle failures'
    print '(a)', 'gs_coupling vacuum (T1) oracles passed'

contains

    subroutine require(ok, what)
        logical, intent(in) :: ok
        character(len=*), intent(in) :: what
        if (.not. ok) then
            failures = failures + 1
            print '(2a)', 'FAIL: ', what
        end if
    end subroutine require

    subroutine solve_case(n_outer, order, panels, errors, iterations)
        !! errors = [relative max node error, relative rms node error,
        !! relative max sample error, relative max trace-node error,
        !! relative residual].
        integer, intent(in) :: n_outer, order, panels
        real(dp), intent(out) :: errors(5)
        integer, intent(out) :: iterations
        type(gs_state_t) :: state
        type(gs_circuit_set_t) :: set
        type(gs_coupling_t) :: coupling
        real(dp), allocatable :: load(:, :), rhs(:), psi(:), exact(:, :), action(:)
        real(dp) :: value, ref(1, 2), scale
        integer :: status, k, found_count
        logical :: found

        call build_state(state, n_outer, order, status)
        call require(status == 0, 'two-loop mesh')
        call build_coils(set)
        call gs_coupling_build(state, coupling, status, panels_per_edge=panels)
        call require(status == GS_COUPLING_OK, 'coupling build')
        call gs_coupling_coil_load(state, coupling, set, load, status)
        call require(status == GS_COUPLING_OK, 'coil load')
        rhs = matmul(load, currents)
        call gs_coupled_solve(coupling, rhs, psi, status, 1.0e-11_dp, iterations)
        call require(status == GS_COUPLING_OK, 'coupled CG converged')
        allocate (exact(size(psi), 2), action(size(psi)))
    call gs_circuit_flux_vectors(set, state % points(1, :), state % points(2, :), exact)
        scale = maxval(abs(matmul(exact, currents)))
        errors(1) = maxval(abs(psi - matmul(exact, currents))) / scale
        errors(2) = sqrt(sum((psi - matmul(exact, currents))**2) / size(psi)) / scale
        errors(4) = maxval(abs(psi(coupling % trace_node) &
                           - matmul(exact(coupling % trace_node, :), currents))) / scale
        errors(3) = 0.0_dp
        found_count = 0
        do k = 1, size(samples, 2)
            call fe_value(state, psi, samples(1, k), samples(2, k), value, found)
            if (.not. found) cycle
            found_count = found_count + 1
            call gs_circuit_flux_vectors(set, samples(1:1, k), samples(2:2, k), ref)
       errors(3) = max(errors(3), abs(value - dot_product(ref(1, :), currents)) / scale)
        end do
        call require(found_count == size(samples, 2), 'sample points located')
        call gs_coupled_apply(coupling, psi, action)
        errors(5) = norm2(rhs - action) / norm2(rhs)
    end subroutine solve_case

    subroutine mesh_convergence(order, levels, panels, min_rate, max_final)
        integer, intent(in) :: order, levels(:), panels
        real(dp), intent(in) :: min_rate, max_final
        real(dp) :: errors(5, size(levels)), rate(size(levels) - 1), t0, t1
        integer :: i, iterations, n

        n = size(levels)
        do i = 1, n
            call cpu_time(t0)
            call solve_case(levels(i), order, panels, errors(:, i), iterations)
            call cpu_time(t1)
            print '(a,i1,a,i4,a,i2,a,5es10.2,a,i5,a,f6.1)', 'P', order, ' edges ', &
                levels(i), ' panels/edge ', panels, '  max rms pt trace res', &
                errors(:, i), '  cg ', iterations, '  s', t1 - t0
            call require(errors(5, i) < 1.0e-10_dp, 'coupled residual')
        end do
        rate = log(errors(2, :n - 1) / errors(2, 2:)) / log(2.0_dp)
        print '(a,i1,a,10f7.2)', 'P', order, ' observed rms rates ', rate
        call require(all(errors(2, 2:) < errors(2, :n - 1)), 'rms error decreases')
        call require(all(errors(1, 2:) < errors(1, :n - 1)), 'max node error decreases')
        call require(all(errors(3, 2:) < errors(3, :n - 1)), 'sample error decreases')
        call require(rate(n - 1) > min_rate, 'rms rate on the finest pair')
        call require(errors(1, n) < max_final, 'accurate on the finest mesh')
    end subroutine mesh_convergence

    subroutine panel_refinement()
        !! Fixed P3 mesh: the BEM trace is the P2 interpolant of the P3 trace,
        !! so refining the panels must reduce the error until the finite
        !! element floor is reached.
        integer, parameter :: panels(3) = [1, 2, 4]
        real(dp) :: errors(5, 3)
        integer :: i, iterations

        do i = 1, 3
            call solve_case(96, 3, panels(i), errors(:, i), iterations)
            print '(a,i2,a,5es10.2)', 'P3 edges 96 panels/edge ', panels(i), &
                '  max rms pt trace res ', errors(:, i)
        end do

       call require(errors(4, 2) < errors(4, 1), 'trace error falls from 1 to 2 panels')
 call require(errors(4, 3) <= 1.02_dp * errors(4, 2), 'trace error settles by 4 panels')
        call require(errors(2, 3) <= 1.02_dp*errors(2, 1), 'rms error stable under panels')
    end subroutine panel_refinement

end program test_gs_coupling_vacuum
