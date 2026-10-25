program test_gs_cuda_fixed_source
    !! Host-assembled fixed-source GS, native CPU/CUDA SPD factors, CPU checks.
    !! This is a host-staged solve prototype, not a resident GPU equilibrium.
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use kin6d_grad_shafranov, only: gs_state_t, initialize_gs, assemble_gs_rhs, &
        apply_gs, gs_cell_data
    use kin6d_gs_direct, only: gs_direct_t, gs_direct_factor, gs_direct_solve, &
        gs_direct_free
    use kin6d_gs_csr, only: gs_csr_spd_matrix
    use kin6d_gs_options, only: gs_set_options
    use kin6d_gs_timers, only: gs_wall_s
    use kin6d_gs_equilibrated, only: gs_equilibrated_majorant, &
        equilibrated_flux_estimate_t
    use kin6d_gs_quadrature_budget, only: gs_checked_equilibrated_majorant
    use kin6d_grad_shafranov_leaves, only: gs_source
    use fortnum_interval, only: interval_t
    use fortfem_triangle_duffy_quadrature, only: triangle_duffy_quadrature
    use nn_solovev_fixture, only: boundary, exact, minor, pp
    implicit none
    integer, parameter :: repeats = 5
    type(gs_state_t) :: state
    type(gs_direct_t) :: cpu, gpu
    type(equilibrated_flux_estimate_t) :: estimate_cpu, estimate_gpu
    type(interval_t) :: audit_cpu, audit_gpu
    real(dp), allocatable :: contour(:, :), rhs(:), xc(:), xg(:), action(:)
    real(dp), allocatable :: values(:), delta(:), tangent(:), goal(:), adjoint(:)
    integer, allocatable :: rowptr(:), colidx(:)
    real(dp) :: max_area, pprime, ffprime, start, assembly_s, checking_s
    real(dp) :: cpu_setup, gpu_setup, cpu_solve(repeats), gpu_solve(repeats)
    real(dp) :: residual_s(2), estimator_s(2), audit_s(2), oracle_s(2), clock
    real(dp) :: rc, rg, difference, energy_difference, reference_energy, pairing
    real(dp) :: cpu_error, gpu_error, estimator_difference, scale
    integer :: n, j, status, device, audit_status_cpu, audit_status_gpu
    logical :: ok
    character(24) :: mode, argument

    mode = 'a3'; n = 16; device = 0
    if (command_argument_count() > 0) call get_command_argument(1, mode)
    if (trim(mode) /= 'a3' .and. trim(mode) /= 'affine') &
        error stop 'supported modes: fixed-source a3 or affine; nonlinear unavailable'
    if (command_argument_count() > 1) then
        call get_command_argument(2, argument)
        read (argument, *, iostat=status) n
        if (status /= 0) error stop 'invalid boundary resolution'
    end if
    if (command_argument_count() > 2) then
        call get_command_argument(3, argument)
        read (argument, *, iostat=status) device
        if (status /= 0) error stop 'invalid CUDA device'
    end if
    if (command_argument_count() > 3) error stop 'unexpected argument'
    if (n < 8 .or. n > 64) error stop 'prototype boundary resolution must be 8..64'
    if (device < 0) error stop 'explicit CUDA device required; no CPU fallback'
    start = gs_wall_s()
    if (trim(mode) == 'a3') then
        allocate (contour(2, n))
        do j = 1, n
            call boundary(real(j - 1, dp)/n, contour(:, j))
        end do
        max_area = (2*acos(-1.0_dp)*minor/n)**2/2
        call initialize_gs(state, contour, max_area, status, 2, curve=boundary)
        pprime = pp; ffprime = 0
    else
        contour = reshape([2.0_dp, -1.0_dp, 4.0_dp, -1.0_dp, &
            4.0_dp, 1.0_dp, 2.0_dp, 1.0_dp], [2, 4])
        max_area = 16.0_dp/real(n*n, dp)
        call initialize_gs(state, contour, max_area, status, 2)
        pprime = 0.1_dp; ffprime = 0.2_dp
    end if
    if (status /= 0) error stop 'GS mesh initialization'
    allocate (rhs(size(state%psi)), xc(size(state%psi)), xg(size(state%psi)), &
        action(size(state%psi)), delta(size(state%psi)))
    call assemble_gs_rhs(state, 1.0_dp, pprime, ffprime, rhs, status)
    if (status /= 0) error stop 'shared GS source assembly'
    call gs_csr_spd_matrix(state%op, state%boundary, rowptr, colidx, values)
    assembly_s = gs_wall_s() - start

    start = gs_wall_s()
    call gs_direct_factor(cpu, size(rhs), rowptr, colidx, values, ok)
    if (.not. ok) error stop 'native CPU GS factor'
    cpu_setup = gs_wall_s() - start
    ! The limit covers explicit backend-owned device payload, not runtime/RSS.
    gpu%solver%options%device_memory_limit_bytes = 1024_int64**3
    start = gs_wall_s()
    call gs_direct_factor(gpu, size(rhs), rowptr, colidx, values, ok, device=device)
    if (.not. ok) then
        print '(a,i0,2a)', 'CUDA_factor_status=', gpu%factor_status%code, &
            ' ', trim(gpu%factor_status%msg)
        error stop 'explicit CUDA GS factor unavailable or failed'
    end if
    gpu_setup = gs_wall_s() - start
    print '(a,i0)', 'gpu_explicit_host_csc_copy_bytes=', &
        int(size(rowptr), int64)*(storage_size(rowptr(1))/8) &
        + int(size(values), int64)*(storage_size(colidx(1))/8 &
        + storage_size(values(1))/8)
    print '(a,i0)', 'gpu_constructor_triplet_row_bytes=', &
        int(size(values), int64)*(storage_size(colidx(1))/8)
    if (gpu%stats%active_device /= device) error stop 'CUDA request fell back'
    if (gpu%stats%device_peak_bytes <= 0) error stop 'missing CUDA payload evidence'
    if (gpu%stats%device_peak_bytes > 1024_int64**3) error stop 'CUDA payload limit'
    do j = 1, repeats
        start = gs_wall_s()
        call gs_direct_solve(cpu, rhs, xc, ok)
        cpu_solve(j) = gs_wall_s() - start
        if (.not. ok) error stop 'native CPU GS solve'
        start = gs_wall_s()
        call gs_direct_solve(gpu, rhs, xg, ok)
        gpu_solve(j) = gs_wall_s() - start
        if (.not. ok) error stop 'native CUDA GS solve'
    end do
    if (.not. all(ieee_is_finite(xg))) error stop 'nonfinite CUDA candidate'

    start = gs_wall_s()
    ! Independent element action: do not reuse the CSR factor's multiplication.
    call gs_set_options(.false., .true.)
    clock = gs_wall_s()
    call apply_gs(state, xc, action)
    rc = norm2(action - rhs)/max(norm2(rhs), tiny(1.0_dp))
    residual_s(1) = gs_wall_s() - clock
    reference_energy = sqrt(max(dot_product(xc, action), 0.0_dp))
    clock = gs_wall_s()
    call apply_gs(state, xg, action)
    rg = norm2(action - rhs)/max(norm2(rhs), tiny(1.0_dp))
    residual_s(2) = gs_wall_s() - clock
    delta = xg - xc
    call apply_gs(state, delta, action)
    energy_difference = sqrt(max(dot_product(delta, action), 0.0_dp)) &
        /max(reference_energy, tiny(1.0_dp))
    scale = max(maxval(abs(xc)), tiny(1.0_dp))
    difference = maxval(abs(delta))/scale
    if (max(rc, rg) > 1.0e-9_dp) error stop 'CPU element residual gate'
    if (max(difference, energy_difference) > 1.0e-9_dp) &
        error stop 'CPU/CUDA GS solution or energy disagreement'

    ! Equation-level tangent/adjoint systems reuse the same GPU factor.
    allocate (tangent(size(rhs)), goal(size(rhs)), adjoint(size(rhs)))
    goal = [(sin(0.37_dp*j), j=1, size(rhs))]
    where (state%boundary) goal = 0
    call gs_direct_solve(gpu, rhs, tangent, ok)
    if (.not. ok) error stop 'implicit tangent solve'
    call gs_direct_solve(gpu, goal, adjoint, ok)
    if (.not. ok) error stop 'implicit adjoint solve'
    pairing = abs(dot_product(goal, tangent) - dot_product(rhs, adjoint)) &
        /max(norm2(goal)*norm2(tangent), tiny(1.0_dp))
    if (pairing > 1.0e-10_dp) error stop 'implicit GPU tangent/adjoint pairing'
    clock = gs_wall_s()
    call gs_equilibrated_majorant(state, xc, source, estimate_cpu, status)
    estimator_s(1) = gs_wall_s() - clock
    if (status /= 0) error stop 'CPU candidate equilibrated estimator'
    clock = gs_wall_s()
    call gs_equilibrated_majorant(state, xg, source, estimate_gpu, status)
    estimator_s(2) = gs_wall_s() - clock
    if (status /= 0) error stop 'GPU candidate CPU equilibrated estimator'
    estimator_difference = abs(estimate_cpu%bound - estimate_gpu%bound) &
        /max(estimate_cpu%bound, tiny(1.0_dp))
    if (estimator_difference > 1.0e-7_dp) error stop 'CPU estimator disagreement'
    clock = gs_wall_s()
    call gs_checked_equilibrated_majorant(state, xc, estimate_cpu, &
        1.0_dp, pprime, ffprime, audit_cpu, audit_status_cpu)
    audit_s(1) = gs_wall_s() - clock
    clock = gs_wall_s()
    call gs_checked_equilibrated_majorant(state, xg, estimate_gpu, &
        1.0_dp, pprime, ffprime, audit_gpu, audit_status_gpu)
    audit_s(2) = gs_wall_s() - clock
    if (trim(mode) == 'affine') then
        if (max(audit_status_cpu, audit_status_gpu) /= 0) &
            error stop 'supported affine strict CPU audit failed'
    else
        if (audit_status_cpu == 0 .or. audit_status_gpu == 0) &
            error stop 'curved A3 audit incorrectly claimed available'
    end if
    cpu_error = -1; gpu_error = -1; oracle_s = 0
    if (trim(mode) == 'a3') then
        clock = gs_wall_s()
        cpu_error = exact_energy_error(xc)
        oracle_s(1) = gs_wall_s() - clock
        clock = gs_wall_s()
        gpu_error = exact_energy_error(xg)
        oracle_s(2) = gs_wall_s() - clock
        if (max(cpu_error, gpu_error) > 0.1_dp*reference_energy) &
            error stop 'public A3 exact-oracle error'
    end if
    checking_s = gs_wall_s() - start
    call gs_set_options(.true., .true.)

    print '(a,a,a,i0,a,i0,a,i0)', 'case=', trim(mode), ' n=', n, &
        ' dofs=', size(rhs), ' degree=', state%element_order
    print '(a)', 'claim=candidate; CPU assembly/checking; host-staged GPU SPD solve'
    print '(a,es16.8)', 'assembly_s=', assembly_s
    print '(a,es16.8)', 'cpu_setup_factor_wall_s=', cpu_setup
    print '(a,es16.8)', 'gpu_setup_factor_wall_s=', gpu_setup
    print '(a,5es16.8)', 'cpu_solve_roundtrip_s=', cpu_solve
    print '(a,5es16.8)', 'gpu_solve_roundtrip_s=', gpu_solve
    print '(a,2es16.8)', 'cpu_residual_evaluation_s_cpu_gpu=', residual_s
    print '(a,2es16.8)', 'cpu_estimator_s_cpu_gpu=', estimator_s
    print '(a,2es16.8)', 'strict_cpu_audit_s_cpu_gpu=', audit_s
    print '(a,2es16.8)', 'exact_oracle_validation_s_cpu_gpu=', oracle_s
    print '(a,es16.8)', 'cpu_first_checked_result_s=', assembly_s + cpu_setup &
        + cpu_solve(1) + residual_s(1) + estimator_s(1) + audit_s(1)
    print '(a,es16.8)', 'gpu_first_checked_result_s=', assembly_s + gpu_setup &
        + gpu_solve(1) + residual_s(2) + estimator_s(2) + audit_s(2)
    print '(a)', 'transfer_timing=unavailable-separately; included-in-roundtrip'
    print '(a,i0)', 'solve_transfer_bytes_per_rhs=', 16_int64*size(rhs)
    print '(a,i0)', 'gpu_active_device=', gpu%stats%active_device
    print '(a,i0)', 'gpu_device_peak_payload_bytes=', gpu%stats%device_peak_bytes
    print '(a,i0)', 'gpu_host_peak_payload_bytes=', gpu%stats%host_peak_bytes
    print '(a,i0)', 'cpu_host_peak_payload_bytes=', cpu%stats%host_peak_bytes
    print '(a,es16.8)', 'cpu_element_residual_relative=', rc
    print '(a,es16.8)', 'gpu_element_residual_relative=', rg
    print '(a,es16.8)', 'flux_cpu_gpu_relative=', difference
    print '(a,es16.8)', 'energy_cpu_gpu_relative=', energy_difference
    print '(a,es16.8)', 'implicit_pairing_relative=', pairing
    print '(a,2es16.8)', 'sampled_energy_estimator_cpu_gpu=', &
        estimate_cpu%bound, estimate_gpu%bound
    print '(a,2i4)', 'strict_cpu_audit_status_cpu_gpu=', &
        audit_status_cpu, audit_status_gpu
    if (audit_status_gpu == 0) then
        print '(a,2es16.8)', 'certified_continuous_energy_radius_cpu_gpu=', &
            audit_cpu%hi, audit_gpu%hi
    else
        print '(a)', 'continuous_energy_certificate=unavailable-curved-geometry'
    end if
    print '(a,2es16.8)', 'A3_exact_energy_error_cpu_gpu_minus1_unavailable=', &
        cpu_error, gpu_error
    print '(a,es16.8)', 'cpu_checking_s=', checking_s
    call gs_direct_free(cpu)
    call gs_direct_free(gpu)
    print '(a)', 'PASS'
contains
    real(dp) function source(cell, r, z, psi) result(value)
        integer, intent(in) :: cell
        real(dp), intent(in) :: r, z, psi
        call gs_source(r, 1.0_dp, pprime, ffprime, value)
        value = value + 0.0_dp*(cell + z + psi)
    end function source

    real(dp) function exact_energy_error(solution) result(error)
        real(dp), intent(in) :: solution(:)
        real(dp), allocatable :: u(:), v(:), w(:)
        real(dp) :: point(2), det, shape(6), gradient(2, 6), hessian(2, 2, 6)
        real(dp) :: psi_exact, b(2), defect(2)
        integer :: e, q, checked
        call triangle_duffy_quadrature(10, u, v, w, checked)
        if (checked /= 0) error stop 'exact-oracle quadrature'
        error = 0
        do e = 1, size(state%triangles, 2)
            do q = 1, size(w)
                call gs_cell_data(state, e, [u(q), v(q)], point, det, shape, &
                    gradient, hessian, checked)
                if (checked /= 0) error stop 'exact-oracle geometry'
                call exact(point, psi_exact, b)
                defect = matmul(gradient, solution(state%triangles(:, e))) &
                    - point(1)*[b(2), -b(1)]
                error = error + w(q)*det*sum(defect**2)/point(1)
            end do
        end do
        error = sqrt(max(error, 0.0_dp))
    end function exact_energy_error
end program test_gs_cuda_fixed_source
