program test_gs_cuda_public_solve
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use kin6d_grad_shafranov, only: gs_state_t, initialize_gs, solve_gs, &
        assemble_gs_rhs, apply_gs, gs_apply_inverse, gs_requested_device
    use kin6d_gs_direct, only: gs_direct_free
    use kin6d_gs_options, only: gs_set_options
    use kin6d_gs_timers, only: gs_wall_s
    use nn_solovev_fixture, only: boundary, minor, pp
    implicit none
    integer(int64), parameter :: budget = 1024_int64**3
    type(gs_state_t) :: cpu, gpu, shared
    real(dp), allocatable :: contour(:, :)
    real(dp) :: max_area, clock, cpu_init, gpu_init
    real(dp) :: cpu_first, gpu_first, gpu_reuse, gpu_reset, residual, energy
    real(dp) :: difference, pairing, reference_energy
    real(dp), allocatable :: rhs(:), action(:), first(:), delta(:)
    real(dp), allocatable :: goal(:), tangent(:), adjoint(:)
    integer :: j, status, n
    logical :: ok
    character(24) :: argument

    n = 16
    if (command_argument_count() > 1) error stop 'expected boundary resolution only'
    if (command_argument_count() == 1) then
        call get_command_argument(1, argument)
        read (argument, *, iostat=status) n
        if (status /= 0) error stop 'invalid boundary resolution'
    end if
    if (n < 8 .or. n > 64) error stop 'boundary resolution must be 8..64'
    allocate (contour(2, n))
    call gs_set_options(.true., .true.)
    do j = 1, n
        call boundary(real(j - 1, dp)/n, contour(:, j))
    end do
    max_area = (2*acos(-1.0_dp)*minor/n)**2/2
    clock = gs_wall_s()
    call initialize_gs(cpu, contour, max_area, status, 2, curve=boundary)
    cpu_init = gs_wall_s() - clock
    call require(status == 0, 'CPU public initialization')
    clock = gs_wall_s()
    call initialize_gs(gpu, contour, max_area, status, 2, curve=boundary, &
        device=0, device_memory_limit_bytes=budget)
    gpu_init = gs_wall_s() - clock
    call require(status == 0, 'GPU public initialization')
    call require(size(cpu%psi) == size(gpu%psi), 'matched state dimension')
    call require(all(shape(cpu%points) == shape(gpu%points)), 'matched point count')
    call require(maxval(abs(cpu%points - gpu%points)) == 0, 'matched physical mesh')
    call require(all(cpu%triangles == gpu%triangles), 'matched topology')

    clock = gs_wall_s()
    call solve_gs(cpu, 1.0_dp, pp, 0.0_dp, status)
    cpu_first = gs_wall_s() - clock
    call require(status == 0, 'CPU public solve')
    clock = gs_wall_s()
    call solve_gs(gpu, 1.0_dp, pp, 0.0_dp, status)
    gpu_first = gs_wall_s() - clock
    call require(status == 0, 'GPU public solve without fallback')
    call require(gs_requested_device(gpu) == 0, 'explicit state device')
    call require(gpu%direct%stats%active_device == 0, 'actual native GPU factor')
    call require(gpu%direct%stats%device_peak_bytes > 0, 'positive GPU payload')
    call require(gpu%direct%stats%device_peak_bytes <= budget, 'GPU payload cap')
    call require(cpu%direct%stats%active_device == -1, 'CPU remains independent')
    call require(.not. associated(cpu%direct, gpu%direct), 'independent caches')
    allocate (rhs(size(gpu%psi)), action(size(gpu%psi)), delta(size(gpu%psi)))
    call assemble_gs_rhs(gpu, 1.0_dp, pp, 0.0_dp, rhs, status)
    call require(status == 0, 'shared source')
    call gs_set_options(.false., .true.)
    call apply_gs(gpu, gpu%psi, action)
    residual = norm2(action - rhs)/norm2(rhs)
    call require(residual <= 1.0e-9_dp, 'independent element residual')
    call apply_gs(cpu, cpu%psi, action)
    reference_energy = sqrt(max(dot_product(cpu%psi, action), tiny(1.0_dp)))
    delta = gpu%psi - cpu%psi
    call apply_gs(gpu, delta, action)
    energy = sqrt(max(dot_product(delta, action), 0.0_dp))/reference_energy
    difference = maxval(abs(delta))/maxval(abs(cpu%psi))
    call require(max(energy, difference) <= 1.0e-9_dp, 'same public solution')

    allocate (goal(size(rhs)), tangent(size(rhs)), adjoint(size(rhs)))
    goal = [(sin(0.37_dp*j), j=1, size(rhs))]
    where (gpu%boundary) goal = 0
    call gs_apply_inverse(gpu, rhs, tangent, ok)
    call require(ok, 'public equation tangent inverse')
    call gs_apply_inverse(gpu, goal, adjoint, ok)
    call require(ok, 'public equation adjoint inverse')
    pairing = abs(dot_product(goal, tangent) - dot_product(rhs, adjoint)) &
        /max(norm2(goal)*norm2(tangent), tiny(1.0_dp))
    call require(pairing <= 1.0e-10_dp, 'equation-level tangent adjoint pairing')

    first = gpu%psi
    shared = gpu
    call require(associated(shared%direct, gpu%direct), 'copy reuses factor')
    clock = gs_wall_s()
    call solve_gs(shared, 1.0_dp, 0.75_dp*pp, 0.0_dp, status)
    gpu_reuse = gs_wall_s() - clock
    call require(status == 0, 'copied state cached solve')
    call require(maxval(abs(shared%psi - 0.75_dp*first)) &
        <= 1.0e-9_dp*maxval(abs(first)), 'independent linear source scaling')
    call require(shared%direct%stats%active_device == 0, 'reuse stays GPU')
    call gs_direct_free(shared%direct)
    clock = gs_wall_s()
    call solve_gs(shared, 1.0_dp, pp, 0.0_dp, status)
    gpu_reset = gs_wall_s() - clock
    call require(status == 0, 'reset public GPU factor')
    call require(shared%direct%stats%active_device == 0, 'reset restores device')
    call require(shared%direct%stats%device_peak_bytes > 0, 'reset GPU payload')
    call require(shared%direct%stats%device_peak_bytes <= budget, 'reset cap')
    call require(maxval(abs(shared%psi - first)) &
        <= 1.0e-9_dp*maxval(abs(first)), 'reset preserves physical solution')
    call gs_set_options(.true., .true.)
    print '(a,i0)', 'public_GS_boundary_resolution=', n
    print '(a,i0)', 'public_GS_dofs=', size(rhs)
    print '(a,2es17.9)', 'public_initialization_cpu_gpu_s=', cpu_init, gpu_init
    print '(a,2es17.9)', 'public_first_solve_cpu_gpu_s=', cpu_first, gpu_first
    print '(a,2es17.9)', 'public_GPU_reuse_reset_s=', gpu_reuse, gpu_reset
    print '(a,3es17.9)', 'public_element_residual_flux_energy=', &
        residual, difference, energy
    print '(a,es17.9)', 'public_implicit_pairing=', pairing
    print '(a,i0)', 'public_GPU_device_payload_bytes=', &
        shared%direct%stats%device_peak_bytes
    print '(a)', 'claim=candidate; fixed-source SPD GS; host assembly/checking'
    print '(a)', 'PASS: public fixed-boundary GPU solve, reuse and reset'
    call gs_direct_free(cpu%direct)
    call gs_direct_free(shared%direct)
contains
    subroutine require(ok, message)
        logical, intent(in) :: ok
        character(*), intent(in) :: message
        if (ok) return
        print '(a)', message
        error stop 'public GS CUDA solve'
    end subroutine
end program test_gs_cuda_public_solve
