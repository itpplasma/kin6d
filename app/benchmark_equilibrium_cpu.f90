program benchmark_equilibrium_cpu
    ! Wall-clock profile of the actual 3D map consumer, with an fp64 residual.
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use kin6d_grad_shafranov, only: gs_state_t, initialize_gs
    use kin6d_equilibrium_chart, only: cylindrical_chart_t
    use kin6d_ideal_map, only: ideal_map_t, ideal_map_init, ideal_map_freeze, &
                               ideal_map_energy, ideal_map_hessian, ideal_map_residual
    use fortfem_element_to_csr, only: csr_matrix_t
    use ideal_map_cases, only: setup_case, boundary_point, pullback_field_t, &
                               pi, a, f_edge
    implicit none
    type(gs_state_t) :: state
    type(ideal_map_t) :: map
    type(cylindrical_chart_t) :: chart
    type(pullback_field_t) :: field
    type(csr_matrix_t) :: hessian
    real(dp), allocatable :: contour(:, :), gradient(:)
    real(dp) :: h, t, energy_s, hessian_s, residual_s, energy, residual, sigma
    integer :: nr, ntor, i, rep, repeats, status
    character(32) :: arg
    nr = 32; ntor = 1; repeats = 3
    if (command_argument_count() >= 1) then
        call get_command_argument(1, arg); read (arg, *) nr
    end if
    if (command_argument_count() >= 2) then
        call get_command_argument(2, arg); read (arg, *) ntor
    end if
    if (command_argument_count() >= 3) then
        call get_command_argument(3, arg); read (arg, *) repeats
    end if
    if (nr < 8 .or. ntor < 0 .or. repeats < 1) error stop 'invalid benchmark size'
    call setup_case(.false.)
    allocate (contour(2, nr))
    do i = 1, nr
        call boundary_point(real(i - 1, dp)/nr, contour(:, i))
    end do
    h = 2*pi*a/nr
    call initialize_gs(state, contour, h*h/2, status, 2, curve=boundary_point)
    if (status /= 0) error stop 'mesh'
    call ideal_map_init(map, state, chart, 3, ntor, status, quadrature_order=7)
    if (status /= 0) error stop 'map init'
    call ideal_map_freeze(map, field, status)
    if (status /= 0) error stop 'map freeze'
    allocate (gradient(map%ndof))
    energy_s = huge(1.0_dp); hessian_s = huge(1.0_dp)
    residual_s = huge(1.0_dp)
    do rep = 0, repeats
        t = wall()
        call ideal_map_energy(map, map%x, energy, status, gradient)
        if (status /= 0) error stop 'energy'
        if (rep > 0) energy_s = min(energy_s, wall() - t)
        t = wall()
        call ideal_map_hessian(map, map%x, hessian)
        if (rep > 0) hessian_s = min(hessian_s, wall() - t)
        t = wall()
        call ideal_map_residual(map, f_edge(), residual, sigma, status)
        if (status /= 0) error stop 'fp64 residual'
        if (rep > 0) residual_s = min(residual_s, wall() - t)
    end do
    print '(a)', '# [V] fixed candidate profile; fp64 residual evaluated; '// &
        'error bound unavailable'
  print '(a)', '# nr,ntor,cells,dof,energy_s,hessian_s,residual_s,energy,residual,sigma'
    print '(i0,3(a,i0),6(a,es23.15))', nr, ',', ntor, ',', map%cells, ',', &
        map%ndof, ',', energy_s, ',', hessian_s, ',', residual_s, ',', energy, &
        ',', residual, ',', sigma
contains
    function wall() result(t)
        real(dp) :: t
        integer(int64) :: count, rate
        call system_clock(count, rate)
        t = real(count, dp)/real(rate, dp)
    end function
end program benchmark_equilibrium_cpu
