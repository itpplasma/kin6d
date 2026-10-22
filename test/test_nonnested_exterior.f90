module nn_exterior_case
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp), r0 = 3.0_dp, minor = 0.5_dp
contains
    subroutine mapping(u, theta, x)
        real(dp), intent(in) :: u, theta
        real(dp), intent(out) :: x(2)
        x = [r0 + minor*u*cos(theta), minor*u*sin(theta)]
    end subroutine
    function source(x) result(s)
        real(dp), intent(in) :: x(2)
        real(dp) :: s
        s = 0*x(1)
    end function
end module

program test_nonnested_exterior
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use nn_exterior_case, only: pi, r0, minor, mapping, source
    use kin6d_mhd_mesh, only: flux_mesh_t, flux_mesh_build
    use kin6d_coils_3d, only: coil_set_t, coil_set_from_polylines, &
                              coil_set_b_cart, coil_mu0
    use kin6d_equilibrium_nonnested, only: free_b_field_t, nn_build, nn_point
    use kin6d_nonnested_coupled, only: nc_problem_t, nc_layout_t, nc_report_t, &
                                 nc_layout, nc_solve, nc_fixed, nc_residual, nc_jacobian
    use kin6d_nonnested_exterior, only: nn_exterior_t, nn_attach_exterior, &
                                        nn_exterior_trace
    use fortfem_parametric_bem_high_order, only: bem_pair_options_t, bem_pair_report_t
    use kin6d_equilibrium_state, only: EQ_VALUES, EQ_BR, EQ_BPHI, EQ_BZ
    use fortsparse, only: csc_t, csc_matvec
    implicit none
    type(flux_mesh_t) :: mesh
    type(coil_set_t) :: coils
    type(free_b_field_t) :: field
    type(nc_problem_t) :: problem
    type(nc_layout_t) :: lay
    type(nc_report_t) :: newton
    type(nn_exterior_t) :: exterior
    type(bem_pair_options_t) :: options
    type(bem_pair_report_t) :: report
    type(csc_t) :: jacobian
    integer, parameter :: npc = 64
    real(dp) :: pts(3, 64*(npc + 1)), currents(npc + 1), phi, theta, rz(2), xq(3), bc(3)
    real(dp) :: values(EQ_VALUES), direct(3), b(3), error, scale, defect, seed_error
    real(dp) :: step, jac_error, start, finish, physical_defect, solve_seconds
    integer(int64) :: solve_tick, end_tick, solve_rate
   real(dp), allocatable :: rings(:), state(:), saved(:), potential(:), flux(:), load(:)
    real(dp), allocatable :: direction(:), rp(:), rm(:)
    logical, allocatable :: fixed(:)
    integer :: nr, modes, nphi, degree, k, j, i, status, counts(npc + 1)
    character(32) :: arg
    nr = 1; modes = 1; nphi = 8; degree = 12
    if (command_argument_count() >= 1) then
        call get_command_argument(1, arg); read (arg, *) nr
    end if
    if (command_argument_count() >= 2) then
        call get_command_argument(2, arg); read (arg, *) nphi
    end if
    if (command_argument_count() >= 3) then
        call get_command_argument(3, arg); read (arg, *) modes
    end if
    allocate (rings(nr)); rings = [(real(i, dp)/nr, i=1, nr)]
    ! Synthetic poloidal coils outside the torus, plus a remote tilted coil.
    do k = 1, npc
        phi = 2*pi*(k - 1)/npc
        do j = 1, 64
            theta = 2*pi*(j - 1)/64
            pts(:, 64*(k - 1) + j) = &
                [(r0 + 1.3_dp*cos(theta))*cos(phi), &
                 (r0 + 1.3_dp*cos(theta))*sin(phi), 1.3_dp*sin(theta)]
        end do
    end do
    do j = 1, 64
        theta = 2*pi*(j - 1)/64
        pts(:, 64*npc + j) = [7.0_dp + 0.2_dp*sin(theta), cos(theta), &
                              1.0_dp + sin(theta)]
    end do
    currents = 160000.0_dp/npc; currents(npc + 1) = 2.e4_dp; counts = 64
    call coil_set_from_polylines(coils, counts, pts, status=status)
    if (status /= 0) error stop 'synthetic coils'
    call flux_mesh_build(rings, 4, 64, 2, mapping, mesh, status)
    if (status /= 0) error stop 'mesh'
    call nn_build(mesh, 2, modes, [r0, 0.0_dp], minor, &
                  -coil_mu0*16*1.e4_dp/(2*pi), coil_mu0, field, status)
    if (status /= 0) error stop 'FEEC build'
    problem%source => source
    options%maximum_degree = 80
    call cpu_time(start)
    call nn_attach_exterior(field, problem, coils, currents, nphi, degree, &
                            options, exterior, report, status)
    print '(a,4i12,es16.7)', 'BEM status/failures/degree/pairs/ratio ', &
        status, report%failed_pairs, report%maximum_degree_used, &
        report%accepted_pairs, report%error_ratio
    if (status /= 0) error stop '3D exterior attachment'
    print '(a,3i12,i18,3es16.7)', 'BEM panels/trace/flux/bytes/assembly/DtN/coupling ', &
        exterior%panels, exterior%trace_dofs, exterior%flux_dofs, &
        exterior%operator_bytes, exterior%assembly_seconds, &
        exterior%dtn_seconds, exterior%coupling_seconds
    lay = nc_layout(field); fixed = nc_fixed(field, problem)
    if (count(fixed(:lay%na)) /= 1) error stop 'magnetic boundary still fixed'
    allocate (state(lay%n), direction(lay%n)); state = 0
    direction = [(1.e-4_dp*sin(real(i, dp)), i=1, lay%n)]
    where (fixed) direction = 0
    ! Nonzero harmonic B keeps the closure regular even at a zero A seed.
    call nc_jacobian(field, problem, state, jacobian, status)
    if (status /= 0) error stop 'Jacobian'

    step = 1.e-5_dp
    call nc_residual(field, problem, state + step*direction, rp, status)
    if (status /= 0) error stop 'plus residual'
    call nc_residual(field, problem, state - step*direction, rm, status)
    if (status /= 0) error stop 'minus residual'
    load = csc_matvec(jacobian, direction)
    where (fixed) load = 0
    jac_error = norm2((rp - rm)/(2*step) - load)/max(norm2(load), tiny(1.0_dp))
    print '(a,es16.7)', 'coupled Jacobian relative defect ', jac_error
    if (jac_error > 1.e-6_dp) error stop 'exterior Jacobian'
    call system_clock(solve_tick, solve_rate)
    call nc_solve(field, problem, state, newton, status)
    call system_clock(end_tick)
    solve_seconds = real(end_tick - solve_tick, dp)/solve_rate
    print '(a,es16.7)', 'Newton wall seconds ', solve_seconds
    print '(a,i4,6es16.7)', 'Newton status/residual/div ', status, &
        newton%residual, newton%divergence
    if (status /= 0) error stop 'free-boundary Newton'
    call nn_exterior_trace(exterior, state(:lay%na), potential, flux, load, &
        defect, status, physical_defect)
    if (status /= 0) error stop 'exterior trace'
    error = 0; scale = 0
    do i = 1, 24
        phi = 2*pi*(i - 0.37_dp)/24
        call mapping(0.7_dp, 2*pi*(i - 0.19_dp)/24, rz)
        call nn_point(field, [rz, phi], values, status)
        if (status /= 0) error stop 'field query'
        xq = [rz(1)*cos(phi), rz(1)*sin(phi), rz(2)]
        call coil_set_b_cart(coils, currents, xq, bc)
        direct = [bc(1)*cos(phi) + bc(2)*sin(phi), &
                  -bc(1)*sin(phi) + bc(2)*cos(phi), bc(3)]
        b = [values(EQ_BR), values(EQ_BPHI), values(EQ_BZ)]
        error = error + sum((b - direct)**2); scale = scale + sum(direct**2)
    end do
    error = sqrt(error/scale)
    print '(a,es16.7)', 'physical normal-flux defect ', physical_defect
    saved = state
    direction(lay%ip + 1:) = 0
    state = state + direction
    call nc_solve(field, problem, state, newton, status)
    if (status /= 0) error stop 'displaced seed Newton'
    seed_error = norm2(state(:lay%na) - saved(:lay%na)) &
                 /max(norm2(saved(:lay%na)), tiny(1.0_dp))
    call cpu_time(finish)
    print '(a,4es16.7)', 'field/weak-interface/seed/cpu ', &
        error, defect, seed_error, finish - start
    if (defect > 1.e-6_dp) error stop 'weak interface'
    if (seed_error > 1.e-3_dp) error stop 'coefficient seed return'
    ! This is a pilot, not G8: virtual casing, continuous traces and a
    ! geometrically displaced seed need independent checks at refined sizes.
    if (error > 1.e-3_dp) error stop 'direct-coil field gate'
end program
