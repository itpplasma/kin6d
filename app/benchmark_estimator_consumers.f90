module estimator_consumer_data
    ! Constant heating is benchmark input; geometry uses the retained NN1 oracle.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortfem_meridian_fourier_scalar, only: meridian_fourier_source_t
    use kin6d_mhd_response_n, only: ring_angles_t, ring_sheet_density
    implicit none
    type, extends(meridian_fourier_source_t) :: heating_t
    contains
        procedure :: source => heating
    end type
    type(ring_angles_t) :: angles
    integer, allocatable :: sheet_cells(:)
    real(dp), allocatable :: sheet_ta(:), sheet_tb(:)
contains
    subroutine heating(self, t, x, f0, f, status)
        class(heating_t), intent(in) :: self
        integer, intent(in) :: t
        real(dp), intent(in) :: x(:, :)
        real(dp), intent(out) :: f0(:), f(:, :)
        integer, intent(out) :: status
        associate (unused => self)
            f0 = 1; f = 0
            status = 0
            if (t < 1 .or. size(x, 1) /= 3) status = 1
        end associate
    end subroutine

    subroutine zero_current(x, j)
        real(dp), intent(in) :: x(2)
        complex(dp), intent(out) :: j(3)
        j = cmplx(0.0_dp*x(1), 0.0_dp, dp)
    end subroutine

    subroutine sheet(cell, edge, tau, g)
        integer, intent(in) :: cell, edge
        real(dp), intent(in) :: tau
        complex(dp), intent(out) :: g
        integer :: i
        g = (0.0_dp, 0.0_dp)
        if (edge /= 1) return
        i = findloc(sheet_cells, cell, dim=1)
        if (i == 0) return
        g = ring_sheet_density(angles, 2, 3, (1.0_dp, 0.0_dp), &
            sheet_ta(i) + tau*(sheet_tb(i) - sheet_ta(i)))
    end subroutine
end module

program benchmark_estimator_consumers
    ! Wall-clock ledgers of real public consumers, without changing their mathematics.
    ! Usage: benchmark_estimator_consumers nn1|nn2|ampere|d3 SIZE REPEATS
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use estimator_consumer_data, only: heating_t, angles, sheet_cells, sheet_ta, &
        sheet_tb, sheet, zero_current
    use nn_solovev_fixture, only: r0, minor, ftor, pp, mapping, coordinates, boundary
    use kin6d_mhd_mesh, only: flux_mesh_t, flux_mesh_build
    use fortfem_meridian_feec_n, only: meridian_feec_t, meridian_feec_build, &
        meridian_point_basis_t, meridian_point_basis_init
    use kin6d_equilibrium_nonnested, only: free_b_field_t, nn_report_t, &
        nn_build, nn_solve
    use kin6d_nonnested_error, only: nn_error_t, nn_axisymmetric_error, &
        nt_axisymmetric_error
    use kin6d_nonnested_transport, only: nt_conductivity_t, nt_report_t, nt_solve
    use kin6d_mhd_ampere_n, only: ampere_n_t, ampere_n_error_t, ampere_n_build, &
        ampere_n_solve, ampere_n_coupled_error, ampere_n_free
    use kin6d_mhd_response_n, only: ring_angles_set, ring_sheet_load, ring_sheet_edges
    use kin6d_grad_shafranov, only: gs_state_t, initialize_gs, solve_gs
    use kin6d_gs_equilibrium, only: gs_equilibrium
    use kin6d_equilibrium_state, only: equilibrium_state_t
    use kin6d_equilibrium_residual, only: equilibrium_residual_t, equilibrium_residual
    implicit none
    type(flux_mesh_t) :: fm
    type(free_b_field_t), target :: field
    type(nn_report_t) :: nn_report
    type(nn_error_t) :: estimate
    type(nt_conductivity_t) :: conductivity
    type(nt_report_t) :: nt_report
    type(heating_t) :: source
    type(meridian_feec_t) :: fsd
    type(meridian_point_basis_t) :: pb
    type(ampere_n_t) :: amp
    type(ampere_n_error_t) :: amp_error
    type(gs_state_t) :: mesh
    type(equilibrium_state_t) :: eq
    type(equilibrium_residual_t) :: residual
    complex(dp), allocatable :: bc(:), seed(:), solution(:), load(:), a(:)
    real(dp), allocatable :: rings(:), p(:), transport_bc(:), contour(:, :)
    real(dp) :: theta(64), clock, solve_s, estimate_s, value, reference, h, reused_s
    integer :: nr, repeats, rep, i, st, ndof
    character(32) :: consumer, arg

    call get_command_argument(1, consumer)
    call get_command_argument(2, arg); read (arg, *) nr
    call get_command_argument(3, arg); read (arg, *) repeats
    if (nr < 2 .or. repeats < 1) error stop 'SIZE >= 2 and REPEATS >= 1 required'
    if (consumer == 'd3') then
        allocate (contour(2, nr))
        do i = 1, nr
            call boundary(real(i - 1, dp)/nr, contour(:, i))
        end do
        h = 2*acos(-1.0_dp)*minor/nr
        call initialize_gs(mesh, contour, h*h/2, st, 2, curve=boundary)
        call check(st, 'GS mesh')
        call solve_gs(mesh, 1.0_dp, pp, 0.0_dp, st)
        call check(st, 'GS solve')
        call gs_equilibrium(mesh, eq, st, ftor, 0.0_dp, pprime=pp, mu0=1.0_dp)
        call check(st, 'GS field adapter')
        ndof = 2*size(mesh%points, 2)
    else
        allocate (rings(nr))
        rings = [(real(i, dp)/nr, i=1, nr)]
        call flux_mesh_build(rings, 8, 4096, 2, mapping, fm, st)
        call check(st, 'flux mesh')
        call nn_build(fm, 2, 0, [r0, 0.0_dp], minor, ftor, 1.0_dp, field, st)
        call check(st, 'NN build')
        field%coordinates => coordinates
        allocate (bc(field%fs%n0), seed(field%fs%n0))
        bc = 0; seed = 0
        ndof = field%fs%n0
        if (consumer == 'nn2') then
            call nn_solve(field, 0, bc, seed, nn_report, st, pprime=pp)
            call check(st, 'NN magnetic solve')
            conductivity%field => field
            conductivity%chi_par = 1.0e3_dp
            conductivity%chi_perp = 1
            allocate (transport_bc(field%fs%n0)); transport_bc = 0
            ! Exercise nt_solve's actual estimate_error route once, then time
            ! that same estimator directly, avoiding subtraction of solve times.
            call nt_solve(field%fs, 0, 1, conductivity, source, transport_bc, p, &
                nt_report, st, estimate_error=.true.)
            call check(st, 'transport with estimate')
            call check(nt_report%error_status, 'transport estimate')
            reference = nt_report%error%bound
        else if (consumer == 'ampere') then
            call meridian_feec_build(fm%mm, 3, fsd, st, quad_degree=10)
            call check(st, 'Ampere potential space')
            call meridian_point_basis_init(field%fs, pb, st)
            call check(st, 'Ampere point basis')
            theta = [(2*acos(-1.0_dp)*(i - 1)/size(theta), i=1, size(theta))]
            call ring_angles_set(1.0_dp, theta, theta, angles)
            call ring_sheet_edges(fm, nr/2, sheet_cells, sheet_ta, sheet_tb)
        else if (consumer /= 'nn1') then
            error stop 'consumer must be nn1, nn2, ampere or d3'
        end if
    end if
    print '(a)', '# consumer size dof repeat solve_s estimator_s bound'
    do rep = 0, repeats
        select case (consumer)
        case ('nn1')
            clock = wall()
            call nn_solve(field, 0, bc, seed, nn_report, st, pprime=pp, &
                solution=solution)
            solve_s = wall() - clock
            call check(st, 'NN1 solve')
            clock = wall()
            call nn_axisymmetric_error(field, solution, estimate, st)
            estimate_s = wall() - clock
            call check(st, 'NN1 estimate')
            value = estimate%bound
        case ('nn2')
            clock = wall()
            call nt_solve(field%fs, 0, 1, conductivity, source, transport_bc, p, &
                nt_report, st)
            solve_s = wall() - clock
            call check(st, 'NN2 solve')
            clock = wall()
            call nt_axisymmetric_error(field%fs, 0, 1, conductivity, &
                conductivity%chi_perp, source, p, estimate, st)
            estimate_s = wall() - clock
            call check(st, 'NN2 estimate')
            value = estimate%bound
            if (abs(value - reference) > 1.0e-9_dp*reference) &
                error stop 'nt_solve and direct estimate differ'
        case ('ampere')
            clock = wall()
            call ampere_n_build(field%fs, 2, amp, st)
            call check(st, 'Ampere build')
            if (allocated(load)) deallocate (load, a)
            allocate (load(amp%nn), a(amp%nn)); load = 0
            call ring_sheet_load(field%fs, pb, fm, nr/2, angles, 2, 3, &
                (1.0_dp, 0.0_dp), load)
            call ampere_n_solve(amp, load, a, st)
            solve_s = wall() - clock
            call check(st, 'Ampere solve')
            ndof = amp%nn
            amp_error = amp%error
            clock = wall()
            call ampere_n_coupled_error(field%fs, fsd, amp, a, zero_current, &
                amp_error, st, sheet_cells, [(1, i=1, size(sheet_cells))], sheet)
            estimate_s = wall() - clock
            call check(st, 'Ampere estimate')
            value = amp_error%coupled_bound
            ! gs_response_3d reuses its per-n factor for multiple sheet RHSs.
            ! Report that denominator as well as setup + the first solve.
            reference = norm2(abs(a))
            clock = wall()
            call ampere_n_solve(amp, load, a, st)
            reused_s = wall() - clock
            call check(st, 'reused Ampere solve')
            if (abs(norm2(abs(a)) - reference) > 1.0e-12_dp*reference) &
                error stop 'reused Ampere solve differs'
            call ampere_n_free(amp)
        case ('d3')
            clock = wall()
            call equilibrium_residual(eq, mesh, residual, st)
            solve_s = wall() - clock
            call check(st, 'D3 Riesz solve')
            clock = wall()
            call equilibrium_residual(eq, mesh, residual, st, error_estimate=.true.)
            estimate_s = wall() - clock - solve_s
            call check(st, 'D3 estimated Riesz solve')
            call check(residual%bound_status, 'D3 estimate')
            value = residual%dual_upper
        end select
        if (rep == 0) cycle ! untimed warm-up
        print '(a,3(1x,i0),3(1x,es24.16))', trim(consumer), nr, ndof, rep, &
            solve_s, estimate_s, value
        if (consumer == 'ampere') then
            print '(a,3(1x,i0),3(1x,es24.16))', 'ampere_reuse', nr, ndof, rep, &
                reused_s, estimate_s, value
        end if
    end do
contains
    real(dp) function wall() result(t)
        integer(int64) :: count, rate
        call system_clock(count, rate)
        t = real(count, dp)/real(rate, dp)
    end function

    subroutine check(status, what)
        integer, intent(in) :: status
        character(*), intent(in) :: what
        if (status == 0) return
        print '(a,1x,i0)', what, status
        error stop 1
    end subroutine
end program
