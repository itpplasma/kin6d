module island_a_fixture
    ! Closure A+D provider case: the NN1 helical-coil vacuum field on the
    ! annulus u in [0.3, 1.2] about the coil (R = 6, Z = 0), the NN2 shell
    ! heat source of test_nonnested_island (u in [0.75, 1.05], p = 0 on
    ! both walls) and a floating inner conductor. The coil's helical
    ! excursion (0.02 m) supplies the resonant n >= 1 boundary trace; it is
    ! split off as the provider's epsilon direction, so epsilon = 1 is the
    ! NN2 island case (island width ~0.19 m at Z = 0) and epsilon = 0 its
    ! axisymmetric, nested counterpart with the same mode-0 data.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use nn_vacuum_fixture, only: r0, ftor, field, active_mode, coil_mu0, &
        make_coil, mean_psi, gauge_potential, annulus_mesh
    use kin6d_mhd_mesh, only: flux_mesh_t
    use fortfem_meridian_feec_n, only: meridian_lagrange_nodes
    use kin6d_equilibrium_nonnested, only: nn_report_t, nn_build, nn_solve, &
        nn_vacuum_mode
    use kin6d_nonnested_coupled, only: nc_problem_t, nc_layout_t, nc_layout
    implicit none
    private
    public :: island_a_setup, shell_source
    real(dp), parameter :: s0 = 3.0e6_dp
contains
    function shell_source(x) result(s)
        real(dp), intent(in) :: x(2)
        real(dp) :: s, u
        u = norm2(x - [r0, 0.0_dp])
        s = s0*max(0.0_dp, 1 - ((u - 0.9_dp)/0.15_dp)**2)**2
    end function

    subroutine island_a_setup(nr, k, modes, base, perturbation, problem)
        ! Builds nn_vacuum_fixture's field. base holds the mode-0 vacuum
        ! state (the epsilon = 0 seed and boundary data), perturbation the
        ! n >= 1 vacuum A (mode n Nedelec A_pol, 2 A_n = c - i s).
        integer, intent(in) :: nr, k, modes
        real(dp), allocatable, intent(out) :: base(:), perturbation(:)
        type(nc_problem_t), intent(out) :: problem
        type(flux_mesh_t) :: fm
        type(nn_report_t) :: report
        type(nc_layout_t) :: lay
        real(dp), allocatable :: nodes(:, :)
        complex(dp), allocatable :: bc(:), seed(:), a(:)
        integer :: st, n, j
        call make_coil()
        call annulus_mesh(nr, k, fm)
        call nn_build(fm, k, modes, [r0, 0.0_dp], 1.0_dp, ftor, coil_mu0, field, st)
        if (st /= 0) error stop 'island_a_setup: build'
        call meridian_lagrange_nodes(field%fs, nodes)
        lay = nc_layout(field)
        allocate (bc(lay%n0), seed(lay%n0), base(lay%n), perturbation(lay%n))
        bc = 0; seed = 0; base = 0; perturbation = 0
        !$omp parallel do default(shared) private(j) schedule(static)
        do j = 1, lay%n0
            if (field%fs%sp(0)%bnd(j)) bc(j) = mean_psi(nodes(:, j))
        end do
        !$omp end parallel do
        call nn_solve(field, 0, bc, seed, report, st, solution=a)
        if (st /= 0) error stop 'island_a_setup: vacuum n0'
        base(lay%ia + field%fs%nn + 1:lay%ia + lay%na1) = real(a, dp)
        do n = 1, modes
            active_mode = n
            call nn_vacuum_mode(field, n, gauge_potential, report, st, solution=a)
            if (st /= 0) error stop 'island_a_setup: vacuum mode'
            perturbation(lay%ia + (2*n - 1)*lay%na1 + 1:lay%ia + (2*n - 1)*lay%na1 &
                + field%fs%nn) = 2*real(a, dp)
            perturbation(lay%ia + 2*n*lay%na1 + 1:lay%ia + 2*n*lay%na1 &
                + field%fs%nn) = -2*aimag(a)
        end do
        problem%floating = field%fs%sp(0)%bnd .and. &
            norm2(nodes - spread([r0, 0.0_dp], 2, lay%n0), 1) < 0.75_dp
        problem%source => shell_source
        problem%chi_perp = 1; problem%chi_par = 1
        problem%g_par = 1; problem%g_perp = 1.0e-4_dp
        problem%current = 0
    end subroutine
end module
