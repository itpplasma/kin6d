module ampere_sheet_fixture
    ! Unit ring sheet of the per-n Ampere problem (ring_sheet_load form).
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_mhd_response_n, only: ring_angles_t, ring_sheet_density
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp)
    integer, parameter :: mode_n = 2, mode_m = 3
    type(ring_angles_t) :: angles
    integer, allocatable :: sheet_cells(:)
    real(dp), allocatable :: sheet_ta(:), sheet_tb(:)
contains
    subroutine circle(u, theta, x)
        real(dp), intent(in) :: u, theta
        real(dp), intent(out) :: x(2)
        x = [3.0_dp + u*cos(theta), u*sin(theta)]
    end subroutine

    subroutine zero_current(x, j)
        real(dp), intent(in) :: x(2)
        complex(dp), intent(out) :: j(3)
        j = cmplx(0.0_dp*x(1), 0.0_dp, dp)
    end subroutine

    subroutine sheet(cell, edge, tau, g)
        ! Density at parameter tau of local edge 1 of a ring cell.
        integer, intent(in) :: cell, edge
        real(dp), intent(in) :: tau
        complex(dp), intent(out) :: g
        integer :: i
        g = (0.0_dp, 0.0_dp)
        if (edge /= 1) return
        i = findloc(sheet_cells, cell, dim=1)
        if (i == 0) return
        g = ring_sheet_density(angles, mode_n, mode_m, (1.0_dp, 0.0_dp), &
            sheet_ta(i) + tau*(sheet_tb(i) - sheet_ta(i)))
    end subroutine
end module

program test_mhd_ampere_sheet_error
    ! Coupled equilibrated bound of a ring-sheet Ampere solve
    ! (ampere_n_coupled_error: sheet jump plus the discrete Steklov trace):
    !  - the Galerkin solution's bound converges at the energy rate (P2
    !    Nedelec, h^2) under ring refinement, which fails for a wrong sheet or
    !    exterior sign;
    !  - for perturbed candidates the bound contains the discrete coupled
    !    energy error (a lower bound of the continuous one by Galerkin
    !    orthogonality) with an efficiency index near 1;
    !  - omitting the sheet loads leaves an O(1) bound;
    !  - claim metadata: practical-discrete-exterior, Steklov term zero by
    !    construction, exterior remainder unavailable.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use ampere_sheet_fixture
    use fortfem_meridian_feec_n, only: meridian_feec_t, meridian_feec_build, &
        meridian_point_basis_t, meridian_point_basis_init
    use kin6d_mhd_mesh, only: flux_mesh_t, flux_mesh_build
    use kin6d_mhd_ampere_n, only: ampere_n_t, ampere_n_error_t, ampere_n_build, &
        ampere_n_solve, ampere_n_error, ampere_n_coupled_error, ampere_n_free
    use kin6d_mhd_response_n, only: ring_angles_set, ring_sheet_load, &
        ring_sheet_edges
    implicit none
    real(dp) :: bound(2), field(2), theta(64)
    integer :: lev, i

    do i = 1, size(theta)
        theta(i) = 2*pi*(i - 1)/size(theta)
    end do
    call ring_angles_set(1.0_dp, theta, theta, angles)
    do lev = 1, 2
        call level(4*2**(lev - 1), bound(lev), field(lev))
    end do
    print '(a, f6.2)', 'observed bound order ', log(bound(1)/bound(2))/log(2.0_dp)
    ! Measured relative bounds 0.14, 0.046, 0.0135 for nr = 4, 8, 16
    ! (orders 1.61, 1.76 towards the P2 energy rate 2).
    call require(bound(2) < 0.45_dp*bound(1), 'bound converges at the energy rate')
    print '(a)', 'PASS'
contains
    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message
        if (.not. condition) then
            print '(a)', 'FAIL: '//message
            error stop 1
        end if
    end subroutine

    subroutine level(nr, galerkin_bound, field_norm)
        integer, intent(in) :: nr
        real(dp), intent(out) :: galerkin_bound, field_norm
        type(flux_mesh_t) :: fm
        type(meridian_feec_t) :: fs, fsd
        type(meridian_point_basis_t) :: pb
        type(ampere_n_t) :: amp
        type(ampere_n_error_t) :: est
        complex(dp), allocatable :: f(:), a(:), delta(:)
        real(dp) :: rings(nr), index_
        integer :: st, j, i

        rings = [(real(i, dp)/nr, i=1, nr)]
        call flux_mesh_build(rings, 8, 64, 2, circle, fm, st)
        call require(st == 0, 'mesh')
        call meridian_feec_build(fm%mm, 2, fs, st, quad_degree=10)
        call require(st == 0, 'primal FEEC')
        call meridian_feec_build(fm%mm, 3, fsd, st, quad_degree=10)
        call require(st == 0, 'potential FEEC')
        call meridian_point_basis_init(fs, pb, st)
        call ampere_n_build(fs, mode_n, amp, st)
        call require(st == 0, 'Ampere build')
        j = nr/2
        call ring_sheet_edges(fm, j, sheet_cells, sheet_ta, sheet_tb)
        allocate (f(amp%nn), a(amp%nn))
        f = (0.0_dp, 0.0_dp)
        call ring_sheet_load(fs, pb, fm, j, angles, mode_n, mode_m, &
            (1.0_dp, 0.0_dp), f)
        call ampere_n_solve(amp, f, a, st)
        call require(st == 0, 'sheet solve')

        est = amp%error
        call ampere_n_coupled_error(fs, fsd, amp, a, zero_current, est, st, &
            sheet_cells, [(1, i=1, size(sheet_cells))], sheet)
        call require(st == 0, 'coupled error')
        galerkin_bound = est%coupled_bound
        field_norm = est%coupled%field_norm
        print '(a, i3, 5(a, es10.3))', 'nr', nr, ' Galerkin bound', &
            galerkin_bound, ' |B_h|', field_norm, ' discrete', &
            amp%error%discrete_energy_error, ' load defect', &
            est%coupled%load_defect, ' mismatch', est%coupled%lift_mismatch
        call require(est%coupled_claim == 'practical-discrete-exterior', 'claim')
        call require(est%steklov_term == 0.0_dp, 'Steklov term zero')
        call require(est%steklov_status == 'zero-by-construction-discrete-S', &
            'Steklov status')
        call require(est%continuous_claim == 'unavailable-exterior-remainder', &
            'exact exterior remainder stays unavailable')
        call require(galerkin_bound < 0.2_dp*field_norm, 'small Galerkin bound')

        ! Perturbed candidates: bound >= discrete coupled energy error.
        allocate (delta(amp%nn))
        do i = 1, amp%nn
            delta(i) = 1.0e-3_dp*cmplx(sin(1.3_dp*i), cos(0.7_dp*i), dp) &
                *maxval(abs(a))
        end do
        call ampere_n_error(amp, f, a + delta, est, st)
        call require(st == 0, 'discrete error')
        call ampere_n_coupled_error(fs, fsd, amp, a + delta, zero_current, est, &
            st, sheet_cells, [(1, i=1, size(sheet_cells))], sheet)
        call require(st == 0, 'perturbed coupled error')
        index_ = est%coupled_bound/est%discrete_energy_error
        print '(a, 2es11.3, a, f7.3)', '  perturbed: bound, discrete error', &
            est%coupled_bound, est%discrete_energy_error, ' index', index_
        call require(est%coupled_bound >= est%discrete_energy_error, &
            'bound contains the discrete coupled energy error')
        ! Measured index 1.32 (nr = 4), 1.01 (nr = 8).
        call require(index_ < 1.5_dp, 'efficiency for a large perturbation')

        ! Without the sheet loads the jump is missing.
        call ampere_n_coupled_error(fs, fsd, amp, a, zero_current, est, st)
        print '(a, es11.3)', '  without sheet loads: bound', est%coupled_bound
        call require(st == 0 .and. est%coupled_bound > 0.2_dp*field_norm, &
            'missing sheet loads are detected')
        call ampere_n_free(amp)
    end subroutine
end program
