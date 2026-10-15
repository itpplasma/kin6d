program test_gs_bem_residual
    !! The same curve, continuous P2 trace and enriched flux test space:
    !! quantify the part of a coarse flux error visible to the fine V matrix.
    !! This checks an algebraic identity, NOT a continuous exterior bound.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortfem_axisymmetric_gs_boundary_operators_2d, only: gs_bem_curve_t, &
        gs_bem_operators_t, gs_bem_curve_from_nodes, gs_bem_assemble, &
        gs_bem_exterior_neumann
    use fortfem_bem_residual_energy, only: bem_residual_factor, &
        bem_residual_energy, BEM_RESIDUAL_OK
    use kin6d_gs_majorant, only: gs_coupled_majorant_t
    implicit none
    type(gs_bem_curve_t) :: curve
    type(gs_bem_operators_t) :: coarse, fine
    type(gs_coupled_majorant_t) :: estimate
    real(dp), allocatable :: nodes(:, :), u(:), tc(:), tf(:), injected(:)
    real(dp), allocatable :: factor(:, :), rhs(:), work(:), delta(:)
    real(dp) :: theta, indicator, exact, zero_indicator
    integer :: p, n, status
    real(dp), parameter :: pi = acos(-1.0_dp)

    n = 8
    allocate (nodes(2, 2*n), u(2*n))
    do p = 1, 2*n
        theta = pi*(p - 1)/n
        nodes(:, p) = [3.0_dp + 0.7_dp*cos(theta), sin(theta)]
        u(p) = 1 + 0.2_dp*cos(2*theta) + 0.1_dp*sin(3*theta)
    end do
    call gs_bem_curve_from_nodes(nodes, curve, status)
    call require(status == 0, 'curve')
    call gs_bem_assemble(curve, 0, coarse, status, hypersingular=.false.)
    call require(status == 0, 'coarse BEM')
    call gs_bem_assemble(curve, 2, fine, status, hypersingular=.false.)
    call require(status == 0, 'fine BEM')
    ! Match the symmetric operator used by gs_bem_exterior_neumann/v_solve.
    coarse%v = 0.5_dp*(coarse%v + transpose(coarse%v))
    fine%v = 0.5_dp*(fine%v + transpose(fine%v))
    call gs_bem_exterior_neumann(coarse, u, tc, status)
    call require(status == 0, 'coarse conormal solve')
    call gs_bem_exterior_neumann(fine, u, tf, status)
    call require(status == 0, 'fine conormal solve')
    allocate (injected(fine%ny), work(fine%ny))
    do p = 1, n
        injected(3*p - 2:3*p) = tc(p)
    end do
    rhs = matmul(fine%k - 0.5_dp*fine%m, u)
    call bem_residual_factor(fine%v, factor, status)
    call require(status == BEM_RESIDUAL_OK, 'factor')
    call bem_residual_energy(fine%v, factor, rhs, injected, work, indicator, status)
    call require(status == BEM_RESIDUAL_OK, 'enriched residual')
    delta = tf - injected
    exact = sqrt(dot_product(delta, matmul(fine%v, delta)))
    call require(indicator > 1.0e-3_dp, 'nontrivial enriched flux error')
    call require(abs(indicator - exact) < 1.0e-12_dp*exact, &
        'reference-space flux energy identity')
    call require(norm2(work - delta) < 1.0e-11_dp*norm2(delta), &
        'residual correction equals independently solved flux difference')
    call bem_residual_energy(fine%v, factor, rhs, tf, work, zero_indicator, status)
    call require(status == BEM_RESIDUAL_OK, 'reference flux residual')
    call require(zero_indicator < 1.0e-12_dp, 'exact reference recovery')
    call require(estimate%exterior_status == 'unavailable-bem-remainder', &
        'a finite residual norm does not supply a continuous remainder')
    print '(a,3es14.5)', 'reference error, residual norm, exact recovery ', &
        exact, indicator, zero_indicator
contains
    subroutine require(ok, message)
        logical, intent(in) :: ok
        character(*), intent(in) :: message
        if (.not. ok) then
            print '(a)', message
            error stop 'GS BEM residual test'
        end if
    end subroutine
end program test_gs_bem_residual
