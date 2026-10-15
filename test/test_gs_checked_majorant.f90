program test_gs_checked_majorant
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_interval, only: interval_t
    use kin6d_grad_shafranov, only: gs_state_t, initialize_gs, solve_gs
    use kin6d_gs_equilibrated, only: gs_equilibrated_majorant, &
                                     equilibrated_flux_estimate_t
    use kin6d_gs_quadrature_budget, only: gs_checked_equilibrated_majorant
    implicit none
    type(gs_state_t) :: state
    type(equilibrated_flux_estimate_t) :: estimate
    type(interval_t) :: bound
    real(dp) :: contour(2, 4), zero, sampled
    integer :: status, order, node

    contour = reshape([2.0_dp, -1.0_dp, 4.0_dp, -1.0_dp, &
                       4.0_dp, 1.0_dp, 2.0_dp, 1.0_dp], [2, 4])
    do order = 1, 2
        call initialize_gs(state, contour, 0.25_dp, status, order)
        if (status /= 0) error stop 'affine test mesh'
        call solve_gs(state, 1.0_dp, 0.1_dp, 0.2_dp, status)
        if (status /= 0) error stop 'constant-profile solve'
        call gs_equilibrated_majorant(state, state%psi, source, estimate, status)
        if (status /= 0) error stop 'candidate RTp reconstruction'
        sampled = estimate%bound
        call gs_checked_equilibrated_majorant(state, state%psi, estimate, &
                                              1.0_dp, 0.1_dp, 0.2_dp, bound, status, 2)
        if (status /= 0) error stop 'continuous interval audit'
        if (bound%lo < 0 .or. bound%hi >= 10.0_dp) error stop 'unusable bound'
        print '(a,i0,2es14.5)', 'P', order, sampled, bound%hi
        ! Exact zero solution independently exercises every accumulated
        ! quadrature and arithmetic budget without an assumed relative slack.
        state%psi = 0
        zero = 0
        estimate%flux_coefficients = 0
        call gs_checked_equilibrated_majorant(state, state%psi, estimate, &
                                              1.0_dp, zero, zero, bound, status)
        if (status /= 0) error stop 'zero audit'
        if (bound%lo > 0 .or. bound%hi > 1.e-100_dp) error stop 'zero containment'
        ! Nonzero boundary data and curved geometry are outside this theorem.
        do node = 1, size(state%psi)
            if (.not. state%boundary(node)) cycle
            state%psi(node) = 1
            exit
        end do
        call gs_checked_equilibrated_majorant(state, state%psi, estimate, &
                                              1.0_dp, zero, zero, bound, status)
        if (status == 0) error stop 'nonhomogeneous boundary was certified'
        state%psi = 0
        state%curved_cell(1) = .true.
        call gs_checked_equilibrated_majorant(state, state%psi, estimate, &
                                             1.0_dp, zero, zero, bound, status)
        if (status == 0) error stop 'curved cell was certified'
    end do
    print '(a)', 'PASS: affine GS recovered-flux continuous interval audit'
contains
    real(dp) function source(cell, r, z, psi) result(f)
        integer, intent(in) :: cell
        real(dp), intent(in) :: r, z, psi
        f = 0.1_dp*r + 0.2_dp/r + 0.0_dp*(cell + z + psi)
    end function source
end program test_gs_checked_majorant
