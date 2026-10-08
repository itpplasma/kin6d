program test_gs_p2_refinement
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_grad_shafranov, only: gs_state_t, initialize_gs, solve_gs, &
        sample_gs, gs_classical_current, gs_element_jet, write_gs_mesh
    implicit none
    type(gs_state_t) :: coarse, fine, affine
    real(dp) :: p2_coarse_field, p2_fine_field, p1_field, coarse_current, fine_current
    real(dp) :: affine_current
    integer :: status
    real(dp) :: boundary(2, 4)
    real(dp), parameter :: c = 0.002_dp, r0 = 6.2_dp
    boundary = reshape([5.8_dp, -0.4_dp, 6.6_dp, -0.4_dp, &
                         6.6_dp, 0.4_dp, 5.8_dp, 0.4_dp], [2, 4])
    call assess(coarse, 2, 0.0064_dp, p2_coarse_field, coarse_current)
    call assess(fine, 2, 0.0016_dp, p2_fine_field, fine_current)
    call assess(affine, 1, 0.0064_dp, p1_field, affine_current)
    print *, 'P2 coarse/fine Bpol RMS [T]:', p2_coarse_field, p2_fine_field
    print *, 'P2 coarse/fine classical current RMS [A/m2]:', coarse_current, fine_current
    print *, 'P1 same-mesh Bpol RMS [T]:', p1_field
    if (p2_fine_field >= 0.65_dp*p2_coarse_field) &
        error stop 'P2 independent physical field refinement failed'
    if (fine_current >= 0.85_dp*coarse_current) &
        error stop 'P2 independent classical current refinement failed'
    if (p2_coarse_field >= 0.35_dp*p1_field) &
        error stop 'P2 same-mesh physical accuracy improvement failed'
    call write_gs_mesh(coarse, 'p2_solovev_coarse_mesh.dat', status)
    if (status /= 0) error stop 'coarse solved P2 export failed'
    call write_gs_mesh(fine, 'p2_solovev_fine_mesh.dat', status)
    if (status /= 0) error stop 'fine solved P2 export failed'
    ! Deliberate unsupported geometry for the enclosure ingress control.
    coarse%points(1,coarse%triangles(4,1)) = &
        coarse%points(1,coarse%triangles(4,1)) + 0.125_dp
    call write_gs_mesh(coarse, 'p2_curved_metadata_mesh.dat', status)
    if (status /= 0) error stop 'invalid geometry fixture export failed'
    print *, 'PASS: native P2 independent polynomial refinement'
contains
    subroutine assess(state, order, max_area, field_error, current_error)
        type(gs_state_t), intent(out) :: state
        integer, intent(in) :: order
        real(dp), intent(in) :: max_area
        real(dp), intent(out) :: field_error, current_error
        real(dp) :: r, z, value, br, bz, bt, current, exact_br, exact_bz
        real(dp) :: grad(2), hessian(2, 2)
        integer :: i, j, e, status
        call initialize_gs(state, boundary, max_area, status, order)
        if (status /= 0) error stop 'P2 refinement mesh failed'
        state%psi = 0
        do i = 1, size(state%psi)
            if (.not. state%boundary(i)) cycle
            r = state%points(1, i); z = state%points(2, i)
            state%psi(i) = c*((r*r-r0*r0)**2+4*r*r*z*z)
        end do
        call solve_gs(state, 1.0_dp, -16*c, 0.0_dp, status)
        if (status /= 0) error stop 'P2 refinement manufactured GS solve failed'
        field_error = 0; current_error = 0
        do i = 1, 17
            r = 5.85_dp + 0.7_dp*i/18
            do j = 1, 19
                z = -0.35_dp + 0.7_dp*j/20
                call sample_gs(state, r, z, 30.0_dp, value, br, bz, bt, status)
                if (status /= 0) error stop 'P2 refinement physical query failed'
                exact_br = -8*c*r*z
                exact_bz = 4*c*(r*r-r0*r0+2*z*z)
                field_error = field_error + (br-exact_br)**2 + (bz-exact_bz)**2
                do e = 1, size(state%triangles, 2)
                    call gs_element_jet(state, e, r, z, value, grad, hessian, status)
                    if (status == 0) exit
                end do
                if (status /= 0) error stop 'P2 current query unavailable'
                call gs_classical_current(state, e, r, z, 1.0_dp, current, status)
                if (status /= 0) error stop 'P2 current sample failed'
                current_error = current_error + (current+16*c*r)**2
            end do
        end do
        field_error = sqrt(field_error/(17*19))
        current_error = sqrt(current_error/(17*19))
    end subroutine assess
end program test_gs_p2_refinement
