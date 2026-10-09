program test_gs_readback_order
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_grad_shafranov, only: gs_state_t, initialize_gs, solve_gs
    use kin6d_gs_readback, only: gs_readback_t, read_gs_equilibrium, gs_readback_q
    implicit none
    real(dp), parameter :: r0 = 6.2_dp, a = 0.62_dp, f = 32.86_dp
    real(dp), parameter :: c = f/(12*r0**3)
    type(gs_state_t) :: state
    type(gs_readback_t) :: result, polynomial
    real(dp), allocatable :: contour(:, :)
    real(dp) :: h, s, angle, radius, exact_q, exact_phi, exact_volume, flux
    real(dp) :: q_error(3), volume_error(3), phi_error(3), axis_error(3)
    real(dp) :: label_error(3), off_grid_error(3), rates(4), query
    integer :: n, mesh, j, i, status, degree
    character(8) :: argument
    degree = 2
    if (command_argument_count() > 0) then
        call get_command_argument(1, argument)
        read (argument, *) degree
    end if
    ! An exported cubic q profile must survive off-grid readback. This oracle
    ! isolates the old linear-interpolation floor from the FE solve error.
    polynomial%s_tor = [0.0_dp, 0.1_dp, 0.4_dp, 0.7_dp, 1.0_dp]
    polynomial%q = 1 + polynomial%s_tor + polynomial%s_tor**2 + polynomial%s_tor**3
    do i = 1, 20
        s = real(i, dp)/21
        if (abs(gs_readback_q(polynomial, s) - (1 + s + s*s + s*s*s)) > 1.0e-13_dp) &
            error stop 'cubic q profile readback'
    end do
    print '(a)', 'n axis_error volume_rel Phi_rel q_rel_max s_tor_abs_max'
    do mesh = 1, 3
        n = 48*2**(mesh - 1)
        allocate (contour(2, n))
        do j = 1, n
            call boundary(real(j - 1, dp)/n, contour(:, j))
        end do
        h = 2*acos(-1.0_dp)*a/n
        call initialize_gs(state, contour, h*h/2, status, degree, curve=boundary)
        if (status /= 0) error stop 'readback mesh'
        call solve_gs(state, 1.0_dp, -16*c, 0.0_dp, status)
        if (status /= 0) error stop 'readback solve'
        call read_gs_equilibrium(state, f, 0.0_dp, 0, 0, result, status)
        if (status /= 0) then
            print *, 'readback status', status
            error stop 'curved readback'
        end if
        exact_volume = 0; exact_phi = 0
        do i = 1, 4096
            angle = 2*acos(-1.0_dp)*(i - 0.5_dp)/4096
            radius = sqrt(r0*r0 + 2*r0*a*cos(angle))
            exact_phi = exact_phi + sin(angle)**2/radius**3
            ! Green's theorem: volume = pi * integral R^2 dZ.
            exact_volume = exact_volume + radius*r0*a*cos(angle) &
                + r0*r0*a*a*sin(angle)**2/radius
        end do
        exact_phi = exact_phi/4096*(2*acos(-1.0_dp)*f*r0*r0*a*a)
        exact_volume = exact_volume/4096*2*acos(-1.0_dp)**2
        axis_error(mesh) = norm2(result%axis - [r0, 0.0_dp])
        volume_error(mesh) = abs(result%volume/exact_volume - 1)
        phi_error(mesh) = abs(result%phi_edge/exact_phi - 1)
        q_error(mesh) = 0; label_error(mesh) = 0; off_grid_error(mesh) = 0
        do j = 2, size(result%q)
            ! Compare at the actual returned physical psi, not a fitted axis label.
            s = 1 + result%psi(j)/(4*c*r0*r0*a*a)
            exact_q = 0; flux = 0
            do i = 1, 4096
                angle = 2*acos(-1.0_dp)*(i - 0.5_dp)/4096
                radius = sqrt(r0*r0 + 2*r0*a*sqrt(s)*cos(angle))
                exact_q = exact_q + 1/radius**3
                flux = flux + sin(angle)**2/radius**3
            end do
            exact_q = exact_q/4096*1.5_dp*r0**3
            flux = flux/4096*(2*acos(-1.0_dp)*f*r0*r0*a*a*s)
            label_error(mesh) = max(label_error(mesh), abs(result%s_tor(j) &
                - flux/exact_phi))
            if (s < 0.05_dp .or. s > 0.98_dp) cycle
            q_error(mesh) = max(q_error(mesh), &
                abs(gs_readback_q(result, flux/exact_phi)/exact_q - 1))
        end do
        do j = 1, 80
            s = 0.05_dp + 0.93_dp*(j - 1)/79
            exact_q = 0; flux = 0
            do i = 1, 4096
                angle = 2*acos(-1.0_dp)*(i - 0.5_dp)/4096
                radius = sqrt(r0*r0 + 2*r0*a*sqrt(s)*cos(angle))
                exact_q = exact_q + 1/radius**3
                flux = flux + sin(angle)**2/radius**3
            end do
            exact_q = exact_q/4096*1.5_dp*r0**3
            flux = flux/4096*(2*acos(-1.0_dp)*f*r0*r0*a*a*s)
            query = flux/exact_phi
            off_grid_error(mesh) = max(off_grid_error(mesh), &
                abs(gs_readback_q(result, query)/exact_q - 1))
        end do
        print '(i4,5es14.5)', n, axis_error(mesh), volume_error(mesh), &
            phi_error(mesh), q_error(mesh), label_error(mesh)
        print '(a,es14.5)', 'off-grid q error: ', off_grid_error(mesh)
        deallocate (contour)
    end do
    do mesh = 2, 3
        rates = log([axis_error(mesh - 1)/axis_error(mesh), &
            volume_error(mesh - 1)/volume_error(mesh), &
            phi_error(mesh - 1)/phi_error(mesh), &
            q_error(mesh - 1)/q_error(mesh)])/log(2.0_dp)
        print '(a,4f10.4)', 'axis, volume, Phi, q rates: ', rates
        if (rates(2) < 3.7_dp .or. rates(3) < 3.7_dp) error stop 'geometry measures'

    end do
    ! Fit the slope over all three meshes; independent triangulations have
    ! mesh-to-mesh variation in the maximum contour error.
    rates(4) = log(q_error(1)/q_error(3))/log(4.0_dp)
    print '(a,f10.4)', 'three-mesh q rate: ', rates(4)
    if (abs(rates(4) - degree) > 0.3_dp) error stop 'q field order'
    if (log(off_grid_error(1)/off_grid_error(3))/log(4.0_dp) < degree - 0.3_dp) &
        error stop 'off-grid q loses field order'
    if (log(label_error(1)/label_error(3))/log(4.0_dp) < degree - 0.3_dp) &
        error stop 'toroidal labels lose field order'
    if (label_error(3) >= label_error(2)) error stop 'toroidal flux refinement'
contains
    subroutine boundary(t, point)
        real(dp), intent(in) :: t
        real(dp), intent(out) :: point(2)
        real(dp) :: theta
        theta = 2*acos(-1.0_dp)*t
        point(1) = sqrt(r0*r0 + 2*r0*a*cos(theta))
        point(2) = r0*a*sin(theta)/point(1)
    end subroutine
end program
