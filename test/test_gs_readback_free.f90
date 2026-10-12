program test_gs_readback_free
    !! Free-boundary readback (psi = psi_b on a plasma sub-mesh) against the
    !! fixed-boundary readback and the exact Solov'ev surfaces. Three checks:
    !!   1. A constant gauge shift with psi_bound = shift on the whole mesh
    !!      reproduces the fixed-boundary q, s_tor and toroidal flux.
    !!   2. An interior surface psi_b = psi(rho_k) taken as the boundary, with
    !!      the plasma restricted to the cells that reach it, gives at s = 1
    !!      the same q as the fixed-boundary readback at rho_k.
    !!   3. On that sub-mesh, q(psi) and Phi(psi) match the exact Solov'ev
    !!      surfaces of the enclosed region.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
    use kin6d_grad_shafranov, only: gs_state_t, initialize_gs, solve_gs
    use kin6d_gs_readback, only: gs_readback_t, read_gs_equilibrium
    implicit none
    real(dp), parameter :: r0 = 6.2_dp, a = 0.62_dp, f = 32.86_dp
    real(dp), parameter :: c = f/(12*r0**3), gauge = 0.37_dp
    type(gs_state_t) :: state, shifted
    type(gs_readback_t) :: fixed, whole, inner
    real(dp), allocatable :: contour(:, :)
    logical, allocatable :: plasma(:)
    real(dp) :: h, psi_b, s_b, phi_b, s, q_exact, phi_exact, q_error, phi_error
    real(dp) :: label_error, fixed_q, fixed_label
    integer :: n, j, k, e, status, degree
    character(8) :: argument
    degree = 2
    if (command_argument_count() > 0) then
        call get_command_argument(1, argument)
        read (argument, *) degree
    end if
    n = 96
    allocate (contour(2, n))
    do j = 1, n
        call boundary(real(j - 1, dp)/n, contour(:, j))
    end do
    h = 2*acos(-1.0_dp)*a/n
    call initialize_gs(state, contour, h*h/2, status, degree, curve=boundary)
    if (status /= 0) error stop 'mesh'
    call solve_gs(state, 1.0_dp, -16*c, 0.0_dp, status)
    if (status /= 0) error stop 'solve'
    call read_gs_equilibrium(state, f, 0.0_dp, 0, 0, fixed, status)
    if (status /= 0) error stop 'fixed readback'

    ! 1. Gauge shift, whole mesh.
    shifted = state
    shifted%psi = state%psi + gauge
    call read_gs_equilibrium(shifted, f, 0.0_dp, 0, 0, whole, status, &
                             psi_bound=gauge)
    if (status /= 0) error stop 'free readback, whole mesh'
    print '(a,3es12.4)', 'whole mesh: q, s_tor, psi differences ', &
        maxval(abs(whole%q/fixed%q - 1)), maxval(abs(whole%s_tor - fixed%s_tor)), &
        maxval(abs(whole%psi - gauge - fixed%psi))
    if (maxval(abs(whole%q/fixed%q - 1)) > 1.0e-11_dp) error stop 'q differs'
    if (maxval(abs(whole%s_tor - fixed%s_tor)) > 1.0e-11_dp) error stop 's_tor differs'
    if (maxval(abs(whole%psi - gauge - fixed%psi)) > 1.0e-12_dp) &
        error stop 'absolute flux levels differ'
    if (abs(whole%psi_axis - gauge - fixed%psi_axis) > 1.0e-12_dp) &
        error stop 'axis flux differs'
    if (norm2(whole%axis - fixed%axis) > 1.0e-10_dp) error stop 'axis differs'
    ! The fixed path reports the area sum; the free path the contour integral.
    if (abs(whole%phi_edge/fixed%phi_edge - 1) > fixed%flux_consistency + 1.0e-12_dp) &
        error stop 'toroidal flux differs'
    if (.not. ieee_is_nan(whole%volume)) error stop 'free volume must be NaN'

    ! 2. An interior surface as the boundary, plasma restricted to a sub-mesh.
    k = (size(fixed%psi) + 1)/2
    psi_b = fixed%psi(k)
    allocate (plasma(size(state%triangles, 2)))
    do e = 1, size(plasma)
        plasma(e) = minval(state%psi(state%triangles(:, e))) &
            <= psi_b - 0.05_dp*fixed%psi_axis
    end do
    if (all(plasma)) error stop 'sub-mesh is the whole mesh'
    call read_gs_equilibrium(shifted, f, 0.0_dp, 0, 0, inner, status, &
                             psi_bound=psi_b + gauge, plasma_cell=plasma)
    if (status /= 0) then
        print *, 'status', status
        error stop 'free readback, sub-mesh'
    end if
    print '(a,es12.4)', 'sub-mesh edge q against fixed q at the same level ', &
        abs(inner%q(size(inner%q))/fixed%q(k) - 1)
    if (abs(inner%q(size(inner%q))/fixed%q(k) - 1) > 1.0e-10_dp) &
        error stop 'edge q differs from the fixed readback'
    if (abs(inner%psi(size(inner%psi)) - gauge - psi_b) > 1.0e-12_dp) &
        error stop 'edge level is not psi_bound'

    ! 3. Exact surfaces of the enclosed region, against the error of the
    ! fixed-boundary readback on the same levels (s >= 0.05, as in
    ! test_gs_readback_order: the P2 axis cell limits q there).
    s_b = 1 + psi_b/(4*c*r0*r0*a*a)
    call exact(s_b, q_exact, phi_b)
    phi_error = abs(inner%phi_edge/phi_b - 1)
    q_error = 0; label_error = 0
    do j = 2, size(inner%q)
        s = 1 + (inner%psi(j) - gauge)/(4*c*r0*r0*a*a)
        call exact(s, q_exact, phi_exact)
        label_error = max(label_error, abs(inner%s_tor(j) - phi_exact/phi_b))
        if (s < 0.05_dp) cycle
        q_error = max(q_error, abs(inner%q(j)/q_exact - 1))
    end do
    fixed_q = 0; fixed_label = 0
    do j = 2, k
        s = 1 + fixed%psi(j)/(4*c*r0*r0*a*a)
        call exact(s, q_exact, phi_exact)
        fixed_label = max(fixed_label, abs(fixed%s_tor(j)*fixed%phi_edge &
            - phi_exact)/phi_b)
        if (s < 0.05_dp) cycle
        fixed_q = max(fixed_q, abs(fixed%q(j)/q_exact - 1))
    end do
    print '(a,f6.3,3es12.4)', 'sub-mesh s_b, q, Phi_edge, s_tor errors ', s_b, &
        q_error, phi_error, label_error
    print '(a,2es12.4)', 'fixed readback q, Phi errors on the same range ', &
        fixed_q, fixed_label
    if (q_error > 2*fixed_q + 1.0e-12_dp) error stop 'sub-mesh q'
    if (phi_error > 2*fixed_label + 1.0e-12_dp) error stop 'sub-mesh toroidal flux'
    if (label_error > 4*fixed_label + 1.0e-12_dp) error stop 'sub-mesh flux labels'
contains
    subroutine exact(s, q, phi)
        real(dp), intent(in) :: s
        real(dp), intent(out) :: q, phi
        real(dp) :: angle, radius
        integer :: i
        q = 0; phi = 0
        do i = 1, 4096
            angle = 2*acos(-1.0_dp)*(i - 0.5_dp)/4096
            radius = sqrt(r0*r0 + 2*r0*a*sqrt(s)*cos(angle))
            q = q + 1/radius**3
            phi = phi + sin(angle)**2/radius**3
        end do
        q = q/4096*1.5_dp*r0**3
        phi = phi/4096*(2*acos(-1.0_dp)*f*r0*r0*a*a*s)
    end subroutine

    subroutine boundary(t, point)
        real(dp), intent(in) :: t
        real(dp), intent(out) :: point(2)
        real(dp) :: theta
        theta = 2*acos(-1.0_dp)*t
        point(1) = sqrt(r0*r0 + 2*r0*a*cos(theta))
        point(2) = r0*a*sin(theta)/point(1)
    end subroutine
end program
