program test_gs_prescribed_q_edge_mode
    ! Edge FF' mode of the prescribed-q inverse, against a known truth.
    ! The q data come from a truth law (pedestal p', smooth quadratic FF') on
    ! a fine mesh; the inverse runs on a coarse mesh with FF' knots graded
    ! toward the edge. With exact q moments the coarse model must reproduce
    ! its own discretization error, and an alternating edge FF' mode that q
    ! barely sees absorbs it (FF' at the edge thousands of times the truth).
    ! The discrepancy-principle regularization (misfit = the truth law's own
    ! coarse-mesh misfit) recovers a smooth FF' close to the truth.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_grad_shafranov, only: gs_state_t, initialize_gs, solve_gs
    use kin6d_gs_prescribed_q, only: solve_gs_prescribed_q, gs_q_table_value
    use kin6d_gs_profile_law, only: gs_profile_law_t, gs_law_sources
    use kin6d_gs_nonlinear, only: solve_gs_cubic_profiles
    use kin6d_gs_readback, only: gs_readback_t, read_gs_equilibrium, &
                                 gs_readback_resolution
    implicit none
    real(dp), parameter :: r0 = 6.2_dp, a = r0/10, mu0 = 1.2566370614359173e-6_dp
    real(dp), parameter :: pi = acos(-1.0_dp), f0 = -8, pedestal = 2
    integer, parameter :: nf = 48, nc = 24, ns = 81
    real(dp), parameter :: knots(11) = [0.0_dp, 0.3_dp, 0.6_dp, 0.8_dp, &
        0.85_dp, 0.9_dp, 0.93_dp, 0.96_dp, 0.98_dp, 0.99_dp, 1.0_dp]
    type(gs_state_t) :: fine, coarse
    type(gs_profile_law_t) :: truth, law
    type(gs_readback_t) :: readback
    real(dp), allocatable :: s(:), q(:), rho(:)
    real(dp) :: p0, phi, pp, ff, delta, sw, sm, x, w, lambda, q_error, error_none
    real(dp) :: error_reg
    integer :: j, status, status_none, nr, na
    p0 = -2*5.3_dp/(1.5_dp*mu0*r0*r0)
    truth%knots = [0.0_dp, 0.8_dp, 0.9_dp, 1.0_dp]
    truth%normalized = .true.
    allocate (truth%pp(4, 3), truth%ff(4, 3))
    truth%pp = 0; truth%pp(1, :) = p0
    truth%pp(3, 2) = p0*pedestal*300; truth%pp(4, 2) = -p0*pedestal*2000
    truth%pp(1, 3) = p0*(1 + pedestal)
    do j = 1, 3
        x = truth%knots(j)
        truth%ff(:, j) = f0*[1 - 0.5_dp*x + 0.8_dp*x*x, -0.5_dp + 1.6_dp*x, &
                             0.8_dp, 0.0_dp]
    end do
    truth%axis = -1; truth%f_edge = 5.3_dp*r0
    ! Data: the truth law's q on the fine mesh.
    call forward(fine, nf, truth)
    rho = [(real(j - 1, dp)/(ns - 1), j=1, ns)]
    call read_gs_equilibrium(fine, truth%f_edge, 0.0_dp, ns, 0, readback, status, &
                             truth, rho_grid=rho)
    if (status /= 0) error stop 'truth readback'
    s = rho**2; q = readback%q; phi = readback%phi_edge
    ! Data/model floor: the truth law on the coarse mesh, weighted like the
    ! solver's samples (uniform in rho).
    law = truth
    call forward(coarse, nc, law)
    nr = 2*nc + 1; na = 0
    call gs_readback_resolution(coarse, nr, na)
    rho = [(real(j - 1, dp)/(nr - 1), j=1, nr)]
    call read_gs_equilibrium(coarse, law%f_edge, 0.0_dp, nr, na, readback, status, &
                             law, rho_grid=rho)
    if (status /= 0) error stop 'coarse truth readback'
    sw = 0; sm = 0
    do j = 1, nr
        x = rho(j)**2
        w = (rho(min(nr, j + 1)) - rho(max(1, j - 1)))/2
        sw = sw + w; sm = sm + w*(readback%q(j)/gs_q_table_value(s, q, x) - 1)**2
    end do
    delta = sqrt(sm/sw)
    print '(a,es10.3)', 'coarse-mesh truth rms q misfit ', delta
    call circle_state(coarse, nc)
    call solve_gs_prescribed_q(coarse, mu0, 0.0_dp, s, q, phi, law, status_none, &
        knots=knots, radial_points=2*nc + 1, tolerance=1.0e-9_dp, &
        q_tolerance=2.0e-2_dp, pressure_law=truth)
    error_none = huge(1.0_dp)
    if (status_none == 0 .or. status_none == 9) error_none = law_error()
    print '(a,i0,a,es10.3)', 'exact moments: status ', status_none, &
        ' max FFprime error ', error_none
    call circle_state(coarse, nc)
    call solve_gs_prescribed_q(coarse, mu0, 0.0_dp, s, q, phi, law, status, &
        knots=knots, radial_points=2*nc + 1, tolerance=1.0e-9_dp, &
        q_tolerance=2.0e-2_dp, pressure_law=truth, regularization='discrepancy', &
        q_accuracy=delta, regularization_used=lambda, q_error=q_error)
    error_reg = huge(1.0_dp)
    if (status == 0) error_reg = law_error()
    print '(a,i0,a,es10.3,a,es10.3,a,es10.3)', 'discrepancy: status ', status, &
        ' lambda ', lambda, ' q error ', q_error, ' max FFprime error ', error_reg
    if (error_none < 1) error stop 'edge mode not excited: the check lost its purpose'
    if (status /= 0) error stop 'regularized inverse failed'
    if (lambda <= 0) error stop 'no regularization weight reported'
    if (error_reg > 0.1_dp) error stop 'regularized FFprime far from the truth'
contains
    real(dp) function law_error() result(worst)
        ! Maximum |FF' - truth| over s, relative to max |truth FF'| = 1.3 |f0|.
        real(dp) :: t, ft
        integer :: k
        worst = 0
        do k = 0, 200
            t = k/200.0_dp
            call gs_law_sources(law, law%axis*(1 - t), pp, ff)
            ft = f0*(1 - 0.5_dp*t + 0.8_dp*t*t)
            worst = max(worst, abs(ff - ft)/(1.3_dp*abs(f0)))
        end do
    end function
    subroutine forward(state, n, source)
        type(gs_state_t), intent(out) :: state
        integer, intent(in) :: n
        type(gs_profile_law_t), intent(inout) :: source
        integer :: st
        call circle_state(state, n)
        call gs_law_sources(source, 0.5_dp*source%axis, pp, ff)
        call solve_gs(state, mu0, pp, ff, st)
        if (st /= 0) error stop 'linear seed'
        call solve_gs_cubic_profiles(state, mu0, source, st, 1.0e-11_dp, &
                                     quadrature_tolerance=1.0e-9_dp)
        if (st /= 0) error stop 'truth forward'
    end subroutine
    subroutine circle_state(state, n)
        type(gs_state_t), intent(out) :: state
        integer, intent(in) :: n
        real(dp) :: boundary(2, n)
        integer :: k, st
        do k = 1, n
            boundary(:, k) = [r0 + a*cos(2*pi*real(k - 1, dp)/n), &
                              a*sin(2*pi*real(k - 1, dp)/n)]
        end do
        call initialize_gs(state, boundary, (2*pi*a/n)**2/2, st, 3, curve=circle)
        if (st /= 0) error stop 'mesh'
    end subroutine
    subroutine circle(t, point)
        real(dp), intent(in) :: t
        real(dp), intent(out) :: point(2)
        point = [r0 + a*cos(2*pi*t), a*sin(2*pi*t)]
    end subroutine
end program
