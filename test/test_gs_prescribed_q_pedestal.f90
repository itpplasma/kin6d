program test_gs_prescribed_q_pedestal
    ! Prescribed-q inverse (constant q = 1.5 on a circular boundary) with a
    ! strong pedestal in p'(s): p' = p0*(1 + pedestal*h(s)), h a C1 smooth step
    ! from 0 to 1 over s in [0.8, 0.9]. The seed field of the solver is built
    ! with the mean p'; for a strong pedestal it can have no readable q profile,
    ! and the outer Newton iteration then starts from a pressure homotopy.
    ! The solver's own independent q check (interleaved radial grid, doubled arc
    ! quadrature) and the delivered p' law are the oracles.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_grad_shafranov, only: gs_state_t, initialize_gs
    use kin6d_gs_prescribed_q, only: solve_gs_prescribed_q
    use kin6d_gs_profile_law, only: gs_profile_law_t, gs_law_sources
    implicit none
    real(dp), parameter :: r0 = 6.2_dp, a = r0/10, mu0 = 1.2566370614359173e-6_dp
    real(dp), parameter :: pedestal = 8.0_dp, pi = acos(-1.0_dp)
    type(gs_state_t) :: state
    type(gs_profile_law_t) :: law, pressure
    real(dp), allocatable :: boundary(:, :)
    real(dp) :: sq(4), qq(4), phi, p0, qerror, t, u, pp, ff
    integer :: n, j, status, evaluations
    sq = [0.0_dp, 1.0_dp/3, 2.0_dp/3, 1.0_dp]; qq = 1.5_dp
    p0 = -2*5.3_dp/(1.5_dp*mu0*r0*r0)
    phi = pi*a*a*5.3_dp
    pressure%knots = [0.0_dp, 0.8_dp, 0.9_dp, 1.0_dp]
    pressure%normalized = .true.
    allocate (pressure%pp(4, 3), pressure%ff(4, 3))
    pressure%ff = 0; pressure%pp = 0
    pressure%pp(1, :) = p0
    ! Step 3u^2 - 2u^3, u = (t - 0.8)/0.1, then 1.
    pressure%pp(3, 2) = p0*pedestal*300
    pressure%pp(4, 2) = -p0*pedestal*2000
    pressure%pp(1, 3) = p0*(1 + pedestal)
    n = 24
    allocate (boundary(2, n))
    do j = 1, n
        boundary(:, j) = [r0 + a*cos(2*pi*real(j - 1, dp)/n), &
                          a*sin(2*pi*real(j - 1, dp)/n)]
    end do
    call initialize_gs(state, boundary, (2*pi*a/n)**2/2, status, 3, curve=circle)
    if (status /= 0) error stop 'mesh initialization'
    call solve_gs_prescribed_q(state, mu0, 0.0_dp, sq, qq, phi, law, status, &
                               intervals=6, radial_points=2*n + 1, q_error=qerror, &
                               evaluations=evaluations, tolerance=1.0e-9_dp, q_tolerance=1.0e-2_dp, &
                               pressure_law=pressure)
    print '(a,i0,a,es10.2,a,i0)', ' pedestal inverse status ', status, ' q error ', &
        qerror, ' evaluations ', evaluations
    if (status /= 0) error stop 'pedestal inverse solve'
    if (qerror > 1.0e-2_dp) error stop 'pedestal inverse q'
    do j = 1, 20
        t = (j - 0.5_dp)/20
        u = min(1.0_dp, max(0.0_dp, (t - 0.8_dp)/0.1_dp))
        call gs_law_sources(law, law%axis*(1 - t), pp, ff)
        if (abs(pp - p0*(1 + pedestal*u*u*(3 - 2*u))) > 1.0e-9_dp*abs(p0)) &
            error stop 'prescribed pressure law altered'
    end do
contains
    subroutine circle(t, point)
        real(dp), intent(in) :: t
        real(dp), intent(out) :: point(2)
        point = [r0 + a*cos(2*pi*t), a*sin(2*pi*t)]
    end subroutine
end program
