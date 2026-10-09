program test_gs_prescribed_q
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use kin6d_grad_shafranov, only: gs_state_t, initialize_gs, gs_cell_data
    use kin6d_gs_prescribed_q, only: solve_gs_prescribed_q
    use kin6d_gs_profile_law, only: gs_profile_law_t, gs_law_primitives, gs_law_sources
    use fortfem_triangle_duffy_quadrature, only: triangle_duffy_quadrature
    use fortnum_quadrature, only: gauss_legendre_ab
    implicit none
    real(dp), parameter :: r0 = 6.2_dp, f = 32.86_dp, a = r0/3
    real(dp), parameter :: c = f/(12*r0**3), level = 4*r0*r0*a*a
    type(gs_state_t) :: state
    type(gs_profile_law_t) :: law
    real(dp), allocatable :: boundary(:, :), x(:), y(:), w(:)
    real(dp) :: s(513), qtarget(513), nodes(128), weights(128), phi
    real(dp) :: errors(4, 3), rates(4), point(2), det, shape(10)
    real(dp) :: gradients(2, 10), hessians(2, 2, 10), grad(2), exact_grad(2)
    real(dp) :: psi, exact, h, elapsed, fp, pp, ff, p, qerror
    integer :: mesh, n, j, e, k, status, evaluations
    integer(int64) :: start, finish, clock_rate
    call triangle_duffy_quadrature(8, x, y, w, status)
    if (status /= 0) error stop 'quadrature'
    do j = 1, size(s)
        s(j) = real(j - 1, dp)/(size(s) - 1)
        qtarget(j) = exact_q(s(j))
    end do
    call gauss_legendre_ab(128, 0.0_dp, 1.0_dp, nodes, weights)
    phi = 0
    do j = 1, size(nodes)
        phi = phi + weights(j)*exact_q(nodes(j))
    end do
    phi = phi*2*acos(-1.0_dp)*c*level
    print '(a)', 'n dof wall_s psi_L2 Bpol_L2 Btor_L2 FFprime_L2 q_max evaluations'
    do mesh = 1, 3
        n = 24*2**(mesh - 1)
        allocate (boundary(2, n))
        call system_clock(start, clock_rate)
        do j = 1, n
            call exact_boundary(real(j - 1, dp)/n, boundary(:, j))
        end do
        h = 2*acos(-1.0_dp)*a/n
        call initialize_gs(state, boundary, h*h/2, status, 3, curve=exact_boundary)
        if (status /= 0) error stop 'mesh initialization'
        call solve_gs_prescribed_q(state, 1.0_dp, -16*c, s, qtarget, phi, law, &
                           status, radial_points=2*n + 1, intervals=1, q_error=qerror, &
                                   evaluations=evaluations, tolerance=1.0e-10_dp)
        print *, 'inverse status:', status
        if (status /= 0) error stop 'inverse Solovev solve'
        call system_clock(finish)
        elapsed = real(finish - start, dp)/clock_rate
        if (elapsed > 300) error stop 'solve time limit'
        errors(:, mesh) = 0
        do e = 1, size(state%triangles, 2)
            do k = 1, size(w)
                call gs_cell_data(state, e, [x(k), y(k)], point, det, &
                                  shape, gradients, hessians, status)
                if (status /= 0) error stop 'physical cell'
                psi = dot_product(shape, state%psi(state%triangles(:, e)))
                grad = matmul(gradients, &
                    state%psi(state%triangles(:, e)) - state%psi(state%triangles(1, e)))
                exact = c*((point(1)**2 - r0*r0)**2 &
                           + 4*point(1)**2*point(2)**2 - level)
                exact_grad = [4*c*point(1)*(point(1)**2 - r0*r0 &
                                            + 2*point(2)**2), 8*c*point(1)**2*point(2)]
                call gs_law_primitives(law, psi, p, fp, status)
                if (status /= 0) error stop 'F recovery'
                call gs_law_sources(law, psi, pp, ff)
                errors(:, mesh) = errors(:, mesh) + w(k)*det &
                          *[(psi - exact)**2, sum((grad - exact_grad)**2)/point(1)**2, &
                                    ((fp - f)/point(1))**2, ff*ff]
            end do
        end do
        errors(:, mesh) = sqrt(errors(:, mesh))
        print '(i4,1x,i7,6(1x,es13.5),1x,i5)', n, &
            count(.not. state%boundary), elapsed, errors(:, mesh), qerror, evaluations
        deallocate (boundary)
    end do
    do mesh = 2, 3
        rates = log(errors(:, mesh - 1)/errors(:, mesh))/log(2.0_dp)
        print '(a,i1,4f10.5)', 'rates pair=', mesh - 1, rates
        if (abs(rates(1) - 4) > 0.5_dp) error stop 'inverse psi P3 rate'
        if (abs(rates(2) - 3) > 0.5_dp) error stop 'inverse Bpol P3 rate'
        if (rates(3) < 2.5_dp) error stop 'inverse F P3 convergence'
        if (rates(4) < 3.5_dp) error stop 'inverse FFprime coefficient convergence'
    end do
contains
    real(dp) function exact_q(s) result(q)
        real(dp), intent(in) :: s
        real(dp) :: theta, radius
        integer :: j
        q = 0
        do j = 1, 512
            theta = 2*acos(-1.0_dp)*(j - 1)/512
            radius = sqrt(r0*r0 + 2*r0*a*sqrt(s)*cos(theta))
            q = q + radius**(-3)
        end do
        q = q*f/(8*c*512)
    end function
    subroutine exact_boundary(t, point)
        real(dp), intent(in) :: t
        real(dp), intent(out) :: point(2)
        real(dp) :: angle
        angle = 2*acos(-1.0_dp)*t
        point(1) = sqrt(r0*r0 + 2*r0*a*cos(angle))
        point(2) = r0*a*sin(angle)/point(1)
    end subroutine
end program
