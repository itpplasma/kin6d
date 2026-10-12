program test_gs_prescribed_q_variable_p
    ! Prescribed-q inverse with a tabulated p'(s): manufactured circular solution.
    !
    ! psi_ex = axis*(1 - rho^2), rho^2 = ((R-r0)^2 + Z^2)/a^2, so the normalized
    ! flux is s = rho^2. The sources pp(s), ff(s) are cubic. The weak equation
    !   -div(grad psi/R) = ff/R + mu0*pp*R + forcing
    ! is satisfied by psi_ex for forcing = -Delta*psi_ex/R - ff/R - mu0*pp*R.
    ! On concentric circles dPhi/ds = pi a^2 F(s)/sqrt(r0^2 - a^2 s) with
    ! F^2 = F_edge^2 + 2*(-axis)*(PFF(s) - PFF(1)), so q(s) and Phi follow in
    ! closed form (Gauss-Legendre for Phi), independently of the solver.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_grad_shafranov, only: gs_state_t, initialize_gs, gs_cell_data
    use kin6d_gs_prescribed_q, only: solve_gs_prescribed_q
    use kin6d_gs_profile_law, only: gs_profile_law_t, gs_law_sources
    use fortfem_triangle_duffy_quadrature, only: triangle_duffy_quadrature
    use fortnum_quadrature, only: gauss_legendre_ab
    implicit none
    real(dp), parameter :: r0 = 6.2_dp, a = r0/10, mu0 = 1.2566370614359173e-6_dp
    real(dp), parameter :: axis = -0.6_dp, f_edge = 32.86_dp, pp_scale = -1.5e5_dp
    real(dp), parameter :: pi = acos(-1.0_dp)
    type(gs_state_t) :: state, other
    type(gs_profile_law_t) :: law, other_law, pressure, knotted
    real(dp), allocatable :: boundary(:, :), x(:), y(:), w(:)
    real(dp) :: s(513), qtarget(513), nodes(64), weights(64), phi
    real(dp) :: errors(4, 3), rates(4), point(2), det, shape(10)
    real(dp) :: gradients(2, 10), hessians(2, 2, 10), psi, h, qerror, pp, ff
    real(dp) :: ffref, ppref, diff
    integer :: mesh, n, j, e, k, status, evaluations
    call triangle_duffy_quadrature(8, x, y, w, status)
    if (status /= 0) error stop 'quadrature'
    do j = 1, size(s)
        s(j) = real(j - 1, dp)/(size(s) - 1)
        qtarget(j) = exact_q(s(j))
    end do
    call gauss_legendre_ab(64, 0.0_dp, 1.0_dp, nodes, weights)
    phi = pi*a*a*sum(weights*[(exact_f(nodes(j))/sqrt(r0*r0 - a*a*nodes(j)), &
                               j=1, size(nodes))])
    pressure = pressure_cubic([0.0_dp, 1.0_dp])
    print '(a)', 'n dof psi_L2 FFprime_L2 Fedge_rel q_max evaluations'
    do mesh = 1, 3
        n = nint(24*1.5_dp**(mesh - 1))
        allocate (boundary(2, n))
        do j = 1, n
            call circle(real(j - 1, dp)/n, boundary(:, j))
        end do
        h = 2*pi*a/n
        call initialize_gs(state, boundary, h*h/2, status, 3, curve=circle)
        if (status /= 0) error stop 'mesh initialization'
        call solve_gs_prescribed_q(state, mu0, 0.0_dp, s, qtarget, phi, law, status, &
                              intervals=1, radial_points=2*n + 1, q_error=qerror, &
                                   evaluations=evaluations, tolerance=1.0e-10_dp, &
                                   pressure_law=pressure, forcing=forcing)
        if (status /= 0) error stop 'variable-pressure inverse solve'
        errors(:, mesh) = 0
        do e = 1, size(state%triangles, 2)
            do k = 1, size(w)
                call gs_cell_data(state, e, [x(k), y(k)], point, det, &
                                  shape, gradients, hessians, status)
                if (status /= 0) error stop 'physical cell'
                psi = dot_product(shape, state%psi(state%triangles(:, e)))
                errors(1, mesh) = errors(1, mesh) + w(k)*det &
                                  *(psi - axis*(1 - rho2(point(1), point(2))))**2
            end do
        end do
        do j = 1, size(nodes)
            call gs_law_sources(law, law%axis*(1 - nodes(j)), pp, ff)
            call exact_sources(nodes(j), ppref, ffref)
            errors(2, mesh) = errors(2, mesh) + weights(j)*(ff - ffref)**2
            errors(4, mesh) = errors(4, mesh) + weights(j)*((pp - ppref)/pp_scale)**2
        end do
        errors(1:2, mesh) = sqrt(errors(1:2, mesh))
        errors(4, mesh) = sqrt(errors(4, mesh))
        errors(3, mesh) = abs(law%f_edge/f_edge - 1)
        print '(i4,1x,i7,4(1x,es13.5),1x,i5)', n, count(.not. state%boundary), &
            errors(1:3, mesh), qerror, evaluations
        if (errors(4, mesh) > 1.0e-12_dp) error stop 'prescribed pressure law altered'
        if (mesh == 1) then
            ! The same cubic written on extra pressure knots must give the same
            ! equilibrium: the union of FFprime and pressure knots is exact.
            call circle_state(other, n)
            knotted = pressure_cubic([0.0_dp, 0.3_dp, 1.0_dp])
            call solve_gs_prescribed_q(other, mu0, 0.0_dp, s, qtarget, phi, &
                                other_law, status, intervals=1, radial_points=2*n + 1, &
                                       tolerance=1.0e-10_dp, forcing=forcing, &
                                       pressure_law=knotted)
            if (status /= 0) error stop 'knotted pressure inverse solve'
            diff = 0
            do j = 1, size(nodes)
                call gs_law_sources(law, law%axis*(1 - nodes(j)), pp, ff)
                call gs_law_sources(other_law, other_law%axis*(1 - nodes(j)), &
                                    ppref, ffref)
                diff = max(diff, abs(ff - ffref), abs(pp - ppref)*mu0*r0*r0)
            end do
            print *, 'knotted vs single-interval pressure law, max source diff:', diff
            if (diff > 1.0e-6_dp*abs(ffref) + 1.0e-8_dp) &
                error stop 'pressure knots changed the solution'
            call check_scalar_compat()
            call check_pressure_guard()
        end if
        deallocate (boundary)
    end do
    do mesh = 2, 3
        rates = log(errors(:, mesh - 1)/errors(:, mesh))/log(1.5_dp)
        print '(a,i1,3f10.5)', 'rates psi FFprime Fedge pair=', mesh - 1, rates(1:3)
    end do
    rates(1:3) = log(errors(1:3, 1)/errors(1:3, 3))/log(2.25_dp)
    print '(a,3f10.5)', 'overall rates psi FFprime Fedge:', rates(1:3)
    if (rates(1) < 3.3_dp) error stop 'variable-p psi P3 rate'
    if (rates(2) < 3.5_dp) error stop 'variable-p FFprime rate'
    if (rates(3) < 3.5_dp) error stop 'variable-p F_edge rate'
contains
    subroutine circle(t, point)
        real(dp), intent(in) :: t
        real(dp), intent(out) :: point(2)
        point = [r0 + a*cos(2*pi*t), a*sin(2*pi*t)]
    end subroutine

    subroutine circle_state(st, n)
        ! Re-initialize a fresh mesh of the same size so no solved state is reused.
        type(gs_state_t), intent(inout) :: st
        integer, intent(in) :: n
        real(dp), allocatable :: b(:, :)
        integer :: i, stat
        allocate (b(2, n))
        do i = 1, n
            call circle(real(i - 1, dp)/n, b(:, i))
        end do
        call initialize_gs(st, b, (2*pi*a/n)**2/2, stat, 3, curve=circle)
        if (stat /= 0) error stop 'mesh initialization'
    end subroutine

    pure real(dp) function rho2(r, z) result(value)
        real(dp), intent(in) :: r, z
        value = ((r - r0)**2 + z*z)/(a*a)
    end function

    pure real(dp) function g_shape(t) result(value)
        real(dp), intent(in) :: t
        value = 1 + 0.4_dp*t - 0.8_dp*t*t + 0.3_dp*t**3
    end function

    pure subroutine exact_sources(t, pp, ff)
        real(dp), intent(in) :: t
        real(dp), intent(out) :: pp, ff
        pp = pp_scale*g_shape(t)
        ff = 20*(1 - 0.5_dp*t + 0.8_dp*t*t - 0.3_dp*t**3)
    end subroutine

    pure real(dp) function exact_pff(t) result(value)
        real(dp), intent(in) :: t
        value = 20*(t - 0.25_dp*t*t + 0.8_dp*t**3/3 - 0.075_dp*t**4)
    end function

    pure real(dp) function exact_f(t) result(value)
        real(dp), intent(in) :: t
        value = sqrt(f_edge**2 + 2*(-axis)*(exact_pff(t) - exact_pff(1.0_dp)))
    end function

    pure real(dp) function exact_q(t) result(value)
        real(dp), intent(in) :: t
        value = a*a*exact_f(t)/(2*(-axis)*sqrt(r0*r0 - a*a*t))
    end function

    function forcing(r, z) result(source)
        real(dp), intent(in) :: r, z
        real(dp) :: source, t, pp, ff, delta_star
        t = rho2(r, z)
        call exact_sources(t, pp, ff)
        delta_star = -4*axis/(a*a) + 2*axis*(r - r0)/(a*a*r)
        source = -delta_star/r - ff/r - mu0*pp*r
    end function

    function pressure_cubic(knots) result(law)
        ! pp(t) = pp_scale*g_shape(t) as exact Taylor re-expansions on each interval.
        real(dp), intent(in) :: knots(:)
        type(gs_profile_law_t) :: law
        integer :: i
        real(dp) :: t
        law%knots = knots; law%normalized = .true.
        allocate (law%pp(4, size(knots) - 1), law%ff(4, size(knots) - 1))
        law%ff = 0
        do i = 1, size(knots) - 1
            t = knots(i)
            law%pp(:, i) = pp_scale*[g_shape(t), 0.4_dp - 1.6_dp*t + 0.9_dp*t*t, &
                                     -0.8_dp + 0.9_dp*t, 0.3_dp]
        end do
    end function

    subroutine check_scalar_compat()
        ! A constant pressure law on [0,1] is the scalar path, bit for bit.
        type(gs_state_t) :: one, two
        type(gs_profile_law_t) :: l1, l2, constant
        real(dp) :: sq(4), qq(4), p0, ph
        integer :: st
        sq = [0.0_dp, 1.0_dp/3, 2.0_dp/3, 1.0_dp]; qq = 1.5_dp
        p0 = -2*5.3_dp/(1.5_dp*mu0*r0*r0)
        ph = pi*a*a*5.3_dp
        call circle_state(one, n); call circle_state(two, n)
        call solve_gs_prescribed_q(one, mu0, p0, sq, qq, ph, l1, st, intervals=1, &
                                   radial_points=2*n + 1, tolerance=1.0e-10_dp)
        if (st /= 0) error stop 'scalar inverse'
        constant%knots = [0.0_dp, 1.0_dp]; constant%normalized = .true.
        allocate (constant%pp(4, 1), constant%ff(4, 1))
        constant%pp = 0; constant%pp(1, 1) = p0; constant%ff = 0
        call solve_gs_prescribed_q(two, mu0, 0.0_dp, sq, qq, ph, l2, st, intervals=1, &
                      radial_points=2*n + 1, tolerance=1.0e-10_dp, &
                                   pressure_law=constant)
        if (st /= 0) error stop 'constant-law inverse'
        if (any(one%psi /= two%psi) .or. any(l1%ff /= l2%ff) .or. &
            any(l1%pp /= l2%pp) .or. l1%f_edge /= l2%f_edge) &
            error stop 'constant law differs from scalar pprime'
        ! Pressure stays nonnegative although p' > 0 on part of the interval.
        constant = pressure_cubic([0.0_dp, 1.0_dp])
        constant%pp = 0
        constant%pp(:, 1) = p0*[1.0_dp, -4.5_dp, 4.5_dp, 0.0_dp]
        call circle_state(two, n)
        call solve_gs_prescribed_q(two, mu0, 0.0_dp, sq, qq, ph, l2, st, intervals=1, &
                      radial_points=2*n + 1, tolerance=1.0e-10_dp, &
                                   pressure_law=constant)
        print *, 'nonmonotone-p'' (p >= 0) status:', st
        if (st /= 0) error stop 'nonnegative pressure with local p'' > 0 rejected'
    end subroutine

    subroutine check_pressure_guard()
        type(gs_state_t) :: trial
        type(gs_profile_law_t) :: bad, ok_law
        real(dp) :: sq(4), qq(4), ph
        integer :: st
        sq = [0.0_dp, 1.0_dp/3, 2.0_dp/3, 1.0_dp]; qq = 1.5_dp
        ph = pi*a*a*5.3_dp
        call circle_state(trial, n)
        ! p(s) = axis*integral_s^1 pp < 0 for pp = pp_scale*(1 - 2 s) near s = 0.
        bad = pressure_cubic([0.0_dp, 1.0_dp])
        bad%pp = 0; bad%pp(1:2, 1) = pp_scale*[1.0_dp, -2.0_dp]
        call solve_gs_prescribed_q(trial, mu0, 0.0_dp, sq, qq, ph, ok_law, st, &
                                   pressure_law=bad)
        if (st == 0) error stop 'negative pressure accepted'
        ! Positive pprime on the constant path is the same negative pressure.
        call solve_gs_prescribed_q(trial, mu0, 1.0e3_dp, sq, qq, ph, ok_law, st)
        if (st == 0) error stop 'positive scalar pprime accepted'
        ! Discontinuous and unnormalized laws are rejected.
        bad = pressure_cubic([0.0_dp, 0.5_dp, 1.0_dp])
        bad%pp(1, 2) = bad%pp(1, 2)*1.1_dp
        call solve_gs_prescribed_q(trial, mu0, 0.0_dp, sq, qq, ph, ok_law, st, &
                                   pressure_law=bad)
        if (st == 0) error stop 'discontinuous pressure law accepted'
        bad = pressure_cubic([0.0_dp, 1.0_dp])
        bad%normalized = .false.
        call solve_gs_prescribed_q(trial, mu0, 0.0_dp, sq, qq, ph, ok_law, st, &
                                   pressure_law=bad)
        if (st == 0) error stop 'unnormalized pressure law accepted'
        bad = pressure_cubic([0.0_dp, 0.8_dp])
        call solve_gs_prescribed_q(trial, mu0, 0.0_dp, sq, qq, ph, ok_law, st, &
                                   pressure_law=bad)
        if (st == 0) error stop 'pressure law not spanning 0..1 accepted'
    end subroutine
end program
