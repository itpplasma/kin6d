program test_gs_prescribed_q_intervals
    ! Prescribed-q inverse with many FF' spline intervals: arbitrary (graded)
    ! knots, equation-level Jacobian and the secant fallback.
    !
    ! The manufactured circular solution is that of the variable-pressure test:
    ! psi_ex = axis*(1 - rho^2) with cubic p'(s), FF'(s) and the matching
    ! forcing; q_ex and Phi_ex follow from independent quadrature of psi_ex.
    ! A cubic FF' lies in every cubic B-spline space, so refining or grading
    ! the knots must keep the discretization orders.
    !   mms K G [R [M]]  orders on M (default 3) meshes; K intervals, G = 0 uniform,
    !               1 graded; R radial q levels per interval (0: 2n+1 per mesh)
    !   jacobian    exact Jacobian against Richardson differences of converged
    !               forward solves, on a small case
    !   secant      exact-Jacobian and secant solutions agree
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
    real(dp) :: s(513), qtarget(513), phi
    character(16) :: mode, text
    integer :: k, graded, per_interval = 0, meshes = 3
    do k = 1, size(s)
        s(k) = real(k - 1, dp)/(size(s) - 1)
        qtarget(k) = exact_q(s(k))
    end do
    phi = exact_phi()
    call get_command_argument(1, mode)
    select case (trim(mode))
    case ('mms')
        call get_command_argument(2, text); read (text, *) k
        call get_command_argument(3, text); read (text, *) graded
        if (command_argument_count() > 3) then
            call get_command_argument(4, text); read (text, *) per_interval
        end if
        if (command_argument_count() > 4) then
            call get_command_argument(5, text); read (text, *) meshes
        end if
        call check_mms(k, graded == 1, per_interval, meshes)
    case ('jacobian')
        call check_jacobian()
    case ('secant')
        call check_secant()
    case default
        error stop 'mode: mms K G | jacobian | secant'
    end select
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


    real(dp) function exact_phi() result(value)
        real(dp) :: nodes(64), weights(64)
        integer :: j
        call gauss_legendre_ab(64, 0.0_dp, 1.0_dp, nodes, weights)
        value = pi*a*a*sum(weights*[(exact_f(nodes(j))/sqrt(r0*r0 - a*a*nodes(j)), &
                                     j=1, size(nodes))])
    end function

    function graded_knots(intervals) result(knots)
        ! Knots refined (at most 2:1) toward the axis and a pedestal at s = 0.94 by
        ! equidistributing a monitor function, as in the TC24 knot study.
        integer, intent(in) :: intervals
        real(dp) :: knots(intervals + 1)
        real(dp), allocatable :: x(:), m(:), cumulative(:)
        integer :: i, j, fine
        fine = 16*intervals
        allocate (x(fine + 1), m(fine + 1), cumulative(fine + 1))
        x = [(real(i, dp)/fine, i=0, fine)]
        m = 1 + (exp(-x/0.01_dp) + exp(-((x - 0.94_dp)/0.02_dp)**2) &
                    + exp(-(1 - x)/0.02_dp))
        cumulative(1) = 0
        do i = 2, fine + 1
            cumulative(i) = cumulative(i - 1) + 0.5_dp*(m(i) + m(i - 1))*(x(i) - x(i - 1))
        end do
        cumulative = cumulative/cumulative(fine + 1)
        knots(1) = 0; knots(intervals + 1) = 1
        j = 1
        do i = 2, intervals
            do while (cumulative(j + 1) < real(i - 1, dp)/intervals)
                j = j + 1
            end do
            knots(i) = x(j) + (x(j + 1) - x(j))*(real(i - 1, dp)/intervals &
                       - cumulative(j))/(cumulative(j + 1) - cumulative(j))
        end do
    end function

    subroutine check_mms(intervals, use_graded, per_interval, meshes)
        integer, intent(in) :: intervals, per_interval, meshes
        logical, intent(in) :: use_graded
        type(gs_state_t) :: state
        type(gs_profile_law_t) :: law, pressure
        real(dp), allocatable :: x(:), y(:), w(:), knots(:)
        real(dp) :: errors(3, 3), rates(3), point(2), det, shape(10)
        real(dp) :: gradients(2, 10), hessians(2, 2, 10), psi, qerror
        real(dp) :: nodes(64), weights(64), pp, ff, ppref, ffref
        integer :: mesh, n, j, e, q, status, evaluations, iterations, levels
        call triangle_duffy_quadrature(8, x, y, w, status)
        if (status /= 0) error stop 'quadrature'
        call gauss_legendre_ab(64, 0.0_dp, 1.0_dp, nodes, weights)
        pressure = pressure_cubic([0.0_dp, 1.0_dp])
        if (use_graded) then
            knots = graded_knots(intervals)
        else
            knots = [(real(j, dp)/intervals, j=0, intervals)]
        end if
        print '(a,i0,a,l1)', 'FF prime intervals', intervals, ' graded ', use_graded
        print '(a)', 'n dof psi_L2 FFprime_L2 Fedge_rel q_max newton forward_solves'
        do mesh = 1, meshes
            n = nint(24*1.5_dp**(mesh - 1))
            levels = 2*n + 1
            if (per_interval > 0) levels = per_interval*intervals + 1
            call circle_state(state, n)
            call solve_gs_prescribed_q(state, mu0, 0.0_dp, s, qtarget, phi, law, &
                status, knots=knots, radial_points=levels, q_error=qerror, &
                evaluations=evaluations, newton_iterations=iterations, &
                tolerance=1.0e-10_dp, pressure_law=pressure, forcing=forcing, &
                verbose=.true.)
            if (status /= 0) then
                print *, 'status', status
                error stop 'many-interval inverse solve'
            end if
            errors(:, mesh) = 0
            do e = 1, size(state%triangles, 2)
                do q = 1, size(w)
                    call gs_cell_data(state, e, [x(q), y(q)], point, det, &
                                      shape, gradients, hessians, status)
                    if (status /= 0) error stop 'physical cell'
                    psi = dot_product(shape, state%psi(state%triangles(:, e)))
                    errors(1, mesh) = errors(1, mesh) + w(q)*det &
                                      *(psi - axis*(1 - rho2(point(1), point(2))))**2
                end do
            end do
            do j = 1, size(nodes)
                call gs_law_sources(law, law%axis*(1 - nodes(j)), pp, ff)
                call exact_sources(nodes(j), ppref, ffref)
                errors(2, mesh) = errors(2, mesh) + weights(j)*(ff - ffref)**2
            end do
            errors(1:2, mesh) = sqrt(errors(1:2, mesh))
            errors(3, mesh) = abs(law%f_edge/f_edge - 1)
            print '(i4,1x,i7,4(1x,es13.5),2(1x,i5))', n, count(.not. state%boundary), &
                errors(:, mesh), qerror, iterations, evaluations
        end do
        if (meshes < 3) return
        rates = log(errors(:, 1)/errors(:, 3))/log(2.25_dp)
        print '(a,3f10.5)', 'overall rates psi FFprime Fedge:', rates
        ! The field keeps the P3 order. The pointwise FF' and F_edge are set by
        ! the q moments, which on a fixed number of knots finer than the mesh
        ! carry a discretization-noise floor: they must still decrease.
        if (rates(1) < 3.3_dp) error stop 'many-interval psi P3 rate'
        if (errors(2, 3) > 0.7_dp*errors(2, 1)) error stop 'FFprime not converging'
        if (errors(3, 3) > 0.7_dp*errors(3, 1)) error stop 'F_edge not converging'
    end subroutine

    subroutine probe_start(seed, st)
        ! The solver's own starting unknowns (flat FF' of the one-interval seed).
        real(dp), intent(out) :: seed(:)
        integer, intent(out) :: st
        type(gs_state_t) :: state
        type(gs_profile_law_t) :: law, pressure
        integer, parameter :: n = 24, intervals = 4
        pressure = pressure_cubic([0.0_dp, 1.0_dp])
        call circle_state(state, n)
        call solve_gs_prescribed_q(state, mu0, 0.0_dp, s, qtarget, phi, law, st, &
            intervals=intervals, radial_points=2*n + 1, tolerance=1.0e-10_dp, &
            pressure_law=pressure, forcing=forcing, probe_seed=seed)
    end subroutine

    subroutine probe(x, st, residual, jacobian)
        ! One converged forward solve at the outer unknowns x on a fresh mesh.
        real(dp), intent(in) :: x(:)
        integer, intent(out) :: st
        real(dp), intent(out), optional :: residual(:), jacobian(:, :)
        type(gs_state_t) :: state
        type(gs_profile_law_t) :: law, pressure
        integer, parameter :: n = 24, intervals = 4
        pressure = pressure_cubic([0.0_dp, 1.0_dp])
        call circle_state(state, n)
        if (present(jacobian)) then
            call solve_gs_prescribed_q(state, mu0, 0.0_dp, s, qtarget, phi, law, st, &
                intervals=intervals, radial_points=2*n + 1, tolerance=1.0e-10_dp, &
                pressure_law=pressure, forcing=forcing, probe_x=x, &
                probe_residual=residual, probe_jacobian=jacobian)
        else
            call solve_gs_prescribed_q(state, mu0, 0.0_dp, s, qtarget, phi, law, st, &
                intervals=intervals, radial_points=2*n + 1, tolerance=1.0e-10_dp, &
                pressure_law=pressure, forcing=forcing, probe_x=x, &
                probe_residual=residual)
        end if
    end subroutine

    subroutine check_jacobian()
        ! The exact Jacobian of the outer residual against Richardson-extrapolated
        ! central differences of converged forward solves (each evaluation
        ! solves the nonlinear weak equation from the seed again).
        integer, parameter :: m = 8
        real(dp) :: seed(m), xp(m), xs(m), r0(m), jac(m, m), fd(m, m)
        real(dp) :: rp(m), rm(m), d1(m), d2(m), h(2), scale, worst
        integer :: j, side, status
        call probe_start(seed, status)
        if (status /= 0) error stop 'probe seed'
        xp = seed + 0.08_dp*[(sin(1.3_dp*j), j=1, m)]
        xp(m) = 0.02_dp
        call probe(xp, status, residual=r0, jacobian=jac)
        if (status /= 0) error stop 'probe jacobian'
        h = [0.04_dp, 0.02_dp]
        do j = 1, m
            do side = 1, 2
                xs = xp; xs(j) = xp(j) + h(side)
                call probe(xs, status, residual=rp)
                if (status /= 0) error stop 'probe plus'
                xs = xp; xs(j) = xp(j) - h(side)
                call probe(xs, status, residual=rm)
                if (status /= 0) error stop 'probe minus'
                if (side == 1) d1 = (rp - rm)/(2*h(1))
                if (side == 2) d2 = (rp - rm)/(2*h(2))
            end do
            fd(:, j) = (4*d2 - d1)/3
        end do
        worst = 0
        do j = 1, m
            scale = maxval(abs(fd(:, j)))
            worst = max(worst, maxval(abs(jac(:, j) - fd(:, j)))/scale)
            print '(a,i2,a,es10.3,a,es10.3)', 'column ', j, ' max|FD| ', scale, &
                ' relative error ', maxval(abs(jac(:, j) - fd(:, j)))/scale
        end do
        print '(a,es10.3)', 'worst column-relative Jacobian error: ', worst
        if (worst > 2.0e-4_dp) error stop 'exact Jacobian differs from Richardson'
    end subroutine

    subroutine check_secant()
        ! Both Jacobians drive the same discrete residual to the same tolerance,
        ! so the recovered fields and laws coincide.
        type(gs_state_t) :: exact_state, secant_state
        type(gs_profile_law_t) :: exact_law, secant_law, pressure
        real(dp) :: nodes(64), weights(64), pp, ff, ppsecant, ffsecant, qe, qs
        real(dp) :: dpsi, dff, dedge
        integer :: n, j, status, ie, is, newton, fe, fs
        n = 24
        pressure = pressure_cubic([0.0_dp, 1.0_dp])
        call gauss_legendre_ab(64, 0.0_dp, 1.0_dp, nodes, weights)
        call circle_state(exact_state, n); call circle_state(secant_state, n)
        call solve_gs_prescribed_q(exact_state, mu0, 0.0_dp, s, qtarget, phi, &
            exact_law, status, intervals=4, radial_points=2*n + 1, q_error=qe, &
            evaluations=fe, newton_iterations=newton, tolerance=1.0e-10_dp, &
            pressure_law=pressure, forcing=forcing, jacobian='exact')
        if (status /= 0) error stop 'exact-Jacobian solve'
        call solve_gs_prescribed_q(secant_state, mu0, 0.0_dp, s, qtarget, phi, &
            secant_law, status, intervals=4, radial_points=2*n + 1, q_error=qs, &
            evaluations=fs, tolerance=1.0e-10_dp, pressure_law=pressure, &
            forcing=forcing, jacobian='secant')
        if (status /= 0) error stop 'secant solve'
        dpsi = maxval(abs(exact_state%psi - secant_state%psi))/abs(axis)
        dff = 0
        do j = 1, size(nodes)
            call gs_law_sources(exact_law, exact_law%axis*(1 - nodes(j)), pp, ff)
            call gs_law_sources(secant_law, secant_law%axis*(1 - nodes(j)), &
                                ppsecant, ffsecant)
            dff = max(dff, abs(ff - ffsecant)/20)
        end do
        dedge = abs(exact_law%f_edge/secant_law%f_edge - 1)
        print '(a,3es11.3,2(a,i0))', 'exact vs secant psi, FFprime, F_edge:', &
            dpsi, dff, dedge, ' forward solves exact ', fe, ' secant ', fs
        if (dpsi > 1.0e-7_dp .or. dff > 1.0e-6_dp .or. dedge > 1.0e-8_dp) &
            error stop 'exact and secant solutions differ'
    end subroutine
end program
