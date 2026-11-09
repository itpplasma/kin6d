program test_gs_symmetry
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use milca_grad_shafranov, only: gs_state_t, initialize_gs, solve_gs, &
        gs_cell_data, gs_element_jet, sample_gs, apply_gs, gs_profile_load
    use milca_gs_symmetry, only: gs_symmetry_t, GS_ROTATION, GS_TRANSLATION, &
        GS_SCREW, GS_CHART_RZ, GS_CHART_CARTESIAN, GS_UNSUPPORTED_SYMMETRY
    use kin6d_gs_symmetry_leaves, only: gs_screw_field, gs_screw_field_jvp, &
        gs_screw_field_vjp, gs_screw_source, gs_screw_source_jvp
    use milca_gs_profile_law, only: gs_profile_law_t
    use milca_gs_nonlinear, only: solve_gs_cubic_profiles, gs_cubic_load, &
        gs_profile_jacobian_t, apply_gs_profile_jacobian, gs_profile_response, &
        gs_normalize_current
    use milca_gs_field_readback, only: sample_gs_many
    use milca_gs_coupling, only: gs_coupling_t, gs_coupling_build
    use milca_gs_majorant, only: gs_functional_majorant
    use milca_gs_inverse, only: gs_inverse_result_t, gs_toroidal_flux
    use fortfem_triangle_duffy_quadrature, only: triangle_duffy_quadrature
    implicit none
    real(dp), parameter :: lambda = 1.3_dp
    real(dp) :: pitch = 0
    logical :: manufactured = .false.
    character(32) :: gate
    call get_command_argument(1, gate)
    select case (trim(gate))
    case ('rotation')
        call rotation()
    case ('translation')
        pitch = 0
        call fixed_source(GS_TRANSLATION)
        call convergence(GS_TRANSLATION)
    case ('screw')
        pitch = 0.7_dp
        call fixed_source(GS_SCREW)
        call check_ck()
        call convergence(GS_SCREW)
    case ('manufactured')
        pitch = -0.8_dp; manufactured = .true.; call convergence(GS_SCREW)
    case ('limit')
        call limit()
    case ('derivatives')
        pitch = 0.7_dp; manufactured = .true.; call derivatives()
    case ('guards')
        call guards()
    case default
        error stop 'unknown RM0 gate'
    end select
contains
    subroutine require(ok, label)
        logical, intent(in) :: ok
        character(*), intent(in) :: label
        if (.not. ok) then
            print '(a)', label
            error stop 'RM0 oracle failed'
        end if
    end subroutine

    subroutine mesh(state, kind, n, order)
        type(gs_state_t), intent(out) :: state
        integer, intent(in) :: kind, n, order
        type(gs_symmetry_t) :: symmetry
        real(dp) :: contour(2, 4*n), t
        integer :: i, status
        do i = 1, n
            t = real(i - 1, dp)/n
            contour(:, i) = [-0.5_dp + t, -0.5_dp]
            contour(:, n + i) = [0.5_dp, -0.5_dp + t]
            contour(:, 2*n + i) = [0.5_dp - t, 0.5_dp]
            contour(:, 3*n + i) = [-0.5_dp, 0.5_dp - t]
        end do
        symmetry = gs_symmetry_t(kind, GS_CHART_CARTESIAN, pitch)
        if (kind == GS_ROTATION) then
            symmetry = gs_symmetry_t()
            contour(1, :) = contour(1, :) + 2
        end if
        call initialize_gs(state, contour, 0.5_dp/n**2, status, &
            element_order=order, symmetry=symmetry)
        call require(status == 0, 'symmetry mesh')
    end subroutine

    subroutine profile(law)
        type(gs_profile_law_t), intent(out) :: law
        law%knots = [-2.0_dp, 2.0_dp]
        allocate (law%pp(4, 1), law%ff(4, 1))
        law%pp = 0; law%ff = 0
        law%f_edge = 1; law%p_edge = 2
        ! F=1+lambda*psi is positive throughout these cases.
        law%ff(1, 1) = lambda - 2*lambda**2
        law%ff(2, 1) = lambda**2
        if (manufactured) then
            law%pp(1, 1) = 0.1_dp - 2*0.03_dp
            law%pp(2, 1) = 0.03_dp
        end if
    end subroutine

    subroutine fixed_source(kind)
        integer, intent(in) :: kind
        type(gs_state_t) :: state
        real(dp) :: alpha, f, source, psi, b(3), expected(3), x, y
        integer :: i, status
        call mesh(state, kind, 8, 2)
        alpha = -0.25_dp; f = 1; source = 1
        if (kind == GS_SCREW) then
            alpha = 0.2_dp; f = 2*alpha/pitch; source = 0
        end if
        state%psi = 0
        do i = 1, size(state%psi)
            if (state%boundary(i)) &
                state%psi(i) = alpha*sum(state%points(:, i)**2)
        end do
        call solve_gs(state, 1.0_dp, source, 0.0_dp, status, f=f)
        call require(status == 0, 'fixed-source solve')
        call require(maxval(abs(state%psi &
            - alpha*sum(state%points**2, dim=1))) < 1.0e-8_dp, &
            'analytic quadratic fixed-source flux')
        x = 0.21_dp; y = -0.17_dp
        call sample_gs(state, x, y, f, psi, b(1), b(2), b(3), status)
        expected = [-y/2, x/2, f]
        if (kind == GS_SCREW) expected = [0.0_dp, 0.0_dp, f]
        call require(status == 0, 'fixed-source field')
        call require(norm2(b - expected) < 1.0e-7_dp, 'analytic fixed-source field')
    end subroutine

    subroutine exact(x, y, psi, b)
        real(dp), intent(in) :: x, y
        real(dp), intent(out) :: psi, b(3)
        real(dp) :: r, c, s, j0, j1, j2, jp, gamma, br, bt, bz, g(2), f
        r = sqrt(x*x + y*y)
        if (manufactured) then
            psi = 0.2_dp*cos(x)*cos(2*y) + 0.1_dp*x*y
            g = [-0.2_dp*sin(x)*cos(2*y) + 0.1_dp*y, &
                -0.4_dp*cos(x)*sin(2*y) + 0.1_dp*x]
            f = 1 + lambda*psi
            ! Independent physical reconstruction at z=0.
            b(1) = (g(2) + pitch**2*x*(-y*g(1) + x*g(2)) - pitch*y*f) &
                /(1 + pitch**2*r*r)
            b(2) = (-g(1) + pitch**2*y*(-y*g(1) + x*g(2)) + pitch*x*f) &
                /(1 + pitch**2*r*r)
            b(3) = (pitch*(x*g(1) + y*g(2)) + f)/(1 + pitch**2*r*r)
            return
        end if
        j0 = bessel_j0(lambda*r); j1 = bessel_j1(lambda*r)
        b = [0.0_dp, 0.0_dp, j0]
        if (r > 1.0e-12_dp) b(:2) = [-y, x]*j1/r
        psi = (j0 + pitch*r*j1 - 1)/lambda
        ! A genuine u=theta-k*z CK mode, added to keep F away from zero.
        gamma = sqrt(lambda**2 - pitch**2)
        if (r < 1.0e-12_dp) then
            b(2) = b(2) - 0.08_dp*gamma*(lambda - pitch)/2
            return
        end if
        c = x/r; s = y/r
        j1 = bessel_j1(gamma*r); j0 = bessel_j0(gamma*r)
        j2 = bessel_jn(2, gamma*r); jp = (j0 - j2)/2
        br = (-lambda*j1/r + pitch*gamma*jp)*s
        bt = (-lambda*gamma*jp + pitch*j1/r)*c
        bz = gamma**2*j1*c
        b = b + 0.08_dp*[br*c - bt*s, br*s + bt*c, bz]
        psi = psi + 0.08_dp*(lambda*j1 - pitch*gamma*r*jp)*c
    end subroutine

    real(dp) function forcing(x, y) result(value)
        real(dp), intent(in) :: x, y
        real(dp) :: psi, b(3), g(2), hxx, hxy, hyy, s, op, pp, ff, f
        call exact(x, y, psi, b)
        g = [-0.2_dp*sin(x)*cos(2*y) + 0.1_dp*y, &
            -0.4_dp*cos(x)*sin(2*y) + 0.1_dp*x]
        hxx = -0.2_dp*cos(x)*cos(2*y)
        hxy = 0.4_dp*sin(x)*sin(2*y) + 0.1_dp
        hyy = -0.8_dp*cos(x)*cos(2*y)
        s = 1 + pitch**2*(x*x + y*y)
        ! Expanded divergence is independent of the generated weak leaf.
        op = ((1 + pitch**2*y*y)*hxx - 2*pitch**2*x*y*hxy &
            + (1 + pitch**2*x*x)*hyy)/s &
            - pitch**2*(s + 2)*(x*g(1) + y*g(2))/s**2
        pp = 0.1_dp + 0.03_dp*psi
        f = 1 + lambda*psi; ff = lambda*f
        value = -op - pp - ff/s + 2*pitch*f/s**2
    end function

    subroutine solve_case(state, kind, n, order)
        type(gs_state_t), intent(out) :: state
        integer, intent(in) :: kind, n, order
        type(gs_profile_law_t) :: law
        real(dp) :: psi, b(3)
        integer :: i, status
        call mesh(state, kind, n, order)
        do i = 1, size(state%psi)
            call exact(state%points(1, i), state%points(2, i), psi, b)
            state%psi(i) = 0
            if (state%boundary(i)) state%psi(i) = psi
        end do
        call profile(law)
        if (manufactured) then
            call solve_gs_cubic_profiles(state, 1.0_dp, law, status, &
                relative_tolerance=1.0e-11_dp, forcing=forcing, &
                quadrature_tolerance=1.0e-11_dp)
        else
            call solve_gs_cubic_profiles(state, 1.0_dp, law, status, &
                relative_tolerance=1.0e-11_dp, quadrature_tolerance=1.0e-11_dp)
        end if
        if (status /= 0) print *, 'solve status, residual: ', status, state%residual
        call require(status == 0, 'symmetry solve')
        call field_consumers(state)
    end subroutine

    subroutine field_consumers(state)
        type(gs_state_t), intent(in) :: state
        real(dp) :: p(2, 3), values(4, 3), scalar(4), f
        integer :: i, status, statuses(3)
        p = reshape([-0.12_dp, 0.13_dp, 0.0_dp, 0.0_dp, 0.21_dp, -0.17_dp], [2, 3])
        f = 1.2_dp
        call sample_gs_many(state, p, f, values, statuses)
        call require(all(statuses == 0), 'batch fields across the cylinder axis')
        do i = 1, 3
            call sample_gs(state, p(1, i), p(2, i), f, scalar(1), &
                scalar(2), scalar(3), scalar(4), status)
            call require(status == 0, 'scalar field across the cylinder axis')
            call require(maxval(abs(values(:, i) - scalar)) < 1.0e-12_dp, &
                'scalar and batch physical fields')
        end do
    end subroutine

    subroutine errors(state, error)
        type(gs_state_t), intent(in) :: state
        real(dp), intent(out) :: error(2)
        real(dp), allocatable :: x(:), y(:), w(:)
        real(dp) :: p(2), det, shape(10), g(2, 10), h(2, 2, 10), psi, b(3), bh(3), ph
        integer :: e, q, nd, status
        error = 0; nd = size(state%triangles, 1)
        call triangle_duffy_quadrature(10, x, y, w, status)
        call require(status == 0, 'error quadrature')
        do e = 1, size(state%area)
            do q = 1, size(w)
                call gs_cell_data(state, e, [x(q), y(q)], p, det, shape(:nd), &
                    g(:, :nd), h(:, :, :nd), status)
                call require(status == 0, 'error geometry')
                call exact(p(1), p(2), psi, b)
                call sample_cell(state, e, p, ph, bh, status)
                call require(status == 0, 'field evaluation')
                error(1) = error(1) + w(q)*det*(ph - psi)**2
                error(2) = error(2) + w(q)*det*sum((bh - b)**2)
            end do
        end do
        error = sqrt(error)
    end subroutine

    subroutine sample_cell(state, e, p, psi, b, status)
        type(gs_state_t), intent(in) :: state
        integer, intent(in) :: e
        real(dp), intent(in) :: p(2)
        real(dp), intent(out) :: psi, b(3)
        integer, intent(out) :: status
        real(dp) :: g(2), h(2, 2)
        call gs_element_jet(state, e, p(1), p(2), psi, g, h, status)
        if (status /= 0) return
        call gs_screw_field(p(1), p(2), pitch, g(1), g(2), 1 + lambda*psi, b(1), b(2), b(3))
    end subroutine

    subroutine convergence(kind)
        integer, intent(in) :: kind
        type(gs_state_t) :: state
        real(dp) :: error(2, 3), rates(2), clock, finish
        integer :: order, level, n
        print '(a)', 'degree n dofs seconds psi_L2 B_L2 (synthetic SI units)'
        do order = 2, 3
            do level = 1, 3
                n = 4*2**level
                call cpu_time(clock)
                call solve_case(state, kind, n, order)
                call errors(state, error(:, level))
                call cpu_time(finish)
                print '(3i8,3es15.6)', order, n, size(state%psi), finish - clock, error(:, level)
                if (level > 1) then
                    rates = log(error(:, level - 1)/error(:, level))/log(2.0_dp)
                    print '(a,2f10.5)', 'psi,B orders: ', rates
                    call require(rates(1) > order + 0.65_dp, 'psi design order')
                    call require(rates(2) > order - 0.3_dp, 'B design order')
                end if
            end do
        end do
    end subroutine

    subroutine rotation()
        type(gs_state_t) :: a, b
        type(gs_symmetry_t) :: symmetry
        real(dp) :: contour(2, 4), psi, br, bz, bt
        integer :: status, i
        contour = reshape([1.0_dp, -0.5_dp, 2.0_dp, -0.5_dp, &
            2.0_dp, 0.5_dp, 1.0_dp, 0.5_dp], [2, 4])
        call initialize_gs(a, contour, 0.01_dp, status, element_order=2)
        call require(status == 0, 'default rotation')
        call initialize_gs(b, contour, 0.01_dp, status, element_order=2, symmetry=symmetry)
        call require(status == 0, 'explicit rotation')
        call require(maxval(abs(a%stiffness - b%stiffness)) <= tiny(1.0_dp), &
            'rotation stiffness regression')
        do i = 1, size(a%psi)
            a%psi(i) = a%points(1, i)**2
        end do
        b%psi = a%psi
        call solve_gs(a, 1.0_dp, 0.0_dp, 0.0_dp, status)
        call require(status == 0, 'rotation vacuum solve')
        call solve_gs(b, 1.0_dp, 0.0_dp, 0.0_dp, status)
        call require(status == 0, 'explicit rotation vacuum solve')
        call require(maxval(abs(a%psi - b%psi)) < 1.0e-13_dp, 'rotation solve regression')
        call sample_gs(a, 1.3_dp, 0.12_dp, 2.0_dp, psi, br, bz, bt, status)
        call require(status == 0, 'rotation sample')
        call require(maxval(abs([psi - 1.3_dp**2, br, bz - 2, bt - 2/1.3_dp])) < 2.0e-8_dp, &
            'independent rotation flux and field')
    end subroutine

    subroutine check_ck()
        real(dp) :: x, y, psi, b(3), plus(3), minus(3), dx(3), dy(3), dz(3), curl(3)
        real(dp), parameter :: step = 1.0e-5_dp
        x = 0.21_dp; y = -0.17_dp
        call exact(x, y, psi, b)
        call exact(x + step, y, psi, plus); call exact(x - step, y, psi, minus)
        dx = (plus - minus)/(2*step)
        call exact(x, y + step, psi, plus); call exact(x, y - step, psi, minus)
        dy = (plus - minus)/(2*step)
        ! Equivariance gives d_z B=k(e_z x B+y d_x B-x d_y B).
        dz = pitch*([-b(2), b(1), 0.0_dp] + y*dx - x*dy)
        curl = [dy(3) - dz(2), dz(1) - dx(3), dx(2) - dy(1)]
        call require(maxval(abs(curl - lambda*b)) < 1.0e-9_dp, 'independent CK curl')
        call require(abs(dx(1) + dy(2) + dz(3)) < 1.0e-9_dp, 'independent CK divergence')
    end subroutine

    subroutine limit()
        type(gs_state_t) :: translation, screw
        real(dp) :: ref, error, previous, b0(3), b(3), psi, field_error, previous_b
        integer :: i, status
        pitch = 0
        call solve_case(translation, GS_TRANSLATION, 12, 2)
        call sample_gs(translation, 0.21_dp, -0.17_dp, 1.0_dp, &
            psi, b0(1), b0(2), b0(3), status)
        call require(status == 0, 'translation limit field')
        previous = huge(1.0_dp)
        previous_b = huge(1.0_dp)
        do i = 1, 4
            pitch = 10.0_dp**(-i - 1)
            call solve_case(screw, GS_SCREW, 12, 2)
            error = norm2(screw%psi - translation%psi)
            call require(error < 0.12_dp*previous, 'zero pitch limit rate')
            previous = error
            call sample_gs(screw, 0.21_dp, -0.17_dp, 1.0_dp, &
                psi, b(1), b(2), b(3), status)
            call require(status == 0, 'screw limit field')
            field_error = norm2(b - b0)
            call require(field_error < 0.12_dp*previous_b, 'zero pitch field rate')
            previous_b = field_error
        end do
        ref = norm2(translation%psi)
        call require(error < 2.0e-4_dp*ref, 'zero pitch within discretization error')
        call require(field_error < 2.0e-4_dp*norm2(b0), 'zero pitch field error')
        pitch = 0
        call solve_case(screw, GS_SCREW, 12, 2)
        call require(norm2(screw%psi - translation%psi) < 2.0e-10_dp*max(1.0_dp, ref), &
            'exact zero pitch translation')
        call sample_gs(screw, 0.21_dp, -0.17_dp, 1.0_dp, &
            psi, b(1), b(2), b(3), status)
        call require(status == 0, 'zero pitch field')
        call require(norm2(b - b0) < 2.0e-10_dp, 'exact zero pitch field equality')
    end subroutine

    subroutine derivatives()
        type(gs_state_t) :: state
        type(gs_profile_law_t) :: law
        type(gs_profile_jacobian_t) :: jac, unused
        real(dp), allocatable :: rhs(:),rp(:),rm(:),v(:),u(:),jv(:),jtu(:),av(:),solution(:)
        real(dp) :: eps, error, psi, b(3)
        integer :: i, status
        call leaf_derivatives()
        call mesh(state, GS_SCREW, 8, 2)
        call profile(law)
        allocate (rhs(size(state%psi)), rp(size(state%psi)), rm(size(state%psi)), &
            v(size(state%psi)), u(size(state%psi)), jv(size(state%psi)), &
            jtu(size(state%psi)), av(size(state%psi)), solution(size(state%psi)))
        do i = 1, size(state%psi)
            call exact(state%points(1, i), state%points(2, i), psi, b)
            state%psi(i) = psi
            v(i) = sin(real(i, dp)); u(i) = cos(real(i, dp))
        end do
        where (state%boundary) v = 0
        where (state%boundary) u = 0
        call gs_cubic_load(state, 1.0_dp, law, rhs, jac, status, tolerance=1.0e-11_dp)
        call require(status == 0, 'derivative load')
        call apply_gs_profile_jacobian(state, jac, v, jv)
        eps = 1.0e-6_dp
        state%psi = state%psi + eps*v
        call gs_cubic_load(state, 1.0_dp, law, rp, unused, status, tolerance=1.0e-11_dp)
        call require(status == 0, 'positive difference load')
        state%psi = state%psi - 2*eps*v
        call gs_cubic_load(state, 1.0_dp, law, rm, unused, status, tolerance=1.0e-11_dp)
        call require(status == 0, 'negative difference load')
        state%psi = state%psi + eps*v
        call apply_gs(state, v, av)
        error = norm2(jv - av + (rp - rm)/(2*eps))/norm2(jv)
        call require(error < 1.0e-8_dp, 'equation Jacobian finite difference')
        call apply_gs_profile_jacobian(state, jac, u, jtu, transpose=.true.)
        call require(abs(dot_product(u, jv) - dot_product(v, jtu)) < 1.0e-10_dp, &
            'equation transpose pairing')
        call gs_profile_response(state, jac, u, solution, status)
        if (status /= 0) print *, 'response status: ', status
        call require(status == 0, 'implicit equation response')
        call apply_gs_profile_jacobian(state, jac, solution, jv)
        call require(norm2(jv - u) < 1.0e-8_dp*norm2(u), 'implicit response residual')
    end subroutine

    subroutine leaf_derivatives()
        real(dp) :: args(6), direction(6), up(6), down(6), db(3), bp(3), bm(3)
        real(dp) :: dual(3), gradient(6), eps, ds, sp, sm
        args = [0.2_dp, -0.3_dp, 0.7_dp, 0.4_dp, -0.5_dp, 1.2_dp]
        direction = [0.3_dp, 0.2_dp, -0.1_dp, 0.7_dp, 0.6_dp, -0.2_dp]
        dual = [0.3_dp, -0.2_dp, 0.8_dp]; eps = 1.0e-6_dp
        up = args + eps*direction; down = args - eps*direction
        call gs_screw_field_jvp(args(1), args(2), args(3), args(4), args(5), args(6), &
            direction(1), direction(2), direction(3), direction(4), direction(5), &
            direction(6), db(1), db(2), db(3))
        call gs_screw_field(up(1), up(2), up(3), up(4), up(5), up(6), bp(1), bp(2), bp(3))
        call gs_screw_field(down(1), down(2), down(3), down(4), down(5), down(6), &
            bm(1), bm(2), bm(3))
        call require(norm2(db - (bp - bm)/(2*eps)) < 1.0e-9_dp, 'field leaf difference')
        call gs_screw_field_vjp(args(1), args(2), args(3), args(4), args(5), args(6), &
            dual(1), dual(2), dual(3), gradient(1), gradient(2), gradient(3), &
            gradient(4), gradient(5), gradient(6))
        call require(abs(dot_product(dual, db) - dot_product(gradient, direction)) &
            < 1.0e-13_dp, 'field leaf transpose')
        call gs_screw_source_jvp(args(1), args(2), args(3), 1.0_dp, 0.3_dp, 0.2_dp, &
            args(6), direction(1), direction(2), direction(3), 0.0_dp, 0.0_dp, &
            0.0_dp, direction(6), ds)
        call gs_screw_source(up(1), up(2), up(3), 1.0_dp, 0.3_dp, 0.2_dp, up(6), sp)
        call gs_screw_source(down(1), down(2), down(3), 1.0_dp, 0.3_dp, 0.2_dp, down(6), sm)
        call require(abs(ds - (sp - sm)/(2*eps)) < 1.0e-9_dp, 'source leaf difference')
    end subroutine

    subroutine guards()
        type(gs_state_t) :: state
        type(gs_coupling_t) :: coupling
        type(gs_symmetry_t) :: symmetry
        type(gs_profile_law_t) :: law
        type(gs_profile_jacobian_t) :: jac
        type(gs_inverse_result_t) :: inverse
        real(dp), allocatable :: rhs(:)
        real(dp) :: m, f, b
        integer :: status
        symmetry = gs_symmetry_t(GS_SCREW, GS_CHART_RZ, 0.5_dp)
        call require(.not. symmetry%valid(), 'reject incompatible chart')
        pitch = 0.5_dp
        call mesh(state, GS_SCREW, 4, 2)
        allocate (rhs(size(state%psi)))
        call gs_coupling_build(state, coupling, status)
        call require(status == GS_UNSUPPORTED_SYMMETRY, 'reject axisymmetric exterior')
        call gs_functional_majorant(state, 1.0_dp, 0.0_dp, 0.0_dp, m, f, b, status)
        call require(status == GS_UNSUPPORTED_SYMMETRY, 'reject axisymmetric certificate')
        call gs_profile_load(state, 1.0_dp, [-1.0_dp, 1.0_dp], &
            [0.0_dp, 0.0_dp], [1.0_dp, 1.0_dp], rhs, status)
        call require(status == GS_UNSUPPORTED_SYMMETRY, 'reject missing F primitive')
        call solve_gs(state, 1.0_dp, 0.0_dp, 0.0_dp, status)
        call require(status == GS_UNSUPPORTED_SYMMETRY, 'reject unspecified fixed F')
        call gs_toroidal_flux(state, inverse, f, status)
        call require(status == GS_UNSUPPORTED_SYMMETRY, 'reject toroidal flux')
        call profile(law)
        call gs_normalize_current(state, 1.0_dp, law, 1.0_dp, f, status)
        call require(status == GS_UNSUPPORTED_SYMMETRY, 'reject toroidal current')
        law%normalized = .true.; law%knots = [0.0_dp, 1.0_dp]
        call gs_cubic_load(state, 1.0_dp, law, rhs, jac, status)
        call require(status == GS_UNSUPPORTED_SYMMETRY, 'reject normalized helical profile')
    end subroutine
end program test_gs_symmetry
