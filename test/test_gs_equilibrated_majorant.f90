program test_gs_equilibrated_majorant
    !! Equilibrated-flux GS majorant against the exact weighted energy error
    !! of a manufactured solution on a straight R x Z box (R0 = 3, a = 1,
    !! elongation 1.5): psi* = sin(alpha(R - R0 + a)) cos(beta Z)(1 + 0.3(R - R0)),
    !! f = -div(grad psi*/R). Checks the upper bound, convergence of the bound
    !! at the energy-error rate h^p and an efficiency index close to one, for
    !! P1, P2 and P3, and compares with the global functional majorant.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_grad_shafranov, only: gs_state_t, initialize_gs, gs_linear_response, &
        gs_cell_data
    use kin6d_gs_majorant, only: gs_functional_majorant
    use kin6d_gs_equilibrated, only: gs_equilibrated_majorant, &
        equilibrated_flux_estimate_t, GS_EQ_OK
    use fortfem_triangle_duffy_quadrature, only: triangle_duffy_quadrature
    implicit none
    real(dp), parameter :: r0 = 3.0_dp, a = 1.0_dp, kap = 1.5_dp, slope = 0.3_dp
    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp), parameter :: alpha = pi/(2*a), beta = pi/(2*kap*a)
    type(gs_state_t) :: state
    type(equilibrated_flux_estimate_t) :: est
    real(dp), allocatable :: contour(:, :), rhs(:), psi(:)
    real(dp) :: error, previous_error, previous_bound, effectivity, rate
    real(dp) :: functional, flux_defect, balance_defect, area
    integer :: order, level, status, levels

    print '(a)', 'p level dofs energy_error equilibrated eff functional eff'
    do order = 1, 3
        levels = merge(2, 3, order == 3)
        previous_error = 0
        previous_bound = 0
        do level = 1, levels
            area = 0.1_dp/4**(level - 1)
            call rectangle(area, contour)
            call initialize_gs(state, contour, area, status, order)
            if (status /= 0) error stop 'mesh'
            call load(rhs)
            call gs_linear_response(state, rhs, psi, status, &
                relative_tolerance=1.0e-12_dp)
            if (status > 1) error stop 'solve'
            error = energy_error(psi)
            call gs_equilibrated_majorant(state, psi, source, est, status)
            if (status /= GS_EQ_OK) error stop 'equilibrated majorant status'
            effectivity = est%bound/error
            functional = 0
            if (order > 1) then
                state%psi = psi
                call gs_functional_majorant(state, 1.0_dp, 0.0_dp, 0.0_dp, &
                    functional, flux_defect, balance_defect, status, &
                    forcing=forcing)
                if (status /= 0) error stop 'functional majorant'
            end if
            print '(i1, i3, i8, 5es11.3)', order, level, count(.not. state%boundary), &
                error, est%bound, effectivity, functional, functional/error
            if (.not. (est%bound >= error)) error stop 'equilibrated bound below error'
            if (effectivity > 1.15_dp) error stop 'equilibrated efficiency index'
            if (level > 1) then
                rate = log(previous_bound/est%bound)/log(2.0_dp)
                if (rate < order - 0.3_dp) error stop 'equilibrated bound rate'
            end if
            previous_error = error
            previous_bound = est%bound
        end do
    end do

contains

    subroutine exact(r, z, value, grad)
        real(dp), intent(in) :: r, z
        real(dp), intent(out) :: value, grad(2)
        real(dp) :: s, c, l

        s = sin(alpha*(r - r0 + a))
        c = cos(beta*z)
        l = 1 + slope*(r - r0)
        value = s*c*l
        grad = [(alpha*cos(alpha*(r - r0 + a))*l + slope*s)*c, -beta*s*l*sin(beta*z)]
    end subroutine exact

    real(dp) function forcing(r, z) result(f)
        !! f = -d/dR(psi_R/R) - psi_ZZ/R with psi = cos(beta Z) g(R).
        real(dp), intent(in) :: r, z
        real(dp) :: s, cs, l, g, g1, g2

        s = sin(alpha*(r - r0 + a))
        cs = cos(alpha*(r - r0 + a))
        l = 1 + slope*(r - r0)
        g = s*l
        g1 = alpha*cs*l + slope*s
        g2 = -alpha**2*s*l + 2*slope*alpha*cs
        f = cos(beta*z)*(-(g2 - beta**2*g)/r + g1/r**2)
    end function forcing

    real(dp) function source(cell, r, z, p) result(f)
        integer, intent(in) :: cell
        real(dp), intent(in) :: r, z, p
        f = forcing(r, z) + 0.0_dp*(cell + p)
    end function source

    subroutine rectangle(max_area, points)
        real(dp), intent(in) :: max_area
        real(dp), allocatable, intent(out) :: points(:, :)
        real(dp) :: h, corner(2, 5)
        integer :: nr, nz, side, j, m, n

        h = sqrt(4*max_area/sqrt(3.0_dp))
        nr = ceiling(2*a/h)
        nz = ceiling(2*kap*a/h)
        corner = reshape([r0 - a, -kap*a, r0 + a, -kap*a, r0 + a, kap*a, &
            r0 - a, kap*a, r0 - a, -kap*a], [2, 5])
        allocate (points(2, 2*(nr + nz)))
        n = 0
        do side = 1, 4
            m = merge(nr, nz, mod(side, 2) == 1)
            do j = 0, m - 1
                n = n + 1
                points(:, n) = corner(:, side) &
                    + (corner(:, side + 1) - corner(:, side))*real(j, dp)/m
            end do
        end do
    end subroutine rectangle

    subroutine load(b)
        real(dp), allocatable, intent(out) :: b(:)
        real(dp), allocatable :: x(:), y(:), w(:)
        real(dp) :: point(2), det, shape(10), g(2, 10), hs(2, 2, 10)
        integer :: e, q, nd, st

        nd = size(state%triangles, 1)
        allocate (b(size(state%psi)))
        b = 0
        call triangle_duffy_quadrature(14, x, y, w, st)
        do e = 1, size(state%triangles, 2)
            do q = 1, size(w)
                call gs_cell_data(state, e, [x(q), y(q)], point, det, shape(:nd), &
                    g(:, :nd), hs(:, :, :nd), st)
                b(state%triangles(:, e)) = b(state%triangles(:, e)) &
                    + w(q)*det*forcing(point(1), point(2))*shape(:nd)
            end do
        end do
        where (state%boundary) b = 0
    end subroutine load

    real(dp) function energy_error(u) result(val)
        real(dp), intent(in) :: u(:)
        real(dp), allocatable :: x(:), y(:), w(:)
        real(dp) :: point(2), det, shape(10), g(2, 10), hs(2, 2, 10), v, gr(2)
        integer :: e, q, nd, st

        nd = size(state%triangles, 1)
        call triangle_duffy_quadrature(16, x, y, w, st)
        val = 0
        do e = 1, size(state%triangles, 2)
            do q = 1, size(w)
                call gs_cell_data(state, e, [x(q), y(q)], point, det, shape(:nd), &
                    g(:, :nd), hs(:, :, :nd), st)
                call exact(point(1), point(2), v, gr)
                val = val + w(q)*det*sum((gr - matmul(g(:, :nd), &
                    u(state%triangles(:, e))))**2)/point(1)
            end do
        end do
        val = sqrt(val)
    end function energy_error

end program test_gs_equilibrated_majorant
