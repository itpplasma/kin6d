program test_gs_checked_source_remainder
    !! Independent manufactured quartic, rational source numerator approximated
    !! by a Taylor polynomial of degree p+3 with an explicit geometric remainder.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_grad_shafranov, only: gs_state_t, initialize_gs, gs_linear_response, &
        gs_cell_data
    use kin6d_gs_equilibrated, only: gs_equilibrated_majorant, &
        equilibrated_flux_estimate_t
    use kin6d_gs_quadrature_budget, only: gs_checked_fixed_source_majorant
    use fortfem_triangle_duffy_quadrature, only: triangle_duffy_quadrature
    use fortfem_verified_bernstein_products, only: verified_bernstein_affine_product, &
        verified_bernstein_elevate
    use fortnum_interval, only: interval_t, interval, operator(+), operator(-), &
        operator(*), operator(/), operator(**)
    use omp_lib, only: omp_get_wtime, omp_set_num_threads
    implicit none
    real(dp), parameter :: r0 = 3.0_dp, a = 1.0_dp, kap = 1.0_dp
    type(gs_state_t) :: state
    type(equilibrated_flux_estimate_t) :: est
    type(interval_t) :: bound
    real(dp), allocatable :: contour(:, :), rhs(:), psi(:)
    real(dp) :: error, previous, area, before, solve_time, audit_time, total_time
    real(dp) :: start, reference_bound, times(3)
    integer :: order, level, status, numerator_degree, first_order, first_level
    integer :: last_order, last_level, threads, repeat
    character(32) :: argument
    logical :: bad_tail

    first_order = 1; last_order = 3; first_level = 1; last_level = 3
    call get_command_argument(1, argument)
    if (trim(argument) == 'study') then
        first_order = 3; last_order = 3; first_level = 4; last_level = 5
        call omp_set_num_threads(4)
    end if
    bad_tail = .false.
    print '(a)', 'p level dofs error sampled checked efficiency ratio solve audit total'
    do order = first_order, last_order
        numerator_degree = order + 3
        previous = 0
        do level = first_level, last_level
            area = 0.1_dp/4**(level - 1)
            call rectangle(area, contour)
            start = omp_get_wtime()
            call initialize_gs(state, contour, area, status, order)
            if (status /= 0) error stop 'mesh'
            call load(rhs)
            before = omp_get_wtime()
            call gs_linear_response(state, rhs, psi, status, &
                relative_tolerance=1.0e-12_dp)
            solve_time = omp_get_wtime() - before
            total_time = omp_get_wtime() - start
            if (status > 1) error stop 'solve'
            error = energy_error(psi)
            call gs_equilibrated_majorant(state, psi, source, est, status)
            if (status /= 0) error stop 'sampled reconstruction'
            before = omp_get_wtime()
            call gs_checked_fixed_source_majorant(state, psi, est, numerator, &
                numerator_degree, 1, bound, status, tail)
            audit_time = omp_get_wtime() - before
            if (status /= 0) error stop 'source remainder audit'
            if (bound%hi < error) error stop 'checked bound below true error'
            if (bound%hi > 1.2_dp*error) error stop 'checked efficiency'
            if (bound%hi > 1.05_dp*est%bound) error stop 'checked sampled ratio'
            if (level > first_level) then
                if (log(previous/bound%hi)/log(2.0_dp) < order - 0.4_dp) &
                    error stop 'checked rate'
            end if
            previous = bound%hi
            print '(i1,i3,i7,8es12.4)', order, level, count(.not. state%boundary), &
                error, est%bound, bound%hi, bound%hi/error, bound%hi/est%bound, &
                solve_time, audit_time, total_time
            if (trim(argument) == 'study') then
                reference_bound = bound%hi
                do threads = 1, 4
                    if (threads == 3) cycle
                    call omp_set_num_threads(threads)
                    do repeat = 1, 3
                        before = omp_get_wtime()
                        call gs_checked_fixed_source_majorant(state, psi, est, &
                            numerator, numerator_degree, 1, bound, status, tail)
                        times(repeat) = omp_get_wtime() - before
                        if (status /= 0) error stop 'scaling audit'
                        if (abs(bound%hi - reference_bound) > &
                            1.e-6_dp*reference_bound) error stop 'thread consistency'
                    end do
                    print '(a,i2,i2,3es13.5)', 'AUDIT ', level, threads, times
                end do
                call omp_set_num_threads(4)
            end if
        end do
    end do
    bad_tail = .true.
    call gs_checked_fixed_source_majorant(state, psi, est, numerator, &
        numerator_degree, 1, bound, status, tail)
    if (status == 0) error stop 'missing source enclosure accepted'
    print '(a)', 'PASS: non-polynomial fixed-source remainder audit'
contains

    subroutine exact(r, z, value, grad)
        real(dp), intent(in) :: r, z
        real(dp), intent(out) :: value, grad(2)
        real(dp) :: x
        x = r - r0
        value = (1 - x*x)*(1 - z*z)
        grad = [-2*x*(1 - z*z), -2*z*(1 - x*x)]
    end subroutine exact

    real(dp) function forcing(r, z) result(f)
        real(dp), intent(in) :: r, z
        f = 6*(1 - z*z)/r**2 + 2*(1 - (r - r0)**2)/r
    end function forcing

    real(dp) function source(cell, r, z, p) result(f)
        integer, intent(in) :: cell
        real(dp), intent(in) :: r, z, p
        f = forcing(r, z) + 0.0_dp*(cell + p)
    end function source

    subroutine numerator(cell, coefficients, ok)
        !! R*f = 6(1-z^2)/R + 2(1-x^2). Expand 1/R at the exact stored
        !! center c, through power numerator_degree-2.
        integer, intent(in) :: cell
        type(interval_t), intent(out) :: coefficients(:)
        logical, intent(out) :: ok
        type(interval_t) :: r(3), z(3), x(3), v(3), c, p(28), q(28), tmp(28)
        real(dp) :: center
        integer :: j, n, power
        r = interval(state%points(1, state%triangles(1:3, cell)))
        z = interval(state%points(2, state%triangles(1:3, cell)))
        x = r - interval(r0)
        center = sum(state%points(1, state%triangles(1:3, cell)))/3
        c = interval(center)
        v = (c - r)/c
        coefficients = interval(0.0_dp)
        do power = 0, numerator_degree - 2
            p = interval(1.0_dp)
            call verified_bernstein_affine_product(0, p, interval(1.0_dp) - z, q)
            call verified_bernstein_affine_product(1, q, interval(1.0_dp) + z, p)
            n = 2
            do j = 1, power
                call verified_bernstein_affine_product(n, p, v, tmp)
                p = tmp; n = n + 1
            end do
            do while (n < numerator_degree)
                call verified_bernstein_elevate(n, p, tmp)
                p = tmp; n = n + 1
            end do
            coefficients = coefficients + interval(6.0_dp)*p(:size(coefficients))/c
        end do
        p = interval(1.0_dp)
        call verified_bernstein_affine_product(0, p, interval(1.0_dp) - x, q)
        call verified_bernstein_affine_product(1, q, interval(1.0_dp) + x, p)
        do n = 2, numerator_degree - 1
            call verified_bernstein_elevate(n, p, tmp)
            p = tmp
        end do
        coefficients = coefficients + interval(2.0_dp)*p(:size(coefficients))
        ok = .true.
    end subroutine numerator

    subroutine tail(cell, radius, ok)
        !! Geometric identity, no transcendental estimate:
        !! |1/R - sum_{j=0}^k v^j/c| = |v|^(k+1)/R, v=(c-R)/c.
        !! On this box |1-z^2| <= 1. The affine R hull is exact.
        integer, intent(in) :: cell
        type(interval_t), intent(out) :: radius
        logical, intent(out) :: ok
        type(interval_t) :: r(3), c, d, v
        real(dp) :: center
        r = interval(state%points(1, state%triangles(1:3, cell)))
        center = sum(state%points(1, state%triangles(1:3, cell)))/3
        c = interval(center)
        d = interval(minval(r%lo), maxval(r%hi)) - c
        v = interval(max(abs(d%lo), abs(d%hi)))/c
        radius = interval(6.0_dp)*v**(numerator_degree - 1)/interval(minval(r%lo))
        ok = .not. bad_tail
    end subroutine tail

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

end program test_gs_checked_source_remainder
