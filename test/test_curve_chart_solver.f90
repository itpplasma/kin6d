module curve_solver_cases
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_ideal_map, only: ideal_frozen_t
    use kin6d_landreman_benchmark, only: landreman_frozen_t
    implicit none
    type, extends(ideal_frozen_t) :: relabelled_landreman_t
    contains
        procedure :: evaluate => relabelled
    end type
    type, extends(ideal_frozen_t) :: tube_field_t
    contains
        procedure :: evaluate => tube
    end type
contains
    subroutine relabelled(field, point, chi, chi_gradient, b3, pressure)
        class(relabelled_landreman_t), intent(in) :: field
        real(dp), intent(in) :: point(2)
        real(dp), intent(out) :: chi, chi_gradient(2), b3, pressure
        type(landreman_frozen_t) :: original
        real(dp) :: gradient(2)
        ! Independent orientation-changing pullback q=(Z,R-1,phi).
        call original%evaluate([1 + point(2), point(1)], chi, gradient, b3, pressure)
        chi = -chi; chi_gradient = -gradient([2, 1])
    end subroutine

    subroutine tube(field, point, chi, chi_gradient, b3, pressure)
        class(tube_field_t), intent(in) :: field
        real(dp), intent(in) :: point(2)
        real(dp), intent(out) :: chi, chi_gradient(2), b3, pressure
        ! Closed reference flux, tangent to the circular wall. Candidate only.
        chi = sum(point**2)/2; chi_gradient = point
        b3 = 1; pressure = 0
    end subroutine
end module

program test_curve_chart_solver
    use curve_solver_cases
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use kin6d_equilibrium_chart, only: curve_chart_t, cylindrical_chart_t
    use kin6d_grad_shafranov, only: gs_state_t, initialize_gs
    use kin6d_landreman_leaves, only: landreman_surface, landreman_axisymmetric
    use kin6d_ideal_map, only: ideal_map_t, ideal_map_init, ideal_map_freeze, &
                          ideal_map_energy, ideal_map_solve, ideal_map_solve_report_t, &
                               ideal_map_reduce, ideal_map_expand, ideal_map_hessian, &
                      ideal_map_eliminate_omega, ideal_map_residual, IDEAL_MAP_BAD_INPUT
    use fortfem_element_to_csr, only: csr_matrix_t
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp), delta = 0.01_dp
    type(curve_chart_t) :: chart
    real(dp) :: c(3, 0:2), s(3, 0:2)
    character(16) :: mode
    integer :: st
    call get_command_argument(1, mode)
    c = 0; s = 0
    select case (trim(mode))
    case ('landreman')
        c(1, 1) = 1; s(2, 1) = 1
        call chart%init(c, s, 1, st)
        if (st /= 0) error stop 'circular chart'
        call circle_solve()
    case ('figure8')
        s(1, 1) = 1; s(2, 2) = 1; c(3, 1) = 0.3_dp
        call chart%init(c, s, 1, st)
        if (st /= 0) error stop 'lifted figure8 chart'
        call figure8_solve()
    case default
        error stop 'chart solver test mode'
    end select
contains
    subroutine wall(t, point)
        real(dp), intent(in) :: t
        real(dp), intent(out) :: point(2)
        real(dp) :: dr, dz
        if (mode == 'landreman') then
            call landreman_surface(2*pi*t, delta, point(1), point(2), dr, dz)
        else
            point = 0.02_dp*[cos(2*pi*t), sin(2*pi*t)]
        end if
    end subroutine

    subroutine mesh_build(mesh)
        type(gs_state_t), intent(out) :: mesh
        real(dp) :: contour(2, 16), h
        integer :: i, status
        if (mode == 'figure8') then
            call normal_plane_fan(mesh)
            return
        end if
        do i = 1, 16
            call wall(real(i - 1, dp)/16, contour(:, i))
        end do
        h = 2*pi*sqrt(delta)/16
        call initialize_gs(mesh, contour, h*h/2, status, 2, curve=wall)
        if (status /= 0) error stop 'chart reference mesh'
    end subroutine

    subroutine normal_plane_fan(mesh)
        ! P2 test fixture: a disk in (q1,q2), with curved boundary edges.
        ! ideal_map_init reads only these three mesh arrays; no GS operator
        ! or positive-cylindrical-radius condition applies to them.
        type(gs_state_t), intent(out) :: mesh
        integer, parameter :: n = 12
        integer :: i, j
        allocate (mesh%points(2, 3*n + 1), mesh%triangles(6, n), &
            mesh%boundary(3*n + 1))
        mesh%points(:, 3*n + 1) = 0
        mesh%boundary = .false.
        mesh%boundary(:n) = .true.; mesh%boundary(2*n + 1:3*n) = .true.
        do i = 1, n
            j = mod(i, n) + 1
            call wall(real(i - 1, dp)/n, mesh%points(:, i))
            mesh%points(:, n + i) = mesh%points(:, i)/2
            call wall((i - 0.5_dp)/n, mesh%points(:, 2*n + i))
            mesh%triangles(:, i) = [i, j, 3*n + 1, 2*n + i, n + j, n + i]
        end do
    end subroutine

    subroutine circle_solve()
        type(gs_state_t) :: mesh, relabelled_mesh
        type(cylindrical_chart_t) :: cylinder
        type(landreman_frozen_t) :: original
        type(relabelled_landreman_t) :: pulled
        type(ideal_map_t) :: a, b
        type(ideal_map_solve_report_t) :: ra, rb
        real(dp) :: wa, wb, err, eb, ep, scale
        integer :: i, status
        call mesh_build(mesh)
        relabelled_mesh = mesh
        relabelled_mesh%points(1, :) = mesh%points(2, :)
        relabelled_mesh%points(2, :) = mesh%points(1, :) - 1
        relabelled_mesh%triangles = mesh%triangles([1, 3, 2, 6, 5, 4], :)
        call ideal_map_init(a, mesh, cylinder, 1, 0, status)
        if (status /= 0) error stop 'cylindrical map init'
        call ideal_map_init(b, relabelled_mesh, chart, 1, 0, status, nphi=16)
        if (status /= 0) error stop 'curve map init'
        call ideal_map_freeze(a, original, status)
        if (status /= 0) error stop 'cylindrical freeze'
        call ideal_map_freeze(b, pulled, status)
        if (status /= 0) error stop 'curve freeze'
        if (a%sigma /= -1 .or. b%sigma /= 1) error stop 'chart orientations'
        ! This comparison fails with the old q1^2 omega gauge, since q1=Z.
        err = maxval(abs(a%omega_gauge - &
                         b%omega_gauge([1, 3, 2, 6, 5, 4], [1, 3, 2, 6, 5, 4], :)))
        scale = maxval(abs(a%omega_gauge))
        print '(a,es12.4)', '[V] circular gauge relative difference ', err/scale
        if (err > 1e-11_dp*scale) error stop 'chart gauge covariance'
        do i = 1, a%nodes
            a%x(1, 1, i) = a%x(1, 1, i) + 0.001_dp*a%axis_t(i)
            b%x(2, 1, i) = b%x(2, 1, i) + 0.001_dp*b%axis_t(i)
        end do
        call ideal_map_energy(a, a%x, wa, status)
        if (status /= 0) error stop 'cylindrical energy'
        call ideal_map_energy(b, b%x, wb, status)
        if (status /= 0) error stop 'curve energy'
        if (abs(wa - wb) > 1e-12_dp*abs(wa)) error stop 'chart energy covariance'
        call ideal_map_solve(a, 1e-6_dp, 80, ra, status)
        if (status /= 0) error stop 'cylindrical Landreman solve'
        call ideal_map_solve(b, 1e-6_dp, 80, rb, status)
        if (status /= 0) error stop 'curve Landreman solve'
        err = maxval(abs(a%x(1, :, :) - b%x(2, :, :) - 1))
        err = max(err, maxval(abs(a%x(2, :, :) - b%x(1, :, :))))
        print '(a,3es12.4)', '[V] solved chart position/energy/gradient differences ', &
            err, abs(ra%energy - rb%energy), abs(ra%gradient - rb%gradient)
        if (err > 1e-5_dp) error stop 'Landreman chart solve covariance'
        call landreman_error(b, eb, ep)
        print '(a,2es12.4)', '[V] Landreman max field/pressure errors ', eb, ep
        if (eb > 0.03_dp*2*sqrt(delta)) error stop 'Landreman field accuracy'
        if (ep > 0.03_dp*2*delta) error stop 'Landreman pressure accuracy'
        print '(a)', '[V] fixed-boundary candidate; continuous bounds unavailable'
    end subroutine

    subroutine landreman_error(map, eb, ep)
        type(ideal_map_t), intent(in) :: map
        real(dp), intent(out) :: eb, ep
        real(dp) :: q(3), qd(3, 3), pos(3), d1(3, 3), d2(3, 3, 3), d3(3, 3, 3, 3)
        real(dp) :: dx(3, 3), jac, bc(3), exact(3), current(3), gp(2), psi, p, r
        integer :: e, g
        eb = 0; ep = 0
        do e = 1, map%cells
            do g = 1, map%nq
                q = matmul(map%x(:, 1, map%conn(:, e)), map%shape(:, g))
                qd = 0; qd(3, 3) = 1
                qd(:, :2) = matmul(map%x(:, 1, map%conn(:, e)), &
                                   transpose(map%dshape(:, :, g)))
                call map%chart%jet(q, pos, d1, d2, d3)
                dx = matmul(d1, qd)
                jac = map%sigma*det(dx)
                bc = matmul(dx, map%b(:, g, e))/jac
                r = norm2(pos(:2))
                call landreman_axisymmetric(r, pos(3), delta, psi, p, &
                     exact(1), exact(2), exact(3), current(1), current(2), current(3), &
                                            gp(1), gp(2))
                eb = max(eb, norm2(bc - [exact(1), exact(2), exact(3)]))
                ep = max(ep, abs(map%p(g, e) - p))
            end do
        end do
    end subroutine

    real(dp) function det(a) result(v)
        real(dp), intent(in) :: a(3, 3)
        v = a(1, 1)*(a(2, 2)*a(3, 3) - a(2, 3)*a(3, 2)) &
            - a(1, 2)*(a(2, 1)*a(3, 3) - a(2, 3)*a(3, 1)) &
            + a(1, 3)*(a(2, 1)*a(3, 2) - a(2, 2)*a(3, 1))
    end function

    subroutine figure8_solve()
        type(gs_state_t) :: mesh
        type(tube_field_t) :: field
        type(ideal_map_t) :: map, fine, zero_mode
        type(ideal_map_solve_report_t) :: report
        real(dp) :: energy, refined, residual, sigma, curvature, radius
        real(dp), allocatable :: grad(:)
        integer :: status
        call chart%tube_radii(2048, curvature, radius, status)
        if (status /= 0) error stop 'figure8 tube admission'
        if (radius <= 0.02_dp) error stop 'figure8 tube radius'
        call mesh_build(mesh)
        call ideal_map_init(map, mesh, chart, 1, 1, status)
        if (status /= IDEAL_MAP_BAD_INPUT) error stop 'missing chart quadrature'
        call ideal_map_init(map, mesh, chart, 2, 1, status, nphi=64)
        if (status /= IDEAL_MAP_BAD_INPUT) error stop 'chart field periods'
        call ideal_map_init(map, mesh, chart, 1, 1, status, nphi=64)
        if (status /= 0) error stop 'figure8 map init'
        call ideal_map_freeze(map, field, status)
        if (status /= 0) error stop 'figure8 freeze'
        call derivatives(map)
        call ideal_map_solve(map, 1e-4_dp, 100, report, status)
        print '(a,i0,a,3es12.4)', '[V] figure8 solve status=', status, &
            ' g0/relative gradient/min J ', report%initial_gradient, &
            report%gradient, report%min_jacobian
        if (status /= 0) error stop 'figure8 solve'
        if (.not. all(ieee_is_finite([report%gradient, report%min_jacobian]))) &
            error stop 'figure8 finite solve'
        if (report%min_jacobian <= 0) error stop 'figure8 positive Jacobian'
        if (report%gradient > 0.01_dp) &
            error stop 'figure8 stationarity'
        call ideal_map_init(fine, mesh, chart, 1, 1, status, nphi=128)
        if (status /= 0) error stop 'refined periodic quadrature'
        call ideal_map_freeze(fine, field, status)
        if (status /= 0) error stop 'refined freeze'
        fine%x = map%x
        call ideal_map_energy(map, map%x, energy, status)
        if (status /= 0) error stop 'figure8 energy'
        call ideal_map_energy(fine, fine%x, refined, status)
        if (status /= 0) error stop 'figure8 refined energy'
        print '(a,es12.4)', '[V] figure8 periodic energy refinement ', &
            abs(energy - refined)/abs(refined)
        if (abs(energy - refined) > 0.01_dp*abs(refined)) &
            error stop 'figure8 periodic quadrature accuracy'
        allocate (grad(map%nred))
        call ideal_map_init(zero_mode, mesh, chart, 1, 0, status, nphi=64)
        if (status /= 0) error stop 'zero-mode curve chart'
        call ideal_map_freeze(zero_mode, field, status)
        if (status /= 0) error stop 'zero-mode freeze'
        call ideal_map_eliminate_omega(zero_mode, zero_mode%x, energy, grad, status)
        if (status /= IDEAL_MAP_BAD_INPUT) error stop 'cylindrical elimination only'
        call ideal_map_residual(map, 1.0_dp, residual, sigma, status)
        if (status /= IDEAL_MAP_BAD_INPUT) error stop 'cylindrical residual only'
        print '(a)', '[V] schematic fixed-wall candidate; continuous bounds unavailable'
    end subroutine

    subroutine derivatives(map)
        type(ideal_map_t), intent(in) :: map
        type(csr_matrix_t) :: hessian
        real(dp), allocatable :: d(:), full(:), gp(:), gm(:), gr(:), hr(:), hd(:)
        real(dp), allocatable :: trial(:, :, :)
        real(dp) :: h, wp, wm, energy, err
        integer :: i, row, k, status
        allocate (d(map%nred), full(map%ndof), gp(map%ndof), gm(map%ndof), &
                  gr(map%nred), hr(map%nred), hd(map%nred))
        do i = 1, map%nred
            d(i) = 0.001_dp*sin(real(i, dp))
        end do
        where (map%fixed) d = 0
        call ideal_map_expand(map, d, full)
        h = 1e-4_dp
        trial = map%x + h*reshape(full, shape(map%x))
        call ideal_map_energy(map, trial, wp, status, gp)
        if (status /= 0) error stop 'figure8 derivative plus'
        trial = map%x - h*reshape(full, shape(map%x))
        call ideal_map_energy(map, trial, wm, status, gm)
        if (status /= 0) error stop 'figure8 derivative minus'
        call ideal_map_reduce(map, (gp - gm)/(2*h), hr)
        call ideal_map_energy(map, map%x, energy, status, gp)
        if (status /= 0) error stop 'figure8 derivative center'
        call ideal_map_reduce(map, gp, gr)
        err = abs((wp - wm)/(2*h) - dot_product(gr, d))
        if (err > 1e-6_dp*max(1e-10_dp, abs(dot_product(gr, d)))) &
            error stop 'figure8 energy gradient'
        call ideal_map_hessian(map, map%x, hessian)
        hd = 0
        do row = 1, map%nred
            do k = hessian%rowptr(row), hessian%rowptr(row + 1) - 1
                hd(row) = hd(row) + hessian%val(k)*d(hessian%colidx(k))
            end do
        end do
        where (map%fixed) hr = 0
        where (map%fixed) hd = 0
        err = norm2(hr - hd)/norm2(hd)
        print '(a,es12.4)', '[V] figure8 Hessian directional FD error ', err
        if (.not. ieee_is_finite(err) .or. err > 1e-5_dp) &
            error stop 'figure8 gradient Hessian'
    end subroutine
end program
