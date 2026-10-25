module fitted_profile_fixture
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_equilibrium_state, only: equilibrium_field_t, EQ_VALUES, &
                                       EQ_BR, EQ_BZ, EQ_BPHI, EQ_P
    use fortfem_meridian_fourier_scalar, only: meridian_fourier_source_t
    use kin6d_nonnested_profile, only: nn_profile_t
    use nn_solovev_fixture, only: r0, pp, ftor, exact, laplacian
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp)
    type, extends(equilibrium_field_t) :: analytic_field_t
        logical :: solovev = .false.
    contains
        procedure :: sample => analytic_sample
    end type
    type, extends(nn_profile_t) :: analytic_profile_t
        logical :: solovev = .false.
        real(dp) :: amplitude = 1
    contains
        procedure :: evaluate => analytic_profile
    end type
    type, extends(meridian_fourier_source_t) :: strong_source_t
    contains
        procedure :: source => strong_source
    end type
contains
    subroutine strong_source(self, t, x, f0, f, status)
        class(strong_source_t), intent(in) :: self
        integer, intent(in) :: t
        real(dp), intent(in) :: x(:, :)
        real(dp), intent(out) :: f0(:), f(:, :)
        integer, intent(out) :: status
        real(dp) :: base(size(x, 2))
        base = exp(x(1, :) - 3)*cos(2*x(2, :))
        f0 = (3 - 1/x(1, :))*base*(1 + 0.2_dp*cos(x(3, :))) &
             + 0.2_dp*base*cos(x(3, :))/x(1, :)**2
        f = 0; status = 0
        if (t < 1) status = 1
    end subroutine
    subroutine circle(u, theta, x)
        real(dp), intent(in) :: u, theta
        real(dp), intent(out) :: x(2)
        x = [3 + u*cos(theta), u*sin(theta)]
    end subroutine
    function old_source(x) result(s)
        real(dp), intent(in) :: x(2)
        real(dp) :: s
        s = -pp*laplacian(x)
    end function
    function old_source_3d(x, phi) result(s)
        real(dp), intent(in) :: x(2), phi
        real(dp) :: s
        s = old_source(x)*(1 + 0.1_dp*cos(phi))
    end function
    subroutine analytic_sample(field, points, values, statuses)
        class(analytic_field_t), intent(in) :: field
        real(dp), intent(in) :: points(:, :)
        real(dp), intent(out) :: values(:, :)
        integer, intent(out) :: statuses(:)
        real(dp) :: psi, b(2)
        integer :: i
        values = 0; statuses = 0
        do i = 1, size(points, 2)
            values(EQ_BPHI, i) = 3/points(1, i)
            if (.not. field%solovev) cycle
            call exact(points(:2, i), psi, b)
            values(EQ_BR, i) = b(1); values(EQ_BZ, i) = b(2)
            values(EQ_BPHI, i) = ftor/points(1, i)
        end do
    end subroutine
    subroutine analytic_profile(self, x, p, gradp, lambda, gradphi, status)
        class(analytic_profile_t), intent(in) :: self
        real(dp), intent(in) :: x(:, :)
        real(dp), intent(out) :: p(:), gradp(:, :), lambda(:), gradphi(:, :)
        integer, intent(out) :: status
        real(dp) :: psi, b(2), base, mode
        integer :: i
        gradphi = 0; gradp = 0
        do i = 1, size(x, 2)
            if (self%solovev) then
                call exact(x(:2, i), psi, b)
                p(i) = pp*psi
                gradp(:, i) = pp*[x(1, i)*b(2), 0.0_dp, -x(1, i)*b(1)]
                lambda(i) = pp*ftor/(sum(b*b) + (ftor/x(1, i))**2)
            else
                base = exp(x(1, i) - 3)*cos(2*x(2, i))
                mode = 1 + 0.2_dp*cos(x(3, i))
                p(i) = base*mode
                gradp(1, i) = p(i)
                gradp(2, i) = -0.2_dp*base*sin(x(3, i))/x(1, i)
                gradp(3, i) = -2*exp(x(1, i) - 3)*sin(2*x(2, i))*mode
                lambda(i) = 0.1_dp*sin(x(1, i))
            end if
        end do
        p = self%amplitude*p; gradp = self%amplitude*gradp
        status = 0
    end subroutine
end module

program test_nonnested_profile
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fitted_profile_fixture
    use nn_solovev_fixture, only: minor, mapping, coordinates
    use kin6d_mhd_mesh, only: flux_mesh_t, flux_mesh_build
    use fortfem_meridian_feec_n, only: meridian_lagrange_nodes, meridian_qpoint
use fortfem_meridian_fourier_scalar, only: meridian_fourier_value, meridian_fourier_load
    use fortsparse, only: csc_t, csc_matvec
    use kin6d_equilibrium_nonnested, only: free_b_field_t, nn_build
    use kin6d_nonnested_profile, only: nn_fitted_source_t, nn_fit_source, &
                                       nn_state_profile_t, nn_profile_defect_t
    use kin6d_nonnested_transport, only: nt_conductivity_t, nt_report_t, nt_solve
    use kin6d_nonnested_coupled, only: nc_problem_t, nc_layout_t, nc_report_t, &
                       nc_layout, nc_solve, nc_fixed, nc_store, nc_jacobian, nc_residual
    implicit none
    character(24) :: mode, arg
    integer :: k, lev
    real(dp) :: errors(4), fields(4), rate
    call get_command_argument(1, mode)
    call get_command_argument(2, arg)
    k = 3
    if (len_trim(arg) > 0) read (arg, *) k
    select case (mode)
    case ('manufactured')
        do lev = 1, 4
            call manufactured(4*2**(lev - 1), k, errors(lev))
        end do
        rate = log(errors(3)/errors(4))/log(2.0_dp)
        print '(a,i0,a,f8.3)', 'MMS P', k, ' pressure order ', rate
        if (rate < k + 0.6_dp) error stop 'manufactured pressure order'
    case ('nested')
        do lev = 1, 4
            call nested(4*2**(lev - 1), k, errors(lev), fields(lev))
        end do
        rate = log(errors(3)/errors(4))/log(2.0_dp)
        print '(a,i0,a,f8.3)', 'nested P', k, ' pressure order ', rate
        if (rate < k + 0.6_dp) error stop 'nested pressure order'
        if (errors(4) > 1.0e-3_dp) error stop 'nested prescribed pressure'
        if (log(fields(3)/fields(4))/log(2.0_dp) < k - 0.4_dp) &
            error stop 'nested field order'
    case ('reproduce')
        call reproduce(k, .false.)
    case ('reproduce3d')
        call reproduce(k, .true.)
    case ('sampling')
        call sampling_check(k)
    case default
        error stop 'mode'
    end select
    print *, 'PASS: frozen fitted-source closure B'
contains
    subroutine manufactured(nr, order, error)
        integer, intent(in) :: nr, order
        real(dp), intent(out) :: error
        type(flux_mesh_t) :: mesh
        type(free_b_field_t) :: space
        type(analytic_field_t), target :: seed
        type(analytic_profile_t) :: target
        type(nn_fitted_source_t) :: fit
        type(nt_conductivity_t) :: tensor
        type(nt_report_t) :: report
        type(strong_source_t) :: independent
        real(dp), allocatable :: nodes(:, :), bc(:), p(:), load(:), strong(:)
        real(dp) :: rings(nr), v(3, 64), xy(2), w, x(3, 1), pe(1), gp(3, 1)
        real(dp) :: lam(1), gf(3, 1), ph, val, err2, norm2p
        real(dp) :: flux(3, 1), drive(1), before(3, 1)
        integer :: i, t, q, j, st
        rings = [(real(i, dp)/nr, i=1, nr)]
        call flux_mesh_build(rings, 8, 1024, order, circle, mesh, st)
        if (st /= 0) error stop 'MMS mesh'
      call nn_build(mesh, order, 1, [3.0_dp, 0.0_dp], 1.0_dp, 3.0_dp, 1.0_dp, space, st)
        if (st /= 0) error stop 'MMS space'
       call nn_fit_source(seed, target, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 0.0_dp, fit, st)
        if (st /= 0) error stop 'MMS fit'
        x(:, 1) = [3.1_dp, 0.2_dp, 0.7_dp]
        call fit%data(x, before, drive, st)
        if (st /= 0) error stop 'fit data'
        target%amplitude = 7
        call fit%data(x, flux, drive, st)
        if (st /= 0 .or. any(flux /= before)) error stop 'target not frozen'
        target%amplitude = 1
        call meridian_lagrange_nodes(space%fs, nodes)
        allocate (bc(3*space%fs%n0)); bc = 0
        do i = 1, space%fs%n0
            bc(i) = exp(nodes(1, i) - 3)*cos(2*nodes(2, i))
            bc(space%fs%n0 + i) = 0.2_dp*bc(i)
        end do
        tensor%field => seed
        call nt_solve(space%fs, 1, 8, tensor, fit, bc, p, report, st)
        if (st /= 0) error stop 'MMS solve'
        if (report%residual > 1.0e-13_dp) error stop 'MMS residual'
        ! Independent strong cylindrical Laplacian, integrated with the same
        ! quadrature, checks the sign and toroidal metric of the weak source.
        call meridian_fourier_load(space%fs, 1, 8, fit, load, st)
        if (st /= 0) error stop 'MMS fitted load'
        call meridian_fourier_load(space%fs, 1, 8, independent, strong, st)
        if (st /= 0) error stop 'MMS independent strong load'
        do i = 0, 2
            where (space%fs%sp(0)%bnd) &
                load(i*space%fs%n0 + 1:(i + 1)*space%fs%n0) = &
                strong(i*space%fs%n0 + 1:(i + 1)*space%fs%n0)
        end do
        val = norm2(load - strong)/norm2(strong)
        print '(a,es14.5)', 'independent strong/weak source defect ', val
        if (val > 1.0e-7_dp) error stop 'strong/weak fitted source'
        err2 = 0; norm2p = 0
        do t = 1, space%fs%mm%topology%n_triangles
            do q = 1, space%fs%nq
                call meridian_qpoint(space%fs, 0, t, q, xy, w, v(:, :space%fs%nl0), st)
                if (st /= 0) error stop 'MMS quadrature'
                do j = 0, 15
                    ph = 2*pi*j/16
                    x(:, 1) = [xy, ph]
                    call target%evaluate(x, pe, gp, lam, gf, st)
                    val = meridian_fourier_value(space%fs%n0, 1, p, &
                                      space%fs%sp(0)%dof(:, t), v(1, :space%fs%nl0), ph)
                    err2 = err2 + w*(val - pe(1))**2/16
                    norm2p = norm2p + w*pe(1)**2/16
                end do
            end do
        end do
        error = sqrt(err2/norm2p)
        print '(a,3i6,es14.5)', 'MMS force-incompatible P/nr/dof/error ', &
            order, nr, size(p), error
    end subroutine

    subroutine setup(nr, order, field, lay, nodes, x, problem, modes)
        integer, intent(in) :: nr, order
        integer, intent(in), optional :: modes
        type(free_b_field_t), intent(out) :: field
        type(nc_layout_t), intent(out) :: lay
        real(dp), allocatable, intent(out) :: nodes(:, :), x(:)
        type(nc_problem_t), intent(out) :: problem
        type(flux_mesh_t) :: mesh
        real(dp) :: rings(nr), psi, b(2), w, xy(2), v(3, 64), area
        integer :: i, t, q, st, nmodes
        nmodes = 0
        if (present(modes)) nmodes = modes
        rings = [(real(i, dp)/nr, i=1, nr)]
        ! Straight cells isolate FE profile order from curved-boundary
        ! geometry interpolation. Boundary values are the analytic trace on
        ! this declared polygon, including nonzero pressure at its corners.
        call flux_mesh_build(rings, 8, 1024, 1, mapping, mesh, st)
        if (st /= 0) error stop 'nested mesh'
        call nn_build(mesh, order, nmodes, [r0, 0.0_dp], minor, ftor, 1.0_dp, field, st)
        if (st /= 0) error stop 'nested space'
        field%coordinates => coordinates
        lay = nc_layout(field)
        call meridian_lagrange_nodes(field%fs, nodes)
        allocate (x(lay%n)); x = 0
        do i = 1, lay%n0
            call exact(nodes(:, i), psi, b)
            x(field%fs%nn + i) = psi
            x(lay%ip + i) = pp*psi
        end do
        x(lay%iv) = 2*pi*r0*r0*pp
        area = 0
        do t = 1, field%fs%mm%topology%n_triangles
            do q = 1, field%fs%nq
                call meridian_qpoint(field%fs, 0, t, q, xy, w, v(:, :field%fs%nl0), st)
                if (st /= 0) error stop 'nested quadrature'
                area = area + w
            end do
        end do
        problem%chi_par = 1; problem%chi_perp = 1
        problem%g_par = 1; problem%g_perp = 1
        problem%current = pp*area
    end subroutine

    subroutine nested(nr, order, error, berror)
        integer, intent(in) :: nr, order
        real(dp), intent(out) :: error, berror
        type(free_b_field_t) :: field
        type(nc_layout_t) :: lay
        type(nc_problem_t) :: problem
        type(nc_report_t) :: report
        type(analytic_field_t) :: seed
        type(analytic_profile_t) :: target
        type(nn_profile_defect_t) :: defects
        type(nn_state_profile_t) :: solved_target
        real(dp), allocatable :: nodes(:, :), x(:)
        real(dp) :: points(3, 64), values(EQ_VALUES, 64), gp(3, 64), gf(3, 64)
        real(dp) :: pe(64), lam(64), psi, b(2), bref(3, 64), weights(64)
        real(dp) :: xy(2), wq, basis(3, 64), ph, ep2, np2
        integer :: i, st, sts(64), t, q
        call setup(nr, order, field, lay, nodes, x, problem)
        seed%solovev = .true.; target%solovev = .true.
        problem%closure%pressure = 'B_frozen_fit'
        allocate (problem%fitted)
        call nn_fit_source(seed, target, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, &
                           x(lay%iv), problem%fitted, st)
        if (st /= 0) error stop 'nested fit'
        ! Offset the pressure and magnetic interior, keeping boundary data.
        where (.not. nc_fixed(field)) x(:lay%n) = 0.9_dp*x
        call nc_solve(field, problem, x, report, st)
        if (st /= 0) error stop 'nested solve'
        if (maxval(report%residual) > 1.0e-11_dp) error stop 'nested residual'
        if (report%divergence > 1.0e-12_dp) error stop 'nested divergence'
        do i = 1, 64
        call mapping(0.2_dp + 0.65_dp*mod(i, 8)/8, 2*pi*(i - 0.31_dp)/64, points(:2, i))
            points(3, i) = 0.7_dp
            call exact(points(:2, i), psi, b)
            bref(:, i) = [b(1), ftor/points(1, i), b(2)]
        end do
        call field%sample(points, values, sts)
        if (any(sts /= 0)) error stop 'nested sample'
        call target%evaluate(points, pe, gp, lam, gf, st)
        ep2 = 0; np2 = 0
        do t = 1, field%fs%mm%topology%n_triangles
            do q = 1, field%fs%nq
                call meridian_qpoint(field%fs, 0, t, q, xy, wq, &
                                     basis(:, :field%fs%nl0), st)
                if (st /= 0) error stop 'nested pressure quadrature'
                call exact(xy, psi, b)
            ph = dot_product(basis(1, :field%fs%nl0), field%p(field%fs%sp(0)%dof(:, t)))
                ep2 = ep2 + wq*(ph - pp*psi)**2
                np2 = np2 + wq*(pp*psi)**2
            end do
        end do
        error = sqrt(ep2/np2)
        berror = sqrt(sum((values(EQ_BR, :) - bref(1, :))**2 &
               + (values(EQ_BZ, :) - bref(3, :))**2)/sum(bref(1, :)**2 + bref(3, :)**2))
        call solved_target%init(field, x(lay%iphi + 1:lay%iphi + lay%nv), &
                                x(lay%iv), problem%g_par, problem%g_perp, st)
        if (st /= 0) error stop 'solved profile adapter'
        call solved_target%evaluate(points, pe, gp, lam, gf, st)
        if (st /= 0) error stop 'solved gradients'
        weights = 1
     call problem%fitted%defects(field, points, weights, gp, gf, x(lay%iv), defects, st)
        if (st /= 0) error stop 'nested defects'
        if (defects%seed_current_fit_rms > 1.0e-12_dp) error stop 'seed current fit'
        print '(a,2i5,4es13.4)', 'nested P/nr p/B/lambda/seed-fit ', order, nr, &
            error, berror, defects%lambda_rms, defects%seed_current_fit_rms
        if (nr == 4) call jacobian_check(field, problem, x)
    end subroutine

    subroutine reproduce(order, three_d)
        integer, intent(in) :: order
        logical, intent(in) :: three_d
        type(free_b_field_t) :: field
        type(nn_state_profile_t) :: target
        type(nc_layout_t) :: lay
        type(nc_problem_t) :: problem
        type(nc_report_t) :: report
        real(dp), allocatable :: nodes(:, :), x(:), reference(:), r(:)
        real(dp) :: defect
        integer :: st
        call setup(merge(2, 4, three_d), order, field, lay, nodes, x, problem, &
                   modes=merge(1, 0, three_d))
        if (three_d) then
            problem%source_3d => old_source_3d
        else
            problem%source => old_source
        end if
        call nc_solve(field, problem, x, report, st, tol=1.0e-13_dp)
        if (st /= 0) error stop 'A+D solve'
        reference = x
        call target%init(field, x(lay%iphi + 1:lay%iphi + lay%nv), x(lay%iv), &
                         1.0_dp, 1.0_dp, st)
        if (st /= 0) error stop 'A+D profile adapter'
        problem%closure%pressure = 'B_frozen_fit'
        allocate (problem%fitted)
        call nn_fit_source(field, target, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, &
                           x(lay%iv), problem%fitted, st)
        if (st /= 0) error stop 'A+D fit'
        nullify (problem%source, problem%source_3d)
        call nc_residual(field, problem, x, r, st)
        if (st /= 0) error stop 'A+D fitted residual'
        print '(a,es14.5)', 'A+D fitted residual at original state ', norm2(r)
        where (.not. nc_fixed(field)) x = 0.97_dp*x
        call nc_solve(field, problem, x, report, st, tol=1.0e-13_dp)
        if (st /= 0) error stop 'A+D fitted solve'
        defect = norm2(x - reference)/norm2(reference)
        print '(a,es14.5)', 'fitted A+D solution return ', defect
        if (defect > 1.0e-8_dp) error stop 'A+D reproduction'
        call jacobian_check(field, problem, x)
        ! Reject ambiguous source declarations instead of silently ignoring S.
        problem%source => old_source
        call nc_solve(field, problem, x, report, st)
        if (st == 0) error stop 'ambiguous fitted source accepted'
    end subroutine

    subroutine sampling_check(order)
        integer, intent(in) :: order
        type(free_b_field_t) :: field
        type(nc_problem_t) :: problem
        type(nc_layout_t) :: lay
        real(dp), allocatable :: nodes(:, :), state(:)
        real(dp) :: point(3, 1), values(EQ_VALUES, 1), xy(2), wq
        real(dp) :: basis(3, 64), b(3), expected, error
        integer :: t, q, st, statuses(1)
        call setup(16, order, field, lay, nodes, state, problem)
        call nc_store(field, state, st)
        if (st /= 0) error stop 'sampling store'
        error = 0
        do t = 1, field%fs%mm%topology%n_triangles
            do q = 1, field%fs%nq
                call meridian_qpoint(field%fs, 2, t, q, xy, wq, &
                    basis(:, :field%fs%sp(2)%nloc), st)
                if (st /= 0) error stop 'sampling element field'
                b = real(matmul(basis(:, :field%fs%sp(2)%nloc), &
                    field%b(field%fs%sp(2)%dof(:, t), 0) &
                    *field%fs%sp(2)%sgn(:, t)), dp)
                point(:, 1) = [xy, 0.37_dp]
                call field%sample(point, values, statuses)
                if (any(statuses /= 0)) error stop 'sampling physical point'
                error = max(error, abs(values(EQ_BR, 1) - b(1)), &
                    abs(values(EQ_BZ, 1) - b(3)))
                call meridian_qpoint(field%fs, 0, t, q, xy, wq, &
                    basis(:, :field%fs%nl0), st)
                if (st /= 0) error stop 'sampling element pressure'
                expected = dot_product(basis(1, :field%fs%nl0), &
                    field%p(field%fs%sp(0)%dof(:, t)))
                error = max(error, abs(values(EQ_P, 1) - expected))
            end do
        end do
        print '(a,es14.5)', 'represented-cell / physical sampling defect ', error
        if (error > 1.0e-9_dp) error stop 'represented-cell sampling'
    end subroutine

    subroutine jacobian_check(field, problem, x)
        type(free_b_field_t), intent(in) :: field
        type(nc_problem_t), intent(in) :: problem
        real(dp), intent(in) :: x(:)
        type(csc_t) :: jac
        real(dp), allocatable :: d(:), rp(:), rm(:), jd(:)
        real(dp) :: h, error
        integer :: i, st
        d = [(sin(1.3_dp*i), i=1, size(x))]
        d = d*maxval(abs(x))*0.01_dp
        where (nc_fixed(field)) d = 0
        h = 1.0e-5_dp
        call nc_jacobian(field, problem, x, jac, st)
        if (st /= 0) error stop 'fitted Jacobian'
        jd = csc_matvec(jac, d)
        where (nc_fixed(field)) jd = 0
        call nc_residual(field, problem, x + h*d, rp, st)
        if (st /= 0) error stop 'fitted residual plus'
        call nc_residual(field, problem, x - h*d, rm, st)
        if (st /= 0) error stop 'fitted residual minus'
        error = norm2((rp - rm)/(2*h) - jd)/norm2(jd)
        print '(a,es14.5)', 'frozen-source Jacobian directional error ', error
        if (error > 1.0e-7_dp) error stop 'fitted Jacobian consistency'
    end subroutine
end program
