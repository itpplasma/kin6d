module nn_anisotropy_fixture
    ! Independent analytical data for the Sovinec diagnostic, not production leaves.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_equilibrium_state, only: equilibrium_field_t, EQ_VALUES, EQ_BR, EQ_BZ
    use fortfem_meridian_fourier_scalar, only: meridian_fourier_source_t
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp)
    type, extends(equilibrium_field_t) :: field_t
        real(dp) :: amplitude = 1
    contains
        procedure :: sample
    end type
    type, extends(meridian_fourier_source_t) :: source_t
        logical :: weak = .false.
    contains
        procedure :: source
    end type
contains
    subroutine sample(field, points, values, statuses)
        class(field_t), intent(in) :: field
        real(dp), intent(in) :: points(:, :)
        real(dp), intent(out) :: values(:, :)
        integer, intent(out) :: statuses(:)
        values = 0
        statuses = 0
        values(EQ_BR, :) = field%amplitude*pi*cos(pi*points(1, :)) &
            *sin(pi*points(2, :))
        values(EQ_BZ, :) = -field%amplitude*pi*sin(pi*points(1, :)) &
            *cos(pi*points(2, :))
    end subroutine

    subroutine source(self, t, x, f0, f, status)
        class(source_t), intent(in) :: self
        integer, intent(in) :: t
        real(dp), intent(in) :: x(:, :)
        real(dp), intent(out) :: f0(:), f(:, :)
        integer, intent(out) :: status
        f0 = 0
        f = 0
        status = 1
        if (t < 1) return
        if (self%weak) then
            ! b.grad(psi) = 0, so K grad(psi) = grad(psi) at every ratio.
            f(1, :) = -pi*sin(pi*x(1, :))*cos(pi*x(2, :))
            f(3, :) = -pi*cos(pi*x(1, :))*sin(pi*x(2, :))
        else
            f0 = 2*pi*pi*cos(pi*x(1, :))*cos(pi*x(2, :))
        end if
        status = 0
    end subroutine

    subroutine square_map(u, theta, x)
        real(dp), intent(in) :: u, theta
        real(dp), intent(out) :: x(2)
        x = 0.5_dp*u*[cos(theta), sin(theta)] &
            /max(abs(cos(theta)), abs(sin(theta)))
    end subroutine
end module

program test_nonnested_anisotropy
    use omp_lib, only: omp_get_wtime
    use nn_anisotropy_fixture, only: dp, pi, field_t, source_t, square_map
    use kin6d_mhd_mesh, only: flux_mesh_t, flux_mesh_build
    use kin6d_nonnested_transport, only: nt_conductivity_t, nt_report_t, nt_solve
    use fortfem_meridian_feec_n, only: meridian_feec_t, meridian_feec_build, &
        meridian_qgrad, meridian_lagrange_nodes
    use fortfem_meridian_curved_mesh, only: meridian_geometry_point
    use fortfem_triangle_duffy_quadrature, only: triangle_duffy_quadrature
    use fortfem_equilibrated_flux_2d, only: equilibrated_flux_estimate_tensor, &
        equilibrated_flux_estimate_t, EQUILIBRATION_OK
    implicit none
    type(meridian_feec_t) :: fs
    type(field_t), target :: field
    type(nt_conductivity_t) :: cond
    type(source_t) :: src
    type(nt_report_t) :: report
    type(equilibrated_flux_estimate_t) :: baseline, corrected
    real(dp), allocatable :: p(:), weak(:), reference(:), bc(:), xy(:, :)
    real(dp), allocatable :: tensor(:, :, :, :), load(:, :), flux(:, :, :)
    real(dp), allocatable :: vertices(:, :), xi(:), eta(:), w(:), lower(:)
    integer, allocatable :: triangles(:, :)
    real(dp) :: ratio, started, solve_time, base_time, correction_time, error
    real(dp) :: x(2), jac(2, 2), val(6), grad(2, 6), k(3, 3, 1), wq, gp(2)
    real(dp) :: energy, work, discrepancy, quad_change
    real(dp) :: previous(4), backward
    integer :: nr, lev, ir, st, nt, centre(1), sweeps
    character(32) :: mode, arg
    logical :: ok

    interface
        subroutine dposv(uplo, n, nrhs, a, lda, rhs, ldb, info)
            import dp
            character, intent(in) :: uplo
            integer, intent(in) :: n, nrhs, lda, ldb
            real(dp), intent(inout) :: a(lda, *), rhs(ldb, *)
            integer, intent(out) :: info
        end subroutine
    end interface
    call get_command_argument(1, mode)
    nr = 0
    sweeps = 8
    call get_command_argument(2, arg)
    if (len_trim(arg) > 0) read (arg, *) nr
    call get_command_argument(3, arg)
    if (len_trim(arg) > 0) read (arg, *) sweeps
    ok = .true.
    cond%field => field
    previous = huge(1.0_dp)
    do lev = 1, merge(1, 2, nr > 0)
        if (nr == 0) then
            call build(4*2**(lev - 1), 8, fs)
        else
            call build(nr, 8, fs)
        end if
        allocate (bc(fs%n0))
        bc = 0
        call meridian_lagrange_nodes(fs, xy)
        centre = minloc(norm2(xy, dim=1))
        do ir = 1, 4
            ratio = 10.0_dp**(3*(ir - 1))
            cond%chi_par = ratio
            src%weak = .false.
            started = omp_get_wtime()
            call nt_solve(fs, 0, 1, cond, src, bc, p, report, st)
            solve_time = omp_get_wtime() - started
            if (st /= 0) error stop 'Sovinec production solve'
            ok = ok .and. report%residual < 1.0e-13_dp
            backward = report%residual
            if (trim(mode) == 'pollution') then
                src%weak = .true.
                call nt_solve(fs, 0, 1, cond, src, bc, weak, report, st)
                if (st /= 0) error stop 'Sovinec weak-load solve'
                call dense_reference(fs, ratio, reference, error, energy, work)
                discrepancy = maxval(abs(p - reference))
                ok = ok .and. discrepancy < 2.0e-9_dp
                ok = ok .and. maxval(abs(p - weak)) < 2.0e-9_dp
                ok = ok .and. abs(error**2 - (pi*pi/2 - 2*work + energy)) &
                    < 2.0e-8_dp
                ! The high-ratio state remains polluted under oversampling.
                call dense_oversampled(fs, ratio, p, quad_change)
                if (ir >= 3) then
                    ok = ok .and. abs(p(centre(1)) - 1) > 0.9_dp
                    ok = ok .and. error > 2.0_dp
                    ok = ok .and. quad_change < 1.0e-2_dp
                end if
                print '(a,i0,a,es9.2,6(a,es11.3))', 'pollution dofs=', fs%n0, &
                    ' ratio=', ratio, ' centre_error=', abs(p(centre(1)) - 1), &
                    ' energy_error=', error, ' dense_delta=', discrepancy, &
                    ' weak_delta=', maxval(abs(p - weak)), ' q24_delta=', &
                    quad_change, ' backward=', backward
            else if (trim(mode) == 'correction' .or. trim(mode) == 'study') then
                call samples(fs, p, tensor, load, flux, vertices, triangles)
                call triangle_duffy_quadrature(8, xi, eta, w, st)
                nt = size(triangles, 2)
                allocate (lower(nt))
                lower = 1
                started = omp_get_wtime()
                call equilibrated_flux_estimate_tensor(vertices, triangles, 2, &
                    xi, eta, w, tensor, load, flux, lower, 1/(sqrt(2.0_dp)*pi), &
                    baseline, st)
                if (st /= EQUILIBRATION_OK) error stop 'baseline reconstruction'
                base_time = omp_get_wtime() - started
                started = omp_get_wtime()
                call equilibrated_flux_estimate_tensor(vertices, triangles, 2, &
                    xi, eta, w, tensor, load, flux, lower, 1/(sqrt(2.0_dp)*pi), &
                    corrected, st, correction_sweeps=sweeps)
                if (st /= EQUILIBRATION_OK) error stop 'corrected reconstruction'
                correction_time = omp_get_wtime() - started
                call energy_error(fs, p, ratio, error)
                ok = ok .and. corrected%bound >= error
                ok = ok .and. corrected%bound <= baseline%bound*(1 + 1.0e-10_dp)
                if (trim(mode) == 'correction') &
                    ok = ok .and. corrected%bound/error < 5.0_dp
                ok = ok .and. corrected%divergence_defect < 1.0e-5_dp
                if (ir <= 2) then
                    if (lev > 1) then
                        ! h halves: bound reduction > 1.8 gives rate > 0.84.
                        ok = ok .and. previous(ir) > 1.8_dp*corrected%bound
                    end if
                end if
                previous(ir) = corrected%bound
                print '(a,i0,a,es9.2,7(a,es11.3))', 'correction dofs=', fs%n0, &
                    ' ratio=', ratio, ' error=', error, ' index=', &
                    corrected%bound/error, ' baseline=', baseline%bound/error, &
                    ' solve_s=', solve_time, ' baseline_s=', base_time, &
                    ' corrected_s=', correction_time, ' cost/solve=', &
                    correction_time/solve_time
                if (trim(mode) == 'study') print '(a,2es24.16)', &
                    'invariance bound/error ', corrected%bound, error
                deallocate (lower)
            else
                error stop 'expected pollution, correction or study'
            end if
        end do
        deallocate (bc)
    end do
    if (.not. ok) error stop 'anisotropy diagnostic failed'
    print '(a)', 'PASS'
contains
    subroutine build(rings_count, qdeg, space)
        integer, intent(in) :: rings_count, qdeg
        type(meridian_feec_t), intent(out) :: space
        type(flux_mesh_t) :: mesh
        real(dp) :: rings(rings_count)
        integer :: i, status
        rings = [(real(i, dp)/rings_count, i=1, rings_count)]
        call flux_mesh_build(rings, 8, 4096, 1, square_map, mesh, status)
        if (status /= 0) error stop 'Sovinec mesh'
        call meridian_feec_build(mesh%mm, 2, space, status, quad_degree=qdeg, &
            cylinder_r0=1.0_dp)
        if (status /= 0) error stop 'Sovinec space'
    end subroutine

    subroutine dense_reference(space, chi, result, err, eh, fh)
        ! Independent directional Gram assembly and LAPACK Cholesky solve.
        ! Handwritten analytical contractions are deliberate test oracles.
        type(meridian_feec_t), intent(in) :: space
        real(dp), intent(in) :: chi
        real(dp), allocatable, intent(out) :: result(:)
        real(dp), intent(out) :: err, eh, fh
        real(dp), allocatable :: matrix(:, :), rhs(:, :), original(:, :), f(:)
        real(dp) :: xx(2), ww, vv(6), gg(2, 6), bb(2), par(6)
        integer :: cell, qp, i, j, ii, jj, info
        allocate (matrix(space%n0, space%n0), rhs(space%n0, 1))
        matrix = 0
        rhs = 0
        do cell = 1, space%mm%topology%n_triangles
            do qp = 1, space%nq
                call meridian_qgrad(space, cell, qp, xx, ww, vv, gg, info)
                if (info /= 0) error stop 'Gram quadrature'
                bb = [cos(pi*xx(1))*sin(pi*xx(2)), &
                    -sin(pi*xx(1))*cos(pi*xx(2))]
                bb = bb/norm2(bb)
                par = matmul(bb, gg)
                do i = 1, 6
                    ii = space%sp(0)%dof(i, cell)
                    rhs(ii, 1) = rhs(ii, 1) &
                        + ww*2*pi*pi*cos(pi*xx(1))*cos(pi*xx(2))*vv(i)
                    do j = 1, 6
                        jj = space%sp(0)%dof(j, cell)
                        matrix(ii, jj) = matrix(ii, jj) + ww &
                            *(dot_product(gg(:, i), gg(:, j)) &
                            + (chi - 1)*par(i)*par(j))
                    end do
                end do
            end do
        end do
        original = matrix
        f = rhs(:, 1)
        do i = 1, space%n0
            if (.not. space%sp(0)%bnd(i)) cycle
            matrix(:, i) = 0
            matrix(i, :) = 0
            matrix(i, i) = 1
            rhs(i, 1) = 0
        end do
        call dposv('L', space%n0, 1, matrix, space%n0, rhs, space%n0, info)
        if (info /= 0) error stop 'independent Cholesky'
        result = rhs(:, 1)
        eh = dot_product(result, matmul(original, result))
        fh = dot_product(f, result)
        call energy_error(space, result, chi, err)
    end subroutine

    subroutine dense_oversampled(space, chi, candidate, change)
        type(meridian_feec_t), intent(in) :: space
        real(dp), intent(in) :: chi, candidate(:)
        real(dp), intent(out) :: change
        type(meridian_feec_t) :: over
        real(dp), allocatable :: result(:)
        real(dp) :: err, eh, fh
        integer :: status
        call meridian_feec_build(space%mm, 2, over, status, quad_degree=24, &
            cylinder_r0=1.0_dp)
        if (status /= 0) error stop 'oversampled space'
        call dense_reference(over, chi, result, err, eh, fh)
        change = maxval(abs(result - candidate))
    end subroutine

    subroutine energy_error(space, candidate, chi, err)
        type(meridian_feec_t), intent(in) :: space
        real(dp), intent(in) :: candidate(:), chi
        real(dp), intent(out) :: err
        real(dp) :: xx(2), ww, vv(6), gg(2, 6), ge(2), bb(2), ee(2), squared
        integer :: cell, qp, info
        squared = 0
        !$omp parallel do private(cell, qp, info, xx, ww, vv, gg, ge, bb, ee) &
        !$omp reduction(+:squared) schedule(static)
        do cell = 1, space%mm%topology%n_triangles
            do qp = 1, space%nq
                call meridian_qgrad(space, cell, qp, xx, ww, vv, gg, info)
                ge = -pi*[sin(pi*xx(1))*cos(pi*xx(2)), &
                    cos(pi*xx(1))*sin(pi*xx(2))]
                bb = [ge(2), -ge(1)]/norm2(ge)
                ee = ge - matmul(gg, candidate(space%sp(0)%dof(:, cell)))
                squared = squared + ww*(dot_product(ee, ee) &
                    + (chi - 1)*dot_product(bb, ee)**2)
            end do
        end do
        !$omp end parallel do
        err = sqrt(squared)
    end subroutine

    subroutine samples(space, candidate, aa, ff, qq, points, cells)
        type(meridian_feec_t), intent(in) :: space
        real(dp), intent(in) :: candidate(:)
        real(dp), allocatable, intent(out) :: aa(:, :, :, :), ff(:, :), qq(:, :, :)
        real(dp), allocatable, intent(out) :: points(:, :)
        integer, allocatable, intent(out) :: cells(:, :)
        real(dp) :: corner(2, 3)
        integer :: cell, j, qp, info, count
        count = space%mm%topology%n_triangles
        cells = space%mm%topology%triangles
        allocate (points(2, space%mm%topology%n_vertices))
        allocate (aa(2, 2, space%nq, count), ff(space%nq, count))
        allocate (qq(2, space%nq, count))
        corner = reshape([0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2, 3])
        do cell = 1, count
            do j = 1, 3
                call meridian_geometry_point(space%mm, cell, corner(1, j), &
                    corner(2, j), x, jac)
                points(:, cells(j, cell)) = x
            end do
        end do
        !$omp parallel do private(cell, qp, info, x, wq, val, grad, k, gp) &
        !$omp schedule(static)
        do cell = 1, count
            do qp = 1, space%nq
                call meridian_qgrad(space, cell, qp, x, wq, val, grad, info)
                call cond%tensor(cell, reshape([x, 0.0_dp], [3, 1]), k, info)
                aa(:, :, qp, cell) = k([1, 3], [1, 3], 1)
                gp = matmul(grad, candidate(space%sp(0)%dof(:, cell)))
                qq(:, qp, cell) = matmul(aa(:, :, qp, cell), gp)
                ff(qp, cell) = 2*pi*pi*cos(pi*x(1))*cos(pi*x(2))
            end do
        end do
        !$omp end parallel do
    end subroutine
end program
