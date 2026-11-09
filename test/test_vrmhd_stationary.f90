module vrmhd_spd_extremes
    !! Extreme eigenpairs of dense SPD matrices without a full decomposition.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_krylov, only: real_linear_operator_t, KRYLOV_OK
    use fortnum_symmetric_eigen, only: block_lanczos_lowest
    implicit none
    private
    public :: spd_extremes
    type, extends(real_linear_operator_t) :: spd_operator_t
        real(dp), allocatable :: a(:, :), factor(:, :)
        logical :: inverse = .false.
    contains
        procedure :: apply => spd_apply
    end type
    interface
        subroutine dpotrf(uplo, n, a, lda, info)
            import :: dp
            character, intent(in) :: uplo
            integer, intent(in) :: n, lda
            real(dp), intent(inout) :: a(lda, *)
            integer, intent(out) :: info
        end subroutine
        subroutine dpotrs(uplo, n, nrhs, a, lda, b, ldb, info)
            import :: dp
            character, intent(in) :: uplo
            integer, intent(in) :: n, nrhs, lda, ldb
            real(dp), intent(in) :: a(lda, *)
            real(dp), intent(inout) :: b(ldb, *)
            integer, intent(out) :: info
        end subroutine
        subroutine dsymv(uplo, n, alpha, a, lda, x, incx, beta, y, incy)
            import :: dp
            character, intent(in) :: uplo
            integer, intent(in) :: n, lda, incx, incy
            real(dp), intent(in) :: alpha, beta, a(lda, *), x(*)
            real(dp), intent(inout) :: y(*)
        end subroutine
    end interface
contains
    subroutine spd_apply(self, input, output)
        class(spd_operator_t), intent(inout) :: self
        real(dp), intent(in) :: input(:)
        real(dp), intent(out) :: output(:)
        integer :: n, info
        n = size(input)
        ! Lanczos returns the lowest pairs, so both maps are negated.
        if (self%inverse) then
            output = -input
            call dpotrs('L', n, 1, self%factor, n, output, n, info)
            if (info /= 0) output = huge(1.0_dp)
        else
            call dsymv('L', n, -1.0_dp, self%a, n, input, 1, 0.0_dp, output, 1)
        end if
    end subroutine

    subroutine spd_extremes(a, lowest, highest, vector, ok)
        !! Lowest eigenpair and highest eigenvalue of the SPD matrix a (lower
        !! triangle). Shift-invert Lanczos gives residual-bounded values;
        !! Cholesky factorizations of a - (1 - margin) lowest I and of
        !! (1 + margin) highest I - a check that no eigenvalue lies outside.
        real(dp), intent(in) :: a(:, :)
        real(dp), intent(out) :: lowest, highest, vector(:)
        logical, intent(out) :: ok
        real(dp), parameter :: tolerance = 1.0e-10_dp, margin = 1.0e-6_dp
        type(spd_operator_t) :: op
        real(dp) :: start(size(a, 1), 1), values(1), vectors(size(a, 1), 1)
        real(dp) :: residuals(1)
        integer :: n, i, info, iterations
        n = size(a, 1); ok = .false.
        op%a = a; op%factor = a
        call dpotrf('L', n, op%factor, n, info)
        if (info /= 0) return
        start(:, 1) = [(sin(real(i, dp)), i=1, n)]
        op%inverse = .true.
        call block_lanczos_lowest(op, start, 1, min(n, 400), tolerance, values, &
            vectors, residuals, info, iterations)
        if (info /= KRYLOV_OK) return
        lowest = -1/values(1); vector = vectors(:, 1)
        op%inverse = .false.
        call block_lanczos_lowest(op, start, 1, min(n, 400), tolerance, values, &
            vectors, residuals, info, iterations)
        if (info /= KRYLOV_OK) return
        highest = -values(1)
        op%factor = a
        do i = 1, n
            op%factor(i, i) = op%factor(i, i) - (1 - margin)*lowest
        end do
        call dpotrf('L', n, op%factor, n, info)
        if (info /= 0) return
        op%factor = -a
        do i = 1, n
            op%factor(i, i) = op%factor(i, i) + (1 + margin)*highest
        end do
        call dpotrf('L', n, op%factor, n, info)
        ok = info == 0
    end subroutine
end module vrmhd_spd_extremes

program test_vrmhd_stationary
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use milca_gs_symmetry, only: gs_symmetry_t, GS_ROTATION, GS_SCREW, &
        GS_CHART_RZ, GS_CHART_CARTESIAN
    use milca_vrmhd_parameters, only: vrmhd_parameters_t
    use milca_vrmhd_section
    use milca_vrmhd_stationary
    use vrmhd_wall_oracle, only: current_wall_dofs
    use milca_vrmhd_boundary, only: vrmhd_current_lift
    use milca_vrmhd_drive, only: vrmhd_drive_build, VRMHD_CURRENT_SOURCE, &
        VRMHD_LOOP_VOLTAGE
    use vrmhd_mms_leaves, only: mms_rotation, mms_screw
    use vrmhd_spd_extremes, only: spd_extremes
    use fortfem_mesh_2d, only: mesh_2d_t
    use fortfem_meridian_curved_mesh, only: meridian_mesh_t, meridian_mesh_straight
    use fortsparse, only: csc_t, csc_from_triplet, csc_matvec, sparse_solve_once, &
        fortsparse_status_t, status_ok
    implicit none
    character(32) :: gate
    integer :: kind, level
    real(dp) :: errors(8, 3), rates(8), curl_errors(2, 3), curl_rates(2)
    real(dp) :: sigma(3), mass_condition(3)
    call get_command_argument(1, gate)
    do kind = 1, 3
        select case (trim(gate))
        case ('boundary')
            call check_boundary(kind)
        case ('derivatives')
            call check_derivatives(kind)
        case ('harmonic')
            call check_harmonic(kind)
        case ('gauge')
            do level = 1, 3
                call check_gauge(kind, 2**(level + 1), sigma(level), &
                    mass_condition(level))
                if (level == 1) cycle
                print '(a,2i3,f10.6)', 'reduced minimum refinement ratio ', &
                    kind, level, sigma(level)/sigma(level - 1)
                call require(sigma(level) >= 0.75_dp*sigma(level - 1), &
                    'reduced stationary minimum collapses under refinement')
                print '(a,2i3,f10.6)', 'scaled current mass condition refinement ratio ', &
                    kind, level, mass_condition(level)/mass_condition(level - 1)
                call require( &
                    mass_condition(level) <= mass_condition(level - 1)/0.75_dp, &
                    'scaled current mass condition grows under refinement')
            end do
        case ('solve')
            do level = 1, 3
                call solve_manufactured(kind, 2**(level + 1), errors(:, level), &
                    curl_errors(:, level))
                print '(a,2i3,8es13.5)', 'solved natural-norm errors ', &
                    kind, level, errors(:, level)
                if (level == 1) cycle
                rates = log(errors(:, level - 1)/errors(:, level))/log(2.0_dp)
                print '(a,2i3,8f8.3)', 'solved natural-norm rates ', kind, level, rates
                curl_rates(1) = log(curl_errors(1, level - 1) &
                    /curl_errors(1, level))/log(2.0_dp)
                print '(a,2i3,f8.3)', 'diagnostic curl J all-cell rate ', &
                    kind, level, curl_rates(1)
                if (level > 2) then
                    curl_rates(2) = log(curl_errors(2, level - 1) &
                        /curl_errors(2, level))/log(2.0_dp)
                    print '(a,2i3,f8.3)', 'diagnostic curl J interior rate ', &
                        kind, level, curl_rates(2)
                    if (curl_rates(2) <= curl_rates(1) + 0.2_dp) &
                        print '(a)', 'OPEN: interior curl J rate is not clearly better'
                end if
                call require(minval(rates) > 0.8_dp, 'P1 solved design order')
            end do
        case default
            error stop 'unknown stationary gate'
        end select
    end do
contains
    subroutine check_gauge(kind, n, sigma_min, mass_condition)
        integer, intent(in) :: kind, n
        real(dp), intent(out) :: sigma_min, mass_condition
        type(vrmhd_section_t) :: op
        type(vrmhd_stationary_t) :: steady
        type(vrmhd_parameters_t) :: p
        type(csc_t) :: jac
        type(fortsparse_status_t) :: st
        real(dp), allocatable :: u(:), load(:), source(:, :, :), basis(:, :)
        real(dp), allocatable :: boundary(:), r(:), zero(:), &
            step(:), unused(:)
        real(dp), allocatable :: trial(:, :), matrix(:, :), minimum_mode(:)
        real(dp), allocatable :: block(:, :), current_mass(:, :), &
            mass_copy(:, :), current_primal(:, :), &
            primal_current(:, :), recovery(:, :), mass_mode(:), &
            probe(:), image(:), mass_diagonal(:)
        integer, allocatable :: free_index(:), gauge_index(:), current_index(:)
        real(dp) :: x(2), weight, value, mode_norms(8), defect, scale
        real(dp) :: raw_condition, lowest, highest
        logical :: certified
        real(dp), allocatable :: gauge_value(:), gauge_gradient(:, :)
        integer :: i, j, t, q, nf, ng, row, col, status, entry, info
        integer :: group, it, nj, current_start, nh
        interface
            subroutine dpotrf(uplo, n, a, lda, info)
                import :: dp
                character, intent(in) :: uplo
                integer, intent(in) :: n, lda
                real(dp), intent(inout) :: a(lda, *)
                integer, intent(out) :: info
            end subroutine
            subroutine dpotrs(uplo, n, nrhs, a, lda, b, ldb, info)
                import :: dp
                character, intent(in) :: uplo
                integer, intent(in) :: n, nrhs, lda, ldb
                real(dp), intent(in) :: a(lda, *)
                real(dp), intent(inout) :: b(ldb, *)
                integer, intent(out) :: info
            end subroutine
            subroutine dtrsm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
                import :: dp
                character, intent(in) :: side, uplo, transa, diag
                integer, intent(in) :: m, n, lda, ldb
                real(dp), intent(in) :: alpha, a(lda, *)
                real(dp), intent(inout) :: b(ldb, *)
            end subroutine
            subroutine dgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, &
                    c, ldc)
                import :: dp
                character, intent(in) :: transa, transb
                integer, intent(in) :: m, n, k, lda, ldb, ldc
                real(dp), intent(in) :: alpha, beta, a(lda, *), b(ldb, *)
                real(dp), intent(inout) :: c(ldc, *)
            end subroutine
            subroutine dsyrk(uplo, trans, n, k, alpha, a, lda, beta, c, ldc)
                import :: dp
                character, intent(in) :: uplo, trans
                integer, intent(in) :: n, k, lda, ldc
                real(dp), intent(in) :: alpha, beta, a(lda, *)
                real(dp), intent(inout) :: c(ldc, *)
            end subroutine
        end interface
        call prepare(kind, n, op, steady, p, u, load)
        allocate (gauge_value(op%fs%nl0), gauge_gradient(3, op%fs%nl0))
        allocate (source(6, op%fs%nq, op%fs%mm%topology%n_triangles), &
            free_index(steady%ndof), gauge_index(op%fs%n0), &
            current_index(steady%ndof))
        source = 0; free_index = 0; gauge_index = 0; current_index = 0
        nf = 0; ng = 0; nj = 0
        current_start = 6*op%scalar_fs%n0 + op%fs%sp(1)%ndof
        allocate (boundary(op%ndof), r(steady%ndof), &
            zero(steady%ndof), step(steady%ndof), unused(steady%ndof))
        boundary = u(:op%ndof); zero = 0

        do it = 1, 10
            call vrmhd_stationary_residual(op, steady, p, u, zero, boundary, &
                source, load, r, unused, status)
            call require(status == 0, 'gauge stationary Newton residual')
            if (norm2(r) < 1.0e-10_dp) exit
            call vrmhd_stationary_jacobian(op, steady, p, u, boundary, source, &
                jac, status)
            call require(status == 0, 'gauge stationary Newton Jacobian')
            call sparse_solve_once(jac, -r, step, st)
            call require(status_ok(st), 'gauge stationary Newton direct solve')
            u = u + step
        end do
        call require(it <= 10, 'gauge probe at converged stationary solution')
        print '(a,2i5,es14.6)', 'gauge probe stationary residual ', kind, n, norm2(r)
        do i = 1, op%ndof
            if (op%fixed(i)) cycle
            if (i > current_start) then
                nj = nj + 1
                current_index(i) = nj
                cycle
            end if
            nf = nf + 1
            free_index(i) = nf
        end do
        do i = 1, op%fs%n0
            if (op%fs%sp(0)%bnd(i)) cycle
            ng = ng + 1
            gauge_index(i) = nf + ng
            free_index(op%ndof + ng) = nf + ng
        end do
        nf = nf + ng
        allocate (trial(nf, nf), matrix(nf, nf), minimum_mode(nf), &
            basis(36, op%nlocal), block(op%nlocal, op%nlocal), &
            current_mass(nj, nj), mass_copy(nj, nj), &
            current_primal(nj, nf), primal_current(nf, nj), &
            recovery(nj, nf), mass_mode(nj), &
            probe(nf), image(nf), mass_diagonal(nj))
        trial = 0; matrix = 0; current_mass = 0
        current_primal = 0; primal_current = 0
        call vrmhd_stationary_jacobian(op, steady, p, u, boundary, source, &
            jac, status)
        call require(status == 0, 'gauge nonsingularity Jacobian')
        do j = 1, steady%ndof
            do entry = jac%col_ptr(j), jac%col_ptr(j + 1) - 1
                i = jac%row_idx(entry)
                row = free_index(i); col = free_index(j)
                if (row > 0 .and. col > 0) matrix(row, col) = jac%val(entry)
                if (row > 0 .and. current_index(j) > 0) &
                    primal_current(row, current_index(j)) = jac%val(entry)
                if (current_index(i) > 0 .and. col > 0) &
                    current_primal(current_index(i), col) = jac%val(entry)
                if (current_index(i) > 0 .and. current_index(j) > 0) &
                    current_mass(current_index(i), current_index(j)) = jac%val(entry)
            end do
        end do
        ! Eliminate algebraic J before measuring stability of the primal fields.
        mass_copy = current_mass
        recovery = current_primal
        call dpotrf('L', nj, mass_copy, nj, info)
        call require(info == 0, 'current mass Cholesky')
        call dpotrs('L', nj, nf, mass_copy, nj, recovery, nj, info)
        call require(info == 0, 'current elimination solve')
        ! trial is free until the Riesz assembly and holds M_J X - C here.
        trial(:nj, :) = current_primal
        call dgemm('N', 'N', nj, nf, nj, 1.0_dp, current_mass, nj, recovery, nj, &
            -1.0_dp, trial, nf)
        call require(maxval(abs(trial(:nj, :))) &
            < 1.0e-10_dp*max(1.0_dp, maxval(abs(current_primal))), &
            'current elimination original-block residual')
        trial = 0
        call dgemm('N', 'N', nf, nf, nj, -1.0_dp, primal_current, nf, recovery, nj, &
            1.0_dp, matrix, nf)
        do i = 1, nf
            probe(i) = sin(real(i, dp))
        end do
        r(:nj) = -matmul(recovery, probe)
        step = 0
        do i = 1, steady%ndof
            if (free_index(i) > 0) step(i) = probe(free_index(i))
            if (current_index(i) > 0) step(i) = r(current_index(i))
        end do
        unused = csc_matvec(jac, step)
        image = matmul(matrix, probe)
        defect = 0
        do i = 1, steady%ndof
            if (free_index(i) > 0) defect = max(defect, &
                abs(unused(i) - image(free_index(i))))
            if (current_index(i) > 0) defect = max(defect, abs(unused(i)))
        end do
        scale = max(1.0_dp, maxval(abs(image)))
        print '(a,2i5,es14.6)', 'reduced versus full Jacobian action ', &
            kind, n, defect/scale
        call require(defect < 1.0e-10_dp*scale, 'reduced action and full sparse oracle')
        ! A symmetric positive mass has its eigenvalues as singular values.
        call require(all(current_mass == transpose(current_mass)), &
            'symmetric current mass')
        call spd_extremes(current_mass/p%mu0, lowest, highest, mass_mode, certified)
        call require(certified, 'positive current mass spectrum')
        raw_condition = highest/lowest
        print '(a,3i6,3es16.8)', 'current mass kind, n, size, min, max, condition ', &
            kind, n, nj, lowest, highest, raw_condition
        nh = op%fs%nn
        print '(a,2i5,2es14.6)', 'current mass minimum mode section, orbit ', &
            kind, n, norm2(mass_mode(:nh)), norm2(mass_mode(nh + 1:))
        print '(a,2i5,4es14.6)', 'current mass section/orbit diagonal min,max ', &
            kind, n, minval([(current_mass(i, i), i=1, nh)]), &
            maxval([(current_mass(i, i), i=1, nh)]), &
            minval([(current_mass(i, i), i=nh + 1, nj)]), &
            maxval([(current_mass(i, i), i=nh + 1, nj)])
        print '(a,2i5,es14.6)', 'current mass section-orbit cross block ', &
            kind, n, maxval(abs(current_mass(:nh, nh + 1:)))
        ! Diagonal scaling removes the different edge and orbit basis amplitudes.
        mass_diagonal = [(current_mass(i, i), i=1, nj)]
        call require(minval(mass_diagonal) > 0, 'positive current mass diagonal')
        mass_diagonal = sqrt(mass_diagonal)
        do j = 1, nj
            mass_copy(:, j) = current_mass(:, j)/(mass_diagonal*mass_diagonal(j))
        end do
        call spd_extremes(mass_copy, lowest, highest, mass_mode, certified)
        call require(certified, 'positive scaled current mass spectrum')
        mass_condition = highest/lowest
        print '(a,3i6,3es16.8)', &
            'scaled current mass kind, n, size, min, max, condition ', &
            kind, n, nj, lowest, highest, mass_condition
        print '(a,2i5,2es14.6)', 'scaled current mass minimum mode section, orbit ', &
            kind, n, norm2(mass_mode(:nh)), norm2(mass_mode(nh + 1:))
        do t = 1, op%fs%mm%topology%n_triangles
            block = 0
                do q = 1, op%fs%nq
                    call vrmhd_section_basis(op, t, q, x, weight, basis, status)
                    call require(status == 0, 'gauge physical Riesz basis')
                    call vrmhd_section_gauge_basis(op, t, q, gauge_value, &
                        gauge_gradient, status)
                    call require(status == 0, 'multiplier physical Riesz basis')
                    block = block + weight*matmul(transpose(basis), basis)
                        do j = 1, op%fs%nl0
                            col = gauge_index(op%fs%sp(0)%dof(j, t))
                            if (col == 0) cycle
                            do i = 1, op%fs%nl0
                                row = gauge_index(op%fs%sp(0)%dof(i, t))
                                if (row == 0) cycle
                                value = weight*(gauge_value(i)*gauge_value(j) &
                                    + dot_product(gauge_gradient(:, i), &
                                    gauge_gradient(:, j)))
                                trial(row, col) = trial(row, col) + value
                            end do
                        end do
                    end do
                    do j = 1, op%nlocal
                        col = free_index(op%map(j, t))
                        if (col == 0) cycle
                        do i = 1, op%nlocal
                            row = free_index(op%map(i, t))
                            if (row == 0) cycle
                            trial(row, col) = trial(row, col) + block(i, j)
                        end do
                    end do
                end do
                call dpotrf('L', nf, trial, nf, info)
                call require(info == 0, 'trial Riesz Cholesky')
                call dtrsm('L', 'L', 'N', 'N', nf, nf, 1.0_dp, trial, nf, matrix, nf)
                call dtrsm('R', 'L', 'T', 'N', nf, nf, 1.0_dp, trial, nf, matrix, nf)
                ! sigma(B)^2 = eig(B^T B); squaring loses log10(kappa(B)^2) digits,
                ! about two for the measured kappa(B) near 11.
                call dsyrk('L', 'T', nf, nf, 1.0_dp, matrix, nf, 0.0_dp, trial, nf)
                call spd_extremes(trial, lowest, highest, minimum_mode, certified)
                call require(certified, 'reduced stationary Gram spectrum')
                lowest = sqrt(lowest); highest = sqrt(highest)
                print '(a,3i6,2es16.8)', 'reduced singular values kind, n, free, min, max ', &
                    kind, n, nf, lowest, highest
                mode_norms = 0
                do i = 1, steady%ndof
                    row = free_index(i)
                    if (row == 0) cycle
                    if (i <= 6*op%scalar_fs%n0) then
                        group = (i - 1)/op%scalar_fs%n0 + 1
                    else if (i <= 6*op%scalar_fs%n0 + op%fs%sp(1)%ndof) then
                        group = 7
                    else
                        group = 8
                    end if
                    mode_norms(group) = mode_norms(group) + minimum_mode(row)**2
                end do
                print '(a,8es12.4)', 'reduced minimum mode rho,vx,vy,vz,Te,Ti,A,lambda ', &
                    sqrt(mode_norms)
                sigma_min = lowest
                call require(lowest > 1.0e-4_dp, 'reduced stationary inf-sup lower bound')
            end subroutine

            subroutine check_harmonic(kind)
                integer, intent(in) :: kind
                type(vrmhd_section_t) :: op
                type(vrmhd_stationary_t) :: steady
                type(vrmhd_parameters_t) :: p
                real(dp), allocatable :: u(:), load(:), h(:), drive(:), expected(:)
                real(dp), allocatable :: source(:, :, :), js(:, :, :), basis(:, :)
                real(dp) :: x(2), weight, field(3), period, defect
                integer :: status, t, q, i
                call prepare(kind, 2, op, steady, p, u, load)
                allocate (h(op%fs%sp(1)%ndof), drive(op%ndof), expected(op%ndof), &
                    basis(36, op%nlocal), &
                    source(6, op%fs%nq, op%fs%mm%topology%n_triangles), &
                    js(3, op%fs%nq, op%fs%mm%topology%n_triangles))
                source = 0.37_dp; js = 0.23_dp
                call vrmhd_drive_build(op, VRMHD_CURRENT_SOURCE, 2.3_dp, 0.0_dp, js, &
                    h, source, drive, status)
                call require(status == 0, 'current-source convention')
                call require(maxval(abs(source(:3, :, :) - js)) < epsilon(1.0_dp), &
                    'prescribed current retained')
                call require(maxval(abs(drive)) < epsilon(1.0_dp), 'zero voltage load')
                call vrmhd_drive_build(op, VRMHD_LOOP_VOLTAGE, 2.3_dp, 0.41_dp, js, &
                    h, source, drive, status)
                call require(status == 0, 'voltage convention')
                call require(maxval(abs(source(:3, :, :))) < epsilon(1.0_dp), &
                    'voltage convention has zero prescribed current')
                call require(maxval(abs(source(4:, :, :) - 0.37_dp)) < epsilon(1.0_dp), &
                    'non-electrical sources retained')
                defect = maxval(abs(csc_matvec(op%fs%curl%a, h)))
                print '(a,i3,es14.6)', 'harmonic cochain curl ', kind, defect
                call require(defect < 64*epsilon(1.0_dp), 'closed discrete harmonic')
                call require(maxval(abs(h(:op%fs%nn))) < epsilon(1.0_dp), &
                    'zero increments along section edges')
                call require(maxval(abs(h(op%fs%nn + 1:) - 1)) < epsilon(1.0_dp), &
                    'unit increments around periodic orbit edges')
                period = 2.3_dp
                if (kind == GS_ROTATION) period = 2*acos(-1.0_dp)
                expected = 0
                do t = 1, op%fs%mm%topology%n_triangles
                    do q = 1, op%fs%nq
                        call vrmhd_section_basis(op, t, q, x, weight, basis, status)
                        call require(status == 0, 'harmonic load basis')
                        field = [0.0_dp, 0.0_dp, 1/period]
                        if (kind == GS_ROTATION) field = [0.0_dp, 1/(period*x(1)), 0.0_dp]
                        do i = 1, op%nlocal
                            expected(op%map(i, t)) = expected(op%map(i, t)) &
                                + 0.41_dp*weight*dot_product(field, basis(5:7, i))
                        end do
                    end do
                end do
                defect = maxval(abs(drive - expected))
                print '(a,i3,es14.6)', 'independent harmonic voltage load ', kind, defect
                call require(defect < 64*epsilon(1.0_dp), 'Cartesian harmonic load and sign')
                call vrmhd_drive_build(op, VRMHD_CURRENT_SOURCE, 2.3_dp, 0.41_dp, js, &
                    h, source, drive, status)
                call require(status /= 0, 'reject voltage in current-source convention')
            end subroutine

            subroutine require(ok, message)
                logical, intent(in) :: ok
                character(*), intent(in) :: message
                if (ok) return
                print '(a)', message
                error stop 'stationary MHD gate'
            end subroutine

            subroutine exact_at(op, x, exact)
                type(vrmhd_section_t), intent(in) :: op
                real(dp), intent(in) :: x(2)
                real(dp), intent(out) :: exact(48)
                if (op%symmetry%kind == GS_ROTATION) then
                    call mms_rotation(x(1), x(2), op%symmetry%k, exact)
                else
                    call mms_screw(x(1), x(2), op%symmetry%k, exact)
                end if
            end subroutine

            subroutine manufactured_current(x, data, current)
                real(dp), intent(in) :: x(2), data(:)
                real(dp), intent(out) :: current(3)
                real(dp) :: exact(48)
                if (nint(data(1)) == GS_ROTATION) then
                    call mms_rotation(x(1), x(2), data(2), exact)
                else
                    call mms_screw(x(1), x(2), data(2), exact)
                end if
                current = exact(8:10)
            end subroutine

            pure function norms(u) result(e)
                real(dp), intent(in) :: u(36)
                real(dp) :: e(8)
                integer :: i
                e(1) = u(1)**2 + sum(u(13:15)**2)
                do i = 2, 4
                    e(i) = u(i)**2 + u(i + 14)**2 + u(i + 17)**2 + u(i + 20)**2
                end do
                e(5) = sum(u(5:7)**2) + sum(u(25:27)**2)
                ! Weak Ampere controls J in L2; no residual term contains curl J.
                e(6) = sum(u(8:10)**2)
                e(7) = u(11)**2 + sum(u(28:30)**2)
                e(8) = u(12)**2 + sum(u(31:33)**2)
            end function

            subroutine prepare(kind, n, op, steady, p, u, load)
                integer, intent(in) :: kind, n
                type(vrmhd_section_t), intent(out) :: op
                type(vrmhd_stationary_t), intent(out) :: steady
                type(vrmhd_parameters_t), intent(out) :: p
                real(dp), allocatable, intent(out) :: u(:), load(:)
                type(mesh_2d_t) :: mesh
                type(meridian_mesh_t) :: mm
                type(gs_symmetry_t) :: symmetry
                real(dp) :: shift
                real(dp), allocatable :: wall(:)
                integer :: status, current_start
                symmetry = gs_symmetry_t(kind, GS_CHART_CARTESIAN, 0.0_dp)
                shift = 0
                if (kind == GS_ROTATION) then
                    shift = 2; symmetry%chart = GS_CHART_RZ
                end if
                if (kind == GS_SCREW) symmetry%k = 0.7_dp
                call mesh%create_rectangular(n + 1, n + 1, shift + 0.1_dp, shift + 1.1_dp, &
                    0.2_dp, 1.2_dp)
                call meridian_mesh_straight(mesh, 1, mm, status)
                call require(status == 0, 'mesh')
                call vrmhd_section_build(mm, symmetry, 1, 1, op, status)
                call require(status == 0, 'section')
                call vrmhd_stationary_build(op, steady, status)
                call require(status == 0, 'Coulomb gauge')
                p%mi = 1; p%mu0 = 1; p%eta0 = 0.2_dp; p%nu0 = 0.2_dp
                p%chi_parallel0 = 0.2_dp; p%chi_perpendicular = 0.05_dp
                p%d_parallel = 0.1_dp; p%d_perpendicular = 0.1_dp; p%epsilon_b = 0.1_dp
                allocate (u(steady%ndof), load(op%ndof))
                u = 0
                call project(op, u(:op%ndof), load)
                allocate (wall(op%fs%sp(1)%ndof))
                call current_wall_dofs(op, manufactured_current, &
                    [real(kind, dp), symmetry%k], wall)
                current_start = 6*op%scalar_fs%n0 + op%fs%sp(1)%ndof
                where (op%fs%sp(1)%bnd) u(current_start + 1:op%ndof) = wall
            end subroutine

            subroutine project(op, u, load)
                type(vrmhd_section_t), intent(in) :: op
                real(dp), intent(out) :: u(op%ndof), load(op%ndof)
                type(csc_t) :: mass
                type(fortsparse_status_t) :: st
                integer, allocatable :: rows(:), cols(:)
                real(dp), allocatable :: values(:), rhs(:)
                real(dp) :: basis(36, op%nlocal), block(op%nlocal, op%nlocal)
                real(dp) :: exact(48), x(2), weight
                integer :: i, j, t, q, pos, status, nt
                nt = op%fs%mm%topology%n_triangles
                allocate (rows(nt*op%nlocal**2), cols(nt*op%nlocal**2), &
                    values(nt*op%nlocal**2), rhs(op%ndof))
                rhs = 0; load = 0
                do t = 1, nt
                    block = 0
                        do q = 1, op%fs%nq
                            call vrmhd_section_basis(op, t, q, x, weight, basis, status)
                            call require(status == 0, 'projection basis')
                            call exact_at(op, x, exact)
                            block = block + weight*matmul(transpose(basis), basis)
                                do i = 1, op%nlocal
                                    j = op%map(i, t)
                                    rhs(j) = rhs(j) + weight*dot_product(basis(:, i), exact(:36))
                                    load(j) = load(j) + weight*dot_product(basis(:12, i), exact(37:48))
                                end do
                            end do
                            pos = (t - 1)*op%nlocal**2
                            do j = 1, op%nlocal
                                do i = 1, op%nlocal
                                    pos = pos + 1
                                    rows(pos) = op%map(i, t); cols(pos) = op%map(j, t)
                                    values(pos) = block(i, j)
                                end do
                            end do
                        end do
                        call csc_from_triplet(op%ndof, op%ndof, rows, cols, values, mass, st)
                        call require(status_ok(st), 'projection assembly')
                        call sparse_solve_once(mass, rhs, u, st)
                        call require(status_ok(st), 'projection solve')
                    end subroutine

                    subroutine check_boundary(kind)
                        integer, intent(in) :: kind
                        type(vrmhd_section_t) :: op
                        type(vrmhd_stationary_t) :: steady
                        type(vrmhd_parameters_t) :: p
                        real(dp), allocatable :: u(:), load(:), wall(:), volume(:)
                        real(dp), allocatable :: basis(:, :)
                        real(dp) :: exact(48), x(2), weight, value
                        real(dp) :: defect
                        integer :: t, q, i, row, status
                        call prepare(kind, 2, op, steady, p, u, load)
                        allocate (wall(op%fs%sp(1)%ndof), volume(op%ndof), basis(36, op%nlocal))
                        call require(all(op%fixed(6*op%scalar_fs%n0 + op%fs%sp(1)%ndof + 1:) &
                            .eqv. op%fs%sp(1)%bnd), 'fixed current trace mask')
                        call vrmhd_current_lift(op, &
                            u(6*op%scalar_fs%n0 + op%fs%sp(1)%ndof + 1:op%ndof), wall)
                        call require(maxval(abs(wall), mask=op%fs%sp(1)%bnd) &
                            > 1.0e-3_dp, 'nonzero manufactured current trace')
                        call require(maxval(abs(wall), mask=.not. op%fs%sp(1)%bnd) &
                            == 0, 'zero interior current lift')
                        volume = 0
                        do t = 1, op%fs%mm%topology%n_triangles
                            do q = 1, op%fs%nq
                                call vrmhd_section_basis(op, t, q, x, weight, basis, status)
                                call exact_at(op, x, exact)
                                do i = 1, op%nlocal
                                    if (op%row(i) /= 4) cycle
                                    value = dot_product(exact(8:10), basis(8:10, i)) &
                                        - dot_product(exact(25:27), basis(34:36, i))
                                    row = op%map(i, t)
                                    volume(row) = volume(row) + weight*value
                                end do
                            end do
                        end do
                        defect = maxval(abs(volume), mask=.not. op%fixed)
                        print '(a,i2,2es14.6)', 'Ampere H0(curl) defect / wall current ', &
                            kind, defect, maxval(abs(wall))
                        call require(defect < 2.0e-8_dp, 'independent Ampere integration by parts')
                    end subroutine

                    subroutine check_derivatives(kind)
                        integer, intent(in) :: kind
                        type(vrmhd_section_t) :: op
                        type(vrmhd_stationary_t) :: steady
                        type(vrmhd_parameters_t) :: p
                        type(csc_t) :: jac
                        real(dp), allocatable :: u(:), load(:), source(:, :, :), du(:), zero(:)
                        real(dp), allocatable :: r(:), jv(:), plus(:), minus(:), unused(:), boundary(:)
                        real(dp) :: eps, error, previous
                        integer :: i, status
                        call prepare(kind, 2, op, steady, p, u, load)
                        allocate (source(6, op%fs%nq, op%fs%mm%topology%n_triangles))
                        allocate (du(steady%ndof), zero(steady%ndof), r(steady%ndof), &
                            jv(steady%ndof), plus(steady%ndof), minus(steady%ndof), &
                            unused(steady%ndof), boundary(op%ndof))
                        source = 0; zero = 0; boundary = u(:op%ndof)
                        do i = 1, steady%ndof
                            du(i) = sin(real(i, dp))/10
                        end do
                        u(op%ndof + 1:) = 0.03_dp
                        call vrmhd_stationary_residual(op, steady, p, u, du, boundary, source, &
                            load, r, jv, status)
                        call require(status == 0, 'gauged residual')
                        call vrmhd_stationary_jacobian(op, steady, p, u, boundary, source, jac, status)
                        call require(status == 0, 'gauged sparse Jacobian')
                        error = maxval(abs(csc_matvec(jac, du) - jv))
                        call require(error < 1.0e-11_dp, 'sparse and matrix-free gauged Jacobian')
                        previous = huge(1.0_dp)
                        do i = 1, 3
                            eps = 0.02_dp/2**i
                            call vrmhd_stationary_residual(op, steady, p, u + eps*du, zero, &
                                boundary, source, load, plus, unused, status)
                            call vrmhd_stationary_residual(op, steady, p, u - eps*du, zero, &
                                boundary, source, load, minus, unused, status)
                            error = maxval(abs((plus - minus)/(2*eps) - jv))
                            print '(a,2i3,es14.6)', 'gauged centered difference ', kind, i, error
                            call require(error < 0.3_dp*previous, 'quadratic centered differences')
                            previous = error
                        end do
                    end subroutine

                    subroutine solve_manufactured(kind, n, error, curl_error)
                        integer, intent(in) :: kind, n
                        real(dp), intent(out) :: error(8), curl_error(2)
                        type(vrmhd_section_t) :: op
                        type(vrmhd_stationary_t) :: steady
                        type(vrmhd_parameters_t) :: p
                        type(csc_t) :: jac
                        type(fortsparse_status_t) :: st
                        real(dp), allocatable :: u(:), load(:), source(:, :, :), boundary(:)
                        real(dp), allocatable :: r(:), unused(:), zero(:), step(:)
                        real(dp), allocatable :: basis(:, :), projected(:)
                        real(dp) :: history(10), x(2), weight, exact(48), field(36), order
                        real(dp) :: current_error(2), projection_error(2), projection(36)
                        real(dp) :: lower(2), upper(2), vertices(2, 3), margin, curl_squared
                        integer :: i, it, status, t, q, interior_count
                        logical :: quadratic, interior
                        call prepare(kind, n, op, steady, p, u, load)
                        allocate (source(6, op%fs%nq, op%fs%mm%topology%n_triangles), &
                            boundary(op%ndof), basis(36, op%nlocal), &
                            projected(op%ndof))
                        allocate (r(steady%ndof), unused(steady%ndof), zero(steady%ndof), &
                            step(steady%ndof))
                        source = 0; zero = 0; boundary = u(:op%ndof)
                        projected = u(:op%ndof)

                        do i = 1, op%ndof
                            if (op%fixed(i)) cycle
                            u(i) = u(i) + 0.02_dp*sin(real(i, dp))
                        end do
                        quadratic = .false.
                        do it = 1, size(history)
                            call vrmhd_stationary_residual(op, steady, p, u, zero, boundary, &
                                source, load, r, unused, status)
                            call require(status == 0, 'Newton residual')
                            history(it) = sqrt(sum(r*r))
                            print '(a,3i4,es16.8)', 'Newton residual ', kind, n, it - 1, history(it)
                            if (it > 2) then
                                if (history(it) > 1.0e-12_dp .and. &
                                    history(it - 1) < 0.2_dp*history(it - 2)) then
                                    order = log(history(it)/history(it - 1)) &
                                        /log(history(it - 1)/history(it - 2))
                                    print '(a,f9.4)', 'Newton observed order ', order
                                    quadratic = quadratic .or. order > 1.8_dp
                                end if
                            end if
                            if (history(it) < 1.0e-10_dp) exit
                            call vrmhd_stationary_jacobian(op, steady, p, u, boundary, source, &
                                jac, status)
                            call require(status == 0, 'Newton assembled sparse Jacobian')
                            call sparse_solve_once(jac, -r, step, st)
                            if (.not. status_ok(st)) print '(a,i0)', 'direct solve status ', st%code
                            call require(status_ok(st), 'Newton direct solve')
                            u = u + step
                        end do
                        call require(it <= size(history), 'manufactured Newton convergence')
                        call require(quadratic, 'quadratic Newton convergence')
                        print '(a,es14.6)', 'Coulomb multiplier maximum ', maxval(abs(u(op%ndof + 1:)))
                        error = 0; current_error = 0; projection_error = 0; curl_error = 0
                        lower = minval(op%fs%mm%topology%vertices, dim=2)
                        upper = maxval(op%fs%mm%topology%vertices, dim=2)
                        margin = 2.0_dp/n - 64*epsilon(1.0_dp)
                        interior_count = 0
                        do t = 1, op%fs%mm%topology%n_triangles
                            vertices = op%fs%mm%topology%vertices(:, &
                                op%fs%mm%topology%triangles(:, t))
                            interior = all(minval(vertices, dim=2) >= lower + margin) .and. &
                                all(maxval(vertices, dim=2) <= upper - margin)
                            if (interior) interior_count = interior_count + 1
                            do q = 1, op%fs%nq
                                call vrmhd_section_basis(op, t, q, x, weight, basis, status)
                                call exact_at(op, x, exact)
                                field = matmul(basis, u(op%map(:, t)))
                                projection = matmul(basis, projected(op%map(:, t)))
                                error = error + weight*norms(field - exact(:36))
                                current_error = current_error + weight*[ &
                                    sum((field(8:10) - exact(8:10))**2), &
                                    sum((field(34:36) - exact(34:36))**2)]
                                curl_squared = weight*sum((field(34:36) - exact(34:36))**2)
                                curl_error(1) = curl_error(1) + curl_squared
                                if (interior) curl_error(2) = curl_error(2) + curl_squared
                                projection_error = projection_error + weight*[ &
                                    sum((projection(8:10) - exact(8:10))**2), &
                                    sum((projection(34:36) - exact(34:36))**2)]
                            end do
                        end do
                        error = sqrt(error)
                        curl_error = sqrt(curl_error)
                        print '(a,2es14.6)', 'solved J L2, curl L2 errors ', sqrt(current_error)
                        print '(a,2es14.6)', 'projected J L2, curl L2 errors ', sqrt(projection_error)
                        print '(a,3i6,2es14.6)', 'diagnostic curl J kind, n, interior cells, errors ', &
                            kind, n, interior_count, curl_error
                        if (interior_count == 0) &
                            print '(a)', 'diagnostic interior curl J unavailable: no interior cells'
                    end subroutine
                end program test_vrmhd_stationary
