module hessian_test_case
    ! Circular cross-section on a large-aspect-ratio torus with a smooth
    ! frozen state; the Hessian assembly needs frozen data, not an equilibrium.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_ideal_map, only: ideal_frozen_t
    implicit none
    private
    public :: hessian_field_t, boundary_point, r0, a
    real(dp), parameter :: r0 = 6.2_dp, a = 1.6_dp
    type, extends(ideal_frozen_t) :: hessian_field_t
    contains
        procedure :: evaluate => field_evaluate
    end type
contains
    subroutine boundary_point(t, point)
        real(dp), intent(in) :: t
        real(dp), intent(out) :: point(2)
        point = [r0 + a*cos(8*atan(1.0_dp)*t), a*sin(8*atan(1.0_dp)*t)]
    end subroutine

    subroutine field_evaluate(field, point, chi, chi_gradient, b3, pressure)
        class(hessian_field_t), intent(in) :: field
        real(dp), intent(in) :: point(2)
        real(dp), intent(out) :: chi, chi_gradient(2), b3, pressure
        chi = ((point(1) - r0)**2 + point(2)**2)/a**2
        chi_gradient = 2*[point(1) - r0, point(2)]/a**2
        b3 = 1.5_dp
        pressure = 0.2_dp*(1 - chi)
    end subroutine
end module hessian_test_case

program test_ideal_map_hessian
    ! The sum-factorised coefficient Hessian of kin6d_ideal_map against the
    ! frozen point-loop assembly (test/support/ideal_map_hessian_oracle.f90) on
    ! random 3D states: full Fourier coupling (ntor > 0), P2 through P6, the plain
    ! transversal basis with its axis expansion and gauge form, the nested
    ! axis basis (vector axis modes) and the position gauge. Both assemble into
    ! the same CSR pattern; the matrices must agree to roundoff relative to the
    ! largest entry, and the production matrix must be exactly symmetric.
    ! The Hessian, energy and gradient must not depend on the OpenMP team
    ! size (bitwise), including teams larger than the number of cells.
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use kin6d_grad_shafranov, only: gs_state_t, initialize_gs
    use kin6d_equilibrium_chart, only: cylindrical_chart_t
    use kin6d_ideal_map, only: ideal_map_t, ideal_map_init, ideal_map_freeze, &
        ideal_map_hessian, ideal_map_energy, ideal_map_expand, ideal_map_reduce, &
        IDEAL_MAP_OK
    use fortfem_element_to_csr, only: csr_matrix_t
    use hessian_test_case, only: hessian_field_t, boundary_point, a
    use ideal_map_hessian_oracle, only: ideal_map_hessian_reference
    !$ use omp_lib, only: omp_get_max_threads, omp_set_num_threads
    implicit none
    integer :: degree, ntor, variant, failures, first_degree, last_degree, ntor_step
    character(16) :: argument
    failures = 0
    first_degree = 2; last_degree = 3; ntor_step = 1
    if (command_argument_count() > 0) then
        call get_command_argument(1, argument)
        read (argument, *) first_degree
        if (first_degree < 4 .or. first_degree > 6) error stop 'expected degree 4..6'
        last_degree = first_degree; ntor_step = 2
    end if
    do degree = first_degree, last_degree
        do ntor = 0, 2, ntor_step
            do variant = 1, 3
                call check(degree, ntor, variant, failures)
            end do
        end do
    end do
    if (failures > 0) error stop 'ideal-map Hessian qualification failed'
contains
    subroutine check(degree, ntor, variant, failures)
        ! variant 1: transversal basis, axis monomials, gauge form;
        ! 2: nested axis basis (omega eliminated, vector axis modes);
        ! 3: position gauge.
        integer, intent(in) :: degree, ntor, variant
        integer, intent(inout) :: failures
        type(gs_state_t) :: state
        type(ideal_map_t) :: map
        type(cylindrical_chart_t) :: chart
        type(csr_matrix_t) :: fast, oracle
        type(hessian_field_t) :: field
        real(dp), allocatable :: contour(:, :), x(:, :, :)
        real(dp) :: scale, diff, energy, asym, amplitude
        integer :: j, status, node, mode, comp, n
        integer(int64) :: seed
        n = 16
        if (degree > 3) n = 8
        allocate (contour(2, n))
        do j = 1, n
            call boundary_point(real(j - 1, dp)/n, contour(:, j))
        end do
        call initialize_gs(state, contour, (8*atan(1.0_dp)*a/n)**2/2, status, degree, &
            curve=boundary_point, geometry_only=degree > 3)
        if (status /= 0) error stop 'mesh'
        call ideal_map_init(map, state, chart, 3, ntor, status, nphi=3*ntor + 2)
        if (status /= IDEAL_MAP_OK) error stop 'map init'
        map%omega_elimination = variant == 2
        map%position_gauge = variant == 3
        call ideal_map_freeze(map, field, status)
        if (status /= IDEAL_MAP_OK) error stop 'freeze'
        ! Deterministic pseudo-random state: every mode and component moves, the
        ! step stays far below the element size so that J > 0 everywhere.
        x = map%x
        amplitude = 4.0e-3_dp*a
        if (degree > 3) amplitude = amplitude*(3.0_dp/degree)**2
        seed = 12345_int64 + 100*degree + 10*ntor + variant
        do node = 1, map%nodes
            if (map%boundary(node)) cycle
            do mode = 1, map%nmode
                do comp = 1, 3
                    seed = mod(seed*1103515245_int64 + 12345_int64, 2147483647_int64)
                    x(comp, mode, node) = x(comp, mode, node) &
                        + amplitude*(real(seed, dp)/2147483647.0_dp - 0.5_dp)
                end do
            end do
        end do
        call ideal_map_energy(map, x, energy, status)
        if (status /= IDEAL_MAP_OK) error stop 'perturbed state is folded'
        call ideal_map_hessian(map, x, fast)
        call ideal_map_hessian_reference(map, x, oracle)
        if (fast%n /= oracle%n) error stop 'matrix dimension mismatch'
        if (size(fast%val) /= size(oracle%val)) error stop 'pattern mismatch'
        if (size(fast%rowptr) /= size(oracle%rowptr)) error stop 'row count mismatch'
        if (any(fast%rowptr /= oracle%rowptr)) error stop 'row pointer mismatch'
        if (any(fast%colidx /= oracle%colidx)) error stop 'column index mismatch'
        scale = maxval(abs(oracle%val))
        diff = maxval(abs(fast%val - oracle%val))
        asym = asymmetry(fast)
        print '(a,3i3,3es12.4)', 'degree ntor variant max|new-old|/max|old|, asym, max|H| ', &
            degree, ntor, variant, diff/scale, asym, scale
        if (.not. diff <= 1.0e-12_dp*scale) failures = failures + 1
        if (asym /= 0) failures = failures + 1
        call thread_invariance(map, x, fast, failures)
        if (degree > 3) then
            print '(a,5i7)', 'nodes cells nq full-local reduced-local ', &
                map%nodes, map%cells, map%nq, 3*map%nmode*map%nloc, &
                2*map%nmode*(map%nloc + map%npoly)
            call block_actions(map, fast, oracle, failures)
            if (ntor == 2) call directional_check(map, x, fast, failures)
        end if
    end subroutine

    subroutine thread_invariance(map, x, reference, failures)
        ! Hessian, energy and gradient with 1, 3 and 7 threads against the
        ! default team: bitwise identical (the work items keep the per-cell
        ! operation order).
        type(ideal_map_t), intent(in) :: map
        real(dp), intent(in) :: x(:, :, :)
        type(csr_matrix_t), intent(in) :: reference
        integer, intent(inout) :: failures
        type(csr_matrix_t) :: other
        real(dp), allocatable :: g0(:), g1(:)
        real(dp) :: e0, e1
        integer :: status, k, saved
        integer, parameter :: teams(3) = [1, 3, 7]
        saved = 1
        !$ saved = omp_get_max_threads()
        allocate (g0(map%ndof), g1(map%ndof))
        call ideal_map_energy(map, x, e0, status, g0)
        do k = 1, size(teams)
            !$ call omp_set_num_threads(teams(k))
            call ideal_map_hessian(map, x, other)
            call ideal_map_energy(map, x, e1, status, g1)
            if (any(other%val /= reference%val) .or. e1 /= e0 .or. any(g1 /= g0)) then
                print '(a,i0,a,i0)', 'thread dependence: team ', teams(k), ' cells ', &
                    map%cells
                failures = failures + 1
            end if
        end do
        !$ call omp_set_num_threads(saved)
    end subroutine

    subroutine multiply(matrix, v, product)
        type(csr_matrix_t), intent(in) :: matrix
        real(dp), intent(in) :: v(:)
        real(dp), intent(out) :: product(:)
        integer :: i
        do i = 1, matrix%n
            product(i) = dot_product(matrix%val(matrix%rowptr(i):matrix%rowptr(i + 1) - 1), &
                v(matrix%colidx(matrix%rowptr(i):matrix%rowptr(i + 1) - 1)))
        end do
    end subroutine

    subroutine block_actions(map, fast, oracle, failures)
        ! Separate nodal-1, nodal-2 and axis inputs expose all nine mixed-block
        ! actions. Each output block has a relative gate plus a roundoff-scale
        ! absolute allowance for blocks whose exact action is zero.
        type(ideal_map_t), intent(in) :: map
        type(csr_matrix_t), intent(in) :: fast, oracle
        integer, intent(inout) :: failures
        real(dp) :: v(map%nred), actual(map%nred), expected(map%nred), scale, worst, budget
        integer :: i, input_block, output_block, nodal_end, block(map%nred)
        nodal_end = 2*map%nmode*map%nodes
        block = 3
        do i = 1, nodal_end
            block(i) = 1 + mod(i - 1, 2)
        end do
        worst = 0
        do input_block = 1, 3
            do i = 1, map%nred
                v(i) = sin(0.37_dp*i)
            end do
            where (block /= input_block) v = 0
            call multiply(fast, v, actual)
            call multiply(oracle, v, expected)
            scale = max(norm2(expected), tiny(1.0_dp))
            do output_block = 1, 3
                budget = 1.0e-12_dp*norm2(pack(expected, block == output_block)) &
                    + 128*epsilon(1.0_dp)*scale
                worst = max(worst, &
                    norm2(pack(actual - expected, block == output_block))/budget)
            end do
        end do
        print '(a,es14.5)', 'mixed-block action discrepancy / tolerance ', worst
        if (.not. worst <= 1.0_dp) failures = failures + 1
    end subroutine

    subroutine directional_check(map, x, hessian, failures)
        ! Differentiate the separately generated gradient with the frozen
        ! reduction P held fixed. This checks H v = d(P^T g)[P v], rather than
        ! symmetry, which both assembled implementations impose by reflection.
        type(ideal_map_t), intent(in) :: map
        real(dp), intent(in) :: x(:, :, :)
        type(csr_matrix_t), intent(in) :: hessian
        integer, intent(inout) :: failures
        real(dp) :: vr(map%nred), v(map%ndof), hv(map%nred), gp(map%ndof), gm(map%ndof)
        real(dp) :: grp(map%nred), grm(map%nred), step, errors(3), energy, scale
        real(dp) :: xp(size(x, 1), size(x, 2), size(x, 3))
        integer :: i, k, status
        do i = 1, map%nred
            vr(i) = cos(0.37_dp*i)
        end do
        where (map%fixed) vr = 0
        call ideal_map_expand(map, vr, v)
        call multiply(hessian, vr, hv)
        scale = max(norm2(hv), tiny(1.0_dp))
        do k = 1, size(errors)
            step = 1.0e-5_dp/2**(k - 1)
            xp = x + step*reshape(v, shape(x))
            call ideal_map_energy(map, xp, energy, status, gp)
            if (status /= IDEAL_MAP_OK) error stop 'positive directional state folded'
            xp = x - step*reshape(v, shape(x))
            call ideal_map_energy(map, xp, energy, status, gm)
            if (status /= IDEAL_MAP_OK) error stop 'negative directional state folded'
            call ideal_map_reduce(map, gp, grp)
            call ideal_map_reduce(map, gm, grm)
            errors(k) = norm2((grp - grm)/(2*step) - hv)/scale
        end do
        print '(a,3es14.5)', 'gradient directional relative errors h=1e-5,5e-6,2.5e-6 ', errors
        ! Same relative threshold as ideal_map_derivatives_and_rotation_null.
        if (.not. errors(size(errors)) <= 1.0e-6_dp) failures = failures + 1
    end subroutine

    real(dp) function asymmetry(matrix) result(worst)
        ! Largest |H(i,j) - H(j,i)| over the stored entries. FortFEM's CSR
        ! contract sorts columns per row; binary search keeps P6 affordable.
        type(csr_matrix_t), intent(in) :: matrix
        integer :: i, k, l, first, last
        worst = 0
        do i = 1, matrix%n
            do k = matrix%rowptr(i), matrix%rowptr(i + 1) - 1
                first = matrix%rowptr(matrix%colidx(k))
                last = matrix%rowptr(matrix%colidx(k) + 1) - 1
                do while (first <= last)
                    l = first + (last - first)/2
                    if (matrix%colidx(l) < i) then
                        first = l + 1
                    else if (matrix%colidx(l) > i) then
                        last = l - 1
                    else
                        exit
                    end if
                end do
                if (first > last) error stop 'asymmetric sparsity pattern'
                worst = max(worst, abs(matrix%val(k) - matrix%val(l)))
            end do
        end do
    end function
end program test_ideal_map_hessian
