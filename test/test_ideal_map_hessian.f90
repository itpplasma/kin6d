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
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use kin6d_grad_shafranov, only: gs_state_t, initialize_gs
    use kin6d_equilibrium_chart, only: cylindrical_chart_t
    use kin6d_ideal_map, only: ideal_map_t, ideal_map_init, ideal_map_freeze, &
        ideal_map_hessian, ideal_map_energy, IDEAL_MAP_OK
    use fortfem_element_to_csr, only: csr_matrix_t
    use hessian_test_case, only: hessian_field_t, boundary_point, a
    use ideal_map_hessian_oracle, only: ideal_map_hessian_reference
    implicit none
    integer :: degree, ntor, variant, failures
    failures = 0
    do degree = 2, 6
        do ntor = 0, 2
            do variant = 1, 3
                call check(degree, ntor, variant, failures)
            end do
        end do
    end do
    if (failures > 0) error stop 'sum-factorised Hessian differs from the oracle'
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
        integer, parameter :: n = 16
        integer :: j, status, node, mode, comp
        integer(int64) :: seed
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
        if (size(fast%val) /= size(oracle%val)) error stop 'pattern mismatch'
        scale = maxval(abs(oracle%val))
        diff = maxval(abs(fast%val - oracle%val))
        asym = asymmetry(fast)
        print '(a,3i3,3es12.4)', 'degree ntor variant max|new-old|/max|old|, asym, max|H| ', &
            degree, ntor, variant, diff/scale, asym, scale
        if (.not. diff <= 1.0e-12_dp*scale) failures = failures + 1
        if (asym /= 0) failures = failures + 1
    end subroutine

    real(dp) function asymmetry(matrix) result(worst)
        ! Largest |H(i,j) - H(j,i)| over the stored entries.
        type(csr_matrix_t), intent(in) :: matrix
        integer :: i, k, l
        worst = 0
        do i = 1, matrix%n
            do k = matrix%rowptr(i), matrix%rowptr(i + 1) - 1
                do l = matrix%rowptr(matrix%colidx(k)), &
                    matrix%rowptr(matrix%colidx(k) + 1) - 1
                    if (matrix%colidx(l) /= i) cycle
                    worst = max(worst, abs(matrix%val(k) - matrix%val(l)))
                end do
            end do
        end do
    end function
end program test_ideal_map_hessian
