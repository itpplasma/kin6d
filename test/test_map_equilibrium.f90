module map_equilibrium_oracle
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_ideal_map, only: ideal_frozen_t
    implicit none
    real(dp), parameter :: r0 = 6.2_dp, minor = 0.62_dp, ftor = 32.86_dp
    real(dp), parameter :: c = ftor/(12*r0**3), shear = 0.08_dp
    type, extends(ideal_frozen_t) :: exact_frozen_t
    contains
        procedure :: evaluate => exact_pullback
    end type
contains
    subroutine exact_pullback(field, point, chi, chi_gradient, b3, pressure)
        class(exact_frozen_t), intent(in) :: field
        real(dp), intent(in) :: point(2)
        real(dp), intent(out) :: chi, chi_gradient(2), b3, pressure
        real(dp) :: x(2), grad(2)
        x = [point(1) + shear*point(2)**2, point(2)]
        chi = c*((x(1)**2 - r0*r0)**2 + 4*x(1)**2*x(2)**2 &
                 - 4*r0*r0*minor*minor)
        grad = c*[4*x(1)*(x(1)**2 - r0*r0) + 8*x(1)*x(2)**2, &
                  8*x(1)**2*x(2)]
        chi_gradient = [grad(1), grad(2) + 2*shear*point(2)*grad(1)]
        b3 = ftor/x(1)
        pressure = -16*c*chi
    end subroutine

    subroutine boundary(t, point)
        real(dp), intent(in) :: t
        real(dp), intent(out) :: point(2)
        real(dp) :: theta, r, z
        theta = 2*acos(-1.0_dp)*t
        r = sqrt(r0*r0 + 2*r0*minor*cos(theta))
        z = r0*minor*sin(theta)/r
        point = [r - shear*z*z, z]
    end subroutine
end module

program test_map_equilibrium
    use map_equilibrium_oracle
    use kin6d_grad_shafranov, only: gs_state_t, initialize_gs
    use kin6d_gs_equilibrium, only: gs_geometry_view
    use kin6d_equilibrium_chart, only: cylindrical_chart_t
    use kin6d_ideal_map, only: ideal_map_t, ideal_map_init, ideal_map_freeze
    use kin6d_map_equilibrium, only: map_equilibrium
    use kin6d_equilibrium_state, only: equilibrium_state_t, EQ_VALUES, &
                 EQ_BR, EQ_BZ, EQ_BPHI, EQ_P, EQ_PSI, EQ_PSI_R, EQ_PSI_Z, EQ_UNAVAILABLE
    use kin6d_equilibrium_residual, only: equilibrium_residual_t, equilibrium_residual
    use kin6d_ideal_map_leaves, only: ideal_axisymmetric_field, &
                              ideal_axisymmetric_field_jvp, ideal_axisymmetric_field_vjp
    use omp_lib, only: omp_set_num_threads, omp_get_wtime
    implicit none
    type(gs_state_t) :: mesh, physical
    type(ideal_map_t) :: map
    type(exact_frozen_t) :: frozen
    type(cylindrical_chart_t) :: chart
    type(equilibrium_state_t) :: eq, rejected
    type(equilibrium_residual_t) :: report
    real(dp), allocatable :: contour(:, :)
    real(dp) :: errors(3), rate
    real(dp) :: point(2), gradient(2), flux, pressure, bt, h
    real(dp) :: points(3, 10000), values(EQ_VALUES, 10000), serial(EQ_VALUES, 10000)
    real(dp) :: start, elapsed, r, z, expected(EQ_VALUES), maximum
    integer :: st, i, j, degree, statuses(10000), threads, level, n, k
    character(16) :: argument
    call check_leaf()
    call get_command_argument(1, argument)
    read (argument, *) degree
    do level = 1, 3
        n = 12*2**level
        allocate (contour(2, n))
        do i = 1, n
            call boundary(real(i - 1, dp)/n, contour(:, i))
        end do
        h = 2*acos(-1.0_dp)*minor/n
        call initialize_gs(mesh, contour, h*h/2, st, degree, curve=boundary)
        if (st /= 0) error stop 'reference mesh'
        call ideal_map_init(map, mesh, chart, 1, 0, st)
        if (st /= 0) error stop 'map init'
        call ideal_map_freeze(map, frozen, st)
        if (st /= 0) error stop 'map freeze'
        ! Quadratic deformation composed with curved reference cells: expect
        ! refinement order K, rather than exactness of the nodal interpolant.
        map%x(1, 1, :) = mesh%points(1, :) + shear*mesh%points(2, :)**2
        call map_equilibrium(map, frozen, mesh, eq, st, ftor)
        if (st /= 0) error stop 'map adapter'
        if (eq%lane /= 'map' .or. eq%claim /= 'candidate') error stop 'provenance'
        call gs_geometry_view(mesh, 0.0_dp, physical)
        physical%points = map%x(:2, 1, :)
        physical%curved_cell = .true.; physical%mapped = .true.
        call equilibrium_residual(eq, physical, report, st)
        if (st /= 0) error stop 'map D3'
        print '(a,i0,2(a,es12.4))', 'map P', degree, ' D3=', report%normalized, &
            ' sigma=', report%sigma
        errors(level) = report%normalized
        deallocate (contour)
    end do
    rate = log(errors(2)/errors(3))/log(2.0_dp)
    print '(a,2f8.3)', 'D3 orders=', log(errors(1)/errors(2))/log(2.0_dp), rate
    if (rate < degree) error stop 'map D3 order below K'
    if (errors(3) > 0.01_dp) error stop 'map D3 exceeds 1 percent'
    if (report%sigma <= 0) error stop 'sigma normalization'
    maximum = 0
    do i = 1, size(points, 2)
        r = r0 + 0.3_dp*sin(real(i, dp))
        z = 0.2_dp*cos(real(i, dp)*1.3_dp)
        points(:, i) = [r, z, real(i, dp)*0.1_dp]
    end do
    do k = 1, 4
        threads = 2**(k - 1)
        call omp_set_num_threads(threads)
        start = omp_get_wtime()
        do j = 1, 5
            call eq%field%sample(points, values, statuses)
        end do
        elapsed = (omp_get_wtime() - start)/5
        if (any(statuses /= 0)) error stop 'map sample status'
        if (threads == 1) serial = values
        if (maxval(abs(values - serial)) > 1.0e-12_dp) error stop 'thread parity'
        print '(a,i0,a,es12.4)', 'sample threads=', threads, ' seconds=', elapsed
    end do
    do i = 1, size(points, 2)
        r = points(1, i); z = points(2, i)
        point = [r - shear*z*z, z]
        call frozen%evaluate(point, flux, gradient, bt, pressure)
        gradient = c*[4*r*(r*r - r0*r0) + 8*r*z*z, 8*r*r*z]
        expected = [-gradient(2)/r, gradient(1)/r, ftor/r, pressure, flux, &
                    gradient(1), gradient(2)]
        maximum = max(maximum, maxval(abs(values(:, i) - expected)))
    end do
    print '(a,es12.4)', 'pointwise oracle error=', maximum
    if (maximum > 1.0e-10_dp) error stop 'Piola oracle'
    ! A nonconstant omega changes the toroidal field, not the flux label.
    map%x(3, 1, :) = 0.02_dp*mesh%points(2, :)
    call map_equilibrium(map, frozen, mesh, eq, st, ftor)
    if (st /= 0) error stop 'omega adapter'
    call eq%field%sample(points, values, statuses)
    if (any(statuses /= 0)) error stop 'omega sample'
    maximum = maxval(abs(values(EQ_BPHI, :) - serial(EQ_BPHI, :) &
                         - 0.02_dp*points(1, :)*serial(EQ_BZ, :)))
    if (maximum > 1.0e-11_dp) error stop 'omega Piola term'
    points(:, 1) = [100.0_dp, 100.0_dp, 0.0_dp]
    call eq%field%sample(points(:, :1), values(:, :1), statuses(:1))
    if (statuses(1) == 0) error stop 'outside domain'
    map%ntor = 1
    call map_equilibrium(map, frozen, mesh, rejected, st, ftor)
    if (st /= EQ_UNAVAILABLE) error stop 'unsupported map admitted'
    if (allocated(rejected%field)) error stop 'failed adapter retained a field'
 print '(a)', 'PASS: map state, Piola, D3, thread parity; candidate, bounds unavailable'

contains
    subroutine check_leaf()
        real(dp) :: v(14), d(14), w(6), plus(6), minus(6), t(6), g(14), h
        integer :: j
        v = [1.1_dp, 0.04_dp, 0.1_dp, 0.95_dp, 1.0_dp, 0.2_dp, &
             -0.1_dp, 0.9_dp, 0.03_dp, -0.02_dp, 0.2_dp, -0.3_dp, 5.0_dp, 6.2_dp]
        d = [(sin(real(j, dp))/10, j=1, 14)]
        w = [(cos(real(j, dp)), j=1, 6)]
        h = 1.0e-5_dp
        call leaf(v + h*d, plus)
        call leaf(v - h*d, minus)
        call ideal_axisymmetric_field_jvp(v(1), v(2), v(3), v(4), v(5), v(6), v(7), &
                v(8), v(9), v(10), v(11), v(12), v(13), v(14), d(1), d(2), d(3), d(4), &
                d(5), d(6), d(7), d(8), d(9), d(10), d(11), d(12), d(13), d(14), t(1), &
                                          t(2), t(3), t(4), t(5), t(6))
        if (maxval(abs(t - (plus - minus)/(2*h))) > 1.0e-8_dp) &
            error stop 'Piola JVP finite difference'
        call ideal_axisymmetric_field_vjp(v(1), v(2), v(3), v(4), v(5), v(6), v(7), &
                v(8), v(9), v(10), v(11), v(12), v(13), v(14), w(1), w(2), w(3), w(4), &
                     w(5), w(6), g(1), g(2), g(3), g(4), g(5), g(6), g(7), g(8), g(9), &
                                          g(10), g(11), g(12), g(13), g(14))
        if (abs(dot_product(t, w) - dot_product(d, g)) > 1.0e-12_dp) &
            error stop 'Piola JVP/VJP pairing'
    end subroutine
    subroutine leaf(v, y)
        real(dp), intent(in) :: v(14)
        real(dp), intent(out) :: y(6)
        call ideal_axisymmetric_field(v(1), v(2), v(3), v(4), v(5), v(6), v(7), &
                v(8), v(9), v(10), v(11), v(12), v(13), v(14), y(1), y(2), y(3), y(4), &
                                      y(5), y(6))
    end subroutine
end program
