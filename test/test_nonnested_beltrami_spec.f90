module test_nonnested_beltrami_spec_callbacks
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_mhd_mesh, only: flux_mesh_t, flux_mesh_build, flux_mesh_locate
    use kin6d_equilibrium_nonnested, only: free_b_field_t, nn_build
    use kin6d_nonnested_mapped, only: nn_map_t, nn_map_build
    use kin6d_nonnested_beltrami, only: beltrami_system_t, beltrami_report_t, &
        beltrami_prepare, beltrami_solve
    use fortfem_meridian_feec_n, only: meridian_point_basis

    use kin6d_nonnested_mapped, only: nn_map_function_t
    implicit none
    private
    public :: run_test_nonnested_beltrami_spec
    real(dp), parameter :: pi = acos(-1.0_dp), r0 = 6.0_dp, mu = 1.1_dp
    type :: run_test_nonnested_beltrami_spec_context_t
        real(dp), pointer :: epsilon => null()
    end type

    type, extends(nn_map_function_t) :: run_test_nonnested_beltrami_spec_deform_t
        type(run_test_nonnested_beltrami_spec_context_t), pointer :: data => null()
    contains
        procedure :: evaluate => run_test_nonnested_beltrami_spec_deform_binding
    end type

contains
    subroutine run_test_nonnested_beltrami_spec()
    ! External synthetic SPEC oracle. Samples contain physical (R,Z,phi,B).
    implicit none
    type(flux_mesh_t) :: mesh
    type(free_b_field_t) :: field
    type(nn_map_t) :: map
    type(beltrami_system_t) :: system
    type(beltrami_report_t) :: report
    real(dp), allocatable :: rings(:), data(:), samples(:, :), values(:, :)
    real(dp), allocatable :: local(:, :)
    integer, allocatable :: cell(:)
    real(dp) , target :: epsilon, error(3), mod_error, scale, ref(2), theta, u, xi, &
        eta, phi
    real(dp) :: jac(2, 2), inv(2, 2), point(2), df(3, 3), det, b(3), bp(3), bm(3), &
        xp(2), xm(2)
    real(dp) :: db(3, 3), jx(2, 3), der(3, 3), curl(3), residual, step, total
    integer :: k, nr, modes, st, i, j, t, n, unit
    character(512) :: arg, path, output
        type(run_test_nonnested_beltrami_spec_context_t), target :: ctx
        type(run_test_nonnested_beltrami_spec_deform_t), target :: deform_operator
    call get_command_argument(6, output)
    call get_command_argument(1, path)
    call get_command_argument(2, arg); read (arg, *) epsilon
    call get_command_argument(3, arg); read (arg, *) nr
    call get_command_argument(4, arg); read (arg, *) k
    call get_command_argument(5, arg); read (arg, *) modes
    open (newunit=unit, file=trim(path), status='old'); read (unit, *) n
    allocate (samples(6, n), values(3, n))
    do i = 1, n
        read (unit, *) samples(:, i)
    end do
    close (unit)
    allocate (rings(nr)); rings = [(real(i, dp)/nr, i=1, nr)]
    call flux_mesh_build(rings, 8, 256, k, mapping, mesh, st)
    if (st /= 0) error stop 'SPEC mesh'
    call nn_build(mesh, k, modes, [r0, 0.0_dp], 1.0_dp, 0.0_dp, 1.0_dp, field, st)
    if (st /= 0) error stop 'SPEC field'
        call bind_run_test_nonnested_beltrami_spec_context()
    call nn_map_build(field%fs, 2, 9, deform_operator, map, st)
    if (st /= 0) error stop 'SPEC map'
    allocate (cell(n), local(2, n))
    do i = 1, n
        phi = samples(3, i)
        call deformation(phi, df, det)
        jac = reshape([df(1, 1), df(3, 1), df(1, 3), df(3, 3)], [2, 2])
        det = jac(1, 1)*jac(2, 2) - jac(1, 2)*jac(2, 1)
        inv = reshape([jac(2, 2), -jac(2, 1), -jac(1, 2), jac(1, 1)], [2, 2])/det
        ref = matmul(inv, samples(:2, i) - [r0, 0.0_dp]) + [r0, 0.0_dp]
        u = norm2(ref - [r0, 0.0_dp]); theta = atan2(ref(2), ref(1) - r0)
        call flux_mesh_locate(field%mesh, ref, u, theta, t, xi, eta, st)
        if (st /= 0) error stop 'SPEC locate'
        if (min(xi, eta, 1 - xi - eta) < 1.0e-7_dp) then
            print *, 'SPEC facet: sample/rings/cell/xi/eta/R/Z/phi ', &
                i, nr, t, xi, eta, samples(:3, i)
            error stop 'SPEC oracle point lies on an FEEC facet'
        end if
        cell(i) = t; local(:, i) = [xi, eta]
    end do
    allocate (data((2*modes + 1)*field%fs%sp(1)%ndof)); data = 0
    call beltrami_prepare(field, 1.0_dp, data, 4*modes + 9, system, st, map=map, nfp=2)
    if (st /= 0) error stop 'SPEC prepare'
    call beltrami_solve(system, mu, field, report)
    if (report%status /= 0) then
        print *, report%status, trim(report%message), report%min_pivot, &
            report%pivot_tolerance
        error stop 'SPEC Beltrami solve'
    end if
    residual = 0; total = 0; step = 1.0e-6_dp
    do i = 1, n
        t = cell(i); xi = local(1, i); eta = local(2, i); phi = samples(3, i)
        call sample(t, xi, eta, phi, b, point)
        values(:, i) = b
        ! Take the derivative probe inside this element. V2 has tangential
        ! jumps, so crossing an edge would not measure its element curl.
        xi = (1 - 1.0e-4_dp)*xi + 1.0e-4_dp/3
        eta = (1 - 1.0e-4_dp)*eta + 1.0e-4_dp/3
        call sample(t, xi, eta, phi, b, point)
        do j = 1, 3
            select case (j)
            case (1)
                call sample(t, xi + step, eta, phi, bp, xp)
                call sample(t, xi - step, eta, phi, bm, xm)
            case (2)
                call sample(t, xi, eta + step, phi, bp, xp)
                call sample(t, xi, eta - step, phi, bm, xm)
            case (3)
                call sample(t, xi, eta, phi + step, bp, xp)
                call sample(t, xi, eta, phi - step, bm, xm)
            end select
            db(:, j) = (bp - bm)/(2*step); jx(:, j) = (xp - xm)/(2*step)
        end do
        det = jx(1, 1)*jx(2, 2) - jx(1, 2)*jx(2, 1)
        inv = reshape([jx(2, 2), -jx(2, 1), -jx(1, 2), jx(1, 1)], [2, 2])/det
        der(:, :2) = matmul(db(:, :2), inv)
        der(:, 3) = db(:, 3) - matmul(der(:, :2), jx(:, 3))
        curl = [der(3, 3)/point(1) - der(2, 2), der(1, 2) - der(3, 1), &
            der(2, 1) + b(2)/point(1) - der(1, 3)/point(1)]
        residual = residual + sum((curl - mu*b)**2); total = total + sum((mu*b)**2)
    end do
    scale = sum(samples(4:6, :)**2)
    error = sqrt(sum((values - samples(4:6, :))**2, dim=2)/scale)
    mod_error = sqrt(sum((sqrt(sum(values**2, dim=1)) &
        - sqrt(sum(samples(4:6, :)**2, dim=1)))**2)/scale)
    print '(a,3i5,a,6es14.6)', 'SPEC rings/order/modes ', nr, k, modes, &
        ' BR/Bphi/BZ/modB/weak/element-curl ', error, mod_error, report%residual, &
        sqrt(residual/total)
    if (len_trim(output) > 0) then
        open (newunit=unit, file=trim(output), status='replace')
        write (unit, '(i0)') n
        do i = 1, n
            write (unit, '(6es25.16)') samples(:3, i), values(:, i)
        end do
        close (unit)
    end if
    if (sqrt(sum(error**2)) > max(0.003_dp, 0.04_dp/nr) &
        .or. report%residual > 1.0e-9_dp) &
        error stop 'SPEC comparison'
contains
    subroutine bind_run_test_nonnested_beltrami_spec_context()
        ctx%epsilon => epsilon
        deform_operator%data => ctx
    end subroutine bind_run_test_nonnested_beltrami_spec_context

    subroutine mapping(u, theta, x)
        real(dp), intent(in) :: u, theta
        real(dp), intent(out) :: x(2)
        x = [r0 + u*cos(theta), u*sin(theta)]
    end subroutine
    subroutine deformation(phi, df, det)
        real(dp), intent(in) :: phi
        real(dp), intent(out) :: df(3, 3), det
        real(dp) :: c, s
        c = epsilon*cos(2*phi); s = epsilon*sin(2*phi)
        df = 0; df(1, 1) = 1 + c; df(1, 3) = s; df(3, 1) = s; df(3, 3) = 1 - c
        df(2, 2) = 1; det = 1 - epsilon**2
    end subroutine
    subroutine sample(t, xi, eta, phi, b, y)
        integer, intent(in) :: t
        real(dp), intent(in) :: xi, eta, phi
        real(dp), intent(out) :: b(3), y(2)
        real(dp) :: x(2), basis(3, field%fs%sp(2)%nloc), hat(3), df(3, 3), det, c, s
        complex(dp) :: bn(3)
        integer :: st, n
        call meridian_point_basis(field%fs, field%pb, 2, t, xi, eta, x, basis, st)
        if (st /= 0) error stop 'SPEC sample'
        hat = 0
        do n = 0, modes
            bn = matmul(basis, field%b(field%fs%sp(2)%dof(:, t), n) &
                *field%fs%sp(2)%sgn(:, t))
            hat = hat + merge(1, 2, n == 0)*real(bn*exp(cmplx(0.0_dp, n*phi, dp)), dp)
        end do
        hat(2) = hat(2) + field%f_vac/x(1)
        call bind_run_test_nonnested_beltrami_spec_context()
        call run_test_nonnested_beltrami_spec_deform(ctx, x, phi, y)
        call deformation(phi, df, det)
        c = epsilon*cos(2*phi); s = epsilon*sin(2*phi)
        df(1, 2) = (-2*s*(x(1) - r0) + 2*c*x(2))/x(1)
        df(3, 2) = (2*c*(x(1) - r0) + 2*s*x(2))/x(1)
        df(2, 2) = y(1)/x(1); det = det*df(2, 2)
        b = matmul(df, hat)/det
    end subroutine

    end subroutine run_test_nonnested_beltrami_spec

    subroutine run_test_nonnested_beltrami_spec_deform(ctx, x, phi, y)
        type(run_test_nonnested_beltrami_spec_context_t), intent(in) :: ctx
        real(dp), intent(in) :: x(2), phi
        real(dp), intent(out) :: y(2)
        real(dp) :: c, s
        c = ctx%epsilon*cos(2*phi); s = ctx%epsilon*sin(2*phi)
        y = [r0 + (1 + c)*(x(1) - r0) + s*x(2), s*(x(1) - r0) + (1 - c)*x(2)]
    end subroutine run_test_nonnested_beltrami_spec_deform

    subroutine run_test_nonnested_beltrami_spec_deform_binding(self, x, phi, y)
        class(run_test_nonnested_beltrami_spec_deform_t), intent(in) :: self
        real(dp), intent(in) :: x(2), phi
        real(dp), intent(out) :: y(2)
        call run_test_nonnested_beltrami_spec_deform(self%data, x, phi, y)
    end subroutine
end module test_nonnested_beltrami_spec_callbacks

program test_nonnested_beltrami_spec
    use test_nonnested_beltrami_spec_callbacks, only: run_test_nonnested_beltrami_spec
    implicit none
    call run_test_nonnested_beltrami_spec()
end program test_nonnested_beltrami_spec
