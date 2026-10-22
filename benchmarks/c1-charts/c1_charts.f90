module c1_problem
    ! Anisotropic diffusion -div(K grad u) = f with K = chi_perp I +
    ! (chi_par - chi_perp) b b^T, b from the nested field and the manufactured
    ! flux function u = cos(pi s / 2). Since b . grad s = 0 exactly,
    ! f = -chi_perp Lap u does not depend on chi_par: every chi_par dependence
    ! of the discrete error is parallel-to-perpendicular pollution.
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use kin6d_nested_chart, only: nested_map_t, nested_point_t, nested_map_eval, &
        nested_map_invert, nested_field_direction, nested_forcing
    implicit none
    private
    public :: map, chi_par, chi_perp, coefficient, forcing, exact, coef_seconds, &
        coef_calls, invert_failures
    type(nested_map_t) :: map
    real(dp) :: chi_par = 1, chi_perp = 1, coef_seconds = 0
    integer(int64) :: coef_calls = 0, invert_failures = 0
    real(dp) :: warm(3) = 0
    real(dp), parameter :: pi = acos(-1.0_dp)
    !$omp threadprivate(warm)
contains
    subroutine locate(x, q, pt)
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: q(3)
        type(nested_point_t), intent(out) :: pt
        integer :: status
        q = warm
        call nested_map_invert(map, x, q, status)
        if (status /= 0) invert_failures = invert_failures + 1
        warm = q
        call nested_map_eval(map, q(1), q(2), q(3), pt)
    end subroutine

    subroutine coefficient(x, k, c)
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: k(:, :), c
        type(nested_point_t) :: pt
        real(dp) :: q(3), b(3)
        integer(int64) :: t0, t1, rate
        integer :: i
        call system_clock(t0, rate)
        call locate(x, q, pt)
        b = nested_field_direction(map, pt, q)
        do i = 1, 3
            k(:, i) = (chi_par - chi_perp)*b*b(i)
            k(i, i) = k(i, i) + chi_perp
        end do
        c = 0
        call system_clock(t1)
        coef_seconds = coef_seconds + real(t1 - t0, dp)/rate
        coef_calls = coef_calls + 1
    end subroutine

    function forcing(x) result(f)
        real(dp), intent(in) :: x(:)
        real(dp) :: f, q(3), grad(3), lap, s
        type(nested_point_t) :: pt
        call locate(x, q, pt)
        call nested_forcing(pt, q, grad, lap)
        s = q(1)**2 + q(2)**2
        f = -chi_perp*(-(pi/2)**2*cos(pi*s/2)*dot_product(grad, grad) &
                       - (pi/2)*sin(pi*s/2)*lap)
    end function

    function exact(x) result(u)
        real(dp), intent(in) :: x(:)
        real(dp) :: u, q(3)
        type(nested_point_t) :: pt
        call locate(x, q, pt)
        u = cos(pi*(q(1)**2 + q(2)**2)/2)
    end function
end module

program c1_charts
    ! Usage: c1_charts MAP {aligned|unaligned} P N NPHI CHI_RATIO [TOL] [RHO_CORE]
    ! Prints one CSV row (header with --header as the only argument).
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use kin6d_nested_chart, only: nested_map_read, build_ogrid_chart, &
        CHART_ALIGNED, CHART_UNALIGNED
    use fortfem_gll_line, only: gll_rule
    use fortfem_tensor_mesh, only: tensor_mesh_t, build_tensor_mesh
    use fortfem_tensor_operator, only: tensor_operator_t, build_tensor_operator, &
        tensor_load, tensor_l2_error
    use fortfem_tensor_solver, only: tensor_solve_info_t, solve_tensor_diffusion
    use c1_problem, only: map, chi_par, chi_perp, coefficient, forcing, exact, &
        coef_seconds, coef_calls, invert_failures
    implicit none
    character(512) :: arg, path
    character(16) :: chart
    integer :: p, n, nphi, kind, status, e, a
    real(dp) :: ratio, tol, rho_core, t_op, t_solve, t_coef, l2, unorm, maxerr, umax
    real(dp), allocatable :: nodes(:), weights(:), geometry(:, :, :), rhs(:), x(:), &
                             zero(:)
    integer, allocatable :: vertices(:, :)
    logical, allocatable :: outer(:)
    type(tensor_mesh_t) :: mesh
    type(tensor_operator_t) :: op
    type(tensor_solve_info_t) :: info
    integer(int64) :: t0, t1, rate
    if (command_argument_count() == 1) then
        print '(a)', 'chart,p,n,nphi,ndof,nelem,chi_ratio,tol,status,iterations,'// &
            'residual,t_load_s,t_solve_s,t_coef_in_solve_s,coef_calls,'// &
            'invert_failures,l2_rel,max_nodal_rel'
        stop
    end if
    if (command_argument_count() < 6) error stop 'usage: c1_charts MAP CHART P N NPHI CHI [TOL] [RHO_CORE]'
    call get_command_argument(1, path)
    call get_command_argument(2, chart)
    p = int_arg(3); n = int_arg(4); nphi = int_arg(5)
    ratio = real_arg(6, -1.0_dp)
    tol = real_arg(7, 1.e-12_dp)
    rho_core = real_arg(8, 0.35_dp)
    select case (trim(chart))
    case ('aligned'); kind = CHART_ALIGNED
    case ('unaligned'); kind = CHART_UNALIGNED
    case default; error stop 'chart must be aligned or unaligned'
    end select
    call nested_map_read(trim(path), map, status)
    if (status /= 0) error stop 'map read failed'
    chi_perp = 1
    chi_par = ratio
    ! Isoparametric: geometry degree q = p on the MF1 GLL nodes.
    call gll_rule(p, nodes, weights, status)
    if (status /= 0) error stop 'gll_rule failed'
    call build_ogrid_chart(map, kind, n, nphi, nodes, rho_core, vertices, geometry, &
                           outer, status)
    if (status /= 0) error stop 'chart build failed'
    call build_tensor_mesh(3, p, p, vertices, geometry, mesh, status)
    if (status /= 0) error stop 'build_tensor_mesh failed'
    ! Dirichlet u = 0 on s = 1: the radial end face of the ring's outer cells.
    do e = 1, mesh%nelem
        if (.not. outer(e)) cycle
        do a = 1, (p + 1)**3
            if (modulo(a - 1, p + 1) == p) mesh%boundary(mesh%dof(a, e)) = .true.
        end do
    end do
    ! Load and error use an over-integrated operator (p + 3 Gauss-Lobatto points).
    call system_clock(t0, rate)
    call build_tensor_operator(mesh, coefficient, op, status, nq=p + 3)
    if (status /= 0) error stop 'build_tensor_operator failed'
    call tensor_load(op, forcing, rhs)
    call system_clock(t1)
    t_op = real(t1 - t0, dp)/rate
    allocate (x(mesh%ndof), zero(mesh%ndof))
    x = 0; zero = 0
    coef_seconds = 0
    call system_clock(t0)
    call solve_tensor_diffusion(mesh, coefficient, rhs, tol, x, info, maxiter=20000)
    call system_clock(t1)
    t_solve = real(t1 - t0, dp)/rate
    t_coef = coef_seconds
    l2 = tensor_l2_error(op, x, exact)
    unorm = tensor_l2_error(op, zero, exact)
    maxerr = 0; umax = 0
    do a = 1, mesh%ndof
        maxerr = max(maxerr, abs(x(a) - exact(mesh%coordinates(:, a))))
        umax = max(umax, abs(exact(mesh%coordinates(:, a))))
    end do
    print '(a,",",i0,",",i0,",",i0,",",i0,",",i0,",",es10.3,",",es9.2,",",i0,",",'// &
        'i0,",",es10.3,",",f0.4,",",f0.4,",",f0.4,",",i0,",",i0,",",es12.5,",",es12.5)', &
        trim(chart), p, n, nphi, mesh%ndof, mesh%nelem, ratio, tol, info%status, &
        info%iterations, info%residual, t_op, t_solve, t_coef, coef_calls, &
        invert_failures, l2/unorm, maxerr/umax
contains
    integer function int_arg(i)
        integer, intent(in) :: i
        call get_command_argument(i, arg)
        read (arg, *) int_arg
    end function

    real(dp) function real_arg(i, default)
        integer, intent(in) :: i
        real(dp), intent(in) :: default
        real_arg = default
        if (command_argument_count() < i) return
        call get_command_argument(i, arg)
        read (arg, *) real_arg
    end function
end program
