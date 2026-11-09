program test_field_shell
    ! milca_field alone: a two-cell P1 shell whose map moves with phi (n = 5)
    ! carries the potential Phi = c Z + f_vac phi, so B = (0, f_vac/R, c)
    ! exactly. Points from pseudo-random master points are located and
    ! evaluated; a point beyond the shell must be reported outside.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use milca_equilibrium_chart, only: cylindrical_chart_t
    use milca_field_shell, only: field_shell_eval, field_shell_locate
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp), c = 0.7_dp, f_vac = 3.0_dp
    integer, parameter :: nfp = 5
    type(cylindrical_chart_t) :: chart
    integer :: conn(3, 2), node, i, e, f, st, nseed, seed
    integer, allocatable :: put(:)
    real(dp) :: reference(2, 4), xs(3, 3, 4), phi(3, 4), u(3), xi(2), xj(2)
    real(dp) :: x(3), y(3), b(3), ph, ph2, miss, worst
    reference = reshape([2.0_dp, -1.0_dp, 4.0_dp, -1.0_dp, 4.0_dp, 1.0_dp, &
        2.0_dp, 1.0_dp], [2, 4])
    conn = reshape([1, 2, 3, 1, 3, 4], [3, 2])
    xs = 0
    do node = 1, 4
        xs(1:2, 1, node) = reference(:, node)
        xs(1, 2, node) = 0.05_dp*reference(2, node)
        xs(2, 3, node) = 0.04_dp*reference(1, node)
    end do
    ! Nodal values of c Z per mode: P1 reproduces c Z of the mapped point.
    phi = c*xs(2, :, :)
    call random_seed(size=nseed)
    put = [(7919*i + 3, i = 1, nseed)]
    call random_seed(put=put)
    miss = 0; worst = 0; seed = 0
    do i = 1, 1000
        call random_number(u)
        e = 1 + int(2*u(1))
        xi = [u(2), u(3)*(1 - u(2))]
        ph = 2*pi*u(3)
        call field_shell_eval(chart, nfp, conn(:, e), reference, xs, phi, f_vac, &
            0.0_dp, [0.0_dp, 0.0_dp], xi, ph, x, b, st)
        if (st /= 0) error stop 'eval'
        call field_shell_locate(chart, nfp, conn, xs, x, seed, f, xj, ph2, st)
        if (st /= 0) error stop 'locate'
        seed = f
        call field_shell_eval(chart, nfp, conn(:, f), reference, xs, phi, f_vac, &
            0.0_dp, [0.0_dp, 0.0_dp], xj, ph2, y, b, st)
        miss = max(miss, maxval(abs(y - x)))
        worst = max(worst, maxval(abs(b - [0.0_dp, f_vac/x(1), c])))
    end do
    print '(a,es10.3,a,es10.3)', 'position miss', miss, '  max |B - B_exact| (T)', worst
    if (miss > 1.0e-12_dp) error stop 'located position'
    if (worst > 1.0e-13_dp) error stop 'field'
    call field_shell_locate(chart, nfp, conn, xs, [6.0_dp, 0.0_dp, 0.3_dp], 0, f, xj, &
        ph2, st)
    if (st == 0) error stop 'point beyond the shell located'
end program test_field_shell
