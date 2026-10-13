program test_gs_steklov_cache
    !! The cached Steklov-Poincare operator of a fixed curve Gamma_c: a repeat
    !! call returns the assembled matrix, and the same curve entered at another
    !! panel (as after a mesh refit) returns the relabelled matrix of a fresh
    !! assembly; a different curve is not served from the cache.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortfem_axisymmetric_gs_boundary_operators_2d, only: gs_bem_curve_t
    use kin6d_gs_coupling, only: gs_coupling_steklov, gs_coupling_steklov_cache_clear
    implicit none

    integer, parameter :: np = 10, nn = 2*np, degree = 2
    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp) :: nodes(2, nn), shifted(2, nn), radius
    real(dp), allocatable :: s_first(:, :), s_hit(:, :), s_shift(:, :), s_fresh(:, :)
    type(gs_bem_curve_t) :: curve
    integer :: stat, j, shift, i
    real(dp) :: err, scale

    call curve_nodes(0.7_dp, nodes)
    call gs_coupling_steklov_cache_clear()
    call gs_coupling_steklov(nodes, degree, curve, s_first, stat)
    call require(stat == 0, 'first assembly')
    call gs_coupling_steklov(nodes, degree, curve, s_hit, stat)
    call require(stat == 0, 'cache hit')
    call require(all(s_hit == s_first), 'a hit returns the cached matrix')

    shift = 6
    do j = 1, nn
        shifted(:, j) = nodes(:, mod(j - 1 + shift, nn) + 1)
    end do
    call gs_coupling_steklov(shifted, degree, curve, s_shift, stat)
    call require(stat == 0, 'shifted entry')
    call gs_coupling_steklov_cache_clear()
    call gs_coupling_steklov(shifted, degree, curve, s_fresh, stat)
    call require(stat == 0, 'fresh assembly of the shifted curve')
    scale = maxval(abs(s_fresh))
    err = 0.0_dp
    do j = 1, nn
        do i = 1, nn
            err = max(err, abs(s_shift(i, j) - s_fresh(i, j)))
        end do
    end do
    call require(err <= 1.0e-11_dp*scale, 'relabelled matrix equals a fresh assembly')

    radius = 0.8_dp
    call curve_nodes(radius, shifted)
    call gs_coupling_steklov(shifted, degree, curve, s_hit, stat)
    call require(stat == 0, 'different curve')
    call require(maxval(abs(s_hit - s_fresh)) > 1.0e-3_dp*scale, &
                 'a different curve is assembled, not served from the cache')
    print '(a)', 'gs_steklov_cache oracles passed'

contains

    subroutine curve_nodes(a, x)
        real(dp), intent(in) :: a
        real(dp), intent(out) :: x(:, :)
        integer :: k
        real(dp) :: t

        do k = 1, nn
            t = 2.0_dp*pi*real(k - 1, dp)/real(nn, dp)
            x(1, k) = 3.0_dp + a*cos(t)
            x(2, k) = 1.1_dp*a*sin(t)
        end do
    end subroutine curve_nodes

    subroutine require(condition, label)
        logical, intent(in) :: condition
        character(*), intent(in) :: label

        if (.not. condition) then
            print '(2a)', 'FAILED: ', label
            error stop 1
        end if
    end subroutine require

end program test_gs_steklov_cache
