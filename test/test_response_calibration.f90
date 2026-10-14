program test_response_calibration
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use kin6d_response_calibration, only: spot_weights, spot_jump, &
        spot_jump_vjp, equivalent_sheet_field, spot_distance
    use kin6d_response_calibration_leaves, only: sheet_field, &
        sheet_field_jvp, sheet_field_vjp, spot_distance_jvp, spot_distance_vjp
    implicit none
    real(dp) :: x(5, 2), w(5, 2), v(7), d(7), vp(7), vm(7), adj(7)
    real(dp) :: yp(2), ym(2), dy(2), seed(2), h, t
    real(dp) :: xs, xp, xm, dxs, location_adj(3)
    complex(dp) :: phi(5, 2), direction(5, 2), pa(5, 2), j, b
    integer :: i, side, st
    call spot_distance(-2, -0.5_dp, 0.1_dp, xs, st)
    call check(st == 0 .and. abs(xs - 0.1_dp) < 1.0e-15_dp, 'spot distance')
    h = 1.0e-5_dp
    call spot_distance(2, 0.5_dp + h*0.13_dp, 0.1_dp - h*0.02_dp, xp, st)
    call spot_distance(2, 0.5_dp - h*0.13_dp, 0.1_dp + h*0.02_dp, xm, st)
    call spot_distance_jvp(2.0_dp, 0.5_dp, 0.1_dp, 0.0_dp, 0.13_dp, &
        -0.02_dp, dxs)
    call spot_distance_vjp(2.0_dp, 0.5_dp, 0.1_dp, 0.3_dp, &
        location_adj(1), location_adj(2), location_adj(3))
    call check(abs(dxs - (xp - xm)/(2*h)) < 1.0e-10_dp, 'spot distance JVP')
    call check(abs(0.3_dp*dxs - 0.13_dp*location_adj(2) + &
        0.02_dp*location_adj(3)) < 1.0e-14_dp, 'spot distance VJP')
    do side = 1, 2
        do i = 1, 5
            x(i, side) = 0.5_dp + merge(-0.1_dp, 0.1_dp, side == 1) &
                + (i - 3)*0.01_dp
            t = x(i, side)
            phi(i, side) = cmplx(t**4, 2*t**3, dp)
            direction(i, side) = cmplx(i + side, i - side, dp)
        end do
    end do
    call spot_weights(x, 0.5_dp, -2, -0.5_dp, 0.1_dp, w, st)
    call check(st == 0, 'signed n and shear produce positive spot distance')
    j = spot_jump(w, phi)
    call check(abs(j - cmplx(4*(0.6_dp**3 - 0.4_dp**3), &
        6*(0.6_dp**2 - 0.4_dp**2), dp)) < 1.0e-12_dp, 'quartic flux derivative')
    call spot_jump_vjp(w, (0.3_dp, -0.7_dp), pa)
    call check(abs(real(conjg(cmplx(0.3_dp, -0.7_dp, dp))* &
        spot_jump(w, direction)) - real(sum(conjg(pa)*direction))) &
        < 1.0e-11_dp, 'fixed-grid Hermitian transpose')
    call spot_weights(x, 0.5_dp, 0, 1.0_dp, 0.1_dp, w, st)
    call check(st /= 0, 'reject n zero')
    call spot_weights(x, 0.5_dp, 2, 0.5_dp, 0.5_dp, w, st)
    call check(st /= 0, 'reject spot outside resolved samples')
    x(5, 1) = 0.51_dp
    call spot_weights(x, 0.5_dp, 2, 0.5_dp, 0.1_dp, w, st)
    call check(st /= 0, 'reject stencil crossing sheet')
    v = [1.2_dp, -0.3_dp, 2.1_dp, 0.4_dp, -0.7_dp, 0.2_dp, 1.8_dp]
    d = [0.1_dp, 0.2_dp, -0.1_dp, 0.3_dp, -0.2_dp, 0.4_dp, 0.5_dp]
    seed = [0.4_dp, -0.8_dp]
    h = 1.0e-5_dp
    vp = v + h*d; vm = v - h*d
    call sheet_field(vp(1), vp(2), vp(3), vp(4), vp(5), vp(6), vp(7), yp(1), yp(2))
    call sheet_field(vm(1), vm(2), vm(3), vm(4), vm(5), vm(6), vm(7), ym(1), ym(2))
    call sheet_field_jvp(v(1), v(2), v(3), v(4), v(5), v(6), v(7), &
        d(1), d(2), d(3), d(4), d(5), d(6), d(7), dy(1), dy(2))
    call sheet_field_vjp(v(1), v(2), v(3), v(4), v(5), v(6), v(7), &
        seed(1), seed(2), adj(1), adj(2), adj(3), adj(4), adj(5), adj(6), adj(7))
    call check(maxval(abs(dy - (yp - ym)/(2*h))) < 1.0e-10_dp, 'normalization JVP')
    call check(abs(sum(seed*dy) - sum(adj*d)) < 1.0e-14_dp, 'normalization VJP')
    call equivalent_sheet_field((2.0_dp, -3.0_dp), (1.0_dp, 0.0_dp), &
        (-4.0_dp, 0.0_dp), 2.0_dp, b, st)
    call check(st == 0 .and. abs(b - (4.0_dp, -6.0_dp)) < 1.0e-14_dp, &
        'peak sign and area against exact unit sheet')
    call equivalent_sheet_field((1.0_dp, 0.0_dp), (0.0_dp, 0.0_dp), &
        (-4.0_dp, 0.0_dp), 2.0_dp, b, st)
    call check(st /= 0, 'reject zero unit jump')
    call equivalent_sheet_field(cmplx(ieee_value(0.0_dp, ieee_quiet_nan), &
        0.0_dp, dp), (1.0_dp, 0.0_dp), (-4.0_dp, 0.0_dp), 2.0_dp, b, st)
    call check(st /= 0, 'reject nonfinite response')
    print '(a)', 'spot calibration algebra: PASS'
contains
    subroutine check(ok, message)
        logical, intent(in) :: ok
        character(*), intent(in) :: message
        if (.not. ok) then
            print '(a)', 'FAIL: '//message
            error stop 1
        end if
    end subroutine check
end program test_response_calibration
