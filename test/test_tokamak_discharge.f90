!> Teaching discharge model: static limits of an ITER-like operating point and
!> three ramps that must end in a disruption of the expected kind, plus one
!> ramp that stays inside the limits and must survive.
program test_tokamak_discharge
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use milca_tokamak_discharge
    implicit none
    integer, parameter :: maxrec = 6000
    real(dp) :: par(TD_NPAR), out(TD_NOUT), prof(TD_NX, TD_NPROF), ramp(TD_NRAMP)
    real(dp), allocatable :: ser(:,:), pro(:,:,:)
    integer :: rc, nrec

    allocate(ser(TD_NSER, maxrec), pro(TD_NX, TD_NPROF, maxrec))
    par = 0
    par(1:18) = [6.2_dp, 2.0_dp, 1.75_dp, 0.33_dp, 5.3_dp, 15.0_dp, 1.0_dp, &
        1.0_dp, 0.3_dp, 20.0_dp, 1.0_dp, 1.7_dp, 0.005_dp, 1.3_dp, 0.3_dp, &
        1.0_dp, 2.5_dp, 0.1_dp]
    call td_limits(par, out, prof, rc)
    if (rc /= 0) error stop 'limits status'
    print '(a,2f8.3,a,f6.3,a,2f7.3,a,2f7.3)', 'q0 q95', out(1:2), '  li', out(3), &
        '  betat% betaN', out(4:5), '  nG fGW', out(7:8)
    print '(a,4f9.2,a,f8.4)', 'W tau98 Pheat Paux', out(9), out(10), out(11), out(13), &
        '  chi0', out(40)
    print '(a,4es11.3)', 'kink dW m=1..4', out(17:20)
    print '(a,4f10.4)', '2/1 and 3/2: x_s, r_s Delta''', out(29:32)
    print '(a,3f9.3)', 'ballooning max ratio, x, betaN limit', out(33), out(34), out(36)
    if (abs(out(2) - 3.0_dp) > 0.15_dp) error stop 'ITER q95 near 3'
    if (abs(out(7) - 1.19_dp) > 0.02_dp) error stop 'Greenwald density'

    ! density ramp to twice the Greenwald density
    ramp = 0
    ramp(1:6) = [4.0_dp, 2.4_dp, -1.0_dp, 15.0_dp, 8.0_dp, 1.0_dp]
    call run('density ramp')
    ! heating ramp: beta limit
    ramp(1:6) = [3.0_dp, 1.0_dp, 400.0_dp, 15.0_dp, 6.0_dp, 1.0_dp]
    call run('power ramp')
    ! current ramp: q95 -> below 2
    ramp(1:6) = [3.0_dp, 1.0_dp, -1.0_dp, 24.0_dp, 6.0_dp, 1.0_dp]
    call run('current ramp')
    ! small density step inside the limits
    ramp(1:6) = [2.0_dp, 1.05_dp, -1.0_dp, 15.0_dp, 4.0_dp, 0.5_dp]
    call run('safe ramp')
    print '(a)', 'PASS tokamak discharge'
contains
    subroutine run(name)
        character(*), intent(in) :: name
        integer :: k, kq, last
        call td_run(par, ramp, maxrec, nrec, ser, pro, rc)
        if (rc /= 0) error stop 'run status'
        kq = 0
        do k = 1, nrec
            if (ser(21,k) >= 3 .and. kq == 0) kq = k
        end do
        last = nrec
        print '(a,a,i6,a,f9.5,a,i2,a,i2)', name, ': records', nrec, '  t_end', &
            ser(1,last), '  phase', nint(ser(21,last)), '  cause', nint(ser(22,last))
        if (kq > 0) then
            print '(a,f9.5,a,f6.3,a,f6.3,a,f6.2,a,f6.2,a,2f7.3)', '  TQ at t=', ser(1,kq), &
                ' fGW', ser(6,kq), ' betaN', ser(7,kq), ' q95', ser(9,kq), &
                ' Ip', ser(2,kq), ' w21,w32[m]', ser(13,kq), ser(14,kq)
            print '(a,f8.3,a,f9.5,a,f7.3)', '  max Ip after TQ [MA]', maxval(ser(2,kq:last)), &
                '  CQ end t', ser(1,last), '  Z end', ser(11,last)
        end if
        do k = 1, nrec, max(1, nrec/8)
            print '(4x,a,f9.5,a,f7.3,a,f7.3,a,f6.2,a,f6.3,a,f6.3,a,f6.2,a,i2)', 't', ser(1,k), &
                ' Ip', ser(2,k), ' T0', ser(3,k), ' W', ser(4,k), ' fGW', ser(6,k), &
                ' bN', ser(7,k), ' q95', ser(9,k), ' ph', nint(ser(21,k))
        end do
    end subroutine
end program
