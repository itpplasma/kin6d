program test_s0_gs_space
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_s0_gs_space_leaves, only: s0_polar_jet, s0_polar_jet_jvp, &
        s0_polar_jet_vjp, s0_polar_geometry, s0_polar_geometry_jvp, &
        s0_polar_geometry_vjp
    implicit none
    real(dp) :: x(9), dx(9), seed(18), y(18), tangent(18), adjoint(9)
    real(dp) :: plus(18), minus(18), step, pr, pz, jac, geom(3), dg(3), bg(4)
    real(dp) :: xp, zp, r, theta, values(6)
    integer :: i, direction
    x = [0.7_dp,0.4_dp,1.0_dp,0.9_dp,-0.3_dp,0.2_dp, &
        0.6_dp,-0.4_dp,0.1_dp]
    dx = sin([(real(i,dp),i=1,9)])
    seed = cos([(real(i,dp),i=1,18)])
    call jet(x,y)
    call s0_polar_jet_jvp(x(1),x(2),x(3),x(4),x(5),x(6),x(7),x(8),x(9), &
        dx(1),dx(2),dx(3),dx(4),dx(5),dx(6),dx(7),dx(8),dx(9), &
        tangent(1),tangent(2),tangent(3),tangent(4),tangent(5),tangent(6), &
        tangent(7),tangent(8),tangent(9),tangent(10),tangent(11),tangent(12), &
        tangent(13),tangent(14),tangent(15),tangent(16),tangent(17),tangent(18))
    call s0_polar_jet_vjp(x(1),x(2),x(3),x(4),x(5),x(6),x(7),x(8),x(9), &
        seed(1),seed(2),seed(3),seed(4),seed(5),seed(6),seed(7),seed(8), &
        seed(9),seed(10),seed(11),seed(12),seed(13),seed(14),seed(15), &
        seed(16),seed(17),seed(18),adjoint(1),adjoint(2),adjoint(3), &
        adjoint(4),adjoint(5),adjoint(6),adjoint(7),adjoint(8),adjoint(9))
    step = 1e-6_dp
    call jet(x+step*dx,plus); call jet(x-step*dx,minus)
    call expect(maxval(abs((plus-minus)/(2*step)-tangent)) < 2e-9_dp, &
        'polar primal/JVP central difference')
    call expect(abs(dot_product(seed,tangent)-dot_product(dx,adjoint)) &
        < 1e-12_dp,'polar JVP/VJP duality')
    ! Independent physical finite difference, without using the emitted jets.
    xp = x(1)*cos(x(2)); zp = x(1)*sin(x(2))
    do direction = 1, 2
        call physical_value(xp+step*merge(1,0,direction==1), &
            zp+step*merge(1,0,direction==2),values)
        plus(:6) = values
        call physical_value(xp-step*merge(1,0,direction==1), &
            zp-step*merge(1,0,direction==2),values)
        minus(:6) = values
        i = 6*direction
        call expect(maxval(abs((plus(:6)-minus(:6))/(2*step)-y(i+1:i+6))) &
            < 2e-9_dp,'physical R/Z polar jet')
    end do
    ! The combined m=1 axis basis is a constant Cartesian transverse vector.
    do i = 1, 8
        x = [1e-8_dp,real(i,dp),1.0_dp,1.0_dp,1.0_dp,0.0_dp, &
            0.0_dp,0.0_dp,0.0_dp]
        call jet(x,y)
        call expect(maxval(abs(y(:6)-[1.0_dp,0.0_dp,0.0_dp, &
            0.0_dp,0.0_dp,1.0_dp])) < 1e-14_dp,'regular m=1 axis value')
        call expect(maxval(abs(y(7:))) < 1e-7_dp,'regular m=1 axis jet')
    end do
    r = 0.7_dp; theta = 0.4_dp
    call s0_polar_geometry(r,theta,20.0_dp,0.2_dp,pr,pz,jac)
    call s0_polar_geometry_jvp(r,theta,20.0_dp,0.2_dp, &
        0.3_dp,-0.2_dp,0.1_dp,0.4_dp,dg(1),dg(2),dg(3))
    call s0_polar_geometry_vjp(r,theta,20.0_dp,0.2_dp, &
        0.2_dp,-0.3_dp,0.4_dp,bg(1),bg(2),bg(3),bg(4))
    call expect(abs(dot_product(dg,[0.2_dp,-0.3_dp,0.4_dp])- &
        dot_product(bg,[0.3_dp,-0.2_dp,0.1_dp,0.4_dp])) < 1e-12_dp, &
        'polar geometry JVP/VJP duality')
    call s0_polar_geometry(r+0.3_dp*step,theta-0.2_dp*step, &
        20.0_dp+0.1_dp*step,0.2_dp+0.4_dp*step,geom(1),geom(2),geom(3))
    call expect(maxval(abs((geom-[pr,pz,jac])/step-dg)) < 2e-6_dp, &
        'polar geometry directional derivative')
    print '(a)', '[V] physical polar jets, regular m=1 axis, primal/JVP/VJP'
contains
    subroutine expect(ok,label)
        logical, intent(in) :: ok
        character(*), intent(in) :: label
        if (.not. ok) then
            print '(a)', 'FAIL: '//label
            error stop 1
        end if
    end subroutine expect

    subroutine jet(input,output)
        real(dp), intent(in) :: input(9)
        real(dp), intent(out) :: output(18)
        call s0_polar_jet(input(1),input(2),input(3),input(4),input(5),input(6), &
            input(7),input(8),input(9),output(1),output(2),output(3),output(4), &
            output(5),output(6),output(7),output(8),output(9),output(10), &
            output(11),output(12),output(13),output(14),output(15),output(16), &
            output(17),output(18))
    end subroutine jet

    subroutine physical_value(radial,vertical,output)
        real(dp), intent(in) :: radial,vertical
        real(dp), intent(out) :: output(6)
        real(dp) :: rr, angle, amplitude(3)
        complex(dp) :: phase, xi(3)
        rr = sqrt(radial**2+vertical**2); angle = atan2(vertical,radial)
        amplitude = x(4:6)+x(7:9)*(rr-x(1))
        phase = exp(cmplx(0.0_dp,x(3)*angle,dp))
        xi = phase*[cmplx(radial*amplitude(1),-vertical*amplitude(2),dp)/rr, &
            cmplx(0.0_dp,amplitude(3),dp), &
            cmplx(vertical*amplitude(1),radial*amplitude(2),dp)/rr]
        output = [real(xi,dp),aimag(xi)]
    end subroutine physical_value
end program test_s0_gs_space
