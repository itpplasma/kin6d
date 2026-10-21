program test_ideal_mhd_coordinates
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use ideal_mhd_energy_generated, only: generated_two_energy, &
        generated_perpendicular_energy, generated_perpendicular_matrix, &
        generated_perpendicular_jvp, generated_perpendicular_vjp
    use kin6d_ideal_xi, only: xi_bilinear
    implicit none
    real(dp) :: w(1), m(2,2), jv(2), vg(8), x(7), dx(7), expected, mx(2)
    complex(dp) :: a(3,3)
    integer :: orientation
    ! Unit uniform-B slab, xi=(u(z),v(z),0). Independent curl/Hessian
    ! gives |u_z|^2+|v_z|^2. The signed chart only changes curl signs.
    do orientation = -1, 1, 2
        call generated_two_energy(1.0_dp,0.0_dp,0.0_dp,0.0_dp, &
            1.0_dp,0.0_dp,real(orientation,dp),1.0_dp,1.0_dp,0.0_dp, &
            0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.4_dp,0.0_dp,0.0_dp, &
            0.7_dp,0.0_dp,-0.3_dp,w)
        call close(w(1),0.7_dp**2+0.3_dp**2,'slab curl energy')
        a=cmplx(0.0_dp,0.0_dp,dp)
        a(1,3)=0.7_dp; a(2,3)=-0.3_dp
        call close(real(xi_bilinear([0.0_dp,0.0_dp,1.0_dp],0.0_dp,1.0_dp,a,a),dp), &
            w(1),'Cartesian Hessian = coordinate leaf')
        x(1)=real(orientation,dp)
        x(2:)=[1.0_dp,1.0_dp,0.0_dp,2.3_dp,0.4_dp,-0.3_dp]
        dx=[0.1_dp,0.2_dp,-0.1_dp,0.3_dp,0.2_dp,-0.2_dp,0.4_dp]
        call generated_perpendicular_energy(x(1),x(2),x(3),x(4),x(5),x(6),x(7),w)
        expected=0.5_dp*x(5)*(x(6)**2+x(7)**2)
        call close(w(1),expected,'Cartesian kinetic norm')
        call generated_perpendicular_matrix(x(1),x(2),x(3),x(4),x(5),m)
        mx=matmul(m,x(6:7))
        call close(dot_product(x(6:7),mx)/2,expected,'mass Hessian')
        call generated_perpendicular_jvp(x(1),x(2),x(3),x(4),x(5),x(6),x(7), &
            dx(1),dx(2),dx(3),dx(4),dx(5),dx(6),dx(7),jv)
        call generated_perpendicular_vjp(x(1),x(2),x(3),x(4),x(5),x(6),x(7),1.0_dp,vg)
        call close(jv(1),expected,'JVP primal')
        call close(jv(2),dot_product(vg(2:),dx),'JVP/VJP duality')
    end do
    print '(a)', 'PASS shared slab delta-W, Cartesian kinetic norm and derivatives'
contains
    subroutine close(actual,reference,label)
        real(dp), intent(in) :: actual,reference
        character(*), intent(in) :: label
        if(abs(actual-reference)>256*epsilon(1.0_dp)*max(1.0_dp,abs(reference))) then
            print *,label,actual,reference
            error stop 1
        endif
    end subroutine close
end program test_ideal_mhd_coordinates
