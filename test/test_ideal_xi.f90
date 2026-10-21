program test_ideal_xi
    use, intrinsic :: iso_fortran_env, only: dp => real64, sp => real32
    use kin6d_ideal_xi, only: xi_bilinear, xi_gradient_n, xi_volume_energy, &
        xi_energy_report_t, XI_OK, XI_BAD_INPUT
    use kin6d_ideal_map_leaves, only: ideal_energy_hessian
    use kin6d_ideal_xi_leaves, only: xi_form_dp, xi_form_jvp_dp, xi_form_vjp_dp
    implicit none
    real(dp) :: h(9,9), packed(45), b(3), mu, p, expected, tolerance
    real(dp) :: x(23), dx(23), g(23), tangent, fd, step, energy
    real(dp) :: w(4), field(3,4), pressure(4), r, angle, qk
    complex(dp) :: a(3,3), c(3,3), value, xi(3), dr(3), dz(3), batch(3,3,4)
    complex(dp) :: gp(3), gm(3), phase_derivative(3), imag_unit
    complex(sp) :: as(3,3), cs(3,3), xs(3), rs(3), zs(3)
    type(xi_energy_report_t) :: report
    integer :: i, j, k, n, status
    tolerance = 512*epsilon(1.0_dp)
    b = [0.7_dp, -1.1_dp, 2.3_dp]; p = 0.4_dp; mu = 1.7_dp
    call ideal_energy_hessian(1.0_dp,0.0_dp,0.0_dp,0.0_dp,1.0_dp,0.0_dp, &
        0.0_dp,0.0_dp,1.0_dp,b(1),b(2),b(3),p,mu,1.0_dp, &
        packed(1),packed(2),packed(3),packed(4),packed(5),packed(6),packed(7), &
        packed(8),packed(9),packed(10),packed(11),packed(12),packed(13), &
        packed(14),packed(15),packed(16),packed(17),packed(18),packed(19), &
        packed(20),packed(21),packed(22),packed(23),packed(24),packed(25), &
        packed(26),packed(27),packed(28),packed(29),packed(30),packed(31), &
        packed(32),packed(33),packed(34),packed(35),packed(36),packed(37), &
        packed(38),packed(39),packed(40),packed(41),packed(42),packed(43), &
        packed(44),packed(45))
    k = 0
    do j = 1, 9
        do i = 1, j
            k = k + 1
            h(i,j) = packed(k); h(j,i) = packed(k)
        end do
    end do
    do j = 1, 9
        do i = 1, 9
            a = 0; c = 0
            a(1+mod(i-1,3),1+(i-1)/3) = 1
            c(1+mod(j-1,3),1+(j-1)/3) = 1
            value = xi_bilinear(b,p,mu,a,c)
            call close(real(value,dp), h(i,j), tolerance, 'emitted map Hessian = L4')
        end do
    end do
    ! General complex jets: Hermitian pairing and an independent tensor formula.
    do j = 1, 3
        do i = 1, 3
            a(i,j) = cmplx(sin(real(i+3*j,dp)),cos(real(2*i-j,dp)),dp)
            c(i,j) = cmplx(cos(real(i-2*j,dp)),sin(real(i+j,dp)),dp)
        end do
    end do
    value = xi_bilinear(b,p,mu,a,c)
    if (abs(value-conjg(xi_bilinear(b,p,mu,c,a))) > tolerance*max(1.0_dp,abs(value))) &
        error stop 'Hermitian pairing'
    imag_unit = cmplx(0,1,dp)
    if (abs(xi_bilinear(b,p,mu,imag_unit*a,c)-imag_unit*value) &
        > tolerance*max(1.0_dp,abs(value))) error stop 'trial complex linearity'
    if (abs(xi_bilinear(b,p,mu,a,imag_unit*c)+imag_unit*value) &
        > tolerance*max(1.0_dp,abs(value))) error stop 'test conjugate linearity'
    expected = oracle(b,p,mu,real(a,dp),real(c,dp)) &
        + oracle(b,p,mu,aimag(a),aimag(c))
    call close(real(value,dp), expected, tolerance, 'independent xi tensor form')
    as = cmplx(a,kind=sp); cs = cmplx(c,kind=sp)
    call close(real(xi_bilinear(real(b,sp),real(p,sp),real(mu,sp),as,cs),dp), &
        real(value,dp), 64.0_dp*epsilon(1.0_sp), 'kind-generic FP32 path')
    ! Independent finite differences and JVP/VJP duality on all 23 arguments.
    do k = 1, 23
        x(k) = sin(real(k,dp)); dx(k) = cos(real(2*k,dp))
    end do
    x(5) = 1.7_dp
    call derivative(x,dx,tangent,g)
    call close(tangent,dot_product(g,dx),tolerance,'JVP/VJP pairing')
    step = 2.0_dp**(-18)
    fd = (form(x+step*dx)-form(x-step*dx))/(2*step)
    call close(tangent,fd,2.0e-8_dp,'finite difference derivative')
    ! exp(+i n phi), with basis rotation checked independently in Cartesian space.
    imag_unit = cmplx(0,1,dp); r = 2.3_dp; n = -3
    xi = [cmplx(1,0.3_dp,dp),cmplx(-0.4_dp,0.7_dp,dp),cmplx(0.2_dp,-0.9_dp,dp)]
    dr = cmplx([0.2_dp,0.4_dp,0.6_dp],0.0_dp,dp)
    dz = cmplx([-0.3_dp,0.5_dp,0.1_dp],0.0_dp,dp)
    call xi_gradient_n(n,r,xi,dr,dz,a,status)
    if (status /= XI_OK) error stop 'cylindrical gradient status'
    angle = 1.0e-5_dp
    gp = rotate(xi*exp(imag_unit*n*angle),angle)
    gm = rotate(xi*exp(-imag_unit*n*angle),-angle)
    phase_derivative = (gp-gm)/(2*angle*r)
    if (maxval(abs(a(:,2)-phase_derivative)) > 2.0e-9_dp) &
        error stop 'Fourier sign or basis rotation'
    if (maxval(abs(a(:,1)-dr)) > tolerance) error stop 'radial gradient'
    if (maxval(abs(a(:,3)-dz)) > tolerance) error stop 'axial gradient'
    xs = cmplx(xi,kind=sp); rs = cmplx(dr,kind=sp); zs = cmplx(dz,kind=sp)
    call xi_gradient_n(n,real(r,sp),xs,rs,zs,as,status)
    if (status /= XI_OK) error stop 'FP32 gradient status'
    if (maxval(abs(cmplx(as,kind=dp)-a)) > 16*epsilon(1.0_sp)) &
        error stop 'FP32 gradient'
    ! Cartesian translation is n=1 in cylindrical components: grad xi=0.
    xi = [cmplx(1,0,dp),imag_unit,cmplx(0,0,dp)]; dr = 0; dz = 0
    call xi_gradient_n(1,r,xi,dr,dz,a,status)
    if (maxval(abs(a)) > tolerance) error stop 'translation null'
    call xi_gradient_n(1,0.0_dp,xi,dr,dz,a,status)
    if (status /= XI_BAD_INPUT) error stop 'axis unsupported'
    ! Screw-pinch rigid helical translation: raw map plasma Hessian = q^2.
    ! B=(-y,x,Bz), xi=(cos kz,sin kz,0); at kz=0 A_yz=k.
    ! This is the plasma-volume piece; the vacuum/boundary term is separate.
    do k = 1, 3
        qk = 0.5_dp*k
        b = [0.3_dp,-0.7_dp,2.0_dp]; a = 0; a(2,3) = qk/b(3)
        value = xi_bilinear(b,0.6_dp,1.0_dp,a,a)
        call close(real(value,dp),qk*qk,tolerance,'cylinder kink plasma limit')
    end do
    ! Axial shear-Alfven limit and the modal 1/2 convention in the batch consumer.
    b = [0.0_dp,0.0_dp,2.0_dp]; a = 0; a(1,3) = imag_unit*0.7_dp
    do k = 1, 4
        w(k) = 0.25_dp; field(:,k) = b; pressure(k) = p; batch(:,:,k) = a
    end do
    call xi_volume_energy(w,field,pressure,mu,batch,energy,report)
    if (report%status /= XI_OK) error stop 'volume status'
    if (report%claim /= 'candidate') error stop 'claim'
    if (report%bound_status /= 'unavailable') error stop 'bound status'
    call close(energy,0.5_dp*sum(w)*4*0.7_dp**2/mu,tolerance,'Alfven energy')
    w(1) = -1
    call xi_volume_energy(w,field,pressure,mu,batch,energy,report)
    if (report%status /= XI_BAD_INPUT) error stop 'negative quadrature weight'
    print *, 'PASS: xi Hessian, tensor oracle, Fourier phases, derivatives, cylinder and Alfven limits'
contains
    subroutine close(actual,wanted,tol,label)
        real(dp), intent(in) :: actual,wanted,tol
        character(*), intent(in) :: label
        if (abs(actual-wanted) > tol*max(1.0_dp,abs(wanted))) then
            print *, label, actual,wanted
            error stop 1
        end if
    end subroutine
    function oracle(b,p,mu,a,c) result(v)
        real(dp), intent(in) :: b(3),p,mu,a(3,3),c(3,3)
        real(dp) :: v,ta,tc
        ta = a(1,1)+a(2,2)+a(3,3); tc = c(1,1)+c(2,2)+c(3,3)
        v = dot_product(matmul(a,b),matmul(c,b))/mu &
            -(dot_product(b,matmul(a,b))*tc+dot_product(b,matmul(c,b))*ta)/mu &
            +(sum(b*b)/(2*mu)-p)*ta*tc &
            +(sum(b*b)/(2*mu)+p)*sum(a*transpose(c))
    end function
    function rotate(x,phi) result(v)
        complex(dp), intent(in) :: x(3)
        real(dp), intent(in) :: phi
        complex(dp) :: v(3)
        v = [cos(phi)*x(1)-sin(phi)*x(2), &
            sin(phi)*x(1)+cos(phi)*x(2),x(3)]
    end function
    function form(x) result(v)
        real(dp), intent(in) :: x(23)
        real(dp) :: v
        call xi_form_dp(x(1),x(2),x(3),x(4),x(5),x(6),x(7),x(8), &
            x(9),x(10),x(11),x(12),x(13),x(14),x(15),x(16),x(17), &
            x(18),x(19),x(20),x(21),x(22),x(23),v)
    end function
    subroutine derivative(x,d,v,g)
        real(dp), intent(in) :: x(23),d(23)
        real(dp), intent(out) :: v,g(23)
        call xi_form_jvp_dp(x(1),x(2),x(3),x(4),x(5),x(6),x(7),x(8), &
            x(9),x(10),x(11),x(12),x(13),x(14),x(15),x(16),x(17), &
            x(18),x(19),x(20),x(21),x(22),x(23),d(1),d(2),d(3),d(4), &
            d(5),d(6),d(7),d(8),d(9),d(10),d(11),d(12),d(13),d(14), &
            d(15),d(16),d(17),d(18),d(19),d(20),d(21),d(22),d(23),v)
        call xi_form_vjp_dp(x(1),x(2),x(3),x(4),x(5),x(6),x(7),x(8), &
            x(9),x(10),x(11),x(12),x(13),x(14),x(15),x(16),x(17), &
            x(18),x(19),x(20),x(21),x(22),x(23),1.0_dp, &
            g(1),g(2),g(3),g(4),g(5),g(6),g(7),g(8),g(9),g(10),g(11), &
            g(12),g(13),g(14),g(15),g(16),g(17),g(18),g(19),g(20),g(21),g(22),g(23))
    end subroutine
end program test_ideal_xi
