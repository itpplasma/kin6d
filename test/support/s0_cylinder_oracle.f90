module s0_cylinder_oracle
    ! Independent Eulerian force/energy form. Does not import S0, its generated
    ! Hessian, geometry, mass or basis leaves. Dense LAPACK is a test oracle only.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_quadrature, only: gauss_legendre_ab
    implicit none
    private
    public :: cylinder_oracle_t, cylinder_solve, cylinder_profile
    type :: cylinder_oracle_t
        real(dp) :: radius = 2, bz_axis = 0.5_dp, k = -0.4_dp
        real(dp) :: twist = 1, shear = 1, omega2 = 0
        integer :: cells = 0, negative = 0
        real(dp), allocatable :: u(:), v(:,:)
    end type cylinder_oracle_t
    interface
        subroutine dsygv(itype,jobz,uplo,n,a,lda,b,ldb,w,work,lwork,info)
            import :: dp
            integer, intent(in) :: itype,n,lda,ldb,lwork
            character, intent(in) :: jobz,uplo
            real(dp), intent(inout) :: a(lda,*),b(ldb,*)
            real(dp), intent(out) :: w(*),work(*)
            integer, intent(out) :: info
        end subroutine dsygv
    end interface
contains
    subroutine cylinder_solve(model, cells)
        type(cylinder_oracle_t), intent(inout) :: model
        integer, intent(in) :: cells
        real(dp), allocatable :: k(:,:), mass(:,:), eigen(:), work(:), x(:), wt(:)
        real(dp) :: phi(3), deriv(3), h, r, bz, bt, dbz, dbt, jp, jz, pp
        real(dp) :: u(6), du(6), v(6), w(6), div(6), qtheta(6), qz(6), qr(6)
        real(dp) :: f, helicity, integrand, factor, norm
        integer :: nd, dof(6), e, g, a, b, i, row, col, stat, node
        model%cells = cells
        h = model%radius/cells; nd = 5*cells-1
        allocate(k(nd,nd),mass(nd,nd),eigen(nd),work(8*nd),x(6),wt(6))
        k = 0; mass = 0
        call gauss_legendre_ab(6,0.0_dp,1.0_dp,x,wt)
        do e = 1, cells
            do a = 1, 3
                node = 2*(e-1)+a
                dof(a) = node
                if (node == 2*cells+1) dof(a) = 0
                dof(a+3) = 2*cells+3*(e-1)+a-1
            end do
            dof(4) = 2*cells+3*(e-1)
            if (e == 1) dof(4) = 0 ! Axis v(0)=u(0), represented by basis 1.
            do g = 1, size(x)
                r = (e-1+x(g))*h
                phi = [(1-x(g))*(1-2*x(g)),4*x(g)*(1-x(g)),x(g)*(2*x(g)-1)]
                deriv = [4*x(g)-3,4-8*x(g),4*x(g)-1]/h
                u = 0; du = 0; v = 0
                u(:3) = phi; du(:3) = deriv; v(4:) = phi
                if (e == 1) then
                    v(1) = phi(1); v(4) = 0
                end if
                bt = model%twist*r; dbt = model%twist
                bz = sqrt(model%bz_axis**2+model%shear*r*r)
                dbz = model%shear*r/bz
                ! mu=rho=1. Radial force balance p'=-Bt(Bt'+Bt/r)-Bz Bz'.
                jp = -dbz; jz = dbt+bt/r
                pp = -bt*jz-bz*dbz
                f = bt/r+model%k*bz
                ! Unit perpendicular direction: no ill-conditioned Btheta/Bz
                ! coefficient norm in the axial-field -> zero limit.
                v = bz/sqrt(bz*bz+bt*bt)*v
                w = -bt/bz*v
                div = du+u/r-v/r-model%k*w
                qr = f*u ! Actual Q_r=i*qr; other components are real.
                do a = 1, 6
                    helicity = v(a)*bz-w(a)*bt
                    qtheta(a) = -model%k*helicity-bt*du(a)-dbt*u(a)
                    qz(a) = helicity/r-bz*du(a)-(bz/r+dbz)*u(a)
                end do
                factor = wt(g)*h*r
                do b = 1, 6
                    col = dof(b)
                    if (col == 0) cycle
                    do a = 1, 6
                        row = dof(a)
                        if (row == 0) cycle
                        ! |curl(xi x B)|^2 - xi*.(J x Q) + p' u* div(xi).
                        integrand = qr(a)*qr(b)+qtheta(a)*qtheta(b)+qz(a)*qz(b) &
                            -u(a)*(jp*qz(b)-jz*qtheta(b)) &
                            -qr(b)*(jz*v(a)-jp*w(a))+pp*div(a)*u(b)
                        k(row,col) = k(row,col)+factor*integrand
                        mass(row,col) = mass(row,col)+factor* &
                            (u(a)*u(b)+v(a)*v(b)+w(a)*w(b))
                    end do
                end do
            end do
        end do
        if (maxval(abs(k-transpose(k))) > 1e-11_dp*maxval(abs(k))) &
            error stop 'independent cylinder force form is not symmetric'
        call dsygv(1,'V','U',nd,k,nd,mass,nd,eigen,work,size(work),stat)
        if (stat /= 0) error stop 'independent cylinder LAPACK eigenproblem'
        model%negative = count(eigen < -1e-9_dp)
        model%omega2 = eigen(1)
        if (allocated(model%u)) deallocate(model%u,model%v)
        allocate(model%u(2*cells+1),model%v(3,cells))
        model%u = 0; model%u(:2*cells) = k(:2*cells,1)
        do e = 1, cells
            do a = 1, 3
                i = 2*cells+3*(e-1)+a-1
                if (i == 2*cells) then
                    model%v(a,e) = model%u(1)
                else
                    model%v(a,e) = k(i,1)
                end if
            end do
        end do
        norm = maxval(abs(model%u))
        model%u = model%u/norm; model%v = model%v/norm
    end subroutine cylinder_solve

    subroutine cylinder_profile(model, r, xi)
        type(cylinder_oracle_t), intent(in) :: model
        real(dp), intent(in) :: r
        real(dp), intent(out) :: xi(3)
        real(dp) :: t, h, phi(3), bz
        integer :: e
        h = model%radius/model%cells
        e = min(model%cells,int(r/h)+1)
        t = r/h-(e-1)
        phi = [(1-t)*(1-2*t),4*t*(1-t),t*(2*t-1)]
        xi(1) = dot_product(phi,model%u(2*e-1:2*e+1))
        bz = sqrt(model%bz_axis**2+model%shear*r*r)
        xi(2) = bz/sqrt(bz*bz+(model%twist*r)**2)* &
            dot_product(phi,model%v(:,e))
        xi(3) = -model%twist*r/bz*xi(2)
    end subroutine cylinder_profile
end module s0_cylinder_oracle
