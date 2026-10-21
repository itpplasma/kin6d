program test_s0
    use, intrinsic :: iso_fortran_env, only: dp => real64, sp => real32, int64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
    use kin6d_s0_volume, only: s0_assemble_volume
    use kin6d_ideal_xi, only: xi_gradient_n, XI_OK
    use kin6d_s0_spectrum, only: s0_pencil_t, s0_report_t, s0_prepare, &
        s0_solve, s0_destroy
    use kin6d_s0_leaves, only: s0_mass_dp, s0_mass_jvp_dp, s0_mass_vjp_dp, &
        s0_mass_sp, s0_alfven_dp, s0_alfven_sp, s0_kink_dp
    use fortsparse, only: csc_z_t, fortsparse_status_t, FORTSPARSE_OK, &
        FORTSPARSE_INERTIA_RELIABLE, FORTSPARSE_INERTIA_UNAVAILABLE, &
        csc_matvec, csc_from_triplet
    use fortnum_eigen, only: eigen_options_t, EIGEN_OK, EIGEN_INVALID, EIGEN_NEAREST
    implicit none
    character(32) :: mode, arg
    real(dp) :: errors(3), ignored
    integer :: i, ne, repeats
    call get_command_argument(1,mode)
    select case(trim(mode))
    case('continuum')
        call continuum()
        call cylinder_continuum()
    case('contracts')
        call contracts()
    case('timing')
        call get_command_argument(2,arg); read(arg,*) ne
        call get_command_argument(3,arg); read(arg,*) repeats
        do i = 1, repeats
            call wave(ne,.false.,ignored)
        end do
    case default
        do i = 1, 3
            call wave(8*2**i,.false.,errors(i))
        end do
        call expect(all(errors(:2)/errors(2:) > 3.8_dp),'observed P1 eigenvalue order')
        call expect(all(errors(:2)/errors(2:) < 4.2_dp),'observed P1 eigenvalue order')
        call wave(32,.true.,ignored)
        print '(a,3es14.5)', 'P1 physical relative errors: ',errors
    end select
contains
    subroutine expect(condition,label)
        logical, intent(in) :: condition
        character(*), intent(in) :: label
        if (.not. condition) then
            print '(a)', 'FAIL: '//label
            error stop 1
        end if
    end subroutine expect

    subroutine wave(ne,complex_shift,error)
        integer, intent(in) :: ne
        logical, intent(in) :: complex_shift
        real(dp), intent(out) :: error
        integer, allocatable :: dof(:,:)
        real(dp), allocatable :: w(:,:), b(:,:,:), p(:,:), rho(:,:)
        complex(dp), allocatable :: value(:,:,:,:), grad(:,:,:,:,:), start(:), &
            vectors(:,:), kx(:), mx(:)
        type(csc_z_t) :: stiffness, mass
        type(fortsparse_status_t) :: status
        type(s0_pencil_t) :: pencil
        type(s0_report_t) :: report
        type(eigen_options_t) :: options
        real(dp) :: h, local, pi, speed2, exact, discrete, mu, density, res(2), &
            tolerance, allowance, a_time, f_time, e_time, begin, end
        complex(dp) :: lambda(2), shift
        integer :: e, q, i, j, n
        integer(int64) :: rate, tick
        n = ne-1; h = 1.0_dp/ne; pi = acos(-1.0_dp)
        density = 1.3_dp; mu = 1.7_dp; speed2 = 4/(mu*density)
        allocate(dof(2,ne),w(2,ne),b(3,2,ne),p(2,ne),rho(2,ne), &
            value(3,2,2,ne),grad(3,3,2,2,ne),start(n),vectors(n,2),kx(n),mx(n))
        do e = 1, ne
            dof(:,e) = [e-1,e]
        end do
        dof(2,ne) = 0
        w = h/2; b = 0; b(3,:,:) = 2; p = 0.7_dp; rho = density
        value = 0; grad = 0
        do e = 1, ne
            do q = 1, 2
                local = (1+(2*q-3)/sqrt(3.0_dp))/2
                value(1,:,q,e) = cmplx([1-local,local],0.0_dp,dp)
                grad(1,3,:,q,e) = cmplx([-1/h,1/h],0.0_dp,dp)
            end do
        end do
        call system_clock(tick,rate); begin = real(tick,dp)/rate
        call s0_assemble_volume(n,dof,w,b,p,rho,mu,value,grad,stiffness,mass,status)
        call system_clock(tick); end = real(tick,dp)/rate; a_time = end-begin
        call expect(status%code == FORTSPARSE_OK,'wave assembly')
        exact = pi*pi*speed2
        shift = cmplx(1.5_dp*exact,0.0_dp,dp)
        if (complex_shift) shift = shift+cmplx(0.0_dp,0.4_dp*exact,dp)
        begin = end
        call s0_prepare(pencil,stiffness,mass,shift,report)
        call system_clock(tick); end = real(tick,dp)/rate; f_time = end-begin
        if (report%status /= EIGEN_OK) print *, report%sparse_status,report%message
        call expect(report%status == EIGEN_OK,'wave factorization')
        if (complex_shift) then
            call expect(report%inertia_status == FORTSPARSE_INERTIA_UNAVAILABLE, &
                'complex shifts have no Sylvester count')
        else
            call expect(report%inertia_status == FORTSPARSE_INERTIA_RELIABLE, &
                'real-shift inertia reliability')
            call expect(report%below == 1 .and. report%at == 0 .and. &
                report%above == n-1,'physical generalized inertia')
        end if
        do i = 1, n
            start(i) = cmplx(sin(real(i,dp)*sqrt(2.0_dp)), &
                cos(real(i,dp)*sqrt(3.0_dp)),dp)
        end do
        options%max_basis = min(n,40); options%max_restarts = 40
        options%which = EIGEN_NEAREST; tolerance = 1.0e-10_dp
        begin = end
        call s0_solve(pencil,start,tolerance,options,lambda,vectors,res,report)
        call system_clock(tick); end = real(tick,dp)/rate; e_time = end-begin
        if (report%status /= EIGEN_OK) print *, 'eigen status/residuals',report%status,res
        call expect(report%status == EIGEN_OK,'frequency-targeted eigen solve')
        call expect(report%claim == 'candidate' .and. &
            report%bound_status == 'unavailable','spectral claim')
        do j = 1, 2
            discrete = 6*speed2/h**2*(1-cos(j*pi/ne))/(2+cos(j*pi/ne))
            kx = csc_matvec(stiffness,vectors(:,j))
            mx = csc_matvec(mass,vectors(:,j))
            call expect(abs(dot_product(vectors(:,j),mx)-1) < 1.0e-9_dp, &
                'physical inertia normalization')
            ! M >= rho*h/3 I. The residual gives distance to the Hermitian
            ! M^-1/2 K M^-1/2 spectrum: ||r||/sqrt(lambda_min(M)), ||x||M=1.
            allowance = res(j)/sqrt(density*h/3)+ &
                128*epsilon(h)*maxval(abs(stiffness%val))/(density*h/3)
            call expect(abs(lambda(j)-discrete) <= allowance,'exact FE sine spectrum')
            call expect(sqrt(sum(abs(kx-lambda(j)*mx)**2)) <= &
                res(j)*(1+1.0e-5_dp)+128*epsilon(h)*maxval(abs(kx)), &
                'original-pencil residual')
        end do
        error = abs(real(lambda(1))-exact)/exact
        call expect(error < 0.01_dp,'one percent physical eigenvalue accuracy')
        call s0_destroy(pencil)
        print '(a,i0,a,l1,a,3es14.5,a,2es17.8)', &
            'wave ne=',ne,' complex=',complex_shift,' assembly/factor/solve s=', &
            a_time,f_time,e_time,' omega2=',real(lambda,dp)
    end subroutine wave

    subroutine continuum()
        ! Straight slab analogue of the cylinder shear continuum: B=Bz(r) ez,
        ! p+Bz^2/(2mu)=constant, xi=ey u(r) exp(ikz), rho constant.
        ! No transverse derivative energy: P0 radial cells discretize the
        ! multiplication operator. They do not create discrete global modes.
        integer, parameter :: ne=32
        integer :: dof(1,ne), e, q
        real(dp) :: w(2,ne), b(3,2,ne), p(2,ne), rho(2,ne), &
            h, r, k, mu, exact, omega2, omega2s, res(2)
        real(sp) :: single
        complex(dp) :: value(3,1,2,ne), grad(3,3,1,2,ne), &
            start(ne), vectors(ne,2), lambda(2), shift
        type(csc_z_t) :: stiffness,mass
        type(fortsparse_status_t) :: status
        type(s0_pencil_t) :: pencil
        type(s0_report_t) :: report
        type(eigen_options_t) :: options
        h = 1.0_dp/ne; k = 0.7_dp; mu = 1.2_dp
        w = h/2; b = 0; value = 0; value(2,1,:,:) = 1
        grad = 0; grad(2,3,1,:,:) = cmplx(0.0_dp,k,dp); rho = 1.4_dp
        do e = 1, ne
            dof(1,e) = e
            start(e) = cmplx(1.0_dp,0.1_dp*e,dp)
            do q = 1, 2
                r = (e-1+(1+(2*q-3)/sqrt(3.0_dp))/2)*h
                b(3,q,e) = 1+r
                p(q,e) = 4-b(3,q,e)**2/(2*mu)
            end do
        end do
        call s0_assemble_volume(ne,dof,w,b,p,rho,mu,value,grad,stiffness,mass,status)
        call expect(status%code == FORTSPARSE_OK,'continuum assembly')
        do e = 1, ne
            r = (e-0.5_dp)*h
            call s0_alfven_dp(k,1+r,mu,rho(1,e),omega2)
            call s0_alfven_sp(real(k,sp),real(1+r,sp),real(mu,sp), &
                real(rho(1,e),sp),single)
            omega2s = real(single,dp)
            call expect(abs(omega2s-omega2) < 16*epsilon(single)*omega2,'FP32 oracle')
            exact = omega2+k*k*h*h/(12*mu*rho(1,e))
            call expect(abs(stiffness%val(e)/mass%val(e)-exact) < &
                64*epsilon(h)*exact,'analytic cell-average Alfven continuum')
        end do
        ! A shift halfway between cells 12 and 13 has exactly 12 lower cells.
        shift = (stiffness%val(12)/mass%val(12)+ &
            stiffness%val(13)/mass%val(13))/2
        call s0_prepare(pencil,stiffness,mass,shift,report)
        if (report%status /= EIGEN_OK) print *, report%sparse_status,report%message
        call expect(report%status == EIGEN_OK,'continuum factor')
        call expect(report%inertia_status == FORTSPARSE_INERTIA_RELIABLE &
            .and. report%below == 12,'continuum counting below shift')
        options%max_basis = 24; options%which = EIGEN_NEAREST
        call s0_solve(pencil,start,1.0e-10_dp,options,lambda,vectors,res,report)
        call expect(report%status == EIGEN_OK,'continuum frequency targeting')
        call expect(maxval(abs(aimag(lambda))) < 1.0e-12_dp,'real ideal spectrum')
        do e = 1, 2
            call expect(minval(abs(stiffness%val/mass%val-lambda(e))) < &
                1.0e-8_dp,'continuum eigenvalue belongs to cell spectrum')
        end do
        call s0_destroy(pencil)
        print '(a,2es17.8)', 'shear-continuum targeted omega2=',real(lambda,dp)
    end subroutine continuum

    subroutine cylinder_continuum()
        ! Torsional n=0 modes on 1<r<2, Bz=1+r, rho constant,
        ! p+Bz^2/(2mu)=C. The cylindrical cross term integrates to
        ! -2*pi*C*[eta*xi] on the radial boundary and vanishes here.
        integer, parameter :: ne=32, n=ne-1
        integer :: dof(2,ne), e, q, i, j, at, row, col, ierr
        real(dp) :: w(3,ne), b(3,3,ne), p(3,ne), rho(3,ne), &
            h, r, k, mu, pi, phi(2), derivative(2), &
            nodes(3), weights(3), reference(n,n), mass_ref(n,n), res(2), &
            lower, upper, tolerance
        complex(dp) :: value(3,2,3,ne), grad(3,3,2,3,ne), &
            xi(3), dr(3), dz(3), start(n), vectors(n,2), lambda(2), shift
        type(csc_z_t) :: stiffness,mass
        type(fortsparse_status_t) :: status
        type(s0_pencil_t) :: pencil
        type(s0_report_t) :: report
        type(eigen_options_t) :: options
        h = 1.0_dp/ne; k = 0.7_dp; mu = 1.2_dp; pi = acos(-1.0_dp)
        nodes = [-sqrt(0.6_dp),0.0_dp,sqrt(0.6_dp)]
        weights = [5.0_dp/9,8.0_dp/9,5.0_dp/9]
        b = 0; rho = 1.4_dp; value = 0; grad = 0
        reference = 0; mass_ref = 0; derivative = [-1/h,1/h]
        do e = 1, ne
            dof(:,e) = [e-1,e]
        end do
        dof(2,ne) = 0
        do e = 1, ne
            do q = 1, 3
                phi = [(1-nodes(q))/2,(1+nodes(q))/2]
                r = 1+(e-1+phi(2))*h
                w(q,e) = weights(q)*h*pi*r
                b(3,q,e) = 1+r; p(q,e) = 5-b(3,q,e)**2/(2*mu)
                do i = 1, 2
                    xi = 0; xi(2) = phi(i)
                    dr = 0; dr(2) = derivative(i)
                    dz = 0; dz(2) = cmplx(0.0_dp,k*phi(i),dp)
                    value(:,i,q,e) = xi
                    call xi_gradient_n(0,r,xi,dr,dz,grad(:,:,i,q,e),ierr)
                    call expect(ierr == XI_OK,'cylindrical physical jet')
                end do
                ! Independent reduced physical multiplication operator;
                ! polynomial degree five, integrated exactly by three points.
                do j = 1, 2
                    col = dof(j,e)
                    if (col == 0) cycle
                    do i = 1, 2
                        row = dof(i,e)
                        if (row == 0) cycle
                        reference(row,col) = reference(row,col)+ &
                            w(q,e)*k*k*(1+r)**2/mu*phi(i)*phi(j)
                        mass_ref(row,col) = mass_ref(row,col)+ &
                            w(q,e)*rho(q,e)*phi(i)*phi(j)
                    end do
                end do
            end do
        end do
        call s0_assemble_volume(n,dof,w,b,p,rho,mu,value,grad,stiffness,mass,status)
        call expect(status%code == FORTSPARSE_OK,'cylindrical continuum assembly')
        tolerance = 1024*epsilon(h)*maxval(abs(reference))
        do col = 1, n
            do at = stiffness%col_ptr(col),stiffness%col_ptr(col+1)-1
                row = stiffness%row_idx(at)
                call expect(abs(stiffness%val(at)-reference(row,col)) < tolerance, &
                    'cylindrical boundary cross-term cancellation')
                call expect(abs(mass%val(at)-mass_ref(row,col)) < tolerance, &
                    'cylindrical physical inertia')
            end do
            start(col) = cmplx(1.0_dp,0.1_dp*col,dp)
        end do
        call s0_alfven_dp(k,2.0_dp,mu,rho(1,1),lower)
        call s0_alfven_dp(k,3.0_dp,mu,rho(1,1),upper)
        call s0_alfven_dp(k,2.5_dp,mu,rho(1,1),r)
        shift = cmplx(r,0.0_dp,dp)
        call s0_prepare(pencil,stiffness,mass,shift,report)
        call expect(report%status == EIGEN_OK,'cylindrical continuum factor')
        call expect(report%inertia_status == FORTSPARSE_INERTIA_RELIABLE, &
            'cylindrical continuum inertia reliability')
        options%max_basis = n; options%which = EIGEN_NEAREST
        call s0_solve(pencil,start,1.0e-10_dp,options,lambda,vectors,res,report)
        call expect(report%status == EIGEN_OK,'cylindrical continuum targeting')
        call expect(all(real(lambda,dp) > lower) .and. all(real(lambda,dp) < upper), &
            'continuum spectral range')
        call expect(maxval(abs(real(lambda,dp)-r)) < 0.04_dp*(upper-lower), &
            'frequency targeting inside cylindrical continuum')
        call s0_destroy(pencil)
        print '(a,2es17.8)', 'cylindrical torsional-continuum omega2=',real(lambda,dp)
    end subroutine cylinder_continuum

    subroutine contracts()
        type(csc_z_t) :: stiffness,mass,bad
        type(fortsparse_status_t) :: status
        type(s0_pencil_t) :: pencil
        type(s0_report_t) :: report
        type(eigen_options_t) :: options
        real(dp) :: x(7), d(7), g(7), v, dv, gamma2, mass_value, res(2), exact(2)
        complex(dp) :: start(2), vectors(2,2), lambda(2)
        real(sp) :: vs
        integer :: i
        do i = 1, 7
            x(i) = 0.2_dp*i; d(i) = cos(real(i,dp))
        end do
        call s0_mass_dp(x(1),x(2),x(3),x(4),x(5),x(6),x(7),v)
        call s0_mass_sp(real(x(1),sp),real(x(2),sp),real(x(3),sp),real(x(4),sp), &
            real(x(5),sp),real(x(6),sp),real(x(7),sp),vs)
        call expect(abs(real(vs,dp)-v) <= 8*epsilon(vs)*abs(v),'FP32 mass')
        call s0_mass_jvp_dp(x(1),x(2),x(3),x(4),x(5),x(6),x(7), &
            d(1),d(2),d(3),d(4),d(5),d(6),d(7),dv)
        call s0_mass_vjp_dp(x(1),x(2),x(3),x(4),x(5),x(6),x(7),1.0_dp, &
            g(1),g(2),g(3),g(4),g(5),g(6),g(7))
        call expect(abs(dv-dot_product(g,d)) < 64*epsilon(v),'mass JVP/VJP')
        mass_value = x(1)*sum(x(2:4)*x(5:7))
        call expect(abs(v-mass_value) < 8*epsilon(v),'independent mass oracle')
        call s0_kink_dp(0.4_dp,1.7_dp,gamma2)
        call expect(abs(gamma2-0.48_dp/1.7_dp) < 8*epsilon(v),'rigid kink growth oracle')
        ! An independently specified 2x2 Hermitian pencil exercises complex
        ! phase, an unstable mode, non-diagonal mass, and a 2x2 LDLH pivot.
        call csc_from_triplet(2,2,[1,2,1,2],[1,1,2,2], &
            [cmplx(0,0,dp),cmplx(0,2,dp),cmplx(0,-2,dp),cmplx(0,0,dp)], &
            stiffness,status)
        call csc_from_triplet(2,2,[1,2,1,2],[1,1,2,2], &
            [cmplx(2,0,dp),cmplx(0,0.2_dp,dp),cmplx(0,-0.2_dp,dp), &
             cmplx(1,0,dp)],mass,status)
        call s0_prepare(pencil,stiffness,mass,cmplx(0,0,dp),report)
        if (report%status /= EIGEN_OK) print *, report%sparse_status,report%message
        call expect(report%status == EIGEN_OK,'complex Hermitian pencil')
        call expect(report%below == 1 .and. report%above == 1,'indefinite inertia')
        options%max_basis = 2; options%which = EIGEN_NEAREST
        start = [cmplx(1,0,dp),cmplx(0.3_dp,0.7_dp,dp)]
        call s0_solve(pencil,start,1.0e-10_dp,options,lambda,vectors,res,report)
        call expect(report%status == EIGEN_OK,'unstable complex Hermitian solve')
        ! det(K-lambda M) = 1.96 lambda^2 + 0.8 lambda - 4.
        exact = [(-0.8_dp-sqrt(32.0_dp))/3.92_dp, &
                 (-0.8_dp+sqrt(32.0_dp))/3.92_dp]
        do i = 1, 2
            call expect(minval(abs(lambda(i)-exact)) < &
                res(i)/sqrt(0.8_dp)+128*epsilon(v),'analytic indefinite pencil roots')
        end do
        bad = mass; bad%val = -bad%val
        call s0_prepare(pencil,stiffness,bad,cmplx(0,0,dp),report)
        call expect(report%status /= EIGEN_OK,'negative physical metric rejected')
        bad = stiffness; bad%val(2) = 1
        call s0_prepare(pencil,bad,mass,cmplx(0,0,dp),report)
        call expect(report%status == EIGEN_INVALID,'non-Hermitian K rejected')
        call s0_prepare(pencil,stiffness,mass, &
            cmplx(ieee_value(1.0_dp,ieee_quiet_nan),0.0_dp,dp),report)
        call expect(report%status == EIGEN_INVALID,'NaN shift rejected')
        call s0_destroy(pencil)
        print '(a)', 'S0 physical metric, phase, inertia and derivative contracts PASS'
    end subroutine contracts
end program test_s0
