program test_s0_tae
    use, intrinsic :: iso_fortran_env, only: dp => real64, sp => real32
    use kin6d_s0_leaves, only: s0_tae_block_dp, s0_tae_block_sp, &
        s0_tae_block_jvp_dp, s0_tae_block_vjp_dp, s0_tae_edges_dp
    use kin6d_s0_spectrum, only: s0_pencil_t, s0_report_t, s0_prepare, &
        s0_solve, s0_destroy
    use fortsparse, only: csc_z_t, csc_from_triplet, csc_matvec, &
        fortsparse_status_t, FORTSPARSE_OK, FORTSPARSE_INERTIA_RELIABLE
    use fortnum_eigen, only: eigen_options_t, EIGEN_OK, EIGEN_NEAREST
    implicit none
    real(dp) :: errors(2,3), order(2)
    integer :: level
    call derivative_contract()
    do level = 1, 3
        call gap(8*2**level,errors(:,level))
    end do
    order = log(errors(:,1)/errors(:,2))/log(2.0_dp)
    call require(all(order > 0.95_dp) .and. all(order < 1.05_dp), &
        'first-order P0 continuum-edge convergence')
    order = log(errors(:,2)/errors(:,3))/log(2.0_dp)
    call require(all(order > 0.95_dp) .and. all(order < 1.05_dp), &
        'asymptotic P0 continuum-edge convergence')
    print '(a,2f10.5)', 'P0 omega edge convergence orders: ',order
contains
    subroutine require(ok,label)
        logical, intent(in) :: ok
        character(*), intent(in) :: label
        if (ok) return
        print *, 'FAIL: '//label
        error stop 1
    end subroutine require

    subroutine block(x,y)
        real(dp), intent(in) :: x(5)
        real(dp), intent(out) :: y(4)
        call s0_tae_block_dp(x(1),x(2),x(3),x(4),x(5), &
            y(1),y(2),y(3),y(4))
    end subroutine block

    subroutine derivative_contract()
        real(dp) :: x(5), dx(5), y(4), dy(4), g(5), seed(4), &
            plus(4), minus(4), step
        real(sp) :: ys(4), xs(5)
        x = [0.1_dp,-0.2_dp,0.08_dp,1.3_dp,2.1_dp]
        dx = [0.03_dp,-0.02_dp,0.01_dp,0.04_dp,0.02_dp]
        seed = [0.2_dp,-0.3_dp,0.5_dp,0.7_dp]
        call block(x,y)
        call s0_tae_block_jvp_dp(x(1),x(2),x(3),x(4),x(5), &
            dx(1),dx(2),dx(3),dx(4),dx(5),dy(1),dy(2),dy(3),dy(4))
        call s0_tae_block_vjp_dp(x(1),x(2),x(3),x(4),x(5), &
            seed(1),seed(2),seed(3),seed(4),g(1),g(2),g(3),g(4),g(5))
        call require(abs(dot_product(seed,dy)-dot_product(g,dx)) < &
            128*epsilon(1.0_dp),'toroidal block JVP/VJP pairing')
        step = epsilon(1.0_dp)**(1.0_dp/3)
        call block(x+step*dx,plus); call block(x-step*dx,minus)
        call require(maxval(abs((plus-minus)/(2*step)-dy)) < &
            8*step**2,'toroidal block independent finite differences')
        xs = real(x,sp)
        call s0_tae_block_sp(xs(1),xs(2),xs(3),xs(4),xs(5), &
            ys(1),ys(2),ys(3),ys(4))
        call require(maxval(abs(real(ys,dp)-y)) < &
            8*epsilon(1.0_sp)*maxval(abs(y)),'kind-generic block')
    end subroutine derivative_contract

    subroutine gap(ne,error)
        integer, intent(in) :: ne
        real(dp), intent(out) :: error(2)
        ! Flat q=3/2, m=1,2 and n=1 in a thin circular annulus.
        ! The retained cold field-line equation has mass weight
        ! 1+2*r/R*cos(theta). No radial derivative: these are continuum
        ! cells. Their two branches bound a gap, not a global TAE mode.
        integer, parameter :: angular=16
        real(dp), parameter :: major=10, density=1.3_dp, speed2=2.1_dp
        real(dp), parameter :: q=1.5_dp, inner=0.5_dp, outer=1.5_dp
        integer :: rows(4*ne), cols(4*ne), e, i, j, t, at, branch, n
        real(dp) :: h, radius, eps, center, lower, upper, edge(2), &
            cell(2), rawk(2,2), rawm(2,2), theta, kpar(2), &
            coeff(4,ne), res(1), tolerance, allowance, norm_mass, measured(2)
        complex(dp) :: kv(4*ne), mv(4*ne), phase(2), deriv(2), &
            start(2*ne), vectors(2*ne,1), lambda(1), shift, mx(2*ne), kx(2*ne)
        type(csc_z_t) :: stiffness,mass
        type(fortsparse_status_t) :: status
        type(s0_pencil_t) :: pencil
        type(s0_report_t) :: report
        type(eigen_options_t) :: options
        n = 2*ne; h = (outer-inner)/ne
        kpar = [(1-1/q)/major,(1-2/q)/major]
        center = speed2/(2*q*major)**2
        call s0_tae_edges_dp(center,inner/major,lower,upper)
        edge = [lower,upper]
        call s0_tae_edges_dp(center,(inner+h/2)/major,cell(1),cell(2))
        ! Independent angular Fourier quadrature of the defining functional.
        ! Structural zeros are retained so K and M have identical CSC pattern.
        !$omp parallel do default(none) shared(ne,h,kpar,rows,cols,kv,mv,coeff) &
        !$omp private(e,radius,eps,rawk,rawm,t,theta,phase,deriv,i,j,at)
        do e = 1, ne
            radius = inner+(e-0.5_dp)*h; eps = radius/major
            rawk = 0; rawm = 0
            do t = 0, angular-1
                theta = 2*acos(-1.0_dp)*t/angular
                phase = [exp(cmplx(0.0_dp,theta,dp)), &
                    exp(cmplx(0.0_dp,2*theta,dp))]
                deriv = cmplx(0.0_dp,1.0_dp,dp)*kpar*phase
                do j = 1, 2
                    do i = 1, 2
                        rawk(i,j) = rawk(i,j)+density*speed2 &
                            *real(conjg(deriv(i))*deriv(j),dp)/angular
                        rawm(i,j) = rawm(i,j)+density*(1+2*eps*cos(theta)) &
                            *real(conjg(phase(i))*phase(j),dp)/angular
                    end do
                end do
            end do
            call s0_tae_block_dp(kpar(1),kpar(2),eps,density,speed2, &
                coeff(1,e),coeff(2,e),coeff(3,e),coeff(4,e))
            do j = 1, 2
                do i = 1, 2
                    at = 4*(e-1)+2*(j-1)+i
                    rows(at) = 2*(e-1)+i; cols(at) = 2*(e-1)+j
                    kv(at) = cmplx(h*rawk(i,j),0.0_dp,dp)
                    mv(at) = cmplx(h*rawm(i,j),0.0_dp,dp)
                end do
            end do
        end do
        !$omp end parallel do
        do e = 1, ne
            at = 4*(e-1)
            call require(abs(kv(at+1)/h-coeff(1,e)) < &
                64*epsilon(h)*coeff(1,e),'Fourier line bending left')
            call require(abs(kv(at+4)/h-coeff(2,e)) < &
                64*epsilon(h)*coeff(2,e),'Fourier line bending right')
            call require(abs(mv(at+2)/h-coeff(4,e)) < &
                64*epsilon(h)*coeff(3,e),'Fourier toroidal coupling')
            call require(abs(mv(at+1)/h-coeff(3,e)) < &
                64*epsilon(h)*coeff(3,e),'Fourier physical inertia')
        end do
        call csc_from_triplet(n,n,rows,cols,kv,stiffness,status)
        call require(status%code == FORTSPARSE_OK,'toroidal stiffness')
        call csc_from_triplet(n,n,rows,cols,mv,mass,status)
        call require(status%code == FORTSPARSE_OK,'toroidal mass')
        do i = 1, n
            start(i) = cmplx(sin(real(i,dp)),cos(sqrt(2.0_dp)*i),dp)
        end do
        options%which = EIGEN_NEAREST; options%max_basis = min(n,48)
        options%max_restarts = 80; tolerance = 1.0e-11_dp
        do branch = 1, 2
            ! Both shifts lie strictly inside the continuous gap.
            shift = cmplx(0.8_dp*edge(branch)+0.2_dp*edge(3-branch),0.0_dp,dp)
            call s0_prepare(pencil,stiffness,mass,shift,report)
            call require(report%status == EIGEN_OK,'gap shift-invert factor')
            call require(report%inertia_status == FORTSPARSE_INERTIA_RELIABLE &
                .and. report%below == ne,'no continuum levels inside gap')
            call s0_solve(pencil,start,tolerance,options,lambda,vectors,res,report)
            call require(report%status == EIGEN_OK,'gap-targeted continuum edge')
            mx = csc_matvec(mass,vectors(:,1))
            kx = csc_matvec(stiffness,vectors(:,1))
            norm_mass = real(dot_product(vectors(:,1),mx),dp)
            call require(abs(norm_mass-1) < 1.0e-9_dp,'physical normalization')
            ! M >= density*h*(1-outer/major) I; use the original residual
            ! to bound distance to the finite-pencil spectrum, plus rounding.
            allowance = res(1)/sqrt(density*h*(1-outer/major)) &
                +128*epsilon(h)*center
            call require(abs(lambda(1)-cell(branch)) < allowance, &
                'exact discrete toroidal gap edge')
            call require(sqrt(sum(abs(kx-lambda(1)*mx)**2)) <= &
                res(1)+128*epsilon(h)*sqrt(sum(abs(kx)**2)), &
                'original equation residual')
            call require(report%claim == 'candidate' .and. &
                report%bound_status == 'unavailable','spectral status')
            measured(branch) = sqrt(real(lambda(1),dp))
            error(branch) = abs(measured(branch)-sqrt(edge(branch)))/sqrt(edge(branch))
            call s0_destroy(pencil)
        end do
        call require(maxval(error) < 0.003_dp,'physical edge frequency accuracy')
        print '(a,i0,a,2es17.8,a,2es12.4)', 'TAE continuum ne=',ne, &
            ' omega edges=',measured,' relative errors=',error
    end subroutine gap
end program test_s0_tae
