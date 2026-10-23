module test_s0_gs_kink_callbacks
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_grad_shafranov, only: gs_state_t, initialize_gs, solve_gs
    use kin6d_gs_equilibrium, only: gs_equilibrium
    use kin6d_equilibrium_state, only: equilibrium_state_t
    use kin6d_s0_gs_space, only: s0_gs_space_t, s0_gs_tables, s0_gs_assemble
    use kin6d_s0_spectrum, only: s0_pencil_t, s0_report_t, s0_prepare, &
        s0_solve, s0_destroy
    use fortsparse, only: csc_z_t, fortsparse_status_t, FORTSPARSE_OK, &
        FORTSPARSE_INERTIA_RELIABLE, csc_matvec
    use fortnum_eigen, only: eigen_options_t, EIGEN_OK, EIGEN_NEAREST
    use s0_cylinder_oracle, only: cylinder_oracle_t, cylinder_solve, cylinder_profile

    use kin6d_curved_triangle, only: boundary_curve_t
    use kin6d_grad_shafranov, only: initialize_gs_object
    implicit none
    private
    public :: run_test_s0_gs_kink
    real(dp), parameter :: radius = 2.0_dp
    type :: run_test_s0_gs_kink_context_t
        real(dp), pointer :: r0 => null()
    end type

    type, extends(boundary_curve_t) :: run_test_s0_gs_kink_boundary_t
        type(run_test_s0_gs_kink_context_t), pointer :: data => null()
    contains
        procedure :: evaluate => run_test_s0_gs_kink_boundary_binding
    end type

contains
    subroutine run_test_s0_gs_kink()
    implicit none
    type(gs_state_t) :: gs
    type(equilibrium_state_t) :: eq
    type(cylinder_oracle_t) :: oracle, analytic
    real(dp) , target :: r0, eigenvalue(4), shape_error(4), reference, rates(2), &
        errors(3)
    real(dp) :: radial_order, radial_estimate, aspect_estimate, quadrature_value
    integer :: i, stat
    character(32) :: mode
        type(run_test_s0_gs_kink_context_t), target :: ctx
        type(run_test_s0_gs_kink_boundary_t), target :: boundary_operator
    call get_command_argument(1,mode)
    select case (trim(mode))
    case ('oracle')
        analytic%twist = 0; analytic%shear = 0
        call cylinder_solve(analytic,16)
        reference = analytic%k**2*analytic%bz_axis**2
        call expect(abs(analytic%omega2/reference-1) < 2e-10_dp, &
            'analytic axial-field Alfven eigenvalue')
        call expect(analytic%negative == 0,'analytic axial-field inertia')
        analytic%twist = 1; analytic%shear = 0; analytic%bz_axis = 1e-8_dp
        do i = 1, 3
            analytic%k = -20.0_dp*2**i
            call cylinder_solve(analytic,64)
            errors(i) = abs(analytic%omega2+3)
            print '(a,f8.1,a,es18.10)', '[V] z-pinch k=',analytic%k, &
                ' omega2 -> -3: ',analytic%omega2
        end do
        rates = log(errors(:2)/errors(2:))/log(2.0_dp)
        call expect(all(rates > 1.8_dp) .and. all(rates < 2.2_dp), &
            'analytic unstable z-pinch large-k limit order')
        call expect(abs(sqrt(-analytic%omega2/3)-1) < 2e-4_dp, &
            'analytic gamma -> sqrt(3) unstable limit')
        do i = 1, 3
            call cylinder_solve(oracle,16*2**i)
            eigenvalue(i) = oracle%omega2
            call expect(oracle%negative == 1,'independent cylinder unstable count')
        end do
        call expect(abs(eigenvalue(3)-eigenvalue(2)) < &
            abs(eigenvalue(2)-eigenvalue(1))/8,'independent P2 eigenvalue refinement')
        call expect(abs(eigenvalue(3)/eigenvalue(2)-1) < 2e-6_dp, &
            'independent cylinder reference budget')
        print '(a,es18.10)', '[V] exact axial-field omega2=',reference
        print '(a,3es18.10)', '[V] independent cylinder P2 omega2=',eigenvalue(:3)
    case ('aspect')
        call cylinder_solve(oracle,128)
        do i = 1, 4
            r0 = 10.0_dp*2**i
            call equilibrium(r0,0.015_dp,eq)
            call spectrum(eq,r0,512,oracle,eigenvalue(i),shape_error(i))
        end do
        rates = log(abs((eigenvalue(:2)-eigenvalue(2:3))/ &
            (eigenvalue(2:3)-eigenvalue(3:4))))/log(2.0_dp)
        print '(a,2f12.6)', '[V] measured aspect-ratio orders=',rates
        call expect(all(rates > 1.8_dp) .and. all(rates < 2.2_dp), &
            'measured second-order cylinder limit')
        aspect_estimate = abs(eigenvalue(4)-eigenvalue(3))/ &
            ((2**rates(2)-1)*2*abs(oracle%omega2))
        call expect(1.25_dp*aspect_estimate < 0.0018_dp, &
            'separate empirical aspect-ratio gamma budget')
        print '(a,es16.7)', '[V] finest relative gamma error=', &
            abs(sqrt(eigenvalue(4)/oracle%omega2)-1)
        print '(a,es16.7)', '[V] aspect gamma estimate with 25% margin=', &
            1.25_dp*aspect_estimate
    case ('equilibrium')
        call cylinder_solve(oracle,128)
        r0 = 160
        do i = 1, 3
            call equilibrium(r0,0.06_dp/4**(i-1),eq)
            call spectrum(eq,r0,512,oracle,eigenvalue(i),shape_error(i))
        end do
        print '(a,3es18.10)', '[V] GS h-refinement omega2=',eigenvalue(:3)
        call expect(abs(eigenvalue(3)-eigenvalue(2)) < &
            abs(eigenvalue(2)-eigenvalue(1)),'GS refinement decreases change')
        call expect(abs(eigenvalue(3)-eigenvalue(2))/(2*abs(oracle%omega2)) &
            < 2e-6_dp,'separate empirical equilibrium gamma budget')
        call spectrum(eq,r0,512,oracle,quadrature_value,reference,angular_points=64)
        call expect(abs(quadrature_value-eigenvalue(3))/(2*abs(oracle%omega2)) &
            < 2e-6_dp,'angular quadrature gamma budget')
    case default
        call cylinder_solve(oracle,128)
        call expect(oracle%negative == 1,'independent cylinder unstable count')
        r0 = 160
        call equilibrium(r0,0.015_dp,eq)
        do i = 1, 3
            call spectrum(eq,r0,64*2**i,oracle,eigenvalue(i),shape_error(i))
        end do
        errors = abs(eigenvalue(:3)-oracle%omega2)
        radial_order = log(abs((eigenvalue(1)-eigenvalue(2))/ &
            (eigenvalue(2)-eigenvalue(3))))/log(2.0_dp)
        call expect(radial_order > 0.9_dp .and. radial_order < 1.2_dp, &
            'observed radial eigenvalue refinement order')
        call expect(all(errors(2:) < errors(:2)),'radial approach to cylinder value')
        radial_estimate = abs(eigenvalue(3)-eigenvalue(2))/ &
            ((2**radial_order-1)*2*abs(oracle%omega2))
        call expect(1.25_dp*radial_estimate < 0.006_dp, &
            'separate empirical radial gamma budget')
        call expect(abs(sqrt(eigenvalue(3)/oracle%omega2)-1) &
            < 0.006_dp+0.0018_dp+4e-6_dp+1e-7_dp,'total empirical gamma budget')
        call expect(all(shape_error(:2)/shape_error(2:3) > 1.3_dp), &
            'physical eigenvector converges to independent radial profile')
        call expect(shape_error(3) < 0.07_dp,'finest physical profile agreement')
        print '(a,f12.6)', '[V] measured radial successive-difference order=', &
            radial_order
        print '(a,3es16.7)', '[V] relative gamma errors=', &
            abs(sqrt(eigenvalue(:3)/oracle%omega2)-1)
        print '(a,3es16.7)', '[V] physical profile distances=',shape_error(:3)
        print '(a,es16.7)', '[V] radial gamma estimate with 25% margin=', &
            1.25_dp*radial_estimate
        print '(a,es16.7)', '[V] cylinder gamma*tau_A (tau_A=1)=', &
            sqrt(-oracle%omega2)
    end select
contains
    subroutine bind_run_test_s0_gs_kink_context()
        ctx%r0 => r0
        boundary_operator%data => ctx
    end subroutine bind_run_test_s0_gs_kink_context

    subroutine expect(condition,label)
        logical, intent(in) :: condition
        character(*), intent(in) :: label
        if (.not. condition) then
            print '(a)', 'FAIL: '//label
            error stop 1
        end if
    end subroutine expect

    subroutine equilibrium(major_radius, area, state)
        real(dp), intent(in) :: major_radius, area
        type(equilibrium_state_t), intent(out) :: state
        real(dp) :: contour(2,128)
        integer :: j
        r0 = major_radius
        do j = 1, size(contour,2)
        call bind_run_test_s0_gs_kink_context()
            call run_test_s0_gs_kink_boundary(ctx, &
                real(j-1,dp)/size(contour,2),contour(:,j))
        end do
        call bind_run_test_s0_gs_kink_context()
        call initialize_gs_object(gs,contour,area,stat,3,curve=boundary_operator)
        call expect(stat == 0,'GS curved P3 mesh')
        call solve_gs(gs,1.0_dp,-3/r0,r0,stat)
        call expect(stat == 0,'GS forward solve')
        call expect(gs%residual < 1e-10_dp,'GS weak equation residual')
        call expect(gs%iterations > 0,'GS solution from a zero seed')
        call gs_equilibrium(gs,state,stat,r0*sqrt(4.25_dp),r0, &
            pprime=-3/r0,mu0=1.0_dp)
        call expect(stat == 0,'GS core state adapter')
        print '(a,f8.1,a,i0,a,es12.4)', '[V] R0=',r0,' GS ndof=', &
            size(gs%psi),' residual=',gs%residual
    end subroutine equilibrium

    subroutine spectrum(state, major_radius, nr, model, omega2, distance, &
            angular_points)
        type(equilibrium_state_t), intent(in) :: state
        real(dp), intent(in) :: major_radius
        integer, intent(in) :: nr
        type(cylinder_oracle_t), intent(in) :: model
        real(dp), intent(out) :: omega2, distance
        integer, intent(in), optional :: angular_points
        type(s0_gs_space_t) :: space
        type(csc_z_t) :: stiffness, mass
        type(fortsparse_status_t) :: sparse
        type(s0_pencil_t) :: pencil
        type(s0_report_t) :: report
        type(eigen_options_t) :: options
        complex(dp), allocatable :: start(:), vectors(:,:), kx(:), mx(:)
        complex(dp) :: lambda(1), xi(3), ref(3), phase, overlap
        real(dp) :: residual(1), scaled_residual, knorm, xnorm, anorm, bnorm
        real(dp) :: r, theta, radial(3), nodes(4), h, c, s
        integer :: ndof, n, e, q, a, t, g, j, nt
        n = nint(model%k*major_radius)
        nt = 32
        if (present(angular_points)) nt = angular_points
        call s0_gs_tables(state,[major_radius,0.0_dp],radius,nr,nt,n, &
            1.0_dp,space,stat)
        call expect(stat == 0,'S0 GS sampler tables')
        call s0_gs_assemble(state,space,stiffness,mass,sparse)
        call expect(sparse%code == FORTSPARSE_OK,'S0 GS volume assembly')
        call s0_prepare(pencil,stiffness,mass,cmplx(0.0_dp,0.0_dp,dp),report)
        if (report%status /= EIGEN_OK) print *, report%message
        call expect(report%status == EIGEN_OK,'zero shift LDLH')
        call expect(report%inertia_status == FORTSPARSE_INERTIA_RELIABLE, &
            'LDLH inertia reliability')
        call expect(report%below == model%negative .and. report%at == 0, &
            'all unstable modes counted independently')
        call expect(report%above == space%ndof-model%negative,'positive mode count')
        call s0_prepare(pencil,stiffness,mass,cmplx(-0.2_dp,0.0_dp,dp),report)
        call expect(report%status == EIGEN_OK,'kink shift factor')
        ndof = space%ndof
        allocate(start(ndof),vectors(ndof,1))
        do j = 1, ndof
            start(j) = cmplx(sin(j*sqrt(2.0_dp)),cos(j*sqrt(3.0_dp)),dp)
        end do
        options%max_basis = min(60,ndof)
        options%max_restarts = 20; options%which = EIGEN_NEAREST
        call s0_solve(pencil,start,1e-10_dp,options,lambda,vectors,residual,report)
        call expect(report%status == EIGEN_OK,'physical kink eigenproblem')
        call expect(report%claim == 'candidate' .and. &
            report%bound_status == 'unavailable','candidate status retained')
        call expect(abs(aimag(lambda(1))) < 1e-13_dp,'real ideal eigenvalue')
        omega2 = real(lambda(1),dp)
        call expect(omega2 < 0,'ideal m=1 kink instability')
        kx = csc_matvec(stiffness,vectors(:,1))
        mx = csc_matvec(mass,vectors(:,1))
        knorm = 0
        do j = 1, ndof
            knorm = max(knorm,sum(abs(stiffness%val(stiffness%col_ptr(j): &
                stiffness%col_ptr(j+1)-1))))
        end do
        xnorm = sqrt(sum(abs(vectors(:,1))**2))
        scaled_residual = sqrt(sum(abs(kx-lambda(1)*mx)**2))/(knorm*xnorm)
        call expect(scaled_residual < 1e-11_dp,'original physical eigen residual')
        call expect(abs(dot_product(vectors(:,1),mx)-1) < 1e-9_dp, &
            'physical M normalization')
        ! Independent radial-profile projection, with the actual physical mass.
        ! Gauss nodes here are written independently of the adapter's leaves.
        nodes = [(1-sqrt((3+2*sqrt(1.2_dp))/7))/2, &
            (1-sqrt((3-2*sqrt(1.2_dp))/7))/2, &
            (1+sqrt((3-2*sqrt(1.2_dp))/7))/2, &
            (1+sqrt((3+2*sqrt(1.2_dp))/7))/2]
        h = radius/nr; anorm = 0; bnorm = 0; overlap = 0
        do e = 1, nr
            do t = 1, nt
                theta = 2*acos(-1.0_dp)*(t-1)/nt
                c = cos(theta); s = sin(theta); phase = cmplx(c,s,dp)
                do g = 1, 4
                    q = 4*(t-1)+g
                    r = h*(e-1+nodes(g))
                    call cylinder_profile(model,r,radial)
                    ref = phase*[cmplx(c*radial(1),-s*radial(2),dp), &
                        cmplx(0.0_dp,radial(3),dp), &
                        cmplx(s*radial(1),c*radial(2),dp)]
                    xi = 0
                    do a = 1, 6
                        j = space%dof(a,e)
                        if (j == 0) cycle
                        xi = xi+space%value(:,a,q,e)*vectors(j,1)
                    end do
                    overlap = overlap+space%weight(q,e)*dot_product(ref,xi)
                    anorm = anorm+space%weight(q,e)*sum(abs(ref)**2)
                    bnorm = bnorm+space%weight(q,e)*sum(abs(xi)**2)
                end do
            end do
        end do
        distance = sqrt(max(0.0_dp,2*(1-abs(overlap)/sqrt(anorm*bnorm))))
        print '(a,f7.1,a,i4,a,es18.10,a,es10.3,a,es10.3)', &
            '[V] R0=',major_radius,' Nr=',nr,' omega2=',omega2, &
            ' residual=',scaled_residual,' profile distance=',distance
        call s0_destroy(pencil)
    end subroutine spectrum


    end subroutine run_test_s0_gs_kink

    subroutine run_test_s0_gs_kink_boundary(ctx, t,point)
        type(run_test_s0_gs_kink_context_t), intent(in) :: ctx
        real(dp), intent(in) :: t
        real(dp), intent(out) :: point(2)
        real(dp) :: theta
        theta = 2*acos(-1.0_dp)*t
        point = [ctx%r0+radius*cos(theta),radius*sin(theta)]
    end subroutine run_test_s0_gs_kink_boundary

    subroutine run_test_s0_gs_kink_boundary_binding(self, t, point)
        class(run_test_s0_gs_kink_boundary_t), intent(in) :: self
        real(dp), intent(in) :: t
        real(dp), intent(out) :: point(2)
        call run_test_s0_gs_kink_boundary(self%data, t, point)
    end subroutine
end module test_s0_gs_kink_callbacks

program test_s0_gs_kink
    use test_s0_gs_kink_callbacks, only: run_test_s0_gs_kink
    implicit none
    call run_test_s0_gs_kink()
end program test_s0_gs_kink
