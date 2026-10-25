module kin6d_equilibrium_3d_driver
    ! D8 consumer of the boundary/profile initializer and native map Newton.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_equilibrium_harness, only: equilibrium_case_t, equilibrium_resolution_t, &
        equilibrium_start_t, equilibrium_run_report_t, D8_OK, D8_BAD_INPUT, D8_BAD_REPORT
    use kin6d_equilibrium_state, only: equilibrium_state_t
    use kin6d_ideal_map, only: ideal_map_t, ideal_map_solve, ideal_map_solve_report_t, &
        ideal_map_energy, ideal_map_reduce
    use kin6d_ideal_seed_3d, only: toroidal_profile_t, toroidal_seed_t, &
        ideal_seed_3d_init, ideal_seed_disk_mesh
    use kin6d_map_3d_equilibrium, only: map_3d_equilibrium, map_3d_sample
    use kin6d_map_3d_readback_leaves, only: map_3d_cartesian
    use omp_lib, only: omp_get_wtime
    use kin6d_grad_shafranov, only: gs_state_t
    use kin6d_curved_triangle, only: triangle_locate
    use kin6d_nested_oracles, only: nested_oracle_t, NESTED_INTEGER
    use kin6d_nested_oracle_leaves, only: nested_boundary_harmonic, nested_setup, &
        nested_integer_averages, nested_boundary_radius
    implicit none
    private
    public :: equilibrium_3d_driver_t, run_boundary_profile_equilibrium
    public :: nested_integer_driver
    ! Measurement only: seed coefficients captured before the Newton solve.
    real(dp), allocatable, public :: prof_x_seed(:, :, :)
    type :: equilibrium_3d_driver_t
        real(dp), allocatable :: boundary(:, :, :)
        type(toroidal_profile_t) :: profile
        real(dp) :: toroidal_flux = 0, tolerance = 1e-5_dp
        integer :: nfp = 2, max_iterations = 160
        type(nested_oracle_t) :: oracle
        type(ideal_map_t) :: map
        type(toroidal_seed_t) :: seed
        type(ideal_map_solve_report_t) :: solve
        real(dp) :: magnetic_error = huge(1.0_dp), pressure_error = huge(1.0_dp)
        real(dp) :: label_error = huge(1.0_dp)
        character(32) :: label_error_status = 'unavailable'
        character(32) :: force_error_status = 'unavailable'
        character(32) :: estimator_status = 'unavailable'
    end type
contains
    subroutine nested_integer_driver(oracle, mpol, ntor, driver, status)
        ! Only boundary and exact profiles/flux enter the initializer.
        ! Interior oracle fields are used exclusively in error diagnostics.
        type(nested_oracle_t), intent(in) :: oracle
        integer, intent(in) :: mpol, ntor
        type(equilibrium_3d_driver_t), intent(out) :: driver
        integer, intent(out) :: status
        real(dp), allocatable :: rz(:, :, :)
        real(dp) :: theta, phi, x(3), term(4), pi, limit, extent, pa1, pa2
        real(dp) :: volume, pmean, bmean, beta, iota
        integer :: angles, j, k, m, n, local, failed
        status = D8_BAD_INPUT
        if (.not. oracle%valid()) return
        if (oracle%family /= NESTED_INTEGER) return
        if (min(mpol, ntor) < 1) return
        driver%oracle = oracle
        pi = acos(-1.0_dp); angles = max(64, 4*max(mpol, ntor) + 1)
        allocate (rz(2, 0:angles - 1, 0:angles - 1))
        allocate (driver%boundary(4, -ntor:ntor, 0:mpol))
        failed = 0
        !$omp parallel do collapse(2) schedule(static) &
        !$omp private(j,k,theta,phi,x,local) reduction(+:failed)
        do j = 0, angles - 1
            do k = 0, angles - 1
                theta = 2*pi*k/angles; phi = pi*j/angles
                call oracle%boundary(theta, phi, x, local)
                failed = failed + local
                call nested_boundary_radius(x(1),x(2),rz(1,k,j))
                rz(2,k,j)=x(3)
            end do
        end do
        !$omp end parallel do
        if (failed /= 0) return
        driver%boundary = 0
        !$omp parallel do collapse(2) schedule(static) private(m,n,k,j,theta,phi,term)
        do m = 0, mpol
            do n = -ntor, ntor
                if (m == 0 .and. n < 0) cycle
                do j = 0, angles - 1
                    do k = 0, angles - 1
                        theta = 2*pi*k/angles; phi = pi*j/angles
                        call nested_boundary_harmonic(rz(1,k,j),rz(2,k,j),theta,phi, &
                            real(m,dp),real(n,dp),term(1),term(2),term(3),term(4))
                        driver%boundary(:,n,m) = driver%boundary(:,n,m) + term
                    end do
                end do
                driver%boundary(:,n,m) = driver%boundary(:,n,m)/angles**2
                if (m /= 0 .or. n /= 0) driver%boundary(:,n,m) = &
                    2*driver%boundary(:,n,m)
            end do
        end do
        !$omp end parallel do
        call nested_setup(oracle%eps,oracle%tau,oracle%delta,oracle%lam, &
            oracle%edge_pressure,limit,extent,pa1,pa2)
        call nested_integer_averages(oracle%eps,oracle%tau,pa1,oracle%delta, &
            volume,pmean,bmean,beta,iota,driver%toroidal_flux)
        driver%profile%iota = [iota]
        driver%profile%pressure = [pa1,oracle%edge_pressure-pa1]
        status = D8_OK
    end subroutine

    subroutine run_boundary_profile_equilibrium(spec,resolution,start,state,report,eq)
        type(equilibrium_case_t), intent(in) :: spec
        type(equilibrium_resolution_t), intent(in) :: resolution
        type(equilibrium_start_t), intent(in) :: start
        class(*), intent(inout) :: state
        type(equilibrium_run_report_t), intent(inout) :: report
        type(equilibrium_state_t), intent(out) :: eq
        type(gs_state_t) :: mesh
        integer :: status, solve_status
        real(dp), allocatable :: gradient(:), reduced(:)
        real(dp) :: energy, pf2_clock, pf2_solve
        report%status = D8_BAD_INPUT
        if (spec%lane /= 'map' .or. start%name /= 'legacy') return
        select type (state)
        type is (equilibrium_3d_driver_t)
            if (.not. allocated(state%boundary)) return
            if (.not. state%oracle%valid()) return
            if (state%tolerance <= 0 .or. state%tolerance >= 1) return
            select case(spec%name)
            case('landreman_3d_nested')
                if (abs(state%oracle%tau) > tiny(1.0_dp)) return
            case('issan_3d_nonsymmetric')
                if (abs(state%oracle%tau) <= tiny(1.0_dp)) return
            case default
                return
            end select
            pf2_clock = omp_get_wtime()
            call ideal_seed_disk_mesh(mesh,resolution%boundary_nodes, &
                resolution%degree,status)
            if (status /= 0) return
            call ideal_seed_3d_init(state%map,state%seed,mesh,state%boundary, &
                state%nfp,state%toroidal_flux,state%profile,status)
            if (status /= 0) return
            print '(a,i0)', 'initialized_free_dofs=',count(.not.state%map%fixed)
            print '(a,i0)', 'initialized_coefficient_dofs=',state%map%ndof
            print '(a,i0)', 'initialized_reduced_dofs=',state%map%nred
            flush(6)
            print *, "pf2_setup_s=", omp_get_wtime()-pf2_clock
            pf2_clock = omp_get_wtime()
            call metrics(state,status)
            if (status /= 0) return
            report%seed_error = state%magnetic_error
            report%status = D8_BAD_REPORT
            print *, "pf2_seed_metrics_s=", omp_get_wtime()-pf2_clock
            prof_x_seed = state%map%x
            pf2_clock = omp_get_wtime()
            call ideal_map_solve(state%map,state%tolerance,state%max_iterations, &
                state%solve,status)
            pf2_solve = omp_get_wtime()-pf2_clock
            print *, "pf2_solver_wall_s=", pf2_solve
            pf2_clock = omp_get_wtime()
            solve_status = status
            report%iterations = state%solve%iterations
            ! The native report changes gradient units after its first
            ! iteration. Check the returned state directly, with the same
            ! generated energy and the declared free increment space.
            allocate(gradient(state%map%ndof),reduced(state%map%nred))
            call ideal_map_energy(state%map,state%map%x,energy,status,gradient)
            if(status/=0) return
            call ideal_map_reduce(state%map,gradient,reduced)
            where(state%map%fixed) reduced=0
            report%residual = norm2(reduced) &
                /max(tiny(1.0_dp),state%solve%initial_gradient)
            call metrics(state,status)
            if (status /= 0) return
            report%oracle_error = state%magnetic_error
            report%residual_norm = 'relative free reduced-gradient l2'
            if (solve_status /= 0) return
            if (report%residual > state%tolerance) return
            call map_3d_equilibrium(state%map,state%seed,eq,status, &
                toroidal_flux=state%toroidal_flux)
            if (status /= 0) return
            print *, "pf2_postsolve_s=", omp_get_wtime()-pf2_clock
            report%status = D8_OK
        end select
    end subroutine

    subroutine metrics(state,status)
        type(equilibrium_3d_driver_t), intent(inout) :: state
        integer, intent(out) :: status
        real(dp) :: xi(2,3,24), grid(3), angles(24), y(2), q(3), values(7), x(3)
        real(dp) :: exact(20), b(3), sums(4), label, pi
        integer :: cells(3,24), ir,it,ip,e,local,failed
        pi=acos(-1.0_dp);grid=[0.25_dp,0.5_dp,0.75_dp]
        do it=1,24
            angles(it)=2*pi*(it-1+0.37_dp)/24
        end do
        status=1
        do it=1,24
            do ir=1,3
                y=[2+sqrt(grid(ir))*cos(angles(it)),sqrt(grid(ir))*sin(angles(it))]
                cells(ir,it)=0
                do e=1,state%map%cells
                    call triangle_locate(state%map%reference(:,state%map%conn(:,e)), &
                        y,xi(:,ir,it),local)
                    if(local/=0) cycle
                    if(minval([xi(:,ir,it),1-sum(xi(:,ir,it))]) < -1e-10_dp) cycle
                    cells(ir,it)=e;exit
                end do
                if(cells(ir,it)==0) return
            end do
        end do
        sums=0;label=0;failed=0
        !$omp parallel do collapse(3) schedule(static) &
        !$omp private(ip,it,ir,e,q,values,local,x,b,exact) &
        !$omp reduction(+:sums,failed) reduction(max:label)
        do ip=1,24
            do it=1,24
                do ir=1,3
                    e=cells(ir,it)
                    call map_3d_sample(state%map,state%seed,e,xi(:,ir,it), &
                        angles(ip)/state%nfp,q,values,local)
                    failed=failed+local
                    call map_3d_cartesian(q(1),q(2),q(3),values(1),values(2),values(3), &
                        x(1),x(2),x(3),b(1),b(2),b(3))
                    call state%oracle%point(x,exact,local)
                    failed=failed+local
                    sums=sums+[sum((b-exact(3:5))**2),sum(exact(3:5)**2), &
                        (values(4)-exact(2))**2,exact(2)**2]
                    if(state%oracle%family==NESTED_INTEGER) &
                        label=max(label,abs(exact(1)/state%oracle%delta-grid(ir)))
                end do
            end do
        end do
        !$omp end parallel do
        if(failed/=0 .or. min(sums(2),sums(4))<=0) return
        state%magnetic_error=sqrt(sums(1)/sums(2))
        state%pressure_error=sqrt(sums(3)/sums(4))
        if(state%oracle%family==NESTED_INTEGER) then
            state%label_error=label
            state%label_error_status='sampled diagnostic'
        end if
        status=0
    end subroutine
end module
