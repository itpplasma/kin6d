program test_nonnested_error
    ! Efficiency of the NN1 N = 0 equilibrated bound against the exact
    ! Solov'ev flux of nn_solovev_fixture (p' only, F F' = 0): the bound must
    ! contain the true energy error |||psi - psi_h|||, stay within a modest
    ! efficiency index and converge with it.
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use nn_solovev_fixture, only: r0, minor, ftor, pp, mapping, coordinates, exact
    use kin6d_mhd_mesh, only: flux_mesh_t, flux_mesh_build
    use fortfem_meridian_feec_n, only: meridian_lagrange_nodes, meridian_qgrad
    use kin6d_equilibrium_nonnested, only: free_b_field_t, nn_report_t, nn_build, &
        nn_solve
    use kin6d_nonnested_error, only: nn_error_t, nn_axisymmetric_error, NN_ERROR_OK
    implicit none
    type(flux_mesh_t) :: fm
    type(free_b_field_t) :: field
    type(nn_report_t) :: report
    type(nn_error_t) :: estimate
    real(dp), allocatable :: rings(:), nodes(:, :)
    complex(dp), allocatable :: seed(:), bc(:), solution(:)
    real(dp) :: true_error(2), bound(2), u, theta, index_, rate
    integer :: k, lev, nr, i, st
    real(dp) :: solve_time, estimate_time
    integer(int64) :: clock0, clock1, clock_rate

    do k = 2, 3
    do lev = 1, 2
        nr = 8*2**(lev - 1)
        if (allocated(rings)) deallocate (rings, seed, bc)
        allocate (rings(nr))
        do i = 1, nr
            rings(i) = real(i, dp)/nr
        end do
        call flux_mesh_build(rings, 8, 1024, k, mapping, fm, st)
        if (st /= 0) error stop 'mesh'
        call nn_build(fm, k, 0, [r0, 0.0_dp], minor, ftor, 1.0_dp, field, st)
        if (st /= 0) error stop 'free-B build'
        field%coordinates => coordinates
        call meridian_lagrange_nodes(field%fs, nodes)
        allocate (seed(field%fs%n0), bc(field%fs%n0))
        bc = 0
        do i = 1, size(seed)
            call coordinates([r0, 0.0_dp], minor, nodes(:, i), u, theta)
            seed(i) = 0.1_dp*(1 - u*u)*cos(theta)*u
        end do
        call system_clock(clock0, clock_rate)
        call nn_solve(field, 0, bc, seed, report, st, pprime=pp, solution=solution)
        if (st /= 0) error stop 'free-B solve'
        call system_clock(clock1)
        solve_time = real(clock1 - clock0, dp)/clock_rate
        call system_clock(clock0)
        call nn_axisymmetric_error(field, solution, estimate, st)
        call system_clock(clock1)
        estimate_time = real(clock1 - clock0, dp)/clock_rate
        print '(a,3es12.4)', 'flux/oscillation/divergence: ', estimate%flux_defect, &
            estimate%oscillation, estimate%divergence_defect
        print '(a,3es12.4)', 'estimator/solve seconds, ratio: ', estimate_time, &
            solve_time, estimate_time/solve_time
        if (st /= NN_ERROR_OK) then
            print '(a,i0)', 'estimator status ', st
            error stop 'equilibrated estimate'
        end if
        if (estimate%claim /= 'practical-fixed-dirichlet-trace') error stop 'claim'
        bound(lev) = estimate%bound
        true_error(lev) = energy_error(real(solution, dp))
        index_ = bound(lev)/true_error(lev)
        print '(a,i0,a,i0,a,es11.4,a,es11.4,a,f7.3)', 'P', k, ' nr=', nr, &
            ' error ', true_error(lev), ' bound ', bound(lev), ' index ', index_
        if (index_ < 1 .or. index_ > 1.2_dp) error stop 'efficiency index'
    end do
    rate = log(bound(1)/bound(2))/log(2.0_dp)
    print '(a,f7.3)', 'bound convergence rate ', rate
    if (rate < k - 0.5_dp) error stop 'bound convergence rate'
    rate = log(true_error(1)/true_error(2))/log(2.0_dp)
    if (rate < k - 0.5_dp) error stop 'energy error convergence rate'
    end do
    call exact_recovery()
    print '(a)', 'PASS: NN1 N=0 equilibrated bound contains the Solovev error'

contains

    subroutine exact_recovery()
        ! psi = R^2 is harmonic for -div(grad psi/R), with nonzero constant
        ! flux (2,0). It is represented exactly by P3 on an affine mesh.
        ! The curved zero-field case also exercises the mapped source path.
        type(flux_mesh_t) :: mesh
        type(free_b_field_t) :: exact_field
        type(nn_error_t) :: res
        real(dp), allocatable :: xy(:, :)
        complex(dp), allocatable :: psi(:)
        integer :: geo, status

        do geo = 1, 3, 2
            call flux_mesh_build([0.5_dp, 1.0_dp], 8, 16, geo, mapping, mesh, status)
            if (status /= 0) error stop 'exact recovery mesh'
            call nn_build(mesh, 3, 0, [r0, 0.0_dp], minor, ftor, 1.0_dp, &
                exact_field, status)
            if (status /= 0) error stop 'exact recovery space'
            call meridian_lagrange_nodes(exact_field%fs, xy)
            psi = cmplx(0.0_dp*xy(1, :), 0.0_dp, dp)
            if (geo == 1) psi = cmplx(xy(1, :)**2, 0.0_dp, dp)
            call nn_axisymmetric_error(exact_field, psi, res, status, pprime=0.0_dp)
            if (status /= NN_ERROR_OK) error stop 'exact recovery estimator'
            print '(a,i0,a,es12.4)', 'exact recovery geometry P', geo, ' bound ', &
                res%bound
            if (res%bound > 1.0e-9_dp) error stop 'exact recovery above rounding'
        end do
    end subroutine exact_recovery

    real(dp) function energy_error(psi) result(error)
        ! int |grad(psi_exact - psi_h)|^2 / R dR dZ on the meridian rule with
        ! the analytic gradient (R B_Z, -R B_R) of the exact flux.
        real(dp), intent(in) :: psi(:)
        real(dp) :: x(2), wq, val(field%fs%nl0), grad(2, field%fs%nl0), g(2)
        real(dp) :: p, b(2)
        integer :: t, q, s
        error = 0
        do t = 1, field%fs%mm%topology%n_triangles
            do q = 1, field%fs%nq
                call meridian_qgrad(field%fs, t, q, x, wq, val, grad, s)
                if (s /= 0) error stop 'quadrature point'
                call exact(x, p, b)
                g = matmul(grad, psi(field%fs%sp(0)%dof(:, t))) &
                    - [x(1)*b(2), -x(1)*b(1)]
                ! wq carries R dR dZ; the energy density is |g|^2 / R.
                error = error + wq*sum(g*g)/x(1)**2
            end do
        end do
        error = sqrt(error)
    end function
end program test_nonnested_error
