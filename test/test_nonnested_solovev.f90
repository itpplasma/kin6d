module nn_solovev_fixture
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    real(dp), parameter :: r0 = 6.2_dp, minor = 0.62_dp, ftor = 32.86_dp
    real(dp), parameter :: c = ftor/(12*r0**3), pp = -16*c
contains
    subroutine mapping(u, theta, x)
        real(dp), intent(in) :: u, theta
        real(dp), intent(out) :: x(2)
        x(1) = sqrt(r0*r0 + 2*r0*minor*u*cos(theta))
        x(2) = r0*minor*u*sin(theta)/x(1)
    end subroutine
    subroutine coordinates(centre, scale, x, u, theta)
        real(dp), intent(in) :: centre(2), scale, x(2)
        real(dp), intent(out) :: u, theta
        real(dp) :: v(2)
        v = [(x(1)**2 - centre(1)**2)/2, x(1)*x(2)]/(centre(1)*scale)
        u = norm2(v); theta = atan2(v(2), v(1))
    end subroutine
    subroutine boundary(t, x)
        real(dp), intent(in) :: t
        real(dp), intent(out) :: x(2)
        call mapping(1.0_dp, 2*acos(-1.0_dp)*t, x)
    end subroutine
    subroutine exact(x, psi, b)
        real(dp), intent(in) :: x(2)
        real(dp), intent(out) :: psi, b(2)
        real(dp) :: r, z
        r = x(1); z = x(2)
        psi = c*((r*r - r0*r0)**2 + 4*r*r*z*z - 4*r0*r0*minor*minor)
        b = [-8*c*r*z, 4*c*(r*r - r0*r0 + 2*z*z)]
    end subroutine
end module

program test_nonnested_solovev
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use nn_solovev_fixture, only: r0, minor, ftor, pp, mapping, coordinates, &
        boundary, exact
    use kin6d_mhd_mesh, only: flux_mesh_t, flux_mesh_build
    use fortfem_meridian_feec_n, only: meridian_lagrange_nodes
    use kin6d_equilibrium_nonnested, only: free_b_field_t, nn_report_t, &
        nn_build, nn_solve, nn_export
    use kin6d_nonnested_residual, only: nn_residual
    use kin6d_equilibrium_residual, only: equilibrium_residual_t
    use kin6d_equilibrium_state, only: equilibrium_state_t, EQ_VALUES, &
        EQ_BR, EQ_BZ, EQ_PSI, EQ_P
    use kin6d_grad_shafranov, only: gs_state_t, initialize_gs, solve_gs
    use kin6d_gs_equilibrium, only: gs_equilibrium
    implicit none
    type(flux_mesh_t) :: fm
    type(free_b_field_t) :: field, varied
    type(nn_report_t) :: solve_report
    type(equilibrium_residual_t) :: stress
    type(equilibrium_state_t) :: free_eq, gs_eq
    type(gs_state_t) :: gs
    real(dp), allocatable :: rings(:), nodes(:, :), contour(:, :), points(:, :)
    complex(dp), allocatable :: seed(:), bc(:)
    real(dp) :: ef(3), ep(3), eg(3), er(3), values(EQ_VALUES, 128)
    real(dp) :: other(EQ_VALUES, 128), b(2), psi, theta, sum_b, sum_p, u, h
    real(dp) :: rate, dvalues(EQ_VALUES, 128), eps, tangent_error
    integer :: k, lev, nr, ntheta, i, j, st, statuses(128)
    character(16) :: arg
    call get_command_argument(1, arg)
    read (arg, *) k
    if (k < 2 .or. k > 3) error stop 'order'
    do lev = 1, 3
        nr = 8*2**(lev - 1); ntheta = 8*nr
        if (allocated(rings)) deallocate (rings, seed, bc, contour, points)
        allocate (rings(nr), contour(2, ntheta), points(3, 128))
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
        call nn_solve(field, 0, bc, seed, solve_report, st, pprime=pp)
        if (st /= 0) error stop 'free-B Newton'
        if (solve_report%energy_after >= solve_report%energy_before) &
            error stop 'offset did not decrease energy'
        if (solve_report%divergence > 1.0e-12_dp) error stop 'closed 2-form'
        do i = 1, ntheta
            call boundary(real(i - 1, dp)/ntheta, contour(:, i))
        end do
        h = 2*acos(-1.0_dp)*minor/ntheta
        call initialize_gs(gs, contour, h*h/2, st, k, curve=boundary)
        if (st /= 0) error stop 'GS mesh'
        call solve_gs(gs, 1.0_dp, pp, 0.0_dp, st)
        if (st /= 0) error stop 'GS solve'
        call gs_equilibrium(gs, gs_eq, st, ftor, 0.0_dp, pprime=pp, mu0=1.0_dp)
        if (st /= 0) error stop 'GS adapter'
        call nn_export(field, 'offset Solovev; fixed boundary', free_eq)
        do i = 1, 4
            do j = 1, 32
                theta = 2*acos(-1.0_dp)*(j - 0.37_dp)/32
                call mapping(0.15_dp + 0.18_dp*i, theta, points(:2, j + 32*(i - 1)))
                points(3, j + 32*(i - 1)) = 0.43_dp
            end do
        end do
        call free_eq%field%sample(points, values, statuses)
        if (any(statuses /= 0)) error stop 'free-B sampling'
        call gs_eq%field%sample(points, other, statuses)
        if (any(statuses /= 0)) error stop 'GS sampling'
        ef(lev) = 0; ep(lev) = 0; sum_b = 0; sum_p = 0
        eg(lev) = sum((values(EQ_BR, :) - other(EQ_BR, :))**2 &
            + (values(EQ_BZ, :) - other(EQ_BZ, :))**2)
        do i = 1, 128
            call exact(points(:2, i), psi, b)
            ef(lev) = ef(lev) + sum(([values(EQ_BR, i), values(EQ_BZ, i)] - b)**2)
            ep(lev) = ep(lev) + (values(EQ_PSI, i) - psi)**2
            sum_b = sum_b + sum(b*b); sum_p = sum_p + psi*psi
            if (abs(values(EQ_P, i) - pp*values(EQ_PSI, i)) > 1.0e-13_dp) &
                error stop 'pressure compatibility'
        end do
        ef(lev) = sqrt(ef(lev)/sum_b); eg(lev) = sqrt(eg(lev)/sum_b)
        ep(lev) = sqrt(ep(lev)/sum_p)
        call nn_residual(field, 0, stress, st)
        if (st /= 0) error stop 'D3'
        er(lev) = stress%normalized
        print '(a,i0,a,i0,a,4es12.4,a,es10.2)', 'P', k, ' nr=', nr, &
            ' field/flux/GS/D3 ', ef(lev), ep(lev), eg(lev), er(lev), &
            ' div ', solve_report%divergence
        if (lev == 1) then
            ! Equation-level sensitivity: solution scales linearly in p'.
            varied = field
            eps = 1.0e-4_dp
            call nn_solve(varied, 0, bc, seed, solve_report, st, pprime=pp*(1 + eps))
            if (st /= 0) error stop 'tangent solve'
            call varied%sample(points, dvalues, statuses)
            tangent_error = maxval(abs((dvalues(EQ_PSI, :) - values(EQ_PSI, :)) &
                /eps - values(EQ_PSI, :)))/maxval(abs(values(EQ_PSI, :)))
            print '(a,es12.4)', 'equation-level tangent error ', tangent_error
            if (tangent_error > 1.0e-9_dp) error stop 'implicit derivative'
        end if
    end do
    rate = log(ef(1)/ef(3))/log(4.0_dp)
    print '(a,4f9.3)', 'orders field/flux/GS/D3 ', rate, &
        log(ep(1)/ep(3))/log(4.0_dp), log(eg(1)/eg(3))/log(4.0_dp), &
        log(er(1)/er(3))/log(4.0_dp)
    if (rate < k - 0.3_dp) error stop 'field order'
    if (log(ep(1)/ep(3))/log(4.0_dp) < k + 0.6_dp) error stop 'flux order'
    if (log(eg(1)/eg(3))/log(4.0_dp) < k - 0.3_dp) error stop 'GS reduction order'
    if (er(3) >= er(1) .or. er(3) > 0.02_dp) error stop 'D3 reduction'
    print '(a)', 'PASS: NN1 independent Solovev reduction and offset return'
end program
