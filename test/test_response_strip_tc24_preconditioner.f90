program test_response_strip_tc24_preconditioner
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_vacuum_surfaces, only: psi_grid_t, eqdsk_read
    use kin6d_response_strip_jets, only: strip_equilibrium_jets_t, &
        strip_equilibrium_grid
    use kin6d_response_fluid, only: response_fluid_equilibrium_t
    use kin6d_response_toroidal_strip, only: toroidal_strip_t, &
        toroidal_strip_work_t, toroidal_strip_status_t, toroidal_strip_init
    use kin6d_response_strip_preconditioner, only: strip_preconditioner_t, &
        strip_preconditioner_init, strip_preconditioned_solve
    implicit none
    integer, parameter :: nr = 7, nt = 25, n = 3
    type(psi_grid_t) :: g
    type(strip_equilibrium_jets_t) :: jets(nr, nt)
    type(response_fluid_equilibrium_t) :: eq(nr, nt)
    type(toroidal_strip_t) :: strip
    type(toroidal_strip_work_t) :: work
    type(toroidal_strip_status_t) :: st
    type(strip_preconditioner_t) :: pc
    real(dp) :: radial(nr), radius(nr, nt), inverse(2, 2, nr, nt)
    real(dp) :: scale(6, nr, nt), eta, momentum_scale, ohm_scale, err
    complex(dp) :: rhs(6, nr, nt), x(6, nr, nt), gamma
    character(512) :: file
    integer :: i, j, level, info
    call get_command_argument(1, file)
    call eqdsk_read(trim(file), g, info)
    call check(info == 0, 'read TC24')
    do i = 1, nr
        radial(i) = sqrt(0.94_dp + 0.02_dp*(i - 1)/(nr - 1))
    end do
    do level = 1, 5
        eta = 1.0e-7_dp*10.0_dp**(-2*(level - 1))
        if (level >= 3) eta = 0
        gamma = (2000.0_dp, 500.0_dp)
        if (level == 4) gamma = 0.3_dp*gamma
        if (level == 5) gamma = 0.1_dp*gamma
        call strip_equilibrium_grid(g, radial, nt, 1.0e-7_dp, eta, jets, info)
        call check(info == 0, 'TC24 map/B/pressure jets')
        do j = 1, nt
            do i = 1, nr
                radius(i, j) = jets(i, j)%x(1)
                inverse(:, :, i, j) = jets(i, j)%inverse
                eq(i, j) = jets(i, j)%fluid
                ohm_scale = max(1.0_dp, abs(gamma)*norm2(eq(i, j)%b0))
                momentum_scale = max(1.0_dp, eq(i, j)%rho*abs(gamma)**2, &
                    eq(i, j)%p, norm2(eq(i, j)%gp))
                scale(:, i, j) = [ohm_scale, ohm_scale, &
                    momentum_scale, momentum_scale, momentum_scale, ohm_scale]
                if (i == 1 .or. i == nr) scale(2:3, i, j) = 1
            end do
        end do
        call check(minval(eq%p) > 0, 'full TC24 pressure retained')
        call toroidal_strip_init(strip, work, radial, radius, inverse, eq, n, info)
        call check(info == 0, 'TC24 strip initialized')
        call strip_preconditioner_init(pc, strip, work, gamma, inverse, scale, info)
        call check(info == 0, 'TC24 radial band factors')
        ! Independent exact solution: A=xi=0, Phi=exp(i*n*phi).
        ! Only -grad Phi is nonzero, in the toroidal Ohm row. The map,
        ! pressure and magnetic field may be arbitrary represented backgrounds.
        rhs = 0
        rhs(6, :, :) = cmplx(0.0_dp, -real(n, dp), dp)/radius
        x = 0
        call strip_preconditioned_solve(pc, strip, work, rhs, x, &
            2.0e-10_dp, 1500, 600, st)
        err = max(maxval(abs(x(1:5, :, :))), maxval(abs(x(6, :, :) - 1)))
        print '(a,2es12.4,a,i0,a,2es12.4)', 'TC24 eta,|gamma|=', eta, abs(gamma), &
            ' iterations=', st%iterations, ' residual,error=', st%residual, err
        call check(st%info == 0, 'TC24 full-operator FGMRES convergence')
        call check(st%residual < 2.0e-10_dp, 'TC24 relative equation residual')
        call check(err < 2.0e-5_dp, 'TC24 analytic potential-drive solution')
        call check(.not. st%certified, 'TC24 inverse bound unavailable')
    end do
contains
    subroutine check(ok, label)
        logical, intent(in) :: ok
        character(*), intent(in) :: label
        if (.not. ok) then
            print *, 'FAIL ', label
            error stop 1
        end if
    end subroutine
end program
