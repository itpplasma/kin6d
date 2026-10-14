program test_response_strip_preconditioner
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64, real128
    use kin6d_mhd_ampere_n, only: mu0
    use kin6d_response_fluid, only: response_fluid_equilibrium_t
    use kin6d_response_toroidal_strip, only: toroidal_strip_t, &
        toroidal_strip_work_t, toroidal_strip_status_t, toroidal_strip_init, &
        toroidal_strip_action, toroidal_strip_solve, toroidal_strip_vjp
    use kin6d_response_strip_preconditioner, only: strip_preconditioner_t, &
        strip_preconditioner_init, strip_preconditioner_apply, &
        strip_preconditioned_solve, strip_preconditioner_bytes
    implicit none
    integer, parameter :: wp = merge(real128, dp, real128 > 0)
    type(toroidal_strip_t) :: strip
    type(toroidal_strip_work_t) :: work
    type(strip_preconditioner_t) :: pc
    type(toroidal_strip_status_t) :: st, baseline
    type(response_fluid_equilibrium_t), allocatable :: eq(:, :)
    real(dp), allocatable :: radial(:), radius(:, :), inverse(:, :, :, :)
    real(dp), allocatable :: scale(:, :, :)
    complex(dp), allocatable :: exact(:, :, :), rhs(:, :, :), x(:, :, :), ax(:, :, :)
    complex(dp), allocatable :: seed(:, :, :), adj(:, :, :), pre(:, :, :), padj(:, :, :)
    complex(dp) :: gamma
    real(dp) :: eta, err, pairing, setup_s, solve_s, apply_s, solve_tol, amplitude
    integer(int64) :: tick, tock, rate
    integer :: nr, nt, level, info, i, j, a, repeat, count, forward_iterations, n
    character(16) :: command
    call get_command_argument(1, command)
    count = 1
    if (command == 'limits') count = 7
    if (command == 'scaling') count = 3
    do level = 1, count
        nr = 7; nt = 9; n = 3
        if (command == 'negative') n = -3
        gamma = (0.7_dp, 0.2_dp); eta = 0.2_dp*mu0
        solve_tol = 2.0e-12_dp
        if (command == 'limits') then
            if (level <= 4) then
                eta = mu0*0.2_dp*10.0_dp**(-2*(level - 1))
            else
                eta = 0.0_dp
                gamma = gamma*10.0_dp**(-(level - 5))
            end if
        end if
        if (command == 'scaling') then
            gamma = (2.0_dp, 0.2_dp)
            nr = 2**(level + 2) + 1; nt = 2**(level + 2) + 1
            solve_tol = 1.2e-12_dp
        end if
        call setup(nr, nt, gamma, eta)
        call system_clock(tick, rate)
        call strip_preconditioner_init(pc, strip, work, gamma, inverse, scale, info)
        call system_clock(tock)
        setup_s = real(tock - tick, dp)/rate
        call check(info == 0, 'band factors')
        x = 0
        call system_clock(tick)
        call strip_preconditioned_solve(pc, strip, work, rhs, x, &
                                        solve_tol, 1000, min(240, 6*nr*nt), st)
        call system_clock(tock)
        solve_s = real(tock - tick, dp)/rate
        err = maxval(abs(x - exact))/maxval(abs(exact))
        print '(a,2i4,a,2es11.3,a,es11.3,a,i5,a,es11.3)', &
            'grid=', nr, nt, ' eta/mu0,gamma=', eta/mu0, abs(gamma), &
            ' residual=', st%residual, ' iterations=', st%iterations, ' error=', err
        print '(a,3es11.3,a,i0)', 'setup,solve,bytes=', setup_s, solve_s, &
            real(strip_preconditioner_bytes(pc), dp), ' info=', st%info
        call check(st%info == 0, 'FGMRES converges')
        call check(st%residual < 2.0e-10_dp, 'independent full equation residual')
        call check(err < 2.0e-5_dp, 'independent analytic forced state')
        forward_iterations = st%iterations
        call check(.not. st%certified, 'candidate, inverse bound unavailable')
        if (command == 'oracles' .or. command == 'negative') then
            do j = 1, nt
                do i = 1, nr
                    do a = 1, 6
                        seed(a, i, j) = cmplx(sin(real(7*a + 3*i + j, dp)), &
                            cos(real(2*a + i + 5*j, dp)), dp)
                    end do
                end do
            end do
            call strip_preconditioner_apply(pc, seed, pre, info)
            call check(info == 0, 'preconditioner action')
            call strip_preconditioner_apply(pc, exact, padj, info, .true.)
            call check(info == 0, 'preconditioner transpose')
            pairing = abs(sum(conjg(exact)*pre) - sum(conjg(padj)*seed)) &
                /max(1.0_dp, abs(sum(conjg(exact)*pre)))
            print '(a,es12.4)', 'preconditioner VJP pairing=', pairing
            call check(pairing < 2.0e-10_dp, 'preconditioner Hermitian transpose')
            ! The inverse JVP is the equation solve for a perturbed drive.
            x = 0
            call strip_preconditioned_solve(pc, strip, work, rhs*(0.13_dp, -0.08_dp), &
                x, 2.0e-10_dp, 1000, 180, st)
            call check(st%info == 0, 'equation-level drive tangent')
            call check(maxval(abs(x - exact*(0.13_dp, -0.08_dp))) &
                /maxval(abs(exact)) < 2.0e-6_dp, 'analytic drive tangent')
            do repeat = 1, 2
                amplitude = 10.0_dp**(24*repeat - 36)
                x = 0.05_dp*amplitude*exact
                call strip_preconditioned_solve(pc, strip, work, amplitude*rhs, &
                    x, 2.0e-10_dp, 1000, 180, st)
                call check(st%info == 0, 'amplitude-independent relative residual')
                call check(maxval(abs(x/amplitude - exact))/maxval(abs(exact)) &
                    < 2.0e-6_dp, 'small/large drive with nonzero initial guess')
            end do
            x = exact
            call strip_preconditioned_solve(pc, strip, work, 0.0_dp*rhs, &
                x, 2.0e-10_dp, 1000, 180, st)
            call check(st%info == 0, 'zero drive')
            call check(maxval(abs(x)) <= 0.0_dp, 'zero drive clears initial guess')
            ! Independently check the adjoint inverse by constructing A^*x.
            call toroidal_strip_vjp(strip, work, gamma, exact/scale, adj, .true.)
            x = 0
            call strip_preconditioned_solve(pc, strip, work, adj, &
                x, 1.0e-13_dp, 1000, 180, st, .true.)
            print '(a,i0,a,es11.3)', 'adjoint info=', st%info, ' residual=', st%residual
            print '(a,es11.3)', 'adjoint state error=', &
                maxval(abs(x - exact))/maxval(abs(exact))
            call check(st%info == 0, 'equation-level adjoint inverse')
            call check(maxval(abs(x - exact))/maxval(abs(exact)) < 2.0e-5_dp, &
                'adjoint inverse against supplied exact state')
            x = 0
            call strip_preconditioned_solve(pc, strip, work, 1.0e-12_dp*adj, &
                x, 1.0e-13_dp, 1000, 180, st, .true.)
            call check(st%info == 0, 'small adjoint drive relative residual')
            call check(maxval(abs(x/1.0e-12_dp - exact))/maxval(abs(exact)) &
                < 2.0e-5_dp, 'small adjoint drive analytic state')
            x = 0
            call toroidal_strip_solve(strip, work, gamma, rhs, x, scale, &
                2.0e-10_dp, 1000, 180, baseline)
            print '(a,i0,a,i0)', 'unpreconditioned iterations=', baseline%iterations, &
                ' info=', baseline%info
            call check(forward_iterations < baseline%iterations, 'iteration reduction')
        end if
        call system_clock(tick)
        do repeat = 1, 30
            call strip_preconditioner_apply(pc, rhs, pre, info)
        end do
        call system_clock(tock)
        apply_s = real(tock - tick, dp)/rate/30
        print '(a,es11.3,a,es16.8)', 'warm apply seconds=', apply_s, &
            ' L1=', sum(abs(pre))
        if (command == 'limits' .and. level == count) then
            call strip_preconditioner_init(pc, strip, work, (0.0_dp, 0.0_dp), &
                inverse, scale, info)
            call check(info /= 0, 'static singular rate explicitly rejected')
        end if
        deallocate (radial, radius, inverse, scale, eq, exact, rhs, x, ax, &
            seed, adj, pre, padj)
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
    subroutine setup(nr, nt, gamma, eta)
        integer, intent(in) :: nr, nt
        complex(dp), intent(in) :: gamma
        real(dp), intent(in) :: eta
        real(dp), parameter :: alpha = 0.03_dp/mu0, pbase = 0.4_dp/mu0
        complex(wp), parameter :: iq = (0.0_wp, 1.0_wp)
        complex(wp) :: div, pressure, gp(3), mom(3), ohm(3), state(6), rateq
        real(wp) :: rq, zq, cq, sq, p0q, g0q, hrq, htq, bq, jq, j1q
        real(dp) :: r, z, c, sn, th, bp, j0
        integer :: i, j
        allocate (radial(nr), radius(nr, nt), inverse(2, 2, nr, nt), &
            scale(6, nr, nt), eq(nr, nt), exact(6, nr, nt), rhs(6, nr, nt), &
            x(6, nr, nt), ax(6, nr, nt), seed(6, nr, nt), adj(6, nr, nt), &
            pre(6, nr, nt), padj(6, nr, nt))
        scale = 1; scale(3:5, :, :) = 1/mu0
        scale(3, 1, :) = 1; scale(3, nr, :) = 1
        do i = 1, nr
            if (command == 'scaling') then
                radial(i) = 0.3_dp + 0.2_dp*(i - 1.0_dp)/(nr - 1)
            else
                radial(i) = 0.3_dp + 0.2_dp*((i - 1.0_dp)/(nr - 1))**1.3_dp
            end if
        end do
        do j = 1, nt
            th = 2*acos(-1.0_dp)*(j - 1)/nt; c = cos(th); sn = sin(th)
            do i = 1, nr
                r = 2 + radial(i)*c; z = radial(i)*sn
                radius(i, j) = r
                inverse(:, 1, i, j) = [c, -sn/radial(i)]
                inverse(:, 2, i, j) = [sn, c/radial(i)]
                bp = sqrt(9/r**2 - mu0*alpha*r**2)
                j0 = -2*alpha*r/bp
                eq(i, j)%b0 = [0.0_dp, bp, 0.0_dp]
                eq(i, j)%j0 = [0.0_dp, 0.0_dp, j0]
                eq(i, j)%p = pbase + alpha*r*r
                eq(i, j)%gp = [2*alpha*r, 0.0_dp, 0.0_dp]
                eq(i, j)%hp = 0
                eq(i, j)%hp(1, 1) = 2*alpha; eq(i, j)%hp(2, 2) = 2*alpha
                eq(i, j)%rho = 1/mu0; eq(i, j)%eta = eta
                ! Independent polynomial oracle: compose the map and state
                ! before rounding, and differentiate using the represented
                ! background jets. Production leaves are not called here.
                cq = real(c, wp); sq = real(sn, wp)
                rq = 2.0_wp + real(radial(i), wp)*cq
                zq = real(radial(i), wp)*sq
                rateq = cmplx(gamma, kind=wp)
                state = cmplx([0.0_wp, rq*rq, rq*rq, rq*zq, zq*zq, rq*zq], &
                              kind=wp)
                exact(:, i, j) = cmplx(state, kind=dp)
                p0q = real(eq(i, j)%p, wp); g0q = real(eq(i, j)%gp(1), wp)
                hrq = real(eq(i, j)%hp(1, 1), wp)
                htq = real(eq(i, j)%hp(2, 2), wp)
                bq = real(bp, wp); jq = real(j0, wp)
                j1q = real(n*n - 4, wp)/real(mu0, wp)
                div = 3*rq + (2 + iq*n)*zq
                pressure = -state(3)*g0q - real(eq(i, j)%adiabatic, wp)*p0q*div
                gp(1) = -2*rq*g0q - state(3)*hrq - real(eq(i, j)%adiabatic, wp) &
                        *(g0q*div + 3*p0q)
                gp(2) = -(iq*n*state(3)/rq - state(4)/rq)*g0q - state(4)*htq &
                        - real(eq(i, j)%adiabatic, wp)*p0q*iq*n*div/rq
                gp(3) = -real(eq(i, j)%adiabatic, wp)*p0q*(2 + iq*n)
                mom = real(eq(i, j)%rho, wp)*rateq**2*state(3:5) + gp
                mom(1) = mom(1) + j1q*bq - 2*rq*jq
                mom(2) = mom(2) - iq*n*rq*jq
                ohm(1) = -zq - rateq*zq*zq*bq
                ohm(2) = -iq*n*zq
                ohm(3) = -rateq*rq*rq - rq + rateq*rq*rq*bq - real(eta, wp)*j1q
                rhs(:, i, j) = cmplx([ohm(1), ohm(3), mom, ohm(2)], kind=dp)
                if (i == 1 .or. i == nr) then
                    rhs(1, i, j) = cmplx(cq*ohm(1) + sq*ohm(3), kind=dp)
                    rhs(2, i, j) = cmplx(-sq*state(1) + cq*state(2), kind=dp)
                    rhs(3, i, j) = cmplx(cq*state(3) + sq*state(5), kind=dp)
                    rhs(5, i, j) = cmplx(-sq*mom(1) + cq*mom(3), kind=dp)
                end if
            end do
        end do
        call toroidal_strip_init(strip, work, radial, radius, inverse, eq, n, info)
        call check(info == 0, 'mapped annulus setup')
        call toroidal_strip_action(strip, work, gamma, exact, ax, .true.)
        call check(maxval(abs((ax - rhs)/scale)) < 2.0e-7_dp, &
            'analytic forcing independent of discrete action')
    end subroutine
end program
