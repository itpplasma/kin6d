program test_response_toroidal_strip
    !! Analytic forced annulus in SI, not a discrete self-manufactured RHS.
    !! A_Z=R^2, xi=(R^2,R Z,Z^2), Phi=R Z, all multiplied by exp(i n phi).
    !! Equilibrium p=pbase+alpha R^2, B_phi^2=C/R^2-mu0 alpha R^2,
    !! J_Z=-2 alpha R/B_phi satisfies J0 x B0=grad p exactly.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_mhd_ampere_n, only: mu0
    use kin6d_response_fluid, only: response_fluid_equilibrium_t
    use kin6d_response_toroidal_strip, only: &
        toroidal_strip_t, toroidal_strip_work_t, toroidal_strip_status_t, &
        toroidal_strip_init, &
        toroidal_strip_action, toroidal_strip_solve, toroidal_strip_vjp, &
        toroidal_strip_defects
    implicit none
    complex(dp), parameter :: iu = (0.0_dp, 1.0_dp), gamma = (0.7_dp, 0.2_dp)
    real(dp), parameter :: pbase = 0.4_dp/mu0, alpha = 0.03_dp/mu0
    type(toroidal_strip_t) :: s
    type(toroidal_strip_work_t) :: work
    type(toroidal_strip_status_t) :: st
    type(response_fluid_equilibrium_t), allocatable :: eq(:, :)
    real(dp), allocatable :: radial(:), radius(:, :), inverse(:, :, :, :)
    real(dp), allocatable :: equation_scale(:, :, :)
    complex(dp), allocatable :: x(:, :, :), y(:, :, :), exact(:, :, :), rhs(:, :, :)
    complex(dp), allocatable :: direction(:, :, :), trial(:, :, :), dy(:, :, :)
    complex(dp) :: div, pressure, gp(3), mom(3), ohm(3), mode, dmode
    real(dp) :: theta, r, z, c, sn, bp, j0, j1, error, scale, eps, div_b, div_j
    complex(dp) :: lhs, rhspair
    integer :: nr, nt, i, j, n, info, m, sign_n, a, boundary_mode
    character(16) :: command

    call get_command_argument(1, command)
    if (command == 'convergence') then
        call radial_convergence()
        stop
    end if
    nr = 7; nt = 9
    if (command == 'solve') then
        nr = 5; nt = 5
    end if
    allocate (radial(nr), radius(nr, nt), inverse(2, 2, nr, nt), eq(nr, nt), &
              x(6, nr, nt), y(6, nr, nt), exact(6, nr, nt), rhs(6, nr, nt), &
              direction(6, nr, nt), trial(6, nr, nt), dy(6, nr, nt), &
              equation_scale(6, nr, nt))
    equation_scale = 1.0_dp
    equation_scale(3:5, :, :) = 1.0_dp/mu0
    equation_scale(3, 1, :) = 1.0_dp
    equation_scale(3, nr, :) = 1.0_dp
    do i = 1, nr
        ! Deliberately nonuniform mesh.
        radial(i) = 0.3_dp + 0.2_dp*((i - 1.0_dp)/(nr - 1))**1.3_dp
    end do
    do sign_n = -1, 1, 2
        n = 3*sign_n
        do j = 1, nt
            theta = 2*acos(-1.0_dp)*(j - 1)/nt
            c = cos(theta); sn = sin(theta)
            do i = 1, nr
                r = 2.0_dp + radial(i)*c; z = radial(i)*sn
                radius(i, j) = r
                inverse(:, 1, i, j) = [c, -sn/radial(i)]
                inverse(:, 2, i, j) = [sn, c/radial(i)]
                bp = sqrt(9.0_dp/r**2 - mu0*alpha*r**2)
                j0 = -2*alpha*r/bp; j1 = (n*n - 4.0_dp)/mu0
                eq(i, j)%b0 = [0.0_dp, bp, 0.0_dp]
                eq(i, j)%j0 = [0.0_dp, 0.0_dp, j0]
                eq(i, j)%p = pbase + alpha*r**2
                eq(i, j)%gp = [2*alpha*r, 0.0_dp, 0.0_dp]
                eq(i, j)%hp = 0.0_dp
                eq(i, j)%hp(1, 1) = 2*alpha
                eq(i, j)%hp(2, 2) = 2*alpha
                eq(i, j)%rho = 1.0_dp/mu0
                eq(i, j)%eta = 0.2_dp*mu0
                x(:, i, j) = cmplx([0.0_dp, r*r, r*r, r*z, z*z, r*z], 0.0_dp, dp)
                div = 3*r + (2 + iu*n)*z
                pressure = -2*alpha*r**3 - eq(i, j)%adiabatic*eq(i, j)%p*div
                gp(1) = -6*alpha*r*r - eq(i, j)%adiabatic &
                        *(2*alpha*r*div + 3*eq(i, j)%p)
                gp(2) = iu*n*pressure/r
                gp(3) = -eq(i, j)%adiabatic*eq(i, j)%p*(2 + iu*n)
                mom = eq(i, j)%rho*gamma**2*x(3:5, i, j) + gp
                mom(1) = mom(1) + j1*bp - 2*r*j0
                mom(2) = mom(2) - iu*n*r*j0
                ohm(1) = -z - gamma*z*z*bp
                ohm(2) = -iu*n*z
                ohm(3) = -gamma*r*r - r + gamma*r*r*bp - eq(i, j)%eta*j1
                exact(:, i, j) = [ohm(1), ohm(3), mom, ohm(2)]
                rhs(:, i, j) = exact(:, i, j)
                if (i == 1 .or. i == nr) then
                    rhs(1, i, j) = c*ohm(1) + sn*ohm(3)
                    rhs(2, i, j) = -sn*x(1, i, j) + c*x(2, i, j)
                    rhs(3, i, j) = c*x(3, i, j) + sn*x(5, i, j)
                    rhs(5, i, j) = -sn*mom(1) + c*mom(3)
                end if
            end do
        end do
        call toroidal_strip_init(s, work, radial, radius, inverse, eq, n, info)
        call check(info == 0, 'valid mapped annulus')
        call toroidal_strip_action(s, work, gamma, x, y, .false.)
        error = maxval(abs(y - exact))/maxval(abs(exact))
        print '(a,i3,a,es12.4)', 'toroidal pressure n=', n, ' action error=', error
        call check(error < 2.0e-10_dp, 'analytic pressure, curvature, current and Phi')
        call toroidal_strip_defects(s, work, div_b, div_j)
        print '(a,2es12.4)', 'nodal div B [T/m], div J [A/m^3]=', div_b, div_j
        call check(div_b < 2.0e-9_dp, 'analytic div B = 0')
        call check(div_j < 2.0e-3_dp, 'analytic div J = 0, SI current units')
        call toroidal_strip_action(s, work, gamma, x, y, .true.)
        error = maxval(abs(y - rhs))/maxval(abs(rhs))
        call check(error < 2.0e-10_dp, 'two prescribed trace rings')

        if (command == 'solve') then
            trial = 0.0_dp
            call toroidal_strip_solve( &
                s, work, gamma, rhs, trial, equation_scale, &
                1.0e-11_dp, 2000, 6*nr*nt, st)
            error = maxval(abs(trial - x))/maxval(abs(x))
            print '(a,i0,a,es12.4,a,es12.4)', 'GMRES info=', st%info, &
                ' residual=', st%residual, ' state error=', error
            call check(st%info == 0, 'full-pressure GMRES convergence')
            call check(error < 2.0e-6_dp, 'forced solution against analytic state')
            call check(.not. st%certified, 'candidate status')
            direction = rhs*(0.13_dp, -0.08_dp)
            trial = 0.0_dp
            call toroidal_strip_solve( &
                s, work, gamma, direction, trial, equation_scale, &
                1.0e-11_dp, 2000, 6*nr*nt, st)
            call check(st%info == 0, 'equation-level tangent solve')
            call check( &
                maxval(abs(trial - x*(0.13_dp, -0.08_dp)))/maxval(abs(x)) &
                < 2.0e-6_dp, 'tangent solution against analytic drive derivative')
        else
            ! Highest retained positive/negative poloidal harmonics, with
            ! electrostatic forcing: a scalar derivative oracle in both charts.
            do m = -(nt - 1)/2, (nt - 1)/2
                trial = 0.0_dp; direction = 0.0_dp
                do j = 1, nt
                    theta = 2*acos(-1.0_dp)*(j - 1)/nt
                    mode = exp(iu*m*theta); dmode = iu*m*mode
                    do i = 1, nr
                        trial(6, i, j) = mode
                        direction(1, i, j) = sin(theta)*dmode/radial(i)
                        direction(2, i, j) = -cos(theta)*dmode/radial(i)
                        direction(6, i, j) = -iu*n*mode/radius(i, j)
                    end do
                end do
                call toroidal_strip_action(s, work, gamma, trial, dy, .false.)
                call check( &
                    maxval(abs(dy - direction)) < 2.0e-10_dp, &
                    'every poloidal harmonic, including m=0 and both signs')
            end do
            ! A fixed-equilibrium perturbation JVP is the linear action itself.
            direction = x*(0.13_dp, -0.08_dp) + trial
            eps = 1.0e-2_dp
            trial = x + eps*direction
            call toroidal_strip_action(s, work, gamma, trial, dy, .true.)
            call toroidal_strip_action(s, work, gamma, direction, exact, .true.)
            scale = maxval(abs(exact))
            call check( &
                maxval(abs((dy - y)/eps - exact))/scale < 2.0e-8_dp, &
                'fixed-equilibrium directional derivative')
            do j = 1, nt
                do i = 1, nr
                    do a = 1, 6
                        direction(a, i, j) = cmplx(sin(0.7_dp*a*i + j), &
                                                   cos(0.6_dp*a*j + i), dp)
                        trial(a, i, j) = cmplx(cos(0.3_dp*a*i - j), &
                                               sin(0.2_dp*a*j - i), dp) &
                                         /equation_scale(a, i, j)
                    end do
                end do
            end do
            do boundary_mode = 0, 1
                call toroidal_strip_action(s, work, gamma, direction, dy, &
                                           boundary_mode == 1)
                call toroidal_strip_vjp(s, work, gamma, trial, exact, &
                                        boundary_mode == 1)
                lhs = sum(conjg(dy)*trial)
                rhspair = sum(conjg(direction)*exact)
                error = abs(lhs - rhspair)/max(abs(lhs), abs(rhspair), 1.0_dp)
                print '(a,i0,a,es12.4)', &
                    'VJP boundary=', boundary_mode, ' error=', error
                call check(error < 2.0e-11_dp, 'complex storage VJP duality')
            end do
        end if
    end do

    ! Reject unsupported gauges/resolution and invalid charts before allocation.
    call toroidal_strip_init(s, work, radial, radius, inverse, eq, 0, info)
    call check(info /= 0, 'axisymmetric gauge is unsupported')
    inverse = 0.0_dp
    call toroidal_strip_init(s, work, radial, radius, inverse, eq, 3, info)
    call check(info /= 0, 'singular chart rejected')
    print '(a)', 'toroidal strip: all checks passed'
contains
    subroutine radial_convergence()
        type(toroidal_strip_t) :: strip
        type(toroidal_strip_work_t) :: scratch
        type(response_fluid_equilibrium_t), allocatable :: equilibrium(:, :)
        real(dp), allocatable :: u(:), major(:, :), inv(:, :, :, :)
        complex(dp), allocatable :: state(:, :, :), result(:, :, :), oracle(:, :, :)
        real(dp) :: errors(3), angle, cosine, sine, h
        complex(dp) :: harmonic, du, dt
        integer :: level, nn, it, ir, status
        do level = 1, 3
            nn = 2**(level + 2) + 1
            allocate (u(nn), major(nn, 17), inv(2, 2, nn, 17), equilibrium(nn, 17), &
                      state(6, nn, 17), result(6, nn, 17), oracle(6, nn, 17))
            h = 0.5_dp/(nn - 1)
            do ir = 1, nn
                u(ir) = 0.3_dp + h*(ir - 1)
            end do
            state = 0.0_dp; oracle = 0.0_dp
            do it = 1, 17
                angle = 2*acos(-1.0_dp)*(it - 1)/17
                cosine = cos(angle); sine = sin(angle)
                harmonic = exp(3*iu*angle)
                do ir = 1, nn
                    major(ir, it) = 2 + u(ir)*cosine
                    inv(:, 1, ir, it) = [cosine, -sine/u(ir)]
                    inv(:, 2, ir, it) = [sine, cosine/u(ir)]
                    equilibrium(ir, it)%p = 1.0_dp
                    state(6, ir, it) = u(ir)**6*harmonic
                    du = 6*u(ir)**5*harmonic
                    dt = 3*iu*u(ir)**6*harmonic
                    oracle(1, ir, it) = -cosine*du + sine*dt/u(ir)
                    oracle(2, ir, it) = -sine*du - cosine*dt/u(ir)
                    oracle(6, ir, it) = -2*iu*state(6, ir, it)/major(ir, it)
                end do
            end do
            call toroidal_strip_init( &
                strip, scratch, u, major, inv, equilibrium, 2, status)
            call check(status == 0, 'convergence chart')
            call toroidal_strip_action(strip, scratch, gamma, state, result, .false.)
            errors(level) = maxval(abs(result - oracle))
            print '(a,i0,a,es12.4)', &
                'radial nodes=', nn, ' analytic error=', errors(level)
            deallocate (u, major, inv, equilibrium, state, result, oracle)
        end do
        call check( &
            all(errors(:2)/errors(2:) > 8.0_dp), &
            'radial analytic error decreases at least cubically, including edges')
    end subroutine radial_convergence

    subroutine check(ok, message)
        logical, intent(in) :: ok
        character(*), intent(in) :: message
        if (.not. ok) then
            print '(a)', 'FAIL: '//message
            error stop 1
        end if
    end subroutine check
end program test_response_toroidal_strip
