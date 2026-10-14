program test_response_fluid
    !! Independent analytic waves and symmetry of a curved equilibrium.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_response_fluid
    use kin6d_response_fluid_leaves, only: fluid_pressure, &
                                  fluid_pressure_jvp, fluid_pressure_vjp, fluid_force, &
                                           fluid_force_jvp, fluid_force_vjp
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp), mu0 = 4.0e-7_dp*pi
    complex(dp), parameter :: iu = (0.0_dp, 1.0_dp)
    type(response_fluid_equilibrium_t) :: eq
    type(response_fluid_jet_t) :: jet
    type(response_fluid_result_t) :: res
    real(dp) :: wave(3), va2, cs2, k2, w2, total, disc, stiff(2, 2), c, radius
    complex(dp) :: gamma, amplitude
    integer :: mode, k
    character(16) :: arg

    eq%rho = 2.0_dp
    eq%b0(3) = sqrt(3.0_dp*mu0)
    eq%p = 0.7_dp
    wave = [0.8_dp, 0.0_dp, 1.3_dp]
    k2 = sum(wave**2)
    va2 = eq%b0(3)**2/(mu0*eq%rho)
    cs2 = eq%adiabatic*eq%p/eq%rho
    total = (va2 + cs2)*k2
    disc = sqrt(total**2 - 4*va2*cs2*k2*wave(3)**2)
    stiff(1, :) = [va2*k2 + cs2*wave(1)**2, cs2*wave(1)*wave(3)]
    stiff(2, :) = [stiff(1, 2), cs2*wave(3)**2]
    amplitude = (0.6_dp, -0.8_dp)
    do mode = 1, 2
        ! Uniform-plasma magnetosonic dispersion, including parallel compression.
        w2 = 0.5_dp*(total + (3 - 2*mode)*disc)
        jet = response_fluid_jet_t()
        jet%xi = amplitude*[stiff(1, 2), 0.0_dp, w2 - stiff(1, 1)]
        call plane_wave(wave)
        gamma = iu*sqrt(w2)
        call response_fluid_residual(eq, jet, gamma, res)
        call check(maxval(abs(res%momentum)) < 2.0e-14_dp, &
                   'fast/slow magnetosonic dispersion at finite pressure')
        call check(maxval(abs(res%ohm)) < 2.0e-14_dp, 'ideal frozen-in Ohm law')
        call check(abs(res%p1 + iu*eq%adiabatic*eq%p*sum(wave*jet%xi)) &
                   < 2.0e-14_dp, 'adiabatic compression pressure amplitude')
    end do

    ! Resistive shear Alfven wave: gamma^2 + eta k^2 gamma/mu0 + v_A^2 k^2 = 0.
    jet = response_fluid_jet_t()
    wave = [0.0_dp, 0.0_dp, 1.3_dp]
    eq%eta = 0.2_dp*mu0
    k2 = wave(3)**2
    gamma = -0.5_dp*eq%eta*k2/mu0 &
            + iu*sqrt(va2*k2 - (0.5_dp*eq%eta*k2/mu0)**2)
    jet%xi(1) = amplitude
    jet%a(2) = eq%rho*gamma**2*amplitude*mu0/(eq%b0(3)*k2)
    jet%b(1) = -iu*wave(3)*jet%a(2)
    jet%j(2) = k2*jet%a(2)/mu0
    jet%dxi(1, 3) = iu*wave(3)*amplitude
    call response_fluid_residual(eq, jet, gamma, res)
    call check(maxval(abs(res%momentum)) < 2.0e-14_dp, &
               'resistive Alfven dispersion: inertia and Lorentz force')
    call check(maxval(abs(res%ohm)) < 2.0e-14_dp, &
               'resistive Alfven dispersion: inductive electric field and eta J')
    call check(real(gamma, dp) < 0.0_dp, 'positive resistivity damps the wave')

    ! Curved equilibrium B = C(-y,x,0), J_z = 2C/mu0,
    ! p = p_axis - C^2(x^2+y^2)/mu0. A rigid translation is a zero mode.
    ! Its force balance needs BOTH J0 x b and the pressure Hessian.
    eq = response_fluid_equilibrium_t()
    jet = response_fluid_jet_t()
    c = sqrt(mu0)*0.4_dp
    radius = 1.7_dp
    eq%b0(2) = c*radius
    eq%j0(3) = 2*c/mu0
    eq%p = 2.0_dp - c**2*radius**2/mu0
    eq%gp(1) = -2*c**2*radius/mu0
    eq%hp(1, 1) = -2*c**2/mu0
    eq%hp(2, 2) = eq%hp(1, 1)
    jet%xi(1) = amplitude
    jet%a(3) = amplitude*c*radius
    jet%b(2) = -amplitude*c
    call response_fluid_residual(eq, jet, (0.0_dp, 0.0_dp), res)
    call check(abs(res%p1 - amplitude*2*c**2*radius/mu0) < 2.0e-14_dp, &
               'curved-equilibrium translation advects the full pressure gradient')
    call check(maxval(abs(res%momentum)) < 2.0e-14_dp, &
               'curved-equilibrium translation: pressure cancels J0 x b')
    call check(abs(res%gp1(1)) > 0.1_dp, 'curvature test has nonzero pressure force')

    ! Product rule checked by differencing a polynomial p1(x), not the kernel.
    call polynomial_pressure()
    call batch_oracle()
    call pressure_derivatives()
    call force_derivatives()
    print '(a)', 'test_response_fluid: all checks passed'
    if (command_argument_count() > 0) then
        call get_command_argument(1, arg)
        if (trim(arg) == 'benchmark') call benchmark()
    end if
contains
    subroutine check(ok, name)
        logical, intent(in) :: ok
        character(*), intent(in) :: name
        if (.not. ok) then
            print '(a)', 'FAIL: '//name
            error stop 1
        end if
        print '(a)', 'PASS: '//name
    end subroutine check

    subroutine plane_wave(kvec)
        real(dp), intent(in) :: kvec(3)
        integer :: i, j
        jet%a = cross(jet%xi, cmplx(eq%b0, 0.0_dp, dp))
        jet%b = iu*cross(cmplx(kvec, 0.0_dp, dp), jet%a)
        jet%j = iu*cross(cmplx(kvec, 0.0_dp, dp), jet%b)/mu0
        do j = 1, 3
            do i = 1, 3
                jet%dxi(i, j) = iu*kvec(j)*jet%xi(i)
            end do
        end do
        jet%gd = -kvec*sum(kvec*jet%xi)
    end subroutine plane_wave

    pure function cross(x, y) result(z)
        complex(dp), intent(in) :: x(3), y(3)
        complex(dp) :: z(3)
        z = [x(2)*y(3) - x(3)*y(2), x(3)*y(1) - x(1)*y(3), &
             x(1)*y(2) - x(2)*y(1)]
    end function cross

    subroutine polynomial_pressure()
        real(dp) :: x(3), xp(3), xm(3), h
        complex(dp) :: reference(3)
        x = [0.3_dp, -0.7_dp, 0.6_dp]
        eq = response_fluid_equilibrium_t()
        jet = response_fluid_jet_t()
        ! p = 3 + x*y + 0.2 z^2, xi = (x^2, yz, z^2).
        eq%p = 3 + x(1)*x(2) + 0.2_dp*x(3)**2
        eq%gp = [x(2), x(1), 0.4_dp*x(3)]
        eq%hp(1, 2) = 1; eq%hp(2, 1) = 1; eq%hp(3, 3) = 0.4_dp
        jet%xi = amplitude*[x(1)**2, x(2)*x(3), x(3)**2]
        jet%dxi(1, 1) = amplitude*2*x(1)
        jet%dxi(2, 2) = amplitude*x(3)
        jet%dxi(2, 3) = amplitude*x(2)
        jet%dxi(3, 3) = amplitude*2*x(3)
        jet%gd = amplitude*[2.0_dp, 0.0_dp, 3.0_dp]
        call response_fluid_residual(eq, jet, (0.0_dp, 0.0_dp), res)
        h = 1.0e-5_dp
        do k = 1, 3
            xp = x; xm = x
            xp(k) = xp(k) + h; xm(k) = xm(k) - h
            reference(k) = (p1(xp) - p1(xm))/(2*h)
        end do
        call check(abs(res%p1 - p1(x)) < 2.0e-14_dp, 'polynomial pressure advection')
        call check(maxval(abs(res%gp1 - reference)) < 2.0e-9_dp, &
                   'pressure gradient against independent polynomial differences')
    end subroutine polynomial_pressure

    function p1(x) result(p)
        real(dp), intent(in) :: x(3)
        complex(dp) :: p
        p = -amplitude*(x(1)**2*x(2) + x(1)*x(2)*x(3) &
                        + 0.4_dp*x(3)**3 + eq%adiabatic &
                        *(3 + x(1)*x(2) + 0.2_dp*x(3)**2)*(2*x(1) + 3*x(3)))
    end function p1

    subroutine batch_oracle()
        type(response_fluid_equilibrium_t) :: e(37)
        type(response_fluid_jet_t) :: v(37)
        type(response_fluid_result_t) :: r(37)
        complex(dp) :: gam(37), scale
        integer :: i
        e = eq
        do i = 1, size(e)
            scale = cmplx(0.3_dp*i, -0.2_dp*i, dp)
            v(i) = jet
            v(i)%xi = scale*jet%xi
            v(i)%dxi = scale*jet%dxi
            v(i)%gd = scale*jet%gd
        end do
        gam = (0.0_dp, 0.0_dp)
        call response_fluid_batch(e, v, gam, r)
        do i = 1, size(e)
            scale = cmplx(0.3_dp*i, -0.2_dp*i, dp)
            if (maxval(abs(r(i)%gp1 - scale*res%gp1)) > 1.0e-12_dp) &
                error stop 'parallel batch pressure mismatch'
            if (maxval(abs(r(i)%momentum - scale*res%momentum)) > 1.0e-12_dp) &
                error stop 'parallel batch force mismatch'
        end do
        call check(.true., 'parallel batch agrees with analytic linearity')
    end subroutine batch_oracle

    subroutine pressure_derivatives()
        real(dp) :: x(29), d(29), seed(4), tangent(4), adjoint(29)
        real(dp) :: plus(4), minus(4), h
        integer :: i
        do i = 1, 29
            x(i) = sin(0.3_dp*i)
            d(i) = cos(0.4_dp*i)
        end do
        seed = [0.2_dp, -0.4_dp, 0.7_dp, 0.3_dp]
        h = 1.0e-5_dp
        call fluid_pressure_jvp(x(1), x(2), x(3), x(4), x(5), x(6), x(7), &
                          x(8), x(9), x(10), x(11), x(12), x(13), x(14), x(15), x(16), &
                        x(17), x(18), x(19), x(20), x(21), x(22), x(23), x(24), x(25), &
                       x(26), x(27), x(28), x(29), d(1), d(2), d(3), d(4), d(5), d(6), &
                           d(7), d(8), d(9), d(10), d(11), d(12), d(13), d(14), d(15), &
                        d(16), d(17), d(18), d(19), d(20), d(21), d(22), d(23), d(24), &
                            d(25), d(26), d(27), d(28), d(29), tangent(1), tangent(2), &
                                tangent(3), tangent(4))
        call fluid_pressure_vjp(x(1), x(2), x(3), x(4), x(5), x(6), x(7), &
                          x(8), x(9), x(10), x(11), x(12), x(13), x(14), x(15), x(16), &
                        x(17), x(18), x(19), x(20), x(21), x(22), x(23), x(24), x(25), &
                       x(26), x(27), x(28), x(29), seed(1), seed(2), seed(3), seed(4), &
                           adjoint(1), adjoint(2), adjoint(3), adjoint(4), adjoint(5), &
                          adjoint(6), adjoint(7), adjoint(8), adjoint(9), adjoint(10), &
                      adjoint(11), adjoint(12), adjoint(13), adjoint(14), adjoint(15), &
                      adjoint(16), adjoint(17), adjoint(18), adjoint(19), adjoint(20), &
                      adjoint(21), adjoint(22), adjoint(23), adjoint(24), adjoint(25), &
                                adjoint(26), adjoint(27), adjoint(28), adjoint(29))
        x = x + h*d
        call fluid_pressure(x(1), x(2), x(3), x(4), x(5), x(6), x(7), x(8), &
                         x(9), x(10), x(11), x(12), x(13), x(14), x(15), x(16), x(17), &
                        x(18), x(19), x(20), x(21), x(22), x(23), x(24), x(25), x(26), &
                            x(27), x(28), x(29), plus(1), plus(2), plus(3), plus(4))
        x = x - 2*h*d
        call fluid_pressure(x(1), x(2), x(3), x(4), x(5), x(6), x(7), x(8), &
                         x(9), x(10), x(11), x(12), x(13), x(14), x(15), x(16), x(17), &
                        x(18), x(19), x(20), x(21), x(22), x(23), x(24), x(25), x(26), &
                            x(27), x(28), x(29), minus(1), minus(2), minus(3), minus(4))
        call check(maxval(abs(tangent - (plus - minus)/(2*h))) < 1.0e-8_dp, &
                   'pressure JVP against central differences of all inputs')
        call check(abs(dot_product(seed, tangent) - dot_product(d, adjoint)) &
                   < 1.0e-12_dp, 'pressure VJP duality with checked JVP')
    end subroutine pressure_derivatives

    subroutine force_derivatives()
        real(dp) :: x(26), d(26), seed(6), tangent(6), adjoint(26)
        real(dp) :: plus(6), minus(6), h
        integer :: i
        do i = 1, 26
            x(i) = sin(0.3_dp*i)
            d(i) = cos(0.4_dp*i)
        end do
        seed = [0.2_dp, -0.4_dp, 0.7_dp, 0.3_dp, -0.6_dp, 0.1_dp]
        h = 1.0e-5_dp
        call fluid_force_jvp(x(1), x(2), x(3), x(4), x(5), x(7), x(8), x(9), &
                        x(11), x(12), x(13), x(15), x(16), x(17), x(19), x(20), x(21), &
                        x(23), x(24), x(25), d(1), d(2), d(3), d(4), d(5), d(6), d(7), &
                          d(8), d(9), d(10), d(11), d(12), d(13), d(14), d(15), d(16), &
                        d(17), d(18), d(19), d(20), d(21), d(22), d(23), d(24), d(25), &
                             d(26), tangent(1), tangent(2), tangent(3), tangent(4), &
                             tangent(5), tangent(6))
        call fluid_force_vjp(x(1), x(2), x(3), x(4), x(5), x(7), x(8), x(9), &
                        x(11), x(12), x(13), x(15), x(16), x(17), x(19), x(20), x(21), &
                     x(23), x(24), x(25), seed(1), seed(2), seed(3), seed(4), seed(5), &
                             seed(6), adjoint(1), adjoint(2), adjoint(3), adjoint(4), &
                           adjoint(5), adjoint(6), adjoint(7), adjoint(8), adjoint(9), &
                      adjoint(10), adjoint(11), adjoint(12), adjoint(13), adjoint(14), &
                      adjoint(15), adjoint(16), adjoint(17), adjoint(18), adjoint(19), &
                      adjoint(20), adjoint(21), adjoint(22), adjoint(23), adjoint(24), &
                             adjoint(25), adjoint(26))
        x = x + h*d
        call fluid_force(x(1), x(2), x(3), x(4), x(5), x(6), x(7), x(8), &
                         x(9), x(10), x(11), x(12), x(13), x(14), x(15), x(16), x(17), &
                        x(18), x(19), x(20), x(21), x(22), x(23), x(24), x(25), x(26), &
                         plus(1), plus(2), plus(3), plus(4), plus(5), plus(6))
        x = x - 2*h*d
        call fluid_force(x(1), x(2), x(3), x(4), x(5), x(6), x(7), x(8), &
                         x(9), x(10), x(11), x(12), x(13), x(14), x(15), x(16), x(17), &
                        x(18), x(19), x(20), x(21), x(22), x(23), x(24), x(25), x(26), &
                         minus(1), minus(2), minus(3), minus(4), minus(5), minus(6))
        call check(maxval(abs(tangent - (plus - minus)/(2*h))) < 1.0e-8_dp, &
                   'momentum/Ohm JVP against central differences of all inputs')
        call check(abs(dot_product(seed, tangent) - dot_product(d, adjoint)) &
                   < 1.0e-12_dp, 'momentum/Ohm VJP duality with checked JVP')
    end subroutine force_derivatives

    subroutine benchmark()
        type(response_fluid_equilibrium_t), allocatable :: e(:)
        type(response_fluid_jet_t), allocatable :: v(:)
        type(response_fluid_result_t), allocatable :: r(:)
        complex(dp), allocatable :: gam(:)
        integer, parameter :: np = 100000, repeats = 100
        integer :: first, last, rate, i
        real(dp) :: norm
        allocate (e(np), v(np), r(np), gam(np))
        e = eq; v = jet; gam = (0.1_dp, 0.7_dp)
        call response_fluid_batch(e, v, gam, r)
        call system_clock(first, rate)
        do i = 1, repeats
            call response_fluid_batch(e, v, gam, r)
        end do
        call system_clock(last)
        norm = sum(abs(r%p1))
        print '(a,f12.6,a,es24.16)', 'batch seconds=', real(last - first, dp)/rate, &
            ' pressure L1=', norm
    end subroutine benchmark
end program test_response_fluid
