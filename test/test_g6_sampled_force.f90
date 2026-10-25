program test_g6_sampled_force
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, &
        ieee_positive_inf, ieee_is_nan
    use g6_sampled_force, only: g6_force_accumulator_t, g6_force_vector, &
        g6_force_add, g6_force_result
    implicit none
    real(dp), parameter :: r0 = 6.2_dp, minor = 0.62_dp, ftor = 32.86_dp
    real(dp), parameter :: c = ftor/(12*r0**3), pp = -16*c
    real(dp), parameter :: zero(3) = 0, current(3) = [-2.0_dp, -3.0_dp, -1.0_dp]
    real(dp), parameter :: gradient(3) = [1.0_dp, 2.0_dp, -3.0_dp]
    real(dp), parameter :: expected(3, 3) = reshape( &
        [-10.0_dp, 4.0_dp, 3.0_dp, -1.0_dp, -3.0_dp, 6.0_dp, &
         1.0_dp, -2.0_dp, -1.0_dp], [3, 3])
    type(g6_force_accumulator_t) :: acc, rotated, scaled, empty
    real(dp) :: xyz(3), b(3), jhat(3), gp(3), force(3), qforce(3), rms
    real(dp) :: r, z, u, theta, phi, br, bz, bphi, gr, gz, tol, h
    integer :: i, status, checks
    checks = 0
    ! Independent Cartesian analytic field: B=(y,2z,3x), curl B=(-2,-3,-1).
    ! mu0*p=10+x+2y-3z gives grad(mu0*p)=(1,2,-3).
    do i = 1, 3
        xyz = 0
        xyz(i) = 1
        b = [xyz(2), 2*xyz(3), 3*xyz(1)]
        call g6_force_vector(b, current, gradient, force, status)
        call require(status == 0, '3D analytic force available')
        call require(maxval(abs(force - expected(:, i))) == 0, &
                     'independent signed Cartesian force vector')
        call g6_force_add(acc, b, current, gradient, status)
        call require(status == 0, 'analytic sample accumulation')
        call g6_force_vector(turn(b), turn(current), turn(gradient), &
                             qforce, status)
        call require(status == 0, 'proper rotation available')
        call require(maxval(abs(qforce - turn(expected(:, i)))) == 0, &
                     'right-handed rotation covariance')
        call g6_force_add(rotated, turn(b), turn(current), turn(gradient), status)
        call require(status == 0, 'rotated sample accumulation')
        ! B*2 and mu0*J*3 require mu0*grad(p)*6, hence force*6.
        call g6_force_vector(2*b, 3*current, 6*gradient, force, status)
        call require(status == 0, 'physical unit scaling available')
        call require(maxval(abs(force - 6*expected(:, i))) == 0, &
                     'no extra mu0 or automatic normalization')
        call g6_force_add(scaled, 2*b, 3*current, 6*gradient, status)
        call require(status == 0, 'scaled sample accumulation')
    end do
    call g6_force_result(acc, rms, status)
    call require(status == 0, 'analytic vector RMS available')
    call require(abs(rms - sqrt(59.0_dp)) < 1.0e-14_dp, &
                 'vector RMS divides by rows, not components')
    call g6_force_result(rotated, rms, status)
    call require(status == 0, 'rotated RMS available')
    call require(abs(rms - sqrt(59.0_dp)) < 1.0e-14_dp, 'RMS rotation invariant')
    call g6_force_result(scaled, rms, status)
    call require(status == 0, 'scaled RMS available')
    call require(abs(rms - 6*sqrt(59.0_dp)) < 1.0e-13_dp, &
                 'absolute T^2/m unit scaling')

    ! Analytic Solovev curl and the independently differentiated pressure.
    ! The test rounding allowance is not an equilibrium acceptance tolerance.
    do i = 1, 3
        u = 0.15_dp + 0.2_dp*i
        theta = 0.31_dp + 0.71_dp*i
        phi = 0.27_dp + 0.83_dp*i
        r = sqrt(r0*r0 + 2*r0*minor*u*cos(theta))
        z = r0*minor*u*sin(theta)/r
        br = -8*c*r*z
        bz = 4*c*(r*r - r0*r0 + 2*z*z)
        bphi = ftor/r
        b = [br*cos(phi) - bphi*sin(phi), &
             br*sin(phi) + bphi*cos(phi), bz]
        jhat = [-pp*r*sin(phi), pp*r*cos(phi), 0.0_dp]
        gr = 4*c*pp*r*(r*r - r0*r0 + 2*z*z)
        gz = 8*c*pp*r*r*z
        gp = [gr*cos(phi), gr*sin(phi), gz]
        call g6_force_vector(b, jhat, gp, force, status)
        call require(status == 0, 'analytic Solovev force available')
        tol = 64*epsilon(1.0_dp)*max(1.0_dp, norm2(gp))
        call require(norm2(force) <= tol, 'independent Solovev force balance')
    end do

    acc = empty
    call g6_force_result(acc, rms, status)
    call require(status /= 0 .and. ieee_is_nan(rms), 'empty metric unavailable')
    call g6_force_add(acc, zero, zero, zero, status)
    call require(status == 0, 'zero force accumulation')
    call g6_force_result(acc, rms, status)
    call require(status == 0 .and. rms == 0, 'zero force is available')
    acc = empty
    call g6_force_add(acc, zero, zero, [-3.0e200_dp, -4.0e200_dp, 0.0_dp], status)
    call require(status == 0, 'large finite force accumulation')
    call g6_force_result(acc, rms, status)
    call require(status == 0, 'scaled norm avoids squared overflow')
    call require(abs(rms/1.0e200_dp - 5) < 1.0e-13_dp, 'large force RMS oracle')
    call g6_force_vector([0.0_dp, 0.0_dp, 1.0e-200_dp], &
                         [0.0_dp, 1.0e200_dp, 0.0_dp], zero, force, status)
    call require(status == 0, 'finite product with disparate scales')
    call require(norm2(force - [1.0_dp, 0.0_dp, 0.0_dp]) < 1.0e-14_dp, &
                 'dimensioned finite product oracle')

    b = zero
    b(1) = ieee_value(0.0_dp, ieee_quiet_nan)
    call g6_force_vector(b, current, gradient, force, status)
    call unavailable('nonfinite B')
    jhat = zero
    jhat(2) = ieee_value(0.0_dp, ieee_positive_inf)
    call g6_force_vector(zero, jhat, zero, force, status)
    call unavailable('nonfinite mu0*J')
    gp = zero
    gp(3) = ieee_value(0.0_dp, ieee_quiet_nan)
    call g6_force_vector(zero, zero, gp, force, status)
    call unavailable('nonfinite mu0*grad(p)')
    h = huge(1.0_dp)
    call g6_force_vector([0.0_dp, 0.0_dp, 2.0_dp], &
                         [0.0_dp, h, 0.0_dp], zero, force, status)
    call unavailable('cross product multiplication overflow')
    call g6_force_vector([0.0_dp, 1.0_dp, 1.0_dp], &
                         [0.0_dp, h, -h], zero, force, status)
    call unavailable('cross product subtraction overflow')
    call g6_force_vector([0.0_dp, 0.0_dp, 1.0_dp], &
                         [0.0_dp, h, 0.0_dp], [-h, 0.0_dp, 0.0_dp], force, status)
    call unavailable('pressure-gradient subtraction overflow')
    acc = empty
    call g6_force_add(acc, [0.0_dp, 0.0_dp, 2.0_dp], &
                      [0.0_dp, h, 0.0_dp], zero, status)
    call require(status /= 0, 'bad force sample unavailable')
    call g6_force_add(acc, zero, zero, zero, status)
    call require(status /= 0, 'unavailable accumulation is sticky')
    call g6_force_result(acc, rms, status)
    call require(status /= 0 .and. ieee_is_nan(rms), 'no partial RMS after failure')
    acc = empty
    call g6_force_add(acc, zero, zero, [-h, -h, -h], status)
    call require(status == 0, 'finite components before overflowing RMS')
    call g6_force_result(acc, rms, status)
    call require(status /= 0 .and. ieee_is_nan(rms), 'overflowing RMS unavailable')
    print '(a,i0,a)', 'PASS: ', checks, ' external sampled force checks'
contains
    pure function turn(v) result(q)
        real(dp), intent(in) :: v(3)
        real(dp) :: q(3)
        q = [-v(2), v(1), v(3)]
    end function
    subroutine require(condition, label)
        logical, intent(in) :: condition
        character(*), intent(in) :: label
        checks = checks + 1
        if (.not. condition) then
            print '(2a)', 'FAIL: ', label
            error stop 1
        end if
    end subroutine
    subroutine unavailable(label)
        character(*), intent(in) :: label
        call require(status /= 0, label)
        call require(all(ieee_is_nan(force)), 'unavailable vector remains NaN')
    end subroutine
end program
