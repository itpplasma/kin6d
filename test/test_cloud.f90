program test_cloud
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_status, only: fortnum_status_t, FORTNUM_OK
    use kin6d_microscopic_cloud, only: cloud_model_t, cloud_force, cloud_step
    use kin6d_cloud_leaves, only: cloud_background_inside, &
                                cloud_background_outside, cloud_background_inside_jvp, &
                                  cloud_background_inside_vjp
    implicit none
    type(cloud_model_t) :: model
    type(fortnum_status_t) :: status
    real(dp) :: state(12), force(6), potential, f(3), u, expected(3), exact
    real(dp) :: tangent(4), adjoint(5), direction(5), weight(4), left, right
    integer :: i
    ! Gauss' law: a unit-radius sphere of total charge 8 gives E_r=8r inside.
    call cloud_background_inside(0.25_dp, -0.5_dp, 0.0_dp, 8.0_dp, 1.0_dp, &
                                 f(1), f(2), f(3), u)
    expected = [-2.0_dp, 4.0_dp, 0.0_dp]
    if (maxval(abs(f - expected)) > 1.0e-14_dp) error stop 'inside Gauss oracle'
    if (abs(u + 10.75_dp) > 1.0e-14_dp) error stop 'inside potential oracle'
    call cloud_background_outside(2.0_dp, 0.0_dp, 0.0_dp, 8.0_dp, 1.0_dp, &
                                  f(1), f(2), f(3), u)
    if (abs(f(1) + 2.0_dp) > 1.0e-14_dp) error stop 'outside Gauss oracle'
    if (abs(u + 4.0_dp) > 1.0e-14_dp) error stop 'outside potential oracle'
    model%sphere_charge = 0.0_dp
    state = 0.0_dp; state(1) = -0.5_dp; state(4) = 0.5_dp
    call cloud_force(state, model, force, potential, status)
    if (status%code /= FORTNUM_OK) error stop 'pair force failure'
    if (abs(force(1) + 1.0_dp) > 1.0e-14_dp) error stop 'Coulomb force oracle'
    if (abs(force(4) - 1.0_dp) > 1.0e-14_dp) error stop 'pair action reaction'
    if (abs(potential - 1.0_dp) > 1.0e-14_dp) error stop 'pair potential oracle'
    ! Single electron inside sphere: x''=-x. Integrator independent closed form.
    state = 0.0_dp; state(1) = 0.5_dp; state(4) = 0.25_dp
    model%sphere_charge = 27.0_dp
    do i = 1, 1000
        call cloud_step(state(:6), 0.001_dp, model, status)
        if (status%code /= FORTNUM_OK) error stop 'harmonic cloud step failure'
    end do
    exact = 0.5_dp*cos(1.0_dp) + 0.25_dp*sin(1.0_dp)
    if (abs(state(1) - exact) > 1.0e-6_dp) error stop 'harmonic orbit oracle'
    direction = [0.1_dp, 0.2_dp, -0.1_dp, 0.3_dp, 0.4_dp]
    weight = [0.5_dp, -0.3_dp, 0.2_dp, 0.1_dp]
    call cloud_background_inside_jvp(0.25_dp, -0.5_dp, 0.4_dp, 8.0_dp, &
                       1.0_dp, direction(1), direction(2), direction(3), direction(4), &
                           direction(5), tangent(1), tangent(2), tangent(3), tangent(4))
    ! Independent directional derivative of -n*x/R^3.
    left = -8.0_dp*direction(1) - 0.25_dp*direction(4) + &
           3*8.0_dp*0.25_dp*direction(5)
    if (abs(tangent(1) - left) > 1.0e-13_dp) error stop 'background JVP oracle'
    call cloud_background_inside_vjp(0.25_dp, -0.5_dp, 0.4_dp, 8.0_dp, &
                       1.0_dp, weight(1), weight(2), weight(3), weight(4), adjoint(1), &
                                     adjoint(2), adjoint(3), adjoint(4), adjoint(5))
    left = dot_product(weight, tangent); right = dot_product(direction, adjoint)
    if (abs(left - right) > 1.0e-13_dp) error stop 'cloud derivative pairing'
    state = 0.0_dp
    call cloud_force(state, model, force, potential, status)
    if (status%code == FORTNUM_OK) error stop 'coincident particles admitted'
    print *, 'cloud: Gauss, Coulomb action/reaction, harmonic orbit and AD pass'
end program
