program test_gs_profile_quadrature
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use kin6d_grad_shafranov, only: gs_state_t, gs_profile_load, &
                                    gs_profile_load_jvp, gs_profile_load_vjp
    implicit none
    type(gs_state_t) :: state, affine
    real(dp) :: flux(3), pp(3), ff(3), load(6), exact(6), step_exact(6)
    real(dp) :: tangent(6), dpsi(6), dflux(3), dpp(3), dff(3), dmu0, estimate
    real(dp) :: dual(6), gpsi(6), gflux(3), gpp(3), gff(3), gmu0, lhs, rhs
    real(dp) :: positive(6), negative(6), original_psi(6), h, weighted, invalid_mu
    integer :: status, passed
    passed = 0
    allocate (state%points(2, 6), state%triangles(6, 1), state%psi(6), &
    state%boundary(6), state%stiffness(6, 6, 1), state%gradient(2, 3, 1), state%area(1))
    ! Independent one-cell physical map R=2+xi, Z=eta, with native edge ordering.
    state%points = reshape([2.0_dp, 0.0_dp, 3.0_dp, 0.0_dp, 2.0_dp, 1.0_dp, &
                            2.5_dp, 0.0_dp, 2.5_dp, 0.5_dp, 2.0_dp, 0.5_dp], [2, 6])
    state%triangles(:, 1) = [1, 2, 3, 4, 5, 6]
    state%element_order = 2
    state%gradient(:, :, 1) = reshape([-1.0_dp, -1.0_dp, 1.0_dp, 0.0_dp, &
                                       0.0_dp, 1.0_dp], [2, 3])
    state%area = 0.5_dp; state%boundary = .false.; state%stiffness = 0
    original_psi = state%points(2, :)**2
    state%psi = original_psi
    flux = [0.0_dp, 0.25_dp, 1.0_dp]; pp = [0.0_dp, 0.0_dp, 0.75_dp]; ff = 0
    ! Closed rational integration of R*N_i*(eta^2-1/4) over eta>1/2.
    ! Integrate xi^i*eta^j by xi=0..1-eta, then eta=1/2..1.
    exact = [-649.0_dp/161280, -37.0_dp/8960, 1229.0_dp/53760, &
             53.0_dp/20160, 13.0_dp/672, 83.0_dp/4480]
    step_exact = [-43.0_dp/1920, -11.0_dp/480, 127.0_dp/1920, &
                  11.0_dp/480, 7.0_dp/60, 53.0_dp/480]
call gs_profile_load(state, 1.0_dp, flux, pp, ff, load, status, error_estimate=estimate)
    call require(status == 0, 'source-cut pressure integration status')
call close (maxval(abs(load - exact)), 0.0_dp, 1.0e-10_dp, 'six exact pressure moments')
call require(estimate >= 0 .and. estimate < 1.0e-10_dp, 'explicit quadrature estimator')
    dpsi = 1; dmu0 = 0; dflux = 0; dpp = 0; dff = 0
    call gs_profile_load_jvp(state, 1.0_dp, flux, pp, ff, dpsi, dmu0, dflux, dpp, dff, &
                             tangent, status)
    call require(status == 0, 'moving-cut continuous-source derivative')
    call close (maxval(abs(tangent - step_exact)), 0.0_dp, 1.0e-10_dp, &
                'six exact pressure psi derivative moments')
    dflux = 1
    call gs_profile_load_jvp(state, 1.0_dp, flux, pp, ff, dpsi, dmu0, dflux, dpp, dff, &
                             tangent, status)
    call require(status == 0, 'simultaneous flux and profile gauge tangent')
 call close (maxval(abs(tangent)), 0.0_dp, 1.0e-10_dp, 'constant-gauge derivative zero')
    state%psi = original_psi + 1048576.0_dp
    flux = flux + 1048576.0_dp
    call gs_profile_load(state, 1.0_dp, flux, pp, ff, load, status)
    call require(status == 0, 'large exactly representable source-profile gauge')
 call close (maxval(abs(load - exact)), 0.0_dp, 1.0e-9_dp, 'gauge-stable pressure load')
    state%psi = original_psi; flux = [0.0_dp, 0.25_dp, 1.0_dp]
    pp = 0; ff = [0.0_dp, 0.0_dp, 0.75_dp]
    call gs_profile_load(state, 1.0_dp, flux, pp, ff, load, status)
    call require(status == 0, 'source-cut FFprime integration status')
    weighted = dot_product(state%points(1, :), load)
    call close (weighted, 5.0_dp/192, 1.0e-10_dp, 'exact FFprime load tested against R')
    dpsi = 1; dflux = 0
    call gs_profile_load_jvp(state, 1.0_dp, flux, pp, ff, dpsi, dmu0, dflux, dpp, dff, &
                             tangent, status)
    call require(status == 0, 'FFprime derivative integration status')
    call close (dot_product(state%points(1, :), tangent), 1.0_dp/8, 1.0e-10_dp, &
                'exact FFprime psi derivative tested against R')
    pp = [0.0_dp, 0.0_dp, 0.75_dp]
    dpsi = [0.13_dp, -0.2_dp, 0.4_dp, 0.1_dp, -0.3_dp, 0.6_dp]
    dflux = [0.04_dp, -0.03_dp, 0.08_dp]
    dpp = [0.3_dp, -0.15_dp, 0.23_dp]; dff = [-0.2_dp, 0.12_dp, 0.31_dp]; dmu0 = 0.17_dp
    dual = [-0.4_dp, 0.2_dp, 0.3_dp, -0.1_dp, 0.5_dp, 0.7_dp]
    call gs_profile_load_jvp(state, 1.0_dp, flux, pp, ff, dpsi, dmu0, dflux, dpp, dff, &
                             tangent, status)
    call require(status == 0, 'combined source state/profile/parameter tangent')
    call gs_profile_load_vjp(state,1.0_dp,flux,pp,ff,dual,gpsi,gmu0,gflux,gpp,gff, &
                             status)
    call require(status == 0, 'source-load VJP status')
    lhs = dot_product(dual, tangent)
    rhs = dot_product(gpsi, dpsi) + gmu0*dmu0 + dot_product(gflux, dflux) + &
          dot_product(gpp, dpp) + dot_product(gff, dff)
    call close (lhs, rhs, 1.0e-10_dp, 'all-active integrated source JVP/VJP pairing')
    h = 1.0e-5_dp
    state%psi = original_psi + h*dpsi
    call gs_profile_load(state, 1 + h*dmu0, flux + h*dflux, pp + h*dpp, ff + h*dff, &
                         positive, status)
    call require(status == 0, 'positive independently rebuilt cut geometry')
    state%psi = original_psi - h*dpsi
    call gs_profile_load(state, 1 - h*dmu0, flux - h*dflux, pp - h*dpp, ff - h*dff, &
                         negative, status)
    call require(status == 0, 'negative independently rebuilt cut geometry')
    call close (maxval(abs((positive - negative)/(2*h) - tangent)), 0.0_dp, 2.0e-7_dp, &
                'independent moving-level source finite difference')
    state%psi = original_psi
    call gs_profile_load(state, 1.0_dp, flux, pp, ff, load, status, max_evaluations=1)
    call require(status /= 0, 'exhausted integration budget fails explicitly')
call close(maxval(abs(load)),0.0_dp,0.0_dp,'failed integration returns no partial load')
    invalid_mu = ieee_value(0.0_dp, ieee_quiet_nan)
    call gs_profile_load(state, invalid_mu, flux, pp, ff, load, status)
    call require(status /= 0, 'nonfinite mu0 domain rejected')
    ! Exact boundary elimination, independently specified stiffness column.
    state%boundary(2) = .true.
    state%stiffness(:, 2, 1) = [0.2_dp, 0.7_dp, -0.1_dp, 0.3_dp, 0.4_dp, -0.2_dp]
    state%stiffness(2, :, 1) = state%stiffness(:, 2, 1)
    state%psi = original_psi + 0.125_dp; flux = flux + 0.125_dp
    pp = [0.0_dp, 0.0_dp, 0.75_dp]; ff = 0
    call gs_profile_load(state, 1.0_dp, flux, pp, ff, load, status)
    call require(status == 0, 'nonzero boundary load status')
    positive = exact - 0.125_dp*state%stiffness(:, 2, 1); positive(2) = 0.125_dp
    call close (maxval(abs(load - positive)), 0.0_dp, 1.0e-10_dp, &
                'independent Dirichlet elimination')
    dpsi = 1; dflux = 0; dpp = 0; dff = 0; dmu0 = 0
    call gs_profile_load_jvp(state, 1.0_dp, flux, pp, ff, dpsi, dmu0, dflux, dpp, dff, &
                             tangent, status)
    call require(status == 0, 'boundary source tangent status')
    positive = step_exact - state%stiffness(:, 2, 1); positive(2) = 1
    call close (maxval(abs(tangent - positive)), 0.0_dp, 1.0e-10_dp, &
                'independent Dirichlet tangent elimination')
    call gs_profile_load_vjp(state,1.0_dp,flux,pp,ff,dual,gpsi,gmu0,gflux,gpp,gff, &
                             status)
    call require(status == 0, 'boundary source VJP status')
    call close (dot_product(dual, tangent), sum(gpsi), 1.0e-10_dp, &
                'Dirichlet source tangent/VJP pairing')
    state%boundary = .false.; state%stiffness = 0
    flux = [0.0_dp, 0.25_dp, 1.0_dp]
    state%psi = 0.25_dp
    call gs_profile_load(state, 1.0_dp, flux, pp, ff, load, status)
    call require(status == 0, 'flat cell primal on continuous source knot')
    call gs_profile_load_jvp(state, 1.0_dp, flux, pp, ff, dpsi, dmu0, dflux, dpp, dff, &
                             tangent, status)
    call require(status == 4, 'whole-cell kink derivative unavailable')
    ! The same physics integration must also repair affine P1 profile cuts.
    allocate (affine%points(2, 3), affine%triangles(3, 1), affine%psi(3), &
              affine%boundary(3), affine%stiffness(3, 3, 1), affine%area(1))
    affine%points = state%points(:, 1:3); affine%triangles(:, 1) = [1, 2, 3]
    affine%element_order = 1; affine%area = 0.5_dp; affine%boundary = .false.
    affine%stiffness = 0; affine%psi = affine%points(2, :)
    flux = [0.0_dp, 0.5_dp, 1.0_dp]; pp = [0.0_dp, 0.0_dp, 0.5_dp]; ff = 0
    call gs_profile_load(affine, 1.0_dp, flux, pp, ff, load(:3), status)
    call require(status == 0, 'affine P1 source-cut status')
    positive(:3) = [7.0_dp/1280, 11.0_dp/1920, 127.0_dp/3840]
    call close (maxval(abs(load(:3) - positive(:3))), 0.0_dp, 1.0e-10_dp, &
                'three exact affine pressure moments')
    pp = 0; ff = [0.0_dp, 0.0_dp, 0.5_dp]
    call gs_profile_load(affine, 1.0_dp, flux, pp, ff, load(:3), status)
    call require(status == 0, 'affine P1 FFprime source-cut status')
    call close (dot_product(affine%points(1, :), load(:3)), 1.0_dp/48, 1.0e-10_dp, &
                'exact affine FFprime load tested against R')
    print '(a,i0)', 'PASS independent source-cut/load/derivative checks: ', passed
contains
    subroutine require(condition, message)
        logical, intent(in)::condition
        character(*), intent(in)::message
        if (.not. condition) then
            print *, trim(message)
            error stop 'profile quadrature oracle'
        end if
        passed = passed + 1
    end subroutine require
    subroutine close (actual, wanted, tolerance, message)
        real(dp), intent(in)::actual, wanted, tolerance
        character(*), intent(in)::message
        if (abs(actual - wanted) > tolerance) then
            print *, trim(message), actual, wanted, tolerance
            error stop 'profile quadrature exact oracle'
        end if
        passed = passed + 1
    end subroutine close
end program test_gs_profile_quadrature
