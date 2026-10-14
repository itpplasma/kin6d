program test_response_strip_readback
    !! Actual R-B1 resolved strip, without a singular fit. Independent ideal
    !! exterior solution for the step equilibrium: psi=0 inside r_s and
    !! psi_ext*(r^m - r_s^(2m) r^-m) outside. Phi=2 i m R0 psi, psi_N=r^2.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_response_strip, only: cyl_step_t, strip_options_t, strip_solution_t, &
        strip_mesh, strip_solve, layer_width, step_rational, step_fprime
    use kin6d_response_calibration, only: spot_weights, spot_jump
    implicit none
    type(cyl_step_t) :: eq
    type(strip_options_t) :: opt
    type(strip_solution_t) :: sol
    real(dp), allocatable :: r(:)
    real(dp), parameter :: r0 = 20.0_dp, spot = 5.0e-4_dp
    real(dp) :: rs, dq, xs, at, nodes(5, 2), w(5, 2), rad, refcurrent
    complex(dp) :: flux(5, 2), jump, ref
    integer :: side, lo, i, st, lev
    eq = cyl_step_t(m=2, n=1, q0=1.2_dp, r1=0.5_dp, psi_ext=(0.7_dp, -0.4_dp))
    rs = step_rational(eq); dq = eq%q0/eq%r1**2
    xs = spot/(eq%n*dq)
    opt = strip_options_t(hcore=0.025_dp, ratio=1.025_dp)
    do lev = 1, 2
        call strip_mesh(rs, layer_width(0.0_dp, cmplx(1.0e-6_dp, 0.0_dp, dp), &
            step_fprime(eq)), opt, r)
        call strip_solve(eq, 0.0_dp, (1.0e-6_dp, 0.0_dp), r, sol)
        if (sol%info /= 0) error stop 'strip solve'
        do side = 1, 2
            at = rs**2 + merge(-xs, xs, side == 1)
            lo = minloc(abs(r*r - at), 1) - 2
            if (lo < 1 .or. lo + 4 > size(r)) error stop 'stencil range'
            nodes(:, side) = r(lo:lo + 4)**2
            do i = 1, 5
                flux(i, side) = cmplx(0.0_dp, 2*eq%m*r0, dp)*sol%psi(lo + i - 1)
            end do
        end do
        call spot_weights(nodes, rs**2, eq%n, dq, spot, w, st)
        if (st /= 0) error stop 'stencil does not resolve spot'
        jump = spot_jump(w, flux)
        rad = sqrt(rs**2 + xs)
        ref = cmplx(0.0_dp, eq%m**2*r0, dp)*eq%psi_ext* &
            (rad**(eq%m - 2) + rs**(2*eq%m)*rad**(-eq%m - 2))
        refcurrent = 2*eq%m*rs**eq%m
        print '(a,i0,a,2es12.3)', 'level=', lev, ' spot/current errors=', &
            abs(jump/ref - 1), abs(sol%current/(eq%psi_ext*refcurrent) - 1)
        if (abs(jump/ref - 1) > 1.0e-3_dp) error stop 'analytic direct spot jump'
        if (abs(sol%current/(eq%psi_ext*refcurrent) - 1) > 1.0e-3_dp) &
            error stop 'analytic current integral'
        opt%ratio = 1.0125_dp; opt%hcore = 0.0125_dp
    end do
    print '(a)', 'resolved single-ring direct spot and current: PASS'
end program test_response_strip_readback
