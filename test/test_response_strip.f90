program test_response_strip
    !! R-B1 small cases: resolved resistive strip in the periodic cylinder,
    !! sharp current step (closed-form Delta', closed-form ideal shielding).
    !!
    !! 1. Constant-psi tearing (m/n = 2/1, q0 = 1.2, r1 = 0.5, r_s Delta' = 5):
    !!    eigenvalue of the coupled strip within 2 % of the resistive-inertial
    !!    dispersion relation Delta' d = Dhat(Q) (FKR/Ara closed form) at
    !!    S = 1e6, 1e7.
    !! 2. Non-constant-psi (m/n = 1/1, q0 = 0.7, Delta' = infinity, resistive
    !!    internal kink): gamma within 2 % of (eta F'^2)^(1/3) at S = 1e6, 1e7.
    !! 3. Ideal limit (eta = 0, gamma = 1e-6): layer current within 1e-3 of the
    !!    ideal shielding current r_s [psi'] = 2m psi_ext r_s^m (S = 2m for the
    !!    step: the outer field of the shielded rational is the vacuum one).
    !! 4. Steady penetration (gamma = 1e-9, S = 1e6): psi(r_s)/psi_vac(r_s)
    !!    within 1 % of P = -S/(r_s Delta') = -2m/(r_s Delta').
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_response_strip, only: cyl_step_t, strip_options_t, &
        strip_solution_t, step_rational, step_fprime, step_delta_prime, &
        strip_growth_rate, strip_mesh, strip_solve, layer_width, &
        ri_dispersion_gamma
    implicit none

    logical :: all_passed
    all_passed = .true.
    call tearing(cyl_step_t(m=2, n=1, q0=1.2_dp, r1=0.5_dp), [1.0e6_dp, 1.0e7_dp])
    call tearing(cyl_step_t(m=1, n=1, q0=0.7_dp, r1=0.5_dp), [1.0e6_dp, 1.0e7_dp])
    call ideal_limit()
    call steady_penetration()
    if (.not. all_passed) error stop 1
    print '(a)', "test_response_strip: all checks passed"

contains

    subroutine record(c, d)
        logical, intent(in) :: c
        character(len=*), intent(in) :: d
        all_passed = all_passed .and. c
        if (c) then
            print '(a,a)', "  PASS ", d
        else
            print '(a,a)', "  FAIL ", d
        end if
    end subroutine record

    subroutine tearing(eq, slist)
        type(cyl_step_t), intent(in) :: eq
        real(dp), intent(in) :: slist(:)
        type(strip_options_t) :: opt
        type(strip_solution_t) :: sol
        complex(dp) :: g
        real(dp) :: eta, rs, dpr, gref
        integer :: k, it, st, st2
        character(len=160) :: msg
        rs = step_rational(eq)
        dpr = step_delta_prime(eq)
        if (dpr < 0.5_dp*huge(1.0_dp)) dpr = dpr/rs
        do k = 1, size(slist)
            eta = 1.0_dp/slist(k)
            call ri_dispersion_gamma(eta, step_fprime(eq), dpr, gref, st2)
            ! deliberately off guess: the root must not depend on it
            call strip_growth_rate(eq, eta, cmplx(1.3_dp*gref, 0.0_dp, dp), opt, &
                g, sol, it, st)
            write (msg, '(a,i0,a,i0,a,es8.1,a,es14.7,a,es14.7,a,i0)') "m/n = ", &
                eq%m, "/", eq%n, " S = ", slist(k), " gamma = ", real(g), &
                " RI = ", gref, " it ", it
            print '(a)', "  "//trim(msg)
            call record(st == 0 .and. st2 == 0 .and. sol%info == 0 .and. &
                abs(g/gref - 1.0_dp) <= 0.02_dp, "growth rate within 2 % of RI layer")
        end do
    end subroutine tearing

    subroutine ideal_limit()
        type(cyl_step_t) :: eq
        type(strip_options_t) :: opt
        type(strip_solution_t) :: sol
        real(dp), allocatable :: r(:)
        complex(dp) :: g
        real(dp) :: rs, ideal
        eq = cyl_step_t(m=2, n=1, q0=1.2_dp, r1=0.5_dp)
        rs = step_rational(eq)
        g = (1.0e-6_dp, 0.0_dp)
        call strip_mesh(rs, layer_width(0.0_dp, g, step_fprime(eq)), opt, r)
        call strip_solve(eq, 0.0_dp, g, r, sol)
        ideal = 2*eq%m*real(eq%psi_ext)*rs**eq%m
        print '(a,2es14.6,a,es10.3)', "  ideal limit: current ", sol%current, &
            " psi_s/psi_vac ", abs(sol%psi_s)/rs**eq%m
        call record(sol%info == 0 .and. abs(sol%current/ideal - 1.0_dp) <= 1.0e-3_dp, &
            "eta = 0 layer current within 1e-3 of ideal shielding 2m psi_vac")
    end subroutine ideal_limit

    subroutine steady_penetration()
        type(cyl_step_t) :: eq
        type(strip_options_t) :: opt
        type(strip_solution_t) :: sol
        real(dp), allocatable :: r(:)
        complex(dp) :: g
        real(dp) :: rs, eta, p_ref
        complex(dp) :: p
        eq = cyl_step_t(m=2, n=1, q0=1.2_dp, r1=0.5_dp)
        rs = step_rational(eq)
        eta = 1.0e-6_dp
        g = (1.0e-9_dp, 0.0_dp)
        call strip_mesh(rs, layer_width(eta, g, step_fprime(eq)), opt, r)
        call strip_solve(eq, eta, g, r, sol)
        p = sol%psi_s/(eq%psi_ext*rs**eq%m)
        p_ref = -2.0_dp*eq%m/step_delta_prime(eq)
        print '(a,2es14.6,a,es14.6)', "  steady: P ", p, " ref ", p_ref
        call record(sol%info == 0 .and. abs(p/p_ref - 1.0_dp) <= 0.01_dp, &
            "steady reconnected flux within 1 % of -S/(r_s Delta')")
    end subroutine steady_penetration

end program test_response_strip
