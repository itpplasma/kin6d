program test_gs_profile_newton
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_grad_shafranov, only: gs_state_t, initialize_gs, solve_gs, apply_gs
    use kin6d_gs_profile_law, only: gs_profile_law_t, gs_law_primitives
    use kin6d_gs_nonlinear, only: gs_profile_jacobian_t, gs_cubic_load, &
                 solve_gs_cubic_profiles, apply_gs_profile_jacobian, gs_profile_response, &
                 gs_normalize_current
    use kin6d_gs_readback, only: gs_readback_t, read_gs_equilibrium
    use kin6d_gs_majorant, only: gs_functional_majorant
    use kin6d_gs_certificate, only: gs_fixed_certifier_t
    use kin6d_equilibrium_certificate, only: equilibrium_certificate_t, certify, &
        EQ_CLAIM_RESIDUAL, EQ_CLAIM_PRACTICAL, EQ_CLAIM_CERTIFIED
    use fortnum_quadrature, only: gauss_legendre_ab
    implicit none
    type(gs_state_t), target :: state
    type(gs_profile_law_t) :: law
    type(gs_profile_jacobian_t) :: jac, scratch
    type(gs_readback_t) :: positive, negative
  real(dp), allocatable :: saved(:), rhs(:), plus(:), minus(:), v(:), w(:), jv(:), jw(:)
    real(dp) :: contour(2, 32), majorant, flux, balance, step, pressure, f, alpha
    real(dp) :: current, original_current, integral
    integer :: degree, status, i, n
    character(16) :: argument
    if (command_argument_count() > 0) then
        call get_command_argument(1, argument)
        if (trim(argument) == 'work') then
            call source_cut_oracle(3)
            call source_cut_oracle(2)
            stop
        end if
    end if
    call rank_one_response_oracle()
    do degree = 2, 3
        do i = 1, 32
            call boundary(real(i - 1, dp)/32, contour(:, i))
        end do
        call initialize_gs(state, contour, 0.045_dp, status, degree, curve=boundary)
        if (status /= 0) error stop 'curved mesh'
        call solve_gs(state, 1.0_dp, -0.3_dp, -0.2_dp, status)
        if (status /= 0) error stop 'Solovev reference'
        saved = state%psi; n = size(saved)
        law%knots = [-10.0_dp, 0.0_dp]
        allocate (law%pp(4, 1), law%ff(4, 1))
        law%pp = 0; law%ff = 0
        law%pp(1, :) = -0.3_dp; law%ff(1, :) = -0.2_dp
        law%f_edge = 5; law%p_edge = 0; law%normalized = .false.
        state%psi = 0
        call solve_gs_cubic_profiles(state, 1.0_dp, law, status)
        if (status /= 0) error stop 'constant cubic law'
        if (norm2(state%psi - saved)/norm2(saved) > 2.0e-8_dp) &
            error stop 'Solovev limit changed'
        call gs_law_primitives(law, -0.5_dp, pressure, f, status)
        if (status /= 0 .or. abs(pressure - 0.15_dp) > 1.0e-13_dp) &
            error stop 'pressure primitive'
        if (abs(f*f - 25.2_dp) > 1.0e-13_dp) error stop 'F primitive'
        call read_gs_equilibrium(state, 5.0_dp, -0.2_dp, 33, 80, &
                                 positive, status, law)
        if (status /= 0) error stop 'negative-axis readback'
        state%psi = -state%psi
        law%pp = -law%pp; law%ff = -law%ff
        call read_gs_equilibrium(state, 5.0_dp, 0.2_dp, 33, 80, &
                                 negative, status, law)
        if (status /= 0) error stop 'positive-axis readback'
        if (maxval(abs(positive%q + negative%q)) > 1.0e-10_dp) &
            error stop 'signed q reversal'
        if (abs(positive%phi_edge - negative%phi_edge) > 1.0e-11_dp) &
            error stop 'toroidal flux changed under poloidal reversal'
        law%f_edge = -5
        call read_gs_equilibrium(state, -5.0_dp, 0.2_dp, 33, 80, &
                                 negative, status, law)
        if (status /= 0) error stop 'negative F readback'
        if (maxval(abs(positive%q - negative%q)) > 1.0e-10_dp) &
            error stop 'signed F and q'
        if (abs(positive%phi_edge + negative%phi_edge) > 1.0e-11_dp) &
            error stop 'signed Phi'
        state%psi = saved
        law%normalized = .true.; law%knots = [0.0_dp, 1.0_dp]
        law%axis = minval(saved); law%f_edge = 5
        ! Independent constant-source integral: J_phi = pprime*R + FFprime/R
        ! is integrated by a boundary-area oracle below, not by the production load.
        law%pp = 0; law%ff = 0; law%pp(1, :) = -0.3_dp
        state%psi = 0
        call solve_gs_cubic_profiles(state, 1.0_dp, law, status)
        if (status /= 0) error stop 'normalized constant law'
        saved = state%psi
        original_current = -0.3_dp * 3.0_dp * acos(-1.0_dp) * 1.4_dp
        ! This fixture boundary is inspected below; the exact current uses its
        ! area and first R moment, rather than a frozen repository output.
        call gs_normalize_current(state, 1.0_dp, law, &
             2*original_current, current, status, 1.0e-10_dp)
        if (status /= 0) error stop 'current normalization'
        if (norm2(state%psi-2*saved)/norm2(saved) > 3.0e-5_dp) &
            error stop 'ellipse current scale disagrees with analytic area moment'
        if (abs(current/(2*original_current)-1) > 1.0e-12_dp) &
            error stop 'prescribed current'
        block
            real(dp), allocatable :: load(:)
            allocate(load(n))
            call gs_cubic_load(state, 1.0_dp, law, load, jac, status, &
                              1.0e-10_dp, source_integral=integral)
        end block
        if (status /= 0 .or. abs(integral/current-1) > 1.0e-9_dp) &
            error stop 'normalized current reintegration'
        ! A negative target must be rejected rather than changing polarity.
        call gs_normalize_current(state, 1.0_dp, law, -current, &
                                  integral, status)
        if (status == 0) error stop 'polarity-changing current accepted'
        state%psi = saved
        law%axis = minval(saved)
        law%pp(:, 1) = [-0.3_dp, -0.08_dp, -0.02_dp, -0.01_dp]
        law%ff = 0
        allocate (rhs(n), plus(n), minus(n), v(n), w(n), jv(n), jw(n))
        do i = 1, n
            v(i) = sin(real(i, dp)); w(i) = cos(real(2*i, dp))
        end do
        where (state%boundary)
            v = 0; w = 0
        end where
        call gs_cubic_load(state, 1.0_dp, law, rhs, jac, status, 1.0e-10_dp)
        if (status /= 0) error stop 'normalized cubic load'
        call apply_gs_profile_jacobian(state, jac, v, jv)
        call apply_gs_profile_jacobian(state, jac, w, jw, transpose=.true.)
        if (abs(dot_product(w, jv) - dot_product(jw, v)) > &
            1.0e-11_dp*norm2(w)*norm2(jv)) error stop 'Jacobian transpose'
        step = 1.0e-5_dp
        state%psi = saved + step*v
        call gs_cubic_load(state, 1.0_dp, law, plus, scratch, status, 1.0e-10_dp)
        if (status /= 0) error stop 'plus load'
        state%psi = saved - step*v
        call gs_cubic_load(state, 1.0_dp, law, minus, scratch, status, 1.0e-10_dp)
        if (status /= 0) error stop 'minus load'
        state%psi = saved
        call apply_gs(state, v, rhs)
        if (norm2((plus - minus)/(2*step) - (rhs - jv)) > &
            2.0e-6_dp*norm2(rhs - jv)) error stop 'normalized load finite difference'
        call gs_profile_response(state, jac, v, plus, status)
        if (status /= 0) error stop 'equation tangent'
        call apply_gs_profile_jacobian(state, jac, plus, rhs)
        if (norm2(rhs - v) > 1.0e-8_dp*norm2(v)) error stop 'tangent residual'
        call gs_profile_response(state, jac, w, minus, status, transpose=.true.)
        if (status /= 0) error stop 'equation adjoint'
        if (abs(dot_product(w, plus) - dot_product(minus, v)) > &
            1.0e-8_dp*norm2(w)*norm2(plus)) error stop 'implicit pairing'
        call gs_functional_majorant(state, 1.0_dp, 0.0_dp, 0.0_dp, majorant, &
                                    flux, balance, status, law=law, &
                                    estimate_stability=.true.)
        if (status /= 4) error stop 'missing nonlinear stability must be explicit'
        call require_claim(4, EQ_CLAIM_RESIDUAL, majorant)
        law%normalized = .false.; law%knots = [-10.0_dp, 0.0_dp]
        law%pp = 0; law%ff(:, 1) = [0.3_dp, 0.05_dp, 0.0_dp, 0.0_dp]
        call gs_functional_majorant(state, 1.0_dp, 0.0_dp, 0.0_dp, majorant, &
            flux, balance, status, law=law, estimate_stability=.true., &
            stability_estimate=alpha)
        if (status /= 5) error stop 'practical nonlinear estimate unavailable'
        if (alpha <= 0 .or. alpha >= 1) error stop 'positive source amplification'
        if (abs(majorant*alpha - flux - balance) > 1.0e-11_dp*majorant) &
            error stop 'linearized majorant scaling'
        call require_claim(5, EQ_CLAIM_PRACTICAL, majorant)
        call require_claim(0, EQ_CLAIM_CERTIFIED, (flux + balance)/0.5_dp, 0.5_dp)
        law%ff(:, 1) = [1.0e4_dp, 1.0e3_dp, 0.0_dp, 0.0_dp]
        call gs_functional_majorant(state, 1.0_dp, 0.0_dp, 0.0_dp, majorant, &
            flux, balance, status, law=law, estimate_stability=.true.)
        if (status /= 4) error stop 'unstable comparison must not report an estimate'
        call require_claim(4, EQ_CLAIM_RESIDUAL, flux + balance)
        deallocate (saved, rhs, plus, minus, v, w, jv, jw, law%pp, law%ff)
        call source_cut_oracle(degree)
    end do
   print *, 'PASS: curved cubic profiles, normalization, implicit derivatives and signs'
contains
    subroutine require_claim(lane_status, claim, bound, lower)
        ! The psi certifier reports the gs_functional_majorant status
        ! unchanged, the core claim for it, and the same bound.
        integer, intent(in) :: lane_status, claim
        real(dp), intent(in) :: bound
        real(dp), intent(in), optional :: lower
        type(gs_fixed_certifier_t) :: cert
        type(equilibrium_certificate_t) :: result
        integer :: stat
        cert%state => state
        cert%mu0 = 1.0_dp
        cert%law = law
        if (present(lower)) cert%stability_lower_bound = lower
        call certify(cert, result, stat)
        if (stat /= 0) error stop 'certificate failed'
        if (result%lane_status /= lane_status) error stop 'certificate lane status'
        if (result%claim /= claim) error stop 'certificate claim'
        if (abs(result%majorant - flux - balance) > 1.0e-12_dp*result%majorant) &
            error stop 'certificate majorant'
        if (abs(result%bound - bound) > 1.0e-12_dp*bound) error stop 'certificate bound'
    end subroutine
    subroutine rank_one_response_oracle()
        type(gs_state_t) :: cell
        type(gs_profile_jacobian_t) :: jacobian
        real(dp) :: answer(3), rhs(3), expected(3)
        integer :: stat
        ! Two active graph-Laplacian DOFs have B=[[2,-1],[-1,2]].
        ! The independent rational matrix J=B-[1/2,1/4]*[1/4,1/2]
        ! sends [2,-1] to [5,-4]; its transpose sends it to [77/16,-35/8].
        cell%ready = .true.
        cell%triangles = reshape([1, 2, 3], [3, 1])
        cell%boundary = [.true., .false., .false.]
        cell%psi = [0.0_dp, 0.0_dp, 0.0_dp]
        allocate (cell%stiffness(3, 3, 1))
        cell%stiffness = -1
        cell%stiffness(1, 1, 1) = 2
        cell%stiffness(2, 2, 1) = 2
        cell%stiffness(3, 3, 1) = 2
        cell%diagonal = [1.0_dp, 2.0_dp, 2.0_dp]
        allocate (jacobian%mass(3, 3, 1))
        jacobian%mass = 0
        jacobian%axis_load = [0.0_dp, 0.5_dp, 0.25_dp]
        jacobian%axis_gradient = [0.0_dp, 0.25_dp, 0.5_dp]
        expected = [0.0_dp, 2.0_dp, -1.0_dp]
        rhs = [0.0_dp, 5.0_dp, -4.0_dp]
        call gs_profile_response(cell, jacobian, rhs, answer, stat)
        if (stat /= 0 .or. norm2(answer - expected) > 1.0e-8_dp) &
            error stop 'rank-one rational tangent'
        rhs = [0.0_dp, 77.0_dp/16, -35.0_dp/8]
        call gs_profile_response(cell, jacobian, rhs, answer, stat, transpose=.true.)
        if (stat /= 0 .or. norm2(answer - expected) > 1.0e-8_dp) &
            error stop 'rank-one rational adjoint'
        ! B=[[1,2],[2,1]] is indefinite despite its positive diagonal.
        ! Negative CG curvature must retain the existing GMRES fallback.
        jacobian%mass(2, 2, 1) = 1; jacobian%mass(3, 3, 1) = 1
        jacobian%mass(2, 3, 1) = -3; jacobian%mass(3, 2, 1) = -3
        rhs = [0.0_dp, 0.0_dp, 3.0_dp]
        call gs_profile_response(cell, jacobian, rhs, answer, stat)
        if (stat /= 0 .or. norm2(answer - expected) > 1.0e-8_dp) &
            error stop 'indefinite rank-one fallback'
        ! J=[[1/2,-1],[-1,2]] is singular and this rhs is inconsistent.
        jacobian%mass = 0
        jacobian%axis_load = [0.0_dp, 1.0_dp, 0.0_dp]
        jacobian%axis_gradient = [0.0_dp, 1.5_dp, 0.0_dp]
        rhs = [0.0_dp, 1.0_dp, 0.0_dp]
        call gs_profile_response(cell, jacobian, rhs, answer, stat)
        if (stat == 0) error stop 'singular rank-one response accepted'
    end subroutine

    subroutine boundary(t, point)
        real(dp), intent(in) :: t
        real(dp), intent(out) :: point(2)
        point = [3 + cos(2*acos(-1.0_dp)*t), 1.4_dp*sin(2*acos(-1.0_dp)*t)]
    end subroutine

    subroutine source_cut_oracle(degree)
        integer, intent(in) :: degree
        type(gs_state_t) :: cell
        type(gs_profile_law_t) :: table
        type(gs_profile_jacobian_t) :: jacobian
        real(dp), parameter :: ref(2, 10) = reshape([0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, &
                                   0.0_dp, 1.0_dp, 1.0_dp/3, 0.0_dp, 2.0_dp/3, 0.0_dp, &
                             2.0_dp/3, 1.0_dp/3, 1.0_dp/3, 2.0_dp/3, 0.0_dp, 2.0_dp/3, &
                                         0.0_dp, 1.0_dp/3, 1.0_dp/3, 1.0_dp/3], [2, 10])
        real(dp) :: reference(2, 10), load(10), x(16), w(16), expected, a, cut
        integer :: i, ndof, status, evaluations
        ndof = (degree + 1)*(degree + 2)/2
        reference = ref
        if (degree == 2) reference(:, 4:6) = reshape([0.5_dp, 0.0_dp, &
                                                0.5_dp, 0.5_dp, 0.0_dp, 0.5_dp], [2, 3])
        cell%element_order = degree; cell%ready = .true.; cell%mapped = .true.
        cell%points = reference(:, :ndof)
        cell%points(1, :) = 1 + reference(1, :ndof)
        cell%points(2, :) = reference(2, :ndof)*(1 + 0.2_dp*reference(1, :ndof))
        cell%triangles = reshape([(i, i=1, ndof)], [ndof, 1])
        cell%psi = reference(1, :ndof)**degree
        cell%boundary = spread(.false., 1, ndof)
        cell%area = [0.5_dp]; cell%curved_cell = [.true.]
        allocate (cell%stiffness(ndof, ndof, 1)); cell%stiffness = 0
        a = 0.374_dp; cut = a**degree
        table%knots = [-1.0_dp, cut, 2.0_dp]
        allocate (table%pp(4, 2), table%ff(4, 2))
        table%pp = 0; table%ff = 0; table%pp(2, 2) = 1
        call gs_cubic_load(cell, 1.0_dp, table, load(:ndof), jacobian, status)
        if (status /= 0) error stop 'cubic source-cut load'
        call gauss_legendre_ab(16, a, 1.0_dp, x, w)
        expected = sum(w*(1 + x)*(x**degree - cut)*(1 + 0.2_dp*x)*(1 - x))
        if (abs(sum(load(:ndof)) - expected) > 2.0e-10_dp) &
            error stop 'true polynomial flux cuts on mapped geometry'
        if (degree == 3) then
            ! Three crossings on one inner slice; integrate the two positive
            ! intervals independently of the FE cut construction.
            cell%psi = (reference(2, :ndof) - 0.2_dp) &
                       *(reference(2, :ndof) - 0.5_dp)*(reference(2, :ndof) - 0.8_dp)
            table%knots = [-1.0_dp, 0.0_dp, 1.0_dp]
            table%pp = 0; table%pp(2, 2) = 1
            call gs_cubic_load(cell, 1.0_dp, table, load(:ndof), jacobian, status)
            if (status /= 0) error stop 'three-root source cut'
            call gauss_legendre_ab(16, 0.2_dp, 0.5_dp, x, w)
            expected = sum(w*(x - 0.2_dp)*(x - 0.5_dp)*(x - 0.8_dp) &
                           *((1 - x) + 0.6_dp*(1 - x)**2 + (0.2_dp/3)*(1 - x)**3))
            call gauss_legendre_ab(16, 0.8_dp, 1.0_dp, x, w)
            expected = expected + sum(w*(x - 0.2_dp)*(x - 0.5_dp)*(x - 0.8_dp) &
                                 *((1 - x) + 0.6_dp*(1 - x)**2 + (0.2_dp/3)*(1 - x)**3))
            if (abs(sum(load(:ndof)) - expected) > 1.0e-10_dp*expected) &
                error stop 'all three cubic crossings required'
        end if
        ! A small closed cut wholly inside the cell: edge crossings alone
        ! cannot locate its two vertical tangencies.
        cell%psi = (reference(1, :ndof) - 0.21_dp)**2 &
                   + (reference(2, :ndof) - 0.23_dp)**2
        cut = 0.011_dp**2
        table%knots = [-1.0_dp, cut, 2.0_dp]
        table%pp = 0
        table%pp(:, 1) = [1 + cut, -1.0_dp, 0.0_dp, 0.0_dp]
        call gs_cubic_load(cell, 1.0_dp, table, load(:ndof), jacobian, status)
        if (status /= 0) error stop 'closed interior source cut'
        ! Integral over the reference disk of (cut-rho^2)*(1+x)*(1+.2*x).
        expected = acos(-1.0_dp)*cut**2/2 &
                   *((1 + 0.21_dp)*(1 + 0.2_dp*0.21_dp) + 0.2_dp*cut/6)
        if (abs(sum(load(:ndof)) - expected) > 1.0e-11_dp*expected) &
            error stop 'missed interior source island'
        table%pp(:, 1) = [(1 + cut)**2, -2*(1 + cut), 1.0_dp, 0.0_dp]
        call gs_cubic_load(cell, 1.0_dp, table, load(:ndof), jacobian, status)
        if (status /= 0) error stop 'C1 interior source cut'
        expected = acos(-1.0_dp)*cut**3/3 &
            *((1 + 0.21_dp)*(1 + 0.2_dp*0.21_dp) + 0.2_dp*cut/8)
        if (abs(sum(load(:ndof)) - expected) > 5.0e-7_dp*expected) &
            error stop 'missed C1 source island'
        table%pp(:, 1) = [1 + cut, -1.0_dp, 0.0_dp, 0.0_dp]
        cell%psi = cut
        call gs_cubic_load(cell, 1.0_dp, table, load(:ndof), jacobian, status)
        if (status /= 4) error stop 'whole-cell kink derivative must be unavailable'
        ! A dense table of one smooth quadratic must not trigger work for
        ! hundreds of artificial cuts. Analytic integrals check both the load
        ! and the Jacobian; the evaluation budget guards the algorithmic cost.
        deallocate (table%pp, table%ff)
        table%knots = [(real(i, dp)/512, i=0, 512)]
        allocate (table%pp(4, 512), table%ff(4, 512))
        do i = 1, 512
            a = table%knots(i)
            table%pp(:, i) = [a*a, 2*a, 1.0_dp, 0.0_dp]
        end do
        table%ff = 0
        cell%psi = reference(1, :ndof)**degree
        call gs_cubic_load(cell, 1.0_dp, table, load(:ndof), jacobian, status, &
            tolerance=1.0e-10_dp, quadrature_evaluations=evaluations)
        if (status /= 0) error stop 'dense smooth source table'
        call gauss_legendre_ab(16, 0.0_dp, 1.0_dp, x, w)
        expected = sum(w*(1 + x)*x**(2*degree)*(1 + 0.2_dp*x)*(1 - x))
        if (abs(sum(load(:ndof)) - expected) > 1.0e-11_dp*expected) &
            error stop 'dense smooth load integral'
        expected = sum(w*(1 + x)*2*x**degree*(1 + 0.2_dp*x)*(1 - x))
        if (abs(sum(jacobian%mass) - expected) > 1.0e-11_dp*expected) &
            error stop 'dense smooth Jacobian integral'
        print *, 'dense smooth source evaluations:', degree, evaluations
        if (evaluations > 5000) error stop 'smooth source quadrature work explosion'
        ! A dense kinked table (kink at a knot, 513 knots crossing the cell) must
        ! integrate band by band: all-levels cut integration exhausted its budget.
        do i = 1, 512
            a = table%knots(i)
            table%pp(:, i) = 0
            if (a >= 0.25_dp) table%pp(:, i) = [a - 0.25_dp, 1.0_dp, 0.0_dp, 0.0_dp]
        end do
        call gs_cubic_load(cell, 1.0_dp, table, load(:ndof), jacobian, status, &
            tolerance=1.0e-10_dp, quadrature_evaluations=evaluations)
        if (status /= 0) error stop 'dense kinked source table'
        call gauss_legendre_ab(16, 0.25_dp**(1.0_dp/degree), 1.0_dp, x, w)
        expected = sum(w*(1 + x)*(x**degree - 0.25_dp)*(1 + 0.2_dp*x)*(1 - x))
        if (abs(sum(load(:ndof)) - expected) > 1.0e-9_dp*expected) &
            error stop 'dense kinked load integral'
        print *, 'dense kinked source evaluations:', degree, evaluations
        if (evaluations > 80000) error stop 'kinked source quadrature work explosion'
    end subroutine
end program
