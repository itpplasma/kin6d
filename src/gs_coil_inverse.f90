module kin6d_gs_coil_inverse
    !! Inverse coil-current problem on the free-boundary solver (DESIGN 3.1 to
    !! 3.4, step P9).
    !!
    !!     min_I  Phi = 1/2 |r(psi(I))|^2 + 1/2 lambda^2 |L (I - I_ref)|^2
    !!     s.t.   F(psi; I) = 0,  lower <= I <= upper,  C I = d,
    !!
    !! with L = diag(1/S_k), S_k the current limit of circuit k. The reference
    !! I_ref is either a scenario (cfg%i_ref allocated) or zero.
    !!
    !! Constraint rows. Every row is a linear functional of the nodal flux at
    !! fixed points plus an optional multiple of the boundary flux psi_b,
    !!
    !!     r_i = ( sum_t c_it E_it(psi) + cb_i psi_b - target_i ) / scale_i,
    !!
    !! with E a point value or a point derivative of the finite element flux
    !! (the same jet as gs_element_jet) and psi_b the boundary flux of the
    !! free-boundary solver, a linear functional g_b . psi (envelope argument
    !! of P7). The row kinds are
    !!
    !!     GS_CI_ISOFLUX_PAIR  psi(p) - psi(q) (target 0 for isoflux)
    !!     GS_CI_XPOINT_BR     B_R(p) = -psi_Z / R   (X point: target 0)
    !!     GS_CI_XPOINT_BZ     B_Z(p) = +psi_R / R
    !!     GS_CI_PSI_POINT     psi(p)
    !!     GS_CI_ISOFLUX_LCFS  psi(p) - psi_b        (psi_b free: boundary flux)
    !!     GS_CI_PSI_BOUND     psi_b                 (psi_b fixed: flux state)
    !!
    !! Sensitivities. F = A psi - b(psi) - sum_k I_k load_k - extra gives
    !! d psi / d I_k = J^-1 load_k exactly, including the plasma response: one
    !! back substitution per circuit with the factorisation of P7 (reused, never
    !! refactored), threaded over the right-hand sides. D = d r / d I follows by
    !! applying the row functionals.
    !!
    !! Gauss-Newton. Each step solves the dense bounded problem
    !!
    !!     min |r + D dI|^2 + lambda^2 |L (I + dI - I_ref)|^2
    !!     s.t. bounds and ties,
    !!
    !! with fortnum BVLS in the variables dI_k / S_k, followed by an Armijo line
    !! search on Phi (each trial is a warm-started free-boundary Newton solve).
    !! The iterates stay feasible: they are convex combinations of feasible
    !! points and are clamped to the bounds.
    !!
    !! FreeGS-emulation mode (cfg%emulate_freegs). D is replaced by the vacuum
    !! response of the circuits (point values of the coil flux; for the psi_b
    !! columns the vacuum flux at the axis and the X point, the limiter point
    !! of highest flux or the mean over the target points, according to the
    !! solver mode), the residual still comes from the full solution, and the
    !! step is taken with a fixed relaxation. That is the FreeGS fixed point
    !! iteration: it stops where D_vac^T r + regulariser = 0 and not where the
    !! gradient of Phi vanishes, so Phi is higher than for Gauss-Newton on the
    !! same problem (DESIGN 3.4).
    !!
    !! GSVD report. The generalised singular values of (D, L) are the singular
    !! values of D diag(S). The coil combination of mode j is column j of the
    !! inverse transpose of the GSVD factor X; a near-uniform flux combination
    !! shows up as the smallest generalised singular value (the flux gauge,
    !! DESIGN 3.3).
    !!
    !! Fit variants (gs_ci_fit): unconstrained (no bounds, no ties), bounded
    !! (bounds and ties), fixed-at-limit (unconstrained refit, pin the circuits
    !! beyond their limit at the limit, refit the rest, repeat).
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_positive_inf, &
                                                                          ieee_is_finite
    use fortnum_bvls, only: bounded_least_squares, BVLS_OK, BVLS_MAX_ITERATIONS
    use fortnum_gsvd, only: generalized_svd
    use kin6d_grad_shafranov, only: gs_state_t, gs_cell_data
    use kin6d_gs_coils, only: gs_circuit_set_t, gs_circuit_flux_vectors
    use kin6d_gs_free_boundary, only: gs_fb_system_t, gs_fb_eval_t, gs_fb_factor_t, &
                                      gs_fb_evaluate, gs_fb_factor, gs_fb_solve_block, &
                                      gs_fb_release, gs_fb_newton, gs_fb_locate, &
                                      GS_FB_OK, GS_FB_XPOINT, GS_FB_LIMITER
    implicit none
    private

    integer, parameter, public :: GS_CI_ISOFLUX_PAIR = 1
    integer, parameter, public :: GS_CI_XPOINT_BR = 2
    integer, parameter, public :: GS_CI_XPOINT_BZ = 3
    integer, parameter, public :: GS_CI_PSI_POINT = 4
    integer, parameter, public :: GS_CI_ISOFLUX_LCFS = 5
    integer, parameter, public :: GS_CI_PSI_BOUND = 6

    integer, parameter, public :: GS_CI_UNCONSTRAINED = 1
    integer, parameter, public :: GS_CI_BOUNDED = 2
    integer, parameter, public :: GS_CI_FIXED_AT_LIMIT = 3

    integer, parameter, public :: GS_CI_OK = 0
    integer, parameter, public :: GS_CI_BAD_INPUT = 1
    integer, parameter, public :: GS_CI_FORWARD_FAILED = 2
    integer, parameter, public :: GS_CI_SENSITIVITY_FAILED = 3
    integer, parameter, public :: GS_CI_BVLS_FAILED = 4
    integer, parameter, public :: GS_CI_NOT_CONVERGED = 5
    integer, parameter, public :: GS_CI_LINE_SEARCH = 6
    integer, parameter, public :: GS_CI_GSVD_FAILED = 7

    type, public :: gs_ci_row_t
        integer :: kind = GS_CI_PSI_POINT
        real(dp) :: p(2) = 0.0_dp          !< first point (R, Z)
        real(dp) :: q(2) = 0.0_dp          !< second point of an isoflux pair
        real(dp) :: target = 0.0_dp        !< value of the functional to reach
        real(dp) :: scale = 1.0_dp         !< tolerance: the row is divided by it
    end type gs_ci_row_t

    type :: gs_ci_point_t
        !! One point evaluation as a nodal functional and its vacuum response.
        integer :: nn
        integer :: node(10)
        real(dp) :: w(10)
        real(dp), allocatable :: vac(:)    !< per unit circuit current
    end type gs_ci_point_t

    type, public :: gs_ci_model_t
        integer :: n_circuits = 0, n_rows = 0
        type(gs_ci_row_t), allocatable :: rows(:)
        type(gs_ci_point_t), allocatable :: pts(:)
        integer, allocatable :: pt_of(:, :)     !< (2, rows), 0 = unused
        real(dp), allocatable :: coef(:, :)     !< (2, rows)
        real(dp), allocatable :: cb(:)          !< psi_b coefficient per row
        type(gs_ci_point_t), allocatable :: lim(:)   !< solver limiter/target points
        type(gs_circuit_set_t) :: set
        real(dp), allocatable :: lower(:), upper(:)  !< circuit limits (A)
        real(dp), allocatable :: scale(:)            !< S_k (A), the limits
    end type gs_ci_model_t

    type, public :: gs_ci_config_t
        real(dp) :: lambda = 0.0_dp
        real(dp), allocatable :: i_ref(:)       !< scenario reference; absent = zero
        logical :: use_bounds = .true.
        real(dp), allocatable :: tie_c(:, :)    !< C (ties, circuits), C I = d
        real(dp), allocatable :: tie_d(:)
        logical :: emulate_freegs = .false.
        real(dp) :: relaxation = 1.0_dp         !< step factor of the emulation
        integer :: max_iterations = 30
        real(dp) :: step_tolerance = 1.0e-8_dp  !< |dI / S|_inf
        real(dp) :: armijo = 1.0e-4_dp
        real(dp) :: function_tolerance = 1.0e-12_dp  !< predicted decrease of Phi
    end type gs_ci_config_t

    type, public :: gs_ci_gsvd_t
        integer :: n = 0
        real(dp), allocatable :: gsv(:)      !< generalised singular values, ascending
        real(dp), allocatable :: modes(:, :) !< coil combination (A), |L v| = 1
        real(dp), allocatable :: xmat(:, :)  !< GSVD factor X, columns in mode order
        integer :: null_mode = 0             !< index of the smallest gsv
        real(dp) :: condition = 0.0_dp       !< gsv(n) / gsv(1)
    end type gs_ci_gsvd_t

    type, public :: gs_ci_result_t
        integer :: status = GS_CI_BAD_INPUT
        integer :: iterations = 0, passes = 0
        real(dp) :: phi = 0.0_dp, phi_data = 0.0_dp, phi_reg = 0.0_dp
        real(dp), allocatable :: currents(:), psi(:)
        real(dp), allocatable :: residual(:)  !< scaled rows at the solution
        real(dp), allocatable :: gradient(:)  !< dPhi/dI (A^-1 units of Phi)
        real(dp), allocatable :: history(:)   !< Phi per iteration
        integer, allocatable :: active(:)     !< last BVLS: -1 lower, 0 free, +1 upper
        logical, allocatable :: pinned(:)     !< fixed-at-limit circuits
        type(gs_ci_gsvd_t) :: gsvd
    end type gs_ci_result_t

    public :: gs_ci_build_model, gs_ci_row_values, gs_ci_set_targets
    public :: gs_ci_flux_sensitivities, gs_ci_row_jacobian, gs_ci_objective
    public :: gs_ci_solve, gs_ci_fit, gs_ci_gsvd_report, gs_ci_mode_coefficients

    interface
        subroutine dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
            import :: dp
            integer, intent(in) :: n, nrhs, lda, ldb
            real(dp), intent(inout) :: a(lda, *), b(ldb, *)
            integer, intent(out) :: ipiv(*), info
        end subroutine dgesv
    end interface

contains

    ! ------------------------------------------------------------------
    ! Model: row functionals and their vacuum response.
    ! ------------------------------------------------------------------

    subroutine gs_ci_build_model(state, sys, set, rows, model, status)
        !! Locate all points and build the nodal functionals. The circuit limits
        !! of `set` become the bounds and scales; edit model%lower, model%upper
        !! and model%scale afterwards to change them.
        type(gs_state_t), intent(in) :: state
        type(gs_fb_system_t), intent(in) :: sys
        type(gs_circuit_set_t), intent(in) :: set
        type(gs_ci_row_t), intent(in) :: rows(:)
        type(gs_ci_model_t), intent(out) :: model
        integer, intent(out) :: status
        integer :: m, i, npt, nlim, k, nc, stat

        status = GS_CI_BAD_INPUT
        m = size(rows)
        nc = size(set%incidence, 2)
        if (m .lt. 1 .or. nc .lt. 1) return
        if (size(sys%coil_load, 2) .ne. nc) return
        model%n_circuits = nc
        model%n_rows = m
        model%rows = rows
        model%set = set
        allocate (model%pt_of(2, m), model%coef(2, m), model%cb(m))
        npt = count(rows%kind .ne. GS_CI_PSI_BOUND &
                    ) + count(rows%kind .eq. GS_CI_ISOFLUX_PAIR)
        allocate (model%pts(npt))
        model%pt_of = 0
        model%coef = 0.0_dp
        model%cb = 0.0_dp
        npt = 0
        do i = 1, m
            select case (rows(i)%kind)
            case (GS_CI_ISOFLUX_PAIR)
                call add_point(rows(i)%p, 0, 1.0_dp, 1)
                call add_point(rows(i)%q, 0, -1.0_dp, 2)
            case (GS_CI_XPOINT_BR)
                call add_point(rows(i)%p, 2, -1.0_dp/rows(i)%p(1), 1)
            case (GS_CI_XPOINT_BZ)
                call add_point(rows(i)%p, 1, 1.0_dp/rows(i)%p(1), 1)
            case (GS_CI_PSI_POINT)
                call add_point(rows(i)%p, 0, 1.0_dp, 1)
            case (GS_CI_ISOFLUX_LCFS)
                call add_point(rows(i)%p, 0, 1.0_dp, 1)
                model%cb(i) = -1.0_dp
            case (GS_CI_PSI_BOUND)
                model%cb(i) = 1.0_dp
            case default
                return
            end select
            if (status .eq. GS_CI_FORWARD_FAILED) return
            if (.not. (rows(i)%scale .gt. 0.0_dp)) then
                status = GS_CI_BAD_INPUT
                return
            end if
        end do
        status = GS_CI_BAD_INPUT
        if (allocated(sys%cfg%points)) then
            nlim = size(sys%cfg%points, 2)
            allocate (model%lim(nlim))
            do k = 1, nlim
                call build_point(state, set, sys%cfg%points(1, k), &
                                 sys%cfg%points(2, k), 0, 1.0_dp, model%lim(k), stat)
                if (stat .ne. 0) return
            end do
        else
            allocate (model%lim(0))
        end if
        allocate (model%lower(nc), model%upper(nc), model%scale(nc))
        model%lower = -ieee_value(1.0_dp, ieee_positive_inf)
        model%upper = ieee_value(1.0_dp, ieee_positive_inf)
        model%scale = 1.0_dp
        if (allocated(set%i_min) .and. allocated(set%i_max)) then
            if (size(set%i_min) .eq. nc .and. size(set%i_max) .eq. nc) then
                model%lower = set%i_min
                model%upper = set%i_max
                do k = 1, nc
                    model%scale(k) = max(abs(set%i_min(k)), abs(set%i_max(k)))
                    if (.not. ieee_is_finite(model%scale(k)) .or. &
                        model%scale(k) .le. 0.0_dp) model%scale(k) = 1.0_dp
                end do
            end if
        end if
        status = GS_CI_OK
    contains
        subroutine add_point(p, kind, coef, slot)
            real(dp), intent(in) :: p(2), coef
            integer, intent(in) :: kind, slot
            integer :: s

            npt = npt + 1
            call build_point(state, set, p(1), p(2), kind, coef, model%pts(npt), s)
            if (s .ne. 0) then
                status = GS_CI_FORWARD_FAILED
                return
            end if
            model%pt_of(slot, i) = npt
            model%coef(slot, i) = 1.0_dp
        end subroutine add_point
    end subroutine gs_ci_build_model

    subroutine build_point(state, set, r, z, kind, factor, pt, status)
        !! Nodal weights of the value (kind 0), R derivative (1) or Z derivative
        !! (2) of the finite element flux at (r, z), times `factor` (the factor
        !! is applied to the weights); vacuum response per unit circuit current.
        type(gs_state_t), intent(in) :: state
        type(gs_circuit_set_t), intent(in) :: set
        real(dp), intent(in) :: r, z, factor
        integer, intent(in) :: kind
        type(gs_ci_point_t), intent(inout) :: pt
        integer, intent(out) :: status
        real(dp) :: xi(2), point(2), det, shape(10), grads(2, 10), hess(2, 2, 10)
        real(dp), allocatable :: g(:, :), gr(:, :), gz(:, :)
        integer :: e, ndof, stat, nc
        logical :: found

        status = 1
        call gs_fb_locate(state, r, z, e, xi, found)
        if (.not. found) return
        ndof = size(state%triangles, 1)
        call gs_cell_data(state, e, xi, point, det, shape(:ndof), grads(:, :ndof), &
                          hess(:, :, :ndof), stat)
        if (stat .ne. 0) return
        pt%nn = ndof
        pt%node = 0
        pt%node(:ndof) = state%triangles(:, e)
        pt%w = 0.0_dp
        select case (kind)
        case (0)
            pt%w(:ndof) = shape(:ndof)
            pt%w(1) = pt%w(1) + 1.0_dp - sum(shape(:ndof))
        case (1)
            pt%w(:ndof) = grads(1, :ndof)
            pt%w(1) = pt%w(1) - sum(grads(1, :ndof))
        case (2)
            pt%w(:ndof) = grads(2, :ndof)
            pt%w(1) = pt%w(1) - sum(grads(2, :ndof))
        case default
            return
        end select
        pt%w = factor*pt%w
        nc = size(set%incidence, 2)
        allocate (g(1, nc), gr(1, nc), gz(1, nc))
        call gs_circuit_flux_vectors(set, [r], [z], g, gr, gz)
        if (allocated(pt%vac)) deallocate (pt%vac)
        allocate (pt%vac(nc))
        select case (kind)
        case (0)
            pt%vac = factor*g(1, :)
        case (1)
            pt%vac = factor*gr(1, :)
        case (2)
            pt%vac = factor*gz(1, :)
        end select
        status = 0
    end subroutine build_point

    pure function point_apply(pt, psi) result(v)
        type(gs_ci_point_t), intent(in) :: pt
        real(dp), intent(in) :: psi(:)
        real(dp) :: v
        integer :: j

        v = 0.0_dp
        do j = 1, pt%nn
            v = v + pt%w(j)*psi(pt%node(j))
        end do
    end function point_apply

    subroutine gs_ci_row_values(model, psi, psi_bound, values)
        !! Unscaled row functionals f_i(psi) (residual is (f - target)/scale).
        !! The coefficient of each point is already folded into its weights.
        type(gs_ci_model_t), intent(in) :: model
        real(dp), intent(in) :: psi(:), psi_bound
        real(dp), intent(out) :: values(:)
        integer :: i, t

        do i = 1, model%n_rows
            values(i) = model%cb(i)*psi_bound
            do t = 1, 2
                if (model%pt_of(t, i) .gt. 0) values(i) = values(i) &
                                        + point_apply(model%pts(model%pt_of(t, i)), psi)
            end do
        end do
    end subroutine gs_ci_row_values

    subroutine gs_ci_set_targets(model, psi, psi_bound)
        !! Targets := the row functionals of the given solution (twin data).
        type(gs_ci_model_t), intent(inout) :: model
        real(dp), intent(in) :: psi(:), psi_bound
        real(dp), allocatable :: f(:)

        allocate (f(model%n_rows))
        call gs_ci_row_values(model, psi, psi_bound, f)
        model%rows%target = f
    end subroutine gs_ci_set_targets

    ! ------------------------------------------------------------------
    ! Sensitivities.
    ! ------------------------------------------------------------------

    subroutine gs_ci_flux_sensitivities(sys, ev, fac, ds, status)
        !! ds(:, k) = d psi / d I_k = J^-1 load_k with the factorisation `fac`
        !! of the converged state; threaded over the circuits.
        type(gs_fb_system_t), intent(in) :: sys
        type(gs_fb_eval_t), intent(in) :: ev
        type(gs_fb_factor_t), intent(inout) :: fac
        real(dp), allocatable, intent(out) :: ds(:, :)
        integer, intent(out) :: status
        integer :: nc, stat

        nc = size(sys%coil_load, 2)
        allocate (ds(sys%n, nc))
        call gs_fb_solve_block(sys, ev, fac, sys%coil_load, ds, stat)
        status = GS_CI_OK
        if (stat /= GS_FB_OK) status = GS_CI_SENSITIVITY_FAILED
    end subroutine gs_ci_flux_sensitivities

    subroutine gs_ci_row_jacobian(model, sys, ev, psi, ds, dmat, emulate)
        !! D(i, k) = d r_i / d I_k from ds (full response) or, with `emulate`,
        !! from the vacuum response of the circuits (ds is then unused).
        type(gs_ci_model_t), intent(in) :: model
        type(gs_fb_system_t), intent(in) :: sys
        type(gs_fb_eval_t), intent(in) :: ev
        real(dp), intent(in) :: psi(:), ds(:, :)
        real(dp), allocatable, intent(out) :: dmat(:, :)
        logical, intent(in), optional :: emulate
        real(dp), allocatable :: val(:, :), vb(:)
        integer :: k, q, i, t, nc
        logical :: vacuum

        vacuum = .false.
        if (present(emulate)) vacuum = emulate
        nc = model%n_circuits
        allocate (dmat(model%n_rows, nc), val(size(model%pts), nc), vb(nc))
        if (vacuum) then
            do q = 1, size(model%pts)
                val(q, :) = model%pts(q)%vac
            end do
            call vacuum_boundary_response(model, sys, ev, psi, vb)
        else
            !$omp parallel do default(shared) private(k, q) schedule(static)
            do k = 1, nc
                do q = 1, size(model%pts)
                    val(q, k) = point_apply(model%pts(q), ds(:, k))
                end do
                vb(k) = dot_product(ev%gb, ds(:, k))
            end do
            !$omp end parallel do
        end if
        do k = 1, nc
            do i = 1, model%n_rows
                dmat(i, k) = model%cb(i)*vb(k)
                do t = 1, 2
                    if (model%pt_of(t, i) .gt. 0) dmat(i, k) = dmat(i, k) &
                                                             + val(model%pt_of(t, i), k)
                end do
                dmat(i, k) = dmat(i, k)/model%rows(i)%scale
            end do
        end do
    end subroutine gs_ci_row_jacobian

    subroutine vacuum_boundary_response(model, sys, ev, psi, vb)
        !! Vacuum response of psi_b: the circuits' flux at the axis and the X
        !! point (cut fraction), at the limiter point of highest flux, or the
        !! mean over the target points.
        type(gs_ci_model_t), intent(in) :: model
        type(gs_fb_system_t), intent(in) :: sys
        type(gs_fb_eval_t), intent(in) :: ev
        real(dp), intent(in) :: psi(:)
        real(dp), intent(out) :: vb(:)
        real(dp), allocatable :: ga(:, :), gx(:, :)
        real(dp) :: best, v, sigma
        integer :: j, jbest

        sigma = real(sys%cfg%axis_sign, dp)
        if (sys%cfg%mode .eq. GS_FB_XPOINT) then
            allocate (ga(1, model%n_circuits), gx(1, model%n_circuits))
            call gs_circuit_flux_vectors(model%set, [ev%axis(1)], [ev%axis(2)], ga)
            call gs_circuit_flux_vectors(model%set, [ev%xpoint(1)], &
                                         [ev%xpoint(2)], gx)
            vb = (1.0_dp - sys%cfg%cut_fraction)*ga(1, :) &
                 + sys%cfg%cut_fraction*gx(1, :)
        else if (sys%cfg%mode .eq. GS_FB_LIMITER) then
            jbest = 1
            best = -huge(1.0_dp)
            do j = 1, size(model%lim)
                v = sigma*point_apply(model%lim(j), psi)
                if (v .gt. best) then
                    best = v
                    jbest = j
                end if
            end do
            vb = model%lim(jbest)%vac
        else
            vb = 0.0_dp
            do j = 1, size(model%lim)
                vb = vb + model%lim(j)%vac/real(size(model%lim), dp)
            end do
        end if
    end subroutine vacuum_boundary_response

    ! ------------------------------------------------------------------
    ! Forward solve, objective and gradient.
    ! ------------------------------------------------------------------

    subroutine forward_light(state, sys, cur, psi, ev, status)
        !! Warm-started free-boundary Newton from psi; ev without Jacobian.
        type(gs_state_t), intent(inout) :: state
        type(gs_fb_system_t), intent(inout) :: sys
        real(dp), intent(in) :: cur(:)
        real(dp), intent(inout) :: psi(:)
        type(gs_fb_eval_t), intent(inout) :: ev
        integer, intent(out) :: status
        integer :: stat

        status = GS_CI_FORWARD_FAILED
        state%psi = psi
        call gs_fb_newton(state, sys, cur, stat)
        if (stat .ne. GS_FB_OK) return
        psi = state%psi
        call gs_fb_evaluate(state, sys, psi, cur, ev, stat, jacobian=.false.)
        if (stat .ne. GS_FB_OK) return
        status = GS_CI_OK
    end subroutine forward_light

    subroutine scaled_residual(model, psi, psib, r)
        type(gs_ci_model_t), intent(in) :: model
        real(dp), intent(in) :: psi(:), psib
        real(dp), intent(out) :: r(:)
        real(dp), allocatable :: f(:)

        allocate (f(model%n_rows))
        call gs_ci_row_values(model, psi, psib, f)
        r = (f - model%rows%target)/model%rows%scale
    end subroutine scaled_residual

    function reference(cfg, n) result(ref)
        type(gs_ci_config_t), intent(in) :: cfg
        integer, intent(in) :: n
        real(dp) :: ref(n)

        ref = 0.0_dp
        if (allocated(cfg%i_ref)) ref = cfg%i_ref
    end function reference

    function objective_value(model, cfg, cur, r, phi_data, phi_reg) result(phi)
        type(gs_ci_model_t), intent(in) :: model
        type(gs_ci_config_t), intent(in) :: cfg
        real(dp), intent(in) :: cur(:), r(:)
        real(dp), intent(out), optional :: phi_data, phi_reg
        real(dp) :: phi, pd, pr

        pd = 0.5_dp*dot_product(r, r)
        pr = 0.5_dp*cfg%lambda**2*sum(((cur - reference(cfg, size(cur)))/ &
                                       model%scale)**2)
        phi = pd + pr
        if (present(phi_data)) phi_data = pd
        if (present(phi_reg)) phi_reg = pr
    end function objective_value

    subroutine evaluate_point(state, sys, model, cfg, cur, psi, fac, r, dmat, &
                              grad, phi, status, phi_data, phi_reg)
        !! Forward solve at `cur` (warm start psi), factorisation at the solution
        !! (kept in fac), residual, Jacobian and gradient of Phi.
        type(gs_state_t), intent(inout) :: state
        type(gs_fb_system_t), intent(inout) :: sys
        type(gs_ci_model_t), intent(in) :: model
        type(gs_ci_config_t), intent(in) :: cfg
        real(dp), intent(in) :: cur(:)
        real(dp), intent(inout) :: psi(:)
        type(gs_fb_factor_t), intent(inout) :: fac
        real(dp), allocatable, intent(out) :: r(:), dmat(:, :), grad(:)
        real(dp), intent(out) :: phi
        integer, intent(out) :: status
        real(dp), intent(out), optional :: phi_data, phi_reg
        type(gs_fb_eval_t) :: ev
        real(dp), allocatable :: ds(:, :)
        integer :: stat

        status = GS_CI_FORWARD_FAILED
        phi = 0.0_dp
        call forward_light(state, sys, cur, psi, ev, stat)
        if (stat .ne. GS_CI_OK) return
        call gs_fb_evaluate(state, sys, psi, cur, ev, stat, jacobian=.true.)
        if (stat .ne. GS_FB_OK) return
        call gs_fb_factor(state, sys, ev, fac, stat)
        if (stat .ne. GS_FB_OK) return
        allocate (r(model%n_rows))
        call scaled_residual(model, psi, ev%psi_bound, r)
        if (cfg%emulate_freegs) then
            allocate (ds(0, 0))
        else
            call gs_ci_flux_sensitivities(sys, ev, fac, ds, stat)
            if (stat .ne. GS_CI_OK) then
                status = GS_CI_SENSITIVITY_FAILED
                return
            end if
        end if
        call gs_ci_row_jacobian(model, sys, ev, psi, ds, dmat, cfg%emulate_freegs)
        allocate (grad(model%n_circuits))
        grad = matmul(transpose(dmat), r) + cfg%lambda**2*(cur - &
                                               reference(cfg, size(cur)))/model%scale**2
        phi = objective_value(model, cfg, cur, r, phi_data, phi_reg)
        status = GS_CI_OK
    end subroutine evaluate_point

    subroutine gs_ci_objective(state, sys, model, cfg, cur, psi, phi, grad, status, &
                               r, dmat)
        !! Objective Phi(I) after a free-boundary solve at the circuit currents
        !! `cur` (psi: warm start in, solution out) and its gradient
        !! dPhi/dI = D^T r + lambda^2 L^2 (I - I_ref), from the exact forward
        !! sensitivities (vacuum columns with cfg%emulate_freegs).
        type(gs_state_t), intent(inout) :: state
        type(gs_fb_system_t), intent(inout) :: sys
        type(gs_ci_model_t), intent(in) :: model
        type(gs_ci_config_t), intent(in) :: cfg
        real(dp), intent(in) :: cur(:)
        real(dp), intent(inout) :: psi(:)
        real(dp), intent(out) :: phi
        real(dp), allocatable, intent(out) :: grad(:)
        integer, intent(out) :: status
        real(dp), allocatable, intent(out), optional :: r(:), dmat(:, :)
        type(gs_fb_factor_t) :: fac
        real(dp), allocatable :: rr(:), dd(:, :)

        status = GS_CI_BAD_INPUT
        if (size(cur) .ne. model%n_circuits .or. size(psi) .ne. sys%n) return
        call evaluate_point(state, sys, model, cfg, cur, psi, fac, rr, dd, grad, &
                            phi, status)
        call gs_fb_release(fac)
        if (status .ne. GS_CI_OK) return
        if (present(r)) call move_alloc(rr, r)
        if (present(dmat)) call move_alloc(dd, dmat)
    end subroutine gs_ci_objective

    ! ------------------------------------------------------------------
    ! Gauss-Newton with BVLS.
    ! ------------------------------------------------------------------

    subroutine gn_step(model, cfg, dmat, r, cur, lo, up, delta, active, info)
        !! Bounded, tied, regularised linear step; delta in A.
        type(gs_ci_model_t), intent(in) :: model
        type(gs_ci_config_t), intent(in) :: cfg
        real(dp), intent(in) :: dmat(:, :), r(:), cur(:), lo(:), up(:)
        real(dp), intent(out) :: delta(:)
        integer, intent(out) :: active(:), info
        real(dp), allocatable :: a(:, :), cs(:, :), ds(:), x(:), xref(:)
        integer :: m, n, k

        m = size(dmat, 1)
        n = size(dmat, 2)
        allocate (a(m, n), x(n), xref(n))
        do k = 1, n
            a(:, k) = dmat(:, k)*model%scale(k)
        end do
        xref = (reference(cfg, n) - cur)/model%scale
        if (allocated(cfg%tie_c)) then
            allocate (cs(size(cfg%tie_c, 1), n), ds(size(cfg%tie_c, 1)))
            do k = 1, n
                cs(:, k) = cfg%tie_c(:, k)*model%scale(k)
            end do
            ds = cfg%tie_d - matmul(cfg%tie_c, cur)
            call bounded_least_squares(a, -r, (lo - cur)/model%scale, &
                                       (up - cur)/model%scale, x, info, &
                                       active=active, c=cs, d=ds, &
                                       lambda=cfg%lambda, x_ref=xref)
        else
            call bounded_least_squares(a, -r, (lo - cur)/model%scale, &
                                       (up - cur)/model%scale, x, info, &
                                       active=active, lambda=cfg%lambda, x_ref=xref)
        end if
        delta = x*model%scale
    end subroutine gn_step

    subroutine project_feasible(model, cfg, lo, up, cur, info)
        !! Nearest (in the L norm) point to `cur` with bounds and ties.
        type(gs_ci_model_t), intent(in) :: model
        type(gs_ci_config_t), intent(in) :: cfg
        real(dp), intent(in) :: lo(:), up(:)
        real(dp), intent(inout) :: cur(:)
        integer, intent(out) :: info
        real(dp), allocatable :: a(:, :), cs(:, :), ds(:), x(:), b(:)
        integer :: n, k

        n = size(cur)
        info = 0
        allocate (a(n, n), x(n), b(n))
        a = 0.0_dp
        do k = 1, n
            a(k, k) = 1.0_dp
        end do
        b = 0.0_dp
        if (allocated(cfg%tie_c)) then
            allocate (cs(size(cfg%tie_c, 1), n), ds(size(cfg%tie_c, 1)))
            do k = 1, n
                cs(:, k) = cfg%tie_c(:, k)*model%scale(k)
            end do
            ds = cfg%tie_d - matmul(cfg%tie_c, cur)
            call bounded_least_squares(a, b, (lo - cur)/model%scale, &
                                       (up - cur)/model%scale, x, info, c=cs, d=ds)
        else
            call bounded_least_squares(a, b, (lo - cur)/model%scale, &
                                       (up - cur)/model%scale, x, info)
        end if
        if (info .ne. BVLS_OK) return
        cur = min(max(cur + x*model%scale, lo), up)
    end subroutine project_feasible

    subroutine gs_ci_solve(state, sys, model, cfg, cur0, psi0, res, status, lower, &
                           upper)
        !! Reduced Gauss-Newton (or the FreeGS fixed-point iteration with
        !! cfg%emulate_freegs) from the circuit currents cur0 and the warm start
        !! psi0. Bounds are model%lower/upper when cfg%use_bounds, else none; the
        !! optional lower/upper replace them (pinned circuits have lower = upper).
        type(gs_state_t), intent(inout) :: state
        type(gs_fb_system_t), intent(inout) :: sys
        type(gs_ci_model_t), intent(in) :: model
        type(gs_ci_config_t), intent(in) :: cfg
        real(dp), intent(in) :: cur0(:), psi0(:)
        type(gs_ci_result_t), intent(out) :: res
        integer, intent(out) :: status
        real(dp), intent(in), optional :: lower(:), upper(:)
        type(gs_fb_factor_t) :: fac
        type(gs_fb_eval_t) :: ev_try
        real(dp), allocatable :: cur(:), psi(:), r(:), dmat(:, :), grad(:), delta(:)
        real(dp), allocatable :: lo(:), up(:), trial(:), psi_try(:), r_try(:), hist(:)
        real(dp) :: phi, phi_try, slope, alpha, axis_save(2), x_save(2), step
        integer, allocatable :: active(:)
        integer :: n, it, ls, info, stat
        logical :: converged

        status = GS_CI_BAD_INPUT
        n = model%n_circuits
        if (size(cur0) .ne. n .or. size(psi0) .ne. sys%n) return
        allocate (cur(n), psi(sys%n), lo(n), up(n), delta(n), active(n), &
                  trial(n), psi_try(sys%n), r_try(model%n_rows), &
                  hist(0:cfg%max_iterations))
        cur = cur0
        psi = psi0
        lo = -ieee_value(1.0_dp, ieee_positive_inf)
        up = ieee_value(1.0_dp, ieee_positive_inf)
        if (cfg%use_bounds) then
            lo = model%lower
            up = model%upper
        end if
        if (present(lower)) lo = lower
        if (present(upper)) up = upper
        if (cfg%use_bounds .or. allocated(cfg%tie_c) .or. present(lower)) then
            call project_feasible(model, cfg, lo, up, cur, info)
            if (info .ne. BVLS_OK) then
                status = GS_CI_BVLS_FAILED
                return
            end if
        end if
        hist = 0.0_dp
        converged = .false.
        status = GS_CI_NOT_CONVERGED
        active = 0
        do it = 0, cfg%max_iterations
            call evaluate_point(state, sys, model, cfg, cur, psi, fac, r, dmat, grad, &
                                phi, stat, res%phi_data, res%phi_reg)
            if (stat .ne. GS_CI_OK) then
                status = stat
                exit
            end if
            hist(it) = phi
            res%iterations = it
            call gn_step(model, cfg, dmat, r, cur, lo, up, delta, active, info)
            if (info .ne. BVLS_OK) then
                status = GS_CI_BVLS_FAILED
                exit
            end if
            step = maxval(abs(delta/model%scale))
            if (step .le. cfg%step_tolerance) then
                converged = .true.
                status = GS_CI_OK
                exit
            end if
            if (it .eq. cfg%max_iterations) exit
            axis_save = sys%axis_hint
            x_save = sys%xpoint_hint
            if (cfg%emulate_freegs) then
                alpha = cfg%relaxation
                trial = min(max(cur + alpha*delta, lo), up)
                cur = trial
                cycle
            end if
            slope = dot_product(grad, delta)
            if (slope .ge. 0.0_dp) then
                converged = -slope .le. cfg%function_tolerance*(1.0_dp + phi)
                status = merge(GS_CI_OK, GS_CI_LINE_SEARCH, converged)
                exit
            end if
            if (-slope .le. cfg%function_tolerance*(1.0_dp + phi)) then
                converged = .true.
                status = GS_CI_OK
                exit
            end if
            alpha = 1.0_dp
            status = GS_CI_LINE_SEARCH
            do ls = 0, 20
                trial = min(max(cur + alpha*delta, lo), up)
                psi_try = psi
                call forward_light(state, sys, trial, psi_try, ev_try, stat)
                if (stat .eq. GS_CI_OK) then
                    call scaled_residual(model, psi_try, ev_try%psi_bound, r_try)
                    phi_try = objective_value(model, cfg, trial, r_try)
                    if (phi_try .le. phi + cfg%armijo*alpha*slope) then
                        status = GS_CI_NOT_CONVERGED
                        exit
                    end if
                end if
                alpha = 0.5_dp*alpha
                sys%axis_hint = axis_save
                sys%xpoint_hint = x_save
            end do
            if (status .eq. GS_CI_LINE_SEARCH) exit
            cur = trial
            psi = psi_try
        end do
        if (status .eq. GS_CI_OK .or. status .eq. GS_CI_NOT_CONVERGED) then
            res%phi = phi
            call move_alloc(r, res%residual)
            call move_alloc(grad, res%gradient)
            call gs_ci_gsvd_report(dmat, model%scale, res%gsvd, stat)
        else if (allocated(r)) then
            res%phi = phi
            call move_alloc(r, res%residual)
            call move_alloc(grad, res%gradient)
        end if
        call gs_fb_release(fac)
        res%status = status
        res%currents = cur
        res%psi = psi
        res%history = hist(0:res%iterations)
        res%active = active
        if (.not. converged .and. status .eq. GS_CI_NOT_CONVERGED) res%status = status
    end subroutine gs_ci_solve

    subroutine gs_ci_fit(state, sys, model, cfg, variant, cur0, psi0, res, status)
        !! The three fit variants of DESIGN 3.3 / the project request:
        !! GS_CI_UNCONSTRAINED (no bounds, no ties), GS_CI_BOUNDED (limits and
        !! ties), GS_CI_FIXED_AT_LIMIT (ties kept; refit unconstrained, pin the
        !! circuits beyond a limit at that limit, refit the others, repeat).
        type(gs_state_t), intent(inout) :: state
        type(gs_fb_system_t), intent(inout) :: sys
        type(gs_ci_model_t), intent(in) :: model
        type(gs_ci_config_t), intent(in) :: cfg
        integer, intent(in) :: variant
        real(dp), intent(in) :: cur0(:), psi0(:)
        type(gs_ci_result_t), intent(out) :: res
        integer, intent(out) :: status
        type(gs_ci_config_t) :: c
        real(dp), allocatable :: lo(:), up(:), cur(:), psi(:)
        logical, allocatable :: pinned(:)
        integer :: n, pass, k
        logical :: violated

        status = GS_CI_BAD_INPUT
        n = model%n_circuits
        c = cfg
        select case (variant)
        case (GS_CI_UNCONSTRAINED)
            c%use_bounds = .false.
            if (allocated(c%tie_c)) deallocate (c%tie_c)
            if (allocated(c%tie_d)) deallocate (c%tie_d)
            call gs_ci_solve(state, sys, model, c, cur0, psi0, res, status)
            allocate (res%pinned(n))
            res%pinned = .false.
        case (GS_CI_BOUNDED)
            c%use_bounds = .true.
            call gs_ci_solve(state, sys, model, c, cur0, psi0, res, status)
            allocate (res%pinned(n))
            res%pinned = .false.
        case (GS_CI_FIXED_AT_LIMIT)
            c%use_bounds = .false.
            call gs_ci_solve(state, sys, model, c, cur0, psi0, res, status)
            allocate (pinned(n), lo(n), up(n))
            pinned = .false.
            c%use_bounds = .true.
            do pass = 1, n + 1
                if (status .ne. GS_CI_OK .and. status .ne. GS_CI_NOT_CONVERGED) exit
                violated = .false.
                do k = 1, n
                    if (pinned(k)) cycle
                    if (res%currents(k) .gt. model%upper(k) .or. &
                        res%currents(k) .lt. model%lower(k)) then
                        pinned(k) = .true.
                        violated = .true.
                    end if
                end do
                if (.not. violated) exit
                lo = -ieee_value(1.0_dp, ieee_positive_inf)
                up = ieee_value(1.0_dp, ieee_positive_inf)
                cur = res%currents
                psi = res%psi
                do k = 1, n
                    if (.not. pinned(k)) cycle
                    cur(k) = min(max(cur(k), model%lower(k)), model%upper(k))
                    lo(k) = cur(k)
                    up(k) = cur(k)
                end do
                call gs_ci_solve(state, sys, model, c, cur, psi, res, status, lo, up)
                res%passes = pass
            end do
            call move_alloc(pinned, res%pinned)
        case default
            return
        end select
    end subroutine gs_ci_fit

    ! ------------------------------------------------------------------
    ! GSVD report.
    ! ------------------------------------------------------------------

    subroutine gs_ci_gsvd_report(dmat, scale, rep, status)
        !! GSVD of (D, L), L = diag(1 / scale): spectrum in ascending order, the
        !! coil combination of each mode normalised to |L v| = 1 (so |D v| =
        !! gsv), and the GSVD factor X with delta = sum_j (x_j . delta) v_j.
        real(dp), intent(in) :: dmat(:, :), scale(:)
        type(gs_ci_gsvd_t), intent(out) :: rep
        integer, intent(out) :: status
        real(dp), allocatable :: b(:, :), alpha(:), beta(:), gsv(:), u(:, :), v(:, :)
        real(dp), allocatable :: xm(:, :), cm(:, :), sm(:, :), xt(:, :), vm(:, :)
        integer, allocatable :: ipiv(:), order(:)
        integer :: n, k, l, info, j, jj
        real(dp) :: s, tmp

        status = GS_CI_GSVD_FAILED
        n = size(dmat, 2)
        allocate (b(n, n))
        b = 0.0_dp
        do j = 1, n
            b(j, j) = 1.0_dp/scale(j)
        end do
        call generalized_svd(dmat, b, alpha, beta, gsv, u, v, xm, cm, sm, k, l, info)
        if (info .ne. 0 .or. k + l .ne. n .or. size(xm, 2) .ne. n) return
        allocate (order(n), xt(n, n), vm(n, n), ipiv(n))
        do j = 1, n
            order(j) = j
        end do
        do j = 2, n          ! insertion sort, ascending gsv
            jj = j
            do while (jj .gt. 1)
                if (gsv(order(jj - 1)) .le. gsv(order(jj))) exit
                info = order(jj)
                order(jj) = order(jj - 1)
                order(jj - 1) = info
                jj = jj - 1
            end do
        end do
        xt = transpose(xm)
        vm = 0.0_dp
        do j = 1, n
            vm(j, j) = 1.0_dp
        end do
        call dgesv(n, n, xt, n, ipiv, vm, n, info)    ! vm = X^-T
        if (info .ne. 0) return
        rep%n = n
        allocate (rep%gsv(n), rep%modes(n, n), rep%xmat(n, n))
        do j = 1, n
            rep%gsv(j) = gsv(order(j))
            rep%modes(:, j) = vm(:, order(j))
            rep%xmat(:, j) = xm(:, order(j))
            s = norm2(rep%modes(:, j)/scale)
            tmp = s
            rep%modes(:, j) = rep%modes(:, j)/tmp
            rep%xmat(:, j) = rep%xmat(:, j)*tmp
        end do
        rep%null_mode = 1
        rep%condition = rep%gsv(n)/max(rep%gsv(1), tiny(1.0_dp))
        status = GS_CI_OK
    end subroutine gs_ci_gsvd_report

    function gs_ci_mode_coefficients(rep, delta) result(c)
        !! Coefficients of a current change `delta` (A) on the GSVD modes.
        type(gs_ci_gsvd_t), intent(in) :: rep
        real(dp), intent(in) :: delta(:)
        real(dp) :: c(rep%n)

        c = 0.0_dp
        if (.not. allocated(rep%xmat)) return
        c = matmul(transpose(rep%xmat), delta)
    end function gs_ci_mode_coefficients

end module kin6d_gs_coil_inverse
