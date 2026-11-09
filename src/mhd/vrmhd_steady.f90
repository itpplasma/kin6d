module milca_vrmhd_steady
    !! Driven stationary states of the invariant section by pseudo-transient
    !! continuation and Newton: (shift M + R_U) du = -R with the shift going
    !! to zero as the residual falls (switched evolution relaxation).
    !! See docs/vrmhd-cylinder.md#stationary-solve.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, &
        ieee_quiet_nan
    use milca_vrmhd_parameters, only: vrmhd_parameters_t
    use milca_vrmhd_section, only: vrmhd_section_t, vrmhd_section_jacobian
    use milca_vrmhd_stationary, only: vrmhd_stationary_t, vrmhd_stationary_residual, &
        vrmhd_stationary_jacobian
    use fortsparse, only: csc_t, csc_matvec, sparse_solver_t, sparse_factor, &
        sparse_solve, sparse_free, sparse_refactor, fortsparse_status_t, status_ok, &
        sparse_factor_stats, sparse_factor_stats_t, &
        FORTSPARSE_LU_SUPERNODAL
    use fortnum_krylov, only: real_linear_operator_t, real_fgmres
    use fortnum_eigen, only: krylov_schur, eigen_options_t, eigen_stats_t, &
        eigen_real_operator_dp_t, EIGEN_LARGEST_MAGNITUDE
    implicit none
    private

    type, public :: vrmhd_steady_options_t
        !> Initial pseudo time step; zero starts with pure Newton.
        real(dp) :: dtau = 0
        real(dp) :: dtau_max = 1.0e30_dp, growth = 10, increase = 10
        !> Convergence of the maximum-norm residual, relative and absolute.
        real(dp) :: rtol = 1.0e-10_dp, atol = 0
        integer :: max_iterations = 60
        !> Largest accepted relative change of either temperature per step.
        real(dp) :: temperature_change = 0.5_dp
        !> Halvings of a step that does not reduce the residual norm.
        integer :: line_search = 12
        !> Relative residual of the FGMRES solve of each linearised step.
        real(dp) :: linear_tolerance = 1.0e-10_dp, forcing = 1.0e-4_dp
        !> Relative field change of a Newton step accepted as converged, the
        !> floor of the field scale, and the pseudo step counted as Newton.
        real(dp) :: step_tolerance = 1.0e-8_dp, step_floor = 1.0e-30_dp
        !> Velocity scale of the step test (the Alfven speed in its units).
        real(dp) :: velocity_scale = 1
        real(dp) :: dtau_newton = 1.0e8_dp
        logical :: verbose = .false.
    end type

    type, public :: vrmhd_steady_report_t
        integer :: iterations = 0, factorizations = 0, krylov = 0, status = 1
        real(dp) :: residual0 = 0, residual = 0, dtau = 0, step = 0
    end type

    public :: vrmhd_steady_solve, vrmhd_steady_spectrum

    !> Assembled Jacobian and its static-pivot factor for FGMRES: the factor
    !> is not refined against the matrix, so perturbed pivots do not force the
    !> serial scalar fallback; the outer Krylov solve restores accuracy.
    type, extends(real_linear_operator_t) :: matrix_op_t
        type(csc_t), pointer :: a => null()
    contains
        procedure :: apply => matrix_apply
    end type
    type, extends(real_linear_operator_t) :: factor_op_t
        type(sparse_solver_t), pointer :: solver => null()
        integer :: failures = 0
    contains
        procedure :: apply => factor_apply
    end type

    !> (R_U + sigma M)^(-1) M with the factor and mass it captures.
    type, extends(eigen_real_operator_dp_t) :: shift_invert_t
        type(csc_t) :: mass
        type(sparse_solver_t) :: solver
        logical, allocatable :: massless(:)
        integer :: ndof = 0
    contains
        procedure :: apply => shift_invert_apply
    end type
contains
    subroutine vrmhd_steady_solve(op, steady, p, boundary, source, load, u, &
            options, report)
        type(vrmhd_section_t), intent(in) :: op
        type(vrmhd_stationary_t), intent(in) :: steady
        type(vrmhd_parameters_t), intent(in) :: p
        real(dp), intent(in) :: boundary(op%ndof)
        real(dp), intent(in) :: source(6, op%fs%nq, op%fs%mm%topology%n_triangles)
        real(dp), intent(in) :: load(op%ndof)
        real(dp), intent(inout) :: u(steady%ndof)
        type(vrmhd_steady_options_t), intent(in) :: options
        type(vrmhd_steady_report_t), intent(out) :: report
        type(fortsparse_status_t) :: st
        real(dp), allocatable :: r(:), jv(:), du(:), step(:), rhs(:), trial(:)
        real(dp) :: rnorm, previous, dtau, shift, scale, ratio, change, floor
        logical :: reused
        type(sparse_factor_stats_t) :: factor_stats
        type(csc_t), target :: jac
        type(sparse_solver_t), target :: solver
        type(matrix_op_t) :: matrix
        type(factor_op_t) :: factor
        integer :: info, krylov
        real(dp) :: linear_residual
        integer :: it, status, lo, hi, s, ls
        allocate (r(steady%ndof), jv(steady%ndof), du(steady%ndof), &
            step(steady%ndof))
        du = 0
        dtau = options%dtau
        previous = -1
        solver%options%lu_method = FORTSPARSE_LU_SUPERNODAL
        solver%options%verify_factor = .false.
        do it = 0, options%max_iterations
            call vrmhd_stationary_residual(op, steady, p, u, du, boundary, source, &
                load, r, jv, status)
            if (status /= 0) exit
            rnorm = maxval(abs(r))
            if (.not. ieee_is_finite(rnorm)) exit
            if (it == 0) report%residual0 = rnorm
            report%iterations = it; report%residual = rnorm; report%dtau = dtau
            if (options%verbose) print '(a,i4,2es13.5)', &
                'steady iteration, residual, dtau: ', it, rnorm, dtau
            if (options%verbose) print '(a,9es9.1)', '  rows ', &
                (maxval(abs(r(op%offset(s) + 1:op%offset(s + 1)))), s=1, 8), &
                maxval(abs(r(op%ndof + 1:)))
            if (rnorm <= max(options%rtol*report%residual0, options%atol)) then
                report%status = 0
                exit
            end if
            if (it == options%max_iterations) exit
            ! Grow the pseudo step unless the residual rises by more than
            ! `increase` (near the roundoff floor it fluctuates), else halve it;
            ! slow viscous and resistive relaxation needs steps beyond 1/rate.
            if (dtau > 0 .and. previous > 0) then
                if (rnorm < options%increase*previous) then
                    dtau = min(options%dtau_max, dtau*options%growth)
                else
                    dtau = dtau/2
                end if
            end if
            shift = 0
            if (dtau > 0) shift = 1/dtau
            call vrmhd_stationary_jacobian(op, steady, p, u, boundary, source, jac, &
                status, shift=shift)
            if (status /= 0) exit
            ! The pattern is fixed, so ordering and matching are analysed once.
            if (report%factorizations > 0) then
                call sparse_refactor(solver, jac, st, reused)
            else
                call sparse_factor(solver, jac, st)
            end if
            if (.not. status_ok(st)) exit
            report%factorizations = report%factorizations + 1
            if (options%verbose) then
                call sparse_factor_stats(solver, factor_stats, st)
                if (factor_stats%lu_fallback) print '(a,a)', '  factor fallback: ', &
                    trim(factor_stats%lu_fallback_reason)
            end if
            ! FGMRES on the scaled system; its tolerance is absolute below one.
            rhs = -r/rnorm
            step = 0
            matrix%a => jac; factor%solver => solver
            call real_fgmres(matrix, factor, rhs, step, &
                options%linear_tolerance, 200, 50, info, krylov, linear_residual)
            if (options%verbose) print '(a,2i5,es11.3,i8)', &
                '  fgmres info, iterations, residual; perturbed pivots: ', info, &
                krylov, linear_residual, factor_stats%perturbed_pivots
            ! Inexact Newton: a stagnated solve is accepted below the forcing term;
            ! the roundoff of an ill-conditioned factor bounds it near 1e-8.
            if (factor%failures > 0) exit
            if (info /= 0 .and. .not. linear_residual <= options%forcing) exit
            step = rnorm*step
            report%krylov = report%krylov + krylov
            ! Damp steps that would change a temperature by more than the bound.
            scale = 1
            do s = 5, 6
                lo = op%offset(s) + 1; hi = op%offset(s + 1)
                ratio = maxval(abs(step(lo:hi))/max(abs(u(lo:hi)), tiny(1.0_dp)))
                if (ratio > options%temperature_change) &
                    scale = min(scale, options%temperature_change/ratio)
            end do
            ! The Ampere rows reach their roundoff floor first, so convergence
            ! is also accepted on a Newton step that no longer changes the fields.
            if (shift <= 1/options%dtau_newton) then
                change = 0
                do s = 1, 8
                    lo = op%offset(s) + 1; hi = op%offset(s + 1)
                    floor = options%step_floor
                    if (s >= 2 .and. s <= 4) floor = options%velocity_scale
                    change = max(change, maxval(abs(step(lo:hi))) &
                        /max(maxval(abs(u(lo:hi))), floor))
                end do
                report%step = change
                if (change <= options%step_tolerance) then
                    u = u + step
                    report%status = 0
                    exit
                end if
            end if
            if (options%verbose) print '(a,9es9.1,f6.2)', '  step ', &
                (maxval(abs(step(op%offset(s) + 1:op%offset(s + 1)))), s=1, 8), &
                maxval(abs(step(op%ndof + 1:))), scale
            ! Backtrack on the residual norm; weakly determined directions (density
            ! at low beta) make full steps leave the region of the linear model.
            ! Pseudo-transient steps may raise the residual on their way.
            do ls = 1, merge(options%line_search, 0, shift <= 0)
                trial = u + scale*step
                call vrmhd_stationary_residual(op, steady, p, trial, du, boundary, &
                    source, load, r, jv, status)
                if (status == 0) then
                    if (maxval(abs(r)) < rnorm) exit
                end if
                scale = scale/2
            end do
            u = u + scale*step
            previous = rnorm
        end do
        if (report%factorizations > 0) call sparse_free(solver)
    end subroutine

    subroutine vrmhd_steady_spectrum(op, steady, p, u, boundary, source, frozen, &
            sigma, values, vectors, residuals, info)
        !! Eigenvalues lambda of the linearisation M d_t x = -R_U x about u
        !! nearest the real shift sigma, by Krylov-Schur on the shift-invert
        !! operator (R_U + sigma M)^(-1) M; lambda = sigma - 1/theta. Fixed and
        !! `frozen` coefficients and the gauge multiplier carry no mass.
        type(vrmhd_section_t), intent(in) :: op
        type(vrmhd_stationary_t), intent(in) :: steady
        type(vrmhd_parameters_t), intent(in) :: p
        real(dp), intent(in) :: u(steady%ndof), boundary(op%ndof)
        real(dp), intent(in) :: source(6, op%fs%nq, op%fs%mm%topology%n_triangles)
        logical, intent(in) :: frozen(op%ndof)
        real(dp), intent(in) :: sigma
        complex(dp), intent(out) :: values(:), vectors(:, :)
        real(dp), intent(out) :: residuals(:)
        integer, intent(out) :: info
        type(vrmhd_section_t) :: held
        type(csc_t) :: shifted
        type(shift_invert_t) :: transform
        type(fortsparse_status_t) :: st
        type(eigen_options_t) :: options
        type(eigen_stats_t) :: stats
        real(dp), allocatable :: start(:), zero(:)
        complex(dp), allocatable :: theta(:)
        integer :: status, i
        info = -1
        held = op
        held%fixed = op%fixed .or. frozen
        call vrmhd_stationary_jacobian(held, steady, p, u, boundary, source, shifted, &
            status, shift=sigma)
        if (status /= 0) return
        allocate (zero(op%ndof), source=0.0_dp)
        call vrmhd_section_jacobian(held, p, u(:op%ndof), zero, source, 1.0_dp, &
            transform%mass, status, mass_only=.true.)
        if (status /= 0) return
        allocate (transform%massless(steady%ndof), source=.true.)
        transform%massless(:op%ndof) = held%fixed
        transform%ndof = op%ndof
        transform%solver%options%lu_method = FORTSPARSE_LU_SUPERNODAL
        call sparse_factor(transform%solver, shifted, st)
        if (.not. status_ok(st)) return
        allocate (start(steady%ndof), theta(size(values)))
        do i = 1, steady%ndof
            start(i) = sin(1.7_dp*i) + 0.3_dp*cos(0.37_dp*i*i)
        end do
        where (transform%massless) start = 0
        options%max_basis = max(40, 4*size(values))
        options%retain = max(2*size(values), 20)
        options%max_restarts = 300
        options%which = EIGEN_LARGEST_MAGNITUDE
        call krylov_schur(transform, start, theta, vectors, residuals, 1.0e-10_dp, &
            options, info, stats)
        call sparse_free(transform%solver)
        values = sigma - 1/theta
    end subroutine

    subroutine shift_invert_apply(self, x, y)
        class(shift_invert_t), intent(inout) :: self
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: y(:)
        real(dp), allocatable :: mx(:)
        type(fortsparse_status_t) :: st
        allocate (mx(size(x)), source=0.0_dp)
        mx(:self%ndof) = csc_matvec(self%mass, merge(0.0_dp, x(:self%ndof), &
            self%massless(:self%ndof)))
        where (self%massless) mx = 0
        call sparse_solve(self%solver, mx, y, st)
        ! A failed solve must not pass as an eigenvector image.
        if (.not. status_ok(st)) y = ieee_value(1.0_dp, ieee_quiet_nan)
    end subroutine

    subroutine matrix_apply(self, input, output)
        class(matrix_op_t), intent(inout) :: self
        real(dp), intent(in) :: input(:)
        real(dp), intent(out) :: output(:)
        output = csc_matvec(self%a, input)
    end subroutine

    subroutine factor_apply(self, input, output)
        class(factor_op_t), intent(inout) :: self
        real(dp), intent(in) :: input(:)
        real(dp), intent(out) :: output(:)
        real(dp), allocatable :: b(:)
        type(fortsparse_status_t) :: st
        b = input
        call sparse_solve(self%solver, b, output, st)
        if (.not. status_ok(st)) then
            self%failures = self%failures + 1
            output = input
        end if
    end subroutine
end module milca_vrmhd_steady
