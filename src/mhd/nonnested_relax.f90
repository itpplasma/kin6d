module kin6d_nonnested_relax
    ! Metriplectic relaxation on the NN1 FEEC Fourier complex, derived in
    ! docs/metriplectic-relaxation.md. The state is the full real-Fourier
    ! V1 potential A (boundary tangential trace fixed), B = curl A + F/R e_phi.
    ! With J = M1^{-1} grad W and B1 = M1^{-1} grad K / 2 (V1 projections of
    ! curl* B and B), the flow is dA/dt = -G(B1) J with a metric G whose null
    ! space contains B1:
    !   RELAX_IDEAL : G = D(B1) = |B1|^2 I - B1 B1^T pointwise (magneto-friction;
    !                 field-line topology is a continuum invariant);
    !   RELAX_TAYLOR: G = I - B1 (B1, .)_M / (B1, B1)_M (global; only helicity
    !                 and flux survive, the fixed point is a Taylor state).
    ! Implicit midpoint keeps the quadratic helicity and the linear fluxes
    ! exactly and decreases the energy by 2 dt (J, G J)_M. Generic sparse and
    ! Krylov pieces stay in FortSparse/FortNum; the pointwise metric is the
    ! FortSym leaf friction_drift (gen/gen_metriplectic.f90).
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_meridian_feec_n, only: meridian_qpoint, meridian_load
    use fortfem_meridian_fourier_curl, only: meridian_fourier_curl, &
        meridian_fourier_pairing, meridian_fourier_diagonal_mass
    use fortsparse, only: csc_t, csc_matvec, csc_matmul, csc_transpose, &
        csc_from_triplet, fortsparse_status_t, status_ok
    use fortnum_krylov, only: real_fgmres_operator, KRYLOV_OK, KRYLOV_STAGNATION
    use kin6d_equilibrium_nonnested, only: free_b_field_t, nn_factor_t, nn_factor, &
        nn_factor_solve_real, nn_factor_free, nn_divergence
    use kin6d_nonnested_mapped, only: nn_map_t, map_metric_t, METRIC_IDENTITY, &
        nn_harmonic_forms
    use kin6d_metriplectic_leaves, only: friction_drift
    implicit none
    private
    public :: relax_system_t, relax_options_t, relax_report_t, relax_state_t
    public :: relax_prepare, relax_free, relax_run, relax_state, relax_store
    public :: relax_pack_modes
    public :: RELAX_IDEAL, RELAX_TAYLOR
    integer, parameter :: RELAX_IDEAL = 1, RELAX_TAYLOR = 2
    real(dp), parameter :: pi = acos(-1.0_dp)

    type :: relax_system_t
        integer :: modes = -1, nb = 0, n1 = 0, n2 = 0, nn = 0, n0 = 0
        integer :: nl1 = 0, nq = 0, nt = 0, nphi = 0
        type(csc_t) :: curl, curl_t, m2, pair, stiff, mass1
        type(nn_factor_t) :: mass
        real(dp), allocatable :: weight(:), ck(:), ch(:), flux_load(:)
        real(dp), allocatable :: basis(:, :, :, :), wq(:, :), trig(:, :)
        integer, allocatable :: dof(:, :), sgn(:, :)
        logical, allocatable :: fixed(:), trace(:)
        real(dp) :: hh = 0, harmonic_flux = 0, mu0 = 1, length = 1, volume = 0
        logical :: ready = .false.
    end type

    type :: relax_options_t
        integer :: kind = RELAX_IDEAL
        integer :: max_steps = 400, max_outer = 40, gmres_max = 200
        ! Time steps cycle geometrically from dt to dt_max (see relax_run).
        real(dp) :: dt = 1.0e-3_dp, dt_max = 1.0e1_dp, growth = 1.5_dp
        real(dp) :: outer_tol = 1.0e-10_dp, gmres_tol = 1.0e-12_dp
        ! Hand-off to Newton when the dimensionless force falls below this.
        real(dp) :: handoff = 1.0e-3_dp
        logical :: verbose = .false.
    end type

    type :: relax_state_t
        ! Physical energy W = int |B|^2/(2 mu0), helicity int A.B, toroidal
        ! flux, dimensionless force (hand-off measure) and Taylor multiplier.
        real(dp) :: energy = 0, helicity = 0, flux = 0, force = 0, lambda = 0
        real(dp), allocatable :: trace(:)
    end type

    type :: relax_report_t
        integer :: status = 1, steps = 0, rejected = 0, outer = 0, krylov = 0
        logical :: handoff = .false.
        type(relax_state_t) :: initial, final
        real(dp) :: helicity_drift = 0, flux_drift = 0, trace_drift = 0
        real(dp) :: energy_increase = 0, identity_defect = 0, dt = 0
        real(dp) :: law_defect = 0
        ! Largest relative size of the scalar helicity projection per step.
        real(dp) :: projection = 0
        character(96) :: message = 'not run'
    end type
contains
    subroutine relax_prepare(field, system, status, nphi)
        ! Call relax_free on a prepared system before preparing it again.
        type(free_b_field_t), intent(in), target :: field
        type(relax_system_t), intent(out) :: system
        integer, intent(out) :: status
        integer, optional, intent(in) :: nphi
        type(nn_map_t), target :: identity
        type(map_metric_t) :: metric
        type(fortsparse_status_t) :: fst
        type(csc_t) :: tmp
        real(dp), allocatable :: g(:)
        complex(dp), allocatable :: load(:)
        real(dp) :: henergy, hvolume, x(2)
        integer :: b, t, q, st, ip, n
        status = 1
        if (.not. allocated(field%b)) return
        system%modes = ubound(field%b, 2); system%nb = 2*system%modes + 1
        system%n1 = field%fs%sp(1)%ndof; system%n2 = field%fs%sp(2)%ndof
        system%nn = field%fs%nn; system%n0 = field%fs%n0
        system%nl1 = field%fs%sp(1)%nloc; system%nq = field%fs%nq
        system%nt = field%fs%mm%topology%n_triangles
        system%mu0 = field%mu0; system%length = field%radial_scale
        ! Products of three fields of degree <= modes are integrated exactly
        ! against a degree-modes test function for nphi >= 4 modes + 1.
        system%nphi = 4*system%modes + 1
        if (present(nphi)) system%nphi = nphi
        if (system%nphi < 2*system%modes + 1) return
        call meridian_fourier_curl(field%fs, system%modes, system%curl, st)
        if (st /= 0) return
        call csc_transpose(system%curl, system%curl_t, fst)
        if (.not. status_ok(fst)) return
        call meridian_fourier_diagonal_mass(field%fs, 2, system%modes, system%m2, st)
        if (st /= 0) return
        call meridian_fourier_diagonal_mass(field%fs, 1, system%modes, &
            system%mass1, st)
        if (st /= 0) return
        call meridian_fourier_pairing(field%fs, system%modes, system%pair, st)
        if (st /= 0) return
        call csc_matmul(system%m2, system%curl, tmp, fst)
        if (.not. status_ok(fst)) return
        call csc_matmul(system%curl_t, tmp, system%stiff, fst)
        if (.not. status_ok(fst)) return
        call nn_factor(field%fs%mass(1), field%fs%sp(1)%bnd, system%mass, st)
        if (st /= 0) return
        allocate (system%weight(0:system%nb - 1), &
            system%fixed(system%nb*system%n1), system%trace(system%nb*system%n2))
        system%weight = 0.5_dp; system%weight(0) = 1
        do b = 0, system%nb - 1
            system%fixed(b*system%n1 + 1:(b + 1)*system%n1) = field%fs%sp(1)%bnd
            system%trace(b*system%n2 + 1:(b + 1)*system%n2) = field%fs%sp(2)%bnd
        end do
        ! Harmonic h = F/R e_phi: <V2, h>, <h, h>, int dR dZ / R (flux).
        metric%fs => field%fs
        identity%modes = 0
        allocate (identity%d(field%fs%n0, 2, 0:0)); identity%d = 0
        metric%map => identity; metric%kind = METRIC_IDENTITY
        call nn_harmonic_forms(field, metric, max(system%nphi, 1), g, henergy, &
            hvolume, st)
        if (st /= 0) return
        system%ck = field%f_vac*csc_matvec(system%curl_t, g)
        system%hh = field%f_vac**2*henergy
        system%harmonic_flux = field%f_vac*hvolume
        call meridian_load(field%fs, 1, harmonic_unit, load, st)
        if (st /= 0) return
        allocate (system%ch(system%nb*system%n1)); system%ch = 0
        system%ch(:system%n1) = field%f_vac*real(load, dp)
        call meridian_load(field%fs, 2, harmonic_unit, load, st)
        if (st /= 0) return
        system%flux_load = real(load, dp)
        ! Cached physical V1 basis at the quadrature points (hot loop input).
        allocate (system%basis(3, system%nl1, system%nq, system%nt), &
            system%wq(system%nq, system%nt))
        system%dof = field%fs%sp(1)%dof; system%sgn = field%fs%sp(1)%sgn
        st = 0
        !$omp parallel do default(shared) private(t, q, x, n) reduction(max:st) &
        !$omp& schedule(static)
        do t = 1, system%nt
            do q = 1, system%nq
                call meridian_qpoint(field%fs, 1, t, q, x, system%wq(q, t), &
                    system%basis(:, :, q, t), n)
                st = max(st, n)
            end do
        end do
        !$omp end parallel do
        if (st /= 0) return
        system%volume = sum(system%wq)
        allocate (system%trig(0:system%nb - 1, system%nphi))
        do ip = 1, system%nphi
            system%trig(0, ip) = 1
            do n = 1, system%modes
                system%trig(2*n - 1, ip) = cos(n*2*pi*(ip - 1)/system%nphi)
                system%trig(2*n, ip) = sin(n*2*pi*(ip - 1)/system%nphi)
            end do
        end do
        if (.not. all(ieee_is_finite(system%ck))) return
        system%ready = .true.
        status = 0
    contains
        subroutine harmonic_unit(x, v)
            real(dp), intent(in) :: x(2)
            complex(dp), intent(out) :: v(3)
            real(dp) :: radius
            radius = x(1)
            if (field%fs%cyl_r0 > 0) radius = field%fs%cyl_r0
            v = 0; v(2) = 1/radius
        end subroutine
    end subroutine

    subroutine relax_free(system)
        type(relax_system_t), intent(inout) :: system
        call nn_factor_free(system%mass)
        system%ready = .false.
    end subroutine

    subroutine relax_pack_modes(system, a0, an, a)
        ! Full real layout from the NN1 per-mode potentials: a0 is the n=0
        ! scalar R A_phi (V0), an(:, n) the n>0 Nedelec potentials in the
        ! R A_phi=0 gauge, with complex amplitude (c - i s)/2 per mode.
        type(relax_system_t), intent(in) :: system
        complex(dp), intent(in) :: a0(:), an(:, :)
        real(dp), allocatable, intent(out) :: a(:)
        integer :: n
        allocate (a(system%nb*system%n1)); a = 0
        a(system%nn + 1:system%n1) = real(a0, dp)
        do n = 1, system%modes
            a((2*n - 1)*system%n1 + 1:(2*n - 1)*system%n1 + system%nn) = &
                2*real(an(:, n), dp)
            a(2*n*system%n1 + 1:2*n*system%n1 + system%nn) = -2*aimag(an(:, n))
        end do
    end subroutine

    subroutine relax_store(system, a, field, status)
        ! B = curl A into the NN1 complex Fourier coefficients of field.
        type(relax_system_t), intent(in) :: system
        real(dp), intent(in) :: a(:)
        type(free_b_field_t), intent(inout) :: field
        integer, intent(out) :: status
        real(dp), allocatable :: bh(:)
        real(dp) :: div
        integer :: n
        allocate (bh(system%nb*system%n2))
        bh = csc_matvec(system%curl, a)
        field%b(:, 0) = cmplx(bh(:system%n2), 0.0_dp, dp)
        do n = 1, system%modes
            field%b(:, n) = cmplx(bh((2*n - 1)*system%n2 + 1:2*n*system%n2), &
                -bh(2*n*system%n2 + 1:(2*n + 1)*system%n2), dp)/2
        end do
        field%psi = a(system%nn + 1:system%n1)
        call nn_divergence(field, div, status)
        if (status == 0 .and. div > 1.0e-11_dp) status = 1
    end subroutine

    subroutine relax_state(system, kind, a, state, status)
        type(relax_system_t), intent(inout) :: system
        integer, intent(in) :: kind
        real(dp), intent(in) :: a(:)
        type(relax_state_t), intent(out) :: state
        integer, intent(out) :: status
        real(dp), allocatable :: ca(:), ge(:), gh(:), j(:), b1(:), load(:)
        real(dp) :: diss, b2
        status = 1
        if (.not. system%ready .or. size(a) /= size(system%fixed)) return
        ca = csc_matvec(system%curl, a)
        ge = gradient_energy(system, a); gh = gradient_helicity(system, a)
        state%energy = pi*(dot_product(a, ge) + dot_product(a, system%ck) &
            + system%hh)/system%mu0
        state%helicity = 2*pi*(dot_product(a, gh) + dot_product(a, system%ch))
        state%flux = dot_product(system%flux_load, ca(:system%n2)) &
            + system%harmonic_flux
        state%trace = pack(ca, system%trace)
        call mass_solve(system, ge, j, status)
        if (status /= 0) return
        call mass_solve(system, gh, b1, status)
        if (status /= 0) return
        state%lambda = dot_product(j, mask(system, gh)) &
            /max(dot_product(b1, mask(system, gh)), tiny(1.0_dp))
        if (kind == RELAX_IDEAL) then
            call ideal_load(system, j, b1, load, diss, b2)
            ! |J x B1|_Q / (|B1|^2_Q / (sqrt(V) L)), dimensionless.
            state%force = sqrt(diss)*sqrt(system%volume)*system%length &
                /max(b2*system%volume, tiny(1.0_dp))
        else
            load = mask(system, ge) - state%lambda*mask(system, gh)
            state%force = sqrt(max(dot_product(j - state%lambda*b1, load), 0.0_dp)) &
                *system%length/max(sqrt(dot_product(b1, mask(system, gh))), &
                tiny(1.0_dp))
        end if
        status = 0
    end subroutine

    subroutine relax_run(system, a, options, report)
        ! Adaptive implicit-midpoint relaxation; a is updated in place and
        ! keeps its fixed (boundary) entries.
        type(relax_system_t), intent(inout) :: system
        real(dp), intent(inout) :: a(:)
        type(relax_options_t), intent(in) :: options
        type(relax_report_t), intent(out) :: report
        type(relax_state_t) :: before, after
        type(nn_factor_t) :: precond
        real(dp), allocatable :: delta(:), trial(:), b1s(:)
        real(dp) :: lambdas, dts
        real(dp) :: dt, beta, factored_dt, factored_beta, dissipation, scale
        logical :: have_factor, ok
        integer :: st, step
        report%status = 1
        if (.not. system%ready .or. size(a) /= size(system%fixed)) return
        if (options%kind /= RELAX_IDEAL .and. options%kind /= RELAX_TAYLOR) return
        if (.not. all(ieee_is_finite(a))) return
        call relax_state(system, options%kind, a, before, st)
        if (st /= 0) return
        report%initial = before
        dt = options%dt; have_factor = .false.; factored_dt = -1; factored_beta = -1
        allocate (delta(size(a))); delta = 0
        step = 0
        do while (step < options%max_steps)
            if (before%force <= options%handoff) then
                report%handoff = .true.
                exit
            end if
            call midpoint_step(dt, delta, ok, dissipation, beta)
            if (ok) then
                trial = a + delta
                call relax_state(system, options%kind, trial, after, st)
                if (st /= 0) return
                scale = max(abs(before%energy), tiny(1.0_dp))
                ok = after%energy <= before%energy + 1.0e-13_dp*scale
            end if
            if (.not. ok) then
                report%rejected = report%rejected + 1
                if (report%rejected > 40) then
                    report%message = 'step rejected repeatedly'
                    exit
                end if
                dt = dt/4; delta = 0
                cycle
            end if
            step = step + 1
            report%energy_increase = max(report%energy_increase, &
                (after%energy - before%energy)/scale)
            ! Discrete energy law W(n+1) - W(n) = -dt (J, G J) at the midpoint:
            ! relative defect on steps whose dissipation exceeds 1e-8 W, and
            ! the absolute defect in units of W on every step.
            if (dissipation > 1.0e-8_dp*scale) report%identity_defect = &
                max(report%identity_defect, &
                abs(after%energy - before%energy + dissipation)/dissipation)
            report%law_defect = max(report%law_defect, &
                abs(after%energy - before%energy + dissipation)/scale)
            if (options%verbose) print '(a,i4,4es11.3,2i7)', 'relax step ', step, dt, &
                after%energy, after%force, after%helicity, report%outer, report%krylov
            a = trial; before = after
            ! Sawtooth time steps: implicit midpoint is A- but not L-stable,
            ! so a mode with rate r is damped strongly only near r dt = 2.
            ! Cycling dt geometrically over [dt, dt_max] damps every mode in
            ! each cycle while each step stays structure preserving.
            if (dt*options%growth > options%dt_max) then
                delta = delta*(options%dt/dt); dt = options%dt
            else
                dt = dt*options%growth; delta = delta*options%growth
            end if
        end do
        call nn_factor_free(precond)
        report%steps = step; report%dt = dt
        report%final = before
        report%helicity_drift = abs(before%helicity - report%initial%helicity) &
            /max(abs(report%initial%helicity), tiny(1.0_dp))
        report%flux_drift = abs(before%flux - report%initial%flux) &
            /max(abs(report%initial%flux), tiny(1.0_dp))
        report%trace_drift = maxval(abs(before%trace - report%initial%trace)) &
            /max(maxval(abs(report%initial%trace)), tiny(1.0_dp))
        if (report%handoff) then
            report%status = 0
            report%message = 'hand-off force reached'
        else if (step >= options%max_steps) then
            report%status = 0
            report%message = 'step budget exhausted before hand-off'
        end if
    contains
        subroutine midpoint_step(dt, delta, ok, dissipation, beta)
            ! Solve M delta = dt Phi(a + delta/2), Phi = -M G J. With B1
            ! (ideal) or lambda (Taylor) frozen the problem is linear and is
            ! solved by FGMRES preconditioned with M + (dt/2) beta K. Ideal:
            ! Picard on B1. Taylor: secant on lambda for the helicity defect
            ! h = delta . M B1(midpoint), which vanishes exactly at the
            ! midpoint solution. A final scalar step along B1(midpoint)
            ! removes the remaining helicity defect to roundoff.
            real(dp), intent(in) :: dt
            real(dp), intent(inout) :: delta(:)
            logical, intent(out) :: ok
            real(dp), intent(out) :: dissipation, beta
            real(dp), allocatable :: abar(:), j(:), b1(:), ge0(:), gh0(:), j0(:)
            real(dp), allocatable :: x(:), load(:), gh(:)
            real(dp) :: lambda, diss, b2, change, l0, l1, h0, h1, hscale
            real(dp) :: qa, qb, qc, sh
            integer :: k, st
            ok = .false.; dissipation = 0; beta = 1; change = huge(1.0_dp)
            ge0 = gradient_energy(system, a); gh0 = gradient_helicity(system, a)
            call mass_solve(system, ge0, j0, st)
            if (st /= 0) return
            abar = a + delta/2
            call frozen(abar, j, b1, lambda, diss, b2, st)
            if (st /= 0) return
            if (options%kind == RELAX_IDEAL) beta = max(b2, tiny(1.0_dp))
            if (.not. have_factor .or. abs(dt - factored_dt) > 0 .or. &
                abs(beta - factored_beta) > 0.3_dp*factored_beta) then
                call nn_factor_free(precond)
                call preconditioner(dt/2*beta, st)
                if (st /= 0) return
                have_factor = .true.; factored_dt = dt; factored_beta = beta
            end if
            if (options%kind == RELAX_IDEAL) then
                do k = 1, options%max_outer
                    b1s = b1
                    call ideal_load(system, j0, b1s, load, diss, b2)
                    call linear_solve(dt*load, delta, x, st)
                    if (st /= 0) return
                    change = norm2(x - delta)/max(norm2(x), tiny(1.0_dp))
                    delta = x
                    report%outer = report%outer + 1
                    if (change <= options%outer_tol) exit
                    call frozen(a + delta/2, j, b1, lambda, diss, b2, st)
                    if (st /= 0) return
                end do
                if (change > options%outer_tol) then
                    report%message = 'ideal midpoint iteration did not converge'
                    return
                end if
            else
                hscale = max(abs(dot_product(a, gh0)), tiny(1.0_dp))
                l0 = lambda
                call taylor_defect(l0, delta, h0, st)
                if (st /= 0) return
                call frozen(a + delta/2, j, b1, l1, diss, b2, st)
                if (st /= 0) return
                do k = 1, options%max_outer
                    if (abs(l1 - l0) <= 1.0e-15_dp*max(1.0_dp, abs(l1))) exit
                    call taylor_defect(l1, delta, h1, st)
                    if (st /= 0) return
                    report%outer = report%outer + 1
                    if (abs(h1) <= 1.0e-15_dp*hscale) exit
                    if (abs(h1 - h0) <= tiny(1.0_dp)) exit
                    lambda = l1 - h1*(l1 - l0)/(h1 - h0)
                    l0 = l1; h0 = h1; l1 = lambda
                end do
            end if
            ! Helicity projection: K(a + delta + s B1bar) = K(a), smaller root.
            abar = a + delta/2
            call frozen(abar, j, b1, lambda, diss, b2, st)
            if (st /= 0) return
            gh = mask(system, gradient_helicity(system, abar))
            qc = 2*dot_product(delta, gh)
            qb = 2*dot_product(b1, mask(system, gradient_helicity(system, &
                a + delta)))
            qa = dot_product(b1, mask(system, csc_matvec(system%pair, &
                csc_matvec(system%curl, b1))))
            if (abs(qb) <= tiny(1.0_dp)) return
            if (qb*qb - 4*qa*qc < 0) return
            sh = -2*qc/(qb + sign(sqrt(qb*qb - 4*qa*qc), qb))
            report%projection = max(report%projection, &
                abs(sh)*norm2(b1)/max(norm2(delta), tiny(1.0_dp)))
            delta = delta + sh*b1
            abar = a + delta/2
            call frozen(abar, j, b1, lambda, diss, b2, st)
            if (st /= 0) return
            if (options%kind == RELAX_IDEAL) then
                dissipation = 2*pi*dt*diss/system%mu0
            else
                load = mask(system, gradient_energy(system, abar)) &
                    - lambda*mask(system, gradient_helicity(system, abar))
                dissipation = 2*pi*dt*dot_product(j - lambda*b1, load)/system%mu0
            end if
            ok = all(ieee_is_finite(delta))
        end subroutine

        subroutine taylor_defect(lambda, delta, h, status)
            ! delta(lambda) of the frozen-lambda midpoint system and the
            ! helicity defect h = delta . M B1(a + delta/2).
            real(dp), intent(in) :: lambda
            real(dp), intent(inout) :: delta(:)
            real(dp), intent(out) :: h
            integer, intent(out) :: status
            real(dp), allocatable :: x(:)
            lambdas = lambda
            call linear_solve(-dt*(mask(system, gradient_energy(system, a)) &
                - lambda*mask(system, gradient_helicity(system, a))), delta, x, status)
            if (status /= 0) return
            delta = x
            h = dot_product(delta, mask(system, gradient_helicity(system, &
                a + delta/2)))
        end subroutine

        subroutine linear_solve(rhs, guess, x, status)
            real(dp), intent(in) :: rhs(:), guess(:)
            real(dp), allocatable, intent(out) :: x(:)
            integer, intent(out) :: status
            real(dp) :: rnorm
            integer :: info, iterations
            status = 1
            x = guess
            dts = dt
            if (norm2(rhs) <= tiny(1.0_dp)) then
                x = 0
            else
                call real_fgmres_operator(apply, precondition, rhs, x, &
                    options%gmres_tol, options%gmres_max, min(options%gmres_max, 60), &
                    info, iterations, rnorm)
                report%krylov = report%krylov + iterations
                if (info /= KRYLOV_OK .and. info /= KRYLOV_STAGNATION) then
                    report%message = 'FGMRES failed in the midpoint solve'
                    return
                end if
                if (.not. all(ieee_is_finite(x))) return
            end if
            where (system%fixed) x = 0
            status = 0
        end subroutine

        subroutine frozen(abar, j, b1, lambda, diss, b2, status)
            real(dp), intent(in) :: abar(:)
            real(dp), allocatable, intent(out) :: j(:), b1(:)
            real(dp), intent(out) :: lambda, diss, b2
            integer, intent(out) :: status
            real(dp), allocatable :: gh(:), load(:)
            call mass_solve(system, gradient_energy(system, abar), j, status)
            if (status /= 0) return
            gh = mask(system, gradient_helicity(system, abar))
            call mass_solve(system, gh, b1, status)
            if (status /= 0) return
            lambda = dot_product(j, gh)/max(dot_product(b1, gh), tiny(1.0_dp))
            diss = 0; b2 = 1
            if (options%kind == RELAX_IDEAL) &
                call ideal_load(system, j, b1, load, diss, b2)
        end subroutine

        subroutine apply(input, output)
            ! (M - dt/2 L) x with L the frozen linear part of M Phi.
            real(dp), intent(in) :: input(:)
            real(dp), intent(out) :: output(:)
            real(dp), allocatable :: x(:), jx(:), load(:)
            real(dp) :: diss, b2
            integer :: st
            x = input
            where (system%fixed) x = 0
            output = mass_apply(system, x)
            if (options%kind == RELAX_IDEAL) then
                call mass_solve(system, csc_matvec(system%stiff, x), jx, st)
                call ideal_load(system, jx, b1s, load, diss, b2)
                output = output - dts/2*load
            else
                output = output + dts/2*(mask(system, csc_matvec(system%stiff, x)) &
                    - lambdas*mask(system, csc_matvec(system%pair, &
                    csc_matvec(system%curl, x))))
            end if
            where (system%fixed) output = input
        end subroutine

        subroutine precondition(input, output)
            real(dp), intent(in) :: input(:)
            real(dp), intent(out) :: output(:)
            integer :: st
            call nn_factor_solve_real(precond, input, output, st)
            where (system%fixed) output = input
        end subroutine

        subroutine preconditioner(s, status)
            real(dp), intent(in) :: s
            integer, intent(out) :: status
            type(csc_t) :: p
            type(fortsparse_status_t) :: fst
            integer :: nm, nk
            nm = system%mass1%nnz; nk = system%stiff%nnz
            call csc_from_triplet(system%mass1%nrow, system%mass1%ncol, &
                [csc_rows(system%mass1), csc_rows(system%stiff)], &
                [csc_cols(system%mass1), csc_cols(system%stiff)], &
                [system%mass1%val(:nm), s*system%stiff%val(:nk)], p, fst)
            status = 1
            if (.not. status_ok(fst)) return
            call nn_factor(p, system%fixed, precond, status)
        end subroutine

    end subroutine

    subroutine ideal_load(system, j, b1, load, diss, b2)
        ! load_i = < phi_i . (J x B1) x B1 >, diss = <|J x B1|^2>,
        ! b2 = <|B1|^2>/V on the cached quadrature, nphi-point phi mean.
        type(relax_system_t), intent(in) :: system
        real(dp), intent(in) :: j(:), b1(:)
        real(dp), allocatable, intent(out) :: load(:)
        real(dp), intent(out) :: diss, b2
        real(dp), allocatable :: local(:, :, :)
        real(dp) :: jb(3, 0:system%nb - 1), bb(3, 0:system%nb - 1)
        real(dp) :: eb(3, 0:system%nb - 1), cj(system%nl1), cb(system%nl1)
        real(dp) :: jx(3), bx(3), e(3), d, w, tr(0:system%nb - 1)
        integer :: t, q, ip, b, nb, n1
        nb = system%nb; n1 = system%n1
        allocate (local(system%nl1, 0:nb - 1, system%nt))
        diss = 0; b2 = 0
        !$omp parallel do default(shared) &
        !$omp& private(t, q, ip, b, jb, bb, eb, cj, cb, jx, bx, e, d, w, tr) &
        !$omp& reduction(+:diss, b2) schedule(static)
        do t = 1, system%nt
            local(:, :, t) = 0
            do q = 1, system%nq
                do b = 0, nb - 1
                    cj = system%sgn(:, t)*j(b*n1 + system%dof(:, t))
                    cb = system%sgn(:, t)*b1(b*n1 + system%dof(:, t))
                    jb(:, b) = matmul(system%basis(:, :, q, t), cj)
                    bb(:, b) = matmul(system%basis(:, :, q, t), cb)
                end do
                eb = 0
                w = system%wq(q, t)/system%nphi
                do ip = 1, system%nphi
                    tr = system%trig(:, ip)
                    jx = matmul(jb, tr); bx = matmul(bb, tr)
                    call friction_drift(jx(1), jx(2), jx(3), bx(1), bx(2), bx(3), &
                        e(1), e(2), e(3), d)
                    diss = diss + w*d
                    b2 = b2 + w*sum(bx**2)
                    do b = 0, nb - 1
                        eb(:, b) = eb(:, b) + (w*tr(b))*e
                    end do
                end do
                do b = 0, nb - 1
                    local(:, b, t) = local(:, b, t) &
                        + matmul(eb(:, b), system%basis(:, :, q, t))
                end do
            end do
            do b = 0, nb - 1
                local(:, b, t) = local(:, b, t)*system%sgn(:, t)
            end do
        end do
        !$omp end parallel do
        allocate (load(nb*n1)); load = 0
        do t = 1, system%nt
            do b = 0, nb - 1
                load(b*n1 + system%dof(:, t)) = load(b*n1 + system%dof(:, t)) &
                    + local(:, b, t)
            end do
        end do
        where (system%fixed) load = 0
        b2 = b2/max(system%volume, tiny(1.0_dp))
    end subroutine

    function gradient_energy(system, a) result(g)
        ! K a + C^T <V2, h>: half the gradient of B.M2 B.
        type(relax_system_t), intent(in) :: system
        real(dp), intent(in) :: a(:)
        real(dp), allocatable :: g(:)
        g = csc_matvec(system%stiff, a) + system%ck
    end function

    function gradient_helicity(system, a) result(g)
        ! P (C a) + <V1, h>: half the gradient of int A.B.
        type(relax_system_t), intent(in) :: system
        real(dp), intent(in) :: a(:)
        real(dp), allocatable :: g(:)
        g = csc_matvec(system%pair, csc_matvec(system%curl, a)) + system%ch
    end function

    function mask(system, v) result(w)
        type(relax_system_t), intent(in) :: system
        real(dp), intent(in) :: v(:)
        real(dp), allocatable :: w(:)
        w = merge(0.0_dp, v, system%fixed)
    end function

    function mass_apply(system, x) result(y)
        type(relax_system_t), intent(in) :: system
        real(dp), intent(in) :: x(:)
        real(dp), allocatable :: y(:)
        y = mask(system, csc_matvec(system%mass1, x))
    end function

    subroutine mass_solve(system, r, x, status)
        ! Interior V1 Fourier mass: one meridian factor, block weights 1, 1/2.
        type(relax_system_t), intent(inout) :: system
        real(dp), intent(in) :: r(:)
        real(dp), allocatable, intent(out) :: x(:)
        integer, intent(out) :: status
        integer :: b, n1
        n1 = system%n1
        allocate (x(size(r))); x = 0
        status = 0
        do b = 0, system%nb - 1
            call nn_factor_solve_real(system%mass, r(b*n1 + 1:(b + 1)*n1), &
                x(b*n1 + 1:(b + 1)*n1), status)
            if (status /= 0) return
            x(b*n1 + 1:(b + 1)*n1) = x(b*n1 + 1:(b + 1)*n1)/system%weight(b)
        end do
    end subroutine

    function csc_rows(m) result(r)
        type(csc_t), intent(in) :: m
        integer, allocatable :: r(:)
        r = m%row_idx(:m%nnz)
    end function

    function csc_cols(m) result(c)
        type(csc_t), intent(in) :: m
        integer, allocatable :: c(:)
        integer :: j
        allocate (c(m%nnz))
        do j = 1, m%ncol
            c(m%col_ptr(j):m%col_ptr(j + 1) - 1) = j
        end do
    end function
end module kin6d_nonnested_relax
