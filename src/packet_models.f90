! Q2: reduced wave-packet models of the Q1 two-electron Coulomb encounter,
! measured against the Q1 partial-wave reference (same physical parameters,
! observables and bins) with a full-Hamiltonian a-posteriori estimator.
!
! Models (one approximation changed at a time; docs/packet-models.md):
!   M2   product of two electron Gaussians (WPMD, canonical, FortSym leaves);
!        the pair Coulomb term erf(q/sqrt(2(s1+s2)))/q ties the CM width to the
!        relative width (alpha_cm = 4 alpha_rel for the symmetric preparation)
!   M2S/M2T  the same electron Gaussians symmetrized (singlet) or
!        antisymmetrized (triplet): two-determinant TDVP with exact overlap and
!        exchange integrals (metric engine, parity tie + CM tie)
!   M3   one thawed Gaussian for the relative coordinate under g/r (canonical);
!        the CM is the exact free Gaussian
!   M3S/M3T  M3 symmetrized/antisymmetrized without the CM tie (metric engine)
!   M4   coherent sum of K relative Gaussians (variational multi-configuration
!        Gaussians, all cross terms, complex amplitudes; K = 1 equals M3);
!        zero-population Gaussians follow their own one-Gaussian TDVP until
!        the variational equations give them weight (basis enrichment).
! The classical Wigner ensemble (M1) is Q1's encounter_classical_ensemble.
!
! Estimator: the model trajectory is projected onto the Q1 basis, a_m = P psi,
! and r_h = i hbar d a_m/dt - H_h a_m is evaluated with the Q1 banded
! Hamiltonian. Since P is linear and time independent, a_m solves the discrete
! equation with source r_h, so ||a_m(t) - a_ref(t)|| <= ||a_m(0) - a_ref(0)||
! + (1/hbar) int ||r_h|| is a bound in the represented space (binary64, time
! integral by trapezoid of sampled ||r_h||). The continuum statement adds the
! reference's discretization error and the model projection defect (both
! reported, neither bounded).
module kin6d_packet_models
    use, intrinsic :: iso_fortran_env, only: dp => real64, i64 => int64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use kin6d_gaussian_pair, only: gaussian_t, pair_rule_t, pair_rule_init, pair_elements, &
                                   gaussian_value, gaussian_norm2
    use kin6d_packet_models_leaves
    use kin6d_coulomb_encounter, only: encounter_t, encounter_state, encounter_project, &
                                       encounter_angular_bins, encounter_norm, in_sector, &
                                       SECTOR_PRODUCT, SECTOR_SINGLET, SECTOR_TRIPLET
    use kin6d_radial_dvr, only: radial_band_hamiltonian, radial_band_matvec
    use fortnum_ode_cash_karp, only: cash_karp_step
    implicit none
    private
    public :: packet_spec_t, packet_sample_t, packet_trace_t
    public :: packet_samples, packet_trace, packet_project, packet_bin_jvp
    public :: model_name, model_sector
    public :: MODEL_M2, MODEL_M2S, MODEL_M2T, MODEL_M3, MODEL_M3S, MODEL_M3T, MODEL_M4

    integer, parameter :: MODEL_M2 = 2, MODEL_M2S = 21, MODEL_M2T = 22, MODEL_M3 = 3, &
                          MODEL_M3S = 31, MODEL_M3T = 32, MODEL_M4 = 4
    real(dp), parameter :: pi = acos(-1.0_dp)
    complex(dp), parameter :: iu = (0.0_dp, 1.0_dp)

    type :: packet_spec_t
        integer :: model = MODEL_M3
        integer :: k = 1                     ! M4 basis size (1..8)
        real(dp) :: kappa = 0                ! harmonic test potential kappa r^2/2
        ! Tikhonov eps relative to the largest equilibrated metric eigenvalue;
        ! <= 0 selects 1e-6 for guided Gaussian sums (M4) and 1e-10 otherwise.
        real(dp) :: regularization = 0
        real(dp) :: rtol = 1.0e-10_dp, atol = 1.0e-12_dp
        real(dp) :: refuse_condition = 1.0e14_dp  ! overlap condition refusal limit
        integer :: boys_nodes = 96
        logical :: metric_engine = .false.   ! M2/M3 through the metric engine (cross-check)
    end type packet_spec_t

    ! One model state at a sample time: normalized relative Gaussian list with
    ! parameter rates, CM data and conditioning status.
    type :: packet_sample_t
        real(dp) :: t = 0
        type(gaussian_t), allocatable :: g(:), rate(:)
        complex(dp) :: cm_overlap = 1        ! <free CM|model CM>
        real(dp) :: cm_residual = 0          ! ||r_C|| (free CM Hamiltonian)
        real(dp) :: energy = 0, norm = 0     ! model expectation values (full model)
        real(dp) :: condition = 1            ! scaled metric eigenvalue ratio
        real(dp) :: overlap_condition = 1    ! normalized Gaussian overlap matrix
        integer :: regularized = 0           ! metric directions below eps
        real(dp) :: width_asymmetry = 0      ! M2: |s1 - s2|/(s1 + s2)
    end type packet_sample_t

    ! Diagnostics of a model trajectory against the Q1 reference.
    type :: packet_trace_t
        real(dp), allocatable :: t(:), error(:), phase_error(:), residual(:), bound(:)
        real(dp), allocatable :: total_error(:), total_bound(:), energy(:), norm(:)
        real(dp), allocatable :: projected_norm(:), dvr_energy(:)
        real(dp) :: initial_error = 0, max_condition = 1, max_overlap_condition = 1
        integer :: max_regularized = 0
        logical :: refused = .false.
        character(:), allocatable :: status
        complex(dp), allocatable :: a(:, :, :)   ! projected model state at the last time
        real(dp) :: seconds_model = 0, seconds_projection = 0
        integer(i64) :: rhs_evaluations = 0
    end type packet_trace_t

    ! Metric (vMCG) engine for Gaussian sums, optional parity and CM ties.
    type :: engine_t
        integer :: k = 1, nred = 8, parity = 0   ! parity: 0 none, +1 singlet, -1 triplet
        logical :: cm_tie = .false., guided = .false.
        real(dp) :: hbar = 1, mu = 1, g = 0, kap = 0, mass_cm = 2, eps = 1.0e-10_dp
        real(dp), allocatable :: jmap(:, :), scale(:)
        type(pair_rule_t) :: rule
    end type engine_t

    type :: canonical_t
        integer :: model = MODEL_M3
        real(dp) :: hbar = 1, mu = 1, m = 1, g = 0, kap = 0
    end type canonical_t

    interface
        subroutine dsyev(jobz, uplo, n, a, lda, w, work, lwork, info)
            import :: dp
            character, intent(in) :: jobz, uplo
            integer, intent(in) :: n, lda, lwork
            real(dp), intent(inout) :: a(lda, *)
            real(dp), intent(out) :: w(*), work(*)
            integer, intent(out) :: info
        end subroutine dsyev
    end interface
contains
    function model_name(spec) result(name)
        type(packet_spec_t), intent(in) :: spec
        character(:), allocatable :: name
        character(8) :: buf
        select case (spec%model)
        case (MODEL_M2); name = 'M2'
        case (MODEL_M2S); name = 'M2S'
        case (MODEL_M2T); name = 'M2T'
        case (MODEL_M3); name = 'M3'
        case (MODEL_M3S); name = 'M3S'
        case (MODEL_M3T); name = 'M3T'
        case default
            write (buf, '(a,i0)') 'M4K', spec%k
            name = trim(buf)
        end select
    end function model_name

    pure integer function model_sector(spec) result(sector)
        type(packet_spec_t), intent(in) :: spec
        select case (spec%model)
        case (MODEL_M2S, MODEL_M3S); sector = SECTOR_SINGLET
        case (MODEL_M2T, MODEL_M3T); sector = SECTOR_TRIPLET
        case default; sector = SECTOR_PRODUCT
        end select
    end function model_sector

    ! Adaptive Cash-Karp integration of y' = rhs(t, y) from t0 to t1.
    subroutine integrate(rhs, ctx, y, t0, t1, rtol, atol, nfev, ok)
        interface
            subroutine rhs(t, y, dydt, ctx)
                import :: dp
                real(dp), intent(in) :: t, y(:)
                real(dp), intent(out) :: dydt(:)
                class(*), intent(in), optional :: ctx
            end subroutine rhs
        end interface
        class(*), intent(in) :: ctx
        real(dp), intent(inout) :: y(:)
        real(dp), intent(in) :: t0, t1, rtol, atol
        integer(i64), intent(inout) :: nfev
        logical, intent(out) :: ok
        real(dp), dimension(size(y)) :: k1, k2, k3, k4, k5, k6, ytmp, y5, yerr
        real(dp) :: t, h, err
        integer :: nf, it
        logical :: have
        ok = .true.
        if (t1 <= t0) return
        t = t0; h = (t1 - t0)/16; have = .false.; nf = 0
        do it = 1, 10000000
            h = min(h, t1 - t)
            call cash_karp_step(rhs, t, y, h, have, k1, k2, k3, k4, k5, k6, ytmp, y5, &
                                yerr, nf, ctx)
            err = maxval(abs(yerr)/(atol + rtol*max(abs(y), abs(y5))))
            if (.not. ieee_is_finite(err)) err = huge(1.0_dp)
            if (err <= 1) then
                t = t + h; y = y5; have = .false.
                if (t >= t1*(1 - 4*epsilon(1.0_dp)) .or. t >= t1) exit
                h = h*min(4.0_dp, 0.9_dp*err**(-0.2_dp))
            else
                have = .true.
                h = h*max(0.1_dp, 0.9_dp*err**(-0.25_dp))
            end if
            if (h < 1.0e-14_dp*max(abs(t1), 1.0_dp)) then
                ok = .false.; exit
            end if
        end do
        nfev = nfev + nf
    end subroutine integrate

    ! ---- Canonical thawed Gaussians (FortSym leaves) -------------------------
    ! M3: y = (qx, qz, px, pz, s, ps, gamma[, tangents 6 each]).
    ! M2: y = (q1x, q1z, q2x, q2z, p1x, p1z, p2x, p2z, s1, s2, ps1, ps2,
    !          gamma_total, gamma_cm).
    subroutine canonical_rhs(t, y, dydt, ctx)
        real(dp), intent(in) :: t, y(:)
        real(dp), intent(out) :: dydt(:)
        class(*), intent(in), optional :: ctx
        real(dp) :: e, sc, psc, ec, qd(2)
        integer :: j, o
        if (.false.) dydt(1) = t
        select type (c => ctx)
        type is (canonical_t)
            if (c%model == MODEL_M3) then
                call m3_rates(c, y(1:6), dydt(1:6), e, dydt(7))
                do j = 1, (size(y) - 7)/6
                    o = 7 + 6*(j - 1)
                    call m3_tangent(c, y(1:6), y(o + 1:o + 6), dydt(o + 1:o + 6))
                end do
            else
                call m2_rates(c, y(1:12), dydt(1:12), e, dydt(13))
                sc = (y(9) + y(10))/4; psc = y(11) + y(12)
                call pm_gaussian_energy(y(5) + y(7), y(6) + y(8), sc, psc, 2*c%m, c%hbar, ec)
                qd = (dydt(1:2) + dydt(3:4))/2
                call pm_df_phase(y(5) + y(7), y(6) + y(8), qd(1), qd(2), sc, &
                                 dydt(11) + dydt(12), ec, c%hbar, dydt(14))
            end if
        end select
    end subroutine canonical_rhs

    ! Gaussian-averaged Coulomb switches to its series form for T = q^2/(2s) < 0.02
    ! (truncation below 1e-16 relative).
    subroutine m3_rates(c, x, r, e, dgamma)
        type(canonical_t), intent(in) :: c
        real(dp), intent(in) :: x(6)
        real(dp), intent(out) :: r(6), e, dgamma
        if ((x(1)**2 + x(2)**2)/(2*x(5)) < 0.02_dp) then
            call pm_m3_rhs_series(x(1), x(2), x(3), x(4), x(5), x(6), c%g, c%kap, c%mu, &
                                  c%hbar, r(1), r(2), r(3), r(4), r(5), r(6), e, dgamma)
        else
            call pm_m3_rhs(x(1), x(2), x(3), x(4), x(5), x(6), c%g, c%kap, c%mu, c%hbar, &
                           r(1), r(2), r(3), r(4), r(5), r(6), e, dgamma)
        end if
    end subroutine m3_rates

    subroutine m3_tangent(c, x, dx, dr)
        type(canonical_t), intent(in) :: c
        real(dp), intent(in) :: x(6), dx(6)
        real(dp), intent(out) :: dr(6)
        real(dp) :: de, dgam
        if ((x(1)**2 + x(2)**2)/(2*x(5)) < 0.02_dp) then
            call pm_m3_rhs_series_jvp(x(1), x(2), x(3), x(4), x(5), x(6), c%g, c%kap, c%mu, &
                                      c%hbar, dx(1), dx(2), dx(3), dx(4), dx(5), dx(6), &
                                      0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, dr(1), dr(2), dr(3), &
                                      dr(4), dr(5), dr(6), de, dgam)
        else
            call pm_m3_rhs_jvp(x(1), x(2), x(3), x(4), x(5), x(6), c%g, c%kap, c%mu, c%hbar, &
                               dx(1), dx(2), dx(3), dx(4), dx(5), dx(6), 0.0_dp, 0.0_dp, &
                               0.0_dp, 0.0_dp, dr(1), dr(2), dr(3), dr(4), dr(5), dr(6), &
                               de, dgam)
        end if
    end subroutine m3_tangent

    subroutine m2_rates(c, x, r, e, dgamma)
        type(canonical_t), intent(in) :: c
        real(dp), intent(in) :: x(12)
        real(dp), intent(out) :: r(12), e, dgamma
        if (((x(1) - x(3))**2 + (x(2) - x(4))**2)/(2*(x(9) + x(10))) < 0.02_dp) then
            call pm_m2_rhs_series(x(1), x(2), x(3), x(4), x(5), x(6), x(7), x(8), x(9), &
                                  x(10), x(11), x(12), c%g, c%m, c%hbar, r(1), r(2), r(3), &
                                  r(4), r(5), r(6), r(7), r(8), r(9), r(10), r(11), r(12), &
                                  e, dgamma)
        else
            call pm_m2_rhs(x(1), x(2), x(3), x(4), x(5), x(6), x(7), x(8), x(9), x(10), &
                           x(11), x(12), c%g, c%m, c%hbar, r(1), r(2), r(3), r(4), r(5), &
                           r(6), r(7), r(8), r(9), r(10), r(11), r(12), e, dgamma)
        end if
    end subroutine m2_rates

    ! Normalized canonical Gaussian (q, p, s, ps, gamma) and the rates (or
    ! tangent) of its complex parameters: alpha = 1/(4s) - i ps/(3 hbar),
    ! d = (2 pi s)^(-3/4) exp(i gamma).
    subroutine canonical_gaussian(x, gam, dx, dgam, hbar, g, rate)
        real(dp), intent(in) :: x(6), gam, dx(6), dgam, hbar
        type(gaussian_t), intent(out) :: g, rate
        real(dp) :: ar, ai, dar, dai
        call pm_width_parameter_jvp(x(5), x(6), hbar, dx(5), dx(6), 0.0_dp, dar, dai)
        call pm_width_parameter(x(5), x(6), hbar, ar, ai)
        g%q = x(1:2); g%kappa = x(3:4)/hbar; g%alpha = cmplx(ar, ai, dp)
        g%d = (2*pi*x(5))**(-0.75_dp)*exp(iu*gam)
        rate%q = dx(1:2); rate%kappa = dx(3:4)/hbar; rate%alpha = cmplx(dar, dai, dp)
        rate%d = g%d*(-0.75_dp*dx(5)/x(5) + iu*dgam)
    end subroutine canonical_gaussian

    ! ---- Metric (variational multi-configuration Gaussian) engine -----------
    ! Real parameters per Gaussian: (qx, qz, px, pz, Re alpha, Im alpha, Re d,
    ! Im d); psi = sum_k d_k g_k. McLachlan form of the Dirac-Frenkel principle
    ! on this complex manifold: Re(M) xdot = Im(h)/hbar, M_ij = <d_i psi|d_j psi>,
    ! h_i = <d_i psi|H|psi>, reduced by the tie map x_full = J x_red and solved
    ! in scaled coordinates with Tikhonov eps; null directions (zero-population
    ! Gaussians) follow the guide velocity.
    subroutine engine_setup(eng, enc, spec, ngauss, parity, cm_tie, guided)
        type(engine_t), intent(out) :: eng
        type(encounter_t), intent(in) :: enc
        type(packet_spec_t), intent(in) :: spec
        integer, intent(in) :: ngauss, parity
        logical, intent(in) :: cm_tie, guided
        real(dp) :: sig, amp
        integer :: k, b
        eng%k = ngauss; eng%parity = parity; eng%cm_tie = cm_tie; eng%guided = guided
        eng%hbar = enc%cfg%hbar; eng%mu = enc%mu; eng%g = enc%cfg%g; eng%kap = spec%kappa
        eng%mass_cm = enc%mtotal; eng%eps = spec%regularization
        if (.not. (eng%eps > 0)) eng%eps = merge(1.0e-6_dp, 1.0e-10_dp, guided)
        call pair_rule_init(eng%rule, spec%boys_nodes)
        if (parity /= 0) then
            eng%nred = 8
            allocate (eng%jmap(16, 8))
            eng%jmap = 0
            do k = 1, 8
                eng%jmap(k, k) = 1
            end do
            eng%jmap(9:12, 1:4) = -reshape([1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1], &
                                           [4, 4])
            eng%jmap(13, 5) = 1; eng%jmap(14, 6) = 1
            eng%jmap(15, 7) = parity; eng%jmap(16, 8) = parity
        else
            eng%nred = 8*ngauss
            allocate (eng%jmap(eng%nred, eng%nred))
            eng%jmap = 0
            do k = 1, eng%nred
                eng%jmap(k, k) = 1
            end do
        end if
        sig = enc%cfg%sigma
        amp = (2*pi*sig**2)**(-0.75_dp)
        allocate (eng%scale(eng%nred))
        do k = 1, eng%nred/8
            b = 8*(k - 1)
            eng%scale(b + 1:b + 8) = [sig, sig, eng%hbar/sig, eng%hbar/sig, 1/sig**2, &
                                      1/sig**2, amp, amp]
        end do
    end subroutine engine_setup

    pure subroutine engine_unpack(eng, x, gs)
        type(engine_t), intent(in) :: eng
        real(dp), intent(in) :: x(:)
        type(gaussian_t), intent(out) :: gs(:)
        real(dp) :: xf(size(eng%jmap, 1))
        integer :: k, b
        xf = matmul(eng%jmap, x(1:eng%nred))
        do k = 1, size(gs)
            b = 8*(k - 1)
            gs(k)%q = xf(b + 1:b + 2); gs(k)%kappa = xf(b + 3:b + 4)/eng%hbar
            gs(k)%alpha = cmplx(xf(b + 5), xf(b + 6), dp)
            gs(k)%d = cmplx(xf(b + 7), xf(b + 8), dp)
        end do
    end subroutine engine_unpack

    subroutine engine_rhs(t, y, dydt, ctx)
        real(dp), intent(in) :: t, y(:)
        real(dp), intent(out) :: dydt(:)
        class(*), intent(in), optional :: ctx
        type(packet_sample_t) :: info
        if (.false.) dydt(1) = t
        select type (e => ctx)
        type is (engine_t)
            call engine_eval(e, y, dydt, info, .false.)
        end select
    end subroutine engine_rhs

    ! Rates of y = (x_red[, phi_cm]); info carries energy, norm and conditioning;
    ! with full = .true. also the overlap condition.
    subroutine engine_eval(eng, y, ydot, info, full)
        type(engine_t), intent(in) :: eng
        real(dp), intent(in) :: y(:)
        real(dp), intent(out) :: ydot(:)
        type(packet_sample_t), intent(inout) :: info
        logical, intent(in) :: full
        integer :: n, nf, nr, k, l, i, j, info_lapack, kmax
        type(gaussian_t) :: gs(size(eng%jmap, 1)/8), gc
        complex(dp) :: s, sig(4, 4), hel(4), dk(8, 4), dl(8, 4), nrm, er, sc, sigc(4, 4), helc(4)
        complex(dp), allocatable :: m(:, :), h(:), o(:), smat(:, :)
        complex(dp) :: kap(5:6)
        real(dp), allocatable :: mr(:, :), mred(:, :), rr(:), rred(:), lam(:), work(:)
        real(dp), allocatable :: guide(:), xdot(:), rhs(:), xf(:), emb(:, :)
        real(dp), allocatable :: scl(:)
        real(dp) :: lmax, eps, r6(6), ecan, dgam, sx, ec, dmax
        kmax = size(gs); nf = 8*kmax; nr = eng%nred; n = nf
        call engine_unpack(eng, y, gs)
        allocate (m(nf, nf), h(nf), o(nf), smat(kmax, kmax))
        m = 0; h = 0; o = 0; nrm = 0; er = 0
        do k = 1, kmax
            call derivative_map(gs(k), eng%hbar, dk)
            do l = 1, kmax
                call pair_elements(eng%rule, gs(k), gs(l), eng%hbar, eng%mu, eng%g, &
                                   eng%kap, s, sig, hel)
                smat(k, l) = s
                call derivative_map(gs(l), eng%hbar, dl)
                m(8*k - 7:8*k, 8*l - 7:8*l) = matmul(conjg(dk), matmul(s*sig, transpose(dl)))
                h(8*k - 7:8*k) = h(8*k - 7:8*k) + gs(l)%d*matmul(conjg(dk), hel)
                o(8*k - 7:8*k) = o(8*k - 7:8*k) + gs(l)%d*matmul(conjg(dk), s*sig(:, 1))
                nrm = nrm + conjg(gs(k)%d)*gs(l)%d*s
                er = er + conjg(gs(k)%d)*gs(l)%d*hel(1)
            end do
        end do
        info%norm = real(nrm); info%energy = real(er)/real(nrm)
        if (eng%cm_tie) then
            gc%q = 0; gc%kappa = 0; gc%alpha = 4*gs(1)%alpha; gc%d = 1
            call pair_elements(eng%rule, gc, gc, eng%hbar, eng%mass_cm, 0.0_dp, 0.0_dp, &
                               sc, sigc, helc)
            kap(5) = -4; kap(6) = -4*iu
            m = sc*m
            do i = 5, 6
                m(i, :) = m(i, :) + conjg(kap(i))*sc*sigc(4, 1)*conjg(o)
                m(:, i) = m(:, i) + kap(i)*sc*sigc(1, 4)*o
            end do
            do i = 5, 6
                do j = 5, 6
                    m(i, j) = m(i, j) + conjg(kap(i))*kap(j)*sc*sigc(4, 4)*nrm
                end do
            end do
            h = helc(1)*o + sc*h
            do i = 5, 6
                h(i) = h(i) + conjg(kap(i))*(helc(4)*nrm + sc*sigc(4, 1)*er)
            end do
            info%norm = real(sc*nrm)
            info%energy = real(helc(1)*nrm + sc*er)/real(sc*nrm)
        end if
        allocate (mr(nf, nf), rr(nf))
        mr = real(m); rr = aimag(h)/eng%hbar
        mred = matmul(transpose(eng%jmap), matmul(mr, eng%jmap))
        rred = matmul(transpose(eng%jmap), rr)
        ! Jacobi equilibration: sc_i = 1/sqrt(M_ii); empty (zero-population)
        ! directions keep the static parameter scale.
        allocate (scl(nr))
        dmax = maxval([(mred(j, j), j=1, nr)])
        do j = 1, nr
            scl(j) = eng%scale(j)
            if (mred(j, j) > 1.0e-14_dp*dmax) scl(j) = 1/sqrt(mred(j, j))
        end do
        do j = 1, nr
            mred(:, j) = mred(:, j)*scl*scl(j)
        end do
        rred = rred*scl
        ! Guide velocity (scaled): each Gaussian's own one-Gaussian TDVP.
        allocate (guide(nr)); guide = 0
        if (eng%guided) then
            do k = 1, kmax
                sx = 1/(4*real(gs(k)%alpha))
                call pm_m3_rhs(gs(k)%q(1), gs(k)%q(2), eng%hbar*gs(k)%kappa(1), &
                               eng%hbar*gs(k)%kappa(2), sx, -3*eng%hbar*aimag(gs(k)%alpha), &
                               eng%g, eng%kap, eng%mu, eng%hbar, r6(1), r6(2), r6(3), r6(4), &
                               r6(5), r6(6), ecan, dgam)
                guide(8*k - 7:8*k) = [r6(1:4), -r6(5)/(4*sx**2), -r6(6)/(3*eng%hbar), &
                                      0.0_dp, 0.0_dp]/scl(8*k - 7:8*k)
            end do
        end if
        allocate (lam(nr), work(max(1, 66*nr)))
        call dsyev('V', 'U', nr, mred, nr, lam, work, size(work), info_lapack)
        lmax = max(maxval(lam), tiny(1.0_dp)); eps = eng%eps*lmax
        allocate (rhs(nr), xdot(nr))
        ! mred now holds eigenvectors V; rhs = V^T (r - M guide).
        rhs = rred - matmul(mred, lam*matmul(transpose(mred), guide))
        rhs = matmul(transpose(mred), rhs)
        xdot = guide + matmul(mred, rhs*lam/(lam**2 + eps**2))
        ydot(1:nr) = xdot*scl
        info%condition = lmax/max(minval(abs(lam)), tiny(1.0_dp))
        info%regularized = count(lam < eps)
        if (eng%cm_tie) then
            xf = matmul(eng%jmap, ydot(1:nr))
            ec = real(helc(1)/sc)
            ydot(nr + 1) = (3*eng%hbar*(1/(4*real(gc%alpha)))*4*xf(6) - ec)/eng%hbar
        end if
        if (.not. full) return
        ! Normalized overlap matrix of the occupied Gaussians (real embedding).
        allocate (emb(2*kmax, 2*kmax))
        do k = 1, kmax
            do l = 1, kmax
                s = smat(k, l)/sqrt(real(smat(k, k))*real(smat(l, l)))
                emb(k, l) = real(s); emb(k + kmax, l + kmax) = real(s)
                emb(k, l + kmax) = -aimag(s); emb(k + kmax, l) = aimag(s)
            end do
        end do
        deallocate (lam, work); allocate (lam(2*kmax), work(66*2*kmax))
        call dsyev('N', 'U', 2*kmax, emb, 2*kmax, lam, work, size(work), info_lapack)
        info%overlap_condition = maxval(lam)/max(minval(lam), tiny(1.0_dp))
        if (.false.) n = info_lapack
    end subroutine engine_eval

    ! d(d_k g_k)/dx_i in the elementary basis (g, u_x g, u_z g, |u|^2 g).
    pure subroutine derivative_map(gk, hbar, dmap)
        type(gaussian_t), intent(in) :: gk
        real(dp), intent(in) :: hbar
        complex(dp), intent(out) :: dmap(8, 4)
        dmap = 0
        dmap(1, 1) = -iu*gk%kappa(1)*gk%d; dmap(1, 2) = 2*gk%alpha*gk%d
        dmap(2, 1) = -iu*gk%kappa(2)*gk%d; dmap(2, 3) = 2*gk%alpha*gk%d
        dmap(3, 2) = iu*gk%d/hbar; dmap(4, 3) = iu*gk%d/hbar
        dmap(5, 4) = -gk%d; dmap(6, 4) = -iu*gk%d
        dmap(7, 1) = 1; dmap(8, 1) = iu
    end subroutine derivative_map

    ! ---- Model trajectories at sample times ----------------------------------
    ! Preparation identical to Q1: relative Gaussian rms sigma per axis centred at
    ! (b, 0, -z0) with wave vector k z and phase exp(i k z); CM Gaussian at rest
    ! with rms sigma_cm (= sigma/2 for two matched electrons).
    subroutine packet_samples(enc, spec, times, samples, ok, message, nfev)
        type(encounter_t), intent(in) :: enc
        type(packet_spec_t), intent(in) :: spec
        real(dp), intent(in) :: times(:)
        type(packet_sample_t), allocatable, intent(out) :: samples(:)
        logical, intent(out) :: ok
        character(:), allocatable, intent(out) :: message
        integer(i64), intent(out) :: nfev
        type(canonical_t) :: can
        type(engine_t) :: eng
        real(dp), allocatable :: y(:)
        real(dp) :: t0, hb, sig, kz
        integer :: it, parity
        logical :: cm, canonical
        allocate (samples(size(times)))
        ok = .false.; message = 'invalid model specification'; nfev = 0
        hb = enc%cfg%hbar; sig = enc%cfg%sigma; kz = enc%cfg%k
        if (spec%k < 1 .or. spec%k > 8) return
        canonical = (spec%model == MODEL_M3 .or. spec%model == MODEL_M2) .and. &
                    .not. spec%metric_engine
        select case (spec%model)
        case (MODEL_M3, MODEL_M2)
            if (.not. canonical) then
                call engine_setup(eng, enc, spec, 1, 0, spec%model == MODEL_M2, .false.)
                call engine_initial(eng, enc, y)
            end if
            can%model = spec%model; can%hbar = hb; can%mu = enc%mu; can%m = enc%cfg%m1
            can%g = enc%cfg%g; can%kap = spec%kappa
            if (.not. canonical) then
                continue
            else if (spec%model == MODEL_M3) then
                y = [enc%cfg%b, -enc%cfg%z0, 0.0_dp, hb*kz, sig**2, 0.0_dp, -kz*enc%cfg%z0]
            else
                message = 'M2 needs equal masses and no harmonic term'
                if (enc%cfg%m1 /= enc%cfg%m2 .or. spec%kappa /= 0) return
                y = [enc%cfg%b/2, -enc%cfg%z0/2, -enc%cfg%b/2, enc%cfg%z0/2, 0.0_dp, hb*kz, &
                     0.0_dp, -hb*kz, sig**2/2, sig**2/2, 0.0_dp, 0.0_dp, -kz*enc%cfg%z0, 0.0_dp]
            end if
        case (MODEL_M2S, MODEL_M2T, MODEL_M3S, MODEL_M3T, MODEL_M4)
            parity = 0
            if (spec%model == MODEL_M2S .or. spec%model == MODEL_M3S) parity = 1
            if (spec%model == MODEL_M2T .or. spec%model == MODEL_M3T) parity = -1
            cm = spec%model == MODEL_M2S .or. spec%model == MODEL_M2T
            if (spec%model == MODEL_M4) then
                call engine_setup(eng, enc, spec, spec%k, 0, .false., .true.)
            else
                call engine_setup(eng, enc, spec, 2, parity, cm, .false.)
            end if
            call engine_initial(eng, enc, y)
        case default
            return
        end select
        message = 'model integration failed'
        t0 = 0
        do it = 1, size(times)
            samples(it)%t = times(it)
            if (canonical) then
                call integrate(canonical_rhs, can, y, t0, times(it), spec%rtol, spec%atol, &
                               nfev, ok)
                if (.not. ok) return
                call canonical_sample(can, enc, y, samples(it))
            else
                call integrate(engine_rhs, eng, y, t0, times(it), spec%rtol, spec%atol, &
                               nfev, ok)
                if (.not. ok) return
                call engine_sample(eng, enc, y, samples(it))
            end if
            t0 = times(it)
        end do
        ok = .true.; message = ''
    end subroutine packet_samples

    ! Initial reduced parameters (and CM phase) of the engine models.
    subroutine engine_initial(eng, enc, y)
        type(engine_t), intent(in) :: eng
        type(encounter_t), intent(in) :: enc
        real(dp), allocatable, intent(out) :: y(:)
        ! Guide offsets for zero-population Gaussians (dx, dz in sigma; dpz in
        ! the momentum width hbar/(2 sigma)); off-axis offsets only when b /= 0.
        real(dp), parameter :: axis(2, 7) = reshape([0, -2, 0, 2, -2, 0, 2, 0, 0, -4, 0, 4, &
                                                     -2, -2], [2, 7])
        real(dp), parameter :: side(3, 7) = reshape([-1, 0, 0, 1, 0, 0, 0, 0, -2, 0, 0, 2, &
                                                     -2, 0, 0, 2, 0, 0, 0, -2, 0], [3, 7])
        type(gaussian_t) :: gs(size(eng%jmap, 1)/8)
        complex(dp) :: s, sig(4, 4), hel(4), d, s12
        real(dp) :: sg, hb, kz, n2, nc, off(3)
        integer :: k, b
        sg = enc%cfg%sigma; hb = enc%cfg%hbar; kz = enc%cfg%k
        allocate (y(eng%nred + merge(1, 0, eng%cm_tie)))
        y = 0
        do k = 1, eng%nred/8
            off = 0
            if (k > 1) then
                if (abs(enc%cfg%b) > 0) then
                    off = side(:, k - 1)
                else
                    off = [0.0_dp, axis(:, k - 1)]
                end if
            end if
            b = 8*(k - 1)
            y(b + 1:b + 6) = [enc%cfg%b + off(1)*sg, -enc%cfg%z0 + off(2)*sg, 0.0_dp, &
                              hb*kz + off(3)*hb/(2*sg), 1/(4*sg**2), 0.0_dp]
        end do
        y(7) = 1; y(8) = 0
        call engine_unpack(eng, y, gs)
        call pair_elements(eng%rule, gs(1), gs(1), hb, 0.0_dp, 0.0_dp, 0.0_dp, s, sig, hel)
        n2 = real(s)
        if (eng%parity /= 0) then
            call pair_elements(eng%rule, gs(1), gs(2), hb, 0.0_dp, 0.0_dp, 0.0_dp, s12, sig, hel)
            n2 = 2*n2 + 2*real(gs(2)%d)*real(s12)
        end if
        nc = 1
        if (eng%cm_tie) nc = (pi/(2*4*real(gs(1)%alpha)))**0.75_dp
        d = exp(-iu*kz*enc%cfg%z0)/(nc*sqrt(n2))
        y(7) = real(d); y(8) = aimag(d)
        do k = 2, eng%nred/8
            y(8*k - 1:8*k) = 0
        end do
    end subroutine engine_initial

    subroutine canonical_sample(can, enc, y, smp)
        type(canonical_t), intent(in) :: can
        type(encounter_t), intent(in) :: enc
        real(dp), intent(in) :: y(:)
        type(packet_sample_t), intent(inout) :: smp
        real(dp) :: dy(size(y)), xr(6), dxr(6), xc(6), dxc(6), e, dg
        type(gaussian_t) :: gc, rc
        call canonical_rhs(0.0_dp, y, dy, can)
        allocate (smp%g(1), smp%rate(1))
        if (can%model == MODEL_M3) then
            call m3_rates(can, y(1:6), dy(1:6), e, dg)
            call canonical_gaussian(y(1:6), y(7), dy(1:6), dy(7), can%hbar, smp%g(1), &
                                    smp%rate(1))
            smp%energy = e; smp%norm = 1
            return
        end if
        call m2_rates(can, y(1:12), dy(1:12), e, dg)
        ! Relative Gaussian (q1-q2, (p1-p2)/2, s1+s2, (ps1+ps2)/4) and CM Gaussian
        ! ((q1+q2)/2, p1+p2, (s1+s2)/4, ps1+ps2) of the symmetric product state.
        xr = [y(1:2) - y(3:4), (y(5:6) - y(7:8))/2, y(9) + y(10), (y(11) + y(12))/4]
        dxr = [dy(1:2) - dy(3:4), (dy(5:6) - dy(7:8))/2, dy(9) + dy(10), (dy(11) + dy(12))/4]
        xc = [(y(1:2) + y(3:4))/2, y(5:6) + y(7:8), (y(9) + y(10))/4, y(11) + y(12)]
        dxc = [(dy(1:2) + dy(3:4))/2, dy(5:6) + dy(7:8), (dy(9) + dy(10))/4, dy(11) + dy(12)]
        call canonical_gaussian(xr, y(13) - y(14), dxr, dy(13) - dy(14), can%hbar, &
                                smp%g(1), smp%rate(1))
        call canonical_gaussian(xc, y(14), dxc, dy(14), can%hbar, gc, rc)
        call cm_diagnostics(enc, gc, xc(5), xc(6), dxc(5), dxc(6), smp)
        smp%energy = e; smp%norm = 1
        smp%width_asymmetry = abs(y(9) - y(10))/(y(9) + y(10))
    end subroutine canonical_sample

    ! Overlap of the model CM Gaussian gc (normalized, phase included) with the
    ! exact free CM Gaussian at time smp%t, and its free-Hamiltonian residual norm.
    subroutine cm_diagnostics(enc, gc, s, ps, ds, dps, smp)
        type(encounter_t), intent(in) :: enc
        type(gaussian_t), intent(in) :: gc
        real(dp), intent(in) :: s, ps, ds, dps
        type(packet_sample_t), intent(inout) :: smp
        type(gaussian_t) :: gf
        type(pair_rule_t) :: rule
        complex(dp) :: a0, tau, sv, sig(4, 4), hel(4)
        call pair_rule_init(rule, 4)
        a0 = 1/(4*enc%sigma_cm**2)
        tau = 2*iu*enc%cfg%hbar*a0*smp%t/enc%mtotal
        gf%q = 0; gf%kappa = 0; gf%alpha = a0/(1 + tau)
        gf%d = (2*a0/pi)**0.75_dp*(1 + tau)**(-1.5_dp)
        call pair_elements(rule, gf, gc, enc%cfg%hbar, 0.0_dp, 0.0_dp, 0.0_dp, sv, sig, hel)
        smp%cm_overlap = conjg(gf%d)*gc%d*sv
        call pm_free_residual(s, ps, ds, dps, enc%mtotal, enc%cfg%hbar, smp%cm_residual)
    end subroutine cm_diagnostics

    subroutine engine_sample(eng, enc, y, smp)
        type(engine_t), intent(in) :: eng
        type(encounter_t), intent(in) :: enc
        real(dp), intent(in) :: y(:)
        type(packet_sample_t), intent(inout) :: smp
        real(dp) :: dy(size(y)), xf(size(eng%jmap, 1)), dxf(size(eng%jmap, 1))
        real(dp) :: arc, aic, darc, daic, phi, dphi, nc, dnc, s
        complex(dp) :: fac, dfac
        type(gaussian_t) :: gc
        integer :: k, b, n
        call engine_eval(eng, y, dy, smp, .true.)
        n = size(eng%jmap, 1)/8
        allocate (smp%g(n), smp%rate(n))
        call engine_unpack(eng, y, smp%g)
        xf = matmul(eng%jmap, y(1:eng%nred)); dxf = matmul(eng%jmap, dy(1:eng%nred))
        do k = 1, n
            b = 8*(k - 1)
            smp%rate(k)%q = dxf(b + 1:b + 2); smp%rate(k)%kappa = dxf(b + 3:b + 4)/eng%hbar
            smp%rate(k)%alpha = cmplx(dxf(b + 5), dxf(b + 6), dp)
            smp%rate(k)%d = cmplx(dxf(b + 7), dxf(b + 8), dp)
        end do
        if (.not. eng%cm_tie) return
        ! Split Psi = C_n (x) R_n with C_n normalized and carrying the CM
        ! Dirac-Frenkel phase phi, so that r = r_C (x) R_n + C_n (x) r_R, r_C _|_ C_n.
        arc = 4*xf(5); aic = 4*xf(6); darc = 4*dxf(5); daic = 4*dxf(6)
        phi = y(eng%nred + 1); dphi = dy(eng%nred + 1)
        nc = (pi/(2*arc))**0.75_dp; dnc = -0.75_dp*nc*darc/arc
        fac = nc*exp(-iu*phi); dfac = (dnc - iu*dphi*nc)*exp(-iu*phi)
        do k = 1, n
            smp%rate(k)%d = smp%rate(k)%d*fac + smp%g(k)%d*dfac
            smp%g(k)%d = smp%g(k)%d*fac
        end do
        gc%q = 0; gc%kappa = 0; gc%alpha = cmplx(arc, aic, dp); gc%d = exp(iu*phi)/nc
        s = 1/(4*arc)
        call cm_diagnostics(enc, gc, s, -3*eng%hbar*aic, -darc/(4*arc**2), &
                            -3*eng%hbar*daic, smp)
    end subroutine engine_sample

    ! ---- Projection onto the Q1 basis and full-Hamiltonian diagnostics -------
    ! a = P psi and adot = P (d psi/dt) for psi = sum_k d_k g_k (rates in rate);
    ! layout a(n, 0:mmax, 0:lmax) as encounter_state.
    subroutine packet_project(enc, g, rate, a, adot)
        type(encounter_t), intent(in) :: enc
        type(gaussian_t), intent(in) :: g(:), rate(:)
        complex(dp), allocatable, intent(out) :: a(:, :, :), adot(:, :, :)
        complex(dp), allocatable :: a4(:, :, :, :)
        real(dp) :: amax, ext
        integer :: k, nph
        logical :: used(size(g))
        amax = 0
        do k = 1, size(g)
            used(k) = abs(g(k)%d) > 0 .or. abs(rate(k)%d) > 0
            if (.not. used(k)) cycle
            ext = sqrt(sum(g(k)%q**2)) + 12/(2*sqrt(real(g(k)%alpha)))
            amax = max(amax, abs(2*g(k)%alpha*g(k)%q(1) + iu*g(k)%kappa(1))*ext)
        end do
        nph = 2*(enc%mmax + ceiling(amax + 6*sqrt(amax)) + 8)
        call encounter_project(enc, enc%lmax, enc%mmax, nph, 2, field, a4)
        allocate (a(enc%mesh%n, 0:enc%mmax, 0:enc%lmax), &
                  adot(enc%mesh%n, 0:enc%mmax, 0:enc%lmax))
        a = a4(:, :, :, 1); adot = a4(:, :, :, 2)
    contains
        subroutine field(rr, x, y, z, psi, active)
            real(dp), intent(in) :: rr, x(:), y(:), z(:)
            complex(dp), intent(out) :: psi(:, :)
            logical, intent(out) :: active
            complex(dp) :: v(size(x)), lin(size(x))
            real(dp) :: ux(size(x)), uz(size(x))
            integer :: j
            active = .false.
            psi = 0
            do j = 1, size(g)
                if (.not. used(j)) cycle
                if ((rr - sqrt(sum(g(j)%q**2)))**2*real(g(j)%alpha) > 45) cycle
                active = .true.
                v = gaussian_value(g(j), x, y, z)
                ux = x - g(j)%q(1); uz = z - g(j)%q(2)
                lin = -rate(j)%alpha*(ux**2 + y**2 + uz**2) &
                      + 2*g(j)%alpha*(ux*rate(j)%q(1) + uz*rate(j)%q(2)) &
                      + iu*(rate(j)%kappa(1)*ux + rate(j)%kappa(2)*uz) &
                      - iu*(g(j)%kappa(1)*rate(j)%q(1) + g(j)%kappa(2)*rate(j)%q(2))
                psi(:, 1) = psi(:, 1) + g(j)%d*v
                psi(:, 2) = psi(:, 2) + (rate(j)%d + g(j)%d*lin)*v
            end do
        end subroutine field
    end subroutine packet_project

    ! Errors, residuals and the Duhamel estimator of a model trajectory against
    ! the Q1 reference in the model's sector (sector-normalized reference).
    ! samples(1)%t must be 0.
    subroutine packet_trace(enc, spec, samples, tr)
        type(encounter_t), intent(in) :: enc
        type(packet_spec_t), intent(in) :: spec
        type(packet_sample_t), intent(in) :: samples(:)
        type(packet_trace_t), intent(out) :: tr
        complex(dp), allocatable :: am(:, :, :), ad(:, :, :), ar(:, :, :), r(:, :)
        real(dp), allocatable :: hb(:, :, :), x(:, :), yv(:, :)
        real(dp) :: ns, hbar, num, eh, ov2
        complex(dp) :: ov
        integer :: ns_t, it, l, nm, n, sector
        integer(i64) :: c0, c1, rate
        sector = model_sector(spec)
        hbar = enc%cfg%hbar; n = enc%mesh%n; ns_t = size(samples)
        allocate (tr%t(ns_t), tr%error(ns_t), tr%phase_error(ns_t), tr%residual(ns_t), &
                  tr%bound(ns_t), tr%total_error(ns_t), tr%total_bound(ns_t), &
                  tr%energy(ns_t), tr%norm(ns_t), tr%projected_norm(ns_t), &
                  tr%dvr_energy(ns_t))
        allocate (hb(enc%mesh%order, n, 0:enc%lmax))
        do l = 0, enc%lmax
            call radial_band_hamiltonian(enc%mesh, l, enc%cfg%g, enc%mu, hbar, hb(:, :, l))
        end do
        allocate (ar(n, 0:enc%mmax, 0:enc%lmax))
        call encounter_state(enc, 0.0_dp, ar)
        ns = encounter_norm(enc, ar, sector)
        tr%status = 'ok'
        call system_clock(c0, rate)
        do it = 1, ns_t
            tr%t(it) = samples(it)%t
            call packet_project(enc, samples(it)%g, samples(it)%rate, am, ad)
            call encounter_state(enc, samples(it)%t, ar)
            num = 0; eh = 0
            do l = 0, enc%lmax
                if (.not. in_sector(l, sector)) then
                    am(:, :, l) = 0; ad(:, :, l) = 0; ar(:, :, l) = 0
                    cycle
                end if
                nm = size(am, 2)
                allocate (x(n, 2*nm), yv(n, 2*nm))
                x(:, 1:nm) = real(am(:, :, l)); x(:, nm + 1:) = aimag(am(:, :, l))
                call radial_band_matvec(hb(:, :, l), x, yv)
                eh = eh + sum(x*yv)
                r = iu*hbar*ad(:, :, l) - cmplx(yv(:, 1:nm), yv(:, nm + 1:), dp)
                num = num + sum(abs(r)**2)
                deallocate (x, yv)
            end do
            ar = ar/sqrt(ns)
            tr%projected_norm(it) = sum(abs(am)**2)
            tr%dvr_energy(it) = eh/tr%projected_norm(it)
            tr%residual(it) = sqrt(num)
            tr%error(it) = sqrt(sum(abs(am - ar)**2))
            ov = sum(conjg(ar)*am)
            tr%phase_error(it) = sqrt(max(0.0_dp, tr%projected_norm(it) + sum(abs(ar)**2) &
                                          - 2*abs(ov)))
            ov2 = tr%projected_norm(it) + sum(abs(ar)**2) - 2*real(samples(it)%cm_overlap*ov)
            tr%total_error(it) = sqrt(max(0.0_dp, ov2))
            tr%energy(it) = samples(it)%energy; tr%norm(it) = samples(it)%norm
            if (it == 1) then
                tr%bound(1) = tr%error(1); tr%total_bound(1) = tr%total_error(1)
                tr%initial_error = tr%error(1)
            else
                tr%bound(it) = tr%bound(it - 1) + (tr%t(it) - tr%t(it - 1))* &
                               (tr%residual(it) + tr%residual(it - 1))/(2*hbar)
                tr%total_bound(it) = tr%total_bound(it - 1) + (tr%t(it) - tr%t(it - 1))* &
                                     (hypot(tr%residual(it), samples(it)%cm_residual) + &
                                      hypot(tr%residual(it - 1), samples(it - 1)%cm_residual)) &
                                     /(2*hbar)
            end if
            tr%max_condition = max(tr%max_condition, samples(it)%condition)
            tr%max_overlap_condition = max(tr%max_overlap_condition, &
                                           samples(it)%overlap_condition)
            tr%max_regularized = max(tr%max_regularized, samples(it)%regularized)
            if (samples(it)%overlap_condition > spec%refuse_condition) then
                tr%refused = .true.; tr%status = 'refused: overlap condition limit'
            end if
        end do
        call system_clock(c1)
        tr%seconds_projection = real(c1 - c0, dp)/real(rate, dp)
        call move_alloc(am, tr%a)
    end subroutine packet_trace

    ! Equation-level JVP of an M3 sector-normalized angular-bin probability with
    ! respect to (b, sigma) along dir, through the tangent of the canonical TDVP
    ! ODE (FortSym JVP leaf). The reference basis (mesh, l/m truncation) is
    ! frozen; finite differences must keep the same encounter_t.
    subroutine packet_bin_jvp(enc, spec, t, sector, rcut, th1, th2, dir, p, dpv, ok)
        type(encounter_t), intent(in) :: enc
        type(packet_spec_t), intent(in) :: spec
        real(dp), intent(in) :: t, rcut, th1, th2, dir(2)
        integer, intent(in) :: sector
        real(dp), intent(out) :: p, dpv
        logical, intent(out) :: ok
        type(canonical_t) :: can
        type(gaussian_t) :: g(1), dg(1)
        complex(dp), allocatable :: a(:, :, :), da(:, :, :)
        real(dp) :: y(13), dy(13), np, nm, pp, pm, pq(1), e, dgam
        integer(i64) :: nfev
        can%model = MODEL_M3; can%hbar = enc%cfg%hbar; can%mu = enc%mu
        can%g = enc%cfg%g; can%kap = spec%kappa
        y(1:7) = [enc%cfg%b, -enc%cfg%z0, 0.0_dp, can%hbar*enc%cfg%k, enc%cfg%sigma**2, &
                  0.0_dp, -enc%cfg%k*enc%cfg%z0]
        y(8:13) = [dir(1), 0.0_dp, 0.0_dp, 0.0_dp, 2*enc%cfg%sigma*dir(2), 0.0_dp]
        nfev = 0
        call integrate(canonical_rhs, can, y, 0.0_dp, t, spec%rtol, spec%atol, nfev, ok)
        p = 0; dpv = 0
        if (.not. ok) return
        call m3_rates(can, y(1:6), dy(1:6), e, dgam)
        call canonical_gaussian(y(1:6), y(7), y(8:13), 0.0_dp, can%hbar, g(1), dg(1))
        call packet_project(enc, g, dg, a, da)
        call encounter_angular_bins(enc, a, sector, rcut, [th1, th2], pq)
        p = pq(1)
        ! Polarization: Re<a|X|da> = (N(a+da) - N(a-da))/4 for quadratic N.
        call encounter_angular_bins(enc, a + da, sector, rcut, [th1, th2], pq)
        np = encounter_norm(enc, a + da, sector); pp = pq(1)*np
        call encounter_angular_bins(enc, a - da, sector, rcut, [th1, th2], pq)
        nm = encounter_norm(enc, a - da, sector); pm = pq(1)*nm
        dpv = ((pp - pm)/2 - p*(np - nm)/2)/encounter_norm(enc, a, sector)
    end subroutine packet_bin_jvp
end module kin6d_packet_models
