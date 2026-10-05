! Q1/Q1b localized two-body quantum Coulomb encounter.
!
! Two nonrelativistic particles with bare Coulomb coupling g/r separate exactly
! into a free centre-of-mass Gaussian (analytic) and the 3D relative motion of
! reduced mass mu. The relative state is expanded in real spherical harmonics,
! psi(r) = sum_{l,m>=0} u_lm(r)/r P_lm(cos th) c_m cos(m phi), because the
! packet (centre (b,0,-z0), mean wave vector k z) is even in y. Every l evolves
! with H_l = -hbar^2/(2mu) d^2/dr^2 + V_l(r) on a radial FEM-DVR (see
! kin6d_radial_dvr); H_l is independent of m. Evolution is exact in the
! represented space: a(t) = V_l exp(-i E t/hbar) V_l^T a(0) with one banded
! symmetric eigendecomposition per l reused for every time. Exchange symmetry:
! Y_lm(-r) = (-1)^l Y_lm(r), so the spatial singlet keeps even l, the triplet
! odd l and the (distinguishable/product) preparation all l.
module kin6d_coulomb_encounter
    use, intrinsic :: iso_fortran_env, only: dp => real64, i64 => int64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use kin6d_radial_dvr, only: radial_mesh_t, build_radial_mesh, mesh_spacing, &
                                radial_band_hamiltonian, radial_band_matvec, &
                                radial_eigen_window, radial_locate
    use kin6d_coulomb_encounter_leaves
    use kin6d_species_leaves, only: species_pair_scales, species_sommerfeld
    use fortnum_quadrature, only: gauss_legendre, gauss_legendre_ab
    use fortnum_interval, only: interval_t, ipoint
    implicit none
    private
    public :: encounter_config_t, encounter_t, encounter_ledger_t
    public :: encounter_init, encounter_state, encounter_norm, encounter_energy
    public :: encounter_angular_bins, encounter_separation_density
    public :: encounter_relative_projection, encounter_lab_projection
    public :: encounter_wall_probability, encounter_state_error
    public :: encounter_bin_jvp, encounter_bin_vjp, encounter_bin_full
    public :: encounter_evaluate_packet, encounter_classical_ensemble
    public :: normalized_legendre, in_sector, encounter_points
    integer, parameter, public :: SECTOR_PRODUCT = 0, SECTOR_SINGLET = 1, &
                                  SECTOR_TRIPLET = 2

    type :: encounter_config_t
        real(dp) :: hbar = 1, g = 1, m1 = 1, m2 = 1      ! physical (declared units)
        real(dp) :: k = 1, sigma = 6, b = 0, z0 = 40      ! relative packet
        real(dp) :: sigma_cm = 0       ! CM rms width; 0 -> matched electrons sigma/2
        real(dp) :: t_final = 0        ! sets the automatic radial domain
        real(dp) :: rmax = 0           ! > 0 overrides the automatic domain
        real(dp) :: hmax = 0           ! > 0 overrides the automatic element size
        integer :: order = 12          ! GLL points per radial element
        real(dp) :: nodes_per_wavelength = 6
        real(dp) :: tolerance = 1.0e-10_dp   ! discarded norm^2 budget (l, m, spectrum)
        integer :: lmax = -1, mmax = -1      ! -1: adaptive
    end type encounter_config_t

    type :: partial_wave_t
        integer :: nm = 0
        real(dp), allocatable :: e(:), v(:, :), resid(:)
        complex(dp), allocatable :: c(:, :)    ! eigen coefficients (nkeep, nm)
        complex(dp), allocatable :: a0(:, :)   ! initial DVR coefficients (n, nm)
    end type partial_wave_t

    ! Error/approximation ledger. Labels: [bound] finite-space bound evaluated in
    ! binary64; [indicator] measured truncation/convergence indicator.
    type :: encounter_ledger_t
        real(dp) :: represented_norm2 = 0    ! sum |a_lm|^2 after l/m truncation
        real(dp) :: projection_defect = 0    ! [indicator] 1 - represented (analytic norm 1)
        real(dp) :: l_discarded2 = 0         ! [indicator] weight in l_est >= l > lmax
        real(dp) :: l_edge2 = 0              ! [indicator] weight in top 3 projected l
        real(dp) :: m_discarded2 = 0
        real(dp) :: m_edge2 = 0
        real(dp) :: spectral_discarded2 = 0  ! [bound] ||a0||^2 - ||V^T a0||^2
        real(dp) :: e0 = 0                   ! [bound] ||a0 - V V^T a0||
        real(dp) :: residual_rate = 0        ! [bound] rho: ||r(t)|| <= rho
        real(dp) :: orthogonality = 0        ! max |V^T V - I|
        real(dp) :: reconstruction_rms = 0   ! [indicator] pointwise off-node error
        real(dp) :: radial_tail2 = 0         ! [bound] analytic packet weight beyond rmax
        real(dp) :: setup_seconds = 0, projection_seconds = 0, eigen_seconds = 0
        integer :: l_estimate = 0, m_estimate = 0, kept_states = 0
        integer(i64) :: stored_bytes = 0
    end type encounter_ledger_t

    type :: encounter_t
        type(encounter_config_t) :: cfg
        type(radial_mesh_t) :: mesh
        real(dp) :: mu = 0, mtotal = 0, eta = 0, lambdabar = 0, r0 = 0
        real(dp) :: energy = 0, speed = 0, sigma_cm = 0, emin = 0, emax = 0, hmax = 0
        integer :: lmax = 0, mmax = 0
        type(partial_wave_t), allocatable :: wave(:)
        type(encounter_ledger_t) :: ledger
        logical :: ready = .false.
    end type encounter_t
contains
    pure logical function in_sector(l, sector)
        integer, intent(in) :: l, sector
        select case (sector)
        case (SECTOR_SINGLET)
            in_sector = mod(l, 2) == 0
        case (SECTOR_TRIPLET)
            in_sector = mod(l, 2) == 1
        case default
            in_sector = .true.
        end select
    end function in_sector

    ! Orthonormal associated Legendre functions: int_{-1}^{1} P_lm^2 dx = 1.
    pure subroutine normalized_legendre(lmax, mmax, x, p)
        integer, intent(in) :: lmax, mmax
        real(dp), intent(in) :: x
        real(dp), intent(out) :: p(0:lmax, 0:mmax)
        real(dp) :: s, a, aprev
        integer :: l, m
        p = 0
        s = sqrt(max(0.0_dp, 1 - x*x))
        p(0, 0) = sqrt(0.5_dp)
        do m = 0, min(mmax, lmax)
            if (m > 0) p(m, m) = sqrt(real(2*m + 1, dp)/real(2*m, dp))*s*p(max(0, m - 1), &
                                                                         max(0, m - 1))
            if (m + 1 <= lmax) p(m + 1, m) = sqrt(real(2*m + 3, dp))*x*p(m, m)
            aprev = sqrt(real(2*m + 3, dp))
            do l = m + 2, lmax
                a = sqrt(real(4*l*l - 1, dp)/real(l*l - m*m, dp))
                p(l, m) = a*(x*p(l - 1, m) - p(l - 2, m)/aprev)
                aprev = a
            end do
        end do
    end subroutine normalized_legendre

    pure real(dp) function cm_norm(m)
        integer, intent(in) :: m
        if (m == 0) then
            cm_norm = 1/sqrt(2*acos(-1.0_dp))
        else
            cm_norm = 1/sqrt(acos(-1.0_dp))
        end if
    end function cm_norm

    ! Analytic relative packet (or its (b,sigma) tangent) at Cartesian points.
    elemental subroutine encounter_evaluate_packet(cfg, x, y, z, db, dsigma, psi)
        type(encounter_config_t), intent(in) :: cfg
        real(dp), intent(in) :: x, y, z, db, dsigma
        complex(dp), intent(out) :: psi
        real(dp) :: re, im
        if (.not. (abs(db) > 0 .or. abs(dsigma) > 0)) then
            call ce_packet(x, y, z, cfg%b, cfg%z0, cfg%sigma, cfg%k, re, im)
        else
            call ce_packet_jvp(x, y, z, cfg%b, cfg%z0, cfg%sigma, cfg%k, 0.0_dp, 0.0_dp, &
                               0.0_dp, db, 0.0_dp, dsigma, 0.0_dp, re, im)
        end if
        psi = cmplx(re, im, dp)
    end subroutine encounter_evaluate_packet

    subroutine encounter_init(enc, cfg, ok, message)
        type(encounter_t), intent(out) :: enc
        type(encounter_config_t), intent(in) :: cfg
        logical, intent(out) :: ok
        character(:), allocatable, intent(out) :: message
        complex(dp), allocatable :: a0(:, :, :)
        real(dp), allocatable :: wl(:), wm(:)
        real(dp) :: sk, khi, vinf, rmax, hmax, h0, growth, st, se, tail, sig_e
        real(dp) :: bohr, rydberg
        integer(i64) :: c0, c1, c2, rate
        integer :: lest, mest, l, mcount
        logical :: built
        ok = .false.; message = 'invalid configuration'
        call system_clock(c0, rate)
        if (.not. valid_config(cfg)) return
        enc%cfg = cfg
        call species_pair_scales(cfg%hbar, cfg%m1, cfg%m2, cfg%g, enc%mu, bohr, rydberg)
        call ce_total_mass(cfg%m1, cfg%m2, enc%mtotal)
        call ce_scales(cfg%g, enc%mu, cfg%hbar, cfg%k, enc%lambdabar, enc%r0, &
                       enc%energy, enc%speed)
        call species_sommerfeld(cfg%g, cfg%hbar, enc%speed, enc%eta)
        sig_e = cfg%sigma/sqrt(2.0_dp)
        call ce_pair_widths(sig_e, se, enc%sigma_cm)
        if (cfg%sigma_cm > 0) enc%sigma_cm = cfg%sigma_cm
        call ce_momentum_width(cfg%sigma, cfg%hbar, sk)
        sk = sk/cfg%hbar
        vinf = sqrt(max(enc%speed**2, enc%speed**2 + 2*cfg%g/(enc%mu*cfg%z0)))
        khi = enc%mu*vinf/cfg%hbar + 8*sk
        call ce_free_width(cfg%sigma, enc%mu, cfg%hbar, cfg%t_final, st)
        rmax = cfg%rmax
        if (rmax <= 0) rmax = max(cfg%z0 + 8*cfg%sigma, vinf*cfg%t_final - cfg%z0 &
                                  + 8*st) + abs(cfg%b)
        hmax = mesh_spacing(khi, cfg%order, cfg%nodes_per_wavelength)
        if (cfg%hmax > 0) hmax = cfg%hmax
        h0 = hmax; growth = 1
        if (cfg%g < 0) then
            h0 = min(hmax, bohr); growth = 1.12_dp
        end if
        enc%hmax = hmax
        call build_radial_mesh(enc%mesh, rmax, cfg%order, hmax, h0, growth, built)
        message = 'radial mesh construction failed'
        if (.not. built) return
        ! Analytic weight of the initial packet outside the ball of radius rmax:
        ! the packet centre lies at distance d; Gaussian tail beyond rmax - d.
        tail = rmax - sqrt(cfg%b**2 + cfg%z0**2)
        enc%ledger%radial_tail2 = 1
        if (tail > 0) call ce_radial_tail(tail/cfg%sigma, enc%ledger%radial_tail2)
        ! Adaptive angular estimate, then trim by measured partial-wave weights.
        lest = ceil_nonneg(khi*(abs(cfg%b) + cfg%z0 + 8*cfg%sigma)) + 8
        if (cfg%lmax >= 0) lest = cfg%lmax
        mest = ceil_nonneg(khi*(abs(cfg%b) + 8*cfg%sigma)) + 8
        if (.not. (abs(cfg%b) > 0)) mest = 0
        if (cfg%mmax >= 0) mest = cfg%mmax
        mest = min(mest, lest)
        enc%ledger%l_estimate = lest; enc%ledger%m_estimate = mest
        call system_clock(c1)
        call project_packet(enc, lest, mest, 0.0_dp, 0.0_dp, a0)
        call system_clock(c2)
        enc%ledger%projection_seconds = real(c2 - c1, dp)/real(rate, dp)
        allocate (wl(0:lest), wm(0:mest))
        do l = 0, lest
            wl(l) = sum(abs(a0(:, :, l))**2)
        end do
        do l = 0, mest
            wm(l) = sum(abs(a0(:, l, :))**2)
        end do
        enc%ledger%l_edge2 = sum(wl(max(0, lest - 2):lest))
        enc%ledger%m_edge2 = sum(wm(max(0, mest - 2):mest))
        enc%lmax = trim_index(wl, cfg%tolerance/3, cfg%lmax >= 0)
        enc%mmax = trim_index(wm, cfg%tolerance/3, cfg%mmax >= 0)
        enc%ledger%l_discarded2 = sum(wl(enc%lmax + 1:lest))
        enc%ledger%m_discarded2 = sum(wm(enc%mmax + 1:mest))
        enc%ledger%represented_norm2 = 0
        do l = 0, enc%lmax
            enc%ledger%represented_norm2 = enc%ledger%represented_norm2 + &
                sum(abs(a0(:, 0:min(l, enc%mmax), l))**2)
        end do
        enc%ledger%projection_defect = 1 - enc%ledger%represented_norm2
        message = 'partial-wave eigen-window failed'
        allocate (enc%wave(0:enc%lmax))
        do l = 0, enc%lmax
            mcount = min(l, enc%mmax) + 1
            enc%wave(l)%nm = mcount
            enc%wave(l)%a0 = a0(:, 0:mcount - 1, l)
        end do
        deallocate (a0)
        call system_clock(c1)
        ! Packet energy window from its momentum band and the Coulomb energy over
        ! its extent; widened adaptively until the discarded weight fits the budget.
        enc%emax = cfg%hbar**2*khi**2/(2*enc%mu) + &
                   max(0.0_dp, cfg%g)/max(cfg%z0 - 6*cfg%sigma, cfg%z0/4)
        enc%emin = cfg%hbar**2*max(0.0_dp, cfg%k - 8*sk)**2/(2*enc%mu) + &
                   min(0.0_dp, cfg%g)/max(cfg%z0 - 6*cfg%sigma, cfg%z0/4)
        do l = 0, enc%lmax
            if (.not. solve_wave(enc, l)) return
        end do
        call system_clock(c2)
        enc%ledger%eigen_seconds = real(c2 - c1, dp)/real(rate, dp)
        call finalize_ledger(enc)
        enc%ledger%reconstruction_rms = reconstruction_error(enc)
        call system_clock(c2)
        enc%ledger%setup_seconds = real(c2 - c0, dp)/real(rate, dp)
        enc%ready = .true.; ok = .true.; message = ''
    end subroutine encounter_init

    pure logical function valid_config(cfg)
        type(encounter_config_t), intent(in) :: cfg
        valid_config = .false.
        if (.not. (ieee_is_finite(cfg%hbar) .and. ieee_is_finite(cfg%g))) return
        if (.not. (ieee_is_finite(cfg%k) .and. ieee_is_finite(cfg%sigma))) return
        if (.not. (ieee_is_finite(cfg%b) .and. ieee_is_finite(cfg%z0))) return
        if (.not. (ieee_is_finite(cfg%t_final) .and. ieee_is_finite(cfg%rmax))) return
        if (cfg%hbar <= 0 .or. cfg%k <= 0 .or. cfg%sigma <= 0) return
        if (cfg%m1 <= 0 .or. cfg%m2 <= 0 .or. cfg%z0 <= 0 .or. cfg%t_final < 0) return
        if (cfg%order < 4 .or. cfg%nodes_per_wavelength < 2) return
        if (.not. (cfg%tolerance > 0 .and. cfg%tolerance < 1)) return
        valid_config = .true.
    end function valid_config

    pure integer function ceil_nonneg(x)
        real(dp), intent(in) :: x
        ceil_nonneg = max(0, ceiling(x))
    end function ceil_nonneg

    ! Smallest L with sum_{l>L} w_l <= budget (fixed: keep all).
    pure integer function trim_index(w, budget, fixed)
        real(dp), intent(in) :: w(0:)
        real(dp), intent(in) :: budget
        logical, intent(in) :: fixed
        real(dp) :: tail
        integer :: j
        trim_index = ubound(w, 1)
        if (fixed) return
        tail = 0
        do j = ubound(w, 1), 1, -1
            if (tail + w(j) > budget) exit
            tail = tail + w(j)
            trim_index = j - 1
        end do
    end function trim_index

    ! a(i,m,l) = r_i sqrt(w_i) int Y_lm psi0 dOmega (Gauss-Legendre in cos th,
    ! periodic trapezoid in phi). With (db,dsigma) /= 0 the packet tangent is
    ! projected instead; the map is linear, so this is the exact tangent of a0.
    subroutine project_packet(enc, lest, mest, db, dsigma, a)
        type(encounter_t), intent(in) :: enc
        integer, intent(in) :: lest, mest
        real(dp), intent(in) :: db, dsigma
        complex(dp), allocatable, intent(out) :: a(:, :, :)
        real(dp), allocatable :: xt(:), wt(:), plm(:, :, :), cphi(:, :), phi(:)
        complex(dp), allocatable :: psi(:, :), fm(:)
        real(dp) :: d, amax, st, x, y, z, rr, scale
        integer :: nth, nph, i, j, q, m, n
        type(encounter_config_t) :: cfg
        cfg = enc%cfg
        n = enc%mesh%n
        nth = lest + 16
        d = sqrt(cfg%b**2 + cfg%z0**2)
        amax = abs(cfg%b)*(d + 12*cfg%sigma)/(2*cfg%sigma**2)
        nph = 2*(mest + ceiling(amax + 6*sqrt(amax)) + 8)
        allocate (xt(nth), wt(nth), plm(0:lest, 0:mest, nth), phi(nph), cphi(nph, 0:mest))
        call gauss_legendre(nth, xt, wt)
        do j = 1, nth
            call normalized_legendre(lest, mest, xt(j), plm(:, :, j))
        end do
        do q = 1, nph
            phi(q) = 2*acos(-1.0_dp)*real(q - 1, dp)/real(nph, dp)
            do m = 0, mest
                cphi(q, m) = cos(m*phi(q))*cm_norm(m)*2*acos(-1.0_dp)/real(nph, dp)
            end do
        end do
        allocate (a(n, 0:mest, 0:lest), psi(nph, nth), fm(nth))
        a = 0
        do i = 1, n
            rr = enc%mesh%r(i)
            if (abs(rr - d) > 13*cfg%sigma) cycle
            do j = 1, nth
                st = sqrt(max(0.0_dp, 1 - xt(j)**2))
                do q = 1, nph
                    call ce_cartesian(rr, xt(j), st, cos(phi(q)), sin(phi(q)), x, y, z)
                    call encounter_evaluate_packet(cfg, x, y, z, db, dsigma, psi(q, j))
                end do
            end do
            scale = rr*sqrt(enc%mesh%w(i))
            do m = 0, mest
                fm = matmul(cphi(:, m), psi)*wt*scale
                a(i, m, m:lest) = matmul(plm(m:lest, m, :), fm)
            end do
        end do
    end subroutine project_packet

    ! Per-l eigen window, adaptively widened until the discarded spectral weight
    ! of the initial state fits the budget; accumulates the ledger.
    logical function solve_wave(enc, l) result(ok)
        type(encounter_t), intent(inout) :: enc
        integer, intent(in) :: l
        real(dp), allocatable :: e(:), v(:, :), resid(:)
        complex(dp), allocatable :: c(:, :)
        real(dp) :: emin, emax, orth, w2, disc, budget, rho2, span
        integer :: attempt, m
        logical :: done
        ok = .false.
        associate (wv => enc%wave(l))
            w2 = sum(abs(wv%a0)**2)
            budget = enc%cfg%tolerance/(3*real(enc%lmax + 1, dp))
            emin = enc%emin; emax = enc%emax; span = enc%emax - enc%emin
            do attempt = 1, 8
                call radial_eigen_window(enc%mesh, l, enc%cfg%g, enc%mu, enc%cfg%hbar, &
                                         emax, e, v, resid, orth, done, emin)
                if (.not. done) return
                c = matmul(transpose(cmplx(v, 0.0_dp, dp)), wv%a0)
                disc = max(0.0_dp, w2 - sum(abs(c)**2))
                if (disc <= budget .or. size(e) == enc%mesh%n) exit
                emin = emin - span; emax = emax + span; span = 2*span
            end do
            if (disc > budget .and. size(e) < enc%mesh%n) return
            call move_alloc(e, wv%e); call move_alloc(v, wv%v)
            call move_alloc(resid, wv%resid); call move_alloc(c, wv%c)
            rho2 = 0
            do m = 1, wv%nm
                rho2 = rho2 + sum(abs(wv%c(:, m))*wv%resid)**2
            end do
            enc%ledger%spectral_discarded2 = enc%ledger%spectral_discarded2 + disc
            enc%ledger%e0 = enc%ledger%e0 + &
                sum(abs(wv%a0 - matmul(cmplx(wv%v, 0.0_dp, dp), wv%c))**2)
            enc%ledger%residual_rate = enc%ledger%residual_rate + rho2
            enc%ledger%orthogonality = max(enc%ledger%orthogonality, orth)
            enc%ledger%kept_states = enc%ledger%kept_states + size(wv%e)
            enc%ledger%stored_bytes = enc%ledger%stored_bytes + &
                int(8*size(wv%v) + 16*size(wv%c) + 16*size(wv%a0) + 16*size(wv%e), i64)
        end associate
        ok = .true.
    end function solve_wave

    subroutine finalize_ledger(enc)
        type(encounter_t), intent(inout) :: enc
        enc%ledger%e0 = sqrt(enc%ledger%e0)
        enc%ledger%residual_rate = sqrt(enc%ledger%residual_rate)
    end subroutine finalize_ledger

    ! Synthesize psi at Cartesian points from DVR coefficients a(n,0:mmax,0:lmax).
    subroutine encounter_points(enc, a, x, y, z, psi)
        type(encounter_t), intent(in) :: enc
        complex(dp), intent(in) :: a(:, 0:, 0:)
        real(dp), intent(in) :: x(:), y(:), z(:)
        complex(dp), intent(out) :: psi(:)
        real(dp) :: plm(0:enc%lmax, 0:enc%mmax), val(enc%mesh%order), rr, ct, phi
        integer :: idx(enc%mesh%order), p, l, m
        complex(dp) :: ul
        do p = 1, size(x)
            rr = sqrt(x(p)**2 + y(p)**2 + z(p)**2)
            psi(p) = 0
            if (rr <= 0 .or. rr >= enc%mesh%rmax) cycle
            ct = z(p)/rr; phi = atan2(y(p), x(p))
            call radial_locate(enc%mesh, rr, idx, val)
            call normalized_legendre(enc%lmax, enc%mmax, ct, plm)
            do l = 0, enc%lmax
                do m = 0, min(l, enc%mmax)
                    ul = sum(val*a(idx, m, l))
                    psi(p) = psi(p) + ul/rr*plm(l, m)*cm_norm(m)*cos(m*phi)
                end do
            end do
        end do
    end subroutine encounter_points

    ! RMS pointwise error of the represented initial packet relative to the RMS
    ! analytic amplitude, on a deterministic lattice within 2 sigma of the centre.
    real(dp) function reconstruction_error(enc) result(err)
        type(encounter_t), intent(in) :: enc
        integer, parameter :: np = 125
        real(dp) :: x(np), y(np), z(np), s
        complex(dp) :: psi(np), ref(np)
        complex(dp), allocatable :: a(:, :, :)
        integer :: p, l
        s = enc%cfg%sigma
        do p = 1, np
            x(p) = enc%cfg%b + s*(4*frac(0.7548776662_dp*p) - 2)
            y(p) = s*(4*frac(0.5698402910_dp*p) - 2)
            z(p) = -enc%cfg%z0 + s*(4*frac(0.4142135624_dp*p) - 2)
        end do
        allocate (a(enc%mesh%n, 0:enc%mmax, 0:enc%lmax))
        a = 0
        do l = 0, enc%lmax
            a(:, 0:enc%wave(l)%nm - 1, l) = enc%wave(l)%a0
        end do
        call encounter_points(enc, a, x, y, z, psi)
        call encounter_evaluate_packet(enc%cfg, x, y, z, 0.0_dp, 0.0_dp, ref)
        err = sqrt(sum(abs(psi - ref)**2)/sum(abs(ref)**2))
    contains
        pure real(dp) function frac(t)
            real(dp), intent(in) :: t
            frac = t - floor(t)
        end function frac
    end function reconstruction_error

    ! Represented relative state a(n,0:mmax,0:lmax) at time t (exact evolution
    ! of the retained eigen coefficients).
    subroutine encounter_state(enc, t, a)
        type(encounter_t), intent(in) :: enc
        real(dp), intent(in) :: t
        complex(dp), intent(out) :: a(:, 0:, 0:)
        complex(dp), allocatable :: ct(:, :)
        real(dp) :: cr, ci
        integer :: l, j
        a = 0
        do l = 0, enc%lmax
            associate (wv => enc%wave(l))
                ct = wv%c
                do j = 1, size(wv%e)
                    call ce_phase(wv%e(j), t, enc%cfg%hbar, cr, ci)
                    ct(j, :) = ct(j, :)*cmplx(cr, ci, dp)
                end do
                a(:, 0:wv%nm - 1, l) = matmul(wv%v, ct)
            end associate
        end do
    end subroutine encounter_state

    real(dp) function encounter_norm(enc, a, sector) result(s)
        type(encounter_t), intent(in) :: enc
        complex(dp), intent(in) :: a(:, 0:, 0:)
        integer, intent(in) :: sector
        integer :: l
        s = 0
        do l = 0, enc%lmax
            if (in_sector(l, sector)) s = s + sum(abs(a(:, :, l))**2)
        end do
    end function encounter_norm

    ! <H> per unit sector norm from an independent banded action (not eigenvalues).
    real(dp) function encounter_energy(enc, a, sector) result(en)
        type(encounter_t), intent(in) :: enc
        complex(dp), intent(in) :: a(:, 0:, 0:)
        integer, intent(in) :: sector
        real(dp), allocatable :: hb(:, :), x(:, :), y(:, :)
        integer :: l, nm
        en = 0
        allocate (hb(enc%mesh%order, enc%mesh%n))
        do l = 0, enc%lmax
            if (.not. in_sector(l, sector)) cycle
            nm = enc%wave(l)%nm
            call radial_band_hamiltonian(enc%mesh, l, enc%cfg%g, enc%mu, enc%cfg%hbar, hb)
            allocate (x(enc%mesh%n, 2*nm), y(enc%mesh%n, 2*nm))
            x(:, 1:nm) = real(a(:, 0:nm - 1, l)); x(:, nm + 1:) = aimag(a(:, 0:nm - 1, l))
            call radial_band_matvec(hb, x, y)
            en = en + sum(x*y)
            deallocate (x, y)
        end do
        en = en/encounter_norm(enc, a, sector)
    end function encounter_energy

    ! Nodal probabilities sum_{l in sector,m} |a_lm(i)|^2.
    subroutine node_weights(enc, a, sector, node)
        type(encounter_t), intent(in) :: enc
        complex(dp), intent(in) :: a(:, 0:, 0:)
        integer, intent(in) :: sector
        real(dp), intent(out) :: node(:)
        integer :: i, l
        node = 0
        do l = 0, enc%lmax
            if (.not. in_sector(l, sector)) cycle
            do i = 1, enc%mesh%n
                node(i) = node(i) + sum(abs(a(i, :, l))**2)
            end do
        end do
    end subroutine node_weights

    ! Probability of |r| in each radial bin, sector-normalized: 8-point
    ! Gauss-Legendre per bin on the exactly interpolated sum_lm |u_lm(r)|^2.
    subroutine encounter_separation_density(enc, a, sector, edges, p)
        type(encounter_t), intent(in) :: enc
        complex(dp), intent(in) :: a(:, 0:, 0:)
        integer, intent(in) :: sector
        real(dp), intent(in) :: edges(0:)
        real(dp), intent(out) :: p(:)
        real(dp) :: xq(8), wq(8), val(enc%mesh%order), dens
        integer :: idx(enc%mesh%order), ib, q, l, m
        p = 0
        do ib = 1, size(p)
            call gauss_legendre_ab(8, edges(ib - 1), min(edges(ib), enc%mesh%rmax), xq, wq)
            do q = 1, 8
                call radial_locate(enc%mesh, xq(q), idx, val)
                dens = 0
                do l = 0, enc%lmax
                    if (.not. in_sector(l, sector)) cycle
                    do m = 0, min(l, enc%mmax)
                        dens = dens + abs(sum(val*a(idx, m, l)))**2
                    end do
                end do
                p(ib) = p(ib) + wq(q)*dens
            end do
        end do
        p = p/encounter_norm(enc, a, sector)
    end subroutine encounter_separation_density

    ! [indicator] probability within width of the hard wall at rmax.
    real(dp) function encounter_wall_probability(enc, a, sector, width) result(p)
        type(encounter_t), intent(in) :: enc
        complex(dp), intent(in) :: a(:, 0:, 0:)
        integer, intent(in) :: sector
        real(dp), intent(in) :: width
        real(dp) :: node(enc%mesh%n)
        call node_weights(enc, a, sector, node)
        p = sum(node, mask=enc%mesh%r > enc%mesh%rmax - width)/sum(node)
    end function encounter_wall_probability

    ! M^m_{ll'} = int_{cos th2}^{cos th1} P_lm P_l'm dx (exact Gauss-Legendre).
    subroutine bin_matrices(enc, th1, th2, mats)
        type(encounter_t), intent(in) :: enc
        real(dp), intent(in) :: th1, th2
        real(dp), intent(out) :: mats(0:enc%lmax, 0:enc%lmax, 0:enc%mmax)
        real(dp) :: xq(enc%lmax + 2), wq(enc%lmax + 2)
        real(dp) :: p(0:enc%lmax, 0:enc%mmax, enc%lmax + 2), pw(0:enc%lmax, enc%lmax + 2)
        integer :: q, m
        call gauss_legendre_ab(enc%lmax + 2, cos(th2), cos(th1), xq, wq)
        do q = 1, size(xq)
            call normalized_legendre(enc%lmax, enc%mmax, xq(q), p(:, :, q))
        end do
        do m = 0, enc%mmax
            do q = 1, size(xq)
                pw(:, q) = p(:, m, q)*wq(q)
            end do
            mats(:, :, m) = matmul(pw, transpose(p(:, m, :)))
        end do
    end subroutine bin_matrices

    ! Sector-normalized probability of the direction theta of r in each bin
    ! (phi integrated) restricted to |r| >= rcut. A finite-time observable.
    subroutine encounter_angular_bins(enc, a, sector, rcut, edges, p)
        type(encounter_t), intent(in) :: enc
        complex(dp), intent(in) :: a(:, 0:, 0:)
        integer, intent(in) :: sector
        real(dp), intent(in) :: rcut, edges(0:)
        real(dp), intent(out) :: p(:)
        real(dp), allocatable :: rr(:, :, :), mats(:, :, :)
        real(dp) :: nrm
        integer :: ib
        call cut_overlaps(enc, a, sector, rcut, rr)
        allocate (mats(0:enc%lmax, 0:enc%lmax, 0:enc%mmax))
        nrm = encounter_norm(enc, a, sector)
        do ib = 1, size(p)
            call bin_matrices(enc, edges(ib - 1), edges(ib), mats)
            p(ib) = sum(rr*mats)/nrm
        end do
    end subroutine encounter_angular_bins

    ! Re R^m_{ll'} = Re sum_{r_i>=rcut} conj(a_lm(i)) a_l'm(i), zero outside sector.
    subroutine cut_overlaps(enc, a, sector, rcut, rr)
        type(encounter_t), intent(in) :: enc
        complex(dp), intent(in) :: a(:, 0:, 0:)
        integer, intent(in) :: sector
        real(dp), intent(in) :: rcut
        real(dp), allocatable, intent(out) :: rr(:, :, :)
        complex(dp), allocatable :: blk(:, :)
        integer :: i0, m, l
        allocate (rr(0:enc%lmax, 0:enc%lmax, 0:enc%mmax))
        i0 = first_node(enc, rcut)
        do m = 0, enc%mmax
            blk = a(i0:, m, :)
            do l = 0, enc%lmax
                if (.not. in_sector(l, sector)) blk(:, l + 1) = 0
            end do
            rr(:, :, m) = real(matmul(conjg(transpose(blk)), blk))
        end do
    end subroutine cut_overlaps

    pure integer function first_node(enc, rcut) result(i0)
        type(encounter_t), intent(in) :: enc
        real(dp), intent(in) :: rcut
        i0 = 1
        do while (i0 < enc%mesh%n .and. enc%mesh%r(i0) < rcut)
            i0 = i0 + 1
        end do
    end function first_node

    ! Finite-space state error bound e(t) <= e0 + t rho/hbar and the radius of a
    ! sector-normalized bounded observable in [0,1] (outward interval leaves;
    ! inputs e0, rho evaluated in binary64).
    subroutine encounter_state_error(enc, t, nrm, p, e, radius)
        type(encounter_t), intent(in) :: enc
        real(dp), intent(in) :: t, nrm, p
        real(dp), intent(out) :: e, radius
        type(interval_t) :: ei, ri
        call ce_state_bound_interval(ipoint(enc%ledger%e0), &
                                     ipoint(enc%ledger%residual_rate), ipoint(t), &
                                     ipoint(enc%cfg%hbar), ei)
        call ce_sector_radius_interval(ipoint(nrm), ipoint(p), ipoint(ei%hi), ri)
        e = ei%hi; radius = ri%hi
    end subroutine encounter_state_error

    ! Relative-coordinate probability per unit area in the (x,z) plane, integrated
    ! over y (trapezoid, ny points on [-rmax,rmax]); image(nz,nx) at pixel centres.
    ! mirror=.true. returns (|psi(r)|^2+|psi(-r)|^2)/2 (spin-unresolved product).
    ! Visualization path: psi is bilinearly interpolated from a uniform (rho,theta)
    ! table of spacing dtab (relative amplitude error ~ (k dtab)^2/8).
    subroutine encounter_relative_projection(enc, a, sector, mirror, xe, ze, ny, dtab, image)
        type(encounter_t), intent(in) :: enc
        complex(dp), intent(in) :: a(:, 0:, 0:)
        integer, intent(in) :: sector, ny
        logical, intent(in) :: mirror
        real(dp), intent(in) :: xe(0:), ze(0:), dtab
        real(dp), intent(out) :: image(:, :)
        complex(dp), allocatable :: u(:, :), tab(:, :), psi(:, :, :, :)
        real(dp), allocatable :: plm(:, :, :), ptab(:, :)
        real(dp) :: val(enc%mesh%order), rmx, dth, x, y, z, rho, th, fr, ft, dy, wy, nrm
        integer :: idx(enc%mesh%order), nr, nt, ia, ib, ix, iz, iy, m, l, s, ns
        nrm = encounter_norm(enc, a, sector)
        rmx = enc%mesh%rmax
        nr = ceiling(rmx/dtab) + 1
        nt = ceiling(acos(-1.0_dp)*rmx/dtab) + 2
        dth = acos(-1.0_dp)/real(nt - 1, dp)
        allocate (plm(0:enc%lmax, 0:enc%mmax, nt), u(nr, 0:enc%lmax), ptab(0:enc%lmax, nt))
        do ib = 1, nt
            call normalized_legendre(enc%lmax, enc%mmax, cos(real(ib - 1, dp)*dth), &
                                     plm(:, :, ib))
        end do
        ns = 1
        if (mirror) ns = 2
        allocate (psi(ny, size(ze) - 1, size(xe) - 1, ns))
        psi = 0
        dy = 2*rmx/real(ny - 1, dp)
        do m = 0, enc%mmax
            u = 0
            do ia = 2, nr
                rho = min(rmx, real(ia - 1, dp)*dtab)
                call radial_locate(enc%mesh, rho, idx, val)
                do l = m, enc%lmax
                    if (in_sector(l, sector)) u(ia, l) = sum(val*a(idx, m, l))/rho
                end do
            end do
            ptab = plm(:, m, :)*cm_norm(m)
            tab = matmul(u, ptab)
            do s = 1, ns
                do ix = 1, size(xe) - 1
                    do iz = 1, size(ze) - 1
                        do iy = 1, ny
                            x = (xe(ix - 1) + xe(ix))/2; z = (ze(iz - 1) + ze(iz))/2
                            y = -rmx + real(iy - 1, dp)*dy
                            if (s == 2) then
                                x = -x; y = -y; z = -z
                            end if
                            rho = sqrt(x*x + y*y + z*z)
                            if (rho >= rmx .or. rho <= 0) cycle
                            th = acos(max(-1.0_dp, min(1.0_dp, z/rho)))
                            fr = rho/dtab; ia = min(nr - 1, int(fr) + 1); fr = fr - (ia - 1)
                            ft = th/dth; ib = min(nt - 1, int(ft) + 1); ft = ft - (ib - 1)
                            psi(iy, iz, ix, s) = psi(iy, iz, ix, s) + cos(m*atan2(y, x))* &
                                ((1 - fr)*((1 - ft)*tab(ia, ib) + ft*tab(ia, ib + 1)) + &
                                 fr*((1 - ft)*tab(ia + 1, ib) + ft*tab(ia + 1, ib + 1)))
                        end do
                    end do
                end do
            end do
        end do
        do ix = 1, size(xe) - 1
            do iz = 1, size(ze) - 1
                image(iz, ix) = 0
                do iy = 1, ny
                    wy = dy
                    if (iy == 1 .or. iy == ny) wy = dy/2
                    image(iz, ix) = image(iz, ix) + wy*sum(abs(psi(iy, iz, ix, :))**2)/ns
                end do
            end do
        end do
        image = image/nrm
    end subroutine encounter_relative_projection

    ! One-body (x,z) probability density of particle `which` (1 or 2) in the
    ! CM rest frame: x1 = X + (m2/M) r, x2 = X - (m1/M) r, convolving the relative
    ! pixel probabilities with the analytic free CM Gaussian at time t.
    subroutine encounter_lab_projection(enc, rel, xr, zr, t, which, xl, zl, image)
        type(encounter_t), intent(in) :: enc
        real(dp), intent(in) :: rel(:, :), xr(0:), zr(0:), t, xl(0:), zl(0:)
        integer, intent(in) :: which
        real(dp), intent(out) :: image(:, :)
        real(dp), allocatable :: ax(:, :), az(:, :), prob(:, :)
        real(dp) :: width, s, lab1, lab2
        integer :: i, j
        call ce_free_width(enc%sigma_cm, enc%mtotal, enc%cfg%hbar, t, width)
        call ce_lab_positions(0.0_dp, 1.0_dp, enc%cfg%m1, enc%cfg%m2, lab1, lab2)
        s = lab1
        if (which == 2) s = lab2
        allocate (ax(size(xl) - 1, size(xr) - 1), az(size(zl) - 1, size(zr) - 1))
        do j = 1, size(xr) - 1
            do i = 1, size(xl) - 1
                call ce_gaussian_density((xl(i - 1) + xl(i))/2, s*(xr(j - 1) + xr(j))/2, &
                                         width, ax(i, j))
            end do
        end do
        do j = 1, size(zr) - 1
            do i = 1, size(zl) - 1
                call ce_gaussian_density((zl(i - 1) + zl(i))/2, s*(zr(j - 1) + zr(j))/2, &
                                         width, az(i, j))
            end do
        end do
        allocate (prob(size(zr) - 1, size(xr) - 1))
        do j = 1, size(xr) - 1
            prob(:, j) = rel(:, j)*(xr(j) - xr(j - 1))*(zr(1:) - zr(:size(zr) - 2))
        end do
        image = matmul(matmul(az, prob), transpose(ax))
    end subroutine encounter_lab_projection

    ! Classical comparison: independent point-pair realizations sampled from the
    ! positive initial Wigner function of the Gaussian (relative and CM parts),
    ! relative motion integrated by the interaction-laboratory Verlet step with
    ! encounter-resolved step control; CM motion is free. Repulsive g > 0 only.
    subroutine encounter_classical_ensemble(enc, nsample, seed, times, rel, cm, drift, ok)
        use fortnum_rng, only: rng_t, rng_seed, rng_normal
        use fortnum_status, only: fortnum_status_t
        use kin6d_interaction_laboratory, only: interaction_state_t, init_interaction, &
            step_interaction, sample_interaction
        use kin6d_interaction_leaves, only: interaction_energy
        use kin6d_species, only: units_t, species_t
        type(encounter_t), intent(in) :: enc
        integer, intent(in) :: nsample
        integer(i64), intent(in) :: seed
        real(dp), intent(in) :: times(:)
        real(dp), intent(out) :: rel(:, :, :), cm(:, :, :), drift
        logical, intent(out) :: ok
        type(rng_t) :: gen
        type(fortnum_status_t) :: status
        type(interaction_state_t) :: st
        type(units_t) :: units
        type(species_t) :: first, second
        type(interval_t) :: rbi(3), ubi(3)
        real(dp) :: r(3), u(3), x(3), p(3), alpha, sk, skc, t, h, e0, e1
        integer :: ns, it, k, d
        logical :: stepped, bounded
        ok = .false.; drift = 0
        if (enc%cfg%g <= 0 .or. nsample < 1) return
        alpha = enc%cfg%g/enc%mu
        ! Declared units with Coulomb constant 1: unit charge on particle 1 and
        ! charge g on particle 2 reproduce the pair coupling g exactly.
        units%hbar = enc%cfg%hbar; units%coulomb = 1
        first = species_t('particle1', 1.0_dp, enc%cfg%m1)
        second = species_t('particle2', enc%cfg%g, enc%cfg%m2)
        call rng_seed(gen, seed, status)
        call ce_momentum_width(enc%cfg%sigma, enc%cfg%hbar, sk)
        call ce_momentum_width(enc%sigma_cm, enc%cfg%hbar, skc)
        do ns = 1, nsample
            do d = 1, 3
                call rng_normal(gen, r(d)); call rng_normal(gen, u(d))
                call rng_normal(gen, x(d)); call rng_normal(gen, p(d))
            end do
            r = r*enc%cfg%sigma + [enc%cfg%b, 0.0_dp, -enc%cfg%z0]
            u = (u*sk + [0.0_dp, 0.0_dp, enc%cfg%hbar*enc%cfg%k])/enc%mu
            x = x*enc%sigma_cm; p = p*skc
            call init_interaction(st, r, u, units, first, second, stepped)
            if (.not. stepped) return
            call interaction_energy(r(1), r(2), r(3), u(1), u(2), u(3), alpha, e0)
            t = 0
            do it = 1, size(times)
                do k = 1, 100000000
                    if (t >= times(it)) exit
                    call sample_interaction(st, r, u, rbi, ubi, stepped, bounded)
                    h = 0.005_dp*min(norm2(r)/max(norm2(u), tiny(1.0_dp)), &
                                     sqrt(norm2(r)**3/alpha))
                    h = min(h, times(it) - t)
                    call step_interaction(st, h, stepped, certify=.false.)
                    if (.not. stepped) return
                    t = t + h
                end do
                call sample_interaction(st, r, u, rbi, ubi, stepped, bounded)
                rel(:, ns, it) = r
                cm(:, ns, it) = x + p*times(it)/enc%mtotal
            end do
            call interaction_energy(r(1), r(2), r(3), u(1), u(2), u(3), alpha, e1)
            drift = max(drift, abs(e1 - e0)/abs(e0))
        end do
        ok = .true.
    end subroutine encounter_classical_ensemble

    ! ---- Equation-level derivatives of a sector-normalized angular-bin
    ! probability P(t) = <a|Pi Q|a>/<a|Q|a>, a(t) = exp(-i H(g) t/hbar) a0(b,sigma),
    ! in the full represented space (complete eigenbasis per l; no window).
    ! Tangent: da = U da0 + dg * V[(V^T (dH/dg) V o F) V^T a0], F the Duhamel
    ! (Daleckii-Krein) kernel. The adjoint uses the transposed contractions.

    subroutine full_eigen(enc, l, e, v)
        type(encounter_t), intent(in) :: enc
        integer, intent(in) :: l
        real(dp), allocatable, intent(out) :: e(:), v(:, :)
        real(dp), allocatable :: resid(:)
        real(dp) :: orth
        logical :: ok
        call radial_eigen_window(enc%mesh, l, enc%cfg%g, enc%mu, enc%cfg%hbar, &
                                 huge(1.0_dp), e, v, resid, orth, ok)
        if (.not. ok) error stop 'coulomb_encounter: full eigendecomposition failed'
    end subroutine full_eigen

    ! dH/dg = diag(dV_l/dg (r_i)) in the eigenbasis, Hadamard with the kernel F.
    subroutine coupling_kernel(enc, l, t, e, v, kern)
        type(encounter_t), intent(in) :: enc
        integer, intent(in) :: l
        real(dp), intent(in) :: t, e(:), v(:, :)
        complex(dp), allocatable, intent(out) :: kern(:, :)
        real(dp), allocatable :: dv(:), wt(:, :)
        real(dp) :: x, sn, re, im
        integer :: i, j, k
        allocate (dv(enc%mesh%n))
        do i = 1, enc%mesh%n
            call ce_radial_potential_jvp(enc%mesh%r(i), real(l, dp), enc%cfg%g, enc%mu, &
                                         enc%cfg%hbar, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, &
                                         0.0_dp, dv(i))
        end do
        wt = matmul(transpose(v), spread(dv, 2, size(e))*v)
        allocate (kern(size(e), size(e)))
        do k = 1, size(e)
            do j = 1, size(e)
                call ce_sinc_argument(e(j), e(k), t, enc%cfg%hbar, x)
                if (abs(x) < 1.0e-3_dp) then
                    call ce_sinc_series(x, sn)
                else
                    call ce_sinc(x, sn)
                end if
                call ce_duhamel(e(j), e(k), t, enc%cfg%hbar, sn, re, im)
                kern(j, k) = wt(j, k)*cmplx(re, im, dp)
            end do
        end do
    end subroutine coupling_kernel

    subroutine phases(enc, e, t, ph)
        type(encounter_t), intent(in) :: enc
        real(dp), intent(in) :: e(:), t
        complex(dp), intent(out) :: ph(:)
        real(dp) :: cr, ci
        integer :: j
        do j = 1, size(e)
            call ce_phase(e(j), t, enc%cfg%hbar, cr, ci)
            ph(j) = cmplx(cr, ci, dp)
        end do
    end subroutine phases

    ! Full-space primal state, probability and covector lam = (Pi - P) Q a / den.
    subroutine bin_primal(enc, t, sector, rcut, th1, th2, a, lam, p)
        type(encounter_t), intent(in) :: enc
        real(dp), intent(in) :: t, rcut, th1, th2
        integer, intent(in) :: sector
        complex(dp), allocatable, intent(out) :: a(:, :, :), lam(:, :, :)
        real(dp), intent(out) :: p
        real(dp), allocatable :: e(:), v(:, :), rr(:, :, :), mats(:, :, :)
        complex(dp), allocatable :: ph(:)
        real(dp) :: den
        integer :: l, m, nm, i0
        allocate (a(enc%mesh%n, 0:enc%mmax, 0:enc%lmax), lam(enc%mesh%n, 0:enc%mmax, &
                                                                0:enc%lmax))
        a = 0; lam = 0
        do l = 0, enc%lmax
            if (.not. in_sector(l, sector)) cycle
            nm = enc%wave(l)%nm
            call full_eigen(enc, l, e, v)
            allocate (ph(size(e)))
            call phases(enc, e, t, ph)
            a(:, 0:nm - 1, l) = matmul(v, spread(ph, 2, nm)*matmul(transpose(v), &
                                                                  enc%wave(l)%a0))
            deallocate (ph)
        end do
        call cut_overlaps(enc, a, sector, rcut, rr)
        allocate (mats(0:enc%lmax, 0:enc%lmax, 0:enc%mmax))
        call bin_matrices(enc, th1, th2, mats)
        den = encounter_norm(enc, a, sector)
        p = sum(rr*mats)/den
        i0 = first_node(enc, rcut)
        do m = 0, enc%mmax
            do l = 0, enc%lmax
                if (.not. in_sector(l, sector)) mats(l, :, m) = 0
                if (.not. in_sector(l, sector)) mats(:, l, m) = 0
            end do
            lam(i0:, m, :) = matmul(a(i0:, m, :), mats(:, :, m))
        end do
        lam = (lam - p*a)/den
    end subroutine bin_primal

    subroutine encounter_bin_full(enc, t, sector, rcut, th1, th2, p)
        type(encounter_t), intent(in) :: enc
        real(dp), intent(in) :: t, rcut, th1, th2
        integer, intent(in) :: sector
        real(dp), intent(out) :: p
        complex(dp), allocatable :: a(:, :, :), lam(:, :, :)
        call bin_primal(enc, t, sector, rcut, th1, th2, a, lam, p)
    end subroutine encounter_bin_full

    ! Forward tangent along direction (db, dsigma, dg).
    subroutine encounter_bin_jvp(enc, t, sector, rcut, th1, th2, dir, p, dp_out)
        type(encounter_t), intent(in) :: enc
        real(dp), intent(in) :: t, rcut, th1, th2, dir(3)
        integer, intent(in) :: sector
        real(dp), intent(out) :: p, dp_out
        complex(dp), allocatable :: a(:, :, :), lam(:, :, :), da0(:, :, :), kern(:, :)
        complex(dp), allocatable :: c(:, :), dc(:, :), ph(:)
        real(dp), allocatable :: e(:), v(:, :)
        integer :: l, nm
        call bin_primal(enc, t, sector, rcut, th1, th2, a, lam, p)
        call project_packet(enc, enc%lmax, enc%mmax, dir(1), dir(2), da0)
        dp_out = 0
        do l = 0, enc%lmax
            if (.not. in_sector(l, sector)) cycle
            nm = enc%wave(l)%nm
            call full_eigen(enc, l, e, v)
            allocate (ph(size(e)))
            call phases(enc, e, t, ph)
            c = matmul(transpose(v), enc%wave(l)%a0)
            dc = spread(ph, 2, nm)*matmul(transpose(v), da0(:, 0:nm - 1, l))
            if (abs(dir(3)) > 0) then
                call coupling_kernel(enc, l, t, e, v, kern)
                dc = dc + dir(3)*matmul(kern, c)
            end if
            dp_out = dp_out + 2*real(sum(conjg(lam(:, 0:nm - 1, l))*matmul(v, dc)))
            deallocate (ph)
        end do
    end subroutine encounter_bin_jvp

    ! Reverse gradient (dP/db, dP/dsigma, dP/dg) from one backward propagation.
    subroutine encounter_bin_vjp(enc, t, sector, rcut, th1, th2, p, grad)
        type(encounter_t), intent(in) :: enc
        real(dp), intent(in) :: t, rcut, th1, th2
        integer, intent(in) :: sector
        real(dp), intent(out) :: p, grad(3)
        complex(dp), allocatable :: a(:, :, :), lam(:, :, :), mu0(:, :, :), d0(:, :, :)
        complex(dp), allocatable :: kern(:, :), lt(:, :), c(:, :), ph(:)
        real(dp), allocatable :: e(:), v(:, :)
        integer :: l, nm
        call bin_primal(enc, t, sector, rcut, th1, th2, a, lam, p)
        allocate (mu0(enc%mesh%n, 0:enc%mmax, 0:enc%lmax))
        mu0 = 0; grad = 0
        do l = 0, enc%lmax
            if (.not. in_sector(l, sector)) cycle
            nm = enc%wave(l)%nm
            call full_eigen(enc, l, e, v)
            allocate (ph(size(e)))
            call phases(enc, e, t, ph)
            lt = matmul(transpose(v), lam(:, 0:nm - 1, l))
            mu0(:, 0:nm - 1, l) = matmul(v, spread(conjg(ph), 2, nm)*lt)
            call coupling_kernel(enc, l, t, e, v, kern)
            c = matmul(transpose(v), enc%wave(l)%a0)
            grad(3) = grad(3) + 2*real(sum(conjg(matmul(conjg(transpose(kern)), lt))*c))
            deallocate (ph)
        end do
        call project_packet(enc, enc%lmax, enc%mmax, 1.0_dp, 0.0_dp, d0)
        grad(1) = 2*real(sum(conjg(mu0)*d0))
        call project_packet(enc, enc%lmax, enc%mmax, 0.0_dp, 1.0_dp, d0)
        grad(2) = 2*real(sum(conjg(mu0)*d0))
    end subroutine encounter_bin_vjp
end module kin6d_coulomb_encounter
