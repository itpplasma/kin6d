! Q1b: attractive hydrogen-like radial (Kepler) wave packet.
!
! A single partial wave l of the relative motion with g < 0 evolves on the same
! radial FEM-DVR and banded eigendecomposition as the scattering encounter.
! The preparation is a radial Gaussian u0(r) at rest at the outer classical
! turning point of the mean level nbar. Breathing at the Kepler period, collapse
! and fractional/full revivals follow from the represented spectrum; the
! independent oracle replaces every identified bound eigenvalue by the exact
! level E_n = -mu g^2/(2 hbar^2 n^2) in the autocorrelation sum.
module kin6d_coulomb_kepler
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use kin6d_radial_dvr, only: radial_mesh_t, build_radial_mesh, mesh_spacing, &
                                radial_eigen_window, radial_locate
    use kin6d_coulomb_encounter_leaves, only: ce_outer_turning, &
        ce_kepler_width, ce_radial_packet, ce_hydrogen_level, ce_kepler, ce_phase, &
        ce_state_bound
    use kin6d_species_leaves, only: species_pair_scales
    implicit none
    private
    public :: kepler_t, kepler_init, kepler_autocorrelation, kepler_mean_radius
    public :: kepler_radial_density, kepler_breathing

    type :: kepler_t
        type(radial_mesh_t) :: mesh
        real(dp) :: g = -1, mu = 1, hbar = 1, nbar = 20, a = 1, rc = 0, s = 0
        real(dp) :: period = 0, revival = 0
        integer :: l = 1, nbound = 0
        real(dp), allocatable :: e(:), v(:, :), resid(:), exact(:)
        complex(dp), allocatable :: c(:)
        real(dp), allocatable :: a0(:)
        ! Ledger: [bound] finite-space e0/rho, [indicator] continuum/window weight,
        ! spectrum agreement of identified bound states.
        real(dp) :: e0 = 0, rho = 0, window_discarded2 = 0, continuum_weight = 0
        real(dp) :: max_level_error = 0, weighted_level_error = 0, norm2 = 0
        real(dp) :: mean_n = 0, rms_n = 0, setup_seconds = 0
    end type kepler_t
contains
    ! Mean level nbar, partial wave l, rms level spread dn (sets the width).
    subroutine kepler_init(kep, g, mu, hbar, l, nbar, dn, order, ok)
        type(kepler_t), intent(out) :: kep
        real(dp), intent(in) :: g, mu, hbar, nbar, dn
        integer, intent(in) :: l, order
        logical, intent(out) :: ok
        real(dp) :: rmax, hmax, nhi, orth, en, w2, u
        integer :: i, j, c0, c1, rate
        logical :: done
        ok = .false.
        if (.not. (g < 0 .and. mu > 0 .and. hbar > 0 .and. dn > 0)) return
        if (.not. (nbar > real(l + 1, dp) + 3*dn)) return
        call system_clock(c0, rate)
        kep%g = g; kep%mu = mu; kep%hbar = hbar; kep%l = l; kep%nbar = nbar
        ! Relative-coordinate scales for reduced mass mu (unit masses split so
        ! that m1 m2/(m1+m2) = mu).
        call species_pair_scales(hbar, 2*mu, 2*mu, g, en, kep%a, w2)
        call ce_outer_turning(nbar, real(l, dp), kep%a, kep%rc)
        call ce_kepler_width(nbar, kep%a, dn, kep%s)
        call ce_kepler(nbar, mu, g, hbar, kep%period, kep%revival)
        nhi = nbar + 8*dn
        rmax = max(kep%rc + 14*kep%s, 3*nhi**2*kep%a)
        ! Local wavenumber of near-threshold states at r is sqrt(2/(a r));
        ! geometric growth from a/2 stays inside that resolution for r < rmax.
        hmax = mesh_spacing(sqrt(2/(kep%a*kep%rc/4)), order, 6.0_dp)
        call build_radial_mesh(kep%mesh, rmax, order, hmax, kep%a/2, 1.08_dp, done)
        if (.not. done) return
        allocate (kep%a0(kep%mesh%n))
        do i = 1, kep%mesh%n
            call ce_radial_packet(kep%mesh%r(i), kep%rc, kep%s, u)
            kep%a0(i) = u*sqrt(kep%mesh%w(i))
        end do
        kep%norm2 = sum(kep%a0**2)
        ! Window: all represented states below the ionization threshold plus a
        ! small continuum margin; the discarded weight is measured.
        call ce_hydrogen_level(nhi + 20, mu, g, hbar, en)
        call radial_eigen_window(kep%mesh, l, g, mu, hbar, abs(en), kep%e, kep%v, &
                                 kep%resid, orth, done)
        if (.not. done) return
        kep%c = matmul(transpose(kep%v), kep%a0)
        w2 = sum(abs(kep%c)**2)
        kep%window_discarded2 = max(0.0_dp, kep%norm2 - w2)
        kep%e0 = norm2(kep%a0 - real(matmul(kep%v, kep%c)))
        kep%rho = sum(abs(kep%c)*kep%resid)
        ! Bound states: E < 0 and spatially inside the box (n^2 a well below rmax).
        allocate (kep%exact(size(kep%e)))
        kep%exact = kep%e
        kep%nbound = 0
        do j = 1, size(kep%e)
            if (kep%e(j) >= 0) exit
            if (2.5_dp*real(l + j, dp)**2*kep%a > rmax) exit
            kep%nbound = j
            call ce_hydrogen_level(real(l + j, dp), mu, g, hbar, kep%exact(j))
        end do
        kep%continuum_weight = sum(abs(kep%c(kep%nbound + 1:))**2)
        kep%max_level_error = 0; kep%weighted_level_error = 0
        kep%mean_n = 0; kep%rms_n = 0
        do j = 1, kep%nbound
            if (abs(kep%c(j))**2 > 1.0e-12_dp) kep%max_level_error = max( &
                kep%max_level_error, abs(kep%e(j) - kep%exact(j))/abs(kep%exact(j)))
            kep%weighted_level_error = kep%weighted_level_error + &
                abs(kep%c(j))**2*abs(kep%e(j) - kep%exact(j))
            kep%mean_n = kep%mean_n + abs(kep%c(j))**2*real(l + j, dp)
            kep%rms_n = kep%rms_n + abs(kep%c(j))**2*real(l + j, dp)**2
        end do
        kep%mean_n = kep%mean_n/w2
        kep%rms_n = sqrt(max(0.0_dp, kep%rms_n/w2 - kep%mean_n**2))
        call system_clock(c1)
        kep%setup_seconds = real(c1 - c0, dp)/real(rate, dp)
        ok = all(ieee_is_finite(kep%e)) .and. kep%nbound > 0
    end subroutine kepler_init

    ! A(t) = <u0|u(t)> from the represented spectrum, and the oracle with every
    ! identified bound eigenvalue replaced by the exact hydrogen level.
    subroutine kepler_autocorrelation(kep, t, numeric, oracle)
        type(kepler_t), intent(in) :: kep
        real(dp), intent(in) :: t
        complex(dp), intent(out) :: numeric, oracle
        real(dp) :: cr, ci
        integer :: j
        numeric = 0; oracle = 0
        do j = 1, size(kep%e)
            call ce_phase(kep%e(j), t, kep%hbar, cr, ci)
            numeric = numeric + abs(kep%c(j))**2*cmplx(cr, ci, dp)
            call ce_phase(kep%exact(j), t, kep%hbar, cr, ci)
            oracle = oracle + abs(kep%c(j))**2*cmplx(cr, ci, dp)
        end do
    end subroutine kepler_autocorrelation

    subroutine kepler_state(kep, t, a)
        type(kepler_t), intent(in) :: kep
        real(dp), intent(in) :: t
        complex(dp), intent(out) :: a(:)
        complex(dp) :: ct(size(kep%e))
        real(dp) :: cr, ci
        integer :: j
        do j = 1, size(kep%e)
            call ce_phase(kep%e(j), t, kep%hbar, cr, ci)
            ct(j) = kep%c(j)*cmplx(cr, ci, dp)
        end do
        a = matmul(kep%v, ct)
    end subroutine kepler_state

    ! <r>(t) and the finite-space state bound e(t) <= e0 + t rho/hbar.
    subroutine kepler_mean_radius(kep, t, rmean, ebound)
        type(kepler_t), intent(in) :: kep
        real(dp), intent(in) :: t
        real(dp), intent(out) :: rmean, ebound
        complex(dp) :: a(kep%mesh%n)
        call kepler_state(kep, t, a)
        rmean = sum(kep%mesh%r*abs(a)**2)/sum(abs(a)**2)
        call ce_state_bound(kep%e0, kep%rho, t, kep%hbar, ebound)
    end subroutine kepler_mean_radius

    ! Breathing amplitude max <r> - min <r> over [t0, t1] (ns+1 samples).
    real(dp) function kepler_breathing(kep, t0, t1, ns) result(amp)
        type(kepler_t), intent(in) :: kep
        real(dp), intent(in) :: t0, t1
        integer, intent(in) :: ns
        real(dp) :: rm, eb, lo, hi
        integer :: i
        lo = huge(1.0_dp); hi = -huge(1.0_dp)
        do i = 0, ns
            call kepler_mean_radius(kep, t0 + (t1 - t0)*real(i, dp)/real(ns, dp), rm, eb)
            lo = min(lo, rm); hi = max(hi, rm)
        end do
        amp = hi - lo
    end function kepler_breathing

    ! Radial probability density |u(r)|^2 at the bin centres, from exact
    ! interpolation of the represented polynomial (no nodal-binning aliasing).
    subroutine kepler_radial_density(kep, t, edges, p)
        type(kepler_t), intent(in) :: kep
        real(dp), intent(in) :: t, edges(0:)
        real(dp), intent(out) :: p(:)
        complex(dp) :: a(kep%mesh%n)
        real(dp) :: val(kep%mesh%order)
        integer :: idx(kep%mesh%order), ib
        call kepler_state(kep, t, a)
        do ib = 1, size(p)
            call radial_locate(kep%mesh, (edges(ib - 1) + edges(ib))/2, idx, val)
            p(ib) = abs(sum(val*a(idx)))**2
        end do
    end subroutine kepler_radial_density
end module kin6d_coulomb_kepler
