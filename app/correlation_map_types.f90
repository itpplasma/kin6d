!> Per-chain data of the H1a correlation map demo (application bookkeeping;
!> no physics). Raw accumulators are kept per batch so that every derived
!> quantity gets a batch/jackknife statistical interval.
module correlation_map_types
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_correlation_hierarchy, only: triplet_histogram_t
    implicit none
    private
    public :: chain_t, nouts, output_names, power_fit
    public :: O_H2, O_C3, O_U3, O_RED, O_DK, O_UG, O_UEX, O_UK, O_UL
    public :: O_EEX, O_EK, O_EL, O_XEX, O_XK, O_XL, O_M2MC, O_M2K, O_M2L
    public :: O_M0MC, O_EHNC, O_EDH, O_XHNC, O_XDH, O_RPEAK, O_GPEAK
    public :: O_C3OUT, O_C3COARSE, O_U3OUT, O_REDOUT

    integer, parameter :: O_H2 = 1, O_C3 = 2, O_U3 = 3, O_RED = 4, O_DK = 5
    integer, parameter :: O_UG = 6, O_UEX = 7, O_UK = 8, O_UL = 9
    integer, parameter :: O_EEX = 10, O_EK = 11, O_EL = 12
    integer, parameter :: O_XEX = 13, O_XK = 14, O_XL = 15
    integer, parameter :: O_M2MC = 16, O_M2K = 17, O_M2L = 18, O_M0MC = 19
    integer, parameter :: O_EHNC = 20, O_EDH = 21, O_XHNC = 22, O_XDH = 23
    integer, parameter :: O_RPEAK = 24, O_GPEAK = 25, O_C3OUT = 26, O_C3COARSE = 27
    integer, parameter :: O_U3OUT = 28, O_REDOUT = 29, nouts = 29
    character(16), parameter :: output_names(nouts) = [character(16) :: &
        'h2_sq', 'c3_sq', 'u3_sq', 'reducible_sq', 'kirkwood_rms', &
        'u_g2_route', 'u_ybg_exact', 'u_ybg_kirkwood', 'u_ybg_linear', &
        'l2_ybg_exact', 'l2_ybg_kirkwood', 'l2_ybg_linear', &
        'max_ybg_exact', 'max_ybg_kirkwood', 'max_ybg_linear', &
        'sl_mc', 'sl_kirkwood', 'sl_linear', 'neutral_mc', &
        'l2_hnc', 'l2_dh', 'max_hnc', 'max_dh', 'peak_r', 'peak_g', &
        'c3_sq_outer', 'c3_sq_coarse', 'u3_sq_outer', 'red_sq_outer']

    type :: chain_t
        real(dp) :: gamma = 0
        integer :: n = 0, nk = 0, bins = 0, shells = 16
        real(dp) :: length = 0, alpha = 0, bound = 0, step = 0, acceptance = 0
        real(dp) :: cpu = 0, rmax_g = 0
        !> beta U/N per production sweep and per equilibration sweep.
        real(dp), allocatable :: energy(:), equil(:)
        !> Per-batch accumulators (last index = batch).
        real(dp), allocatable :: pair(:, :), pair_force(:, :), force_sum(:, :)
        real(dp), allocatable :: triplet(:, :, :, :), sk(:, :)
        integer, allocatable :: snaps(:), fsnaps(:), sk_vectors(:)
        type(triplet_histogram_t) :: tri, coarse
        ! Analysis results.
        real(dp) :: u = 0, u_err = 0, u_batch_err = 0, tau = 0, ess = 0
        real(dp) :: var = 0, var_err = 0, dudg = 0, dudg_err = 0, u_tail = 0
        real(dp) :: u_hnc = 0, u_dh = 0, bin_indicator = 0, energy_drift = 0
        real(dp), allocatable :: r(:), g(:), g_err(:), g_exact(:), g_kirk(:)
        real(dp), allocatable :: g_lin(:), g_hnc(:), g_dh(:), g_dhl(:)
        real(dp), allocatable :: eq_r(:), eq_g3(:), eq_g3_err(:), eq_kirk(:), eq_lin(:)
        real(dp), allocatable :: s_k(:), s_mc(:), s_err(:), s_dh(:), s_hnc(:)
        real(dp) :: out(nouts) = 0, out_err(nouts) = 0
    end type
contains
    !> Weighted least-squares fit y = A x^p in log-log over points with
    !> y > 2 sigma (statistically resolved); sigma_p from the fit covariance.
    subroutine power_fit(x, y, sigma, p, amplitude, p_err, used)
        real(dp), intent(in) :: x(:), y(:), sigma(:)
        real(dp), intent(out) :: p, amplitude, p_err
        integer, intent(out) :: used
        real(dp) :: w(size(x)), lx(size(x)), ly(size(x)), sw, sx, sy, sxx, sxy, d
        logical :: ok(size(x))
        ok = y > 2*sigma .and. y > 0
        used = count(ok); p = 0; amplitude = 0; p_err = huge(1.0_dp)
        if (used < 2) return
        where (ok)
            lx = log(x); ly = log(y); w = (y/max(sigma, 1.0e-300_dp))**2
        elsewhere
            lx = 0; ly = 0; w = 0
        end where
        sw = sum(w); sx = sum(w*lx); sy = sum(w*ly); sxx = sum(w*lx*lx)
        sxy = sum(w*lx*ly); d = sw*sxx - sx*sx
        p = (sw*sxy - sx*sy)/d; amplitude = exp((sy - p*sx)/sw)
        p_err = sqrt(sw/d)
    end subroutine
end module correlation_map_types
