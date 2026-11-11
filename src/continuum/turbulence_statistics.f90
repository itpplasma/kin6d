!> Fluctuation statistics of a statistically stationary, homogeneous run on a
!> periodic box: the output contract of docs/drift-turbulence.md.
!>
!> Input: equally spaced samples of the Fourier coefficients phi(kx, ky) and
!> n(kx, ky), ky >= 0, of two real fields, normalised so that the box average
!> of phi^2 is the sum of |phi(k)|^2 over all k. Time dependence is counted as
!> exp(-i omega t).
!>
!> Spectra are Welch estimates: Hann-windowed segments with half overlap,
!> normalised by the window power, so that the sum of S_phi_phi over all
!> cells times dkx dky domega is the window-weighted variance of phi.
!> Correlation functions are the transforms of the spectrum divided by the
!> autocorrelation of the window. Every interval is statistical: one standard
!> error from batches of contiguous segments, or from the integrated
!> autocorrelation time of a series; it is never a deterministic enclosure.
!> A negative error means "not available" (fewer than two complete batches),
!> a negative correlation scale "1/e not reached inside the lag range".
module milca_turbulence_statistics
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_fft, only: fortnum_fft_plan_t, fft_c2c_plan_init, fft_c2c
    use milca_series_statistics, only: series_autocorrelation
    implicit none
    private
    public :: turbulence_scalar_t, turbulence_spectrum_t, turbulence_accumulator_t
    public :: turbulence_series_scalar, turbulence_accumulator_init, &
        turbulence_accumulate, turbulence_spectrum

    real(dp), parameter :: pi = 3.141592653589793238462643383279502884_dp
    integer, parameter :: batch_slots = 16

    !> Estimate with one standard error; tau is the integrated autocorrelation
    !> time of the series it came from (zero for batch estimates).
    type :: turbulence_scalar_t
        real(dp) :: value = 0, error = -1, tau = 0
    end type

    type :: turbulence_spectrum_t
        !> Axes: kx and ky ascending and symmetric about zero, omega ascending.
        real(dp), allocatable :: kx(:), ky(:), omega(:)
        !> Power spectral density of phi and cross spectral density
        !> <n conjg(phi)> on (kx, ky, omega), per unit dkx dky domega.
        real(dp), allocatable :: s_phi_phi(:, :, :), s_n_phi_re(:, :, :), s_n_phi_im(:, :, :)
        !> Cuts of the correlation function <phi(x + r, t + s) phi(x, t)>
        !> through the origin.
        real(dp), allocatable :: lag_t(:), lag_x(:), lag_y(:)
        real(dp), allocatable :: c_t(:), c_x(:), c_y(:)
        !> Variance of phi inside the retained modes.
        type(turbulence_scalar_t) :: variance
        !> First 1/e crossing of the single-point autocorrelation.
        type(turbulence_scalar_t) :: tau_euler
        !> First 1/e crossing of sum_k |<phi_k(t + s) conjg(phi_k(t))>|: decay
        !> of the mode amplitudes without their phase rotation.
        type(turbulence_scalar_t) :: tau_mode
        !> First 1/e crossings of the correlation cuts in x and y.
        type(turbulence_scalar_t) :: ell_x, ell_y
        !> tau_euler sqrt(<(d phi/dy)^2>)/ell_x.
        type(turbulence_scalar_t) :: kubo
        real(dp) :: dt = 0
        integer :: segment_length = 0, segments = 0, batches = 0
    end type

    type :: turbulence_accumulator_t
        integer :: nseg = 0, nkx = 0, nky = 0, fill = 0
        integer :: segments = 0, slots = 0, per_slot = 1, in_slot = 0
        real(dp) :: dt = 0, window_power = 0
        real(dp), allocatable :: kx(:), ky(:), window(:)
        complex(dp), allocatable :: buffer(:, :, :, :)
        !> Summed periodograms of phi per batch slot and of the cross term.
        real(dp), allocatable :: power(:, :, :, :)
        complex(dp), allocatable :: cross(:, :, :)
        type(fortnum_fft_plan_t) :: plan
    end type
contains
    !> Mean of a stationary series sampled at spacing dt, with the standard
    !> error from its integrated autocorrelation time.
    function turbulence_series_scalar(x, dt) result(scalar)
        real(dp), intent(in) :: x(:), dt
        type(turbulence_scalar_t) :: scalar
        real(dp) :: tau, ess

        call series_autocorrelation(x, scalar%value, scalar%error, tau, ess)
        scalar%tau = tau*dt
        if (size(x) < 4) scalar%error = -1
    end function

    !> kx: ascending, symmetric about zero, equally spaced. ky: 0, dky, 2 dky, ...
    !> nseg: even segment length in samples of spacing dt.
    subroutine turbulence_accumulator_init(acc, kx, ky, dt, nseg)
        type(turbulence_accumulator_t), intent(out) :: acc
        real(dp), intent(in) :: kx(:), ky(:), dt
        integer, intent(in) :: nseg
        integer :: t

        if (nseg < 8 .or. mod(nseg, 4) /= 0) error stop 'turbulence: nseg must be a multiple of 4'
        if (size(kx) < 3 .or. mod(size(kx), 2) /= 1 .or. size(ky) < 2) &
            error stop 'turbulence: need an odd number of kx and at least two ky'
        if (abs(kx(1) + kx(size(kx))) > 1e-12_dp*abs(kx(1)) .or. abs(ky(1)) > 0) &
            error stop 'turbulence: kx must be symmetric and ky must start at zero'
        if (.not. dt > 0) error stop 'turbulence: need dt > 0'
        acc%nseg = nseg; acc%nkx = size(kx); acc%nky = size(ky); acc%dt = dt
        acc%kx = kx; acc%ky = ky
        allocate (acc%window(nseg), acc%buffer(nseg, acc%nkx, acc%nky, 2))
        allocate (acc%power(nseg, acc%nkx, acc%nky, batch_slots))
        allocate (acc%cross(nseg, acc%nkx, acc%nky))
        do t = 1, nseg
            acc%window(t) = sin(pi*(t - 1)/nseg)**2
        end do
        acc%window_power = sum(acc%window**2)
        acc%power = 0; acc%cross = 0
        call fft_c2c_plan_init(acc%plan, nseg)
    end subroutine

    subroutine turbulence_accumulate(acc, phi, n)
        type(turbulence_accumulator_t), intent(inout) :: acc
        complex(dp), intent(in) :: phi(:, :), n(:, :)
        integer :: half

        acc%fill = acc%fill + 1
        acc%buffer(acc%fill, :, :, 1) = phi
        acc%buffer(acc%fill, :, :, 2) = n
        if (acc%fill < acc%nseg) return
        call add_segment(acc)
        half = acc%nseg/2
        acc%buffer(1:half, :, :, :) = acc%buffer(half + 1:, :, :, :)
        acc%fill = half
    end subroutine

    subroutine add_segment(acc)
        type(turbulence_accumulator_t), intent(inout) :: acc
        complex(dp) :: x(acc%nseg), y(acc%nseg)
        real(dp) :: scale
        integer :: i, j, slot, half

        if (acc%in_slot == 0) then
            if (acc%slots == batch_slots) then
                do slot = 1, batch_slots/2
                    acc%power(:, :, :, slot) = acc%power(:, :, :, 2*slot - 1) &
                                               + acc%power(:, :, :, 2*slot)
                end do
                acc%power(:, :, :, batch_slots/2 + 1:) = 0
                acc%slots = batch_slots/2
                acc%per_slot = 2*acc%per_slot
            end if
            acc%slots = acc%slots + 1
        end if
        scale = 1/(acc%nseg*acc%window_power)
        half = acc%nseg/2
        !$omp parallel do collapse(2) private(i, j, x, y) schedule(static)
        do j = 1, acc%nky
            do i = 1, acc%nkx
                ! Backward transform: bin omega multiplies exp(-i omega t).
                x = acc%window*acc%buffer(:, i, j, 1)
                y = acc%window*acc%buffer(:, i, j, 2)
                call fft_c2c(x, 1, acc%plan)
                call fft_c2c(y, 1, acc%plan)
                x = cshift(x, half); y = cshift(y, half)
                acc%power(:, i, j, acc%slots) = acc%power(:, i, j, acc%slots) &
                                                + scale*(real(x, dp)**2 + aimag(x)**2)
                acc%cross(:, i, j) = acc%cross(:, i, j) + scale*y*conjg(x)
            end do
        end do
        !$omp end parallel do
        acc%segments = acc%segments + 1
        acc%in_slot = acc%in_slot + 1
        if (acc%in_slot == acc%per_slot) acc%in_slot = 0
    end subroutine

    subroutine turbulence_spectrum(acc, spectrum)
        type(turbulence_accumulator_t), intent(in) :: acc
        type(turbulence_spectrum_t), intent(out) :: spectrum
        real(dp), allocatable :: total(:, :, :), batch(:, :)
        real(dp) :: cell, scalars(6), dkx, dky, domega
        integer :: nseg, nkx, nky, nlag, i, j, w, jm, im, wm, complete, slot

        nseg = acc%nseg; nkx = acc%nkx; nky = acc%nky
        if (acc%segments < 1) error stop 'turbulence_spectrum: no complete segment'
        dkx = acc%kx(2) - acc%kx(1); dky = acc%ky(2) - acc%ky(1)
        domega = 2*pi/(nseg*acc%dt)
        cell = dkx*dky*domega
        spectrum%dt = acc%dt
        spectrum%segment_length = nseg
        spectrum%segments = acc%segments
        spectrum%kx = acc%kx
        spectrum%ky = [(-acc%ky(nky + 1 - j), j=1, nky - 1), acc%ky]
        spectrum%omega = [(domega*(w - 1 - nseg/2), w=1, nseg)]
        nlag = nseg/4
        spectrum%lag_t = [(acc%dt*w, w=0, nlag)]
        spectrum%lag_x = [(pi/dkx*i/(2*nkx), i=0, 2*nkx)]
        spectrum%lag_y = [(pi/dky*j/(4*nky), j=0, 4*nky)]

        total = sum(acc%power(:, :, :, 1:acc%slots), dim=4)/acc%segments
        allocate (spectrum%s_phi_phi(nkx, 2*nky - 1, nseg))
        allocate (spectrum%s_n_phi_re(nkx, 2*nky - 1, nseg))
        allocate (spectrum%s_n_phi_im(nkx, 2*nky - 1, nseg))
        do j = 1, nky
            do i = 1, nkx
                spectrum%s_phi_phi(i, nky - 1 + j, :) = total(:, i, j)/cell
                spectrum%s_n_phi_re(i, nky - 1 + j, :) = &
                    real(acc%cross(:, i, j), dp)/(acc%segments*cell)
                spectrum%s_n_phi_im(i, nky - 1 + j, :) = &
                    aimag(acc%cross(:, i, j))/(acc%segments*cell)
            end do
        end do
        ! Real fields: X(-k, -omega) = conjg(X(k, omega)).
        do j = 1, nky - 1
            jm = 2*nky - j
            do i = 1, nkx
                im = nkx + 1 - i
                do w = 1, nseg
                    wm = modulo(nseg + 1 - w, nseg) + 1
                    spectrum%s_phi_phi(i, j, w) = spectrum%s_phi_phi(im, jm, wm)
                    spectrum%s_n_phi_re(i, j, w) = spectrum%s_n_phi_re(im, jm, wm)
                    spectrum%s_n_phi_im(i, j, w) = -spectrum%s_n_phi_im(im, jm, wm)
                end do
            end do
        end do

        allocate (spectrum%c_t(0:nlag), spectrum%c_x(0:2*nkx), spectrum%c_y(0:4*nky))
        call correlation_scalars(acc, spectrum, total, scalars, spectrum%c_t, &
                                 spectrum%c_x, spectrum%c_y)
        spectrum%variance%value = scalars(1)
        spectrum%tau_euler%value = scalars(2)
        spectrum%tau_mode%value = scalars(3)
        spectrum%ell_x%value = scalars(4)
        spectrum%ell_y%value = scalars(5)
        spectrum%kubo%value = scalars(6)

        complete = acc%slots
        if (acc%in_slot /= 0) complete = complete - 1
        spectrum%batches = complete
        if (complete < 2) return
        allocate (batch(6, complete))
        block
            real(dp) :: ct(0:nlag), cx(0:2*nkx), cy(0:4*nky)
            do slot = 1, complete
                total = acc%power(:, :, :, slot)/acc%per_slot
                call correlation_scalars(acc, spectrum, total, batch(:, slot), ct, cx, cy)
            end do
        end block
        spectrum%variance%error = batch_error(batch(1, :))
        spectrum%tau_euler%error = batch_error(batch(2, :))
        spectrum%tau_mode%error = batch_error(batch(3, :))
        spectrum%ell_x%error = batch_error(batch(4, :))
        spectrum%ell_y%error = batch_error(batch(5, :))
        spectrum%kubo%error = batch_error(batch(6, :))
    end subroutine

    !> Standard error of the mean of independent batch values; unavailable
    !> when a batch did not reach its 1/e crossing.
    pure real(dp) function batch_error(x)
        real(dp), intent(in) :: x(:)
        real(dp) :: mean
        batch_error = -1
        if (any(x < 0)) return
        mean = sum(x)/size(x)
        batch_error = sqrt(sum((x - mean)**2)/(size(x)*(size(x) - 1)))
    end function

    !> scalars: variance, tau_euler, tau_mode, ell_x, ell_y, kubo from the
    !> mode power spectrum power(omega, kx, ky >= 0) (not yet per cell).
    subroutine correlation_scalars(acc, spectrum, power, scalars, c_t, c_x, c_y)
        type(turbulence_accumulator_t), intent(in) :: acc
        type(turbulence_spectrum_t), intent(in) :: spectrum
        real(dp), intent(in) :: power(:, :, :)
        real(dp), intent(out) :: scalars(6), c_t(0:), c_x(0:), c_y(0:)
        complex(dp) :: mode(0:ubound(c_t, 1))
        complex(dp), allocatable :: rotation(:, :)
        real(dp) :: amplitude(0:ubound(c_t, 1)), weight, energy, taper, velocity2
        integer :: i, j, w, lag

        allocate (rotation(acc%nseg, 0:ubound(c_t, 1)))
        do lag = 0, ubound(c_t, 1)
            do w = 1, acc%nseg
                rotation(w, lag) = exp(cmplx(0, -spectrum%omega(w)*spectrum%lag_t(lag + 1), dp))
            end do
        end do
        c_t = 0; c_x = 0; c_y = 0; amplitude = 0; velocity2 = 0
        do j = 1, acc%nky
            ! The half plane ky > 0 stands for its mirror image as well.
            weight = merge(1.0_dp, 2.0_dp, j == 1)
            do i = 1, acc%nkx
                energy = sum(power(:, i, j))
                do lag = 0, ubound(c_t, 1)
                    mode(lag) = sum(power(:, i, j)*rotation(:, lag))
                end do
                c_t = c_t + weight*real(mode, dp)
                amplitude = amplitude + weight*abs(mode)
                c_x = c_x + weight*energy*cos(acc%kx(i)*spectrum%lag_x)
                c_y = c_y + weight*energy*cos(acc%ky(j)*spectrum%lag_y)
                velocity2 = velocity2 + weight*energy*acc%ky(j)**2
            end do
        end do
        ! Undo the taper that the window puts on the lagged products.
        do lag = 0, ubound(c_t, 1)
            taper = sum(acc%window(1:acc%nseg - lag)*acc%window(1 + lag:acc%nseg)) &
                    /acc%window_power
            c_t(lag) = c_t(lag)/taper
            amplitude(lag) = amplitude(lag)/taper
        end do
        scalars(1) = c_t(0)
        scalars(2) = first_crossing(spectrum%lag_t, c_t)
        scalars(3) = first_crossing(spectrum%lag_t, amplitude)
        scalars(4) = first_crossing(spectrum%lag_x, c_x)
        scalars(5) = first_crossing(spectrum%lag_y, c_y)
        scalars(6) = -1
        if (scalars(2) > 0 .and. scalars(4) > 0) &
            scalars(6) = scalars(2)*sqrt(velocity2)/scalars(4)
    end subroutine

    !> Lag at which c/c(0) first falls below 1/e, by linear interpolation.
    pure real(dp) function first_crossing(lag, c)
        real(dp), intent(in) :: lag(:), c(:)
        real(dp) :: level
        integer :: k

        first_crossing = -1
        if (.not. c(1) > 0) return
        level = c(1)*exp(-1.0_dp)
        do k = 2, size(c)
            if (c(k) < level) then
                first_crossing = lag(k - 1) + (lag(k) - lag(k - 1))*(c(k - 1) - level) &
                                 /(c(k - 1) - c(k))
                return
            end if
        end do
    end function
end module milca_turbulence_statistics
