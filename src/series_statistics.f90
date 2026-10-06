!> Statistical summaries of correlated Markov-chain time series (H1a).
!> Intervals from these routines are statistical (one standard error), never
!> deterministic enclosures.
!>
!> - series_autocorrelation: integrated autocorrelation time
!>   tau = 1/2 + sum_{t>=1} rho(t) with Sokal's automatic window (smallest
!>   W with W >= c tau(W), c = 6), effective sample size n/(2 tau) and the
!>   standard error sqrt(2 tau var/n) of the mean.
!> - series_batch_means: mean and standard error from nb contiguous batches
!>   (labelled statistical); valid when the batch length >> tau.
module kin6d_series_statistics
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    private
    public :: series_autocorrelation, series_batch_means
contains
    subroutine series_autocorrelation(x, mean, error, tau, ess)
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: mean, error, tau, ess
        real(dp) :: var, c
        integer :: n, t
        n = size(x)
        mean = sum(x)/max(1, n)
        var = sum((x - mean)**2)/max(1, n)
        tau = 0.5_dp; error = 0; ess = n
        if (n < 4 .or. var <= 0) return
        do t = 1, n/2
            c = sum((x(1:n - t) - mean)*(x(1 + t:n) - mean))/((n - t)*var)
            tau = tau + c
            if (t >= 6*tau) exit
        end do
        tau = max(0.5_dp, tau)
        ess = n/(2*tau)
        error = sqrt(2*tau*var/n)
    end subroutine

    subroutine series_batch_means(x, nb, mean, error)
        real(dp), intent(in) :: x(:)
        integer, intent(in) :: nb
        real(dp), intent(out) :: mean, error
        real(dp) :: b(nb)
        integer :: i, len
        len = size(x)/nb
        mean = sum(x)/max(1, size(x)); error = 0
        if (len < 1 .or. nb < 2) return
        do i = 1, nb
            b(i) = sum(x((i - 1)*len + 1:i*len))/len
        end do
        mean = sum(b)/nb
        error = sqrt(sum((b - mean)**2)/(nb*(nb - 1)))
    end subroutine
end module kin6d_series_statistics
