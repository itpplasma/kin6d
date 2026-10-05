!> Static pair-correlation estimators of a periodic one-component snapshot
!> (data for C5/H1): the radial distribution g(r) and the structure factor
!> S(k) = |rho_k|^2/N, rho_k = sum_i exp(-i k.x_i), k = 2 pi n/L, n /= 0.
!>
!> g(r): minimum-image pair counts in shells of width rmax/bins (rmax <= L/2,
!> linked cells of edge >= rmax, each pair once) normalized by the ideal-gas
!> count N (N-1)/2 * shell volume/V. S(k): all lattice vectors with
!> 0 < |n| <= nmax, binned by |n| (exact integer shells n^2). Both are
!> per-snapshot estimators; averaging and intervals belong to the caller.
module kin6d_pair_correlation
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    private
    public :: pair_counts, pair_distribution, structure_factor

    real(dp), parameter :: pi = 3.14159265358979323846_dp, two_pi = 2*pi
contains
    !> Add minimum-image pair counts of the first n particles of x(3n).
    subroutine pair_counts(x, n, length, rmax, counts)
        real(dp), intent(in) :: x(:), length, rmax
        integer, intent(in) :: n
        real(dp), intent(inout) :: counts(:)
        integer, allocatable :: head(:), next(:)
        integer :: nc, i, j, c, c2, a, b, e, oa, ob, oe, q(3), bins
        real(dp) :: d(3), r2, scale
        bins = size(counts); scale = bins/rmax
        nc = int(length/rmax)
        if (nc < 3) then
            do i = 1, n
                do j = i + 1, n
                    call add(i, j)
                end do
            end do
            return
        end if
        allocate (head(nc**3), next(n))
        head = 0
        do i = 1, n
            do a = 1, 3
                q(a) = min(nc - 1, int((x(3*i - 3 + a)/length - &
                                        floor(x(3*i - 3 + a)/length))*nc))
            end do
            c = 1 + q(1) + nc*(q(2) + nc*q(3))
            next(i) = head(c); head(c) = i
        end do
        do e = 0, nc - 1
            do b = 0, nc - 1
                do a = 0, nc - 1
                    c = 1 + a + nc*(b + nc*e)
                    do oe = -1, 1
                        do ob = -1, 1
                            do oa = -1, 1
                                c2 = 1 + modulo(a + oa, nc) + nc*(modulo(b + ob, nc) + &
                                     nc*modulo(e + oe, nc))
                                if (c2 < c) cycle
                                i = head(c)
                                do while (i > 0)
                                    j = head(c2)
                                    if (c2 == c) j = next(i)
                                    do while (j > 0)
                                        call add(i, j)
                                        j = next(j)
                                    end do
                                    i = next(i)
                                end do
                            end do
                        end do
                    end do
                end do
            end do
        end do
    contains
        subroutine add(i, j)
            integer, intent(in) :: i, j
            integer :: k
            d = x(3*i - 2:3*i) - x(3*j - 2:3*j)
            d = d - length*anint(d/length)
            r2 = d(1)*d(1) + d(2)*d(2) + d(3)*d(3)
            if (r2 >= rmax*rmax) return
            k = 1 + int(sqrt(r2)*scale)
            if (k <= bins) counts(k) = counts(k) + 1
        end subroutine
    end subroutine

    !> g at bin centres from accumulated counts over `snapshots` snapshots.
    subroutine pair_distribution(counts, n, length, rmax, snapshots, r, g)
        real(dp), intent(in) :: counts(:), length, rmax
        integer, intent(in) :: n, snapshots
        real(dp), intent(out) :: r(:), g(:)
        real(dp) :: dr, shell, ideal
        integer :: k
        dr = rmax/size(counts)
        do k = 1, size(counts)
            r(k) = (k - 0.5_dp)*dr
            shell = 4*pi/3*((k*dr)**3 - ((k - 1)*dr)**3)
            ideal = real(n, dp)*(n - 1)/2*shell/length**3*snapshots
            g(k) = counts(k)/ideal
        end do
    end subroutine

    !> Shell-averaged S over lattice shells n^2 = 1..nmax^2: s(m) for the
    !> shell m = n^2 (0 where the shell is empty), and its vector count.
    subroutine structure_factor(x, n, length, nmax, s, vectors)
        real(dp), intent(in) :: x(:), length
        integer, intent(in) :: n, nmax
        real(dp), intent(out) :: s(:)
        integer, intent(out) :: vectors(:)
        complex(dp), allocatable :: phase(:, :, :)
        complex(dp) :: rho
        integer :: i, a, m, h, k, l, n2
        real(dp) :: theta
        s = 0; vectors = 0
        allocate (phase(-nmax:nmax, 3, n))
        do i = 1, n
            do a = 1, 3
                theta = two_pi*x(3*i - 3 + a)/length
                do m = -nmax, nmax
                    phase(m, a, i) = cmplx(cos(m*theta), -sin(m*theta), dp)
                end do
            end do
        end do
        do l = -nmax, nmax
            do k = -nmax, nmax
                do h = 0, nmax
                    ! Half space: rho_{-k} = conj(rho_k) has the same |rho|^2.
                    if (h == 0 .and. (k < 0 .or. (k == 0 .and. l <= 0))) cycle
                    n2 = h*h + k*k + l*l
                    if (n2 > min(size(s), nmax*nmax)) cycle
                    rho = (0.0_dp, 0.0_dp)
                    do i = 1, n
                        rho = rho + phase(h, 1, i)*phase(k, 2, i)*phase(l, 3, i)
                    end do
                    s(n2) = s(n2) + (real(rho)**2 + aimag(rho)**2)/n
                    vectors(n2) = vectors(n2) + 1
                end do
            end do
        end do
        where (vectors > 0) s = s/vectors
    end subroutine
end module kin6d_pair_correlation
