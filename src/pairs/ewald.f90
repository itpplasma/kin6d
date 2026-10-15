!> Ewald parameters of the periodic Coulomb sum, chosen from a requested
!> truncation tolerance rather than fixed numbers.
!>
!> Convention (kin6d_symbolic_physics): cube of edge L, neutralizing uniform
!> background, tin-foil boundary, k = 2 pi n/L with n /= 0. The real-space sum
!> runs over images with max-norm index <= real_cutoff after wrapping the
!> separation into [-L/2,L/2]^3; the reciprocal sum over 0 < |n|_inf <=
!> reciprocal_cutoff. Selection minimizes a simple direct-sum cost model.
!>
!> Truncation bound: per unit |g|, the omitted real images and reciprocal
!> vectors are bounded shell by shell (ewald_*_shell_bound leaves, derived in
!> gen/gen_coulomb.f90). Shells are summed until their contribution is below
!> 1e-20 of the running bound; the shell terms decay like exp(-(alpha L m)^2)
!> and exp(-(pi m/(alpha L))^2), so the neglected remainder is far below the
!> reported bound. The bound holds for the pair potential (potential_bound,
!> units |g|/L) and each Cartesian force component (force_bound, |g|/L^2).
module kin6d_ewald
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use kin6d_coulomb_leaves, only: ewald_real_shell_bound, &
                                    ewald_reciprocal_shell_bound
    implicit none
    private
    public :: ewald_parameters_t, ewald_select, ewald_truncation_bound

    type :: ewald_parameters_t
        real(dp) :: length = 0.0_dp, alpha = 0.0_dp
        integer :: real_cutoff = 0, reciprocal_cutoff = 0
        !> Requested relative tolerance and the achieved truncation bounds.
        real(dp) :: tolerance = 0.0_dp
        real(dp) :: potential_bound = 0.0_dp, force_bound = 0.0_dp
    end type
contains
    !> For each real cutoff 1..3: the smallest alpha, then the smallest
    !> reciprocal cutoff, whose truncation bounds are each <= tolerance/2 in
    !> units of |g|/L (potential) and |g|/L^2 (force). Keeps the cheapest by
    !> the pair-count cost model n(n-1)/2 (2 nr+1)^3 + 2 n modes for
    !> `particles` charges (default 2). `real_cutoff` (0..3) restricts the
    !> search to that one cutoff; 0 is the minimum-image real sum (all
    !> omitted images have distance >= L/2), used by single-particle Monte
    !> Carlo updates.
    subroutine ewald_select(length, tolerance, p, ok, particles, real_cutoff)
        real(dp), intent(in) :: length, tolerance
        type(ewald_parameters_t), intent(out) :: p
        logical, intent(out) :: ok
        integer, intent(in), optional :: particles, real_cutoff
        type(ewald_parameters_t) :: trial
        real(dp) :: cost, best
        integer :: n, nr, first, last
        ok = .false.
        if (.not. ieee_is_finite(length) .or. .not. ieee_is_finite(tolerance)) return
        if (length <= 0 .or. tolerance <= 1.0e-15_dp .or. tolerance >= 1) return
        n = 2
        if (present(particles)) n = max(1, particles)
        best = huge(1.0_dp)
        first = 1; last = 3
        if (present(real_cutoff)) then
            if (real_cutoff < 0 .or. real_cutoff > 3) return
            first = real_cutoff; last = real_cutoff
        end if
        do nr = first, last
            call select_for_cutoff(length, tolerance, nr, trial)
            cost = 5*real(n, dp)*(n - 1)/2*real(2*nr + 1, dp)**3 + &
                   real(n, dp)*(real(2*trial%reciprocal_cutoff + 1, dp)**3 - 1)/2
            if (cost < best .and. trial%potential_bound*length <= tolerance .and. &
                trial%force_bound*length**2 <= tolerance) then
                best = cost; p = trial; ok = .true.
            end if
        end do
    end subroutine

    subroutine select_for_cutoff(length, tolerance, nr, p)
        real(dp), intent(in) :: length, tolerance
        integer, intent(in) :: nr
        type(ewald_parameters_t), intent(out) :: p
        real(dp) :: lo, hi, mid, potential, force
        integer :: it
        p%length = length; p%tolerance = tolerance; p%real_cutoff = nr
        lo = 0.1_dp/length; hi = 50.0_dp/length
        do it = 1, 200
            mid = (lo + hi)/2
            call real_tail(mid, length, nr, potential, force)
            if (potential*length <= tolerance/2 .and. &
                force*length**2 <= tolerance/2) then
                hi = mid
            else
                lo = mid
            end if
        end do
        p%alpha = hi
        do it = 1, 64
            p%reciprocal_cutoff = it
            call reciprocal_tail(p%alpha, length, it, potential, force)
            if (potential*length <= tolerance/2 .and. &
                force*length**2 <= tolerance/2) exit
        end do
        call ewald_truncation_bound(p)
    end subroutine

    !> Recompute the truncation bounds (per unit |g|) of given parameters.
    subroutine ewald_truncation_bound(p)
        type(ewald_parameters_t), intent(inout) :: p
        real(dp) :: rp, rf, kp, kf
        call real_tail(p%alpha, p%length, p%real_cutoff, rp, rf)
        call reciprocal_tail(p%alpha, p%length, p%reciprocal_cutoff, kp, kf)
        p%potential_bound = rp + kp
        p%force_bound = rf + kf
    end subroutine

    subroutine real_tail(alpha, length, cutoff, potential, force)
        real(dp), intent(in) :: alpha, length
        integer, intent(in) :: cutoff
        real(dp), intent(out) :: potential, force
        real(dp) :: u, f
        integer :: m
        potential = 0; force = 0
        do m = cutoff + 1, cutoff + 400
            call ewald_real_shell_bound(real(m, dp), alpha, length, u, f)
            potential = potential + u; force = force + f
            if (u <= 1.0e-20_dp*potential .and. f <= 1.0e-20_dp*force) exit
        end do
    end subroutine

    subroutine reciprocal_tail(alpha, length, cutoff, potential, force)
        real(dp), intent(in) :: alpha, length
        integer, intent(in) :: cutoff
        real(dp), intent(out) :: potential, force
        real(dp) :: u, f
        integer :: m
        potential = 0; force = 0
        do m = cutoff + 1, cutoff + 400
            call ewald_reciprocal_shell_bound(real(m, dp), alpha, length, u, f)
            potential = potential + u; force = force + f
            if (u <= 1.0e-20_dp*potential .and. f <= 1.0e-20_dp*force) exit
        end do
    end subroutine
end module kin6d_ewald
