!> NTV field line of an axisymmetric Boozer `.bc` surface (kin6d_boozer_file)
!! as a dk_ntv provider: |B|, d|B|/dtheta and d|B|/ds along theta_B, and the
!! flux functions, through the FortSym leaf ntv_boozer_point. The kernel's
!! precession prefactor for this line is m v^2/(2 e) [V], b = |B|/B_ref with
!! B_ref the (0, 0) harmonic. Only n = 0 rows enter; a file with n /= 0 rows
!! is refused.
module kin6d_ntv_boozer_line
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_dk_ntv, only: ntv_line_t
    use kin6d_boozer_file, only: boozer_surface_t
    use kin6d_ntv_boozer_leaves, only: ntv_boozer_point
    implicit none
    private
    public :: ntv_boozer_line_t, ntv_boozer_line, NTV_LINE_OK, NTV_LINE_REFUSED

    integer, parameter :: NTV_LINE_OK = 0, NTV_LINE_REFUSED = 1

    type, extends(ntv_line_t) :: ntv_boozer_line_t
        type(boozer_surface_t) :: surf
        real(dp) :: bref = 0
    contains
        procedure :: sample => boozer_sample
    end type

contains

    subroutine ntv_boozer_line(surf, line, status)
        type(boozer_surface_t), intent(in) :: surf
        type(ntv_boozer_line_t), intent(out) :: line
        integer, intent(out) :: status
        status = NTV_LINE_REFUSED
        if (any(surf%n /= 0)) return
        if (.not. any(surf%m == 0)) return
        line%surf = surf
        line%bref = sum(surf%bmn, mask=surf%m == 0)
        status = NTV_LINE_OK
    end subroutine ntv_boozer_line

    subroutine boozer_sample(self, theta, stream, mirror, drift, drift2, source, measure)
        class(ntv_boozer_line_t), intent(in) :: self
        real(dp), intent(in) :: theta(:)
        real(dp), intent(out) :: stream(:), mirror(:), drift(:), drift2(:), &
            source(:), measure(:)
        real(dp) :: c, sn, bmod, bmod_t, bmod_s
        integer :: it, k
        associate (sf => self%surf)
            do it = 1, size(theta)
                bmod = 0; bmod_t = 0; bmod_s = 0
                do k = 1, size(sf%m)
                    c = cos(sf%m(k)*theta(it)); sn = sin(sf%m(k)*theta(it))
                    bmod = bmod + sf%bmn(k)*c + sf%bma(k)*sn
                    bmod_t = bmod_t + sf%m(k)*(sf%bma(k)*c - sf%bmn(k)*sn)
                    bmod_s = bmod_s + sf%bmn_s(k)*c + sf%bma_s(k)*sn
                end do
                call ntv_boozer_point(bmod, bmod_t, bmod_s, self%bref, sf%iota, &
                    sf%diota_ds, sf%g_cov, sf%i_cov, sf%g_cov_s, sf%i_cov_s, &
                    sf%psi_prime, stream(it), mirror(it), drift(it), drift2(it), &
                    source(it), measure(it))
            end do
        end associate
    end subroutine boozer_sample

end module kin6d_ntv_boozer_line
