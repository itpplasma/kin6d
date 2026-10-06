!> Short-range semiclassical pair term of the periodic electron cell
!> (docs/quantum-classical-correlations.md): the declared addition to the
!> one Ewald Coulomb pair energy, so the long-range part keeps a single
!> definition (kin6d_ewald/kin6d_coulomb_leaves).
!>
!> pair_coulomb       no addition (bare Ewald Coulomb);
!> pair_kelbg         periodic image sum of the Kelbg diffraction excess
!>                    Kelbg(r) - g/r;
!> pair_kelbg_pauli   Kelbg excess plus the spin-averaged Pauli potential of
!>                    the periodic leading-order exchange factor.
!>
!> Every image |d + nL| beyond the image range carries Kelbg/exchange terms
!> below exp(-36) relative to their contact values. Leaves come from
!> kin6d_quantum_classical_leaves (one symbolic source).
module kin6d_pair_law
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_quantum_classical_leaves, only: qcm_kelbg, qcm_kelbg_series, &
        qcm_kelbg_excess, qcm_exchange, qcm_pauli
    implicit none
    private
    public :: pair_law_t, pair_law_init, pair_law_excess, pair_law_exchange_sum
    integer, parameter, public :: pair_coulomb = 1, pair_kelbg = 2, &
        pair_kelbg_pauli = 3

    type :: pair_law_t
        integer :: model = pair_coulomb, images = 0
        !> Cell edge, Coulomb coupling g, pair thermal length lambda, kT and
        !> the same-spin probability of a pair (Pauli term).
        real(dp) :: length = 0, coupling = 1, lambda = 0, kt = 0, same_spin = 0.5_dp
    end type
contains
    !> Rejects unknown models and nonpositive scales of the requested model.
    subroutine pair_law_init(law, model, length, coupling, lambda, kt, same_spin, ok)
        type(pair_law_t), intent(out) :: law
        integer, intent(in) :: model
        real(dp), intent(in) :: length, coupling, lambda, kt, same_spin
        logical, intent(out) :: ok
        ok = .false.
        if (model < pair_coulomb .or. model > pair_kelbg_pauli) return
        if (.not. (length > 0)) return
        if (model /= pair_coulomb .and. .not. (lambda > 0 .and. kt > 0 .and. &
                                               same_spin >= 0 .and. same_spin <= 1)) return
        law%model = model; law%length = length; law%coupling = coupling
        law%lambda = lambda; law%kt = kt; law%same_spin = same_spin
        law%images = 0
        if (lambda > 0) law%images = ceiling((6*lambda + length/2)/length)
        ok = .true.
    end subroutine

    !> Pair addition for the minimum-image separation d (|d_a| <= L/2): the
    !> minimum image uses Kelbg(r) - g/r (series near contact), every other
    !> image the Kelbg excess; then the Pauli term of the exchange factor.
    real(dp) function pair_law_excess(law, d) result(u)
        type(pair_law_t), intent(in) :: law
        real(dp), intent(in) :: d(3)
        real(dp) :: r, v, xs, e(3)
        integer :: a, b, c
        u = 0
        if (law%model == pair_coulomb) return
        r = norm2(d)
        if (r < 0.1_dp*law%lambda) then
            call qcm_kelbg_series(r, law%lambda, law%coupling, v)
        else
            call qcm_kelbg(r, law%lambda, law%coupling, v)
        end if
        u = v - law%coupling/r
        xs = 0
        do a = -law%images, law%images
            do b = -law%images, law%images
                do c = -law%images, law%images
                    e = d + law%length*[a, b, c]
                    r = norm2(e)
                    if (a /= 0 .or. b /= 0 .or. c /= 0) then
                        call qcm_kelbg_excess(r, law%lambda, law%coupling, v)
                        u = u + v
                    end if
                    call qcm_exchange(r, law%lambda, v)
                    xs = xs + v
                end do
            end do
        end do
        if (law%model == pair_kelbg_pauli) then
            call qcm_pauli(xs, law%kt, law%same_spin, v)
            u = u + v
        end if
    end function

    !> Periodic leading-order exchange factor sum_n exp(-|d + nL|^2/lambda^2)
    !> for the minimum-image separation d.
    real(dp) function pair_law_exchange_sum(law, d) result(xs)
        type(pair_law_t), intent(in) :: law
        real(dp), intent(in) :: d(3)
        real(dp) :: v
        integer :: a, b, c
        xs = 0
        do a = -law%images, law%images
            do b = -law%images, law%images
                do c = -law%images, law%images
                    call qcm_exchange(norm2(d + law%length*[a, b, c]), law%lambda, v)
                    xs = xs + v
                end do
            end do
        end do
    end function
end module kin6d_pair_law
