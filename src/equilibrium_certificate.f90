module kin6d_equilibrium_certificate
    ! One certificate interface for every lane (docs/equilibrium-core.md, D4,
    ! migration M3):
    !     ||U - U_h|| <~ M(U_h) / beta_h.
    ! A lane implements M, a computable majorant of the D3 residual in the
    ! dual norm of its displacement space, and beta_h, an inf-sup (or
    ! coercivity) constant of its Jacobian on the gauge complement, together
    ! with the kind of beta_h: none, a practical estimate, or a lower bound.
    ! The core owns the composition bound = M/beta and the claim label:
    !   candidate           M unavailable (the lane failed);
    !   residual verified   M only (no beta_h);
    !   practical estimate  M divided by an estimated beta_h;
    !   certified           M divided by a lower bound of beta_h.
    ! The rigour of M itself (quadrature, linearisation) is the lane's
    ! declared scope; the label does not upgrade it.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use kin6d_equilibrium_state, only: equilibrium_state_t
    implicit none
    private
    public :: equilibrium_certifier_t, equilibrium_certificate_t
    public :: certify, claim_name

    integer, parameter, public :: EQ_CLAIM_CANDIDATE = 0, EQ_CLAIM_RESIDUAL = 1, &
        EQ_CLAIM_PRACTICAL = 2, EQ_CLAIM_CERTIFIED = 3
    integer, parameter, public :: EQ_STABILITY_NONE = 0, EQ_STABILITY_ESTIMATE = 1, &
        EQ_STABILITY_BOUND = 2

    type :: equilibrium_certificate_t
        real(dp) :: majorant = huge(1.0_dp)  ! M, dual norm of the residual
        real(dp) :: stability = 0            ! beta_h (0 when none)
        real(dp) :: bound = huge(1.0_dp)     ! M/beta_h, or M when no beta_h
        integer :: stability_kind = EQ_STABILITY_NONE
        integer :: claim = EQ_CLAIM_CANDIDATE
        integer :: lane_status = 0           ! the lane's own status code
    end type

    ! A lane's implementation of M and beta_h at its current state.
    type, abstract :: equilibrium_certifier_t
    contains
        procedure(certifier_estimate), deferred :: estimate
    end type
    abstract interface
        ! majorant: M; stability, kind: beta_h and its kind (kind NONE leaves
        ! stability unused); status /= 0: M unavailable; lane_status: the
        ! lane's own status code, reported unchanged.
        subroutine certifier_estimate(certifier, majorant, stability, kind, &
                                      status, lane_status)
            import :: equilibrium_certifier_t, dp
            class(equilibrium_certifier_t), intent(inout) :: certifier
            real(dp), intent(out) :: majorant, stability
            integer, intent(out) :: kind, status, lane_status
        end subroutine
    end interface
contains
    subroutine certify(certifier, cert, status, eq)
        ! Evaluate the lane's M and beta_h and assign the claim; with eq, the
        ! claim label is written to its provenance.
        class(equilibrium_certifier_t), intent(inout) :: certifier
        type(equilibrium_certificate_t), intent(out) :: cert
        integer, intent(out) :: status
        type(equilibrium_state_t), intent(inout), optional :: eq
        real(dp) :: m, beta
        integer :: kind
        call certifier%estimate(m, beta, kind, status, cert%lane_status)
        if (status == 0 .and. .not. ieee_is_finite(m)) status = 1
        if (status == 0) then
            cert%majorant = m
            cert%bound = m
            cert%claim = EQ_CLAIM_RESIDUAL
            if (kind /= EQ_STABILITY_NONE .and. ieee_is_finite(beta)) then
                if (beta > 0) then
                    cert%stability = beta
                    cert%stability_kind = kind
                    cert%bound = m/beta
                    cert%claim = EQ_CLAIM_PRACTICAL
                    if (kind == EQ_STABILITY_BOUND) cert%claim = EQ_CLAIM_CERTIFIED
                end if
            end if
        end if
        if (present(eq)) eq%claim = claim_name(cert%claim)
    end subroutine

    pure function claim_name(claim) result(name)
        integer, intent(in) :: claim
        character(32) :: name
        select case (claim)
        case (EQ_CLAIM_RESIDUAL)
            name = 'residual verified'
        case (EQ_CLAIM_PRACTICAL)
            name = 'practical estimate'
        case (EQ_CLAIM_CERTIFIED)
            name = 'certified'
        case default
            name = 'candidate'
        end select
    end function
end module kin6d_equilibrium_certificate
