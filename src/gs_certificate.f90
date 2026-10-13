module kin6d_gs_certificate
    ! The psi lane's implementations of the D4 certificate interface
    ! (docs/equilibrium-core.md, migration M3). Both wrap kin6d_gs_majorant
    ! unchanged and report its majorant and stability constant separately;
    ! the core composes the bound and assigns the claim.
    !
    ! Fixed boundary (gs_functional_majorant, Repin D + C_F B):
    !   status 0 without law  linear problem, the majorant bounds the energy
    !                         error directly: beta = 1, a bound (certified);
    !   status 0 with law     caller lower bound of the branch stability:
    !                         beta = that bound (certified, conditional on
    !                         its stated assumptions);
    !   status 4              frozen-source defect only (residual verified);
    !   status 5              estimated coercivity alpha (practical estimate).
    ! Free boundary (gs_fb_majorant, P10): M = D + B_T + N of the frozen
    ! source, beta the power-iteration inf-sup estimate of the Jacobian in
    ! the coupled energy pairing (practical estimate).
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_equilibrium_certificate, only: equilibrium_certifier_t, &
        EQ_STABILITY_NONE, EQ_STABILITY_ESTIMATE, EQ_STABILITY_BOUND
    use kin6d_grad_shafranov, only: gs_state_t
    use kin6d_gs_profile_law, only: gs_profile_law_t
    use kin6d_gs_nonlinear, only: gs_forcing_t
    use kin6d_gs_free_boundary, only: gs_fb_system_t
    use kin6d_gs_majorant, only: gs_functional_majorant, gs_fb_majorant, &
        gs_coupled_estimator_t, gs_fb_majorant_t, GS_CM_OK
    implicit none
    private
    public :: gs_fixed_certifier_t, gs_fb_certifier_t

    type, extends(equilibrium_certifier_t) :: gs_fixed_certifier_t
        type(gs_state_t), pointer :: state => null()
        real(dp) :: mu0 = 0, pprime = 0, ffprime = 0
        ! Profile law (unallocated: the linear problem with pprime, ffprime).
        type(gs_profile_law_t), allocatable :: law
        procedure(gs_forcing_t), pointer, nopass :: forcing => null()
        ! Caller lower bound of the branch stability (<= 0: none); otherwise
        ! the coercivity is estimated when estimate_stability is set.
        real(dp) :: stability_lower_bound = 0
        logical :: estimate_stability = .true.
    contains
        procedure :: estimate => fixed_estimate
    end type

    type, extends(equilibrium_certifier_t) :: gs_fb_certifier_t
        type(gs_state_t), pointer :: state => null()
        type(gs_fb_system_t), pointer :: sys => null()
        type(gs_coupled_estimator_t), pointer :: est => null()
        real(dp), allocatable :: currents(:)
        procedure(gs_forcing_t), pointer, nopass :: forcing => null()
        integer :: max_iterations = 200
        ! The lane's full result of the last estimate.
        type(gs_fb_majorant_t) :: last
    contains
        procedure :: estimate => fb_estimate
    end type
contains
    subroutine fixed_estimate(certifier, majorant, stability, kind, status, &
                              lane_status)
        class(gs_fixed_certifier_t), intent(inout) :: certifier
        real(dp), intent(out) :: majorant, stability
        integer, intent(out) :: kind, status, lane_status
        real(dp) :: bounded, flux, balance, alpha, lower
        logical :: practical
        majorant = huge(1.0_dp); stability = 0; kind = EQ_STABILITY_NONE
        status = 1; lane_status = 1
        if (.not. associated(certifier%state)) return
        lower = certifier%stability_lower_bound
        practical = certifier%estimate_stability .and. .not. lower > 0
        alpha = 0
        if (lower > 0) then
            if (associated(certifier%forcing)) then
                call gs_functional_majorant(certifier%state, certifier%mu0, &
                    certifier%pprime, certifier%ffprime, bounded, flux, balance, &
                    lane_status, law=certifier%law, stability_lower_bound=lower, &
                    forcing=certifier%forcing)
            else
                call gs_functional_majorant(certifier%state, certifier%mu0, &
                    certifier%pprime, certifier%ffprime, bounded, flux, balance, &
                    lane_status, law=certifier%law, stability_lower_bound=lower)
            end if
        else if (associated(certifier%forcing)) then
            call gs_functional_majorant(certifier%state, certifier%mu0, &
                certifier%pprime, certifier%ffprime, bounded, flux, balance, &
                lane_status, law=certifier%law, forcing=certifier%forcing, &
                estimate_stability=practical, stability_estimate=alpha)
        else
            call gs_functional_majorant(certifier%state, certifier%mu0, &
                certifier%pprime, certifier%ffprime, bounded, flux, balance, &
                lane_status, law=certifier%law, estimate_stability=practical, &
                stability_estimate=alpha)
        end if
        ! The majorant before any division, as gs_functional_majorant forms it.
        majorant = flux + balance
        select case (lane_status)
        case (0)
            status = 0
            kind = EQ_STABILITY_BOUND
            stability = 1
            if (allocated(certifier%law)) stability = lower
        case (4)
            status = 0
        case (5)
            status = 0
            kind = EQ_STABILITY_ESTIMATE
            stability = alpha
        end select
    end subroutine

    subroutine fb_estimate(certifier, majorant, stability, kind, status, &
                           lane_status)
        class(gs_fb_certifier_t), intent(inout) :: certifier
        real(dp), intent(out) :: majorant, stability
        integer, intent(out) :: kind, status, lane_status
        majorant = huge(1.0_dp); stability = 0; kind = EQ_STABILITY_NONE
        status = 1; lane_status = 1
        if (.not. (associated(certifier%state) .and. associated(certifier%sys))) return
        if (.not. (associated(certifier%est) .and. allocated(certifier%currents))) return
        if (associated(certifier%forcing)) then
            call gs_fb_majorant(certifier%state, certifier%sys, certifier%est, &
                certifier%currents, certifier%last, lane_status, &
                certifier%forcing, certifier%max_iterations)
        else
            call gs_fb_majorant(certifier%state, certifier%sys, certifier%est, &
                certifier%currents, certifier%last, lane_status, &
                max_iterations=certifier%max_iterations)
        end if
        if (lane_status /= GS_CM_OK) return
        status = 0
        majorant = certifier%last%frozen%majorant
        stability = certifier%last%beta
        kind = EQ_STABILITY_ESTIMATE
    end subroutine
end module kin6d_gs_certificate
