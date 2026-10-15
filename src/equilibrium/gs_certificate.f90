module kin6d_gs_certificate
    ! Fixed-boundary estimates prefer local RT equilibration. The legacy
    ! functional recovery remains a fallback, and neither sampled quadrature
    ! nor floating-point geometry is silently promoted to a certificate.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use kin6d_equilibrium_certificate, only: equilibrium_certifier_t, &
        EQ_STABILITY_NONE, EQ_STABILITY_ESTIMATE, EQ_STABILITY_BOUND
    use kin6d_grad_shafranov, only: gs_state_t
    use kin6d_grad_shafranov_leaves, only: gs_source
    use kin6d_gs_profile_law, only: gs_law_sources
    use kin6d_gs_equilibrated, only: gs_equilibrated_majorant, &
        equilibrated_flux_estimate_t, GS_EQ_OK
    use kin6d_gs_error_control, only: gs_stability_t, gs_fixed_stability, &
        gs_fixed_source_stability, GS_ERROR_OK
    use kin6d_gs_quadrature_budget, only: gs_checked_equilibrated_majorant, &
        GS_CHECKED_OK
    use fortnum_interval, only: interval_t
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
        ! Optional caller stability lower bound (<= 0: none). Sampled
        ! residual quadrature remains practical even when this bound is known.
        real(dp) :: stability_lower_bound = 0
        logical :: estimate_stability = .true. ! enables automatic nonlinear stability
        character(32) :: method = 'unavailable'
        character(32) :: quadrature_status = 'approximate'
        type(equilibrated_flux_estimate_t) :: last
        type(gs_stability_t) :: stability_result
        real(dp) :: flux_defect = 0, balance_defect = 0
        logical :: recovery_converged = .false.
        integer :: recovery_iterations = 0
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
        real(dp) :: bounded
        type(interval_t) :: energy_bound
        integer :: eq_status, stability_status, checked_status
        majorant = huge(1.0_dp); stability = 0; kind = EQ_STABILITY_NONE
        status = 1; lane_status = 1
        certifier%method = 'unavailable'
        certifier%quadrature_status = 'approximate'
        certifier%recovery_converged = .false.
        certifier%recovery_iterations = 0
        if (.not. associated(certifier%state)) return
        if (certifier%mu0 <= 0) return
        if (.not. all(ieee_is_finite([certifier%mu0, certifier%pprime, &
            certifier%ffprime]))) return
        if (allocated(certifier%law)) then
            if (.not. certifier%law%valid()) return
        end if
        certifier%stability_result = gs_stability_t()
        call gs_equilibrated_majorant(certifier%state, certifier%state%psi, &
            source, certifier%last, eq_status)
        if (eq_status == GS_EQ_OK) then
            certifier%method = 'equilibrated-rt'
            certifier%flux_defect = certifier%last%flux_defect
            certifier%balance_defect = certifier%last%oscillation &
                + certifier%last%divergence_defect
            majorant = certifier%last%bound
            certifier%recovery_converged = .true.
        else
            ! Unsupported cells/sources keep the established recovery path.
            call gs_functional_majorant(certifier%state, certifier%mu0, &
                certifier%pprime, certifier%ffprime, bounded, &
                certifier%flux_defect, certifier%balance_defect, lane_status, &
                law=certifier%law, forcing=certifier%forcing, &
                recovery_converged=certifier%recovery_converged, &
                recovery_iterations=certifier%recovery_iterations, &
                estimate_stability=.false.)
            if (lane_status /= 0 .and. lane_status /= 4) return
            majorant = certifier%flux_defect + certifier%balance_defect
            certifier%method = 'functional-recovery'
        end if
        if (allocated(certifier%law)) then
            stability_status = 1
            if (certifier%estimate_stability) then
                call gs_fixed_stability(certifier%state, certifier%law, certifier%mu0, &
                    certifier%stability_result, stability_status)
            end if
        else
            call gs_fixed_source_stability(certifier%state, &
                certifier%stability_result, stability_status)
        end if
        status = 0
        lane_status = 4 ! frozen-source residual when nonlinear stability is absent
        if (stability_status == GS_ERROR_OK) then
            stability = certifier%stability_result%alpha
            kind = EQ_STABILITY_ESTIMATE
            lane_status = 5
        end if
        if (certifier%stability_lower_bound > 0) then
            stability = certifier%stability_lower_bound
            kind = EQ_STABILITY_ESTIMATE
            lane_status = 5
        end if
        ! Stability alone cannot upgrade sampled residual norms. An independent
        ! outward audit may establish the full continuous polygon-model bound
        ! for affine P1/P2 and the fixed constant-profile source only.
        if (eq_status == GS_EQ_OK .and. .not. allocated(certifier%law)) then
            if (.not. associated(certifier%forcing)) then
                call gs_checked_equilibrated_majorant(certifier%state, &
                    certifier%state%psi, certifier%last, certifier%mu0, &
                    certifier%pprime, certifier%ffprime, energy_bound, checked_status)
                if (checked_status == GS_CHECKED_OK) then
                    majorant = energy_bound%hi
                    stability = 1
                    kind = EQ_STABILITY_BOUND
                    lane_status = 0
                    certifier%method = 'equilibrated-rt0-audit'
                    certifier%quadrature_status = 'rigorous-rt0-audit'
                end if
            end if
        end if
    contains
        function source(cell, r, z, psi) result(f)
            integer, intent(in) :: cell
            real(dp), intent(in) :: r, z, psi
            real(dp) :: f, pp, ff
            if (cell < 1) then
                f = huge(1.0_dp)
                return
            end if
            pp = certifier%pprime; ff = certifier%ffprime
            if (allocated(certifier%law)) then
                call gs_law_sources(certifier%law, psi, pp, ff)
            end if
            call gs_source(r, certifier%mu0, pp, ff, f)
            if (associated(certifier%forcing)) f = f + certifier%forcing(r, z)
        end function source
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
