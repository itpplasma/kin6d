!> Electromagnetic wave at a plasma edge (docs/em-slab-reflection.md).
!> Modes (one CTest each): oracles | vacuum | reflection | transmission.
!> The linear-theory references are computed here independently with complex
!> arithmetic; the generated oracle leaves are checked against them.
program test_em_slab
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_em_slab, only: em_slab_config_t, em_slab_result_t, em_slab_run, &
        em_slab_t, em_slab_init, em_slab_step, em_slab_gauss_residual, &
        em_slab_reflection, em_slab_decay_length, em_slab_mean_amplitude, &
        em_slab_wavenumber
    use kin6d_em_slab_leaves, only: em_slab_boris, em_slab_boris_jvp, &
        em_slab_boris_vjp, em_slab_fresnel, em_slab_evanescent
    implicit none
    character(len=32) :: mode
    real(dp), parameter :: pi = 3.14159265358979323846_dp
    if (command_argument_count() /= 1) error stop 'usage: test_em_slab MODE'
    call get_command_argument(1, mode)
    select case (trim(mode))
    case ('oracles'); call check_oracles()
    case ('vacuum'); call check_vacuum()
    case ('reflection'); call check_reflection()
    case ('transmission'); call check_transmission()
    case ('refinement'); call check_refinement()
    case default; error stop 'unknown mode'
    end select
    print '(2a)', 'em_slab ', trim(mode)//' passed'

contains

    !> Default case geometry with the measurement window used by the demo.
    subroutine run_case(omega, plasma, start, periods, result, dx)
        real(dp), intent(in) :: omega, start
        logical, intent(in) :: plasma
        integer, intent(in) :: periods
        type(em_slab_result_t), intent(out) :: result
        real(dp), intent(in), optional :: dx
        type(em_slab_config_t) :: config
        config%omega = omega; config%plasma = plasma
        if (present(dx)) config%dx = dx
        call em_slab_run(config, start, periods, -16.0_dp, -12.0_dp, result)
    end subroutine

    subroutine check_oracles()
        real(dp) :: ux, uy, a, b, n, r, t, k, lambda, kappa, delta, rre, rim, phase
        real(dp) :: tabs, d(6), g(6), du(2), gu(2), plus(2), minus(2), eps, lhs, rhs
        complex(dp) :: nc, rc, tc
        real(dp) :: omegas(3) = [0.3_dp, 0.5_dp, 0.9_dp]
        integer :: i
        type(em_slab_t) :: s
        type(em_slab_config_t) :: config
        ! Pure magnetic rotation preserves |u| (hence gamma) to rounding.
        call em_slab_boris(0.3_dp, -0.4_dp, 0.0_dp, 0.0_dp, 2.7_dp, 0.05_dp, a, b)
        call require(abs(hypot(a, b) - 0.5_dp) < 4*epsilon(1.0_dp), 'Boris |u|')
        ! Uniform E_y kick from rest, B = 0: u_y = -E_y dt exactly.
        call em_slab_boris(0.0_dp, 0.0_dp, 0.0_dp, 0.2_dp, 0.0_dp, 0.05_dp, a, b)
        call require(abs(b + 0.02_dp) < 4*epsilon(1.0_dp) .and. a == 0, 'Boris kick')
        ! JVP against central differences, VJP against JVP (Euclidean pairing).
        ux = 0.11_dp; uy = -0.07_dp
        d = [0.3_dp, -0.2_dp, 0.5_dp, 0.1_dp, -0.4_dp, 0.25_dp]
        g = [0.17_dp, -0.31_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
        call em_slab_boris_jvp(ux, uy, 0.2_dp, -0.3_dp, 1.4_dp, 0.05_dp, &
                               d(1), d(2), d(3), d(4), d(5), d(6), du(1), du(2))
        eps = 1.0e-6_dp
        call em_slab_boris(ux + eps*d(1), uy + eps*d(2), 0.2_dp + eps*d(3), &
                           -0.3_dp + eps*d(4), 1.4_dp + eps*d(5), 0.05_dp + eps*d(6), &
                           plus(1), plus(2))
        call em_slab_boris(ux - eps*d(1), uy - eps*d(2), 0.2_dp - eps*d(3), &
                           -0.3_dp - eps*d(4), 1.4_dp - eps*d(5), 0.05_dp - eps*d(6), &
                           minus(1), minus(2))
        call require(maxval(abs((plus - minus)/(2*eps) - du)) < 1.0e-8_dp, 'Boris JVP')
        call em_slab_boris_vjp(ux, uy, 0.2_dp, -0.3_dp, 1.4_dp, 0.05_dp, g(1), g(2), &
                               gu(1), gu(2), a, b, n, r)
        lhs = g(1)*du(1) + g(2)*du(2)
        rhs = gu(1)*d(1) + gu(2)*d(2) + a*d(3) + b*d(4) + n*d(5) + r*d(6)
        call require(abs(lhs - rhs) < 1.0e-14_dp, 'Boris VJP pairing')
        ! Generated linear-theory leaves against complex Fresnel matching.
        call em_slab_fresnel(2.0_dp, n, r, t, k, lambda)
        nc = sqrt(cmplx(1 - 1/2.0_dp**2, 0, dp))
        call require(abs(r - real((1 - nc)/(1 + nc))) < 1.0e-15_dp, 'Fresnel R')
        call require(abs(t - real(2/(1 + nc))) < 1.0e-15_dp, 'Fresnel T')
        call require(abs(r**2 + n*t**2 - 1) < 1.0e-15_dp, 'Fresnel flux balance')
        call require(abs(lambda - 2*pi/sqrt(3.0_dp)) < 1.0e-14_dp, 'wavelength')
        do i = 1, size(omegas)
            call em_slab_evanescent(omegas(i), kappa, delta, rre, rim, phase, tabs)
            nc = sqrt(cmplx(1 - 1/omegas(i)**2, 0, dp))     ! i kappa
            rc = (1 - nc)/(1 + nc); tc = 1 + rc
            call require(abs(rre - real(rc)) + abs(rim - aimag(rc)) < 1.0e-14_dp, &
                         'evanescent R')
            call require(abs(phase - atan2(aimag(rc), real(rc))) < 1.0e-14_dp, &
                         'evanescent phase')
            call require(abs(tabs - abs(tc)) < 1.0e-14_dp, 'evanescent |T|')
            call require(abs(delta - 1/(omegas(i)*aimag(nc))) < 1.0e-12_dp, &
                         'skin depth')
        end do
        ! Charge-conserving deposit keeps the discrete Gauss law at rounding
        ! level through a strongly driven short run.
        config%omega = 0.8_dp; config%amplitude = 0.3_dp; config%ramp = 2
        config%x_left = -4; config%x_source = -3; config%x_edge = 0
        config%x_right = 12
        call em_slab_init(s, config)
        do i = 1, 400
            call em_slab_step(s)
        end do
        print '(a,es10.3)', 'Gauss residual after 400 driven steps: ', &
            em_slab_gauss_residual(s)
        call require(em_slab_gauss_residual(s) < 1.0e-11_dp, 'discrete Gauss law')
        call require(maxval(abs(s%uy)) > 0.05_dp, 'driven run moved electrons')
    end subroutine

    !> No plasma: the forward wave travels at c and nothing comes back.
    subroutine check_vacuum()
        type(em_slab_result_t) :: res
        real(dp) :: dphi, speed
        call run_case(2.0_dp, .false., 50.0_dp, 4, res)
        dphi = modulo(atan2(aimag(res%forward_far/res%forward), &
                            real(res%forward_far/res%forward)), 2*pi)
        ! Probes 4 apart: omega D/c = 8 rad = 2 pi + 1.717 rad.
        speed = 2.0_dp*(res%probe_far - res%probe)/(dphi + 2*pi)
        print '(a,f12.8,a,es10.3,a,es10.3)', 'vacuum phase speed c = ', speed, &
            ', |backward/forward| = ', abs(res%backward/res%forward), &
            ', energy balance = ', res%balance
        call require(abs(speed - 1) < 2.0e-4_dp, 'vacuum phase speed')
        call require(abs(abs(res%forward)/5.0e-3_dp - 1) < 1.0e-2_dp, &
                     'incident amplitude')
        call require(abs(res%backward/res%forward) < 2.0e-3_dp, 'no vacuum reflection')
        call require(abs(res%balance) < 2.0e-3_dp, 'vacuum energy balance')
    end subroutine

    !> omega = 0.5 omega_pe: |R| = 1, evanescent skin depth, energy balance.
    subroutine check_reflection()
        type(em_slab_result_t) :: res
        real(dp) :: kappa, delta, rre, rim, phase, tabs, measured_delta, rmeas
        real(dp) :: phase_meas
        complex(dp) :: ratio
        call em_slab_evanescent(0.5_dp, kappa, delta, rre, rim, phase, tabs)
        call run_case(0.5_dp, .true., 75.0_dp, 4, res)
        ratio = em_slab_reflection(res, 0.5_dp, 0.0_dp)
        rmeas = abs(ratio); phase_meas = atan2(aimag(ratio), real(ratio))
        measured_delta = em_slab_decay_length(res, 0.25_dp, 3.0_dp)
        print '(a,f10.6,a,f10.6,a,f10.6)', '|R| = ', rmeas, '  arg R = ', phase_meas, &
            '  (theory ', phase, ')'
        print '(a,f10.6,a,f10.6)', 'skin depth = ', measured_delta, '  theory ', delta
        print '(a,f10.6,a,f10.6)', '|E(edge)|/A = ', &
            abs(res%amplitude(lbound(res%x, 1) + nint(20/0.04_dp)))/abs(res%forward), &
            '  theory ', tabs
        call report_energy(res)
        call require(abs(rmeas - 1) < 5.0e-3_dp, '|R| = 1 below cutoff')
        call require(abs(phase_meas - phase) < 5.0e-3_dp, 'reflection phase')
        call require(abs(measured_delta/delta - 1) < 2.0e-3_dp, 'skin depth')
        call require(abs(res%balance) < 2.0e-3_dp, 'energy balance below cutoff')
        call require(res%gauss_residual < 1.0e-11_dp, 'Gauss law')
    end subroutine

    !> omega = 2 omega_pe: Fresnel R and T, transmitted wavelength, energy.
    subroutine check_transmission()
        type(em_slab_result_t) :: res
        real(dp) :: n, r, t, k, lambda, rmeas, tmeas, kmeas, phase_meas
        complex(dp) :: ratio
        call em_slab_fresnel(2.0_dp, n, r, t, k, lambda)
        call run_case(2.0_dp, .true., 62.0_dp, 8, res)
        ratio = em_slab_reflection(res, 2.0_dp, 0.0_dp)
        rmeas = abs(ratio); phase_meas = atan2(aimag(ratio), real(ratio))
        tmeas = em_slab_mean_amplitude(res, 1.0_dp, 10.0_dp)/abs(res%forward)
        kmeas = em_slab_wavenumber(res, 1.0_dp, 10.0_dp)
        print '(a,f10.6,a,f10.6,a,f10.6)', 'R = ', rmeas, '  theory ', r, &
            '  arg R = ', phase_meas
        print '(a,f10.6,a,f10.6)', 'T = ', tmeas, '  theory ', t
        print '(a,f10.6,a,f10.6)', 'wavelength = ', 2*pi/kmeas, '  theory ', lambda
        print '(a,es10.3)', 'R^2 + N T^2 - 1 (measured) = ', rmeas**2 + n*tmeas**2 - 1
        call report_energy(res)
        call require(abs(rmeas - r) < 2.0e-3_dp, 'Fresnel |R|')
        call require(abs(phase_meas) < 2.0e-2_dp, 'Fresnel R real positive')
        call require(abs(tmeas/t - 1) < 2.0e-3_dp, 'Fresnel T')
        call require(abs(2*pi/kmeas/lambda - 1) < 1.0e-4_dp, 'transmitted wavelength')
        call require(abs(res%balance) < 2.0e-3_dp, 'energy balance above cutoff')
        call require(res%gauss_residual < 1.0e-11_dp, 'Gauss law')
    end subroutine

    !> Halving dx (and dt) must shrink the Fresnel errors: an independent
    !> indicator of the discretization part of the reported deviations.
    subroutine check_refinement()
        type(em_slab_result_t) :: coarse, fine
        real(dp) :: n, r, t, k, lambda, err(2, 2)
        call em_slab_fresnel(2.0_dp, n, r, t, k, lambda)
        call run_case(2.0_dp, .true., 62.0_dp, 8, coarse, 0.08_dp)
        call run_case(2.0_dp, .true., 62.0_dp, 8, fine, 0.04_dp)
        err(:, 1) = [abs(coarse%backward/coarse%forward) - r, &
                     em_slab_mean_amplitude(coarse, 1.0_dp, 10.0_dp)/ &
                     abs(coarse%forward)/t - 1]
        err(:, 2) = [abs(fine%backward/fine%forward) - r, &
                     em_slab_mean_amplitude(fine, 1.0_dp, 10.0_dp)/ &
                     abs(fine%forward)/t - 1]
        print '(a,2es11.3)', 'dx = 0.08: |R| - R, T/T_th - 1 = ', err(:, 1)
        print '(a,2es11.3)', 'dx = 0.04: |R| - R, T/T_th - 1 = ', err(:, 2)
        call require(all(abs(err(:, 2)) < abs(err(:, 1))), 'errors shrink with dx')
    end subroutine

    subroutine report_energy(res)
        type(em_slab_result_t), intent(in) :: res
        print '(a,4es12.4)', 'in, reflected, field, kinetic: ', res%incident, &
            res%reflected, res%field_change, res%kinetic_change
        print '(a,es12.4,a,es10.3,a,es10.3)', 'outflow ', res%outflow, &
            '  balance ', res%balance, '  Gauss ', res%gauss_residual
    end subroutine

    subroutine require(condition, what)
        logical, intent(in) :: condition
        character(*), intent(in) :: what
        if (.not. condition) then
            print '(2a)', 'FAILED: ', what
            error stop 1
        end if
    end subroutine
end program test_em_slab
