!> Independent-method cross-check of the H1a Monte Carlo chains: molecular
!> dynamics of the same OCP (kin6d_ocp: PME + encounter-resolved splitting,
!> plasma units). Preparation: uniform positions, Maxwellian velocities,
!> isokinetic rescaling every 5 steps during equilibration only (declared
!> thermostat); production is thermostat-free NVE. Gamma maps to the
!> plasma-unit density by Lambda = n lambda_D^3 = 1/(sqrt(48) pi Gamma^1.5).
module correlation_map_md
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use fortnum_status, only: fortnum_status_t, FORTNUM_OK
    use kin6d_ocp, only: ocp_model_t, ocp_system, ocp_prepare, ocp_thermostat, &
        ocp_advance, ocp_kinetic_temperature
    use kin6d_nbody, only: nbody_system_t, nbody_energy
    use kin6d_pair_correlation, only: pair_counts, pair_distribution
    use kin6d_series_statistics, only: series_autocorrelation
    implicit none
    private
    public :: md_result_t, md_run

    type :: md_result_t
        real(dp) :: gamma = 0, dt = 0, temperature = 0, temperature_err = 0
        real(dp) :: u = 0, u_err = 0, tau = 0, drift = 0, wall = 0, a = 0
        integer :: n = 0, steps = 0, failed_step = 0
        character(128) :: failure = ''
        real(dp), allocatable :: r(:), g(:)
    end type
contains
    subroutine md_run(gamma, n, t_equil, t_prod, width, seed, md)
        real(dp), intent(in) :: gamma, t_equil, t_prod, width
        integer, intent(in) :: n
        integer(int64), intent(in) :: seed
        type(md_result_t), intent(out) :: md
        type(ocp_model_t) :: model
        type(nbody_system_t) :: system
        type(fortnum_status_t) :: status
        real(dp), allocatable :: state(:), u(:), temp(:), counts(:)
        real(dp) :: pi, lambda, e0, e, ke, err, ess, rmax
        integer :: rejected, i, s, ns, bins
        integer(int64) :: c0, c1, rate
        logical :: ok
        call system_clock(c0, rate)
        pi = acos(-1.0_dp)
        lambda = 1/(sqrt(48.0_dp)*pi*gamma**1.5_dp)
        model%debye_volume = lambda; model%box = (n/lambda)**(1.0_dp/3)
        ! A quarter of the default outer step 0.01 (100/Lambda)^(1/3): at
        ! Gamma = 1 the default step fails closed in the cluster check within
        ! ~1e4 steps and dt/2 drifts by ~1e-2 k_B T per particle per 100/omega_p;
        ! dt/4 drifts by < 1e-3 (measured).
        model%dt = 0.0025_dp*(100/lambda)**(1.0_dp/3)
        call ocp_system(model, system, ok)
        if (.not. ok .or. abs(model%scales%gamma/gamma - 1) > 1.0e-6_dp) &
            error stop 'MD: OCP system'
        md%gamma = gamma; md%n = model%electrons; md%dt = model%dt
        md%a = model%scales%wigner_seitz
        allocate (state(6*model%particles))
        call ocp_prepare(model, seed, state, rejected, status)
        do i = 1, nint(t_equil/model%dt)
            call ocp_advance(system, state, model%dt, status)
            if (status%code /= FORTNUM_OK) error stop 'MD equilibration step'
            if (mod(i, 5) == 0) call ocp_thermostat(model, state)
        end do
        call ocp_thermostat(model, state)
        md%steps = nint(t_prod/model%dt); ns = md%steps/5
        allocate (u(ns), temp(ns))
        bins = int(model%box/(2*md%a)/width); rmax = bins*width*md%a
        allocate (counts(bins), md%r(bins), md%g(bins))
        counts = 0
        call nbody_energy(system, state, e0, status)
        s = 0
        do i = 1, md%steps
            call ocp_advance(system, state, model%dt, status)
            if (status%code /= FORTNUM_OK) then
                ! Fail closed: keep the samples before the failed step.
                md%failed_step = i; md%failure = status%msg
                exit
            end if
            if (mod(i, 5) /= 0) cycle
            s = s + 1
            call nbody_energy(system, state, e, status)
            ke = sum(state(3*model%particles + 1:)**2)/2
            ! Plasma units: k_B T_ref = 1, so beta U/N = U/N at T_ref.
            u(s) = (e - ke)/md%n
            temp(s) = ocp_kinetic_temperature(model, state)
            call pair_counts(state(1:3*md%n), md%n, model%box, rmax, counts)
        end do
        ns = s
        md%drift = (e - e0)/md%n
        call series_autocorrelation(u(:ns), md%u, md%u_err, md%tau, ess)
        call series_autocorrelation(temp(:ns), md%temperature, md%temperature_err, err, ess)
        call pair_distribution(counts, md%n, model%box, rmax, ns, md%r, md%g)
        ! Same canonical 1/N normalization as the MC tables (S(0) = 0).
        md%g = md%g*(md%n - 1)/real(md%n, dp)
        md%r = md%r/md%a
        call system_clock(c1)
        md%wall = real(c1 - c0, dp)/rate
        print '(a,f6.2,a,i0,a,f8.1,a)', 'MD Gamma=', gamma, ' N=', md%n, ' wall ', &
            md%wall, ' s'
    end subroutine
end module correlation_map_md
