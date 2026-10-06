!> Classical electron one-component plasma (OCP) in a periodic 3D3V cube with
!> a uniform neutralizing background and bare Coulomb interaction (C3, C4).
!>
!> Plasma units of the responsive species (electrons only): k_B T = m = 1,
!> lambda_D = 1, omega_p = 1, so n = Lambda = n lambda_D^3 and the Coulomb
!> constant is 1/(4 pi Lambda). The values are checked against the shared
!> kin6d_species scales (one definition). Optional fixed repulsive impurity of
!> charge Q (electron charge -1, same sign) at the origin, compensated in the
!> background; kappa = |Q|/(4 pi Lambda) = r_nl/lambda_D is its nonlinear
!> radius |qQ|/(4 pi eps0 k_B T) in Debye lengths. Interactions: the shared
!> N-body system with the particle-mesh Ewald action (kin6d_pme) and the
!> encounter-resolved near/far step (kin6d_near_encounter). No softening,
!> screened potential or collision operator is imposed.
module kin6d_ocp
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortnum_status, only: fortnum_status_t, FORTNUM_OK, FORTNUM_DOMAIN_ERROR, &
                              status_set
    use fortnum_rng, only: rng_t, rng_seed, rng_uniform, rng_normal, rng_split
    use kin6d_species, only: units_t, species_t, reduced_units, electron, &
        plasma_scales, plasma_scales_t, electron_mass_si, elementary_charge_si
    use kin6d_pme, only: pme_parameters_t, pme_select
    use kin6d_ocp_leaves, only: debye_image, debye_zero_mode, debye_shell_bound
    use kin6d_nbody, only: nbody_system_t, nbody_periodic_mesh, nbody_fix, &
        nbody_step, nbody_step_jvp, nbody_set_encounter
    implicit none
    private
    public :: ocp_model_t, ocp_units, ocp_system, ocp_prepare, ocp_thermostat
    public :: ocp_advance, ocp_advance_jvp, ocp_kinetic_temperature, ocp_mode
    public :: ocp_displace, ocp_displace_tangent, ocp_radial_histogram
    public :: ocp_debye_green, ocp_debye_shell

    real(dp), parameter :: pi = 3.14159265358979323846_dp, two_pi = 2*pi

    type :: ocp_model_t
        !> Lambda = n lambda_D^3 and box edge L/lambda_D.
        real(dp) :: debye_volume = 100.0_dp, box = 4.0_dp
        !> Requested rms PME force error relative to |c| e^2/a^2.
        real(dp) :: tolerance = 1.0e-4_dp
        !> Impurity nonlinear radius r_nl/lambda_D (0: no impurity).
        real(dp) :: impurity = 0.0_dp
        !> Outer step (omega_p^-1; <= 0 selects 0.01 (100/Lambda)^(1/3), which
        !> keeps the switch radius r2 ~ 0.4 a) and outer/inner encounter step
        !> factors; reference relative speed sqrt(2) v_th for the switch radii.
        real(dp) :: dt = 0.0_dp, eta_outer = 0.5_dp, eta_inner = 0.05_dp
        ! Resolved at ocp_system:
        integer :: electrons = 0, particles = 0
        real(dp) :: coulomb = 0.0_dp, impurity_charge = 0.0_dp
        type(plasma_scales_t) :: scales
        type(pme_parameters_t) :: pme
    end type
contains
    !> Plasma units: hbar is irrelevant (classical) and set to 1; the SI map
    !> takes the electron mass and charge as mass and charge units.
    function ocp_units(model) result(units)
        type(ocp_model_t), intent(in) :: model
        type(units_t) :: units
        units = reduced_units(1.0_dp, 1/(4*pi*model%debye_volume), electron_mass_si, &
                              elementary_charge_si)
    end function

    subroutine ocp_system(model, system, ok)
        type(ocp_model_t), intent(inout) :: model
        type(nbody_system_t), intent(out) :: system
        logical, intent(out) :: ok
        type(units_t) :: units
        type(species_t), allocatable :: particles(:)
        ok = .false.
        if (.not. (model%debye_volume > 0 .and. model%box > 0)) return
        if (.not. (model%impurity >= 0)) return
        units = ocp_units(model)
        model%coulomb = units%coulomb
        model%electrons = nint(model%debye_volume*model%box**3)
        model%particles = model%electrons
        if (model%impurity > 0) model%particles = model%electrons + 1
        model%scales = plasma_scales(units, electron(), real(model%electrons, dp)/ &
                                     model%box**3, 1.0_dp)
        ! The declared plasma units must agree with the shared definitions
        ! (up to the rounding of N to an integer).
        if (abs(model%scales%debye_length - 1) > 1/real(model%electrons, dp) .or. &
            abs(model%scales%plasma_frequency - 1) > 1/real(model%electrons, dp)) return
        allocate (particles(model%particles))
        particles = electron()
        model%impurity_charge = 0
        if (model%impurity > 0) then
            model%impurity_charge = -4*pi*model%debye_volume*model%impurity
            particles(model%particles)%charge = model%impurity_charge
        end if
        call pme_select(model%box, model%electrons, model%coulomb, model%tolerance, &
                        model%pme, ok)
        if (.not. ok) return
        call nbody_periodic_mesh(system, units, particles, model%pme, ok)
        if (.not. ok) return
        if (model%impurity > 0) call nbody_fix(system, model%particles, ok)
        if (.not. ok) return
        if (.not. (model%dt > 0)) model%dt = 0.01_dp*(100/model%debye_volume)**(1.0_dp/3)
        call nbody_set_encounter(system, model%dt, sqrt(2.0_dp), ok, &
                                 model%eta_outer, model%eta_inner)
    end subroutine

    !> Prepared ensemble member: electrons uniform in the cube (with the
    !> impurity, accepted with the bare-impurity Boltzmann factor
    !> exp(-kappa/r), r the minimum-image distance; repulsive, so no bound
    !> states), Maxwellian velocities at T = 1 with the mean velocity removed
    !> (zero total momentum; a uniform translation is not an oscillation).
    !> The impurity sits at the origin at rest. Returns rejected draws.
    subroutine ocp_prepare(model, seed, state, rejected, status)
        type(ocp_model_t), intent(in) :: model
        integer(int64), intent(in) :: seed
        real(dp), intent(out) :: state(:)
        integer, intent(out) :: rejected
        type(fortnum_status_t), intent(out) :: status
        type(rng_t) :: g
        real(dp) :: u(3), r, accept, mean(3)
        integer :: i, n, n3, a
        n = model%particles; n3 = 3*n; rejected = 0
        call status_set(status, FORTNUM_DOMAIN_ERROR, 'state size must be 6 N')
        if (size(state) /= 6*n) return
        call rng_seed(g, seed, status)
        if (status%code /= FORTNUM_OK) return
        state = 0
        do i = 1, model%electrons
            do
                do a = 1, 3
                    call rng_uniform(g, u(a))
                end do
                state(3*i - 2:3*i) = model%box*u
                if (model%impurity <= 0) exit
                u = state(3*i - 2:3*i) - model%box*anint(state(3*i - 2:3*i)/model%box)
                r = sqrt(sum(u**2))
                call rng_uniform(g, accept)
                if (r > 0) then
                    if (accept < exp(-model%impurity/r)) exit
                end if
                rejected = rejected + 1
            end do
            do a = 1, 3
                call rng_normal(g, state(n3 + 3*i - 3 + a))
            end do
        end do
        do a = 1, 3
            mean(a) = sum(state(n3 + a:n3 + 3*model%electrons:3))/model%electrons
            state(n3 + a:n3 + 3*model%electrons:3) = &
                state(n3 + a:n3 + 3*model%electrons:3) - mean(a)
        end do
        call status_set(status, FORTNUM_OK, '')
    end subroutine

    !> Kinetic temperature sum m v^2/(3 N - 3) of the electrons (m = 1).
    real(dp) function ocp_kinetic_temperature(model, state) result(t)
        type(ocp_model_t), intent(in) :: model
        real(dp), intent(in) :: state(:)
        integer :: n3
        n3 = 3*model%particles
        t = sum(state(n3 + 1:n3 + 3*model%electrons)**2)/(3*model%electrons - 3)
    end function

    !> Isokinetic velocity rescaling to T = 1. This is a different (non-
    !> Hamiltonian) dynamics, used only in a declared equilibration stage.
    subroutine ocp_thermostat(model, state)
        type(ocp_model_t), intent(in) :: model
        real(dp), intent(inout) :: state(:)
        integer :: n3
        n3 = 3*model%particles
        state(n3 + 1:n3 + 3*model%electrons) = state(n3 + 1:n3 + 3*model%electrons)/ &
            sqrt(ocp_kinetic_temperature(model, state))
    end subroutine

    !> One step dt of the Hamiltonian dynamics: FortNum Strang composition
    !> of far kicks and the encounter-resolved near-pair flow
    !> (kin6d_near_encounter); `clusters` reports the near clusters.
    subroutine ocp_advance(system, state, dt, status)
        type(nbody_system_t), intent(inout) :: system
        real(dp), intent(inout) :: state(:)
        real(dp), intent(in) :: dt
        type(fortnum_status_t), intent(out) :: status
        call nbody_step(system, state, dt, status)
    end subroutine

    !> Step with its tangent; cluster structure and substep counts are held
    !> at the primal values (piecewise-smooth derivative).
    subroutine ocp_advance_jvp(system, state, tangent, dt, status)
        type(nbody_system_t), intent(inout) :: system
        real(dp), intent(inout) :: state(:), tangent(:)
        real(dp), intent(in) :: dt
        type(fortnum_status_t), intent(out) :: status
        call nbody_step_jvp(system, state, tangent, dt, status)
    end subroutine

    !> Electron density mode rho_k = sum_i exp(-i k.x_i).
    complex(dp) function ocp_mode(model, state, k) result(rho)
        type(ocp_model_t), intent(in) :: model
        real(dp), intent(in) :: state(:), k(3)
        integer :: i
        real(dp) :: phase
        rho = (0.0_dp, 0.0_dp)
        do i = 1, model%electrons
            phase = k(1)*state(3*i - 2) + k(2)*state(3*i - 1) + k(3)*state(3*i)
            rho = rho + cmplx(cos(phase), -sin(phase), dp)
        end do
    end function

    !> Spatial-mode displacement x_a -> x_a + (eps/k) sin(k x_a), k = 2 pi/L,
    !> on the selected axes: density 1 - eps cos(k x_a) to first order with
    !> unchanged Maxwellian velocities.
    subroutine ocp_displace(model, state, eps, axes)
        type(ocp_model_t), intent(in) :: model
        real(dp), intent(inout) :: state(:)
        real(dp), intent(in) :: eps
        logical, intent(in) :: axes(3)
        real(dp) :: k
        integer :: i, a
        k = two_pi/model%box
        do i = 1, model%electrons
            do a = 1, 3
                if (axes(a)) state(3*i - 3 + a) = state(3*i - 3 + a) + &
                    eps/k*sin(k*state(3*i - 3 + a))
            end do
        end do
    end subroutine

    !> Tangent of ocp_displace with respect to eps at the undisplaced state.
    subroutine ocp_displace_tangent(model, state, axes, tangent)
        type(ocp_model_t), intent(in) :: model
        real(dp), intent(in) :: state(:)
        logical, intent(in) :: axes(3)
        real(dp), intent(out) :: tangent(:)
        real(dp) :: k
        integer :: i, a
        k = two_pi/model%box
        tangent = 0
        do i = 1, model%electrons
            do a = 1, 3
                if (axes(a)) tangent(3*i - 3 + a) = sin(k*state(3*i - 3 + a))/k
            end do
        end do
    end subroutine

    !> Add electron counts in radial bins of width rmax/size(counts) around
    !> the minimum image of `center`.
    subroutine ocp_radial_histogram(model, state, center, rmax, counts)
        type(ocp_model_t), intent(in) :: model
        real(dp), intent(in) :: state(:), center(3), rmax
        real(dp), intent(inout) :: counts(:)
        real(dp) :: d(3), r
        integer :: i, b
        do i = 1, model%electrons
            d = state(3*i - 2:3*i) - center
            d = d - model%box*anint(d/model%box)
            r = sqrt(sum(d**2))
            if (r >= rmax) cycle
            b = 1 + int(r/rmax*size(counts))
            counts(b) = counts(b) + 1
        end do
    end subroutine

    !> Periodic screened-Poisson Green function of unit coupling in the
    !> neutral cube, G(d) = sum_n exp(-kappa |d + n L|)/|d + n L| - 4 pi/(V kappa^2)
    !> (k /= 0 modes of the Debye-Hueckel linear response; leaves from
    !> gen/gen_ocp.f90). Image shells by max-norm until the shell bound is
    !> below 1e-17 of the sum; `bound` is the omitted-shell bound.
    real(dp) function ocp_debye_green(d, kappa, length, bound) result(g)
        real(dp), intent(in) :: d(3), kappa, length
        real(dp), intent(out) :: bound
        real(dp) :: w(3), u, b
        integer :: m, a, c, e, top
        w = d - length*anint(d/length)
        top = 1
        do m = 1, 200
            call debye_shell_bound(real(m, dp), kappa, length, b)
            if (b < 1.0e-17_dp) exit
            top = m
        end do
        g = 0
        do e = -top, top
            do c = -top, top
                do a = -top, top
                    call debye_image(w(1) + a*length, w(2) + c*length, w(3) + e*length, &
                                     1.0_dp, kappa, u)
                    g = g + u
                end do
            end do
        end do
        call debye_zero_mode(1.0_dp, kappa, length**3, u)
        g = g - u
        bound = 0
        do m = top + 1, top + 50
            call debye_shell_bound(real(m, dp), kappa, length, b)
            bound = bound + b
        end do
    end function

    !> Spherical average of ocp_debye_green at radius r around the origin
    !> (r <= L/2) over `directions` Fibonacci-sphere directions.
    real(dp) function ocp_debye_shell(r, kappa, length, directions) result(g)
        real(dp), intent(in) :: r, kappa, length
        integer, intent(in) :: directions
        real(dp) :: z, phi, rho, b
        integer :: j
        g = 0
        do j = 1, directions
            z = 1 - (2*j - 1.0_dp)/directions
            rho = sqrt(1 - z*z); phi = pi*(3 - sqrt(5.0_dp))*j
            g = g + ocp_debye_green(r*[rho*cos(phi), rho*sin(phi), z], kappa, length, b)
        end do
        g = g/directions
    end function
end module kin6d_ocp
