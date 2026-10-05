!> Classical point-charge N-body Hamiltonian H = sum p^2/(2m) + U shared by
!> the microscopic cases (C0 relative encounter, C1 periodic pair, C2 cloud).
!>
!> State layout: state(1:3n) positions, state(3n+1:6n) velocities.
!> Interactions come only from kin6d_coulomb_leaves (one Coulomb definition):
!> - isolated: bare pair Coulomb, optional fixed point sources and a fixed
!>   uniform sphere centred at the origin;
!> - periodic: Ewald sum (kin6d_ewald parameters) with neutralizing
!>   background; the energy includes self-image, self and background terms.
!> - periodic mesh: particle-mesh Ewald (kin6d_pme: cell-list real space and
!>   B-spline mesh), the scalable action verified against the direct sum.
!> Particles may be fixed (infinite mass: an impurity held at rest).
!> Pair work is a tiled direct sum (or structure-factor reciprocal sum); the
!> logical all-pairs matrix is never stored. There is no softening: a
!> coincident pair or non-finite force fails closed with a domain error.
!> One splitting path: FortNum geometric Strang (kick-drift-kick).
module kin6d_nbody
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortnum_status, only: fortnum_status_t, FORTNUM_OK, FORTNUM_DOMAIN_ERROR, &
                              status_set
    use fortnum_ode_geometric, only: geometric_strang_step
    use kin6d_coulomb_leaves, only: coulomb_pair, coulomb_pair_jvp, &
        coulomb_sphere_inside, coulomb_sphere_inside_jvp, ewald_real_pair, &
        ewald_real_pair_jvp, ewald_reciprocal_weight, ewald_self, ewald_background
    use kin6d_ewald, only: ewald_parameters_t
    use kin6d_pme, only: pme_parameters_t, pme_t, pme_setup, pme_force, &
                         pme_force_jvp, pme_constant
    use kin6d_species, only: units_t, species_t
    implicit none
    private
    public :: nbody_system_t, nbody_isolated, nbody_periodic, nbody_add_source
    public :: nbody_set_sphere, nbody_force, nbody_force_jvp, nbody_energy
    public :: nbody_step, nbody_step_jvp, nbody_encounter_step, nbody_wrap
    public :: nbody_periodic_mesh, nbody_fix, nbody_source_encounter_step

    integer, parameter :: tile = 64
    real(dp), parameter :: two_pi = 6.2831853071795864769_dp

    type :: nbody_system_t
        integer :: n = 0
        real(dp) :: coulomb = 1.0_dp
        real(dp), allocatable :: charge(:), mass(:)
        !> Fixed point sources (isolated systems only).
        real(dp), allocatable :: source_position(:, :), source_charge(:)
        !> Fixed uniform sphere centred at the origin (isolated systems only).
        real(dp) :: sphere_charge = 0.0_dp, sphere_radius = 0.0_dp
        logical :: periodic = .false., mesh = .false.
        type(ewald_parameters_t) :: ewald
        type(pme_t) :: pme
        !> Fixed particles: kicks keep their velocity zero.
        logical, allocatable :: fixed(:)
        integer :: modes = 0
        integer, allocatable :: wave(:, :)
        !> Half-space reciprocal weights, doubled for the +/-k pair, g=1.
        real(dp), allocatable :: weight(:)
        ! Scratch sized once at setup; reused by every force evaluation.
        real(dp), allocatable :: force_cache(:), position_cache(:)
        real(dp), allocatable :: augmented(:), augmented_force(:)
        complex(dp), allocatable :: phase(:, :, :), mode_phase(:)
        !> Periodic self-image, self and background energy (set at setup).
        real(dp) :: constant_energy = 0.0_dp
        real(dp) :: potential_cache = 0.0_dp
        logical :: cached = .false.
    end type

    !> Context handed to the FortNum flows; scratch is reached by pointer.
    type :: flow_context_t
        type(nbody_system_t), pointer :: system => null()
    end type
contains
    subroutine nbody_isolated(system, units, particles, ok)
        type(nbody_system_t), intent(out) :: system
        type(units_t), intent(in) :: units
        type(species_t), intent(in) :: particles(:)
        logical, intent(out) :: ok
        call setup(system, units, particles, ok)
    end subroutine

    subroutine nbody_periodic(system, units, particles, ewald, ok)
        type(nbody_system_t), intent(out) :: system
        type(units_t), intent(in) :: units
        type(species_t), intent(in) :: particles(:)
        type(ewald_parameters_t), intent(in) :: ewald
        logical, intent(out) :: ok
        integer :: i, j, k, m, nk
        real(dp) :: volume, k2
        call setup(system, units, particles, ok)
        if (.not. ok) return
        ok = .false.
        if (ewald%length <= 0 .or. ewald%alpha <= 0) return
        if (ewald%real_cutoff < 0 .or. ewald%reciprocal_cutoff < 1) return
        system%periodic = .true.; system%ewald = ewald
        nk = ewald%reciprocal_cutoff
        m = ((2*nk + 1)**3 - 1)/2
        allocate (system%wave(3, m), system%weight(m))
        allocate (system%phase(0:nk, 3, system%n), system%mode_phase(system%n))
        volume = ewald%length**3
        m = 0
        do i = 0, nk
            do j = -nk, nk
                do k = -nk, nk
                    ! One representative of each +/-n pair; n = 0 omitted.
                    if (i == 0 .and. (j < 0 .or. (j == 0 .and. k <= 0))) cycle
                    m = m + 1
                    system%wave(:, m) = [i, j, k]
                    k2 = (two_pi/ewald%length)**2*real(i*i + j*j + k*k, dp)
                    call ewald_reciprocal_weight(k2, ewald%alpha, volume, 1.0_dp, &
                                                 system%weight(m))
                    system%weight(m) = 2*system%weight(m)
                end do
            end do
        end do
        system%modes = m
        call periodic_constant(system)
        ok = ieee_is_finite(system%constant_energy)
    end subroutine

    !> Periodic system with the particle-mesh Ewald action of parameters p
    !> (kin6d_pme); ewald%length carries the box edge for minimum images.
    subroutine nbody_periodic_mesh(system, units, particles, p, ok)
        type(nbody_system_t), intent(out) :: system
        type(units_t), intent(in) :: units
        type(species_t), intent(in) :: particles(:)
        type(pme_parameters_t), intent(in) :: p
        logical, intent(out) :: ok
        call setup(system, units, particles, ok)
        if (.not. ok) return
        call pme_setup(system%pme, p, system%n, ok)
        if (.not. ok) return
        system%periodic = .true.; system%mesh = .true.
        system%ewald%length = p%length; system%ewald%alpha = p%alpha
        system%constant_energy = pme_constant(system%pme, system%coulomb, system%charge)
        ok = ieee_is_finite(system%constant_energy)
    end subroutine

    !> Hold particle i at rest (infinite mass); its velocity must be zero.
    subroutine nbody_fix(system, i, ok)
        type(nbody_system_t), intent(inout) :: system
        integer, intent(in) :: i
        logical, intent(out) :: ok
        ok = i >= 1 .and. i <= system%n
        if (ok) system%fixed(i) = .true.
    end subroutine

    subroutine setup(system, units, particles, ok)
        type(nbody_system_t), intent(out) :: system
        type(units_t), intent(in) :: units
        type(species_t), intent(in) :: particles(:)
        logical, intent(out) :: ok
        ok = .false.
        if (size(particles) < 1) return
        if (.not. ieee_is_finite(units%coulomb)) return
        if (.not. all(ieee_is_finite(particles%charge))) return
        if (.not. all(ieee_is_finite(particles%mass))) return
        if (any(particles%mass <= 0)) return
        system%n = size(particles); system%coulomb = units%coulomb
        system%charge = particles%charge; system%mass = particles%mass
        allocate (system%source_position(3, 0), system%source_charge(0))
        allocate (system%force_cache(3*system%n), system%position_cache(3*system%n))
        allocate (system%fixed(system%n))
        system%fixed = .false.
        ok = .true.
    end subroutine

    !> Add a fixed point charge (isolated systems). Setup only, not hot.
    subroutine nbody_add_source(system, position, charge, ok)
        type(nbody_system_t), intent(inout) :: system
        real(dp), intent(in) :: position(3), charge
        logical, intent(out) :: ok
        ok = .not. system%periodic .and. all(ieee_is_finite(position)) .and. &
             ieee_is_finite(charge)
        if (.not. ok) return
        system%source_position = reshape([reshape(system%source_position, &
            [size(system%source_position)]), position], &
            [3, size(system%source_charge) + 1])
        system%source_charge = [system%source_charge, charge]
        system%cached = .false.
    end subroutine

    subroutine nbody_set_sphere(system, charge, radius, ok)
        type(nbody_system_t), intent(inout) :: system
        real(dp), intent(in) :: charge, radius
        logical, intent(out) :: ok
        ok = .not. system%periodic .and. ieee_is_finite(charge) .and. &
             ieee_is_finite(radius) .and. radius > 0
        if (.not. ok) return
        system%sphere_charge = charge; system%sphere_radius = radius
        system%cached = .false.
    end subroutine

    !> Minimum-image wrap of a separation into [-L/2, L/2].
    elemental real(dp) function nbody_wrap(d, length) result(w)
        real(dp), intent(in) :: d, length
        w = d - length*anint(d/length)
    end function

    !> Forces F = -grad U on every particle and the potential energy U.
    subroutine nbody_force(system, state, force, potential, status)
        type(nbody_system_t), intent(inout) :: system
        real(dp), intent(in) :: state(:)
        real(dp), intent(out) :: force(:), potential
        type(fortnum_status_t), intent(out) :: status
        logical :: ok
        integer :: n3
        call status_set(status, FORTNUM_DOMAIN_ERROR, 'invalid N-body state')
        force = 0; potential = 0
        n3 = 3*system%n
        if (system%n < 1 .or. size(state) < n3 .or. size(force) /= n3) return
        if (.not. all(ieee_is_finite(state(:n3)))) return
        if (system%periodic) then
            call periodic_force(system, state(:n3), force, potential, ok)
        else
            call isolated_force(system, state(:n3), force, potential, ok)
        end if
        if (ok) ok = all(ieee_is_finite(force)) .and. ieee_is_finite(potential)
        if (.not. ok) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, &
                'coincident charges or non-finite force; no softening applied')
            return
        end if
        call status_set(status, FORTNUM_OK, '')
    end subroutine

    subroutine isolated_force(system, x, force, potential, ok)
        type(nbody_system_t), intent(in) :: system
        real(dp), intent(in) :: x(:)
        real(dp), intent(inout) :: force(:), potential
        logical, intent(out) :: ok
        integer :: ib, ie, jb, je, i, j
        ok = .true.
        do ib = 1, system%n, tile
            ie = min(system%n, ib + tile - 1)
            do i = ib, ie
                call external_force(system, x, i, force, potential, ok)
                if (.not. ok) return
                do j = i + 1, ie
                    call bare_pair(system, x, i, j, force, potential, ok)
                    if (.not. ok) return
                end do
            end do
            do jb = ie + 1, system%n, tile
                je = min(system%n, jb + tile - 1)
                do i = ib, ie
                    do j = jb, je
                        call bare_pair(system, x, i, j, force, potential, ok)
                        if (.not. ok) return
                    end do
                end do
            end do
        end do
    end subroutine

    subroutine bare_pair(system, x, i, j, force, potential, ok)
        type(nbody_system_t), intent(in) :: system
        real(dp), intent(in) :: x(:)
        integer, intent(in) :: i, j
        real(dp), intent(inout) :: force(:), potential
        logical, intent(out) :: ok
        real(dp) :: d(3), f(3), u
        d = x(3*i - 2:3*i) - x(3*j - 2:3*j)
        ok = maxval(abs(d)) > 0
        if (.not. ok) return
        call coulomb_pair(d(1), d(2), d(3), &
                          system%coulomb*system%charge(i)*system%charge(j), &
                          f(1), f(2), f(3), u)
        force(3*i - 2:3*i) = force(3*i - 2:3*i) + f
        force(3*j - 2:3*j) = force(3*j - 2:3*j) - f
        potential = potential + u
    end subroutine

    subroutine external_force(system, x, i, force, potential, ok)
        type(nbody_system_t), intent(in) :: system
        real(dp), intent(in) :: x(:)
        integer, intent(in) :: i
        real(dp), intent(inout) :: force(:), potential
        logical, intent(out) :: ok
        real(dp) :: p(3), d(3), f(3), u, g
        integer :: s
        ok = .true.
        p = x(3*i - 2:3*i)
        if (system%sphere_radius > 0) then
            g = system%coulomb*system%charge(i)*system%sphere_charge
            if (sqrt(p(1)*p(1) + p(2)*p(2) + p(3)*p(3)) <= system%sphere_radius) then
                call coulomb_sphere_inside(p(1), p(2), p(3), g, system%sphere_radius, &
                                           f(1), f(2), f(3), u)
            else
                call coulomb_pair(p(1), p(2), p(3), g, f(1), f(2), f(3), u)
            end if
            force(3*i - 2:3*i) = force(3*i - 2:3*i) + f
            potential = potential + u
        end if
        do s = 1, size(system%source_charge)
            d = p - system%source_position(:, s)
            ok = maxval(abs(d)) > 0
            if (.not. ok) return
            g = system%coulomb*system%charge(i)*system%source_charge(s)
            call coulomb_pair(d(1), d(2), d(3), g, f(1), f(2), f(3), u)
            force(3*i - 2:3*i) = force(3*i - 2:3*i) + f
            potential = potential + u
        end do
    end subroutine

    subroutine periodic_force(system, x, force, potential, ok)
        type(nbody_system_t), intent(inout) :: system
        real(dp), intent(in) :: x(:)
        real(dp), intent(inout) :: force(:), potential
        logical, intent(out) :: ok
        integer :: ib, ie, jb, je, i, j, m
        real(dp) :: k(3), t, q
        complex(dp) :: s
        ok = .true.
        if (system%mesh) then
            call pme_force(system%pme, system%coulomb, system%charge, x, force, &
                           potential, ok)
            potential = potential + system%constant_energy
            return
        end if
        do ib = 1, system%n, tile
            ie = min(system%n, ib + tile - 1)
            do jb = ib, system%n, tile
                je = min(system%n, jb + tile - 1)
                do i = ib, ie
                    do j = max(jb, i + 1), je
                        call ewald_pair(system, x, i, j, force, potential, ok)
                        if (.not. ok) return
                    end do
                end do
            end do
        end do
        potential = potential + system%constant_energy
        call fill_phases(system, x)
        do m = 1, system%modes
            k = two_pi/system%ewald%length*real(system%wave(:, m), dp)
            s = (0.0_dp, 0.0_dp)
            do i = 1, system%n
                system%mode_phase(i) = mode_phase(system, m, i)
                s = s + system%charge(i)*system%mode_phase(i)
            end do
            potential = potential + system%coulomb*system%weight(m)* &
                        (real(s)**2 + aimag(s)**2)/2
            do i = 1, system%n
                q = system%coulomb*system%charge(i)*system%weight(m)
                t = aimag(conjg(system%mode_phase(i))*s)
                force(3*i - 2:3*i) = force(3*i - 2:3*i) - q*t*k
            end do
        end do
    end subroutine

    !> Real-space Ewald interaction of particles i<j over the image cube.
    subroutine ewald_pair(system, x, i, j, force, potential, ok)
        type(nbody_system_t), intent(in) :: system
        real(dp), intent(in) :: x(:)
        integer, intent(in) :: i, j
        real(dp), intent(inout) :: force(:), potential
        logical, intent(out) :: ok
        real(dp) :: d(3), e(3), f(3), u, g, length
        integer :: a, b, c, nr
        length = system%ewald%length; nr = system%ewald%real_cutoff
        g = system%coulomb*system%charge(i)*system%charge(j)
        d = nbody_wrap(x(3*i - 2:3*i) - x(3*j - 2:3*j), length)
        ok = .true.
        do a = -nr, nr
            do b = -nr, nr
                do c = -nr, nr
                    e(1) = d(1) + length*a; e(2) = d(2) + length*b; e(3) = d(3) + length*c
                    ok = maxval(abs(e)) > 0
                    if (.not. ok) return
                    call ewald_real_pair(e(1), e(2), e(3), system%ewald%alpha, g, &
                                         f(1), f(2), f(3), u)
                    force(3*i - 2:3*i) = force(3*i - 2:3*i) + f
                    force(3*j - 2:3*j) = force(3*j - 2:3*j) - f
                    potential = potential + u
                end do
            end do
        end do
    end subroutine

    !> Configuration-independent periodic energy: each charge with its own
    !> images (real part), the Gaussian self term and the neutralizing
    !> background. The reciprocal self part is inside |S(k)|^2.
    subroutine periodic_constant(system)
        type(nbody_system_t), intent(inout) :: system
        real(dp) :: g, u, f(3), length, volume
        integer :: i, a, b, c, nr
        length = system%ewald%length; nr = system%ewald%real_cutoff
        volume = length**3
        system%constant_energy = 0
        do i = 1, system%n
            g = system%coulomb*system%charge(i)**2
            do a = -nr, nr
                do b = -nr, nr
                    do c = -nr, nr
                        if (a == 0 .and. b == 0 .and. c == 0) cycle
                        call ewald_real_pair(length*a, length*b, length*c, &
                                             system%ewald%alpha, g, f(1), f(2), f(3), u)
                        system%constant_energy = system%constant_energy + u/2
                    end do
                end do
            end do
            call ewald_self(system%ewald%alpha, g, u)
            system%constant_energy = system%constant_energy + u
        end do
        call ewald_background(system%ewald%alpha, volume, &
                              system%coulomb*sum(system%charge)**2, u)
        system%constant_energy = system%constant_energy + u/2
    end subroutine

    subroutine fill_phases(system, x)
        type(nbody_system_t), intent(inout) :: system
        real(dp), intent(in) :: x(:)
        integer :: i, a, m
        real(dp) :: theta
        do i = 1, system%n
            do a = 1, 3
                theta = two_pi*x(3*(i - 1) + a)/system%ewald%length
                system%phase(0, a, i) = (1.0_dp, 0.0_dp)
                do m = 1, system%ewald%reciprocal_cutoff
                    system%phase(m, a, i) = cmplx(cos(m*theta), sin(m*theta), dp)
                end do
            end do
        end do
    end subroutine

    !> exp(i k.x_i) for half-space mode m from the per-axis phase table.
    pure complex(dp) function mode_phase(system, m, i) result(e)
        type(nbody_system_t), intent(in) :: system
        integer, intent(in) :: m, i
        integer :: a, w
        e = (1.0_dp, 0.0_dp)
        do a = 1, 3
            w = system%wave(a, m)
            if (w >= 0) then
                e = e*system%phase(w, a, i)
            else
                e = e*conjg(system%phase(-w, a, i))
            end if
        end do
    end function

    !> Total energy sum m v^2/2 + U.
    subroutine nbody_energy(system, state, energy, status)
        type(nbody_system_t), intent(inout) :: system
        real(dp), intent(in) :: state(:)
        real(dp), intent(out) :: energy
        type(fortnum_status_t), intent(out) :: status
        integer :: i, j
        call cached_force(system, state, status)
        energy = system%potential_cache
        if (status%code /= FORTNUM_OK) return
        if (size(state) < 6*system%n) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, 'state lacks velocities')
            return
        end if
        do i = 1, system%n
            j = 3*system%n + 3*i - 2
            energy = energy + system%mass(i)*(state(j)*state(j) + &
                     state(j + 1)*state(j + 1) + state(j + 2)*state(j + 2))/2
        end do
    end subroutine

    !> Force at the current positions, reused while positions are unchanged
    !> (the KDK closing kick and the next opening kick share one evaluation).
    subroutine cached_force(system, state, status)
        type(nbody_system_t), intent(inout) :: system
        real(dp), intent(in) :: state(:)
        type(fortnum_status_t), intent(out) :: status
        integer :: n3
        n3 = 3*system%n
        if (system%cached .and. size(state) >= n3) then
            if (maxval(abs(state(:n3) - system%position_cache)) <= 0) then
                call status_set(status, FORTNUM_OK, '')
                return
            end if
        end if
        system%cached = .false.
        call nbody_force(system, state, system%force_cache, system%potential_cache, &
                         status)
        if (status%code /= FORTNUM_OK) return
        system%position_cache = state(:n3)
        system%cached = .true.
    end subroutine

    !> One kick-drift-kick step of size dt (FortNum geometric Strang).
    subroutine nbody_step(system, state, dt, status)
        type(nbody_system_t), intent(inout), target :: system
        real(dp), intent(inout) :: state(:)
        real(dp), intent(in) :: dt
        type(fortnum_status_t), intent(out) :: status
        type(flow_context_t) :: context
        call status_set(status, FORTNUM_DOMAIN_ERROR, 'state size must be 6n')
        if (size(state) /= 6*system%n) return
        context%system => system
        call geometric_strang_step(state, dt, kick, drift, status, context)
    end subroutine

    !> Step the state and its tangent [dx, dv] together: the tangent of a
    !> Strang step is the Strang step of the tangent-augmented flows.
    subroutine nbody_step_jvp(system, state, tangent, dt, status)
        type(nbody_system_t), intent(inout), target :: system
        real(dp), intent(inout) :: state(:), tangent(:)
        real(dp), intent(in) :: dt
        type(fortnum_status_t), intent(out) :: status
        type(flow_context_t) :: context
        integer :: n6
        n6 = 6*system%n
        call status_set(status, FORTNUM_DOMAIN_ERROR, 'state/tangent size must be 6n')
        if (size(state) /= n6 .or. size(tangent) /= n6) return
        if (.not. allocated(system%augmented)) then
            allocate (system%augmented(2*n6), system%augmented_force(3*system%n))
        end if
        system%augmented(:n6) = state; system%augmented(n6 + 1:) = tangent
        context%system => system
        call geometric_strang_step(system%augmented, dt, kick, drift, status, context)
        if (status%code /= FORTNUM_OK) return
        state = system%augmented(:n6); tangent = system%augmented(n6 + 1:)
    end subroutine

    subroutine kick(state, dt, status, context)
        real(dp), intent(inout) :: state(:)
        real(dp), intent(in) :: dt
        type(fortnum_status_t), intent(out) :: status
        class(*), intent(in), optional :: context
        type(nbody_system_t), pointer :: s
        integer :: i, n3, n6
        call status_set(status, FORTNUM_DOMAIN_ERROR, 'N-body context required')
        if (.not. present(context)) return
        select type (context)
        type is (flow_context_t)
            s => context%system
            n3 = 3*s%n; n6 = 6*s%n
            call cached_force(s, state, status)
            if (status%code /= FORTNUM_OK) return
            do i = 1, s%n
                if (s%fixed(i)) then
                    state(n3 + 3*i - 2:n3 + 3*i) = 0; cycle
                end if
                state(n3 + 3*i - 2:n3 + 3*i) = state(n3 + 3*i - 2:n3 + 3*i) + &
                    dt*(s%force_cache(3*i - 2:3*i)/s%mass(i))
            end do
            if (size(state) == 2*n6) then
                ! Tangent: dv += dt dF(x; dx)/m, using the augmented scratch.
                call nbody_force_jvp(s, state(:n3), state(n6 + 1:n6 + n3), &
                                     s%augmented_force, status)
                if (status%code /= FORTNUM_OK) return
                do i = 1, s%n
                    if (s%fixed(i)) then
                        state(n6 + n3 + 3*i - 2:n6 + n3 + 3*i) = 0; cycle
                    end if
                    state(n6 + n3 + 3*i - 2:n6 + n3 + 3*i) = &
                        state(n6 + n3 + 3*i - 2:n6 + n3 + 3*i) + &
                        dt*(s%augmented_force(3*i - 2:3*i)/s%mass(i))
                end do
            end if
        end select
    end subroutine

    subroutine drift(state, dt, status, context)
        real(dp), intent(inout) :: state(:)
        real(dp), intent(in) :: dt
        type(fortnum_status_t), intent(out) :: status
        class(*), intent(in), optional :: context
        integer :: n3, n6
        call status_set(status, FORTNUM_DOMAIN_ERROR, 'N-body context required')
        if (.not. present(context)) return
        select type (context)
        type is (flow_context_t)
            n3 = 3*context%system%n; n6 = 6*context%system%n
            state(:n3) = state(:n3) + dt*state(n3 + 1:n6)
            if (size(state) == 2*n6) state(n6 + 1:n6 + n3) = state(n6 + 1:n6 + n3) + &
                dt*state(n6 + n3 + 1:)
            call status_set(status, FORTNUM_OK, '')
        end select
    end subroutine

    !> Directional derivative dF = (dF/dx) dx of the forces (Hessian action).
    subroutine nbody_force_jvp(system, position, dposition, dforce, status)
        type(nbody_system_t), intent(inout) :: system
        real(dp), intent(in) :: position(:), dposition(:)
        real(dp), intent(out) :: dforce(:)
        type(fortnum_status_t), intent(out) :: status
        integer :: i, j, n3
        logical :: ok
        call status_set(status, FORTNUM_DOMAIN_ERROR, 'invalid N-body tangent')
        dforce = 0
        n3 = 3*system%n
        if (size(position) /= n3 .or. size(dposition) /= n3 .or. &
            size(dforce) /= n3) return
        ok = .true.
        if (system%mesh) then
            call pme_force_jvp(system%pme, system%coulomb, system%charge, position, &
                               dposition, dforce, ok)
        else
            call direct_jvp(ok)
        end if
        if (.not. ok .or. .not. all(ieee_is_finite(dforce))) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, &
                'coincident charges or non-finite force tangent')
            return
        end if
        call status_set(status, FORTNUM_OK, '')
    contains
        subroutine direct_jvp(ok)
            logical, intent(out) :: ok
            ok = .true.
            do i = 1, system%n
                if (.not. system%periodic) call external_jvp(system, position, &
                                                             dposition, i, dforce, ok)
                if (.not. ok) exit
                do j = i + 1, system%n
                    call pair_jvp(system, position, dposition, i, j, dforce, ok)
                    if (.not. ok) exit
                end do
                if (.not. ok) exit
            end do
            if (ok .and. system%periodic) call reciprocal_jvp(system, position, &
                                                              dposition, dforce)
        end subroutine
    end subroutine

    subroutine pair_jvp(system, x, dx, i, j, dforce, ok)
        type(nbody_system_t), intent(in) :: system
        real(dp), intent(in) :: x(:), dx(:)
        integer, intent(in) :: i, j
        real(dp), intent(inout) :: dforce(:)
        logical, intent(out) :: ok
        real(dp) :: d(3), e(3), dd(3), df(3), du, g
        integer :: a, b, c, nr
        g = system%coulomb*system%charge(i)*system%charge(j)
        d = x(3*i - 2:3*i) - x(3*j - 2:3*j)
        dd = dx(3*i - 2:3*i) - dx(3*j - 2:3*j)
        ok = .true.
        if (.not. system%periodic) then
            ok = maxval(abs(d)) > 0
            if (.not. ok) return
            call coulomb_pair_jvp(d(1), d(2), d(3), g, dd(1), dd(2), dd(3), 0.0_dp, &
                                  df(1), df(2), df(3), du)
            dforce(3*i - 2:3*i) = dforce(3*i - 2:3*i) + df
            dforce(3*j - 2:3*j) = dforce(3*j - 2:3*j) - df
            return
        end if
        d = nbody_wrap(d, system%ewald%length); nr = system%ewald%real_cutoff
        do a = -nr, nr
            do b = -nr, nr
                do c = -nr, nr
                    e(1) = d(1) + system%ewald%length*a
                    e(2) = d(2) + system%ewald%length*b
                    e(3) = d(3) + system%ewald%length*c
                    ok = maxval(abs(e)) > 0
                    if (.not. ok) return
                    call ewald_real_pair_jvp(e(1), e(2), e(3), system%ewald%alpha, g, &
                                             dd(1), dd(2), dd(3), 0.0_dp, 0.0_dp, &
                                             df(1), df(2), df(3), du)
                    dforce(3*i - 2:3*i) = dforce(3*i - 2:3*i) + df
                    dforce(3*j - 2:3*j) = dforce(3*j - 2:3*j) - df
                end do
            end do
        end do
    end subroutine

    subroutine external_jvp(system, x, dx, i, dforce, ok)
        type(nbody_system_t), intent(in) :: system
        real(dp), intent(in) :: x(:), dx(:)
        integer, intent(in) :: i
        real(dp), intent(inout) :: dforce(:)
        logical, intent(out) :: ok
        real(dp) :: p(3), dp3(3), d(3), df(3), du, g
        integer :: s
        ok = .true.
        p = x(3*i - 2:3*i); dp3 = dx(3*i - 2:3*i)
        if (system%sphere_radius > 0) then
            g = system%coulomb*system%charge(i)*system%sphere_charge
            if (sqrt(p(1)*p(1) + p(2)*p(2) + p(3)*p(3)) <= system%sphere_radius) then
                call coulomb_sphere_inside_jvp(p(1), p(2), p(3), g, &
                    system%sphere_radius, dp3(1), dp3(2), dp3(3), 0.0_dp, 0.0_dp, &
                    df(1), df(2), df(3), du)
            else
                call coulomb_pair_jvp(p(1), p(2), p(3), g, dp3(1), dp3(2), dp3(3), &
                                      0.0_dp, df(1), df(2), df(3), du)
            end if
            dforce(3*i - 2:3*i) = dforce(3*i - 2:3*i) + df
        end if
        do s = 1, size(system%source_charge)
            d = p - system%source_position(:, s)
            ok = maxval(abs(d)) > 0
            if (.not. ok) return
            g = system%coulomb*system%charge(i)*system%source_charge(s)
            call coulomb_pair_jvp(d(1), d(2), d(3), g, dp3(1), dp3(2), dp3(3), &
                                  0.0_dp, df(1), df(2), df(3), du)
            dforce(3*i - 2:3*i) = dforce(3*i - 2:3*i) + df
        end do
    end subroutine

    !> Tangent of F_i = -k q_i sum_m W_m k_m Im(conj(e_i) S_m).
    subroutine reciprocal_jvp(system, x, dx, dforce)
        type(nbody_system_t), intent(inout) :: system
        real(dp), intent(in) :: x(:), dx(:)
        real(dp), intent(inout) :: dforce(:)
        integer :: i, m
        real(dp) :: k(3), q, kdx
        complex(dp) :: s, ds, z
        call fill_phases(system, x)
        do m = 1, system%modes
            k = two_pi/system%ewald%length*real(system%wave(:, m), dp)
            s = (0.0_dp, 0.0_dp); ds = (0.0_dp, 0.0_dp)
            do i = 1, system%n
                system%mode_phase(i) = mode_phase(system, m, i)
                kdx = dot_product(k, dx(3*i - 2:3*i))
                s = s + system%charge(i)*system%mode_phase(i)
                ds = ds + system%charge(i)*kdx*(0.0_dp, 1.0_dp)*system%mode_phase(i)
            end do
            do i = 1, system%n
                q = system%coulomb*system%charge(i)*system%weight(m)
                kdx = dot_product(k, dx(3*i - 2:3*i))
                z = conjg(system%mode_phase(i))*((0.0_dp, -1.0_dp)*kdx*s + ds)
                dforce(3*i - 2:3*i) = dforce(3*i - 2:3*i) - q*aimag(z)*k
            end do
        end do
    end subroutine

    !> Encounter-resolved step: dt = eta * min over interacting pairs (and
    !> fixed sources) of the Kepler time sqrt(r^3 mu/|g|) and the crossing
    !> time r/|v_rel|; a uniform sphere adds its harmonic time
    !> sqrt(R^3 m/|g|). Periodic pairs use the minimum-image separation.
    !> A coincident pair fails closed; nothing is softened.
    subroutine nbody_encounter_step(system, state, eta, dt, status)
        type(nbody_system_t), intent(in) :: system
        real(dp), intent(in) :: state(:), eta
        real(dp), intent(out) :: dt
        type(fortnum_status_t), intent(out) :: status
        integer :: i, j, s, n3
        real(dp) :: d(3), v(3), mu, g
        call status_set(status, FORTNUM_DOMAIN_ERROR, 'invalid encounter-step input')
        dt = huge(1.0_dp)
        n3 = 3*system%n
        if (size(state) /= 2*n3 .or. .not. (eta > 0)) return
        if (.not. all(ieee_is_finite(state))) return
        do i = 1, system%n
            do j = i + 1, system%n
                d = state(3*i - 2:3*i) - state(3*j - 2:3*j)
                if (system%periodic) d = nbody_wrap(d, system%ewald%length)
                v = state(n3 + 3*i - 2:n3 + 3*i) - state(n3 + 3*j - 2:n3 + 3*j)
                mu = system%mass(i)*system%mass(j)/(system%mass(i) + system%mass(j))
                g = system%coulomb*system%charge(i)*system%charge(j)
                call limit(d, v, mu, g, dt)
            end do
            if (system%periodic) cycle
            do s = 1, size(system%source_charge)
                d = state(3*i - 2:3*i) - system%source_position(:, s)
                v = state(n3 + 3*i - 2:n3 + 3*i)
                g = system%coulomb*system%charge(i)*system%source_charge(s)
                call limit(d, v, system%mass(i), g, dt)
            end do
            g = system%coulomb*system%charge(i)*system%sphere_charge
            if (system%sphere_radius > 0 .and. abs(g) > 0) dt = min(dt, &
                sqrt(system%sphere_radius**3*system%mass(i)/abs(g)))
        end do
        if (dt <= 0) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, &
                'coincident charges: encounter unresolved, no softening applied')
            return
        end if
        dt = eta*dt
        call status_set(status, FORTNUM_OK, '')
    contains
        pure subroutine limit(d, v, mu, g, dt)
            real(dp), intent(in) :: d(3), v(3), mu, g
            real(dp), intent(inout) :: dt
            real(dp) :: r, speed
            r = sqrt(d(1)*d(1) + d(2)*d(2) + d(3)*d(3))
            if (r <= 0) then
                dt = 0; return
            end if
            if (abs(g) > 0) dt = min(dt, sqrt(r**3*mu/abs(g)))
            speed = sqrt(v(1)*v(1) + v(2)*v(2) + v(3)*v(3))
            if (speed > 0) dt = min(dt, r/speed)
        end subroutine
    end subroutine

    !> O(N) encounter step against the fixed particles only (an impurity):
    !> dt = eta * min over mobile i and fixed s of the Kepler time
    !> sqrt(r^3 m_i/|g|) and the crossing time r/|v_i| (minimum image when
    !> periodic). Mobile-mobile encounters are not examined here.
    subroutine nbody_source_encounter_step(system, state, eta, dt, status)
        type(nbody_system_t), intent(in) :: system
        real(dp), intent(in) :: state(:), eta
        real(dp), intent(out) :: dt
        type(fortnum_status_t), intent(out) :: status
        integer :: i, s, n3
        real(dp) :: d(3), v(3), g, r, speed
        call status_set(status, FORTNUM_DOMAIN_ERROR, 'invalid encounter-step input')
        dt = huge(1.0_dp)
        n3 = 3*system%n
        if (size(state) /= 2*n3 .or. .not. (eta > 0)) return
        do s = 1, system%n
            if (.not. system%fixed(s)) cycle
            do i = 1, system%n
                if (system%fixed(i)) cycle
                d = state(3*i - 2:3*i) - state(3*s - 2:3*s)
                if (system%periodic) d = nbody_wrap(d, system%ewald%length)
                v = state(n3 + 3*i - 2:n3 + 3*i)
                r = sqrt(d(1)*d(1) + d(2)*d(2) + d(3)*d(3))
                if (.not. (r > 0)) then
                    call status_set(status, FORTNUM_DOMAIN_ERROR, &
                        'charge on a fixed source: no softening applied')
                    return
                end if
                g = system%coulomb*system%charge(i)*system%charge(s)
                if (abs(g) > 0) dt = min(dt, sqrt(r**3*system%mass(i)/abs(g)))
                speed = sqrt(v(1)*v(1) + v(2)*v(2) + v(3)*v(3))
                if (speed > 0) dt = min(dt, r/speed)
            end do
        end do
        dt = eta*dt
        call status_set(status, FORTNUM_OK, '')
    end subroutine
end module kin6d_nbody
