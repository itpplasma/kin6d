!> Quasilinear NTV torque density of every species on a set of flux surfaces
!! of a Boozer equilibrium with one toroidal harmonic n of a Boozer |B|
!! perturbation, from the monoenergetic kernel of milca_dk_ntv with the
!! NEO-2-QL switch set (weight (1 + xi^2)/b, Boozer streaming, both #192
!! sources) and pitch-angle scattering on all species.
!!
!! The perturbation file holds delta B = Re sum_m (bmnc - i bmns)
!! e^{i(m theta + n phi)} (NEO-2 convention), passed to the kernel as
!! c_{-m} = (bmnc + i bmns)/B_00 with its rotation argument -Omega_tE, where
!! Omega_tE is NEO-2's E x B rotation frequency. Per surface, species a and
!! speed v = v_a sqrt(x), v_a = sqrt(2 T_a/m_a), the kernel gives R(x), and
!!   K_1j = (T/(Z e))^2/(4 sqrt(pi) J_0) sum_k w_k x_k^(j-1) R(x_k),
!! with Gauss-Laguerre nodes and weights of int dx x^(5/2) e^-x and
!! J_0 = int dtheta (B_00/B)^2. The poloidal-flux particle flux is
!! Gamma_psi = -n (K_11 A_1 + K_12 A_2), A_1 = n'/n + (Z e/T) Omega_tE
!! - 3/2 A_2, A_2 = T'/T, ' = d/dpsi_pol with psi_pol = iota psi_tor per
!! radian, and the torque density on the species, positive along phi, is
!! Z e Gamma_psi. With r the effective radius of <|grad r|> = 1, NEO-2's
!! D_1j^NA = K_1j/(dpsi_pol/dr)^2.
!!
!! The deflection frequency of species a on all species b is
!!   nu_D(x) = sum_b nu_ab (erf(y_b) - G(y_b))/x^(3/2), y_b = v/v_b,
!!   nu_ab = (3 sqrt(pi)/8) v_a kappa_a n_b Z_b^2/(n_a Z_a^2),
!! with kappa_a = 2/(v_a tau_aa) the NEO-2 collisionality of the species.
!!
!! With orbit_width the monoenergetic R(x) comes from the finite-width
!! kernel of milca_ntv_orbit_geometry instead: the orbits of psi* = s of the
!! surface at the node's energy (canonical Maxwellian, exact guiding-centre
!! orbits with the E x B rotation Omega_tE of the profile in their
!! potential), the perturbation at the orbit point, pitch-angle scattering at
!! fixed psi*, on the hp mesh of the resolution. It has neither of the two
!! #192 sources: their sum is the first-order difference between the local
!! and the canonical drive, which this kernel replaces, so its thin limit is
!! the kernel with sources off.
!!
!! Work units (species, energy node, chunk of consecutive surfaces) run
!! under OpenMP with dynamic scheduling. Each unit factorises its first
!! surface and solves the following ones by GMRES preconditioned with that
!! factorisation, refactorising when the iteration count exceeds a bound.
module milca_ntv_device
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    !$ use omp_lib, only: omp_get_max_threads
    use fortnum_quadrature, only: gauss_gen_laguerre, gauss_legendre
    use milca_boozer_file, only: boozer_file_t, boozer_surface_t, boozer_surface, &
        boozer_average, BOOZER_OK
    use milca_ntv_boozer_line, only: ntv_boozer_line_t, ntv_boozer_line, NTV_LINE_OK
    use milca_dk_ntv, only: ntv_surface_t, ntv_kernel_t, ntv_factor_t, ntv_kernel, &
        ntv_f10, ntv_energy_coefficients, ntv_axisymmetric_coefficients, NTV_OK, &
            NTV_NO_CONVERGENCE, NTV_WEIGHT_FULL
    use milca_boozer_orbit, only: orbit_particle_t, orbit_set_potential, orbit_flux_p
    use milca_dk_hp, only: hp_operator_t, hp_mesh_t, hp_default_mesh, hp_ntv_kernel, HP_OK, &
        HP_SINGULAR, hp_ntv_estimate, hp_refine, hp_geometry_t
    use milca_ntv_hp_line, only: hp_line_field_t, hp_line_field
    use milca_ntv_orbit_geometry, only: hp_orbit_geometry_t, hp_orbit_geometry, &
        hp_orbit_perturbation, hp_orbit_confined
    use milca_dk_collop, only: collop_species_t, collop_grid_t, collop_make_grid, &
        collop_nuhat, collop_vth, COLLOP_OK
    implicit none
    private
    public :: ntv_species_t, ntv_profile_t, ntv_resolution_t, ntv_species_result_t
    public :: ntv_device_profile, deflection_frequency, ntv_read_input, &
        ntv_device_axisymmetric, ntv_axisymmetric_result_t, ntv_omega_te_from_rotation, &
        ntv_indicator_t, ntv_device_indicators
    integer, parameter, public :: DEVICE_OK = 0, DEVICE_BAD_INPUT = 1, DEVICE_SOLVE = 2, &
        DEVICE_NO_CONVERGENCE = 3, DEVICE_ORBIT_LOST = 4
    ! Mesh refinement modes of one kernel evaluation (adapt_kernel).
    integer, parameter :: MESH_FIXED = 0, MESH_SETTLE = 1, MESH_CHECK = 2

    !> End of one mesh adaptation: |eta|/|R + eta|, the unknowns and whether
    !> the tolerance was met before the size or round limit.
    type :: adapt_info_t
        real(dp) :: eta_rel = 0
        integer :: n_dof = 0
        logical :: settled = .true.
    end type

    real(dp), parameter :: pi = 4*atan(1.0_dp), elementary = 1.602176634e-19_dp

    !> Species constants: charge number and mass [kg].
    type :: ntv_species_t
        real(dp) :: z = 0, mass = 0
    end type

    !> Plasma on the surfaces: s, Omega_tE [rad/s], and per species (columns)
    !> density [m^-3], temperature [J], their s derivatives and the NEO-2
    !> collisionality kappa = 2/(v_th tau_aa) [1/m].
    type :: ntv_profile_t
        real(dp), allocatable :: s(:), omega_te(:)
        real(dp), allocatable :: density(:, :), temperature(:, :), density_s(:, :), &
            temperature_s(:, :), kappa(:, :)
    end type

    type :: ntv_resolution_t
        integer :: n_theta = 0, n_leg = 0, n_x = 0
        ! Surfaces per preconditioned sweep and the GMRES iteration count
        ! above which the next surface is factorised afresh.
        integer :: chunk = 8, max_iterations = 60
        logical :: sources = .true.
        ! 0: monoenergetic pitch-angle scattering on Gauss-Laguerre nodes;
        ! otherwise a milca_dk_collop model on the speed grid of nelem
        ! elements of the given order on [0, x_max] (energy-coupled kernel).
        integer :: collisions = 0, nelem = 4, order = 4
        real(dp) :: x_max = 5
        ! Monoenergetic energy integral: relative tolerance of adaptive
        ! composite Gauss-Legendre panels in x = v^2/v_th^2 on [0, 25]
        ! (0: n_x Gauss-Laguerre nodes). The E x B and precession resonance
        ! puts kinks into the integrand that fixed rules do not resolve: on
        ! AUG 30835 at T/10 the 24-node rule was off by a factor of 2.
        real(dp) :: x_tol = 1e-3_dp
        ! Finite orbit width (monoenergetic paths only) on the hp
        ! mesh of hp_default_mesh: element order, columns, bulk rows and the
        ! relative width of the layer at the separatrix.
        logical :: orbit_width = .false.
        ! Thin orbits: the spectral-element kernel on the Boozer line, on the
        ! mesh below (default), or the Legendre kernel of milca_dk_ntv.
        logical :: thin_legendre = .false.
        integer :: hp_order = 4, hp_columns = 24, hp_rows = 8
        real(dp) :: hp_delta = 1e-3_dp
        ! Goal-oriented refinement of that mesh (hp_ntv_estimate) at seed
        ! energies per surface and species, until the corrected kernel R + eta
        ! changes by at most hp_adapt relative between two refinements; every
        ! energy is then solved on the union of its two neighbouring seed
        ! meshes. 0 keeps the mesh.
        real(dp) :: hp_adapt = 0
        ! Refinement produces no mesh with more unknowns than this.
        integer :: hp_max_unknowns = 250000
    end type

    !> Per surface and species: K_1j and K_2j [m^2/s Wb^2/m^2], dpsi_pol/dr [Wb/m],
    !> Gamma_psi [Wb/(m^3 s)], the energy flux -n T (K_21 A_1 + K_22 A_2)
    !> [J Wb/(m^3 s)], the torque density [N m/m^3], the offset rotation
    !> [rad/s], and the factorisations and GMRES iterations spent on it. The
    !> offset rotation is the Omega_tE at which Gamma_psi of the species
    !> vanishes with the coefficients held at their values for the given
    !> Omega_tE; in the 1/nu regime, where they do not depend on Omega_tE, it is
    !> the rotation the torque relaxes to. status is DEVICE_OK, or the failure
    !> of this surface and species, whose other fields are then zero; the
    !> profile status is DEVICE_SOLVE when any of them failed. With finite
    !> orbit width on Gauss-Laguerre nodes, energies whose orbits leave the
    !> file (lost to the wall or across the axis) carry no NTV in this model and
    !> contribute zero; lost_weight is their share of the K_11 weight. With
    !> hp_adapt, mesh_error is the largest final |eta|/|R + eta| over the seed
    !> and check energies of the mesh adaptation, and unsettled counts those
    !> whose refinement stopped at its size or round limit before the
    !> tolerance was met; mesh_unknowns is the largest mesh solved, adapted or
    !> on a fixed band mesh.
    type :: ntv_species_result_t
        real(dp) :: k11 = 0, k12 = 0, k21 = 0, k22 = 0, dpsi_dr = 0, gamma_psi = 0, &
            heat_flux = 0, torque = 0, offset_rotation = 0, lost_weight = 0, mesh_error = 0
        integer :: factorisations = 0, iterations = 0, status = DEVICE_OK, unsettled = 0, &
            mesh_unknowns = 0
    end type

    !> Axisymmetric transport matrix of one surface: d(a, b, i, j) of
    !> ntv_axisymmetric_coefficients (responses i = particle, energy, parallel
    !> flow of species a to the drives j = A_1, A_2, parallel field of species
    !> b), b.grad theta/b (kappa), q, dpsi_pol/dr [Wb/m], dpsi_pol/ds [Wb], the
    !> covariant toroidal field G [T m] of the file and G_phi, the same along
    !> the toroidal direction of Omega_tE, <B^2> [T^2], the parallel flow <u_par B> of each
    !> species [T m/s], the bootstrap current <j_par B> [T A/m^2], and the
    !> toroidal rotation V_phi/R = omega + k_pol/R^2 [rad/s] of each species with
    !> omega = Omega_tE + p'/(Z e n) and k_pol = G_phi (<u_par B> - omega G_phi)
    !> /<B^2> [m^2/s], and the status. With D_3j = -2 kappa q G d(a, b, 3, j) [m^2/s T per unit A_j],
    !> <u_par B>_a = -sum_b (D_31 A_1 + D_32 A_2) of species b; momentum-conserving
    !> self-collisions give D_31 = G T/(Z e) for a single species: A_1 drives
    !> the rigid rotation alone.
    type :: ntv_axisymmetric_result_t
        real(dp), allocatable :: d(:, :, :, :), flow(:), omega(:), k_pol(:)
        real(dp) :: kappa = 0, q = 0, dpsi_dr = 0, psi_pol_s = 0, g_cov = 0, g_phi = 0, &
            b2 = 0, bootstrap = 0
        integer :: status = DEVICE_OK
    end type

    !> Validity indicators of the thin-orbit model on one surface, per species
    !> at the thermal speed v_a = sqrt(2 T_a/m_a): the inverse aspect ratio
    !> eps = (B_max - B_min)/(B_max + B_min), the gradient length
    !> L = min(n/|dn/dr|, T/|dT/dr|) [m], the radial excursion of a passing
    !> orbit with v_par = v_a, dr = m v_a |G|/(Z e sqrt(<B^2>) |dpsi_pol/dr|) [m]
    !> (rho_L B/B_theta in the circular limit), and the banana width
    !> 2 sqrt(2 eps) dr of a deeply trapped orbit. The thin-orbit expansion
    !> needs dr/L << 1; the excursion grows like the speed.
    type :: ntv_indicator_t
        real(dp) :: eps = 0
        real(dp), allocatable :: gradient_length(:), excursion(:), banana(:)
        integer :: status = DEVICE_OK
    end type

contains

    !> nu_D of species a at x = v^2/v_a^2 [1/s].
    pure real(dp) function deflection_frequency(species, density, temperature, kappa, &
            a, x) result(nu)
        type(ntv_species_t), intent(in) :: species(:)
        real(dp), intent(in) :: density(:), temperature(:), kappa(:), x
        integer, intent(in) :: a
        real(dp) :: va, y, g, nu_ab
        integer :: b
        va = sqrt(2*temperature(a)/species(a)%mass)
        nu = 0
        do b = 1, size(species)
            y = sqrt(x*temperature(a)*species(b)%mass/(temperature(b)*species(a)%mass))
            g = (erf(y) - 2*y*exp(-y*y)/sqrt(pi))/(2*y*y)
            nu_ab = 3*sqrt(pi)/8*va*kappa(a)*density(b)*species(b)%z**2 &
                /(density(a)*species(a)%z**2)
            nu = nu + nu_ab*(erf(y) - g)/x**1.5_dp
        end do
    end function deflection_frequency

    !> Orbits outside the model: they leave the files, or the orbit-class
    !> mesh cannot be fitted to them (potato and stagnation orbits near the
    !> axis at high energy, where the system is singular). They carry no NTV;
    !> lost_weight reports their share.
    pure logical function orbit_outside(st)
        integer, intent(in) :: st
        orbit_outside = st == DEVICE_ORBIT_LOST .or. st == HP_SINGULAR
    end function orbit_outside

    !> Whether the thin-orbit kernel is the Legendre one (on request).
    pure logical function thin_legendre(res)
        type(ntv_resolution_t), intent(in) :: res
        thin_legendre = res%thin_legendre
    end function thin_legendre

    subroutine ntv_device_profile(eq, pert, n_tor, species, prof, res, result, status)
        type(boozer_file_t), target, intent(in) :: eq, pert
        integer, intent(in) :: n_tor
        type(ntv_species_t), intent(in) :: species(:)
        type(ntv_profile_t), intent(in) :: prof
        type(ntv_resolution_t), intent(in) :: res
        type(ntv_species_result_t), allocatable, intent(out) :: result(:, :)
        integer, intent(out) :: status
        type(ntv_surface_t), allocatable :: sf(:)
        real(dp), allocatable :: x(:), w(:), r(:, :, :), j0(:), dpsi_dr(:), psi_pol_s(:), &
            b_ref(:), phi(:)
        integer, allocatable :: units(:, :), nfac(:, :, :), nits(:, :, :), ustat(:)
        logical, allocatable :: lost(:, :, :)
        ! Energies x of the resonance kinks per surface and species, from the
        ! panels the adaptive quadrature refined.
        integer, parameter :: max_kinks = 16
        real(dp), allocatable :: kink(:, :, :)
        integer, allocatable :: nkink(:, :)
        ! With hp_adapt, per surface and species, meshes adapted at seed
        ! energies, normalised to theta - theta_x and lambda/lambda_c; band k
        ! joins seeds k and k + 1, and every energy is solved on its band.
        integer, parameter :: max_bands = 16
        type(hp_mesh_t), allocatable :: band(:, :, :)
        real(dp), allocatable :: band_seed(:, :, :)
        integer, allocatable :: nband(:, :)
        integer :: ns, nsp, nchunk, nunit, i, a, k, c, u

        ns = size(prof%s)
        nsp = size(species)
        status = DEVICE_BAD_INPUT
        if (ns < 1 .or. nsp < 1 .or. res%n_x < 1 .or. res%chunk < 1) return
        if (res%orbit_width .and. res%collisions /= 0) return
        allocate (sf(ns), j0(ns), dpsi_dr(ns), psi_pol_s(ns), b_ref(ns), phi(ns))
        do i = 1, ns
            call surface(eq, pert, n_tor, prof%s(i), res, sf(i), j0(i), dpsi_dr(i), &
                psi_pol_s(i), status, b_ref=b_ref(i))
            if (status /= DEVICE_OK) return
        end do
        ! Potential of the profile, dPhi/dpsi_pol = Omega_tE along NEO-2's
        ! toroidal direction (-phi of the file), by the trapezoidal rule.
        phi(1) = 0
        do i = 2, ns
            phi(i) = phi(i - 1) + 0.5_dp*(prof%omega_te(i - 1) + prof%omega_te(i)) &
                *orbit_flux_p(eq, prof%s(i - 1), prof%s(i))
        end do
        allocate (band(max_bands, ns, nsp), band_seed(max_bands + 1, ns, nsp), nband(ns, nsp))
        nband = 0
        if (res%x_tol > 0) then
            allocate (result(ns, nsp), ustat(ns*nsp), kink(max_kinks, ns, nsp), &
                nkink(ns, nsp))
            nkink = 0
            ! Surfaces and species fill the threads when there are enough of
            ! them; otherwise they run in turn and the energies of each panel
            ! take the threads (panel_eval).
            !$omp parallel do schedule(dynamic, 1) default(shared) private(u) &
            !$omp if(ns*nsp >= omp_get_max_threads())
            do u = 1, ns*nsp
                call adaptive((u - 1)/nsp + 1, mod(u - 1, nsp) + 1, ustat(u))
            end do
            !$omp end parallel do
            ! The full operator puts speed-element edges at the kinks found here.
            if (res%collisions /= 0) call energy_profile(eq, n_tor, species, prof, res, sf, &
                psi_pol_s, dpsi_dr, result, status, kink, nkink)
            do u = 1, ns*nsp
                associate (r => result((u - 1)/nsp + 1, mod(u - 1, nsp) + 1))
                    if (ustat(u) /= DEVICE_OK) r = ntv_species_result_t(status=ustat(u))
                end associate
            end do
            if (any(ustat /= DEVICE_OK)) status = DEVICE_SOLVE
            if (res%collisions == 0 .and. all(ustat == DEVICE_OK)) status = DEVICE_OK
            return
        end if
        if (res%collisions /= 0) then
            call energy_profile(eq, n_tor, species, prof, res, sf, psi_pol_s, dpsi_dr, &
                result, status)
            return
        end if
        allocate (x(res%n_x), w(res%n_x), result(ns, nsp))
        call gauss_gen_laguerre(res%n_x, 2.5_dp, x, w)
        if (res%hp_adapt > 0 .and. .not. thin_legendre(res)) then
            do a = 1, nsp
                do i = 1, ns
                    ! With few nodes, the nodes themselves are the seeds.
                    if (res%n_x <= 4) then
                        call seed_bands(i, a, x)
                    else
                        call seed_bands(i, a)
                    end if
                end do
            end do
        end if

        ! Units: (species, node, first surface of the chunk).
        nchunk = (ns + res%chunk - 1)/res%chunk
        nunit = nsp*res%n_x*nchunk
        allocate (units(3, nunit), ustat(nunit), r(ns, res%n_x, nsp), &
            nfac(ns, res%n_x, nsp), nits(ns, res%n_x, nsp), lost(ns, res%n_x, nsp))
        lost = .false.
        u = 0
        do a = 1, nsp
            do k = 1, res%n_x
                do c = 1, nchunk
                    u = u + 1
                    units(:, u) = [a, k, (c - 1)*res%chunk + 1]
                end do
            end do
        end do
        r = 0
        nfac = 0
        nits = 0
        ! Units are independent energies; the kernels' own threading is then
        ! nested and serial, which uses the cores far better than one
        ! factorisation at a time.
        !$omp parallel do schedule(dynamic, 1) default(shared) private(u)
        do u = 1, nunit
            call sweep(units(1, u), units(2, u), units(3, u), &
                min(units(3, u) + res%chunk - 1, ns), ustat(u))
        end do
        !$omp end parallel do

        do a = 1, nsp
            do i = 1, ns
                call assemble(i, a)
            end do
        end do
        ! A failed unit voids its species on every surface of its chunk.
        do u = 1, nunit
            if (ustat(u) == DEVICE_OK) cycle
            do i = units(3, u), min(units(3, u) + res%chunk - 1, ns)
                result(i, units(1, u)) = ntv_species_result_t(status=ustat(u))
            end do
        end do
        status = merge(DEVICE_OK, DEVICE_SOLVE, all(ustat == DEVICE_OK))

    contains

        subroutine sweep(a, k, first, last, st)
            integer, intent(in) :: a, k, first, last
            integer, intent(out) :: st
            type(ntv_factor_t) :: factor
            type(ntv_kernel_t) :: kern
            real(dp), allocatable :: f10(:, :)
            real(dp) :: v, nu
            integer :: i, ks
            logical :: fresh
            st = DEVICE_SOLVE
            if (res%orbit_width) then
                do i = first, last
                    ! The node sum stops where the adaptive panels do: beyond
                    ! x = 25 lie 1e-8 of the K11 and 4e-7 of the K22 weight.
                    if (x(k) > 25) then
                        r(i, k, a) = 0
                        cycle
                    end if
                    if (nband(i, a) > 0) then
                        r(i, k, a) = orbit_r(i, a, x(k), ks, &
                            start=band(band_of(i, a, x(k)), i, a))
                    else
                        r(i, k, a) = orbit_r(i, a, x(k), ks)
                    end if
                    lost(i, k, a) = orbit_outside(ks)
                    if (lost(i, k, a)) then
                        r(i, k, a) = 0
                        cycle
                    end if
                    if (ks /= HP_OK) return
                    nfac(i, k, a) = 1
                end do
                st = DEVICE_OK
                return
            end if
            if (.not. thin_legendre(res)) then
                do i = first, last
                    if (nband(i, a) > 0) then
                        r(i, k, a) = thin_r(i, a, x(k), ks, &
                            start=band(band_of(i, a, x(k)), i, a))
                    else
                        r(i, k, a) = thin_r(i, a, x(k), ks)
                    end if
                    if (ks /= HP_OK) return
                    nfac(i, k, a) = 1
                end do
                st = DEVICE_OK
                return
            end if
            fresh = .true.
            do i = first, last
                v = sqrt(2*prof%temperature(i, a)/species(a)%mass*x(k))
                nu = deflection_frequency(species, prof%density(i, :), &
                    prof%temperature(i, :), prof%kappa(i, :), a, x(k))
                call ntv_f10(sf(i), v, nu, res%n_theta, res%n_leg, f10, ks)
                if (ks /= NTV_OK) return
                if (.not. fresh) then
                    call ntv_kernel(sf(i), v, nu, -prof%omega_te(i), &
                        species(a)%mass*v*v/(2*species(a)%z*elementary), res%n_theta, &
                        res%n_leg, kern, ks, indicators=.false., f10=f10, &
                        precondition=factor)
                    fresh = ks /= NTV_OK .or. kern%iterations > res%max_iterations
                    nits(i, k, a) = kern%iterations
                end if
                if (fresh) then
                    call ntv_kernel(sf(i), v, nu, -prof%omega_te(i), &
                        species(a)%mass*v*v/(2*species(a)%z*elementary), res%n_theta, &
                        res%n_leg, kern, ks, indicators=.false., f10=f10, keep=factor)
                    if (ks /= NTV_OK) return
                    nfac(i, k, a) = 1
                    fresh = .false.
                end if
                r(i, k, a) = kern%r
            end do
            st = DEVICE_OK
        end subroutine sweep

        !> R of species a on surface i at x = v^2/v_th^2, thin orbits, on the
        !> spectral elements of the Boozer line.
        real(dp) function thin_r(i, a, xv, st, start, used, mode, info) result(rv)
            integer, intent(in) :: i, a
            real(dp), intent(in) :: xv
            integer, intent(out) :: st
            type(hp_mesh_t), intent(in), optional :: start
            type(hp_mesh_t), intent(out), optional :: used
            integer, intent(in), optional :: mode
            type(adapt_info_t), intent(out), optional :: info
            type(hp_line_field_t) :: g
            type(hp_operator_t) :: op
            type(hp_mesh_t) :: mesh
            real(dp) :: v, kappa
            rv = 0
            v = sqrt(2*prof%temperature(i, a)/species(a)%mass*xv)
            call hp_line_field(sf(i)%line, 4*res%n_theta + 1, g, kappa, st)
            if (st /= HP_OK) return
            op%k = kappa
            op%k_power = 1
            op%q = sf(i)%q
            op%n = n_tor
            op%m = sf(i)%m
            op%c = sf(i)%c
            op%v = v
            op%wsign = 1
            op%omega_e = -prof%omega_te(i)
            op%pre = species(a)%mass*v*v/(2*species(a)%z*elementary)
            op%nu = deflection_frequency(species, prof%density(i, :), &
                prof%temperature(i, :), prof%kappa(i, :), a, xv)/v
            op%source_dpar = res%sources
            op%source_dlf10 = res%sources
            call start_mesh(g, mesh, start)
            call adapt_kernel(g, mesh, op, rv, st, kernel_mode(start, mode), info)
            if (present(used)) call normalised(g, mesh, used)
        end function thin_r

        !> The refinement mode of a kernel: as given, else MESH_FIXED on a start
        !> mesh and MESH_SETTLE from the mesh of the resolution.
        pure integer function kernel_mode(start, mode) result(m)
            type(hp_mesh_t), intent(in), optional :: start
            integer, intent(in), optional :: mode
            m = MESH_SETTLE
            if (present(start)) m = MESH_FIXED
            if (present(mode)) m = mode
        end function kernel_mode

        !> The mesh of the resolution, or start carried to g's X point and
        !> separatrix level.
        subroutine start_mesh(g, mesh, start)
            class(hp_geometry_t), intent(in) :: g
            type(hp_mesh_t), intent(out) :: mesh
            type(hp_mesh_t), intent(in), optional :: start
            if (.not. present(start)) then
                ! Adaptation starts from at least 16 columns and rows: a p + 1
                ! estimate does not see a resonance narrower than its element,
                ! and from 4 rows the refinement settled 8 % off on AUG 30835.
                if (res%hp_adapt > 0) then
                    call hp_default_mesh(g, res%hp_order, max(res%hp_columns, 16), &
                        max(res%hp_rows, 16), res%hp_delta, 0.3_dp, mesh)
                else
                    call hp_default_mesh(g, res%hp_order, res%hp_columns, res%hp_rows, &
                        res%hp_delta, 0.3_dp, mesh)
                end if
                return
            end if
            mesh%p = start%p
            mesh%theta = start%theta + g%theta_x
            mesh%lam = start%lam*g%lambda_c
            mesh%zeta = start%zeta
        end subroutine start_mesh

        subroutine normalised(g, mesh, out)
            class(hp_geometry_t), intent(in) :: g
            type(hp_mesh_t), intent(in) :: mesh
            type(hp_mesh_t), intent(out) :: out
            out%p = mesh%p
            out%theta = mesh%theta - g%theta_x
            out%lam = mesh%lam/g%lambda_c
            out%zeta = mesh%zeta
        end subroutine normalised

        !> hp_ntv_kernel on mesh, or with res%hp_adapt the corrected kernel
        !> R + eta of goal-oriented refinement from it (at most 8 estimates),
        !> refined until R + eta has settled and eta is within the tolerance.
        subroutine adapt_kernel(g, mesh, op, rv, st, mode, info)
            class(hp_geometry_t), intent(in) :: g
            type(hp_mesh_t), intent(inout) :: mesh
            type(hp_operator_t), intent(in) :: op
            real(dp), intent(out) :: rv
            integer, intent(out) :: st
            ! MESH_FIXED: solve on mesh. MESH_CHECK: refine until eta is
            ! within hp_adapt of R + eta, so that R on the mesh is.
            ! MESH_SETTLE: in addition until R + eta changes by at most hp_adapt.
            integer, intent(in) :: mode
            type(adapt_info_t), intent(out), optional :: info
            real(dp) :: eta, val, prev, fraction
            logical :: settled, full
            real(dp), allocatable :: eta_el(:)
            type(hp_mesh_t) :: next
            integer :: it, n_dof
            if (mode == MESH_FIXED .or. .not. res%hp_adapt > 0) then
                call hp_ntv_kernel(g, mesh, op, rv, st, n_dof)
                if (present(info)) info = adapt_info_t(0.0_dp, n_dof, .true.)
                return
            end if
            prev = huge(1.0_dp)
            settled = .false.
            do it = 1, 8
                call hp_ntv_estimate(g, mesh, op, rv, eta, eta_el, st, n_dof)
                if (st /= HP_OK) return
                val = rv + eta
                ! The other energies of a band are solved on this mesh without
                ! the correction eta, so eta has to be small in both modes.
                settled = abs(eta) <= res%hp_adapt*abs(val)
                if (mode == MESH_SETTLE) settled = settled &
                    .and. abs(val - prev) <= res%hp_adapt*abs(val)
                if (settled) exit
                prev = val
                ! Mark less while the refined mesh would exceed the size limit.
                fraction = 0.5_dp
                do
                    call hp_refine(mesh, eta_el, fraction, next)
                    full = real(n_dof, dp)*elements(next) &
                        > real(res%hp_max_unknowns, dp)*elements(mesh)
                    if (.not. full .or. fraction < 0.05_dp) exit
                    fraction = fraction/2
                end do
                if (full) exit
                mesh = next
            end do
            rv = val
            if (present(info)) info = adapt_info_t(abs(eta)/max(abs(val), tiny(val)), &
                n_dof, settled)
        end subroutine adapt_kernel

        pure real(dp) function elements(m)
            type(hp_mesh_t), intent(in) :: m
            elements = real(size(m%theta) - 1, dp)*(2*(size(m%lam) - 1) + size(m%zeta) - 1)
        end function elements

        !> R of species a on surface i at x = v^2/v_th^2 with finite orbit width.
        real(dp) function orbit_r(i, a, xv, st, start, used, mode, info) result(rv)
            integer, intent(in) :: i, a
            real(dp), intent(in) :: xv
            integer, intent(out) :: st
            type(hp_mesh_t), intent(in), optional :: start
            type(hp_mesh_t), intent(out), optional :: used
            integer, intent(in), optional :: mode
            type(adapt_info_t), intent(out), optional :: info
            type(orbit_particle_t) :: p
            type(hp_orbit_geometry_t) :: g
            type(hp_operator_t) :: op
            type(hp_mesh_t) :: mesh
            real(dp) :: v
            rv = 0
            v = sqrt(2*prof%temperature(i, a)/species(a)%mass*xv)
            p%mass = species(a)%mass
            p%charge = species(a)%z*elementary
            p%energy = p%mass*v*v/2
            p%s_ref = prof%s(i)
            p%omega_e = -prof%omega_te(i)
            if (ns >= 4) then
                call orbit_set_potential(p, prof%s, phi - phi(i), st)
                if (st /= 0) return
            end if
            call hp_orbit_geometry(eq, p, prof%s(i), b_ref(i), g, st)
            if (st /= 0) return
            call hp_orbit_perturbation(g, pert, n_tor, b_ref(i))
            ! Orbits leaving the file are outside the model and carry no NTV.
            st = DEVICE_ORBIT_LOST
            if (.not. hp_orbit_confined(g, 32, 32)) return
            op%k = g%kappa
            op%k_power = 1
            op%q = sf(i)%q
            op%n = n_tor
            op%m = sf(i)%m
            op%c = sf(i)%c
            op%v = v
            op%wsign = 1
            op%nu = deflection_frequency(species, prof%density(i, :), &
                prof%temperature(i, :), prof%kappa(i, :), a, xv)/v
            call start_mesh(g, mesh, start)
            call adapt_kernel(g, mesh, op, rv, st, kernel_mode(start, mode), info)
            if (present(used)) call normalised(g, mesh, used)
        end function orbit_r

        !> R of species a on surface i at x = v^2/v_th^2.
        real(dp) function kernel_r(i, a, xv, st, start, used, mode, info) result(rv)
            integer, intent(in) :: i, a
            real(dp), intent(in) :: xv
            integer, intent(out) :: st
            ! Adapted meshes (normalised) in and out and the refinement mode
            ! (adapt_kernel), hp kernels only.
            type(hp_mesh_t), intent(in), optional :: start
            type(hp_mesh_t), intent(out), optional :: used
            integer, intent(in), optional :: mode
            type(adapt_info_t), intent(out), optional :: info
            type(ntv_kernel_t) :: kern
            type(adapt_info_t) :: inf
            real(dp), allocatable :: f10(:, :)
            real(dp) :: v, nu
            if (res%orbit_width) then
                rv = orbit_r(i, a, xv, st, start, used, mode, inf)
                call note_unknowns(i, a, st, inf)
                if (present(info)) info = inf
                if (orbit_outside(st)) then
                    rv = 0
                    st = DEVICE_ORBIT_LOST
                end if
                return
            end if
            if (.not. thin_legendre(res)) then
                rv = thin_r(i, a, xv, st, start, used, mode, inf)
                call note_unknowns(i, a, st, inf)
                if (present(info)) info = inf
                if (st == HP_OK) st = NTV_OK
                return
            end if
            rv = 0
            v = sqrt(2*prof%temperature(i, a)/species(a)%mass*xv)
            nu = deflection_frequency(species, prof%density(i, :), prof%temperature(i, :), &
                prof%kappa(i, :), a, xv)
            call ntv_f10(sf(i), v, nu, res%n_theta, res%n_leg, f10, st)
            if (st /= NTV_OK) return
            call ntv_kernel(sf(i), v, nu, -prof%omega_te(i), &
                species(a)%mass*v*v/(2*species(a)%z*elementary), res%n_theta, res%n_leg, &
                kern, st, indicators=.false., f10=f10)
            rv = kern%r
        end function kernel_r

        !> Raises mesh_unknowns of surface i and species a to the mesh of a
        !> successful hp kernel, also one solved on a fixed band mesh.
        subroutine note_unknowns(i, a, st, info)
            integer, intent(in) :: i, a, st
            type(adapt_info_t), intent(in) :: info
            if (st /= HP_OK .or. .not. allocated(result)) return
            !$omp critical (ntv_unknowns)
            result(i, a)%mesh_unknowns = max(result(i, a)%mesh_unknowns, info%n_dof)
            !$omp end critical (ntv_unknowns)
        end subroutine note_unknowns

        !> int x^(5/2) e^-x (1, x) R dx by adaptive composite Gauss-Legendre
        !> panels: a panel is bisected while its 6- and 12-point rules differ
        !> by more than x_tol of the running total.
        subroutine adaptive(i, a, st)
            integer, intent(in) :: i, a
            integer, intent(out) :: st
            integer, parameter :: max_panels = 400
            real(dp), parameter :: start(11) = [0.0_dp, 0.5_dp, 1.0_dp, 2.0_dp, 3.0_dp, &
                4.5_dp, 6.5_dp, 9.0_dp, 13.0_dp, 18.0_dp, 25.0_dp]
            real(dp) :: lo(max_panels), hi(max_panels), val(3, max_panels), err(max_panels), &
                wt(2, max_panels)
            real(dp) :: x6(6), w6(6), x12(12), w12(12), total(3), t, scale, a1, a2, &
                pv(3, 2), pe(2), pw(2, 2)
            integer :: np, k, ks, nev
            st = DEVICE_SOLVE
            call gauss_legendre(6, x6, w6)
            call gauss_legendre(12, x12, w12)
            if (res%hp_adapt > 0 .and. .not. thin_legendre(res)) call seed_bands(i, a)
            np = size(start) - 1
            lo(1:np) = start(1:np)
            hi(1:np) = start(2:)
            call panel_eval(i, a, lo(1:np), hi(1:np), x6, w6, x12, w12, val(:, 1:np), &
                err(1:np), ks, wt(:, 1:np))
            if (ks /= NTV_OK) return
            nev = 18*np
            do
                total = sum(val(:, 1:np), 2)
                k = maxloc(err(1:np), 1)
                if (err(k) <= res%x_tol*maxval(abs(total))) exit
                if (np + 1 > max_panels) return
                np = np + 1
                lo(np) = (lo(k) + hi(k))/2
                hi(np) = hi(k)
                hi(k) = lo(np)
                call panel_eval(i, a, [lo(k), lo(np)], [hi(k), hi(np)], x6, w6, x12, w12, &
                    pv, pe, ks, pw)
                if (ks /= NTV_OK) return
                val(:, k) = pv(:, 1)
                val(:, np) = pv(:, 2)
                err(k) = pe(1)
                err(np) = pe(2)
                wt(:, k) = pw(:, 1)
                wt(:, np) = pw(:, 2)
                nev = nev + 36
            end do
            t = prof%temperature(i, a)
            scale = (t/(species(a)%z*elementary))**2/(4*sqrt(pi)*j0(i))
            result(i, a)%k11 = scale*total(1)
            result(i, a)%k12 = scale*total(2)
            result(i, a)%k21 = result(i, a)%k12
            result(i, a)%k22 = scale*total(3)
            a2 = prof%temperature_s(i, a)/(t*psi_pol_s(i))
            a1 = prof%density_s(i, a)/(prof%density(i, a)*psi_pol_s(i)) &
                + species(a)%z*elementary/t*prof%omega_te(i) - 1.5_dp*a2
            result(i, a)%dpsi_dr = dpsi_dr(i)
            result(i, a)%gamma_psi = -prof%density(i, a)*(result(i, a)%k11*a1 &
                + result(i, a)%k12*a2)
            result(i, a)%heat_flux = -prof%density(i, a)*t*(result(i, a)%k21*a1 &
                + result(i, a)%k22*a2)
            result(i, a)%torque = species(a)%z*elementary*result(i, a)%gamma_psi
            result(i, a)%offset_rotation = prof%omega_te(i) + result(i, a)%gamma_psi &
                /(prof%density(i, a)*result(i, a)%k11*species(a)%z*elementary/t)
            result(i, a)%factorisations = nev
            result(i, a)%iterations = np
            result(i, a)%lost_weight = sum(wt(2, 1:np))/sum(wt(1, 1:np))
            call record_kinks(i, a, lo(1:np), hi(1:np), start)
            st = DEVICE_OK
        end subroutine adaptive

        !> Meshes adapted at seed energies spanning the Maxwellian weight (or at
        !> the nodes), joined pairwise into bands between neighbouring seeds and
        !> checked at the band's geometric midpoint; energies whose orbits leave
        !> the file are skipped.
        subroutine seed_bands(i, a, nodes)
            integer, intent(in) :: i, a
            real(dp), intent(in), optional :: nodes(:)
            real(dp), allocatable :: seeds(:), checks(:), ok_x(:)
            type(hp_mesh_t), allocatable :: used(:), ok_mesh(:)
            type(adapt_info_t), allocatable :: info(:)
            integer, allocatable :: sts(:)
            integer :: k, n
            real(dp) :: r, dof
            if (present(nodes)) then
                seeds = nodes
            else
                seeds = [0.4_dp, 0.8_dp, 1.5_dp, 2.5_dp, 4.0_dp, 6.0_dp, 9.0_dp, 14.0_dp]
            end if
            allocate (used(size(seeds)), sts(size(seeds)), info(size(seeds)))
            !$omp parallel do schedule(dynamic, 1) default(shared) private(k, r)
            do k = 1, size(seeds)
                r = kernel_r(i, a, seeds(k), sts(k), used=used(k), info=info(k))
            end do
            !$omp end parallel do
            call note_adaptation(i, a, seeds, sts, info)
            ok_x = pack(seeds, sts == NTV_OK)
            allocate (ok_mesh(size(ok_x)))
            n = 0
            do k = 1, size(seeds)
                if (sts(k) /= NTV_OK) cycle
                n = n + 1
                ok_mesh(n) = used(k)
            end do
            if (n == 0) return
            dof = per_element(used, sts, info)
            n = min(n, max_bands + 1)
            ! One seed: one band, checked at the seed itself.
            nband(i, a) = max(1, n - 1)
            band_seed(1:n, i, a) = ok_x(1:n)
            if (n == 1) then
                band(1, i, a) = ok_mesh(1)
                checks = ok_x(1:1)
            else
                do k = 1, n - 1
                    band(k, i, a) = ok_mesh(k)
                    call join_within(band(k, i, a), ok_mesh(k + 1), dof)
                end do
                checks = sqrt(ok_x(1:n - 1)*ok_x(2:n))
            end if
            deallocate (used, sts, info)
            allocate (used(size(checks)), sts(size(checks)), info(size(checks)))
            !$omp parallel do schedule(dynamic, 1) default(shared) private(k, r)
            do k = 1, size(checks)
                r = kernel_r(i, a, checks(k), sts(k), start=band(k, i, a), used=used(k), &
                    mode=MESH_CHECK, info=info(k))
            end do
            !$omp end parallel do
            call note_adaptation(i, a, checks, sts, info)
            dof = max(dof, per_element(used, sts, info))
            do k = 1, size(checks)
                if (sts(k) == NTV_OK .and. allocated(used(k)%theta)) &
                    call join_within(band(k, i, a), used(k), dof)
            end do
        end subroutine seed_bands

        !> The largest number of unknowns per element of the meshes used.
        real(dp) function per_element(used, sts, info) result(dof)
            type(hp_mesh_t), intent(in) :: used(:)
            integer, intent(in) :: sts(:)
            type(adapt_info_t), intent(in) :: info(:)
            integer :: k
            dof = 0
            do k = 1, size(used)
                if (sts(k) == NTV_OK .and. allocated(used(k)%theta)) &
                    dof = max(dof, info(k)%n_dof/elements(used(k)))
            end do
        end function per_element

        !> Joins m into mesh unless the union, at dof unknowns per element,
        !> would exceed hp_max_unknowns. The union of two meshes has up to the
        !> sum of their columns times the sum of their rows, and the energies
        !> of a band are solved on it in parallel without a size check; past
        !> the limit the finer mesh is kept and the band's check energy
        !> refines it within the limit.
        subroutine join_within(mesh, m, dof)
            type(hp_mesh_t), intent(inout) :: mesh
            type(hp_mesh_t), intent(in) :: m
            real(dp), intent(in) :: dof
            type(hp_mesh_t) :: joined
            joined = mesh
            call join(joined, m)
            if (dof*elements(joined) <= real(res%hp_max_unknowns, dp)) then
                mesh = joined
            else if (elements(m) > elements(mesh)) then
                mesh = m
            end if
        end subroutine join_within

        !> Records the adaptation at energies x in mesh_error and unsettled of
        !> surface i and species a, and reports each unsettled energy.
        subroutine note_adaptation(i, a, x, sts, info)
            integer, intent(in) :: i, a, sts(:)
            real(dp), intent(in) :: x(:)
            type(adapt_info_t), intent(in) :: info(:)
            integer :: k
            do k = 1, size(x)
                if (sts(k) /= NTV_OK) cycle
                result(i, a)%mesh_error = max(result(i, a)%mesh_error, info(k)%eta_rel)
                if (info(k)%settled) cycle
                result(i, a)%unsettled = result(i, a)%unsettled + 1
                write (error_unit, '(a,f10.6,a,i0,a,f8.4,a,es10.3,a,i0)') &
                    'mesh not settled: s ', prof%s(i), ' species ', a, ' x ', x(k), &
                    ' |eta|/|R + eta| ', info(k)%eta_rel, ' unknowns ', info(k)%n_dof
            end do
        end subroutine note_adaptation

        !> The band of energy x: band k spans seeds k to k + 1, the outer bands
        !> extend to all lower and higher energies.
        pure integer function band_of(i, a, x) result(k)
            integer, intent(in) :: i, a
            real(dp), intent(in) :: x
            k = 1
            if (nband(i, a) < 2) return
            k = count(band_seed(2:nband(i, a), i, a) <= x) + 1
            k = min(max(k, 1), nband(i, a))
        end function band_of

        !> Adds the breakpoints of m to mesh.
        pure subroutine join(mesh, m)
            type(hp_mesh_t), intent(inout) :: mesh
            type(hp_mesh_t), intent(in) :: m
            mesh%theta = merged(mesh%theta, m%theta)
            mesh%lam = merged(mesh%lam, m%lam)
            mesh%zeta = merged(mesh%zeta, m%zeta)
        end subroutine join

        !> Sorted union of two increasing breakpoint lists; points closer than
        !> 1e-12 count as one.
        pure function merged(u, v) result(w)
            real(dp), intent(in) :: u(:), v(:)
            real(dp), allocatable :: w(:)
            real(dp) :: tmp(size(u) + size(v)), x
            integer :: i, j, n
            i = 1
            j = 1
            n = 0
            do while (i <= size(u) .or. j <= size(v))
                if (j > size(v)) then
                    x = u(i)
                    i = i + 1
                else if (i > size(u)) then
                    x = v(j)
                    j = j + 1
                else if (u(i) <= v(j)) then
                    x = u(i)
                    i = i + 1
                else
                    x = v(j)
                    j = j + 1
                end if
                if (n > 0) then
                    if (abs(x - tmp(n)) <= 1e-12_dp) cycle
                end if
                n = n + 1
                tmp(n) = x
            end do
            w = tmp(1:n)
        end function merged

        !> A kink lies in each run of adjacent panels the quadrature bisected at
        !> least once; its energy is the run's centre.
        subroutine record_kinks(i, a, lo, hi, start)
            integer, intent(in) :: i, a
            real(dp), intent(in) :: lo(:), hi(:), start(:)
            real(dp) :: l(size(lo)), h(size(lo)), w0
            integer :: order(size(lo)), k, m, j
            logical :: fine(size(lo))
            order = [(k, k=1, size(lo))]
            do k = 2, size(lo)
                m = order(k)
                j = k - 1
                do while (j >= 1)
                    if (lo(order(j)) <= lo(m)) exit
                    order(j + 1) = order(j)
                    j = j - 1
                end do
                order(j + 1) = m
            end do
            l = lo(order)
            h = hi(order)
            do k = 1, size(l)
                j = count(start <= l(k) + 1e-12_dp)
                w0 = start(min(j + 1, size(start))) - start(max(j, 1))
                fine(k) = h(k) - l(k) < 0.75_dp*w0
            end do
            k = 1
            do while (k <= size(l))
                if (.not. fine(k)) then
                    k = k + 1
                    cycle
                end if
                m = k
                do while (m < size(l))
                    if (.not. fine(m + 1)) exit
                    m = m + 1
                end do
                if (nkink(i, a) < max_kinks) then
                    nkink(i, a) = nkink(i, a) + 1
                    kink(nkink(i, a), i, a) = (l(k) + h(m))/2
                end if
                k = m + 1
            end do
        end subroutine record_kinks

        !> 6- and 12-point Gauss-Legendre values of int x^(5/2) e^-x (1, x) R
        !> on the panels [lo, hi] and their differences.
        subroutine panel_eval(i, a, lo, hi, x6, w6, x12, w12, fine, err, ks, weight)
            integer, intent(in) :: i, a
            real(dp), intent(in) :: lo(:), hi(:), x6(:), w6(:), x12(:), w12(:)
            real(dp), intent(out) :: fine(:, :), err(:)
            integer, intent(out) :: ks
            ! K_11 weight of each panel and its part on lost orbits.
            real(dp), intent(out) :: weight(:, :)
            real(dp) :: xs(size(x6) + size(x12), size(lo)), rs(size(xs, 1), size(lo))
            integer :: q, j, nq, sts(size(xs, 1), size(lo))
            nq = size(xs, 1)
            do j = 1, size(lo)
                xs(:, j) = (lo(j) + hi(j))/2 + (hi(j) - lo(j))/2*[x6, x12]
            end do
            ! The energies of all panels are independent and take the threads
            ! together (one surface and species at a time).
            !$omp parallel do schedule(dynamic, 1) default(shared) private(q, j)
            do q = 1, nq*size(lo)
                j = (q - 1)/nq + 1
                if (nband(i, a) > 0) then
                    rs(q - (j - 1)*nq, j) = kernel_r(i, a, xs(q - (j - 1)*nq, j), &
                        sts(q - (j - 1)*nq, j), &
                        start=band(band_of(i, a, xs(q - (j - 1)*nq, j)), i, a))
                else
                    rs(q - (j - 1)*nq, j) = kernel_r(i, a, xs(q - (j - 1)*nq, j), &
                        sts(q - (j - 1)*nq, j))
                end if
            end do
            !$omp end parallel do
            do j = 1, size(lo)
                call panel_sum(lo(j), hi(j), xs(:, j), rs(:, j), sts(:, j), w6, w12, &
                    fine(:, j), err(j), ks, weight(:, j))
                if (ks /= NTV_OK) return
            end do
        end subroutine panel_eval

        subroutine panel_sum(lo, hi, xs, rs, sts, w6, w12, fine, err, ks, weight)
            real(dp), intent(in) :: lo, hi, xs(:), rs(:), w6(:), w12(:)
            integer, intent(in) :: sts(:)
            real(dp), intent(out) :: fine(3), err, weight(2)
            integer, intent(out) :: ks
            real(dp) :: h, coarse(3), g
            integer :: q, n6
            h = (hi - lo)/2
            n6 = size(w6)
            coarse = 0
            fine = 0
            err = 0
            weight = 0
            do q = 1, size(xs)
                ks = sts(q)
                if (q > n6) then
                    weight(1) = weight(1) + h*w12(q - n6)*xs(q)**2.5_dp*exp(-xs(q))
                    if (ks == DEVICE_ORBIT_LOST) weight(2) = weight(2) &
                        + h*w12(q - n6)*xs(q)**2.5_dp*exp(-xs(q))
                end if
                if (ks == DEVICE_ORBIT_LOST) ks = NTV_OK
                if (ks /= NTV_OK) return
                if (q <= n6) then
                    g = h*w6(q)*xs(q)**2.5_dp*exp(-xs(q))*rs(q)
                    coarse = coarse + g*[1.0_dp, xs(q), xs(q)**2]
                else
                    g = h*w12(q - n6)*xs(q)**2.5_dp*exp(-xs(q))*rs(q)
                    fine = fine + g*[1.0_dp, xs(q), xs(q)**2]
                end if
            end do
            err = maxval(abs(fine - coarse))
        end subroutine panel_sum

        subroutine assemble(i, a)
            integer, intent(in) :: i, a
            real(dp) :: t, scale, a2, a1
            t = prof%temperature(i, a)
            scale = (t/(species(a)%z*elementary))**2/(4*sqrt(pi)*j0(i))
            result(i, a)%k11 = scale*sum(w*r(i, :, a))
            result(i, a)%k12 = scale*sum(w*x*r(i, :, a))
            result(i, a)%k21 = result(i, a)%k12
            result(i, a)%k22 = scale*sum(w*x**2*r(i, :, a))
            a2 = prof%temperature_s(i, a)/(t*psi_pol_s(i))
            a1 = prof%density_s(i, a)/(prof%density(i, a)*psi_pol_s(i)) &
                + species(a)%z*elementary/t*prof%omega_te(i) - 1.5_dp*a2
            result(i, a)%dpsi_dr = dpsi_dr(i)
            result(i, a)%gamma_psi = -prof%density(i, a)*(result(i, a)%k11*a1 &
                + result(i, a)%k12*a2)
            result(i, a)%heat_flux = -prof%density(i, a)*t*(result(i, a)%k21*a1 &
                + result(i, a)%k22*a2)
            result(i, a)%torque = species(a)%z*elementary*result(i, a)%gamma_psi
            result(i, a)%offset_rotation = prof%omega_te(i) + result(i, a)%gamma_psi &
                /(prof%density(i, a)*result(i, a)%k11*species(a)%z*elementary/t)
            result(i, a)%lost_weight = sum(w, mask=lost(i, :, a))/sum(w)
            result(i, a)%factorisations = sum(nfac(i, :, a))
            result(i, a)%iterations = sum(nits(i, :, a))
        end subroutine assemble

    end subroutine ntv_device_profile

    !> Reads the ntv_profile input: equilibrium and perturbation paths, then
    !> n_tor n_theta n_leg n_x chunk max_iterations sources(0|1), an optional
    !> line 'orbit_width [order columns rows delta adapt]' for the finite-width
    !> kernel, 'hp [order columns rows delta adapt]' for the mesh of the thin kernel
    !> (the default) or 'legendre' for the Legendre thin kernel (n_leg applies
    !> to it only), an optional
    !> line collision model x_max elements order [x_tol], the species (Z,
    !> mass [kg]) and per surface s Omega_tE [rad/s] with one line per
    !> species: density [m^-3] T [J] dn/ds dT/ds kappa [1/m].
    subroutine ntv_read_input(path, eq_path, pert_path, n_tor, species, prof, res)
        character(*), intent(in) :: path
        character(*), intent(out) :: eq_path, pert_path
        integer, intent(out) :: n_tor
        type(ntv_species_t), allocatable, intent(out) :: species(:)
        type(ntv_profile_t), intent(out) :: prof
        type(ntv_resolution_t), intent(out) :: res
        character(1024) :: line
        integer :: u, nsp, ns, i, a, status, sources
        open (newunit=u, file=path, status='old', action='read')
        read (u, '(a)') eq_path
        read (u, '(a)') pert_path
        read (u, *) n_tor, res%n_theta, res%n_leg, res%n_x, res%chunk, res%max_iterations, &
            sources
        res%sources = sources /= 0
        read (u, '(a)') line
        if (index(adjustl(line), 'orbit_width') == 1) then
            res%orbit_width = .true.
            call read_mesh(line(index(line, 'orbit_width') + 11:), res)
            read (u, '(a)') line
        else if (index(adjustl(line), 'hp') == 1) then
            call read_mesh(line(index(line, 'hp') + 2:), res)
            read (u, '(a)') line
        else if (index(adjustl(line), 'legendre') == 1) then
            res%thin_legendre = .true.
            read (u, '(a)') line
        end if
        read (line, *, iostat=status) res%collisions, res%x_max, res%nelem, res%order, &
            res%x_tol
        if (status /= 0) read (line, *, iostat=status) res%collisions, res%x_max, &
            res%nelem, res%order
        if (status == 0) then
            read (u, *) nsp
        else
            res%collisions = 0
            read (line, *) nsp
        end if
        allocate (species(nsp))
        do a = 1, nsp
            read (u, *) species(a)%z, species(a)%mass
        end do
        read (u, *) ns
        allocate (prof%s(ns), prof%omega_te(ns), prof%density(ns, nsp), &
            prof%temperature(ns, nsp), prof%density_s(ns, nsp), &
            prof%temperature_s(ns, nsp), prof%kappa(ns, nsp))
        do i = 1, ns
            read (u, *) prof%s(i), prof%omega_te(i)
            do a = 1, nsp
                read (u, *) prof%density(i, a), prof%temperature(i, a), &
                    prof%density_s(i, a), prof%temperature_s(i, a), prof%kappa(i, a)
            end do
        end do
        close (u)
    end subroutine ntv_read_input

    !> Optional mesh values 'order columns rows delta [adapt]' of an input line.
    subroutine read_mesh(text, res)
        character(*), intent(in) :: text
        type(ntv_resolution_t), intent(inout) :: res
        integer :: order, columns, rows, status
        real(dp) :: delta, adapt
        read (text, *, iostat=status) order, columns, rows, delta
        if (status /= 0) return
        res%hp_order = order
        res%hp_columns = columns
        res%hp_rows = rows
        res%hp_delta = delta
        read (text, *, iostat=status) order, columns, rows, delta, adapt
        if (status == 0) res%hp_adapt = adapt
    end subroutine read_mesh

    !> Axisymmetric transport matrix on every surface of the profile with
    !> the collision model res%collisions (nonzero) on the speed grid of res;
    !> the perturbation file only completes the surface geometry. One OpenMP
    !> unit per surface; a failed surface keeps its status and the others.
    subroutine ntv_device_axisymmetric(eq, pert, n_tor, species, prof, res, result, status)
        type(boozer_file_t), intent(in) :: eq, pert
        integer, intent(in) :: n_tor
        type(ntv_species_t), intent(in) :: species(:)
        type(ntv_profile_t), intent(in) :: prof
        type(ntv_resolution_t), intent(in) :: res
        type(ntv_axisymmetric_result_t), allocatable, intent(out) :: result(:)
        integer, intent(out) :: status
        type(ntv_surface_t), allocatable :: sf(:)
        type(collop_grid_t) :: grid
        real(dp) :: j0
        integer :: ns, nsp, i, st
        ns = size(prof%s)
        nsp = size(species)
        status = DEVICE_BAD_INPUT
        if (ns < 1 .or. nsp < 1 .or. res%collisions == 0) return
        call collop_make_grid(res%x_max, res%nelem, res%order, grid, st)
        if (st /= COLLOP_OK) return
        allocate (sf(ns), result(ns))
        do i = 1, ns
            call surface(eq, pert, n_tor, prof%s(i), res, sf(i), j0, result(i)%dpsi_dr, &
                result(i)%psi_pol_s, status, result(i)%g_cov, result(i)%b2, result(i)%g_phi)
            if (status /= DEVICE_OK) return
        end do
        !$omp parallel do schedule(dynamic, 1) default(shared) private(i)
        do i = 1, ns
            call axisymmetric_surface(i)
        end do
        !$omp end parallel do
        status = merge(DEVICE_OK, DEVICE_SOLVE, all(result%status == DEVICE_OK))
    contains
        subroutine axisymmetric_surface(i)
            integer, intent(in) :: i
            type(collop_species_t) :: sp(nsp)
            real(dp) :: lnl(nsp, nsp)
            integer :: a, its, ks
            do a = 1, nsp
                sp(a) = collop_species_t(species(a)%mass, species(a)%z, prof%density(i, a), &
                    prof%temperature(i, a))
            end do
            do a = 1, nsp
                lnl(a, :) = 3*sqrt(pi)/8*collop_vth(sp(a))*prof%kappa(i, a) &
                    /collop_nuhat(sp(a), sp(a), 1.0_dp)
            end do
            allocate (result(i)%d(nsp, nsp, 3, 3))
            result(i)%q = sf(i)%q
            call ntv_axisymmetric_coefficients(sf(i), sp, lnl, grid, res%collisions, &
                res%n_theta, res%n_leg, result(i)%d, its, ks, result(i)%kappa)
            if (ks /= NTV_OK) then
                result(i)%status = DEVICE_SOLVE
                result(i)%d = 0
                return
            end if
            call flows(result(i), i)
        end subroutine axisymmetric_surface

        subroutine flows(r, i)
            type(ntv_axisymmetric_result_t), intent(inout) :: r
            integer, intent(in) :: i
            real(dp) :: a1(nsp), a2(nsp), conv
            integer :: a
            do a = 1, nsp
                a2(a) = prof%temperature_s(i, a)/(prof%temperature(i, a)*r%psi_pol_s)
                a1(a) = prof%density_s(i, a)/(prof%density(i, a)*r%psi_pol_s) &
                    + species(a)%z*elementary/prof%temperature(i, a)*prof%omega_te(i) &
                    - 1.5_dp*a2(a)
            end do
            conv = -2*r%kappa*r%q*r%g_cov
            allocate (r%flow(nsp))
            do a = 1, nsp
                r%flow(a) = -conv*(sum(r%d(a, :, 3, 1)*a1) + sum(r%d(a, :, 3, 2)*a2))
            end do
            r%bootstrap = sum(species%z*elementary*prof%density(i, :)*r%flow)
            allocate (r%omega(nsp), r%k_pol(nsp))
            do a = 1, nsp
                r%omega(a) = prof%omega_te(i) + prof%temperature(i, a) &
                    /(species(a)%z*elementary)*(prof%density_s(i, a)/prof%density(i, a) &
                    + prof%temperature_s(i, a)/prof%temperature(i, a))/r%psi_pol_s
                r%k_pol(a) = r%g_phi*(r%flow(a) - r%omega(a)*r%g_phi)/r%b2
            end do
        end subroutine flows
    end subroutine ntv_device_axisymmetric

    !> Orbit-width indicators (ntv_indicator_t) on every surface of the profile.
    subroutine ntv_device_indicators(eq, pert, n_tor, species, prof, res, ind, status)
        type(boozer_file_t), intent(in) :: eq, pert
        integer, intent(in) :: n_tor
        type(ntv_species_t), intent(in) :: species(:)
        type(ntv_profile_t), intent(in) :: prof
        type(ntv_resolution_t), intent(in) :: res
        type(ntv_indicator_t), allocatable, intent(out) :: ind(:)
        integer, intent(out) :: status
        type(ntv_surface_t) :: sf
        real(dp) :: j0, dpsi_dr, psi_pol_s, g_cov, b2, ds_dr, ln, lt, va
        integer :: i, a, nsp
        nsp = size(species)
        allocate (ind(size(prof%s)))
        do i = 1, size(prof%s)
            call surface(eq, pert, n_tor, prof%s(i), res, sf, j0, dpsi_dr, psi_pol_s, status, &
                g_cov=g_cov, b2_avg=b2, eps=ind(i)%eps)
            if (status /= DEVICE_OK) return
            ds_dr = dpsi_dr/psi_pol_s
            allocate (ind(i)%gradient_length(nsp), ind(i)%excursion(nsp), ind(i)%banana(nsp))
            do a = 1, nsp
                ln = prof%density(i, a)/max(abs(prof%density_s(i, a)*ds_dr), tiny(1.0_dp))
                lt = prof%temperature(i, a)/max(abs(prof%temperature_s(i, a)*ds_dr), &
                    tiny(1.0_dp))
                ind(i)%gradient_length(a) = min(ln, lt)
                va = sqrt(2*prof%temperature(i, a)/species(a)%mass)
                ind(i)%excursion(a) = species(a)%mass*va*abs(g_cov) &
                    /(abs(species(a)%z)*elementary*sqrt(b2)*abs(dpsi_dr))
                ind(i)%banana(a) = 2*sqrt(2*ind(i)%eps)*ind(i)%excursion(a)
            end do
        end do
        status = DEVICE_OK
    end subroutine ntv_device_indicators

    !> E x B rotation Omega_tE [rad/s] of surface i from the radial force
    !> balance, given the toroidal rotation V_phi/R [rad/s] of species a at
    !> major radius r_major [m] (NEO-2 isw_calc_er = 1). V_phi/R is linear in
    !> Omega_tE through omega_a and through A_1 of every species in the flow,
    !> so two evaluations of the result's model give the root.
    function ntv_omega_te_from_rotation(r, species, prof, i, a, r_major, rotation) &
            result(omega_te)
        type(ntv_axisymmetric_result_t), intent(in) :: r
        type(ntv_species_t), intent(in) :: species(:)
        type(ntv_profile_t), intent(in) :: prof
        integer, intent(in) :: i, a
        real(dp), intent(in) :: r_major, rotation
        real(dp) :: omega_te, v0, v1
        v0 = model(0.0_dp)
        v1 = model(1.0_dp)
        omega_te = (rotation - v0)/(v1 - v0)
    contains
        real(dp) function model(om)
            real(dp), intent(in) :: om
            real(dp) :: a1(size(species)), a2(size(species)), flow, omega
            integer :: b
            do b = 1, size(species)
                a2(b) = prof%temperature_s(i, b)/(prof%temperature(i, b)*r%psi_pol_s)
                a1(b) = prof%density_s(i, b)/(prof%density(i, b)*r%psi_pol_s) &
                    + species(b)%z*elementary/prof%temperature(i, b)*om - 1.5_dp*a2(b)
            end do
            flow = 2*r%kappa*r%q*r%g_cov*(sum(r%d(a, :, 3, 1)*a1) + sum(r%d(a, :, 3, 2)*a2))
            omega = om + prof%temperature(i, a)/(species(a)%z*elementary) &
                *(prof%density_s(i, a)/prof%density(i, a) &
                + prof%temperature_s(i, a)/prof%temperature(i, a))/r%psi_pol_s
            model = omega + r%g_phi*(flow - omega*r%g_phi)/(r%b2*r_major**2)
        end function model
    end function ntv_omega_te_from_rotation

    !> Energy-coupled kernel per surface (one OpenMP unit per surface), with
    !> the species coupled through the collision operator:
    !> Gamma_psi of a = -n_a sum_b (K(a, b, 1) A_1 + K(a, b, 2) A_2) of b.
    !> The common Coulomb logarithm follows from kappa of each species.
    subroutine energy_profile(eq, n_tor, species, prof, res, sf, psi_pol_s, dpsi_dr, &
            result, status, kink, nkink)
        type(boozer_file_t), intent(in) :: eq
        integer, intent(in) :: n_tor
        type(ntv_species_t), intent(in) :: species(:)
        type(ntv_profile_t), intent(in) :: prof
        type(ntv_resolution_t), intent(in) :: res
        type(ntv_surface_t), intent(in) :: sf(:)
        real(dp), intent(in) :: psi_pol_s(:), dpsi_dr(:)
        type(ntv_species_result_t), allocatable, intent(inout) :: result(:, :)
        integer, intent(out) :: status
        real(dp), intent(in), optional :: kink(:, :, :)
        integer, intent(in), optional :: nkink(:, :)
        type(collop_grid_t) :: grid
        integer, allocatable :: ustat(:)
        integer :: i, ns, nsp
        ns = size(prof%s)
        nsp = size(species)
        status = DEVICE_BAD_INPUT
        if (n_tor == 0 .or. .not. allocated(eq%spline%x)) return
        call collop_make_grid(res%x_max, res%nelem, res%order, grid, i)
        if (i /= COLLOP_OK) return
        if (allocated(result)) deallocate (result)
        allocate (result(ns, nsp), ustat(ns))
        !$omp parallel do schedule(dynamic, 1) default(shared) private(i)
        do i = 1, ns
            call energy_surface(i, ustat(i))
        end do
        !$omp end parallel do
        do i = 1, ns
            if (ustat(i) /= DEVICE_OK) result(i, :) = ntv_species_result_t(status=ustat(i))
        end do
        status = merge(DEVICE_OK, DEVICE_SOLVE, all(ustat == DEVICE_OK))
    contains
        subroutine energy_surface(i, st)
            integer, intent(in) :: i
            integer, intent(out) :: st
            type(collop_species_t) :: sp(nsp)
            type(ntv_surface_t) :: local
            type(collop_grid_t) :: local_grid
            real(dp) :: lnl(nsp, nsp), kk(nsp, nsp, 2), qq(nsp, nsp, 2), a1(nsp), a2(nsp)
            integer :: a, b, its, ks
            st = DEVICE_SOLVE
            do a = 1, nsp
                sp(a) = collop_species_t(species(a)%mass, species(a)%z, prof%density(i, a), &
                    prof%temperature(i, a))
            end do
            do a = 1, nsp
                lnl(a, :) = 3*sqrt(pi)/8*collop_vth(sp(a))*prof%kappa(i, a) &
                    /collop_nuhat(sp(a), sp(a), 1.0_dp)
            end do
            local = sf(i)
            if (present(kink)) then
                call kink_grid(i, local_grid, ks)
                if (ks /= COLLOP_OK) return
            else
                local_grid = grid
            end if
            call ntv_energy_coefficients(local, sp, lnl, local_grid, res%collisions, &
                -prof%omega_te(i), res%n_theta, res%n_leg, kk, its, ks, qq)
            if (ks == NTV_NO_CONVERGENCE) st = DEVICE_NO_CONVERGENCE
            if (ks /= NTV_OK) return
            do a = 1, nsp
                a2(a) = prof%temperature_s(i, a)/(prof%temperature(i, a)*psi_pol_s(i))
                a1(a) = prof%density_s(i, a)/(prof%density(i, a)*psi_pol_s(i)) &
                    + species(a)%z*elementary/prof%temperature(i, a)*prof%omega_te(i) &
                    - 1.5_dp*a2(a)
            end do
            do a = 1, nsp
                result(i, a)%k11 = kk(a, a, 1)
                result(i, a)%k12 = kk(a, a, 2)
                result(i, a)%k21 = qq(a, a, 1)
                result(i, a)%k22 = qq(a, a, 2)
                result(i, a)%dpsi_dr = dpsi_dr(i)
                result(i, a)%gamma_psi = 0
                result(i, a)%heat_flux = 0
                do b = 1, nsp
                    result(i, a)%gamma_psi = result(i, a)%gamma_psi &
                        - prof%density(i, a)*(kk(a, b, 1)*a1(b) + kk(a, b, 2)*a2(b))
                    result(i, a)%heat_flux = result(i, a)%heat_flux &
                        - prof%density(i, a)*prof%temperature(i, a) &
                        *(qq(a, b, 1)*a1(b) + qq(a, b, 2)*a2(b))
                end do
                result(i, a)%torque = species(a)%z*elementary*result(i, a)%gamma_psi
                result(i, a)%offset_rotation = prof%omega_te(i) + result(i, a)%gamma_psi &
                    /(prof%density(i, a)*sum(kk(a, :, 1)*species%z*elementary &
                    /prof%temperature(i, :)))
                result(i, a)%factorisations = 1
                result(i, a)%iterations = its
            end do
            st = DEVICE_OK
        end subroutine energy_surface

        !> Uniform speed elements with extra edges at the kink speeds sqrt(x) of
        !> every species on surface i; uniform edges close to a kink give way.
        subroutine kink_grid(i, g, st)
            integer, intent(in) :: i
            type(collop_grid_t), intent(out) :: g
            integer, intent(out) :: st
            real(dp), parameter :: gap = 0.12_dp
            real(dp) :: y(res%nelem + 1 + 16*nsp), e(res%nelem + 1 + 16*nsp), t
            integer :: a, k, n, m
            n = 0
            do a = 1, nsp
                do k = 1, nkink(i, a)
                    t = sqrt(kink(k, i, a))
                    if (t > gap .and. t < res%x_max - gap) then
                        n = n + 1
                        y(n) = t
                    end if
                end do
            end do
            do k = 0, res%nelem
                t = res%x_max*k/res%nelem
                if (k == 0 .or. k == res%nelem .or. all(abs(y(1:n) - t) > gap)) then
                    n = n + 1
                    y(n) = t
                end if
            end do
            do k = 2, n
                t = y(k)
                m = k - 1
                do while (m >= 1)
                    if (y(m) <= t) exit
                    y(m + 1) = y(m)
                    m = m - 1
                end do
                y(m + 1) = t
            end do
            m = 1
            e(1) = y(1)
            do k = 2, n
                if (y(k) - e(m) > gap/2) then
                    m = m + 1
                    e(m) = y(k)
                end if
            end do
            e(m) = res%x_max
            call collop_make_grid(res%x_max, m - 1, res%order, g, st, e(1:m))
        end subroutine kink_grid
    end subroutine energy_profile

    !> Kernel surface of label s: the equilibrium line, the perturbation
    !> harmonics of toroidal number n_tor resolved by the theta grid, and the
    !> flux-function factors of the torque.
    subroutine surface(eq, pert, n_tor, s, res, sf, j0, dpsi_dr, psi_pol_s, status, g_cov, &
            b2_avg, g_phi, eps, b_ref)
        type(boozer_file_t), intent(in) :: eq, pert
        integer, intent(in) :: n_tor
        real(dp), intent(in) :: s
        type(ntv_resolution_t), intent(in) :: res
        type(ntv_surface_t), intent(out) :: sf
        real(dp), intent(out) :: j0, dpsi_dr, psi_pol_s
        integer, intent(out) :: status
        real(dp), intent(out), optional :: g_cov, b2_avg, g_phi, eps, b_ref
        type(boozer_surface_t) :: es, ps
        type(ntv_boozer_line_t) :: line
        real(dp) :: theta, bmod, grad_s, b2, bmin, bmax
        logical, allocatable :: keep(:)
        integer :: st, it, k

        status = DEVICE_BAD_INPUT
        call boozer_surface(eq, s, es, st)
        if (st /= BOOZER_OK) return
        call ntv_boozer_line(es, line, st)
        if (st /= NTV_LINE_OK) return
        call boozer_surface(pert, s, ps, st)
        if (st /= BOOZER_OK) return
        ! Toroidal mode numbers of the file count per field period.
        keep = ps%n*ps%nper == n_tor .and. abs(ps%m) <= (res%n_theta - 1)/2 .and. &
            (ps%bmn /= 0 .or. ps%bma /= 0)
        if (.not. any(keep)) return
        ! The file holds delta B = Re sum_m (bmnc - i bmns) e^{i(m theta + n phi)}
        ! (NEO-2 convention), the kernel Re sum_m c_m e^{i(m theta - n phi)}.
        sf%m = -pack(ps%m, keep)
        sf%c = cmplx(pack(ps%bmn, keep), pack(ps%bma, keep), dp)/line%bref
        sf%n = n_tor
        sf%q = 1/es%iota
        sf%weight = NTV_WEIGHT_FULL
        sf%boozer_stream = .true.
        sf%source_dpar = res%sources
        sf%source_dlf10 = res%sources
        allocate (sf%line, source=line)

        j0 = 0
        bmin = huge(1.0_dp)
        bmax = 0
        do it = 0, res%n_theta - 1
            theta = 2*pi*it/res%n_theta
            bmod = 0
            do k = 1, size(es%m)
                bmod = bmod + es%bmn(k)*cos(es%m(k)*theta) + es%bma(k)*sin(es%m(k)*theta)
            end do
            j0 = j0 + (line%bref/bmod)**2
        end do
        if (present(eps)) then
            do it = 0, 16*res%n_theta - 1
                theta = 2*pi*it/(16*res%n_theta)
                bmod = 0
                do k = 1, size(es%m)
                    bmod = bmod + es%bmn(k)*cos(es%m(k)*theta) + es%bma(k)*sin(es%m(k)*theta)
                end do
                bmin = min(bmin, abs(bmod))
                bmax = max(bmax, abs(bmod))
            end do
            eps = (bmax - bmin)/(bmax + bmin)
        end if
        j0 = 2*pi*j0/res%n_theta
        psi_pol_s = es%iota*es%psi_prime
        call boozer_average(es, 2*res%n_theta, 1, grad_s, b2)
        dpsi_dr = psi_pol_s*grad_s
        if (present(g_cov)) g_cov = es%g_cov
        ! NEO-2 counts Omega_tE and V_phi along the toroidal direction in which
        ! the covariant field is -G of the file (ITER TC24 and AUG 30835).
        if (present(g_phi)) g_phi = -es%g_cov
        if (present(b2_avg)) b2_avg = b2
        if (present(b_ref)) b_ref = line%bref
        status = DEVICE_OK
    end subroutine surface

end module milca_ntv_device
