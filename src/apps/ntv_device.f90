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
!! Work units (species, energy node, chunk of consecutive surfaces) run
!! under OpenMP with dynamic scheduling. Each unit factorises its first
!! surface and solves the following ones by GMRES preconditioned with that
!! factorisation, refactorising when the iteration count exceeds a bound.
module milca_ntv_device
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_quadrature, only: gauss_gen_laguerre, gauss_legendre
    use milca_boozer_file, only: boozer_file_t, boozer_surface_t, boozer_surface, &
        boozer_average, BOOZER_OK
    use milca_ntv_boozer_line, only: ntv_boozer_line_t, ntv_boozer_line, NTV_LINE_OK
    use milca_dk_ntv, only: ntv_surface_t, ntv_kernel_t, ntv_factor_t, ntv_kernel, &
        ntv_f10, ntv_energy_coefficients, ntv_axisymmetric_coefficients, NTV_OK, &
            NTV_NO_CONVERGENCE, NTV_WEIGHT_FULL
    use milca_dk_collop, only: collop_species_t, collop_grid_t, collop_make_grid, &
        collop_nuhat, collop_vth, COLLOP_OK
    implicit none
    private
    public :: ntv_species_t, ntv_profile_t, ntv_resolution_t, ntv_species_result_t
    public :: ntv_device_profile, deflection_frequency, ntv_read_input, &
        ntv_device_axisymmetric, ntv_axisymmetric_result_t
    integer, parameter, public :: DEVICE_OK = 0, DEVICE_BAD_INPUT = 1, DEVICE_SOLVE = 2, &
        DEVICE_NO_CONVERGENCE = 3

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
        ! puts kinks into the integrand that fixed rules do not resolve.
        real(dp) :: x_tol = 0
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
    !> profile status is DEVICE_SOLVE when any of them failed.
    type :: ntv_species_result_t
        real(dp) :: k11 = 0, k12 = 0, k21 = 0, k22 = 0, dpsi_dr = 0, gamma_psi = 0, &
            heat_flux = 0, torque = 0, offset_rotation = 0
        integer :: factorisations = 0, iterations = 0, status = DEVICE_OK
    end type

    !> Axisymmetric transport matrix of one surface: d(a, b, i, j) of
    !> ntv_axisymmetric_coefficients (responses i = particle, energy, parallel
    !> flow of species a to the drives j = A_1, A_2, parallel field of species
    !> b), b.grad theta/b (kappa), q, dpsi_pol/dr [Wb/m] and the status.
    type :: ntv_axisymmetric_result_t
        real(dp), allocatable :: d(:, :, :, :)
        real(dp) :: kappa = 0, q = 0, dpsi_dr = 0
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

    subroutine ntv_device_profile(eq, pert, n_tor, species, prof, res, result, status)
        type(boozer_file_t), intent(in) :: eq, pert
        integer, intent(in) :: n_tor
        type(ntv_species_t), intent(in) :: species(:)
        type(ntv_profile_t), intent(in) :: prof
        type(ntv_resolution_t), intent(in) :: res
        type(ntv_species_result_t), allocatable, intent(out) :: result(:, :)
        integer, intent(out) :: status
        type(ntv_surface_t), allocatable :: sf(:)
        real(dp), allocatable :: x(:), w(:), r(:, :, :), j0(:), dpsi_dr(:), psi_pol_s(:)
        integer, allocatable :: units(:, :), nfac(:, :, :), nits(:, :, :), ustat(:)
        ! Energies x of the resonance kinks per surface and species, from the
        ! panels the adaptive quadrature refined.
        integer, parameter :: max_kinks = 16
        real(dp), allocatable :: kink(:, :, :)
        integer, allocatable :: nkink(:, :)
        integer :: ns, nsp, nchunk, nunit, i, a, k, c, u

        ns = size(prof%s)
        nsp = size(species)
        status = DEVICE_BAD_INPUT
        if (ns < 1 .or. nsp < 1 .or. res%n_x < 1 .or. res%chunk < 1) return
        allocate (sf(ns), j0(ns), dpsi_dr(ns), psi_pol_s(ns))
        do i = 1, ns
            call surface(eq, pert, n_tor, prof%s(i), res, sf(i), j0(i), dpsi_dr(i), &
                psi_pol_s(i), status)
            if (status /= DEVICE_OK) return
        end do
        if (res%x_tol > 0) then
            allocate (result(ns, nsp), ustat(ns*nsp), kink(max_kinks, ns, nsp), &
                nkink(ns, nsp))
            nkink = 0
            !$omp parallel do schedule(dynamic, 1) default(shared) private(u)
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
        allocate (x(res%n_x), w(res%n_x))
        call gauss_gen_laguerre(res%n_x, 2.5_dp, x, w)

        ! Units: (species, node, first surface of the chunk).
        nchunk = (ns + res%chunk - 1)/res%chunk
        nunit = nsp*res%n_x*nchunk
        allocate (units(3, nunit), ustat(nunit), r(ns, res%n_x, nsp), &
            nfac(ns, res%n_x, nsp), nits(ns, res%n_x, nsp))
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
        !$omp parallel do schedule(dynamic, 1) default(shared) private(u)
        do u = 1, nunit
            call sweep(units(1, u), units(2, u), units(3, u), &
                min(units(3, u) + res%chunk - 1, ns), ustat(u))
        end do
        !$omp end parallel do

        allocate (result(ns, nsp))
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

        !> R of species a on surface i at x = v^2/v_th^2.
        real(dp) function kernel_r(i, a, xv, st) result(rv)
            integer, intent(in) :: i, a
            real(dp), intent(in) :: xv
            integer, intent(out) :: st
            type(ntv_kernel_t) :: kern
            real(dp), allocatable :: f10(:, :)
            real(dp) :: v, nu
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

        !> int x^(5/2) e^-x (1, x) R dx by adaptive composite Gauss-Legendre
        !> panels: a panel is bisected while its 6- and 12-point rules differ
        !> by more than x_tol of the running total.
        subroutine adaptive(i, a, st)
            integer, intent(in) :: i, a
            integer, intent(out) :: st
            integer, parameter :: max_panels = 400
            real(dp), parameter :: start(11) = [0.0_dp, 0.5_dp, 1.0_dp, 2.0_dp, 3.0_dp, &
                4.5_dp, 6.5_dp, 9.0_dp, 13.0_dp, 18.0_dp, 25.0_dp]
            real(dp) :: lo(max_panels), hi(max_panels), val(3, max_panels), err(max_panels)
            real(dp) :: x6(6), w6(6), x12(12), w12(12), total(3), t, scale, a1, a2
            integer :: np, k, ks, nev
            st = DEVICE_SOLVE
            call gauss_legendre(6, x6, w6)
            call gauss_legendre(12, x12, w12)
            np = size(start) - 1
            lo(1:np) = start(1:np)
            hi(1:np) = start(2:)
            nev = 0
            do k = 1, np
                call panel_eval(i, a, lo(k), hi(k), x6, w6, x12, w12, val(:, k), err(k), ks)
                if (ks /= NTV_OK) return
                nev = nev + 18
            end do
            do
                total = sum(val(:, 1:np), 2)
                k = maxloc(err(1:np), 1)
                if (err(k) <= res%x_tol*maxval(abs(total))) exit
                if (np + 1 > max_panels) return
                np = np + 1
                lo(np) = (lo(k) + hi(k))/2
                hi(np) = hi(k)
                hi(k) = lo(np)
                call panel_eval(i, a, lo(k), hi(k), x6, w6, x12, w12, val(:, k), err(k), ks)
                if (ks /= NTV_OK) return
                call panel_eval(i, a, lo(np), hi(np), x6, w6, x12, w12, val(:, np), err(np), &
                    ks)
                if (ks /= NTV_OK) return
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
            call record_kinks(i, a, lo(1:np), hi(1:np), start)
            st = DEVICE_OK
        end subroutine adaptive

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
        !> on [lo, hi] and their difference.
        subroutine panel_eval(i, a, lo, hi, x6, w6, x12, w12, fine, err, ks)
            integer, intent(in) :: i, a
            real(dp), intent(in) :: lo, hi, x6(:), w6(:), x12(:), w12(:)
            real(dp), intent(out) :: fine(3), err
            integer, intent(out) :: ks
            real(dp) :: c, h, xv, rv, coarse(3), g
            integer :: q
            c = (lo + hi)/2
            h = (hi - lo)/2
            coarse = 0
            fine = 0
            err = 0
            do q = 1, size(x6)
                xv = c + h*x6(q)
                rv = kernel_r(i, a, xv, ks)
                if (ks /= NTV_OK) return
                g = h*w6(q)*xv**2.5_dp*exp(-xv)*rv
                coarse = coarse + g*[1.0_dp, xv, xv**2]
            end do
            do q = 1, size(x12)
                xv = c + h*x12(q)
                rv = kernel_r(i, a, xv, ks)
                if (ks /= NTV_OK) return
                g = h*w12(q)*xv**2.5_dp*exp(-xv)*rv
                fine = fine + g*[1.0_dp, xv, xv**2]
            end do
            err = maxval(abs(fine - coarse))
        end subroutine panel_eval

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
            result(i, a)%factorisations = sum(nfac(i, :, a))
            result(i, a)%iterations = sum(nits(i, :, a))
        end subroutine assemble

    end subroutine ntv_device_profile

    !> Reads the ntv_profile input: equilibrium and perturbation paths, then
    !> n_tor n_theta n_leg n_x chunk max_iterations sources(0|1), an optional
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
        real(dp) :: j0, psi_pol_s
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
                psi_pol_s, status)
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
            end if
        end subroutine axisymmetric_surface
    end subroutine ntv_device_axisymmetric

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
    subroutine surface(eq, pert, n_tor, s, res, sf, j0, dpsi_dr, psi_pol_s, status)
        type(boozer_file_t), intent(in) :: eq, pert
        integer, intent(in) :: n_tor
        real(dp), intent(in) :: s
        type(ntv_resolution_t), intent(in) :: res
        type(ntv_surface_t), intent(out) :: sf
        real(dp), intent(out) :: j0, dpsi_dr, psi_pol_s
        integer, intent(out) :: status
        type(boozer_surface_t) :: es, ps
        type(ntv_boozer_line_t) :: line
        real(dp) :: theta, bmod, grad_s, b2
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
        do it = 0, res%n_theta - 1
            theta = 2*pi*it/res%n_theta
            bmod = 0
            do k = 1, size(es%m)
                bmod = bmod + es%bmn(k)*cos(es%m(k)*theta) + es%bma(k)*sin(es%m(k)*theta)
            end do
            j0 = j0 + (line%bref/bmod)**2
        end do
        j0 = 2*pi*j0/res%n_theta
        psi_pol_s = es%iota*es%psi_prime
        call boozer_average(es, 2*res%n_theta, 1, grad_s, b2)
        dpsi_dr = psi_pol_s*grad_s
        status = DEVICE_OK
    end subroutine surface

end module milca_ntv_device
