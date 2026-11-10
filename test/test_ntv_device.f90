!> NTV torque profile on Boozer files (milca_ntv_device): the deflection
!! frequency in its Lorentz limit, and sweeps over neighbouring surfaces
!! preconditioned by one factorisation agree with a factorisation per surface,
!! a second solve at the reported offset rotation carries no flux, and a
!! surface that fails leaves the other surfaces of the profile intact, and the
!! radial force balance returns the Omega_tE of the rotation it was solved for;
!! the orbit-width indicators match the concentric-circle geometry of the file,
!! the thin kernel refined from a coarse mesh by its error estimate reaches the
!! converged Legendre kernel, and the finite-width kernel reproduces the thin
!! kernel for thin orbits and counts the energies whose orbits leave the file.
program test_ntv_device
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use milca_boozer_file, only: boozer_file_t, load_boozer_file, BOOZER_OK
    use milca_ntv_device, only: ntv_species_t, ntv_profile_t, ntv_resolution_t, &
        ntv_species_result_t, ntv_device_profile, deflection_frequency, DEVICE_OK, &
        DEVICE_SOLVE, ntv_axisymmetric_result_t, ntv_device_axisymmetric, &
        ntv_omega_te_from_rotation, ntv_indicator_t, ntv_device_indicators
    use milca_boozer_orbit, only: orbit_particle_t
    use milca_dk_hp, only: hp_operator_t
    use milca_ntv_orbit_geometry, only: hp_orbit_geometry_t, hp_orbit_geometry, &
        hp_orbit_perturbation, hp_orbit_confined
    implicit none
    real(dp), parameter :: r0 = 1.7_dp, a = 0.5_dp, b0 = 2.3_dp, eps = 0.1_dp
    real(dp), parameter :: jpol = -1.2e7_dp, itor_s = -3.0e5_dp, flux = -1.3_dp
    real(dp), parameter :: pi = 4*atan(1.0_dp), kev = 1.602176634e-16_dp
    integer :: failures = 0

    call lorentz_limit()
    call preconditioned_sweep()
    call energy_coupled()
    call offset_rotation()
    call failed_surface()
    call force_balance()
    call indicators()
    call hp_thin_matches_legendre()
    call hp_adaptive()
    call orbit_width_thin()
    call orbit_width_lost()
    if (failures > 0) error stop 'FAIL'
    print '(a)', 'All tests passed!'

contains

    !> Electrons on ions far slower than them: erf(y) = 1 and G(y) = 1/(2 y^2)
    !> to double precision at y = 148, so nu_D is the electron-ion rate
    !> (3 sqrt(pi)/8) v_e kappa_e n_i Z_i^2/n_e (1 - 1/(2 y^2))/x^(3/2).
    subroutine lorentz_limit()
        type(ntv_species_t) :: sp(2)
        real(dp) :: n(2), t(2), kappa(2), x, ve, want, nu_ee, y
        sp(1) = ntv_species_t(-1.0_dp, 9.1093837e-31_dp)
        sp(2) = ntv_species_t(2.0_dp, 6.6446573e-27_dp)
        n = [2e19_dp, 1e19_dp]
        t = [5*kev, 5*kev]
        kappa = [1e-4_dp, 0.0_dp]
        x = 3.0_dp
        ve = sqrt(2*t(1)/sp(1)%mass)
        nu_ee = deflection_frequency(sp(1:1), n(1:1), t(1:1), kappa(1:1), 1, x)
        y = sqrt(x*sp(2)%mass/sp(1)%mass)
        want = 3*sqrt(pi)/8*ve*kappa(1)*n(2)*4/n(1)*(1 - 1/(2*y*y))/x**1.5_dp
        call close_to('electron-ion deflection', &
            deflection_frequency(sp, n, t, kappa, 1, x) - nu_ee, want, 1e-12_dp)
    end subroutine lorentz_limit

    subroutine preconditioned_sweep()
        character(*), parameter :: eq_name = 'ntv_device_eq.bc', pert_name = 'ntv_device_pert.bc'
        type(boozer_file_t) :: eq, pert
        type(ntv_species_t) :: sp(2)
        type(ntv_profile_t) :: prof
        type(ntv_resolution_t) :: res
        type(ntv_species_result_t), allocatable :: direct(:, :), swept(:, :)
        integer :: status, i, ns
        call write_file(eq_name, .false.)
        call write_file(pert_name, .true.)
        call load_boozer_file(eq_name, eq, status)
        call check('equilibrium file', status == BOOZER_OK)
        call load_boozer_file(pert_name, pert, status)
        call check('perturbation file', status == BOOZER_OK)
        sp(1) = ntv_species_t(1.0_dp, 3.3435837768e-27_dp)
        sp(2) = ntv_species_t(-1.0_dp, 9.1093837e-31_dp)
        ns = 5
        prof%s = [(0.40_dp + 0.01_dp*i, i=1, ns)]
        prof%omega_te = [(2e3_dp*(1 - 0.1_dp*i), i=1, ns)]
        allocate (prof%density(ns, 2), prof%temperature(ns, 2), prof%density_s(ns, 2), &
            prof%temperature_s(ns, 2), prof%kappa(ns, 2))
        prof%density(:, 1) = 3e19_dp
        prof%density(:, 2) = 3e19_dp
        prof%temperature = 2*kev
        prof%density_s = -1e19_dp
        prof%temperature_s = -1*kev
        prof%kappa(:, 1) = 3e-4_dp
        prof%kappa(:, 2) = 6e-4_dp
        res%x_tol = 0
        res%n_theta = 33
        res%n_leg = 60
        res%n_x = 6
        res%chunk = 1
        ! The preconditioned sweep belongs to the Legendre kernel.
        res%thin_legendre = .true.
        call ntv_device_profile(eq, pert, 3, sp, prof, res, direct, status)
        call check('direct profile', status == DEVICE_OK)
        res%chunk = ns
        call ntv_device_profile(eq, pert, 3, sp, prof, res, swept, status)
        call check('swept profile', status == DEVICE_OK)
        if (failures > 0) return
        call check('one factorisation per sweep', &
            sum(swept%factorisations) < sum(direct%factorisations) .and. &
            sum(swept%iterations) > 0)
        do i = 1, ns
            call close_to('torque ions', swept(i, 1)%torque, direct(i, 1)%torque, 1e-8_dp)
            call close_to('torque electrons', swept(i, 2)%torque, direct(i, 2)%torque, &
                1e-8_dp)
            call close_to('K12 ions', swept(i, 1)%k12, direct(i, 1)%k12, 1e-8_dp)
        end do
    end subroutine preconditioned_sweep

    !> The energy-coupled kernel with Lorentz collisions reproduces the
    !> monoenergetic coefficients (lumped speed nodes against Gauss-Laguerre
    !> nodes), and the full operator, whose field-particle part couples the
    !> species, converges to coefficients that differ from pitch-angle
    !> scattering alone.
    subroutine energy_coupled()
        character(*), parameter :: eq_name = 'ntv_energy_eq.bc', pert_name = 'ntv_energy_pert.bc'
        type(boozer_file_t) :: eq, pert
        type(ntv_species_t) :: sp(2)
        type(ntv_profile_t) :: prof
        type(ntv_resolution_t) :: res
        type(ntv_species_result_t), allocatable :: mono(:, :), lorentz(:, :), full(:, :)
        integer :: status, a
        call write_file(eq_name, .false.)
        call write_file(pert_name, .true.)
        call load_boozer_file(eq_name, eq, status)
        call load_boozer_file(pert_name, pert, status)
        sp(1) = ntv_species_t(1.0_dp, 3.3435837768e-27_dp)
        sp(2) = ntv_species_t(-1.0_dp, 9.1093837e-31_dp)
        prof%s = [0.45_dp]
        prof%omega_te = [1e3_dp]
        allocate (prof%density(1, 2), prof%temperature(1, 2), prof%density_s(1, 2), &
            prof%temperature_s(1, 2), prof%kappa(1, 2))
        prof%density = 3e19_dp
        prof%temperature = 2*kev
        prof%density_s = -1e19_dp
        prof%temperature_s = -1*kev
        prof%kappa(1, :) = [3e-4_dp, 6e-4_dp]
        res%x_tol = 0
        res%n_theta = 17
        res%n_leg = 24
        res%n_x = 16
        res%sources = .false.
        ! The energy-coupled operator is Legendre in pitch; the same pitch
        ! basis isolates the speed treatment (24 degrees are far from
        ! converged here).
        res%thin_legendre = .true.
        call ntv_device_profile(eq, pert, 3, sp, prof, res, mono, status)
        call check('monoenergetic', status == DEVICE_OK)
        res%collisions = 1
        res%nelem = 6
        call ntv_device_profile(eq, pert, 3, sp, prof, res, lorentz, status)
        call check('energy-coupled Lorentz', status == DEVICE_OK)
        res%collisions = 3
        res%sources = .true.
        call ntv_device_profile(eq, pert, 3, sp, prof, res, full, status)
        call check('energy-coupled full operator', status == DEVICE_OK)
        if (failures > 0) return
        do a = 1, 2
            call close_to('Lorentz K11', lorentz(1, a)%k11, mono(1, a)%k11, 2e-2_dp)
            call close_to('Lorentz K12', lorentz(1, a)%k12, mono(1, a)%k12, 2e-2_dp)
            call close_to('Lorentz K22', lorentz(1, a)%k22, mono(1, a)%k22, 2e-2_dp)
            ! Onsager symmetry: without energy scattering each speed node
            ! carries both moments, so K_21 = K_12 to the iterative tolerance.
            call close_to('Lorentz K21 = K12', lorentz(1, a)%k21, lorentz(1, a)%k12, &
                1e-6_dp)
        end do
        call check('field-particle part acts', &
            abs(full(1, 1)%k11/lorentz(1, 1)%k11 - 1) > 1e-3_dp .and. full(1, 1)%iterations > 0)
    end subroutine energy_coupled

    !> Weak gradients keep the offset rotation small against the collision
    !> frequency, so the coefficients barely change between Omega_tE = 0 and
    !> the offset, and the second solve there must carry almost no flux.
    subroutine offset_rotation()
        character(*), parameter :: eq_name = 'ntv_offset_eq.bc', pert_name = 'ntv_offset_pert.bc'
        type(boozer_file_t), target :: eq, pert
        type(ntv_species_t) :: sp(1)
        type(ntv_profile_t) :: prof
        type(ntv_resolution_t) :: res
        type(ntv_species_result_t), allocatable :: rest(:, :), shifted(:, :)
        integer :: status
        call write_file(eq_name, .false.)
        call write_file(pert_name, .true.)
        call load_boozer_file(eq_name, eq, status)
        call load_boozer_file(pert_name, pert, status)
        sp(1) = ntv_species_t(1.0_dp, 3.3435837768e-27_dp)
        prof%s = [0.45_dp]
        prof%omega_te = [0.0_dp]
        allocate (prof%density(1, 1), prof%temperature(1, 1), prof%density_s(1, 1), &
            prof%temperature_s(1, 1), prof%kappa(1, 1))
        prof%density = 3e19_dp
        prof%temperature = 2*kev
        prof%density_s = -1e16_dp
        prof%temperature_s = -1e-3_dp*kev
        prof%kappa = 3e-3_dp
        res%x_tol = 0
        res%n_theta = 17
        res%n_leg = 24
        res%n_x = 16
        res%sources = .false.
        call ntv_device_profile(eq, pert, 3, sp, prof, res, rest, status)
        call check('offset: solve at rest', status == DEVICE_OK)
        if (failures > 0) return
        prof%omega_te = [rest(1, 1)%offset_rotation]
        call ntv_device_profile(eq, pert, 3, sp, prof, res, shifted, status)
        call check('offset: solve at the offset', status == DEVICE_OK)
        if (failures > 0) return
        print '(a,3es12.4)', 'offset', rest(1, 1)%offset_rotation, rest(1, 1)%gamma_psi, &
            shifted(1, 1)%gamma_psi
        call check('flux vanishes at the offset rotation', rest(1, 1)%gamma_psi /= 0 .and. &
            abs(shifted(1, 1)%gamma_psi) < 1e-3_dp*abs(rest(1, 1)%gamma_psi))
    end subroutine offset_rotation

    !> No density on the second surface makes its collision operator invalid;
    !> the first surface must still be solved and reported.
    subroutine failed_surface()
        character(*), parameter :: eq_name = 'ntv_failed_eq.bc', pert_name = 'ntv_failed_pert.bc'
        type(boozer_file_t) :: eq, pert
        type(ntv_species_t) :: sp(2)
        type(ntv_profile_t) :: prof
        type(ntv_resolution_t) :: res
        type(ntv_species_result_t), allocatable :: both(:, :), alone(:, :)
        integer :: status
        call write_file(eq_name, .false.)
        call write_file(pert_name, .true.)
        call load_boozer_file(eq_name, eq, status)
        call load_boozer_file(pert_name, pert, status)
        sp(1) = ntv_species_t(1.0_dp, 3.3435837768e-27_dp)
        sp(2) = ntv_species_t(-1.0_dp, 9.1093837e-31_dp)
        prof%s = [0.45_dp, 0.44_dp]
        prof%omega_te = [1e3_dp, 1e3_dp]
        allocate (prof%density(2, 2), prof%temperature(2, 2), prof%density_s(2, 2), &
            prof%temperature_s(2, 2), prof%kappa(2, 2))
        prof%density = 3e19_dp
        prof%density(2, :) = 0
        prof%temperature = 2*kev
        prof%density_s = -1e19_dp
        prof%temperature_s = -1*kev
        prof%kappa = 3e-4_dp
        res%x_tol = 0
        res%n_theta = 17
        res%n_leg = 24
        res%n_x = 16
        res%sources = .false.
        res%collisions = 1
        res%nelem = 4
        call ntv_device_profile(eq, pert, 3, sp, prof, res, both, status)
        call check('failed surface: profile reports the failure', status == DEVICE_SOLVE)
        if (failures > 0) return
        call check('failed surface: its rows carry the failure', &
            all(both(2, :)%status == DEVICE_SOLVE) .and. all(both(2, :)%k11 == 0))
        prof%s = prof%s(1:1)
        prof%omega_te = prof%omega_te(1:1)
        prof%density = prof%density(1:1, :)
        prof%temperature = prof%temperature(1:1, :)
        prof%density_s = prof%density_s(1:1, :)
        prof%temperature_s = prof%temperature_s(1:1, :)
        prof%kappa = prof%kappa(1:1, :)
        call ntv_device_profile(eq, pert, 3, sp, prof, res, alone, status)
        call check('failed surface: the good surface alone', status == DEVICE_OK)
        if (failures > 0) return
        call check('failed surface: the good surface is kept', &
            all(both(1, :)%status == DEVICE_OK) .and. all(both(1, :)%k11 == alone(1, :)%k11))
    end subroutine failed_surface

    !> V_phi/R of the deuterons at R = 1.9 m from the axisymmetric solve at
    !> Omega_tE = 2e3 rad/s; the force balance must return that Omega_tE.
    subroutine force_balance()
        character(*), parameter :: eq_name = 'ntv_force_eq.bc', pert_name = 'ntv_force_pert.bc'
        type(boozer_file_t) :: eq, pert
        type(ntv_species_t) :: sp(2)
        type(ntv_profile_t) :: prof
        type(ntv_resolution_t) :: res
        type(ntv_axisymmetric_result_t), allocatable :: r(:)
        real(dp), parameter :: radius = 1.9_dp, omega_te = 2e3_dp
        real(dp) :: rotation
        integer :: status
        call write_file(eq_name, .false.)
        call write_file(pert_name, .true.)
        call load_boozer_file(eq_name, eq, status)
        call load_boozer_file(pert_name, pert, status)
        sp(1) = ntv_species_t(1.0_dp, 3.3435837768e-27_dp)
        sp(2) = ntv_species_t(-1.0_dp, 9.1093837e-31_dp)
        prof%s = [0.45_dp]
        prof%omega_te = [omega_te]
        allocate (prof%density(1, 2), prof%temperature(1, 2), prof%density_s(1, 2), &
            prof%temperature_s(1, 2), prof%kappa(1, 2))
        prof%density = 3e19_dp
        prof%temperature = 2*kev
        prof%density_s = -1e19_dp
        prof%temperature_s = -1*kev
        prof%kappa(1, :) = [3e-4_dp, 6e-4_dp]
        res%x_tol = 0
        res%n_theta = 17
        res%n_leg = 24
        res%collisions = 3
        res%nelem = 4
        call ntv_device_axisymmetric(eq, pert, 3, sp, prof, res, r, status)
        call check('force balance: axisymmetric solve', status == DEVICE_OK)
        if (failures > 0) return
        rotation = r(1)%omega(1) + r(1)%k_pol(1)/radius**2
        call close_to('force balance: Omega_tE', ntv_omega_te_from_rotation(r(1), sp, prof, &
            1, 1, radius, rotation), omega_te, 1e-10_dp)
    end subroutine force_balance

    !> The file's surfaces are concentric circles r = a sqrt(s) with
    !> B = b0 (1 - 0.1 sqrt(s) cos(theta)): eps = 0.1 sqrt(s), |grad r| = 1,
    !> dpsi_pol/dr = iota (flux/2 pi) 2 r/a^2 and G = mu0 I_pol/(2 pi).
    subroutine indicators()
        character(*), parameter :: eq_name = 'ntv_ind_eq.bc', pert_name = 'ntv_ind_pert.bc'
        real(dp), parameter :: mu0 = 1.25663706212e-6_dp, elementary = 1.602176634e-19_dp
        type(boozer_file_t), target :: eq, pert
        type(ntv_species_t) :: sp(1)
        type(ntv_profile_t) :: prof
        type(ntv_resolution_t) :: res
        type(ntv_indicator_t), allocatable :: ind(:)
        real(dp) :: s, r, iota, dpsi_dr, ds_dr, va, want
        integer :: status
        call write_file(eq_name, .false.)
        call write_file(pert_name, .true.)
        call load_boozer_file(eq_name, eq, status)
        call load_boozer_file(pert_name, pert, status)
        sp(1) = ntv_species_t(1.0_dp, 3.3435837768e-27_dp)
        s = 0.45_dp
        prof%s = [s]
        prof%omega_te = [0.0_dp]
        allocate (prof%density(1, 1), prof%temperature(1, 1), prof%density_s(1, 1), &
            prof%temperature_s(1, 1), prof%kappa(1, 1))
        prof%density = 3e19_dp
        prof%temperature = 2*kev
        prof%density_s = -1e19_dp
        prof%temperature_s = -4*kev
        prof%kappa = 3e-4_dp
        res%x_tol = 0
        res%n_theta = 33
        call ntv_device_indicators(eq, pert, 3, sp, prof, res, ind, status)
        call check('indicators: status', status == DEVICE_OK)
        if (failures > 0) return
        r = a*sqrt(s)
        iota = 0.9_dp - 0.3_dp*s
        ds_dr = 2*r/a**2
        dpsi_dr = abs(iota*flux/(2*pi)*ds_dr)
        va = sqrt(2*prof%temperature(1, 1)/sp(1)%mass)
        call close_to('indicators: eps', ind(1)%eps, eps*sqrt(s), 1e-12_dp)
        call close_to('indicators: gradient length', ind(1)%gradient_length(1), &
            min(3e19_dp/(1e19_dp*ds_dr), 2.0_dp/(4.0_dp*ds_dr)), 1e-6_dp)
        ! <B^2> with the geometric weight R = r0 + r cos(theta) of the circles:
        ! b0^2 (1 + e^2/2 - e r/r0) for B = b0 (1 - e cos(theta)).
        want = sp(1)%mass*va*mu0*abs(jpol)/(2*pi)/(elementary*dpsi_dr*b0 &
            *sqrt(1 + (eps*sqrt(s))**2/2 - eps*sqrt(s)*r/r0))
        call close_to('indicators: orbit excursion', ind(1)%excursion(1), want, 1e-3_dp)
    end subroutine indicators

    !> The default thin kernel, spectral elements on the Boozer line, against
    !> the Legendre kernel at a pitch resolution where it has converged.
    subroutine hp_thin_matches_legendre()
        character(*), parameter :: eq_name = 'ntv_hpthin_eq.bc', pert_name = 'ntv_hpthin_pert.bc'
        type(boozer_file_t) :: eq, pert
        type(ntv_species_t) :: sp(1)
        type(ntv_profile_t) :: prof
        type(ntv_resolution_t) :: res
        type(ntv_species_result_t), allocatable :: hp(:, :), legendre(:, :)
        integer :: status
        call write_file(eq_name, .false.)
        call write_file(pert_name, .true.)
        call load_boozer_file(eq_name, eq, status)
        call load_boozer_file(pert_name, pert, status)
        sp(1) = ntv_species_t(1.0_dp, 3.3435837768e-27_dp)
        prof%s = [0.45_dp]
        prof%omega_te = [1e3_dp]
        allocate (prof%density(1, 1), prof%temperature(1, 1), prof%density_s(1, 1), &
            prof%temperature_s(1, 1), prof%kappa(1, 1))
        prof%density = 3e19_dp
        prof%temperature = 2*kev
        prof%density_s = -1e19_dp
        prof%temperature_s = -1*kev
        prof%kappa = 3e-4_dp
        res%x_tol = 0
        res%n_theta = 33
        res%n_leg = 400
        res%n_x = 1
        res%sources = .false.
        call ntv_device_profile(eq, pert, 3, sp, prof, res, hp, status)
        call check('hp thin kernel', status == DEVICE_OK)
        res%thin_legendre = .true.
        call ntv_device_profile(eq, pert, 3, sp, prof, res, legendre, status)
        call check('Legendre thin kernel', status == DEVICE_OK)
        if (failures > 0) return
        call close_to('hp thin K11', hp(1, 1)%k11, legendre(1, 1)%k11, 2e-3_dp)
        call close_to('hp thin K12', hp(1, 1)%k12, legendre(1, 1)%k12, 2e-3_dp)
        ! With the two NEO-2-QL sources, which move K11 by 1.2 % here.
        res%sources = .true.
        res%n_leg = 800
        res%n_theta = 129
        call ntv_device_profile(eq, pert, 3, sp, prof, res, legendre, status)
        call check('Legendre thin kernel with sources', status == DEVICE_OK)
        res%thin_legendre = .false.
        res%n_theta = 33
        call ntv_device_profile(eq, pert, 3, sp, prof, res, hp, status)
        call check('hp thin kernel with sources', status == DEVICE_OK)
        if (failures > 0) return
        call close_to('hp thin K11, sources', hp(1, 1)%k11, legendre(1, 1)%k11, 5e-4_dp)
        call close_to('hp thin K12, sources', hp(1, 1)%k12, legendre(1, 1)%k12, 5e-4_dp)
    end subroutine hp_thin_matches_legendre

    !> The hp thin kernel on a mesh too coarse for 2e-3 (order 4, 8 columns,
    !> 4 rows) reaches the converged Legendre kernel within 2e-3 when refined
    !> by its dual-weighted residual estimate to hp_adapt = 1e-3.
    subroutine hp_adaptive()
        character(*), parameter :: eq_name = 'ntv_hpadapt_eq.bc', pert_name = 'ntv_hpadapt_pert.bc'
        type(boozer_file_t) :: eq, pert
        type(ntv_species_t) :: sp(1)
        type(ntv_profile_t) :: prof
        type(ntv_resolution_t) :: res
        type(ntv_species_result_t), allocatable :: coarse(:, :), adapted(:, :), legendre(:, :)
        integer :: status
        call write_file(eq_name, .false.)
        call write_file(pert_name, .true.)
        call load_boozer_file(eq_name, eq, status)
        call load_boozer_file(pert_name, pert, status)
        sp(1) = ntv_species_t(1.0_dp, 3.3435837768e-27_dp)
        prof%s = [0.45_dp]
        prof%omega_te = [1e3_dp]
        allocate (prof%density(1, 1), prof%temperature(1, 1), prof%density_s(1, 1), &
            prof%temperature_s(1, 1), prof%kappa(1, 1))
        prof%density = 3e19_dp
        prof%temperature = 2*kev
        prof%density_s = -1e19_dp
        prof%temperature_s = -1*kev
        prof%kappa = 3e-4_dp
        res%x_tol = 0
        res%n_theta = 33
        res%n_leg = 400
        res%n_x = 1
        res%sources = .false.
        res%thin_legendre = .true.
        call ntv_device_profile(eq, pert, 3, sp, prof, res, legendre, status)
        call check('adaptive: Legendre reference', status == DEVICE_OK)
        res%thin_legendre = .false.
        res%hp_columns = 8
        res%hp_rows = 4
        call ntv_device_profile(eq, pert, 3, sp, prof, res, coarse, status)
        call check('adaptive: coarse mesh', status == DEVICE_OK)
        res%hp_adapt = 1e-3_dp
        call ntv_device_profile(eq, pert, 3, sp, prof, res, adapted, status)
        call check('adaptive: refined mesh', status == DEVICE_OK)
        if (failures > 0) return
        call check('adaptive: the coarse mesh misses 2e-3', &
            abs(coarse(1, 1)%k11/legendre(1, 1)%k11 - 1) > 2e-3_dp)
        call close_to('adaptive K11', adapted(1, 1)%k11, legendre(1, 1)%k11, 2e-3_dp)
        call close_to('adaptive K12', adapted(1, 1)%k12, legendre(1, 1)%k12, 2e-3_dp)
    end subroutine hp_adaptive

    !> Circular surfaces with a toroidal current and shear; the perturbation
    !> file carries harmonics m = 1..4 of toroidal number 3.
    !> Deuterons with the temperature and the rotation scaled by 1e-4 (orbit
    !> width and poloidal Mach number by 1e-2, precession and Doppler shift
    !> in the same ratio): the finite-width profile equals the thin one
    !> without the two sources, whose sum the canonical drive replaces, to
    !> first order in the width. The rotation enters the thin kernel as
    !> -Omega_tE and the orbits through their potential, so a sign error
    !> between the two shows here.
    subroutine orbit_width_thin()
        character(*), parameter :: eq_name = 'ntv_orbit_eq.bc', pert_name = 'ntv_orbit_pert.bc'
        type(boozer_file_t), target :: eq, pert
        type(ntv_species_t) :: sp(1)
        type(ntv_profile_t) :: prof
        type(ntv_resolution_t) :: res
        type(ntv_species_result_t), allocatable :: thin(:, :), wide(:, :)
        real(dp), parameter :: factor = 1e-4_dp
        integer :: status, i
        call write_file(eq_name, .false.)
        call write_file(pert_name, .true.)
        call load_boozer_file(eq_name, eq, status)
        call load_boozer_file(pert_name, pert, status)
        sp(1) = ntv_species_t(1.0_dp, 3.3435837768e-27_dp)
        prof%s = [0.45_dp, 0.5_dp]
        prof%omega_te = [3e4_dp, 2e4_dp]*factor
        allocate (prof%density(2, 1), prof%temperature(2, 1), prof%density_s(2, 1), &
            prof%temperature_s(2, 1), prof%kappa(2, 1))
        prof%density = 3e19_dp
        prof%temperature = 2*kev*factor
        prof%density_s = -1e19_dp
        prof%temperature_s = -1*kev*factor
        prof%kappa = 3e-4_dp
        res%x_tol = 0
        res%n_theta = 65
        res%n_leg = 300
        res%n_x = 2
        res%sources = .false.
        call ntv_device_profile(eq, pert, 3, sp, prof, res, thin, status)
        call check('thin profile', status == DEVICE_OK)
        res%orbit_width = .true.
        res%hp_columns = 16
        res%hp_rows = 6
        call ntv_device_profile(eq, pert, 3, sp, prof, res, wide, status)
        call check('finite-width profile', status == DEVICE_OK)
        if (failures > 0) return
        do i = 1, 2
            call close_to('K11 thin orbits', wide(i, 1)%k11, thin(i, 1)%k11, 2e-3_dp)
            call close_to('torque thin orbits', wide(i, 1)%torque, thin(i, 1)%torque, 2e-3_dp)
        end do
    end subroutine orbit_width_thin

    !> Physical deuterons on s = 0.12 of a file that starts at s = 0.1: their
    !> orbits leave the file, which the finite-width model cannot represent.
    !> Those energies carry no NTV and are counted in lost_weight instead of
    !> being evaluated outside the file.
    subroutine orbit_width_lost()
        character(*), parameter :: eq_name = 'ntv_lost_eq.bc', pert_name = 'ntv_lost_pert.bc'
        type(boozer_file_t), target :: eq, pert
        type(ntv_species_t) :: sp(1)
        type(ntv_profile_t) :: prof
        type(ntv_resolution_t) :: res
        type(ntv_species_result_t), allocatable :: wide(:, :)
        integer :: status
        call write_file(eq_name, .false.)
        call write_file(pert_name, .true.)
        call load_boozer_file(eq_name, eq, status)
        call load_boozer_file(pert_name, pert, status)
        sp(1) = ntv_species_t(1.0_dp, 3.3435837768e-27_dp)
        prof%s = [0.12_dp]
        prof%omega_te = [0.0_dp]
        allocate (prof%density(1, 1), prof%temperature(1, 1), prof%density_s(1, 1), &
            prof%temperature_s(1, 1), prof%kappa(1, 1))
        prof%density = 3e19_dp
        prof%temperature = 20*kev
        prof%density_s = -1e19_dp
        prof%temperature_s = -1*kev
        prof%kappa = 3e-4_dp
        res%x_tol = 0
        res%n_theta = 33
        res%n_leg = 80
        res%n_x = 2
        res%sources = .false.
        res%orbit_width = .true.
        call ntv_device_profile(eq, pert, 3, sp, prof, res, wide, status)
        call check('lost energies keep the profile', status == DEVICE_OK)
        call check('lost energies are counted', wide(1, 1)%lost_weight > 0 .and. &
            wide(1, 1)%lost_weight <= 1)
        ! The adaptive energy panels count them the same way.
        res%x_tol = 0.3_dp
        res%hp_columns = 8
        res%hp_rows = 4
        call ntv_device_profile(eq, pert, 3, sp, prof, res, wide, status)
        call check('lost energies keep the adaptive profile', status == DEVICE_OK)
        call check('adaptive panels count lost energies', wide(1, 1)%lost_weight > 0 .and. &
            wide(1, 1)%lost_weight <= 1)
        ! The drive at points whose orbits leave the file is zero, not read
        ! from outside it, and the geometry reports the loss.
        block
            type(orbit_particle_t) :: p
            type(hp_orbit_geometry_t) :: g
            type(hp_operator_t) :: op
            complex(dp) :: src
            real(dp) :: worst
            integer :: j, k
            p%mass = sp(1)%mass
            p%charge = 1.602176634e-19_dp
            p%energy = 20*kev
            p%s_ref = 0.12_dp
            call hp_orbit_geometry(eq, p, 0.12_dp, b0, g, status)
            call hp_orbit_perturbation(g, pert, 3, b0)
            op%n = 3
            worst = 0
            do j = 0, 31
                do k = 0, 31
                    src = g%perturbation(2*pi*j/32, -1 + 2*(k + 0.5_dp)/32, op)
                    worst = max(worst, abs(src))
                end do
            end do
            call check('drive stays finite', worst < huge(1.0_dp))
            call check('lost orbits detected', .not. hp_orbit_confined(g, 32, 32))
        end block
    end subroutine orbit_width_lost

    subroutine write_file(name, perturbation)
        character(*), intent(in) :: name
        logical, intent(in) :: perturbation
        integer :: unit, k, m
        real(dp) :: s
        open (newunit=unit, file=name, status='replace', action='write')
        write (unit, '(a)') 'CC synthetic test file'
        write (unit, '(a)') ' m0b n0b nsurf nper flux/[Tm^2] a/[m] R/[m]'
        write (unit, '(a,i0,a,es24.16,a)') ' ', merge(4, 1, perturbation), ' 0 9 1 ', &
            flux, ' 0.5 1.7'
        do k = 1, 9
            s = 0.1_dp*k
            write (unit, '(a)') ' s iota curr_pol/nper curr_tor pprime sqrt g(0,0)'
            write (unit, '(a)') ' [A] [A] dp/ds,[Pa] (dV/ds)/nper'
            write (unit, '(6es24.16)') s, 0.9_dp - 0.3_dp*s, jpol, itor_s*s, 0.0_dp, 0.0_dp
            write (unit, '(a)') ' m n r/[m] z/[m] (phib-phi)*nper/twopi bmn/[T]'
            if (perturbation) then
                write (unit, '(2i5,4es24.16)') 0, 3, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp
                do m = 1, 4
                    write (unit, '(2i5,4es24.16)') m, 3, 0.0_dp, 0.0_dp, 0.0_dp, &
                        1e-3_dp*b0*s**(0.5_dp*m)/m
                end do
            else
                write (unit, '(2i5,4es24.16)') 0, 0, r0, 0.0_dp, 0.0_dp, b0
                write (unit, '(2i5,4es24.16)') 1, 0, a*sqrt(s), a*sqrt(s), 0.0_dp, &
                    -eps*b0*sqrt(s)
            end if
        end do
        close (unit)
    end subroutine write_file

    subroutine close_to(label, got, want, rtol)
        character(*), intent(in) :: label
        real(dp), intent(in) :: got, want, rtol
        call check(label, abs(got - want) <= rtol*abs(want))
        if (abs(got - want) > rtol*abs(want)) print '(a,2es24.15)', '  got, want ', got, want
    end subroutine close_to

    subroutine check(label, ok)
        character(*), intent(in) :: label
        logical, intent(in) :: ok
        if (ok) return
        failures = failures + 1
        print '(a)', 'FAIL '//label
    end subroutine check

end program test_ntv_device
