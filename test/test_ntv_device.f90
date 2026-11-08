!> NTV torque profile on Boozer files (kin6d_ntv_device): the deflection
!! frequency in its Lorentz limit, and sweeps over neighbouring surfaces
!! preconditioned by one factorisation agree with a factorisation per surface,
!! and a second solve at the reported offset rotation carries no flux.
program test_ntv_device
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_boozer_file, only: boozer_file_t, load_boozer_file, BOOZER_OK
    use kin6d_ntv_device, only: ntv_species_t, ntv_profile_t, ntv_resolution_t, &
        ntv_species_result_t, ntv_device_profile, deflection_frequency, DEVICE_OK
    implicit none
    real(dp), parameter :: r0 = 1.7_dp, a = 0.5_dp, b0 = 2.3_dp, eps = 0.1_dp
    real(dp), parameter :: jpol = -1.2e7_dp, itor_s = -3.0e5_dp, flux = -1.3_dp
    real(dp), parameter :: pi = 4*atan(1.0_dp), kev = 1.602176634e-16_dp
    integer :: failures = 0

    call lorentz_limit()
    call preconditioned_sweep()
    call energy_coupled()
    call offset_rotation()
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
        res%n_theta = 33
        res%n_leg = 60
        res%n_x = 6
        res%chunk = 1
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
        res%n_theta = 17
        res%n_leg = 24
        res%n_x = 16
        res%sources = .false.
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
        type(boozer_file_t) :: eq, pert
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

    !> Circular surfaces with a toroidal current and shear; the perturbation
    !> file carries harmonics m = 1..4 of toroidal number 3.
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
