!> Quasilinear NTV torque profile on a Boozer equilibrium (milca_ntv_device).
!! Usage: ntv_profile INPUT. INPUT, SI units, free format:
!!   equilibrium .bc path
!!   perturbation .bc path
!!   n_tor n_theta n_leg n_x chunk max_iterations sources(0|1)
!!   n_species, then per species: Z mass[kg]
!!   n_surfaces, then per surface: s Omega_tE[rad/s], followed by one line per
!!   species: density[m^-3] T[J] dn/ds dT/ds kappa[1/m]
!! Prints per surface and species: s, species, K_11, K_12, dpsi_pol/dr,
!! Gamma_psi, torque density [N m/m^3], factorisations, GMRES iterations,
!! K_21, K_22, the energy flux [J Wb/(m^3 s)], the offset rotation [rad/s], the
!! share of the K_11 weight on energies whose finite-width orbits leave the file
!! (zero without orbit_width) and the status (0 solved, 2 solve failed, 3 GMRES
!! did not converge, 4 orbits leave the file with adaptive energy panels; the
!! row is then zero). The program ends with an error when any status is not 0.
program ntv_profile
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use milca_boozer_file, only: boozer_file_t, load_boozer_file, BOOZER_OK
    use milca_ntv_device, only: ntv_species_t, ntv_profile_t, ntv_resolution_t, &
        ntv_species_result_t, ntv_device_profile, ntv_read_input, DEVICE_OK, DEVICE_SOLVE
    implicit none
    type(boozer_file_t), target :: eq, pert
    type(ntv_species_t), allocatable :: species(:)
    type(ntv_profile_t) :: prof
    type(ntv_resolution_t) :: res
    type(ntv_species_result_t), allocatable :: result(:, :)
    character(1024) :: path, eq_path, pert_path
    integer :: n_tor, nsp, ns, i, a, status
    integer(8) :: c0, c1, rate

    if (command_argument_count() /= 1) error stop 'usage: ntv_profile INPUT'
    call get_command_argument(1, path)
    call ntv_read_input(trim(path), eq_path, pert_path, n_tor, species, prof, res)
    nsp = size(species)
    ns = size(prof%s)

    call load_boozer_file(trim(eq_path), eq, status)
    if (status /= BOOZER_OK) error stop 'cannot read the equilibrium file'
    call load_boozer_file(trim(pert_path), pert, status)
    if (status /= BOOZER_OK) error stop 'cannot read the perturbation file'
    call system_clock(c0, rate)
    call ntv_device_profile(eq, pert, n_tor, species, prof, res, result, status)
    call system_clock(c1)
    if (status /= DEVICE_OK .and. status /= DEVICE_SOLVE) then
        write (error_unit, '(a,i0)') 'ntv_device_profile failed with status ', status
        error stop 1
    end if
    do i = 1, ns
        do a = 1, nsp
            print '(f12.8,i3,5es18.9,2i8,5es18.9,i3)', prof%s(i), a, result(i, a)%k11, &
                result(i, a)%k12, result(i, a)%dpsi_dr, result(i, a)%gamma_psi, &
                result(i, a)%torque, result(i, a)%factorisations, result(i, a)%iterations, &
                result(i, a)%k21, result(i, a)%k22, result(i, a)%heat_flux, &
                result(i, a)%offset_rotation, result(i, a)%lost_weight, result(i, a)%status
        end do
    end do
    write (error_unit, '(a,f10.2)') 'wall s ', real(c1 - c0, dp)/rate
    if (status /= DEVICE_OK) then
        write (error_unit, '(a,i0,a)') 'ntv_device_profile: ', count(result%status /= DEVICE_OK), &
            ' surface-species pairs failed; see the status column'
        error stop 1
    end if
end program ntv_profile
