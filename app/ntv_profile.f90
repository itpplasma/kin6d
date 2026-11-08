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
!! K_21, K_22, the energy flux [J Wb/(m^3 s)] and the offset rotation [rad/s].
program ntv_profile
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use milca_boozer_file, only: boozer_file_t, load_boozer_file, BOOZER_OK
    use milca_ntv_device, only: ntv_species_t, ntv_profile_t, ntv_resolution_t, &
        ntv_species_result_t, ntv_device_profile, DEVICE_OK
    implicit none
    type(boozer_file_t) :: eq, pert
    type(ntv_species_t), allocatable :: species(:)
    type(ntv_profile_t) :: prof
    type(ntv_resolution_t) :: res
    type(ntv_species_result_t), allocatable :: result(:, :)
    character(1024) :: path, eq_path, pert_path
    integer :: u, n_tor, nsp, ns, i, a, status, sources
    integer(8) :: c0, c1, rate

    if (command_argument_count() /= 1) error stop 'usage: ntv_profile INPUT'
    call get_command_argument(1, path)
    open (newunit=u, file=trim(path), status='old', action='read')
    read (u, '(a)') eq_path
    read (u, '(a)') pert_path
    read (u, *) n_tor, res%n_theta, res%n_leg, res%n_x, res%chunk, res%max_iterations, &
        sources
    res%sources = sources /= 0
    ! Optional line: collision model (0 monoenergetic, else dk_collop model),
    ! x_max, elements, order of the speed grid.
    read (u, '(a)') path
    read (path, *, iostat=status) res%collisions, res%x_max, res%nelem, res%order, res%x_tol
    if (status /= 0) read (path, *, iostat=status) res%collisions, res%x_max, res%nelem, &
        res%order
    if (status == 0) then
        read (u, *) nsp
    else
        res%collisions = 0
        read (path, *) nsp
    end if
    allocate (species(nsp))
    do a = 1, nsp
        read (u, *) species(a)%z, species(a)%mass
    end do
    read (u, *) ns
    allocate (prof%s(ns), prof%omega_te(ns), prof%density(ns, nsp), &
        prof%temperature(ns, nsp), prof%density_s(ns, nsp), prof%temperature_s(ns, nsp), &
        prof%kappa(ns, nsp))
    do i = 1, ns
        read (u, *) prof%s(i), prof%omega_te(i)
        do a = 1, nsp
            read (u, *) prof%density(i, a), prof%temperature(i, a), prof%density_s(i, a), &
                prof%temperature_s(i, a), prof%kappa(i, a)
        end do
    end do
    close (u)

    call load_boozer_file(trim(eq_path), eq, status)
    if (status /= BOOZER_OK) error stop 'cannot read the equilibrium file'
    call load_boozer_file(trim(pert_path), pert, status)
    if (status /= BOOZER_OK) error stop 'cannot read the perturbation file'
    call system_clock(c0, rate)
    call ntv_device_profile(eq, pert, n_tor, species, prof, res, result, status)
    call system_clock(c1)
    if (status /= DEVICE_OK) then
        write (error_unit, '(a,i0)') 'ntv_device_profile failed with status ', status
        error stop 1
    end if
    do i = 1, ns
        do a = 1, nsp
            print '(f12.8,i3,5es18.9,2i8,4es18.9)', prof%s(i), a, result(i, a)%k11, &
                result(i, a)%k12, result(i, a)%dpsi_dr, result(i, a)%gamma_psi, &
                result(i, a)%torque, result(i, a)%factorisations, result(i, a)%iterations, &
                result(i, a)%k21, result(i, a)%k22, result(i, a)%heat_flux, &
                result(i, a)%offset_rotation
        end do
    end do
    write (error_unit, '(a,f10.2)') 'wall s ', real(c1 - c0, dp)/rate
end program ntv_profile
