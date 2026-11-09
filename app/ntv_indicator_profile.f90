!> Orbit-width indicators of the thin-orbit model on the surfaces of an
!! ntv_profile input (milca_ntv_device, ntv_indicator_t); no kinetic solve.
!! Usage: ntv_indicator_profile INPUT. Prints per surface and species: s,
!! species, eps, the gradient length L [m], the passing excursion dr [m] and
!! the banana width [m] at the thermal speed, dr/L and banana width/L.
program ntv_indicator_profile
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use milca_boozer_file, only: boozer_file_t, load_boozer_file, BOOZER_OK
    use milca_ntv_device, only: ntv_species_t, ntv_profile_t, ntv_resolution_t, &
        ntv_indicator_t, ntv_device_indicators, ntv_read_input, DEVICE_OK
    implicit none
    type(boozer_file_t) :: eq, pert
    type(ntv_species_t), allocatable :: species(:)
    type(ntv_profile_t) :: prof
    type(ntv_resolution_t) :: res
    type(ntv_indicator_t), allocatable :: ind(:)
    character(1024) :: path, eq_path, pert_path
    integer :: n_tor, i, a, status

    if (command_argument_count() /= 1) error stop 'usage: ntv_indicator_profile INPUT'
    call get_command_argument(1, path)
    call ntv_read_input(trim(path), eq_path, pert_path, n_tor, species, prof, res)
    call load_boozer_file(trim(eq_path), eq, status)
    if (status /= BOOZER_OK) error stop 'cannot read the equilibrium file'
    call load_boozer_file(trim(pert_path), pert, status)
    if (status /= BOOZER_OK) error stop 'cannot read the perturbation file'
    call ntv_device_indicators(eq, pert, n_tor, species, prof, res, ind, status)
    if (status /= DEVICE_OK) then
        write (error_unit, '(a,i0)') 'ntv_device_indicators failed with status ', status
        error stop 1
    end if
    do i = 1, size(prof%s)
        do a = 1, size(species)
            print '(f12.8,i3,7es16.7)', prof%s(i), a, ind(i)%eps, ind(i)%gradient_length(a), &
                ind(i)%excursion(a), ind(i)%banana(a), &
                ind(i)%excursion(a)/ind(i)%gradient_length(a), &
                ind(i)%banana(a)/ind(i)%gradient_length(a)
        end do
    end do
end program ntv_indicator_profile
