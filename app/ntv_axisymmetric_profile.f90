!> Axisymmetric transport matrix on the surfaces of an ntv_profile input
!! (milca_ntv_device). Usage: ntv_axisymmetric_profile INPUT, with the input of
!! ntv_profile; the collision line is required (model 1 pitch-angle
!! scattering, 3 full operator), the perturbation file only completes the
!! surface geometry and Omega_tE is not used.
!! Prints per surface and species pair (a, b): s, a, b, d(a, b, i, j) for
!! i = 1..3 (particle and energy coefficients in the units of the NTV K_1j,
!! K_2j, then the parallel flow) and j = 1..3 (A_1, A_2, parallel field of
!! species b), b.grad theta/b, q, dpsi_pol/dr [Wb/m], dpsi_pol/ds [Wb], G [T m],
!! <B^2> [T^2] and the status; then per surface a line 'flows' with s, the
!! parallel flow <u_par B> of every species [T m/s] and the bootstrap current
!! <j_par B> [T A/m^2] for the gradients and Omega_tE of the input.
program ntv_axisymmetric_profile
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use milca_boozer_file, only: boozer_file_t, load_boozer_file, BOOZER_OK
    use milca_ntv_device, only: ntv_species_t, ntv_profile_t, ntv_resolution_t, &
        ntv_axisymmetric_result_t, ntv_device_axisymmetric, ntv_read_input, DEVICE_OK, &
        DEVICE_SOLVE
    implicit none
    type(boozer_file_t) :: eq, pert
    type(ntv_species_t), allocatable :: species(:)
    type(ntv_profile_t) :: prof
    type(ntv_resolution_t) :: res
    type(ntv_axisymmetric_result_t), allocatable :: result(:)
    character(1024) :: path, eq_path, pert_path
    integer :: n_tor, i, a, b, status
    integer(8) :: c0, c1, rate

    if (command_argument_count() /= 1) error stop 'usage: ntv_axisymmetric_profile INPUT'
    call get_command_argument(1, path)
    call ntv_read_input(trim(path), eq_path, pert_path, n_tor, species, prof, res)
    call load_boozer_file(trim(eq_path), eq, status)
    if (status /= BOOZER_OK) error stop 'cannot read the equilibrium file'
    call load_boozer_file(trim(pert_path), pert, status)
    if (status /= BOOZER_OK) error stop 'cannot read the perturbation file'
    call system_clock(c0, rate)
    call ntv_device_axisymmetric(eq, pert, n_tor, species, prof, res, result, status)
    call system_clock(c1)
    if (status /= DEVICE_OK .and. status /= DEVICE_SOLVE) then
        write (error_unit, '(a,i0)') 'ntv_device_axisymmetric failed with status ', status
        error stop 1
    end if
    do i = 1, size(prof%s)
        do a = 1, size(species)
            do b = 1, size(species)
                if (result(i)%status == DEVICE_OK) then
                    print '(f12.8,2i3,15es18.9,i3)', prof%s(i), a, b, &
                        result(i)%d(a, b, :, :), result(i)%kappa, result(i)%q, &
                        result(i)%dpsi_dr, result(i)%psi_pol_s, result(i)%g_cov, &
                        result(i)%b2, result(i)%status
                else
                    print '(f12.8,2i3,15es18.9,i3)', prof%s(i), a, b, spread(0.0_dp, 1, 15), &
                        result(i)%status
                end if
            end do
        end do
        if (result(i)%status == DEVICE_OK) print '(a,f12.8,*(es18.9))', 'flows ', prof%s(i), &
            result(i)%flow, result(i)%bootstrap
    end do
    write (error_unit, '(a,f10.2)') 'wall s ', real(c1 - c0, dp)/rate
    if (status /= DEVICE_OK) error stop 1
end program ntv_axisymmetric_profile
