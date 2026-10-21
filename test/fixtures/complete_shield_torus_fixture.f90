program complete_shield_torus_fixture
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use geqdsk_tools, only: geqdsk_t, geqdsk_write
    use kin6d_vacuum_surfaces, only: psi_grid_t, eqdsk_read, psi_eval, &
                              profiles_of_psin, surface_t, surface_build, flux_ray_point
    use kin6d_grad_shafranov, only: gs_axis_safety_factor
    implicit none
    type(psi_grid_t) :: grid
    type(geqdsk_t) :: g
    type(surface_t) :: surf
    real(dp) :: theta, rho, f, df, p, pp, value, pr, pz, prr, prz, pzz
    integer :: i, status, original_unit, generated_unit, output_unit, ios
    logical :: ok
    character(1024) :: path, output
    character(80) :: original_record, generated_record
    call get_command_argument(1, path); call get_command_argument(2, output)
    call eqdsk_read(trim(path), grid, status)
    if (status /= 0) error stop 'legacy prefix read'
    g%nw = grid%nr; g%nh = grid%nz; g%nbbbs = 128; g%limitr = 0
    g%header = 'CHEASE A=10 q0=1.55; circular shielding test'
    ! Only the writer-produced trailer is retained. All preceding records in
    ! the committed fixture are preserved, including original pp/ff profiles.
    g%rdim = 1; g%zdim = 1; g%rcentr = 1; g%rleft = 1; g%zmid = 0
    g%rmaxis = 1; g%zmaxis = 0; g%simag = 0; g%sibry = 1; g%bcentr = 1; g%current = 1
    allocate (g%fpol(65), g%pres(65), g%ffprim(65), g%pprime(65), &
              g%qpsi(65), g%psirz(65, 65), g%rbbbs(128), g%zbbbs(128), &
              g%rlim(0), g%zlim(0))
    g%fpol = 1; g%pres = 0; g%ffprim = 0; g%pprime = 0; g%psirz = 0
    call profiles_of_psin(grid, 0.0_dp, f, df, p, pp)
    call psi_eval(grid, grid%r_axis, grid%z_axis, value, pr, pz, prr, prz, pzz)
    call gs_axis_safety_factor(grid%r_axis, f, &
                               reshape([prr, prz, prz, pzz], [2, 2]), g%qpsi(1), status)
    if (status /= 0) error stop 'axis q'
    do i = 2, 65
        call surface_build(grid, real(i - 1, dp)/64, 512, surf, status)
        if (status /= 0) error stop 'surface q'
        g%qpsi(i) = -surf%q ! geometric clockwise field -> positive COCOS 1 q
    end do
    do i = 1, 128
        theta = 2*acos(-1.0_dp)*(i - 1)/128
        call flux_ray_point(grid, 1.0_dp, theta, rho, ok)
        if (.not. ok) error stop 'LCFS ray'
        g%rbbbs(i) = (grid%r_axis + rho*cos(theta))*100
        g%zbbbs(i) = (grid%z_axis + rho*sin(theta))*100
    end do
    print *, 'reconstructed source q axis / edge', g%qpsi(1), g%qpsi(65)
    if (trim(path) == trim(output)) error stop 'use a separate output file'
    call geqdsk_write(g, trim(output)//'.libneo')
    open (newunit=original_unit, file=trim(path), status='old', action='read')
    open (newunit=generated_unit, file=trim(output)//'.libneo', status='old')
    open (newunit=output_unit, file=trim(output), status='replace', action='write')
    ! Keep the standard libneo header and the original 901 numeric records.
    do i = 1, 902
        read (original_unit, '(a)') original_record
        read (generated_unit, '(a)') generated_record
        if (i == 1) original_record = generated_record
        write (output_unit, '(a)') trim(original_record)
    end do
    do
        read (generated_unit, '(a)', iostat=ios) generated_record
        if (ios /= 0) exit
        write (output_unit, '(a)') trim(generated_record)
    end do
    close (original_unit); close (generated_unit, status='delete'); close (output_unit)
end program
