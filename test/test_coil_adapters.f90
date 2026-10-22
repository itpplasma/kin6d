program test_coil_adapters
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_coils_3d, only: coil_set_t, coil_set_read_makegrid, &
                  coil_set_from_polylines, coil_set_b_cart, coil_w7x_normalized_currents
    implicit none
    type(coil_set_t) :: coils, reference
    real(dp), allocatable :: currents(:)
    integer, allocatable :: groups(:)
    real(dp) :: points(3, 6), direct(3), imported(3), ratios(7), expected(7, 6)
    character(len=3), parameter :: configurations(6) = &
                                   ['EIM', 'DMM', 'FTM', 'SVM', 'WNJ', 'YFK']
    integer :: unit, status, k, variant

    ! Two synthetic triangles, with opposite signed currents and shared group.
    points = reshape([2.0_dp, 0.0_dp, 0.0_dp, 2.0_dp, 1.0_dp, 0.0_dp, &
                      2.0_dp, 0.0_dp, 1.0_dp, -2.0_dp, 0.0_dp, 0.0_dp, &
                      -2.0_dp, 1.0_dp, 0.0_dp, -2.0_dp, 0.0_dp, 1.0_dp], [3, 6])
    call coil_set_from_polylines(reference, [3, 3], points, status=status)
    if (status /= 0) error stop 'reference coils'
    call coil_set_b_cart(reference, [3.0_dp, -7.0_dp], &
                         [0.3_dp, 0.4_dp, 0.6_dp], direct)
    do variant = 0, 7
        open (newunit=unit, file='synthetic_makegrid.dat', status='replace')
        write (unit, '(a)') '  # Synthetic test geometry, metres and A-turns', &
            'periods 5', 'begin filament'
        if (variant == 1) then
            write (unit, '(a)') 'mirror stellarator'
        else
            write (unit, '(a)') 'mirror NUL'
        end if
        do k = 1, 3
            write (unit, *) points(:, k), merge(0.0_dp, 3.0_dp, variant == 7) + &
                merge(1.0_dp, 0.0_dp, variant == 2 .and. k == 2)
        end do
        if (variant == 3) then
            write (unit, *) points(:, 2), 0.0_dp, 2, ' tilted'
        else if (variant == 6) then
            write (unit, '(a)') '2.0 0.0 0.0 NaN 2 tilted'
        else
            write (unit, *) points(:, 1), 0.0_dp, 2, ' tilted'
        end if
        do k = 4, 6
            write (unit, *) points(:, k), -7.0_dp
        end do
        write (unit, *) points(:, 4), 0.0_dp, 2, ' reversed'
        if (variant /= 4) write (unit, '(a)') 'end'
        if (variant == 5) write (unit, '(a)') 'unexpected trailing record'
        close (unit)
        call coil_set_read_makegrid(coils, 'synthetic_makegrid.dat', &
                                    currents, groups, status)
        if (variant == 0 .or. variant == 7) then
            if (status /= 0) error stop 'MAKEGRID read'
            if (coils%ncoil /= 2 .or. coils%nseg /= 6) error stop 'explicit loops'
            if (any(groups /= 2)) error stop 'group identity'
            if (any(abs(currents - &
                [merge(0.0_dp, 3.0_dp, variant == 7), -7.0_dp]) > 0)) &
                error stop 'current signs'
            call coil_set_b_cart(reference, currents, &
                [0.3_dp, 0.4_dp, 0.6_dp], direct)
            if (coils%name(1) /= 'tilted') error stop 'coil name'
            call coil_set_b_cart(coils, currents, [0.3_dp, 0.4_dp, 0.6_dp], imported)
            if (norm2(imported - direct) > 1.e-14_dp*norm2(direct)) &
                error stop 'imported direct field'
        else
            if (status == 0) error stop 'invalid MAKEGRID accepted'
            if (allocated(currents) .or. allocated(groups)) &
                error stop 'failed reader exposed currents'
        end if
    end do
    open (newunit=unit, file='synthetic_makegrid.dat', status='old')
    close (unit, status='delete')
    expected(:, 1) = [1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 0.0_dp, 0.0_dp]
    expected(:, 2) = [1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 0.25_dp, 0.25_dp]
    expected(:, 3) = [1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, -0.23_dp, -0.23_dp]
    expected(:, 4) = [1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 0.6_dp, -0.23_dp, -0.23_dp]
    expected(:, 5) = [1.0_dp, 1.0_dp, 0.5_dp, 1.0_dp, 0.23_dp, -0.08_dp, 0.08_dp]
    expected(:, 6) = [1.0_dp, 1.0_dp, 0.5_dp, 1.0_dp, 0.1_dp, 0.13_dp, 0.25_dp]
    do k = 1, 6
        call coil_w7x_normalized_currents(configurations(k), ratios, status)
        if (status /= 0) error stop 'public configuration'
        if (maxval(abs(ratios - expected(:, k))) > 1.e-15_dp) &
            error stop 'public filament-current ratios'
    end do
    call coil_w7x_normalized_currents('unknown', ratios, status)
    if (status == 0) error stop 'unknown configuration accepted'
    print '(a)', 'PASS: MAKEGRID signed fields, invalid inputs, public currents'
end program
