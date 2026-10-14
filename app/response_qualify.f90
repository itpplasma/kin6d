program response_qualify
    !! One frozen root, 72 cases x two phases x two target rationals.
    !! Optional adapters only serialize inputs; all gates are owned here.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_result_io, only: run_t, result_file_t
    use kin6d_response_qualification, only: response_record_t, qualification_t, &
        qualify_response
    implicit none
    integer, parameter :: maxrows = 288
    integer :: nrows = 0, cut(maxrows) = 0, phase(maxrows) = 0, target(maxrows) = 0
    integer :: axis(maxrows) = 0, level(maxrows) = 0
    complex(dp) :: b(maxrows) = 0, b_vac(maxrows) = 0, total_flux(maxrows) = 0
    complex(dp) :: current(maxrows) = 0, gpec(2, 2, 3) = 0
    real(dp) :: residual(maxrows) = huge(1.0_dp), balance(maxrows) = huge(1.0_dp)
    real(dp) :: extent(2, maxrows) = 0, mars_amplitude(2, 2, 3) = 0
    logical :: available(maxrows) = .false., coupled(maxrows) = .false.
    logical :: calibrated(maxrows) = .false., full_pressure(maxrows) = .false.
    logical :: all_rationals(maxrows) = .false., interface_test(maxrows) = .false.
    logical :: whole_edge(maxrows) = .false., mars_admitted(2, 2, 3) = .false.
    namelist /qualification/ nrows, cut, phase, target, axis, level, b, b_vac, &
        total_flux, current, gpec, residual, balance, extent, mars_amplitude, &
        available, coupled, calibrated, full_pressure, all_rationals, &
        interface_test, whole_edge, mars_admitted
    type(response_record_t) :: rows(3, 8, 2, 2, 3)
    type(qualification_t) :: q(2, 2, 3)
    logical :: seen(3, 8, 2, 2, 3)
    type(run_t) :: run
    type(result_file_t) :: output
    integer :: u, ios, j, c, p, t, a, l, failures
    character(512) :: message
    real(dp) :: gates(8, 12), changes(2, 8, 12), ratios(3, 12), passed(12), verdict(12)
    call run%start('response_qualify')
    if (.not. run%has_config()) error stop 'response_qualify --config qualification.nml'
    call run%open_config(u)
    read (u, nml=qualification, iostat=ios, iomsg=message)
    call run%close_config(u, 'qualification', ios, message)
    if (ios /= 0) error stop 'qualification group required'
    if (nrows < 0 .or. nrows > maxrows) error stop 'invalid row count'
    seen = .false.
    do j = 1, nrows
        c = cut(j); p = phase(j); t = target(j); a = axis(j); l = level(j)
        if (c < 1 .or. c > 3 .or. p < 1 .or. p > 2 .or. t < 1 .or. t > 2 &
            .or. a < 1 .or. a > 8 .or. l < 1 .or. l > 3) error stop 'row identity'
        if (seen(l, a, t, p, c)) error stop 'duplicate row'
        seen(l, a, t, p, c) = .true.
        rows(l, a, t, p, c) = response_record_t(b=b(j), b_vac=b_vac(j), &
            total_flux=total_flux(j), current=current(j), residual=residual(j), &
            balance=balance(j), extent=extent(:, j), available=available(j), &
         coupled=coupled(j), calibrated=calibrated(j), full_pressure=full_pressure(j), &
            all_rationals=all_rationals(j), interface_test=interface_test(j), &
            whole_edge=whole_edge(j))
    end do
    call run%open_resolved(u)
    write (u, nml=qualification)
    call run%close_resolved(u)
    failures = 0; j = 0
    do c = 1, 3
        do p = 1, 2
            do t = 1, 2
                j = j + 1
                call qualify_response(rows(:, :, t, p, c), gpec(t, p, c), &
                    mars_amplitude(t, p, c), mars_admitted(t, p, c), q(t, p, c))
                if (.not. q(t, p, c)%passed) failures = failures + 1
                gates(:, j) = merge(1.0_dp, 0.0_dp, q(t, p, c)%gates)
                changes(:, :, j) = q(t, p, c)%changes
                ratios(:, j) = [q(t, p, c)%gpec_ratio, q(t, p, c)%gpec_phase, &
                    q(t, p, c)%mars_ratio]
                passed(j) = merge(1.0_dp, 0.0_dp, q(t, p, c)%passed)
                verdict(j) = q(t, p, c)%verdict
            end do
        end do
    end do
    if (failures > 6) then
        verdict = 0
    end if
    print '(a)', 'cut phase target gates_passed gpec_ratio phase_rad verdict'
    j = 0
    do c = 1, 3
        do p = 1, 2
            do t = 1, 2
                j = j + 1
                print '(3i4,i5,2es14.5,i4)', c, p, t, count(q(t, p, c)%gates), &
                    ratios(1:2, j), nint(verdict(j))
            end do
        end do
    end do
    if (failures > 6) then
        print '(a,i0,a)', 'STOP_NO_TIEBREAK: ', failures, '/12 rows fail gates'
    else
        print '(a,i0,a)', 'QUALIFICATION: ', 12 - failures, '/12 rows pass gates'
    end if
    call execute_command_line('mkdir -p output', exitstat=ios)
    if (ios /= 0) error stop 'output directory'
    call output%create(run)
    call output%attribute('claim', &
        'candidate; numerical/model/input bounds unavailable')
    call output%attribute('row_order', 'cut c0995,c0997,full; phase 000,010; q 3,10/3')
    call output%attribute('axes', &
        'radial core grading harmonic poloidal extent S gamma')
    call output%write('gates', gates, '1', 'empirical qualification gates', &
        'axis,row', 'measurement')
    call output%write('changes', changes, '1', 'complex ladder relative changes', &
        'change,axis,row', 'measurement')
    call output%write('comparison', ratios, '1', &
        'GPEC ratio, phase radians, MARS ratio', &
        'diagnostic,row', 'measurement')
    call output%write('passed', passed, '1', 'all gates passed', 'row', 'measurement')
    call output%write('verdict', verdict, '1', '0 none, 1 GPEC, 2 MARS', &
        'row', 'measurement')
    call output%close()
end program response_qualify
