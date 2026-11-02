program test_island_vacuum_oracle
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use nn_vacuum_fixture, only: make_coil, coils, current, r0, ftor, field, &
                                 direct_coil, direct_field, solve_vacuum, orbit
    use kin6d_coil_quadrature_field, only: coil_quadrature_field_t, &
                                           coil_quadrature_init
    use kin6d_equilibrium_state, only: EQ_VALUES, EQ_BR, EQ_BZ, EQ_BPHI
    use kin6d_island_diagnostics, only: island_trace_options_t, island_point_t, &
                    island_width_t, island_periodic_point, island_width, island_trace, &
                                        island_stability
    use fortnum_multiroot, only: multiroot_hybrids
    use fortnum_status, only: fortnum_status_t
    implicit none
    type(coil_quadrature_field_t) :: quad, fine
    type(island_trace_options_t) :: opt, reference_opt
    type(island_point_t) :: op, xp, ref_o, ref_x
    type(island_width_t) :: width, ref_width
    type(fortnum_status_t) :: stat
    real(dp) :: points(3, 64), values(EQ_VALUES, 64), finer(EQ_VALUES, 64)
    real(dp) :: exact(3), err, quad_err, phase_error(4), width_error(4), rate, phase_rate
    real(dp) :: map512(2), map1024(2), span, root(2), residual(2), jc(2, 2), jf(2, 2)
    real(dp) :: closure, floor_error
    integer :: st, statuses(64), j, k, base, lev, levels, n
    logical :: bad_points
    character(32) :: arg, mode
    call get_command_argument(1, mode)
    call make_coil()
    call coil_quadrature_init(quad, coils, [current], 8, ftor, st)
    if (st /= 0) error stop 'quadrature field'
    call coil_quadrature_init(fine, coils, [current], 12, ftor, st)
    if (st /= 0) error stop 'fine quadrature field'
    do j = 1, 64
        points(:, j) = [r0 + 0.4_dp*cos(j*0.73_dp), &
                        0.4_dp*sin(j*0.73_dp), j*0.37_dp]
    end do
    call quad%sample(points, values, statuses)
    if (any(statuses /= 0)) error stop 'quadrature samples'
    call fine%sample(points, finer, statuses)
    if (any(statuses /= 0)) error stop 'fine samples'
    err = 0
    do j = 1, 64
        call direct_field(points(:2, j), points(3, j), exact)
        err = max(err, norm2([values(EQ_BR, j), values(EQ_BPHI, j), &
                              values(EQ_BZ, j)] - exact))
    end do
    quad_err = maxval(abs(values - finer))
    print '(a,2es14.6)', 'Biot-Savart quadrature vs segments / 8 vs 12 ', err, quad_err
    if (max(err, quad_err) > 1.0e-8_dp) error stop 'Biot-Savart quadrature'
    if (trim(mode) == 'quadrature') stop
    call get_command_argument(2, arg); read (arg, *) k
    call get_command_argument(3, arg); read (arg, *) base
    if (trim(mode) == 'seed') then
        call solve_vacuum(base, k, 3)
        print '(a)', 'PASS: high-order vacuum seed meets energy and residual gates'
        stop
    end if
    if (trim(mode) == 'coarse') then
        call solve_vacuum(base, k, 3)
        ! Reproduce the old root algorithm against the actual FEEC consumer.
        call multiroot_hybrids(old_residual, 2, [5.52_dp, 0.0_dp], root, stat, &
                               xtol=1.0e-8_dp, ftol=2.0e-6_dp, max_iter=30)
        call old_residual(root, residual)
        print '(a,i0,a,3es14.6)', 'legacy status ', stat%code, ' root/closure ', &
            root, norm2(residual)
        call orbit(root, 1.0_dp, 512, .false., map512, span)
        call orbit(root, 1.0_dp, 1024, .false., map1024, span)
        floor_error = norm2(map512 - map1024)
        opt%fixed_steps = 512; opt%map_delta = 6.055454452393343e-6_dp*r0
        call island_stability(field, root, 0.0_dp, 1.0_dp, opt, jc, closure, st)
        opt%fixed_steps = 1024
        call island_stability(field, root, 0.0_dp, 1.0_dp, opt, jf, closure, st)
        print '(a,2es14.6)', 'map step-doubling / Jacobian schedule change ', &
            floor_error, maxval(abs(jc - jf))
        opt = island_trace_options_t()
        opt%rtol = 1.0e-9_dp; opt%atol = 1.0e-11_dp
        opt%root_tol = 2.0e-6_dp
        call measure_points()
        if (op%status /= 0 .or. xp%status /= 0) error stop 'unresolved FEEC O/X'
        print '(a)', 'PASS: coarse FEEC map roots with controlled tracing'
        stop
    end if
    call get_command_argument(4, arg); read (arg, *) levels
    if (levels < 2 .or. levels > 4) error stop 'ladder levels'
    opt%fixed_steps = 4096; opt%root_tol = 2.0e-5_dp
    opt%rtol = 1.0e-10_dp; opt%atol = 1.0e-12_dp
    opt%max_step = 2*acos(-1.0_dp)/256; opt%max_steps = 1000000
    reference_opt = opt; reference_opt%fixed_steps = 0
    call island_periodic_point(quad, [5.52_dp, 0.0_dp], 0.0_dp, 1, &
                               [r0, 0.0_dp], 1, -1, reference_opt, ref_o)
    call island_periodic_point(quad, [6.7_dp, 0.0_dp], 0.0_dp, 1, &
                               [r0, 0.0_dp], 1, -1, reference_opt, ref_x)
    if (ref_o%status /= 0 .or. ref_x%status /= 0) error stop 'oracle O/X'
    call island_width(direct_coil, [0.36_dp, 0.63_dp], 0.48_dp, [r0, 0.0_dp], &
                      [-1.0_dp, 0.0_dp], 0.0_dp, 32.0_dp, 1, -1, opt, radius, &
                      'u = distance from (6,0) [m]', ref_width, scans=16)
    if (ref_width%status /= 0) error stop 'oracle width'
    print '(a,5es14.6)', 'oracle width / trace / duration / scan / phase ', &
        ref_width%metres, ref_width%trace_error, ref_width%duration_error, &
        ref_width%scan_error, ref_o%phase
    bad_points = .false.
    do lev = 1, levels
        n = base*2**(lev - 1)
        call solve_vacuum(n, k, 3)
        call measure_points()
        if (op%status /= 0 .or. xp%status /= 0) bad_points = .true.
        call island_width(field, [0.36_dp, 0.63_dp], 0.48_dp, [r0, 0.0_dp], &
                          [-1.0_dp, 0.0_dp], 0.0_dp, 32.0_dp, 1, -1, opt, radius, &
                          'u = distance from (6,0) [m]', width, scans=16)
        if (width%status /= 0) then
            print '(a,3i10)', 'width failure degree/nr/status ', k, n, width%status
            error stop 'FEEC width'
        end if
        width_error(lev) = abs(width%metres - ref_width%metres)
        phase_error(lev) = abs(atan2(sin(op%phase - ref_o%phase), &
                                     cos(op%phase - ref_o%phase)))
        print '(a,i0,a,i0,a,7es14.6)', 'P', k, ' nr=', n, &
            ' h / width / error / trace / duration / phase error / phase indicator ', &
            0.9_dp/n, width%metres, width_error(lev), width%trace_error, &
            width%duration_error, phase_error(lev), op%phase_error
    end do
    rate = log(max(width_error(1), tiny(1.0_dp)) &
               /max(width_error(levels), tiny(1.0_dp)))/log(2.0_dp**(levels - 1))
    print '(a,f10.5)', 'finite-time island-width error order ', rate
    phase_rate = log(max(phase_error(1), tiny(1.0_dp)) &
        /max(phase_error(levels), tiny(1.0_dp)))/log(2.0_dp**(levels - 1))
    print '(a,f10.5)', 'helical-phase error order ', phase_rate
    if (bad_points) error stop 'FEEC O/X unresolved on ladder'
    if (rate <= 0) error stop 'width does not converge'
    if (phase_error(levels) > phase_error(1) + op%phase_error + ref_o%phase_error) &
        error stop 'phase does not converge within tracing indicator'
    print '(a)', 'PASS: finite-time vacuum width/phase comparison; enclosure unavailable'
contains
    function radius(x) result(u)
        real(dp), intent(in) :: x(2)
        real(dp) :: u
        u = norm2(x - [r0, 0.0_dp])
    end function
    subroutine old_residual(x, f, ctx)
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: f(:)
        class(*), intent(in), optional :: ctx
        real(dp) :: finish(2), range
        call orbit(x, 1.0_dp, 512, .false., finish, range)
        f = finish - x
    end subroutine
    subroutine measure_points()
        type(island_trace_options_t) :: point_opt
        point_opt = opt
        point_opt%fixed_steps = 0
        point_opt%rtol = min(point_opt%rtol, 1.0e-11_dp)
        point_opt%atol = min(point_opt%atol, 1.0e-13_dp)
        point_opt%max_step = min(point_opt%max_step, 2*acos(-1.0_dp)/512)
        point_opt%map_delta = 1.0e-4_dp
        point_opt%root_tol = 2.0e-6_dp
        call island_periodic_point(field, [5.52_dp, 0.0_dp], 0.0_dp, 1, &
                                 [r0, 0.0_dp], 1, -1, point_opt, op, [5.49_dp, 5.54_dp])
        call island_periodic_point(field, [6.7_dp, 0.0_dp], 0.0_dp, 1, &
                                 [r0, 0.0_dp], 1, -1, point_opt, xp, [6.69_dp, 6.72_dp])
        print '(a,2i3,a,4es13.5,a,2l2)', 'O/X statuses ', op%status, xp%status, &
            ' closure/residue ', op%closure, xp%closure, op%residue, xp%residue, &
            ' symmetry ', op%used_symmetry, xp%used_symmetry
        if (op%status == 0) then
            if (op%residue <= 0 .or. op%residue >= 1) error stop 'O classification'
        end if
        if (xp%status == 0) then
            if (xp%residue >= 0) error stop 'X classification'
        end if
    end subroutine
end program
