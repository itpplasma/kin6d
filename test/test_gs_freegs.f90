program test_gs_freegs
    !! Public end-to-end comparison with FreeGS 01-freeboundary.py v0.8.2.
    !! The reference fixes currents and the (1-s)^2 profile shape. Native
    !! solves impose the published pressure/current constraints as well.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_axisymmetric_gs_kernel, only: gs_mu0
    use fortfem_triangle_duffy_quadrature, only: triangle_duffy_quadrature
    use fortnum_polynomial, only: lagrange_weights
    use fortnum_multiroot, only: multiroot_hybrids
    use fortnum_status, only: fortnum_status_t, status_ok
    use kin6d_grad_shafranov, only: gs_state_t, initialize_gs, gs_cell_data, sample_gs
    use kin6d_gs_coils, only: gs_circuit_set_t, gs_read_coil_table
    use kin6d_gs_profile_law, only: gs_profile_law_t, gs_law_primitives
    use kin6d_gs_free_boundary, only: gs_fb_config_t, gs_fb_system_t, &
                                      gs_fb_eval_t, gs_fb_history_t, &
                                      gs_fb_setup, gs_fb_newton, &
                                      gs_fb_evaluate, gs_fb_lcfs_levels, &
                                      gs_fb_trace_lcfs, GS_FB_XPOINT
    use kin6d_gs_equilibrium, only: gs_equilibrium
    use kin6d_equilibrium_state, only: equilibrium_state_t, labelled_field_t
    implicit none

    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp), parameter :: levels(4) = [0.01_dp, 0.5_dp, 0.9_dp, 0.95_dp]
    character(20), parameter :: names(16) = [character(20) :: &
                                             'psi_axis', 'psi_boundary', &
                                             'axis_R', 'axis_Z', &
                                             'Xlower_R', 'Xlower_Z', &
                                             'Xupper_R', 'Xupper_Z', &
                                             'q(.01), signed', 'q(.5), signed', &
                                             'q(.9), signed', 'q(.95), signed', &
                                             'beta_p', 'li2', 'Ip', 'volume']
    type(gs_state_t) :: state
    type(gs_circuit_set_t) :: coils
    type(gs_profile_law_t) :: law
    type(gs_fb_system_t) :: sys
    type(gs_fb_config_t) :: cfg
    type(gs_fb_eval_t) :: ev
    type(gs_fb_history_t) :: history
    type(equilibrium_state_t) :: eq
    real(dp) :: reference(23, 3), seed(33, 25), currents(4), values(16)
    real(dp) :: pa, pb, axis(2), xp(2), f, pressure, coarea(5), flux_levels(5)
    real(dp), allocatable :: contour(:, :), ga(:), gb(:), chain(:, :), xpoints(:, :)
    real(dp) :: p_integral, b_integral, volume, r_geo, paxis
    real(dp) :: pt(2), br, bz, bt, psi_probe, gap, max_gap, seed_error
    real(dp) :: candidates(16, 2), pressures(2), gaps(2), expected(16), budgets(16)
    real(dp) :: grid_error(16), q_error(4), q_rows(7, 4), gap_row(2), ref_gap
    real(dp) :: scales(16), caps(16), pressure_budget, gap_budget
    real(dp) :: critical_row(9), volume_rows(4, 3), volume_error(3)
    real(dp) :: trace_distances(2), trace_gaps(2), distance, trace_gap
    integer :: unit, i, j, n, status, cell, nside, rung, failures
    character(1024) :: path, argument, line

    call get_command_argument(1, path)
    nside = 24
    if (command_argument_count() > 1) then
        call get_command_argument(2, argument)
        read (argument, *) nside
    end if
    open (newunit=unit, file=trim(path)//'/reference.dat', status='old')
    do i = 1, 3
        call data_line(unit, line)
        read (line, *) reference(:, i)
    end do
    close (unit)
    open (newunit=unit, file=trim(path)//'/critical-refinement.dat', status='old')
    do i = 1, 3
        call data_line(unit, line)
        read (line, *) critical_row
    end do
    close (unit)
    open (newunit=unit, file=trim(path)//'/volume-convergence.dat', status='old')
    do i = 1, 3
        call data_line(unit, line)
        read (line, *) volume_rows(:, i)
    end do
    close (unit)
    volume_error = max(abs(volume_rows(2:4, 3) - volume_rows(2:4, 2)), &
                       abs(volume_rows(2:4, 2) - volume_rows(2:4, 1)))
    open (newunit=unit, file=trim(path)//'/seed.dat', status='old')
    do i = 1, 25
        call data_line(unit, line)
        read (line, *) seed(:, i)
    end do
    close (unit)
    call gs_read_coil_table(trim(path)//'/coils.dat', coils, status, line)
    if (status /= 0) print *, trim(line)
    call check(status, 'coil table')
    failures = 0
    q_error = 0
    open (newunit=unit, file=trim(path)//'/q-convergence.dat', status='old')
    do i = 1, 3
        do j = 1, 4
            call data_line(unit, line)
            read (line, *) q_rows(:, j)
        end do
        if (i /= 3) cycle
        q_error = max(abs(q_rows(4:7, 4) - q_rows(4:7, 3)), &
                      abs(q_rows(4:7, 3) - q_rows(4:7, 2)))
    end do
    close (unit)
    ref_gap = 0
    open (newunit=unit, file=trim(path)//'/lcfs-convergence.dat', status='old')
    do i = 1, 3
        call data_line(unit, line)
        read (line, *) gap_row
        if (i > 1) ref_gap = max(ref_gap, gap_row(2))
    end do
    close (unit)
    currents = reference(2:5, 3)
    allocate (law%knots(2), law%pp(4, 1), law%ff(4, 1))
    law%normalized = .true.
    law%knots = [0.0_dp, 1.0_dp]
    law%axis = reference(8, 3) - reference(9, 3)
    law%f_edge = 2.0_dp
    law%pp(:, 1) = reference(6, 3)*[1.0_dp, -2.0_dp, 1.0_dp, 0.0_dp]
    law%ff(:, 1) = reference(7, 3)*[1.0_dp, -2.0_dp, 1.0_dp, 0.0_dp]

    do rung = 1, 2
        call run_case(nside*2**(rung - 1))
        candidates(:, rung) = values
        pressures(rung) = paxis
        gaps(rung) = max_gap
    end do
    expected = reference(8:23, 3)
    expected(9:12) = -expected(9:12)
    grid_error = max(abs(reference(8:23, 3) - reference(8:23, 2)), &
                     abs(reference(8:23, 2) - reference(8:23, 1)))
    grid_error(9:12) = grid_error(9:12) + q_error
    grid_error(:8) = max(grid_error(:8), abs(critical_row(2:9)))
    grid_error(13:14) = grid_error(13:14) + volume_error(:2)
    grid_error(16) = grid_error(16) + volume_error(3)
    scales = max(abs(expected), 1.0e-6_dp)
    scales(3:8) = 1.0_dp
    ! Empirical refinement indicators, never continuous enclosures. No assumed
    ! order: retain the larger FreeGS step because critical/mask outputs can
    ! stagnate or oscillate. A physical cap rejects an unresolved comparison.
    budgets = 2*(grid_error + abs(candidates(:, 2) - candidates(:, 1))) &
              + 1.0e-8_dp*scales
    caps = 0.01_dp*scales
    caps(3:8) = 0.005_dp
    do i = 1, 16
        print '(a20,4es16.7)', names(i), candidates(i, 2), expected(i), &
            abs(candidates(i, 2) - expected(i)), budgets(i)
        call require(budgets(i) <= caps(i), trim(names(i))//' budget resolved')
        call require(abs(candidates(i, 2) - expected(i)) <= budgets(i), &
                     trim(names(i))//' cross-code agreement')
    end do
    pressure_budget = 2*abs(pressures(2) - pressures(1)) + 1.0e-5_dp
    call require(pressure_budget <= 10.0_dp, 'axis pressure budget resolved')
    call require(abs(pressures(2) - 1000.0_dp) <= pressure_budget, &
                 'same 1 kPa axis pressure constraint')
    gap_budget = 2*(ref_gap + abs(critical_row(3)) &
                    + abs(gaps(2) - gaps(1))) + 1.0e-10_dp
    call require(gap_budget <= 0.01_dp*abs(law%axis), 'LCFS budget resolved')
    call require(gaps(2) <= gap_budget, 'LCFS flux gaps at FreeGS samples')
    print '(a,3es16.7)', 'paxis: candidate, discrepancy, budget ', &
        pressures(2), abs(pressures(2) - 1000.0_dp), pressure_budget
    print '(a,2es16.7)', 'LCFS: flux gap, budget (Wb/rad) ', gaps(2), gap_budget
    call require(trace_gaps(2) <= gap_budget, &
                 'native traced LCFS at unchanged FreeGS flux budget')
    call require(trace_distances(2) <= minval(budgets(5:6)), &
                 'native traced LCFS at unchanged X-point position budget')
    print '(a,3es16.7)', 'native LCFS: distance (m), flux-distance, flux budget ', &
        trace_distances(2), trace_gaps(2), gap_budget
    print '(a)', 'candidate comparison; continuous bounds unavailable'
    if (failures /= 0) error stop 'FreeGS parity budgets failed'
contains
    subroutine run_case(resolution)
        integer, intent(in) :: resolution
        type(fortnum_status_t) :: root_status
        real(dp) :: amplitudes(2), constraint_residual(2)

        ! Artificial vacuum curve only: all four physical coils remain outside.
        ! No reference flux/LCFS value is used as a Dirichlet condition.
        n = 4*resolution
        allocate (contour(2, n))
        do i = 1, n
            call outer_curve(real(i - 1, dp)/n, contour(:, i))
        end do
        call initialize_gs(state, contour, 0.5_dp*(1.9_dp/resolution)**2, status, &
                           element_order=3, curved=.true., curve=outer_curve)
        call check(status, 'mesh')
        cfg%mode = GS_FB_XPOINT
        cfg%axis_sign = 1
        cfg%axis_guess = reference(10:11, 3)
        cfg%xpoint_guess = reference(12:13, 3)
        cfg%search_radius = 0.2_dp
        call gs_fb_setup(state, sys, law, cfg, coils, status)
        call check(status, 'FEM-BEM setup')
        do i = 1, size(state%psi)
            state%psi(i) = seed_value(state%points(:, i))
        end do
        call gs_fb_newton(state, sys, currents, status, history)
        call check(status, 'native free-boundary Newton')
        print *, 'mesh ', resolution, size(state%psi), ' Newton ', history%iterations
        print *, 'residuals ', history%residual
        ! ConstrainPaxisIp fixes two amplitudes, not just a frozen derivative
        ! table. Use the existing cubic law and native equilibrium/readback
        ! equations as a test consumer of FortNum's two-variable root solve.
        ! Finite differences here belong to that independent setup driver;
        ! production free-boundary derivatives remain equation-level.
        call multiroot_hybrids(profile_constraints, 2, [1.0_dp, 1.0_dp], &
                               amplitudes, root_status, ftol=1.0e-9_dp, max_iter=6)
        if (.not. status_ok(root_status)) then
            print *, root_status%msg
            error stop 'pressure/current normalization failed'
        end if
        call profile_constraints(amplitudes, constraint_residual)
        call require(maxval(abs(constraint_residual)) < 1.0e-9_dp, &
                     'same pressure/current constraints')
        law = sys%law
        print '(a,2es19.10)', 'pprime(0), FFprime(0) ', law%pp(1, 1), law%ff(1, 1)
        seed_error = maxval(abs(state%psi - [(seed_value(state%points(:, i)), &
                                              i=1, size(state%psi))]))
        print *, 'change from seed ', seed_error
        call gs_fb_evaluate(state, sys, state%psi, currents, ev, status, &
                            jacobian=.false.)
        call check(status, 'final residual/current')
        values(1:6) = [ev%psi_axis, ev%psi_bound, ev%axis, ev%xpoint]
        call gs_fb_trace_lcfs(state, state%psi, ev%axis, ev%psi_bound, &
                             1.0_dp, 1024, chain, status, xpoints=xpoints)
        call check(status, 'native separatrix trace')
        call require(size(xpoints, 2) == 1, 'one bounding lower X-point')
        if (size(xpoints, 2) == 1) then
            call require(norm2(xpoints(:, 1)-ev%xpoint) < 2e-14_dp, &
                         'trace returns the solved bounding X-point')
        end if
        call require(minval(norm2(chain-spread(ev%xpoint, 2, size(chain, 2)), &
                     dim=1)) < 2e-14_dp, 'chain retains the exact X-point')
        trace_gap = 0
        do i = 1, size(chain, 2)
            call sample_gs(state, chain(1, i), chain(2, i), 2.0_dp, &
                           psi_probe, br, bz, bt, status)
            call check(status, 'native LCFS point sample')
            trace_gap = max(trace_gap, abs(psi_probe-ev%psi_bound))
        end do
        call require(trace_gap < 2e-13_dp, 'native LCFS flux roots near roundoff')
        print '(a,es16.7)', 'native LCFS maximum root flux gap ', trace_gap
        cfg%xpoint_guess = reference(14:15, 3)
        call gs_fb_lcfs_levels(state, state%psi, cfg, ev%axis, cfg%xpoint_guess, &
                               pa, pb, axis, xp, cell, ga, gb, status)
        call check(status, 'upper X point')
        values(7:8) = xp
        law%axis = ev%psi_axis - ev%psi_bound
        call gs_law_primitives(law, law%axis, paxis, f, status)
        call check(status, 'axis pressure')
        call gs_equilibrium(state, eq, status, law%f_edge, 0.0_dp, law=law, &
                            psi_bound=ev%psi_bound, plasma_cell=ev%plasma_cell)
        call check(status, 'core evaluator')
        flux_levels(1) = ev%psi_axis
        flux_levels(2:) = ev%psi_axis + levels*(ev%psi_bound - ev%psi_axis)
        select type (field => eq%field)
        class is (labelled_field_t)
            call field%label_axis(axis, flux_levels(1), status)
            call check(status, 'readback axis')
            call field%label_coarea(flux_levels, 96, coarea, status)
        class default
            error stop 'missing flux label capability'
        end select
        call check(status, 'coarea readback')
        do i = 1, 4
            call gs_law_primitives(law, flux_levels(i + 1) - ev%psi_bound, &
                                   pressure, f, status)
            call check(status, 'F readback')
            ! FreeGS q is positive for this f>0, Ip>0 example. KIN6D's COCOS 3
            ! signed q has the opposite sign (sigma=sign(psi_b-psi_a)).
            values(8 + i) = sign(1.0_dp, ev%psi_bound - ev%psi_axis) &
                            *f*coarea(i + 1)/(2*pi)
        end do
        call integrals(volume, p_integral, b_integral)
        open (newunit=unit, file=trim(path)//'/lcfs.dat', status='old')
        max_gap = 0
        trace_distances(rung) = 0
        trace_gaps(rung) = 0
        r_geo = 0
        do j = 1, 64
            call data_line(unit, line)
            read (line, *) pt
            call sample_gs(state, pt(1), pt(2), 2.0_dp, psi_probe, br, bz, bt, status)
            call check(status, 'LCFS point sample')
            gap = abs(psi_probe - ev%psi_bound)
            max_gap = max(gap, max_gap)
            distance = polyline_distance(pt, chain)
            trace_distances(rung) = max(distance, trace_distances(rung))
            trace_gaps(rung) = max(trace_gaps(rung), distance*pt(1)*hypot(br, bz))
        end do
        close (unit)
        ! FreeGS internalInductance2: the normalization uses magnetic-axis R.
        r_geo = ev%axis(1)
        values(13) = 2*gs_mu0*p_integral/b_integral
        values(14) = 2*b_integral/(gs_mu0**2*r_geo*ev%plasma_current**2)
        values(15:16) = [ev%plasma_current, volume]
        print '(a,16es19.10)', 'observables ', values
        if (.not. all(ieee_is_finite(values))) error stop 'nonfinite readback'
        deallocate (contour)
    end subroutine

    function polyline_distance(point, chain) result(distance)
        ! Independent geometric comparison with the returned native contour.
        real(dp), intent(in) :: point(2), chain(:, :)
        real(dp) :: distance, edge(2), offset(2), t, length2
        integer :: i, next
        distance = huge(1.0_dp)
        do i = 1, size(chain, 2)
            next = mod(i, size(chain, 2))+1
            edge = chain(:, next)-chain(:, i)
            offset = point-chain(:, i)
            length2 = dot_product(edge, edge)
            t = 0
            if (length2 > 0) t = max(0.0_dp, min(1.0_dp, &
                dot_product(offset, edge)/length2))
            distance = min(distance, norm2(offset-t*edge))
        end do
    end function polyline_distance

    subroutine profile_constraints(amplitudes, residual, ctx)
        real(dp), intent(in) :: amplitudes(:)
        real(dp), intent(out) :: residual(:)
        class(*), intent(in), optional :: ctx
        type(gs_profile_law_t) :: executed
        real(dp) :: axis_pressure, axis_f
        integer :: stat
        sys%law%pp(:, 1) = reference(6, 3)*amplitudes(1) &
                           *[1.0_dp, -2.0_dp, 1.0_dp, 0.0_dp]
        sys%law%ff(:, 1) = reference(7, 3)*amplitudes(2) &
                           *[1.0_dp, -2.0_dp, 1.0_dp, 0.0_dp]
        call gs_fb_newton(state, sys, currents, stat)
        call check(stat, 'normalized profile forward solve')
        call gs_fb_evaluate(state, sys, state%psi, currents, ev, stat, &
                            jacobian=.false.)
        call check(stat, 'normalized profile current')
        executed = sys%law
        executed%axis = ev%psi_axis - ev%psi_bound
        call gs_law_primitives(executed, executed%axis, axis_pressure, axis_f, stat)
        call check(stat, 'normalized profile pressure')
        residual = [axis_pressure/1000.0_dp - 1, ev%plasma_current/2.0e5_dp - 1]
    end subroutine

    subroutine require(ok, what)
        logical, intent(in) :: ok
        character(*), intent(in) :: what
        if (ok) return
        failures = failures + 1
        print *, 'FAIL: ', what
    end subroutine

    subroutine outer_curve(t, point)
        real(dp), intent(in) :: t
        real(dp), intent(out) :: point(2)
        point = [1.0_dp + 0.8_dp*cos(2*pi*t), 0.95_dp*sin(2*pi*t)]
    end subroutine

    subroutine check(code, what)
        integer, intent(in) :: code
        character(*), intent(in) :: what
        if (code == 0) return
        print *, 'FAIL: ', what, code
        error stop 'FreeGS parity setup/solve/readback failed'
    end subroutine

    subroutine data_line(u, text)
        integer, intent(in) :: u
        character(*), intent(out) :: text
        do
            read (u, '(a)') text
            if (len_trim(text) == 0) cycle
            if (text(1:1) /= '#') return
        end do
    end subroutine

    real(dp) function seed_value(point) result(value)
        real(dp), intent(in) :: point(2)
        real(dp) :: x, y, wx(4), wy(4)
        integer :: ix, iy, k
        x = (point(1) - 0.2_dp)*24/1.45_dp
        y = (point(2) + 0.95_dp)*32/1.9_dp
        ix = max(1, min(22, floor(x)))
        iy = max(1, min(30, floor(y)))
        call lagrange_weights(4, x - (ix - 1), &
                              [0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp], wx)
        call lagrange_weights(4, y - (iy - 1), &
                              [0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp], wy)
        value = 0
        do k = 1, 4
            value = value + wx(k)*dot_product(wy, seed(iy:iy + 3, ix + k - 1))
        end do
    end function

    subroutine integrals(vol, pint, bint)
        real(dp), intent(out) :: vol, pint, bint
        real(dp), allocatable :: x(:), y(:), w(:)
        real(dp) :: shape(10), grad(2, 10), hess(2, 2, 10), point(2), det
        real(dp) :: ph, g(2), p, field, dv
        integer :: e, q, stat, bad
        call triangle_duffy_quadrature(16, x, y, w, stat)
        call check(stat, 'observable quadrature')
        vol = 0; pint = 0; bint = 0; bad = 0
        !$omp parallel do default(shared) private(e,q,shape,grad,hess,point,det, &
        !$omp ph,g,p,field,dv,stat) reduction(+:vol,pint,bint,bad) schedule(static)
        do e = 1, size(state%triangles, 2)
            if (.not. ev%plasma_cell(e)) cycle
            do q = 1, size(w)
                call gs_cell_data(state, e, [x(q), y(q)], point, det, shape, grad, &
                                  hess, stat)
                if (stat /= 0) then
                    bad = bad + 1
                    cycle
                end if
                ph = dot_product(shape, state%psi(state%triangles(:, e))) &
                     - ev%psi_bound
                if (ph <= 0) cycle
                call gs_law_primitives(law, ph, p, field, stat)
                if (stat /= 0) bad = bad + 1
                g = matmul(grad, state%psi(state%triangles(:, e)))
                dv = 2*pi*point(1)*w(q)*det
                vol = vol + dv
                pint = pint + dv*p
                bint = bint + dv*sum(g**2)/point(1)**2
            end do
        end do
        !$omp end parallel do
        call check(bad, 'observable integrands')
    end subroutine
end program
