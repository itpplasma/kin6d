module island_a_study_trace
    ! Field-line tracing on the converged provider field that marks lines
    ! leaving the mesh instead of stopping (open lines reach the walls).
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_equilibrium_nonnested, only: free_b_field_t, nn_velocity
    implicit none
    type(free_b_field_t), pointer :: traced => null()
    logical :: left_mesh = .false.
    !$omp threadprivate(left_mesh)
contains
    subroutine rhs(phi, y, dy, ctx)
        real(dp), intent(in) :: phi, y(:)
        real(dp), intent(out) :: dy(:)
        class(*), intent(in), optional :: ctx
        integer :: st
        call nn_velocity(traced, [y, phi], dy, st)
        if (st /= 0) then
            left_mesh = .true.; dy = 0
        end if
    end subroutine
end module

program study_island_closure_a
    ! Closure A+D island study on the helical-coil annulus (not a CTest):
    !   scan k nr outdir eps...    flattening G, island width and residuals
    !                              against chi_par/chi_perp (ISLAND_SCAN_LOG10)
    !   poincare k nr outdir eps ratio   Poincare section and p(R, Z) at phi = 0
    !   connect k nr outdir eps ratio    wall connection (turns, both directions)
    !                              of lines seeded on the Z = 0 diameter at
    !                              phi = 0, with p/p_max at each seed
    ! ISLAND_MODES sets the toroidal truncation n = 0..N (default 1).
    ! scan also writes plane.dat (min p/p_max and the area fractions with
    ! p < 0 and p < 0.01 p_max on the phi = 0 plane, solve wall time) and
    ! diameter.dat (p/p_max on the whole Z = 0 diameter); eps = 0 skips the
    ! island tracing (G, W, u_O are NaN).
    ! G is the O-point slope over +-W/6 divided by the chord slope over
    ! +-2W/3 on Z = 0, phi = 0 (fractions of the measured width W, frozen
    ! across the comparison; docs/island-equilibrium.md).
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use island_a_fixture, only: island_a_setup, shell_source
    use nn_vacuum_fixture, only: field, r0, pi
    use island_a_study_trace, only: traced, left_mesh, rhs
    use fortnum_ode_cash_karp, only: cash_karp_step
    use kin6d_nonnested_coupled, only: nc_problem_t
    use kin6d_island_provider, only: island_provider_t, island_provider_init, &
        island_provider_solve, island_provider_sample, IP_OK, IP_CHI_PAR, &
        island_provider_provenance
    implicit none
    type(island_provider_t), target :: prov
    character(256) :: mode, arg, outdir
    integer :: k, nr, i
    real(dp), allocatable :: epss(:)
    real(dp) :: ratio
    call get_command_argument(1, mode)
    call get_command_argument(2, arg); read (arg, *) k
    call get_command_argument(3, arg); read (arg, *) nr
    call get_command_argument(4, outdir)
    select case (trim(mode))
    case ('scan')
        allocate (epss(command_argument_count() - 4))
        do i = 1, size(epss)
            call get_command_argument(4 + i, arg); read (arg, *) epss(i)
        end do
        call scan()
    case ('connect')
        call get_command_argument(5, arg); allocate (epss(1)); read (arg, *) epss(1)
        call get_command_argument(6, arg); read (arg, *) ratio
        call connect(epss(1), ratio)
    case ('poincare')
        call get_command_argument(5, arg); allocate (epss(1)); read (arg, *) epss(1)
        call get_command_argument(6, arg); read (arg, *) ratio
        call poincare(epss(1), ratio)
    case default
        error stop 'mode'
    end select
contains
    subroutine start(eps, ratio0)
        real(dp), intent(in) :: eps, ratio0
        type(nc_problem_t) :: problem
        real(dp), allocatable :: base(:), pert(:)
        integer :: st, modes, ln
        character(16) :: env
        ! Toroidal modes n = 0..ISLAND_MODES (default 1) of A, p and Phi.
        modes = 1
        call get_environment_variable('ISLAND_MODES', env, ln)
        if (ln > 0) read (env, *) modes
        call island_a_setup(nr, k, modes, base, pert, problem)
        call island_provider_init(prov, field, problem, base, pert, shell_source, &
            [1.0_dp, eps, 0.0_dp, ratio0], 0.9_dp, 1.0e5_dp, &
            'helical-coil-annulus', st)
        if (st /= IP_OK) error stop 'init'
    end subroutine

    real(dp) function p_at(rr) result(p)
        real(dp), intent(in) :: rr
        real(dp) :: v(4, 1)
        integer :: st
        call island_provider_sample(prov, reshape([rr, 0.0_dp, 0.0_dp], [3, 1]), &
            v, st)
        if (st /= IP_OK) error stop 'sample'
        p = v(4, 1)
    end function

    subroutine scan()
!$      use omp_lib, only: omp_get_wtime
        real(dp) :: eps, ratio, ro, rx, uo, w, w1, e1, e2, g, a, b, pmax, rr
        real(dp) :: t0, wall, pmin, fneg, fzero, nan
        integer :: ie, ir, st, us, uc, up, ud, j, nratio, ln
        real(dp) :: grid(3)
        character(64) :: env
        ! log10 ratio grid: ISLAND_SCAN_LOG10="first last step" (default 0 8 0.5);
        ! the provider warm-starts each ratio from the previous one.
        grid = [0.0_dp, 8.0_dp, 0.5_dp]
        call get_environment_variable('ISLAND_SCAN_LOG10', env, ln)
        if (ln > 0) read (env, *) grid
        nratio = nint((grid(2) - grid(1))/grid(3)) + 1
        open (newunit=us, file=trim(outdir)//'/scan.dat', status='replace')
        open (newunit=uc, file=trim(outdir)//'/chords.dat', status='replace')
        write (us, '(a)') '# eps ratio W(traced at first/last ratio) u_O G p_O/p_max p_X/p_max p_max eps_parallel '// &
            'pressure_rows D3 newton stages'
        write (uc, '(a)') '# eps ratio R p/p_max (Z = 0, phi = 0, inner side)'
        open (newunit=up, file=trim(outdir)//'/plane.dat', status='replace')
        write (up, '(a)') '# eps ratio min(p)/p_max frac(p<0) frac(p<0.01 p_max) '// &
            'solve_wall_s (phi = 0 plane, 241 x 241 grid in the annulus)'
        open (newunit=ud, file=trim(outdir)//'/diameter.dat', status='replace')
        write (ud, '(a)') '# eps ratio R p/p_max (Z = 0, phi = 0, whole diameter)'
        nan = ieee_value(0.0_dp, ieee_quiet_nan)
        do ie = 1, size(epss)
            eps = epss(ie)
            call start(eps, 1.0_dp)
            do ir = 1, nratio
                ratio = 10.0_dp**(grid(1) + grid(3)*(ir - 1))
                prov%theta(IP_CHI_PAR) = ratio
                t0 = 0
!$              t0 = omp_get_wtime()
                call island_provider_solve(prov, prov%theta, st)
                wall = 0
!$              wall = omp_get_wtime() - t0
                if (st /= IP_OK) then
                    print '(a,es9.2,a,es9.2,a)', 'eps ', eps, ' ratio ', ratio, &
                        ' solve FAILED'
                    exit
                end if
                field = prov%field
                traced => prov%field
                call plane_metrics(pmin, fneg, fzero)
                write (up, '(2es14.6,3es14.6,f12.2)') eps, ratio, pmin, fneg, fzero, wall
                flush (up)
                if (eps > 0) then
                    ro = half_turn_root(5.47_dp, 5.56_dp)
                    rx = half_turn_root(6.66_dp, 6.75_dp)
                    uo = r0 - ro
                    ! W is traced at the first and last ratio; G uses the
                    ! first-ratio W for the whole scan (frozen offsets).
                    if (ir == 1 .or. ir == nratio) then
                        e1 = edge(max(0.31_dp, uo - 0.2_dp), uo)
                        e2 = edge(uo + 0.2_dp, uo)
                        w = e2 - e1
                        if (ir == 1) w1 = w
                    end if
                    a = w1/6; b = 2*w1/3
                    g = abs(p_at(ro + a) - p_at(ro - a))/(2*a) &
                        /(abs(p_at(ro + b) - p_at(ro - b))/(2*b))
                else
                    ro = r0 - 0.48_dp; rx = r0 + 0.70_dp
                    uo = nan; w = nan; g = nan
                end if
                pmax = 0
                do j = 0, 180
                    rr = r0 - 1.18_dp + 0.86_dp*j/180
                    pmax = max(pmax, p_at(rr))
                end do
                do j = 0, 180
                    rr = r0 - 1.18_dp + 0.86_dp*j/180
                    write (uc, '(4es16.8)') eps, ratio, rr, p_at(rr)/pmax
                end do
                write (uc, *)
                do j = 0, 480
                    rr = r0 - 1.199_dp + 2.398_dp*j/480
                    if (abs(rr - r0) < 0.301_dp) cycle
                    write (ud, '(4es16.8)') eps, ratio, rr, p_at(rr)/pmax
                end do
                write (ud, *)
                flush (ud)
                write (us, '(4es14.6,4es14.6,3es12.4,2i5)') eps, ratio, w, uo, g, &
                    p_at(ro)/pmax, p_at(rx)/pmax, pmax, prov%report%epsilon_parallel, &
                    prov%report%pressure, prov%report%force%normalized, &
                    prov%solve_report%newton_steps, prov%continued_stages
                flush (us); flush (uc)
                print '(a,es9.2,a,es9.2,a,f8.4,a,f8.4,a,f8.4,a,es9.2,a,es9.2,a,i0)', &
                    'eps ', eps, ' ratio ', ratio, ' W ', w, ' u_O ', uo, ' G ', g, &
                    ' eps_par ', prov%report%epsilon_parallel, ' D3 ', &
                    prov%report%force%normalized, ' stages ', prov%continued_stages
            end do
            print '(a)', island_provider_provenance(prov)
        end do
        close (us); close (uc); close (up); close (ud)
    end subroutine

    subroutine plane_metrics(pmin, fneg, fzero)
        ! p on the phi = 0 plane (the poincare mode's 241 x 241 grid restricted
        ! to the annulus): min p/p_max and area fractions of p < 0 and
        ! p < 0.01 p_max (equal-area cells).
        real(dp), intent(out) :: pmin, fneg, fzero
        integer, parameter :: ng = 241
        real(dp), allocatable :: xyz(:, :), v(:, :)
        real(dp) :: x(2)
        integer :: i, j, n, st
        allocate (xyz(3, ng*ng))
        n = 0
        do i = 1, ng
            do j = 1, ng
                x = [r0 - 1.25_dp + 2.5_dp*(j - 1)/(ng - 1), &
                    -1.25_dp + 2.5_dp*(i - 1)/(ng - 1)]
                if (norm2(x - [r0, 0.0_dp]) <= 0.301_dp) cycle
                if (norm2(x - [r0, 0.0_dp]) >= 1.199_dp) cycle
                n = n + 1
                xyz(:, n) = [x(1), 0.0_dp, x(2)]
            end do
        end do
        allocate (v(4, n))
        call island_provider_sample(prov, xyz(:, :n), v, st)
        if (st /= IP_OK) error stop 'plane metrics'
        pmin = minval(v(4, :))/maxval(v(4, :))
        fneg = count(v(4, :) < 0)/real(n, dp)
        fzero = count(v(4, :) < 0.01_dp*maxval(v(4, :)))/real(n, dp)
    end subroutine

    subroutine connect(eps, ratio)
        ! Wall connection of field lines seeded on the Z = 0 diameter at
        ! phi = 0: toroidal turns until the line leaves the annulus, forward
        ! and backward (nturn when it stays inside). A line that reaches a
        ! wall carries the Dirichlet value p = 0 in the chi_par -> infinity
        ! limit; a line that stays inside does not.
        real(dp), intent(in) :: eps, ratio
        integer, parameter :: nside = 240, nturn = 100, steps = 128
        real(dp) :: rs(2*nside), pseed(2*nside), pmax, u
        integer :: turns(2, 2*nside), s, uo, st
        real(dp), allocatable :: xyz(:, :), v(:, :)
        call start(eps, 1.0_dp)
        prov%theta(IP_CHI_PAR) = ratio
        call island_provider_solve(prov, prov%theta, st)
        if (st /= IP_OK) error stop 'connect solve'
        print '(a)', island_provider_provenance(prov)
        traced => prov%field
        do s = 1, nside
            u = 0.302_dp + 0.896_dp*(s - 1)/(nside - 1)
            rs(s) = r0 - u; rs(nside + s) = r0 + u
        end do
        allocate (xyz(3, 2*nside), v(4, 2*nside))
        xyz(1, :) = rs; xyz(2:3, :) = 0
        call island_provider_sample(prov, xyz, v, st)
        if (st /= IP_OK) error stop 'connect sample'
        pmax = maxval(v(4, :)); pseed = v(4, :)/pmax
        !$omp parallel do default(shared) private(s) schedule(dynamic)
        do s = 1, 2*nside
            turns(1, s) = turns_to_wall(rs(s), 1.0_dp, nturn, steps)
            turns(2, s) = turns_to_wall(rs(s), -1.0_dp, nturn, steps)
        end do
        !$omp end parallel do
        open (newunit=uo, file=trim(outdir)//'/connect.dat', status='replace')
        write (uo, '(a,i0,a)') '# R p/p_max turns_forward turns_backward (Z = 0, '// &
            'phi = 0; ', nturn, ' = no wall contact)'
        do s = 1, 2*nside
            write (uo, '(2es16.8,2i6)') rs(s), pseed(s), turns(:, s)
        end do
        close (uo)
    end subroutine

    integer function turns_to_wall(rr, direction, nturn, steps) result(t)
        ! Toroidal turns until the line from (rr, 0, phi = 0) leaves the
        ! annulus; nturn when it stays inside.
        real(dp), intent(in) :: rr, direction
        integer, intent(in) :: nturn, steps
        real(dp) :: y(2), k1(2), k2(2), k3(2), k4(2), k5(2), k6(2), tmp(2)
        real(dp) :: nxt(2), err(2), h
        integer :: i, j, nfev
        y = [rr, 0.0_dp]; left_mesh = .false.; nfev = 0
        h = direction*2*pi/steps
        do t = 0, nturn - 1
            do j = 1, steps
                i = t*steps + j
                call cash_karp_step(rhs, h*(i - 1), y, h, .false., k1, k2, k3, &
                    k4, k5, k6, tmp, nxt, err, nfev)
                y = nxt
                if (left_mesh) return
                if (norm2(y - [r0, 0.0_dp]) <= 0.3_dp .or. &
                    norm2(y - [r0, 0.0_dp]) >= 1.2_dp) return
            end do
        end do
        t = nturn
    end function

    real(dp) function half_turn_root(left, right) result(root)
        ! Period-1 point on Z = 0 at phi = 0 (the coil is stellarator
        ! symmetric): sign change of Z after half a toroidal turn, bracketed
        ! by three passes of nprobe parallel traces (nn_vacuum_fixture's
        ! periodic_radius predicate, bisection replaced by probes).
        real(dp), intent(in) :: left, right
        integer, parameter :: nprobe = 24
        real(dp) :: a, b, r(0:nprobe), z(0:nprobe)
        integer :: pass, j
        a = left; b = right
        do pass = 1, 3
            r = [(a + (b - a)*j/real(nprobe, dp), j=0, nprobe)]
            !$omp parallel do default(shared) private(j) schedule(dynamic)
            do j = 0, nprobe
                z(j) = half_turn_z(r(j))
            end do
            !$omp end parallel do
            do j = 1, nprobe
                if (z(j - 1)*z(j) <= 0) exit
            end do
            if (j > nprobe) error stop 'periodic orbit not bracketed'
            a = r(j - 1); b = r(j)
        end do
        root = (a + b)/2
    end function

    real(dp) function half_turn_z(rr) result(z)
        real(dp), intent(in) :: rr
        integer, parameter :: steps = 128
        real(dp) :: y(2), k1(2), k2(2), k3(2), k4(2), k5(2), k6(2), tmp(2), nxt(2)
        real(dp) :: err(2), h
        integer :: i, nfev
        y = [rr, 0.0_dp]; left_mesh = .false.; nfev = 0
        h = pi/steps
        do i = 1, steps
            call cash_karp_step(rhs, h*(i - 1), y, h, .false., k1, k2, k3, k4, &
                k5, k6, tmp, nxt, err, nfev)
            y = nxt
        end do
        if (left_mesh) error stop 'periodic trace left the mesh'
        z = y(2)
    end function

    real(dp) function edge(outside, inside) result(e)
        ! Island separatrix crossing on the inner chord (R = r0 - u, Z = 0)
        ! between a circulating start (outside) and the O point (inside):
        ! two passes of nprobe parallel 16-turn traces, each narrowing the
        ! bracket to the probe spacing. A line circulates when its poloidal
        ! angle about the coil, relative to the toroidal angle, spans more
        ! than 2 pi (nn_vacuum_fixture's orbit predicate).
        real(dp), intent(in) :: outside, inside
        integer, parameter :: nprobe = 24
        real(dp) :: a, b, u(0:nprobe)
        logical :: circ(0:nprobe)
        integer :: pass, j
        a = outside; b = inside
        do pass = 1, 2
            u = [(a + (b - a)*j/real(nprobe, dp), j=0, nprobe)]
            !$omp parallel do default(shared) private(j) schedule(dynamic)
            do j = 0, nprobe
                circ(j) = circulates(u(j))
            end do
            !$omp end parallel do
            if (.not. circ(0) .or. circ(nprobe)) error stop 'edge bracket'
            j = findloc(circ, .false., 1) - 1
            a = u(j - 1); b = u(j)
        end do
        e = (a + b)/2
    end function

    logical function circulates(u) result(c)
        real(dp), intent(in) :: u
        integer, parameter :: turns = 16, steps = 128
        real(dp) :: y(2), k1(2), k2(2), k3(2), k4(2), k5(2), k6(2), tmp(2), nxt(2)
        real(dp) :: err(2), h, theta, oldtheta, angle, lo, hi
        integer :: i, nfev
        y = [r0 - u, 0.0_dp]; left_mesh = .false.; nfev = 0
        oldtheta = atan2(y(2), y(1) - r0); angle = oldtheta; lo = angle; hi = angle
        h = 2*pi/steps
        c = .false.
        do i = 1, turns*steps
            call cash_karp_step(rhs, h*(i - 1), y, h, .false., k1, k2, k3, k4, &
                k5, k6, tmp, nxt, err, nfev)
            y = nxt
            if (left_mesh) error stop 'width trace left the mesh'
            theta = atan2(y(2), y(1) - r0)
            angle = angle + atan2(sin(theta - oldtheta), cos(theta - oldtheta))
            oldtheta = theta
            lo = min(lo, angle + i*h); hi = max(hi, angle + i*h)
            if (hi - lo > 2*pi) then
                c = .true.
                return
            end if
        end do
    end function

    subroutine poincare(eps, ratio)
        real(dp), intent(in) :: eps, ratio
        integer, parameter :: nseed = 72, nturn = 400, steps = 256
        integer, parameter :: ng = 241
        real(dp), allocatable :: pts(:, :, :), xyz(:, :), v(:, :)
        integer, allocatable :: reached(:), idx(:)
        logical, allocatable :: inside(:)
        real(dp) :: y(2), k1(2), k2(2), k3(2), k4(2), k5(2), k6(2), tmp(2), nxt(2)
        real(dp) :: err(2), h, u, nan
        integer :: s, t, j, nfev, st, uo, i
        real(dp), allocatable :: vi(:, :)
        call start(eps, 1.0_dp)
        prov%theta(IP_CHI_PAR) = ratio
        call island_provider_solve(prov, prov%theta, st)
        if (st /= IP_OK) error stop 'poincare solve'
        print '(a)', island_provider_provenance(prov)
        traced => prov%field
        allocate (pts(2, nturn, nseed), reached(nseed))
        h = 2*pi/steps
        !$omp parallel do default(shared) private(s, t, j, y, u, k1, k2, k3, k4, &
        !$omp& k5, k6, tmp, nxt, err, nfev) schedule(dynamic)
        do s = 1, nseed
            ! 48 seeds alternating between both sides of Z = 0, then 24 on
            ! the inner chord across the island (u in [0.30, 0.70]).
            if (s <= 48) then
                u = 0.32_dp + 0.86_dp*(s - 1)/47
                y = [r0 - u, 0.0_dp]
                if (mod(s, 2) == 0) y = [r0 + u, 0.0_dp]
            else
                u = 0.30_dp + 0.40_dp*(s - 49)/23
                y = [r0 - u, 0.0_dp]
            end if
            left_mesh = .false.; reached(s) = 0; nfev = 0
            do t = 1, nturn
                do j = 1, steps
                    call cash_karp_step(rhs, h*(j - 1), y, h, .false., k1, k2, k3, &
                        k4, k5, k6, tmp, nxt, err, nfev)
                    y = nxt
                    if (left_mesh) exit
                end do
                if (left_mesh) exit
                pts(:, t, s) = y; reached(s) = t
            end do
        end do
        !$omp end parallel do
        open (newunit=uo, file=trim(outdir)//'/poincare.dat', status='replace')
        write (uo, '(a)') '# seed R Z (phi = 0, one point per toroidal turn)'
        do s = 1, nseed
            do t = 1, reached(s)
                write (uo, '(i4,2es16.8)') s, pts(:, t, s)
            end do
        end do
        close (uo)
        allocate (xyz(3, ng*ng), v(4, ng*ng))
        do i = 1, ng
            do j = 1, ng
                xyz(:, (i - 1)*ng + j) = [r0 - 1.25_dp + 2.5_dp*(j - 1)/(ng - 1), &
                    0.0_dp, -1.25_dp + 2.5_dp*(i - 1)/(ng - 1)]
            end do
        end do
        nan = ieee_value(0.0_dp, ieee_quiet_nan)
        inside = [(norm2(xyz([1, 3], j) - [r0, 0.0_dp]) > 0.301_dp .and. &
            norm2(xyz([1, 3], j) - [r0, 0.0_dp]) < 1.199_dp, j=1, ng*ng)]
        idx = pack([(j, j=1, ng*ng)], inside)
        v = nan
        allocate (vi(4, size(idx)))
        call island_provider_sample(prov, xyz(:, idx), vi, st)
        if (st /= IP_OK) error stop 'pressure grid'
        v(:, idx) = vi
        open (newunit=uo, file=trim(outdir)//'/pressure.dat', status='replace')
        write (uo, '(a,i0)') '# R Z p (phi = 0; NaN outside the annulus) n=', ng
        do j = 1, ng*ng
            write (uo, '(3es16.8)') xyz(1, j), xyz(3, j), v(4, j)
        end do
        close (uo)
    end subroutine
end program
