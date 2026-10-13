module nn_island_case
    ! NN2 finite-pressure island: the NN1 helical-coil vacuum field on the
    ! annulus u in [0.3, 1.2] (docs/equilibrium-nonnested-nn2.md, gate 4)
    ! with a pressure source in the shell u in [0.75, 1.05] (u = distance
    ! from the coil axis R = 6, Z = 0) and p = 0 on both walls. The island
    ! (u ~ 0.41..0.57 at phi = 0, Z = 0) then lies on the inner pressure
    ! gradient, where the source vanishes. The inner wall is a floating
    ! conductor (zero net current into it, one potential).
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    real(dp), parameter :: s0 = 3.0e6_dp ! source peak per unit chi_perp
contains
    function source(x) result(s)
        real(dp), intent(in) :: x(2)
        real(dp) :: s, u
        u = norm2(x - [6.0_dp, 0.0_dp])
        s = s0*max(0.0_dp, 1 - ((u - 0.9_dp)/0.15_dp)**2)**2
    end function
end module

program test_nonnested_island
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use nn_island_case, only: source
    use nn_vacuum_fixture, only: pi, r0, ftor, field, active_mode, coil_mu0, &
        make_coil, mean_psi, gauge_potential, periodic_radius, &
        stability, island_edge, annulus_mesh
    use kin6d_mhd_mesh, only: flux_mesh_t
    use fortfem_meridian_feec_n, only: meridian_lagrange_nodes
    use kin6d_equilibrium_nonnested, only: nn_report_t, nn_build, nn_solve, &
        nn_vacuum_mode, nn_point
    use kin6d_equilibrium_state, only: EQ_VALUES, EQ_P
    use kin6d_nonnested_residual, only: nn_residual
    use kin6d_equilibrium_residual, only: equilibrium_residual_t
    use kin6d_nonnested_coupled, only: nc_problem_t, nc_layout_t, nc_report_t, &
        nc_layout, nc_solve
    implicit none
    character(16) :: mode, arg
    real(dp), allocatable :: xvac(:), x(:), ratios(:)
    real(dp) :: flat(8), oxdiff(8), d3(4), ratio, edges(2)
    integer :: k, modes, nr, i, nlev
    type(nc_layout_t) :: lay
    type(nc_problem_t) :: problem
    call get_command_argument(1, mode)
    call get_command_argument(2, arg)
    read (arg, *) k
    call get_command_argument(3, arg)
    read (arg, *) nr
    call get_command_argument(4, arg)
    read (arg, *) modes
    if (k < 2 .or. k > 4 .or. nr < 2 .or. modes < 1) error stop 'arguments'
    call make_coil()
    select case (trim(mode))
    case ('ratio')
        ! Flattening against chi_par/chi_perp at a fixed mesh, each ratio
        ! solved from the vacuum (zero-pressure) start.
        ratios = [1.0_dp, 1.0e2_dp, 1.0e4_dp, 1.0e6_dp]
        if (command_argument_count() >= 5) then
            call get_command_argument(5, arg)
            read (arg, *) ratio
            do while (ratios(size(ratios)) < ratio/2)
                ratios = [ratios, 100*ratios(size(ratios))]
            end do
        end if
        call setup(nr)
        do i = 1, size(ratios)
            call solve(ratios(i), flat(i), d3(1), i == size(ratios), oxdiff(i), &
                .true.)
        end do
        print '(a,8es11.3)', 'O-point gradient ratio G ', flat(:size(ratios))
        print '(a,8es11.3)', '|p_O - p_X|/p_max (signed crossing, no gate) ', &
            oxdiff(:size(ratios))
        ! G = 1 for a profile the island does not affect, 0 when flat at O.
        if (any(flat(2:size(ratios)) > flat(:size(ratios) - 1) + 0.01_dp)) &
            error stop 'flattening not monotone in the ratio'
        if (flat(size(ratios)) > 0.5_dp*flat(1)) error stop 'no flattening'
        print '(a)', 'PASS: NN2 island pressure flattening against the ratio'
    case ('ladder')
        ! D3 against resolution at a fixed ratio, nr, 2 nr, ...
        call get_command_argument(5, arg)
        read (arg, *) ratio
        nlev = 3
        if (command_argument_count() >= 6) then
            call get_command_argument(6, arg)
            read (arg, *) nlev
        end if
        do i = 1, nlev
            call setup(nr*2**(i - 1))
            ! The O/X bisection needs nr >= 8 (piecewise-smooth FEEC map).
            call solve(ratio, flat(i), d3(i), .false., oxdiff(i), &
                nr*2**(i - 1) >= 8)
        end do
        print '(a,4es11.3)', 'D3 ladder ', d3(:nlev)
        print '(a,4es11.3)', 'O-point gradient ratio G ladder ', flat(:nlev)
        print '(a,f8.3)', 'D3 order ', log(d3(1)/d3(nlev))/log(2.0_dp**(nlev - 1))
        if (any(d3(2:nlev) >= d3(:nlev - 1))) error stop 'D3 does not fall'
        print '(a)', 'PASS: NN2 island D3 falls with resolution'
    case default
        error stop 'mode'
    end select
contains
    subroutine setup(nrad)
        ! Vacuum state (NN1 solves) mapped to the coupled block layout:
        ! mode 0 is R A_phi, mode n the Nedelec A_pol with 2 A_n = c - i s.
        integer, intent(in) :: nrad
        type(flux_mesh_t) :: fm
        type(nn_report_t) :: report
        real(dp), allocatable :: nodes(:, :)
        complex(dp), allocatable :: bc(:), seed(:), a(:)
        integer :: st, n, j
        if (allocated(xvac)) deallocate (xvac)
        call annulus_mesh(nrad, k, fm)
        call nn_build(fm, k, modes, [r0, 0.0_dp], 1.0_dp, ftor, coil_mu0, field, st)
        if (st /= 0) error stop 'build'
        call meridian_lagrange_nodes(field%fs, nodes)
        lay = nc_layout(field)
        allocate (bc(lay%n0), seed(lay%n0), xvac(lay%n))
        bc = 0; seed = 0; xvac = 0
        !$omp parallel do default(shared) private(j) schedule(static)
        do j = 1, lay%n0
            if (field%fs%sp(0)%bnd(j)) bc(j) = mean_psi(nodes(:, j))
        end do
        !$omp end parallel do
        call nn_solve(field, 0, bc, seed, report, st, solution=a)
        if (st /= 0) error stop 'vacuum n0'
        xvac(lay%ia + field%fs%nn + 1:lay%ia + lay%na1) = real(a, dp)
        do n = 1, modes
            active_mode = n
            call nn_vacuum_mode(field, n, gauge_potential, report, st, solution=a)
            if (st /= 0) error stop 'vacuum mode'
            xvac(lay%ia + (2*n - 1)*lay%na1 + 1:lay%ia + (2*n - 1)*lay%na1 &
                + field%fs%nn) = 2*real(a, dp)
            xvac(lay%ia + 2*n*lay%na1 + 1:lay%ia + 2*n*lay%na1 + field%fs%nn) = &
                -2*aimag(a)
        end do
        ! Inner wall nodes form the floating conductor.
        problem%floating = field%fs%sp(0)%bnd .and. &
            norm2(nodes - spread([r0, 0.0_dp], 2, lay%n0), 1) < 0.75_dp
        problem%source => source
        print '(a,i0,a,i0,a,i0,a,i0)', 'P', k, ' nr=', nrad, ' N=', modes, &
            ' unknowns ', lay%n
    end subroutine

    subroutine solve(ratio_i, flattening, dual, edges_wanted, ox, measure)
        ! flattening: G = |dp/dR| at the O point (central difference over
        ! +-0.03 m on Z = 0, phi = 0) over the mean slope of the chord
        ! R_O +- 0.12 m, which spans the island (width ~ 0.18 m).
        real(dp), intent(in) :: ratio_i
        real(dp), intent(out) :: flattening, dual, ox
        logical, intent(in) :: edges_wanted, measure
        type(nc_report_t) :: report
        type(equilibrium_residual_t) :: stress
        real(dp) :: oh, xh, opoint(2), xpoint(2), trace_o, trace_x, defect, closure
        real(dp) :: val(EQ_VALUES), p_o, p_x, p_max, chord(9), rr
        integer :: st, j
        problem%chi_perp = 1; problem%chi_par = ratio_i
        problem%g_par = 1; problem%g_perp = 1/ratio_i
        problem%current = 0
        x = xvac
        call nc_solve(field, problem, x, report, st)
        print '(a,es9.1,a,15es9.1)', 'ratio ', ratio_i, ' Newton ', &
            report%history(:report%newton_steps)
        if (st /= 0) error stop 'coupled Newton'
        ! div curl A = 0 is algebraic; its rounding grows with |A| ~ 20 Wb
        ! and the P_k derivative scale, unlike the unit-scale Solov'ev case.
        print '(a,es9.1,a,es10.2)', 'ratio ', ratio_i, ' div B rounding ', &
            report%divergence
        if (report%divergence > 1.0e-10_dp) error stop 'closed 2-form'
        ! lambda balances <dA, J> ~ |A|/mu0 (Ampere), so it is compared there.
        if (report%multiplier > 1.0e-10_dp*maxval(abs(x(:lay%na)))/coil_mu0) &
            error stop 'gauge multiplier'
        call nn_residual(field, 1, stress, st)
        if (st /= 0) error stop 'D3'
        dual = stress%normalized
        flattening = -1; ox = 0
        print '(a,es9.1,a,es11.3)', 'ratio ', ratio_i, ' D3 ', dual
        if (.not. measure) return
        oh = periodic_radius(5.49_dp, 5.54_dp, .false., 512)
        xh = periodic_radius(6.69_dp, 6.72_dp, .false., 512)
        ! The coil is stellarator symmetric (phi, Z) -> (-phi, -Z), so the
        ! period-1 points at phi = 0 lie on Z = 0 and the half-turn bisection
        ! finds them; the FEEC map is only piecewise smooth, so no 2D Newton
        ! refinement. The one-turn closure is reported.
        opoint = [oh, 0.0_dp]; xpoint = [xh, 0.0_dp]
        call stability(opoint, .false., trace_o, defect)
        if (abs(trace_o) >= 2) error stop 'O point'
        call stability(xpoint, .false., trace_x, closure)
        if (trace_x <= 2) error stop 'X point'
        print '(a,2es10.2,a,2f8.3)', 'closure O, X ', defect, closure, &
            ' traces ', trace_o, trace_x
        call nn_point(field, [opoint, 0.0_dp], val, st)
        if (st /= 0) error stop 'p at O'
        p_o = val(EQ_P)
        call nn_point(field, [xpoint, 0.0_dp], val, st)
        if (st /= 0) error stop 'p at X'
        p_x = val(EQ_P)
        p_max = 0
        do j = 1, 181
            rr = r0 - 1.2_dp + 2.4_dp*(j - 1)/180
            if (abs(rr - r0) < 0.3_dp .or. abs(rr - r0) > 1.2_dp) cycle
            call nn_point(field, [rr, 0.0_dp, 0.0_dp], val, st)
            if (st /= 0) error stop 'p profile'
            p_max = max(p_max, val(EQ_P))
        end do
        ox = (p_o - p_x)/p_max
        ! p/p_max on the chord Z = 0, phi = 0 through the O point.
        do j = 1, 9
            call nn_point(field, [opoint(1) + 0.03_dp*(j - 5), 0.0_dp, 0.0_dp], &
                val, st)
            if (st /= 0) error stop 'chord'
            chord(j) = val(EQ_P)/p_max
        end do
        flattening = abs(chord(6) - chord(4))/0.06_dp &
            /(abs(chord(9) - chord(1))/0.24_dp)
        print '(a,es9.1,a,2f9.5,a,2f9.5,a,es11.3)', 'ratio ', ratio_i, ' O ', &
            opoint, ' X ', xpoint, ' p_max ', p_max
        print '(a,es9.1,a,3es11.3,a,es10.2,a,es10.2)', 'ratio ', ratio_i, &
            ' p_O p_X G ', p_o/p_max, p_x/p_max, flattening, &
            ' D3 ', dual, ' V ', report%voltage
        print '(a,es9.1,a,es10.2)', 'ratio ', ratio_i, ' lambda mu0/|A| ', &
            report%multiplier*coil_mu0/maxval(abs(x(:lay%na)))
        print '(a,9f8.4)', 'chord p/p_max (dR = 0.03) ', chord
        if (edges_wanted) then
            edges(1) = island_edge(0.36_dp, 0.48_dp, .false., 32.0_dp, 256)
            edges(2) = island_edge(0.63_dp, 0.48_dp, .false., 32.0_dp, 256)
            print '(a,3es12.4)', 'island edges (u at Z=0, inner side) and width ', &
                edges, edges(2) - edges(1)
        end if
    end subroutine
end program
