program test_ideal_map_axisymmetric
    ! N1 axisymmetric reduction of the 3D frozen-in energy map (DESIGN 7.1/7.5).
    ! Exact Solovev states: A3 (circular-type, FF'=0) and a shaped one
    ! (kappa 1.61, delta 0.18, FF' /= 0). The frozen labels are the exact state
    ! pulled back through a smooth interior shift Y (Y = id on the boundary), so
    ! the seed X0 = id is not an equilibrium; the exact minimiser is X = Y up to
    ! relabelling. Oracles: (i) the exact flux, through the gauge-invariant
    ! surface defect psi_exact(X(y)) - chi(y); (ii) KIN6D's 2D GS solution on the
    ! same mesh, psi_2D(X(y)) - chi(y). Modes:
    !   derivatives      FD gradient/Hessian, rigid-rotation null, N=1 path
    !   leakage          N=1, nfp=3, n=3 interior seed: n/=0 content decays
    !   solovev <k>      A3 ladder, degree k, observed orders
    !   shaped <k>       shaped ladder, degree k, observed orders
    !   seed <k>, seedshaped <k>  frozen data from the nested-seed builder
    !                    (boundary + exact q(s), p(s), Phi_edge) with the seed
    !                    axis at the boundary-box centre (VMEC-like), offset from
    !                    the magnetic axis; the solver must return it
    !   gauge <k>        W_red curvature of field-line slides at the exact map
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64, output_unit
    use kin6d_grad_shafranov, only: gs_state_t, initialize_gs, solve_gs, &
        gs_element_jet
    use kin6d_equilibrium_chart, only: cylindrical_chart_t
    use kin6d_ideal_map, only: ideal_map_t, ideal_map_init, &
        ideal_map_freeze, ideal_map_energy, ideal_map_hessian, ideal_map_solve, &
        ideal_map_point, ideal_map_rdof, ideal_map_solve_report_t, &
        ideal_map_reduce, ideal_map_expand, ideal_map_prolong, ideal_map_residual, &
        ideal_map_eliminate_omega, IDEAL_MAP_BAD_INPUT, ideal_map_gauge_project
    use fortfem_element_to_csr, only: csr_matrix_t
    use kin6d_triangle_basis, only: triangle_basis_jet
    use fortfem_gauss_quadrature_2d, only: gauss_quadrature_triangle_t
    use kin6d_ideal_map, only: ideal_frozen_t
    use kin6d_ideal_seed, only: nested_seed_t, ideal_seed_init
    use ideal_map_cases, only: setup_case, psi_exact, boundary_point, label, &
        phi_edge, pullback_field_t, solovev_profile_t, pi, a, raxis, psi_axis, r0, &
        f_edge, seed_map_point, &
        pprime, ffprime, shaped
    use kin6d_equilibrium_harness, only: equilibrium_case_t, &
        equilibrium_resolution_t, equilibrium_start_t, equilibrium_run_report_t, &
        run_equilibrium_case, write_equilibrium_report
    implicit none
    type :: map_case_state_t
        type(gs_state_t) :: gs
        type(ideal_map_t) :: map
        type(ideal_map_t), allocatable :: previous
        class(ideal_frozen_t), allocatable :: field
        type(ideal_map_solve_report_t) :: solve
        real(dp) :: errors(4), elapsed
        logical :: sequenced = .false.
    end type
    type(pullback_field_t) :: pullback
    character(16) :: mode, argument
    integer :: degree
    logical :: verbose_solve = .false.
    integer :: quadrature_order = 7
    integer :: axis_degree = -1, axis_power = 1
    real(dp) :: initial_damping = 1.0e-6_dp, gauge_weight = -1
    real(dp) :: seed_tolerance = 1.0e-11_dp
    logical :: gauge_hessian_only = .false.
    ! Nested solve: omega eliminated exactly, Newton on W_red(nu) (no gauge form).
    logical :: omega_elimination = .false.
    integer :: max_newton = 600
    logical :: position_gauge = .false.
    real(dp) :: solve_tolerance = 1.0e-6_dp
    mode = 'solovev'; degree = 2
    if (command_argument_count() > 0) call get_command_argument(1, mode)
    if (command_argument_count() > 1) then
        call get_command_argument(2, argument)
        read (argument, *) degree
    end if
    if (command_argument_count() > 2) then
        call get_command_argument(3, argument)
        verbose_solve = trim(argument) /= 'p'
        position_gauge = trim(argument) == 'p' .or. trim(argument) == 'pv'
    end if
    if (command_argument_count() > 3) then
        call get_command_argument(4, argument)
        read (argument, *) quadrature_order
    end if
    if (command_argument_count() > 5) then
        call get_command_argument(5, argument)
        read (argument, *) axis_degree
        call get_command_argument(6, argument)
        read (argument, *) axis_power
    end if
    if (command_argument_count() > 6) then
        call get_command_argument(7, argument)
        read (argument, *) initial_damping
    end if
    if (command_argument_count() > 7) then
        call get_command_argument(8, argument)
        read (argument, *) gauge_weight
    end if
    if (command_argument_count() > 8) then
        call get_command_argument(9, argument)
        read (argument, *) seed_tolerance
    end if
    if (command_argument_count() > 9) then
        call get_command_argument(10, argument)
        gauge_hessian_only = trim(argument) == 'h'
        omega_elimination = trim(argument) == 'n'
    end if
    call setup_case(trim(mode) == 'shaped' .or. trim(mode) == 'seedshaped')
    select case (trim(mode))
    case ('derivatives')
        call derivative_checks()
    case ('reduced')
        call reduced_derivative_checks()
    case ('leakage')
        call leakage_check()
    case ('solovev', 'shaped')
        call ladder(pullback)
    case ('seed', 'seedshaped')
        call seed_ladder()
    case ('seeddiag')
        call seed_diagnosis()
    case ('gauge')
        call gauge_probe()
    case ('residual')
        call residual_check()
    case default
        error stop 'unknown mode'
    end select
contains
    subroutine make_mesh(n, k, state)
        integer, intent(in) :: n, k
        type(gs_state_t), intent(out) :: state
        real(dp), allocatable :: contour(:, :)
        real(dp) :: h
        integer :: j, status
        allocate (contour(2, n))
        do j = 1, n
            call boundary_point(real(j - 1, dp)/n, contour(:, j))
        end do
        h = 2*pi*a/n
        if (shaped) h = h*1.35_dp
        call initialize_gs(state, contour, h*h/2, status, k, curve=boundary_point)
        if (status /= 0) error stop 'mesh'
    end subroutine

    subroutine build_map(state, nfp, ntor, map, field)
        type(gs_state_t), intent(in) :: state
        class(ideal_frozen_t), intent(in), optional :: field
        integer, intent(in) :: nfp, ntor
        type(ideal_map_t), intent(out) :: map
        type(cylindrical_chart_t) :: chart
        integer :: status
        call ideal_map_init(map, state, chart, nfp, ntor, status, &
            quadrature_order=quadrature_order, axis_degree=axis_degree)
        if (status /= 0) error stop 'map init'
        if (map%sigma /= -1) error stop 'cylindrical orientation'
        map%axis_power = axis_power
        map%initial_damping = initial_damping
        if (gauge_weight > 0) map%gauge_weight = gauge_weight
        map%position_gauge = position_gauge
        map%gauge_hessian_only = gauge_hessian_only
        map%omega_elimination = omega_elimination
        if (present(field)) then
            call ideal_map_freeze(map, field, status)
        else
            call ideal_map_freeze(map, pullback, status)
        end if
        if (status /= 0) error stop 'freeze'
    end subroutine

    subroutine derivative_checks()
        ! Full-coefficient gradient against central differences of W; reduced
        ! Hessian against differences of the reduced gradient; the rigid
        ! toroidal rotation is an exact null direction of W.
        type(gs_state_t) :: state
        type(ideal_map_t) :: map
        type(csr_matrix_t) :: hess
        real(dp), allocatable :: x(:, :, :), v(:), vr(:), g(:), gp(:), gm(:), hv(:)
        real(dp), allocatable :: gr(:), grp(:), grm(:), xp(:, :, :)
        real(dp) :: e, ep, em, step, fd, an, scale
        integer :: status, i, ntor
        call make_mesh(16, 2, state)
        do ntor = 0, 1
            call build_map(state, 3, ntor, map)
            allocate (v(map%ndof), g(map%ndof), gp(map%ndof), gm(map%ndof), &
                vr(map%nred), gr(map%nred), grp(map%nred), grm(map%nred), &
                hv(map%nred))
            x = map%x
            do i = 1, map%nred
                vr(i) = sin(1.7_dp*i)*1.0e-2_dp
            end do
            where (map%fixed) vr = 0
            call ideal_map_expand(map, vr, v)
            ! Perturbed base point (all modes) so the n /= 0 chain is exercised.
            x = x + reshape(v, shape(x))
            do i = 1, map%ndof
                v(i) = cos(0.37_dp*i)
            end do
            call ideal_map_energy(map, x, e, status, g)
            if (status /= 0) error stop 'base energy'
            step = 1.0e-5_dp
            xp = x + step*reshape(v, shape(x))
            call ideal_map_energy(map, xp, ep, status)
            xp = x - step*reshape(v, shape(x))
            call ideal_map_energy(map, xp, em, status)
            fd = (ep - em)/(2*step); an = dot_product(g, v)
            print '(a,i1,3es14.5)', 'ntor gradient fd analytic rel ', ntor, fd, an, &
                abs(fd - an)/abs(an)
            if (abs(fd - an) > 1.0e-6_dp*abs(an)) error stop 'energy gradient'
            do i = 1, map%nred
                vr(i) = cos(0.37_dp*i)
            end do
            call ideal_map_expand(map, vr, v)
            xp = x + step*reshape(v, shape(x))
            call ideal_map_energy(map, xp, ep, status, gp)
            xp = x - step*reshape(v, shape(x))
            call ideal_map_energy(map, xp, em, status, gm)
            call ideal_map_reduce(map, gp, grp)
            call ideal_map_reduce(map, gm, grm)
            call ideal_map_hessian(map, x, hess)
            call multiply(hess, vr, hv)
            grp = (grp - grm)/(2*step)
            scale = norm2(hv)
            print '(a,i1,es14.5)', 'ntor reduced hessian rel ', ntor, &
                norm2(grp - hv)/scale
            if (norm2(grp - hv) > 1.0e-6_dp*scale) error stop 'energy hessian'
            vr = 0
            do i = 1, map%nodes
                vr(ideal_map_rdof(map, 2, 1, i)) = 1
            end do
            call ideal_map_reduce(map, g, gr)
            call multiply(hess, vr, hv)
            print '(a,i1,2es14.5)', 'ntor rotation null |g.v|, |Hv|/(max|H| |v|) ', &
                ntor, abs(dot_product(gr, vr))/(norm2(gr)*norm2(vr)), &
                norm2(hv)/(maxval(abs(hess%val))*norm2(vr))
            if (abs(dot_product(gr, vr)) > 1.0e-11_dp*norm2(gr)*norm2(vr)) &
                error stop 'rotation gradient null'
            if (norm2(hv) > 1.0e-10_dp*maxval(abs(hess%val))*norm2(vr)) &
                error stop 'rotation hessian null'
            deallocate (v, g, gp, gm, vr, gr, grp, grm, hv)
        end do
    end subroutine

    subroutine reduced_derivative_checks()
        ! Nested omega elimination: after it the omega gradient vanishes;
        ! the gradient of W_red(nu) = min_omega W against central differences
        ! of W_red (envelope property: the nu part of the reduced gradient);
        ! the Schur complement S = H_tt - H_to H_oo^-1 H_ot against
        ! differences of that gradient, S d = (H [d; y])_t with y = d omega*
        ! (differenced), and (H [d; y])_o = 0.
        type(gs_state_t) :: state
        type(ideal_map_t) :: map
        type(csr_matrix_t) :: hess
        real(dp), allocatable :: x(:, :, :), xp(:, :, :), xm(:, :, :), v(:)
        real(dp), allocatable :: vr(:), gr(:), grp(:), grm(:), z(:), hz(:)
        logical, allocatable :: omask(:)
        real(dp) :: e, ep, em, step, fd, an, scale
        integer :: status, i, ntor, a, m
        omega_elimination = .true.
        call make_mesh(16, 2, state)
        ! ntor > 0 is rejected: J depends on omega there.
        call build_map(state, 3, 1, map)
        x = map%x
        allocate (gr(map%nred))
        call ideal_map_eliminate_omega(map, x, e, gr, status)
        if (status /= IDEAL_MAP_BAD_INPUT) error stop 'ntor > 0 elimination'
        deallocate (gr)
        do ntor = 0, 0
            call build_map(state, 3, ntor, map)
            allocate (v(map%ndof), vr(map%nred), gr(map%nred), grp(map%nred), &
                grm(map%nred), z(map%nred), hz(map%nred), omask(map%nred))
            omask = .false.
            do a = 1, map%nodes
                do m = 1, map%nmode
                    omask(ideal_map_rdof(map, 2, m, a)) = .true.
                end do
            end do
            omask = omask .and. .not. map%fixed
            x = map%x
            do i = 1, map%nred
                vr(i) = sin(1.7_dp*i)*1.0e-2_dp
            end do
            where (map%fixed) vr = 0
            call ideal_map_expand(map, vr, v)
            x = x + reshape(v, shape(x))
            call ideal_map_eliminate_omega(map, x, e, gr, status)
            if (status /= 0) error stop 'elimination'
            print '(a,i1,es14.5)', 'ntor |g_omega|/|g| after elimination ', ntor, &
                norm2(pack(gr, omask))/norm2(gr)
            if (norm2(pack(gr, omask)) > 1.0e-8_dp*norm2(gr)) &
                error stop 'omega gradient after elimination'
            ! Transversal direction (omega part zero).
            do i = 1, map%nred
                vr(i) = cos(0.37_dp*i)
            end do
            where (map%fixed .or. omask) vr = 0
            call ideal_map_expand(map, vr, v)
            step = 1.0e-5_dp
            xp = x + step*reshape(v, shape(x))
            call ideal_map_eliminate_omega(map, xp, ep, grp, status)
            if (status /= 0) error stop 'elimination +'
            xm = x - step*reshape(v, shape(x))
            call ideal_map_eliminate_omega(map, xm, em, grm, status)
            if (status /= 0) error stop 'elimination -'
            fd = (ep - em)/(2*step); an = dot_product(gr, vr)
            print '(a,i1,3es14.5)', 'ntor W_red gradient fd analytic rel ', ntor, fd, &
                an, abs(fd - an)/abs(an)
            if (abs(fd - an) > 1.0e-6_dp*abs(an)) error stop 'reduced gradient'
            ! y = d omega*/d nu . d by differences; z = [d; y].
            z = vr
            do a = 1, map%nodes
                do m = 1, map%nmode
                    i = ideal_map_rdof(map, 2, m, a)
                    if (omask(i)) z(i) = (xp(3, m, a) - xm(3, m, a))/(2*step)
                end do
            end do
            call ideal_map_hessian(map, x, hess)
            call multiply(hess, z, hz)
            where (map%fixed) hz = 0
            grp = (grp - grm)/(2*step)
            where (omask) grp = 0
            scale = norm2(pack(hz, .not. omask))
            print '(a,i1,2es14.5)', 'ntor Schur fd rel, |(Hz)_o|/|(Hz)_t| ', ntor, &
                norm2(grp - merge(0.0_dp, hz, omask))/scale, &
                norm2(pack(hz, omask))/scale
            if (norm2(grp - merge(0.0_dp, hz, omask)) > 1.0e-5_dp*scale) &
                error stop 'Schur complement'
            ! Constrained omega (gauge moments): (Hz)_o = C lambda, zero after
            ! projection onto ker C^T.
            call ideal_map_gauge_project(map, hz)
            print '(a,es14.5)', 'ntor projected |(Hz)_o|/|(Hz)_t| ', &
                norm2(pack(hz, omask))/scale
            if (norm2(pack(hz, omask)) > 1.0e-5_dp*scale) error stop 'omega stationarity'
            deallocate (v, vr, gr, grp, grm, z, hz, omask)
        end do
    end subroutine

    subroutine multiply(matrix, v, y)
        type(csr_matrix_t), intent(in) :: matrix
        real(dp), intent(in) :: v(:)
        real(dp), intent(out) :: y(:)
        integer :: i, k
        y = 0
        do i = 1, matrix%n
            do k = matrix%rowptr(i), matrix%rowptr(i + 1) - 1
                y(i) = y(i) + matrix%val(k)*v(matrix%colidx(k))
            end do
        end do
    end subroutine

    subroutine solve_map(map, report, verbose, reference_gradient)
        type(ideal_map_t), intent(inout) :: map
        type(ideal_map_solve_report_t), intent(out) :: report
        logical, intent(in) :: verbose
        real(dp), intent(in), optional :: reference_gradient
        integer :: status
        call ideal_map_solve(map, solve_tolerance, max_newton, report, status, verbose, &
            reference_gradient)
        if (status /= 0) then
            print *, 'solve status', status, report%iterations, report%gradient
            error stop 'ideal map solve'
        end if
    end subroutine

    subroutine defects(map, state2d, field, errors)
        ! errors(1) ||psi_exact(X) - chi||, (2) ||psi_2D(X) - chi||,
        ! (3) seed ||psi_exact(y) - chi||, (4) ||psi_2D - psi_exact|| (2D alone),
        ! all over the reference with area weights, divided by |psi_axis|*sqrt(A).
        type(ideal_map_t), intent(in) :: map
        type(gs_state_t), intent(in) :: state2d
        class(ideal_frozen_t), intent(in) :: field
        real(dp), intent(out) :: errors(4)
        type(gauss_quadrature_triangle_t) :: rule
        real(dp) :: q(3), y(2), chi, ex, two, ex_y, g(2), area, w, b3, p
        real(dp) :: shape(10), dshape(2, 10), dg(2, 2), det
        integer :: e, k, status, owner
        call rule%init(7)
        errors = 0; area = 0; owner = 1
        do e = 1, map%cells
            do k = 1, rule%n_points
                call ideal_map_point(map, e, [rule%xi(k), rule%eta(k)], 0.0_dp, q, y)
                call triangle_basis_jet([rule%xi(k), rule%eta(k)], &
                    shape(:map%nloc), dshape(:, :map%nloc))
                dg = matmul(map%reference(:, map%conn(:, e)), &
                    transpose(dshape(:, :map%nloc)))
                det = dg(1, 1)*dg(2, 2) - dg(1, 2)*dg(2, 1)
                w = rule%weights(k)*det
                call field%evaluate(y, chi, g, b3, p)
                call psi_exact(q(1), q(2), ex, g)
                call psi_exact(y(1), y(2), ex_y, g)
                call locate_sample(state2d, q(1), q(2), owner, two, status)
                if (status /= 0) error stop '2D sample'
                errors(1) = errors(1) + w*(ex - chi)**2
                errors(2) = errors(2) + w*(two - chi)**2
                errors(3) = errors(3) + w*(ex_y - chi)**2
                errors(4) = errors(4) + w*(two - ex)**2
                area = area + w
            end do
        end do
        errors = sqrt(errors/area)/abs(psi_axis)
    end subroutine

    subroutine locate_sample(state, r, z, owner, value, status)
        ! psi_2D at (r, z): walk from the previous owner over all cells.
        type(gs_state_t), intent(in) :: state
        real(dp), intent(in) :: r, z
        integer, intent(inout) :: owner
        real(dp), intent(out) :: value
        integer, intent(out) :: status
        real(dp) :: grad(2), hessian(2, 2)
        integer :: k, e, n
        n = size(state%triangles, 2)
        do k = 0, n - 1
            e = 1 + mod(owner - 1 + merge(k/2, n - (k + 1)/2, mod(k, 2) == 0), n)
            call gs_element_jet(state, e, r, z, value, grad, hessian, status)
            if (status == 0) then
                owner = e; return
            end if
        end do
    end subroutine

    subroutine residual_check()
        ! Normalised weak residual (equilibrium-core D3) on the A3 pullback case:
        ! seed X0 = id, the interpolated exact map X = Y, and the discrete
        ! solution, on three meshes. sigma excludes the vacuum toroidal field.
        type(gs_state_t) :: state
        type(ideal_map_t) :: map
        type(ideal_map_solve_report_t) :: report
        real(dp) :: table(3, 3), sigma, x(2), dy(2, 2), bp2
        integer :: mesh, n, i, status
        print '(a)', 'n seed exact_interpolant solution sigma'
        do mesh = 1, 3
            n = 32*2**(mesh - 1)
            call make_mesh(n, degree, state)
            call build_map(state, 1, 0, map)
            call ideal_map_residual(map, f_edge(), table(1, mesh), sigma, status)
            if (status /= 0) error stop 'residual (seed)'
            do i = 1, map%nodes
                call seed_map_point(map%reference(:, i), x, dy)
                map%x(1:2, 1, i) = x
            end do
            call ideal_map_residual(map, f_edge(), table(2, mesh), sigma, status)
            if (status /= 0) error stop 'residual (exact)'
            call build_map(state, 1, 0, map)
            call solve_map(map, report, .false.)
            call ideal_map_residual(map, f_edge(), table(3, mesh), sigma, status)
            if (status /= 0) error stop 'residual (solution)'
            bp2 = sigma
            print '(i4,4es12.4)', n, table(:, mesh), bp2
        end do
        print '(a,3f8.3)', 'orders pair 1 ', log(table(:, 1)/table(:, 2))/log(2.0_dp)
        print '(a,3f8.3)', 'orders pair 2 ', log(table(:, 2)/table(:, 3))/log(2.0_dp)
        ! The discrete solution is stationary in the gauge space only; the
        ! tangential (relabelling) components remain at discretisation level.
        ! L2 force of the exact interpolant: order k - 1; of the solution: k.
        if (any(table(3, :) > table(2, :))) error stop 'solution above interpolant'
        if (log(table(2, 2)/table(2, 3))/log(2.0_dp) < degree - 1 - 0.3_dp) &
            error stop 'exact-interpolant residual order'
        if (log(table(3, 2)/table(3, 3))/log(2.0_dp) < degree - 0.3_dp) &
            error stop 'solution residual order'
    end subroutine

    real(dp) function reduced_gradient_norm(map) result(norm)
        ! Free reduced-gradient norm of the current map.
        type(ideal_map_t), intent(in) :: map
        real(dp), allocatable :: g(:), gr(:)
        real(dp) :: energy
        integer :: status
        allocate (g(map%ndof), gr(map%nred))
        call ideal_map_energy(map, map%x, energy, status, g)
        call ideal_map_reduce(map, g, gr)
        where (map%fixed) gr = 0
        norm = norm2(gr)
    end function

    subroutine seed_diagnosis()
        ! A3 only. Separates the map space from the solve: X* = P(xi(y)) with
        ! P(xi) = (sqrt(R0^2 + 2 R0 a xi1), R0 a xi2/R) is a smooth exact map
        ! for the nested-seed labels (psi_exact(X*) = chi). Per mesh: (1) its
        ! P_k nodal interpolant, unsolved; (2) the solve started from it;
        ! (3) the solve started from the seed X0 = id.
        type(nested_seed_t) :: seed
        type(solovev_profile_t) :: profile
        type(gs_state_t) :: state
        type(ideal_map_t) :: map
        type(ideal_map_solve_report_t) :: report
        real(dp) :: errors(4), table(3, 3), xi(2)
        real(dp), allocatable :: start(:, :, :)
        integer :: mesh, n, i, status
        logical :: ok
        call ideal_seed_init(seed, boundary_point, phi_edge(), profile, status)
        if (status /= 0) error stop 'seed init'
        solve_tolerance = 1.0e-8_dp
        print '(a)', 'n interpolant solve_from_exact solve_from_seed'
        do mesh = 1, 3
            n = 32*2**(mesh - 1)
            call make_mesh(n, degree, state)
            call solve_gs(state, 1.0_dp, pprime, ffprime, status)
            call build_map(state, 1, 0, map, seed)
            do i = 1, map%nodes
                if (map%boundary(i)) cycle
                call seed%inverse(map%reference(:, i), xi, ok)
                if (.not. ok) error stop 'seed inverse'
                map%x(1, 1, i) = sqrt(r0*r0 + 2*r0*a*xi(1))
                map%x(2, 1, i) = r0*a*xi(2)/map%x(1, 1, i)
            end do
            call defects(map, state, seed, errors)
            table(1, mesh) = errors(1)
            start = map%x
            call solve_map(map, report, verbose_solve)
            call defects(map, state, seed, errors)
            table(2, mesh) = errors(1)
            print '(a,3es12.4)', '  moved from X*: max |dq1,dq2| (m), max |domega|, h (m) ', &
                maxval(norm2(map%x(1:2, 1, :) - start(1:2, 1, :), dim=1)), &
                maxval(abs(map%x(3, 1, :) - start(3, 1, :))), 2*pi*a/n
            call build_map(state, 1, 0, map, seed)
            call solve_map(map, report, verbose_solve)
            call defects(map, state, seed, errors)
            table(3, mesh) = errors(1)
            print '(i4,3es12.4)', n, table(:, mesh)
        end do
        print '(a,3f8.3)', 'orders pair 1 ', log(table(:, 1)/table(:, 2))/log(2.0_dp)
        print '(a,3f8.3)', 'orders pair 2 ', log(table(:, 2)/table(:, 3))/log(2.0_dp)
    end subroutine

    subroutine gauge_probe()
        ! A3, nested. W_red curvature along probe directions at the exact
        ! map X* (as seed_diagnosis), relative to a physical axis shift
        ! (probe 6), and of W along omega = f(chi) (probes 7, 8).
        ! Relabellings y -> y + eps v (X* o (id + eps v), exact
        ! through the seed inverse): field-line slides v = g b_pol with
        ! g = T, T u1, T u2 (probes 1-3, continuum-exact nulls of W_red: the
        ! ratio must vanish with h), the tangential part of the r^2-weighted
        ! translation (4, = -u2 r t: a smooth slide only for circular
        ! surfaces) and of the plain translation (5, not a smooth slide).
        type(nested_seed_t) :: seed
        type(solovev_profile_t) :: profile
        type(gs_state_t) :: state
        type(ideal_map_t) :: map
        real(dp) :: xi(2), chi, grad(2), b3, p, t(2), e0, ep, eps
        real(dp) :: kappa(8, 2), d2(8), s(2), floor, w(-1:1)
        real(dp), allocatable :: v(:, :), x0(:, :, :), x(:, :, :), gr(:)
        integer :: mesh, n, i, k, status, sgn, nk
        logical :: ok
        call ideal_seed_init(seed, boundary_point, phi_edge(), profile, status)
        if (status /= 0) error stop 'seed init'
        omega_elimination = .true.
        print '(a)', 'n  W_red curvature / axis-shift curvature, probes 1-6, floor'
        do mesh = 1, 2
            n = merge(24, 32, degree == 3)*2**(mesh - 1)
            call make_mesh(n, degree, state)
            call build_map(state, 1, 0, map, seed)
            do i = 1, map%nodes
                if (map%boundary(i)) cycle
                call seed%inverse(map%reference(:, i), xi, ok)
                if (.not. ok) error stop 'seed inverse'
                map%x(1:2, 1, i) = exact(xi)
            end do
            x0 = map%x
            allocate (v(2, map%nodes), gr(map%nred))
            call ideal_map_eliminate_omega(map, x0, e0, gr, status)
            if (status /= 0) error stop 'elimination'
            do k = 1, 6
                ! The quadratic axis monomials for probe 4 exist only at P3.
                ! P2 already reports this inapplicable probe as zero below.
                if (k == 4 .and. map%npoly < 6) then
                    d2(k) = 0
                    cycle
                end if
                do i = 1, map%nodes
                    call seed%evaluate(map%reference(:, i), chi, grad, b3, p)
                    t = [-map%direction(2, i), map%direction(1, i)]
                    select case (k)
                    case (1:3)
                        v(:, i) = map%axis_basis(k, i)*[-grad(2), grad(1)]
                    case (4)
                        v(:, i) = t*t(1)*(map%axis_basis(4, i) + map%axis_basis(6, i))
                    case (5)
                        v(:, i) = t*t(1)*map%axis_basis(1, i)
                    case (6)
                        v(:, i) = [map%axis_basis(1, i), 0.0_dp]
                    end select
                end do
                ! Second differences at eps and 2 eps (amplitude 1e-3 a),
                ! Richardson-combined: the linear relabelling id + eps v differs
                ! from the slide flow at O(eps^2), which the combination removes.
                eps = 1.0e-3_dp*a/maxval(norm2(v, dim=1))
                s = -2*e0
                do sgn = -2, 2
                    if (sgn == 0) cycle
                    x = x0
                    do i = 1, map%nodes
                        if (map%boundary(i)) cycle
                        if (k == 6) then
                            x(1:2, 1, i) = x0(1:2, 1, i) + sgn*eps*v(:, i)
                        else
                            call seed%inverse(map%reference(:, i) + sgn*eps*v(:, i), &
                                xi, ok)
                            if (.not. ok) error stop 'seed inverse (probe)'
                            x(1:2, 1, i) = exact(xi)
                        end if
                    end do
                    call ideal_map_eliminate_omega(map, x, ep, gr, status)
                    if (status /= 0) error stop 'elimination (probe)'
                    s(abs(sgn)) = s(abs(sgn)) + ep
                end do
                d2(k) = (4*s(1) - s(2)/4)/3
            end do
            ! Omega gauge f(chi) (probes 7, 8: f = T, sin(2 pi T)) at fixed
            ! transversal dofs: W is quadratic in omega, so the second
            ! difference is the exact curvature f.H_oo f.
            do k = 7, 8
                do sgn = -1, 1
                    x = x0
                    do i = 1, map%nodes
                        eps = map%axis_t(i)
                        if (k == 8) eps = sin(2*pi*eps)
                        x(3, 1, i) = x0(3, 1, i) + sgn*1.0e-3_dp*eps
                    end do
                    call ideal_map_energy(map, x, w(sgn), status)
                    if (status /= 0) error stop 'energy (omega probe)'
                end do
                ! Same amplitude scale as the axis shift: |eps v| = 1e-3 a.
                d2(k) = (w(1) + w(-1) - 2*w(0))*a**2
            end do
            kappa(:, mesh) = d2/d2(6)
            ! Roundoff floor of the differences (W_red to ~64 ulp).
            floor = 64*epsilon(1.0_dp)*abs(e0)/abs(d2(6))
            print '(i4,9es11.3)', n, kappa(:, mesh), floor
            deallocate (v, gr)
        end do
        d2 = 0
        do k = 1, 8
            if (abs(kappa(k, 2)) > 0) &
                d2(k) = log(abs(kappa(k, 1)/kappa(k, 2)))/log(2.0_dp)
        end do
        print '(a,8f8.2)', 'orders ', d2
        ! Field-line slides (1-3) and, at P3, the tangential part of the
        ! r^2-weighted translation (4, = -u2 r t, a slide of the same kind):
        ! curvature far below the physical one (either sign) and decreasing
        ! at least at order 1 with h (measured 2.7-6.8). The tangential
        ! translation (5) is not a smooth slide: O(1).
        nk = merge(4, 3, degree >= 3)
        if (any(abs(kappa(1:nk, 2)) > 1.0e-3_dp)) error stop 'slide not near-null'
        if (any(abs(kappa(1:nk, 2)) > max(abs(kappa(1:nk, 1))/2, floor))) &
            error stop 'slide curvature not vanishing with h'
        if (abs(kappa(5, 2)) < 0.1_dp) &
            error stop 'translation slide unexpectedly null'
        ! Omega gauge f(chi) (7: f = T, near-null; 8: f = sin(2 pi T), one
        ! radial wave): curvature vanishing with h (measured orders 4.3, 2.8
        ! at P2; at P3 probe 7 reaches the roundoff floor on the finer mesh).
        if (abs(kappa(7, 2)) > 1.0e-3_dp) error stop 'omega gauge not near-null'
        if (any(abs(kappa(7:8, 2)) > max(abs(kappa(7:8, 1))/2, floor))) &
            error stop 'omega gauge curvature not vanishing with h'
    end subroutine

    function exact(xi) result(q)
        ! Smooth exact A3 map P(xi) of the nested-seed labels (seed_diagnosis).
        real(dp), intent(in) :: xi(2)
        real(dp) :: q(2)
        q(1) = sqrt(r0*r0 + 2*r0*a*xi(1))
        q(2) = r0*a*xi(2)/q(1)
    end function

    subroutine seed_ladder()
        type(nested_seed_t) :: seed
        type(solovev_profile_t) :: profile
        integer :: status
        ! Seed axis: the harmonic-extension centre Y(0) (boundary mean in its
        ! parameter); the magnetic axis is elsewhere and must be found.
        call ideal_seed_init(seed, boundary_point, phi_edge(), profile, status)
        print '(a,2f9.4)', 'seed axis offset from magnetic axis (m) ', &
            real(seed%a(0)) - raxis, aimag(seed%a(0))
        if (status /= 0) error stop 'seed init'
        if (abs(seed%psi_axis/psi_axis - 1) > 1.0e-12_dp) error stop 'seed axis flux'
        ! Stop at seed_tolerance x |g(X0 = id)| on each mesh. P3 uses the
        ! stronger transversal gauge weight 1e-2, for which the end point is
        ! independent of the initial damping and the factoriser.
        solve_tolerance = seed_tolerance
        if (gauge_weight <= 0 .and. degree >= 3) gauge_weight = 1.0e-2_dp
        call ladder(seed, sequenced=.true.)
    end subroutine

    subroutine ladder(field, sequenced)
        ! sequenced: each finer map starts from the previous solution (grid
        ! sequencing); otherwise every mesh starts from the seed X0 = id.
        class(ideal_frozen_t), intent(in) :: field
        logical, intent(in), optional :: sequenced
        type(map_case_state_t) :: result
        type(equilibrium_run_report_t) :: verification
        real(dp) :: errors(4, 3), rates(4)
        integer :: mesh, n, base
        base = 32
        if (degree == 3 .and. .not. shaped) base = 24
        allocate (result%field, source=field)
        if (present(sequenced)) result%sequenced = sequenced
        print '(a)', 'case k n dof newton wall_s d3D_exact d3D_2D seed d2D_exact'
        do mesh = 1, 3
            n = base*2**(mesh - 1)
            call run_equilibrium_case(trim(mode), 'map', &
                equilibrium_resolution_t(degree, n), 'legacy', run_map_case, &
                result, verification)
            call write_equilibrium_report(output_unit, verification)
            if (verification%status /= 0) error stop 'D8 map case'
            errors(:, mesh) = result%errors
            print '(a,3es24.16)', '  end point: d3D_exact, |x|_2, W ', &
                errors(1, mesh), norm2(result%map%x), result%solve%energy
            print '(a8,i2,i5,i8,i4,f9.2,4es12.4)', trim(mode), degree, n, &
                count(.not. result%map%fixed), result%solve%iterations, &
                result%elapsed, errors(:, mesh)
            print '(a,5f9.2)', &
                '  wall energy hessian factor+solve (analysis numeric) ', &
                result%solve%energy_s, result%solve%hessian_s, &
                result%solve%factor_s, result%solve%analysis_s, result%solve%numeric_s
        end do
        do mesh = 2, 3
            rates = log(errors(:, mesh - 1)/errors(:, mesh))/log(2.0_dp)
            print '(a,i1,4f8.3)', 'orders pair ', mesh - 1, rates
        end do
        rates = log(errors(:, 2)/errors(:, 3))/log(2.0_dp)
        if (present(sequenced)) then
            ! Nested seed: return from the offset seed axis (and the k + 1 order
            ! gates below).
            if (errors(1, 3) > 1.0e-3_dp*errors(3, 3)) error stop 'seed not returned'
        end if
        if (rates(1) < degree + 1 - 0.4_dp) error stop '3D N=0 exact-flux order'
        if (rates(2) < degree + 1 - 0.4_dp) error stop '3D N=0 vs 2D order'
        if (errors(1, 3) > 1.0e-2_dp*errors(3, 3)) error stop 'seed not corrected'
    end subroutine

    subroutine run_map_case(spec, resolution, start, state, verification)
        type(equilibrium_case_t), intent(in) :: spec
        type(equilibrium_resolution_t), intent(in) :: resolution
        type(equilibrium_start_t), intent(in) :: start
        class(*), intent(inout) :: state
        type(equilibrium_run_report_t), intent(inout) :: verification
        real(dp) :: seed_gradient
        integer :: status
        integer(int64) :: begin, finish, rate
        if (spec%lane /= 'map' .or. start%name /= 'legacy') return
        select type (state)
        type is (map_case_state_t)
            call make_mesh(resolution%boundary_nodes, resolution%degree, state%gs)
            call system_clock(begin, rate)
            call build_map(state%gs, 1, 0, state%map, state%field)
            if (state%sequenced) then
                seed_gradient = reduced_gradient_norm(state%map)
                if (allocated(state%previous)) then
                    call ideal_map_prolong(state%map, state%previous, status)
                    if (status /= 0) then
                        verification%status = status
                        return
                    end if
                end if
                call solve_map(state%map, state%solve, verbose_solve, seed_gradient)
            else
                call solve_map(state%map, state%solve, verbose_solve)
            end if
            state%previous = state%map
            call system_clock(finish)
            state%elapsed = real(finish - begin, dp)/rate
            call solve_gs(state%gs, 1.0_dp, pprime, ffprime, status)
            if (status /= 0) then
                verification%status = status
                return
            end if
            call defects(state%map, state%gs, state%field, state%errors)
            verification%oracle_error = state%errors(1)
            verification%seed_error = state%errors(3)
            verification%iterations = state%solve%iterations
            call ideal_map_residual(state%map, f_edge(), verification%residual, &
                verification%sigma, verification%status)
            verification%residual_norm = 'D3: L2 force * L / (sigma sqrt(V))'
        end select
    end subroutine

    subroutine leakage_check()
        ! N = 1, nfp = 3: an n = 3 interior seed (axisymmetric boundary) must
        ! decay. n /= 0 map coefficients are not gauge invariant (relabelling
        ! along B), so leakage is measured on the surfaces: the phi-variation of
        ! the defect psi_exact(X(y, phi)) - chi(y) over one period. Relabelling
        ! changes the discretisation error pattern, so the variation must fall
        ! to the axisymmetric defect level and converge with it.
        type(gs_state_t) :: state
        type(ideal_map_t) :: map, axisym
        type(ideal_map_solve_report_t) :: report, axisym_report
        real(dp) :: seeded, final(2), defect(2), difference, v, g(2), rate
        integer :: i, mesh
        do mesh = 1, 2
            call make_mesh(16*mesh, 2, state)
            call build_map(state, 1, 0, axisym)
            call solve_map(axisym, axisym_report, .false.)
            call toroidal_variation(axisym, v, defect(mesh))
            call build_map(state, 3, 1, map)
            do i = 1, map%nodes
                if (map%boundary(i)) cycle
                call psi_exact(map%reference(1, i), map%reference(2, i), v, g)
                map%x(1, 2, i) = 0.04_dp*a*v/psi_axis
                map%x(2, 3, i) = 0.03_dp*a*v/psi_axis
                map%x(3, 2, i) = 0.02_dp*v/psi_axis
            end do
            call toroidal_variation(map, seeded, v)
            call solve_map(map, report, verbose_solve)
            call toroidal_variation(map, final(mesh), v)
            difference = abs(3*report%energy - axisym_report%energy) &
                /abs(axisym_report%energy)
            print '(a,i3,4es12.4)', 'n seeded final N=0-defect dW ', 16*mesh, &
                seeded, final(mesh), defect(mesh), difference
            if (final(mesh) > defect(mesh)) error stop 'n=3 above axisymmetric defect'
            if (difference > 1.0e-8_dp) error stop 'energy differs from N=0 solve'
        end do
        rate = log(final(1)/final(2))/log(2.0_dp)
        print '(a,f8.3)', 'n=3 surface variation order ', rate
        ! Gauge relabelling (near-null) sets the residual pattern, so only
        ! convergence is required here, not the k+1 order (measured 1.6).
        if (rate < 1) error stop 'n=3 variation does not converge'
    end subroutine

    subroutine toroidal_variation(map, variation, defect)
        ! RMS over the reference of the phi-variation (8 angles per period) of
        ! d = psi_exact(X(y, phi)) - chi(y), and RMS of d; both / |psi_axis|.
        type(ideal_map_t), intent(in) :: map
        real(dp), intent(out) :: variation, defect
        type(gauss_quadrature_triangle_t) :: rule
        real(dp) :: q(3), y(2), ex, g(2), d(8), area, w
        integer :: e, k, j
        call rule%init(7)
        variation = 0; defect = 0; area = 0
        do e = 1, map%cells
            do k = 1, rule%n_points
                do j = 1, 8
                    call ideal_map_point(map, e, [rule%xi(k), rule%eta(k)], &
                        2*pi*(j - 1)/(8*map%nfp), q, y)
                    call psi_exact(q(1), q(2), ex, g)
                    d(j) = ex - label(y)
                end do
                w = rule%weights(k)
                variation = variation + w*sum((d - sum(d)/8)**2)/8
                defect = defect + w*sum(d**2)/8
                area = area + w
            end do
        end do
        variation = sqrt(variation/area)/abs(psi_axis)
        defect = sqrt(defect/area)/abs(psi_axis)
    end subroutine
end program
