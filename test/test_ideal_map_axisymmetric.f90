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
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use kin6d_grad_shafranov, only: gs_state_t, initialize_gs, solve_gs, &
        gs_element_jet
    use kin6d_ideal_map, only: ideal_map_t, cylindrical_chart_t, ideal_map_init, &
        ideal_map_freeze, ideal_map_energy, ideal_map_hessian, ideal_map_solve, &
        ideal_map_point, ideal_map_rdof, ideal_map_solve_report_t, &
        ideal_map_reduce, ideal_map_expand, ideal_map_prolong
    use fortfem_element_to_csr, only: csr_matrix_t
    use kin6d_triangle_basis, only: triangle_basis_jet
    use fortfem_gauss_quadrature_2d, only: gauss_quadrature_triangle_t
    use kin6d_ideal_map, only: ideal_frozen_t
    use kin6d_ideal_seed, only: nested_seed_t, ideal_seed_init
    use ideal_map_cases, only: setup_case, psi_exact, boundary_point, label, &
        phi_edge, pullback_field_t, solovev_profile_t, pi, a, raxis, psi_axis, &
        pprime, ffprime, shaped
    implicit none
    type(pullback_field_t) :: pullback
    character(16) :: mode, argument
    integer :: degree
    logical :: verbose_solve = .false.
    real(dp) :: solve_tolerance = 1.0e-6_dp
    mode = 'solovev'; degree = 2
    if (command_argument_count() > 0) call get_command_argument(1, mode)
    if (command_argument_count() > 1) then
        call get_command_argument(2, argument)
        read (argument, *) degree
    end if
    if (command_argument_count() > 2) verbose_solve = .true.
    call setup_case(trim(mode) == 'shaped' .or. trim(mode) == 'seedshaped')
    select case (trim(mode))
    case ('derivatives')
        call derivative_checks()
    case ('leakage')
        call leakage_check()
    case ('solovev', 'shaped')
        call ladder(pullback)
    case ('seed', 'seedshaped')
        call seed_ladder()
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
        call ideal_map_init(map, state, chart, nfp, ntor, status)
        if (status /= 0) error stop 'map init'
        if (map%sigma /= -1) error stop 'cylindrical orientation'
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

    subroutine solve_map(map, report, verbose)
        type(ideal_map_t), intent(inout) :: map
        type(ideal_map_solve_report_t), intent(out) :: report
        logical, intent(in) :: verbose
        integer :: status
        call ideal_map_solve(map, solve_tolerance, 100, report, status, verbose)
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
        ! The seed is far from equilibrium (large initial gradient), so the
        ! relative stopping tolerance is tightened accordingly.
        solve_tolerance = 1.0e-9_dp
        call ladder(seed, sequenced=.true.)
    end subroutine

    subroutine ladder(field, sequenced)
        ! sequenced: each finer map starts from the previous solution (grid
        ! sequencing); otherwise every mesh starts from the seed X0 = id.
        class(ideal_frozen_t), intent(in) :: field
        logical, intent(in), optional :: sequenced
        type(gs_state_t) :: state
        type(ideal_map_t) :: map, previous
        type(ideal_map_solve_report_t) :: report
        real(dp) :: errors(4, 3), rates(4), elapsed
        integer :: mesh, n, status, base
        integer(int64) :: start, finish, rate
        base = 32
        if (degree == 3 .and. .not. shaped) base = 24
        print '(a)', 'case k n dof newton wall_s d3D_exact d3D_2D seed d2D_exact'
        do mesh = 1, 3
            n = base*2**(mesh - 1)
            call make_mesh(n, degree, state)
            call system_clock(start, rate)
            call build_map(state, 1, 0, map, field)
            if (present(sequenced) .and. mesh > 1) then
                if (sequenced) then
                    call ideal_map_prolong(map, previous, status)
                    if (status /= 0) error stop 'prolong'
                end if
            end if
            call solve_map(map, report, verbose_solve)
            previous = map
            call system_clock(finish)
            elapsed = real(finish - start, dp)/rate
            call solve_gs(state, 1.0_dp, pprime, ffprime, status)
            if (status /= 0) error stop '2D GS solve'
            call defects(map, state, field, errors(:, mesh))
            print '(a8,i2,i5,i8,i4,f9.2,4es12.4)', trim(mode), degree, n, &
                count(.not. map%fixed), report%iterations, elapsed, errors(:, mesh)
            print '(a,3f9.2)', '  wall energy hessian factor+solve ', &
                report%energy_s, report%hessian_s, report%factor_s
        end do
        do mesh = 2, 3
            rates = log(errors(:, mesh - 1)/errors(:, mesh))/log(2.0_dp)
            print '(a,i1,4f8.3)', 'orders pair ', mesh - 1, rates
        end do
        rates = log(errors(:, 2)/errors(:, 3))/log(2.0_dp)
        if (rates(1) < degree + 1 - 0.4_dp) error stop '3D N=0 exact-flux order'
        if (rates(2) < degree + 1 - 0.4_dp) error stop '3D N=0 vs 2D order'
        if (errors(1, 3) > 1.0e-2_dp*errors(3, 3)) error stop 'seed not corrected'
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
