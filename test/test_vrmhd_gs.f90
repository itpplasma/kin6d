program test_vrmhd_gs
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_gs_symmetry, only: gs_symmetry_t, GS_ROTATION, GS_SCREW, &
        GS_CHART_RZ, GS_CHART_CARTESIAN
    use kin6d_grad_shafranov, only: gs_state_t, initialize_gs, gs_cell_data, &
        gs_element_jet, gs_linear_response
    use kin6d_vrmhd_parameters, only: vrmhd_parameters_t
    use kin6d_vrmhd_section, only: vrmhd_section_t, vrmhd_section_build, &
        vrmhd_section_basis
    use kin6d_vrmhd_stationary, only: vrmhd_stationary_t, vrmhd_stationary_build, &
        vrmhd_stationary_residual
    use kin6d_vrmhd_boundary, only: vrmhd_ampere_boundary
    use kin6d_vrmhd_drive, only: vrmhd_drive_build, VRMHD_CURRENT_SOURCE, &
        VRMHD_LOOP_VOLTAGE
    use kin6d_vrmhd_gs_leaves, only: gs_lift_rotation, gs_lift_screw, &
        vrmhd_gs_flux_rotation, vrmhd_gs_flux_screw, &
        vrmhd_gs_heat_rotation, vrmhd_gs_heat_screw
    use fortfem_mesh_2d, only: mesh_2d_t
    use fortfem_meridian_curved_mesh, only: meridian_mesh_t, meridian_mesh_straight
    use fortfem_triangle_duffy_quadrature, only: triangle_duffy_quadrature
    use fortsparse, only: csc_t, csc_from_triplet, sparse_solve_once, &
        fortsparse_status_t, status_ok
    implicit none
    integer :: kind, convention, level, trace_kind
    real(dp) :: defect(3), error(3), rate
    do convention = 0, 2
        do kind = GS_ROTATION, GS_SCREW, 2
            do level = 1, 3
                call consistency(kind, convention, 2**(level + 1), &
                    defect(level), error(level))
                if (level == 1) cycle
                rate = log(error(level - 1)/error(level))/log(2.0_dp)
                print '(a,3i4,f10.5)', 'GS lifted field convergence rate ', &
                    kind, convention, level, rate
                call require(rate > 0.8_dp, 'GS lifted P1 design order')
                rate = log(defect(level - 1)/defect(level))/log(2.0_dp)
                print '(a,3i4,f10.5)', 'GS stationary dual residual rate ', &
                    kind, convention, level, rate
                call require(rate > 0.8_dp, 'GS stationary residual design order')
            end do
        end do
    end do
contains
    subroutine require(ok, message)
        logical, intent(in) :: ok
        character(*), intent(in) :: message
        if (ok) return
        print '(a)', message
        error stop 'GS stationary consistency'
    end subroutine

    subroutine flux(kind, x, psi, grad, source)
        integer, intent(in) :: kind
        real(dp), intent(in) :: x(2)
        real(dp), intent(out) :: psi, grad(2), source
        if (kind == GS_ROTATION) then
            call vrmhd_gs_flux_rotation(x(1), x(2), 0.0_dp, psi, &
                grad(1), grad(2), source)
        else
            call vrmhd_gs_flux_screw(x(1), x(2), 0.7_dp, psi, &
                grad(1), grad(2), source)
        end if
    end subroutine

    subroutine lift(kind, x, psi, grad, values)
        integer, intent(in) :: kind
        real(dp), intent(in) :: x(2), psi, grad(2)
        real(dp), intent(out) :: values(38)
        if (kind == GS_ROTATION) then
            call gs_lift_rotation(x(1), x(2), 0.0_dp, psi, grad(1), grad(2), values)
        else
            call gs_lift_screw(x(1), x(2), 0.7_dp, psi, grad(1), grad(2), values)
        end if
    end subroutine

    subroutine oracle(kind, x, values)
        integer, intent(in) :: kind
        real(dp), intent(in) :: x(2)
        real(dp), intent(out) :: values(38)
        real(dp) :: r, z, psi
        r = x(1); z = x(2); values = 0
        values(11:12) = 1
        if (kind == GS_ROTATION) then
            psi = 0.1_dp*r*r*(log(r) - 0.5_dp)
            values(1) = 3 - 0.01_dp*log(r)**2
            values(6) = psi/r
            values(9) = -0.2_dp/r
            values(13) = -0.02_dp*log(r)/r
            values(27) = 0.2_dp*log(r)
            values(37) = psi; values(38) = -0.2_dp/r
        else
            values(1) = 3 - 0.02_dp*(r*r + z*z)
            values(5:7) = [-z/2, r/2, -0.1_dp*(r*r + z*z)]
            values(10) = 0.4_dp
            values(13:14) = -0.04_dp*x
            values(25:27) = [-0.2_dp*z, 0.2_dp*r, 1.0_dp]
            values(37) = 0.25_dp*(r*r + z*z)
            values(38) = -1/(1 + 0.49_dp*(r*r + z*z))**2
        end if
    end subroutine

    subroutine magnetic_trace(x, b)
        real(dp), intent(in) :: x(2)
        real(dp), intent(out) :: b(3)
        real(dp) :: exact(38)
        call oracle(trace_kind, x, exact)
        b = exact(25:27)
    end subroutine

    subroutine equilibrium(kind, n, gs, mm)
        integer, intent(in) :: kind, n
        type(gs_state_t), intent(out) :: gs
        type(meridian_mesh_t), intent(out) :: mm
        type(gs_symmetry_t) :: symmetry
        type(mesh_2d_t) :: mesh
        real(dp), allocatable :: rhs(:), solution(:), xi(:), eta(:), qw(:)
        real(dp) :: contour(2, 4), shift, psi, grad(2), source, exact(38), values(38)
        real(dp) :: x(2), det, shape(6), gradients(2, 6), hessians(2, 2, 6)
        integer :: i, j, e, q, status, iterations
        shift = 0
        symmetry = gs_symmetry_t(kind, GS_CHART_CARTESIAN, 0.7_dp)
        if (kind == GS_ROTATION) then
            shift = 2; symmetry = gs_symmetry_t(kind, GS_CHART_RZ, 0.0_dp)
        end if
        contour(:, 1) = [shift + 0.1_dp, 0.2_dp]
        contour(:, 2) = [shift + 1.1_dp, 0.2_dp]
        contour(:, 3) = [shift + 1.1_dp, 1.2_dp]
        contour(:, 4) = [shift + 0.1_dp, 1.2_dp]
        call initialize_gs(gs, contour, 0.5_dp/n**2, status, &
            element_order=2, symmetry=symmetry)
        call require(status == 0, 'RM0 equilibrium mesh')
        allocate (rhs(size(gs%psi)))
        gs%psi = 0
        do i = 1, size(gs%psi)
            call flux(kind, gs%points(:, i), psi, grad, source)
            call oracle(kind, gs%points(:, i), exact)
            call require(abs(psi - exact(37)) < 1.0e-13_dp, 'independent GS flux')
            call lift(kind, gs%points(:, i), psi, grad, values)
            call require(maxval(abs(values - exact)) < 1.0e-13_dp, &
                'independent Cartesian GS equilibrium and lift')
            if (gs%boundary(i)) gs%psi(i) = psi
        end do
        call triangle_duffy_quadrature(10, xi, eta, qw, status)
        call require(status == 0, 'GS pressure-source quadrature')
        rhs = 0
        do e = 1, size(gs%triangles, 2)
            do q = 1, size(qw)
                call gs_cell_data(gs, e, [xi(q), eta(q)], x, det, shape, &
                    gradients, hessians, status)
                call require(status == 0, 'RM0 source geometry')
                call flux(kind, x, psi, grad, source)
                do i = 1, 6
                    j = gs%triangles(i, e)
                    rhs(j) = rhs(j) + det*qw(q)*source*shape(i)
                end do
            end do
            do j = 1, 6
                if (.not. gs%boundary(gs%triangles(j, e))) cycle
                do i = 1, 6
                    rhs(gs%triangles(i, e)) = rhs(gs%triangles(i, e)) &
                        - gs%stiffness(i, j, e)*gs%psi(gs%triangles(j, e))
                end do
            end do
        end do
        where (gs%boundary) rhs = gs%psi
        ! The spatial source is the analytic GS profile evaluated on its monotone flux.
        call gs_linear_response(gs, rhs, solution, status, iterations, &
            relative_tolerance=1.0e-12_dp, residual_norm=gs%residual)
        call require(status == 0, 'RM0 prescribed-profile equilibrium solve')
        call require(gs%residual < 1.0e-10_dp, 'RM0 original-operator residual')
        gs%psi = solution
        mesh%n_vertices = maxval(gs%triangles(1:3, :))
        mesh%n_triangles = size(gs%triangles, 2); mesh%has_triangles = .true.
        mesh%vertices = gs%points(:, :mesh%n_vertices)
        mesh%triangles = gs%triangles(1:3, :)
        call meridian_mesh_straight(mesh, 1, mm, status)
        call require(status == 0, 'GS to invariant-section mesh')
    end subroutine

    subroutine consistency(kind, convention, n, defect, error)
        integer, intent(in) :: kind, convention, n
        real(dp), intent(out) :: defect, error
        type(gs_state_t) :: gs
        type(meridian_mesh_t) :: mm
        type(vrmhd_section_t) :: op
        type(vrmhd_stationary_t) :: steady
        type(vrmhd_parameters_t) :: p
        type(csc_t) :: gram, dual_gram, current_gram
        type(fortsparse_status_t) :: st
        integer, allocatable :: rows(:), cols(:), index(:)
        real(dp), allocatable :: vals(:), rhs(:), u(:), load(:), surface(:), r(:)
        real(dp), allocatable :: zero(:), unused(:), basis(:, :), metric(:, :)
        real(dp), allocatable :: source(:, :, :), js(:, :, :), h(:), dual(:), recovered(:)
        real(dp) :: x(2), weight, psi, grad(2), hessian(2, 2), values(38), exact(38)
        real(dp) :: field(36), errors(3), parts(7), voltage, gs_error
        real(dp), allocatable :: local(:, :)
        integer :: t, q, i, j, row, col, pos, status, nt, nl, ng, nc, entry
        call equilibrium(kind, n, gs, mm)
        call vrmhd_section_build(mm, gs%symmetry, 1, 1, op, status)
        call require(status == 0, 'GS section')
        call vrmhd_stationary_build(op, steady, status)
        call require(status == 0, 'GS Coulomb gauge')
        p%mi = 1; p%mu0 = 1
        p%eta0 = 0
        if (convention > 0) p%eta0 = 0.2_dp
        nt = mm%topology%n_triangles; nl = op%nlocal
        nc = 6*op%fs%n0 + op%fs%sp(1)%ndof
        allocate (basis(36, nl), metric(36, nl), local(nl, nl), rhs(op%ndof), &
            rows(nt*(nl**2 + op%fs%nq*op%fs%nl0**2)), &
            cols(nt*(nl**2 + op%fs%nq*op%fs%nl0**2)), &
            vals(nt*(nl**2 + op%fs%nq*op%fs%nl0**2)), u(steady%ndof), &
            load(op%ndof), surface(op%ndof), r(steady%ndof), &
            zero(steady%ndof), unused(steady%ndof), source(6, op%fs%nq, nt), &
            js(3, op%fs%nq, nt), h(op%fs%sp(1)%ndof), index(op%fs%n0), &
            recovered(op%ndof), dual(steady%ndof))
        rhs = 0; pos = 0; u = 0; source = 0; js = 0; gs_error = 0
        do t = 1, nt
            local = 0
            do q = 1, op%fs%nq
                call vrmhd_section_basis(op, t, q, x, weight, basis, status)
                call require(status == 0, 'GS projection basis')
                call gs_element_jet(gs, t, x(1), x(2), psi, grad, hessian, status)
                call require(status == 0, 'solved RM0 flux readback')
                call lift(kind, x, psi, grad, values)
                call oracle(kind, x, exact)
                gs_error = gs_error + weight*(psi - exact(37))**2
                metric = basis; metric(34:36, :) = 0
                local = local + weight*matmul(transpose(metric), metric)
                do i = 1, nl
                    row = op%map(i, t)
                    rhs(row) = rhs(row) + weight*dot_product(metric(:, i), values(:36))
                end do
                js(:, q, t) = exact(8:10)
                if (convention == 2) then
                    ! Analytic cooling balances Ohmic heat left by the voltage drive.
                    if (kind == GS_ROTATION) then
                        call vrmhd_gs_heat_rotation(x(1), p%eta0, source(5, q, t))
                    else
                        call vrmhd_gs_heat_screw(p%eta0, source(5, q, t))
                    end if
                    call require(abs(source(5, q, t) &
                        + p%eta0*dot_product(exact(8:10), exact(8:10))) < 1.0e-14_dp, &
                        'independent analytic GS heat sink')
                end if
            end do
            do j = 1, nl
                do i = 1, nl
                    pos = pos + 1
                    rows(pos) = op%map(i, t); cols(pos) = op%map(j, t)
                    vals(pos) = local(i, j)
                end do
            end do
        end do
        call csc_from_triplet(op%ndof, op%ndof, rows(:pos), cols(:pos), &
            vals(:pos), gram, st)
        call require(status_ok(st), 'GS physical projection Gram matrix')
        call sparse_solve_once(gram, rhs, recovered, st)
        call require(status_ok(st), 'GS physical projection')
        u(:op%ndof) = recovered
        zero = 0; trace_kind = kind
        call vrmhd_ampere_boundary(op, magnetic_trace, surface, status)
        call require(status == 0, 'GS Ampere trace')
        u(nc + 1:op%ndof) = 0
        call vrmhd_stationary_residual(op, steady, p, u, zero, u(:op%ndof), &
            source, -surface, r, unused, status)
        call require(status == 0, 'GS weak-current recovery residual')
        pos = 0
        do j = nc + 1, op%ndof
            do entry = gram%col_ptr(j), gram%col_ptr(j + 1) - 1
                row = gram%row_idx(entry)
                if (row <= nc) cycle
                pos = pos + 1; rows(pos) = row - nc; cols(pos) = j - nc
                vals(pos) = gram%val(entry)
            end do
        end do
        call csc_from_triplet(op%ndof - nc, op%ndof - nc, rows(:pos), &
            cols(:pos), vals(:pos), current_gram, st)
        call require(status_ok(st), 'GS current L2 Riesz matrix')
        deallocate (recovered); allocate (recovered(op%ndof - nc))
        call sparse_solve_once(current_gram, -r(nc + 1:op%ndof), recovered, st)
        call require(status_ok(st), 'GS weak Ampere solve')
        u(nc + 1:op%ndof) = recovered
        load = 0
        if (convention == 1) then
            call vrmhd_drive_build(op, VRMHD_CURRENT_SOURCE, 2.3_dp, 0.0_dp, &
                js, h, source, load, status)
            call require(status == 0, 'GS current-source drive')
        else if (convention == 2) then
            voltage = 0.2_dp*0.4_dp*2.3_dp
            if (kind == GS_ROTATION) voltage = -0.2_dp*0.2_dp*2*acos(-1.0_dp)
            call vrmhd_drive_build(op, VRMHD_LOOP_VOLTAGE, 2.3_dp, voltage, &
                js, h, source, load, status)
            call require(status == 0, 'GS harmonic voltage drive')
        end if
        load = load - surface
        call vrmhd_stationary_residual(op, steady, p, u, zero, u(:op%ndof), &
            source, load, r, unused, status)
        call require(status == 0, 'GS gauged stationary residual')
        pos = 0
        do j = 1, op%ndof
            if (op%fixed(j)) then
                pos = pos + 1; rows(pos) = j; cols(pos) = j; vals(pos) = 1
                cycle
            end if
            do entry = gram%col_ptr(j), gram%col_ptr(j + 1) - 1
                row = gram%row_idx(entry)
                if (op%fixed(row)) cycle
                pos = pos + 1; rows(pos) = row; cols(pos) = j
                vals(pos) = gram%val(entry)
            end do
        end do
        ng = 0; index = 0
        do i = 1, op%fs%n0
            if (op%fs%sp(0)%bnd(i)) cycle
            ng = ng + 1; index(i) = op%ndof + ng
        end do
        do t = 1, nt
            do q = 1, op%fs%nq
                call vrmhd_section_basis(op, t, q, x, weight, basis, status)
                call require(status == 0, 'GS gauge Riesz basis')
                do j = 1, op%fs%nl0
                    col = index(op%fs%sp(0)%dof(j, t))
                    if (col == 0) cycle
                    do i = 1, op%fs%nl0
                        row = index(op%fs%sp(0)%dof(i, t))
                        if (row == 0) cycle
                        pos = pos + 1; rows(pos) = row; cols(pos) = col
                        vals(pos) = weight*dot_product(basis(:, i), basis(:, j))
                    end do
                end do
            end do
        end do
        call csc_from_triplet(steady%ndof, steady%ndof, rows(:pos), &
            cols(:pos), vals(:pos), dual_gram, st)
        call require(status_ok(st), 'GS constrained dual Riesz matrix')
        call sparse_solve_once(dual_gram, r, dual, st)
        call require(status_ok(st), 'GS dual residual solve')
        defect = sqrt(max(0.0_dp, dot_product(r, dual)))
        parts = 0
        do i = 1, steady%ndof
            if (i <= op%fs%n0) then
                j = 1
            else if (i <= 4*op%fs%n0) then
                j = 2
            else if (i <= 5*op%fs%n0) then
                j = 5
            else if (i <= 6*op%fs%n0) then
                j = 6
            else if (i <= nc) then
                j = 3
            else if (i <= op%ndof) then
                j = 4
            else
                j = 7
            end if
            parts(j) = parts(j) + r(i)*dual(i)
        end do
        errors = 0
        do t = 1, nt
            do q = 1, op%fs%nq
                call vrmhd_section_basis(op, t, q, x, weight, basis, status)
                call require(status == 0, 'GS error basis')
                call oracle(kind, x, exact)
                field = matmul(basis, u(op%map(:, t))) - exact(:36)
                errors(1) = errors(1) + weight*(field(1)**2 + sum(field(13:15)**2))
                errors(2) = errors(2) + weight*(sum(field(5:7)**2) &
                    + sum(field(25:27)**2))
                errors(3) = errors(3) + weight*sum(field(8:10)**2)
            end do
        end do
        errors = sqrt(errors); error = norm2(errors(:2))
        print '(a,3i5,es14.6)', 'RM0 flux L2 error kind, convention, n ', &
            kind, convention, n, sqrt(gs_error)
        print '(a,3i5,3es14.6)', 'GS errors rho H1, A Hcurl, J L2 ', &
            kind, convention, n, errors
        print '(a,3i5,es14.6)', 'GS stationary dual defect ', kind, convention, n, defect
        print '(a,7es14.6)', 'GS dual parts rho,v,A,J,Te,Ti,lambda ', &
            sqrt(max(0.0_dp, parts))
    end subroutine
end program
