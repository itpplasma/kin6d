module gll_benchmark_context
    ! A single benchmark invocation owns this context. Module callbacks avoid
    ! executable-stack trampolines on GNU 12; this is not a solver library.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_ideal_gll_spike
    use fortfem_tensor_solver, only: tensor_multigrid_t, build_tensor_multigrid, &
        tensor_vcycle, free_tensor_multigrid
    use fortfem_tensor_mesh, only: tensor_mesh_t
    use fortnum_krylov, only: real_conjugate_gradient_operator, KRYLOV_OK
    use omp_lib, only: omp_get_wtime, omp_get_max_threads
    implicit none
    type(gll_map_t) :: map
    type(gll_block_jacobi_t) :: node_pc, scalar_pc
    type(tensor_multigrid_t) :: mg, omega_mg
    real(dp) :: shift = 0, hv_seconds = 0, pc_seconds = 0
    integer :: hv_calls = 0, pc_calls = 0, pc_status = 0
    integer :: pc_mode = 1
contains
    subroutine apply(v, y)
        real(dp), intent(in) :: v(:)
        real(dp), intent(out) :: y(:)
        real(dp) :: started
        started = omp_get_wtime()
        call apply_gll_hessian(map, v, y)
        where (.not. map%fixed) y = y + shift*v
        hv_seconds = hv_seconds + omp_get_wtime() - started
        hv_calls = hv_calls + 1
    end subroutine

    subroutine precondition(v, y)
        real(dp), intent(in) :: v(:)
        real(dp), intent(out) :: y(:)
        real(dp) :: started
        integer :: c, lo, hi, st
        started = omp_get_wtime()
        select case (pc_mode)
        case (1)
            do c = 1, 3
                lo = 1 + (c - 1)*map%nodes; hi = c*map%nodes
                if (c < 3) then
                    call tensor_vcycle(mg, 1, v(lo:hi), y(lo:hi), st)
                else
                    call tensor_vcycle(omega_mg, 1, v(lo:hi), y(lo:hi), st)
                end if
                if (st /= 0) pc_status = st
            end do
        case (3)
            call apply_gll_block_jacobi(scalar_pc, v, y, st)
            if (st /= 0) pc_status = st
        case (4)
            call apply_gll_block_jacobi(node_pc, v, y, st)
            if (st /= 0) pc_status = st
        case default
            y = v
            pc_status = 1
        end select
        where (map%fixed) y = v
        pc_seconds = pc_seconds + omp_get_wtime() - started
        pc_calls = pc_calls + 1
    end subroutine
end module

program benchmark_gll
    use gll_benchmark_context
    use kin6d_landreman_benchmark, only: landreman_frozen_t
    use kin6d_landreman_leaves, only: landreman_surface
    use kin6d_nested_oracles, only: nested_oracle_t
    use kin6d_grad_shafranov, only: gs_state_t, initialize_gs
    use kin6d_equilibrium_chart, only: cylindrical_chart_t
    use kin6d_ideal_map, only: ideal_map_t, ideal_map_init, ideal_map_freeze, &
        ideal_map_energy
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    real(dp), parameter :: delta = 1.0_dp/64, pi = acos(-1.0_dp)
    type(landreman_frozen_t) :: field
    type(nested_oracle_t) :: oracle
    type(tensor_mesh_t) :: omega_mesh
    real(dp), allocatable :: x(:), base(:), g(:), gp(:), gm(:), v(:), hv(:), &
        w(:), hw(:), solution(:), rhs(:), check(:)
    real(dp) :: energy, exact, plus, minus, plus2, minus2, h, fdg, fdg4, fdh, &
        errg, errh, err2, &
        fd2, started, setup, mgsetup, blocksetup, diagsetup, pcsetup, &
        linearization, evaluation, hvtime, residual, &
        normg, cgtime, eafter, symmetry, pairing, s, r, z, moments(4)
    integer :: p, n, status, i, j, iterations, nd, rep, repetitions, pc
    character(32) :: argument
    p = 4; n = 1; repetitions = 20; shift = 1.e-2_dp
    call get_command_argument(1, argument)
    if (len_trim(argument) > 0) read(argument, *) p
    call get_command_argument(2, argument)
    if (len_trim(argument) > 0) read(argument, *) n
    call get_command_argument(3, argument)
    if (len_trim(argument) > 0) read(argument, *) repetitions
    call get_command_argument(4, argument)
    if (len_trim(argument) > 0) read(argument, *) shift
    if (shift <= 0) error stop 'positive damping required'
    if (p < 2 .or. p > 8 .or. n < 1 .or. repetitions < 1) error stop 'arguments'
    field%delta = delta
    oracle%eps = 0; oracle%tau = 0; oracle%delta = delta
    call oracle%averages(32, moments)
    exact = moments(1)*(moments(3)/2 - moments(2))
    started = omp_get_wtime()
    call build_gll_map(p, n, delta, field, map, x, status)
    if (status /= 0) error stop 'GLL setup'
    setup = omp_get_wtime() - started
    nd = size(x)
    allocate(g(nd), gp(nd), gm(nd), v(nd), hv(nd), w(nd), hw(nd), &
        solution(nd), rhs(nd), check(nd))
    base = x
    call evaluate_gll_map(map, x, energy, g, status, .false.)
    if (status /= 0) error stop 'identity energy'
    if (abs(energy - exact)/abs(exact) > 1.e-5_dp) error stop 'GLL identity energy'
    print '(a)', 'case=L_axisymmetric_eps0_tau0_delta1over64; certification=unavailable'
    print '(a,4(i0,1x))', 'p,n,dofs,threads=', p, n, nd, omp_get_max_threads()
    print '(a,3(es24.16,1x))', 'identity_energy,exact,rel_error=', &
        energy, exact, abs(energy - exact)/abs(exact)
    if (p == 4 .and. n == 1) call triangle_comparison(exact)
    ! Smooth admissible perturbation away from the identity equilibrium.
    do i = 1, map%nodes
        r = base(i); z = base(map%nodes + i)
        s = ((r*r - 1)**2/4 + z*z)/delta
        x(i) = x(i) + 1.e-4_dp*(1 - s)*cos(17*z)
        x(map%nodes + i) = x(map%nodes + i) + 1.e-4_dp*(1 - s)*sin(11*r)
        x(2*map%nodes + i) = 1.e-3_dp*(1 - s)*sin(13*r + 9*z)
    end do
    where (map%fixed) x = base
    do i = 1, nd
        v(i) = sin(real(17*i, dp)); w(i) = cos(real(13*i, dp))
    end do
    where (map%fixed)
        v = 0; w = 0
    end where
    v = v/norm2(v); w = w/norm2(w)
    started = omp_get_wtime()
    call evaluate_gll_map(map, x, energy, g, status, .true.)
    if (status /= 0) error stop 'linearization'
    linearization = omp_get_wtime() - started
    call apply_gll_hessian(map, v, hv)
    call apply_gll_hessian(map, w, hw)
    pairing = dot_product(v, hv)
    symmetry = abs(dot_product(w, hv) - dot_product(v, hw)) &
        /max(1.0_dp, norm2(hv), norm2(hw))
    errg = huge(h); errh = huge(h); err2 = huge(h)
    do j = 2, 6
        h = 10.0_dp**(-j)
        call evaluate_gll_map(map, x + h*v, plus, gp, status, .false.)
        if (status /= 0) cycle
        call evaluate_gll_map(map, x - h*v, minus, gm, status, .false.)
        if (status /= 0) cycle
        fdg = (plus - minus)/(2*h)
        call evaluate_gll_map(map, x + 2*h*v, plus2, check, status, .false.)
        if (status /= 0) cycle
        call evaluate_gll_map(map, x - 2*h*v, minus2, check, status, .false.)
        if (status /= 0) cycle
        fdg4 = (-plus2 + 8*plus - 8*minus + minus2)/(12*h)
        fdh = norm2((gp - gm)/(2*h) - hv)/max(1.0_dp, norm2(hv))
        fd2 = abs((plus - 2*energy + minus)/h**2 - pairing)/max(1.0_dp, abs(pairing))
        errg = min(errg, abs(fdg4 - dot_product(g, v)) &
            /max(1.e-8_dp, abs(dot_product(g, v))))
        errh = min(errh, fdh); err2 = min(err2, fd2)
        print '(a,5(es14.6,1x))', 'fd h,g2_rel,g4_rel,hvp_rel,energy_hessian_rel=', &
            h, abs(fdg - dot_product(g, v))/max(1.e-8_dp, abs(dot_product(g, v))), &
            abs(fdg4 - dot_product(g, v))/max(1.e-8_dp, abs(dot_product(g, v))), &
            fdh, fd2
    end do
    print '(a,es14.6)', 'hessian_pairing_symmetry=', symmetry
    if (errg > 2.e-5_dp .or. errh > 2.e-5_dp .or. err2 > 2.e-4_dp &
        .or. symmetry > 1.e-12_dp) error stop 'derivative oracle'
    ! A nonzero perturbation on fixed columns must not leak into free rows.
    check = 0
    where (map%fixed) check = 1
    call apply_gll_hessian(map, check, solution)
    where (map%fixed) solution = solution - check
    if (maxval(abs(solution)) > 1.e-12_dp) error stop 'fixed row/column mask'
    started = omp_get_wtime()
    do rep = 1, repetitions
        call apply_gll_hessian(map, v, hv)
    end do
    hvtime = (omp_get_wtime() - started)/repetitions
    started = omp_get_wtime()
    do rep = 1, repetitions
        call evaluate_gll_map(map, x, plus, gp, status, .false.)
    end do
    evaluation = (omp_get_wtime() - started)/repetitions
    started = omp_get_wtime()
    call build_tensor_multigrid(map%reference%mesh, diffusion, mg, status)
    if (status /= 0) error stop 'pMG RZ setup'
    omega_mesh = map%reference%mesh
    omega_mesh%boundary = .false.
    call build_tensor_multigrid(omega_mesh, regularized_diffusion, omega_mg, status)
    if (status /= 0) error stop 'pMG omega setup'
    mgsetup = omp_get_wtime() - started
    started = omp_get_wtime()
    call build_gll_block_jacobi(map, shift, scalar_pc, status, diagonal_only=.true.)
    if (status /= 0) error stop 'unchanged-shift diagonal is not positive'
    diagsetup = omp_get_wtime() - started
    started = omp_get_wtime()
    call build_gll_block_jacobi(map, shift, node_pc, status)
    if (status /= 0) error stop 'unchanged-shift nodal block is not SPD'
    blocksetup = omp_get_wtime() - started
    normg = norm2(g); rhs = -g/normg
    print '(a,5(es14.6,1x))', 'setup,mg_setup,linearize,energy_gradient,hvp_s=', &
        setup, mgsetup, linearization, evaluation, hvtime
    print '(a,es14.6)', 'coupled_blocks_setup_s=', blocksetup
    print '(a,es14.6)', 'scalar_Jacobi_setup_s=', diagsetup
    print '(a)', 'pc1=scalar_pMG;pc2=identity;pc3=scalar_Jacobi;pc4=coupled_node3'
    ! One damped Newton system at a common user-selected shift across p.
    do pc = 1, 4
        pc_mode = pc
        solution = 0; hv_seconds = 0; pc_seconds = 0; hv_calls = 0; pc_calls = 0
        pc_status = 0
        started = omp_get_wtime()
        if (pc == 2) then
            call real_conjugate_gradient_operator(apply, rhs, solution, 1.e-8_dp, &
                5000, status, iterations, residual)
        else
            call real_conjugate_gradient_operator(apply, rhs, solution, 1.e-8_dp, &
                5000, status, iterations, residual, precondition)
        end if
        cgtime = omp_get_wtime() - started
        if (status /= KRYLOV_OK .or. pc_status /= 0) error stop 'damped Newton CG'
        call apply(solution, check)
        residual = norm2(check - rhs)
        call evaluate_gll_map(map, x + normg*solution, eafter, gp, i, .false.)
        print '(a,4(i0,1x),7(es14.6,1x))', &
            'cg pc,status,it,hv_calls;shift,res,total,hv,pc,dE,gdotstep=', &
            pc, status, iterations, hv_calls, shift, residual, cgtime, hv_seconds, &
            pc_seconds, eafter - energy, normg*dot_product(g, solution)
        if (residual > 2.e-8_dp .or. .not. ieee_is_finite(residual) &
            .or. .not. ieee_is_finite(eafter)) &
            error stop 'true Newton residual'
        if (i /= 0 .or. eafter >= energy .or. dot_product(g, solution) >= 0) &
            error stop 'Newton step must be unfolded and decrease energy'
        pcsetup = 0
        if (pc == 1) pcsetup = mgsetup
        if (pc == 3) pcsetup = diagsetup
        if (pc == 4) pcsetup = blocksetup
        print '(a,i0,a,4(es14.6,1x))', 'cg_cost pc=', pc, &
            ';pc_setup,solve,pc_setup_plus_solve,common_setup_plus_solve_s=', &
            pcsetup, cgtime, pcsetup + cgtime, &
            setup + linearization + pcsetup + cgtime
    end do
    call free_tensor_multigrid(mg)
    call free_tensor_multigrid(omega_mg)
contains
    subroutine wall(t, point)
        real(dp), intent(in) :: t
        real(dp), intent(out) :: point(2)
        real(dp) :: dr, dz
        call landreman_surface(2*pi*t, delta, point(1), point(2), dr, dz)
    end subroutine

    subroutine triangle_comparison(exact_energy)
        real(dp), intent(in) :: exact_energy
        type(gs_state_t) :: mesh
        type(ideal_map_t) :: triangle
        type(cylindrical_chart_t) :: chart
        real(dp) :: contour(2, 32), et, max_area
        integer :: bnodes, a, degree, st
        do degree = 2, 6, 2
            do bnodes = 8, 32, 8
                do a = 1, bnodes
                    call wall(real(a - 1, dp)/bnodes, contour(:, a))
                end do
                max_area = (2*pi*sqrt(delta)/bnodes)**2/2
                call initialize_gs(mesh, contour(:, :bnodes), max_area, st, &
                    degree, curve=wall, geometry_only=.true.)
                if (st /= 0) error stop 'triangle mesh'
                call ideal_map_init(triangle, mesh, chart, 1, 0, st, &
                    quadrature_order=2*degree + 7)
                if (st /= 0) error stop 'triangle map'
                triangle%gauge_weight = 0
                call ideal_map_freeze(triangle, field, st)
                if (st /= 0) error stop 'triangle freeze'
                call ideal_map_energy(triangle, triangle%x, et, st)
                if (st /= 0) error stop 'triangle energy'
                print '(a,3(i0,1x),2(es24.16,1x))', &
                    'triangle p,bnodes,nodes,E,rel_error=', &
                    degree, bnodes, triangle%nodes, et, &
                    abs(et - exact_energy)/abs(exact_energy)
                if (bnodes == 32 .and. degree == 6) then
                    if (abs(et - exact_energy)/abs(exact_energy) > 1.e-8_dp) &
                        error stop 'triangle refinement energy'
                end if
            end do
        end do
    end subroutine
end program
