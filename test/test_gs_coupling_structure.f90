program test_gs_coupling_structure
    !! Structure of the two-loop mesh and of the coupled operator:
    !!  - the `boundary` flag marks exactly the nodes of Gamma_c; the Gamma_p
    !!    chain edges are curved (mid-edge nodes on the spline, off the chord)
    !!    and shared by two cells;
    !!  - the coupled matrix is symmetric, positive definite, annihilates
    !!    constants only through its stiffness part, and agrees with apply;
    !!  - energy: the exterior energy <S u, u> of the decaying extension of the
    !!    trace of a ring-current flux G(., y0), y0 inside Gamma_c, equals
    !!    -int_Gamma_c u d_nu u ds, with panel and mesh refinement; the interior
    !!    energy of the vacuum solution equals int_Gamma_c psi d_nu psi ds.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortfem_axisymmetric_gs_kernel, only: gs_green_gradient
    use kin6d_grad_shafranov, only: gs_state_t
    use kin6d_gs_coils, only: gs_circuit_set_t, gs_circuit_flux_vectors
    use kin6d_gs_coupling, only: gs_coupling_t, gs_coupling_build, gs_coupled_apply, &
                           gs_coupled_matrix, gs_coupling_coil_load, gs_coupled_solve, &
                           gs_coupling_energies, gs_coupling_edge_points, GS_COUPLING_OK
    use fortnum_quadrature, only: gauss_legendre
 use gs_coupling_fixture, only: gc_r0, gc_a, gc_b, gp_r0, gp_a, gp_b, pi, build_state, &
                                   build_coils
    implicit none

    interface
        subroutine dpotrf(uplo, n, a, lda, info)
            import :: dp
            character(len=1), intent(in) :: uplo
            integer, intent(in) :: n, lda
            real(dp), intent(inout) :: a(lda, *)
            integer, intent(out) :: info
        end subroutine dpotrf
    end interface

    real(dp), parameter :: ring(2) = [3.2_dp, 0.3_dp]
    real(dp), parameter :: currents(2) = [1.2e5_dp, -0.8e5_dp]
    integer :: failures

    failures = 0
    call test_mesh_flags(2)
    call test_mesh_flags(3)
    call test_operator(2)
    call test_operator(3)
    call test_exterior_energy()
    call test_interior_energy()
    if (failures /= 0) error stop 'gs_coupling structure oracle failures'
    print '(a)', 'gs_coupling structure and energy oracles passed'

contains

    subroutine require(ok, what)
        logical, intent(in) :: ok
        character(len=*), intent(in) :: what
        if (.not. ok) then
            failures = failures + 1
            print '(2a)', 'FAIL: ', what
        end if
    end subroutine require

    subroutine test_mesh_flags(order)
        integer, intent(in) :: order
        integer, parameter :: n_outer = 24
        type(gs_state_t) :: state
        real(dp) :: chain(2, n_outer / 2), chord(2), on_curve, off_chord, x(2)
        integer :: status, j, e, k, a, b, ja, jb, ids(2), curved_incidences, nc, i
        real(dp) :: rho2c

        call build_state(state, n_outer, order, status)
        call require(status == 0, 'two-loop mesh builds')
        nc = n_outer / 2
        do j = 1, nc
            chain(:, j) = [gp_r0 + gp_a * cos(2.0_dp * pi * (j - 1) / nc), &
                           gp_b * sin(2.0_dp * pi * (j - 1) / nc)]
        end do
        ! boundary flag: exactly the P2/P3 nodes of Gamma_c, all on the curve.
   call require(count(state%boundary) == n_outer*order, 'boundary flag only on Gamma_c')
        do i = 1, size(state % points, 2)
            if (.not. state % boundary(i)) cycle
   rho2c = ((state % points(1, i) - gc_r0) / gc_a)**2 + (state % points(2, i) / gc_b)**2
            if (abs(rho2c - 1.0_dp) > 1.0e-12_dp) then
                call require(.false., 'boundary node on Gamma_c')
                exit
            end if
        end do
        ! chain edges: both chain vertices exist, the edge is shared, and its
        ! interior nodes lie on the periodic spline, off the straight chord.
        curved_incidences = 0
        do e = 1, size(state % triangles, 2)
            do k = 1, 3
                a = state % triangles(k, e)
                b = state % triangles(mod(k, 3) + 1, e)
                ja = chain_index(state % points(:, a), chain)
                jb = chain_index(state % points(:, b), chain)
                if (ja == 0 .or. jb == 0) cycle
                if (mod(ja, nc) + 1 /= jb .and. mod(jb, nc) + 1 /= ja) cycle
                call require(state % neighbors(k, e) > 0, 'chain edge is interior')
                call require(state % curved_cell(e), 'cell on a chain edge is curved')
                ids(1) = state % triangles(3 + (k - 1) * (order - 1) + 1, e)
                x = state % points(:, ids(1))
                chord = (1.0_dp - real(1, dp) / order) * state % points(:, a) &
                        + real(1, dp) / order * state % points(:, b)
                off_chord = norm2(x - chord)
                on_curve = abs(((x(1) - gp_r0) / gp_a)**2 + (x(2) / gp_b)**2 - 1.0_dp)
                call require(off_chord > 1.0e-4_dp, 'chain edge node is off the chord')
                call require(on_curve < 5.0e-3_dp, 'chain edge node is on Gamma_p')
                curved_incidences = curved_incidences + 1
            end do
        end do
call require(curved_incidences == 2 * nc, 'each chain edge is a mesh edge of two cells')
    end subroutine test_mesh_flags

    pure function chain_index(x, chain) result(j)
        real(dp), intent(in) :: x(2), chain(:, :)
        integer :: j, i

        j = 0
        do i = 1, size(chain, 2)
            if (maxval(abs(x - chain(:, i))) < 1.0e-14_dp) j = i
        end do
    end function chain_index

    subroutine test_operator(order)
        integer, intent(in) :: order
        type(gs_state_t) :: state
        type(gs_coupling_t) :: coupling
        real(dp), allocatable :: a(:, :), x(:), y(:), z(:), chol(:, :), ones(:)
        real(dp) :: scale, interior, exterior
        integer :: status, info, i, n

        call build_state(state, 24, order, status)
        call gs_coupling_build(state, coupling, status, panels_per_edge=2)
        call require(status == GS_COUPLING_OK, 'coupling build')
        call gs_coupled_matrix(coupling, a)
        n = size(a, 1)
        scale = maxval(abs(a))
        call require(maxval(abs(a - transpose(a))) <= 1.0e-14_dp * scale, &
                     'coupled matrix is symmetric')
        call require(maxval(abs(coupling % s_trace - transpose(coupling % s_trace))) &
          <= 1.0e-14_dp * maxval(abs(coupling % s_trace)), 'Gamma_c block is symmetric')
        allocate (x(n), y(n), z(n), ones(n))
        do i = 1, n
            x(i) = sin(0.37_dp * i) + 0.1_dp * cos(1.3_dp * i)
        end do
        call gs_coupled_apply(coupling, x, y)
        z = matmul(a, x)
        call require(maxval(abs(y - z)) <= 1.0e-12_dp * maxval(abs(z)), &
                     'apply agrees with assembled matrix')
        ! Constants: the stiffness annihilates them, S does not.
        ones = 1.0_dp
        call gs_coupling_energies(coupling, ones, interior, exterior)
 call require(abs(interior) <= 1.0e-12_dp * exterior, 'stiffness annihilates constants')
        call require(exterior > 0.0_dp, 'exterior energy of a constant is positive')
        ! Positive definiteness through a Cholesky factorisation.
        chol = a
        call dpotrf('L', n, chol, n, info)
        call require(info == 0, 'coupled matrix is positive definite')
        print '(a,i1,a,i5,a,es10.2,a,es10.2)', 'P', order, ' nodes ', n, &
            ' asymmetry ', maxval(abs(a - transpose(a))) / scale, &
            ' constant exterior energy ', exterior
    end subroutine test_operator

    subroutine ring_flux(x, u, du_dr, du_dz)
        real(dp), intent(in) :: x(2)
        real(dp), intent(out) :: u, du_dr, du_dz
        real(dp) :: grad(4)

        call gs_green_gradient(x(1), x(2), ring(1), ring(2), u, grad)
        du_dr = grad(1)
        du_dz = grad(2)
    end subroutine ring_flux

    function exact_exterior_energy(state, coupling) result(energy)
        !! -int_Gamma_c u d_nu u ds for u = G(., ring), by Gauss quadrature.
        type(gs_state_t), intent(in) :: state
        type(gs_coupling_t), intent(in) :: coupling
        real(dp) :: energy
        real(dp) :: gx(12), gw(12), x(2), tangent(2), u, ur, uz, speed
        integer :: i, q

        call gauss_legendre(12, gx, gw)
        gx = 0.5_dp * (gx + 1.0_dp)
        gw = 0.5_dp * gw
        energy = 0.0_dp
        do i = 1, coupling % n_edges
            do q = 1, 12
                call gs_coupling_edge_points(state, coupling, i, gx(q), x, tangent)
                speed = norm2(tangent)
                call ring_flux(x, u, ur, uz)
             energy = energy - gw(q) * speed * u * (tangent(2) * ur - tangent(1) * uz) &
                         / (speed * x(1))
            end do
        end do
    end function exact_exterior_energy

    subroutine test_exterior_energy()
        integer, parameter :: meshes(3) = [24, 48, 96]
        integer, parameter :: panels(4) = [1, 2, 4, 8]
        type(gs_state_t) :: state
        type(gs_coupling_t) :: coupling
        real(dp) :: errors(4), mesh_errors(3), exact, u, ur, uz
        real(dp), allocatable :: trace(:)
        integer :: i, j, status

        ! Panel refinement on a fixed P3 trace.
        call build_state(state, 48, 3, status)
        do j = 1, size(panels)
            call gs_coupling_build(state, coupling, status, panels_per_edge=panels(j))
            call require(status == GS_COUPLING_OK, 'coupling build')
            if (j == 1) exact = exact_exterior_energy(state, coupling)
            allocate (trace(coupling % n_trace))
            do i = 1, coupling % n_trace
                call ring_flux(state % points(:, coupling % trace_node(i)), u, ur, uz)
                trace(i) = u
            end do
  errors(j) = abs(dot_product(trace, matmul(coupling % s_trace, trace)) - exact) / exact
            deallocate (trace)
        end do
     print '(a,4es11.3)', 'exterior energy error, P3 edges 48, panels 1 2 4 8: ', errors
        ! One panel per edge is BEM-limited; from two panels on the error is
        ! the (panel-independent) error of the P3 trace on this mesh.
        call require(errors(2) < 0.2_dp * errors(1), &
                     'exterior energy falls under panel refinement')
        call require(abs(errors(4) - errors(3)) < 0.05_dp * errors(3), &
                     'exterior energy panel error settles')
        ! Mesh refinement at four panels per edge.
        do j = 1, size(meshes)
            call build_state(state, meshes(j), 3, status)
            call gs_coupling_build(state, coupling, status, panels_per_edge=4)
            exact = exact_exterior_energy(state, coupling)
            allocate (trace(coupling % n_trace))
            do i = 1, coupling % n_trace
                call ring_flux(state % points(:, coupling % trace_node(i)), u, ur, uz)
                trace(i) = u
            end do
            mesh_errors(j) = abs(dot_product(trace, matmul(coupling % s_trace, trace)) &
                                 - exact) / exact
            deallocate (trace)
        end do
print '(a,3es11.3)', 'exterior energy error, P3 edges 24 48 96, panels 4: ', mesh_errors
        call require(all(mesh_errors(2:) < mesh_errors(:2)), &
                     'exterior energy falls under mesh refinement')
        call require(log(mesh_errors(2) / mesh_errors(3)) / log(2.0_dp) > 3.5_dp, &
                     'exterior energy order on the finest pair')
   call require(mesh_errors(3) < 5.0e-8_dp, 'exterior energy accurate on the fine mesh')
    end subroutine test_exterior_energy

    subroutine test_interior_energy()
        !! Vacuum solution: x^T K x = int_Gamma_c psi d_nu psi ds (L psi = 0).
        integer, parameter :: meshes(3) = [24, 48, 96]
        type(gs_state_t) :: state
        type(gs_circuit_set_t) :: set
        type(gs_coupling_t) :: coupling
        real(dp), allocatable :: load(:, :), psi(:)
        real(dp) :: gx(12), gw(12), x(2), tangent(2), speed, interior, exterior
        real(dp) :: exact, errors(3), flux(1, 2), dflux_r(1, 2), dflux_z(1, 2)
        integer :: j, i, q, status, iterations

        call build_coils(set)
        call gauss_legendre(12, gx, gw)
        gx = 0.5_dp * (gx + 1.0_dp)
        gw = 0.5_dp * gw
        do j = 1, size(meshes)
            call build_state(state, meshes(j), 3, status)
            call gs_coupling_build(state, coupling, status, panels_per_edge=2)
            call gs_coupling_coil_load(state, coupling, set, load, status)
            call gs_coupled_solve(coupling, matmul(load, currents), psi, status, &
                                  1.0e-11_dp, iterations)
            call require(status == GS_COUPLING_OK, 'coupled CG converged')
            call gs_coupling_energies(coupling, psi, interior, exterior)
            exact = 0.0_dp
            do i = 1, coupling % n_edges
                do q = 1, 12
                    call gs_coupling_edge_points(state, coupling, i, gx(q), x, tangent)
                    speed = norm2(tangent)
               call gs_circuit_flux_vectors(set, x(1:1), x(2:2), flux, dflux_r, dflux_z)
               exact = exact + gw(q) * dot_product(flux(1, :), currents) * (tangent(2) &
                                   * dot_product(dflux_r(1, :), currents) - tangent(1) &
                                          * dot_product(dflux_z(1, :), currents)) / x(1)
                end do
            end do
            errors(j) = abs(interior - exact) / exact
        end do
     print '(a,3es11.3)', 'interior energy error, P3 edges 24 48 96, panels 2: ', errors
    call require(all(errors(2:) < errors(:2)), 'interior energy falls under refinement')
        call require(log(errors(2) / errors(3)) / log(2.0_dp) > 3.5_dp, &
                     'interior energy order on the finest pair')
        call require(errors(3) < 5.0e-8_dp, 'interior energy accurate on the fine mesh')
    end subroutine test_interior_energy

end program test_gs_coupling_structure
