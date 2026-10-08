program test_gs_p2
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_grad_shafranov, only: gs_state_t, initialize_gs, solve_gs, &
        sample_gs, gs_element_jet, gs_axis_candidate, gs_axis_safety_factor, &
        gs_classical_current, gs_current_jump, gs_linear_response, write_gs_mesh
    use kin6d_grad_shafranov, only: solve_gs_profiles
    use kin6d_grad_shafranov_leaves, only: gs_smooth_axis_q_jvp, &
        gs_smooth_axis_q_vjp, gs_current_density_jvp, gs_current_density_vjp
    implicit none
    type(gs_state_t) :: state, vertex_case, p1
    real(dp) :: boundary(2, 4), triangle(2, 3), axis(2), wanted_axis(2), h(2, 2)
    real(dp) :: point(2), grad(2), value, br, bz, bt, current, q, expected, hc(2, 2)
    real(dp) :: left_center(2), right_center(2), edge_origin(2), normal(2), height
    real(dp) :: jump, reverse, lhs, rhs, dq, ar, af, ahrr, ahrz, ahzz
    real(dp) :: aj, amu, apr, dhrr, dhzz
    real(dp), allocatable :: source(:), goal(:), tangent(:), adjoint(:)
    integer :: status, i, e, j, neighbor, edge, opposite, axis_element, passed
    passed = 0
    boundary = reshape([5.8_dp, -0.4_dp, 6.6_dp, -0.4_dp, &
                         6.6_dp, 0.4_dp, 5.8_dp, 0.4_dp], [2, 4])
    call initialize_gs(state, boundary, 0.012_dp, status, 2)
    call require(status == 0, 'native P2 initialization')
    call require(size(state%triangles, 1) == 6, 'six local scalar DOFs')
    do e = 1, size(state%triangles, 2)
        do j = 1, 3
            point = 0.5_dp*(state%points(:, state%triangles(j, e)) + &
                state%points(:, state%triangles(mod(j, 3)+1, e)))
            call close(norm2(point-state%points(:, state%triangles(j+3, e))), &
                       0.0_dp, 1.0e-13_dp, 'physical midpoint identity')
            neighbor = state%neighbors(j, e)
            call require(state%boundary(state%triangles(j+3, e)) .eqv. &
                (neighbor < 1), 'boundary midpoint topology')
            if (neighbor < 1) cycle
            call require(any(state%triangles(4:6, neighbor) == &
                state%triangles(j+3, e)), 'global edge DOF shared with neighbor')
        end do
    end do
    ! Independent full quadratic, deliberately centered strictly inside an element.
    wanted_axis = sum(state%points(:, state%triangles(1:3, 1)), dim=2)/3
    hc = reshape([2.0_dp, 0.5_dp, 0.5_dp, 3.0_dp], [2, 2])
    do i = 1, size(state%psi)
        point = state%points(:, i)-wanted_axis
        state%psi(i) = 1.0_dp + point(1)**2 + 0.5_dp*point(1)*point(2) + &
            1.5_dp*point(2)**2
    end do
    do e = 1, size(state%triangles, 2)
        point = sum(state%points(:, state%triangles(1:3, e)), dim=2)/3
        call gs_element_jet(state, e, point(1), point(2), value, grad, h, status)
        call require(status == 0, 'P2 jet exists at each physical centroid')
        call close(value, 1+(point(1)-wanted_axis(1))**2 + &
            0.5_dp*(point(1)-wanted_axis(1))*(point(2)-wanted_axis(2)) + &
            1.5_dp*(point(2)-wanted_axis(2))**2, 3.0e-12_dp, 'quadratic value')
        call close(maxval(abs(grad-matmul(hc, point-wanted_axis))), &
            0.0_dp, 3.0e-11_dp, 'independent quadratic physical gradient')
        call close(maxval(abs(h-hc)), 0.0_dp, 5.0e-10_dp, 'quadratic physical Hessian')
        call gs_classical_current(state, e, point(1), point(2), 2.0_dp, current, status)
        expected = -(5.0_dp-(2*(point(1)-wanted_axis(1)) + &
            0.5_dp*(point(2)-wanted_axis(2)))/point(1))/(2*point(1))
        call close(current, expected, 4.0e-11_dp, 'classical curl current')
    end do
    call gs_axis_candidate(state, axis, value, h, axis_element, status)
    call require(status == 0, 'regular strict interior stationary axis')
    call close(norm2(axis-wanted_axis), 0.0_dp, 3.0e-11_dp, 'independent exact axis')
    call close(value, 1.0_dp, 2.0e-12_dp, 'axis value uses P2 flux')
    call gs_axis_safety_factor(axis(1), -30.0_dp, h, q, status)
    call require(status == 0, 'signed elliptic axis-q candidate')
    call close(q, independent_coarea(axis(1), -30.0_dp), &
        3.0e-9_dp, 'independent elliptic contour coarea axis-q limit')
    ! All interior traces of the global quadratic have zero sheet current.
    edge = 0; e = 0
    outer: do i = 1, size(state%triangles, 2)
        do j = 1, 3
            if (state%neighbors(j, i) < 1) cycle
            e = i; edge = j
            exit outer
        end do
    end do outer
    call require(edge /= 0, 'interior face available')
    neighbor = state%neighbors(edge, e)
    edge_origin = state%points(:, state%triangles(edge, e))
    point = 0.5_dp*(edge_origin + &
        state%points(:, state%triangles(mod(edge, 3)+1, e)))
    normal = state%points(:, state%triangles(mod(edge, 3)+1, e))-edge_origin
    normal = [normal(2), -normal(1)]/norm2(normal)
    left_center = sum(state%points(:, state%triangles(1:3, e)), dim=2)/3
    right_center = sum(state%points(:, state%triangles(1:3, neighbor)), dim=2)/3
    if (dot_product(normal, right_center-left_center) < 0) normal = -normal
    call gs_current_jump(state, e, neighbor, point(1), point(2), &
                         normal, 2.0_dp, jump, status)
    call require(status == 0, 'shared-face current diagnostic')
    call close(jump, 0.0_dp, 1.0e-11_dp, 'global quadratic has no sheet current')
    opposite = mod(edge+1, 3)+1
    height = dot_product(normal, edge_origin-state%points(:, state%triangles(opposite, e)))
    i = state%triangles(opposite, e)
    state%psi(i) = state%psi(i)+1
    call gs_current_jump(state, e, neighbor, point(1), point(2), &
                         normal, 2.0_dp, jump, status)
    call require(status == 0, 'perturbed continuous P2 face diagnostic')
    call close(jump, 1/(height*2*point(1)), 3.0e-10_dp, &
        'independent opposite-vertex polynomial sheet current')
    call gs_current_jump(state, neighbor, e, point(1), point(2), &
                         -normal, 2.0_dp, reverse, status)
    call close(jump, reverse, 3.0e-11_dp, 'physical sheet independent of face orientation')
    call gs_current_jump(state, e, neighbor, point(1), point(2), &
                         -normal, 2.0_dp, reverse, status)
    call require(status /= 0, 'wrong normal orientation is rejected')
    ! Exact manufactured GS solution: psi=R^2+Z has Delta*psi=0.
    ! Initialize only boundary data; no supplied interior answer.
    state%psi = 0
    where (state%boundary) state%psi = state%points(1, :)**2 + state%points(2, :)
    call solve_gs(state, 2.0_dp, 0.0_dp, 0.0_dp, status)
    call require(status == 0, 'P2 manufactured homogeneous weak solve')
    call close(maxval(abs(state%psi-state%points(1, :)**2-state%points(2, :))), &
               0.0_dp, 2.0e-8_dp, 'independent R^2+Z solution from boundary only')
    call sample_gs(state, 6.21_dp, 0.053_dp, -30.0_dp, value, br, bz, bt, status)
    call require(status == 0, 'consistent P2 physical sample')
    call close(br, -1/6.21_dp, 2.0e-8_dp, 'independent signed radial field')
    call close(bz, 2.0_dp, 2.0e-8_dp, 'independent vertical field')
    call close(bt, -30/6.21_dp, 1.0e-12_dp, 'independent negative toroidal field')
    ! psi=Z^2 has Delta*psi=2, hence FFprime=-2 and j_phi=-2/(mu0*R).
    state%psi = 0
    where (state%boundary) state%psi = state%points(2, :)**2
    call solve_gs(state, 2.0_dp, 0.0_dp, -2.0_dp, status)
    call require(status == 0, 'P2 manufactured nonzero-source weak solve')
    call close(maxval(abs(state%psi-state%points(2, :)**2)), &
               0.0_dp, 2.0e-8_dp, 'independent Z^2 weak current source')
    do e = 1, size(state%triangles, 2)
        point = sum(state%points(:, state%triangles(1:3, e)), dim=2)/3
        call gs_classical_current(state, e, point(1), point(2), 2.0_dp, current, status)
        call close(current, -1/point(1), 2.0e-7_dp, 'manufactured signed source current')
    end do
    state%psi = 0
    where (state%boundary) state%psi = state%points(1, :)**2
    call solve_gs_profiles(state, 1.0_dp, [0.0_dp, 50.0_dp], &
        [-0.01_dp, -0.01_dp], [0.0_dp, 0.5_dp], status, 1.0e-11_dp)
    call require(status == 0, 'P2 flux-dependent weak equation')
    call close(maxval(abs(state%psi-state%points(1, :)**2)), &
        0.0_dp, 2.0e-7_dp, 'independent nonlinear R^2 source cancellation')
    allocate (source(size(state%psi)), goal(size(state%psi)))
    do i = 1, size(source)
        source(i) = sin(real(i, dp)); goal(i) = cos(real(i, dp))
    end do
    where (state%boundary) source = 0
    where (state%boundary) goal = 0
    call gs_linear_response(state, source, tangent, status)
    call require(status == 0, 'P2 equation tangent')
    call gs_linear_response(state, goal, adjoint, status)
    call require(status == 0, 'P2 equation adjoint')
    lhs = dot_product(goal, tangent); rhs = dot_product(source, adjoint)
    call close(lhs, rhs, 3.0e-8_dp*max(1.0_dp, abs(lhs)), 'P2 implicit pairing')
    call gs_smooth_axis_q_jvp(6.0_dp, -30.0_dp, 2.0_dp, 0.5_dp, 3.0_dp, &
        0.3_dp, -0.2_dp, 0.1_dp, -0.4_dp, 0.2_dp, dq)
    call gs_smooth_axis_q_vjp(6.0_dp, -30.0_dp, 2.0_dp, 0.5_dp, 3.0_dp, &
        1.3_dp, ar, af, ahrr, ahrz, ahzz)
    call close(1.3_dp*dq, 0.3_dp*ar-0.2_dp*af+0.1_dp*ahrr-0.4_dp*ahrz+ &
        0.2_dp*ahzz, 1.0e-12_dp, 'generated axis-q JVP/VJP pairing')
    call gs_current_density_jvp(6.0_dp, 2.0_dp, 1.0_dp, 2.0_dp, 3.0_dp, &
        0.3_dp, -0.2_dp, 0.1_dp, -0.4_dp, 0.2_dp, aj)
    call gs_current_density_vjp(6.0_dp, 2.0_dp, 1.0_dp, 2.0_dp, 3.0_dp, &
        1.3_dp, ar, amu, apr, dhrr, dhzz)
    call close(1.3_dp*aj, 0.3_dp*ar-0.2_dp*amu+0.1_dp*apr-0.4_dp*dhrr+ &
        0.2_dp*dhzz, 1.0e-12_dp, 'generated classical-current JVP/VJP pairing')
    ! Positive Hessian alone must not turn a vertex/edge minimum into a smooth axis.
    triangle = reshape([6.0_dp, 0.0_dp, 7.0_dp, 0.0_dp, &
                          6.0_dp, 1.0_dp], [2, 3])
    call initialize_gs(vertex_case, triangle, 10.0_dp, status, 2)
    call require(status == 0, 'single-element rejection control')
    call require(size(vertex_case%triangles, 2) == 1, 'single-element geometry')
    vertex_case%psi = (vertex_case%points(1, :)-6)**2 + &
        vertex_case%points(2, :)**2 + vertex_case%points(1, :)-6 + &
        vertex_case%points(2, :)
    call gs_axis_candidate(vertex_case, axis, value, h, axis_element, status)
    call require(status == 4, 'vertex cusp is not assigned a finite smooth axis-q')
    vertex_case%psi = (vertex_case%points(1, :)-6.5_dp)**2 + &
        vertex_case%points(2, :)**2
    call gs_axis_candidate(vertex_case, axis, value, h, axis_element, status)
    call require(status == 4, 'edge stationary minimum is explicitly unqualified')
    call initialize_gs(p1, triangle, 10.0_dp, status)
    call gs_axis_candidate(p1, axis, value, h, axis_element, status)
    call require(status == 3, 'P1 Hessian is not fabricated')
    call write_gs_mesh(p1, 'p1_native_oracle_mesh.dat', status)
    call require(status == 0, 'native P1 mesh-v1 export remains supported')
    call initialize_gs(vertex_case, triangle, 10.0_dp, status, 3)
    call require(status /= 0, 'unsupported polynomial order rejected')
    call write_gs_mesh(state, 'p2_native_oracle_mesh.dat', status)
    call require(status == 0, 'native six-DOF mesh-v2 export')
    print *, 'PASS: native P2 independent physical oracles', passed
contains
    subroutine require(condition, name)
        logical, intent(in) :: condition
        character(*), intent(in) :: name
        if (.not. condition) then
            print *, 'FAIL: ', name
            error stop 1
        end if
        passed = passed+1
    end subroutine require
    subroutine close(actual, wanted, tolerance, name)
        real(dp), intent(in) :: actual, wanted, tolerance
        character(*), intent(in) :: name
        if (abs(actual-wanted) > tolerance) then
            print *, 'FAIL: ', name, actual, wanted, tolerance
            error stop 1
        end if
        passed = passed+1
    end subroutine close
    function independent_coarea(r_axis, f_axis) result(q_axis)
        real(dp), intent(in) :: r_axis, f_axis
        real(dp) :: q_axis, theta, x, z, dx, dz, gr, gz, rho
        integer :: k
        integer, parameter :: n = 20000
        ! Explicit rotated ellipse: 2*x^2+x*z+3*z^2=rho^2.
        ! Direct dl/(R*|grad psi|), no production Hessian determinant formula.
        rho = 1.0e-5_dp
        q_axis = 0
        do k = 0, n-1
            theta = 2*acos(-1.0_dp)*k/n
            z = rho*sin(theta)/sqrt(2.875_dp)
            x = (rho*cos(theta)-0.5_dp*z/sqrt(2.0_dp))/sqrt(2.0_dp)
            dz = rho*cos(theta)/sqrt(2.875_dp)
            dx = (-rho*sin(theta)-0.5_dp*dz/sqrt(2.0_dp))/sqrt(2.0_dp)
            gr = 2*x+0.5_dp*z; gz = 0.5_dp*x+3*z
            q_axis = q_axis + norm2([dx, dz])/((r_axis+x)*norm2([gr, gz]))
        end do
        q_axis = f_axis*q_axis/n
    end function independent_coarea
end program test_gs_p2
