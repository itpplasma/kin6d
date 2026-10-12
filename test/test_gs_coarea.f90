program test_gs_coarea
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_grad_shafranov, only: gs_state_t, &
                                    initialize_gs, gs_axis_candidate, &
                                    gs_element_jet
    use kin6d_gs_coarea, only: gs_level_coarea
    use kin6d_gs_coarea_leaves, &
        only: gs_contour_root_response
    implicit none
    type(gs_state_t) :: state, plus_state, minus_state
    real(dp) :: contour(2, 4), axis(2), wanted_axis(2), &
                centroid(2), distance
    real(dp) :: levels(17), coarse(17), fine(17), &
                radius2, expected, f, q
    real(dp) :: response, closest, two_pi, span, &
                tangent_coarea(17), tangent_levels(17)
    real(dp), allocatable::tangent_psi(:)
    real(dp)::hessian(2, 2), gradient(2), axis_flux, &
               perturbation_step, difference_error
    real(dp)::plus_levels(17), minus_levels(17), &
               plus_coarea(17), minus_coarea(17), tvalue
    integer :: status, e, i, element
    contour = reshape([5.8_dp, -0.4_dp, 6.6_dp, -0.4_dp, &
                       6.6_dp, 0.4_dp, 5.8_dp, 0.4_dp], [2, 4])
    call initialize_gs(state, contour, 0.008_dp, status, 2)
    if (status /= 0) error stop 'coarea polynomial mesh failed'
    closest = huge(1.0_dp)
    do e = 1, size(state%triangles, 2)
        centroid = sum(state%points(:, state%triangles(1:3, &
                                                       e)), dim=2)/3
        distance = norm2(centroid - [6.2_dp, 0.0_dp])
        if (distance >= closest) cycle
        closest = distance; wanted_axis = centroid
    end do
    span = 0.0625_dp
    state%psi = (state%points(1, :) - wanted_axis(1))**2 + &
(state%points(2, :) - wanted_axis(2))**2 - span
    do i = 1, size(levels)
        levels(i) = -span + span*(i - 1)/(size(levels) - 1)
    end do
    call gs_level_coarea(state, levels, 128, coarse, &
                         axis, status)
    if (status /= 0) then
        print *, 'coarea status:', status
        error stop 'native P2 coarea rejected an exact regular quadratic'
    end if
    if (norm2(axis - &
              wanted_axis) > 1.0e-11_dp) error stop 'coarea changed the native axis'
    call gs_level_coarea(state, levels, 256, fine, axis, &
                         status)
    if (status /= 0) error stop 'coarea angular refinement failed'
    two_pi = 2*acos(-1.0_dp)
    do i = 1, size(levels)
        radius2 = levels(i) + span
        ! psi=radius^2-span => |grad psi|=2radius. Exact finite-R contour integral.
        expected = two_pi/(2*sqrt(wanted_axis(1)**2 - radius2))
        if (abs(fine(i)/expected - 1) > 2.0e-9_dp) &
            error stop 'independent finite-R circular coarea oracle'
        q = -30.0_dp/(2*sqrt(wanted_axis(1)**2 - radius2))
        f = q*two_pi/fine(i)
        if (abs(f + &
                30) > 7.0e-8_dp) error stop 'signed F/q coarea correspondence'
    end do
    if (maxval(abs(fine/coarse - 1)) > 2.0e-9_dp) &
        error stop 'independent periodic angular quadrature control'
    tangent_psi = state%psi; tangent_levels = levels
    call gs_level_coarea(state, levels, 256, fine, axis, &
                         status, tangent_psi, &
                         tangent_levels, tangent_coarea)
    if (status /= 0) then
        print *, 'scale tangent status:', status
        error stop 'native coarea scale tangent rejected'
    end if
    if (maxval(abs(tangent_coarea/fine + 1)) > 2.0e-9_dp) &
        error stop 'independent inverse-amplitude coarea derivative'
    tangent_psi = -2*(state%points(1, :) - &
                      wanted_axis(1)); tangent_levels = 0
    call gs_level_coarea(state, levels, 256, fine, axis, &
                         status, tangent_psi, &
                         tangent_levels, tangent_coarea)
    if (status /= 0) error stop 'native coarea moving-axis tangent rejected'
    do i = 1, size(levels)
        radius2 = levels(i) + span
        expected = -fine(i)*wanted_axis(1)/(wanted_axis(1)**2 - radius2)
        if (abs(tangent_coarea(i) - expected) > 2.0e-9_dp) &
            error stop 'independent finite-R translated-axis coarea derivative'
    end do
    tangent_psi = 1; tangent_levels = 1
    call gs_level_coarea(state, levels, 256, fine, axis, &
                         status, tangent_psi, &
                         tangent_levels, tangent_coarea)
    if (status /= 0 .or. &
        maxval(abs(tangent_coarea)) > 2.0e-11_dp) &
        error stop 'coarea constant gauge derivative'
    ! A C0 P2 polynomial with genuine inter-cell gradient jumps. Its
    ! finite-difference oracle moves face intersections independently;
    ! omitting the Leibniz jump terms must fail this control.
    state%psi = state%psi + &
                0.004_dp*sin(9*(state%points(1, :) - &
                                wanted_axis(1)))* &
                sin(13*(state%points(2, :) - wanted_axis(2)))
    call gs_axis_candidate(state, axis, axis_flux, &
                           hessian, element, status)
    if (status /= 0) error stop 'jump oracle axis unqualified'
    levels = axis_flux*(1 - [(real(i - 1, &
                                   dp)/(size(levels) - 1), i=1, size(levels))])
    tangent_psi = 0.01_dp*cos(11*(state%points(1, :) - &
                                  wanted_axis(1)))* &
                  sin(7*(state%points(2, :) - wanted_axis(2)))
    call gs_element_jet_with_tangent_value()
    tangent_levels = tvalue*(1 - [(real(i - 1, &
                                        dp)/(size(levels) - 1), i=1, size(levels))])
    call gs_level_coarea(state, levels, 256, fine, axis, &
                         status, tangent_psi, &
                         tangent_levels, tangent_coarea)
    if (status /= 0) then
        print *, 'jump derivative status:', status
        error stop 'jump coarea tangent rejected'
    end if
    perturbation_step = 1.0e-5_dp
    plus_state = state; plus_state%psi = state%psi + &
 perturbation_step*tangent_psi
    minus_state = state; minus_state%psi = state%psi - &
 perturbation_step*tangent_psi
    call gs_axis_candidate(plus_state, axis, axis_flux, &
                           hessian, element, status)
    plus_levels = axis_flux*(1 - [(real(i - 1, &
                                        dp)/(size(levels) - 1), i=1, size(levels))])
    call gs_axis_candidate(minus_state, axis, axis_flux, &
                           hessian, element, status)
    minus_levels = axis_flux*(1 - [(real(i - 1, &
                                         dp)/(size(levels) - 1), i=1, size(levels))])
    call gs_level_coarea(plus_state, plus_levels, 256, &
                         plus_coarea, axis, status)
    if (status /= 0) error stop 'jump coarea positive FD rejected'
    call gs_level_coarea(minus_state, minus_levels, 256, &
                         minus_coarea, axis, status)
    if (status /= 0) error stop 'jump coarea negative FD rejected'
    difference_error = maxval(abs((plus_coarea - &
                                   minus_coarea)/(2*perturbation_step) &
                                  - tangent_coarea))
    print *, &
        'native gradient-jump coarea independent FD error:', &
        difference_error
    if (difference_error > 2.0e-7_dp) &
        error stop 'moving-face coarea derivative failed FD oracle'
    call gs_contour_root_response(3.0_dp, 5.0_dp, 2.0_dp, &
                                  response)
    if (abs(response - &
            1) > 1.0e-14_dp) error stop 'equation-level contour root tangent'
    call gs_level_coarea(state, levels, 2, fine, axis, &
                         status)
    if (status == 0) error stop 'underresolved angular control accepted'
    levels(1) = levels(1) - 0.01_dp
    call gs_level_coarea(state, levels, 256, fine, axis, &
                         status)
    if (status == 0) error stop 'non-native axis flux accepted'
    call refined_face_endpoint_oracle()
    call near_vertex_contour_oracle(2)
    call near_vertex_contour_oracle(3)
    print *, &
        'PASS: native P2 exact finite-R coarea and signed source correspondence'
contains
    subroutine refined_face_endpoint_oracle()
        type(gs_state_t)::refined
        real(dp)::boundary(2, 63), center(2), wanted(2), computed_axis(2), best, d, theta
        real(dp)::ell(65), integral(65), derivative(65), exact, radius_squared
        integer::i, cell, stat, amplitude_index, count_index, n
        real(dp)::amplitude
        do i = 1, 63
            theta = 2*acos(-1.0_dp)*(i - 1)/63
            boundary(:, i) = [4.0_dp + 0.4_dp*cos(theta), 0.4_dp*sin(theta)]
        end do
        call initialize_gs(refined, boundary, 0.0003125_dp, stat, 2)
        if (stat /= 0) error stop 'refined polynomial mesh failed'
        best = huge(1.0_dp)
        do cell = 1, size(refined%triangles, 2)
            center = sum(refined%points(:, refined%triangles(1:3, cell)), dim=2)/3
            d = norm2(center - [4.0_dp, 0.0_dp])
            if (d >= best) cycle
            best = d
            wanted = center
        end do
        ! Exactly representable degree two, with all tested circles strictly
        ! inside the polygon. No Grad-Shafranov PDE solve is involved.
        do amplitude_index = 0, 3
            amplitude = 10.0_dp**amplitude_index
            refined%psi = amplitude*(0.25_dp*((refined%points(1, :) - wanted(1))**2 + &
                                             (refined%points(2, :) - wanted(2))**2) - 0.02_dp)
            do count_index = 1, 4
                n = 2**(count_index + 2) + 1
                ell(:n) = [((-0.02_dp + 0.02_dp*real(i - 1, dp)/(n - 1))*amplitude, i=1, n)]
                call gs_level_coarea(refined, ell(:n), 128, integral(:n), computed_axis, stat)
                if (stat /= 0) error stop 'refined exact polynomial primal rejected'
                call gs_level_coarea(refined, ell(:n), 128, integral(:n), computed_axis, stat, &
                                     refined%psi, ell(:n), derivative(:n))
                if (stat /= 0) error stop 'known face endpoint was re-solved outside its cell'
                do i = 1, n
                    radius_squared = 4*(ell(i)/amplitude + 0.02_dp)
                    exact = 4*acos(-1.0_dp)/(amplitude*sqrt(wanted(1)**2 - radius_squared))
                    if (abs(integral(i)/exact - 1) > 2.0e-9_dp) &
                        error stop 'refined exact circular coarea oracle'
                end do
                if (maxval(abs(derivative(:n)/integral(:n) + 1)) > 2.0e-9_dp) &
                    error stop 'refined exact amplitude derivative oracle'
            end do
        end do
    end subroutine refined_face_endpoint_oracle

    subroutine near_vertex_contour_oracle(degree)
        ! A level passing 1e-11 m from an interior mesh vertex crosses the
        ! incident faces at distinct, qualified edge parameters, but within
        ! 1e-10 rad of each other (TC24 readback status 18). The contour is
        ! regular there; its coarea (and P2 tangent) must match the exact circle.
        integer, intent(in)::degree
        type(gs_state_t)::mesh
        real(dp)::boundary(2, 63), center(2), wanted(2), computed_axis(2)
        real(dp)::ell(2), integral(2), derivative(2), theta, best, d, radius
        real(dp)::offsets(4), exact, point(2)
        integer::i, k, cell, vertex, stat, target
        offsets = [1.0e-11_dp, 2.0e-11_dp, -1.0e-11_dp, -2.0e-11_dp]
        do i = 1, 63
            theta = 2*acos(-1.0_dp)*(i - 1)/63
            boundary(:, i) = [4.0_dp + 0.4_dp*cos(theta), 0.4_dp*sin(theta)]
        end do
        call initialize_gs(mesh, boundary, 0.0003125_dp, stat, degree)
        if (stat /= 0) error stop 'near-vertex polynomial mesh failed'
        best = huge(1.0_dp)
        do cell = 1, size(mesh%triangles, 2)
            center = sum(mesh%points(:, mesh%triangles(1:3, cell)), dim=2)/3
            d = norm2(center - [4.0_dp, 0.0_dp])
            if (d >= best) cycle
            best = d
            wanted = center
        end do
        ! Keep the P3 axis off the shared centroid nodes of the cubic cells.
        if (degree == 3) wanted = wanted + [1.0e-3_dp, 7.0e-4_dp]
        mesh%psi = 0.25_dp*((mesh%points(1, :) - wanted(1))**2 + &
                            (mesh%points(2, :) - wanted(2))**2) - 0.02_dp
        ! A generic vertex, and one just above polar angle zero so that the
        ! nearly coincident events open the angle-sorted event list.
        do target = 1, 2
            best = huge(1.0_dp)
            vertex = 0
            do cell = 1, size(mesh%triangles, 2)
                do k = 1, 3
                    point = mesh%points(:, mesh%triangles(k, cell)) - wanted
                    if (target == 1) then
                        d = norm2(point - 0.2_dp*[sqrt(0.5_dp), sqrt(0.5_dp)])
                    else
                        if (point(2) <= 0) cycle
                        d = norm2(point - [0.2_dp, 0.0_dp])
                    end if
                    if (d >= best) cycle
                    best = d
                    vertex = mesh%triangles(k, cell)
                end do
            end do
            radius = norm2(mesh%points(:, vertex) - wanted)
            do i = 1, size(offsets)
                ell = [-0.02_dp, 0.25_dp*(radius + offsets(i))**2 - 0.02_dp]
                call gs_level_coarea(mesh, ell, 128, integral, computed_axis, stat)
                if (stat /= 0) then
                    print *, 'near-vertex coarea status:', stat, ' degree:', &
                        degree, ' vertex:', target, ' offset:', offsets(i)
                    error stop 'regular contour near a mesh vertex rejected'
                end if
                exact = 4*acos(-1.0_dp)/sqrt(wanted(1)**2 - (radius + offsets(i))**2)
                if (abs(integral(2)/exact - 1) > 2.0e-9_dp) &
                    error stop 'near-vertex exact circular coarea oracle'
                if (degree /= 2) cycle
                call gs_level_coarea(mesh, ell, 128, integral, computed_axis, stat, &
                                     mesh%psi, ell, derivative)
                if (stat /= 0) error stop 'near-vertex coarea tangent rejected'
                if (abs(derivative(2)/integral(2) + 1) > 2.0e-9_dp) &
                    error stop 'near-vertex exact amplitude derivative oracle'
            end do
        end do
    end subroutine near_vertex_contour_oracle

    subroutine gs_element_jet_with_tangent_value()
        type(gs_state_t)::tangent_state
        tangent_state = state; tangent_state%psi = tangent_psi
        call gs_element_jet(tangent_state, element, axis(1), &
                            axis(2), tvalue, gradient, &
                            hessian, status)
        if (status /= 0) error stop 'jump tangent jet unqualified'
    end subroutine
end program test_gs_coarea
