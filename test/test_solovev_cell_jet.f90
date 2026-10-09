program test_solovev_cell_jet
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use kin6d_solovev_geometry, only: solovev_geometry_t, geometry_sample_t, &
        initialize_solovev_geometry, geometry_ok, geometry_invalid, &
        geometry_interface, geometry_arithmetic
    use kin6d_solovev_cell_jet, only: solovev_cell_scalar_jet
    use fortfem_geometry_jet_2d, only: transform_scalar_jet_2d
    implicit none
    type(solovev_geometry_t) :: map, uninitialized
    type(geometry_sample_t) :: geometry
    real(dp) :: polygon(2, 8), gamma, vertex(2), r0, level, tau
    real(dp) :: vertices(2, 3), m(2, 2), xi(2), point(2), want_point(2)
    real(dp) :: jg(2, 2), hg(2, 2, 2), jt(2, 2), ht(2, 2, 2)
    real(dp) :: gradient_reference(2), hessian_reference(2, 2)
    real(dp) :: gradient(2), hessian(2, 2), want_gradient(2), want_hessian(2, 2)
    real(dp) :: determinant, inverse(2, 2), bad_gradient(2), bad_hessian(2, 2)
    real(dp) :: reordered(2, 2, 2), wrong_order(2, 2, 2), bad_det, bad_inverse(2, 2)
    real(dp) :: br, bz, current, want_current, mu0, swapped(2, 3)
    integer :: status, i, a, j, k, passed
    passed = 0
    r0 = 3; level = 5.3145_dp; tau = 0.75_dp
    gamma = 1.44_dp*(sqrt(2.0_dp)-1)
    do i = 1, 8
        select case (i)
        case (1); vertex = [-gamma, -1.44_dp]
        case (2); vertex = [gamma, -1.44_dp]
        case (3); vertex = [1.44_dp, -gamma]
        case (4); vertex = [1.44_dp, gamma]
        case (5); vertex = [gamma, 1.44_dp]
        case (6); vertex = [-gamma, 1.44_dp]
        case (7); vertex = [-1.44_dp, gamma]
        case (8); vertex = [-1.44_dp, -gamma]
        end select
        polygon(:, i) = [r0, 0.0_dp] + &
            [4*vertex(1)-3*vertex(2), 3*vertex(1)+4*vertex(2)]/25
    end do
    call initialize_solovev_geometry(map, r0, level, tau, polygon, status)
    call require(status == geometry_ok, 'noncircular fixture initialized')
    xi = [0.25_dp, 0.125_dp]
    m = reshape([1.0_dp/64, 1.0_dp/256, 1.0_dp/96, 3.0_dp/256], [2, 2])
    point = [r0+128.0_dp/625, 256.0_dp/1875]
    want_point = [r0+0.24_dp, 0.16_dp]
    call cell_at(point)
    ! Exact independent formal-series oracle for G at this quartic/facet point.
    jg = reshape([6885.0_dp/4864, 395.0_dp/2432, &
        2145.0_dp/9728, 6415.0_dp/4864], [2, 2])
    hg(1, :, :) = reshape([8399001875.0_dp/3034202112.0_dp, &
        12750974375.0_dp/4045602816.0_dp, 12750974375.0_dp/4045602816.0_dp, &
        -286488125.0_dp/674267136], [2, 2])
    hg(2, :, :) = reshape([73758125.0_dp/284456448, &
        15614155625.0_dp/6068404224.0_dp, 15614155625.0_dp/6068404224.0_dp, &
        3782714375.0_dp/2022801408.0_dp], [2, 2])
    jt = matmul(jg, m)
    do a = 1, 2
        ht(a, :, :) = matmul(transpose(m), matmul(hg(a, :, :), m))
    end do
    call physical_polynomial(want_point, want_gradient, want_hessian)
    gradient_reference = matmul(transpose(jt), want_gradient)
    hessian_reference = matmul(transpose(jt), matmul(want_hessian, jt)) &
        + want_gradient(1)*ht(1, :, :) + want_gradient(2)*ht(2, :, :)
    call evaluate(map)
    call require(status == geometry_ok, 'annular scalar-jet status')
    call require(geometry%facet == 3, 'facet tag retained')
    call require(.not. geometry%identity, 'annular identity flag retained')
    call require(.not. geometry%certified, 'composition remains uncertified')
    call close(maxval(abs(geometry%value-want_point)), 0.0_dp, 3.0e-13_dp, &
        'composed physical point')
    call close(maxval(abs(geometry%jacobian-jt)), 0.0_dp, 1.0e-13_dp, &
        'nonsymmetric affine-cell Jacobian composition')
    call close(maxval(abs(geometry%hessian-ht)), 0.0_dp, 1.0e-13_dp, &
        'component-first affine-cell Hessian composition')
    call close(maxval(abs(gradient-want_gradient)), 0.0_dp, 3.0e-11_dp, &
        'physical polynomial gradient')
    call close(maxval(abs(hessian-want_hessian)), 0.0_dp, 3.0e-9_dp, &
        'physical polynomial Hessian with connection')
    call close(maxval(abs(matmul(geometry%jacobian, inverse) &
        -reshape([1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2, 2]))), &
        0.0_dp, 3.0e-14_dp, 'inverse convention')
    mu0 = 1.2566370614359173e-6_dp
    br = -gradient(2)/geometry%value(1); bz = gradient(1)/geometry%value(1)
    current = -(hessian(1, 1)+hessian(2, 2)-gradient(1)/geometry%value(1)) &
        /(mu0*geometry%value(1))
    want_current = -(6-want_gradient(1)/want_point(1))/(mu0*want_point(1))
    call close(br, -want_gradient(2)/want_point(1), 2.0e-11_dp, 'signed B_R')
    call close(bz, want_gradient(1)/want_point(1), 2.0e-11_dp, 'signed B_Z')
    call close(current, want_current, 2.0e-3_dp, 'signed cylindrical current')
    ! Behavioral falsifiers: wrong layout and an omitted connection term.
    do a = 1, 2
        do j = 1, 2
            do k = 1, 2
                reordered(j, k, a) = ht(a, j, k)
                wrong_order(j, k, a) = ht(j, k, a)
            end do
        end do
    end do
    call transform_scalar_jet_2d(want_point, jt, 0*reordered, &
        gradient_reference, hessian_reference, bad_det, bad_inverse, &
        bad_gradient, bad_hessian, status)
    call require(status == 0, 'omitted connection is finite but wrong')
    call require(maxval(abs(bad_hessian-want_hessian)) > 0.1_dp, &
        'nonaffine curvature oracle rejects omitted connection')
    call transform_scalar_jet_2d(want_point, jt, wrong_order, &
        gradient_reference, hessian_reference, bad_det, bad_inverse, &
        bad_gradient, bad_hessian, status)
    call require(status == 0, 'wrong tensor permutation is finite but wrong')
    call require(maxval(abs(bad_hessian-want_hessian)) > 0.1_dp, &
        'physical polynomial rejects wrong component order')
    ! Identity G composed with a nonsymmetric affine cell remains affine.
    call cell_at([r0, 0.0_dp])
    call physical_polynomial([r0, 0.0_dp], want_gradient, want_hessian)
    gradient_reference = matmul(transpose(m), want_gradient)
    hessian_reference = matmul(transpose(m), matmul(want_hessian, m))
    call evaluate(map)
    call require(status == geometry_ok .and. geometry%identity, 'axis identity branch')
    call close(maxval(abs(geometry%jacobian-m)), 0.0_dp, 2.0e-15_dp, &
        'identity G gives cell Jacobian, not identity matrix')
    call close(sum(abs(geometry%hessian)), 0.0_dp, 0.0_dp, 'identity curvature zero')
    call close(maxval(abs(hessian-want_hessian)), 0.0_dp, 2.0e-11_dp, &
        'identity branch physical Hessian')
    call cell_at([r0+0.075_dp, 0.05_dp])
    call evaluate(map)
    call require(status == geometry_interface, 'inner cut status propagated')
    call require(geometry%inner_interface, 'inner cut tag retained')
    call cleared('inner cut')
    call cell_at([r0, 0.0_dp]+0.98_dp*(polygon(:, 4)-[r0, 0.0_dp]))
    call evaluate(map)
    call require(status == geometry_interface, 'sector cut status propagated')
    call require(geometry%sector_interface, 'sector cut tag retained')
    call cleared('sector cut')
    call cell_at([r0, 0.0_dp])
    swapped = vertices; vertices(:, 2) = swapped(:, 3); vertices(:, 3) = swapped(:, 2)
    call evaluate(map)
    call require(status == geometry_invalid, 'reversed cell rejected')
    call cleared('reversed cell')
    vertices = swapped; vertices(:, 3) = vertices(:, 2)
    call evaluate(map)
    call require(status == geometry_invalid, 'degenerate cell rejected')
    call cleared('degenerate cell')
    call cell_at([r0, 0.0_dp]); xi = [-0.1_dp, 0.1_dp]
    call evaluate(map)
    call require(status == geometry_invalid, 'outside reference triangle rejected')
    call cleared('outside reference')
    xi = [0.25_dp, 0.125_dp]
    gradient_reference(1) = ieee_value(0.0_dp, ieee_quiet_nan)
    call evaluate(map)
    call require(status == geometry_invalid, 'nonfinite scalar jet rejected')
    call cleared('nonfinite scalar jet')
    gradient_reference = huge(1.0_dp); hessian_reference = 0
    call evaluate(map)
    call require(status == geometry_arithmetic, &
        'finite scalar overflow is arithmetic failure')
    call cleared('scalar arithmetic failure')
    gradient_reference = 0; hessian_reference = 0
    call evaluate(uninitialized)
    call require(status == geometry_invalid, 'uninitialized geometry propagated')
    call cleared('uninitialized map')
    print *, 'Solovev affine-cell scalar jet checks:', passed
contains
    subroutine cell_at(center)
        real(dp), intent(in) :: center(2)
        vertices(:, 1) = center-matmul(m, xi)
        vertices(:, 2) = vertices(:, 1)+m(:, 1)
        vertices(:, 3) = vertices(:, 1)+m(:, 2)
    end subroutine
    subroutine evaluate(selected_map)
        type(solovev_geometry_t), intent(in) :: selected_map
        call solovev_cell_scalar_jet(selected_map, vertices, xi, &
            gradient_reference, hessian_reference, geometry, determinant, &
            inverse, gradient, hessian, status)
    end subroutine
    subroutine physical_polynomial(x, g, h)
        real(dp), intent(in) :: x(2)
        real(dp), intent(out) :: g(2), h(2, 2)
        ! u=R^2+3RZ+2Z^2+5R-7Z+11, independent Cartesian polynomial.
        g = [2*x(1)+3*x(2)+5, 3*x(1)+4*x(2)-7]
        h = reshape([2.0_dp, 3.0_dp, 3.0_dp, 4.0_dp], [2, 2])
    end subroutine
    subroutine cleared(label)
        character(*), intent(in) :: label
        call close(sum(abs(geometry%value))+sum(abs(geometry%jacobian)) &
            +sum(abs(geometry%hessian))+abs(determinant)+sum(abs(inverse)) &
            +sum(abs(gradient))+sum(abs(hessian)), 0.0_dp, 0.0_dp, label//' clears')
        call require(geometry%status == status, label//' status agrees')
        call require(.not. geometry%certified, label//' remains uncertified')
    end subroutine
    subroutine require(ok, label)
        logical, intent(in) :: ok
        character(*), intent(in) :: label
        if (.not. ok) then
            print *, 'FAIL: ', label
            error stop 1
        end if
        passed = passed+1
    end subroutine
    subroutine close(actual, expected, tolerance, label)
        real(dp), intent(in) :: actual, expected, tolerance
        character(*), intent(in) :: label
        call require(abs(actual-expected) <= tolerance, label)
    end subroutine
end program
