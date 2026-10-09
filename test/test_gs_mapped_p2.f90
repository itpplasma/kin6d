program test_gs_mapped_p2
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use kin6d_grad_shafranov, only: gs_state_t, initialize_gs, solve_gs, &
                       gs_linear_response, gs_element_jet, sample_gs, gs_profile_load, &
                          gs_profile_load_jvp, gs_profile_load_vjp, solve_gs_profiles, &
              gs_axis_candidate, gs_classical_current, gs_current_jump, write_gs_mesh, &
                                    GS_UNSUPPORTED_GEOMETRY
    use kin6d_gs_coarea, only: gs_level_coarea
    use fortfem_gauss_quadrature_2d, only: gauss_quadrature_triangle_t
    implicit none
    type(gs_state_t) :: state, bad
    type(gauss_quadrature_triangle_t) :: rule
    real(dp) :: polygon(2, 8), gamma, corner(2), normals(2, 8), beta(8)
    real(dp), parameter :: r0 = 3, level = 5.3145_dp, tau = 0.75_dp
    real(dp), allocatable :: matrix(:, :), load(:), vector(:), response(:), before(:)
    real(dp), allocatable :: action(:), chord_load(:), radius_load(:)
    real(dp) :: vertices(2, 3), affine(2, 2), inverse(2, 2), det, point(2), physical(2)
    real(dp) :: phi(6), dphi(2, 6), gradients(2, 6), jf(2, 2), hf(2, 2, 2)
    real(dp) :: local_matrix(6, 6), local_load(6), weight, source, value, grad(2)
    real(dp) :: hessian(2, 2), want_grad(2), want_h(2, 2), ref_grad(2), ref_h(2, 2)
    real(dp) :: br, bz, bt, axis(2), axis_h(2, 2), current, coarea(2), tangent(2)
    real(dp) :: wrong_gradient(2, 6), wrong_matrix(6, 6), transpose_gap, connection_gap
    real(dp) :: gmu, flux(2), pp(2), ff(2), gflux(2), gpp(2), gff(2), xi(2)
    real(dp) :: x, z, mu0, pprime, ffprime, transformed(2, 6), product(2, 2)
    integer :: i, j, e, k, l, q, status, passed, n, axis_element, unit, ios
    character(32) :: sentinel
    passed = 0
    gamma = 1.44_dp*(sqrt(2.0_dp) - 1)
    do i = 1, 8
        select case (i)
        case (1); corner(1) = -gamma; corner(2) = -1.44_dp
        case (2); corner(1) = gamma; corner(2) = -1.44_dp
        case (3); corner(1) = 1.44_dp; corner(2) = -gamma
        case (4); corner(1) = 1.44_dp; corner(2) = gamma
        case (5); corner(1) = gamma; corner(2) = 1.44_dp
        case (6); corner(1) = -gamma; corner(2) = 1.44_dp
        case (7); corner(1) = -1.44_dp; corner(2) = gamma
        case (8); corner(1) = -1.44_dp; corner(2) = -gamma
        end select
        polygon(1, i) = r0 + (4*corner(1) - 3*corner(2))/25
        polygon(2, i) = (3*corner(1) + 4*corner(2))/25
    end do
    do i = 1, 8
        corner = polygon(:, mod(i, 8) + 1) - polygon(:, i)
        normals(1, i) = corner(2); normals(2, i) = -corner(1)
        beta(i) = dot_product(normals(:, i), polygon(:, i) - [r0, 0.0_dp])
    end do
    call initialize_gs(state, polygon, 0.01_dp, status, 2, r0, level, tau)
    call require(status == 0, 'mapped state initialization')
    call require(state%ready .and. state%mapped, 'mapped ready contract')
    call require(.not. allocated(state%gradient), 'no affine physical gradient cache')
   call require(.not. allocated(state%area), 'no chord area masquerades as curved area')
    call require(size(state%triangles, 2) > 1, 'multiple computational cells')
    call require(maxval(abs(state%points - state%reference_points)) > 0.01_dp, &
                 'distinct computational and physical positions')
    n = size(state%psi)
    allocate (matrix(n, n), load(n), vector(n), before(n), action(n), &
              chord_load(n), radius_load(n))
    matrix = 0; load = 0; chord_load = 0; radius_load = 0
    transpose_gap = 0; connection_gap = 0
    mu0 = 2; pprime = -7.0_dp/5; ffprime = 3.0_dp/7
    call rule%init(5)
    ! Independent forward F implicit derivatives; never call generated G/PR6/PR7.
    ! Quadrature coordinates/weights are the declared native finite-node contract.
    do e = 1, size(state%triangles, 2)
        do j = 1, 3
            vertices(:, j) = state%reference_points(:, state%triangles(j, e))
        end do
        affine(:, 1) = vertices(:, 2) - vertices(:, 1)
        affine(:, 2) = vertices(:, 3) - vertices(:, 1)
        det = affine(1, 1)*affine(2, 2) - affine(1, 2)*affine(2, 1)
        inverse(1, 1) = affine(2, 2)/det
        inverse(2, 1) = -affine(2, 1)/det
        inverse(1, 2) = -affine(1, 2)/det
        inverse(2, 2) = affine(1, 1)/det
        local_matrix = 0; local_load = 0; wrong_matrix = 0
        do q = 1, 7
            xi(1) = rule%xi(q); xi(2) = rule%eta(q)
            point = matmul(affine, xi)
            point = vertices(:, 1) + point
            call inverse_forward(point, physical, jf, hf)
            call shapes(xi, phi, dphi)
            transformed = matmul(transpose(inverse), dphi)
            gradients = matmul(transpose(jf), transformed)
            weight = rule%weights(q)*det &
                     /(jf(1, 1)*jf(2, 2) - jf(1, 2)*jf(2, 1))
            source = mu0*pprime*physical(1) + ffprime/physical(1)
            product = matmul(inverse, jf)
            wrong_gradient = matmul(product, dphi)
            do l = 1, 6
                j = state%triangles(l, e)
                radius_load(j) = radius_load(j) + weight &
                                 *(mu0*pprime*point(1) + ffprime/point(1))*phi(l)
                chord_load(j) = chord_load(j) + rule%weights(q)*det &
                                *(mu0*pprime*point(1) + ffprime/point(1))*phi(l)
                local_load(l) = local_load(l) + weight*source*phi(l)
                do k = 1, 6
                    local_matrix(k, l) = local_matrix(k, l) + weight &
                              /physical(1)*dot_product(gradients(:, k), gradients(:, l))
                    wrong_matrix(k, l) = wrong_matrix(k, l) + weight/physical(1) &
                                *dot_product(wrong_gradient(:, k), wrong_gradient(:, l))
                end do
            end do
        end do
        transpose_gap = max(transpose_gap, maxval(abs(local_matrix - wrong_matrix)))
        call close (maxval(abs(local_matrix - state%stiffness(:, :, e))), &
                    0.0_dp, 3.0e-10_dp, 'independent full cell finite-sum stiffness')
        do l = 1, 6
            j = state%triangles(l, e)
            load(j) = load(j) + local_load(l)
            do k = 1, 6
                i = state%triangles(k, e)
                matrix(i, j) = matrix(i, j) + local_matrix(k, l)
            end do
        end do
    end do
    call require(transpose_gap > 0.01_dp, 'wrong inverse-transpose falsifier')
    call require(maxval(abs(load - radius_load)) > 1.0e-4_dp, &
                 'computational radius in source falsifier')
    call require(maxval(abs(load - chord_load)) > 1.0e-3_dp, &
                 'retained polygon measure and load falsifier')
    vector = 1
    action = matmul(matrix, vector)
    call close (maxval(abs(action)), 0.0_dp, 3.0e-12_dp, &
                'global constants annihilated')
    do i = 1, n
        vector(i) = polynomial(state%reference_points(:, i))
    end do
    before = vector
    where (state%boundary) vector = 0
    action = matmul(matrix, vector)
    where (state%boundary) action = 0
    call gs_linear_response(state, action, response, status)
    call require(status == 0, 'existing CG shared-DOF finite algebra')
    call close (maxval(abs(response - vector)), 0.0_dp, 2.0e-9_dp, &
                'independent assembled action inversion')
    where (state%boundary) action = before
    call gs_linear_response(state, action, response, status)
    call require(status == 0, 'nonzero boundary equation response')
    call require(all(response == before .or. .not. state%boundary), &
                 'equation response preserves exact boundary coefficients')
    state%psi = before
    call solve_gs(state, mu0, pprime, ffprime, status)
    call require(status == 0, 'fixed-source discrete unit solve')
    vector = matmul(matrix, state%psi)
    vector = vector - load
    where (state%boundary) vector = 0
    call close (norm2(vector), 0.0_dp, 3.0e-9_dp, &
                'mapped signed load and complete boundary lifting residual')
    call close (maxval(abs(state%psi - before), mask=state%boundary), &
                0.0_dp, 1.0e-14_dp, 'Dirichlet coefficients retained')
    ! Complete scalar readback of an independently prescribed computational polynomial.
    state%psi = before
    do e = 1, size(state%triangles, 2)
        do j = 1, 3
            vertices(:, j) = state%reference_points(:, state%triangles(j, e))
        end do
        point = 0.21_dp*vertices(:, 1) + 0.32_dp*vertices(:, 2) + 0.47_dp*vertices(:, 3)
        call inverse_forward(point, physical, jf, hf)
        x = point(1) - r0; z = point(2)
        ref_grad(1) = 1 + z; ref_grad(2) = 2 + x + 2*z
        ref_h = reshape([0.0_dp, 1.0_dp, 1.0_dp, 2.0_dp], [2, 2])
        want_grad = matmul(transpose(jf), ref_grad)
        connection_gap = max(connection_gap, &
                         maxval(abs(ref_grad(1)*hf(1, :, :) + ref_grad(2)*hf(2, :, :))))
        product = matmul(ref_h, jf)
        want_h = matmul(transpose(jf), product)
        want_h = want_h + ref_grad(1)*hf(1, :, :) + ref_grad(2)*hf(2, :, :)
        call gs_element_jet(state, e, physical(1), physical(2), value, grad, &
                            hessian, status)
        call require(status == 0, 'physical-addressed full element jet')
        call close (value, polynomial(point), 3.0e-11_dp, 'complete P2 value')
        call close (maxval(abs(grad - want_grad)), 0.0_dp, 3.0e-9_dp, &
                    'forward-chart gradient oracle')
        call close (maxval(abs(hessian - want_h)), 0.0_dp, 4.0e-8_dp, &
                    'forward-chart Hessian connection oracle')
     call sample_gs(state, physical(1), physical(2), -9.0_dp, value, br, bz, bt, status)
        call require(status == 0, 'physical location and field sample')
        call close (br, -want_grad(2)/physical(1), 2.0e-9_dp, 'signed radial field')
        call close (bz, want_grad(1)/physical(1), 2.0e-9_dp, 'signed vertical field')
        call close (bt, -9.0_dp/physical(1), 2.0e-13_dp, 'signed toroidal field')
    end do
    call require(connection_gap > 0.1_dp, 'omitted Hessian connection falsifier')
    ! Independent exact-ray position: physical fraction4/5 on ray(3/10,1/5).
    ! Its ideal rational computational position is (r0+128/625,256/1875).
    point(1) = r0+128.0_dp/625; point(2) = 256.0_dp/1875
    physical(1) = r0+0.24_dp; physical(2) = 0.16_dp
    call sample_gs(state, physical(1), physical(2), -9.0_dp, &
        value, br, bz, bt, status)
    call require(status == 0, 'prescribed rational-ray physical location')
    call close(value, polynomial(point), 3.0e-11_dp, &
        'physical location has independent prescribed computational value')
    ! A global cell-edge query cannot silently choose an affine trace.
    point = 0.5_dp*(state%reference_points(:, state%triangles(1, 1)) &
                    + state%reference_points(:, state%triangles(2, 1)))
    call inverse_forward(point, physical, jf, hf)
    call sample_gs(state, physical(1), physical(2), 1.0_dp, value, br, bz, bt, status)
    call require(status == 2, 'global cell-edge readback unavailable')
    call require(value == 0 .and. br == 0 .and. bz == 0 .and. bt == 0, &
                 'failed sample clears outputs')
    call gs_element_jet(state, 1, ieee_value(0.0_dp, ieee_quiet_nan), &
                        0.0_dp, value, grad, hessian, status)
    call require(status /= 0, 'nonfinite location rejected')
    call require(value == 0 .and. maxval(abs(grad)) == 0 .and. &
                 maxval(abs(hessian)) == 0, 'failed jet clears all components')
    before = state%psi
    call solve_gs(state, 2.0_dp, huge(1.0_dp), 0.0_dp, status)
    call require(status == 3, 'private mapped load arithmetic status propagated')
   call require(all(state%psi == before), 'failed load cannot publish partial solution')
   call require(state%residual == huge(1.0_dp), 'failed solve has unavailable residual')
    flux = [0.0_dp, 1.0_dp]; pp = 1; ff = 0
    call gs_profile_load(state, 1.0_dp, flux, pp, ff, load, status)
    call unavailable(status, 'profile quadrature')
    call require(all(load == 0), 'unavailable profile load cleared')
    call gs_profile_load_jvp(state, 1.0_dp, flux, pp, ff, before, 0.0_dp, &
                             flux, pp, ff, load, status)
    call unavailable(status, 'profile JVP')
    call require(all(load == 0), 'unavailable profile tangent cleared')
    call gs_profile_load_vjp(state, 1.0_dp, flux, pp, ff, before, &
                             load, gmu, gflux, gpp, gff, status)
    call unavailable(status, 'profile VJP')
    call require(all(load == 0) .and. gmu == 0 .and. &
  all(gflux == 0) .and. all(gpp == 0) .and. all(gff == 0), 'profile cotangents cleared')
    call solve_gs_profiles(state, 1.0_dp, flux, pp, ff, status)
    call unavailable(status, 'profile nonlinear solve')
    call gs_axis_candidate(state, axis, value, axis_h, axis_element, status)
    call unavailable(status, 'affine axis search')
    call require(all(axis == 0) .and. axis_element == 0 .and. &
                 value == huge(1.0_dp), 'unavailable axis output policy')
    call gs_classical_current(state, 1, r0, 0.0_dp, 1.0_dp, current, status)
    call unavailable(status, 'classical current')
    call require(current == 0, 'unavailable current cleared')
    call gs_current_jump(state, 1, 2, r0, 0.0_dp, [1.0_dp, 0.0_dp], &
                         1.0_dp, current, status)
    call unavailable(status, 'sheet current straight-edge consumer')
    call gs_level_coarea(state, [0.0_dp, 1.0_dp], 32, coarea, axis, status, &
                         coarea_tangent=tangent)
    call unavailable(status, 'coarea and inverse prerequisite')
    call require(all(coarea == 0) .and. all(tangent == 0), 'coarea outputs cleared')
    open (newunit=unit, file='mapped_export_guard.dat', status='replace')
    write (unit, '(a)') 'must remain unchanged'
    close (unit)
    call write_gs_mesh(state, 'mapped_export_guard.dat', status)
    call unavailable(status, 'affine mesh export')
    open (newunit=unit, file='mapped_export_guard.dat', status='old')
    read (unit, '(a)', iostat=ios) sentinel
    close (unit, status='delete')
    call require(ios == 0 .and. sentinel == 'must remain unchanged', &
                 'export rejected before replacing existing file')
    call initialize_gs(bad, polygon, 0.01_dp, status, 2, mapped_r0=r0)
    call require(status /= 0 .and. .not. bad%ready, 'partial map parameters rejected')
    call gs_linear_response(bad, before, response, status)
    call require(status /= 0 .and. .not. allocated(response), &
                 'failed initialization guarded before solver caches')
    call initialize_gs(bad, polygon, 0.01_dp, status, 1, r0, level, tau)
    call require(status /= 0 .and. .not. bad%ready, 'mapped P1 unsupported')
    call initialize_gs(bad, polygon(:, 8:1:-1), 0.01_dp, status, 2, r0, level, tau)
    call require(status /= 0 .and. .not. bad%ready, 'reversed boundary cycle rejected')
    call require(.not. allocated(bad%psi), 'cycle failure precedes state caches')
    write (*, '(a,i0,a)') 'PASS: ', passed, ' mapped finite-form checks'
contains
    real(dp) function polynomial(point) result(value)
        real(dp), intent(in) :: point(2)
        real(dp) :: x, z
        x = point(1) - r0; z = point(2)
        value = 1 + x + 2*z + x*z + z*z
    end function

    subroutine shapes(xi, phi, gradient)
        real(dp), intent(in) :: xi(2)
        real(dp), intent(out) :: phi(6), gradient(2, 6)
        real(dp) :: bary(3), gb(2, 3)
        integer :: i, j
        bary(1) = 1 - sum(xi); bary(2:3) = xi
        gb = reshape([-1.0_dp, -1.0_dp, 1.0_dp, 0.0_dp, &
                      0.0_dp, 1.0_dp], [2, 3])
        do i = 1, 3
            phi(i) = bary(i)*(2*bary(i) - 1)
            gradient(:, i) = (4*bary(i) - 1)*gb(:, i)
            j = mod(i, 3) + 1
            phi(i + 3) = 4*bary(i)*bary(j)
            gradient(:, i + 3) = 4*(bary(i)*gb(:, j) + bary(j)*gb(:, i))
        end do
    end subroutine

    real(dp) function boundary_level(v) result(value)
        real(dp), intent(in) :: v(2)
        real(dp) :: radius
        radius = r0 + v(1)
        value = (radius*radius - r0*r0)**2 + 4*radius*radius*v(2)**2
    end function

    real(dp) function ray_boundary(v) result(h)
        real(dp), intent(in) :: v(2)
        real(dp) :: lo, hi, mid, probe(2)
        integer :: i
        lo = 0; hi = 1
        do
            probe = hi*v
            if (boundary_level(probe) >= level) exit
            hi = 2*hi
        end do
        do i = 1, 100
            mid = (lo + hi)/2
            probe = mid*v
            if (boundary_level(probe) < level) then
                lo = mid
            else
                hi = mid
            end if
        end do
        h = (lo + hi)/2
    end function

    subroutine inverse_forward(point, physical, jf, hf)
        real(dp), intent(in) :: point(2)
        real(dp), intent(out) :: physical(2), jf(2, 2), hf(2, 2, 2)
        real(dp) :: v(2), h, p, lo, hi, t, scale, probe(2), scores(8)
        integer :: i
        v = point - [r0, 0.0_dp]
        probe = v/(1 - tau)
        if (boundary_level(probe) <= level) then
            physical = point
        else
            h = ray_boundary(v)
            scores = matmul(transpose(normals), v)
            p = 1/maxval(scores/beta)
            lo = 1; hi = h
            ! F(t v)=v, solved by independent scalar bisection, not inverse G leaf.
            do i = 1, 100
                t = (lo + hi)/2
                scale = (1 - t/h + p*t/h**2 - (1 - tau)*p/h)/tau
                if (t*scale < 1) then
                    lo = t
                else
                    hi = t
                end if
            end do
            physical = [r0, 0.0_dp] + ((lo + hi)/2)*v
        end if
        call forward_jet(physical, jf, hf)
    end subroutine

    subroutine forward_jet(physical, jacobian, hessian)
        real(dp), intent(in) :: physical(2)
        real(dp), intent(out) :: jacobian(2, 2), hessian(2, 2, 2)
        real(dp) :: v(2), boundary(2), radius, h, p, g(2), k(2, 2), denominator
        real(dp) :: hi(2), hij(2, 2), pi(2), pij(2, 2), wi(2, 2), normal(2)
        real(dp) :: sh, sp, shh, shp, scale, si(2), sij(2, 2), q
        real(dp) :: probe(2), scores(8), kw(2)
        integer :: a, i, j, facet
        v = physical - [r0, 0.0_dp]
        jacobian = 0; hessian = 0
        jacobian(1, 1) = 1; jacobian(2, 2) = 1
        q = 1 - tau
        probe = v/q
        if (boundary_level(probe) <= level) return
        h = ray_boundary(v)
        boundary = [r0, 0.0_dp] + h*v
        radius = boundary(1)
        g(1) = 4*radius*(radius*radius - r0*r0 + 2*boundary(2)**2)
        g(2) = 8*radius*radius*boundary(2)
        k(1, 1) = 12*radius*radius - 4*r0*r0 + 8*boundary(2)**2
        k(1, 2) = 16*radius*boundary(2); k(2, 1) = k(1, 2)
        k(2, 2) = 8*radius*radius
        denominator = dot_product(g, v)
        hi = -h*g/denominator
        do i = 1, 2
            wi(:, i) = v*hi(i)
            wi(i, i) = wi(i, i) + h
        end do
        do j = 1, 2
            do i = 1, 2
                kw = matmul(k, wi(:, j))
                hij(i, j) = -(dot_product(wi(:, i), kw) &
                              + g(i)*hi(j) + g(j)*hi(i))/denominator
            end do
        end do
        scores = matmul(transpose(normals), v)
        facet = maxloc(scores/beta, dim=1)
        normal = normals(:, facet)/beta(facet)
        p = 1/dot_product(normal, v)
        pi = -p*p*normal
        do j = 1, 2
            pij(:, j) = 2*p**3*normal*normal(j)
        end do
        scale = (1 - 1/h + p/h**2 - q*p/h)/tau
        sh = (1/h**2 - 2*p/h**3 + q*p/h**2)/tau
        sp = (1/h**2 - q/h)/tau
        shh = (-2/h**3 + 6*p/h**4 - 2*q*p/h**3)/tau
        shp = (-2/h**3 + q/h**2)/tau
        si = sh*hi + sp*pi
        do j = 1, 2
            do i = 1, 2
                sij(i, j) = sh*hij(i, j) + sp*pij(i, j) + shh*hi(i)*hi(j) &
                            + shp*(hi(i)*pi(j) + pi(i)*hi(j))
            end do
        end do
        jacobian = 0
        do a = 1, 2
            jacobian(a, :) = v(a)*si
            jacobian(a, a) = jacobian(a, a) + scale
            hessian(a, :, :) = v(a)*sij
            hessian(a, a, :) = hessian(a, a, :) + si
            hessian(a, :, a) = hessian(a, :, a) + si
        end do
    end subroutine

    subroutine unavailable(status, name)
        integer, intent(in) :: status
        character(*), intent(in) :: name
        call require(status == GS_UNSUPPORTED_GEOMETRY, name)
    end subroutine

    subroutine close (value, expected, tolerance, name)
        real(dp), intent(in) :: value, expected, tolerance
        character(*), intent(in) :: name
        call require(abs(value - expected) <= tolerance, name)
    end subroutine

    subroutine require(condition, name)
        logical, intent(in) :: condition
        character(*), intent(in) :: name
        if (.not. condition) then
            print *, 'FAIL: ', name
            error stop 1
        end if
        passed = passed + 1
    end subroutine
end program
