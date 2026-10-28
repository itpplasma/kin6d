program test_ideal_map_field_angle
    ! Field-line angle of the frozen field (docs/ideal-seed-3d.md):
    ! 1. generated leaves against finite differences and JVP/VJP pairing,
    ! 2. map energy gradient and Hessian against energy/gradient differences
    !    in node and angle directions, at the seed and after Newton steps,
    ! 3. regression: the axisymmetric sheared oracle boundary is supplied in
    !    its chart angle, which is not the straight-field-line angle. With the
    !    wall pinned in that angle the field is wrong at the percent level;
    !    with the field angle it agrees with the exact field.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_ideal_field_angle_leaves, only: field_angle_mode, field_angle_mode_jvp, &
        field_angle_mode_vjp, field_angle_form, field_angle_form_jvp, &
        field_angle_form_vjp
    use kin6d_ideal_map_leaves, only: ideal_point_gradient_jets, &
        ideal_point_field_gradient_jets, ideal_point_field_jets, &
        ideal_axisymmetric_field
    use kin6d_ideal_map, only: ideal_map_t, ideal_map_energy, ideal_map_hessian, &
        ideal_map_reduce, ideal_map_expand, ideal_map_solve, &
        ideal_map_solve_report_t, ideal_map_point, ideal_map_angle_rdof, IDEAL_MAP_OK
    use kin6d_ideal_seed_3d, only: toroidal_seed_t, toroidal_profile_t, &
        ideal_seed_3d_init, ideal_seed_disk_mesh, ideal_seed_3d_angle_density
    use kin6d_nested_oracles, only: nested_oracle_t, NESTED_SHEARED
    use kin6d_nested_seed_input, only: nested_seed_boundary, nested_seed_profile, &
        NESTED_SEED_OK
    use kin6d_grad_shafranov, only: gs_state_t
    use kin6d_triangle_basis, only: triangle_basis_jet, TRIANGLE_MAX_NODES
    use fortfem_element_to_csr, only: csr_matrix_t
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp)
    type(nested_oracle_t) :: oracle
    type(toroidal_profile_t) :: profile
    real(dp), allocatable :: boundary(:, :, :)
    real(dp) :: flux, fit_error, av(4), exact_energy, error_off, error_on
    real(dp) :: energy_off, energy_on
    integer :: status
    call check_leaves()
    oracle%family = NESTED_SHEARED; oracle%eps = 0; oracle%delta = 0.14_dp
    call nested_seed_profile(oracle, 12, 1e-9_dp, profile, flux, fit_error, status)
    if (status /= NESTED_SEED_OK) error stop 'sheared profile'
    call nested_seed_boundary(oracle, 8, 1, boundary, status)
    if (status /= NESTED_SEED_OK) error stop 'sheared boundary'
    call check_map_derivatives()
    call oracle%averages(48, av)
    exact_energy = av(1)*(av(3)/2 - av(2))/2
    call regression(-1, energy_off, error_off)
    call regression(4, energy_on, error_on)
    print '(a,2es12.4)', 'chart-angle wall, B rms error pinned / field angle:', &
        error_off, error_on
    print '(a,2es12.4)', 'relative energy above exact, pinned / field angle:', &
        energy_off/exact_energy - 1, energy_on/exact_energy - 1
    if (.not. error_off > 2e-2_dp) error stop 'pinned chart-angle wall should be wrong'
    if (.not. error_on < 2e-3_dp) error stop 'field angle: B error'
    ! Conforming: every discrete state is an admissible ideal state.
    if (.not. energy_on > exact_energy*(1 - 1e-7_dp)) error stop 'energy below exact'
    if (.not. energy_on/exact_energy - 1 < (energy_off/exact_energy - 1)/10) &
        error stop 'field angle: energy'
    print '(a)', 'PASS ideal map field angle'
contains
    subroutine check_leaves()
        real(dp) :: x(8), d(8), y(3), yp(3), ym(3), dy(3), wt(3), g(8), h
        real(dp) :: f(6), df(6), z(3), zp(3), zm(3), dz(3), wz(3), gf(6)
        real(dp) :: d1(3, 3), d2(3, 3, 3), d3(3, 3, 3, 3), qd(3, 3), b(3), db(3)
        real(dp) :: w, wp, wm, jac, gv(12), gvp(12), gvm(12), gb(3), gbp(3), gbm(3)
        real(dp) :: hbb(6), hvb(12, 3), hm(3, 3), worst
        integer :: k
        h = 1e-6_dp
        x = [0.31_dp, -0.47_dp, 3.0_dp, 2.0_dp, cos(0.7_dp), sin(0.7_dp), 0.23_dp, &
            -0.11_dp]
        d = [0.2_dp, -0.3_dp, 0.0_dp, 0.0_dp, -0.4_dp, 0.5_dp, 0.7_dp, 0.9_dp]
        call field_angle_mode(x(1), x(2), x(3), x(4), x(5), x(6), x(7), x(8), &
            y(1), y(2), y(3))
        call field_angle_mode(x(1) + h*d(1), x(2) + h*d(2), x(3), x(4), x(5) + h*d(5), &
            x(6) + h*d(6), x(7) + h*d(7), x(8) + h*d(8), yp(1), yp(2), yp(3))
        call field_angle_mode(x(1) - h*d(1), x(2) - h*d(2), x(3), x(4), x(5) - h*d(5), &
            x(6) - h*d(6), x(7) - h*d(7), x(8) - h*d(8), ym(1), ym(2), ym(3))
        call field_angle_mode_jvp(x(1), x(2), x(3), x(4), x(5), x(6), x(7), x(8), &
            d(1), d(2), d(3), d(4), d(5), d(6), d(7), d(8), dy(1), dy(2), dy(3))
        if (maxval(abs((yp - ym)/(2*h) - dy)) > 1e-8_dp) error stop 'mode JVP'
        wt = [0.3_dp, -1.1_dp, 0.6_dp]
        call field_angle_mode_vjp(x(1), x(2), x(3), x(4), x(5), x(6), x(7), x(8), &
            wt(1), wt(2), wt(3), g(1), g(2), g(3), g(4), g(5), g(6), g(7), g(8))
        if (abs(dot_product(g, d) - dot_product(wt, dy)) > 1e-12_dp) error stop 'mode VJP'
        f = [0.4_dp, -0.2_dp, 0.3_dp, -0.6_dp, 1.7_dp, pi]
        df = [0.5_dp, 0.1_dp, -0.8_dp, 0.3_dp, 0.0_dp, 0.0_dp]
        call field_angle_form(f(1), f(2), f(3), f(4), f(5), f(6), z(1), z(2), z(3))
        call field_angle_form(f(1) + h*df(1), f(2) + h*df(2), f(3) + h*df(3), &
            f(4) + h*df(4), f(5), f(6), zp(1), zp(2), zp(3))
        call field_angle_form(f(1) - h*df(1), f(2) - h*df(2), f(3) - h*df(3), &
            f(4) - h*df(4), f(5), f(6), zm(1), zm(2), zm(3))
        call field_angle_form_jvp(f(1), f(2), f(3), f(4), f(5), f(6), df(1), df(2), &
            df(3), df(4), df(5), df(6), dz(1), dz(2), dz(3))
        if (maxval(abs((zp - zm)/(2*h) - dz)) > 1e-8_dp) error stop 'form JVP'
        wz = [-0.7_dp, 0.2_dp, 1.3_dp]
        call field_angle_form_vjp(f(1), f(2), f(3), f(4), f(5), f(6), wz(1), wz(2), &
            wz(3), gf(1), gf(2), gf(3), gf(4), gf(5), gf(6))
        if (abs(dot_product(gf, df) - dot_product(wz, dz)) > 1e-12_dp) error stop 'form VJP'
        ! Energy derivatives in b on a cylindrical chart jet with a sheared Q.
        d1 = reshape([1.1_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.9_dp, 0.0_dp, &
            2.3_dp, 0.0_dp], [3, 3])
        d2 = 0; d2(2, 1, 3) = 1; d2(2, 3, 1) = 1; d2(1, 3, 3) = -2.3_dp
        d3 = 0; d3(2, 3, 3, 3) = -2.3_dp
        qd = reshape([1.2_dp, 0.1_dp, 0.0_dp, -0.2_dp, 0.8_dp, 0.0_dp, 0.05_dp, &
            0.02_dp, 1.0_dp], [3, 3])
        b = [0.3_dp, -0.2_dp, 1.4_dp]
        call ideal_point_field_jets(d1, d2, d3, qd, b, 0.2_dp, 1.0_dp, 1.0_dp, gb, &
            hbb, hvb)
        hm = reshape([hbb(1), hbb(2), hbb(4), hbb(2), hbb(3), hbb(5), hbb(4), hbb(5), &
            hbb(6)], [3, 3])
        worst = 0
        do k = 1, 3
            db = 0; db(k) = h
            call ideal_point_gradient_jets(d1, d2, d3, qd, b + db, 0.2_dp, 1.0_dp, &
                1.0_dp, wp, jac, gvp)
            call ideal_point_gradient_jets(d1, d2, d3, qd, b - db, 0.2_dp, 1.0_dp, &
                1.0_dp, wm, jac, gvm)
            call ideal_point_field_gradient_jets(d1, d2, d3, qd, b + db, 0.2_dp, &
                1.0_dp, 1.0_dp, gbp)
            call ideal_point_field_gradient_jets(d1, d2, d3, qd, b - db, 0.2_dp, &
                1.0_dp, 1.0_dp, gbm)
            worst = max(worst, abs((wp - wm)/(2*h) - gb(k)), &
                maxval(abs((gbp - gbm)/(2*h) - hm(:, k))), &
                maxval(abs((gvp - gvm)/(2*h) - hvb(:, k))))
        end do
        call ideal_point_gradient_jets(d1, d2, d3, qd, b, 0.2_dp, 1.0_dp, 1.0_dp, w, &
            jac, gv)
        if (.not. abs(w) > 0) error stop 'leaf energy'
        print '(a,es11.3)', 'field leaves FD worst:', worst
        if (worst > 1e-7_dp) error stop 'field leaves FD'
    end subroutine

    subroutine setup(map, seed, angle_mpol, ntor)
        type(ideal_map_t), intent(out) :: map
        type(toroidal_seed_t), intent(out) :: seed
        integer, intent(in) :: angle_mpol, ntor
        type(gs_state_t) :: mesh
        integer :: st
        call ideal_seed_disk_mesh(mesh, 8, 3, st)
        if (st /= 0) error stop 'mesh'
        call ideal_seed_3d_init(map, seed, mesh, boundary(:, -ntor:ntor, :), 2, &
            flux, profile, st, angle_mpol=angle_mpol, angle_ntor=ntor)
        if (st /= IDEAL_MAP_OK) error stop 'seed'
    end subroutine

    subroutine check_map_derivatives()
        ! Directional checks in the reduced space (node and angle dofs).
        type(ideal_map_t) :: map
        type(toroidal_seed_t) :: seed
        type(ideal_map_solve_report_t) :: report
        integer :: pass, st
        call setup(map, seed, 3, 1)
        if (map%nangle /= 2*(1 + 3*3)) error stop 'angle dof count'
        do pass = 1, 2
            call directional(map, pass)
            if (pass == 2) exit
            ! A few Newton steps (not converged) move off the seed.
            call ideal_map_solve(map, 1e-3_dp, 3, report, st)
            if (maxval(abs(map%angle)) < 1e-3_dp) error stop 'angle did not move'
        end do
    end subroutine

    subroutine directional(map, pass)
        type(ideal_map_t), intent(inout) :: map
        integer, intent(in) :: pass
        type(csr_matrix_t) :: matrix
        real(dp), allocatable :: gred(:), d(:), full(:), g(:), gp(:), gm(:), hd(:)
        real(dp), allocatable :: x0(:, :, :), a0(:), rp(:), rm(:)
        real(dp) :: e0, ep, em, h, slope, curvature
        integer :: st, k, i
        allocate (g(map%ndof), gp(map%ndof), gm(map%ndof), full(map%ndof), &
            gred(map%nred), d(map%nred), hd(map%nred), rp(map%nred), rm(map%nred))
        x0 = map%x; a0 = map%angle
        do k = 1, map%nred
            d(k) = sin(1.3_dp*k + pass)
        end do
        where (map%fixed) d = 0
        ! Angle directions carry radians, node directions reference lengths.
        do k = 1, map%nangle
            i = ideal_map_angle_rdof(map, k)
            d(i) = 20*d(i)
        end do
        call ideal_map_expand(map, d, full)
        call ideal_map_energy(map, map%x, e0, st, g)
        if (st /= IDEAL_MAP_OK) error stop 'energy'
        call ideal_map_reduce(map, g, gred)
        call ideal_map_hessian(map, map%x, matrix)
        hd = 0
        do k = 1, map%nred
            do i = matrix%rowptr(k), matrix%rowptr(k + 1) - 1
                hd(k) = hd(k) + matrix%val(i)*d(matrix%colidx(i))
            end do
        end do
        ! Central differences: the error is O(h^2) (1.4e-6 at h = 1e-5, 1.4e-8 at 1e-6).
        h = 1e-6_dp
        call shifted(map, x0, a0, full, h, ep, gp)
        call shifted(map, x0, a0, full, -h, em, gm)
        call ideal_map_reduce(map, gp, rp)
        call ideal_map_reduce(map, gm, rm)
        slope = abs((ep - em)/(2*h) - dot_product(gred, d))/abs(dot_product(gred, d))
        curvature = maxval(abs((rp - rm)/(2*h) - hd), mask=.not. map%fixed) &
            /maxval(abs(hd), mask=.not. map%fixed)
        print '(a,i0,2es11.3)', 'map derivatives (pass, slope, Hessian action): ', &
            pass, slope, curvature
        if (slope > 1e-6_dp) error stop 'gradient vs energy difference'
        if (curvature > 1e-6_dp) error stop 'Hessian vs gradient difference'
        map%x = x0; map%angle = a0
    end subroutine

    subroutine shifted(map, x0, a0, full, t, energy, grad)
        type(ideal_map_t), intent(inout) :: map
        real(dp), intent(in) :: x0(:, :, :), a0(:), full(:), t
        real(dp), intent(out) :: energy, grad(:)
        integer :: st
        map%x = x0 + t*reshape(full(1:size(x0)), shape(x0))
        map%angle = a0 + t*full(size(x0) + 1:)
        call ideal_map_energy(map, map%x, energy, st, grad)
        if (st /= IDEAL_MAP_OK) error stop 'shifted energy'
    end subroutine

    subroutine regression(angle_mpol, energy, error)
        ! Axisymmetric (ntor 0) solve on the P3 disk mesh with 8 boundary nodes,
        ! B sampled at interior points of every cell against the exact field.
        integer, intent(in) :: angle_mpol
        real(dp), intent(out) :: energy, error
        type(ideal_map_t) :: map
        type(toroidal_seed_t) :: seed
        type(ideal_map_solve_report_t) :: report
        real(dp) :: xi(2, 3), shape(TRIANGLE_MAX_NODES), grad(2, TRIANGLE_MAX_NODES)
        real(dp) :: q(3), y(2), dg(2, 2), a(2, 2), chi, cg(2), b3, p, dcg(2), db3
        real(dp) :: bfield(3), jac, psi(2), exact(20), sums(2)
        real(dp), allocatable :: g(:)
        integer :: st, e, k, n
        call setup(map, seed, angle_mpol, 0)
        call ideal_map_solve(map, 1e-8_dp, 60, report, st)
        if (st /= IDEAL_MAP_OK) error stop 'regression solve'
        allocate (g(map%ndof))
        call ideal_map_energy(map, map%x, energy, st, g)
        xi = reshape([1/3.0_dp, 1/3.0_dp, 0.2_dp, 0.6_dp, 0.6_dp, 0.2_dp], [2, 3])
        n = map%nloc; sums = 0
        do e = 1, map%cells
            do k = 1, 3
                call triangle_basis_jet(xi(:, k), shape(:n), grad(:, :n))
                call ideal_map_point(map, e, xi(:, k), 0.0_dp, q, y)
                dg = matmul(map%reference(:, map%conn(:, e)), transpose(grad(:, :n)))
                a = matmul(map%x(1:2, 1, map%conn(:, e)), transpose(grad(:, :n)))
                call seed%evaluate(y, chi, cg, b3, p)
                call ideal_seed_3d_angle_density(map, seed, y, 0.0_dp, dcg, db3)
                call ideal_axisymmetric_field(a(1, 1), a(2, 1), a(1, 2), a(2, 2), &
                    dg(1, 1), dg(2, 1), dg(1, 2), dg(2, 2), 0.0_dp, 0.0_dp, &
                    cg(1) + dcg(1), cg(2) + dcg(2), b3 + db3, q(1), bfield(1), &
                    bfield(2), bfield(3), psi(1), psi(2), jac)
                ! The leaf divides by R det(dq/ds), while J = sigma det(DX) and
                ! the cylindrical (R, Z, phi) chart has det = -R.
                bfield = -map%sigma*bfield
                call oracle%point([q(1), 0.0_dp, q(2)], exact, st)
                if (st /= 0) error stop 'oracle point'
                sums = sums + [sum((bfield - exact([3, 5, 4]))**2), sum(exact(3:5)**2)]
            end do
        end do
        error = sqrt(sums(1)/sums(2))
    end subroutine
end program
