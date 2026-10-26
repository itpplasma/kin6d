program test_ideal_map_wall_slide
    ! Sliding fixed boundary (ideal_map_wall_slide): the generated wall-mode
    ! leaf against central differences and its own JVP/VJP pairing, and the
    ! reduced gradient and Hessian (wall Jacobian plus curvature term) of a
    ! 3D shaped seed against differences of the energy along the retraction
    ! that the Newton solve uses (interior linear, wall nodes placed on the
    ! wall at the stepped angles).
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_ideal_seed_3d_leaves, only: seed_wall_mode, seed_wall_mode_jvp, &
        seed_wall_mode_vjp
    use kin6d_ideal_map, only: ideal_map_t, ideal_map_energy, ideal_map_reduce, &
        ideal_map_expand, ideal_map_hessian, ideal_map_wall_set, ideal_map_rdof, &
        IDEAL_MAP_OK
    use kin6d_ideal_seed_3d, only: toroidal_profile_t, toroidal_seed_t, &
        ideal_seed_3d_init, ideal_seed_disk_mesh
    use kin6d_grad_shafranov, only: gs_state_t
    use fortfem_element_to_csr, only: csr_matrix_t
    implicit none
    call check_leaf()
    call check_reduced_derivatives()
    print '(a)', 'PASS ideal map wall slide'
contains
    subroutine check_leaf()
        real(dp) :: a(9), d(9), o(6), op(6), om(6), dout(6), bar(6), gin(9)
        real(dp), parameter :: h = 1.0e-5_dp
        real(dp) :: err
        a = [0.37_dp, 0.21_dp, 3.0_dp, 1.0_dp, 2.0_dp, 0.3_dp, -0.2_dp, 0.15_dp, &
            0.4_dp]
        d = [0.3_dp, -0.7_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.5_dp, 0.2_dp, -0.4_dp, 0.9_dp]
        call leaf(a, o)
        ! theta-derivatives against central differences.
        call leaf(a + h*unit(1), op); call leaf(a - h*unit(1), om)
        err = maxval(abs([(op(1) - om(1))/(2*h) - o(2), (op(2) - om(2))/(2*h) - o(3), &
            (op(4) - om(4))/(2*h) - o(5), (op(5) - om(5))/(2*h) - o(6)]))
        print '(a,es10.2)', 'wall leaf theta-derivative FD error:', err
        if (err > 1.0e-8_dp) error stop 'wall leaf theta derivatives'
        ! JVP against a directional central difference.
        call leaf(a + h*d, op); call leaf(a - h*d, om)
        call seed_wall_mode_jvp(a(1), a(2), a(3), a(4), a(5), a(6), a(7), a(8), &
            a(9), d(1), d(2), d(3), d(4), d(5), d(6), d(7), d(8), d(9), dout(1), &
            dout(2), dout(3), dout(4), dout(5), dout(6))
        err = maxval(abs((op - om)/(2*h) - dout))
        print '(a,es10.2)', 'wall leaf JVP FD error:', err
        if (err > 1.0e-8_dp) error stop 'wall leaf JVP'
        ! VJP pairing <bar, J d> = <J^T bar, d>.
        bar = [0.7_dp, -0.3_dp, 0.2_dp, 0.5_dp, -0.8_dp, 0.1_dp]
        call seed_wall_mode_vjp(a(1), a(2), a(3), a(4), a(5), a(6), a(7), a(8), &
            a(9), bar(1), bar(2), bar(3), bar(4), bar(5), bar(6), gin(1), gin(2), &
            gin(3), gin(4), gin(5), gin(6), gin(7), gin(8), gin(9))
        err = abs(dot_product(bar, dout) - dot_product(gin, d))
        print '(a,es10.2)', 'wall leaf VJP pairing error:', err
        if (err > 1.0e-13_dp) error stop 'wall leaf VJP'
    end subroutine

    subroutine leaf(a, o)
        real(dp), intent(in) :: a(9)
        real(dp), intent(out) :: o(6)
        call seed_wall_mode(a(1), a(2), a(3), a(4), a(5), a(6), a(7), a(8), a(9), &
            o(1), o(2), o(3), o(4), o(5), o(6))
    end subroutine

    pure function unit(k) result(e)
        integer, intent(in) :: k
        real(dp) :: e(9)
        e = 0; e(k) = 1
    end function

    subroutine check_reduced_derivatives()
        type(gs_state_t) :: mesh
        type(ideal_map_t) :: map
        type(toroidal_seed_t) :: seed
        type(toroidal_profile_t) :: profile
        type(csr_matrix_t) :: matrix
        real(dp) :: c(4, 3, 0:2), e0, ep, em, slope, curvature, hdd, eps
        real(dp), allocatable :: g(:), gr(:), d(:), hd(:), g0(:, :), x0(:, :, :)
        integer :: status, k, i, p, free
        ! Shaped stellarator-symmetric boundary, nfp 2, ntor 1, in the
        ! convention of kin6d_ideal_seed_3d (rc, rs, zc, zs) per (n, m).
        c = 0
        c(1, 2, 0) = 3.0_dp
        c(1, 2, 1) = 0.8_dp; c(4, 2, 1) = 1.1_dp
        c(1, 2, 2) = 0.12_dp
        c(1, 3, 1) = 0.06_dp; c(4, 3, 1) = 0.05_dp
        profile%iota = [0.55_dp, 0.25_dp]
        profile%pressure = [0.02_dp, -0.018_dp]
        call ideal_seed_disk_mesh(mesh, 8, 2, status)
        if (status /= 0) error stop 'disk mesh'
        call ideal_seed_3d_init(map, seed, mesh, c, 2, 1.0_dp, profile, status, &
            sliding_wall=.true.)
        if (status /= IDEAL_MAP_OK) error stop 'sliding seed'
        if (map%nwall /= count(map%boundary)) error stop 'wall nodes'
        free = count(.not. map%fixed)
        ! Off the seed: smooth wall angles in every mode (the pinned mean of
        ! one node kept) and an interior perturbation.
        allocate (g0(map%nmode, map%nwall))
        do k = 1, map%nwall
            do i = 1, map%nmode
                g0(i, k) = 0.04_dp*sin(i + 1.3_dp*map%wall_angle(k)) &
                    + 0.02_dp*cos(2*map%wall_angle(k) - i)
            end do
        end do
        g0(1, map%wall_pin) = 0
        call ideal_map_wall_set(map, g0)
        allocate (g(map%ndof), gr(map%nred), d(map%nred), hd(map%nred))
        do k = 1, map%nred
            d(k) = sin(0.7_dp*k) + 0.3_dp*cos(1.9_dp*k)
        end do
        where (map%fixed) d = 0
        ! Interior displacement scale ~ 1e-2, wall angles ~ 1e-2 rad.
        d = 1.0e-2_dp*d
        call ideal_map_expand(map, d, g)
        x0 = map%x
        call ideal_map_energy(map, map%x, e0, status, g)
        if (status /= IDEAL_MAP_OK) error stop 'energy'
        call ideal_map_reduce(map, g, gr)
        where (map%fixed) gr = 0
        call ideal_map_hessian(map, map%x, matrix)
        hd = 0
        do i = 1, map%nred
            if (map%fixed(i)) cycle
            do p = matrix%rowptr(i), matrix%rowptr(i + 1) - 1
                if (map%fixed(matrix%colidx(p))) cycle
                hd(i) = hd(i) + matrix%val(p)*d(matrix%colidx(p))
            end do
        end do
        eps = 1.0e-3_dp
        call along(map, x0, g0, d, eps, ep)
        call along(map, x0, g0, d, -eps, em)
        slope = (ep - em)/(2*eps)
        curvature = (ep - 2*e0 + em)/eps**2
        hdd = dot_product(d, hd)
        print '(a,i0,a,i0)', 'sliding wall free dofs ', free, ' wall nodes ', map%nwall
        print '(a,2es14.6,es10.2)', 'directional slope reduced/FD:', &
            dot_product(gr, d), slope, abs(slope/dot_product(gr, d) - 1)
        print '(a,2es14.6,es10.2)', 'directional curvature Hessian/FD:', hdd, &
            curvature, abs(curvature/hdd - 1)
        if (abs(slope/dot_product(gr, d) - 1) > 1.0e-5_dp) error stop 'reduced gradient'
        if (abs(curvature/hdd - 1) > 1.0e-4_dp) error stop 'reduced Hessian'
        ! Wall-only direction: isolates the wall Jacobian and curvature term.
        do i = 1, map%nred
            if (mod(i, 2) == 0) d(i) = 0
        end do
        do k = 1, map%nodes
            if (map%boundary(k)) cycle
            do i = 1, map%nmode
                d(ideal_map_rdof(map, 1, i, k)) = 0
            end do
        end do
        hd = 0
        do i = 1, map%nred
            if (map%fixed(i)) cycle
            do p = matrix%rowptr(i), matrix%rowptr(i + 1) - 1
                if (map%fixed(matrix%colidx(p))) cycle
                hd(i) = hd(i) + matrix%val(p)*d(matrix%colidx(p))
            end do
        end do
        call along(map, x0, g0, d, eps, ep)
        call along(map, x0, g0, d, -eps, em)
        slope = (ep - em)/(2*eps)
        curvature = (ep - 2*e0 + em)/eps**2
        hdd = dot_product(d, hd)
        print '(a,2es14.6,es10.2)', 'wall-only slope reduced/FD:', &
            dot_product(gr, d), slope, abs(slope/dot_product(gr, d) - 1)
        print '(a,2es14.6,es10.2)', 'wall-only curvature Hessian/FD:', hdd, &
            curvature, abs(curvature/hdd - 1)
        if (abs(slope/dot_product(gr, d) - 1) > 1.0e-5_dp) error stop 'wall gradient'
        if (abs(curvature/hdd - 1) > 1.0e-4_dp) error stop 'wall Hessian'
    end subroutine

    subroutine along(map, x0, g0, d, t, energy)
        ! Retraction of the solve: x + t expand(d), wall placed at g + t d.
        type(ideal_map_t), intent(in) :: map
        real(dp), intent(in) :: x0(:, :, :), g0(:, :), d(:), t
        real(dp), intent(out) :: energy
        type(ideal_map_t) :: work
        real(dp), allocatable :: v(:), gw(:, :)
        integer :: kk, ii, st
        work = map
        allocate (v(map%ndof), gw(map%nmode, map%nwall))
        call ideal_map_expand(map, t*d, v)
        work%x = x0 + reshape(v, shape(x0))
        do kk = 1, map%nwall
            do ii = 1, map%nmode
                gw(ii, kk) = g0(ii, kk) &
                    + t*d(ideal_map_rdof(map, 1, ii, map%wall_node(kk)))
            end do
        end do
        call ideal_map_wall_set(work, gw)
        call ideal_map_energy(work, work%x, energy, st)
        if (st /= IDEAL_MAP_OK) error stop 'energy along the retraction'
    end subroutine
end program
