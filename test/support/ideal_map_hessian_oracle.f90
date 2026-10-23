module ideal_map_hessian_oracle
    ! Frozen copy of the coefficient-Hessian assembly of kin6d_ideal_map before
    ! its sum factorisation (point loop with per-point rank updates, then
    ! P^T H P and reflection). It is the independent oracle of the production
    ! assembly in test_ideal_map_hessian: same definition, different loop
    ! structure. Do not optimise or share code with src/mhd/ideal_map.f90.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_ideal_map, only: ideal_map_t
    use kin6d_ideal_map_leaves, only: ideal_energy_density, ideal_energy_gradient, &
        ideal_energy_hessian
    use fortfem_element_to_csr, only: csr_matrix_t, csr_build_pattern, &
        csr_assemble_blocks
    implicit none
    private
    public :: ideal_map_hessian_reference, point_jet
contains
    subroutine point_jet(map, e, g, j, x, q, qd)
        ! q(3) and Q(k, i) = dq_k/ds_i at quadrature point g, plane j.
        type(ideal_map_t), intent(in) :: map
        real(dp), intent(in) :: x(:, :, :)
        integer, intent(in) :: e, g, j
        real(dp), intent(out) :: q(3), qd(3, 3)
        real(dp) :: v
        integer :: a, m, node
        q = 0; qd = 0
        do a = 1, map%nloc
            node = map%conn(a, e)
            do m = 1, map%nmode
                v = map%trig(m, j)
                q = q + x(:, m, node)*map%shape(a, g)*v
                qd(:, 1) = qd(:, 1) + x(:, m, node)*map%dshape(1, a, g)*v
                qd(:, 2) = qd(:, 2) + x(:, m, node)*map%dshape(2, a, g)*v
                qd(:, 3) = qd(:, 3) + x(:, m, node)*map%shape(a, g)*map%dtrig(m, j)
            end do
        end do
        q(3) = q(3) + map%phi(j)
        qd(3, 3) = qd(3, 3) + 1
    end subroutine

    subroutine local_derivatives(map, e, g, j, x, w, jac, gv, hv, want_hessian)
        ! Energy density w, Jacobian and its first/second derivatives with
        ! respect to v = (q(1:3), Q(1:3, 1:3)) (index 3 + k + 3*(i-1) for Q(k,i)).
        type(ideal_map_t), intent(in) :: map
        integer, intent(in) :: e, g, j
        real(dp), intent(in) :: x(:, :, :)
        real(dp), intent(out) :: w, jac, gv(12), hv(12, 12)
        logical, intent(in) :: want_hessian
        real(dp) :: q(3), qd(3, 3), pos(3), d1(3, 3), d2(3, 3, 3), d3(3, 3, 3, 3)
        real(dp) :: dx(3, 3), gd(3, 3), hd(9, 9), dvx(9, 12), hflat(45)
        integer :: a, i, k, m, n, r, c
        call point_jet(map, e, g, j, x, q, qd)
        call map%chart%jet(q, pos, d1, d2, d3)
        dx = matmul(d1, qd)
        call ideal_energy_density(dx(1, 1), dx(2, 1), dx(3, 1), dx(1, 2), &
            dx(2, 2), dx(3, 2), dx(1, 3), dx(2, 3), dx(3, 3), map%b(1, g, e), &
            map%b(2, g, e), map%b(3, g, e), map%p(g, e), map%mu0, &
            map%sigma, w, jac)
        call ideal_energy_gradient(dx(1, 1), dx(2, 1), dx(3, 1), dx(1, 2), &
            dx(2, 2), dx(3, 2), dx(1, 3), dx(2, 3), dx(3, 3), map%b(1, g, e), &
            map%b(2, g, e), map%b(3, g, e), map%p(g, e), map%mu0, &
            map%sigma, gd(1, 1), gd(2, 1), gd(3, 1), gd(1, 2), &
            gd(2, 2), gd(3, 2), gd(1, 3), gd(2, 3), gd(3, 3))
        ! dDX(a,i)/dv: q_m -> sum_k d2(a,k,m) Q(k,i); Q(k,i') -> d1(a,k) delta.
        dvx = 0
        do i = 1, 3
            do a = 1, 3
                r = a + 3*(i - 1)
                do m = 1, 3
                    dvx(r, m) = dot_product(d2(a, :, m), qd(:, i))
                end do
                do k = 1, 3
                    dvx(r, 3 + k + 3*(i - 1)) = d1(a, k)
                end do
            end do
        end do
        gv = matmul(reshape(gd, [9]), dvx)
        hv = 0
        if (.not. want_hessian) return
        call hessian_leaf(hflat)
        n = 0
        do c = 1, 9
            do r = 1, c
                n = n + 1
                hd(r, c) = hflat(n); hd(c, r) = hflat(n)
            end do
        end do
        hv = matmul(transpose(dvx), matmul(hd, dvx))
        do m = 1, 3
            do n = 1, 3
                do i = 1, 3
                    do a = 1, 3
                        hv(m, n) = hv(m, n) + gd(a, i)*dot_product(d3(a, :, m, n), &
                            qd(:, i))
                    end do
                end do
            end do
            do i = 1, 3
                do k = 1, 3
                    r = 3 + k + 3*(i - 1)
                    hv(m, r) = hv(m, r) + dot_product(gd(:, i), d2(:, k, m))
                    hv(r, m) = hv(m, r)
                end do
            end do
        end do
    contains
        subroutine hessian_leaf(t)
            real(dp), intent(out) :: t(45)
            call ideal_energy_hessian(dx(1, 1), dx(2, 1), dx(3, 1), dx(1, 2), &
                dx(2, 2), dx(3, 2), dx(1, 3), dx(2, 3), dx(3, 3), map%b(1, g, e), &
                map%b(2, g, e), map%b(3, g, e), map%p(g, e), map%mu0, &
                map%sigma, t(1), t(2), t(3), t(4), t(5), t(6), t(7), &
                t(8), t(9), t(10), t(11), t(12), t(13), t(14), t(15), t(16), t(17), &
                t(18), t(19), t(20), t(21), t(22), t(23), t(24), t(25), t(26), &
                t(27), t(28), t(29), t(30), t(31), t(32), t(33), t(34), t(35), &
                t(36), t(37), t(38), t(39), t(40), t(41), t(42), t(43), t(44), t(45))
        end subroutine
    end subroutine

    pure subroutine test_functions(map, g, j, a, m, f)
        ! f(0:3) = (N T, dN/ds1 T, dN/ds2 T, N dT/ds3) of node a, mode m.
        type(ideal_map_t), intent(in) :: map
        integer, intent(in) :: g, j, a, m
        real(dp), intent(out) :: f(0:3)
        f(0) = map%shape(a, g)*map%trig(m, j)
        f(1) = map%dshape(1, a, g)*map%trig(m, j)
        f(2) = map%dshape(2, a, g)*map%trig(m, j)
        f(3) = map%shape(a, g)*map%dtrig(m, j)
    end subroutine

    subroutine ideal_map_hessian_reference(map, x, matrix)
        ! Assembled coefficient Hessian (all Fourier modes coupled) as CSR. The
        ! sparsity pattern of a matrix already built for this map is reused.
        type(ideal_map_t), intent(in) :: map
        real(dp), intent(in) :: x(:, :, :)
        type(csr_matrix_t), intent(inout) :: matrix
        real(dp), allocatable :: blocks(:, :, :), fa(:, :), fb(:, :), full(:, :)
        real(dp), allocatable :: half(:, :)
        real(dp) :: w, jac, gv(12), hv(12, 12), dphi, weight, s
        integer :: e, g, j, a, m, c, r, l, lb, nl, c2, r2, nr
        nl = 3*map%nmode*map%nloc
        nr = 2*map%nmode*map%nloc + 2*map%npoly*map%nmode
        allocate (blocks(nr, nr, map%cells))
        dphi = 2*acos(-1.0_dp)/(map%nfp*map%nphi)
        !$omp parallel do default(none) shared(map, x, blocks, dphi, nl, nr) &
        !$omp private(e, g, j, a, m, c, r, l, lb, w, jac, gv, hv, weight, fa, fb, &
        !$omp s, c2, r2, full, half) schedule(static)
        do e = 1, map%cells
            allocate (fa(0:3, map%nloc*map%nmode), fb(nl, 12), full(nl, nl), &
                half(nl, nr))
            full = 0
            do j = 1, map%nphi
                do g = 1, map%nq
                    call local_derivatives(map, e, g, j, x, w, jac, gv, hv, .true.)
                    weight = map%weight(g)*dphi
                    l = 0
                    do a = 1, map%nloc
                        do m = 1, map%nmode
                            l = l + 1
                            call test_functions(map, g, j, a, m, fa(:, l))
                        end do
                    end do
                    ! fb(row, v) = sum_r f_r(row) hv(idx(c,r), v)
                    do l = 1, map%nloc*map%nmode
                        do c = 1, 3
                            lb = 3*(l - 1) + c
                            fb(lb, :) = 0
                            do r = 0, 3
                                fb(lb, :) = fb(lb, :) + fa(r, l)*hv(c + 3*r, :)
                            end do
                        end do
                    end do
                    do l = 1, map%nloc*map%nmode
                        do c2 = 1, 3
                            lb = 3*(l - 1) + c2
                            do r2 = 0, 3
                                s = weight*fa(r2, l)
                                full(:, lb) = full(:, lb) + s*fb(:, c2 + 3*r2)
                            end do
                        end do
                    end do
                end do
            end do
            do a = 1, map%nloc
                do c = 1, map%nloc
                    l = 3*(1 + map%nmode*(a - 1))
                    r = 3*(1 + map%nmode*(c - 1))
                    full(l, r) = full(l, r) + map%omega_gauge(a, c, e)
                end do
            end do
            ! Reduce to (nu, omega): P maps (nu, omega) of local node-mode l to
            ! (d1 nu, d2 nu, omega).
            do l = 1, map%nloc*map%nmode
                a = map%conn(1 + (l - 1)/map%nmode, e)
                half(:, 2*l - 1) = matmul(full(:, 3*l - 2:3*l), map%basis(:, 1, a))
                half(:, 2*l) = matmul(full(:, 3*l - 2:3*l), map%basis(:, 2, a))
            end do
            do m = 1, map%nmode
                do j = 1, map%npoly
                    do c = 1, 2
                        r = 2*map%nmode*map%nloc + c + 2*(j - 1) + 2*map%npoly*(m - 1)
                        half(:, r) = 0
                        do a = 1, map%nloc
                            l = m + map%nmode*(a - 1)
                            if (allocated(map%axis_vec)) then
                                half(:, r) = half(:, r) &
                                    + matmul(full(:, 3*l - 2:3*l - 1), &
                                    map%axis_vec(:, c + 2*(j - 1), map%conn(a, e)))
                                cycle
                            end if
                            half(:, r) = half(:, r) &
                                + map%axis_basis(j, map%conn(a, e))*full(:, 3*l - 3 + c)
                        end do
                    end do
                end do
            end do
            do l = 1, map%nloc*map%nmode
                a = map%conn(1 + (l - 1)/map%nmode, e)
                blocks(2*l - 1, :, e) = matmul(map%basis(:, 1, a), half(3*l - 2:3*l, :))
                blocks(2*l, :, e) = matmul(map%basis(:, 2, a), half(3*l - 2:3*l, :))
            end do
            do m = 1, map%nmode
                do j = 1, map%npoly
                    do c = 1, 2
                        r = 2*map%nmode*map%nloc + c + 2*(j - 1) + 2*map%npoly*(m - 1)
                        blocks(r, :, e) = 0
                        do a = 1, map%nloc
                            l = m + map%nmode*(a - 1)
                            if (allocated(map%axis_vec)) then
                                blocks(r, :, e) = blocks(r, :, e) + matmul( &
                                    map%axis_vec(:, c + 2*(j - 1), map%conn(a, e)), &
                                    half(3*l - 2:3*l - 1, :))
                                cycle
                            end if
                            blocks(r, :, e) = blocks(r, :, e) &
                                + map%axis_basis(j, map%conn(a, e))*half(3*l - 3 + c, :)
                        end do
                    end do
                end do
            end do
            ! The scalar-energy Hessian has one authoritative triangle.
            ! Reflect after P^T H P: cancellation in independently evaluated
            ! contractions must not create an antisymmetric Newton operator.
            do c = 1, nr
                do r = 1, c - 1
                    blocks(r, c, e) = blocks(c, r, e)
                end do
            end do
            deallocate (fa, fb, full, half)
        end do
        !$omp end parallel do
        if (matrix%n /= map%nred .or. matrix%nloc /= nr) &
            call csr_build_pattern(map%nred, map%redconn, .false., matrix)
        call csr_assemble_blocks(matrix, map%redconn, blocks)
    end subroutine ideal_map_hessian_reference
end module ideal_map_hessian_oracle
