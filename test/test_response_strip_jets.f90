program test_response_strip_jets
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_vacuum_surfaces, only: psi_grid_t, eqdsk_read
    use kin6d_response_strip_jets, only: strip_equilibrium_jets_t, &
        strip_equilibrium_point, strip_equilibrium_grid
    implicit none
    type(psi_grid_t) :: g
    type(strip_equilibrium_jets_t) :: jet, plus, minus, edge
    type(strip_equilibrium_jets_t) :: grid(5, 9)
    real(dp), parameter :: r0 = 3.0_dp, kappa = 1.4_dp, minor = 0.8_dp
    real(dp), parameter :: p0 = 7.0_dp, p1 = -3.0_dp, p2 = 1.2_dp
    real(dp) :: theta, u, h, ray, c, sn, rr, zz, x, depth, pr, pz, pp
    real(dp) :: gp(3), hp(3, 3), identity(2, 2), err, radial(5)
    integer :: i, j, status
    character(512) :: file
    call get_command_argument(1, file)
    if (len_trim(file) > 0) then
        call eqdsk_read(trim(file), g, status)
        call check(status == 0, 'read TC24')
        radial = sqrt([0.15_dp, 0.4_dp, 0.65_dp, 0.85_dp, 0.999_dp])
        call strip_equilibrium_grid(g, radial, 9, 1.0e-6_dp, &
            1.0e-7_dp, grid, status)
        call check(status == 0, 'TC24 grid jets')
        print '(a,es12.4)', 'TC24 max label defect=', maxval(grid%label_defect)
        print '(a,es12.4)', 'TC24 max Ampere defect [A/m^2]=', &
            maxval(grid%ampere_defect)
        print '(a,2es12.4)', 'TC24 retained outer pressure [Pa]=', &
            minval(grid(5, :)%fluid%p), maxval(grid(5, :)%fluid%p)
        call check(maxval(grid%label_defect) < 2.0e-12_dp, 'TC24 level equation')
        call check(minval(grid(5, :)%fluid%p) > 0, 'TC24 edge pressure retained')
        h = 2.0e-6_dp
        err = 0
        do j = 1, 9
            theta = 2*acos(-1.0_dp)*(j - 1)/9
            do i = 1, 5
                u = radial(i)
                jet = grid(i, j)
                call strip_equilibrium_point(g, u + h, theta, 1.0e-6_dp, &
                    1.0e-7_dp, plus, status)
                call check(status == 0, 'TC24 plus')
                call strip_equilibrium_point(g, u - h, theta, 1.0e-6_dp, &
                    1.0e-7_dp, minus, status)
                call check(status == 0, 'TC24 minus')
                gp = matmul(jet%db(:, [1, 3]), jet%jac(:, 1))
                err = max(err, norm2((plus%fluid%b0 - minus%fluid%b0)/(2*h) - gp) &
                    /max(1.0_dp, norm2(gp)))
                gp = matmul(jet%fluid%hp(:, [1, 3]), jet%jac(:, 1))
                call check(norm2((plus%fluid%gp - minus%fluid%gp)/(2*h) - gp) &
                    /max(1.0_dp, norm2(gp)) < 2.0e-5_dp, 'TC24 pressure Hessian')
                call check(maxval(abs((plus%x - minus%x)/(2*h) - jet%jac(:, 1))) &
                    < 2.0e-5_dp, 'TC24 implicit map tangent')
                call check(maxval(abs((plus%jac - minus%jac)/(2*h) &
                    - jet%hmap(:, 1:2)))/max(1.0_dp, maxval(abs(jet%hmap(:, 1:2)))) &
                    < 2.0e-5_dp, 'TC24 implicit map Hessian radial column')
                pp = dot_product(jet%fluid%gp([1, 3]), jet%jac(:, 1))
                call check(abs((plus%fluid%p - minus%fluid%p)/(2*h) - pp) &
                    /max(1.0_dp, abs(pp)) < 2.0e-5_dp, 'TC24 pressure tangent')
                call strip_equilibrium_point(g, u, theta + h, 1.0e-6_dp, &
                    1.0e-7_dp, plus, status)
                call check(status == 0, 'TC24 theta plus')
                call strip_equilibrium_point(g, u, theta - h, 1.0e-6_dp, &
                    1.0e-7_dp, minus, status)
                call check(status == 0, 'TC24 theta minus')
                gp = matmul(jet%db(:, [1, 3]), jet%jac(:, 2))
                err = max(err, norm2((plus%fluid%b0 - minus%fluid%b0)/(2*h) - gp) &
                    /max(1.0_dp, norm2(gp)))
                call check(maxval(abs((plus%x - minus%x)/(2*h) - jet%jac(:, 2))) &
                    < 2.0e-5_dp, 'TC24 implicit map angular tangent')
                call check(maxval(abs((plus%jac - minus%jac)/(2*h) &
                    - jet%hmap(:, 2:3)))/max(1.0_dp, maxval(abs(jet%hmap(:, 2:3)))) &
                    < 2.0e-5_dp, 'TC24 implicit map Hessian angular column')
                gp = matmul(jet%fluid%hp(:, [1, 3]), jet%jac(:, 2))
                call check(norm2((plus%fluid%gp - minus%fluid%gp)/(2*h) - gp) &
                    /max(1.0_dp, norm2(gp)) < 2.0e-5_dp, &
                    'TC24 angular pressure Hessian')
            end do
        end do
        print '(a,es12.4)', 'TC24 B-jet difference error=', err
        call check(err < 2.0e-5_dp, 'TC24 B jets independent differences')
        stop
    end if
    call manufactured_grid(g)
    depth = g%psi_edge - g%psi_axis
    identity = 0; identity(1, 1) = 1; identity(2, 2) = 1
    u = 0.73_dp; h = 2.0e-5_dp
    do j = 1, 17
        theta = 2*acos(-1.0_dp)*(j - 1)/17
        c = cos(theta); sn = sin(theta)
        ray = minor*u/sqrt(c*c + sn*sn/kappa**2)
        call strip_equilibrium_point(g, u, theta, 0.4_dp, 1.0e-7_dp, jet, status)
        call check(status == 0, 'manufactured jets admitted')
        err = maxval(abs(jet%x - [r0 + ray*c, ray*sn]))
        call check(err < 2.0e-12_dp, 'elliptical level map')
        call check(maxval(abs(matmul(jet%inverse, jet%jac) - identity)) &
            < 1.0e-12_dp, 'inverse map composition')
        rr = jet%x(1); zz = jet%x(2); x = u*u
        pr = rr - r0; pz = zz/kappa**2
        gp = [(p1 + 2*p2*x)*pr/depth, 0.0_dp, (p1 + 2*p2*x)*pz/depth]
        hp = 0
        hp(1, 1) = 2*p2*pr**2/depth**2 + (p1 + 2*p2*x)/depth
        hp(1, 3) = 2*p2*pr*pz/depth**2; hp(3, 1) = hp(1, 3)
        hp(3, 3) = 2*p2*pz**2/depth**2 + (p1 + 2*p2*x)/(depth*kappa**2)
        hp(2, 2) = gp(1)/rr
        call check(abs(jet%fluid%p - (p0 + p1*x + p2*x*x)) < 1.0e-12_dp, &
            'quadratic retained pressure')
        call check(maxval(abs(jet%fluid%gp - gp)) < 2.0e-12_dp, 'pressure gradient')
        call check(maxval(abs(jet%fluid%hp - hp)) < 2.0e-11_dp, &
            'Cartesian pressure Hessian including toroidal basis')
        call check(abs(jet%db(1, 1) + jet%db(2, 2) + jet%db(3, 3)) &
            < 2.0e-13_dp, 'div B=0')
        call check(maxval(abs(jet%db(:, 2) &
            - [-jet%fluid%b0(2), jet%fluid%b0(1), 0.0_dp]/rr)) &
            < 2.0e-13_dp, 'Cartesian B basis derivative')
        call strip_equilibrium_point(g, u + h, theta, 0.4_dp, &
            1.0e-7_dp, plus, status)
        call check(status == 0, 'plus jet')
        call strip_equilibrium_point(g, u - h, theta, 0.4_dp, &
            1.0e-7_dp, minus, status)
        call check(status == 0, 'minus jet')
        call check(maxval(abs((plus%x - minus%x)/(2*h) - jet%jac(:, 1))) &
            < 2.0e-8_dp, 'implicit map tangent')
        call check(maxval(abs((plus%jac - minus%jac)/(2*h) &
            - jet%hmap(:, 1:2))) < 2.0e-8_dp, 'map Hessian u column')
        call check(maxval(abs((plus%fluid%b0 - minus%fluid%b0)/(2*h) &
            - matmul(jet%db(:, [1, 3]), jet%jac(:, 1)))) &
            < 2.0e-8_dp, 'B spatial jet chain rule')
        call strip_equilibrium_point(g, u, theta + h, 0.4_dp, &
            1.0e-7_dp, plus, status)
        call check(status == 0, 'theta plus')
        call strip_equilibrium_point(g, u, theta - h, 0.4_dp, &
            1.0e-7_dp, minus, status)
        call check(status == 0, 'theta minus')
        call check(maxval(abs((plus%jac - minus%jac)/(2*h) &
            - jet%hmap(:, 2:3))) < 2.0e-8_dp, 'map Hessian theta column')
    end do
    call strip_equilibrium_point(g, 1.0_dp, 0.0_dp, 0.4_dp, 0.0_dp, edge, status)
    call check(status == 0, 'edge jets')
    call check(abs(edge%fluid%p - (p0 + p1 + p2)) < 1.0e-12_dp, &
        'nonzero edge pressure retained')
    call strip_equilibrium_point(g, 0.0_dp, 0.0_dp, 0.4_dp, 0.0_dp, edge, status)
    call check(status /= 0, 'singular axis chart rejected')
    call products()
    call map_products()
    print *, 'PASS manufactured map/B/pressure jets and products'
contains
    subroutine check(ok, label)
        logical, intent(in) :: ok
        character(*), intent(in) :: label
        if (.not. ok) then
            print *, 'FAIL ', label
            error stop 1
        end if
    end subroutine
    subroutine manufactured_grid(g)
        type(psi_grid_t), intent(out) :: g
        real(dp) :: rr, zz, t, spacing
        integer :: i, j
        g%nr = 25; g%nz = 29; g%nf = 21
        g%rleft = r0 - 1.3_dp; g%dr = 2.6_dp/(g%nr - 1)
        g%zbot = -1.8_dp; g%dz = 3.6_dp/(g%nz - 1)
        g%r_axis = r0; g%z_axis = 0; g%psi_axis = 0
        g%psi_edge = minor**2/2; g%ip = 1; g%s_pol = 1; g%sgn_dpsi = 1
        allocate (g%coef(0:g%nr + 1, 0:g%nz + 1), &
            g%fcoef(0:g%nf + 1), g%pcoef(0:g%nf + 1))
        do j = 0, g%nz + 1
            zz = g%zbot + (j - 1)*g%dz
            do i = 0, g%nr + 1
                rr = g%rleft + (i - 1)*g%dr - r0
                g%coef(i, j) = (rr**2 + zz**2/kappa**2 &
                    - (g%dr**2 + g%dz**2/kappa**2)/3)/2
            end do
        end do
        spacing = 1.0_dp/(g%nf - 1)
        do i = 0, g%nf + 1
            t = (i - 1)*spacing
            g%fcoef(i) = 4.0_dp + 0.3_dp*t
            g%pcoef(i) = p0 + p1*t + p2*(t*t - spacing**2/3)
        end do
    end subroutine
    subroutine products()
        use kin6d_response_strip_jets_leaves, only: strip_field_jets, &
            strip_field_jets_jvp, strip_field_jets_vjp
        real(dp) :: v(13), d(13), seed(24), y(24), yp(24), ym(24), dy(24), adj(13)
        real(dp) :: step, pairing
        integer :: k
        v = [3.2_dp, 0.5_dp, 1.0_dp, 0.1_dp, 0.3_dp, 1.0_dp, &
            0.2_dp, 2.0_dp, 4.0_dp, 0.4_dp, -3.0_dp, 0.8_dp, 0.9_dp]
        do k = 1, 13
            d(k) = sin(real(k, dp))
        end do
        do k = 1, 24
            seed(k) = cos(real(k, dp))
        end do
        step = 1.0e-6_dp
        call strip_field_jets(v(1), v(2), v(3), v(4), v(5), v(6), v(7), v(8), &
            v(9), v(10), v(11), v(12), v(13), y(1), y(2), y(3), y(4), y(5), &
            y(6), y(7), y(8), y(9), y(10), y(11), y(12), y(13), y(14), y(15), &
            y(16), y(17), y(18), y(19), y(20), y(21), y(22), y(23), y(24))
        call strip_field_jets_jvp(v(1), v(2), v(3), v(4), v(5), v(6), v(7), &
            v(8), v(9), v(10), v(11), v(12), v(13), d(1), d(2), d(3), d(4), &
            d(5), d(6), d(7), d(8), d(9), d(10), d(11), d(12), d(13), dy(1), &
            dy(2), dy(3), dy(4), dy(5), dy(6), dy(7), dy(8), dy(9), dy(10), &
            dy(11), dy(12), dy(13), dy(14), dy(15), dy(16), dy(17), dy(18), &
            dy(19), dy(20), dy(21), dy(22), dy(23), dy(24))
        call strip_field_jets_vjp(v(1), v(2), v(3), v(4), v(5), v(6), v(7), &
            v(8), v(9), v(10), v(11), v(12), v(13), seed(1), seed(2), seed(3), &
            seed(4), seed(5), seed(6), seed(7), seed(8), seed(9), seed(10), &
            seed(11), seed(12), seed(13), seed(14), seed(15), seed(16), &
            seed(17), seed(18), seed(19), seed(20), seed(21), seed(22), &
            seed(23), seed(24), adj(1), adj(2), adj(3), adj(4), adj(5), adj(6), &
            adj(7), adj(8), adj(9), adj(10), adj(11), adj(12), adj(13))
        v = v + step*d
        call strip_field_jets(v(1), v(2), v(3), v(4), v(5), v(6), v(7), v(8), &
            v(9), v(10), v(11), v(12), v(13), yp(1), yp(2), yp(3), yp(4), &
            yp(5), yp(6), yp(7), yp(8), yp(9), yp(10), yp(11), yp(12), yp(13), &
            yp(14), yp(15), yp(16), yp(17), yp(18), yp(19), yp(20), yp(21), &
            yp(22), yp(23), yp(24))
        v = v - 2*step*d
        call strip_field_jets(v(1), v(2), v(3), v(4), v(5), v(6), v(7), v(8), &
            v(9), v(10), v(11), v(12), v(13), ym(1), ym(2), ym(3), ym(4), &
            ym(5), ym(6), ym(7), ym(8), ym(9), ym(10), ym(11), ym(12), ym(13), &
            ym(14), ym(15), ym(16), ym(17), ym(18), ym(19), ym(20), ym(21), &
            ym(22), ym(23), ym(24))
        call check(maxval(abs((yp - ym)/(2*step) - dy)) < 2.0e-8_dp, &
            'generated field jet JVP independent difference')
        pairing = abs(dot_product(seed, dy) - dot_product(adj, d))
        call check(pairing < 2.0e-12_dp, 'generated field jet VJP pairing')
    end subroutine
    subroutine map_products()
        use kin6d_response_strip_jets_leaves, only: strip_map_jets, &
            strip_map_jets_jvp, strip_map_jets_vjp
        real(dp) :: v(10), d(10), seed(14), yp(14), ym(14), dy(14), adj(10)
        real(dp) :: step
        integer :: k
        v = [0.7_dp, 0.5_dp, 0.3_dp, 0.7_dp, 0.6_dp, 0.1_dp, &
            0.3_dp, 1.0_dp, 0.2_dp, 2.0_dp]
        do k = 1, 10
            d(k) = 0.1_dp*sin(real(k, dp))
        end do
        do k = 1, 14
            seed(k) = cos(real(k, dp))
        end do
        step = 1.0e-6_dp
        call strip_map_jets_jvp(v(1), v(2), v(3), v(4), v(5), v(6), v(7), v(8), &
            v(9), v(10), d(1), d(2), d(3), d(4), d(5), d(6), d(7), d(8), d(9), &
            d(10), dy(1), dy(2), dy(3), dy(4), dy(5), dy(6), dy(7), dy(8), &
            dy(9), dy(10), dy(11), dy(12), dy(13), dy(14))
        call strip_map_jets_vjp(v(1), v(2), v(3), v(4), v(5), v(6), v(7), v(8), &
            v(9), v(10), seed(1), seed(2), seed(3), seed(4), seed(5), seed(6), &
            seed(7), seed(8), seed(9), seed(10), seed(11), seed(12), seed(13), &
            seed(14), adj(1), adj(2), adj(3), adj(4), adj(5), adj(6), adj(7), &
            adj(8), adj(9), adj(10))
        v = v + step*d
        call strip_map_jets(v(1), v(2), v(3), v(4), v(5), v(6), v(7), v(8), &
            v(9), v(10), yp(1), yp(2), yp(3), yp(4), yp(5), yp(6), yp(7), &
            yp(8), yp(9), yp(10), yp(11), yp(12), yp(13), yp(14))
        v = v - 2*step*d
        call strip_map_jets(v(1), v(2), v(3), v(4), v(5), v(6), v(7), v(8), &
            v(9), v(10), ym(1), ym(2), ym(3), ym(4), ym(5), ym(6), ym(7), &
            ym(8), ym(9), ym(10), ym(11), ym(12), ym(13), ym(14))
        call check(maxval(abs((yp - ym)/(2*step) - dy)) < 2.0e-8_dp, &
            'generated map JVP independent difference')
        call check(abs(dot_product(seed, dy) - dot_product(adj, d)) < 2.0e-12_dp, &
            'generated map VJP pairing')
    end subroutine
end program
