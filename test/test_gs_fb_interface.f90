program test_gs_fb_interface
    !! Interface term of the moving source cut (DESIGN 1.2): the flow
    !! parametrisation sigma int g/|grad psi| dl with g = f_edge v_i v_j
    !! against 1D Gauss quadrature along an analytically known cut, on an
    !! affine and on a curved cell, for a straight and a circular cut.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_quadrature, only: gauss_legendre
    use kin6d_triangle_basis, only: triangle_basis_jet
    use kin6d_gs_free_boundary, only: gs_fb_interface_mass
    implicit none

    real(dp), parameter :: mu0 = 1.2566370614359172e-6_dp
    integer :: failures

    failures = 0
    call straight_cut(2, .true.)
    call straight_cut(3, .true.)
    call straight_cut(2, .false.)
    call circular_cut(2, .true.)
    call circular_cut(3, .false.)
    if (failures /= 0) error stop 'gs_fb_interface failures'
    print '(a)', 'gs_fb_interface oracles passed'

contains

    subroutine require(ok, what)
        logical, intent(in) :: ok
        character(len=*), intent(in) :: what
        if (.not. ok) then
            failures = failures + 1
            print '(2a)', 'FAIL: ', what
        end if
    end subroutine require

    subroutine nodes(order, curved, pts)
        !! Reference node positions of the P-order Lagrange triangle mapped by
        !! a straight or a smoothly curved map.
        integer, intent(in) :: order
        logical, intent(in) :: curved
        real(dp), intent(out) :: pts(2, 10)
        real(dp) :: xi(2, 10), x(2)
        integer :: j, ndof

        ndof = (order + 1)*(order + 2)/2
        call reference_nodes(order, xi(:, :ndof))
        do j = 1, ndof
            x = [3.0_dp + 0.6_dp*xi(1, j) - 0.1_dp*xi(2, j), &
                 0.7_dp*xi(2, j) + 0.1_dp*xi(1, j)]
            if (curved) x = x + [0.05_dp*xi(1, j)*xi(2, j), &
                                 0.04_dp*xi(1, j)*(1.0_dp - xi(1, j) - xi(2, j))]
            pts(:, j) = x
        end do
    end subroutine nodes

    subroutine reference_nodes(order, xi)
        integer, intent(in) :: order
        real(dp), intent(out) :: xi(:, :)
        integer :: k, j
        real(dp), parameter :: v(2, 3) = reshape([0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, &
                                                  0.0_dp, 1.0_dp], [2, 3])

        xi(:, 1:3) = v
        k = 3
        do j = 1, 3
            do k = 1, order - 1
                xi(:, 3 + (j - 1)*(order - 1) + k) = v(:, j) + (v(:, mod(j, 3) + 1) &
                                                            - v(:, j))*real(k, dp)/order
            end do
        end do
        if (order == 3) xi(:, 10) = [1.0_dp/3, 1.0_dp/3]
    end subroutine reference_nodes

    subroutine map_point(order, pts, xi, x, det)
        integer, intent(in) :: order
        real(dp), intent(in) :: pts(2, 10), xi(2)
        real(dp), intent(out) :: x(2), det
        real(dp) :: n(10), g(2, 10), jac(2, 2)
        integer :: ndof

        ndof = (order + 1)*(order + 2)/2
        call triangle_basis_jet(xi, n(:ndof), g(:, :ndof))
        x = matmul(pts(:, :ndof), n(:ndof))
        jac = matmul(pts(:, :ndof), transpose(g(:, :ndof)))
        det = jac(1, 1)*jac(2, 2) - jac(1, 2)*jac(2, 1)
    end subroutine map_point

    subroutine reference_mass(order, pts, phi_fn, a, b, speed_fn, pp, ff, sigma, mass)
        !! 1D Gauss along the reference arc xi(s), s in [a, b].
        integer, intent(in) :: order
        real(dp), intent(in) :: pts(2, 10), a, b, pp, ff, sigma
        interface
            function phi_fn(s) result(xi)
                import :: dp
                real(dp), intent(in) :: s
                real(dp) :: xi(2)
            end function phi_fn
            function speed_fn(s) result(w)
                import :: dp
                real(dp), intent(in) :: s
                real(dp) :: w
            end function speed_fn
        end interface
        real(dp), intent(out) :: mass(:, :)
        real(dp) :: gx(40), gw(40), s, x(2), det, n(10), g(2, 10), xi(2)
        integer :: q, ndof, i, j

        ndof = (order + 1)*(order + 2)/2
        call gauss_legendre(40, gx, gw)
        mass = 0.0_dp
        do q = 1, 40
            s = a + (b - a)*0.5_dp*(gx(q) + 1.0_dp)
            xi = phi_fn(s)
            call map_point(order, pts, xi, x, det)
            call triangle_basis_jet(xi, n(:ndof), g(:, :ndof))
            do j = 1, ndof
                do i = 1, ndof
                    mass(i, j) = mass(i, j) + 0.5_dp*(b - a)*gw(q)*speed_fn(s) &
                                 *(ff/x(1) + mu0*pp*x(1))*det*n(i)*n(j)
                end do
            end do
        end do
        mass = sigma*mass
    end subroutine reference_mass

    subroutine straight_cut(order, affine)
        integer, intent(in) :: order
        logical, intent(in) :: affine
        real(dp) :: pts(2, 10), phi(10), xi_n(2, 10), mass(10, 10), ref(10, 10), dropped
        integer :: j, ndof, stat

        ndof = (order + 1)*(order + 2)/2
        call nodes(order,.not. affine, pts)
        call reference_nodes(order, xi_n(:, :ndof))
        ! phi = xi + eta/2 - 2/5, the cut runs from (2/5, 0) to (0, 4/5).
        do j = 1, ndof
            phi(j) = xi_n(1, j) + 0.5_dp*xi_n(2, j) - 0.4_dp
        end do
        call gs_fb_interface_mass(pts(:, :ndof), affine, phi(:ndof), mu0, 2.0e5_dp, &
                                  1.5_dp, -1.0_dp, mass(:ndof, :ndof), dropped, j)
        call require(j == 0 .and. dropped == 0.0_dp, 'straight cut traced')
        ! Reference arc length is sqrt(0.4**2 + 0.8**2); the flow speed is
        ! |grad phi| = sqrt(1.25), so dt = d sigma / sqrt(1.25).
        call reference_mass(order, pts, line, 0.0_dp, 1.0_dp, line_speed, 2.0e5_dp, &
                            1.5_dp, -1.0_dp, ref(:ndof, :ndof))
       print '(a,i1,l2,es10.2)', 'straight cut order/affine/max diff ', order, affine, &
            maxval(abs(mass(:ndof, :ndof) - ref(:ndof, :ndof)))
        call require(maxval(abs(mass(:ndof, :ndof) - ref(:ndof, :ndof))) &
                     <= 1.0e-10_dp*maxval(abs(ref(:ndof, :ndof))), &
                     'straight cut interface mass')
    end subroutine straight_cut

    subroutine circular_cut(order, affine)
        integer, intent(in) :: order
        logical, intent(in) :: affine
        real(dp) :: pts(2, 10), phi(10), xi_n(2, 10), mass(10, 10), ref(10, 10), dropped
        real(dp) :: theta1, theta2
        integer :: j, ndof

        ndof = (order + 1)*(order + 2)/2
        call nodes(order,.not. affine, pts)
        call reference_nodes(order, xi_n(:, :ndof))
        do j = 1, ndof
            phi(j) = (xi_n(1, j) - 0.15_dp)**2 + (xi_n(2, j) - 0.15_dp)**2 - 0.09_dp
        end do
        theta1 = atan2(-0.15_dp, sqrt(0.09_dp - 0.0225_dp))
        theta2 = atan2(sqrt(0.09_dp - 0.0225_dp), -0.15_dp)
        call gs_fb_interface_mass(pts(:, :ndof), affine, phi(:ndof), mu0, 0.0_dp, &
                                  1.0_dp, 1.0_dp, mass(:ndof, :ndof), dropped, j)
        call require(j == 0 .and. dropped == 0.0_dp, 'circular cut traced')
        call reference_mass(order, pts, arc, theta1, theta2, arc_speed, 0.0_dp, &
                            1.0_dp, 1.0_dp, ref(:ndof, :ndof))
       print '(a,i1,l2,es10.2)', 'circular cut order/affine/max diff ', order, affine, &
            maxval(abs(mass(:ndof, :ndof) - ref(:ndof, :ndof)))
        call require(maxval(abs(mass(:ndof, :ndof) - ref(:ndof, :ndof))) &
                     <= 1.0e-10_dp*maxval(abs(ref(:ndof, :ndof))), &
                     'circular cut interface mass')
    end subroutine circular_cut

    function line(s) result(xi)
        real(dp), intent(in) :: s
        real(dp) :: xi(2)
        xi = [0.4_dp*(1.0_dp - s), 0.8_dp*s]
    end function line
    function line_speed(s) result(w)
        real(dp), intent(in) :: s
        real(dp) :: w
        ! d sigma = sqrt(0.4^2 + 0.8^2) ds, dt = d sigma / |grad phi|.
        w = sqrt(0.16_dp + 0.64_dp)/sqrt(1.25_dp)
        if (s < 0.0_dp) w = 0.0_dp
    end function line_speed

    function arc(s) result(xi)
        real(dp), intent(in) :: s
        real(dp) :: xi(2)
        xi = [0.15_dp + 0.3_dp*cos(s), 0.15_dp + 0.3_dp*sin(s)]
    end function arc
    function arc_speed(s) result(w)
        real(dp), intent(in) :: s
        real(dp) :: w
        ! d sigma = r d theta, |grad phi| = 2 r.
        w = 0.3_dp/0.6_dp
        if (s > 1.0e9_dp) w = 0.0_dp
    end function arc_speed

end program test_gs_fb_interface
