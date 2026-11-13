module milca_ideal_flux_tokamak
    ! Axisymmetric equilibrium from a poloidal-flux grid (G-EQDSK) as a chart
    ! for milca_ideal_flux (docs/stability.md, P1e.3). rho = sqrt(psin/psin_max)
    ! with the wall on the surface psin_max; theta is the straight-field-line
    ! angle (PEST or Hamada) of milca_vacuum_surfaces and zeta the straight
    ! toroidal angle, phi = zeta + nu - q theta; with angle = 0 theta is the
    ! geometric angle about the axis, zeta = phi and the stream function
    ! lambda = theta_PEST - theta makes theta + lambda straight. R, Z, the
    ! shift nu - q theta and lambda (unscaled, finite on the axis)
    ! are Fourier series in theta whose coefficients, divided by rho^k with
    ! k = min(m, 2 - mod(m, 2)) (lambda: k = mod(m, 2)), are fitted by quintic
    ! B-splines in s = rho^2 (least squares over sampled surfaces).
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_bspline, only: bspline_workspace_t, bspline_init, &
        bspline_set_knots, bspline_eval_deriv
    use fortnum_status, only: fortnum_status_t, FORTNUM_OK
    use milca_ideal_flux, only: flux_equilibrium_t, flux_point_t
    use milca_ideal_flux_leaves, only: flux_toroidal_chart
    use milca_vacuum_surfaces, only: psi_grid_t, surface_t, surface_build, &
        profiles_of_psin, psi_eval
    implicit none
    private
    public :: tokamak_chart_t, tokamak_chart_init, tokamak_grid_field, &
        tokamak_position

    real(dp), parameter :: twopi = 2*acos(-1.0_dp)

    type, extends(flux_equilibrium_t) :: tokamak_chart_t
        type(psi_grid_t) :: grid
        integer :: mmax = 0
        real(dp) :: psin_max = 0, dpsi = 0, density = 1
        real(dp) :: psi_sign = 1, iota_sign = 1
        type(bspline_workspace_t) :: ws
        !> Columns: R cos, R sin, Z cos, Z sin, shift cos, shift sin, lambda cos,
        !> lambda sin for m = 0..mmax, then |q|.
        real(dp), allocatable :: coef(:, :)
    contains
        procedure :: evaluate => tokamak_evaluate
    end type tokamak_chart_t
contains
    subroutine tokamak_chart_init(chart, grid, psin_max, angle, mmax, nbreak, &
            nsample, ntheta, density, status)
        type(tokamak_chart_t), intent(out) :: chart
        type(psi_grid_t), intent(in) :: grid
        real(dp), intent(in) :: psin_max, density
        integer, intent(in) :: angle, mmax, nbreak, nsample, ntheta
        integer, intent(out) :: status
        type(surface_t) :: surf
        type(fortnum_status_t) :: fs
        type(flux_point_t) :: pt
        real(dp), allocatable :: a(:, :), y(:, :), basis(:, :), work(:)
        real(dp) :: s, rho, ts, w, scale, dual(3, 3), jac, bvec(3), lam, shift
        integer :: i, j, m, ncol, info, ok
        status = 1
        if (psin_max <= 0 .or. psin_max >= 1 .or. mmax < 0) return
        if (nbreak < 2 .or. ntheta < 8) return
        chart%grid = grid; chart%mmax = mmax; chart%psin_max = psin_max
        chart%dpsi = abs(grid%psi_edge - grid%psi_axis); chart%density = density
        call bspline_init(chart%ws, 6, nbreak, fs)
        if (fs%code /= FORTNUM_OK) return
        call bspline_set_knots(chart%ws, [(real(i, dp)/(nbreak - 1), i = 0, &
            nbreak - 1)], fs)
        if (fs%code /= FORTNUM_OK) return
        if (nsample < chart%ws%ncoef) return
        ncol = 8*(mmax + 1) + 1
        allocate(a(nsample, chart%ws%ncoef), y(nsample, ncol), &
            basis(0:0, chart%ws%ncoef))
        do i = 1, nsample
            s = (i - 0.5_dp)/nsample; rho = sqrt(s)
            call bspline_eval_deriv(chart%ws, s, 0, basis, fs)
            a(i, :) = basis(0, :)
            call surface_build(grid, psin_max*s, ntheta, surf, ok)
            if (ok /= 0) return
            y(i, :) = 0
            do j = 1, ntheta
                if (angle == 0) then
                    ts = surf%theta_g(j)
                    w = 2.0_dp/ntheta
                    lam = surf%theta_s(j, 1) - ts
                    shift = 0
                else
                    ts = surf%theta_s(j, angle)
                    w = surf%dtheta_s(j, angle)*twopi/ntheta/acos(-1.0_dp)
                    lam = 0
                    shift = surf%nu(j) - surf%q*ts
                end if
                do m = 0, mmax
                    y(i, 7 + 6*mmax + m) = y(i, 7 + 6*mmax + m) &
                        + merge(w/2, w, m == 0)*lam*cos(m*ts)/rho**mod(m, 2)
                    y(i, 8 + 7*mmax + m) = y(i, 8 + 7*mmax + m) &
                        + w*lam*sin(m*ts)/rho**mod(m, 2)
                    scale = w/rho**power(m)
                    if (m == 0) scale = scale/2
                    y(i, 1 + m) = y(i, 1 + m) + scale*surf%r(j)*cos(m*ts)
                    y(i, 2 + mmax + m) = y(i, 2 + mmax + m) + scale*surf%r(j)*sin(m*ts)
                    y(i, 3 + 2*mmax + m) = y(i, 3 + 2*mmax + m) &
                        + scale*surf%z(j)*cos(m*ts)
                    y(i, 4 + 3*mmax + m) = y(i, 4 + 3*mmax + m) &
                        + scale*surf%z(j)*sin(m*ts)
                    y(i, 5 + 4*mmax + m) = y(i, 5 + 4*mmax + m) &
                        + scale*shift*cos(m*ts)
                    y(i, 6 + 5*mmax + m) = y(i, 6 + 5*mmax + m) &
                        + scale*shift*sin(m*ts)
                end do
            end do
            y(i, ncol) = abs(surf%q)
        end do
        allocate(work(64*(nsample + ncol)))
        call dgels('N', nsample, chart%ws%ncoef, ncol, a, nsample, y, nsample, &
            work, size(work), info)
        if (info /= 0) return
        chart%coef = y(1:chart%ws%ncoef, :)
        ! Orientation: match sqrt(g) B^zeta and B^theta/B^zeta with the grid
        ! field B = (s/R)(psi_Z, -psi_R) + F/R e_phi at one point.
        chart%psi_sign = 1; chart%iota_sign = 1
        call chart%evaluate(0.5_dp, 0.3_dp, 0.0_dp, pt)
        jac = dot_product(pt%xu(:, 1), cross(pt%xu(:, 2), pt%xu(:, 3)))
        dual(:, 1) = cross(pt%xu(:, 2), pt%xu(:, 3))/jac
        dual(:, 2) = cross(pt%xu(:, 3), pt%xu(:, 1))/jac
        dual(:, 3) = cross(pt%xu(:, 1), pt%xu(:, 2))/jac
        bvec = tokamak_grid_field(chart, 0.5_dp, 0.3_dp, 0.0_dp)
        chart%psi_sign = sign(1.0_dp, jac*dot_product(bvec, dual(:, 3)))
        chart%iota_sign = sign(1.0_dp, dot_product(bvec, dual(:, 2)) &
            /dot_product(bvec, dual(:, 3)))
        status = 0
    end subroutine tokamak_chart_init

    subroutine tokamak_evaluate(self, rho, theta, zeta, point)
        class(tokamak_chart_t), intent(in) :: self
        real(dp), intent(in) :: rho, theta, zeta
        type(flux_point_t), intent(out) :: point
        real(dp) :: col(0:2, size(self%coef, 2)), s
        real(dp) :: rj(0:9), zj(0:9), fj(0:9), lj(0:9), q, qr, f, df, p, dpx, &
            psin, dps
        s = rho*rho
        call synthesize(self, rho, theta, zeta, rj, zj, fj, col, lj)
        call flux_toroidal_chart(rj(0), rj(1), rj(2), rj(3), rj(4), rj(5), rj(6), &
            rj(7), rj(8), rj(9), zj(1), zj(2), zj(3), zj(4), zj(5), zj(6), zj(7), &
            zj(8), zj(9), fj(0), fj(1), fj(2), fj(3), fj(4), fj(5), fj(6), fj(7), &
            fj(8), fj(9), &
            point%xu(1,1), point%xu(2,1), point%xu(3,1), point%xu(1,2), &
            point%xu(2,2), point%xu(3,2), point%xu(1,3), point%xu(2,3), &
            point%xu(3,3), &
            point%xuu(1,1,1), point%xuu(2,1,1), point%xuu(3,1,1), &
            point%xuu(1,1,2), point%xuu(2,1,2), point%xuu(3,1,2), &
            point%xuu(1,2,2), point%xuu(2,2,2), point%xuu(3,2,2), &
            point%xuu(1,1,3), point%xuu(2,1,3), point%xuu(3,1,3), &
            point%xuu(1,2,3), point%xuu(2,2,3), point%xuu(3,2,3), &
            point%xuu(1,3,3), point%xuu(2,3,3), point%xuu(3,3,3))
        call symmetrize(point%xuu)
        ! |q| and its rho derivatives from the s fit.
        q = col(0, size(col, 2)); qr = 2*rho*col(1, size(col, 2))
        dps = 2*self%psin_max*self%dpsi
        point%psi1 = self%psi_sign*dps*rho*q
        point%psi2 = self%psi_sign*dps*(q + rho*qr)
        point%iota = self%iota_sign/q
        point%iota1 = -self%iota_sign*qr/q**2
        ! lambda_theta, lambda_zeta, d(rho, theta, zeta) lambda_theta,
        ! d(rho, zeta) lambda_zeta; lambda does not depend on zeta here.
        point%lam = [lj(2), 0.0_dp, lj(5), lj(6), 0.0_dp, 0.0_dp, 0.0_dp]
        psin = self%psin_max*s
        call profiles_of_psin(self%grid, psin, f, df, p, dpx)
        point%p = p; point%pp = dpx*2*rho*self%psin_max
        point%density = self%density
    end subroutine tokamak_evaluate

    subroutine synthesize(self, rho, theta, zeta, rj, zj, fj, col, lj)
        ! Jets of R, Z and phi: index 0 value, 1:3 d(rho, theta, zeta), 4:9 the
        ! second derivatives (11, 12, 22, 13, 23, 33).
        class(tokamak_chart_t), intent(in) :: self
        real(dp), intent(in) :: rho, theta, zeta
        real(dp), intent(out) :: rj(0:9), zj(0:9), fj(0:9), col(0:, :)
        real(dp), intent(out), optional :: lj(0:9)
        type(fortnum_status_t) :: fs
        real(dp) :: lt(0:9)
        real(dp) :: basis(0:2, self%ws%ncoef), s
        integer :: m
        s = rho*rho
        call bspline_eval_deriv(self%ws, s, 2, basis, fs)
        col = matmul(basis, self%coef)
        rj = 0; zj = 0; fj = 0; lt = 0
        do m = 0, self%mmax
            call add_mode(rj, m, power(m), rho, theta, col(:, 1 + m), &
                col(:, 2 + self%mmax + m))
            call add_mode(zj, m, power(m), rho, theta, col(:, 3 + 2*self%mmax + m), &
                col(:, 4 + 3*self%mmax + m))
            call add_mode(fj, m, power(m), rho, theta, col(:, 5 + 4*self%mmax + m), &
                col(:, 6 + 5*self%mmax + m))
            call add_mode(lt, m, mod(m, 2), rho, theta, col(:, 7 + 6*self%mmax + m), &
                col(:, 8 + 7*self%mmax + m))
        end do
        if (present(lj)) lj = lt
        fj(0) = fj(0) + zeta; fj(3) = fj(3) + 1
    end subroutine synthesize

    function tokamak_position(self, rho, theta, zeta) result(x)
        ! (R, Z, phi) of the chart point.
        class(tokamak_chart_t), intent(in) :: self
        real(dp), intent(in) :: rho, theta, zeta
        real(dp) :: x(3), rj(0:9), zj(0:9), fj(0:9), col(0:2, size(self%coef, 2))
        call synthesize(self, rho, theta, zeta, rj, zj, fj, col)
        x = [rj(0), zj(0), fj(0)]
    end function tokamak_position

    function tokamak_grid_field(self, rho, theta, zeta) result(b)
        ! Cartesian B of the flux grid at the chart point:
        ! B = (s/R)(psi_Z e_R - psi_R e_Z) + (F/R) e_phi.
        class(tokamak_chart_t), intent(in) :: self
        real(dp), intent(in) :: rho, theta, zeta
        real(dp) :: b(3), rj(0:9), zj(0:9), fj(0:9), col(0:2, size(self%coef, 2))
        real(dp) :: psi, pr, pz, f, df, p, dpx, r
        call synthesize(self, rho, theta, zeta, rj, zj, fj, col)
        r = rj(0)
        call psi_eval(self%grid, r, zj(0), psi, pr, pz)
        call profiles_of_psin(self%grid, (psi - self%grid%psi_axis) &
            /(self%grid%psi_edge - self%grid%psi_axis), f, df, p, dpx)
        b = self%grid%s_pol*pz/r*[cos(fj(0)), sin(fj(0)), 0.0_dp] &
            + f/r*[-sin(fj(0)), cos(fj(0)), 0.0_dp] &
            + [0.0_dp, 0.0_dp, -self%grid%s_pol*pr/r]
    end function tokamak_grid_field

    pure subroutine add_mode(x, m, k, rho, theta, c, sn)
        ! Adds rho^k [a_c(s) cos(m theta) + a_s(s) sin(m theta)], s = rho^2, to
        ! the jet x(0:9); c and sn hold a(s), a'(s), a''(s).
        real(dp), intent(inout) :: x(0:9)
        integer, intent(in) :: m, k
        real(dp), intent(in) :: rho, theta, c(0:2), sn(0:2)
        real(dp) :: f0, f1, f2, cs, ss
        cs = cos(m*theta); ss = sin(m*theta)
        call radial(k, rho, c, f0, f1, f2)
        x(0) = x(0) + f0*cs; x(1) = x(1) + f1*cs; x(4) = x(4) + f2*cs
        x(2) = x(2) - m*f0*ss; x(5) = x(5) - m*f1*ss; x(6) = x(6) - m*m*f0*cs
        call radial(k, rho, sn, f0, f1, f2)
        x(0) = x(0) + f0*ss; x(1) = x(1) + f1*ss; x(4) = x(4) + f2*ss
        x(2) = x(2) + m*f0*cs; x(5) = x(5) + m*f1*cs; x(6) = x(6) - m*m*f0*ss
    end subroutine add_mode

    pure integer function power(m)
        ! A harmonic m of a regular function near the axis vanishes at least
        ! like rho^power with the parity of m, in every regular angle, also the
        ! geometric angle about an elliptic axis (odd m like rho, even m like
        ! rho^2). Dividing by more would amplify rounding noise.
        integer, intent(in) :: m
        power = min(m, 2 - mod(m, 2))
    end function power

    pure subroutine radial(k, rho, a, v0, v1, v2)
        ! f(rho) = rho^k a(rho^2) and its first two rho derivatives.
        integer, intent(in) :: k
        real(dp), intent(in) :: rho, a(0:2)
        real(dp), intent(out) :: v0, v1, v2
        real(dp) :: rm, rm1, rm2
        rm = rho**k
        rm1 = 0
        if (k >= 1) rm1 = k*rho**(k - 1)
        rm2 = 0
        if (k >= 2) rm2 = k*(k - 1)*rho**(k - 2)
        v0 = rm*a(0)
        v1 = rm1*a(0) + 2*rho*rm*a(1)
        v2 = rm2*a(0) + 4*rho*rm1*a(1) + 2*rm*a(1) + 4*rho*rho*rm*a(2)
    end subroutine radial

    pure subroutine symmetrize(xuu)
        real(dp), intent(inout) :: xuu(3, 3, 3)
        integer :: u, v
        do v = 1, 3
            do u = v + 1, 3
                xuu(:, u, v) = xuu(:, v, u)
            end do
        end do
    end subroutine symmetrize

    pure function cross(p, q) result(w)
        real(dp), intent(in) :: p(3), q(3)
        real(dp) :: w(3)
        w = [p(2)*q(3) - p(3)*q(2), p(3)*q(1) - p(1)*q(3), p(1)*q(2) - p(2)*q(1)]
    end function cross
end module milca_ideal_flux_tokamak
