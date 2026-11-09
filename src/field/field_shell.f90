module milca_field_shell
    ! Field of the free-boundary vacuum shell at physical points. The shell is
    ! a curved triangle mesh in the reference plane times a real Fourier
    ! series in phi_hat (mode 1 constant, then cos, sin per harmonic k, angle
    ! k nfp phi_hat). Its map xs gives the chart point q, and the potential
    !   Phi = phi_s + f_vac q3 + circulation theta/(2 pi)
    ! gives B = grad Phi. The kernels take the arrays explicitly so the solver
    ! state (milca_ideal_free_boundary) and a client's copy share them.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use milca_equilibrium_chart, only: ideal_chart_t, cylindrical_chart_t
    use milca_triangle_basis, only: triangle_basis_jet
    implicit none
    private
    public :: field_shell_jet, field_shell_eval, field_shell_locate
    real(dp), parameter :: pi = acos(-1.0_dp)
contains
    subroutine field_shell_jet(nfp, conn, xs, phi, xi, phi_hat, q, qd, gt)
        ! Chart point q of the cell with nodes conn at (xi, phi_hat), its
        ! derivatives qd(:, j) along (xi1, xi2, phi_hat) and the master
        ! gradient gt of the single-valued potential phi_s.
        integer, intent(in) :: nfp, conn(:)
        real(dp), intent(in) :: xs(:, :, :), phi(:, :), xi(2), phi_hat
        real(dp), intent(out) :: q(3), qd(3, 3), gt(3)
        real(dp) :: shape(10), grad(2, 10), t, dt, nn
        integer :: a, m, node, nloc
        nloc = size(conn)
        call triangle_basis_jet(xi, shape(:nloc), grad(:, :nloc))
        q = 0; qd = 0; gt = 0
        do a = 1, nloc
            node = conn(a)
            do m = 1, size(xs, 2)
                nn = (m/2)*nfp
                if (m == 1) then
                    t = 1; dt = 0
                else if (mod(m, 2) == 0) then
                    t = cos(nn*phi_hat); dt = -nn*sin(nn*phi_hat)
                else
                    t = sin(nn*phi_hat); dt = nn*cos(nn*phi_hat)
                end if
                q = q + xs(:, m, node)*shape(a)*t
                qd(:, 1) = qd(:, 1) + xs(:, m, node)*grad(1, a)*t
                qd(:, 2) = qd(:, 2) + xs(:, m, node)*grad(2, a)*t
                qd(:, 3) = qd(:, 3) + xs(:, m, node)*shape(a)*dt
                gt(1) = gt(1) + phi(m, node)*grad(1, a)*t
                gt(2) = gt(2) + phi(m, node)*grad(2, a)*t
                gt(3) = gt(3) + phi(m, node)*shape(a)*dt
            end do
        end do
        q(3) = q(3) + phi_hat; qd(3, 3) = qd(3, 3) + 1
    end subroutine

    subroutine field_shell_eval(chart, nfp, conn, reference, xs, phi, f_vac, &
            circulation, theta_centre, xi, phi_hat, x, b, status)
        ! Physical point x = (R, Z, phi) and cylindrical field (B_R, B_phi,
        ! B_Z) at master point xi, phi_hat of the cell with nodes conn.
        class(ideal_chart_t), intent(in) :: chart
        integer, intent(in) :: nfp, conn(:)
        real(dp), intent(in) :: reference(:, :), xs(:, :, :), phi(:, :)
        real(dp), intent(in) :: f_vac, circulation, theta_centre(2), xi(2), phi_hat
        real(dp), intent(out) :: x(3), b(3)
        integer, intent(out) :: status
        real(dp) :: shape(10), grad(2, 10), q(3), qd(3, 3), gt(3), d(2), tg(2)
        real(dp) :: dg(2, 2), y(2), gg(3), cart(3), d1(3, 3), d2(3, 3, 3), dx(3, 3)
        real(dp) :: d3(3, 3, 3, 3), inv(3, 3), bc(3), det
        integer :: nloc
        status = 1; x = 0; b = 0
        nloc = size(conn)
        call field_shell_jet(nfp, conn, xs, phi, xi, phi_hat, q, qd, gt)
        call triangle_basis_jet(xi, shape(:nloc), grad(:, :nloc))
        y = matmul(reference(:, conn), shape(:nloc))
        dg = matmul(reference(:, conn), transpose(grad(:, :nloc)))
        d = y - theta_centre
        tg = [-d(2), d(1)]/(2*pi*sum(d**2))
        gg = gt + f_vac*qd(3, :) + circulation*[matmul(transpose(dg), tg), 0.0_dp]
        call chart%jet(q, cart, d1, d2, d3)
        dx = matmul(d1, qd)
        call inverse3(dx, inv, det)
        if (.not. abs(det) > 0) return
        ! Cartesian B = DX^{-T} g, then cylindrical components.
        bc = matmul(transpose(inv), gg)
        select type (chart)
        type is (cylindrical_chart_t)
            x = q
        class default
            ! Cylindrical coordinates of the Cartesian point, the angle
            ! unwrapped to the branch nearest the chart angle q3.
            x(3) = atan2(cart(2), cart(1))
            x(3) = x(3) + 2*pi*nint((q(3) - x(3))/(2*pi))
            x(1) = hypot(cart(1), cart(2)); x(2) = cart(3)
        end select
        b(1) = bc(1)*cos(x(3)) + bc(2)*sin(x(3))
        b(2) = -bc(1)*sin(x(3)) + bc(2)*cos(x(3))
        b(3) = bc(3)
        status = 0
    end subroutine

    subroutine field_shell_locate(chart, nfp, conn, xs, x, seed, e, xi, phi_hat, status)
        ! Cell e and master coordinates (xi, phi_hat) of the physical point
        ! x = (R, Z, phi). Candidates are all cells whose node images at phi
        ! span a box (padded by half its extent, for curved edges) containing
        ! x, tried by Newton on the Cartesian position in order of centroid
        ! distance, seed first. status 0 found, 1 outside the shell (plasma,
        ! beyond Gamma_c or a fold).
        class(ideal_chart_t), intent(in) :: chart
        integer, intent(in) :: nfp, conn(:, :), seed
        real(dp), intent(in) :: xs(:, :, :), x(3)
        integer, intent(out) :: e, status
        real(dp), intent(out) :: xi(2), phi_hat
        real(dp) :: target(3), q(3), qd(3, 3), gt(3), cart(3), d1(3, 3)
        real(dp) :: d2(3, 3, 3), d3(3, 3, 3, 3), dx(3, 3), inv(3, 3), det, t, nn
        real(dp) :: residual(3), step(3), scale, lo(3), hi(3), pad(3)
        real(dp), allocatable :: dist(:), node_cart(:, :), dummy(:, :)
        integer :: iteration, c, ncell, node, m
        status = 1; e = 0; xi = 0; phi_hat = x(3)
        ncell = size(conn, 2)
        if (ncell < 1 .or. .not. x(1) > 0) return
        target = [x(1)*cos(x(3)), x(1)*sin(x(3)), x(2)]
        scale = max(1.0_dp, x(1), abs(x(2)))
        allocate (dummy(size(xs, 2), size(xs, 3)), source=0.0_dp)
        allocate (node_cart(3, size(xs, 3)), dist(ncell))
        do node = 1, size(xs, 3)
            q = xs(:, 1, node)
            do m = 2, size(xs, 2)
                nn = (m/2)*nfp
                t = merge(cos(nn*x(3)), sin(nn*x(3)), mod(m, 2) == 0)
                q = q + xs(:, m, node)*t
            end do
            q(3) = q(3) + x(3)
            call chart%jet(q, node_cart(:, node), d1, d2, d3)
        end do
        do c = 1, ncell
            lo = minval(node_cart(:, conn(:, c)), 2)
            hi = maxval(node_cart(:, conn(:, c)), 2)
            pad = 0.5_dp*maxval(hi - lo) + 1.0e-12_dp*scale
            if (any(target < lo - pad) .or. any(target > hi + pad)) then
                dist(c) = huge(1.0_dp)
            else
                dist(c) = sum((sum(node_cart(:, conn(:, c)), 2)/size(conn, 1) &
                    - target)**2)
            end if
        end do
        if (seed >= 1 .and. seed <= ncell) then
            if (dist(seed) < huge(1.0_dp)) dist(seed) = -1
        end if
        do
            c = minloc(dist, 1)
            if (.not. dist(c) < huge(1.0_dp)) exit
            dist(c) = huge(1.0_dp)
            xi = 1.0_dp/3; phi_hat = x(3)
            do iteration = 1, 30
                call field_shell_jet(nfp, conn(:, c), xs, dummy, xi, phi_hat, q, qd, gt)
                call chart%jet(q, cart, d1, d2, d3)
                residual = cart - target
                if (maxval(abs(residual)) < 1.0e-12_dp*scale) then
                    if (minval([xi, 1 - sum(xi)]) >= -1.0e-10_dp) then
                        e = c; status = 0; return
                    end if
                    exit
                end if
                dx = matmul(d1, qd)
                call inverse3(dx, inv, det)
                if (.not. abs(det) > 0) exit
                step = matmul(inv, residual)
                if (.not. all(ieee_is_finite(step))) exit
                step = step/max(1.0_dp, maxval(abs(step(:2))))
                xi = xi - step(:2); phi_hat = phi_hat - step(3)
                if (minval(xi) < -2 .or. maxval(xi) > 3) exit
            end do
        end do
        e = 0
    end subroutine

    pure subroutine inverse3(a, inv, det)
        real(dp), intent(in) :: a(3, 3)
        real(dp), intent(out) :: inv(3, 3), det
        inv(1, 1) = a(2, 2)*a(3, 3) - a(2, 3)*a(3, 2)
        inv(1, 2) = a(1, 3)*a(3, 2) - a(1, 2)*a(3, 3)
        inv(1, 3) = a(1, 2)*a(2, 3) - a(1, 3)*a(2, 2)
        inv(2, 1) = a(2, 3)*a(3, 1) - a(2, 1)*a(3, 3)
        inv(2, 2) = a(1, 1)*a(3, 3) - a(1, 3)*a(3, 1)
        inv(2, 3) = a(1, 3)*a(2, 1) - a(1, 1)*a(2, 3)
        inv(3, 1) = a(2, 1)*a(3, 2) - a(2, 2)*a(3, 1)
        inv(3, 2) = a(1, 2)*a(3, 1) - a(1, 1)*a(3, 2)
        inv(3, 3) = a(1, 1)*a(2, 2) - a(1, 2)*a(2, 1)
        det = a(1, 1)*inv(1, 1) + a(1, 2)*inv(2, 1) + a(1, 3)*inv(3, 1)
        if (abs(det) > 0) inv = inv/det
    end subroutine
end module milca_field_shell
