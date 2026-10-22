module kin6d_nested_chart
    ! Optional IH2 accelerator pilot: curve-following O-grid charts built from a
    ! nested solution, emitted as plain MF1 mesh arrays (vertices, geometry).
    !
    ! The nested map is a stellarator-symmetric Fourier series in the
    ! axis-regular disk coordinates (u, v) = rho (cos theta, sin theta):
    !   R + i Z = sum_j [a_j(s) Re + i c_j(s) Im] (u + i v)^m_j exp(-i n_j zeta),
    ! s = u^2 + v^2 and polynomial a_j, c_j, so every field is smooth through
    ! the magnetic axis. lambda uses the Z form. zeta is the cylindrical angle.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    private
    public :: nested_map_t, nested_point_t, nested_map_read, nested_map_eval
    public :: nested_map_invert, nested_field_direction, nested_forcing
    public :: build_ogrid_chart, CHART_ALIGNED, CHART_UNALIGNED, CHART_FIELDLINE
    integer, parameter :: CHART_ALIGNED = 1, CHART_UNALIGNED = 2, CHART_FIELDLINE = 3
    real(dp), parameter :: pi = acos(-1.0_dp)
    type :: nested_map_t
        integer :: nfp = 1, nmode = 0, nk = 0, sign_theta = 1
        integer, allocatable :: m(:), n(:)
        real(dp), allocatable :: rc(:, :), zs(:, :), ls(:, :), iota(:)
    end type
    type :: nested_point_t
        ! Cartesian position x, Jacobian a(:, q) = dx/dq for q = (u, v, zeta),
        ! Hessians h(:, :, i) of x_i and lambda jets in (u, v, zeta).
        real(dp) :: r = 0, z = 0, x(3) = 0, a(3, 3) = 0, h(3, 3, 3) = 0
        real(dp) :: rq(3) = 0, zq(3) = 0, lam(3) = 0
    end type
contains
    subroutine nested_map_read(path, map, status)
        character(*), intent(in) :: path
        type(nested_map_t), intent(out) :: map
        integer, intent(out) :: status
        character(256) :: header
        integer :: unit, j
        open (newunit=unit, file=path, status='old', action='read', iostat=status)
        if (status /= 0) return
        read (unit, '(a)', iostat=status) header
        if (status == 0 .and. index(header, 'kin6d nested map v1') == 0) status = 2
        if (status == 0) read (unit, *, iostat=status) map%nfp, map%nmode, map%nk
        if (status == 0) then
            if (map%nfp < 1 .or. map%nmode < 1 .or. map%nk < 1) status = 3
        end if
        if (status == 0) then
            allocate (map%m(map%nmode), map%n(map%nmode), map%iota(map%nk), &
                      map%rc(map%nk, map%nmode), map%zs(map%nk, map%nmode), &
                      map%ls(map%nk, map%nmode))
            read (unit, *, iostat=status) map%iota
        end if
        do j = 1, map%nmode
            if (status /= 0) exit
            read (unit, *, iostat=status) map%m(j), map%n(j), map%rc(:, j), &
                map%zs(:, j), map%ls(:, j)
        end do
        close (unit)
        if (status /= 0) return
        if (any(map%m < 0)) status = 4
        if (status == 0) call orient_theta(map)
    end subroutine

    subroutine orient_theta(map)
        ! Choose the logical poloidal direction so that the chart reference
        ! frame (radial, tangent, zeta) is positively oriented in (x, y, z).
        type(nested_map_t), intent(inout) :: map
        type(nested_point_t) :: pt
        real(dp) :: drho(2), dth(2)
        map%sign_theta = 1
        call nested_map_eval(map, 0.5_dp, 0.0_dp, 0.0_dp, pt)
        drho = [pt%rq(1), pt%zq(1)]
        dth = [pt%rq(2), pt%zq(2)]*0.5_dp
        ! det d(R,Z)/d(rho,theta) must be negative (see build_ogrid_chart).
        if (drho(1)*dth(2) - drho(2)*dth(1) > 0) map%sign_theta = -1
    end subroutine

    pure subroutine nested_map_eval(map, u, v, zeta, pt)
        type(nested_map_t), intent(in) :: map
        real(dp), intent(in) :: u, v, zeta
        type(nested_point_t), intent(out) :: pt
        complex(dp) :: w, w1, w2, wz, e, c, fr(10), fz(10), fl(10)
        real(dp) :: s, sk, ar(0:2), az(0:2), al(0:2), cz, sz
        real(dp) :: r1(3), r2(3, 3)
        integer :: j, k, mm, nn
        s = u*u + v*v
        c = cmplx(u, v, dp)
        fr = 0; fz = 0; fl = 0
        do j = 1, map%nmode
            mm = map%m(j)
            nn = map%n(j)
            e = exp(cmplx(0.0_dp, -nn*zeta, dp))
            w = c**mm*e
            w1 = 0; w2 = 0
            if (mm >= 1) w1 = mm*c**(mm - 1)*e
            if (mm >= 2) w2 = mm*(mm - 1)*c**(mm - 2)*e
            wz = cmplx(0.0_dp, -real(nn, dp), dp)
            ar = 0; az = 0; al = 0
            sk = 1
            do k = 1, map%nk
                ! Value and first two s-derivatives of each radial polynomial.
                ar(0) = ar(0) + map%rc(k, j)*sk
                az(0) = az(0) + map%zs(k, j)*sk
                al(0) = al(0) + map%ls(k, j)*sk
                if (k >= 2) then
                    ar(1) = ar(1) + (k - 1)*map%rc(k, j)*s**(k - 2)
                    az(1) = az(1) + (k - 1)*map%zs(k, j)*s**(k - 2)
                    al(1) = al(1) + (k - 1)*map%ls(k, j)*s**(k - 2)
                end if
                if (k >= 3) then
                    ar(2) = ar(2) + (k - 1)*(k - 2)*map%rc(k, j)*s**(k - 3)
                    az(2) = az(2) + (k - 1)*(k - 2)*map%zs(k, j)*s**(k - 3)
                    al(2) = al(2) + (k - 1)*(k - 2)*map%ls(k, j)*s**(k - 3)
                end if
                sk = sk*s
            end do
            call accumulate(ar, fr)
            call accumulate(az, fz)
            call accumulate(al, fl)
        end do
        ! fr holds R jets in Re, fz/fl hold Z/lambda jets in Im:
        ! order (val, u, v, zeta, uu, vv, uv, uz, vz, zz).
        pt%r = real(fr(1), dp)
        pt%z = aimag(fz(1))
        pt%rq = real(fr(2:4), dp)
        pt%zq = aimag(fz(2:4))
        pt%lam = aimag(fl(2:4))
        r1 = pt%rq
        r2(1, 1) = real(fr(5), dp); r2(2, 2) = real(fr(6), dp)
        r2(1, 2) = real(fr(7), dp); r2(2, 1) = r2(1, 2)
        r2(1, 3) = real(fr(8), dp); r2(3, 1) = r2(1, 3)
        r2(2, 3) = real(fr(9), dp); r2(3, 2) = r2(2, 3)
        r2(3, 3) = real(fr(10), dp)
        cz = cos(zeta); sz = sin(zeta)
        pt%x = [pt%r*cz, pt%r*sz, pt%z]
        pt%a(1, :) = r1*cz
        pt%a(2, :) = r1*sz
        pt%a(1, 3) = pt%a(1, 3) - pt%r*sz
        pt%a(2, 3) = pt%a(2, 3) + pt%r*cz
        pt%a(3, :) = pt%zq
        pt%h(:, :, 1) = r2*cz
        pt%h(:, :, 2) = r2*sz
        pt%h(1:2, 3, 1) = pt%h(1:2, 3, 1) - r1(1:2)*sz
        pt%h(3, 1:2, 1) = pt%h(1:2, 3, 1)
        pt%h(1:2, 3, 2) = pt%h(1:2, 3, 2) + r1(1:2)*cz
        pt%h(3, 1:2, 2) = pt%h(1:2, 3, 2)
        pt%h(3, 3, 1) = r2(3, 3)*cz - 2*r1(3)*sz - pt%r*cz
        pt%h(3, 3, 2) = r2(3, 3)*sz + 2*r1(3)*cz - pt%r*sz
        pt%h(1, 1, 3) = aimag(fz(5)); pt%h(2, 2, 3) = aimag(fz(6))
        pt%h(1, 2, 3) = aimag(fz(7)); pt%h(2, 1, 3) = pt%h(1, 2, 3)
        pt%h(1, 3, 3) = aimag(fz(8)); pt%h(3, 1, 3) = pt%h(1, 3, 3)
        pt%h(2, 3, 3) = aimag(fz(9)); pt%h(3, 2, 3) = pt%h(2, 3, 3)
        pt%h(3, 3, 3) = aimag(fz(10))
    contains
        pure subroutine accumulate(a, f)
            real(dp), intent(in) :: a(0:2)
            complex(dp), intent(inout) :: f(10)
            f(1) = f(1) + a(0)*w
            f(2) = f(2) + 2*u*a(1)*w + a(0)*w1
            f(3) = f(3) + 2*v*a(1)*w + a(0)*cmplx(0, 1, dp)*w1
            f(4) = f(4) + a(0)*wz*w
            f(5) = f(5) + 4*u*u*a(2)*w + 2*a(1)*w + 4*u*a(1)*w1 + a(0)*w2
            f(6) = f(6) + 4*v*v*a(2)*w + 2*a(1)*w &
                   + 4*v*a(1)*cmplx(0, 1, dp)*w1 - a(0)*w2
            f(7) = f(7) + 4*u*v*a(2)*w + 2*u*a(1)*cmplx(0, 1, dp)*w1 &
                   + 2*v*a(1)*w1 + a(0)*cmplx(0, 1, dp)*w2
            f(8) = f(8) + 2*u*a(1)*wz*w + a(0)*wz*w1
            f(9) = f(9) + 2*v*a(1)*wz*w + a(0)*wz*cmplx(0, 1, dp)*w1
            f(10) = f(10) + a(0)*wz*wz*w
        end subroutine
    end subroutine

    subroutine nested_map_invert(map, x, q, status)
        ! Inverse (u, v) of (R, Z)(u, v; zeta) = target, zeta from x. q holds a
        ! warm start on entry and (u, v, zeta) on exit. The fitted map folds
        ! outside s = 1, so a root is accepted only for s <= 1 + margin (points
        ! of an isoparametric boundary may lie slightly outside). Otherwise
        ! Newton is continued from the axis along the straight segment.
        type(nested_map_t), intent(in) :: map
        real(dp), intent(in) :: x(3)
        real(dp), intent(inout) :: q(3)
        integer, intent(out) :: status
        real(dp), parameter :: margin = 0.05_dp
        type(nested_point_t) :: pt
        real(dp) :: target(2), axis(2), guess(2), t, dt
        integer :: steps
        q(3) = atan2(x(2), x(1))
        target = [hypot(x(1), x(2)), x(3)]
        if (q(1)**2 + q(2)**2 <= 1 + margin) then
            guess = q(1:2)
            call newton(target, guess, status)
            if (status == 0 .and. sum(guess**2) <= 1 + margin) then
                q(1:2) = guess
                return
            end if
        end if
        call nested_map_eval(map, 0.0_dp, 0.0_dp, q(3), pt)
        axis = [pt%r, pt%z]
        guess = 0
        t = 0
        dt = 0.125_dp
        steps = 0
        do while (t < 1 .and. steps < 200)
            steps = steps + 1
            q(1:2) = guess
            call newton(axis + min(1.0_dp, t + dt)*(target - axis), q(1:2), status)
            if (status == 0 .and. norm2(q(1:2) - guess) < 0.25_dp) then
                t = min(1.0_dp, t + dt)
                guess = q(1:2)
                dt = min(0.25_dp, 2*dt)
            else
                dt = dt/4
                if (dt < 1.e-6_dp) exit
            end if
        end do
        q(1:2) = guess
        status = merge(0, 1, t >= 1 .and. sum(guess**2) <= 1 + margin)
    contains
        subroutine newton(goal, w, info)
            real(dp), intent(in) :: goal(2)
            real(dp), intent(inout) :: w(2)
            integer, intent(out) :: info
            type(nested_point_t) :: cur, trial
            real(dp) :: res(2), j(2, 2), det, step(2), f0, f1, lam, tol
            integer :: it, ls
            info = 1
            tol = 1.e-14_dp*max(1.0_dp, goal(1))
            call nested_map_eval(map, w(1), w(2), q(3), cur)
            res = [cur%r, cur%z] - goal
            f0 = norm2(res)
            do it = 1, 40
                if (f0 <= tol) then
                    info = 0
                    return
                end if
                j(1, :) = cur%rq(1:2)
                j(2, :) = cur%zq(1:2)
                det = j(1, 1)*j(2, 2) - j(1, 2)*j(2, 1)
                if (abs(det) <= tiny(1.0_dp)) return
                step(1) = (j(2, 2)*res(1) - j(1, 2)*res(2))/det
                step(2) = (-j(2, 1)*res(1) + j(1, 1)*res(2))/det
                lam = 1
                f1 = huge(1.0_dp)
                do ls = 1, 30
                    call nested_map_eval(map, w(1) - lam*step(1), &
                                         w(2) - lam*step(2), q(3), trial)
                    f1 = norm2([trial%r, trial%z] - goal)
                    if (f1 < f0) exit
                    lam = lam/2
                end do
                if (.not. (f1 < f0)) then
                    if (f0 <= 100*tol) info = 0
                    return
                end if
                w = w - lam*step
                cur = trial
                res = [cur%r, cur%z] - goal
                f0 = f1
            end do
        end subroutine
    end subroutine

    pure function nested_field_direction(map, pt, q) result(b)
        ! Unit vector of B ~ (iota - lambda_zeta) e_theta + (1 + lambda_theta)
        ! e_zeta, tangent to the nested surfaces by construction.
        type(nested_map_t), intent(in) :: map
        type(nested_point_t), intent(in) :: pt
        real(dp), intent(in) :: q(3)
        real(dp) :: b(3), etheta(3), iota, ltheta, s
        integer :: k
        s = q(1)**2 + q(2)**2
        iota = 0
        do k = map%nk, 1, -1
            iota = iota*s + map%iota(k)
        end do
        etheta = -q(2)*pt%a(:, 1) + q(1)*pt%a(:, 2)
        ltheta = -q(2)*pt%lam(1) + q(1)*pt%lam(2)
        b = (iota - pt%lam(3))*etheta + (1 + ltheta)*pt%a(:, 3)
        b = b/norm2(b)
    end function

    pure subroutine nested_forcing(pt, q, gradient, laplacian)
        ! Physical gradient and Laplacian of s = u^2 + v^2 from the inverse
        ! function rule H_x s = A^-T (H_q s - sum_i (grad_x s)_i H_q x_i) A^-1.
        type(nested_point_t), intent(in) :: pt
        real(dp), intent(in) :: q(3)
        real(dp), intent(out) :: gradient(3), laplacian
        real(dp) :: ainv(3, 3), hq(3, 3), gq(3)
        integer :: i
        call invert3(pt%a, ainv)
        gq = [2*q(1), 2*q(2), 0.0_dp]
        gradient = matmul(transpose(ainv), gq)
        hq = 0
        hq(1, 1) = 2; hq(2, 2) = 2
        do i = 1, 3
            hq = hq - gradient(i)*pt%h(:, :, i)
        end do
        laplacian = 0
        do i = 1, 3
            laplacian = laplacian + dot_product(ainv(:, i), matmul(hq, ainv(:, i)))
        end do
    end subroutine

    pure subroutine invert3(a, b)
        real(dp), intent(in) :: a(3, 3)
        real(dp), intent(out) :: b(3, 3)
        real(dp) :: det
        b(1, 1) = a(2, 2)*a(3, 3) - a(2, 3)*a(3, 2)
        b(1, 2) = a(1, 3)*a(3, 2) - a(1, 2)*a(3, 3)
        b(1, 3) = a(1, 2)*a(2, 3) - a(1, 3)*a(2, 2)
        b(2, 1) = a(2, 3)*a(3, 1) - a(2, 1)*a(3, 3)
        b(2, 2) = a(1, 1)*a(3, 3) - a(1, 3)*a(3, 1)
        b(2, 3) = a(1, 3)*a(2, 1) - a(1, 1)*a(2, 3)
        b(3, 1) = a(2, 1)*a(3, 2) - a(2, 2)*a(3, 1)
        b(3, 2) = a(1, 2)*a(3, 1) - a(1, 1)*a(3, 2)
        b(3, 3) = a(1, 1)*a(2, 2) - a(1, 2)*a(2, 1)
        det = a(1, 1)*b(1, 1) + a(1, 2)*b(2, 1) + a(1, 3)*b(3, 1)
        b = b/det
    end subroutine

    subroutine build_ogrid_chart(map, kind, n, nphi, nodes, rho_core, vertices, &
                                 geometry, outer, status)
        ! Five-block O-grid x zeta over one field period, with the MF1 block
        ! layout (core square + four rotated ring blocks; reference direction 1
        ! radial, 2 tangent, 3 zeta). Ring blocks follow the curve family
        ! C(rho, theta, zeta): nested surfaces for CHART_ALIGNED, scaled copies
        ! of the boundary about its m=0 centre for CHART_UNALIGNED. For
        ! CHART_FIELDLINE the logical angle is the straight-field-line angle
        ! theta + lambda = theta_logical + iota zeta, so zeta lines are field
        ! lines; plain periodic vertex IDs then require iota/nfp integer. The core is
        ! a Coons patch of the rho_core curve in (R, Z), so the axis is an
        ! ordinary interior point. Vertex IDs are combinatorial and periodic in
        ! zeta, which identifies the field-period end planes. nodes are the
        ! geometry interpolation nodes on [-1, 1] (MF1: GLL of degree q).
        type(nested_map_t), intent(in) :: map
        integer, intent(in) :: kind, n, nphi
        real(dp), intent(in) :: nodes(:), rho_core
        integer, allocatable, intent(out) :: vertices(:, :)
        real(dp), allocatable, intent(out) :: geometry(:, :, :)
        logical, allocatable, intent(out) :: outer(:)
        integer, intent(out) :: status
        real(dp) :: period, ref(3), xi, eta, zeta, rz(2), d(3, 3)
        integer :: nq, ng, ne, e, iz, block, ix, iy, a, c, bits(3), n2, &
                   corner(3), idx(3), k, lo(3)
        status = 1
        if (kind < CHART_ALIGNED .or. kind > CHART_FIELDLINE) return
        if (kind == CHART_FIELDLINE) then
            if (map%nk > 1) then
                if (any(abs(map%iota(2:)) > 1.e-10_dp)) return
            end if
            if (abs(map%iota(1)/map%nfp - nint(map%iota(1)/map%nfp)) > 1.e-10_dp) return
        end if
        if (n < 1 .or. nphi < 1 .or. size(nodes) < 2) return
        if (rho_core <= 0 .or. rho_core >= 1) return
        nq = size(nodes)
        ng = nq**3
        ne = 5*n*n*nphi
        n2 = (n + 1)**2 + 4*n*n
        period = 2*pi/map%nfp
        allocate (vertices(8, ne), geometry(3, ng, ne), outer(ne))
        e = 0
        do iz = 0, nphi - 1
            do block = 0, 4
                do iy = 0, n - 1
                    do ix = 0, n - 1
                        e = e + 1
                        outer(e) = block > 0 .and. ix == n - 1
                        do c = 1, 8
                            bits = [modulo(c - 1, 2), modulo((c - 1)/2, 2), (c - 1)/4]
                            corner = [ix, iy, iz] + bits
                            vertices(c, e) = vertex_id(block, corner(1), corner(2)) &
                                             + n2*modulo(corner(3), nphi)
                        end do
                        do a = 1, ng
                            idx = [modulo(a - 1, nq), modulo((a - 1)/nq, nq), &
                                   (a - 1)/nq**2] + 1
                            ref = (nodes(idx) + 1)/2
                            xi = (ix + ref(1))/n
                            eta = (iy + ref(2))/n
                            zeta = period*(iz + ref(3))/nphi
                            rz = chart_rz(block, xi, eta, zeta)
                            geometry(:, a, e) = [rz(1)*cos(zeta), rz(1)*sin(zeta), &
                                                 rz(2)]
                        end do
                        ! Orientation check at every corner from node differences.
                        do c = 1, 8
                            bits = [modulo(c - 1, 2), modulo((c - 1)/2, 2), (c - 1)/4]
                            lo = merge(nq - 1, 1, bits == 1)
                            do k = 1, 3
                                idx = merge(nq, 1, bits == 1)
                                idx(k) = lo(k) + 1
                                d(:, k) = geometry(:, flat(idx), e)
                                idx(k) = lo(k)
                                d(:, k) = d(:, k) - geometry(:, flat(idx), e)
                            end do
                            if (det3(d) <= 0) then
                                status = 5
                                return
                            end if
                        end do
                    end do
                end do
            end do
        end do
        status = 0
    contains
        pure integer function flat(j)
            integer, intent(in) :: j(3)
            flat = j(1) + nq*(j(2) - 1) + nq*nq*(j(3) - 1)
        end function

        pure integer function vertex_id(blk, i, j)
            ! Logical 2D vertex: core grid, then ring blocks without their
            ! inner (core) edge and without the edge shared with the next block.
            integer, intent(in) :: blk, i, j
            integer :: ci, cj, b, jj
            if (blk == 0 .or. i == 0) then
                ci = i; cj = j
                if (blk == 1) then
                    ci = n; cj = j
                else if (blk == 2) then
                    ci = n - j; cj = n
                else if (blk == 3) then
                    ci = 0; cj = n - j
                else if (blk == 4) then
                    ci = j; cj = 0
                end if
                vertex_id = 1 + ci + (n + 1)*cj
                return
            end if
            b = blk
            jj = j
            if (jj == n) then
                b = modulo(blk, 4) + 1
                jj = 0
            end if
            vertex_id = (n + 1)**2 + 1 + ((b - 1)*n + (i - 1))*n + jj
        end function

        function chart_rz(blk, xr, xt, zt) result(rz)
            integer, intent(in) :: blk
            real(dp), intent(in) :: xr, xt, zt
            real(dp) :: rz(2), p00(2), p10(2), p01(2), p11(2)
            if (blk > 0) then
                rz = curve(rho_core + (1 - rho_core)*xr, &
                           (blk - 1)*pi/2 + pi*(2*xt - 1)/4, zt)
                return
            end if
            ! Core edges carry the ring blocks' inner curves (see vertex_id).
            p00 = curve(rho_core, 5*pi/4, zt)
            p10 = curve(rho_core, 7*pi/4, zt)
            p01 = curve(rho_core, 3*pi/4, zt)
            p11 = curve(rho_core, pi/4, zt)
            rz = (1 - xr)*curve(rho_core, pi - pi*(2*xt - 1)/4, zt) &
                 + xr*curve(rho_core, pi*(2*xt - 1)/4, zt) &
                 + (1 - xt)*curve(rho_core, 3*pi/2 + pi*(2*xr - 1)/4, zt) &
                 + xt*curve(rho_core, pi/2 - pi*(2*xr - 1)/4, zt) &
                 - (1 - xr)*(1 - xt)*p00 - xr*(1 - xt)*p10 &
                 - (1 - xr)*xt*p01 - xr*xt*p11
        end function

        function curve(rho, theta_logical, zt) result(rz)
            real(dp), intent(in) :: rho, theta_logical, zt
            real(dp) :: rz(2), theta, centre(2)
            type(nested_point_t) :: pt
            integer :: j
            theta = map%sign_theta*theta_logical
            if (kind == CHART_FIELDLINE) theta = sfl_theta(rho, theta + map%iota(1)*zt, zt)
            if (kind /= CHART_UNALIGNED) then
                call nested_map_eval(map, rho*cos(theta), rho*sin(theta), zt, pt)
                rz = [pt%r, pt%z]
                return
            end if
            ! Boundary-only chart: the centre is the boundary's m=0 part, i.e.
            ! its exact angular mean at s = 1.
            centre = 0
            do j = 1, map%nmode
                if (map%m(j) /= 0) cycle
                centre = centre + [sum(map%rc(:, j))*cos(map%n(j)*zt), &
                                   -sum(map%zs(:, j))*sin(map%n(j)*zt)]
            end do
            call nested_map_eval(map, cos(theta), sin(theta), zt, pt)
            rz = centre + rho*([pt%r, pt%z] - centre)
        end function

        function sfl_theta(rho, target, zt) result(theta)
            ! Solve theta + lambda(rho, theta, zeta) = target by Newton.
            real(dp), intent(in) :: rho, target, zt
            real(dp) :: theta, f, df
            type(nested_point_t) :: pt
            integer :: it
            theta = target
            do it = 1, 50
                call nested_map_eval(map, rho*cos(theta), rho*sin(theta), zt, pt)
                ! The point jets hold lambda derivatives only, so sum lambda here.
                f = theta + lambda_value(rho, theta, zt) - target
                df = 1 - rho*sin(theta)*pt%lam(1) + rho*cos(theta)*pt%lam(2)
                theta = theta - f/df
                if (abs(f) < 1.e-15_dp) exit
            end do
        end function

        real(dp) function lambda_value(rho, theta, zt)
            real(dp), intent(in) :: rho, theta, zt
            integer :: j
            lambda_value = 0
            do j = 1, map%nmode
                lambda_value = lambda_value + rho**map%m(j) &
                    *poly(map%ls(:, j), rho*rho)*sin(map%m(j)*theta - map%n(j)*zt)
            end do
        end function

        real(dp) function poly(c, x)
            real(dp), intent(in) :: c(:), x
            integer :: k
            poly = 0
            do k = size(c), 1, -1
                poly = poly*x + c(k)
            end do
        end function
    end subroutine

    pure real(dp) function det3(a)
        real(dp), intent(in) :: a(3, 3)
        det3 = a(1, 1)*(a(2, 2)*a(3, 3) - a(2, 3)*a(3, 2)) &
               - a(1, 2)*(a(2, 1)*a(3, 3) - a(2, 3)*a(3, 1)) &
               + a(1, 3)*(a(2, 1)*a(3, 2) - a(2, 2)*a(3, 1))
    end function
end module
