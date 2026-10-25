program test_ideal_fb_surface
    !! Public synthetic chart/pullback oracles, without a BEM assembly.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use kin6d_equilibrium_chart, only: cylindrical_chart_t, curve_chart_t
    use kin6d_ideal_free_boundary, only: ideal_fb_t
    use kin6d_ideal_free_boundary_bem3d, only: coupling_surface_t, ideal_fb_surface
    use fortfem_parametric_bem_high_order, only: build_p2_trace_connectivity
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp), radius = .15_dp
    real(dp), parameter :: nodes(2, 6) = reshape( &
        [0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, &
         .5_dp, 0.0_dp, .5_dp, .5_dp, 0.0_dp, .5_dp], [2, 6])
    character(32) :: mode, arg
    integer :: order, family
    call get_command_argument(1, mode)
    call get_command_argument(2, arg)
    if (trim(mode) == 'invalid') then
        call invalid()
    else
        read (arg, *) order
        select case (trim(mode))
        case ('cylinder'); family = 0
        case ('circle'); family = 1
        case ('helix'); family = 2
        case default; error stop 'unknown surface test'
        end select
        call geometry(family, order)
    end if
contains
    subroutine require(ok, label)
        logical, intent(in) :: ok
        character(*), intent(in) :: label
        if (ok) return
        write (*, '(a)') label
        error stop 'chart surface gate'
    end subroutine

    subroutine fixture(fb, family, order)
        type(ideal_fb_t), intent(out) :: fb
        integer, intent(in) :: family, order
        type(curve_chart_t) :: chart
        real(dp) :: c(3, 0:4), s(3, 0:4), theta
        integer :: e, k, i, status
        fb%order = order; fb%nedge = 8; fb%nt = order*fb%nedge; fb%ns = fb%nt
        fb%nfp = 1
        if (family == 2) fb%nfp = 3
        fb%ntor = 1; fb%nmode = 3
        allocate (fb%reference(2, fb%ns), fb%trace(fb%nt), fb%edge(0:order, fb%nedge))
        do e = 1, fb%nedge
            do k = 0, order
                i = mod((e - 1)*order + k, fb%ns) + 1
                theta = 2*pi*real((e - 1)*order + k, dp)/fb%ns
                fb%reference(:, i) = radius*[cos(theta), sin(theta)]
                fb%edge(k, e) = i
            end do
        end do
        fb%trace = [(i, i=1, fb%nt)]
        if (family == 0) then
            fb%reference(1, :) = fb%reference(1, :) + 3
            allocate (cylindrical_chart_t :: fb%plasma%chart)
            fb%plasma%sigma = -1
        else
            c = 0; s = 0; c(1, 1) = 3; s(2, 1) = 3
            if (family == 2) then
                c(1, 2) = .1_dp; s(2, 2) = -.1_dp
                c(1, 4) = .1_dp; s(2, 4) = .1_dp; s(3, 3) = .2_dp
            end if
            call chart%init(c, s, fb%nfp, status)
            call require(status == 0, 'existing public generated chart initializes')
            allocate (fb%plasma%chart, source=chart)
            fb%plasma%sigma = 1
        end if
    end subroutine

    subroutine edge_point(fb, e, t, q, tangent)
        !! Independent Lagrange products, rather than production triangle leaves.
        type(ideal_fb_t), intent(in) :: fb
        integer, intent(in) :: e
        real(dp), intent(in) :: t
        real(dp), intent(out) :: q(2), tangent(2)
        real(dp) :: w, dw, term, tk, tl
        integer :: k, l, j
        q = 0; tangent = 0
        do k = 0, fb%order
            tk = real(k, dp)/fb%order; w = 1; dw = 0
            do l = 0, fb%order
                if (l == k) cycle
                tl = real(l, dp)/fb%order
                w = w*(t - tl)/(tk - tl)
                term = 1/(tk - tl)
                do j = 0, fb%order
                    if (j == k .or. j == l) cycle
                    term = term*(t - real(j, dp)/fb%order)/(tk - real(j, dp)/fb%order)
                end do
                dw = dw + term
            end do
            q = q + w*fb%reference(:, fb%edge(k, e))
            tangent = tangent + dw*fb%reference(:, fb%edge(k, e))
        end do
    end subroutine

    function cross(a, b) result(c)
        real(dp), intent(in) :: a(3), b(3)
        real(dp) :: c(3)
        c = [a(2)*b(3) - a(3)*b(2), a(3)*b(1) - a(1)*b(3), a(1)*b(2) - a(2)*b(1)]
    end function

    function inverse(a) result(b)
        real(dp), intent(in) :: a(3, 3)
        real(dp) :: b(3, 3), det
        b(1, :) = cross(a(:, 2), a(:, 3))
        b(2, :) = cross(a(:, 3), a(:, 1))
        b(3, :) = cross(a(:, 1), a(:, 2))
        det = dot_product(a(:, 1), b(1, :)); b = b/det
    end function

    subroutine topology(faces)
        integer, intent(in) :: faces(:, :)
        integer, allocatable :: count(:, :), orient(:, :)
        integer :: p, i, a, b, nv, ne
        nv = maxval(faces); allocate (count(nv, nv), orient(nv, nv)); count = 0; orient = 0
        do p = 1, size(faces, 2)
            do i = 1, 3
                a = faces(i, p); b = faces(mod(i, 3) + 1, p)
                count(min(a, b), max(a, b)) = count(min(a, b), max(a, b)) + 1
                orient(min(a, b), max(a, b)) = orient(min(a, b), max(a, b)) + merge(1, -1, a < b)
            end do
        end do
        ne = count_nonzero(count)
        call require(all(count == 0 .or. count == 2), 'every closed edge has two faces')
        call require(all(orient == 0), 'paired faces traverse each edge oppositely')
        call require(nv - ne + size(faces, 2) == 0, 'closed periodic genus-one topology')
    end subroutine

    integer function count_nonzero(values)
        integer, intent(in) :: values(:, :)
        count_nonzero = count(values /= 0)
    end function

    subroutine geometry(family, order)
        integer, intent(in) :: family, order
        type(ideal_fb_t) :: fb
        type(coupling_surface_t) :: surface
        integer, allocatable :: faces(:, :), td(:, :), fd(:, :)
        real(dp), allocatable :: restriction(:, :), coefficients(:), trace(:)
        real(dp) :: uv(2), parameter(2), q(3), tangent(2), x(3), a(3), b(3)
        real(dp) :: d1(3, 3), d2(3, 3, 3), d3(3, 3, 3, 3), ji(3, 3)
        real(dp) :: normal(3), area(3), refnormal(3), field(3), bref(3), det, parameter_det
        real(dp) :: xp(3), xm(3), ap(3), bp(3), fd1(3), fd2(3), h, expected, error_trace
        real(dp) :: error_flux, error_tangent, error_period, outward, period, t0, phi0
        integer :: refine, p, i, e, status, nphi
        call fixture(fb, family, order)
        nphi = 16; bref = [.7_dp, -.2_dp, .3_dp]; h = 1e-5_dp
        allocate (coefficients(fb%nt*fb%nmode))
        coefficients(:fb%nt) = fb%reference(1, :)
        coefficients(fb%nt + 1:2*fb%nt) = .3_dp*fb%reference(2, :)
        coefficients(2*fb%nt + 1:) = .2_dp*fb%reference(1, :)
        error_trace = 0; error_flux = 0; error_tangent = 0; error_period = 0; outward = 1
        do refine = 1, 3
            call ideal_fb_surface(fb, nphi, surface, faces, restriction, status, refine)
            call require(status == 0, 'factorized chart surface constructor')
            call require(size(faces, 2) == 2*fb%nedge*nphi*refine, 'independent trace subdivision')
            call topology(faces)
            if (family == 0 .and. refine == 1) then
                call require(all(faces(:, 1) == [1, 9, 2]), 'default cylindrical first face unchanged')
                call require(all(faces(:, 2) == [9, 10, 2]), 'default cylindrical second face unchanged')
            end if
            call build_p2_trace_connectivity(faces, td, fd, status)
            call require(status == 0, 'continuous P2 trace connectivity')
            trace = matmul(restriction, coefficients)
            do p = 1, size(faces, 2)
                e = surface%panel_edge(p)
                do i = 1, 6
                    uv = nodes(:, i)
                    parameter = surface%parameters(:, 1, p) &
                        + uv(1)*(surface%parameters(:, 2, p) - surface%parameters(:, 1, p)) &
                        + uv(2)*(surface%parameters(:, 3, p) - surface%parameters(:, 1, p))
                    call edge_point(fb, e, parameter(1), q(:2), tangent)
                    q(3) = parameter(2)
                    expected = q(1) + .3_dp*q(2)*cos(fb%nfp*q(3)) + .2_dp*q(1)*sin(fb%nfp*q(3))
                    error_trace = max(error_trace, abs(trace(td(i, p)) - expected))
                end do
                uv = [.23_dp, .31_dp]
                parameter = surface%parameters(:, 1, p) &
                    + uv(1)*(surface%parameters(:, 2, p) - surface%parameters(:, 1, p)) &
                    + uv(2)*(surface%parameters(:, 3, p) - surface%parameters(:, 1, p))
                call edge_point(fb, e, parameter(1), q(:2), tangent); q(3) = parameter(2)
                call surface%evaluate(p, uv, x, a, b, status)
                call require(status == 0, 'regular physical surface point')
                call fb%plasma%chart%jet(q, xp, d1, d2, d3)
                error_tangent = max(error_tangent, maxval(abs(x - xp)))
                ji = inverse(d1)
                det = dot_product(d1(:, 1), cross(d1(:, 2), d1(:, 3)))
                refnormal = [tangent(2), -tangent(1), 0.0_dp]
                normal = matmul(transpose(ji), refnormal); area = cross(a, b)
                outward = min(outward, dot_product(area, normal)/(norm2(area)*norm2(normal)))
                parameter_det = abs((surface%parameters(1, 2, p) - surface%parameters(1, 1, p))* &
                    (surface%parameters(2, 3, p) - surface%parameters(2, 1, p)) - &
                    (surface%parameters(1, 3, p) - surface%parameters(1, 1, p))* &
                    (surface%parameters(2, 2, p) - surface%parameters(2, 1, p)))
                field = matmul(d1, bref)/det
                expected = fb%plasma%sigma*dot_product(bref, refnormal)*parameter_det
                error_flux = max(error_flux, abs(dot_product(field, area) - expected))
                call surface%evaluate(p, uv + [h, 0.0_dp], xp, ap, bp, status)
                call surface%evaluate(p, uv - [h, 0.0_dp], xm, ap, bp, status)
                fd1 = (xp - xm)/(2*h)
                call surface%evaluate(p, uv + [0.0_dp, h], xp, ap, bp, status)
                call surface%evaluate(p, uv - [0.0_dp, h], xm, ap, bp, status)
                fd2 = (xp - xm)/(2*h)
                error_tangent = max(error_tangent, norm2(a - fd1), norm2(b - fd2))
            end do
            ! The pulled-back one-form dphi has a full-turn period 2*pi.
            period = 0; t0 = .29_dp
            call edge_point(fb, 1, t0, q(:2), tangent)
            do i = 1, 32
                phi0 = 2*pi*(real(i, dp) - .5_dp)/32; q(3) = phi0
                call fb%plasma%chart%jet(q, x, d1, d2, d3); ji = inverse(d1)
                q(3) = phi0 + h; call fb%plasma%chart%jet(q, xp, d1, d2, d3)
                q(3) = phi0 - h; call fb%plasma%chart%jet(q, xm, d1, d2, d3)
                period = period + dot_product(ji(3, :), (xp - xm)/(2*h))*2*pi/32
            end do
            error_period = max(error_period, abs(period - 2*pi))
            q(3) = 0; call fb%plasma%chart%jet(q, xp, d1, d2, d3)
            q(3) = 2*pi; call fb%plasma%chart%jet(q, xm, d1, d2, d3)
            call require(norm2(xp - xm) < 1e-13_dp, 'full-turn physical seam closes')
        end do
        write (*, '(a,4es15.7)') 'trace/flux/tangent/period errors ', &
            error_trace, error_flux, error_tangent, error_period
        write (*, '(a,es15.7)') 'minimum outward normal cosine ', outward
        call require(outward > 1 - 1e-12_dp, 'physical normals are outward for both map orientations')
        call require(error_trace < 2e-13_dp, 'independent polynomial/Fourier trace pullback')
        call require(error_flux < 2e-13_dp, 'signed contravariant Piola normal flux')
        call require(error_tangent < 2e-9_dp, 'independent finite-difference tangents')
        call require(error_period < 2e-9_dp, 'one-form full-turn period without nfp rescaling')
    end subroutine

    subroutine invalid()
        type(ideal_fb_t) :: fb
        type(coupling_surface_t) :: surface
        integer, allocatable :: faces(:, :)
        real(dp), allocatable :: restriction(:, :)
        integer :: status
        call fixture(fb, 1, 2)
        call ideal_fb_surface(fb, 16, surface, faces, restriction, status, 0)
        call require(status /= 0, 'reject zero subdivision')
        call ideal_fb_surface(fb, 2, surface, faces, restriction, status)
        call require(status /= 0, 'reject unresolved azimuthal trace')
        fb%edge(2, 1) = fb%edge(0, 1)
        call ideal_fb_surface(fb, 16, surface, faces, restriction, status)
        call require(status /= 0, 'reject a disconnected reference chain')
        call fixture(fb, 1, 2); fb%trace(2) = fb%trace(1)
        call ideal_fb_surface(fb, 16, surface, faces, restriction, status)
        call require(status /= 0, 'reject repeated trace node')
        call fixture(fb, 1, 2); fb%reference = 0
        call ideal_fb_surface(fb, 16, surface, faces, restriction, status)
        call require(status /= 0, 'reject degenerate surface tangent')
        call fixture(fb, 1, 2); fb%reference(1, 1) = ieee_value(0.0_dp, ieee_quiet_nan)
        call ideal_fb_surface(fb, 16, surface, faces, restriction, status)
        call require(status /= 0, 'reject nonfinite reference geometry')
    end subroutine
end program
