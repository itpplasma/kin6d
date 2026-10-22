program test_nested_chart
    ! Independent checks of the nested map jets, inverse, field direction,
    ! manufactured Laplacian and both O-grid charts. A synthetic shifted,
    ! rotating-ellipse map is built here; an optional argument adds a fitted
    ! wout map file.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_nested_chart, only: nested_map_t, nested_point_t, nested_map_read, &
        nested_map_eval, nested_map_invert, nested_field_direction, nested_forcing, &
        build_ogrid_chart, CHART_ALIGNED, CHART_UNALIGNED, CHART_FIELDLINE
    implicit none
    type(nested_map_t) :: map
    character(512) :: path
    integer :: failures, status
    failures = 0
    call synthetic(map)
    call check_map(map, 'synthetic')
    if (command_argument_count() >= 1) then
        call get_command_argument(1, path)
        call nested_map_read(trim(path), map, status)
        call expect(status == 0, 'read '//trim(path))
        if (status == 0) call check_map(map, 'file')
    end if
    if (failures /= 0) then
        print '(a,i0)', 'FAILED checks: ', failures
        error stop 1
    end if
    print '(a)', 'test_nested_chart passed'
contains
    subroutine expect(ok, label)
        logical, intent(in) :: ok
        character(*), intent(in) :: label
        if (.not. ok) then
            failures = failures + 1
            print '(2a)', 'FAIL: ', label
        end if
    end subroutine

    subroutine synthetic(map)
        type(nested_map_t), intent(out) :: map
        character(*), parameter :: file = 'synthetic_nested_map.txt'
        integer :: unit, status
        ! R = 1 + 0.03 s + 0.2 rho cos(t) + 0.02 rho cos(t - 2 z) + 0.01 s rho^2 cos(2t)
        ! Z = 0.25 rho sin(t) + 0.02 rho sin(t - 2 z) - 0.04 sin(-2 z), lambda small.
        open (newunit=unit, file=file, status='replace', action='write')
        write (unit, '(a)') '# kin6d nested map v1 synthetic'
        write (unit, '(a)') '2 4 2'
        write (unit, '(a)') '0.4 0.1'
        write (unit, '(a)') '0 0 1.0 0.03 0.0 0.0 0.0 0.0'
        write (unit, '(a)') '0 2 0.0 0.0 0.04 0.0 0.0 0.0'
        write (unit, '(a)') '1 0 0.2 0.0 0.25 0.0 0.01 0.0'
        write (unit, '(a)') '1 2 0.02 0.0 0.02 0.0 0.0 0.005'
        close (unit)
        call nested_map_read(file, map, status)
        call expect(status == 0, 'synthetic read')
    end subroutine

    subroutine check_map(map, label)
        type(nested_map_t), intent(in) :: map
        character(*), intent(in) :: label
        type(nested_point_t) :: pt, pp, pm
        real(dp) :: q(3), qs(3), dq(3), err1, err2, b(3), grad(3), lap, x(3), &
                    h, lapfd, s0, sp, sm, gradfd(3), scale
        integer :: i, k, trial, status
        real(dp), parameter :: samples(3, 4) = reshape([ &
            0.3_dp, -0.2_dp, 0.4_dp, 0.0_dp, 0.05_dp, 1.3_dp, &
            -0.6_dp, 0.5_dp, 2.0_dp, 0.1_dp, -0.85_dp, 2.9_dp], [3, 4])
        err1 = 0; err2 = 0
        h = 1.e-5_dp
        do trial = 1, 4
            q = samples(:, trial)
            call nested_map_eval(map, q(1), q(2), q(3), pt)
            scale = maxval(abs(pt%x))
            do k = 1, 3
                dq = 0; dq(k) = h
                call nested_map_eval(map, q(1) + dq(1), q(2) + dq(2), q(3) + dq(3), pp)
                call nested_map_eval(map, q(1) - dq(1), q(2) - dq(2), q(3) - dq(3), pm)
                err1 = max(err1, maxval(abs((pp%x - pm%x)/(2*h) - pt%a(:, k)))/scale)
                do i = 1, 3
                    err2 = max(err2, maxval(abs((pp%a(i, :) - pm%a(i, :))/(2*h) &
                                                - pt%h(:, k, i)))/scale)
                end do
            end do
            ! Inverse round trip.
            qs = [0.0_dp, 0.0_dp, 0.0_dp]
            call nested_map_invert(map, pt%x, qs, status)
            call expect(status == 0, label//' invert status')
            call expect(norm2(qs(1:2) - q(1:2)) < 1.e-11_dp .and. &
                        abs(modulo(qs(3) - q(3) + 3.1_dp, 6.2831853071795864_dp) &
                            - 3.1_dp) < 1.e-12_dp, label//' invert round trip')
            ! Field direction is tangent to the surfaces and the manufactured
            ! Laplacian matches second differences of s(x).
            b = nested_field_direction(map, pt, q)
            call nested_forcing(pt, q, grad, lap)
            call expect(abs(dot_product(b, grad)) < 1.e-12_dp*norm2(grad) + 1.e-14_dp, &
                        label//' b . grad s = 0')
            s0 = q(1)**2 + q(2)**2
            lapfd = 0
            h = 2.e-4_dp*scale
            do k = 1, 3
                x = pt%x; x(k) = x(k) + h
                qs = q
                call nested_map_invert(map, x, qs, status)
                sp = qs(1)**2 + qs(2)**2
                x = pt%x; x(k) = x(k) - h
                qs = q
                call nested_map_invert(map, x, qs, status)
                sm = qs(1)**2 + qs(2)**2
                gradfd(k) = (sp - sm)/(2*h)
                lapfd = lapfd + (sp - 2*s0 + sm)/h**2
            end do
            ! Central differences of an iterative inverse: O(h^2) ~ 1e-5 near
            ! the axis, where s varies on the scale of rho.
            call expect(norm2(gradfd - grad) < 1.e-4_dp*norm2(grad), &
                        label//' grad s vs differences')
            call expect(abs(lapfd - lap) < 1.e-4_dp*max(1.0_dp, abs(lap)), &
                        label//' laplacian s vs differences')
            if (abs(lapfd - lap) >= 1.e-4_dp*max(1.0_dp, abs(lap))) &
                print '(a,2es14.6)', '  laplacian analytic/fd ', lap, lapfd
            h = 1.e-5_dp
        end do
        call expect(err1 < 1.e-8_dp, label//' first-derivative jets')
        call expect(err2 < 1.e-7_dp, label//' second-derivative jets')
        print '(a,a,2es10.2)', label, ' jet errors ', err1, err2
        call check_chart(map, CHART_ALIGNED, label//' aligned')
        call check_chart(map, CHART_UNALIGNED, label//' unaligned')
        call check_fieldline(map, label)
    end subroutine

    subroutine check_chart(map, kind, label)
        type(nested_map_t), intent(in) :: map
        integer, intent(in) :: kind
        character(*), intent(in) :: label
        integer, allocatable :: vertices(:, :)
        real(dp), allocatable :: geometry(:, :, :), seen(:, :)
        logical, allocatable :: visited(:)
        logical, allocatable :: outer(:)
        real(dp) :: nodes(4), mismatch, q(3), surface
        integer :: status, e, c, a, v, nq, corner_node(8), status_inv, i, j, k
        real(dp), parameter :: pi = acos(-1.0_dp)
        nq = size(nodes)
        nodes = [(-cos(pi*i/(nq - 1)), i=0, nq - 1)]
        call build_ogrid_chart(map, kind, 3, 3, nodes, 0.35_dp, vertices, &
                               geometry, outer, status)
        call expect(status == 0, label//' build (positive corner Jacobians)')
        if (status /= 0) return
        call expect(size(vertices, 2) == 5*9*3, label//' element count')
        call expect(count(outer) == 4*3*3, label//' outer element count')
        call expect(maxval(vertices) == 3*((3 + 1)**2 + 4*9), label//' vertex count')
        do c = 1, 8
            i = modulo(c - 1, 2); j = modulo((c - 1)/2, 2); k = (c - 1)/4
            corner_node(c) = 1 + i*(nq - 1) + nq*j*(nq - 1) + nq*nq*k*(nq - 1)
        end do
        ! Shared vertices must carry identical geometry (conformity), except the
        ! periodic end planes, which differ by the field-period rotation.
        allocate (seen(3, maxval(vertices)), visited(maxval(vertices)))
        visited = .false.
        mismatch = 0
        do e = 1, size(vertices, 2)
            do c = 1, 8
                v = vertices(c, e)
                if (c > 4 .and. e > size(vertices, 2)*2/3) cycle
                if (.not. visited(v)) then
                    visited(v) = .true.
                    seen(:, v) = geometry(:, corner_node(c), e)
                else
                    mismatch = max(mismatch, norm2(seen(:, v) - geometry(:, corner_node(c), e)))
                end if
            end do
        end do
        call expect(mismatch < 1.e-12_dp, label//' shared vertex geometry')
        ! Aligned ring nodes lie on their nested surface s = rho^2.
        if (kind == CHART_ALIGNED) then
            surface = 0
            do e = 1, size(vertices, 2)
                if (.not. outer(e)) cycle
                do a = 1, nq**3
                    if (modulo(a - 1, nq) /= nq - 1) cycle
                    q = 0
                    call nested_map_invert(map, geometry(:, a, e), q, status_inv)
                    surface = max(surface, abs(q(1)**2 + q(2)**2 - 1))
                    if (abs(q(1)**2 + q(2)**2 - 1) > 1.e-10_dp .or. status_inv /= 0) &
                        print '(a,i0,i3,es12.4)', '  boundary node ', e, status_inv, &
                        q(1)**2 + q(2)**2 - 1
                end do
            end do
            call expect(surface < 1.e-10_dp, label//' boundary nodes on s = 1')
        end if
    end subroutine

    subroutine check_fieldline(map, label)
        ! Zeta lines of the field-line chart keep theta + lambda - iota zeta
        ! fixed; a non-integer iota/nfp must be rejected.
        type(nested_map_t), intent(in) :: map
        character(*), intent(in) :: label
        integer, allocatable :: vertices(:, :)
        real(dp), allocatable :: geometry(:, :, :)
        logical, allocatable :: outer(:)
        real(dp) :: nodes(4), q(3), angle(4), drift, lam, rho, theta
        integer :: status, e, i, j, t, a, nq
        real(dp), parameter :: pi = acos(-1.0_dp)
        nq = size(nodes)
        nodes = [(-cos(pi*i/(nq - 1)), i=0, nq - 1)]
        call build_ogrid_chart(map, CHART_FIELDLINE, 3, 6, nodes, 0.35_dp, &
                               vertices, geometry, outer, status)
        if (abs(map%iota(1)/map%nfp - nint(map%iota(1)/map%nfp)) > 1.e-10_dp) then
            call expect(status /= 0, label//' fieldline rejects non-integer iota/nfp')
            return
        end if
        call expect(status == 0, label//' fieldline build')
        if (status /= 0) return
        drift = 0
        do e = 1, size(vertices, 2)
            if (.not. outer(e)) cycle
            do j = 1, nq
                do i = 2, nq
                    do t = 1, nq
                        a = i + nq*(j - 1) + nq*nq*(t - 1)
                        q = 0
                        call nested_map_invert(map, geometry(:, a, e), q, status)
                        rho = hypot(q(1), q(2))
                        theta = atan2(q(2), q(1))
                        lam = lambda_at(map, rho, theta, q(3))
                        angle(t) = theta + lam - map%iota(1)*q(3)
                    end do
                    angle = modulo(angle - angle(1) + pi, 2*pi) - pi
                    drift = max(drift, maxval(abs(angle)))
                end do
            end do
        end do
        call expect(drift < 1.e-9_dp, label//' fieldline zeta lines follow b')
        print '(a,a,es10.2)', label, ' fieldline drift ', drift
    end subroutine

    real(dp) function lambda_at(map, rho, theta, zeta)
        type(nested_map_t), intent(in) :: map
        real(dp), intent(in) :: rho, theta, zeta
        integer :: j, k
        lambda_at = 0
        do j = 1, map%nmode
            do k = 1, map%nk
                lambda_at = lambda_at + map%ls(k, j)*rho**(map%m(j) + 2*(k - 1)) &
                    *sin(map%m(j)*theta - map%n(j)*zeta)
            end do
        end do
    end function
end program
