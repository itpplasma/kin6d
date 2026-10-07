program test_gs_p1_weak_current
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_quadrature, only: gauss_legendre_ab
    use kin6d_grad_shafranov_leaves, only: gs_field
    implicit none
    integer, parameter :: nw = 5
    real(dp), parameter :: xy(2, 4) = reshape([1._dp, 0._dp, 2._dp, 0._dp, &
        2._dp, 1._dp, 1._dp, 1._dp], [2, 4])
    real(dp) :: values(4), cell(nw), jump(nw), weak(nw), boundary(nw), lhs(nw)
    real(dp) :: expected(4, nw), l2
    integer :: triangles(3, 2), family, ordering, i
    character(32) :: mutation

    call get_command_argument(1, mutation)
    l2 = log(2._dp)
    do family = 1, 2
        if (family == 1) then
            values = xy(1, :) + 2*xy(2, :)
            expected(1, :) = [0.5_dp, l2, 0.25_dp, l2/2, (3*l2 - 2)/6]
            expected(2, :) = 0
            expected(3, :) = [0._dp, l2, 2*l2, 2 + l2/2, (3*l2 - 2)/6]
            expected(4, :) = [-0.5_dp, 0._dp, 2*l2 - 0.25_dp, 2._dp, 0._dp]
        else
            values = xy(1, :)*xy(2, :)
            expected(1, :) = [1 - l2, 2*l2 - 1, l2 - 0.5_dp, 0.25_dp, &
                (25 - 36*l2)/6]
            expected(2, :) = [2*l2, 2._dp, 2 - 2*l2, 1._dp, 8*l2 - 5.5_dp]
            expected(3, :) = [0._dp, 2*l2 - 1, 1._dp, 1.75_dp, 2*l2 - 4._dp/3]
            expected(4, :) = [-1 - l2, -2._dp, l2 - 0.5_dp, 0.5_dp, 0._dp]
        end if
        do ordering = 1, 4
            select case (ordering)
            case (1, 4)
                triangles = reshape([1, 2, 3, 1, 3, 4], [3, 2])
            case (2)
                triangles = reshape([3, 2, 1, 4, 3, 1], [3, 2])
            case (3)
                triangles = reshape([4, 1, 3, 3, 1, 2], [3, 2])
            end select
            call integrate(triangles, values, ordering == 4, cell, jump, weak, boundary)
            do i = 1, nw
                call require_close(cell(i), expected(1, i), 'cell-current oracle')
                call require_close(jump(i), expected(2, i), 'edge-current oracle')
                call require_close(weak(i), expected(3, i), 'weak-volume oracle')
                call require_close(boundary(i), expected(4, i), 'boundary-flux oracle')
            end do
            lhs = cell + jump
            select case (trim(mutation))
            case ('omit-jump')
                lhs = cell
            case ('omit-cylindrical')
                lhs = jump
            case ('wrong-current-sign')
                lhs = -lhs
            case ('')
            case default
                error stop 'Unknown current mutation'
            end select
            do i = 1, nw
                call require_close(lhs(i), weak(i) - boundary(i), &
                    'distributional current must include cells, jumps and boundary')
            end do
        end do
    end do
    print *, 'PASS: native signed field and P1 cell/edge/boundary reorder oracles'
contains
    subroutine require_close(actual, oracle, label)
        real(dp), intent(in) :: actual, oracle
        character(*), intent(in) :: label
        if (abs(actual - oracle) <= 3.e-12_dp) return
        print *, label, ': actual/oracle/difference', actual, oracle, actual - oracle
        error stop 1
    end subroutine require_close

    subroutine weights(r, z, v, vr, vz)
        real(dp), intent(in) :: r, z
        real(dp), intent(out) :: v(nw), vr(nw), vz(nw)
        ! Independent test polynomials; the fifth vanishes on the whole boundary.
        v = [1._dp, r, z, r*z, (r - 1)*(2 - r)*z*(1 - z)]
        vr = [0._dp, 1._dp, 0._dp, z, (3 - 2*r)*z*(1 - z)]
        vz = [0._dp, 0._dp, 1._dp, r, (r - 1)*(2 - r)*(1 - 2*z)]
    end subroutine weights

    subroutine native_gradient(r, g, gn, mu0_jcell)
        real(dp), intent(in) :: r, g(2)
        real(dp), intent(out) :: gn(2), mu0_jcell
        real(dp) :: br, bz, bt
        call gs_field(r, g(1), g(2), 7._dp, br, bz, bt)
        call require_close(br, -g(2)/r, 'independent native BR sign/units')
        call require_close(bz, g(1)/r, 'independent native BZ sign/units')
        call require_close(bt, 7/r, 'independent native Bphi sign/units')
        gn = [r*bz, -r*br]
        ! Even with zero affine Hessian, -d(BZ)/dR = g_R/R^2 survives.
        mu0_jcell = bz/r
    end subroutine native_gradient

    subroutine integrate(tri, val, reverse_edge, cell, jump, weak, boundary)
        integer, intent(in) :: tri(3, 2)
        real(dp), intent(in) :: val(4)
        logical, intent(in) :: reverse_edge
        real(dp), intent(out) :: cell(nw), jump(nw), weak(nw), boundary(nw)
        real(dp) :: q(32), w(32)
        real(dp) :: g(2, 2), a(2), b(2), c(2), u(2), v(2), det, centroid(2)
        real(dp) :: point(2), normal(2), gn(2), gn2(2), jcell, dummy
        real(dp) :: vv(nw), vr(nw), vz(nw), factor
        integer :: e, k, i, j, ia, ib, other
        cell = 0; jump = 0; weak = 0; boundary = 0
        call gauss_legendre_ab(32, 0._dp, 1._dp, q, w)
        do e = 1, 2
            a = xy(:, tri(1, e)); u = xy(:, tri(2, e)) - a
            v = xy(:, tri(3, e)) - a
            det = u(1)*v(2) - u(2)*v(1)
            g(1, e) = ((val(tri(2, e)) - val(tri(1, e)))*v(2) - &
                (val(tri(3, e)) - val(tri(1, e)))*u(2))/det
            g(2, e) = (u(1)*(val(tri(3, e)) - val(tri(1, e))) - &
                v(1)*(val(tri(2, e)) - val(tri(1, e))))/det
            do j = 1, size(q)
                do i = 1, size(q)
                    point = a + q(i)*u + (1 - q(i))*q(j)*v
                    factor = abs(det)*(1 - q(i))*w(i)*w(j)
                    call weights(point(1), point(2), vv, vr, vz)
                    call native_gradient(point(1), g(:, e), gn, jcell)
                    cell = cell + factor*jcell*vv
                    weak = weak + factor*(gn(1)*vr + gn(2)*vz)/point(1)
                end do
            end do
        end do
        do e = 1, 2
            centroid = sum(xy(:, tri(:, e)), dim=2)/3
            other = 3 - e
            do k = 1, 3
                ia = tri(k, e); ib = tri(mod(k, 3) + 1, e)
                a = xy(:, ia); b = xy(:, ib)
                normal = [b(2) - a(2), a(1) - b(1)]
                if (dot_product(normal, (a + b)/2 - centroid) < 0) normal = -normal
                if (reverse_edge) then
                    c = a; a = b; b = c
                end if
                if (any(tri(:, other) == ia)) then
                    if (any(tri(:, other) == ib)) then
                        if (e == 2) cycle
                        do i = 1, size(q)
                            point = a + q(i)*(b - a)
                            call weights(point(1), point(2), vv, vr, vz)
                            call native_gradient(point(1), g(:, e), gn, dummy)
                            call native_gradient(point(1), g(:, other), gn2, dummy)
                            jump = jump - w(i)*dot_product(gn2 - gn, normal)* &
                                vv/point(1)
                        end do
                        cycle
                    end if
                end if
                do i = 1, size(q)
                    point = a + q(i)*(b - a)
                    call weights(point(1), point(2), vv, vr, vz)
                    call native_gradient(point(1), g(:, e), gn, dummy)
                    boundary = boundary + w(i)*dot_product(gn, normal)*vv/point(1)
                end do
            end do
        end do
    end subroutine integrate
end program test_gs_p1_weak_current
