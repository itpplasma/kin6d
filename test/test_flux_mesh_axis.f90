program test_flux_mesh_axis
    !! flux_mesh_t geometry near the magnetic axis.  A P_k nodal interpolant
    !! of a smooth function on the isoparametric mesh must converge at order
    !! k + 1 over the whole disk, axis fan included; a non-polynomial element
    !! map in the fan triangles caps it at a lower order.  The domain area
    !! checks the geometry itself (the circle is exact on the outer ring).
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_mhd_mesh, only: flux_mesh_t, flux_mesh_build, flux_mesh_geometry
    use fortfem_triangle_lagrange_arbitrary_order, only: &
        triangle_lagrange_nodes, evaluate_triangle_lagrange_basis
    implicit none

    real(dp), parameter :: pi = 3.14159265358979323846_dp
    integer, parameter :: nlev = 5
    real(dp) :: err(nlev), area_(nlev), rate
    integer :: geo, lev, nr
    logical :: ok

    ok = .true.
    do geo = 1, 3
        do lev = 1, nlev
            nr = 4*2**(lev - 1)
            call interp_error(nr, geo, err(lev), area_(lev))
        end do
        rate = log(err(nlev - 1)/err(nlev))/log(2.0_dp)
        print '(a, i0, a, 5es11.3, a, f6.2)', 'geo ', geo, ' L2 err', err, &
            ' order', rate
        ok = ok .and. rate > real(geo + 1, dp) - 0.35_dp
        if (geo >= 2) ok = ok .and. abs(area_(nlev) - pi) < 1.0e-4_dp
    end do
    if (.not. ok) error stop 'test_flux_mesh_axis FAILED'
    print '(a)', 'test_flux_mesh_axis passed'

contains

    subroutine circle_map(u, theta, x)
        real(dp), intent(in) :: u, theta
        real(dp), intent(out) :: x(2)
        x = u*[cos(theta), sin(theta)]
    end subroutine circle_map

    pure real(dp) function f(x)
        real(dp), intent(in) :: x(2)
        f = exp(x(1))*cos(2.0_dp*x(2)) + x(1)*x(2)
    end function f

    subroutine interp_error(nr, geo, e, area)
        integer, intent(in) :: nr, geo
        real(dp), intent(out) :: e, area
        type(flux_mesh_t) :: fm
        real(dp), allocatable :: u(:), lat(:, :), gv(:), gg(:, :), fn(:)
        real(dp) :: jac(2, 2), pos(2), det, w, xi, eta, uh
        integer :: st, t, k, iq, jq, nn, nq
        real(dp), parameter :: gx(4) = [0.06943184420297371_dp, &
            0.33000947820757187_dp, 0.66999052179242813_dp, 0.93056815579702629_dp]
        real(dp), parameter :: gw(4) = [0.17392742256872693_dp, &
            0.32607257743127307_dp, 0.32607257743127307_dp, 0.17392742256872693_dp]

        allocate (u(nr))
        do k = 1, nr
            u(k) = real(k, dp)/nr
        end do
        call flux_mesh_build(u, 8, 4096, geo, circle_map, fm, st)
        if (st /= 0) error stop 'mesh build failed'
        nn = size(fm%mm%nodes, 2)
        allocate (fn(nn), gv(nn), gg(2, nn))
        e = 0.0_dp
        area = 0.0_dp
        nq = size(gx)
        do t = 1, fm%mm%topology%n_triangles
            do k = 1, nn
                fn(k) = f(fm%mm%nodes(:, k, t))
            end do
            ! Gauss-Legendre product rule collapsed onto the triangle
            do iq = 1, nq
                do jq = 1, nq
                    xi = gx(iq)*(1.0_dp - gx(jq))
                    eta = gx(jq)
                    w = gw(iq)*gw(jq)*(1.0_dp - gx(jq))
                    call flux_mesh_geometry(fm, t, xi, eta, pos, jac)
                    det = jac(1, 1)*jac(2, 2) - jac(1, 2)*jac(2, 1)
                    call evaluate_triangle_lagrange_basis(fm%gb, xi, eta, gv, gg, st)
                    uh = dot_product(fn, gv)
                    e = e + w*abs(det)*(f(pos) - uh)**2
                    area = area + w*abs(det)
                end do
            end do
        end do
        e = sqrt(e)
    end subroutine interp_error

end program test_flux_mesh_axis
