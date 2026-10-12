module kin6d_mhd_mesh
    !! K-M1 (first stage): flux-aligned meridional meshes for the per-n linear
    !! response.  The mesh is a nest of rings around the magnetic axis.
    !!
    !! Topology.  Ring j = 1 .. nring sits at the radial label u_j
    !! (0 < u_1 < ... < u_nring) and carries nth_j vertices at the angles
    !! theta_i = 2 pi (i - 1)/nth_j; vertex 1 is the axis.  nth_j is n_core
    !! times a power of two, grows at most by a factor two from ring to ring
    !! (the outer vertices of a doubling strip include the inner angles) and
    !! keeps the cells of a uniform u ladder near unit aspect ratio.  The fan
    !! of ring 1 around the axis and the strips between consecutive rings are
    !! triangulated counter-clockwise.  Every ring edge (theta_i -> theta_i+1)
    !! is local edge 1 of the triangle on its inner side, so ring edges,
    !! traversed with increasing theta, have that triangle on their left; for
    !! the last ring this is the counter-clockwise domain boundary.
    !!
    !! Geometry.  The physical map x(u, theta) is supplied by the caller (a
    !! flux-surface map for G-EQDSK input, or an analytic map for tests).  The
    !! isoparametric nodes of order geo are the images of the equispaced
    !! lattice of each triangle in the polar parameters (u, theta), so every
    !! node of a ring edge lies on the ring curve x(u_j, .).  In the fan
    !! triangles the axis has no angle; their lattice uses the collapsed polar
    !! parameters u = u_1 (l_2 + l_3), theta = (l_2 theta_a + l_3 theta_b)/
    !! (l_2 + l_3), which agree on shared edges with both neighbours.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortfem_mesh_2d, only: mesh_2d_t
    use fortfem_meridian_curved_mesh, only: meridian_mesh_t
    use fortfem_triangle_lagrange_arbitrary_order, only: &
        initialize_triangle_lagrange_basis, triangle_lagrange_basis_t, &
        triangle_lagrange_nodes, evaluate_triangle_lagrange_basis
    implicit none
    private

    real(dp), parameter :: pi = 3.14159265358979323846_dp

    abstract interface
        subroutine polar_map_i(u, theta, x)
            !! Physical point x = (R, Z) of the radial label u and angle
            !! theta; x(0, theta) is the axis for every theta.
            import :: dp
            real(dp), intent(in) :: u, theta
            real(dp), intent(out) :: x(2)
        end subroutine polar_map_i
    end interface

    type, public :: flux_mesh_t
        type(meridian_mesh_t) :: mm
        integer :: nring = 0
        real(dp), allocatable :: u(:)          !! ring labels
        integer, allocatable :: nth(:)         !! vertices per ring
        integer, allocatable :: vfirst(:)      !! first vertex of ring j
        integer, allocatable :: tlast(:)       !! last triangle of strip j
        integer, allocatable :: roff(:)        !! ring edge offset of ring j
        integer, allocatable :: ring_tri(:)    !! inner triangle of ring edge
        type(triangle_lagrange_basis_t) :: gb  !! geometry basis
    end type flux_mesh_t

    public :: polar_map_i
    public :: flux_mesh_build
    public :: flux_mesh_ring_edge
    public :: flux_mesh_locate
    public :: flux_mesh_geometry

contains

    subroutine flux_mesh_build(u, n_core, nth_max, geo, map, fm, status)
        real(dp), intent(in) :: u(:)
        integer, intent(in) :: n_core, nth_max, geo
        procedure(polar_map_i) :: map
        type(flux_mesh_t), intent(out) :: fm
        integer, intent(out) :: status

        type(mesh_2d_t) :: topo
        integer :: nr, j, i, nv, nt, ratio, a, b, c, t, ne
        real(dp) :: x(2)

        status = 1
        nr = size(u)
        if (nr < 1 .or. n_core < 3 .or. geo < 1 .or. geo > 3) return
        if (u(1) <= 0.0_dp) return
        do j = 2, nr
            if (u(j) <= u(j - 1)) return
        end do
        fm%nring = nr
        fm%u = u
        allocate (fm%nth(nr), fm%vfirst(nr), fm%tlast(nr), fm%roff(nr + 1))
        fm%nth(1) = n_core
        do j = 2, nr
            fm%nth(j) = fm%nth(j - 1)
            if (2*fm%nth(j - 1) <= nth_max .and. &
                2.0_dp*pi*u(j)/fm%nth(j - 1) > 1.5_dp*(u(j) - u(j - 1))) &
                fm%nth(j) = 2*fm%nth(j - 1)
        end do
        nv = 1
        do j = 1, nr
            fm%vfirst(j) = nv + 1
            nv = nv + fm%nth(j)
        end do
        nt = fm%nth(1)
        do j = 2, nr
            nt = nt + fm%nth(j) + fm%nth(j - 1)
        end do
        topo%n_vertices = nv
        topo%n_triangles = nt
        topo%has_triangles = .true.
        allocate (topo%vertices(2, nv), topo%triangles(3, nt))
        call map(0.0_dp, 0.0_dp, x)
        topo%vertices(:, 1) = x
        do j = 1, nr
            do i = 1, fm%nth(j)
                call map(u(j), theta_of(fm, j, i), x)
                topo%vertices(:, fm%vfirst(j) + i - 1) = x
            end do
        end do
        ne = 0
        do j = 1, nr
            fm%roff(j) = ne
            ne = ne + fm%nth(j)
        end do
        fm%roff(nr + 1) = ne
        allocate (fm%ring_tri(ne))
        t = 0
        ! axis fan
        do i = 1, fm%nth(1)
            t = t + 1
            topo%triangles(:, t) = [vid(fm, 1, i), vid(fm, 1, i + 1), 1]
            fm%ring_tri(fm%roff(1) + i) = t
        end do
        fm%tlast(1) = t
        do j = 2, nr
            ratio = fm%nth(j)/fm%nth(j - 1)
            do i = 1, fm%nth(j - 1)
                if (ratio == 1) then
                    a = vid(fm, j, i)
                    b = vid(fm, j, i + 1)
                    c = vid(fm, j - 1, i)
                    t = t + 1
                    topo%triangles(:, t) = [a, b, c]
                    fm%ring_tri(fm%roff(j) + i) = t
                    t = t + 1
                    topo%triangles(:, t) = [c, b, vid(fm, j - 1, i + 1)]
                else
                    a = vid(fm, j, 2*i - 1)
                    b = vid(fm, j, 2*i)
                    c = vid(fm, j, 2*i + 1)
                    t = t + 1
                    topo%triangles(:, t) = [a, b, vid(fm, j - 1, i)]
                    fm%ring_tri(fm%roff(j) + 2*i - 1) = t
                    t = t + 1
                    topo%triangles(:, t) = [b, c, vid(fm, j - 1, i + 1)]
                    fm%ring_tri(fm%roff(j) + 2*i) = t
                    t = t + 1
                    topo%triangles(:, t) = [vid(fm, j - 1, i), b, &
                        vid(fm, j - 1, i + 1)]
                end if
            end do
            fm%tlast(j) = t
        end do
        if (t /= nt) return
        call build_nodes(fm, topo, geo, map, status)
    end subroutine flux_mesh_build

    pure real(dp) function theta_of(fm, j, i)
        type(flux_mesh_t), intent(in) :: fm
        integer, intent(in) :: j, i
        theta_of = 2.0_dp*pi*real(i - 1, dp)/real(fm%nth(j), dp)
    end function theta_of

    pure integer function vid(fm, j, i)
        !! Vertex of ring j, angle index i (periodic).
        type(flux_mesh_t), intent(in) :: fm
        integer, intent(in) :: j, i
        vid = fm%vfirst(j) + modulo(i - 1, fm%nth(j))
    end function vid

    subroutine vertex_polar(fm, v, ur, th)
        !! Polar parameters of vertex v (axis: u = 0, th = 0).
        type(flux_mesh_t), intent(in) :: fm
        integer, intent(in) :: v
        real(dp), intent(out) :: ur, th
        integer :: j

        ur = 0.0_dp
        th = 0.0_dp
        if (v == 1) return
        do j = fm%nring, 1, -1
            if (v >= fm%vfirst(j)) exit
        end do
        ur = fm%u(j)
        th = theta_of(fm, j, v - fm%vfirst(j) + 1)
    end subroutine vertex_polar

    subroutine build_nodes(fm, topo, geo, map, status)
        type(flux_mesh_t), intent(inout) :: fm
        type(mesh_2d_t), intent(in) :: topo
        integer, intent(in) :: geo
        procedure(polar_map_i) :: map
        integer, intent(out) :: status

        type(triangle_lagrange_basis_t) :: basis
        real(dp), allocatable :: lat(:, :)
        real(dp) :: up(3), tp(3), l(3), uu, tt, s
        logical :: ring(3)
        integer :: t, i, k

        call initialize_triangle_lagrange_basis(geo, basis, status)
        if (status /= 0) return
        fm%gb = basis
        call triangle_lagrange_nodes(basis, lat)
        fm%mm%topology = topo
        fm%mm%geo_order = geo
        allocate (fm%mm%nodes(2, size(lat, 2), topo%n_triangles))
        do t = 1, topo%n_triangles
            do k = 1, 3
                call vertex_polar(fm, topo%triangles(k, t), up(k), tp(k))
            end do
            ! unwrap the angles of the triangle to one branch (the axis
            ! vertex, vertex 1, has no angle)
            ring = topo%triangles(:, t) /= 1
            do k = 1, 3
                if (.not. ring(k)) cycle
                if (tp(k) < maxval(tp, mask=ring) - pi) tp(k) = tp(k) + 2.0_dp*pi
            end do
            do i = 1, size(lat, 2)
                l = [1.0_dp - lat(1, i) - lat(2, i), lat(1, i), lat(2, i)]
                if (topo%triangles(3, t) == 1) then
                    s = l(1) + l(2)
                    uu = up(1)*s
                    if (s > 0.0_dp) then
                        tt = (l(1)*tp(1) + l(2)*tp(2))/s
                    else
                        tt = 0.0_dp
                    end if
                else
                    uu = sum(l*up)
                    tt = sum(l*tp)
                end if
                call map(uu, tt, fm%mm%nodes(:, i, t))
            end do
        end do
        status = 0
    end subroutine build_nodes

    pure subroutine flux_mesh_ring_edge(fm, j, i, t, theta_a, theta_b)
        !! Triangle t of ring edge i of ring j (its local edge 1) and the
        !! angles at its ends (theta_b > theta_a, unwrapped).
        type(flux_mesh_t), intent(in) :: fm
        integer, intent(in) :: j, i
        integer, intent(out) :: t
        real(dp), intent(out) :: theta_a, theta_b

        t = fm%ring_tri(fm%roff(j) + i)
        theta_a = theta_of(fm, j, i)
        theta_b = theta_a + 2.0_dp*pi/real(fm%nth(j), dp)
    end subroutine flux_mesh_ring_edge

    subroutine flux_mesh_locate(fm, x, ulab, theta, t, xi, eta, status)
        !! Triangle t and reference point (xi, eta) of the physical point x
        !! with polar parameters (ulab, theta) (ulab may be any label inside
        !! the strip of x, theta its angle parameter).  Newton inversion of
        !! the isoparametric map over the candidate triangles of the strip,
        !! projected onto the reference triangle; a point slightly outside the
        !! polygonal boundary of a straight-sided mesh (distance below 0.05
        !! cell diameters) is moved to the nearest point of the nearest cell.
        type(flux_mesh_t), intent(in) :: fm
        real(dp), intent(in) :: x(2), ulab, theta
        integer, intent(out) :: t
        real(dp), intent(out) :: xi, eta
        integer, intent(out) :: status

        integer :: j, nin, per, i0, di, c, tc, it
        real(dp) :: th, v(2), jac(2, 2), r(2), det, d(2), tol, best, viol, xq, eq

        status = 1
        do j = 1, fm%nring
            if (fm%u(j) >= ulab) exit
        end do
        if (j > fm%nring) return
        th = modulo(theta, 2.0_dp*pi)
        if (j == 1) then
            nin = fm%nth(1)
            per = 1
        else
            nin = fm%nth(j - 1)
            per = 2 + fm%nth(j)/fm%nth(j - 1) - 1
        end if
        i0 = int(th/(2.0_dp*pi)*nin)
        tol = 1.0e-11_dp
        best = huge(1.0_dp)
        do di = 0, 2
            do c = 1, per
                if (j == 1) then
                    if (c > 1) cycle
                    tc = modulo(i0 + di - 1 + nin, nin) + 1
                else
                    tc = fm%tlast(j - 1) + per*modulo(i0 + di - 1 + nin, nin) + c
                end if
                xq = 1.0_dp/3.0_dp
                eq = 1.0_dp/3.0_dp
                do it = 1, 40
                    call flux_mesh_geometry(fm, tc, xq, eq, v, jac)
                    r = v - x
                    det = jac(1, 1)*jac(2, 2) - jac(1, 2)*jac(2, 1)
                    if (abs(det) <= tiny(1.0_dp)) exit
                    d(1) = (jac(2, 2)*r(1) - jac(1, 2)*r(2))/det
                    d(2) = (-jac(2, 1)*r(1) + jac(1, 1)*r(2))/det
                    ! Newton step projected onto the closed reference triangle
                    xq = max(0.0_dp, xq - d(1))
                    eq = max(0.0_dp, eq - d(2))
                    if (xq + eq > 1.0_dp) then
                        xq = xq/(xq + eq)
                        eq = 1.0_dp - xq
                    end if
                    if (abs(d(1)) + abs(d(2)) < 1.0e-14_dp) exit
                end do
                call flux_mesh_geometry(fm, tc, xq, eq, v, jac)
                viol = norm2(v - x)/diameter(fm, tc)
                if (viol < best) then
                    best = viol
                    t = tc
                    xi = xq
                    eta = eq
                end if
                if (viol <= tol) then
                    status = 0
                    return
                end if
            end do
        end do
        ! a point just outside a straight-sided boundary cell (geometry order
        ! 1 on a curved boundary) is moved to the nearest boundary point
        if (best < 0.05_dp) status = 0
    contains
        real(dp) function diameter(fm, t)
            type(flux_mesh_t), intent(in) :: fm
            integer, intent(in) :: t
            integer :: a(3)
            a = fm%mm%topology%triangles(:, t)
            associate (v => fm%mm%topology%vertices)
                diameter = max(norm2(v(:, a(1)) - v(:, a(2))), &
                    norm2(v(:, a(2)) - v(:, a(3))), norm2(v(:, a(3)) - v(:, a(1))))
            end associate
        end function diameter
    end subroutine flux_mesh_locate

    subroutine flux_mesh_geometry(fm, t, xi, eta, x, jac)
        !! Physical point and Jacobian jac(:, c) = dx/dxi_c of triangle t at
        !! the reference point (xi, eta) (inside the closed reference
        !! triangle).
        type(flux_mesh_t), intent(in) :: fm
        integer, intent(in) :: t
        real(dp), intent(in) :: xi, eta
        real(dp), intent(out) :: x(2), jac(2, 2)
        real(dp) :: gv(size(fm%mm%nodes, 2)), gg(2, size(fm%mm%nodes, 2))
        integer :: k, st

        call evaluate_triangle_lagrange_basis(fm%gb, xi, eta, gv, gg, st)
        x = 0.0_dp
        jac = 0.0_dp
        do k = 1, size(gv)
            x = x + gv(k)*fm%mm%nodes(:, k, t)
            jac(:, 1) = jac(:, 1) + gg(1, k)*fm%mm%nodes(:, k, t)
            jac(:, 2) = jac(:, 2) + gg(2, k)*fm%mm%nodes(:, k, t)
        end do
    end subroutine flux_mesh_geometry

end module kin6d_mhd_mesh
